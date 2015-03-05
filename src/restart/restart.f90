#include "boltzplatz.h"

MODULE MOD_Restart
!===================================================================================================================================
! Module to handle Boltzplatz's restart
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitRestart
  MODULE PROCEDURE InitRestart
END INTERFACE

INTERFACE Restart
  MODULE PROCEDURE Restart
END INTERFACE

INTERFACE FinalizeRestart
  MODULE PROCEDURE FinalizeRestart
END INTERFACE

PUBLIC :: InitRestart,FinalizeRestart
PUBLIC :: Restart
!===================================================================================================================================

CONTAINS

SUBROUTINE InitRestart()
!===================================================================================================================================
! Initialize all necessary information to perform the restart
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Interpolation_Vars, ONLY: xGP,InterpolationInitIsDone
USE MOD_Restart_Vars
USE MOD_HDF5_Input,ONLY:OpenDataFile,CloseDataFile,GetDataProps,ReadAttribute,File_ID
USE MOD_ReadInTools,ONLY:GETLOGICAL,GETREALARRAY,ReadInDone 
#ifdef PARTICLES
USE MOD_DSMC_Vars,ONLY: UseDSMC
#endif /*PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: nArgs
! Special handling for DSMC chemistry data-file
INTEGER            :: maxNArgs
!===================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.RestartInitIsDone)THEN
   CALL abort(__STAMP__,'InitRestart not ready to be called or already called.',999,999.)
   RETURN
END IF

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT RESTART...'

#ifdef PARTICLES
! DSMC handling:
useDSMC=GETLOGICAL('UseDSMC','.FALSE.')
IF (useDSMC) THEN
  ReadInDone = .FALSE.
  maxNArgs = 3
ELSE
  maxNArgs = 2
END IF
#else
maxNArgs=2
#endif /*PARTICLES*/


! Check if we want to perform a restart
nArgs=COMMAND_ARGUMENT_COUNT()
IF (nArgs .EQ. maxNArgs) THEN
  ! Read in the state file we want to restart from
  CALL GETARG(maxNArgs,RestartFile)
  SWRITE(UNIT_StdOut,'(A,A,A)')' | Restarting from file "',TRIM(RestartFile),'":'
  DoRestart = .TRUE.
#ifdef MPI
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.)
#else
  CALL OpenDataFile(RestartFile,create=.FALSE.)
#endif
  CALL GetDataProps(nVar_Restart,N_Restart,nElems_Restart,NodeType_Restart)
  ! Read in time from restart file
  CALL ReadAttribute(File_ID,'Time',1,RealScalar=RestartTime)
  CALL CloseDataFile() 
  PrimScaling=GETLOGICAL('PrimScaling','.FALSE.')
  IF(PrimScaling) PrimScale=GETREALARRAY('PrimScale',PP_nVar)
ELSE
  RestartTime = 0.
  SWRITE(UNIT_StdOut,'(A)')' | No restart wanted, doing a fresh computation!'
END IF
IF(DoRestart .AND. (N_Restart .NE. PP_N))THEN
  BuildNewMesh       =.TRUE.
  WriteNewMesh       =.TRUE.
  InterpolateSolution=.TRUE.
END IF

IF(InterpolateSolution)THEN
  CALL initRestartBasis(PP_N,N_Restart,xGP)
END IF

RestartInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT RESTART DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitRestart



SUBROUTINE InitRestartBasis(N_in,N_Restart_in,xGP)
!===================================================================================================================================
! Initialize all necessary information to perform the restart
!===================================================================================================================================
! MODULES
USE MOD_Restart_Vars, ONLY:Vdm_GaussNRestart_GaussN
USE MOD_Basis,        ONLY:LegendreGaussNodesAndWeights,LegGaussLobNodesAndWeights,ChebyGaussLobNodesAndWeights
USE MOD_Basis,        ONLY:BarycentricWeights,InitializeVandermonde
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                  :: N_in,N_Restart_in
REAL,INTENT(IN),DIMENSION(0:N_in)   :: xGP
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(0:N_Restart_in)  :: xGP_Restart,wBary_Restart
!===================================================================================================================================
  ALLOCATE(Vdm_GaussNRestart_GaussN(0:N_in,0:N_Restart_in))
#if (PP_NodeType==1)
  CALL LegendreGaussNodesAndWeights(N_Restart_in,xGP_Restart)
#elif (PP_NodeType==2)
  CALL LegGaussLobNodesAndWeights(N_Restart_in,xGP_Restart)
#elif (PP_NodeType==3)
  CALL ChebyGaussLobNodesAndWeights(N_Restart_in,xGP_Restart)
#endif
  CALL BarycentricWeights(N_Restart_in,xGP_Restart,wBary_Restart)
  CALL InitializeVandermonde(N_Restart_in,N_in,wBary_Restart,xGP_Restart,xGP,Vdm_GaussNRestart_GaussN)
END SUBROUTINE InitRestartBasis



SUBROUTINE Restart()
!===================================================================================================================================
! Read in mesh (if available) and state, set time for restart
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,         ONLY:U
USE MOD_Mesh_Vars,       ONLY:offsetElem
USE MOD_Restart_Vars,    ONLY:Vdm_GaussNRestart_GaussN
USE MOD_Restart_Vars,    ONLY:DoRestart,N_Restart,RestartFile,RestartTime,InterpolateSolution
USE MOD_ChangeBasis,     ONLY:ChangeBasis3D
USE MOD_HDF5_input ,     ONLY:OpenDataFile,CloseDataFile,ReadArray,ReadAttribute
USE MOD_HDF5_Output,     ONLY:FlushHDF5
USE MOD_PML_Vars,        ONLY:DoPML,PMLToElem,U2,nPMLElems
#ifdef PP_POIS
USE MOD_Equation_Vars,   ONLY:Phi
#endif
#ifdef PARTICLES
USE MOD_Particle_Vars,   ONLY:PartState, PartSpecies, PEM, PDM, Species, nSpecies, usevMPF, PartMPF
USE MOD_part_tools,      ONLY: UpdateNextFreePosition
USE MOD_DSMC_Vars,       ONLY: UseDSMC, CollisMode,PartStateIntEn, DSMC
!USE MOD_BoundaryTools,   ONLY : SingleParticleToExactElement, ParticleInsideQuad3D
!USE MOD_BoundaryTools,   ONLY: ParticleInsideQuad3D
USE MOD_Eval_XYZ,        ONLY: EVal_xyz_ElemCheck
USE MOD_Particle_Mesh,   ONLY:SingleParticleToExactElement
#ifdef MPI
USE MOD_Particle_MPI_Vars,ONLY:PartMPI
#endif
#endif /*PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE         :: U_local(:,:,:,:,:)
REAL,ALLOCATABLE         :: U_local2(:,:,:,:,:)
INTEGER                  :: iElem,iPML
#ifdef PARTICLES
INTEGER                  :: FirstElemInd,LastelemInd,i,iInit
INTEGER,ALLOCATABLE      :: PartInt(:,:)
INTEGER,PARAMETER        :: PartIntSize=2        !number of entries in each line of PartInt
INTEGER                  :: PartDataSize         !number of entries in each line of PartData
INTEGER                  :: locnPart,offsetnPart
INTEGER,PARAMETER        :: ELEM_FirstPartInd=1
INTEGER,PARAMETER        :: ELEM_LastPartInd=2
REAL,ALLOCATABLE         :: PartData(:,:)
REAL                     :: xi(3)
LOGICAL                  :: InElementCheck
REAL                     :: det(16)
INTEGER                  :: COUNTER, COUNTER2
REAL,PARAMETER           :: eps=1e-8 ! same value as in eval_xyz_elem
REAL                     :: epsOne
#ifdef MPI
REAL, ALLOCATABLE        :: SendBuff(:), RecBuff(:)
INTEGER                  :: LostParts(0:PartMPI%nProcs-1), Displace(0:PartMPI%nProcs-1),CurrentPartNum
INTEGER                  :: NbrOfFoundParts, CompleteNbrOfFound, RecCount(0:PartMPI%nProcs-1)
#endif /*MPI*/
#endif /*PARTICLES*/
!===================================================================================================================================
IF(DoRestart)THEN
SWRITE(UNIT_stdOut,*)'Restarting from File:',TRIM(RestartFile)
#ifdef MPI
   CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.)
#else
   CALL OpenDataFile(RestartFile,create=.FALSE.)
#endif
  ! Read in time from restart file
  !CALL ReadAttribute(File_ID,'Time',1,RealScalar=RestartTime)
  ! Read in state
  IF(.NOT. InterpolateSolution)THEN
    ! No interpolation needed, read solution directly from file
#ifdef PP_POIS
#if (PP_nVar==8)
    CALL ReadArray('DG_SolutionE',5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),OffsetElem,5,RealArray=U)
    CALL ReadArray('DG_SolutionPhi',5,(/4,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),OffsetElem,5,RealArray=Phi)
#else
    CALL ReadArray('DG_SolutionE',5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),OffsetElem,5,RealArray=U)
    CALL ReadArray('DG_SolutionPhi',5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),OffsetElem,5,RealArray=Phi)
#endif
#else
    CALL ReadArray('DG_Solution',5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),OffsetElem,5,RealArray=U)
    IF(DoPML)THEN
      ALLOCATE(U_local(6,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
      CALL ReadArray('PML_Solution',5,(/6,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),OffsetElem,5,RealArray=U_local)
      DO iPML=1,nPMLElems
        U2(:,:,:,:,iPML) = U_local(:,:,:,:,PMLToElem(iPML))
      END DO ! iPML
      DEALLOCATE(U_local)
    END IF ! DoPML
#endif
    !CALL ReadState(RestartFile,PP_nVar,PP_N,PP_nElems,U)
  ELSE
    ! We need to interpolate the solution to the new computational grid
    SWRITE(UNIT_stdOut,*)'Interpolating solution from restart grid with N=',N_restart,' to computational grid with N=',PP_N
#ifdef PP_POIS
#if (PP_nVar==8)
    ALLOCATE(U_local(PP_nVar,0:N_Restart,0:N_Restart,0:N_Restart,PP_nElems))
    CALL ReadArray('DG_SolutionE',5,(/PP_nVar,N_Restart+1,N_Restart+1,N_Restart+1,PP_nElems/),OffsetElem,5,RealArray=U_local)
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(PP_nVar,N_Restart,PP_N,Vdm_GaussNRestart_GaussN,U_local(:,:,:,:,iElem),U(:,:,:,:,iElem))
    END DO
    DEALLOCATE(U_local)

    ALLOCATE(U_local(4,0:N_Restart,0:N_Restart,0:N_Restart,PP_nElems))
    CALL ReadArray('DG_SolutionPhi',5,(/4,N_Restart+1,N_Restart+1,N_Restart+1,PP_nElems/),OffsetElem,5,RealArray=U_local)
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(4,N_Restart,PP_N,Vdm_GaussNRestart_GaussN,U_local(:,:,:,:,iElem),Phi(:,:,:,:,iElem))
    END DO
    DEALLOCATE(U_local)
#else
    ALLOCATE(U_local(PP_nVar,0:N_Restart,0:N_Restart,0:N_Restart,PP_nElems))
    CALL ReadArray('DG_SolutionE',5,(/PP_nVar,N_Restart+1,N_Restart+1,N_Restart+1,PP_nElems/),OffsetElem,5,RealArray=U_local)
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(PP_nVar,N_Restart,PP_N,Vdm_GaussNRestart_GaussN,U_local(:,:,:,:,iElem),U(:,:,:,:,iElem))
    END DO
    CALL ReadArray('DG_SolutionPhi',5,(/PP_nVar,N_Restart+1,N_Restart+1,N_Restart+1,PP_nElems/),OffsetElem,5,RealArray=U_local)
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(PP_nVar,N_Restart,PP_N,Vdm_GaussNRestart_GaussN,U_local(:,:,:,:,iElem),Phi(:,:,:,:,iElem))
    END DO
    DEALLOCATE(U_local)
#endif
#else
    ALLOCATE(U_local(PP_nVar,0:N_Restart,0:N_Restart,0:N_Restart,PP_nElems))
    CALL ReadArray('DG_Solution',5,(/PP_nVar,N_Restart+1,N_Restart+1,N_Restart+1,PP_nElems/),OffsetElem,5,RealArray=U_local)
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(PP_nVar,N_Restart,PP_N,Vdm_GaussNRestart_GaussN,U_local(:,:,:,:,iElem),U(:,:,:,:,iElem))
    END DO
    DEALLOCATE(U_local)
    IF(DoPML)THEN
      ALLOCATE(U_local(6,0:N_Restart,0:N_Restart,0:N_Restart,PP_nElems))
      ALLOCATE(U_local2(6,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
      CALL ReadArray('PML_Solution',5,(/6,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),OffsetElem,5,RealArray=U_local)
      DO iElem=1,PP_nElems
        CALL ChangeBasis3D(PP_nVar,N_Restart,PP_N,Vdm_GaussNRestart_GaussN,U_local(:,:,:,:,iElem),U_local2(:,:,:,:,iElem))
      END DO
      DO iPML=1,nPMLElems
        U2(:,:,:,:,iPML) = U_local2(:,:,:,:,PMLToElem(iPML))
      END DO ! iPML
      DEALLOCATE(U_local,U_local2)
    END IF ! DoPML
#endif
    SWRITE(UNIT_stdOut,*)'DONE!'
  END IF

#ifdef PARTICLES
  epsOne=1.0+eps
  IF (useDSMC) THEN
    IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicState)) THEN !int ener + 3, vmpf +1
      PartDataSize=11
    ELSE IF ((CollisMode.GT.1).AND.((usevMPF) .OR. (DSMC%ElectronicState)) ) THEN ! int ener + 2 and vmpf +1
                                                                                  ! or int ener + 3 and no vmpf
      PartDataSize=10
    ELSE IF (CollisMode.GT.1) THEN
      PartDataSize=9 !int ener + 2
    ELSE IF ( usevMPF) THEN 
      PartDataSize=8 ! vMPF + 1
    ELSE 
      PartDataSize=7 !+ 0
    END IF
  ELSE IF (usevMPF) THEN
    PartDataSize=8 !vmpf +1
  ELSE
    PartDataSize=7
  END IF  

  SWRITE(UNIT_stdOut,*)'Reading Particles from Restartfile...' 
  !read local ElemInfo from HDF5
  FirstElemInd=offsetElem+1
  LastElemInd=offsetElem+PP_nElems
  ! read local ParticleInfo from HDF5
  ALLOCATE(PartInt(FirstElemInd:LastElemInd,PartIntSize))
  CALL ReadArray('PartInt',2,(/PP_nElems,PartIntSize/),offsetElem,1,IntegerArray=PartInt)
  ! read local Particle Data from HDF5
  locnPart=PartInt(LastElemInd,ELEM_LastPartInd)-PartInt(FirstElemInd,ELEM_FirstPartInd)
  offsetnPart=PartInt(FirstElemInd,ELEM_FirstPartInd)
  ALLOCATE(PartData(offsetnPart+1:offsetnPart+locnPart,PartDataSize))
  CALL ReadArray('PartData',2,(/locnPart,PartDataSize/),offsetnPart,1,RealArray=PartData)!,&
                         !xfer_mode_independent=.TRUE.)
  IF (locnPart.GT.0) THEN
    PartState(1:locnPart,1)   = PartData(offsetnPart+1:offsetnPart+locnPart,1)
    PartState(1:locnPart,2)   = PartData(offsetnPart+1:offsetnPart+locnPart,2)
    PartState(1:locnPart,3)   = PartData(offsetnPart+1:offsetnPart+locnPart,3)
    PartState(1:locnPart,4)   = PartData(offsetnPart+1:offsetnPart+locnPart,4)
    PartState(1:locnPart,5)   = PartData(offsetnPart+1:offsetnPart+locnPart,5)
    PartState(1:locnPart,6)   = PartData(offsetnPart+1:offsetnPart+locnPart,6)
    PartSpecies(1:locnPart)= INT(PartData(offsetnPart+1:offsetnPart+locnPart,7))
    IF (useDSMC) THEN
      IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicState)) THEN
        PartStateIntEn(1:locnPart,1)=PartData(offsetnPart+1:offsetnPart+locnPart,8)
        PartStateIntEn(1:locnPart,2)=PartData(offsetnPart+1:offsetnPart+locnPart,9)
        PartStateIntEn(1:locnPart,3)=PartData(offsetnPart+1:offsetnPart+locnPart,10)
        PartMPF(1:locnPart)=PartData(offsetnPart+1:offsetnPart+locnPart,11)
      ELSE IF ((CollisMode.GT.1).AND. (usevMPF)) THEN
        PartStateIntEn(1:locnPart,1)=PartData(offsetnPart+1:offsetnPart+locnPart,8)
        PartStateIntEn(1:locnPart,2)=PartData(offsetnPart+1:offsetnPart+locnPart,9)
        PartMPF(1:locnPart)=PartData(offsetnPart+1:offsetnPart+locnPart,10)
      ELSE IF ((CollisMode.GT.1).AND. (DSMC%ElectronicState)) THEN
        PartStateIntEn(1:locnPart,1)=PartData(offsetnPart+1:offsetnPart+locnPart,8)
        PartStateIntEn(1:locnPart,2)=PartData(offsetnPart+1:offsetnPart+locnPart,9)
        PartStateIntEn(1:locnPart,3)=PartData(offsetnPart+1:offsetnPart+locnPart,10)
      ELSE IF (CollisMode.GT.1) THEN
        PartStateIntEn(1:locnPart,1)=PartData(offsetnPart+1:offsetnPart+locnPart,8)
        PartStateIntEn(1:locnPart,2)=PartData(offsetnPart+1:offsetnPart+locnPart,9)
      ELSE IF (usevMPF) THEN
        PartMPF(1:locnPart)=PartData(offsetnPart+1:offsetnPart+locnPart,8)        
      END IF
    ELSE IF (usevMPF) THEN
      PartMPF(1:locnPart)=PartData(offsetnPart+1:offsetnPart+locnPart,8)
    END IF
    DO iElem=FirstElemInd,LastElemInd
      IF (PartInt(iElem,ELEM_LastPartInd).GT.PartInt(iElem,ELEM_FirstPartInd)) THEN
        PEM%Element(PartInt(iElem,ELEM_FirstPartInd)-offsetnPart+1 : &
                               PartInt(iElem,ELEM_LastPartInd) -offsetnPart)  = iElem-offsetElem
        PEM%LastElement(PartInt(iElem,ELEM_FirstPartInd)-offsetnPart+1 : &
                               PartInt(iElem,ELEM_LastPartInd) -offsetnPart)  = iElem-offsetElem
      END IF
    END DO
    PDM%ParticleInside(1:locnPart) = .TRUE.
  END IF
  DEALLOCATE(PartInt,PartData)
  PDM%ParticleVecLength = PDM%ParticleVecLength + locnPart
  CALL UpdateNextFreePosition()
  SWRITE(UNIT_stdOut,*)'DONE!' 
  DO i=1,nSpecies
    DO iInit = Species(i)%StartnumberOfInits, Species(i)%NumberOfInits
      Species(i)%Init(iInit)%InsertedParticle = INT(Species(i)%Init(iInit)%ParticleEmission * RestartTime,8)
    END DO
  END DO
  ! if ParticleVecLength GT maxParticleNumber: Stop
  IF (PDM%ParticleVecLength.GT.PDM%maxParticleNumber) THEN
    WRITE(*,*) 'ERROR: Number of Particles in Restart higher than MaxParticleNumber'
    STOP
  END IF

  ! Since the elementside-local node number are NOT persistant and dependent on the location
  ! of the MPI borders, all particle-element mappings need to be checked after a restart
  ! Step 1: Identify particles that are not in the element in which they were before the restart
  COUNTER = 0
  COUNTER2 = 0
  DO i = 1,PDM%ParticleVecLength
    CALL Eval_xyz_ElemCheck(PartState(i,1:3),Xi,PEM%Element(i))
    IF(ALL(ABS(Xi).LE.EpsOne)) THEN ! particle inside
      InElementCheck=.TRUE.
    ELSE
      InElementCheck=.FALSE.
    END IF
    IF (.NOT.InElementCheck) THEN  ! try to find them within MyProc
      COUNTER = COUNTER + 1
      CALL SingleParticleToExactElement(i)
      IF (.NOT.PDM%ParticleInside(i)) THEN
        COUNTER2 = COUNTER2 + 1
      ELSE
        PEM%LastElement(i) = PEM%Element(i)
      END IF
    END IF
  END DO
#ifdef MPI
  ! Step 2: All particles that are not found withing MyProc need to be communicated to the others and located there
  ! Combine number of lost particles of all processes and allocate variables
  CALL MPI_ALLGATHER(COUNTER2, 1, MPI_INTEGER, LostParts, 1, MPI_INTEGER, PartMPI%COMM, IERROR)
  IF (SUM(LostParts).GT.0) THEN
    ALLOCATE(SendBuff(1:COUNTER2*PartDataSize))
    ALLOCATE(RecBuff(1:SUM(LostParts)*PartDataSize))
    ! Fill SendBuffer
    COUNTER = 0
    DO i = 1, PDM%ParticleVecLength
      IF (.NOT.PDM%ParticleInside(i)) THEN
        SendBuff(COUNTER+1:COUNTER+6) = PartState(i,1:6)
        SendBuff(COUNTER+7)           = REAL(PartSpecies(i))
        IF (useDSMC) THEN
          IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicState)) THEN
            SendBuff(COUNTER+8)  = PartStateIntEn(i,1)
            SendBuff(COUNTER+9)  = PartStateIntEn(i,2)
            SendBuff(COUNTER+10) = PartMPF(i)
            SendBuff(COUNTER+11) = PartStateIntEn(i,3)
          ELSE IF ((CollisMode.GT.1).AND. (usevMPF)) THEN
            SendBuff(COUNTER+8)  = PartStateIntEn(i,1)
            SendBuff(COUNTER+9)  = PartStateIntEn(i,2)
            SendBuff(COUNTER+10) = PartMPF(i)
          ELSE IF ((CollisMode.GT.1).AND. (DSMC%ElectronicState)) THEN
            SendBuff(COUNTER+8)  = PartStateIntEn(i,1)
            SendBuff(COUNTER+9)  = PartStateIntEn(i,2)
            SendBuff(COUNTER+10) = PartStateIntEn(i,3)
          ELSE IF (CollisMode.GT.1) THEN
            SendBuff(COUNTER+8)  = PartStateIntEn(i,1)
            SendBuff(COUNTER+9)  = PartStateIntEn(i,2)
          ELSE IF (usevMPF) THEN
            SendBuff(COUNTER+8) = PartMPF(i)
          END IF
        ELSE IF (usevMPF) THEN
          SendBuff(COUNTER+8) = PartMPF(i)
        END IF
        COUNTER = COUNTER + PartDataSize
      END IF
    END DO
    ! Distribute lost particles to all procs
    COUNTER = 0
    DO i = 0, PartMPI%nProcs-1
      RecCount(i) = LostParts(i) * PartDataSize
      Displace(i) = COUNTER
      COUNTER = COUNTER + LostParts(i)*PartDataSize
    END DO
    CALL MPI_ALLGATHERV(SendBuff, PartDataSize*LostParts(PartMPI%MyRank), MPI_DOUBLE_PRECISION, &
         RecBuff, RecCount, Displace, MPI_DOUBLE_PRECISION, PartMPI%COMM, IERROR)
    ! Add them to particle list and check if they are in MyProcs domain
    NbrOfFoundParts = 0
    CurrentPartNum = PDM%ParticleVecLength+1
    COUNTER = 0
    DO i = 1, SUM(LostParts)
      PartState(CurrentPartNum,1:6) = RecBuff(COUNTER+1:COUNTER+6)
      PDM%ParticleInside(CurrentPartNum) = .true.
      CALL SingleParticleToExactElement(CurrentPartNum)
      IF (PDM%ParticleInside(CurrentPartNum)) THEN
        PEM%LastElement(CurrentPartNum) = PEM%Element(CurrentPartNum)
        NbrOfFoundParts = NbrOfFoundParts + 1
        PartSpecies(CurrentPartNum) = INT(RecBuff(COUNTER+7))
        IF (useDSMC) THEN
          IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicState)) THEN
            PartStateIntEn(CurrentPartNum,1) = RecBuff(COUNTER+8)
            PartStateIntEn(CurrentPartNum,2) = RecBuff(COUNTER+9)
            PartStateIntEn(CurrentPartNum,3) = RecBuff(COUNTER+11)
            PartMPF(CurrentPartNum)          = RecBuff(COUNTER+10)
          ELSE IF ((CollisMode.GT.1).AND. (usevMPF)) THEN
            PartStateIntEn(CurrentPartNum,1) = RecBuff(COUNTER+8)
            PartStateIntEn(CurrentPartNum,2) = RecBuff(COUNTER+9)
            PartMPF(CurrentPartNum)          = RecBuff(COUNTER+10)
          ELSE IF ((CollisMode.GT.1).AND. (DSMC%ElectronicState)) THEN
            PartStateIntEn(CurrentPartNum,1) = RecBuff(COUNTER+8)
            PartStateIntEn(CurrentPartNum,2) = RecBuff(COUNTER+9)
            PartStateIntEn(CurrentPartNum,3) = RecBuff(COUNTER+10)
          ELSE IF (CollisMode.GT.1) THEN
            PartStateIntEn(CurrentPartNum,1) = RecBuff(COUNTER+8)
            PartStateIntEn(CurrentPartNum,2) = RecBuff(COUNTER+9)
          ELSE IF (usevMPF) THEN
            PartMPF(CurrentPartNum)          = RecBuff(COUNTER+8)
          END IF
        ELSE IF (usevMPF) THEN
          PartMPF(CurrentPartNum)          = RecBuff(COUNTER+8)
        END IF
        CurrentPartNum = CurrentPartNum + 1
      END IF
      COUNTER = COUNTER + PartDataSize
    END DO
    PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfFoundParts
    ! Combine number of found particles to make sure none are lost completely
    CALL MPI_ALLREDUCE(NbrOfFoundParts, CompleteNbrOfFound, 1, MPI_INTEGER, MPI_SUM, PartMPI%COMM, IERROR)
    SWRITE(UNIT_stdOut,*) SUM(LostParts),'were not in the correct proc after restart.'
    SWRITE(UNIT_stdOut,*) CompleteNbrOfFound,'of these were found in other procs.'
    SWRITE(UNIT_stdOut,*) SUM(LostParts)-CompleteNbrOfFound,'were not found and have been removed.'
  END IF
#else
  IF (COUNTER.NE.0) WRITE(*,*) COUNTER,'Particles are in different element after restart!'
  IF (COUNTER2.NE.0) WRITE(*,*) COUNTER2,'of which could not be found and are removed!'
#endif
  CALL UpdateNextFreePosition()

#endif /*PARTICLES*/

  CALL CloseDataFile() 
  ! Delete all files that will be rewritten
  CALL FlushHDF5(RestartTime)
  SWRITE(UNIT_stdOut,*)'Restart DONE!' 
ELSE
  ! Delete all files since we are doing a fresh start
  CALL FlushHDF5()
END IF
END SUBROUTINE Restart



SUBROUTINE FinalizeRestart()
!===================================================================================================================================
! Finalizes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Restart_Vars,ONLY:Vdm_GaussNRestart_GaussN,RestartInitIsDone
! IMPLICIT VARIABLE HANDLINGDGInitIsDone
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(Vdm_GaussNRestart_GaussN)
RestartInitIsDone = .FALSE.
END SUBROUTINE FinalizeRestart

END MODULE MOD_Restart
