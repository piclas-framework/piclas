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
USE MOD_HDF5_Input,         ONLY:OpenDataFile,CloseDataFile,GetDataProps,ReadAttribute,File_ID
USE MOD_ReadInTools,        ONLY:GETLOGICAL,GETREALARRAY,ReadInDone 
#ifdef PARTICLES
USE MOD_DSMC_Vars,          ONLY: UseDSMC
USE MOD_LD_Vars,            ONLY: UseLD
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
   CALL abort(&
__STAMP__&
,'InitRestart not ready to be called or already called.',999,999.)
   RETURN
END IF

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT RESTART...'

#ifdef PARTICLES
! DSMC handling:
useDSMC=GETLOGICAL('UseDSMC','.FALSE.')
useLD=GETLOGICAL('UseLD','.FALSE.')
IF(useLD) useDSMC=.TRUE.
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
USE MOD_DG_Vars,                 ONLY:U
USE MOD_Mesh_Vars,               ONLY:offsetElem,DoWriteStateToHDF5
#ifdef PP_HDG
USE MOD_Mesh_Vars,               ONLY:offsetSide,nSides,nMPISides_YOUR, offsetSide
#else
USE MOD_Restart_Vars,            ONLY:Vdm_GaussNRestart_GaussN
#endif /*PP_HDG*/
USE MOD_Restart_Vars,            ONLY:DoRestart,N_Restart,RestartFile,RestartTime,InterpolateSolution
USE MOD_ChangeBasis,             ONLY:ChangeBasis3D
USE MOD_HDF5_input ,             ONLY:OpenDataFile,CloseDataFile,ReadArray,ReadAttribute
USE MOD_HDF5_Output,             ONLY:FlushHDF5
#ifndef PP_HDG
USE MOD_PML_Vars,                ONLY:DoPML,PMLToElem,U2,nPMLElems
#endif /*not PP_HDG*/
#ifdef PP_POIS
USE MOD_Equation_Vars,           ONLY:Phi
#endif /*PP_POIS*/
#ifdef PARTICLES
USE MOD_Particle_Vars,           ONLY:PartState, PartSpecies, PEM, PDM, Species, nSpecies, usevMPF, PartMPF,PartPosRef
USE MOD_part_tools,              ONLY: UpdateNextFreePosition
USE MOD_DSMC_Vars,               ONLY: UseDSMC, CollisMode,PartStateIntEn, DSMC
!USE MOD_BoundaryTools,           ONLY : SingleParticleToExactElement, ParticleInsideQuad3D
!USE MOD_BoundaryTools,           ONLY: ParticleInsideQuad3D
USE MOD_Eval_XYZ,                ONLY: EVal_xyz_ElemCheck
USE MOD_Particle_Mesh,           ONLY:SingleParticleToExactElement,SingleParticleToExactElementNoMap
USE MOD_Particle_Mesh_Vars,      ONLY:epsOneCell
USE MOD_Particle_Tracking_Vars,  ONLY:DoRefMapping
#ifdef MPI
USE MOD_Particle_MPI_Vars,       ONLY:PartMPI
#endif /*MPI*/
#ifdef PP_HDG
USE MOD_HDG_Vars,             ONLY: lambda, nGP_face
USE MOD_HDG,                  ONLY: RestartHDG
USE MOD_HDF5_Input,           ONLY: File_ID,DatasetExists
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
!USE MOD_TimeDisc,             ONLY: TimeStepPoissonByLSERK
#endif
#if (PP_TimeDiscMethod==500)
!USE MOD_TimeDisc,             ONLY:  TimeStepPoisson
#endif
USE MOD_TimeDisc_Vars,           ONLY: dt
#endif /*PP_HDG*/
#endif /*PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#ifndef PP_HDG
REAL,ALLOCATABLE         :: U_local(:,:,:,:,:)
REAL,ALLOCATABLE         :: U_local2(:,:,:,:,:)
INTEGER                  :: iPML
#else
LOGICAL                  :: DG_SolutionLambdaExists,DG_SolutionUExists
INTEGER(KIND=8)          :: iter
#endif /*not PP_HDG*/
INTEGER                  :: iElem
#ifdef MPI
REAL                     :: StartT,EndT
#endif /*MPI*/
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
INTEGER                  :: COUNTER, COUNTER2
#ifdef MPI
REAL, ALLOCATABLE        :: SendBuff(:), RecBuff(:)
INTEGER                  :: LostParts(0:PartMPI%nProcs-1), Displace(0:PartMPI%nProcs-1),CurrentPartNum
INTEGER                  :: NbrOfFoundParts, CompleteNbrOfFound, RecCount(0:PartMPI%nProcs-1)
#endif /*MPI*/
REAL                     :: VFR_total
#endif /*PARTICLES*/
!===================================================================================================================================
IF(DoRestart)THEN
SWRITE(UNIT_stdOut,*)'Restarting from File:',TRIM(RestartFile)
#ifdef MPI
  StartT=MPI_WTIME()
#endif
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
#elif defined PP_HDG
    CALL DatasetExists(File_ID,'DG_SolutionU',DG_SolutionUExists)
    IF(DG_SolutionUExists)THEN
      CALL ReadArray('DG_SolutionU',5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),OffsetElem,5,RealArray=U)
    ELSE
      CALL ReadArray('DG_Solution' ,5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),OffsetElem,5,RealArray=U)
    END IF
    CALL DatasetExists(File_ID,'DG_SolutionLambda',DG_SolutionLambdaExists)
    IF(DG_SolutionLambdaExists)THEN
      CALL ReadArray('DG_SolutionLambda',3,(/PP_nVar,nGP_face,nSides-nMPISides_YOUR/),offsetSide,3,RealArray=lambda)
      CALL RestartHDG(U) !
    ELSE
      lambda=0.
    END IF



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
#elif defined PP_HDG
CALL abort(&
__STAMP__&
,' Restart with changed polynomial degree not implemented for HDG!')
!    ALLOCATE(U_local(PP_nVar,0:N_Restart,0:N_Restart,0:N_Restart,PP_nElems))
!    CALL ReadArray('DG_SolutionLambda',5,(/PP_nVar,N_Restart+1,N_Restart+1,N_Restart+1,PP_nElems/),OffsetElem,5,RealArray=U_local)
!    DO iElem=1,PP_nElems
!      CALL ChangeBasis3D(PP_nVar,N_Restart,PP_N,Vdm_GaussNRestart_GaussN,U_local(:,:,:,:,iElem),U(:,:,:,:,iElem))
!    END DO
!    DEALLOCATE(U_local)
    !CALL RestartHDG(U)     
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
    SWRITE(UNIT_stdOut,*)' DONE!'
  END IF

#ifdef PARTICLES
  IF (useDSMC) THEN
    IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicModel)) THEN !int ener + 3, vmpf +1
      PartDataSize=11
    ELSE IF ((CollisMode.GT.1).AND.((usevMPF) .OR. (DSMC%ElectronicModel)) ) THEN ! int ener + 2 and vmpf +1
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
      IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicModel)) THEN
        PartStateIntEn(1:locnPart,1)=PartData(offsetnPart+1:offsetnPart+locnPart,8)
        PartStateIntEn(1:locnPart,2)=PartData(offsetnPart+1:offsetnPart+locnPart,9)
        PartStateIntEn(1:locnPart,3)=PartData(offsetnPart+1:offsetnPart+locnPart,10)
        PartMPF(1:locnPart)=PartData(offsetnPart+1:offsetnPart+locnPart,11)
      ELSE IF ((CollisMode.GT.1).AND. (usevMPF)) THEN
        PartStateIntEn(1:locnPart,1)=PartData(offsetnPart+1:offsetnPart+locnPart,8)
        PartStateIntEn(1:locnPart,2)=PartData(offsetnPart+1:offsetnPart+locnPart,9)
        PartMPF(1:locnPart)=PartData(offsetnPart+1:offsetnPart+locnPart,10)
      ELSE IF ((CollisMode.GT.1).AND. (DSMC%ElectronicModel)) THEN
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
  SWRITE(UNIT_stdOut,*)' DONE!' 
  DO i=1,nSpecies
    DO iInit = Species(i)%StartnumberOfInits, Species(i)%NumberOfInits
      Species(i)%Init(iInit)%InsertedParticle = INT(Species(i)%Init(iInit)%ParticleEmission * RestartTime,8)
    END DO
    DO iInit = 1, Species(i)%nSurfacefluxBCs
      IF (Species(i)%Surfaceflux(iInit)%ReduceNoise) THEN
        VFR_total = Species(i)%Surfaceflux(iInit)%VFR_total_allProcsTotal !proc global total (for non-root: dummy!!!)
      ELSE
        VFR_total = Species(i)%Surfaceflux(iInit)%VFR_total               !proc local total
      END IF
      Species(i)%Surfaceflux(iInit)%InsertedParticle = INT(Species(i)%Surfaceflux(iInit)%PartDensity * RestartTime &
        / Species(i)%MacroParticleFactor * VFR_total,8)
    END DO
  END DO
  ! if ParticleVecLength GT maxParticleNumber: Stop
  IF (PDM%ParticleVecLength.GT.PDM%maxParticleNumber) THEN
    CALL abort(&
__STAMP__&
,' Number of Particles in Restart higher than MaxParticleNumber!')
  END IF

  ! Since the elementside-local node number are NOT persistant and dependent on the location
  ! of the MPI borders, all particle-element mappings need to be checked after a restart
  ! Step 1: Identify particles that are not in the element in which they were before the restart
  COUNTER = 0
  COUNTER2 = 0
  IF(DoRefMapping) THEN
    DO i = 1,PDM%ParticleVecLength
      CALL Eval_xyz_ElemCheck(PartState(i,1:3),Xi,PEM%Element(i))
      IF(ALL(ABS(Xi).LE.EpsOneCell(PEM%Element(i)))) THEN ! particle inside
        InElementCheck=.TRUE.
        PartPosRef(1:3,i)=Xi
      ELSE
        InElementCheck=.FALSE.
      END IF
      IF (.NOT.InElementCheck) THEN  ! try to find them within MyProc
        COUNTER = COUNTER + 1
        !CALL SingleParticleToExactElement(i)
        CALL SingleParticleToExactElement(i,doHALO=.FALSE.)
        IF (.NOT.PDM%ParticleInside(i)) THEN
          COUNTER2 = COUNTER2 + 1
          PartPosRef(1:3,i) = -888.
        ELSE
          PEM%LastElement(i) = PEM%Element(i)
        END IF
      END IF
    END DO
  ELSE ! no Ref Mapping
    DO i = 1,PDM%ParticleVecLength
      CALL Eval_xyz_ElemCheck(PartState(i,1:3),Xi,PEM%Element(i))
      IF(ALL(ABS(Xi).LE.1.0)) THEN ! particle inside
        InElementCheck=.TRUE.
        IF(ALLOCATED(PartPosRef)) PartPosRef(1:3,i)=Xi
      ELSE
        InElementCheck=.FALSE.
      END IF
      IF (.NOT.InElementCheck) THEN  ! try to find them within MyProc
        COUNTER = COUNTER + 1
        !CALL SingleParticleToExactElement(i)
        CALL SingleParticleToExactElementNoMap(i,doHALO=.FALSE.)
        IF (.NOT.PDM%ParticleInside(i)) THEN
          COUNTER2 = COUNTER2 + 1
        ELSE
          PEM%LastElement(i) = PEM%Element(i)
        END IF
      END IF
    END DO
  END IF
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
          IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicModel)) THEN
            SendBuff(COUNTER+8)  = PartStateIntEn(i,1)
            SendBuff(COUNTER+9)  = PartStateIntEn(i,2)
            SendBuff(COUNTER+10) = PartMPF(i)
            SendBuff(COUNTER+11) = PartStateIntEn(i,3)
          ELSE IF ((CollisMode.GT.1).AND. (usevMPF)) THEN
            SendBuff(COUNTER+8)  = PartStateIntEn(i,1)
            SendBuff(COUNTER+9)  = PartStateIntEn(i,2)
            SendBuff(COUNTER+10) = PartMPF(i)
          ELSE IF ((CollisMode.GT.1).AND. (DSMC%ElectronicModel)) THEN
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
      IF(DoRefMapping)THEN
        CALL SingleParticleToExactElement(CurrentPartNum,doHALO=.FALSE.)
      ELSE
        CALL SingleParticleToExactElementNoMap(CurrentPartNum,doHALO=.FALSE.)
      END IF
      !CALL SingleParticleToExactElement(CurrentPartNum)
      IF (PDM%ParticleInside(CurrentPartNum)) THEN
        PEM%LastElement(CurrentPartNum) = PEM%Element(CurrentPartNum)
        NbrOfFoundParts = NbrOfFoundParts + 1
        PartSpecies(CurrentPartNum) = INT(RecBuff(COUNTER+7))
        IF (useDSMC) THEN
          IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicModel)) THEN
            PartStateIntEn(CurrentPartNum,1) = RecBuff(COUNTER+8)
            PartStateIntEn(CurrentPartNum,2) = RecBuff(COUNTER+9)
            PartStateIntEn(CurrentPartNum,3) = RecBuff(COUNTER+11)
            PartMPF(CurrentPartNum)          = RecBuff(COUNTER+10)
          ELSE IF ((CollisMode.GT.1).AND. (usevMPF)) THEN
            PartStateIntEn(CurrentPartNum,1) = RecBuff(COUNTER+8)
            PartStateIntEn(CurrentPartNum,2) = RecBuff(COUNTER+9)
            PartMPF(CurrentPartNum)          = RecBuff(COUNTER+10)
          ELSE IF ((CollisMode.GT.1).AND. (DSMC%ElectronicModel)) THEN
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



#ifdef PP_HDG
!print*,RestartTime
iter=0
!print*,iter
IF(RestartTime.GT.0.0)THEN
  dt=RestartTime/1e6
ELSE
  dt=1e-19
END IF
!print*,dt
#if (PP_TimeDiscMethod==500)
    CALL RestartTimeStepPoisson(RestartTime) ! Euler Explicit, Poisson
#else
    CALL RestartTimeStepPoissonByLSERK(RestartTime,iter,0.)  !Runge Kutta Explicit, Poisson
#endif
#endif /*PP_HDG*/




  ! Delete all files that will be rewritten
  IF(DoWriteStateToHDF5) CALL FlushHDF5(RestartTime)
#ifdef MPI
  EndT=MPI_WTIME()
  SWRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')' Restart took  [',EndT-StartT,'s] for readin.'
  SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' Restart DONE!' 
#else
  SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' Restart DONE!' 
#endif
ELSE
  ! Delete all files since we are doing a fresh start
  IF(DoWriteStateToHDF5) CALL FlushHDF5()
END IF
END SUBROUTINE Restart


#ifdef PP_HDG
#if (PP_TimeDiscMethod==500)
SUBROUTINE RestartTimeStepPoisson(t)
!===================================================================================================================================
! Hesthaven book, page 64
! Low-Storage Runge-Kutta integration of degree 4 with 5 stages.
! This procedure takes the current time t, the time step dt and the solution at
! the current time U(t) and returns the solution at the next time level.
!===================================================================================================================================
! MODULES
USE MOD_DG_Vars,                 ONLY: U
USE MOD_PreProc
USE MOD_TimeDisc_Vars,           ONLY: dt,iter
USE MOD_HDG,                     ONLY: HDG
USE MOD_Particle_Tracking_vars,  ONLY: DoRefMapping!,MeasureTrackTime
#ifdef PARTICLES
USE MOD_PICDepo,                 ONLY: Deposition
USE MOD_PICInterpolation,        ONLY: InterpolateFieldToParticle
USE MOD_Particle_Vars,           ONLY: PartState, Pt, LastPartPos,PEM, PDM, usevMPF, doParticleMerge, DelayTime, PartPressureCell
!USE MOD_Particle_Vars,           ONLY : Time
USE MOD_Particle_Vars,           ONLY: DoSurfaceFlux
USE MOD_part_RHS,                ONLY: CalcPartRHS
!USE MOD_part_boundary,           ONLY : ParticleBoundary
USE MOD_part_emission,           ONLY: ParticleInserting, ParticleSurfaceflux
USE MOD_DSMC,                    ONLY: DSMC_main
USE MOD_DSMC_Vars,               ONLY: useDSMC, DSMC_RHS
USE MOD_part_MPFtools,           ONLY: StartParticleMerge
USE MOD_PIC_Analyze,             ONLY: VerifyDepositedCharge
USE MOD_Particle_Analyze_Vars,   ONLY: DoVerifyCharge
#ifdef MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY: PartMPIExchange
#endif
!USE MOD_PIC_Analyze,      ONLY: CalcDepositedCharge
USE MOD_part_tools,              ONLY: UpdateNextFreePosition
USE MOD_Particle_Tracking_vars,  ONLY: tTracking,tLocalization,DoRefMapping!,MeasureTrackTime
USE MOD_Particle_Tracking,       ONLY: ParticleTrackingCurved,ParticleRefTracking
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)       :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iPart
REAL    :: RandVal, dtFrac
!===================================================================================================================================
!Time = t

IF ((t.GE.DelayTime).OR.(iter.EQ.0)) THEN
  CALL Deposition(doInnerParts=.TRUE.) ! because of emmision and UpdateParticlePosition
#ifdef MPI
  ! here: finish deposition with delta kernal
  !       maps source terms in physical space
  ! ALWAYS require
  PartMPIExchange%nMPIParticles=0
#endif /*MPI*/
  CALL Deposition(doInnerParts=.FALSE.) ! needed for closing communication
  IF(DoVerifyCharge) CALL VerifyDepositedCharge()
END IF

! EM field
CALL HDG(t,U,iter)

IF (t.GE.DelayTime) THEN
  CALL InterpolateFieldToParticle(doInnerParts=.TRUE.)
  !CALL InterpolateFieldToParticle(doInnerParts=.FALSE.) ! only needed when MPI communation changes the number of parts
  CALL CalcPartRHS()
END IF

IF (DoSurfaceFlux) THEN
  LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
  LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
  LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
  PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
  IF (t.GE.DelayTime) THEN
    CALL ParticleSurfaceflux() !dtFracPush (SurfFlux): LastPartPos and LastElem already set!
  END IF
  IF (t.GE.DelayTime) THEN ! Euler-Explicit only for Particles
    DO iPart=1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        IF (.NOT.PDM%dtFracPush(iPart)) THEN
          PartState(iPart,1) = PartState(iPart,1) + PartState(iPart,4) * dt
          PartState(iPart,2) = PartState(iPart,2) + PartState(iPart,5) * dt
          PartState(iPart,3) = PartState(iPart,3) + PartState(iPart,6) * dt
          PartState(iPart,4) = PartState(iPart,4) + Pt(iPart,1) * dt
          PartState(iPart,5) = PartState(iPart,5) + Pt(iPart,2) * dt
          PartState(iPart,6) = PartState(iPart,6) + Pt(iPart,3) * dt
        ELSE
          CALL RANDOM_NUMBER(RandVal)
          dtFrac = dt * RandVal
          PartState(iPart,1) = PartState(iPart,1) + PartState(iPart,4) * dtFrac
          PartState(iPart,2) = PartState(iPart,2) + PartState(iPart,5) * dtFrac
          PartState(iPart,3) = PartState(iPart,3) + PartState(iPart,6) * dtFrac
          PDM%dtFracPush(iPart) = .FALSE.
        END IF
      END IF
    END DO
  END IF
ELSE
  LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
  LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
  LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
  PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
  IF (t.GE.DelayTime) THEN ! Euler-Explicit only for Particles
    PartState(1:PDM%ParticleVecLength,1) = PartState(1:PDM%ParticleVecLength,1) + PartState(1:PDM%ParticleVecLength,4) * dt
    PartState(1:PDM%ParticleVecLength,2) = PartState(1:PDM%ParticleVecLength,2) + PartState(1:PDM%ParticleVecLength,5) * dt
    PartState(1:PDM%ParticleVecLength,3) = PartState(1:PDM%ParticleVecLength,3) + PartState(1:PDM%ParticleVecLength,6) * dt
    PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) + Pt(1:PDM%ParticleVecLength,1) * dt
    PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) + Pt(1:PDM%ParticleVecLength,2) * dt
    PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) + Pt(1:PDM%ParticleVecLength,3) * dt
  END IF
END IF

IF ((t.GE.DelayTime).OR.(iter.EQ.0)) THEN
#ifdef MPI
  CALL IRecvNbofParticles() ! open receive buffer for number of particles
#endif
  IF(DoRefMapping)THEN
    CALL ParticleRefTracking()
  ELSE
    CALL ParticleTrackingCurved()
  END IF
#ifdef MPI
  CALL SendNbOfParticles() ! send number of particles
  CALL MPIParticleSend()  ! finish communication of number of particles and send particles
  CALL MPIParticleRecv()  ! finish communication
#endif
END IF

IF (t.GE.DelayTime) CALL ParticleInserting()

IF (doParticleMerge) THEN
  IF (.NOT.(useDSMC.OR.PartPressureCell)) THEN
    ALLOCATE(PEM%pStart(1:PP_nElems)           , &
             PEM%pNumber(1:PP_nElems)          , &
             PEM%pNext(1:PDM%maxParticleNumber), &
             PEM%pEnd(1:PP_nElems) )
  END IF
END IF

IF ((t.GE.DelayTime).OR.(iter.EQ.0)) THEN
  CALL UpdateNextFreePosition()
END IF

IF (doParticleMerge) THEN
  CALL StartParticleMerge()  
  IF (.NOT.(useDSMC.OR.PartPressureCell)) THEN
    DEALLOCATE(PEM%pStart , &
               PEM%pNumber, &
               PEM%pNext  , &
               PEM%pEnd   )
  END IF
  CALL UpdateNextFreePosition()
END IF

IF (useDSMC) THEN
  IF (t.GE.DelayTime) THEN
    CALL DSMC_main()
    PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) &
                                           + DSMC_RHS(1:PDM%ParticleVecLength,1)
    PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) &
                                           + DSMC_RHS(1:PDM%ParticleVecLength,2)
    PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) &
                                           + DSMC_RHS(1:PDM%ParticleVecLength,3)
  END IF
END IF
END SUBROUTINE RestartTimeStepPoisson
#endif /*(PP_TimeDiscMethod==500)*/

#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
SUBROUTINE RestartTimeStepPoissonByLSERK(t,iter,tEndDiff)
!===================================================================================================================================
! Hesthaven book, page 64
! Low-Storage Runge-Kutta integration of degree 4 with 5 stages.
! This procedure takes the current time t, the time step dt and the solution at
! the current time U(t) and returns the solution at the next time level.
!===================================================================================================================================
! MODULES
USE MOD_Globals,               ONLY: Abort
USE MOD_PreProc
USE MOD_Analyze,               ONLY: PerformAnalyze
USE MOD_TimeDisc_Vars,         ONLY: dt,iStage,RKdtFrac,RKdtFracTotal
USE MOD_TimeDisc_Vars,         ONLY: RK_a,RK_b,RK_c,nRKStages
USE MOD_DG_Vars,               ONLY: U
USE MOD_Particle_Tracking_vars,ONLY: DoRefMapping
#ifdef PARTICLES
USE MOD_PICDepo,               ONLY: Deposition
USE MOD_PICInterpolation,      ONLY: InterpolateFieldToParticle
USE MOD_Particle_Vars,         ONLY: PartState, Pt, Pt_temp, LastPartPos, DelayTime,  PEM, PDM, & 
                                     doParticleMerge,PartPressureCell,DoSurfaceFlux
USE MOD_part_RHS,              ONLY: CalcPartRHS
USE MOD_part_emission,         ONLY: ParticleInserting, ParticleSurfaceflux
USE MOD_DSMC,                  ONLY: DSMC_main
USE MOD_DSMC_Vars,             ONLY: useDSMC, DSMC_RHS
USE MOD_part_MPFtools,         ONLY: StartParticleMerge
USE MOD_PIC_Analyze,           ONLY: VerifyDepositedCharge
USE MOD_Particle_Analyze_Vars,   ONLY: DoVerifyCharge
#ifdef MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY: PartMPIExchange
#endif
USE MOD_Particle_Tracking_vars, ONLY: DoRefMapping
USE MOD_part_tools,             ONLY: UpdateNextFreePosition
USE MOD_Particle_Tracking,      ONLY: ParticleTrackingCurved,ParticleRefTracking
#endif /*PARTICLES*/
USE MOD_HDG           ,ONLY: HDG
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: tendDiff
REAL,INTENT(INOUT)            :: t
INTEGER(KIND=8),INTENT(INOUT) :: iter
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: tStage,b_dt(1:nRKStages),RK_a_rebuilt(1:nRKStages)
INTEGER                       :: iPart              
REAL                          :: RandVal, dtFrac
!===================================================================================================================================

DO iStage=1,nRKStages
  ! RK coefficients
  b_dt(iStage)=RK_b(iStage)*dt
  ! Rebuild Pt_tmp(1:3)-coefficients assuming F=0 and const. v in previous stages: Pt_tmp(i)=(/v(i)*RK_a_rebuilt(i),a(i)/)
  IF (iStage.EQ.1) THEN
    RK_a_rebuilt(:)=1.
  ELSE
    RK_a_rebuilt(iStage) = 1. - RK_a(iStage)*RK_a_rebuilt(iStage-1)
  END IF
END DO
iStage=1

! first RK step
#ifdef PARTICLES
!Time=t
RKdtFrac = RK_c(2)
RKdtFracTotal=RKdtFrac

IF ((t.GE.DelayTime).OR.(iter.EQ.0)) THEN
  CALL Deposition(doInnerParts=.TRUE.) ! because of emmision and UpdateParticlePosition
#ifdef MPI
  ! here: finish deposition with delta kernal
  !       maps source terms in physical space
  ! ALWAYS require
  PartMPIExchange%nMPIParticles=0
#endif /*MPI*/
  CALL Deposition(doInnerParts=.FALSE.) ! needed for closing communication
  IF(DoVerifyCharge) CALL VerifyDepositedCharge()
END IF
#endif /*PARTICLES*/

CALL HDG(t,U,iter)

#ifdef PARTICLES
IF (t.GE.DelayTime) THEN
  CALL InterpolateFieldToParticle(doInnerParts=.TRUE.)   ! forces on particles
  !CALL InterpolateFieldToParticle(doInnerParts=.FALSE.) ! only needed when MPI communation changes the number of parts
  CALL CalcPartRHS()
END IF
#endif /*PARTICLES*/

! calling the analyze routines
CALL PerformAnalyze(t,iter,tendDiff,forceAnalyze=.FALSE.,OutPut=.FALSE.)

#ifdef PARTICLES
! particles
LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
IF (DoSurfaceFlux .AND. (t.GE.DelayTime)) THEN
  CALL ParticleSurfaceflux() !dtFracPush (SurfFlux): LastPartPos and LastElem already set!
END IF
IF (t.GE.DelayTime) THEN
  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
      !-- Pt is not known only for new Surfaceflux-Parts -> change IsNewPart back to F for other Parts
      IF (.NOT.DoSurfaceFlux) THEN
        PDM%IsNewPart(iPart)=.FALSE.
      ELSE
        IF (.NOT.PDM%dtFracPush(iPart)) PDM%IsNewPart(iPart)=.FALSE.
      END IF
      !-- Particle Push
      IF (.NOT.PDM%IsNewPart(iPart)) THEN
        Pt_temp(iPart,1) = PartState(iPart,4)
        Pt_temp(iPart,2) = PartState(iPart,5)
        Pt_temp(iPart,3) = PartState(iPart,6)
        Pt_temp(iPart,4) = Pt(iPart,1)
        Pt_temp(iPart,5) = Pt(iPart,2)
        Pt_temp(iPart,6) = Pt(iPart,3)
        PartState(iPart,1) = PartState(iPart,1) + PartState(iPart,4)*b_dt(1)
        PartState(iPart,2) = PartState(iPart,2) + PartState(iPart,5)*b_dt(1)
        PartState(iPart,3) = PartState(iPart,3) + PartState(iPart,6)*b_dt(1)
        PartState(iPart,4) = PartState(iPart,4) + Pt(iPart,1)*b_dt(1)
        PartState(iPart,5) = PartState(iPart,5) + Pt(iPart,2)*b_dt(1)
        PartState(iPart,6) = PartState(iPart,6) + Pt(iPart,3)*b_dt(1)
      ELSE !IsNewPart
        IF (DoSurfaceFlux .AND. PDM%dtFracPush(iPart)) THEN !SF, new in current RKStage (no forces assumed in this stage)
          CALL RANDOM_NUMBER(RandVal)
          dtFrac = dt * RKdtFrac * RandVal
          PDM%dtFracPush(iPart) = .FALSE.
        ELSE
          CALL abort(&
          __STAMP__&
          ,'Error in LSERK-HDG-Timedisc: This case should be impossible...')
        END IF
        PartState(iPart,1) = PartState(iPart,1) + PartState(iPart,4) * dtFrac
        PartState(iPart,2) = PartState(iPart,2) + PartState(iPart,5) * dtFrac
        PartState(iPart,3) = PartState(iPart,3) + PartState(iPart,6) * dtFrac
        !!!PDM%IsNewPart(iPart) = .FALSE. !do not(!) change to false: flag needed for next RKStages
      END IF
    END IF
  END DO
END IF

IF ((t.GE.DelayTime).OR.(iter.EQ.0)) THEN
#ifdef MPI
  CALL IRecvNbofParticles() ! open receive buffer for number of particles
#endif
  IF(DoRefMapping)THEN
    CALL ParticleRefTracking()
  ELSE
    CALL ParticleTrackingCurved()
  END IF
#ifdef MPI
  CALL SendNbOfParticles() ! send number of particles
  CALL MPIParticleSend()   ! finish communication of number of particles and send particles
  CALL MPIParticleRecv()   ! finish communication
#endif
END IF

IF (t.GE.DelayTime) CALL ParticleInserting()
#endif /*PARTICLES*/

! perform RK steps
DO iStage=2,nRKStages
  tStage=t+dt*RK_c(iStage)
#ifdef PARTICLES
  IF (iStage.NE.nRKStages) THEN
    RKdtFrac = RK_c(iStage+1)-RK_c(iStage)
    RKdtFracTotal=RKdtFracTotal+RKdtFrac
  ELSE
    RKdtFrac = 1.-RK_c(nRKStages)
    RKdtFracTotal=1.
  END IF

  ! deposition 
  IF (t.GE.DelayTime) THEN
    CALL Deposition(doInnerParts=.TRUE.) ! because of emmision and UpdateParticlePosition
#ifdef MPI
    ! here: finish deposition with delta kernal
    !       maps source terms in physical space
    ! ALWAYS require
    PartMPIExchange%nMPIParticles=0
#endif /*MPI*/
    CALL Deposition(doInnerParts=.FALSE.) ! needed for closing communication
    IF(DoVerifyCharge) CALL VerifyDepositedCharge()
  END IF
#endif /*PARTICLES*/

  CALL HDG(tStage,U,iter)

#ifdef PARTICLES
  IF (t.GE.DelayTime) THEN
    ! forces on particle
    CALL InterpolateFieldToParticle(doInnerParts=.TRUE.)   ! forces on particles
    !CALL InterpolateFieldToParticle(doInnerParts=.FALSE.) ! only needed when MPI communation changes the number of parts
    CALL CalcPartRHS()

    ! particle step
    LastPartPos(1:PDM%ParticleVecLength,1)=PartState(1:PDM%ParticleVecLength,1)
    LastPartPos(1:PDM%ParticleVecLength,2)=PartState(1:PDM%ParticleVecLength,2)
    LastPartPos(1:PDM%ParticleVecLength,3)=PartState(1:PDM%ParticleVecLength,3)
    PEM%lastElement(1:PDM%ParticleVecLength)=PEM%Element(1:PDM%ParticleVecLength)
    IF (DoSurfaceFlux) CALL ParticleSurfaceflux() !dtFracPush (SurfFlux): LastPartPos and LastElem already set!
    DO iPart=1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        IF (.NOT.PDM%IsNewPart(iPart)) THEN
          Pt_temp(iPart,1) = PartState(iPart,4) - RK_a(iStage) * Pt_temp(iPart,1)
          Pt_temp(iPart,2) = PartState(iPart,5) - RK_a(iStage) * Pt_temp(iPart,2)
          Pt_temp(iPart,3) = PartState(iPart,6) - RK_a(iStage) * Pt_temp(iPart,3)
          Pt_temp(iPart,4) = Pt(iPart,1) - RK_a(iStage) * Pt_temp(iPart,4)
          Pt_temp(iPart,5) = Pt(iPart,2) - RK_a(iStage) * Pt_temp(iPart,5)
          Pt_temp(iPart,6) = Pt(iPart,3) - RK_a(iStage) * Pt_temp(iPart,6)
          PartState(iPart,1) = PartState(iPart,1) + Pt_temp(iPart,1)*b_dt(iStage)
          PartState(iPart,2) = PartState(iPart,2) + Pt_temp(iPart,2)*b_dt(iStage)
          PartState(iPart,3) = PartState(iPart,3) + Pt_temp(iPart,3)*b_dt(iStage)
          PartState(iPart,4) = PartState(iPart,4) + Pt_temp(iPart,4)*b_dt(iStage)
          PartState(iPart,5) = PartState(iPart,5) + Pt_temp(iPart,5)*b_dt(iStage)
          PartState(iPart,6) = PartState(iPart,6) + Pt_temp(iPart,6)*b_dt(iStage)
        ELSE !IsNewPart: no Pt_temp history available!
          IF (DoSurfaceFlux .AND. PDM%dtFracPush(iPart)) THEN !SF, new in current RKStage (no forces assumed in this stage)
            CALL RANDOM_NUMBER(RandVal)
            dtFrac = dt * RKdtFrac * RandVal
            PartState(iPart,1) = PartState(iPart,1) + PartState(iPart,4) * dtFrac
            PartState(iPart,2) = PartState(iPart,2) + PartState(iPart,5) * dtFrac
            PartState(iPart,3) = PartState(iPart,3) + PartState(iPart,6) * dtFrac
            PDM%dtFracPush(iPart) = .FALSE.
            !!!PDM%IsNewPart(iPart) = .FALSE. !do not(!) change to false: flag needed for next RKStages
          ELSE !new but without SF in current RKStage (i.e., from ParticleInserting or diffusive wall reflection)
               ! -> assuming F=0 and const. v in previous stages with RK_a_rebuilt (see above)
            Pt_temp(iPart,1) = PartState(iPart,4) * RK_a_rebuilt(iStage)
            Pt_temp(iPart,2) = PartState(iPart,5) * RK_a_rebuilt(iStage)
            Pt_temp(iPart,3) = PartState(iPart,6) * RK_a_rebuilt(iStage)
            Pt_temp(iPart,4) = Pt(iPart,1)
            Pt_temp(iPart,5) = Pt(iPart,2)
            Pt_temp(iPart,6) = Pt(iPart,3)
            PartState(iPart,1) = PartState(iPart,1) + Pt_temp(iPart,1)*b_dt(iStage)
            PartState(iPart,2) = PartState(iPart,2) + Pt_temp(iPart,2)*b_dt(iStage)
            PartState(iPart,3) = PartState(iPart,3) + Pt_temp(iPart,3)*b_dt(iStage)
            PartState(iPart,4) = PartState(iPart,4) + Pt_temp(iPart,4)*b_dt(iStage)
            PartState(iPart,5) = PartState(iPart,5) + Pt_temp(iPart,5)*b_dt(iStage)
            PartState(iPart,6) = PartState(iPart,6) + Pt_temp(iPart,6)*b_dt(iStage)
            PDM%IsNewPart(iPart) = .FALSE. !"normal" part in next iter
          END IF
        END IF
      END IF
    END DO

    ! particle tracking
#ifdef MPI
    CALL IRecvNbofParticles() ! open receive buffer for number of particles
#endif
    IF(DoRefMapping)THEN
      CALL ParticleRefTracking()
    ELSE
      CALL ParticleTrackingCurved()
    END IF
#ifdef MPI
    CALL SendNbOfParticles() ! send number of particles
    CALL MPIParticleSend()   ! finish communication of number of particles and send particles
    CALL MPIParticleRecv()   ! finish communication
#endif
    CALL ParticleInserting()
  END IF
#endif /*PARTICLES*/
END DO

#ifdef PARTICLES
IF (doParticleMerge) THEN
  IF (.NOT.(useDSMC.OR.PartPressureCell)) THEN
    ALLOCATE(PEM%pStart(1:PP_nElems)           , &
             PEM%pNumber(1:PP_nElems)          , &
             PEM%pNext(1:PDM%maxParticleNumber), &
             PEM%pEnd(1:PP_nElems) )
  END IF
END IF

IF ((t.GE.DelayTime).OR.(iter.EQ.0)) THEN
  CALL UpdateNextFreePosition()
END IF

IF (doParticleMerge) THEN
  CALL StartParticleMerge()  
  IF (.NOT.(useDSMC.OR.PartPressureCell)) THEN
    DEALLOCATE(PEM%pStart , &
               PEM%pNumber, &
               PEM%pNext  , &
               PEM%pEnd   )
  END IF
  CALL UpdateNextFreePosition()
END IF

IF (useDSMC) THEN
  IF (t.GE.DelayTime) THEN
    CALL DSMC_main()
    PartState(1:PDM%ParticleVecLength,4) = PartState(1:PDM%ParticleVecLength,4) &
                                           + DSMC_RHS(1:PDM%ParticleVecLength,1)
    PartState(1:PDM%ParticleVecLength,5) = PartState(1:PDM%ParticleVecLength,5) &
                                           + DSMC_RHS(1:PDM%ParticleVecLength,2)
    PartState(1:PDM%ParticleVecLength,6) = PartState(1:PDM%ParticleVecLength,6) &
                                           + DSMC_RHS(1:PDM%ParticleVecLength,3)
  END IF
END IF
#endif /*PARTICLES*/

END SUBROUTINE RestartTimeStepPoissonByLSERK
#endif /*(PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)*/
#endif /*PP_HDG*/

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
