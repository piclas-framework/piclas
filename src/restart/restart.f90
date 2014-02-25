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
USE MOD_HDF5_Input,ONLY:OpenDataFile,CloseDataFile,GetDataProps,ReadAttribute
USE MOD_ReadInTools,ONLY:GETLOGICAL,GETREALARRAY,ReadInDone 
USE MOD_DSMC_Vars,ONLY: UseDSMC
USE MOD_LD_Vars,ONLY: UseLD
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

! DSMC handling:
useDSMC=GETLOGICAL('UseDSMC','.FALSE.')
IF (useDSMC) THEN
  ReadInDone = .FALSE.
  maxNArgs = 3
ELSE
  maxNArgs = 2
END IF
! LD handling:
useLD=GETLOGICAL('UseLD','.FALSE.')
IF (useLD) THEN
  useDSMC=.TRUE.
  ReadInDone = .FALSE.
  maxNArgs = 3
END IF

! Check if we want to perform a restart
nArgs=COMMAND_ARGUMENT_COUNT()
IF (nArgs .EQ. maxNArgs) THEN
  ! Read in the state file we want to restart from
  CALL GETARG(maxNArgs,RestartFile)
  SWRITE(UNIT_StdOut,'(A,A,A)')' | Restarting from file "',TRIM(RestartFile),'":'
  DoRestart = .TRUE.
  CALL OpenDataFile(RestartFile,.FALSE.)
  CALL GetDataProps(nVar_Restart,N_Restart,nElems_Restart)
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
USE MOD_DG_Vars,ONLY:U
USE MOD_Mesh_Vars,ONLY:offsetElem
USE MOD_Restart_Vars,ONLY:Vdm_GaussNRestart_GaussN
USE MOD_Restart_Vars,ONLY:DoRestart,N_Restart,RestartFile,RestartTime,InterpolateSolution
USE MOD_Output_Vars,ONLY:ProjectName
USE MOD_ChangeBasis,ONLY:ChangeBasis3D
USE MOD_HDF5_input ,ONLY:OpenDataFile,CloseDataFile,File_ID,ReadArray,ReadAttribute
USE MOD_HDF5_Output,ONLY:FlushHDF5
#ifdef PP_POIS
USE MOD_Equation_Vars,ONLY:Phi
#endif
#ifdef PARTICLES
USE MOD_Particle_Vars, ONLY:PartState, PartSpecies, PEM, PDM, Species, nSpecies, usevMPF, PartMPF
USE MOD_part_tools, ONLY: UpdateNextFreePosition
USE MOD_DSMC_Vars,ONLY: UseDSMC, CollisMode,PartStateIntEn, DSMC
USE MOD_LD_Vars,       ONLY: UseLD, PartStateBulkValues
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
INTEGER                  :: iElem
#ifdef PARTICLES
INTEGER                  :: FirstElemInd,LastelemInd,iVar,i
INTEGER,ALLOCATABLE      :: PartInt(:,:)
INTEGER,PARAMETER        :: PartIntSize=2        !number of entries in each line of PartInt
INTEGER                  :: PartDataSize         !number of entries in each line of PartData
INTEGER                  :: locnPart,offsetnPart
INTEGER,PARAMETER        :: ELEM_FirstPartInd=1
INTEGER,PARAMETER        :: ELEM_LastPartInd=2
REAL,ALLOCATABLE         :: PartData(:,:)
#endif /*PARTICLES*/
!===================================================================================================================================
IF(DoRestart)THEN
SWRITE(UNIT_stdOut,*)'Restarting from File:',TRIM(RestartFile)
  CALL OpenDataFile(RestartFile,.FALSE.)
  ! Read in time from restart file
  CALL ReadAttribute(File_ID,'Time',1,RealScalar=RestartTime)
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
#endif
    SWRITE(UNIT_stdOut,*)'DONE!'
  END IF

#ifdef PARTICLES
  IF (useDSMC.AND.(.NOT.(useLD))) THEN
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
  ELSE IF (useLD) THEN
    IF ((CollisMode.GT.1).AND.(usevMPF) .AND. DSMC%ElectronicState ) THEN !int ener + 3, vmpf +1
      PartDataSize=16
    ELSE IF ((CollisMode.GT.1).AND.( (usevMPF) .OR. DSMC%ElectronicState ) ) THEN !int ener + 2 and vmpf + 1
                                                                             ! or int energ +3 but no vmpf +1
      PartDataSize=15
    ELSE IF (CollisMode.GT.1) THEN
      PartDataSize=14!int ener + 2
    ELSE IF (usevMPF) THEN
      PartDataSize=13!+ 1 vmpf
    ELSE
      PartDataSize=12 !+ 0
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
    IF (useDSMC.AND.(.NOT.(useLD))) THEN
      IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicState)) THEN
        PartStateIntEn(1:locnPart,1)=PartData(offsetnPart+1:offsetnPart+locnPart,8)
        PartStateIntEn(1:locnPart,2)=PartData(offsetnPart+1:offsetnPart+locnPart,9)
        PartMPF(1:locnPart)=PartData(offsetnPart+1:offsetnPart+locnPart,10)
        PartStateIntEn(1:locnPart,3)=PartData(offsetnPart+1:offsetnPart+locnPart,11)
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
    ELSE IF (useLD) THEN
      IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicState)) THEN
        PartStateIntEn(1:locnPart,1)=PartData(offsetnPart+1:offsetnPart+locnPart,8)
        PartStateIntEn(1:locnPart,2)=PartData(offsetnPart+1:offsetnPart+locnPart,9)
        PartMPF(1:locnPart)=PartData(offsetnPart+1:offsetnPart+locnPart,10)
        PartStateIntEn(1:locnPart,3)=PartData(offsetnPart+1:offsetnPart+locnPart,11)
        PartStateBulkValues(1:locnPart,1)=PartData(offsetnPart+1:offsetnPart+locnPart,12)
        PartStateBulkValues(1:locnPart,2)=PartData(offsetnPart+1:offsetnPart+locnPart,13)
        PartStateBulkValues(1:locnPart,3)=PartData(offsetnPart+1:offsetnPart+locnPart,14)
        PartStateBulkValues(1:locnPart,4)=PartData(offsetnPart+1:offsetnPart+locnPart,15)
        PartStateBulkValues(1:locnPart,5)=PartData(offsetnPart+1:offsetnPart+locnPart,16)
      ELSE IF ((CollisMode.GT.1).AND. (usevMPF)) THEN
        PartStateIntEn(1:locnPart,1)=PartData(offsetnPart+1:offsetnPart+locnPart,8)
        PartStateIntEn(1:locnPart,2)=PartData(offsetnPart+1:offsetnPart+locnPart,9)
        PartMPF(1:locnPart)=PartData(offsetnPart+1:offsetnPart+locnPart,10)
        PartStateBulkValues(1:locnPart,1)=PartData(offsetnPart+1:offsetnPart+locnPart,11)
        PartStateBulkValues(1:locnPart,2)=PartData(offsetnPart+1:offsetnPart+locnPart,12)
        PartStateBulkValues(1:locnPart,3)=PartData(offsetnPart+1:offsetnPart+locnPart,13)
        PartStateBulkValues(1:locnPart,4)=PartData(offsetnPart+1:offsetnPart+locnPart,14)
        PartStateBulkValues(1:locnPart,5)=PartData(offsetnPart+1:offsetnPart+locnPart,15)
      ELSE IF ((CollisMode.GT.1).AND. (DSMC%ElectronicState)) THEN
        PartStateIntEn(1:locnPart,1)=PartData(offsetnPart+1:offsetnPart+locnPart,8)
        PartStateIntEn(1:locnPart,2)=PartData(offsetnPart+1:offsetnPart+locnPart,9)
        PartStateIntEn(1:locnPart,3)=PartData(offsetnPart+1:offsetnPart+locnPart,10)
        PartStateBulkValues(1:locnPart,1)=PartData(offsetnPart+1:offsetnPart+locnPart,11)
        PartStateBulkValues(1:locnPart,2)=PartData(offsetnPart+1:offsetnPart+locnPart,12)
        PartStateBulkValues(1:locnPart,3)=PartData(offsetnPart+1:offsetnPart+locnPart,13)
        PartStateBulkValues(1:locnPart,4)=PartData(offsetnPart+1:offsetnPart+locnPart,14)
        PartStateBulkValues(1:locnPart,5)=PartData(offsetnPart+1:offsetnPart+locnPart,15)
      ELSE IF (CollisMode.GT.1) THEN
        PartStateIntEn(1:locnPart,1)=PartData(offsetnPart+1:offsetnPart+locnPart,8)
        PartStateIntEn(1:locnPart,2)=PartData(offsetnPart+1:offsetnPart+locnPart,9)
        PartStateBulkValues(1:locnPart,1)=PartData(offsetnPart+1:offsetnPart+locnPart,10)
        PartStateBulkValues(1:locnPart,2)=PartData(offsetnPart+1:offsetnPart+locnPart,11)
        PartStateBulkValues(1:locnPart,3)=PartData(offsetnPart+1:offsetnPart+locnPart,12)
        PartStateBulkValues(1:locnPart,4)=PartData(offsetnPart+1:offsetnPart+locnPart,13)
        PartStateBulkValues(1:locnPart,5)=PartData(offsetnPart+1:offsetnPart+locnPart,14)
      ELSE IF (usevMPF) THEN
        PartMPF(1:locnPart)=PartData(offsetnPart+1:offsetnPart+locnPart,8)
        PartStateBulkValues(1:locnPart,1)=PartData(offsetnPart+1:offsetnPart+locnPart,9)
        PartStateBulkValues(1:locnPart,2)=PartData(offsetnPart+1:offsetnPart+locnPart,10)
        PartStateBulkValues(1:locnPart,3)=PartData(offsetnPart+1:offsetnPart+locnPart,11)
        PartStateBulkValues(1:locnPart,4)=PartData(offsetnPart+1:offsetnPart+locnPart,12)
        PartStateBulkValues(1:locnPart,5)=PartData(offsetnPart+1:offsetnPart+locnPart,13)
      ELSE
        PartStateBulkValues(1:locnPart,1)=PartData(offsetnPart+1:offsetnPart+locnPart,8)
        PartStateBulkValues(1:locnPart,2)=PartData(offsetnPart+1:offsetnPart+locnPart,9)
        PartStateBulkValues(1:locnPart,3)=PartData(offsetnPart+1:offsetnPart+locnPart,10)
        PartStateBulkValues(1:locnPart,4)=PartData(offsetnPart+1:offsetnPart+locnPart,11)
        PartStateBulkValues(1:locnPart,5)=PartData(offsetnPart+1:offsetnPart+locnPart,12)  
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
    Species(i)%InsertedParticle = INT(Species(i)%ParticleEmission * RestartTime)
  END DO
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
