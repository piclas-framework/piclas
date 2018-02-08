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

PUBLIC :: DefineParametersRestart
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters.
!==================================================================================================================================
SUBROUTINE DefineParametersRestart()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Restart")
CALL prms%CreateLogicalOption('ResetTime', "Override solution time to t=0 on restart.", '.FALSE.')
END SUBROUTINE DefineParametersRestart


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
USE MOD_ReadInTools,        ONLY:GETLOGICAL,GETREALARRAY
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
  !ReadInDone = .FALSE.
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
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
#else
  CALL OpenDataFile(RestartFile,create=.FALSE.,readOnly=.TRUE.)
#endif
  CALL GetDataProps(nVar_Restart,N_Restart,nElems_Restart,NodeType_Restart)
  ! Read in time from restart file
  CALL ReadAttribute(File_ID,'Time',1,RealScalar=RestartTime)
  CALL CloseDataFile() 
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
USE MOD_HDF5_Input,              ONLY:DatasetExists
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
USE MOD_PML_Vars,                ONLY:DoPML,PMLToElem,U2,nPMLElems,PMLnVar
#endif /*not PP_HDG*/
#ifdef PP_POIS
USE MOD_Equation_Vars,           ONLY:Phi
#endif /*PP_POIS*/
#ifdef PARTICLES
USE MOD_Particle_Vars,           ONLY:PartState, PartSpecies, PEM, PDM, Species, nSpecies, usevMPF, PartMPF,PartPosRef
USE MOD_part_tools,              ONLY:UpdateNextFreePosition
USE MOD_DSMC_Vars,               ONLY:UseDSMC, CollisMode,PartStateIntEn, DSMC
USE MOD_Eval_XYZ,                ONLY:EVal_xyz_ElemCheck
USE MOD_Particle_Mesh,           ONLY:SingleParticleToExactElement,SingleParticleToExactElementNoMap
USE MOD_Particle_Mesh_Vars,      ONLY:epsOneCell
USE MOD_Particle_Tracking_Vars,  ONLY:DoRefMapping
#ifdef MPI
USE MOD_Particle_MPI_Vars,       ONLY:PartMPI
#endif /*MPI*/
USE MOD_Particle_Tracking,       ONLY:ParticleCollectCharges
USE MOD_TTM_Vars,                ONLY:DoImportTTMFile,TTM_DG
#endif /*PARTICLES*/
#ifdef PP_HDG
USE MOD_HDG_Vars,                ONLY:lambda, nGP_face
USE MOD_HDG,                     ONLY:RestartHDG
#endif /*PP_HDG*/
#if USE_QDS_DG
USE MOD_QDS_DG_Vars,             ONLY:DoQDS,QDSMacroValues,nQDSElems,QDSSpeciesMass
#endif /*USE_QDS_DG*/
USE MOD_HDF5_Input,              ONLY:File_ID,DatasetExists
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
INTEGER                  :: FirstElemInd,LastelemInd,iInit
INTEGER,ALLOCATABLE      :: PartInt(:,:)
INTEGER,PARAMETER        :: PartIntSize=2        !number of entries in each line of PartInt
INTEGER                  :: PartDataSize         !number of entries in each line of PartData
INTEGER                  :: locnPart,offsetnPart
INTEGER,PARAMETER        :: ELEM_FirstPartInd=1
INTEGER,PARAMETER        :: ELEM_LastPartInd=2
REAL,ALLOCATABLE         :: PartData(:,:)
REAL                     :: xi(3)
LOGICAL                  :: InElementCheck,PartIntExists
INTEGER                  :: COUNTER, COUNTER2
#ifdef MPI
REAL, ALLOCATABLE        :: SendBuff(:), RecBuff(:)
INTEGER                  :: LostParts(0:PartMPI%nProcs-1), Displace(0:PartMPI%nProcs-1),CurrentPartNum
INTEGER                  :: NbrOfFoundParts, CompleteNbrOfFound, RecCount(0:PartMPI%nProcs-1)
#endif /*MPI*/
REAL                     :: VFR_total
CHARACTER(255)           :: TTMRestartFile !> TTM Data file for restart
LOGICAL                  :: TTM_DG_SolutionExists
#endif /*PARTICLES*/
INTEGER                  :: IndNum         !> auxiliary variable containing the index number of a substring within a string
#if USE_QDS_DG
CHARACTER(255)           :: QDSRestartFile !> QDS Data file for restart
LOGICAL                  :: QDS_DG_SolutionExists
INTEGER                  :: j,k
#endif /*USE_QDS_DG*/
#if (USE_QDS_DG) || (PARTICLES)
INTEGER                  :: i
#endif
!===================================================================================================================================
IF(DoRestart)THEN
  SWRITE(UNIT_stdOut,*)'Restarting from File:',TRIM(RestartFile)
#ifdef MPI
  StartT=MPI_WTIME()
#endif

#ifdef PARTICLES
  IF(DoImportTTMFile)THEN ! read TTM data from "XXXXX_TTM_000.0XXXXXXXXXXX.h5"
    IF(.NOT. InterpolateSolution)THEN! No interpolation needed, read solution directly from file
      IndNum=INDEX(RestartFile,'State')
      IF(IndNum.LE.0)CALL abort(&
        __STAMP__&
        ,' Restart file does not contain "State" character string?! Supply restart file for reading TTM data')
      TTMRestartFile=TRIM(RestartFile(1:IndNum-1))//'TTM'//TRIM(RestartFile(IndNum+5:LEN(RestartFile)))
#ifdef MPI
      CALL OpenDataFile(TTMRestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
#else
      CALL OpenDataFile(TTMRestartFile,create=.FALSE.,readOnly=.TRUE.)
#endif
      ALLOCATE( TTM_DG(1:11,0:PP_N,0:PP_N,0:PP_N,PP_nElems) )
      CALL DatasetExists(File_ID,'DG_Solution',TTM_DG_SolutionExists)
      IF(TTM_DG_SolutionExists)THEN
        CALL ReadArray('DG_Solution',5,(/11,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),OffsetElem,5,RealArray=TTM_DG)
      ELSE
        CALL abort(&
          __STAMP__&
          ,' Restart file does not contain "DG_Solution" in restart file for reading TTM data')
      END IF
      CALL CloseDataFile() 
    ELSE! We need to interpolate the solution to the new computational grid
      CALL abort(&
        __STAMP__&
        ,' Restart with changed polynomial degree not implemented for TTM!')
    END IF
  END IF
#endif /*PARTICLES*/

#if USE_QDS_DG
  IF(DoQDS)THEN ! read QDS data from "XXXXX_QDS_000.0XXXXXXXXXXX.h5"
    IndNum=INDEX(RestartFile,'State')
    IF(IndNum.LE.0)CALL abort(&
      __STAMP__&
      ,' Restart file does not contain "State" character string?! Supply restart file for reading QDS data')
    QDSRestartFile=TRIM(RestartFile(1:IndNum-1))//'QDS'//TRIM(RestartFile(IndNum+5:LEN(RestartFile)))
#ifdef MPI
    CALL OpenDataFile(QDSRestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
#else
    CALL OpenDataFile(QDSRestartFile,create=.FALSE.,readOnly=.TRUE.)
#endif
    !ALLOCATE( QDSMacroValues(1:6,0:PP_N,0:PP_N,0:PP_N,PP_nElems) )
    CALL DatasetExists(File_ID,'DG_Solution',QDS_DG_SolutionExists)

    IF(.NOT.QDS_DG_SolutionExists)THEN
      CALL abort(&
        __STAMP__&
        ,' Restart file does not contain "DG_Solution" in restart file for reading QDS data')
    END IF

    IF(.NOT. InterpolateSolution)THEN! No interpolation needed, read solution directly from file
      CALL ReadArray('DG_Solution',5,(/6,PP_N+1,PP_N+1,PP_N+1,nQDSElems/),OffsetElem,5,RealArray=QDSMacroValues)
      DO iElem =1, nQDSElems
        DO k=0, PP_N; DO j=0, PP_N; DO i=0, PP_N
          QDSMacroValues(1  ,i,j,k,iElem) = QDSMacroValues(1  ,i,j,k,iElem)*QDSSpeciesMass!Species(QDS_Species)%MassIC
          QDSMacroValues(2:4,i,j,k,iElem) = QDSMacroValues(2:4,i,j,k,iElem) * QDSMacroValues(1,i,j,k,iElem)
        END DO; END DO; END DO
      END DO
      CALL CloseDataFile() 
    ELSE! We need to interpolate the solution to the new computational grid
      ALLOCATE(U_local(6,0:N_Restart,0:N_Restart,0:N_Restart,nQDSElems))
      CALL ReadArray('DG_Solution',5,(/6,N_Restart+1,N_Restart+1,N_Restart+1,PP_nElems/),OffsetElem,5,RealArray=U_local)
      DO iElem =1, nQDSElems
        DO k=0, N_Restart; DO j=0, N_Restart; DO i=0, N_Restart
          U_local(1  ,i,j,k,iElem) = U_local(1  ,i,j,k,iElem)*QDSSpeciesMass!Species(QDS_Species)%MassIC
          U_local(2:4,i,j,k,iElem) = U_local(2:4,i,j,k,iElem) * U_local(1,i,j,k,iElem)
        END DO; END DO; END DO
        CALL ChangeBasis3D(6,N_Restart,PP_N,Vdm_GaussNRestart_GaussN,U_local(:,:,:,:,iElem),QDSMacroValues(:,:,:,:,iElem))
      END DO
      DEALLOCATE(U_local)
    END IF
  END IF
#endif /*USE_QDS_DG*/

#ifdef MPI
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
#else
  CALL OpenDataFile(RestartFile,create=.FALSE.,readOnly=.TRUE.)
#endif
  ! Read in time from restart file
  !CALL ReadAttribute(File_ID,'Time',1,RealScalar=RestartTime)
  ! Read in state
  IF(.NOT. InterpolateSolution)THEN! No interpolation needed, read solution directly from file
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
      CALL RestartHDG(U) ! calls PostProcessGradient for calculate the derivative, e.g., the electric field E
    ELSE
      lambda=0.
    END IF
#else
    CALL ReadArray('DG_Solution',5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),OffsetElem,5,RealArray=U)
    IF(DoPML)THEN
      ALLOCATE(U_local(PMLnVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
      CALL ReadArray('PML_Solution',5,(/PMLnVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),OffsetElem,5,RealArray=U_local)
      DO iPML=1,nPMLElems
        U2(:,:,:,:,iPML) = U_local(:,:,:,:,PMLToElem(iPML))
      END DO ! iPML
      DEALLOCATE(U_local)
    END IF ! DoPML
#endif
    !CALL ReadState(RestartFile,PP_nVar,PP_N,PP_nElems,U)
  ELSE! We need to interpolate the solution to the new computational grid
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
      ALLOCATE(U_local(PMLnVar,0:N_Restart,0:N_Restart,0:N_Restart,PP_nElems))
      ALLOCATE(U_local2(PMLnVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
      CALL ReadArray('PML_Solution',5,(/PMLnVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),OffsetElem,5,RealArray=U_local)
      DO iElem=1,PP_nElems
        CALL ChangeBasis3D(PMLnVar,N_Restart,PP_N,Vdm_GaussNRestart_GaussN,U_local(:,:,:,:,iElem),U_local2(:,:,:,:,iElem))
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
  CALL DatasetExists(File_ID,'PartInt',PartIntExists)
  IF(PartIntExists)THEN
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
          CALL SingleParticleToExactElement(i,doHALO=.FALSE.,initFix=.FALSE.,doRelocate=.FALSE.)
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
          CALL SingleParticleToExactElementNoMap(i,doHALO=.FALSE.,doRelocate=.FALSE.)
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
          CALL SingleParticleToExactElement(CurrentPartNum,doHALO=.FALSE.,initFix=.FALSE.,doRelocate=.FALSE.)
        ELSE
          CALL SingleParticleToExactElementNoMap(CurrentPartNum,doHALO=.FALSE.,doRelocate=.FALSE.)
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
  ELSE
      SWRITE(UNIT_stdOut,*)'PartData does not exists in restart file'
  END IF ! PartIntExists
#endif /*PARTICLES*/

  CALL CloseDataFile() 

#ifdef PARTICLES
  ! include initially collected particles for first call of field-solver (e.g. in RecomputeLambda)
  CALL ParticleCollectCharges(initialCall_opt=.TRUE.)
#endif /*PARTICLES*/
#ifdef PP_HDG
  iter=0
  ! INSTEAD OF ALL THIS **** DO
  ! 1) MPI-Communication for shape-function particles
  ! 2) Deposition
  ! 3) ONE HDG solve
  CALL  RecomputeLambda(RestartTime)
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
#ifdef PARTICLES
  ! include initially collected particles for first call of field-solver (here because of consistency, but not used until timedisc)
  CALL ParticleCollectCharges(initialCall_opt=.TRUE.)
#endif /*PARTICLES*/
  ! Delete all files since we are doing a fresh start
  IF(DoWriteStateToHDF5) CALL FlushHDF5()
END IF !IF(DoRestart)
END SUBROUTINE Restart

#ifdef PP_HDG
SUBROUTINE RecomputeLambda(t)
!===================================================================================================================================
! The lambda-solution is stored per side, however, the side-list is computed with the OLD domain-decomposition. To allow for
! a change in the load-distribution, number of used cores, etc,... lambda has to be recomputed ONCE
!===================================================================================================================================
! MODULES
USE MOD_DG_Vars,                 ONLY: U
USE MOD_PreProc
USE MOD_HDG,                     ONLY: HDG
USE MOD_TimeDisc_Vars,           ONLY: iter
#ifdef PARTICLES
USE MOD_PICDepo,                 ONLY: Deposition
#ifdef MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY: PartMPIExchange,DoExternalParts
#endif /*MPI*/
#endif /*PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)       :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

#ifdef PARTICLES
#ifdef MPI
IF(DoExternalParts)THEN
  ! communication of shape-function particles, YEAH.
  CALL IRecvNbofParticles() ! open receive buffer for number of particles
  CALL SendNbOfParticles() ! send number of particles
  CALL MPIParticleSend()  ! finish communication of number of particles and send particles
  CALL MPIParticleRecv()  ! finish communication
END IF
#endif

! Deposition of particles
CALL Deposition(doInnerParts=.TRUE.) 
#ifdef MPI
! here: finish deposition with delta kernal
!       maps source terms in physical space
! ALWAYS require
PartMPIExchange%nMPIParticles=0
#endif /*MPI*/
CALL Deposition(doInnerParts=.FALSE.) 
#endif /*PARTICLES*/

! recompute fields
! EM field
CALL HDG(t,U,iter)

END SUBROUTINE RecomputeLambda
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
