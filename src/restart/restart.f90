!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (gitlab.com/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

MODULE MOD_Restart
!===================================================================================================================================
! Module to handle PICLas's restart
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
!CALL prms%CreateLogicalOption('ResetTime', "Override solution time to t=0 on restart.", '.FALSE.')
#if USE_LOADBALANCE
CALL prms%CreateLogicalOption('DoInitialAutoRestart',&
                               "Set Flag for doing automatic initial restart with loadbalancing routines "// &
                               "after first 'InitialAutoRestartSample'-number of iterations.\n"// &
                               "Restart is done if Imbalance > 'Load-DeviationThreshold'."&
                               , '.FALSE.')
CALL prms%CreateIntOption('InitialAutoRestartSample',&
                               "Define number of iterations at simulation start used for ElemTime "// &
                               "sampling before performing automatic initial restart.\n"// &
                               "IF 0 than one iteration is sampled and statefile written has zero timeflag.\n"// &
                               " DEFAULT: LoadBalanceSample.")
CALL prms%CreateLogicalOption( 'InitialAutoRestart-PartWeightLoadBalance', &
                               "Set flag for doing initial auto restart with partMPIWeight instead of"//&
                               " ElemTimes. ElemTime array in state file is filled with nParts*PartMPIWeight for each Elem. "//&
                               " If Flag [TRUE] InitialAutoRestartSample is set to 0 and vice versa.", '.FALSE.')
#endif /*USE_LOADBALANCE*/
CALL prms%CreateLogicalOption( 'RestartNullifySolution', &
                               "Set the DG solution to zero (ignore the DG solution in the state file)",&
                               '.FALSE.')
CALL prms%CreateLogicalOption('Particles-MacroscopicRestart', &
                              "TO-DO",&
                              '.FALSE.')
CALL prms%CreateStringOption( 'Particles-MacroscopicRestart-Filename', &
                              'TO-DO')
END SUBROUTINE DefineParametersRestart


SUBROUTINE InitRestart()
!===================================================================================================================================
! Initialize all necessary information to perform the restart
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: FileVersionHDF5
USE MOD_PreProc
USE MOD_ReadInTools        ,ONLY: GETLOGICAL,GETSTR
#if USE_LOADBALANCE
USE MOD_ReadInTools        ,ONLY: GETINT
USE MOD_LoadBalance_Vars   ,ONLY: LoadBalanceSample
USE MOD_ReadInTools        ,ONLY: PrintOption
#endif /*USE_LOADBALANCE*/
USE MOD_Interpolation_Vars ,ONLY: xGP,InterpolationInitIsDone
USE MOD_Restart_Vars
USE MOD_HDF5_Input         ,ONLY: OpenDataFile,CloseDataFile,GetDataProps,ReadAttribute,File_ID
USE MOD_HDF5_Input         ,ONLY: DatasetExists
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_LOADBALANCE
CHARACTER(20)               :: hilf
#endif /*USE_LOADBALANCE*/
#if USE_HDG
LOGICAL                     :: DG_SolutionUExists
#endif /*USE_HDG*/
LOGICAL                     :: FileVersionExists
!===================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.RestartInitIsDone)THEN
   CALL abort(&
__STAMP__&
,'InitRestart not ready to be called or already called.',999,999.)
   RETURN
END IF

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT RESTART...'

! Set the DG solution to zero (ignore the DG solution in the state file)
RestartNullifySolution = GETLOGICAL('RestartNullifySolution','F')

! Macroscopic restart
DoMacroscopicRestart = GETLOGICAL('Particles-MacroscopicRestart')
IF(DoMacroscopicRestart) MacroRestartFileName = GETSTR('Particles-MacroscopicRestart-Filename')

! Check if we want to perform a restart
IF (LEN_TRIM(RestartFile).GT.0) THEN
  ! Read in the state file we want to restart from
  DoRestart = .TRUE.
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
#ifdef PP_POIS
#if (PP_nVar==8)
  !The following arrays are read from the file
  !CALL ReadArray('DG_SolutionE',5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),OffsetElem,5,RealArray=U)
  !CALL ReadArray('DG_SolutionPhi',5,(/4,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),OffsetElem,5,RealArray=Phi)
  CALL abort(&
      __STAMP__&
      ,'InitRestart: This case is not implemented here. Fix this!')
#else
  !The following arrays are read from the file
  !CALL ReadArray('DG_SolutionE',5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),OffsetElem,5,RealArray=U)
  !CALL ReadArray('DG_SolutionPhi',5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),OffsetElem,5,RealArray=Phi)
  CALL abort(&
      __STAMP__&
      ,'InitRestart: This case is not implemented here. Fix this!')
#endif
#elif USE_HDG
  CALL DatasetExists(File_ID,'DG_SolutionU',DG_SolutionUExists)
  IF(DG_SolutionUExists)THEN
    CALL GetDataProps('DG_SolutionU',nVar_Restart,N_Restart,nElems_Restart,NodeType_Restart)
  END IF
#else
  CALL GetDataProps('DG_Solution',nVar_Restart,N_Restart,nElems_Restart,NodeType_Restart)
#endif
  IF(RestartNullifySolution)THEN ! Open the restart file and neglect the DG solution (only read particles if present)
    SWRITE(UNIT_stdOut,*)' | Restarting from File: "',TRIM(RestartFile),'" (but without reading the DG solution)'
  ELSE ! Use the solution in the restart file
    SWRITE(UNIT_stdOut,*)' | Restarting from File: "',TRIM(RestartFile),'"'
    IF(PP_nVar.NE.nVar_Restart)THEN
      SWRITE(UNIT_StdOut,'(A,I5)')"     PP_nVar =", PP_nVar
      SWRITE(UNIT_StdOut,'(A,I5)')"nVar_Restart =", nVar_Restart
      CALL abort(&
          __STAMP__&
          ,'InitRestart: PP_nVar.NE.nVar_Restart (number of variables in restat file does no match the compiled equation system).')
    END IF
  END IF
  ! Read in time from restart file
  CALL ReadAttribute(File_ID,'Time',1,RealScalar=RestartTime)
  ! check file version
  CALL DatasetExists(File_ID,'File_Version',FileVersionExists,attrib=.TRUE.)
  IF (FileVersionExists) THEN
    CALL ReadAttribute(File_ID,'File_Version',1,RealScalar=FileVersionHDF5)
  ELSE
    CALL abort(&
        __STAMP__&
        ,'Error in InitRestart(): Attribute "File_Version" does not exist!')
  END IF
  IF(FileVersionHDF5.LT.1.5)THEN
    SWRITE(UNIT_StdOut,'(A)')' '
    SWRITE(UNIT_StdOut,'(A)')' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
    SWRITE(UNIT_StdOut,'(A)')' '
    SWRITE(UNIT_StdOut,'(A)')' Restart file is too old! "File_Version" in restart file < 1.5!'
    SWRITE(UNIT_StdOut,'(A)')' The format used in the restart file is not compatible with this version of PICLas.'
    SWRITE(UNIT_StdOut,'(A)')' Among others, the particle format (PartData) has changed.'
    SWRITE(UNIT_StdOut,'(A)')' Run python script '
    SWRITE(UNIT_StdOut,'(A)')' '
    SWRITE(UNIT_StdOut,'(A)')'     python  ./tools/flip_PartState/flip_PartState.py  --help'
    SWRITE(UNIT_StdOut,'(A)')' '
    SWRITE(UNIT_StdOut,'(A)')' for info regarding the usage and run the script against the restart file, e.g., '
    SWRITE(UNIT_StdOut,'(A)')' '
    SWRITE(UNIT_StdOut,'(A)')'     python  ./tools/flip_PartState/flip_PartState.py  ProjectName_State_000.0000xxxxxx.h5'
    SWRITE(UNIT_StdOut,'(A)')' '
    SWRITE(UNIT_StdOut,'(A)')' to update the format and file version number.'
    SWRITE(UNIT_StdOut,'(A)')' Note that the format can be changed back to the old one by running the script a second time.'
    SWRITE(UNIT_StdOut,'(A)')' '
    SWRITE(UNIT_StdOut,'(A)')' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
    CALL abort(&
    __STAMP__&
    ,'Error in InitRestart(): "File_Version" in restart file < 1.5. See error message above to fix. File version in restart file =',&
    RealInfoOpt=FileVersionHDF5)
  END IF ! FileVersionHDF5.LT.1.5
  CALL CloseDataFile()
ELSE
  RestartTime = 0.
  SWRITE(UNIT_StdOut,'(A)')' | No restart wanted, doing a fresh computation!'
END IF

! Automatically do a load balance step at the beginning of a new simulation or a user-restarted simulation
#if USE_LOADBALANCE
DoInitialAutoRestart = GETLOGICAL('DoInitialAutoRestart')
IF(nProcessors.LT.2) DoInitialAutoRestart = .FALSE.
WRITE(UNIT=hilf,FMT='(I0)') LoadBalanceSample
InitialAutoRestartSample = GETINT('InitialAutoRestartSample',TRIM(hilf))
IAR_PerformPartWeightLB = GETLOGICAL('InitialAutoRestart-PartWeightLoadBalance','F')
IF (IAR_PerformPartWeightLB) THEN
  InitialAutoRestartSample = 0 ! deactivate loadbalance sampling of ElemTimes if balancing with partweight is enabled
  CALL PrintOption('InitialAutoRestart-PartWeightLoadBalance = T : InitialAutoRestartSample','INFO',IntOpt=InitialAutoRestartSample)
ELSE IF (InitialAutoRestartSample.EQ.0) THEN
  IAR_PerformPartWeightLB = .TRUE. ! loadbalance (ElemTimes) is done with partmpiweight if loadbalancesampling is set to zero
  CALL PrintOption('InitialAutoRestart-PartWeightLoadBalance','INFO',LogOpt=IAR_PerformPartWeightLB)
END IF
#endif /*USE_LOADBALANCE*/

! Set wall time to the beginning of the simulation or when a restart is performed to the current wall time
RestartWallTime=PICLASTIME()

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
USE MOD_IO_HDF5
USE MOD_DG_Vars                ,ONLY: U
USE MOD_Mesh_Vars              ,ONLY: OffsetElem,DoWriteStateToHDF5
#if USE_HDG
USE MOD_Mesh_Vars              ,ONLY: offsetSide,nSides,nMPISides_YOUR, offsetSide
#endif
#if (USE_QDS_DG) || ! (USE_HDG)
USE MOD_Restart_Vars           ,ONLY: Vdm_GaussNRestart_GaussN
#endif /*USE_HDG*/
USE MOD_Restart_Vars           ,ONLY: DoRestart,N_Restart,RestartFile,RestartTime,InterpolateSolution,RestartNullifySolution
USE MOD_ChangeBasis            ,ONLY: ChangeBasis3D
USE MOD_HDF5_input             ,ONLY: OpenDataFile,CloseDataFile,ReadArray,ReadAttribute,GetDataSize
USE MOD_HDF5_Output            ,ONLY: FlushHDF5
#if ! (USE_HDG)
USE MOD_PML_Vars               ,ONLY: DoPML,PMLToElem,U2,nPMLElems,PMLnVar
#endif /*not USE_HDG*/
#ifdef PP_POIS
USE MOD_Equation_Vars          ,ONLY: Phi
#endif /*PP_POIS*/
#ifdef PARTICLES
USE MOD_Restart_Tools          ,ONLY: ReadNodeSourceExtFromHDF5
USE MOD_Restart_Vars           ,ONLY: DoMacroscopicRestart
USE MOD_Particle_Vars          ,ONLY: PartState, PartSpecies, PEM, PDM, nSpecies, usevMPF, PartMPF,PartPosRef, SpecReset, Species
USE MOD_part_tools             ,ONLY: UpdateNextFreePosition
USE MOD_DSMC_Vars              ,ONLY: UseDSMC,CollisMode,PartStateIntEn,DSMC,VibQuantsPar,PolyatomMolDSMC,SpecDSMC,RadialWeighting
USE MOD_Eval_XYZ               ,ONLY: GetPositionInRefElem
USE MOD_Particle_Mesh          ,ONLY: SingleParticleToExactElement,SingleParticleToExactElementNoMap
USE MOD_Particle_Mesh_Tools    ,ONLY: ParticleInsideQuad3D
USE MOD_Particle_Mesh_Vars     ,ONLY: epsOneCell
USE MOD_Particle_Tracking_Vars ,ONLY: DoRefMapping, TriaTracking
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_SurfaceModel_Vars      ,ONLY: SurfDistInfo, Adsorption
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfBC
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample,SurfMesh,offSetSurfSide,PartBound,nPartBound
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*USE_MPI*/
USE MOD_Particle_Tracking      ,ONLY: ParticleCollectCharges
USE MOD_PICDepo_Vars           ,ONLY: DoDeposition, RelaxDeposition, PartSourceOld
USE MOD_Dielectric_Vars        ,ONLY: DoDielectricSurfaceCharge
#endif /*PARTICLES*/
#if USE_HDG
USE MOD_HDG_Vars               ,ONLY: lambda, nGP_face
USE MOD_HDG                    ,ONLY: RestartHDG
#endif /*USE_HDG*/
#if USE_QDS_DG
USE MOD_QDS_DG_Vars            ,ONLY: DoQDS,QDSMacroValues,nQDSElems,QDSSpeciesMass
#endif /*USE_QDS_DG*/
#if (USE_QDS_DG) || (PARTICLES) || (USE_HDG)
USE MOD_HDF5_Input             ,ONLY: File_ID,DatasetExists,nDims,HSize
#endif

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if (USE_QDS_DG) || !(USE_HDG)
REAL,ALLOCATABLE                   :: U_local(:,:,:,:,:)
REAL,ALLOCATABLE                   :: U_local2(:,:,:,:,:)
INTEGER                            :: iPML
#endif
#if USE_HDG
LOGICAL                            :: DG_SolutionLambdaExists,DG_SolutionUExists
INTEGER(KIND=8)                    :: iter
#endif /*USE_HDG*/
INTEGER                            :: iElem
#if USE_MPI
REAL                               :: StartT,EndT
#endif /*USE_MPI*/
#ifdef PARTICLES
CHARACTER(LEN=255),ALLOCATABLE     :: StrVarNames(:)
CHARACTER(LEN=255),ALLOCATABLE     :: StrVarNames_HDF5(:)
INTEGER                            :: FirstElemInd,LastelemInd,j,k
INTEGER(KIND=IK),ALLOCATABLE       :: PartInt(:,:)
INTEGER,PARAMETER                  :: PartIntSize=2                  ! number of entries in each line of PartInt
INTEGER                            :: PartDataSize,PartDataSize_HDF5 ! number of entries in each line of PartData
INTEGER(KIND=IK)                   :: locnPart,offsetnPart,iLoop
INTEGER,PARAMETER                  :: ELEM_FirstPartInd=1
INTEGER,PARAMETER                  :: ELEM_LastPartInd=2
REAL,ALLOCATABLE                   :: PartData(:,:)
REAL                               :: xi(3)
LOGICAL                            :: InElementCheck,PartIntExists,PartDataExists,VibQuantDataExists,changedVars,DGSourceExists
REAL                               :: det(6,2)
INTEGER                            :: COUNTER, COUNTER2, CounterPoly
INTEGER, ALLOCATABLE               :: VibQuantData(:,:)
INTEGER                            :: MaxQuantNum, iPolyatMole, iSpec, iPart, iVar
! 2D Symmetry RadialWeighting
LOGICAL                            :: CloneExists
#if USE_MPI
REAL, ALLOCATABLE                  :: SendBuff(:), RecBuff(:)
INTEGER                            :: LostParts(0:PartMPI%nProcs-1), Displace(0:PartMPI%nProcs-1),CurrentPartNum
INTEGER                            :: NbrOfFoundParts, CompleteNbrOfFound, RecCount(0:PartMPI%nProcs-1)
INTEGER, ALLOCATABLE               :: SendBuffPoly(:), RecBuffPoly(:)
INTEGER                            :: LostPartsPoly(0:PartMPI%nProcs-1), DisplacePoly(0:PartMPI%nProcs-1)
#endif /*USE_MPI*/
INTEGER                            :: locnSurfPart,offsetnSurfPart
INTEGER,ALLOCATABLE                :: SurfPartInt(:,:,:,:,:)
INTEGER,ALLOCATABLE                :: SurfPartData(:,:)
REAL,ALLOCATABLE                   :: SurfCalcData(:,:,:,:,:)
REAL,ALLOCATABLE                   :: PartSource_HDF5(:,:,:,:,:)
INTEGER                            :: Coordinations, SurfPartIntSize, SurfPartDataSize
INTEGER                            :: UsedSiteMapPos, nVar, nfreeArrayindeces, lastfreeIndx, current
INTEGER                            :: xpos, ypos, firstpart, lastpart, PartBoundID, SideID
INTEGER                            :: iCoord, SpecID, iSurfSide, isubsurf, jsubsurf, iInterAtom
INTEGER                            :: nSpecies_HDF5, nSurfSample_HDF5, nSurfBC_HDF5, nSurfSides_HDF5
LOGICAL                            :: SurfCalcDataExists, SurfPartIntExists, SurfPartDataExists, MoveToLastFree, implemented
LOGICAL,ALLOCATABLE                :: readVarFromState(:)
LOGICAL                            :: WallmodelExists(1:nPartBound), SurfModelTypeExists
INTEGER,ALLOCATABLE                :: SurfModelType(:)
#endif /*PARTICLES*/
#if USE_QDS_DG
CHARACTER(255)                     :: QDSRestartFile        ! > QDS Data file for restart
LOGICAL                            :: QDS_DG_SolutionExists
INTEGER                            :: j,k
INTEGER                            :: IndNum    ! > auxiliary variable containing the index number of a substring within a string
#endif /*USE_QDS_DG*/
#if (USE_QDS_DG) || (PARTICLES)
INTEGER                            :: i
#endif
INTEGER(KIND=IK)                   :: PP_NTmp,OffsetElemTmp,PP_nVarTmp,PP_nElemsTmp,N_RestartTmp
#if !(USE_HDG)
INTEGER(KIND=IK)                   :: PMLnVarTmp
#endif /*not USE_HDG*/
!===================================================================================================================================
IF(DoRestart)THEN
#if USE_MPI
  StartT=MPI_WTIME()
#endif
  ! Temp. vars for integer KIND=8 possibility
  PP_NTmp       = INT(PP_N,IK)
  OffsetElemTmp = INT(OffsetElem,IK)
  PP_nVarTmp    = INT(PP_nVar,IK)
  PP_nElemsTmp  = INT(PP_nElems,IK)
  N_RestartTmp  = INT(N_Restart,IK)
#if !(USE_HDG)
  PMLnVarTmp    = INT(PMLnVar,IK)
#endif /*not USE_HDG*/
  ! ===========================================================================
  ! 1.) Read the field solution
  ! ===========================================================================
  RestartNullifySolution=.False.
  IF(RestartNullifySolution)THEN ! Open the restart file and neglect the DG solution (only read particles if present)
    SWRITE(UNIT_stdOut,*)'Restarting from File: ',TRIM(RestartFile),' (but without reading the DG solution)'
    CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)

  ELSE ! Use the solution in the restart file
    SWRITE(UNIT_stdOut,*)'Restarting from File: ',TRIM(RestartFile)
#if USE_QDS_DG
    IF(DoQDS)THEN ! read QDS data from "XXXXX_QDS_000.0XXXXXXXXXXX.h5"
      IndNum=INDEX(RestartFile,'State')
      IF(IndNum.LE.0)CALL abort(&
          __STAMP__&
          ,' Restart file does not contain "State" character string?! Supply restart file for reading QDS data')
      QDSRestartFile=TRIM(RestartFile(1:IndNum-1))//'QDS'//TRIM(RestartFile(IndNum+5:LEN(RestartFile)))
      CALL OpenDataFile(QDSRestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
      !ALLOCATE( QDSMacroValues(1:6,0:PP_N,0:PP_N,0:PP_N,PP_nElems) )
      CALL DatasetExists(File_ID,'DG_Solution',QDS_DG_SolutionExists)

      IF(.NOT.QDS_DG_SolutionExists)THEN
        CALL abort(&
            __STAMP__&
            ,' Restart file does not contain "DG_Solution" in restart file for reading QDS data')
      END IF

      IF(.NOT. InterpolateSolution)THEN! No interpolation needed, read solution directly from file
        CALL ReadArray('DG_Solution',5,(/6_IK,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_NTmp+1_IK,INT(nQDSElems,IK)/),&
            OffsetElemTmp,5,RealArray=QDSMacroValues)
        DO iElem =1, nQDSElems
          DO k=0, PP_N; DO j=0, PP_N; DO i=0, PP_N
            QDSMacroValues(1  ,i,j,k,iElem) = QDSMacroValues(1  ,i,j,k,iElem)*QDSSpeciesMass!Species(QDS_Species)%MassIC
            QDSMacroValues(2:4,i,j,k,iElem) = QDSMacroValues(2:4,i,j,k,iElem) * QDSMacroValues(1,i,j,k,iElem)
          END DO; END DO; END DO
        END DO
      ELSE! We need to interpolate the solution to the new computational grid
        ALLOCATE(U_local(6,0:N_Restart,0:N_Restart,0:N_Restart,nQDSElems))
        CALL ReadArray('DG_Solution',5,(/6_IK,N_RestartTmp+1_IK,N_RestartTmp+1_IK,N_RestartTmp+1_IK,INT(nQDSElems,IK)/),&
            OffsetElemTmp,5,RealArray=U_local)
        DO iElem =1, nQDSElems
          DO k=0, N_Restart; DO j=0, N_Restart; DO i=0, N_Restart
            U_local(1  ,i,j,k,iElem) = U_local(1  ,i,j,k,iElem)*QDSSpeciesMass!Species(QDS_Species)%MassIC
            U_local(2:4,i,j,k,iElem) = U_local(2:4,i,j,k,iElem) * U_local(1,i,j,k,iElem)
          END DO; END DO; END DO
          CALL ChangeBasis3D(6,N_Restart,PP_N,Vdm_GaussNRestart_GaussN,U_local(:,:,:,:,iElem),QDSMacroValues(:,:,:,:,iElem))
        END DO
        DEALLOCATE(U_local)
      END IF
      CALL CloseDataFile()
    END IF ! DoQDS
#endif /*USE_QDS_DG*/

    CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
#ifdef PARTICLES
    !-- read PartSource if relaxation is performed (might be needed for RestartHDG)
    IF (DoDeposition .AND. RelaxDeposition) THEN
      CALL DatasetExists(File_ID,'DG_Source',DGSourceExists)
      IF(DGSourceExists)THEN
        IF(.NOT.InterpolateSolution)THEN! No interpolation needed, read solution directly from file
          ALLOCATE(PartSource_HDF5(1:4,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
          CALL ReadArray('DG_Source' ,5,(/4_IK,PP_NTmp+1,PP_NTmp+1,PP_NTmp+1,PP_nElemsTmp/),OffsetElemTmp,5,RealArray=PartSource_HDF5)
          DO iElem =1, PP_nElems
            DO k=0, PP_N; DO j=0, PP_N; DO i=0, PP_N
#if ((USE_HDG) && (PP_nVar==1))
              PartSourceOld(1,1,i,j,k,iElem) = PartSource_HDF5(4,i,j,k,iElem)
              PartSourceOld(1,2,i,j,k,iElem) = PartSource_HDF5(4,i,j,k,iElem)
#else
              PartSourceOld(1:4,1,i,j,k,iElem) = PartSource_HDF5(1:4,i,j,k,iElem)
              PartSourceOld(1:4,2,i,j,k,iElem) = PartSource_HDF5(1:4,i,j,k,iElem)
#endif
            END DO; END DO; END DO
          END DO
          DEALLOCATE(PartSource_HDF5)
        ELSE! We need to interpolate the solution to the new computational grid
          CALL abort(&
              __STAMP__&
              ,' Restart with changed polynomial degree not implemented for DG_Source!')
        END IF
      END IF
    END IF
#endif /*PARTICLES*/
    ! Read in time from restart file
    !CALL ReadAttribute(File_ID,'Time',1,RealScalar=RestartTime)
    ! Read in state
    IF(.NOT. InterpolateSolution)THEN! No interpolation needed, read solution directly from file
#ifdef PP_POIS
#if (PP_nVar==8)
      CALL ReadArray('DG_SolutionE',5,(/PP_nVarTmp,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_nElemsTmp/),OffsetElemTmp,5,RealArray=U)
      CALL ReadArray('DG_SolutionPhi',5,(/4_IK,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_nElemsTmp/),OffsetElemTmp,5,RealArray=Phi)
#else
      CALL ReadArray('DG_SolutionE',5,(/PP_nVarTmp,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_nElemsTmp/),OffsetElemTmp,5,RealArray=U)
      CALL ReadArray('DG_SolutionPhi',5,(/PP_nVarTmp,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_nElemsTmp/),OffsetElemTmp,5,RealArray=Phi)
#endif
#elif USE_HDG
      CALL DatasetExists(File_ID,'DG_SolutionU',DG_SolutionUExists)
      IF(DG_SolutionUExists)THEN
        CALL ReadArray('DG_SolutionU',5,(/PP_nVarTmp,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_nElemsTmp/),OffsetElemTmp,5,RealArray=U)
      ELSE
        ! CALL abort(&
        !     __STAMP__&
        !     ,' DG_SolutionU does not exist in restart-file!')
        ! !DG_Solution contains a 4er-/3er-/7er-array, not PP_nVar!!!
        CALL ReadArray('DG_Solution' ,5,(/PP_nVarTmp,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_nElemsTmp/),OffsetElemTmp,5,RealArray=U)
      END IF
      CALL DatasetExists(File_ID,'DG_SolutionLambda',DG_SolutionLambdaExists)
      IF(DG_SolutionLambdaExists)THEN
        CALL ReadArray('DG_SolutionLambda',3,(/PP_nVarTmp,nGP_face,nSides-nMPISides_YOUR/),offsetSide,3,RealArray=lambda)
        CALL RestartHDG(U) ! calls PostProcessGradient for calculate the derivative, e.g., the electric field E
      ELSE
        lambda=0.
      END IF
#else
      CALL ReadArray('DG_Solution',5,(/PP_nVarTmp,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_nElemsTmp/),OffsetElemTmp,5,RealArray=U)
      IF(DoPML)THEN
        ALLOCATE(U_local(PMLnVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
        CALL ReadArray('PML_Solution',5,(/INT(PMLnVar,IK),PP_NTmp+1_IK,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_nElemsTmp/),&
            OffsetElemTmp,5,RealArray=U_local)
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
      CALL ReadArray('DG_SolutionE',5,(/PP_nVarTmp,N_RestartTmp+1_IK,N_RestartTmp+1_IK,N_RestartTmp+1_IK,PP_nElemsTmp/),&
          OffsetElemTmp,5,RealArray=U_local)
      DO iElem=1,PP_nElems
        CALL ChangeBasis3D(PP_nVar,N_Restart,PP_N,Vdm_GaussNRestart_GaussN,U_local(:,:,:,:,iElem),U(:,:,:,:,iElem))
      END DO
      DEALLOCATE(U_local)

      ALLOCATE(U_local(4,0:N_Restart,0:N_Restart,0:N_Restart,PP_nElems))
      CALL ReadArray('DG_SolutionPhi',5,(/4_IK,N_RestartTmp+1_IK,N_RestartTmp+1_IK,N_RestartTmp+1_IK,PP_nElemsTmp/),&
          OffsetElemTmp,5,RealArray=U_local)
      DO iElem=1,PP_nElems
        CALL ChangeBasis3D(4,N_Restart,PP_N,Vdm_GaussNRestart_GaussN,U_local(:,:,:,:,iElem),Phi(:,:,:,:,iElem))
      END DO
      DEALLOCATE(U_local)
#else
      ALLOCATE(U_local(PP_nVar,0:N_Restart,0:N_Restart,0:N_Restart,PP_nElems))
      CALL ReadArray('DG_SolutionE',5,(/PP_nVarTmp,N_RestartTmp+1_IK,N_RestartTmp+1_IK,N_RestartTmp+1_IK,PP_nElemsTmp/),&
          OffsetElemTmp,5,RealArray=U_local)
      DO iElem=1,PP_nElems
        CALL ChangeBasis3D(PP_nVar,N_Restart,PP_N,Vdm_GaussNRestart_GaussN,U_local(:,:,:,:,iElem),U(:,:,:,:,iElem))
      END DO
      CALL ReadArray('DG_SolutionPhi',5,(/PP_nVarTmp,N_RestartTmp+1_IK,N_RestartTmp+1_IK,N_RestartTmp+1_IK,PP_nElemsTmp/),&
          OffsetElemTmp,5,RealArray=U_local)
      DO iElem=1,PP_nElems
        CALL ChangeBasis3D(PP_nVar,N_Restart,PP_N,Vdm_GaussNRestart_GaussN,U_local(:,:,:,:,iElem),Phi(:,:,:,:,iElem))
      END DO
      DEALLOCATE(U_local)
#endif
#elif USE_HDG
      CALL abort(&
          __STAMP__&
          ,' Restart with changed polynomial degree not implemented for HDG!')
      !    ALLOCATE(U_local(PP_nVar,0:N_Restart,0:N_Restart,0:N_Restart,PP_nElems))
      !    CALL ReadArray('DG_SolutionLambda',5,(/PP_nVar,N_RestartTmp+1_IK,N_RestartTmp+1_IK,N_RestartTmp+1_IK,PP_nElemsTmp/),OffsetElem,5,RealArray=U_local)
      !    DO iElem=1,PP_nElems
      !      CALL ChangeBasis3D(PP_nVar,N_Restart,PP_N,Vdm_GaussNRestart_GaussN,U_local(:,:,:,:,iElem),U(:,:,:,:,iElem))
      !    END DO
      !    DEALLOCATE(U_local)
      !CALL RestartHDG(U)
#else
      ALLOCATE(U_local(PP_nVar,0:N_Restart,0:N_Restart,0:N_Restart,PP_nElems))
      CALL ReadArray('DG_Solution',5,(/PP_nVarTmp,N_RestartTmp+1_IK,N_RestartTmp+1_IK,N_RestartTmp+1_IK,PP_nElemsTmp/),&
          OffsetElemTmp,5,RealArray=U_local)
      DO iElem=1,PP_nElems
        CALL ChangeBasis3D(PP_nVar,N_Restart,PP_N,Vdm_GaussNRestart_GaussN,U_local(:,:,:,:,iElem),U(:,:,:,:,iElem))
      END DO
      DEALLOCATE(U_local)
      IF(DoPML)THEN
        ALLOCATE(U_local(PMLnVar,0:N_Restart,0:N_Restart,0:N_Restart,PP_nElems))
        ALLOCATE(U_local2(PMLnVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
        CALL ReadArray('PML_Solution',5,(/INT(PMLnVar,IK),PP_NTmp+1_IK,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_nElemsTmp/),&
            OffsetElemTmp,5,RealArray=U_local)
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
    END IF ! IF(.NOT. InterpolateSolution)
  END IF ! IF(.NOT. RestartNullifySolution)

  ! ------------------------------------------------
  ! NodeSourceExt (external/additional charge source terms)
  ! ------------------------------------------------
#if PARTICLES
  IF(DoDielectricSurfaceCharge) CALL ReadNodeSourceExtFromHDF5()
#endif /*PARTICLES*/


#ifdef PARTICLES
  ! ===========================================================================
  ! 2.) Read the particle solution
  ! ===========================================================================
  implemented=.FALSE.
  IF(.NOT.DoMacroscopicRestart) THEN
    IF(useDSMC)THEN
      IF((CollisMode.GT.1).AND.(usevMPF).AND.(DSMC%ElectronicModel))THEN
        PartDataSize=11
        ALLOCATE(StrVarNames(PartDataSize))
        StrVarNames( 8)='Vibrational'
        StrVarNames( 9)='Rotational'
        StrVarNames(10)='Electronic'
        StrVarNames(11)='MPF'
        implemented = .TRUE.
      ELSE IF ( (CollisMode .GT. 1) .AND. (usevMPF) ) THEN
        PartDataSize=10
        ALLOCATE(StrVarNames(PartDataSize))
        StrVarNames( 8)='Vibrational'
        StrVarNames( 9)='Rotational'
        StrVarNames(10)='MPF'
        implemented = .TRUE.
      ELSE IF ( (CollisMode .GT. 1) .AND. (DSMC%ElectronicModel) ) THEN
        PartDataSize=10
        ALLOCATE(StrVarNames(PartDataSize))
        StrVarNames( 8)='Vibrational'
        StrVarNames( 9)='Rotational'
        StrVarNames(10)='Electronic'
      ELSE IF (CollisMode.GT.1) THEN
        implemented=.TRUE.
        PartDataSize=9 !int ener + 2
        ALLOCATE(StrVarNames(PartDataSize))
        StrVarNames( 8)='Vibrational'
        StrVarNames( 9)='Rotational'
      ELSE IF (usevMPF) THEN
        PartDataSize=8 !+ 1 vmpf
        ALLOCATE(StrVarNames(PartDataSize))
        StrVarNames( 8)='MPF'
        implemented=.TRUE.
      ELSE
        PartDataSize=7 !+ 0
        ALLOCATE(StrVarNames(PartDataSize))
      END IF
    ELSE IF (usevMPF) THEN
      PartDataSize=8 !vmpf +1
      ALLOCATE(StrVarNames(PartDataSize))
      StrVarNames( 8)='MPF'
    ELSE
      PartDataSize=7
      ALLOCATE(StrVarNames(PartDataSize))
    END IF ! UseDSMC
    StrVarNames(1)='ParticlePositionX'
    StrVarNames(2)='ParticlePositionY'
    StrVarNames(3)='ParticlePositionZ'
    StrVarNames(4)='VelocityX'
    StrVarNames(5)='VelocityY'
    StrVarNames(6)='VelocityZ'
    StrVarNames(7)='Species'
    ALLOCATE(readVarFromState(PartDataSize))
    readVarFromState=.TRUE.

    IF (useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
      MaxQuantNum = 0
      DO iSpec = 1, nSpecies
        IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
          iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
          IF (PolyatomMolDSMC(iPolyatMole)%VibDOF.GT.MaxQuantNum) MaxQuantNum = PolyatomMolDSMC(iPolyatMole)%VibDOF
        END IF ! SpecDSMC(iSpec)%PolyatomicMol
      END DO ! iSpec = 1, nSpecies
    END IF ! useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)

    SWRITE(UNIT_stdOut,*)'Reading Particles from Restartfile...'
    !read local ElemInfo from HDF5
    FirstElemInd=offsetElem+1
    LastElemInd=offsetElem+PP_nElems
    ! read local ParticleInfo from HDF5
    CALL DatasetExists(File_ID,'PartInt',PartIntExists)
    IF(PartIntExists)THEN
      ALLOCATE(PartInt(FirstElemInd:LastElemInd,PartIntSize))

      ! Associate construct for integer KIND=8 possibility
      ASSOCIATE (&
            PP_nElems   => INT(PP_nElems,IK)   ,&
            PartIntSize => INT(PartIntSize,IK) ,&
            offsetElem  => INT(offsetElem,IK)   )
        CALL ReadArray('PartInt',2,(/PP_nElems,PartIntSize/),offsetElem,1,IntegerArray=PartInt)
      END ASSOCIATE
      ! read local Particle Data from HDF5
      locnPart=PartInt(LastElemInd,ELEM_LastPartInd)-PartInt(FirstElemInd,ELEM_FirstPartInd)
      offsetnPart=PartInt(FirstElemInd,ELEM_FirstPartInd)
      CALL DatasetExists(File_ID,'PartData',PartDataExists)
      IF(PartDataExists)THEN
        ! Read in parameters from the State file
        CALL GetDataSize(File_ID,'VarNamesParticles',nDims,HSize,attrib=.TRUE.)
        PartDataSize_HDF5 = INT(HSize(1),4)
        ALLOCATE(StrVarNames_HDF5(PartDataSize_HDF5))
        CALL ReadAttribute(File_ID,'VarNamesParticles',PartDataSize_HDF5,StrArray=StrVarNames_HDF5)
        IF (PartDataSize_HDF5.NE.PartDataSize) THEN
          changedVars=.TRUE.
        ELSE IF (.NOT.ALL(StrVarNames_HDF5.EQ.StrVarNames)) THEN
          changedVars=.TRUE.
        ELSE
          changedVars=.FALSE.
        END IF ! PartDataSize_HDF5.NE.PartDataSize
        IF (changedVars) THEN
          SWRITE(*,*) 'WARNING: VarNamesParticles have changed from restart-file!!!'
          IF (.NOT.implemented) CALL Abort(&
              __STAMP__&
              ,"not implemented yet!")
          readVarFromState=.FALSE.
          DO iVar=1,PartDataSize_HDF5
            IF (TRIM(StrVarNames(iVar)).EQ.TRIM(StrVarNames_HDF5(iVar))) THEN
              readVarFromState(iVar)=.TRUE.
            ELSE
              CALL Abort(&
                  __STAMP__&
                  ,"not associated VarNamesParticles in HDF5!")
            END IF
          END DO ! iVar=1,PartDataSize_HDF5
          DO iVar=1,PartDataSize
            IF (.NOT.readVarFromState(iVar)) THEN
              IF (TRIM(StrVarNames(iVar)).EQ.'Vibrational' .OR. TRIM(StrVarNames(iVar)).EQ.'Rotational') THEN
                SWRITE(*,*) 'WARNING: The following VarNamesParticles will be set to zero: '//TRIM(StrVarNames(iVar))
              ELSE IF(TRIM(StrVarNames(iVar)).EQ.'MPF') THEN
                SWRITE(*,*) 'WARNING: The particle weighting factor will be initialized with the given global weighting factor!'
              ELSE
                CALL Abort(&
                    __STAMP__&
                    ,"not associated VarNamesParticles to be reset!")
              END IF ! TRIM(StrVarNames(iVar)).EQ.'Vibrational' .OR. TRIM(StrVarNames(iVar)).EQ.'Rotational'
            END IF ! .NOT.readVarFromState(iVar)
          END DO ! iVar=1,PartDataSize
        END IF ! changedVars
        ALLOCATE(PartData(PartDataSize_HDF5,offsetnPart+1_IK:offsetnPart+locnPart))

        CALL ReadArray('PartData',2,(/INT(PartDataSize_HDF5,IK),locnPart/),offsetnPart,2,RealArray=PartData)!,&
        !xfer_mode_independent=.TRUE.)

        IF (useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
          CALL DatasetExists(File_ID,'VibQuantData',VibQuantDataExists)
          IF (.NOT.VibQuantDataExists) CALL abort(&
              __STAMP__&
              ,' Restart file does not contain "VibQuantData" in restart file for reading of polyatomic data')
          ALLOCATE(VibQuantData(MaxQuantNum,offsetnPart+1_IK:offsetnPart+locnPart))

          CALL ReadArray('VibQuantData',2,(/INT(MaxQuantNum,IK),locnPart/),offsetnPart,2,IntegerArray_i4=VibQuantData)
          !+1 is real number of necessary vib quants for the particle
        END IF ! useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)

        iPart=0
        DO iLoop = 1_IK,locnPart
          IF(SpecReset(INT(PartData(7,offsetnPart+iLoop),4))) CYCLE
          iPart = iPart + 1
          PartState(1,iPart)   = PartData(1,offsetnPart+iLoop)
          PartState(2,iPart)   = PartData(2,offsetnPart+iLoop)
          PartState(3,iPart)   = PartData(3,offsetnPart+iLoop)
          PartState(4,iPart)   = PartData(4,offsetnPart+iLoop)
          PartState(5,iPart)   = PartData(5,offsetnPart+iLoop)
          PartState(6,iPart)   = PartData(6,offsetnPart+iLoop)
          PartSpecies(iPart)= INT(PartData(7,offsetnPart+iLoop),4)
          IF (useDSMC) THEN
            IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicModel)) THEN
              PartStateIntEn(1,iPart)=PartData(8,offsetnPart+iLoop)
              PartStateIntEn(2,iPart)=PartData(9,offsetnPart+iLoop)
              PartStateIntEn(3,iPart)=PartData(10,offsetnPart+iLoop)
              PartMPF(iPart)=PartData(11,offsetnPart+iLoop)
            ELSE IF ((CollisMode.GT.1).AND. (usevMPF)) THEN
              PartStateIntEn(1,iPart)=PartData(8,offsetnPart+iLoop)
              PartStateIntEn(2,iPart)=PartData(9,offsetnPart+iLoop)
              PartMPF(iPart)=PartData(10,offsetnPart+iLoop)
            ELSE IF ((CollisMode.GT.1).AND. (DSMC%ElectronicModel)) THEN
              PartStateIntEn(1,iPart)=PartData(8,offsetnPart+iLoop)
              PartStateIntEn(2,iPart)=PartData(9,offsetnPart+iLoop)
              PartStateIntEn(3,iPart)=PartData(10,offsetnPart+iLoop)
            ELSE IF (CollisMode.GT.1) THEN
              IF (readVarFromState(8).AND.readVarFromState(9)) THEN
                PartStateIntEn(1,iPart)=PartData(8,offsetnPart+iLoop)
                PartStateIntEn(2,iPart)=PartData(9,offsetnPart+iLoop)
              ELSE IF ((SpecDSMC(PartSpecies(iPart))%InterID.EQ.1).OR.&
                       (SpecDSMC(PartSpecies(iPart))%InterID.EQ.10).OR.&
                       (SpecDSMC(PartSpecies(iPart))%InterID.EQ.15)) THEN
                !- setting inner DOF to 0 for atoms
                PartStateIntEn(1,iPart)=0.
                PartStateIntEn(2,iPart)=0.
              ELSE
                CALL Abort(&
                    __STAMP__&
                    ,"resetting inner DOF for molecules is not implemented yet!"&
                ,SpecDSMC(PartSpecies(iPart))%InterID , PartData(7,offsetnPart+iLoop))
              END IF ! readVarFromState(8).AND.readVarFromState(9)
            ELSE IF (usevMPF) THEN
              PartMPF(iPart)=PartData(8,offsetnPart+iLoop)
            END IF ! (CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicModel)
          ELSE IF (usevMPF) THEN
            PartMPF(iPart)=PartData(8,offsetnPart+iLoop)
          END IF ! UseDSMC

          IF (useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
            IF (SpecDSMC(PartSpecies(iPart))%PolyatomicMol) THEN
              iPolyatMole = SpecDSMC(PartSpecies(iPart))%SpecToPolyArray
              IF(ALLOCATED(VibQuantsPar(iPart)%Quants)) DEALLOCATE(VibQuantsPar(iPart)%Quants)
              ALLOCATE(VibQuantsPar(iPart)%Quants(PolyatomMolDSMC(iPolyatMole)%VibDOF))
              VibQuantsPar(iPart)%Quants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF)= &
                  VibQuantData(1:PolyatomMolDSMC(iPolyatMole)%VibDOF,offsetnPart+iLoop)
            END IF ! SpecDSMC(PartSpecies(iPart))%PolyatomicMol
          END IF ! useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)

          PDM%ParticleInside(iPart) = .TRUE.
        END DO ! iLoop = 1_IK,locnPart
        iPart = 0
        DO iElem=FirstElemInd,LastElemInd
          IF (PartInt(iElem,ELEM_LastPartInd).GT.PartInt(iElem,ELEM_FirstPartInd)) THEN
            DO iLoop = PartInt(iElem,ELEM_FirstPartInd)-offsetnPart+1_IK , PartInt(iElem,ELEM_LastPartInd)- offsetnPart
              IF(SpecReset(INT(PartData(7,offsetnPart+iLoop),4))) CYCLE
              iPart = iPart +1
              PEM%Element(iPart)  = iElem-offsetElem
              PEM%LastElement(iPart)  = iElem-offsetElem
            END DO ! iLoop
          END IF ! PartInt(iElem,ELEM_LastPartInd).GT.PartInt(iElem,ELEM_FirstPartInd)
        END DO ! iElem=FirstElemInd,LastElemInd
        DEALLOCATE(PartData)
        IF (useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
          DEALLOCATE(VibQuantData)
        END IF
      ELSE ! not PartDataExists
        SWRITE(UNIT_stdOut,*)'PartData does not exists in restart file'
      END IF ! PartDataExists
      DEALLOCATE(PartInt)

      PDM%ParticleVecLength = PDM%ParticleVecLength + iPart
      CALL UpdateNextFreePosition()
      SWRITE(UNIT_stdOut,*)' DONE!'

      ! if ParticleVecLength GT maxParticleNumber: Stop
      IF (PDM%ParticleVecLength.GT.PDM%maxParticleNumber) THEN
        SWRITE (UNIT_stdOut,*) "PDM%ParticleVecLength =", PDM%ParticleVecLength
        SWRITE (UNIT_stdOut,*) "PDM%maxParticleNumber =", PDM%maxParticleNumber
        CALL abort(__STAMP__&
            ,' Number of Particles in Restart file is higher than MaxParticleNumber! Increase MaxParticleNumber!')
      END IF ! PDM%ParticleVecLength.GT.PDM%maxParticleNumber
      ! Since the elementside-local node number are NOT persistant and dependent on the location
      ! of the MPI borders, all particle-element mappings need to be checked after a restart
      ! Step 1: Identify particles that are not in the element in which they were before the restart
      COUNTER = 0
      COUNTER2 = 0
      CounterPoly = 0

      IF(DoRefMapping) THEN
        DO i = 1,PDM%ParticleVecLength
          CALL GetPositionInRefElem(PartState(1:3,i),Xi,PEM%Element(i))
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
              IF(useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
                IF(SpecDSMC(PartSpecies(i))%PolyatomicMol) THEN
                  iPolyatMole = SpecDSMC(PartSpecies(i))%SpecToPolyArray
                  CounterPoly = CounterPoly + PolyatomMolDSMC(iPolyatMole)%VibDOF
                END IF
              END IF
              PartPosRef(1:3,i) = -888.
            ELSE
              PEM%LastElement(i) = PEM%Element(i)
            END IF
          END IF
        END DO ! i = 1,PDM%ParticleVecLength
      ELSE ! no Ref Mapping
        IF (TriaTracking) THEN
          DO i = 1,PDM%ParticleVecLength
            CALL ParticleInsideQuad3D(PartState(1:3,i),PEM%Element(i),InElementCheck,det)
            IF (.NOT.InElementCheck) THEN  ! try to find them within MyProc
              COUNTER = COUNTER + 1
              CALL SingleParticleToExactElement(i,doHALO=.FALSE.,initFix=.FALSE.,doRelocate=.FALSE.)
              IF (.NOT.PDM%ParticleInside(i)) THEN
                COUNTER2 = COUNTER2 + 1
                IF(useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
                  IF(SpecDSMC(PartSpecies(i))%PolyatomicMol) THEN
                    iPolyatMole = SpecDSMC(PartSpecies(i))%SpecToPolyArray
                    CounterPoly = CounterPoly + PolyatomMolDSMC(iPolyatMole)%VibDOF
                  END IF
                END IF
              ELSE
                PEM%LastElement(i) = PEM%Element(i)
              END IF
            END IF
          END DO ! i = 1,PDM%ParticleVecLength
        ELSE ! not TriaTracking
          DO i = 1,PDM%ParticleVecLength
            CALL GetPositionInRefElem(PartState(1:3,i),Xi,PEM%Element(i))
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
                IF(useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
                  IF(SpecDSMC(PartSpecies(i))%PolyatomicMol) THEN
                    iPolyatMole = SpecDSMC(PartSpecies(i))%SpecToPolyArray
                    CounterPoly = CounterPoly + PolyatomMolDSMC(iPolyatMole)%VibDOF
                  END IF
                END IF
              ELSE
                PEM%LastElement(i) = PEM%Element(i)
              END IF ! .NOT.PDM%ParticleInside(i)
            END IF ! .NOT.InElementCheck
          END DO ! i = 1,PDM%ParticleVecLength
        END IF ! TriaTracking
      END IF ! DoRefMapping
#if USE_MPI
      ! Step 2: All particles that are not found withing MyProc need to be communicated to the others and located there
      ! Combine number of lost particles of all processes and allocate variables
      CALL MPI_ALLGATHER(COUNTER2, 1, MPI_INTEGER, LostParts, 1, MPI_INTEGER, PartMPI%COMM, IERROR)
      IF(useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) &
          CALL MPI_ALLGATHER(CounterPoly, 1, MPI_INTEGER, LostPartsPoly, 1, MPI_INTEGER, PartMPI%COMM, IERROR)
      IF (SUM(LostParts).GT.0) THEN
        ALLOCATE(SendBuff(1:COUNTER2*PartDataSize))
        ALLOCATE(RecBuff(1:SUM(LostParts)*PartDataSize))
        IF(useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
          ALLOCATE(SendBuffPoly(1:CounterPoly))
          ALLOCATE(RecBuffPoly(1:SUM(LostPartsPoly)))
        END IF
        ! Fill SendBuffer
        COUNTER = 0
        CounterPoly = 0
        DO i = 1, PDM%ParticleVecLength
          IF (.NOT.PDM%ParticleInside(i)) THEN
            SendBuff(COUNTER+1:COUNTER+6) = PartState(1:6,i)
            SendBuff(COUNTER+7)           = REAL(PartSpecies(i))
            IF (useDSMC) THEN
              IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicModel)) THEN
                SendBuff(COUNTER+8)  = PartStateIntEn(1,i)
                SendBuff(COUNTER+9)  = PartStateIntEn(2,i)
                SendBuff(COUNTER+10) = PartMPF(i)
                SendBuff(COUNTER+11) = PartStateIntEn(3,i)
              ELSE IF ((CollisMode.GT.1).AND. (usevMPF)) THEN
                SendBuff(COUNTER+8)  = PartStateIntEn(1,i)
                SendBuff(COUNTER+9)  = PartStateIntEn(2,i)
                SendBuff(COUNTER+10) = PartMPF(i)
              ELSE IF ((CollisMode.GT.1).AND. (DSMC%ElectronicModel)) THEN
                SendBuff(COUNTER+8)  = PartStateIntEn(1,i)
                SendBuff(COUNTER+9)  = PartStateIntEn(2,i)
                SendBuff(COUNTER+10) = PartStateIntEn(3,i)
              ELSE IF (CollisMode.GT.1) THEN
                SendBuff(COUNTER+8)  = PartStateIntEn(1,i)
                SendBuff(COUNTER+9)  = PartStateIntEn(2,i)
              ELSE IF (usevMPF) THEN
                SendBuff(COUNTER+8) = PartMPF(i)
              END IF
            ELSE IF (usevMPF) THEN
              SendBuff(COUNTER+8) = PartMPF(i)
            END IF
            COUNTER = COUNTER + PartDataSize

            !--- receive the polyatomic vibquants per particle at the end of the message
            IF(useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
              IF(SpecDSMC(PartSpecies(i))%PolyatomicMol) THEN
                iPolyatMole = SpecDSMC(PartSpecies(i))%SpecToPolyArray
                SendBuffPoly(CounterPoly+1:CounterPoly+PolyatomMolDSMC(iPolyatMole)%VibDOF) &
                    = VibQuantsPar(i)%Quants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF)
                CounterPoly = CounterPoly + PolyatomMolDSMC(iPolyatMole)%VibDOF
              END IF ! SpecDSMC(PartSpecies(i))%PolyatomicMol
            END IF ! useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)
          END IF ! .NOT.PDM%ParticleInside(i)
        END DO ! i = 1, PDM%ParticleVecLength
        ! Distribute lost particles to all procs
        COUNTER = 0
        CounterPoly = 0
        DO i = 0, PartMPI%nProcs-1
          RecCount(i) = LostParts(i) * PartDataSize
          Displace(i) = COUNTER
          COUNTER = COUNTER + LostParts(i)*PartDataSize
          IF(useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
            DisplacePoly(i) = CounterPoly
            CounterPoly = CounterPoly + LostPartsPoly(i)
          END IF
        END DO ! i = 0, PartMPI%nProcs-1
        CALL MPI_ALLGATHERV(SendBuff, PartDataSize*LostParts(PartMPI%MyRank), MPI_DOUBLE_PRECISION, &
            RecBuff, RecCount, Displace, MPI_DOUBLE_PRECISION, PartMPI%COMM, IERROR)
        IF(useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) CALL MPI_ALLGATHERV(SendBuffPoly, LostPartsPoly(PartMPI%MyRank), MPI_INTEGER, &
            RecBuffPoly, LostPartsPoly, DisplacePoly, MPI_INTEGER, PartMPI%COMM, IERROR)
        ! Add them to particle list and check if they are in MyProcs domain
        NbrOfFoundParts = 0
        CurrentPartNum = PDM%ParticleVecLength+1
        COUNTER = 0
        CounterPoly = 0
        DO i = 1, SUM(LostParts)
          PartState(1:6,CurrentPartNum) = RecBuff(COUNTER+1:COUNTER+6)
          PDM%ParticleInside(CurrentPartNum) = .true.
          IF(DoRefMapping.OR.TriaTracking)THEN
            CALL SingleParticleToExactElement(CurrentPartNum,doHALO=.FALSE.,initFix=.FALSE.,doRelocate=.FALSE.)
          ELSE
            CALL SingleParticleToExactElementNoMap(CurrentPartNum,doHALO=.FALSE.,doRelocate=.FALSE.)
          END IF
          IF (PDM%ParticleInside(CurrentPartNum)) THEN
            PEM%LastElement(CurrentPartNum) = PEM%Element(CurrentPartNum)
            NbrOfFoundParts = NbrOfFoundParts + 1
            PartSpecies(CurrentPartNum) = INT(RecBuff(COUNTER+7))
            IF (useDSMC) THEN
              IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicModel)) THEN
                PartStateIntEn(1,CurrentPartNum) = RecBuff(COUNTER+8)
                PartStateIntEn(2,CurrentPartNum) = RecBuff(COUNTER+9)
                PartStateIntEn(3,CurrentPartNum) = RecBuff(COUNTER+11)
                PartMPF(CurrentPartNum)          = RecBuff(COUNTER+10)
              ELSE IF ((CollisMode.GT.1).AND. (usevMPF)) THEN
                PartStateIntEn(1,CurrentPartNum) = RecBuff(COUNTER+8)
                PartStateIntEn(2,CurrentPartNum) = RecBuff(COUNTER+9)
                PartMPF(CurrentPartNum)          = RecBuff(COUNTER+10)
              ELSE IF ((CollisMode.GT.1).AND. (DSMC%ElectronicModel)) THEN
                PartStateIntEn(1,CurrentPartNum) = RecBuff(COUNTER+8)
                PartStateIntEn(2,CurrentPartNum) = RecBuff(COUNTER+9)
                PartStateIntEn(3,CurrentPartNum) = RecBuff(COUNTER+10)
              ELSE IF (CollisMode.GT.1) THEN
                PartStateIntEn(1,CurrentPartNum) = RecBuff(COUNTER+8)
                PartStateIntEn(2,CurrentPartNum) = RecBuff(COUNTER+9)
              ELSE IF (usevMPF) THEN
                PartMPF(CurrentPartNum)          = RecBuff(COUNTER+8)
              END IF
            ELSE IF (usevMPF) THEN
              PartMPF(CurrentPartNum)          = RecBuff(COUNTER+8)
            END IF
            IF(useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
              IF(SpecDSMC(PartSpecies(CurrentPartNum))%PolyatomicMol) THEN
                iPolyatMole = SpecDSMC(PartSpecies(CurrentPartNum))%SpecToPolyArray
                IF(ALLOCATED(VibQuantsPar(CurrentPartNum)%Quants)) DEALLOCATE(VibQuantsPar(CurrentPartNum)%Quants)
                ALLOCATE(VibQuantsPar(CurrentPartNum)%Quants(PolyatomMolDSMC(iPolyatMole)%VibDOF))
                VibQuantsPar(CurrentPartNum)%Quants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF) &
                    = RecBuffPoly(CounterPoly+1:CounterPoly+PolyatomMolDSMC(iPolyatMole)%VibDOF)
                CounterPoly = CounterPoly + PolyatomMolDSMC(iPolyatMole)%VibDOF
              END IF
            END IF
            CurrentPartNum = CurrentPartNum + 1
          END IF
          COUNTER = COUNTER + PartDataSize
        END DO ! i = 1, SUM(LostParts)
        PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfFoundParts
        ! Combine number of found particles to make sure none are lost completely
        CALL MPI_ALLREDUCE(NbrOfFoundParts, CompleteNbrOfFound, 1, MPI_INTEGER, MPI_SUM, PartMPI%COMM, IERROR)
        SWRITE(UNIT_stdOut,*) SUM(LostParts),'were not in the correct proc after restart.'
        SWRITE(UNIT_stdOut,*) CompleteNbrOfFound,'of these were found in other procs.'
        SWRITE(UNIT_stdOut,*) SUM(LostParts)-CompleteNbrOfFound,'were not found and have been removed.'
      END IF ! SUM(LostParts).GT.0
#else
      IF (COUNTER.NE.0) WRITE(*,*) COUNTER,'Particles are in different element after restart!'
      IF (COUNTER2.NE.0) WRITE(*,*) COUNTER2,'of which could not be found and are removed!'
#endif

      CALL UpdateNextFreePosition()

      IF (RadialWeighting%PerformCloning) THEN
        CALL DatasetExists(File_ID,'CloneData',CloneExists)
        IF(CloneExists) THEN
          CALL RestartClones()
        ELSE
          SWRITE(*,*) 'No clone data found! Restart without cloning.'
          IF(RadialWeighting%CloneMode.EQ.1) THEN
            RadialWeighting%CloneDelayDiff = 1
          ELSEIF (RadialWeighting%CloneMode.EQ.2) THEN
            RadialWeighting%CloneDelayDiff = 0
          END IF ! RadialWeighting%CloneMode.EQ.1
        END IF ! CloneExists
      END IF ! RadialWeighting%PerformCloning
    ELSE ! not PartIntExists
      SWRITE(UNIT_stdOut,*)'PartInt does not exists in restart file'
    END IF ! PartIntExists
  ELSE ! DoMacroscopicRestart
    CALL CloseDataFile()
    CALL MacroscopicRestart()
    CALL UpdateNextFreePosition()
  END IF ! .NOT.DoMacroscopicRestart

  ! ------------------------------------------------
  ! Reactive Surfaces
  ! ------------------------------------------------
  WallModelExists(:)=.FALSE.
  IF (ANY(PartBound%Reactive)) THEN
    ! check if datasets for restarting of surface model from state exists in state file used for restart
    SurfModelTypeExists=.FALSE.
    CALL DatasetExists(File_ID,'SurfaceModelType',SurfModelTypeExists)
    IF (SurfModelTypeExists) THEN
      SWRITE(UNIT_stdOut,'(A,A)')' GET NUMBER OF SURFACE-SIDES IN RESTART FILE... '
      CALL GetDataSize(File_ID,'SurfaceModelType',nDims,HSize,attrib=.FALSE.)
      nSurfSides_HDF5 = INT(HSize(1),4)
      IF (nSurfSides_HDF5.NE.SurfMesh%nGlobalSides) THEN
        SWRITE(UNIT_stdOut,'(A,A)') ' NUMBER OF SURFACE-SIDES IN RESTART FILE NOT EQUAL TO CALCULATION ... RESTARTING FROM INI'
      ELSE
        ALLOCATE(SurfModelType(SurfMesh%nOutputSides))
        ! Associate construct for integer KIND=8 possibility
        ASSOCIATE (&
              nSides          => INT(SurfMesh%nOutputSides,IK) ,&
              offsetSurfSide  => INT(offsetSurfSide,IK) )
          CALL ReadArray('SurfaceModelType',1,(/nSides/) ,&
              offsetSurfSide,1,IntegerArray_i4=SurfModelType)
        END ASSOCIATE
        WallModelExists(:)=.TRUE.
        DO iSurfSide=1,SurfMesh%nOutputSides
          SideID = SurfMesh%SurfIDToSideID(iSurfSide)
          PartboundID = PartBound%MapToPartBC(BC(SideID))
          IF (PartBound%SurfaceModel(PartboundID).NE.SurfModelType(iSurfSide)) THEN
            WallModelExists(PartBoundID)=.FALSE.
            EXIT
          END IF
        END DO
      END IF ! nsurfsides_hdf5 != nglobalsurfsides
    END IF
    IF (ANY(WallModelExists(:))) THEN
      SWRITE(UNIT_stdOut,*)'Reading surface calculation infos from Restartfile...'
      ! do sanity checks of data in h5 file before proceeding
      CALL GetDataSize(File_ID,'Surface_BCs',nDims,HSize,attrib=.TRUE.)
      nSurfBC_HDF5 = INT(HSize(1),4)
      IF (nSurfBC_HDF5.NE.nSurfBC) CALL abort(&
          __STAMP__&
          ,'Error in surface restart: number of surface boundaries in HDF5-file does not match!')
      CALL ReadAttribute(File_ID,'nSurfSample',1,IntegerScalar=nSurfSample_HDF5)
      IF (nSurfSample_HDF5.NE.nSurfSample) CALL abort(&
          __STAMP__&
          ,'Error in surface restart: number of surface subsides (nSurfSample) in HDF5-file does not match!')
      CALL ReadAttribute(File_ID,'nSpecies',1,IntegerScalar=nSpecies_HDF5)
      IF (nSpecies_HDF5.NE.nSpecies) CALL abort(&
          __STAMP__&
          ,'Error in surface restart: number of Species in HDF5-file does not match!')

      SurfCalcDataExists=.FALSE.
      CALL DatasetExists(File_ID,'SurfCalcData',SurfCalcDataExists)
      IF (SurfCalcDataExists) THEN
        nVar = 4
        ALLOCATE(SurfCalcData(nVar,nSurfSample,nSurfSample,SurfMesh%nOutputSides,nSpecies))

        ! Associate construct for integer KIND=8 possibility
        ASSOCIATE (&
              nVar            => INT(nVar,IK) ,&
              nSurfSample     => INT(nSurfSample,IK) ,&
              nSides          => INT(SurfMesh%nOutputSides,IK) ,&
              nSpecies        => INT(nSpecies,IK) ,&
              offsetSurfSide  => INT(offsetSurfSide,IK) )
          CALL ReadArray('SurfCalcData',5,(/nVar,nSurfSample,nSurfSample,nSides,nSpecies/) ,&
              offsetSurfSide,4,RealArray=SurfCalcData)
        END ASSOCIATE
        DO iSurfSide = 1,SurfMesh%nOutputSides
          SideID = SurfMesh%SurfIDToSideID(iSurfSide)
          PartboundID = PartBound%MapToPartBC(BC(SideID))
          IF (PartBound%Reactive(PartboundID).AND.WallModelExists(PartBoundID)) THEN
            DO jsubsurf = 1,nSurfSample
              DO isubsurf = 1,nSurfSample
                Adsorption%Coverage(iSubSurf,jSubSurf,iSurfSide,:) = SurfCalcData(1,iSubSurf,jSubSurf,iSurfSide,:)
                IF (PartBound%SurfaceModel(PartBoundID).EQ.3) THEN
                  SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%adsorbnum_tmp(:) = SurfCalcData(2,iSubSurf,jSubSurf,iSurfSide,:)
                  SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%desorbnum_tmp(:) = SurfCalcData(3,iSubSurf,jSubSurf,iSurfSide,:)
                  SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%reactnum_tmp(:)  = SurfCalcData(4,iSubSurf,jSubSurf,iSurfSide,:)
                END IF
              END DO
            END DO
          END IF
        END DO
        DEALLOCATE(SurfCalcData)
        ! read additional data for wallmodel 3
        !IF (PartSurfaceModel.EQ.3) THEN
        Coordinations    = 3
        SurfPartIntSize  = 3
        SurfPartDataSize = 2
        ! check if surfpartint exists
        SurfPartIntExists=.FALSE.
        CALL DatasetExists(File_ID,'SurfPartInt',SurfPartIntExists)
        IF(SurfPartIntExists)THEN
          ALLOCATE(SurfPartInt(offsetSurfSide+1:offsetSurfSide+SurfMesh%nOutputSides &
              ,nSurfSample,nSurfSample,Coordinations,SurfPartIntSize))

          ! Associate construct for integer KIND=8 possibility
          ASSOCIATE (&
                nSides          => INT(SurfMesh%nOutputSides,IK) ,&
                nSurfSample     => INT(nSurfSample,IK)     ,&
                Coordinations   => INT(Coordinations,IK)   ,&
                SurfPartIntSize => INT(SurfPartIntSize,IK) ,&
                offsetSurfSide  => INT(offsetSurfSide,IK)   )
            ! read local Surface Particle indexing from HDF5
            CALL ReadArray('SurfPartInt',5,(/nSides,nSurfSample,nSurfSample,Coordinations,SurfPartIntSize/) &
                ,offsetSurfSide,1,IntegerArray_i4=SurfPartInt)
          END ASSOCIATE
          ! check if surfpartdata exists
          SurfPartDataExists=.FALSE.
          CALL DatasetExists(File_ID,'SurfPartData',SurfPartDataExists)
          IF(SurfPartDataExists)THEN
            IF (SurfMesh%nOutputSides.GT.0) THEN
              locnSurfPart = SurfPartInt(offsetSurfSide+SurfMesh%nOutputSides,nSurfSample,nSurfSample,Coordinations,3) &
                  - SurfPartInt(offsetSurfSide+1,1,1,1,2)
              offsetnSurfPart=SurfPartInt(offsetSurfSide+1,1,1,1,2)
            ELSE
              locnSurfPart = 0
              offsetnSurfPart = 0
            END IF
            ALLOCATE(SurfPartData(SurfPartDataSize,offsetnSurfPart+1:offsetnSurfPart+locnSurfPart))
            ! read local Surface Particle Data from HDF5

            ! Associate construct for integer KIND=8 possibility
            ASSOCIATE (&
                  locnSurfPart      => INT(locnSurfPart,IK)      ,&
                  SurfPartDataSize  => INT(SurfPartDataSize,IK)  ,&
                  offsetnSurfPart   => INT(offsetnSurfPart,IK)   )
              CALL ReadArray('SurfPartData',2,(/SurfPartDataSize,locnSurfPart/),offsetnSurfPart,2,IntegerArray_i4=SurfPartData)
            END ASSOCIATE
            IF (locnSurfPart.GT.0) THEN
              DO iSurfSide = 1,SurfMesh%nOutputSides
                SideID = SurfMesh%SurfIDToSideID(iSurfSide)
                PartboundID = PartBound%MapToPartBC(BC(SideID))
                IF (WallModelExists(PartBoundID).AND.PartBound%SurfaceModel(PartboundID).EQ.3) THEN
                  DO jsubsurf = 1,nSurfSample
                    DO isubsurf = 1,nSurfSample
                      DO iCoord = 1,Coordinations
                        firstpart = SurfPartInt(offsetSurfSide+iSurfSide,isubsurf,jsubsurf,iCoord,2) + 1
                        lastpart  = SurfPartInt(offsetSurfSide+iSurfSide,isubsurf,jsubsurf,iCoord,3)
                        ! set the surfpartdata array values
                        DO iPart = firstpart, lastpart
                          UsedSiteMapPos = SurfPartData(1,iPart)
                          SpecID         = SurfPartData(2,ipart)
                          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(iCoord)%Species(UsedSiteMapPos) = SpecID
                          ! assign bond order of respective surface atoms in the surface lattice
                          DO iInterAtom = 1,SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(iCoord)%nInterAtom
                            xpos = SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(iCoord)%BondAtomIndx( &
                                UsedSiteMapPos,iInterAtom)
                            ypos = SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(iCoord)%BondAtomIndy( &
                                UsedSiteMapPos,iInterAtom)
                            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%SurfAtomBondOrder(SpecID,xpos,ypos) = &
                                SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%SurfAtomBondOrder(SpecID,xpos,ypos) + 1
                          END DO
                        END DO ! iPart = firstpart,lastpart
                        ! sort and rearrange UsedSiteMap-Surfpos-array
                        ! structure of UsedSiteMap array for one coordination
                        !               [<---------------nSites---------------------------------------->]
                        ! Name        :  nfreeArrayindeces   (lastfreeIndx)                   Adsorbates
                        ! current     :  1 2                   3           |      4  5  6  7  8  9 10 11
                        ! UsedSiteMap :  1 7                   8           |     11  9 10  3  4  5  2  6
                        lastfreeIndx = SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(iCoord)
                        nfreeArrayindeces = lastfreeIndx - ( lastpart - (firstpart-1) )
                        DO current = 1,SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(iCoord)
                          UsedSiteMapPos =  SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(iCoord)%UsedSiteMap(current)
                          IF (SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(iCoord)%Species(UsedSiteMapPos).GT.0) THEN
                            MoveToLastFree = .TRUE.
                            ! move value to end of array and end of array to current array index
                            DO WHILE (MoveToLastFree)
                              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(iCoord)%UsedSiteMap(current) = &
                                  SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(iCoord)%UsedSiteMap(lastfreeIndx)
                              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(iCoord)%UsedSiteMap(lastfreeIndx) = &
                                  UsedSiteMapPos
                              UsedSiteMapPos =  SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(iCoord)%UsedSiteMap(current)
                              IF (SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(iCoord)%Species(UsedSiteMapPos).EQ.0) THEN
                                MoveToLastFree = .FALSE.
                              END IF
                              lastfreeIndx = lastfreeIndx - 1
                              IF (lastfreeIndx .EQ. nfreeArrayindeces) EXIT
                            END DO ! current Indx not empty
                          END IF
                          IF (lastfreeIndx .EQ. nfreeArrayindeces) EXIT
                        END DO ! current = 1,nSites
                        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%SitesRemain(iCoord) = nfreeArrayindeces
                      END DO
                    END DO
                  END DO
                END IF
              END DO
              DEALLOCATE(SurfPartData)
            END IF
          END IF ! SurfPartDataExists
          DEALLOCATE(SurfPartInt)
        END IF ! SurfPartIntExists
        !END IF ! PartSurfaceModel.EQ.3
      END IF ! SurfCalcDataExists
    ELSE
      SWRITE(UNIT_stdOut,*)'Data for current wallmodel does not exists in restart file'
    END IF ! WallModel_HDF5.NE.PartSurfaceModel
  END IF ! ANY(PartBound%Reactive)
#endif /*PARTICLES*/

CALL CloseDataFile()

#ifdef PARTICLES
  ! include initially collected particles for first call of field-solver (e.g. in RecomputeLambda)
  CALL ParticleCollectCharges(initialCall_opt=.TRUE.)
#endif /*PARTICLES*/
#if USE_HDG
  iter=0
  ! INSTEAD OF ALL THIS STUFF DO
  ! 1) MPI-Communication for shape-function particles
  ! 2) Deposition
  ! 3) ONE HDG solve
  CALL  RecomputeLambda(RestartTime)
#endif /*USE_HDG*/


  ! Delete all files that will be rewritten
  IF(DoWriteStateToHDF5) CALL FlushHDF5(RestartTime)
#if USE_MPI
  EndT=MPI_WTIME()
  SWRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')' Restart took  [',EndT-StartT,'s] for readin.'
  SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' Restart DONE!'
#else
  SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' Restart DONE!'
#endif
ELSE ! no restart
#ifdef PARTICLES
  ! include initially collected particles for first call of field-solver (here because of consistency, but not used until timedisc)
  CALL ParticleCollectCharges(initialCall_opt=.TRUE.)
#endif /*PARTICLES*/
  ! Delete all files since we are doing a fresh start
  IF(DoWriteStateToHDF5) CALL FlushHDF5()
END IF !IF(DoRestart)
END SUBROUTINE Restart

#ifdef PARTICLES
SUBROUTINE RestartClones()
!===================================================================================================================================
! Axisymmetric 2D simulation with particle weighting: Read-in of clone particles saved during output of particle data
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_HDF5_input
USE MOD_io_hdf5
USE MOD_Mesh_Vars,                ONLY : offsetElem, nElems
USE MOD_DSMC_Vars,                ONLY : UseDSMC, CollisMode, DSMC, PolyatomMolDSMC, SpecDSMC
USE MOD_DSMC_Vars,                ONLY : RadialWeighting, ClonedParticles
USE MOD_Particle_Vars,            ONLY : nSpecies, usevMPF
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                           :: nDimsClone, CloneDataSize, ClonePartNum, iPart, iDelay, maxDelay, iElem, tempDelay
  INTEGER(HSIZE_T), POINTER         :: SizeClone(:)
  REAL,ALLOCATABLE                  :: CloneData(:,:)
  INTEGER                           :: iPolyatmole, MaxQuantNum, iSpec, compareDelay
  INTEGER,ALLOCATABLE               :: pcount(:), VibQuantData(:,:)
!===================================================================================================================================

  CALL GetDataSize(File_ID,'CloneData',nDimsClone,SizeClone)

  CloneDataSize = INT(SizeClone(1),4)
  ClonePartNum = INT(SizeClone(2),4)
  DEALLOCATE(SizeClone)

  IF(ClonePartNum.GT.0) THEN
    ALLOCATE(CloneData(1:CloneDataSize,1:ClonePartNum))
    ASSOCIATE(ClonePartNum  => INT(ClonePartNum,IK)  ,&
              CloneDataSize => INT(CloneDataSize,IK) )
      CALL ReadArray('CloneData',2,(/CloneDataSize,ClonePartNum/),0_IK,2,RealArray=CloneData)
    END ASSOCIATE
    SWRITE(*,*) 'Read-in of cloned particles complete. Total clone number: ', ClonePartNum
    ! Determing the old clone delay
    maxDelay = INT(MAXVAL(CloneData(9,:)))
    IF(RadialWeighting%CloneMode.EQ.1) THEN
      ! Array is allocated from 0 to maxDelay
      compareDelay = maxDelay + 1
    ELSE
      compareDelay = maxDelay
    END IF
    IF(compareDelay.GT.RadialWeighting%CloneInputDelay) THEN
      SWRITE(*,*) 'Old clone delay is greater than the new delay. Old delay:', compareDelay
      RadialWeighting%CloneDelayDiff = RadialWeighting%CloneInputDelay + 1
    ELSEIF(compareDelay.EQ.RadialWeighting%CloneInputDelay) THEN
      SWRITE(*,*) 'The clone delay has not been changed.'
      RadialWeighting%CloneDelayDiff = RadialWeighting%CloneInputDelay + 1
    ELSE
      SWRITE(*,*) 'New clone delay is greater than the old delay. Old delay:', compareDelay
      RadialWeighting%CloneDelayDiff = compareDelay + 1
    END IF
    IF(RadialWeighting%CloneMode.EQ.1) THEN
      tempDelay = RadialWeighting%CloneInputDelay - 1
    ELSE
      tempDelay = RadialWeighting%CloneInputDelay
    END IF
    ALLOCATE(pcount(0:tempDelay))
    pcount(0:tempDelay) = 0
    ! Polyatomic clones: determining the size of the VibQuant array
    IF (UseDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
      MaxQuantNum = 0
      DO iSpec = 1, nSpecies
        IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
          iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
          IF (PolyatomMolDSMC(iPolyatMole)%VibDOF.GT.MaxQuantNum) MaxQuantNum = PolyatomMolDSMC(iPolyatMole)%VibDOF
        END IF
      END DO
      ALLOCATE(VibQuantData(1:MaxQuantNum,1:ClonePartNum))
      ASSOCIATE(ClonePartNum => INT(ClonePartNum,IK),MaxQuantNum => INT(MaxQuantNum,IK))
        CALL ReadArray('CloneVibQuantData',2,(/MaxQuantNum,ClonePartNum/),0_IK,2,IntegerArray_i4=VibQuantData)
      END ASSOCIATE
    END IF
    ! Copying particles into ClonedParticles array
    DO iPart = 1, ClonePartNum
      iDelay = INT(CloneData(9,iPart))
      iElem = INT(CloneData(8,iPart)) - offsetElem
      IF((iElem.LE.nElems).AND.(iElem.GT.0)) THEN
        IF(iDelay.LE.tempDelay) THEN
          pcount(iDelay) = pcount(iDelay) + 1
          RadialWeighting%ClonePartNum(iDelay) = pcount(iDelay)
          ClonedParticles(pcount(iDelay),iDelay)%PartState(1) = CloneData(1,iPart)
          ClonedParticles(pcount(iDelay),iDelay)%PartState(2) = CloneData(2,iPart)
          ClonedParticles(pcount(iDelay),iDelay)%PartState(3) = CloneData(3,iPart)
          ClonedParticles(pcount(iDelay),iDelay)%PartState(4) = CloneData(4,iPart)
          ClonedParticles(pcount(iDelay),iDelay)%PartState(5) = CloneData(5,iPart)
          ClonedParticles(pcount(iDelay),iDelay)%PartState(6) = CloneData(6,iPart)
          ClonedParticles(pcount(iDelay),iDelay)%Species = INT(CloneData(7,iPart))
          ClonedParticles(pcount(iDelay),iDelay)%Element = iElem
          ClonedParticles(pcount(iDelay),iDelay)%lastPartPos(1:3) = CloneData(1:3,iPart)
          IF (UseDSMC) THEN
            IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicModel) ) THEN
              ClonedParticles(pcount(iDelay),iDelay)%PartStateIntEn(1) = CloneData(10,iPart)
              ClonedParticles(pcount(iDelay),iDelay)%PartStateIntEn(2) = CloneData(11,iPart)
              ClonedParticles(pcount(iDelay),iDelay)%PartStateIntEn(3) = CloneData(12,iPart)
              ClonedParticles(pcount(iDelay),iDelay)%WeightingFactor   = CloneData(13,iPart)
            ELSE IF ( (CollisMode .GT. 1) .AND. (usevMPF) ) THEN
              ClonedParticles(pcount(iDelay),iDelay)%PartStateIntEn(1) = CloneData(10,iPart)
              ClonedParticles(pcount(iDelay),iDelay)%PartStateIntEn(2) = CloneData(11,iPart)
              ClonedParticles(pcount(iDelay),iDelay)%WeightingFactor   = CloneData(12,iPart)
            ELSE IF ( (CollisMode .GT. 1) .AND. (DSMC%ElectronicModel) ) THEN
              ClonedParticles(pcount(iDelay),iDelay)%PartStateIntEn(1) = CloneData(10,iPart)
              ClonedParticles(pcount(iDelay),iDelay)%PartStateIntEn(2) = CloneData(11,iPart)
              ClonedParticles(pcount(iDelay),iDelay)%PartStateIntEn(3) = CloneData(12,iPart)
            ELSE IF (CollisMode.GT.1) THEN
              ClonedParticles(pcount(iDelay),iDelay)%PartStateIntEn(1) = CloneData(10,iPart)
              ClonedParticles(pcount(iDelay),iDelay)%PartStateIntEn(2) = CloneData(11,iPart)
            ELSE IF (usevMPF) THEN
              ClonedParticles(pcount(iDelay),iDelay)%WeightingFactor = CloneData(10,iPart)
            END IF
          ELSE IF (usevMPF) THEN
              ClonedParticles(pcount(iDelay),iDelay)%WeightingFactor = CloneData(10,iPart)
          END IF
          IF (UseDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
            IF (SpecDSMC(ClonedParticles(pcount(iDelay),iDelay)%Species)%PolyatomicMol) THEN
              iPolyatMole = SpecDSMC(ClonedParticles(pcount(iDelay),iDelay)%Species)%SpecToPolyArray
              ALLOCATE(ClonedParticles(pcount(iDelay),iDelay)%VibQuants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF))
              ClonedParticles(pcount(iDelay),iDelay)%VibQuants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF) &
                = VibQuantData(1:PolyatomMolDSMC(iPolyatMole)%VibDOF,iPart)
            ELSE
               VibQuantData(:,iPart) = 0
            END IF
          END IF
        END IF
      END IF
    END DO
  ELSE
    SWRITE(*,*) 'Read-in of cloned particles complete. No clones detected.'
  END IF

END SUBROUTINE RestartClones


SUBROUTINE MacroscopicRestart()
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_io_hdf5
USE MOD_HDF5_Input    ,ONLY: OpenDataFile,CloseDataFile,ReadArray,GetDataSize
USE MOD_HDF5_Input    ,ONLY: nDims,HSize,File_ID
USE MOD_Restart_Vars  ,ONLY: MacroRestartFileName, MacroRestartValues
USE MOD_Mesh_Vars     ,ONLY: offsetElem, nElems
USE MOD_Particle_Vars ,ONLY: nSpecies
USE MOD_macro_restart ,ONLY: MacroRestart_InsertParticles
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: nVar_HDF5, iVar, iSpec, iElem
REAL, ALLOCATABLE                 :: ElemData_HDF5(:,:)
!===================================================================================================================================

SWRITE(UNIT_stdOut,*) 'Using macroscopic values from file: ',TRIM(MacroRestartFileName)

CALL OpenDataFile(MacroRestartFileName,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)

CALL GetDataSize(File_ID,'ElemData',nDims,HSize,attrib=.FALSE.)
nVar_HDF5=INT(HSize(1),4)

ALLOCATE(MacroRestartValues(1:nElems,1:nSpecies+1,1:DSMC_NVARS))
MacroRestartValues = 0.

ALLOCATE(ElemData_HDF5(1:nVar_HDF5,1:nElems))
! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
  nVar_HDF5  => INT(nVar_HDF5,IK) ,&
  offsetElem => INT(offsetElem,IK),&
  nElems     => INT(nElems,IK)    )
  CALL ReadArray('ElemData',2,(/nVar_HDF5,nElems/),offsetElem,2,RealArray=ElemData_HDF5(:,:))
END ASSOCIATE

iVar = 1
DO iSpec = 1, nSpecies
  DO iElem = 1, nElems
    MacroRestartValues(iElem,iSpec,:) = ElemData_HDF5(iVar:iVar-1+DSMC_NVARS,iElem)
  END DO
  iVar = iVar + DSMC_NVARS
END DO

CALL MacroRestart_InsertParticles()

DEALLOCATE(MacroRestartValues)
DEALLOCATE(ElemData_HDF5)

END SUBROUTINE MacroscopicRestart
#endif /*PARTICLES*/


#if USE_HDG
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
#if USE_MPI
USE MOD_Particle_MPI,            ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY: PartMPIExchange,DoExternalParts
#endif /*USE_MPI*/
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
#if USE_MPI
IF(DoExternalParts)THEN
  ! communication of shape-function particles, YEAH.
  CALL IRecvNbofParticles() ! open receive buffer for number of particles
  CALL SendNbOfParticles() ! send number of particles
  CALL MPIParticleSend()  ! finish communication of number of particles and send particles
  CALL MPIParticleRecv()  ! finish communication
END IF
#endif

! Deposition of particles
CALL Deposition(DoInnerParts=.TRUE.)
#if USE_MPI
! here: finish deposition with delta kernal
!       maps source terms in physical space
! ALWAYS require
PartMPIExchange%nMPIParticles=0
#endif /*USE_MPI*/
CALL Deposition(DoInnerParts=.FALSE.)
#endif /*PARTICLES*/

! recompute fields
! EM field
CALL HDG(t,U,iter)

END SUBROUTINE RecomputeLambda
#endif /*USE_HDG*/

SUBROUTINE FinalizeRestart()
!===================================================================================================================================
! Finalizes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Restart_Vars,ONLY:Vdm_GaussNRestart_GaussN,RestartInitIsDone,DoMacroscopicRestart
! IMPLICIT VARIABLE HANDLINGDGInitIsDone
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(Vdm_GaussNRestart_GaussN)
RestartInitIsDone = .FALSE.
! Avoid performing a macroscopic restart during an automatic load balance restart
DoMacroscopicRestart = .FALSE.
END SUBROUTINE FinalizeRestart

END MODULE MOD_Restart
