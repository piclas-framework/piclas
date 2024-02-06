!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (piclas.boltzplatz.eu/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
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

PUBLIC :: InitRestart
PUBLIC :: Restart
PUBLIC :: FinalizeRestart

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
                              'Utilize a macroscopic result to restart the simulation. Particles are inserted based on the '//&
                              'cell local species-specific number density, temperature and velocity from a DSMCState file.', &
                              '.FALSE.')
CALL prms%CreateStringOption( 'Particles-MacroscopicRestart-Filename', &
                              'File name of the DSMCState to be utilized as the input for the particle insertion.')
CALL prms%CreateLogicalOption('FlushInitialState',&
                              'Check whether (during restart) the statefile from which the restart is performed should be deleted.'&
                            , '.FALSE.')
END SUBROUTINE DefineParametersRestart


SUBROUTINE InitRestart()
!===================================================================================================================================
! Initialize all necessary information to perform the restart
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: FileVersionHDF5Real
USE MOD_PreProc
USE MOD_ReadInTools            ,ONLY: GETLOGICAL,GETSTR
#if USE_LOADBALANCE
USE MOD_ReadInTools            ,ONLY: GETINT
USE MOD_LoadBalance_Vars       ,ONLY: LoadBalanceSample
USE MOD_ReadInTools            ,ONLY: PrintOption
#endif /*USE_LOADBALANCE*/
USE MOD_Interpolation_Vars     ,ONLY: xGP,InterpolationInitIsDone
USE MOD_Restart_Vars
USE MOD_HDF5_Input             ,ONLY: OpenDataFile,CloseDataFile,GetDataProps,ReadAttribute,File_ID
USE MOD_HDF5_Input             ,ONLY: DatasetExists
#ifdef PARTICLES
USE MOD_Particle_Tracking_Vars ,ONLY: TotalNbrOfMissingParticlesSum
#if USE_HDG
USE MOD_Part_BR_Elecron_Fluid  ,ONLY: InitSwitchBRElectronModel
#endif /*USE_HDG*/
#endif /*PARTICLES*/
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
INTEGER                     :: FileVersionHDF5Int
!===================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.RestartInitIsDone)THEN
   CALL abort(__STAMP__,'InitRestart not ready to be called or already called.',999,999.)
   RETURN
END IF

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT RESTART...'

#ifdef PARTICLES
! Set counter for particle that go missing during restart (if they are not located within their host element during restart)
TotalNbrOfMissingParticlesSum = 0
#if USE_HDG
  ! Initialize variables (only once, never during load balance restart) for switching between BR electron fluid model and fully
  ! kinetic model in HDG simulations
  CALL InitSwitchBRElectronModel()
#endif /*USE_HDG*/
#endif /*PARTICLES*/

! Set the DG solution to zero (ignore the DG solution in the state file)
RestartNullifySolution = GETLOGICAL('RestartNullifySolution','F')

! Check if we want to perform a restart
IF (LEN_TRIM(RestartFile).GT.0) THEN
  ! Read in the state file we want to restart from
  DoRestart = .TRUE.
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)
#ifdef PP_POIS
#if (PP_nVar==8)
  !The following arrays are read from the file
  !CALL ReadArray('DG_SolutionE',5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),OffsetElem,5,RealArray=U)
  !CALL ReadArray('DG_SolutionPhi',5,(/4,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),OffsetElem,5,RealArray=Phi)
  CALL abort(__STAMP__,'InitRestart: This case is not implemented here. Fix this!')
#else
  !The following arrays are read from the file
  !CALL ReadArray('DG_SolutionE',5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),OffsetElem,5,RealArray=U)
  !CALL ReadArray('DG_SolutionPhi',5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),OffsetElem,5,RealArray=Phi)
  CALL abort(__STAMP__,'InitRestart: This case is not implemented here. Fix this!')
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
#if USE_HDG
      SWRITE(UNIT_StdOut,'(A)')" HDG: Restarting from a different equation system."
#else
      CALL abort(__STAMP__,'PP_nVar!=nVar_Restart (Number of variables in restart file does no match the compiled equation system).')
#endif /*USE_HDG*/
    END IF
  END IF
  ! Read in time from restart file
  CALL ReadAttribute(File_ID,'Time',1,RealScalar=RestartTime)
  ! check file version
  CALL DatasetExists(File_ID,'File_Version',FileVersionExists,attrib=.TRUE.)
  IF (FileVersionExists) THEN
    CALL ReadAttribute(File_ID,'File_Version',1,RealScalar=FileVersionHDF5Real)

    IF(FileVersionHDF5Real.LT.1.5)THEN
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
      CALL abort(__STAMP__,&
      'Error in InitRestart(): "File_Version" in restart file < 1.5. See error message above to fix. File version in restart file =',&
      RealInfoOpt=FileVersionHDF5Real)
    END IF ! FileVersionHDF5Real.LT.1.5
  ELSE
    CALL DatasetExists(File_ID,'Piclas_VersionInt',FileVersionExists,attrib=.TRUE.)
    IF (FileVersionExists) THEN
      CALL ReadAttribute(File_ID,'Piclas_VersionInt',1,IntScalar=FileVersionHDF5Int)
    ELSE
      CALL abort(__STAMP__,'Error in InitRestart(): Attribute "Piclas_VersionInt" does not exist!')
    END IF
  END IF
  CALL CloseDataFile()
ELSE
  RestartTime = 0.
  SWRITE(UNIT_StdOut,'(A)')' | No restart wanted, doing a fresh computation!'
END IF

! Macroscopic restart
IF(DoRestart) THEN
  DoMacroscopicRestart = GETLOGICAL('Particles-MacroscopicRestart')
  IF(DoMacroscopicRestart) MacroRestartFileName = GETSTR('Particles-MacroscopicRestart-Filename')
END IF

! Automatically do a load balance step at the beginning of a new simulation or a user-restarted simulation
#if USE_LOADBALANCE
DoInitialAutoRestart = GETLOGICAL('DoInitialAutoRestart')
IF(nProcessors.LT.2) DoInitialAutoRestart = .FALSE.
WRITE(UNIT=hilf,FMT='(I0)') LoadBalanceSample
InitialAutoRestartSample = GETINT('InitialAutoRestartSample',TRIM(hilf))
InitialAutoRestartPartWeight = GETLOGICAL('InitialAutoRestart-PartWeightLoadBalance','F')
IF (InitialAutoRestartPartWeight) THEN
  InitialAutoRestartSample = 0 ! deactivate loadbalance sampling of ElemTimes if balancing with partweight is enabled
  CALL PrintOption('InitialAutoRestart-PartWeightLoadBalance = T : InitialAutoRestartSample','INFO',IntOpt=InitialAutoRestartSample)
ELSE IF (InitialAutoRestartSample.EQ.0) THEN
  InitialAutoRestartPartWeight = .TRUE. ! loadbalance (ElemTimes) is done with partmpiweight if loadbalancesampling is set to zero
  CALL PrintOption('InitialAutoRestart-PartWeightLoadBalance','INFO',LogOpt=InitialAutoRestartPartWeight)
END IF
#endif /*USE_LOADBALANCE*/

IF(DoRestart .AND. (N_Restart .NE. PP_N))THEN
  InterpolateSolution=.TRUE.
END IF

IF(InterpolateSolution)THEN
  CALL initRestartBasis(PP_N,N_Restart,xGP)
END IF

! Check whether (during restart) the statefile from which the restart is performed should be deleted
FlushInitialState = GETLOGICAL('FlushInitialState')

! Set wall time to the beginning of the simulation or when a restart is performed to the current wall time
RestartWallTime=PICLASTIME()
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
USE MOD_Restart_Vars           ,ONLY: DoRestart,RestartTime
USE MOD_ChangeBasis            ,ONLY: ChangeBasis3D
USE MOD_HDF5_Output            ,ONLY: FlushHDF5
#ifdef PP_POIS
USE MOD_Equation_Vars          ,ONLY: Phi
#endif /*PP_POIS*/
#if defined(PARTICLES)
USE MOD_Particle_Restart       ,ONLY: ParticleRestart
USE MOD_RayTracing             ,ONLY: RayTracing
#endif /*defined(PARTICLES)*/
#if USE_HDG
USE MOD_Restart_Tools          ,ONLY: RecomputeLambda
USE MOD_HDG                    ,ONLY: RestartHDG
#endif /*USE_HDG*/
USE MOD_Restart_Field          ,ONLY: FieldRestart
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                               :: StartT,EndT
#if USE_HDG
INTEGER(KIND=8)                    :: iter
#endif /*USE_HDG*/
!===================================================================================================================================
IF(DoRestart)THEN
  GETTIME(StartT)

  ! Restart field arrays
  CALL FieldRestart()

#ifdef PARTICLES
  ! Restart particle arrays
  CALL ParticleRestart()
  ! Get ray tracing volume and surface data
  CALL RayTracing()
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
  CALL FlushHDF5(RestartTime)

  GETTIME(EndT)
  CALL DisplayMessageAndTime(EndT-StartT, 'DONE!', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
ELSE ! no restart
  ! Delete all files since we are doing a fresh start
  CALL FlushHDF5()
END IF !IF(DoRestart)
END SUBROUTINE Restart


SUBROUTINE FinalizeRestart()
!===================================================================================================================================
! Finalizes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Restart_Vars     ,ONLY: Vdm_GaussNRestart_GaussN,RestartInitIsDone,DoMacroscopicRestart
#if defined(PARTICLES)
USE MOD_Particle_Restart ,ONLY: FinalizeParticleRestart
#endif /*defined(PARTICLES)*/
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

#if defined(PARTICLES)
CALL FinalizeParticleRestart()
#endif /*defined(PARTICLES)*/

END SUBROUTINE FinalizeRestart

END MODULE MOD_Restart
