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

MODULE MOD_Piclas_Init
!===================================================================================================================================
! contains global init and finalize
!===================================================================================================================================

INTERFACE InitPiclas
  MODULE PROCEDURE InitPiclas
END INTERFACE

INTERFACE FinalizePiclas
   MODULE PROCEDURE FinalizePiclas
END INTERFACE

PUBLIC:: FinalizePiclas
PUBLIC:: InitPiclas
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters.
!==================================================================================================================================
SUBROUTINE DefineParametersPiclas()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Piclas Initialization")

CALL prms%CreateIntOption(      'TimeStampLength' , 'Length of the floating number time stamp' , '21')
#ifdef PARTICLES
CALL prms%CreateLogicalOption(  'UseDSMC'         , "Flag for using DSMC methods"              , '.FALSE.')
CALL prms%CreateLogicalOption(  'UseRayTracing'   , "Flag for using ray tracing tools"         , '.FALSE.')
#endif

END SUBROUTINE DefineParametersPiclas


SUBROUTINE InitPiclas(IsLoadBalance)
!----------------------------------------------------------------------------------------------------------------------------------!
! init Piclas data structure
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Globals_Vars         ,ONLY: TimeStampLenStr,TimeStampLenStr2
USE MOD_Preproc
USE MOD_ReadInTools          ,ONLY: prms
USE MOD_Interpolation_Vars   ,ONLY: InterpolationInitIsDone
USE MOD_Restart_Vars         ,ONLY: RestartInitIsDone
USE MOD_Restart              ,ONLY: InitRestart
USE MOD_Restart_Vars         ,ONLY: DoRestart
USE MOD_Mesh                 ,ONLY: InitMesh
USE MOD_Mesh_Vars            ,ONLY: GetMeshMinMaxBoundariesIsDone
USE MOD_Equation             ,ONLY: InitEquation
USE MOD_GetBoundaryFlux      ,ONLY: InitBC
USE MOD_DG                   ,ONLY: InitDG
USE MOD_Mortar               ,ONLY: InitMortar
#if ! (USE_HDG)
USE MOD_PML                  ,ONLY: InitPML
#endif /*USE_HDG*/
USE MOD_Dielectric           ,ONLY: InitDielectric
USE MOD_Analyze              ,ONLY: InitAnalyze
USE MOD_RecordPoints         ,ONLY: InitRecordPoints
#if defined(ROS) || defined(IMPA)
USE MOD_LinearSolver         ,ONLY: InitLinearSolver
#endif /*ROS or IMPA*/
USE MOD_Restart_Vars         ,ONLY: N_Restart,InterpolateSolution,RestartNullifySolution
#if USE_MPI
USE MOD_MPI                  ,ONLY: InitMPIvars
#endif /*USE_MPI*/
#ifdef PARTICLES
USE MOD_DSMC_Vars            ,ONLY: UseDSMC
USE MOD_ParticleInit         ,ONLY: InitParticleGlobals,InitParticles
USE MOD_TTMInit              ,ONLY: InitTTM,InitIMD_TTM_Coupling
USE MOD_TTM_Vars             ,ONLY: DoImportTTMFile
USE MOD_Particle_Analyze     ,ONLY: InitParticleAnalyze
USE MOD_SurfaceModel_Analyze ,ONLY: InitSurfModelAnalyze
USE MOD_Particle_MPI         ,ONLY: InitParticleMPI
USE MOD_DSMC_Symmetry        ,ONLY: Init_Symmetry
#if USE_MPI
USE mod_readIMD              ,ONLY: initReadIMDdata,read_IMD_results
#endif /* USE_MPI */
#if defined(IMPA) || defined(ROS)
USE MOD_ParticleSolver       ,ONLY: InitPartSolver
#endif
#endif
#if USE_HDG
USE MOD_HDG                  ,ONLY: InitHDG
#endif
#if (PP_TimeDiscMethod==600)
USE MOD_RadiationTrans_Init        ,ONLY: InitRadiationTransport
USE MOD_Radiation_Init             ,ONLY: InitRadiation
#endif
USE MOD_Interfaces           ,ONLY: InitInterfaces
USE MOD_ReadInTools          ,ONLY: GETLOGICAL,GETREALARRAY,GETINT
USE MOD_TimeDisc_Vars        ,ONLY: TEnd
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
LOGICAL,INTENT(IN)      :: IsLoadBalance
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: TimeStampLength
!===================================================================================================================================
! Get length of the floating number time stamp
TimeStampLength = GETINT('TimeStampLength')
IF((TimeStampLength.LT.4).OR.(TimeStampLength.GT.30)) CALL abort(__STAMP__&
    ,'TimeStampLength cannot be smaller than 4 and not larger than 30')
WRITE(UNIT=TimeStampLenStr2,FMT='(I0)') TimeStampLength-4
! Check if TEnd overflows the output floating format
IF(TEnd.GE.1000.)THEN
  ! Add at least 1 digit
  TimeStampLength = TimeStampLength + MAX(1,CEILING(LOG10(TEnd))-3)
END IF
WRITE(UNIT=TimeStampLenStr ,FMT='(I0)') TimeStampLength

#ifdef PARTICLES
! DSMC handling:
useDSMC=GETLOGICAL('UseDSMC')

CALL Init_Symmetry()

#endif /*PARTICLES*/

! Initialization
IF(IsLoadBalance)THEN
  DoRestart=.TRUE.
  RestartInitIsDone=.TRUE.
  InterpolationInitIsDone=.TRUE.
  RestartNullifySolution=.FALSE.
  InterpolateSolution=.FALSE.
  N_Restart=PP_N
  CALL InitMortar()
ELSE
  CALL InitMortar()
  CALL InitRestart()
  GetMeshMinMaxBoundariesIsDone = .FALSE. ! Initialize this logical only once (assume that the mesh sizes does not change during load balance restarts)
END IF

#ifdef PARTICLES
CALL InitParticleGlobals(IsLoadBalance)
#endif

CALL InitMesh(2)
#if USE_MPI
CALL InitMPIvars()
#endif /*USE_MPI*/
CALL InitEquation()
CALL InitBC()
#if !(USE_HDG)
CALL InitPML() ! Perfectly Matched Layer (PML): electromagnetic-wave-absorbing layer
#endif /*USE_HDG*/
CALL InitDielectric() ! Dielectric media
CALL InitDG()
#if defined(ROS) || defined(IMPA)
CALL InitLinearSolver()
#endif /*ROS /IMEX*/
#ifdef PARTICLES
CALL InitParticleMPI
CALL InitParticles()
#if defined(IMPA) || defined(ROS)
CALL InitPartSolver()
#endif
#endif
CALL InitAnalyze()
CALL InitRecordPoints()
#ifdef PARTICLES
CALL InitParticleAnalyze()
CALL InitSurfModelAnalyze()
#endif

#if USE_HDG
CALL InitHDG() ! Hybridizable Discontinuous Galerkin Method (HDGSEM)
#endif

#if (PP_TimeDiscMethod==600)
CALL InitRadiation()
CALL InitRadiationTransport()
#endif

#ifdef PARTICLES
! Old IMD format
  CALL InitTTM() ! FD grid based data from a Two-Temperature Model (TTM) from Molecular Dynamics (MD) Code IMD
IF(DoImportTTMFile)THEN
  CALL InitIMD_TTM_Coupling() ! use MD and TTM data to distribute the cell averaged charge to the atoms/ions
END IF
#if USE_MPI
! New IMD binary format (not TTM needed as this information is stored on the atoms)
CALL initReadIMDdata()
CALL read_IMD_results()
#endif /* USE_MPI */
#endif /*PARTICLES*/

CALL InitInterfaces() ! set Riemann solver identifier for face connectivity (vacuum, dielectric, PML ...)

! !do this last
! write out parameters that are not used and remove multiple and unused, that are not needed to do restart if no parameter.ini is
! read in
IF (.NOT.IsLoadBalance) THEN
  CALL prms%WriteUnused()
  ! CALL prms%RemoveUnnecessary() ! We need to keep the numberedMulti in case we want the default value
END IF


END SUBROUTINE InitPiclas


SUBROUTINE FinalizePiclas(IsLoadBalance)
!----------------------------------------------------------------------------------------------------------------------------------!
! finalize Piclas data structure
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Commandline_Arguments      ,ONLY: FinalizeCommandlineArguments
USE MOD_Globals
USE MOD_ReadInTools                ,ONLY: prms,FinalizeParameters
USE MOD_Restart                    ,ONLY: FinalizeRestart
USE MOD_Interpolation              ,ONLY: FinalizeInterpolation
USE MOD_Mesh                       ,ONLY: FinalizeMesh
USE MOD_Equation                   ,ONLY: FinalizeEquation
USE MOD_Interfaces                 ,ONLY: FinalizeInterfaces
USE MOD_GetBoundaryFlux            ,ONLY: FinalizeBC
USE MOD_DG                         ,ONLY: FinalizeDG
USE MOD_Mortar                     ,ONLY: FinalizeMortar
USE MOD_Dielectric                 ,ONLY: FinalizeDielectric
#if ! (USE_HDG)
USE MOD_PML                        ,ONLY: FinalizePML
#else
USE MOD_HDG                        ,ONLY: FinalizeHDG
#endif /*USE_HDG*/
USE MOD_Analyze                    ,ONLY: FinalizeAnalyze
USE MOD_RecordPoints               ,ONLY: FinalizeRecordPoints
USE MOD_RecordPoints_Vars          ,ONLY: RP_Data
#if defined(ROS) || defined(IMPA)
USE MOD_LinearSolver               ,ONLY: FinalizeLinearSolver
#endif /*IMEX*/
#if USE_MPI
USE MOD_MPI                        ,ONLY: FinalizeMPI
USE MOD_MPI_Shared                 ,ONLY: FinalizeMPIShared
#endif /*USE_MPI*/
#ifdef PARTICLES
USE MOD_RayTracing_Init            ,ONLY: FinalizeRayTracing
USE MOD_Particle_Surfaces          ,ONLY: FinalizeParticleSurfaces
USE MOD_InitializeBackgroundField  ,ONLY: FinalizeBackGroundField
USE MOD_SuperB_Init                ,ONLY: FinalizeSuperB
USE MOD_Particle_Mesh              ,ONLY: FinalizeParticleMesh
USE MOD_Particle_Analyze           ,ONLY: FinalizeParticleAnalyze
USE MOD_PICDepo                    ,ONLY: FinalizeDeposition
USE MOD_PICInterpolation           ,ONLY: FinalizePICInterpolation
USE MOD_ParticleInit               ,ONLY: FinalizeParticles
USE MOD_Particle_Sampling_Adapt    ,ONLY: FinalizeParticleSamplingAdaptive
USE MOD_Particle_Boundary_Init     ,ONLY: FinalizeParticleBoundary
USE MOD_TTMInit                    ,ONLY: FinalizeTTM
USE MOD_DSMC_Init                  ,ONLY: FinalizeDSMC
USE MOD_MCC_Init                   ,ONLY: FinalizeMCC
USE MOD_SurfaceModel_Porous        ,ONLY: FinalizePorousBoundaryCondition
#if (PP_TimeDiscMethod==300)
USE MOD_FPFlow_Init                ,ONLY: FinalizeFPFlow
#endif
#if (PP_TimeDiscMethod==400)
USE MOD_BGK_Init                   ,ONLY: FinalizeBGK
#endif
USE MOD_SurfaceModel_Init          ,ONLY: FinalizeSurfaceModel
USE MOD_SurfaceModel_Analyze       ,ONLY: FinalizeSurfaceModelAnalyze
USE MOD_Particle_Boundary_Sampling ,ONLY: FinalizeParticleBoundarySampling
USE MOD_Particle_Vars              ,ONLY: ParticlesInitIsDone
USE MOD_PIC_Vars                   ,ONLY: PICInitIsDone
#if USE_MPI
USE MOD_Particle_MPI               ,ONLY: FinalizeParticleMPI
USE MOD_Particle_MPI_Vars          ,ONLY: ParticleMPIInitisdone
#if defined(MEASURE_MPI_WAIT)
USE MOD_MPI                        ,ONLY: OutputMPIW8Time
#endif /*defined(MEASURE_MPI_WAIT)*/
#endif /*USE_MPI*/
#endif /*PARTICLES*/
USE MOD_IO_HDF5                    ,ONLY: FinalizeElemData,ElementOut
USE MOD_TimeDiscInit               ,ONLY: FinalizeTimeDisc
#if (PP_TimeDiscMethod==600)
USE MOD_Radiation_Init             ,ONLY: FinalizeRadiation
USE MOD_RadiationTrans_Init        ,ONLY: FinalizeRadiationTransport
USE MOD_Photon_Tracking            ,ONLY: FinalizePhotonSurfSample
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
LOGICAL,INTENT(IN)      :: IsLoadBalance
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Time
!===================================================================================================================================
CALL FinalizeElemData(ElementOut)
!Finalize
CALL FinalizeRecordPoints()
CALL FinalizeAnalyze()
CALL FinalizeDG()
#if defined(IMPA) || defined(ROS)
!CALL FinalizeCSR()
CALL FinalizeLinearSolver()
#endif /*IMEX*/
CALL FinalizeDielectric()
#if !(USE_HDG)
CALL FinalizePML()
#else
CALL FinalizeHDG()
#endif /*USE_HDG*/
CALL FinalizeEquation()
CALL FinalizeBC()
IF(.NOT.IsLoadBalance) CALL FinalizeInterpolation()
CALL FinalizeRestart()
CALL FinalizeMesh()
CALL FinalizeMortar()
#ifdef PARTICLES
CALL FinalizeRayTracing()
CALL FinalizeSurfaceModel()
CALL FinalizeSurfaceModelAnalyze()
CALL FinalizeParticleBoundarySampling()
CALL FinalizePorousBoundaryCondition()
CALL FinalizeParticleSurfaces()
CALL FinalizeParticleMesh()
CALL FinalizeParticleAnalyze()
CALL FinalizeDeposition()
CALL FinalizePICInterpolation()
#if USE_MPI
CALL FinalizeParticleMPI()
#endif /*USE_MPI*/
CALL FinalizeDSMC()
CALL FinalizeMCC()
#if (PP_TimeDiscMethod==300)
CALL FinalizeFPFlow()
#endif
#if (PP_TimeDiscMethod==400)
CALL FinalizeBGK()
#endif
CALL FinalizeParticles()
CALL FinalizeParticleSamplingAdaptive(IsLoadBalance)
CALL FinalizeParticleBoundary()
CALL FinalizeBackGroundField()
CALL FinalizeSuperB()
#endif /*PARTICLES*/
#if USE_MPI
CALL FinalizeMPI()
#endif /*USE_MPI*/

#ifdef PARTICLES
ParticlesInitIsDone = .FALSE.
PICInitIsDone = .FALSE.
#if USE_MPI
ParticleMPIInitIsDone=.FALSE.
#endif /*USE_MPI*/

#if (PP_TimeDiscMethod==600)
CALL FinalizeRadiation()
CALL FinalizeRadiationTransport()
CALL FinalizePhotonSurfSample()
#endif

CALL FinalizeTTM() ! FD grid based data from a Two-Temperature Model (TTM) from Molecular Dynamics (MD) Code IMD
#endif /*PARTICLES*/

CALL FinalizeInterfaces()
CALL prms%finalize(IsLoadBalance) ! is the same as CALL FinalizeParameters(), but considers load balancing
CALL FinalizeCommandlineArguments()

CALL FinalizeTimeDisc()
! mssing arrays to deallocate
SDEALLOCATE(RP_Data)

! Before program termination: Finalize load balance
! Measure simulation duration
Time=PICLASTIME()
CALL FinalizeLoadBalance(IsLoadBalance)
IF(.NOT.IsLoadBalance)THEN
  CALL DisplaySimulationTime(Time, StartTime, 'FINISHED')
#if USE_MPI
  ! Free the communicators!
  CALL FinalizeMPIShared()
#if defined(MEASURE_MPI_WAIT)
  ! Collect the MPI_WAIT() over all procs and output
  IF(nProcessors.GT.1) CALL OutputMPIW8Time()
#endif /*defined(MEASURE_MPI_WAIT)*/
  ! Free the last communicator after OutputMPIW8Time
  CALL MPI_BARRIER  (MPI_COMM_PICLAS,iError)
  IF(MPI_COMM_PICLAS.NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(MPI_COMM_PICLAS,IERROR)
#endif /*USE_MPI*/
ELSE
  CALL DisplaySimulationTime(Time, StartTime, 'RUNNING')
END IF ! .NOT.IsLoadBalance
SWRITE(UNIT_stdOut,'(132("="))')

END SUBROUTINE FinalizePiclas


SUBROUTINE FinalizeLoadBalance(IsLoadBalance)
!===================================================================================================================================
! Deallocate arrays
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_LoadBalance_Vars
#if USE_LOADBALANCE
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars   ,ONLY: MPI_COMM_SHARED
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)  :: IsLoadBalance
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(nPartsPerElem)
SDEALLOCATE(nDeposPerElem)
SDEALLOCATE(nTracksPerElem)
SDEALLOCATE(nSurfacefluxPerElem)
SDEALLOCATE(nPartsPerBCElem)
SDEALLOCATE(nSurfacePartsPerElem)

SDEALLOCATE(LoadDistri)
SDEALLOCATE(PartDistri)
SDEALLOCATE(ElemGlobalTime)
SDEALLOCATE(ElemHDGSides)
SDEALLOCATE(ElemTime_tmp)
!SDEALLOCATE(ElemTime)

IF(.NOT.IsLoadBalance) THEN
#if USE_LOADBALANCE
  IF (ASSOCIATED(ElemInfoRank_Shared)) THEN
    ! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
    CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
    CALL UNLOCK_AND_FREE(ElemInfoRank_Shared_Win)
    CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

    ! Then, free the pointers or arrays
    NULLIFY(ElemInfoRank_Shared)
  END IF

  SDEALLOCATE(tCurrent)
  InitLoadBalanceIsDone = .FALSE.

  SDEALLOCATE(offsetElemMPIOld)
#endif /*USE_LOADBALANCE*/
END IF

END SUBROUTINE FinalizeLoadBalance


END MODULE MOD_Piclas_Init
