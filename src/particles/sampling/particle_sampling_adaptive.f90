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

MODULE MOD_Particle_Sampling_Adapt
!===================================================================================================================================
!> Subroutines required for the sampling of macroscopic properties at elements with a boundary for the porous and adaptive boundary
!> conditions
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: DefineParametersParticleSamplingAdaptive, InitAdaptiveBCSampling, AdaptiveBCSampling, FinalizeParticleSamplingAdaptive
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for the sampling for the adaptive boundary conditions
!==================================================================================================================================
SUBROUTINE DefineParametersParticleSamplingAdaptive()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%CreateLogicalOption('AdaptiveBC-AverageValuesOverBC',  'Flag to enable the usage of average macroscopic values'//&
                                                                    ' across the whole boundary.', '.FALSE.')
CALL prms%CreateRealOption(   'AdaptiveBC-RelaxationFactor',  'Relaxation factor for weighting of current'//&
                                                              ' values with those of previous iterations.', '0.001')
CALL prms%CreateIntOption(    'AdaptiveBC-SamplingIteration', 'Number of iterations the values will be sampled before updating'//&
                                                              ' the utilized value. Alternative is to constantly update'//&
                                                              ' with a relaxation factor', '0')
CALL prms%CreateLogicalOption('AdaptiveBC-TruncateRunningAverage',  'Flag to enable a truncated running average for'//&
                                                                    ' the last number of SamplingIteration')
END SUBROUTINE DefineParametersParticleSamplingAdaptive


SUBROUTINE InitAdaptiveBCSampling()
!===================================================================================================================================
!> Routine to initialize the sampling of elements adjacent to adaptive surface flux and/or porous BCs
!> 1) Count the number of sample elements and create mapping from ElemID to SampleElemID
!>  1a) Add elements for the adaptive surface flux
!>  1b) Add elements for the porous BCs
!> 2) Allocate the sampling arrays and create mapping from SampleElemID to ElemID
!>  2a) Initializing the array with the given velocity vector and magnitude
!> 3) Read-in of the additional variables for sampling and and array allocation
!> 4) Sampling of near adaptive boundary element values in the first time step to get initial distribution
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Particle_Sampling_Vars
USE MOD_HDF5_INPUT              ,ONLY: ReadArray, ReadAttribute, DatasetExists, GetDataSize
USE MOD_Mesh_Vars               ,ONLY: offsetElem, nElems, SideToElem
USE MOD_Particle_Vars           ,ONLY: Species, nSpecies, UseCircularInflow
USE MOD_Particle_Surfaces_Vars  ,ONLY: BCdata_auxSF, SurfFluxSideSize, SurfMeshSubSideData
USE MOD_Restart_Vars            ,ONLY: DoRestart
USE MOD_SurfaceModel_Vars       ,ONLY: nPorousBC
USE MOD_Particle_Boundary_Vars  ,ONLY: nPorousSides, PorousBCInfo_Shared, SurfSide2GlobalSide
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared, ElemVolume_Shared
USE MOD_LoadBalance_Vars        ,ONLY: DoLoadBalance, PerformLoadBalance
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars        ,ONLY: UseH5IOLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                           :: UseAdaptiveType4
INTEGER                           :: iElem, iSpec, iSF, iSide, ElemID, SampleElemID, GlobalSideID, GlobalElemID, currentBC
INTEGER                           :: jSample, iSample, BCSideID
INTEGER                           :: nSurfacefluxBCs
#if USE_MPI
INTEGER                           :: offSetElemAdaptBCSampleMPI(0:nProcessors-1)
#endif
!===================================================================================================================================

LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT SAMPLING FOR ADAPTIVE BC ...'

#if USE_LOADBALANCE
IF(DoLoadBalance.AND.(.NOT.UseH5IOLoadBalance)) THEN
  CALL abort(__STAMP__,'ERROR: Adaptive BCs currently only support a load balance using an HDF5 output (UseH5IOLoadBalance = T)!')
END IF
#endif /*USE_LOADBALANCE*/

UseAdaptiveType4 = .FALSE.
AdaptBCSampleElemNum = 0
ALLOCATE(AdaptBCMapElemToSample(nElems))
AdaptBCMapElemToSample = 0
ALLOCATE(AdaptiveData(1:7*nSpecies,1:nElems))
AdaptiveData = 0.

nSurfacefluxBCs = MAXVAL(Species(:)%nSurfacefluxBCs)

IF(UseCircularInflow) THEN
  ALLOCATE(AdaptBCAreaSurfaceFlux(1:nSpecies,1:nSurfacefluxBCs))
  AdaptBCAreaSurfaceFlux = 0.
END IF

ALLOCATE(AdaptBCVolSurfaceFlux(1:nSpecies,1:nSurfacefluxBCs))
AdaptBCVolSurfaceFlux = 0.

! 1) Count the number of sample elements and create mapping from ElemID to SampleElemID
! 1a) Add elements for the adaptive surface flux
DO iSpec=1,nSpecies
  DO iSF=1,Species(iSpec)%nSurfacefluxBCs
    currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
    ! Set flag for every processor as read-in of the state file is performed in parallel to avoid deadlock
    IF(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType.EQ.4) UseAdaptiveType4 = .TRUE.
    ! Skip processors without a surface flux
    IF (BCdata_auxSF(currentBC)%SideNumber.EQ.0) CYCLE
    ! Skip a regular surface flux
    IF (.NOT.Species(iSpec)%Surfaceflux(iSF)%Adaptive) CYCLE
    ! Loop over sides on the surface flux
    DO iSide=1,BCdata_auxSF(currentBC)%SideNumber
      BCSideID = BCdata_auxSF(currentBC)%SideList(iSide)
      ElemID = SideToElem(S2E_ELEM_ID,BCSideID)
      ! Skip elements outside of the circular inflow
      IF(Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
        IF(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide).EQ.1) CYCLE
        ! Determine the area of the surface flux
        DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
          AdaptBCAreaSurfaceFlux(iSpec,iSF)=AdaptBCAreaSurfaceFlux(iSpec,iSF)+SurfMeshSubSideData(iSample,jSample,BCSideID)%area
        END DO; END DO
      END IF
      ! Calculate the sum of the cell volumes adjacent to the surface flux
      AdaptBCVolSurfaceFlux(iSpec,iSF) = AdaptBCVolSurfaceFlux(iSpec,iSF) + ElemVolume_Shared(GetCNElemID(ElemID+offsetElem))
      ! Only add elements once
      IF(AdaptBCMapElemToSample(ElemID).NE.0) CYCLE
      ! Add the element to the BC sampling
      AdaptBCSampleElemNum = AdaptBCSampleElemNum + 1
      AdaptBCMapElemToSample(ElemID) = AdaptBCSampleElemNum
    END DO
  END DO
END DO

#if USE_MPI
IF(UseCircularInflow) THEN
  IF(MPIRoot)THEN
    CALL MPI_REDUCE(MPI_IN_PLACE,AdaptBCAreaSurfaceFlux,nSpecies*nSurfacefluxBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS,iError)
  ELSE
    CALL MPI_REDUCE(AdaptBCAreaSurfaceFlux,0,nSpecies*nSurfacefluxBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS,iError)
  END IF
END IF
#endif /*USE_MPI*/

#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,AdaptBCVolSurfaceFlux,nSpecies*nSurfacefluxBCs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/

! 1b) Add elements for the porous BCs
IF(nPorousBC.GT.0) THEN
  DO iSide = 1, nPorousSides
    GlobalSideID = SurfSide2GlobalSide(SURF_SIDEID,PorousBCInfo_Shared(2,iSide))
    GlobalElemID = SideInfo_Shared(SIDE_ELEMID,GlobalSideID)
    ! Only treat your proc-local elements
    IF ((GlobalElemID.LT.1+offSetElem).OR.(GlobalElemID.GT.nElems+offSetElem)) CYCLE
    ElemID = GlobalElemID - offsetElem
    ! Only add elements once
    IF(AdaptBCMapElemToSample(ElemID).NE.0) CYCLE
    ! Add the element to the BC sampling
    AdaptBCSampleElemNum = AdaptBCSampleElemNum + 1
    AdaptBCMapElemToSample(ElemID) = AdaptBCSampleElemNum
  END DO
END IF

! 2) Allocate the sampling arrays and create mapping from SampleElemID to ElemID
SampleElemID = 0
ALLOCATE(AdaptBCMapSampleToElem(AdaptBCSampleElemNum))
DO iElem = 1,nElems
  IF(AdaptBCMapElemToSample(iElem).NE.0) THEN
    AdaptBCMapSampleToElem(AdaptBCMapElemToSample(iElem)) = iElem
    SampleElemID = SampleElemID + 1
  END IF
END DO

! Sanity check
IF(SampleElemID.NE.AdaptBCSampleElemNum) THEN
  IPWRITE(*,*) 'Number of elements counted in the mapping array: ', SampleElemID
  IPWRITE(*,*) 'Number of elements counted during surface flux and porous BC checks: ', AdaptBCSampleElemNum
  CALL abort(__STAMP__,&
      'ERROR: Mapping between sampling elements and regular elements is not complete!')
END IF

#if USE_MPI
! Gather the number of sampling elements per proc
CALL MPI_GATHER(AdaptBCSampleElemNum,1,MPI_INTEGER_INT_KIND,offSetElemAdaptBCSampleMPI,1,MPI_INTEGER_INT_KIND,0,MPI_COMM_PICLAS,iError)
! Distribute the number of elements per proc to each each proc
CALL MPI_BCAST(offSetElemAdaptBCSampleMPI,nProcessors,MPI_INTEGER,0,MPI_COMM_PICLAS,iERROR)
! Determine the offset for the sampling elements
IF(myRank.EQ.0) THEN
  offSetElemAdaptBCSample = 0
ELSE
  offSetElemAdaptBCSample = SUM(offSetElemAdaptBCSampleMPI(0:myRank-1))
END IF
! Determine the total number of sampling elements (global)
AdaptBCSampleElemNumGlobal = SUM(offSetElemAdaptBCSampleMPI)
#else
AdaptBCSampleElemNumGlobal = AdaptBCSampleElemNum
#endif /*USE_MPI*/

ALLOCATE(AdaptBCMacroVal(1:7,1:AdaptBCSampleElemNum,1:nSpecies))
AdaptBCMacroVal(:,:,:) = 0.0
ALLOCATE(AdaptBCSample(1:8,1:AdaptBCSampleElemNum,1:nSpecies))
AdaptBCSample = 0.0

! 2a) Initializing an array with the given velocity vector and magnitude. It is used as a fallback, when the sampled velocity
! in the cell is zero (e.g. when starting a simulation with zero particles)
ALLOCATE(AdaptBCBackupVelocity(1:3,1:AdaptBCSampleElemNum,1:nSpecies))
AdaptBCBackupVelocity = 0.0
DO iSpec=1,nSpecies
  DO iSF=1,Species(iSpec)%nSurfacefluxBCs
    currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
    ! Skip processors without a surface flux
    IF (BCdata_auxSF(currentBC)%SideNumber.EQ.0) CYCLE
    ! Loop over sides on the surface flux
    DO iSide=1,BCdata_auxSF(currentBC)%SideNumber
      ! Skip elements outside of the circular inflow
      IF(Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
        IF(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide).EQ.1) CYCLE
      END IF
      ElemID = SideToElem(S2E_ELEM_ID,BCdata_auxSF(currentBC)%SideList(iSide))
      SampleElemID = AdaptBCMapElemToSample(ElemID)
      IF(SampleElemID.GT.0) THEN
        AdaptBCBackupVelocity(1:3,SampleElemID,iSpec) = Species(iSpec)%Surfaceflux(iSF)%VeloIC*Species(iSpec)%Surfaceflux(iSF)%VeloVecIC(1:3)
      END IF
    END DO
  END DO
END DO

! 3) Read-in of the additional variables for sampling and array allocation
AdaptBCAverageValBC = GETLOGICAL('AdaptiveBC-AverageValuesOverBC')
AdaptBCRelaxFactor = GETREAL('AdaptiveBC-RelaxationFactor')
AdaptBCSampIter = GETINT('AdaptiveBC-SamplingIteration')

! Type=4 requires the number of particles that left through that BC in the previous time step
! Only allocate and nullify, if new simulation or restart, keep values during load balance step to avoid wrong mass flow during
IF(.NOT.PerformLoadBalance) THEN
  IF(UseAdaptiveType4) THEN
    ALLOCATE(AdaptBCPartNumOut(1:nSpecies,1:nSurfacefluxBCs))
    AdaptBCPartNumOut = 0
  END IF
  IF(AdaptBCAverageValBC) THEN
    ALLOCATE(AdaptBCAverageMacroVal(1:3,1:nSpecies,1:nSurfacefluxBCs))
    AdaptBCAverageMacroVal = 0
  END IF
END IF

IF(AdaptBCSampIter.GT.0) THEN
  AdaptBCTruncAverage = GETLOGICAL('AdaptiveBC-TruncateRunningAverage','.TRUE.')
  IF(AdaptBCTruncAverage) THEN
    ALLOCATE(AdaptBCAverage(1:8,1:AdaptBCSampIter,1:AdaptBCSampleElemNum,1:nSpecies))
    AdaptBCAverage = 0.0
  END IF
ELSE
  AdaptBCTruncAverage = GETLOGICAL('AdaptiveBC-TruncateRunningAverage','.FALSE.')
  IF(AdaptBCTruncAverage) THEN
    AdaptBCTruncAverage = .FALSE.
    LBWRITE(*,*) '| WARNING: TruncateRunningAverage is set to true, but SamplingIteration is zero. RelaxationFactor is utilized instead.'
  END IF
END IF

IF(.NOT.DoRestart.AND..NOT.PerformLoadBalance) THEN
! 4) Initialize values from parameter file as a fallback
  DO iSpec=1,nSpecies
    DO iSF=1,Species(iSpec)%nSurfacefluxBCs
      currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
      ! Skip processors without a surface flux
      IF (BCdata_auxSF(currentBC)%SideNumber.EQ.0) CYCLE
      ! Loop over sides on the surface flux
      DO iSide=1,BCdata_auxSF(currentBC)%SideNumber
        ! Skip elements outside of the circular inflow
        IF(Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
          IF(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide).EQ.1) CYCLE
        END IF
        ElemID = SideToElem(S2E_ELEM_ID,BCdata_auxSF(currentBC)%SideList(iSide))
        SampleElemID = AdaptBCMapElemToSample(ElemID)
        IF(SampleElemID.GT.0) THEN
          AdaptBCMacroVal(1:3,SampleElemID,iSpec) = Species(iSpec)%Surfaceflux(iSF)%VeloIC*Species(iSpec)%Surfaceflux(iSF)%VeloVecIC(1:3)
          AdaptBCMacroVal(4,SampleElemID,iSpec) = Species(iSpec)%Surfaceflux(iSF)%PartDensity
        END IF
      END DO
    END DO
  END DO
  LBWRITE(*,*) '| Macroscopic values have been initialized with the input parameters for velocity and number density (fallback).'
END IF

END SUBROUTINE InitAdaptiveBCSampling


SUBROUTINE AdaptiveBCSampling(initSampling_opt,initTruncAverage_opt)
!===================================================================================================================================
!> Sampling of variables (part-density and velocity) for adaptive BC and porous BC elements
!> 0) Determination whether a calculation of macroscopic values is required and set the counter for AdaptBCTruncAverage
!> 1) Loop over all particles, while cycling elements which are not at a boundary through AdaptBCMapElemToSample
!> 2) AdaptBCAverageValBC: Calculate the macroscopic values averaged over the whole BC
!>  2a) SamplingIteration: use the truncated running average (e.g. sample over the last 100 iterations)
!>      OR
!>  2a) RelaxationFactor: update the current values by a certain percentage
!> 3) Populate the element local arrays for AdaptBCAverageValBC
!>  3a) AdaptBCAverageValBC: Populate the element local arrays
!> 4) Calculation of the cell-local values (if not AdaptBCAverageValBC)
!>  4a) SamplingIteration
!>      OR
!>  4b) RelaxationFactor
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_Particle_Sampling_Vars
USE MOD_DSMC_Vars              ,ONLY: RadialWeighting
USE MOD_Mesh_Vars              ,ONLY: offsetElem, SideToElem
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID
USE MOD_Part_Tools             ,ONLY: GetParticleWeight
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemVolume_Shared
USE MOD_Particle_Vars          ,ONLY: PartState, PDM, PartSpecies, Species, nSpecies, PEM, usevMPF
USE MOD_Timedisc_Vars          ,ONLY: iter
USE MOD_Particle_Surfaces_Vars ,ONLY: BCdata_auxSF
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers     ,ONLY: LBStartTime, LBElemSplitTime, LBPauseTime
USE MOD_LoadBalance_vars       ,ONLY: nPartsPerBCElem
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL, INTENT(IN), OPTIONAL   :: initSampling_opt
LOGICAL, INTENT(IN), OPTIONAL   :: initTruncAverage_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: ElemID, CNElemID, SampleElemID, iPart, iSpec, SamplingIteration, TruncIter, ModIter
INTEGER                         :: RestartSampIter, nSurfacefluxBCs
INTEGER                         :: BCSideID, currentBC, iSF, iSide
REAL                            :: partWeight, TTrans_TempFac, RelaxationFactor
REAL,ALLOCATABLE                :: AdaptBCMeanValues(:,:,:)
LOGICAL                         :: initSampling, CalcValues, initTruncAverage
#if USE_LOADBALANCE
REAL                            :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================

! Optional flag for the utilization of the routine for an initial sampling of the density and pressure distribution before simstart
IF(PRESENT(initSampling_opt)) THEN
  initSampling = initSampling_opt
ELSE
  initSampling = .FALSE.
END IF

! Optional flag for the utilization of the routine for an initial calculation of the density and pressure distribution from the
! read-in truncated average sample values before simulation start
IF(PRESENT(initTruncAverage_opt)) THEN
  initTruncAverage = initTruncAverage_opt
ELSE
  initTruncAverage = .FALSE.
END IF

CalcValues = .FALSE.
RestartSampIter = 0
nSurfacefluxBCs = MAXVAL(Species(:)%nSurfacefluxBCs)

! If no particles are present during the initial sampling, leave the routine, otherwise initial variables for the
! adaptive inlet surface flux will be overwritten by zero's.
IF (PDM%ParticleVecLength.LT.1.AND..NOT.initTruncAverage) RETURN

! 0) Determination whether a calculation of macroscopic values is required and set the counter for AdaptBCTruncAverage
! Calculate the counter for the truncated moving average
IF(AdaptBCTruncAverage.AND..NOT.initSampling) THEN
  IF (AdaptBCSampleElemNum.GT.0) THEN
    RestartSampIter =  INT(iter,4)+AdaptBCSampIterReadIn
    ModIter = MOD(RestartSampIter,AdaptBCSampIter)
    TruncIter = MERGE(AdaptBCSampIter,ModIter,ModIter.EQ.0)
    ! Delete the oldest sample (required, otherwise it would be added to the new sample)
    IF(.NOT.initTruncAverage) AdaptBCAverage(1:8,TruncIter,1:AdaptBCSampleElemNum,1:nSpecies) = 0.
  END IF
END IF

! Determine whether a calculation of macroscopic values is required and set the appropriate counters
IF(initSampling) THEN
  RelaxationFactor = 1
  SamplingIteration = 1
  CalcValues = .TRUE.
ELSE
  RelaxationFactor = AdaptBCRelaxFactor
  IF(AdaptBCSampIter.GT.0) THEN
    IF(AdaptBCTruncAverage.AND.(RestartSampIter.GT.0).AND.(RestartSampIter.LT.AdaptBCSampIter)) THEN
      ! Truncated average: get the correct number of samples to calculate the average number density while the
      ! sampling array is populated
      SamplingIteration = RestartSampIter
    ELSE
      SamplingIteration = AdaptBCSampIter
    END IF
    ! Determine whether the macroscopic values shall be calculated from the sample (e.g. every 100 steps)
    CalcValues = MOD(INT(iter,4),SamplingIteration).EQ.0
  END IF
END IF

! 1) Loop over all particles, while cycling elements which are not at a boundary through AdaptBCMapElemToSample
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
IF(.NOT.initTruncAverage) THEN
  DO iPart = 1, PDM%ParticleVecLength
    IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
    ElemID = PEM%LocalElemID(iPart)
    SampleElemID = AdaptBCMapElemToSample(ElemID)
    ! Cycle particles inside non-sampling elements
    IF(SampleElemID.LE.0) CYCLE
#if USE_LOADBALANCE
    nPartsPerBCElem(ElemID) = nPartsPerBCElem(ElemID) + 1
#endif /*USE_LOADBALANCE*/
    ! Sample the particle properties
    iSpec = PartSpecies(iPart)
    partWeight = GetParticleWeight(iPart)
    IF(AdaptBCTruncAverage.AND..NOT.initSampling) THEN
      ! Store the samples of the last AdaptBCSampIter and replace the oldest with the newest sample
      AdaptBCAverage(1:3,TruncIter,SampleElemID, iSpec) = AdaptBCAverage(1:3,TruncIter,SampleElemID,iSpec) + PartState(4:6,iPart) * partWeight
      AdaptBCAverage(4:6,TruncIter,SampleElemID, iSpec) = AdaptBCAverage(4:6,TruncIter,SampleElemID,iSpec) + PartState(4:6,iPart)**2 * partWeight
      AdaptBCAverage(7,  TruncIter,SampleElemID, iSpec) = AdaptBCAverage(7,  TruncIter,SampleElemID,iSpec) + 1.0  ! simulation particle number
      AdaptBCAverage(8,  TruncIter,SampleElemID, iSpec) = AdaptBCAverage(8,  TruncIter,SampleElemID,iSpec) + partWeight
    ELSE
      AdaptBCSample(1:3,SampleElemID, iSpec) = AdaptBCSample(1:3,SampleElemID,iSpec) + PartState(4:6,iPart) * partWeight
      AdaptBCSample(4:6,SampleElemID, iSpec) = AdaptBCSample(4:6,SampleElemID,iSpec) + PartState(4:6,iPart)**2 * partWeight
      AdaptBCSample(7,SampleElemID, iSpec) = AdaptBCSample(7,SampleElemID, iSpec) + 1.0  ! simulation particle number
      AdaptBCSample(8,SampleElemID, iSpec) = AdaptBCSample(8,SampleElemID, iSpec) + partWeight
    END IF
  END DO
END IF
#if USE_LOADBALANCE
CALL LBPauseTime(LB_ADAPTIVE,tLBStart)
#endif /*USE_LOADBALANCE*/

IF(AdaptBCTruncAverage.AND..NOT.initSampling.AND.AdaptBCSampleElemNum.GT.0) THEN
  ! Sum-up the complete sample over the number of sampling iterations
  AdaptBCSample(1:8,1:AdaptBCSampleElemNum,1:nSpecies) = SUM(AdaptBCAverage(1:8,:,1:AdaptBCSampleElemNum,1:nSpecies),2)
END IF

! 2) AdaptBCAverageValBC: Calculate the macroscopic values averaged over the whole BC
IF(AdaptBCAverageValBC) THEN
  IF(CalcValues.OR.AdaptBCTruncAverage.OR.(AdaptBCSampIter.EQ.0)) THEN
    ALLOCATE(AdaptBCMeanValues(1:8,1:nSpecies,1:nSurfacefluxBCs))
    AdaptBCMeanValues = 0.
    ! Sum up values per species and surface flux
    DO iSpec=1,nSpecies
      DO iSF=1,Species(iSpec)%nSurfacefluxBCs
        currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
        ! Loop over sides on the surface flux
        DO iSide=1,BCdata_auxSF(currentBC)%SideNumber
          BCSideID = BCdata_auxSF(currentBC)%SideList(iSide)
          ElemID = SideToElem(S2E_ELEM_ID,BCSideID)
          SampleElemID = AdaptBCMapElemToSample(ElemID)
          IF (SampleElemID.GT.0) THEN
            ! Determine the mean flow velocity
            AdaptBCMeanValues(1:8,iSpec,iSF) = AdaptBCMeanValues(1:8,iSpec,iSF) + AdaptBCSample(1:8,SampleElemID,iSpec)
            ! Resetting sampled values
            AdaptBCSample(1:8,SampleElemID,iSpec) = 0.
          END IF
        END DO
      END DO
    END DO
    ! MPI Communication
#if USE_MPI
    IF(MPIRoot)THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,AdaptBCMeanValues,8*nSpecies*nSurfacefluxBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS,iError)
    ELSE
      CALL MPI_REDUCE(AdaptBCMeanValues,0.,8*nSpecies*nSurfacefluxBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS,iError)
    END IF
#endif /*USE_MPI*/
    IF(MPIRoot) THEN
      DO iSpec=1,nSpecies
        DO iSF=1,Species(iSpec)%nSurfacefluxBCs
          IF(AdaptBCMeanValues(7,iSpec,iSF).GT.0.) THEN
            ! Average the values across the BC using the weight sum
            AdaptBCMeanValues(1:6,iSpec,iSF) = AdaptBCMeanValues(1:6,iSpec,iSF) / AdaptBCMeanValues(8,iSpec,iSF)
            ! ================================================================
            ! Sampling iteration: sampling for AdaptBCSampIter iterations, calculating the macro values and resetting sample OR
            ! AdaptBCTruncAverage: continuous average of the last AdaptBCSampIter iterations
            IF(AdaptBCSampIter.GT.0) THEN
              ! Calculate the average number density
              IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
                AdaptBCAverageMacroVal(1,iSpec,iSF) = AdaptBCMeanValues(8,iSpec,iSF) / REAL(SamplingIteration) &
                                                      / AdaptBCVolSurfaceFlux(iSpec,iSF)
              ELSE
                AdaptBCAverageMacroVal(1,iSpec,iSF) = AdaptBCMeanValues(8,iSpec,iSF) / REAL(SamplingIteration) &
                                                      / AdaptBCVolSurfaceFlux(iSpec,iSF) * Species(iSpec)%MacroParticleFactor
              END IF
              IF(AdaptBCMeanValues(7,iSpec,iSF).GT.1) THEN
                ! Calculate the average temperature
                AdaptBCAverageMacroVal(2,iSpec,iSF) = (AdaptBCMeanValues(7,iSpec,iSF)/(AdaptBCMeanValues(7,iSpec,iSF)-1.0)) &
                  *Species(iSpec)%MassIC / BoltzmannConst*(AdaptBCMeanValues(4,iSpec,iSF) - AdaptBCMeanValues(1,iSpec,iSF)**2  &
                                                         + AdaptBCMeanValues(5,iSpec,iSF) - AdaptBCMeanValues(2,iSpec,iSF)**2  &
                                                         + AdaptBCMeanValues(6,iSpec,iSF) - AdaptBCMeanValues(3,iSpec,iSF)**2)/3.
                ! Calculate the average pressure
                AdaptBCAverageMacroVal(3,iSpec,iSF) = AdaptBCAverageMacroVal(1,iSpec,iSF) * BoltzmannConst &
                                                      * AdaptBCAverageMacroVal(2,iSpec,iSF)
              END IF
            ELSE
            ! ================================================================
            ! Relaxation factor: updating the macro values with a certain percentage of the current sampled value
              ! Calculate the average number density
              IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
                AdaptBCAverageMacroVal(1,iSpec,iSF) = (1-RelaxationFactor) * AdaptBCAverageMacroVal(1,iSpec,iSF) &
                  + RelaxationFactor * AdaptBCMeanValues(8,iSpec,iSF)/REAL(SamplingIteration)/AdaptBCVolSurfaceFlux(iSpec,iSF)
              ELSE
                AdaptBCAverageMacroVal(1,iSpec,iSF) = (1-RelaxationFactor) * AdaptBCAverageMacroVal(1,iSpec,iSF) &
                  + RelaxationFactor * AdaptBCMeanValues(8,iSpec,iSF)/REAL(SamplingIteration)/AdaptBCVolSurfaceFlux(iSpec,iSF) &
                    * Species(iSpec)%MacroParticleFactor
              END IF
              IF(AdaptBCMeanValues(7,iSpec,iSF).GT.1) THEN
                ! Calculate the average temperature
                AdaptBCAverageMacroVal(2,iSpec,iSF) = (AdaptBCMeanValues(7,iSpec,iSF)/(AdaptBCMeanValues(7,iSpec,iSF)-1.0)) &
                  *Species(iSpec)%MassIC / BoltzmannConst*(AdaptBCMeanValues(4,iSpec,iSF) - AdaptBCMeanValues(1,iSpec,iSF)**2  &
                                                         + AdaptBCMeanValues(5,iSpec,iSF) - AdaptBCMeanValues(2,iSpec,iSF)**2  &
                                                         + AdaptBCMeanValues(6,iSpec,iSF) - AdaptBCMeanValues(3,iSpec,iSF)**2)/3.
                ! Calculate the average pressure
                AdaptBCAverageMacroVal(3,iSpec,iSF) = (1-RelaxationFactor) * AdaptBCAverageMacroVal(3,iSpec,iSF) &
                  + RelaxationFactor * AdaptBCAverageMacroVal(1,iSpec,iSF)*BoltzmannConst*AdaptBCAverageMacroVal(2,iSpec,iSF)
              END IF
            END IF
          END IF
        END DO
      END DO
    END IF
#if USE_MPI
    CALL MPI_BCAST(AdaptBCMeanValues,8*nSpecies*nSurfacefluxBCs, MPI_DOUBLE_PRECISION,0,MPI_COMM_PICLAS,iERROR)
    CALL MPI_BCAST(AdaptBCAverageMacroVal,3*nSpecies*nSurfacefluxBCs, MPI_DOUBLE_PRECISION,0,MPI_COMM_PICLAS,iERROR)
#endif /*USE_MPI*/
  END IF
END IF

IF(AdaptBCAverageValBC) THEN
! ================================================================
! 3) Populate the element local arrays
! BC-averaged values
  IF(CalcValues.OR.AdaptBCTruncAverage) THEN
    DO iSpec = 1,nSpecies
      DO iSF = 1,Species(iSpec)%nSurfacefluxBCs
        currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
        ! Copy values onto the element-wise array
        IF (BCdata_auxSF(currentBC)%SideNumber.GT.0) THEN
          IF(.NOT.initSampling) THEN
            DO iSide = 1, BCdata_auxSF(currentBC)%SideNumber
              BCSideID = BCdata_auxSF(currentBC)%SideList(iSide)
              ElemID = SideToElem(S2E_ELEM_ID,BCSideID)
              SampleElemID = AdaptBCMapElemToSample(ElemID)
              IF(SampleElemID.GT.0) AdaptBCMacroVal(1:3,SampleElemID,iSpec) = AdaptBCMeanValues(1:3,iSpec,iSF)
            END DO
          END IF
          DO iSide = 1, BCdata_auxSF(currentBC)%SideNumber
            BCSideID = BCdata_auxSF(currentBC)%SideList(iSide)
            ElemID = SideToElem(S2E_ELEM_ID,BCSideID)
            SampleElemID = AdaptBCMapElemToSample(ElemID)
            IF(SampleElemID.GT.0) THEN
              AdaptBCMacroVal(4,SampleElemID,iSpec) = AdaptBCAverageMacroVal(1,iSpec,iSF)
              AdaptBCMacroVal(6,SampleElemID,iSpec) = AdaptBCAverageMacroVal(3,iSpec,iSF)
            END IF
          END DO
        END IF    ! BCdata_auxSF(currentBC)%SideNumber.GT.0
      END DO      ! iSF = 1,Species(iSpec)%nSurfacefluxBCs
    END DO        ! iSpec = 1,nSpecies
  END IF          ! CalcValues.OR.AdaptBCTruncAverage
ELSE              ! .NOT. AdaptBCAverageValBC
! ================================================================
! 4) Calculation of the cell-local values
  IF(AdaptBCSampIter.GT.0) THEN
  ! ================================================================
  ! Sampling iteration: sampling for AdaptBCSampIter iterations, calculating the macro values and resetting sample OR
  ! AdaptBCTruncAverage: continuous average of the last AdaptBCSampIter iterations
    IF(CalcValues.OR.AdaptBCTruncAverage) THEN
      DO SampleElemID = 1,AdaptBCSampleElemNum
#if USE_LOADBALANCE
        CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
        ElemID = AdaptBCMapSampleToElem(SampleElemID)
        CNElemID = GetCNElemID(ElemID+offsetElem)
        DO iSpec = 1,nSpecies
          IF (AdaptBCSample(7,SampleElemID,iSpec).GT.0.0) THEN
            ! Calculate the average velocities
            AdaptBCSample(1:6,SampleElemID,iSpec) = AdaptBCSample(1:6,SampleElemID,iSpec) / AdaptBCSample(8,SampleElemID,iSpec)
            IF(.NOT.initSampling) THEN
              ! Compute flow velocity (during computation, not for the initial distribution, where the velocity from the ini is used)
              AdaptBCMacroVal(1:3,SampleElemID,iSpec) = AdaptBCSample(1:3,SampleElemID, iSpec)
            END IF
            ! number density
            IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
              AdaptBCMacroVal(4,SampleElemID,iSpec) = AdaptBCSample(8,SampleElemID,iSpec) / REAL(SamplingIteration) / ElemVolume_Shared(CNElemID)
            ELSE
              AdaptBCMacroVal(4,SampleElemID,iSpec) = AdaptBCSample(8,SampleElemID,iSpec) / REAL(SamplingIteration) / ElemVolume_Shared(CNElemID) &
                                                * Species(iSpec)%MacroParticleFactor
            END IF
            ! Calculation of the pressure
            IF(AdaptBCSample(7,SampleElemID,iSpec).GT.1) THEN
              ! instantaneous temperature WITHOUT 1/BoltzmannConst
              TTrans_TempFac = (AdaptBCSample(7,SampleElemID,iSpec)/(AdaptBCSample(7,SampleElemID,iSpec)-1.0)) &
                  *Species(iSpec)%MassIC*(AdaptBCSample(4,SampleElemID,iSpec) - AdaptBCSample(1,SampleElemID,iSpec)**2   &
                                        + AdaptBCSample(5,SampleElemID,iSpec) - AdaptBCSample(2,SampleElemID,iSpec)**2   &
                                        + AdaptBCSample(6,SampleElemID,iSpec) - AdaptBCSample(3,SampleElemID,iSpec)**2) / 3.
              ! pressure (BoltzmannConstant canceled out in temperature calculation)
              AdaptBCMacroVal(6,SampleElemID,iSpec)=AdaptBCMacroVal(4,SampleElemID,iSpec)*TTrans_TempFac
            END IF
          END IF    ! AdaptBCSample(7,SampleElemID,iSpec).GT.0.0
          ! Resetting sampled values
          AdaptBCSample(1:8,SampleElemID,iSpec) = 0.
        END DO      ! iSpec = 1,nSpecies
#if USE_LOADBALANCE
        CALL LBElemSplitTime(ElemID,tLBStart)
#endif /*USE_LOADBALANCE*/
      END DO        ! SampleElemID = 1,AdaptBCSampleElemNum
    END IF          ! CalcValues.OR.AdaptBCTruncAverage
  ELSE              ! AdaptBCSampIter.LE.0
  ! ================================================================
  ! Relaxation factor: updating the macro values with a certain percentage of the current sampled value
    DO SampleElemID = 1,AdaptBCSampleElemNum
#if USE_LOADBALANCE
      CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
      ElemID = AdaptBCMapSampleToElem(SampleElemID)
      CNElemID = GetCNElemID(ElemID+offsetElem)
      DO iSpec = 1,nSpecies
        IF (AdaptBCSample(7,SampleElemID,iSpec).GT.0.0) THEN
          ! Calculate the average velocities
          IF(.NOT.AdaptBCAverageValBC) AdaptBCSample(1:6,SampleElemID,iSpec) = AdaptBCSample(1:6,SampleElemID,iSpec) / AdaptBCSample(8,SampleElemID,iSpec)
          IF(.NOT.initSampling) THEN
            ! compute flow velocity (during computation, not for the initial distribution, where the velocity from the ini is used)
            AdaptBCMacroVal(1:3,SampleElemID,iSpec) = (1-RelaxationFactor)*AdaptBCMacroVal(1:3,SampleElemID,iSpec) &
                                                + RelaxationFactor*AdaptBCSample(1:3,SampleElemID, iSpec)
          END IF
          ! Calculation of the number density
          IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
            AdaptBCMacroVal(4,SampleElemID,iSpec) = (1-RelaxationFactor)*AdaptBCMacroVal(4,SampleElemID,iSpec) &
              + RelaxationFactor*AdaptBCSample(8,SampleElemID,iSpec) / ElemVolume_Shared(CNElemID)
          ELSE
            AdaptBCMacroVal(4,SampleElemID,iSpec) = (1-RelaxationFactor)*AdaptBCMacroVal(4,SampleElemID,iSpec) &
              + RelaxationFactor*AdaptBCSample(8,SampleElemID,iSpec) / ElemVolume_Shared(CNElemID)*Species(iSpec)%MacroParticleFactor
          END IF
          ! Calculation of the pressure
          IF (AdaptBCSample(7,SampleElemID,iSpec).GT.1.0) THEN
            ! Compute instantaneous temperature WITHOUT 1/BoltzmannConst
            TTrans_TempFac = (AdaptBCSample(7,SampleElemID,iSpec)/(AdaptBCSample(7,SampleElemID,iSpec)-1.0)) &
                              * Species(iSpec)%MassIC * (AdaptBCSample(4,SampleElemID,iSpec) - AdaptBCSample(1,SampleElemID,iSpec)**2 &
                                                        + AdaptBCSample(5,SampleElemID,iSpec) - AdaptBCSample(2,SampleElemID,iSpec)**2 &
                                                        + AdaptBCSample(6,SampleElemID,iSpec) - AdaptBCSample(3,SampleElemID,iSpec)**2) / 3.
            AdaptBCMacroVal(6,SampleElemID,iSpec) = (1-RelaxationFactor)*AdaptBCMacroVal(6,SampleElemID,iSpec) &
                                              + RelaxationFactor*AdaptBCMacroVal(4,SampleElemID,iSpec)*TTrans_TempFac
          END IF
        ELSE
          ! Relax the values towards zero
          AdaptBCMacroVal(1:7,SampleElemID,iSpec) = (1-RelaxationFactor)*AdaptBCMacroVal(1:7,SampleElemID,iSpec)
        END IF    ! AdaptBCSample(7,SampleElemID,iSpec).GT.0.0
        ! Resetting sampled values
        AdaptBCSample(1:8,SampleElemID,iSpec) = 0.
      END DO      ! iSpec = 1,nSpecies
#if USE_LOADBALANCE
      CALL LBElemSplitTime(ElemID,tLBStart)
#endif /*USE_LOADBALANCE*/
    END DO        ! SampleElemID = 1,AdaptBCSampleElemNum
  END IF          !  AdaptBCSampIter.GT.0
END IF            ! AdaptBCAverageValBC

END SUBROUTINE AdaptiveBCSampling


SUBROUTINE FinalizeParticleSamplingAdaptive(IsLoadBalance)
!----------------------------------------------------------------------------------------------------------------------------------!
!>
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Sampling_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
LOGICAL,INTENT(IN)              :: IsLoadBalance
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(AdaptBCAverage)
SDEALLOCATE(AdaptBCMacroVal)
SDEALLOCATE(AdaptBCSample)
SDEALLOCATE(AdaptBCMapSampleToElem)
SDEALLOCATE(AdaptBCMapElemToSample)
SDEALLOCATE(AdaptiveData)
SDEALLOCATE(AdaptBCAreaSurfaceFlux)
SDEALLOCATE(AdaptBCVolSurfaceFlux)
SDEALLOCATE(AdaptBCBackupVelocity)
IF(.NOT.IsLoadBalance) THEN
  SDEALLOCATE(AdaptBCPartNumOut)
  SDEALLOCATE(AdaptBCAverageMacroVal)
END IF

END SUBROUTINE FinalizeParticleSamplingAdaptive

END MODULE MOD_Particle_Sampling_Adapt
