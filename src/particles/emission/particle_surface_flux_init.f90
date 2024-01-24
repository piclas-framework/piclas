!==================================================================================================================================
! Copyright (c) 2021 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_Particle_SurfFlux_Init
!===================================================================================================================================
!> Module for the initialization of particle emission through the surface flux
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: DefineParametersParticleSurfaceFlux, InitializeParticleSurfaceflux
!===================================================================================================================================
CONTAINS

!==================================================================================================================================
!> Define parameters for particle surface flux
!==================================================================================================================================
SUBROUTINE DefineParametersParticleSurfaceFlux()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Particle Surface Flux")

CALL prms%CreateIntOption(      'Part-Species[$]-nSurfacefluxBCs', &
                                'Number of surface flux emissions per species', '0', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-Surfaceflux[$]-BC', &
                                'PartBound to be emitted from', numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Part-Species[$]-Surfaceflux[$]-velocityDistribution', &
                                'Specifying keyword for velocity distribution' , 'constant', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-VeloIC', &
                                'Velocity magnitude in meter/second', '0.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Species[$]-Surfaceflux[$]-VeloIsNormal', &
                                'VeloIC is normal to the surface and not a specific vector (=VeloVecIC)', &
                                '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-Surfaceflux[$]-VeloVecIC', &
                                'Normalized velocity vector' , '0.0 , 0.0 , 0.0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-MWTemperatureIC', &
                                'Temperature for Maxwell Distribution of particle velocities', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-PartDensity', &
                                'Number density (real particles per cubic meter)', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-EmissionCurrent', &
                                'Current over the emission surface in ampere (as an alternative to PartDensity for charged ' //&
                                'species). Velocity magnitude can be zero or above.', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-Massflow', &
                                'Mass flow over surface flux surface (as an alternative to PartDensity e.g. for outgassing. ' //&
                                'Velocity magnitude can be zero or above.', '0.', numberedmulti=.TRUE.)
! === Unclear/Deprecated
CALL prms%CreateLogicalOption(  'Part-Species[$]-Surfaceflux[$]-ReduceNoise' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Reduce stat. noise by global calc. of PartIns', '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Species[$]-Surfaceflux[$]-AcceptReject' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  ' Perform ARM for skewness of RefMap-positioning', '.TRUE.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-Surfaceflux[$]-ARM_DmaxSampleN' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Number of sample intervals in xi/eta for Dmax-calc.', '1', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'DoForceFreeSurfaceFlux' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Flag if the stage reconstruction uses a force' , '.FALSE.')
! === Circular inflow
CALL prms%CreateLogicalOption(  'Part-Species[$]-Surfaceflux[$]-CircularInflow', &
                                'Enables the utilization of a circular region as a surface flux on the selected boundary. '//&
                                'Only possible on surfaces, which are in xy, xz, or yz-planes.', '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-Surfaceflux[$]-axialDir', &
                                'Normal direction of the surface, where the circular inflow is defined: x = 1, y = 2, z = 3', &
                                numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-Surfaceflux[$]-origin', &
                                'Origin of circular inflow on the surface, where the coordinates depend on the axialDir:\n' //&
                                'x (=1): (y,z); y (=2): (z,x); z (=3): (x,y)', numberedmulti=.TRUE., no=2)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-rmax', &
                                'Maximum radius of the circular inflow to define a circle (rmin undefined) or a ring (rmin ' //&
                                'defined)', '1e21', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-rmin', &
                                'Minimal radius of the circular inflow to define a ring (rmax defined) or exclude an inner ' //&
                                'circle (rmax undefined)', '0.', numberedmulti=.TRUE.)
! === Adaptive surface flux types
CALL prms%CreateLogicalOption(  'Part-Species[$]-Surfaceflux[$]-Adaptive' &
                                      , 'Flag for the definition of adaptive boundary conditions', '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-Surfaceflux[$]-Adaptive-Type' &
                                , 'Define the type of the adaptive boundary condition. Options:\n' //&
                                  '(1) Const. static pressure inlet after Farbar & Boyd 2014 (Type 1)\n' //&
                                  '(2) Const. static pressure outlet after Farbar & Boyd 2014 (Type 1)\n' //&
                                  '(3) Const. massflow inlet after Farbar & Boyd 2014 (Type 2): Number of particles for regular '//&
                                  'surface flux is calculated with velocity and given mass flow. Requires an open BC.' //&
                                  '(4) Const. massflow inlet after Lei 2017 (cf_3): N_in = N_mdot + N_out (counting particles, '//&
                                  'which exist the domain through the adaptive BC).' &
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-Adaptive-Pressure' &
                                , 'Static pressure in [Pa] for the adaptive boundary conditions of type 1 and 2.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-Adaptive-Massflow' &
                                , 'Massflow in [kg/s] for the adaptive boundary conditions of type 3 and 4.', numberedmulti=.TRUE.)
! === Thermionic emission
CALL prms%CreateLogicalOption(  'Part-Species[$]-Surfaceflux[$]-ThermionicEmission', &
                                'Flag for the definition of a thermionic emission', '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Species[$]-Surfaceflux[$]-ThermionicEmission-SchottkyEffect', &
                                'Flag for the consideration of the Schottky effect: reduction of the work function (lowering the '//&
                                'escape barrier) due to an electric field (requires PIC HDG simulation)', '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-ThermionicEmission-WorkFunction', &
                                'Material specific work function W for the thermionic emission [eV]', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-ThermionicEmission-RichardsonConstant', &
                                'Material specific constant A*b for the thermionic emission [A / (cm^2 K^2)]', numberedmulti=.TRUE.)
END SUBROUTINE DefineParametersParticleSurfaceFlux


SUBROUTINE InitializeParticleSurfaceflux()
!===================================================================================================================================
! Init Particle Inserting via Surface Flux
!===================================================================================================================================
! Modules
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Mesh_Vars              ,ONLY: nBCSides, SideToElem, offsetElem, NGeo
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces_Vars ,ONLY: BCdata_auxSF, BezierSampleN, SurfMeshSubSideData, SurfMeshSideAreas, tBCdata_auxSFRadWeight
USE MOD_Particle_Surfaces_Vars ,ONLY: SurfFluxSideSize, TriaSurfaceFlux
USE MOD_Particle_Surfaces      ,ONLY: GetBezierSampledAreas
USE MOD_Particle_Vars          ,ONLY: Species, nSpecies, DoSurfaceFlux
USE MOD_Particle_Vars          ,ONLY: UseCircularInflow, DoForceFreeSurfaceFlux
USE MOD_Particle_Vars          ,ONLY: VarTimeStep
USE MOD_Particle_Sampling_Vars ,ONLY: UseAdaptiveBC
USE MOD_Restart_Vars           ,ONLY: DoRestart, RestartTime
USE MOD_DSMC_Vars              ,ONLY: AmbiPolarSFMapping, DSMC, useDSMC
USE MOD_Particle_Vars          ,ONLY: Symmetry
#if USE_MPI
USE MOD_Particle_Vars          ,ONLY: DoPoissonRounding, DoTimeDepInflow
#endif /*USE_MPI*/
#ifdef CODE_ANALYZE
USE MOD_Particle_Vars          ,ONLY: CountCircInflowType
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Local variable declaration
INTEGER               :: iSpec,iSF,SideID,BCSideID,iSide,ElemID,iLocSide,iSample,jSample,currentBC, MaxSF, iSFElec
INTEGER               :: iCopy1, iCopy2, iCopy3, MaxSurfacefluxBCs,nDataBC
REAL                  :: tmp_SubSideDmax(SurfFluxSideSize(1),SurfFluxSideSize(2))
REAL                  :: tmp_SubSideAreas(SurfFluxSideSize(1),SurfFluxSideSize(2))
REAL                  :: tmp_BezierControlPoints2D(2,0:NGeo,0:NGeo,SurfFluxSideSize(1),SurfFluxSideSize(2))
REAL                  :: VFR_total, RestartTimeVar
TYPE(tBCdata_auxSFRadWeight),ALLOCATABLE          :: BCdata_auxSFTemp(:)
#if USE_MPI
REAL                  :: totalAreaSF_global
#endif
!===================================================================================================================================
ALLOCATE(SurfMeshSubSideData(SurfFluxSideSize(1),SurfFluxSideSize(2),1:nBCSides),SurfMeshSideAreas(1:nBCSides))
SurfMeshSideAreas=0.
! global calculations for sampling the faces for area and vector calculations (checks the integration with CODE_ANALYZE)
CALL BCSurfMeshSideAreasandNormals()

UseCircularInflow=.FALSE.
UseAdaptiveBC=.FALSE.
MaxSurfacefluxBCs=0
nDataBC=0
DoSurfaceFlux=.FALSE.

!-- 1.: read/prepare parameters and determine nec. BCs
CALL ReadInAndPrepareSurfaceFlux(MaxSurfacefluxBCs, nDataBC)

#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoPoissonRounding,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_PICLAS,iError) !set T if this is for all procs
CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoTimeDepInflow,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_PICLAS,iError) !set T if this is for all procs
#endif /*USE_MPI*/

CALL CreateSideListAndFinalizeAreasSurfFlux(nDataBC, BCdata_auxSFTemp)

#ifdef CODE_ANALYZE
IF (UseCircularInflow) THEN
  ALLOCATE(CountCircInflowType(1:3,1:MaxSurfacefluxBCs,1:nSpecies))
  CountCircInflowType = 0
END IF
#endif

!-- 3.: initialize Surfaceflux-specific data
DO iSpec=1,nSpecies
  DO iSF=1,Species(iSpec)%nSurfacefluxBCs
    currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
    IF (BCdata_auxSF(currentBC)%SideNumber.EQ.-1) THEN
      CALL abort(__STAMP__,'ERROR in ParticleSurfaceflux: Someting is wrong with SideNumber of BC ',currentBC)
    END IF
    ! Loop over sides on the surface flux
    DO iSide=1,BCdata_auxSF(currentBC)%SideNumber
      BCSideID=BCdata_auxSF(currentBC)%SideList(iSide)
      ElemID = SideToElem(S2E_ELEM_ID,BCSideID)
      iLocSide = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
      SideID=GetGlobalNonUniqueSideID(offsetElem+ElemID,iLocSide)
      ! Calculate the total area of the surface flux
      IF (Species(iSpec)%Surfaceflux(iSF)%AcceptReject) THEN
        CALL GetBezierSampledAreas(SideID=SideID,BezierSampleN=BezierSampleN &
          ,BezierSurfFluxProjection_opt=.NOT.Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal &
          ,SurfMeshSubSideAreas=tmp_SubSideAreas,DmaxSampleN_opt=Species(iSpec)%Surfaceflux(iSF)%ARM_DmaxSampleN &
          ,Dmax_opt=tmp_SubSideDmax,BezierControlPoints2D_opt=tmp_BezierControlPoints2D)
      ELSE IF (.NOT.TriaSurfaceFlux) THEN
        CALL GetBezierSampledAreas(SideID=SideID,BezierSampleN=BezierSampleN &
          ,BezierSurfFluxProjection_opt=.NOT.Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal,SurfMeshSubSideAreas=tmp_SubSideAreas)
      ELSE ! TriaSurfaceFlux
        DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
          tmp_SubSideAreas(iSample,jSample)=SurfMeshSubSideData(iSample,jSample,BCSideID)%area
          IF(.NOT.Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
            Species(iSpec)%Surfaceflux(iSF)%totalAreaSF = Species(iSpec)%Surfaceflux(iSF)%totalAreaSF &
                                                          + SurfMeshSubSideData(iSample,jSample,BCSideID)%area
          END IF
        END DO; END DO
      END IF
      ! Initialize circular inflow (determine if elements are (partially) inside/outside)
      IF (Species(iSpec)%Surfaceflux(iSF)%CircularInflow) CALL DefineCircInflowRejectType(iSpec, iSF, iSide)
      ! Initialize the volume flow rate
      CALL InitVolumeFlowRate(iSpec, iSF, iSide, tmp_SubSideAreas, BCdata_auxSFTemp)
      ! Initialize acceptance-rejection on SF
      IF (Species(iSpec)%Surfaceflux(iSF)%AcceptReject) THEN
        DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
          Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Dmax = tmp_SubSideDmax(iSample,jSample)
          IF (.NOT.Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
            ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample &
                                                                        ,iSide)%BezierControlPoints2D(1:2,0:NGeo,0:NGeo))
            DO iCopy1=0,NGeo; DO iCopy2=0,NGeo; DO iCopy3=1,2
              Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample &
                                                                  ,iSide)%BezierControlPoints2D(iCopy3,iCopy2,iCopy1) &
                = tmp_BezierControlPoints2D(iCopy3,iCopy2,iCopy1,iSample,jSample)
            END DO; END DO; END DO
          END IF !.NOT.VeloIsNormal
        END DO; END DO !jSample=1,SurfFluxSideSize(2); iSample=1,SurfFluxSideSize(1)
      END IF
    END DO ! iSide
    !--- 3b: ReduceNoise initialization
    IF (Species(iSpec)%Surfaceflux(iSF)%ReduceNoise) CALL InitReduceNoiseSF(iSpec, iSF)
    ! Calculate the total area per surface flux
#if USE_MPI
    IF(.NOT.Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
      totalAreaSF_global = 0.0
      CALL MPI_ALLREDUCE(Species(iSpec)%Surfaceflux(iSF)%totalAreaSF,totalAreaSF_global,1, &
                          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_PICLAS,IERROR)
      Species(iSpec)%Surfaceflux(iSF)%totalAreaSF = totalAreaSF_global
    END IF
#endif
    ! Inserting particles through a rate instead of particle density. This assumes that the volume flow rate (VFR) has been replaced
    ! by the local area (in InitVolumeFlowRate).
    IF(Species(iSpec)%Surfaceflux(iSF)%UseEmissionCurrent) THEN
      ! Store the current as particles per second per square meter in the particle density variable.
      IF(Species(iSpec)%Surfaceflux(iSF)%ThermionicEmission) THEN
        ! Thermionic emission: Richardson-Dushman equation gives directly the current density [A/m2]
        Species(iSpec)%Surfaceflux(iSF)%PartDensity = Species(iSpec)%Surfaceflux(iSF)%EmissionCurrent / ABS(Species(iSpec)%ChargeIC)
      ELSE
        Species(iSpec)%Surfaceflux(iSF)%PartDensity = Species(iSpec)%Surfaceflux(iSF)%EmissionCurrent &
                                                      / (ABS(Species(iSpec)%ChargeIC) * Species(iSpec)%Surfaceflux(iSF)%totalAreaSF)
      END IF
    END IF
    IF(Species(iSpec)%Surfaceflux(iSF)%UseMassflow) THEN
      ! Store the mass flow as particles per second per square meter in the particle density variable.
      Species(iSpec)%Surfaceflux(iSF)%PartDensity = Species(iSpec)%Surfaceflux(iSF)%Massflow &
                                                    / (Species(iSpec)%MassIC * Species(iSpec)%Surfaceflux(iSF)%totalAreaSF)
    END IF
    ! Output of the number of sides for circular inflow (only if code was compiled with CODE_ANALYZE = TRUE
#ifdef CODE_ANALYZE
    IF (BCdata_auxSF(currentBC)%SideNumber.GT.0 .AND. Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
      IPWRITE(*,'(I4,A,2(1X,I0),A,3(1X,I0))') ' For Surfaceflux/Spec',iSF,iSpec,' are nType0,1,2: ', &
        CountCircInflowType(1,iSF,iSpec),CountCircInflowType(2, iSF,iSpec), CountCircInflowType(3, iSF,iSpec)
    END IF
#endif /*CODE_ANALYZE*/
  END DO !iSF
END DO !iSpec

#ifdef CODE_ANALYZE
SDEALLOCATE(CountCircInflowType)
#endif
! Deallocate auxiliary variable container (no pointers used inside container)
SDEALLOCATE(BCdata_auxSFTemp)

! Setting variables required after a restart
IF(DoRestart) THEN
  DO iSpec=1,nSpecies
    ! Species-specific time step
    IF(VarTimeStep%UseSpeciesSpecific) THEN
      RestartTimeVar = RestartTime * Species(iSpec)%TimeStepFactor
    ELSE
      RestartTimeVar = RestartTime
    END IF
    DO iSF = 1, Species(iSpec)%NumberOfInits
      Species(iSpec)%Init(iSF)%InsertedParticle = INT(Species(iSpec)%Init(iSF)%ParticleNumber * RestartTimeVar,8)
    END DO
    DO iSF = 1, Species(iSpec)%nSurfacefluxBCs
      IF (Species(iSpec)%Surfaceflux(iSF)%ReduceNoise) THEN
        VFR_total = Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcsTotal !proc global total (for non-root: dummy!!!)
      ELSE
        VFR_total = Species(iSpec)%Surfaceflux(iSF)%VFR_total               !proc local total
      END IF

      Species(iSpec)%Surfaceflux(iSF)%InsertedParticle = INT(Species(iSpec)%Surfaceflux(iSF)%PartDensity * RestartTimeVar &
        / Species(iSpec)%MacroParticleFactor * VFR_total,8)
    END DO
  END DO
END IF

! Initialize the adaptive surface flux of type 4 (performed at the end, requires the totalAreaSF for zero mass flow weighting)
DO iSpec=1,nSpecies
  DO iSF=1,Species(iSpec)%nSurfacefluxBCs
    IF(Species(iSpec)%Surfaceflux(iSF)%Adaptive) THEN
      IF(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType.EQ.4) THEN
        currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
        IF(Symmetry%Axisymmetric) CALL abort(__STAMP__, 'ERROR: AdaptiveType = 4 is not implemented with axisymmetric simulations!')
        ! Weighting factor to account for a
        ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%ConstMassflowWeight(1:SurfFluxSideSize(1),1:SurfFluxSideSize(2), &
                  1:BCdata_auxSF(currentBC)%SideNumber))
        Species(iSpec)%Surfaceflux(iSF)%ConstMassflowWeight = 0.0
        ! In case the mass flow is set to zero, the weighting factor can be determined once initially based on area ratio
        IF(ALMOSTEQUAL(Species(iSpec)%Surfaceflux(iSF)%AdaptiveMassflow,0.)) CALL CalcConstMassflowWeightForZeroMassFlow(iSpec,iSF)
        ! Circular inflow in combination with AdaptiveType = 4 requires the partial circle area per tria side
        IF(Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
          IF(Symmetry%Axisymmetric) CALL abort(__STAMP__, 'ERROR: Circular inflow is not implemented with axisymmetric simulations!')
          ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%CircleAreaPerTriaSide(1:SurfFluxSideSize(1),1:SurfFluxSideSize(2), &
                  1:BCdata_auxSF(currentBC)%SideNumber))
          Species(iSpec)%Surfaceflux(iSF)%CircleAreaPerTriaSide = 0.0
          CALL CalcCircInflowAreaPerTria(iSpec,iSF)
        END IF
      END IF
    END IF
  END DO  ! iSF=1,Species(iSpec)%nSurfacefluxBCs
END DO    ! iSpec=1,nSpecies

#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoSurfaceFlux,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_PICLAS,iError) !set T if at least 1 proc have SFs
#endif  /*USE_MPI*/
IF (.NOT.DoSurfaceFlux) THEN !-- no SFs defined
  LBWRITE(*,*) 'WARNING: No Sides for SurfacefluxBCs found! DoSurfaceFlux is now disabled!'
END IF
DoForceFreeSurfaceFlux = GETLOGICAL('DoForceFreeSurfaceFlux','.FALSE.')

IF (useDSMC) THEN
  IF (DSMC%DoAmbipolarDiff) THEN
    MaxSF = 0
    DO iSpec = 1,nSpecies
      MaxSF = MAX(MaxSF,Species(iSpec)%nSurfacefluxBCs)
    END DO

    ALLOCATE(AmbiPolarSFMapping(nSpecies,MaxSF))
    AmbiPolarSFMapping = -1
    DO iSpec = 1,nSpecies
      IF (Species(iSpec)%ChargeIC.LE.0.0) CYCLE
      DO iSF = 1,Species(iSpec)%nSurfacefluxBCs
        DO iSFElec = 1,Species(DSMC%AmbiDiffElecSpec)%nSurfacefluxBCs
          IF(Species(iSpec)%Surfaceflux(iSF)%BC.EQ.Species(DSMC%AmbiDiffElecSpec)%Surfaceflux(iSFElec)%BC) THEN
            AmbiPolarSFMapping(iSpec,iSF) = iSFElec
          END IF
        END DO
        IF(AmbiPolarSFMapping(iSpec,iSF).EQ.-1) CALL abort(__STAMP__,&
            'ERROR: No corresponding Electron Surface Flux found for Species: ',IntInfoOpt=iSpec)
      END DO
    END DO
  END IF
END IF

END SUBROUTINE InitializeParticleSurfaceflux


SUBROUTINE ReadInAndPrepareSurfaceFlux(MaxSurfacefluxBCs, nDataBC)
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst, Pi
USE MOD_Particle_Vars          ,ONLY: nSpecies, Species, UseVarTimeStep, VarTimeStep, DoPoissonRounding, DoTimeDepInflow
USE MOD_Particle_Vars          ,ONLY: UseCircularInflow
USE MOD_Particle_Sampling_Vars ,ONLY: UseAdaptiveBC
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound,nPartBound
USE MOD_DSMC_Vars              ,ONLY: useDSMC, BGGas, RadialWeighting
USE MOD_Particle_Surfaces_Vars ,ONLY: BCdata_auxSF, BezierSampleN, TriaSurfaceFlux
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_Mesh_Vars              ,ONLY: NGeo
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(INOUT) :: MaxSurfacefluxBCs, nDataBC
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(42)         :: hilf, hilf2, hilf3
INTEGER               :: iSpec, iSF
!===================================================================================================================================
DO iSpec=1,nSpecies
  WRITE(UNIT=hilf,FMT='(I0)') iSpec
  ! Read-in the number of surface flux BCs per species
  Species(iSpec)%nSurfacefluxBCs = GETINT('Part-Species'//TRIM(hilf)//'-nSurfacefluxBCs')
  ! Skip the remainder of the loop, if no surface fluxes have been defined
  IF (Species(iSpec)%nSurfacefluxBCs.EQ.0) CYCLE
  ! Sanity check: Background species cannot have a surface flux emission
  IF (useDSMC.AND.(BGGas%NumberOfSpecies.GT.0)) THEN
    IF (BGGas%BackgroundSpecies(iSpec)) THEN
      CALL abort(__STAMP__,'ERROR: Surface flux is not implemented for background gas species!')
    END IF
  END IF
  ! Determine the maximum number of surface flux BCs for all species
  MaxSurfacefluxBCs=MAX(MaxSurfacefluxBCs,Species(iSpec)%nSurfacefluxBCs)
  ! Allocate the species-specific surface flux array
  ALLOCATE(Species(iSpec)%Surfaceflux(1:Species(iSpec)%nSurfacefluxBCs))
  ! Initialize Surfaceflux to BC mapping
  Species(iSpec)%Surfaceflux(:)%BC=-1
  ! Loop over the surface flux BCs
  DO iSF=1,Species(iSpec)%nSurfacefluxBCs
    WRITE(UNIT=hilf2,FMT='(I0)') iSF
    hilf2=TRIM(hilf)//'-Surfaceflux'//TRIM(hilf2)
    ASSOCIATE(SF => Species(iSpec)%Surfaceflux(iSF))
    SF%BC = GETINT('Part-Species'//TRIM(hilf2)//'-BC')
    ! Sanity check: BC index must be within the number of defined particle boundary conditions
    IF (SF%BC.LT.1 .OR. SF%BC.GT.nPartBound) CALL abort(__STAMP__, 'SurfacefluxBCs must be between 1 and nPartBound!')
    ! Default surface flux type (constant particle number per side)
    SF%Type = 0
    SF%AdaptiveType = 0
    ! Initialize SF-specific variables
    WRITE(UNIT=hilf2,FMT='(I0)') iSF
    hilf2=TRIM(hilf)//'-Surfaceflux'//TRIM(hilf2)
    SF%InsertedParticle = 0
    SF%InsertedParticleSurplus = 0
    SF%VFR_total = 0
    SF%VFR_total_allProcsTotal = 0
    SF%totalAreaSF = 0.
    ! SideNumber has not been set yet
    IF (BCdata_auxSF(SF%BC)%SideNumber.EQ. -1) THEN
      BCdata_auxSF(SF%BC)%SideNumber=0
      nDataBC=nDataBC+1
    END IF
    SF%velocityDistribution  = TRIM(GETSTR('Part-Species'//TRIM(hilf2)//'-velocityDistribution','constant'))
    IF (TRIM(SF%velocityDistribution).NE.'constant' .AND. TRIM(SF%velocityDistribution).NE.'maxwell' .AND. &
        TRIM(SF%velocityDistribution).NE.'maxwell_lpn') THEN
      CALL abort(__STAMP__,'Only constant or maxwell-like velocity distributions implemented for surface flux!')
    END IF
    SF%VeloIC                = GETREAL('Part-Species'//TRIM(hilf2)//'-VeloIC')
    SF%VeloIsNormal          = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-VeloIsNormal')
    IF (SF%VeloIsNormal) THEN
      SF%CircularInflow=.FALSE.
    ELSE
      SF%VeloVecIC          = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-VeloVecIC',3)
      SF%CircularInflow     = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-CircularInflow')
      IF(SF%CircularInflow) THEN
        UseCircularInflow = .TRUE.
        SF%dir(1)         = GETINT('Part-Species'//TRIM(hilf2)//'-axialDir')
        IF (SF%dir(1).EQ.1) THEN
          SF%dir(2)=2
          SF%dir(3)=3
        ELSE IF (SF%dir(1).EQ.2) THEN
          SF%dir(2)=3
          SF%dir(3)=1
        ELSE IF (SF%dir(1).EQ.3) THEN
          SF%dir(2)=1
          SF%dir(3)=2
        ELSE
          CALL abort(__STAMP__,'ERROR in Surface Flux: axialDir for circular inflow must be between 1 and 3!')
        END IF
        SF%origin   = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-origin',2)
        WRITE(UNIT=hilf3,FMT='(E16.8)') HUGE(SF%rmax)
        SF%rmax     = GETREAL('Part-Species'//TRIM(hilf2)//'-rmax',TRIM(hilf3))
        SF%rmin     = GETREAL('Part-Species'//TRIM(hilf2)//'-rmin')
        ! Total area of surface flux
        SF%totalAreaSF = Pi*(SF%rmax*SF%rmax - SF%rmin*SF%rmin)
      END IF
    END IF !.NOT.VeloIsNormal
    IF (.NOT.SF%VeloIsNormal) THEN
      !--- normalize VeloVecIC
      IF (.NOT. ALL(SF%VeloVecIC(:).EQ.0.)) SF%VeloVecIC = SF%VeloVecIC / SQRT(DOT_PRODUCT(SF%VeloVecIC,SF%VeloVecIC))
    END IF
    SF%MWTemperatureIC       = GETREAL('Part-Species'//TRIM(hilf2)//'-MWTemperatureIC')
    SF%PartDensity           = GETREAL('Part-Species'//TRIM(hilf2)//'-PartDensity')
    SF%EmissionCurrent       = GETREAL('Part-Species'//TRIM(hilf2)//'-EmissionCurrent')
    SF%Massflow              = GETREAL('Part-Species'//TRIM(hilf2)//'-Massflow')
    SF%SampledMassflow       = 0.
    ! === Sanity checks & compatibility
    IF(SF%PartDensity.GT.0.) THEN
      IF(SF%EmissionCurrent.GT.0..OR.SF%Massflow.GT.0.) THEN
        CALL abort(__STAMP__,'ERROR in Surface Flux: PartDensity and EmissionCurrent/Massflow cannot be both above 0!')
      END IF
    END IF
    IF(SF%EmissionCurrent.GT.0.) THEN
      SF%UseEmissionCurrent = .TRUE.
      IF(SF%Massflow.GT.0.) THEN
        CALL abort(__STAMP__,'ERROR in Surface Flux: Mass flow and emission current cannot be defined at the same time!')
      END IF
      IF(Species(iSpec)%ChargeIC.EQ.0.) THEN
        CALL abort(__STAMP__,'ERROR in Surface Flux: Using the emission current is only possible for charged species!')
      END IF
      IF(TRIM(SF%velocityDistribution).EQ.'constant') THEN
        CALL abort(__STAMP__,'ERROR in Surface Flux: Constant velocity distribution is not supported for the emission current!')
      END IF
    ELSE
      SF%UseEmissionCurrent = .FALSE.
    END IF
    IF(SF%Massflow.GT.0.) THEN
      SF%UseMassflow = .TRUE.
      IF(TRIM(SF%velocityDistribution).EQ.'constant') THEN
        CALL abort(__STAMP__,'ERROR in Surface Flux: Constant velocity distribution is not supported for the mass flow!')
      END IF
    ELSE
      SF%UseMassflow = .FALSE.
    END IF
    ! === ReduceNoise & AcceptReject ===============================================================================================
    SF%ReduceNoise           = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-ReduceNoise')
    IF (DoPoissonRounding .AND. SF%ReduceNoise) THEN
      LBWRITE(*,*)'WARNING: Poisson sampling not possible for noise reduction of surfacefluxes:'
      LBWRITE(*,*)'switching now to Random rounding...'
      DoPoissonRounding   = .FALSE.
    END IF
    IF (DoTimeDepInflow .AND. SF%ReduceNoise) THEN
      LBWRITE(*,*)'WARNING: Time-dependent inflow is not possible for noise reduction of surfacefluxes:'
      LBWRITE(*,*)'switching now to Random rounding...'
      DoTimeDepInflow   = .FALSE.
    END IF
    IF (TriaSurfaceFlux) THEN
      SF%AcceptReject = .FALSE.
    ELSE
      SF%AcceptReject = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-AcceptReject','.TRUE.')
    END IF
    IF (SF%AcceptReject .AND. BezierSampleN.GT.1) THEN
      LBWRITE(*,*)'WARNING: BezierSampleN > 0 may not be necessary as ARM is used for SurfaceFlux!'
    ELSE IF (.NOT.SF%AcceptReject .AND. BezierSampleN.LE.NGeo .AND. .NOT.TriaSurfaceFlux) THEN
      LBWRITE(*,*)'WARNING: The choosen small BezierSampleN (def.: NGeo) might result in inhom. SurfFluxes without ARM!'
    END IF
    IF (SF%AcceptReject) THEN
      WRITE( hilf3, '(I0.2)') NGeo*NGeo*NGeo !1 for linear elements, this is an arbitray estimation for higher N!
      SF%ARM_DmaxSampleN = GETINT('Part-Species'//TRIM(hilf2)//'-ARM_DmaxSampleN',hilf3)
    ELSE
      SF%ARM_DmaxSampleN = 0
    END IF
    ! === Set the surface flux type
    IF(DoPoissonRounding) SF%Type = 3
    IF(DoTimeDepInflow)   SF%Type = 4
    IF(RadialWeighting%DoRadialWeighting) SF%Type = 2
    ! === ADAPTIVE BC ==============================================================================================================
    SF%Adaptive         = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-Adaptive')
    IF(SF%Adaptive) THEN
      UseAdaptiveBC  = .TRUE.
      SF%Type = 1
      IF(TrackingMethod.EQ.REFMAPPING) THEN
        CALL abort(__STAMP__,'ERROR: Adaptive surface flux boundary conditions are not implemented with RefMapping!')
      END IF
      IF(UseVarTimeStep.AND..NOT.VarTimeStep%UseDistribution) THEN
        CALL abort(__STAMP__,'ERROR: Adaptive surface flux boundary conditions are not implemented with variable time step!')
      END IF
      SF%AdaptiveType         = GETINT('Part-Species'//TRIM(hilf2)//'-Adaptive-Type')
      SELECT CASE(SF%AdaptiveType)
      ! Farbar2014 - Case 1: Inlet Type 1, constant pressure and temperature
      !              Case 2: Outlet Type 1, constant pressure
      !              Case 3: Inlet Type 2, constant mass flow and temperature
      ! Lei2017    - Case 4: cf_3, constant mass flow and temperature N through mass flow and particles out
      CASE(1,2)
        SF%AdaptivePressure  = GETREAL('Part-Species'//TRIM(hilf2)//'-Adaptive-Pressure')
        SF%PartDensity       = SF%AdaptivePressure / (BoltzmannConst * SF%MWTemperatureIC)
      CASE(3,4)
        SF%AdaptiveMassflow  = GETREAL('Part-Species'//TRIM(hilf2)//'-Adaptive-Massflow')
        IF(ALMOSTEQUAL(SF%AdaptiveMassflow,0.).AND.SF%AdaptiveMassflow.NE.0.) THEN
          CALL abort(__STAMP__,'ERROR in adaptive inlet: given mass flow is within machine tolerance!')
        END IF
        IF(SF%VeloIC.LE.0.0) THEN
          CALL abort(__STAMP__,'ERROR in adaptive inlet: positive initial guess of velocity for Type 3/Type 4 condition required!')
        END IF
      END SELECT
      ! Sanity check: regular surface flux must be on an open BC
      IF (PartBound%TargetBoundCond(SF%BC).EQ.PartBound%ReflectiveBC) THEN
        IF(.NOT.SF%CircularInflow.AND.(SF%AdaptiveType.NE.4)) THEN
          CALL abort(__STAMP__&
            ,'ERROR in adaptive surface flux: using a reflective BC without circularInflow is only allowed for Type 4!')
        END IF
      END IF
    END IF
    ! === THERMIONIC EMISSION ======================================================================================================
    SF%ThermionicEmission = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-ThermionicEmission')
    IF(SF%ThermionicEmission) CALL ReadInThermionicEmission(iSpec,iSF)
    ! ==============================================================================================================================
    END ASSOCIATE
  END DO !iSF
END DO ! iSpec

END SUBROUTINE ReadInAndPrepareSurfaceFlux


SUBROUTINE BCSurfMeshSideAreasandNormals()
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Surfaces_Vars ,ONLY: SurfFluxSideSize, SurfMeshSubSideData, BezierSampleN, SurfMeshSideAreas, TriaSurfaceFlux
USE MOD_Mesh_Vars              ,ONLY: nBCSides, offsetElem, SideToElem
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces      ,ONLY: GetBezierSampledAreas, CalcNormAndTangTriangle
#if CODE_ANALYZE && USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*CODE_ANALYZE && USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: totalArea
INTEGER                 :: BCSideID, ElemID, iLocSide, SideID, jSample, iSample
REAL                    :: tmp_SubSideAreas(SurfFluxSideSize(1),SurfFluxSideSize(2))
REAL                    :: tmp_Vec_nOut(3,SurfFluxSideSize(1),SurfFluxSideSize(2))
REAL                    :: tmp_Vec_t1(3,SurfFluxSideSize(1),SurfFluxSideSize(2))
REAL                    :: tmp_Vec_t2(3,SurfFluxSideSize(1),SurfFluxSideSize(2))
!===================================================================================================================================
totalArea=0.
DO BCSideID=1,nBCSides
  ElemID = SideToElem(S2E_ELEM_ID,BCSideID)
  iLocSide = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
  SideID=GetGlobalNonUniqueSideID(offsetElem+ElemID,iLocSide)
  IF (TriaSurfaceFlux) THEN
    IF (SurfFluxSideSize(1).NE.1 .OR. SurfFluxSideSize(2).NE.2) CALL abort(__STAMP__,&
      'SurfFluxSideSize must be 1,2 for TriaSurfaceFlux!')
    DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
      CALL CalcNormAndTangTriangle(SideID=SideID,nVec=tmp_Vec_nOut(:,iSample,jSample) &
        ,tang1=tmp_Vec_t1(:,iSample,jSample) &
        ,tang2=tmp_Vec_t2(:,iSample,jSample) &
        ,area=tmp_SubSideAreas(iSample,jSample) &
        ,TriNum=jSample)
      SurfMeshSideAreas(BCSideID)=SurfMeshSideAreas(BCSideID)+tmp_SubSideAreas(iSample,jSample)
    END DO; END DO
  ELSE
    IF (ANY(SurfFluxSideSize.NE.BezierSampleN)) CALL abort(__STAMP__,&
      'SurfFluxSideSize must be BezierSampleN,BezierSampleN for .NOT.TriaSurfaceFlux!')
    CALL GetBezierSampledAreas(SideID=SideID &
      ,BezierSampleN=BezierSampleN &
      ,SurfMeshSubSideAreas=tmp_SubSideAreas &
      ,SurfMeshSideArea_opt=SurfMeshSideAreas(BCSideID) &
      ,SurfMeshSubSideVec_nOut_opt=tmp_Vec_nOut &
      ,SurfMeshSubSideVec_t1_opt=tmp_Vec_t1 &
      ,SurfMeshSubSideVec_t2_opt=tmp_Vec_t2)
  END IF
  totalArea=totalArea+SurfMeshSideAreas(BCSideID)
  DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
    SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn=-tmp_Vec_nOut(:,iSample,jSample)
    SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t1=tmp_Vec_t1(:,iSample,jSample)
    SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t2=tmp_Vec_t2(:,iSample,jSample)
    SurfMeshSubSideData(iSample,jSample,BCSideID)%area=tmp_SubSideAreas(iSample,jSample)
  END DO; END DO
END DO
#ifdef CODE_ANALYZE
#if USE_LOADBALANCE
IF(.NOT.PerformLoadBalance)THEN
#endif /*USE_LOADBALANCE*/
  IPWRITE(*,*)" ===== TOTAL AREA (all BCsides) ====="
  IPWRITE(*,*)"totalArea       = ",totalArea
  IPWRITE(*,*)"totalArea/(pi) = ",totalArea/(ACOS(-1.))
  IPWRITE(*,*)" ===== TOTAL AREA (all BCsides) ====="
#if USE_LOADBALANCE
END IF ! PerformLoadBalance
#endif /*USE_LOADBALANCE*/
#endif /*CODE_ANALYZE*/
END SUBROUTINE BCSurfMeshSideAreasandNormals

SUBROUTINE CreateSideListAndFinalizeAreasSurfFlux(nDataBC, BCdata_auxSFTemp)
!===================================================================================================================================
! SideList for SurfaceFlux in BCdata_auxSF is created. Furthermore, the side areas are corrected for the 1D/2D case and finally
! communicated.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Mesh_Tools             ,ONLY: GetCNSideID
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces_Vars ,ONLY: BCdata_auxSF, SurfMeshSubSideData, SurfFluxSideSize, TriaSurfaceFlux, tBCdata_auxSFRadWeight
USE MOD_Particle_Surfaces_Vars ,ONLY: SideType
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound,nPartBound
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_Mesh_Vars              ,ONLY: nBCSides, offsetElem, BC, SideToElem
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO, ElemMidPoint_Shared, SideInfo_Shared
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID
USE MOD_Particle_Vars          ,ONLY: UseCircularInflow, Species, DoSurfaceFlux, nSpecies, Symmetry
USE MOD_DSMC_Symmetry          ,ONLY: DSMC_1D_CalcSymmetryArea, DSMC_2D_CalcSymmetryArea, DSMC_2D_CalcSymmetryAreaSubSides
USE MOD_DSMC_Vars              ,ONLY: RadialWeighting
USE MOD_Particle_Surfaces      ,ONLY: CalcNormAndTangTriangle
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                           :: nDataBC
TYPE(tBCdata_auxSFRadWeight), ALLOCATABLE, INTENT(INOUT)        :: BCdata_auxSFTemp(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: TmpMapToBC(1:nDataBC), TmpSideStart(1:nDataBC), TmpSideNumber(1:nDataBC), TmpSideEnd(1:nDataBC)
! PartBC, Start of Linked List for Sides in SurfacefluxBC, Number of Particles in Sides in SurfacefluxBC, End of Linked List for Sides in SurfacefluxBC
INTEGER               :: TmpSideNext(1:nBCSides) !Next: Sides of diff. BCs ar not overlapping!
INTEGER               :: countDataBC,iBC,BCSideID,currentBC,iSF,iCount,iLocSide,SideID,CNSideID,iPartBound
INTEGER               :: ElemID,CNElemID,GlobalElemID
INTEGER               :: iSample,jSample,iSpec,iSub
#if USE_MPI
REAL, ALLOCATABLE     :: areasLoc(:),areasGlob(:)
#endif /*USE_MPI*/
REAL                  :: ymax,ymin,yMaxTemp,yMinTemp
!===================================================================================================================================
!-- 2.: create Side lists for applicable BCs
!--- 2a: temporary (linked) lists
TmpMapToBC = 0; TmpSideStart = 0; TmpSideNumber = 0; TmpSideEnd = 0; TmpSideNext = 0
countDataBC=0
DO iBC=1,nPartBound
  IF (BCdata_auxSF(iBC)%SideNumber.EQ. -1) CYCLE !not set for SFs
  countDataBC=countDataBC+1
  TmpMapToBC(countDataBC)=iBC
END DO
DO BCSideID=1,nBCSides
  currentBC=0
  DO iBC=1,countDataBC
    IF (PartBound%MapToPartBC(BC(BCSideID)) .EQ. TmpMapToBC(iBC)) currentBC=iBC
  END DO
  IF (currentBC.EQ.0) CYCLE
  IF (TmpSideNumber(currentBC).EQ.0) THEN
    TmpSideStart(currentBC) = BCSideID ! Start of Linked List for Sides
  ELSE
    TmpSideNext(TmpSideEnd(currentBC)) = BCSideID ! Next Side
  END IF
  !-- prepare for next entry in list
  TmpSideEnd(currentBC) = BCSideID
  TmpSideNumber(currentBC) = TmpSideNumber(currentBC) + 1  ! Number of Sides
END DO ! BCSideID

IF(RadialWeighting%DoRadialWeighting) THEN
  ALLOCATE(BCdata_auxSFTemp(1:nPartBound))
END IF

!--- 2b: save sequential lists in BCdata_auxSF
DO iBC=1,countDataBC
  BCdata_auxSF(TmpMapToBC(iBC))%SideNumber=TmpSideNumber(iBC)
  IF (TmpSideNumber(iBC).EQ.0) CYCLE
  ALLOCATE(BCdata_auxSF(TmpMapToBC(iBC))%SideList(1:TmpSideNumber(iBC)))
  IF (TriaSurfaceFlux) THEN
    ALLOCATE(BCdata_auxSF(TmpMapToBC(iBC))%TriaSwapGeo(SurfFluxSideSize(1),SurfFluxSideSize(2),1:TmpSideNumber(iBC)))
    ALLOCATE(BCdata_auxSF(TmpMapToBC(iBC))%TriaSideGeo(1:TmpSideNumber(iBC)))
  END IF
  IF(RadialWeighting%DoRadialWeighting) THEN
    ALLOCATE(BCdata_auxSFTemp(TmpMapToBC(iBC))%WeightingFactor(1:TmpSideNumber(iBC)))
    BCdata_auxSFTemp(TmpMapToBC(iBC))%WeightingFactor = 1.
    ALLOCATE(BCdata_auxSFTemp(TmpMapToBC(iBC))%SubSideWeight(1:TmpSideNumber(iBC),1:RadialWeighting%nSubSides))
    BCdata_auxSFTemp(TmpMapToBC(iBC))%SubSideWeight = 0.
    ALLOCATE(BCdata_auxSFTemp(TmpMapToBC(iBC))%SubSideArea(1:TmpSideNumber(iBC),1:RadialWeighting%nSubSides))
    BCdata_auxSFTemp(TmpMapToBC(iBC))%SubSideArea = 0.
  END IF
  DO iSpec=1,nSpecies
    DO iSF=1,Species(iSpec)%nSurfacefluxBCs
      IF (TmpMapToBC(iBC).EQ.Species(iSpec)%Surfaceflux(iSF)%BC) THEN !only surfacefluxes with iBC
        ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(SurfFluxSideSize(1),SurfFluxSideSize(2),1:TmpSideNumber(iBC)))
        IF (UseCircularInflow .AND. (iSF .LE. Species(iSpec)%nSurfacefluxBCs)) THEN
          ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(1:TmpSideNumber(iBC)) )
        END IF
        IF(RadialWeighting%DoRadialWeighting) THEN
          ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%nVFRSub(1:TmpSideNumber(iBC),1:RadialWeighting%nSubSides))
        END IF
      END IF
    END DO
  END DO
  BCSideID=TmpSideStart(iBC)
  iCount=0
  DO !follow BCSideID list seq. with iCount
    iCount=iCount+1
    BCdata_auxSF(TmpMapToBC(iBC))%SideList(iCount)=BCSideID
    IF (TriaSurfaceFlux) THEN
      ElemID = SideToElem(S2E_ELEM_ID,BCSideID)
      iLocSide = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
      SideID=GetGlobalNonUniqueSideID(offsetElem+ElemID,iLocSide)
      !----- symmetry specific area calculation start
      IF(Symmetry%Order.EQ.2) THEN
        GlobalElemID = SideInfo_Shared(SIDE_ELEMID,SideID)
        CNElemID     = GetCNElemID(GlobalElemID)
        iLocSide = SideInfo_Shared(SIDE_LOCALID,SideID)
        IF(Symmetry%Axisymmetric) THEN
          ! Calculate the correct area for the axisymmetric (ring area) and 2D (length) and get ymin and ymax for element
          SurfMeshSubSideData(1,1,BCSideID)%area = DSMC_2D_CalcSymmetryArea(iLocSide,CNElemID, ymin, ymax)
          SurfMeshSubSideData(1,2,BCSideID)%area = 0.0
          ! Determination of the mean radial weighting factor for calculation of the number of particles to be inserted
          IF (RadialWeighting%DoRadialWeighting) THEN
            IF((ymax - ymin).GT.0.0) THEN
              ! Surfaces that are NOT parallel to the YZ-plane
              IF(RadialWeighting%CellLocalWeighting) THEN
                ! Cell local weighting
                BCdata_auxSFTemp(TmpMapToBC(iBC))%WeightingFactor(iCount) = (1. + ElemMidPoint_Shared(2,CNElemID) &
                                                                        / GEO%ymaxglob*(RadialWeighting%PartScaleFactor-1.))
              ELSE
                BCdata_auxSFTemp(TmpMapToBC(iBC))%WeightingFactor(iCount) = 1.
                DO iSub = 1, RadialWeighting%nSubSides
                  yMinTemp = ymin + (iSub-1) * (ymax - ymin) / RadialWeighting%nSubSides
                  yMaxTemp = ymin + iSub * (ymax - ymin) / RadialWeighting%nSubSides
                  BCdata_auxSFTemp(TmpMapToBC(iBC))%SubSideWeight(iCount,iSub) = 1.          &
                    + (yMaxTemp**2/(GEO%ymaxglob*2.)*(RadialWeighting%PartScaleFactor-1.) &
                    -  yMinTemp**2/(GEO%ymaxglob*2.)*(RadialWeighting%PartScaleFactor-1.))/(yMaxTemp - yMinTemp)
                END DO
                BCdata_auxSFTemp(TmpMapToBC(iBC))%SubSideArea(iCount,:) = DSMC_2D_CalcSymmetryAreaSubSides(iLocSide,CNElemID)
              END IF
            ELSE ! surfaces parallel to the x-axis (ymax = ymin)
              BCdata_auxSFTemp(TmpMapToBC(iBC))%WeightingFactor(iCount) = 1. &
                                                                        + ymax/(GEO%ymaxglob)*(RadialWeighting%PartScaleFactor-1.)
            END IF
          END IF
        ELSE
          SurfMeshSubSideData(1,1:2,BCSideID)%area = DSMC_2D_CalcSymmetryArea(iLocSide,CNElemID) / 2.
        END IF
      ELSE IF(Symmetry%Order.EQ.1) THEN
        SurfMeshSubSideData(1,1:2,BCSideID)%area = DSMC_1D_CalcSymmetryArea(iLocSide,ElemID) / 2.
      END IF
      !----- symmetry specific area calculation end
      IF (TrackingMethod.NE.TRIATRACKING) THEN !check that all sides are planar if TriaSurfaceFlux is used for tracing or refmapping
        CNSideID = GetCNSideID(SideID)
        IF (SideType(CNSideID).NE.PLANAR_RECT .AND. SideType(CNSideID).NE.PLANAR_NONRECT) CALL abort(__STAMP__,&
          'every surfaceflux-sides must be planar if TriaSurfaceFlux is used for tracing or refmapping!!!')
      END IF ! TrackingMethod.NE.TRIATRACKING

      DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
        CALL CalcNormAndTangTriangle(SideID=SideID &
          ,midpoint=BCdata_auxSF(TmpMapToBC(iBC))%TriaSwapGeo(iSample,jSample,iCount)%midpoint &
          ,ndist=BCdata_auxSF(TmpMapToBC(iBC))%TriaSwapGeo(iSample,jSample,iCount)%ndist &
          ,xyzNod=BCdata_auxSF(TmpMapToBC(iBC))%TriaSideGeo(iCount)%xyzNod &
          ,Vectors=BCdata_auxSF(TmpMapToBC(iBC))%TriaSideGeo(iCount)%Vectors &
          ,TriNum=jSample)
      END DO; END DO
    END IF !TriaSurfaceFlux

    !-- BC-list specific data
    DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
      BCdata_auxSF(TmpMapToBC(iBC))%LocalArea = BCdata_auxSF(TmpMapToBC(iBC))%LocalArea &
        + SurfMeshSubSideData(iSample,jSample,BCSideID)%area
    END DO; END DO

    !-- next Side
    IF (BCSideID .EQ. TmpSideEnd(iBC)) THEN
      IF (TmpSideNumber(iBC).NE.iCount) THEN
        CALL abort(__STAMP__,'Someting is wrong with TmpSideNumber of iBC',iBC,999.)
      ELSE
#ifdef CODE_ANALYZE
        IPWRITE(*,'(I4,I7,A53,I0)') iCount,' Sides have been found for Surfaceflux-linked PartBC ',TmpMapToBC(iBC)
#endif /*CODE_ANALYZE*/
        DoSurfaceFlux=.TRUE.
        EXIT
      END IF
    END IF
    BCSideID=TmpSideNext(BCSideID)
  END DO ! BCSideID (iCount)
END DO !iBC
!-- communicate areas
#if USE_MPI
   ALLOCATE( areasLoc(1:nPartBound) , areasGlob(1:nPartBound) )
   areasLoc=0.
   areasGlob=0.
   DO iPartBound=1,nPartBound
     areasLoc(iPartBound)=BCdata_auxSF(iPartBound)%LocalArea
   END DO
   CALL MPI_ALLREDUCE(areasLoc,areasGlob,nPartBound,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_PICLAS,IERROR)
#endif
   DO iPartBound=1,nPartBound
#if USE_MPI
     BCdata_auxSF(iPartBound)%GlobalArea=areasGlob(iPartBound)
#else
     BCdata_auxSF(iPartBound)%GlobalArea=BCdata_auxSF(iPartBound)%LocalArea
#endif
   END DO
#if USE_MPI
   DEALLOCATE(areasLoc,areasGlob)
#endif

END SUBROUTINE CreateSideListAndFinalizeAreasSurfFlux


SUBROUTINE DefineCircInflowRejectType(iSpec, iSF, iSide)
!===================================================================================================================================
!> Check where the sides are located relative to rmax (based on corner nodes of bounding box)
!> RejectType=0: complete side is inside valid bounds
!> RejectType=1: complete side is outside of valid bounds
!> RejectType=2: side is partly inside valid bounds
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars              ,ONLY: offsetElem, SideToElem
USE MOD_Particle_Surfaces_Vars ,ONLY: BCdata_auxSF
USE MOD_Particle_Vars          ,ONLY: Species
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Boundary_Tools,ONLY: GetRadialDistance2D
#ifdef CODE_ANALYZE
USE MOD_Particle_Vars          ,ONLY: CountCircInflowType
#endif
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)   :: iSpec, iSF, iSide
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: rmax, rmin
INTEGER               :: currentBC, BCSideID, ElemID, iLocSide, GlobalSideID
!===================================================================================================================================
currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
BCSideID=BCdata_auxSF(currentBC)%SideList(iSide)
ElemID = SideToElem(S2E_ELEM_ID,BCSideID)
iLocSide = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
GlobalSideID=GetGlobalNonUniqueSideID(offsetElem+ElemID,iLocSide)

CALL GetRadialDistance2D(GlobalSideID,Species(iSpec)%Surfaceflux(iSF)%dir,Species(iSpec)%Surfaceflux(iSF)%origin,rmin,rmax)

! define rejecttype
IF ( (rmin .GT. Species(iSpec)%Surfaceflux(iSF)%rmax) .OR. (rmax .LT. Species(iSpec)%Surfaceflux(iSF)%rmin) ) THEN
  Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide)=1
#ifdef CODE_ANALYZE
  CountCircInflowType(2,iSF,iSpec)=CountCircInflowType(2,iSF,iSpec)+1
#endif
ELSE IF ( (rmax .LE. Species(iSpec)%Surfaceflux(iSF)%rmax) .AND. (rmin .GE. Species(iSpec)%Surfaceflux(iSF)%rmin) ) THEN
  Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide)=0
#ifdef CODE_ANALYZE
  CountCircInflowType(1,iSF,iSpec)=CountCircInflowType(1,iSF,iSpec)+1
#endif
ELSE
  Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide)=2
#ifdef CODE_ANALYZE
  CountCircInflowType(3,iSF,iSpec)=CountCircInflowType(3,iSF,iSpec)+1
#endif
END IF !  (rmin > Surfaceflux-rmax) .OR. (rmax < Surfaceflux-rmin)
END SUBROUTINE DefineCircInflowRejectType


SUBROUTINE InitVolumeFlowRate(iSpec, iSF, iSide, tmp_SubSideAreas, BCdata_auxSFTemp)
!===================================================================================================================================
!> Calculate the volume flow rate (VFR) and store projected & tangential velocity vectors
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: BoltzmannConst, PI
USE MOD_Particle_Surfaces_Vars  ,ONLY: SurfFluxSideSize, SurfMeshSubSideData, tBCdata_auxSFRadWeight, BCdata_auxSF
USE MOD_Particle_Vars           ,ONLY: Species
USE MOD_DSMC_Vars               ,ONLY: RadialWeighting
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)   :: iSpec, iSF, iSide
REAL, INTENT(IN)      :: tmp_SubSideAreas(SurfFluxSideSize(1),SurfFluxSideSize(2))
TYPE(tBCdata_auxSFRadWeight), ALLOCATABLE, INTENT(IN)        :: BCdata_auxSFTemp(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: jSample, iSample, iSub, currentBC, BCSideID
REAL                  :: vec_nIn(3), nVFR, vec_t1(3), vec_t2(3), projFak, v_thermal, a, vSF
!===================================================================================================================================
currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
BCSideID=BCdata_auxSF(currentBC)%SideList(iSide)
DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
  vec_nIn = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn
  vec_t1 = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t1
  vec_t2 = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t2
  IF (.NOT.Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
    projFak = DOT_PRODUCT(vec_nIn,Species(iSpec)%Surfaceflux(iSF)%VeloVecIC) ! VeloVecIC projected to inwards normal
  ELSE
    projFak = 1.
  END IF
  v_thermal = SQRT(2.*BoltzmannConst*Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC/Species(iSpec)%MassIC) ! thermal speed
  a = 0. !dummy for projected speed ratio in constant v-distri
  !-- compute total volume flow rate through surface
  SELECT CASE(TRIM(Species(iSpec)%Surfaceflux(iSF)%velocityDistribution))
  CASE('constant')
    vSF = MAX(Species(iSpec)%Surfaceflux(iSF)%VeloIC * projFak,0.) ! VFR proj. to inwards normal (only positive parts!)
  CASE('maxwell','maxwell_lpn')
    IF ( ALMOSTEQUAL(v_thermal,0.)) THEN
      CALL abort(__STAMP__,' ERROR in SurfaceFlux: Calculated thermal velocity is zero! Temperature input might be missing (-MWTemperatureIC) ')
    END IF
    IF(Species(iSpec)%Surfaceflux(iSF)%VeloIC.GT.0.) THEN
      a = Species(iSpec)%Surfaceflux(iSF)%VeloIC * projFak / v_thermal !speed ratio proj. to inwards n (can be negative!)
      vSF = v_thermal / (2.0*SQRT(PI)) * ( EXP(-(a*a)) + a*SQRT(PI)*(1+ERF(a)) ) ! mean flux velocity through normal sub-face
    ELSE
      vSF = v_thermal / (2.0*SQRT(PI))  ! mean flux velocity through normal sub-face
    END IF
  CASE DEFAULT
    CALL abort(__STAMP__, 'ERROR in SurfaceFlux: Wrong velocity distribution!')
  END SELECT
  ! Calculate the volume flow rate. In case of an emission current and mass flow, it contains only the area
  IF(Species(iSpec)%Surfaceflux(iSF)%UseEmissionCurrent .OR. Species(iSpec)%Surfaceflux(iSF)%UseMassflow .OR. &
      Species(iSpec)%Surfaceflux(iSF)%Adaptive) THEN
    nVFR = tmp_SubSideAreas(iSample,jSample) ! Area
    ! vSF set to 1 to allow the utilization with radial weighting (untested)
    vSF = 1.
  ELSE
    nVFR = tmp_SubSideAreas(iSample,jSample) * vSF ! VFR projected to inwards normal of sub-side
  END IF
  IF(RadialWeighting%DoRadialWeighting) THEN
    nVFR = nVFR / BCdata_auxSFTemp(currentBC)%WeightingFactor(iSide)
    DO iSub = 1, RadialWeighting%nSubSides
      IF(ABS(BCdata_auxSFTemp(currentBC)%SubSideWeight(iSide,iSub)).GT.0.)THEN
        Species(iSpec)%Surfaceflux(iSF)%nVFRSub(iSide,iSub) = BCdata_auxSFTemp(currentBC)%SubSideArea(iSide,iSub) &
                                                            * vSF / BCdata_auxSFTemp(currentBC)%SubSideWeight(iSide,iSub)
      END IF
    END DO
  END IF
  IF (Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
    ! Check whether cell is completely outside of the circular inflow region and set the volume flow rate to zero
    IF (Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide).EQ.1) nVFR = 0.
  END IF
  Species(iSpec)%Surfaceflux(iSF)%VFR_total = Species(iSpec)%Surfaceflux(iSF)%VFR_total + nVFR
  !-- store SF-specific SubSide data in SurfFluxSubSideData (incl. projected velos)
  Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%nVFR = nVFR
  Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%projFak = projFak
  Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%a_nIn = a
  IF (.NOT.Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
    Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t1 &
      = Species(iSpec)%Surfaceflux(iSF)%VeloIC &
      * DOT_PRODUCT(vec_t1,Species(iSpec)%Surfaceflux(iSF)%VeloVecIC) !v in t1-dir
    Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t2 &
      = Species(iSpec)%Surfaceflux(iSF)%VeloIC &
      * DOT_PRODUCT(vec_t2,Species(iSpec)%Surfaceflux(iSF)%VeloVecIC) !v in t2-dir
  ELSE
    Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t1 = 0. !v in t1-dir
    Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t2 = 0. !v in t2-dir
  END IF! .NOT.VeloIsNormal
END DO; END DO !jSample=1,SurfFluxSideSize(2); iSample=1,SurfFluxSideSize(1)

END SUBROUTINE InitVolumeFlowRate


SUBROUTINE InitReduceNoiseSF(iSpec, iSF)
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: Species
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)   :: iSpec, iSF
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
INTEGER                :: iProc
#endif  /*USE_MPI*/
!===================================================================================================================================
IF(MPIroot)THEN
  ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcs(0:nProcessors-1))
  Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcs=0.
ELSE
  ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcs(1)) !dummy for debug
END IF !MPIroot
#if USE_MPI
CALL MPI_GATHER(Species(iSpec)%Surfaceflux(iSF)%VFR_total,1,MPI_DOUBLE_PRECISION &
  ,Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcs,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_PICLAS,iError)
IF(MPIroot)THEN
  DO iProc=0,nProcessors-1
    Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcsTotal = Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcsTotal &
      + Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcs(iProc)
  END DO
END IF
#else  /*USE_MPI*/
Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcs=Species(iSpec)%Surfaceflux(iSF)%VFR_total
Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcsTotal=Species(iSpec)%Surfaceflux(iSF)%VFR_total
#endif  /*USE_MPI*/
END SUBROUTINE InitReduceNoiseSF


SUBROUTINE CalcCircInflowAreaPerTria(iSpec,iSF)
!===================================================================================================================================
!> Routine calculates the partial area of the circular infow per triangle for AdaptiveType=4 (Monte Carlo integration)
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Vars               ,ONLY:Species
USE MOD_Particle_Surfaces_Vars      ,ONLY:SurfMeshSubSideData, BCdata_auxSF, SurfFluxSideSize, TriaSurfaceFlux
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: iSpec, iSF
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSide, BCSideID, iMC, iSample, jSample, dir(3), currentBC, MCVar, Node1, Node2, counter
REAL                            :: CircleArea, point(2), origin(2), radius, Vector1(3), Vector2(3), RandVal2(2), Particle_pos(3)
REAL                            :: PartDistance, TriaNode(1:3), midpoint(1:3), ndist(1:3)
!===================================================================================================================================

MCVar = 1000000

IF (.NOT.TriaSurfaceFlux) THEN
  CALL abort(&
__STAMP__&
,'ERROR: CircularInflow only with TriaSurfaceFlux!')
END IF

currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
dir       = Species(iSpec)%Surfaceflux(iSF)%dir
origin    = Species(iSpec)%Surfaceflux(iSF)%origin

DO iSide=1,BCdata_auxSF(currentBC)%SideNumber
  BCSideID=BCdata_auxSF(currentBC)%SideList(iSide)
  IF(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide).EQ.1) CYCLE
  TriaNode(1:3) = BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%xyzNod(1:3)
  DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
    !-- compute parallelogram of triangle
    Node1 = jSample+1     ! normal = cross product of 1-2 and 1-3 for first triangle
    Node2 = jSample+2     !          and 1-3 and 1-4 for second triangle
    Vector1 = BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%Vectors(:,Node1-1)
    Vector2 = BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%Vectors(:,Node2-1)
    midpoint(1:3) = BCdata_auxSF(currentBC)%TriaSwapGeo(iSample,jSample,iSide)%midpoint(1:3)
    ndist(1:3) = BCdata_auxSF(currentBC)%TriaSwapGeo(iSample,jSample,iSide)%ndist(1:3)
    IF(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide).EQ.0) THEN
      CircleArea = SurfMeshSubSideData(iSample,jSample,BCSideID)%area
    ELSE
      CircleArea = 0.
      counter = 0
      DO iMC = 1,MCVar
        CALL RANDOM_NUMBER(RandVal2)
        Particle_pos(1:3) = TriaNode(1:3) + Vector1(1:3) * RandVal2(1)
        Particle_pos(1:3) = Particle_pos(1:3) + Vector2(1:3) * RandVal2(2)
        PartDistance = ndist(1)*(Particle_pos(1)-midpoint(1)) & !Distance from v1-v2
                      + ndist(2)*(Particle_pos(2)-midpoint(2)) &
                      + ndist(3)*(Particle_pos(3)-midpoint(3))
        IF (PartDistance.GT.0.) THEN !flip into right triangle if outside
          Particle_pos(1:3) = 2*midpoint(1:3)-Particle_pos(1:3)
        END IF
        point(1)=Particle_pos(dir(2))-origin(1)
        point(2)=Particle_pos(dir(3))-origin(2)
        radius=SQRT( (point(1))**2+(point(2))**2 )
        IF ((radius.LE.Species(iSpec)%Surfaceflux(iSF)%rmax).AND.(radius.GE.Species(iSpec)%Surfaceflux(iSF)%rmin)) THEN
          CircleArea = CircleArea + 1./REAL(MCVar)
          counter = counter + 1
        END IF
      END DO
      CircleArea = CircleArea * SurfMeshSubSideData(iSample,jSample,BCSideID)%area
    END IF
    Species(iSpec)%Surfaceflux(iSF)%CircleAreaPerTriaSide(iSample,jSample,iSide) = CircleArea
  END DO; END DO
END DO

END SUBROUTINE CalcCircInflowAreaPerTria


SUBROUTINE CalcConstMassflowWeightForZeroMassFlow(iSpec,iSF)
!===================================================================================================================================
!> Routine calculates the ConstMassflowWeight for AdaptiveType=4 in case of zero massflow rate
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Vars           ,ONLY: Species
USE MOD_Particle_Surfaces_Vars  ,ONLY: SurfMeshSubSideData, BCdata_auxSF, SurfFluxSideSize
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: iSpec, iSF
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSample, jSample, currentBC, iSide, BCSideID
!===================================================================================================================================

currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
DO iSide=1,BCdata_auxSF(currentBC)%SideNumber
  BCSideID=BCdata_auxSF(currentBC)%SideList(iSide)
  DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
    Species(iSpec)%Surfaceflux(iSF)%ConstMassflowWeight(iSample,jSample,iSide)=SurfMeshSubSideData(iSample,jSample,BCSideID)%area &
                                                                                / Species(iSpec)%Surfaceflux(iSF)%totalAreaSF
  END DO; END DO
END DO

END SUBROUTINE CalcConstMassflowWeightForZeroMassFlow


SUBROUTINE ReadInThermionicEmission(iSpec,iSF)
!===================================================================================================================================
!> Thermionic emission model: Read-in of work function and Richardson constant and calculation of the current density
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst, Pi, eV2Joule
USE MOD_Particle_Vars          ,ONLY: Species
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER               :: iSpec, iSF
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(42)         :: help, help2
REAL                  :: WallTemp
!===================================================================================================================================

ASSOCIATE(SF => Species(iSpec)%Surfaceflux(iSF))

! Sanity checks
IF(SF%UseEmissionCurrent) THEN
  CALL abort(__STAMP__,'ERROR in Surface Flux: Thermionic emission modelling and emission current cannot be used at the same time!')
ELSE
  ! Utilize the same calculation of particles per iteration
  SF%UseEmissionCurrent = .TRUE.
END IF
IF(PartBound%TargetBoundCond(SF%BC) .NE. PartBound%ReflectiveBC) THEN
  CALL abort(__STAMP__,'ERROR in Surface Flux: Thermionic emission modelling requires a reflective boundary condition!')
END IF

! Read-in of variables
WRITE(UNIT=help,FMT='(I0)') iSpec
WRITE(UNIT=help2,FMT='(I0)') iSF
help2=TRIM(help)//'-Surfaceflux'//TRIM(help2)

! Consider influence of electric field on the material work function (Schottky effect)
SF%SchottkyEffectTE = GETLOGICAL('Part-Species'//TRIM(help2)//'-ThermionicEmission-SchottkyEffect')
IF(SF%SchottkyEffectTE) THEN
  SF%Type = 5
#if !(USE_HDG)
  CALL abort(__STAMP__,'ERROR in Surface Flux: Thermionic emission with Schottky effect requires an electric field!')
#endif
END IF
! Material-specific work function read-in in eV and converted to K
SF%WorkFunctionTE = GETREAL('Part-Species'//TRIM(help2)//'-ThermionicEmission-WorkFunction') * eV2Joule
! Material-specific constant read-in in A/(cm^2 K^2) and converted to m^2
SF%RichardsonConstant = GETREAL('Part-Species'//TRIM(help2)//'-ThermionicEmission-RichardsonConstant') * 1E4

! Richardson-Dushman equation to calculate the curren density [A/m2]
! (stored in the EmissionCurrent variable, different treatment in InitializeParticleSurfaceflux)
WallTemp = PartBound%WallTemp(SF%BC)
IF(WallTemp.GT.0) THEN
  SF%EmissionCurrent = SF%RichardsonConstant * WallTemp**2 * EXP(-SF%WorkFunctionTE / (BoltzmannConst * WallTemp))
ELSE
  CALL abort(__STAMP__,'ERROR in Surface Flux: Thermionic emission modelling requires a wall temperature for BC: ',IntInfoOpt=SF%BC)
END IF

END ASSOCIATE

END SUBROUTINE ReadInThermionicEmission

END MODULE MOD_Particle_SurfFlux_Init
