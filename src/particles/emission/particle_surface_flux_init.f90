!==================================================================================================================================
! Copyright (c) 2021 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-VeloIC' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Velocity for inital Data', '0.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Species[$]-Surfaceflux[$]-VeloIsNormal' &
                                , 'TODO-DEFINE-PARAMETER VeloIC is in Surf-Normal instead of VeloVecIC' &
                                , '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-Surfaceflux[$]-VeloVecIC' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Normalized velocity vector' , '0.0 , 0.0 , 0.0', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Species[$]-Surfaceflux[$]-CircularInflow' &
                                , 'Enables the utilization of a circular region as a surface flux on the selected boundary. '//&
                                  'Only possible on surfaces, which are in xy, xz, and yz-planes.' &
                                , '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-Surfaceflux[$]-axialDir' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Axial direction of coordinates in polar system', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-Surfaceflux[$]-origin' &
                                , 'TODO-DEFINE-PARAMETER Origin in orth(ogonal?) coordinates of polar system' , '0.0 , 0.0'&
                                ,  numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-rmax' &
                                , 'TODO-DEFINE-PARAMETER Max radius of to-be inserted particles', '1e21', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-rmin' &
                                , 'TODO-DEFINE-PARAMETER Min radius of to-be inserted particles', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-MWTemperatureIC' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Temperature for Maxwell Distribution', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-PartDensity' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'PartDensity (real particles per m^3) or  (vpi_)cub./cyl. as alternative  to'//&
                                  ' Part.Emis. in Type1'  , '0.', numberedmulti=.TRUE.)
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
! Parameters for adaptive boundary conditions
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
USE MOD_Particle_Sampling_Vars ,ONLY: UseAdaptive
USE MOD_SurfaceModel_Chemistry
USE MOD_SurfaceModel_Vars
!USE MOD_Particle_SurfChemFlux_Init
USE MOD_Restart_Vars           ,ONLY: DoRestart, RestartTime
#if USE_MPI
USE MOD_Particle_Vars          ,ONLY: DoPoissonRounding, DoTimeDepInflow
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*USE_MPI*/
#ifdef CODE_ANALYZE
USE MOD_Particle_Vars          ,ONLY: CountCircInflowType
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Local variable declaration
INTEGER               :: iReac, SurfNumReac
INTEGER               :: iSpec,iSF,SideID,BCSideID,iSide,ElemID,iLocSide,iSample,jSample,currentBC
INTEGER               :: iCopy1, iCopy2, iCopy3, MaxSurfacefluxBCs,nDataBC
REAL                  :: tmp_SubSideDmax(SurfFluxSideSize(1),SurfFluxSideSize(2))
REAL                  :: tmp_SubSideAreas(SurfFluxSideSize(1),SurfFluxSideSize(2))
REAL                  :: tmp_BezierControlPoints2D(2,0:NGeo,0:NGeo,SurfFluxSideSize(1),SurfFluxSideSize(2))
REAL                  :: VFR_total
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
UseAdaptive=.FALSE.
MaxSurfacefluxBCs=0
nDataBC=0
DoSurfaceFlux=.FALSE.

!-- 1.: read/prepare parameters and determine nec. BCs
CALL ReadInAndPrepareSurfaceFlux(MaxSurfacefluxBCs, nDataBC)

! Call to the reactive surfaces
CALL ReadInAndPrepareSurfChemFlux(nDataBC)

#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoPoissonRounding,1,MPI_LOGICAL,MPI_LAND,PartMPI%COMM,iError) !set T if this is for all procs
CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoTimeDepInflow,1,MPI_LOGICAL,MPI_LAND,PartMPI%COMM,iError) !set T if this is for all procs
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
    IF (BCdata_auxSF(currentBC)%SideNumber.GT.0) THEN
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
            IF(Species(iSpec)%Surfaceflux(iSF)%Adaptive) THEN
              IF(.NOT.Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
                Species(iSpec)%Surfaceflux(iSF)%totalAreaSF = Species(iSpec)%Surfaceflux(iSF)%totalAreaSF &
                                                              + SurfMeshSubSideData(iSample,jSample,BCSideID)%area
              END IF
            END IF
          END DO; END DO
        END IF
        ! Initialize circular inflow (determine if elements are (partially) inside/outside)
        IF (Species(iSpec)%Surfaceflux(iSF)%CircularInflow) CALL DefineCircInflowRejectType(iSpec, iSF, iSide)
        ! Initialize surface flux
        CALL InitSurfFlux(iSpec, iSF, iSide, tmp_SubSideAreas, BCdata_auxSFTemp)
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
        !Init adaptive SF
      END DO ! iSide
    ELSE IF (BCdata_auxSF(currentBC)%SideNumber.EQ.-1) THEN
      CALL abort(__STAMP__&
        ,'ERROR in ParticleSurfaceflux: Someting is wrong with SideNumber of BC ',currentBC)
    END IF
    !--- 3b: ReduceNoise initialization
    IF (Species(iSpec)%Surfaceflux(iSF)%ReduceNoise) CALL InitReduceNoiseSF(iSpec, iSF)
    !--- Finalize adaptive SF
#if USE_MPI
    IF(Species(iSpec)%Surfaceflux(iSF)%Adaptive) THEN
      IF(.NOT.Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
        totalAreaSF_global = 0.0
        CALL MPI_ALLREDUCE(Species(iSpec)%Surfaceflux(iSF)%totalAreaSF,totalAreaSF_global,1, &
                            MPI_DOUBLE_PRECISION,MPI_SUM,PartMPI%COMM,IERROR)
        Species(iSpec)%Surfaceflux(iSF)%totalAreaSF = totalAreaSF_global
      END IF
    END IF
#endif

#ifdef CODE_ANALYZE
    IF (BCdata_auxSF(currentBC)%SideNumber.GT.0 .AND. Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
      IPWRITE(*,'(I4,A,2(1X,I0),A,3(1X,I0))') ' For Surfaceflux/Spec',iSF,iSpec,' are nType0,1,2: ' &
                                            , CountCircInflowType(1,iSF,iSpec),CountCircInflowType(2, iSF,iSpec) &
                                            , CountCircInflowType(3, iSF,iSpec)
    END IF
#endif /*CODE_ANALYZE*/
  END DO !iSF
END DO !iSpec

SurfNumReac = SurfChemReac%NumOfReact
!initialize Surfaceflux-specific data
DO iReac=1,SurfNumReac
  DO iSF=1,SurfChemReac%NumOfBounds(iReac)
    currentBC = SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%BC
    IF (BCdata_auxSF(currentBC)%SideNumber.GT.0) THEN
      
      ! Loop over sides on the surface flux
      DO iSide=1,BCdata_auxSF(currentBC)%SideNumber
        BCSideID=BCdata_auxSF(currentBC)%SideList(iSide)
        ElemID = SideToElem(S2E_ELEM_ID,BCSideID)
        iLocSide = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
        SideID=GetGlobalNonUniqueSideID(offsetElem+ElemID,iLocSide)
        IF (SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%AcceptReject) THEN
          CALL GetBezierSampledAreas(SideID=SideID,BezierSampleN=BezierSampleN &
            ,BezierSurfFluxProjection_opt=.NOT.SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloIsNormal  &
            ,SurfMeshSubSideAreas=tmp_SubSideAreas,DmaxSampleN_opt=SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%ARM_DmaxSampleN &
            ,Dmax_opt=tmp_SubSideDmax,BezierControlPoints2D_opt=tmp_BezierControlPoints2D)
        ELSE IF (.NOT.TriaSurfaceFlux) THEN
          CALL GetBezierSampledAreas(SideID=SideID,BezierSampleN=BezierSampleN &
            ,BezierSurfFluxProjection_opt=.NOT.SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloIsNormal,SurfMeshSubSideAreas=tmp_SubSideAreas)
        ELSE ! TriaSurfaceFlux
          DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
            tmp_SubSideAreas(iSample,jSample)=SurfMeshSubSideData(iSample,jSample,BCSideID)%area
          END DO; END DO
        END IF
        ! Initialize surface flux
        CALL InitSurfChemFlux(iReac, iSF, iSide, tmp_SubSideAreas, BCdata_auxSFTemp)
        ! Initialize acceptance-rejection on SF
        IF (SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%AcceptReject) THEN
          DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
            SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Dmax = tmp_SubSideDmax(iSample,jSample)
          !  IF (.NOT.SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloIsNormal) THEN
         !     ALLOCATE(SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample &
           !                                                               ,iSide)%BezierControlPoints2D(1:2,0:NGeo,0:NGeo))
            !  DO iCopy1=0,NGeo; DO iCopy2=0,NGeo; DO iCopy3=1,2
             !   SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample &
              !                                                     ,iSide)%BezierControlPoints2D(iCopy3,iCopy2,iCopy1) &
               !   = tmp_BezierControlPoints2D(iCopy3,iCopy2,iCopy1,iSample,jSample)
             ! END DO; END DO; END DO
            !END IF !.NOT.VeloIsNormal
          END DO; END DO !jSample=1,SurfFluxSideSize(2); iSample=1,SurfFluxSideSize(1)
        END IF
        !Init adaptive SF
      END DO ! iSide
    ELSE IF (BCdata_auxSF(currentBC)%SideNumber.EQ.-1) THEN
      CALL abort(__STAMP__&
        ,'ERROR in ParticleSurfaceflux: Someting is wrong with SideNumber of BC ',currentBC)
    END IF
    !--- Finalize adaptive SF
#if USE_MPI
    IF(SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%Adaptive) THEN
      totalAreaSF_global = 0.0
      CALL MPI_ALLREDUCE(SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%totalAreaSF,totalAreaSF_global,1, &
                          MPI_DOUBLE_PRECISION,MPI_SUM,PartMPI%COMM,IERROR)
      SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%totalAreaSF = totalAreaSF_global
    END IF
#endif

#ifdef CODE_ANALYZE
    IF (BCdata_auxSF(currentBC)%SideNumber.GT.0 .AND. SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%CircularInflow) THEN
      IPWRITE(*,'(I4,A,2(1X,I0),A,3(1X,I0))') ' For Surfaceflux/Spec',iSF,iSpec,' are nType0,1,2: ' &
                                            , CountCircInflowType(1,iSF,iSpec),CountCircInflowType(2, iSF,iSpec) &
                                            , CountCircInflowType(3, iSF,iSpec)
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
  DO iReac=1,SurfNumReac
    DO iSF=1,SurfChemReac%NumOfBounds(iReac)
      VFR_total = SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VFR_total               !proc local total
      !Species(iSpec)%Surfaceflux(iSF)%InsertedParticle = INT(Species(iSpec)%Surfaceflux(iSF)%PartDensity * RestartTime &
      ! / Species(iSpec)%MacroParticleFactor * VFR_total,8)
    END DO
  END DO
END IF

IF(DoRestart) THEN
  DO iSpec=1,nSpecies
    DO iSF = 1, Species(iSpec)%NumberOfInits
      Species(iSpec)%Init(iSF)%InsertedParticle = INT(Species(iSpec)%Init(iSF)%ParticleNumber * RestartTime,8)
    END DO
    DO iSF = 1, Species(iSpec)%nSurfacefluxBCs
      IF (Species(iSpec)%Surfaceflux(iSF)%ReduceNoise) THEN
        VFR_total = Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcsTotal !proc global total (for non-root: dummy!!!)
      ELSE
        VFR_total = Species(iSpec)%Surfaceflux(iSF)%VFR_total               !proc local total
      END IF
      Species(iSpec)%Surfaceflux(iSF)%InsertedParticle = INT(Species(iSpec)%Surfaceflux(iSF)%PartDensity * RestartTime &
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
        ! Weighting factor to account for a
        ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%ConstMassflowWeight(1:SurfFluxSideSize(1),1:SurfFluxSideSize(2), &
                  1:BCdata_auxSF(currentBC)%SideNumber))
        Species(iSpec)%Surfaceflux(iSF)%ConstMassflowWeight = 0.0
        ! In case the mass flow is set to zero, the weighting factor can be determined once initially based on area ratio
        IF(ALMOSTEQUAL(Species(iSpec)%Surfaceflux(iSF)%AdaptiveMassflow,0.)) CALL CalcConstMassflowWeightForZeroMassFlow(iSpec,iSF)
        ! Circular inflow in combination with AdaptiveType = 4 requires the partial circle area per tria side
        IF(Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
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
CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoSurfaceFlux,1,MPI_LOGICAL,MPI_LOR,PartMPI%COMM,iError) !set T if at least 1 proc have SFs
#endif  /*USE_MPI*/
IF ((.NOT.DoSurfaceFlux).AND. (.NOT.DoChemSurface)) THEN !-- no SFs defined
  SWRITE(*,*) 'WARNING: No Sides for SurfacefluxBCs found! DoSurfaceFlux is now disabled!'
END IF
DoForceFreeSurfaceFlux = GETLOGICAL('DoForceFreeSurfaceFlux','.FALSE.')

END SUBROUTINE InitializeParticleSurfaceflux


SUBROUTINE ReadInAndPrepareSurfaceFlux(MaxSurfacefluxBCs, nDataBC)
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst, Pi
USE MOD_Particle_Vars          ,ONLY: nSpecies, Species, VarTimeStep, DoPoissonRounding, DoTimeDepInflow
USE MOD_Particle_Vars          ,ONLY: Symmetry, UseCircularInflow
USE MOD_Particle_Sampling_Vars ,ONLY: UseAdaptive
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound,nPartBound
USE MOD_DSMC_Vars              ,ONLY: useDSMC, BGGas
USE MOD_Particle_Surfaces_Vars ,ONLY: BCdata_auxSF, BezierSampleN, TriaSurfaceFlux
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_Mesh_Vars              ,ONLY: NGeo
USE MOD_SurfaceModel_Chemistry
USE MOD_SurfaceModel_Vars
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
  Species(iSpec)%nSurfacefluxBCs = GETINT('Part-Species'//TRIM(hilf)//'-nSurfacefluxBCs')
  IF (useDSMC.AND.(BGGas%NumberOfSpecies.GT.0)) THEN
    IF (BGGas%BackgroundSpecies(iSpec)) THEN
      IF (Species(iSpec)%nSurfacefluxBCs.GT.0) CALL abort(__STAMP__,&
        'SurfaceFlux is not implemented for the BGG-species!')
    END IF
  END IF
  IF (Species(iSpec)%nSurfacefluxBCs.EQ.0) THEN
    CYCLE
  ELSE
    ALLOCATE(Species(iSpec)%Surfaceflux(1:Species(iSpec)%nSurfacefluxBCs))
    ! Initialize Surfaceflux to BC mapping
    Species(iSpec)%Surfaceflux(:)%BC=-1
    DO iSF=1,Species(iSpec)%nSurfacefluxBCs
      WRITE(UNIT=hilf2,FMT='(I0)') iSF
      hilf2=TRIM(hilf)//'-Surfaceflux'//TRIM(hilf2)
      Species(iSpec)%Surfaceflux(iSF)%BC = GETINT('Part-Species'//TRIM(hilf2)//'-BC')
    END DO
  END IF

  MaxSurfacefluxBCs=MAX(MaxSurfacefluxBCs,Species(iSpec)%nSurfacefluxBCs)
  DO iSF=1,Species(iSpec)%nSurfacefluxBCs
    WRITE(UNIT=hilf2,FMT='(I0)') iSF
    hilf2=TRIM(hilf)//'-Surfaceflux'//TRIM(hilf2)
    Species(iSpec)%Surfaceflux(iSF)%InsertedParticle = 0
    Species(iSpec)%Surfaceflux(iSF)%InsertedParticleSurplus = 0
    Species(iSpec)%Surfaceflux(iSF)%VFR_total = 0
    Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcsTotal = 0
    ! get surfaceflux data
    IF (Species(iSpec)%Surfaceflux(iSF)%BC.LT.1 .OR. Species(iSpec)%Surfaceflux(iSF)%BC.GT.nPartBound) THEN
      CALL abort(&
__STAMP__&
, 'SurfacefluxBCs must be between 1 and nPartBound!')
    ELSE IF (BCdata_auxSF(Species(iSpec)%Surfaceflux(iSF)%BC)%SideNumber.EQ. -1) THEN !not set yet
      BCdata_auxSF(Species(iSpec)%Surfaceflux(iSF)%BC)%SideNumber=0
      nDataBC=nDataBC+1
    END IF
    Species(iSpec)%Surfaceflux(iSF)%velocityDistribution  = &
        TRIM(GETSTR('Part-Species'//TRIM(hilf2)//'-velocityDistribution','constant'))
    IF (TRIM(Species(iSpec)%Surfaceflux(iSF)%velocityDistribution).NE.'constant' .AND. &
        TRIM(Species(iSpec)%Surfaceflux(iSF)%velocityDistribution).NE.'maxwell' .AND. &
        TRIM(Species(iSpec)%Surfaceflux(iSF)%velocityDistribution).NE.'maxwell_lpn') THEN
      CALL abort(&
__STAMP__&
, 'Only constant or maxwell-like velodistri implemented for surfaceflux!')
    END IF
    Species(iSpec)%Surfaceflux(iSF)%VeloIC                = GETREAL('Part-Species'//TRIM(hilf2)//'-VeloIC')
    Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal          = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-VeloIsNormal')
    IF (Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
      Species(iSpec)%Surfaceflux(iSF)%CircularInflow=.FALSE.
    ELSE
      Species(iSpec)%Surfaceflux(iSF)%VeloVecIC          =GETREALARRAY('Part-Species'//TRIM(hilf2)//'-VeloVecIC',3)
      Species(iSpec)%Surfaceflux(iSF)%CircularInflow=GETLOGICAL('Part-Species'//TRIM(hilf2)//'-CircularInflow')
      IF(Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
        UseCircularInflow=.TRUE.
        Species(iSpec)%Surfaceflux(iSF)%dir(1)       = GETINT('Part-Species'//TRIM(hilf2)//'-axialDir')
        IF (Species(iSpec)%Surfaceflux(iSF)%dir(1).EQ.1) THEN
          Species(iSpec)%Surfaceflux(iSF)%dir(2)=2
          Species(iSpec)%Surfaceflux(iSF)%dir(3)=3
        ELSE IF (Species(iSpec)%Surfaceflux(iSF)%dir(1).EQ.2) THEN
          Species(iSpec)%Surfaceflux(iSF)%dir(2)=3
          Species(iSpec)%Surfaceflux(iSF)%dir(3)=1
        ELSE IF (Species(iSpec)%Surfaceflux(iSF)%dir(1).EQ.3) THEN
          Species(iSpec)%Surfaceflux(iSF)%dir(2)=1
          Species(iSpec)%Surfaceflux(iSF)%dir(3)=2
        ELSE
          CALL abort(__STAMP__&
            ,'ERROR in init: axialDir for SFradial must be between 1 and 3!')
        END IF
        IF ( Species(iSpec)%Surfaceflux(iSF)%VeloVecIC(Species(iSpec)%Surfaceflux(iSF)%dir(2)).NE.0. .OR. &
             Species(iSpec)%Surfaceflux(iSF)%VeloVecIC(Species(iSpec)%Surfaceflux(iSF)%dir(3)).NE.0. ) THEN
          CALL abort(__STAMP__&
            ,'ERROR in init: axialDir for SFradial do not correspond to VeloVecIC!')
        END IF
        Species(iSpec)%Surfaceflux(iSF)%origin       = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-origin',2,'0. , 0.')
        WRITE(UNIT=hilf3,FMT='(E16.8)') HUGE(Species(iSpec)%Surfaceflux(iSF)%rmax)
        Species(iSpec)%Surfaceflux(iSF)%rmax     = GETREAL('Part-Species'//TRIM(hilf2)//'-rmax',TRIM(hilf3))
        Species(iSpec)%Surfaceflux(iSF)%rmin     = GETREAL('Part-Species'//TRIM(hilf2)//'-rmin','0.')
      END IF
    END IF !.NOT.VeloIsNormal
    IF (.NOT.Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
      !--- normalize VeloVecIC
      IF (.NOT. ALL(Species(iSpec)%Surfaceflux(iSF)%VeloVecIC(:).eq.0.)) THEN
        Species(iSpec)%Surfaceflux(iSF)%VeloVecIC = Species(iSpec)%Surfaceflux(iSF)%VeloVecIC &
          /SQRT(DOT_PRODUCT(Species(iSpec)%Surfaceflux(iSF)%VeloVecIC,Species(iSpec)%Surfaceflux(iSF)%VeloVecIC))
      END IF
    END IF
    Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC       = GETREAL('Part-Species'//TRIM(hilf2)//'-MWTemperatureIC','0.')
    Species(iSpec)%Surfaceflux(iSF)%PartDensity           = GETREAL('Part-Species'//TRIM(hilf2)//'-PartDensity','0.')
    Species(iSpec)%Surfaceflux(iSF)%ReduceNoise           = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-ReduceNoise','.FALSE.')
    IF (DoPoissonRounding .AND. Species(iSpec)%Surfaceflux(iSF)%ReduceNoise) THEN
      SWRITE(*,*)'WARNING: Poisson sampling not possible for noise reduction of surfacefluxes:'
      SWRITE(*,*)'switching now to Random rounding...'
      DoPoissonRounding   = .FALSE.
    END IF
    IF (DoTimeDepInflow .AND. Species(iSpec)%Surfaceflux(iSF)%ReduceNoise) THEN
      SWRITE(*,*)'WARNING: Time-dependent inflow is not possible for noise reduction of surfacefluxes:'
      SWRITE(*,*)'switching now to Random rounding...'
      DoTimeDepInflow   = .FALSE.
    END IF
    IF (TriaSurfaceFlux) THEN
      Species(iSpec)%Surfaceflux(iSF)%AcceptReject=.FALSE.
    ELSE
      Species(iSpec)%Surfaceflux(iSF)%AcceptReject          = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-AcceptReject','.TRUE.')
    END IF
    IF (Species(iSpec)%Surfaceflux(iSF)%AcceptReject .AND. BezierSampleN.GT.1) THEN
      SWRITE(*,*)'WARNING: BezierSampleN > 0 may not be necessary as ARM is used for SurfaceFlux!'
    ELSE IF (.NOT.Species(iSpec)%Surfaceflux(iSF)%AcceptReject .AND. BezierSampleN.LE.NGeo .AND. .NOT.TriaSurfaceFlux) THEN
      SWRITE(*,*)'WARNING: The choosen small BezierSampleN (def.: NGeo) might result in inhom. SurfFluxes without ARM!'
    END IF
    IF (Species(iSpec)%Surfaceflux(iSF)%AcceptReject) THEN
      WRITE( hilf3, '(I0.2)') NGeo*NGeo*NGeo !1 for linear elements, this is an arbitray estimation for higher N!
      Species(iSpec)%Surfaceflux(iSF)%ARM_DmaxSampleN = GETINT('Part-Species'//TRIM(hilf2)//'-ARM_DmaxSampleN',hilf3)
    ELSE
      Species(iSpec)%Surfaceflux(iSF)%ARM_DmaxSampleN = 0
    END IF
    ! ================================= ADAPTIVE BC READ IN START =================================================================!
    Species(iSpec)%Surfaceflux(iSF)%Adaptive         = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-Adaptive','.FALSE.')
    IF(Species(iSpec)%Surfaceflux(iSF)%Adaptive) THEN
      DoPoissonRounding = .TRUE.
      UseAdaptive  = .TRUE.
      IF(TrackingMethod.EQ.REFMAPPING) THEN
        CALL abort(__STAMP__&
            ,'ERROR: Adaptive surface flux boundary conditions are not implemented with RefMapping!')
      END IF
      IF((Symmetry%Order.LE.2).OR.VarTimeStep%UseVariableTimeStep) THEN
        CALL abort(__STAMP__&
            ,'ERROR: Adaptive surface flux boundary conditions are not implemented with 2D/axisymmetric or variable time step!')
      END IF
      ! Total area of surface flux
      IF(Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
        Species(iSpec)%Surfaceflux(iSF)%totalAreaSF = Pi*(Species(iSpec)%Surfaceflux(iSF)%rmax &
          *Species(iSpec)%Surfaceflux(iSF)%rmax - Species(iSpec)%Surfaceflux(iSF)%rmin*Species(iSpec)%Surfaceflux(iSF)%rmin)
      ELSE
        Species(iSpec)%Surfaceflux(iSF)%totalAreaSF = 0.
      END IF
      Species(iSpec)%Surfaceflux(iSF)%AdaptiveType         = GETINT('Part-Species'//TRIM(hilf2)//'-Adaptive-Type')
      SELECT CASE(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType)
      ! Farbar2014 - Case 1: Inlet Type 1, constant pressure and temperature
      !              Case 2: Outlet Type 1, constant pressure
      !              Case 3: Inlet Type 2, constant mass flow and temperature
      ! Lei2017    - Case 4: cf_3, constant mass flow and temperature N through mass flow and particles out
      CASE(1,2)
        Species(iSpec)%Surfaceflux(iSF)%AdaptivePressure  = GETREAL('Part-Species'//TRIM(hilf2)//'-Adaptive-Pressure')
        Species(iSpec)%Surfaceflux(iSF)%PartDensity       = Species(iSpec)%Surfaceflux(iSF)%AdaptivePressure &
                                                            / (BoltzmannConst * Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC)
      CASE(3,4)
        Species(iSpec)%Surfaceflux(iSF)%AdaptiveMassflow     = GETREAL('Part-Species'//TRIM(hilf2)//'-Adaptive-Massflow')
        IF(Species(iSpec)%Surfaceflux(iSF)%VeloIC.LE.0.0) THEN
          CALL abort(__STAMP__&
            ,'ERROR in init of adaptive inlet: positive initial guess of velocity for Type 3/Type 4 condition required!')
        END IF
      END SELECT
      ! Sanity check: regular surface flux must be on an open BC
      IF (PartBound%TargetBoundCond(Species(iSpec)%Surfaceflux(iSF)%BC).EQ.PartBound%ReflectiveBC) THEN
        IF(.NOT.Species(iSpec)%Surfaceflux(iSF)%CircularInflow.AND.(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType.NE.4)) THEN
          CALL abort(__STAMP__&
            ,'ERROR in adaptive surface flux: using a reflective BC without circularInflow is only allowed for Type 4!')
        END IF
      END IF
    ELSE
      Species(iSpec)%Surfaceflux(iSF)%AdaptiveType = 0
    END IF
    ! ================================= ADAPTIVE BC READ IN END ===================================================================!
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
USE MOD_SurfaceModel_Chemistry
USE MOD_SurfaceModel_Vars
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
    IF (SurfFluxSideSize(1).NE.1 .OR. SurfFluxSideSize(2).NE.2) CALL abort(&
__STAMP__&
, 'SurfFluxSideSize must be 1,2 for TriaSurfaceFlux!')
    DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
      CALL CalcNormAndTangTriangle(SideID=SideID,nVec=tmp_Vec_nOut(:,iSample,jSample) &
        ,tang1=tmp_Vec_t1(:,iSample,jSample) &
        ,tang2=tmp_Vec_t2(:,iSample,jSample) &
        ,area=tmp_SubSideAreas(iSample,jSample) &
        ,TriNum=jSample)
      SurfMeshSideAreas(BCSideID)=SurfMeshSideAreas(BCSideID)+tmp_SubSideAreas(iSample,jSample)
    END DO; END DO
  ELSE
    IF (ANY(SurfFluxSideSize.NE.BezierSampleN)) CALL abort(&
__STAMP__&
, 'SurfFluxSideSize must be BezierSampleN,BezierSampleN for .NOT.TriaSurfaceFlux!')
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
IPWRITE(*,*)" ===== TOTAL AREA (all BCsides) ====="
IPWRITE(*,*)"totalArea       = ",totalArea
IPWRITE(*,*)"totalArea/(pi) = ",totalArea/(ACOS(-1.))
IPWRITE(*,*)" ===== TOTAL AREA (all BCsides) ====="
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
USE MOD_SurfaceModel_Chemistry
USE MOD_SurfaceModel_Vars
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*USE_MPI*/
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
INTEGER               :: iReac, SurfNumReac
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

  SurfNumReac = SurfChemReac%NumOfReact
  DO iReac=1,SurfNumReac
    DO iSF=1,SurfChemReac%NumOfBounds(iReac)
      IF (TmpMapToBC(iBC).EQ.SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%BC) THEN !only surfacefluxes with iBC
        ALLOCATE(SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%SurfFluxSubSideData(SurfFluxSideSize(1),SurfFluxSideSize(2),1:TmpSideNumber(iBC)))
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
        IF (SideType(CNSideID).NE.PLANAR_RECT .AND. SideType(CNSideID).NE.PLANAR_NONRECT) CALL abort(&
__STAMP__&
,'every surfaceflux-sides must be planar if TriaSurfaceFlux is used for tracing or refmapping!!!')
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
        CALL abort(&
__STAMP__&
,'Someting is wrong with TmpSideNumber of iBC',iBC,999.)
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
   CALL MPI_ALLREDUCE(areasLoc,areasGlob,nPartBound,MPI_DOUBLE_PRECISION,MPI_SUM,PartMPI%COMM,IERROR)
#endif
   DO iPartBound=1,nPartBound
#if USE_MPI
     BCdata_auxSF(iPartBound)%GlobalArea=areasGlob(iPartBound)
#else
     BCdata_auxSF(iPartBound)%GlobalArea=BCdata_auxSF(iPartBound)%LocalArea
#endif
!     IPWRITE(*,'(I4,A,I4,2(x,E16.8))') 'areas:-', &
!       iPartBound,BCdata_auxSF(iPartBound)%GlobalArea,BCdata_auxSF(iPartBound)%LocalArea
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


SUBROUTINE InitSurfFlux(iSpec, iSF, iSide, tmp_SubSideAreas, BCdata_auxSFTemp)
!===================================================================================================================================
!> Initialize surface flux variables in SurfFluxSubSideData type
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
    projFak = DOT_PRODUCT(vec_nIn,Species(iSpec)%Surfaceflux(iSF)%VeloVecIC) !VeloVecIC projected to inwards normal
  ELSE
    projFak = 1.
  END IF
  v_thermal = SQRT(2.*BoltzmannConst*Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC/Species(iSpec)%MassIC) !thermal speed
  a = 0 !dummy for projected speed ratio in constant v-distri
  !-- compute total volume flow rate through surface
  SELECT CASE(TRIM(Species(iSpec)%Surfaceflux(iSF)%velocityDistribution))
  CASE('constant')
    vSF = Species(iSpec)%Surfaceflux(iSF)%VeloIC * projFak !Velo proj. to inwards normal
    nVFR = MAX(tmp_SubSideAreas(iSample,jSample) * vSF,0.) !VFR proj. to inwards normal (only positive parts!)
  CASE('maxwell','maxwell_lpn')
    IF ( ALMOSTEQUAL(v_thermal,0.)) THEN
      CALL abort(__STAMP__,' ERROR in SurfaceFlux: Calculated thermal velocity is zero! Temperature input might be missing (-MWTemperatureIC) ')
    END IF
    a = Species(iSpec)%Surfaceflux(iSF)%VeloIC * projFak / v_thermal !speed ratio proj. to inwards n (can be negative!)
    vSF = v_thermal / (2.0*SQRT(PI)) * ( EXP(-(a*a)) + a*SQRT(PI)*(1+ERF(a)) ) !mean flux velocity through normal sub-face
    nVFR = tmp_SubSideAreas(iSample,jSample) * vSF !VFR projected to inwards normal of sub-side
    IF(RadialWeighting%DoRadialWeighting) THEN
      nVFR = nVFR / BCdata_auxSFTemp(currentBC)%WeightingFactor(iSide)
      DO iSub = 1, RadialWeighting%nSubSides
        IF(ABS(BCdata_auxSFTemp(currentBC)%SubSideWeight(iSide,iSub)).GT.0.)THEN
          Species(iSpec)%Surfaceflux(iSF)%nVFRSub(iSide,iSub) = BCdata_auxSFTemp(currentBC)%SubSideArea(iSide,iSub) &
                                                             * vSF / BCdata_auxSFTemp(currentBC)%SubSideWeight(iSide,iSub)
        END IF
      END DO
    END IF
  CASE DEFAULT
    CALL abort(__STAMP__,&
      'ERROR in SurfaceFlux: Wrong velocity distribution!')
  END SELECT
  IF (Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN !check rmax-rejection
    IF (Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide).EQ.1) THEN ! complete side is outside of valid bounds
      nVFR = 0.
    END IF
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

END SUBROUTINE InitSurfFlux


SUBROUTINE InitReduceNoiseSF(iSpec, iSF)
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: Species
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*USE_MPI*/
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
  ,Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcs,1,MPI_DOUBLE_PRECISION,0,PartMPI%COMM,iError)
IF(MPIroot)THEN
  DO iProc=0,PartMPI%nProcs-1
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

SUBROUTINE InitSurfChemFlux(iReac, iSF, iSide, tmp_SubSideAreas, BCdata_auxSFTemp)
  !===================================================================================================================================
  !> Initialize surface flux variables in SurfFluxSubSideData type
  !===================================================================================================================================
  ! MODULES
  USE MOD_Globals
  USE MOD_Globals_Vars            ,ONLY: BoltzmannConst, PI
  USE MOD_Particle_Surfaces_Vars  ,ONLY: SurfFluxSideSize, SurfMeshSubSideData, tBCdata_auxSFRadWeight, BCdata_auxSF
  USE MOD_Particle_Vars           ,ONLY: Species, nSpecies
  USE MOD_SurfaceModel_Vars
  USE MOD_DSMC_Vars               ,ONLY: RadialWeighting
  ! IMPLICIT VARIABLE HANDLING
   IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  INTEGER, INTENT(IN)   :: iReac, iSF, iSide
  REAL, INTENT(IN)      :: tmp_SubSideAreas(SurfFluxSideSize(1),SurfFluxSideSize(2))
  TYPE(tBCdata_auxSFRadWeight), ALLOCATABLE, INTENT(IN)        :: BCdata_auxSFTemp(:)
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
  INTEGER               ::iSpec
  INTEGER               :: jSample, iSample, iSub, currentBC, BCSideID
  REAL                  :: vec_nIn(3), nVFR, vec_t1(3), vec_t2(3), projFak, v_thermal, a, vSF
  !===================================================================================================================================
  currentBC = SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%BC
  BCSideID=BCdata_auxSF(currentBC)%SideList(iSide)
  DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
    vec_nIn = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn
    vec_t1 = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t1
    vec_t2 = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t2
    IF (.NOT.SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloIsNormal) THEN
      projFak = DOT_PRODUCT(vec_nIn,SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloVecIC) !VeloVecIC projected to inwards normal
    ELSE
      projFak = 1.
    END IF

    DO iSpec=1,nSpecies  
      IF(ANY(SurfChemReac%Products(iReac,:).EQ.iSpec)) THEN
        v_thermal = SQRT(2.*BoltzmannConst*SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%MWTemperatureIC/Species(iSpec)%MassIC) !thermal speed
      ELSE 
        v_thermal = 0.
      END IF
    END DO
  
    a = 0 !dummy for projected speed ratio in constant v-distri
    !-- compute total volume flow rate through surface
    SELECT CASE(TRIM(SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%velocityDistribution))
    CASE('constant')
      vSF = SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloIC * projFak !Velo proj. to inwards normal
      nVFR = MAX(tmp_SubSideAreas(iSample,jSample) * vSF,0.) !VFR proj. to inwards normal (only positive parts!)
    CASE('maxwell','maxwell_lpn')
      IF(v_thermal.NE.0.) THEN
      a = SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloIC * projFak / v_thermal !speed ratio proj. to inwards n (can be negative!)
    ELSE 
      a = 0.
    END IF
      vSF = v_thermal / (2.0*SQRT(PI)) * ( EXP(-(a*a)) + a*SQRT(PI)*(1+ERF(a)) ) !mean flux velocity through normal sub-face
      nVFR = tmp_SubSideAreas(iSample,jSample) * vSF !VFR projected to inwards normal of sub-side
      IF(RadialWeighting%DoRadialWeighting) THEN
        nVFR = nVFR / BCdata_auxSFTemp(currentBC)%WeightingFactor(iSide)
        DO iSub = 1, RadialWeighting%nSubSides
          IF(ABS(BCdata_auxSFTemp(currentBC)%SubSideWeight(iSide,iSub)).GT.0.)THEN
            SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%nVFRSub(iSide,iSub) = BCdata_auxSFTemp(currentBC)%SubSideArea(iSide,iSub) &
                                                               * vSF / BCdata_auxSFTemp(currentBC)%SubSideWeight(iSide,iSub)
          END IF
        END DO
      END IF
    CASE DEFAULT
      CALL abort(__STAMP__,&
        'ERROR in SurfaceFlux: Wrong velocity distribution!')
    END SELECT
    
    SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VFR_total = SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VFR_total + nVFR
    !-- store SF-specific SubSide data in SurfFluxSubSideData (incl. projected velos)
    SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%nVFR = nVFR
    SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%projFak = projFak
    SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%a_nIn = a
    IF (.NOT.SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloIsNormal) THEN
      SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t1 &
        = SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloIC &
        * DOT_PRODUCT(vec_t1,SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloVecIC) !v in t1-dir
      SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t2 &
        = SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloIC &
        * DOT_PRODUCT(vec_t2,SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloVecIC) !v in t2-dir
    ELSE
      SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t1 = 0. !v in t1-dir
      SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t2 = 0. !v in t2-dir
    END IF! .NOT.VeloIsNormal
  END DO; END DO !jSample=1,SurfFluxSideSize(2); iSample=1,SurfFluxSideSize(1)
  
  END SUBROUTINE InitSurfChemFlux

  SUBROUTINE ReadInAndPrepareSurfChemFlux(nDataBC)
    !===================================================================================================================================
    ! Initialize the variables first
    !===================================================================================================================================
    ! MODULES
    USE MOD_Globals
    USE MOD_ReadInTools
    USE MOD_Globals_Vars           ,ONLY: BoltzmannConst, Pi
    USE MOD_Particle_Boundary_Vars ,ONLY: PartBound,nPartBound
    USE MOD_SurfaceModel_Tools     ,ONLY: GetWallTemperature
    USE MOD_SurfaceModel_Vars
    USE MOD_Particle_Surfaces_Vars ,ONLY: BCdata_auxSF, BezierSampleN, SurfMeshSubSideData, SurfMeshSideAreas, tBCdata_auxSFRadWeight
    USE MOD_Particle_Surfaces_Vars ,ONLY: SurfFluxSideSize, TriaSurfaceFlux
    USE MOD_Particle_Surfaces      ,ONLY: GetBezierSampledAreas
    USE MOD_Particle_Vars          ,ONLY: Species, nSpecies, DoSurfaceFlux
    USE MOD_Particle_Vars          ,ONLY: UseCircularInflow, DoForceFreeSurfaceFlux
    USE MOD_Particle_Sampling_Vars ,ONLY: UseAdaptive
    USE MOD_Restart_Vars           ,ONLY: DoRestart, RestartTime
    ! IMPLICIT VARIABLE HANDLING
     IMPLICIT NONE
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    INTEGER, INTENT(INOUT) :: nDataBC
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    INTEGER               :: iSF
    INTEGER               :: iReac, SurfNumReac
    !===================================================================================================================================
    SurfNumReac = SurfChemReac%NumOfReact
    ALLOCATE(SurfChemReac%SFMap(SurfChemReac%NumOfReact))
    DO iReac=1,SurfNumReac
      IF (SurfChemReac%NumOfBounds(iReac).EQ.0) THEN
        CYCLE
      ELSE
        ALLOCATE(SurfChemReac%SFMap(iReac)%Surfaceflux(1:SurfChemReac%NumOfBounds(iReac)))
        ! Initialize Surfaceflux to BC mapping
        SurfChemReac%SFMap(iReac)%Surfaceflux(:)%BC=-1
        DO iSF=1,SurfChemReac%NumOfBounds(iReac)
          SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%BC = SurfChemReac%BoundMap(iReac)%Boundaries(iSF)
        END DO
      END IF
    
      DO iSF=1,SurfChemReac%NumOfBounds(iReac)
        IF (TriaSurfaceFlux) THEN
          SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%AcceptReject=.FALSE.
        END IF

        SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%InsertedParticle = 0
        SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VFR_total = 0
        SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VFR_total_allProcsTotal = 0
    
        ! get surfaceflux data
        IF (SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%BC.LT.1 .OR. SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%BC.GT.nPartBound) THEN
          CALL abort(&
    __STAMP__&
    , 'SurfacefluxBCs must be between 1 and nPartBound!')
        ELSE IF (BCdata_auxSF(SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%BC)%SideNumber.EQ. -1) THEN !not set yet
          BCdata_auxSF(SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%BC)%SideNumber=0
          nDataBC=nDataBC+1
        END IF
    
        SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%velocityDistribution = 'maxwell_lpn'
        SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloIC = 0.
        SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloIsNormal = .FALSE.
        SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%MWTemperatureIC = PartBound%WallTemp(SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%BC)
      
      END DO !iSF
    END DO ! iReac
    
    END SUBROUTINE ReadInAndPrepareSurfChemFlux

END MODULE MOD_Particle_SurfFlux_Init
