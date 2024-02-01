!==================================================================================================================================
! Copyright (c) 2010 - 2019 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_Particle_SurfFlux
!===================================================================================================================================
!> Module for particle insertion through the surface flux
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: ParticleSurfaceflux
!===================================================================================================================================
CONTAINS

!===================================================================================================================================
!> Particle Inserting via Surface Flux and (if present) adaptiveBC (Surface Flux adapting part density, velocity or temperature)
!===================================================================================================================================
SUBROUTINE ParticleSurfaceflux()
! Modules
USE MOD_Globals
USE MOD_Particle_Vars
USE MOD_part_tools              ,ONLY: CalcRadWeightMPF, IncreaseMaxParticleNumber
USE MOD_DSMC_Vars               ,ONLY: useDSMC, CollisMode, RadialWeighting, DSMC
USE MOD_Eval_xyz                ,ONLY: GetPositionInRefElem
USE MOD_Mesh_Vars               ,ONLY: SideToElem, offsetElem
USE MOD_Part_Tools              ,ONLY: GetParticleWeight, GetNextFreePosition
USE MOD_Part_Emission_Tools     ,ONLY: SetParticleChargeAndMass, SetParticleMPF
USE MOD_Particle_Analyze_Vars   ,ONLY: CalcPartBalance, CalcSurfFluxInfo, nPartIn, PartEkinIn
USE MOD_Particle_Analyze_Tools  ,ONLY: CalcEkinPart
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Sampling_Vars  ,ONLY: AdaptBCPartNumOut
USE MOD_Particle_Surfaces_Vars  ,ONLY: SurfFluxSideSize, TriaSurfaceFlux, BCdata_auxSF
USE MOD_Particle_TimeStep       ,ONLY: GetParticleTimeStep
USE MOD_Timedisc_Vars           ,ONLY: RKdtFrac, dt
USE MOD_DSMC_AmbipolarDiffusion ,ONLY: AD_SetSFElectronVelo
#if defined(IMPA) || defined(ROS)
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
#endif /*IMPA*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars        ,ONLY: nSurfacefluxPerElem
USE MOD_LoadBalance_Timers      ,ONLY: LBStartTime, LBElemSplitTime, LBPauseTime
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
INTEGER                     :: iSpec , PositionNbr, iSF, iSide, currentBC, SideID, NbrOfParticle, ExtraParts, ParticleIndexNbr
INTEGER                     :: BCSideID, ElemID, iLocSide, iSample, jSample, PartInsSubSide, iPart, iPartTotal
INTEGER                     :: allowedRejections, PartsEmitted, Node1, Node2, globElemId
INTEGER                     :: PartInsSideRadWeight(1:RadialWeighting%nSubSides)
REAL                        :: Particle_pos(3), RandVal1,  xyzNod(3), RVec(2), minPos(2), xi(2), Vector1(3), Vector2(3)
REAL                        :: ndist(3), midpoint(3)
LOGICAL                     :: AcceptPos
REAL,ALLOCATABLE            :: particle_positions(:), particle_xis(:)
INTEGER,ALLOCATABLE         :: PartInsSubSides(:,:,:)
#if USE_LOADBALANCE
REAL                        :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
DO iSpec=1,nSpecies
  IF(useDSMC) THEN
    IF (DSMC%DoAmbipolarDiff) THEN
      IF (iSpec.EQ.DSMC%AmbiDiffElecSpec) CYCLE
    END IF
  END IF
  DO iSF=1,Species(iSpec)%nSurfacefluxBCs
    PartsEmitted = 0
    ASSOCIATE(SF => Species(iSpec)%Surfaceflux(iSF))
    currentBC = SF%BC
    NbrOfParticle = 0 ! calculated within (sub)side-Loops!
    iPartTotal=0
    ! Adaptive BC, Type = 4 (Const. massflow): Sum-up the global number of particles exiting through BC and calculate new weights
    IF(SF%AdaptiveType.EQ.4) THEN
#if USE_MPI
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,AdaptBCPartNumOut(iSpec,iSF),1,MPI_INTEGER,MPI_SUM,MPI_COMM_PICLAS,IERROR)
#endif
      IF(.NOT.ALMOSTEQUAL(SF%AdaptiveMassflow,0.)) CALL CalcConstMassflowWeight(iSpec,iSF)
    END IF
    ! Calc Particles for insertion in standard case
    IF (SF%Type.EQ.0) CALL CalcPartInsSubSidesStandardCase(iSpec,iSF, PartInsSubSides)

!----- 0.: go through (sub)sides if present in proc
    IF (BCdata_auxSF(currentBC)%SideNumber.EQ.0) THEN
      IF(SF%AdaptiveType.EQ.4) AdaptBCPartNumOut(iSpec,iSF) = 0
      CYCLE
    ELSE IF (BCdata_auxSF(currentBC)%SideNumber.EQ.-1) THEN
      CALL abort(__STAMP__,'ERROR in ParticleSurfaceflux: Someting is wrong with SideNumber of BC ',currentBC)
    END IF
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
    DO iSide=1,BCdata_auxSF(currentBC)%SideNumber
      IF (SF%CircularInflow) THEN
        IF(SF%SurfFluxSideRejectType(iSide).EQ.1) CYCLE
      END IF
      BCSideID=BCdata_auxSF(currentBC)%SideList(iSide)
      ElemID = SideToElem(S2E_ELEM_ID,BCSideID)
      iLocSide = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
      globElemId = ElemID + offSetElem
      SideID=GetGlobalNonUniqueSideID(globElemId,iLocSide)
      IF (TriaSurfaceFlux) xyzNod(1:3) = BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%xyzNod(1:3)
      DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
        ExtraParts = 0 !set here number of additional to-be-inserted particles in current BCSideID/subsides (e.g. desorption)
        IF(Symmetry%Axisymmetric.AND.(jSample.EQ.2)) CYCLE
        IF (TriaSurfaceFlux) THEN
          !-- compute parallelogram of triangle
          Node1 = jSample+1     ! normal = cross product of 1-2 and 1-3 for first triangle
          Node2 = jSample+2     !          and 1-3 and 1-4 for second triangle
          Vector1 = BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%Vectors(:,Node1-1)
          Vector2 = BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%Vectors(:,Node2-1)
          midpoint(1:3) = BCdata_auxSF(currentBC)%TriaSwapGeo(iSample,jSample,iSide)%midpoint(1:3)
          ndist(1:3) = BCdata_auxSF(currentBC)%TriaSwapGeo(iSample,jSample,iSide)%ndist(1:3)
        END IF

        ! REQUIRED LATER FOR THE POSITION START
        IF(Symmetry%Axisymmetric) CALL DefineSideDirectVec2D(SideID, xyzNod, minPos, RVec)

        !-- compute number of to be inserted particles
        SELECT CASE(SF%Type)
        CASE(0)
          ! Standard surface flux with fixed number of particles per side
          PartInsSubSide=PartInsSubSides(iSample,jSample,iSide)
        CASE(1)
          ! Adaptive surface flux: Number of particles depends on velocity/temperature/mass flow (includes treatment for RadialWeighting)
          CALL CalcPartInsAdaptive(iSpec, iSF, BCSideID, iSide, iSample, jSample, minPos, RVec, PartInsSubSide, PartInsSideRadWeight)
        CASE(2)
          ! Radial weighting: Number of particles depends on the modified area and includes insertion over subsides
          CALL CalcPartInsRadWeight(iSpec, iSF, iSample, jSample, iSide, minPos, RVec, PartInsSubSide, PartInsSideRadWeight)
        CASE(3)
          ! DoPoissonRounding .AND. .NOT.DoTimeDepInflow
          CALL CalcPartInsPoissonDistr(iSpec, iSF, iSample, jSample, iSide, PartInsSubSide)
        CASE(4)
          ! DoTimeDepInflow
          CALL RANDOM_NUMBER(RandVal1)
          PartInsSubSide = INT(SF%PartDensity / Species(iSpec)%MacroParticleFactor &
                          * dt*RKdtFrac * SF%SurfFluxSubSideData(iSample,jSample,iSide)%nVFR+RandVal1)
#if (USE_HDG)
        CASE(5)
          ! SF%ThermionicEmission.AND.SF%SchottkyEffectTE
          CALL CalcPartInsThermionicEmissionSchottky(iSpec,iSF,iSample,jSample,iSide,iLocSide,ElemID,PartInsSubSide)
#endif
        CASE DEFAULT
          CALL abort(__STAMP__,'ERROR in ParticleSurfaceflux: Given surface flux type is not defined!')
        END SELECT

        !-- proceed with calculated to be inserted particles
        IF (PartInsSubSide.LT.0) THEN
          IPWRITE(*,*) 'ERROR in ParticleSurfaceflux: Calculated number of particles to insert below zero! PartInsSubSide: ', PartInsSubSide
          IPWRITE(*,*) 'Species Index: ', iSpec, 'SurfaceFlux Index: ', iSF, 'iSide: ', iSide
          CALL abort(__STAMP__,'ERROR in ParticleSurfaceflux: PartInsSubSide.LT.0!')
        ELSE IF (PartInsSubSide + ExtraParts.LE.0) THEN
          CYCLE
        END IF
        PartInsSubSide = PartInsSubSide + ExtraParts
        NbrOfParticle = NbrOfParticle + PartInsSubSide
        ALLOCATE(particle_positions(1:PartInsSubSide*3))
        IF (SF%VeloIsNormal .AND. .NOT.TriaSurfaceFlux) THEN
          ALLOCATE( particle_xis(1:PartInsSubSide*2))
        END IF !VeloIsNormal
        !-- put particles in subside (rejections are used if contraint reduces actual inserted number)
        iPart=1
        allowedRejections=0

        !-- Set Positions
        IF(Symmetry%Axisymmetric) THEN
          CALL CalcPartPosAxisym(iSpec, iSF, iSide, minPos, RVec, PartInsSubSide, PartInsSideRadWeight, particle_positions, allowedRejections)
        ELSE
          DO WHILE (iPart+allowedRejections .LE. PartInsSubSide)
            IF (TriaSurfaceFlux) THEN
              Particle_pos(1:3) = CalcPartPosTriaSurface(xyzNod, Vector1, Vector2, ndist, midpoint)
            ELSE !.NOT.TriaSurfaceFlux
              Particle_pos(1:3) = CalcPartPosBezier(iSpec, iSF, iSample, jSample, iSide, SideID, xi)
            END IF !TriaSurfaceFlux

            AcceptPos=.TRUE.
            IF (SF%CircularInflow) THEN !check rmax-rejection
              IF (.NOT.InSideCircularInflow(iSpec, iSF, iSide, Particle_pos)) AcceptPos=.FALSE.
            END IF ! CircularInflow
            !-- save position if accepted:
            IF (AcceptPos) THEN
              particle_positions(iPart*3-2) = Particle_pos(1)
              particle_positions(iPart*3-1) = Particle_pos(2)
              particle_positions(iPart*3  ) = Particle_pos(3)
              IF (SF%VeloIsNormal .AND. .NOT.TriaSurfaceFlux) THEN
                particle_xis(iPart*2-1) = xi(1)
                particle_xis(iPart*2  ) = xi(2)
              END IF !VeloIsNormal
              iPart=iPart+1
            ELSE
              IF (SF%CircularInflow) THEN !check rmax-rejection
                allowedRejections=allowedRejections+1
              END IF
            END IF
          END DO !put particles in subside: WHILE(iPart+allowedRejections .LE. PartInsSubSide)
        END IF
        PartInsSubSide = PartInsSubSide - allowedRejections
        NbrOfParticle = NbrOfParticle - allowedRejections

        !-- Fill Particle Information (PartState, Partelem, etc.)
        ParticleIndexNbr = 1
        DO iPart=1,PartInsSubSide
          IF ((iPart.EQ.1).OR.PDM%ParticleInside(ParticleIndexNbr)) &
              ParticleIndexNbr = GetNextFreePosition(iPartTotal+1)
          PartState(1:3,ParticleIndexNbr) = particle_positions(3*(iPart-1)+1:3*(iPart-1)+3)
          IF (SF%VeloIsNormal.AND.(.NOT.TriaSurfaceFlux)) THEN
            PartState(4:5,ParticleIndexNbr) = particle_xis(2*(iPart-1)+1:2*(iPart-1)+2) !use velo as dummy-storage for xi!
          END IF
          LastPartPos(1:3,ParticleIndexNbr)=PartState(1:3,ParticleIndexNbr)
#if defined(IMPA) || defined(ROS)
          IF(TrackingMethod.EQ.REFMAPPING) CALL GetPositionInRefElem(PartState(1:3,ParticleIndexNbr),PartPosRef(1:3,ParticleIndexNbr),globElemId)
#endif /*IMPA*/
          PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
          PDM%dtFracPush(ParticleIndexNbr) = .TRUE.
          PDM%IsNewPart(ParticleIndexNbr) = .TRUE.
          PEM%GlobalElemID(ParticleIndexNbr) = globElemId
          PEM%LastGlobalElemID(ParticleIndexNbr) = globElemId !needed when ParticlePush is not executed, e.g. "delay"
          iPartTotal = iPartTotal + 1
          IF (UseVarTimeStep) THEN
            PartTimeStep(ParticleIndexNbr) = GetParticleTimeStep(PartState(1,ParticleIndexNbr),PartState(2,ParticleIndexNbr), &
                                                              PEM%LocalElemID(ParticleIndexNbr))
          END IF
          IF (RadialWeighting%DoRadialWeighting) THEN
            PartMPF(ParticleIndexNbr) = CalcRadWeightMPF(PartState(2,ParticleIndexNbr), iSpec,ParticleIndexNbr)
          END IF
          IF(CalcSurfFluxInfo) THEN
            SF%SampledMassflow = SF%SampledMassflow + GetParticleWeight(ParticleIndexNbr)
          END IF
#ifdef CODE_ANALYZE
          CALL AnalyzePartPos(ParticleIndexNbr)
#endif /*CODE_ANALYZE*/
        END DO
!----- 2a.: set velocities if special for each subside
        CALL SetSurfacefluxVelocities(iSpec,iSF,iSample,jSample,iSide,BCSideID,SideID,ElemID,NbrOfParticle,PartInsSubSide)

        PartsEmitted = PartsEmitted + PartInsSubSide

        IF (useDSMC) THEN
          IF (DSMC%DoAmbipolarDiff) CALL AD_SetSFElectronVelo(iSpec,iSF,iSample,jSample,iSide,BCSideID,SideID,ElemID,NbrOfParticle,PartInsSubSide,particle_xis)
        END IF

        IF (SF%VeloIsNormal .AND. .NOT.TriaSurfaceFlux) DEALLOCATE(particle_xis)
        DEALLOCATE(particle_positions)
#if USE_LOADBALANCE
        !used for calculating LoadBalance of tCurrent(LB_SURFFLUX) ==> "2b.: set remaining properties"
        nSurfacefluxPerElem(ElemID)=nSurfacefluxPerElem(ElemID)+PartInsSubSide
#endif /*USE_LOADBALANCE*/

      END DO; END DO !jSample=1,SurfFluxSideSize(2); iSample=1,SurfFluxSideSize(1)
#if USE_LOADBALANCE
      CALL LBElemSplitTime(ElemID,tLBStart)
#endif /*USE_LOADBALANCE*/
    END DO ! iSide

    IF(SF%Adaptive) THEN
      IF(SF%AdaptiveType.EQ.4) AdaptBCPartNumOut(iSpec,iSF) = 0
    END IF
    IF (NbrOfParticle.NE.iPartTotal) CALL abort(__STAMP__, 'Error ParticleSurfaceflux: Mismatch between the determined and inserted number of particles!')
!----- 2b.: set remaining properties
    CALL SetParticleChargeAndMass(iSpec,NbrOfParticle)
    IF (usevMPF.AND.(.NOT.RadialWeighting%DoRadialWeighting)) CALL SetParticleMPF(iSpec,-1,NbrOfParticle)
    ! define molecule stuff
    IF (useDSMC.AND.(CollisMode.GT.1)) CALL SetInnerEnergies(iSpec, iSF, NbrOfParticle)
    IF(CalcPartBalance) THEN
    ! Compute number of input particles and energy
      nPartIn(iSpec)=nPartIn(iSpec) + NbrOfParticle
      DO iPart=1,NbrOfParticle
        PositionNbr = GetNextFreePosition(iPart)
        PartEkinIn(iSpec) = PartEkinIn(iSpec)+CalcEkinPart(PositionNbr)
      END DO ! iPart
    END IF ! CalcPartBalance
    ! instead of an UpdateNextfreePosition we update the particleVecLength only - enough ?!?
    IF(iPartTotal.GT.0) THEN
      PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + NbrOfParticle
      PDM%ParticleVecLength = MAX(PDM%ParticleVecLength,GetNextFreePosition(0))
    END IF
#ifdef CODE_ANALYZE
    IF(PDM%ParticleVecLength.GT.PDM%maxParticleNumber) CALL Abort(__STAMP__,'PDM%ParticleVeclength exceeds PDM%maxParticleNumber, Difference:',IntInfoOpt=PDM%ParticleVeclength-PDM%maxParticleNumber)
    DO iPart=PDM%ParticleVecLength+1,PDM%maxParticleNumber
      IF (PDM%ParticleInside(iPart)) THEN
        IPWRITE(*,*) iPart,PDM%ParticleVecLength,PDM%maxParticleNumber
        CALL Abort(__STAMP__,'Particle outside PDM%ParticleVeclength',IntInfoOpt=iPart)
      END IF
    END DO
#endif
#if USE_LOADBALANCE
    CALL LBPauseTime(LB_SURFFLUX,tLBStart)
#endif /*USE_LOADBALANCE*/
    ! Sample Energies on Surfaces when particles are emitted from them
    IF (NbrOfParticle.NE.PartsEmitted) THEN
      ! should be equal for including the following lines in tSurfaceFlux
      CALL abort(__STAMP__,'ERROR in ParticleSurfaceflux: NbrOfParticle.NE.PartsEmitted')
    END IF
    END ASSOCIATE
  END DO !iSF
END DO !iSpec

END SUBROUTINE ParticleSurfaceflux


!===================================================================================================================================
!> Calculate the particle number per side for the case of a regular emission
!===================================================================================================================================
SUBROUTINE CalcPartInsSubSidesStandardCase(iSpec, iSF, PartInsSubSides)
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars           ,ONLY: Species, VarTimeStep
USE MOD_TimeDisc_Vars           ,ONLY: dt, RKdtFrac, RKdtFracTotal, Time
USE MOD_Particle_Surfaces_Vars  ,ONLY: SurfFluxSideSize, BCdata_auxSF
USE MOD_Part_Emission_Tools     ,ONLY: IntegerDivide, SamplePoissonDistri
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                 :: iSpec, iSF
INTEGER, INTENT(OUT), ALLOCATABLE   :: PartInsSubSides(:,:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(KIND=8)        :: inserted_Particle_iter,inserted_Particle_time,inserted_Particle_diff
INTEGER                :: currentBC, PartInsSF, IntSample
REAL                   :: VFR_total, PartIns, RandVal1
INTEGER, ALLOCATABLE   :: PartInsProc(:)
REAL                   :: dtVar, TimeVar
!===================================================================================================================================

ASSOCIATE(SF => Species(iSpec)%Surfaceflux(iSF))
!--- Noise reduction (both ReduceNoise=T (with comm.) and F (proc local), but not for DoPoissonRounding)
currentBC = SF%BC
IF (SF%ReduceNoise) THEN
  !-- calc global to-be-inserted number of parts and distribute to procs (root)
  ALLOCATE(PartInsProc(0:nProcessors-1))
  PartInsProc=0
END IF !ReduceNoise
IF (.NOT.SF%ReduceNoise .OR. MPIroot) THEN !ReduceNoise: root only
  IF (SF%ReduceNoise) THEN
    VFR_total = SF%VFR_total_allProcsTotal !proc global total
  ELSE
    VFR_total = SF%VFR_total               !proc local total
  END IF
  ! Species-specific time step
  IF(VarTimeStep%UseSpeciesSpecific) THEN
    dtVar = dt * Species(iSpec)%TimeStepFactor
    TimeVar = Time * Species(iSpec)%TimeStepFactor
  ELSE
    dtVar = dt
    TimeVar = Time
  END IF
  PartIns = SF%PartDensity / Species(iSpec)%MacroParticleFactor * dtVar*RKdtFrac * VFR_total
  inserted_Particle_iter = INT(PartIns,8)
  PartIns = SF%PartDensity / Species(iSpec)%MacroParticleFactor * (TimeVar + dtVar*RKdtFracTotal) * VFR_total
  !-- random-round the inserted_Particle_time for preventing periodicity
  IF (inserted_Particle_iter.GE.1) THEN
    CALL RANDOM_NUMBER(RandVal1)
    inserted_Particle_time = INT(PartIns+RandVal1,8)
  ELSE IF (inserted_Particle_iter.GE.0) THEN
    ! Needed, since InsertedParticleSurplus can increase and _iter>1 needs to be possible for preventing periodicity
    IF (ALMOSTEQUAL(PartIns,0.)) THEN !dummy for procs without SFs (needed for mpi-comm, are cycled later)
      inserted_Particle_time = INT(PartIns,8)
    ELSE !poisson-distri of PartIns-INT(PartIns)
      CALL SamplePoissonDistri( PartIns-INT(PartIns) , IntSample )
      inserted_Particle_time = INT(INT(PartIns)+IntSample,8) !INT(PartIns) + POISDISTRI( PartIns-INT(PartIns) )
    END IF
  ELSE !dummy for procs without SFs (needed for mpi-comm, are cycled later)
    inserted_Particle_time = INT(PartIns,8)
  END IF
  !-- evaluate inserted_Particle_time and inserted_Particle_iter
  inserted_Particle_diff = inserted_Particle_time - SF%InsertedParticle - inserted_Particle_iter - SF%InsertedParticleSurplus
  SF%InsertedParticleSurplus = ABS(MIN(inserted_Particle_iter + inserted_Particle_diff,0_8))
  PartInsSF = MAX(INT(inserted_Particle_iter + inserted_Particle_diff,4),0)
  SF%InsertedParticle = SF%InsertedParticle + INT(PartInsSF,8)
  IF (SF%ReduceNoise) THEN
#if USE_MPI
    CALL IntegerDivide(PartInsSF,nProcessors,SF%VFR_total_allProcs(0:nProcessors-1),PartInsProc(0:nProcessors-1))
#else  /*USE_MPI*/
    PartInsProc=PartInsSF
#endif  /*USE_MPI*/
  END IF !ReduceNoise
END IF !ReduceNoise, MPIroot
#if USE_MPI
IF (SF%ReduceNoise) THEN !scatter PartInsProc into PartInsSF of procs
  CALL MPI_SCATTER(PartInsProc(0:nProcessors-1),1,MPI_INTEGER,PartInsSF,1,MPI_INTEGER,0,MPI_COMM_PICLAS,IERROR)
END IF !ReduceNoise
#endif  /*USE_MPI*/
!-- calc global to-be-inserted number of parts and distribute to SubSides (proc local)
SDEALLOCATE(PartInsSubSides)
ALLOCATE(PartInsSubSides(SurfFluxSideSize(1),SurfFluxSideSize(2),1:BCdata_auxSF(currentBC)%SideNumber))
PartInsSubSides=0
IF (BCdata_auxSF(currentBC)%SideNumber.LT.1) THEN
  IF (PartInsSF.NE.0) CALL abort(__STAMP__,'ERROR in ParticleSurfaceflux: Someting is wrong with PartInsSF of BC ',currentBC)
ELSE
  CALL IntegerDivide(PartInsSF,BCdata_auxSF(currentBC)%SideNumber*SurfFluxSideSize(1)*SurfFluxSideSize(2) &
    ,SF%SurfFluxSubSideData(1:SurfFluxSideSize(1),1:SurfFluxSideSize(2),1:BCdata_auxSF(currentBC)%SideNumber)%nVFR &
    ,PartInsSubSides(1:SurfFluxSideSize(1),1:SurfFluxSideSize(2),1:BCdata_auxSF(currentBC)%SideNumber) )
END IF
END ASSOCIATE

END SUBROUTINE CalcPartInsSubSidesStandardCase


!===================================================================================================================================
!>
!===================================================================================================================================
SUBROUTINE DefineSideDirectVec2D(SideID, xyzNod, minPos, RVec)
! MODULES
USE MOD_Globals
USE MOD_Particle_Mesh_Vars        ,ONLY: NodeCoords_Shared, ElemSideNodeID_Shared, SideInfo_Shared
USE MOD_Mesh_Tools                ,ONLY: GetCNElemID
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                 :: SideID
REAL, INTENT(IN)                    :: xyzNod(3)
REAL, INTENT(OUT)                   :: minPos(2), RVec(2)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iLocSide, CNElemID, Node1, Node2, minVec
REAL                    :: Vector1(3), Vector2(3), Vector2D(2)
!===================================================================================================================================
iLocSide = SideInfo_Shared(SIDE_LOCALID,SideID)
CNElemID = GetCNElemID(SideInfo_Shared(SIDE_ELEMID,SideID))
!-- compute parallelogram of triangle (only simple 2 value adds/subs, other from init)
Node1 = 2     ! normal = cross product of 1-2 and 1-3 for first triangle
Node2 = 4     !          and 1-3 and 1-4 for second triangle
Vector1(1:3) = NodeCoords_Shared(1:3,ElemSideNodeID_Shared(Node1,iLocSide,CNElemID)+1) - xyzNod(1:3)
Vector2(1:3) = NodeCoords_Shared(1:3,ElemSideNodeID_Shared(Node2,iLocSide,CNElemID)+1) - xyzNod(1:3)
IF (ABS(Vector1(3)).GT.ABS(Vector2(3))) THEN
  Vector2D(1:2) = Vector2(1:2)
ELSE
  Vector2D(1:2) = Vector1(1:2)
END IF
minVec = MINLOC((/xyzNod(2), xyzNod(2)+Vector2D(2)/),1)
SELECT CASE(minVec)
CASE(1)
  minPos(1:2) = xyzNod(1:2)
  RVec(1:2) =  Vector2D(1:2)
CASE(2)
  minPos(1:2) = xyzNod(1:2) + Vector2D(1:2)
  RVec(1:2) = - Vector2D(1:2)
END SELECT
END SUBROUTINE DefineSideDirectVec2D


#ifdef CODE_ANALYZE
!===================================================================================================================================
!>
!===================================================================================================================================
SUBROUTINE AnalyzePartPos(ParticleIndexNbr)
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars      ,ONLY: LastPartPos, PDM
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
#ifdef IMPA
USE MOD_Particle_Vars      ,ONLY: PartDtFrac,PartIsImplicit
#endif /*IMPA*/
#if  defined(IMPA) || defined(ROS)
USE MOD_Timedisc_Vars      ,ONLY: iStage
#endif
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                        :: ParticleIndexNbr
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(   (LastPartPos(1,ParticleIndexNbr).GT.GEO%xmaxglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(1,ParticleIndexNbr),GEO%xmaxglob) &
  .OR.(LastPartPos(1,ParticleIndexNbr).LT.GEO%xminglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(1,ParticleIndexNbr),GEO%xminglob) &
  .OR.(LastPartPos(2,ParticleIndexNbr).GT.GEO%ymaxglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(2,ParticleIndexNbr),GEO%ymaxglob) &
  .OR.(LastPartPos(2,ParticleIndexNbr).LT.GEO%yminglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(2,ParticleIndexNbr),GEO%yminglob) &
  .OR.(LastPartPos(3,ParticleIndexNbr).GT.GEO%zmaxglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(3,ParticleIndexNbr),GEO%zmaxglob) &
  .OR.(LastPartPos(3,ParticleIndexNbr).LT.GEO%zminglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(3,ParticleIndexNbr),GEO%zminglob) ) THEN
  IPWRITE(UNIt_stdOut,'(I0,A18,L1)')                            ' ParticleInside ',PDM%ParticleInside(ParticleIndexNbr)
#ifdef IMPA
  IPWRITE(UNIt_stdOut,'(I0,A18,L1)')                            ' PartIsImplicit ', PartIsImplicit(ParticleIndexNbr)
  IPWRITE(UNIt_stdOut,'(I0,A18,ES25.14)')                       ' PartDtFrac ', PartDtFrac(ParticleIndexNbr)
#endif /*IMPA*/
  IPWRITE(UNIt_stdOut,'(I0,A18,L1)')                            ' PDM%IsNewPart ', PDM%IsNewPart(ParticleIndexNbr)
  IPWRITE(UNIt_stdOut,'(I0,A18,1X,A18,1X,A18)')                  '    min ', ' value ', ' max '
  IPWRITE(UNIt_stdOut,'(I0,A2,1X,ES25.14,1X,ES25.14,1X,ES25.14)') ' x', GEO%xminglob, LastPartPos(1,ParticleIndexNbr) &
                                                                , GEO%xmaxglob
  IPWRITE(UNIt_stdOut,'(I0,A2,1X,ES25.14,1X,ES25.14,1X,ES25.14)') ' y', GEO%yminglob, LastPartPos(2,ParticleIndexNbr) &
                                                                , GEO%ymaxglob
  IPWRITE(UNIt_stdOut,'(I0,A2,1X,ES25.14,1X,ES25.14,1X,ES25.14)') ' z', GEO%zminglob, LastPartPos(3,ParticleIndexNbr) &
                                                                , GEO%zmaxglob
  CALL abort(&
     __STAMP__ &
#if  defined(IMPA) || defined(ROS)
     ,' LastPartPos outside of mesh. iPart=, iStage',ParticleIndexNbr,REAL(iStage))
#else
     ,' LastPartPos outside of mesh. iPart=',ParticleIndexNbr)
#endif
END IF
END SUBROUTINE AnalyzePartPos
#endif /*CODE_ANALYZE*/


!===================================================================================================================================
!>
!===================================================================================================================================
SUBROUTINE SetInnerEnergies(iSpec, iSF, NbrOfParticle)
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars               ,ONLY: SpecDSMC
USE MOD_DSMC_PolyAtomicModel    ,ONLY: DSMC_SetInternalEnr_Poly
USE MOD_part_emission_tools     ,ONLY: DSMC_SetInternalEnr_LauxVFD
USE MOD_Part_Tools              ,ONLY: GetNextFreePosition
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                        :: iSpec, iSF, NbrOfParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iPart, PositionNbr
!===================================================================================================================================
iPart = 1
DO iPart=1,NbrOfParticle
  PositionNbr = GetNextFreePosition(iPart)
    IF (SpecDSMC(iSpec)%PolyatomicMol) THEN
      CALL DSMC_SetInternalEnr_Poly(iSpec,iSF,PositionNbr,2)
    ELSE
      CALL DSMC_SetInternalEnr_LauxVFD(iSpec, iSF, PositionNbr,2)
    END IF
END DO
END SUBROUTINE SetInnerEnergies


!===================================================================================================================================
!> Circular inflow: check whether particle position is inside or outside of the defined circle
!===================================================================================================================================
FUNCTION InSideCircularInflow(iSpec, iSF, iSide, Particle_pos)
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars           ,ONLY: Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES)
INTEGER, INTENT(IN)             :: iSpec, iSF, iSide
REAL, INTENT(IN)                :: Particle_pos(3)
LOGICAL                         :: InSideCircularInflow
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                        :: point(2), radius, origin(2)
!===================================================================================================================================
InSideCircularInflow=.FALSE. ! Initialize
origin=Species(iSpec)%Surfaceflux(iSF)%origin
SELECT CASE(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide))
CASE(0) !- RejectType=0 : complete side is inside valid bounds
  InSideCircularInflow=.TRUE.
CASE(1) !- RejectType=1 : complete side is outside of valid bounds
  CALL abort(__STAMP__,'side outside of valid bounds was considered although nVFR=0...?!')
CASE(2) !- RejectType=2 : side is partly inside valid bounds
  point(1)=Particle_pos(Species(iSpec)%Surfaceflux(iSF)%dir(2))-origin(1)
  point(2)=Particle_pos(Species(iSpec)%Surfaceflux(iSF)%dir(3))-origin(2)
  radius=SQRT( (point(1))**2+(point(2))**2 )
  IF ((radius.LE.Species(iSpec)%Surfaceflux(iSF)%rmax).AND.(radius.GE.Species(iSpec)%Surfaceflux(iSF)%rmin)) THEN
    InSideCircularInflow=.TRUE.
  ELSE
    InSideCircularInflow=.FALSE.
  END IF
CASE DEFAULT
  CALL abort(__STAMP__,'wrong SurfFluxSideRejectType!')
END SELECT !SurfFluxSideRejectType

END FUNCTION InSideCircularInflow


!===================================================================================================================================
!> Calculate random normalized vector in 3D (unit space)
!===================================================================================================================================
FUNCTION CalcPartPosBezier(iSpec, iSF, iSample, jSample, iSide, SideID, xi)
! MODULES
USE MOD_Globals
USE MOD_Particle_Surfaces_Vars  ,ONLY: BezierControlPoints3D,BezierSampleXi
USE MOD_Particle_Surfaces       ,ONLY: EvaluateBezierPolynomialAndGradient
USE MOD_Particle_Vars           ,ONLY: Species
USE MOD_Mesh_Vars               ,ONLY: NGeo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES)
INTEGER, INTENT(IN)         :: iSpec, iSF, iSide, SideID, iSample, jSample
REAL, INTENT(OUT)           :: xi(2)
REAL                        :: CalcPartPosBezier(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                        :: RandVal2(2), xiab(1:2,1:2), E, F, G, D, gradXiEta2D(1:2,1:2),gradXiEta3D(1:2,1:3), RandVal1
INTEGER                     :: iLoop
!===================================================================================================================================
iLoop=0
DO !ARM for xi considering the dA of the Subside in RefSpace
  iLoop = iLoop+1
  CALL RANDOM_NUMBER(RandVal2)
  xiab(1,1:2)=(/BezierSampleXi(iSample-1),BezierSampleXi(iSample)/) !correct order?!?
  xiab(2,1:2)=(/BezierSampleXi(JSample-1),BezierSampleXi(JSample)/) !correct order?!?
  xi=(xiab(:,2)-xiab(:,1))*RandVal2+xiab(:,1)
  IF (Species(iSpec)%Surfaceflux(iSF)%AcceptReject) THEN
    IF (.NOT.Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
      CALL EvaluateBezierPolynomialAndGradient(xi,NGeo,2 &
        ,Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample &
        ,iSide)%BezierControlPoints2D(1:2,0:NGeo,0:NGeo) &
        ,Gradient=gradXiEta2D)
      E=DOT_PRODUCT(gradXiEta2D(1,1:2),gradXiEta2D(1,1:2))
      F=DOT_PRODUCT(gradXiEta2D(1,1:2),gradXiEta2D(2,1:2))
      G=DOT_PRODUCT(gradXiEta2D(2,1:2),gradXiEta2D(2,1:2))
    ELSE
      CALL EvaluateBezierPolynomialAndGradient(xi,NGeo,3,BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID) &
        ,Gradient=gradXiEta3D)
      E=DOT_PRODUCT(gradXiEta3D(1,1:3),gradXiEta3D(1,1:3))
      F=DOT_PRODUCT(gradXiEta3D(1,1:3),gradXiEta3D(2,1:3))
      G=DOT_PRODUCT(gradXiEta3D(2,1:3),gradXiEta3D(2,1:3))
    END IF !.NOT.VeloIsNormal
    D=SQRT(E*G-F*F)
    D=D/Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Dmax !scaled Jacobian of xi
    IF (D .GT. 1.01) THEN !arbitrary warning threshold
      IPWRITE(*,'(I4,1X,A28,I0,A9,I0,A22,I0)') &
        'WARNING: ARM of SurfaceFlux ',iSF,' of Spec ',iSpec,' has inaccurate Dmax! ',D
    END IF
    CALL RANDOM_NUMBER(RandVal1)
    IF (RandVal1.LE.D) THEN
      EXIT !accept xi
    ELSE
      IF (MOD(iLoop,100).EQ.0) THEN !arbitrary warning threshold
        IPWRITE(*,'(I4,1X,A28,I0,A9,I0,A18,I0)') &
          'WARNING: ARM of SurfaceFlux ',iSF,' of Spec ',iSpec,' has reached loop ',iLoop
        IPWRITE(*,'(I4,1X,A19,2(1X,E16.8))') &
          '         R, D/Dmax:',RandVal1,D
      END IF
    END IF
  ELSE !no ARM -> accept xi
    EXIT
  END IF
END DO !Jacobian-based ARM-loop
IF(MINVAL(XI).LT.-1.)THEN
  IPWRITE(UNIT_StdOut,'(I0,A,E16.8)') ' Xi<-1',XI
END IF
IF(MAXVAL(XI).GT.1.)THEN
  IPWRITE(UNIT_StdOut,'(I0,A,E16.8)') ' Xi>1',XI
END IF
CALL EvaluateBezierPolynomialAndGradient(xi,NGeo,3,BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID),Point=CalcPartPosBezier)

END FUNCTION CalcPartPosBezier


!===================================================================================================================================
!> Calculate a random particle position for the case of regular surface flux with TriaTracking
!===================================================================================================================================
FUNCTION CalcPartPosTriaSurface(xyzNod, Vector1, Vector2, ndist, midpoint)
! MODULES
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)            :: xyzNod(3), Vector1(3), Vector2(3), ndist(3), midpoint(3)
REAL                        :: CalcPartPosTriaSurface(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                        :: RandVal2(2), PartDistance
REAL, PARAMETER             :: eps_nontria=1.0E-6
!===================================================================================================================================
  CALL RANDOM_NUMBER(RandVal2)
  IF (TrackingMethod.NE.TRIATRACKING) THEN !prevent inconsistency with non-triatracking by bilinear-routine (tol. might be increased)
    RandVal2 = RandVal2 + eps_nontria*(1. - 2.*RandVal2) !shift randVal off from 0 and 1
    DO WHILE (ABS(RandVal2(1)+RandVal2(2)-1.0).LT.eps_nontria) !sum must not be 1, since this corresponds to third egde
      CALL RANDOM_NUMBER(RandVal2)
      RandVal2 = RandVal2 + eps_nontria*(1. - 2.*RandVal2)
    END DO
  END IF
  CalcPartPosTriaSurface = xyzNod + Vector1 * RandVal2(1)
  CalcPartPosTriaSurface = CalcPartPosTriaSurface + Vector2 * RandVal2(2)
  PartDistance = ndist(1)*(CalcPartPosTriaSurface(1)-midpoint(1)) & !Distance from v1-v2
              + ndist(2)*(CalcPartPosTriaSurface(2)-midpoint(2)) &
              + ndist(3)*(CalcPartPosTriaSurface(3)-midpoint(3))
  IF (PartDistance.GT.0.) THEN !flip into right triangle if outside
    CalcPartPosTriaSurface(1:3) = 2.*midpoint(1:3)-CalcPartPosTriaSurface(1:3)
  END IF

END FUNCTION CalcPartPosTriaSurface


!===================================================================================================================================
!> Calculate a random particle position for the case of radial weighting (2D axisymmetric)
!===================================================================================================================================
SUBROUTINE CalcPartPosAxisym(iSpec,iSF,iSide,minPos,RVec,PartInsSubSide,PartInsSideRadWeight,particle_positions,allowedRejections)
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_Globals
USE MOD_Particle_Vars           ,ONLY: Species
USE MOD_DSMC_Vars               ,ONLY: RadialWeighting
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)         :: iSpec, iSF, iSide
INTEGER, INTENT(IN)         :: PartInsSubSide, PartInsSideRadWeight(:)
REAL, INTENT(IN)            :: minPos(2), RVec(2)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(OUT)        :: allowedRejections
REAL, INTENT(OUT)           :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                        :: RandVal1, PminTemp, PmaxTemp, Particle_pos(2)
INTEGER                     :: iSub, iPart, iPartSub, allowedRejectionsSub
!===================================================================================================================================
iPart=1
allowedRejections = 0
IF (RadialWeighting%DoRadialWeighting.AND.(.NOT.(ALMOSTEQUAL(minPos(2),minPos(2)+RVec(2))))) THEN
  IF(RadialWeighting%CellLocalWeighting) THEN
    DO WHILE (iPart+allowedRejections.LE.PartInsSubSide)
      CALL RANDOM_NUMBER(RandVal1)
      Particle_pos(2) = minPos(2) + RandVal1 * RVec(2)
      ! x-position depending on the y-location
      Particle_pos(1) = minPos(1) + (Particle_pos(2)-minPos(2)) * RVec(1) / RVec(2)
      particle_positions(iPart*3-2) = Particle_pos(1)
      particle_positions(iPart*3-1) = Particle_pos(2)
      particle_positions(iPart*3  ) = 0.
      IF (Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
        IF (.NOT.InSideCircularInflow(iSpec, iSF, iSide, (/Particle_pos(1),Particle_pos(2),0.0/))) THEN
          allowedRejections = allowedRejections + 1
          CYCLE
        END IF
      END IF ! CircularInflow
      iPart = iPart + 1
    END DO
  ELSE
    DO iSub = 1, RadialWeighting%nSubSides
      iPartSub = 1
      allowedRejectionsSub = 0
      DO WHILE (iPartSub+allowedRejectionsSub.LE.PartInsSideRadWeight(iSub))
        CALL RANDOM_NUMBER(RandVal1)
        PminTemp = minPos(2) + RVec(2)/RadialWeighting%nSubSides*(iSub-1.)
        PmaxTemp = minPos(2) + RVec(2)/RadialWeighting%nSubSides*iSub
        Particle_pos(2) = PminTemp + RandVal1 * (PmaxTemp - PminTemp)
        ! x-position depending on the y-location
        Particle_pos(1) = minPos(1) + (Particle_pos(2)-minPos(2)) * RVec(1) / RVec(2)
        particle_positions(iPart*3-2) = Particle_pos(1)
        particle_positions(iPart*3-1) = Particle_pos(2)
        particle_positions(iPart*3  ) = 0.
        IF (Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
          IF (.NOT.InSideCircularInflow(iSpec, iSF, iSide, (/Particle_pos(1),Particle_pos(2),0.0/))) THEN
            allowedRejectionsSub = allowedRejectionsSub + 1
            CYCLE
          END IF
        END IF ! CircularInflow
        iPart = iPart + 1
        iPartSub = iPartSub + 1
      END DO
      allowedRejections = allowedRejections + allowedRejectionsSub
    END DO
  END IF
ELSE
  DO WHILE (iPart+allowedRejections.LE.PartInsSubSide)
    CALL RANDOM_NUMBER(RandVal1)
    IF (ALMOSTEQUAL(minPos(2),minPos(2)+RVec(2))) THEN
      ! y_min = y_max, faces parallel to x-direction, constant distribution
      Particle_pos(1:2) = minPos(1:2) + RVec(1:2) * RandVal1
    ELSE
    ! No RadialWeighting, regular linear distribution of particle positions
      Particle_pos(1:2) = minPos(1:2) + RVec(1:2) &
          * ( SQRT(RandVal1*((minPos(2) + RVec(2))**2-minPos(2)**2)+minPos(2)**2) - minPos(2) ) / (RVec(2))
    END IF
    particle_positions(iPart*3-2) = Particle_pos(1)
    particle_positions(iPart*3-1) = Particle_pos(2)
    particle_positions(iPart*3  ) = 0.
    IF (Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
      IF (.NOT.InSideCircularInflow(iSpec, iSF, iSide, (/Particle_pos(1),Particle_pos(2),0.0/))) THEN
        allowedRejections = allowedRejections + 1
        CYCLE
      END IF
    END IF ! CircularInflow
    iPart = iPart + 1
  END DO
END IF

END SUBROUTINE CalcPartPosAxisym


!===================================================================================================================================
!> Calculate the particle number per side for the case of radial weighting (2D axisymmetric)
!===================================================================================================================================
SUBROUTINE CalcPartInsRadWeight(iSpec, iSF, iSample, jSample, iSide, minPos, RVec, PartInsSubSide, PartInsSideRadWeight)
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_Globals
USE MOD_TimeDisc_Vars           ,ONLY: dt,RKdtFrac
USE MOD_Particle_Vars           ,ONLY: Species, VarTimeStep
USE MOD_DSMC_Vars               ,ONLY: RadialWeighting
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)         :: iSpec, iSF, iSample, jSample, iSide
REAL, INTENT(IN)            :: minPos(2), RVec(2)
INTEGER, INTENT(OUT)        :: PartInsSubSide, PartInsSideRadWeight(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                        :: RandVal1, dtVar
INTEGER                     :: iSub
!===================================================================================================================================

! Species-specific time step
IF(VarTimeStep%UseSpeciesSpecific) THEN
  dtVar = dt * Species(iSpec)%TimeStepFactor
ELSE
  dtVar = dt
END IF

CALL RANDOM_NUMBER(RandVal1)
PartInsSubSide = INT(Species(iSpec)%Surfaceflux(iSF)%PartDensity / Species(iSpec)%MacroParticleFactor &
  * dtVar*RKdtFrac * Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%nVFR + RandVal1)
IF(.NOT.RadialWeighting%CellLocalWeighting) THEN
  IF(.NOT.ALMOSTEQUAL(minPos(2),minPos(2)+RVec(2))) THEN
    PartInsSubSide = 0
    DO iSub = 1, RadialWeighting%nSubSides
      CALL RANDOM_NUMBER(RandVal1)
      PartInsSideRadWeight(iSub) = INT(Species(iSpec)%Surfaceflux(iSF)%PartDensity / Species(iSpec)%MacroParticleFactor &
              * dtVar*RKdtFrac * Species(iSpec)%Surfaceflux(iSF)%nVFRSub(iSide,iSub)+ RandVal1)
      PartInsSubSide = PartInsSubSide + PartInsSideRadWeight(iSub)
    END DO
  END IF
END IF

END SUBROUTINE CalcPartInsRadWeight


!===================================================================================================================================
!> Calculate the particle number per side for the case of a Poisson distribution
!===================================================================================================================================
SUBROUTINE CalcPartInsPoissonDistr(iSpec, iSF, iSample, jSample, iSide, PartInsSubSide)
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_Globals
USE MOD_TimeDisc_Vars           ,ONLY: dt,RKdtFrac
USE MOD_Part_Emission_Tools     ,ONLY: SamplePoissonDistri
USE MOD_Particle_Vars           ,ONLY: Species, VarTimeStep
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)         :: iSpec, iSF, iSample, jSample, iSide
INTEGER, INTENT(OUT)        :: PartInsSubSide
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                        :: PartIns, dtVar
!===================================================================================================================================

! Species-specific time step
IF(VarTimeStep%UseSpeciesSpecific) THEN
  dtVar = dt * Species(iSpec)%TimeStepFactor
ELSE
  dtVar = dt
END IF

PartIns = Species(iSpec)%Surfaceflux(iSF)%PartDensity / Species(iSpec)%MacroParticleFactor &
                      * dtVar*RKdtFrac * Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%nVFR
IF (EXP(-PartIns).LE.TINY(PartIns)) THEN
  CALL abort(__STAMP__,'ERROR in ParticleSurfaceflux: flux is too large for poisson sampling!')
ELSE !poisson-sampling instead of random rounding (reduces numerical non-equlibrium effects [Tysanner and Garcia 2004]
  CALL SamplePoissonDistri( PartIns , PartInsSubSide )
END IF

END SUBROUTINE CalcPartInsPoissonDistr


!===================================================================================================================================
!> Calculate the particle number per side for the case of adaptive surface flux BCs
!===================================================================================================================================
SUBROUTINE CalcPartInsAdaptive(iSpec, iSF, BCSideID, iSide, iSample, jSample, minPos, RVec, PartInsSubSide, PartInsSideRadWeight)
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: BoltzmannConst, Pi
USE MOD_TimeDisc_Vars           ,ONLY: dt,RKdtFrac
USE MOD_Particle_Vars           ,ONLY: Species, VarTimeStep
USE MOD_Particle_Sampling_Vars  ,ONLY: AdaptBCMacroVal, AdaptBCMapElemToSample, AdaptBCBackupVelocity, AdaptBCPartNumOut
USE MOD_Particle_Surfaces_Vars  ,ONLY: SurfMeshSubSideData
USE MOD_Mesh_Vars               ,ONLY: SideToElem
USE MOD_DSMC_Vars               ,ONLY: RadialWeighting
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)         :: iSpec, iSF, BCSideID, iSide, jSample, iSample
REAL, INTENT(IN)            :: minPos(2), RVec(2)
INTEGER, INTENT(OUT)        :: PartInsSubSide, PartInsSideRadWeight(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                        :: ElemPartDensity, T, pressure, VeloVec(3), vec_nIn(3), veloNormal, VeloIC, VeloVecIC(3)
REAL                        :: projFak, a, v_thermal, vSF, nVFR, RandVal1, dtVar
INTEGER                     :: iSub, ElemID, SampleElemID
!===================================================================================================================================

! Species-specific time step
IF(VarTimeStep%UseSpeciesSpecific) THEN
  dtVar = dt * RKdtFrac * Species(iSpec)%TimeStepFactor
ELSE
  dtVar = dt * RKdtFrac
END IF

! 1) Calculate the velocity, density and/or temperature based on the selected BC type
ElemID = SideToElem(1,BCSideID)
SampleElemID = AdaptBCMapElemToSample(ElemID)
SELECT CASE(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType)
CASE(1) ! Pressure inlet (pressure, temperature const)
  ElemPartDensity = Species(iSpec)%Surfaceflux(iSF)%PartDensity
  T =  Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC
CASE(2) ! adaptive Outlet/freestream
  ElemPartDensity = AdaptBCMacroVal(4,SampleElemID,iSpec)
  pressure = Species(iSpec)%Surfaceflux(iSF)%AdaptivePressure
  T = pressure / (BoltzmannConst * ElemPartDensity)
CASE(3) ! Mass flow, temperature constant
  VeloVec(1:3) = AdaptBCMacroVal(1:3,SampleElemID,iSpec)
  vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
  veloNormal = VeloVec(1)*vec_nIn(1) + VeloVec(2)*vec_nIn(2) + VeloVec(3)*vec_nIn(3)
  IF(veloNormal.GT.0.0) THEN
    ElemPartDensity = Species(iSpec)%Surfaceflux(iSF)%AdaptiveMassflow &
                      / (veloNormal * Species(iSpec)%Surfaceflux(iSF)%totalAreaSF * Species(iSpec)%MassIC)
    AdaptBCBackupVelocity(1:3,SampleElemID,iSpec) = VeloVec(1:3)
  ELSE
    ! Using the old velocity vector, overwriting the sampled value with the old one
    AdaptBCMacroVal(1:3,SampleElemID,iSpec) = AdaptBCBackupVelocity(1:3,SampleElemID,iSpec)
    VeloVec(1:3) = AdaptBCMacroVal(1:3,SampleElemID,iSpec)
    vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
    veloNormal = VeloVec(1)*vec_nIn(1) + VeloVec(2)*vec_nIn(2) + VeloVec(3)*vec_nIn(3)
    IF(veloNormal.GT.0.0) THEN
      ElemPartDensity = Species(iSpec)%Surfaceflux(iSF)%AdaptiveMassflow &
        / (veloNormal * Species(iSpec)%Surfaceflux(iSF)%totalAreaSF * Species(iSpec)%MassIC)
    ELSE
      SWRITE(*,*) 'WARNING: No particles inserted!'
      SWRITE(*,*) 'WARNING: Possibly different adaptive BCs of Type3/4 have been defined next to each other.'
      SWRITE(*,*) 'WARNING: Adaptive BCs sharing a mesh element is currently not supported -> wrong velocity vector!'
      ElemPartDensity = 0
    END IF
  END IF
  T =  Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC
CASE(4) !Const. massflow inlet after Lei 2017
  T =  Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC
CASE DEFAULT
  SWRITE(*,*) 'Selected adaptive boundary condition type: ', Species(iSpec)%Surfaceflux(iSF)%AdaptiveType
  CALL abort(__STAMP__,'ERROR Adaptive Inlet: Wrong adaptive type for Surfaceflux!')
END SELECT
VeloVec(1:3) = AdaptBCMacroVal(1:3,SampleElemID,iSpec)
VeloIC = SQRT(DOT_PRODUCT(VeloVec,VeloVec))
IF (ABS(VeloIC).GT.0.) THEN
  VeloVecIC = VeloVec / VeloIC
ELSE
  VeloVecIC = (/1.,0.,0./)
END IF
vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
projFak = DOT_PRODUCT(vec_nIn,VeloVecIC) !VeloVecIC projected to inwards normal
v_thermal = SQRT(2.*BoltzmannConst*T/Species(iSpec)%MassIC) !thermal speed
a = 0 !dummy for projected speed ratio in constant v-distri

! 2) Compute total volume flow rate through surface depending on the area and velocity/temperature
SELECT CASE(TRIM(Species(iSpec)%Surfaceflux(iSF)%velocityDistribution))
CASE('constant')
  vSF = VeloIC * projFak !Velo proj. to inwards normal
  vSF = MAX(vSF,0.) !VFR proj. to inwards normal (only positive parts!)
CASE('maxwell','maxwell_lpn')
  IF ( ALMOSTEQUAL(v_thermal,0.)) THEN
    v_thermal = 1.
  END IF
  a = VeloIC * projFak / v_thermal !speed ratio proj. to inwards n (can be negative!)
  vSF = v_thermal / (2.0*SQRT(PI)) * ( EXP(-(a*a)) + a*SQRT(PI)*(1+ERF(a)) ) !mean flux velocity through normal sub-face
CASE DEFAULT
  CALL abort(__STAMP__,'ERROR in CalcPartInsAdaptive: Wrong velocity distribution!')
END SELECT

IF(RadialWeighting%DoRadialWeighting) THEN
  ! In case of adaptive SF, nVFR is initialized as only a weighted area
  nVFR = Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%nVFR * vSF
ELSE
  nVFR = SurfMeshSubSideData(iSample,jSample,BCSideID)%area * vSF
END IF

! 3) Calculate the actual number of particles per side
IF(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType.EQ.4) THEN
  ! TODO: RADIAL WEIGHTING TREATMENT
  CALL RANDOM_NUMBER(RandVal1)
  PartInsSubSide = INT(Species(iSpec)%Surfaceflux(iSF)%ConstMassflowWeight(iSample,jSample,iSide)     &
                          * (Species(iSpec)%Surfaceflux(iSF)%AdaptiveMassflow * dtVar    &
                              / (Species(iSpec)%MassIC * Species(iSpec)%MacroParticleFactor)  &
                              + REAL(AdaptBCPartNumOut(iSpec,iSF))) +RandVal1)
ELSE
  CALL RANDOM_NUMBER(RandVal1)
  PartInsSubSide = INT(ElemPartDensity / Species(iSpec)%MacroParticleFactor * dtVar * nVFR + RandVal1)
  ! Radial weighting: subdivide the side into smaller subsides to improve distribution
  IF(RadialWeighting%DoRadialWeighting.AND..NOT.RadialWeighting%CellLocalWeighting) THEN
    ! Skip sides parallel to rotational axis
    IF(.NOT.ALMOSTEQUAL(minPos(2),minPos(2)+RVec(2))) THEN
      PartInsSubSide = 0
      DO iSub = 1, RadialWeighting%nSubSides
        CALL RANDOM_NUMBER(RandVal1)
        PartInsSideRadWeight(iSub) = INT(ElemPartDensity / Species(iSpec)%MacroParticleFactor &
                * dtVar * Species(iSpec)%Surfaceflux(iSF)%nVFRSub(iSide,iSub) * vSF + RandVal1)
        PartInsSubSide = PartInsSubSide + PartInsSideRadWeight(iSub)
      END DO
    END IF
  END IF
END IF

! DEBUG OUTPUT
IF (PartInsSubSide.LT.0) THEN
  IPWRITE(*,*) 'ERROR in CalcPartInsAdaptive: Calculated number of particles to insert below zero! PartInsSubSide: ', PartInsSubSide
  IPWRITE(*,*) 'Species Index: ', iSpec, 'SurfaceFlux Index: ', iSF, 'iSide: ', iSide
  IPWRITE(*,*) 'Utilized values for adaptive surface flux:'
  IPWRITE(*,*) 'VeloVec(1:3)', VeloVec(1:3), 'Temperature', T
  IPWRITE(*,*) 'MassIC: ', Species(iSpec)%MassIC, 'MPF: ', Species(iSpec)%MacroParticleFactor, 'dt', dtVar
  SELECT CASE(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType)
  CASE(1)
    IPWRITE(*,*) 'ElemPartDensity', ElemPartDensity
  CASE(2)
    IPWRITE(*,*) 'ElemPartDensity', ElemPartDensity, 'pressure', pressure
  CASE(3)
    IPWRITE(*,*) 'ElemPartDensity', ElemPartDensity
    IPWRITE(*,*) 'AdaptBCBackupVelocity(1:3,SampleElemID,iSpec)', AdaptBCBackupVelocity(1:3,SampleElemID,iSpec)
  CASE(4)
    IPWRITE(*,*) 'ConstMassflowWeight(iSample,jSample,iSide)', Species(iSpec)%Surfaceflux(iSF)%ConstMassflowWeight(iSample,jSample,iSide)
    IPWRITE(*,*) 'AdaptBCPartNumOut(iSpec,iSF)', AdaptBCPartNumOut(iSpec,iSF)
    IPWRITE(*,*) 'AdaptiveMassflow', Species(iSpec)%Surfaceflux(iSF)%AdaptiveMassflow
  END SELECT
  CALL abort(__STAMP__,'ERROR in CalcPartInsAdaptive: PartInsSubSide.LT.0!')
END IF

END SUBROUTINE CalcPartInsAdaptive


!===================================================================================================================================
!> Routine calculates the weights of the triangles for AdaptiveType=4 depending on the side-specific volume flow rate
!===================================================================================================================================
SUBROUTINE CalcConstMassflowWeight(iSpec,iSF)
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst, Pi
USE MOD_Particle_Vars          ,ONLY: Species
USE MOD_Particle_Sampling_Vars ,ONLY: AdaptBCMacroVal, AdaptBCMapElemToSample, AdaptBCBackupVelocity
USE MOD_Particle_Surfaces_Vars ,ONLY: SurfMeshSubSideData, BCdata_auxSF, SurfFluxSideSize
USE MOD_Mesh_Vars              ,ONLY: SideToElem, offsetElem
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: iSpec, iSF
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSide, BCSideID, ElemID, SideID, currentBC, iSample, jSample, SampleElemID
REAL                            :: VeloVec(1:3), vec_nIn(1:3), nVFRTotal, VeloIC, VeloVecIC(1:3), projFak
REAL                            :: v_thermal, a, vSF, nVFR, area
!===================================================================================================================================

ASSOCIATE(SF => Species(iSpec)%Surfaceflux(iSF))

currentBC = SF%BC

nVFRTotal = 0.
SF%ConstMassflowWeight = 0.

DO iSide=1,BCdata_auxSF(currentBC)%SideNumber
  ! Skip sides outside of the circular inflow region
  IF (SF%CircularInflow) THEN
    IF(SF%SurfFluxSideRejectType(iSide).EQ.1) CYCLE
  END IF
  BCSideID=BCdata_auxSF(currentBC)%SideList(iSide)
  ElemID = SideToElem(S2E_ELEM_ID,BCSideID)
  SampleElemID = AdaptBCMapElemToSample(ElemID)
  SideID=GetGlobalNonUniqueSideID(offsetElem+ElemID,SideToElem(S2E_LOC_SIDE_ID,BCSideID))
  ! Get the sampled velocity vector
  VeloVec(1:3) = AdaptBCMacroVal(1:3,SampleElemID,iSpec)
  ! Determine the velocity magnitude
  VeloIC = SQRT(DOT_PRODUCT(VeloVec,VeloVec))
  IF (ABS(VeloIC).GT.0.) THEN
    ! Calculate the normalized velocity vector
    VeloVecIC = VeloVec / VeloIC
    ! Store the vector as backup for low particle numbers
    AdaptBCBackupVelocity(1:3,SampleElemID,iSpec) = VeloVec(1:3)
  ELSE
    ! Using the old velocity vector, overwriting the sampled value with the old one
    VeloVec(1:3) = AdaptBCBackupVelocity(1:3,SampleElemID,iSpec)
    AdaptBCMacroVal(1:3,SampleElemID,iSpec) = AdaptBCBackupVelocity(1:3,SampleElemID,iSpec)
    VeloIC = SQRT(DOT_PRODUCT(VeloVec,VeloVec))
    IF(ABS(VeloIC).GT.0.) THEN
      VeloVecIC = VeloVec / VeloIC
    ELSE
      ! Dummy value, for maxwell only the thermal velocity will be considered
      VeloVecIC = (/1.,0.,0./)
    END IF
  END IF
  ! Loop over the triangles
  DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
    ! Set the area of the side, different area for circular inflow
    IF(SF%CircularInflow) THEN
      area = SF%CircleAreaPerTriaSide(iSample,jSample,iSide)
    ELSE
      area = SurfMeshSubSideData(iSample,jSample,BCSideID)%area
    END IF
    ! VeloVecIC projected to inwards normal
    vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
    projFak = DOT_PRODUCT(vec_nIn,VeloVecIC)
    ! Compute total volume flow rate through surface
    SELECT CASE(TRIM(SF%velocityDistribution))
    CASE('constant')
      vSF = VeloIC * projFak !Velo proj. to inwards normal
      nVFR = MAX(area * vSF,0.) !VFR proj. to inwards normal (only positive parts!)
    CASE('maxwell','maxwell_lpn')
      ! Thermal velocity
      v_thermal = SQRT(2.*BoltzmannConst*SF%MWTemperatureIC/Species(iSpec)%MassIC)
      IF ( ALMOSTEQUAL(v_thermal,0.)) THEN
        v_thermal = 1.
      END IF
      a = VeloIC * projFak / v_thermal !speed ratio proj. to inwards n (can be negative!)
      vSF = v_thermal / (2.0*SQRT(PI)) * ( EXP(-(a*a)) + a*SQRT(PI)*(1+ERF(a)) ) !mean flux velocity through normal sub-face
      nVFR = area * vSF !VFR projected to inwards normal of sub-side
    CASE DEFAULT
      CALL abort(__STAMP__,'ERROR in CalcConstMassflowWeight: Wrong velocity distribution!')
    END SELECT
    ! Skip the side if a negative/zero volume flow rate has been determined
    IF(nVFR.LE.0.0) CYCLE
    ! Calculate the volume flow rate per side
    SF%ConstMassflowWeight(iSample,jSample,iSide) = nVFR
    ! Calculate the total volume flow rate
    nVFRTotal = nVFRTotal + nVFR
  END DO; END DO
END DO

! Calculate the total volume flow rate over the whole BC across all processors
#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,nVFRTotal,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_PICLAS,IERROR)
#endif

! Determine the weight of each side compared to the total volume flow rate
IF(nVFRTotal.GT.0.) THEN
  SF%ConstMassflowWeight(:,:,:) = SF%ConstMassflowWeight(:,:,:) / REAL(nVFRTotal)
ELSE
  SF%ConstMassflowWeight(:,:,:) = 0.
END IF

IF(SF%CircularInflow) THEN
  ! Scaling up the number of particles to be inserted on the triaside
  DO iSide=1,BCdata_auxSF(currentBC)%SideNumber
    BCSideID=BCdata_auxSF(currentBC)%SideList(iSide)
    DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
      IF(SF%CircleAreaPerTriaSide(iSample,jSample,iSide).GT.0.0) THEN
        SF%ConstMassflowWeight(iSample,jSample,iSide) = SF%ConstMassflowWeight(iSample,jSample,iSide) &
            * SurfMeshSubSideData(iSample,jSample,BCSideID)%area / SF%CircleAreaPerTriaSide(iSample,jSample,iSide)
      ELSE
        SF%ConstMassflowWeight(iSample,jSample,iSide) = 0.0
      END IF
    END DO; END DO
  END DO
END IF

END ASSOCIATE

END SUBROUTINE CalcConstMassflowWeight


!===================================================================================================================================
!> Determine the particle velocity of each inserted particle
!===================================================================================================================================
SUBROUTINE SetSurfacefluxVelocities(iSpec,iSF,iSample,jSample,iSide,BCSideID,SideID,ElemID,NbrOfParticle,PartIns)
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: PI, BoltzmannConst
USE MOD_Particle_Vars
USE MOD_Particle_Surfaces_Vars  ,ONLY: SurfMeshSubSideData, TriaSurfaceFlux
USE MOD_Particle_Surfaces       ,ONLY: CalcNormAndTangBezier
USE MOD_Particle_Sampling_Vars  ,ONLY: AdaptBCMapElemToSample, AdaptBCMacroVal
USE MOD_Part_Tools              ,ONLY: InRotRefFrameCheck, GetNextFreePosition
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: iSpec,iSF,iSample,jSample,iSide,BCSideID,SideID,ElemID,NbrOfParticle,PartIns
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: i,PositionNbr,envelope,currentBC,SampleElemID
REAL                             :: Vec3D(3), vec_nIn(1:3), vec_t1(1:3), vec_t2(1:3)
REAL                             :: a,zstar,RandVal1,RandVal2(2),RandVal3(3),u,RandN,RandN_save,Velo1,Velo2,Velosq,T,beta,z
LOGICAL                          :: RandN_in_Mem
REAL                             :: projFak                          ! VeloVecIC projected to inwards normal of tria
REAL                             :: Velo_t1                          ! Velo comp. of first orth. vector in tria
REAL                             :: Velo_t2                          ! Velo comp. of second orth. vector in tria
REAL                             :: VeloIC
REAL                             :: VeloVec(1:3)
REAL                             :: VeloVecIC(1:3),v_thermal, pressure
!===================================================================================================================================

IF(PartIns.LT.1) RETURN

RandN_in_Mem=.FALSE.
envelope=-1
currentBC = Species(iSpec)%Surfaceflux(iSF)%BC

IF (.NOT.Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
  vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
  vec_t1(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t1(1:3)
  vec_t2(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t2(1:3)
END IF !.NOT.VeloIsNormal

IF(.NOT.Species(iSpec)%Surfaceflux(iSF)%Adaptive) THEN
  VeloIC = Species(iSpec)%Surfaceflux(iSF)%VeloIC
  T = Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC
  a = Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%a_nIn
  projFak = Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%projFak
  Velo_t1 = Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t1
  Velo_t2 = Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t2
ELSE !Species(iSpec)%Surfaceflux(iSF)%Adaptive
  SampleElemID = AdaptBCMapElemToSample(ElemID)
  SELECT CASE(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType)
  CASE(1,3,4) ! Pressure and massflow inlet (pressure/massflow, temperature const)
    T =  Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC
  CASE(2) ! adaptive Outlet/freestream
    pressure = Species(iSpec)%Surfaceflux(iSF)%AdaptivePressure
    T = pressure / (BoltzmannConst * AdaptBCMacroVal(4,SampleElemID,iSpec))
  CASE DEFAULT
    CALL abort(__STAMP__,'ERROR in SurfaceFlux: Wrong adaptive type for Surfaceflux velocities!')
  END SELECT
  VeloVec(1) = AdaptBCMacroVal(DSMC_VELOX,SampleElemID,iSpec)
  VeloVec(2) = AdaptBCMacroVal(DSMC_VELOY,SampleElemID,iSpec)
  VeloVec(3) = AdaptBCMacroVal(DSMC_VELOZ,SampleElemID,iSpec)
  vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
  VeloVec(1:3) = DOT_PRODUCT(VeloVec,vec_nIn)*vec_nIn(1:3)
  VeloIC = SQRT(DOT_PRODUCT(VeloVec,VeloVec))
  IF (ABS(VeloIC).GT.0.) THEN
    VeloVecIC = VeloVec / VeloIC
  ELSE
    VeloVecIC = (/1.,0.,0./)
  END IF
  projFak = DOT_PRODUCT(vec_nIn,VeloVecIC) !VeloVecIC projected to inwards normal
  v_thermal = SQRT(2.*BoltzmannConst*T/Species(iSpec)%MassIC) !thermal speed
  IF ( ALMOSTEQUAL(v_thermal,0.)) THEN
    v_thermal = 1.
  END IF
  a = VeloIC * projFak / v_thermal !speed ratio proj. to inwards n (can be negative!)
  Velo_t1 = VeloIC * DOT_PRODUCT(vec_t1,VeloVecIC) !v in t1-dir
  Velo_t2 = VeloIC * DOT_PRODUCT(vec_t2,VeloVecIC) !v in t2-dir
END IF !Adaptive SurfaceFlux

! Set velocities
SELECT CASE(TRIM(Species(iSpec)%Surfaceflux(iSF)%velocityDistribution))
CASE('constant')
  IF (.NOT.Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
    VeloVecIC(1:3) = Species(iSpec)%Surfaceflux(iSF)%VeloVecIC(1:3)
    VeloVecIC(1:3) = VeloVecIC(1:3) / VECNORM(VeloVecIC(1:3))
  END IF
  DO i = NbrOfParticle-PartIns+1,NbrOfParticle
    PositionNbr = GetNextFreePosition(i)
    ! In case of side-normal velocities: calc n-vector at particle position, xi was saved in PartState(4:5)
    IF (Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal .AND. TriaSurfaceFlux) THEN
      vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
      vec_t1(1:3) = 0. !dummy
      vec_t2(1:3) = 0. !dummy
    ELSE IF (Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
      CALL CalcNormAndTangBezier( nVec=vec_nIn(1:3),xi=PartState(4,PositionNbr),eta=PartState(5,PositionNbr),SideID=SideID )
      vec_nIn(1:3) = -vec_nIn(1:3)
      vec_t1(1:3) = 0. !dummy
      vec_t2(1:3) = 0. !dummy
    ELSE
      vec_nIn(1:3) = VeloVecIC(1:3)
    END IF !VeloIsNormal
    ! Build complete velo-vector
    Vec3D(1:3) = vec_nIn(1:3) * Species(iSpec)%Surfaceflux(iSF)%VeloIC
    PartState(4:6,PositionNbr) = Vec3D(1:3)
  END DO !i = ...NbrOfParticle
CASE('maxwell','maxwell_lpn')
  !-- determine envelope for most efficient ARM [Garcia and Wagner 2006, JCP217-2]
  IF (ALMOSTZERO(VeloIC*projFak)) THEN
    ! Rayleigh distri
    envelope = 0
  ELSE IF (-0.4.LT.a .AND. a.LT.1.3) THEN
    ! low speed flow
    IF (a.LE.0.) THEN
      envelope = 1
    ELSE
      envelope = 3
    END IF !choose envelope based on flow direction
  ELSE
    ! high speed / general flow
    IF (a.LT.0.) THEN
      envelope = 2
    ELSE
      envelope = 4
    END IF !choose envelope based on flow direction
  END IF !low speed / high speed / rayleigh flow

  DO i = NbrOfParticle-PartIns+1,NbrOfParticle
    PositionNbr = GetNextFreePosition(i)
    !-- 0a.: In case of side-normal velocities: calc n-/t-vectors at particle position, xi was saved in PartState(4:5)
    IF (Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal .AND. TriaSurfaceFlux) THEN
      vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
      vec_t1(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t1(1:3)
      vec_t2(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t2(1:3)
    ELSE IF (Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
      CALL CalcNormAndTangBezier( nVec=vec_nIn(1:3),tang1=vec_t1(1:3),tang2=vec_t2(1:3) &
        ,xi=PartState(4,PositionNbr),eta=PartState(5,PositionNbr),SideID=SideID )
      vec_nIn(1:3) = -vec_nIn(1:3)
    END IF !VeloIsNormal
    !-- 1.: determine zstar (initial generation of potentially too many RVu is for needed indentities of RVu used multiple times!
    SELECT CASE(envelope)
    CASE(0)
      CALL RANDOM_NUMBER(RandVal1)
      zstar = -SQRT(-LOG(RandVal1))
    CASE(1)
      DO
        CALL RANDOM_NUMBER(RandVal2)
        zstar = -SQRT(a*a-LOG(RandVal2(1)))
        IF ( -(a-zstar)/zstar .GT. RandVal2(2)) THEN
          EXIT
        END IF
      END DO
    CASE(2)
      z = 0.5*(a-SQRT(a*a+2.))
      beta  = a-(1.0-a)*(a-z)
      DO
        CALL RANDOM_NUMBER(RandVal3)
        IF (EXP(-(beta*beta))/(EXP(-(beta*beta))+2.0*(a-z)*(a-beta)*EXP(-(z*z))).GT.RandVal3(1)) THEN
          zstar=-SQRT(beta*beta-LOG(RandVal3(2)))
          IF ( -(a-zstar)/zstar .GT. RandVal3(3)) THEN
            EXIT
          END IF
        ELSE
          zstar=beta+(a-beta)*RandVal3(2)
          IF ( (a-zstar)/(a-z)*EXP(z*z-(zstar*zstar)) .GT. RandVal3(3)) THEN
            EXIT
          END IF
        END IF
      END DO
    CASE(3)
      DO
        CALL RANDOM_NUMBER(RandVal3)
        u = RandVal3(1)
        IF ( a*SQRT(PI)/(a*SQRT(PI)+1+a*a) .GT. u) THEN
!            IF (.NOT.DoZigguratSampling) THEN !polar method
            IF (RandN_in_Mem) THEN !reusing second RandN form previous polar method
              RandN = RandN_save
              RandN_in_Mem=.FALSE.
            ELSE
              Velosq = 2
              DO WHILE ((Velosq .GE. 1.) .OR. (Velosq .EQ. 0.))
                CALL RANDOM_NUMBER(RandVal2)
                Velo1 = 2.*RandVal2(1) - 1.
                Velo2 = 2.*RandVal2(2) - 1.
                Velosq = Velo1**2 + Velo2**2
              END DO
              RandN = Velo1*SQRT(-2*LOG(Velosq)/Velosq)
              RandN_save = Velo2*SQRT(-2*LOG(Velosq)/Velosq)
              RandN_in_Mem=.TRUE.
            END IF
!            ELSE !ziggurat method
!              RandN=rnor()
!            END IF
          zstar = -1./SQRT(2.)*ABS(RandN)
          EXIT
        ELSE IF ( (a*SQRT(PI)+1.)/(a*SQRT(PI)+1+a*a) .GT. u) THEN
          zstar = -SQRT(-LOG(RandVal3(2)))
          EXIT
        ELSE
          zstar = (1.0-SQRT(RandVal3(2)))*a
          IF (EXP(-(zstar*zstar)).GT.RandVal3(3)) THEN
            EXIT
          END IF
        END IF
      END DO
    CASE(4)
      DO
        CALL RANDOM_NUMBER(RandVal3)
        IF (1.0/(2.0*a*SQRT(PI)+1.0).GT.RandVal3(1)) THEN
          zstar=-SQRT(-LOG(RandVal3(2)))
        ELSE
!            IF (.NOT.DoZigguratSampling) THEN !polar method
            IF (RandN_in_Mem) THEN !reusing second RandN form previous polar method
              RandN = RandN_save
              RandN_in_Mem=.FALSE.
            ELSE
              Velosq = 2
              DO WHILE ((Velosq .GE. 1.) .OR. (Velosq .EQ. 0.))
                CALL RANDOM_NUMBER(RandVal2)
                Velo1 = 2.*RandVal2(1) - 1.
                Velo2 = 2.*RandVal2(2) - 1.
                Velosq = Velo1**2 + Velo2**2
              END DO
              RandN = Velo1*SQRT(-2*LOG(Velosq)/Velosq)
              RandN_save = Velo2*SQRT(-2*LOG(Velosq)/Velosq)
              RandN_in_Mem=.TRUE.
            END IF
!            ELSE !ziggurat method
!              RandN=rnor()
!            END IF
          zstar = 1./SQRT(2.)*RandN
        END IF
        IF ( (a-zstar)/a .GT. RandVal3(3)) THEN
          EXIT
        END IF
      END DO
    CASE DEFAULT
      CALL abort(__STAMP__,'ERROR in SurfaceFlux: Wrong envelope in SetSurfacefluxVelocities!')
    END SELECT
    !-- 2.: sample normal directions and build complete velo-vector
    Vec3D(1:3) = vec_nIn(1:3) * SQRT(2.*BoltzmannConst*T/Species(iSpec)%MassIC)*(a-zstar)
!      IF (.NOT.DoZigguratSampling) THEN !polar method
      Velosq = 2
      DO WHILE ((Velosq .GE. 1.) .OR. (Velosq .EQ. 0.))
        CALL RANDOM_NUMBER(RandVal2)
        Velo1 = 2.*RandVal2(1) - 1.
        Velo2 = 2.*RandVal2(2) - 1.
        Velosq = Velo1**2 + Velo2**2
      END DO
      Velo1 = Velo1*SQRT(-2*LOG(Velosq)/Velosq)
      Velo2 = Velo2*SQRT(-2*LOG(Velosq)/Velosq)
!      ELSE !ziggurat method
!        Velo1=rnor()
!        Velo2=rnor()
!      END IF
    Vec3D(1:3) = Vec3D(1:3) + vec_t1(1:3) * ( Velo_t1+Velo1*SQRT(BoltzmannConst*T/Species(iSpec)%MassIC) )
    Vec3D(1:3) = Vec3D(1:3) + vec_t2(1:3) * ( Velo_t2+Velo2*SQRT(BoltzmannConst*T/Species(iSpec)%MassIC) )
    PartState(4:6,PositionNbr) = Vec3D(1:3)
  END DO !i = ...NbrOfParticle
CASE DEFAULT
  CALL abort(__STAMP__,'ERROR in SurfaceFlux: Wrong velocity distribution!')
END SELECT

IF(UseRotRefFrame) THEN
  DO i = NbrOfParticle-PartIns+1,NbrOfParticle
    PositionNbr = GetNextFreePosition(i)
    ! Detect if particle is within a RotRefDomain
    PDM%InRotRefFrame(PositionNbr) = InRotRefFrameCheck(PositionNbr)
    ! Initialize velocity in the rotational frame of reference
    IF(PDM%InRotRefFrame(PositionNbr)) THEN
      PartVeloRotRef(1:3,PositionNbr) = PartState(4:6,PositionNbr) - CROSS(RotRefFrameOmega(1:3),PartState(1:3,PositionNbr))
    END IF
  END DO
END IF

END SUBROUTINE SetSurfacefluxVelocities

#if (USE_HDG)
!===================================================================================================================================
!> Thermionic emission: Calculation of the work function reduction due to an electric field (Schottky effect)
!===================================================================================================================================
SUBROUTINE CalcPartInsThermionicEmissionSchottky(iSpec,iSF,iSample,jSample,iSide,iLocSide,ElemID,PartInsSubSide)
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_Globals
USE MOD_PreProc
! VARIABLES
USE MOD_Globals_Vars            ,ONLY: BoltzmannConst, Pi, ElementaryCharge, eps0
USE MOD_TimeDisc_Vars           ,ONLY: dt,RKdtFrac
USE MOD_Particle_Vars           ,ONLY: Species, VarTimeStep
USE MOD_Equation_Vars           ,ONLY: E
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound
! ROUTINES
USE MOD_ProlongToFace           ,ONLY: ProlongToFace_Elementlocal
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)         :: iSpec,iSF,iSample,jSample,iSide,iLocSide,ElemID
INTEGER, INTENT(OUT)        :: PartInsSubSide
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                        :: EFieldFace(1:3,0:PP_N,0:PP_N), EFaceMag, CurrentDensity, WallTemp, WorkFunction, RandVal1, dtVar
!===================================================================================================================================
ASSOCIATE(SF => Species(iSpec)%Surfaceflux(iSF))

! Species-specific time step
IF(VarTimeStep%UseSpeciesSpecific) THEN
  dtVar = dt * Species(iSpec)%TimeStepFactor
ELSE
  dtVar = dt
END IF

! 1) Determine the electric field at the surface
CALL ProlongToFace_Elementlocal(nVar=3,locSideID=iLocSide,Uvol=E(:,:,:,:,ElemID),Uface=EFieldFace)

! 2) Average the e-field vector and calculate magnitude
EFaceMag = (PP_N+1)**2
EFaceMag = VECNORM((/SUM(EFieldFace(1,:,:))/EFaceMag, SUM(EFieldFace(2,:,:))/EFaceMag, SUM(EFieldFace(3,:,:))/EFaceMag/))

! 3) Calculate the work function with the Schottky effect and the new current density [A/m2]
WallTemp = PartBound%WallTemp(SF%BC)
WorkFunction = SF%WorkFunctionTE- SQRT((EFaceMag*ElementaryCharge**3)/(4*Pi*eps0))

CurrentDensity = SF%RichardsonConstant * WallTemp**2 * EXP(-WorkFunction / (BoltzmannConst * WallTemp))
CurrentDensity = CurrentDensity / ABS(Species(iSpec)%ChargeIC)

! 4) Determine the number of particles per side
CALL RANDOM_NUMBER(RandVal1)
PartInsSubSide = INT(CurrentDensity / Species(iSpec)%MacroParticleFactor * dtVar*RKdtFrac &
                     * SF%SurfFluxSubSideData(iSample,jSample,iSide)%nVFR + RandVal1)

END ASSOCIATE

END SUBROUTINE CalcPartInsThermionicEmissionSchottky
#endif /*(USE_HDG)*/

END MODULE MOD_Particle_SurfFlux
