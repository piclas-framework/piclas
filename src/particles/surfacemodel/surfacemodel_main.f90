!==================================================================================================================================
! Copyright (c) 2015 - 2019 Wladimir Reschke
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

MODULE MOD_SurfaceModel
!===================================================================================================================================
!> Main Module for surface model
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE ReactiveSurfaceTreatment
  MODULE PROCEDURE ReactiveSurfaceTreatment
END INTERFACE

PUBLIC :: SurfaceModel_main
PUBLIC :: UpdateSurfModelVars
PUBLIC :: ReactiveSurfaceTreatment
!===================================================================================================================================

CONTAINS

SUBROUTINE SurfaceModel_main()
!===================================================================================================================================
!> Main Routine treating all surface and calculating desorbing / evaporating number of particles
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound
USE MOD_SMCR                   ,ONLY: SMCR_PartDesorb, SMCR_Diffusion
USE MOD_SurfaceModel_Tools     ,ONLY: CalcEvapPartNum
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
IF (.NOT.ANY(PartBound%Reactive)) RETURN
! evaporating particles for surface models using mean values on surfaces
CALL CalcEvapPartNum()
! desorbing particles for surface model 3
CALL SMCR_Diffusion()
CALL SMCR_PartDesorb()

END SUBROUTINE SurfaceModel_main


SUBROUTINE UpdateSurfModelVars()
!===================================================================================================================================
!> Update and sample surface values for adsorption, desorption and reactions on surfaces (heterogeneous reactions)
!> Only procs with surface enter function, all other exit routine
!> 1. Communicate number of particles that were absorbed by halo-sides of neighbour procs (SumAdsorbPart and SumERDesorbed),
!>    so that own proc has number of total adsorption particles for each species on own surfaces
!> 2. After communication go through all own sides (no mpi halo sides) and adjust coverage resulting from changes through
!>    adsorption and reactions for all surfacemodels
!> 2.1  For surfacemodel=3 calculate number of adsorbate change (surfacempf!=gasmpf) and if changed Call AdjustReconstructMapNum
!>      to adjust number of adsorbates on reconstructed Monte Carlo surface
!> 2.2  Sample macroscopic surface coverage values
!> 3. Reset/Adjust surface model sum counters
!> 4. Calculate global mean probabilities for surface models
!> 5. Send Coverages and Distribution info of own sides to halo sides of other procs (other procs have own sides as halo sides
!>    and need coverage info for adsorption calculation (mpi routines know which sides to communicate)
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Particle_Vars          ,ONLY: WriteMacroSurfaceValues, KeepWallParticles, Species, nSpecies
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, SampWall, PartBound
USE MOD_TimeDisc_Vars          ,ONLY: tend,time
#if USE_LOADBALANCE
USE MOD_LoadBalance_tools      ,ONLY: LBStartTime,LBSplitTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption, SurfDistInfo, SurfModel
USE MOD_SurfaceModel_Tools     ,ONLY: CalcAdsorbProb, CalcDesorbProb
USE MOD_SurfaceModel_Tools     ,ONLY: SMCR_AdjustMapNum, IsReactiveSurface, SurfaceHasModelNum
#ifdef MPI
USE MOD_SurfaceModel_MPI       ,ONLY: ExchangeAdsorbNum, ExchangeCoverageInfo, ExchangeSurfDistInfo
#endif /*MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                          :: iSpec, iSurfSide, p, q, new_adsorbates, numSites
REAL                             :: maxPart
REAL                             :: coverage_tmp, coverage_corrected
#if USE_LOADBALANCE
REAL                             :: tLBStart
#endif /*USE_LOADBALANCE*/
#if (PP_TimeDiscMethod==42)
REAL                             :: desorbnum_covreduce
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
IF (.NOT.SurfMesh%SurfOnProc) RETURN
IF (.NOT.ANY(Partbound%Reactive)) RETURN

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
IF (.NOT.KeepWallParticles) THEN
#ifdef MPI
!----- 1.
  CALL ExchangeAdsorbNum()
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_SURFCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*MPI*/

  ! adjust coverages of all species on surfaces
  DO iSpec = 1,nSpecies
#if (PP_TimeDiscMethod==42)
    IF (DSMC%ReservoirRateStatistic) SurfModel%Info(iSpec)%WallSpecNumCount = 0
#endif
!----- 2.
    DO iSurfSide = 1,SurfMesh%nSides
      IF (.NOT.IsReactiveSurface(iSurfSide)) CYCLE
      DO q = 1,nSurfSample
        DO p = 1,nSurfSample
#if (PP_TimeDiscMethod==42)
          ! write number of adsorbates into info array
          IF (DSMC%ReservoirRateStatistic) THEN
            maxPart = REAL(INT(Adsorption%DensSurfAtoms(iSurfSide) * SurfMesh%SurfaceArea(p,q,iSurfSide),8))
            SurfModel%Info(iSpec)%WallSpecNumCount = SurfModel%Info(iSpec)%WallSpecNumCount &
                                                          + INT( Adsorption%Coverage(p,q,iSurfSide,iSpec) &
                                                          * maxPart/Species(iSpec)%MacroParticleFactor)
          END IF
#endif
          SELECT CASE (SurfaceHasModelNum(iSurfSide))
          CASE(1)
            maxPart = Adsorption%DensSurfAtoms(iSurfSide) * SurfMesh%SurfaceArea(p,q,iSurfSide)
            Adsorption%Coverage(p,q,iSurfSide,iSpec) = Adsorption%Coverage(p,q,iSurfSide,iSpec) &
                + ( SurfModel%SumAdsorbPart(p,q,iSurfSide,iSpec) &
                - (SurfModel%SumDesorbPart(p,q,iSurfSide,iSpec) - SurfModel%SumReactPart(p,q,iSurfSide,iSpec)) ) &
                * Species(iSpec)%MacroParticleFactor / maxPart
          CASE(2)
            Adsorption%Coverage(p,q,iSurfSide,iSpec) = Adsorption%Coverage(p,q,iSurfSide,iSpec) &
                + SurfModel%SumAdsorbPart(p,q,iSurfSide,iSpec)
          CASE(3)
#if (PP_TimeDiscMethod==42)
            IF (Adsorption%CoverageReduction) THEN
              ! calculate number of desorbed particles for each species
              desorbnum_covreduce = (REAL(Adsorption%CovReductionStep(iSpec))/REAL(SurfDistInfo(p,q,iSurfSide)%nSites(3))) &
                  * REAL(INT(Adsorption%DensSurfAtoms(iSurfSide) &
                  * SurfMesh%SurfaceArea(p,q,iSurfSide),8)) / Species(iSpec)%MacroParticleFactor
              coverage_tmp = Adsorption%Coverage(p,q,iSurfSide,iSpec) &
                  - REAL(desorbnum_covreduce) * Species(iSpec)%MacroParticleFactor &
                  / REAL(INT(Adsorption%DensSurfAtoms(iSurfSide) * SurfMesh%SurfaceArea(p,q,iSurfSide),8))
              IF(coverage_tmp.GT.0.) THEN
                Adsorption%Coverage(p,q,iSurfSide,iSpec) = coverage_tmp
                new_adsorbates = -1*Adsorption%CovReductionStep(iSpec)
                CALL SMCR_AdjustMapNum(p,q,iSurfSide,new_adsorbates,iSpec)
              END IF
            END IF
#endif
            maxPart = REAL(INT(Adsorption%DensSurfAtoms(iSurfSide) * SurfMesh%SurfaceArea(p,q,iSurfSide),8))
            coverage_tmp = Adsorption%Coverage(p,q,iSurfSide,iSpec) &
                + REAL(( SurfModel%SumAdsorbPart(p,q,iSurfSide,iSpec) &
                - (SurfModel%SumDesorbPart(p,q,iSurfSide,iSpec) - SurfModel%SumReactPart(p,q,iSurfSide,iSpec)) ) &
                * Species(iSpec)%MacroParticleFactor) / maxPart
            IF (coverage_tmp.LT.0.) THEN ! can only happen for (ER + desorption) or (ER + MPF_surf<MPF_gas) --> SumAdsorbPart<0
              coverage_corrected = 0. !Adsorption%Coverage(p,q,iSurfSide,iSpec) &
                  !+ REAL(( - (SurfModel%SumDesorbPart(p,q,iSurfSide,iSpec) - SurfModel%SumReactPart(p,q,iSurfSide,iSpec)) ) &
                  !* Species(iSpec)%MacroParticleFactor) / maxPart
            ELSE
              coverage_corrected = coverage_tmp
            END IF
!----- 2.1
            IF (SurfModel%SumAdsorbPart(p,q,iSurfSide,iSpec).NE.0) THEN
              ! calculate number of adsorbed particles on reconstructed surface for each species
              numSites = SurfDistInfo(p,q,iSurfSide)%nSites(3) !number of simulated surface atoms
              SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec) = SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec) &
                    + (REAL(SurfModel%SumAdsorbPart(p,q,iSurfSide,iSpec)) * Species(iSpec)%MacroParticleFactor &
                    / maxPart) * REAL(numSites)
            END IF ! SumAdsorbPart!=0
            ! convert to integer adsorbates
            new_adsorbates = INT(SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec))
            IF (new_adsorbates.NE.0) THEN
              ! Adjust tracking of adsorbing background particles
              CALL SMCR_AdjustMapNum(p,q,iSurfSide,new_adsorbates,iSpec,SampleFlag=.TRUE.)
              SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec) = SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec) &
                                                                - new_adsorbates
            END IF
            Adsorption%Coverage(p,q,iSurfSide,iSpec) = coverage_corrected
          END SELECT
!----- 2.2
          IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd))&
              .OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
            SampWall(iSurfSide)%SurfModelState(5+iSpec,p,q) = SampWall(iSurfSide)%SurfModelState(5+iSpec,p,q) &
                                                        + Adsorption%Coverage(p,q,iSurfSide,iSpec)
          END IF
        END DO
      END DO
    END DO
  END DO
!----- 3.
  ! SumEvapPart is nullified in particle emission (surface flux) after inserting particles at corresponding surfaces
  SurfModel%SumEvapPart(:,:,:,:) = SurfModel%SumEvapPart(:,:,:,:) + SurfModel%SumERDesorbed(:,:,1:SurfMesh%nSides,:)
  SurfModel%SumERDesorbed(:,:,:,:) = 0
  SurfModel%SumDesorbPart(:,:,:,:) = 0
  SurfModel%SumAdsorbPart(:,:,:,:) = 0
  SurfModel%SumReactPart(:,:,:,:) = 0
END IF
Adsorption%NumCovSamples = Adsorption%NumCovSamples + 1
!----- 4.
CALL CalcDesorbProb()
CALL CalcAdsorbProb()
#if USE_LOADBALANCE
CALL LBPauseTime(LB_SURF,tLBStart)
#endif /*USE_LOADBALANCE*/

!----- 5.
#ifdef MPI
! communicate coverage and probabilities to halo sides of neighbour procs
CALL ExchangeCoverageInfo()
! communicate distribution to halo-sides of neighbour procs
CALL ExchangeSurfDistInfo()
#endif

END SUBROUTINE UpdateSurfModelVars


SUBROUTINE ReactiveSurfaceTreatment(PartTrajectory,LengthPartTrajectory,alpha,xi,eta,PartID,sideID_IN,flip,IsSpeciesSwap,&
                              adsindex,BCSideID,Opt_Reflected,TriNum)
!===================================================================================================================================
!> Routine for Selection of Surface interaction
!===================================================================================================================================
USE MOD_Globals                ,ONLY: CROSSNORM,UNITVECTOR
USE MOD_Globals_Vars           ,ONLY: Pi
USE MOD_Particle_Tracking_Vars ,ONLY: TriaTracking
USE MOD_Part_Tools             ,ONLY: VELOFROMDISTRIBUTION, CreateParticle
USE MOD_DSMC_Analyze           ,ONLY: CalcWallSample
USE MOD_Particle_Vars          ,ONLY: WriteMacroSurfaceValues
USE MOD_Particle_Vars          ,ONLY: PartState,Species,PartSpecies
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
USE MOD_Particle_Vars          ,ONLY: LastPartPos, PEM
USE MOD_Mesh_Vars              ,ONLY: BC,NGeo
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_Particle_Boundary_Tools,ONLY: PartEnergyToSurface,SurfaceToPartEnergy
USE MOD_Particle_Boundary_Tools,ONLY: TSURUTACONDENSCOEFF
USE MOD_Particle_Boundary_Vars ,ONLY: SurfMesh, dXiEQ_SurfSample, Partbound, SampWall
USE MOD_TimeDisc_Vars          ,ONLY: TEnd, time, dt, RKdtFrac
USE MOD_Particle_Surfaces_vars ,ONLY: SideNormVec,SideType,BezierControlPoints3D
USE MOD_Particle_Surfaces      ,ONLY: CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption, ModelERSpecular, SurfModel
USE MOD_SMCR                   ,ONLY: SMCR_PartAdsorb
USE MOD_SEE                    ,ONLY: SecondaryElectronEmission
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(INOUT)       :: adsindex
REAL,INTENT(INOUT)          :: PartTrajectory(1:3), LengthPartTrajectory, alpha
REAL,INTENT(IN)             :: xi, eta
INTEGER,INTENT(IN)          :: PartID
INTEGER,INTENT(IN)          :: sideID_IN
INTEGER,INTENT(IN)          :: flip
LOGICAL,INTENT(IN)          :: IsSpeciesSwap
INTEGER,INTENT(IN),OPTIONAL :: BCSideID
INTEGER,INTENT(IN),OPTIONAL :: TriNum
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT),OPTIONAL :: Opt_Reflected
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: ProductSpec(2)   ! 1: additional species added or removed from surface
                                                     ! 2: product species of incident particle (also used for simple reflection)
REAL                             :: RanNum
REAL                             :: Xitild,EtaTild
INTEGER                          :: p,q
REAL                             :: n_loc(1:3), tang1(1:3),tang2(1:3)
REAL                             :: Adsorption_prob, Recombination_prob
INTEGER                          :: interactionCase
INTEGER                          :: SurfSideID, SpecID
REAL                             :: Norm_velo!, Norm_Ec
! variables for Energy sampling
REAL                             :: TransArray(1:6),IntArray(1:6), reactionEnthalpie
REAL                             :: oldVelo(1:3)
INTEGER                          :: locBCID
REAL                             :: VeloReal, EtraOld
REAL                             :: EtraWall, EtraNew
REAL                             :: WallVelo(1:3), WallTemp
REAL                             :: TransACC!, VibACC, RotACC
! Polyatomic Molecules
INTEGER                          :: iReact, RecombReactID
REAL                             :: VeloCrad, Fak_D, NewVelo(3)
REAL                             :: Phi, Cmr, VeloCx, VeloCy, VeloCz
REAL                             :: POI_fak, TildTrajectory(3)
CHARACTER(30)                    :: velocityDistribution             ! specifying keyword for velocity distribution
!===================================================================================================================================

! find normal vector two perpendicular tangential vectors (normal_vector points outwards !!!)
IF(PRESENT(BCSideID))THEN
  SELECT CASE(SideType(BCSideID))
  CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
    n_loc=SideNormVec(1:3,BCSideID)
    tang1=UNITVECTOR(BezierControlPoints3D(:,NGeo,0,BCSideID)-BezierControlPoints3D(:,0,0,BCSideID))
    tang2=CROSSNORM(n_loc,tang1)
  CASE(BILINEAR)
    CALL CalcNormAndTangBilinear(n_loc,tang1,tang2,xi,eta,BCSideID)
  CASE(CURVED)
    CALL CalcNormAndTangBezier(n_loc,tang1,tang2,xi,eta,BCSideID)
  END SELECT
ELSE
  IF (TriaTracking) THEN
    CALL CalcNormAndTangTriangle(nVec=n_loc,tang1=tang1,tang2=tang2,TriNum=TriNum,SideID=sideID_IN)
  ELSE
    SELECT CASE(SideType(sideID_IN))
    CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
      n_loc=SideNormVec(1:3,sideID_IN)
        tang1=UNITVECTOR(BezierControlPoints3D(:,NGeo,0,sideID_IN)-BezierControlPoints3D(:,0,0,sideID_IN))
        tang2=CROSSNORM(n_loc,tang1)
    CASE(BILINEAR)
      CALL CalcNormAndTangBilinear(n_loc,tang1,tang2,xi,eta,sideID_IN)
    CASE(CURVED)
      CALL CalcNormAndTangBezier(n_loc,tang1,tang2,xi,eta,sideID_IN)
    END SELECT
    IF(flip.NE.0) n_loc=-n_loc
  END IF
END IF

! check if BC was already crossed
IF(DOT_PRODUCT(n_loc,PartTrajectory).LT.0.)  THEN
  IF(PRESENT(opt_Reflected)) opt_Reflected=.FALSE.
  RETURN
ELSE
  IF(PRESENT(opt_Reflected)) opt_Reflected=.TRUE.
END IF

! additional states
locBCID=PartBound%MapToPartBC(BC(sideID_IN))
! get BC values
WallVelo     = PartBound%WallVelo(1:3,locBCID)
WallTemp     = PartBound%WallTemp(locBCID)

! initialize sampling arrays
TransArray(:) = 0.0
IntArray(:) = 0.0

! compute p and q
! correction of xi and eta, can only be applied if xi & eta are not used later!
IF (TriaTracking) THEN
  p=1 ; q=1
ELSE
  Xitild =MIN(MAX(-1.,xi ),0.99)
  Etatild=MIN(MAX(-1.,eta),0.99)
  p=INT((Xitild +1.0)/dXiEQ_SurfSample)+1
  q=INT((Etatild+1.0)/dXiEQ_SurfSample)+1
END IF

SurfSideID = SurfMesh%SideIDToSurfID(sideID_IN)
SpecID = PartSpecies(PartID)
#if (PP_TimeDiscMethod==42)
! Update wallcollision counter
SurfModel%Info(SpecID)%WallCollCount = SurfModel%Info(SpecID)%WallCollCount + 1
IF (PartBound%SurfaceModel(locBCID).EQ.1) THEN
  SurfModel%Info(SpecID)%Accomodation = SurfModel%Info(SpecID)%Accomodation &
      + (PartBound%TransACC(locBCID) + PartBound%VibACC(locBCID)+ PartBound%RotACC(locBCID))/3.
END IF
#endif

interactionCase = 0
reactionEnthalpie = 0. ! negative at evaporation and positive at condensation
SELECT CASE(PartBound%SurfaceModel(locBCID))
CASE (1)
  Adsorption_prob = Adsorption%ProbAds(p,q,SurfSideID,SpecID)
  CALL RANDOM_NUMBER(RanNum)
  IF ( (Adsorption_prob.GE.RanNum) .AND. &
     (Adsorption%Coverage(p,q,SurfSideID,SpecID).LT.Adsorption%MaxCoverage(SurfSideID,SpecID)) ) THEN
    ProductSpec(2) = SpecID
    interactionCase = 1
  END IF
CASE (2)
  ! Set probabilities
  Adsorption_prob = Adsorption%ProbAds(p,q,SurfSideID,SpecID)
  Recombination_prob = Adsorption%ProbDes(p,q,SurfSideID,SpecID)
  ! check if still enough saved particles on surface
  IF (Adsorption%Coverage(p,q,SurfSideID,SpecID).LE.(-SurfModel%SumAdsorbPart(p,q,SurfSideID,SpecID))) THEN
    Adsorption_prob = Adsorption_prob + Recombination_prob
    Recombination_prob = 0.
  END IF
  ! Decide what happens to colliding particle
  CALL RANDOM_NUMBER(RanNum)
  IF ((Adsorption_prob+Recombination_prob).GE.RanNum) THEN
    CALL RANDOM_NUMBER(RanNum)
    IF ((Adsorption_prob/(Adsorption_prob+Recombination_prob)).GE.RanNum) THEN
      interactionCase = 1
      ProductSpec(1) = 0
      ProductSpec(2) = SpecID
    ELSE
      interactionCase = 3
      DO iReact = Adsorption%DissNum+1,(Adsorption%ReactNum)
        RecombReactID = iReact-Adsorption%DissNum
        IF (Adsorption%RecombReact(2,RecombReactID,SpecID).EQ.Adsorption%ResultSpec(locBCID,SpecID)) THEN
          EXIT
        END IF
      END DO
      ProductSpec(1) = Adsorption%RecombReact(1,RecombReactID,SpecID)
      ProductSpec(2) = Adsorption%RecombReact(2,RecombReactID,SpecID)
      reactionEnthalpie = - Adsorption%EDissBond(iReact,SpecID) * Adsorption%ReactAccomodation(locBCID,SpecID) * BoltzmannConst
    END IF
  END IF
CASE (3)
  Norm_velo = DOT_PRODUCT(PartState(PartID,4:6),n_loc(1:3))
  !Norm_Ec = 0.5 * Species(SpecID)%MassIC * Norm_velo**2 + PartStateIntEn(PartID,1) + PartStateIntEn(PartID,2)
  CALL SMCR_PartAdsorb(p,q,SurfSideID,PartID,Norm_velo,interactionCase,ProductSpec,reactionEnthalpie)
CASE (4)
  ! TODO
CASE (5,6) ! Copied from CASE(1) and adjusted for secondary e- emission (SEE)
           ! 5: SEE by Levko2015
           ! 6: SEE by Pagonakis2016 (originally from Harrower1956)
  ! Get electron emission probability
  CALL SecondaryElectronEmission(PartBound%SurfaceModel(locBCID),PartID,Adsorption_prob,interactionCase,ProductSpec)
  !Adsorption_prob = 1. !Adsorption%ProbAds(p,q,SurfSideID,SpecID)
  ! CALL RANDOM_NUMBER(RanNum)
  ! IF(Adsorption_prob.GE.RanNum)THEN
  !    !  .AND. &
  !    !(Adsorption%Coverage(p,q,SurfSideID,SpecID).LT.Adsorption%MaxCoverage(SurfSideID,SpecID)) ) THEN
  !   ProductSpec(1) = SpecID
  !   ProductSpec(2) = 4!SpecID ! electron
  !   interactionCase = -2 ! perfect elastic scattering + particle creation
  !   WRITE (*,*) "SpecID,ProductSpec(1),ProductSpec(2),interactionCase =", SpecID,ProductSpec(1),ProductSpec(2),interactionCase
  ! END IF
CASE (101) ! constant condensation coefficient
  ProductSpec(2) = SpecID
  Adsorption_prob = Adsorption%ProbAds(p,q,SurfSideID,SpecID)
  CALL RANDOM_NUMBER(RanNum)
  IF ( (Adsorption_prob.GE.RanNum) ) THEN
    interactionCase = 1
  END IF
CASE (102) ! calculate condensation probability by tsuruta2005 and reflection distribution function
  ProductSpec(2) = SpecID
  interactionCase = 4
  velocityDistribution='liquid_refl'
  Norm_velo = DOT_PRODUCT(PartState(PartID,4:6),n_loc(1:3))
  CALL RANDOM_NUMBER(RanNum)
  IF ( (TSURUTACONDENSCOEFF(SpecID,Norm_velo,WallTemp).GE.RanNum) ) THEN
    interactionCase = 1
  END IF
END SELECT

SELECT CASE(interactionCase)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(-4) ! Remove bombarding particle
!-----------------------------------------------------------------------------------------------------------------------------------
  adsindex = 1
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(-3) ! Remove bombarding particle + electron creation
!-----------------------------------------------------------------------------------------------------------------------------------
  SurfModel%SumERDesorbed(p,q,SurfSideID,ProductSpec(2)) = SurfModel%SumERDesorbed(p,q,SurfSideID,ProductSpec(2)) + 1
  adsindex = 1
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(-2) ! Perfect elastic scattering + electron creation
!-----------------------------------------------------------------------------------------------------------------------------------
  SurfModel%SumERDesorbed(p,q,SurfSideID,ProductSpec(2)) = SurfModel%SumERDesorbed(p,q,SurfSideID,ProductSpec(2)) + 1
  adsindex = -1
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(-1) ! Perfect elastic scattering
!-----------------------------------------------------------------------------------------------------------------------------------
  adsindex = -1
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(1) ! Molecular adsorption
!-----------------------------------------------------------------------------------------------------------------------------------
  SurfModel%SumAdsorbPart(p,q,SurfSideID,ProductSpec(2)) = SurfModel%SumAdsorbPart(p,q,SurfSideID,ProductSpec(2)) + 1
  adsindex = 1
#if (PP_TimeDiscMethod==42)
  SurfModel%Info(SpecID)%NumOfAds = SurfModel%Info(SpecID)%NumOfAds + 1
#endif
  CALL PartEnergyToSurface(PartID,SpecID,Transarray,IntArray)
  CALL CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,PartTrajectory,alpha,IsSpeciesSwap,locBCID)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(2) ! dissociative adsorption (particle dissociates on adsorption)
!-----------------------------------------------------------------------------------------------------------------------------------
  SurfModel%SumAdsorbPart(p,q,SurfSideID,ProductSpec(1)) = SurfModel%SumAdsorbPart(p,q,SurfSideID,ProductSpec(1)) + 1
  SurfModel%SumAdsorbPart(p,q,SurfSideID,ProductSpec(2)) = SurfModel%SumAdsorbPart(p,q,SurfSideID,ProductSpec(2)) + 1
  adsindex = 1
#if (PP_TimeDiscMethod==42)
  SurfModel%Info(SpecID)%NumOfAds = SurfModel%Info(SpecID)%NumOfAds + 1
  SurfModel%ProperInfo(SpecID)%HeatFlux(1) = SurfModel%ProperInfo(SpecID)%HeatFlux(1) &
     + reactionEnthalpie/BoltzmannConst
#endif
  IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
    !----  Sampling of reactionEnthalpie
    SampWall(SurfSideID)%SurfModelState(4,p,q) = SampWall(SurfSideID)%SurfModelState(4,p,q) &
                                           + reactionEnthalpie * Species(SpecID)%MacroParticleFactor
    ! Sample reaction counter
    DO iReact = 1,Adsorption%DissNum
      IF (Adsorption%DissocReact(1,iReact,SpecID).EQ.ProductSpec(1) .AND. Adsorption%DissocReact(2,iReact,SpecID).EQ.ProductSpec(2))THEN
        SampWall(SurfSideID)%SurfModelReactCount(iReact,SpecID,p,q)=SampWall(SurfSideID)%SurfModelReactCount(iReact,SpecID,p,q) + 1
      END IF
    END DO
  END IF
  CALL PartEnergyToSurface(PartID,SpecID,Transarray,IntArray)
  CALL CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,PartTrajectory,alpha,IsSpeciesSwap,locBCID)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(3) ! Eley-Rideal reaction (reflecting particle and changes species at contact and reaction partner removed from surface)
!-----------------------------------------------------------------------------------------------------------------------------------
  !SurfModel%SumERDesorbed(p,q,SurfSideID,ProductSpec(2)) = SurfModel%SumERDesorbed(p,q,SurfSideID,ProductSpec(2)) + 1
  SurfModel%SumAdsorbPart(p,q,SurfSideID,ProductSpec(1)) = SurfModel%SumAdsorbPart(p,q,SurfSideID,ProductSpec(1)) - 1
  adsindex = 2
  ! --------
  ! sampling and analyze stuff for heat flux and reaction rates
#if (PP_TimeDiscMethod==42)
  SurfModel%Info(ProductSpec(1))%NumOfDes = SurfModel%Info(ProductSpec(1))%NumOfDes + 1
  SurfModel%Info(ProductSpec(2))%NumOfDes = SurfModel%Info(ProductSpec(2))%NumOfDes + 1
  SurfModel%ProperInfo(SpecID)%HeatFlux(1) = SurfModel%ProperInfo(SpecID)%HeatFlux(1) &
     + reactionEnthalpie/BoltzmannConst
#endif
  IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
    ! Sample recombination reaction counter
    DO iReact = 1,Adsorption%RecombNum
      IF (Adsorption%RecombReact(2,iReact,SpecID).EQ.ProductSpec(2))THEN
        SampWall(SurfSideID)%SurfModelReactCount(Adsorption%DissNum+iReact,SpecID,p,q) = &
            SampWall(SurfSideID)%SurfModelReactCount(Adsorption%DissNum+iReact,SpecID,p,q) + 1
      END IF
    END DO
    !----  Sampling of reactionEnthalpie
    reactionEnthalpie = reactionEnthalpie * Adsorption%ReactAccomodation(locBCID,SpecID)
    SampWall(SurfSideID)%SurfModelState(3,p,q) = SampWall(SurfSideID)%SurfModelState(3,p,q) &
                                           + reactionEnthalpie * Species(SpecID)%MacroParticleFactor
  END IF
  ! --------
  ! reflect particle and change its species
  oldVelo(1:3) = PartState(PartID,4:6)
  CALL PartEnergyToSurface(PartID,SpecID,Transarray,IntArray)
  CALL CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,PartTrajectory,alpha,IsSpeciesSwap,locBCID)

  IF (ModelERSpecular) THEN
    ! perfect velocity reflection
    NewVelo(1:3) = oldVelo(1:3) - 2.*DOT_PRODUCT(oldVelo(1:3),n_loc)*n_loc
    ! mass changes, therefore velocity is scaled because impuls remains the same
    NewVelo(1:3) = NewVelo(1:3) * (Species(ProductSpec(2))%MassIC/Species(SpecID)%MassIC)
  ELSE
    ! diffuse reflection
    TransACC   = PartBound%TransACC(locBCID)
    !VibACC     = PartBound%VibACC(locBCID)
    !RotACC     = PartBound%RotACC(locBCID)
    CALL RANDOM_NUMBER(RanNum)
    VeloCrad    = SQRT(-LOG(RanNum))
    CALL RANDOM_NUMBER(RanNum)
    VeloCz      = SQRT(-LOG(RanNum))
    Fak_D       = VeloCrad**2 + VeloCz**2
    EtraWall    = BoltzmannConst * WallTemp * Fak_D
    VeloReal    = SQRT(DOT_PRODUCT(oldVelo,oldVelo))
    EtraOld     = 0.5 * Species(PartSpecies(PartID))%MassIC * VeloReal**2
    EtraNew     = EtraOld + TransACC * (EtraWall - EtraOld)
    Cmr         = SQRT(2.0 * EtraNew / (Species(ProductSpec(2))%MassIC * Fak_D))
    CALL RANDOM_NUMBER(RanNum)
    Phi     = 2.0 * PI * RanNum
    VeloCx  = Cmr * VeloCrad * COS(Phi) ! tang1
    VeloCy  = Cmr * VeloCrad * SIN(Phi) ! tang2
    VeloCz  = Cmr * VeloCz
    NewVelo = VeloCx*tang1-tang2*VeloCy-VeloCz*n_loc
  END IF
  ! intersection point with surface
  LastPartPos(PartID,1:3) = LastPartPos(PartID,1:3) + PartTrajectory(1:3)*alpha
  ! recompute initial position and ignoring preceding reflections and trajectory between current position and recomputed position
  TildTrajectory=dt*RKdtFrac*oldVelo(1:3)
  POI_fak=1.- (lengthPartTrajectory-alpha)/SQRT(DOT_PRODUCT(TildTrajectory,TildTrajectory))
  ! travel rest of particle vector
  IF (PartBound%Resample(locBCID)) CALL RANDOM_NUMBER(POI_fak) !Resample Equilibirum Distribution
  PartState(PartID,1:3)   = LastPartPos(PartID,1:3) + (1.0 - POI_fak) * dt*RKdtFrac * NewVelo(1:3)
  !----  saving new particle velocity
  PartState(PartID,4:6)   = NewVelo(1:3) + WallVelo(1:3)

  ! recompute trajectory etc
  PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
  lengthPartTrajectory=SQRT(DOT_PRODUCT(PartTrajectory,PartTrajectory))
  PartTrajectory=PartTrajectory/lengthPartTrajectory
  ! set new species
  PartSpecies(PartID) = ProductSpec(2)

  CALL SurfaceToPartEnergy(PartID,ProductSpec(2),WallTemp,Transarray,IntArray)
  CALL CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,PartTrajectory,alpha,IsSpeciesSwap,locBCID,emission_opt=.TRUE.)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(4) ! Distribution function reflection  (reflecting particle according to defined distribution function)
!-----------------------------------------------------------------------------------------------------------------------------------
  adsindex = 2
  ! sampling and analyze stuff for heat flux and reaction rates
#if (PP_TimeDiscMethod==42)
  SurfModel%ProperInfo(SpecID)%HeatFlux(1) = SurfModel%ProperInfo(SpecID)%HeatFlux(1) &
     + reactionEnthalpie/BoltzmannConst
#endif
  IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
    SampWall(SurfSideID)%SurfModelState(5,p,q) = SampWall(SurfSideID)%SurfModelState(5,p,q) &
                                           + ReactionEnthalpie * Species(SpecID)%MacroParticleFactor
  END IF
  oldVelo(1:3) = PartState(PartID,4:6)
  CALL PartEnergyToSurface(PartID,SpecID,Transarray,IntArray)
  CALL CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,PartTrajectory,alpha,IsSpeciesSwap,locBCID)

  ! sample new velocity
  NewVelo(1:3) = VELOFROMDISTRIBUTION(velocityDistribution,SpecID,WallTemp)
  ! important: n_loc points outwards
  PartState(PartID,4:6) = tang1(1:3)*NewVelo(1) + tang2(1:3)*NewVelo(2) - n_Loc(1:3)*NewVelo(3) + WallVelo(1:3)

  ! intersection point with surface
  LastPartPos(PartID,1:3) = LastPartPos(PartID,1:3) + PartTrajectory(1:3)*alpha
  ! recompute initial position and ignoring preceding reflections and trajectory between current position and recomputed position
  TildTrajectory=dt*RKdtFrac*oldVelo(1:3)
  POI_fak=1.- (lengthPartTrajectory-alpha)/SQRT(DOT_PRODUCT(TildTrajectory,TildTrajectory))
  ! travel rest of particle vector
  IF (PartBound%Resample(locBCID)) CALL RANDOM_NUMBER(POI_fak) !Resample Equilibirum Distribution

  ! recompute trajectory etc
  PartState(PartID,1:3)   = LastPartPos(PartID,1:3) + (1.0 - POI_fak) * dt*RKdtFrac * PartState(PartID,4:6)
  PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
  lengthPartTrajectory=SQRT(DOT_PRODUCT(PartTrajectory,PartTrajectory))
  PartTrajectory=PartTrajectory/lengthPartTrajectory

  CALL SurfaceToPartEnergy(PartID,ProductSpec(2),WallTemp,Transarray,IntArray)
  CALL CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,PartTrajectory,alpha,IsSpeciesSwap,locBCID,emission_opt=.TRUE.)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(5) ! reflect incident particle according to distribution function (variable "velocityDistribution" has to be defined)
!-----------------------------------------------------------------------------------------------------------------------------------
  adsindex = 2
  ! sampling and analyze stuff for heat flux and reaction rates
#if (PP_TimeDiscMethod==42)
  SurfModel%Info(ProductSpec(1))%NumOfDes = SurfModel%Info(ProductSpec(1))%NumOfDes + 1
  SurfModel%ProperInfo(SpecID)%HeatFlux(1) = SurfModel%ProperInfo(SpecID)%HeatFlux(1) &
     + reactionEnthalpie/BoltzmannConst
#endif
  IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
    SampWall(SurfSideID)%SurfModelState(5,p,q) = SampWall(SurfSideID)%SurfModelState(5,p,q) &
                                           + ReactionEnthalpie * Species(SpecID)%MacroParticleFactor
  END IF
  oldVelo(1:3) = PartState(PartID,4:6)
  CALL PartEnergyToSurface(PartID,SpecID,Transarray,IntArray)
  CALL CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,PartTrajectory,alpha,IsSpeciesSwap,locBCID)

  ! sample new velocity
  NewVelo(1:3) = VELOFROMDISTRIBUTION(velocityDistribution,SpecID,WallTemp)
  ! important: n_loc points outwards
  PartState(PartID,4:6) = tang1(1:3)*NewVelo(1) + tang2(1:3)*NewVelo(2) - n_Loc(1:3)*NewVelo(3) + WallVelo(1:3)

  ! intersection point with surface
  LastPartPos(PartID,1:3) = LastPartPos(PartID,1:3) + PartTrajectory(1:3)*alpha
  ! recompute initial position and ignoring preceding reflections and trajectory between current position and recomputed position
  TildTrajectory=dt*RKdtFrac*oldVelo(1:3)
  POI_fak=1.- (lengthPartTrajectory-alpha)/SQRT(DOT_PRODUCT(TildTrajectory,TildTrajectory))
  ! travel rest of particle vector
  IF (PartBound%Resample(locBCID)) CALL RANDOM_NUMBER(POI_fak) !Resample Equilibirum Distribution

  ! recompute trajectory etc
  PartState(PartID,1:3)   = LastPartPos(PartID,1:3) + (1.0 - POI_fak) * dt*RKdtFrac * PartState(PartID,4:6)
  PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
  lengthPartTrajectory=SQRT(DOT_PRODUCT(PartTrajectory,PartTrajectory))
  PartTrajectory=PartTrajectory/lengthPartTrajectory

  ! set new species
  PartSpecies(PartID) = ProductSpec(2)

  CALL SurfaceToPartEnergy(PartID,ProductSpec(2),WallTemp,Transarray,IntArray)
  CALL CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,PartTrajectory,alpha,IsSpeciesSwap,locBCID,emission_opt=.TRUE.)

  ! create new particle and assign correct energies
  ! sample new velocity
  NewVelo(1:3) = VELOFROMDISTRIBUTION(velocityDistribution,SpecID,WallTemp)
  ! important: n_loc points outwards
  NewVelo(1:3) = tang1(1:3)*NewVelo(1) + tang2(1:3)*NewVelo(2) - n_Loc(1:3)*NewVelo(3) + WallVelo(1:3)
  CALL CreateParticle(ProductSpec(1),LastPartPos(PartID,1:3),PEM%Element(PartID),NewVelo(1:3),0.,0.,0.)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE DEFAULT ! diffuse reflection
!-----------------------------------------------------------------------------------------------------------------------------------
  adsindex = 0
END SELECT

END SUBROUTINE ReactiveSurfaceTreatment


END MODULE MOD_SurfaceModel
