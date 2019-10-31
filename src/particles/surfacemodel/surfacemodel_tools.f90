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

MODULE MOD_SurfaceModel_Tools
!===================================================================================================================================
! Module for surface model tools
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
INTERFACE CalcEvappartNum
  MODULE PROCEDURE CalcEvappartNum
END INTERFACE

INTERFACE CalcAdsorbProb
  MODULE PROCEDURE CalcAdsorbProb
END INTERFACE

INTERFACE CalcDesorbProb
  MODULE PROCEDURE CalcDesorbProb
END INTERFACE

INTERFACE Calc_Adsorb_Heat
  MODULE PROCEDURE Calc_Adsorb_Heat
END INTERFACE

INTERFACE CalcDissRecombActEnergy
  MODULE PROCEDURE CalcDissRecombActEnergy
END INTERFACE

INTERFACE Calc_E_Act
  MODULE PROCEDURE Calc_E_Act
END INTERFACE

INTERFACE CalcAdsorbReactProb
  MODULE PROCEDURE CalcAdsorbReactProb
END INTERFACE

INTERFACE SpaceOccupied
  MODULE PROCEDURE SpaceOccupied
END INTERFACE

INTERFACE UpdateSurfPos
  MODULE PROCEDURE UpdateSurfPos
END INTERFACE

INTERFACE SampleAdsorptionHeat
  MODULE PROCEDURE SampleAdsorptionHeat
END INTERFACE

INTERFACE IsReactiveSurface
  MODULE PROCEDURE IsReactiveSurface
END INTERFACE

INTERFACE SurfaceHasModelNum
  MODULE PROCEDURE SurfaceHasModelNum
END INTERFACE

INTERFACE SMCR_AdjustMapNum
  MODULE PROCEDURE SMCR_AdjustMapNum
END INTERFACE

PUBLIC :: CalcEvapPartNum
PUBLIC :: CalcAdsorbProb
PUBLIC :: CalcDesorbProb
PUBLIC :: Calc_Adsorb_Heat
PUBLIC :: Calc_E_Act
PUBLIC :: CalcDissRecombActEnergy
PUBLIC :: CalcAdsorbReactProb
PUBLIC :: SpaceOccupied
PUBLIC :: UpdateSurfPos
PUBLIC :: SampleAdsorptionHeat
PUBLIC :: IsReactiveSurface
PUBLIC :: SurfaceHasModelNum
PUBLIC :: SMCR_AdjustMapNum
!===================================================================================================================================

CONTAINS

SUBROUTINE CalcEvapPartNum()
!===================================================================================================================================
!> calculation of number of evaporating/desorbing particles when mean surface densities are used (mean probabilities)
!===================================================================================================================================
USE MOD_Globals_Vars           ,ONLY: PI, BoltzmannConst
USE MOD_Particle_Vars          ,ONLY: nSpecies, Species
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption, surfmodel, SpecSurf
USE MOD_Part_Tools             ,ONLY: TSURUTACONDENSCOEFF
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, PartBound
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_TimeDisc_Vars          ,ONLY: dt
#if USE_LOADBALANCE
USE MOD_Particle_Mesh_Vars     ,ONLY: PartSideToElem
USE MOD_LoadBalance_Vars       ,ONLY: nSurfacePartsPerElem, PerformLBSample
#endif /*USE_LOADBALANCE*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                          :: iSurfSide, iSpec, p, q, NPois, PartEvapInfo
REAL                             :: WallPartNum, PartEvap, RanNum, Tpois, LiquidSurfTemp
REAL                             :: A, B, C, pressureVapor
INTEGER                          :: iPart, PartEvap_temp
REAl                             :: veloPart, sigma
#if USE_LOADBALANCE
INTEGER                          :: ElemID
#endif /*USE_LOADBALANCE*/
!----------------------------------------------------------------------------------------------------------------------------------!
IF (.NOT.SurfMesh%SurfOnProc) RETURN

DO iSpec = 1,nSpecies
  DO iSurfSide = 1,SurfMesh%nMasterSides

#if USE_LOADBALANCE
    IF(PerformLBSample) THEN
      ElemID = PartSideToElem(S2E_ELEM_ID,SurfMesh%SurfIDToSideID(iSurfSide))
      nSurfacePartsPerElem(ElemID) = nSurfacePartsPerElem(ElemID) + 1
    END IF
#endif /*USE_LOADBALANCE*/
    DO q = 1,nSurfSample
      DO p = 1,nSurfSample
        SELECT CASE (SurfaceHasModelNum(iSurfSide))
        CASE(1)
          IF (Adsorption%Coverage(p,q,iSurfSide,iSpec).GT.0) THEN
            WallPartNum = Adsorption%Coverage(p,q,iSurfSide,iSpec) * Adsorption%DensSurfAtoms(iSurfSide) &
                      * SurfMesh%SurfaceArea(p,q,iSurfSide) &
                      / Species(iSpec)%MacroParticleFactor
            PartEvap = WallPartNum * Adsorption%ProbDes(p,q,iSurfSide,iSpec)
          ELSE
            PartEvap = 0.
          END IF
        CASE(101,102)
          IF (Adsorption%SurfaceSpec(PartBound%MapToPartBC(BC( SurfMesh%SurfIDToSideID(iSurfSide) )),iSpec)) THEN
            LiquidSurfTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC( SurfMesh%SurfIDToSideID(iSurfSide) )))
            ! Antoine parameters defined in ini file are chosen so pressure given is in bar
            A = SpecSurf(iSpec)%ParamAntoine(1)
            B = SpecSurf(iSpec)%ParamAntoine(2)
            C = SpecSurf(iSpec)%ParamAntoine(3)
            ! Use Antoine Eq. to calculate pressure vapor
            pressureVapor = 10 ** (A- B/(C+LiquidSurfTemp)) * 1e5 !transformation bar -> Pa
            ! Use Hertz-Knudsen equation to calculate number of evaporating liquid particles from surface
            PartEvap = pressureVapor / ( 2*PI*Species(iSpec)%MassIC*BoltzmannConst*LiquidSurfTemp)**0.5 &
                     * SurfMesh%SurfaceArea(p,q,iSurfSide) / Species(iSpec)%MacroParticleFactor * dt &
                     * Adsorption%ProbDes(p,q,iSurfSide,iSpec)
            IF (SurfaceHasModelNum(iSurfSide).EQ.102) THEN
              sigma=SQRT(BoltzmannConst*LiquidSurfTemp/Species(iSpec)%MassIC)
              PartEvap_temp = INT(PartEvap,4)
              PartEvap = 0
              DO iPart = 1,PartEvap_temp
                CALL RANDOM_NUMBER(RanNum)
                veloPart = SQRT(-2*LOG(RanNum))*sigma
                CALL RANDOM_NUMBER(RanNum)
                IF ( (TSURUTACONDENSCOEFF(iSpec,veloPart,LiquidSurfTemp).GE.RanNum) ) PartEvap=PartEvap+1
              END DO
            END IF
            WallPartNum = PartEvap
          ELSE
            WallPartNum = 0
          END IF
        CASE DEFAULT
          WallPartNum = 0
        END SELECT

        IF (WallPartNum.GT.0) THEN
          IF (PartEvap.GT.WallPartNum .OR. PartEvap.LT.0.) THEN
            PartEvapInfo = INT(WallPartNum)
          ELSE
            CALL RANDOM_NUMBER(RanNum)
            IF (EXP(-PartEvap).LE.TINY(PartEvap) &
#if (PP_TimeDiscMethod==42)
                .OR. Adsorption%TPD &
#endif
                ) THEN
              PartEvapInfo = INT(PartEvap + RanNum)
            ELSE !poisson-sampling instead of random rounding (reduces numeric non-equlibrium effects [Tysanner and Garcia 2004]
              Npois=0
              Tpois=1.0
              DO
                Tpois=RanNum*Tpois
                IF (Tpois.LT.TINY(Tpois)) THEN
                  PartEvapInfo = INT(PartEvap + RanNum)
                  EXIT
                END IF
                IF (Tpois.GT.EXP(-PartEvap)) THEN
                  Npois=Npois+1
                  CALL RANDOM_NUMBER(RanNum)
                ELSE
                  PartEvapInfo = Npois
                  EXIT
                END IF
              END DO
            END IF
          END IF !PartEvap.GT.WallPartNum
          SurfModel%SumEvapPart(p,q,iSurfSide,iSpec) = SurfModel%SumEvapPart(p,q,iSurfSide,iSpec) + INT(PartEvapInfo)
        END IF !WallPartNum.GT.0
      END DO
    END DO
  END DO
END DO

END SUBROUTINE CalcEvapPartNum

SUBROUTINE CalcAdsorbProb()
!===================================================================================================================================
!> Calculcation of adsorption probability for different model (wallmodel 1 and 2)
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption, SurfModel
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, PartBound
#if (PP_TimeDiscMethod==42)
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! Local variable declaration
INTEGER                          :: SurfSide, iSpec, p, q
REAL                             :: Theta_req, Kfactor, S_0
INTEGER                          :: PartBoundID
!===================================================================================================================================
DO iSpec=1,nSpecies
  DO SurfSide=1,SurfMesh%nMasterSides
    PartBoundID = PartBound%MapToPartBC(BC(SurfMesh%SurfIDToSideID(SurfSide)))
    IF (.NOT.PartBound%Reactive(PartboundID)) CYCLE
    DO q = 1,nSurfSample
      DO p = 1,nSurfSample
        SELECT CASE (PartBound%SurfaceModel(PartboundID))
!----------------------------------------------------------------------------------------------------------------------------------!
        CASE (1) ! Kisluik Sticking Model from Kolasinski's Surface Science (book)
!----------------------------------------------------------------------------------------------------------------------------------!
          ! enhance later to co-adsorption
          Theta_req = (1.0 - Adsorption%Coverage(p,q,SurfSide,iSpec)/Adsorption%MaxCoverage(SurfSide,iSpec)) &
                    **Adsorption%Adsorbexp(SurfSide,iSpec)
          !----- kann später auf von Wandtemperatur abhängige Werte erweitert werden
          Kfactor = Adsorption%PrefactorStick(SurfSide,iSpec)
          S_0 = Adsorption%InitStick(SurfSide,iSpec)
          !-----
          IF (Theta_req.EQ.0) THEN
            Adsorption%ProbAds(p,q,SurfSide,iSpec) = 0.
          ELSE
            Adsorption%ProbAds(p,q,SurfSide,iSpec) = S_0 / (1.0 + Kfactor * ( 1.0/Theta_req - 1.0))
          END IF
!----------------------------------------------------------------------------------------------------------------------------------!
        CASE (2) ! Recombination Model described by Laux
!----------------------------------------------------------------------------------------------------------------------------------!
          Adsorption%ProbAds(p,q,SurfSide,iSpec) = Adsorption%ReactCoeff(PartBoundID,iSpec)-Adsorption%ProbDes(p,q,SurfSide,iSpec)
!----------------------------------------------------------------------------------------------------------------------------------!
        CASE (101) ! simple condensation coefficient
!----------------------------------------------------------------------------------------------------------------------------------!
          Adsorption%ProbAds(p,q,SurfSide,iSpec) = Adsorption%ReactCoeff(PartBoundID,iSpec)
        END SELECT
!----------------------------------------------------------------------------------------------------------------------------------!
        SurfModel%Info(iSpec)%MeanProbAds = SurfModel%Info(iSpec)%MeanProbAds+Adsorption%ProbAds(p,q,SurfSide,iSpec)
        SurfModel%Info(iSpec)%MeanProbAdsCount = SurfModel%Info(iSpec)%MeanProbAdsCount + 1
      END DO
    END DO
  END DO
END DO
END SUBROUTINE CalcAdsorbProb


SUBROUTINE CalcDesorbProb()
!===================================================================================================================================
!> Calculcation of desorption probability for different model (wallmodel 1 and 2)
!===================================================================================================================================
USE MOD_Globals_Vars           ,ONLY: PlanckConst, BoltzmannConst
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption, SurfModel
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, PartBound
USE MOD_TimeDisc_Vars          ,ONLY: dt
#if (PP_TimeDiscMethod==42)
USE MOD_TimeDisc_Vars          ,ONLY: iter
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! Local variable declaration
INTEGER                          :: SurfSide, iSpec, p, q
REAL                             :: Theta, nu_des, rate, WallTemp
REAL                             :: E_des
INTEGER                          :: PartBoundID, iReactNum, RecombReactID, jSpec, kSpec
!===================================================================================================================================
DO SurfSide=1,SurfMesh%nMasterSides
  PartBoundID = PartBound%MapToPartBC(BC(SurfMesh%SurfIDToSideID(SurfSide)))
  IF (.NOT.PartBound%Reactive(PartboundID)) CYCLE
! special TPD (temperature programmed desorption) temperature adjustment routine
#if (PP_TimeDiscMethod==42)
  IF (Adsorption%TPD) THEN
    WallTemp = PartBound%WallTemp(PartBoundID) + (Adsorption%TPD_beta * iter * dt)
    Adsorption%TPD_Temp = Walltemp
  ELSE
    WallTemp = PartBound%WallTemp(PartBoundID)
  END IF
#else
  WallTemp = PartBound%WallTemp(PartBoundID)
#endif

  DO iSpec = 1,nSpecies
    DO q = 1,nSurfSample
      DO p = 1,nSurfSample
        SELECT CASE (PartBound%SurfaceModel(PartboundID))
!----------------------------------------------------------------------------------------------------------------------------------!
        CASE (1) ! Polanyi-Wigner-eq. from Kolasinski's Surface Science (book)
!----------------------------------------------------------------------------------------------------------------------------------!
          Theta = Adsorption%Coverage(p,q,SurfSide,iSpec)! / Adsorption%MaxCoverage(SurfSide,iSpec)
          !----- kann später auf von Wandtemperatur/Translationsenergie abhängige Werte erweitert werden
          E_des = Adsorption%DesorbEnergy(SurfSide,iSpec) + Adsorption%Intensification(SurfSide,iSpec) * Theta
          nu_des = 10**(Adsorption%Nu_a(SurfSide,iSpec) + Adsorption%Nu_b(SurfSide,iSpec) * Theta)!/10000
          !-----
          rate = nu_des &!*(Adsorption%DensSurfAtoms(SurfSide)**(Adsorption%Adsorbexp(SurfSide,iSpec)-1)) &
                        * (Theta**Adsorption%Adsorbexp(SurfSide,iSpec)) * exp(-E_des/WallTemp)
          IF (Theta.GT.0) THEN
            Adsorption%ProbDes(p,q,SurfSide,iSpec) = rate * dt /Theta
          ELSE
            Adsorption%ProbDes(p,q,SurfSide,iSpec) = 0.0
          END IF
!----------------------------------------------------------------------------------------------------------------------------------!
        CASE (2) ! Recombination Model described by Fasoulas/Laux
!----------------------------------------------------------------------------------------------------------------------------------!
          jSpec = 0 ! initialize reaction partner with zero
          DO iReactNum = Adsorption%DissNum+1,(Adsorption%ReactNum)
            RecombReactID = iReactNum-Adsorption%DissNum
            ! resulting species
            kSpec = Adsorption%RecombReact(2,RecombReactID,iSpec)
            IF (kSpec.EQ.Adsorption%ResultSpec(PartBoundID,iSpec)) THEN
              ! reaction partner
              jSpec = Adsorption%RecombReact(1,RecombReactID,iSpec)
              EXIT
            END IF
          END DO
          IF (jSpec.LE.0) THEN
            Adsorption%ProbDes(p,q,SurfSide,iSpec) = 0.
          ELSE
            IF (Adsorption%Coverage(p,q,SurfSide,jSpec).LE.0) THEN
              Adsorption%ProbDes(p,q,SurfSide,iSpec) = 0.
            ELSE
              Adsorption%ProbDes(p,q,SurfSide,iSpec) = Adsorption%ReactCoeff(PartBoundID,iSpec) &
                  * ( 1 - exp( - Adsorption%Coverage(p,q,SurfSide,jSpec) ) )
            END IF
          END IF
!----------------------------------------------------------------------------------------------------------------------------------!
        CASE (101,102) ! simple condensation coefficient
!----------------------------------------------------------------------------------------------------------------------------------!
          Adsorption%ProbDes(p,q,SurfSide,iSpec) = Adsorption%ReactCoeff(PartBoundID,iSpec)
        END SELECT
        SurfModel%Info(iSpec)%MeanProbDes = SurfModel%Info(iSpec)%MeanProbDes + Adsorption%ProbDes(p,q,SurfSide,iSpec)
        SurfModel%Info(iSpec)%MeanProbDesCount = SurfModel%Info(iSpec)%MeanProbDesCount + 1
      END DO
    END DO
  END DO
END DO
END SUBROUTINE CalcDesorbProb


REAL FUNCTION Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Species,Surfpos,IsAdsorption)
!===================================================================================================================================
!> Calculates the Heat of adsorption for given species and given surface position
!> Uses UBI-QEP model approach with Surface Monte Carlo Reconstruction
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_Particle_Boundary_vars ,ONLY: PartBound, SurfMesh
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_DSMC_Vars              ,ONLY: SpecDSMC
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption, SurfDistInfo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)            :: subsurfxi, subsurfeta, SurfSideID
INTEGER, INTENT(IN)            :: Species, Surfpos
LOGICAL, INTENT(IN)            :: IsAdsorption
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: Coordination, i, j, Indx, Indy, PartBoundID
REAL , ALLOCATABLE             :: x(:)
INTEGER , ALLOCATABLE          :: m(:)
INTEGER                        :: bondorder
REAL                           :: D_AB, D_AX, D_BX
REAL                           :: Heat_A, Heat_B
REAL                           :: A, B, sigma, sigma_m
! additional parameters for attraction between associative neighbours
REAL , ALLOCATABLE             :: D_AL(:)
REAL , ALLOCATABLE             :: attractBondOrder(:,:)
INTEGER , ALLOCATABLE          :: Neigh_bondorder(:)
REAL                           :: HeatAttraction
INTEGER                        :: neighSpec, neighSpec2, Coord2, Coord3, iRecombReact, ReactNum, nNeigh_interactions
INTEGER                        :: l, k, NeighPos
!===================================================================================================================================
PartBoundID = PartBound%MapToPartBC(BC(SurfMesh%SurfIDToSideID(SurfSideID)))
IF (.NOT.PartBound%Reactive(PartboundID)) CALL Abort(&
__STAMP__,&
'Calc_Adsorb_Heat_ERROR: Given SurfSideID is not reactive',SurfSideID)
Coordination = Adsorption%Coordination(PartBoundID,Species)
ALLOCATE( x(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom) )
ALLOCATE( m(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom) )
x(:) = 1. ! averaged bond-index of adsorbate with respective surface atom
m(:) = 1  ! number of adsorbates belonging to the respective surface atom
Calc_Adsorb_Heat = 0.
sigma = 0.
IF (Surfpos.GT.0) THEN
  DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom
    Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%BondAtomIndx(Surfpos,j)
    Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%BondAtomIndy(Surfpos,j)
    bondorder = 0
    DO i = 1,nSpecies
      bondorder = bondorder + SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(i,Indx,Indy)
    END DO
    IF (IsAdsorption) THEN
      ! calculate bond order for heat of adsorption (in case of adsorption treatment)
      m(j) = (bondorder + 1) !adsorbing particle itself has to be added
    ELSE
      ! calculate bond order for heat of adsorption (in case of desorption treatment)
      m(j) = bondorder
    END IF
    IF (m(j).LT.1) THEN !should never occur except calculating desorb heat for empty site (IsAdsorption=FALSE)
      CALL Abort(&
__STAMP__,&
'Calc_Adsorb_Heat_ERROR: Calculating Heat of adsorbtion not possible for surface position',Surfpos)
    END IF
  END DO
END IF

! calculate scaling factor for M-A bond weakening due to lateral interactions
#if (PP_TimeDiscMethod==42)
IF (Adsorption%LateralInactive) THEN
  sigma_m = 1.
ELSE
#endif
  sigma_m = 0.
  ! calculate local scaling factor for chosen surface site
  DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom
  !     x(j) = x(j) / REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom)
  !     sigma = sigma + (2.*x(j) - x(j)**2.) * (2.*(1./REAL(m(j))) - (1./REAL(m(j)))**2.)
    sigma_m = sigma_m + (2.*(1./REAL(m(j))) - (1./REAL(m(j)))**2) &
                    / REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom)
  END DO
#if (PP_TimeDiscMethod==42)
END IF
#endif

! calculate additional heat of adsorption for direct A-A interaction (attraction of associating adsorbates)
nNeigh_interactions = 0
HeatAttraction = 0.
IF (Adsorption%EnableAdsAttraction) THEN
  IF ((Adsorption%RecombNum.GT.0) .AND. (Surfpos.GT.0) ) THEN
    ALLOCATE(attractBondOrder(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom, &
             1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nNeighbours))
    attractBondOrder(:,:) = 0.
    ALLOCATE(D_AL(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nNeighbours),&
             Neigh_bondorder(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nNeighbours))
    D_AL(:) = 0.
    Neigh_bondorder(:) = 0
    nNeigh_interactions = 0
    ! define dissociation bond energies of neighbours and count interacting neighbours
    DO l = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nNeighbours
      IF (.NOT.SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%IsNearestNeigh(SurfPos,l)) CYCLE
      Coord2 = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%NeighSite(Surfpos,l)
      IF (Coord2.GT.0) THEN
        NeighPos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%NeighPos(Surfpos,l)
        neighSpec = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%Species(NeighPos)
        IF (neighSpec.NE.0) THEN
          DO iRecombReact = 1,Adsorption%RecombNum
            ReactNum = iRecombReact + Adsorption%DissNum
            IF ( neighSpec.EQ.Adsorption%RecombReact(1,iRecombReact,Species) .AND. &
                 (Adsorption%RecombReact(2,iRecombReact,Species).NE.0)) THEN
              DO i = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom
                Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%BondAtomIndx(Surfpos,i)
                Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%BondAtomIndy(Surfpos,i)
                DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%nInterAtom
                  IF (Indx.EQ.SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%BondAtomIndx(NeighPos,j) &
                      .AND. Indy.EQ.SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%BondAtomIndy(NeighPos,j) ) THEN
                    attractBondOrder(i,l) = attractBondOrder(i,l) + 1
                  END IF
                END DO
              END DO
              D_AL(l) = Adsorption%EDissBond(ReactNum,Species)
              nNeigh_interactions = nNeigh_interactions + 1
              CYCLE
            END IF
          END DO
          ! estimtate bondorder of neighbours
          DO k = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%nNeighbours
            IF (.NOT.SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%IsNearestNeigh(NeighPos,k)) CYCLE
            Coord3 = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%NeighSite(NeighPos,k)
            IF (Coord3.GT.0) THEN
              neighSpec2 = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord3)%Species( &
                           SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%NeighPos(NeighPos,k))
              IF ( (neighSpec2.NE.0) ) THEN
                DO iRecombReact = 1,Adsorption%RecombNum
                  IF ( (neighSpec2.EQ.Adsorption%RecombReact(1,iRecombReact,neighSpec)) .AND. &
                       (Adsorption%RecombReact(2,iRecombReact,neighSpec).NE.0)) THEN
                    Neigh_bondorder(l) = Neigh_bondorder(l) + 1
                    CYCLE
                  END IF
                END DO
              END IF
            END IF
          END DO
          IF (IsAdsorption) THEN
            Neigh_bondorder(:) = Neigh_bondorder(:) + 1
          END IF
        END IF
      END IF
    END DO
    ! normalize sum of attractbondorder for each surface atomto one
    DO i = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom
      IF (SUM(attractBondOrder(i,:)).GT.1) THEN
        attractBondOrder(i,:) = attractBondOrder(i,:) / SUM(attractBondOrder(i,:))
      END IF
    END DO
    ! calculate interaction energy between adsorbate and neighbours
    IF (nNeigh_interactions.GT.0) THEN
      DO l = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nNeighbours
        IF (Neigh_bondorder(l).GT.0) THEN
          D_AL(l) = D_AL(l) * ( 2. - 1./REAL(Neigh_bondorder(l)) ) / REAL(Neigh_bondorder(l))
        END IF
        DO i = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom
          HeatAttraction = HeatAttraction + 0.5*D_AL(l) * (2*attractbondOrder(i,l) &
            / REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom) - (attractBondOrder(i,l) &
            / REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom))**2)
        END DO
      END DO
    END IF

    DEALLOCATE(D_AL)
    DEALLOCATE(Neigh_bondorder)
    !DEALLOCATE(delta)
    IF (nNeigh_interactions.EQ.0) DEALLOCATE(attractBondOrder)
  END IF
END IF

! caluclate bond index for M-A interaction
DO i = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom
  IF (nNeigh_interactions.GT.0) THEN
    x(i) = (1.-SUM(attractBondOrder(i,:))) / REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom)
  ELSE
    x(i) = 1. / REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom)
  END IF
  sigma = sigma + 2.*x(i) - x(i)**2
END DO

! Testing if the adsorption particle is an atom or molecule, if molecule: is it polyatomic?
! and calculate right heat of adsorption to surface atoms
Heat_A = Adsorption%HeatOfAdsZero(PartBoundID,Species)
D_AB = Adsorption%EDissBond(0,Species)
IF(SpecDSMC(Species)%InterID.EQ.2) THEN
  ! Cases for binding type
  SELECT CASE(Adsorption%DiCoord(PartBoundID,Species))
  CASE(1) ! strong bonding
    Calc_Adsorb_Heat = (Heat_A*sigma)**2/(D_AB+Heat_A*sigma) * sigma_m
  CASE(2) ! weak bonding
    Calc_Adsorb_Heat = Heat_A**2/(D_AB+Heat_A/REAL(1./(2-sigma))) * sigma_m
  CASE(3) ! intermediate binding (something between strong and weak)
    Calc_Adsorb_Heat = ( (Heat_A*sigma)**2/(D_AB+Heat_A*sigma) + Heat_A**2/(D_AB+Heat_A/REAL(1./(2-sigma))) )/2. * sigma_m
  CASE(4) ! parallel to surface, each molecule atom is bound to one surface atom (bridge site, acceptor adsorbate)
    IF(SpecDSMC(Species)%PolyatomicMol) THEN
      ! dicoordination e.g. (HCOOH --> M--(HC)O-O(H)--M) (M--O bond)
      !D_AB = Adsorption%EDissBond(0,Species) ! Bond O-O
      D_AX = Adsorption%EDissBondAdsorbPoly(0,Species) ! Bond HC--O
      D_BX = Adsorption%EDissBondAdsorbPoly(1,Species) ! Bond O--H
      Heat_A = Adsorption%HeatOfAdsZero(PartBoundID,Species)
      A = Heat_A**2./(D_AX+D_AB+Heat_A)
      Heat_B = Adsorption%HeatOfAdsZero(PartBoundID,Species)
      B = Heat_B**2./(D_BX+D_AB+Heat_B)
      Calc_Adsorb_Heat = ( A*B*( A + B ) + D_AB*( A - B )**2. ) / ( A*B + D_AB*( A + B ) ) * sigma_m
    ELSE
      Heat_B = Adsorption%HeatOfAdsZero(PartBoundID,Species)
      A = Heat_A**2 * ( Heat_A + 2.*Heat_B ) / ( Heat_A + Heat_B )**2
      B = Heat_B**2 * ( Heat_B + 2.*Heat_A ) / ( Heat_A + Heat_B )**2
      Calc_Adsorb_Heat = ( A*B*( A + B ) + D_AB*( A - B )**2 ) / ( A*B + D_AB*( A + B ) ) * sigma_m
    END IF
  CASE(5) ! parallel to surface, each molecule atom is bound to one surface atom (on top site, donor adsorbate)
    IF(SpecDSMC(Species)%PolyatomicMol) THEN
      Heat_B = Adsorption%HeatOfAdsZero(PartBoundID,Species)
      D_AX = Adsorption%EDissBondAdsorbPoly(0,Species) ! Bond HC--O
      D_BX = Adsorption%EDissBondAdsorbPoly(1,Species) ! Bond O--H
      Heat_A = Heat_A * 3./4.
      Heat_B = Heat_B * 3./4.
      A = Heat_A**2./(D_AX+Heat_A)
      B = Heat_B**2./(D_BX+Heat_B)
      Calc_Adsorb_Heat = ( A*B*( A + B ) + D_AB*( A - B )**2. ) / ( A*B + D_AB*( A + B ) ) * sigma_m
    ELSE
      Heat_B = Adsorption%HeatOfAdsZero(PartBoundID,Species)
      A = Heat_A**2 * ( Heat_A + 2.*Heat_B ) / ( Heat_A + Heat_B )**2
      B = Heat_B**2 * ( Heat_B + 2.*Heat_A ) / ( Heat_A + Heat_B )**2
      Calc_Adsorb_Heat = ( A*B*( A + B ) + D_AB*( A - B )**2 ) / ( A*B + D_AB*( A + B ) ) * sigma_m
    END IF
  CASE(6) ! parallel to surface, each molecule atom is bound to both surface atoms (bridge site, donor adsorbate)
    Heat_B = Adsorption%HeatOfAdsZero(PartBoundID,Species)
    A = Heat_A *3./4.
    B = Heat_B *3./4.
    Calc_Adsorb_Heat = 2*( A*B*( A + B ) + 2*D_AB*( A - B )**2 ) / ( A*B + 2*D_AB*( A + B ) ) * sigma_m
  CASE(7) ! chelating bridge, e.g. (NO2 --> M--O-N-O--M) no direct bonding between adsorbate ends
    IF(SpecDSMC(Species)%PolyatomicMol) THEN
      D_AX = Adsorption%EDissBondAdsorbPoly(0,Species) ! Bond O--N
      D_BX = Adsorption%EDissBondAdsorbPoly(1,Species) ! Bond N--O
      Heat_A = Adsorption%HeatOfAdsZero(PartBoundID,Species)
      Heat_A = Heat_A**2/(D_AX+Heat_A)
      Heat_B = Adsorption%HeatOfAdsZero(PartBoundID,Species)
      Heat_B = Heat_B**2/(D_BX+Heat_B)
      A = Heat_A**2. * ( Heat_A + 2.*Heat_B ) / ( Heat_A + Heat_B )**2.
      B = Heat_B**2. * ( Heat_B + 2.*Heat_A ) / ( Heat_A + Heat_B )**2.
      Calc_Adsorb_Heat = (A + B) * sigma_m
    END IF
  CASE DEFAULT
    CALL abort(&
__STAMP__&
,"ERROR in Calc_Adsorb_Heat: wrong dicoord for species:",Species)
  END SELECT
ELSE
  Calc_Adsorb_Heat = (Heat_A*sigma) * sigma_m
END IF

! calculate total adsorption heat
IF(nNeigh_interactions.GT.0) THEN
  Calc_Adsorb_Heat = Calc_Adsorb_Heat + HeatAttraction
  DEALLOCATE(attractBondOrder)
END IF

DEALLOCATE(x,m)

END FUNCTION Calc_Adsorb_Heat


REAL FUNCTION Calc_E_Act(Heat_Product_A,Heat_Product_B,Heat_Reactant_A,Heat_Reactant_B,&
                         D_Product_A,D_Product_B,D_Reactant_A,D_Reactant_B)
!===================================================================================================================================
!> Calculates the Activation energy for a given reaction
!> A_Reactant_ads + B_Reactant_ads --> A_Product_ads + B_Product_ads
!> Adsorption --> forward reaction
!> Forward reaction is defined by D_Educt > D_Products
!> Examples:
!> (1)
!> O2 desorbed directly to gasphase from reaction of two O (O_ads + O_ads -> O2_g):
!> ==> forward reaction: O2_g + (-)_ads -> O_ads + O_ads
!> ==> IsAdsorption = .FALSE.
!> ==> Heat_Reactant_A = Heat_O2_g = 0. | Heat_Product_A_ads = Heat_Product_B_ads = Heat_O_ads
!> (2)
!> adsorbed CH radical reacts with adsorbed O-atom to adsorbed C-atom and OH-radical (CH_ads + O_ads -> C_ads + OH_ads):
!> ==> forward reaction: CH_ads + O_ads -> C_ads + OH_ads
!> ==> IsAdsorption = .TRUE.
!> ==> Heat_Reactant_A = Heat_CH_ads | Heat_Reactant_B = Heat_O_ads | Heat_Product_A = Heat_C_ads | Heat_Product_B = Heat_OH_ads
!> (3)
!> adsorbed OH radical reacts with adsorbed C-atom to adsorbed O-atom and gasphase CH-radical (OH_ads + C_ads -> O_ads + CH_g):
!> ==> forward reaction: CH_g + O_ads -> C_ads + OH_ads
!> ==> IsAdsorption = .FALSE.
!> ==> Heat_Reactant_A = Heat_CH_g = 0. | Heat_Reactant_B = Heat_O_ads | Heat_Product_A = Heat_C_ads | Heat_Product_B = Heat_OH_ads
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)               :: Heat_Product_A, Heat_Product_B, Heat_Reactant_A, Heat_Reactant_B
REAL, INTENT(IN)               :: D_Product_A, D_Product_B, D_Reactant_A, D_Reactant_B
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Delta_H
LOGICAL                        :: Forward
!===================================================================================================================================
! decide if forward or reverse reaction
Forward = .FALSE.
IF ( (D_Reactant_A +D_Reactant_B -D_Product_A -D_Product_B).GT.0. ) Forward = .TRUE.

IF (Forward) THEN
  Delta_H = ( Heat_Reactant_A +Heat_Reactant_B -Heat_Product_A -Heat_Product_B ) &
          + ( D_Reactant_A +D_Reactant_B -D_Product_A -D_Product_B )
  Calc_E_Act = 0.5 * ( Delta_H + (Heat_Product_A*Heat_Product_B / (Heat_Product_A+Heat_Product_B)) )
ELSE
  Delta_H = ( Heat_Product_A +Heat_Product_B -Heat_Reactant_A -Heat_Reactant_B ) &
          + ( +D_Product_A +D_Product_B -D_Reactant_A -D_Reactant_B )
  Calc_E_Act = 0.5 * ( Delta_H + (Heat_Reactant_A*Heat_Reactant_B / (Heat_Reactant_A+Heat_Reactant_B)) )
  IF (Calc_E_Act.LT.0.) Calc_E_Act = 0.
  Calc_E_Act = Calc_E_Act - Delta_H
END IF
IF (Calc_E_Act.LT.0.) Calc_E_Act = 0.

END FUNCTION Calc_E_Act


REAL FUNCTION CalcDissRecombActEnergy(HeatProduct,HeatReactantA,HeatReactantB,DProduct,ProdSpec)
!===================================================================================================================================
!> Calculates the Activation energy for a dissociation or recombination reaction
!> specify heat and dissociation bond energy of product (recombined) species and heat of adosrptions of educt (dissociated) species
!>   DISSOCIATION: O2_ads + (-)_ads -> O_ads + O_ads
!>   HeatProduct = Heat_O2 | HeatReactantA = HeatReactantB = Heat_O
!>   RECOMBINATION: O_ads + O_ads -> O2_gas
!>   HeatProduct = 0. | HeatReactantA = HeatReactantB = Heat_O
!> ProdSpec is recombined species, which is used for selection of linear polyatomic case
!>   DISSOCIATION: CO2_ads + (-)_ads -> CO_ads + O_ads
!>   HeatProduct = Heat_CO2 | HeatReactantA = CO | HeatReactantB = Heat_O  | product -> polyatomic-linear -> E_a = 2*E_a
!>   RECOMBINATION: CO_ads + O_ads -> CO2_gas
!>   HeatProduct = 0. | HeatReactantA = CO | HeatReactantB = Heat_O  | product -> polyatomic-linear -> E_a = 2*E_a
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars     ,ONLY: SpecDSMC, PolyatomMolDSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)               :: HeatProduct, HeatReactantA, HeatReactantB
REAL, INTENT(IN)               :: DProduct
INTEGER,INTENT(IN)             :: ProdSpec
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Delta_H
LOGICAL                        :: Dissociation
!===================================================================================================================================
! decide if dissociation or recombination reaction
Dissociation = .FALSE.
IF (DProduct.GT.0.) Dissociation = .TRUE.

Delta_H = HeatProduct +ABS(DProduct) -HeatReactantA -HeatReactantB
CalcDissRecombActEnergy = 0.5 * ( (HeatReactantA*HeatReactantB / (HeatReactantA+HeatReactantB)) + Delta_H )

IF (Dissociation) THEN
  IF (SpecDSMC(ProdSpec)%PolyatomicMol) THEN
    IF(PolyatomMolDSMC(SpecDSMC(ProdSpec)%SpecToPolyArray)%LinearMolec) THEN
      CalcDissRecombActEnergy = 2.*CalcDissRecombActEnergy
    END IF
  END IF
ELSE
  IF (CalcDissRecombActEnergy.LT.0.) CalcDissRecombActEnergy = 0.
  IF (SpecDSMC(ProdSpec)%PolyatomicMol) THEN
    IF(PolyatomMolDSMC(SpecDSMC(ProdSpec)%SpecToPolyArray)%LinearMolec) THEN
      CalcDissRecombActEnergy = 2.*CalcDissRecombActEnergy - Delta_H
    ELSE
      CalcDissRecombActEnergy = CalcDissRecombActEnergy - Delta_H
    END IF
  ELSE
    CalcDissRecombActEnergy = CalcDissRecombActEnergy - Delta_H
  END IF
END IF
IF (CalcDissRecombActEnergy.LT.0.) CalcDissRecombActEnergy = 0.

END FUNCTION CalcDissRecombActEnergy


#if (PP_TimeDiscMethod==42)
REAL FUNCTION CalcAdsorbReactProb(ReactionCase,ReactNum,PartID,SurfID,NormalVelo,E_Activation,E_Activation_max,loc_ActE,loc_nu)
#else
REAL FUNCTION CalcAdsorbReactProb(ReactionCase,ReactNum,PartID,SurfID,NormalVelo,E_Activation,E_Activation_max)
#endif
!===================================================================================================================================
!> Calculates the Probability for Adsorption with TCE Model
!>   if automatic TST is enabled, then mean probability from rate expression with particle temperature is used
!> 1: molecular adsorption
!> 2: dissociative adsorption
!> 3: eley-rideal reaction
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars           ,ONLY: PlanckConst, BoltzmannConst, PI
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: PartSpecies, Species !, PartState
USE MOD_DSMC_Vars              ,ONLY: DSMC, SpecDSMC, PartStateIntEn, PolyatomMolDSMC
USE MOD_DSMC_Analyze           ,ONLY: CalcTVib, CalcTVibPoly
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption
#if (PP_TimeDiscMethod==42)
USE MOD_SurfaceModel_Vars      ,ONLY: SurfModel
#endif
USE MOD_SurfaceModel_PartFunc  ,ONLY: PartitionFuncActAdsorb, PartitionFuncSurf, PartitionFuncGas
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)          :: ReactionCase
INTEGER, INTENT(IN)          :: ReactNum
INTEGER, INTENT(IN)          :: PartID
INTEGER, INTENT(IN)          :: SurfID
REAL, INTENT(IN)             :: NormalVelo
REAL, INTENT(IN)             :: E_Activation
REAL, INTENT(IN)             :: E_Activation_max
!INTEGER, INTENT(IN),OPTIONAL :: PartnerSpecies
#if (PP_TimeDiscMethod==42)
REAL, INTENT(INOUT),OPTIONAL :: loc_ActE
REAL, INTENT(INOUT),OPTIONAL :: loc_nu
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: EZeroPoint_Educt, Xi_Rot, Xi_Vib, Xi_Total, Norm_Ec, phi_1, phi_2!, PartVelo, MeanNormalVelo
REAL    :: Beta, a_f, b_f, c_f, AdsorptionTemp
INTEGER :: SpecID
INTEGER :: DissocNum, AssocNum
!INTEGER :: iDof, iPolyAtMole
!INTEGER :: iQuant
!REAL    :: RanNum
#if (PP_TimeDiscMethod==42)
INTEGER :: iSampleReact
#endif
!===================================================================================================================================
!IF(ReactionCase.EQ.3.AND. (.NOT.PRESENT(PartnerSpecies)))THEN
!  CALL abort(&
!__STAMP__&
!,"CalcAdsorbReactProb can't be calculated for Eley-Rideal without Partnerspecies")
!END IF
SpecID = PartSpecies(PartID)
#if (PP_TimeDiscMethod==42)
a_f = 0.
b_f = 0.
c_f = 0.
#endif

! set DOF
! Testing if the adsorption particle is an atom or molecule, if molecule: is it polyatomic?
EZeroPoint_Educt = 0.
Xi_Rot = 0
IF(SpecDSMC(SpecID)%InterID.EQ.2) THEN
  IF(SpecDSMC(SpecID)%PolyatomicMol) THEN
    EZeroPoint_Educt = EZeroPoint_Educt + SpecDSMC(SpecID)%EZeroPoint
    ! Calculation of the vibrational degree of freedom for the particle
    IF (PartStateIntEn(1,PartID).GT.SpecDSMC(SpecID)%EZeroPoint) THEN
      Xi_vib = 2.*(PartStateIntEn(1,PartID)-SpecDSMC(SpecID)%EZeroPoint) &
              / (BoltzmannConst*CalcTVibPoly(PartStateIntEn(1,PartID), SpecID))
    ELSE
      Xi_vib = 0.0
    END IF
    IF(PolyatomMolDSMC(SpecDSMC(SpecID)%SpecToPolyArray)%LinearMolec) THEN
      Xi_Rot = 3
    ELSE
      Xi_Rot = 2
    END IF
  ELSE
    EZeroPoint_Educt = EZeroPoint_Educt + DSMC%GammaQuant*BoltzmannConst*SpecDSMC(SpecID)%CharaTVib
    IF((PartStateIntEn(1,PartID)-DSMC%GammaQuant*BoltzmannConst*SpecDSMC(SpecID)%CharaTVib).GT.0.0) THEN
!           IF(ChemReac%MeanEVibQua_PerIter(SpecID).GT.0.0) THEN
      Xi_vib = 2.*(PartStateIntEn(1,PartID)-DSMC%GammaQuant*BoltzmannConst*SpecDSMC(SpecID)%CharaTVib) &
              / (BoltzmannConst*CalcTVib(SpecDSMC(SpecID)%CharaTVib, PartStateIntEn(1,PartID), SpecDSMC(SpecID)%MaxVibQuant))
!             Xi_vib = 2.0*ChemReac%MeanEVibQua_PerIter(SpecID) &
!                     * LOG(1.0/ChemReac%MeanEVibQua_PerIter(SpecID) + 1.0)
    ELSE
      Xi_vib = 0.0
    END IF
    Xi_Rot = 2
  END IF
ELSE
  Xi_vib = 0.0
END IF


CalcAdsorbReactProb = 0.0
IF (Adsorption%TST_Calc(ReactNum,SpecID)) THEN
  !PartVelo = VECNORM(PartState(4:6,PartID))
  !Norm_Ec = PartVelo**2 * 0.5*Species(SpecID)%MassIC !+ PartStateIntEn(2,PartID) + PartStateIntEn(1,PartID) - EZeroPoint_Educt
  !Xi_Total = 1.! Xi_vib + Xi_rot + 3.
  !AdsorptionTemp=2.*Norm_Ec/Xi_Total/BoltzmannConst
  !MeanNormalVelo = SQRT((BoltzmannConst*AdsorptionTemp) / (2*PI*Species(SpecID)%MassIC)) ! equilibrium mean thermal velo for AdsorptionTemp
  AdsorptionTemp = (Species(SpecID)%MassIC*NormalVelo**2 / BoltzmannConst) * 2./PI
#if (PP_TimeDiscMethod==42)
  a_f = (BoltzmannConst*AdsorptionTemp/PlanckConst) &
        *(PartitionFuncActAdsorb(SpecID, AdsorptionTemp)/PartitionFuncGas(SpecID, AdsorptionTemp))
#endif
  SELECT CASE(ReactionCase)
  CASE(1) ! adsorption
    CalcAdsorbReactProb = 1. - EXP(-E_Activation_max/(BoltzmannConst*AdsorptionTemp))
  CASE(2,3) ! dissociation or eley-rideal
    CalcAdsorbReactProb = EXP(-E_activation/(BoltzmannConst*AdsorptionTemp))
  END SELECT
ELSE
  Beta = 0.0
  Norm_Ec = NormalVelo**2 * 0.5*Species(SpecID)%MassIC + PartStateIntEn(2,PartID) + PartStateIntEn(1,PartID) - EZeroPoint_Educt
  Xi_Total = Xi_vib + Xi_rot + 1.
  SELECT CASE(ReactionCase)
  CASE(1) ! adsorption
    IF ((Norm_Ec.GE.E_Activation) .AND. (Norm_Ec.LT.E_Activation_max)) THEN
      a_f = Adsorption%Ads_Prefactor(SpecID)
      b_f = Adsorption%Ads_Powerfactor(SpecID)
      phi_1 = b_f - 1. + Xi_Total/2.
      phi_2 = 1. - Xi_Total/2.
      IF((phi_1+1).GT.0.0) THEN
        c_f = BoltzmannConst/PlanckConst &
            * REAL(Adsorption%DensSurfAtoms(SurfID)*Adsorption%AreaIncrease(SurfID)) &
            / ( (BoltzmannConst / (2*Pi*Species(SpecID)%MassIC))**0.5 )
        Beta = a_f * c_f * BoltzmannConst**(-b_f) * GAMMA(Xi_Total/2.) / (GAMMA(phi_1+1))
      END IF
      CalcAdsorbReactProb = Beta * ((Norm_Ec) - E_Activation)**phi_1 * (Norm_Ec) ** phi_2
    END IF
  CASE(2,3) ! dissociation or eley-rideal
    IF ((Norm_Ec.GE.E_Activation) ) THEN
      SELECT CASE(ReactionCase)
      CASE(2)
        a_f = Adsorption%Diss_Prefactor(DissocNum,SpecID)
        b_f = Adsorption%Diss_Powerfactor(DissocNum,SpecID)
      CASE(3)
        a_f = Adsorption%ER_Prefactor(AssocNum,SpecID)
        b_f = Adsorption%ER_Powerfactor(AssocNum,SpecID)
      END SELECT
      phi_1 = b_f - 1. + Xi_Total/2.
      phi_2 = 1. - Xi_Total/2.
      IF((phi_1+1).GT.0.0) THEN
        c_f = BoltzmannConst/PlanckConst &
            * REAL(Adsorption%DensSurfAtoms(SurfID)*Adsorption%AreaIncrease(SurfID)) &
            / ( (BoltzmannConst / (2*Pi*Species(SpecID)%MassIC))**0.5 )
        Beta = a_f * c_f * BoltzmannConst**(-b_f) * GAMMA(Xi_Total/2.) / (GAMMA(phi_1+1))
      END IF
      CalcAdsorbReactProb = Beta * ((Norm_Ec) - E_Activation)**phi_1 * (Norm_Ec) ** phi_2
    END IF
  END SELECT
END IF

#if (PP_TimeDiscMethod==42)
iSampleReact = 1 + ReactNum
IF ((.NOT.DSMC%ReservoirRateStatistic).AND.(CalcAdsorbReactProb.GT.0.)) THEN
  !IF (calcAdsorbReactProb.GT.1) THEN
  !  SurfModel%ProperInfo(SpecID)%NumAdsReact(iSampleReact) = &
  !      SurfModel%ProperInfo(SpecID)%NumAdsReact(iSampleReact) + 1.
  !ELSE
    SurfModel%ProperInfo(SpecID)%NumAdsReact(iSampleReact) = &
        SurfModel%ProperInfo(SpecID)%NumAdsReact(iSampleReact) + CalcAdsorbReactProb
  !END IF
END IF
SurfModel%ProperInfo(SpecID)%MeanAdsActE(iSampleReact) = &
    SurfModel%ProperInfo(SpecID)%MeanAdsActE(iSampleReact) + (E_Activation - E_Activation_max) /BoltzmannConst
loc_ActE = (E_Activation - E_Activation_max) /BoltzmannConst
IF (Adsorption%TST_Calc(ReactNum,SpecID)) THEN
  SurfModel%ProperInfo(SpecID)%MeanAdsnu(iSampleReact) = &
      SurfModel%ProperInfo(SpecID)%MeanAdsnu(iSampleReact) + a_f
  loc_nu = a_f
ELSE
  AdsorptionTemp=2.*Norm_Ec/Xi_Total/BoltzmannConst
  SurfModel%ProperInfo(SpecID)%MeanAdsnu(iSampleReact) = &
      SurfModel%ProperInfo(SpecID)%MeanAdsnu(iSampleReact) + a_f*c_f*AdsorptionTemp**b_f
  loc_nu = a_f*c_f*AdsorptionTemp**b_f
END IF
SurfModel%ProperInfo(SpecID)%AdsReactCount(iSampleReact) = &
    SurfModel%ProperInfo(SpecID)%AdsReactCount(iSampleReact) + 1
#endif

END FUNCTION CalcAdsorbReactProb


LOGICAL FUNCTION SpaceOccupied(SurfID,subsurfxi,subsurfeta,Coordination,SurfPos)
!===================================================================================================================================
!> Check if particle has enough space on given SurfPos
!>  cycle through all neighbours and check if nearest (valid) neighbour is occupied and blocks considered current position
!===================================================================================================================================
! MODULES
USE MOD_SurfaceModel_Vars ,ONLY: SurfDistInfo, BlockingNeigh
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: SurfID, SubSurfxi, SubSurfeta, Coordination, SurfPos
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iNeigh, NeighCoord
!===================================================================================================================================
SpaceOccupied = .FALSE.
IF ( ANY(BlockingNeigh(Coordination,1:3)) ) THEN
  DO iNeigh = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%nNeighbours
    IF (.NOT.SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%IsNearestNeigh(SurfPos,iNeigh)) CYCLE
    NeighCoord = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%NeighSite(SurfPos,iNeigh)
    IF (NeighCoord.EQ.0) CYCLE ! outside of boundaries of lattice
    IF ( .NOT.BlockingNeigh(Coordination,NeighCoord) ) CYCLE
    ASSOCIATE (NeighSpec => SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(NeighCoord)%Species( &
                            SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%NeighPos(SurfPos,iNeigh)))
      IF ( (NeighSpec.NE.0) ) SpaceOccupied = .TRUE.
    END ASSOCIATE
  END DO
END IF

END FUNCTION SpaceOccupied


SUBROUTINE UpdateSurfPos(SurfID,subsurfxi,subsurfeta,Coordination,SurfPos,Species,removeFlag,relaxation)
!===================================================================================================================================
!> updates bond order for surfpos and species (if removeflag=True then remove species from space else add it)
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_SurfaceModel_Vars      ,ONLY: SurfDistInfo
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound, SurfMesh
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars              ,ONLY: SpecDSMC, DSMC, PolyatomMolDSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)          :: SurfID, SubSurfxi, SubSurfeta, Coordination, SurfPos, Species
LOGICAL, INTENT(IN)          :: removeFlag
LOGICAL, INTENT(IN),OPTIONAL :: relaxation
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iInterAtom, Indx, Indy
INTEGER :: BondOrderAddon, iQuant
REAL    :: iRan, WallTemp
INTEGER :: iPolyatMole, iDOF
!===================================================================================================================================
IF (removeFlag) THEN
  SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%Species(SurfPos) = 0
  BondOrderAddon = -1
  IF (PRESENT(relaxation)) SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%EVib(SurfPos) = 0.0
ELSE
  SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%Species(SurfPos) = Species
  BondOrderAddon = 1
  IF (PRESENT(relaxation)) THEN
    ! set vibrational energy of adsorbate
    IF ((SpecDSMC(Species)%InterID.EQ.2).OR.(SpecDSMC(Species)%InterID.EQ.20)) THEN
      WallTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC(SurfMesh%SurfIDToSideID(SurfID))))
      IF(SpecDSMC(Species)%PolyatomicMol) THEN
        iPolyatMole = SpecDSMC(Species)%SpecToPolyArray
        DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
          CALL RANDOM_NUMBER(iRan)
          iQuant = INT(-LOG(iRan)*WallTemp/PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF))
          DO WHILE (iQuant.GE.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(iDOF))
            CALL RANDOM_NUMBER(iRan)
            iQuant = INT(-LOG(iRan)*WallTemp/PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF))
          END DO
          SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%EVib(SurfPos) = &
              SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%EVib(SurfPos) &
              + (iQuant + DSMC%GammaQuant)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
        END DO
      ELSE
        CALL RANDOM_NUMBER(iRan)
        iQuant = INT(-LOG(iRan)*WallTemp/SpecDSMC(Species)%CharaTVib)
        DO WHILE (iQuant.GE.SpecDSMC(Species)%MaxVibQuant)
          CALL RANDOM_NUMBER(iRan)
          iQuant = INT(-LOG(iRan)*WallTemp/SpecDSMC(Species)%CharaTVib)
        END DO
        SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%EVib(SurfPos) = &
            (iQuant + DSMC%GammaQuant)*SpecDSMC(Species)%CharaTVib*BoltzmannConst
      END IF
    ELSE
      SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%EVib(SurfPos) = 0.0
    END IF
  END IF
END IF

DO iInterAtom = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%nInterAtom
  Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%BondAtomIndx(SurfPos,iInterAtom)
  Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%BondAtomIndy(SurfPos,iInterAtom)
  SurfDistInfo(subsurfxi,subsurfeta,SurfID)%SurfAtomBondOrder(Species,Indx,Indy) = &
      SurfDistInfo(subsurfxi,subsurfeta,SurfID)%SurfAtomBondOrder(Species,Indx,Indy) + BondOrderAddon
END DO

END SUBROUTINE UpdateSurfPos


REAL FUNCTION SampleAdsorptionHeat(SurfID,iSubSurf,jSubSurf)
!===================================================================================================================================
!> Sums up the current heat of adsorption on the specified SMCR surface
!> If SurfID is non catalytic, adsorptionheat is zero
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
USE MOD_SurfaceModel_Vars      ,ONLY: SurfDistInfo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: SurfID, iSubSurf, jSubSurf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: SurfPos, SpecID, AdsorbID, Coord
!===================================================================================================================================
SampleAdsorptionHeat = 0.0
IF (.NOT.IsReactiveSurface(SurfID)) RETURN

ASSOCIATE ( nSites => SurfDistInfo(iSubSurf,jSubSurf,SurfID)%nSites(:) ,&
            nSitesRemain => SurfDistInfo(iSubSurf,jSubSurf,SurfID)%SitesRemain(:) )
  DO Coord = 1,3
    DO AdsorbID = 1,nSites(Coord)-nSitesRemain(Coord)
      Surfpos = SurfDistInfo(iSubSurf,jSubSurf,SurfID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+AdsorbID)
      SpecID = SurfDistInfo(iSubSurf,jSubSurf,SurfID)%AdsMap(Coord)%Species(Surfpos)
      SampleAdsorptionHeat = SampleAdsorptionHeat + Calc_Adsorb_Heat(iSubSurf,jSubSurf,SurfID,SpecID,Surfpos,.FALSE.) &
                           + SurfDistInfo(iSubSurf,jSubSurf,SurfID)%AdsMap(Coord)%EVib(Surfpos)/BoltzmannConst
    END DO
  END DO
END ASSOCIATE

END FUNCTION SampleAdsorptionHeat


LOGICAL FUNCTION IsReactiveSurface(SurfID)
!===================================================================================================================================
!> Checks if SurfID has reactive boundary flag
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound, SurfMesh
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: SurfID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IsReactiveSurface = .FALSE.
IF (PartBound%Reactive(PartBound%MapToPartBC(BC(SurfMesh%SurfIDToSideID(SurfID))))) IsReactiveSurface = .TRUE.
END FUNCTION IsReactiveSurface

INTEGER FUNCTION SurfaceHasModelNum(SurfID)
!===================================================================================================================================
!> Checks if SurfID has reactive boundary flag
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound, SurfMesh
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: SurfID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SurfaceHasModelNum = PartBound%SurfaceModel(PartBound%MapToPartBC(BC(SurfMesh%SurfIDToSideID(SurfID))))
END FUNCTION SurfaceHasModelNum


SUBROUTINE SMCR_AdjustMapNum(subsurfxi,subsurfeta,SurfSideID,adsorbates_num,SpecID,SampleFlag)
!===================================================================================================================================
!> Routine for adjusting the number of Adsorbates for the adsorbate background distribution (wallmodel 3)
!> in case adsorption took place in SMCR_PartAdsorb and Coverage changed sufficiently (depending on particle weighting)
!> if more particles are adsorbed than space left on surface then adsoprtion number is adjusted
!> same for the case where particles are removed from surface
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals_Vars           ,ONLY: PlanckConst, BoltzmannConst
USE MOD_Globals                ,ONLY: abort
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound, SurfMesh, SampWall
USE MOD_Particle_Vars          ,ONLY: Species, WriteMacroSurfaceValues
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption, SurfDistInfo
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_TimeDisc_Vars          ,ONLY: TEnd, time
!----------------------------------------------------------------------------------------------------------------------------------!
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)               :: subsurfxi,subsurfeta,SurfSideID,SpecID
INTEGER,INTENT(INOUT)            :: adsorbates_num
LOGICAL,INTENT(IN),OPTIONAL      :: SampleFlag
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                          :: dist, PartBoundID
INTEGER                          :: Coord, nEmptySites, nAdsorbates, Surfpos, UsedSiteMapPos, nSites, nSitesRemain
INTEGER                          :: ntreated, SitesRemain
REAL                             :: RanNum
LOGICAL                          :: LocSampleFlag
!===================================================================================================================================
IF (PRESENT(SampleFlag)) THEN
  LocSampleFlag = SampleFlag
ELSE
  LocSampleFlag = .FALSE.
END IF
IF (LocSampleFlag) THEN
  IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
    SampWall(SurfSideID)%SurfModelState(5,subsurfxi,subsurfeta) = SampWall(SurfSideID)%SurfModelState(5,subsurfxi,subsurfeta) &
        + (SampleAdsorptionHeat(SurfSideID,subsurfxi,subsurfeta) * BoltzmannConst &
        / REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(3))) &
        * REAL(INT(Adsorption%DensSurfAtoms(SurfSideID) &
        * SurfMesh%SurfaceArea(subsurfxi,subsurfeta,SurfSideID),8)) / Species(1)%MacroParticleFactor
  END IF
END IF

PartBoundID = PartBound%MapToPartBC(BC(SurfMesh%SurfIDToSideID(SurfSideID)))
IF (adsorbates_num.GT.0) THEN
  ! distribute adsorbates randomly on the surface on the correct site and assign surface atom bond order
  ! do this only if chosen surface position is not occupied
  Coord = Adsorption%Coordination(PartBoundID,SpecID)
  SitesRemain = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord)
  ! check if new_adsorbates greater than number of empty adsorbate-sites on surface and correct to be added number
  ! the remaining number of to be added particles is still kept in tmp array
  IF ((SitesRemain - adsorbates_num).LT.0) THEN
    adsorbates_num = SitesRemain
  END IF
  ntreated = 0
  dist = 0
  nEmptySites = SitesRemain
  DO WHILE (ntreated.LT.adsorbates_num .AND. nEmptySites.GT.0)
    CALL RANDOM_NUMBER(RanNum)
    Surfpos = 1 + INT(nEmptySites * RanNum)
    UsedSiteMapPos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(Surfpos)
    IF (.NOT.SpaceOccupied(SurfSideID,SubSurfxi,SubSurfeta,Coord,UsedSiteMapPos)) THEN
      CALL UpdateSurfPos(SurfSideID,SubSurfxi,SubSurfeta,Coord,UsedSiteMapPos,SpecID,.FALSE.,relaxation=.TRUE.)
      ntreated = ntreated + 1
    END IF
    ! rearrange UsedSiteMap-Surfpos-array
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(Surfpos) = &
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nEmptySites)
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nEmptySites) = UsedSiteMapPos
    nEmptySites = nEmptySites - 1
    dist = dist + 1
  END DO
  IF (ntreated.LT.dist) THEN
    ! sort usedsitemap array
    nEmptySites = SitesRemain
    SurfPos = nEmptySites
    DO dist=1,SitesRemain
      UsedSiteMapPos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(Surfpos)
      IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(UsedSiteMapPos).NE.0) THEN
        ! rearrange UsedSiteMap-Surfpos-array
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(Surfpos) = &
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nEmptySites)
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nEmptySites) = UsedSiteMapPos
        nEmptySites = nEmptySites - 1
      END IF
      SurfPos = SurfPos - 1
    END DO
  END IF
  adsorbates_num = ntreated
  SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord) = nEmptySites
ELSE IF (adsorbates_num.LT.0) THEN
  ! remove adsorbates randomly on the surface on the correct site and assign surface atom bond order
  dist = -1
  Coord = Adsorption%Coordination(PartBoundID,SpecID)
  nSites = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Coord)
  nSitesRemain = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord)
  nAdsorbates = nSites - nSitesRemain
  ! check if new_adsorbates lower than number of adsorbates tracked on surface and correct to be removed number
  ! the remaining number of to be removed particles is still kept in tmp array
  IF ((nAdsorbates - ABS(adsorbates_num)).LT.0) THEN
    adsorbates_num = -nAdsorbates
  END IF
  DO WHILE (dist.GE.adsorbates_num)
    CALL RANDOM_NUMBER(RanNum)
    Surfpos = 1 + INT(nAdsorbates * RanNum)
    UsedSiteMapPos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain+Surfpos)
    CALL UpdateSurfPos(SurfSideID,SubSurfxi,SubSurfeta,Coord,UsedSiteMapPos,SpecID,.TRUE.,relaxation=.TRUE.)
    ! rearrange UsedSiteMap-Surfpos-array
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain+Surfpos) = &
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain+1)
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain+1) = UsedSiteMapPos
    nAdsorbates = nAdsorbates - 1
    nSitesRemain = nSitesRemain + 1
    dist = dist - 1
  END DO
  SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord) = nSitesRemain
END IF

IF (LocSampleFlag) THEN
  IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
    SampWall(SurfSideID)%SurfModelState(5,subsurfxi,subsurfeta) = SampWall(SurfSideID)%SurfModelState(5,subsurfxi,subsurfeta) &
        - (SampleAdsorptionHeat(SurfSideID,subsurfxi,subsurfeta) * BoltzmannConst &
        / REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(3))) &
        * REAL(INT(Adsorption%DensSurfAtoms(SurfSideID) &
        * SurfMesh%SurfaceArea(subsurfxi,subsurfeta,SurfSideID),8)) / Species(1)%MacroParticleFactor
  END IF
END IF

END SUBROUTINE SMCR_AdjustMapNum


END MODULE MOD_SurfaceModel_Tools
