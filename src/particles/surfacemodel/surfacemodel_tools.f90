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
INTERFACE CalcAdsorbProb
  MODULE PROCEDURE CalcAdsorbProb
END INTERFACE

INTERFACE CalcDesorbProb
  MODULE PROCEDURE CalcDesorbProb
END INTERFACE

INTERFACE Calc_Adsorb_Heat
  MODULE PROCEDURE Calc_Adsorb_Heat
END INTERFACE

INTERFACE Calc_E_Act
  MODULE PROCEDURE Calc_E_Act
END INTERFACE

INTERFACE Set_TST_Factors
  MODULE PROCEDURE Set_TST_Factors
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

INTERFACE SMCR_AdjustMapNum
  MODULE PROCEDURE SMCR_AdjustMapNum
END INTERFACE

PUBLIC :: CalcAdsorbProb
PUBLIC :: CalcDesorbProb
PUBLIC :: Calc_Adsorb_Heat
PUBLIC :: Calc_E_Act
PUBLIC :: Set_TST_Factors
PUBLIC :: CalcAdsorbReactProb
PUBLIC :: SpaceOccupied
PUBLIC :: UpdateSurfPos
PUBLIC :: SMCR_AdjustMapNum
!===================================================================================================================================

CONTAINS

SUBROUTINE CalcAdsorbProb()
!===================================================================================================================================
!> Calculcation of adsorption probability for different model (wallmodel 1 and 2)
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Particle_Vars          ,ONLY: nSpecies, PartSurfaceModel
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, PartBound
#if (PP_TimeDiscMethod==42)
USE MOD_DSMC_Vars              ,ONLY: DSMC
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
  DO SurfSide=1,SurfMesh%nSides
    PartBoundID = PartBound%MapToPartBC(BC(Adsorption%SurfSideToGlobSideMap(SurfSide)))
    DO q = 1,nSurfSample
      DO p = 1,nSurfSample
!----------------------------------------------------------------------------------------------------------------------------------!
        IF (PartSurfaceModel.EQ.1) THEN
!----------------------------------------------------------------------------------------------------------------------------------!
!   Kisluik Sticking Model from Kolasinski's Surface Science (book)
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
        ELSE IF (PartSurfaceModel.EQ.2) THEN
!----------------------------------------------------------------------------------------------------------------------------------!
! Recombination Model described by Laux
          Adsorption%ProbAds(p,q,SurfSide,iSpec) = Adsorption%RecombCoeff(PartBoundID,iSpec)-Adsorption%ProbDes(p,q,SurfSide,iSpec)
        END IF
!----------------------------------------------------------------------------------------------------------------------------------!
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          Adsorption%AdsorpInfo(iSpec)%MeanProbAds = Adsorption%AdsorpInfo(iSpec)%MeanProbAds+Adsorption%ProbAds(p,q,SurfSide,iSpec)
        END IF
#endif
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
USE MOD_Particle_Vars          ,ONLY: nSpecies, PartSurfaceModel
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, PartBound
USE MOD_TimeDisc_Vars          ,ONLY: dt
#if (PP_TimeDiscMethod==42)
USE MOD_TimeDisc_Vars          ,ONLY: iter
USE MOD_DSMC_Vars              ,ONLY: DSMC
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! Local variable declaration
INTEGER                          :: SurfSide, iSpec, p, q
REAL                             :: Theta, nu_des, rate, WallTemp
REAL                             :: E_des
INTEGER                          :: PartBoundID
!===================================================================================================================================
! CALL CalcSurfDistInteraction()
DO SurfSide=1,SurfMesh%nSides
  PartBoundID = PartBound%MapToPartBC(BC(Adsorption%SurfSideToGlobSideMap(SurfSide)))
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
!----------------------------------------------------------------------------------------------------------------------------------!
        IF (PartSurfaceModel.EQ.1) THEN
!----------------------------------------------------------------------------------------------------------------------------------!
!   Polanyi-Wigner-eq. from Kolasinski's Surface Science (book)
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
#if (PP_TimeDiscMethod==42)
          IF (.NOT.DSMC%ReservoirRateStatistic) THEN
            Adsorption%AdsorpInfo(iSpec)%MeanProbDes = Adsorption%AdsorpInfo(iSpec)%MeanProbDes &
                                                     + Adsorption%ProbDes(p,q,SurfSide,iSpec)
          END IF
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
        ELSE IF (PartSurfaceModel.EQ.2) THEN
!----------------------------------------------------------------------------------------------------------------------------------!
! Recombination Model described by Laux
!----------------------------------------------------------------------------------------------------------------------------------!
          IF (Adsorption%RecombData(1,iSpec).LE.0) THEN
            Adsorption%ProbDes(p,q,SurfSide,iSpec) = 0.
          ELSE
            IF (Adsorption%Coverage(p,q,SurfSide,Adsorption%RecombData(1,iSpec)).LE.0) THEN
              Adsorption%ProbDes(p,q,SurfSide,iSpec) = 0.
            ELSE
              Adsorption%ProbDes(p,q,SurfSide,iSpec) = Adsorption%RecombCoeff(PartBoundID,iSpec) &
                  * ( 1 - exp( - Adsorption%Coverage(p,q,SurfSide, Adsorption%RecombData(1,iSpec) ) ) )
            END IF
          END IF
        END IF ! PartSurfaceModel
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
USE MOD_Particle_Boundary_vars ,ONLY: PartBound
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
REAL , ALLOCATABLE             :: x(:)!, D_AL(:), delta(:)
INTEGER , ALLOCATABLE          :: m(:)!, Neigh_bondorder(:)
INTEGER                        :: bondorder
REAL                           :: D_AB, D_AX, D_BX
REAL                           :: Heat_A, Heat_B
REAL                           :: A, B, sigma, sigma_m
!REAL                           :: Heat_D_AL
!INTEGER                        :: neighSpec, neighSpec2, Coord2, Coord3, ReactNum, nNeigh_interactions
!===================================================================================================================================
PartBoundID = PartBound%MapToPartBC(BC(Adsorption%SurfSideToGlobSideMap(SurfSideID)))
Coordination = Adsorption%Coordination(PartBoundID,Species)
ALLOCATE( x(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom) )
!   ALLOCATE( z(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom) )
ALLOCATE( m(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom) )
x(:) = 1. ! averaged bond-index for surface atoms
m(:) = 1  ! number of adsorbates belonging to surface atom
Calc_Adsorb_Heat = 0.
sigma = 0.
sigma_m = 0.
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
    x(j) = 1.
  END DO
END IF
! calculate local scaling factor for chosen surface site
DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom
!     x(j) = x(j) / REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom)
!     sigma = sigma + (2.*x(j) - x(j)**2.) * (2.*(1./REAL(m(j))) - (1./REAL(m(j)))**2.)
  sigma_m = sigma_m + (2.*(1./REAL(m(j))) - (1./REAL(m(j)))**2.) &
                  / REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom)
END DO
IF (Coordination.EQ.1) THEN
  sigma = (2 - 1. / REAL(Adsorption%CrystalIndx(SurfSideID)) )
ELSE
  sigma = (2 - 1. / REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom) )
END IF

! Testing if the adsorption particle is an atom or molecule, if molecule: is it polyatomic?
! and calculate right heat of adsorption to surface atoms
IF(SpecDSMC(Species)%InterID.EQ.2) THEN
  IF(SpecDSMC(Species)%PolyatomicMol) THEN
    D_AB = Adsorption%EDissBond(0,Species)
    SELECT CASE(Coordination)
    CASE(1) ! hollow (radical with localized electron like NH2)
      Heat_A = Adsorption%HeatOfAdsZero(PartBoundID,Species)
      Calc_Adsorb_Heat = Heat_A**2/(D_AB+Heat_A/REAL(Adsorption%CrystalIndx(SurfSideID))) * sigma_m
    CASE(2) ! bridge
      IF (Adsorption%DiCoord(PartBoundID,Species).EQ.1) THEN ! dicoordination (HCOOH --> M--(HC)O-O(H)--M) (M--O bond)
        D_AX = Adsorption%EDissBondAdsorbPoly(0,Species) ! Bond HC--O
        Heat_A = Adsorption%HeatOfAdsZero(PartBoundID,Species) * (2.*(1./REAL(m(1))) - (1./REAL(m(1)))**2.)
        A = Heat_A**2./(D_AX+Heat_A)
        D_BX = Adsorption%EDissBondAdsorbPoly(1,Species) ! Bond O--H
        Heat_B = Adsorption%HeatOfAdsZero(PartBoundID,Species) * (2.*(1./REAL(m(2))) - (1./REAL(m(2)))**2.)
        B = Heat_B**2./(D_BX+Heat_B)
        Calc_Adsorb_Heat = ( A*B*( A + B ) + D_AB*( A - B )**2. ) / ( A*B + D_AB*( A + B ) ) * sigma_m

!        Heat_A = Adsorption%HeatOfAdsZero(PartBoundID,Species) * (2.*x(1) - x(1)**2)
!        Heat_M = Adsorption%HeatOfAdsZeroM(Species) * (2.*x(1) - x(1)**2)
!        A = Heat_A * (1- ( (m*Heat_B/( m*Heat_A + Heat_B ))**2) )
!
!        Heat_A = Adsorption%HeatOfAdsZero2(Species) * (2.*x(2) - x(2)**2)
!        Heat_M = Adsorption%HeatOfAdsZeroR(Species) * (2.*x(2) - x(2)**2)
!        B = Heat_A * (1- ( (mtilde*Heat_B/( mtilde*Heat_A + Heat_B ))**2) )

      ELSE IF (Adsorption%DiCoord(PartBoundID,Species).EQ.2) THEN ! chelate binding (NO2 --> M--O-N-O--M)
        D_AX = Adsorption%EDissBondAdsorbPoly(0,Species) ! Bond O--N
        Heat_A = Adsorption%HeatOfAdsZero(PartBoundID,Species)
        Heat_A = Heat_A**2/(D_AX+Heat_A)
        D_BX = Adsorption%EDissBondAdsorbPoly(1,Species) ! Bond N--O
        Heat_B = Adsorption%HeatOfAdsZero(PartBoundID,Species)
        Heat_B = Heat_B**2/(D_BX+Heat_B)
        A = Heat_A**2. * ( Heat_A + 2.*Heat_B ) / ( Heat_A + Heat_B )**2.
        B = Heat_B**2. * ( Heat_B + 2.*Heat_A ) / ( Heat_A + Heat_B )**2.
        Calc_Adsorb_Heat = (A + B) * sigma_m
      ELSE IF (Adsorption%DiCoord(PartBoundID,Species).EQ.3) THEN !weak binding
        Heat_A = Adsorption%HeatOfAdsZero(PartBoundID,Species)
        Calc_Adsorb_Heat = Heat_A**2./(D_AB+Heat_A/2.) * sigma_m
      ELSE
        Heat_A = Adsorption%HeatOfAdsZero(PartBoundID,Species) * sigma * sigma_m
        Calc_Adsorb_Heat = Heat_A**2./(D_AB+Heat_A)
      END IF
    CASE(3) ! on-top (closed shell or open shell with unlocalized electron like CO)
      Heat_A = Adsorption%HeatOfAdsZero(PartBoundID,Species)
      Calc_Adsorb_Heat = Heat_A**2./(D_AB+Heat_A) * sigma_m
    END SELECT
  ELSE
    D_AB = Adsorption%EDissBond(0,Species)
    SELECT CASE(Coordination)
    CASE(1) ! hollow (radical with localized electron like C-H)
      Heat_A = Adsorption%HeatOfAdsZero(PartBoundID,Species)
      Calc_Adsorb_Heat = (Heat_A**2)/(D_AB+Heat_A/Adsorption%CrystalIndx(SurfSideID)) * sigma_m
    CASE(2) ! bridge (closed shell like O2)
      IF (Adsorption%DiCoord(PartBoundID,Species).EQ.1) THEN
        Heat_A = Adsorption%HeatOfAdsZero(PartBoundID,Species)
        Heat_B = Adsorption%HeatOfAdsZero(PartBoundID,Species)
        A = Heat_A**2. * ( Heat_A + 2.*Heat_B ) / ( Heat_A + Heat_B )**2.
        B = Heat_B**2. * ( Heat_B + 2.*Heat_A ) / ( Heat_A + Heat_B )**2.
        Calc_Adsorb_Heat = ( A*B*( A + B ) + D_AB*( A - B )**2. ) / ( A*B + D_AB*( A + B ) ) * sigma_m
      ELSE IF (Adsorption%DiCoord(PartBoundID,Species).EQ.3) THEN !weak binding
        Heat_A = Adsorption%HeatOfAdsZero(PartBoundID,Species)
        Calc_Adsorb_Heat = Heat_A**2./(D_AB+Heat_A/2.) * sigma_m
      ELSE
        Heat_A = Adsorption%HeatOfAdsZero(PartBoundID,Species) * sigma * sigma_m
        Calc_Adsorb_Heat = Heat_A**2./(D_AB+Heat_A)
      END IF
    CASE(3) ! on-top (closed shell or open shell with unlocalized electron like CO)
      Heat_A = Adsorption%HeatOfAdsZero(PartBoundID,Species)
      Calc_Adsorb_Heat = Heat_A**2./(D_AB+Heat_A) * sigma * sigma_m
    END SELECT
  END IF
ELSE
  Calc_Adsorb_Heat = Adsorption%HeatOfAdsZero(PartBoundID,Species) * sigma * sigma_m
END IF

! routine wird nicht benutz, da höheres Rauschen und größerer Rechenaufwand aber aufgehoben für spätere Einsicht
!   ! calculate additional heat of adsorption for direct interaction (attraction of associating adsorbates)
!   IF ((Adsorption%ReactNum.GT.0) .AND. (Surfpos.GT.0)) THEN
!     ALLOCATE(D_AL(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nNeighbours),&
!              delta(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nNeighbours),&
!              Neigh_bondorder(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nNeighbours))
!     D_AL(:) = 0.
!     delta(:) = 0.
!     Neigh_bondorder(:) = 0
!     nNeigh_interactions = 0
!     ! define dissociation bond energies of neighbours and count interacting neighbours
!     DO l = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nNeighbours
!       Coord2 = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%NeighSite(Surfpos,l)
!       IF (Coord2.GT.0) THEN
!         NeighPos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%NeighPos(Surfpos,l)
!         neighSpec = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%Species(NeighPos)
!         IF ( (neighSpec.NE.0) ) THEN
!           DO ReactNum = 1,(Adsorption%ReactNum-Adsorption%DissNum)
!             IF ( neighSpec.EQ.Adsorption%AssocReact(1,ReactNum,Species) .AND. &
!                  (Adsorption%AssocReact(2,ReactNum,Species).NE.0)) THEN
!               D_AL(l) = Adsorption%EDissBond((Adsorption%DissNum+ReactNum),Species)
!               nNeigh_interactions = nNeigh_interactions + 1
!               CYCLE
!             END IF
!           END DO
!           ! associate bondorder of neighbours
!           DO k = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%nNeighbours
!             Coord3 = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%NeighSite(NeighPos,k)
!             IF (Coord3.GT.0) THEN
!               neighSpec2 = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord3)%Species( &
!                       SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%NeighPos(NeighPos,k))
!               IF ( (neighSpec2.NE.0) ) THEN
!                 DO ReactNum = 1,(Adsorption%ReactNum-Adsorption%DissNum)
!                   IF ( (neighSpec2.EQ.Adsorption%AssocReact(1,ReactNum,neighSpec)) .AND. &
!                        (Adsorption%AssocReact(2,ReactNum,neighSpec).NE.0)) THEN
!                     Neigh_bondorder(l) = Neigh_bondorder(l) + 1
!                     CYCLE
!                   END IF
!                 END DO
!               END IF
!             END IF
!           END DO
!         END IF
!       END IF
!     END DO
!     ! calculate interaction energy between adsorbate and neighbours
!     IF (nNeigh_interactions.NE.0) THEN
!       Heat_D_AL = 0.
!       DO l = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nNeighbours
!         IF (Neigh_bondorder(l).EQ.0) THEN
!           delta(l) = 0.
!         ELSE
!           delta(l) = 1 / REAL(Neigh_bondorder(l))
!         END IF
!         Heat_D_AL = Heat_D_AL + 0.5*D_AL(l) * (2*delta(l)-delta(l)**2)
!       END DO
!     nInterAtom = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom
!     Calc_Adsorb_Heat = (Calc_Adsorb_Heat*nInterAtom + Heat_D_AL*nNeigh_interactions) /(nInterAtom + nNeigh_interactions)
!     END IF
!
!     DEALLOCATE(D_AL,delta,Neigh_bondorder)
!   END IF

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


SUBROUTINE Set_TST_Factors(ReactionCase,a_f,b_f,PartID,ReactNum)!,PartBoundID)
!===================================================================================================================================
!> set the Transition state theory (TST) factors for surface chemistry
!===================================================================================================================================
! MODULES
!USE MOD_Basis         ,ONLY: GetInverse
USE MOD_Globals_Vars      ,ONLY: PlanckConst, BoltzmannConst
USE MOD_SurfaceModel_Vars ,ONLY: Adsorption
USE MOD_Particle_Vars     ,ONLY: PEM, PartSpecies
!===================================================================================================================================
IMPLICIT NONE
!===================================================================================================================================
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: ReactionCase, ReactNum
INTEGER, INTENT(IN)           :: PartID!, PartBoundID
!===================================================================================================================================
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: a_f, b_f
!===================================================================================================================================
! LOCAL VARIABLES
INTEGER                       :: SpecID, ElemID
!REAL                          :: GasTemp
!===================================================================================================================================
SpecID = PartSpecies(PartID)
ElemID = PEM%Element(PartID)
!GasTemp = Adsorption%PartitionTemp(ElemID,SpecID)
SELECT CASE(ReactionCase)
Case(1)
  IF (Adsorption%TST_Calc(0,SpecID)) THEN
    a_f = 0.0
    b_f = 0.0
  ELSE
    a_f = Adsorption%Ads_Prefactor(SpecID)
    b_f = Adsorption%Ads_Powerfactor(SpecID)
  END IF
CASE(2)
  IF (Adsorption%TST_Calc(ReactNum,SpecID)) THEN
    a_f = 0.0
    b_f = 0.0
  ELSE
    a_f = Adsorption%Diss_Prefactor(ReactNum,SpecID)
    b_f = Adsorption%Diss_Powerfactor(ReactNum,SpecID)
  END IF
CASE(3)
  IF (Adsorption%TST_Calc(Adsorption%DissNum+ReactNum,SpecID)) THEN
    a_f = 0.0
    b_f = 0.0
  ELSE
    a_f = Adsorption%ER_Prefactor(ReactNum,SpecID)
    b_f = Adsorption%ER_Powerfactor(ReactNum,SpecID)
  END IF
END SELECT

END SUBROUTINE Set_TST_Factors


REAL FUNCTION CalcAdsorbReactProb(ReactionCase,PartID,NormalVelo,E_Activation,E_Activation_max,a_f,b_f,c_f &
                                 ,PartnerSpecies,CharaTemp,WallTemp)
!===================================================================================================================================
!> Calculates the Probability for Adsorption with TCE Model
!> 1: molecular adsorption
!> 2: dissociative adsorption
!> 3: eley-rideal reaction
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars   ,ONLY: PlanckConst, BoltzmannConst
USE MOD_Globals
USE MOD_Particle_Vars  ,ONLY: PartSpecies, Species !,PartState
USE MOD_DSMC_Vars      ,ONLY: DSMC, SpecDSMC, PartStateIntEn, PolyatomMolDSMC
USE MOD_DSMC_Analyze   ,ONLY: CalcTVib, CalcTVibPoly
!USE MOD_DSMC_ChemReact ,ONLY: gammainc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)          :: ReactionCase
INTEGER, INTENT(IN)          :: PartID
REAL, INTENT(IN)             :: NormalVelo, E_Activation, E_Activation_max
REAL, INTENT(IN)             :: a_f, b_f, c_f
INTEGER, INTENT(IN),OPTIONAL :: PartnerSpecies
REAL, INTENT(IN),OPTIONAL    :: CharaTemp, WallTemp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: EZeroPoint_Educt, Xi_Rot, Xi_Vib, Xi_Total, Norm_Ec, phi_1, phi_2
REAL                           :: SurfPartIntE, SurfPartVibE, RanNum, Beta!, Velocity
INTEGER                        :: iSpec, iQuant, iPolyAtMole, iDOF
!===================================================================================================================================
IF(ReactionCase.EQ.3.AND. (.NOT.PRESENT(PartnerSpecies)))THEN
  CALL abort(&
__STAMP__&
,"CalcAdsorbReactProb can't be calculated for Eley-Rideal without Partnerspecies")
END IF
iSpec = PartSpecies(PartID)
!Velocity = SQRT(PartState(PartID,4)**2+PartState(PartID,5)**2+PartState(PartID,6)**2)
! Testing if the adsorption particle is an atom or molecule, if molecule: is it polyatomic?
EZeroPoint_Educt = 0.
Xi_Rot = 0
IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
  IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
    EZeroPoint_Educt = EZeroPoint_Educt + SpecDSMC(iSpec)%EZeroPoint
    ! Calculation of the vibrational degree of freedom for the particle
    IF (PartStateIntEn(PartID,1).GT.SpecDSMC(iSpec)%EZeroPoint) THEN
      Xi_vib = 2.*(PartStateIntEn(PartID,1)-SpecDSMC(iSpec)%EZeroPoint) &
              / (BoltzmannConst*CalcTVibPoly(PartStateIntEn(PartID,1), iSpec))
    ELSE
      Xi_vib = 0.0
    END IF
    IF(PolyatomMolDSMC(SpecDSMC(iSpec)%SpecToPolyArray)%LinearMolec) THEN
      Xi_Rot = 3
    ELSE
      Xi_Rot = 2
    END IF
  ELSE
    EZeroPoint_Educt = EZeroPoint_Educt + DSMC%GammaQuant*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib
    IF((PartStateIntEn(PartID,1)-DSMC%GammaQuant*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib).GT.0.0) THEN
!           IF(ChemReac%MeanEVibQua_PerIter(iSpec).GT.0.0) THEN
      Xi_vib = 2.*(PartStateIntEn(PartID,1)-DSMC%GammaQuant*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib) &
              / (BoltzmannConst*CalcTVib(SpecDSMC(iSpec)%CharaTVib, PartStateIntEn(PartID,1), SpecDSMC(iSpec)%MaxVibQuant))
!             Xi_vib = 2.0*ChemReac%MeanEVibQua_PerIter(iSpec) &
!                     * LOG(1.0/ChemReac%MeanEVibQua_PerIter(iSpec) + 1.0)
    ELSE
      Xi_vib = 0.0
    END IF
    Xi_Rot = 2
  END IF
ELSE
  Xi_vib = 0.0
END IF

CalcAdsorbReactProb = 0.0
Beta = 0.0
!-----------------------------------------------------------------------------------------------------------------------------------
SELECT CASE(ReactionCase)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(1) ! adsorption
!-----------------------------------------------------------------------------------------------------------------------------------
  Norm_Ec = NormalVelo**2. * 0.5*Species(iSpec)%MassIC + PartStateIntEn(PartID,2) + PartStateIntEn(PartID,1) - EZeroPoint_Educt
  IF ((Norm_Ec.GE.E_Activation) .AND. (Norm_Ec.LT.E_Activation_max)) THEN
    Xi_Total = Xi_vib + Xi_rot + 1.
    phi_1 = b_f - 1. + Xi_Total/2.
    phi_2 = 1. - Xi_Total/2.
    IF((phi_1+1).GT.0.0) THEN
      Beta = a_f * c_f * BoltzmannConst**(-b_f) * GAMMA(Xi_Total/2.) / (GAMMA(phi_1+1))
    END IF
    CalcAdsorbReactProb = Beta * ((Norm_Ec) - E_Activation)**phi_1 * (Norm_Ec) ** phi_2
  END IF
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(2) ! dissociation
!-----------------------------------------------------------------------------------------------------------------------------------
  Norm_Ec = NormalVelo**2 * 0.5*Species(iSpec)%MassIC + PartStateIntEn(PartID,2) + PartStateIntEn(PartID,1) - EZeroPoint_Educt
  IF ((Norm_Ec.GE.E_Activation) ) THEN
    Xi_Total = Xi_vib + Xi_rot + 1.
    phi_1 = b_f - 1. + Xi_Total/2.
    phi_2 = 1. - Xi_Total/2.
    IF((phi_1+1).GT.0.0) THEN
      Beta = a_f * c_f * BoltzmannConst**(-b_f) * GAMMA(Xi_Total/2.) / (GAMMA(phi_1+1))
    END IF
    CalcAdsorbReactProb = Beta * ((Norm_Ec) - E_Activation)**phi_1 * (Norm_Ec) ** phi_2
  END IF
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(3) ! eley-rideal
!-----------------------------------------------------------------------------------------------------------------------------------
  EZeroPoint_Educt = EZeroPoint_Educt + DSMC%GammaQuant*BoltzmannConst*CharaTemp
  SurfPartIntE = 0.
  SurfPartVibE = 0.

  ! Set surface2particle vibrational energy
  CALL RANDOM_NUMBER(RanNum)
  iQuant = INT(-LOG(RanNum)*WallTemp/CharaTemp)
  DO WHILE (iQuant.GE.200)
    CALL RANDOM_NUMBER(RanNum)
    iQuant = INT(-LOG(RanNum)*WallTemp/CharaTemp)
  END DO
  SurfPartIntE = SurfPartIntE + (iQuant + DSMC%GammaQuant)*CharaTemp*BoltzmannConst!*Adsorption%CrystalIndx(SurfSideID)

  ! set vibrational energy of particle
  IF(SpecDSMC(PartnerSpecies)%InterID.EQ.2) THEN
    IF(SpecDSMC(PartnerSpecies)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(PartnerSpecies)%SpecToPolyArray
      DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
        CALL RANDOM_NUMBER(RanNum)
        iQuant = INT(-LOG(RanNum)*WallTemp/PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF))
        DO WHILE (iQuant.GE.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(iDOF))
          CALL RANDOM_NUMBER(RanNum)
          iQuant = INT(-LOG(RanNum)*WallTemp/PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF))
        END DO
        SurfPartVibE = SurfPartVibE + (iQuant + DSMC%GammaQuant)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
      END DO
    ELSE
      CALL RANDOM_NUMBER(RanNum)
      iQuant = INT(-LOG(RanNum)*WallTemp/SpecDSMC(PartnerSpecies)%CharaTVib)
      DO WHILE (iQuant.GE.SpecDSMC(PartnerSpecies)%MaxVibQuant)
        CALL RANDOM_NUMBER(RanNum)
        iQuant = INT(-LOG(RanNum)*Walltemp/SpecDSMC(PartnerSpecies)%CharaTVib)
      END DO
      SurfPartVibE = SurfPartVibE + (iQuant + DSMC%GammaQuant)*SpecDSMC(PartnerSpecies)%CharaTVib*BoltzmannConst
    END IF
  END IF

  IF(SpecDSMC(PartnerSpecies)%InterID.EQ.2) THEN
    IF(SpecDSMC(PartnerSpecies)%PolyatomicMol) THEN
      EZeroPoint_Educt = EZeroPoint_Educt + SpecDSMC(PartnerSpecies)%EZeroPoint
      ! Calculation of the vibrational degree of freedom for the particle
      IF (SurfPartVibE.GT.SpecDSMC(PartnerSpecies)%EZeroPoint) THEN
        Xi_vib = Xi_vib + 2.*(SurfPartVibE-SpecDSMC(PartnerSpecies)%EZeroPoint) &
                / (BoltzmannConst*CalcTVibPoly(SurfPartVibE, PartnerSpecies))
      END IF
    ELSE
      EZeroPoint_Educt = EZeroPoint_Educt + DSMC%GammaQuant*BoltzmannConst*SpecDSMC(PartnerSpecies)%CharaTVib
      IF((SurfPartVibE-DSMC%GammaQuant*BoltzmannConst*SpecDSMC(PartnerSpecies)%CharaTVib).GT.0.0) THEN
        Xi_vib = 2.*(SurfPartVibE-DSMC%GammaQuant*BoltzmannConst*SpecDSMC(PartnerSpecies)%CharaTVib) &
                / (BoltzmannConst*CalcTVib(SpecDSMC(PartnerSpecies)%CharaTVib, SurfPartIntE, SpecDSMC(PartnerSpecies)%MaxVibQuant))
      END IF
    END IF
  END IF
  Norm_Ec = NormalVelo**2. * 0.5*Species(iSpec)%MassIC + PartStateIntEn(PartID,2) + PartStateIntEn(PartID,1) - EZeroPoint_Educt &
          + SurfPartIntE + SurfPartVibE
  IF ((Norm_Ec.GE.E_Activation) ) THEN
    Xi_Total = Xi_vib + Xi_rot + 2.
    phi_1 = b_f - 1. + Xi_Total/2.
    phi_2 = 1. - Xi_Total/2.
    IF((phi_1+1).GT.0.0) THEN
      !Beta = a_f * c_f * BoltzmannConst**(-b_f) * GAMMA(Xi_Total/2.) / (GAMMA(phi_1+1)*((gammainc((/phi_1+1,E_Activation/)))))
      Beta = a_f * c_f * BoltzmannConst**(-b_f) * GAMMA(Xi_Total/2.) / (GAMMA(phi_1+1))
    END IF
    CalcAdsorbReactProb = Beta * ((Norm_Ec) - E_Activation)**phi_1 * (Norm_Ec) ** phi_2
  END IF
!-----------------------------------------------------------------------------------------------------------------------------------
END SELECT
!-----------------------------------------------------------------------------------------------------------------------------------

END FUNCTION CalcAdsorbReactProb


LOGICAL FUNCTION SpaceOccupied(SurfID,subsurfxi,subsurfeta,Coordination,SurfPos)
!===================================================================================================================================
!> Check if particle has enough space on given SurfPos
!===================================================================================================================================
! MODULES
USE MOD_SurfaceModel_Vars ,ONLY: SurfDistInfo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: SurfID, SubSurfxi, SubSurfeta, Coordination, SurfPos
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iNeigh, iCoord
!===================================================================================================================================
SpaceOccupied = .FALSE.
DO iNeigh = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%nNeighbours
  IF (.NOT.SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%IsNearestNeigh(SurfPos,iNeigh)) CYCLE
  IF ((Coordination.EQ.1) .AND. (SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%NeighSite(SurfPos,iNeigh).EQ.1)) &
      CYCLE
  DO iCoord = 1,3
    IF (iCoord.EQ.SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%NeighSite(SurfPos,iNeigh)) THEN
      IF ( (SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(iCoord)%Species( &
            SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%NeighPos(SurfPos,iNeigh)) &
            .NE.0) ) THEN
            SpaceOccupied = .TRUE.
      END IF
    END IF
  END DO
END DO

END FUNCTION SpaceOccupied


SUBROUTINE UpdateSurfPos(SurfID,subsurfxi,subsurfeta,Coordination,SurfPos,Species,removeFlag)
!===================================================================================================================================
!> updates bond order for surfpos and species (if removeflag=True then remove species from space else add it)
!===================================================================================================================================
! MODULES
USE MOD_SurfaceModel_Vars ,ONLY: SurfDistInfo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: SurfID, SubSurfxi, SubSurfeta, Coordination, SurfPos, Species
LOGICAL, INTENT(IN) :: removeFlag
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iInterAtom, Indx, Indy
INTEGER :: BondOrderAddon
!===================================================================================================================================
IF (removeFlag) THEN
  SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%Species(SurfPos) = 0
  BondOrderAddon = -1
ELSE
  SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%Species(SurfPos) = Species
  BondOrderAddon = 1
END IF

DO iInterAtom = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%nInterAtom
  Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%BondAtomIndx(SurfPos,iInterAtom)
  Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coordination)%BondAtomIndy(SurfPos,iInterAtom)
  SurfDistInfo(subsurfxi,subsurfeta,SurfID)%SurfAtomBondOrder(Species,Indx,Indy) = &
      SurfDistInfo(subsurfxi,subsurfeta,SurfID)%SurfAtomBondOrder(Species,Indx,Indy) + BondOrderAddon
END DO

END SUBROUTINE UpdateSurfPos


SUBROUTINE SMCR_AdjustMapNum(subsurfxi,subsurfeta,SurfSideID,adsorbates_num,SpecID)
!===================================================================================================================================
!> Routine for adjusting the number of Adsorbates of adsorbates background distribution (wallmodel 3)
!> if adsorption took place in SMCR_PartAdsorb and Coverage changed sufficiently (depending on particle weighting)
!> if more particles are adsorbed than space left on surface then adsoprtion number is adjusted
!> same for case if too many particles are removed from surface
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals_Vars           ,ONLY: PlanckConst, BoltzmannConst
USE MOD_Globals                ,ONLY: abort
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption, SurfDistInfo
!----------------------------------------------------------------------------------------------------------------------------------!
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)               :: subsurfxi,subsurfeta,SurfSideID,SpecID
INTEGER,INTENT(INOUT)            :: adsorbates_num
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                          :: dist, PartBoundID
INTEGER                          :: Coord, Surfnum, Surfpos, UsedSiteMapPos, nSites, nSitesRemain
INTEGER                          :: iInterAtom, xpos, ypos
REAL                             :: RanNum
!===================================================================================================================================
PartBoundID = PartBound%MapToPartBC(BC(Adsorption%SurfSideToGlobSideMap(SurfSideID)))
IF (adsorbates_num.GT.0) THEN
  ! distribute adsorbates randomly on the surface on the correct site and assign surface atom bond order
  dist = 1
  Coord = Adsorption%Coordination(PartBoundID,SpecID)
  Surfnum = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord)
  ! check if new_adsorbates greater than number of empty adsorbate-sites on surface and correct to be added number
  ! the remaining number of to be added particles is still kept in tmp array
  IF ((SurfNum - adsorbates_num).LT.0) THEN
!    CALL abort(&
!__STAMP__&
!,'Error in AdjustReconstructMapNum: Too many new Adsorbates! not enough Sites for Coordination:' &
!,Adsorption%Coordination(PartBoundID,SpecID))
    adsorbates_num = SurfNum
  END IF
  DO WHILE (dist.LE.adsorbates_num)
    CALL RANDOM_NUMBER(RanNum)
    Surfpos = 1 + INT(Surfnum * RanNum)
    UsedSiteMapPos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(Surfpos)
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(UsedSiteMapPos) = SpecID
    ! assign bond order of respective surface atoms in the surfacelattice
    DO iInterAtom = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
      xpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(UsedSiteMapPos,iInterAtom)
      ypos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(UsedSiteMapPos,iInterAtom)
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(SpecID,xpos,ypos) = &
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(SpecID,xpos,ypos) + 1
    END DO
    ! rearrange UsedSiteMap-Surfpos-array
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(Surfpos) = &
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(Surfnum)
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(Surfnum) = UsedSiteMapPos
    Surfnum = Surfnum - 1
    dist = dist + 1
  END DO
  SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord) = Surfnum
ELSE IF (adsorbates_num.LT.0) THEN
  ! remove adsorbates randomly on the surface on the correct site and assign surface atom bond order
  dist = -1
  Coord = Adsorption%Coordination(PartBoundID,SpecID)
  nSites = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Coord)
  nSitesRemain = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord)
  Surfnum = nSites - nSitesRemain
  ! check if new_adsorbates lower than number of adsorbates tracked on surface and correct to be removed number
  ! the remaining number of to be removed particles is still kept in tmp array
  IF ((SurfNum - ABS(adsorbates_num)).LT.0) THEN
!    CALL abort(&
!__STAMP__&
!,'Error in AdjustReconstructMapNum: Too few Adsorbates on surface to remove:' &
!,Adsorption%Coordination(PartBoundID,SpecID))
    adsorbates_num = -SurfNum
  END IF
  DO WHILE (dist.GE.adsorbates_num)
    CALL RANDOM_NUMBER(RanNum)
    Surfpos = 1 + INT(Surfnum * RanNum)
    UsedSiteMapPos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain+Surfpos)
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(UsedSiteMapPos) = 0
    ! assign bond order of respective surface atoms in the surfacelattice
    DO iInterAtom = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
      xpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(UsedSiteMapPos,iInterAtom)
      ypos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(UsedSiteMapPos,iInterAtom)
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(SpecID,xpos,ypos) = &
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(SpecID,xpos,ypos) - 1
    END DO
    ! rearrange UsedSiteMap-Surfpos-array
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain+Surfpos) = &
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain+1)
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain+1) = UsedSiteMapPos
    Surfnum = Surfnum - 1
    nSitesRemain = nSitesRemain + 1
    dist = dist - 1
  END DO
  SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord) = nSitesRemain
END IF

END SUBROUTINE SMCR_AdjustMapNum

END MODULE MOD_SurfaceModel_Tools
