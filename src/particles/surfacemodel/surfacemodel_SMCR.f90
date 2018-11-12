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

MODULE MOD_SMCR
!===================================================================================================================================
!> Main Routines of Surface Approximation Monte Carlo
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
PUBLIC :: SMCR_PartAdsorb
PUBLIC :: SMCR_PartDesorb
PUBLIC :: SMCR_Diffusion
!===================================================================================================================================

CONTAINS

SUBROUTINE SMCR_PartAdsorb(subsurfxi,subsurfeta,SurfID,PartID,Norm_Velo,adsorption_case,outSpec,AdsorptionEnthalpie)
!===================================================================================================================================
!> Particle adsorption probability calculation for one impinging particle using a surface reconstruction (SMCR) (surfacemodel = 3)
!===================================================================================================================================
USE MOD_Globals_Vars           ,ONLY: PlanckConst, BoltzmannConst, PI
USE MOD_Particle_Vars          ,ONLY: PartSpecies, nSpecies, Species, WriteMacroSurfaceValues
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_DSMC_Vars              ,ONLY: DSMC, SpecDSMC
USE MOD_SurfaceModel_Vars      ,ONLY: SurfDistInfo, Adsorption
USE MOD_SurfaceModel_Tools     ,ONLY: Calc_Adsorb_Heat, Calc_E_Act, Set_TST_Factors
USE MOD_SurfaceModel_Tools     ,ONLY: CalcAdsorbReactProb, SpaceOccupied, UpdateSurfPos
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound, SampWall, SurfMesh
USE MOD_TimeDisc_Vars          ,ONLY: TEnd, time
#if (PP_TimeDiscMethod==42)
USE MOD_TimeDisc_Vars          ,ONLY: iter, dt
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)               :: subsurfxi,subsurfeta,SurfID,PartID
INTEGER,INTENT(OUT)              :: adsorption_case
INTEGER,INTENT(OUT)              :: outSpec(2)
REAL   ,INTENT(OUT)              :: AdsorptionEnthalpie
REAL   ,INTENT(IN)               :: Norm_Velo
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                          :: iSpec, globSide, PartBoundID
INTEGER                          :: Coord, Coord2, i, AdsorbID, IDRearrange
INTEGER                          :: nSites, nSitesRemain
REAL                             :: WallTemp, RanNum
REAL                             :: Prob_ads
REAL , ALLOCATABLE               :: P_Eley_Rideal(:), Prob_diss(:)
INTEGER                          :: Surfpos, ReactNum
INTEGER                          :: jSpec, kSpec, jCoord, kCoord
REAL                             :: sum_probabilities
INTEGER , ALLOCATABLE            :: NeighbourID(:,:), NeighSpec(:)
INTEGER                          :: SiteSpec, nNeigh_trap, Neighpos_j, Neighpos_k, chosen_Neigh_j, chosen_Neigh_k
INTEGER                          :: n_empty_Neigh(3), n_Neigh(3), adsorbates(nSpecies)
REAL                             :: E_a, c_f
REAL                             :: Heat_A, Heat_B, Heat_AB, D_AB, D_A, D_B
REAL                             :: vel_norm, vel_coll, potential_pot, a_const, mu, surfmass, trapping_prob
LOGICAL                          :: Cell_Occupied
REAL                             :: CharaTemp
INTEGER                          :: DissocReactID, AssocReactID
REAL                             :: a_f, b_f, E_d
REAL                             :: coverage_check
#if (PP_TimeDiscMethod==42)
INTEGER                          :: iSampleReact
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
! special TPD (temperature programmed desorption) surface temperature adjustment part
globSide = Adsorption%SurfSideToGlobSideMap(SurfID)
PartBoundID = PartBound%MapToPartBC(BC(globSide))
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
! calculate number of adsorbates for each species (already on surface)
adsorbates(:) = 0
DO Coord = 1,3
nSites = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%nSites(Coord)
nSitesRemain = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%SitesRemain(Coord)
DO AdsorbID = 1,nSites-nSitesRemain
  Surfpos = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain+AdsorbID)
  iSpec = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%Species(Surfpos)
  adsorbates(iSpec) = adsorbates(iSpec) + 1
END DO
END DO
! initialize variables
Prob_ads = 0.
ALLOCATE( P_Eley_Rideal(1:Adsorption%ReactNum),&
          Prob_diss(1:Adsorption%ReactNum))
P_Eley_Rideal(:) = 0.
Prob_diss(:) = 0.
Cell_Occupied = .FALSE.

! Choose Random surface site with species coordination
CALL RANDOM_NUMBER(RanNum)
iSpec = PartSpecies(PartID)
Coord = Adsorption%Coordination(PartBoundID,iSpec)
AdsorbID = 1 + INT(SurfDistInfo(subsurfxi,subsurfeta,SurfID)%nSites(Coord)*RanNum)
Surfpos = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%UsedSiteMap(AdsorbID)
SiteSpec = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%Species(Surfpos)
!-----------------------------------------------------------------------------------------------------------------------------------
! calculate trapping probability (using hard cube collision with surface atom or adsorbate)
!-----------------------------------------------------------------------------------------------------------------------------------
! if site is empty nearest neighbour site can be occupied and this influences the collision cube mass
IF (SiteSpec.EQ.0) THEN
  ALLOCATE(NeighSpec(1:SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%nNeighbours))
  NeighSpec = 0
  nNeigh_trap = 0
  DO i = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%nNeighbours
    IF (.NOT.SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%IsNearestNeigh(Surfpos,i)) CYCLE
    IF ((Coord.EQ.1) .AND. (SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%NeighSite(Surfpos,i).EQ.1)) CYCLE
    DO Coord2 = 1,3
      IF (Coord2.EQ.SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%NeighSite(Surfpos,i)) THEN
        IF ( (SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord2)%Species( &
              SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%NeighPos(Surfpos,i)) &
              .NE.0) ) THEN
              Cell_Occupied = .TRUE.
              nNeigh_trap = nNeigh_trap + 1
              NeighSpec(nNeigh_trap) = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord2)%Species( &
              SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%NeighPos(Surfpos,i))
        END IF
      END IF
    END DO
  END DO
END IF
potential_pot = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfID,iSpec,Surfpos,.TRUE.)*BoltzmannConst
vel_norm = - (( 2*(potential_pot)/Species(iSpec)%MassIC)**(0.5) + Norm_Velo)
a_const = (Species(iSpec)%MassIC/(2*BoltzmannConst*WallTemp))**(0.5)
IF (SiteSpec.EQ.0) THEN
  IF (nNeigh_trap.EQ.0) THEN
    surfmass = PartBound%SolidMassIC(PartBoundID) !Adsorption%SurfMassIC(SurfID)
  ELSE
    surfmass = 0.
    DO i = 1,nNeigh_trap
      surfmass = surfmass + Species(NeighSpec(i))%MassIC
    END DO
    surfmass = surfmass / nNeigh_trap
  END IF
  SDEALLOCATE(NeighSpec)
  mu = Species(iSpec)%MassIC / surfmass
ELSE
  mu = Species(iSpec)%MassIC / (Species(SiteSpec)%MassIC)
END IF
vel_coll = 0.5 * ( (1+mu)*(2*(potential_pot/Species(iSpec)%MassIC))**(0.5) + (1-mu)*vel_norm )
trapping_prob = abs(0.5 + 0.5*ERF(a_const*vel_coll) + ( EXP(-(a_const**2)*(vel_coll**2)) / (2*a_const*(vel_norm*PI**0.5)) ))
IF (trapping_prob.GT.1.) trapping_prob = 1.
#if (PP_TimeDiscMethod==42)
Adsorption%AdsorpInfo(iSpec)%Accomodation = Adsorption%AdsorpInfo(iSpec)%Accomodation + trapping_prob
#endif
IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
  SampWall(SurfID)%Accomodation(iSpec,subsurfxi,subsurfeta) = SampWall(SurfID)%Accomodation(iSpec,subsurfxi,subsurfeta) &
                                                                + trapping_prob
END IF
! adaptive accomodation
!IF (Adaptive_ACC_FLAG) THEN
  CALL RANDOM_NUMBER(RanNum)
  IF(RanNum.GE.PartBound%MomentumACC(PartBoundID)) THEN
    outSpec(1) = iSpec
    outSpec(2) = 0
    AdsorptionEnthalpie = 0.
    adsorption_case = -1
    RETURN
  END IF
!END IF

!! if no trapping return and perform elastic reflection
!CALL RANDOM_NUMBER(RanNum)
!IF (RanNum.GT.trapping_prob) THEN
!  outSpec(1) = iSpec
!  outSpec(2) = 0
!  AdsorptionEnthalpie = 0.
!  adsorption_case = -1
!  RETURN
!END IF
c_f = BoltzmannConst/PlanckConst &
    !* REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfID)%nSites(3))/SurfMesh%SurfaceArea(subsurfxi,subsurfeta,SurfID) &
    * REAL(Adsorption%DensSurfAtoms(SurfID)*Adsorption%AreaIncrease(SurfID)) &
    / ( (BoltzmannConst / (2*Pi*Species(iSpec)%MassIC))**0.5 )
!-----------------------------------------------------------------------------------------------------------------------------------
! calculate probability for molecular adsorption
!-----------------------------------------------------------------------------------------------------------------------------------
IF ( (SiteSpec.EQ.0) .AND. (.NOT.Cell_Occupied) &
    .AND. (INT(SurfDistInfo(subsurfxi,subsurfeta,SurfID)%adsorbnum_tmp(iSpec)).LE.0) ) THEN
  ! calculation of molecular adsorption probability with TCE
  E_a = 0.
  E_d = 0.1 * Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfID,iSpec,Surfpos,.TRUE.) * Boltzmannconst
  CALL Set_TST_Factors(1,a_f,b_f,PartID,0)!,PartBoundID)
  Prob_ads = CalcAdsorbReactProb(1,PartID,Norm_velo,E_a,E_d,a_f,b_f,c_f)
#if (PP_TimeDiscMethod==42)
  iSampleReact = 1
  IF ((.NOT.DSMC%ReservoirRateStatistic).AND.(Prob_Ads.GT.0.)) THEN
    !IF (Prob_Ads.GT.1) THEN
    !  Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(iSampleReact) = &
    !      Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(iSampleReact) + 1.
    !  Adsorption%AdsorpInfo(iSpec)%MeanProbAds = Adsorption%AdsorpInfo(iSpec)%MeanProbAds + 1.
    !ELSE
      Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(iSampleReact) = &
          Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(iSampleReact) + Prob_ads
      Adsorption%AdsorpInfo(iSpec)%MeanProbAds = Adsorption%AdsorpInfo(iSpec)%MeanProbAds + Prob_ads
    !END IF
  END IF
  Adsorption%AdsorpReactInfo(iSpec)%AdsReactCount(iSampleReact) = &
      Adsorption%AdsorpReactInfo(iSpec)%AdsReactCount(iSampleReact) + 1
#endif
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! sort Neighbours to coordinations for search of two random neighbour positions from impact position for dissociative adsorption
!-----------------------------------------------------------------------------------------------------------------------------------
n_Neigh(:) = 0
n_empty_Neigh(:) = 0
ALLOCATE(NeighbourID(1:3,1:SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%nNeighbours))

DO i = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%nNeighbours
DO Coord2 = 1,3
  IF (Coord2.EQ.SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%NeighSite(Surfpos,i) &
      .AND.SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%IsNearestNeigh(Surfpos,i)) THEN
    IF ( (SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord2)%Species( &
          SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%NeighPos(Surfpos,i)) &
          .EQ.0) ) THEN
      n_empty_Neigh(Coord2) = n_empty_Neigh(Coord2) + 1
    n_Neigh(Coord2) = n_Neigh(Coord2) + 1
    NeighbourID(Coord2,n_Neigh(Coord2)) = i
    END IF
  END IF
END DO
END DO
!---------------------------------------------------------------------------------------------------------------------------------
! calculate probability for dissociative adsorption
!---------------------------------------------------------------------------------------------------------------------------------
IF ((SpecDSMC(iSpec)%InterID.EQ.2)) THEN ! particle has to be at least di-atomic
DO ReactNum = 1,(Adsorption%DissNum)
  jSpec = Adsorption%DissocReact(1,ReactNum,iSpec)
  kSpec = Adsorption%DissocReact(2,ReactNum,iSpec)
  IF ((jSpec.NE.0) .AND. (kSpec.NE.0)) THEN !if 2 resulting species, dissociation possible
    ! in one of the last iterations the surface-coverage of one resulting species was saturated
    IF (INT(SurfDistInfo(subsurfxi,subsurfeta,SurfID)%adsorbnum_tmp(jSpec)).GT.0 .OR. &
        INT(SurfDistInfo(subsurfxi,subsurfeta,SurfID)%adsorbnum_tmp(kSpec)).GT.0 ) THEN
      CYCLE
    END IF
    jCoord = Adsorption%Coordination(PartBoundID,jSpec)
    kCoord = Adsorption%Coordination(PartBoundID,kSpec)
    Neighpos_j = 0
    Neighpos_k = 0
    !At least one of resulting species has same coordination as surfaceposition of impact
    IF ( (jCoord.EQ.Coord) .OR. (kCoord.EQ.Coord) ) THEN
      IF (jCoord.EQ.Coord) THEN
        Neighpos_j = Surfpos
        IF (n_empty_Neigh(kCoord).GT.0) THEN
          ! assign availiable neighbour position for k
          CALL RANDOM_NUMBER(RanNum)
          chosen_Neigh_k = 1 + INT(REAL(n_Neigh(kCoord))*RanNum)
          Neighpos_k = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%NeighPos(Surfpos,NeighbourID(kCoord,chosen_Neigh_k))
          IDRearrange = NeighbourID(kCoord,chosen_Neigh_k)
          NeighbourID(kCoord,chosen_Neigh_k) = NeighbourID(kCoord,n_Neigh(kCoord))
          NeighbourID(kCoord,n_Neigh(kCoord)) = IDRearrange
        END IF
      ELSE
        Neighpos_k = Surfpos
        IF (n_empty_Neigh(jCoord).GT.0) THEN
          ! assign availiable neighbour position for j
          CALL RANDOM_NUMBER(RanNum)
          chosen_Neigh_j = 1 + INT(REAL(n_Neigh(jCoord))*RanNum)
          Neighpos_j = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%NeighPos(Surfpos,NeighbourID(jCoord,chosen_Neigh_j))
          IDRearrange = NeighbourID(jCoord,chosen_Neigh_j)
          NeighbourID(jCoord,chosen_Neigh_j) = NeighbourID(jCoord,n_Neigh(jCoord))
          NeighbourID(jCoord,n_Neigh(jCoord)) = IDRearrange
        END IF
      END IF
    !NEITHER of resulting species have same coordination as surfaceposition of impact
    ELSE
      IF ( (jCoord.EQ.kCoord) .AND. (n_empty_Neigh(jCoord).GT.1) ) THEN !resulting species are equal
        ! assign availiable neighbour positions
        CALL RANDOM_NUMBER(RanNum)
        chosen_Neigh_j = 1 + INT(REAL(n_Neigh(jCoord))*RanNum)
        Neighpos_j = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%NeighPos(Surfpos,NeighbourID(jCoord,chosen_Neigh_j))
        IDRearrange = NeighbourID(jCoord,chosen_Neigh_j)
        NeighbourID(jCoord,chosen_Neigh_j) = NeighbourID(jCoord,n_Neigh(jCoord))
        NeighbourID(jCoord,n_Neigh(jCoord)) = IDRearrange
        CALL RANDOM_NUMBER(RanNum)
        chosen_Neigh_k = 1 + INT(REAL(n_Neigh(kCoord)-1)*RanNum)
        Neighpos_k = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%NeighPos(Surfpos,NeighbourID(kCoord,chosen_Neigh_k))
        IDRearrange = NeighbourID(kCoord,chosen_Neigh_k)
        NeighbourID(kCoord,chosen_Neigh_k) = NeighbourID(kCoord,n_Neigh(kCoord)-1)
        NeighbourID(kCoord,n_Neigh(kCoord)-1) = IDRearrange
      ELSE IF ( (jCoord.NE.kCoord) .AND. (n_empty_Neigh(jCoord).GT.0) .AND. (n_empty_Neigh(kCoord).GT.0) ) THEN
        ! assign availiable neighbour positions
        CALL RANDOM_NUMBER(RanNum)
        chosen_Neigh_j = 1 + INT(REAL(n_Neigh(jCoord))*RanNum)
        Neighpos_j = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%NeighPos(Surfpos,NeighbourID(jCoord,chosen_Neigh_j))
        IDRearrange = NeighbourID(jCoord,chosen_Neigh_j)
        NeighbourID(jCoord,chosen_Neigh_j) = NeighbourID(jCoord,n_Neigh(jCoord))
        NeighbourID(jCoord,n_Neigh(jCoord)) = IDRearrange
        CALL RANDOM_NUMBER(RanNum)
        chosen_Neigh_k = 1 + INT(REAL(n_Neigh(kCoord))*RanNum)
        Neighpos_k = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%NeighPos(Surfpos,NeighbourID(kCoord,chosen_Neigh_k))
        IDRearrange = NeighbourID(kCoord,chosen_Neigh_j)
        NeighbourID(kCoord,chosen_Neigh_k) = NeighbourID(kCoord,n_Neigh(kCoord))
        NeighbourID(kCoord,n_Neigh(kCoord)) = IDRearrange
      END IF
    END IF
    IF ( (Neighpos_j.GT.0) .AND. (Neighpos_k.GT.0) ) THEN !both neighbour positions assigned
    !both assigned neighbour positions empty
    IF ( (SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(jCoord)%Species(Neighpos_j).EQ.0) .AND. &
        (SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(kCoord)%Species(Neighpos_k).EQ.0) ) THEN
      ! check for occupation of neirest Neighbours of both positions and Cycle if space of at least one position already occupied
      IF ( SpaceOccupied(SurfID,subsurfxi,subsurfeta,jCoord,Neighpos_j) &
          .OR. SpaceOccupied(SurfID,subsurfxi,subsurfeta,kCoord,Neighpos_k) ) CYCLE
      ! assign bond order of respective surface atoms in the surfacelattice for molecule
      CALL UpdateSurfPos(SurfID,subsurfxi,subsurfeta,Coord,Surfpos,iSpec,.FALSE.)
      ! calculation of activation energy
      Heat_A = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfID,jSpec,Neighpos_j,.TRUE.)
      Heat_B = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfID,kSpec,Neighpos_k,.TRUE.)
      ! remove bond order of respective surface atoms in the surfacelattice for molecule again
      CALL UpdateSurfPos(SurfID,subsurfxi,subsurfeta,Coord,Surfpos,iSpec,.TRUE.)
      Heat_AB = 0. ! direct dissociative adsorption
      D_AB = Adsorption%EDissBond(ReactNum,iSpec)
      D_A = 0.
      D_B = 0.
      E_a = Calc_E_Act(Heat_A,Heat_B,Heat_AB,0.,D_A,D_B,D_AB,0.) * BoltzmannConst
      E_d = 0.0  !0.1 * Calc_E_Act(Heat_AB,0.,Heat_A,Heat_B,D_AB,0.,D_A,D_B) * BoltzmannConst
      ! calculation of dissociative adsorption probability with TCE
      CALL Set_TST_Factors(2,a_f,b_f,PartID,ReactNum)!,PartBoundID)
      Prob_diss(ReactNum) = CalcAdsorbReactProb(2,PartID,Norm_Velo,E_a,E_d,a_f,b_f,c_f)
#if (PP_TimeDiscMethod==42)
      iSampleReact = 1 + ReactNum
      IF ((.NOT.DSMC%ReservoirRateStatistic).AND.(Prob_diss(ReactNum).GT.0.)) THEN
        !IF (Prob_diss(ReactNum).GT.1) THEN
        !  Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(iSampleReact) = &
        !      Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(iSampleReact) + 1.
        !  Adsorption%AdsorpInfo(iSpec)%MeanProbAds = Adsorption%AdsorpInfo(iSpec)%MeanProbAds + 1.
        !ELSE
          Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(iSampleReact) = &
              Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(iSampleReact) + Prob_diss(ReactNum)
          Adsorption%AdsorpInfo(iSpec)%MeanProbAds = Adsorption%AdsorpInfo(iSpec)%MeanProbAds + Prob_diss(ReactNum)
        !END IF
      END IF
      Adsorption%AdsorpReactInfo(iSpec)%MeanAdsActE(iSampleReact-1) = &
          Adsorption%AdsorpReactInfo(iSpec)%MeanAdsActE(iSampleReact-1) + E_a/BoltzmannConst
      Adsorption%AdsorpReactInfo(iSpec)%AdsReactCount(iSampleReact) = &
          Adsorption%AdsorpReactInfo(iSpec)%AdsReactCount(iSampleReact) + 1
#endif
    END IF
    END IF !both neighbour positions assigned
  END IF
END DO
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! calculate probability for Eley-Rideal reaction
!-----------------------------------------------------------------------------------------------------------------------------------
DO ReactNum = 1,(Adsorption%RecombNum)
  ! reaction partner
  jSpec = Adsorption%AssocReact(1,ReactNum,iSpec)
  Neighpos_j = 0
  IF (jSpec.EQ.0) CYCLE
  IF (INT(SurfDistInfo(subsurfxi,subsurfeta,SurfID)%adsorbnum_tmp(jSpec)).LT.0 ) CYCLE
  coverage_check = Adsorption%Coverage(subsurfxi,subsurfeta,SurfID,jSpec) &
                 + Adsorption%SumAdsorbPart(subsurfxi,subsurfeta,SurfID,jSpec) &
                 * Species(jSpec)%MacroParticleFactor &
                 / REAL(INT(Adsorption%DensSurfAtoms(SurfID) * SurfMesh%SurfaceArea(subsurfxi,subsurfeta,SurfID),8))
  IF (coverage_check.LE.0.) CYCLE
  ! reaction results
  kSpec = Adsorption%AssocReact(2,ReactNum,iSpec)
  jCoord = Adsorption%Coordination(PartBoundID,jSpec)
  ! Choose Random surface site with reactionpartner coordination
  CALL RANDOM_NUMBER(RanNum)
  AdsorbID = 1 + INT(SurfDistInfo(subsurfxi,subsurfeta,SurfID)%nSites(jCoord)*RanNum)
  Neighpos_j = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(jCoord)%UsedSiteMap(AdsorbID)
  !IF ( n_Neigh(jCoord)-n_empty_Neigh(jCoord).GT.0 ) THEN
  !  CALL RANDOM_NUMBER(RanNum)
  !  chosen_Neigh_j = 1 + INT(n_Neigh(jCoord)*RanNum)
  !  Neighpos_j = &
  !          SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(jCoord)%NeighPos(Surfpos,NeighbourID(jCoord,chosen_Neigh_j))
  !END IF
  IF ( (Neighpos_j.GT.0) ) THEN
  IF ( (SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(jCoord)%Species(Neighpos_j).EQ.jSpec) ) THEN
    ! calculation of activation energy
    Heat_A = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfID,jSpec,Neighpos_j,.TRUE.)
    Heat_B = 0. !Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfID,kSpec,Neighpos_k,.TRUE.)
    Heat_AB = 0. ! direct associative desorption
    D_AB = Adsorption%EDissBond(ReactNum+Adsorption%DissNum,iSpec)
    D_A = 0.
    D_B = 0.
    E_a = Calc_E_Act(Heat_AB,0.,Heat_A,Heat_B,D_AB,0.,D_A,D_B) * BoltzmannConst
    E_d = 0.0 !0.1 * Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfID,jSpec,Neighpos_j,.TRUE.)
    ! estimate characteristic vibrational temperature of surface-particle bond
    IF (Coord.EQ.1) THEN
      CharaTemp = Heat_A / 200. / (2 - 1./Adsorption%CrystalIndx(SurfID))
    ELSE
      CharaTemp = Heat_A / 200. / (2 - 1./SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%nInterAtom)
    END IF
    ! calculation of ER-reaction probability with TCE
    CALL Set_TST_Factors(3,a_f,b_f,PartID,ReactNum)!,PartBoundID)
    P_Eley_Rideal(ReactNum+Adsorption%DissNum) = CalcAdsorbReactProb(3,PartID,Norm_Velo,E_a,E_d,a_f,b_f,c_f &
                             ,PartnerSpecies=jSpec,CharaTemp=CharaTemp,WallTemp=WallTemp)
#if (PP_TimeDiscMethod==42)
    iSampleReact = 1 + Adsorption%DissNum + ReactNum
    IF ((.NOT.DSMC%ReservoirRateStatistic).AND.(P_Eley_Rideal(ReactNum+Adsorption%DissNum).GT.0.)) THEN
      !IF (P_Eley_Rideal(ReactNum).GT.1) THEN
      !  Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(iSampleReact) = &
      !      Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(iSampleReact) + 1.
      !ELSE
        Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(iSampleReact) = &
            Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(iSampleReact) + P_Eley_Rideal(ReactNum+Adsorption%DissNum)
      !END IF
    END IF
    Adsorption%AdsorpReactInfo(iSpec)%MeanAdsActE(iSampleReact-1) = &
        Adsorption%AdsorpReactInfo(iSpec)%MeanAdsActE(iSampleReact-1) + E_a/BoltzmannConst
    Adsorption%AdsorpReactInfo(iSpec)%AdsReactCount(iSampleReact) = &
        Adsorption%AdsorpReactInfo(iSpec)%AdsReactCount(iSampleReact) + 1
#endif
  END IF
  END IF
END DO
!-----------------------------------------------------------------------------------------------------------------------------------

sum_probabilities = Prob_ads
DO ReactNum = 1,Adsorption%ReactNum
  sum_probabilities = sum_probabilities + Prob_diss(ReactNum) + P_Eley_Rideal(ReactNum)
END DO
! initialize adsorption case
adsorption_case = 0
AdsorptionEnthalpie = 0.
!-----------------------------------------------------------------------------------------------------------------------------------
! choose which adsorption type takes place
!-----------------------------------------------------------------------------------------------------------------------------------
IF (DSMC%ReservoirSurfaceRate) sum_probabilities = 0. !only probabilities are calculated without adsorption taking place

CALL RANDOM_NUMBER(RanNum)
IF (sum_probabilities .GT. RanNum) THEN
  ! chose surface reaction case (0=inelastic scattering, 1=adsorption, 2=reaction (dissociation), 3=reaction (Eley-Rideal))
  DO ReactNum = 1,(Adsorption%ReactNum)
    CALL RANDOM_NUMBER(RanNum)
    IF ((P_Eley_Rideal(ReactNum)/sum_probabilities).GT.RanNum) THEN
      ! if ER-reaction set output parameters
      adsorption_case = 3
      AssocReactID = ReactNum - Adsorption%DissNum
      outSpec(1) = Adsorption%AssocReact(1,AssocReactID,iSpec)
      outSpec(2) = Adsorption%AssocReact(2,AssocReactID,iSpec)
      ! calculate adsorption Enthalpie
      Heat_A = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfID,outSpec(1),-1,.TRUE.)
      Heat_B = 0.
      Heat_AB = 0.
      D_AB = Adsorption%EDissBond(ReactNum,iSpec)
      D_A = 0.
      D_B = 0.
      AdsorptionEnthalpie = -(( Heat_AB -Heat_A -Heat_B ) + ( D_AB -D_A -D_B )) * BoltzmannConst
#if (PP_TimeDiscMethod==42)
      iSampleReact = 1 + ReactNum
      IF (DSMC%ReservoirRateStatistic) THEN
        Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(iSampleReact) = &
            Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(iSampleReact) + 1.
      END IF
#endif
      EXIT
    END IF
    sum_probabilities = sum_probabilities - P_Eley_Rideal(ReactNum)
    CALL RANDOM_NUMBER(RanNum)
    IF ((Prob_diss(ReactNum)/sum_probabilities).GT.RanNum) THEN
      ! if dissocciative adsorption set output parameters
      adsorption_case = 2
      DissocReactID = ReactNum
      outSpec(1) = Adsorption%DissocReact(1,DissocReactID,iSpec)
      outSpec(2) = Adsorption%DissocReact(2,DissocReactID,iSpec)
      ! calculate adsorption Enthalpie
      Heat_A = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfID,outSpec(1),-1,.TRUE.)
      Heat_B = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfID,outSpec(2),-1,.TRUE.)
      Heat_AB = 0.
      D_AB = Adsorption%EDissBond(ReactNum,iSpec)
      D_A = 0.
      D_B = 0.
      AdsorptionEnthalpie = (( Heat_AB -Heat_A -Heat_B ) + ( D_AB -D_A -D_B )) * BoltzmannConst
#if (PP_TimeDiscMethod==42)
      iSampleReact = 1 + ReactNum
      IF (DSMC%ReservoirRateStatistic) THEN
        Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(iSampleReact) = &
            Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(iSampleReact) + 1.
        Adsorption%AdsorpInfo(iSpec)%MeanProbAds = Adsorption%AdsorpInfo(iSpec)%MeanProbAds + 1.
      END IF
#endif
      EXIT
    END IF
    sum_probabilities = sum_probabilities - Prob_diss(ReactNum)
  END DO
  IF (adsorption_case.EQ.0) THEN
    CALL RANDOM_NUMBER(RanNum)
    IF ((Prob_ads/sum_probabilities).GT.RanNum) THEN
      ! if molecular adsorption set output parameters
      adsorption_case = 1
      outSpec(1) = iSpec
      outSpec(2) = 0
      ! calculate adsorption Enthalpie
      AdsorptionEnthalpie = - Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfID,outSpec(1),-1,.TRUE.) * BoltzmannConst
#if (PP_TimeDiscMethod==42)
      iSampleReact = 1 + ReactNum
      IF (DSMC%ReservoirRateStatistic) THEN
        Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(iSampleReact) = &
            Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(iSampleReact) + 1.
        Adsorption%AdsorpInfo(iSpec)%MeanProbAds = Adsorption%AdsorpInfo(iSpec)%MeanProbAds + 1.
      END IF
#endif
    END IF
  END IF
END IF

DEALLOCATE(P_Eley_Rideal,Prob_diss,NeighbourID)

END SUBROUTINE SMCR_PartAdsorb

SUBROUTINE SMCR_PartDesorb()
!===================================================================================================================================
!> Calculation of number of desorbing particles using surface reconstruction (surfacemodel = 3)
!> Performing Surface Monte Carlo step (MCS)
!> 1. diffusion into equilibrium (Quasi Chemical Approximation - QCA) is performed for particles on surface
!> 2. Loop over all particles on surfaces and calulate desorption/reaction probabities
!>   - choose respective process
!>   - update surface distribution
!===================================================================================================================================
USE MOD_Globals_Vars           ,ONLY: PlanckConst, BoltzmannConst
USE MOD_Particle_Vars          ,ONLY: nSpecies, Species, WriteMacroSurfaceValues
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_DSMC_Vars              ,ONLY: DSMC, SpecDSMC, PolyatomMolDSMC
USE MOD_SurfaceModel_Vars      ,ONLY: SurfDistInfo, Adsorption
USE MOD_SurfaceModel_Tools     ,ONLY: Calc_Adsorb_Heat, Calc_E_Act
USE MOD_SurfaceModel_Tools     ,ONLY: CalcAdsorbReactProb, SpaceOccupied, UpdateSurfPos
USE MOD_SurfaceModel_PartFunc  ,ONLY: PartitionFuncAct, PartitionFuncAct_dissoc
USE MOD_SurfaceModel_PartFunc  ,ONLY: PartitionFuncSurf, PartitionFuncAct_recomb, PartitionFuncAct_exch
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, PartBound, SampWall
USE MOD_TimeDisc_Vars          ,ONLY: dt
USE MOD_TimeDisc_Vars          ,ONLY: TEnd, time
#if (PP_TimeDiscMethod==42)
USE MOD_TimeDisc_Vars          ,ONLY: iter
#endif
#if USE_LOADBALANCE
USE MOD_Particle_Mesh_Vars     ,ONLY: PartSideToElem
USE MOD_LoadBalance_Vars       ,ONLY: nSurfacePartsPerElem, PerformLBSample
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
IMPLICIT NONE
!===================================================================================================================================
! Local variable declaration
INTEGER                           :: iSurf, iSpec, globSide, jSubSurf, iSubSurf, PartBoundID
INTEGER                           :: Coord, Coord3, i, j, AdsorbID, numSites
REAL                              :: WallTemp, RanNum
REAL                              :: Heat_A, Heat_B, nu_des, rate, P_actual_des
INTEGER                           :: Indx, Indy, Surfpos
INTEGER                           :: trace, traceNum, ReactNum_run
INTEGER , ALLOCATABLE             :: desorbnum(:), reactdesorbnum(:)
INTEGER , ALLOCATABLE             :: adsorbnum(:)
INTEGER , ALLOCATABLE             :: nSites(:), nSitesRemain(:), remainNum(:), adsorbates(:)
LOGICAL                           :: Cell_Occupied
!---------- reaction variables
INTEGER                           :: n_Neigh(3), n_empty_Neigh(3), jSpec, iReact, PartnerID!,react_Neigh
INTEGER                           :: surf_react_case, NeighID
REAL                              :: E_d, nu_react, E_diff
REAL                              :: Heat_AB, D_AB, sum_probabilities
REAL                              :: VarPartitionFuncWall1, VarPartitionFuncWall2, VarPartitionFuncAct
REAL                              :: CharaTemp
INTEGER , ALLOCATABLE             :: Pos_ReactP(:), Coord_ReactP(:), Pos_Product(:,:), Coord_Product(:,:)
INTEGER , ALLOCATABLE             :: React_NeighbourID(:,:), NeighbourID(:,:)
REAL , ALLOCATABLE                :: P_actual_react(:),P_react_forward(:),P_react_back(:)
REAL                              :: AdsorptionEnthalpie
!variables used for sampling of vibrational energies of reacting/desorbing particles
INTEGER                           :: SampleParts, iPart
INTEGER                           :: iPolyatMole, iDOF
INTEGER                           :: VibQuantWall
INTEGER, ALLOCATABLE              :: VibQuantWallPoly(:)
REAL, ALLOCATABLE                 :: RanNumPoly(:)
REAL                              :: EvibOld, EvibWall, EVibNew
! variables used for exchange reactions
REAL                              :: D_A, D_ReactP, D_Product1, D_Product2, Heat_ReactP, Heat_Product1, Heat_Product2
INTEGER                           :: Prod_Spec1, Prod_Spec2
LOGICAL                           :: ReactDirForward, exch_react_possible
INTEGER                           :: DissocNum, ExchNum, jCoord, kCoord, Neighpos_j, Neighpos_k, chosen_Neigh_j, chosen_Neigh_k
INTEGER                           :: IDRearrange
#if (PP_TimeDiscMethod==42)
! reservoir sample variables
INTEGER                           :: iSampleReact
REAL                              :: loc_SurfActE(1:Adsorption%ReactNum+Adsorption%NumOfExchReact+1)
#endif
INTEGER                           :: NumDesorbLH(1:nSpecies,1:Adsorption%RecombNum)
#if USE_LOADBALANCE
INTEGER                           :: ElemID
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
IF (.NOT.SurfMesh%SurfOnProc) RETURN

! diffusion into equilibrium distribution
CALL SMCR_Diffusion()

! allocate number of reactions and desorptions in current MCS
ALLOCATE (desorbnum(1:nSpecies),&
          reactdesorbnum(1:nSpecies))
! allocate arrays for surface positions
ALLOCATE( adsorbnum(1:4),&
          nSites(1:4),&
          nSitesRemain(1:4),&
          remainNum(1:4),&
          adsorbates(1:nSpecies))
! allocate reaction probability arrays
ALLOCATE( P_react_forward(1:Adsorption%nDisPropReactions),&
          P_react_back(1:Adsorption%nDisPropReactions),&
          P_actual_react(1:Adsorption%ReactNum+1+Adsorption%nDisPropReactions),&
          Coord_ReactP(1:Adsorption%ReactNum+Adsorption%nDisPropReactions),&
          Coord_Product(1:2,1:Adsorption%ReactNum+Adsorption%nDisPropReactions),&
          Pos_ReactP(1:Adsorption%ReactNum+1+Adsorption%nDisPropReactions),&
          Pos_Product(1:2,1:Adsorption%ReactNum+Adsorption%nDisPropReactions))

! loop over all surfaces and decide if catalytic boundary
DO iSurf = 1,SurfMesh%nSides
  globSide = Adsorption%SurfSideToGlobSideMap(iSurf)
  PartBoundID = PartBound%MapToPartBC(BC(globSide))
  IF (.NOT.PartBound%SolidCatalytic(PartboundID)) CYCLE
#if USE_LOADBALANCE
  IF(PerformLBSample) ElemID = PartSideToElem(S2E_ELEM_ID,globSide)
#endif /*USE_LOADBALANCE*/
! special TPD (temperature programmed desorption) surface temperature adjustment part
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

! loop over all subsurfaces
DO jSubSurf = 1,nSurfSample
DO iSubSurf = 1,nSurfSample
  desorbnum(:) = 0
  reactdesorbnum(:) = 0
  nSites(:) = 0
  nSitesRemain(:) = 0
  remainNum(:) = 0
  P_react_back(:) = 0.

  NumDesorbLH(:,:) = 0

  DO Coord=1,3
    nSites(Coord) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%nSites(Coord)
    nSitesRemain(Coord) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%SitesRemain(Coord)
    nSites(4) = nSites(4) + nSites(Coord)
    nSitesRemain(4) = nSitesRemain(4) + nSitesRemain(Coord)
  END DO
  adsorbnum(:) = nSites(:) - nSitesRemain(:)
  traceNum = adsorbnum(4) ! number of all surface particles to be considered
  trace = 1 ! trace number for tracking of already considered surface particles (new adsorbates resulting from reaction not needed)
  adsorbates(:) = 0

  ! structure of UsedSiteMap array for one coordination
  !           [<---------------nSites-------------->]
  ! Name    :  nSitesRemain                remainNum
  ! AdsorbID:  1 2 3          4  5 6  7    8 9 10 11
  ! SurfPos :  1 7 8          11 9 10 3    4 5 2  6

  ! calculate number of adsorbates for each species (for analyze)
  DO Coord = 1,3
  DO AdsorbID = 1,nSites(Coord)-nSitesRemain(Coord)
    Surfpos = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+AdsorbID)
    iSpec = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%Species(Surfpos)
    adsorbates(iSpec) = adsorbates(iSpec) + 1
#if USE_LOADBALANCE
    IF(PerformLBSample) nSurfacePartsPerElem(ElemID) = nSurfacePartsPerElem(ElemID) + 1
#endif /*USE_LOADBALANCE*/
  END DO
  END DO

  IF (adsorbnum(4) .GT. 0) THEN
  DO WHILE (trace.LE.traceNum)
    ! choose random adsorption-ID
    CALL RANDOM_NUMBER(RanNum)
    AdsorbID = nSitesRemain(4) + 1 + INT((adsorbnum(4)-remainNum(4))*RanNum)
    ! choose right coordination and adjust adsorbate ID
    IF ( AdsorbID.GT.(adsorbnum(1)+adsorbnum(2)+nSitesRemain(4)-remainNum(1)-remainNum(2)) ) THEN
      AdsorbID = AdsorbID - (adsorbnum(1) + adsorbnum(2) + nSitesRemain(1) + nSitesRemain(2) - remainNum(1) - remainNum(2))
!       AdsorbID = AdsorbID - (nSitesRemain(4) + adsorbnum(1) + adsorbnum(2) - remainNum(1) - remainNum(2)) + nSitesRemain(3)
      Coord = 3
    ELSE IF ( (AdsorbID.GT.(adsorbnum(1)+nSitesRemain(4)-remainNum(1))) &
             .AND. (AdsorbID.LE.(adsorbnum(1)+adsorbnum(2)+nSitesRemain(4))-remainNum(1)-remainNum(2)) ) THEN
      AdsorbID = AdsorbID - (adsorbnum(1) + nSitesRemain(1) + nSitesRemain(3) - remainNum(1))
      Coord = 2
    ELSE IF ((AdsorbID.GT.nSitesRemain(4)) .AND. (AdsorbID.LE.(adsorbnum(1)+nSitesRemain(4)-remainNum(1))) ) THEN
      AdsorbID = AdsorbID - (nSitesRemain(2) + nSitesRemain(3))
      Coord = 1
    ELSE
      CYCLE
    END IF

#if (PP_TimeDiscMethod==42)
    loc_SurfActE = 0.
#endif

    Surfpos = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(AdsorbID)
    iSpec = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%Species(Surfpos)

    ! move considered particle to end of array, at beginning of the already considered and on surface remaining particles
    SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = &
        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(nSites(Coord)-remainNum(Coord))
    SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(nSites(Coord)-remainNum(Coord)) = Surfpos
    AdsorbID = nSites(Coord)-remainNum(Coord)

    ! calculate heat of adsorption for actual site
    Heat_A = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,iSpec,Surfpos,.FALSE.)
    E_diff = 0.

    P_actual_react(:) = 0.
    Coord_ReactP(:) = -1
    Coord_Product(:,:) = -1
    Pos_ReactP(:) = 0
    Pos_Product(:,:) = 0
    ReactNum_run = 0
    ! Choose Random surface positions appropriate for each possible Reaction and calculate reaction probability
    !-------------------------------------------------------------------------------------------------------------------------------
    ! calculate probability for associative reactions
    !-------------------------------------------------------------------------------------------------------------------------------
    ! product of associative reactions desorb after reaction as in most cases exothermic (LH-reaction)
    DO iReact = 1,Adsorption%RecombNum
      IF (Adsorption%AssocReact(1,iReact,iSpec).LT.1) CYCLE ! no partner for this associative reaction
      Coord_ReactP(iReact) = Adsorption%Coordination(PartBoundID,Adsorption%AssocReact(1,iReact,iSpec))
      CALL RANDOM_NUMBER(RanNum)
      IF (Coord.EQ.Coord_ReactP(iReact)) THEN
        NeighID = 1 + INT((nSites(Coord_ReactP(iReact))-remainNum(Coord_ReactP(iReact))-1)*RanNum)
      ELSE
        NeighID = 1 + INT((nSites(Coord_ReactP(iReact))-remainNum(Coord_ReactP(iReact)))*RanNum)
      END IF
      Pos_ReactP(iReact) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_ReactP(iReact))%UsedSiteMap(NeighID)
      jSpec = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_ReactP(iReact))%Species(Pos_ReactP(iReact))

      IF ( jSpec.EQ.Adsorption%AssocReact(1,iReact,iSpec) .AND. (jSpec.GT.0)) THEN
        Prod_Spec1 = Adsorption%AssocReact(2,iReact,iSpec)
!         Coord_Product(1,iReact) = Adsorption%Coordination(PartBoundID,Prod_Spec1)
!         CALL RANDOM_NUMBER(RanNum)
!         NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,iReact))))*RanNum)
!         Pos_Product(1,iReact) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
!                                           Coord_Product(1,iReact))%UsedSiteMap(NeighID)
!         Heat_AB = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Prod_Spec1,Pos_Product(1,iReact),.TRUE.)
        ! calculate heats of adsorption
        Heat_AB = 0. ! immediate desorption
        Heat_B = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,jSpec,Pos_ReactP(iReact),.FALSE.)
        D_AB = Adsorption%EDissBond((Adsorption%DissNum+iReact),iSpec)
        E_d = Calc_E_Act(Heat_AB,0.,Heat_A,Heat_B,D_AB,0.,0.,0.)
        ! check diffusion barrier
        IF ((Coord.NE.3).OR.(Coord_ReactP(iReact).NE.3)) THEN
          IF ((Coord.NE.3).AND.(Coord_ReactP(iReact).NE.3)) THEN
            E_diff = (2* (0.1*Heat_A*0.1*Heat_B))/(0.1*Heat_A + 0.1*Heat_B)
          ELSE IF ((Coord.NE.3).AND.(Coord_ReactP(iReact).EQ.3)) THEN
            E_diff = 0.1*Heat_A
          ELSE IF ((Coord.EQ.3).AND.(Coord_ReactP(iReact).NE.3)) THEN
            E_diff = 0.1*Heat_B
          ELSE
            E_diff = 0.
          END IF
        END IF
        ! if diffusion barrier is greater then reaction activation, then it is limiting --> activation energy
        IF (E_Diff.GT.E_d) E_d = E_diff
        ! estimate vibrational temperatures of surface-particle bond
        IF (Coord.EQ.1) THEN
          CharaTemp = Heat_A / 200. / (2 - 1./Adsorption%CrystalIndx(iSurf))
        ELSE
          CharaTemp = Heat_A / 200. / (2 - 1./SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%nInterAtom)
        END IF
        ! calculate partition function of first particle bound on surface
        CALL PartitionFuncSurf(iSpec, WallTemp, VarPartitionFuncWall1,CharaTemp,PartBoundID)
        IF (Coord_ReactP(iReact).EQ.1) THEN
          CharaTemp = Heat_B / 200. / (2 - 1./Adsorption%CrystalIndx(iSurf))
        ELSE
          CharaTemp = Heat_B / 200. / (2 - 1./SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_ReactP(iReact))%nInterAtom)
        END IF
        ! calculate partition function of second particle bound on surface
        CALL PartitionFuncSurf(jSpec, WallTemp, VarPartitionFuncWall2,CharaTemp,PartBoundID)
        ! estimate partition function of activated complex
        CALL PartitionFuncAct_recomb(iSpec,jSpec, WallTemp, VarPartitionFuncAct, &
                    Adsorption%DensSurfAtoms(iSurf)/Adsorption%AreaIncrease(iSurf))
        ! transition state theory to estimate pre-exponential factor
        nu_react = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct)/(VarPartitionFuncWall1*VarPartitionFuncWall2)
        rate = nu_react * exp(-E_d/(WallTemp))
        P_actual_react(iReact) = rate * dt
#if (PP_TimeDiscMethod==42)
! sample reservoir reaction probabilities
        iSampleReact = iReact + Adsorption%DissNum + 1
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          IF (rate*dt.GT.1) THEN
            Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(iSampleReact) = &
                Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(iSampleReact) + 1.
!            Adsorption%AdsorpInfo(ProdSpec1)%MeanProbDes = Adsorption%AdsorpInfo(ProdSpec1)%MeanProbDes + 1.
          ELSE
            Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(iSampleReact) = &
                Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(iSampleReact) + P_actual_react(iReact)
!            Adsorption%AdsorpInfo(ProdSpec1)%MeanProbDes = Adsorption%AdsorpInfo(ProdSpec1)%MeanProbDes + P_actual_react(iReact)
          END IF
        END IF
        loc_SurfActE(iSampleReact) = E_d
        Adsorption%AdsorpReactInfo(iSpec)%MeanSurfActE(iSampleReact) = &
            Adsorption%AdsorpReactInfo(iSpec)%MeanSurfActE(iSampleReact) + E_d
        Adsorption%AdsorpReactInfo(iSpec)%SurfReactCount(iSampleReact) = &
            Adsorption%AdsorpReactInfo(iSpec)%SurfReactCount(iSampleReact) + 1
#endif
      END IF
    END DO
    ReactNum_run = ReactNum_run + Adsorption%RecombNum

    !-------------------------------------------------------------------------------------------------------------------------------
    ! sort Neighbours to coordinations for search of two random neighbour positions from particle position
    ! for dissociation and diffusion
    !-------------------------------------------------------------------------------------------------------------------------------
    IF (Adsorption%DissNum.GT.0 .AND. SpecDSMC(iSpec)%InterID.EQ.2) THEN
      n_Neigh(:) = 0
      n_empty_Neigh(:) = 0
      ALLOCATE(NeighbourID(1:3,1:SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%nNeighbours))
      ALLOCATE(React_NeighbourID(1:3,1:SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%nNeighbours))

      DO i = 1,SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%nNeighbours
        Coord3 = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%NeighSite(Surfpos,i)
        IF ((Coord3.GT.0).AND.SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%IsNearestNeigh(Surfpos,i)) THEN
          IF ( (SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord3)%Species( &
                SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%NeighPos(Surfpos,i)) &
                .EQ.0) ) THEN
            n_empty_Neigh(Coord3) = n_empty_Neigh(Coord3) + 1
            NeighbourID(Coord3,n_empty_Neigh(Coord3)) = i
          END IF
          n_Neigh(Coord3) = n_Neigh(Coord3) + 1
          React_NeighbourID(Coord3,n_Neigh(Coord3)) = i
        END IF
      END DO
    END IF
    !-------------------------------------------------------------------------------------------------------------------------------
    ! calculate probability for dissociation
    !-------------------------------------------------------------------------------------------------------------------------------
    ! (both products stay adsorbed --> endothermic)
    ! This routine looks if neirest neighbours are empty and then calculates dissociation probability
    IF (SpecDSMC(iSpec)%InterID.EQ.2) THEN
    DO iReact = 1,Adsorption%DissNum
      Prod_Spec1 = Adsorption%DissocReact(1,iReact,iSpec)
      Prod_Spec2 = Adsorption%DissocReact(2,iReact,iSpec)
      IF ((Prod_Spec1.NE.0) .AND. (Prod_Spec2.NE.0)) THEN !if 2 resulting species, dissociation possible
        Coord_Product(1,iReact+ReactNum_run) = Adsorption%Coordination(PartBoundID,Prod_Spec1)
        Coord_Product(2,iReact+ReactNum_run) = Adsorption%Coordination(PartBoundID,Prod_Spec2)
        IF ((Coord_Product(1,iReact+ReactNum_run).EQ.Coord_Product(2,iReact+ReactNum_run))&
              .AND.(nSitesRemain(Coord_Product(1,iReact+ReactNum_run)).LT.2)) CYCLE
        IF ((Coord_Product(1,iReact+ReactNum_run).NE.Coord_Product(2,iReact+ReactNum_run))&
              .AND.((nSitesRemain(Coord_Product(1,iReact+ReactNum_run)).LT.1)&
              .OR.(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)).LT.1))) CYCLE
        jCoord = Adsorption%Coordination(PartBoundID,Prod_Spec1)
        kCoord = Adsorption%Coordination(PartBoundID,Prod_Spec2)
        Neighpos_j = 0
        Neighpos_k = 0
        !At least one of resulting species has same coordination as surfaceposition of dissociating particle
        IF ( (jCoord.EQ.Coord) .OR. (kCoord.EQ.Coord) ) THEN
          IF (jCoord.EQ.Coord) THEN
            Neighpos_j = Surfpos
            IF (n_empty_Neigh(kCoord).GT.0) THEN
              ! assign availiable neighbour position for k
              CALL RANDOM_NUMBER(RanNum)
              chosen_Neigh_k = 1 + INT(REAL(n_empty_Neigh(kCoord))*RanNum)
              Neighpos_k = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%NeighPos(Surfpos,&
                                                                                      NeighbourID(kCoord,chosen_Neigh_k))
              IDRearrange = NeighbourID(kCoord,chosen_Neigh_k)
              NeighbourID(kCoord,chosen_Neigh_k) = NeighbourID(kCoord,n_empty_Neigh(kCoord))
              NeighbourID(kCoord,n_empty_Neigh(kCoord)) = IDRearrange
            END IF
          ELSE
            Neighpos_k = Surfpos
            IF (n_empty_Neigh(jCoord).GT.0) THEN
              ! assign availiable neighbour position for j
              CALL RANDOM_NUMBER(RanNum)
              chosen_Neigh_j = 1 + INT(REAL(n_empty_Neigh(jCoord))*RanNum)
              Neighpos_j = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%NeighPos(Surfpos,&
                                                                                      NeighbourID(jCoord,chosen_Neigh_j))
              IDRearrange = NeighbourID(jCoord,chosen_Neigh_j)
              NeighbourID(jCoord,chosen_Neigh_j) = NeighbourID(jCoord,n_empty_Neigh(jCoord))
              NeighbourID(jCoord,n_empty_Neigh(jCoord)) = IDRearrange
            END IF
          END IF
        !NEITHER of resulting species have same coordination as surfaceposition of dissociating particle
        ELSE
          IF ( (jCoord.EQ.kCoord) .AND. (n_empty_Neigh(jCoord).GT.1) ) THEN !resulting species are equal
            ! assign availiable neighbour positions
            CALL RANDOM_NUMBER(RanNum)
            chosen_Neigh_j = 1 + INT(REAL(n_empty_Neigh(jCoord))*RanNum)
            Neighpos_j = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%NeighPos(Surfpos,&
                                                                                    NeighbourID(jCoord,chosen_Neigh_j))
            IDRearrange = NeighbourID(jCoord,chosen_Neigh_j)
            NeighbourID(jCoord,chosen_Neigh_j) = NeighbourID(jCoord,n_empty_Neigh(jCoord))
            NeighbourID(jCoord,n_empty_Neigh(jCoord)) = IDRearrange
            CALL RANDOM_NUMBER(RanNum)
            chosen_Neigh_k = 1 + INT(REAL(n_empty_Neigh(kCoord)-1)*RanNum)
            Neighpos_k = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%NeighPos(Surfpos,&
                                                                                    NeighbourID(kCoord,chosen_Neigh_k))
            IDRearrange = NeighbourID(kCoord,chosen_Neigh_j)
            NeighbourID(kCoord,chosen_Neigh_k) = NeighbourID(kCoord,n_empty_Neigh(kCoord)-1)
            NeighbourID(kCoord,n_empty_Neigh(kCoord)-1) = IDRearrange
          ELSE IF ( (jCoord.NE.kCoord) .AND. (n_empty_Neigh(jCoord).GT.0) .AND. (n_empty_Neigh(kCoord).GT.0) ) THEN
            ! assign availiable neighbour positions
            CALL RANDOM_NUMBER(RanNum)
            chosen_Neigh_j = 1 + INT(REAL(n_empty_Neigh(jCoord))*RanNum)
            Neighpos_j = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%NeighPos(Surfpos,&
                                                                                    NeighbourID(jCoord,chosen_Neigh_j))
            IDRearrange = NeighbourID(jCoord,chosen_Neigh_j)
            NeighbourID(jCoord,chosen_Neigh_j) = NeighbourID(jCoord,n_empty_Neigh(jCoord))
            NeighbourID(jCoord,n_empty_Neigh(jCoord)) = IDRearrange
            CALL RANDOM_NUMBER(RanNum)
            chosen_Neigh_k = 1 + INT(REAL(n_empty_Neigh(kCoord))*RanNum)
            Neighpos_k = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%NeighPos(Surfpos,&
                                                                                    NeighbourID(kCoord,chosen_Neigh_k))
            IDRearrange = NeighbourID(kCoord,chosen_Neigh_j)
            NeighbourID(kCoord,chosen_Neigh_k) = NeighbourID(kCoord,n_empty_Neigh(kCoord))
            NeighbourID(kCoord,n_empty_Neigh(kCoord)) = IDRearrange
          END IF
        END IF
        Pos_Product(1,iReact+ReactNum_run) = Neighpos_j
        Pos_Product(2,iReact+ReactNum_run) = Neighpos_k
        IF ( (Neighpos_j.GT.0) .AND. (Neighpos_k.GT.0) ) THEN !both neighbour positions assigned
          ! remove adsorbate temporarily for ocupancy-check
          SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%Species(Surfpos) = 0
          ! both assigned neighbour positions empty
          IF ( (SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%Species(Neighpos_j).EQ.0) .AND. &
              (SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(kCoord)%Species(Neighpos_k).EQ.0) ) THEN
            Cell_Occupied = .FALSE.
            ! check for occupation with neirest Neighbours of both cells
            IF ( SpaceOccupied(iSurf,iSubSurf,jSubSurf,jCoord,Neighpos_j) &
                .OR. SpaceOccupied(iSurf,iSubSurf,jSubSurf,kCoord,Neighpos_k))  Cell_Occupied = .TRUE.
            ! removed adsorbate temporarily for occupancy-check
            SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%Species(Surfpos) = iSpec
            IF (Cell_Occupied) CYCLE ! Cycle if space is already occupied
            ! calculation of activation energy
            Heat_Product1 = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Prod_Spec1,Pos_Product(1,iReact+ReactNum_run),.TRUE.)
            Heat_Product2 = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Prod_Spec2,Pos_Product(2,iReact+ReactNum_run),.TRUE.)
            D_A = Adsorption%EDissBond(iReact,iSpec)
            D_Product1 = 0.
            D_Product2 = 0.
            E_d = Calc_E_Act(Heat_Product1,Heat_Product2,Heat_A,0.,D_Product1,D_Product2,D_A,0.)
            ! estimate vibrational temperatures of surface-particle bond
            IF (Coord.EQ.1) THEN
              CharaTemp = Heat_A / 200. / (2 - 1./Adsorption%CrystalIndx(iSurf))
            ELSE
              CharaTemp = Heat_A / 200. / (2 - 1./SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%nInterAtom)
            END IF
            ! calculate partition function of particles bound on surface
            CALL PartitionFuncSurf(iSpec, WallTemp, VarPartitionFuncWall1,CharaTemp,PartBoundID)
            ! estimate partition function of activated complex
            CALL PartitionFuncAct_dissoc(Prod_Spec1,Prod_Spec2, WallTemp, VarPartitionFuncAct)
            !CALL PartitionFuncAct_dissoc(iSpec,Prod_Spec1,Prod_Spec2, WallTemp, VarPartitionFuncAct, &
            !                      Adsorption%DensSurfAtoms(iSurf)/Adsorption%AreaIncrease(iSurf))
            ! transition state theory to estimate pre-exponential factor
            nu_react = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct/VarPartitionFuncWall1)
            rate = nu_react * exp(-E_d/(WallTemp))
            P_actual_react(ReactNum_run+iReact) = rate * dt
#if (PP_TimeDiscMethod==42)
! sample reservoir reaction probabilities
            iSampleReact = iReact + 1
            IF (.NOT.DSMC%ReservoirRateStatistic) THEN
              IF (rate*dt.GT.1) THEN
                Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(iSampleReact) = &
                    Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(iSampleReact) + 1.
              ELSE
                Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(iSampleReact) = &
                    Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(iSampleReact) + P_actual_react(ReactNum_run+iReact)
              END IF
            END IF
            loc_SurfActE(iSampleReact) = E_d
            Adsorption%AdsorpReactInfo(iSpec)%MeanSurfActE(iSampleReact) = &
                Adsorption%AdsorpReactInfo(iSpec)%MeanSurfActE(iSampleReact) + E_d
            Adsorption%AdsorpReactInfo(iSpec)%SurfReactCount(iSampleReact) = &
                Adsorption%AdsorpReactInfo(iSpec)%SurfReactCount(iSampleReact) + 1
#endif
          ELSE
            ! removed adsorbate temporarily for occupancy-check
            SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%Species(Surfpos) = iSpec
          END IF
        END IF
      END IF
    END DO
    END IF
    ReactNum_run = ReactNum_run + Adsorption%DissNum

    !-------------------------------------------------------------------------------------------------------------------------------
    ! Calculate propabilities for exchange reactions
    !-------------------------------------------------------------------------------------------------------------------------------
    ! maybe try later (if exothermic --> Delta_H < 0, product with lower adsorbheat desorbs else both stay adsorbed)
    exch_react_possible = .FALSE.
    DO iReact = 1,Adsorption%nDisPropReactions
      ! check if considered particle is one of defined reactants
      IF ( (Adsorption%ChemReactant(1,iReact+Adsorption%nDissocReactions).NE.iSpec) &
           .AND. (Adsorption%ChemReactant(2,iReact+Adsorption%nDissocReactions).NE.iSpec) ) THEN
        ! check if considered particle is one of defined products
        IF ( (Adsorption%ChemProduct(1,iReact+Adsorption%nDissocReactions).NE.iSpec) &
           .AND. (Adsorption%ChemProduct(2,iReact+Adsorption%nDissocReactions).NE.iSpec) ) CYCLE
      END IF
      ! Choose which reaction partner considered particle is
      ! -> defines which species is second reaction partner and if forward or reverse reaction
      ! ----------------------------------------------------------------------------------------------------------------------------
      IF (Adsorption%ChemReactant(1,iReact+Adsorption%nDissocReactions).EQ.iSpec) THEN ! -> first forward
      ! ----------------------------------------------------------------------------------------------------------------------------
        ReactDirForward = .TRUE.
        Coord_ReactP(iReact+ReactNum_run) = Adsorption%Coordination(PartBoundID,&
            Adsorption%ChemReactant(2,iReact+Adsorption%nDissocReactions))
        CALL RANDOM_NUMBER(RanNum)
        IF (Coord.EQ.Coord_ReactP(iReact+ReactNum_run)) THEN
          NeighID = 1 + INT((nSites(Coord_ReactP(iReact+ReactNum_run))-remainNum(Coord_ReactP(iReact+ReactNum_run))-1)*RanNum)
        ELSE
          NeighID = 1 + INT((nSites(Coord_ReactP(iReact+ReactNum_run))-remainNum(Coord_ReactP(iReact+ReactNum_run)))*RanNum)
        END IF
        Pos_ReactP(iReact+ReactNum_run) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                                          Coord_ReactP(iReact+ReactNum_run))%UsedSiteMap(NeighID)
        jSpec = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                Coord_ReactP(iReact+ReactNum_run))%Species(Pos_ReactP(iReact+ReactNum_run))
        ! continue only if species on chosen site is appropriate partner
        IF ((jSpec.EQ.Adsorption%ChemReactant(2,iReact+Adsorption%nDissocReactions)).AND.(jSpec.GT.0)) THEN
          exch_react_possible = .TRUE.
          Prod_Spec1 = Adsorption%ChemProduct(1,iReact+Adsorption%nDissocReactions)
          Coord_Product(1,iReact+ReactNum_run) = Adsorption%Coordination(PartBoundID,Prod_Spec1)
          Prod_Spec2 = Adsorption%ChemProduct(2,iReact+Adsorption%nDissocReactions)
          Coord_Product(2,iReact+ReactNum_run) = Adsorption%Coordination(PartBoundID,Prod_Spec2)
          IF ((Coord_Product(1,iReact+ReactNum_run).EQ.Coord_Product(2,iReact+ReactNum_run))&
               .AND.(nSitesRemain(Coord_Product(1,iReact+ReactNum_run)).LT.2)) CYCLE
          IF ((Coord_Product(1,iReact+ReactNum_run).NE.Coord_Product(2,iReact+ReactNum_run))&
               .AND.((nSitesRemain(Coord_Product(1,iReact+ReactNum_run)).LT.1)&
               .OR.(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)).LT.1))) CYCLE
          IF (Coord_Product(1,iReact+ReactNum_run).EQ.Coord_Product(2,iReact+ReactNum_run)) THEN
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,iReact+ReactNum_run)))/2.)*RanNum)
            Pos_Product(1,iReact+ReactNum_run) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                                              Coord_Product(1,iReact+ReactNum_run))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)))/2.)*RanNum &
                        + (REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)))/2.))
            Pos_Product(2,iReact+ReactNum_run) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                                              Coord_Product(2,iReact+ReactNum_run))%UsedSiteMap(NeighID)
          ELSE
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,iReact+ReactNum_run))))*RanNum)
            Pos_Product(1,iReact+ReactNum_run) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                                              Coord_Product(1,iReact+ReactNum_run))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run))))*RanNum)
            Pos_Product(2,iReact+ReactNum_run) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                                              Coord_Product(2,iReact+ReactNum_run))%UsedSiteMap(NeighID)
          END IF
          ! define dissociation bonds
          D_A = Adsorption%Reactant_DissBond_K(1,iReact+Adsorption%nDissocReactions)
          D_ReactP = Adsorption%Reactant_DissBond_K(2,iReact+Adsorption%nDissocReactions)
          D_Product1 = Adsorption%Product_DissBond_K(1,iReact+Adsorption%nDissocReactions)
          D_Product2 = Adsorption%Product_DissBond_K(2,iReact+Adsorption%nDissocReactions)
        END IF
      ! ----------------------------------------------------------------------------------------------------------------------------
      ELSE IF (Adsorption%ChemReactant(2,iReact+Adsorption%nDissocReactions).EQ.iSpec) THEN ! -> second forward
      ! ----------------------------------------------------------------------------------------------------------------------------
        ReactDirForward = .TRUE.
        Coord_ReactP(iReact+ReactNum_run) = Adsorption%Coordination(PartBoundID,&
            Adsorption%ChemReactant(1,iReact+Adsorption%nDissocReactions))
        CALL RANDOM_NUMBER(RanNum)
        IF (Coord.EQ.Coord_ReactP(iReact+ReactNum_run)) THEN
          NeighID = 1 + INT((nSites(Coord_ReactP(iReact+ReactNum_run))-remainNum(Coord_ReactP(iReact+ReactNum_run))-1)*RanNum)
        ELSE
          NeighID = 1 + INT((nSites(Coord_ReactP(iReact+ReactNum_run))-remainNum(Coord_ReactP(iReact+ReactNum_run)))*RanNum)
        END IF
        Pos_ReactP(iReact+ReactNum_run) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                                          Coord_ReactP(iReact+ReactNum_run))%UsedSiteMap(NeighID)
        jSpec = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                Coord_ReactP(iReact+ReactNum_run))%Species(Pos_ReactP(iReact+ReactNum_run))
        ! continue only if species on chosen site is appropriate partner
        IF ((jSpec.EQ.Adsorption%ChemReactant(1,iReact+Adsorption%nDissocReactions)).AND.(jSpec.GT.0)) THEN
          exch_react_possible = .TRUE.
          Prod_Spec1 = Adsorption%ChemProduct(1,iReact+Adsorption%nDissocReactions)
          Coord_Product(1,iReact+ReactNum_run) = Adsorption%Coordination(PartBoundID,Prod_Spec1)
          Prod_Spec2 = Adsorption%ChemProduct(2,iReact+Adsorption%nDissocReactions)
          Coord_Product(2,iReact+ReactNum_run) = Adsorption%Coordination(PartBoundID,Prod_Spec2)
          IF ((Coord_Product(1,iReact+ReactNum_run).EQ.Coord_Product(2,iReact+ReactNum_run))&
               .AND.(nSitesRemain(Coord_Product(1,iReact+ReactNum_run)).LT.2)) CYCLE
          IF ((Coord_Product(1,iReact+ReactNum_run).NE.Coord_Product(2,iReact+ReactNum_run))&
               .AND.((nSitesRemain(Coord_Product(1,iReact+ReactNum_run)).LT.1)&
               .OR.(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)).LT.1))) CYCLE
          IF (Coord_Product(1,iReact+ReactNum_run).EQ.Coord_Product(2,iReact+ReactNum_run)) THEN
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,iReact+ReactNum_run)))/2.)*RanNum)
            Pos_Product(1,iReact+ReactNum_run) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                                              Coord_Product(1,iReact+ReactNum_run))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)))/2.)*RanNum &
                        + (REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)))/2.))
            Pos_Product(2,iReact+ReactNum_run) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                                              Coord_Product(2,iReact+ReactNum_run))%UsedSiteMap(NeighID)
          ELSE
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,iReact+ReactNum_run))))*RanNum)
            Pos_Product(1,iReact+ReactNum_run) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                                              Coord_Product(1,iReact+ReactNum_run))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run))))*RanNum)
            Pos_Product(2,iReact+ReactNum_run) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                                              Coord_Product(2,iReact+ReactNum_run))%UsedSiteMap(NeighID)
          END IF
          ! define dissociation bonds
          D_A = Adsorption%Reactant_DissBond_K(2,iReact+Adsorption%nDissocReactions)
          D_ReactP = Adsorption%Reactant_DissBond_K(1,iReact+Adsorption%nDissocReactions)
          D_Product1 = Adsorption%Product_DissBond_K(1,iReact+Adsorption%nDissocReactions)
          D_Product2 = Adsorption%Product_DissBond_K(2,iReact+Adsorption%nDissocReactions)
        END IF
      ! ----------------------------------------------------------------------------------------------------------------------------
      ELSE IF (Adsorption%ChemProduct(1,iReact+Adsorption%nDissocReactions).EQ.iSpec) THEN ! -> first reverse
      ! ----------------------------------------------------------------------------------------------------------------------------
        ReactDirForward = .FALSE.
        Coord_ReactP(iReact+ReactNum_run) = Adsorption%Coordination(PartBoundID,&
            Adsorption%ChemProduct(2,iReact+Adsorption%nDissocReactions))
        CALL RANDOM_NUMBER(RanNum)
        IF (Coord.EQ.Coord_ReactP(iReact+ReactNum_run)) THEN
          NeighID = 1 + INT((nSites(Coord_ReactP(iReact+ReactNum_run))-remainNum(Coord_ReactP(iReact+ReactNum_run))-1)*RanNum)
        ELSE
          NeighID = 1 + INT((nSites(Coord_ReactP(iReact+ReactNum_run))-remainNum(Coord_ReactP(iReact+ReactNum_run)))*RanNum)
        END IF
        Pos_ReactP(iReact+ReactNum_run) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                                          Coord_ReactP(iReact+ReactNum_run))%UsedSiteMap(NeighID)
        jSpec = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                Coord_ReactP(iReact+ReactNum_run))%Species(Pos_ReactP(iReact+ReactNum_run))
        ! continue only if species on chosen site is appropriate partner
        IF ((jSpec.EQ.Adsorption%ChemProduct(2,iReact+Adsorption%nDissocReactions)).AND.(jSpec.GT.0)) THEN
          exch_react_possible = .TRUE.
          Prod_Spec1 = Adsorption%ChemReactant(1,iReact+Adsorption%nDissocReactions)
          Coord_Product(1,iReact+ReactNum_run) = Adsorption%Coordination(PartBoundID,Prod_Spec1)
          Prod_Spec2 = Adsorption%ChemReactant(2,iReact+Adsorption%nDissocReactions)
          Coord_Product(2,iReact+ReactNum_run) = Adsorption%Coordination(PartBoundID,Prod_Spec2)
          IF ((Coord_Product(1,iReact+ReactNum_run).EQ.Coord_Product(2,iReact+ReactNum_run))&
               .AND.(nSitesRemain(Coord_Product(1,iReact+ReactNum_run)).LT.2)) CYCLE
          IF ((Coord_Product(1,iReact+ReactNum_run).NE.Coord_Product(2,iReact+ReactNum_run))&
               .AND.((nSitesRemain(Coord_Product(1,iReact+ReactNum_run)).LT.1)&
               .OR.(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)).LT.1))) CYCLE
          IF (Coord_Product(1,iReact+ReactNum_run).EQ.Coord_Product(2,iReact+ReactNum_run)) THEN
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,iReact+ReactNum_run)))/2.)*RanNum)
            Pos_Product(1,iReact+ReactNum_run) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                                              Coord_Product(1,iReact+ReactNum_run))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)))/2.)*RanNum &
                        + (REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)))/2.))
            Pos_Product(2,iReact+ReactNum_run) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                                              Coord_Product(2,iReact+ReactNum_run))%UsedSiteMap(NeighID)
          ELSE
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,iReact+ReactNum_run))))*RanNum)
            Pos_Product(1,iReact+ReactNum_run) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                                              Coord_Product(1,iReact+ReactNum_run))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run))))*RanNum)
            Pos_Product(2,iReact+ReactNum_run) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                                              Coord_Product(2,iReact+ReactNum_run))%UsedSiteMap(NeighID)
          END IF
          ! define dissociation bonds
          D_A = Adsorption%Product_DissBond_K(1,iReact+Adsorption%nDissocReactions)
          D_ReactP = Adsorption%Product_DissBond_K(2,iReact+Adsorption%nDissocReactions)
          D_Product1 = Adsorption%Reactant_DissBond_K(1,iReact+Adsorption%nDissocReactions)
          D_Product2 = Adsorption%Reactant_DissBond_K(2,iReact+Adsorption%nDissocReactions)
        END IF
      ! ----------------------------------------------------------------------------------------------------------------------------
      ELSE IF (Adsorption%ChemProduct(2,iReact+Adsorption%nDissocReactions).EQ.iSpec) THEN ! -> second reverse
      ! ----------------------------------------------------------------------------------------------------------------------------
        ReactDirForward = .FALSE.
        Coord_ReactP(iReact+ReactNum_run) = Adsorption%Coordination(PartBoundID,&
            Adsorption%ChemProduct(1,iReact+Adsorption%nDissocReactions))
        CALL RANDOM_NUMBER(RanNum)
        IF (Coord.EQ.Coord_ReactP(iReact+ReactNum_run)) THEN
          NeighID = 1 + INT((nSites(Coord_ReactP(iReact+ReactNum_run))-remainNum(Coord_ReactP(iReact+ReactNum_run))-1)*RanNum)
        ELSE
          NeighID = 1 + INT((nSites(Coord_ReactP(iReact+ReactNum_run))-remainNum(Coord_ReactP(iReact+ReactNum_run)))*RanNum)
        END IF
        Pos_ReactP(iReact+ReactNum_run) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                                          Coord_ReactP(iReact+ReactNum_run))%UsedSiteMap(NeighID)
        jSpec = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                Coord_ReactP(iReact+ReactNum_run))%Species(Pos_ReactP(iReact+ReactNum_run))
        ! continue only if species on chosen site is appropriate partner
        IF ((jSpec.EQ.Adsorption%ChemProduct(1,iReact+Adsorption%nDissocReactions)).AND.(jSpec.GT.0)) THEN
          exch_react_possible = .TRUE.
          Prod_Spec1 = Adsorption%ChemReactant(1,iReact+Adsorption%nDissocReactions)
          Coord_Product(1,iReact+ReactNum_run) = Adsorption%Coordination(PartBoundID,Prod_Spec1)
          Prod_Spec2 = Adsorption%ChemReactant(2,iReact+Adsorption%nDissocReactions)
          Coord_Product(2,iReact+ReactNum_run) = Adsorption%Coordination(PartBoundID,Prod_Spec2)
          IF ((Coord_Product(1,iReact+ReactNum_run).EQ.Coord_Product(2,iReact+ReactNum_run))&
               .AND.(nSitesRemain(Coord_Product(1,iReact+ReactNum_run)).LT.2)) CYCLE
          IF ((Coord_Product(1,iReact+ReactNum_run).NE.Coord_Product(2,iReact+ReactNum_run))&
               .AND.((nSitesRemain(Coord_Product(1,iReact+ReactNum_run)).LT.1)&
               .OR.(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)).LT.1))) CYCLE
          IF (Coord_Product(1,iReact+ReactNum_run).EQ.Coord_Product(2,iReact+ReactNum_run)) THEN
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,iReact+ReactNum_run)))/2.)*RanNum)
            Pos_Product(1,iReact+ReactNum_run) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                                              Coord_Product(1,iReact+ReactNum_run))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)))/2.)*RanNum &
                        + (REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)))/2.))
            Pos_Product(2,iReact+ReactNum_run) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                                              Coord_Product(2,iReact+ReactNum_run))%UsedSiteMap(NeighID)
          ELSE
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,iReact+ReactNum_run))))*RanNum)
            Pos_Product(1,iReact+ReactNum_run) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                                              Coord_Product(1,iReact+ReactNum_run))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run))))*RanNum)
            Pos_Product(2,iReact+ReactNum_run) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                                              Coord_Product(2,iReact+ReactNum_run))%UsedSiteMap(NeighID)
          END IF
          ! define dissociation bonds
          D_A = Adsorption%Product_DissBond_K(2,iReact+Adsorption%nDissocReactions)
          D_ReactP = Adsorption%Product_DissBond_K(1,iReact+Adsorption%nDissocReactions)
          D_Product1 = Adsorption%Reactant_DissBond_K(1,iReact+Adsorption%nDissocReactions)
          D_Product2 = Adsorption%Reactant_DissBond_K(2,iReact+Adsorption%nDissocReactions)
        END IF
      END IF
      ! calculate probability if reaction possible
      IF (exch_react_possible) THEN
        ! define adsorption heats
        Heat_ReactP = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,jSpec,Pos_ReactP(iReact+ReactNum_run),.FALSE.)
        Heat_Product1 = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Prod_Spec1,Pos_Product(1,iReact+ReactNum_run),.TRUE.)
        Heat_Product2 = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Prod_Spec2,Pos_Product(2,iReact+ReactNum_run),.TRUE.)
        ! calculate activation energy
        E_d = Calc_E_Act(Heat_Product1,Heat_Product2,Heat_A,Heat_ReactP,D_Product1,D_Product2,D_A,D_ReactP)
        ! check diffusion barrier
        IF ((Coord.NE.3).OR.(Coord_ReactP(iReact+ReactNum_run).NE.3)) THEN
          IF ((Coord.NE.3).AND.(Coord_ReactP(iReact+ReactNum_run).NE.3)) THEN
            E_diff = (2* (0.1*Heat_A*0.1*Heat_ReactP))/(0.1*Heat_A + 0.1*Heat_ReactP)
          ELSE IF ((Coord.NE.3).AND.(Coord_ReactP(iReact+ReactNum_run).EQ.3)) THEN
            E_diff = 0.1*Heat_A
          ELSE IF ((Coord.EQ.3).AND.(Coord_ReactP(iReact+ReactNum_run).NE.3)) THEN
            E_diff = 0.1*Heat_ReactP
          ELSE
            E_diff = 0.
          END IF
        END IF
        ! if diffusion barrier is greater then reaction activation, then it is limiting --> activation energy
        IF (E_Diff.GT.E_d) E_d = E_diff
        ! estimate vibrational temperatures of surface-particle bond
        IF (Coord.EQ.1) THEN
          CharaTemp = Heat_A / 200. / (2 - 1./Adsorption%CrystalIndx(iSurf))
        ELSE
          CharaTemp = Heat_A / 200. / (2 - 1./SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%nInterAtom)
        END IF
        ! calculate partition function of first particle bound on surface
        CALL PartitionFuncSurf(iSpec, WallTemp, VarPartitionFuncWall1,CharaTemp,PartBoundID)
        IF (Coord_ReactP(iReact+ReactNum_run).EQ.1) THEN
          CharaTemp = Heat_ReactP / 200. / (2 - 1./Adsorption%CrystalIndx(iSurf))
        ELSE
          CharaTemp = Heat_ReactP / 200. &
                    / (2 - 1./SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_ReactP(iReact+ReactNum_run))%nInterAtom)
        END IF
        ! calculate partition function of second particle bound on surface
        CALL PartitionFuncSurf(jSpec, WallTemp, VarPartitionFuncWall2,CharaTemp,PartBoundID)
        ! estimate partition function of activated complex
        CALL PartitionFuncAct_exch(iSpec,jSpec, WallTemp, VarPartitionFuncAct, &
                    Adsorption%DensSurfAtoms(iSurf)/Adsorption%AreaIncrease(iSurf))
        ! transition state theory to estimate pre-exponential factor
        nu_react = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct/(VarPartitionFuncWall1*VarPartitionFuncWall2))
        rate = nu_react * exp(-E_d/(WallTemp))
        P_actual_react(iReact+ReactNum_run) = rate * dt
        IF (rate*dt.GT.1) THEN
          IF (ReactDirForward) THEN
            P_react_forward(iReact) = P_react_forward(iReact) + 1
          ELSE
            P_react_back(iReact) = P_react_back(iReact) + 1
          END IF
        ELSE
          IF (ReactDirForward) THEN
            P_react_forward(iReact) = P_react_forward(iReact) + P_actual_react(iReact+ReactNum_run)
          ELSE
            P_react_back(iReact) = P_react_back(iReact) + P_actual_react(iReact+ReactNum_run)
          END IF
        END IF
      END IF
    END DO
    ReactNum_run = ReactNum_run + Adsorption%nDisPropReactions

!    !-------------------------------------------------------------------------------------------------------------------------------
!    ! if empty neighbour sites are available choose one and calculate diffusion probability
!    !-------------------------------------------------------------------------------------------------------------------------------
!    IF (n_empty_Neigh(Coord).GT.0) THEN
!      ! choose random empty neighbour site from relevant neighbours
!      CALL RANDOM_NUMBER(RanNum)
!      react_Neigh = 1 + INT(n_empty_Neigh(Coord) * RanNum)
!      Coord3 = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%NeighSite(Surfpos,NeighbourID(Coord,react_Neigh))
!      Pos_ReactP(ReactNum_run+1) = &
!          SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%NeighPos(Surfpos,NeighbourID(Coord,react_Neigh))
!      ! calculate heat of adsorption for actual site
!      Heat_A = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,iSpec,Surfpos,.FALSE.)
!      ! remove considered particle from bondorder
!      CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,Coord,Surfpos,iSpec,.TRUE.)
!      !calculate heat of adsorption for diffusion site
!      Heat_B = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,iSpec,Pos_ReactP(ReactNum_run+1),.TRUE.)
!      ! add considered particle to bondorder again
!      CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,Coord,Surfpos,iSpec,.FALSE.)
!      ! calculate diffusion probability
!      P_actual_react(ReactNum_run+1) = exp(-(Heat_A - Heat_B)/WallTemp) / (1+exp(-(Heat_A - Heat_B)/Walltemp))
!    END IF

    SDEALLOCATE(NeighbourID)
    SDEALLOCATE(React_NeighbourID)

    !-------------------------------------------------------------------------------------------------------------------------------
    ! calculate molecular desorption probability
    !-------------------------------------------------------------------------------------------------------------------------------
    ! estimate vibrational temperature of surface-particle bond
    IF (Coord.EQ.1) THEN
      CharaTemp = Heat_A / 200. / (2 - 1./Adsorption%CrystalIndx(iSurf))
    ELSE
      CharaTemp = Heat_A / 200. / (2 - 1./SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%nInterAtom)
    END IF
    ! calculate partition function of particles bound on surface
    CALL PartitionFuncSurf(iSpec, WallTemp, VarPartitionFuncWall1,CharaTemp,PartBoundID)
    ! estimate partition function of activated complex
    CALL PartitionFuncAct(iSpec, WallTemp, VarPartitionFuncAct)!, Adsorption%DensSurfAtoms(iSurf)/Adsorption%AreaIncrease(iSurf))
    ! transition state theory to estimate pre-exponential factor
    nu_des = ((BoltzmannConst*Walltemp)/PlanckConst) * (VarPartitionFuncAct / VarPartitionFuncWall1)
    E_d = Calc_E_Act(0.,0.,Heat_A,0.,0.,0.,0.,0.)
    rate = nu_des * exp(-E_d/WallTemp)
    P_actual_des = rate * dt
#if (PP_TimeDiscMethod==42)
! sample reservoir reaction probabilities
    IF (.NOT.DSMC%ReservoirRateStatistic) THEN
      IF (rate*dt.GT.1) THEN
        Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(1) = Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(1) + 1.
        Adsorption%AdsorpInfo(iSpec)%MeanProbDes = Adsorption%AdsorpInfo(iSpec)%MeanProbDes + 1.
      ELSE
        Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(1) = Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(1) &
                                                               + P_actual_des
        Adsorption%AdsorpInfo(iSpec)%MeanProbDes = Adsorption%AdsorpInfo(iSpec)%MeanProbDes + P_actual_des
      END IF
    END IF
    loc_SurfActE(1) = E_d
    Adsorption%AdsorpReactInfo(iSpec)%MeanSurfActE(1) = Adsorption%AdsorpReactInfo(iSpec)%MeanSurfActE(1) + E_d
    Adsorption%AdsorpReactInfo(iSpec)%SurfReactCount(1) = Adsorption%AdsorpReactInfo(iSpec)%SurfReactCount(1) + 1
#endif
    ! initialize sum of all probabilities
    sum_probabilities = P_actual_des
    DO iReact = 1,ReactNum_run+1
      sum_probabilities = sum_probabilities + P_actual_react(iReact)
    END DO
    !-------------------------------------------------------------------------------------------------------------------------------
    ! choose which surface reaction takes place
    !-------------------------------------------------------------------------------------------------------------------------------
    IF (DSMC%ReservoirSurfaceRate) sum_probabilities = 0. !only probabilities are calculated without desorption taking place
    CALL RANDOM_NUMBER(RanNum)
    IF (sum_probabilities .GT. RanNum) THEN
      DO iReact = 1,ReactNum_run+1
        CALL RANDOM_NUMBER(RanNum)
        IF ((P_actual_react(iReact)/sum_probabilities).GT.RanNum) THEN
          IF (iReact.LE.(Adsorption%RecombNum)) THEN
            surf_react_case = 1
          ELSE IF ((iReact.GT.(Adsorption%RecombNum)).AND.(iReact.LE.(Adsorption%ReactNum))) THEN
            surf_react_case = 2
          ELSE IF ((iReact.GT.(Adsorption%ReactNum)).AND.(iReact.LT.ReactNum_run)) THEN
            surf_react_case = 3
          ELSE IF (iReact.EQ.(ReactNum_run+1)) THEN
            surf_react_case = 4
          END IF
          EXIT
        ELSE
          surf_react_case = 5
          sum_probabilities = sum_probabilities - P_actual_react(iReact)
        END IF
      END DO
      !-----------------------------------------------------------------------------------------------------------------------------
      SELECT CASE(surf_react_case)
      !-----------------------------------------------------------------------------------------------------------------------------
      CASE(1) ! association
      !-----------------------------------------------------------------------------------------------------------------------------
        ! update number of desorptions
        desorbnum(Adsorption%AssocReact(2,iReact,iSpec)) = desorbnum(Adsorption%AssocReact(2,iReact,iSpec)) + 1
        reactdesorbnum(Adsorption%AssocReact(2,iReact,iSpec)) = reactdesorbnum(Adsorption%AssocReact(2,iReact,iSpec)) + 1
        reactdesorbnum(Adsorption%AssocReact(1,iReact,iSpec)) = reactdesorbnum(Adsorption%AssocReact(1,iReact,iSpec)) - 1
        reactdesorbnum(iSpec) = reactdesorbnum(iSpec) - 1
#if (PP_TimeDiscMethod==42)
        IF (DSMC%ReservoirRateStatistic) THEN
          Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(iReact+Adsorption%DissNum+1) = &
              Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(iReact+Adsorption%DissNum+1) + 1
        END IF
        Adsorption%AdsorpReactInfo(iSpec)%ProperSurfActE(iReact+Adsorption%DissNum+1) = &
            Adsorption%AdsorpReactInfo(iSpec)%ProperSurfActE(iReact+Adsorption%DissNum+1)+loc_SurfActE(iReact+Adsorption%DissNum+1)
        Adsorption%AdsorpReactInfo(iSpec)%ProperSurfReactCount(iReact+Adsorption%DissNum+1) = &
            Adsorption%AdsorpReactInfo(iSpec)%ProperSurfReactCount(iReact+Adsorption%DissNum+1) + 1
#endif
        IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd))&
            .OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
          ! calculate Enthalpie of desorption and sample
          Heat_A = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,iSpec,-1,.FALSE.)
          Heat_B = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Adsorption%AssocReact(1,iReact,iSpec),-1,.FALSE.)
          Heat_AB = 0.
          D_AB = Adsorption%EDissBond(Adsorption%DissNum+iReact,iSpec)
          AdsorptionEnthalpie = (-(( Heat_AB -Heat_A -Heat_B ) + D_AB) * BoltzmannConst &
                          / REAL(SurfDistInfo(iSubSurf,jSubSurf,iSurf)%nSites(3))) &
                          * REAL(INT(Adsorption%DensSurfAtoms(iSurf) &
                          * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurf),8)) / Species(iSpec)%MacroParticleFactor
          !----  Sampling of energies
          SampWall(iSurf)%Adsorption(1,iSubSurf,jSubSurf) = SampWall(iSurf)%Adsorption(1,iSubSurf,jSubSurf) &
                                                                + AdsorptionEnthalpie * Species(iSpec)%MacroParticleFactor
          !NumDesorbLH(Adsorption%AssocReact(2,iReact,iSpec),iReact) = &
          !    NumDesorbLH(Adsorption%AssocReact(2,iReact,iSpec),iReact) + 1
          NumDesorbLH(iSpec,iReact) = NumDesorbLH(iSpec,iReact) + 1
        END IF
        ! remove adsorbate and update map
        CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,Coord,Surfpos,iSpec,.TRUE.)
        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = &
            SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1)
        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1) = Surfpos
        nSitesRemain(Coord) = nSitesRemain(Coord) + 1
        nSitesRemain(4) = nSitesRemain(4) + 1
        ! remove adsorbate reaction partner and update map
        DO PartnerID = nSitesRemain(Coord_ReactP(iReact))+1,nSites(Coord_ReactP(iReact))
        IF (SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_ReactP(iReact))%UsedSiteMap(PartnerID) &
            .EQ.Pos_ReactP(iReact)) THEN
          CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,Coord_ReactP(iReact), &
                             Pos_ReactP(iReact),Adsorption%AssocReact(1,iReact,iSpec),.TRUE.)
          SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_ReactP(iReact))%UsedSiteMap(PartnerID) = &
              SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_ReactP(iReact))%UsedSiteMap( &
              nSitesRemain(Coord_ReactP(iReact))+1)
          SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_ReactP(iReact))%UsedSiteMap( &
              nSitesRemain(Coord_ReactP(iReact))+1) = Pos_ReactP(iReact)
          nSitesRemain(Coord_ReactP(iReact)) = nSitesRemain(Coord_ReactP(iReact)) + 1
          nSitesRemain(4) = nSitesRemain(4) + 1
          ! additional increment to trace number
          trace = trace + 1
        END IF
        END DO !PartnerID
      !-----------------------------------------------------------------------------------------------------------------------------
      CASE(2) ! dissociation
      !-----------------------------------------------------------------------------------------------------------------------------
        ! update number of reactions
        DissocNum = iReact - Adsorption%RecombNum
        reactdesorbnum(Adsorption%DissocReact(1,DissocNum,iSpec)) = &
                reactdesorbnum(Adsorption%DissocReact(1,DissocNum,iSpec)) + 1
        reactdesorbnum(Adsorption%DissocReact(2,DissocNum,iSpec)) = &
                reactdesorbnum(Adsorption%DissocReact(2,DissocNum,iSpec)) + 1
        reactdesorbnum(iSpec) = reactdesorbnum(iSpec) - 1
#if (PP_TimeDiscMethod==42)
        IF (DSMC%ReservoirRateStatistic) THEN
          Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(DissocNum+1) = &
              Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(DissocNum+1) + 1
        END IF
        Adsorption%AdsorpReactInfo(iSpec)%ProperSurfActE(DissocNum+1) = &
            Adsorption%AdsorpReactInfo(iSpec)%ProperSurfActE(DissocNum+1) + loc_SurfActE(DissocNum+1)
        Adsorption%AdsorpReactInfo(iSpec)%ProperSurfReactCount(DissocNum+1) = &
            Adsorption%AdsorpReactInfo(iSpec)%ProperSurfReactCount(DissocNum+1) + 1
#endif
        IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd))&
            .OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
          ! calculate Enthalpie of desorption and sample
          Heat_A = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,iSpec,-1,.FALSE.)
          Heat_Product1 = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Adsorption%DissocReact(1,DissocNum,iSpec),-1,.FALSE.)
          Heat_Product2 = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Adsorption%DissocReact(2,DissocNum,iSpec),-1,.FALSE.)
          D_A = Adsorption%EDissBond(DissocNum,iSpec)
          AdsorptionEnthalpie = ((( Heat_A -Heat_Product1 -Heat_Product2 ) + D_A) * BoltzmannConst &
                          / REAL(SurfDistInfo(iSubSurf,jSubSurf,iSurf)%nSites(3))) &
                          * REAL(INT(Adsorption%DensSurfAtoms(iSurf) &
                          * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurf),8)) / Species(iSpec)%MacroParticleFactor
          !----  Sampling of energies
          SampWall(iSurf)%Adsorption(1,iSubSurf,jSubSurf) = SampWall(iSurf)%Adsorption(1,iSubSurf,jSubSurf) &
                                                                + AdsorptionEnthalpie * Species(iSpec)%MacroParticleFactor
        END IF
        ! remove adsorbate and update map
        CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,Coord,SurfPos,iSpec,.TRUE.)
        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = &
            SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1)
        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1) = Surfpos
        nSitesRemain(Coord) = nSitesRemain(Coord) + 1
        nSitesRemain(4) = nSitesRemain(4) + 1
        ! add resulting adsorbates from dissociation and update map
        ! first product
        jCoord = Coord_Product(1,iReact)
        DO i = 1,nSites(jCoord)
        IF (SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(i).EQ.Pos_Product(1,iReact)) THEN
          jSpec = Adsorption%DissocReact(1,DissocNum,iSpec)
          CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,jCoord,Pos_Product(1,iReact),jSpec,.FALSE.)
          SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(i) = &
              SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord))
!           SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord)) = Pos_Product(1,iReact)
          ! move Surfpos to MapID in remainNum segment
          SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord)) = &
              SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(nSites(jCoord)-remainNum(jCoord))
          SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(nSites(jCoord)-remainNum(jCoord)) = &
              Pos_Product(1,iReact)
          remainNum(jCoord) = remainNum(jCoord) + 1
          remainNum(4) = remainNum(4) + 1
          nSitesRemain(jCoord) = nSitesRemain(jCoord) - 1
          nSitesRemain(4) = nSitesRemain(4) - 1
          EXIT
        END IF
        END DO !i
        ! second product
        jCoord = Coord_Product(2,iReact)
        DO i = 1,nSites(jCoord)
        IF (SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(i).EQ.Pos_Product(2,iReact)) THEN
          jSpec = Adsorption%DissocReact(2,DissocNum,iSpec)
          CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,jCoord,Pos_Product(2,iReact),jSpec,.FALSE.)
          SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(i) = &
              SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord))
!           SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord)) = Pos_Product(2,iReact)
          ! move Surfpos to MapID in remainNum segment
          SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord)) = &
              SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(nSites(jCoord)-remainNum(jCoord))
          SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(nSites(jCoord)-remainNum(jCoord)) = &
              Pos_Product(2,iReact)
          remainNum(jCoord) = remainNum(jCoord) + 1
          remainNum(4) = remainNum(4) + 1
          nSitesRemain(jCoord) = nSitesRemain(jCoord) - 1
          nSitesRemain(4) = nSitesRemain(4) - 1
          EXIT
        END IF
        END DO !i
      !-----------------------------------------------------------------------------------------------------------------------------
      CASE(3) ! exchange reactions (disproportionation)
      !-----------------------------------------------------------------------------------------------------------------------------
        ExchNum = iReact - Adsorption%ReactNum + Adsorption%nDissocReactions
        ! choose if forward or reverse reaction and update number of reactions
        IF ( (Adsorption%ChemReactant(1,ExchNum).EQ.iSpec) .OR. (Adsorption%ChemReactant(2,ExchNum).EQ.iSpec) ) THEN
          reactdesorbnum(Adsorption%ChemReactant(1,ExchNum)) = reactdesorbnum(Adsorption%ChemReactant(1,ExchNum)) - 1
          reactdesorbnum(Adsorption%ChemReactant(2,ExchNum)) = reactdesorbnum(Adsorption%ChemReactant(2,ExchNum)) - 1
          reactdesorbnum(Adsorption%ChemProduct(1,ExchNum)) = reactdesorbnum(Adsorption%ChemProduct(1,ExchNum)) + 1
          reactdesorbnum(Adsorption%ChemProduct(2,ExchNum)) = reactdesorbnum(Adsorption%ChemProduct(2,ExchNum)) + 1
          ReactDirForward = .TRUE.
        ELSE
          reactdesorbnum(Adsorption%ChemReactant(1,ExchNum)) = reactdesorbnum(Adsorption%ChemReactant(1,ExchNum)) + 1
          reactdesorbnum(Adsorption%ChemReactant(2,ExchNum)) = reactdesorbnum(Adsorption%ChemReactant(2,ExchNum)) + 1
          reactdesorbnum(Adsorption%ChemProduct(1,ExchNum)) = reactdesorbnum(Adsorption%ChemProduct(1,ExchNum)) - 1
          reactdesorbnum(Adsorption%ChemProduct(2,ExchNum)) = reactdesorbnum(Adsorption%ChemProduct(2,ExchNum)) - 1
          ReactDirForward = .FALSE.
        END IF
        IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd))&
            .OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
          ! calculate Enthalpie of reaction and sample
          Heat_A = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Adsorption%ChemReactant(1,ExchNum),-1,.FALSE.)
          Heat_ReactP = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Adsorption%ChemReactant(2,ExchNum),-1,.FALSE.)
          Heat_Product1 = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Adsorption%ChemProduct(1,ExchNum),-1,.FALSE.)
          Heat_Product2 = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Adsorption%ChemProduct(2,ExchNum),-1,.FALSE.)
          D_A = Adsorption%Reactant_DissBond_K(1,ExchNum)
          D_ReactP = Adsorption%Reactant_DissBond_K(2,ExchNum)
          D_Product1 = Adsorption%Product_DissBond_K(1,ExchNum)
          D_Product2 = Adsorption%Product_DissBond_K(2,ExchNum)
          AdsorptionEnthalpie = (-(( Heat_A +Heat_ReactP -Heat_Product1 -Heat_Product2 ) &
                          + (D_A +D_ReactP -D_Product1 -D_Product2)) * BoltzmannConst &
                          / REAL(SurfDistInfo(iSubSurf,jSubSurf,iSurf)%nSites(3))) &
                          * REAL(INT(Adsorption%DensSurfAtoms(iSurf) &
                          * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurf),8)) / Species(iSpec)%MacroParticleFactor
          IF (ReactDirForward) AdsorptionEnthalpie = -AdsorptionEnthalpie
          !----  Sampling of energies
          SampWall(iSurf)%Adsorption(1,iSubSurf,jSubSurf) = SampWall(iSurf)%Adsorption(1,iSubSurf,jSubSurf) &
                                                                + AdsorptionEnthalpie * Species(iSpec)%MacroParticleFactor
        END IF
        ! remove first adsorbate and update map
        CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,Coord,Surfpos,iSpec,.TRUE.)
        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = &
            SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1)
        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1) = Surfpos
        nSitesRemain(Coord) = nSitesRemain(Coord) + 1
        nSitesRemain(4) = nSitesRemain(4) + 1
        ! remove second adsorbate (reaction partner) and update map
        IF (ReactDirForward) THEN
          IF (Adsorption%ChemReactant(1,ExchNum).EQ.iSpec) jSpec = Adsorption%ChemReactant(2,ExchNum)
          IF (Adsorption%ChemReactant(2,ExchNum).EQ.iSpec) jSpec = Adsorption%ChemReactant(1,ExchNum)
        ELSE
          IF (Adsorption%ChemProduct(1,ExchNum).EQ.iSpec) jSpec = Adsorption%ChemProduct(2,ExchNum)
          IF (Adsorption%ChemProduct(2,ExchNum).EQ.iSpec) jSpec = Adsorption%ChemProduct(1,ExchNum)
        END IF
        jCoord = Coord_ReactP(iReact)
        CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,jCoord,Pos_ReactP(iReact),jSpec,.TRUE.)
        DO i = nSitesRemain(jCoord)+1,nSites(jCoord)
        IF (SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(i).EQ.Pos_ReactP(iReact)) THEN
          SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(i) = SurfDistInfo( &
              iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord)+1)
          SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord) + 1) = Pos_ReactP(iReact)
          nSitesRemain(jCoord) = nSitesRemain(jCoord) + 1
          nSitesRemain(4) = nSitesRemain(4) + 1
          ! additional increment to trace number
          trace = trace + 1
          EXIT
        END IF
        END DO
        ! add resulting adsorbates from dissociation and update map
        ! first product
        jCoord = Coord_Product(1,iReact)
        DO i = 1,nSites(jCoord)
        IF (SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(i).EQ.Pos_Product(1,iReact)) THEN
          IF (ReactDirForward) THEN
            jSpec = Adsorption%ChemProduct(1,ExchNum)
          ELSE
            jSpec = Adsorption%ChemReactant(1,ExchNum)
          END IF
          CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,jCoord,Pos_Product(1,iReact),jSpec,.FALSE.)
          SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(i) = &
              SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord))
!           SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord)) = Pos_Product(1,iReact)
          ! move Surfpos to MapID in remainNum segment
          SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord)) = &
              SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(nSites(jCoord)-remainNum(jCoord))
          SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(nSites(jCoord)-remainNum(jCoord)) = &
              Pos_Product(1,iReact)
          remainNum(jCoord) = remainNum(jCoord) + 1
          remainNum(4) = remainNum(4) + 1
          nSitesRemain(jCoord) = nSitesRemain(jCoord) - 1
          nSitesRemain(4) = nSitesRemain(4) - 1
          EXIT
        END IF
        END DO !i
        ! second product
        jCoord = Coord_Product(2,iReact)
        DO i = 1,nSites(jCoord)
        IF (SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(i).EQ.Pos_Product(2,iReact)) THEN
          IF (ReactDirForward) THEN
            jSpec = Adsorption%ChemProduct(2,ExchNum)
          ELSE
            jSpec = Adsorption%ChemReactant(2,ExchNum)
          END IF
          CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,jCoord,Pos_Product(2,iReact),jSpec,.FALSE.)
          SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(i) = &
              SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord))
!           SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord)) = Pos_Product(2,iReact)
          ! move Surfpos to MapID in remainNum segment
          SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord)) = &
              SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(nSites(jCoord)-remainNum(jCoord))
          SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(jCoord)%UsedSiteMap(nSites(jCoord)-remainNum(jCoord)) = &
              Pos_Product(2,iReact)
          remainNum(jCoord) = remainNum(jCoord) + 1
          remainNum(4) = remainNum(4) + 1
          nSitesRemain(jCoord) = nSitesRemain(jCoord) - 1
          nSitesRemain(4) = nSitesRemain(4) - 1
          EXIT
        END IF
        END DO !i
      !-----------------------------------------------------------------------------------------------------------------------------
      CASE(4) !diffusion
      !-----------------------------------------------------------------------------------------------------------------------------
        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%Species(Surfpos) = 0
        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%Species(Pos_ReactP(ReactNum_run+1)) = iSpec
        DO j = 1,SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%nInterAtom
          Indx = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%BondAtomIndx(Surfpos,j)
          Indy = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%BondAtomIndy(Surfpos,j)
          SurfDistInfo(iSubSurf,jSubSurf,iSurf)%SurfAtomBondOrder(iSpec,Indx,Indy) = &
              SurfDistInfo(iSubSurf,jSubSurf,iSurf)%SurfAtomBondOrder(iSpec,Indx,Indy) - 1
          Indx = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%BondAtomIndx(Pos_ReactP(ReactNum_run+1),j)
          Indy = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%BondAtomIndy(Pos_ReactP(ReactNum_run+1),j)
          SurfDistInfo(iSubSurf,jSubSurf,iSurf)%SurfAtomBondOrder(iSpec,Indx,Indy) = &
              SurfDistInfo(iSubSurf,jSubSurf,iSurf)%SurfAtomBondOrder(iSpec,Indx,Indy) + 1
        END DO
        ! move adsorbate to empty site and update map
        DO i = 1,nSitesRemain(Coord)
          IF (SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(i) &
              .EQ.Pos_ReactP(ReactNum_run+1)) THEN
            SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(i) = Surfpos
!             SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = Pos_ReactP(ReactNum_run+1)
            ! move Surfpos to MapID in remainNum segment
            SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = &
                SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(nSites(Coord)-remainNum(Coord))
            SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(nSites(Coord)-remainNum(Coord)) &
                = Pos_ReactP(ReactNum_run+1)
            EXIT
          END IF
        END DO
        remainNum(Coord) = remainNum(Coord) + 1
        remainNum(4) = remainNum(4) + 1
        ! additional increment to trace number
        trace = trace + 1
      !-----------------------------------------------------------------------------------------------------------------------------
      CASE(5) !(desorption)
      !-----------------------------------------------------------------------------------------------------------------------------
        desorbnum(iSpec) = desorbnum(iSpec) + 1
#if (PP_TimeDiscMethod==42)
        IF (DSMC%ReservoirRateStatistic) THEN
          Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(1) = Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(1) + 1
        END IF
        Adsorption%AdsorpReactInfo(iSpec)%ProperSurfActE(1) = &
            Adsorption%AdsorpReactInfo(iSpec)%ProperSurfActE(1) + loc_SurfActE(1)
        Adsorption%AdsorpReactInfo(iSpec)%ProperSurfReactCount(1) = &
            Adsorption%AdsorpReactInfo(iSpec)%ProperSurfReactCount(1) + 1
#endif
        IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd))&
            .OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
          ! calculate Enthalpie of desorption and sample
          AdsorptionEnthalpie = (Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,iSpec,-1,.TRUE.) * BoltzmannConst &
                          / REAL(SurfDistInfo(iSubSurf,jSubSurf,iSurf)%nSites(3))) &
                          * REAL(INT(Adsorption%DensSurfAtoms(iSurf) &
                          * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurf),8)) / Species(iSpec)%MacroParticleFactor
          !----  Sampling of energies
          SampWall(iSurf)%Adsorption(1,iSubSurf,jSubSurf) = SampWall(iSurf)%Adsorption(1,iSubSurf,jSubSurf) &
                                                                  + AdsorptionEnthalpie * Species(iSpec)%MacroParticleFactor
        END IF

        CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,Coord,Surfpos,iSpec,.TRUE.)
        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = &
            SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1)
        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1) = Surfpos
        nSitesRemain(Coord) = nSitesRemain(Coord) + 1
        nSitesRemain(4) = nSitesRemain(4) + 1
      END SELECT
      !-----------------------------------------------------------------------------------------------------------------------------
    ELSE ! adsorbate remains on surface on same site
!       SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = &
!           SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(nSites(Coord)-remainNum(Coord))
!       SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(nSites(Coord)-remainNum(Coord)) = Surfpos
      remainNum(Coord) = remainNum(Coord) + 1
      remainNum(4) = remainNum(4) + 1
    END IF
    ! update number of adsorbates
    adsorbnum(:) = nSites(:) - nSitesRemain(:)
    SurfDistInfo(iSubSurf,jSubSurf,iSurf)%SitesRemain(:) = nSitesRemain(1:3)
    ! increment number of considered surface particles
    trace = trace + 1
  END DO ! trace < traceNum
  END IF ! number of adsobred particles on surface > 0
  !---------------------------------------------------------------------------------------------------------------------------------
  ! update and sample relevant values (desorbed particles, adsorption heat, ...)
  !---------------------------------------------------------------------------------------------------------------------------------
  DO iSpec = 1,nSpecies
    numSites = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%nSites(3)
    ! calculate number of desorbed particles for each species
    SurfDistInfo(iSubSurf,jSubSurf,iSurf)%desorbnum_tmp(iSpec) = &
                        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%desorbnum_tmp(iSpec) &
                        + ((REAL(desorbnum(iSpec)) / REAL(numSites)) &
                        * REAL(INT(Adsorption%DensSurfAtoms(iSurf) &
                        * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurf),8)) / Species(iSpec)%MacroParticleFactor)
    ! calculate number of desorbing simulation particles (round to integer)
    Adsorption%SumDesorbPart(iSubSurf,jSubSurf,iSurf,iSpec) = &
                        INT(SurfDistInfo(iSubSurf,jSubSurf,iSurf)%desorbnum_tmp(iSpec))
    ! Adjust tracking desorbing simulation particles
    SurfDistInfo(iSubSurf,jSubSurf,iSurf)%desorbnum_tmp(iSpec) = &
                        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%desorbnum_tmp(iSpec) &
                        - REAL(Adsorption%SumDesorbPart(iSubSurf,jSubSurf,iSurf,iSpec))
    ! calculate number of reacted particles for each species
    SurfDistInfo(iSubSurf,jSubSurf,iSurf)%reactnum_tmp(iSpec) = &
                        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%reactnum_tmp(iSpec) &
                        + ((REAL(reactdesorbnum(iSpec)) / REAL(numSites)) &
                        * REAL(INT(Adsorption%DensSurfAtoms(iSurf) &
                        * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurf),8)) / Species(iSpec)%MacroParticleFactor)
    ! calculate number of reacting simulation particles on surface (round to integer)
    Adsorption%SumReactPart(iSubSurf,jSubSurf,iSurf,iSpec) = &
                        INT(SurfDistInfo(iSubSurf,jSubSurf,iSurf)%reactnum_tmp(iSpec))
    ! Adjust tracking reacting simulation particles
    SurfDistInfo(iSubSurf,jSubSurf,iSurf)%reactnum_tmp(iSpec) = &
                        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%reactnum_tmp(iSpec) &
                        - REAL(Adsorption%SumReactPart(iSubSurf,jSubSurf,iSurf,iSpec))
    ! Sample vibrational energies
    ! due to reaction a part of energies can be transformed into other vibrational groundstates and the change is sampled
    ! the energies of the emitted particles are sampled in surfflux routine (particle_emission.f90) because calculated there
    IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
      IF (SpecDSMC(iSpec)%InterID.EQ.2) THEN
        SampleParts = Adsorption%SumDesorbPart(iSubSurf,jSubSurf,iSurf,iSpec) &
                    - Adsorption%SumReactPart(iSubSurf,jSubSurf,iSurf,iSpec)
        IF (SampleParts.GT.0) THEN
          IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
            iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
            ALLOCATE(RanNumPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF),VibQuantWallPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF))
            CALL RANDOM_NUMBER(RanNumPoly)
            VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
            DO WHILE (ALL(VibQuantWallPoly.GE.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF))
              CALL RANDOM_NUMBER(RanNumPoly)
              VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
            END DO
          END IF
          EvibNew = 0. ! later added in particle emission when particle is emitted
          DO iPart = 1,SampleParts
            EvibOld = 0.
            EvibWall = 0.
            IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
              DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
                EvibOld = EvibOld + (VibQuantWallPoly(iDOF) + DSMC%GammaQuant) * BoltzmannConst &
                          * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) * Species(iSpec)%MacroParticleFactor
                EvibWall = EvibWall + VibQuantWallPoly(iDOF) * BoltzmannConst &
                          * SpecDSMC(iSpec)%CharaTVib * Species(iSpec)%MacroParticleFactor
              END DO
            ELSE
              VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(iSpec)%CharaTVib)
              DO WHILE (VibQuantWall.GE.SpecDSMC(iSpec)%MaxVibQuant)
                CALL RANDOM_NUMBER(RanNum)
                VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(iSpec)%CharaTVib)
              END DO
              EvibOld = (VibQuantWall + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(iSpec)%CharaTVib
              EvibWall = VibQuantWall * BoltzmannConst * SpecDSMC(iSpec)%CharaTVib
            END IF
            !----  Sampling for internal (vibrational) energy accommodation at walls for desorbed simulation gas particles
            SampWall(iSurf)%State(7,iSubSurf,jSubSurf) = SampWall(iSurf)%State(7,iSubSurf,jSubSurf) &
                                              + EvibOld * Species(iSpec)%MacroParticleFactor
            SampWall(iSurf)%State(8,iSubSurf,jSubSurf) = SampWall(iSurf)%State(8,iSubSurf,jSubSurf) &
                                              + EvibWall * Species(iSpec)%MacroParticleFactor
            SampWall(iSurf)%State(9,iSubSurf,jSubSurf) = SampWall(iSurf)%State(9,iSubSurf,jSubSurf) &
                                              + EvibNew * Species(iSpec)%MacroParticleFactor
          END DO
        END IF
        SDEALLOCATE(RanNumPoly)
        SDEALLOCATE(VibQuantWallPoly)
      END IF
      DO iReact = 1,Adsorption%RecombNum
        IF (NumDesorbLH(iSpec,iReact).GT.0) THEN
          SampWall(iSurf)%Reaction(iReact,iSpec,iSubSurf,jSubSurf) = &
              SampWall(iSurf)%Reaction(iReact,iSpec,iSubSurf,jSubSurf) &
              + ((REAL(NumDesorbLH(iSpec,iReact)) / REAL(numSites)) * REAL(INT(Adsorption%DensSurfAtoms(iSurf) &
              * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurf),8)) / Species(iSpec)%MacroParticleFactor)
        END IF
      END DO
    END IF
    !-------------------------------------------------------------------------------------------------------------------------------
    ! analyze rate data
    !-------------------------------------------------------------------------------------------------------------------------------
#if (PP_TimeDiscMethod==42)
    Adsorption%AdsorpInfo(iSpec)%NumOfDes = Adsorption%AdsorpInfo(iSpec)%NumOfDes &
                                          + Adsorption%SumDesorbPart(iSubSurf,jSubSurf,iSurf,iSpec)
    !Adsorption%AdsorpInfo(iSpec)%MeanProbDes = Adsorption%AdsorpInfo(iSpec)%MeanProbDes &
    !                                          + Adsorption%ProbDes(iSubSurf,jSubSurf,iSurf,iSpec)
#endif
  END DO ! nSpecies (analyze)
END DO ! nSurfSample
END DO ! nSurfSample
END DO ! SurfMesh%nSides

DEALLOCATE(desorbnum,adsorbnum,nSites,nSitesRemain,remainNum,adsorbates)
DEALLOCATE(P_actual_react,P_react_forward,P_react_back)
DEALLOCATE(Coord_ReactP,Pos_ReactP)

END SUBROUTINE SMCR_PartDesorb


SUBROUTINE SMCR_Diffusion()
!===================================================================================================================================
!> Calculation of diffusion on reconstructed surface with assumption of Quasi Chemical Approximation (QCA)
!===================================================================================================================================
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption, SurfDistInfo
USE MOD_SurfaceModel_Tools     ,ONLY: UpdateSurfPos, Calc_Adsorb_Heat
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, PartBound
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                          :: iSurf, SpecID, globSide, jSubSurf, iSubSurf
INTEGER                          :: Coord, nSites, nSitesRemain, i, AdsorbID
REAL                             :: WallTemp, Prob_diff, RanNum
REAL                             :: Heat_i, Heat_j, Heat_temp
INTEGER                          :: n_equal_site_Neigh, Surfpos, newpos
INTEGER , ALLOCATABLE            :: free_Neigh_pos(:)
!----------------------------------------------------------------------------------------------------------------------------------!

DO iSurf = 1,SurfMesh%nSides
DO jSubSurf = 1,nSurfSample
DO iSubSurf = 1,nSurfSample

  DO Coord = 1,3
    nSites = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%nSites(Coord)
    nSitesRemain = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%SitesRemain(Coord)

    ALLOCATE ( free_Neigh_pos(1:SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%nNeighbours))
    DO AdsorbID = nSitesRemain+1,nSites,1
      Surfpos = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(AdsorbID)
      SpecID = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%Species(Surfpos)
      n_equal_site_Neigh = 0
      free_Neigh_pos(:) = 0

      ! find free Neighbour positions of the same site-coordination
      DO i = 1,SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%nNeighbours
        IF (SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%NeighSite(Surfpos,i) .EQ. Coord) THEN
          IF ( (SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%Species( &
              SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%NeighPos(Surfpos,i)).EQ.0) ) THEN
            n_equal_site_Neigh = n_equal_site_Neigh + 1
            free_Neigh_pos(n_equal_site_Neigh) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%NeighPos(Surfpos,i)
          END IF
        END IF
      END DO

      ! calculate heat of adsorption for actual site then reduce bond order of bondatoms
      Heat_i = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,SpecID,Surfpos,.FALSE.)
      ! update surfatom bond order and species map
      CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,Coord,Surfpos,SpecID,.TRUE.)

      ! choose Neighbour position with highest heat of adsorption if adsorbate would move there
      Heat_j = 0.
      DO i = 1,n_equal_site_Neigh
        Heat_temp = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,SpecID,free_Neigh_pos(i),.TRUE.)
        IF (Heat_temp .GT. Heat_j) THEN
          Heat_j = Heat_temp
          newpos = free_Neigh_pos(i)
        END IF
      END DO

      ! only try to diffuse particle if unoccupied sites available
      IF (n_equal_site_Neigh .GE. 1) THEN
        globSide = Adsorption%SurfSideToGlobSideMap(iSurf)
        WallTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC(globSide)))
        Prob_diff = exp(-(Heat_i - Heat_j)/WallTemp) / (1+exp(-(Heat_i - Heat_j)/Walltemp))
        CALL RANDOM_NUMBER(RanNum)
        IF (Prob_diff.GT.RanNum) THEN
        ! move particle to new position and update map
          DO i = 1,nSitesRemain
            IF (SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(i).EQ.newpos) THEN
              SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(i) = Surfpos
              SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = newpos
            END IF
          END DO
        ELSE
          newpos = Surfpos
        END IF
      ELSE
        newpos = Surfpos
      END IF ! end if (n_equal_site_Neigh >= 1)

      ! update surfatom bond order and species map
      CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,Coord,newpos,SpecID,.FALSE.)
    END DO
    DEALLOCATE(free_Neigh_pos)
  END DO
END DO  !iSubSurf = 1,nSurfSample
END DO  !jSubSurf = 1,nSurfSample
END DO  !iSurf = 1,SurfMesh%nSides

END SUBROUTINE SMCR_Diffusion

END MODULE MOD_SMCR
