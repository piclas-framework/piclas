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

MODULE MOD_SMCR
!===================================================================================================================================
!> Main Routines of Surface Monte Carlo Replication
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

SUBROUTINE SMCR_PartAdsorb(subsurfxi,subsurfeta,SurfID,PartID,Norm_Velo,adsorption_case,ProductSpec,AdsorptionEnthalpy)
!===================================================================================================================================
!> Particle adsorption probability calculation for one impinging particle using a surface replication (SMCR) (surfacemodel = 3)
!===================================================================================================================================
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst!, PI
USE MOD_Particle_Vars          ,ONLY: PartSpecies, nSpecies, Species!, WriteMacroSurfaceValues
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_DSMC_Vars              ,ONLY: DSMC, SpecDSMC
USE MOD_SurfaceModel_Vars      ,ONLY: SurfDistInfo, Adsorption, SurfModel
USE MOD_SurfaceModel_Tools     ,ONLY: Calc_Adsorb_Heat, CalcDissRecombActEnergy
USE MOD_SurfaceModel_Tools     ,ONLY: CalcAdsorbReactProb, SpaceOccupied, UpdateSurfPos
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound, SurfMesh !, SampWall
!USE MOD_TimeDisc_Vars          ,ONLY: TEnd, time
#if (PP_TimeDiscMethod==42)
USE MOD_TimeDisc_Vars          ,ONLY: iter, dt
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)               :: subsurfxi,subsurfeta,SurfID,PartID
INTEGER,INTENT(OUT)              :: adsorption_case
INTEGER,INTENT(OUT)              :: ProductSpec(2)      !< ProductSpec(1) new ID of impacting particle (the old one can change)
                                                        !< ProductSpec(2) new ID of created or consumed partner
REAL   ,INTENT(OUT)              :: AdsorptionEnthalpy
REAL   ,INTENT(IN)               :: Norm_Velo
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                          :: SpecID, globSide, PartBoundID
INTEGER                          :: Coord, Coord2, i, AdsorbID, IDRearrange
INTEGER                          :: nSites, nSitesRemain
REAL                             :: WallTemp, RanNum
REAL , ALLOCATABLE               :: ProbAds(:)
INTEGER                          :: Surfpos, ReactNum
INTEGER                          :: jSpec, kSpec, jCoord, kCoord
REAL                             :: sum_probabilities
INTEGER , ALLOCATABLE            :: NeighbourID(:,:)!, NeighSpec(:)
INTEGER                          :: SiteSpec, Neighpos_j, Neighpos_k, chosen_Neigh_j, chosen_Neigh_k
INTEGER                          :: n_empty_Neigh(3), n_Neigh(3), adsorbates(nSpecies)
REAL                             :: E_a
REAL                             :: Heat_A, Heat_B, Heat_AB, D_AB
!REAL                             :: vel_norm, vel_coll, potential_pot, a_const, mu, surfmass, trapping_prob
!INTEGER                          :: nNeigh_trap
LOGICAL                          :: Cell_Occupied
REAL                             :: CharaTemp
INTEGER                          :: DissocReactID, RecombReactID
REAL                             :: E_d
REAL                             :: coverage_check
#if (PP_TimeDiscMethod==42)
! reservoir sample variables
INTEGER                          :: iSampleReact
REAL                             :: loc_AdsActE(1:Adsorption%ReactNum+1)
REAL                             :: loc_Adsnu(1:Adsorption%ReactNum+1)
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
! special TPD (temperature programmed desorption) surface temperature adjustment part
globSide = SurfMesh%SurfIDToSideID(SurfID)
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
  SpecID = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%Species(Surfpos)
  adsorbates(SpecID) = adsorbates(SpecID) + 1
END DO
END DO
! initialize probabilites
ALLOCATE( ProbAds(0:Adsorption%ReactNum))
ProbAds(:) = 0.
Cell_Occupied = .FALSE.
#if (PP_TimeDiscMethod==42)
loc_AdsActE = 0.
loc_Adsnu = 0.
#endif

! Choose Random surface site with species coordination
CALL RANDOM_NUMBER(RanNum)
SpecID = PartSpecies(PartID)
Coord = Adsorption%Coordination(PartBoundID,SpecID)
AdsorbID = 1 + INT(SurfDistInfo(subsurfxi,subsurfeta,SurfID)%nSites(Coord)*RanNum)
Surfpos = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%UsedSiteMap(AdsorbID)
SiteSpec = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%Species(Surfpos)

! sample the normal velo at surfaces for the incident temperature that is needed in adsorption probability calculation
Adsorption%SurfaceNormalVelo(subsurfxi,subsurfeta,SurfID,SpecID) = &
    Adsorption%SurfaceNormalVelo(subsurfxi,subsurfeta,SurfID,SpecID) + Norm_Velo
Adsorption%CollSpecPartNum(subsurfxi,subsurfeta,SurfID,SpecID) = Adsorption%CollSpecPartNum(subsurfxi,subsurfeta,SurfID,SpecID) + 1

!!-----------------------------------------------------------------------------------------------------------------------------------
!! calculate trapping probability (using hard cube collision with surface atom or adsorbate)
!!-----------------------------------------------------------------------------------------------------------------------------------
!! if site is empty nearest neighbour site can be occupied and this influences the collision cube mass
!IF (SiteSpec.EQ.0) THEN
!  ALLOCATE(NeighSpec(1:SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%nNeighbours))
!  NeighSpec = 0
!  nNeigh_trap = 0
!  DO i = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%nNeighbours
!    IF (.NOT.SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%IsNearestNeigh(Surfpos,i)) CYCLE
!    IF ((Coord.EQ.1) .AND. (SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%NeighSite(Surfpos,i).EQ.1)) CYCLE
!    DO Coord2 = 1,3
!      IF (Coord2.EQ.SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%NeighSite(Surfpos,i)) THEN
!        IF ( (SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord2)%Species( &
!              SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%NeighPos(Surfpos,i)) &
!              .NE.0) ) THEN
!              Cell_Occupied = .TRUE.
!              nNeigh_trap = nNeigh_trap + 1
!              NeighSpec(nNeigh_trap) = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord2)%Species( &
!              SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%NeighPos(Surfpos,i))
!        END IF
!      END IF
!    END DO
!  END DO
!END IF
!potential_pot = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfID,SpecID,Surfpos,.TRUE.)*BoltzmannConst
!vel_norm = - (( 2*(potential_pot)/Species(SpecID)%MassIC)**(0.5) + Norm_Velo)
!a_const = (Species(SpecID)%MassIC/(2*BoltzmannConst*WallTemp))**(0.5)
!IF (SiteSpec.EQ.0) THEN
!  IF (nNeigh_trap.EQ.0) THEN
!    surfmass = PartBound%SolidMassIC(PartBoundID) !Adsorption%SurfMassIC(SurfID)
!  ELSE
!    surfmass = 0.
!    DO i = 1,nNeigh_trap
!      surfmass = surfmass + Species(NeighSpec(i))%MassIC
!    END DO
!    surfmass = surfmass / nNeigh_trap
!  END IF
!  SDEALLOCATE(NeighSpec)
!  mu = Species(SpecID)%MassIC / surfmass
!ELSE
!  mu = Species(SpecID)%MassIC / (Species(SiteSpec)%MassIC)
!END IF
!vel_coll = 0.5 * ( (1+mu)*(2*(potential_pot/Species(SpecID)%MassIC))**(0.5) + (1-mu)*vel_norm )
!trapping_prob = abs(0.5 + 0.5*ERF(a_const*vel_coll) + ( EXP(-(a_const**2)*(vel_coll**2)) / (2*a_const*(vel_norm*PI**0.5)) ))
!IF (trapping_prob.GT.1.) trapping_prob = 1.
!#if (PP_TimeDiscMethod==42)
!SurfModel%Info(SpecID)%Accomodation = SurfModel%Info(SpecID)%Accomodation + trapping_prob
!#endif
!IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
!  SampWall(SurfID)%Accomodation(SpecID,subsurfxi,subsurfeta) = SampWall(SurfID)%Accomodation(SpecID,subsurfxi,subsurfeta) &
!                                                                + trapping_prob
!END IF

ProductSpec(1) = SpecID
ProductSpec(2) = 0
AdsorptionEnthalpy = 0.
CALL RANDOM_NUMBER(RanNum)
IF(RanNum.GE.PartBound%MomentumACC(PartBoundID)) THEN
  adsorption_case = 1
  RETURN
END IF
adsorption_case = 2

!-----------------------------------------------------------------------------------------------------------------------------------
! calculate probability for molecular adsorption
!-----------------------------------------------------------------------------------------------------------------------------------
ReactNum = 0
! check for occupation of nearest Neighbours or wether
IF ( (SiteSpec.EQ.0) .AND. (.NOT.SpaceOccupied(SurfID,subsurfxi,subsurfeta,Coord,SurfPos)) &
    .AND. (INT(SurfDistInfo(subsurfxi,subsurfeta,SurfID)%adsorbnum_tmp(SpecID)).LT.1)) THEN
  ! calculation of molecular adsorption probability
  E_a = 0.
  E_d = 0.1 * Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfID,SpecID,Surfpos,.TRUE.) * Boltzmannconst
#if (PP_TimeDiscMethod==42)
  iSampleReact = 1 + ReactNum
  ProbAds(ReactNum) = CalcAdsorbReactProb(1,ReactNum,PartID,SurfID &
      ,Adsorption%IncidentNormalVeloAtSurf(subsurfxi,subsurfeta,SurfID,SpecID),E_a,E_d &
      ,loc_ActE=loc_AdsActE(iSampleReact),loc_nu=loc_Adsnu(iSampleReact))
#else
  ProbAds(Reactnum) = CalcAdsorbReactProb(1,ReactNum,PartID,SurfID &
      ,Adsorption%IncidentNormalVeloAtSurf(subsurfxi,subsurfeta,SurfID,SpecID),E_a,E_d)
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
IF ((SpecDSMC(SpecID)%InterID.EQ.2)) THEN ! particle has to be at least di-atomic
DO ReactNum = 1,(Adsorption%DissNum)
  DissocReactID = ReactNum
  jSpec = Adsorption%DissocReact(1,DissocReactID,SpecID)
  kSpec = Adsorption%DissocReact(2,DissocReactID,SpecID)
  IF ((jSpec.NE.0) .AND. (kSpec.NE.0)) THEN !if 2 resulting species, dissociation possible
    ! in one of the last iterations the surface-coverage of one resulting species was saturated
    IF (INT(SurfDistInfo(subsurfxi,subsurfeta,SurfID)%adsorbnum_tmp(jSpec)).GE.1. .OR. &
        INT(SurfDistInfo(subsurfxi,subsurfeta,SurfID)%adsorbnum_tmp(kSpec)).GE.1. ) THEN
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
        IF (SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(Coord)%Species(SurfPos).NE.0) THEN
          ! calculation of activation energy
          Heat_A = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfID,jSpec,Neighpos_j,.TRUE.)
          Heat_B = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfID,kSpec,Neighpos_k,.TRUE.)
        ELSE
          CALL UpdateSurfPos(SurfID,subsurfxi,subsurfeta,Coord,Surfpos,SpecID,.FALSE.)
          ! calculation of activation energy
          Heat_A = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfID,jSpec,Neighpos_j,.TRUE.)
          Heat_B = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfID,kSpec,Neighpos_k,.TRUE.)
          ! remove bond order of respective surface atoms in the surfacelattice for molecule again
          CALL UpdateSurfPos(SurfID,subsurfxi,subsurfeta,Coord,Surfpos,SpecID,.TRUE.)
        END IF
        Heat_AB = 0. ! direct dissociative adsorption
        D_AB = Adsorption%EDissBond(ReactNum,SpecID)
        E_a = CalcDissRecombActEnergy(Heat_AB,Heat_A,Heat_B,D_AB,SpecID) * BoltzmannConst
        E_d = 0.0
        ! calculation of dissociative adsorption probability
#if (PP_TimeDiscMethod==42)
        iSampleReact = 1 + ReactNum
        ProbAds(ReactNum) = CalcAdsorbReactProb(2,ReactNum,PartID,SurfID &
            ,Adsorption%IncidentNormalVeloAtSurf(subsurfxi,subsurfeta,SurfID,SpecID),E_a,E_d&
            ,loc_ActE=loc_AdsActE(iSampleReact),loc_nu=loc_Adsnu(iSampleReact))
#else
        ProbAds(ReactNum) = CalcAdsorbReactProb(2,ReactNum,PartID,SurfID &
            ,Adsorption%IncidentNormalVeloAtSurf(subsurfxi,subsurfeta,SurfID,SpecID),E_a,E_d)
#endif
      END IF
    END IF !both neighbour positions assigned
  END IF
END DO
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! calculate probability for Eley-Rideal reaction
!-----------------------------------------------------------------------------------------------------------------------------------
DO ReactNum = Adsorption%DissNum+1,(Adsorption%ReactNum)
  RecombReactID = ReactNum-Adsorption%DissNum
  ! reaction partner
  jSpec = Adsorption%RecombReact(1,RecombReactID,SpecID)
  Neighpos_j = 0
  IF (jSpec.EQ.0) CYCLE
  IF (INT(SurfDistInfo(subsurfxi,subsurfeta,SurfID)%adsorbnum_tmp(jSpec)).LT.0 ) CYCLE
  coverage_check = Adsorption%Coverage(subsurfxi,subsurfeta,SurfID,jSpec) &
                 + SurfModel%SumAdsorbPart(subsurfxi,subsurfeta,SurfID,jSpec) &
                 * Species(jSpec)%MacroParticleFactor &
                 / REAL(INT(Adsorption%DensSurfAtoms(SurfID) * SurfMesh%SurfaceArea(subsurfxi,subsurfeta,SurfID),8))
  IF (coverage_check.LE.0.) CYCLE
  ! reaction results
  kSpec = Adsorption%RecombReact(2,RecombReactID,SpecID)
  ! Choose surface site with reactionpartner coordination
  jCoord = Adsorption%Coordination(PartBoundID,jSpec)
  IF (jCoord.EQ.Coord) THEN
    Neighpos_j = SurfPos
  ELSE
  !  CALL RANDOM_NUMBER(RanNum)
  !  AdsorbID = 1 + INT(SurfDistInfo(subsurfxi,subsurfeta,SurfID)%nSites(jCoord)*RanNum)
  !  Neighpos_j = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(jCoord)%UsedSiteMap(AdsorbID)
    IF ( n_Neigh(jCoord)-n_empty_Neigh(jCoord).GT.0 ) THEN
      CALL RANDOM_NUMBER(RanNum)
      chosen_Neigh_j = 1 + INT(n_Neigh(jCoord)*RanNum)
      Neighpos_j = SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(jCoord)%NeighPos(Surfpos,NeighbourID(jCoord,chosen_Neigh_j))
    END IF
  END IF
  IF ( (Neighpos_j.GT.0) ) THEN
    IF ( (SurfDistInfo(subsurfxi,subsurfeta,SurfID)%AdsMap(jCoord)%Species(Neighpos_j).EQ.jSpec) ) THEN
      ! calculation of activation energy
      Heat_A = 0.1*Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfID,SpecID,SurfPos,.TRUE.)
      Heat_B = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfID,jSpec,Neighpos_j,.FALSE.)
      Heat_AB = 0. ! direct associative desorption (recombination)
      D_AB = Adsorption%EDissBond(ReactNum,SpecID)
      E_a = CalcDissRecombActEnergy(Heat_AB,Heat_A,Heat_B,-D_AB,kSpec) * BoltzmannConst
      E_d = 0.
      ! estimate characteristic vibrational temperature of surface-particle bond
      CharaTemp = Heat_A / 200.
      ! calculation of ER-reaction probability
#if (PP_TimeDiscMethod==42)
      iSampleReact = 1 + ReactNum
      ProbAds(ReactNum) = CalcAdsorbReactProb(3,ReactNum,PartID,SurfID &
          ,Adsorption%IncidentNormalVeloAtSurf(subsurfxi,subsurfeta,SurfID,SpecID),E_a,E_d &
          ,loc_ActE=loc_AdsActE(iSampleReact),loc_nu=loc_Adsnu(iSampleReact))
#else
      ProbAds(ReactNum) = CalcAdsorbReactProb(3,ReactNum,PartID,SurfID &
          ,Adsorption%IncidentNormalVeloAtSurf(subsurfxi,subsurfeta,SurfID,SpecID),E_a,E_d)
#endif
    END IF
  END IF
END DO
!-----------------------------------------------------------------------------------------------------------------------------------

sum_probabilities = 0.
DO ReactNum = 0,Adsorption%ReactNum
  sum_probabilities = sum_probabilities + ProbAds(ReactNum)
END DO
SurfModel%Info(SpecID)%MeanProbAds = SurfModel%Info(SpecID)%MeanProbAds + sum_probabilities
SurfModel%Info(SpecID)%MeanProbAdsCount = SurfModel%Info(SpecID)%MeanProbAdsCount + 1

!-----------------------------------------------------------------------------------------------------------------------------------
! choose which adsorption type takes place
!-----------------------------------------------------------------------------------------------------------------------------------
CALL RANDOM_NUMBER(RanNum)
IF (sum_probabilities .GT. RanNum) THEN
  ! chose surface reaction case (0=inelastic scattering, 1=adsorption, 2=reaction (dissociation), 3=reaction (Eley-Rideal))
  adsorption_case = 3
  DO ReactNum = 0,(Adsorption%ReactNum)
    CALL RANDOM_NUMBER(RanNum)
    IF ((ProbAds(ReactNum)/sum_probabilities).GT.RanNum) THEN
      IF (ReactNum.EQ.0) THEN
        ! if molecular adsorption set output parameters
        ProductSpec(1) = -SpecID
        ProductSpec(2) = 0
      ELSE IF (ReactNum.GT.0 .AND. ReactNum.LE.Adsorption%DissNum) THEN
        ! if dissocciative adsorption set output parameters
        DissocReactID = ReactNum
        ProductSpec(1) = -Adsorption%DissocReact(1,DissocReactID,SpecID)
        ProductSpec(2) = -Adsorption%DissocReact(2,DissocReactID,SpecID)
        ! calculate reaction Enthalpy
        D_AB = Adsorption%EDissBond(ReactNum,SpecID)
        AdsorptionEnthalpy = D_AB * BoltzmannConst
      ELSE IF (ReactNum.GT.0 .AND. ReactNum.GT.Adsorption%DissNum) THEN
        ! if ER-reaction set output parameters
        RecombReactID = ReactNum - Adsorption%DissNum
        ProductSpec(1) = Adsorption%RecombReact(2,RecombReactID,SpecID)
        ProductSpec(2) = Adsorption%RecombReact(1,RecombReactID,SpecID)
        ! calculate adsorption Enthalpy
        D_AB = Adsorption%EDissBond(ReactNum,SpecID)
        AdsorptionEnthalpy = -D_AB * BoltzmannConst
      END IF
      EXIT
    END IF
    sum_probabilities = sum_probabilities - ProbAds(ReactNum)
  END DO
#if (PP_TimeDiscMethod==42)
  IF (adsorption_case.GT.2) THEN
    iSampleReact = 1 + ReactNum
    IF (DSMC%ReservoirRateStatistic) THEN
      SurfModel%ProperInfo(SpecID)%NumAdsReact(iSampleReact) = &
          SurfModel%ProperInfo(SpecID)%NumAdsReact(iSampleReact) + 1.
    END IF
    SurfModel%ProperInfo(SpecID)%ProperAdsActE(iSampleReact) = &
        SurfModel%ProperInfo(SpecID)%ProperAdsActE(iSampleReact) + loc_AdsActE(iSampleReact)
    SurfModel%ProperInfo(SpecID)%ProperAdsnu(iSampleReact) = &
        SurfModel%ProperInfo(SpecID)%ProperAdsnu(iSampleReact) + loc_Adsnu(iSampleReact)
    SurfModel%ProperInfo(SpecID)%ProperAdsReactCount(iSampleReact) = &
        SurfModel%ProperInfo(SpecID)%ProperAdsReactCount(iSampleReact) + 1
    SurfModel%ProperInfo(SpecID)%HeatFluxAdsCount(iSampleReact) = &
        SurfModel%ProperInfo(SpecID)%HeatFluxAdsCount(iSampleReact) + 1.
  END IF
#endif
END IF

IF (DSMC%ReservoirSurfaceRate) adsorption_case = 1

DEALLOCATE(ProbAds,NeighbourID)

END SUBROUTINE SMCR_PartAdsorb

SUBROUTINE SMCR_PartDesorb()
!===================================================================================================================================
!> Calculation of number of desorbing particles using surface replication (surfacemodel = 3)
!> Performing Surface Monte Carlo step (MCS)
!>  Loop over all particles on surfaces and calulate desorption/reaction probabities
!>  - choose respective process
!>  - update surface distribution
!===================================================================================================================================
USE MOD_Globals_Vars           ,ONLY: PlanckConst, BoltzmannConst
USE MOD_Particle_Vars          ,ONLY: nSpecies, Species, WriteMacroSurfaceValues
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_DSMC_Vars              ,ONLY: DSMC, SpecDSMC
USE MOD_SurfaceModel_Vars      ,ONLY: SurfDistInfo, Adsorption, surfmodel
USE MOD_SurfaceModel_Tools     ,ONLY: Calc_Adsorb_Heat, Calc_E_Act, CalcDissRecombActEnergy, SampleAdsorptionHeat
USE MOD_SurfaceModel_Tools     ,ONLY: SpaceOccupied, UpdateSurfPos, SurfaceHasModelNum
USE MOD_SurfaceModel_PartFunc  ,ONLY: PartitionFuncActDesorb, PartitionFuncSurf
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, PartBound, SampWall
USE MOD_TimeDisc_Vars          ,ONLY: dt
USE MOD_TimeDisc_Vars          ,ONLY: TEnd, time
#if (PP_TimeDiscMethod==42)
USE MOD_SurfaceModel_Analyze   ,ONLY: AnalyzeSurfRates
USE MOD_TimeDisc_Vars          ,ONLY: iter
#endif
#if USE_LOADBALANCE
USE MOD_Particle_Mesh_Vars     ,ONLY: PartSideToElem
USE MOD_LoadBalance_Vars       ,ONLY: nSurfacePartsPerElem, PerformLBSample
#endif /*USE_LOADBALANCE*/
#if USE_MPI
USE MOD_SurfaceModel_MPI       ,ONLY: ExchangeSurfDistInfo
#endif /*USE_MPI*/
!===================================================================================================================================
IMPLICIT NONE
!===================================================================================================================================
! Local variable declaration
INTEGER                           :: iSurf, iSpec, SpecID, globSide, jSubSurf, iSubSurf, PartBoundID
INTEGER                           :: iDissocReact, iRecombReact, ReactID
INTEGER                           :: Coord, Coord3, i, AdsorbID, numSites
REAL                              :: WallTemp, RanNum
REAL                              :: Heat_A, Heat_B, rate
INTEGER                           :: Surfpos
INTEGER                           :: trace, traceNum
INTEGER , ALLOCATABLE             :: desorbnum(:), reactdesorbnum(:)
INTEGER , ALLOCATABLE             :: adsorbnum(:)
INTEGER , ALLOCATABLE             :: nSites(:), nSitesRemain(:), remainNum(:), adsorbates(:)
LOGICAL                           :: Cell_Occupied
!---------- reaction variables
INTEGER                           :: n_Neigh(3), n_empty_Neigh(3), jSpec, iReact, PartnerID
INTEGER                           :: surf_react_case, NeighID
REAL                              :: E_d, nu_react, E_diff
REAL                              :: Heat_AB, D_AB, sum_probabilities
REAL                              :: VarPartitionFuncWall1, VarPartitionFuncWall2, VarPartitionFuncAct
REAL                              :: CharaTemp
INTEGER , ALLOCATABLE             :: Pos_ReactP(:), Coord_ReactP(:), Pos_Product(:,:), Coord_Product(:,:)
INTEGER , ALLOCATABLE             :: React_NeighbourID(:,:), NeighbourID(:,:)
REAL , ALLOCATABLE                :: ProbDes(:),P_react_forward(:),P_react_back(:)
REAL                              :: AdsorptionEnthalpy
!variables used for sampling of vibrational energies of reacting/desorbing particles
!INTEGER                           :: SampleParts, iPart
!INTEGER                           :: iPolyatMole, iDOF
!INTEGER                           :: VibQuantWall
!INTEGER, ALLOCATABLE              :: VibQuantWallPoly(:)
!REAL, ALLOCATABLE                 :: RanNumPoly(:)
!REAL                              :: EvibOld, EvibWall, EVibNew
! variables used for exchange reactions
REAL                              :: D_A, D_ReactP, D_Product1, D_Product2, Heat_ReactP, Heat_Product1, Heat_Product2
INTEGER                           :: Prod_Spec1, Prod_Spec2
LOGICAL                           :: ReactDirForward, exch_react_possible
INTEGER                           :: ExchNum, jCoord, kCoord, Neighpos_j, Neighpos_k, chosen_Neigh_j, chosen_Neigh_k
INTEGER                           :: IDRearrange
#if (PP_TimeDiscMethod==42)
! reservoir sample variables
INTEGER                           :: iSampleReact
REAL                              :: loc_SurfActE(0:Adsorption%ReactNum+Adsorption%NumOfExchReact)
REAL                              :: loc_Surfnu(0:Adsorption%ReactNum+Adsorption%NumOfExchReact)
#endif
INTEGER                           :: NumSurfReact(1:nSpecies,1:Adsorption%ReactNum)
#if USE_LOADBALANCE
INTEGER                           :: ElemID
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
IF (.NOT.SurfMesh%SurfOnProc) RETURN

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
ALLOCATE( P_react_forward(1:Adsorption%nExchReactions),&
          P_react_back(1:Adsorption%nExchReactions),&
          ProbDes(0:Adsorption%ReactNum+Adsorption%nExchReactions),&
          Coord_ReactP(1:Adsorption%ReactNum+Adsorption%nExchReactions),&
          Coord_Product(1:2,1:Adsorption%ReactNum+Adsorption%nExchReactions),&
          Pos_ReactP(1:Adsorption%ReactNum+Adsorption%nExchReactions),&
          Pos_Product(1:2,1:Adsorption%ReactNum+Adsorption%nExchReactions))

! sample energy of surfaces before desorption treatment
DO iSurf = 1,SurfMesh%nOutputSides
  IF (SurfaceHasModelNum(iSurf).NE.3) CYCLE
  DO jSubSurf = 1,nSurfSample ; DO iSubSurf = 1,nSurfSample
    IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
      SampWall(iSurf)%SurfModelState(5,iSubSurf,jSubSurf) = SampWall(iSurf)%SurfModelState(5,iSubSurf,jSubSurf) &
          + (SampleAdsorptionHeat(iSurf,iSubSurf,jSubSurf) * BoltzmannConst &
          / REAL(SurfDistInfo(iSubSurf,jSubSurf,iSurf)%nSites(3))) &
          * REAL(INT(Adsorption%DensSurfAtoms(iSurf) &
          * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurf),8)) / Species(1)%MacroParticleFactor
    END IF
  END DO ; END DO
END DO

! loop over all surfaces and decide if catalytic boundary of modeltype 3
DO iSurf = 1,SurfMesh%nOutputSides
  IF (SurfaceHasModelNum(iSurf).NE.3) CYCLE
  globSide = SurfMesh%SurfIDToSideID(iSurf)
  PartBoundID = PartBound%MapToPartBC(BC(globSide))
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
DO jSubSurf = 1,nSurfSample ; DO iSubSurf = 1,nSurfSample
  desorbnum(:) = 0
  reactdesorbnum(:) = 0
  nSites(:) = 0
  nSitesRemain(:) = 0
  remainNum(:) = 0
  P_react_back(:) = 0.

  NumSurfReact(:,:) = 0

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
    SpecID = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%Species(Surfpos)
    adsorbates(SpecID) = adsorbates(SpecID) + 1
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
    loc_Surfnu = 0.
#endif

    Surfpos = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(AdsorbID)
    SpecID = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%Species(Surfpos)

    ! move considered particle to end of array, at beginning of the already considered and on surface remaining particles
    SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = &
        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(nSites(Coord)-remainNum(Coord))
    SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(nSites(Coord)-remainNum(Coord)) = Surfpos
    AdsorbID = nSites(Coord)-remainNum(Coord)

    ! calculate heat of adsorption for actual site
    Heat_A = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,SpecID,Surfpos,.FALSE.)
    E_diff = 0.

    ProbDes(:) = 0.
    Coord_ReactP(:) = -1
    Coord_Product(:,:) = -1
    Pos_ReactP(:) = 0
    Pos_Product(:,:) = 0
    ! Choose Random neighbour surface positions appropriate for each possible Reaction and calculate reaction probability

    !-------------------------------------------------------------------------------------------------------------------------------
    ! calculate molecular desorption probability
    !-------------------------------------------------------------------------------------------------------------------------------
    ReactID = 0
    ! estimate vibrational temperature of surface-particle bond
    CharaTemp = Heat_A / 200.
    ! calculate partition function of particles bound on surface
    VarPartitionFuncWall1 = PartitionFuncSurf(SpecID, WallTemp,CharaTemp)
    ! estimate partition function of activated complex
    VarPartitionFuncAct = PartitionFuncActDesorb(SpecID, WallTemp, Adsorption%DensSurfAtoms(iSurf))
    ! transition state theory to estimate pre-exponential factor
    nu_react = ((BoltzmannConst*Walltemp)/PlanckConst) * (VarPartitionFuncAct / VarPartitionFuncWall1)
    E_d = Heat_A
    rate = nu_react * exp(-E_d/WallTemp)
    ProbDes(ReactID) = rate *dt
#if (PP_TimeDiscMethod==42)
    loc_SurfActE(ReactID) = E_d
    loc_Surfnu(ReactID) = nu_react
    CALL AnalyzeSurfRates(1,SpecID,ReactID,E_d,nu_react,ProbDes(ReactID))
#endif
    !-------------------------------------------------------------------------------------------------------------------------------
    ! sort Neighbours to coordinations for search of two random neighbour positions from particle position for dissociation
    !-------------------------------------------------------------------------------------------------------------------------------
    IF (Adsorption%DissNum.GT.0 .OR. Adsorption%RecombNum.GT.0) THEN ! .AND. SpecDSMC(SpecID)%InterID.EQ.2) THEN
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
    IF (SpecDSMC(SpecID)%InterID.EQ.2) THEN
    DO iDissocReact = 1,Adsorption%DissNum
      ReactID = iDissocReact
      Prod_Spec1 = Adsorption%DissocReact(1,iDissocReact,SpecID)
      Prod_Spec2 = Adsorption%DissocReact(2,iDissocReact,SpecID)
      IF ((Prod_Spec1.NE.0) .AND. (Prod_Spec2.NE.0)) THEN !if 2 resulting species, dissociation possible
        Coord_Product(1,ReactID) = Adsorption%Coordination(PartBoundID,Prod_Spec1)
        Coord_Product(2,ReactID) = Adsorption%Coordination(PartBoundID,Prod_Spec2)
        IF ((Coord_Product(1,ReactID).EQ.Coord_Product(2,ReactID))&
              .AND.(nSitesRemain(Coord_Product(1,ReactID)).LT.2)) CYCLE
        IF ((Coord_Product(1,ReactID).NE.Coord_Product(2,ReactID))&
              .AND.((nSitesRemain(Coord_Product(1,ReactID)).LT.1)&
              .OR.(nSitesRemain(Coord_Product(2,ReactID)).LT.1))) CYCLE
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
        Pos_Product(1,ReactID) = Neighpos_j
        Pos_Product(2,ReactID) = Neighpos_k
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
            SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%Species(Surfpos) = SpecID
            IF (Cell_Occupied) CYCLE ! Cycle if space is already occupied
            ! calculation of activation energy
            Heat_Product1 = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Prod_Spec1,Pos_Product(1,ReactID),.TRUE.)
            Heat_Product2 = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Prod_Spec2,Pos_Product(2,ReactID),.TRUE.)
            D_A = Adsorption%EDissBond(iDissocReact,SpecID)
            E_d = CalcDissRecombActEnergy(Heat_A,Heat_Product1,Heat_Product2,D_A,SpecID)
            ! estimate vibrational temperatures of surface-particle bond
            CharaTemp = Heat_A / 200.
            ! calculate partition function of particles bound on surface
            VarPartitionFuncWall1 = PartitionFuncSurf(SpecID, WallTemp,CharaTemp)
            ! estimate partition function of activated complex
            VarPartitionFuncAct = PartitionFuncActDesorb(SpecID,WallTemp,Adsorption%DensSurfAtoms(iSurf))
            ! transition state theory to estimate pre-exponential factor
            nu_react = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct/VarPartitionFuncWall1)
            rate = nu_react * exp(-E_d/(WallTemp))
            ProbDes(ReactID) = rate * dt
#if (PP_TimeDiscMethod==42)
            loc_SurfActE(ReactID) = E_d
            loc_Surfnu(ReactID) = nu_react
            CALL AnalyzeSurfRates(1,SpecID,ReactID,E_d,nu_react,ProbDes(ReactID))
#endif
          ELSE
            ! removed adsorbate temporarily for occupancy-check
            SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%Species(Surfpos) = SpecID
          END IF
        END IF
      END IF
    END DO
    END IF
    !-------------------------------------------------------------------------------------------------------------------------------
    ! calculate probability for associative reactions (recombination)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! product of associative reactions desorb after reaction as in most cases exothermic (LH-reaction)
    DO iRecombReact = 1,Adsorption%RecombNum
      ReactID = iRecombReact + Adsorption%DissNum
      IF (Adsorption%RecombReact(1,iRecombReact,SpecID).LT.1) CYCLE ! no partner for this associative reaction
      Coord_ReactP(ReactID) = Adsorption%Coordination(PartBoundID,Adsorption%RecombReact(1,iRecombReact,SpecID))

      IF (Adsorption%EnableAdsAttraction) THEN
        CALL RANDOM_NUMBER(RanNum)
        NeighID = 1 + INT(REAL(n_Neigh(Coord_ReactP(ReactID)))*RanNum)
        Pos_ReactP(ReactID) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%NeighPos(Surfpos, &
                                  React_NeighbourID(Coord_ReactP(ReactID),NeighID))
      ELSE
        CALL RANDOM_NUMBER(RanNum)
        IF (Coord.EQ.Coord_ReactP(ReactID)) THEN
          NeighID = 1 + INT((nSites(Coord_ReactP(ReactID))-remainNum(Coord_ReactP(ReactID))-1)*RanNum)
        ELSE
          NeighID = 1 + INT((nSites(Coord_ReactP(ReactID))-remainNum(Coord_ReactP(ReactID)))*RanNum)
        END IF
        Pos_ReactP(ReactID) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_ReactP(ReactID))%UsedSiteMap(NeighID)
      END IF

      jSpec = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_ReactP(ReactID))%Species(Pos_ReactP(ReactID))

      IF ( jSpec.EQ.Adsorption%RecombReact(1,iRecombReact,SpecID) .AND. (jSpec.GT.0)) THEN
        Prod_Spec1 = Adsorption%RecombReact(2,iRecombReact,SpecID)
!         Coord_Product(1,ReactID) = Adsorption%Coordination(PartBoundID,Prod_Spec1)
!         CALL RANDOM_NUMBER(RanNum)
!         NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,ReactID))))*RanNum)
!         Pos_Product(1,ReactID) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
!                                           Coord_Product(1,ReactID))%UsedSiteMap(NeighID)
!         Heat_AB = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Prod_Spec1,Pos_Product(1,ReactID),.TRUE.)
        ! calculate heats of adsorption
        Heat_AB = 0. ! immediate desorption
        Heat_B = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,jSpec,Pos_ReactP(ReactID),.FALSE.)
        D_AB = Adsorption%EDissBond(ReactID,SpecID)
        E_d = CalcDissRecombActEnergy(Heat_AB,Heat_A,Heat_B,-D_AB,Prod_Spec1)
        ! check diffusion barrier
        IF ((Coord.NE.3).OR.(Coord_ReactP(ReactID).NE.3)) THEN
          IF ((Coord.NE.3).AND.(Coord_ReactP(ReactID).NE.3)) THEN
            E_diff = (2* (0.1*Heat_A*0.1*Heat_B))/(0.1*Heat_A + 0.1*Heat_B)
          ELSE IF ((Coord.NE.3).AND.(Coord_ReactP(ReactID).EQ.3)) THEN
            E_diff = 0.1*Heat_A
          ELSE IF ((Coord.EQ.3).AND.(Coord_ReactP(ReactID).NE.3)) THEN
            E_diff = 0.1*Heat_B
          ELSE
            E_diff = 0.
          END IF
        END IF
        ! if diffusion barrier is greater then reaction activation, then it is limiting --> activation energy
        IF (E_Diff.GT.E_d) E_d = E_diff
        ! estimate vibrational temperatures of surface-particle bond
        CharaTemp = Heat_A / 200.
        ! calculate partition function of first particle bound on surface
        VarPartitionFuncWall1 = PartitionFuncSurf(SpecID, WallTemp,CharaTemp)
        CharaTemp = Heat_B / 200.
        ! calculate partition function of second particle bound on surface
        VarPartitionFuncWall2 = PartitionFuncSurf(jSpec,WallTemp,CharaTemp)
        ! estimate partition function of activated complex
        VarPartitionFuncAct = PartitionFuncActDesorb(SpecID,WallTemp,Adsorption%DensSurfAtoms(iSurf)) &
                            * PartitionFuncActDesorb(jSpec,WallTemp,Adsorption%DensSurfAtoms(iSurf))
        ! transition state theory to estimate pre-exponential factor
        nu_react = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct)/(VarPartitionFuncWall1*VarPartitionFuncWall2)
        rate = nu_react * exp(-E_d/(WallTemp))
        ProbDes(ReactID) = rate * dt
#if (PP_TimeDiscMethod==42)
        loc_SurfActE(ReactID) = E_d
        loc_Surfnu(ReactID) = nu_react
        CALL AnalyzeSurfRates(1,SpecID,ReactID,E_d,nu_react,ProbDes(ReactID))
#endif
      END IF
    END DO

    !-------------------------------------------------------------------------------------------------------------------------------
    ! Calculate propabilities for exchange reactions
    !-------------------------------------------------------------------------------------------------------------------------------
    ! maybe try later (if exothermic --> Delta_H < 0, product with lower adsorbheat desorbs else both stay adsorbed)
    exch_react_possible = .FALSE.
    DO iReact = 1,Adsorption%nExchReactions
      ReactID = iReact + Adsorption%ReactNum
      ! check if considered particle is one of defined reactants
      IF ( (Adsorption%ChemReactant(1,iReact+Adsorption%nDissocReactions).NE.SpecID) &
           .AND. (Adsorption%ChemReactant(2,iReact+Adsorption%nDissocReactions).NE.SpecID) ) THEN
        ! check if considered particle is one of defined products
        IF ( (Adsorption%ChemProduct(1,iReact+Adsorption%nDissocReactions).NE.SpecID) &
           .AND. (Adsorption%ChemProduct(2,iReact+Adsorption%nDissocReactions).NE.SpecID) ) CYCLE
      END IF
      ! Choose which reaction partner considered particle is
      ! -> defines which species is second reaction partner and if forward or reverse reaction
      ! ----------------------------------------------------------------------------------------------------------------------------
      IF (Adsorption%ChemReactant(1,iReact+Adsorption%nDissocReactions).EQ.SpecID) THEN ! -> first forward
      ! ----------------------------------------------------------------------------------------------------------------------------
        ReactDirForward = .TRUE.
        Coord_ReactP(ReactID) = Adsorption%Coordination(PartBoundID,&
            Adsorption%ChemReactant(2,iReact+Adsorption%nDissocReactions))
        CALL RANDOM_NUMBER(RanNum)
        IF (Coord.EQ.Coord_ReactP(ReactID)) THEN
          NeighID = 1 + INT((nSites(Coord_ReactP(ReactID))-remainNum(Coord_ReactP(ReactID))-1)*RanNum)
        ELSE
          NeighID = 1 + INT((nSites(Coord_ReactP(ReactID))-remainNum(Coord_ReactP(ReactID)))*RanNum)
        END IF
        Pos_ReactP(ReactID) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_ReactP(ReactID))%UsedSiteMap(NeighID)
        jSpec = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_ReactP(ReactID))%Species(Pos_ReactP(ReactID))
        ! continue only if species on chosen site is appropriate partner
        IF ((jSpec.EQ.Adsorption%ChemReactant(2,iReact+Adsorption%nDissocReactions)).AND.(jSpec.GT.0)) THEN
          exch_react_possible = .TRUE.
          Prod_Spec1 = Adsorption%ChemProduct(1,iReact+Adsorption%nDissocReactions)
          Coord_Product(1,ReactID) = Adsorption%Coordination(PartBoundID,Prod_Spec1)
          Prod_Spec2 = Adsorption%ChemProduct(2,iReact+Adsorption%nDissocReactions)
          Coord_Product(2,ReactID) = Adsorption%Coordination(PartBoundID,Prod_Spec2)
          IF ((Coord_Product(1,ReactID).EQ.Coord_Product(2,ReactID))&
               .AND.(nSitesRemain(Coord_Product(1,ReactID)).LT.2)) CYCLE
          IF ((Coord_Product(1,ReactID).NE.Coord_Product(2,ReactID))&
               .AND.((nSitesRemain(Coord_Product(1,ReactID)).LT.1)&
               .OR.(nSitesRemain(Coord_Product(2,ReactID)).LT.1))) CYCLE
          IF (Coord_Product(1,ReactID).EQ.Coord_Product(2,ReactID)) THEN
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,ReactID)))/2.)*RanNum)
            Pos_Product(1,ReactID) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_Product(1,ReactID))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,ReactID)))/2.)*RanNum &
                        + (REAL(nSitesRemain(Coord_Product(2,ReactID)))/2.))
            Pos_Product(2,ReactID) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_Product(2,ReactID))%UsedSiteMap(NeighID)
          ELSE
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,ReactID))))*RanNum)
            Pos_Product(1,ReactID) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_Product(1,ReactID))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,ReactID))))*RanNum)
            Pos_Product(2,ReactID) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_Product(2,ReactID))%UsedSiteMap(NeighID)
          END IF
          ! define dissociation bonds
          D_A = Adsorption%Reactant_DissBond_K(1,iReact+Adsorption%nDissocReactions)
          D_ReactP = Adsorption%Reactant_DissBond_K(2,iReact+Adsorption%nDissocReactions)
          D_Product1 = Adsorption%Product_DissBond_K(1,iReact+Adsorption%nDissocReactions)
          D_Product2 = Adsorption%Product_DissBond_K(2,iReact+Adsorption%nDissocReactions)
        END IF
      ! ----------------------------------------------------------------------------------------------------------------------------
      ELSE IF (Adsorption%ChemReactant(2,iReact+Adsorption%nDissocReactions).EQ.SpecID) THEN ! -> second forward
      ! ----------------------------------------------------------------------------------------------------------------------------
        ReactDirForward = .TRUE.
        Coord_ReactP(ReactID) = Adsorption%Coordination(PartBoundID,&
            Adsorption%ChemReactant(1,iReact+Adsorption%nDissocReactions))
        CALL RANDOM_NUMBER(RanNum)
        IF (Coord.EQ.Coord_ReactP(ReactID)) THEN
          NeighID = 1 + INT((nSites(Coord_ReactP(ReactID))-remainNum(Coord_ReactP(ReactID))-1)*RanNum)
        ELSE
          NeighID = 1 + INT((nSites(Coord_ReactP(ReactID))-remainNum(Coord_ReactP(ReactID)))*RanNum)
        END IF
        Pos_ReactP(ReactID) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                                          Coord_ReactP(ReactID))%UsedSiteMap(NeighID)
        jSpec = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(&
                Coord_ReactP(ReactID))%Species(Pos_ReactP(ReactID))
        ! continue only if species on chosen site is appropriate partner
        IF ((jSpec.EQ.Adsorption%ChemReactant(1,iReact+Adsorption%nDissocReactions)).AND.(jSpec.GT.0)) THEN
          exch_react_possible = .TRUE.
          Prod_Spec1 = Adsorption%ChemProduct(1,iReact+Adsorption%nDissocReactions)
          Coord_Product(1,ReactID) = Adsorption%Coordination(PartBoundID,Prod_Spec1)
          Prod_Spec2 = Adsorption%ChemProduct(2,iReact+Adsorption%nDissocReactions)
          Coord_Product(2,ReactID) = Adsorption%Coordination(PartBoundID,Prod_Spec2)
          IF ((Coord_Product(1,ReactID).EQ.Coord_Product(2,ReactID)).AND.(nSitesRemain(Coord_Product(1,ReactID)).LT.2)) CYCLE
          IF ((Coord_Product(1,ReactID).NE.Coord_Product(2,ReactID))&
               .AND.((nSitesRemain(Coord_Product(1,ReactID)).LT.1)&
               .OR.(nSitesRemain(Coord_Product(2,ReactID)).LT.1))) CYCLE
          IF (Coord_Product(1,ReactID).EQ.Coord_Product(2,ReactID)) THEN
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,ReactID)))/2.)*RanNum)
            Pos_Product(1,ReactID) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_Product(1,ReactID))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,ReactID)))/2.)*RanNum &
                        + (REAL(nSitesRemain(Coord_Product(2,ReactID)))/2.))
            Pos_Product(2,ReactID) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_Product(2,ReactID))%UsedSiteMap(NeighID)
          ELSE
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,ReactID))))*RanNum)
            Pos_Product(1,ReactID) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_Product(1,ReactID))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,ReactID))))*RanNum)
            Pos_Product(2,ReactID) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_Product(2,ReactID))%UsedSiteMap(NeighID)
          END IF
          ! define dissociation bonds
          D_A = Adsorption%Reactant_DissBond_K(2,iReact+Adsorption%nDissocReactions)
          D_ReactP = Adsorption%Reactant_DissBond_K(1,iReact+Adsorption%nDissocReactions)
          D_Product1 = Adsorption%Product_DissBond_K(1,iReact+Adsorption%nDissocReactions)
          D_Product2 = Adsorption%Product_DissBond_K(2,iReact+Adsorption%nDissocReactions)
        END IF
      ! ----------------------------------------------------------------------------------------------------------------------------
      ELSE IF (Adsorption%ChemProduct(1,iReact+Adsorption%nDissocReactions).EQ.SpecID) THEN ! -> first reverse
      ! ----------------------------------------------------------------------------------------------------------------------------
        ReactDirForward = .FALSE.
        Coord_ReactP(ReactID) = Adsorption%Coordination(PartBoundID,&
            Adsorption%ChemProduct(2,iReact+Adsorption%nDissocReactions))
        CALL RANDOM_NUMBER(RanNum)
        IF (Coord.EQ.Coord_ReactP(ReactID)) THEN
          NeighID = 1 + INT((nSites(Coord_ReactP(ReactID))-remainNum(Coord_ReactP(ReactID))-1)*RanNum)
        ELSE
          NeighID = 1 + INT((nSites(Coord_ReactP(ReactID))-remainNum(Coord_ReactP(ReactID)))*RanNum)
        END IF
        Pos_ReactP(ReactID) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_ReactP(ReactID))%UsedSiteMap(NeighID)
        jSpec = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_ReactP(ReactID))%Species(Pos_ReactP(ReactID))
        ! continue only if species on chosen site is appropriate partner
        IF ((jSpec.EQ.Adsorption%ChemProduct(2,iReact+Adsorption%nDissocReactions)).AND.(jSpec.GT.0)) THEN
          exch_react_possible = .TRUE.
          Prod_Spec1 = Adsorption%ChemReactant(1,iReact+Adsorption%nDissocReactions)
          Coord_Product(1,ReactID) = Adsorption%Coordination(PartBoundID,Prod_Spec1)
          Prod_Spec2 = Adsorption%ChemReactant(2,iReact+Adsorption%nDissocReactions)
          Coord_Product(2,ReactID) = Adsorption%Coordination(PartBoundID,Prod_Spec2)
          IF ((Coord_Product(1,ReactID).EQ.Coord_Product(2,ReactID)).AND.(nSitesRemain(Coord_Product(1,ReactID)).LT.2)) CYCLE
          IF ((Coord_Product(1,ReactID).NE.Coord_Product(2,ReactID))&
               .AND.((nSitesRemain(Coord_Product(1,ReactID)).LT.1)&
               .OR.(nSitesRemain(Coord_Product(2,ReactID)).LT.1))) CYCLE
          IF (Coord_Product(1,ReactID).EQ.Coord_Product(2,ReactID)) THEN
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,ReactID)))/2.)*RanNum)
            Pos_Product(1,ReactID) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_Product(1,ReactID))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,ReactID)))/2.)*RanNum &
                        + (REAL(nSitesRemain(Coord_Product(2,ReactID)))/2.))
            Pos_Product(2,ReactID) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_Product(2,ReactID))%UsedSiteMap(NeighID)
          ELSE
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,ReactID))))*RanNum)
            Pos_Product(1,ReactID) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_Product(1,ReactID))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,ReactID))))*RanNum)
            Pos_Product(2,ReactID) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_Product(2,ReactID))%UsedSiteMap(NeighID)
          END IF
          ! define dissociation bonds
          D_A = Adsorption%Product_DissBond_K(1,iReact+Adsorption%nDissocReactions)
          D_ReactP = Adsorption%Product_DissBond_K(2,iReact+Adsorption%nDissocReactions)
          D_Product1 = Adsorption%Reactant_DissBond_K(1,iReact+Adsorption%nDissocReactions)
          D_Product2 = Adsorption%Reactant_DissBond_K(2,iReact+Adsorption%nDissocReactions)
        END IF
      ! ----------------------------------------------------------------------------------------------------------------------------
      ELSE IF (Adsorption%ChemProduct(2,iReact+Adsorption%nDissocReactions).EQ.SpecID) THEN ! -> second reverse
      ! ----------------------------------------------------------------------------------------------------------------------------
        ReactDirForward = .FALSE.
        Coord_ReactP(ReactID) = Adsorption%Coordination(PartBoundID,&
            Adsorption%ChemProduct(1,iReact+Adsorption%nDissocReactions))
        CALL RANDOM_NUMBER(RanNum)
        IF (Coord.EQ.Coord_ReactP(ReactID)) THEN
          NeighID = 1 + INT((nSites(Coord_ReactP(ReactID))-remainNum(Coord_ReactP(ReactID))-1)*RanNum)
        ELSE
          NeighID = 1 + INT((nSites(Coord_ReactP(ReactID))-remainNum(Coord_ReactP(ReactID)))*RanNum)
        END IF
        Pos_ReactP(ReactID) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_ReactP(ReactID))%UsedSiteMap(NeighID)
        jSpec = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_ReactP(ReactID))%Species(Pos_ReactP(ReactID))
        ! continue only if species on chosen site is appropriate partner
        IF ((jSpec.EQ.Adsorption%ChemProduct(1,iReact+Adsorption%nDissocReactions)).AND.(jSpec.GT.0)) THEN
          exch_react_possible = .TRUE.
          Prod_Spec1 = Adsorption%ChemReactant(1,iReact+Adsorption%nDissocReactions)
          Coord_Product(1,ReactID) = Adsorption%Coordination(PartBoundID,Prod_Spec1)
          Prod_Spec2 = Adsorption%ChemReactant(2,iReact+Adsorption%nDissocReactions)
          Coord_Product(2,ReactID) = Adsorption%Coordination(PartBoundID,Prod_Spec2)
          IF ((Coord_Product(1,ReactID).EQ.Coord_Product(2,ReactID)).AND.(nSitesRemain(Coord_Product(1,ReactID)).LT.2)) CYCLE
          IF ((Coord_Product(1,ReactID).NE.Coord_Product(2,ReactID))&
               .AND.((nSitesRemain(Coord_Product(1,ReactID)).LT.1)&
               .OR.(nSitesRemain(Coord_Product(2,ReactID)).LT.1))) CYCLE
          IF (Coord_Product(1,ReactID).EQ.Coord_Product(2,ReactID)) THEN
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,ReactID)))/2.)*RanNum)
            Pos_Product(1,ReactID) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_Product(1,ReactID))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,ReactID)))/2.)*RanNum &
                        + (REAL(nSitesRemain(Coord_Product(2,ReactID)))/2.))
            Pos_Product(2,ReactID) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_Product(2,ReactID))%UsedSiteMap(NeighID)
          ELSE
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,ReactID))))*RanNum)
            Pos_Product(1,ReactID) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_Product(1,ReactID))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,ReactID))))*RanNum)
            Pos_Product(2,ReactID) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_Product(2,ReactID))%UsedSiteMap(NeighID)
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
        Heat_ReactP = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,jSpec,Pos_ReactP(ReactID),.FALSE.)
        Heat_Product1 = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Prod_Spec1,Pos_Product(1,ReactID),.TRUE.)
        Heat_Product2 = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Prod_Spec2,Pos_Product(2,ReactID),.TRUE.)
        ! calculate activation energy
        E_d = Calc_E_Act(Heat_Product1,Heat_Product2,Heat_A,Heat_ReactP,D_Product1,D_Product2,D_A,D_ReactP)
        ! check diffusion barrier
        IF ((Coord.NE.3).OR.(Coord_ReactP(ReactID).NE.3)) THEN
          IF ((Coord.NE.3).AND.(Coord_ReactP(ReactID).NE.3)) THEN
            E_diff = (2* (0.1*Heat_A*0.1*Heat_ReactP))/(0.1*Heat_A + 0.1*Heat_ReactP)
          ELSE IF ((Coord.NE.3).AND.(Coord_ReactP(ReactID).EQ.3)) THEN
            E_diff = 0.1*Heat_A
          ELSE IF ((Coord.EQ.3).AND.(Coord_ReactP(ReactID).NE.3)) THEN
            E_diff = 0.1*Heat_ReactP
          ELSE
            E_diff = 0.
          END IF
        END IF
        ! if diffusion barrier is greater then reaction activation, then it is limiting --> activation energy
        IF (E_Diff.GT.E_d) E_d = E_diff
        ! estimate vibrational temperatures of surface-particle bond
        CharaTemp = Heat_A / 200.
        ! calculate partition function of first particle bound on surface
        VarPartitionFuncWall1 = PartitionFuncSurf(SpecID, WallTemp,CharaTemp)
        CharaTemp = Heat_ReactP / 200.
        ! calculate partition function of second particle bound on surface
        VarPartitionFuncWall2 = PartitionFuncSurf(jSpec, WallTemp,CharaTemp)
        ! estimate partition function of activated complex
        VarPartitionFuncAct = PartitionFuncActDesorb(SpecID,WallTemp,Adsorption%DensSurfAtoms(iSurf)) &
                            * PartitionFuncActDesorb(jSpec,WallTemp,Adsorption%DensSurfAtoms(iSurf))
        ! transition state theory to estimate pre-exponential factor
        nu_react = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct/(VarPartitionFuncWall1*VarPartitionFuncWall2))
        rate = nu_react * exp(-E_d/(WallTemp))
        ProbDes(ReactID) = rate * dt
        IF (rate*dt.GT.1) THEN
          IF (ReactDirForward) THEN
            P_react_forward(iReact) = P_react_forward(iReact) + 1
          ELSE
            P_react_back(iReact) = P_react_back(iReact) + 1
          END IF
        ELSE
          IF (ReactDirForward) THEN
            P_react_forward(iReact) = P_react_forward(iReact) + ProbDes(ReactID)
          ELSE
            P_react_back(iReact) = P_react_back(iReact) + ProbDes(ReactID)
          END IF
        END IF
      END IF
    END DO

    SDEALLOCATE(NeighbourID)
    SDEALLOCATE(React_NeighbourID)

    ! initialize sum of all probabilities
    sum_probabilities = 0.
    DO iReact = 0,Adsorption%ReactNum+Adsorption%nExchReactions
      sum_probabilities = sum_probabilities + ProbDes(iReact)
    END DO
    SurfModel%Info(SpecID)%MeanProbDes = SurfModel%Info(SpecID)%MeanProbDes + sum_probabilities
    SurfModel%Info(SpecID)%MeanProbDesCount = SurfModel%Info(SpecID)%MeanProbDesCount + 1
    !-------------------------------------------------------------------------------------------------------------------------------
    ! choose which surface reaction takes place
    !-------------------------------------------------------------------------------------------------------------------------------
    CALL RANDOM_NUMBER(RanNum)
    IF (sum_probabilities .GT. RanNum) THEN
      DO iReact = 0,Adsorption%ReactNum+Adsorption%nExchReactions
        CALL RANDOM_NUMBER(RanNum)
        IF ((ProbDes(iReact)/sum_probabilities).GT.RanNum) THEN
          IF (iReact.EQ.0) THEN
            surf_react_case = 1 ! desorption
          ELSE IF ( iReact.GT.0 .AND. iReact.LE.Adsorption%DissNum .AND. iReact.LE.Adsorption%ReactNum ) THEN
            surf_react_case = 2 ! dissociation
          ELSE IF ( iReact.GT.0 .AND. iReact.GT.Adsorption%DissNum .AND. iReact.LE.Adsorption%ReactNum ) THEN
            surf_react_case = 3 ! recombination
          ELSE IF ( iReact.GT.Adsorption%ReactNum ) THEN
            surf_react_case = 4 ! exchange reaction
          END IF
          EXIT
        ELSE
          sum_probabilities = sum_probabilities - ProbDes(iReact)
        END IF
      END DO
#if (PP_TimeDiscMethod==42)
      CALL AnalyzeSurfRates(2,SpecID,iReact,loc_SurfActE(iReact),loc_Surfnu(iReact),ProbDes(iReact))
#endif
      IF (DSMC%ReservoirSurfaceRate) surf_react_case = 0 !only probabilities and analyze are calculated without actual desorption
      !-----------------------------------------------------------------------------------------------------------------------------
      SELECT CASE(surf_react_case)
      !-----------------------------------------------------------------------------------------------------------------------------
      CASE(1) !(desorption)
      !-----------------------------------------------------------------------------------------------------------------------------
        desorbnum(SpecID) = desorbnum(SpecID) + 1
        CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,Coord,Surfpos,SpecID,.TRUE.,relaxation=.TRUE.)
        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = &
            SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1)
        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1) = Surfpos
        nSitesRemain(Coord) = nSitesRemain(Coord) + 1
        nSitesRemain(4) = nSitesRemain(4) + 1
      !-----------------------------------------------------------------------------------------------------------------------------
      CASE(2) ! dissociation
      !-----------------------------------------------------------------------------------------------------------------------------
        iDissocReact = iReact
        ! update number of reactions
        reactdesorbnum(Adsorption%DissocReact(1,iDissocReact,SpecID)) = &
                reactdesorbnum(Adsorption%DissocReact(1,iDissocReact,SpecID)) + 1
        reactdesorbnum(Adsorption%DissocReact(2,iDissocReact,SpecID)) = &
                reactdesorbnum(Adsorption%DissocReact(2,iDissocReact,SpecID)) + 1
        reactdesorbnum(SpecID) = reactdesorbnum(SpecID) - 1

        IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd))&
            .OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
        !  ! calculate Enthalpy of desorption and sample
        !  Heat_A = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,SpecID,-1,.FALSE.)
        !  Heat_Product1 = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Adsorption%DissocReact(1,iDissocReact,SpecID),-1,.FALSE.)
        !  Heat_Product2 = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Adsorption%DissocReact(2,iDissocReact,SpecID),-1,.FALSE.)
          Heat_A = 0.0
          Heat_Product1 = 0.0
          Heat_Product2 = 0.0
          D_A = Adsorption%EDissBond(iReact,SpecID)
          AdsorptionEnthalpy = ((( Heat_A -Heat_Product1 -Heat_Product2 ) + D_A) * BoltzmannConst &
                          / REAL(SurfDistInfo(iSubSurf,jSubSurf,iSurf)%nSites(3))) &
                          * REAL(INT(Adsorption%DensSurfAtoms(iSurf) &
                          * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurf),8)) / Species(SpecID)%MacroParticleFactor
          !----  Sampling of energies
          SampWall(iSurf)%SurfModelState(2,iSubSurf,jSubSurf) = SampWall(iSurf)%SurfModelState(2,iSubSurf,jSubSurf) &
                                                                + AdsorptionEnthalpy * Species(SpecID)%MacroParticleFactor
          NumSurfReact(SpecID,iReact) = NumSurfReact(SpecID,iReact) + 1
        END IF
#if (PP_TimeDiscMethod==42)
        D_A = Adsorption%EDissBond(iReact,SpecID)
        AdsorptionEnthalpy = ((( Heat_A -Heat_Product1 -Heat_Product2 ) + D_A) * BoltzmannConst &
                        / REAL(SurfDistInfo(iSubSurf,jSubSurf,iSurf)%nSites(3))) &
                        * REAL(INT(Adsorption%DensSurfAtoms(iSurf) &
                        * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurf),8)) / Species(SpecID)%MacroParticleFactor
        SurfModel%ProperInfo(SpecID)%HeatFlux(2) = SurfModel%ProperInfo(SpecID)%HeatFlux(2) &
            + AdsorptionEnthalpy/BoltzmannConst
        iSampleReact = iReact + 1
        SurfModel%ProperInfo(SpecID)%HeatFluxDesCount(iSampleReact) = &
            SurfModel%ProperInfo(SpecID)%HeatFluxDesCount(iSampleReact) &
                        + (1. / REAL(SurfDistInfo(iSubSurf,jSubSurf,iSurf)%nSites(3))) &
                        * REAL(INT(Adsorption%DensSurfAtoms(iSurf) * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurf),8))
#endif
        ! remove adsorbate and update map
        CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,Coord,SurfPos,SpecID,.TRUE.,relaxation=.TRUE.)
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
          jSpec = Adsorption%DissocReact(1,iDissocReact,SpecID)
          CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,jCoord,Pos_Product(1,iReact),jSpec,.FALSE.,relaxation=.TRUE.)
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
            jSpec = Adsorption%DissocReact(2,iDissocReact,SpecID)
            CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,jCoord,Pos_Product(2,iReact),jSpec,.FALSE.,relaxation=.TRUE.)
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
      CASE(3) ! recombination
      !-----------------------------------------------------------------------------------------------------------------------------
        iRecombReact = iReact - Adsorption%DissNum
        ! update number of desorptions
        Prod_Spec1 = Adsorption%RecombReact(1,iRecombReact,SpecID) ! Partner species
        Prod_Spec2 = Adsorption%RecombReact(2,iRecombReact,SpecID) ! Results species
        desorbnum(Prod_Spec2) = desorbnum(Prod_Spec2) + 1
        reactdesorbnum(Prod_Spec2) = reactdesorbnum(Prod_Spec2) + 1
        reactdesorbnum(Prod_Spec1) = reactdesorbnum(Prod_Spec1) - 1
        reactdesorbnum(SpecID) = reactdesorbnum(SpecID) - 1
        IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd))&
            .OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
        !  ! calculate Enthalpy of desorption and sample
        !  Heat_A = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,SpecID,-1,.FALSE.)
        !  Heat_B = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Prod_Spec1,-1,.FALSE.)
        !  Heat_AB = 0.
          Heat_A = 0.
          Heat_B = 0.
          Heat_AB = 0.
          D_AB = Adsorption%EDissBond(iReact,SpecID)
          AdsorptionEnthalpy = (-(( Heat_AB -Heat_A -Heat_B ) + D_AB) * BoltzmannConst &
                          / REAL(SurfDistInfo(iSubSurf,jSubSurf,iSurf)%nSites(3))) &
                          * REAL(INT(Adsorption%DensSurfAtoms(iSurf) &
                          * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurf),8)) / Species(SpecID)%MacroParticleFactor
          !----  Sampling of energies
          SampWall(iSurf)%SurfModelState(1,iSubSurf,jSubSurf) = SampWall(iSurf)%SurfModelState(1,iSubSurf,jSubSurf) &
                                                                + AdsorptionEnthalpy * Species(SpecID)%MacroParticleFactor
          NumSurfReact(SpecID,iReact) = NumSurfReact(SpecID,iReact) + 1
        END IF
#if (PP_TimeDiscMethod==42)
        D_AB = Adsorption%EDissBond(iReact,SpecID)
        AdsorptionEnthalpy = (-(( Heat_AB -Heat_A -Heat_B ) + D_AB) * BoltzmannConst &
                        / REAL(SurfDistInfo(iSubSurf,jSubSurf,iSurf)%nSites(3))) &
                        * REAL(INT(Adsorption%DensSurfAtoms(iSurf) &
                        * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurf),8)) / Species(SpecID)%MacroParticleFactor
        SurfModel%ProperInfo(SpecID)%HeatFlux(2) = SurfModel%ProperInfo(SpecID)%HeatFlux(2) &
            + AdsorptionEnthalpy/BoltzmannConst
        iSampleReact = iReact + 1
        SurfModel%ProperInfo(SpecID)%HeatFluxDesCount(iSampleReact) = &
            SurfModel%ProperInfo(SpecID)%HeatFluxDesCount(iSampleReact) &
                        + (1. / REAL(SurfDistInfo(iSubSurf,jSubSurf,iSurf)%nSites(3))) &
                        * REAL(INT(Adsorption%DensSurfAtoms(iSurf) * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurf),8))
#endif
        ! remove adsorbate and update map
        CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,Coord,Surfpos,SpecID,.TRUE.,relaxation=.TRUE.)
        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = &
            SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1)
        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1) = Surfpos
        nSitesRemain(Coord) = nSitesRemain(Coord) + 1
        nSitesRemain(4) = nSitesRemain(4) + 1
        ! remove adsorbate reaction partner and update map
        DO PartnerID = nSitesRemain(Coord_ReactP(iReact))+1,nSites(Coord_ReactP(iReact))
          IF (SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord_ReactP(iReact))%UsedSiteMap(PartnerID) &
              .EQ.Pos_ReactP(iReact)) THEN
            CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,Coord_ReactP(iReact),Pos_ReactP(iReact),Prod_Spec1,.TRUE.,relaxation=.TRUE.)
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
      CASE(4) ! exchange reactions
      !-----------------------------------------------------------------------------------------------------------------------------
        ExchNum = iReact - Adsorption%ReactNum + Adsorption%nDissocReactions
        ! choose if forward or reverse reaction and update number of reactions
        IF ( (Adsorption%ChemReactant(1,ExchNum).EQ.SpecID) .OR. (Adsorption%ChemReactant(2,ExchNum).EQ.SpecID) ) THEN
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
        !  ! calculate Enthalpy of reaction and sample
        !  Heat_A = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Adsorption%ChemReactant(1,ExchNum),-1,.FALSE.)
        !  Heat_ReactP = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Adsorption%ChemReactant(2,ExchNum),-1,.FALSE.)
        !  Heat_Product1 = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Adsorption%ChemProduct(1,ExchNum),-1,.FALSE.)
        !  Heat_Product2 = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,Adsorption%ChemProduct(2,ExchNum),-1,.FALSE.)
          Heat_A = 0.0
          Heat_ReactP = 0.0
          Heat_Product1 = 0.0
          Heat_Product2 = 0.0
          D_A = Adsorption%Reactant_DissBond_K(1,ExchNum)
          D_ReactP = Adsorption%Reactant_DissBond_K(2,ExchNum)
          D_Product1 = Adsorption%Product_DissBond_K(1,ExchNum)
          D_Product2 = Adsorption%Product_DissBond_K(2,ExchNum)
          AdsorptionEnthalpy = (-(( Heat_A +Heat_ReactP -Heat_Product1 -Heat_Product2 ) &
                          + (D_A +D_ReactP -D_Product1 -D_Product2)) * BoltzmannConst &
                          / REAL(SurfDistInfo(iSubSurf,jSubSurf,iSurf)%nSites(3))) &
                          * REAL(INT(Adsorption%DensSurfAtoms(iSurf) &
                          * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurf),8)) / Species(SpecID)%MacroParticleFactor
          IF (ReactDirForward) AdsorptionEnthalpy = -AdsorptionEnthalpy
          !----  Sampling of energies
          SampWall(iSurf)%SurfModelState(1,iSubSurf,jSubSurf) = SampWall(iSurf)%SurfModelState(1,iSubSurf,jSubSurf) &
                                                                + AdsorptionEnthalpy * Species(SpecID)%MacroParticleFactor
        END IF
        ! remove first adsorbate and update map
        CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,Coord,Surfpos,SpecID,.TRUE.,relaxation=.TRUE.)
        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = &
            SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1)
        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1) = Surfpos
        nSitesRemain(Coord) = nSitesRemain(Coord) + 1
        nSitesRemain(4) = nSitesRemain(4) + 1
        ! remove second adsorbate (reaction partner) and update map
        IF (ReactDirForward) THEN
          IF (Adsorption%ChemReactant(1,ExchNum).EQ.SpecID) jSpec = Adsorption%ChemReactant(2,ExchNum)
          IF (Adsorption%ChemReactant(2,ExchNum).EQ.SpecID) jSpec = Adsorption%ChemReactant(1,ExchNum)
        ELSE
          IF (Adsorption%ChemProduct(1,ExchNum).EQ.SpecID) jSpec = Adsorption%ChemProduct(2,ExchNum)
          IF (Adsorption%ChemProduct(2,ExchNum).EQ.SpecID) jSpec = Adsorption%ChemProduct(1,ExchNum)
        END IF
        jCoord = Coord_ReactP(iReact)
        CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,jCoord,Pos_ReactP(iReact),jSpec,.TRUE.,relaxation=.TRUE.)
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
          CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,jCoord,Pos_Product(1,iReact),jSpec,.FALSE.,relaxation=.TRUE.)
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
            CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,jCoord,Pos_Product(2,iReact),jSpec,.FALSE.,relaxation=.TRUE.)
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
    SurfModel%SumDesorbPart(iSubSurf,jSubSurf,iSurf,iSpec) = &
                        INT(SurfDistInfo(iSubSurf,jSubSurf,iSurf)%desorbnum_tmp(iSpec))
    ! add desorbed particles to be inserted at emission
    SurfModel%SumEvapPart(iSubSurf,jSubSurf,iSurf,iSpec) = SurfModel%SumEvapPart(iSubSurf,jSubSurf,iSurf,iSpec) &
                        + SurfModel%SumDesorbPart(iSubSurf,jSubSurf,iSurf,iSpec)
    ! Adjust tracking desorbing simulation particles
    SurfDistInfo(iSubSurf,jSubSurf,iSurf)%desorbnum_tmp(iSpec) = &
                        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%desorbnum_tmp(iSpec) &
                        - REAL(SurfModel%SumDesorbPart(iSubSurf,jSubSurf,iSurf,iSpec))
    ! calculate number of reacted particles for each species
    SurfDistInfo(iSubSurf,jSubSurf,iSurf)%reactnum_tmp(iSpec) = &
                        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%reactnum_tmp(iSpec) &
                        + ((REAL(reactdesorbnum(iSpec)) / REAL(numSites)) &
                        * REAL(INT(Adsorption%DensSurfAtoms(iSurf) &
                        * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurf),8)) / Species(iSpec)%MacroParticleFactor)
    ! calculate number of reacting simulation particles on surface (round to integer)
    SurfModel%SumReactPart(iSubSurf,jSubSurf,iSurf,iSpec) = &
                        INT(SurfDistInfo(iSubSurf,jSubSurf,iSurf)%reactnum_tmp(iSpec))
    ! Adjust tracking reacting simulation particles
    SurfDistInfo(iSubSurf,jSubSurf,iSurf)%reactnum_tmp(iSpec) = &
                        SurfDistInfo(iSubSurf,jSubSurf,iSurf)%reactnum_tmp(iSpec) &
                        - REAL(SurfModel%SumReactPart(iSubSurf,jSubSurf,iSurf,iSpec))
    ! Sample vibrational energies
    ! due to reaction a part of energies can be transformed into other vibrational groundstates and the change is sampled
    ! the energies of the emitted particles are sampled in surfflux routine (particle_emission.f90) because calculated there
    IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
      DO iReact = 1,Adsorption%ReactNum
        IF (NumSurfReact(iSpec,iReact).GT.0) THEN
          SampWall(iSurf)%SurfModelReactCount(iReact+Adsorption%ReactNum,iSpec,iSubSurf,jSubSurf) = &
              SampWall(iSurf)%SurfModelReactCount(iReact+Adsorption%ReactNum,iSpec,iSubSurf,jSubSurf) &
              + ((REAL(NumSurfReact(iSpec,iReact)) / REAL(numSites)) * REAL(INT(Adsorption%DensSurfAtoms(iSurf) &
              * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurf),8)) / Species(iSpec)%MacroParticleFactor)
        END IF
      END DO
    END IF
#if (PP_TimeDiscMethod==42)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! analyze desorption data in case TPD is enabled
    !-------------------------------------------------------------------------------------------------------------------------------
    IF (DSMC%ReservoirSimu) THEN
      SurfModel%Info(iSpec)%NumOfDes = SurfModel%Info(iSpec)%NumOfDes + SurfModel%SumDesorbPart(iSubSurf,jSubSurf,iSurf,iSpec)
    END IF
#endif
  END DO ! nSpecies (analyze)
END DO ; END DO ! nSurfSample
END DO ! SurfMesh%nOutputSides

! sample energy of surfaces after desorption treatment
DO iSurf = 1,SurfMesh%nOutputSides
  IF (SurfaceHasModelNum(iSurf).NE.3) CYCLE
  DO jSubSurf = 1,nSurfSample ; DO iSubSurf = 1,nSurfSample
    IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
      SampWall(iSurf)%SurfModelState(5,iSubSurf,jSubSurf) = SampWall(iSurf)%SurfModelState(5,iSubSurf,jSubSurf) &
          - (SampleAdsorptionHeat(iSurf,iSubSurf,jSubSurf) * BoltzmannConst &
          / REAL(SurfDistInfo(iSubSurf,jSubSurf,iSurf)%nSites(3))) &
          * REAL(INT(Adsorption%DensSurfAtoms(iSurf) &
          * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurf),8)) / Species(1)%MacroParticleFactor
    END IF
  END DO ; END DO
END DO

DEALLOCATE(desorbnum,adsorbnum,nSites,nSitesRemain,remainNum,adsorbates)
DEALLOCATE(ProbDes,P_react_forward,P_react_back)
DEALLOCATE(Coord_ReactP,Pos_ReactP)

! 6. communicate surface state to halo sides of neighbours
#if USE_MPI
! communicate distribution to halo-sides of neighbour procs
CALL ExchangeSurfDistInfo()
#endif

END SUBROUTINE SMCR_PartDesorb


SUBROUTINE SMCR_Diffusion()
!===================================================================================================================================
!> Calculation of diffusion on reconstructed surface with assumption of Quasi Chemical Approximation (QCA)
!>   diffusion into equilibrium (Quasi Chemical Approximation - QCA) is performed for particles on surface
!===================================================================================================================================
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_Particle_Vars          ,ONLY: Species, WriteMacroSurfaceValues
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption, SurfDistInfo
USE MOD_SurfaceModel_Tools     ,ONLY: UpdateSurfPos, Calc_Adsorb_Heat, SpaceOccupied, SampleAdsorptionHeat, SurfaceHasModelNum
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, PartBound, SampWall
USE MOD_TimeDisc_Vars          ,ONLY: dt, TEnd, time
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_SurfaceModel_PartFunc  ,ONLY: PartitionFuncActDesorb, PartitionFuncSurf
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
IF (.NOT.SurfMesh%SurfOnProc) RETURN
IF (Adsorption%NoDiffusion) RETURN

! diffusion into equilibrium distribution
DO iSurf=1,SurfMesh%nOutputSides
  IF (SurfaceHasModelNum(iSurf).NE.3) CYCLE
  DO jSubSurf=1,nSurfSample ; DO iSubSurf=1,nSurfSample

    IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
      SampWall(iSurf)%SurfModelState(5,iSubSurf,jSubSurf) = SampWall(iSurf)%SurfModelState(5,iSubSurf,jSubSurf) &
          + (SampleAdsorptionHeat(iSurf,iSubSurf,jSubSurf) * BoltzmannConst &
          / REAL(SurfDistInfo(iSubSurf,jSubSurf,iSurf)%nSites(3))) &
          * REAL(INT(Adsorption%DensSurfAtoms(iSurf) &
          * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurf),8)) / Species(1)%MacroParticleFactor
    END IF

    DO Coord = 1,3
      nSites = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%nSites(Coord)
      nSitesRemain = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%SitesRemain(Coord)

      ALLOCATE ( free_Neigh_pos(1:SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%nNeighbours))
      DO AdsorbID = nSitesRemain+1,nSites,1
        Surfpos = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%UsedSiteMap(AdsorbID)
        SpecID = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%Species(Surfpos)

        ! choose Random vacant neighbour position
        n_equal_site_Neigh = 0
        free_Neigh_pos(:) = 0

        ! find free Neighbour positions of the same site-coordination
        DO i = 1,SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%nNeighbours
          IF (SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%NeighSite(Surfpos,i) .EQ. Coord) THEN
            IF ( (SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%Species( &
                SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%NeighPos(Surfpos,i)).EQ.0) ) THEN
              ! check for occupation with nearest Neighbours of position
              IF (SpaceOccupied(iSurf,iSubSurf,jSubSurf,Coord &
                  ,SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%NeighPos(Surfpos,i))) CYCLE
              n_equal_site_Neigh = n_equal_site_Neigh + 1
              free_Neigh_pos(n_equal_site_Neigh) = SurfDistInfo(iSubSurf,jSubSurf,iSurf)%AdsMap(Coord)%NeighPos(Surfpos,i)
            END IF
          END IF
        END DO

        ! calculate heat of adsorption for actual site
        Heat_i = Calc_Adsorb_Heat(iSubSurf,jSubSurf,iSurf,SpecID,Surfpos,.FALSE.)
        ! update surfatom bond order and species map
        CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,Coord,Surfpos,SpecID,.TRUE.,relaxation=.TRUE.)

        ! choose Neighbour position with highest heat of adsorption
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
          globSide = SurfMesh%SurfIDToSideID(iSurf)
          WallTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC(globSide)))
          Prob_diff = exp(-(Heat_i - Heat_j)/WallTemp) / (1+exp(-(Heat_i - Heat_j)/Walltemp)) ! QCA
          CALL RANDOM_NUMBER(RanNum)
          IF (dt.LT.1e-1) Prob_diff = Prob_diff*dt
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
        CALL UpdateSurfPos(iSurf,iSubSurf,jSubSurf,Coord,newpos,SpecID,.FALSE.,relaxation=.TRUE.)
      END DO
      DEALLOCATE(free_Neigh_pos)
    END DO
    IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
      SampWall(iSurf)%SurfModelState(5,iSubSurf,jSubSurf) = SampWall(iSurf)%SurfModelState(5,iSubSurf,jSubSurf) &
          - (SampleAdsorptionHeat(iSurf,iSubSurf,jSubSurf) * BoltzmannConst &
          / REAL(SurfDistInfo(iSubSurf,jSubSurf,iSurf)%nSites(3))) &
          * REAL(INT(Adsorption%DensSurfAtoms(iSurf) &
          * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurf),8)) / Species(1)%MacroParticleFactor
    END IF
  END DO ; END DO !iSubSurf = 1,nSurfSample; jSubSurf = 1,nSurfSample
END DO !iSurf = 1,SurfMesh%nOutputSides

END SUBROUTINE SMCR_Diffusion

END MODULE MOD_SMCR
