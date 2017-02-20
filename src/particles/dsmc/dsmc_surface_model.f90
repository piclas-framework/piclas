#include "boltzplatz.h"

MODULE MOD_DSMC_SurfModel_Tools
!===================================================================================================================================
! Initialization of DSMC
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
PUBLIC                       :: DSMC_Update_Wall_Vars
! PUBLIC                       :: Particle_Wall_Adsorb
! PUBLIC                       :: CalcPartAdsorb
PUBLIC                       :: Calc_PartNum_Wall_Desorb
PUBLIC                       :: CalcBackgndPartAdsorb
PUBLIC                       :: CalcBackgndPartDesorb
PUBLIC                       :: CalcAdsorbProb
PUBLIC                       :: CalcDesorbProb
PUBLIC                       :: CalcDiffusion
#ifdef MPI
PUBLIC                       :: ExchangeSurfDistInfo
PUBLIC                       :: ExchangeSurfDistSize
#endif /*MPI*/
!===================================================================================================================================

CONTAINS

SUBROUTINE DSMC_Update_Wall_Vars()   
!===================================================================================================================================
! Update and sample DSMC-values for adsorption, desorption and reactions on surfaces
!===================================================================================================================================
  USE MOD_PARTICLE_Vars,          ONLY : WriteMacroValues, nSpecies
  USE MOD_PARTICLE_Vars,          ONLY : KeepWallParticles, Species
  USE MOD_DSMC_Vars,              ONLY : DSMC, Adsorption, SurfDistInfo
  USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh, SampWall
!===================================================================================================================================
  IMPLICIT NONE
!===================================================================================================================================
! Local variable declaration 
  INTEGER                          :: iSpec, iSurfSide, p, q, new_adsorbates, numSites
  REAL                             :: maxPart
!===================================================================================================================================

  IF (DSMC%WallModel.GT.0) THEN    
    IF (.NOT.KeepWallParticles) THEN
#ifdef MPI
      ! communicate number of particles that were adsorbed on halo-sides of neighbour procs
      CALL ExchangeAdsorbNum()
#endif
      ! adjust coverages of all species on surfaces
      DO iSpec = 1,nSpecies
#if (PP_TimeDiscMethod==42)
        IF (DSMC%ReservoirRateStatistic) Adsorption%AdsorpInfo(iSpec)%WallSpecNumCount = 0
#endif
        DO iSurfSide = 1,SurfMesh%nSides
          DO q = 1,nSurfSample
            DO p = 1,nSurfSample
#if (PP_TimeDiscMethod==42)
              ! write number of adsorbates into info array
              IF (DSMC%ReservoirRateStatistic) THEN
                maxPart = REAL(INT(Adsorption%DensSurfAtoms(iSurfSide) * SurfMesh%SurfaceArea(p,q,iSurfSide),8))
                Adsorption%AdsorpInfo(iSpec)%WallSpecNumCount = Adsorption%AdsorpInfo(iSpec)%WallSpecNumCount &
                                                              + INT( Adsorption%Coverage(p,q,iSurfSide,iSpec) &
                                                              * maxPart/Species(iSpec)%MacroParticleFactor)
              END IF
#endif
              IF (DSMC%WallModel.EQ.1) THEN
                maxPart = Adsorption%DensSurfAtoms(iSurfSide) * SurfMesh%SurfaceArea(p,q,iSurfSide)
                Adsorption%Coverage(p,q,iSurfSide,iSpec) = Adsorption%Coverage(p,q,iSurfSide,iSpec) &
                    + ( Adsorption%SumAdsorbPart(p,q,iSurfSide,iSpec) &
                    - (Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) - Adsorption%SumReactPart(p,q,iSurfSide,iSpec)) ) &
                    * Species(iSpec)%MacroParticleFactor / maxPart
              ELSE IF (DSMC%WallModel.GT.1) THEN
                maxPart = REAL(INT(Adsorption%DensSurfAtoms(iSurfSide) * SurfMesh%SurfaceArea(p,q,iSurfSide),8))
                Adsorption%Coverage(p,q,iSurfSide,iSpec) = Adsorption%Coverage(p,q,iSurfSide,iSpec) &
                    + ( Adsorption%SumAdsorbPart(p,q,iSurfSide,iSpec) &
                    - (Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) - Adsorption%SumReactPart(p,q,iSurfSide,iSpec)) ) &
                    * Species(iSpec)%MacroParticleFactor / maxPart
                ! adjust number of background mapping adsorbates if SumAdsorbPart > 0
                IF (Adsorption%SumAdsorbPart(p,q,iSurfSide,iSpec).GT.0) THEN
                  ! calculate number of adsorbed particles on background for each species
                  numSites = SurfDistInfo(p,q,iSurfSide)%nSites(3) !number of simulated surface atoms
                  SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec) = SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec) &
                        + (REAL(Adsorption%SumAdsorbPart(p,q,iSurfSide,iSpec)) * Species(iSpec)%MacroParticleFactor &
                        / maxPart) * REAL(numSites)
                  ! convert to integer adsorbates
                  new_adsorbates = INT(SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec))
                  IF (new_adsorbates.GT.0) THEN
                    ! Adjust tracking of adsorbing background particles
                    SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec) = SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec) &
                                                                      - new_adsorbates
                    CALL AdjustBackgndAdsNum(p,q,iSurfSide,new_adsorbates,iSpec)
                  END IF
                END IF
              END IF
              ! sample adsorption coverage
              IF (WriteMacroValues) THEN
                SampWall(iSurfSide)%Adsorption(1+iSpec,p,q) = SampWall(iSurfSide)%Adsorption(1+iSpec,p,q) &
                                                            + Adsorption%Coverage(p,q,iSurfSide,iSpec)
              END IF
            END DO
          END DO
        END DO
      END DO
      Adsorption%SumDesorbPart(:,:,:,:) = 0
      Adsorption%SumAdsorbPart(:,:,:,:) = 0
      Adsorption%SumReactPart(:,:,:,:) = 0
    END IF
    
#ifdef MPI
    IF (DSMC%WallModel.EQ.3) THEN
      ! communicate distribution to halo-sides of neighbour procs
      CALL ExchangeSurfDistInfo()
    END IF
#endif

    IF (DSMC%WallModel.EQ.2) THEN
      CALL CalcDistNumChange()
      CALL CalcDiffusion()
      CALL CalcSurfDistInteraction()
    END IF
    
    IF (DSMC%WallModel.EQ.1 .AND. DSMC%WallModel.EQ.2) THEN
      CALL CalcAdsorbProb()
      IF (KeepWallParticles) CALL CalcDesorbprob()
    END IF
  END IF

END SUBROUTINE DSMC_Update_Wall_Vars  


SUBROUTINE Calc_PartNum_Wall_Desorb()
!===================================================================================================================================
! calculation of desorbing particle number when particles deleted at adsorption and inserted at desorption
!===================================================================================================================================
USE MOD_Particle_Vars,          ONLY : nSpecies, Species
USE MOD_DSMC_Vars,              ONLY : Adsorption, DSMC
USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh
!===================================================================================================================================
  IMPLICIT NONE
!===================================================================================================================================
! Local variable declaration
!    INTEGER                          :: i, PartAds2
   INTEGER                          :: iSurfSide, iSpec, p, q, NPois, WallPartNum, iProbSigma
   REAL                             :: PartAds, PartDes, RanNum, Tpois
!===================================================================================================================================
#if (PP_TimeDiscMethod==42)
  Adsorption%AdsorpInfo(:)%MeanProbDes = 0.
  Adsorption%AdsorpInfo(:)%NumOfDes = 0
#endif
  CALL CalcDesorbProb()
  
  DO iSpec = 1,nSpecies
    DO iSurfSide = 1,SurfMesh%nSides
      DO q = 1,nSurfSample
        DO p = 1,nSurfSample
          IF (Adsorption%Coverage(p,q,iSurfSide,iSpec).GT.0) THEN
            WallPartNum = INT(Adsorption%Coverage(p,q,iSurfSide,iSpec) * Adsorption%DensSurfAtoms(iSurfSide) &
                      * SurfMesh%SurfaceArea(p,q,iSurfSide) &
                      / Species(iSpec)%MacroParticleFactor)
            IF (WallPartNum .GT. 0) THEN
            PartAds = Adsorption%Coverage(p,q,iSurfSide,iSpec) * Adsorption%DensSurfAtoms(iSurfSide) &
                        * SurfMesh%SurfaceArea(p,q,iSurfSide) &
                        / Species(iSpec)%MacroParticleFactor
            
            IF ((DSMC%WallModel.EQ.1) .OR. (DSMC%WallModel.EQ.3)) THEN
              PartDes = PartAds * Adsorption%ProbDes(p,q,iSurfSide,iSpec)
            ELSE IF (DSMC%WallModel.EQ.2) THEN 
              PartDes = 0.
              DO iProbSigma = (1+(iSpec-1)*36),(36*iSpec)
                PartDes = PartDes + PartAds * Adsorption%ProbSigma(p,q,iSurfSide,iSpec,iProbSigma) &
                                            * Adsorption%ProbSigDes(p,q,iSurfSide,iSpec,iProbSigma) 
              END DO
            END IF
              IF (PartDes.GT.PartAds) THEN
                Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) = Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) + INT(PartAds)
              ELSE
                CALL RANDOM_NUMBER(RanNum)            
                IF (EXP(-PartDes).LE.TINY(PartDes) &
#if (PP_TimeDiscMethod==42)
                .OR. Adsorption%TPD &
#endif
                ) THEN
                  Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) = Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec)&
                  +INT(PartDes + RanNum)
                ELSE !poisson-sampling instead of random rounding (reduces numeric non-equlibrium effects [Tysanner and Garcia 2004]
                  Npois=0
                  Tpois=1.0
                  DO
                    Tpois=RanNum*Tpois
                    IF (Tpois.LT.TINY(Tpois)) THEN
                      Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) = Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec)&
                      +INT(PartDes + RanNum)
                      EXIT
                    END IF
                    IF (Tpois.GT.EXP(-PartDes)) THEN
                      Npois=Npois+1
                      CALL RANDOM_NUMBER(RanNum)
                    ELSE
                      Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) = Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec)+Npois
                      EXIT
                    END IF
                  END DO
                END IF
              END IF !PartDes.GT.WallPartNum
              IF (Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec).GT.WallPartNum  &
              .OR.Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec).LT.0) THEN
                Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) = WallPartNum
              END IF
            ELSE !not PartAds.GT.0
              Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) = 0
            END IF !PartAds.GT.0
#if (PP_TimeDiscMethod==42)
            Adsorption%AdsorpInfo(iSpec)%NumOfDes = Adsorption%AdsorpInfo(iSpec)%NumOfDes &
                                                  + Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec)
#endif
          END IF
        END DO  
      END DO
    END DO
  END DO
END SUBROUTINE Calc_PartNum_Wall_Desorb

SUBROUTINE CalcDiffusion()
!===================================================================================================================================
! Model for diffusion calculation
!===================================================================================================================================
  USE MOD_Mesh_Vars,              ONLY : BC
  USE MOD_DSMC_Vars,              ONLY : Adsorption, SurfDistInfo
  USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh, PartBound
!===================================================================================================================================
   IMPLICIT NONE
!===================================================================================================================================
! Local variable declaration
   INTEGER                          :: SurfSideID, iSpec, globSide, subsurfeta, subsurfxi
   INTEGER                          :: Coord, nSites, nSitesRemain, i, j, AdsorbID
   REAL                             :: WallTemp, Prob_diff, RanNum
   REAL                             :: Heat_i, Heat_j, Heat_temp
   INTEGER                          :: n_equal_site_Neigh, Indx, Indy, Surfpos, newpos
   INTEGER , ALLOCATABLE            :: free_Neigh_pos(:)
!===================================================================================================================================

DO SurfSideID = 1,SurfMesh%nSides
DO subsurfeta = 1,nSurfSample
DO subsurfxi = 1,nSurfSample

  DO Coord = 1,3
    nSites = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Coord)
    nSitesRemain = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord)
    
    ALLOCATE ( free_Neigh_pos(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nNeighbours))
    DO AdsorbID = nSitesRemain+1,nSites,1
      Surfpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID)
      iSpec = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(Surfpos)
      n_equal_site_Neigh = 0
      free_Neigh_pos(:) = 0
      
      ! find free Neighbour positions of the same site-coordination
      DO i = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nNeighbours
        IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighSite(Surfpos,i) .EQ. Coord) THEN
          IF ( (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species( &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,i)).EQ.0) ) THEN
            n_equal_site_Neigh = n_equal_site_Neigh + 1
            free_Neigh_pos(n_equal_site_Neigh) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,i)
          END IF
        END IF
      END DO
      
      ! calculate heat of adsorption for actual site and reduce bond order of bondatoms
      Heat_i = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,iSpec,Surfpos,.FALSE.)
      DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
        Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(Surfpos,j)
        Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(Surfpos,j)
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) = &
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) - 1
      END DO
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(Surfpos) = 0
      
      ! choose Neighbour position with highest heat of adsorption if adsorbate would move there
      Heat_j = 0.
      DO i = 1,n_equal_site_Neigh
        Heat_temp = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,iSpec,free_Neigh_pos(i),.TRUE.)
        IF (Heat_temp .GT. Heat_j) THEN
          Heat_j = Heat_temp
          newpos = free_Neigh_pos(i)
        END IF
      END DO
      
      ! only try to diffuse particle if unoccupied sites available
      IF (n_equal_site_Neigh .GE. 1) THEN
        globSide = Adsorption%SurfSideToGlobSideMap(SurfSideID)
        WallTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC(globSide)))
        Prob_diff = exp(-(Heat_i - Heat_j)/WallTemp) / (1+exp(-(Heat_i - Heat_j)/Walltemp))
        CALL RANDOM_NUMBER(RanNum)
        IF (Prob_diff.GT.RanNum) THEN
        ! move particle to new position and update map
          DO i = 1,nSitesRemain
            IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(i).EQ.newpos) THEN
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(i) = Surfpos
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = newpos
            END IF
          END DO
        ELSE 
          newpos = Surfpos
        END IF
      ELSE
        newpos = Surfpos        
      END IF ! end if (n_equal_site_Neigh >= 1)
      
      ! update surfatom bond order and species map
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(newpos) = iSpec
      DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
        Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(newpos,j)
        Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(newpos,j) 
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) = &
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) + 1
      END DO
    END DO
    DEALLOCATE(free_Neigh_pos)
    
  END DO
  
END DO
END DO
END DO

END SUBROUTINE CalcDiffusion


SUBROUTINE CalcBackgndPartAdsorb(subsurfxi,subsurfeta,SurfSideID,PartID,Norm_Ec,adsorption_case,outSpec,AdsorptionEnthalpie)
!===================================================================================================================================
! Particle Adsorption probability calculation for wallmodel 3
!===================================================================================================================================
  USE MOD_Globals_Vars,           ONLY : PlanckConst
  USE MOD_Particle_Vars,          ONLY : PartSpecies, nSpecies, Species, BoltzmannConst, WriteMacroValues
  USE MOD_Mesh_Vars,              ONLY : BC
  USE MOD_DSMC_Vars,              ONLY : Adsorption, DSMC, SurfDistInfo, SpecDSMC, PartStateIntEn
  USE MOD_Particle_Boundary_Vars, ONLY : PartBound, SampWall
  USE MOD_DSMC_Analyze,           ONLY : CalcTVib, CalcTVibPoly
  USE MOD_TimeDisc_Vars,          ONLY : TEnd, time
#if (PP_TimeDiscMethod==42)  
  USE MOD_TimeDisc_Vars,          ONLY : iter, dt
#endif
!===================================================================================================================================
  IMPLICIT NONE
!===================================================================================================================================
! argument list declaration
  INTEGER,INTENT(IN)               :: subsurfxi,subsurfeta,SurfSideID,PartID
  REAL   ,INTENT(IN)               :: Norm_Ec
  INTEGER,INTENT(OUT)              :: adsorption_case
  INTEGER,INTENT(OUT)              :: outSpec(2)
  REAL   ,INTENT(OUT)              :: AdsorptionEnthalpie
! LOCAL VARIABLES
  INTEGER                          :: iSpec, globSide
  INTEGER                          :: Coord, Coord2, i, AdsorbID, IDRearrange
  INTEGER                          :: nSites, nSitesRemain
  REAL                             :: WallTemp, RanNum
  REAL                             :: Prob_ads
  REAL , ALLOCATABLE               :: P_Eley_Rideal(:), Prob_diss(:)
  INTEGER                          :: Surfpos, ReactNum
  REAL, PARAMETER                  :: Pi=3.14159265358979323846_8
  INTEGER                          :: jSpec, kSpec, jCoord, kCoord
  REAL                             :: sum_probabilities
  INTEGER , ALLOCATABLE            :: NeighbourID(:,:), NeighSpec(:)
  INTEGER                          :: SiteSpec, nNeigh_trap, Neighpos_j, Neighpos_k, chosen_Neigh_j, chosen_Neigh_k
  INTEGER                          :: n_empty_Neigh(3), n_Neigh(3), adsorbates(nSpecies)
  REAL                             :: E_a, c_f, EZeroPoint_Educt, phi_1, phi_2, Xi_Total, Xi_vib
  REAL                             :: Heat_A, Heat_B, Heat_AB, D_AB, D_A, D_B
  
  REAL                             :: vel_norm, vel_coll, potential_pot, a_const, mu, surfmass, trapping_prob
!===================================================================================================================================
! special TPD (temperature programmed desorption) surface temperature adjustment part
  globSide = Adsorption%SurfSideToGlobSideMap(SurfSideID)
#if (PP_TimeDiscMethod==42)
  IF (Adsorption%TPD) THEN
    WallTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC(globSide))) + (Adsorption%TPD_beta * iter * dt)
    Adsorption%TPD_Temp = Walltemp
  ELSE
    WallTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC(globSide)))
  END IF
#else
  WallTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC(globSide)))
#endif
  ! calculate number of adsorbates for each species (already on surface)
  adsorbates(:) = 0
  DO Coord = 1,3
  nSites = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Coord)
  nSitesRemain = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord)
  DO AdsorbID = 1,nSites-nSitesRemain
    Surfpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain+AdsorbID)
    iSpec = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(Surfpos)
    adsorbates(iSpec) = adsorbates(iSpec) + 1
  END DO
  END DO
  ! initialize variables
  Prob_ads = 0.
  ALLOCATE( P_Eley_Rideal(1:Adsorption%ReactNum),&
            Prob_diss(1:Adsorption%ReactNum))
  P_Eley_Rideal(:) = 0.
  Prob_diss(:) = 0.
  
  ! Choose Random surface site with species coordination
  CALL RANDOM_NUMBER(RanNum)
  iSpec = PartSpecies(PartID)
  Coord = Adsorption%Coordination(iSpec)
  AdsorbID = 1 + INT(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Coord)*RanNum)
  Surfpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID)
  SiteSpec = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(Surfpos)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! calculate trapping probability (using hard cube collision with surface atom or adsorbate)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! if site is empty nearest neighbour site can be occupied and this influences the collision cube mass
  IF (SiteSpec.EQ.0) THEN
    ALLOCATE(NeighSpec(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nNeighbours))
    NeighSpec = 0
    nNeigh_trap = 0
    DO i = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nNeighbours
      IF (.NOT.SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%IsNearestNeigh(Surfpos,i)) CYCLE
    DO Coord2 = 1,3
      IF (Coord2.EQ.SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighSite(Surfpos,i)) THEN
        IF ( (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%Species( & 
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,i)) & 
              .NE.0) ) THEN
              nNeigh_trap = nNeigh_trap + 1
              NeighSpec(nNeigh_trap) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%Species( & 
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,i))
        END IF
      END IF
    END DO
    END DO
  END IF
  potential_pot = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,iSpec,Surfpos,.TRUE.)*BoltzmannConst
  vel_norm = - ( 2*(Norm_Ec+potential_pot)/Species(iSpec)%MassIC)**(0.5)
  a_const = (Species(iSpec)%MassIC/(2*BoltzmannConst*WallTemp))**(0.5)
  IF (SiteSpec.EQ.0) THEN
    IF (nNeigh_trap.EQ.0) THEN
      surfmass = Adsorption%SurfMassIC
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
  IF ((DSMC%CalcSurfaceVal.AND.(Time.ge.(1-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroValues)) THEN
    SampWall(SurfSideID)%Accomodation(iSpec,subsurfxi,subsurfeta) = SampWall(SurfSideID)%Accomodation(iSpec,subsurfxi,subsurfeta) &
                                                                  + trapping_prob
  END IF
  ! if no trapping return and perform elastic reflection
  CALL RANDOM_NUMBER(RanNum)
  IF (RanNum.GT.trapping_prob) THEN
    outSpec(1) = iSpec
    outSpec(2) = 0
    AdsorptionEnthalpie = 0.
    adsorption_case = -1
    RETURN
  END IF
  !---------------------------------------------------------------------------------------------------------------------------------
  ! calculate probability for molecular adsorption
  !---------------------------------------------------------------------------------------------------------------------------------
  IF (SiteSpec.EQ.0) THEN
    ! calculation of molecular adsorption probability with TCE
    EZeroPoint_Educt = 0.
    E_a = 0.
    c_f = Adsorption%DensSurfAtoms(SurfSideID) / ( (BoltzmannConst / (2*Pi*Species(iSpec)%MassIC))**0.5 )
    ! Testing if the adsorption particle is an atom or molecule, if molecule: is it polyatomic?
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
      END IF
    ELSE
      Xi_vib = 0.0
    END IF
    Xi_Total = Xi_vib + 2 + 3 !Xi_rot + 3
    phi_1 = Adsorption%Ads_Powerfactor(iSpec) - 3./2. + Xi_Total
    phi_2 = 1 - Xi_Total
    Prob_ads = Calc_Beta_Adsorb(iSpec,Xi_Total,c_f) * (Norm_Ec - E_a)**phi_1 * Norm_Ec ** phi_2
  ELSE
    Prob_ads = 0.
  END IF
  !---------------------------------------------------------------------------------------------------------------------------------
  ! sort Neighbours to coordinations for search of two random neighbour positions from impact position for dissociative adsorption
  !---------------------------------------------------------------------------------------------------------------------------------
  n_Neigh(:) = 0
  n_empty_Neigh(:) = 0
  ALLOCATE(NeighbourID(1:3,1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nNeighbours))

  DO i = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nNeighbours
  DO Coord2 = 1,3
    IF (Coord2.EQ.SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighSite(Surfpos,i)) THEN
      IF ( (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%Species( & 
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,i)) & 
            .EQ.0) ) THEN
        n_empty_Neigh(Coord2) = n_empty_Neigh(Coord2) + 1
      END IF
      n_Neigh(Coord2) = n_Neigh(Coord2) + 1
      NeighbourID(Coord2,n_Neigh(Coord2)) = i
    END IF
  END DO
  END DO
  !---------------------------------------------------------------------------------------------------------------------------------
  ! calculate probability for dissociative adsorption
  !---------------------------------------------------------------------------------------------------------------------------------
  IF (SpecDSMC(iSpec)%InterID.EQ.2) THEN
  DO ReactNum = 1,(Adsorption%DissNum)
    jSpec = Adsorption%DissocReact(1,ReactNum,iSpec)
    kSpec = Adsorption%DissocReact(2,ReactNum,iSpec)
    IF ((jSpec.NE.0) .AND. (kSpec.NE.0)) THEN !if 2 resulting species, dissociation possible
      jCoord = Adsorption%Coordination(jSpec)
      kCoord = Adsorption%Coordination(kSpec)
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
            Neighpos_k = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,&
                                                                                              NeighbourID(kCoord,chosen_Neigh_k))
            IDRearrange = NeighbourID(jCoord,chosen_Neigh_k)
            NeighbourID(jCoord,chosen_Neigh_k) = NeighbourID(kCoord,n_Neigh(kCoord))
            NeighbourID(jCoord,n_Neigh(jCoord)) = IDRearrange
          END IF
        ELSE
          Neighpos_k = Surfpos
          IF (n_empty_Neigh(jCoord).GT.0) THEN
            ! assign availiable neighbour position for j
            CALL RANDOM_NUMBER(RanNum)
            chosen_Neigh_j = 1 + INT(REAL(n_Neigh(jCoord))*RanNum)
            Neighpos_j = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,&
                                                                                              NeighbourID(jCoord,chosen_Neigh_j))
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
          Neighpos_j = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,&
                                                                                            NeighbourID(jCoord,chosen_Neigh_j))
          IDRearrange = NeighbourID(jCoord,chosen_Neigh_j)
          NeighbourID(jCoord,chosen_Neigh_j) = NeighbourID(jCoord,n_Neigh(jCoord))
          NeighbourID(jCoord,n_Neigh(jCoord)) = IDRearrange
          CALL RANDOM_NUMBER(RanNum)
          chosen_Neigh_k = 1 + INT(REAL(n_Neigh(kCoord)-1)*RanNum)
          Neighpos_k = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,&
                                                                                            NeighbourID(kCoord,chosen_Neigh_k))
          IDRearrange = NeighbourID(jCoord,chosen_Neigh_j)
          NeighbourID(jCoord,chosen_Neigh_k) = NeighbourID(kCoord,n_Neigh(kCoord)-1)
          NeighbourID(jCoord,n_Neigh(jCoord)-1) = IDRearrange
        ELSE IF ( (jCoord.NE.kCoord) .AND. (n_empty_Neigh(jCoord).GT.0) .AND. (n_empty_Neigh(kCoord).GT.0) ) THEN
          ! assign availiable neighbour positions
          CALL RANDOM_NUMBER(RanNum)
          chosen_Neigh_j = 1 + INT(REAL(n_Neigh(jCoord))*RanNum)
          Neighpos_j = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,&
                                                                                            NeighbourID(jCoord,chosen_Neigh_j))
          IDRearrange = NeighbourID(jCoord,chosen_Neigh_j)
          NeighbourID(jCoord,chosen_Neigh_j) = NeighbourID(jCoord,n_Neigh(jCoord))
          NeighbourID(jCoord,n_Neigh(jCoord)) = IDRearrange
          CALL RANDOM_NUMBER(RanNum)
          chosen_Neigh_k = 1 + INT(REAL(n_Neigh(kCoord))*RanNum)
          Neighpos_k = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,&
                                                                                            NeighbourID(kCoord,chosen_Neigh_k))
          IDRearrange = NeighbourID(jCoord,chosen_Neigh_j)
          NeighbourID(jCoord,chosen_Neigh_k) = NeighbourID(kCoord,n_Neigh(kCoord))
          NeighbourID(jCoord,n_Neigh(jCoord)) = IDRearrange
        END IF
      END IF
      IF ( (Neighpos_j.GT.0) .AND. (Neighpos_k.GT.0) ) THEN !both neighbour positions assigned
      !both assigned neighbour positions empty
      IF ( (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%Species(Neighpos_j).EQ.0) .AND. &
          (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(kCoord)%Species(Neighpos_k).EQ.0) ) THEN
        ! calculation of activation energy
        Heat_A = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,jSpec,Neighpos_j,.TRUE.)
        Heat_B = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,kSpec,Neighpos_k,.TRUE.)
        Heat_AB = 0. ! direct dissociative adsorption
        D_AB = Adsorption%EDissBond(ReactNum,iSpec)
        ! find dissociation bond energies of products if at least one of them is diatomic
!         IF ((SpecDSMC(jSpec)%InterID.EQ.2) .OR. (SpecDSMC(kSpec)%InterID.EQ.2)) THEN
!           DO iReactNum2 = 1,(Adsorption%ReactNum-Adsorption%DissNum)
!             IF ((Adsorption%AssocReact(2,iReactNum2,jSpec).EQ.iSpec) .OR. (Adsorption%AssocReact(2,iReactNum2,kSpec).EQ.iSpec) THEN
!             IF ((Adsorption%AssocReact(1,iReactNum2,jSpec).EQ.kSpec) .OR. (Adsorption%AssocReact(1,iReactNum2,kSpec).EQ.jSpec) THEN
!               D_A = Adsorption%EDissBond(Adsorption%ReactNum+iReactNum2,jSpec)
!               D_B = Adsorption%EDissBond(Adsorption%ReactNum+iReactNum2,kSpec)
!             END IF
!             END IF
!           END DO
!         ELSE
          D_A = 0.
          D_B = 0.
!         END IF
        E_a = Calc_E_Act(Heat_A,Heat_B,Heat_AB,0.,D_A,D_B,D_AB,0.) * BoltzmannConst
        ! calculation of dissociative adsorption probability with TCE
        EZeroPoint_Educt = 0.
        c_f = Adsorption%DensSurfAtoms(SurfSideID) / ( (BoltzmannConst / (2*Pi*Species(iSpec)%MassIC))**0.5 )
        ! Testing if the adsorption particle is an atom or molecule, if molecule: is it polyatomic?
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
          END IF
        ELSE
          Xi_vib = 0.0
        END IF
        IF ((Norm_Ec).GT.E_a) THEN
          Xi_Total = Xi_vib + 2 + 3 !Xi_rot + 3
          phi_1 = Adsorption%Diss_Powerfactor(ReactNum,iSpec) - 3./2. + Xi_Total
          phi_2 = 1 - Xi_Total
          Prob_diss(ReactNum) = Calc_Beta_Diss(ReactNum,iSpec,Xi_Total,c_f) * (Norm_Ec - E_a)**phi_1 * Norm_Ec ** phi_2
        END IF
      END IF
      END IF
    END IF
  END DO
  END IF
  !---------------------------------------------------------------------------------------------------------------------------------
  ! calculate probability for Eley-Rideal reaction (not ready yet)
  !---------------------------------------------------------------------------------------------------------------------------------
  DO ReactNum = 1,(Adsorption%ReactNum-Adsorption%DissNum)
    ! reaction partner
    jSpec = Adsorption%AssocReact(1,ReactNum,iSpec)
    Neighpos_j = 0
    IF (jSpec.EQ.0) CYCLE
    ! reaction results
    kSpec = Adsorption%AssocReact(2,ReactNum,iSpec)
    jCoord = Adsorption%Coordination(jSpec)
    IF (jCoord.EQ.Coord) THEN
      Neighpos_j = Surfpos
    ELSE
      IF ( n_Neigh(jCoord)-n_empty_Neigh(jCoord).GT.0 ) THEN
        CALL RANDOM_NUMBER(RanNum)
        chosen_Neigh_j = 1 + INT(n_Neigh(jCoord)*RanNum)
        Neighpos_j = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighPos(Surfpos,NeighbourID(jCoord,chosen_Neigh_j))
      END IF
    END IF
    IF ( (Neighpos_j.GT.0) ) THEN
    IF ( (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%Species(Neighpos_j).EQ.jSpec) ) THEN
      ! calculation of activation energy
      ! calculation of reaction probability with TCE
    END IF
    END IF
  END DO
  !---------------------------------------------------------------------------------------------------------------------------------
  
  sum_probabilities = Prob_ads
  DO ReactNum = 1,Adsorption%ReactNum
    sum_probabilities = sum_probabilities + Prob_diss(ReactNum) + P_Eley_Rideal(ReactNum)
  END DO
  ! choose adsorption case
  adsorption_case = 0
  AdsorptionEnthalpie = 0.
  CALL RANDOM_NUMBER(RanNum)
  IF (sum_probabilities .GT. RanNum) THEN
    ! chose surface reaction case (0=inelastic scattering, 1=adsorption, 2=reaction (dissociation), 3=reaction (Eley-Rideal))
    DO ReactNum = 1,(Adsorption%ReactNum)
      CALL RANDOM_NUMBER(RanNum)
      IF ((P_Eley_Rideal(ReactNum)/sum_probabilities).GT.RanNum) THEN
        ! if ER-reaction set output parameters
        adsorption_case = 3
        outSpec(1) = Adsorption%AssocReact(1,ReactNum,iSpec)
        outSpec(2) = Adsorption%AssocReact(2,ReactNum,iSpec)
        ! calculate adsorption Enthalpie
        Heat_A = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,outSpec(1),-1,.TRUE.)
        Heat_B = 0.
        Heat_AB = 0.
        D_AB = Adsorption%EDissBond(ReactNum,iSpec)
        D_A = 0.
        D_B = 0.
        AdsorptionEnthalpie = (( Heat_AB -Heat_A -Heat_B ) + ( D_AB -D_A -D_B )) * BoltzmannConst
        EXIT
      END IF
      sum_probabilities = sum_probabilities - P_Eley_Rideal(ReactNum)
      CALL RANDOM_NUMBER(RanNum)
      IF ((Prob_diss(ReactNum)/sum_probabilities).GT.RanNum) THEN
        ! if dissocciative adsorption set output parameters
        adsorption_case = 2
        outSpec(1) = Adsorption%DissocReact(1,ReactNum,iSpec)
        outSpec(2) = Adsorption%DissocReact(2,ReactNum,iSpec)
        ! calculate adsorption Enthalpie
        Heat_A = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,outSpec(1),-1,.TRUE.)
        Heat_B = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,outSpec(2),-1,.TRUE.)
        Heat_AB = 0.
        D_AB = Adsorption%EDissBond(ReactNum,iSpec)
        D_A = 0.
        D_B = 0.
        AdsorptionEnthalpie = (( Heat_AB -Heat_A -Heat_B ) + ( D_AB -D_A -D_B )) * BoltzmannConst
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
        AdsorptionEnthalpie = - Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,outSpec(1),-1,.TRUE.) * BoltzmannConst
      END IF
    END IF
  END IF
#if (PP_TimeDiscMethod==42)
  IF (.NOT.DSMC%ReservoirRateStatistic) THEN
    Adsorption%AdsorpInfo(iSpec)%MeanProbAds = Adsorption%AdsorpInfo(iSpec)%MeanProbDes + Prob_ads
  END IF
#endif

  DEALLOCATE(P_Eley_Rideal,Prob_diss,NeighbourID)

END SUBROUTINE CalcBackgndPartAdsorb

SUBROUTINE CalcBackgndPartDesorb()
!===================================================================================================================================
! Routine for calculation of number of desorbing particles from background surface distribution
!===================================================================================================================================
  USE MOD_Globals_Vars,           ONLY : PlanckConst
  USE MOD_Particle_Vars,          ONLY : WriteMacroValues, nSpecies, Species, BoltzmannConst
  USE MOD_Mesh_Vars,              ONLY : BC
  USE MOD_DSMC_Vars,              ONLY : Adsorption, DSMC, SurfDistInfo, SpecDSMC, PolyatomMolDSMC
  USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh, PartBound, SampWall
  USE MOD_TimeDisc_Vars,          ONLY : dt
  USE MOD_TimeDisc_Vars,          ONLY : TEnd, time
#if (PP_TimeDiscMethod==42)  
  USE MOD_TimeDisc_Vars,          ONLY : iter
#endif
!===================================================================================================================================
   IMPLICIT NONE
!===================================================================================================================================
! Local variable declaration
   INTEGER                          :: SurfSideID, iSpec, globSide, subsurfeta, subsurfxi
   INTEGER                          :: Coord, Coord3, i, j, AdsorbID, numSites
   REAL                             :: WallTemp, RanNum
   REAL                             :: Heat_A, Heat_B, nu_des, rate, P_actual_des
   REAL , ALLOCATABLE               :: P_des(:)
#if (PP_TimeDiscMethod==42)
   REAL , ALLOCATABLE               :: Energy(:)
#endif
   INTEGER                          :: Indx, Indy, Surfpos
   REAL                             :: VarPartitionFuncAct, VarPartitionFuncWall
   INTEGER                          :: trace, traceNum, ReactNum_run
   INTEGER , ALLOCATABLE            :: desorbnum(:), reactdesorbnum(:)
   INTEGER , ALLOCATABLE            :: adsorbnum(:)
   INTEGER , ALLOCATABLE            :: nSites(:), nSitesRemain(:), remainNum(:), adsorbates(:)
!---------- reaction variables
   INTEGER                          :: react_Neigh, n_Neigh(3), n_empty_Neigh(3), jSpec, iReact, PartnerID
   INTEGER                          :: surf_react_case, NeighID
   REAL                             :: E_d, nu_react
   REAL                             :: Heat_AB, D_AB, sum_probabilities
   REAL                             :: VarPartitionFuncWall1, VarPartitionFuncWall2, VarPartitionFuncAct1, VarPartitionFuncAct2
   INTEGER , ALLOCATABLE            :: Pos_ReactP(:), Coord_ReactP(:), Pos_Product(:,:), Coord_Product(:,:)
   INTEGER , ALLOCATABLE            :: React_NeighbourID(:,:), NeighbourID(:)
   REAL , ALLOCATABLE               :: P_react(:), P_actual_react(:)
   REAL                             :: AdsorptionEnthalpie
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
  LOGICAL                           :: ReactDirForward
  INTEGER                           :: DissocNum, ExchNum, jCoord
!===================================================================================================================================
ALLOCATE (desorbnum(1:nSpecies),&
          reactdesorbnum(1:nSpecies),&
          adsorbnum(1:4),&
          nSites(1:4),&
          nSitesRemain(1:4),&
          remainNum(1:4),&
          adsorbates(1:nSpecies),&
          P_des(1:nSpecies),&
          P_react(1:Adsorption%ReactNum+Adsorption%nDisPropReactions),&
          P_actual_react(1:Adsorption%ReactNum+1+Adsorption%nDisPropReactions),&
          Coord_ReactP(1:Adsorption%ReactNum+Adsorption%nDisPropReactions),&
          Coord_Product(1:2,1:Adsorption%ReactNum+Adsorption%nDisPropReactions),&
          Pos_ReactP(1:Adsorption%ReactNum+1+Adsorption%nDisPropReactions),&
          Pos_Product(1:2,1:Adsorption%ReactNum+Adsorption%nDisPropReactions))
#if (PP_TimeDiscMethod==42)
IF (Adsorption%TPD) THEN
  ALLOCATE(Energy(1:nSpecies))
  Adsorption%AdsorpInfo(:)%MeanEAds = 0.
END IF
! Adsorption%AdsorpInfo(:)%MeanProbDes = 0.
Adsorption%AdsorpInfo(:)%NumOfDes = 0
Adsorption%AdsorpInfo(:)%NumOfAds = 0
#endif

DO SurfSideID = 1,SurfMesh%nSides
  globSide = Adsorption%SurfSideToGlobSideMap(SurfSideID)
  
! special TPD (temperature programmed desorption) temperature adjustment part
#if (PP_TimeDiscMethod==42)
  IF (Adsorption%TPD) THEN
    WallTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC(globSide))) + (Adsorption%TPD_beta * iter * dt)
    Adsorption%TPD_Temp = Walltemp
  ELSE
    WallTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC(globSide)))
  END IF
#else
  WallTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC(globSide)))
#endif

DO subsurfeta = 1,nSurfSample
DO subsurfxi = 1,nSurfSample
  desorbnum(:) = 0
  reactdesorbnum(:) = 0
  nSites(:) = 0
  nSitesRemain(:) = 0
  remainNum(:) = 0
  P_des(:) = 0.
  P_react(:) = 0.
#if (PP_TimeDiscMethod==42)
  IF (Adsorption%TPD) Energy(:) = 0.
#endif

  DO Coord=1,3
    nSites(Coord) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Coord)
    nSitesRemain(Coord) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord)
    nSites(4) = nSites(4) + nSites(Coord)
    nSitesRemain(4) = nSitesRemain(4) + nSitesRemain(Coord)
  END DO
  adsorbnum(:) = nSites(:) - nSitesRemain(:)
  traceNum = adsorbnum(4) ! number of all surface particles to be considered
  trace = 1 ! trace number for tracking of already considered surface particles (new adsorbates resulting from reaction not needed)
  adsorbates(:) = 0

  ! calculate number of adsorbates for each species (for analyze)
  DO Coord = 1,3
  DO AdsorbID = 1,nSites(Coord)-nSitesRemain(Coord)
    Surfpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+AdsorbID)
    iSpec = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(Surfpos)
    adsorbates(iSpec) = adsorbates(iSpec) + 1
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
    
    Surfpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID)
    iSpec = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(Surfpos)
    
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = &
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSites(Coord)-remainNum(Coord))
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSites(Coord)-remainNum(Coord)) = Surfpos
    AdsorbID = nSites(Coord)-remainNum(Coord)
    
    ! calculate heat of adsorption for actual site
    Heat_A = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,iSpec,Surfpos,.FALSE.)

    P_actual_react(:) = 0.
    Coord_ReactP(:) = -1
    Coord_Product(:,:) = -1
    Pos_ReactP(:) = 0
    Pos_Product(:,:) = 0
    ReactNum_run = 0
    ! Choose Random surface positions appropriate for each possible Reaction and calculate reaction probability
    ! First do associative reactions
    ! product of associative reactions desorb after reaction (LH-reaction)
    DO iReact = 1,(Adsorption%ReactNum-Adsorption%DissNum)
      IF (Adsorption%AssocReact(1,iReact,iSpec).LT.1) CYCLE ! no partner for this associative reaction
      Coord_ReactP(iReact) = Adsorption%Coordination(Adsorption%AssocReact(1,iReact,iSpec))
      CALL RANDOM_NUMBER(RanNum)
      IF (Coord.EQ.Coord_ReactP(iReact)) THEN
        NeighID = 1 + INT((nSites(Coord_ReactP(iReact))-remainNum(Coord_ReactP(iReact))-1)*RanNum)
      ELSE
        NeighID = 1 + INT((nSites(Coord_ReactP(iReact))-remainNum(Coord_ReactP(iReact)))*RanNum)
      END IF
      Pos_ReactP(iReact) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord_ReactP(iReact))%UsedSiteMap(NeighID)
      jSpec = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord_ReactP(iReact))%Species(Pos_ReactP(iReact))
      
      IF ( jSpec.EQ.Adsorption%AssocReact(1,iReact,iSpec) .AND. (jSpec.GT.0)) THEN
        Prod_Spec1 = Adsorption%AssocReact(2,iReact,iSpec)
        ! calculate heats of adsorption
!         Heat_AB = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Prod_Spec1,-1,.FALSE.)
        Heat_AB = 0. ! immediate desorption
        Heat_B = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,jSpec,Pos_ReactP(iReact),.FALSE.)
        D_AB = Adsorption%EDissBond((Adsorption%DissNum+iReact),iSpec)
        ! calculate LH reaction probability
        E_d = Calc_E_Act(Heat_AB,0.,Heat_A,Heat_B,D_AB,0.,0.,0.)
        CALL PartitionFuncAct(Prod_Spec1, WallTemp, VarPartitionFuncAct, Adsorption%DensSurfAtoms(SurfSideID)/Adsorption%AreaIncrease)
        CALL PartitionFuncSurf(iSpec, WallTemp, VarPartitionFuncWall1)
        CALL PartitionFuncSurf(jSpec, WallTemp, VarPartitionFuncWall2)
        nu_react = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct &
                    / (VarPartitionFuncWall1 * VarPartitionFuncWall2))
        rate = nu_react * exp(-E_d/(WallTemp))
        P_actual_react(iReact) = rate * dt
        IF (rate*dt.GT.1) THEN
          P_react(iReact) = P_react(iReact) + 1
          P_des(Prod_Spec1) = P_des(Prod_Spec1) + 1
        ELSE
          P_react(iReact) = P_react(iReact) + P_actual_react(iReact)
          P_des(Prod_Spec1) = P_des(Prod_Spec1) + P_actual_react(iReact)
        END IF
      END IF
    END DO
    ReactNum_run = ReactNum_run + Adsorption%ReactNum-Adsorption%DissNum
    ! Then do dissociative reactions (both products stay adsorbed)
    
    
!     ! find neighbour positions with empty site of same coord for diffusion
!     n_empty_Neigh = 0
!     ALLOCATE(NeighbourID(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nNeighbours))
!     DO i = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nNeighbours
!       Coord3 = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighSite(Surfpos,i)
!       IF (Coord3.GT.0) THEN
!         IF ( (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord3)%Species( &
!               SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,i)) &
!               .EQ.0) ) THEN
!           IF (Coord3.EQ.Coord) THEN
!             n_empty_Neigh = n_empty_Neigh + 1
!             NeighbourID(n_empty_Neigh) = i
!           END IF
!         END IF
!       END IF
!     END DO
    
    !-------------------------------------------------------------------------------------------------------------------------------
    ! sort Neighbours to coordinations for search of two random neighbour positions from particle position
    ! for dissociation and diffusion
    !-------------------------------------------------------------------------------------------------------------------------------
    n_Neigh(:) = 0
    n_empty_Neigh(:) = 0
    ALLOCATE(NeighbourID(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nNeighbours))
    ALLOCATE(React_NeighbourID(1:3,1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nNeighbours))

    DO i = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nNeighbours
      Coord3 = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighSite(Surfpos,i)
      IF (Coord3.GT.0) THEN
        IF ( (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord3)%Species( &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,i)) &
              .EQ.0) ) THEN
          n_empty_Neigh(Coord3) = n_empty_Neigh(Coord3) + 1
          IF (Coord3.EQ.Coord) THEN
            NeighbourID(n_empty_Neigh(Coord3)) = i
          END IF
        END IF
        n_Neigh(Coord3) = n_Neigh(Coord3) + 1
        React_NeighbourID(Coord3,n_Neigh(Coord3)) = i
      END IF
    END DO
    !-------------------------------------------------------------------------------------------------------------------------------
    ! calculate probability for dissociation
    !-------------------------------------------------------------------------------------------------------------------------------
    ! This routine chooses random empty sites from remaining sites array and then calculates dissociation probability
    IF (SpecDSMC(iSpec)%InterID.EQ.2) THEN
    DO iReact = 1,Adsorption%DissNum
      Prod_Spec1 = Adsorption%DissocReact(1,iReact,iSpec)
      Prod_Spec2 = Adsorption%DissocReact(2,iReact,iSpec)
      IF ((Prod_Spec1.NE.0) .AND. (Prod_Spec2.NE.0)) THEN !if 2 resulting species, dissociation possible
        Coord_Product(1,iReact+ReactNum_run) = Adsorption%Coordination(Prod_Spec1)
        Coord_Product(2,iReact+ReactNum_run) = Adsorption%Coordination(Prod_Spec2)
        
        IF ((Coord_Product(1,iReact+ReactNum_run).EQ.Coord_Product(2,iReact+ReactNum_run))&
              .AND.(nSitesRemain(Coord_Product(1,iReact+ReactNum_run)).LT.2)) CYCLE
        IF ((Coord_Product(1,iReact+ReactNum_run).NE.Coord_Product(2,iReact+ReactNum_run))&
              .AND.((nSitesRemain(Coord_Product(1,iReact+ReactNum_run)).LT.1)&
              .OR.(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)).LT.1))) CYCLE
        IF (Coord_Product(1,iReact+ReactNum_run).EQ.Coord_Product(2,iReact+ReactNum_run)) THEN
          CALL RANDOM_NUMBER(RanNum)
          NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,iReact+ReactNum_run)))/2.)*RanNum)
          Pos_Product(1,iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                            Coord_Product(1,iReact+ReactNum_run))%UsedSiteMap(NeighID)
          CALL RANDOM_NUMBER(RanNum)
          NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)))/2.)*RanNum &
                      + (REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)))/2.))
          Pos_Product(2,iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                            Coord_Product(2,iReact+ReactNum_run))%UsedSiteMap(NeighID)
        ELSE
          CALL RANDOM_NUMBER(RanNum)
          NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,iReact+ReactNum_run))))*RanNum)
          Pos_Product(1,iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                            Coord_Product(1,iReact+ReactNum_run))%UsedSiteMap(NeighID)
          CALL RANDOM_NUMBER(RanNum)
          NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run))))*RanNum)
          Pos_Product(2,iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                            Coord_Product(2,iReact+ReactNum_run))%UsedSiteMap(NeighID)
        END IF
        ! calculation of activation energy
        Heat_Product1 = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Prod_Spec1,Pos_Product(1,iReact+ReactNum_run),.TRUE.)
        Heat_Product2 = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Prod_Spec2,Pos_Product(2,iReact+ReactNum_run),.TRUE.)
        D_A = Adsorption%EDissBond(iReact,iSpec)
        D_Product1 = 0.
        D_Product2 = 0.
        E_d = Calc_E_Act(Heat_Product1,Heat_Product2,Heat_A,0.,D_Product1,D_Product2,D_A,0.)
        CALL PartitionFuncSurf(iSpec, WallTemp, VarPartitionFuncWall)
        CALL PartitionFuncAct(Prod_Spec1, WallTemp, VarPartitionFuncAct1, &
                              Adsorption%DensSurfAtoms(SurfSideID)/Adsorption%AreaIncrease)
        CALL PartitionFuncAct(Prod_Spec2, WallTemp, VarPartitionFuncAct2, &
                              Adsorption%DensSurfAtoms(SurfSideID)/Adsorption%AreaIncrease)
        nu_react = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct1 * VarPartitionFuncAct2 &
                    / (VarPartitionFuncWall))
        rate = nu_react * exp(-E_d/(WallTemp))
        P_actual_react(ReactNum_run+iReact) = rate * dt
        IF (rate*dt.GT.1) THEN
          P_react(ReactNum_run+iReact) = P_react(ReactNum_run+iReact) + 1.
        ELSE
          P_react(ReactNum_run+iReact) = P_react(ReactNum_run+iReact) + P_actual_react(ReactNum_run+iReact)
        END IF
      END IF
    END DO
    END IF
    ReactNum_run = ReactNum_run + Adsorption%DissNum

! This routine looks if neirest neighbours are empty and then calculates dissociation probability
!     IF (SpecDSMC(iSpec)%InterID.EQ.2) THEN
!     DO iReact = 1,(Adsorption%DissNum)
!       Prod_Spec1 = Adsorption%DissocReact(1,iReact,iSpec)
!       Prod_Spec2 = Adsorption%DissocReact(2,iReact,iSpec)
!         jCoord = Adsorption%Coordination(Prod_Spec1)
!         kCoord = Adsorption%Coordination(Prod_Spec2)
!         Neighpos_j = 0
!         Neighpos_k = 0
!         !At least one of resulting species has same coordination as surfaceposition of dissociating particle
!         IF ( (jCoord.EQ.Coord) .OR. (kCoord.EQ.Coord) ) THEN 
!           IF (jCoord.EQ.Coord) THEN
!             Neighpos_j = Surfpos
!             IF (n_empty_Neigh(kCoord).GT.0) THEN
!               ! assign availiable neighbour position for k
!               CALL RANDOM_NUMBER(RanNum)
!               chosen_Neigh_k = 1 + INT(REAL(n_Neigh(kCoord))*RanNum)
!               Neighpos_k = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,&
!                                                                                       React_NeighbourID(kCoord,chosen_Neigh_k))
!               IDRearrange = React_NeighbourID(jCoord,chosen_Neigh_k)
!               React_NeighbourID(jCoord,chosen_Neigh_k) = React_NeighbourID(kCoord,n_Neigh(kCoord))
!               React_NeighbourID(jCoord,n_Neigh(jCoord)) = IDRearrange
!             END IF
!           ELSE
!             Neighpos_k = Surfpos
!             IF (n_empty_Neigh(jCoord).GT.0) THEN
!               ! assign availiable neighbour position for j
!               CALL RANDOM_NUMBER(RanNum)
!               chosen_Neigh_j = 1 + INT(REAL(n_Neigh(jCoord))*RanNum)
!               Neighpos_j = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,&
!                                                                                       React_NeighbourID(jCoord,chosen_Neigh_j))
!               IDRearrange = React_NeighbourID(jCoord,chosen_Neigh_j)
!               React_NeighbourID(jCoord,chosen_Neigh_j) = React_NeighbourID(jCoord,n_Neigh(jCoord))
!               React_NeighbourID(jCoord,n_Neigh(jCoord)) = IDRearrange
!             END IF
!           END IF
!         !NEITHER of resulting species have same coordination as surfaceposition of dissociating particle
!         ELSE
!           IF ( (jCoord.EQ.kCoord) .AND. (n_empty_Neigh(jCoord).GT.1) ) THEN !resulting species are equal
!             ! assign availiable neighbour positions
!             CALL RANDOM_NUMBER(RanNum)
!             chosen_Neigh_j = 1 + INT(REAL(n_Neigh(jCoord))*RanNum)
!             Neighpos_j = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,&
!                                                                                     React_NeighbourID(jCoord,chosen_Neigh_j))
!             IDRearrange = React_NeighbourID(jCoord,chosen_Neigh_j)
!             React_NeighbourID(jCoord,chosen_Neigh_j) = React_NeighbourID(jCoord,n_Neigh(jCoord))
!             React_NeighbourID(jCoord,n_Neigh(jCoord)) = IDRearrange
!             CALL RANDOM_NUMBER(RanNum)
!             chosen_Neigh_k = 1 + INT(REAL(n_Neigh(kCoord)-1)*RanNum)
!             Neighpos_k = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,&
!                                                                                     React_NeighbourID(kCoord,chosen_Neigh_k))
!             IDRearrange = React_NeighbourID(jCoord,chosen_Neigh_j)
!             React_NeighbourID(jCoord,chosen_Neigh_k) = React_NeighbourID(kCoord,n_Neigh(kCoord)-1)
!             React_NeighbourID(jCoord,n_Neigh(jCoord)-1) = IDRearrange
!           ELSE IF ( (jCoord.NE.kCoord) .AND. (n_empty_Neigh(jCoord).GT.0) .AND. (n_empty_Neigh(kCoord).GT.0) ) THEN
!             ! assign availiable neighbour positions
!             CALL RANDOM_NUMBER(RanNum)
!             chosen_Neigh_j = 1 + INT(REAL(n_Neigh(jCoord))*RanNum)
!             Neighpos_j = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,&
!                                                                                     React_NeighbourID(jCoord,chosen_Neigh_j))
!             IDRearrange = React_NeighbourID(jCoord,chosen_Neigh_j)
!             React_NeighbourID(jCoord,chosen_Neigh_j) = React_NeighbourID(jCoord,n_Neigh(jCoord))
!             React_NeighbourID(jCoord,n_Neigh(jCoord)) = IDRearrange
!             CALL RANDOM_NUMBER(RanNum)
!             chosen_Neigh_k = 1 + INT(REAL(n_Neigh(kCoord))*RanNum)
!             Neighpos_k = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,&
!                                                                                     React_NeighbourID(kCoord,chosen_Neigh_k))
!             IDRearrange = React_NeighbourID(jCoord,chosen_Neigh_j)
!             React_NeighbourID(jCoord,chosen_Neigh_k) = React_NeighbourID(kCoord,n_Neigh(kCoord))
!             React_NeighbourID(jCoord,n_Neigh(jCoord)) = IDRearrange
!           END IF
!         END IF
!         IF ( (Neighpos_j.GT.0) .AND. (Neighpos_k.GT.0) ) THEN !both neighbour positions assigned
!         !both assigned neighbour positions empty
!         IF ( (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%Species(Neighpos_j).EQ.0) .AND. &
!             (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(kCoord)%Species(Neighpos_k).EQ.0) ) THEN
!           ! calculation of activation energy
!           IF (Neighpos_j.EQ.Surfpos) THEN
!             Heat_A = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Prod_Spec1,Neighpos_j,.FALSE.)
!             Heat_B = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Prod_Spec2,Neighpos_k,.TRUE.)
!           ELSE IF (Neighpos_k.EQ.Surfpos) THEN
!             Heat_A = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Prod_Spec1,Neighpos_j,.TRUE.)
!             Heat_B = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Prod_Spec2,Neighpos_k,.FALSE.)
!           ELSE
!             Heat_A = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Prod_Spec1,Neighpos_j,.TRUE.)
!             Heat_B = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Prod_Spec2,Neighpos_k,.TRUE.)
!           END IF
!           Heat_AB = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,iSpec,Surfpos,.FALSE.)
!           D_AB = Adsorption%EDissBond(iReact,iSpec)
!           D_A = 0.
!           D_B = 0.
!           E_d = Calc_E_Act(Heat_A,Heat_B,Heat_AB,0.,D_A,D_B,D_AB,0.)
!           CALL PartitionFuncSurf(iSpec, WallTemp, VarPartitionFuncWall)
!           CALL PartitionFuncAct(Prod_Spec1, WallTemp, VarPartitionFuncAct1, Adsorption%DensSurfAtoms(SurfSideID)/Adsorption%AreaIncrease)
!           CALL PartitionFuncAct(Prod_Spec2, WallTemp, VarPartitionFuncAct2, Adsorption%DensSurfAtoms(SurfSideID)/Adsorption%AreaIncrease)
!           nu_react = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct1 * VarPartitionFuncAct2 &
!                       / (VarPartitionFuncWall))
!           rate = nu_react * exp(-E_d/(WallTemp))
!           P_actual_react(ReactNum_run+iReact) = rate * dt
!           IF (rate*dt.GT.1) THEN
!             P_react(ReactNum_run+iReact) = P_react(ReactNum_run+iReact) + 1.
!           ELSE
!             P_react(ReactNum_run+iReact) = P_react(ReactNum_run+iReact) + P_actual_react(ReactNum_run+iReact)
!           END IF
!         END IF
!         END IF
!       END IF
!     END DO
!     END IF
!     ReactNum_run = ReactNum_run + iReact - 1
    
    ! Calculate propabilities for exchange reactions 
    ! maybe try later (if exothermic --> Delta_H < 0, product with lower adsorbheat desorbs else both stay adsorbed)
    DO iReact = 1,Adsorption%nDisPropReactions
      IF ( (Adsorption%ChemReactant(1,iReact+Adsorption%nDissocReactions).NE.iSpec) &
           .AND. (Adsorption%ChemReactant(2,iReact+Adsorption%nDissocReactions).NE.iSpec) ) THEN
        IF ( (Adsorption%ChemProduct(1,iReact+Adsorption%nDissocReactions).NE.iSpec) &
           .AND. (Adsorption%ChemProduct(2,iReact+Adsorption%nDissocReactions).NE.iSpec) ) CYCLE
      END IF
      ! Choose which reaction partner considered particle is 
      ! -> defines which species is second reaction partner and if forward or reverse reaction
      ! ----------------------------------------------------------------------------------------------------------------------------
      IF (Adsorption%ChemReactant(1,iReact+Adsorption%nDissocReactions).EQ.iSpec) THEN ! -> first forward
      ! ----------------------------------------------------------------------------------------------------------------------------
        Coord_ReactP(iReact+ReactNum_run) = Adsorption%Coordination(Adsorption%ChemReactant(2,iReact+Adsorption%nDissocReactions))
        CALL RANDOM_NUMBER(RanNum)
        IF (Coord.EQ.Coord_ReactP(iReact)) THEN
          NeighID = 1 + INT((nSites(Coord_ReactP(iReact+ReactNum_run))-remainNum(Coord_ReactP(iReact+ReactNum_run))-1)*RanNum)
        ELSE
          NeighID = 1 + INT((nSites(Coord_ReactP(iReact+ReactNum_run))-remainNum(Coord_ReactP(iReact+ReactNum_run)))*RanNum)
        END IF
        Pos_ReactP(iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                          Coord_ReactP(iReact+ReactNum_run))%UsedSiteMap(NeighID)
        jSpec = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                Coord_ReactP(iReact+ReactNum_run))%Species(Pos_ReactP(iReact+ReactNum_run))
        ! continue only if species on chosen site is appropriate partner
        IF ((jSpec.EQ.Adsorption%ChemReactant(2,iReact+Adsorption%nDissocReactions)).AND.(jSpec.GT.0)) THEN
          Prod_Spec1 = Adsorption%ChemProduct(1,iReact+Adsorption%nDissocReactions)
          Coord_Product(1,iReact+ReactNum_run) = Adsorption%Coordination(Prod_Spec1)
          Prod_Spec2 = Adsorption%ChemProduct(2,iReact+Adsorption%nDissocReactions)
          Coord_Product(2,iReact+ReactNum_run) = Adsorption%Coordination(Prod_Spec2)
          IF ((Coord_Product(1,iReact+ReactNum_run).EQ.Coord_Product(2,iReact+ReactNum_run))&
               .AND.(nSitesRemain(Coord_Product(1,iReact+ReactNum_run)).LT.2)) CYCLE
          IF ((Coord_Product(1,iReact+ReactNum_run).NE.Coord_Product(2,iReact+ReactNum_run))&
               .AND.((nSitesRemain(Coord_Product(1,iReact+ReactNum_run)).LT.1)&
               .OR.(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)).LT.1))) CYCLE
          IF (Coord_Product(1,iReact+ReactNum_run).EQ.Coord_Product(2,iReact+ReactNum_run)) THEN
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,iReact+ReactNum_run)))/2.)*RanNum)
            Pos_Product(1,iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                              Coord_Product(1,iReact+ReactNum_run))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)))/2.)*RanNum &
                        + (REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)))/2.))
            Pos_Product(2,iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                              Coord_Product(2,iReact+ReactNum_run))%UsedSiteMap(NeighID)
          ELSE
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,iReact+ReactNum_run))))*RanNum)
            Pos_Product(1,iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                              Coord_Product(1,iReact+ReactNum_run))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run))))*RanNum)
            Pos_Product(2,iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                              Coord_Product(2,iReact+ReactNum_run))%UsedSiteMap(NeighID)
          END IF
          ! calculate exchange reaction probability
          Heat_ReactP = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,jSpec,Pos_ReactP(iReact+ReactNum_run),.FALSE.)
          Heat_Product1 = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Prod_Spec1,Pos_Product(1,iReact+ReactNum_run),.TRUE.)
          Heat_Product2 = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Prod_Spec2,Pos_Product(2,iReact+ReactNum_run),.TRUE.)
          D_A = Adsorption%Reactant_DissBond_K(1,iReact+Adsorption%nDissocReactions)
          D_ReactP = Adsorption%Reactant_DissBond_K(2,iReact+Adsorption%nDissocReactions)
          D_Product1 = Adsorption%Product_DissBond_K(1,iReact+Adsorption%nDissocReactions)
          D_Product2 = Adsorption%Product_DissBond_K(2,iReact+Adsorption%nDissocReactions)
          E_d = Calc_E_Act(Heat_Product1,Heat_Product2,Heat_A,Heat_ReactP,D_Product1,D_Product2,D_A,D_ReactP)
          CALL PartitionFuncAct(Prod_Spec1, WallTemp, VarPartitionFuncAct1, &
                                Adsorption%DensSurfAtoms(SurfSideID)/Adsorption%AreaIncrease)
          CALL PartitionFuncAct(Prod_Spec2, WallTemp, VarPartitionFuncAct2, &
                                Adsorption%DensSurfAtoms(SurfSideID)/Adsorption%AreaIncrease)
          CALL PartitionFuncSurf(iSpec, WallTemp, VarPartitionFuncWall1)
          CALL PartitionFuncSurf(jSpec, WallTemp, VarPartitionFuncWall2)
          nu_react = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct1 * VarPartitionFuncAct2 &
                      / (VarPartitionFuncWall1 * VarPartitionFuncWall2))
          rate = nu_react * exp(-E_d/(WallTemp))
          P_actual_react(iReact+ReactNum_run) = rate * dt
          IF (rate*dt.GT.1) THEN
            P_react(iReact+ReactNum_run) = P_react(iReact+ReactNum_run) + 1
          ELSE
            P_react(iReact+ReactNum_run) = P_react(iReact+ReactNum_run) + P_actual_react(iReact+ReactNum_run)
          END IF
        END IF
      ! ----------------------------------------------------------------------------------------------------------------------------
      ELSE IF (Adsorption%ChemReactant(2,iReact+Adsorption%nDissocReactions).EQ.iSpec) THEN ! -> second forward
      ! ----------------------------------------------------------------------------------------------------------------------------
        Coord_ReactP(iReact+ReactNum_run) = Adsorption%Coordination(Adsorption%ChemReactant(1,iReact+Adsorption%nDissocReactions))
        CALL RANDOM_NUMBER(RanNum)
        IF (Coord.EQ.Coord_ReactP(iReact)) THEN
          NeighID = 1 + INT((nSites(Coord_ReactP(iReact+ReactNum_run))-remainNum(Coord_ReactP(iReact+ReactNum_run))-1)*RanNum)
        ELSE
          NeighID = 1 + INT((nSites(Coord_ReactP(iReact+ReactNum_run))-remainNum(Coord_ReactP(iReact+ReactNum_run)))*RanNum)
        END IF
        Pos_ReactP(iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                          Coord_ReactP(iReact+ReactNum_run))%UsedSiteMap(NeighID)
        jSpec = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                Coord_ReactP(iReact+ReactNum_run))%Species(Pos_ReactP(iReact+ReactNum_run))
        ! continue only if species on chosen site is appropriate partner
        IF ((jSpec.EQ.Adsorption%ChemReactant(1,iReact+Adsorption%nDissocReactions)).AND.(jSpec.GT.0)) THEN
          Prod_Spec1 = Adsorption%ChemProduct(1,iReact+Adsorption%nDissocReactions)
          Coord_Product(1,iReact+ReactNum_run) = Adsorption%Coordination(Prod_Spec1)
          Prod_Spec2 = Adsorption%ChemProduct(2,iReact+Adsorption%nDissocReactions)
          Coord_Product(2,iReact+ReactNum_run) = Adsorption%Coordination(Prod_Spec2)
          IF ((Coord_Product(1,iReact+ReactNum_run).EQ.Coord_Product(2,iReact+ReactNum_run))&
               .AND.(nSitesRemain(Coord_Product(1,iReact+ReactNum_run)).LT.2)) CYCLE
          IF ((Coord_Product(1,iReact+ReactNum_run).NE.Coord_Product(2,iReact+ReactNum_run))&
               .AND.((nSitesRemain(Coord_Product(1,iReact+ReactNum_run)).LT.1)&
               .OR.(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)).LT.1))) CYCLE
          IF (Coord_Product(1,iReact+ReactNum_run).EQ.Coord_Product(2,iReact+ReactNum_run)) THEN
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,iReact+ReactNum_run)))/2.)*RanNum)
            Pos_Product(1,iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                              Coord_Product(1,iReact+ReactNum_run))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)))/2.)*RanNum &
                        + (REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)))/2.))
            Pos_Product(2,iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                              Coord_Product(2,iReact+ReactNum_run))%UsedSiteMap(NeighID)
          ELSE
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,iReact+ReactNum_run))))*RanNum)
            Pos_Product(1,iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                              Coord_Product(1,iReact+ReactNum_run))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run))))*RanNum)
            Pos_Product(2,iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                              Coord_Product(2,iReact+ReactNum_run))%UsedSiteMap(NeighID)
          END IF
          ! calculate exchange reaction probability
          Heat_ReactP = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,jSpec,Pos_ReactP(iReact+ReactNum_run),.FALSE.)
          Heat_Product1 = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Prod_Spec1,Pos_Product(1,iReact+ReactNum_run),.TRUE.)
          Heat_Product2 = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Prod_Spec2,Pos_Product(2,iReact+ReactNum_run),.TRUE.)
          D_A = Adsorption%Reactant_DissBond_K(2,iReact+Adsorption%nDissocReactions)
          D_ReactP = Adsorption%Reactant_DissBond_K(1,iReact+Adsorption%nDissocReactions)
          D_Product1 = Adsorption%Product_DissBond_K(1,iReact+Adsorption%nDissocReactions)
          D_Product2 = Adsorption%Product_DissBond_K(2,iReact+Adsorption%nDissocReactions)
          E_d = Calc_E_Act(Heat_Product1,Heat_Product2,Heat_A,Heat_ReactP,D_Product1,D_Product2,D_A,D_ReactP)
          CALL PartitionFuncAct(Prod_Spec1, WallTemp, VarPartitionFuncAct1, &
                                Adsorption%DensSurfAtoms(SurfSideID)/Adsorption%AreaIncrease)
          CALL PartitionFuncAct(Prod_Spec2, WallTemp, VarPartitionFuncAct2, &
                                Adsorption%DensSurfAtoms(SurfSideID)/Adsorption%AreaIncrease)
          CALL PartitionFuncSurf(jSpec, WallTemp, VarPartitionFuncWall1)
          CALL PartitionFuncSurf(iSpec, WallTemp, VarPartitionFuncWall2)
          nu_react = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct1 * VarPartitionFuncAct2 &
                      / (VarPartitionFuncWall1 * VarPartitionFuncWall2))
          rate = nu_react * exp(-E_d/(WallTemp))
          P_actual_react(iReact+ReactNum_run) = rate * dt
          IF (rate*dt.GT.1) THEN
            P_react(iReact+ReactNum_run) = P_react(iReact+ReactNum_run) + 1
          ELSE
            P_react(iReact+ReactNum_run) = P_react(iReact+ReactNum_run) + P_actual_react(iReact+ReactNum_run)
          END IF
        END IF
      ! ----------------------------------------------------------------------------------------------------------------------------
      ELSE IF (Adsorption%ChemProduct(1,iReact+Adsorption%nDissocReactions).EQ.iSpec) THEN ! -> first reverse
      ! ----------------------------------------------------------------------------------------------------------------------------
        Coord_ReactP(iReact+ReactNum_run) = Adsorption%Coordination(Adsorption%ChemProduct(2,iReact+Adsorption%nDissocReactions))
        CALL RANDOM_NUMBER(RanNum)
        IF (Coord.EQ.Coord_ReactP(iReact)) THEN
          NeighID = 1 + INT((nSites(Coord_ReactP(iReact+ReactNum_run))-remainNum(Coord_ReactP(iReact+ReactNum_run))-1)*RanNum)
        ELSE
          NeighID = 1 + INT((nSites(Coord_ReactP(iReact+ReactNum_run))-remainNum(Coord_ReactP(iReact+ReactNum_run)))*RanNum)
        END IF
        Pos_ReactP(iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                          Coord_ReactP(iReact+ReactNum_run))%UsedSiteMap(NeighID)
        jSpec = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                Coord_ReactP(iReact+ReactNum_run))%Species(Pos_ReactP(iReact+ReactNum_run))
        ! continue only if species on chosen site is appropriate partner
        IF ((jSpec.EQ.Adsorption%ChemProduct(2,iReact+Adsorption%nDissocReactions)).AND.(jSpec.GT.0)) THEN
          Prod_Spec1 = Adsorption%ChemReactant(1,iReact+Adsorption%nDissocReactions)
          Coord_Product(1,iReact+ReactNum_run) = Adsorption%Coordination(Prod_Spec1)
          Prod_Spec2 = Adsorption%ChemReactant(2,iReact+Adsorption%nDissocReactions)
          Coord_Product(2,iReact+ReactNum_run) = Adsorption%Coordination(Prod_Spec2)
          IF ((Coord_Product(1,iReact+ReactNum_run).EQ.Coord_Product(2,iReact+ReactNum_run))&
               .AND.(nSitesRemain(Coord_Product(1,iReact+ReactNum_run)).LT.2)) CYCLE
          IF ((Coord_Product(1,iReact+ReactNum_run).NE.Coord_Product(2,iReact+ReactNum_run))&
               .AND.((nSitesRemain(Coord_Product(1,iReact+ReactNum_run)).LT.1)&
               .OR.(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)).LT.1))) CYCLE
          IF (Coord_Product(1,iReact+ReactNum_run).EQ.Coord_Product(2,iReact+ReactNum_run)) THEN
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,iReact+ReactNum_run)))/2.)*RanNum)
            Pos_Product(1,iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                              Coord_Product(1,iReact+ReactNum_run))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)))/2.)*RanNum &
                        + (REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)))/2.))
            Pos_Product(2,iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                              Coord_Product(2,iReact+ReactNum_run))%UsedSiteMap(NeighID)
          ELSE
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,iReact+ReactNum_run))))*RanNum)
            Pos_Product(1,iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                              Coord_Product(1,iReact+ReactNum_run))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run))))*RanNum)
            Pos_Product(2,iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                              Coord_Product(2,iReact+ReactNum_run))%UsedSiteMap(NeighID)
          END IF
          ! calculate exchange reaction probability
          Heat_ReactP = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,jSpec,Pos_ReactP(iReact+ReactNum_run),.FALSE.)
          Heat_Product1 = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Prod_Spec1,Pos_Product(1,iReact+ReactNum_run),.TRUE.)
          Heat_Product2 = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Prod_Spec2,Pos_Product(2,iReact+ReactNum_run),.TRUE.)
          D_A = Adsorption%Product_DissBond_K(1,iReact+Adsorption%nDissocReactions)
          D_ReactP = Adsorption%Product_DissBond_K(2,iReact+Adsorption%nDissocReactions)
          D_Product1 = Adsorption%Reactant_DissBond_K(1,iReact+Adsorption%nDissocReactions)
          D_Product2 = Adsorption%Reactant_DissBond_K(2,iReact+Adsorption%nDissocReactions)
          E_d = Calc_E_Act(Heat_Product1,Heat_Product2,Heat_A,Heat_ReactP,D_Product1,D_Product2,D_A,D_ReactP)
          CALL PartitionFuncAct(Prod_Spec1, WallTemp, VarPartitionFuncAct1, &
                                Adsorption%DensSurfAtoms(SurfSideID)/Adsorption%AreaIncrease)
          CALL PartitionFuncAct(Prod_Spec2, WallTemp, VarPartitionFuncAct2, &
                                Adsorption%DensSurfAtoms(SurfSideID)/Adsorption%AreaIncrease)
          CALL PartitionFuncSurf(jSpec, WallTemp, VarPartitionFuncWall1)
          CALL PartitionFuncSurf(iSpec, WallTemp, VarPartitionFuncWall2)
          nu_react = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct1 * VarPartitionFuncAct2 &
                      / (VarPartitionFuncWall1 * VarPartitionFuncWall2))
          rate = nu_react * exp(-E_d/(WallTemp))
          P_actual_react(iReact+ReactNum_run) = rate * dt
          IF (rate*dt.GT.1) THEN
            P_react(iReact+ReactNum_run) = P_react(iReact+ReactNum_run) + 1
          ELSE
            P_react(iReact+ReactNum_run) = P_react(iReact+ReactNum_run) + P_actual_react(iReact+ReactNum_run)
          END IF
        END IF
      ! ----------------------------------------------------------------------------------------------------------------------------
      ELSE IF (Adsorption%ChemProduct(2,iReact+Adsorption%nDissocReactions).EQ.iSpec) THEN ! -> second reverse
      ! ----------------------------------------------------------------------------------------------------------------------------
        Coord_ReactP(iReact+ReactNum_run) = Adsorption%Coordination(Adsorption%ChemProduct(1,iReact+Adsorption%nDissocReactions))
        CALL RANDOM_NUMBER(RanNum)
        IF (Coord.EQ.Coord_ReactP(iReact)) THEN
          NeighID = 1 + INT((nSites(Coord_ReactP(iReact+ReactNum_run))-remainNum(Coord_ReactP(iReact+ReactNum_run))-1)*RanNum)
        ELSE
          NeighID = 1 + INT((nSites(Coord_ReactP(iReact+ReactNum_run))-remainNum(Coord_ReactP(iReact+ReactNum_run)))*RanNum)
        END IF
        Pos_ReactP(iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                          Coord_ReactP(iReact+ReactNum_run))%UsedSiteMap(NeighID)
        jSpec = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                Coord_ReactP(iReact+ReactNum_run))%Species(Pos_ReactP(iReact+ReactNum_run))
        ! continue only if species on chosen site is appropriate partner
        IF ((jSpec.EQ.Adsorption%ChemProduct(1,iReact+Adsorption%nDissocReactions)).AND.(jSpec.GT.0)) THEN
          Prod_Spec1 = Adsorption%ChemReactant(1,iReact+Adsorption%nDissocReactions)
          Coord_Product(1,iReact+ReactNum_run) = Adsorption%Coordination(Prod_Spec1)
          Prod_Spec2 = Adsorption%ChemReactant(2,iReact+Adsorption%nDissocReactions)
          Coord_Product(2,iReact+ReactNum_run) = Adsorption%Coordination(Prod_Spec2)
          IF ((Coord_Product(1,iReact+ReactNum_run).EQ.Coord_Product(2,iReact+ReactNum_run))&
               .AND.(nSitesRemain(Coord_Product(1,iReact+ReactNum_run)).LT.2)) CYCLE
          IF ((Coord_Product(1,iReact+ReactNum_run).NE.Coord_Product(2,iReact+ReactNum_run))&
               .AND.((nSitesRemain(Coord_Product(1,iReact+ReactNum_run)).LT.1)&
               .OR.(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)).LT.1))) CYCLE
          IF (Coord_Product(1,iReact+ReactNum_run).EQ.Coord_Product(2,iReact+ReactNum_run)) THEN
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,iReact+ReactNum_run)))/2.)*RanNum)
            Pos_Product(1,iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                              Coord_Product(1,iReact+ReactNum_run))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)))/2.)*RanNum &
                        + (REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run)))/2.))
            Pos_Product(2,iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                              Coord_Product(2,iReact+ReactNum_run))%UsedSiteMap(NeighID)
          ELSE
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(1,iReact+ReactNum_run))))*RanNum)
            Pos_Product(1,iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                              Coord_Product(1,iReact+ReactNum_run))%UsedSiteMap(NeighID)
            CALL RANDOM_NUMBER(RanNum)
            NeighID = 1 + INT((REAL(nSitesRemain(Coord_Product(2,iReact+ReactNum_run))))*RanNum)
            Pos_Product(2,iReact+ReactNum_run) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(&
                                              Coord_Product(2,iReact+ReactNum_run))%UsedSiteMap(NeighID)
          END IF
          ! calculate exchange reaction probability
          Heat_ReactP = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,jSpec,Pos_ReactP(iReact+ReactNum_run),.FALSE.)
          Heat_Product1 = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Prod_Spec1,Pos_Product(1,iReact+ReactNum_run),.TRUE.)
          Heat_Product2 = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Prod_Spec2,Pos_Product(2,iReact+ReactNum_run),.TRUE.)
          D_A = Adsorption%Product_DissBond_K(2,iReact+Adsorption%nDissocReactions)
          D_ReactP = Adsorption%Product_DissBond_K(1,iReact+Adsorption%nDissocReactions)
          D_Product1 = Adsorption%Reactant_DissBond_K(1,iReact+Adsorption%nDissocReactions)
          D_Product2 = Adsorption%Reactant_DissBond_K(2,iReact+Adsorption%nDissocReactions)
          E_d = Calc_E_Act(Heat_Product1,Heat_Product2,Heat_A,Heat_ReactP,D_Product1,D_Product2,D_A,D_ReactP)
          CALL PartitionFuncAct(Prod_Spec1, WallTemp, VarPartitionFuncAct1, &
                                Adsorption%DensSurfAtoms(SurfSideID)/Adsorption%AreaIncrease)
          CALL PartitionFuncAct(Prod_Spec2, WallTemp, VarPartitionFuncAct2, &
                                Adsorption%DensSurfAtoms(SurfSideID)/Adsorption%AreaIncrease)
          CALL PartitionFuncSurf(jSpec, WallTemp, VarPartitionFuncWall1)
          CALL PartitionFuncSurf(iSpec, WallTemp, VarPartitionFuncWall2)
          nu_react = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct1 * VarPartitionFuncAct2 &
                      / (VarPartitionFuncWall1 * VarPartitionFuncWall2))
          rate = nu_react * exp(-E_d/(WallTemp))
          P_actual_react(iReact+ReactNum_run) = rate * dt
          IF (rate*dt.GT.1) THEN
            P_react(iReact+ReactNum_run) = P_react(iReact+ReactNum_run) + 1
          ELSE
            P_react(iReact+ReactNum_run) = P_react(iReact+ReactNum_run) + P_actual_react(iReact+ReactNum_run)
          END IF
        END IF
      END IF
    END DO
    ReactNum_run = ReactNum_run + Adsorption%nDisPropReactions
    
    ! if empty neighbour sites are available choose one and calculate diffusion probability
    IF (n_empty_Neigh(Coord).GT.0) THEN
      ! choose random empty neighbour site from relevant neighbours
      CALL RANDOM_NUMBER(RanNum)
      react_Neigh = 1 + INT(n_empty_Neigh(Coord) * RanNum)
      Coord3 = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighSite(Surfpos,NeighbourID(react_Neigh))
      Pos_ReactP(ReactNum_run+1) = &
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,NeighbourID(react_Neigh))
      ! calculate heat of adsorption for actual site
      Heat_A = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,iSpec,Surfpos,.FALSE.)
      ! remove considered particle from bondorder
      DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
        Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(Surfpos,j)
        Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(Surfpos,j)
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) = &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) - 1
      END DO
      !calculate heat of adsorption for diffusion site
      Heat_B = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,iSpec,Pos_ReactP(ReactNum_run+1),.TRUE.)
      ! add considered particle to bondorder again
      DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
        Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(Surfpos,j)
        Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(Surfpos,j)
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) = &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) + 1
      END DO
      ! calculate diffusion probability
      P_actual_react(ReactNum_run+1) = exp(-(Heat_A - Heat_B)/WallTemp) / (1+exp(-(Heat_A - Heat_B)/Walltemp))
    END IF
    
    SDEALLOCATE(NeighbourID)
    SDEALLOCATE(React_NeighbourID)
    
    ! calculate molecular desorption probability
    CALL PartitionFuncAct(iSpec, WallTemp, VarPartitionFuncAct, Adsorption%DensSurfAtoms(SurfSideID)/Adsorption%AreaIncrease)
    CALL PartitionFuncSurf(iSpec, WallTemp, VarPartitionFuncWall)
    nu_des = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct / VarPartitionFuncWall)
    E_d = Calc_E_Act(Heat_A, 0.,0.,0.,0.,0.,0.,0.)
    rate = nu_des * exp(-E_d/WallTemp)
    P_actual_des = rate * dt
    IF (rate*dt.GT.1) THEN
      P_des(iSpec) = P_des(iSpec) + 1
    ELSE
      P_des(iSpec) = P_des(iSpec) + P_actual_des
    END IF
    
    ! initialize sum of all probabilities
    sum_probabilities = P_actual_des
    DO iReact = 1,ReactNum_run+1
      sum_probabilities = sum_probabilities + P_actual_react(iReact)
    END DO
    
    ! choose which surface reaction takes place
    CALL RANDOM_NUMBER(RanNum)
    IF (sum_probabilities .GT. RanNum) THEN
      DO iReact = 1,ReactNum_run+1
        CALL RANDOM_NUMBER(RanNum)
        IF ((P_actual_react(iReact)/sum_probabilities).GT.RanNum) THEN
          IF (iReact.LE.(Adsorption%ReactNum-Adsorption%DissNum)) THEN
            surf_react_case = 1
          ELSE IF ((iReact.GT.(Adsorption%ReactNum-Adsorption%DissNum)).AND.(iReact.LE.(Adsorption%ReactNum))) THEN
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
        IF (WriteMacroValues) THEN
          ! calculate Enthalpie of desorption and sample
          Heat_A = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,iSpec,-1,.FALSE.)
          Heat_B = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Adsorption%AssocReact(1,iReact,iSpec),-1,.FALSE.)
          Heat_AB = 0.
          D_AB = Adsorption%EDissBond(Adsorption%DissNum+iReact,iSpec)
          AdsorptionEnthalpie = (-(( Heat_AB -Heat_A -Heat_B ) + D_AB) * BoltzmannConst &
                          / REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(3))) &
                          * REAL(INT(Adsorption%DensSurfAtoms(SurfSideID) &
                          * SurfMesh%SurfaceArea(subsurfxi,subsurfeta,SurfSideID),8)) / Species(iSpec)%MacroParticleFactor
          !----  Sampling of energies
          IF ((DSMC%CalcSurfaceVal.AND.(Time.ge.(1-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroValues)) THEN
            SampWall(SurfSideID)%Adsorption(1,subsurfxi,subsurfeta) = SampWall(SurfSideID)%Adsorption(1,subsurfxi,subsurfeta) &
                                                                  + AdsorptionEnthalpie * Species(iSpec)%MacroParticleFactor
          END IF
        END IF
        ! remove adsorbate and update map
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(Surfpos) = 0
        DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
          Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(Surfpos,j)
          Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(Surfpos,j)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) = &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) - 1
        END DO
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = &
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1)
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1) = Surfpos
        nSitesRemain(Coord) = nSitesRemain(Coord) + 1
        nSitesRemain(4) = nSitesRemain(4) + 1
        ! remove adsorbate reaction partner and update map
        DO PartnerID = nSitesRemain(Coord_ReactP(iReact))+1,nSites(Coord_ReactP(iReact))
        IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord_ReactP(iReact))%UsedSiteMap(PartnerID) &
            .EQ.Pos_ReactP(iReact)) THEN
          jSpec = Adsorption%AssocReact(1,iReact,iSpec)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord_ReactP(iReact))%Species(Pos_ReactP(iReact)) = 0
          DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord_ReactP(iReact))%nInterAtom
            Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord_ReactP(iReact))%BondAtomIndx(Pos_ReactP(iReact),j)
            Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord_ReactP(iReact))%BondAtomIndy(Pos_ReactP(iReact),j)
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(jSpec,Indx,Indy) = &
                SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(jSpec,Indx,Indy) - 1
          END DO
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord_ReactP(iReact))%UsedSiteMap(PartnerID) = &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord_ReactP(iReact))%UsedSiteMap( &
              nSitesRemain(Coord_ReactP(iReact))+1)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord_ReactP(iReact))%UsedSiteMap( &
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
        DissocNum = iReact -(Adsorption%ReactNum-Adsorption%DissNum)
        reactdesorbnum(Adsorption%DissocReact(1,DissocNum,iSpec)) = &
                reactdesorbnum(Adsorption%DissocReact(1,DissocNum,iSpec)) + 1
        reactdesorbnum(Adsorption%DissocReact(2,DissocNum,iSpec)) = &
                reactdesorbnum(Adsorption%DissocReact(2,DissocNum,iSpec)) + 1
        reactdesorbnum(iSpec) = reactdesorbnum(iSpec) - 1
        IF (WriteMacroValues) THEN
          ! calculate Enthalpie of desorption and sample
          Heat_A = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,iSpec,-1,.FALSE.)
          Heat_Product1 = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Adsorption%DissocReact(1,DissocNum,iSpec),-1,.FALSE.)
          Heat_Product2 = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Adsorption%DissocReact(2,DissocNum,iSpec),-1,.FALSE.)
          D_A = Adsorption%EDissBond(Adsorption%DissNum+iReact,iSpec)
          AdsorptionEnthalpie = ((( Heat_A -Heat_Product1 -Heat_Product2 ) + D_A) * BoltzmannConst &
                          / REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(3))) &
                          * REAL(INT(Adsorption%DensSurfAtoms(SurfSideID) &
                          * SurfMesh%SurfaceArea(subsurfxi,subsurfeta,SurfSideID),8)) / Species(iSpec)%MacroParticleFactor
          !----  Sampling of energies
          IF ((DSMC%CalcSurfaceVal.AND.(Time.ge.(1-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroValues)) THEN
            SampWall(SurfSideID)%Adsorption(1,subsurfxi,subsurfeta) = SampWall(SurfSideID)%Adsorption(1,subsurfxi,subsurfeta) &
                                                                  + AdsorptionEnthalpie * Species(iSpec)%MacroParticleFactor
          END IF
        END IF
        ! remove adsorbate and update map
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(Surfpos) = 0
        DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
          Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(Surfpos,j)
          Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(Surfpos,j)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) = &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) - 1
        END DO
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = &
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1)
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1) = Surfpos
        nSitesRemain(Coord) = nSitesRemain(Coord) + 1
        nSitesRemain(4) = nSitesRemain(4) + 1
        ! add resulting adsorbates from dissociation and update map
        ! first product
        jCoord = Coord_Product(1,iReact)
        DO i = 1,nSites(jCoord)
        IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(i).EQ.Pos_Product(1,iReact)) THEN
          jSpec = Adsorption%DissocReact(1,DissocNum,iSpec)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%Species(Pos_Product(1,iReact)) = jSpec
          DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%nInterAtom
            Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%BondAtomIndx(Pos_Product(1,iReact),j)
            Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%BondAtomIndy(Pos_Product(1,iReact),j)
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(jSpec,Indx,Indy) = &
                SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(jSpec,Indx,Indy) + 1
          END DO
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(i) = &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord))
!           SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord)) = Pos_Product(1,iReact)
          ! move Surfpos to MapID in remainNum segment
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord)) = &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(nSites(jCoord)-remainNum(jCoord))
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(nSites(jCoord)-remainNum(jCoord)) = &
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
        IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(i).EQ.Pos_Product(2,iReact)) THEN
          jSpec = Adsorption%DissocReact(2,DissocNum,iSpec)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%Species(Pos_Product(2,iReact)) = jSpec
          DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%nInterAtom
            Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%BondAtomIndx(Pos_Product(2,iReact),j)
            Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%BondAtomIndy(Pos_Product(2,iReact),j)
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(jSpec,Indx,Indy) = &
                SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(jSpec,Indx,Indy) + 1
          END DO
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(i) = &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord))
!           SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord)) = Pos_Product(2,iReact)
          ! move Surfpos to MapID in remainNum segment
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord)) = &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(nSites(jCoord)-remainNum(jCoord))
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(nSites(jCoord)-remainNum(jCoord)) = &
              Pos_Product(2,iReact)
          remainNum(jCoord) = remainNum(jCoord) + 1
          remainNum(4) = remainNum(4) + 1
          nSitesRemain(jCoord) = nSitesRemain(jCoord) - 1
          nSitesRemain(4) = nSitesRemain(4) - 1
          EXIT
        END IF
        END DO !i
      !-----------------------------------------------------------------------------------------------------------------------------
      CASE(3) ! disproportionation
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
        IF (WriteMacroValues) THEN
          ! calculate Enthalpie of reaction and sample
          Heat_A = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Adsorption%ChemReactant(1,ExchNum),-1,.FALSE.)
          Heat_ReactP = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Adsorption%ChemReactant(2,ExchNum),-1,.FALSE.)
          Heat_Product1 = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Adsorption%ChemProduct(1,ExchNum),-1,.FALSE.)
          Heat_Product2 = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,Adsorption%ChemProduct(2,ExchNum),-1,.FALSE.)
          D_A = Adsorption%Reactant_DissBond_K(1,ExchNum)
          D_ReactP = Adsorption%Reactant_DissBond_K(2,ExchNum)
          D_Product1 = Adsorption%Product_DissBond_K(1,ExchNum)
          D_Product2 = Adsorption%Product_DissBond_K(2,ExchNum)
          AdsorptionEnthalpie = (-(( Heat_A +Heat_ReactP -Heat_Product1 -Heat_Product2 ) &
                          + (D_A +D_ReactP -D_Product1 -D_Product2)) * BoltzmannConst &
                          / REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(3))) &
                          * REAL(INT(Adsorption%DensSurfAtoms(SurfSideID) &
                          * SurfMesh%SurfaceArea(subsurfxi,subsurfeta,SurfSideID),8)) / Species(iSpec)%MacroParticleFactor
          IF (ReactDirForward) AdsorptionEnthalpie = -AdsorptionEnthalpie
          !----  Sampling of energies
          IF ((DSMC%CalcSurfaceVal.AND.(Time.ge.(1-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroValues)) THEN
            SampWall(SurfSideID)%Adsorption(1,subsurfxi,subsurfeta) = SampWall(SurfSideID)%Adsorption(1,subsurfxi,subsurfeta) &
                                                                  + AdsorptionEnthalpie * Species(iSpec)%MacroParticleFactor
          END IF
        END IF
        ! remove first adsorbate and update map
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(Surfpos) = 0
        DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
          Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(Surfpos,j)
          Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(Surfpos,j)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) = &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) - 1
        END DO
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = &
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1)
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1) = Surfpos
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
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%Species(Pos_ReactP(iReact)) = 0
        DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%nInterAtom
          Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%BondAtomIndx(Pos_ReactP(iReact),j)
          Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%BondAtomIndy(Pos_ReactP(iReact),j)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(jSpec,Indx,Indy) = &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(jSpec,Indx,Indy) - 1
        END DO
        DO i = nSitesRemain(jCoord)+1,nSites(jCoord)
        IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(i).EQ.Pos_ReactP(iReact)) THEN
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(i) = SurfDistInfo( &
              subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord)+1)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord) + 1) = Pos_ReactP(iReact)
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
        IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(i).EQ.Pos_Product(1,iReact)) THEN
          IF (ReactDirForward) THEN
            jSpec = Adsorption%ChemProduct(1,ExchNum)
          ELSE
            jSpec = Adsorption%ChemReactant(1,ExchNum)
          END IF
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%Species(Pos_Product(1,iReact)) = jSpec
          DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%nInterAtom
            Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%BondAtomIndx(Pos_Product(1,iReact),j)
            Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%BondAtomIndy(Pos_Product(1,iReact),j)
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(jSpec,Indx,Indy) = &
                SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(jSpec,Indx,Indy) + 1
          END DO
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(i) = &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord))
!           SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord)) = Pos_Product(1,iReact)
          ! move Surfpos to MapID in remainNum segment
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord)) = &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(nSites(jCoord)-remainNum(jCoord))
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(nSites(jCoord)-remainNum(jCoord)) = &
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
        IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(i).EQ.Pos_Product(2,iReact)) THEN
          IF (ReactDirForward) THEN
            jSpec = Adsorption%ChemProduct(2,ExchNum)
          ELSE
            jSpec = Adsorption%ChemReactant(2,ExchNum)
          END IF
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%Species(Pos_Product(2,iReact)) = jSpec
          DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%nInterAtom
            Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%BondAtomIndx(Pos_Product(2,iReact),j)
            Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%BondAtomIndy(Pos_Product(2,iReact),j)
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(jSpec,Indx,Indy) = &
                SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(jSpec,Indx,Indy) + 1
          END DO
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(i) = &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord))
!           SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord)) = Pos_Product(2,iReact)
          ! move Surfpos to MapID in remainNum segment
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(nSitesRemain(jCoord)) = &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(nSites(jCoord)-remainNum(jCoord))
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%UsedSiteMap(nSites(jCoord)-remainNum(jCoord)) = &
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
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(Surfpos) = 0
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(Pos_ReactP(ReactNum_run+1)) = iSpec
        DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
          Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(Surfpos,j)
          Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(Surfpos,j)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) = &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) - 1
          Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(Pos_ReactP(ReactNum_run+1),j)
          Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(Pos_ReactP(ReactNum_run+1),j)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) = &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) + 1
        END DO
        ! move adsorbate to empty site and update map
        DO i = 1,nSitesRemain(Coord)
          IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(i) &
              .EQ.Pos_ReactP(ReactNum_run+1)) THEN
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(i) = Surfpos
!             SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = Pos_ReactP(ReactNum_run+1)
            ! move Surfpos to MapID in remainNum segment
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = &
                SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSites(Coord)-remainNum(Coord))
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSites(Coord)-remainNum(Coord)) &
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
        IF (Adsorption%TPD) Energy(iSpec) = Energy(iSpec) + Heat_A
#endif
        IF (WriteMacroValues) THEN
          ! calculate Enthalpie of desorption and sample
          AdsorptionEnthalpie = (Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,iSpec,-1,.TRUE.) * BoltzmannConst &
                          / REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(3))) &
                          * REAL(INT(Adsorption%DensSurfAtoms(SurfSideID) &
                          * SurfMesh%SurfaceArea(subsurfxi,subsurfeta,SurfSideID),8)) / Species(iSpec)%MacroParticleFactor
          !----  Sampling of energies
          IF ((DSMC%CalcSurfaceVal.AND.(Time.ge.(1-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroValues)) THEN
            SampWall(SurfSideID)%Adsorption(1,subsurfxi,subsurfeta) = SampWall(SurfSideID)%Adsorption(1,subsurfxi,subsurfeta) &
                                                                    + AdsorptionEnthalpie * Species(iSpec)%MacroParticleFactor
          END IF
        END IF
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(Surfpos) = 0
        DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
          Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(Surfpos,j)
          Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(Surfpos,j)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) = &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) - 1
        END DO
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = &
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1)
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1) = Surfpos
        nSitesRemain(Coord) = nSitesRemain(Coord) + 1
        nSitesRemain(4) = nSitesRemain(4) + 1
      END SELECT
      !-----------------------------------------------------------------------------------------------------------------------------
    ELSE ! adsorbate remains on surface on same site
!       SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = &
!           SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSites(Coord)-remainNum(Coord))
!       SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSites(Coord)-remainNum(Coord)) = Surfpos
      remainNum(Coord) = remainNum(Coord) + 1
      remainNum(4) = remainNum(4) + 1
    END IF
    ! update number of adsorbates
    adsorbnum(:) = nSites(:) - nSitesRemain(:)
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(:) = nSitesRemain(:)
    ! increment number of considered surface particles
    trace = trace + 1
  END DO
  END IF
  DO iSpec = 1,nSpecies
    numSites = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(3)
    ! calculate number of desorbed particles for each species
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%desorbnum_tmp(iSpec) = &
                        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%desorbnum_tmp(iSpec) &
                        + ((REAL(desorbnum(iSpec)) / REAL(numSites)) &
                        * REAL(INT(Adsorption%DensSurfAtoms(SurfSideID) &
                        * SurfMesh%SurfaceArea(subsurfxi,subsurfeta,SurfSideID),8)) / Species(iSpec)%MacroParticleFactor)
    ! calculate number of desorbing simulation particles (round to integer)
    Adsorption%SumDesorbPart(subsurfxi,subsurfeta,SurfSideID,iSpec) = &
                        INT(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%desorbnum_tmp(iSpec))
    ! Adjust tracking desorbing simulation particles
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%desorbnum_tmp(iSpec) = &
                        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%desorbnum_tmp(iSpec) &
                        - Adsorption%SumDesorbPart(subsurfxi,subsurfeta,SurfSideID,iSpec)
    ! calculate number of reacted particles for each species
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%reactnum_tmp(iSpec) = &
                        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%reactnum_tmp(iSpec) &
                        + ((REAL(reactdesorbnum(iSpec)) / REAL(numSites)) &
                        * REAL(INT(Adsorption%DensSurfAtoms(SurfSideID) &
                        * SurfMesh%SurfaceArea(subsurfxi,subsurfeta,SurfSideID),8)) / Species(iSpec)%MacroParticleFactor)
    ! calculate number of reacting simulation particles on surface (round to integer)
    Adsorption%SumReactPart(subsurfxi,subsurfeta,SurfSideID,iSpec) = &
                        INT(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%reactnum_tmp(iSpec))
    ! Adjust tracking reacting simulation particles
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%reactnum_tmp(iSpec) = &
                        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%reactnum_tmp(iSpec) &
                        - Adsorption%SumReactPart(subsurfxi,subsurfeta,SurfSideID,iSpec)
    ! Sample vibrational energies 
    ! due to reaction a part of energies can be transformed into other vibrational groundstates and the change is sampled 
    ! the energies of the emitted particles are sampled in surfflux routine (particle_emission.f90) because calculated there
    IF ((DSMC%CalcSurfaceVal.AND.(Time.ge.(1-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroValues)) THEN
      IF (SpecDSMC(iSpec)%InterID.EQ.2) THEN
        SampleParts = Adsorption%SumDesorbPart(subsurfxi,subsurfeta,SurfSideID,iSpec) &
                    - Adsorption%SumReactPart(subsurfxi,subsurfeta,SurfSideID,iSpec)
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
          EvibNew = 0. ! later added in particle emission when particle ist emitted
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
            SampWall(SurfSideID)%State(7,subsurfxi,subsurfeta) = SampWall(SurfSideID)%State(7,subsurfxi,subsurfeta) &
                                              + EvibOld * Species(iSpec)%MacroParticleFactor
            SampWall(SurfSideID)%State(8,subsurfxi,subsurfeta) = SampWall(SurfSideID)%State(8,subsurfxi,subsurfeta) &
                                              + EvibWall * Species(iSpec)%MacroParticleFactor
            SampWall(SurfSideID)%State(9,subsurfxi,subsurfeta) = SampWall(SurfSideID)%State(9,subsurfxi,subsurfeta) &
                                              + EvibNew * Species(iSpec)%MacroParticleFactor
          END DO
        END IF
        SDEALLOCATE(RanNumPoly)
        SDEALLOCATE(VibQuantWallPoly)
      END IF
    END IF
    ! analyze rate data
    IF (adsorbates(iSpec).EQ.0) THEN
      Adsorption%ProbDes(subsurfxi,subsurfeta,SurfSideID,iSpec) = 0
    ELSE
      Adsorption%ProbDes(subsurfxi,subsurfeta,SurfSideID,iSpec) = P_des(iSpec) / REAL(adsorbates(iSpec))
    END IF
#if (PP_TimeDiscMethod==42)
    IF (Adsorption%TPD) THEN
      IF ((desorbnum(iSpec)-reactdesorbnum(iSpec)).LE.0) THEN
        Energy(iSpec) = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,iSpec,-1,.FALSE.)
      ELSE
        Energy(iSpec) = Energy(iSpec) /REAL(desorbnum(iSpec)-reactdesorbnum(iSpec))
      END IF
    END IF
    Adsorption%AdsorpInfo(iSpec)%NumOfDes = Adsorption%AdsorpInfo(iSpec)%NumOfDes &
                                          + Adsorption%SumDesorbPart(subsurfxi,subsurfeta,SurfSideID,iSpec)
    Adsorption%AdsorpInfo(iSpec)%MeanProbDes = Adsorption%AdsorpInfo(iSpec)%MeanProbDes &
                                              + Adsorption%ProbDes(subsurfxi,subsurfeta,SurfSideID,iSpec)
    IF (Adsorption%TPD) Adsorption%AdsorpInfo(iSpec)%MeanEAds = Adsorption%AdsorpInfo(iSpec)%MeanEAds + Energy(iSpec)
#endif
  END DO
END DO
END DO
END DO
#if (PP_TimeDiscMethod==42)
IF (.NOT.DSMC%ReservoirRateStatistic) THEN
  Adsorption%AdsorpInfo(:)%MeanProbDes = Adsorption%AdsorpInfo(:)%MeanProbDes / (nSurfSample * nSurfSample * SurfMesh%nSides)
  IF (Adsorption%TPD) Adsorption%AdsorpInfo(:)%MeanEAds = Adsorption%AdsorpInfo(:)%MeanEAds &
                                                        / (nSurfSample * nSurfSample * SurfMesh%nSides)
END IF
IF (Adsorption%TPD) DEALLOCATE(Energy)
#endif

DEALLOCATE(desorbnum,adsorbnum,nSites,nSitesRemain,remainNum,P_des,adsorbates)
DEALLOCATE(P_react,P_actual_react)
DEALLOCATE(Coord_ReactP,Pos_ReactP)

END SUBROUTINE CalcBackgndPartDesorb

SUBROUTINE AdjustBackgndAdsNum(subsurfxi,subsurfeta,SurfSideID,adsorbates_num,iSpec)
!===================================================================================================================================
! Routine for adjusting the number of Adsorbates for background surface, if Adsorption took place in CalcBackgndPartAdsorb 
! and Coverage changed strong enough
!===================================================================================================================================
  USE MOD_Globals_Vars,           ONLY : PlanckConst
  USE MOD_Globals,                ONLY : abort
  USE MOD_Particle_Vars,          ONLY : BoltzmannConst
  USE MOD_DSMC_Vars,              ONLY : Adsorption, SurfDistInfo
!===================================================================================================================================
   IMPLICIT NONE
!===================================================================================================================================
! Local variable declaration
  INTEGER,INTENT(IN)               :: subsurfxi,subsurfeta,SurfSideID,adsorbates_num,iSpec
! LOCAL VARIABLES
  INTEGER                          :: dist
  INTEGER                          :: Coord, Surfnum, Surfpos, UsedSiteMapPos, nSites, nSitesRemain
  INTEGER                          :: iInterAtom, xpos, ypos
  REAL                             :: RanNum
  
!===================================================================================================================================
  IF (adsorbates_num.GT.0) THEN
    ! distribute adsorbates randomly on the surface on the correct site and assign surface atom bond order
    dist = 1
    Coord = Adsorption%Coordination(iSpec)
    Surfnum = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord)
    IF ((SurfNum - adsorbates_num).LT.0) THEN
      CALL abort(&
__STAMP__&
,'Error in AdjustBackgndAdsNum: Too many new Adsorbates! not enough Sites for Coordination:',Adsorption%Coordination(iSpec))
    END IF
    DO WHILE (dist.LE.adsorbates_num) 
      CALL RANDOM_NUMBER(RanNum)
      Surfpos = 1 + INT(Surfnum * RanNum)
      UsedSiteMapPos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(Surfpos)
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(UsedSiteMapPos) = iSpec
      ! assign bond order of respective surface atoms in the surfacelattice
      DO iInterAtom = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
        xpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(UsedSiteMapPos,iInterAtom)
        ypos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(UsedSiteMapPos,iInterAtom)
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xpos,ypos) = &
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xpos,ypos) + 1
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
    Coord = Adsorption%Coordination(iSpec)
    nSites = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Coord)
    nSitesRemain = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord)
    Surfnum = nSites - nSitesRemain
    DO WHILE (dist.GE.adsorbates_num) 
      CALL RANDOM_NUMBER(RanNum)
      Surfpos = 1 + INT(Surfnum * RanNum)
      UsedSiteMapPos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain+Surfpos)
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(UsedSiteMapPos) = 0
      ! assign bond order of respective surface atoms in the surfacelattice
      DO iInterAtom = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
        xpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(UsedSiteMapPos,iInterAtom)
        ypos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(UsedSiteMapPos,iInterAtom)
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xpos,ypos) = &
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xpos,ypos) - 1
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
    
END SUBROUTINE AdjustBackgndAdsNum

#ifdef MPI
SUBROUTINE ExchangeAdsorbNum() 
!===================================================================================================================================
! exchange the number of adsorbing particles on halo surface 
! only processes with samling sides in their halo region and the original process participate on the communication
! structure is similar to surface sampling/particle communication
! each process sends his halo-information directly to the origin process by use of a list, containing the surfsideids for sending
! the receiving process adds the new data to his own sides
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Vars               ,ONLY:nSpecies
USE MOD_DSMC_Vars                   ,ONLY:Adsorption
USE MOD_Particle_Boundary_Vars      ,ONLY:SurfComm,nSurfSample
USE MOD_Particle_MPI_Vars           ,ONLY:AdsorbSendBuf,AdsorbRecvBuf,SurfExchange
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: MessageSize,nValues,iSurfSide,SurfSideID
INTEGER                         :: iPos,p,q,iProc
INTEGER                         :: recv_status_list(1:MPI_STATUS_SIZE,1:SurfCOMM%nMPINeighbors)
!===================================================================================================================================

nValues=nSpecies*(nSurfSample)**2

! open receive buffer
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesRecv(iProc).EQ.0) CYCLE
  MessageSize=SurfExchange%nSidesRecv(iProc)*nValues
  CALL MPI_IRECV( AdsorbRecvBuf(iProc)%content_int             &
                , MessageSize                                  &
                , MPI_INT                                      &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID     &
                , 1010                                         &
                , SurfCOMM%COMM                                &
                , SurfExchange%RecvRequest(iProc)              & 
                , IERROR )
END DO ! iProc

! build message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesSend(iProc).EQ.0) CYCLE
  iPos=0
  DO iSurfSide=1,SurfExchange%nSidesSend(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%SendList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        AdsorbSendBuf(iProc)%content_int(iPos+1:iPos+nSpecies)= Adsorption%SumAdsorbPart(p,q,SurfSideID,:)
        iPos=iPos+nSpecies
      END DO ! p=0,nSurfSample
    END DO ! q=0,nSurfSample
!     Adsorption%SumAdsorbPart(:,:,SurfSideID,:)=0.
  END DO ! iSurfSide=1,nSurfExchange%nSidesSend(iProc)
END DO

! send message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesSend(iProc).EQ.0) CYCLE
  MessageSize=SurfExchange%nSidesSend(iProc)*nValues
  CALL MPI_ISEND( AdsorbSendBuf(iProc)%content_int         &
                , MessageSize                              & 
                , MPI_INT                                  &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID & 
                , 1010                                     &
                , SurfCOMM%COMM                            &   
                , SurfExchange%SendRequest(iProc)          &
                , IERROR )                                     
END DO ! iProc                                                

! 4) Finish Received number of particles
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesSend(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfExchange%SendRequest(iProc),MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error', IERROR)
  END IF
  IF(SurfExchange%nSidesRecv(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfExchange%RecvRequest(iProc),recv_status_list(:,iProc),IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error', IERROR)
  END IF
END DO ! iProc

! add data do my list
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesRecv(iProc).EQ.0) CYCLE
  MessageSize=SurfExchange%nSidesSend(iProc)*nValues
  iPos=0
  DO iSurfSide=1,SurfExchange%nSidesRecv(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%RecvList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        Adsorption%SumAdsorbPart(p,q,SurfSideID,:)=Adsorption%SumAdsorbPart(p,q,SurfSideID,:) &
                                         +AdsorbRecvBuf(iProc)%content_int(iPos+1:iPos+nSpecies)
        iPos=iPos+nSpecies
      END DO ! p=0,nSurfSample
    END DO ! q=0,nSurfSample
  END DO ! iSurfSide=1,nSurfExchange%nSidesSend(iProc)
  AdsorbRecvBuf(iProc)%content_int = 0
  AdsorbSendBuf(iProc)%content_int = 0
END DO ! iProc

END SUBROUTINE ExchangeAdsorbNum

SUBROUTINE ExchangeSurfDistSize()
!===================================================================================================================================
! exchange the number of surface distribution sites to communicate to halosides with neighbours
! only processes with surface sides in their halo region and the original process participate on the communication
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Boundary_Vars      ,ONLY:nSurfSample,SurfComm
USE MOD_Particle_MPI_Vars           ,ONLY:SurfDistSendBuf,SurfDistRecvBuf,SurfExchange
USE MOD_DSMC_Vars                   ,ONLY:SurfDistInfo
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: MessageSize,iSurfSide,SurfSideID
INTEGER                         :: p,q,iProc
INTEGER                         :: recv_status_list(1:MPI_STATUS_SIZE,1:SurfCOMM%nMPINeighbors)
INTEGER                         :: iCoord
!===================================================================================================================================

! open receive buffer
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSurfDistSidesRecv(iProc).EQ.0) CYCLE
  MessageSize=SurfExchange%nSurfDistSidesRecv(iProc)
  CALL MPI_IRECV( SurfExchange%NbrOfPos(iProc)%nPosRecv(:)     &
                , MessageSize                                  &
                , MPI_INT                                      &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID     &
                , 1011                                         &
                , SurfCOMM%COMM                                &
                , SurfExchange%SurfDistRecvRequest(1,iProc)    & 
                , IERROR )
END DO ! iProc

DO iProc=1,SurfCOMM%nMPINeighbors
DO iSurfSide=1,SurfExchange%nSurfDistSidesSend(iProc)
SurfSideID=SurfCOMM%MPINeighbor(iProc)%SurfDistSendList(iSurfSide)
DO q=1,nSurfSample
DO p=1,nSurfSample
DO iCoord = 1,3
    SurfExchange%NbrOfPos(iProc)%nPosSend(iSurfSide) = SurfExchange%NbrOfPos(iProc)%nPosSend(iSurfSide) &
                                                     + 3 + 2*SurfDistInfo(p,q,SurfSideID)%nSites(iCoord)
END DO
END DO
END DO
END DO
END DO

! send message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSurfDistSidesSend(iProc).EQ.0) CYCLE
  MessageSize=SurfExchange%nSurfDistSidesSend(iProc)
  CALL MPI_ISEND( SurfExchange%NbrOfPos(iProc)%nPosSend(:)     &
                , MessageSize                                  & 
                , MPI_INT                                      &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID     & 
                , 1011                                         &
                , SurfCOMM%COMM                                &   
                , SurfExchange%SurfDistSendRequest(1,iProc)    &
                , IERROR )
END DO ! iProc

! 4) Finish receiving commsizes
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSurfDistSidesSend(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfExchange%SurfDistSendRequest(1,iProc),MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error in number of surface distribution sites (send)', IERROR)
  END IF
  IF(SurfExchange%nSurfDistSidesRecv(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfExchange%SurfDistRecvRequest(1,iProc),recv_status_list(:,iProc),IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error in number of surface distribution sites (receive)', IERROR)
  END IF
END DO ! iProc

DO iProc=1,SurfCOMM%nMPINeighbors
ALLOCATE(SurfDistSendBuf(iProc)%content_int(1:SUM(SurfExchange%NbrOfPos(iProc)%nPosSend(:))))
ALLOCATE(SurfDistRecvBuf(iProc)%content_int(1:SUM(SurfExchange%NbrOfPos(iProc)%nPosRecv(:))))
END DO

END SUBROUTINE ExchangeSurfDistSize


SUBROUTINE ExchangeSurfDistInfo() 
!===================================================================================================================================
! exchange the surface distribution to halosides of neighbours
! only processes with surface sides in their halo region and the original process participate on the communication
! structure is similar to surface sampling/particle communication
! each process sends his halo-information directly to the origin process by use of a list, containing the surfsideids for sending
! the receiving process adds the new data to his own sides
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Boundary_Vars      ,ONLY:SurfMesh,SurfComm,nSurfSample
USE MOD_Particle_MPI_Vars           ,ONLY:SurfDistSendBuf,SurfDistRecvBuf,SurfExchange
USE MOD_DSMC_Vars                   ,ONLY:SurfDistInfo
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: MessageSize, iSurfSide, SurfSideID, iSurf
INTEGER                         :: iPos,p,q,iProc
INTEGER                         :: recv_status_list(1:MPI_STATUS_SIZE,1:SurfCOMM%nMPINeighbors)
INTEGER                         :: iCoord,nSites,nSitesRemain,iSite,iInteratom,UsedSiteMapPos,iSpec,xpos,ypos
!===================================================================================================================================

! open receive buffer
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSurfDistSidesRecv(iProc).EQ.0) CYCLE
  Messagesize = 0
  DO iSurf = 1,SurfExchange%nSurfDistSidesRecv(iProc)
    MessageSize = MessageSize + SurfExchange%NbrOfPos(iProc)%nPosRecv(iSurf)
  END DO
  CALL MPI_IRECV( SurfDistRecvBuf(iProc)%content_int           &
                , MessageSize                                  &
                , MPI_INT                                      &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID     &
                , 1012                                         &
                , SurfCOMM%COMM                                &
                , SurfExchange%SurfDistRecvRequest(2,iProc)      & 
                , IERROR )
END DO ! iProc

! build message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSurfDistSidesSend(iProc).EQ.0) CYCLE
  iPos=0
  DO iSurfSide=1,SurfExchange%nSurfDistSidesSend(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%SurfDistSendList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        DO iCoord = 1,3
          SurfDistSendBuf(iProc)%content_int(iPos+1) = SurfDistInfo(p,q,SurfSideID)%SitesRemain(iCoord)
          iPos=iPos+1
          SurfDistSendBuf(iProc)%content_int(iPos+1:iPos+SurfDistInfo(p,q,SurfSideID)%nSites(iCoord)) = &
              SurfDistInfo(p,q,SurfSideID)%AdsMap(iCoord)%Species(:)
          iPos=iPos+SurfDistInfo(p,q,SurfSideID)%nSites(iCoord)
          SurfDistSendBuf(iProc)%content_int(iPos+1:iPos+SurfDistInfo(p,q,SurfSideID)%nSites(iCoord)) = &
              SurfDistInfo(p,q,SurfSideID)%AdsMap(iCoord)%UsedSiteMap(:)
          iPos=iPos+SurfDistInfo(p,q,SurfSideID)%nSites(iCoord)
        END DO
      END DO ! p=0,nSurfSample
    END DO ! q=0,nSurfSample
  END DO ! iSurfSide=1,nSurfExchange%nSurfDistSidesSend(iProc)
END DO

! send message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSurfDistSidesSend(iProc).EQ.0) CYCLE
  Messagesize = 0
  DO iSurf = 1,SurfExchange%nSurfDistSidesSend(iProc)
    MessageSize = MessageSize + SurfExchange%NbrOfPos(iProc)%nPosSend(iSurf)
  END DO
  CALL MPI_ISEND( SurfDistSendBuf(iProc)%content_int       &
                , MessageSize                              & 
                , MPI_INT                                  &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID & 
                , 1012                                     &
                , SurfCOMM%COMM                            &   
                , SurfExchange%SurfDistSendRequest(2,iProc)  &
                , IERROR )
END DO ! iProc

! 4) Finish received surface distribution
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSurfDistSidesSend(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfExchange%SurfDistSendRequest(2,iProc),MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error in surface distribution (send)', IERROR)
  END IF
  IF(SurfExchange%nSurfDistSidesRecv(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfExchange%SurfDistRecvRequest(2,iProc),recv_status_list(:,iProc),IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error in Surface distribution (receive)', IERROR)
  END IF
END DO ! iProc

! add data do my list
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSurfDistSidesRecv(iProc).EQ.0) CYCLE
  iPos=0
  DO iSurfSide=1,SurfExchange%nSurfDistSidesRecv(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%SurfDistRecvList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        DO iCoord = 1,3
          SurfDistInfo(p,q,SurfSideID)%SitesRemain(iCoord) = SurfDistRecvBuf(iProc)%content_int(iPos+1)
          iPos=iPos+1
          SurfDistInfo(p,q,SurfSideID)%AdsMap(iCoord)%Species(:) = &
              SurfDistRecvBuf(iProc)%content_int(iPos+1:iPos+SurfDistInfo(p,q,SurfSideID)%nSites(iCoord))
          iPos=iPos+SurfDistInfo(p,q,SurfSideID)%nSites(iCoord)
          SurfDistInfo(p,q,SurfSideID)%AdsMap(iCoord)%UsedSiteMap(:) = &
              SurfDistRecvBuf(iProc)%content_int(iPos+1:iPos+SurfDistInfo(p,q,SurfSideID)%nSites(iCoord))
          iPos=iPos+SurfDistInfo(p,q,SurfSideID)%nSites(iCoord)
        END DO
      END DO ! p=0,nSurfSample
    END DO ! q=0,nSurfSample.
  END DO ! iSurfSide=1,nSurfExchange%nSidesSend(iProc)
  SurfDistRecvBuf(iProc)%content_int = 0
  SurfDistSendBuf(iProc)%content_int = 0
END DO ! iProc

! assign bond order to surface atoms in the surfacelattice for halo sides
DO iSurfSide = SurfMesh%nSides+1,SurfMesh%nTotalSides
  DO q=1,nSurfSample
    DO p=1,nSurfSample
      SurfDistInfo(p,q,iSurfSide)%SurfAtomBondOrder(:,:,:) = 0
      DO iCoord = 1,3
        nSitesRemain = SurfDistInfo(p,q,iSurfSide)%SitesRemain(iCoord)
        nSites = SurfDistInfo(p,q,iSurfSide)%nSites(iCoord)
        DO iSite = nSitesRemain+1,nSites
          UsedSiteMapPos = SurfDistInfo(p,q,iSurfSide)%AdsMap(iCoord)%UsedSiteMap(iSite)
          iSpec = SurfDistInfo(p,q,iSurfSide)%AdsMap(iCoord)%Species(UsedSiteMapPos)
          DO iInterAtom = 1,SurfDistInfo(p,q,iSurfSide)%AdsMap(iCoord)%nInterAtom
            xpos = SurfDistInfo(p,q,iSurfSide)%AdsMap(iCoord)%BondAtomIndx(UsedSiteMapPos,iInterAtom)
            ypos = SurfDistInfo(p,q,iSurfSide)%AdsMap(iCoord)%BondAtomIndy(UsedSiteMapPos,iInterAtom)
            SurfDistInfo(p,q,iSurfSide)%SurfAtomBondOrder(iSpec,xpos,ypos) = &
              SurfDistInfo(p,q,iSurfSide)%SurfAtomBondOrder(iSpec,xpos,ypos) + 1
          END DO
        END DO
      END DO
    END DO
  END DO
END DO

END SUBROUTINE ExchangeSurfDistInfo
#endif /*MPI*/


SUBROUTINE CalcDistNumChange()
!===================================================================================================================================
! updates mapping of adsorbed particles if particles desorbed (sumdesorbpart.GT.0)
!===================================================================================================================================
  USE MOD_Particle_Vars,          ONLY : nSpecies
  USE MOD_DSMC_Vars,              ONLY : Adsorption, SurfDistInfo
  USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh
!===================================================================================================================================
   IMPLICIT NONE
!===================================================================================================================================
! Local variable declaration
   INTEGER                          :: SurfSideID, iSpec, subsurfeta, subsurfxi, Adsorbates
   INTEGER                          :: Coord, nSites, nSitesRemain, iInterAtom, xpos, ypos
   REAL                             :: RanNum
   INTEGER                          :: Surfpos, UsedSiteMapPos, dist, Surfnum, newAdsorbates
!===================================================================================================================================

DO SurfSideID = 1,SurfMesh%nSides
DO subsurfeta = 1,nSurfSample
DO subsurfxi = 1,nSurfSample

  DO iSpec = 1,nSpecies
    Coord = Adsorption%Coordination(iSpec)
    nSites = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Coord)
    nSitesRemain = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord)
    Adsorbates = INT(Adsorption%Coverage(subsurfxi,subsurfeta,SurfSideID,iSpec) &
                * SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Adsorption%Coordination(iSpec)))
    IF (nSites-nSitesRemain.GT.Adsorbates) THEN
      newAdsorbates = (nSites-nSitesRemain) - Adsorbates
      ! remove adsorbates randomly from the surface from the correct site and assign surface atom bond order
      dist = 1
      Surfnum = nSites - nSitesRemain
      DO WHILE (dist.LE.newAdsorbates) 
        CALL RANDOM_NUMBER(RanNum)
        Surfpos = 1 + INT(Surfnum * RanNum)
        UsedSiteMapPos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain+Surfpos)
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(UsedSiteMapPos) = 0
        ! assign bond order of respective surface atoms in the surfacelattice
        DO iInterAtom = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
          xpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(UsedSiteMapPos,iInterAtom)
          ypos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(UsedSiteMapPos,iInterAtom)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xpos,ypos) = &
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xpos,ypos) - 1
        END DO
        ! rearrange UsedSiteMap-Surfpos-array
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain+Surfpos) = &
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain+1)
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain+1) = UsedSiteMapPos
        Surfnum = Surfnum - 1
        nSitesRemain = nSitesRemain + 1
        dist = dist + 1
      END DO
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord) = nSitesRemain
    ELSE IF (nSites-nSitesRemain.LT.Adsorbates) THEN
      newAdsorbates = Adsorbates - (nSites-nSitesRemain)
      ! add new Adsorbates to macro adsorbate map
      dist = 1
      Surfnum = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord)
      DO WHILE (dist.LE.newAdsorbates) 
        CALL RANDOM_NUMBER(RanNum)
        Surfpos = 1 + INT(Surfnum * RanNum)
        UsedSiteMapPos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(Surfpos)
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(UsedSiteMapPos) = iSpec
        ! assign bond order of respective surface atoms in the surfacelattice
        DO iInterAtom = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
          xpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(UsedSiteMapPos,iInterAtom)
          ypos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(UsedSiteMapPos,iInterAtom)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xpos,ypos) = &
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xpos,ypos) + 1
        END DO
        ! rearrange UsedSiteMap-Surfpos-array
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(Surfpos) = &
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(Surfnum)
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(Surfnum) = UsedSiteMapPos
        Surfnum = Surfnum - 1
        dist = dist + 1
      END DO
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord) = Surfnum
    ELSE
      CYCLE
    END IF
  END DO

END DO
END DO
END DO

END SUBROUTINE CalcDistNumChange

SUBROUTINE CalcSurfDistInteraction()
!===================================================================================================================================
! Model for calculating surface distibution and effects of coverage on heat of adsorption
!===================================================================================================================================
  USE MOD_Particle_Vars,          ONLY : nSpecies
  USE MOD_DSMC_Vars,              ONLY : Adsorption, SurfDistInfo
  USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh
!===================================================================================================================================
  IMPLICIT NONE
!=================================================================================================================================== 
! Local variable declaration
  INTEGER                          :: SurfSideID, subsurfxi, subsurfeta, iSpec, counter
  REAL                             :: sigmasub 
  REAL, ALLOCATABLE                :: sigma(:,:,:,:), ProbSigma(:,:,:,:)
  INTEGER                          :: firstval, secondval, thirdval, fourthval, pos
  INTEGER                          :: Surfpos, nSites, nSitesRemain, xpos, ypos, AdsorbID, iInterAtom, bondorder, Coord
!===================================================================================================================================
ALLOCATE( sigma(1:36*nSpecies,1:36*nSpecies,1:36*nSpecies,1:36*nSpecies),&
          ProbSigma(1:36*nSpecies,1:36*nSpecies,1:36*nSpecies,1:36*nSpecies) )

DO SurfSideID=1,SurfMesh%nSides
DO subsurfeta = 1,nSurfSample
DO subsurfxi = 1,nSurfSample

  DO iSpec=1,nSpecies
    Coord = Adsorption%Coordination(iSpec)
    nSites = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Coord)
    nSitesRemain = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord)
    
    pos = 1
    sigma = 1.
    ProbSigma = 0.
    firstval = 1
    secondval = 1
    thirdval = 1
    fourthval = 1
    counter = 0
      
    DO AdsorbID = nSitesRemain+1,nSites,1
      Surfpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID)

      sigmasub = 0.
      DO iInterAtom = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
        xpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(Surfpos,iInterAtom)
        ypos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(Surfpos,iInterAtom)
        bondorder = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xpos,ypos)
        sigmasub = sigmasub + (1./SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom &
                  * 1./bondorder * (2.-1./bondorder))
        IF (bondorder.GE.firstval) THEN
          fourthval = thirdval
          thirdval = secondval
          secondval = firstval
          firstval = bondorder
        ELSE IF (bondorder.GE.secondval) THEN
          fourthval = thirdval
          thirdval = secondval
          secondval = bondorder
        ELSE IF (bondorder.GE.thirdval) THEN
          fourthval = thirdval
          thirdval = bondorder
        ELSE
          fourthval = bondorder
        END IF
      END DO
      sigma(firstval,secondval,thirdval,fourthval) = sigmasub
      ProbSigma(firstval,secondval,thirdval,fourthval) = ProbSigma(firstval,secondval,thirdval,fourthval) + 1.
      counter = counter + 1
    END DO
    
    IF (counter.GT.0) THEN
      firstval = 1
      secondval = 1
      thirdval = 1
      fourthval = 1
      DO WHILE ((firstval.LE.4).AND.(secondval.LE.4).AND.(thirdval.LE.4).AND.(fourthval.LE.4))
        DO firstval = 1,4
        DO secondval = 1,firstval
        DO thirdval = 1,secondval
        DO fourthval = 1,thirdval
          Adsorption%Sigma(subsurfxi,subsurfeta,SurfSideID,iSpec,pos) = sigma(firstval,secondval,thirdval,fourthval)
          Adsorption%ProbSigma(subsurfxi,subsurfeta,SurfSideID,iSpec,pos) = ProbSigma(firstval,secondval,thirdval,fourthval)/counter
          pos = pos + 1
        END DO
        END DO
        END DO
        END DO
      END DO
    ELSE
      Adsorption%Sigma(subsurfxi,subsurfeta,SurfSideID,iSpec,(1+(iSpec-1)*36):(36*iSpec)) = 1.
      Adsorption%ProbSigma(subsurfxi,subsurfeta,SurfSideID,iSpec,(1+(iSpec-1)*36):(36*iSpec)) = 0.
    END IF
    
  END DO! iSpec = 1,nSpecies
    
END DO
END DO
END DO

DEALLOCATE(sigma,ProbSigma)

END SUBROUTINE CalcSurfDistInteraction


SUBROUTINE CalcAdsorbProb()
!===================================================================================================================================
! Models for adsorption probability calculation
!===================================================================================================================================
  USE MOD_Particle_Vars,          ONLY : nSpecies
  USE MOD_DSMC_Vars,              ONLY : Adsorption, DSMC
  USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh
!===================================================================================================================================
  IMPLICIT NONE
!===================================================================================================================================
! Local variable declaration
  INTEGER                          :: SurfSide, iSpec, p, q
  REAL                             :: Theta_req, Kfactor, S_0
!=================================================================================================================================== 
  DO iSpec=1,nSpecies
  DO SurfSide=1,SurfMesh%nSides
  DO q = 1,nSurfSample
  DO p = 1,nSurfSample
!===================================================================================================================================  
  IF (DSMC%WallModel.EQ.1) THEN
!===================================================================================================================================
! Kisluik Sticking Model from Kolasinski's Surface Science (book)
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
!===================================================================================================================================    
  ELSE IF (DSMC%WallModel.EQ.2) THEN
!=================================================================================================================================== 
! transition state theory(TST) and unity bond index-quadratic exponential potential (UBI-QEP) 

!     Theta = Adsorption%Coverage(p,q,SurfSide,iSpec)  
!     IF (Theta.GT.0) THEN
!       Q_0A = 26312.
!       D_A  = 59870.
!       n = 2
!       m = Theta * n
!       IF (m.LT.1) THEN
!   !       Q = (9*Q_0A**2)/(6*Q_0A+16*D_A)
!         Q = (Q_0A**2)/(Q_0A/n+D_A)
!       ELSE 
!         sigma = 1/m*(2-1/m)
!   !       Q = (9*Q_0A**2)/(6*Q_0A+16*D_A) *sigma
!         Q = (Q_0A**2)/(Q_0A/n+D_A) * sigma
!       END IF
!       
!       CALL PartitionFuncAct(iSpec, WallTemp, VarPartitionFuncAct)
!       CALL PartitionFuncGas(iSpec, WallTemp, VarPartitionFuncGas)
!       
!       E_des = Q
!       nu_des = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct / VarPartitionFuncWall)
!       rate = nu_des * exp(-E_des/WallTemp)
!       Adsorption%AdsorpInfo(iSpec)%ProbDes(p,q,SurfSide) = rate * dt
!     ELSE
!       Adsorption%AdsorpInfo(iSpec)%ProbDes(p,q,SurfSide) = 0.0
!     END IF
!===================================================================================================================================
  END IF
#if (PP_TimeDiscMethod==42)
    IF (.NOT.DSMC%ReservoirRateStatistic) THEN
      Adsorption%AdsorpInfo(iSpec)%MeanProbAds = Adsorption%AdsorpInfo(iSpec)%MeanProbAds + Adsorption%ProbAds(p,q,SurfSide,iSpec)
    END IF
#endif
  END DO
  END DO
  END DO
  END DO
#if (PP_TimeDiscMethod==42)
IF (.NOT.DSMC%ReservoirRateStatistic) THEN
  Adsorption%AdsorpInfo(:)%MeanProbAds = Adsorption%AdsorpInfo(:)%MeanProbAds / (nSurfSample * nSurfSample * SurfMesh%nSides)
END IF
#endif
END SUBROUTINE CalcAdsorbProb

SUBROUTINE CalcDesorbProb()
!===================================================================================================================================
! Models for desorption probability calculation
!===================================================================================================================================
  USE MOD_Globals_Vars,           ONLY : PlanckConst
  USE MOD_Particle_Vars,          ONLY : nSpecies, BoltzmannConst
  USE MOD_Mesh_Vars,              ONLY : BC
  USE MOD_DSMC_Vars,              ONLY : Adsorption, DSMC
  USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh, PartBound
  USE MOD_TimeDisc_Vars,          ONLY : dt
#if (PP_TimeDiscMethod==42)  
  USE MOD_TimeDisc_Vars,          ONLY : iter
  USE MOD_DSMC_Vars,              ONLY : DSMC
#endif
!===================================================================================================================================
   IMPLICIT NONE
!===================================================================================================================================
! Local variable declaration
   INTEGER                          :: SurfSide, iSpec, globSide, p, q
   REAL                             :: Theta, nu_des, rate, WallTemp
   REAL                             :: Q_0A, D_A, Heat, E_des, sigma
!    REAL                :: sigma(10)
   INTEGER                          :: n, iProbSigma
   REAL                             :: VarPartitionFuncAct, VarPartitionFuncWall
!===================================================================================================================================
! CALL CalcSurfDistInteraction()
DO SurfSide=1,SurfMesh%nSides
  globSide = Adsorption%SurfSideToGlobSideMap(SurfSide)
! special TPD (temperature programmed desorption) temperature adjustment routine    
#if (PP_TimeDiscMethod==42)
  IF (Adsorption%TPD) THEN
    WallTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC(globSide))) + (Adsorption%TPD_beta * iter * dt)
    Adsorption%TPD_Temp = Walltemp
  ELSE
    WallTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC(globSide)))
  END IF
#else
  WallTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC(globSide)))
#endif

  DO iSpec = 1,nSpecies 
  DO q = 1,nSurfSample
  DO p = 1,nSurfSample

!===================================================================================================================================
  IF (DSMC%WallModel.EQ.1) THEN
!===================================================================================================================================
! Polanyi-Wigner-eq. from Kolasinski's Surface Science (book) 
!===================================================================================================================================
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
      Adsorption%AdsorpInfo(iSpec)%MeanProbDes = Adsorption%AdsorpInfo(iSpec)%MeanProbDes + Adsorption%ProbDes(p,q,SurfSide,iSpec)
    END IF
#endif
!===================================================================================================================================
  ELSE IF (DSMC%WallModel.EQ.2) THEN
!===================================================================================================================================
! transition state theory(TST) and unity bond index-quadratic exponential potential (UBI-QEP) --> lattice gas model 
    Theta = Adsorption%Coverage(p,q,SurfSide,iSpec)  
    IF (Theta.GT.0) THEN
      Q_0A = Adsorption%HeatOfAdsZero(iSpec)
      D_A  = 59870.
      
      rate = 0.
!       iProbSigma = 1
      DO iProbSigma = (1+(iSpec-1)*36),(36*iSpec)
        SELECT CASE(Adsorption%Coordination(iSpec))
        CASE(1)
          n = 4
          sigma = Adsorption%Sigma(p,q,SurfSide,iSpec,iProbSigma)
          Heat = (Q_0A**2)/(Q_0A/n+D_A) * sigma
        CASE(2)
          n = 2
!         m = Theta * n
!         IF (m.LT.1) THEN
!           Heat = (9*Q_0A**2)/(6*Q_0A+16*D_A)
! !           Heat = (Q_0A**2)/(Q_0A/n+D_A)
!         ELSE 
!           sigma = 1/m*(2-1/m)
          sigma = Adsorption%Sigma(p,q,SurfSide,iSpec,iProbSigma)
          Heat = (9*Q_0A**2)/(6*Q_0A+16*D_A) * sigma
  !         Heat = (Q_0A**2)/(Q_0A/n+D_A) * sigma
!         END IF
        CASE DEFAULT
        END SELECT
        E_des = Heat
!         rate = exp(-E_des/WallTemp)
!         rate = rate + Adsorption%ProbSigma(p,q,SurfSide,iSpec,iProbSigma) * exp(-E_des/WallTemp)
!       END DO
      
        CALL PartitionFuncAct(iSpec, WallTemp, VarPartitionFuncAct, Adsorption%DensSurfAtoms(SurfSide)/Adsorption%AreaIncrease)
        CALL PartitionFuncSurf(iSpec, WallTemp, VarPartitionFuncWall)
        nu_des = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct / VarPartitionFuncWall)
        rate = nu_des * exp(-E_des/WallTemp)
        
        Adsorption%ProbSigDes(p,q,SurfSide,iSpec,iProbSigma) = rate * dt
        IF (Adsorption%ProbSigDes(p,q,SurfSide,iSpec,iProbSigma).GT.1.) THEN
          Adsorption%ProbSigDes(p,q,SurfSide,iSpec,iProbSigma) = 1.
        END IF
      END DO
    ELSE
      Adsorption%ProbSigDes(p,q,SurfSide,iSpec,(1+(iSpec-1)*36):(36*iSpec)) = 0.0
    END IF
#if (PP_TimeDiscMethod==42)
    IF (.NOT.DSMC%ReservoirRateStatistic) THEN
      DO iProbSigma = (1+(iSpec-1)*36),(36*iSpec)
        Adsorption%AdsorpInfo(iSpec)%MeanProbDes = Adsorption%AdsorpInfo(iSpec)%MeanProbDes &
                                                  + Adsorption%ProbSigDes(p,q,SurfSide,iSpec,iProbSigma) &
                                                  * Adsorption%ProbSigma(p,q,SurfSide,iSpec,iProbSigma)
      END DO
    END IF
#endif
!===================================================================================================================================
  END IF ! DSMC%WallModel
! (QCA-quasi-chemical-approximation not considered)
  
  END DO
  END DO
  END DO
END DO
#if (PP_TimeDiscMethod==42)
IF (.NOT.DSMC%ReservoirRateStatistic) THEN
  Adsorption%AdsorpInfo(:)%MeanProbDes = Adsorption%AdsorpInfo(:)%MeanProbDes / (nSurfSample * nSurfSample * SurfMesh%nSides)
END IF
#endif
END SUBROUTINE CalcDesorbProb

SUBROUTINE PartitionFuncGas(iSpec, Temp, VarPartitionFuncGas)
!===================================================================================================================================
! Init of DSMC Vars
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,       ONLY: PlanckConst
USE MOD_DSMC_Vars,          ONLY: SpecDSMC, PolyatomMolDSMC
USE MOD_Particle_Vars,      ONLY: BoltzmannConst, Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
REAL, INTENT(IN)              :: Temp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: VarPartitionFuncGas
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL, PARAMETER               :: Pi=3.14159265358979323846_8
INTEGER                       :: iPolyatMole, iDOF
REAL                          :: Qtra, Qrot, Qvib!, Qelec
!===================================================================================================================================
  Qtra = (2. * Pi * Species(iSpec)%MassIC * BoltzmannConst * Temp / (PlanckConst**2))**(1.5)
  IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
    IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
      IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec) THEN
        Qrot = Temp / (SpecDSMC(iSpec)%SymmetryFactor * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))
      ELSE
        Qrot = SQRT(Pi) / SpecDSMC(iSpec)%SymmetryFactor * SQRT(Temp**3/( PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)    &
                                                                        * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(2)    &
                                                                        * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(3)))
      END IF
      Qvib = 1.
      DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
        Qvib = Qvib * EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / (2. * Temp)) &
                / (1. - EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / Temp))
      END DO
    ELSE
      Qrot = Temp / (SpecDSMC(iSpec)%SymmetryFactor * SpecDSMC(iSpec)%CharaTRot)
      Qvib = EXP(-SpecDSMC(iSpec)%CharaTVib / (2. * Temp)) / (1. - EXP(-SpecDSMC(iSpec)%CharaTVib / Temp))
    END IF
  ELSE
    Qrot = 1.
    Qvib = 1.
  END IF
!   Qelec = 0.
!   DO iDOF=1, SpecDSMC(iSpec)%NumElecLevels
!     Qelec = Qelec + SpecDSMC(iSpec)%ElectronicState(1,iDOF) * EXP(-SpecDSMC(iSpec)%ElectronicState(2,iDOF) / Temp)
!   END DO
  VarPartitionFuncGas = Qtra * Qrot * Qvib! * Qelec

END SUBROUTINE PartitionFuncGas

SUBROUTINE PartitionFuncAct(iSpec, Temp, VarPartitionFuncAct, Surfdensity)
!===================================================================================================================================
! Partitionfunction of activated complex
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars,       ONLY : PlanckConst
USE MOD_DSMC_Vars,          ONLY : SpecDSMC, PolyatomMolDSMC
USE MOD_Particle_Vars,      ONLY : BoltzmannConst, Species
!===================================================================================================================================
IMPLICIT NONE
!===================================================================================================================================
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
REAL, INTENT(IN)              :: Temp
REAL, INTENT(IN)              :: Surfdensity
!===================================================================================================================================
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: VarPartitionFuncAct
!===================================================================================================================================
! LOCAL VARIABLES
REAL, PARAMETER               :: Pi=3.14159265358979323846_8
INTEGER                       :: iPolyatMole, iDOF
REAL                          :: Qtra, Qrot, Qvib
!===================================================================================================================================
  Qtra = (2. * Pi * Species(iSpec)%MassIC * BoltzmannConst * Temp / (PlanckConst**2))/Surfdensity
  IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
    IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
      IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec) THEN
        Qrot = Temp / (SpecDSMC(iSpec)%SymmetryFactor * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))
      ELSE
        Qrot = SQRT(Pi) / SpecDSMC(iSpec)%SymmetryFactor * SQRT(Temp**3/( PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)    &
                                                                        * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(2)    &
                                                                        * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(3)))
      END IF
      Qvib = 1.
      DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
        Qvib = Qvib * EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / (2. * Temp)) &
                / (1. - EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / Temp))
      END DO
    ELSE
      Qrot = Temp / (SpecDSMC(iSpec)%SymmetryFactor * SpecDSMC(iSpec)%CharaTRot)
      Qvib = EXP(-SpecDSMC(iSpec)%CharaTVib / (2. * Temp)) / (1. - EXP(-SpecDSMC(iSpec)%CharaTVib / Temp))
    END IF
  ELSE
    Qrot = 1.
    Qvib = 1.
  END IF
  VarPartitionFuncAct = Qtra * Qrot * Qvib

END SUBROUTINE PartitionFuncAct

SUBROUTINE PartitionFuncSurf(iSpec, Temp, VarPartitionFuncSurf)
!===================================================================================================================================
! partition function of adsorbates
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars,       ONLY : PlanckConst
USE MOD_DSMC_Vars,          ONLY : SpecDSMC, PolyatomMolDSMC
USE MOD_Particle_Vars,      ONLY : BoltzmannConst
!===================================================================================================================================
IMPLICIT NONE
!===================================================================================================================================
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
REAL, INTENT(IN)              :: Temp
!===================================================================================================================================
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: VarPartitionFuncSurf
!===================================================================================================================================
! LOCAL VARIABLES
REAL, PARAMETER               :: Pi=3.14159265358979323846_8
INTEGER                       :: iPolyatMole, iDOF
REAL                          :: Qtra, Qrot, Qvib
!===================================================================================================================================
  Qtra = 1.
  IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
    IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
      Qvib = 1.
      DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
        Qvib = Qvib * EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / (2. * Temp)) &
                / (1. - EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / Temp))
      END DO
    ELSE
      Qvib = EXP(-SpecDSMC(iSpec)%CharaTVib / (2. * Temp)) / (1. - EXP(-SpecDSMC(iSpec)%CharaTVib / Temp))
    END IF
  ELSE
    Qvib = 1.
  END IF
  Qrot = 1.
  VarPartitionFuncSurf = Qtra * Qrot * Qvib

END SUBROUTINE PartitionFuncSurf

REAL FUNCTION Calc_Beta_Adsorb(iSpec,Xi_Total,Constant_Adsorb)
!===================================================================================================================================
! Calculates the Beta coefficient for molecular (non-associative) adsorption reaction
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_PARTICLE_Vars,      ONLY : BoltzmannConst
  USE MOD_DSMC_Vars,          ONLY : Adsorption
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)            :: iSpec
  REAL, INTENT(IN)               :: Xi_Total
  REAL, INTENT(IN)               :: Constant_Adsorb
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

  IF((Adsorption%Ads_Powerfactor(iSpec) - 0.5 + Xi_Total).GT.0.0) THEN
    Calc_Beta_Adsorb = Adsorption%Ads_Prefactor(iSpec) * Constant_Adsorb &
      * BoltzmannConst**(0.5 - Adsorption%Ads_Powerfactor(iSpec) ) * GAMMA(Xi_Total) &
      / GAMMA(Adsorption%Ads_Powerfactor(iSpec) - 0.5 + Xi_Total) 
  ELSE
    Calc_Beta_Adsorb = 0.0
  END IF

END FUNCTION Calc_Beta_Adsorb

REAL FUNCTION Calc_Beta_Diss(ReactNum,iSpec,Xi_Total,Constant_Adsorb)
!===================================================================================================================================
! Calculates the Beta coefficient for molecular (non-associative) adsorption reaction
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_PARTICLE_Vars,      ONLY : BoltzmannConst
  USE MOD_DSMC_Vars,          ONLY : Adsorption
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)            :: ReactNum
  INTEGER, INTENT(IN)            :: iSpec
  REAL, INTENT(IN)               :: Xi_Total
  REAL, INTENT(IN)               :: Constant_Adsorb
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

  IF((Adsorption%Diss_Powerfactor(ReactNum,iSpec) - 0.5 + Xi_Total).GT.0.0) THEN
    Calc_Beta_Diss = Adsorption%Diss_Prefactor(ReactNum,iSpec) * Constant_Adsorb &
      * BoltzmannConst**(0.5 - Adsorption%Diss_Powerfactor(ReactNum,iSpec) ) * GAMMA(Xi_Total) &
      / GAMMA(Adsorption%Diss_Powerfactor(ReactNum,iSpec) - 0.5 + Xi_Total) 
  ELSE
    Calc_Beta_Diss = 0.0
  END IF

END FUNCTION Calc_Beta_Diss

REAL FUNCTION Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,iSpec,Surfpos,IsAdsorption)
!===================================================================================================================================
! Calculates the Heat of adsorption for given species and given surface position
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_PARTICLE_Vars,      ONLY : BoltzmannConst, nSpecies
  USE MOD_DSMC_Vars,          ONLY : Adsorption, SurfDistInfo, SpecDSMC
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)            :: subsurfxi, subsurfeta, SurfSideID
  INTEGER, INTENT(IN)            :: iSpec, Surfpos
  LOGICAL, INTENT(IN)            :: IsAdsorption
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                        :: Coordination, i, j, Indx, Indy
  REAL , ALLOCATABLE             :: x(:)!, D_AL(:), delta(:) !, z(:)
  INTEGER , ALLOCATABLE          :: m(:)!, Neigh_bondorder(:)
  INTEGER                        :: bondorder
  REAL                           :: D_AB, D_AX, D_BX
  REAL                           :: Heat_A, Heat_B
  REAL                           :: A, B, sigma
!   REAL                           :: Heat_D_AL
!   INTEGER                        :: neighSpec, neighSpec2, Coord2, Coord3, ReactNum, nNeigh_interactions
!===================================================================================================================================
  Coordination = Adsorption%Coordination(iSpec)
  ALLOCATE( x(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom) )
!   ALLOCATE( z(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom) )
  ALLOCATE( m(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom) )
  x(:) = 1. ! averaged bond-index for surface atoms 
  m(:) = 1  ! number of adsorbates belonging to surface atom
!   z(:) = 1.
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
      x(j) = 1. / REAL(m(j))
    END DO
  END IF
!   z(:) = x(:)
  ! calculate local scaling factor for chosen surface site
  DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom
    x(j) = x(j) / REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom)
    sigma = sigma + (2.*x(j) - x(j)**2.)
  END DO
!   sigma = sigma / REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom)
  
  ! Testing if the adsorption particle is an atom or molecule, if molecule: is it polyatomic?
  ! and calculate right heat of adsorption to surface atoms
  IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
    IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
      D_AB = Adsorption%EDissBond(0,iSpec)
      SELECT CASE(Coordination)
      CASE(1) ! hollow (radical with localized electron like NH2)
        Heat_A = Adsorption%HeatOfAdsZero(iSpec) * sigma
        Calc_Adsorb_Heat = Heat_A**2/(D_AB+Heat_A)
      CASE(2) ! bridge
        IF (Adsorption%DiCoord(iSpec).EQ.1) THEN ! dicoordination (HCOOH --> M--(HC)O-O(H)--M) (M--O bond)
          D_AX = Adsorption%EDissBondAdsorbPoly(0,iSpec) ! Bond HC--O
          Heat_A = Adsorption%HeatOfAdsZero(iSpec)! * (2.*z(1) - z(1)**2)
          A = Heat_A**2./(D_AX+Heat_A)
          D_BX = Adsorption%EDissBondAdsorbPoly(1,iSpec) ! Bond O--H
          Heat_B = Adsorption%HeatOfAdsZero(iSpec)! * (2.*z(2) - z(2)**2)
          B = Heat_B**2./(D_BX+Heat_B)
          Calc_Adsorb_Heat = ( A*B*( A + B ) + D_AB*( A - B )**2. ) / ( A*B + D_AB*( A + B ) ) * sigma
        
!           Heat_A = Adsorption%HeatOfAdsZero(iSpec) * (2.*x(1) - x(1)**2)
!           Heat_M = Adsorption%HeatOfAdsZeroM(iSpec) * (2.*x(1) - x(1)**2)
!           A = Heat_A * (1- ( (m*Heat_B/( m*Heat_A + Heat_B ))**2) )
!           
!           Heat_A = Adsorption%HeatOfAdsZero2(iSpec) * (2.*x(2) - x(2)**2)
!           Heat_M = Adsorption%HeatOfAdsZeroR(iSpec) * (2.*x(2) - x(2)**2)
!           B = Heat_A * (1- ( (mtilde*Heat_B/( mtilde*Heat_A + Heat_B ))**2) )
          
        ELSE IF (Adsorption%DiCoord(iSpec).EQ.2) THEN ! chelate binding (NO2 --> M--O-N-O--M)
          D_AX = Adsorption%EDissBondAdsorbPoly(0,iSpec) ! Bond O--N
          Heat_A = Adsorption%HeatOfAdsZero(iSpec)! * (2.*z(1) - z(1)**2)
          Heat_A = Heat_A**2/(D_AX+Heat_A)
          D_BX = Adsorption%EDissBondAdsorbPoly(1,iSpec) ! Bond N--O
          Heat_B = Adsorption%HeatOfAdsZero(iSpec)! * (2.*z(2) - z(2)**2)
          Heat_B = Heat_B**2/(D_BX+Heat_B)
          A = Heat_A**2. * ( Heat_A + 2.*Heat_B ) / ( Heat_A + Heat_B )**2. * sigma
          B = Heat_B**2. * ( Heat_B + 2.*Heat_A ) / ( Heat_A + Heat_B )**2. * sigma
          Calc_Adsorb_Heat = A + B
        END IF
      CASE(3) ! on-top (closed shell or open shell with unlocalized electron like CO)
        Heat_A = Adsorption%HeatOfAdsZero(iSpec) * sigma
        Calc_Adsorb_Heat = Heat_A**2./(D_AB+Heat_A)
      END SELECT
    ELSE 
      D_AB = Adsorption%EDissBond(0,iSpec)
      SELECT CASE(Coordination)
      CASE(1) ! hollow (radical with localized electron like C-H)
        Heat_A = Adsorption%HeatOfAdsZero(iSpec) * sigma
        Calc_Adsorb_Heat = (Heat_A**2)/(D_AB+Heat_A)
      CASE(2) ! bridge (closed shell like O2)
        IF (Adsorption%DiCoord(iSpec).EQ.1) THEN
          Heat_A = Adsorption%HeatOfAdsZero(iSpec) !* (2.*z(1) - z(1)**2)
          Heat_B = Adsorption%HeatOfAdsZero(iSpec) !* (2.*z(2) - z(2)**2)
          A = Heat_A**2. * ( Heat_A + 2.*Heat_B ) / ( Heat_A + Heat_B )**2. 
          B = Heat_B**2. * ( Heat_B + 2.*Heat_A ) / ( Heat_A + Heat_B )**2.
          Calc_Adsorb_Heat = ( A*B*( A + B ) + D_AB*( A - B )**2. ) / ( A*B + D_AB*( A + B ) ) * sigma
        ELSE
          Heat_A = Adsorption%HeatOfAdsZero(iSpec) * sigma
          Calc_Adsorb_Heat = Heat_A**2./(D_AB+Heat_A)
        END IF
      CASE(3) ! on-top (closed shell or open shell with unlocalized electron like CO)
        Heat_A = Adsorption%HeatOfAdsZero(iSpec) * sigma
        Calc_Adsorb_Heat = Heat_A**2./(D_AB+Heat_A)
      END SELECT
    END IF
  ELSE
    Calc_Adsorb_Heat = Adsorption%HeatOfAdsZero(iSpec) * sigma
  END IF
  
  ! routine wird nicht benutz, da höheres Rauschen und größerer Rechenaufwand aber aufgehoben für spätere Einsicht
!   ! calculate additional heat of adsorption for direct interaction (attraction of associating adsorabtes)
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
!             IF ( neighSpec.EQ.Adsorption%AssocReact(1,ReactNum,iSpec) .AND. &
!                  (Adsorption%AssocReact(2,ReactNum,iSpec).NE.0)) THEN
!               D_AL(l) = Adsorption%EDissBond((Adsorption%DissNum+ReactNum),iSpec)
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
! Calculates the Activation energy for a given reaction
! A_Reactant_ads + B_Reactant_ads --> A_Product_ads + B_Product_ads
! Adsorption --> forward reaction
! Forward reaction is defined by D_Educt > D_Products
! Examples:
! (1)
! O2 desorbed directly to gasphase from reaction of two O (O_ads + O_ads -> O2_g): 
! ==> forward reaction: O2_g + (-)_ads -> O_ads + O_ads
! ==> IsAdsorption = .FALSE.
! ==> Heat_Reactant_A = Heat_O2_g = 0. | Heat_Product_A_ads = Heat_Product_B_ads = Heat_O_ads
! (2)
! adsorbed CH radical reacts with adsorbed O-atom to adsorbed C-atom and OH-radical (CH_ads + O_ads -> C_ads + OH_ads): 
! ==> forward reaction: CH_ads + O_ads -> C_ads + OH_ads
! ==> IsAdsorption = .TRUE.
! ==> Heat_Reactant_A = Heat_CH_ads | Heat_Reactant_B = Heat_O_ads | Heat_Product_A = Heat_C_ads | Heat_Product_B = Heat_OH_ads
! (3)
! adsorbed OH radical reacts with adsorbed C-atom to adsorbed O-atom and gasphase CH-radical (OH_ads + C_ads -> O_ads + CH_g): 
! ==> forward reaction: CH_g + O_ads -> C_ads + OH_ads
! ==> IsAdsorption = .FALSE.
! ==> Heat_Reactant_A = Heat_CH_g = 0. | Heat_Reactant_B = Heat_O_ads | Heat_Product_A = Heat_C_ads | Heat_Product_B = Heat_OH_ads
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
  IF ( (D_Reactant_A +D_Reactant_B -D_Product_A -D_Product_B).GE.0. ) Forward = .TRUE.
  
  Delta_H = ( Heat_Reactant_A +Heat_Reactant_B -Heat_Product_A -Heat_Product_B ) &
          + ( D_Reactant_A +D_Reactant_B -D_Product_A -D_Product_B )
  IF (Forward) THEN
    Calc_E_Act = 0.5 * ( Delta_H + (Heat_Product_A*Heat_Product_B / (Heat_Product_A+Heat_Product_B)) )
  ELSE
    Calc_E_Act = 0.5 * ( Delta_H + (Heat_Reactant_A*Heat_Reactant_B / (Heat_Reactant_A+Heat_Reactant_B)) )
  END IF
  IF (Calc_E_Act.LT.0.) Calc_E_Act = 0.

END FUNCTION Calc_E_Act

END MODULE MOD_DSMC_SurfModel_Tools
