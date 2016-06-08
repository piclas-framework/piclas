#include "boltzplatz.h"

MODULE MOD_DSMC_SurfModelInit
!===================================================================================================================================
! Initialization of DSMC
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE InitDSMCSurfModel
  MODULE PROCEDURE InitDSMCSurfModel
END INTERFACE

INTERFACE FinalizeDSMCSurfModel
  MODULE PROCEDURE FinalizeDSMCSurfModel
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC                       :: InitDSMCSurfModel
PUBLIC                       :: FinalizeDSMCSurfModel
!===================================================================================================================================

CONTAINS


SUBROUTINE InitDSMCSurfModel()
!===================================================================================================================================
! Init of DSMC Vars
!===================================================================================================================================
! MODULES
USE MOD_Globals,                ONLY : abort
USE MOD_Mesh_Vars,              ONLY : nElems, nBCSides, BC
USE MOD_DSMC_Vars,              ONLY : Adsorption
USE MOD_PARTICLE_Vars,          ONLY : nSpecies, PDM
USE MOD_PARTICLE_Vars,          ONLY : KeepWallParticles, PEM
USE MOD_ReadInTools
USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh, PartBound
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                                      !
  CHARACTER(32)                    :: hilf
  INTEGER                          :: iSpec, iSide, iSurf, IDcounter
!===================================================================================================================================
KeepWallParticles = GETLOGICAL('Particles-KeepWallParticles','.FALSE.')
IF (KeepWallParticles) THEN
  ALLOCATE(PDM%ParticleAtWall(1:PDM%maxParticleNumber)  , &
          PDM%PartAdsorbSideIndx(1:3,1:PDM%maxParticleNumber))
  PDM%ParticleAtWall(1:PDM%maxParticleNumber) = .FALSE.
  ALLOCATE(PEM%wNumber(1:nElems))
END IF
! allocate info and constants
#if (PP_TimeDiscMethod==42)
ALLOCATE( Adsorption%AdsorpInfo(1:nSpecies))
#endif
ALLOCATE( Adsorption%InitStick(1:SurfMesh%nSides,1:nSpecies),& 
          Adsorption%PrefactorStick(1:SurfMesh%nSides,1:nSpecies),& 
          Adsorption%Adsorbexp(1:SurfMesh%nSides,1:nSpecies),& 
          Adsorption%Nu_a(1:SurfMesh%nSides,1:nSpecies),& 
          Adsorption%Nu_b(1:SurfMesh%nSides,1:nSpecies),& 
          Adsorption%DesorbEnergy(1:SurfMesh%nSides,1:nSpecies),& 
          Adsorption%Intensification(1:SurfMesh%nSides,1:nSpecies))
ALLOCATE( Adsorption%HeatOfAdsZero(1:nSpecies),& 
!           Adsorption%DissResultsSpecNum(1:nSpecies),&
!           Adsorption%DissResultsSpec(1:2,1:(2*(nSpecies-1)-1),1:nSpecies),&
!           Adsorption%EDissBond(1:nSpecies,1:nSpecies,1:nSpecies),& 
          Adsorption%Coordination(1:nSpecies))
! initialize info and constants
DO iSpec = 1,nSpecies
#if (PP_TimeDiscMethod==42)
  Adsorption%AdsorpInfo(iSpec)%MeanProbAds = 0.
  Adsorption%AdsorpInfo(iSpec)%MeanProbDes = 0.
  Adsorption%AdsorpInfo(iSpec)%WallCollCount = 0
  Adsorption%AdsorpInfo(iSpec)%NumOfAds = 0
  Adsorption%AdsorpInfo(iSpec)%NumOfDes = 0
#endif
  WRITE(UNIT=hilf,FMT='(I2)') iSpec
  Adsorption%InitStick(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-InitialStick','0.')
  Adsorption%PrefactorStick(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-PrefactorStick','0.')
  Adsorption%Adsorbexp(:,iSpec) = GETINT('Part-Species'//TRIM(hilf)//'-Adsorbexp','1')
  Adsorption%Nu_a(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Nu-a','0.')
  Adsorption%Nu_b(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Nu-b','0.')
  Adsorption%DesorbEnergy(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Desorption-Energy-K','1.')
  Adsorption%Intensification(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Intensification-K','0.')
  Adsorption%Coordination(iSpec) = GETINT('Part-Species'//TRIM(hilf)//'-Coordination','0')
  Adsorption%HeatOfAdsZero(iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-HeatOfAdsorption-K','0.')
END DO
#if (PP_TimeDiscMethod==42)
  Adsorption%TPD = GETLOGICAL('Particles-DSMC-Adsorption-doTPD','.FALSE.')
  Adsorption%TPD_beta = GETREAL('Particles-DSMC-Adsorption-TPD-Beta','0.')
  Adsorption%TPD_Temp = 0.
#endif
! allocate and initialize adsorption variables
ALLOCATE( Adsorption%MaxCoverage(1:SurfMesh%nSides,1:nSpecies),& 
          Adsorption%Coverage(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides,1:nSpecies),&
          Adsorption%ProbAds(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides,1:nSpecies),&
          Adsorption%ProbDes(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides,1:nSpecies,1:36*nSpecies),&
          Adsorption%SumDesorbPart(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides,1:nSpecies),&
          Adsorption%SumAdsorbPart(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides,1:nSpecies),&
          Adsorption%SurfSideToGlobSideMap(1:SurfMesh%nSides),&
          Adsorption%Sigma(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides,1:nSpecies,1:36*nSpecies),&
          Adsorption%ProbSigma(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides,1:nSpecies,1:36*nSpecies),&
          Adsorption%DensSurfAtoms(1:SurfMesh%nSides))
IDcounter = 0         
DO iSide = 1,nBCSides 
  IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(BC(iSide))).EQ.PartBound%ReflectiveBC) THEN
    IDcounter = IDcounter + 1
    Adsorption%SurfSideToGlobSideMap(IDcounter) = iSide
  END IF
END DO
DO iSurf = 1,SurfMesh%nSides
  WRITE(UNIT=hilf,FMT='(I2)') iSurf
  Adsorption%DensSurfAtoms(iSurf) = GETREAL('Particles-Surface'//TRIM(hilf)//'-AtomsDensity','1.5E+19')
END DO
DO iSpec = 1,nSpecies
  WRITE(UNIT=hilf,FMT='(I2)') iSpec
  Adsorption%MaxCoverage(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-MaximumCoverage','0.')
  Adsorption%Coverage(:,:,:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-InitialCoverage','0.')
END DO
Adsorption%ProbAds(:,:,:,:) = 0.
Adsorption%ProbDes(:,:,:,:,:) = 0.
Adsorption%SumDesorbPart(:,:,:,:) = 0
Adsorption%SumAdsorbPart(:,:,:,:) = 0
Adsorption%Sigma(:,:,:,:,:) = 1.
Adsorption%ProbSigma(:,:,:,:,:) = 0.

CALL Init_SurfDist()

END SUBROUTINE InitDSMCSurfModel

SUBROUTINE Init_SurfDist()
!===================================================================================================================================
! Initializing surface distibution Model for calculating of coverage effects on heat of adsorption
!===================================================================================================================================
  USE MOD_Globals,                ONLY : abort, MPIRoot, UNIT_StdOut
  USE MOD_ReadInTools
  USE MOD_Particle_Vars,          ONLY : nSpecies, Species, KeepWallParticles
  USE MOD_DSMC_Vars,              ONLY : Adsorption, SurfDistInfo
  USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh
!===================================================================================================================================
  IMPLICIT NONE
!=================================================================================================================================== 
! Local variable declaration
  INTEGER                          :: SurfSideID, subsurfxi, subsurfeta, iSpec, iInterAtom
  INTEGER                          :: surfsquare, dist, Adsorbates
  INTEGER                          :: Surfpos, Surfnum, Indx, Indy, UsedSiteMapPos
  REAL                             :: RanNum
  INTEGER                          :: xpos, ypos
  INTEGER                          :: Coord, nSites, nInterAtom, nNeighbours
!===================================================================================================================================
! position of binding sites in the surface lattice (rectangular lattice)
!------------[        surfsquare       ]--------------
!             |       |       |       |
!         3---2---3---2---3---2---3---2---3
!         |       |       |       |       |
!         2   1   2   1   2   1   2   1   2
!         |       |       |       |       |
!         3---2---3---2---3---2---3---2---3
!         |       |       |       |       |
!         2   1   2   1   2   1   2   1   2
!         |       |       |       |       |
!         3---2---3---2---3---2---3---2---3
!         |       |       |       |       |
!         2   1   2   1   2   1   2   1   2
!         |       |       |       |       |
!         3---2---3---2---3---2---3---2---3
!         |       |       |       |       |
!         2   1   2   1   2   1   2   1   2
!         |       |       |       |       |
!         3---2---3---2---3---2---3---2---3
! For now:
! Neighbours are all sites, that have the same binding surface atom. 
! Except for top sites(3) they also interact with the next top site.
SWRITE(UNIT_stdOut,'(A)')' INIT SURFACE DISTRIBUTION...'
ALLOCATE(SurfDistInfo(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides))
DO SurfSideID=1,SurfMesh%nSides
  DO subsurfeta = 1,nSurfSample
  DO subsurfxi = 1,nSurfSample
    ALLOCATE( SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(1:3),&
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(1:3))
  END DO
  END DO
END DO
IF (.NOT.KeepWallParticles) THEN
  surfsquare = GETINT('Particles-DSMC-AdsorptionSites','1000')
  surfsquare = INT(SQRT(REAL(surfsquare))) - 1
END IF

DO SurfSideID=1,SurfMesh%nSides
  DO subsurfeta = 1,nSurfSample
  DO subsurfxi = 1,nSurfSample
    IF (KeepWallParticles) THEN ! does not work with vMPF
      surfsquare = INT(Adsorption%DensSurfAtoms(SurfSideID) &
                    * SurfMesh%SurfaceArea(subsurfxi,subsurfeta,SurfSideID) &
                    / Species(1)%MacroParticleFactor)
      surfsquare = INT(SQRT(REAL(surfsquare))) - 1
    END IF
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(1) = INT(surfsquare**2)
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(2) = INT( 2*(surfsquare*(surfsquare+1)) )
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(3) = INT((surfsquare+1)**2)
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(:) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(:)
    
    ALLOCATE( SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(1:nSpecies,1:surfsquare+1,1:surfsquare+1))
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(:,:,:) = 0
    
    ALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1:3))
    DO Coord = 1,3    
      SELECT CASE (Coord)
      CASE(1)
        nSites = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Coord)
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom = 4
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nNeighbours = 16
      CASE(2)
        nSites = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Coord)
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom = 2
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nNeighbours = 14
      CASE(3)
        nSites = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Coord)
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom = 1
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nNeighbours = 12
      END SELECT
        
      nInterAtom = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
      nNeighbours = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nNeighbours
      ALLOCATE( SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(1:nSites,nInterAtom),&
                SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(1:nSites,nInterAtom),&
                SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(1:nSites,1:nNeighbours),&
                SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighSite(1:nSites,1:nNeighbours))
      ALLOCATE( SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(1:nSites),&
                SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(1:nSites))
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(:) = 0
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(:) = 0
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(:,:) = 0
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(:,:) = 0
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(:,:) = 0
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighSite(:,:) = 0
    END DO
  END DO
  END DO
END DO

DO SurfSideID = 1,SurfMesh%nSides
DO subsurfeta = 1,nSurfSample
DO subsurfxi = 1,nSurfSample
  ! surfsquare chosen from nSite(1) for correct SurfIndx definitions
  surfsquare = INT(SQRT(REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(1))))
  ! allocate and define surface indexes for adsorbate distribution and build mapping of respective bondatoms and neighbours      
  Indx = 1
  Indy = 1
  DO Surfpos = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(1)
    IF (Indx.GT.surfsquare) THEN
      Indx = 1
      Indy = Indy + 1
    END IF
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%UsedSiteMap(Surfpos) = Surfpos
    ! mapping respective neighbours first hollow then bridge then top
    ! hollow
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighPos(Surfpos,1) = Surfpos - surfsquare - 1
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighPos(Surfpos,2) = Surfpos - surfsquare
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighPos(Surfpos,3) = Surfpos - surfsquare + 1
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighPos(Surfpos,4) = Surfpos - 1
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighPos(Surfpos,5) = Surfpos + 1
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighPos(Surfpos,6) = Surfpos + surfsquare - 1
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighPos(Surfpos,7) = Surfpos + surfsquare
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighPos(Surfpos,8) = Surfpos + surfsquare + 1
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighSite(Surfpos,1:8) = 1
    ! bridge
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighPos(Surfpos,9) = Surfpos +(surfsquare+1)*(Indy-1)
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighPos(Surfpos,10) = Surfpos +surfsquare +(surfsquare+1)*(Indy-1)
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighPos(Surfpos,11) = Surfpos +surfsquare +(surfsquare+1)*(Indy-1) +1
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighPos(Surfpos,12) = Surfpos +surfsquare +(surfsquare+1)*(Indy)
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighSite(Surfpos,9:12) = 2
    ! top
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighPos(Surfpos,13) = Surfpos + (Indy-1)
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighPos(Surfpos,14) = Surfpos + 1 + (Indy-1)
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighPos(Surfpos,15) = Surfpos + (surfsquare+1) + (Indy-1)
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighPos(Surfpos,16) = Surfpos + (surfsquare+1) + 1 + (Indy-1)
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighSite(Surfpos,13:16) = 3
    ! account for empty edges
    IF (Indy .EQ. 1) SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighSite(Surfpos,1:3) = 0
    IF (Indy .EQ. surfsquare) SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighSite(Surfpos,6:8) = 0
    IF (Indx .EQ. 1) THEN
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighSite(Surfpos,1) = 0
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighSite(Surfpos,4) = 0
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighSite(Surfpos,6) = 0
    END IF
    IF (Indx .EQ. surfsquare) THEN
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighSite(Surfpos,3) = 0
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighSite(Surfpos,5) = 0
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighSite(Surfpos,8) = 0
    END IF
    ! mapping respective bond atoms
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%BondAtomIndx(Surfpos,1) = Indx
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%BondAtomIndy(Surfpos,1) = Indy
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%BondAtomIndx(Surfpos,2) = Indx+1
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%BondAtomIndy(Surfpos,2) = Indy
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%BondAtomIndx(Surfpos,3) = Indx
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%BondAtomIndy(Surfpos,3) = Indy+1
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%BondAtomIndx(Surfpos,4) = Indx+1
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%BondAtomIndy(Surfpos,4) = Indy+1
    Indx = Indx + 1
  END DO
  Indx = 1
  Indy = 1
  DO Surfpos = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(2)
    IF (Indx.GT.(2*surfsquare+1)) THEN
      Indx = 1
      Indy = Indy + 1
    END IF
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%UsedSiteMap(Surfpos) = Surfpos
    IF (Indx .LE. surfsquare) THEN ! surface atoms are LEFT an RIGHT of adsorbate site
      ! mapping respective neighbours first hollow then bridge then top
      ! hollow
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,1) = Surfpos-surfsquare -(surfsquare+1)*(Indy-1) -1
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,2) = Surfpos-surfsquare -(surfsquare+1)*(Indy-1)
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,3) = Surfpos-surfsquare -(surfsquare+1)*(Indy-1) +1
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,4) = Surfpos -(surfsquare+1)*(Indy-1) -1
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,5) = Surfpos -(surfsquare+1)*(Indy-1)
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,6) = Surfpos -(surfsquare+1)*(Indy-1) +1
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,1:6) = 1
      ! bridge
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,7) = Surfpos - (surfsquare+1)
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,8) = Surfpos - (surfsquare+1) + 1
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,9) = Surfpos - 1
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,10) = Surfpos + 1
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,11) = Surfpos + surfsquare
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,12) = Surfpos + surfsquare + 1
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,7:12) = 2
      ! top
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,13) = Surfpos -(surfsquare)*(Indy-1)
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,14) = Surfpos +1 -(surfsquare)*(Indy-1)
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,13:14) = 3
      ! account for empty edges
      IF (Indy .EQ. 1) THEN
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,1:3) = 0
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,7:8) = 0
      END IF
      IF (Indy .EQ. surfsquare+1) THEN
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,4:6) = 0
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,11:12) = 0
      END IF
      IF (Indx .EQ. 1) THEN
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,1) = 0
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,4) = 0
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,9) = 0
      END IF
      IF (Indx .EQ. surfsquare) THEN
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,3) = 0
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,6) = 0
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,10) = 0
      END IF
      ! mapping respective bond atoms
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%BondAtomIndx(Surfpos,1) = Indx
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%BondAtomIndy(Surfpos,1) = Indy
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%BondAtomIndx(Surfpos,2) = Indx+1
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%BondAtomIndy(Surfpos,2) = Indy
    ELSE ! surface atoms are TOP and DOWN of adsorbate site
      ! mapping respective neighbours first hollow then bridge then top
      ! hollow
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,1) = Surfpos -(2*surfsquare) &
                                                                                  -(surfsquare+1)*(Indy-1) -1
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,2) = Surfpos -(2*surfsquare)&
                                                                                  -(surfsquare+1)*(Indy-1)
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,3) = Surfpos -surfsquare &
                                                                                  -(surfsquare+1)*(Indy-1) -1
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,4) = Surfpos - surfsquare -(surfsquare+1)*(Indy-1)
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,5) = Surfpos -(surfsquare+1)*(Indy-2) -1
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,6) = Surfpos -(surfsquare+1)*(Indy-2)
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,1:6) = 1
      ! bridge
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,7) = Surfpos - surfsquare - (surfsquare+1)
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,8) = Surfpos - surfsquare - 1
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,9) = Surfpos - surfsquare
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,10) = Surfpos + (surfsquare+1) - 1
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,11) = Surfpos + (surfsquare+1)
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,12) = Surfpos + surfsquare + (surfsquare+1)
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,7:12) = 2
      ! top
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,13) = Surfpos -surfsquare*(Indy)
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighPos(Surfpos,14) = Surfpos -surfsquare*(Indy) +(surfsquare+1)
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,13:14) = 3
      ! account for empty edges
      IF (Indy .EQ. 1) THEN
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,1:2) = 0
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,7) = 0
      END IF
      IF (Indy .EQ. surfsquare) THEN
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,5:6) = 0
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,12) = 0
      END IF
      IF (Indx .EQ. (surfsquare+1)) THEN
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,1) = 0
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,3) = 0
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,5) = 0
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,8) = 0
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,10) = 0
      END IF
      IF (Indx .EQ. 2*surfsquare+1) THEN
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,2) = 0
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,4) = 0
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,6) = 0
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,9) = 0
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%NeighSite(Surfpos,11) = 0
      END IF
      ! mapping respective bond atoms
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%BondAtomIndx(Surfpos,1) = Indx - surfsquare
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%BondAtomIndy(Surfpos,1) = Indy 
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%BondAtomIndx(Surfpos,2) = Indx - surfsquare
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(2)%BondAtomIndy(Surfpos,2) = Indy+1
    END IF
    Indx = Indx + 1
  END DO
  Indx = 1
  Indy = 1
  DO Surfpos = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(3)
    IF (Indx.GT.surfsquare+1) THEN
      Indx = 1
      Indy = Indy + 1
    END IF
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(3)%UsedSiteMap(Surfpos) = Surfpos
    ! mapping respective neighbours first hollow then bridge then top
    ! hollow
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(3)%NeighPos(Surfpos,1) = Surfpos - (surfsquare+1) - 1
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(3)%NeighPos(Surfpos,2) = Surfpos - (surfsquare+1)
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(3)%NeighPos(Surfpos,3) = Surfpos - surfsquare + 1
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(3)%NeighPos(Surfpos,4) = Surfpos - 1
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(3)%NeighSite(Surfpos,1:4) = 1
    ! bridge
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(3)%NeighPos(Surfpos,5) = Surfpos + 1
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(3)%NeighPos(Surfpos,6) = Surfpos + surfsquare - 1
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(3)%NeighPos(Surfpos,7) = Surfpos + surfsquare
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(3)%NeighPos(Surfpos,8) = Surfpos + surfsquare + 1
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(3)%NeighSite(Surfpos,5:8) = 2
    ! top
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(3)%NeighPos(Surfpos,9) = Surfpos - surfsquare - 1
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(3)%NeighPos(Surfpos,10) = Surfpos - surfsquare
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(3)%NeighPos(Surfpos,11) = Surfpos - surfsquare + 1
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(3)%NeighPos(Surfpos,12) = Surfpos - surfsquare + 1
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(3)%NeighSite(Surfpos,9:12) = 3
    ! mapping respective bond atoms
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(3)%BondAtomIndx(Surfpos,1) = Indx
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(3)%BondAtomIndy(Surfpos,1) = Indy
    Indx = Indx + 1
  END DO
    
  DO iSpec = 1,nSpecies
    Adsorbates = INT(Adsorption%Coverage(subsurfxi,subsurfeta,SurfSideID,iSpec) &
                * SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Adsorption%Coordination(iSpec)))
    IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Adsorption%Coordination(iSpec)).LT.Adsorbates) THEN
      CALL abort(&
      __STAMP__&
      ,'Error in Init_SurfDist: Too many Adsorbates! - Choose lower Coverages')
    END IF
    
    ! distribute adsorbates randomly on the surface on the correct site and assign surface atom bond order
    dist = 1
    Coord = Adsorption%Coordination(iSpec)
    Surfnum = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord)
    DO WHILE (dist.LE.Adsorbates) 
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
  END DO
  
END DO
END DO
END DO
SWRITE(UNIT_stdOut,'(A)')' INIT SURFACE DISTRIBUTION DONE!'
    
END SUBROUTINE Init_SurfDist

SUBROUTINE FinalizeDSMCSurfModel()
!===================================================================================================================================
! Init of DSMC Vars
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars,              ONLY : Adsorption, SurfDistInfo
USE MOD_PARTICLE_Vars,          ONLY : PDM, PEM
USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: subsurfxi,subsurfeta,SurfSideID,Coord
!===================================================================================================================================
SDEALLOCATE(PDM%ParticleAtWall)
SDEALLOCATE(PDM%PartAdsorbSideIndx)
SDEALLOCATE(PEM%wNumber)
#if (PP_TimeDiscMethod==42)
SDEALLOCATE( Adsorption%AdsorpInfo)
#endif
SDEALLOCATE(Adsorption%InitStick)
SDEALLOCATE(Adsorption%PrefactorStick)
SDEALLOCATE(Adsorption%Adsorbexp)
SDEALLOCATE(Adsorption%Nu_a)
SDEALLOCATE(Adsorption%Nu_b)
SDEALLOCATE(Adsorption%DesorbEnergy)
SDEALLOCATE(Adsorption%Intensification)
SDEALLOCATE(Adsorption%HeatOfAdsZero)
SDEALLOCATE(Adsorption%DissResultsSpecNum)
SDEALLOCATE(Adsorption%DissResultsSpec)
SDEALLOCATE(Adsorption%EDissBond)
SDEALLOCATE(Adsorption%Coordination)

SDEALLOCATE(Adsorption%MaxCoverage)
SDEALLOCATE(Adsorption%Coverage)
SDEALLOCATE(Adsorption%ProbAds)
SDEALLOCATE(Adsorption%ProbDes)
SDEALLOCATE(Adsorption%SumDesorbPart)
SDEALLOCATE(Adsorption%SumAdsorbPart)
SDEALLOCATE(Adsorption%SurfSideToGlobSideMap)
SDEALLOCATE(Adsorption%Sigma)
SDEALLOCATE(Adsorption%ProbSigma)
SDEALLOCATE(Adsorption%DensSurfAtoms)

IF (ALLOCATED(SurfDistInfo)) THEN
DO SurfSideID=1,SurfMesh%nSides
  DO subsurfeta = 1,nSurfSample
  DO subsurfxi = 1,nSurfSample
!     DEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%HollowSite)
!     DEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%BridgeSite)
!     DEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%TopSite)
    DEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder)
    DEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain)
    DEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites)
    IF (ALLOCATED(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap)) THEN
      DO Coord = 1,3
        DEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap)
        DEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species)
        DEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx)
        DEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy)
        DEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos)
        DEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighSite)
      END DO
    END IF
  END DO
  END DO
END DO
DEALLOCATE(SurfDistInfo)
END IF

END SUBROUTINE FinalizeDSMCSurfModel

END MODULE MOD_DSMC_SurfModelInit