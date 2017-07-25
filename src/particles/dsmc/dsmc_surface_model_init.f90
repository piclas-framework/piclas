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
USE MOD_Mesh_Vars,              ONLY : nElems, BC
USE MOD_DSMC_Vars,              ONLY : Adsorption, DSMC!, CollisMode
USE MOD_PARTICLE_Vars,          ONLY : nSpecies, PDM
USE MOD_PARTICLE_Vars,          ONLY : KeepWallParticles, PEM
USE MOD_Particle_Mesh_Vars,     ONLY : nTotalSides
USE MOD_ReadInTools
USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh, nPartBound, PartBound
#if (PP_TimeDiscMethod==42)
USE MOD_DSMC_SurfModel_Tools,   ONLY : CalcDiffusion
#endif
#ifdef MPI
USE MOD_Particle_Boundary_Vars, ONLY : SurfCOMM
USE MOD_Particle_MPI_Vars     , ONLY : AdsorbSendBuf,AdsorbRecvBuf,SurfExchange
#endif
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  CHARACTER(32)                    :: hilf, hilf2
  INTEGER                          :: iSpec, iSide, iPartBound
  INTEGER                          :: SideID, PartBoundID
  REAL , ALLOCATABLE               :: Coverage_tmp(:,:)
#if (PP_TimeDiscMethod==42)
  INTEGER                          :: idiff
#endif
#ifdef MPI
  INTEGER                          :: iProc
#endif
!===================================================================================================================================
! do not initialize variables if processor does not have any surfaces on in domain
IF(.NOT.SurfMesh%SurfOnProc) RETURN
! initialize surface chemistry
KeepWallParticles = GETLOGICAL('Particles-KeepWallParticles','.FALSE.')
IF (KeepWallParticles) THEN
  ALLOCATE(PDM%ParticleAtWall(1:PDM%maxParticleNumber)  , &
          PDM%PartAdsorbSideIndx(1:3,1:PDM%maxParticleNumber))
  PDM%ParticleAtWall(1:PDM%maxParticleNumber) = .FALSE.
  ALLOCATE(PEM%wNumber(1:nElems))
END IF
! allocate info and constants
ALLOCATE( Adsorption%AdsorpInfo(1:nSpecies))
IF (DSMC%WallModel.EQ.1) THEN 
  ALLOCATE( Adsorption%MaxCoverage(1:SurfMesh%nSides,1:nSpecies),& 
            Adsorption%InitStick(1:SurfMesh%nSides,1:nSpecies),& 
            Adsorption%PrefactorStick(1:SurfMesh%nSides,1:nSpecies),& 
            Adsorption%Adsorbexp(1:SurfMesh%nSides,1:nSpecies),& 
            Adsorption%Nu_a(1:SurfMesh%nSides,1:nSpecies),& 
            Adsorption%Nu_b(1:SurfMesh%nSides,1:nSpecies),& 
            Adsorption%DesorbEnergy(1:SurfMesh%nSides,1:nSpecies),& 
            Adsorption%Intensification(1:SurfMesh%nSides,1:nSpecies))
ELSE IF (DSMC%WallModel.EQ.3) THEN 
  ALLOCATE( Adsorption%HeatOfAdsZero(1:nPartBound,1:nSpecies),&
            Adsorption%Coordination(1:nPartBound,1:nSpecies),&
            Adsorption%DiCoord(1:nPartBound,1:nSpecies))
  ! initialize info and constants
END IF
DO iSpec = 1,nSpecies
#if (PP_TimeDiscMethod==42)
  Adsorption%AdsorpInfo(iSpec)%MeanProbAds = 0.
  Adsorption%AdsorpInfo(iSpec)%MeanProbDiss = 0.
  Adsorption%AdsorpInfo(iSpec)%MeanProbDes = 0.
  Adsorption%AdsorpInfo(iSpec)%MeanEads = 0.
  Adsorption%AdsorpInfo(iSpec)%WallCollCount = 0
  Adsorption%AdsorpInfo(iSpec)%WallSpecNumCount = 0
  Adsorption%AdsorpInfo(iSpec)%NumOfAds = 0
  Adsorption%AdsorpInfo(iSpec)%NumOfDes = 0
  Adsorption%AdsorpInfo(iSpec)%Accomodation = 0
#endif
  WRITE(UNIT=hilf,FMT='(I2)') iSpec
  IF (DSMC%WallModel.EQ.1) THEN
    Adsorption%MaxCoverage(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-MaximumCoverage','0.')
    Adsorption%InitStick(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-InitialStick','0.')
    Adsorption%PrefactorStick(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-PrefactorStick','0.')
    Adsorption%Adsorbexp(:,iSpec) = GETINT('Part-Species'//TRIM(hilf)//'-Adsorbexp','1')
    Adsorption%Nu_a(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Nu-a','0.')
    Adsorption%Nu_b(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Nu-b','0.')
    Adsorption%DesorbEnergy(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Desorption-Energy-K','1.')
    Adsorption%Intensification(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Intensification-K','0.')
  ELSE IF (DSMC%WallModel.EQ.3) THEN 
    DO iPartBound=1,nPartBound
      IF(PartBound%SolidCatalytic(iPartBound))THEN
        WRITE(UNIT=hilf2,FMT='(I2)') iPartBound 
        hilf2=TRIM(hilf)//'-PartBound'//TRIM(hilf2)
        Adsorption%Coordination(iPartBound,iSpec) = GETINT('Part-Species'//TRIM(hilf2)//'-Coordination','0')
        Adsorption%DiCoord(iPartBound,iSpec) = GETINT('Part-Species'//TRIM(hilf2)//'-DiCoordination','0')
        Adsorption%HeatOfAdsZero(iPartbound,iSpec) = GETREAL('Part-Species'//TRIM(hilf2)//'-HeatOfAdsorption-K','0.')
        IF (Adsorption%Coordination(iPartBound,iSpec).EQ.0)THEN
        WRITE(UNIT=hilf2,FMT='(I2)') iPartBound 
          CALL abort(&
__STAMP__,&
'Coordination of Species '//TRIM(hilf)//' for catalytic particle boundary '//TRIM(hilf2)//' not defined')
        END IF
      END IF
    END DO
  END IF
END DO
! initialize temperature programmed desorption specific variables
#if (PP_TimeDiscMethod==42)
  Adsorption%TPD = GETLOGICAL('Particles-DSMC-Adsorption-doTPD','.FALSE.')
  Adsorption%TPD_beta = GETREAL('Particles-DSMC-Adsorption-TPD-Beta','0.')
  Adsorption%TPD_Temp = 0.
#endif
! allocate and initialize adsorption variables
ALLOCATE( Adsorption%Coverage(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides,1:nSpecies),&
          Adsorption%ProbAds(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides,1:nSpecies),&
          Adsorption%ProbDes(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides,1:nSpecies),&
          Adsorption%SumDesorbPart(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides,1:nSpecies),&
          Adsorption%SumReactPart(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides,1:nSpecies),&
          Adsorption%SumAdsorbPart(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides,1:nSpecies),&
          Adsorption%SurfSideToGlobSideMap(1:SurfMesh%nTotalSides),&
          Adsorption%DensSurfAtoms(1:SurfMesh%nTotalSides),&
          Adsorption%AreaIncrease(1:SurfMesh%nTotalSides),&
          Adsorption%CrystalIndx(1:SurfMesh%nTotalSides))

Adsorption%SurfSideToGlobSideMap(:) = -1
DO iSide = 1,nTotalSides
  IF (SurfMesh%SideIDToSurfID(iSide).LE.0) CYCLE
  Adsorption%SurfSideToGlobSideMap(SurfMesh%SideIDToSurfID(iSide)) = iSide
END DO
! initialize surface properties from particle boundary values
DO iSide=1,SurfMesh%nTotalSides
  SideID = Adsorption%SurfSideToGlobSideMap(iSide)
  PartboundID = PartBound%MapToPartBC(BC(SideID))
  !IF (DSMC%WallModel.EQ.3) Adsorption%SurfMassIC(iSide) = PartBound%SolidMassIC(PartBoundID)
  Adsorption%DensSurfAtoms(iSide) = PartBound%SolidPartDens(PartBoundID)
  Adsorption%AreaIncrease(iSide) = PartBound%SolidAreaIncrease(PartBoundID)
  Adsorption%CrystalIndx(iSide) = PartBound%SolidCrystalIndx(PartBoundID)
  Adsorption%DensSurfAtoms(iSide) = Adsorption%DensSurfAtoms(iSide)*Adsorption%AreaIncrease(iSide)
END DO
! initialize temporary coverage array for each particle boundary
ALLOCATE(Coverage_tmp(1:nPartBound,1:nSpecies))
Coverage_tmp(:,:) = 0.
DO iSpec = 1,nSpecies
  WRITE(UNIT=hilf,FMT='(I2)') iSpec
  DO iPartBound=1,nPartBound
    WRITE(UNIT=hilf2,FMT='(I2)') iPartBound 
    hilf2=TRIM(hilf)//'-PartBound'//TRIM(hilf2)
    Coverage_tmp(iPartBound,iSpec) = GETREAL('Part-Species'//TRIM(hilf2)//'-InitialCoverage','0.')
  END DO
END DO
! write temporary coverage values into global coverage array for each surface
DO iSide=1,SurfMesh%nSides
  SideID = Adsorption%SurfSideToGlobSideMap(iSide)
  PartboundID = PartBound%MapToPartBC(BC(SideID))
  DO iSpec=1,nSpecies
    Adsorption%Coverage(:,:,iSide,iSpec) = Coverage_tmp(PartBoundID,iSpec)
  END DO
END DO
SDEALLOCATE(Coverage_tmp)

Adsorption%ProbAds(:,:,:,:) = 0.
Adsorption%ProbDes(:,:,:,:) = 0.
Adsorption%SumDesorbPart(:,:,:,:) = 0
Adsorption%SumAdsorbPart(:,:,:,:) = 0
Adsorption%SumReactPart(:,:,:,:) = 0

#ifdef MPI
! allocate send and receive buffer
ALLOCATE(AdsorbSendBuf(SurfCOMM%nMPINeighbors))
ALLOCATE(AdsorbRecvBuf(SurfCOMM%nMPINeighbors))
DO iProc=1,SurfCOMM%nMPINeighbors
  ALLOCATE(AdsorbSendBuf(iProc)%content_int(nSpecies*(nSurfSample**2)*SurfExchange%nSidesSend(iProc)))
  ALLOCATE(AdsorbRecvBuf(iProc)%content_int(nSpecies*(nSurfSample**2)*SurfExchange%nSidesRecv(iProc)))
  AdsorbSendBuf(iProc)%content_int=0
  AdsorbRecvBuf(iProc)%content_int=0
END DO ! iProc
#endif /*MPI*/

IF (DSMC%WallModel.EQ.3) THEN
  CALL Init_SurfDist()
  CALL Init_SurfChem()
#if (PP_TimeDiscMethod==42)
! initial distribution into equilibrium distribution (needed for better TPD results)
  DO idiff=1,3
    CALL CalcDiffusion()
  END DO
#endif
END IF

END SUBROUTINE InitDSMCSurfModel

SUBROUTINE Init_SurfDist()
!===================================================================================================================================
! Initializing surface distibution Model for calculating of coverage effects on heat of adsorption
!===================================================================================================================================
  USE MOD_Globals
  USE MOD_ReadInTools
  USE MOD_Mesh_Vars,              ONLY : BC
  USE MOD_Particle_Vars,          ONLY : nSpecies, Species!, KeepWallParticles
  USE MOD_DSMC_Vars,              ONLY : Adsorption, SurfDistInfo
  USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh, PartBound
#ifdef MPI
  USE MOD_Particle_Boundary_Vars, ONLY : SurfCOMM
  USE MOD_Particle_MPI_Vars,      ONLY : SurfDistSendBuf,SurfDistRecvBuf,SurfExchange
  USE MOD_DSMC_SurfModel_Tools,   ONLY : ExchangeSurfDistInfo, ExchangeSurfDistSize
#endif /*MPI*/
!===================================================================================================================================
  IMPLICIT NONE
!=================================================================================================================================== 
! Local variable declaration
  CHARACTER(64)                    :: particle_mpf
  REAL                             :: surface_mpf
  INTEGER                          :: Max_Surfsites_num
  INTEGER                          :: Max_Surfsites_own
  INTEGER                          :: Max_Surfsites_halo
  INTEGER                          :: iSurfSide, subsurfxi, subsurfeta, iSpec, iInterAtom
  INTEGER                          :: SideID, PartBoundID
  INTEGER                          :: surfsquare, dist, Adsorbates
  INTEGER                          :: Surfpos, Surfnum, Indx, Indy, UsedSiteMapPos
  REAL                             :: RanNum
  INTEGER                          :: xpos, ypos
  INTEGER                          :: Coord, nSites, nInterAtom, nNeighbours
#ifdef MPI
  INTEGER                          :: iProc
#endif
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
ALLOCATE(SurfDistInfo(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides))
DO iSurfSide = 1,SurfMesh%nTotalSides
  DO subsurfeta = 1,nSurfSample
  DO subsurfxi = 1,nSurfSample
    ALLOCATE( SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%nSites(1:3),&
              SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%SitesRemain(1:3),&
              SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%adsorbnum_tmp(1:nSpecies),&
              SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%desorbnum_tmp(1:nSpecies),&
              SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%reactnum_tmp(1:nSpecies))
  END DO
  END DO
END DO
! IF (.NOT.KeepWallParticles) THEN
!   surfsquare = GETINT('Particles-DSMC-AdsorptionSites','10000')
!   surfsquare = INT(SQRT(REAL(surfsquare))) - 1
! END IF
WRITE(UNIT=particle_mpf,FMT='(E11.3)') Species(1)%MacroParticleFactor
surface_mpf = GETREAL('Particles-Surface-MacroParticleFactor',TRIM(particle_mpf))
Max_Surfsites_num = 0
Max_Surfsites_own = 0
Max_Surfsites_halo = 0

! Allocate and initializes number of surface sites and neighbours 
DO iSurfSide = 1,SurfMesh%nTotalSides
  DO subsurfeta = 1,nSurfSample
  DO subsurfxi = 1,nSurfSample
!     IF (KeepWallParticles) THEN ! does not work with vMPF
!       surfsquare = INT(Adsorption%DensSurfAtoms(iSurfSide) &
!                     * SurfMesh%SurfaceArea(subsurfxi,subsurfeta,iSurfSide) &
!                     / Species(1)%MacroParticleFactor)
!       surfsquare = INT(SQRT(REAL(surfsquare))) - 1
!     END IF
    surfsquare = INT(Adsorption%DensSurfAtoms(iSurfSide) &
                  * SurfMesh%SurfaceArea(subsurfxi,subsurfeta,iSurfSide) &
                  / surface_mpf)
    surfsquare = INT(SQRT(REAL(surfsquare))) - 1
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%nSites(1) = INT(surfsquare**2)
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%nSites(2) = INT( 2*(surfsquare*(surfsquare+1)) )
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%nSites(3) = INT((surfsquare+1)**2)
    
    Max_Surfsites_num = Max_Surfsites_num + SUM(SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%nSites(:))
    IF (iSurfSide.LE.SurfMesh%nSides) THEN
      Max_Surfsites_own = Max_Surfsites_own + SUM(SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%nSites(:))
    ELSE
      Max_Surfsites_halo = Max_Surfsites_halo + SUM(SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%nSites(:))
    END IF
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%SitesRemain(:) = SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%nSites(:)
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%adsorbnum_tmp(1:nSpecies) = 0.
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%desorbnum_tmp(1:nSpecies) = 0.
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%reactnum_tmp(1:nSpecies) = 0.
    
    ALLOCATE( SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%SurfAtomBondOrder(1:nSpecies,1:surfsquare+1,1:surfsquare+1))
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%SurfAtomBondOrder(:,:,:) = 0
    
    ALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1:3))
    DO Coord = 1,3    
      SELECT CASE (Coord)
      CASE(1)
        nSites = SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%nSites(Coord)
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%nInterAtom = 4
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%nNeighbours = 16
      CASE(2)
        nSites = SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%nSites(Coord)
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%nInterAtom = 2
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%nNeighbours = 14
      CASE(3)
        nSites = SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%nSites(Coord)
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%nInterAtom = 1
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%nNeighbours = 12
      END SELECT
        
      nInterAtom = SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%nInterAtom
      nNeighbours = SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%nNeighbours
      ALLOCATE( SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%BondAtomIndx(1:nSites,nInterAtom),&
                SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%BondAtomIndy(1:nSites,nInterAtom),&
                SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%NeighPos(1:nSites,1:nNeighbours),&
                SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%NeighSite(1:nSites,1:nNeighbours),&
                SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%IsNearestNeigh(1:nSites,1:nNeighbours))
      ALLOCATE( SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%UsedSiteMap(1:nSites),&
                SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%Species(1:nSites))
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%UsedSiteMap(:) = 0
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%Species(:) = 0
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%BondAtomIndx(:,:) = 0
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%BondAtomIndy(:,:) = 0
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%NeighPos(:,:) = 0
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%NeighSite(:,:) = 0
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%IsNearestNeigh(:,:) = .FALSE.
    END DO
  END DO
  END DO
END DO

DO iSurfSide = 1,SurfMesh%nTotalSides
DO subsurfeta = 1,nSurfSample
DO subsurfxi = 1,nSurfSample
  ! surfsquare chosen from nSite(1) for correct SurfIndx definitions
  surfsquare = INT(SQRT(REAL(SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%nSites(1))))
  ! allocate and define surface indexes for adsorbate distribution and build mapping of respective bondatoms and neighbours      
  Indx = 1
  Indy = 1
  DO Surfpos = 1,SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%nSites(1)
    IF (Indx.GT.surfsquare) THEN
      Indx = 1
      Indy = Indy + 1
    END IF
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%UsedSiteMap(Surfpos) = Surfpos
    ! mapping respective neighbours first hollow then bridge then top
    ! hollow
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,1) = Surfpos - surfsquare - 1
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,2) = Surfpos - surfsquare
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,3) = Surfpos - surfsquare + 1
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,4) = Surfpos - 1
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,5) = Surfpos + 1
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,6) = Surfpos + surfsquare - 1
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,7) = Surfpos + surfsquare
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,8) = Surfpos + surfsquare + 1
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,1:8) = 1
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,2) = .TRUE.
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,4) = .TRUE.
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,5) = .TRUE.
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,7) = .TRUE.
    ! bridge
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,9) = Surfpos +(surfsquare+1)*(Indy-1)
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,10) = Surfpos +surfsquare +(surfsquare+1)*(Indy-1)
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,11) = Surfpos +surfsquare +(surfsquare+1)*(Indy-1) +1
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,12) = Surfpos +surfsquare +(surfsquare+1)*(Indy)
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,9:12) = 2
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,9:12) = .TRUE.
    ! top
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,13) = Surfpos + (Indy-1)
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,14) = Surfpos + 1 + (Indy-1)
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,15) = Surfpos + (surfsquare+1) + (Indy-1)
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,16) = Surfpos + (surfsquare+1) + 1 + (Indy-1)
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,13:16) = 3
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,13:16) = .TRUE.
    ! account for empty edges
    IF (Indy .EQ. 1) SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,1:3) = 0
    IF (Indy .EQ. surfsquare) SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,6:8) = 0
    IF (Indx .EQ. 1) THEN
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,1) = 0
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,4) = 0
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,6) = 0
    END IF
    IF (Indx .EQ. surfsquare) THEN
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,3) = 0
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,5) = 0
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,8) = 0
    END IF
    ! mapping respective bond atoms
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%BondAtomIndx(Surfpos,1) = Indx
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%BondAtomIndy(Surfpos,1) = Indy
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%BondAtomIndx(Surfpos,2) = Indx+1
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%BondAtomIndy(Surfpos,2) = Indy
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%BondAtomIndx(Surfpos,3) = Indx
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%BondAtomIndy(Surfpos,3) = Indy+1
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%BondAtomIndx(Surfpos,4) = Indx+1
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(1)%BondAtomIndy(Surfpos,4) = Indy+1
    Indx = Indx + 1
  END DO
  Indx = 1
  Indy = 1
  DO Surfpos = 1,SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%nSites(2)
    IF (Indx.GT.(2*surfsquare+1)) THEN
      Indx = 1
      Indy = Indy + 1
    END IF
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%UsedSiteMap(Surfpos) = Surfpos
    IF (Indx .LE. surfsquare) THEN ! surface atoms are LEFT an RIGHT of adsorbate site
      ! mapping respective neighbours first hollow then bridge then top
      ! hollow
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,1) = Surfpos-surfsquare -(surfsquare+1)*(Indy-1) -1
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,2) = Surfpos-surfsquare -(surfsquare+1)*(Indy-1)
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,3) = Surfpos-surfsquare -(surfsquare+1)*(Indy-1) +1
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,4) = Surfpos -(surfsquare+1)*(Indy-1) -1
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,5) = Surfpos -(surfsquare+1)*(Indy-1)
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,6) = Surfpos -(surfsquare+1)*(Indy-1) +1
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,1:6) = 1
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,2) = .TRUE.
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,5) = .TRUE.
      ! bridge
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,7) = Surfpos - (surfsquare+1)
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,8) = Surfpos - (surfsquare+1) + 1
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,9) = Surfpos - 1
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,10) = Surfpos + 1
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,11) = Surfpos + surfsquare
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,12) = Surfpos + surfsquare + 1
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,7:12) = 2
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,8) = .TRUE.
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,11) = .TRUE.
      ! top
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,13) = Surfpos -(surfsquare)*(Indy-1)
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,14) = Surfpos +1 -(surfsquare)*(Indy-1)
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,13:14) = 3
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,13:14) = .TRUE.
      ! account for empty edges
      IF (Indy .EQ. 1) THEN
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,1:3) = 0
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,7:8) = 0
      END IF
      IF (Indy .EQ. surfsquare+1) THEN
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,4:6) = 0
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,11:12) = 0
      END IF
      IF (Indx .EQ. 1) THEN
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,1) = 0
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,4) = 0
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,9) = 0
      END IF
      IF (Indx .EQ. surfsquare) THEN
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,3) = 0
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,6) = 0
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,10) = 0
      END IF
      ! mapping respective bond atoms
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%BondAtomIndx(Surfpos,1) = Indx
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%BondAtomIndy(Surfpos,1) = Indy
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%BondAtomIndx(Surfpos,2) = Indx+1
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%BondAtomIndy(Surfpos,2) = Indy
    ELSE ! surface atoms are TOP and DOWN of adsorbate site
      ! mapping respective neighbours first hollow then bridge then top
      ! hollow
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,1) = Surfpos -(2*surfsquare) &
                                                                                  -(surfsquare+1)*(Indy-1) -1
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,2) = Surfpos -(2*surfsquare)&
                                                                                  -(surfsquare+1)*(Indy-1)
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,3) = Surfpos -surfsquare &
                                                                                  -(surfsquare+1)*(Indy-1) -1
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,4) = Surfpos - surfsquare -(surfsquare+1)*(Indy-1)
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,5) = Surfpos -(surfsquare+1)*(Indy-1) -1
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,6) = Surfpos -(surfsquare+1)*(Indy-1)
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,1:6) = 1
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,3:4) = .TRUE.
      ! bridge
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,7) = Surfpos - surfsquare - (surfsquare+1)
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,8) = Surfpos - surfsquare - 1
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,9) = Surfpos - surfsquare
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,10) = Surfpos + (surfsquare+1) - 1
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,11) = Surfpos + (surfsquare+1)
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,12) = Surfpos + surfsquare + (surfsquare+1)
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,7:12) = 2
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,8:11) = .TRUE.
      ! top
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,13) = Surfpos -surfsquare*(Indy)
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,14) = Surfpos -surfsquare*(Indy) +(surfsquare+1)
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,13:14) = 3
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,13:14) = .TRUE.
      ! account for empty edges
      IF (Indy .EQ. 1) THEN
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,1:2) = 0
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,7) = 0
      END IF
      IF (Indy .EQ. surfsquare) THEN
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,5:6) = 0
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,12) = 0
      END IF
      IF (Indx .EQ. (surfsquare+1)) THEN
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,1) = 0
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,3) = 0
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,5) = 0
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,8) = 0
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,10) = 0
      END IF
      IF (Indx .EQ. 2*surfsquare+1) THEN
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,2) = 0
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,4) = 0
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,6) = 0
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,9) = 0
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,11) = 0
      END IF
      ! mapping respective bond atoms
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%BondAtomIndx(Surfpos,1) = Indx - surfsquare
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%BondAtomIndy(Surfpos,1) = Indy 
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%BondAtomIndx(Surfpos,2) = Indx - surfsquare
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(2)%BondAtomIndy(Surfpos,2) = Indy+1
    END IF
    Indx = Indx + 1
  END DO
  Indx = 1
  Indy = 1
  DO Surfpos = 1,SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%nSites(3)
    IF (Indx.GT.surfsquare+1) THEN
      Indx = 1
      Indy = Indy + 1
    END IF
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%UsedSiteMap(Surfpos) = Surfpos
    ! mapping respective neighbours first hollow then bridge then top
    ! hollow
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,1) = Surfpos - (surfsquare) - 1 -(Indy-1)
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,2) = Surfpos - (surfsquare) - (Indy-1)
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,3) = Surfpos - 1 - (Indy-1)
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,4) = Surfpos - (Indy-1)    
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,1:4) = 1
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%IsNearestNeigh(Surfpos,1:4) = .TRUE.
    ! bridge
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,5) = Surfpos + surfsquare*(Indy-1) -(surfsquare+1)
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,6) = Surfpos -1 +(surfsquare)*(Indy-1)
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,7) = Surfpos +(surfsquare)*(Indy-1)
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,8) = Surfpos + surfsquare*(Indy)
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,5:8) = 2
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%IsNearestNeigh(Surfpos,5:8) = .TRUE.
    ! top
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,9) = Surfpos - (surfsquare+1)
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,10) = Surfpos - 1
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,11) = Surfpos + 1
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,12) = Surfpos + (surfsquare+1)
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,9:12) = 3
    ! account for empty edges
    IF (Indy .EQ. 1) THEN
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,1:2) = 0
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,5) = 0
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,9) = 0
    END IF
    IF (Indx .EQ. 1) THEN
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,1) = 0
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,3) = 0
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,6) = 0
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,10) = 0
    END IF
    IF (Indy .EQ. (surfsquare+1)) THEN
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,3:4) = 0
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,8) = 0
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,12) = 0
    END IF
    IF (Indx .EQ. (surfsquare+1)) THEN
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,2) = 0
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,4) = 0
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,7) = 0
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,11) = 0
    END IF
    ! mapping respective bond atoms
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%BondAtomIndx(Surfpos,1) = Indx
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(3)%BondAtomIndy(Surfpos,1) = Indy
    Indx = Indx + 1
  END DO
  
END DO
END DO
END DO

DO iSurfSide = 1,SurfMesh%nSides
SideID = Adsorption%SurfSideToGlobSideMap(iSurfSide)
PartboundID = PartBound%MapToPartBC(BC(SideID))
DO subsurfeta = 1,nSurfSample
DO subsurfxi = 1,nSurfSample    
  DO iSpec = 1,nSpecies
    ! adjust coverage to actual discret value
    Adsorbates = INT(Adsorption%Coverage(subsurfxi,subsurfeta,iSurfSide,iSpec) &
                * SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%nSites(Adsorption%Coordination(PartboundID,iSpec)))
    Adsorption%Coverage(subsurfxi,subsurfeta,iSurfSide,iSpec) = REAL(Adsorbates) &
        / REAL(SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%nSites(3))
    IF (SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%SitesRemain(Adsorption%Coordination(PartboundID,iSpec)).LT.Adsorbates) THEN
      CALL abort(&
__STAMP__&
,'Error in Init_SurfDist: Too many Adsorbates! - Choose lower Coverages for coordination:', &
Adsorption%Coordination(PartboundID,iSpec))
    END IF
    ! distribute adsorbates randomly on the surface on the correct site and assign surface atom bond order
    dist = 1
    Coord = Adsorption%Coordination(PartboundID,iSpec)
    Surfnum = SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%SitesRemain(Coord)
    DO WHILE (dist.LE.Adsorbates) 
      CALL RANDOM_NUMBER(RanNum)
      Surfpos = 1 + INT(Surfnum * RanNum)
      UsedSiteMapPos = SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%UsedSiteMap(Surfpos)
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%Species(UsedSiteMapPos) = iSpec
      ! assign bond order of respective surface atoms in the surfacelattice
      DO iInterAtom = 1,SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%nInterAtom
        xpos = SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%BondAtomIndx(UsedSiteMapPos,iInterAtom)
        ypos = SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%BondAtomIndy(UsedSiteMapPos,iInterAtom)
        SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%SurfAtomBondOrder(iSpec,xpos,ypos) = &
          SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%SurfAtomBondOrder(iSpec,xpos,ypos) + 1
      END DO
      ! rearrange UsedSiteMap-Surfpos-array
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%UsedSiteMap(Surfpos) = &
          SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%UsedSiteMap(Surfnum)
      SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(Coord)%UsedSiteMap(Surfnum) = UsedSiteMapPos
      Surfnum = Surfnum - 1
      dist = dist + 1
    END DO
    SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%SitesRemain(Coord) = Surfnum
  END DO
END DO
END DO
END DO


#ifdef MPI
#ifdef CODE_ANALYZE
! write out the number of sites on all surface of the proc, that are considered for adsorption
WRITE(UNIT_stdOut,'(A,I3,I13,A,I13,A,I13)')' | Maximum number of surface sites on proc: ',myRank,Max_Surfsites_num,&
  ' | own: ',Max_Surfsites_own,' | halo: ',Max_Surfsites_halo
#endif /*CODE_ANALYZE*/
  

ALLOCATE(SurfExchange%nSurfDistSidesSend(1:SurfCOMM%nMPINeighbors) &
        ,SurfExchange%nSurfDistSidesRecv(1:SurfCOMM%nMPINeighbors) &
        ,SurfExchange%SurfDistSendRequest(1:2,1:SurfCOMM%nMPINeighbors)  &
        ,SurfExchange%SurfDistRecvRequest(1:2,1:SurfCOMM%nMPINeighbors)  &
        ,SurfExchange%NbrOfPos(1:SurfCOMM%nMPINeighbors))
! allocate send and receive buffer
ALLOCATE(SurfDistSendBuf(SurfCOMM%nMPINeighbors))
ALLOCATE(SurfDistRecvBuf(SurfCOMM%nMPINeighbors))
DO iProc=1,SurfCOMM%nMPINeighbors
  SurfExchange%nSurfDistSidesSend(iProc) = SurfExchange%nSidesRecv(iProc)
  SurfExchange%nSurfDistSidesRecv(iProc) = SurfExchange%nSidesSend(iProc)
  ALLOCATE(SurfCOMM%MPINeighbor(iProc)%SurfDistSendList(SurfExchange%nSurfDistSidesSend(iProc)))
  ALLOCATE(SurfCOMM%MPINeighbor(iProc)%SurfDistRecvList(SurfExchange%nSurfDistSidesRecv(iProc)))
  SurfCOMM%MPINeighbor(iProc)%SurfDistRecvList(:)=SurfCOMM%MPINeighbor(iProc)%SendList(:)
  SurfCOMM%MPINeighbor(iProc)%SurfDistSendList(:)=SurfCOMM%MPINeighbor(iProc)%RecvList(:)
  ALLOCATE(SurfExchange%NbrOfPos(iProc)%nPosSend(1:SurfExchange%nSurfDistSidesSend(iProc)))
  ALLOCATE(SurfExchange%NbrOfPos(iProc)%nPosRecv(1:SurfExchange%nSurfDistSidesRecv(iProc)))
  SurfExchange%NbrOfPos(iProc)%nPosSend = 0
  SurfExchange%NbrOfPos(iProc)%nPosRecv = 0
END DO ! iProc

! communicate the number of surface sites for surfdist communication
CALL ExchangeSurfDistSize()
! fill halo surface distribution through mpi communication
CALL ExchangeSurfDistInfo()
#endif /*MPI*/

SWRITE(UNIT_stdOut,'(A)')' INIT SURFACE DISTRIBUTION DONE!'
    
END SUBROUTINE Init_SurfDist

SUBROUTINE Init_SurfChem()
!===================================================================================================================================
! Initializing surface reaction variables
!===================================================================================================================================
! MODULES
USE MOD_Globals,                ONLY : abort, MPIRoot, UNIT_StdOut
USE MOD_DSMC_Vars,              ONLY : Adsorption, SpecDSMC
USE MOD_PARTICLE_Vars,          ONLY : nSpecies
USE MOD_ReadInTools
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  CHARACTER(32)                    :: hilf, hilf2
  INTEGER                          :: iSpec, iSpec2, iReactNum, iReactNum2, iReactant
  INTEGER                          :: ReactNum
  INTEGER                          :: MaxDissNum, MaxReactNum, MaxAssocNum
  INTEGER , ALLOCATABLE            :: nAssocReact(:)
  INTEGER                          :: nDissoc, nDisProp
  INTEGER                          :: CalcTST_Case
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)')' INIT SURFACE CHEMISTRY...'

Adsorption%NumOfDissocReact = 0
Adsorption%NumOfAssocReact = 0
Adsorption%NumOfExchReact = 0

! Adsorption constants
ALLOCATE( Adsorption%Ads_Powerfactor(1:nSpecies),&
          Adsorption%Ads_Prefactor(1:nSpecies))!,&
!           Adsorption%ER_Powerfactor(1:nSpecies),&
!           Adsorption%ER_Prefactor(1:nSpecies))
DO iSpec = 1,nSpecies            
  WRITE(UNIT=hilf,FMT='(I2)') iSpec
  Adsorption%Ads_Powerfactor(iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Adsorption-Powerfactor','0.')
  Adsorption%Ads_Prefactor(iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Adsorption-Prefactor','0.')
  !Adsorption%ER_Powerfactor(iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-ER-Adsorption-Powerfactor','0.') 
  !Adsorption%ER_Prefactor(iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-ER-Adsorption-Prefactor','0.')
END DO

MaxDissNum = GETINT('Part-Species-MaxDissNum','0')
MaxAssocNum = MaxDissNum

! allocate and initialize dissociative and associative reactions species map
IF ( (MaxDissNum.GT.0) .OR. (MaxAssocNum.GT.0) ) THEN
  ALLOCATE( Adsorption%DissocReact(1:2,1:MaxDissNum,1:nSpecies),&
            Adsorption%Diss_Powerfactor(1:MaxDissNum,1:nSpecies),&
            Adsorption%Diss_Prefactor(1:MaxDissNum,1:nSpecies))
  ! Read in dissociative reactions and respective dissociation bond energies
  DO iSpec = 1,nSpecies            
    WRITE(UNIT=hilf,FMT='(I2)') iSpec
    DO iReactNum = 1,MaxDissNum
      WRITE(UNIT=hilf2,FMT='(I2)') iReactNum
      Adsorption%DissocReact(:,iReactNum,iSpec) = &
                                       GETINTARRAY('Part-Species'//TRIM(hilf)//'-SurfDiss'//TRIM(hilf2)//'-Products',2,'0,0')
      IF ((Adsorption%DissocReact(1,iReactNum,iSpec).GT.nSpecies).OR.(Adsorption%DissocReact(2,iReactNum,iSpec).GT.nSpecies) ) THEN
        CALL abort(&
__STAMP__&
,'Error in Init_SurfChem: Product species for reaction '//TRIM(hilf2)//' not defined!')
      END IF
      Adsorption%Diss_Powerfactor(iReactNum,iSpec) = &
                                       GETREAL('Part-Species'//TRIM(hilf)//'-SurfDiss'//TRIM(hilf2)//'-Powerfactor','0.')
      Adsorption%Diss_Prefactor(iReactNum,iSpec) = &
                                       GETREAL('Part-Species'//TRIM(hilf)//'-SurfDiss'//TRIM(hilf2)//'-Prefactor','0.')
    END DO
  END DO
  
  ! find max number of associative reactions for each species from dissociations
  ALLOCATE(nAssocReact(1:nSpecies))
  nAssocReact(:) = 0
  DO iSpec = 1,nSpecies
    DO iSpec2 = 1,nSpecies
    DO iReactNum = 1,MaxDissNum
      IF ((Adsorption%DissocReact(1,iReactNum,iSpec2).EQ.iSpec).OR.(Adsorption%DissocReact(2,iReactNum,iSpec2).EQ.iSpec) ) THEN
        nAssocReact(iSpec) = nAssocReact(iSpec) + 1
        Adsorption%NumOfAssocReact = Adsorption%NumOfDissocReact + 1
      END IF
    END DO
    END DO
  END DO
  MaxAssocNum = MAXVAL(nAssocReact)
  DEALLOCATE(nAssocReact)
  
  ! fill associative reactions species map from defined dissociative reactions
  MaxReactNum = MaxDissNum + MaxAssocNum
  ALLOCATE( Adsorption%AssocReact(1:2,1:MaxAssocNum,1:nSpecies),&
            Adsorption%EDissBond(0:MaxReactNum,1:nSpecies),&
            Adsorption%EDissBondAdsorbPoly(0:1,1:nSpecies))
  Adsorption%EDissBond(0:MaxReactNum,1:nSpecies) = 0.
  Adsorption%EDissBondAdsorbPoly(0:1,1:nSpecies) = 0.
  DO iSpec = 1,nSpecies            
    WRITE(UNIT=hilf,FMT='(I2)') iSpec
    DO iReactNum = 1,MaxDissNum
      WRITE(UNIT=hilf2,FMT='(I2)') iReactNum
      Adsorption%EDissBond(iReactNum,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-SurfDiss'//TRIM(hilf2)//'-EDissBond','0.')
    END DO
  END DO
  DO iSpec = 1,nSpecies            
    WRITE(UNIT=hilf,FMT='(I2)') iSpec
    IF (SpecDSMC(iSpec)%InterID.EQ.2) THEN
      Adsorption%EDissBond(0,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Adsorption-EDissBond','0.')
      IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
        Adsorption%EDissBondAdsorbPoly(0,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Adsorption-EDissBondPoly1','0.')
        Adsorption%EDissBondAdsorbPoly(1,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Adsorption-EDissBondPoly2','0.')
        IF (( MAXVAL(Adsorption%DiCoord(:,iSpec)).NE.0) .AND. (Adsorption%EDissBondAdsorbPoly(0,iSpec).EQ.0)) THEN
          CALL abort(&
__STAMP__&
,'Error in Init_SurfChem: Adsorption dissociation bond energy of dicoord for species '//TRIM(hilf)//' not defined!')
        END IF
      END IF
      IF ((.NOT.SpecDSMC(iSpec)%PolyatomicMol).AND.(MAXVAL(Adsorption%DiCoord(:,iSpec)).EQ.0) &
          .AND.(Adsorption%EDissBond(0,iSpec).EQ.0.))THEN
        CALL abort(&
__STAMP__&
,'Error in Init_SurfChem: Adsorption dissociation bond energy for species '//TRIM(hilf)//' not defined!')
      END IF
    END IF
  END DO
  DO iSpec = 1,nSpecies
    ReactNum = 1
    DO iSpec2 = 1,nSpecies
    DO iReactNum2 = 1,MaxDissNum
      IF (Adsorption%DissocReact(1,iReactNum2,iSpec2).EQ.iSpec) THEN
        Adsorption%AssocReact(1,ReactNum,iSpec) = Adsorption%DissocReact(2,iReactNum2,iSpec2)
        Adsorption%AssocReact(2,ReactNum,iSpec) = iSpec2
        Adsorption%EDissBond((MaxDissNum+ReactNum),iSpec) = Adsorption%EDissBond(iReactNum2,iSpec2)
        ReactNum = ReactNum + 1
      ELSE IF (Adsorption%DissocReact(2,iReactNum2,iSpec2).EQ.iSpec) THEN
        Adsorption%AssocReact(1,ReactNum,iSpec) = Adsorption%DissocReact(1,iReactNum2,iSpec2)
        Adsorption%AssocReact(2,ReactNum,iSpec) = iSpec2
        Adsorption%EDissBond((MaxDissNum+ReactNum),iSpec) = Adsorption%EDissBond(iReactNum2,iSpec2)
        ReactNum = ReactNum + 1
      ELSE
        CYCLE
      END IF
    END DO
    END DO
    IF (ReactNum.LE.(MaxAssocNum)) THEN
      Adsorption%AssocReact(:,ReactNum:(MaxReactNum-MaxDissNum),iSpec) = 0
    END IF
  END DO
  
  nDissoc = GETINT('Part-SurfChem-Nbr-DissocReactions','0')
  IF ((nDissoc.GT.0) .AND. (nDissoc.NE.nSpecies)) THEN
    CALL abort(&
__STAMP__&
,'Error in Init_SurfChem: given number of dissociation reactions in INI-File differs of number from given parameters!')
  END IF
  DO iSpec=1,nSpecies
    DO iReactNum=1,MaxDissNum
      IF (Adsorption%DissocReact(1,iReactNum,iSpec).NE.0)THEN
        Adsorption%NumOfDissocReact = Adsorption%NumOfDissocReact + 1
      END IF
    END DO
  END DO
  nDissoc =  Adsorption%NumOfDissocReact
  nDisProp = GETINT('Part-SurfChem-Nbr-ExchangeReactions','0')
  ! Allocate and fill one array for all types of reactions (dissociation, association, disproportionation/exchange reaction)
  ALLOCATE( Adsorption%ChemReactant(1:2,1:nDissoc+nDisProp),&
            Adsorption%ChemProduct(1:2,1:nDissoc+nDisProp),&
            Adsorption%Reactant_DissBond_K(1:2,1:nDissoc+nDisProp),&
            Adsorption%Product_DissBond_K(1:2,1:nDissoc+nDisProp))
  ! Initialize allocated variables
  Adsorption%ChemReactant(1:2,1:nDissoc+nDisProp)=0
  Adsorption%ChemProduct(1:2,1:nDissoc+nDisProp)=0
  Adsorption%Reactant_DissBond_K(1:2,1:nDissoc+nDisProp)=0.
  Adsorption%Product_DissBond_K(1:2,1:nDissoc+nDisProp)=0.
  ! fill dissociation reactions (can also be used for association)
  ReactNum = 0
  DO iSpec=1,nSpecies
    DO iReactNum=1,MaxDissNum
      IF (Adsorption%DissocReact(1,iReactNum,iSpec).NE.0)THEN
        ReactNum = ReactNum + 1
        Adsorption%ChemReactant(1,ReactNum) = iSpec
        Adsorption%Reactant_DissBond_K(1,ReactNum) = Adsorption%EDissBond(iReactNum,iSpec)
        Adsorption%ChemProduct(1,ReactNum) = Adsorption%DissocReact(1,iReactNum,iSpec)
        Adsorption%ChemProduct(2,ReactNum) = Adsorption%DissocReact(2,iReactNum,iSpec)
!         DO iReactant=1,2
!           IF (SpecDSMC(Adsorption%ChemProduct(iReactant,ReactNum))%InterID.EQ.2) THEN
!             IF(SpecDSMC(Adsorption%ChemProduct(iReactant,ReactNum))%PolyatomicMol) THEN
!               !---------------------------------------------------------------------------------------------------------------------
!               SELECT CASE(MaxDissNum)
!               !---------------------------------------------------------------------------------------------------------------------
!               CASE(1) !only one possible dissociation per species
!               !---------------------------------------------------------------------------------------------------------------------
!                 Adsorption%Product_DissBond_K(iReactant,ReactNum) = &
!                         Adsorption%EDissBond(1,Adsorption%ChemProduct(iReactant,ReactNum))
!               !---------------------------------------------------------------------------------------------------------------------
!               CASE DEFAULT !more than one dissociation possible per species (special case for some polyatomic)
!               !---------------------------------------------------------------------------------------------------------------------
!                 DO iReactNum2=1,MaxDissNum
!                   IF (Adsorption%EDissBond(iReactNum2,Adsorption%ChemProduct(iReactant,ReactNum)).LE.0.) CYCLE
!                   IF ( (Adsorption%Product_DissBond_K(iReactant,ReactNum).GT.&
!                         Adsorption%EDissBond(iReactNum2,Adsorption%ChemProduct(iReactant,ReactNum))) .AND. &
!                        (Adsorption%Product_DissBond_K(iReactant,ReactNum).GT.0.) ) THEN
!                     Adsorption%Product_DissBond_K(iReactant,ReactNum) = &
!                             Adsorption%EDissBond(iReactNum2,Adsorption%ChemProduct(iReactant,ReactNum))
!                   END IF
!                 END DO
!               END SELECT
!             ELSE
!               Adsorption%Product_DissBond_K(iReactant,ReactNum) = Adsorption%EDissBond(1,Adsorption%ChemProduct(iReactant,ReactNum))
!             END IF
!           END IF
!         END DO
      END IF
    END DO
  END DO
  ! fill disproportionation reactions (fancy stuff)
  DO iReactNum = 1,nDisProp
    WRITE(UNIT=hilf,FMT='(I2)') iReactNum
    Adsorption%ChemReactant(:,iReactNum+nDissoc) = &
                                      GETINTARRAY('Part-SurfChem-Disprop'//TRIM(hilf)//'-Reactants',2,'0,0')
    Adsorption%ChemProduct(:,iReactNum+nDissoc) = &
                                      GETINTARRAY('Part-SurfChem-Disprop'//TRIM(hilf)//'-Products',2,'0,0')
    ! Error output if species in reaction not defined
    IF ((Adsorption%ChemReactant(1,iReactNum+nDissoc).GT.nSpecies).OR.&
        (Adsorption%ChemReactant(2,iReactNum+nDissoc).GT.nSpecies).OR.&
        (Adsorption%ChemProduct(1,iReactNum+nDissoc).GT.nSpecies).OR.&
        (Adsorption%ChemProduct(2,iReactNum+nDissoc).GT.nSpecies) ) THEN
      CALL abort(&
__STAMP__&
,'Error in Init_SurfChem: one reaction species for disproportionation reaction '//TRIM(hilf)//' not defined!')
    END IF
    IF ((Adsorption%ChemProduct(1,iReactNum+nDissoc).EQ.0).OR.&
        (Adsorption%ChemProduct(2,iReactNum+nDissoc).EQ.0) ) THEN
      CALL abort(&
__STAMP__&
,'Error in Init_SurfChem: Incorrect definition of disproportionation reaction '//TRIM(hilf)//'!')
    END IF
    ! Check if diatomic / polyatomic species in reactants and products have dissociation defined
    ! if not exit with error
    DO iReactant=1,2
      IF ( (SpecDSMC(Adsorption%ChemReactant(iReactant,iReactNum+nDissoc))%InterID.EQ.2) .AND. &
            (Adsorption%DissocReact(1,1,Adsorption%ChemReactant(iReactant,iReactNum+nDissoc)).EQ.0) ) THEN
      WRITE(UNIT=hilf2,FMT='(I2)') Adsorption%ChemReactant(iReactant,iReactNum+nDissoc)
        CALL abort(&
__STAMP__&
,'Error in Init_SurfChem Disproportionation: Dissociation for reactant species '//TRIM(hilf2)//' not defined!')
      END IF
      IF (SpecDSMC(Adsorption%ChemProduct(iReactant,iReactNum+nDissoc))%InterID.EQ.2 .AND. &
            (Adsorption%DissocReact(1,1,Adsorption%ChemProduct(iReactant,iReactNum+nDissoc)).EQ.0) ) THEN
      WRITE(UNIT=hilf2,FMT='(I2)') Adsorption%ChemProduct(iReactant,iReactNum+nDissoc)
        CALL abort(&
__STAMP__&
,'Error in Init_SurfChem Disproportionation: Dissociation for product species '//TRIM(hilf2)//' not defined!')
      END IF
    END DO
    ! Read dissociation bond energies of reactants and products
    Adsorption%Reactant_DissBond_K(:,iReactNum+nDissoc) = &
            GETREALARRAY('Part-SurfChem-Disprop'//TRIM(hilf)//'-DissBond_K-Reactants',2,'0.,0.')
    Adsorption%Product_DissBond_K(:,iReactNum+nDissoc) = &
            GETREALARRAY('Part-SurfChem-Disprop'//TRIM(hilf)//'-DissBond_K-Products',2,'0.,0.')
    ! Check if dissociation bond energies of reactants and products are all defined if they are at least diatomic
    ! If they are not defined take them from the appropriate dissociation reactions
    DO iReactant=1,2
      IF ( (SpecDSMC(Adsorption%ChemReactant(iReactant,iReactNum+nDissoc))%InterID.EQ.2) .AND. &
            (Adsorption%Reactant_DissBond_K(iReactant,iReactNum+nDissoc).EQ.0.) ) THEN
        IF(SpecDSMC(Adsorption%ChemReactant(iReactant,iReactNum+nDissoc))%PolyatomicMol) THEN
          !-------------------------------------------------------------------------------------------------------------------------
          SELECT CASE(MaxDissNum)
          !-------------------------------------------------------------------------------------------------------------------------
          CASE(1) !only one possible dissociation per species
          !-------------------------------------------------------------------------------------------------------------------------
            Adsorption%Reactant_DissBond_K(iReactant,iReactNum+nDissoc) = &
                    Adsorption%EDissBond(1,Adsorption%ChemReactant(iReactant,iReactNum+nDissoc))
          !-------------------------------------------------------------------------------------------------------------------------
          CASE DEFAULT !more than one dissociation possible per species (special case for some polyatomic)
          !-------------------------------------------------------------------------------------------------------------------------
          WRITE(UNIT=hilf2,FMT='(I2)') iReactNum
            CALL abort(&
__STAMP__&
,'Error in Init_SurfChem: Dissocation bond energy in Disproportionation reaction'//TRIM(hilf2)//' not defined!')
          END SELECT
        ELSE
          Adsorption%Reactant_DissBond_K(iReactant,iReactNum+nDissoc) = &
                  Adsorption%EDissBond(1,Adsorption%ChemReactant(iReactant,iReactNum+nDissoc))
        END IF
      END IF
      IF (SpecDSMC(Adsorption%ChemProduct(iReactant,iReactNum+nDissoc))%InterID.EQ.2 .AND. &
            (Adsorption%Product_DissBond_K(iReactant,iReactNum+nDissoc).EQ.0.) ) THEN
        IF(SpecDSMC(Adsorption%ChemProduct(iReactant,iReactNum+nDissoc))%PolyatomicMol) THEN
          !-------------------------------------------------------------------------------------------------------------------------
          SELECT CASE(MaxDissNum)
          !-------------------------------------------------------------------------------------------------------------------------
          CASE(1) !only one possible dissociation per species
          !-------------------------------------------------------------------------------------------------------------------------
            Adsorption%Product_DissBond_K(iReactant,iReactNum+nDissoc) = &
                    Adsorption%EDissBond(1,Adsorption%ChemProduct(iReactant,iReactNum+nDissoc))
          !-------------------------------------------------------------------------------------------------------------------------
          CASE DEFAULT !more than one dissociation possible per species (special case for some polyatomic)
          !-------------------------------------------------------------------------------------------------------------------------
          WRITE(UNIT=hilf2,FMT='(I2)') iReactNum
            CALL abort(&
__STAMP__&
,'Error in Init_SurfChem: Dissocation bond energy in Disproportionation reaction'//TRIM(hilf2)//' not defined!')
          END SELECT
        ELSE
          Adsorption%Product_DissBond_K(iReactant,iReactNum+nDissoc) = &
                  Adsorption%EDissBond(1,Adsorption%ChemProduct(iReactant,iReactNum+nDissoc))
        END IF
      END IF
    END DO
  END DO
ELSE !MaxDissNum = 0
  nDissoc = 0
  nDisProp = 0
  MaxReactNum = 0
  ALLOCATE(Adsorption%EDissBond(0:1,1:nSpecies))
  ALLOCATE(Adsorption%EDissBondAdsorbPoly(0:1,1:nSpecies))
  Adsorption%EDissBond(:,:)=0.
  Adsorption%EDissBondAdsorbPoly(:,:) = 0.
  DO iSpec = 1,nSpecies
    WRITE(UNIT=hilf,FMT='(I2)') iSpec
    IF (SpecDSMC(iSpec)%InterID.EQ.2) THEN
      Adsorption%EDissBond(0,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Adsorption-EDissBond','0.')
      IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
        Adsorption%EDissBondAdsorbPoly(0,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Adsorption-EDissBondPoly1','0.')
        Adsorption%EDissBondAdsorbPoly(1,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Adsorption-EDissBondPoly2','0.')
        IF (( MAXVAL(Adsorption%DiCoord(:,iSpec)).NE.0) .AND. (Adsorption%EDissBondAdsorbPoly(0,iSpec).EQ.0)) THEN
          CALL abort(&
__STAMP__&
,'Error in Init_SurfChem: Adsorption dissociation bond energy of dicoord for species '//TRIM(hilf)//' not defined!')
        END IF
      END IF
      IF ((.NOT.SpecDSMC(iSpec)%PolyatomicMol).AND.(MAXVAL(Adsorption%DiCoord(:,iSpec)).EQ.0) &
          .AND.(Adsorption%EDissBond(0,iSpec).EQ.0.))THEN
        CALL abort(&
__STAMP__&
,'Error in Init_SurfChem: Adsorption dissociation bond energy for species '//TRIM(hilf)//' not defined!')
      END IF
    END IF
  END DO
END IF !MaxDissNum > 0
! save defined number of surface reactions
Adsorption%DissNum = MaxDissNum
Adsorption%ReactNum = MaxReactNum
Adsorption%nDissocReactions = nDissoc
Adsorption%nDisPropReactions = nDisProp
Adsorption%NumOfExchReact = nDisProp

CalcTST_Case = GETINT('Particles-DSMC-Adsorption-CalcTST','0')
ALLOCATE(Adsorption%TST_Calc(0:Adsorption%ReactNum,1:nSpecies))
Adsorption%TST_Calc(:,:) = .FALSE.
IF (CalcTST_Case.GT.0) CALL Init_TST_Coeff(CalcTST_Case)

SWRITE(UNIT_stdOut,'(A)')' INIT SURFACE CHEMISTRY DONE!'

END SUBROUTINE Init_SurfChem


SUBROUTINE Init_TST_Coeff(TST_Case)
!===================================================================================================================================
! Initializing surface reaction variables
!===================================================================================================================================
! MODULES
USE MOD_Globals,                ONLY : abort, MPIRoot, UNIT_StdOut
USE MOD_Mesh_Vars,              ONLY : nElems
USE MOD_DSMC_Vars,              ONLY : Adsorption, SpecDSMC
USE MOD_PARTICLE_Vars,          ONLY : nSpecies
USE MOD_ReadInTools
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER , INTENT(IN)            :: TST_Case
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: PartitionArraySize, iSpec, iReactNum
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)')' INIT SURFACE TST REACTION COEFFICIENTS!'

Adsorption%PartitionMaxTemp = GETREAL('Particles-DSMC-AdsorptionTST-PartitionMaxTemp','10000.')
Adsorption%PartitionInterval = GETREAL('Particles-DSMC-AdsorptionTST-PartitionInterval','20.')
ALLOCATE(Adsorption%PartitionTemp(1:nElems,1:nSpecies))

!ALLOCATE(Adsorption%TST_Calc(0:Adsorption%ReactNum,1:nSpecies))
!Adsorption%TST_Calc(:,:) = .FALSE.
IF (TST_Case.EQ.1) THEN
  DO iSpec=1,nSpecies
    DO iReactNum=0,Adsorption%ReactNum
      IF ((iReactNum.EQ.0) .AND. (Adsorption%Ads_Prefactor(iSpec).EQ.0.) .AND. (Adsorption%Ads_Powerfactor(iSpec).EQ.0.)) THEN
        Adsorption%TST_Calc(iReactNum,iSpec) = .TRUE.
      ELSE IF ((iReactNum.GT.0) .AND. (iReactNum.LT.Adsorption%DissNum)) THEN
        IF ((Adsorption%Diss_Prefactor(iReactNum,iSpec).EQ.0.) .AND. (Adsorption%Diss_Powerfactor(iReactNum,iSpec).EQ.0.)) THEN
          Adsorption%TST_Calc(iReactNum,iSpec) = .TRUE.
        END IF
      ELSE IF ((iReactNum.GT.0) .AND. (iReactNum.GT.Adsorption%DissNum)) THEN
        !IF ((Adsorption%ER_Prefactor(iReactNum,iSpec).EQ.0.) .AND. (Adsorption%ER_Powerfactor(iReactNum,iSpec).EQ.0.)) THEN
          Adsorption%TST_Calc(iReactNum,iSpec) = .TRUE.
        !END IF
      END IF
    END DO
  END DO
ELSE
  Adsorption%TST_Calc(:,:) = .TRUE. 
END IF

!IF(MOD(Adsorption%PartitionMaxTemp,Adsorption%PartitionInterval).EQ.0.0) THEN
!  PartitionArraySize = INT(Adsorption%PartitionMaxTemp / Adsorption%PartitionInterval)
!ELSE
!  CALL abort(&
!__STAMP__&
!,'ERROR INIT_TST_FACTORS: Partition temperature limit must be multiple of partition interval!')
!END IF  
!
!! calculate array of partitionfunction values
!ALLOCATE(Adsorption%TST_Factors(1:2,0:Adsorption%ReactNum,1:nSpecies)
!DO iSpec=1,nSpecies
!  DO iReactNum=0,Adsorption%ReactNum
!    IF(Adsorption%TST_Calc(iReactNum,iSpec))THEN
!      
!    END IF
!  END DO
!END DO

SWRITE(UNIT_stdOut,'(A)')' INIT SURFACE TST REACTION COEFFICIENTS DONE!'
END SUBROUTINE Init_TST_Coeff


SUBROUTINE FinalizeDSMCSurfModel()
!===================================================================================================================================
! Deallocate Surface model vars
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars,              ONLY : Adsorption, SurfDistInfo
USE MOD_PARTICLE_Vars,          ONLY : PDM, PEM
USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh
#ifdef MPI
USE MOD_Particle_Boundary_Vars, ONLY : SurfCOMM
USE MOD_Particle_MPI_Vars,      ONLY : SurfExchange
USE MOD_Particle_MPI_Vars,      ONLY : AdsorbSendBuf,AdsorbRecvBuf,SurfDistSendBuf,SurfDistRecvBuf
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: subsurfxi,subsurfeta,iSurfSide,iCoord
#ifdef MPI
INTEGER                      :: iProc
#endif /*MPI*/
!===================================================================================================================================
! variables used if particles are kept after adsorption (currentyl not working)
SDEALLOCATE(PDM%ParticleAtWall)
SDEALLOCATE(PDM%PartAdsorbSideIndx)
SDEALLOCATE(PEM%wNumber)
! generaly used adsorption variables
SDEALLOCATE(Adsorption%AdsorpInfo)
SDEALLOCATE(Adsorption%Coverage)
SDEALLOCATE(Adsorption%ProbAds)
SDEALLOCATE(Adsorption%ProbDes)
SDEALLOCATE(Adsorption%SumDesorbPart)
SDEALLOCATE(Adsorption%SumReactPart)
SDEALLOCATE(Adsorption%SumAdsorbPart)
SDEALLOCATE(Adsorption%DensSurfAtoms)
SDEALLOCATE(Adsorption%AreaIncrease)
SDEALLOCATE(Adsorption%CrystalIndx)
SDEALLOCATE(Adsorption%SurfSideToGlobSideMap)
! parameters for Kisliuk and Polanyi Wigner model (wallmodel=1)
SDEALLOCATE(Adsorption%MaxCoverage)
SDEALLOCATE(Adsorption%InitStick)
SDEALLOCATE(Adsorption%PrefactorStick)
SDEALLOCATE(Adsorption%Adsorbexp)
SDEALLOCATE(Adsorption%Nu_a)
SDEALLOCATE(Adsorption%Nu_b)
SDEALLOCATE(Adsorption%DesorbEnergy)
SDEALLOCATE(Adsorption%Intensification)
! parameters for UBI-QEP model (wallmodel=3)
SDEALLOCATE(Adsorption%HeatOfAdsZero)
SDEALLOCATE(Adsorption%DissocReact)
SDEALLOCATE(Adsorption%Diss_Powerfactor)
SDEALLOCATE(Adsorption%Diss_Prefactor)
SDEALLOCATE(Adsorption%EDissBond)
SDEALLOCATE(Adsorption%EDissBondAdsorbPoly)
SDEALLOCATE(Adsorption%AssocReact)
SDEALLOCATE(Adsorption%ChemReactant)
SDEALLOCATE(Adsorption%ChemProduct)
SDEALLOCATE(Adsorption%Reactant_DissBond_K)
SDEALLOCATE(Adsorption%Product_DissBond_K)
SDEALLOCATE(Adsorption%Coordination)
SDEALLOCATE(Adsorption%DiCoord)
SDEALLOCATE(Adsorption%Ads_Powerfactor)
SDEALLOCATE(Adsorption%Ads_Prefactor)
! surfaces distribution variables (currently wallmodel=3)
IF (ALLOCATED(SurfDistInfo)) THEN
DO iSurfSide=1,SurfMesh%nSides
  DO subsurfeta = 1,nSurfSample
  DO subsurfxi = 1,nSurfSample
#ifdef MPI
    SDEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%Nbr_changed)
#endif /*MPI*/
    SDEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%nSites)
    SDEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%SitesRemain)
    SDEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%SurfAtomBondOrder)
    SDEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%desorbnum_tmp)
    SDEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%adsorbnum_tmp)
    SDEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%reactnum_tmp)
    IF (ALLOCATED(SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap)) THEN
      DO iCoord = 1,3
        SDEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(iCoord)%UsedSiteMap)
        SDEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(iCoord)%Species)
        SDEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(iCoord)%BondAtomIndx)
        SDEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(iCoord)%BondAtomIndy)
        SDEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(iCoord)%NeighPos)
        SDEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(iCoord)%NeighSite)
        SDEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap(iCoord)%IsNearestNeigh)
      END DO
      DEALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,iSurfSide)%AdsMap)
    END IF
  END DO
  END DO
END DO
DEALLOCATE(SurfDistInfo)
END IF

#ifdef MPI
IF (ALLOCATED(AdsorbSendBuf)) THEN
  DO iProc=1,SurfCOMM%nMPINeighbors
    SDEALLOCATE(AdsorbSendBuf(iProc)%content_int)
  END DO
  DEALLOCATE(AdsorbSendBuf)
END IF
IF (ALLOCATED(AdsorbRecvBuf)) THEN
  DO iProc=1,SurfCOMM%nMPINeighbors
    SDEALLOCATE(AdsorbRecvBuf(iProc)%content_int)
  END DO
  DEALLOCATE(AdsorbRecvBuf)
END IF
IF (ALLOCATED(SurfCOMM%MPINeighbor)) THEN
  DO iProc=1,SurfCOMM%nMPINeighbors
  SDEALLOCATE(SurfCOMM%MPINeighbor(iProc)%SurfDistSendList)
  SDEALLOCATE(SurfCOMM%MPINeighbor(iProc)%SurfDistRecvList)
  END DO
END IF
SDEALLOCATE(SurfExchange%nSurfDistSidesSend)
SDEALLOCATE(SurfExchange%nSurfDistSidesRecv)
SDEALLOCATE(SurfExchange%SurfDistSendRequest)
SDEALLOCATE(SurfExchange%SurfDistRecvRequest)
IF (ALLOCATED(SurfDistSendBuf)) THEN
  DO iProc=1,SurfCOMM%nMPINeighbors
    SDEALLOCATE(SurfDistSendBuf(iProc)%content_int)
  END DO
  DEALLOCATE(SurfDistSendBuf)
END IF
IF (ALLOCATED(SurfDistRecvBuf)) THEN
  DO iProc=1,SurfCOMM%nMPINeighbors
    SDEALLOCATE(SurfDistRecvBuf(iProc)%content_int)
  END DO
  DEALLOCATE(SurfDistRecvBuf)
END IF
#endif /*MPI*/

END SUBROUTINE FinalizeDSMCSurfModel

END MODULE MOD_DSMC_SurfModelInit
