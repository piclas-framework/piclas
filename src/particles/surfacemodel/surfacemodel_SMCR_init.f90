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

MODULE MOD_SMCR_Init
!===================================================================================================================================
!> Module for initialization of surface models
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
PUBLIC :: InitSMCR
PUBLIC :: InitSMCR_Chem
!===================================================================================================================================
CONTAINS

SUBROUTINE InitSMCR()
!===================================================================================================================================
!> Initializing surface distribution reconstruction model for calculating of coverage effects on heat of adsorption
!> Positions of binding sites in the surface lattice (always rectangular surface lattice assumed)
!> ------------[        surfsquare       ]--------------
!>              |       |       |       |
!>          3---2---3---2---3---2---3---2---3
!>          |       |       |       |       |
!>          2   1   2   1   2   1   2   1   2
!>          |       |       |       |       |
!>          3---2---3---2---3---2---3---2---3
!>          |       |       |       |       |
!>          2   1   2   1   2   1   2   1   2
!>          |       |       |       |       |
!>          3---2---3---2---3---2---3---2---3
!>          |       |       |       |       |
!>          2   1   2   1   2   1   2   1   2
!>          |       |       |       |       |
!>          3---2---3---2---3---2---3---2---3
!>          |       |       |       |       |
!>          2   1   2   1   2   1   2   1   2
!>          |       |       |       |       |
!>          3---2---3---2---3---2---3---2---3
!> For now:
!> Neighbours are all sites, that have the same binding surface atom.
!> Except for top sites(3) they also interact with the next top site.
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_ReadInTools            ,ONLY: GETREAL
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_Particle_Vars          ,ONLY: nSpecies, Species
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption, SurfDistInfo
USE MOD_SurfaceModel_Tools     ,ONLY: UpdateSurfPos
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, PartBound
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
USE MOD_SurfaceModel_MPI       ,ONLY: InitSMCR_MPI
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
CHARACTER(64)                    :: particle_mpf
REAL                             :: surface_mpf
INTEGER                          :: Max_Surfsites_num
INTEGER                          :: Max_Surfsites_own
INTEGER                          :: Max_Surfsites_halo
INTEGER                          :: iSurfSide, iSubSurf, jSubSurf, iSpec
INTEGER                          :: SideID, PartBoundID
INTEGER                          :: surfsquare, dist, Adsorbates
INTEGER                          :: Surfpos, Surfnum, Indx, Indy, UsedSiteMapPos
REAL                             :: RanNum
INTEGER                          :: Coord, nSites, nInterAtom, nNeighbours
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)')' INIT SURFACE DISTRIBUTION...'
ALLOCATE(SurfDistInfo(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides))
DO iSurfSide = 1,SurfMesh%nTotalSides
  DO iSubSurf = 1,nSurfSample
    DO jSubSurf = 1,nSurfSample
      ALLOCATE( SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(1:3),&
                SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%SitesRemain(1:3),&
                SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%adsorbnum_tmp(1:nSpecies),&
                SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%desorbnum_tmp(1:nSpecies),&
                SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%reactnum_tmp(1:nSpecies))
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
  SideID = Adsorption%SurfSideToGlobSideMap(iSurfSide)
  PartboundID = PartBound%MapToPartBC(BC(SideID))
  DO iSubSurf = 1,nSurfSample
    DO jSubSurf = 1,nSurfSample
      IF (PartBound%SolidReactive(PartboundID)) THEN
  !     IF (KeepWallParticles) THEN ! does not work with vMPF
  !       surfsquare = INT(Adsorption%DensSurfAtoms(iSurfSide) &
  !                     * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurfSide) &
  !                     / Species(1)%MacroParticleFactor)
  !       surfsquare = INT(SQRT(REAL(surfsquare))) - 1
  !     END IF
        surfsquare = INT(Adsorption%DensSurfAtoms(iSurfSide) &
                      * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurfSide) &
                      / surface_mpf)
        surfsquare = INT(SQRT(REAL(surfsquare))) - 1
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(1) = INT(surfsquare**2)
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(2) = INT( 2*(surfsquare*(surfsquare+1)) )
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(3) = INT((surfsquare+1)**2)
        IF (surfsquare.LT.2)THEN
          CALL abort(&
            __STAMP__&
            ,'not enough surface spaces for distribution. Surface MacroParticleFactor to to high',surfsquare)
        END IF

        Max_Surfsites_num = Max_Surfsites_num + SUM(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(:))
        IF (iSurfSide.LE.SurfMesh%nSides) THEN
          Max_Surfsites_own = Max_Surfsites_own + SUM(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(:))
        ELSE
          Max_Surfsites_halo = Max_Surfsites_halo + SUM(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(:))
        END IF
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%SitesRemain(:) = SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(:)
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%adsorbnum_tmp(1:nSpecies) = 0.
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%desorbnum_tmp(1:nSpecies) = 0.
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%reactnum_tmp(1:nSpecies) = 0.

        ALLOCATE( SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%SurfAtomBondOrder(1:nSpecies,1:surfsquare+1,1:surfsquare+1))
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%SurfAtomBondOrder(:,:,:) = 0

        ALLOCATE(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1:3))
        DO Coord = 1,3
          SELECT CASE (Coord)
          CASE(1)
            nSites = SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(Coord)
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%nInterAtom = 4
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%nNeighbours = 16
          CASE(2)
            nSites = SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(Coord)
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%nInterAtom = 2
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%nNeighbours = 14
          CASE(3)
            nSites = SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(Coord)
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%nInterAtom = 1
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%nNeighbours = 12
          END SELECT

          nInterAtom = SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%nInterAtom
          nNeighbours = SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%nNeighbours
          ALLOCATE( SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%BondAtomIndx(1:nSites,nInterAtom),&
                    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%BondAtomIndy(1:nSites,nInterAtom),&
                    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%NeighPos(1:nSites,1:nNeighbours),&
                    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%NeighSite(1:nSites,1:nNeighbours),&
                    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%IsNearestNeigh(1:nSites,1:nNeighbours))
          ALLOCATE( SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%UsedSiteMap(1:nSites),&
                    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%Species(1:nSites))
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%UsedSiteMap(:) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%Species(:) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%BondAtomIndx(:,:) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%BondAtomIndy(:,:) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%NeighPos(:,:) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%NeighSite(:,:) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%IsNearestNeigh(:,:) = .FALSE.
        END DO
      ELSE !PartBound%SolidReactive(PartboundID)
        nSites=1 !dummy for correct allocation
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(1:3)=nSites
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%SitesRemain(1:3)=0
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%adsorbnum_tmp(1:nSpecies)=0
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%desorbnum_tmp(1:nSpecies)=0
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%reactnum_tmp(1:nSpecies)=0
        ALLOCATE(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1:3))
        DO Coord = 1,3
          ALLOCATE( SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%UsedSiteMap(1:nSites),&
                    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%Species(1:nSites))
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%UsedSiteMap(:) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%Species(:) = 0
        END DO
      END IF
    END DO
  END DO
END DO

DO iSurfSide = 1,SurfMesh%nTotalSides
SideID = Adsorption%SurfSideToGlobSideMap(iSurfSide)
PartboundID = PartBound%MapToPartBC(BC(SideID))
IF (.NOT.PartBound%SolidReactive(PartboundID)) CYCLE
DO iSubSurf = 1,nSurfSample
DO jSubSurf = 1,nSurfSample
  ! surfsquare chosen from nSite(1) for correct SurfIndx definitions
  surfsquare = INT(SQRT(REAL(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(1))))
  ! allocate and define surface indexes for adsorbate distribution and build mapping of respective bondatoms and neighbours
  Indx = 1
  Indy = 1
  DO Surfpos = 1,SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(1)
    IF (Indx.GT.surfsquare) THEN
      Indx = 1
      Indy = Indy + 1
    END IF
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%UsedSiteMap(Surfpos) = Surfpos
    ! mapping respective neighbours first hollow then bridge then top
    ! hollow
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,1) = Surfpos - surfsquare - 1
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,2) = Surfpos - surfsquare
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,3) = Surfpos - surfsquare + 1
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,4) = Surfpos - 1
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,5) = Surfpos + 1
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,6) = Surfpos + surfsquare - 1
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,7) = Surfpos + surfsquare
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,8) = Surfpos + surfsquare + 1
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,1:8) = 1
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,2) = .TRUE.
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,4) = .TRUE.
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,5) = .TRUE.
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,7) = .TRUE.
    ! bridge
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,9) = Surfpos +(surfsquare+1)*(Indy-1)
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,10) = Surfpos +surfsquare +(surfsquare+1)*(Indy-1)
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,11) = Surfpos +surfsquare +(surfsquare+1)*(Indy-1) +1
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,12) = Surfpos +surfsquare +(surfsquare+1)*(Indy)
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,9:12) = 2
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,9:12) = .TRUE.
    ! top
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,13) = Surfpos + (Indy-1)
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,14) = Surfpos + 1 + (Indy-1)
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,15) = Surfpos + (surfsquare+1) + (Indy-1)
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,16) = Surfpos + (surfsquare+1) + 1 + (Indy-1)
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,13:16) = 3
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,13:16) = .TRUE.
    ! account for empty edges
    IF (Indy .EQ. 1) SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,1:3) = 0
    IF (Indy .EQ. surfsquare) SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,6:8) = 0
    IF (Indx .EQ. 1) THEN
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,1) = 0
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,4) = 0
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,6) = 0
    END IF
    IF (Indx .EQ. surfsquare) THEN
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,3) = 0
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,5) = 0
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,8) = 0
    END IF
    ! mapping respective bond atoms
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%BondAtomIndx(Surfpos,1) = Indx
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%BondAtomIndy(Surfpos,1) = Indy
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%BondAtomIndx(Surfpos,2) = Indx+1
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%BondAtomIndy(Surfpos,2) = Indy
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%BondAtomIndx(Surfpos,3) = Indx
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%BondAtomIndy(Surfpos,3) = Indy+1
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%BondAtomIndx(Surfpos,4) = Indx+1
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%BondAtomIndy(Surfpos,4) = Indy+1
    Indx = Indx + 1
  END DO
  Indx = 1
  Indy = 1
  DO Surfpos = 1,SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(2)
    IF (Indx.GT.(2*surfsquare+1)) THEN
      Indx = 1
      Indy = Indy + 1
    END IF
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%UsedSiteMap(Surfpos) = Surfpos
    IF (Indx .LE. surfsquare) THEN ! surface atoms are LEFT an RIGHT of adsorbate site
      ! mapping respective neighbours first hollow then bridge then top
      ! hollow
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,1) = Surfpos-surfsquare -(surfsquare+1)*(Indy-1) -1
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,2) = Surfpos-surfsquare -(surfsquare+1)*(Indy-1)
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,3) = Surfpos-surfsquare -(surfsquare+1)*(Indy-1) +1
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,4) = Surfpos -(surfsquare+1)*(Indy-1) -1
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,5) = Surfpos -(surfsquare+1)*(Indy-1)
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,6) = Surfpos -(surfsquare+1)*(Indy-1) +1
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,1:6) = 1
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,2) = .TRUE.
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,5) = .TRUE.
      ! bridge
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,7) = Surfpos - (surfsquare+1)
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,8) = Surfpos - (surfsquare+1) + 1
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,9) = Surfpos - 1
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,10) = Surfpos + 1
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,11) = Surfpos + surfsquare
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,12) = Surfpos + surfsquare + 1
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,7:12) = 2
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,8) = .TRUE.
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,11) = .TRUE.
      ! top
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,13) = Surfpos -(surfsquare)*(Indy-1)
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,14) = Surfpos +1 -(surfsquare)*(Indy-1)
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,13:14) = 3
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,13:14) = .TRUE.
      ! account for empty edges
      IF (Indy .EQ. 1) THEN
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,1:3) = 0
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,7:8) = 0
      END IF
      IF (Indy .EQ. surfsquare+1) THEN
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,4:6) = 0
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,11:12) = 0
      END IF
      IF (Indx .EQ. 1) THEN
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,1) = 0
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,4) = 0
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,9) = 0
      END IF
      IF (Indx .EQ. surfsquare) THEN
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,3) = 0
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,6) = 0
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,10) = 0
      END IF
      ! mapping respective bond atoms
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%BondAtomIndx(Surfpos,1) = Indx
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%BondAtomIndy(Surfpos,1) = Indy
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%BondAtomIndx(Surfpos,2) = Indx+1
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%BondAtomIndy(Surfpos,2) = Indy
    ELSE ! surface atoms are TOP and DOWN of adsorbate site
      ! mapping respective neighbours first hollow then bridge then top
      ! hollow
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,1) = Surfpos -(2*surfsquare) &
                                                                                  -(surfsquare+1)*(Indy-1) -1
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,2) = Surfpos -(2*surfsquare)&
                                                                                  -(surfsquare+1)*(Indy-1)
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,3) = Surfpos -surfsquare &
                                                                                  -(surfsquare+1)*(Indy-1) -1
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,4) = Surfpos - surfsquare -(surfsquare+1)*(Indy-1)
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,5) = Surfpos -(surfsquare+1)*(Indy-1) -1
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,6) = Surfpos -(surfsquare+1)*(Indy-1)
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,1:6) = 1
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,3:4) = .TRUE.
      ! bridge
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,7) = Surfpos - surfsquare - (surfsquare+1)
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,8) = Surfpos - surfsquare - 1
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,9) = Surfpos - surfsquare
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,10) = Surfpos + (surfsquare+1) - 1
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,11) = Surfpos + (surfsquare+1)
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,12) = Surfpos + surfsquare + (surfsquare+1)
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,7:12) = 2
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,8:11) = .TRUE.
      ! top
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,13) = Surfpos -surfsquare*(Indy)
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,14) = Surfpos -surfsquare*(Indy) +(surfsquare+1)
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,13:14) = 3
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,13:14) = .TRUE.
      ! account for empty edges
      IF (Indy .EQ. 1) THEN
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,1:2) = 0
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,7) = 0
      END IF
      IF (Indy .EQ. surfsquare) THEN
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,5:6) = 0
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,12) = 0
      END IF
      IF (Indx .EQ. (surfsquare+1)) THEN
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,1) = 0
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,3) = 0
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,5) = 0
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,8) = 0
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,10) = 0
      END IF
      IF (Indx .EQ. 2*surfsquare+1) THEN
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,2) = 0
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,4) = 0
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,6) = 0
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,9) = 0
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,11) = 0
      END IF
      ! mapping respective bond atoms
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%BondAtomIndx(Surfpos,1) = Indx - surfsquare
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%BondAtomIndy(Surfpos,1) = Indy
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%BondAtomIndx(Surfpos,2) = Indx - surfsquare
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%BondAtomIndy(Surfpos,2) = Indy+1
    END IF
    Indx = Indx + 1
  END DO
  Indx = 1
  Indy = 1
  DO Surfpos = 1,SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(3)
    IF (Indx.GT.surfsquare+1) THEN
      Indx = 1
      Indy = Indy + 1
    END IF
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%UsedSiteMap(Surfpos) = Surfpos
    ! mapping respective neighbours first hollow then bridge then top
    ! hollow
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,1) = Surfpos - (surfsquare) - 1 -(Indy-1)
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,2) = Surfpos - (surfsquare) - (Indy-1)
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,3) = Surfpos - 1 - (Indy-1)
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,4) = Surfpos - (Indy-1)
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,1:4) = 1
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%IsNearestNeigh(Surfpos,1:4) = .TRUE.
    ! bridge
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,5) = Surfpos + surfsquare*(Indy-1) -(surfsquare+1)
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,6) = Surfpos -1 +(surfsquare)*(Indy-1)
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,7) = Surfpos +(surfsquare)*(Indy-1)
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,8) = Surfpos + surfsquare*(Indy)
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,5:8) = 2
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%IsNearestNeigh(Surfpos,5:8) = .TRUE.
    ! top
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,9) = Surfpos - (surfsquare+1)
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,10) = Surfpos - 1
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,11) = Surfpos + 1
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,12) = Surfpos + (surfsquare+1)
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,9:12) = 3
    ! account for empty edges
    IF (Indy .EQ. 1) THEN
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,1:2) = 0
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,5) = 0
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,9) = 0
    END IF
    IF (Indx .EQ. 1) THEN
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,1) = 0
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,3) = 0
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,6) = 0
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,10) = 0
    END IF
    IF (Indy .EQ. (surfsquare+1)) THEN
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,3:4) = 0
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,8) = 0
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,12) = 0
    END IF
    IF (Indx .EQ. (surfsquare+1)) THEN
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,2) = 0
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,4) = 0
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,7) = 0
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,11) = 0
    END IF
    ! mapping respective bond atoms
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%BondAtomIndx(Surfpos,1) = Indx
    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%BondAtomIndy(Surfpos,1) = Indy
    Indx = Indx + 1
  END DO

END DO
END DO
END DO

! Use Coverage information to distribute adsorbates randomly on surface
IF (MAXVAL(Adsorption%Coverage(:,:,:,:)).GT.0) THEN
  DO iSurfSide = 1,SurfMesh%nSides
  SideID = Adsorption%SurfSideToGlobSideMap(iSurfSide)
  PartboundID = PartBound%MapToPartBC(BC(SideID))
  IF (.NOT.PartBound%SolidReactive(PartboundID)) CYCLE
  DO iSubSurf = 1,nSurfSample
  DO jSubSurf = 1,nSurfSample
    DO iSpec = 1,nSpecies
      ! adjust coverage to actual discrete value
      Adsorbates = INT(Adsorption%Coverage(iSubSurf,jSubSurf,iSurfSide,iSpec) &
                  * SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(Adsorption%Coordination(PartboundID,iSpec)))
      Adsorption%Coverage(iSubSurf,jSubSurf,iSurfSide,iSpec) = REAL(Adsorbates) &
          / REAL(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(3))
      IF (SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%SitesRemain(Adsorption%Coordination(PartboundID,iSpec)).LT.Adsorbates) THEN
        CALL abort(&
  __STAMP__&
  ,'Error in Init_SurfDist: Too many Adsorbates! - Choose lower Coverages for coordination:', &
  Adsorption%Coordination(PartboundID,iSpec))
      END IF
      ! distribute adsorbates randomly on the surface on the correct site and assign surface atom bond order
      dist = 1
      Coord = Adsorption%Coordination(PartboundID,iSpec)
      Surfnum = SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%SitesRemain(Coord)
      DO WHILE (dist.LE.Adsorbates)
        CALL RANDOM_NUMBER(RanNum)
        Surfpos = 1 + INT(Surfnum * RanNum)
        UsedSiteMapPos = SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%UsedSiteMap(Surfpos)
        ! add species to position and assign bond order of respective surface atoms in the surface lattice
        CALL UpdateSurfPos(iSurfSide,iSubSurf,jSubSurf,Coord,UsedSiteMapPos,iSpec,.FALSE.)
        ! rearrange UsedSiteMap-Surfpos-array
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%UsedSiteMap(Surfpos) = &
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%UsedSiteMap(Surfnum)
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%UsedSiteMap(Surfnum) = UsedSiteMapPos
        Surfnum = Surfnum - 1
        dist = dist + 1
      END DO
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%SitesRemain(Coord) = Surfnum
    END DO
  END DO
  END DO
  END DO
END IF

#if USE_MPI
#ifdef CODE_ANALYZE
! write out the number of sites on all surface of the proc, that are considered for adsorption
WRITE(UNIT_stdOut,'(A,I3,I13,A,I13,A,I13)')' | Maximum number of surface sites on proc: ',myRank,Max_Surfsites_num,&
  ' | own: ',Max_Surfsites_own,' | halo: ',Max_Surfsites_halo
#endif /*CODE_ANALYZE*/
CALL MPI_ALLREDUCE(MPI_IN_PLACE,Max_Surfsites_own,1,MPI_INTEGER,MPI_SUM,PartMPI%COMM,iError) ! write only if mpiroot of all comms
SWRITE(UNIT_stdOut,'(A3,A,I0)') ' > ','Surface sites for all catalytic boundaries: ', Max_SurfSites_own

IF (SurfMesh%SurfOnProc) THEN
  CALL InitSMCR_MPI()
END IF
#endif /*USE_MPI*/

SWRITE(UNIT_stdOut,'(A)')' INIT SURFACE DISTRIBUTION DONE!'

END SUBROUTINE InitSMCR

SUBROUTINE InitSMCR_Chem()
!===================================================================================================================================
!> Initializing surface reaction variables
!===================================================================================================================================
! MODULES
USE MOD_Globals                ,ONLY: abort, MPIRoot, UNIT_StdOut
USE MOD_DSMC_Vars              ,ONLY: SpecDSMC
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_ReadInTools            ,ONLY: GETREAL, GETINT, GETREALARRAY, GETINTARRAY
#if !(USE_LOADBALANCE)
USE MOD_Particle_Boundary_Vars ,ONLY: SurfMesh
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

#if !(USE_LOADBALANCE)
IF (SurfMesh%SurfOnProc .OR. MPIRoot) THEN
#endif
  ! Adsorption constants
  ALLOCATE( Adsorption%Ads_Powerfactor(1:nSpecies),&
            Adsorption%Ads_Prefactor(1:nSpecies))
  DO iSpec = 1,nSpecies
    WRITE(UNIT=hilf,FMT='(I0)') iSpec
    Adsorption%Ads_Powerfactor(iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Adsorption-Powerfactor','0.')
    Adsorption%Ads_Prefactor(iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Adsorption-Prefactor','0.')
  END DO

#if (PP_TimeDiscMethod==42)
  ALLOCATE( Adsorption%AdsorpReactInfo(1:nSpecies))
#endif

  MaxDissNum = GETINT('Surface-MaxDissNum','0')
  MaxAssocNum = MaxDissNum

  ! allocate and initialize dissociative and associative reactions species map
  IF ( (MaxDissNum.GT.0) .OR. (MaxAssocNum.GT.0) ) THEN
    ALLOCATE( Adsorption%DissocReact(1:2,1:MaxDissNum,1:nSpecies),&
              Adsorption%Diss_Powerfactor(1:MaxDissNum,1:nSpecies),&
              Adsorption%Diss_Prefactor(1:MaxDissNum,1:nSpecies))
    ! Read in dissociative reactions and respective dissociation bond energies
    DO iSpec = 1,nSpecies
      WRITE(UNIT=hilf,FMT='(I0)') iSpec
      DO iReactNum = 1,MaxDissNum
        WRITE(UNIT=hilf2,FMT='(I0)') iReactNum
        Adsorption%DissocReact(:,iReactNum,iSpec) = &
                                         GETINTARRAY('Part-Species'//TRIM(hilf)//'-SurfDiss'//TRIM(hilf2)//'-Products',2,'0,0')
        IF ((Adsorption%DissocReact(1,iReactNum,iSpec).GT.nSpecies) &
          .OR.(Adsorption%DissocReact(2,iReactNum,iSpec).GT.nSpecies) ) THEN
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
        END IF
      END DO
      END DO
    END DO
    MaxAssocNum = MAXVAL(nAssocReact)
    Adsorption%NumOfAssocReact = SUM(nAssocReact(:))
    Adsorption%nAssocReactions = SUM(nAssocReact(:))
    Adsorption%RecombNum = MaxAssocNum
    DEALLOCATE(nAssocReact)

    ! fill associative reactions species map from defined dissociative reactions
    MaxReactNum = MaxDissNum + MaxAssocNum
    ALLOCATE( Adsorption%AssocReact(1:2,1:MaxAssocNum,1:nSpecies),&
              Adsorption%ER_Powerfactor(1:MaxAssocNum,1:nSpecies),&
              Adsorption%ER_Prefactor(1:MaxAssocNum,1:nSpecies),&
              Adsorption%EDissBond(0:MaxReactNum,1:nSpecies),&
              Adsorption%EDissBondAdsorbPoly(0:1,1:nSpecies))
    Adsorption%EDissBond(0:MaxReactNum,1:nSpecies) = 0.
    Adsorption%EDissBondAdsorbPoly(0:1,1:nSpecies) = 0.
    DO iSpec = 1,nSpecies
      WRITE(UNIT=hilf,FMT='(I0)') iSpec
      DO iReactNum = 1,MaxDissNum
        WRITE(UNIT=hilf2,FMT='(I0)') iReactNum
        Adsorption%EDissBond(iReactNum,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-SurfDiss'//TRIM(hilf2)//'-EDissBond','0.')
      END DO
    END DO
    DO iSpec = 1,nSpecies
      WRITE(UNIT=hilf,FMT='(I0)') iSpec
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
      WRITE(UNIT=hilf,FMT='(I0)') iSpec
      ReactNum = 1
      DO iSpec2 = 1,nSpecies
      DO iReactNum2 = 1,MaxDissNum
        IF (Adsorption%DissocReact(1,iReactNum2,iSpec2).EQ.iSpec) THEN
          Adsorption%AssocReact(1,ReactNum,iSpec) = Adsorption%DissocReact(2,iReactNum2,iSpec2)
          Adsorption%AssocReact(2,ReactNum,iSpec) = iSpec2
          Adsorption%EDissBond((MaxDissNum+ReactNum),iSpec) = Adsorption%EDissBond(iReactNum2,iSpec2)
          WRITE(UNIT=hilf2,FMT='(I0)') ReactNum
          Adsorption%ER_Powerfactor(ReactNum,iSpec) = &
              GETREAL('Part-Species'//TRIM(hilf)//'-Surf-ER'//TRIM(hilf2)//'-Powerfactor','0.')
          Adsorption%ER_Prefactor(ReactNum,iSpec) = &
              GETREAL('Part-Species'//TRIM(hilf)//'-Surf-ER'//TRIM(hilf2)//'-Prefactor','0.')
          ReactNum = ReactNum + 1
        ELSE IF (Adsorption%DissocReact(2,iReactNum2,iSpec2).EQ.iSpec) THEN
          Adsorption%AssocReact(1,ReactNum,iSpec) = Adsorption%DissocReact(1,iReactNum2,iSpec2)
          Adsorption%AssocReact(2,ReactNum,iSpec) = iSpec2
          Adsorption%EDissBond((MaxDissNum+ReactNum),iSpec) = Adsorption%EDissBond(iReactNum2,iSpec2)
          WRITE(UNIT=hilf2,FMT='(I0)') ReactNum
          Adsorption%ER_Powerfactor(ReactNum,iSpec) = &
              GETREAL('Part-Species'//TRIM(hilf)//'-Surf-ER'//TRIM(hilf2)//'-Powerfactor','0.')
          Adsorption%ER_Prefactor(ReactNum,iSpec) = &
              GETREAL('Part-Species'//TRIM(hilf)//'-Surf-ER'//TRIM(hilf2)//'-Prefactor','0.')
          ReactNum = ReactNum + 1
        ELSE
          CYCLE
        END IF
      END DO
      END DO
      IF (ReactNum.LE.(MaxAssocNum)) THEN
        Adsorption%AssocReact(:,ReactNum:(MaxReactNum-MaxDissNum),iSpec) = 0
        Adsorption%ER_Powerfactor(ReactNum:(MaxReactNum-MaxDissNum),iSpec) = 0.
        Adsorption%ER_Prefactor(ReactNum:(MaxReactNum-MaxDissNum),iSpec) = 0.
      END IF
    END DO

    nDissoc = GETINT('Surface-Nbr-DissocReactions','0')
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
    nDisProp = GETINT('Surface-Nbr-ExchangeReactions','0')
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
  !        DO iReactant=1,2
  !          IF (SpecDSMC(Adsorption%ChemProduct(iReactant,ReactNum))%InterID.EQ.2) THEN
  !            IF(SpecDSMC(Adsorption%ChemProduct(iReactant,ReactNum))%PolyatomicMol) THEN
  !              !------------------------------------------------------------------------------------------------------------------
  !              SELECT CASE(MaxDissNum)
  !              !------------------------------------------------------------------------------------------------------------------
  !              CASE(1) !only one possible dissociation per species
  !              !------------------------------------------------------------------------------------------------------------------
  !                Adsorption%Product_DissBond_K(iReactant,ReactNum) = &
  !                        Adsorption%EDissBond(1,Adsorption%ChemProduct(iReactant,ReactNum))
  !              !------------------------------------------------------------------------------------------------------------------
  !              CASE DEFAULT !more than one dissociation possible per species (special case for some polyatomic)
  !              !------------------------------------------------------------------------------------------------------------------
  !                DO iReactNum2=1,MaxDissNum
  !                  IF (Adsorption%EDissBond(iReactNum2,Adsorption%ChemProduct(iReactant,ReactNum)).LE.0.) CYCLE
  !                  IF ( (Adsorption%Product_DissBond_K(iReactant,ReactNum).GT.&
  !                        Adsorption%EDissBond(iReactNum2,Adsorption%ChemProduct(iReactant,ReactNum))) .AND. &
  !                       (Adsorption%Product_DissBond_K(iReactant,ReactNum).GT.0.) ) THEN
  !                    Adsorption%Product_DissBond_K(iReactant,ReactNum) = &
  !                            Adsorption%EDissBond(iReactNum2,Adsorption%ChemProduct(iReactant,ReactNum))
  !                  END IF
  !                END DO
  !              END SELECT
  !            ELSE
  !              Adsorption%Product_DissBond_K(iReactant,ReactNum) = &
  !                  Adsorption%EDissBond(1,Adsorption%ChemProduct(iReactant,ReactNum))
  !            END IF
  !          END IF
  !        END DO
        END IF
      END DO
    END DO
    ! fill disproportionation reactions (fancy stuff)
    DO iReactNum = 1,nDisProp
      WRITE(UNIT=hilf,FMT='(I0)') iReactNum
      Adsorption%ChemReactant(:,iReactNum+nDissoc) = &
                                        GETINTARRAY('Surface-ExchReact'//TRIM(hilf)//'-Reactants',2,'0,0')
      Adsorption%ChemProduct(:,iReactNum+nDissoc) = &
                                        GETINTARRAY('Surface-ExchReact'//TRIM(hilf)//'-Products',2,'0,0')
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
        WRITE(UNIT=hilf2,FMT='(I0)') Adsorption%ChemReactant(iReactant,iReactNum+nDissoc)
          CALL abort(&
__STAMP__&
,'Error in Init_SurfChem Disproportionation: Dissociation for reactant species '//TRIM(hilf2)//' not defined!')
        END IF
        IF (SpecDSMC(Adsorption%ChemProduct(iReactant,iReactNum+nDissoc))%InterID.EQ.2 .AND. &
              (Adsorption%DissocReact(1,1,Adsorption%ChemProduct(iReactant,iReactNum+nDissoc)).EQ.0) ) THEN
        WRITE(UNIT=hilf2,FMT='(I0)') Adsorption%ChemProduct(iReactant,iReactNum+nDissoc)
          CALL abort(&
__STAMP__&
,'Error in Init_SurfChem Disproportionation: Dissociation for product species '//TRIM(hilf2)//' not defined!')
        END IF
      END DO
      ! Read dissociation bond energies of reactants and products
      Adsorption%Reactant_DissBond_K(:,iReactNum+nDissoc) = &
              GETREALARRAY('Surface-ExchReact'//TRIM(hilf)//'-DissBond_K-Reactants',2,'0.,0.')
      Adsorption%Product_DissBond_K(:,iReactNum+nDissoc) = &
              GETREALARRAY('Surface-ExchReact'//TRIM(hilf)//'-DissBond_K-Products',2,'0.,0.')
      ! Check if dissociation bond energies of reactants and products are all defined if they are at least diatomic
      ! If they are not defined take them from the appropriate dissociation reactions
      DO iReactant=1,2
        IF ( (SpecDSMC(Adsorption%ChemReactant(iReactant,iReactNum+nDissoc))%InterID.EQ.2) .AND. &
              (Adsorption%Reactant_DissBond_K(iReactant,iReactNum+nDissoc).EQ.0.) ) THEN
          IF(SpecDSMC(Adsorption%ChemReactant(iReactant,iReactNum+nDissoc))%PolyatomicMol) THEN
            !-----------------------------------------------------------------------------------------------------------------------
            SELECT CASE(MaxDissNum)
            !-----------------------------------------------------------------------------------------------------------------------
            CASE(1) !only one possible dissociation per species
            !-----------------------------------------------------------------------------------------------------------------------
              Adsorption%Reactant_DissBond_K(iReactant,iReactNum+nDissoc) = &
                      Adsorption%EDissBond(1,Adsorption%ChemReactant(iReactant,iReactNum+nDissoc))
            !-----------------------------------------------------------------------------------------------------------------------
            CASE DEFAULT !more than one dissociation possible per species (special case for some polyatomic)
            !-----------------------------------------------------------------------------------------------------------------------
            WRITE(UNIT=hilf2,FMT='(I0)') iReactNum
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
            !-----------------------------------------------------------------------------------------------------------------------
            SELECT CASE(MaxDissNum)
            !-----------------------------------------------------------------------------------------------------------------------
            CASE(1) !only one possible dissociation per species
            !-----------------------------------------------------------------------------------------------------------------------
              Adsorption%Product_DissBond_K(iReactant,iReactNum+nDissoc) = &
                      Adsorption%EDissBond(1,Adsorption%ChemProduct(iReactant,iReactNum+nDissoc))
            !-----------------------------------------------------------------------------------------------------------------------
            CASE DEFAULT !more than one dissociation possible per species (special case for some polyatomic)
            !-----------------------------------------------------------------------------------------------------------------------
            WRITE(UNIT=hilf2,FMT='(I0)') iReactNum
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
      WRITE(UNIT=hilf,FMT='(I0)') iSpec
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
  Adsorption%RecombNum = MaxAssocNum
  Adsorption%ReactNum = MaxReactNum
  Adsorption%nDissocReactions = nDissoc
  Adsorption%nDisPropReactions = nDisProp
  Adsorption%NumOfExchReact = nDisProp

#if (PP_TimeDiscMethod==42)
  DO iSpec=1,nSpecies
    ALLOCATE( Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(1:Adsorption%ReactNum+1),&
              Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(1:Adsorption%ReactNum+Adsorption%NumOfExchReact+1),&
              Adsorption%AdsorpReactInfo(iSpec)%MeanAdsActE(1:Adsorption%ReactNum),&
              Adsorption%AdsorpReactInfo(iSpec)%MeanSurfActE(1:Adsorption%ReactNum+Adsorption%NumOfExchReact+1),&
              Adsorption%AdsorpReactInfo(iSpec)%ProperSurfActE(1:Adsorption%ReactNum+Adsorption%NumOfExchReact+1),&
              Adsorption%AdsorpReactInfo(iSpec)%AdsReactCount(1:Adsorption%ReactNum+1),&
              Adsorption%AdsorpReactInfo(iSpec)%SurfReactCount(1:Adsorption%ReactNum+Adsorption%NumOfExchReact+1),&
              Adsorption%AdsorpReactInfo(iSpec)%ProperSurfReactCount(1:Adsorption%ReactNum+Adsorption%NumOfExchReact+1))
    Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(:) = 0.
    Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(:) = 0.
    Adsorption%AdsorpReactInfo(iSpec)%MeanAdsActE(:) = 0.
    Adsorption%AdsorpReactInfo(iSpec)%MeanSurfActE(:) = 0.
    Adsorption%AdsorpReactInfo(iSpec)%ProperSurfActE(:) = 0.
    Adsorption%AdsorpReactInfo(iSpec)%AdsReactCount(:) = 0
    Adsorption%AdsorpReactInfo(iSpec)%SurfReactCount(:) = 0
    Adsorption%AdsorpReactInfo(iSpec)%ProperSurfReactCount(:) = 0
  END DO
#endif

  CalcTST_Case = GETINT('Surface-Adsorption-CalcTST','0')
  ALLOCATE(Adsorption%TST_Calc(0:Adsorption%ReactNum,1:nSpecies))
  Adsorption%TST_Calc(:,:) = .FALSE.
  IF (CalcTST_Case.GT.0) CALL Init_TST_Coeff(CalcTST_Case)
#if !(USE_LOADBALANCE)
END IF ! SurfMesh%SurfOnProc .OR. MPIRoot
#endif

SWRITE(UNIT_stdOut,'(A)')' INIT SURFACE CHEMISTRY DONE!'

END SUBROUTINE InitSMCR_Chem


SUBROUTINE Init_TST_Coeff(TST_Case)
!===================================================================================================================================
!> Initializion of the Transition state theory (TST) factors for surface chemistry
!===================================================================================================================================
! MODULES
USE MOD_Globals                ,ONLY: abort,UNIT_StdOut
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_Particle_Boundary_Vars ,ONLY: SurfMesh
USE MOD_ReadInTools            ,ONLY: GETREAL
#if USE_MPI
USE MOD_Globals                ,ONLY: MPIRoot
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER , INTENT(IN)            :: TST_Case
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL                            :: PartitionArraySize
INTEGER                         :: iSpec, iReactNum
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)')' INIT SURFACE TST REACTION COEFFICIENTS!'

Adsorption%PartitionMaxTemp = GETREAL('Surface-AdsorptionTST-PartitionMaxTemp','10000.')
Adsorption%PartitionInterval = GETREAL('Surface-AdsorptionTST-PartitionInterval','20.')
IF (SurfMesh%SurfOnProc) THEN
ALLOCATE(Adsorption%PartitionTemp(1:nElems,1:nSpecies))

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
        IF ((Adsorption%ER_Prefactor(iReactNum,iSpec).EQ.0.) .AND. (Adsorption%ER_Powerfactor(iReactNum,iSpec).EQ.0.)) THEN
          Adsorption%TST_Calc(iReactNum,iSpec) = .TRUE.
        END IF
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
END IF ! SurfMesh%SurfOnProc

SWRITE(UNIT_stdOut,'(A)')' INIT SURFACE TST REACTION COEFFICIENTS DONE!'
END SUBROUTINE Init_TST_Coeff

END MODULE MOD_SMCR_Init
