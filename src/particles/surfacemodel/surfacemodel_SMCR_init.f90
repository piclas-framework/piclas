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
!===================================================================================================================================
CONTAINS

SUBROUTINE InitSMCR()
!===================================================================================================================================
!> Initializing surface distibution reconstruction model for calculating of coverage effects on heat of adsorption
!> For now:
!> Neighbours are all sites, that have the same binding surface atom.
!> Except for top sites(3) they also interact with the next top site.
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_ReadInTools            ,ONLY: GETREAL, GETLOGICAL, GETINT, GETINTARRAY
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_Particle_Vars          ,ONLY: nSpecies, Species
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption, SurfDistInfo, BlockingNeigh
USE MOD_SurfaceModel_Tools     ,ONLY: UpdateSurfPos, SMCR_AdjustMapNum
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, PartBound
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
USE MOD_SurfaceModel_MPI       ,ONLY: InitSMCR_MPI
#endif /*USE_MPI*/
#if (PP_TimeDiscMethod==42)
USE MOD_Particle_Vars          ,ONLY: ManualTimeStep
USE MOD_TimeDisc_Vars          ,ONLY: tend
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
CHARACTER(32)                    :: hilf
CHARACTER(64)                    :: particle_mpf
REAL                             :: surface_mpf
INTEGER                          :: Max_Surfsites_num
INTEGER                          :: Max_Surfsites_own
INTEGER                          :: Max_Surfsites_halo
INTEGER                          :: iSurfSide, iSubSurf, jSubSurf, iSpec
INTEGER                          :: SideID, PartBoundID
INTEGER                          :: surfsquare, Adsorbates
INTEGER                          :: Coord, nSites, nInterAtom, nNeighbours
INTEGER                          :: DistSquareNum
LOGICAL                          :: DistNumCase
INTEGER                          :: BlockedNeightmp(3), i
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
DistNumCase = GETLOGICAL('Particles-Surface-DistNumCase')
IF (DistNumCase) THEN
  DistSquareNum = GETINT('Particles-Surface-DistSquareNumber')
ELSE
  WRITE(UNIT=particle_mpf,FMT='(E11.3)') Species(1)%MacroParticleFactor
  surface_mpf = GETREAL('Particles-Surface-MacroParticleFactor',TRIM(particle_mpf))
#if (PP_TimeDiscMethod==42)
  IF (Adsorption%CoverageReduction) CALL abort(&
__STAMP__&
,'Do not use coverage reduction flag with different surface weightings')
#endif
END IF
Max_Surfsites_num = 0
Max_Surfsites_own = 0
Max_Surfsites_halo = 0

! Allocate and initializes number of surface sites and neighbours
DO iSurfSide = 1,SurfMesh%nTotalSides
  SideID = SurfMesh%SurfIDToSideID(iSurfSide)
  PartboundID = PartBound%MapToPartBC(BC(SideID))
  DO iSubSurf = 1,nSurfSample
    DO jSubSurf = 1,nSurfSample
      IF (PartBound%SurfaceModel(PartboundID).EQ.3) THEN
  !     IF (KeepWallParticles) THEN ! does not work with vMPF
  !       surfsquare = INT(Adsorption%DensSurfAtoms(iSurfSide) &
  !                     * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurfSide) &
  !                     / Species(1)%MacroParticleFactor)
  !       surfsquare = INT(SQRT(REAL(surfsquare))) - 1
  !     END IF
        IF (DistNumCase) THEN
          surfsquare = DistSquareNum - 1
        ELSE
          surfsquare = INT(Adsorption%DensSurfAtoms(iSurfSide) &
                        * SurfMesh%SurfaceArea(iSubSurf,jSubSurf,iSurfSide) &
                        / surface_mpf)
          surfsquare = INT(SQRT(REAL(surfsquare))) - 1
        END IF
        SELECT CASE (PartBound%SolidStructure(PartBoundID))
        CASE(1) !fcc(100)
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(1) = INT(surfsquare**2)
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(2) = INT( 2*(surfsquare*(surfsquare+1)) )
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(3) = INT((surfsquare+1)**2)
        CASE(2) !fcc(111)
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(1) = INT(2*surfsquare**2)
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(2) = INT( surfsquare*(3*(surfsquare+1)-1) )
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(3) = INT((surfsquare+1)**2)
        END SELECT
        IF (surfsquare.LT.9)THEN
          CALL abort(&
            __STAMP__&
            ,'not enough surface spaces for distribution. Surface-MacroParticleFactor too high or DistSquareNumber too low'&
            ,surfsquare)
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
          SELECT CASE (PartBound%SolidStructure(PartBoundID))
          CASE(1) !fcc(100)
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
          CASE(2) !fcc(111)
            SELECT CASE (Coord)
            CASE(1)
              nSites = SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(Coord)
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%nInterAtom = 3
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%nNeighbours = 30
            CASE(2)
              nSites = SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(Coord)
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%nInterAtom = 2
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%nNeighbours = 24
            CASE(3)
              nSites = SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(Coord)
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%nInterAtom = 1
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%nNeighbours = 24
            END SELECT
          END SELECT

          nInterAtom = SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%nInterAtom
          nNeighbours = SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%nNeighbours
          ALLOCATE( SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%BondAtomIndx(1:nSites,nInterAtom),&
                    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%BondAtomIndy(1:nSites,nInterAtom),&
                    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%NeighPos(1:nSites,1:nNeighbours),&
                    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%NeighSite(1:nSites,1:nNeighbours),&
                    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%IsNearestNeigh(1:nSites,1:nNeighbours))
          ALLOCATE( SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%UsedSiteMap(1:nSites),&
                    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%Species(1:nSites),&
                    SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%EVib(1:nSites))
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%UsedSiteMap(:) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%Species(:) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%EVib(:) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%BondAtomIndx(:,:) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%BondAtomIndy(:,:) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%NeighPos(:,:) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%NeighSite(:,:) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%IsNearestNeigh(:,:) = .FALSE.
        END DO
      ELSE !PartBound%Reactive(PartboundID)
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

! which coordinations of nearest neighbour sites block current coordination sites
BlockingNeigh(:,:) = .FALSE.
DO Coord=1,3
  WRITE(UNIT=hilf,FMT='(I0)') Coord
  BlockedNeightmp(1:3) = GETINTARRAY('Surface-Coordination'//TRIM(hilf)//'-BlockingNeigh',3,'0,0,0')
  DO i=1,3
    IF (BlockedNeightmp(i).GT.3 .OR. BlockedNeightmp(i).LT.0 ) THEN
      CALL abort(&
__STAMP__&
,'ERROR: BlockingNeigh has to be 0 =< n =< 3 for Coordination:',Coord)
    END IF
    IF (BlockedNeightmp(i).GT.0) THEN
      BlockingNeigh(Coord,BlockedNeightmp(i)) = .TRUE.
    END IF
  END DO
END DO

CALL Initfcc100Mapping()
CALL Initfcc111Mapping()

#if (PP_TimeDiscMethod==42)
IF(Adsorption%CoverageReduction) ALLOCATE(Adsorption%CovReductionStep(1:nSpecies))
#endif
! Use Coverage information to distribute adsorbates randomly on surface
IF (MAXVAL(Adsorption%Coverage(:,:,:,:)).GT.0) THEN
  DO iSurfSide = 1,SurfMesh%nSides
  SideID = SurfMesh%SurfIDToSideID(iSurfSide)
  PartboundID = PartBound%MapToPartBC(BC(SideID))
  IF (PartBound%SurfaceModel(PartboundID).NE.3) CYCLE
  DO iSubSurf = 1,nSurfSample ;  DO jSubSurf = 1,nSurfSample
    DO iSpec = 1,nSpecies
      ! adjust coverage to discrete integer value
      Adsorbates = INT(Adsorption%Coverage(iSubSurf,jSubSurf,iSurfSide,iSpec)*SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(3))
#if (PP_TimeDiscMethod==42)
      IF(Adsorption%CoverageReduction) THEN
        Adsorption%CovReductionStep(iSpec) = NINT( REAL(Adsorbates) / (tend/ManualTimeStep))
        IF (Adsorption%CovReductionStep(iSpec).LE.0) Adsorption%CovReductionStep(iSpec) = 1
      END IF
#endif
      IF (SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%SitesRemain(Adsorption%Coordination(PartboundID,iSpec)).LT.Adsorbates) THEN
        CALL abort(&
__STAMP__&
,'Error in Init_SurfDist: Too many Adsorbates! - Choose lower Coverages for coordination:', &
Adsorption%Coordination(PartboundID,iSpec))
      END IF
      CALL SMCR_AdjustMapNum(iSubSurf,jSubSurf,iSurfSide,Adsorbates,iSpec)
      Adsorption%Coverage(iSubSurf,jSubSurf,iSurfSide,iSpec) = REAL(Adsorbates) &
          / REAL(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(3))
    END DO
  END DO ; END DO
  END DO
END IF

#ifdef MPI
#ifdef CODE_ANALYZE
! write out the number of sites on all surface of the proc, that are considered for adsorption
WRITE(UNIT_stdOut,'(A,I3,I13,A,I13,A,I13)')' | Maximum number of surface sites on proc: ',myRank,Max_Surfsites_num,&
  ' | own: ',Max_Surfsites_own,' | halo: ',Max_Surfsites_halo
#endif /*CODE_ANALYZE*/
CALL MPI_ALLREDUCE(MPI_IN_PLACE,Max_Surfsites_own,1,MPI_INTEGER,MPI_SUM,PartMPI%COMM,iError) ! write only if mpiroot of all comms
SWRITE(UNIT_stdOut,'(A3,A,I0)') ' > ','Surface sites for all catalytic boundaries: ', Max_SurfSites_own

IF (SurfMesh%SurfOnProc) CALL InitSMCR_MPI()
#endif /*MPI*/

SWRITE(UNIT_stdOut,'(A)')' INIT SURFACE DISTRIBUTION DONE!'

END SUBROUTINE InitSMCR


SUBROUTINE Initfcc100Mapping()
!===================================================================================================================================
!> Initializion of the fcc100 mapping for neighbours and neighrest neighbours site type and position
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
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_SurfaceModel_Vars      ,ONLY: SurfDistInfo
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, PartBound
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: iSurfSide, iSubSurf, jSubSurf
INTEGER                          :: SideID, PartBoundID
INTEGER                          :: surfsquare
INTEGER                          :: Surfpos, Indx, Indy
!===================================================================================================================================
DO iSurfSide = 1,SurfMesh%nTotalSides
SideID = SurfMesh%SurfIDToSideID(iSurfSide)
PartboundID = PartBound%MapToPartBC(BC(SideID))
IF (PartBound%SurfaceModel(PartboundID).NE.3 .OR. PartBound%SolidStructure(PartBoundID).NE.1) CYCLE
DO iSubSurf = 1,nSurfSample
DO jSubSurf = 1,nSurfSample
  ! surfsquare chosen from nSite(1) for correct SurfIndx definitions (Nx-1)
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
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,7:12) = .TRUE.
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,9) = .FALSE.
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,10) = .FALSE.
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

END SUBROUTINE Initfcc100Mapping


SUBROUTINE Initfcc111Mapping()
!===================================================================================================================================
!> Initializion of the fcc111 mapping for neighbours and neighrest neighbours site type and position
!> Positions of binding sites in the surface lattice (always rectangular surface lattice assumed)
!>          ------------[        surfsquare       ]--------------
!>                       |       |       |       |
!>                   3---2---3---2---3---2---3---2---3
!>                  / \  1  / \  1  / \  1  / \  1  /
!>                 2   2   2   2   2   2   2   2   2
!>                /  1  \ /  1  \ /  1  \ /  1  \ /
!>               3---2---3---2---3---2---3---2---3
!>              / \  1  / \  1  / \  1  / \  1  /
!>             2   2   2   2   2   2   2   2   2
!>            /  1  \ /  1  \ /  1  \ /  1  \ /
!>           3---2---3---2---3---2---3---2---3
!>          / \  1  / \  1  / \  1  / \  1  /
!>         2   2   2   2   2   2   2   2   2
!>        /  1  \ /  1  \ /  1  \ /  1  \ /
!>       3---2---3---2---3---2---3---2---3
!>      / \  1  / \  1  / \  1  / \  1  /
!>     2   2   2   2   2   2   2   2   2
!>    /  1  \ /  1  \ /  1  \ /  1  \ /
!>   3---2---3---2---3---2---3---2---3
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_SurfaceModel_Vars      ,ONLY: SurfDistInfo
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, PartBound
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: iSurfSide, iSubSurf, jSubSurf
INTEGER                          :: SideID, PartBoundID
INTEGER                          :: surfsquare
INTEGER                          :: Surfpos, Indx, Indy
!===================================================================================================================================

DO iSurfSide = 1,SurfMesh%nTotalSides
  SideID = SurfMesh%SurfIDToSideID(iSurfSide)
  PartboundID = PartBound%MapToPartBC(BC(SideID))
  IF (PartBound%SurfaceModel(PartboundID).NE.3 .OR. PartBound%SolidStructure(PartBoundID).NE.2) CYCLE
  DO iSubSurf = 1,nSurfSample ; DO jSubSurf = 1,nSurfSample
    ! surfsquare chosen from nSite(3) for correct SurfIndx definitions
    surfsquare = NINT(SQRT(REAL(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(1))/2.))
    ! allocate and define surface indexes for adsorbate distribution and build mapping of respective bondatoms and neighbours
    !sq => surfsquare ,&
    ASSOCIATE ( hstep =>2, &
                b1step=>1, &
                b2step=>2, &
                b3step=>2, &
                tstep =>1 &
               )
    Indx = 1
    Indy = 1
    DO Surfpos = 1,SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(1)
      IF (Indx.GT.(2*surfsquare)) THEN
        Indx = 1
        Indy = Indy + 1
      END IF
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%UsedSiteMap(Surfpos) = Surfpos
      ASSOCIATE ( &! transform sites to other types with corresponding y-index
                  minusyHollow=>2*surfsquare*(indy-2), &
                  plusyHollow=>2*surfsquare*(indy), &
                  mainyHollow=>2*surfsquare*(indy-1), &
                  minusyBridge=>(3*surfsquare+1)*(indy-2), &
                  plusyBridge=>(3*surfsquare+1)*(indy), &
                  mainyBridge=>(3*surfsquare+1)*(indy-1), &
                  minusyTop=>(surfsquare+1)*(indy-2), &
                  plusyTop=>(surfsquare+1)*(indy), &
                  mainyTop=>(surfsquare+1)*(indy-1), &
                  ! transform Hollow sites to other types of same x-index
                  h1toh2=>Indx+1, &
                  h2toh1=>Indx-1, &
                  h1tob1=>(Indx+1)/2, &
                  h1tob2=>Indx+surfsquare, &
                  h1tob3=>Indx+surfsquare+1, &
                  h2tob1=>(Indx)/2, &
                  h2tob2=>Indx+surfsquare-1, &
                  h2tob3=>Indx+surfsquare, &
                  h1tot =>(Indx+1)/2, &
                  h2tot =>Indx/2 &
                  )
        IF (MOD(Indx,2).NE.0) THEN ! first do the triangles with corner at top
          ! mapping respective neighbours first hollow then bridge then top
          ! hollow
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,1)  = minusyHollow +Indx   -hstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,2)  = minusyHollow +h1toh2 -hstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,3)  = minusyHollow +Indx
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,4)  = mainyHollow  +Indx   -hstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,5)  = mainyHollow  +h1toh2 -hstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,6)  = mainyHollow  +h1toh2
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,7)  = mainyHollow  +Indx   +hstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,8)  = plusyHollow  +h1toh2 -hstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,9)  = plusyHollow  +Indx
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,10) = plusyHollow  +h1toh2
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,11) = plusyHollow  +Indx   +hstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,12) = plusyHollow  +h1toh2 +hstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,1:12) = 1
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,5) = .TRUE.
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,6) = .TRUE.
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,10)= .TRUE.
          ! bridge
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,13) = minusyBridge +h1tob3 -b3step
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,14) = minusyBridge +h1tob2
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,15) = mainyBridge  +h1tob1 -b1step
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,16) = mainyBridge  +h1tob1
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,17) = mainyBridge  +h1tob3 -b3step
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,18) = mainyBridge  +h1tob2
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,19) = mainyBridge  +h1tob3
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,20) = mainyBridge  +h1tob2 +b2step
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,21) = plusyBridge  +h1tob1 -b1step
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,22) = plusyBridge  +h1tob1
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,23) = plusyBridge  +h1tob1 +b1step
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,24) = plusyBridge  +h1tob2
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,25) = plusyBridge  +h1tob3
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,26) = plusyBridge  +h1tob2 +b2step
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,27) = plusyBridge  +h1tob3 +b3step
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,13:27) = 2
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,15) = .TRUE.
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,17) = .TRUE.
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,22) = .TRUE.
          ! top
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,28) = mainyTop + h1tot
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,29) = plusyTop + h1tot
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,30) = plusyTop + h1tot +tstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,28:30) = 3
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,28:30) = .TRUE.
          ! account for empty edges
          IF (Indy .EQ. 1) THEN
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,1:3) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,13:14) = 0
          END IF
          IF (Indy .GE. surfsquare) THEN
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,8:12) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,24:27) = 0
          END IF
          IF (Indx .EQ. 1) THEN
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,1:2) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,4:5) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,8) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,13) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,15) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,17) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,21) = 0
          END IF
          IF (Indx .EQ. 2*surfsquare-1) THEN
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,7) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,11:12) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,23) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,27) = 0
          END IF
          ! mapping respective bond atoms
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%BondAtomIndx(Surfpos,1) = h1tot
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%BondAtomIndy(Surfpos,1) = Indy
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%BondAtomIndx(Surfpos,2) = h1tot
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%BondAtomIndy(Surfpos,2) = Indy+1
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%BondAtomIndx(Surfpos,3) = h1tot +tstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%BondAtomIndy(Surfpos,3) = Indy+1
        ELSE ! then do the traingles with corner at bottom
          ! mapping respective neighbours first hollow then bridge then top
          ! hollow
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,1)  = MinusyHollow +h2toh1 -hstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,2)  = MinusyHollow +Indx   -hstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,3)  = MinusyHollow +h2toh1
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,4)  = MinusyHollow +Indx
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,5)  = MinusyHollow +h2toh1 +hstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,6)  = MainyHollow  +Indx   -hstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,7)  = MainyHollow  +h2toh1
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,8)  = MainyHollow  +h2toh1 +hstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,9)  = MainyHollow  +Indx   +hstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,10) = PlusyHollow  +Indx
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,11) = PlusyHollow  +h2toh1 +hstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,12) = PlusyHollow  +Indx   +hstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,1:12) = 1
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,3) = .TRUE.
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,7) = .TRUE.
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,8) = .TRUE.
          ! bridge
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,13) = minusyBridge +h2tob3 -b3step
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,14) = minusyBridge +h2tob2
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,15) = minusyBridge +h2tob3
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,16) = minusyBridge +h2tob2 +b2step
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,17) = mainyBridge  +h2tob1 -b1step
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,18) = mainyBridge  +h2tob1
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,19) = mainyBridge  +h2tob1 +b1step
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,20) = mainyBridge  +h2tob2
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,21) = mainyBridge  +h2tob3
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,22) = mainyBridge  +h2tob2 +b2step
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,23) = mainyBridge  +h2tob3
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,24) = plusyBridge  +h2tob1
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,25) = plusyBridge  +h2tob1 +b1step
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,26) = plusyBridge  +h2tob2 +b2step
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,27) = plusyBridge  +h2tob3 +b3step
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,13:27) = 2
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,18) = .TRUE.
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,21) = .TRUE.
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,22) = .TRUE.
          ! top
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,28) = mainyTop + h2tot
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,29) = mainyTop + h2tot +tstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighPos(Surfpos,30) = plusyTop + h2tot +tstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,28:30) = 3
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%IsNearestNeigh(Surfpos,28:30) = .TRUE.
          ! account for empty edges
          IF (Indy .EQ. 1) THEN
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,1:5) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,13:16) = 0
          END IF
          IF (Indy .GE. surfsquare) THEN
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,10:12) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,26:27) = 0
          END IF
          IF (Indx .EQ. 2) THEN
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,1:2) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,6) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,13) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,17) = 0
          END IF
          IF (Indx .EQ. 2*surfsquare) THEN
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,5) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,8:9) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,11:12) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,19) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,23) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,25) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%NeighSite(Surfpos,27) = 0
          END IF
          ! mapping respective bond atoms
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%BondAtomIndx(Surfpos,1) = h2tot
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%BondAtomIndy(Surfpos,1) = Indy
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%BondAtomIndx(Surfpos,2) = h2tot +tstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%BondAtomIndy(Surfpos,2) = Indy
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%BondAtomIndx(Surfpos,3) = h2tot +tstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(1)%BondAtomIndy(Surfpos,3) = Indy+1
        END IF
      END ASSOCIATE
      Indx = Indx + 1
    END DO
    Indx = 1
    Indy = 1
    DO Surfpos = 1,SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites(2)
      IF (Indx.GT.(2*surfsquare+surfsquare+1)) THEN
        Indx = 1
        Indy = Indy + 1
      END IF
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%UsedSiteMap(Surfpos) = Surfpos
      ASSOCIATE ( &! transform sites to other types with corresponding y-index
                  minusyHollow=>2*surfsquare*(indy-2), &
                  plusyHollow=>2*surfsquare*(indy), &
                  mainyHollow=>2*surfsquare*(indy-1), &
                  minusyBridge=>(3*surfsquare+1)*(indy-2), &
                  plusyBridge=>(3*surfsquare+1)*(indy), &
                  mainyBridge=>(3*surfsquare+1)*(indy-1), &
                  minusyTop=>(surfsquare+1)*(indy-2), &
                  plusyTop=>(surfsquare+1)*(indy), &
                  mainyTop=>(surfsquare+1)*(indy-1), &
                  ! transform Bridge sites to other types of same x-index
                  b1toh1=>Indx*2-1, &
                  b1toh2=>Indx*2, &
                  b1tob2=>Indx*2+surfsquare-1, &
                  b1tob3=>(Indx*2)+surfsquare, &
                  b1tot =>Indx, &
                  b2toh1=>Indx-surfsquare, &
                  b2toh2=>Indx-surfsquare+1, &
                  b2tob1=>(Indx+1-surfsquare)/2, &
                  b2tob3=>Indx+1, &
                  b2tot =>(Indx+1-surfsquare)/2, &
                  b3toh1=>Indx-surfsquare-1, &
                  b3toh2=>Indx-surfsquare, &
                  b3tob1=>(Indx-surfsquare)/2, &
                  b3tob2=>Indx-1, &
                  b3tot =>(indx-surfsquare)/2 &
                 )
        IF (Indx .LE. surfsquare) THEN ! surface atoms are LEFT an RIGHT of adsorbate site
          ! mapping respective neighbours first hollow then bridge then top
          ! hollow
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,1)  = MinusyHollow +b1toh1 -hstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,2)  = MinusyHollow +b1toh2 -hstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,3)  = MinusyHollow +b1toh1
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,4)  = MinusyHollow +b1toh2
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,5)  = MinusyHollow +b1toh1 +hstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,6)  = MainyHollow  +b1toh2 -hstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,7)  = MainyHollow  +b1toh1
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,8)  = MainyHollow  +b1toh2
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,9)  = MainyHollow  +b1toh1 +hstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,10) = MainyHollow  +b1toh2 +hstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,1:10) = 1
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,3) = .TRUE.
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,8) = .TRUE.
          ! bridge
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,11) = minusyBridge +b1tob3 -b3step
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,12) = minusyBridge +b1tob2
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,13) = minusyBridge +b1tob3
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,14) = minusyBridge +b1tob2 +b2step
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,15) = mainyBridge  +Indx   -b1step
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,16) = mainyBridge  +Indx   +b1step
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,17) = mainyBridge  +b1tob2
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,18) = mainyBridge  +b1tob3
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,19) = mainyBridge  +b1tob2 +b2step
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,20) = mainyBridge  +b1tob3 +b3step
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,11:20) = 2
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,12:13) = .TRUE.
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,18:19) = .TRUE.
          ! top
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,21) = minusyTop +b1tot
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,22) = mainyTop  +b1tot
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,23) = mainyTop  +b1tot +tstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,24) = plusyTop  +b1tot +tstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,21:24) = 3
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,22:23) = .TRUE.
          ! account for empty edges
          IF (Indy .EQ. 1) THEN
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,1:5) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,11:14) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,21) = 0
          END IF
          IF (Indy .EQ. surfsquare+1) THEN
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,6:10) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,17:20) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,24) = 0
          END IF
          IF (Indx .EQ. 1) THEN
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,1:2) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,6) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,11) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,15) = 0
          END IF
          IF (Indx .EQ. surfsquare) THEN
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,5) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,9:10) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,16) = 0
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,20) = 0
          END IF
          ! mapping respective bond atoms
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%BondAtomIndx(Surfpos,1) = b1tot
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%BondAtomIndy(Surfpos,1) = Indy
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%BondAtomIndx(Surfpos,2) = b1tot +tstep
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%BondAtomIndy(Surfpos,2) = Indy
        ELSE ! surface atoms are TOP and DOWN of adsorbate site
          IF (MOD((Indx-surfsquare),2).NE.0) THEN ! first do / direction
            ! mapping respective neighbours first hollow then bridge then top
            ! hollow
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,1)  = MinusyHollow +b2toh1 -hstep
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,2)  = MinusyHollow +b2toh2 -hstep
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,3)  = MinusyHollow +b2toh1
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,4)  = MainyHollow  +b2toh1 -hstep
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,5)  = MainyHollow  +b2toh1 -hstep
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,6)  = MainyHollow  +b2toh1
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,7)  = MainyHollow  +b2toh2
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,8)  = PlusyHollow  +b2toh2 -hstep
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,9)  = PlusyHollow  +b2toh1
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,10) = PlusyHollow  +b2toh2
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,1:10) = 1
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,5:6) = .TRUE.
            ! bridge
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,11) = minusyBridge +b2tob3 -b3step
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,12) = minusyBridge +Indx
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,13) = mainyBridge  +b2tob1 -b1step
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,14) = mainyBridge  +b2tob1
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,15) = mainyBridge  +b2tob3 -b3step
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,16) = mainyBridge  +b2tob3
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,17) = plusyBridge  +b2tob1 -b1step
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,18) = plusyBridge  +b2tob1
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,19) = plusyBridge  +Indx
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,20) = plusyBridge  +b2tob3
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,11:20) = 2
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,13) = .TRUE.
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,15:16) = .TRUE.
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,18) = .TRUE.
            ! top
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,21) = mainyTop + b2tot -tstep
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,22) = mainyTop + b2tot
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,23) = plusyTop + b2tot
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,24) = plusyTop + b2tot +tstep
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,21:24) = 3
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,22:23) = .TRUE.
            ! account for empty edges
            IF (Indy .EQ. 1) THEN
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,1:3) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,11:12) = 0
            END IF
            IF (Indy .EQ. surfsquare) THEN
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,8:10) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,19:20) = 0
            END IF
            IF ((Indx-surfsquare) .EQ. 1) THEN
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,1:2) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,4:5) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,8) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,11) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,13) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,15) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,17) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,21) = 0
            END IF
            IF ((Indx-surfsquare) .EQ. 2*surfsquare+1) THEN
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,3) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,6:7) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,9:10) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,14) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,16) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,18) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,20) = 0
            END IF
            ! mapping respective bond atoms
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%BondAtomIndx(Surfpos,1) = b2tot
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%BondAtomIndy(Surfpos,1) = Indy
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%BondAtomIndx(Surfpos,2) = b2tot
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%BondAtomIndy(Surfpos,2) = Indy+1
          ELSE ! then do \ direction
            ! mapping respective neighbours first hollow then bridge then top
            ! hollow
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,1)  = MinusyHollow +b3toh1 -hstep
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,2)  = MinusyHollow +b3toh2 -hstep
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,3)  = MinusyHollow +b3toh1
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,4)  = MainyHollow  +b3toh2 -hstep
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,5)  = MainyHollow  +b3toh1
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,6)  = MainyHollow  +b3toh2
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,7)  = MainyHollow  +b3toh1 +hstep
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,8)  = PlusyHollow  +b3toh2
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,9)  = PlusyHollow  +b3toh1 +hstep
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,10) = PlusyHollow  +b3toh2 +hstep
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,1:10) = 1
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,5:6) = .TRUE.
            ! bridge
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,11) = minusyBridge +Indx   -b3step
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,12) = minusyBridge +b3tob2
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,13) = mainyBridge  +b3tob1 -b1step
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,14) = mainyBridge  +b3tob1
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,15) = mainyBridge  +b3tob2
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,16) = mainyBridge  +b3tob2 +b2step
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,17) = plusyBridge  +b3tob1
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,18) = plusyBridge  +b3tob1 +b1step
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,19) = plusyBridge  +b3tob2 +b2step
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,20) = plusyBridge  +Indx   +b2step
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,11:20) = 2
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,14:17) = .TRUE.
            ! top
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,21) = mainyTop +b3tot
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,22) = mainyTop +b3tot +tstep
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,23) = plusyTop +b3tot
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighPos(Surfpos,24) = plusyTop +b3tot +tstep
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,21:24) = 3
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%IsNearestNeigh(Surfpos,22:23) = .TRUE.
            ! account for empty edges
            IF (Indy .EQ. 1) THEN
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,1:3) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,4) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,11:12) = 0
            END IF
            IF (Indy .EQ. surfsquare) THEN
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,8:10) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,19:20) = 0
            END IF
            IF ((Indx-surfsquare) .EQ. 2) THEN
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,1:2) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,4) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,11) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,13) = 0
            END IF
            IF ((Indx-surfsquare) .EQ. 2*surfsquare) THEN
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,7) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,9:10) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,18) = 0
              SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%NeighSite(Surfpos,20) = 0
            END IF
            ! mapping respective bond atoms
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%BondAtomIndx(Surfpos,1) = b3tot
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%BondAtomIndy(Surfpos,1) = Indy
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%BondAtomIndx(Surfpos,2) = b3tot +tstep
            SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(2)%BondAtomIndy(Surfpos,2) = Indy+1
          END IF
        END IF
      END ASSOCIATE
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
      ASSOCIATE ( &! transform sites to other types with corresponding y-index
                  minusyHollow=>2*surfsquare*(indy-2), &
                  plusyHollow=>2*surfsquare*(indy), &
                  mainyHollow=>2*surfsquare*(indy-1), &
                  minusyBridge=>(3*surfsquare+1)*(indy-2), &
                  plusyBridge=>(3*surfsquare+1)*(indy), &
                  mainyBridge=>(3*surfsquare+1)*(indy-1), &
                  minusyTop=>(surfsquare+1)*(indy-2), &
                  plusyTop=>(surfsquare+1)*(indy), &
                  mainyTop=>(surfsquare+1)*(indy-1), &
                  ! transform Top sites to other types of same x-index
                  ttoh1=>Indx*2-1, &
                  ttoh2=>Indx*2, &
                  ttob1=>Indx, &
                  ttob2=>Indx*2-1+surfsquare, &
                  ttob3=>Indx*2+surfsquare &
                 )
        ! mapping respective neighbours first hollow then bridge then top
        ! hollow
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,1) = minusyHollow +ttoh1 -hstep
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,2) = minusyHollow +ttoh2 -hstep
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,3) = minusyHollow +ttoh1
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,4) = mainyHollow  +ttoh2 -hstep
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,5) = mainyHollow  +ttoh1
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,6) = mainyHollow  +ttoh2
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,1:6) = 1
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%IsNearestNeigh(Surfpos,1:6) = .TRUE.
        ! bridge
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,7)  = minusyBridge +ttob1 -b1step
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,8)  = minusyBridge +ttob2 -b2step
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,9)  = minusyBridge +ttob3 -b3step
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,10) = minusyBridge +ttob2
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,11) = minusyBridge +ttob3
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,12) = mainyBridge +ttob1 -b1step
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,13) = mainyBridge  +ttob1
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,14) = mainyBridge  +ttob3 -b3step
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,15) = mainyBridge  +ttob2
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,16) = mainyBridge  +ttob3
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,17) = mainyBridge  +ttob2 +b2step
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,18) = plusyBridge  +ttob1
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,7:18) = 2
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%IsNearestNeigh(Surfpos,9:10) = .TRUE.
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%IsNearestNeigh(Surfpos,12:13) = .TRUE.
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%IsNearestNeigh(Surfpos,15:16) = .TRUE.
        ! top
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,19) = minusyTop +Indx -1
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,20) = minusyTop +Indx
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,21) = mainyTop  +Indx -1
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,22) = mainyTop  +Indx +1
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,23) = plusyTop  +Indx
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighPos(Surfpos,24) = plusyTop  +Indx +1
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,19:24) = 3
        ! account for empty edges
        IF (Indy .EQ. 1) THEN
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,1:3) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,7:11) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,19:20) = 0
        END IF
        IF (Indx .EQ. 1) THEN
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,1:2) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,4) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,7:9) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,12) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,14) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,19) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,21) = 0
        END IF
        IF (Indy .EQ. (surfsquare+1)) THEN
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,4:6) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,14:18) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,23:24) = 0
        END IF
        IF (Indx .EQ. (surfsquare+1)) THEN
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,3) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,5:6) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,11) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,13) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,16:18) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,22) = 0
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%NeighSite(Surfpos,24) = 0
        END IF
        ! mapping respective bond atoms
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%BondAtomIndx(Surfpos,1) = Indx
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(3)%BondAtomIndy(Surfpos,1) = Indy
      END ASSOCIATE
      Indx = Indx + 1
    END DO
    END ASSOCIATE
  END DO; END DO
END DO

!print*,'site type 1 ========================================'
!DO iSubSurf=1,SurfDistInfo(1,1,1)%nSites(1)
!  DO iSurfSide=1,30
!    print*,'SurfPos: ',iSubSurf,'        Neighbour: ',iSurfSide,'      neighpos: '&
!      ,SurfDistInfo(1,1,1)%AdsMap(1)%NeighPos(iSubSurf,iSurfSide),&
!      '           neighsite: ',SurfDistInfo(1,1,1)%AdsMap(1)%NeighSite(iSubSurf,iSurfSide)
!  END DO
!  print*,'----------------------------------------------------'
!END DO
!print*,'site type 2 ========================================'
!DO iSubSurf=1,SurfDistInfo(1,1,1)%nSites(2)
!  DO iSurfSide=1,24
!    print*,'SurfPos: ',iSubSurf,'        Neighbour: ',iSurfSide,'      neighpos: '&
!      ,SurfDistInfo(1,1,1)%AdsMap(2)%NeighPos(iSubSurf,iSurfSide),&
!      '           neighsite: ',SurfDistInfo(1,1,1)%AdsMap(2)%NeighSite(iSubSurf,iSurfSide)
!  END DO
!  print*,'----------------------------------------------------'
!END DO
!print*,'site type 3 ========================================'
!DO iSubSurf=1,SurfDistInfo(1,1,1)%nSites(3)
!  DO iSurfSide=1,24
!    print*,'SurfPos: ',iSubSurf,'        Neighbour: ',iSurfSide,'      neighpos: '&
!      ,SurfDistInfo(1,1,1)%AdsMap(3)%NeighPos(iSubSurf,iSurfSide),&
!      '           neighsite: ',SurfDistInfo(1,1,1)%AdsMap(3)%NeighSite(iSubSurf,iSurfSide)
!  END DO
!  print*,'----------------------------------------------------'
!END DO

END SUBROUTINE Initfcc111Mapping


END MODULE MOD_SMCR_Init
