!==================================================================================================================================
! Copyright (c) 2015-2019 Wladimir Reschke
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

MODULE MOD_SurfaceModel_Init
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
PUBLIC :: DefineParametersSurfModel
PUBLIC :: InitSurfaceModel
PUBLIC :: FinalizeSurfaceModel
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for surface model
!==================================================================================================================================
SUBROUTINE DefineParametersSurfModel()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("SurfaceModel")

!CALL prms%CreateLogicalOption(  'Particles-KeepWallParticles'&
!  , 'Flag to track adsorbed Particles on surface if they are adsorbed. Currently only [FALSE] implemented','.FALSE.')
CALL prms%CreateLogicalOption(  'Surface-ModelERSpecular'&
  , 'Flag for specular reflection for ER-reaction particles [FALSE] diffuse reflection','.FALSE.')
#if (PP_TimeDiscMethod==42)
CALL prms%CreateLogicalOption(  'Surface-Adsorption-LateralInactive'&
  , 'Flag for disabling lateral innteractions. Only for TD=42 (RESERVOIR)','.FALSE.')
CALL prms%CreateLogicalOption(  'Surface-Adsorption-CoverageReduction'&
  , 'Flag for reducing coverage by time dependant interval. Interval=coverage_start/num_iter. Only for TD=42 (RESERVOIR)','.FALSE.')
CALL prms%CreateLogicalOption(  'Surface-Adsorption-doTPD'&
  , 'Flag for TPD spectrum calculation. Only for TD=42 (RESERVOIR)','.FALSE.')
CALL prms%CreateRealOption(     'Surface-Adsorption-TPD-Beta'&
  , 'Temperature slope used for TPD calculation [K/s]. Surface temperature used for desorption probability is increased'//&
    'each iteration','0.')
#endif
CALL prms%CreateStringOption(  'Particles-SurfCoverageFile'&
  , 'Give relative path to Surface-state-file ([Projectname]_DSMCSurfState***.h5) used for coverage init. '//&
    'File must be of same format as calculation (same mesh, same amount of species, cells and surfacemodel).\n'//&
    'If no file specified:\n'// &
    'Define Part-Species[$]-PartBound[$]-InitialCoverage for initial coverage different than 0.' &
  , 'none')

CALL prms%CreateRealOption(     'Part-Species[$]-PartBound[$]-InitialCoverage'&
  , 'Initial coverage used for species [$] and surfaces of boundary [$] in case of no surface-state-file init' &
  , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-PartBound[$]-ReactionAccomodation'&
  , 'Define energy accomodation coefficient.\n'//&
    'Describes the percentage of reaction enthalpy of surface reaction transmitted to surface.','1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Adsorption-EDissBond', 'Bond dissociation energy (K) needed for calculation of'//&
                                           'adsorption enthalpy in bond order methodology' &
                                         , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Adsorption-EDissBondPoly1'&
                                         , 'Dissociation energy (K) of first bond needed for calculation of'//&
                                           'adsorption enthalpy in bond order methodology for Poly atomic species' &
                                         , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Adsorption-EDissBondPoly2'&
                                         , 'Dissociation energy (K) of second bond needed for calculation of'//&
                                           'adsorption enthalpy in bond order methodology for Poly atomic species' &
                                         , '0.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Surface-Adsorption-CalcTST'&
                                          , 'Define, which rate coefficient for adsorption are used:\n'//&
                                          ' 0: only those, specified in INI\n'//&
                                          ' 1: those not define in INI are calculated with TST\n'//&
                                          ' 2: calculate every coefficient with TST','2')
CALL prms%CreateIntOption(      'Surface-MaxDissNum', 'Define maximum number of Dissociations per species','0')
CALL prms%CreateIntArrayOption( 'Part-Species[$]-SurfDiss[$]-Products'&
                                          , 'Define product species for surface dissociation number of considered species' &
                                          , '0 , 0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-SurfDiss[$]-EDissBond'&
                                         , 'Bond dissociation energy (K) for dissociation into product species of'//&
                                          ' surface dissociation','0.', numberedmulti=.TRUE.)

CALL prms%CreateRealOption(     'Part-Species[$]-Adsorption-Powerfactor'&
                                          , 'TODO-DEFINE-PARAMETER','-1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Adsorption-Prefactor'&
                                          , 'TODO-DEFINE-PARAMETER','-1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-SurfDiss[$]-Powerfactor'&
                                         , 'TODO-DEFINE-PARAMETER','-1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-SurfDiss[$]-Prefactor'&
                                         , 'TODO-DEFINE-PARAMETER','-1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surf-ER[$]-Powerfactor'&
                                          , 'TODO-DEFINE-PARAMETER','-1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surf-ER[$]-Prefactor'&
                                          , 'TODO-DEFINE-PARAMETER','-1.', numberedmulti=.TRUE.)

CALL prms%CreateIntOption(      'Surface-Nbr-DissocReactions'&
                                          , 'TODO-DEFINE-PARAMETER\n'//&
                                          'Resulting species for given dissoc (2,MaxDissNum,nSpecies)','0')
CALL prms%CreateIntOption(      'Surface-Nbr-ExchangeReactions'&
                                          , 'TODO-DEFINE-PARAMETER','0')
CALL prms%CreateIntArrayOption( 'Surface-ExchReact[$]-Reactants'&
                                          , 'TODO-DEFINE-PARAMETER','0 , 0', numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'Surface-ExchReact[$]-Products'&
                                          , 'TODO-DEFINE-PARAMETER','0 , 0', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Surface-ExchReact[$]-DissBond_K-Reactants'&
                                          , 'TODO-DEFINE-PARAMETER','0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Surface-ExchReact[$]-DissBond_K-Products'&
                                         , 'TODO-DEFINE-PARAMETER','0. , 0.', numberedmulti=.TRUE.)

CALL prms%SetSection("Surfacemodel1 [currently not working]")

CALL prms%CreateRealOption(     'Part-Species[$]-MaximumCoverage'&
  , 'Maximum coverage of surfaces with species [$] (used for surfacemodel=1)','0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-InitialStick'&
  , 'Initial sticking coefficient (S_0) of species [$] for surfaces (used for Kisliuk model, surfacemodel=1)' &
  , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-PrefactorStick'&
  , 'Prefactor of sticking coefficient of species [$] for surfaces (used for Kisliuk model, surfacemodel=1)' &
  , '0.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-Adsorbexp'&
  , 'Adsorption exponent of species [$] for surfaces (used for Kisliuk model, surfacemodel=1)' &
  , '1', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Nu-a'&
                                         , 'TODO-DEFINE-PARAMETER\n'//&
                                              'Nu exponent a for surface n','0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Nu-b'&
                                         , 'TODO-DEFINE-PARAMETER\n'//&
                                              'Nu exponent b for surface n','0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Desorption-Energy-K'&
                                          , 'TODO-DEFINE-PARAMETER\n'//&
                                            'Desorption energy (K) for surface n','1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Intensification-K'&
                                          , 'TODO-DEFINE-PARAMETER\n'//&
                                             'Intensification energy (K) for surface n','0.', numberedmulti=.TRUE.)


CALL prms%SetSection("Surfacemodel2")

CALL prms%CreateIntOption(     'Part-Species[$]-PartBound[$]-ResultSpec'&
                               ,'Resulting recombination species (one of nSpecies)','-1', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-PartBound[$]-RecombinationCoeff'&
                                          , 'TODO-DEFINE-PARAMETER','0.', numberedmulti=.TRUE.)


CALL prms%SetSection("SMCR")

CALL prms%CreateLogicalOption(     'Surface-Adsorption-EnableAttraction'&
  , 'Activates attracive forces between adsorbates. \n[TRUE] or [FALSE].','.FALSE.')
CALL prms%CreateLogicalOption(     'Surface-Adsorption-NoDiffusion'&
  , 'Disables diffusion into equlibrium on the surface. \n[TRUE] or [FALSE].','.FALSE.')
CALL prms%CreateLogicalOption(     'Particles-Surface-DistNumCase'&
  , 'Sets the surface distribution case. \nSpecific number[TRUE] or weighting[FALSE].','.TRUE.')
CALL prms%CreateIntOption(     'Particles-Surface-DistSquareNumber'&
  , 'Defines the number of simulated surface atoms for every surface element (DistSquareNumber x DistSquareNumber) '//&
    '[surfacemodel=3]. \nIf one surface contains less then 10x10 surface atoms program abort is called.\n', '20')
CALL prms%CreateRealOption(     'Particles-Surface-MacroParticleFactor'&
  , 'Weighting factor used for particles adsorbed on surface in case of reconstruction [surfacemodel=3].\n'//&
    'If one surface contains less then 5x5 surface atoms program abort is called.\n'//&
    'Default: Species(1)%MPF: Uses macro particle factor of species1.')
CALL prms%CreateIntArrayOption( 'Surface-Coordination[$]-BlockingNeigh'&
                              , 'Define which neighbour coordinations can block considered adsorption position' &
                              , '0 , 0 , 0', numberedmulti=.TRUE.)

CALL prms%CreateIntOption(      'Part-Species[$]-PartBound[$]-Coordination'&
  , 'Coordination at which particle of species [$] is bound on surface of boundary [$].\n'//&
    '1=hollow\n'//&
    '2=bridge\n'//&
    '3=on-top'//&
    '[surfacemodel=3]','0', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-PartBound[$]-DiCoordination'&
  , 'Define energy interaction type for particle of species [$] at Boundary [$] (di-, polyatomic).\n'//&
    '1: strong, erect\n'//&
    '2: weak, erect\n'//&
    '3: intermediate (strong+weak)/2 \n'//&
    '4: parallel, bridge span, acceptor \n'//&
    '5: on top, parallel to one surface atom \n'//&
    '6: parallel, across bridge, donor \n'//&
    '7: chelating, similar to 4 for poly \n'//&
    '[surfacemodel=3]','1', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-PartBound[$]-HeatOfAdsorption-K'&
  , 'Define heat of adsorption [K] on clear surface for binding atom of species [$] on boundary [$].\n'// &
    '[Assumption of on-top side bind, surfacemodel=3]','0.', numberedmulti=.TRUE.)

CALL prms%SetSection("Liquid-Surface")

CALL prms%CreateLogicalOption(     'Part-Species[$]-PartBound[$]-LiquidSpec'&
  , 'Sets the flag for being treated as surface species. If marked, antoine parameters must to be defined'&
  ,'.FALSE.', numberedmulti=.TRUE.)

CALL prms%CreateRealArrayOption(  'Part-Species[$]-ParamAntoine'  &
                                , 'Parameters for Antoine Eq (vapor pressure)', '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-condensCase'&
                                          , 'TODO-DEFINE-PARAMETER','1', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-liquidAlpha'&
                                          , 'TODO-DEFINE-PARAMETER','1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-liquidBeta'&
                                          , 'TODO-DEFINE-PARAMETER','0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption(  'Part-Species[$]-liquidBetaCoeff'  &
                                          , 'TODO-DEFINE-PARAMETER','0. , 0. , 0. , 0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-liquidTkrit'&
                                          , 'TODO-DEFINE-PARAMETER','100', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-liquidTmelt'&
                                          , 'TODO-DEFINE-PARAMETER','50', numberedmulti=.TRUE.)

END SUBROUTINE DefineParametersSurfModel


SUBROUTINE InitSurfaceModel()
!===================================================================================================================================
!> Initialize surface model variables
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars               ,ONLY: BoltzmannConst, PI
USE MOD_Globals
USE MOD_Mesh_Vars                  ,ONLY: nElems, BC
USE MOD_DSMC_Vars                  ,ONLY: DSMC, CollisMode, SpecDSMC
USE MOD_Particle_Vars              ,ONLY: nSpecies, PDM, WriteMacroSurfaceValues
USE MOD_Particle_Vars              ,ONLY: KeepWallParticles, PEM, Species
USE MOD_ReadInTools                ,ONLY: GETINT,GETREAL,GETLOGICAL,GETREALARRAY
USE MOD_Particle_Boundary_Vars     ,ONLY: nSurfSample, SurfMesh, nPartBound, PartBound
USE MOD_Particle_Boundary_Sampling ,ONLY: InitParticleBoundarySampling
USE MOD_SurfaceModel_Vars          ,ONLY: Adsorption, ModelERSpecular, SurfModel, SpecSurf
USE MOD_SurfaceModel_Tools         ,ONLY: CalcAdsorbProb, CalcDesorbProb
USE MOD_SMCR_Init                  ,ONLY: InitSMCR
#if USE_MPI
USE MOD_SurfaceModel_MPI           ,ONLY: InitSurfModel_MPI
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
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
!===================================================================================================================================
IF (.NOT.(ANY(PartBound%Reactive))) RETURN
IF (CollisMode.LE.1) THEN
  CALL abort(&
__STAMP__&
,'Error in InitSurfaceModel - wrong collismode! needs to be >1')
END IF
! initialize variables only for processors that have any surfaces in own domain else they are skipped or not allocated
!KeepWallParticles = GETLOGICAL('Particles-KeepWallParticles','.FALSE.')
ModelERSpecular = GETLOGICAL('Surface-ModelERSpecular')
Adsorption%EnableAdsAttraction = GETLOGICAL('Surface-Adsorption-EnableAttraction')
Adsorption%NoDiffusion = GETLOGICAL('Surface-Adsorption-NoDiffusion')
KeepWallParticles = .FALSE.
IF (KeepWallParticles) THEN
  IF(SurfMesh%SurfOnProc) THEN
    ALLOCATE(PDM%ParticleAtWall(1:PDM%maxParticleNumber)  , &
            PDM%PartAdsorbSideIndx(1:3,1:PDM%maxParticleNumber))
    PDM%ParticleAtWall(1:PDM%maxParticleNumber) = .FALSE.
  END IF !SurfMesh%SurfOnProc
  ALLOCATE(PEM%wNumber(1:nElems))
END IF
! allocate info and constants
ALLOCATE( SurfModel%Info(1:nSpecies))
DO iSpec = 1,nSpecies
  SurfModel%Info(iSpec)%WallCollCount = 0
  SurfModel%Info(iSpec)%NumOfAds = 0
  SurfModel%Info(iSpec)%NumOfDes = 0
  SurfModel%Info(iSpec)%MeanProbAds  = 0.
  SurfModel%Info(iSpec)%MeanProbAdsCount  = 0
  SurfModel%Info(iSpec)%MeanProbDes  = 0.
  SurfModel%Info(iSpec)%MeanProbDesCount  = 0
END DO
#if (PP_TimeDiscMethod==42)
DO iSpec = 1,nSpecies
  SurfModel%Info(iSpec)%WallSpecNumCount = 0
  SurfModel%Info(iSpec)%Accomodation = 0
END DO
! initialize specific variables for analyze functionality like TPD and rate analysis
Adsorption%LateralInactive = GETLOGICAL('Surface-Adsorption-LateralInactive','.FALSE.')
Adsorption%CoverageReduction = GETLOGICAL('Surface-Adsorption-CoverageReduction','.FALSE.')
Adsorption%TPD = GETLOGICAL('Surface-Adsorption-doTPD','.FALSE.')
Adsorption%TPD_beta = GETREAL('Surface-Adsorption-TPD-Beta','0.')
Adsorption%TPD_Temp = 0.
#endif

! initialize species data
IF (.NOT.ALLOCATED(SpecSurf)) ALLOCATE(SpecSurf(nSpecies))
DO iSpec = 1,nSpecies
  WRITE(UNIT=hilf,FMT='(I0)') iSpec
  SpecSurf(iSpec)%ParamAntoine(1:3) = GETREALARRAY('Part-Species'//TRIM(hilf)//'-ParamAntoine',3,'0. , 0. , 0.')
  SpecSurf(iSpec)%condensCase = GETINT('Part-Species'//TRIM(hilf)//'-condensCase')
  SpecSurf(iSpec)%liquidTkrit = GETREAL('Part-Species'//TRIM(hilf)//'-liquidTkrit')
  SpecSurf(iSpec)%liquidTmelt = GETREAL('Part-Species'//TRIM(hilf)//'-liquidTmelt')
  SELECT CASE (SpecSurf(iSpec)%condensCase)
  CASE (1)
    SpecSurf(iSpec)%liquidAlpha = GETREAL('Part-Species'//TRIM(hilf)//'-liquidAlpha')
    SpecSurf(iSpec)%liquidBeta = GETREAL('Part-Species'//TRIM(hilf)//'-liquidBeta')
  CASE (2)
    SpecSurf(iSpec)%liquidBetaCoeff = GETREALARRAY('Part-Species'//TRIM(hilf)//'-liquidBetaCoeff',6,'0.,0.,0.,0.,0.,0.')
  CASE DEFAULT
    CALL abort(&
__STAMP__&
,'condensation case must be 1 or 2 for species: ',iSpec)
  END SELECT
END DO

! initialize model specific variables
DO iSpec = 1,nSpecies
  WRITE(UNIT=hilf,FMT='(I0)') iSpec
  DO iPartBound=1,nPartBound
    IF(.NOT.PartBound%Reactive(iPartBound)) CYCLE
    WRITE(UNIT=hilf2,FMT='(I0)') iPartBound
    hilf2=TRIM(hilf)//'-PartBound'//TRIM(hilf2)
    SELECT CASE(PartBound%SurfaceModel(iPartBound))
!-----------------------------------------------------------------------------------------------------------------------------------
    CASE(1)
!-----------------------------------------------------------------------------------------------------------------------------------
      CALL abort(&
__STAMP__&
,'surfacemode=1 not working!')
    !  ALLOCATE( Adsorption%MaxCoverage(1:SurfMesh%nMasterSides,1:nSpecies),&
    !            Adsorption%InitStick(1:SurfMesh%nMasterSides,1:nSpecies),&
    !            Adsorption%PrefactorStick(1:SurfMesh%nMasterSides,1:nSpecies),&
    !            Adsorption%Adsorbexp(1:SurfMesh%nMasterSides,1:nSpecies),&
    !            Adsorption%Nu_a(1:SurfMesh%nMasterSides,1:nSpecies),&
    !            Adsorption%Nu_b(1:SurfMesh%nMasterSides,1:nSpecies),&
    !            Adsorption%DesorbEnergy(1:SurfMesh%nMasterSides,1:nSpecies),&
    !            Adsorption%Intensification(1:SurfMesh%nMasterSides,1:nSpecies))
    !DO iSpec = 1,nSpecies
    !  WRITE(UNIT=hilf,FMT='(I0)') iSpec
    !  Adsorption%MaxCoverage(:,iSpec)     = GETREAL('Part-Species'//TRIM(hilf)//'-MaximumCoverage','0.')
    !  Adsorption%InitStick(:,iSpec)       = GETREAL('Part-Species'//TRIM(hilf)//'-InitialStick','0.')
    !  Adsorption%PrefactorStick(:,iSpec)  = GETREAL('Part-Species'//TRIM(hilf)//'-PrefactorStick','0.')
    !  Adsorption%Adsorbexp(:,iSpec)       = GETINT('Part-Species'//TRIM(hilf)//'-Adsorbexp','1')
    !  Adsorption%Nu_a(:,iSpec)            = GETREAL('Part-Species'//TRIM(hilf)//'-Nu-a','0.')
    !  Adsorption%Nu_b(:,iSpec)            = GETREAL('Part-Species'//TRIM(hilf)//'-Nu-b','0.')
    !  Adsorption%DesorbEnergy(:,iSpec)    = GETREAL('Part-Species'//TRIM(hilf)//'-Desorption-Energy-K','1.')
    !  Adsorption%Intensification(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Intensification-K','0.')
    !END DO
!-----------------------------------------------------------------------------------------------------------------------------------
    CASE(2)
!-----------------------------------------------------------------------------------------------------------------------------------
      IF (.NOT.ALLOCATED(Adsorption%ResultSpec)) ALLOCATE( Adsorption%ResultSpec(1:nPartBound,1:nSpecies))
      IF (.NOT.ALLOCATED(Adsorption%ReactCoeff)) ALLOCATE( Adsorption%ReactCoeff(1:nPartBound,1:nSpecies))
      IF (.NOT.ALLOCATED(Adsorption%ReactAccomodation)) ALLOCATE( Adsorption%ReactAccomodation(1:nPartBound,1:nSpecies))
      Adsorption%ResultSpec(iPartBound,iSpec)        = GETINT('Part-Species'//TRIM(hilf2)//'-ResultSpec')
      Adsorption%ReactCoeff(iPartBound,iSpec)        = GETREAL('Part-Species'//TRIM(hilf2)//'-RecombinationCoeff')
      Adsorption%ReactAccomodation(iPartBound,iSpec) = GETREAL('Part-Species'//TRIM(hilf2)//'-ReactionAccomodation')
      IF ((Adsorption%ResultSpec(iPartBound,iSpec).EQ.-1).AND.(Adsorption%ReactCoeff(iPartBound,iSpec).NE.0.)) THEN
        CALL abort(&
__STAMP__,&
'Resulting species for species '//TRIM(hilf)//' not defined although recombination coefficient .GT. 0')
      END IF
      IF (Adsorption%ResultSpec(iPartBound,iSpec).EQ.iSpec) THEN
        CALL abort(&
__STAMP__,&
'Resulting species for species '//TRIM(hilf)//' equal to incident species not possible for surfacemodel=2')
      END IF
!-----------------------------------------------------------------------------------------------------------------------------------
    CASE(3)
!-----------------------------------------------------------------------------------------------------------------------------------
      IF (.NOT.ALLOCATED(Adsorption%HeatOfAdsZero)) ALLOCATE( Adsorption%HeatOfAdsZero(1:nPartBound,1:nSpecies))
      IF (.NOT.ALLOCATED(Adsorption%Coordination)) ALLOCATE( Adsorption%Coordination(1:nPartBound,1:nSpecies))
      IF (.NOT.ALLOCATED(Adsorption%DiCoord)) ALLOCATE( Adsorption%DiCoord(1:nPartBound,1:nSpecies))
      IF (.NOT.ALLOCATED(Adsorption%ReactAccomodation)) ALLOCATE( Adsorption%ReactAccomodation(1:nPartBound,1:nSpecies))
      Adsorption%Coordination(iPartBound,iSpec) = GETINT('Part-Species'//TRIM(hilf2)//'-Coordination')
      Adsorption%DiCoord(iPartBound,iSpec) = GETINT('Part-Species'//TRIM(hilf2)//'-DiCoordination')
      ! check posibilities of coodrination of dicoord pairing. some pairing unphysical
      IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
        SELECT CASE(Adsorption%Coordination(iPartBound,iSpec))
        CASE(1)
          SELECT CASE(Adsorption%DiCoord(iPartBound,iSpec))
          CASE(1,2,3)
          CASE DEFAULT
            CALL abort(&
__STAMP__&
,"ERROR in INIT: for molecule at coordination=1 only dicoord 1,2,3 possible. Wrong dicoord for species:",iSpec)
          END SELECT
        CASE(2)
          IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
            SELECT CASE(Adsorption%DiCoord(iPartBound,iSpec))
            CASE(1,2,3,4,7)
            CASE DEFAULT
              CALL abort(&
__STAMP__&
,"ERROR in INIT: for molecule at coordination=2 only dicoord 1,2,3,4,7 possible. Wrong dicoord for species:",iSpec)
            END SELECT
          ELSE
            SELECT CASE(Adsorption%DiCoord(iPartBound,iSpec))
            CASE(1,2,3,4,6)
            CASE DEFAULT
              CALL abort(&
__STAMP__&
,"ERROR in INIT: for molecule at coordination=2 only dicoord 1,2,3,4,6 possible. Wrong dicoord for species:",iSpec)
            END SELECT
          END IF
        CASE(3)
          SELECT CASE(Adsorption%DiCoord(iPartBound,iSpec))
          CASE(1,2,3,5)
          CASE DEFAULT
            CALL abort(&
__STAMP__&
,"ERROR in INIT: for molecule at coordination=3 only dicoord 1,2,3,5 possible. Wrong dicoord for species:",iSpec)
          END SELECT
        END SELECT
      END IF
      Adsorption%HeatOfAdsZero(iPartbound,iSpec) = GETREAL('Part-Species'//TRIM(hilf2)//'-HeatOfAdsorption-K')
      Adsorption%ReactAccomodation(iPartBound,iSpec) = GETREAL('Part-Species'//TRIM(hilf2)//'-ReactionAccomodation')
      IF (Adsorption%Coordination(iPartBound,iSpec).EQ.0)THEN
        WRITE(UNIT=hilf2,FMT='(I0)') iPartBound
        CALL abort(&
__STAMP__,&
'Coordination of Species '//TRIM(hilf)//' for catalytic particle boundary '//TRIM(hilf2)//' not defined')
      END IF

!-----------------------------------------------------------------------------------------------------------------------------------
    CASE(5,6)
!-----------------------------------------------------------------------------------------------------------------------------------
      IF (.NOT.ALLOCATED(Adsorption%ResultSpec)) ALLOCATE( Adsorption%ResultSpec(1:nPartBound,1:nSpecies))
      Adsorption%ResultSpec(iPartBound,iSpec) = GETINT('Part-Species'//TRIM(hilf2)//'-ResultSpec')

!-----------------------------------------------------------------------------------------------------------------------------------
    CASE(101,102)
!-----------------------------------------------------------------------------------------------------------------------------------
      IF (.NOT.ALLOCATED(Adsorption%ResultSpec)) ALLOCATE( Adsorption%ResultSpec(1:nPartBound,1:nSpecies))
      IF (.NOT.ALLOCATED(Adsorption%ReactCoeff)) ALLOCATE( Adsorption%ReactCoeff(1:nPartBound,1:nSpecies))
      IF (.NOT.ALLOCATED(Adsorption%ReactAccomodation)) ALLOCATE( Adsorption%ReactAccomodation(1:nPartBound,1:nSpecies))
      IF (.NOT.ALLOCATED(Adsorption%ReactEnergy)) ALLOCATE( Adsorption%ReactEnergy(1:nPartBound,1:nSpecies))
      Adsorption%ResultSpec(iPartBound,iSpec) = 0
      Adsorption%ReactEnergy(iPartBound,iSpec) = 0
      IF (PartBound%SurfaceModel(iPartBound).EQ.101) THEN
        Adsorption%ReactCoeff(iPartBound,iSpec) = GETREAL('Part-Species'//TRIM(hilf2)//'-CondensationCoeff')
      ELSE
        Adsorption%ReactCoeff(iPartBound,iSpec) = 1.
      END IF
      Adsorption%ReactAccomodation(iPartBound,iSpec) = 1.
      IF (.NOT.ALLOCATED(Adsorption%SurfaceSpec)) ALLOCATE(Adsorption%SurfaceSpec(1:nPartBound,1:nSpecies))
      Adsorption%SurfaceSpec(iPartBound,iSpec) = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-LiquidSpec')
      ! check parameters used for evaporation pressure of Antoine Eq
      IF (Adsorption%SurfaceSpec(iPartBound,iSpec) .AND. (ALMOSTZERO(SpecSurf(iSpec)%ParamAntoine(1))) &
           .AND. (ALMOSTZERO(SpecSurf(iSpec)%ParamAntoine(2))) &
           .AND. (ALMOSTZERO(SpecSurf(iSpec)%ParamAntoine(3))) ) THEN
        CALL abort(&
__STAMP__&
,'Antoine Parameters not defined for species: ',iPartBound)
      END IF
!-----------------------------------------------------------------------------------------------------------------------------------
    END SELECT
!-----------------------------------------------------------------------------------------------------------------------------------
  END DO
END DO


! allocate and initialize adsorption variables
ALLOCATE( SurfModel%SumEvapPart(1:nSurfSample,1:nSurfSample,1:SurfMesh%nMasterSides,1:nSpecies),&
          SurfModel%SumDesorbPart(1:nSurfSample,1:nSurfSample,1:SurfMesh%nMasterSides,1:nSpecies),&
          SurfModel%SumReactPart(1:nSurfSample,1:nSurfSample,1:SurfMesh%nMasterSides,1:nSpecies),&
          SurfModel%SumAdsorbPart(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides,1:nSpecies),&
          SurfModel%SumERDesorbed(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides,1:nSpecies))
ALLOCATE( Adsorption%ProbAds(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides,1:nSpecies),&
          Adsorption%ProbDes(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides,1:nSpecies),&
          Adsorption%Coverage(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides,1:nSpecies))
ALLOCATE( Adsorption%DensSurfAtoms(1:SurfMesh%nTotalSides),&
          Adsorption%AreaIncrease(1:SurfMesh%nTotalSides),&
          Adsorption%CrystalIndx(1:SurfMesh%nTotalSides))

IF (SurfMesh%SurfOnProc) THEN
  SurfModel%SumEvapPart(:,:,:,:) = 0
  SurfModel%SumDesorbPart(:,:,:,:) = 0
  SurfModel%SumAdsorbPart(:,:,:,:) = 0
  SurfModel%SumReactPart(:,:,:,:)  = 0
  SurfModel%SumERDesorbed(:,:,:,:) = 0
  Adsorption%ProbAds(:,:,:,:) = 0.
  Adsorption%ProbDes(:,:,:,:) = 0.
END IF ! SurfMesh%SurfOnProc

! Initialize surface properties from particle boundary values
Adsorption%DensSurfAtoms(:) = 0
Adsorption%AreaIncrease(:)  = 0
Adsorption%CrystalIndx(:)   = 0
DO iSide=1,SurfMesh%nTotalSides
  SideID = SurfMesh%SurfIDToSideID(iSide)
  PartboundID = PartBound%MapToPartBC(BC(SideID))
  IF (.NOT.PartBound%Reactive(PartboundID)) CYCLE
  IF (PartBound%SolidState(PartBoundID)) THEN
    Adsorption%AreaIncrease(iSide)  = PartBound%SolidAreaIncrease(PartBoundID)
    Adsorption%CrystalIndx(iSide)   = PartBound%SolidCrystalIndx(PartBoundID)
    Adsorption%DensSurfAtoms(iSide) = PartBound%SolidPartDens(PartBoundID)*Adsorption%AreaIncrease(iSide)
  ELSE
  END IF
END DO
Adsorption%NumCovSamples = 0

ALLOCATE ( Adsorption%IncidentNormalVeloAtSurf(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides,1:nSpecies),&
           Adsorption%SurfaceNormalVelo(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides,1:nSpecies),&
           Adsorption%CollSpecPartNum(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides,1:nSpecies) )
DO iSpec = 1,nSpecies
  ! Expacted value for an assumed Rayleigh distribution
  Adsorption%IncidentNormalVeloAtSurf(:,:,:,iSpec) =  &
      SQRT(BoltzmannConst*Species(iSpec)%Init(0)%MWTemperatureIC/Species(iSpec)%MassIC) *SQRT(PI/2.)
END DO
Adsorption%SurfaceNormalVelo(:,:,:,:)  = 0.
Adsorption%CollSpecPartNum(:,:,:,:)    = 0

! Initialize surface coverage
CALL InitSurfCoverage()

#if USE_MPI
IF (SurfMesh%SurfOnProc) CALL InitSurfModel_MPI()
#endif /*USE_MPI*/

CALL InitSMCR()
CALL InitSurfChem()
IF (WriteMacroSurfaceValues.OR.DSMC%CalcSurfaceVal) THEN
  IF(SurfMesh%SurfOnProc) CALL Init_SurfChemistrySampling()
END IF

END SUBROUTINE InitSurfaceModel


SUBROUTINE InitSurfChem()
!===================================================================================================================================
!> Initializing surface reaction variables
!===================================================================================================================================
! MODULES
USE MOD_Globals                ,ONLY: abort, MPIRoot, UNIT_StdOut
USE MOD_DSMC_Vars              ,ONLY: SpecDSMC
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_ReadInTools            ,ONLY: GETREAL, GETINT, GETREALARRAY, GETINTARRAY, PrintOption
#if (PP_TimeDiscMethod==42)
USE MOD_SurfaceModel_Vars      ,ONLY: SurfModel
#endif
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
INTEGER                          :: MaxDissNum, MaxRecombNum, MaxReactNum
INTEGER , ALLOCATABLE            :: SpecNRecombReact(:)
INTEGER                          :: nDissoc, nExch
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)')' INIT SURFACE CHEMISTRY...'

! total number of all dossiciation reactions
Adsorption%NumOfDissocReact = 0
Adsorption%NumOfRecombReact = 0
Adsorption%NumOfExchReact = 0

#if !(USE_LOADBALANCE)
IF (SurfMesh%SurfOnProc .OR. MPIRoot) THEN
#endif

#if (PP_TimeDiscMethod==42)
  ALLOCATE( SurfModel%ProperInfo(1:nSpecies))
#endif

  MaxDissNum = GETINT('Surface-MaxDissNum','0')
  MaxRecombNum = MaxDissNum

  ! allocate and initialize dissociative and associative reactions species map
  IF ( (MaxDissNum.GT.0) .OR. (MaxRecombNum.GT.0) ) THEN
    ALLOCATE( Adsorption%DissocReact(1:2,1:MaxDissNum,1:nSpecies))
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
,'Error in Init_SurfChem: wrong Product species for reaction '//TRIM(hilf2)//'! Species > nSpecies!')
        END IF
        IF ((Adsorption%DissocReact(1,iReactNum,iSpec).LT.0) &
          .OR.(Adsorption%DissocReact(2,iReactNum,iSpec).LT.0) ) THEN
          CALL abort(&
__STAMP__&
,'Error in Init_SurfChem: Product species for reaction '//TRIM(hilf2)//' negative!')
        END IF
      END DO
    END DO

    ! find the max number of recombination reactions for each species from dissociation reaction mapping
    ALLOCATE(SpecNRecombReact(1:nSpecies))
    SpecNRecombReact(:) = 0
    DO iSpec = 1,nSpecies
      DO iSpec2 = 1,nSpecies ; DO iReactNum = 1,MaxDissNum
        IF ((Adsorption%DissocReact(1,iReactNum,iSpec2).EQ.iSpec).OR.(Adsorption%DissocReact(2,iReactNum,iSpec2).EQ.iSpec) ) THEN
          SpecNRecombReact(iSpec) = SpecNRecombReact(iSpec) + 1
        END IF
      END DO ; END DO
    END DO
    MaxRecombNum = MAXVAL(SpecNRecombReact)
    Adsorption%NumOfRecombReact = SUM(SpecNRecombReact(:))
    Adsorption%RecombNum = MaxRecombNum
    DEALLOCATE(SpecNRecombReact)

    ! fill recombination reactions species map from defined dissociative reactions
    MaxReactNum = MaxDissNum + MaxRecombNum
    ALLOCATE( Adsorption%RecombReact(1:2,1:MaxRecombNum,1:nSpecies),&
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
          IF (ALLOCATED(Adsorption%DiCoord) .AND. (Adsorption%EDissBondAdsorbPoly(0,iSpec).EQ.0)) THEN
            CALL abort(&
__STAMP__&
,'Error in Init_SurfChem: Adsorption dissociation bond energy of dicoord for species '//TRIM(hilf)//' not defined!')
          END IF
        ELSE
          IF (ALLOCATED(Adsorption%DiCoord) .AND.(Adsorption%EDissBond(0,iSpec).EQ.0.))THEN
            CALL abort(&
__STAMP__&
,'Error in Init_SurfChem: Adsorption dissociation bond energy for species '//TRIM(hilf)//' not defined!')
          END IF
        END IF
      END IF
    END DO
    DO iSpec = 1,nSpecies
      ReactNum = 1
      DO iSpec2 = 1,nSpecies
        DO iReactNum2 = 1,MaxDissNum
          IF (Adsorption%DissocReact(1,iReactNum2,iSpec2).EQ.iSpec) THEN
            Adsorption%RecombReact(1,ReactNum,iSpec) = Adsorption%DissocReact(2,iReactNum2,iSpec2)
            Adsorption%RecombReact(2,ReactNum,iSpec) = iSpec2
            Adsorption%EDissBond((MaxDissNum+ReactNum),iSpec) = Adsorption%EDissBond(iReactNum2,iSpec2)
            ReactNum = ReactNum + 1
          ELSE IF (Adsorption%DissocReact(2,iReactNum2,iSpec2).EQ.iSpec) THEN
            Adsorption%RecombReact(1,ReactNum,iSpec) = Adsorption%DissocReact(1,iReactNum2,iSpec2)
            Adsorption%RecombReact(2,ReactNum,iSpec) = iSpec2
            Adsorption%EDissBond((MaxDissNum+ReactNum),iSpec) = Adsorption%EDissBond(iReactNum2,iSpec2)
            ReactNum = ReactNum + 1
          ELSE
            CYCLE
          END IF
          CALL PrintOption('Recombination reaction: ','OUTPUT',IntOpt=ReactNum-1)
          CALL PrintOption('For species: ','OUTPUT',IntOpt=iSpec)
          CALL PrintOption('PartnerSpec: ','OUTPUT',IntOpt=Adsorption%RecombReact(1,ReactNum-1,iSpec))
          CALL PrintOption('ResultSpec: ','OUTPUT',IntOpt=Adsorption%RecombReact(2,ReactNum-1,iSpec))
          CALL PrintOption('DissBondEnergy: ','OUTPUT',RealOpt=Adsorption%EDissBond((MaxDissNum+ReactNum-1),iSpec))
        END DO
      END DO
      IF (ReactNum.LE.(MaxRecombNum)) THEN
        Adsorption%RecombReact(:,ReactNum:(MaxReactNum-MaxDissNum),iSpec) = 0
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
    nExch = GETINT('Surface-Nbr-ExchangeReactions','0')
    ! Allocate and fill one array for all types of reactions (dissociation, recombination, exchange reaction)
    ALLOCATE( Adsorption%ChemReactant(1:2,1:nDissoc+nExch),&
              Adsorption%ChemProduct(1:2,1:nDissoc+nExch),&
              Adsorption%Reactant_DissBond_K(1:2,1:nDissoc+nExch),&
              Adsorption%Product_DissBond_K(1:2,1:nDissoc+nExch))
    ! Initialize allocated variables
    Adsorption%ChemReactant(1:2,1:nDissoc+nExch)=0
    Adsorption%ChemProduct(1:2,1:nDissoc+nExch)=0
    Adsorption%Reactant_DissBond_K(1:2,1:nDissoc+nExch)=0.
    Adsorption%Product_DissBond_K(1:2,1:nDissoc+nExch)=0.
    ! fill dissociation reactions (can also be used for recombination)
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
    DO iReactNum = 1,nExch
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
    nExch = 0
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
          IF (ALLOCATED(Adsorption%DiCoord) .AND. (Adsorption%EDissBondAdsorbPoly(0,iSpec).EQ.0)) THEN
            CALL abort(&
__STAMP__&
,'Error in Init_SurfChem: Adsorption dissociation bond energy of dicoord for species '//TRIM(hilf)//' not defined!')
          END IF
        ELSE
          IF (ALLOCATED(Adsorption%DiCoord) .AND. (Adsorption%EDissBond(0,iSpec).EQ.0.)) THEN
            CALL abort(&
__STAMP__&
,'Error in Init_SurfChem: Adsorption dissociation bond energy for species '//TRIM(hilf)//' not defined!')
          END IF
        END IF
      END IF
    END DO
  END IF !MaxDissNum > 0
  ! save defined number of surface reactions
  Adsorption%DissNum = MaxDissNum
  Adsorption%RecombNum = MaxRecombNum
  Adsorption%ReactNum = MaxReactNum
  Adsorption%nDissocReactions = nDissoc
  Adsorption%nExchReactions = nExch
  Adsorption%NumOfExchReact = nExch

#if (PP_TimeDiscMethod==42)
! allocate and initialize analyze arrays for surface reaction rate/processes analyze
  DO iSpec=1,nSpecies
    ALLOCATE( SurfModel%ProperInfo(iSpec)%NumAdsReact(1:Adsorption%ReactNum+1),&
              SurfModel%ProperInfo(iSpec)%MeanAdsActE(1:Adsorption%ReactNum+1),&
              SurfModel%ProperInfo(iSpec)%ProperAdsActE(1:Adsorption%ReactNum+1),&
              SurfModel%ProperInfo(iSpec)%MeanAdsnu(1:Adsorption%ReactNum+1),&
              SurfModel%ProperInfo(iSpec)%ProperAdsnu(1:Adsorption%ReactNum+1),&
              SurfModel%ProperInfo(iSpec)%AdsReactCount(1:Adsorption%ReactNum+1),&
              SurfModel%ProperInfo(iSpec)%ProperAdsReactCount(1:Adsorption%ReactNum+1),&
              SurfModel%ProperInfo(iSpec)%NumSurfReact(1:Adsorption%ReactNum+Adsorption%NumOfExchReact+1),&
              SurfModel%ProperInfo(iSpec)%MeanSurfActE(1:Adsorption%ReactNum+Adsorption%NumOfExchReact+1),&
              SurfModel%ProperInfo(iSpec)%ProperSurfActE(1:Adsorption%ReactNum+Adsorption%NumOfExchReact+1),&
              SurfModel%ProperInfo(iSpec)%MeanSurfnu(1:Adsorption%ReactNum+Adsorption%NumOfExchReact+1),&
              SurfModel%ProperInfo(iSpec)%ProperSurfnu(1:Adsorption%ReactNum+Adsorption%NumOfExchReact+1),&
              SurfModel%ProperInfo(iSpec)%SurfReactCount(1:Adsorption%ReactNum+Adsorption%NumOfExchReact+1),&
              SurfModel%ProperInfo(iSpec)%ProperSurfReactCount(1:Adsorption%ReactNum+Adsorption%NumOfExchReact+1),&
              SurfModel%ProperInfo(iSpec)%HeatFluxAdsCount(1:Adsorption%ReactNum+1),&
              SurfModel%ProperInfo(iSpec)%HeatFluxDesCount(1:Adsorption%ReactNum+Adsorption%NumOfExchReact+1),&
              SurfModel%ProperInfo(iSpec)%HeatFlux(1:2))
    SurfModel%ProperInfo(iSpec)%NumAdsReact(:)         = 0.
    SurfModel%ProperInfo(iSpec)%MeanAdsActE(:)         = 0.
    SurfModel%ProperInfo(iSpec)%ProperAdsActE(:)       = 0.
    SurfModel%ProperInfo(iSpec)%MeanAdsnu(:)           = 0.
    SurfModel%ProperInfo(iSpec)%ProperAdsnu(:)         = 0.
    SurfModel%ProperInfo(iSpec)%AdsReactCount(:)       = 0
    SurfModel%ProperInfo(iSpec)%ProperAdsReactCount(:) = 0
    SurfModel%ProperInfo(iSpec)%NumSurfReact(:)         = 0.
    SurfModel%ProperInfo(iSpec)%MeanSurfActE(:)         = 0.
    SurfModel%ProperInfo(iSpec)%ProperSurfActE(:)       = 0.
    SurfModel%ProperInfo(iSpec)%MeanSurfnu(:)           = 0.
    SurfModel%ProperInfo(iSpec)%ProperSurfnu(:)         = 0.
    SurfModel%ProperInfo(iSpec)%SurfReactCount(:)       = 0
    SurfModel%ProperInfo(iSpec)%ProperSurfReactCount(:) = 0
    SurfModel%ProperInfo(iSpec)%HeatFluxAdsCount(:) = 0.
    SurfModel%ProperInfo(iSpec)%HeatFluxDesCount(:) = 0.
    SurfModel%ProperInfo(iSpec)%HeatFlux(:)        = 0.
  END DO
#endif
  CALL Init_TST_Coeff()
#if !(USE_LOADBALANCE)
END IF ! SurfMesh%SurfOnProc .OR. MPIRoot
#endif

SWRITE(UNIT_stdOut,'(A)')' INIT SURFACE CHEMISTRY DONE!'

END SUBROUTINE InitSurfChem


SUBROUTINE Init_TST_Coeff()
!===================================================================================================================================
!> Initializion of the Transition state theory (TST) factors for surface chemistry
!> nu_react_1 = prefactor * T^(powerfactor)
!> nu_react_2 => calculated with partitionfunctions of activated state Q* and ground state Q_a => Q*/Q_a
!===================================================================================================================================
! MODULES
USE MOD_Globals                ,ONLY: abort, UNIT_StdOut
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_ReadInTools            ,ONLY: GETREAL, GETINT
#if USE_MPI
USE MOD_Globals                ,ONLY: MPIRoot
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL                            :: PartitionArraySize
CHARACTER(32)                   :: hilf, hilf2
INTEGER                         :: TST_Case
INTEGER                         :: iSpec, iReactNum, DissocID, RecombID
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)')' INIT SURFACE TST REACTION COEFFICIENTS!'

ALLOCATE(Adsorption%TST_Calc(0:Adsorption%ReactNum,1:nSpecies))
Adsorption%TST_Calc(:,:) = .FALSE.
TST_Case = GETINT('Surface-Adsorption-CalcTST')

SELECT CASE (TST_Case)
CASE (0,1)
  ! reaction rate coefficients
  ALLOCATE( Adsorption%Ads_Powerfactor(1:nSpecies),&
            Adsorption%Ads_Prefactor(1:nSpecies),&
            Adsorption%Diss_Powerfactor(1:Adsorption%DissNum,1:nSpecies),&
            Adsorption%Diss_Prefactor(1:Adsorption%DissNum,1:nSpecies),&
            Adsorption%ER_Powerfactor(1:Adsorption%RecombNum,1:nSpecies),&
            Adsorption%ER_Prefactor(1:Adsorption%RecombNum,1:nSpecies))
  DO iSpec = 1,nSpecies
    WRITE(UNIT=hilf,FMT='(I0)') iSpec
    DO iReactNum=0,Adsorption%ReactNum
      IF (iReactNum.EQ.0) THEN
        Adsorption%Ads_Powerfactor(iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Adsorption-Powerfactor')
        Adsorption%Ads_Prefactor(iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Adsorption-Prefactor')
        IF ((TST_Case.EQ.1) .AND. (Adsorption%Ads_Prefactor(iSpec).EQ.-1.) .AND. (Adsorption%Ads_Powerfactor(iSpec).EQ.-1.)) THEN
          Adsorption%TST_Calc(iReactNum,iSpec) = .TRUE.
        END IF
      ELSE IF ((iReactNum.GT.0) .AND. (iReactNum.LE.Adsorption%DissNum)) THEN
        DissocID = iReactNum
        WRITE(UNIT=hilf2,FMT='(I0)') DissocID
        Adsorption%Diss_Powerfactor(DissocID,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-SurfDiss'//TRIM(hilf2)//'-Powerfactor')
        Adsorption%Diss_Prefactor(DissocID,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-SurfDiss'//TRIM(hilf2)//'-Prefactor')
        IF (TST_Case.EQ.1) THEN
          IF ((Adsorption%Diss_Prefactor(DissocID,iSpec).EQ.-1.) .AND. (Adsorption%Diss_Powerfactor(DissocID,iSpec).EQ.-1.)) THEN
            Adsorption%TST_Calc(iReactNum,iSpec) = .TRUE.
          END IF
        END IF
      ELSE IF ((iReactNum.GT.0) .AND. (iReactNum.GT.Adsorption%DissNum)) THEN
        RecombID = iReactNum-Adsorption%DissNum
        WRITE(UNIT=hilf2,FMT='(I0)') RecombID
        Adsorption%ER_Powerfactor(RecombID,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Surf-ER'//TRIM(hilf2)//'-Powerfactor')
        Adsorption%ER_Prefactor(RecombID,iSpec) =   GETREAL('Part-Species'//TRIM(hilf)//'-Surf-ER'//TRIM(hilf2)//'-Prefactor')
        IF (TST_Case.EQ.1) THEN
          IF ((Adsorption%ER_Prefactor(RecombID,iSpec).EQ.-1.) .AND. (Adsorption%ER_Powerfactor(RecombID,iSpec).EQ.-1.)) THEN
            Adsorption%TST_Calc(iReactNum,iSpec) = .TRUE.
          END IF
        END IF
      END IF
    END DO
  END DO
CASE (2)
  Adsorption%TST_Calc(:,:) = .TRUE.
CASE DEFAULT
  CALL abort(&
__STAMP__&
,'Surface SMCR init ERROR: TST coefficient has to be 0,1 or 2:',TST_Case)
END SELECT

SWRITE(UNIT_stdOut,'(A)')' INIT SURFACE TST REACTION COEFFICIENTS DONE!'
END SUBROUTINE Init_TST_Coeff


SUBROUTINE InitSurfCoverage()
!===================================================================================================================================
!> check if restart is done and if surface data is given in state file
!> if not then
!> Surface coverage is initialized either from given surface init file or from constant values for each given particle boundary.
!> Therefore, the ini file is checked if state file is specified. If not coverage is initialized from parameter else
!> .h5 file is read and used for init. If file has wrong entries programm aborts.
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_HDF5_INPUT             ,ONLY: DatasetExists,GetDataProps,ReadAttribute,ReadArray,GetDataSize
USE MOD_Restart_Vars           ,ONLY: DoRestart,RestartFile
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_ReadInTools            ,ONLY: GETSTR,GETREAL
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, nPartBound, PartBound
USE MOD_Particle_Boundary_Vars ,ONLY: offSetSurfSide, nSurfBC
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
USE MOD_ReadInTools            ,ONLY: PrintOption
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
CHARACTER(LEN=255)               :: SurfaceFileName, Type_HDF5, NodeType_HDF5
CHARACTER(32)                    :: hilf, hilf2
INTEGER                          :: iSpec, iSurfSide, iPartBound, iSubSurf, jSubSurf, iName, iVar
INTEGER                          :: SideID, PartBoundID
REAL , ALLOCATABLE               :: Coverage_tmp(:,:)
REAL , ALLOCATABLE               :: SurfState_HDF5(:,:,:,:)
LOGICAL                          :: exists
INTEGER                          :: nSpecies_HDF5, nVarSurf_HDF5, nSurfSides_HDF5, N_HDF5, nSurfBC_HDF5
INTEGER                          :: nVar2D, nVar2D_Spec, nVar2D_Total
CHARACTER(LEN=255),ALLOCATABLE   :: SurfBCName_HDF5(:)
REAL                             :: Version_HDF5
LOGICAL                          :: WallmodelExists(1:nPartBound), SurfModelTypeExists
INTEGER,ALLOCATABLE              :: SurfModelType(:)
!===================================================================================================================================

IF (DoRestart) THEN
  ! check if datasets for restarting from state exists in state file used for restart
  WallModelExists(:)=.TRUE.
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=PartMPI%COMM)!MPI_COMM_WORLD)
  SurfModelTypeExists=.FALSE.
  CALL DatasetExists(File_ID,'SurfModelType',SurfModelTypeExists)
  IF (SurfModelTypeExists) THEN
    SWRITE(UNIT_stdOut,'(A,A)')' GET NUMBER OF SURFACE-SIDES IN RESTART FILE... '
    CALL GetDataProps('SurfModelType',nVarSurf_HDF5,N_HDF5,nSurfSides_HDF5,NodeType_HDF5)
    SWRITE(UNIT_stdOut,'(A)')' DONE!'
    IF (nSurfSides_HDF5.NE.SurfMesh%nGlobalSides) THEN
      SWRITE(UNIT_stdOut,'(A,A)') ' NUMBER OF SURFACE-SIDES IN RESTART FILE NOT EQUAL TO CALCULATION ... RESTARTING FROM INI'
      WallModelExists(:)=.FALSE.
    ELSE
      ALLOCATE(SurfModelType(SurfMesh%nMasterSides))
      ! Associate construct for integer KIND=8 possibility
      ASSOCIATE (&
            nSides          => INT(SurfMesh%nMasterSides,IK) ,&
            nSpecies        => INT(nSpecies,IK) ,&
            offsetSurfSide  => INT(offsetSurfSide,IK) )
        CALL ReadArray('SurfaceModelType',1,(/nSides/) ,&
                       offsetSurfSide,1,IntegerArray_i4=SurfModelType)
      END ASSOCIATE
      DO iSurfSide=1,SurfMesh%nMasterSides
        SideID = SurfMesh%SurfIDToSideID(iSurfSide)
        PartboundID = PartBound%MapToPartBC(BC(SideID))
        IF (PartBound%SurfaceModel(PartboundID).NE.SurfModelType(iSurfSide)) THEN
          WallModelExists(PartBoundID)=.FALSE.
          EXIT
        END IF
      END DO
    END IF
  ELSE
    WallModelExists(:)=.FALSE.
  END IF
  CALL CloseDataFile()
  !DO iPartBound=1,nPartBound
  !  IF (WallModelExists(iPartBound)) IPWRITE(UNIT_StdOut, '(A,I0,A)')' Surfaces for particle boundary: ' &
  !    ,iPartBound', will be initialized from Restartfile'
  !END DO
ELSE
  WallModelExists(:)=.FALSE.
END IF ! DoRestart

IF (SurfMesh%SurfOnProc) THEN
  ! write zeros into global coverage array for each surface
  Adsorption%Coverage(:,:,:,:) = 0.
END IF

SurfaceFileName=GETSTR('Particles-SurfCoverageFile')
! If no surface file is given, initialize from ini parameters else use values from file
IF (TRIM(SurfaceFileName).EQ.'none') THEN
  SWRITE(UNIT_StdOut, '(A)')' INIT SURFACE COVERAGE FROM INI FILE VALUES...'

  ! initialize temporary coverage array for each particle boundary in case no surfacestatefile is used
  ALLOCATE(Coverage_tmp(1:nPartBound,1:nSpecies))
  Coverage_tmp(:,:) = 0.
  DO iSpec = 1,nSpecies
    WRITE(UNIT=hilf,FMT='(I0)') iSpec
    DO iPartBound=1,nPartBound
      WRITE(UNIT=hilf2,FMT='(I0)') iPartBound
      hilf2=TRIM(hilf)//'-PartBound'//TRIM(hilf2)
      Coverage_tmp(iPartBound,iSpec) = GETREAL('Part-Species'//TRIM(hilf2)//'-InitialCoverage','0.')
    END DO
  END DO
  IF (SurfMesh%SurfOnProc) THEN
    ! write temporary coverage values into global coverage array for each surface
    DO iSurfSide=1,SurfMesh%nMasterSides
      SideID = SurfMesh%SurfIDToSideID(iSurfSide)
      PartboundID = PartBound%MapToPartBC(BC(SideID))
      IF (.NOT.WallModelExists(PartBoundID)) THEN
        DO iSpec=1,nSpecies
          Adsorption%Coverage(:,:,iSurfSide,iSpec) = Coverage_tmp(PartBoundID,iSpec)
        END DO
      END IF
    END DO
  END IF
  SDEALLOCATE(Coverage_tmp)
  SWRITE(UNIT_StdOut, '(A)')' INIT SURFACE COVERAGE FROM INI FILE VALUES DONE!'
ELSE ! SurfaceFileName.EQ.'none'
  SWRITE(UNIT_StdOut, '(A)')' INIT SURFACE COVERAGE FROM '//TRIM(SurfaceFileName)//' ...'
  CALL OpenDataFile(TRIM(SurfaceFileName),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=PartMPI%COMM)!MPI_COMM_WORLD)
  exists=.FALSE.
  ! check if given file is of type 'DSMCSurfState'
  CALL DatasetExists(File_ID,'File_Type',exists,attrib=.TRUE.)
  IF (exists) THEN
    CALL ReadAttribute(File_ID,'File_Type',1,StrScalar=Type_HDF5)
    IF (TRIM(Type_HDF5).NE.'DSMCSurfState') CALL abort(&
__STAMP__&
,'Error in surface coverage init: filetype of given surfcoverage-initfile does not match!')
  ELSE
    CALL abort(&
__STAMP__&
,'Error in surface coverage init: attribute "filetype" does not exist!')
  END IF
  ! check if number of species is equal
  CALL DatasetExists(File_ID,'DSMC_nSpecies',exists,attrib=.TRUE.)
  IF (exists) THEN
    CALL ReadAttribute(File_ID,'DSMC_nSpecies',1,IntegerScalar=nSpecies_HDF5)
    IF (nSpecies_HDF5.NE.nSpecies) CALL abort(&
__STAMP__&
,'Error in surface coverage init: number of Species does not match!')
  ELSE
    CALL abort(&
__STAMP__&
,'Error in surface coverage init: attribute "nSpecies" does not exist!')
  END IF
  ! check for given file version
  CALL DatasetExists(File_ID,'File_Version',exists,attrib=.TRUE.)
  IF (exists) THEN
    CALL ReadAttribute(File_ID,'File_Version',1,RealScalar=Version_HDF5)
  ELSE
    CALL abort(&
__STAMP__&
,'Error in surface coverage init: attribute "fileversion" does not exist!')
  END IF

  ! check if Dataset SurfaceData exists and read from container
  CALL DatasetExists(File_ID,'SurfaceData',exists)
  IF (exists) THEN
    CALL GetDataProps('SurfaceData',nVarSurf_HDF5,N_HDF5,nSurfSides_HDF5,NodeType_HDF5)
    SWRITE(UNIT_stdOut,'(A,A)')' GET NUMBER AND NAMES OF SURFACE-BC-SIDES IN COVERAGE-INIT FILE... '
    CALL GetDataSize(File_ID,'BC_Surf',nDims,HSize,attrib=.TRUE.)
    nSurfBC_HDF5 = INT(HSize(1),4)
    CALL PrintOption('Number of Surface BCs','HDF5',IntOpt=nSurfBC_HDF5)
    IF (MPIROOT) THEN
      ALLOCATE(SurfBCName_HDF5(nSurfBC_HDF5))
      CALL ReadAttribute(File_ID,'BC_Surf',nSurfBC_HDF5,StrArray=SurfBCName_HDF5)
      DO iName = 1,nSurfBC_HDF5
        WRITE(UNIT=hilf,FMT='(I0)') iName
        CALL PrintOption('BC'//ADJUSTL(TRIM(hilf))//'Name','HDF5',StrOpt=SurfBCName_HDF5(iName))
      END DO
    END IF
    SWRITE(UNIT_stdOut,'(A)')' DONE!'
    IF (nSurfBC_HDF5.NE.nSurfBC) CALL abort(&
__STAMP__&
,'Error in surface coverage init: number of surface boundaries in HDF5-file does not match!')
    IF (nSurfSides_HDF5.NE.SurfMesh%nGlobalSides) CALL abort(&
__STAMP__&
,'Error in surface coverage init: number of surface sides in HDF5-file does not match!')
    ! number comes from boundary sampling, if surfacemodel is used then nVar2d_Spec=4
    IF (Version_HDF5.GT.0.1)THEN
      nVar2D = 10
    ELSE
      nVar2D = 5
    END IF
    nVar2D_Spec = 4
    nVar2D_Total = nVar2D + nVar2D_Spec*nSpecies
    IF (nVarSurf_HDF5.NE.nVar2D_Total) CALL abort(&
__STAMP__&
,'Error in surface coverage init: number of variables in HDF5-file does not match!')
    SDEALLOCATE(SurfState_HDF5)
    ALLOCATE(SurfState_HDF5(1:nVarSurf_HDF5,1:nSurfSample,1:nSurfSample,SurfMesh%nMasterSides))

    ! Associate construct for integer KIND=8 possibility
    ASSOCIATE (&
          nVarSurf_HDF5   => INT(nVarSurf_HDF5,IK)   ,&
          nSurfSample     => INT(nSurfSample,IK)     ,&
          nSides          => INT(SurfMesh%nMasterSides,IK) ,&
          offsetSurfSide  => INT(offsetSurfSide,IK)  )
      CALL ReadArray(  'SurfaceData',4,(/nVarSurf_HDF5,nSurfSample,nSurfSample,nSides/)&
                     ,offsetSurfSide,4,RealArray=SurfState_HDF5)
    END ASSOCIATE
    iVar = 3
    DO iSpec = 1, nSpecies
      DO iSurfSide = 1, SurfMesh%nMasterSides
        SideID = SurfMesh%SurfIDToSideID(iSurfSide)
        PartboundID = PartBound%MapToPartBC(BC(SideID))
        IF (PartBound%Reactive(PartboundID)) THEN
          DO jSubSurf = 1, nSurfSample
            DO iSubSurf = 1, nSurfSample
              Adsorption%Coverage(iSubSurf,jSubSurf,iSurfSide,iSpec) = SurfState_HDF5(iVar,iSubSurf,jSubSurf,iSurfSide)
            END DO ! iSubSurf = 1, nSurfSample
          END DO ! jSubSurf = 1, nSurfSample
        ELSE
          Adsorption%Coverage(:,:,iSurfSide,iSpec) = 0.
        END IF
      END DO ! iSurfSide = 1, SurfMesh%nMasterSides
      iVar = iVar + nVar2D_Spec
    END DO ! iSpec = 1, nSpecies
    SDEALLOCATE(SurfState_HDF5)
  ELSE
    CALL abort(&
__STAMP__&
,'Error in surface coverage init: dataset "SurfaceData" does not exist!')
  END IF
  CALL CloseDataFile()
  SWRITE(UNIT_StdOut, '(A)')' INIT SURFACE COVERAGE FROM '//TRIM(SurfaceFileName)//' DONE!'
END IF ! SurfaceFileName.EQ.'none'

END SUBROUTINE InitSurfCoverage


SUBROUTINE Init_SurfChemistrySampling()
!===================================================================================================================================
!> Initialization of additional surface sampling if any Partbound#-SurfaceModel>0
!===================================================================================================================================
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, SampWall
#if USE_MPI
USE MOD_Particle_Boundary_Vars ,ONLY: SurfCOMM
USE MOD_Particle_MPI_Vars      ,ONLY: SurfSendBuf,SurfRecvBuf,SurfExchange
#endif /*USE_MPI*/
!===================================================================================================================================
IMPLICIT NONE
!===================================================================================================================================
! Local variable declaration
INTEGER                          :: iSide
#if USE_MPI
INTEGER                          :: iProc, SendArraySize, RecvArraySize
#endif
!===================================================================================================================================

SurfMesh%ReactiveSampSize=(5+nSpecies+nSpecies+2*(Adsorption%ReactNum*nSpecies))
DO iSide=1,SurfMesh%nTotalSides ! caution: iSurfSideID
  ALLOCATE(SampWall(iSide)%SurfModelState(1:5+nSpecies,1:nSurfSample,1:nSurfSample))
  SampWall(iSide)%SurfModelState=0.
  ALLOCATE(SampWall(iSide)%Accomodation(1:nSpecies,1:nSurfSample,1:nSurfSample))
  SampWall(iSide)%Accomodation=0.
  ALLOCATE(SampWall(iSide)%SurfModelReactCount(1:2*Adsorption%ReactNum,1:nSpecies,1:nSurfSample,1:nSurfSample))
  SampWall(iSide)%SurfModelReactCount=0.
END DO
#if USE_MPI
! Reallocate buffer for mpi communication of sampling
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesSend(iProc).GT.0) THEN
    SendArraySize = SIZE(SurfSendBuf(iProc)%content,DIM=1,KIND=4)
    SDEALLOCATE(SurfSendBuf(iProc)%content)
    ALLOCATE(SurfSendBuf(iProc)%content(SendArraySize+SurfMesh%ReactiveSampSize*(nSurfSample**2)*SurfExchange%nSidesSend(iProc)))
    SurfSendBuf(iProc)%content=0.
  END IF
  IF(SurfExchange%nSidesRecv(iProc).GT.0) THEN
    RecvArraySize = SIZE(SurfRecvBuf(iProc)%content,DIM=1,KIND=4)
    SDEALLOCATE(SurfRecvBuf(iProc)%content)
    ALLOCATE(SurfRecvBuf(iProc)%content(RecvArraySize+SurfMesh%ReactiveSampSize*(nSurfSample**2)*SurfExchange%nSidesRecv(iProc)))
    SurfRecvBuf(iProc)%content=0.
  END IF
END DO ! iProc
#endif

END SUBROUTINE Init_SurfChemistrySampling


SUBROUTINE FinalizeSurfaceModel()
!===================================================================================================================================
!> Deallocate surface model vars
!===================================================================================================================================
! MODULES
USE MOD_SurfaceModel_Vars         ,ONLY: Adsorption, SurfDistInfo, SurfModel, SpecSurf
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: SurfModelAnalyzeInitIsDone
USE MOD_Particle_Vars             ,ONLY: PDM, PEM
USE MOD_Particle_Boundary_Vars    ,ONLY: nSurfSample, SurfMesh
#if USE_MPI
USE MOD_Particle_Boundary_Vars    ,ONLY: SurfCOMM
USE MOD_SurfaceModel_MPI_Vars     ,ONLY: SurfModelExchange
USE MOD_SurfaceModel_MPI_Vars     ,ONLY: SurfDistSendBuf,SurfDistRecvBuf
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iSubSurf,jSubSurf,iSurfSide,iCoord
#if USE_MPI
INTEGER                      :: iProc
#endif /*USE_MPI*/
!===================================================================================================================================
SurfModelAnalyzeInitIsDone=.FALSE.
! variables used if particles are kept after adsorption (currentyl not working)
SDEALLOCATE(PDM%ParticleAtWall)
SDEALLOCATE(PDM%PartAdsorbSideIndx)
SDEALLOCATE(PEM%wNumber)
! generaly used adsorption variables
#if (PP_TimeDiscMethod==42)
SDEALLOCATE(SurfModel%Info)
SDEALLOCATE(SurfModel%ProperInfo)
#endif
SDEALLOCATE(SpecSurf)
SDEALLOCATE(Adsorption%Coverage)
SDEALLOCATE(Adsorption%ProbAds)
SDEALLOCATE(Adsorption%ProbDes)
SDEALLOCATE(SurfModel%SumEvapPart)
SDEALLOCATE(SurfModel%SumDesorbPart)
SDEALLOCATE(SurfModel%SumReactPart)
SDEALLOCATE(SurfModel%SumAdsorbPart)
SDEALLOCATE(SurfModel%SumERDesorbed)
SDEALLOCATE(Adsorption%DensSurfAtoms)
SDEALLOCATE(Adsorption%AreaIncrease)
SDEALLOCATE(Adsorption%CrystalIndx)
SDEALLOCATE(Adsorption%ReactCoeff)
SDEALLOCATE(Adsorption%ReactEnergy)
SDEALLOCATE(Adsorption%ReactAccomodation)
SDEALLOCATE(Adsorption%SurfaceSpec)
! parameters for Kisliuk and Polanyi Wigner model (surfacemodel=1)
SDEALLOCATE(Adsorption%MaxCoverage)
SDEALLOCATE(Adsorption%InitStick)
SDEALLOCATE(Adsorption%PrefactorStick)
SDEALLOCATE(Adsorption%Adsorbexp)
SDEALLOCATE(Adsorption%Nu_a)
SDEALLOCATE(Adsorption%Nu_b)
SDEALLOCATE(Adsorption%DesorbEnergy)
SDEALLOCATE(Adsorption%Intensification)
! parameters for UBI-QEP model (surfacemodel=3)
SDEALLOCATE(Adsorption%HeatOfAdsZero)
SDEALLOCATE(Adsorption%DissocReact)
SDEALLOCATE(Adsorption%Diss_Powerfactor)
SDEALLOCATE(Adsorption%Diss_Prefactor)
SDEALLOCATE(Adsorption%ER_Powerfactor)
SDEALLOCATE(Adsorption%ER_Prefactor)
SDEALLOCATE(Adsorption%EDissBond)
SDEALLOCATE(Adsorption%EDissBondAdsorbPoly)
SDEALLOCATE(Adsorption%RecombReact)
SDEALLOCATE(Adsorption%ChemReactant)
SDEALLOCATE(Adsorption%ChemProduct)
SDEALLOCATE(Adsorption%Reactant_DissBond_K)
SDEALLOCATE(Adsorption%Product_DissBond_K)
SDEALLOCATE(Adsorption%Coordination)
SDEALLOCATE(Adsorption%DiCoord)
SDEALLOCATE(Adsorption%Ads_Powerfactor)
SDEALLOCATE(Adsorption%Ads_Prefactor)
SDEALLOCATE(Adsorption%TST_calc)
! surfaces distribution variables (currently surfacemodel=3)
IF (ALLOCATED(SurfDistInfo)) THEN
DO iSurfSide=1,SurfMesh%nMasterSides
  DO iSubSurf = 1,nSurfSample
  DO jSubSurf = 1,nSurfSample
#if USE_MPI
    SDEALLOCATE(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%Nbr_changed)
#endif /*USE_MPI*/
    SDEALLOCATE(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%nSites)
    SDEALLOCATE(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%SitesRemain)
    SDEALLOCATE(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%SurfAtomBondOrder)
    SDEALLOCATE(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%desorbnum_tmp)
    SDEALLOCATE(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%adsorbnum_tmp)
    SDEALLOCATE(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%reactnum_tmp)
    IF (ALLOCATED(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap)) THEN
      DO iCoord = 1,3
        SDEALLOCATE(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(iCoord)%UsedSiteMap)
        SDEALLOCATE(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(iCoord)%Species)
        SDEALLOCATE(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(iCoord)%BondAtomIndx)
        SDEALLOCATE(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(iCoord)%BondAtomIndy)
        SDEALLOCATE(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(iCoord)%NeighPos)
        SDEALLOCATE(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(iCoord)%NeighSite)
        SDEALLOCATE(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(iCoord)%IsNearestNeigh)
      END DO
      DEALLOCATE(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap)
    END IF
  END DO
  END DO
END DO
DEALLOCATE(SurfDistInfo)
END IF

#if USE_MPI
IF (ALLOCATED(SurfCOMM%MPINeighbor)) THEN
  DO iProc=1,SurfCOMM%nMPINeighbors
    SDEALLOCATE(SurfCOMM%MPINeighbor(iProc)%SurfDistSendList)
    SDEALLOCATE(SurfCOMM%MPINeighbor(iProc)%SurfDistRecvList)
  END DO
END IF
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
SDEALLOCATE(SurfModelExchange%nSurfDistSidesSend)
SDEALLOCATE(SurfModelExchange%nSurfDistSidesRecv)
SDEALLOCATE(SurfModelExchange%SurfDistSendRequest)
SDEALLOCATE(SurfModelExchange%SurfDistRecvRequest)
SDEALLOCATE(SurfModelExchange%NbrOfPos)

IF (ALLOCATED(SurfCOMM%MPINeighbor)) THEN
  DO iProc=1,SurfCOMM%nMPINeighbors
    SDEALLOCATE(SurfCOMM%MPINeighbor(iProc)%H2OSendList)
    SDEALLOCATE(SurfCOMM%MPINeighbor(iProc)%H2ORecvList)
    SDEALLOCATE(SurfCOMM%MPINeighbor(iProc)%O2HSendList)
    SDEALLOCATE(SurfCOMM%MPINeighbor(iProc)%O2HRecvList)
  END DO
END IF
IF (ALLOCATED(SurfModelExchange%H2OSendBuf)) THEN
  DO iProc=1,SurfCOMM%nMPINeighbors
    SDEALLOCATE(SurfModelExchange%H2OSendBuf(iProc)%content_int)
    SDEALLOCATE(SurfModelExchange%H2OSendBuf(iProc)%content)
  END DO
  DEALLOCATE(SurfModelExchange%H2OSendBuf)
END IF
IF (ALLOCATED(SurfModelExchange%H2ORecvBuf)) THEN
  DO iProc=1,SurfCOMM%nMPINeighbors
    SDEALLOCATE(SurfModelExchange%H2ORecvBuf(iProc)%content_int)
    SDEALLOCATE(SurfModelExchange%H2ORecvBuf(iProc)%content)
  END DO
  DEALLOCATE(SurfModelExchange%H2ORecvBuf)
END IF
IF (ALLOCATED(SurfModelExchange%O2HSendBuf)) THEN
  DO iProc=1,SurfCOMM%nMPINeighbors
    SDEALLOCATE(SurfModelExchange%O2HSendBuf(iProc)%content_int)
    SDEALLOCATE(SurfModelExchange%O2HSendBuf(iProc)%content)
  END DO
  DEALLOCATE(SurfModelExchange%O2HSendBuf)
END IF
IF (ALLOCATED(SurfModelExchange%O2HRecvBuf)) THEN
  DO iProc=1,SurfCOMM%nMPINeighbors
    SDEALLOCATE(SurfModelExchange%O2HRecvBuf(iProc)%content_int)
    SDEALLOCATE(SurfModelExchange%O2HRecvBuf(iProc)%content)
  END DO
  DEALLOCATE(SurfModelExchange%O2HRecvBuf)
END IF
SDEALLOCATE(SurfModelExchange%nH2OSidesSend)
SDEALLOCATE(SurfModelExchange%nH2OSidesRecv)
SDEALLOCATE(SurfModelExchange%nO2HSidesSend)
SDEALLOCATE(SurfModelExchange%nO2HSidesRecv)
SDEALLOCATE(SurfModelExchange%SendRequest)
SDEALLOCATE(SurfModelExchange%RecvRequest)
#endif /*USE_MPI*/

END SUBROUTINE FinalizeSurfaceModel

END MODULE MOD_SurfaceModel_Init
