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
PUBLIC :: InitLiquidSurfaceModel
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
#if (PP_TimeDiscMethod==42)
CALL prms%CreateLogicalOption(  'Surface-Adsorption-doTPD'&
  , 'Flag for TPD spectrum calculation. Only for TD=42 (RESERVOIR)','.FALSE.')
CALL prms%CreateRealOption(     'Surface-Adsorption-TPD-Beta'&
  , 'Temperature slope used for TPD calculation [K/s]. Surface temperature used for desorption probability is increased'//&
    'each iteration','0.')
#endif

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

CALL prms%CreateIntOption(      'Part-Species[$]-Recomb-PartnerSpec'&
                                          , 'TODO-DEFINE-PARAMETER\n'//&
                                            'Partner recombination species (nSpecies)','-1', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-Recomb-ResultSpec'&
                                          , 'TODO-DEFINE-PARAMETER\n'//&
                                            'Resulting recombination species (nSpecies)','-1', numberedmulti=.TRUE.)

CALL prms%CreateRealOption(     'Part-Species[$]-PartBound[$]-RecombinationCoeff'&
                                          , 'TODO-DEFINE-PARAMETER','0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-PartBound[$]-RecombinationEnergy'&
                                         , 'TODO-DEFINE-PARAMETER\n'//&
                                            'Energy transformed by reaction (nPartBound,nSpecies)','0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-PartBound[$]-RecombinationAccomodation'&
  , 'Define energy accomodation coefficient.\n'//&
    'Describes the percentage of reaction enthalpy of surface reaction transmitted to surface.','1.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-PartBound[$]-Coordination'&
  , 'Coordination at which particle of species [$] is bound on surface of boundary [$].\n'//&
    '1=hollow\n'//&
    '2=bridge\n'//&
    '3=on-top'//&
    '[surfacemodel=3]','0', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-PartBound[$]-DiCoordination'&
  , 'If particles of species [$] are di-, polyatomic and bind with additional coordination at Boundary [$].\n'//&
    '0: no DiCoordination\n'//&
    '1: bound via bridge bonding\n'//&
    '2: chelating binding\n'//&
    '[surfacemodel=3','0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-PartBound[$]-HeatOfAdsorption-K'&
  , 'Define heat of adsorption [K] on clear surface for binding atom of species [$] on boundary [$].\n'// &
    '[Assumption of on-top side bind, surfacemodel=3]','0.', numberedmulti=.TRUE.)

CALL prms%CreateStringOption(  'Particles-SurfCoverageFile'&
  , 'Give relative path to Surface-state-file ([Projectname]_DSMCSurfState***.h5) used for coverage init. '//&
    'File must be of same format as calculation (same mesh, same amount of species, cells and surfacemodel).\n'//&
    'If no file specified:\n'// &
    'Define Part-Species[$]-PartBound[$]-InitialCoverage for initial coverage different than 0.' &
  , 'none')
CALL prms%CreateRealOption(     'Part-Species[$]-PartBound[$]-InitialCoverage'&
  , 'Initial coverage used for species [$] and surfaces of boundary [$] in case of no surface-state-file init' &
  , '0.', numberedmulti=.TRUE.)


CALL prms%CreateRealOption(     'Particles-Surface-MacroParticleFactor'&
  , 'Weighting factor used for particles adsorbed on surface in case of reconstruction [surfacemodel=3].\n'//&
    'If one surface contains less then 10 surface atoms program abort is called.\n'//&
    'Default: Species(1)%MPF: Uses macro particle factor of species1.')
CALL prms%CreateIntOption(      'Surface-MaxDissNum'&
                                         , 'TODO-DEFINE-PARAMETER','0')
CALL prms%CreateIntOption(      'Surface-Nbr-DissocReactions'&
                                          , 'TODO-DEFINE-PARAMETER\n'//&
                                          'Resulting species for given dissoc (2,MaxDissNum,nSpecies)','0')
CALL prms%CreateIntOption(      'Surface-Nbr-ExchangeReactions'&
                                          , 'TODO-DEFINE-PARAMETER','0')


CALL prms%CreateRealOption(     'Part-Species[$]-Adsorption-Powerfactor'&
                                          , 'TODO-DEFINE-PARAMETER','0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Adsorption-Prefactor'&
                                          , 'TODO-DEFINE-PARAMETER','0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Adsorption-EDissBond'&
                                          , 'TODO-DEFINE-PARAMETER\n'//&
                                             'Bond dissociation energy (K) for diss into resulting'//&
                                          'species (ReactNum,nspecies)?','0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Adsorption-EDissBondPoly1'&
                                          , 'TODO-DEFINE-PARAMETER','0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Adsorption-EDissBondPoly2'&
                                         , 'TODO-DEFINE-PARAMETER','0.', numberedmulti=.TRUE.)

CALL prms%CreateIntArrayOption( 'Part-Species[$]-SurfDiss[$]-Products'&
                                          , 'TODO-DEFINE-PARAMETER','0 , 0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-SurfDiss[$]-Powerfactor'&
                                         , 'TODO-DEFINE-PARAMETER','0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-SurfDiss[$]-Prefactor'&
                                         , 'TODO-DEFINE-PARAMETER','0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-SurfDiss[$]-EDissBond'&
                                         , 'TODO-DEFINE-PARAMETER\n'//&
                                             'Bond dissociation energy (K) for diss into resulting'//&
                                          'species (ReactNum,nspecies)?','0.', numberedmulti=.TRUE.)

CALL prms%CreateRealOption(     'Part-Species[$]-Surf-ER[$]-Powerfactor'&
                                          , 'TODO-DEFINE-PARAMETER','0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surf-ER[$]-Prefactor'&
                                          , 'TODO-DEFINE-PARAMETER','0.', numberedmulti=.TRUE.)


CALL prms%CreateIntArrayOption( 'Surface-ExchReact[$]-Reactants'&
                                          , 'TODO-DEFINE-PARAMETER','0 , 0', numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'Surface-ExchReact[$]-Products'&
                                          , 'TODO-DEFINE-PARAMETER','0 , 0', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Surface-ExchReact[$]-DissBond_K-Reactants'&
                                          , 'TODO-DEFINE-PARAMETER','0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Surface-ExchReact[$]-DissBond_K-Products'&
                                         , 'TODO-DEFINE-PARAMETER','0. , 0.', numberedmulti=.TRUE.)

CALL prms%CreateIntOption(      'Surface-Adsorption-CalcTST'&
                                          , 'TODO-DEFINE-PARAMETER','0')
CALL prms%CreateRealOption(     'Surface-AdsorptionTST-PartitionMaxTemp'&
                                          , 'TODO-DEFINE-PARAMETER\n'//&
                                             'Temperature limit for pre-stored partition function (DEF: 20 000K)','10000.')
CALL prms%CreateRealOption(     'Surface-AdsorptionTST-PartitionInterval'&
                                          , 'TODO-DEFINE-PARAMETER\n'//&
                                             'Temperature interval for pre-stored partition function (DEF: 10K)','20.')

END SUBROUTINE DefineParametersSurfModel


SUBROUTINE InitSurfaceModel()
!===================================================================================================================================
!> Initialize surface model variables
!===================================================================================================================================
! MODULES
USE MOD_Globals                    ,ONLY: abort
USE MOD_Mesh_Vars                  ,ONLY: nElems, BC
USE MOD_DSMC_Vars                  ,ONLY: DSMC, CollisMode
USE MOD_Particle_Vars              ,ONLY: nSpecies, PDM, WriteMacroSurfaceValues, PartSurfaceModel
USE MOD_Particle_Vars              ,ONLY: KeepWallParticles, PEM
USE MOD_Particle_Mesh_Vars         ,ONLY: nTotalSides
USE MOD_ReadInTools                ,ONLY: GETINT,GETREAL,GETLOGICAL
USE MOD_Particle_Boundary_Vars     ,ONLY: nSurfSample, SurfMesh, nPartBound, PartBound
USE MOD_Particle_Boundary_Sampling ,ONLY: InitParticleBoundarySampling
USE MOD_SurfaceModel_Vars          ,ONLY: Adsorption
USE MOD_SurfaceModel_Tools         ,ONLY: CalcAdsorbProb, CalcDesorbProb
USE MOD_SMCR_Init                  ,ONLY: InitSMCR, InitSMCR_Chem
#if USE_MPI
USE MOD_SurfaceModel_MPI           ,ONLY: InitSurfModel_MPI, ExchangeCoverageInfo
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
IF (CollisMode.GT.1) THEN
  IF (PartSurfaceModel.EQ.1) CALL abort(&
__STAMP__&
,'Error in InitSurfaceModel: SurfaceModel 1 not working!')
  IF (PartSurfaceModel.GT.6 .OR. PartSurfaceModel.LT.0) CALL abort(&
__STAMP__&
,'Error in InitSurfaceModel: SurfaceModel must be 0,1,2 or 3!')
ELSE IF (CollisMode.LE.1) THEN
  CALL abort(&
__STAMP__&
,'Error in InitSurfaceModel - wrong collismode!')
END IF
! initialize variables only for processors that have any surfaces in own domain else they are skipped or not allocated
! initialize surface chemistry
!KeepWallParticles = GETLOGICAL('Particles-KeepWallParticles','.FALSE.')
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
#if (PP_TimeDiscMethod==42)
ALLOCATE( Adsorption%AdsorpInfo(1:nSpecies))
#endif
SELECT CASE(PartSurfaceModel)
CASE(1)
  ALLOCATE( Adsorption%MaxCoverage(1:SurfMesh%nSides,1:nSpecies),&
            Adsorption%InitStick(1:SurfMesh%nSides,1:nSpecies),&
            Adsorption%PrefactorStick(1:SurfMesh%nSides,1:nSpecies),&
            Adsorption%Adsorbexp(1:SurfMesh%nSides,1:nSpecies),&
            Adsorption%Nu_a(1:SurfMesh%nSides,1:nSpecies),&
            Adsorption%Nu_b(1:SurfMesh%nSides,1:nSpecies),&
            Adsorption%DesorbEnergy(1:SurfMesh%nSides,1:nSpecies),&
            Adsorption%Intensification(1:SurfMesh%nSides,1:nSpecies))
CASE(2)
  ALLOCATE( Adsorption%RecombCoeff(1:nPartBound,1:nSpecies),&
            Adsorption%RecombEnergy(1:nPartBound,1:nSpecies),&
            Adsorption%RecombAccomodation(1:nPartBound,1:nSpecies),&
            Adsorption%RecombData(1:2,1:nSpecies))
CASE(3)
  ALLOCATE( Adsorption%HeatOfAdsZero(1:nPartBound,1:nSpecies),&
            Adsorption%Coordination(1:nPartBound,1:nSpecies),&
            Adsorption%DiCoord(1:nPartBound,1:nSpecies))
END SELECT

! Initialize info and constants
DO iSpec = 1,nSpecies
#if (PP_TimeDiscMethod==42)
  Adsorption%AdsorpInfo(iSpec)%MeanProbAds  = 0.
  Adsorption%AdsorpInfo(iSpec)%MeanProbDes  = 0.
  Adsorption%AdsorpInfo(iSpec)%WallCollCount = 0
  Adsorption%AdsorpInfo(iSpec)%WallSpecNumCount = 0
  Adsorption%AdsorpInfo(iSpec)%NumOfAds = 0
  Adsorption%AdsorpInfo(iSpec)%NumOfDes = 0
  Adsorption%AdsorpInfo(iSpec)%Accomodation = 0
#endif
  WRITE(UNIT=hilf,FMT='(I0)') iSpec
  SELECT CASE(PartSurfaceModel)
  CASE(1)
    Adsorption%MaxCoverage(:,iSpec)     = GETREAL('Part-Species'//TRIM(hilf)//'-MaximumCoverage','0.')
    Adsorption%InitStick(:,iSpec)       = GETREAL('Part-Species'//TRIM(hilf)//'-InitialStick','0.')
    Adsorption%PrefactorStick(:,iSpec)  = GETREAL('Part-Species'//TRIM(hilf)//'-PrefactorStick','0.')
    Adsorption%Adsorbexp(:,iSpec)       = GETINT('Part-Species'//TRIM(hilf)//'-Adsorbexp','1')
    Adsorption%Nu_a(:,iSpec)            = GETREAL('Part-Species'//TRIM(hilf)//'-Nu-a','0.')
    Adsorption%Nu_b(:,iSpec)            = GETREAL('Part-Species'//TRIM(hilf)//'-Nu-b','0.')
    Adsorption%DesorbEnergy(:,iSpec)    = GETREAL('Part-Species'//TRIM(hilf)//'-Desorption-Energy-K','1.')
    Adsorption%Intensification(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Intensification-K','0.')
  CASE(2)
    Adsorption%RecombData(1,iSpec) = GETINT('Part-Species'//TRIM(hilf)//'-Recomb-PartnerSpec','-1')
    Adsorption%RecombData(2,iSpec) = GETINT('Part-Species'//TRIM(hilf)//'-Recomb-ResultSpec','-1')
    DO iPartBound=1,nPartBound
      IF((PartBound%TargetBoundCond(iPartBound).EQ.PartBound%ReflectiveBC).AND.PartBound%SolidState(iPartBound))THEN
        IF(PartBound%SolidReactive(iPartBound))THEN
          WRITE(UNIT=hilf2,FMT='(I0)') iPartBound
          hilf2=TRIM(hilf)//'-PartBound'//TRIM(hilf2)
          Adsorption%RecombCoeff(iPartBound,iSpec)        = GETREAL('Part-Species'//TRIM(hilf2)//'-RecombinationCoeff','0.')
          Adsorption%RecombEnergy(iPartBound,iSpec)       = GETREAL('Part-Species'//TRIM(hilf2)//'-RecombinationEnergy','0.')
          Adsorption%RecombAccomodation(iPartBound,iSpec) = GETREAL('Part-Species'//TRIM(hilf2)//'-RecombinationAccomodation','1.')
          IF ((Adsorption%RecombData(2,iSpec).EQ.-1).AND.(Adsorption%RecombCoeff(iPartBound,iSpec).NE.0.)) THEN
            CALL abort(&
                __STAMP__,&
                'Resulting species for species '//TRIM(hilf)//' not defined although recombination coefficient .GT. 0')
          END IF
        END IF
      END IF
    END DO
  CASE(3)
    DO iPartBound=1,nPartBound
      IF((PartBound%TargetBoundCond(iPartBound).EQ.PartBound%ReflectiveBC).AND.PartBound%SolidState(iPartBound))THEN
        IF(PartBound%SolidReactive(iPartBound))THEN
          WRITE(UNIT=hilf2,FMT='(I0)') iPartBound
          hilf2=TRIM(hilf)//'-PartBound'//TRIM(hilf2)
          Adsorption%Coordination(iPartBound,iSpec)  = GETINT('Part-Species'//TRIM(hilf2)//'-Coordination','0')
          Adsorption%DiCoord(iPartBound,iSpec)       = GETINT('Part-Species'//TRIM(hilf2)//'-DiCoordination','0')
          Adsorption%HeatOfAdsZero(iPartbound,iSpec) = GETREAL('Part-Species'//TRIM(hilf2)//'-HeatOfAdsorption-K','0.')
          IF (Adsorption%Coordination(iPartBound,iSpec).EQ.0)THEN
            WRITE(UNIT=hilf2,FMT='(I0)') iPartBound
            CALL abort(&
                __STAMP__,&
                'Coordination of Species '//TRIM(hilf)//' for catalytic particle boundary '//TRIM(hilf2)//' not defined')
          END IF
        END IF
      END IF
    END DO
  END SELECT
END DO
! Initialize temperature programmed desorption specific variables
#if (PP_TimeDiscMethod==42)
  Adsorption%TPD = GETLOGICAL('Surface-Adsorption-doTPD','.FALSE.')
  Adsorption%TPD_beta = GETREAL('Surface-Adsorption-TPD-Beta','0.')
  Adsorption%TPD_Temp = 0.
#endif
! Allocate and initialize adsorption variables
ALLOCATE( Adsorption%Coverage(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides,1:nSpecies),&
          Adsorption%ProbAds(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides,1:nSpecies),&
          Adsorption%ProbDes(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides,1:nSpecies),&
          Adsorption%SumDesorbPart(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides,1:nSpecies),&
          Adsorption%SumReactPart(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides,1:nSpecies),&
          Adsorption%SumAdsorbPart(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides,1:nSpecies),&
          Adsorption%SumERDesorbed(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides,1:nSpecies),&
          Adsorption%SurfSideToGlobSideMap(1:SurfMesh%nTotalSides),&
          Adsorption%DensSurfAtoms(1:SurfMesh%nTotalSides),&
          Adsorption%AreaIncrease(1:SurfMesh%nTotalSides),&
          Adsorption%CrystalIndx(1:SurfMesh%nTotalSides))

Adsorption%SurfSideToGlobSideMap(:) = -1
DO iSide = 1,nTotalSides
  IF (SurfMesh%SideIDToSurfID(iSide).LE.0) CYCLE
  Adsorption%SurfSideToGlobSideMap(SurfMesh%SideIDToSurfID(iSide)) = iSide
END DO
! Initialize surface properties from particle boundary values
DO iSide=1,SurfMesh%nTotalSides
  SideID = Adsorption%SurfSideToGlobSideMap(iSide)
  PartboundID = PartBound%MapToPartBC(BC(SideID))
  IF (PartBound%SolidReactive(PartboundID)) THEN
    !IF (PartSurfaceModel.EQ.3) Adsorption%SurfMassIC(iSide) = PartBound%SolidMassIC(PartBoundID)
    Adsorption%DensSurfAtoms(iSide) = PartBound%SolidPartDens(PartBoundID)
    Adsorption%AreaIncrease(iSide)  = PartBound%SolidAreaIncrease(PartBoundID)
    Adsorption%CrystalIndx(iSide)   = PartBound%SolidCrystalIndx(PartBoundID)
    Adsorption%DensSurfAtoms(iSide) = Adsorption%DensSurfAtoms(iSide)*Adsorption%AreaIncrease(iSide)
  ELSE
    Adsorption%DensSurfAtoms(iSide) = 0
    Adsorption%AreaIncrease(iSide)  = 0
    Adsorption%CrystalIndx(iSide)   = 0
    Adsorption%DensSurfAtoms(iSide) = 0
  END IF
END DO

! Initialize surface coverage
CALL InitSurfCoverage()

IF (SurfMesh%SurfOnProc) THEN
  Adsorption%ProbAds(:,:,:,:)       = 0.
  Adsorption%ProbDes(:,:,:,:)       = 0.
  Adsorption%SumDesorbPart(:,:,:,:) = 0
  Adsorption%SumAdsorbPart(:,:,:,:) = 0
  Adsorption%SumReactPart(:,:,:,:)  = 0
  Adsorption%SumERDesorbed(:,:,:,:) = 0
#if USE_MPI
  CALL InitSurfModel_MPI()
#endif /*USE_MPI*/
END IF ! SurfMesh%SurfOnProc

SELECT CASE(PartSurfaceModel)
CASE(1)
  ! Define number of possible recombination reactions per species needed for sampling
  Adsorption%RecombNum = 0
CASE(2)
  IF (SurfMesh%SurfOnProc) THEN
    CALL CalcDesorbProb()
    CALL CalcAdsorbProb()
#if USE_MPI
    CALL ExchangeCoverageInfo()
#endif /*USE_MPI*/
  END IF
  ! Define number of possible recombination reactions per species needed for sampling
  Adsorption%RecombNum = 1
CASE(3)
  CALL InitSMCR()
  CALL InitSMCR_Chem()
END SELECT

IF (WriteMacroSurfaceValues.OR.DSMC%CalcSurfaceVal) THEN
  IF(SurfMesh%SurfOnProc) CALL Init_SurfChemistrySampling()
END IF

END SUBROUTINE InitSurfaceModel


SUBROUTINE InitLiquidSurfaceModel()
!===================================================================================================================================
!> init of particle liquid interfaces
!> 1) mark sides for liquid boundary consideration
!> 2) init MPI communication
!===================================================================================================================================
USE MOD_Particle_Vars          ,ONLY: nSpecies, WriteMacroSurfaceValues
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_SurfaceModel_Vars      ,ONLY: Liquid
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, SampWall
#if USE_MPI
USE MOD_Particle_Boundary_Vars ,ONLY: SurfCOMM
USE MOD_Particle_MPI_Vars      ,ONLY: SurfSendBuf, SurfRecvBuf, SurfExchange
#endif /*USE_MPI*/
!===================================================================================================================================
IMPLICIT NONE
!===================================================================================================================================
! INPUT VARIABLES
!===================================================================================================================================
! OUTPUT VARIABLES
!===================================================================================================================================
! LOCAL VARIABLES
INTEGER                          :: iSpec, iSide
#if USE_MPI
INTEGER                          :: iProc, SendArraySize, RecvArraySize
#endif
!===================================================================================================================================
! allocate info and constants
#if (PP_TimeDiscMethod==42)
ALLOCATE( Liquid%Info(1:nSpecies))
#endif
! initialize info and constants
DO iSpec = 1,nSpecies
#if (PP_TimeDiscMethod==42)
  Liquid%Info(iSpec)%MeanProbAds = 0.
  Liquid%Info(iSpec)%MeanProbDes = 0.
  Liquid%Info(iSpec)%WallCollCount = 0
  Liquid%Info(iSpec)%NumOfAds = 0
  Liquid%Info(iSpec)%NumOfDes = 0
#endif
END DO
IF (.NOT.SurfMesh%SurfOnProc) RETURN
! allocate and initialize liquid variables
ALLOCATE( Liquid%ProbCondens(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides,1:nSpecies),&
          Liquid%ProbEvap(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides,1:nSpecies),&
          Liquid%SumCondensPart(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides,1:nSpecies),&
          Liquid%SumEvapPart(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides,1:nSpecies) )
Liquid%ProbCondens(:,:,:,:) = 1. !0.
Liquid%ProbEvap(:,:,:,:) = 0.
Liquid%SumCondensPart(:,:,:,:) = 0
Liquid%SumEvapPart(:,:,:,:) = 0

IF (WriteMacroSurfaceValues.OR.DSMC%CalcSurfaceVal) THEN
  DO iSide=1,SurfMesh%nTotalSides ! caution: iSurfSideID
    ALLOCATE(SampWall(iSide)%Evaporation(1:nSpecies+1,1:nSurfSample,1:nSurfSample))
    SampWall(iSide)%Evaporation=0.
  END DO

#if USE_MPI
  ! Reallocate buffer for mpi communication of sampling
  DO iProc=1,SurfCOMM%nMPINeighbors
    SendArraySize = SIZE(SurfSendBuf(iProc)%content,DIM=1,KIND=4)
    RecvArraySize = SIZE(SurfRecvBuf(iProc)%content,DIM=1,KIND=4)
    SDEALLOCATE(SurfSendBuf(iProc)%content)
    SDEALLOCATE(SurfRecvBuf(iProc)%content)
    IF(SurfExchange%nSidesSend(iProc).GT.0) THEN
      ALLOCATE(SurfSendBuf(iProc)%content(SendArraySize+(nSpecies+1)&
                                          *(nSurfSample**2)*SurfExchange%nSidesSend(iProc)))
      SurfSendBuf(iProc)%content=0.
    END IF
    IF(SurfExchange%nSidesRecv(iProc).GT.0) THEN
      ALLOCATE(SurfRecvBuf(iProc)%content(RecvArraySize+(nSpecies+1)&
                                          *(nSurfSample**2)*SurfExchange%nSidesRecv(iProc)))
      SurfRecvBuf(iProc)%content=0.
    END IF
  END DO ! iProc
#endif
END IF

END SUBROUTINE InitLiquidSurfaceModel


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
USE MOD_Particle_Vars          ,ONLY: nSpecies, PartSurfaceModel
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
LOGICAL                          :: SurfCalcDataExists, WallmodelExists
INTEGER                          :: WallModel_HDF5
REAL                             :: Version_HDF5
!===================================================================================================================================

IF (DoRestart) THEN
  ! check if datasets for restarting from state exists in state file used for restart
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=PartMPI%COMM)!MPI_COMM_WORLD)
  WallmodelExists=.FALSE.
  CALL DatasetExists(File_ID,'WallModel',WallmodelExists,attrib=.TRUE.)
  IF (WallmodelExists) THEN
    CALL ReadAttribute(File_ID,'WallModel',1,IntegerScalar=WallModel_HDF5)
    IF (WallModel_HDF5.NE.PartSurfaceModel) WallmodelExists=.FALSE.
  END IF
  SurfCalcDataExists=.FALSE.
  CALL DatasetExists(File_ID,'SurfCalcData',SurfCalcDataExists)
  CALL CloseDataFile()
  IF (SurfCalcDataExists.AND.WallmodelExists) THEN
    SWRITE(UNIT_StdOut, '(A)')' Surface will be initialized from Restartfile'
    IF (SurfMesh%SurfOnProc) THEN
      ! write zeros into global coverage array for each surface
      Adsorption%Coverage(:,:,:,:) = 0.
    END IF
    ! return as coverage is set in restart.f90
    RETURN
  END IF
END IF ! DoRestart

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
    IF (PartSurfaceModel.EQ.2) THEN
      Adsorption%Coverage(:,:,:,:) = 0.
    ELSE
      DO iSurfSide=1,SurfMesh%nSides
        SideID = Adsorption%SurfSideToGlobSideMap(iSurfSide)
        PartboundID = PartBound%MapToPartBC(BC(SideID))
        DO iSpec=1,nSpecies
          Adsorption%Coverage(:,:,iSurfSide,iSpec) = Coverage_tmp(PartBoundID,iSpec)
        END DO
      END DO
    END IF
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
      nVar2D = 6
    ELSE
      nVar2D = 5
    END IF
    nVar2D_Spec = 4
    nVar2D_Total = nVar2D + nVar2D_Spec*nSpecies
    IF (nVarSurf_HDF5.NE.nVar2D_Total) CALL abort(&
__STAMP__&
,'Error in surface coverage init: number of variables in HDF5-file does not match!')
    SDEALLOCATE(SurfState_HDF5)
    ALLOCATE(SurfState_HDF5(1:nVarSurf_HDF5,1:nSurfSample,1:nSurfSample,SurfMesh%nSides))

    ! Associate construct for integer KIND=8 possibility
    ASSOCIATE (&
          nVarSurf_HDF5   => INT(nVarSurf_HDF5,IK)   ,&
          nSurfSample     => INT(nSurfSample,IK)     ,&
          nSides          => INT(SurfMesh%nSides,IK) ,&
          offsetSurfSide  => INT(offsetSurfSide,IK)  )
      CALL ReadArray(  'SurfaceData',4,(/nVarSurf_HDF5,nSurfSample,nSurfSample,nSides/)&
                     ,offsetSurfSide,4,RealArray=SurfState_HDF5)
    END ASSOCIATE
    iVar = 3
    DO iSpec = 1, nSpecies
      DO iSurfSide = 1, SurfMesh%nSides
        SideID = Adsorption%SurfSideToGlobSideMap(iSurfSide)
        PartboundID = PartBound%MapToPartBC(BC(SideID))
        IF (PartBound%SolidReactive(PartboundID)) THEN
          DO jSubSurf = 1, nSurfSample
            DO iSubSurf = 1, nSurfSample
              Adsorption%Coverage(iSubSurf,jSubSurf,iSurfSide,iSpec) = SurfState_HDF5(iVar,iSubSurf,jSubSurf,iSurfSide)
            END DO ! iSubSurf = 1, nSurfSample
          END DO ! jSubSurf = 1, nSurfSample
        ELSE
          Adsorption%Coverage(:,:,iSurfSide,iSpec) = 0.
        END IF
      END DO ! iSurfSide = 1, SurfMesh%nSides
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
!> Initializion of additional surface sampling if PartWallModel>0
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

DO iSide=1,SurfMesh%nTotalSides ! caution: iSurfSideID
  ALLOCATE(SampWall(iSide)%Adsorption(1:(nSpecies+1),1:nSurfSample,1:nSurfSample))
  SampWall(iSide)%Adsorption=0.
  ALLOCATE(SampWall(iSide)%Accomodation(1:nSpecies,1:nSurfSample,1:nSurfSample))
  SampWall(iSide)%Accomodation=0.
  ALLOCATE(SampWall(iSide)%Reaction(1:Adsorption%RecombNum,1:nSpecies,1:nSurfSample,1:nSurfSample))
  SampWall(iSide)%Reaction=0.
END DO
#if USE_MPI
! Reallocate buffer for mpi communication of sampling
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesSend(iProc).GT.0) THEN
    SendArraySize = SIZE(SurfSendBuf(iProc)%content,DIM=1,KIND=4)
    SDEALLOCATE(SurfSendBuf(iProc)%content)
    ALLOCATE(SurfSendBuf(iProc)%content(SendArraySize+(2*nSpecies+1+(Adsorption%RecombNum*nSpecies))&
                                        *(nSurfSample**2)*SurfExchange%nSidesSend(iProc)))
    !ALLOCATE(SurfSendBuf(iProc)%content((2*nSpecies+1+SurfMesh%SampSize+(Adsorption%RecombNum*nSpecies))&
    !                                    *(nSurfSample**2)*SurfExchange%nSidesSend(iProc)))
    SurfSendBuf(iProc)%content=0.
  END IF
  IF(SurfExchange%nSidesRecv(iProc).GT.0) THEN
    RecvArraySize = SIZE(SurfRecvBuf(iProc)%content,DIM=1,KIND=4)
    SDEALLOCATE(SurfRecvBuf(iProc)%content)
    ALLOCATE(SurfRecvBuf(iProc)%content(RecvArraySize+(2*nSpecies+1+(Adsorption%RecombNum*nSpecies))&
                                        *(nSurfSample**2)*SurfExchange%nSidesRecv(iProc)))
    !ALLOCATE(SurfRecvBuf(iProc)%content((2*nSpecies+1+SurfMesh%SampSize+(Adsorption%RecombNum*nSpecies))&
    !                                    *(nSurfSample**2)*SurfExchange%nSidesRecv(iProc)))
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
USE MOD_SurfaceModel_Vars         ,ONLY: Adsorption, SurfDistInfo, Liquid
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: SurfModelAnalyzeInitIsDone
USE MOD_Particle_Vars             ,ONLY: PDM, PEM
USE MOD_Particle_Boundary_Vars    ,ONLY: nSurfSample, SurfMesh
#if USE_MPI
USE MOD_Particle_Boundary_Vars    ,ONLY: SurfCOMM
USE MOD_SurfaceModel_MPI_Vars     ,ONLY: SurfModelExchange
USE MOD_SurfaceModel_MPI_Vars     ,ONLY: AdsorbSendBuf,AdsorbRecvBuf,SurfDistSendBuf,SurfDistRecvBuf
USE MOD_SurfaceModel_MPI_Vars     ,ONLY: SurfCoverageSendBuf,SurfCoverageRecvBuf
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
SDEALLOCATE(Adsorption%AdsorpInfo)
SDEALLOCATE(Adsorption%AdsorpReactInfo)
#endif
SDEALLOCATE(Adsorption%Coverage)
SDEALLOCATE(Adsorption%ProbAds)
SDEALLOCATE(Adsorption%ProbDes)
SDEALLOCATE(Adsorption%SumDesorbPart)
SDEALLOCATE(Adsorption%SumReactPart)
SDEALLOCATE(Adsorption%SumAdsorbPart)
SDEALLOCATE(Adsorption%SumERDesorbed)
SDEALLOCATE(Adsorption%DensSurfAtoms)
SDEALLOCATE(Adsorption%AreaIncrease)
SDEALLOCATE(Adsorption%CrystalIndx)
SDEALLOCATE(Adsorption%SurfSideToGlobSideMap)
! parameters for Kisliuk and Polanyi Wigner model (surfacemodel=1)
SDEALLOCATE(Adsorption%MaxCoverage)
SDEALLOCATE(Adsorption%InitStick)
SDEALLOCATE(Adsorption%PrefactorStick)
SDEALLOCATE(Adsorption%Adsorbexp)
SDEALLOCATE(Adsorption%Nu_a)
SDEALLOCATE(Adsorption%Nu_b)
SDEALLOCATE(Adsorption%DesorbEnergy)
SDEALLOCATE(Adsorption%Intensification)
! parameters for surface recombination model (surfacemodel=2)
SDEALLOCATE(Adsorption%RecombCoeff)
SDEALLOCATE(Adsorption%RecombAccomodation)
SDEALLOCATE(Adsorption%RecombEnergy)
SDEALLOCATE(Adsorption%RecombData)
! parameters for UBI-QEP model (surfacemodel=3)
SDEALLOCATE(Adsorption%HeatOfAdsZero)
SDEALLOCATE(Adsorption%DissocReact)
SDEALLOCATE(Adsorption%Diss_Powerfactor)
SDEALLOCATE(Adsorption%Diss_Prefactor)
SDEALLOCATE(Adsorption%ER_Powerfactor)
SDEALLOCATE(Adsorption%ER_Prefactor)
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
SDEALLOCATE(Adsorption%TST_calc)
! surfaces distribution variables (currently surfacemodel=3)
IF (ALLOCATED(SurfDistInfo)) THEN
DO iSurfSide=1,SurfMesh%nSides
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
  SDEALLOCATE(SurfCOMM%MPINeighbor(iProc)%CoverageSendList)
  SDEALLOCATE(SurfCOMM%MPINeighbor(iProc)%CoverageRecvList)
  END DO
END IF
SDEALLOCATE(SurfModelExchange%nSurfDistSidesSend)
SDEALLOCATE(SurfModelExchange%nSurfDistSidesRecv)
SDEALLOCATE(SurfModelExchange%SurfDistSendRequest)
SDEALLOCATE(SurfModelExchange%SurfDistRecvRequest)
SDEALLOCATE(SurfModelExchange%NbrOfPos)
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
SDEALLOCATE(SurfModelExchange%nCoverageSidesSend)
SDEALLOCATE(SurfModelExchange%nCoverageSidesRecv)
IF (ALLOCATED(SurfCoverageSendBuf)) THEN
  DO iProc=1,SurfCOMM%nMPINeighbors
    SDEALLOCATE(SurfCoverageSendBuf(iProc)%content)
  END DO
  DEALLOCATE(SurfCoverageSendBuf)
END IF
IF (ALLOCATED(SurfCoverageRecvBuf)) THEN
  DO iProc=1,SurfCOMM%nMPINeighbors
    SDEALLOCATE(SurfCoverageRecvBuf(iProc)%content)
  END DO
  DEALLOCATE(SurfCoverageRecvBuf)
END IF
SDEALLOCATE(SurfModelExchange%nSidesSend)
SDEALLOCATE(SurfModelExchange%nSidesRecv)
SDEALLOCATE(SurfModelExchange%SendRequest)
SDEALLOCATE(SurfModelExchange%RecvRequest)
#endif /*USE_MPI*/

#if (PP_TimeDiscMethod==42)
SDEALLOCATE(Liquid%Info)
#endif

SDEALLOCATE(Liquid%ProbCondens)
SDEALLOCATE(Liquid%ProbEvap)
SDEALLOCATE(Liquid%SumCondensPart)
SDEALLOCATE(Liquid%SumEvapPart)

END SUBROUTINE FinalizeSurfaceModel

END MODULE MOD_SurfaceModel_Init
