#include "boltzplatz.h"

MODULE MOD_DSMC_SurfModelInit
!===================================================================================================================================
!> Module for initialization of surface models
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
PUBLIC::DefineParametersSurfModel

CONTAINS

!==================================================================================================================================
!> Define parameters for surface model in DSMC
!==================================================================================================================================
SUBROUTINE DefineParametersSurfModel()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("DSMC Surfmodel")

CALL prms%CreateLogicalOption(  'Particles-KeepWallParticles'&
  , 'Flag to track adsorbed Particles on surface if they are adsorbed. Currently only [FALSE] implemented','.FALSE.')
CALL prms%CreateLogicalOption(  'Particles--DSMC-Adsorption-doTPD'&
  , 'Flag for TPD spectrum calculation. Only for TD=42 (RESERVOIR)','.FALSE.')
CALL prms%CreateRealOption(     'Particles-DSMC-Adsorption-TPD-Beta'&
  , 'Temperature slope used for TPD calculation [K/s]. Surface temperature used for desorption probability is increased'//&
    'each iteration','0.')

CALL prms%CreateRealOption(     'Part-Species[$]-MaximumCoverage'&
  , 'Maximum coverage of surfaces with species [$] (used for wallmodel=1)','0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-InitialStick'&
  , 'Initial sticking coefficient (S_0) of species [$] for surfaces (used for Kisliuk model, wallmodel=1)' &
  , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-PrefactorStick'&
  , 'Prefactor of sticking coefficient of species [$] for surfaces (used for Kisliuk model, wallmodel=1)' &
  , '0.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-Adsorbexp'&
  , 'Adsorption exponent of species [$] for surfaces (used for Kisliuk model, wallmodel=1)' &
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
    '[wallmodel=3]','0', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-PartBound[$]-DiCoordination'&
  , 'If particles of species [$] are di-, polyatomic and bind with additional coordination at Boundary [$].\n'//&
    '0: no DiCoordination\n'//&
    '1: bound via bridge bonding\n'//&
    '2: chelating binding\n'//&
    '[wallmodel=3','0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-PartBound[$]-HeatOfAdsorption-K'&
  , 'Define heat of adsorption [K] on clear surface for binding atom of species [$] on boundary [$].\n'// &
    '[Assumption of on-top side bind, wallmodel=3]','0.', numberedmulti=.TRUE.)

CALL prms%CreateStringOption(  'Particles-SurfCoverageFile'&
  , 'Give relative path to Surface-state-file ([Projectname]_DSMCSurfState***.h5) used for coverage init. '//&
    'File must be of same format as calculation (same mesh, same amount of species, cells and wallmodel).\n'//&
    'If no file specified:\n'// &
    'Define Part-Species[$]-PartBound[$]-InitialCoverage for initial coverage different than 0.' &
  , 'none')
CALL prms%CreateRealOption(     'Part-Species[$]-PartBound[$]-InitialCoverage'&
  , 'Initial coverage used for species [$] and surfaces of boundary [$] in case of no surface-state-file init' &
  , '0.', numberedmulti=.TRUE.)


CALL prms%CreateRealOption(     'Particles-Surface-MacroParticleFactor'&
  , 'Weighting factor used for particles adsorbed on surface in case of reconstruction [wallmodel=3].\n'//&
    'If one surface contains less then 10 surface atoms program abort is called.\n'//&
    'Default: Species(1)%MPF: Uses macro particle factor of species1.')
CALL prms%CreateIntOption(      'Part-Species-MaxDissNum'&
                                         , 'TODO-DEFINE-PARAMETER','0')
CALL prms%CreateIntOption(      'Part-SurfChem-Nbr-DissocReactions'&
                                          , 'TODO-DEFINE-PARAMETER\n'//&
                                          'Resulting species for given dissoc (2,MaxDissNum,nSpecies)','0')
CALL prms%CreateIntOption(      'Part-SurfChem-Nbr-ExchangeReactions'&
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


CALL prms%CreateIntArrayOption( 'Part-SurfChem-Disprop[$]-Reactants'&
                                          , 'TODO-DEFINE-PARAMETER','0 , 0', numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'Part-SurfChem-Disprop[$]-Products'&
                                          , 'TODO-DEFINE-PARAMETER','0 , 0', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-SurfChem-Disprop[$]-DissBond_K-Reactants'&
                                          , 'TODO-DEFINE-PARAMETER','0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-SurfChem-Disprop[$]-DissBond_K-Products'&
                                         , 'TODO-DEFINE-PARAMETER','0. , 0.', numberedmulti=.TRUE.)

CALL prms%CreateIntOption(      'Particles-DSMC-Adsorption-CalcTST'&
                                          , 'TODO-DEFINE-PARAMETER','0')
CALL prms%CreateRealOption(     'Particles-DSMC-AdsorptionTST-PartitionMaxTemp'&
                                          , 'TODO-DEFINE-PARAMETER\n'//&
                                             'Temperature limit for pre-stored partition function (DEF: 20 000K)','10000.')
CALL prms%CreateRealOption(     'Particles-DSMC-AdsorptionTST-PartitionInterval'&
                                          , 'TODO-DEFINE-PARAMETER\n'//&
                                             'Temperature interval for pre-stored partition function (DEF: 10K)','20.')

END SUBROUTINE DefineParametersSurfModel


SUBROUTINE InitDSMCSurfModel()
!===================================================================================================================================
!> Init of surface model variables
!===================================================================================================================================
! MODULES
USE MOD_Globals                ,ONLY: abort
USE MOD_Mesh_Vars              ,ONLY: nElems, BC
USE MOD_DSMC_Vars              ,ONLY: Adsorption, DSMC!, CollisMode
USE MOD_PARTICLE_Vars          ,ONLY: nSpecies, PDM
USE MOD_PARTICLE_Vars          ,ONLY: KeepWallParticles, PEM
USE MOD_Particle_Mesh_Vars     ,ONLY: nTotalSides
USE MOD_ReadInTools
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, nPartBound, PartBound
USE MOD_DSMC_SurfModel_Tools   ,ONLY: CalcAdsorbProb, CalcDesorbProb
#if (PP_TimeDiscMethod==42)
USE MOD_DSMC_SurfModel_Tools   ,ONLY: CalcDiffusion
#endif
#ifdef MPI
USE MOD_DSMC_SurfModel_Tools   ,ONLY: ExchangeCoverageInfo
USE MOD_Particle_Boundary_Vars ,ONLY: SurfCOMM
USE MOD_Particle_MPI_Vars      ,ONLY: AdsorbSendBuf,AdsorbRecvBuf,SurfExchange,SurfCoverageSendBuf,SurfCoverageRecvBuf
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
CHARACTER(LEN=255)               :: SurfaceFileName
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
! initialize variables only for processors that have any surfaces in own domain else they are skipped or not allocated
! initialize surface chemistry
KeepWallParticles = GETLOGICAL('Particles-KeepWallParticles','.FALSE.')
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
IF (DSMC%WallModel.EQ.1) THEN
  ALLOCATE( Adsorption%MaxCoverage(1:SurfMesh%nSides,1:nSpecies),&
            Adsorption%InitStick(1:SurfMesh%nSides,1:nSpecies),&
            Adsorption%PrefactorStick(1:SurfMesh%nSides,1:nSpecies),&
            Adsorption%Adsorbexp(1:SurfMesh%nSides,1:nSpecies),&
            Adsorption%Nu_a(1:SurfMesh%nSides,1:nSpecies),&
            Adsorption%Nu_b(1:SurfMesh%nSides,1:nSpecies),&
            Adsorption%DesorbEnergy(1:SurfMesh%nSides,1:nSpecies),&
            Adsorption%Intensification(1:SurfMesh%nSides,1:nSpecies))
ELSE IF (DSMC%WallModel.EQ.2) THEN
  ALLOCATE( Adsorption%RecombCoeff(1:nPartBound,1:nSpecies),&
            Adsorption%RecombEnergy(1:nPartBound,1:nSpecies),&
            Adsorption%RecombAccomodation(1:nPartBound,1:nSpecies),&
            Adsorption%RecombData(1:2,1:nSpecies))
ELSE IF (DSMC%WallModel.EQ.3) THEN
  ALLOCATE( Adsorption%HeatOfAdsZero(1:nPartBound,1:nSpecies),&
            Adsorption%Coordination(1:nPartBound,1:nSpecies),&
            Adsorption%DiCoord(1:nPartBound,1:nSpecies))
END IF

! initialize info and constants
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
  IF (DSMC%WallModel.EQ.1) THEN
    Adsorption%MaxCoverage(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-MaximumCoverage','0.')
    Adsorption%InitStick(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-InitialStick','0.')
    Adsorption%PrefactorStick(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-PrefactorStick','0.')
    Adsorption%Adsorbexp(:,iSpec) = GETINT('Part-Species'//TRIM(hilf)//'-Adsorbexp','1')
    Adsorption%Nu_a(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Nu-a','0.')
    Adsorption%Nu_b(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Nu-b','0.')
    Adsorption%DesorbEnergy(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Desorption-Energy-K','1.')
    Adsorption%Intensification(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Intensification-K','0.')
  ELSE IF (DSMC%WallModel.EQ.2) THEN
    Adsorption%RecombData(1,iSpec) = GETINT('Part-Species'//TRIM(hilf)//'-Recomb-PartnerSpec','-1')
    Adsorption%RecombData(2,iSpec) = GETINT('Part-Species'//TRIM(hilf)//'-Recomb-ResultSpec','-1')
    DO iPartBound=1,nPartBound
      IF((PartBound%TargetBoundCond(iPartBound).EQ.PartBound%ReflectiveBC).AND.PartBound%SolidState(iPartBound))THEN
        IF(PartBound%SolidCatalytic(iPartBound))THEN
          WRITE(UNIT=hilf2,FMT='(I0)') iPartBound
          hilf2=TRIM(hilf)//'-PartBound'//TRIM(hilf2)
          Adsorption%RecombCoeff(iPartBound,iSpec) = GETREAL('Part-Species'//TRIM(hilf2)//'-RecombinationCoeff','0.')
          Adsorption%RecombEnergy(iPartBound,iSpec) = GETREAL('Part-Species'//TRIM(hilf2)//'-RecombinationEnergy','0.')
          Adsorption%RecombAccomodation(iPartBound,iSpec) = GETREAL('Part-Species'//TRIM(hilf2)//'-RecombinationAccomodation','1.')
          IF ((Adsorption%RecombData(2,iSpec).EQ.-1).AND.(Adsorption%RecombCoeff(iPartBound,iSpec).NE.0.)) THEN
            CALL abort(&
__STAMP__,&
'Resulting species for species '//TRIM(hilf)//' not defined although recombination coefficient .GT. 0')
          END IF
        END IF
      END IF
    END DO
  ELSE IF (DSMC%WallModel.EQ.3) THEN
    DO iPartBound=1,nPartBound
      IF((PartBound%TargetBoundCond(iPartBound).EQ.PartBound%ReflectiveBC).AND.PartBound%SolidState(iPartBound))THEN
        IF(PartBound%SolidCatalytic(iPartBound))THEN
          WRITE(UNIT=hilf2,FMT='(I0)') iPartBound
          hilf2=TRIM(hilf)//'-PartBound'//TRIM(hilf2)
          Adsorption%Coordination(iPartBound,iSpec) = GETINT('Part-Species'//TRIM(hilf2)//'-Coordination','0')
          Adsorption%DiCoord(iPartBound,iSpec) = GETINT('Part-Species'//TRIM(hilf2)//'-DiCoordination','0')
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
  END IF
END DO
! initialize temperature programmed desorption specific variables
#if (PP_TimeDiscMethod==42)
  Adsorption%TPD = GETLOGICAL('Particles-DSMC-Adsorption-doTPD','.FALSE.')
  Adsorption%TPD_beta = GETREAL('Particles-DSMC-Adsorption-TPD-Beta','0.')
  Adsorption%TPD_Temp = 0.
#endif
! allocate and initialize adsorption variables
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
! initialize surface properties from particle boundary values
DO iSide=1,SurfMesh%nTotalSides
  SideID = Adsorption%SurfSideToGlobSideMap(iSide)
  PartboundID = PartBound%MapToPartBC(BC(SideID))
  IF (PartBound%SolidCatalytic(PartboundID)) THEN
    !IF (DSMC%WallModel.EQ.3) Adsorption%SurfMassIC(iSide) = PartBound%SolidMassIC(PartBoundID)
    Adsorption%DensSurfAtoms(iSide) = PartBound%SolidPartDens(PartBoundID)
    Adsorption%AreaIncrease(iSide) = PartBound%SolidAreaIncrease(PartBoundID)
    Adsorption%CrystalIndx(iSide) = PartBound%SolidCrystalIndx(PartBoundID)
    Adsorption%DensSurfAtoms(iSide) = Adsorption%DensSurfAtoms(iSide)*Adsorption%AreaIncrease(iSide)
  ELSE
    Adsorption%DensSurfAtoms(iSide) = 0
    Adsorption%AreaIncrease(iSide) = 0
    Adsorption%CrystalIndx(iSide) = 0
    Adsorption%DensSurfAtoms(iSide) = 0
  END IF
END DO

! initialize surface coverage
CALL InitSurfCoverage()

IF (SurfMesh%SurfOnProc) THEN
  Adsorption%ProbAds(:,:,:,:) = 0.
  Adsorption%ProbDes(:,:,:,:) = 0.
  Adsorption%SumDesorbPart(:,:,:,:) = 0
  Adsorption%SumAdsorbPart(:,:,:,:) = 0
  Adsorption%SumReactPart(:,:,:,:) = 0
  Adsorption%SumERDesorbed(:,:,:,:) = 0

#ifdef MPI
  ! allocate send and receive buffer
  ALLOCATE(AdsorbSendBuf(SurfCOMM%nMPINeighbors))
  ALLOCATE(AdsorbRecvBuf(SurfCOMM%nMPINeighbors))
  DO iProc=1,SurfCOMM%nMPINeighbors
    ALLOCATE(AdsorbSendBuf(iProc)%content_int(2*nSpecies*(nSurfSample**2)*SurfExchange%nSidesSend(iProc)))
    ALLOCATE(AdsorbRecvBuf(iProc)%content_int(2*nSpecies*(nSurfSample**2)*SurfExchange%nSidesRecv(iProc)))
    AdsorbSendBuf(iProc)%content_int=0
    AdsorbRecvBuf(iProc)%content_int=0
  END DO ! iProc
  IF (DSMC%WallModel.EQ.2) THEN
    ALLOCATE(SurfCoverageSendBuf(SurfCOMM%nMPINeighbors))
    ALLOCATE(SurfCoverageRecvBuf(SurfCOMM%nMPINeighbors))
    ALLOCATE(SurfExchange%nCoverageSidesSend(SurfCOMM%nMPINeighbors))
    ALLOCATE(SurfExchange%nCoverageSidesRecv(SurfCOMM%nMPINeighbors))
    DO iProc=1,SurfCOMM%nMPINeighbors
      SurfExchange%nCoverageSidesSend(iProc) = SurfExchange%nSidesRecv(iProc)
      SurfExchange%nCoverageSidesRecv(iProc) = SurfExchange%nSidesSend(iProc)
      ALLOCATE(SurfCoverageSendBuf(iProc)%content(1:3*nSpecies*(nSurfSample**2)*SurfExchange%nCoverageSidesSend(iProc)))
      ALLOCATE(SurfCoverageRecvBuf(iProc)%content(1:3*nSpecies*(nSurfSample**2)*SurfExchange%nCoverageSidesRecv(iProc)))
      ALLOCATE(SurfCOMM%MPINeighbor(iProc)%CoverageSendList(SurfExchange%nCoverageSidesSend(iProc)))
      ALLOCATE(SurfCOMM%MPINeighbor(iProc)%CoverageRecvList(SurfExchange%nCoverageSidesRecv(iProc)))
      SurfCOMM%MPINeighbor(iProc)%CoverageRecvList(:)=SurfCOMM%MPINeighbor(iProc)%SendList(:)
      SurfCOMM%MPINeighbor(iProc)%CoverageSendList(:)=SurfCOMM%MPINeighbor(iProc)%RecvList(:)
      SurfCoverageSendBuf(iProc)%content = 0.
      SurfCoverageRecvBuf(iProc)%content = 0.
    END DO
  END IF
#endif /*MPI*/

  IF (DSMC%WallModel.EQ.2) THEN
    CALL CalcDesorbProb()
    CALL CalcAdsorbProb()
#ifdef MPI
    CALL ExchangeCoverageInfo()
#endif /*MPI*/
  END IF
END IF ! SurfMesh%SurfOnProc

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

! define number of possible recombination reactions per species needed for sampling
IF (DSMC%WallModel.EQ.1) Adsorption%RecombNum = 0
IF (DSMC%WallModel.EQ.2) Adsorption%RecombNum = 1
CALL Init_ChemistrySampling()

END SUBROUTINE InitDSMCSurfModel


SUBROUTINE InitSurfCoverage()
!===================================================================================================================================
!> Surface coverage is initialized either from given surface init file or from constant values for each given particle boundary.
!> Therefore, the ini file is checked if state file is specified. If not coverage is initialized from parameter else
!> .h5 file is read and used for init. If file has wrong entries programm aborts.
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_HDF5_INPUT             ,ONLY: DatasetExists,GetDataProps,ReadAttribute,ReadArray,GetDataSize
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_DSMC_Vars              ,ONLY: Adsorption, DSMC!, CollisMode
USE MOD_PARTICLE_Vars          ,ONLY: nSpecies, PDM
USE MOD_ReadInTools            ,ONLY: GETSTR,GETREAL
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, nPartBound, PartBound
USE MOD_Particle_Boundary_Vars ,ONLY: offSetSurfSide, nSurfBC
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
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
!===================================================================================================================================

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
    IF (DSMC%WallModel.EQ.2) THEN
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

  ! check if Dataset SurfaceData exists and read from container
  CALL DatasetExists(File_ID,'SurfaceData',exists)
  IF (exists) THEN
    CALL GetDataProps('SurfaceData',nVarSurf_HDF5,N_HDF5,nSurfSides_HDF5,NodeType_HDF5)
    SWRITE(UNIT_stdOut,'(A,A)')' GET NUMBER AND NAMES OF SURFACE-BC-SIDES IN COVERAGE-INIT FILE... '
    CALL GetDataSize(File_ID,'BC_Surf',nDims,HSize,attrib=.TRUE.)
    nSurfBC_HDF5 = INT(HSize(1),4)
    SWRITE(UNIT_stdOut,'(A3,A45,A3,I33,A13)')' | ','Number of Surface BCs',' | ',nSurfBC_HDF5,' | HDF5    | '
    IF (MPIROOT) THEN
      ALLOCATE(SurfBCName_HDF5(nSurfBC_HDF5))
      CALL ReadAttribute(File_ID,'BC_Surf',nSurfBC_HDF5,StrArray=SurfBCName_HDF5)
      DO iName = 1,nSurfBC_HDF5
        SWRITE(UNIT_stdOut,'(A3,A38,I2.1,A5,A3,A33,A13)')' | ','BC',iName,'Name',' | ',TRIM(SurfBCName_HDF5(iName)),' | HDF5    | '
      END DO
    END IF
    SWRITE(UNIT_stdOut,'(A)')' DONE!'
    IF (nSurfBC_HDF5.NE.nSurfBC) CALL abort(&
__STAMP__&
,'Error in surface coverage init: number of surface boundaries in HDF5-file does not match!')
    IF (nSurfSides_HDF5.NE.SurfMesh%nGlobalSides) CALL abort(&
__STAMP__&
,'Error in surface coverage init: number of surface sides in HDF5-file does not match!')
    ! number comes from boundary sampling, if wallmodel is used then nVar2d_Spec=4
    nVar2D = 5
    nVar2D_Spec = 4
    nVar2D_Total = nVar2D + nVar2D_Spec*nSpecies
    IF (nVarSurf_HDF5.NE.nVar2D_Total) CALL abort(&
__STAMP__&
,'Error in surface coverage init: number of variables in HDF5-file does not match (5+4*nSpecies)!')
    SDEALLOCATE(SurfState_HDF5)
    ALLOCATE(SurfState_HDF5(1:nVarSurf_HDF5,1:nSurfSample,1:nSurfSample,SurfMesh%nSides))
    CALL ReadArray('SurfaceData',4,(/nVarSurf_HDF5,nSurfSample,nSurfSample,SurfMesh%nSides/)&
                   ,offsetSurfSide,4,RealArray=SurfState_HDF5)
    iVar = 3
    DO iSpec = 1, nSpecies
      DO iSurfSide = 1, SurfMesh%nSides
        SideID = Adsorption%SurfSideToGlobSideMap(iSurfSide)
        PartboundID = PartBound%MapToPartBC(BC(SideID))
        IF (PartBound%SolidCatalytic(PartboundID)) THEN
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


SUBROUTINE Init_SurfDist()
!===================================================================================================================================
!> Initializing surface distibution reconstruction model for calculating of coverage effects on heat of adsorption
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
USE MOD_Particle_Vars          ,ONLY: nSpecies, Species!, KeepWallParticles
USE MOD_DSMC_Vars              ,ONLY: Adsorption, SurfDistInfo
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, PartBound
#ifdef MPI
USE MOD_Particle_Boundary_Vars ,ONLY: SurfCOMM
USE MOD_Particle_MPI_Vars      ,ONLY: SurfDistSendBuf, SurfDistRecvBuf, SurfExchange, PartMPI
USE MOD_DSMC_SurfModel_Tools   ,ONLY: ExchangeSurfDistInfo, ExchangeSurfDistSize
#endif /*MPI*/
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
INTEGER                          :: iSurfSide, iSubSurf, jSubSurf, iSpec, iInterAtom
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
      IF (PartBound%SolidCatalytic(PartboundID)) THEN
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
      ELSE !PartBound%SolidCatalytic(PartboundID)
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
IF (.NOT.PartBound%SolidCatalytic(PartboundID)) CYCLE
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
DO iSurfSide = 1,SurfMesh%nSides
SideID = Adsorption%SurfSideToGlobSideMap(iSurfSide)
PartboundID = PartBound%MapToPartBC(BC(SideID))
IF (.NOT.PartBound%SolidCatalytic(PartboundID)) CYCLE
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
      SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%Species(UsedSiteMapPos) = iSpec
      ! assign bond order of respective surface atoms in the surface lattice
      DO iInterAtom = 1,SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%nInterAtom
        xpos = SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%BondAtomIndx(UsedSiteMapPos,iInterAtom)
        ypos = SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(Coord)%BondAtomIndy(UsedSiteMapPos,iInterAtom)
        SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%SurfAtomBondOrder(iSpec,xpos,ypos) = &
          SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%SurfAtomBondOrder(iSpec,xpos,ypos) + 1
      END DO
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

#ifdef MPI
#ifdef CODE_ANALYZE
! write out the number of sites on all surface of the proc, that are considered for adsorption
WRITE(UNIT_stdOut,'(A,I3,I13,A,I13,A,I13)')' | Maximum number of surface sites on proc: ',myRank,Max_Surfsites_num,&
  ' | own: ',Max_Surfsites_own,' | halo: ',Max_Surfsites_halo
#endif /*CODE_ANALYZE*/
CALL MPI_ALLREDUCE(MPI_IN_PLACE,Max_Surfsites_own,1,MPI_INTEGER,MPI_SUM,PartMPI%COMM,iError) ! write only if mpiroot of all comms
SWRITE(UNIT_stdOut,'(A3,A,I0)') ' > ','Surface sites for all catalytic boundaries: ', Max_SurfSites_own

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
  IF(SurfExchange%nSurfDistSidesRecv(iProc).NE.0) THEN
    ALLOCATE(SurfCOMM%MPINeighbor(iProc)%SurfDistRecvList(SurfExchange%nSurfDistSidesRecv(iProc)))
    SurfCOMM%MPINeighbor(iProc)%SurfDistRecvList(:)=SurfCOMM%MPINeighbor(iProc)%SendList(:)
    ALLOCATE(SurfExchange%NbrOfPos(iProc)%nPosRecv(1:SurfExchange%nSurfDistSidesRecv(iProc)))
    SurfExchange%NbrOfPos(iProc)%nPosRecv = 0
  END IF
  IF(SurfExchange%nSurfDistSidesSend(iProc).NE.0) THEN
    ALLOCATE(SurfCOMM%MPINeighbor(iProc)%SurfDistSendList(SurfExchange%nSurfDistSidesSend(iProc)))
    SurfCOMM%MPINeighbor(iProc)%SurfDistSendList(:)=SurfCOMM%MPINeighbor(iProc)%RecvList(:)
    ALLOCATE(SurfExchange%NbrOfPos(iProc)%nPosSend(1:SurfExchange%nSurfDistSidesSend(iProc)))
    SurfExchange%NbrOfPos(iProc)%nPosSend = 0
  END IF
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
!> Initializing surface reaction variables
!===================================================================================================================================
! MODULES
USE MOD_Globals                ,ONLY: abort, MPIRoot, UNIT_StdOut
USE MOD_DSMC_Vars              ,ONLY: Adsorption, SpecDSMC
USE MOD_PARTICLE_Vars          ,ONLY: nSpecies
USE MOD_Particle_Boundary_Vars ,ONLY: SurfMesh
USE MOD_ReadInTools            ,ONLY: GETREAL, GETINT, GETREALARRAY, GETINTARRAY
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

IF (SurfMesh%SurfOnProc .OR. MPIRoot) THEN
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

  MaxDissNum = GETINT('Part-Species-MaxDissNum','0')
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
              GETREALARRAY('Part-SurfChem-Disprop'//TRIM(hilf)//'-DissBond_K-Reactants',2,'0.,0.')
      Adsorption%Product_DissBond_K(:,iReactNum+nDissoc) = &
              GETREALARRAY('Part-SurfChem-Disprop'//TRIM(hilf)//'-DissBond_K-Products',2,'0.,0.')
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
              Adsorption%AdsorpReactInfo(iSpec)%AdsReactCount(1:Adsorption%ReactNum+1),&
              Adsorption%AdsorpReactInfo(iSpec)%SurfReactCount(1:Adsorption%ReactNum+Adsorption%NumOfExchReact+1))
    Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(:) = 0.
    Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(:) = 0.
    Adsorption%AdsorpReactInfo(iSpec)%MeanAdsActE(:) = 0.
    Adsorption%AdsorpReactInfo(iSpec)%MeanSurfActE(:) = 0.
    Adsorption%AdsorpReactInfo(iSpec)%AdsReactCount(:) = 0
    Adsorption%AdsorpReactInfo(iSpec)%SurfReactCount(:) = 0
  END DO
#endif

  CalcTST_Case = GETINT('Particles-DSMC-Adsorption-CalcTST','0')
  ALLOCATE(Adsorption%TST_Calc(0:Adsorption%ReactNum,1:nSpecies))
  Adsorption%TST_Calc(:,:) = .FALSE.
  IF (CalcTST_Case.GT.0) CALL Init_TST_Coeff(CalcTST_Case)
END IF ! SurfMesh%SurfOnProc .OR. MPIRoot

SWRITE(UNIT_stdOut,'(A)')' INIT SURFACE CHEMISTRY DONE!'

END SUBROUTINE Init_SurfChem

SUBROUTINE Init_ChemistrySampling()
!===================================================================================================================================
!> Initializion of additional surface sampling if catalysis is enabled
!===================================================================================================================================
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_DSMC_Vars              ,ONLY: Adsorption
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, SampWall
#ifdef MPI
USE MOD_Particle_Boundary_Vars ,ONLY: SurfCOMM
USE MOD_Particle_MPI_Vars      ,ONLY: SurfSendBuf,SurfRecvBuf,SurfExchange
USE MOD_DSMC_SurfModel_Tools   ,ONLY: ExchangeSurfDistInfo, ExchangeSurfDistSize
#endif /*MPI*/
!===================================================================================================================================
IMPLICIT NONE
!=================================================================================================================================== 
! Local variable declaration
INTEGER                          :: iSide
#ifdef MPI
INTEGER                          :: iProc
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
#ifdef MPI
! Reallocate buffer for mpi communication of sampling
DO iProc=1,SurfCOMM%nMPINeighbors
  SDEALLOCATE(SurfSendBuf(iProc)%content)
  SDEALLOCATE(SurfRecvBuf(iProc)%content)
  IF(SurfExchange%nSidesSend(iProc).GT.0) THEN
    ALLOCATE(SurfSendBuf(iProc)%content((2*nSpecies+1+SurfMesh%SampSize+(Adsorption%RecombNum*nSpecies))&
                                        *(nSurfSample**2)*SurfExchange%nSidesSend(iProc)))
    SurfSendBuf(iProc)%content=0.
  END IF
  IF(SurfExchange%nSidesRecv(iProc).GT.0) THEN
    ALLOCATE(SurfRecvBuf(iProc)%content((2*nSpecies+1+SurfMesh%SampSize+(Adsorption%RecombNum*nSpecies))&
                                        *(nSurfSample**2)*SurfExchange%nSidesRecv(iProc)))
    SurfRecvBuf(iProc)%content=0.
  END IF
END DO ! iProc
#endif

END SUBROUTINE Init_ChemistrySampling

SUBROUTINE Init_TST_Coeff(TST_Case)
!===================================================================================================================================
!> Initializion of the Transition state theory (TST) factors for surface chemistry
!===================================================================================================================================
! MODULES
USE MOD_Globals                ,ONLY: abort, MPIRoot, UNIT_StdOut
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_DSMC_Vars              ,ONLY: Adsorption
USE MOD_PARTICLE_Vars          ,ONLY: nSpecies
USE MOD_Particle_Boundary_Vars ,ONLY: SurfMesh
USE MOD_ReadInTools            ,ONLY: GETREAL
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

Adsorption%PartitionMaxTemp = GETREAL('Particles-DSMC-AdsorptionTST-PartitionMaxTemp','10000.')
Adsorption%PartitionInterval = GETREAL('Particles-DSMC-AdsorptionTST-PartitionInterval','20.')
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


SUBROUTINE FinalizeDSMCSurfModel()
!===================================================================================================================================
!> Deallocate surface model vars
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars              ,ONLY: Adsorption, SurfDistInfo
USE MOD_PARTICLE_Vars          ,ONLY: PDM, PEM
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh
#ifdef MPI
USE MOD_Particle_Boundary_Vars ,ONLY: SurfCOMM
USE MOD_Particle_MPI_Vars      ,ONLY: SurfExchange
USE MOD_Particle_MPI_Vars      ,ONLY: AdsorbSendBuf,AdsorbRecvBuf,SurfDistSendBuf,SurfDistRecvBuf
USE MOD_Particle_MPI_Vars      ,ONLY: SurfCoverageSendBuf,SurfCoverageRecvBuf
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iSubSurf,jSubSurf,iSurfSide,iCoord
#ifdef MPI
INTEGER                      :: iProc
#endif /*MPI*/
!===================================================================================================================================
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
! parameters for Kisliuk and Polanyi Wigner model (wallmodel=1)
SDEALLOCATE(Adsorption%MaxCoverage)
SDEALLOCATE(Adsorption%InitStick)
SDEALLOCATE(Adsorption%PrefactorStick)
SDEALLOCATE(Adsorption%Adsorbexp)
SDEALLOCATE(Adsorption%Nu_a)
SDEALLOCATE(Adsorption%Nu_b)
SDEALLOCATE(Adsorption%DesorbEnergy)
SDEALLOCATE(Adsorption%Intensification)
! parameters for surface recombination model (wallmodel=2)
SDEALLOCATE(Adsorption%RecombCoeff)
SDEALLOCATE(Adsorption%RecombAccomodation)
SDEALLOCATE(Adsorption%RecombEnergy)
SDEALLOCATE(Adsorption%RecombData)
! parameters for UBI-QEP model (wallmodel=3)
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
! surfaces distribution variables (currently wallmodel=3)
IF (ALLOCATED(SurfDistInfo)) THEN
DO iSurfSide=1,SurfMesh%nSides
  DO iSubSurf = 1,nSurfSample
  DO jSubSurf = 1,nSurfSample
#ifdef MPI
    SDEALLOCATE(SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%Nbr_changed)
#endif /*MPI*/
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
  SDEALLOCATE(SurfCOMM%MPINeighbor(iProc)%CoverageSendList)
  SDEALLOCATE(SurfCOMM%MPINeighbor(iProc)%CoverageRecvList)
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
SDEALLOCATE(SurfExchange%nCoverageSidesSend)
SDEALLOCATE(SurfExchange%nCoverageSidesRecv)
IF (ALLOCATED(SurfCoverageSendBuf)) THEN
  DO iProc=1,SurfCOMM%nMPINeighbors
    SDEALLOCATE(SurfCoverageSendBuf(iProc)%content)
  END DO
  DEALLOCATE(SurfDistSendBuf)
END IF
IF (ALLOCATED(SurfCoverageRecvBuf)) THEN
  DO iProc=1,SurfCOMM%nMPINeighbors
    SDEALLOCATE(SurfCoverageRecvBuf(iProc)%content)
  END DO
  DEALLOCATE(SurfDistRecvBuf)
END IF
#endif /*MPI*/

END SUBROUTINE FinalizeDSMCSurfModel

END MODULE MOD_DSMC_SurfModelInit
