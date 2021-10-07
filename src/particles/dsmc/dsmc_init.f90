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

MODULE MOD_DSMC_Init
!===================================================================================================================================
! Initialization of DSMC
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE InitDSMC
  MODULE PROCEDURE InitDSMC
END INTERFACE

INTERFACE FinalizeDSMC
  MODULE PROCEDURE FinalizeDSMC
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: InitDSMC, FinalizeDSMC
!===================================================================================================================================
PUBLIC::DefineParametersDSMC
CONTAINS

!==================================================================================================================================
!> Define parameters for DSMC (Direct Simulation Monte Carlo)
!==================================================================================================================================
SUBROUTINE DefineParametersDSMC()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("DSMC")

CALL prms%CreateIntOption(      'Particles-DSMC-CollisMode'      &
                                        , 'Define mode of collision handling in DSMC.\n'//&
                                            '0: No Collisions (=free molecular flow with DSMC-Sampling-Routines).\n'//&
                                            '1: Elastic Collision \n'//&
                                            '2: Relaxation + Elastic Collision \n'//&
                                            '3: Mode 2 + Chemical Reactions.', '1')
CALL prms%CreateIntOption(      'Particles-DSMC-SelectionProcedure'     &
                                        , 'Mode of Selection Procedure\n'//&
                                          '1: Laux\n'//&
                                          '2: Gimelsheim.', '1')
CALL prms%CreateRealOption(     'Particles-DSMC-RotRelaxProb'&
                                          , 'Define the rotational relaxation probability upon collision of molecules\n'//&
                                          '0-1: constant\n'//&
                                          '2: variable, Boyd)', '0.2')
CALL prms%CreateRealOption(     'Particles-DSMC-VibRelaxProb'&
                                          , 'Define the vibrational relaxation probability upon collision of molecules', '0.004')
CALL prms%CreateRealOption(     'Particles-DSMC-ElecRelaxProb'&
                                          , 'Define the electronic relaxation probability upon collision of molecules', '0.01')
CALL prms%CreateRealOption(     'Particles-DSMC-GammaQuant'&
                                          , 'Set the GammaQuant for zero point energy in Evib (perhaps also Erot) should be'//&
                                          ' 0.5 or 0.', '0.5')
CALL prms%CreateLogicalOption(  'Particles-DSMC-AmbipolarDiffusion', &
                                          'Enables the ambipolar diffusion modelling of electrons, which are attached to the '//&
                                          'ions, however, retain their own velocity vector to participate in collision events.',&
                                          '.FALSE.')
!-----------------------------------------------------------------------------------
CALL prms%CreateLogicalOption(  'Particles-DSMC-CalcQualityFactors', &
                                          'Enables [TRUE] / disables [FALSE] the calculation and output of:\n'//&
                                          'Maximal collision probability\n'//&
                                          'Time-averaged mean collision probability\n'//&
                                          'Mean collision separation distance over mean free path' , '.FALSE.')
CALL prms%CreateLogicalOption(  'Particles-DSMCReservoirSim', &
                                          'Only TD=Reservoir (42).\n'//&
                                          'Set [TRUE] to disable particle movement. Use for reservoir simulations.' , '.FALSE.')
CALL prms%CreateLogicalOption(  'Particles-DSMCReservoirSimRate', &
                                          'Only TD=Reservoir (42).\n'//&
                                          'Set [TRUE] to disable particle reactions. Only probabilities (rates) are calculated.', &
                                          '.FALSE.')
CALL prms%CreateLogicalOption(  'Particles-DSMCReservoirStatistic', &
                                          'Only TD=Reservoir (42).\n'//&
                                          'Probabilities (rates) are calculated\n'//&
                                          ' [TRUE] counting reacting particles.\n'//&
                                          ' [FALSE] summing reaction probabilities (does not work with Q-K).' , '.FALSE.')
CALL prms%CreateLogicalOption(  'Particles-DSMC-TEVR-Relaxation'&
                                          , 'Flag for Translational-Vibrational-Electric-Rotational relaxation (T-V-E-R)\n'//&
                                          '[TRUE] or more simple T-V-R T-E-R\n'//&
                                          '[FALSE] relaxation.' , '.FALSE.')
CALL prms%CreateIntOption(  'Particles-DSMC-ElectronicModel', &
                                          'Select model for the electronic states of atoms and molecules:\n'//&
                                          '0: No electronic energy treatment [default]\n'//&
                                          '1: Model by Liechty, each particle has a specific electronic state\n'//&
                                          '2: Model by Burt, each particle has an electronic distribution function', '0')
CALL prms%CreateStringOption(   'Particles-DSMCElectronicDatabase'&
                                          , 'If electronic model is used give (relative) path to (h5) Name of Electronic State'//&
                                          ' Database', 'none')
CALL prms%CreateRealOption(     'EpsMergeElectronicState'&
                                         , 'Percentage parameter of electronic energy level merging.' , '1E-4')
CALL prms%CreateLogicalOption(  'Particles-DSMC-PolyRelaxSingleMode'&
                                         , 'Set [TRUE] for separate relaxation of each vibrational mode of a polyatomic in a '//&
                                           'loop over all vibrational modes.\n'//&
                                           'Every mode has its own corrected relaxation probability, comparison with the '//&
                                           'same random number while the previous probability is added to the next', '.FALSE.')
CALL prms%CreateLogicalOption(  'Particles-DSMC-CompareLandauTeller'&
                                         ,'Only TD=Reservoir (42). ', '.FALSE.')
CALL prms%CreateLogicalOption(  'Particles-DSMC-UseOctree'&
                                         ,'Use octree method for dynamic grid resolution based on the current mean free path '//&
                                          'and the particle number', '.FALSE.')
CALL prms%CreateIntOption(      'Particles-OctreePartNumNode'&
                                         ,'Resolve grid until the maximum number of particles in a subcell equals'//&
                                          ' OctreePartNumNode')
CALL prms%CreateIntOption(      'Particles-OctreePartNumNodeMin'&
                                         ,'Allow grid division until the minimum number of particles in a subcell is above '//&
                                          'OctreePartNumNodeMin')
CALL prms%CreateLogicalOption(  'Particles-DSMC-UseNearestNeighbour'&
                                         ,'Enable/disable the nearest neighbour search algorithm','.FALSE.')
CALL prms%CreateLogicalOption(  'Particles-DSMC-ProhibitDoubleCollisions'&
                                         ,'2D/Axisymmetric only: Prohibit the occurrence of repeated collisions between the '//&
                                          'same particle pairs in order to reduce the statistical dependence')
CALL prms%CreateLogicalOption(  'Particles-DSMC-MergeSubcells'&
                                         ,'2D/Axisymmetric only: Merge subcells divided by the quadtree algorithm to satisfy '//&
                                          'the minimum particle per subcell requirement', '.FALSE.')

CALL prms%SetSection("DSMC Collision")
CALL prms%CreateIntOption(      'Particles-DSMC-crossSectionConstantMode'  &
                                            ,'Flags which cross section (sigma) constant Cab calculation mode is used.\n'//&
                                            ' sigma=Cab * cr^(-2 omega)\n'//&
                                            ' 0 : single omega for the computational domain\n'//&
                                            '     Part-Species1=omega will be set \n     for all and Cab will be calculated\n'//&
                                            '     via species-specific factor A_j\n'//&
                                            '     Cab=(A_1+A_2)^2*m_red^-omega\n     (see Bird 1981 and Laux 1996)\n'//&
                                            ' 1 : Cab will be calculated\n     via species-specific factor A_j\n'//&
                                            ' 2 : Cab will be calculated directly\n     (see Bird 1981 eq (9):constants)', '0')
CALL prms%CreateLogicalOption(   'Particles-DSMC-averagedCollisionParameters'  &
                                           ,'Flags if collision parameters are specific to colliding species'//&
                                            ' or species specific and averaged for the collision itself.\n'//&
                                            ' T: Part-Species[$]-omega,-Tref,-dref,-alphaVSS\n'//&
                                            '    species-specific parameters(T)\n    can be found in tables e.g. in\n'//&
                                            '    VHS: bird1994 (table A1 and A2)\n    VSS: bird1994 (table A1 and A3)\n'//&
                                            '         weaver2015\n         (https://doi.org/10.1063/1.4921245)\n'//&
                                            ' F: Part-Collision[$]-omega,-Tref,-dref,-alphaVSS\n'//&
                                            '    collision-specific parameters(F)\n    can be found in tables e.g. in\n'//&
                                            '    VHS/VSS: krishnan2015\n             (https://doi.org/10.2514/6.2015-3373)\n'//&
                                            '    VHS/VSS: krishnan2016\n             (https://doi.org/10.1063/1.4939719)', 'T')
CALL prms%CreateIntOption(   'Part-Collision[$]-partnerSpecies[$]'  &
                                           ,'Colliding partnerSpecies(1,2) equal to SpeciesID from ini.\n'//&
                                            'e.g: Collision ID=1 Ar+NO\n     Part-Species2(Ar) + Part-Species3(NO)\n'//&
                                            '     write Part-Collision1-partnerSpecies1=2\n'//&
                                            '     and   Part-Collision1-partnerSpecies2=3', '0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Collision[$]-Tref'  &
                                           ,'collision parameter: collision-specific reference temperature for VHS/VSS model.'&
                                           , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Collision[$]-dref'  &
                                           ,'collision parameter: collision-specific reference diameter for VHS/VSS model. '&
                                           , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Collision[$]-omega'  &
                                           ,'collision parameter: collision-specific\n'//&
                                            'temperature exponent omega=2/(eta-1)\nfor VHS/VSS model.\n'//&
                                            'omega=0 + alpha=1\nreproduces HS model\n' //&
                                            'CAUTION: omega=omega_bird1994-0.5', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Collision[$]-alphaVSS'  &
                                           ,'collision parameter: collision-specific scattering\n'//&
                                            'exponent alpha for VSS model\n'//&
                                            'default alphaVSS=1\nreproduces VHS model\n'//&
                                            '(see Bird 1994 p.42 for more information.)'&
                                           , '1.', numberedmulti=.TRUE.)

CALL prms%SetSection("DSMC Species")

CALL prms%CreateStringOption(  'Part-Species[$]-SpeciesName'  &
                                         ,'Species name of Species[$]', 'none', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(     'Part-Species[$]-InteractionID' , 'ID for identification of particles \n'//&
                                                                 '  1: Atom\n'//&
                                                                 '  2: Molecule\n'//&
                                                                 '  4: Electron\n'//&
                                                                 ' 10: Atomic Ion\n'//&
                                                                 ' 20: Molecular Ion\n'//&
                                                                 ' 40: Excited Atom\n'//&
                                                                 '100: Excited Atomic Ion\n'//&
                                                                 '200: Excited Molecule\n'//&
                                                                 '400: Excited Molecular Ion)', '0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Tref'  &
                                           ,'collision parameter: species-specific reference temperature for VHS/VSS model.' &
                                           , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-dref' &
                                           ,'collision parameter: species-specific reference diameter for VHS/VSS model.'&
                                           , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-omega'  &
                                           ,'collision parameter: species-specific temperature exponent omega = 2 / (eta - 1)'//&
                                            ' for VHS/VSS model.\nCAUTION: omega = omega_bird1994 - 0.5'&
                                           , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-alphaVSS'&
                                           ,'collision parameter: species-specific scattering exponent alpha for VSS model\n'//&
                                            'default alphaVSS=1 reproduces VHS model\n'//&
                                            '(see Bird 1994 p.42 for more information)', '1.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-CharaTempVib','Characteristic vibrational temperature.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-CharaTempRot'  &
                                           ,'Characteristic rotational temperature', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Ediss_eV','Energy of Dissoziation in [eV].', numberedmulti=.TRUE.)
! ----------------------------------------------------------------------------------------------------------------------------------
CALL prms%CreateLogicalOption(  'Particles-DSMC-useRelaxProbCorrFactor'&
                                           ,'Use the relaxation probability correction factor of Lumpkin', '.FALSE.')
CALL prms%CreateRealOption(     'Part-Species[$]-CollNumRotInf'  &
                                           ,'Collision number for rotational relaxation according to Parker or'//&
                                            'Zhang, ini_2 -> model dependent!', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-TempRefRot'  &
                                           ,'Reference temperature for rotational relaxation according to Parker or '//&
                                            'Zhang, ini_2 -> model dependent!', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-MWConstA-[$]-[$]'  &
                                           ,'Millikan-White constant A for variable vibrational relaxation probability, ini_2' &
                                           ,'0.0',  numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-MWConstB-[$]-[$]'  &
                                           ,'Millikan-White constant B for variable vibrational relaxation probability, ini_2' &
                                           ,'0.0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Particles-DSMC-alpha'  &
                                           ,'Relaxation factor of vibration probability for VibRelaxProb = 2' &
                                           ,'0.99')
CALL prms%CreateRealOption(     'Part-Species[$]-VibCrossSection'  &
                                           , 'Vibrational collision cross-section to Boyd, ini_2','1.E-19', numberedmulti=.TRUE.)
! ----------------------------------------------------------------------------------------------------------------------------------
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-TempVib'  &
                                           ,'Vibrational temperature.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-TempRot'  &
                                           ,'Rotational temperature.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-TempElec'  &
                                           ,'Electronic temperature.', numberedmulti=.TRUE.)

CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-TempVib'  &
                                           ,'Vibrational temperature.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-TempRot'  &
                                           ,'Rotational temperature.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-TempElec'  &
                                           ,'Electronic temperature.', numberedmulti=.TRUE.)

CALL prms%CreateRealOption(     'Part-Species[$]-HeatOfFormation_K'  &
                                           ,'Heat of formation of the respective species [Kelvin]'&
                                           , numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-PreviousState'  &
                                           ,'Species number of the previous state (e.g. N for NIon) ', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Species[$]-FullyIonized'  &
                                           ,'Flag if the species is fully ionized, e.g., C^6+ ', '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-NextIonizationSpecies'  &
                                           , 'SpeciesID of the next higher ionization level (required for field ionization)'&
                                           , '0', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-NumElectronicLevels'  &
                                           ,'Max elec quantum number + 1 ', '0', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-ElectronicDegeneracy-Level[$]'  &
                                           ,'Electronic degeneracy level of respective species', '0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-ElectronicEnergyLevel-Level[$]'  &
                                           ,'Electronic energy level of respective species', '0.', numberedmulti=.TRUE.)

CALL prms%SetSection("DSMC Species Polyatomic")
CALL prms%CreateLogicalOption(  'Part-Species[$]-PolyatomicMol'  &
                                           ,'Allows the usage of polyatomic molecules (3 or more atoms).', '.FALSE.' &
                                           , numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Species[$]-LinearMolec'  &
                                           ,'Flag if the polyatomic molecule is a linear molecule (e.g. CO2 is linear, while '//&
                                            'H2O is not.)', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-NumOfAtoms'  &
                                           ,'Number of atoms in the molecule (e.g. CH4 -> 5 atoms).', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-CharaTempVib[$]'  &
                                           ,'Characteristic vibrational temperature [K], given per mode. Degenerate modes should '//&
                                            'simply be given repeatedly, corresponding to the degeneracy.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-CharaTempRot[$]'  &
                                           ,'Characteristic rotational temperature [K]. Linear molecules require only a single '//&
                                            'input, while non-linear molecules require three.', '0.', numberedmulti=.TRUE.)

CALL prms%CreateLogicalOption(  'Part-Species[$]-UseCollXSec'  &
                                           ,'Utilize collision cross sections for the determination of collision probabilities' &
                                           ,'.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Species[$]-UseVibXSec'  &
                                           ,'Utilize vibrational cross sections for the determination of relaxation probabilities' &
                                           ,'.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Particles-CollXSec-Database', 'File name for the collision cross section database. Container '//&
                                                               'should be named with species pair (e.g. "Ar-electron"). The '//&
                                                               'first column shall contain the energy in eV and the second '//&
                                                               'column the cross-section in m^2', 'none')
CALL prms%CreateLogicalOption(  'Particles-CollXSec-NullCollision'  &
                                  ,'Utilize the null collision method for the determination of the number of pairs '//&
                                  'based on the maximum collision frequency and time step (only with a background gas)' &
                                  ,'.TRUE.')

END SUBROUTINE DefineParametersDSMC

SUBROUTINE InitDSMC()
!===================================================================================================================================
! Init of DSMC Vars
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_ReadInTools
USE MOD_DSMC_Vars
USE MOD_Mesh_Vars              ,ONLY: nElems, NGEo
USE MOD_Globals_Vars           ,ONLY: Pi, BoltzmannConst, ElementaryCharge
USE MOD_Particle_Vars          ,ONLY: nSpecies, Species, PDM, PartSpecies, Symmetry, VarTimeStep
USE MOD_Particle_Vars          ,ONLY: DoFieldIonization
USE MOD_DSMC_ParticlePairing   ,ONLY: DSMC_init_octree
USE MOD_DSMC_ChemInit          ,ONLY: DSMC_chemical_init
USE MOD_DSMC_PolyAtomicModel   ,ONLY: InitPolyAtomicMolecs, DSMC_SetInternalEnr_Poly
USE MOD_DSMC_SpecXSec          ,ONLY: MCC_Init
USE MOD_DSMC_CollisVec         ,ONLY: DiceDeflectedVelocityVector4Coll, DiceVelocityVector4Coll, PostCollVec
USE MOD_part_emission_tools    ,ONLY: DSMC_SetInternalEnr_LauxVFD
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(32)         :: hilf , hilf2
INTEGER               :: iCase, iSpec, jSpec, nCase, iPart, iInit, iDOF, VarNum
INTEGER               :: iColl, jColl, pColl, nCollision ! for collision parameter read in
REAL                  :: A1, A2     ! species constant for cross section (p. 24 Laux)
LOGICAL               :: PostCollPointerSet
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' DSMC INIT ...'

! Initialize counter (Count the number of ReactionProb>1)
ReactionProbGTUnityCounter = 0

! reading and reset general DSMC values
CollisMode = GETINT('Particles-DSMC-CollisMode','1') !0: no collis, 1:elastic col, 2:elast+rela, 3:chem
SelectionProc = GETINT('Particles-DSMC-SelectionProcedure','1') !1: Laux, 2:Gimelsheim
IF(RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
  IF(SelectionProc.NE.1) THEN
    CALL abort(__STAMP__&
        ,'ERROR: Radial weighting or variable time step is not implemented with the chosen SelectionProcedure: ' &
        ,IntInfoOpt=SelectionProc)
  END IF
END IF

DSMC%MergeSubcells = GETLOGICAL('Particles-DSMC-MergeSubcells')
IF(DSMC%MergeSubcells.AND.(Symmetry%Order.NE.2)) THEN
  CALL abort(__STAMP__&
      ,'ERROR: Merging of subcells only supported within a 2D/axisymmetric simulation!')
END IF
  IF(CollisMode.GE.2) THEN
    DSMC%RotRelaxProb = GETREAL('Particles-DSMC-RotRelaxProb')
    DSMC%VibRelaxProb = GETREAL('Particles-DSMC-VibRelaxProb')
  ELSE
    DSMC%RotRelaxProb = 0.
    DSMC%VibRelaxProb = 0.
  END IF
DSMC%ElecRelaxProb = GETREAL('Particles-DSMC-ElecRelaxProb')
DSMC%GammaQuant   = GETREAL('Particles-DSMC-GammaQuant')
!-----------------------------------------------------------------------------------
DSMC%CalcQualityFactors = GETLOGICAL('Particles-DSMC-CalcQualityFactors')
DSMC%ReservoirSimu = GETLOGICAL('Particles-DSMCReservoirSim')
IF (DSMC%CalcQualityFactors.AND.(CollisMode.LT.1)) THEN
  CALL abort(&
      __STAMP__&
      ,'ERROR: Do not use DSMC%CalcQualityFactors for CollisMode < 1')
END IF ! DSMC%CalcQualityFactors.AND.(CollisMode.LT.1)
DSMC%ReservoirSimuRate       = GETLOGICAL('Particles-DSMCReservoirSimRate')
DSMC%ReservoirRateStatistic  = GETLOGICAL('Particles-DSMCReservoirStatistic')
DSMC%DoTEVRRelaxation        = GETLOGICAL('Particles-DSMC-TEVR-Relaxation')
IF(RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
  IF(DSMC%DoTEVRRelaxation) THEN
    CALL abort(__STAMP__&
        ,'ERROR: Radial weighting or variable time step is not implemented with T-E-V-R relaxation!')
  END IF
END IF
DSMC%ElectronicModel         = GETINT('Particles-DSMC-ElectronicModel')
IF (DSMC%ElectronicModel.EQ.2) THEN
  IF(.NOT.ALLOCATED(ElectronicDistriPart)) ALLOCATE(ElectronicDistriPart(PDM%maxParticleNumber))
END IF
DSMC%ElectronicModelDatabase = TRIM(GETSTR('Particles-DSMCElectronicDatabase','none'))
IF ((DSMC%ElectronicModelDatabase .NE. 'none').AND.&
    ((CollisMode .GT. 1).OR.(CollisMode .EQ. 0))) THEN ! CollisMode=0 is for use of in PIC simulation without collisions
  DSMC%EpsElecBin = GETREAL('EpsMergeElectronicState','1E-4')
ELSEIF(DSMC%ElectronicModel.GT.0) THEN
  CALL Abort(&
      __STAMP__,&
      'ERROR: Electronic model requires a electronic levels database and CollisMode > 1!')
END IF
DSMC%NumPolyatomMolecs = 0
SamplingActive = .FALSE.

IF(DSMC%CalcQualityFactors) THEN
  ! 1: Maximal collision probability per cell/subcells (octree)
  ! 2: Mean collision probability within cell
  ! 3: Mean collision separation distance over mean free path
  ! 4: Counter (is not simply the number of iterations in case of a coupled BGK/FP-DSMC simulation)
  VarNum = 4
  ! VarNum + 1: Number of cloned particles per cell
  ! VarNum + 2: Number of identical particles (no relative velocity)
  IF(RadialWeighting%PerformCloning) VarNum = VarNum + 2
  ALLOCATE(DSMC%QualityFacSamp(nElems,VarNum))
  DSMC%QualityFacSamp(1:nElems,1:VarNum) = 0.0
END IF

! definition of DSMC particle values
ALLOCATE(DSMC_RHS(1:3,1:PDM%maxParticleNumber))
DSMC_RHS = 0

IF (nSpecies.LE.0) THEN
  CALL Abort(&
      __STAMP__&
      ,"ERROR: nSpecies .LE. 0:", nSpecies)
END IF

!-----------------------------------------------------------------------------------------------------------------------------------
! reading in collision model variables
!-----------------------------------------------------------------------------------------------------------------------------------
  ! Either CollisMode.GT.0 or without chemical reactions due to collisions but with field ionization
  IF(DoFieldIonization.OR.CollisMode.NE.0) THEN
    ! Flags for collision parameters
    CollInf%averagedCollisionParameters     = GETLOGICAL('Particles-DSMC-averagedCollisionParameters')
    CollInf%crossSectionConstantMode        = GETINT('Particles-DSMC-crossSectionConstantMode','0')
    ALLOCATE(SpecDSMC(nSpecies))
    DO iSpec = 1, nSpecies
      WRITE(UNIT=hilf,FMT='(I0)') iSpec
      SpecDSMC(iSpec)%Name    = TRIM(GETSTR('Part-Species'//TRIM(hilf)//'-SpeciesName','none'))
      SpecDSMC(iSpec)%InterID = GETINT('Part-Species'//TRIM(hilf)//'-InteractionID','0')
      ! averagedCollisionParameters set true: species-specific collision parameters get read in
      IF(CollInf%averagedCollisionParameters) THEN
        SpecDSMC(iSpec)%Tref         = GETREAL('Part-Species'//TRIM(hilf)//'-Tref'     )
        SpecDSMC(iSpec)%dref         = GETREAL('Part-Species'//TRIM(hilf)//'-dref'     )
        SpecDSMC(iSpec)%omega    = GETREAL('Part-Species'//TRIM(hilf)//'-omega'    )
        SpecDSMC(iSpec)%alphaVSS     = GETREAL('Part-Species'//TRIM(hilf)//'-alphaVSS' )
        ! check for faulty parameters
        IF((SpecDSMC(iSpec)%InterID * SpecDSMC(iSpec)%Tref * SpecDSMC(iSpec)%dref * SpecDSMC(iSpec)%alphaVSS) .EQ. 0) THEN
          CALL Abort(&
          __STAMP__&
          ,'ERROR in species data: check collision parameters in ini \n'//&
           'Part-Species'//TRIM(hilf)//'-(InterID * Tref * dref * alphaVSS) .EQ. 0'//&
           ' - but must not be 0')
        END IF ! (Tref * dref * alphaVSS) .EQ. 0
        IF ((SpecDSMC(iSpec)%alphaVSS.LT.0.0) .OR. (SpecDSMC(iSpec)%alphaVSS.GT.2.0)) THEN
          CALL Abort(&
          __STAMP__&
          ,'ERROR: Check set parameter Part-Species'//TRIM(hilf)//'-alphaVSS must not be lower 0 or greater 2')
        END IF ! alphaVSS parameter check
      END IF ! averagedCollisionParameters
      SpecDSMC(iSpec)%FullyIonized  = GETLOGICAL('Part-Species'//TRIM(hilf)//'-FullyIonized')
      IF(SpecDSMC(iSpec)%InterID.EQ.4) THEN
        DSMC%ElectronSpecies = iSpec
      END IF ! interID .EQ.4
      ! reading electronic state informations from HDF5 file
      IF((DSMC%ElectronicModelDatabase.NE.'none').AND.(SpecDSMC(iSpec)%InterID.NE.4)) CALL SetElectronicModel(iSpec)
    END DO ! iSpec = nSpecies

    ! determine number of different species combinations and allocate collidingSpecies array
    nCollision=0
    DO iColl=1,nSpecies
      DO jColl=iColl,nSpecies
        nCollision=nCollision+1
      END DO !jColl = nSpecies
    END DO !iColl = nSpecies
    ALLOCATE(CollInf%collidingSpecies(nCollision,2))
    CollInf%collidingSpecies(:,:) = 0 ! default value to determine if collidingSpecies are all set

    IF(CollInf%averagedCollisionParameters) THEN ! partnerSpecies for collidingSpecies are set
      iColl = 0
      DO iSpec = 1 , nSpecies
        DO jSpec = iSpec , nSpecies
          iColl = iColl + 1 ! sum up total number of collisions
          CollInf%collidingSpecies(iColl,1) = iSpec ! assign partnerSpecies for iColl's collision
          CollInf%collidingSpecies(iColl,2) = jSpec
        END DO ! jSpec = nSpecies
      END DO ! iSpec = nSpecies
    ELSE ! .NOT. averagedCollisionParameters     : partnerSpecies for collidingSpecies per collision are read in
      DO iColl = 1, nCollision
        WRITE(UNIT=hilf,FMT='(I0)')  iColl
        CollInf%collidingSpecies(iColl,:) = GETINTARRAY('Part-Collision'//TRIM(hilf)//'-partnerSpecies',2,'0,0')
      END DO ! iColl = nCollision
    END IF ! averagedCollisionParameters
    DO iColl = 1, nCollision ! check if any collidingSpecies pair is set multiple times
      WRITE(UNIT=hilf,FMT='(I0)') iColl
      DO pColl = 1,2 ! collision partner
        WRITE (UNIT = hilf2,FMT = '(I0)') pColl
        IF (CollInf%collidingSpecies(iColl,pColl) .EQ. 0) THEN
            CALL Abort(&
            __STAMP__&
            ,'ERROR: Partner species '//TRIM(hilf2)//' for Collision'//TRIM(hilf)//' not defined. '// &
            'Part-Collision'//TRIM(hilf)//'-partnerSpecies required ')
        END IF ! collidingSpecies .EQ. 0
        IF (CollInf%collidingSpecies(iColl,pColl).GT.nSpecies) THEN
          CALL Abort(&
              __STAMP__&
              ,'ERROR: Partner species '//TRIM(hilf2)//' for Collision'//TRIM(hilf)//' .GT. nSpecies')
        END IF
      END DO ! pColl = 2
      DO jColl=1, nCollision
        WRITE(UNIT=hilf2,FMT='(I0)') jColl
        IF ((CollInf%collidingSpecies(iColl,1) .EQ. CollInf%collidingSpecies(jColl,2))  .AND. &
            (CollInf%collidingSpecies(iColl,2) .EQ. CollInf%collidingSpecies(jColl,1))) THEN
          IF (iColl.NE.jColl) THEN
              CALL Abort(&
              __STAMP__&
              ,'ERROR: Partner species for Collision'//TRIM(hilf)//' .EQ. Collision'//TRIM(hilf2))
          END IF ! iColl .EQ. jColl
        END IF ! check for redundant collision partner combination
      END DO !jColl = nColl
    END DO ! iColl = nColl

    ! allocate and initialize collision parameter arrays
    ALLOCATE(CollInf%Tref(nSpecies,nSpecies))
    ALLOCATE(CollInf%dref(nSpecies,nSpecies))
    ALLOCATE(CollInf%omega(nSpecies,nSpecies))
    ALLOCATE(CollInf%alphaVSS(nSpecies,nSpecies))

    ! read collision parameters in and check if all are set
    DO iColl = 1, nCollision
      iSpec = MINVAL (CollInf%collidingSpecies(iColl,:)) ! sorting for filling upper
      jSpec = MAXVAL (CollInf%collidingSpecies(iColl,:)) ! triangular matrix
      WRITE(UNIT=hilf,FMT='(I0)')  iColl
      IF(CollInf%averagedCollisionParameters) THEN ! collision-averaged parameters
        CollInf%Tref      (iSpec,jSpec) = 0.5 * (SpecDSMC(iSpec)%Tref      + SpecDSMC(jSpec)%Tref)
        CollInf%dref      (iSpec,jSpec) = 0.5 * (SpecDSMC(iSpec)%dref      + SpecDSMC(jSpec)%dref)
        CollInf%omega (iSpec,jSpec) = 0.5 * (SpecDSMC(iSpec)%omega + SpecDSMC(jSpec)%omega)
        CollInf%alphaVSS  (iSpec,jSpec) = 0.5 * (SpecDSMC(iSpec)%alphaVSS  + SpecDSMC(jSpec)%alphaVSS)
      ELSE ! collision-specific parameters
        CollInf%Tref      (iSpec,jSpec) = GETREAL('Part-Collision'//TRIM(hilf)//'-Tref'     )
        CollInf%dref      (iSpec,jSpec) = GETREAL('Part-Collision'//TRIM(hilf)//'-dref'     )
        CollInf%omega (iSpec,jSpec) = GETREAL('Part-Collision'//TRIM(hilf)//'-omega'    )
        CollInf%alphaVSS  (iSpec,jSpec) = GETREAL('Part-Collision'//TRIM(hilf)//'-alphaVSS' )
      END IF ! averagedCollisionParameters
      IF (iSpec.NE.jSpec) THEN ! fill lower triangular matrix
        CollInf%Tref      (jSpec,iSpec) = CollInf%Tref      (iSpec,jSpec)
        CollInf%dref      (jSpec,iSpec) = CollInf%dref      (iSpec,jSpec)
        CollInf%omega (jSpec,iSpec) = CollInf%omega (iSpec,jSpec)
        CollInf%alphaVSS  (jSpec,iSpec) = CollInf%alphaVSS  (iSpec,jSpec)
      END IF ! filled lower triangular matrix
      IF(CollInf%dref(iSpec,jSpec) * CollInf%Tref(iSpec,jSpec) * CollInf%alphaVSS(iSpec,jSpec) .EQ. 0) THEN
        CALL Abort(&
        __STAMP__&
        ,'ERROR: Check collision parameters! (Part-Collision'//TRIM(hilf)//'-Tref * dref * alphaVSS) .EQ. 0 - but must not be 0)')
      END IF ! check if collision parameters are set
      IF ((CollInf%alphaVSS(iSpec,jSpec).LT.1) .OR. (CollInf%alphaVSS(iSpec,jSpec).GT.2)) THEN
        CALL Abort(&
        __STAMP__&
        ,'ERROR: Check set parameter Part-Collision'//TRIM(hilf)//'-alphaVSS must not be lower 1 or greater 2')
      END IF ! alphaVSS parameter check
    END DO ! iColl=nColl

    PostCollVec => DiceVelocityVector4Coll
    PostCollPointerSet = .FALSE.
    DO iSpec = 1 , nSpecies
      DO jSpec = iSpec , nSpecies
        IF ( CollInf%alphaVSS(iSpec,jSpec).NE.1.0 ) THEN
          PostCollVec => DiceDeflectedVelocityVector4Coll
          PostCollPointerSet = .TRUE.
          EXIT
        END IF
      END DO ! jSpec = nSpecies
      IF (PostCollPointerSet) EXIT
    END DO ! iSpec = nSpecies

    IF(CollInf%crossSectionConstantMode.EQ.0) THEN ! one omega for all - DEFAULT
      CollInf%omega(:,:)=CollInf%omega(1,1)
    END IF ! CollInf%crossSectionConstantMode=0
  END IF ! DoFieldIonization.OR.CollisMode.NE.0

! allocate internal energy arrays
IF (DSMC%ElectronicModel.GT.0) THEN
  ALLOCATE(PartStateIntEn(1:3,PDM%maxParticleNumber))
ELSE
  ALLOCATE(PartStateIntEn(1:2,PDM%maxParticleNumber))
ENDIF
PartStateIntEn = 0. ! nullify

IF (CollisMode.EQ.0) THEN
#if (PP_TimeDiscMethod==42)
  CALL Abort(&
      __STAMP__&
      , "Free Molecular Flow (CollisMode=0) is not supported for reservoir!")
#endif
ELSE !CollisMode.GT.0
  ! species and case assignment arrays
  ALLOCATE(CollInf%Coll_Case(nSpecies,nSpecies))
  iCase = 0
  DO iSpec = 1, nSpecies
    DO jSpec = iSpec, nSpecies
      iCase = iCase + 1
      CollInf%Coll_Case(iSpec,jSpec) = iCase
      CollInf%Coll_Case(jSpec,iSpec) = iCase
    END DO
  END DO
  nCase = iCase
  CollInf%NumCase = nCase
  ALLOCATE(DSMC%NumColl(nCase +1))
  DSMC%NumColl = 0.
  ALLOCATE(CollInf%Coll_CaseNum(nCase))
  CollInf%Coll_CaseNum = 0
  ALLOCATE(CollInf%Coll_SpecPartNum(nSpecies))
  CollInf%Coll_SpecPartNum = 0.
  ALLOCATE(CollInf%MeanMPF(nCase))
  CollInf%MeanMPF = 0.
  ALLOCATE(CollInf%FracMassCent(nSpecies, nCase)) ! Calculation of mx/(mx+my) and reduced mass
  CollInf%FracMassCent = 0
  ALLOCATE(CollInf%MassRed(nCase))
  CollInf%MassRed = 0
  DO iSpec = 1, nSpecies
    DO jSpec = iSpec, nSpecies
      iCase = CollInf%Coll_Case(iSpec,jSpec)
      CollInf%FracMassCent(iSpec, iCase)  = Species(iSpec)%MassIC &
          / (Species(iSpec)%MassIC + Species(jSpec)%MassIC)
      CollInf%FracMassCent(jSpec, iCase)  = Species(jSpec)%MassIC &
          / (Species(iSpec)%MassIC + Species(jSpec)%MassIC)
      CollInf%MassRed(iCase) = (Species(iSpec)%MassIC*Species(jSpec)%MassIC) &
          / (Species(iSpec)%MassIC+Species(jSpec)%MassIC)
    END DO
  END DO
!-----------------------------------------------------------------------------------------------------------------------------------
! Factor calculation for particle collision
!-----------------------------------------------------------------------------------------------------------------------------------
  ALLOCATE(CollInf%Cab(nCase))
  ALLOCATE(CollInf%KronDelta(nCase))
  CollInf%Cab = 0
  CollInf%KronDelta = 0

  DO iSpec = 1, nSpecies
    DO jSpec = iSpec, nSpecies
      iCase = CollInf%Coll_Case(iSpec,jSpec)
      IF (iSpec.eq.jSpec) THEN
        CollInf%KronDelta(iCase) = 1
      ELSE !ispec.NE.jspec
        CollInf%KronDelta(iCase) = 0
      END IF !kronecker delta
      ! Laux (2.37) prefactor Cab calculation depending on omega coll-averaged or -specific
      SELECT CASE (CollInf%crossSectionConstantMode) ! sigma=Cab * cr^(-2 omega) , see Bird1981 for details
      CASE (0,1)
        ! A1,A2 species constants, see laux1996 (2.38),(2.39)
        A1 = 0.5 * SQRT(Pi) * CollInf%dref(iSpec,iSpec) * (2 * BoltzmannConst * CollInf%Tref(iSpec , iSpec)) ** &
           (CollInf%omega(iSpec , iSpec) * 0.5) / SQRT (GAMMA(2.0 - CollInf%omega(iSpec , iSpec)))
        A2 = 0.5 * SQRT(Pi) * CollInf%dref(jSpec , jSpec)*(2 * BoltzmannConst * CollInf%Tref(jSpec , jSpec)) ** &
           (CollInf%omega(jSpec , jSpec)*0.5) / SQRT (GAMMA(2.0 - CollInf%omega(jSpec , jSpec)))
        CollInf%Cab(iCase) = (A1 + A2) ** 2 * CollInf%MassRed(iCase) ** ( - CollInf%omega(iSpec , jSpec))

      CASE (2) ! cross section constant Cab without Laux simplification (needed if multiple omegas are used) see bird1981 (9)
        CollInf%Cab(iCase) = (SQRT (Pi) * CollInf%dref(iSpec , jSpec) *                                      &
                                  (2 * BoltzmannConst * CollInf%Tref(iSpec , jSpec)) ** (CollInf%omega(iSpec , jSpec) * 0.5) &
                                  / SQRT (GAMMA (2.0 - CollInf%omega(iSpec , jSpec)))) ** 2                                  &
                                  * CollInf%MassRed(iCase)** ( - CollInf%omega(iSpec , jSpec))
      END SELECT !crossSectionConstant
    END DO !jspec=nspecies
  END DO ! ispec=nspecies
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! Collision cross-sections (required for MCC treatment of particles)
  !-----------------------------------------------------------------------------------------------------------------------------------
  DO iSpec = 1, nSpecies
    WRITE(UNIT=hilf,FMT='(I0)') iSpec
    SpecDSMC(iSpec)%UseCollXSec=GETLOGICAL('Part-Species'//TRIM(hilf)//'-UseCollXSec')
    SpecDSMC(iSpec)%UseVibXSec=GETLOGICAL('Part-Species'//TRIM(hilf)//'-UseVibXSec')
    IF(SpecDSMC(iSpec)%UseCollXSec.AND.BGGas%BackgroundSpecies(iSpec)) THEN
      CALL Abort(&
          __STAMP__&
          ,'ERROR: Please supply the collision cross-section data for the particle species and NOT the background species!')
    END IF
  END DO
  IF(ANY(SpecDSMC(:)%UseCollXSec).OR.ANY(SpecDSMC(:)%UseVibXSec)) THEN
    UseMCC = .TRUE.
    CALL MCC_Init()
  ELSE
    UseMCC = .FALSE.
    XSec_NullCollision =.FALSE.
    XSec_Relaxation = .FALSE.
  END IF
  ! Ambipolar diffusion is not implemented with the regular background gas, only with MCC
  IF(DSMC%DoAmbipolarDiff) THEN
    IF((BGGas%NumberOfSpecies.GT.0).AND.(.NOT.UseMCC)) THEN
      CALL abort(__STAMP__&
          ,'ERROR: Ambipolar diffusion is not implemented with the regular background gas!')
    END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! reading/writing molecular stuff
  !-----------------------------------------------------------------------------------------------------------------------------------
  IF ((CollisMode.EQ.2).OR.(CollisMode.EQ.3)) THEN ! perform relaxation (molecular) reactions
    ! reading molecular stuff
    SpecDSMC(1:nSpecies)%Xi_Rot = 0
    SpecDSMC(1:nSpecies)%MaxVibQuant = 0
    SpecDSMC(1:nSpecies)%CharaTVib = 0
    SpecDSMC(1:nSpecies)%EZeroPoint = 0.0
    SpecDSMC(1:nSpecies)%PolyatomicMol=.false.
    SpecDSMC(1:nSpecies)%SpecToPolyArray = 0
    useRelaxProbCorrFactor=GETLOGICAL('Particles-DSMC-useRelaxProbCorrFactor','.FALSE.')
    DO iSpec = 1, nSpecies
      IF(SpecDSMC(iSpec)%InterID.NE.4) THEN
        WRITE(UNIT=hilf,FMT='(I0)') iSpec
        SpecDSMC(iSpec)%PolyatomicMol=GETLOGICAL('Part-Species'//TRIM(hilf)//'-PolyatomicMol','.FALSE.')
        IF(SpecDSMC(iSpec)%PolyatomicMol.AND.DSMC%DoTEVRRelaxation)  THEN
          CALL Abort(&
              __STAMP__&
              ,'! Simulation of Polyatomic Molecules and T-E-V-R relaxation not possible yet!!!')
        END IF
        IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
          DSMC%NumPolyatomMolecs = DSMC%NumPolyatomMolecs + 1
          SpecDSMC(iSpec)%SpecToPolyArray = DSMC%NumPolyatomMolecs
        ELSEIF ((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
          SpecDSMC(iSpec)%Xi_Rot     = 2
          SpecDSMC(iSpec)%CharaTVib  = GETREAL('Part-Species'//TRIM(hilf)//'-CharaTempVib')
          SpecDSMC(iSpec)%CharaTRot  = GETREAL('Part-Species'//TRIM(hilf)//'-CharaTempRot','0')
          SpecDSMC(iSpec)%Ediss_eV   = GETREAL('Part-Species'//TRIM(hilf)//'-Ediss_eV')
          SpecDSMC(iSpec)%MaxVibQuant = 200
          ! Calculation of the zero-point energy
          SpecDSMC(iSpec)%EZeroPoint = DSMC%GammaQuant * BoltzmannConst * SpecDSMC(iSpec)%CharaTVib
          ! Calculation of the dissociation quantum number (used for QK chemistry)
          SpecDSMC(iSpec)%DissQuant = INT(SpecDSMC(iSpec)%Ediss_eV*ElementaryCharge/(BoltzmannConst*SpecDSMC(iSpec)%CharaTVib))
        END IF
        ! Read in species values for rotational relaxation models of Boyd/Zhang if necessary
        IF(DSMC%RotRelaxProb.GT.1.0.AND.((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20))) THEN
          SpecDSMC(iSpec)%CollNumRotInf = GETREAL('Part-Species'//TRIM(hilf)//'-CollNumRotInf')
          SpecDSMC(iSpec)%TempRefRot    = GETREAL('Part-Species'//TRIM(hilf)//'-TempRefRot')
          IF(SpecDSMC(iSpec)%CollNumRotInf*SpecDSMC(iSpec)%TempRefRot.EQ.0) THEN
            CALL Abort(&
            __STAMP__&
            ,'Error! CollNumRotRef or TempRefRot is equal to zero for species:', iSpec)
          END IF
        END IF
        ! Read in species values for vibrational relaxation models of Milikan-White if necessary
        IF(DSMC%VibRelaxProb.EQ.2.0) THEN
          ! Only molecules or charged molecules
          IF(((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20))) THEN
            ALLOCATE(SpecDSMC(iSpec)%MW_ConstA(1:nSpecies))
            ALLOCATE(SpecDSMC(iSpec)%MW_ConstB(1:nSpecies))
            DO jSpec = 1, nSpecies
              WRITE(UNIT=hilf2,FMT='(I0)') jSpec
              hilf2=TRIM(hilf)//'-'//TRIM(hilf2)
              SpecDSMC(iSpec)%MW_ConstA(jSpec)     = GETREAL('Part-Species'//TRIM(hilf)//'-MWConstA-'//TRIM(hilf2))
              SpecDSMC(iSpec)%MW_ConstB(jSpec)     = GETREAL('Part-Species'//TRIM(hilf)//'-MWConstB-'//TRIM(hilf2))

              IF(SpecDSMC(iSpec)%MW_ConstA(jSpec).EQ.0) THEN
                CALL Abort(&
                __STAMP__&
                ,'Error! MW_ConstA is equal to zero for species:', iSpec)
              END IF
              IF(SpecDSMC(iSpec)%MW_ConstB(jSpec).EQ.0) THEN
                CALL Abort(&
                __STAMP__&
                ,'Error! MW_ConstB is equal to zero for species:', iSpec)
              END IF
            END DO
          END IF
          SpecDSMC(iSpec)%VibCrossSec    = GETREAL('Part-Species'//TRIM(hilf)//'-VibCrossSection')
          ! Only molecules or charged molecules
          IF((SpecDSMC(iSpec)%VibCrossSec.EQ.0).AND.((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20))) THEN
            CALL Abort(&
            __STAMP__&
            ,'Error! VibCrossSec is equal to zero for species:', iSpec)
          END IF
        END IF
        ! Setting the values of Rot-/Vib-RelaxProb to a fix value
        SpecDSMC(iSpec)%RotRelaxProb  = DSMC%RotRelaxProb
        SpecDSMC(iSpec)%VibRelaxProb  = DSMC%VibRelaxProb     ! 0.004
        SpecDSMC(iSpec)%ElecRelaxProb = DSMC%ElecRelaxProb    ! or 0.02 | Bird: somewhere in range 0.01 .. 0.02
        ! multi init stuff
        ALLOCATE(SpecDSMC(iSpec)%Init(0:Species(iSpec)%NumberOfInits))
        DO iInit = 1, Species(iSpec)%NumberOfInits
          WRITE(UNIT=hilf2,FMT='(I0)') iInit
          hilf2=TRIM(hilf)//'-Init'//TRIM(hilf2)
          IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
            SpecDSMC(iSpec)%Init(iInit)%TVib      = GETREAL('Part-Species'//TRIM(hilf2)//'-TempVib')
            SpecDSMC(iSpec)%Init(iInit)%TRot      = GETREAL('Part-Species'//TRIM(hilf2)//'-TempRot')
          END IF
          ! read electronic temperature
          IF (DSMC%ElectronicModel.GT.0) THEN
            SpecDSMC(iSpec)%Init(iInit)%Telec   = GETREAL('Part-Species'//TRIM(hilf2)//'-TempElec')
          END IF ! electronic model
        END DO !Inits
        ALLOCATE(SpecDSMC(iSpec)%Surfaceflux(1:Species(iSpec)%nSurfacefluxBCs))
        DO iInit = 1, Species(iSpec)%nSurfacefluxBCs
          WRITE(UNIT=hilf2,FMT='(I0)') iInit
          hilf2=TRIM(hilf)//'-Surfaceflux'//TRIM(hilf2)
          IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
            SpecDSMC(iSpec)%Surfaceflux(iInit)%TVib      = GETREAL('Part-Species'//TRIM(hilf2)//'-TempVib','0.')
            SpecDSMC(iSpec)%Surfaceflux(iInit)%TRot      = GETREAL('Part-Species'//TRIM(hilf2)//'-TempRot','0.')
            IF (SpecDSMC(iSpec)%Surfaceflux(iInit)%TRot*SpecDSMC(iSpec)%Surfaceflux(iInit)%TVib.EQ.0.) THEN
              CALL Abort(&
                  __STAMP__&
                  ,'Error! TVib and TRot not def. in Part-SpeciesXX-SurfacefluxXX-TempVib/TempRot for iSpec, iInit',iSpec,REAL(iInit))
            END IF
          END IF
          ! read electronic temperature
          IF (DSMC%ElectronicModel.GT.0) THEN
            SpecDSMC(iSpec)%Surfaceflux(iInit)%Telec   = GETREAL('Part-Species'//TRIM(hilf2)//'-TempElec','0.')
            IF (SpecDSMC(iSpec)%Surfaceflux(iInit)%Telec.EQ.0.) THEN
              CALL Abort(&
                  __STAMP__&
                  ,' Error! Telec not defined in Part-SpeciesXX-SurfacefluxXX-Tempelec for iSpec, iInit',iSpec,REAL(iInit))
            END IF
          END IF
        END DO !SurfaceFluxBCs
      END IF ! not electron
    END DO !Species

    ! Initialization of polyatomic species and burn-in phase (Metropolis-Hastings) per initialization region
    IF(DSMC%NumPolyatomMolecs.GT.0) THEN
      DSMC%PolySingleMode = GETLOGICAL('Particles-DSMC-PolyRelaxSingleMode','.FALSE.')
      IF((SelectionProc.NE.2).AND.DSMC%PolySingleMode) THEN
        ! Single-mode relaxation of vibrational modes of polyatomic molecules only possible when the prohibiting double relaxation
        ! method is used (Gimelshein, SelectionProc = 2)
        CALL Abort(&
            __STAMP__&
            ,'ERROR: No single-mode polyatomic relaxation possible with chosen selection procedure! SelectionProc:', SelectionProc)
      END IF
      IF(.NOT.ALLOCATED(VibQuantsPar)) ALLOCATE(VibQuantsPar(PDM%maxParticleNumber))
      ALLOCATE(PolyatomMolDSMC(DSMC%NumPolyatomMolecs))
      DO iSpec = 1, nSpecies
        IF (SpecDSMC(iSpec)%PolyatomicMol) THEN
          CALL InitPolyAtomicMolecs(iSpec)
        END IF
      END DO
    END IF
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! Comparison with Landau-Teller equation, including a different selection procedure (restricts relaxation to a single mode)
    ! and a correctional factor, both in dsmc_collis_mode.f90, also the translational temperature is fixed, in timedisc.f90
    !-----------------------------------------------------------------------------------------------------------------------------------
#if (PP_TimeDiscMethod==42)
    DSMC%CompareLandauTeller = GETLOGICAL('Particles-DSMC-CompareLandauTeller','.FALSE.')
    IF(DSMC%CompareLandauTeller) THEN
      IF(CollisMode.NE.2) THEN
        CALL abort(&
            __STAMP__&
            ,'ERROR: Comparison with Landau-Teller only available in CollisMode = 2, CollisMode:', CollisMode)
      END IF
      IF(nSpecies.GT.1) THEN
        CALL abort(&
            __STAMP__&
            ,'ERROR: Comparison with Landau-Teller only available for a single species, nSpecies:', nSpecies)
      END IF
    END IF
#endif
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! Setting the internal energy value of every particle
    DO iPart = 1, PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        IF (Species(PartSpecies(iPart))%NumberOfInits.GT.0) THEN
          iInit = PDM%PartInit(iPart)
          IF (SpecDSMC(PartSpecies(iPart))%PolyatomicMol) THEN
            CALL DSMC_SetInternalEnr_Poly(PartSpecies(iPart),iInit,iPart,1)
          ELSE
            CALL DSMC_SetInternalEnr_LauxVFD(PartSpecies(iPart),iInit,iPart,1)
          END IF
        END IF
      END IF
    END DO
    ! Array not required anymore after the initialization is completed
    DEALLOCATE(PDM%PartInit)
  END IF ! CollisMode .EQ. 2 or 3
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! Define chemical reactions (including ionization and backward reaction rate)
  !-----------------------------------------------------------------------------------------------------------------------------------
  IF (CollisMode.EQ.3) THEN ! perform chemical reactions
    DO iSpec = 1, nSpecies
      WRITE(UNIT=hilf,FMT='(I0)') iSpec
      ! Read-in of heat of formation, ions are treated later using the heat of formation of their ground state and data from the
      ! from the electronic state database to ensure consistent energies across chemical reactions of QK and Arrhenius type.
      IF((SpecDSMC(iSpec)%InterID.EQ.10).OR.(SpecDSMC(iSpec)%InterID.EQ.20).OR.(SpecDSMC(iSpec)%InterID.EQ.4)) THEN
        SpecDSMC(iSpec)%HeatOfFormation = 0.0
      ELSE
        SpecDSMC(iSpec)%HeatOfFormation = GETREAL('Part-Species'//TRIM(hilf)//'-HeatOfFormation_K')
        SpecDSMC(iSpec)%HeatOfFormation = SpecDSMC(iSpec)%HeatOfFormation * BoltzmannConst
      END IF
      ! Heat of formation of ionized species is modified with the ionization energy directly from read-in electronic energy levels
      ! of the ground/previous state of the respective species (Input requires a species number (eg species number of N for NIon1))
      IF((SpecDSMC(iSpec)%InterID.EQ.10).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
        SpecDSMC(iSpec)%PreviousState = GETINT('Part-Species'//TRIM(hilf)//'-PreviousState')
      ELSE
        SpecDSMC(iSpec)%PreviousState = 0
      END IF
      ! Read-in of species for field ionization (only required if it cannot be determined automatically)
      IF(SpecDSMC(iSpec)%InterID.NE.4) THEN
        SpecDSMC(iSpec)%NextIonizationSpecies = GETINT('Part-Species'//TRIM(hilf)//'-NextIonizationSpecies')
      ELSE
        SpecDSMC(iSpec)%NextIonizationSpecies = 0
        DSMC%ElectronSpecies = iSpec
      END IF
      ! Read-in of electronic levels for QK and backward reaction rate -------------------------------------------------------------
      IF (DSMC%ElectronicModelDatabase .EQ.'none') THEN
        IF ((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
          SpecDSMC(iSpec)%MaxElecQuant               = GETINT('Part-Species'//TRIM(hilf)//'-NumElectronicLevels','0')
          IF(SpecDSMC(iSpec)%MaxElecQuant.GT.0) THEN
            ALLOCATE(SpecDSMC(iSpec)%ElectronicState(2,0:SpecDSMC(iSpec)%MaxElecQuant-1))
            DO iDOF=1, SpecDSMC(iSpec)%MaxElecQuant
              WRITE(UNIT=hilf2,FMT='(I0)') iDOF
              SpecDSMC(iSpec)%ElectronicState(1,iDOF-1) &
                  = GETINT('Part-Species'//TRIM(hilf)//'-ElectronicDegeneracy-Level'//TRIM(hilf2),'0')
              SpecDSMC(iSpec)%ElectronicState(2,iDOF-1) &
                  = GETREAL('Part-Species'//TRIM(hilf)//'-ElectronicEnergyLevel-Level'//TRIM(hilf2),'0')
            END DO
          END IF
        END IF
      END IF
    END DO

    ! Calculating the heat of formation for ionized species (including higher ionization levels)
    ! Requires the completed read-in of species data
    CALL CalcHeatOfFormation()

    ! Set "NextIonizationSpecies" information for field ionization from "PreviousState" info
    ! NextIonizationSpecies => SpeciesID of the next higher ionization level
    CALL SetNextIonizationSpecies()

    CALL DSMC_chemical_init()
  ELSE
    DSMC%BackwardReacRate = .FALSE.
  END IF

    ! Check whether calculation of instantaneous translational temperature is required
  IF(((CollisMode.GT.1).AND.(SelectionProc.EQ.2)).OR.DSMC%BackwardReacRate.OR.DSMC%CalcQualityFactors &
            .OR.(DSMC%VibRelaxProb.EQ.2).OR.(DSMC%ElectronicModel.EQ.2)) THEN
    ! 1. Case: Inelastic collisions and chemical reactions with the Gimelshein relaxation procedure and variable vibrational
    !           relaxation probability (CalcGammaVib)
    ! 2. Case: Backward reaction rates
    ! 3. Case: Temperature required for the mean free path with the VHS model
    ALLOCATE(DSMC%InstantTransTemp(nSpecies+1))
    DSMC%InstantTransTemp = 0.0
    IF(DSMC%ElectronicModel.EQ.2) THEN
      ALLOCATE(DSMC%InstantTXiElec(2,nSpecies))
      DSMC%InstantTXiElec = 0.0
    END IF
  END IF

  !-----------------------------------------------------------------------------------------------------------------------------------
  ! Octree + Nearest Neighbour: octree cell splitting is performed until the mean free path is resolved or a certain number of
  ! particles is reached. Afterwards a nearest neighbour search for the collision partner is performed.
  ! Source: Pfeiffer, M., Mirza, A. and Fasoulas, S. (2013). A grid-independent particle pairing strategy for DSMC.
  ! Journal of Computational Physics 246, 28–36. doi:10.1016/j.jcp.2013.03.018
  !-----------------------------------------------------------------------------------------------------------------------------------
  DSMC%UseOctree = GETLOGICAL('Particles-DSMC-UseOctree')
  DSMC%UseNearestNeighbour = GETLOGICAL('Particles-DSMC-UseNearestNeighbour')
  IF(DSMC%UseOctree) THEN
    IF(NGeo.GT.PP_N) CALL abort(__STAMP__,' Set PP_N to NGeo, else, the volume is not computed correctly.')
    CALL DSMC_init_octree()
  END IF
  IF(Symmetry%Order.LE.2) THEN
    CollInf%ProhibitDoubleColl = GETLOGICAL('Particles-DSMC-ProhibitDoubleCollisions','.TRUE.')
    IF (CollInf%ProhibitDoubleColl) THEN
      IF(.NOT.ALLOCATED(CollInf%OldCollPartner)) ALLOCATE(CollInf%OldCollPartner(1:PDM%maxParticleNumber))
      CollInf%OldCollPartner = 0
    END IF
  ELSE
    CollInf%ProhibitDoubleColl = GETLOGICAL('Particles-DSMC-ProhibitDoubleCollisions','.FALSE.')
    IF (CollInf%ProhibitDoubleColl) CALL abort(__STAMP__,&
      'ERROR: Prohibiting double collisions is only supported within a 2D/axisymmetric simulation!')
  END IF

  !-----------------------------------------------------------------------------------------------------------------------------------
  ! Background gas: Check compatibility with other features
  !-----------------------------------------------------------------------------------------------------------------------------------
  IF (BGGas%NumberOfSpecies.GT.0) THEN
    IF (DSMC%UseOctree) THEN
      CALL abort(__STAMP__,&
          'ERROR: Utilization of the octree and nearest neighbour scheme not possible with the background gas!')
    END IF
    DO iSpec = 1, nSpecies
      IF(BGGas%BackgroundSpecies(iSpec)) THEN
        IF(SpecDSMC(iSpec)%InterID.EQ.4) THEN
          CALL abort(__STAMP__,&
            'ERROR in BGGas: Electrons as background gas are not yet available!')
        END IF
      END IF
    END DO
  END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! Calculate vib collision numbers and characteristic velocity, according to Abe
!-----------------------------------------------------------------------------------------------------------------------------------
  ! (i) dref changed from dref = 0.5 * (dref_1+dref_2) 
  !                  to   dref(iSpec,jSpec) which is identical to old definition (for averagedCollisionParameters=TRUE (DEFAULT))
  ! in case of averagedCollisionParameter=FALSE dref(iSpec,jSpec) contains collision specific dref see --help for details
  IF((DSMC%VibRelaxProb.EQ.2).AND.(CollisMode.GE.2)) THEN
    VarVibRelaxProb%alpha = GETREAL('Particles-DSMC-alpha','0.99')
    IF ((VarVibRelaxProb%alpha.LT.0).OR.(VarVibRelaxProb%alpha.GE.1)) THEN
      CALL abort(&
      __STAMP__&
      ,'ERROR: Particles-DSMC-alpha has to be in the range between 0 and 1')
    END IF
    DO iSpec = 1, nSpecies
      IF(.NOT.((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20))) CYCLE
      ALLOCATE(SpecDSMC(iSpec)%CharaVelo(1:nSpecies))
      ALLOCATE(SpecDSMC(iSpec)%CollNumVib(1:nSpecies))
      DO jSpec = 1, nSpecies
        iCase = CollInf%Coll_Case(iSpec,jSpec)
        ! Calculation Of CharaVelo (g^*) according to Abe (doi:10.1063/1.868094)
        SpecDSMC(iSpec)%CharaVelo(jSpec) = SQRT(  BoltzmannConst / CollInf%MassRed(iCase) &
                                         * (2./3.*SpecDSMC(iSpec)%MW_ConstA(jSpec))**3.)
        ! Calculation Of CollNumVib (Z_0) according to Abe (doi:10.1063/1.868094)
        SpecDSMC(iSpec)%CollNumVib(jSpec) = (2.* (SpecDSMC(iSpec)%CharaVelo(jSpec)**2) * EXP(SpecDSMC(iSpec)%MW_ConstB(jSpec)) &
                                  / (SQRT(3.) * CollInf%MassRed(iCase)) ) &
                                  * pi/4.*(CollInf%dref(iSpec,jSpec))**2 &
                                  * ( 2.* (2.-CollInf%omega(iSpec,iSpec)) * BoltzmannConst * CollInf%Tref(iSpec,iSpec) &
                                  / CollInf%MassRed(iCase) ) ** CollInf%omega(iSpec,iSpec)
      END DO ! jSpec
      DEALLOCATE(SpecDSMC(iSpec)%MW_ConstA)
      DEALLOCATE(SpecDSMC(iSpec)%MW_ConstB)
    END DO ! iSpec
    IF(DSMC%CalcQualityFactors) THEN
      IF(nSpecies.GT.1) THEN
        ALLOCATE(DSMC%QualityFacSampVib(1:nElems,1:nSpecies+1,1:2))
        ALLOCATE(DSMC%QualityFacSampVibSamp(1:nElems,1:nSpecies+1,2))
      ELSE
        ALLOCATE(DSMC%QualityFacSampVib(1:nElems,1,1:2))
        ALLOCATE(DSMC%QualityFacSampVibSamp(1:nElems,1,2))
      END IF
      DSMC%QualityFacSampVib = 0.
      DSMC%QualityFacSampVibSamp = 0
    END IF
    CALL SetVarVibProb2Elems()
    ! CHeck if DSMC%InstantTransTemp is still needed
    IF(.NOT.(((CollisMode.GT.1).AND.(SelectionProc.EQ.2)).OR.DSMC%BackwardReacRate.OR.DSMC%CalcQualityFactors)) THEN
      SDEALLOCATE(DSMC%InstantTransTemp)
    END IF
  END IF ! VibRelaxProb = 2

  IF((DSMC%RotRelaxProb.GE.2).AND.DSMC%CalcQualityFactors.AND.(CollisMode.GE.2)) THEN
    IF(nSpecies.GT.1) THEN
      ALLOCATE(DSMC%QualityFacSampRot(1:nElems,1:nSpecies+1,1:2))
      ALLOCATE(DSMC%QualityFacSampRotSamp(1:nElems,1:nSpecies+1))
    ELSE
      ALLOCATE(DSMC%QualityFacSampRot(1:nElems,1,1:2))
      ALLOCATE(DSMC%QualityFacSampRotSamp(1:nElems,1))
    END IF
    ALLOCATE(DSMC%CalcRotProb(1:nSpecies,1:3))
    DSMC%QualityFacSampRot = 0.
    DSMC%QualityFacSampRotSamp = 0
    DSMC%CalcRotProb = 0.
  END IF
  IF(DSMC%CalcQualityFactors) THEN
    ALLOCATE(DSMC%CalcVibProb(1:nSpecies,1:3))
    DSMC%CalcVibProb = 0.
  END IF
  IF(XSec_Relaxation) SpecXSec(:)%VibCount = 0.
END IF !CollisMode.GT.0

! If field ionization is used without chemical reactions due to collisions (DSMC chemistry)
IF(DoFieldIonization.AND.(CollisMode.NE.3))THEN
  DO iSpec = 1, nSpecies
    WRITE(UNIT=hilf,FMT='(I0)') iSpec
    ! Heat of formation of ionized species is modified with the ionization energy directly from read-in electronic energy levels
    ! of the ground/previous state of the respective species (Input requires a species number (eg species number of N for NIon1))
    IF((SpecDSMC(iSpec)%InterID.EQ.10).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
      SpecDSMC(iSpec)%PreviousState = GETINT('Part-Species'//TRIM(hilf)//'-PreviousState')
    ELSE
      SpecDSMC(iSpec)%PreviousState = 0
    END IF ! (SpecDSMC(iSpec)%InterID.EQ.10).OR.(SpecDSMC(iSpec)%InterID.EQ.20)

    ! Read-in of species for field ionization (only required if it cannot be determined automatically)
    IF ((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
      SpecDSMC(iSpec)%NextIonizationSpecies = GETINT('Part-Species'//TRIM(hilf)//'-NextIonizationSpecies')
    ELSE
      SpecDSMC(iSpec)%NextIonizationSpecies = 0
    END IF
  END DO ! iSpec = 1, nSpecies

  ! Set "NextIonizationSpecies" information for field ionization from "PreviousState" info
  ! NextIonizationSpecies => SpeciesID of the next higher ionization level
  CALL SetNextIonizationSpecies()
END IF

SWRITE(UNIT_stdOut,'(A)')' INIT DSMC DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitDSMC


SUBROUTINE SetElectronicModel(iSpec)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals              ,ONLY: abort
USE MOD_DSMC_Vars            ,ONLY: SpecDSMC
USE MOD_DSMC_ElectronicModel ,ONLY: ReadSpeciesLevel
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: iSpec
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(SpecDSMC(iSpec)%Name.EQ.'none') THEN
  CALL Abort(&
      __STAMP__,&
      "Read-in from electronic database requires the definition of species name! Species:",iSpec)
END IF
IF(.NOT.SpecDSMC(iSpec)%FullyIonized)THEN
  CALL ReadSpeciesLevel(SpecDSMC(iSpec)%Name,iSpec)
END IF
END SUBROUTINE SetElectronicModel


SUBROUTINE CalcHeatOfFormation()
!===================================================================================================================================
! Calculating the heat of formation for ionized species (including higher ionization levels)
! Requires the completed read-in of species data
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_ReadInTools
USE MOD_Globals       ,ONLY: abort,UNIT_stdOut
#if USE_MPI
USE MOD_Globals       ,ONLY: mpiroot
#endif
USE MOD_Globals_Vars  ,ONLY: BoltzmannConst
USE MOD_PARTICLE_Vars ,ONLY: nSpecies
USE MOD_DSMC_Vars     ,ONLY: SpecDSMC
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(32) :: hilf2
LOGICAL       :: AutoDetect
INTEGER       :: iSpec,jSpec,counter,MaxElecQua
!===================================================================================================================================
AutoDetect=.TRUE.
DO iSpec = 1, nSpecies
  counter = 0
  IF((SpecDSMC(iSpec)%InterID.EQ.10).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    IF(SpecDSMC(iSpec)%PreviousState.EQ.0) THEN
      WRITE(UNIT=hilf2,FMT='(I0)') iSpec
      SpecDSMC(iSpec)%HeatOfFormation = GETREAL('Part-Species'//TRIM(hilf2)//'-HeatOfFormation_K') * BoltzmannConst
    ELSE
      IF(SpecDSMC(SpecDSMC(iSpec)%PreviousState)%MaxElecQuant.GT.0) THEN
        jSpec = SpecDSMC(iSpec)%PreviousState
        DO
          MaxElecQua = SpecDSMC(jSpec)%MaxElecQuant - 1
          SpecDSMC(iSpec)%HeatOfFormation = SpecDSMC(iSpec)%HeatOfFormation &
              + SpecDSMC(jSpec)%ElectronicState(2,MaxElecQua)*BoltzmannConst
          IF(SpecDSMC(jSpec)%PreviousState.EQ.0) EXIT
          jSpec = SpecDSMC(jSpec)%PreviousState
          ! Fail-safe, abort after 100 iterations
          counter = counter + 1
          IF(counter.GT.100) THEN
            CALL abort(&
                __STAMP__&
                ,'ERROR: Nbr. of ionization lvls per spec limited to 100. More likely wrong input in PreviuosState of spec:', iSpec)
          END IF
        END DO
        IF(AutoDetect)THEN
          SWRITE(UNIT_stdOut,'(A)')' Automatically determined HeatOfFormation:'
          AutoDetect=.FALSE.
        END IF
        ! Add the heat of formation of the ground state
        SpecDSMC(iSpec)%HeatOfFormation = SpecDSMC(iSpec)%HeatOfFormation + SpecDSMC(jSpec)%HeatOfFormation
        WRITE(UNIT=hilf2,FMT='(I0)') iSpec
        CALL PrintOption('part-species'//TRIM(hilf2)//'-heatofformation_k','CALCUL.',&
            RealOpt=SpecDSMC(iSpec)%HeatOfFormation/BoltzmannConst)
      ELSE
        CALL abort(&
            __STAMP__&
            ,'ERROR: Chemical reactions with ionized species require an input of electronic energy level(s)!', iSpec)
      END IF
    END IF
  END IF
END DO
END SUBROUTINE CalcHeatOfFormation


SUBROUTINE SetNextIonizationSpecies()
!===================================================================================================================================
! Set "NextIonizationSpecies" information for field ionization from "PreviousState" info
! NextIonizationSpecies => SpeciesID of the next higher ionization level
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals       ,ONLY: UNIT_stdOut
#if USE_MPI
USE MOD_Globals       ,ONLY: mpiroot
#endif
USE MOD_PARTICLE_Vars ,ONLY: nSpecies
USE MOD_DSMC_Vars     ,ONLY: SpecDSMC
USE MOD_ReadInTools   ,ONLY: PrintOption
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(32) :: hilf2
LOGICAL       :: AutoDetect
INTEGER       :: iSpec
!===================================================================================================================================
AutoDetect=.FALSE.
DO iSpec = 1, nSpecies
  ! loop all species, except electrons (also loop over fully ionized species!)
  IF(SpecDSMC(iSpec)%InterID.NE.4) THEN
    IF(SpecDSMC(iSpec)%PreviousState.NE.0)THEN
      SpecDSMC(SpecDSMC(iSpec)%PreviousState)%NextIonizationSpecies = iSpec
      AutoDetect=.TRUE.
    END IF
  ELSE
    SpecDSMC(iSpec)%NextIonizationSpecies = 0
  END IF
END DO
IF(AutoDetect)THEN
  SWRITE(UNIT_stdOut,'(A)')' Automatically determined NextIonizationSpecies:'
  DO iSpec = 1, nSpecies
    WRITE(UNIT=hilf2,FMT='(I0)') iSpec
    CALL PrintOption('iSpec='//TRIM(hilf2)//': NextIonizationSpecies','CALCUL.',IntOpt=SpecDSMC(iSpec)%NextIonizationSpecies)
  END DO
END IF
END SUBROUTINE SetNextIonizationSpecies


SUBROUTINE SetVarVibProb2Elems()
!===================================================================================================================================
! Set initial vibrational relaxation probability to all elements
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals                ,ONLY: abort, IK, MPI_COMM_WORLD
USE MOD_PARTICLE_Vars          ,ONLY: nSpecies, Species
USE MOD_Restart_Vars           ,ONLY: DoRestart,RestartFile
USE MOD_Particle_Vars          ,ONLY: nSpecies, PartSpecies
USE MOD_part_tools             ,ONLY: GetParticleWeight
USE MOD_HDF5_INPUT             ,ONLY: DatasetExists,ReadAttribute,ReadArray
USE MOD_IO_HDF5
USE MOD_DSMC_Vars              ,ONLY: VarVibRelaxProb, CollInf, SpecDSMC, Coll_pData, DSMC
USE MOD_Mesh_Vars              ,ONLY: nElems, offsetElem
USE MOD_DSMC_Analyze           ,ONLY: CalcInstantTransTemp
USE MOD_Particle_Vars          ,ONLY: PEM
USE MOD_DSMC_Relaxation        ,ONLY: DSMC_calc_var_P_vib
USE MOD_part_emission_tools    ,ONLY: CalcVelocity_maxwell_lpn
#if USE_MPI
USE MOD_Globals                ,ONLY: MPIRoot
#endif /*USE_MPI*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: iSpec, jSpec, iPart, iElem, dim, n
  REAL                  :: VibProb, Ti, Tj, CRela2, Velo1(3), Velo2(3)
  INTEGER               :: nPart, iLoop, nLoop
  INTEGER, ALLOCATABLE  :: iPartIndx(:), nPerSpec(:)
  LOGICAL               :: VibProbInitDone, VibProbDataExists

!===================================================================================================================================
  ALLOCATE(VarVibRelaxProb%ProbVibAv(1:nElems,1:nSpecies))
  ALLOCATE(VarVibRelaxProb%ProbVibAvNew(1:nSpecies))
  ALLOCATE(VarVibRelaxProb%nCollis(1:nSpecies))
  VarVibRelaxProb%ProbVibAv = 0
  VibProbInitDone = .FALSE.
  IF (DoRestart) THEN
    CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
    ! read local ParticleInfo from HDF5
    CALL DatasetExists(File_ID,'VibProbInfo',VibProbDataExists)
    IF(VibProbDataExists)THEN
      VibProbInitDone = .TRUE.
      ! Associate construct for integer KIND=8 possibility
      SWRITE(*,*) 'Set variable vibrational relaxation probability from restart file'
      ASSOCIATE (&
            nSpecies   => INT(nSpecies,IK) ,&
            offsetElem => INT(offsetElem,IK),&
            nElems     => INT(nElems,IK)    )
        CALL ReadArray('VibProbInfo',2,(/nElems, nSpecies/),offsetElem,1,RealArray=VarVibRelaxProb%ProbVibAv(:,:))
      END ASSOCIATE
    END IF ! If 'VibProbInfo' exists
    CALL DatasetExists(File_ID,'VibProbConstInfo',VibProbDataExists,attrib=.TRUE.)
    IF((.NOT.VibProbInitDone).AND.VibProbDataExists) THEN
      VibProbInitDone = .TRUE.
      CALL ReadAttribute(File_ID,'VibProbConstInfo',1,RealScalar=VibProb)
      ! Set vibrational relaxation probability to former value
      SWRITE(*,*) 'Set uniform vibrational relaxation probability from restart file'
      DO iElem = 1, nElems
        DO iSpec = 1, nSpecies
          VarVibRelaxProb%ProbVibAv(iElem,iSpec) = VibProb
        END DO
      END DO
    END IF ! If 'VibProbConstInfo' exists
    IF(.NOT.VibProbInitDone) THEN
      ! Set vibrational relaxation probability to default value
      SWRITE(*,*) 'No vibrational relaxation probability data in restart file\n', &
                  'Set uniform vibrational relaxation probability of', 0.004
      DO iElem = 1, nElems
        DO iSpec = 1, nSpecies
          VarVibRelaxProb%ProbVibAv(iElem,iSpec) = 0.004
        END DO
      END DO
    END IF ! No restart information exist
    CALL CloseDataFile()
  ELSE ! If not DoRestart
    ALLOCATE(Coll_pData(1))
    ALLOCATE(nPerSpec(nSpecies))
    SWRITE(*,*) 'Set vibrational relaxation probability based on temperature in the cell'
    DO iElem = 1, nElems
      nPerSpec = 0
      nPart = PEM%pNumber(iElem)
      ! List of particles in the cell neccessary
      ALLOCATE(iPartIndx(nPart))
      iPartIndx(1:nPart) = 0
      ! create particle index list
      iPart = PEM%pStart(iElem)
      DO iLoop = 1, nPart
        iPartIndx(iLoop) = iPart
        iPart = PEM%pNext(iPart)
      END DO
      CollInf%Coll_SpecPartNum = 0
      DO iPart = 1, nPart
        CollInf%Coll_SpecPartNum(PartSpecies(iPartIndx(iPart))) = &
                  CollInf%Coll_SpecPartNum(PartSpecies(iPartIndx(iPart))) + GetParticleWeight(iPartIndx(iPart))
        nPerSpec(PartSpecies(iPartIndx(iPart))) = nPerSpec(PartSpecies(iPartIndx(iPart))) + 1
      END DO
      CALL CalcInstantTransTemp(iPartIndx,nPart)
      DO iSpec = 1, nSpecies
        IF(.NOT.((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)))CYCLE
        IF((DSMC%InstantTransTemp(iSpec).NE.0).AND.(nPerSpec(iSpec).GE.5)) THEN
          Ti = DSMC%InstantTransTemp(iSpec)
        ELSE
          IF(DSMC%InstantTransTemp(nSpecies + 1).NE.0) THEN
            Ti = DSMC%InstantTransTemp(nSpecies + 1)
          ELSE
            Ti = Species(iSpec)%Init(1)%MWTemperatureIC
          END IF
        END IF
        n = 1
        DO jSpec = 1, nSpecies
          IF((DSMC%InstantTransTemp(jSpec).NE.0).AND.(nPerSpec(jSpec).GE.5)) THEN
            Tj = DSMC%InstantTransTemp(jSpec)
          ELSE
            IF(DSMC%InstantTransTemp(nSpecies + 1).NE.0) THEN
              Tj = DSMC%InstantTransTemp(nSpecies + 1)
            ELSE
              Tj = Species(jSpec)%Init(1)%MWTemperatureIC
            END IF
          END IF
          Coll_pData(1)%PairType = CollInf%Coll_Case(iSpec, jSpec)
          ! Calculate number of samples vor each collision pair dependent to alpha and the mole fraction in the cell
          IF(nPart.GT.0) THEN
            nLoop = INT( 1. / (1.-VarVibRelaxProb%alpha) * nPerSpec(jSpec)/nPart )
          ELSE
            nLoop = INT( 1. / (1.-VarVibRelaxProb%alpha) / nSpecies )
          END IF
          DO iLoop = 1,nLoop
            ! Calculate random relative velocity
            CRela2 = 0
            CALL CalcVelocity_maxwell_lpn(iSpec, Velo1, Temperature=Ti)
            CALL CalcVelocity_maxwell_lpn(jSpec, Velo2, Temperature=Tj)
            DO dim = 1,3
              CRela2 = CRela2 + (Velo1(dim)-Velo2(dim))**2
            END DO ! dim = 3
            Coll_pData(1)%CRela2 = CRela2
            CALL DSMC_calc_var_P_vib(iSpec,jSpec,1,VibProb)
            VarVibRelaxProb%ProbVibAv(iElem,iSpec) = VarVibRelaxProb%ProbVibAv(iElem,iSpec) &
                                                   + (VibProb - VarVibRelaxProb%ProbVibAv(iElem,iSpec)) / n
            n = n + 1
          END DO ! iLoop = nLoop
        END DO ! jSpec = nSpecies
      END DO ! iSpec = nSpecies
      SDEALLOCATE(iPartIndx)
    END DO ! iElem = nElems
    SDEALLOCATE(Coll_pData)
    SDEALLOCATE(nPerSpec)
  END IF

END SUBROUTINE SetVarVibProb2Elems


SUBROUTINE FinalizeDSMC()
!----------------------------------------------------------------------------------------------------------------------------------!
! finalize dsmc variables
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_DSMC_Vars
USE MOD_Particle_Vars, ONLY:PDM
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER       :: iPoly
!===================================================================================================================================
SDEALLOCATE(VarVibRelaxProb%ProbVibAv)
SDEALLOCATE(VarVibRelaxProb%ProbVibAvNew)
SDEALLOCATE(VarVibRelaxProb%nCollis)
SDEALLOCATE(DSMC%QualityFacSampRot)
SDEALLOCATE(DSMC%QualityFacSampVib)
SDEALLOCATE(DSMC%QualityFacSampRotSamp)
SDEALLOCATE(DSMC%QualityFacSampVibSamp)
SDEALLOCATE(DSMC%CalcVibProb)
SDEALLOCATE(DSMC%CalcRotProb)
SDEALLOCATE(SampDSMC)
SDEALLOCATE(DSMC_RHS)
SDEALLOCATE(PartStateIntEn)
SDEALLOCATE(SpecDSMC)
IF(DSMC%NumPolyatomMolecs.GT.0) THEN
  DO iPoly=1,DSMC%NumPolyatomMolecs
    SDEALLOCATE(PolyatomMolDSMC(iPoly)%CharaTVibDOF)
    SDEALLOCATE(PolyatomMolDSMC(iPoly)%LastVibQuantNums)
    SDEALLOCATE(PolyatomMolDSMC(iPoly)%MaxVibQuantDOF)
    SDEALLOCATE(PolyatomMolDSMC(iPoly)%GammaVib)
    SDEALLOCATE(PolyatomMolDSMC(iPoly)%VibRelaxProb)
    SDEALLOCATE(PolyatomMolDSMC(iPoly)%CharaTRotDOF)
  END DO
  SDEALLOCATE(PolyatomMolDSMC)
END IF
SDEALLOCATE(DSMC%NumColl)
SDEALLOCATE(DSMC%InstantTransTemp)
IF(DSMC%CalcQualityFactors) THEN
  SDEALLOCATE(DSMC%QualityFacSamp)
END IF
SDEALLOCATE(PDM%PartInit)
SDEALLOCATE(Coll_pData)
SDEALLOCATE(SampDSMC)
SDEALLOCATE(MacroDSMC)
SDEALLOCATE(QKChemistry)

SDEALLOCATE(ChemReac%ReactModel)
SDEALLOCATE(ChemReac%BackwardReac)
SDEALLOCATE(ChemReac%QKRColl)
SDEALLOCATE(ChemReac%QKTCollCorrFac)
SDEALLOCATE(ChemReac%NumReac)
SDEALLOCATE(ChemReac%ReacCount)
SDEALLOCATE(ChemReac%ReacCollMean)
SDEALLOCATE(ChemReac%NumReac)
SDEALLOCATE(ChemReac%ReactType)
SDEALLOCATE(ChemReac%Reactants)
SDEALLOCATE(ChemReac%Products)
SDEALLOCATE(ChemReac%ReactCase)
SDEALLOCATE(ChemReac%ReactNum)
SDEALLOCATE(ChemReac%Arrhenius_Prefactor)
SDEALLOCATE(ChemReac%Arrhenius_Powerfactor)
SDEALLOCATE(ChemReac%EActiv)
SDEALLOCATE(ChemReac%EForm)
SDEALLOCATE(ChemReac%MeanEVib_PerIter)
SDEALLOCATE(ChemReac%MeanEVibQua_PerIter)
SDEALLOCATE(ChemReac%MeanXiVib_PerIter)
SDEALLOCATE(ChemReac%CEXa)
SDEALLOCATE(ChemReac%CEXb)
SDEALLOCATE(ChemReac%MEXa)
SDEALLOCATE(ChemReac%MEXb)
SDEALLOCATE(ChemReac%ELa)
SDEALLOCATE(ChemReac%ELb)
SDEALLOCATE(ChemReac%DoScat)
SDEALLOCATE(ChemReac%ReactInfo)
SDEALLOCATE(ChemReac%TLU_FileName)
SDEALLOCATE(ChemReac%CrossSection)
SDEALLOCATE(ChemReac%ReactNumRecomb)
SDEALLOCATE(ChemReac%Hab)
SDEALLOCATE(ChemReac%DeleteProductsList)
SDEALLOCATE(ChemReac%CollCaseInfo)

SDEALLOCATE(CollInf%collidingSpecies)
SDEALLOCATE(CollInf%Coll_Case)
SDEALLOCATE(CollInf%Coll_CaseNum)
SDEALLOCATE(CollInf%Coll_SpecPartNum)
SDEALLOCATE(CollInf%Cab)
SDEALLOCATE(CollInf%KronDelta)
SDEALLOCATE(CollInf%FracMassCent)
SDEALLOCATE(CollInf%MassRed)
SDEALLOCATE(CollInf%MeanMPF)
SDEALLOCATE(CollInf%alphaVSS)
SDEALLOCATE(CollInf%omega)
SDEALLOCATE(CollInf%dref)
SDEALLOCATE(CollInf%Tref)
!SDEALLOCATE(SampWall)
SDEALLOCATE(MacroSurfaceVal)
!SDEALLOCATE(VibQuantsPar)
! SDEALLOCATE(XiEq_Surf)
SDEALLOCATE(DSMC_Solution)
CALL DeleteElemNodeVol()
SDEALLOCATE(BGGas%PairingPartner)
SDEALLOCATE(BGGas%BackgroundSpecies)
SDEALLOCATE(BGGas%MapSpecToBGSpec)
SDEALLOCATE(BGGas%MapBGSpecToSpec)
SDEALLOCATE(BGGas%SpeciesFraction)
SDEALLOCATE(BGGas%NumberDensity)
SDEALLOCATE(RadialWeighting%ClonePartNum)
SDEALLOCATE(ClonedParticles)
SDEALLOCATE(SymmetrySide)
SDEALLOCATE(SpecXSec)
END SUBROUTINE FinalizeDSMC


SUBROUTINE DeleteElemNodeVol()
!----------------------------------------------------------------------------------------------------------------------------------!
! Delete the pointer tree ElemNodeVol
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_DSMC_Vars
USE MOD_Mesh_Vars              ,ONLY: nElems
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iElem
!===================================================================================================================================
IF(.NOT.DSMC%UseOctree) RETURN
DO iElem=1,nElems
  CALL DeleteNodeVolume(ElemNodeVol(iElem)%Root)
  DEALLOCATE(ElemNodeVol(iElem)%Root)
END DO
DEALLOCATE(ElemNodeVol)
END SUBROUTINE DeleteElemNodeVol


RECURSIVE SUBROUTINE DeleteNodeVolume(Node)
!----------------------------------------------------------------------------------------------------------------------------------!
! Check if the Node has subnodes and delete them
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_DSMC_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
TYPE (tNodeVolume), INTENT(IN), POINTER  :: Node
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(ASSOCIATED(Node%SubNode1)) THEN 
  CALL DeleteNodeVolume(Node%SubNode1)
  DEALLOCATE(Node%SubNode1)
END IF
IF(ASSOCIATED(Node%SubNode2)) THEN
  CALL DeleteNodeVolume(Node%SubNode2)
  DEALLOCATE(Node%SubNode2)
END IF
IF(ASSOCIATED(Node%SubNode3)) THEN
  CALL DeleteNodeVolume(Node%SubNode3)
  DEALLOCATE(Node%SubNode3)
END IF
IF(ASSOCIATED(Node%SubNode4)) THEN
  CALL DeleteNodeVolume(Node%SubNode4)
  DEALLOCATE(Node%SubNode4)
END IF
IF(ASSOCIATED(Node%SubNode5)) THEN
  CALL DeleteNodeVolume(Node%SubNode5)
  DEALLOCATE(Node%SubNode5)
END IF
IF(ASSOCIATED(Node%SubNode6)) THEN
  CALL DeleteNodeVolume(Node%SubNode6)
  DEALLOCATE(Node%SubNode6)
END IF
IF(ASSOCIATED(Node%SubNode7)) THEN
  CALL DeleteNodeVolume(Node%SubNode7)
  DEALLOCATE(Node%SubNode7)
END IF
IF(ASSOCIATED(Node%SubNode8)) THEN
  CALL DeleteNodeVolume(Node%SubNode8)
  DEALLOCATE(Node%SubNode8)
END IF
END SUBROUTINE DeleteNodeVolume


END MODULE MOD_DSMC_Init
