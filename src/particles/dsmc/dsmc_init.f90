!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (piclas.boltzplatz.eu/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
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

PUBLIC :: DefineParametersDSMC, InitDSMC, FinalizeDSMC
PUBLIC :: SetVarVibProb2Elems
!===================================================================================================================================
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
                                            '0: No collisions (=free molecular flow with DSMC-Sampling-Routines).\n'//&
                                            '1: Elastic collisions \n'//&
                                            '2: Relaxation + elastic collisions \n'//&
                                            '3: Mode 2 + chemical reactions.', '1')
CALL prms%CreateIntOption(      'Particles-DSMC-SelectionProcedure'     &
                                        , 'Selection procedure for relaxation process\n'//&
                                          '1: Laux, multiple relaxation processes possible\n'//&
                                          '2: Gimelshein, prohibiting double/multiple relaxation', '1')
CALL prms%CreateRealOption(     'Particles-DSMC-RotRelaxProb'&
                                          , 'Define the rotational relaxation probability upon collision of molecules\n'//&
                                          '0-1: constant probability\n'//&
                                          '2: variable, Boyd)', '0.2')
CALL prms%CreateIntOption(  'Particles-DSMC-RotationalRelaxModel' &
                                        , 'Define model of rotational relaxation handling in DSMC \n'//&
                                          '0: Continuous treatment \n'//&
                                          '1: Analytic model + rotational energy levels from unified'//&
                                              'species database only for asymmetric top molecules \n'//&
                                          '2: rotational energy levels for all species from unified species database','0')
CALL prms%CreateIntOption(  'Particles-DSMC-VibrationalRelaxModel' &
                                        , 'Define model of vibrational relaxation handling in DSMC \n'//&
                                          '0: Simple harmonic oscillator (SHO) \n'//&
                                          '1: Anharmonic oscillator (AHO, Morse potential)','0')
CALL prms%CreateRealOption(     'Particles-DSMC-VibRelaxProb'&
                                          , 'Define the vibrational relaxation probability upon collision of molecules', '0.004')
CALL prms%CreateRealOption(     'Part-Species[$]-ElecRelaxProb'  &
                                           ,'Define the electronic relaxation probability per species','0.01',numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Particles-DSMC-Vib-Anharmonic-Species[$]-Data'&
                                          , 'If vibrational AHO model is used give (relative) path to .csv of database for each'//&
                                          'species', 'none', numberedmulti=.TRUE.)
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
                                          'Set [TRUE] to disable particle movement. Use for reservoir simulations.' , '.FALSE.')
CALL prms%CreateLogicalOption(  'Particles-DSMCReservoirSimRate', &
                                          'Only with Particles-DSMCReservoirSim = T\n'//&
                                          'Set [TRUE] to disable particle reactions. Only probabilities (rates) are calculated.', &
                                          '.FALSE.')
CALL prms%CreateLogicalOption(  'Particles-DSMCReservoirStatistic', &
                                          'Only with Particles-DSMCReservoirSim = T\n'//&
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
                                          '2: Model by Burt, each particle has an electronic distribution function\n'//&
                                          '3: MCC model, utilizing cross-section data for specific levels\n'//&
                                          '4: Landau-Teller based model, relaxation of given distribution function', '0')
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
                                         ,'Allows the comparison with Landau-Teller equation. Only with Particles-DSMCReservoirSim = T.',&
                                          '.FALSE.')
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
CALL prms%CreateRealOption(     'Part-Species[$]-MomentOfInertia','Moment of inertia in [kg*m²]', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Ediss_eV','Energy of Dissoziation in [eV].', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Vib-Anharmonic-omegaE' &
                                           ,'Spectroscopy constant omegaE for anharmonic oscillator model, in 1/m'&
                                           , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Vib-Anharmonic-chiE' &
                                           ,'Spectroscopy constant chiE for anharmonic oscillator model'&
                                           , numberedmulti=.TRUE.)
! ----------------------------------------------------------------------------------------------------------------------------------
CALL prms%CreateLogicalOption(  'Particles-DSMC-useRelaxProbCorrFactor'&
                                           ,'Use the relaxation probability correction factor of Lumpkin', '.FALSE.')
CALL prms%CreateRealOption(     'Part-Species[$]-[$]-VibConstant-C1'  &
                                           ,'Constant for vibrational collision number according to Bird' &
                                           , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-[$]-VibConstant-C2'  &
                                           ,'Constant for vibrational collision number according to Bird' &
                                           , numberedmulti=.TRUE.)
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
CALL prms%CreateRealOption(     'Part-Species[$]-MomentOfInertia[$]'  &
                                            ,'MomentOfInertia [kg*m^2], linear molecules require only a single input'//&
                                            'while non-linear molecules require three''0.', numberedmulti=.TRUE.)
! Granular species
CALL prms%CreateRealOption(     'Part-Species[$]-GranularPartCsp' &
                                            ,'Solid particle specific heat [J/(kg*K)].'&
                                            ,'0.0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-GranularPartTau' &
                                            ,'Thermal accommodation coefficient of gas species during granular particle interaction.'&
                                            ,'1.0', numberedmulti=.TRUE.)
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
USE MOD_Globals_Vars           ,ONLY: Pi, BoltzmannConst, ElementaryCharge, PlanckConst
USE MOD_Particle_Vars          ,ONLY: nSpecies, Species, PDM, UseVarTimeStep, usevMPF, UseGranularSpecies
USE MOD_Particle_Vars          ,ONLY: DoFieldIonization, SpeciesDatabase,  SampleElecExcitation
USE MOD_part_tools             ,ONLY: RotInitPolyRoutineFuncPTR, CalcERotQuant_particle, CalcERot_particle, CalcERotDataset_particle
USE MOD_Symmetry_Vars          ,ONLY: Symmetry
USE MOD_DSMC_ParticlePairing   ,ONLY: DSMC_init_octree
USE MOD_DSMC_ChemInit          ,ONLY: DSMC_chemical_init
USE MOD_DSMC_PolyAtomicModel   ,ONLY: InitPolyAtomicMolecs, DSMC_RotRelaxDatabasePoly, DSMC_RotRelaxQuantPoly, DSMC_RotRelaxPoly
USE MOD_DSMC_PolyAtomicModel   ,ONLY: RotRelaxPolyRoutineFuncPTR
USE MOD_DSMC_Relaxation        ,ONLY: DSMC_RotRelaxDiaContinuous,DSMC_RotRelaxDiaQuant, RotRelaxDiaRoutineFuncPTR
USE MOD_DSMC_CollisVec         ,ONLY: DiceDeflectedVelocityVector4Coll, DiceVelocityVector4Coll, PostCollVec
USE MOD_DSMC_BGGas             ,ONLY: BGGas_RegionsSetInternalTemp
USE MOD_io_hdf5
USE MOD_HDF5_input             ,ONLY: ReadAttribute, DatasetExists, AttributeExists
USE MOD_TimeDisc_Vars          ,ONLY: ManualTimeStep, TEnd
USE MOD_Restart_Vars           ,ONLY: DoRestart
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(32)         :: hilf , hilf2
INTEGER               :: iCase, iSpec, jSpec, iInit, iDOF, VarNum, err
INTEGER               :: iColl, jColl, pColl  ! for collision parameter read in
REAL                  :: A1, A2, delta_ij     ! species constant for cross section (p. 24 Laux)
LOGICAL               :: PostCollPointerSet
LOGICAL               :: AttrExists
CHARACTER(LEN=64)     :: dsetname
INTEGER(HID_T)        :: file_id_specdb                       ! File identifier
INTEGER               :: IntToLog
CHARACTER(LEN=255)    :: hilfname, FileNameAHO(nSpecies)
!===================================================================================================================================
LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' DSMC INIT ...'

! Initialize variables
ReactionProbGTUnityCounter = 0
DSMC%NumPolyatomMolecs = 0
SamplingActive = .FALSE.
DSMC%VibAHO = .FALSE.
DSMC%GammaQuant = 0.5

!-----------------------------------------------------------------------------------
DSMC%ReservoirSimu = GETLOGICAL('Particles-DSMCReservoirSim')
DSMC%ReservoirSimuRate       = GETLOGICAL('Particles-DSMCReservoirSimRate')
DSMC%ReservoirRateStatistic  = GETLOGICAL('Particles-DSMCReservoirStatistic')
!-----------------------------------------------------------------------------------
! Initialization of collision mode and internal energy modelling (rotational, vibrational and electronic relaxation)
!-----------------------------------------------------------------------------------
CollisMode = GETINT('Particles-DSMC-CollisMode')
SelectionProc = GETINT('Particles-DSMC-SelectionProcedure','1') ! 1: Laux, 2: Gimelsheim
IF(CollisMode.GE.2) THEN
  DSMC%RotRelaxProb = GETREAL('Particles-DSMC-RotRelaxProb')
  DSMC%VibRelaxProb = GETREAL('Particles-DSMC-VibRelaxProb')
  DSMC%RotRelaxModel = GETINT('Particles-DSMC-RotationalRelaxModel')
  DSMC%VibRelaxModel = GETINT('Particles-DSMC-VibrationalRelaxModel')
  ! Rotational relaxation: set function pointer for rotational energy calculation (no pointer for initial diatomic particle insertion)
  IF(DSMC%RotRelaxModel.EQ.1)THEN
    RotRelaxPolyRoutineFuncPTR => DSMC_RotRelaxQuantPoly
    RotInitPolyRoutineFuncPTR  => CalcERotQuant_particle
    RotRelaxDiaRoutineFuncPTR  => DSMC_RotRelaxDiaQuant
  ELSE IF(DSMC%RotRelaxModel.EQ.2)THEN
    CALL ABORT(__STAMP__,'Rotational Relaxation with database of energy levels and degeneracies not tested yet!')
    RotRelaxPolyRoutineFuncPTR => DSMC_RotRelaxDatabasePoly
    RotInitPolyRoutineFuncPTR  => CalcERotDataset_particle
    RotRelaxDiaRoutineFuncPTR  => DSMC_RotRelaxDatabasePoly
  ELSE
    RotRelaxPolyRoutineFuncPTR => DSMC_RotRelaxPoly
    RotInitPolyRoutineFuncPTR  => CalcERot_particle
    RotRelaxDiaRoutineFuncPTR  => DSMC_RotRelaxDiaContinuous
  END IF
  ! Vibrational relaxation
  IF (DSMC%VibRelaxModel.EQ.1) THEN
    DSMC%VibAHO = .TRUE.
    ALLOCATE(AHO%NumVibLevels(nSpecies))
    ALLOCATE(AHO%omegaE(nSpecies))
    ALLOCATE(AHO%chiE(nSpecies))
    IF (DSMC%VibRelaxProb.EQ.2.) THEN
      CALL Abort(__STAMP__,'ERROR: Variable vibration relaxation probabilities according to Boyd not supported with AHO!')
    END IF
  END IF
ELSE
  DSMC%RotRelaxProb = 0.
  DSMC%VibRelaxProb = 0.
END IF
#if ((PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
IF(DSMC%RotRelaxProb.GT.1.0.OR.DSMC%VibRelaxProb.GT.1.0)THEN
  CALL Abort(__STAMP__,'ERROR: Rotational and vibrational relaxation probabilities must be between 0 and 1 for BGK or FP!')
END IF
#endif /*((PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/
! If granular species exist, the internal degree of freedom needs to be used (restart and communication)
! as long as the vibration energy is used as the bulk temperature value for granular species
IF((CollisMode.LT.2).AND.ANY(Species(:)%InterID.EQ.100)) THEN
  CollisMode = 2
END IF
DSMC%ElectronicModel         = GETINT('Particles-DSMC-ElectronicModel')
IF(SampleElecExcitation.AND.(DSMC%ElectronicModel.NE.3)) CALL CollectiveStop(__STAMP__,&
    'Part-SampleElectronicExcitation = T requires Particles-DSMC-ElectronicModel = 3')
IF (DSMC%ElectronicModel.GT.0) THEN
  ! Allocate internal energy array WITH electronic energy
  ALLOCATE(PartStateIntEn(1:3,PDM%maxParticleNumber))
  ! Allocate model-specific variables
  SELECT CASE(DSMC%ElectronicModel)
    CASE(1) ! Model by Liechty, each particle has a specific electronic state
    CASE(2) ! Model by Burt, each particle has an electronic distribution function
      IF(.NOT.ALLOCATED(ElectronicDistriPart)) ALLOCATE(ElectronicDistriPart(PDM%maxParticleNumber))
    CASE(3) ! MCC model, utilizing cross-section data for specific levels
      IF(SelectionProc.EQ.2) THEN
        CALL Abort(__STAMP__,'ERROR: Selected combination of -SelectionProcedure (2) and -ElectronicModel (3) not supported!')
      END IF
    CASE(4) ! Landau-Teller based model, relaxation of given distribution function
      ALLOCATE(ElecRelaxPart(1:PDM%maxParticleNumber))
      ElecRelaxPart = .TRUE.
    CASE DEFAULT
      CALL Abort(__STAMP__,'ERROR: Please select an electronic model between 1 and 4!')
  END SELECT
ELSE
  ! Allocate internal energy array WITHOUT electronic energy
  ALLOCATE(PartStateIntEn(1:2,PDM%maxParticleNumber))
ENDIF

PartStateIntEn = 0. ! nullify

DSMC%ElectronicModelDatabase = TRIM(GETSTR('Particles-DSMCElectronicDatabase','none'))
IF (SpeciesDatabase.EQ.'none') THEN
  IF ((DSMC%ElectronicModelDatabase .NE. 'none').AND.((CollisMode .GT. 1).OR.(CollisMode .EQ. 0))) THEN
    ! CollisMode=0 is for use of in PIC simulation without collisions
    DSMC%EpsElecBin = GETREAL('EpsMergeElectronicState','1E-4')
  ELSEIF(DSMC%ElectronicModel.EQ.1.OR.DSMC%ElectronicModel.EQ.2.OR.DSMC%ElectronicModel.EQ.4) THEN
    CALL Abort(__STAMP__,'ERROR: Electronic models 1, 2 & 4 require an electronic levels database and CollisMode > 1!')
  END IF
ELSE
  IF ((CollisMode .GT. 1).OR.(CollisMode .EQ. 0)) THEN
    DSMC%EpsElecBin = GETREAL('EpsMergeElectronicState','1E-4')
  END IF
END IF

DSMC%DoTEVRRelaxation        = GETLOGICAL('Particles-DSMC-TEVR-Relaxation')
IF(UseVarTimeStep.OR.usevMPF) THEN
  IF(DSMC%DoTEVRRelaxation) THEN
    CALL abort(__STAMP__,'ERROR: Radial weighting or variable time step is not implemented with T-E-V-R relaxation!')
  END IF
END IF

!-----------------------------------------------------------------------------------
! Initialization of quality factors: collision probability, mean collision separation over mean free path, relaxation factor (BGK)
!-----------------------------------------------------------------------------------
DSMC%CalcQualityFactors = GETLOGICAL('Particles-DSMC-CalcQualityFactors')
IF(DSMC%CalcQualityFactors) THEN
  IF (CollisMode.LT.1) THEN
    CALL abort(__STAMP__,'ERROR: Do not use DSMC%CalcQualityFactors for CollisMode < 1')
  END IF
  ! 1: Maximal collision probability per cell/subcells (octree)
  ! 2: Mean collision probability within cell
  ! 3: Mean collision separation distance over mean free path
  ! 4: Counter (is not simply the number of iterations in case of a coupled BGK/FP-DSMC simulation)
  VarNum = 4
  ! VarNum + 1: Number of cloned particles per cell
  ! VarNum + 2: Number of identical particles (no relative velocity)
  IF(ParticleWeighting%PerformCloning) VarNum = VarNum + 2
  ALLOCATE(DSMC%QualityFacSamp(nElems,VarNum))
  DSMC%QualityFacSamp(1:nElems,1:VarNum) = 0.0
END IF

IF (nSpecies.LE.0) THEN
  CALL Abort(__STAMP__,"ERROR: Number of simulation species must be greater zero! nSpecies: ", nSpecies)
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
  ! Read-in from species database based on SpeciesName
  IF(SpeciesDatabase.NE.'none') THEN
    ! Initialize FORTRAN interface.
    CALL H5OPEN_F(err)
    CALL H5FOPEN_F (TRIM(SpeciesDatabase), H5F_ACC_RDONLY_F, file_id_specdb, err)
    DO iSpec = 1, nSpecies
      ! Skip species that are not read-in from the database
      IF(Species(iSpec)%DoOverwriteParameters) CYCLE
      ! averagedCollisionParameters set true: species-specific collision parameters get read in
      IF(CollInf%averagedCollisionParameters) THEN
        LBWRITE (UNIT_stdOut,'(68(". "))')
        CALL PrintOption('Species Name','INFO',StrOpt=TRIM(Species(iSpec)%Name))
        dsetname = TRIM('/Species/'//TRIM(Species(iSpec)%Name))
        ! Reference temperature
        CALL ReadAttribute(file_id_specdb,'Tref',1,DatasetName = dsetname,RealScalar=SpecDSMC(iSpec)%Tref,ReadFromGroup=.TRUE.)
        CALL PrintOption('Tref','DB',RealOpt=SpecDSMC(iSpec)%Tref)
        ! Reference diameter
        CALL ReadAttribute(file_id_specdb,'dref',1,DatasetName = dsetname,RealScalar=SpecDSMC(iSpec)%dref,ReadFromGroup=.TRUE.)
        CALL PrintOption('dref','DB',RealOpt=SpecDSMC(iSpec)%dref)
        ! Viscosity exponent
        CALL ReadAttribute(file_id_specdb,'omega',1,DatasetName = dsetname,RealScalar=SpecDSMC(iSpec)%omega,ReadFromGroup=.TRUE.)
        CALL PrintOption('omega','DB',RealOpt=SpecDSMC(iSpec)%omega)
        ! Alpha parameter (VSS only)
        CALL AttributeExists(file_id_specdb,'alphaVSS',TRIM(dsetname), AttrExists=AttrExists,ReadFromGroup=.TRUE.)
        IF (AttrExists) THEN
          CALL ReadAttribute(file_id_specdb,'alphaVSS',1,DatasetName = dsetname,RealScalar=SpecDSMC(iSpec)%alphaVSS,ReadFromGroup=.TRUE.)
        ELSE
          SpecDSMC(iSpec)%alphaVSS = 1.0
        END IF
        CALL PrintOption('alphaVSS','DB',RealOpt=SpecDSMC(iSpec)%alphaVSS)
      END IF ! averagedCollisionParameters
      ! Flag to identify fully ionized species
      CALL AttributeExists(file_id_specdb,'FullyIonized',TRIM(dsetname), AttrExists=AttrExists,ReadFromGroup=.TRUE.)
      IF (AttrExists) THEN
        CALL ReadAttribute(file_id_specdb,'FullyIonized',1,DatasetName = dsetname,LogicalScalar=SpecDSMC(iSpec)%FullyIonized,ReadFromGroup=.TRUE.)
      ELSE
        SpecDSMC(iSpec)%FullyIonized = .FALSE.
      END IF
    END DO ! iSpec = nSpecies
    ! Close the file.
    CALL H5FCLOSE_F(file_id_specdb, err)
    ! Close FORTRAN interface.
    CALL H5CLOSE_F(err)
  END IF !database

  ! Read-in from parameter file
  DO iSpec = 1, nSpecies
    IF(Species(iSpec)%DoOverwriteParameters) THEN
      LBWRITE (UNIT_stdOut,'(68(". "))')
      WRITE(UNIT=hilf,FMT='(I0)') iSpec
      SpecDSMC(iSpec)%Tref     = GETREAL('Part-Species'//TRIM(hilf)//'-Tref'     )
      SpecDSMC(iSpec)%dref     = GETREAL('Part-Species'//TRIM(hilf)//'-dref'     )
      SpecDSMC(iSpec)%omega    = GETREAL('Part-Species'//TRIM(hilf)//'-omega'    )
      SpecDSMC(iSpec)%alphaVSS = GETREAL('Part-Species'//TRIM(hilf)//'-alphaVSS' )
      SpecDSMC(iSpec)%FullyIonized  = GETLOGICAL('Part-Species'//TRIM(hilf)//'-FullyIonized')
    END IF
  END DO !iSpec

  ! Sanity check of VHS/VSS parameters
  DO iSpec = 1, nSpecies
    ! Skip sanity check for granular species
    IF(Species(iSpec)%InterID.EQ.100) CYCLE
    WRITE(UNIT=hilf,FMT='(I0)') iSpec
    ! check for faulty parameters
    IF((Species(iSpec)%InterID * SpecDSMC(iSpec)%Tref * SpecDSMC(iSpec)%dref * SpecDSMC(iSpec)%alphaVSS) .EQ. 0) THEN
      CALL Abort(__STAMP__,'ERROR in species data: check collision parameters \n'//&
        'Part-Species'//TRIM(hilf)//'-(InterID * Tref * dref * alphaVSS) .EQ. 0 - but must not be 0!')
    END IF ! (Tref * dref * alphaVSS) .EQ. 0
    ! omega is defined between 0 (= hard sphere) and 0.5 (= Maxwell molecule), CAUTION: omega_PICLas = omega_Bird1994 - 0.5
    IF ((SpecDSMC(iSpec)%omega.LT.0.0) .OR. (SpecDSMC(iSpec)%omega.GT.0.5)) THEN
      CALL Abort(__STAMP__,'ERROR: Check set parameter Part-Species'//TRIM(hilf)//'-omega, which must be between 0 and 0.5 (CAUTION: omega_PICLas = omega_Bird1994 - 0.5)!')
    END IF
    ! alphaVSS is defined between 0 and 2
    IF ((SpecDSMC(iSpec)%alphaVSS.LT.0.0) .OR. (SpecDSMC(iSpec)%alphaVSS.GT.2.0)) THEN
      CALL Abort(__STAMP__,'ERROR: Check set parameter Part-Species'//TRIM(hilf)//'-alphaVSS, which must not be lower 0 or greater 2!')
    END IF ! alphaVSS parameter check
  END DO

  ! Granular species
  SpecDSMC(:)%ThermalACCGranularPart         = 0.0
  SpecDSMC(:)%SpecificHeatSolid              = 0.0
  IF(UseGranularSpecies) THEN
    DO iSpec = 1, nSpecies
      WRITE(UNIT=hilf,FMT='(I0)') iSpec
      IF (Species(iSpec)%InterID.EQ.100) THEN
        SpecDSMC(iSpec)%SpecificHeatSolid              = GETREAL('Part-Species'//TRIM(hilf)//'-GranularPartCsp' )
        IF(ALMOSTZERO(SpecDSMC(iSpec)%SpecificHeatSolid)) THEN
          CALL Abort(__STAMP__,'ERROR in species data: check speciﬁc heat [J/(kg*K)] of granular species. It must not be 0')
        END IF
      ELSE
        SpecDSMC(iSpec)%ThermalACCGranularPart         = GETREAL('Part-Species'//TRIM(hilf)//'-GranularPartTau' )
      END IF
    END DO
  END IF

  DO iSpec=1, nSpecies
    ! Save the electron species into a global variable
    IF(Species(iSpec)%InterID.EQ.4) DSMC%ElectronSpecies = iSpec
    ! reading electronic state information from HDF5 file
    IF(((DSMC%ElectronicModelDatabase.NE.'none').OR.(SpeciesDatabase.NE.'none')).AND.(Species(iSpec)%InterID.NE.4).AND.(Species(iSpec)%InterID.NE.100)) THEN
      CALL SetElectronicModel(iSpec)
    END IF
  END DO

  ! determine number of different species combinations and allocate collidingSpecies array
  CollInf%NumCase = 0
  ALLOCATE(CollInf%Coll_Case(nSpecies,nSpecies))
  DO iSpec = 1, nSpecies
    DO jSpec = iSpec, nSpecies
      CollInf%NumCase = CollInf%NumCase + 1
      CollInf%Coll_Case(iSpec,jSpec) = CollInf%NumCase
      CollInf%Coll_Case(jSpec,iSpec) = CollInf%NumCase
    END DO
  END DO
  ALLOCATE(CollInf%collidingSpecies(CollInf%NumCase,2))
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
    DO iColl = 1, CollInf%NumCase
      WRITE(UNIT=hilf,FMT='(I0)')  iColl
      CollInf%collidingSpecies(iColl,:) = GETINTARRAY('Part-Collision'//TRIM(hilf)//'-partnerSpecies',2,'0,0')
    END DO ! iColl = CollInf%NumCase
  END IF ! averagedCollisionParameters
  DO iColl = 1, CollInf%NumCase ! check if any collidingSpecies pair is set multiple times
    WRITE(UNIT=hilf,FMT='(I0)') iColl
    DO pColl = 1,2 ! collision partner
      WRITE (UNIT = hilf2,FMT = '(I0)') pColl
      IF (CollInf%collidingSpecies(iColl,pColl) .EQ. 0) THEN
          CALL Abort(__STAMP__,'ERROR: Partner species '//TRIM(hilf2)//' for Collision'//TRIM(hilf)//' not defined. '// &
                    'Part-Collision'//TRIM(hilf)//'-partnerSpecies required ')
      END IF ! collidingSpecies .EQ. 0
      IF (CollInf%collidingSpecies(iColl,pColl).GT.nSpecies) THEN
        CALL Abort(__STAMP__,'ERROR: Partner species '//TRIM(hilf2)//' for Collision'//TRIM(hilf)//' .GT. nSpecies')
      END IF
    END DO ! pColl = 2
    DO jColl=1, CollInf%NumCase
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
  DO iColl = 1, CollInf%NumCase
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
      CALL Abort(__STAMP__,&
        'ERROR: Check collision parameters! (Part-Collision'//TRIM(hilf)//'-Tref * dref * alphaVSS) .EQ. 0 - but must not be 0)')
    END IF ! check if collision parameters are set
    ! omega is defined between 0 (= hard sphere) and 0.5 (= Maxwell molecule), CAUTION: omega_PICLas = omega_Bird1994 - 0.5
    IF ((CollInf%omega(iSpec,jSpec).LT.0.0) .OR. (CollInf%omega(iSpec,jSpec).GT.0.5)) THEN
      CALL Abort(__STAMP__,'ERROR: Check set parameter Part-Collision'//TRIM(hilf)//'-omega, which must be between 0 and 0.5 (CAUTION: omega_PICLas = omega_Bird1994 - 0.5)!')
    END IF
    IF ((CollInf%alphaVSS(iSpec,jSpec).LT.1) .OR. (CollInf%alphaVSS(iSpec,jSpec).GT.2)) THEN
      CALL Abort(__STAMP__,&
        'ERROR: Check set parameter Part-Collision'//TRIM(hilf)//'-alphaVSS must not be lower 1 or greater 2')
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

IF (CollisMode.EQ.0) THEN
  IF (DSMC%ReservoirSimu) THEN
    CALL Abort(__STAMP__, "Free Molecular Flow (CollisMode=0) is not supported for reservoir!")
  END IF
ELSE !CollisMode.GT.0
  ! species and case assignment arrays
  ALLOCATE(DSMC%NumColl(CollInf%NumCase +1))
  DSMC%NumColl = 0.
  ALLOCATE(CollInf%Coll_CaseNum(CollInf%NumCase))
  CollInf%Coll_CaseNum = 0
  ALLOCATE(CollInf%Coll_SpecPartNum(nSpecies))
  CollInf%Coll_SpecPartNum = 0.
  ALLOCATE(CollInf%SumPairMPF(CollInf%NumCase))
  CollInf%SumPairMPF = 0.
  ALLOCATE(CollInf%FracMassCent(nSpecies, CollInf%NumCase)) ! Calculation of mx/(mx+my) and reduced mass
  CollInf%FracMassCent = 0
  ALLOCATE(CollInf%MassRed(CollInf%NumCase))
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
  ALLOCATE(CollInf%Cab(CollInf%NumCase))
  ALLOCATE(CollInf%KronDelta(CollInf%NumCase))
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
  ! reading/writing molecular stuff
  !-----------------------------------------------------------------------------------------------------------------------------------
  IF ((CollisMode.EQ.2).OR.(CollisMode.EQ.3)) THEN ! perform relaxation (molecular) reactions
    ! reading molecular stuff
    SpecDSMC(1:nSpecies)%Xi_Rot = 0
    SpecDSMC(1:nSpecies)%MaxVibQuant = 0
    SpecDSMC(1:nSpecies)%CharaTVib = 0
    SpecDSMC(1:nSpecies)%EZeroPoint = 0.0
    SpecDSMC(1:nSpecies)%PolyatomicMol = .FALSE.
    SpecDSMC(1:nSpecies)%SpecToPolyArray = 0
    useRelaxProbCorrFactor=GETLOGICAL('Particles-DSMC-useRelaxProbCorrFactor','.FALSE.')
    IF(DSMC%VibAHO.AND.useRelaxProbCorrFactor) THEN
      CALL Abort(__STAMP__,'ERROR: Utilization of vibrational relaxation probability correction factor not possible with AHO!')
    END IF
    IF (DSMC%VibRelaxProb.EQ.3.0) THEN
      ALLOCATE(SpecDSMC(nSpecies)%C1(nSpecies))
      ALLOCATE(SpecDSMC(nSpecies)%C2(nSpecies))
    END IF

    IF(SpeciesDatabase.NE.'none') THEN
      ! Initialize FORTRAN interface.
      CALL H5OPEN_F(err)
      CALL H5FOPEN_F (TRIM(SpeciesDatabase), H5F_ACC_RDONLY_F, file_id_specdb, err)
      DO iSpec = 1, nSpecies
        IF(.NOT.Species(iSpec)%DoOverwriteParameters) THEN
          dsetname = TRIM('/Species/'//TRIM(Species(iSpec)%Name))
          WRITE(UNIT=hilf,FMT='(I0)') iSpec
          IF((Species(iSpec)%InterID.NE.4).AND.(Species(iSpec)%InterID.NE.100)) THEN
            dsetname = TRIM('/Species/'//TRIM(Species(iSpec)%Name))
            CALL AttributeExists(file_id_specdb,'PolyatomicMol',TRIM(dsetname),AttrExists=AttrExists,ReadFromGroup=.TRUE.)
            IF (AttrExists) THEN
              CALL ReadAttribute(file_id_specdb,'PolyatomicMol',1,DatasetName = dsetname,IntScalar=IntToLog, &
                ReadFromGroup=.TRUE.)
              IF(IntToLog.EQ.1) SpecDSMC(iSpec)%PolyatomicMol = .TRUE.
              CALL PrintOption('Species Name','INFO',StrOpt=TRIM(Species(iSpec)%Name))
              CALL PrintOption('PolyatomicMol','DB',LogOpt=SpecDSMC(iSpec)%PolyatomicMol)
            ELSE
              SpecDSMC(iSpec)%PolyatomicMol = .FALSE.
            END IF
            IF(SpecDSMC(iSpec)%PolyatomicMol.AND.DSMC%DoTEVRRelaxation)  THEN
              CALL Abort(__STAMP__,'! Simulation of Polyatomic Molecules and T-E-V-R relaxation not possible yet!!!')
            END IF
            IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
              IF(DSMC%VibAHO) THEN
                CALL Abort(__STAMP__,'ERROR: The anharmonic model is not implemented for polyatomic species yet!')
              END IF
              DSMC%NumPolyatomMolecs = DSMC%NumPolyatomMolecs + 1
              SpecDSMC(iSpec)%SpecToPolyArray = DSMC%NumPolyatomMolecs
            ELSE IF ((Species(iSpec)%InterID.EQ.2).OR.(Species(iSpec)%InterID.EQ.20)) THEN
              SpecDSMC(iSpec)%Xi_Rot     = 2
              CALL ReadAttribute(file_id_specdb,'CharaTempVib',1,DatasetName = dsetname, &
                RealScalar=SpecDSMC(iSpec)%CharaTVib,ReadFromGroup=.TRUE.)
              CALL PrintOption('CharaTempVib','DB',RealOpt=SpecDSMC(iSpec)%CharaTVib)
              IF(DSMC%RotRelaxModel.EQ.1)THEN
                CALL AttributeExists(file_id_specdb,'MomentOfInertia',TRIM(dsetname), AttrExists=AttrExists, &
                  ReadFromGroup=.TRUE.)
                IF (AttrExists) THEN
                  CALL ReadAttribute(file_id_specdb,'MomentOfInertia',1,DatasetName = dsetname, &
                    RealScalar=SpecDSMC(iSpec)%MomentOfInertia,ReadFromGroup=.TRUE.)
                  SpecDSMC(iSpec)%CharaTRot = PlanckConst**2 / (8 * PI**2 * SpecDSMC(iSpec)%MomentOfInertia * BoltzmannConst)
                  CALL PrintOption('MomentOfInertia','DB',RealOpt=SpecDSMC(iSpec)%MomentOfInertia)
                ELSE
                  CALL abort(__STAMP__,'Moment of inertia necessary for quantized rotational energy and is not set for species '&
                    //(Species(iSpec)%Name))
                END IF
              ELSE
                CALL AttributeExists(file_id_specdb,'CharaTempRot',TRIM(dsetname), AttrExists=AttrExists, &
                  ReadFromGroup=.TRUE.)
                IF(AttrExists)THEN  ! check if CharaTempRot is set without Moment of Inertia
                  CALL ReadAttribute(file_id_specdb,'CharaTempRot',1,DatasetName = dsetname, &
                    RealScalar=SpecDSMC(iSpec)%CharaTRot,ReadFromGroup=.TRUE.)
                ELSE
                  CALL AttributeExists(file_id_specdb,'MomentOfInertia',TRIM(dsetname), AttrExists=AttrExists, &
                    ReadFromGroup=.TRUE.)
                  IF (AttrExists) THEN
                    CALL ReadAttribute(file_id_specdb,'MomentOfInertia',1,DatasetName = dsetname, &
                      RealScalar=SpecDSMC(iSpec)%MomentOfInertia,ReadFromGroup=.TRUE.)
                    SpecDSMC(iSpec)%CharaTRot = PlanckConst**2 / (8 * PI**2 * SpecDSMC(iSpec)%MomentOfInertia * BoltzmannConst)
                  ELSE
                    SpecDSMC(iSpec)%CharaTRot = 0.0
                  END IF
                END IF
              END IF
              CALL PrintOption('CharaTempRot','DB',RealOpt=SpecDSMC(iSpec)%CharaTRot)
              CALL ReadAttribute(file_id_specdb,'Ediss_eV',1,DatasetName = dsetname,RealScalar=SpecDSMC(iSpec)%Ediss_eV, &
                ReadFromGroup=.TRUE.)
              CALL PrintOption('Ediss_eV','DB',RealOpt=SpecDSMC(iSpec)%Ediss_eV)
              ! Anharmonic Oscillator Model
              IF(DSMC%VibAHO) THEN
                CALL ReadAttribute(file_id_specdb,'Vib-OmegaE',1,DatasetName = dsetname,RealScalar=AHO%omegaE(iSpec), &
                  ReadFromGroup=.TRUE.)
                CALL PrintOption('Vib-OmegaE','DB',RealOpt=AHO%omegaE(iSpec))
                CALL AttributeExists(file_id_specdb,'Vib-ChiE ',TRIM(dsetname), AttrExists=AttrExists,ReadFromGroup=.TRUE.)
                CALL ReadAttribute(file_id_specdb,'Vib-ChiE',1,DatasetName = dsetname,RealScalar=AHO%chiE(iSpec), &
                  ReadFromGroup=.TRUE.)
                CALL PrintOption('Vib-ChiE','DB',RealOpt=AHO%chiE(iSpec))
              ELSE ! SHO model
                CALL ReadAttribute(file_id_specdb,'CharaTempVib',1,DatasetName = dsetname,RealScalar=SpecDSMC(iSpec)%CharaTVib, &
                  ReadFromGroup=.TRUE.)
                CALL PrintOption('CharaTempVib','DB',RealOpt=SpecDSMC(iSpec)%CharaTVib)
                SpecDSMC(iSpec)%MaxVibQuant = 200
                ! Calculation of the dissociation quantum number (used for QK chemistry)
                SpecDSMC(iSpec)%DissQuant = INT(SpecDSMC(iSpec)%Ediss_eV*ElementaryCharge/(BoltzmannConst*SpecDSMC(iSpec)%CharaTVib))
                ! Calculation of the zero-point energy
                SpecDSMC(iSpec)%EZeroPoint = DSMC%GammaQuant * BoltzmannConst * SpecDSMC(iSpec)%CharaTVib
              END IF
            END IF
            ! Read in species values for rotational relaxation models of Boyd/Zhang if necessary
            IF(DSMC%RotRelaxProb.GT.1.0.AND.DSMC%RotRelaxProb.LE.3.0.AND.((Species(iSpec)%InterID.EQ.2).OR. &
              (Species(iSpec)%InterID.EQ.20))) THEN
              SpecDSMC(iSpec)%CollNumRotInf = GETREAL('Part-Species'//TRIM(hilf)//'-CollNumRotInf')
              SpecDSMC(iSpec)%TempRefRot    = GETREAL('Part-Species'//TRIM(hilf)//'-TempRefRot')
              IF(SpecDSMC(iSpec)%CollNumRotInf*SpecDSMC(iSpec)%TempRefRot.EQ.0) THEN
                CALL Abort(__STAMP__,'Error! CollNumRotRef or TempRefRot is equal to zero for species:', iSpec)
              END IF
            END IF
            ! Read in species values for vibrational relaxation models of Milikan-White if necessary
            IF(DSMC%VibRelaxProb.EQ.2.0) THEN
              ! Only molecules or charged molecules
              IF(((Species(iSpec)%InterID.EQ.2).OR.(Species(iSpec)%InterID.EQ.20))) THEN
                ALLOCATE(SpecDSMC(iSpec)%MW_ConstA(1:nSpecies))
                ALLOCATE(SpecDSMC(iSpec)%MW_ConstB(1:nSpecies))
                DO jSpec = 1, nSpecies
                  WRITE(UNIT=hilf2,FMT='(I0)') jSpec
                  hilf2=TRIM(hilf)//'-'//TRIM(hilf2)
                  SpecDSMC(iSpec)%MW_ConstA(jSpec)     = GETREAL('Part-Species'//TRIM(hilf)//'-MWConstA-'//TRIM(hilf2))
                  SpecDSMC(iSpec)%MW_ConstB(jSpec)     = GETREAL('Part-Species'//TRIM(hilf)//'-MWConstB-'//TRIM(hilf2))

                  IF(SpecDSMC(iSpec)%MW_ConstA(jSpec).EQ.0) THEN
                    CALL Abort(__STAMP__,'Error! MW_ConstA is equal to zero for species:', iSpec)
                  END IF
                  IF(SpecDSMC(iSpec)%MW_ConstB(jSpec).EQ.0) THEN
                    CALL Abort(__STAMP__,'Error! MW_ConstB is equal to zero for species:', iSpec)
                  END IF
                END DO
              END IF
              SpecDSMC(iSpec)%VibCrossSec    = GETREAL('Part-Species'//TRIM(hilf)//'-VibCrossSection')
              ! Only molecules or charged molecules
              IF((SpecDSMC(iSpec)%VibCrossSec.EQ.0).AND.((Species(iSpec)%InterID.EQ.2).OR.(Species(iSpec)%InterID.EQ.20))) THEN
                CALL Abort(__STAMP__,'Error! VibCrossSec is equal to zero for species:', iSpec)
              END IF
            END IF
            ! Setting the values of Rot-/Vib-RelaxProb to a fix value (electronic: species-specific values are possible)
            SpecDSMC(iSpec)%RotRelaxProb  = DSMC%RotRelaxProb
            SpecDSMC(iSpec)%VibRelaxProb  = DSMC%VibRelaxProb
            SpecDSMC(iSpec)%ElecRelaxProb = GETREAL('Part-Species'//TRIM(hilf)//'-ElecRelaxProb')
            ! multi init stuff
            ALLOCATE(SpecDSMC(iSpec)%Init(0:Species(iSpec)%NumberOfInits))
            ! Skip the read-in of temperatures if a background gas distribution is used but not if background gas regions are used
            IF(BGGas%NumberOfSpecies.GT.0) THEN
              IF(BGGas%BackgroundSpecies(iSpec).AND.BGGas%UseDistribution.AND.(.NOT.BGGas%UseRegions)) THEN
                SpecDSMC(iSpec)%Init(1)%TVib  = 0.
                SpecDSMC(iSpec)%Init(1)%TRot  = 0.
                SpecDSMC(iSpec)%Init(1)%Telec = 0.
                CYCLE
              END IF
            END IF
            DO iInit = 1, Species(iSpec)%NumberOfInits
              WRITE(UNIT=hilf2,FMT='(I0)') iInit
              hilf2=TRIM(hilf)//'-Init'//TRIM(hilf2)
              IF(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'EmissionDistribution')THEN
                SpecDSMC(iSpec)%Init(iInit)%TVib      = GETREAL('Part-Species'//TRIM(hilf2)//'-TempVib','300.0')
                SpecDSMC(iSpec)%Init(iInit)%TRot      = GETREAL('Part-Species'//TRIM(hilf2)//'-TempRot','300.0')
              ELSEIF((Species(iSpec)%InterID.EQ.2).OR.(Species(iSpec)%InterID.EQ.20)) THEN
                SpecDSMC(iSpec)%Init(iInit)%TVib      = GETREAL('Part-Species'//TRIM(hilf2)//'-TempVib')
                SpecDSMC(iSpec)%Init(iInit)%TRot      = GETREAL('Part-Species'//TRIM(hilf2)//'-TempRot')
              END IF
              ! read electronic temperature
              IF(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'EmissionDistribution')THEN
                SpecDSMC(iSpec)%Init(iInit)%Telec   = GETREAL('Part-Species'//TRIM(hilf2)//'-TempElec','300.0')
              ELSEIF (DSMC%ElectronicModel.GT.0) THEN
                SpecDSMC(iSpec)%Init(iInit)%Telec   = GETREAL('Part-Species'//TRIM(hilf2)//'-TempElec')
              END IF ! electronic model
            END DO !Inits
            ALLOCATE(SpecDSMC(iSpec)%Surfaceflux(1:Species(iSpec)%nSurfacefluxBCs))
            DO iInit = 1, Species(iSpec)%nSurfacefluxBCs
              WRITE(UNIT=hilf2,FMT='(I0)') iInit
              hilf2=TRIM(hilf)//'-Surfaceflux'//TRIM(hilf2)
              IF((Species(iSpec)%InterID.EQ.2).OR.(Species(iSpec)%InterID.EQ.20)) THEN
                SpecDSMC(iSpec)%Surfaceflux(iInit)%TVib      = GETREAL('Part-Species'//TRIM(hilf2)//'-TempVib')
                SpecDSMC(iSpec)%Surfaceflux(iInit)%TRot      = GETREAL('Part-Species'//TRIM(hilf2)//'-TempRot')
                IF (SpecDSMC(iSpec)%Surfaceflux(iInit)%TRot*SpecDSMC(iSpec)%Surfaceflux(iInit)%TVib.EQ.0.) THEN
                  CALL Abort(__STAMP__,'Error! TVib and TRot not def. in Part-SpeciesXX-SurfacefluxXX-TempVib/TempRot for iSpec, iInit',iSpec,REAL(iInit))
                END IF
              END IF
              ! read electronic temperature
              IF (DSMC%ElectronicModel.GT.0) THEN
                SpecDSMC(iSpec)%Surfaceflux(iInit)%Telec   = GETREAL('Part-Species'//TRIM(hilf2)//'-TempElec')
                IF (SpecDSMC(iSpec)%Surfaceflux(iInit)%Telec.EQ.0.) THEN
                  CALL Abort(__STAMP__,' Error! Telec not defined in Part-SpeciesXX-SurfacefluxXX-Tempelec for iSpec, iInit',iSpec,REAL(iInit))
                END IF
              END IF
            END DO !SurfaceFluxBCs
            ! read rotational levels from database
            IF(DSMC%RotRelaxModel.EQ.2) CALL ReadRotationalSpeciesLevel(iSpec)
          END IF ! not electron
        END IF
      END DO !Species

      ! Close the file.
      CALL H5FCLOSE_F(file_id_specdb, err)
      ! Close FORTRAN interface.
      CALL H5CLOSE_F(err)

      ! Finalize read-in for vibrational energy levels of anharmonic oscillator model
      IF(DSMC%VibAHO) THEN
        DO iSpec = 1, nSpecies
          IF(.NOT.Species(iSpec)%DoOverwriteParameters) THEN
            IF ((Species(iSpec)%InterID.EQ.2).OR.(Species(iSpec)%InterID.EQ.20)) THEN
              CALL ReadVibrationalSpeciesLevel(iSpec)
              SpecDSMC(iSpec)%EZeroPoint = AHO%VibEnergy(iSpec,1)
            END IF
          END IF
        END DO
      END IF

    END IF !database

    DO iSpec = 1, nSpecies
      IF(Species(iSpec)%DoOverwriteParameters) THEN
        IF((Species(iSpec)%InterID.NE.4).AND.(Species(iSpec)%InterID.NE.100)) THEN
          LBWRITE (UNIT_stdOut,'(68(". "))')
          WRITE(UNIT=hilf,FMT='(I0)') iSpec
          SpecDSMC(iSpec)%PolyatomicMol=GETLOGICAL('Part-Species'//TRIM(hilf)//'-PolyatomicMol','.FALSE.')
          IF(DSMC%VibAHO.AND.SpecDSMC(iSpec)%PolyatomicMol) THEN
            CALL Abort(__STAMP__,'ERROR: The anharmonic oscillator model is only available for diatomic species!')
          END IF
          IF(SpecDSMC(iSpec)%PolyatomicMol.AND.DSMC%DoTEVRRelaxation)  THEN
            CALL Abort(__STAMP__,'! Simulation of Polyatomic Molecules and T-E-V-R relaxation not possible yet!!!')
          END IF
          IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
            DSMC%NumPolyatomMolecs = DSMC%NumPolyatomMolecs + 1
            SpecDSMC(iSpec)%SpecToPolyArray = DSMC%NumPolyatomMolecs
          ELSEIF ((Species(iSpec)%InterID.EQ.2).OR.(Species(iSpec)%InterID.EQ.20)) THEN
            SpecDSMC(iSpec)%Xi_Rot          = 2
            IF(DSMC%RotRelaxModel.EQ.1)THEN  ! If MomentOfInertia needed calc CharaTempRot with Moment
              SpecDSMC(iSpec)%MomentOfInertia = GETREAL('Part-Species'//TRIM(hilf)//'-MomentOfInertia')
              SpecDSMC(iSpec)%CharaTRot       = PlanckConst**2 / (8 * PI**2 * SpecDSMC(iSpec)%MomentOfInertia * BoltzmannConst)
            ELSE  ! read in CharaTRot
              SpecDSMC(iSpec)%CharaTRot       = GETREAL('Part-Species'//TRIM(hilf)//'-CharaTempRot')
            END IF
            SpecDSMC(iSpec)%Ediss_eV   = GETREAL('Part-Species'//TRIM(hilf)//'-Ediss_eV')
            ! Anharmonic Oscillator Model
            IF(DSMC%VibAHO) THEN
              WRITE(UNIT=hilf,FMT='(I0)') iSpec
              AHO%omegaE(iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Vib-Anharmonic-omegaE')
              AHO%chiE(iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Vib-Anharmonic-chiE')
              FileNameAHO(iSpec) = GETSTR('Particles-DSMC-Vib-Anharmonic-Species'//TRIM(hilf)//'-Data')
              CALL ReadAHOLevelsFromCSV(FileNameAHO(iSpec),iSpec)
            ELSE
              SpecDSMC(iSpec)%CharaTVib  = GETREAL('Part-Species'//TRIM(hilf)//'-CharaTempVib')
              SpecDSMC(iSpec)%MaxVibQuant = 200
              ! Calculation of the dissociation quantum number (used for QK chemistry)
              SpecDSMC(iSpec)%DissQuant = INT(SpecDSMC(iSpec)%Ediss_eV*ElementaryCharge/(BoltzmannConst*SpecDSMC(iSpec)%CharaTVib))
              ! Calculation of the zero-point energy
              SpecDSMC(iSpec)%EZeroPoint = DSMC%GammaQuant * BoltzmannConst * SpecDSMC(iSpec)%CharaTVib
            END IF
          END IF
          ! Read in species values for rotational relaxation models of Boyd/Zhang if necessary
          IF(DSMC%RotRelaxProb.GT.1.0.AND.DSMC%RotRelaxProb.LE.3.0 &
            .AND.((Species(iSpec)%InterID.EQ.2).OR.(Species(iSpec)%InterID.EQ.20))) THEN
            SpecDSMC(iSpec)%CollNumRotInf = GETREAL('Part-Species'//TRIM(hilf)//'-CollNumRotInf')
            SpecDSMC(iSpec)%TempRefRot    = GETREAL('Part-Species'//TRIM(hilf)//'-TempRefRot')
            IF(SpecDSMC(iSpec)%CollNumRotInf*SpecDSMC(iSpec)%TempRefRot.EQ.0) THEN
              CALL Abort(__STAMP__,'Error! CollNumRotRef or TempRefRot is equal to zero for species:', iSpec)
            END IF
          END IF
          ! Read in species values for vibrational relaxation models of Milikan-White if necessary
          IF(DSMC%VibRelaxProb.EQ.2.0) THEN
            ! Only molecules or charged molecules
            IF(((Species(iSpec)%InterID.EQ.2).OR.(Species(iSpec)%InterID.EQ.20))) THEN
              ALLOCATE(SpecDSMC(iSpec)%MW_ConstA(1:nSpecies))
              ALLOCATE(SpecDSMC(iSpec)%MW_ConstB(1:nSpecies))
              DO jSpec = 1, nSpecies
                WRITE(UNIT=hilf2,FMT='(I0)') jSpec
                hilf2=TRIM(hilf)//'-'//TRIM(hilf2)
                SpecDSMC(iSpec)%MW_ConstA(jSpec)     = GETREAL('Part-Species'//TRIM(hilf)//'-MWConstA-'//TRIM(hilf2))
                SpecDSMC(iSpec)%MW_ConstB(jSpec)     = GETREAL('Part-Species'//TRIM(hilf)//'-MWConstB-'//TRIM(hilf2))
                IF(SpecDSMC(iSpec)%MW_ConstA(jSpec).EQ.0) THEN
                  CALL Abort(__STAMP__,'Error! MW_ConstA is equal to zero for species:', iSpec)
                END IF
                IF(SpecDSMC(iSpec)%MW_ConstB(jSpec).EQ.0) THEN
                  CALL Abort(__STAMP__,'Error! MW_ConstB is equal to zero for species:', iSpec)
                END IF
              END DO
              SpecDSMC(iSpec)%VibCrossSec    = GETREAL('Part-Species'//TRIM(hilf)//'-VibCrossSection')
              IF(SpecDSMC(iSpec)%VibCrossSec.EQ.0) CALL Abort(__STAMP__,'Error! VibCrossSec is equal to zero for species:', iSpec)
            END IF
          ELSE IF (DSMC%VibRelaxProb.EQ.3.0) THEN
            DO jSpec = 1, nSpecies
              WRITE(UNIT=hilf2,FMT='(I0)') jSpec
              SpecDSMC(iSpec)%C1(jSpec)     = GETREAL('Part-Species'//TRIM(hilf)//'-'//TRIM(hilf2)//'-VibConstant-C1')
              SpecDSMC(iSpec)%C2(jSpec)     = GETREAL('Part-Species'//TRIM(hilf)//'-'//TRIM(hilf2)//'-VibConstant-C2')
            END DO
          END IF
          ! Setting the values of Rot-/Vib-RelaxProb to a fix value (electronic: species-specific values are possible)
          SpecDSMC(iSpec)%RotRelaxProb  = DSMC%RotRelaxProb
          SpecDSMC(iSpec)%VibRelaxProb  = DSMC%VibRelaxProb
          SpecDSMC(iSpec)%ElecRelaxProb = GETREAL('Part-Species'//TRIM(hilf)//'-ElecRelaxProb')
          ! multi init stuff
          ALLOCATE(SpecDSMC(iSpec)%Init(0:Species(iSpec)%NumberOfInits))
          ! Skip the read-in of temperatures if a background gas distribution is used but not if background gas regions are used
          IF(BGGas%NumberOfSpecies.GT.0) THEN
            IF(BGGas%BackgroundSpecies(iSpec).AND.BGGas%UseDistribution.AND.(.NOT.BGGas%UseRegions)) THEN
              SpecDSMC(iSpec)%Init(1)%TVib  = 0.
              SpecDSMC(iSpec)%Init(1)%TRot  = 0.
              SpecDSMC(iSpec)%Init(1)%Telec = 0.
              CYCLE
            END IF
          END IF
          DO iInit = 1, Species(iSpec)%NumberOfInits
            WRITE(UNIT=hilf2,FMT='(I0)') iInit
            hilf2=TRIM(hilf)//'-Init'//TRIM(hilf2)
            IF((Species(iSpec)%InterID.EQ.2).OR.(Species(iSpec)%InterID.EQ.20)) THEN
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
            IF((Species(iSpec)%InterID.EQ.2).OR.(Species(iSpec)%InterID.EQ.20)) THEN
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
      END IF !iSpec overwrite parameters
    END DO !Species
    ! Finalize read-in for vibrational energy levels of anharmonic oscillator model
    IF(DSMC%VibAHO) THEN
      IF(SpeciesDatabase.EQ.'None') THEN
        DO iSpec = 1, nSpecies
          IF(Species(iSpec)%DoOverwriteParameters) THEN
            IF ((Species(iSpec)%InterID.EQ.2).OR.(Species(iSpec)%InterID.EQ.20)) THEN
              CALL ReadAHOEnergiesFromCSV(FileNameAHO(iSpec),iSpec)
              SpecDSMC(iSpec)%EZeroPoint = AHO%VibEnergy(iSpec,1)
            END IF
          END IF
        END DO
      END IF
    END IF

    ! Initialization of polyatomic species and burn-in phase (Metropolis-Hastings) per initialization region
    IF(DSMC%NumPolyatomMolecs.GT.0) THEN
      DSMC%PolySingleMode = GETLOGICAL('Particles-DSMC-PolyRelaxSingleMode','.FALSE.')
      IF((SelectionProc.NE.2).AND.DSMC%PolySingleMode) THEN
        ! Single-mode relaxation of vibrational modes of polyatomic molecules only possible when the prohibiting double relaxation
        ! method is used (Gimelshein, SelectionProc = 2)
        CALL Abort(__STAMP__,'ERROR: No single-mode polyatomic relaxation possible with chosen selection procedure! SelectionProc:', SelectionProc)
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
    IF (DSMC%ReservoirSimu) THEN
      DSMC%CompareLandauTeller = GETLOGICAL('Particles-DSMC-CompareLandauTeller','.FALSE.')
      IF(DSMC%CompareLandauTeller) THEN
        IF(CollisMode.NE.2) THEN
          CALL abort(__STAMP__,'ERROR: Comparison with Landau-Teller only available in CollisMode = 2, CollisMode:', CollisMode)
        END IF
        IF(nSpecies.GT.1) THEN
          CALL abort(__STAMP__,'ERROR: Comparison with Landau-Teller only available for a single species, nSpecies:', nSpecies)
        END IF
      END IF
    END IF
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! Setting the internal energy value of every particle
    IF(BGGas%UseRegions) THEN
      CALL BGGas_RegionsSetInternalTemp()
    END IF
  END IF ! CollisMode .EQ. 2 or 3
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! Define chemical reactions (including ionization and backward reaction rate)
  !-----------------------------------------------------------------------------------------------------------------------------------
  IF (CollisMode.EQ.3) THEN ! perform chemical reactions
    IF(SpeciesDatabase.NE.'none') THEN
      ! Initialize FORTRAN interface.
      CALL H5OPEN_F(err)
      CALL H5FOPEN_F (TRIM(SpeciesDatabase), H5F_ACC_RDONLY_F, file_id_specdb, err)
      DO iSpec = 1, nSpecies
        dsetname = TRIM('/Species/'//TRIM(Species(iSpec)%Name))
        WRITE(UNIT=hilf,FMT='(I0)') iSpec
        ! Read-in of heat of formation, ions are treated later using the heat of formation of their ground state and data from the
        ! from the electronic state database to ensure consistent energies across chemical reactions of QK and Arrhenius type.
        IF((Species(iSpec)%InterID.EQ.10).OR.(Species(iSpec)%InterID.EQ.20).OR.(Species(iSpec)%InterID.EQ.4)) THEN
          SpecDSMC(iSpec)%HeatOfFormation = 0.0
        ELSE
          CALL ReadAttribute(file_id_specdb,'HeatOfFormation_K',1,DatasetName = dsetname,RealScalar=SpecDSMC(iSpec)%HeatOfFormation,ReadFromGroup=.TRUE.)
          SpecDSMC(iSpec)%HeatOfFormation = SpecDSMC(iSpec)%HeatOfFormation * BoltzmannConst
        ENDIF
        CALL PrintOption('HeatOfFormation [J], '//TRIM(Species(iSpec)%Name),'DB',RealOpt=SpecDSMC(iSpec)%HeatOfFormation)
      END DO ! iSpec = nSpecies
      ! Close the file.
      CALL H5FCLOSE_F(file_id_specdb, err)
      ! Close FORTRAN interface.
      CALL H5CLOSE_F(err)
    END IF !database

    DO iSpec = 1, nSpecies
      WRITE(UNIT=hilf,FMT='(I0)') iSpec
      IF(Species(iSpec)%DoOverwriteParameters) THEN
      ! Read-in of heat of formation, ions are treated later using the heat of formation of their ground state and data from the
      ! from the electronic state database to ensure consistent energies across chemical reactions of QK and Arrhenius type.
        IF((Species(iSpec)%InterID.EQ.10).OR.(Species(iSpec)%InterID.EQ.20).OR.(Species(iSpec)%InterID.EQ.4)) THEN
          SpecDSMC(iSpec)%HeatOfFormation = 0.0
        ELSE
          SpecDSMC(iSpec)%HeatOfFormation = GETREAL('Part-Species'//TRIM(hilf)//'-HeatOfFormation_K')
          SpecDSMC(iSpec)%HeatOfFormation = SpecDSMC(iSpec)%HeatOfFormation * BoltzmannConst
        END IF
      END IF
    END DO
      ! Heat of formation of ionized species is modified with the ionization energy directly from read-in electronic energy levels
      ! of the ground/previous state of the respective species (Input requires a species number (eg species number of N for NIon1))
    IF(SpeciesDatabase.NE.'none') THEN
      ! Initialize FORTRAN interface.
      CALL H5OPEN_F(err)
      CALL H5FOPEN_F (TRIM(SpeciesDatabase), H5F_ACC_RDONLY_F, file_id_specdb, err)
      DO iSpec=1,nSpecies
        dsetname = TRIM('/Species/'//TRIM(Species(iSpec)%Name))
        CALL AttributeExists(file_id_specdb,'PreviousState',TRIM(dsetname), AttrExists=AttrExists,ReadFromGroup=.TRUE.)
        IF (AttrExists) THEN
          CALL ReadAttribute(file_id_specdb,'PreviousState',1,DatasetName = dsetname,StrScalar=hilfname,ReadFromGroup=.TRUE.)
          DO jSpec=1,nSpecies
            IF(TRIM(hilfname).EQ.TRIM(Species(jSpec)%Name)) THEN
              SpecDSMC(iSpec)%PreviousState=jSpec
            END IF
          END DO
        ELSE
          SpecDSMC(iSpec)%PreviousState = 0
        END IF
        IF(SpecDSMC(iSpec)%PreviousState.GT.0) THEN
          CALL PrintOption('Previous species of '//TRIM(Species(iSpec)%Name)//' (ground state or previous ionized level)', &
            'DB',StrOpt=TRIM(Species(SpecDSMC(iSpec)%PreviousState)%Name))
        END IF
      END DO
      ! Close the file.
      CALL H5FCLOSE_F(file_id_specdb, err)
      ! Close FORTRAN interface.
      CALL H5CLOSE_F(err)
    END IF !database

    DO iSpec=1,nSpecies
      WRITE(UNIT=hilf,FMT='(I0)') iSpec
      IF(Species(iSpec)%DoOverwriteParameters) THEN
        IF((Species(iSpec)%InterID.EQ.10).OR.(Species(iSpec)%InterID.EQ.20)) THEN
          SpecDSMC(iSpec)%PreviousState = GETINT('Part-Species'//TRIM(hilf)//'-PreviousState')
        ELSE
          SpecDSMC(iSpec)%PreviousState = 0
        END IF
      END IF
    END DO

    DO iSpec=1,nSpecies
      WRITE(UNIT=hilf,FMT='(I0)') iSpec
      ! Read-in of species for field ionization (only required if it cannot be determined automatically)
      IF((Species(iSpec)%InterID.NE.4).AND.(Species(iSpec)%InterID.NE.100)) THEN
        SpecDSMC(iSpec)%NextIonizationSpecies = GETINT('Part-Species'//TRIM(hilf)//'-NextIonizationSpecies')
      ELSE
        SpecDSMC(iSpec)%NextIonizationSpecies = 0
        DSMC%ElectronSpecies = iSpec
      END IF
      ! Read-in of electronic levels for QK and backward reaction rate -------------------------------------------------------------
      IF ((DSMC%ElectronicModelDatabase .EQ.'none').AND.(SpeciesDatabase.EQ.'none')) THEN
        IF ((Species(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
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
    CALL CalcHeatOfFormationIons()

    ! Set "NextIonizationSpecies" information for field ionization from "PreviousState" info
    ! NextIonizationSpecies => SpeciesID of the next higher ionization level
    CALL SetNextIonizationSpecies()

    CALL DSMC_chemical_init()
  ELSE
    DSMC%BackwardReacRate = .FALSE.
  END IF

    ! Check whether calculation of instantaneous translational temperature is required
  IF(((CollisMode.GT.1).AND.(SelectionProc.EQ.2)).OR.DSMC%BackwardReacRate.OR.DSMC%CalcQualityFactors.OR.useRelaxProbCorrFactor &
            .OR.(DSMC%VibRelaxProb.EQ.2).OR.(DSMC%ElectronicModel.EQ.2).OR.(DSMC%ElectronicModel.EQ.4)) THEN
    ! 1. Case: Inelastic collisions and chemical reactions with the Gimelshein relaxation procedure and variable vibrational
    !           relaxation probability (CalcGammaVib)
    ! 2. Case: Backward reaction rates
    ! 3. Case: Temperature required for the mean free path with the VHS model
    ALLOCATE(DSMC%InstantTransTemp(nSpecies+1))
    DSMC%InstantTransTemp = 0.0
    IF((DSMC%ElectronicModel.EQ.2).OR.useRelaxProbCorrFactor) THEN
      ALLOCATE(DSMC%InstantTXiElec(2,nSpecies))
      DSMC%InstantTXiElec = 0.0
    END IF
    IF (useRelaxProbCorrFactor.AND.(DSMC%ElectronicModel.EQ.1)) THEN
      DO iSpec = 1, nSpecies
        ALLOCATE(SpecDSMC(iSpec)%ElecRelaxCorrectFac(nSpecies))
      END DO
    END IF
  END IF

  !-----------------------------------------------------------------------------------------------------------------------------------
  ! Octree + Nearest Neighbour: octree cell splitting is performed until the mean free path is resolved or a certain number of
  ! particles is reached. Afterwards a nearest neighbour search for the collision partner is performed.
  ! Source: Pfeiffer, M., Mirza, A. and Fasoulas, S. (2013). A grid-independent particle pairing strategy for DSMC.
  ! Journal of Computational Physics 246, 28–36. doi:10.1016/j.jcp.2013.03.018
  !-----------------------------------------------------------------------------------------------------------------------------------
  DSMC%UseOctree = GETLOGICAL('Particles-DSMC-UseOctree')
  IF(DSMC%ReservoirSimu.AND.DSMC%UseOctree) CALL abort(__STAMP__,'Particles-DSMC-UseOctree = T not allowed for RESERVOIR simulations!')
  DSMC%UseNearestNeighbour = GETLOGICAL('Particles-DSMC-UseNearestNeighbour')
  IF(DSMC%ReservoirSimu.AND.DSMC%UseNearestNeighbour.AND.(DoRestart.OR.(NINT(TEnd/ManualTimeStep).GT.1))) THEN
    CALL abort(__STAMP__,'Particles-DSMC-UseNearestNeighbour = T not allowed for RESERVOIR simulations, if you simulate more than one time step!')
  END IF
  IF(DSMC%UseOctree) THEN
    DO iSpec = 1, nSpecies
      DO iInit = 1, Species(iSpec)%NumberOfInits
        IF (TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'point') THEN
          CALL abort(__STAMP__,'ERROR: No combination of octree and SpaceIC=point possible!')
        END IF
      END DO
    END DO
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
  DSMC%MergeSubcells = GETLOGICAL('Particles-DSMC-MergeSubcells')
  IF(DSMC%MergeSubcells.AND.(Symmetry%Order.NE.2)) THEN
    CALL abort(__STAMP__&
        ,'ERROR: Merging of subcells only supported within a 2D/axisymmetric simulation!')
  END IF
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! Background gas: Check compatibility with other features
  !-----------------------------------------------------------------------------------------------------------------------------------
  IF (BGGas%NumberOfSpecies.GT.0) THEN
    IF (DSMC%UseOctree) THEN
      CALL abort(__STAMP__,'ERROR: Utilization of the octree and nearest neighbour scheme not possible with the background gas!')
    END IF
    DO iSpec = 1, nSpecies
      IF(BGGas%BackgroundSpecies(iSpec)) THEN
        IF(Species(iSpec)%InterID.EQ.4) CALL abort(__STAMP__,'ERROR in BGGas: Electrons as background gas are not yet available!')
      END IF
    END DO
  ELSE
    IF(usevMPF.AND..NOT.(DoRadialWeighting.OR.DoLinearWeighting.OR.DoCellLocalWeighting)) &
      CALL abort(__STAMP__,'ERROR in DSMC: Variable weighting factors are only available with a background gas!')
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
      CALL abort(__STAMP__,'ERROR: Particles-DSMC-alpha has to be in the range between 0 and 1')
    END IF
    DO iSpec = 1, nSpecies
      IF(.NOT.((Species(iSpec)%InterID.EQ.2).OR.(Species(iSpec)%InterID.EQ.20))) CYCLE
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
END IF !CollisMode.GT.0

! If field ionization is used without chemical reactions due to collisions (DSMC chemistry)
IF(DoFieldIonization.AND.(CollisMode.NE.3))THEN
  DO iSpec = 1, nSpecies
    WRITE(UNIT=hilf,FMT='(I0)') iSpec
    ! Heat of formation of ionized species is modified with the ionization energy directly from read-in electronic energy levels
    ! of the ground/previous state of the respective species (Input requires a species number (eg species number of N for NIon1))
    IF((Species(iSpec)%InterID.EQ.10).OR.(Species(iSpec)%InterID.EQ.20)) THEN
      SpecDSMC(iSpec)%PreviousState = GETINT('Part-Species'//TRIM(hilf)//'-PreviousState')
    ELSE
      SpecDSMC(iSpec)%PreviousState = 0
    END IF ! (Species(iSpec)%InterID.EQ.10).OR.(Species(iSpec)%InterID.EQ.20)

    ! Read-in of species for field ionization (only required if it cannot be determined automatically)
    IF ((Species(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
      SpecDSMC(iSpec)%NextIonizationSpecies = GETINT('Part-Species'//TRIM(hilf)//'-NextIonizationSpecies')
    ELSE
      SpecDSMC(iSpec)%NextIonizationSpecies = 0
    END IF
  END DO ! iSpec = 1, nSpecies

  ! Set "NextIonizationSpecies" information for field ionization from "PreviousState" info
  ! NextIonizationSpecies => SpeciesID of the next higher ionization level
  CALL SetNextIonizationSpecies()
END IF

IF (DSMC%ElectronicModel.EQ.4) THEN
  DO iSpec=1, nSpecies
    ALLOCATE(SpecDSMC(iSpec)%CollFreqPreFactor(nSpecies))
    DO jSpec=1, nSpecies
      IF (iSpec.EQ.jSpec) THEN
        delta_ij = 1.0
      ELSE
        delta_ij = 1.0
      END IF
      SpecDSMC(iSpec)%CollFreqPreFactor(jSpec)= 4.*(2.-delta_ij)*CollInf%dref(iSpec,jSpec)**2.0 &
          * SQRT(Pi*BoltzmannConst*CollInf%Tref(iSpec,jSpec)*(Species(iSpec)%MassIC + Species(jSpec)%MassIC) &
          /(2.*(Species(iSpec)%MassIC * Species(jSpec)%MassIC)))/CollInf%Tref(iSpec,jSpec)**(-CollInf%omega(iSpec,jSpec) +0.5)
    END DO
  END DO
END IF

LBWRITE(UNIT_stdOut,'(A)')' INIT DSMC DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitDSMC

SUBROUTINE ReadAHOLevelsFromCSV(FileNameAHO,iSpec)
!===================================================================================================================================
! ATTENTION: needs to be defined on equidistant data-points as .csv
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars                 ,ONLY: AHO
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255), INTENT(IN)    :: FileNameAHO
INTEGER, INTENT(IN)               :: iSpec
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: ioUnit, err, ncounts
REAL                  :: dummy
!===================================================================================================================================
! Read from csv file
OPEN(NEWUNIT=ioUnit,FILE=FileNameAHO,STATUS='OLD')

! Count number of columns
err = 0
ncounts = 0
DO WHILE (err.EQ.0)
  READ(ioUnit,*,IOSTAT = err) dummy
  IF (err.EQ.-1) THEN
    EXIT
  END IF
  err = 0
  ncounts = ncounts + 1
END DO

! Set number of vib energy levels for this species
AHO%NumVibLevels(iSpec) = ncounts

CLOSE (ioUnit)

END SUBROUTINE ReadAHOLevelsFromCSV


SUBROUTINE ReadAHOEnergiesFromCSV(FileNameAHO,iSpec)
!===================================================================================================================================
! ATTENTION: needs to be defined on equidistant data-points as .csv
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars                 ,ONLY: AHO
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255), INTENT(IN)    :: FileNameAHO
INTEGER, INTENT(IN)               :: iSpec
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: ioUnit, columns
!===================================================================================================================================
! Read from csv file
OPEN(NEWUNIT=ioUnit,FILE=FileNameAHO,STATUS='OLD')

! Write to array
DO columns = 1, AHO%NumVibLevels(iSpec)
  read(ioUnit,*) AHO%VibEnergy(iSpec,columns)
END DO
CLOSE (ioUnit)

END SUBROUTINE ReadAHOEnergiesFromCSV


SUBROUTINE SetElectronicModel(iSpec)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals              ,ONLY: abort
USE MOD_DSMC_Vars            ,ONLY: SpecDSMC
USE MOD_Particle_Vars        ,ONLY: Species
USE MOD_DSMC_ElectronicModel ,ONLY: ReadSpeciesLevel
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: iSpec
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(Species(iSpec)%Name.EQ.'none') CALL Abort(__STAMP__,&
    "Read-in from electronic database requires the definition of species name! Species:",IntInfoOpt=iSpec)
IF(.NOT.SpecDSMC(iSpec)%FullyIonized) CALL ReadSpeciesLevel(Species(iSpec)%Name,iSpec)
END SUBROUTINE SetElectronicModel


SUBROUTINE CalcHeatOfFormationIons()
!===================================================================================================================================
! Calculating the heat of formation for ionized species (including higher ionization levels)
! Requires the completed read-in of species data
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_ReadInTools
USE MOD_Globals          ,ONLY: abort,UNIT_stdOut
#if USE_MPI
USE MOD_Globals          ,ONLY: mpiroot
#endif
USE MOD_Globals_Vars     ,ONLY: BoltzmannConst,Joule2eV
USE MOD_PARTICLE_Vars    ,ONLY: nSpecies, Species, SpeciesDatabase
USE MOD_DSMC_Vars        ,ONLY: SpecDSMC
USE MOD_io_hdf5
USE MOD_HDF5_input       ,ONLY:ReadAttribute
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(32)         :: hilf2
LOGICAL               :: AutoDetect
INTEGER               :: iSpec,jSpec,counter,MaxElecQua,err
CHARACTER(LEN=64)     :: dsetname
INTEGER(HID_T)        :: file_id_specdb                       ! File identifier
!===================================================================================================================================
AutoDetect=.TRUE.
DO iSpec = 1, nSpecies
  counter = 0
  IF((Species(iSpec)%InterID.EQ.10).OR.(Species(iSpec)%InterID.EQ.20)) THEN
    IF(SpecDSMC(iSpec)%PreviousState.EQ.0) THEN
      ! Fall-back for "fake" ions (e.g. the previous state is not available)
      IF(SpeciesDatabase.NE.'none') THEN
        ! Initialize FORTRAN interface.
        CALL H5OPEN_F(err)
        CALL H5FOPEN_F (TRIM(SpeciesDatabase), H5F_ACC_RDONLY_F, file_id_specdb, err)
        dsetname = TRIM('/Species/'//TRIM(Species(iSpec)%Name))
        CALL ReadAttribute(file_id_specdb,'HeatOfFormation_K',1,DatasetName = dsetname,RealScalar=SpecDSMC(iSpec)%HeatOfFormation,ReadFromGroup=.TRUE.)
        SpecDSMC(iSpec)%HeatOfFormation = SpecDSMC(iSpec)%HeatOfFormation * BoltzmannConst
        CALL PrintOption('HeatOfFormation [J], '//TRIM(Species(iSpec)%Name),'DB',RealOpt=SpecDSMC(iSpec)%HeatOfFormation)
        ! Close the file.
        CALL H5FCLOSE_F(file_id_specdb, err)
        ! Close FORTRAN interface.
        CALL H5CLOSE_F(err)
      ELSE
        WRITE(UNIT=hilf2,FMT='(I0)') iSpec
        SpecDSMC(iSpec)%HeatOfFormation = GETREAL('Part-Species'//TRIM(hilf2)//'-HeatOfFormation_K') * BoltzmannConst
      END IF
    ELSE
      ! Get the last electronic state level from the previous state (e.g. for NIon get the last electronic state for N)
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
          IF(counter.GT.100) CALL abort(__STAMP__,&
              'ERROR: Nbr. of ionization lvls per spec limited to 100. More likely wrong input in PreviuosState of spec:', iSpec)
        END DO
        IF(AutoDetect)THEN
          LBWRITE(UNIT_stdOut,'(A)')' | Automatically determined HeatOfFormation:'
          AutoDetect=.FALSE.
        END IF
        ! Add the heat of formation of the ground state
        SpecDSMC(iSpec)%HeatOfFormation = SpecDSMC(iSpec)%HeatOfFormation + SpecDSMC(jSpec)%HeatOfFormation
        WRITE(UNIT=hilf2,FMT='(I0)') iSpec
        CALL PrintOption('Part-Species'//TRIM(hilf2)//'-HeatOfFormation_K  [K]','CALCUL.',&
            RealOpt=SpecDSMC(iSpec)%HeatOfFormation/BoltzmannConst)
        CALL PrintOption('converted to [eV]','CALCUL.',&
            RealOpt=SpecDSMC(iSpec)%HeatOfFormation*Joule2eV)
      ELSE
        CALL abort(__STAMP__,'Chemical reactions with ionized species require an input of electronic energy level(s)!', iSpec)
      END IF
    END IF
  END IF
END DO

END SUBROUTINE CalcHeatOfFormationIons


SUBROUTINE SetNextIonizationSpecies()
!===================================================================================================================================
! Set "NextIonizationSpecies" information for field ionization from "PreviousState" info
! NextIonizationSpecies => SpeciesID of the next higher ionization level
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals          ,ONLY: UNIT_stdOut
#if USE_MPI
USE MOD_Globals          ,ONLY: mpiroot
#endif
USE MOD_PARTICLE_Vars    ,ONLY: nSpecies, Species
USE MOD_DSMC_Vars        ,ONLY: SpecDSMC
USE MOD_ReadInTools      ,ONLY: PrintOption
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
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
  IF((Species(iSpec)%InterID.NE.4).AND.(Species(iSpec)%InterID.NE.100)) THEN
    IF(SpecDSMC(iSpec)%PreviousState.NE.0)THEN
      SpecDSMC(SpecDSMC(iSpec)%PreviousState)%NextIonizationSpecies = iSpec
      AutoDetect=.TRUE.
    END IF
  ELSE
    SpecDSMC(iSpec)%NextIonizationSpecies = 0
  END IF
END DO
IF(AutoDetect)THEN
  LBWRITE(UNIT_stdOut,'(A)')' | Automatically determined NextIonizationSpecies:'
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
USE MOD_Globals                ,ONLY: abort, IK, MPI_COMM_PICLAS
USE MOD_PARTICLE_Vars          ,ONLY: nSpecies, Species
USE MOD_Restart_Vars           ,ONLY: DoRestart,RestartFile
USE MOD_Particle_Vars          ,ONLY: nSpecies, PartSpecies, Species
USE MOD_part_tools             ,ONLY: GetParticleWeight
USE MOD_HDF5_INPUT             ,ONLY: DatasetExists,ReadAttribute,ReadArray
USE MOD_IO_HDF5
USE MOD_DSMC_Vars              ,ONLY: VarVibRelaxProb, CollInf, Coll_pData, DSMC
USE MOD_Mesh_Vars              ,ONLY: nElems, offsetElem
USE MOD_DSMC_Analyze           ,ONLY: CalcInstantTransTemp
USE MOD_Particle_Vars          ,ONLY: PEM
USE MOD_DSMC_Relaxation        ,ONLY: DSMC_calc_var_P_vib
USE MOD_part_emission_tools    ,ONLY: CalcVelocity_maxwell_lpn
#if USE_MPI
USE MOD_Globals                ,ONLY: MPIRoot
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSpec, jSpec, iPart, iElem, dim, n
REAL                  :: VibProb, Ti, Tj, CRela2, Velo1(3), Velo2(3), MPF
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
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)
  ! read local ParticleInfo from HDF5
  CALL DatasetExists(File_ID,'VibProbInfo',VibProbDataExists)
  IF(VibProbDataExists)THEN
    VibProbInitDone = .TRUE.
    ! Associate construct for integer KIND=8 possibility
    LBWRITE(*,*) 'Set variable vibrational relaxation probability from restart file'
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
    LBWRITE(*,*) 'Set uniform vibrational relaxation probability from restart file'
    DO iElem = 1, nElems
      DO iSpec = 1, nSpecies
        VarVibRelaxProb%ProbVibAv(iElem,iSpec) = VibProb
      END DO
    END DO
  END IF ! If 'VibProbConstInfo' exists
  IF(.NOT.VibProbInitDone) THEN
    ! Set vibrational relaxation probability to default value
    LBWRITE(*,*) 'No vibrational relaxation probability data in restart file\n', &
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
  LBWRITE(*,*) '| Set vibrational relaxation probability based on temperature in the cell'
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
      MPF = GetParticleWeight(iPartIndx(iPart))
      CollInf%Coll_SpecPartNum(PartSpecies(iPartIndx(iPart))) = CollInf%Coll_SpecPartNum(PartSpecies(iPartIndx(iPart))) + MPF
      nPerSpec(PartSpecies(iPartIndx(iPart))) = nPerSpec(PartSpecies(iPartIndx(iPart))) + 1
    END DO
    CALL CalcInstantTransTemp(iPartIndx,nPart)
    DO iSpec = 1, nSpecies
      IF(.NOT.((Species(iSpec)%InterID.EQ.2).OR.(Species(iSpec)%InterID.EQ.20)))CYCLE
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


SUBROUTINE ReadRotationalSpeciesLevel (iSpec)
!===================================================================================================================================
! Subroutine to read the rotational levels from SpeciesDatabase.h5
!===================================================================================================================================
! use module
USE MOD_io_hdf5
USE MOD_Globals
USE MOD_DSMC_Vars             ,ONLY: SpecDSMC
USE MOD_HDF5_Input            ,ONLY: DatasetExists
USE MOD_Particle_Vars         ,ONLY: Species, SpeciesDatabase
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars      ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                    :: iSpec
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                               :: err
! HDF5 specifier taken from extractParticles
CHARACTER(LEN=256)                                    :: ElLevelDatabase
CHARACTER(LEN=64)                                     :: datasetname
INTEGER(HSIZE_T), DIMENSION(2)                        :: dims,sizeMax
INTEGER(HID_T)                                        :: file_id_dsmc                       ! File identifier
INTEGER(HID_T)                                        :: dset_id_dsmc                       ! Dataset identifier
INTEGER(HID_T)                                        :: filespace                          ! filespace identifier
REAL,ALLOCATABLE                                      :: RotationalState(:,:)
LOGICAL                                               :: DataSetFound
!===================================================================================================================================
datasetname = TRIM('/Species/'//TRIM(Species(iSpec)%Name)//'/RotationalLevel')
ElLevelDatabase = TRIM(SpeciesDatabase)
LBWRITE(UNIT_StdOut,'(A)') ' | Read rotational level entries '//TRIM(datasetname)//' from '//TRIM(ElLevelDatabase)
! Initialize FORTRAN interface.
CALL H5OPEN_F(err)
! Open the file.
CALL H5FOPEN_F (TRIM(ElLevelDatabase), H5F_ACC_RDONLY_F, file_id_dsmc, err)
CALL DatasetExists(File_ID_DSMC,TRIM(datasetname),DataSetFound)
IF(.NOT.DataSetFound)THEN
  CALL abort(__STAMP__,'DataSet not found: ['//TRIM(datasetname)//'] ['//TRIM(ElLevelDatabase)//']')
END IF
! Open the  dataset.
CALL H5DOPEN_F(file_id_dsmc, datasetname, dset_id_dsmc, err)
! Get the file space of the dataset.
CALL H5DGET_SPACE_F(dset_id_dsmc, FileSpace, err)
! get size
CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, dims, SizeMax, err)
! Allocate rotational_state
ALLOCATE (RotationalState( 1:dims(1), 0:dims(2)-1 ) )
! read data
CALL H5dread_f(dset_id_dsmc, H5T_NATIVE_DOUBLE, RotationalState, dims, err)

ALLOCATE ( SpecDSMC(iSpec)%RotationalState( 1:dims(1), 0:dims(2)-1 ) )
SpecDSMC(iSpec)%RotationalState = RotationalState
SpecDSMC(iSpec)%MaxRotQuant  = SIZE(SpecDSMC(iSpec)%RotationalState,2)
! Close the file.
CALL H5FCLOSE_F(file_id_dsmc, err)
! Close FORTRAN interface.
CALL H5CLOSE_F(err)

END SUBROUTINE ReadRotationalSpeciesLevel


SUBROUTINE ReadVibrationalSpeciesLevel (iSpec)
!===================================================================================================================================
! Subroutine to read the vibrational levels from SpeciesDatabase.h5
!===================================================================================================================================
! use module
USE MOD_io_hdf5
USE MOD_Globals
USE MOD_DSMC_Vars             ,ONLY: AHO
USE MOD_HDF5_Input            ,ONLY: DatasetExists
USE MOD_Particle_Vars         ,ONLY: Species, SpeciesDatabase, nSpecies
USE MOD_Array_Operations      ,ONLY: ChangeSizeArray
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars      ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                    :: iSpec
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                               :: err
! HDF5 specifier taken from extractParticles
CHARACTER(LEN=256)                                    :: LevelDatabase
CHARACTER(LEN=64)                                     :: datasetname
INTEGER(HSIZE_T), DIMENSION(2)                        :: dims,sizeMax
INTEGER(HID_T)                                        :: file_id_dsmc                       ! File identifier
INTEGER(HID_T)                                        :: dset_id_dsmc                       ! Dataset identifier
INTEGER(HID_T)                                        :: filespace                          ! filespace identifier
LOGICAL                                               :: DataSetFound
!===================================================================================================================================
datasetname = TRIM('/Species/'//TRIM(Species(iSpec)%Name)//'/VibrationalLevel')
LevelDatabase = TRIM(SpeciesDatabase)
LBWRITE(UNIT_StdOut,'(A)') ' | Read vibrational level entries '//TRIM(datasetname)//' from '//TRIM(LevelDatabase)
! Initialize FORTRAN interface.
CALL H5OPEN_F(err)
! Open the file.
CALL H5FOPEN_F (TRIM(LevelDatabase), H5F_ACC_RDONLY_F, file_id_dsmc, err)
CALL DatasetExists(File_ID_DSMC,TRIM(datasetname),DataSetFound)
IF(.NOT.DataSetFound)THEN
  CALL abort(__STAMP__,'DataSet not found: ['//TRIM(datasetname)//'] ['//TRIM(LevelDatabase)//']')
END IF
! Open the  dataset.
CALL H5DOPEN_F(file_id_dsmc, datasetname, dset_id_dsmc, err)
! Get the file space of the dataset.
CALL H5DGET_SPACE_F(dset_id_dsmc, FileSpace, err)
! get size
CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, dims, SizeMax, err)
! Set number of vibrational levels
AHO%NumVibLevels(iSpec) = INT(dims(2),4)
IF(.NOT.ALLOCATED(AHO%VibEnergy))THEN
  ALLOCATE(AHO%VibEnergy(nSpecies,MAXVAL(AHO%NumVibLevels)))
  AHO%VibEnergy = 0.
END IF
! Adapt allocation if vibrational levels of ispec is larger than the current maximum
IF(AHO%NumVibLevels(iSpec).GT.MAXVAL(AHO%NumVibLevels))THEN
  CALL ChangeSizeArray(AHO%VibEnergy,MaxVal(AHO%NumVibLevels(1:iSpec-1)),AHO%NumVibLevels(iSpec))
END IF

! read data
CALL H5dread_f(dset_id_dsmc, H5T_NATIVE_DOUBLE, AHO%VibEnergy(iSpec, 1:AHO%NumVibLevels(iSpec)), dims, err)

! Close the file.
CALL H5FCLOSE_F(file_id_dsmc, err)
! Close FORTRAN interface.
CALL H5CLOSE_F(err)

END SUBROUTINE ReadVibrationalSpeciesLevel


SUBROUTINE FinalizeDSMC()
!----------------------------------------------------------------------------------------------------------------------------------!
! finalize dsmc variables
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_DSMC_Vars
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
IF(DSMC%VibAHO) THEN
  SDEALLOCATE(AHO%NumVibLevels)
  SDEALLOCATE(AHO%VibEnergy)
  SDEALLOCATE(AHO%omegaE)
  SDEALLOCATE(AHO%chiE)
END IF
SDEALLOCATE(DSMC%QualityFacSampRot)
SDEALLOCATE(DSMC%QualityFacSampVib)
SDEALLOCATE(DSMC%QualityFacSampRotSamp)
SDEALLOCATE(DSMC%QualityFacSampVibSamp)
SDEALLOCATE(DSMC%CalcVibProb)
SDEALLOCATE(DSMC%CalcRotProb)
SDEALLOCATE(DSMC%InstantTXiElec)
SDEALLOCATE(SampDSMC)
SDEALLOCATE(PartStateIntEn)
SDEALLOCATE(ElecRelaxPart)
SDEALLOCATE(SpecDSMC)
IF(DSMC%NumPolyatomMolecs.GT.0) THEN
  DO iPoly=1,DSMC%NumPolyatomMolecs
    SDEALLOCATE(PolyatomMolDSMC(iPoly)%CharaTVibDOF)
    SDEALLOCATE(PolyatomMolDSMC(iPoly)%LastVibQuantNums)
    SDEALLOCATE(PolyatomMolDSMC(iPoly)%MaxVibQuantDOF)
    SDEALLOCATE(PolyatomMolDSMC(iPoly)%GammaVib)
    SDEALLOCATE(PolyatomMolDSMC(iPoly)%VibRelaxProb)
    SDEALLOCATE(PolyatomMolDSMC(iPoly)%CharaTRotDOF)
    SDEALLOCATE(PolyatomMolDSMC(iPoly)%MomentOfInertia)
  END DO
  SDEALLOCATE(PolyatomMolDSMC)
END IF
SDEALLOCATE(DSMC%NumColl)
SDEALLOCATE(DSMC%InstantTransTemp)
IF(DSMC%CalcQualityFactors) THEN
  SDEALLOCATE(DSMC%QualityFacSamp)
END IF
SDEALLOCATE(Coll_pData)
SDEALLOCATE(SampDSMC)
SDEALLOCATE(MacroDSMC)
SDEALLOCATE(QKChemistry)

SDEALLOCATE(ChemReac%ReactModel)
SDEALLOCATE(ChemReac%BackwardReacForwardIndx)
SDEALLOCATE(ChemReac%BackwardReac)
SDEALLOCATE(ChemReac%QKRColl)
SDEALLOCATE(ChemReac%QKTCollCorrFac)
SDEALLOCATE(ChemReac%NumReac)
SDEALLOCATE(ChemReac%ReacCount)
SDEALLOCATE(ChemReac%ReacCollMean)
SDEALLOCATE(ChemReac%NumReac)
SDEALLOCATE(ChemReac%ReactType)
SDEALLOCATE(ChemReac%ReactionName)
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
SDEALLOCATE(ChemReac%totalReacToModel)

SDEALLOCATE(CollInf%collidingSpecies)
SDEALLOCATE(CollInf%Coll_Case)
SDEALLOCATE(CollInf%Coll_CaseNum)
SDEALLOCATE(CollInf%Coll_SpecPartNum)
SDEALLOCATE(CollInf%Cab)
SDEALLOCATE(CollInf%KronDelta)
SDEALLOCATE(CollInf%FracMassCent)
SDEALLOCATE(CollInf%MassRed)
SDEALLOCATE(CollInf%SumPairMPF)
SDEALLOCATE(CollInf%alphaVSS)
SDEALLOCATE(CollInf%omega)
SDEALLOCATE(CollInf%dref)
SDEALLOCATE(CollInf%Tref)
SDEALLOCATE(CollInf%OldCollPartner)
CollInf%ProhibitDoubleColl=.FALSE.
!SDEALLOCATE(VibQuantsPar)
! SDEALLOCATE(XiEq_Surf)
SDEALLOCATE(DSMC_Solution)
SDEALLOCATE(DSMC_SolutionPressTens)
CALL DeleteElemNodeVol()

SDEALLOCATE(BGGas%PairingPartner)
SDEALLOCATE(BGGas%BackgroundSpecies)
SDEALLOCATE(BGGas%TraceSpecies)
SDEALLOCATE(BGGas%MapSpecToBGSpec)
SDEALLOCATE(BGGas%MapBGSpecToSpec)
SDEALLOCATE(BGGas%SpeciesFraction)
SDEALLOCATE(BGGas%SpeciesFractionElem)
SDEALLOCATE(BGGas%NumberDensity)
SDEALLOCATE(BGGas%DistributionSpeciesIndex)
SDEALLOCATE(BGGas%Distribution)
SDEALLOCATE(BGGas%DistributionNumDens)
SDEALLOCATE(BGGas%Region)
SDEALLOCATE(BGGas%RegionElemType)
#ifdef drift_diffusion
SDEALLOCATE(BGGas%ReducedTownsendCoefficient)
SDEALLOCATE(BGGas%DriftDiffusionCoefficient)
SDEALLOCATE(BGGas%ElectronMobility)
#endif

SDEALLOCATE(ParticleWeighting%ClonePartNum)
SDEALLOCATE(ParticleWeighting%PartInsSide)
SDEALLOCATE(LinearWeighting%ScalePoint)
SDEALLOCATE(LinearWeighting%VarMPF)
SDEALLOCATE(ClonedParticles)
SDEALLOCATE(AmbiPolarSFMapping)
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
USE MOD_Symmetry_Vars         ,ONLY: Symmetry
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
TYPE (tNodeVolume), INTENT(INOUT)  :: Node
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER     ::  iLoop, nLoop
!===================================================================================================================================
nLoop = 2**Symmetry%Order
IF(ASSOCIATED(Node%SubNode)) THEN
  DO iLoop = 1, nLoop
    CALL DeleteNodeVolume(Node%SubNode(iLoop))
  END DO
  DEALLOCATE(Node%SubNode)
END IF

END SUBROUTINE DeleteNodeVolume


END MODULE MOD_DSMC_Init
