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

INTERFACE DSMC_SetInternalEnr_LauxVFD
  MODULE PROCEDURE DSMC_SetInternalEnr_LauxVFD
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: InitDSMC, DSMC_SetInternalEnr_LauxVFD, FinalizeDSMC
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

CALL prms%CreateLogicalOption(  'Particles-DSMC-OutputMeshInit'      &
                                        ,  'not working currently | Writeoutput mesh for constant pressure BC at initialization.'&
                                        , '.FALSE.')

CALL prms%CreateLogicalOption(  'Particles-DSMC-OutputMeshSamp'      &
                                        , 'not working currently | Write output mesh for constant pressure BC with sampling'//&
                                          'values at t_analyze.' , '.FALSE.')

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
                                              '(HALOWIKI:)Choice of the vibrational relaxation probability calculation\n'//&
                                          '0-1: constant\n'//&
                                          '2: variable, Boyd)', '0.2')
CALL prms%CreateRealOption(     'Particles-DSMC-VibRelaxProb'&
                                          , 'Define the vibrational relaxation probability upon collision of molecules', '0.004')
CALL prms%CreateRealOption(     'Particles-DSMC-ElecRelaxProb'&
                                          , 'Define the elextronic relaxation probability upon collision of molecules', '0.01')
CALL prms%CreateRealOption(     'Particles-DSMC-GammaQuant'&
                                          , 'Set the GammaQuant for zero point energy in Evib (perhaps also Erot) should be'//&
                                          ' 0.5 or 0.', '0.5')
CALL prms%CreateLogicalOption(  'Particles-DSMC-BackwardReacRate'&
                                         , 'Set [TRUE] to enable the automatic calculation of the backward reaction rate '//&
                                           'coefficientusing the equilibrium constant calculated by partition functions\n'//&
                                           '[FALSE] if they are defined as separate reactions.' , '.FALSE.')
CALL prms%CreateRealOption(     'Particles-DSMC-PartitionMaxTemp'&
                                          , 'Define temperature limit for pre-stored partition function that are used for '//&
                                          'calculation of backwards rates', '20000.0')
CALL prms%CreateRealOption(     'Particles-DSMC-PartitionInterval'&
                                          , 'Define temperature interval for pre-stored partition functions that are used for '//&
                                          'calculation of backwards rates', '10.0')
CALL prms%CreateRealOption(     'Particles-DSMC-veloMinColl-Spec[$]' , 'min velo magn. for spec allowed to perform collision' , '0.' &
                                          , numberedmulti=.TRUE.)
!-----------------------------------------------------------------------------------
CALL prms%CreateLogicalOption(  'Particles-DSMC-CalcQualityFactors'&
                                          , 'Enables [TRUE] / disables [FALSE] the calculation and output of flow-field variable.\n'//&
                                           'Maximal collision probability\n'//&
                                          'Time-averaged mean collision probability\n'//&
                                          'Mean collision separation distance over mean free path ' , '.FALSE.')
CALL prms%CreateLogicalOption(  'Particles-DSMCReservoirSim'&
                                            , 'Only TD=Reservoir (42).\n'//&
                                          'Set [TRUE] to disable particle movement. Use for reservoir simulations.' , '.FALSE.')
CALL prms%CreateLogicalOption(  'Particles-DSMCReservoirSimRate'&
                                          , 'Only TD=Reservoir (42).\n'//&
                                          'Set [TRUE] to disable particle reactions.Only probabilities (rates) are calculated.' &
                                        , '.FALSE.')
CALL prms%CreateLogicalOption(  'Particles-DSMCReservoirStatistic'&
                                         , 'Only TD=Reservoir (42).\n'//&
                                          'Probabilities (rates) are calculated\n'//&
                                          ' [TRUE] counting reacting particles.\n'//&
                                          ' [FALSE] summing reaction probabilities (does not work with Q-K).' , '.FALSE.')
CALL prms%CreateLogicalOption(  'Particles-DSMCReservoirSurfaceRate'&
                                          , 'Only TD=Reservoir (42).\n'//&
                                          'Set [TRUE] to disable particle adsorption and desorption and keep surface coverage '//&
                                            'constant. Only probabilities (rates) are calculated.' , '.FALSE.')
CALL prms%CreateIntOption(      'Particles-ModelForVibrationEnergy'&
                                          , 'Define model used for vibrational degrees of freedom.\n'//&
                                          ' 0: SHO  simple harmonic oscillator \n'//&
                                          ' 1: TSHO truncated simple harmonic oscillator .', '0')
CALL prms%CreateLogicalOption(  'Particles-DSMC-TEVR-Relaxation'&
                                          , 'Flag for Translational-Vibrational-Electric-Rotational relaxation (T-V-E-R)\n'//&
                                          '[TRUE] or more simple T-V-R T-E-R\n'//&
                                          '[FALSE] relaxation.' , '.FALSE.')
CALL prms%CreateLogicalOption(  'Particles-DSMC-ElectronicModel'&
                                          , 'Set [TRUE] to model electronic states of atoms and molecules.' , '.FALSE.')
CALL prms%CreateStringOption(   'Particles-DSMCElectronicDatabase'&
                                          , 'If electronic model is used give (relative) path to (h5) Name of Electronic State'//&
                                          ' Database', 'none')
CALL prms%CreateRealOption(     'EpsMergeElectronicState'&
                                         , 'Percentage parameter of electronic energy level merging.' , '1E-4')
CALL prms%CreateLogicalOption(  'Particles-DSMC-UseQCrit'&
                                         , 'Set [TRUE] to enable steady state detection and sampling start using Q-criterion'//&
                                           ' (Burt/Boyd).', '.FALSE.')
CALL prms%CreateLogicalOption(  'Particles-DSMC-UseSSD'&
                                         , 'Set [TRUE] to enable steady state detection and sampling start using 3SD routines.' &
                                         , '.FALSE.')
CALL prms%CreateIntOption(      'Particles-DSMCBackgroundGas'&
                                         , 'Define Species number that is used as background gas species', '0')
CALL prms%CreateRealOption(     'Particles-DSMCBackgroundGasDensity'&
                                         , 'Define Species number density for background gas', '0.')
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
                                         ,'Allows to enable/disable the nearest neighbour search algorithm within the ocrtree '//&
                                          'cell refinement','.TRUE.')
CALL prms%CreateLogicalOption(  'Particles-DSMC-ProhibitDoubleCollisions'&
                                         ,'2D/Axisymmetric only: Prohibit the occurrence of repeated collisions between the '//&
                                          'same particle pairs in order to reduce the statistical dependence')
CALL prms%CreateLogicalOption(  'Particles-DSMC-MergeSubcells'&
                                         ,'2D/Axisymmetric only: Merge subcells divided by the quadtree algorithm to satisfy '//&
                                          'the minimum particle per subcell requirement', '.FALSE.')

CALL prms%SetSection("DSMC Collision")
CALL prms%CreateIntOption(      'crossSectionConstantMode'  &
                                            ,' Flags which cross section Cab Mode is used. sigma=Cab * cr^(-2 omega).\n ' //&
                                            ' 0 : single omega for the computational domain(Bird 1981) -'                 //&
                                            '     Part-Species1=temperatureExponentOmega will be set for all and \n'      //&
                                            '     Cab will be calculated via species-specific factor A_j \n '             //&
                                            '     Cab=(A_1+A_2)^2*m_red^-omegaLaux'                                       //&
                                            ' 1 : Cab will be calculated via species-specific factor A_j \n '             //&
                                            ' 2 : Cab will be calculated directly (see bird1981 eq (9),constant values)', '0')
CALL prms%CreateLogicalOption(   'averagedCollisionParameters'  &
                                           ,' Flags if collision parameters are specific to colliding species'        //&
                                            ' or species specific and averaged for the collision itself.'             //&
                                            ' T: Part-Species[$]-temperatureExponentOmega,-referenceTemperature,'     //&
                                            '    -referenceDiameter,-alphaVSS to be set.\n '                          //&
                                            '    species-specific parameters(T) can be found in tables e.g. in\n'     //& 
                                            '    VHS/VSS bird1994 VHS: table A1 and A2/VSS: table A1 and A3'          //&
                                            '    VSS     weaver2014(https://doi.org/10.1063/1.4921245)'               //&
                                            ' F: Part-Collision[$]-[$]-temperatureExponentOmega,-referenceTemperature,'  //&
                                            '    -referenceDiameter,-alphaVSS to be set.\n'                              //&
                                            '    collision-specific parameters(F) can be found in tables e.g. in\n'      //& 
                                            '    VHS/VSS krishnan2015(https://doi.org/10.2514/6.2015-3373),\n'           //&
                                            '    VHS/VSS krishnan2016(https://doi.org/10.1063/1.4939719)', 'T')
CALL prms%CreateRealOption(     'Part-Collision[$]-[$]-referenceTemperature'  &
                                           ,' collision parameter: collision-specific reference temperature for VHS/VSS model.'//&
                                            ' (index=speciesi-speciesj)' , '0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Collision[$]-[$]-referenceDiameter'  &
                                           ,' collision parameter: collision-specific reference diameter for VHS/VSS model. '&
                                           , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Collision[$]-[$]-temperatureExponentOmega'  &
                                           ,' collision parameter: collision-specific temperature exponent omega=2/(eta-1) '//&
                                            ' for VHS/VSS model. default omega=0 plus default alpha=1 reproduces HS model ' //&
                                            ' CAUTION: temperatureExponentOmega=omega_bird1994-0.5','0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Collision[$]-[$]-scatteringExponentAlphaVSS'  &
                                           ,' collision parameter: collision-specific scattering exponent alpha for VSS model' //&
                                            ' default alphaVSS=1 reproduces VHS model'                                         //&
                                            ' See Bird 1994 p.42 for more information.'                                          &
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
CALL prms%CreateRealOption(     'Part-Species[$]-referenceTemperature'  &
                                           ,' collision parameter: species-specific reference temperature for VHS/VSS model.' &
                                           , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-referenceDiameter' &
                                           ,'collision parameter: species-specific reference diameter for VHS/VSS model.'&
                                           , '0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-temperatureExponentOmega'  &
                                           ,' collision parameter: species-specific temperature exponent omega=2/(eta-1) '//&
                                            ' for VHS/VSS model. CAUTION: temperatureExponentOmega=omega_bird1994-0.5'      &
                                           , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-scatteringExponentalphaVSS'&
                                           ,' collision parameter: species-specific scattering exponent alpha for VSS model'  //&
                                            ' default alphaVSS=1 reproduces VHS model'                                         //&
                                            ' See Bird 1994 p.42 for more information.', '1.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-CharaTempVib','Characteristic vibrational temperature.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-CharaTempRot'  &
                                           ,'Characteristic rotational temperature', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Ediss_eV','Energy of Dissoziation in [eV].', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-VFDPhi3'  &
                                           ,'Factor of Phi3 in VFD Method: Phi3 = 0 => VFD', '0.'&
                                           , numberedmulti=.TRUE.)
! ----------------------------------------------------------------------------------------------------------------------------------
CALL prms%CreateLogicalOption(  'Particles-DSMC-useRelaxProbCorrFactor'&
                                           ,'Use the relaxation probability correction factor of Lumpkin', '.FALSE.')
CALL prms%CreateRealOption(     'Part-Species[$]-CollNumRotInf'  &
                                           ,'Collision number for rotational relaxation according to Parker or'//&
                                            'Zhang, ini_2 -> model dependent!', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-TempRefRot'  &
                                           ,'Referece temperature for rotational relaxation according to Parker or '//&
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
CALL prms%CreateRealOption(     'Part-Species[$]-TempVib'  &
                                           ,'Vibrational temperature.', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-TempRot'  &
                                           ,'Rotational temperature.', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-TempElec'  &
                                           ,'Electronic temperature.', '0.', numberedmulti=.TRUE.)

CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-TempVib'  &
                                           ,'Vibrational temperature.', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-TempRot'  &
                                           ,'Rotational temperature.', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Init[$]-TempElec'  &
                                           ,'Electronic temperature.', '0.', numberedmulti=.TRUE.)

CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-TempVib'  &
                                           ,'Vibrational temperature.', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-TempRot'  &
                                           ,'Rotational temperature.', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-TempElec'  &
                                           ,'Electronic temperature.', '0.', numberedmulti=.TRUE.)

CALL prms%CreateRealOption(     'Part-Species[$]-HeatOfFormation_K'  &
                                           ,'Heat of formation of the respective species [Kelvin]'&
                                           , numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-PreviousState'  &
                                           ,'Species number of the previous state (e.g. N for NIon) ', '0', numberedmulti=.TRUE.)
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

CALL prms%CreateIntOption(      'Part-Species[$]-SymmetryFactor'  &
                                           , 'TODO-DEFINE-PARAMETER', '0', numberedmulti=.TRUE.)

CALL prms%CreateRealOption(     'Part-Species[$]-IonizationEn_eV'  &
                                           ,'Energy of Ionization in [eV].', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-RelPolarizability'  &
                                           ,'Relative Polarizability', '0.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-NumEquivElecOutShell'  &
                                           ,'Number of equivalent electrons in outer shells', '0'&
                                           , numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-NumOfProtons'  &
                                           ,'Number of protons for respective species.', '0', numberedmulti=.TRUE.)


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

CALL prms%SetSection("DSMC Chemistry")
CALL prms%CreateIntOption(      'DSMC-NumOfReactions'  &
                                           ,'Number of reactions.', '0')
CALL prms%CreateIntOption(      'DSMC-Reaction[$]-NumberOfNonReactives'  &
                                           ,'TODO-DEFINE-PARAMETER', '0', numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'DSMC-Reaction[$]-NonReactiveSpecies'  &
                                           ,'Array with the non-reactive collision partners for dissociation'&
                                           ,numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'DSMC-Reaction[$]-ReactionType'  &
                                           ,'Used reaction type\n'//&
                                            'I: electron impact ionization\n'//&
                                            'R: molecular recombination\n'//&
                                            'D: molecular dissociation\n'//&
                                            'E: molecular exchange reaction\n'//&
                                            'X: simple charge exchange reaction)', 'none', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'DSMC-Reaction[$]-QKProcedure'  &
                                           ,'Flag to use quantum-kinetic model', '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'DSMC-Reaction[$]-QK-Method'  &
                                           ,'Recombination Method for Q-K model\n'//&
                                            '1: by Bird\n'//&
                                            '2: by Gallis)\n'//&
                                            'If using bird, define the variables:\n'//&
                                            'DSMC-Reaction[$]-QK-Coeff1\n'//&
                                            'DSMC-Reaction[$]-QK-Coeff2 ', '0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-QK-Coeff1'  &
                                           ,'First Q-K coefficient for Birds method.', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-QK-Coeff2'  &
                                           ,'Second Q-K coefficient for Birds method.', '0.', numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'DSMC-Reaction[$]-Reactants'  &
                                           ,'Reactants of Reaction[$]\n'//&
                                            '(SpecNumOfReactant1,\n'//&
                                            'SpecNumOfReactant2,\n'//&
                                            'SpecNumOfReactant3)', '0 , 0 , 0' , numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'DSMC-Reaction[$]-Products'  &
                                           ,'Products of Reaction[j] (Product1, Product2, Product3)', '0 , 0 , 0' &
                                           , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-Arrhenius-Prefactor'  &
                                           , 'TODO-DEFINE-PARAMETER ', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-Arrhenius-Powerfactor'  &
                                           , 'TODO-DEFINE-PARAMETER', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-Activation-Energy_K'  &
                                           , 'Activation energy (relativ to k_Boltzmann) for Reaction[$].', '0.' &
                                           , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-CEXa'  &
                                , 'CEX log-factor '//&
                                '(g-dep. cross section in Angstrom, def.: value for Xe+)', '-27.2' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-CEXb'  &
                                , 'CEX const. factor '//&
                                '(g-dep. cross section in Angstrom, def.: value for Xe+)', '175.269' , numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'DSMC-Reaction[$]-DoScat'  &
                                , 'Perform scattering-based charge-exchange instead of isotropic '//&
                                '(model of Samuel Araki by lookup table)', '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-ELa'  &
                                , 'with DoScat=T: EL log-factor '//&
                                '(g&cut-off-angle-dep. cs in Angstrom, def.: value for Xe+)', '-26.8' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-ELb'  &
                                , 'with DoScat=T: EL const. factor '//&
                                '(g&cut-off-angle-dep. cs in Angstrom, def.: value for Xe+)', '148.975' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-MEXa'  &
                                , 'with DoScat=F: MEX log-factor '//&
                                '(g-dep. cross section in Angstrom, def.: value for Xe+)', '-27.2' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'DSMC-Reaction[$]-MEXb'  &
                                , 'with DoScat=F: MEX const. factor '//&
                                '(g-dep. cross section in Angstrom, def.: value for Xe+)', '175.269' , numberedmulti=.TRUE.)
CALL prms%CreateStringOption(     'DSMC-Reaction[$]-TLU_FileName'  &
                                , 'with DoScat=F: No TLU-File needed '//&
                                '(def.: )', '0' , numberedmulti=.TRUE.)

END SUBROUTINE DefineParametersDSMC

SUBROUTINE InitDSMC()
!===================================================================================================================================
! Init of DSMC Vars
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc                ,ONLY: PP_N
USE MOD_Mesh_Vars              ,ONLY: nElems, NGEo, SideToElem
USE MOD_Globals_Vars           ,ONLY: Pi, BoltzmannConst, ElementaryCharge
USE MOD_ReadInTools
USE MOD_DSMC_Vars
USE MOD_Particle_Vars          ,ONLY: nSpecies, Species, PDM, PartSpecies, Adaptive_MacroVal, Symmetry2D, VarTimeStep
USE MOD_Particle_Vars          ,ONLY: DoFieldIonization
USE MOD_DSMC_Analyze           ,ONLY: InitHODSMC
USE MOD_DSMC_ParticlePairing   ,ONLY: DSMC_init_octree
USE MOD_DSMC_SteadyState       ,ONLY: DSMC_SteadyStateInit
USE MOD_DSMC_ChemInit          ,ONLY: DSMC_chemical_init
USE MOD_DSMC_PolyAtomicModel   ,ONLY: InitPolyAtomicMolecs, DSMC_FindFirstVibPick, DSMC_SetInternalEnr_Poly
USE MOD_Particle_Boundary_Vars ,ONLY: nAdaptiveBC, PartBound
USE MOD_Particle_Surfaces_Vars ,ONLY: BCdata_auxSF
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(32)         :: hilf , hilf2
INTEGER               :: iCase, iSpec, jSpec, nCase, iPart, iInit, iPolyatMole, iDOF
REAL                  :: A1, A2     ! species constant for cross section (p. 24 Laux)
REAL                  :: BGGasEVib
INTEGER               :: currentBC, ElemID, iSide, BCSideID, VarNum
#if ( PP_TimeDiscMethod ==42 )
#ifdef CODE_ANALYZE
CHARACTER(LEN=64)     :: DebugElectronicStateFilename
INTEGER               :: ii
#endif
#endif
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' DSMC INIT ...'

! Initialize counter (Count the number of ReactionProb>1)
ReactionProbGTUnityCounter = 0

! reading/writing OutputMesh stuff
DSMC%OutputMeshInit = GETLOGICAL('Particles-DSMC-OutputMeshInit','.FALSE.')
DSMC%OutputMeshSamp = GETLOGICAL('Particles-DSMC-OutputMeshSamp','.FALSE.')
!  IF (DSMC%OutputMeshInit) THEN
!    SWRITE(UNIT_stdOut,'(A)')' WRITING OUTPUT-MESH...'
!    CALL WriteOutputMesh()
!    SWRITE(UNIT_stdOut,'(A)')' WRITING OUTPUT-MESH DONE!'
!  END IF
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
DSMC%MergeSubcells = GETLOGICAL('Particles-DSMC-MergeSubcells','.FALSE.')
IF(DSMC%MergeSubcells.AND.(.NOT.Symmetry2D)) THEN
  CALL abort(__STAMP__&
      ,'ERROR: Merging of subcells only supported within a 2D/axisymmetric simulation!')
END IF
  IF(CollisMode.GE.2) THEN
    DSMC%RotRelaxProb = GETREAL('Particles-DSMC-RotRelaxProb','0.2')
    DSMC%VibRelaxProb = GETREAL('Particles-DSMC-VibRelaxProb','0.004')
  ELSE
    DSMC%RotRelaxProb = 0.
    DSMC%VibRelaxProb = 0.
  END IF
DSMC%ElecRelaxProb = GETREAL('Particles-DSMC-ElecRelaxProb','0.01')
DSMC%GammaQuant   = GETREAL('Particles-DSMC-GammaQuant', '0.5')
ALLOCATE(DSMC%veloMinColl(nSpecies))
DO iSpec = 1, nSpecies
  WRITE(UNIT=hilf,FMT='(I0)') iSpec
  DSMC%veloMinColl(iSpec) = GETREAL('Particles-DSMC-veloMinColl-Spec'//TRIM(hilf),'0.')
END DO
!-----------------------------------------------------------------------------------
! Flag for the automatic calculation of the backward reaction rate with the partition functions and equilibrium constant.
DSMC%BackwardReacRate  = GETLOGICAL('Particles-DSMC-BackwardReacRate','.FALSE.')
! Partition functions are calculated for each species during initialization and stored for values starting with the
! DSMC%PartitionInterval up to DSMC%PartitionMaxTemp, interpolation between the stored values (also used for analytic QK reactions)
DSMC%PartitionMaxTemp  = GETREAL('Particles-DSMC-PartitionMaxTemp','20000')
DSMC%PartitionInterval = GETREAL('Particles-DSMC-PartitionInterval','10')
!-----------------------------------------------------------------------------------
DSMC%CalcQualityFactors = GETLOGICAL('Particles-DSMC-CalcQualityFactors','.FALSE.')
DSMC%ReservoirSimu = GETLOGICAL('Particles-DSMCReservoirSim','.FALSE.')
IF (DSMC%CalcQualityFactors.AND.(CollisMode.LT.1)) THEN
  CALL abort(&
      __STAMP__&
      ,'ERROR: Do not use DSMC%CalcQualityFactors for CollisMode < 1')
END IF ! DSMC%CalcQualityFactors.AND.(CollisMode.LT.1)
DSMC%ReservoirSimuRate       = GETLOGICAL('Particles-DSMCReservoirSimRate','.FALSE.')
DSMC%ReservoirSurfaceRate    = GETLOGICAL('Particles-DSMCReservoirSurfaceRate','.FALSE.')
DSMC%ReservoirRateStatistic  = GETLOGICAL('Particles-DSMCReservoirStatistic','.FALSE.')
DSMC%VibEnergyModel          = GETINT('Particles-ModelForVibrationEnergy','0')
DSMC%DoTEVRRelaxation        = GETLOGICAL('Particles-DSMC-TEVR-Relaxation','.FALSE.')
IF(RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
  IF(DSMC%DoTEVRRelaxation) THEN
    CALL abort(__STAMP__&
        ,'ERROR: Radial weighting or variable time step is not implemented with T-E-V-R relaxation!')
  END IF
END IF
DSMC%ElectronicModel         = GETLOGICAL('Particles-DSMC-ElectronicModel','.FALSE.')
DSMC%ElectronicModelDatabase = TRIM(GETSTR('Particles-DSMCElectronicDatabase','none'))
IF ((DSMC%ElectronicModelDatabase .NE. 'none').AND.&
    ((CollisMode .GT. 1).OR.(CollisMode .EQ. 0))) THEN ! CollisMode=0 is for use of in PIC simulation without collisions
  DSMC%EpsElecBin = GETREAL('EpsMergeElectronicState','1E-4')
ELSEIF(DSMC%ElectronicModel) THEN
  CALL Abort(&
      __STAMP__,&
      'ERROR: Electronic model requires a electronic levels database and CollisMode > 1!')
END IF
IF ((DSMC%VibEnergyModel.EQ.1).AND.(CollisMode.EQ.3)) THEN
  CALL Abort(&
      __STAMP__&
      ,'TSHO Model is not working with Chemical Reactions !!!')
ELSE IF ((DSMC%VibEnergyModel.GT.1).OR.(DSMC%VibEnergyModel.LT.0)) THEN
  CALL Abort(&
      __STAMP__&
      ,'ERROR in ModelForVibrationEnergy Flag!')
END IF
DSMC%NumPolyatomMolecs = 0
! Steady - State Detection: Use Q-Criterion or SSD-Alogrithm?
SamplingActive = .FALSE.
UseQCrit = GETLOGICAL('Particles-DSMC-UseQCrit','.FALSE.')
UseSSD = GETLOGICAL('Particles-DSMC-UseSSD','.FALSE.')
IF(UseQCrit.OR.UseSSD) CALL DSMC_SteadyStateInit()

ALLOCATE(HValue(nElems))
HValue(1:nElems) = 0.0

IF(DSMC%CalcQualityFactors) THEN
  ! 1: Maximal collision probability per cell/subcells (octree)
  ! 2: Mean collision probability within cell
  ! 3: Mean collision separation distance over mean free path
  ! 4: Counter (is not simply the number of iterations in case of a coupled BGK/FP-DSMC simulation)
  VarNum = 4
  ! VarNum + 1: Number of cloned particles per cell
  ! VarNum + 2: Number of identical particles (no relative velocity)
  IF(RadialWeighting%DoRadialWeighting) VarNum = VarNum + 2
  ALLOCATE(DSMC%QualityFacSamp(nElems,VarNum))
  DSMC%QualityFacSamp(1:nElems,1:VarNum) = 0.0
END IF

! definition of DSMC particle values
ALLOCATE(DSMC_RHS(PDM%maxParticleNumber,3))
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
  IF(DoFieldIonization.OR.CollisMode.NE.0)THEN
    CollInf%averagedCollisionParameters     = GETLOGICAL('averagedCollisionParameters','.TRUE.')
    CollInf%crossSectionConstantMode        = GETINT('crossSectionConstantMode','0') 
    ALLOCATE(SpecDSMC(nSpecies))
    DO iSpec = 1, nSpecies
      WRITE(UNIT=hilf,FMT='(I0)') iSpec
      SpecDSMC(iSpec)%Name    = TRIM(GETSTR('Part-Species'//TRIM(hilf)//'-SpeciesName','none'))
      SpecDSMC(iSpec)%InterID = GETINT('Part-Species'//TRIM(hilf)//'-InteractionID','0')
      IF(CollInf%averagedCollisionParameters) THEN
        SpecDSMC(iSpec)%Tref         = GETREAL('Part-Species'//TRIM(hilf)//'-referenceTemperature'      ,'0')
        SpecDSMC(iSpec)%dref         = GETREAL('Part-Species'//TRIM(hilf)//'-referenceDiameter'         ,'0')
        SpecDSMC(iSpec)%omegaLaux    = GETREAL('Part-Species'//TRIM(hilf)//'-temperatureExponentOmega'  ,'0') 
        SpecDSMC(iSpec)%alphaVSS     = GETREAL('Part-Species'//TRIM(hilf)//'-scatteringExponentalphaVSS','1')
        IF((SpecDSMC(iSpec)%InterID*SpecDSMC(iSpec)%Tref*SpecDSMC(iSpec)%dref).EQ.0) THEN !check if species parameters are set
          CALL Abort(&
          __STAMP__&
          ,"ERROR in species data: check collision parameters in ini_2"//&
           "Part-Species-?-(InterID*referenceTemperature*referenceDiameter) is zero")
        END IF 
      END IF ! averagediCollisionParameters 
      SpecDSMC(iSpec)%FullyIonized  = GETLOGICAL('Part-Species'//TRIM(hilf)//'-FullyIonized')
      IF(SpecDSMC(iSpec)%InterID.EQ.4) THEN
        DSMC%ElectronSpecies = iSpec
      END IF ! interID .EQ.4
      ! reading electronic state informations from HDF5 file
      IF((DSMC%ElectronicModelDatabase.NE.'none').AND.(SpecDSMC(iSpec)%InterID.NE.4)) CALL SetElectronicModel(iSpec)
    END DO ! nSpecies

    ! allocate and initialize collision parameter arrays
    ALLOCATE(CollInf%Tref(nSpecies,nSpecies))
    ALLOCATE(CollInf%dref(nSpecies,nSpecies))
    ALLOCATE(CollInf%omegaLaux(nSpecies,nSpecies))
    ALLOCATE(CollInf%alphaVSS(nSpecies,nSpecies))

    DO iSpec=1,nSpecies
      DO jSpec=iSpec,nSpecies
        WRITE(UNIT=hilf,FMT='(I0)')  iSpec
        WRITE(UNIT=hilf2,FMT='(I0)') jSpec
        IF(CollInf%averagedCollisionParameters) THEN ! collision-averaged parameters
          CollInf%Tref(iSpec,jSpec)         = 0.5 * (SpecDSMC(iSpec)%Tref + SpecDSMC(jSpec)%Tref)
          CollInf%dref(iSpec,jSpec)         = 0.5 * (SpecDSMC(iSpec)%dref + SpecDSMC(jSpec)%dref)
          CollInf%omegaLaux(iSpec,jSpec)    = 0.5 * (SpecDSMC(iSpec)%omegaLaux + SpecDSMC(jSpec)%omegaLaux)
          CollInf%alphaVSS(iSpec,jSpec)     = 0.5 * (SpecDSMC(iSpec)%alphaVSS + SpecDSMC(jSpec)%alphaVSS)
        ELSE !  collision-specific parameters
          CollInf%Tref(iSpec,jSpec)     =GETREAL('Part-Collision'//TRIM(hilf)//'-'//TRIM(hilf2)//'-referenceTemperature'      ,'0')
          CollInf%dref(iSpec,jSpec)     =GETREAL('Part-Collision'//TRIM(hilf)//'-'//TRIM(hilf2)//'-referenceDiameter'         ,'0')
          CollInf%omegaLaux(iSpec,jSpec)=GETREAL('Part-Collision'//TRIM(hilf)//'-'//TRIM(hilf2)//'-temperatureExponentOmega'  ,'0')
          CollInf%alphaVSS(iSpec,jSpec) =GETREAL('Part-Collision'//TRIM(hilf)//'-'//TRIM(hilf2)//'-scatteringExponentalphaVSS','1')
        END IF ! averagedCollisionParameters
        IF (iSpec.NE.jSpec) THEN ! fill minor diagonal of collision-info matrices
          CollInf%Tref(jSpec,iSpec)      = CollInf%Tref(iSpec,jSpec)
          CollInf%dref(jSpec,iSpec)      = CollInf%dref(iSpec,jSpec)
          CollInf%omegaLaux(jSpec,iSpec) = CollInf%omegaLaux(iSpec,jSpec)
          CollInf%alphaVSS(jSpec,iSpec)  = CollInf%alphaVSS(iSpec,jSpec)
        END IF ! fill minor diagonal
        IF(CollInf%dref(iSpec,jSpec)*CollInf%Tref(iSpec,jSpec).EQ.0) THEN
            CALL Abort(&
            __STAMP__&
            ,'ERROR: Check collision parameters! (Part-Collision-?-referenceTemperature * referenceDiameter) is zero)')
        END IF ! check if collision parameters are set
      END DO ! jspec=nspec
    END DO ! ispec=nspec
    IF(CollInf%crossSectionConstantMode.EQ.0) THEN
      CollInf%omegaLaux(:,:)=CollInf%omegaLaux(1,1)
    END IF ! CollInf%crossSectionConstantMode=0
  END IF ! DoFieldIonization.OR.CollisMode.NE.0

! allocate internal energy arrays
IF ( DSMC%ElectronicModel ) THEN
  ALLOCATE(PartStateIntEn(PDM%maxParticleNumber,3))
ELSE
  ALLOCATE(PartStateIntEn(PDM%maxParticleNumber,2))
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
  DSMC%NumColl = 0
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
  ALLOCATE(CollInf%crossSectionConstantCab(nCase))
  ALLOCATE(CollInf%KronDelta(nCase))
  CollInf%crossSectionConstantCab = 0
  CollInf%KronDelta = 0

  DO iSpec = 1, nSpecies
    DO jSpec = iSpec, nSpecies
      iCase = CollInf%Coll_Case(iSpec,jSpec)
      IF (iSpec.eq.jSpec) THEN
        CollInf%KronDelta(iCase) = 1
      ELSE !ispec.NE.jspec
        CollInf%KronDelta(iCase) = 0
      END IF !kronecker delta
      ! Laux (2.37) prefactor crossSectionConstantCab calculation depending on omegaLaux coll-averaged or -specific
      SELECT CASE (CollInf%crossSectionConstantMode) ! sigma=Cab * cr^(-2 omegaLaux)
      CASE (0,1) 
        ! A1,A2 species constants, see laux1996 (2.38),(2.39) 
        A1 = 0.5 * SQRT(Pi) * CollInf%dref(iSpec,iSpec)*(2*BoltzmannConst* CollInf%Tref(iSpec,iSpec))** &
           (CollInf%omegaLaux(iSpec,iSpec)*0.5) /SQRT(GAMMA(2.0 - CollInf%omegaLaux(iSpec,iSpec)))
        A2 = 0.5 * SQRT(Pi) * CollInf%dref(jSpec,jSpec)*(2*BoltzmannConst*CollInf%Tref(jSpec,jSpec))** &
           (CollInf%omegaLaux(jSpec,jSpec)*0.5) /SQRT(GAMMA(2.0 - CollInf%omegaLaux(jSpec,jSpec)))
        CollInf%crossSectionConstantCab(iCase) = (A1 + A2)**2 * CollInf%MassRed(iCase)** ( - CollInf%omegaLaux(iSpec,jSpec))
      CASE (2) ! cross section constant Cab without Laux simplification (needed if multiple omegas are used) see bird1981 (9)
        CollInf%crossSectionConstantCab(iCase) = (SQRT(Pi) * CollInf%dref(iSpec,jSpec)*                                           &
                                               (2*BoltzmannConst*CollInf%Tref(iSpec,jSpec))**(CollInf%omegaLaux(iSpec,jSpec)*0.5) &
                                               /SQRT(GAMMA(2.0 - CollInf%omegaLaux(iSpec,jSpec))))**2                             &
                                               * CollInf%MassRed(iCase)** ( - CollInf%omegaLaux(iSpec,jSpec))
      END SELECT !crossSectionConstant
    END DO !jspec=nspecies
  END DO ! ispec=nspecies
!-----------------------------------------------------------------------------------------------------------------------------------
! reading BG Gas stuff (required for the temperature definition in iInit=0)
!-----------------------------------------------------------------------------------------------------------------------------------
  !...moved to InitializeVariables!!!

!-----------------------------------------------------------------------------------------------------------------------------------
! reading/writing molecular stuff
!-----------------------------------------------------------------------------------------------------------------------------------
  ! Check whether calculation of instantaneous translational temperature is required
  IF(((CollisMode.GT.1).AND.(SelectionProc.EQ.2)).OR.((CollisMode.EQ.3).AND.DSMC%BackwardReacRate).OR.DSMC%CalcQualityFactors &
            .OR.(DSMC%VibRelaxProb.EQ.2)) THEN
    ! 1. Case: Inelastic collisions and chemical reactions with the Gimelshein relaxation procedure and variable vibrational
    !           relaxation probability (CalcGammaVib)
    ! 2. Case: Chemical reactions and backward rate require cell temperature for the partition function and equilibrium constant
    ! 3. Case: Temperature required for the mean free path with the VHS model
    ALLOCATE(DSMC%InstantTransTemp(nSpecies+1))
    DSMC%InstantTransTemp = 0.0
  END IF

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
          IF (DSMC%VibEnergyModel.EQ.0) THEN
            SpecDSMC(iSpec)%MaxVibQuant = 200
          ELSE
            SpecDSMC(iSpec)%MaxVibQuant = INT(SpecDSMC(iSpec)%Ediss_eV*ElementaryCharge/&
                (BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)) + 1
          END IF
          ! Calculation of the zero-point energy
          SpecDSMC(iSpec)%EZeroPoint = DSMC%GammaQuant * BoltzmannConst * SpecDSMC(iSpec)%CharaTVib
          ! Calculation of the dissociation quantum number (used for QK chemistry)
          SpecDSMC(iSpec)%DissQuant = INT(SpecDSMC(iSpec)%Ediss_eV*ElementaryCharge/(BoltzmannConst*SpecDSMC(iSpec)%CharaTVib))
        END IF
        SpecDSMC(iSpec)%VFD_Phi3_Factor = GETREAL('Part-Species'//TRIM(hilf)//'-VFDPhi3','0.')
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
        DO iInit = 0, Species(iSpec)%NumberOfInits
          IF (iInit .EQ. 0) THEN ! 0. entry := old style parameter def. (default values if not def., some values might be needed)
            hilf2=TRIM(hilf)
          ELSE ! iInit >0
            WRITE(UNIT=hilf2,FMT='(I0)') iInit
            hilf2=TRIM(hilf)//'-Init'//TRIM(hilf2)
          END IF ! iInit
          IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
            IF (Species(iSpec)%Init(iInit)%ElemTVibFileID.EQ.0) THEN
              SpecDSMC(iSpec)%Init(iInit)%TVib      = GETREAL('Part-Species'//TRIM(hilf2)//'-TempVib','0.')
              IF (SpecDSMC(iSpec)%Init(iInit)%TVib.EQ.0.) THEN
                IF (iInit.EQ.0)THEN
                  IF (Species(iSpec)%StartnumberOfInits.EQ.0)THEN
                    CALL Abort(&
                        __STAMP__&
                        ,'Error! TVib needs to be defined in Part-SpeciesXX-TempVib for iSpec',iSpec)
                  ELSE IF (BGGas%BGGasSpecies.EQ.iSpec) THEN !cases which need values of fixed iInit=0 (indep. from Startnr.OfInits)
                    CALL Abort(&
                        __STAMP__&
                        ,'Error! TVib needs to be defined in Part-SpeciesXX-TempVib for BGGas')
                  END IF
                ELSE ! iInit >0
                  CALL Abort(&
                      __STAMP__&
                      ,'Error! TVib needs to be defined in Part-SpeciesXX-InitXX-TempVib for iSpec, iInit'&
                  ,iSpec,REAL(iInit))
                END IF
              END IF
            END IF !ElemMacroRestart TVib
            IF (Species(iSpec)%Init(iInit)%ElemTRotFileID.EQ.0) THEN
              SpecDSMC(iSpec)%Init(iInit)%TRot      = GETREAL('Part-Species'//TRIM(hilf2)//'-TempRot','0.')
              IF (SpecDSMC(iSpec)%Init(iInit)%TRot.EQ.0.) THEN
                IF (iInit.EQ.0)THEN
                  IF (Species(iSpec)%StartnumberOfInits.EQ.0)THEN
                    CALL Abort(&
                        __STAMP__&
                        ,'Error! TRot needs to be defined in Part-SpeciesXX-TempRot for iSpec',iSpec)
                  ELSE IF (BGGas%BGGasSpecies.EQ.iSpec) THEN !cases which need values of fixed iInit=0 (indep. from Startnr.OfInits)
                    CALL Abort(&
                        __STAMP__&
                        ,'Error! TRot needs to be defined in Part-SpeciesXX-TempRot for BGGas')
                  END IF
                ELSE ! iInit >0
                  CALL Abort(&
                      __STAMP__&
                      ,'Error! TRot needs to be defined in Part-SpeciesXX-InitXX-TempRot for iSpec, iInit'&
                  ,iSpec,REAL(iInit))
                END IF
              END IF
            END IF
          END IF ! ElemMacroRestart TRot
          ! read electronic temperature
          IF ( DSMC%ElectronicModel ) THEN
            IF (Species(iSpec)%Init(iInit)%ElemTElecFileID.EQ.0) THEN
              SpecDSMC(iSpec)%Init(iInit)%Telec   = GETREAL('Part-Species'//TRIM(hilf2)//'-TempElec','0.')
              IF (SpecDSMC(iSpec)%Init(iInit)%Telec.EQ.0.) THEN
                IF (iInit.EQ.0)THEN
                  IF (Species(iSpec)%StartnumberOfInits.EQ.0)THEN
                    CALL Abort(&
                        __STAMP__&
                        ,' Error! Telec needs to defined in Part-SpeciesXX-Tempelec for Species',iSpec)
                  ELSE IF (BGGas%BGGasSpecies.EQ.iSpec) THEN !cases which need values of fixed iInit=0 (indep. from Startnr.OfInits)
                    CALL Abort(&
                        __STAMP__&
                        ,' Error! Telec needs to defined in Part-SpeciesXX-Tempelec for BGGas')
                  END IF
                ELSE ! iInit >0
                  CALL Abort(&
                      __STAMP__&
                      ,' Error! Telec needs to defined in Part-SpeciesXX-InitXX-Tempelc for iSpec, iInit',iSpec,REAL(iInit))
                END IF
              END IF
            END IF !ElemMacroRestart TElec
          END IF ! electronic model
        END DO !Inits
        ALLOCATE(SpecDSMC(iSpec)%Surfaceflux(1:Species(iSpec)%nSurfacefluxBCs+nAdaptiveBC))
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
          IF ( DSMC%ElectronicModel ) THEN
            SpecDSMC(iSpec)%Surfaceflux(iInit)%Telec   = GETREAL('Part-Species'//TRIM(hilf2)//'-TempElec','0.')
            IF (SpecDSMC(iSpec)%Surfaceflux(iInit)%Telec.EQ.0.) THEN
              CALL Abort(&
                  __STAMP__&
                  ,' Error! Telec not defined in Part-SpeciesXX-SurfacefluxXX-Tempelec for iSpec, iInit',iSpec,REAL(iInit))
            END IF
          END IF
        END DO !SurfaceFluxBCs
        ! add Adaptive boundaries
        ! initialize  rot, vib and elec temperature of macrovalues
        DO iInit = (Species(iSpec)%nSurfacefluxBCs+1),(Species(iSpec)%nSurfacefluxBCs+nAdaptiveBC)
          ! read rot and vib temperatures
          IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
            SpecDSMC(iSpec)%Surfaceflux(iInit)%TVib      = Partbound%AdaptiveTemp(Species(iSpec)%Surfaceflux(iInit)%BC)
            !SpecDSMC(iSpec)%Init(0)%TVib
            SpecDSMC(iSpec)%Surfaceflux(iInit)%TRot      = Partbound%AdaptiveTemp(Species(iSpec)%Surfaceflux(iInit)%BC)
            !SpecDSMC(iSpec)%Init(0)%TRot
            IF (SpecDSMC(iSpec)%Surfaceflux(iInit)%TRot*SpecDSMC(iSpec)%Surfaceflux(iInit)%TVib.EQ.0.) THEN
              CALL Abort(&
                  __STAMP__&
                  ,'Error! TVib and TRot not def. in Part-SpeciesXX-SurfacefluxXX-TempVib/TempRot for iSpec, iInit',iSpec,REAL(iInit))
            END IF
            currentBC = Species(iSpec)%Surfaceflux(iInit)%BC !go through sides if present in proc...
            IF (BCdata_auxSF(currentBC)%SideNumber.GT.0) THEN
              DO iSide=1,BCdata_auxSF(currentBC)%SideNumber
                BCSideID=BCdata_auxSF(currentBC)%SideList(iSide)
                ElemID = SideToElem(1,BCSideID)
                IF (ElemID.LT.1) THEN !not sure if necessary
                  ElemID = SideToElem(2,BCSideID)
                END IF
                Adaptive_MacroVal(8,ElemID,iSpec) = SpecDSMC(iSpec)%Surfaceflux(iInit)%TVib
                Adaptive_MacroVal(9,ElemID,iSpec) = SpecDSMC(iSpec)%Surfaceflux(iInit)%TRot
              END DO
            ELSE IF (BCdata_auxSF(currentBC)%SideNumber.EQ.-1) THEN
              CALL abort(&
                  __STAMP__&
                  ,'ERROR in DSMC_init of rot, vib and elec_shell: Someting is wrong with SideNumber of BC ',currentBC)
            END IF
          END IF
          ! read electronic temperature
          IF ( DSMC%ElectronicModel ) THEN
            SpecDSMC(iSpec)%Surfaceflux(iInit)%Telec   =  Partbound%AdaptiveTemp(Species(iSpec)%Surfaceflux(iInit)%BC)
            !SpecDSMC(iSpec)%Init(0)%Telec
            IF (SpecDSMC(iSpec)%Surfaceflux(iInit)%Telec.EQ.0.) THEN
              CALL Abort(&
                  __STAMP__&
                  ,' Error! Telec not defined in Part-SpeciesXX-SurfacefluxXX-Tempelec for iSpec, iInit',iSpec,REAL(iInit))
            END IF
            currentBC = Species(iSpec)%Surfaceflux(iInit)%BC !go through sides if present in proc...
            IF (BCdata_auxSF(currentBC)%SideNumber.GT.0) THEN
              DO iSide=1,BCdata_auxSF(currentBC)%SideNumber
                BCSideID=BCdata_auxSF(currentBC)%SideList(iSide)
                ElemID = SideToElem(1,BCSideID)
                IF (ElemID.LT.1) THEN !not sure if necessary
                  ElemID = SideToElem(2,BCSideID)
                END IF
                Adaptive_MacroVal(10,ElemID,iSpec) = SpecDSMC(iSpec)%Surfaceflux(iInit)%Telec
              END DO
            ELSE IF (BCdata_auxSF(currentBC)%SideNumber.EQ.-1) THEN
              CALL abort(&
                  __STAMP__&
                  ,'ERROR in DSMC_init of elec_temperatur: Someting is wrong with SideNumber of BC ',currentBC)
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
          !          ! Required if the Metropolis-Hastings random-walk is utilized for the initialization of polyatomic molecules (sampling
          !          ! of all modes at once)
          !          iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
          !          ALLOCATE( &
          !              PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(1:PolyatomMolDSMC(iPolyatMole)%VibDOF, &
          !                                                             0:Species(iSpec)%NumberOfInits+Species(iSpec)%nSurfacefluxBCs))
          !          DO iInit = Species(iSpec)%StartnumberOfInits, Species(iSpec)%NumberOfInits
          !            CALL DSMC_FindFirstVibPick(iInit, iSpec, 1)
          !          END DO
          !          DO iInit = 1,Species(iSpec)%nSurfacefluxBCs
          !            CALL DSMC_FindFirstVibPick(iInit, iSpec, 2)
          !          END DO
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
#if (PP_TimeDiscMethod==42)
#ifdef CODE_ANALYZE
    IF ( DSMC%ElectronicModel ) THEN
      DO iSpec = 1, nSpecies
        IF ( (SpecDSMC(iSpec)%InterID .eq. 4).OR.SpecDSMC(iSpec)%FullyIonized) THEN
          SpecDSMC(iSpec)%MaxElecQuant = 0
        ELSE
          ALLOCATE( SpecDSMC(iSpec)%levelcounter         ( 0:size(SpecDSMC(iSpec)%ElectronicState,2)-1) , &
                    SpecDSMC(iSpec)%dtlevelcounter       ( 0:size(SpecDSMC(ispec)%ElectronicState,2)-1) , &
                    SpecDSMC(iSpec)%ElectronicTransition ( 1:nSpecies                                   , &
                                                           0:size(SpecDSMC(ispec)%ElectronicState,2)-1,   &
                                                           0:size(SpecDSMC(ispec)%ElectronicState,2)-1)   )
          SpecDSMC(iSpec)%levelcounter         = 0
          SpecDSMC(iSpec)%dtlevelcounter       = 0
          SpecDSMC(iSpec)%ElectronicTransition = 0
        END IF
      END DO
    END IF
#endif
#endif
    ! Setting the internal energy value of every particle
    DO iPart = 1, PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        IF (Species(PartSpecies(iPart))%NumberOfInits.EQ.0) THEN
          IF (SpecDSMC(PartSpecies(iPart))%PolyatomicMol) THEN
            CALL DSMC_SetInternalEnr_Poly(PartSpecies(iPart),0,iPart,1)
          ELSE
            CALL DSMC_SetInternalEnr_LauxVFD(PartSpecies(iPart),0,iPart,1)
          END IF
        ELSE
          iInit = PDM%PartInit(iPart)
          IF (SpecDSMC(PartSpecies(iPart))%PolyatomicMol) THEN
            CALL DSMC_SetInternalEnr_Poly(PartSpecies(iPart),iInit,iPart,1)
          ELSE
            CALL DSMC_SetInternalEnr_LauxVFD(PartSpecies(iPart),iInit,iPart,1)
          END IF
        END IF
      END IF
    END DO

#if (PP_TimeDiscMethod==42)
#ifdef CODE_ANALYZE
    ! Debug Output for initialized electronic state
    IF ( DSMC%ElectronicModel ) THEN
      DO iSpec = 1, nSpecies
        print*,SpecDSMC(iSpec)%InterID
        IF ((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
          IF (  SpecDSMC(iSpec)%levelcounter(0) .ne. 0) THEN
            WRITE(DebugElectronicStateFilename,'(I2.2)') iSpec
            DebugElectronicStateFilename = 'Initial_Electronic_State_Species_'//trim(DebugElectronicStateFilename)//'.dat'
            open(unit=483,file=DebugElectronicStateFilename,form='formatted',status='unknown')
            DO ii = 0, SpecDSMC(iSpec)%MaxElecQuant - 1
              WRITE(483,'(I3.1,3x,F12.7)') ii, REAL( SpecDSMC(iSpec)%levelcounter(ii) ) / &
                                               REAL( Species(iSpec)%Init(0)%initialParticleNumber )
            END DO
            close(unit=483)
          END IF
        END IF
      END DO
    END IF
#endif
#endif

#if (PP_TimeDiscMethod!=300)
    DEALLOCATE(PDM%PartInit)
#endif
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
      SpecDSMC(iSpec)%PreviousState = GETINT('Part-Species'//TRIM(hilf)//'-PreviousState','0')
      IF((SpecDSMC(iSpec)%InterID.EQ.10).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
        IF(SpecDSMC(iSpec)%PreviousState.EQ.0) THEN
          CALL abort(&
              __STAMP__&
              ,'ERROR: Please specify the previous state of the ion species:', iSpec)
        END IF
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
      ! Read-in of species parameters for the partition function calculation -------------------------------------------------------
      IF(DSMC%BackwardReacRate) THEN
        IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
          SpecDSMC(iSpec)%SymmetryFactor              = GETINT('Part-Species'//TRIM(hilf)//'-SymmetryFactor','0')
          IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
            iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
            IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec) THEN
              IF(PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)*SpecDSMC(iSpec)%SymmetryFactor.EQ.0) THEN
                CALL abort(&
                    __STAMP__&
                    ,'ERROR: Char. rotational temperature or symmetry factor not defined properly for backward rate!', iSpec)
              END IF
            ELSE
              IF(PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)*PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(2)  &
                  * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(3)*SpecDSMC(iSpec)%SymmetryFactor.EQ.0) THEN
                CALL abort(&
                    __STAMP__&
                    ,'ERROR: Char. rotational temperature or symmetry factor not defined properly for backward rate!', iSpec)
              END IF
            END IF
          ELSE
            IF(SpecDSMC(iSpec)%CharaTRot*SpecDSMC(iSpec)%SymmetryFactor.EQ.0) THEN
              CALL abort(&
                  __STAMP__&
                  ,'ERROR: Char. rotational temperature or symmetry factor not defined properly for backward rate!', iSpec)
            END IF
          END IF
        END IF
        IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
          IF(.NOT.ALLOCATED(SpecDSMC(iSpec)%ElectronicState)) THEN
            CALL abort(&
                __STAMP__&
                ,'ERROR: Electronic energy levels required for the calculation of backward reaction rate!',iSpec)
          END IF
        END IF
      END IF
      !-----------------------------------------------------------------------------------------------------------------------------
      SpecDSMC(iSpec)%Eion_eV               = GETREAL('Part-Species'//TRIM(hilf)//'-IonizationEn_eV','0')
      SpecDSMC(iSpec)%RelPolarizability     = GETREAL('Part-Species'//TRIM(hilf)//'-RelPolarizability','0')
      SpecDSMC(iSpec)%NumEquivElecOutShell  = GETINT('Part-Species'//TRIM(hilf)//'-NumEquivElecOutShell','0')
      SpecDSMC(iSpec)%NumOfPro              = GETINT('Part-Species'//TRIM(hilf)//'-NumOfProtons','0')
      IF((SpecDSMC(iSpec)%Eion_eV*SpecDSMC(iSpec)%RelPolarizability*SpecDSMC(iSpec)%NumEquivElecOutShell &
          *SpecDSMC(iSpec)%NumOfPro).EQ.0) THEN
        SWRITE(*,*) "Ionization parameters are not defined for species:", iSpec
      END IF
    END DO

    ! Calculating the heat of formation for ionized species (including higher ionization levels)
    ! Requires the completed read-in of species data
    CALL CalcHeatOfFormation()

    ! Set "NextIonizationSpecies" information for field ionization from "PreviousState" info
    ! NextIonizationSpecies => SpeciesID of the next higher ionization level
    CALL SetNextIonizationSpecies()

    CALL DSMC_chemical_init()
  ELSE IF (ANY(PartBound%Reactive) .AND. CollisMode.GT.1) THEN
    DO iSpec = 1, nSpecies
      WRITE(UNIT=hilf,FMT='(I0)') iSpec
      IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
        SpecDSMC(iSpec)%SymmetryFactor              = GETINT('Part-Species'//TRIM(hilf)//'-SymmetryFactor','0')
        IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
          iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
          IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec) THEN
            IF(PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)*SpecDSMC(iSpec)%SymmetryFactor.EQ.0) THEN
              CALL abort(&
                  __STAMP__&
                  ,'ERROR: Char. rotational temperature or symmetry factor not defined properly for Adsorptionmodel!', iSpec)
            END IF
          ELSE
            IF(PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)*PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(2)  &
                * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(3)*SpecDSMC(iSpec)%SymmetryFactor.EQ.0) THEN
              CALL abort(&
                  __STAMP__&
                  ,'ERROR: Char. rotational temperature or symmetry factor not defined properly for Adsorptionmodel!', iSpec)
            END IF
          END IF
        ELSE
          IF(SpecDSMC(iSpec)%CharaTRot*SpecDSMC(iSpec)%SymmetryFactor.EQ.0) THEN
            CALL abort(&
                __STAMP__&
                ,'ERROR: Char. rotational temperature or symmetry factor not defined properly for Adsorptionmodel!', iSpec)
          END IF
        END IF
      END IF
    END DO
  END IF

  !-----------------------------------------------------------------------------------------------------------------------------------
  ! Octree + Nearest Neighbour: octree cell splitting is performed until the mean free path is resolved or a certain number of
  ! particles is reached. Afterwards a nearest neighbour search for the collision partner is performed.
  ! Source: Pfeiffer, M., Mirza, A. and Fasoulas, S. (2013). A grid-independent particle pairing strategy for DSMC.
  ! Journal of Computational Physics 246, 28–36. doi:10.1016/j.jcp.2013.03.018
  !-----------------------------------------------------------------------------------------------------------------------------------
  DSMC%UseOctree = GETLOGICAL('Particles-DSMC-UseOctree')
  IF(DSMC%UseOctree) THEN
    DSMC%UseNearestNeighbour = GETLOGICAL('Particles-DSMC-UseNearestNeighbour')
    IF((.NOT.Symmetry2D).AND.(.NOT.DSMC%UseNearestNeighbour)) THEN
      CALL abort(&
          __STAMP__&
          ,'Statistical Pairing with Octree not yet supported in 3D!')
    END IF
  END IF
  ! If number of particles is greater than OctreePartNumNode, cell is going to be divided for performance of nearest neighbour
  IF(Symmetry2D) THEN
    DSMC%PartNumOctreeNode = GETINT('Particles-OctreePartNumNode','40')
  ELSE
    DSMC%PartNumOctreeNode = GETINT('Particles-OctreePartNumNode','80')
  END IF
  ! If number of particles is less than OctreePartNumNodeMin, cell is NOT going to be split even if mean free path is not resolved
  ! 3D: 50/8; 2D: 28/4 -> ca. 6-7 particles per cell
  IF(Symmetry2D) THEN
    DSMC%PartNumOctreeNodeMin = GETINT('Particles-OctreePartNumNodeMin','28')
  ELSE
    DSMC%PartNumOctreeNodeMin = GETINT('Particles-OctreePartNumNodeMin','50')
  END IF
  IF (DSMC%PartNumOctreeNodeMin.LT.20) THEN
    CALL abort(&
        __STAMP__&
        ,'ERROR: Given Particles-OctreePartNumNodeMin is less than 20!')
  END IF
  IF(DSMC%UseOctree) THEN
    IF(NGeo.GT.PP_N) CALL abort(&
        __STAMP__&
        ,' Set PP_N to NGeo, else, the volume is not computed correctly.')
    CALL DSMC_init_octree()
  END IF
  IF(Symmetry2D) THEN
    CollInf%ProhibitDoubleColl = GETLOGICAL('Particles-DSMC-ProhibitDoubleCollisions','.TRUE.')
    IF (CollInf%ProhibitDoubleColl) THEN
      IF(.NOT.ALLOCATED(CollInf%OldCollPartner)) ALLOCATE(CollInf%OldCollPartner(1:PDM%maxParticleNumber))
      CollInf%OldCollPartner = 0
    END IF
  ELSE
    IF (CollInf%ProhibitDoubleColl) THEN
      CollInf%ProhibitDoubleColl = GETLOGICAL('Particles-DSMC-ProhibitDoubleCollisions','.FALSE.')
      CALL abort(__STAMP__,&
          'ERROR: Prohibiting double collisions is only supported within a 2D/axisymmetric simulation!')
    END IF
  END IF

  !-----------------------------------------------------------------------------------------------------------------------------------
  ! Set mean VibQua of BGGas for dissoc reaction
  !-----------------------------------------------------------------------------------------------------------------------------------
  IF (BGGas%BGGasSpecies.NE.0) THEN
    IF (DSMC%UseOctree) THEN
      CALL abort(__STAMP__,&
          'ERROR: Utilization of the octree and nearest neighbour scheme not possible with the background gas')
    END IF
    IF(SpecDSMC(BGGas%BGGasSpecies)%InterID.EQ.4) THEN
      CALL abort(__STAMP__,&
          'ERROR: Electrons as background gas are not yet available!!')
    END IF
    IF((SpecDSMC(BGGas%BGGasSpecies)%InterID.EQ.2).OR.(SpecDSMC(BGGas%BGGasSpecies)%InterID.EQ.20)) THEN
      IF(SpecDSMC(BGGas%BGGasSpecies)%PolyatomicMol) THEN
        CALL abort(&
            __STAMP__&
            ,'ERROR: Polyatomic species as background gas are not yet available!')
      ELSE
        BGGasEVib = DSMC%GammaQuant * BoltzmannConst * SpecDSMC(BGGas%BGGasSpecies)%CharaTVib &
            + BoltzmannConst * SpecDSMC(BGGas%BGGasSpecies)%CharaTVib  &
            /  (EXP(SpecDSMC(BGGas%BGGasSpecies)%CharaTVib / SpecDSMC(BGGas%BGGasSpecies)%Init(0)%TVib) - 1) &
            - BoltzmannConst * SpecDSMC(BGGas%BGGasSpecies)%CharaTVib * SpecDSMC(BGGas%BGGasSpecies)%MaxVibQuant &
            / (EXP(SpecDSMC(BGGas%BGGasSpecies)%CharaTVib * SpecDSMC(BGGas%BGGasSpecies)%MaxVibQuant &
            / SpecDSMC(BGGas%BGGasSpecies)%Init(0)%TVib) - 1)
        BGGasEVib = BGGasEVib/(BoltzmannConst*SpecDSMC(BGGas%BGGasSpecies)%CharaTVib) - DSMC%GammaQuant
        BGGas%BGMeanEVibQua = MIN(INT(BGGasEVib) + 1, SpecDSMC(BGGas%BGGasSpecies)%MaxVibQuant)
      END IF
    ELSE
      BGGas%BGMeanEVibQua = 0
    END IF
  END IF

!-----------------------------------------------------------------------------------------------------------------------------------
! Calculate vib collision numbers and characteristic velocity, according to Abe
!-----------------------------------------------------------------------------------------------------------------------------------
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
                                  * ( 2.* (2.-SpecDSMC(iSpec)%omegaVHS) * BoltzmannConst * SpecDSMC(iSpec)%TrefVHS &
                                  / CollInf%MassRed(iCase) ) ** SpecDSMC(iSpec)%omegaVHS
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
      ALLOCATE(DSMC%CalcVibProb(1:nSpecies,1:3))
      DSMC%QualityFacSampVib = 0.
      DSMC%QualityFacSampVibSamp = 0
      DSMC%CalcVibProb = 0.
    END IF
    CALL SetVarVibProb2Elems()
    ! CHeck if DSMC%InstantTransTemp is still needed
    IF(.NOT.(((CollisMode.GT.1).AND.(SelectionProc.EQ.2)).OR.((CollisMode.EQ.3).AND.DSMC%BackwardReacRate) &
            .OR.DSMC%CalcQualityFactors)) THEN
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

END IF ! CollisMode.GT.0

! If field ionization is used without chemical reactions due to collisions (DSMC chemistry)
IF(DoFieldIonization.AND.(CollisMode.NE.3))THEN
  DO iSpec = 1, nSpecies
    WRITE(UNIT=hilf,FMT='(I0)') iSpec
    ! Heat of formation of ionized species is modified with the ionization energy directly from read-in electronic energy levels
    ! of the ground/previous state of the respective species (Input requires a species number (eg species number of N for NIon1))
    SpecDSMC(iSpec)%PreviousState = GETINT('Part-Species'//TRIM(hilf)//'-PreviousState','0')
    IF((SpecDSMC(iSpec)%InterID.EQ.10).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
      IF(SpecDSMC(iSpec)%PreviousState.EQ.0) THEN
        CALL abort(&
            __STAMP__&
            ,'ERROR: Please specify the previous state of the ion species:', iSpec)
      END IF ! SpecDSMC(iSpec)%PreviousState.EQ.0
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
USE MOD_Globals       ,ONLY: abort,UNIT_stdOut
#if USE_MPI
USE MOD_Globals       ,ONLY: mpiroot
#endif
USE MOD_Globals_Vars  ,ONLY: BoltzmannConst
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
INTEGER       :: iSpec,jSpec,counter,MaxElecQua
!===================================================================================================================================
AutoDetect=.TRUE.
DO iSpec = 1, nSpecies
  counter = 0
  IF((SpecDSMC(iSpec)%InterID.EQ.10).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
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


SUBROUTINE DSMC_SetInternalEnr_LauxVFD(iSpecies, iInit, iPart, init_or_sf)
!===================================================================================================================================
! Energy distribution according to dissertation of Laux (diatomic)
!===================================================================================================================================
! MODULES
  USE MOD_Globals,               ONLY : abort
  USE MOD_Globals_Vars,          ONLY : BoltzmannConst
  USE MOD_DSMC_Vars,             ONLY : PartStateIntEn, SpecDSMC, DSMC
  USE MOD_Particle_Vars,         ONLY : Species, PEM, Adaptive_MacroVal
  USE MOD_Particle_Boundary_Vars,ONLY: PartBound
  USE MOD_DSMC_ElectronicModel,  ONLY : InitElectronShell
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iSpecies, iInit, iPart, init_or_sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                          :: iRan
  INTEGER                       :: iQuant
  REAL                        :: TVib                       ! vibrational temperature
  REAL                        :: TRot                       ! rotational temperature
  INTEGER                     :: ElemID
  REAL                        :: pressure
!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
! Set internal energies (vibrational and rotational)
!-----------------------------------------------------------------------------------------------------------------------------------
  IF ((SpecDSMC(iSpecies)%InterID.EQ.2).OR.(SpecDSMC(iSpecies)%InterID.EQ.20)) THEN
    SELECT CASE (init_or_sf)
    CASE(1) !iInit
      IF (Species(iSpecies)%Init(iInit)%ElemTVibFileID.EQ.0) THEN
        TVib=SpecDSMC(iSpecies)%Init(iInit)%TVib
      ELSE
        TVib=Species(iSpecies)%Init(iInit)%ElemTVib(PEM%Element(iPart))
      END IF
      IF (Species(iSpecies)%Init(iInit)%ElemTRotFileID.EQ.0) THEN
        TRot=SpecDSMC(iSpecies)%Init(iInit)%TRot
      ELSE
        TRot=Species(iSpecies)%Init(iInit)%ElemTRot(PEM%Element(iPart))
      END IF
    CASE(2) !SurfaceFlux
      IF(iInit.GT.Species(iSpecies)%nSurfacefluxBCs)THEN
        !-- compute number of to be inserted particles
        SELECT CASE(PartBound%AdaptiveType(Species(iSpecies)%Surfaceflux(iInit)%BC))
        CASE(1) ! Pressure inlet (pressure, temperature const)
          TVib=SpecDSMC(iSpecies)%SurfaceFlux(iInit)%TVib
          TRot=SpecDSMC(iSpecies)%SurfaceFlux(iInit)%TRot
        CASE(2) ! adaptive Outlet/freestream
          ElemID = PEM%Element(iPart)
          pressure = PartBound%AdaptivePressure(Species(iSpecies)%Surfaceflux(iInit)%BC)
          TVib = pressure / (BoltzmannConst * SUM(Adaptive_MacroVal(7,ElemID,:)))
          TRot = TVib
          !TVib=Adaptive_MacroVal(8,ElemID,iSpecies)
          !TRot=Adaptive_MacroVal(9,ElemID,iSpecies)
        CASE(3) ! pressure outlet (pressure defined)
        CASE DEFAULT
          CALL abort(&
__STAMP__&
,'wrong adaptive type for Surfaceflux in int_energy -> lauxVDF!')
        END SELECT
      ELSE
        IF(Species(iSpecies)%Surfaceflux(iInit)%Adaptive) THEN
          SELECT CASE(Species(iSpecies)%Surfaceflux(iInit)%AdaptiveType)
            CASE(1,3,4) ! Pressure and massflow inlet (pressure/massflow, temperature const)
              TVib=SpecDSMC(iSpecies)%Surfaceflux(iInit)%TVib
              TRot=SpecDSMC(iSpecies)%Surfaceflux(iInit)%TRot
            CASE(2) ! adaptive Outlet/freestream
              ElemID = PEM%Element(iPart)
              TVib = Species(iSpecies)%Surfaceflux(iInit)%AdaptivePressure &
                      / (BoltzmannConst * Adaptive_MacroVal(DSMC_NUMDENS,ElemID,iSpecies))
              TRot = TVib
            CASE DEFAULT
              CALL abort(&
              __STAMP__&
              ,'Wrong adaptive type for Surfaceflux in int_energy -> lauxVDF!')
          END SELECT
        ELSE
          TVib=SpecDSMC(iSpecies)%Surfaceflux(iInit)%TVib
          TRot=SpecDSMC(iSpecies)%Surfaceflux(iInit)%TRot
        END IF
      END IF
    CASE DEFAULT
      CALL abort(&
      __STAMP__&
      ,'neither iInit nor Surfaceflux defined as reference!')
    END SELECT
    ! Set vibrational energy
    CALL RANDOM_NUMBER(iRan)
    iQuant = INT(-LOG(iRan)*TVib/SpecDSMC(iSpecies)%CharaTVib)
    DO WHILE (iQuant.GE.SpecDSMC(iSpecies)%MaxVibQuant)
      CALL RANDOM_NUMBER(iRan)
      iQuant = INT(-LOG(iRan)*TVib/SpecDSMC(iSpecies)%CharaTVib)
    END DO
    !evtl muss partstateinten nochmal geändert werden, mpi, resize etc..
    PartStateIntEn(iPart, 1) = (iQuant + DSMC%GammaQuant)*SpecDSMC(iSpecies)%CharaTVib*BoltzmannConst
    ! Set rotational energy
    CALL RANDOM_NUMBER(iRan)
    PartStateIntEn(iPart, 2) = -BoltzmannConst*TRot*LOG(iRan)
  ELSE
    ! Nullify energy for atomic species
    PartStateIntEn(iPart, 1) = 0
    PartStateIntEn(iPart, 2) = 0
  END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! Set electronic energy
!-----------------------------------------------------------------------------------------------------------------------------------
  IF (DSMC%ElectronicModel) THEN
    IF((SpecDSMC(iSpecies)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpecies)%FullyIonized)) THEN
      CALL InitElectronShell(iSpecies,iPart,iInit,init_or_sf)
    ELSE
      PartStateIntEn(iPart, 3) = 0.
    END IF
  ENDIF

END SUBROUTINE DSMC_SetInternalEnr_LauxVFD


SUBROUTINE SetVarVibProb2Elems()
!===================================================================================================================================
! Set initial vibrational relaxation probability to all elements
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals                ,ONLY: abort,UNIT_stdOut, IK, MPI_COMM_WORLD, MPIRoot
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
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
USE MOD_DSMC_Collis            ,ONLY: DSMC_calc_var_P_vib
USE MOD_part_emission_tools    ,ONLY: CalcVelocity_maxwell_lpn
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
            Ti = Species(iSpec)%Init(0)%MWTemperatureIC
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
              Tj = Species(jSpec)%Init(0)%MWTemperatureIC
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
IF(DSMC%VibRelaxProb.EQ.2) THEN
  SDEALLOCATE(VarVibRelaxProb%ProbVibAv)
  SDEALLOCATE(VarVibRelaxProb%ProbVibAvNew)
  SDEALLOCATE(VarVibRelaxProb%nCollis)
  SDEALLOCATE(DSMC%QualityFacSampRot)
  SDEALLOCATE(DSMC%QualityFacSampVib)
  SDEALLOCATE(DSMC%QualityFacSampRotSamp)
  SDEALLOCATE(DSMC%QualityFacSampVibSamp)
  SDEALLOCATE(DSMC%CalcVibProb)
  SDEALLOCATE(DSMC%CalcRotProb)
END IF
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
SDEALLOCATE(DSMC%veloMinColl)
SDEALLOCATE(DSMC%NumColl)
SDEALLOCATE(DSMC%InstantTransTemp)
IF(DSMC%CalcQualityFactors) THEN
  SDEALLOCATE(DSMC%QualityFacSamp)
END IF
SDEALLOCATE(PDM%PartInit)
SDEALLOCATE(Coll_pData)
SDEALLOCATE(SampDSMC)
SDEALLOCATE(MacroDSMC)
SDEALLOCATE(QKAnalytic)
SDEALLOCATE(ChemReac%QKProcedure)
SDEALLOCATE(ChemReac%QKMethod)
SDEALLOCATE(ChemReac%QKCoeff)
SDEALLOCATE(ChemReac%NumReac)
SDEALLOCATE(ChemReac%ReacCount)
SDEALLOCATE(ChemReac%ReacCollMean)
SDEALLOCATE(ChemReac%ReacCollMeanCount)
SDEALLOCATE(ChemReac%NumReac)
SDEALLOCATE(ChemReac%ReactType)
SDEALLOCATE(ChemReac%DefinedReact)
SDEALLOCATE(ChemReac%ReactCase)
SDEALLOCATE(ChemReac%ReactNum)
SDEALLOCATE(ChemReac%Arrhenius_Prefactor)
SDEALLOCATE(ChemReac%Arrhenius_Powerfactor)
SDEALLOCATE(ChemReac%EActiv)
SDEALLOCATE(ChemReac%EForm)
SDEALLOCATE(ChemReac%MeanEVib_PerIter)
SDEALLOCATE(ChemReac%MeanEVibQua_PerIter)
SDEALLOCATE(ChemReac%CEXa)
SDEALLOCATE(ChemReac%CEXb)
SDEALLOCATE(ChemReac%MEXa)
SDEALLOCATE(ChemReac%MEXb)
SDEALLOCATE(ChemReac%ELa)
SDEALLOCATE(ChemReac%ELb)
SDEALLOCATE(ChemReac%DoScat)
SDEALLOCATE(ChemReac%ReactInfo)
SDEALLOCATE(ChemReac%TLU_FileName)
SDEALLOCATE(ChemReac%ReactNumRecomb)
SDEALLOCATE(ChemReac%Hab)
SDEALLOCATE(CollInf%Coll_Case)
SDEALLOCATE(CollInf%Coll_CaseNum)
SDEALLOCATE(CollInf%Coll_SpecPartNum)
SDEALLOCATE(CollInf%crossSectionConstantCab)
SDEALLOCATE(CollInf%KronDelta)
SDEALLOCATE(CollInf%FracMassCent)
SDEALLOCATE(CollInf%MassRed)
SDEALLOCATE(CollInf%MeanMPF)
SDEALLOCATE(CollInf%alphaVSS) 
SDEALLOCATE(CollInf%omegaLaux)    
SDEALLOCATE(CollInf%dref)    
SDEALLOCATE(CollInf%Tref)    
SDEALLOCATE(HValue)
!SDEALLOCATE(SampWall)
SDEALLOCATE(MacroSurfaceVal)
!SDEALLOCATE(VibQuantsPar)
! SDEALLOCATE(XiEq_Surf)
SDEALLOCATE(DSMC_HOSolution)
SDEALLOCATE(ElemNodeVol)
SDEALLOCATE(BGGas%PairingPartner)
SDEALLOCATE(RadialWeighting%ClonePartNum)
SDEALLOCATE(ClonedParticles)
SDEALLOCATE(SymmetrySide)
END SUBROUTINE FinalizeDSMC


END MODULE MOD_DSMC_Init
