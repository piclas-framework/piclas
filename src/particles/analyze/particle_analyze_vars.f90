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
MODULE MOD_Particle_Analyze_Vars
!===================================================================================================================================
!> Contains global variables used by the Analyze modules.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                       :: ParticleAnalyzeInitIsDone = .FALSE.
LOGICAL                       :: CalcNumSpec                         !< Calculate the number of simulated particles per species
LOGICAL                       :: CalcNumDens                         !< Calculate the number density per species within the domain
LOGICAL                       :: CalcCollRates                       !< Calculate the collision rates per collision pair
LOGICAL                       :: CalcReacRates                       !< Calculate the reaction rate per reaction
LOGICAL                       :: CalcEkin                            !< Compute the kinetic energy of each species
LOGICAL                       :: CalcEtot                            !< Compute the total energy as sum of potential and kin eng
LOGICAL                       :: CalcEint                            !< Compute the internal energy of each species
LOGICAL                       :: CalcTemp                            !< Computation of the temperature (trans, rot, vib, total)
LOGICAL                       :: CalcCoupledPower                    !< Computation of the power that is coupled into plasma
REAL                          :: PCoupl                              !< Power that is coupled into plasma
REAL                          :: PCouplAverage                       !< Power that is coupled into plasma (moving average)
LOGICAL                       :: CalcPartBalance                     !< Particle Power Balance - input and outflow energy of all
                                                                     !< particles
LOGICAL                       :: CalcVelos                           !< Computes the drift and thermal velocity of each species
LOGICAL                       :: VeloDirs(4)                         !< Select the direction for velocity computation
LOGICAL                       :: TrackParticlePosition               !< Track the particle movement
                                                                     !< Stored in .csv format, debug only, no MPI
INTEGER                       :: nSpecAnalyze                        !< Number of analyzed species 1 or nSpecies+1
LOGICAL                       :: IsRestart                           !< Check if restart, add data to Database
LOGICAL                       :: ChargeCalcDone                      !< Check flag
LOGICAL                       :: CalcShapeEfficiency                 !< Efficiency of shape function
CHARACTER(LEN=256)            :: CalcShapeEfficiencyMethod           !< Explanations in particle_analyze.f90
INTEGER                       :: ShapeEfficiencyNumber               !< Explanations in particle_analyze.f90
INTEGER                       :: FieldAnalyzeStep                    !< Analyze is performed each Nth time step
LOGICAL                       :: DoPartAnalyze                       !< perform analyze
INTEGER                       :: PartAnalyzeStep                     !< Analyze is performed each Nth time step
INTEGER,ALLOCATABLE           :: nPartIn(:)                          !< Number of entry and leaving particles
INTEGER,ALLOCATABLE           :: nPartOut(:)                         !< Number of entry and leaving particles
INTEGER,ALLOCATABLE           :: nPartInTmp(:)                       !< Number of entry and leaving particles
REAL,ALLOCATABLE              :: PartEkinIn(:)                       !< Energy and temperature of input particle
REAL,ALLOCATABLE              :: PartEkinOut(:)                      !< Energy and temperature of input particle
REAL,ALLOCATABLE              :: PartEKinInTmp(:)                    !< Energy and temperature of input particle

! get derived particle properties (for IMD/TTM initialization these values are calculated from the TTM grid values)
LOGICAL                       :: CalcDebyeLength                     !< Compute the Debye length (min and max) in each cell
LOGICAL                       :: CalcPICTimeStep                     !< Compute the PIC time step (min and max) in each cell
LOGICAL                       :: CalcElectronIonDensity              !< Compute the electron density in each cell
LOGICAL                       :: CalcElectronTemperature             !< Compute the electron temperature in each cell
LOGICAL                       :: CalcPlasmaParameter                 !< Compute the plasma parameter in each cell
!LOGICAL                       :: ElectronTemperatureIsMaxwell        ! Assumption of Maxwell-Boltzmann or undistributed electrons
LOGICAL                       :: CalcPlasmaFrequency                 !< Compute the electron frequency in each cell
LOGICAL                       :: CalcPointsPerDebyeLength            !< Compute the points per Debye length:
!                                                                    !< PPD=(p+1)lambda_D/L_cell
LOGICAL                       :: CalcPointsPerShapeFunction          !< Compute the points per shape function sphere
!                                                                    !< PPS = DOF_cell*VolumeShapeFunction/Volume_cell

LOGICAL                       :: CalcIonizationDegree                !< Compute the ionization degree and quasi neutrality
!                                                                    !< in each cell
LOGICAL                       :: CalcLaserInteraction                !< Compute laser-plasma interaction properties such as maximum
REAL                          :: LaserInteractionEkinMaxRadius       !< maximum radius (x- and y-dir) of particle to be considered
!                                                                    !< for Ekin maximum calculation (default is HUGE)
!                                                                    !< OR LaserInteractionEkinMaxZPosMin
REAL                          :: LaserInteractionEkinMaxZPosMin      !< minimum z-position of particle to be considered for Ekin
!                                                                    !< maximum calculation (default is -1.*HUGE)
!                                                                    !< OR LaserInteractionEkinMaxRadius
!                                                                    !<particle energy per species. Default=.FALSE.
REAL,ALLOCATABLE              :: IonizationCell(:)                   !< Ionization degree cell value
REAL,ALLOCATABLE              :: QuasiNeutralityCell(:)              !< QuasiNeutrality degree cell value
REAL,ALLOCATABLE              :: PPDCell(:)                          !< Points per Debye length (cell mean value)
REAL,ALLOCATABLE              :: PPDCellX(:)                         !< Points per Debye length in X (cell mean value)
REAL,ALLOCATABLE              :: PPDCellY(:)                         !< Points per Debye length in Y (cell mean value)
REAL,ALLOCATABLE              :: PPDCellZ(:)                         !< Points per Debye length in Z (cell mean value)
REAL,ALLOCATABLE              :: PPSCell(:)                          !< Points per shape function sphere (cell mean value):
                                                                     !<   calculate cell local number excluding neighbor DOFs
REAL,ALLOCATABLE              :: PPSCellEqui(:)                      !< Points per shape function sphere (cell mean value):
                                                                     !<   assume Cartesian grid and calculate to total number
                                                                     !<   including neighbor DOFs
REAL,ALLOCATABLE              :: DebyeLengthCell(:)                  !< Debye length (cell mean value)
REAL,ALLOCATABLE              :: PICTimeStepCell(:)                  !< Approximated PIC Time Step (mean cell value)
REAL,ALLOCATABLE              :: PlasmaParameterCell(:)              !< Approximated PIC Time Step (mean cell value)
REAL,ALLOCATABLE              :: ElectronDensityCell(:)              !< Electron density (cell mean value)
REAL,ALLOCATABLE              :: IonDensityCell(:)                   !< Ion density (cell mean value)
REAL,ALLOCATABLE              :: NeutralDensityCell(:)               !< Neutral density (cell mean value)
REAL,ALLOCATABLE              :: ChargeNumberCell(:)                 !< Charge number (cell mean value)
REAL,ALLOCATABLE              :: ElectronTemperatureCell(:)          !< Electron temperature (cell mean value)
REAL,ALLOCATABLE              :: PlasmaFrequencyCell(:)              !< Plasma electron frequency (cell mean value)

LOGICAL                       :: CalcCharge                          !< Compute the whole deposited charge and abs and relative
                                                                     !< Charge error
LOGICAL                       :: DoVerifyCharge                      !< Validate the charge after each deposition and produces
                                                                     !< an output in std.out
REAL                          :: PartCharge(3)                       !< Contains the whole deposited charge and its absolute
                                                                     !< and relative error
LOGICAL                       :: printDiff                           !< TODO
REAL                          :: printDiffTime                       !< TODO
REAL                          :: printDiffVec(6)                     !< TODO
REAL                          :: ChemEnergySum                       !< TODO
LOGICAL                       :: CalcPorousBCInfo                    !< Calculate output for porous BCs (averaged over whole BC)
!===================================================================================================================================
END MODULE MOD_Particle_Analyze_Vars
