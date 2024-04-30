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
REAL                          :: ParticleAnalyzeSampleTime           !< Accumulated simulation time between two outputs to ParticleAnalyze.csv
LOGICAL                       :: CalcSimNumSpec                      !< Calculate the number of simulated particles per species
LOGICAL                       :: CalcNumDens                         !< Calculate the number density per species within the domain
LOGICAL                       :: CalcSurfFluxInfo                    !< Calculate the current/mass flow through or pressure (adaptive/subsonic BC) at the surface flux boundaries
LOGICAL                       :: CalcCollRates                       !< Calculate the collision rates per collision pair
LOGICAL                       :: CalcReacRates                       !< Calculate the reaction rate per reaction
LOGICAL                       :: CalcRelaxProb                       !< Calculate relaxation probabilities
LOGICAL                       :: CalcEkin                            !< Compute the kinetic energy of each species
LOGICAL                       :: CalcEtot                            !< Compute the total energy as sum of potential and kin eng
LOGICAL                       :: CalcEint(2)                         !< Compute the internal energy of each species [1: Calculate, 2: Output]
LOGICAL                       :: CalcTemp(2)                         !< Computation of the temperature (trans, rot, vib, total)
LOGICAL                       :: CalcCoupledPower                    !< Computation of the power that is coupled into plasma
LOGICAL                       :: DisplayCoupledPower                 !< Display coupled power in UNIT_stdOut
REAL                          :: EDiff                               !< Difference in kinetic energy before and after the particle
                                                                     !< push (only charged particles)
REAL                          :: PCoupl                              !< Power that is coupled into plasma in [W]
REAL                          :: PCouplAverage                       !< Power that is coupled into plasma (moving average) in [W]
REAL                          :: PCouplIntAverage                    !< Power that is coupled into plasma (moving integrated average) in [W]
REAL                          :: PCouplAverageOld                    !< Power that is coupled into plasma (moving integrated average) in [W] - old value from last call to PartAnalyze()
TYPE tPCoupl
  REAL,ALLOCATABLE            :: DensityAvgElem(:)                   !< Power per volume that is coupled into plasma (moving average
                                                                     !< for each element) in [W/m^3]
END TYPE
TYPE(tPCoupl),ALLOCATABLE     :: PCouplSpec(:)                       !< DensityAvgElem array for each species



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
LOGICAL                       :: DoPartAnalyze                       !< perform analyze
INTEGER(KIND=8)               :: PartAnalyzeStep                     !< Analyze is performed each Nth time step
INTEGER,ALLOCATABLE           :: nPartIn(:)                          !< Number of entry and leaving particles
INTEGER,ALLOCATABLE           :: nPartOut(:)                         !< Number of entry and leaving particles
REAL,ALLOCATABLE              :: PartEkinIn(:)                       !< Energy and temperature of input particle
REAL,ALLOCATABLE              :: PartEkinOut(:)                      !< Energy and temperature of input particle

! get derived particle properties (for IMD/TTM initialization these values are calculated from the TTM grid values)
LOGICAL                       :: CalcDebyeLength                     !< Compute the Debye length (min and max) in each cell
LOGICAL                       :: CalcPICTimeStep                     !< Compute the PIC time step from plasma frequency (min and max) in each cell
LOGICAL                       :: CalcPICTimeStepCyclotron            !< Compute the PIC time step from cyclotron motion (min and max) in each cell
LOGICAL                       :: CalcElectronIonDensity              !< Compute the electron density in each cell
LOGICAL                       :: CalcElectronTemperature             !< Compute the electron temperature in each cell
LOGICAL                       :: CalcElectronEnergy                  !< Compute the electron min/max/average energy in each cell
LOGICAL                       :: CalcPlasmaParameter                 !< Compute the plasma parameter in each cell
!LOGICAL                       :: ElectronTemperatureIsMaxwell        ! Assumption of Maxwell-Boltzmann or undistributed electrons
LOGICAL                       :: CalcPlasmaFrequency                 !< Compute the electron plasma frequency in each cell
LOGICAL                       :: CalcCyclotronFrequency              !< Compute the electron cyclotron frequency in each cell
!                                                                    !< (requires a magnetic field)
!LOGICAL                       :: CalcGyroradius                      !< Compute the electron cyclotron radius from the particle velocity and yclotron frequency in each cell
LOGICAL                       :: CalcPointsPerDebyeLength            !< Compute the points per Debye length:
LOGICAL                       :: CalcPICCFLCondition                 !< Compute a PIC CFL condition for each cell
!                                                                    !< in terms of cell lengths in X, Y and Z for each cell
!                                                                    !< PPD=(p+1)lambda_D/L_cell
LOGICAL                       :: CalcMaxPartDisplacement             !< Compute the maximum displacement of the fastest particle
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
INTEGER                       :: PPDCellResolved(4)                  !> Number of elements with PPD>=1
INTEGER                       :: PICTimeCellResolved                 !> Number of cells where the time step is resolved
REAL,ALLOCATABLE              :: PPDCell(:)                          !< Points per Debye length (cell mean value)
REAL,ALLOCATABLE              :: PPDCellX(:)                         !< Points per Debye length in X (cell mean value)
REAL,ALLOCATABLE              :: PPDCellY(:)                         !< Points per Debye length in Y (cell mean value)
REAL,ALLOCATABLE              :: PPDCellZ(:)                         !< Points per Debye length in Z (cell mean value)
REAL,ALLOCATABLE              :: PICCFLCell(:)                       !< PIC CFL Condition (cell mean value)
REAL,ALLOCATABLE              :: PICCFLCellX(:)                      !< PIC CFL Condition in X (cell mean value)
REAL,ALLOCATABLE              :: PICCFLCellY(:)                      !< PIC CFL Condition in Y (cell mean value)
REAL,ALLOCATABLE              :: PICCFLCellZ(:)                      !< PIC CFL Condition in Z (cell mean value)
REAL,ALLOCATABLE              :: MaxPartDisplacementCell(:)          !< Maximum particle displacement (cell mean value)
REAL,ALLOCATABLE              :: MaxPartDisplacementCellX(:)         !< Maximum particle displacement in X (cell mean value)
REAL,ALLOCATABLE              :: MaxPartDisplacementCellY(:)         !< Maximum particle displacement in Y (cell mean value)
REAL,ALLOCATABLE              :: MaxPartDisplacementCellZ(:)         !< Maximum particle displacement in Z (cell mean value)
REAL,ALLOCATABLE              :: PPSCell(:)                          !< Points per shape function sphere (cell mean value):
                                                                     !<   calculate cell local number excluding neighbor DOFs
REAL,ALLOCATABLE              :: PPSCellCartesian(:)                 !< Points per shape function sphere (cell mean value):
                                                                     !<   assume Cartesian grid and calculate to total number
                                                                     !<   including neighbor DOFs
REAL,ALLOCATABLE              :: ShapeFunctionRadius(:)              !< Additional array (shape function radius is already stored in
                                                                     !< the shared array) for output to .h5 (debugging)
REAL,ALLOCATABLE              :: ShapeFunctionFraction(:)            !< Element to shape function volume ratio
REAL,ALLOCATABLE              :: DebyeLengthCell(:)                  !< Debye length (cell mean value)
REAL,ALLOCATABLE              :: PICTimeStepCell(:)                  !< Approximated PIC Time Step due to plasma frequency (mean cell value)
REAL,ALLOCATABLE              :: PICTimeStepCyclotronCell(:)         !< Approximated PIC Time Step due to cyclotron frequency (mean cell value)
REAL,ALLOCATABLE              :: PlasmaParameterCell(:)              !< Plasma parameter (cell mean value)
REAL,ALLOCATABLE              :: ElectronDensityCell(:)              !< Electron density (cell mean value)
REAL,ALLOCATABLE              :: IonDensityCell(:)                   !< Ion density (cell mean value)
REAL,ALLOCATABLE              :: NeutralDensityCell(:)               !< Neutral density (cell mean value)
REAL,ALLOCATABLE              :: ChargeNumberCell(:)                 !< Charge number (cell mean value)
INTEGER,ALLOCATABLE           :: PICValidPlasmaCell(:)               !< Check that quasi-neutrality is above 0.5 and at least 20 particles are inside the element
INTEGER                       :: PICValidPlasmaCellSum               !< Global number of elements that have quasi-neutrality above 0.5 and at least 20 particles are inside the element
REAL,ALLOCATABLE              :: ElectronTemperatureCell(:)          !< Electron temperature (cell mean value)
REAL,ALLOCATABLE              :: ElectronMinEnergyCell(:)            !< Electron minimum cell energy [eV]
REAL,ALLOCATABLE              :: ElectronMaxEnergyCell(:)            !< Electron maximum cell energy [eV]
REAL,ALLOCATABLE              :: ElectronAverageEnergyCell(:)        !< Electron average cell energy [eV]
REAL,ALLOCATABLE              :: PlasmaFrequencyCell(:)              !< Plasma electron frequency (cell mean value)
REAL,ALLOCATABLE              :: CyclotronFrequencyMaxCell(:)        !< Electron cyclotron frequency (cell MAX value)
REAL,ALLOCATABLE              :: CyclotronFrequencyMinCell(:)        !< Electron cyclotron frequency (cell MIN value)
REAL,ALLOCATABLE              :: GyroradiusMinCell(:)                !< Electron gyroradius (cyclotron or Larmor radius), from frequency (cell MIN value)
REAL,ALLOCATABLE              :: GyroradiusMaxCell(:)                !< Electron gyroradius (cyclotron or Larmor radius), from frequency (cell MAX value)

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
REAL,ALLOCATABLE              :: FlowRateSurfFlux(:,:)               !< Particle balance per surface flux BC, utilized to calculate mass flog or current
REAL,ALLOCATABLE              :: PressureAdaptiveBC(:,:)
LOGICAL                       :: CalcEMFieldOutput                   !< Output the electro-magnetic fields on each DOF to .h5 calculated by PIC interpolation external fields and from field solver
!===================================================================================================================================
END MODULE MOD_Particle_Analyze_Vars
