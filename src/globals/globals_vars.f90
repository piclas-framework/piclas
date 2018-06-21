#include "boltzplatz.h"

!===================================================================================================================================
!> Provides parameters, used globally (please use EXTREMLY carefully!) 
!===================================================================================================================================
MODULE MOD_Globals_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
REAL                         :: SimulationTime                        !> Wall time needed by a simulation (is not reset by 
                                                                      !> performing a load balance step, only by user restart)
REAL                         :: InitializationWallTime                !> Wall time needed to initialize a simulation (or
                                                                      !> re-initialize a simulation by performing a load balance 
                                                                      !>  step)
REAL                         :: SimulationEfficiency                  !> relates the simulated time to the used CPUh (SIMULATION TIME PER
                                                                      !> CALCULATION in [s]/[CPUh])
REAL                         :: PID                                   !> Performance index: (CalcTimeEnd-CalcTimeStart)*nProcessors/
                                                                      !> (nGlobalElems*(PP_N+1)**3*iter_loc)
REAL                         :: PI                                    !> the number pi ~= 3.14
REAL                         :: sPI                                   !> inverse of pi
REAL                         :: epsMach,TwoepsMach                    !> TODO-DEFINE-PARAMETER
REAL,PARAMETER               :: EuMas          = 0.577215664901533_8  !> Euler-Mascheroni constant
REAL,PARAMETER               :: PlanckConst    = 6.62606957E-34       !> Planck constant [J s] SI-Unit!
REAL,PARAMETER               :: ElectronCharge = 1.60217653e-19       !> charge of an electron
REAL,PARAMETER               :: ElectronMass   = 9.1093826E-31        !> mass of an electron
REAL,PARAMETER               :: FileVersion    = 0.11                 !> FileVersion number saved in each hdf5 file with hdf5 header
CHARACTER(LEN=255),PARAMETER :: ProgramName    = 'Boltzplatz'         !> name of this program
CHARACTER(LEN=255)           :: ProjectName                           !> TODO-DEFINE-PARAMETER
CHARACTER(LEN=255)           :: ParameterFile                         !> filename of the parameter file
CHARACTER(LEN=255)           :: ParameterDSMCFile                     !> filename of the parameterDSMC file
#ifndef PARTICLES
REAL, PARAMETER              :: BoltzmannConst=1.380648813E-23        !> Boltzmann constant [J/K] SI-Unit! in m^2/(s^2*K)
#endif
!===================================================================================================================================

!CONTAINS

END MODULE MOD_Globals_Vars
