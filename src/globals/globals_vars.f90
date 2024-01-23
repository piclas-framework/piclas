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

!===================================================================================================================================
!> Provides parameters, used globally (please use EXTREMELY carefully!)
!===================================================================================================================================
MODULE MOD_Globals_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=6),PARAMETER :: ProgramName  = 'PICLas'              !> name of this program
INTEGER,PARAMETER          :: MajorVersion = 3                     !> FileVersion number saved in each hdf5 file with hdf5 header
INTEGER,PARAMETER          :: MinorVersion = 1                     !> FileVersion number saved in each hdf5 file with hdf5 header
INTEGER,PARAMETER          :: PatchVersion = 0                     !> FileVersion number saved in each hdf5 file with hdf5 header
REAL,PARAMETER             :: FileVersionReal  = REAL(MajorVersion,8)+REAL(MinorVersion,8)/10.+REAL(PatchVersion,8)/100.
                                                                   !> OLD number saved in each hdf5 file with hdf5 header
INTEGER,PARAMETER          :: FileVersionInt = PatchVersion+MinorVersion*100+MajorVersion*10000
                                                                   !>  NEWnumber saved in each hdf5 file with hdf5 header
CHARACTER(LEN=10)          :: PiclasVersionStr                     !> PiclasVersionStrnumber saved in each hdf5 file with hdf5 header
REAL                       :: FileVersionHDF5Real                  !> OLD FileVersion number read from hdf5 restart file
REAL                       :: FileVersionHDF5Int                   !> NEW FileVersion number read from hdf5 restart file
REAL                       :: WallTime                             !> Wall time needed by a simulation (is not reset by
                                                                   !> performing a load balance step, only by user restart)
REAL                       :: InitializationWallTime               !> Wall time needed to initialize a simulation (or
                                                                   !> re-initialize a simulation by performing a load balance
                                                                   !>  step)
REAL                       :: ReadMeshWallTime                     !> Wall time needed to read the mesh (SUBROUTINE ReadMesh)
REAL                       :: DomainDecompositionWallTime          !> Wall time needed for domain decomposition
REAL                       :: CommMeshReadinWallTime               !> Shared memory mesh communication
REAL                       :: SimulationEfficiency                 !> relates the simulated time to the used CPUh (SIMULATION TIME PER
                                                                   !> CALCULATION in [s]/[CPUh])
REAL                       :: StartT                               !> Timer start
REAL                       :: PID                                  !> Performance index: (CalcTimeEnd-CalcTimeStart)*nProcessors/(nGlobalElems*(PP_N+1)**3*iter_loc)
REAL                       :: memory(1:4)                          !> RAM: used, available, total and initial (total at the beginning of the simulation)
REAL,PARAMETER             :: PI=ACOS(-1.0)                         !> the number pi ~= 3.14
REAL,PARAMETER             :: sPI=1.0/PI                            !> inverse of pi
REAL,PARAMETER             :: epsMach=EPSILON(0.0)                  !> Machine accuracy
REAL,PARAMETER             :: TwoepsMach=2.0d0*epsMach              !> twice the machine accuracy
REAL,PARAMETER             :: EuMas          = 0.577215664901533_8  !> Euler-Mascheroni constant
REAL,PARAMETER             :: PlanckConst    = 6.62606957e-34       !> Planck constant [J s] SI-Unit!
REAL,PARAMETER             :: ElementaryCharge = 1.602176634e-19    !> redefinition of SI base units in 2018-2019,
                                                                    !> => negative charge of an electron, eV to Joule, ...
REAL,PARAMETER             :: StefanBoltzmannConst = 5.670374419E-8
REAL,PARAMETER             :: ElectronMass   = 9.1093826e-31        !> mass of an electron
CHARACTER(LEN=255)         :: ProjectName                           !> TODO-DEFINE-PARAMETER
CHARACTER(LEN=255)         :: ParameterFile                         !> filename of the parameter file
CHARACTER(LEN=255)         :: ParameterDSMCFile                     !> filename of the parameterDSMC file
REAL, PARAMETER            :: BoltzmannConst=1.380648813e-23        !> Boltzmann constant [J/K] SI-Unit! in m^2/(s^2*K)
REAL, PARAMETER            :: Kelvin2eV=8.61732814974056e-5         !> Conversion factor [K]  -> [eV] (Kelvin to electron volt)
REAL, PARAMETER            :: eV2Kelvin=1.16045250061657e4          !> Conversion factor [eV] -> [K]  (electron volt to Kelvin)
REAL, PARAMETER            :: Joule2eV=6.241506363094e+18           !> Conversion factor [J]  -> [eV] (Joule to electron volt)
REAL, PARAMETER            :: eV2Joule=1.60217734e-19               !> Conversion factor [eV] -> [J]  (electron volt to Joule)
CHARACTER(LEN=5)           :: TimeStampLenStr,TimeStampLenStr2      !> Strings for timestamp format of time
REAL,PARAMETER             :: BohrRadius     = 5.2917721067E-11     !> Radius, 1st Bohr orbit for H (a0) [m]
REAL,PARAMETER             :: AtomicMassUnit = 1.660539040E-27      !> Atomic mass unit [kg]

REAL,PARAMETER             :: maxEXP= LOG(HUGE(maxexp))
! Set variables (natural constants and derived quantities) from user input or hard coded
! depending on compile flag (PICLAS_READIN_CONSTANTS=ON)
#if USE_READIN_CONSTANTS
REAL           :: eps0                        !> permittivity eps0
REAL           :: mu0                         !> permeability mu0
REAL           :: smu0                        !> 1/mu0
REAL           :: c                           !> speed of light c
REAL           :: c2                          !> c^2
REAL           :: c2_inv                      !> 1/c^2
REAL           :: c_inv                       !> 1/c
REAL           :: RelativisticLimit           !> for comparison with velocity^2 and is automatically set to 0.3% speed of light
#else
REAL,PARAMETER :: eps0   = 8.8541878176e-12   !> permittivity eps0 of vacuum [F/m]
REAL,PARAMETER :: mu0    = 1.2566370614e-6    !> permeability mu0 of vacuum [H/m]
REAL,PARAMETER :: smu0   = 1/mu0              !> 1/mu0 = 7.9577471548222157e5 [m/H]
REAL,PARAMETER :: c      = 299792458.0        !> speed of light c in vacuum [m/s]
REAL,PARAMETER :: c2     = c**2               !> c^2   = 8.9875517873681764e16 [m^2/s^2]
REAL,PARAMETER :: c2_inv = 1/c2               !> 1/c^2 = 1.1126500560536184e-17 [s^2/m^2]
REAL,PARAMETER :: c_inv  = 1/c                !> 1/c   = 3.3356409519815204e-9 [s/m]
REAL,PARAMETER :: RelativisticLimit = 1e12    !> for comparison with velocity^2 when speed of light is 299792458. Corresponds to 0.3% speed of light
#endif /*USE_READIN_CONSTANTS*/
!===================================================================================================================================

!CONTAINS

END MODULE MOD_Globals_Vars
