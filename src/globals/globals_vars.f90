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
CHARACTER(LEN=6),PARAMETER   :: ProgramName    = 'PICLas'             !> name of this program
REAL,PARAMETER               :: FileVersion    = 1.3                  !> FileVersion number saved in each hdf5 file with hdf5 header
REAL                         :: WallTime                              !> Wall time needed by a simulation (is not reset by
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
REAL,PARAMETER               :: ElementaryCharge = 1.602176634e-19    !> redefinition of SI base units in 2018-2019,
                                                                      !> => negative charge of an electron, joule to eV, ...
REAL,PARAMETER               :: ElectronMass   = 9.1093826E-31        !> mass of an electron
CHARACTER(LEN=255)           :: ProjectName                           !> TODO-DEFINE-PARAMETER
CHARACTER(LEN=255)           :: ParameterFile                         !> filename of the parameter file
CHARACTER(LEN=255)           :: ParameterDSMCFile                     !> filename of the parameterDSMC file
REAL, PARAMETER              :: BoltzmannConst=1.380648813E-23        !> Boltzmann constant [J/K] SI-Unit! in m^2/(s^2*K)
CHARACTER(LEN=5)             :: TimeStampLenStr,TimeStampLenStr2      !> Strings for timestamp format of time
!===================================================================================================================================

!CONTAINS

END MODULE MOD_Globals_Vars
