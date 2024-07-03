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


MODULE MOD_Transport_Data
! MODULES
IMPLICIT NONE
PRIVATE

!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC::CalcDriftDiffusionCoeff, InterpolateCoefficient
!==================================================================================================================================
CONTAINS

SUBROUTINE CalcDriftDiffusionCoeff(ElectricField,Density,mu,D)
!==================================================================================================================================
!> Calculate the transport (drift & diffusion) coefficients for the drift-diffusion electron fluid model
!==================================================================================================================================
! MODULES
USE MOD_PreProc ! PP_N
USE MOD_DSMC_Vars          ,ONLY: BGGas
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                                  :: ElectricField, Density
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                                 :: mu, D
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                             :: ReducedElectricField
!===================================================================================================================================

IF (Density.GT.0.) THEN
    ReducedElectricField=1.e21*ElectricField/Density ! E/n in Townsend as defined in LXCAT database
    mu = InterpolateCoefficient(BGGas%ElectronMobility,ReducedElectricField)
    D = InterpolateCoefficient(BGGas%DriftDiffusionCoefficient,ReducedElectricField)
    ! (actually not energy but kT/e so eV is already SI here)
END IF

END SUBROUTINE CalcDriftDiffusionCoeff


PPURE REAL FUNCTION InterpolateCoefficient(CoeffData,ReducedElectricField)
!===================================================================================================================================
!> Interpolate the collision coefficient [m^2] from the available data at the given electric field [J]
!> Collision energies below and above the given data will be set at the first and last level of the data set
!> Note: Requires the data to be sorted by ascending energy values
!> Assumption: First species given is the particle species, second species input is the background gas species
!===================================================================================================================================
! MODULES
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,INTENT(IN)               :: CoeffData(:,:)                    !< Array for the interpolation [1:2,1:MaxDOF]
REAL,INTENT(IN)               :: ReducedElectricField                 !< electric field in [J]
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iDOF, MaxDOF
!===================================================================================================================================

InterpolateCoefficient = 0.
MaxDOF = SIZE(CoeffData,2)

IF(ReducedElectricField.GT.CoeffData(1,MaxDOF)) THEN
  ! If the electric field is greater than the maximal value, extrapolate from the last two values
  IF((MaxDOF.LT.2).OR.(CoeffData(2,MaxDOF).LE.0.))THEN
    ! If only one value is given or the last coefficient is zero
    InterpolateCoefficient = CoeffData(2,MaxDOF)
  ELSE
    ! Extrapolate
    InterpolateCoefficient = CoeffData(2,MaxDOF-1) + (ReducedElectricField - CoeffData(1,MaxDOF-1)) &
                            / (CoeffData(1,MaxDOF) - CoeffData(1,MaxDOF-1)) * (CoeffData(2,MaxDOF) - CoeffData(2,MaxDOF-1))
    ! Check if extrapolation drops under zero
    IF(InterpolateCoefficient.LE.0.) InterpolateCoefficient = 0.
  END IF ! (MaxDOF.LT.2).OR.(CoeffData(2,MaxDOF).LE.0.)
  ! Leave routine
  RETURN
ELSE IF(ReducedElectricField.LE.CoeffData(1,1)) THEN
  ! If electric field is below the minimal value, get the coefficient of the first level and leave routine
  InterpolateCoefficient = CoeffData(2,1)
  ! Leave routine
  RETURN
END IF

DO iDOF = 1, MaxDOF
  ! Check if the stored energy value is above the electric field
  IF(CoeffData(1,iDOF).GT.ReducedElectricField) THEN
    ! Interpolate the coefficient from the data set using the current and the energy level below
    InterpolateCoefficient = CoeffData(2,iDOF-1) + (ReducedElectricField  - CoeffData(1,iDOF-1)) &
                              / (CoeffData(1,iDOF) - CoeffData(1,iDOF-1)) * (CoeffData(2,iDOF) - CoeffData(2,iDOF-1))
    ! Leave routine and do not finish DO loop
    RETURN
  ELSE IF(ReducedElectricField.EQ.CoeffData(1,iDOF)) THEN
    ! In case the electric field is exactly the stored value
    InterpolateCoefficient = CoeffData(2,iDOF)
  END IF
END DO

END FUNCTION InterpolateCoefficient


END MODULE MOD_Transport_Data
