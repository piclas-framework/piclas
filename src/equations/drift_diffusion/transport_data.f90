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

PUBLIC::CalcDriftDiffusionCoeff
!==================================================================================================================================
CONTAINS

SUBROUTINE CalcDriftDiffusionCoeff(ElectricField,Density,mu,D)
!==================================================================================================================================
!> Calculate the transport (drift & diffusion) coefficients for the drift-diffusion electron fluid model
!> Now using data for electrons in H2 gas from LXCat database (Klaus Berkhan, PhD thesis, 1994, University Heidelberg)
!==================================================================================================================================
! MODULES
USE MOD_PreProc ! PP_N
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
REAL                                             :: EnergyTable(2,26), muTable(2,56), ReducedElectricField
!===================================================================================================================================
mu= 0.
D = 0.

EnergyTable(1,:) = (/0.04, &
0.05, &
0.06, &
0.07, &
0.08, &
0.1, &
0.12, &
0.14, &
0.17, &
0.2, &
0.25, &
0.3, &
0.35, &
0.4, &
0.5, &
0.6, &
0.7, &
0.8, &
1., &
1.2, &
1.4, &
1.7, &
2., &
2.5, &
3., &
3.5/)

EnergyTable(2,:) =(/0.026, &
0.0258, &
0.0256, &
0.0261, &
0.026, &
0.0263, &
0.0268, &
0.0272, &
0.0278, &
0.0287, &
0.03, &
0.0314, &
0.0334, &
0.0342, &
0.0366, &
0.0388, &
0.0417, &
0.045, &
0.0488, &
0.0546, &
0.0596, &
0.0678, &
0.0763, &
0.0923, &
0.1071, &
0.1208/)

muTable(1,:) =(/0.008, &
0.01, &
0.012, &
0.014, &
0.017, &
0.02, &
0.025, &
0.03, &
0.035, &
0.04, &
0.05, &
0.06, &
0.07, &
0.08, &
0.1, &
0.12, &
0.14, &
0.17, &
0.2, &
0.25, &
0.3, &
0.35, &
0.4, &
0.5, &
0.6, &
0.7, &
0.8, &
1., &
1.2, &
1.4, &
1.7, &
2., &
2.5, &
3., &
3.5, &
4., &
5., &
6., &
7., &
8., &
10., &
12., &
14., &
17., &
20., &
25., &
30., &
35., &
40., &
45., &
50., &
55., &
60., &
65., &
70., &
71./)

muTable(2,:) =(/1.53625E+025, &
1.536E+025, &
1.5325E+025, &
1.53E+025, &
1.52412E+025, &
1.5185E+025, &
1.5096E+025, &
1.49833E+025, &
1.48886E+025, &
1.473E+025, &
1.4518E+025, &
1.42517E+025, &
1.40114E+025, &
1.3775E+025, &
1.332E+025, &
1.29083E+025, &
1.25357E+025, &
1.20235E+025, &
1.157E+025, &
1.0908E+025, &
1.03367E+025, &
9.83714E+024, &
9.39E+024, &
8.626E+024, &
7.985E+024, &
7.44571E+024, &
6.98E+024, &
6.222E+024, &
5.62917E+024, &
5.15357E+024, &
4.60118E+024, &
4.1825E+024, &
3.682E+024, &
3.32933E+024, &
3.07143E+024, &
2.857E+024, &
2.5678E+024, &
2.33E+024, &
2.18229E+024, &
2.06025E+024, &
1.8658E+024, &
1.705E+024, &
1.60179E+024, &
1.47988E+024, &
1.396E+024, &
1.274E+024, &
1.18923E+024, &
1.1072E+024, &
1.1146E+024, &
1.0898E+024, &
1.09736E+024, &
1.09038E+024, &
1.13938E+024, &
1.16058E+024, &
1.26347E+024, &
1.26347E+024/)

IF (Density.GT.0.) THEN
    ReducedElectricField=1.e21*ElectricField/Density ! E/n in Townsend as defined in LXCAT database
    mu = InterpolateCoefficient(muTable,ReducedElectricField)/Density ! table gives mu/n
    D = InterpolateCoefficient(EnergyTable,ReducedElectricField)*mu ! table gives energy=D/mu
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
