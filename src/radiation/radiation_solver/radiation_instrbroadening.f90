!==================================================================================================================================
! Copyright (c) 2010 - 2019 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_Radiation_InstrBroadening
!===================================================================================================================================
! Module for Radiation
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE radiation_instrbroadening
  MODULE PROCEDURE radiation_instrbroadening
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: radiation_instrbroadening
!===================================================================================================================================

CONTAINS


SUBROUTINE radiation_instrbroadening(iElem)
!===================================================================================================================================
! Accounting for trapezoidal slit-function (instrumental broadening)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Radiation_Vars,       ONLY   : RadiationParameter, Radiation_Emission_spec
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)             :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

  REAL, PARAMETER   :: basewidth = 3.42E-10*5.
  REAL, PARAMETER   :: topwidth  = 1.79E-10*5.
  REAL              :: basewidth_half, topwidth_half, slope
  REAL              :: wavelength_min_base, wavelength_max_base, wavelength_min_top, wavelength_max_top
  INTEGER           :: io_error, w
  INTEGER           :: iWave, iWave_min, i
  REAL, ALLOCATABLE :: Radiation_Emission_spec_conv(:)
  REAL              :: fractionl, fractionr, delta_base, delta_top

!===================================================================================================================================
  
  ALLOCATE(Radiation_Emission_spec_conv(RadiationParameter%WaveLenDiscr))

  iWave_min = 1!0

  basewidth_half = 0.5 * basewidth
  topwidth_half  = 0.5 * topwidth
  slope          = 1. / (basewidth_half-topwidth_half)

  DO iWave=1, RadiationParameter%WaveLenDiscr
    wavelength_min_base = RadiationParameter%WaveLen(iWave) - basewidth_half
    wavelength_max_base = RadiationParameter%WaveLen(iWave) + basewidth_half
    wavelength_min_top  = RadiationParameter%WaveLen(iWave) - topwidth_half
    wavelength_max_top  = RadiationParameter%WaveLen(iWave) + topwidth_half

! --- start index determination
    DO WHILE(RadiationParameter%WaveLen(iWave_min+1) .LT. wavelength_min_base)
      iWave_min = iWave_min + 1
    END DO

! --- slit function
    DO i = iWave_min, RadiationParameter%WaveLenDiscr-1
      IF(RadiationParameter%WaveLen(i) .LT. wavelength_min_base) THEN
        fractionl = 0.
        IF(RadiationParameter%WaveLen(i+1) .GT. wavelength_min_top) THEN
          STOP 'slit function: step width is too big!'
        END IF
        fractionr  = slope * (RadiationParameter%WaveLen(i+1) - RadiationParameter%WaveLen(iWave) + basewidth_half)
        delta_base = RadiationParameter%WaveLen(i+1) - wavelength_min_base
        delta_top  = 0.
      ELSEIF(RadiationParameter%WaveLen(i+1) .LT. wavelength_min_top) THEN
        fractionl  = slope * (RadiationParameter%WaveLen(i  ) - RadiationParameter%WaveLen(iWave) + basewidth_half)
        fractionr  = slope * (RadiationParameter%WaveLen(i+1) - RadiationParameter%WaveLen(iWave) + basewidth_half)
        delta_base = RadiationParameter%WaveLenIncr
        delta_top  = 0.
      ELSEIF(RadiationParameter%WaveLen(i  ) .LT. wavelength_min_top) THEN
        fractionl  = slope * (RadiationParameter%WaveLen(i  ) - RadiationParameter%WaveLen(iWave) + basewidth_half)
        fractionr  = 1.
        delta_base = wavelength_min_top - RadiationParameter%WaveLen(i)
        delta_top  = RadiationParameter%WaveLen(i+1) - wavelength_min_top
      ELSEIF(RadiationParameter%WaveLen(i+1) .LT. wavelength_max_top) THEN
        delta_base = 0.
        delta_top  = RadiationParameter%WaveLenIncr
      ELSEIF(RadiationParameter%WaveLen(i  ) .LT. wavelength_max_top) THEN
        fractionl  = 1.
        fractionr  = - slope * (RadiationParameter%WaveLen(i+1) - RadiationParameter%WaveLen(iWave) - basewidth_half)
        delta_base = RadiationParameter%WaveLen(i+1) - wavelength_max_top
        delta_top  = wavelength_max_top - RadiationParameter%WaveLen(i)
      ELSEIF(RadiationParameter%WaveLen(i+1) .LT. wavelength_max_base) THEN
        fractionl  = - slope * (RadiationParameter%WaveLen(i  ) - RadiationParameter%WaveLen(iWave) - basewidth_half)
        fractionr  = - slope * (RadiationParameter%WaveLen(i+1) - RadiationParameter%WaveLen(iWave) - basewidth_half)
        delta_base = RadiationParameter%WaveLenIncr
        delta_top  = 0.
      ELSEIF(RadiationParameter%WaveLen(i) .LT. wavelength_max_base) THEN
        fractionl  = - slope * (RadiationParameter%WaveLen(i  ) - RadiationParameter%WaveLen(iWave) - basewidth_half)
        fractionr  = 0.
        delta_base = wavelength_max_base - RadiationParameter%WaveLen(i)
        delta_top  = 0.
      ELSE
        EXIT
      END IF

      Radiation_Emission_spec_conv(iWave) = Radiation_Emission_spec_conv(iWave) &
          + ((fractionl+fractionr)*.5*delta_base+delta_top) &
          * Radiation_Emission_spec(i+1,iElem)


    END DO

! --- transformation to Laux's units
!    Radiation_Emission_spec_conv(iWave) = Radiation_Emission_spec_conv(iWave)/1.D11

  END DO

  OPEN(unit=30,file='Broadening.dat',status='replace',action='write', iostat=io_error)
  DO w=1, RadiationParameter%WaveLenDiscr
    WRITE(30,*) RadiationParameter%WaveLen(w)*1.E10 , Radiation_Emission_spec_conv(w)!*1.E10!/3.E12
  END DO
  CLOSE(unit=30)

  WRITE(*,*) '***  INSTRUMENTAL BROADENING SUCCESSFULLY DONE  ***'
  WRITE(*,*) ''

END SUBROUTINE radiation_instrbroadening


END MODULE MOD_Radiation_InstrBroadening
