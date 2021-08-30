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

MODULE MOD_Radiation_ExportSpectrum
!===================================================================================================================================
! Module for Radiation
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE radiation_exportspectrum
  MODULE PROCEDURE radiation_exportspectrum
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: radiation_exportspectrum
!===================================================================================================================================

CONTAINS


SUBROUTINE radiation_exportspectrum(iElement, output_format)
!===================================================================================================================================
! exports spectral emission and absorption to *.dat-files
! 
! Radiation_Emission_Absorption_xx_iGlobalElement(:,3) -> wavelength, emission (W/m3/str/m), absorption (1/m))
! 
! further formats can be added to output_format, currently available:
! 1: gnuplot
! 2: tikz
!===================================================================================================================================
! MODULES
  USE MOD_Radiation_Vars,       ONLY: RadiationParameter, Radiation_Emission_spec, Radiation_Absorption_spec
  USE MOD_Mesh_Vars,            ONLY: offsetElem
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)             :: iElement, output_format
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER           :: w, io_error, iGlobalElement
  CHARACTER(32)     :: hilf
!===================================================================================================================================
  
  iGlobalElement = iElement + offsetElem
  WRITE(UNIT=hilf,FMT='(I0)') iGlobalElement

  SELECT CASE(output_format)

    CASE(1) !Output for gnuplot
      OPEN(unit=20,file='Radiation_Emission_Absorption_'//TRIM(hilf)//'.dat', status='replace',action='write',iostat=io_error)
      DO w=1, RadiationParameter%WaveLenDiscr
        WRITE(20,*) RadiationParameter%WaveLen(w)*1.E10, Radiation_Emission_spec(w,iElement), Radiation_Absorption_spec(w,iElement)
      END DO
      CLOSE(unit=20)

    CASE(2)!Output for pgfplots
      OPEN(unit=20,file='Radiation_Emission_Absorption_tikz_'//TRIM(hilf)//'.dat', status='replace',action='write',iostat=io_error)
      WRITE(20,*) 'x,y1,y2'
      DO w=1, RadiationParameter%WaveLenDiscr
        WRITE(20,*) RadiationParameter%WaveLen(w)*1.E10,',',Radiation_Emission_spec(w,iElement),',', &
                    Radiation_Absorption_spec(w,iElement)
      END DO
      CLOSE(unit=20)

    CASE DEFAULT
      WRITE(*,*) 'output format is not defined'

  END SELECT

END SUBROUTINE radiation_exportspectrum


END MODULE MOD_Radiation_ExportSpectrum
