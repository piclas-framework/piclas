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

MODULE MOD_Radiation_ReadIn
!===================================================================================================================================
! Module for Radiation
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE Radiation_readin_atoms
  MODULE PROCEDURE Radiation_readin_atoms
END INTERFACE

INTERFACE Radiation_readin_molecules
  MODULE PROCEDURE Radiation_readin_molecules
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: Radiation_readin_atoms, Radiation_readin_molecules
!===================================================================================================================================

CONTAINS


SUBROUTINE Radiation_readin_atoms(iSpec)
!===================================================================================================================================
! Reads-in species constants for atomic radiation
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,    ONLY     : PlanckConst, c
USE MOD_Radiation_Vars,  ONLY     : SpeciesRadiation, RadiationInput
USE MOD_ReadInTools
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)              :: iSpec
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  CHARACTER(32)         :: hilf
  INTEGER :: errtemp, iLine, iLevel, EnLevelIndex, NumOfLevels
  REAL    :: dump  
!===================================================================================================================================
  WRITE(UNIT=hilf,FMT='(I0)') iSpec
  RadiationInput(iSpec)%RadiationSpectraFileName = GETSTR('Radiation-Species'//TRIM(hilf)//'-SpectraFileName','none')
  IF (RadiationInput(iSpec)%RadiationSpectraFileName.EQ.'none') THEN
    SpeciesRadiation(iSpec)%nLevels = 0
    SpeciesRadiation(iSpec)%nLines  = 0
    ! STOP
  END IF

  IF (RadiationInput(iSpec)%RadiationSpectraFileName.NE.'none') THEN
    OPEN( UNIT= 42, file = RadiationInput(iSpec)%RadiationSpectraFileName, status = 'old', form = 'formatted')
    errtemp=0
    READ(42,*,IOSTAT = errtemp) SpeciesRadiation(iSpec)%nLevels
  !  ALLOCATE(SpeciesRadiation(iSpec)%Level(SpeciesRadiation(iSpec)%nLevels,5))
    NumOfLevels= 0
    DO iLevel =1, SpeciesRadiation(iSpec)%nLevels
      READ(42,*,IOSTAT = errtemp) dump,dump, EnLevelIndex, dump, dump
      NumOfLevels = MAX(NumOfLevels, EnLevelIndex)
    END DO
    REWIND(42)
    READ(42,*,IOSTAT = errtemp) dump
    ALLOCATE(SpeciesRadiation(iSpec)%Level(NumOfLevels,5)) ! TODO change to max no of energy levels
    SpeciesRadiation(iSpec)%Level = 0.0
    ALLOCATE(SpeciesRadiation(iSpec)%NumDensExc(SpeciesRadiation(iSpec)%nLevels))
    DO iLevel =1, SpeciesRadiation(iSpec)%nLevels
      READ(42,*,IOSTAT = errtemp) SpeciesRadiation(iSpec)%Level(iLevel,1), SpeciesRadiation(iSpec)%Level(iLevel,2), &
        EnLevelIndex, SpeciesRadiation(iSpec)%Level(iLevel,4), SpeciesRadiation(iSpec)%Level(iLevel,5)
      SpeciesRadiation(iSpec)%Level(iLevel,2) = SpeciesRadiation(iSpec)%Level(iLevel,2)*PlanckConst*c*100.
      SpeciesRadiation(iSpec)%Level(EnLevelIndex,3) = iLevel
      IF (SpeciesRadiation(iSpec)%Level(iLevel,4) .EQ. 0.) SpeciesRadiation(iSpec)%Level(iLevel,4) = 3.
    END DO
    SWRITE(*,'(A6,I6,A17)') ' Found ',SpeciesRadiation(iSpec)%nLevels,' Level entries in File:'
    SWRITE(*,*) RadiationInput(iSpec)%RadiationSpectraFileName

    READ(42,*,IOSTAT = errtemp) SpeciesRadiation(iSpec)%nLines
    ALLOCATE(SpeciesRadiation(iSpec)%LinesReal(SpeciesRadiation(iSpec)%nLines,3))
    ALLOCATE(SpeciesRadiation(iSpec)%LinesInt(SpeciesRadiation(iSpec)%nLines,4))
    DO iLine =1, SpeciesRadiation(iSpec)%nLines
      READ(42,*,IOSTAT = errtemp) SpeciesRadiation(iSpec)%LinesReal(iLine,1), SpeciesRadiation(iSpec)%LinesInt(iLine,1),&
        SpeciesRadiation(iSpec)%LinesInt(iLine,2),&
        SpeciesRadiation(iSpec)%LinesInt(iLine,3),SpeciesRadiation(iSpec)%LinesInt(iLine,4), &
        SpeciesRadiation(iSpec)%LinesReal(iLine,2), & 
        SpeciesRadiation(iSpec)%LinesReal(iLine,3), dump
        SpeciesRadiation(iSpec)%LinesReal(iLine,1) = SpeciesRadiation(iSpec)%LinesReal(iLine,1)*1.E-10 !in m instead of Angstrom
    END DO
    SWRITE(*,'(A6,I6,A17)') ' Found ',SpeciesRadiation(iSpec)%nLines,' Line entries in File:'
    SWRITE(*,*) RadiationInput(iSpec)%RadiationSpectraFileName

    CLOSE(unit=42)
  END IF
END SUBROUTINE Radiation_readin_atoms





SUBROUTINE Radiation_readin_molecules(iSpec)
!===================================================================================================================================
! Reads-in species constants for molecular radiation
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,    ONLY     : PlanckConst, c
USE MOD_Radiation_Vars,  ONLY     : SpeciesRadiation, RadiationInput
USE MOD_ReadInTools
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)              :: iSpec
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  CHARACTER(66)         :: hilf
  INTEGER :: dump, length, errtemp, iLoop, i, charlen, numLines
  REAL    :: realdump
!===================================================================================================================================
  
  WRITE(UNIT=hilf,FMT='(I0)') iSpec
  RadiationInput(iSpec)%RadiationSpectraFileName = GETSTR('Radiation-Species'//TRIM(hilf)//'-SpectraFileName','none')
  IF (RadiationInput(iSpec)%RadiationSpectraFileName.EQ.'none') THEN
    SpeciesRadiation(iSpec)%nBands = 0
    ! STOP
  END IF

  IF (RadiationInput(iSpec)%RadiationSpectraFileName.NE.'none') THEN
    OPEN( UNIT= 304, file = RadiationInput(iSpec)%RadiationSpectraFileName, status = 'old', form = 'formatted')
    errtemp=0  
    
    SpeciesRadiation(iSpec)%nBands = 0
    DO 
      READ(304,*,IOSTAT = errtemp) hilf
      IF (errtemp.NE.0) EXIT
      IF (hilf(1:5).EQ.'[band') SpeciesRadiation(iSpec)%nBands = SpeciesRadiation(iSpec)%nBands + 1
    END DO
    REWIND(304)
    ALLOCATE(SpeciesRadiation(iSpec)%BandName(SpeciesRadiation(iSpec)%nBands))
    ALLOCATE(SpeciesRadiation(iSpec)%NumMolecularTransitions(SpeciesRadiation(iSpec)%nBands))
    
    READ(304,*,IOSTAT = errtemp) hilf
    DO WHILE (hilf(1:7) .NE. '[Energy') 
      READ(304,*,IOSTAT = errtemp) hilf
    END DO
    READ(304,*,IOSTAT = errtemp) length, dump, dump, dump, dump
    ALLOCATE(SpeciesRadiation(iSpec)%EnergyLevelName(length))
    DO iLoop=1, length
      READ(304,*,IOSTAT = errtemp) hilf 
      DO WHILE (hilf(1:1).EQ.'c')
        READ(304,*,IOSTAT = errtemp) hilf    
      END DO
      SpeciesRadiation(iSpec)%EnergyLevelName(iLoop) = TRIM(hilf)
    END DO

    ALLOCATE(SpeciesRadiation(iSpec)%BandProperties(SpeciesRadiation(iSpec)%nBands,4))
    DO iLoop=1, SpeciesRadiation(iSpec)%nBands
      DO WHILE (hilf(1:5).NE.'[band')
        READ(304,*,IOSTAT = errtemp) hilf 
      END DO
      charlen = LEN_TRIM(hilf)
      hilf(1:charlen-7) = hilf(7:charlen-1)
      hilf(charlen-6:charlen) = ' '
      hilf=TRIM(hilf)
      SpeciesRadiation(iSpec)%BandName(iLoop) = hilf
      
      READ(304,*,IOSTAT = errtemp) hilf 
      
      READ(304,*,IOSTAT = errtemp) hilf 
      DO WHILE (hilf(1:1).EQ.'c')
        READ(304,*,IOSTAT = errtemp) hilf    
      END DO
      BACKSPACE(304)
      READ(304,*,IOSTAT = errtemp) SpeciesRadiation(iSpec)%BandProperties(iLoop,1), & 
          SpeciesRadiation(iSpec)%BandProperties(iLoop,2), & 
          SpeciesRadiation(iSpec)%BandProperties(iLoop,3), realdump
      SpeciesRadiation(iSpec)%BandProperties(iLoop,4) = NINT(realdump)
    END DO

    ALLOCATE(SpeciesRadiation(iSpec)%EnergyLevelProperties(length,17))
    DO iLoop=1, length
      REWIND(304)
      DO WHILE (hilf .NE. '['//TRIM(SpeciesRadiation(iSpec)%EnergyLevelName(iLoop))//']')
        READ(304,*,IOSTAT = errtemp) hilf
      IF (hilf(1:1).EQ.'c') CYCLE
      END DO
      READ(304,*,IOSTAT = errtemp) hilf
      i=0
      DO WHILE (i.NE.17)
        READ(304,*,IOSTAT = errtemp) hilf 
        IF (hilf(1:1).EQ.'c') CYCLE
        BACKSPACE(304)
        i = i + 1
        READ(304,*,IOSTAT = errtemp) SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,i)      
      END DO
    END DO

    DO iLoop=1, length
      SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,15)=SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,15)*0.01
      SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,2)=SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,2)*PlanckConst*c*100.
      SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,3)=SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,3)*PlanckConst*c*100.
      SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,4)=SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,4)*PlanckConst*c*100.
      SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,5)=SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,5)*PlanckConst*c*100.
      SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,6)=SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,6)*PlanckConst*c*100.
      SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,7)=SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,7)*PlanckConst*c*100.
      SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,8)=SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,8)*PlanckConst*c*100.
      SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,9)=SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,9)*PlanckConst*c*100.
      SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,13)=SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,13) &
        * PlanckConst*c*100.
      SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,14)=SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,14) & 
        * PlanckConst*c*100.
      
      SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,16)=SpeciesRadiation(iSpec)%EnergyLevelProperties(iLoop,16) & 
        * PlanckConst*c*100.
    END DO

    ALLOCATE(SpeciesRadiation(iSpec)%Bands(SpeciesRadiation(iSpec)%nBands))  
    DO iLoop=1, SpeciesRadiation(iSpec)%nBands
      REWIND(304)
      DO WHILE (hilf .NE. '[transition-'//TRIM(SpeciesRadiation(iSpec)%BandName(iLoop))//']')
        READ(304,*,IOSTAT = errtemp) hilf
      END DO
      READ(304,*,IOSTAT = errtemp) numLines, dump, dump, dump, dump
      SpeciesRadiation(iSpec)%NumMolecularTransitions(iLoop) = numLines
      ALLOCATE(SpeciesRadiation(iSpec)%Bands(iLoop)%MolTransLines(numLines,5))
      READ(304,*,IOSTAT = errtemp) hilf
      DO i = 1, numLines
        READ(304,*,IOSTAT = errtemp) SpeciesRadiation(iSpec)%Bands(iLoop)%MolTransLines(i,1), &
          SpeciesRadiation(iSpec)%Bands(iLoop)%MolTransLines(i,2), SpeciesRadiation(iSpec)%Bands(iLoop)%MolTransLines(i,3), &
          SpeciesRadiation(iSpec)%Bands(iLoop)%MolTransLines(i,4)
  !      PRINT*, SpeciesRadiation(iSpec)%Bands(iLoop)%MolTransLines(i,1), &
  !        SpeciesRadiation(iSpec)%Bands(iLoop)%MolTransLines(i,2), SpeciesRadiation(iSpec)%Bands(iLoop)%MolTransLines(i,3), &
  !        SpeciesRadiation(iSpec)%Bands(iLoop)%MolTransLines(i,4)!, SpeciesRadiation(iSpec)%Bands(iLoop)%MolTransLines(i,5)
      END DO 
    END DO

    SpeciesRadiation(iSpec)%nLevels = length ! TODO: check!

    SWRITE(*,'(A6,I6,A17)') ' Found ',length,' energy levels in File:'
    SWRITE(*,*) TRIM(RadiationInput(iSpec)%RadiationSpectraFileName)
    DO iLoop = 1,length
      SWRITE(*,*) TRIM(SpeciesRadiation(iSpec)%EnergyLevelName(iLoop))    
    END DO

    SWRITE(*,'(A6,I6,A17)') ' Found ',SpeciesRadiation(iSpec)%nBands,' transition bands in File:'
    SWRITE(*,*) TRIM(RadiationInput(iSpec)%RadiationSpectraFileName)
    DO iLoop = 1,SpeciesRadiation(iSpec)%nBands
      SWRITE(*,*) TRIM(SpeciesRadiation(iSpec)%BandName(iLoop)), ' (', &
        SpeciesRadiation(iSpec)%NumMolecularTransitions(iLoop), 'lines)'    
    END DO

    CLOSE(unit=304)
  END IF
END SUBROUTINE Radiation_readin_molecules


END MODULE MOD_Radiation_ReadIn
