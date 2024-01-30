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

MODULE MOD_Radiation_Excitation
!===================================================================================================================================
! Module for calculation of the excited state density for radiative transitions
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE radiation_excitation
  MODULE PROCEDURE radiation_excitation
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: radiation_excitation
!===================================================================================================================================

CONTAINS


SUBROUTINE radiation_excitation()
!===================================================================================================================================
! Main routine of populating the excited state
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,      ONLY   : ElementaryCharge, BoltzmannConst
USE MOD_Radiation_Vars,    ONLY   : RadiationInput, SpeciesRadiation, NumDensElectrons
USE MOD_PARTICLE_Vars,     ONLY   : nSpecies
USE MOD_DSMC_Vars,         ONLY   : SpecDSMC

! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

  REAL              :: low_IonizationPot                     !lowering of ionization Potenzial [J]
!  REAL              :: ElectronicPartFunc                    !Partition function for  Maxwell distribution of excited state
  INTEGER           :: iLevel, iLevel_considered             ! loop indices
  INTEGER           :: nLevels_considered                    ! actual number of considered levels           
  INTEGER           :: iSpec
  REAL              :: BoltzmannFactor                       ! degeneracy*EXP(-E/kB T)
  REAL              :: RotVibPartFunc                        ! combined rotational-vibrational partition function
  REAL              :: Gvib, Gvib_prev                       ! normalized vibrational energy (G=Delta EVib / hc)
  INTEGER           :: v                                     ! vibrational quantum number [-]
  REAL              :: deltaqv
!===================================================================================================================================

  DO iSpec = 1, nSpecies

    IF(.NOT.RadiationInput(iSpec)%DoRadiation) CYCLE

! --- atoms (1) and atomic ions (10)
    IF((SpecDSMC(iSpec)%InterID .EQ. 1) .OR. (SpecDSMC(iSpec)%InterID .EQ. 10)) THEN
      IF (RadiationInput(iSpec)%Telec.LE.0.0) CYCLE
      IF (SpeciesRadiation(iSpec)%nLevels.EQ.0) CYCLE
      low_IonizationPot = 2.9E-8*SQRT(NumDensElectrons/1.E6/MAX(1.,RadiationInput(iSpec)%Telec))*ElementaryCharge

      nLevels_considered  = SpeciesRadiation(iSpec)%nLevels

      IF (low_IonizationPot .NE. 0.0) THEN
        DO iLevel_considered = 1, SpeciesRadiation(iSpec)%nLevels
          IF ( SpeciesRadiation(iSpec)%Level(iLevel_considered,2) .LT. (RadiationInput(iSpec)%IonizationEn-low_IonizationPot) ) THEN
            nLevels_considered = iLevel_considered
          END IF
        END DO
      END IF

      SpeciesRadiation(iSpec)%PartFunc = 0.0

      DO iLevel = 1, nLevels_considered !SpeciesRadiation(1)%nLevels
        SpeciesRadiation(iSpec)%PartFunc = SpeciesRadiation(iSpec)%PartFunc + SpeciesRadiation(iSpec)%Level(iLevel,1) &
          * EXP(-SpeciesRadiation(iSpec)%Level(iLevel,2)/(BoltzmannConst*RadiationInput(iSpec)%Telec))
      END DO

      DO iLevel = 1, nLevels_considered !SpeciesRadiation(1)%nLevels
        SpeciesRadiation(iSpec)%NumDensExc(iLevel) = RadiationInput(iSpec)%NumDens*SpeciesRadiation(iSpec)%Level(iLevel,1) &
          * EXP(-SpeciesRadiation(iSpec)%Level(iLevel,2)/(BoltzmannConst*RadiationInput(iSpec)%Telec))&
          /SpeciesRadiation(iSpec)%PartFunc
      END DO 
   
! --- diatomic molecules (2) and diatomic molecular ions (20)
    ELSEIF((SpecDSMC(iSpec)%InterID .EQ. 2) .OR. (SpecDSMC(iSpec)%InterID .EQ. 20)) THEN 

!! --- Initialization
      SpeciesRadiation(iSpec)%PartFunc = 0.0
      IF ((RadiationInput(iSpec)%Telec.LE.0.0).OR.(RadiationInput(iSpec)%Tvib.LE.0.0))CYCLE !!!!!!!TODO!!!
      DO iLevel = 1, SpeciesRadiation(iSpec)%nLevels
! --- Initialization
        Gvib_prev = 0.0
        v         = 0
        RotVibPartFunc = 0.0

! --- calculation of Boltzmann Factor (ge*EXP(-E/kBT))
        BoltzmannFactor = SpeciesRadiation(iSpec)%EnergyLevelProperties(iLevel,1) & 
          * EXP(MIN(7.d2, -SpeciesRadiation(iSpec)%EnergyLevelProperties(iLevel,2) / (BoltzmannConst*RadiationInput(iSpec)%Telec)))

        DO !WHILE (.TRUE.)
! --- vibrational energy of quantum number v, powers of the vibrational quantum number + .5
! --- G = omega_w(v+1/2) - omega_w x_e(v+1/2)**2 + omega_w y_e(v+1/2)**3 + omega_w z_e(v+1/2)**4 + ...
          Gvib = SpeciesRadiation(iSpec)%EnergyLevelProperties(iLevel,4) * (REAL(v)+0.5) &
            - SpeciesRadiation(iSpec)%EnergyLevelProperties(iLevel,5) * (REAL(v)+0.5)**2 &
            + SpeciesRadiation(iSpec)%EnergyLevelProperties(iLevel,6) * (REAL(v)+0.5)**3 &
            + SpeciesRadiation(iSpec)%EnergyLevelProperties(iLevel,7) * (REAL(v)+0.5)**4

! --- contribution of vibrational quantum number v to the partition function
          deltaqv = EXP(-Gvib/(BoltzmannConst*RadiationInput(iSpec)%Tvib))

! --- cutoff criteria
          IF(Gvib .GT. SpeciesRadiation(iSpec)%EnergyLevelProperties(iLevel,3)) THEN
            EXIT !Vibrational energy exceeds dissociation energy
          ELSEIF(Gvib .LT. Gvib_prev) THEN
            EXIT !Vibrational energy reached fictitious peak?
          END IF

! --- contribution of rotational excitation
          IF((SpeciesRadiation(iSpec)%EnergyLevelProperties(iLevel,8) &
            - SpeciesRadiation(iSpec)%EnergyLevelProperties(iLevel,9)*(REAL(v)+0.5) ) .GT. 0.0) THEN

! --- Analytic expression with cutoff at dissociation energy
            RotVibPartFunc = RotVibPartFunc + RadiationInput(iSpec)%Trot * deltaqv &
              / ( 1.0 / BoltzmannConst * ( SpeciesRadiation(iSpec)%EnergyLevelProperties(iLevel,8) &
              - SpeciesRadiation(iSpec)%EnergyLevelProperties(iLevel,9)*(REAL(v)+0.5) ) ) & ! + gamma_e*(REAL(v)+0.5)**2
              * ( 1.d0 - exp( ( Gvib - SpeciesRadiation(iSpec)%EnergyLevelProperties(iLevel,3) ) &
              / BoltzmannConst / RadiationInput(iSpec)%Trot) )

          END IF

          Gvib_prev = Gvib
          v=v+1

        END DO

        SpeciesRadiation(iSpec)%PartFunc = SpeciesRadiation(iSpec)%PartFunc + BoltzmannFactor * RotVibPartFunc
        
      END DO

    ELSE

      CYCLE

    END IF

  END DO

END SUBROUTINE radiation_excitation

END MODULE MOD_Radiation_Excitation
