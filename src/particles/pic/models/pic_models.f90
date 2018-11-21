!==================================================================================================================================
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

MODULE MOD_PICModels
!===================================================================================================================================
! Physical models for PIC
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE FieldIonization
  MODULE PROCEDURE FieldIonization
END INTERFACE
PUBLIC::FieldIonization
!===================================================================================================================================

CONTAINS

SUBROUTINE FieldIonization()
!===================================================================================================================================
! Field Ionization:
! * Ammosov-Delone-Krainov (ADK) model (tunnel ionization)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst, ElementaryCharge
USE MOD_TimeDisc_Vars         ,ONLY: dt
USE MOD_Particle_Vars         ,ONLY: PDM, Species, PartSpecies,  usevMPF, PartState, PEM, PartMPF
USE MOD_DSMC_Vars             ,ONLY: DSMC, SpecDSMC
USE MOD_PICInterpolation_Vars ,ONLY: FieldAtParticle
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iPart, MaxElecQua, ChargedNum, SumOfFormedParticles, ElectronIndex
REAL                    :: FieldStrength_GV, IonizationEnergy_eV, iRan, QuantumTunnelProb, EffQuantNum
REAL                    :: CriticalValue_GV
!===================================================================================================================================
SumOfFormedParticles = 0

DO iPart = 1, PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart)) THEN
    IF(SpecDSMC(PartSpecies(iPart))%NextIonizationSpecies.NE.0) THEN
      ! Saving field strength in giga volt
      FieldStrength_GV = SQRT(FieldAtParticle(iPart,1)**2 + FieldAtParticle(iPart,2)**2 + FieldAtParticle(iPart,3)**2) / 1E9
      ! Ionization energy (same as in QK model)
      MaxElecQua=SpecDSMC(PartSpecies(iPart))%MaxElecQuant - 1
      IonizationEnergy_eV=SpecDSMC(PartSpecies(iPart))%ElectronicState(2,MaxElecQua)*BoltzmannConst / ElementaryCharge
      ! Checking applicability of ADK model (normalized variables)
      ! Particle-in-cell simulations of tunneling ionization effects in plasma-based accelerators, David L. Bruhwiler
      CriticalValue_GV = (SQRT(2.) - 1.) * (IonizationEnergy_eV / 27.2)**(1.5) * 5.14E+2
      IF(FieldStrength_GV.GT.CriticalValue_GV) THEN
        CALL abort(&
          __STAMP__&
          ,'ERROR FieldIonization: ADK model is not applicable for electric fields > critical value!', PartSpecies(iPart))
      END IF
      ! Z (ChargedNum): Charge number of the atom/ion AFTER the ionization (thus + 1)
      ChargedNum = NINT(Species(PartSpecies(iPart))%ChargeIC/ElementaryCharge) + 1
      EffQuantNum = 3.69*REAL(ChargedNum) / SQRT(IonizationEnergy_eV)
      QuantumTunnelProb = 1.52E+15 * 4.**(EffQuantNum)*IonizationEnergy_eV / (EffQuantNum*GAMMA(2.*EffQuantNum)) &
                        * (20.5*IonizationEnergy_eV**(3./2.)/FieldStrength_GV)**(2.*(EffQuantNum-1.)) &
                        * EXP(-6.83*IonizationEnergy_eV**(3./2.)/FieldStrength_GV) * dt
      CALL RANDOM_NUMBER(iRan)
      IF(QuantumTunnelProb.GT.iRan) THEN
        !.... Get free particle index for the 3rd particle produced
        SumOfFormedParticles = SumOfFormedParticles + 1
        ElectronIndex = PDM%nextFreePosition(SumOfFormedParticles+PDM%CurrentNextFreePosition)
        IF (ElectronIndex.EQ.0) THEN
          CALL abort(__STAMP__,&
          'New Particle Number greater max Part Num in Field Ionization.')
        END IF
        !Set new Species of new particle
        PDM%ParticleInside(ElectronIndex) = .TRUE.
        PartSpecies(ElectronIndex) = DSMC%ElectronSpecies
        PartState(ElectronIndex,1:3) = PartState(iPart,1:3)
        PartState(ElectronIndex,4:6) = Species(DSMC%ElectronSpecies)%MassIC / Species(PartSpecies(iPart))%MassIC &
                                        * PartState(iPart,4:6)
        PartState(iPart,4:6) = Species(SpecDSMC(PartSpecies(iPart))%NextIonizationSpecies)%MassIC &
                                / Species(PartSpecies(iPart))%MassIC * PartState(iPart,4:6)
        PEM%Element(ElectronIndex) = PEM%Element(iPart)
        ! Setting the species of the ionized particle
        PartSpecies(iPart) = SpecDSMC(PartSpecies(iPart))%NextIonizationSpecies
        IF(usevMPF) PartMPF(ElectronIndex) = PartMPF(iPart)
        ! Setting the field for the new particle for the following integration
        FieldAtParticle(ElectronIndex,1:6) = FieldAtParticle(iPart,1:6)
      END IF
    END IF
  END IF
END DO

PDM%ParticleVecLength = PDM%ParticleVecLength + SumOfFormedParticles
PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + SumOfFormedParticles

END SUBROUTINE FieldIonization


END MODULE MOD_PICModels
