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
! * Ammosov-Delone-Krainov (ADK) model (only tunnel ionization no BSI)
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars,ONLY:FieldIonizationModel
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SELECT CASE(FieldIonizationModel)
CASE(1)
  CALL ADK_Bruhwiler2003() ! Bruhwiler 2003: requires E<E_crit (without BSI)
CASE(2)
  CALL ADK_Yu2018() ! Yu 2018: used for tunneling/BSI regardless of E_crit
END SELECT
END SUBROUTINE FieldIonization


SUBROUTINE ADK_Bruhwiler2003()
!===================================================================================================================================
! Field Ionization:
! * Ammosov-Delone-Krainov (ADK) model (only tunnel ionization no BSI)
! * from Bruhwiler, Particle-in-cell simulations of tunneling ionization effects in plasma-based accelerators, 2003
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
REAL              :: FieldStrength
!===================================================================================================================================
SumOfFormedParticles = 0

DO iPart = 1, PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart)) THEN
    ASSOCIATE (& 
          oldSpec => PartSpecies(iPart) ,&
          newSpec => SpecDSMC(PartSpecies(iPart))%NextIonizationSpecies )
      IF(newSpec.EQ.0) CYCLE

      ! Saving field strength in giga volt
      FieldStrength_GV = SQRT(FieldAtParticle(iPart,1)**2 + FieldAtParticle(iPart,2)**2 + FieldAtParticle(iPart,3)**2) / 1E9

      !write(*,*) "E="
      !read(*,*)   FieldStrength
      !FieldStrength_GV=FieldStrength/1e9
      
      ! Ionization energy (same as in QK model)
      MaxElecQua=SpecDSMC(oldSpec)%MaxElecQuant - 1
      IonizationEnergy_eV=SpecDSMC(oldSpec)%ElectronicState(2,MaxElecQua)*BoltzmannConst / ElementaryCharge
      ! Checking applicability of ADK model (normalized variables)
      ! Particle-in-cell simulations of tunneling ionization effects in plasma-based accelerators, David L. Bruhwiler
      CriticalValue_GV = (SQRT(2.) - 1.) * (IonizationEnergy_eV / 27.2)**(1.5) * 5.14E+2
      IF(FieldStrength_GV.GT.CriticalValue_GV) THEN
        WRITE (*,*) "IonizationEnergy_eV =", IonizationEnergy_eV
        WRITE (*,*) "CriticalValue_GV    =", CriticalValue_GV   
        WRITE (*,*) "FieldStrength_GV    =", FieldStrength_GV   
        CALL abort(&
            __STAMP__&
            ,'ERROR FieldIonization: ADK model is not applicable for electric fields > critical value!', oldSpec)
      END IF
      ! Z (ChargedNum): Charge number of the atom/ion AFTER the ionization (thus + 1)
      ChargedNum = NINT(Species(oldSpec)%ChargeIC/ElementaryCharge) + 1
      EffQuantNum = 3.69*REAL(ChargedNum) / SQRT(IonizationEnergy_eV)
      QuantumTunnelProb = 1.52E+15 * 4.**(EffQuantNum)*IonizationEnergy_eV / (EffQuantNum*GAMMA(2.*EffQuantNum)) &
          * (20.5*IonizationEnergy_eV**(3./2.)/FieldStrength_GV)**(2.*(EffQuantNum-1.)) &
          * EXP(-6.83*IonizationEnergy_eV**(3./2.)/FieldStrength_GV) * dt
      !WRITE (*,'(ES25.14E3,ES25.14E3)') FieldStrength,QuantumTunnelProb/dt
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
        PartSpecies(ElectronIndex)        = DSMC%ElectronSpecies
        PartState(ElectronIndex,1:3)      = PartState(iPart,1:3)
        PartState(ElectronIndex,4:6)      = Species(DSMC%ElectronSpecies)%MassIC / Species(oldSpec)%MassIC * PartState(iPart,4:6)
        PartState(iPart,4:6)              = Species(newSpec)%MassIC              / Species(oldSpec)%MassIC * PartState(iPart,4:6)
        PEM%Element(ElectronIndex)        = PEM%Element(iPart)
        ! Setting the species of the ionized particle
        oldSpec = newSpec
        IF(usevMPF) PartMPF(ElectronIndex) = PartMPF(iPart)
        ! Setting the field for the new particle for the following integration
        FieldAtParticle(ElectronIndex,1:6) = FieldAtParticle(iPart,1:6)
      END IF
    END ASSOCIATE
  END IF
END DO

PDM%ParticleVecLength = PDM%ParticleVecLength + SumOfFormedParticles
PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + SumOfFormedParticles

END SUBROUTINE ADK_Bruhwiler2003


SUBROUTINE ADK_Yu2018()
!===================================================================================================================================
! Field Ionization:
! * Ammosov-Delone-Krainov (ADK) model
! * from Yu, Shaping of ion energy spectrum due to ionization in ion acceleration driven by an ultra-short pulse laser, 2018
!   (which is originally from Penetrante, Residual energy in plasmas produced by intense subpicosecond lasers, 1991)
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
REAL                    :: FieldStrength, IonizationEnergy_eV, iRan, QuantumTunnelProb, EffQuantNum
REAL              :: n
REAL                    :: CriticalValue
!===================================================================================================================================
SumOfFormedParticles = 0

DO iPart = 1, PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart)) THEN
    ASSOCIATE ( oldSpec => PartSpecies(iPart) ,&
                newSpec => SpecDSMC(PartSpecies(iPart))%NextIonizationSpecies )
      IF(newSpec.EQ.0) CYCLE

      !write(*,*) "E="
      !read(*,*)   FieldStrength
      !WRITE (*,*) "FieldStrength =", FieldStrength
      ASSOCIATE (&
            E_au     => 5.1e11 ,& ! [V/m] atomic unit field strength
            omega_au => 4.1e16 ,& ! [1/s] atomic unit frequency strength
            Z        => NINT(Species(oldSpec)%ChargeIC/ElementaryCharge) + 1 ,& ! Charge number of the ion AFTER the ionization (+1)
            !E        => FieldStrength&
            E        => SQRT(FieldAtParticle(iPart,1)**2 + FieldAtParticle(iPart,2)**2 + FieldAtParticle(iPart,3)**2)& ! [V/m]
            )
        ! Ionization energy (same as in QK model)
        MaxElecQua=SpecDSMC(oldSpec)%MaxElecQuant - 1
        IonizationEnergy_eV=SpecDSMC(oldSpec)%ElectronicState(2,MaxElecQua)*BoltzmannConst / ElementaryCharge
        ! Checking applicability of ADK model (normalized variables)
        ! Particle-in-cell simulations of tunneling ionization effects in plasma-based accelerators, David L. Bruhwiler
        CriticalValue = (SQRT(2.) - 1.) * (IonizationEnergy_eV / 27.2)**(1.5) * 5.14E+2 * 1e9
        !IF(E.GT.CriticalValue) THEN
          !WRITE (*,*) "IonizationEnergy_eV =", IonizationEnergy_eV
          !WRITE (*,*) "CriticalValue    =", CriticalValue   
          !WRITE (*,*) "FieldStrength    =", E   
          !CALL abort(&
          !__STAMP__&
          !,'ERROR FieldIonization: ADK model is not applicable for electric fields > critical value!', oldSpec)
        !END IF
        n = 3.69*REAL(Z) / SQRT(IonizationEnergy_eV)
        QuantumTunnelProb = 1.61 * omega_au * Z**2 * n**(-9./2.)  &
            * ((10.87 * Z**3 * E_au / (n**4 * E))**(2*n-3./2.))   &
            * EXP(-2*Z**3*E_au / (3*n**3*E))                      &
            * dt
      END ASSOCIATE
      !WRITE (*,*) TRIM(SpecDSMC(oldSpec)%Name)," ==> ",TRIM(SpecDSMC(newSpec)%Name)," + e^-"
      !WRITE (*,WRITEFORMAT) QuantumTunnelProb/dt
      !WRITE (*,'(ES25.14E3,ES25.14E3)') FieldStrength,QuantumTunnelProb/dt
      CALL RANDOM_NUMBER(iRan)
      !WRITE (*,*) "QuantumTunnelProb =", QuantumTunnelProb
      !WRITE (*,*) " "  
      IF(QuantumTunnelProb.GT.iRan) THEN
        !WRITE (*,*) TRIM(SpecDSMC(oldSpec)%Name)," ==> ",TRIM(SpecDSMC(newSpec)%Name)," + e^-"
        !.... Get free particle index for the 3rd particle produced
        SumOfFormedParticles = SumOfFormedParticles + 1
        ElectronIndex = PDM%nextFreePosition(SumOfFormedParticles+PDM%CurrentNextFreePosition)
        IF (ElectronIndex.EQ.0) THEN
          CALL abort(__STAMP__,&
              'New Particle Number greater max Part Num in Field Ionization.')
        END IF
        !Set new Species of new particle
        PDM%ParticleInside(ElectronIndex) = .TRUE.
        PartSpecies(ElectronIndex)        = DSMC%ElectronSpecies
        PartState(ElectronIndex,1:3)      = PartState(iPart,1:3)
        PartState(ElectronIndex,4:6)      = Species(DSMC%ElectronSpecies)%MassIC / Species(oldSpec)%MassIC * PartState(iPart,4:6)
        PartState(iPart,4:6)              = Species(newSpec)%MassIC              / Species(oldSpec)%MassIC * PartState(iPart,4:6)
        PEM%Element(ElectronIndex)        = PEM%Element(iPart)
        ! Setting the species of the ionized particle
        oldSpec = newSpec
        IF(usevMPF) PartMPF(ElectronIndex) = PartMPF(iPart)
        ! Setting the field for the new particle for the following integration
        FieldAtParticle(ElectronIndex,1:6) = FieldAtParticle(iPart,1:6)
      END IF
      !write(*,*) " "
      !write(*,*) " "
      !write(*,*) " "
    END ASSOCIATE
  END IF
END DO

PDM%ParticleVecLength = PDM%ParticleVecLength + SumOfFormedParticles
PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + SumOfFormedParticles

END SUBROUTINE ADK_Yu2018


END MODULE MOD_PICModels
