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
USE MOD_Part_Tools            ,ONLY: GetNextFreePosition
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iPart, MaxElecQua, ChargedNum, ElectronIndex
REAL                    :: IonizationEnergy_eV, iRan, QuantumTunnelProb, EffQuantNum
REAL                    :: CriticalValue_GV
REAL              :: E_GV
#ifdef CODE_ANALYZE
INTEGER           :: ii,jj
INTEGER,PARAMETER :: NN=16
INTEGER,PARAMETER :: KK=9
REAL(KIND=8)      :: a(NN) = (/(10.0**ii, ii=1,NN, 1)/)
REAL(KIND=8)      :: b(KK) = (/(ii, ii=1,KK, 1)/)
#endif /* CODE_ANALYZE */
!===================================================================================================================================

DO iPart = 1, PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart)) THEN
    ASSOCIATE (&
          oldSpec => PartSpecies(iPart) ,&
          newSpec => SpecDSMC(PartSpecies(iPart))%NextIonizationSpecies )
      IF(newSpec.EQ.0) CYCLE ! skip species that cannot be ionized (e.g. electrons or fully ionized species)
#ifdef CODE_ANALYZE
        DO ii = 1, NN
          DO jj = 1, KK
            E_GV = b(jj)*a(ii) / 1E9  ! [GV/m] for analyis
#else
            E_GV = SQRT(FieldAtParticle(1,iPart)**2 + FieldAtParticle(2,iPart)**2 + FieldAtParticle(3,iPart)**2) / 1E9 ! [GV/m]
#endif /* CODE_ANALYZE */
            ! Ionization energy (same as in QK model)
            MaxElecQua=SpecDSMC(oldSpec)%MaxElecQuant - 1
            IonizationEnergy_eV=SpecDSMC(oldSpec)%ElectronicState(2,MaxElecQua)*BoltzmannConst / ElementaryCharge
            ! Checking applicability of ADK model (normalized variables)
            ! Particle-in-cell simulations of tunneling ionization effects in plasma-based accelerators, David L. Bruhwiler
            CriticalValue_GV = (SQRT(2.) - 1.) * (IonizationEnergy_eV / 27.2)**(1.5) * 5.14E+2
            IF(E_GV.GT.CriticalValue_GV) THEN
              WRITE(UNIT_stdOut,'(A,ES25.14E3)') "IonizationEnergy_eV =", IonizationEnergy_eV
              WRITE(UNIT_stdOut,'(A,ES25.14E3)') "CriticalValue_GV    =", CriticalValue_GV
              WRITE(UNIT_stdOut,'(A,ES25.14E3)') "E_GV    =", E_GV
#ifdef CODE_ANALYZE
              WRITE(UNIT_stdOut,'(A)') "ERROR FieldIonization: ADK model is not applicable for electric fields > critical value!"
#else
              CALL abort(&
                  __STAMP__&
                  ,'ERROR FieldIonization: ADK model is not applicable for electric fields > critical value!', oldSpec)
#endif /* CODE_ANALYZE */
            END IF
            ! Z (ChargedNum): Charge number of the atom/ion AFTER the ionization (thus + 1)
            ChargedNum = NINT(Species(oldSpec)%ChargeIC/ElementaryCharge) + 1
            EffQuantNum = 3.69*REAL(ChargedNum) / SQRT(IonizationEnergy_eV)
            QuantumTunnelProb = 1.52E+15 * 4.**(EffQuantNum)*IonizationEnergy_eV / (EffQuantNum*GAMMA(2.*EffQuantNum)) &
                * (20.5*IonizationEnergy_eV**(3./2.)/E_GV)**(2.*(EffQuantNum-1.)) &
                * EXP(-6.83*IonizationEnergy_eV**(3./2.)/E_GV) * dt
#ifdef CODE_ANALYZE
            CALL WriteFieldIonizationRate(E_GV*1e9,QuantumTunnelProb/dt)
          END DO ! jj = 1, KK
        END DO ! ii = 1, NN
        WRITE (*,*) "\n\n Ionization output for \n ",TRIM(SpecDSMC(oldSpec)%Name)," ==> ",TRIM(SpecDSMC(newSpec)%Name)," + e-\n\n "
        RETURN
#endif /* CODE_ANALYZE */
      CALL RANDOM_NUMBER(iRan)
      IF(QuantumTunnelProb.GT.iRan) THEN
        !.... Get free particle index for the 3rd particle produced
        ElectronIndex = GetNextFreePosition()
        !Set new Species of new particle
        PDM%ParticleInside(ElectronIndex) = .TRUE.
        PDM%isNewPart(ElectronIndex)      = .TRUE.
        PartSpecies(ElectronIndex)        = DSMC%ElectronSpecies
        PartState(1:3,ElectronIndex)      = PartState(1:3,iPart)
        PartState(4:6,ElectronIndex)      = Species(DSMC%ElectronSpecies)%MassIC / Species(oldSpec)%MassIC * PartState(4:6,iPart)
        PartState(4:6,iPart)              = Species(newSpec)%MassIC              / Species(oldSpec)%MassIC * PartState(4:6,iPart)
        PEM%GlobalElemID(ElectronIndex)        = PEM%GlobalElemID(iPart)
        ! Setting the species of the ionized particle
        oldSpec = newSpec
        IF(usevMPF) PartMPF(ElectronIndex) = PartMPF(iPart)
        ! Setting the field for the new particle for the following integration
        FieldAtParticle(1:6,ElectronIndex) = FieldAtParticle(1:6,iPart)
      END IF
    END ASSOCIATE
  END IF
END DO

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
USE MOD_Part_Tools            ,ONLY: GetNextFreePosition
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iPart, MaxElecQua, ElectronIndex
REAL                :: IonizationEnergy_eV, iRan, QuantumTunnelProb
REAL                :: n
#ifdef CODE_ANALYZE
INTEGER             :: ii,jj
INTEGER,PARAMETER   :: NN=16
INTEGER,PARAMETER   :: KK=9
REAL(KIND=8)        :: a(NN) = (/(10.0**ii, ii=1,NN, 1)/)
REAL(KIND=8)        :: b(KK) = (/(ii, ii=1,KK, 1)/)
REAL                :: E
#endif /* CODE_ANALYZE */
!===================================================================================================================================

DO iPart = 1, PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart)) THEN
    ASSOCIATE ( oldSpec => PartSpecies(iPart) ,&
          newSpec => SpecDSMC(PartSpecies(iPart))%NextIonizationSpecies )
      IF(newSpec.EQ.0) CYCLE
      ASSOCIATE (&
            E_au     => 5.1e11 ,& ! [V/m] atomic unit field strength
            omega_au => 4.1e16 ,& ! [1/s] atomic unit frequency strength
            Z        => NINT(Species(oldSpec)%ChargeIC/ElementaryCharge) + 1 & ! Charge number of the ion AFTER the ionization (+1)
#ifdef CODE_ANALYZE
             ) ! [V/m] for analyis
        DO ii = 1, NN
          DO jj = 1, KK
            E = b(jj)*a(ii) ! [V/m] for analyis
#else
           ,E        => SQRT(FieldAtParticle(1,iPart)**2 + FieldAtParticle(2,iPart)**2 + FieldAtParticle(3,iPart)**2) ) ! [V/m]
#endif /* CODE_ANALYZE */
            ! Ionization energy (same as in QK model)
            MaxElecQua=SpecDSMC(oldSpec)%MaxElecQuant - 1
            IonizationEnergy_eV=SpecDSMC(oldSpec)%ElectronicState(2,MaxElecQua)*BoltzmannConst / ElementaryCharge
            n = 3.69*REAL(Z) / SQRT(IonizationEnergy_eV)
            QuantumTunnelProb = 1.61 * omega_au * Z**2 * n**(-9./2.)  &
                * ((10.87 * Z**3 * E_au / (n**4 * E))**(2*n-3./2.))   &
                * EXP(-2*Z**3*E_au / (3*n**3*E))                      &
                * dt
#ifdef CODE_ANALYZE
            CALL WriteFieldIonizationRate(E,QuantumTunnelProb/dt)
          END DO ! jj = 1, KK
        END DO ! ii = 1, NN
        WRITE (*,*) "\n\n Ionization output for \n ",TRIM(SpecDSMC(oldSpec)%Name)," ==> ",TRIM(SpecDSMC(newSpec)%Name)," + e-\n\n "
        RETURN
#endif /* CODE_ANALYZE */
      END ASSOCIATE
      CALL RANDOM_NUMBER(iRan)
      IF(QuantumTunnelProb.GT.iRan) THEN
        !.... Get free particle index for the 3rd particle produced
        ElectronIndex = GetNextFreePosition()
        IF (ElectronIndex.EQ.0) THEN
          CALL abort(__STAMP__,&
              'New Particle Number greater max Part Num in Field Ionization.')
        END IF
        !Set new Species of new particle
        PDM%ParticleInside(ElectronIndex) = .TRUE.
        PDM%isNewPart(ElectronIndex)      = .TRUE.
        PartSpecies(ElectronIndex)        = DSMC%ElectronSpecies
        PartState(1:3,ElectronIndex)      = PartState(1:3,iPart)
        PartState(4:6,ElectronIndex)      = Species(DSMC%ElectronSpecies)%MassIC / Species(oldSpec)%MassIC * PartState(4:6,iPart)
        PartState(4:6,iPart)              = Species(newSpec)%MassIC              / Species(oldSpec)%MassIC * PartState(4:6,iPart)
        PEM%GlobalElemID(ElectronIndex)        = PEM%GlobalElemID(iPart)
        ! Setting the species of the ionized particle
        oldSpec = newSpec
        IF(usevMPF) PartMPF(ElectronIndex) = PartMPF(iPart)
        ! Setting the field for the new particle for the following integration
        FieldAtParticle(1:6,ElectronIndex) = FieldAtParticle(1:6,iPart)
      END IF
    END ASSOCIATE
  END IF
END DO

END SUBROUTINE ADK_Yu2018


#ifdef CODE_ANALYZE
!----------------------------------------------------------------------------------------------------------------------------------!
!> Write electric field strength and ionization rate to FieldIonizationRate.csv file
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE WriteFieldIonizationRate(E,W)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals               ,ONLY: MPIRoot,FILEEXISTS,unit_stdout
USE MOD_Restart_Vars          ,ONLY: DoRestart
USE MOD_Globals               ,ONLY: abort
USE MOD_TimeDisc_Vars         ,ONLY: iter
!----------------------------------------------------------------------------------------------------------------------------------!                                                                    ! ----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)                  :: E ! Electric Field Strength [V/m]
REAL,INTENT(IN)                  :: W ! Ionization Rate [1/s]
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=28),PARAMETER              :: outfile='FieldIonizationRate.csv'
INTEGER                                  :: ioUnit,I
CHARACTER(LEN=150)                       :: formatStr
INTEGER,PARAMETER                        :: nOutputVar=2
CHARACTER(LEN=255),DIMENSION(nOutputVar) :: StrVarNames(nOutputVar)=(/ CHARACTER(LEN=255) :: &
    'E[V/m]',     &
    'W[1/s]'      &
    /)
CHARACTER(LEN=255),DIMENSION(nOutputVar) :: tmpStr ! needed because PerformAnalyze is called multiple times at the beginning
CHARACTER(LEN=1000)                      :: tmpStr2
CHARACTER(LEN=1),PARAMETER               :: delimiter=","
LOGICAL                                  :: FileExist,CreateFile
!===================================================================================================================================
! only the root shall write this file
IF(.NOT.MPIRoot)THEN
  CALL abort(__STAMP__,&
      'Field Ionization: The output of field ionization data can only be run single core (MPI not implemented)!')
END IF

! check if file is to be created
CreateFile=.TRUE.
IF(iter.GT.0)CreateFile=.FALSE.                             ! don't create new file if this is not the first iteration
IF((DoRestart).AND.(FILEEXISTS(outfile)))CreateFile=.FALSE. ! don't create new file if this is a restart and the file already exists
!                                                           ! assume continued simulation and old load balance data is still needed

! check if new file with header is to be created
INQUIRE(FILE = outfile, EXIST=FileExist)
IF(.NOT.FileExist)CreateFile=.TRUE.                         ! if no file exists, create one

! create file with header
IF(CreateFile) THEN
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfile),STATUS="UNKNOWN")
  tmpStr=""
  DO I=1,nOutputVar
    WRITE(tmpStr(I),'(A)')delimiter//'"'//TRIM(StrVarNames(I))//'"'
  END DO
  WRITE(formatStr,'(A1)')'('
  DO I=1,nOutputVar
    IF(I.EQ.nOutputVar)THEN ! skip writing "," and the end of the line
      WRITE(formatStr,'(A,A1,I2)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(I))
    ELSE
      WRITE(formatStr,'(A,A1,I2,A1)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(I)),','
    END IF
  END DO

  WRITE(formatStr,'(A,A1)')TRIM(formatStr),')' ! finish the format
  WRITE(tmpStr2,formatStr)tmpStr               ! use the format and write the header names to a temporary string
  tmpStr2(1:1) = " "                           ! remove possible relimiter at the beginning (e.g. a comma)
  WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2))    ! clip away the front and rear white spaces of the temporary string

  CLOSE(ioUnit)
  iter=iter+1
END IF

! Print info to file
IF(FILEEXISTS(outfile))THEN
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfile),POSITION="APPEND",STATUS="OLD")
      WRITE(formatStr,'(A2,I2,A14,A1)')'(',nOutputVar,CSVFORMAT,')'
  WRITE(tmpStr2,formatStr)&
      " ",E, &           ! Electric field strength
      delimiter,W        ! Ionization rate
  WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2)) ! clip away the front and rear white spaces of the data line
  CLOSE(ioUnit)
ELSE
  SWRITE(UNIT_StdOut,'(A)')TRIM(outfile)//" does not exist. Cannot write particle tracking info!"
END IF

END SUBROUTINE WriteFieldIonizationRate
#endif /* CODE_ANALYZE */

END MODULE MOD_PICModels
