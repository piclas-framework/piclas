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

MODULE MOD_MCC_Init
!===================================================================================================================================
!> Initialization of the Monte Carlo Collision module
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

PUBLIC :: DefineParametersMCC, InitMCC, MCC_Chemistry_Init, FinalizeMCC
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for MCC
!==================================================================================================================================
SUBROUTINE DefineParametersMCC()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("MCC")

CALL prms%CreateStringOption(   'Particles-CollXSec-Database', 'File name for the collision cross section database. Container '//&
                                                               'should be named with species pair (e.g. "Ar-electron"). The '//&
                                                               'first column shall contain the energy in eV and the second '//&
                                                               'column the cross-section in m^2', 'none')
CALL prms%CreateLogicalOption(  'Particles-CollXSec-NullCollision'  &
                                  ,'Utilize the null collision method for the determination of the number of pairs '//&
                                  'based on the maximum collision frequency and time step (only with a background gas)' &
                                  ,'.TRUE.')

CALL prms%CreateLogicalOption(  'Part-Species[$]-UseCollXSec'  &
                                           ,'Utilize collision cross sections for the determination of collision probabilities' &
                                           ,'.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Species[$]-UseVibXSec'  &
                                           ,'Utilize vibrational cross sections for the determination of relaxation probabilities' &
                                           ,'.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Species[$]-UseElecXSec'  &
                                           ,'Utilize electronic cross sections, only in combination with ElectronicModel = 3' &
                                           ,'.FALSE.', numberedmulti=.TRUE.)

END SUBROUTINE DefineParametersMCC


SUBROUTINE InitMCC()
!===================================================================================================================================
!> Initialization & read-in of the cross-section database for:
!> - Collisional, vibrational and electronic relaxation
!> - Initialization of the null collision method
!> - Chemical reactions and photo-ionization
!> NOTE: Must be called after InitDSMC, which contains the initialization of the chemistry module.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Globals_Vars  ,ONLY: ElementaryCharge
USE MOD_Particle_Vars ,ONLY: nSpecies, SampleElecExcitation, ExcitationLevelCounter, ExcitationSampleData, ExcitationLevelMapping
USE MOD_Particle_Vars ,ONLY: VarTimeStep
USE MOD_Mesh_Vars     ,ONLY: nElems
USE MOD_DSMC_Vars     ,ONLY: BGGas, SpecDSMC, CollInf, DSMC, ChemReac, CollisMode
USE MOD_MCC_Vars      ,ONLY: UseMCC, XSec_Database, SpecXSec, XSec_NullCollision, XSec_Relaxation
USE MOD_MCC_Vars      ,ONLY: NbrOfPhotonXsecReactions
USE MOD_MCC_XSec      ,ONLY: ReadCollXSec, ReadVibXSec, InterpolateCrossSection, ReadElecXSec
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
#if defined(PARTICLES) && USE_HDG
USE MOD_HDG_Vars      ,ONLY: UseBRElectronFluid,BRNullCollisionDefault
USE MOD_ReadInTools   ,ONLY: PrintOption
#endif /*defined(PARTICLES) && USE_HDG*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=3)      :: hilf
INTEGER               :: iSpec, jSpec, iCase, partSpec
REAL                  :: TotalProb(nSpecies), CrossSection, MaxEnergyColl
INTEGER               :: iLevel, nVib, iStep, MaxDim, MaxDimLevel, NumLevel, MaxLevelIndex
INTEGER, ALLOCATABLE  :: CounterLevelAboveMax(:)
REAL, ALLOCATABLE     :: CollXSecDataTemp(:,:)
!===================================================================================================================================

UseMCC             = .FALSE.
XSec_NullCollision = .FALSE.
XSec_Relaxation    = .FALSE.

IF(CollisMode.EQ.0) RETURN

!-----------------------------------------------------------------------------------------------------------------------------------
! Determine whether cross-section data is utilized
!-----------------------------------------------------------------------------------------------------------------------------------
DO iSpec = 1, nSpecies
  WRITE(UNIT=hilf,FMT='(I0)') iSpec
  SpecDSMC(iSpec)%UseCollXSec=GETLOGICAL('Part-Species'//TRIM(hilf)//'-UseCollXSec')
  SpecDSMC(iSpec)%UseVibXSec=GETLOGICAL('Part-Species'//TRIM(hilf)//'-UseVibXSec')
  SpecDSMC(iSpec)%UseElecXSec=GETLOGICAL('Part-Species'//TRIM(hilf)//'-UseElecXSec')
  IF(SpecDSMC(iSpec)%UseCollXSec.AND.BGGas%BackgroundSpecies(iSpec)) THEN
    CALL Abort(__STAMP__,'ERROR: Please supply the collision cross-section flag for the particle species and NOT the background species!')
  END IF
  IF(SpecDSMC(iSpec)%UseElecXSec.AND.SpecDSMC(iSpec)%InterID.EQ.4) THEN
    CALL Abort(__STAMP__,'ERROR: Electronic relaxation should be enabled for the respective heavy species, not the electrons!')
  END IF
  IF(.NOT.DSMC%ReservoirSimu) THEN
    IF(SpecDSMC(iSpec)%UseElecXSec.AND.(.NOT.BGGas%BackgroundSpecies(iSpec))) THEN
      SWRITE(*,*) 'NOTE: Electronic relaxation via cross-sections for regular DSMC is currently only enabled using the flag'
      SWRITE(*,*) 'NOTE: Particles-DSMCReservoirSim = T to test the calculation of the probabilities during regression testing.'
      SWRITE(*,*) 'NOTE: For DSMC, a de-excitation model should be implemented. Regression test: WEK_Reservoir/MCC_N2_XSec_Elec'
      CALL Abort(__STAMP__,'ERROR: Electronic relaxation via cross-section (-UseElecXSec) is only supported with a background gas!')
    END IF
  END IF
END DO

! Enable MCC if any of the flags was set
IF(ANY(SpecDSMC(:)%UseCollXSec).OR.ANY(SpecDSMC(:)%UseVibXSec).OR.ANY(SpecDSMC(:)%UseElecXSec).OR.ChemReac%AnyXSecReaction &
    .OR.(NbrOfPhotonXsecReactions.GT.0)) THEN
  UseMCC = .TRUE.
END IF

! Ambipolar diffusion is not implemented with the regular background gas, only with MCC
IF(DSMC%DoAmbipolarDiff) THEN
  IF((BGGas%NumberOfSpecies.GT.0).AND.(.NOT.UseMCC)) CALL abort(__STAMP__,&
      'ERROR: Ambipolar diffusion is not implemented with the regular background gas!')
END IF

! MCC and variable vibrational relaxation probability is not supported
IF(UseMCC.AND.(DSMC%VibRelaxProb.EQ.2.0)) CALL abort(__STAMP__&
      ,'ERROR: Monte Carlo Collisions and variable vibrational relaxation probability (DSMC-based) are not compatible!')

IF(.NOT.UseMCC.AND.VarTimeStep%UseSpeciesSpecific) CALL abort(__STAMP__,&
      'ERROR: Only MCC is implemented with a species-specific time step!')

! Disable SampleElecExcitation if electronic excitation has not been enabled for at least one species
IF(SampleElecExcitation.AND.(.NOT.ANY(SpecDSMC(:)%UseElecXSec))) THEN
  SampleElecExcitation = .FALSE.
  LBWRITE(*,*) '| WARNING: Part-SampleElectronicExcitation has been disabled as no electronic excitation has been enabled through -UseElecXSec = T!'
END IF
! Leave the routine
IF(.NOT.UseMCC) RETURN

!-----------------------------------------------------------------------------------------------------------------------------------
! Initialize & read-in of cross-section data
!-----------------------------------------------------------------------------------------------------------------------------------
! Read-in of the cross-section database
XSec_Database = GETSTR('Particles-CollXSec-Database')
! Checks/read-in depending on background gas
IF(BGGas%NumberOfSpecies.GT.0) THEN
  ! Null collision method only works with a background gas
  XSec_NullCollision = GETLOGICAL('Particles-CollXSec-NullCollision')
  ! Disable variable time step for MCC but keep it for the push (parameter is defined in DefineParametersVariableTimeStep)
  IF(VarTimeStep%UseSpeciesSpecific) VarTimeStep%DisableForMCC = GETLOGICAL('Part-VariableTimeStep-DisableForMCC')
ELSE
  ! Variable time step is only implemented with a background gas
  IF(VarTimeStep%UseSpeciesSpecific) CALL abort(__STAMP__,'ERROR: Only MCC is implemented with a species-specific time step!')
END IF

ALLOCATE(SpecXSec(CollInf%NumCase))
SpecXSec(:)%UseCollXSec = .FALSE.
SpecXSec(:)%UseVibXSec = .FALSE.
SpecXSec(:)%UseElecXSec = .FALSE.
SpecXSec(:)%NumElecLevel = 0
SpecXSec(:)%CollXSec_Effective = .FALSE.
SpecXSec(:)%SpeciesToRelax = 0
SpecXSec(:)%ProbNull = 0.
SpecXSec(:)%NumElecLevel = 0
TotalProb = 0.

DO iSpec = 1, nSpecies
  DO jSpec = iSpec, nSpecies
    iCase = CollInf%Coll_Case(iSpec,jSpec)
    IF(XSec_NullCollision)THEN
      IF(BGGas%UseDistribution)THEN
        ALLOCATE(SpecXSec(iCase)%ProbNullElem(1:nElems))
        SpecXSec(iCase)%ProbNullElem = 0.
      END IF ! BGGas%UseDistribution
    END IF ! XSec_NullCollision
    ! Skip species, which shall not be treated with collision cross-sections
    IF(.NOT.SpecDSMC(iSpec)%UseCollXSec.AND..NOT.SpecDSMC(jSpec)%UseCollXSec.AND. &
       .NOT.SpecDSMC(iSpec)%UseVibXSec.AND..NOT.SpecDSMC(jSpec)%UseVibXSec.AND. &
       .NOT.SpecDSMC(iSpec)%UseElecXSec.AND..NOT.SpecDSMC(jSpec)%UseElecXSec) CYCLE
    ! Skip pairing with itself and pairing with other particle species, if background gas is active
    IF(BGGas%NumberOfSpecies.GT.0) THEN
      IF(iSpec.EQ.jSpec) CYCLE
      IF(.NOT.BGGas%BackgroundSpecies(iSpec).AND..NOT.BGGas%BackgroundSpecies(jSpec)) CYCLE
    END IF
    ! Read-in cross-section data for collisions of particles, allocating CollXSecData within the following routine
    IF(SpecDSMC(iSpec)%UseCollXSec.OR.SpecDSMC(jSpec)%UseCollXSec) CALL ReadCollXSec(iCase, iSpec, jSpec)
    ! Check if both species were given the UseCollXSec flag and store the energy value in Joule
    IF(SpecXSec(iCase)%UseCollXSec) THEN
      IF(SpecDSMC(iSpec)%UseCollXSec.AND.SpecDSMC(jSpec)%UseCollXSec) CALL abort(__STAMP__&
          ,'ERROR: Both species defined to use collisional cross-section, define only the source species with UseCollXSec!')
      ! Store the energy value in J (read-in was in eV)
      SpecXSec(iCase)%CollXSecData(1,:) = SpecXSec(iCase)%CollXSecData(1,:) * ElementaryCharge
    END IF
    ! Read-in vibrational cross sections
    IF(SpecDSMC(iSpec)%UseVibXSec.OR.SpecDSMC(jSpec)%UseVibXSec) CALL ReadVibXSec(iCase, iSpec, jSpec)
    ! Vibrational relaxation probabilities: Interpolate and store the probability at the collision cross-section levels
    IF(SpecXSec(iCase)%UseVibXSec) THEN
      IF(SpecDSMC(iSpec)%UseVibXSec.AND.SpecDSMC(jSpec)%UseVibXSec) CALL abort(__STAMP__&
          ,'ERROR: Both species defined to use vib. cross-section, define only the source species with UseVibXSec!')
      ! Save which species shall use the vibrational cross-section data for relaxation probabilities
      ! If the species which was given the UseVibXSec flag is diatomic/polyatomic, use the cross-section for that species
      ! If the species is an atom/electron, use the cross-section for the other collision partner (the background species)
      IF(SpecDSMC(iSpec)%UseVibXSec) THEN
        IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
          SpecXSec(iCase)%SpeciesToRelax = iSpec
        ELSE
          SpecXSec(iCase)%SpeciesToRelax = jSpec
        END IF
      ELSE
        IF((SpecDSMC(jSpec)%InterID.EQ.2).OR.(SpecDSMC(jSpec)%InterID.EQ.20)) THEN
          SpecXSec(iCase)%SpeciesToRelax = jSpec
        ELSE
          SpecXSec(iCase)%SpeciesToRelax = iSpec
        END IF
      END IF
      XSec_Relaxation = .TRUE.
      SpecXSec(iCase)%VibCount = 0.
      nVib = SIZE(SpecXSec(iCase)%VibMode)
      DO iLevel = 1, nVib
        ! Store the energy value in J (read-in was in eV)
        SpecXSec(iCase)%VibMode(iLevel)%XSecData(1,:) = SpecXSec(iCase)%VibMode(iLevel)%XSecData(1,:) * ElementaryCharge
        SpecXSec(iCase)%VibMode(iLevel)%Threshold = SpecXSec(iCase)%VibMode(iLevel)%Threshold * ElementaryCharge
      END DO
      IF(SpecXSec(iCase)%UseCollXSec) THEN
        ! Collision cross-sections are available
        ! Array bounds of collisional cross-section data
        MaxDim = UBOUND(SpecXSec(iCase)%CollXSecData,dim=2)
        ! Maximum energy value of cross-section data
        MaxEnergyColl = SpecXSec(iCase)%CollXSecData(1,MaxDim)
        ALLOCATE(CounterLevelAboveMax(nVib))
        CounterLevelAboveMax = 0
        DO iLevel = 1, nVib
          ! Determine whether the maximal energy level of vibrational levels is greater than the provided collisional level
          CounterLevelAboveMax(iLevel) = COUNT(SpecXSec(iCase)%VibMode(iLevel)%XSecData(1,:).GT.MaxEnergyColl)
        END DO
        ! Abort in case of EFFECTIVE or resizing the CollXSecData array in case of ELASTIC, if vibrational energy levels are above
        IF(ANY(CounterLevelAboveMax.GT.0)) THEN
          IF(SpecXSec(iCase)%CollXSec_Effective) THEN
            CALL abort(__STAMP__,'ERROR: Effective cross-section has a smaller energy range than the vibrational level cross-section!')
          ELSE
            ! Get the index of vibrational level with the most energy values above
            MaxLevelIndex = MAXLOC(CounterLevelAboveMax,DIM=1)
            NumLevel = MAXVAL(CounterLevelAboveMax)
            ! Allocate a temporary array for the new collisional cross-section data
            ALLOCATE(CollXSecDataTemp(1:2,1:MaxDim+NumLevel))
            ! Store the old array in the new temporary array
            CollXSecDataTemp(1:2,1:MaxDim) = SpecXSec(iCase)%CollXSecData(1:2,1:MaxDim)
            ! Determine the bounds of vibrational energy levels to be added to the collision array
            MaxDimLevel = UBOUND(SpecXSec(iCase)%VibMode(MaxLevelIndex)%XSecData,DIM=2)
            ! Add the energy levels but not the cross-section, it will added in the next step
            CollXSecDataTemp(1,MaxDim+1:MaxDim+NumLevel) = SpecXSec(iCase)%VibMode(MaxLevelIndex)%XSecData(1,MaxDimLevel-NumLevel+1:MaxDimLevel)
            CollXSecDataTemp(2,MaxDim+1:MaxDim+NumLevel) = 0.
            CALL MOVE_ALLOC(CollXSecDataTemp,SpecXSec(iCase)%CollXSecData)
            LBWRITE(*,*) 'Resizing CollXSecData array from ', MaxDim, ' to ', UBOUND(SpecXSec(iCase)%CollXSecData,dim=2)
            MaxDim = UBOUND(SpecXSec(iCase)%CollXSecData,dim=2)
          END IF
        END IF
        ! Interpolate the vibrational cross section at the energy levels of the collision collision cross section and sum-up the
        ! vibrational probability (vibrational cross-section divided by the effective)
        DO iStep = 1, MaxDim
          DO iLevel = 1, nVib
            CrossSection = InterpolateCrossSection(SpecXSec(iCase)%VibMode(iLevel)%XSecData,SpecXSec(iCase)%CollXSecData(1,iStep))
            ! When no effective cross-section is available, the vibrational cross-section has to be added to the collisional
            IF(SpecXSec(iCase)%CollXSec_Effective) THEN
              IF(CrossSection.GT.SpecXSec(iCase)%CollXSecData(2,iStep)) THEN
                SWRITE(*,*) 'Current vibrational energy level [eV]: ', SpecXSec(iCase)%VibMode(iLevel)%Threshold / ElementaryCharge
                SWRITE(*,*) 'Vibrational cross-section: ', CrossSection
                SWRITE(*,*) 'Effective cross-section: ', SpecXSec(iCase)%CollXSecData(2,iStep)
                SWRITE(*,*) 'Effective cross-section should be greater as the vibrational is supposed to be part of the effective cross-section.'
                SWRITE(*,*) 'Check the last value of the vibrational data, the cross-section should be zero, otherwise the last value will be taken for energies outside the vibrational data.'
                CALL abort(__STAMP__,'ERROR: Effective cross-section is smaller than the interpolated vibrational level cross-section!')
              END IF
            ELSE
              SpecXSec(iCase)%CollXSecData(2,iStep) = SpecXSec(iCase)%CollXSecData(2,iStep) + CrossSection
            END IF
          END DO
        END DO
        DEALLOCATE(CounterLevelAboveMax)
      END IF    ! SpecXSec(iCase)%UseCollXSec
    END IF      ! SpecXSec(iCase)%UseVibXSec
    ! Read-in electronic cross-section data
    IF(DSMC%ElectronicModel.EQ.3) THEN
      ! Read-in electronic level cross-section data if flags have been defined
      IF(SpecDSMC(iSpec)%UseElecXSec.OR.SpecDSMC(jSpec)%UseElecXSec) CALL ReadElecXSec(iCase, iSpec, jSpec)
      IF(SpecXSec(iCase)%UseElecXSec) THEN
        ! Check if only heavy-species - electron combinations were found
        IF(SpecDSMC(iSpec)%UseElecXSec.AND.SpecDSMC(jSpec)%UseElecXSec) THEN
          CALL abort(__STAMP__,'ERROR: Electronic excitation using cross-section data is currently only supported through electron collisions!')
        END IF
        DO iLevel = 1, SpecXSec(iCase)%NumElecLevel
          ! Store the energy value in J (read-in was in eV)
          SpecXSec(iCase)%ElecLevel(iLevel)%XSecData(1,:) = SpecXSec(iCase)%ElecLevel(iLevel)%XSecData(1,:) * ElementaryCharge
          SpecXSec(iCase)%ElecLevel(iLevel)%Threshold = SpecXSec(iCase)%ElecLevel(iLevel)%Threshold * ElementaryCharge
        END DO
        ! Interpolate and store levels at the collision cross-section intervals
        IF(SpecXSec(iCase)%UseCollXSec) THEN
          IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(SpecDSMC(jSpec)%InterID.NE.4)) THEN
            ! Special treatment required if both collision partners have electronic energy levels (ie. one is not an electron)
            CALL abort(__STAMP__,'ERROR: Electronic relaxation with cross-section is only possible for electron collisions!')
          END IF
          ! Collision cross-sections are available
          ! Array bounds of collisional cross-section data
          MaxDim = UBOUND(SpecXSec(iCase)%CollXSecData,dim=2)
          ! Maximum energy value of cross-section data
          MaxEnergyColl = SpecXSec(iCase)%CollXSecData(1,MaxDim)
          ALLOCATE(CounterLevelAboveMax(SpecXSec(iCase)%NumElecLevel))
          CounterLevelAboveMax = 0
          DO iLevel = 1, SpecXSec(iCase)%NumElecLevel
            ! Determine whether the maximal energy level of electronic levels is greater than the provided collisional level
            CounterLevelAboveMax(iLevel) = COUNT(SpecXSec(iCase)%ElecLevel(iLevel)%XSecData(1,:).GT.MaxEnergyColl)
          END DO
          ! Abort in case of EFFECTIVE or resizing the CollXSecData array in case of ELASTIC, if electronic energy levels are above
          IF(ANY(CounterLevelAboveMax.GT.0)) THEN
            IF(SpecXSec(iCase)%CollXSec_Effective) THEN
              CALL abort(__STAMP__,'ERROR: Effective cross-section has a smaller energy range than the electronic level cross-section!')
            ELSE
              ! Get the index of electronic level with the most energy values above
              MaxLevelIndex = MAXLOC(CounterLevelAboveMax,DIM=1)
              NumLevel = MAXVAL(CounterLevelAboveMax)
              ! Allocate a temporary array for the new collisional cross-section data
              ALLOCATE(CollXSecDataTemp(1:2,1:MaxDim+NumLevel))
              ! Store the old array in the new temporary array
              CollXSecDataTemp(1:2,1:MaxDim) = SpecXSec(iCase)%CollXSecData(1:2,1:MaxDim)
              ! Determine the bounds of electronic energy levels to be added to the collision array
              MaxDimLevel = UBOUND(SpecXSec(iCase)%ElecLevel(MaxLevelIndex)%XSecData,DIM=2)
              ! Add the energy levels but not the cross-section, it will added in the next step
              CollXSecDataTemp(1,MaxDim+1:MaxDim+NumLevel) = SpecXSec(iCase)%ElecLevel(MaxLevelIndex)%XSecData(1,MaxDimLevel-NumLevel+1:MaxDimLevel)
              CollXSecDataTemp(2,MaxDim+1:MaxDim+NumLevel) = 0.
              CALL MOVE_ALLOC(CollXSecDataTemp,SpecXSec(iCase)%CollXSecData)
              LBWRITE(*,*) 'Resizing CollXSecData array from ', MaxDim, ' to ', UBOUND(SpecXSec(iCase)%CollXSecData,dim=2)
              MaxDim = UBOUND(SpecXSec(iCase)%CollXSecData,dim=2)
            END IF
          END IF
          ! Interpolate the electronic cross section at the energy levels of the collision collision cross section and sum-up the
          ! electronic probability (electronic cross-section divided by the effective)
          DO iStep = 1, MaxDim
            DO iLevel = 1, SpecXSec(iCase)%NumElecLevel
              CrossSection = InterpolateCrossSection(SpecXSec(iCase)%ElecLevel(iLevel)%XSecData,SpecXSec(iCase)%CollXSecData(1,iStep))
              ! When no effective cross-section is available, the electronic cross-section has to be added to the collisional
              IF(SpecXSec(iCase)%CollXSec_Effective) THEN
                IF(CrossSection.GT.SpecXSec(iCase)%CollXSecData(2,iStep)) THEN
                  SWRITE(*,*) 'Current electronic energy level [eV]: ', SpecXSec(iCase)%ElecLevel(iLevel)%Threshold / ElementaryCharge
                  SWRITE(*,*) 'Electronic cross-section: ', CrossSection
                  SWRITE(*,*) 'Effective cross-section: ', SpecXSec(iCase)%CollXSecData(2,iStep)
                  SWRITE(*,*) 'Effective cross-section should be greater as the electronic is supposed to be part of the effective cross-section.'
                  SWRITE(*,*) 'Check the last value of the electronic data, it should be zero, otherwise the last value will be taken for energies outside the electronic data.'
                  CALL abort(__STAMP__,'ERROR: Effective cross-section is smaller than the interpolated electronic level cross-section!')
                END IF
              ELSE
                SpecXSec(iCase)%CollXSecData(2,iStep) = SpecXSec(iCase)%CollXSecData(2,iStep) + CrossSection
              END IF
            END DO
          END DO
          DEALLOCATE(CounterLevelAboveMax)
        END IF    ! SpecXSec(iCase)%UseCollXSec
      END IF      ! SpecXSec(iCase)%UseElecXSec
    END IF
    ! Determine the maximum collision frequency for the null collision method
    IF(SpecXSec(iCase)%UseCollXSec) THEN
      IF(XSec_NullCollision) THEN
        CALL DetermineNullCollProb(iCase,iSpec,jSpec)
        ! Select the particle species in order to sum-up the total null collision probability per particle species
        IF(BGGas%BackgroundSpecies(iSpec)) THEN
          partSpec = jSpec
        ELSE
          partSpec = iSpec
        END IF
        IF(BGGas%UseDistribution)THEN
          TotalProb(partSpec) = TotalProb(partSpec) + MAXVAL(SpecXSec(iCase)%ProbNullElem)
        ELSE
          TotalProb(partSpec) = TotalProb(partSpec) + SpecXSec(iCase)%ProbNull
        END IF ! BGGas%UseDistribution
        ! Sum of null collision probability per particle species should be lower than 1, otherwise not enough collision pairs
        IF(TotalProb(partSpec).GT.1.0) THEN
          CALL abort(__STAMP__,'Total null collision probability is above unity. Please reduce the time step! Probability is: '&
          ,RealInfoOpt=TotalProb(partSpec))
        ELSEIF(TotalProb(partSpec).GT.0.1) THEN
          LBWRITE(*,*) 'Total null collision probability is above 0.1. A value of 1E-2 is recommended in literature!'
          LBWRITE(*,*) 'Particle Species: ', TRIM(SpecDSMC(partSpec)%Name), ' Probability: ', TotalProb(partSpec)
        END IF ! TotalProb(partSpec).GT.1.0
      END IF ! XSec_NullCollision
    END IF ! SpecXSec(iCase)%UseCollXSec
  END DO ! jSpec = iSpec, nSpecies
END DO ! iSpec = 1, nSpecies

! Read-in of the reaction cross-section database and re-calculation of the null collision probability
IF (CollisMode.EQ.3) THEN
  ! Dissociation, exchange, and ionization through cross-section data
  IF(ChemReac%AnyXSecReaction) THEN
    CALL MCC_Chemistry_Init()
  END IF
  ! Photo-ionization through cross-section data
  IF(NbrOfPhotonXsecReactions.GT.0) THEN
    CALL InitPhotoionizationXSec()
    CALL CheckPhotoionizationXSec()
  END IF
END IF

! Allocation of electronic excitation array
IF(SampleElecExcitation) THEN
  ExcitationLevelCounter = 0
  ALLOCATE(ExcitationLevelMapping(CollInf%NumCase,MAXVAL(SpecXSec(:)%NumElecLevel)))
  ! Count the number of excitation levels for all collision cases
  DO iCase = 1, CollInf%NumCase
    IF(.NOT.SpecXSec(iCase)%UseElecXSec) CYCLE
    DO iLevel = 1, SpecXSec(iCase)%NumElecLevel
      ExcitationLevelCounter = ExcitationLevelCounter + 1
      ExcitationLevelMapping(iCase,iLevel) = ExcitationLevelCounter
    END DO
  END DO
  ! Sanity check
  IF(ExcitationLevelCounter.NE.SUM(SpecXSec(:)%NumElecLevel)) THEN
    IPWRITE(UNIT_StdOut,*) "ExcitationLevelCounter        =", ExcitationLevelCounter
    IPWRITE(UNIT_StdOut,*) "SUM(SpecXSec(:)%NumElecLevel) =", SUM(SpecXSec(:)%NumElecLevel)
    CALL abort(__STAMP__,'Electronic excitation sampling: Wrong level counter!')
  END IF
  ALLOCATE(ExcitationSampleData(ExcitationLevelCounter,nElems))
  ExcitationSampleData = 0.
END IF

#if defined(PARTICLES) && USE_HDG
BRNullCollisionDefault = XSec_NullCollision ! Backup read-in parameter value (for switching null collision on/off)
IF(XSec_NullCollision.AND.UseBRElectronFluid)THEN
  XSec_NullCollision = .FALSE. ! Deactivate null collision when using BR electrons due to (possibly) increased time step
  CALL PrintOption('Using BR electron fluid model: Particles-CollXSec-NullCollision','INFO',LogOpt=XSec_NullCollision)
END IF
#endif /*defined(PARTICLES) && USE_HDG*/

END SUBROUTINE InitMCC


SUBROUTINE DetermineNullCollProb(iCase,iSpec,jSpec)
!===================================================================================================================================
!> Routine for the MCC method: calculates the maximal collision frequency for a given species and the collision probability
!> Utilizing the ManualTimeStep as dt is not yet defined.
!===================================================================================================================================
! MODULES
USE MOD_ReadInTools
USE MOD_Globals_Vars          ,ONLY: Pi
USE MOD_Mesh_Vars             ,ONLY: nElems
USE MOD_Particle_Vars         ,ONLY: Species, VarTimeStep
USE MOD_TimeDisc_Vars         ,ONLY: ManualTimeStep
USE MOD_DSMC_Vars             ,ONLY: BGGas
USE MOD_MCC_Vars              ,ONLY: SpecXSec
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)            :: iCase                            !< Case index
INTEGER,INTENT(IN)            :: iSpec
INTEGER,INTENT(IN)            :: jSpec
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: MaxDOF, bggSpec, iElem, partSpec
REAL                          :: MaxCollFreq, MaxCollFreqElem, Mass, dtVar
REAL,ALLOCATABLE              :: Velocity(:)
!===================================================================================================================================

! Select the background species as the target cloud and use the mass of particle species
IF(BGGas%BackgroundSpecies(iSpec)) THEN
  bggSpec = BGGas%MapSpecToBGSpec(iSpec)
  Mass = Species(jSpec)%MassIC
  partSpec = jSpec
ELSE
  bggSpec = BGGas%MapSpecToBGSpec(jSpec)
  Mass = Species(iSpec)%MassIC
  partSpec = iSpec
END IF

! Set the correct time step for the particle species
IF(VarTimeStep%UseSpeciesSpecific.AND..NOT.VarTimeStep%DisableForMCC) THEN
  dtVar = ManualTimeStep * Species(partSpec)%TimeStepFactor
ELSE
  dtVar = ManualTimeStep
END IF

MaxDOF = SIZE(SpecXSec(iCase)%CollXSecData,2)
ALLOCATE(Velocity(MaxDOF))

! Determine the mean relative velocity at the given energy level
Velocity(1:MaxDOF) = SQRT(2.) * SQRT(8.*SpecXSec(iCase)%CollXSecData(1,1:MaxDOF)/(Pi*Mass))

! Calculate the maximal collision frequency: Step 1
MaxCollFreq = MAXVAL(Velocity(1:MaxDOF) * SpecXSec(iCase)%CollXSecData(2,1:MaxDOF))

IF(BGGas%UseDistribution) THEN
  DO iElem = 1, nElems
    ! Calculate the maximal collision frequency: Step 2
    MaxCollFreqElem = MaxCollFreq * BGGas%Distribution(bggSpec,7,iElem)
    ! Determine the collision probability
    SpecXSec(iCase)%ProbNullElem(iElem) = 1. - EXP(-MaxCollFreqElem*dtVar)
  END DO ! iElem = 1, nElems
ELSE
  ! Calculate the maximal collision frequency: Step 2
  MaxCollFreq = MaxCollFreq * BGGas%NumberDensity(bggSpec)
  ! Determine the collision probability
  SpecXSec(iCase)%ProbNull = 1. - EXP(-MaxCollFreq*dtVar)
END IF

DEALLOCATE(Velocity)

END SUBROUTINE DetermineNullCollProb


SUBROUTINE MCC_Chemistry_Init()
!===================================================================================================================================
!> Read-in of the reaction cross-section database and re-calculation of the null collision probability
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_MCC_XSec      ,ONLY: ReadReacXSec, InterpolateCrossSection
USE MOD_PARTICLE_Vars ,ONLY: nSpecies
USE MOD_DSMC_Vars     ,ONLY: BGGas, CollInf, ChemReac
USE MOD_MCC_Vars      ,ONLY: SpecXSec, XSec_NullCollision
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSpec, jSpec, iCase, iReac
REAL                  :: TotalProb, ReactionCrossSection, MaxEnergyColl
INTEGER               :: iStep, MaxDim, MaxLevelIndex, NumLevel, MaxDimLevel
INTEGER               :: iPath, NumPaths
INTEGER, ALLOCATABLE  :: CounterLevelAboveMax(:)
REAL, ALLOCATABLE     :: CollXSecDataTemp(:,:)
!===================================================================================================================================

IF(BGGas%NumberOfSpecies.LE.0) CALL abort(__STAMP__,'Cross-section-based chemistry without background gas has not been tested yet!')

! 1.) Read-in of cross-section data for chemical reactions
DO iCase = 1, CollInf%NumCase
  NumPaths = ChemReac%CollCaseInfo(iCase)%NumOfReactionPaths
  IF(ChemReac%CollCaseInfo(iCase)%HasXSecReaction) ALLOCATE(SpecXSec(iCase)%ReactionPath(1:NumPaths))
  DO iPath = 1, NumPaths
    iReac = ChemReac%CollCaseInfo(iCase)%ReactionIndex(iPath)
    IF(TRIM(ChemReac%ReactModel(iReac)).EQ.'XSec') THEN
      CALL ReadReacXSec(iCase,iPath)
    END IF
  END DO
END DO

! 2.) Add the chemical reaction cross-section to the total collision cross-section
DO iCase = 1, CollInf%NumCase
  ! Collision cross-sections are available
  IF(SpecXSec(iCase)%UseCollXSec) THEN
    ! Skip cases without reactions
    IF(.NOT.ChemReac%CollCaseInfo(iCase)%HasXSecReaction) CYCLE
    ! Array bounds of collisional cross-section data
    MaxDim = UBOUND(SpecXSec(iCase)%CollXSecData,dim=2)
    ! Maximum energy value of cross-section data
    MaxEnergyColl = SpecXSec(iCase)%CollXSecData(1,MaxDim)
    NumPaths = ChemReac%CollCaseInfo(iCase)%NumOfReactionPaths
    ALLOCATE(CounterLevelAboveMax(NumPaths))
    CounterLevelAboveMax = 0
    DO iPath = 1, NumPaths
      ! Determine whether the maximal energy level of vibrational levels is greater than the provided collisional level
      CounterLevelAboveMax(iPath) = COUNT(SpecXSec(iCase)%ReactionPath(iPath)%XSecData(1,:).GT.MaxEnergyColl)
    END DO
    ! Abort in case of EFFECTIVE or resizing the CollXSecData array in case of ELASTIC, if vibrational energy levels are above
    IF(ANY(CounterLevelAboveMax.GT.0)) THEN
      IF(SpecXSec(iCase)%CollXSec_Effective) THEN
        CALL abort(__STAMP__,'ERROR: Effective cross-section has a smaller energy range than the vibrational level cross-section!')
      ELSE
        ! Get the index of vibrational level with the most energy values above
        MaxLevelIndex = MAXLOC(CounterLevelAboveMax,DIM=1)
        NumLevel = MAXVAL(CounterLevelAboveMax)
        ! Allocate a temporary array for the new collisional cross-section data
        ALLOCATE(CollXSecDataTemp(1:2,1:MaxDim+NumLevel))
        ! Store the old array in the new temporary array
        CollXSecDataTemp(1:2,1:MaxDim) = SpecXSec(iCase)%CollXSecData(1:2,1:MaxDim)
        ! Determine the bounds of vibrational energy levels to be added to the collision array
        MaxDimLevel = UBOUND(SpecXSec(iCase)%ReactionPath(MaxLevelIndex)%XSecData,DIM=2)
        ! Add the energy levels but not the cross-section, it will added in the next step
        CollXSecDataTemp(1,MaxDim+1:MaxDim+NumLevel) = SpecXSec(iCase)%ReactionPath(MaxLevelIndex)%XSecData(1,MaxDimLevel-NumLevel+1:MaxDimLevel)
        CollXSecDataTemp(2,MaxDim+1:MaxDim+NumLevel) = 0.
        CALL MOVE_ALLOC(CollXSecDataTemp,SpecXSec(iCase)%CollXSecData)
        LBWRITE(*,*) 'Resizing CollXSecData array from ', MaxDim, ' to ', UBOUND(SpecXSec(iCase)%CollXSecData,dim=2), &
                        'due to reaction #', ChemReac%CollCaseInfo(iCase)%ReactionIndex(MaxLevelIndex)
        MaxDim = UBOUND(SpecXSec(iCase)%CollXSecData,dim=2)
      END IF
    END IF
    DEALLOCATE(CounterLevelAboveMax)
    ! When no effective cross-section is available, the total cross-section has to be determined
    IF(.NOT.SpecXSec(iCase)%CollXSec_Effective) THEN
      MaxDim = SIZE(SpecXSec(iCase)%CollXSecData,2)
      ! Interpolate the reaction cross section at the energy levels of the collision collision cross section
      DO iPath = 1, NumPaths
        DO iStep = 1, MaxDim
          ReactionCrossSection = InterpolateCrossSection(SpecXSec(iCase)%ReactionPath(iPath)%XSecData,SpecXSec(iCase)%CollXSecData(1,iStep))
          SpecXSec(iCase)%CollXSecData(2,iStep) = SpecXSec(iCase)%CollXSecData(2,iStep) + ReactionCrossSection
        END DO
      END DO
    END IF  ! SpecXSec(iCase)%CollXSec_Effective
  END IF    ! SpecXSec(iCase)%UseCollXSec
END DO

! 3.) Recalculate the null collision probability with the new total cross-section
IF(XSec_NullCollision) THEN
  DO iSpec = 1, nSpecies
    TotalProb = 0.
    DO jSpec = iSpec, nSpecies
      iCase = CollInf%Coll_Case(iSpec,jSpec)
      IF(SpecXSec(iCase)%UseCollXSec) THEN
        CALL DetermineNullCollProb(iCase,iSpec,jSpec)
        IF(BGGas%UseDistribution)THEN
          TotalProb = TotalProb + MAXVAL(SpecXSec(iCase)%ProbNullElem)
        ELSE
          TotalProb = TotalProb + SpecXSec(iCase)%ProbNull
        END IF ! BGGas%UseDistribution
        IF(TotalProb.GT.1.0) CALL abort(__STAMP__,&
      'ERROR: Total null collision probability is above unity. Please reduce the time step! Probability is: ',RealInfoOpt=TotalProb)
      END IF
    END DO
  END DO
END IF

END SUBROUTINE MCC_Chemistry_Init


!===================================================================================================================================
!> Initialize Photoionization cross-section (XSec) by reading the data from .h5
!> 1. Check which reactions are 'phIonXSec'
!> 2. Check the educts (photon+X), also switch the ordering of the names.e.g, N2-photon or photon-N2 of the container
!> 3. Load the photon energy spectrum from the container "SPECTRUM" under "N2-photon"
!> 4. Load the chemical reaction cross-sections but do not vary the order of the container name, e.g., H-HIon1-electron
!>    if the name is not correctly ordered it will abort
!> 5. Interpolate the cross-section data to the wavelength spectrum and discard the original data. Keep only the interpolated data.
!> 6. Find the first and last wavelength for which cross-sections are available, ignore wavelengths outside in the simulation
!===================================================================================================================================
SUBROUTINE InitPhotoionizationXSec()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_MCC_Vars  ,ONLY: NbrOfPhotonXsecReactions,SpecPhotonXSec,PhotoReacToReac,NbrOfPhotonXsecLines
USE MOD_MCC_Vars  ,ONLY: SpecPhotonXSecInterpolated,PhotoIonFirstLine,PhotoIonLastLine,PhotonDistribution,ReacToPhotoReac
USE MOD_MCC_Vars  ,ONLY: PhotonSpectrum,PhotonEnergies,MaxPhotonXSec
USE MOD_MCC_XSec  ,ONLY: ReadReacPhotonXSec,ReadReacPhotonSpectrum
USE MOD_DSMC_Vars ,ONLY: SpecDSMC,ChemReac
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=50)     :: EductPair,EductPairOld
INTEGER               :: iReac,iPhotoReac
REAL                  :: dx,dy,TotalEnergyFraction
INTEGER               :: ReadInNumOfReact,dims(2),iLine,location(3)
!===================================================================================================================================

ReadInNumOfReact = ChemReac%NumOfReact
! Set Mapping
ALLOCATE(PhotoReacToReac(NbrOfPhotonXsecReactions))
PhotoReacToReac = -1 ! Initialize
ALLOCATE(ReacToPhotoReac(ReadInNumOfReact))
ReacToPhotoReac = -1 ! Initialize
iPhotoReac = 0
!> 1. Check which reactions are 'phIonXSec'
DO iReac = 1, ReadInNumOfReact
  IF(TRIM(ChemReac%ReactModel(iReac)).EQ.'phIonXSec')THEN
    ! Create mappings
    iPhotoReac                  = iPhotoReac + 1
    PhotoReacToReac(iPhotoReac) = iReac
    ReacToPhotoReac(iReac)      = iPhotoReac
!> 2. Check the educts (photon+X), also switch the ordering of the names.e.g, N2-photon or photon-N2 of the container
    EductPair = TRIM(SpecDSMC(ChemReac%Reactants(iReac,1))%Name)//'-photon'
    IF(iPhotoReac.GT.1)THEN
      IF(TRIM(EductPair).NE.TRIM(EductPairOld)) CALL abort(__STAMP__,'Currently only one photo reaction is implemented')
    END IF ! iPhotoReac.GT.1
    EductPairOld = EductPair
  END IF ! TRIM(ChemReac%ReactModel(iReac)).EQ.'phIonXSec'
END DO ! iReac = 1, ReadInNumOfReact

!> 3. Load the photon energy spectrum from the container "SPECTRUM" under "N2-photon"
CALL ReadReacPhotonSpectrum(1) ! 1 corresponds to iPhotoReac=1 and for all reactions, the Educt currently must be the same

!> 4. Load the chemical reaction cross-sections but do not vary the order of the container name, e.g., H-HIon1-electron
!>    if the name is not correctly ordered it will abort
ALLOCATE(SpecPhotonXSec(NbrOfPhotonXsecReactions))
DO iPhotoReac = 1, NbrOfPhotonXsecReactions
  CALL ReadReacPhotonXSec(iPhotoReac)
END DO ! iPhotoReac = 1, NbrOfPhotonXsecReactions

! Normalize energy fraction
TotalEnergyFraction = SUM(PhotonSpectrum(2,:))
IF(.NOT.ALMOSTEQUALRELATIVE(TotalEnergyFraction,1.0,1e-3)) CALL abort(__STAMP__,'Sum of energy fraction is not 1.0')
PhotonSpectrum(2,:) = PhotonSpectrum(2,:) / TotalEnergyFraction

!> 5. Interpolate the cross-section data to the wavelength spectrum and discard the original data. Keep only the interpolated data.
dims=SHAPE(PhotonSpectrum)
NbrOfPhotonXsecLines = dims(2)
ALLOCATE(SpecPhotonXSecInterpolated(NbrOfPhotonXsecLines,NbrOfPhotonXsecReactions+3))
SpecPhotonXSecInterpolated = 0.
ALLOCATE(PhotonDistribution(NbrOfPhotonXsecLines))
DO iLine = 1, NbrOfPhotonXsecLines
  SpecPhotonXSecInterpolated(iLine,1:2) = PhotonSpectrum(:,iLine)
END DO ! iLine = 1, dims(2)
DEALLOCATE(PhotonSpectrum)

! Interpolate the cross-section data to to photon energy spectrum
DO iPhotoReac = 1, NbrOfPhotonXsecReactions
  ASSOCIATE( MINeV => MINVAL(SpecPhotonXSec(iPhotoReac)%XSecData(1,:)), MAXeV => MAXVAL(SpecPhotonXSec(iPhotoReac)%XSecData(1,:)))
    DO iLine = 1, NbrOfPhotonXsecLines
      ASSOCIATE( energy => SpecPhotonXSecInterpolated(iLine,1) )
        IF((energy.GE.MINeV).AND.(energy.LE.MAXeV))THEN
          ! Get the location of the element in the array with min value
          location(3) = MINLOC(ABS(SpecPhotonXSec(iPhotoReac)%XSecData(1,:)-energy),1)

          ! Get lower index of photon energy in XSec array
          IF(energy.LE.SpecPhotonXSec(iPhotoReac)%XSecData(1,location(3)))THEN
            location(2) = location(3)
            location(1) = location(3)-1
          ELSEIF(energy.GE.SpecPhotonXSec(iPhotoReac)%XSecData(1,location(3)))THEN
            location(2) = location(3)+1
            location(1) = location(3)
          ELSE
            CALL abort(__STAMP__,'Photon spectrum interpolation failed.')
          END IF ! energy.LE.SpecPhotonXSec(iPhotoReac)%XSecData(1,location(3))

          ! Check if last element is exactly matched
          IF(location(1).EQ.dims(2))THEN
            SpecPhotonXSecInterpolated(iLine,2+iPhotoReac) = SpecPhotonXSec(iPhotoReac)%XSecData(location(1),2)
          ELSE
            dx = SpecPhotonXSec(iPhotoReac)%XSecData(1,location(2))-SpecPhotonXSec(iPhotoReac)%XSecData(1,location(1))
            dy = SpecPhotonXSec(iPhotoReac)%XSecData(2,location(2))-SpecPhotonXSec(iPhotoReac)%XSecData(2,location(1))
            IF(ABS(dx).GT.0.)THEN
              ! linear interpolation
              SpecPhotonXSecInterpolated(iLine,2+iPhotoReac) = SpecPhotonXSec(iPhotoReac)%XSecData(2,location(1)) &
                  + dy/dx * (energy-SpecPhotonXSec(iPhotoReac)%XSecData(1,location(1)))
            ELSE
              SpecPhotonXSecInterpolated(iLine,2+iPhotoReac) = SpecPhotonXSec(iPhotoReac)%XSecData(2,location(1))
            END IF ! dy.G
          END IF ! location.EQ.dims(2)

          ! Calculate total cross section
          SpecPhotonXSecInterpolated(iLine,NbrOfPhotonXsecReactions+3) = &
              SpecPhotonXSecInterpolated(iLine,NbrOfPhotonXsecReactions+3) + SpecPhotonXSecInterpolated(iLine,2+iPhotoReac)
        END IF ! (SpecPhotonXSecInterpolated(iLine,1).GE.MINeV).AND.()
      END ASSOCIATE
    END DO ! iLine = 1, NbrOfPhotonXsecLines
  END ASSOCIATE
END DO ! iPhotoReac = 1, NbrOfPhotonXsecReactions
DEALLOCATE(SpecPhotonXSec)

!> 6. Find the first and last wavelength for which cross-sections are available, ignore wavelengths outside in the simulation
PhotoIonFirstLine = HUGE(1)
PhotoIonLastLine = 0
DO iLine = 1, NbrOfPhotonXsecLines
  IF(SpecPhotonXSecInterpolated(iLine,NbrOfPhotonXsecReactions+3).GT.0.)THEN
    PhotoIonFirstLine = MIN(iLine,PhotoIonFirstLine)
    PhotoIonLastLine  = MAX(iLine,PhotoIonLastLine)
  END IF ! SpecPhotonXSecInterpolated(iLine,NbrOfPhotonXsecReactions+3).GT.0.
END DO ! iLine = 1, NbrOfPhotonXsecLines
IF(PhotoIonLastLine.LT.PhotoIonFirstLine) CALL abort(__STAMP__,'Photoionization XSec read-in failed. No lines interpolated.')
MaxPhotonXSec = MAXVAL(SpecPhotonXSecInterpolated(PhotoIonFirstLine:PhotoIonLastLine,2))
ALLOCATE(PhotonEnergies(PhotoIonFirstLine:PhotoIonLastLine,1+NbrOfPhotonXsecReactions))

END SUBROUTINE InitPhotoionizationXSec


!===================================================================================================================================
!> 1. Check whether the photon energy is sufficient to trigger the chemical reaction
!> 2. Check for intermediate lines with zero cross sections (not allowed)
!> 3. Determine the lost energy (wavelengths outside of the range of the cross-section data)
!===================================================================================================================================
SUBROUTINE CheckPhotoionizationXSec()
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars  ,ONLY: eV2Joule
USE MOD_MCC_Vars      ,ONLY: NbrOfPhotonXsecReactions,NbrOfPhotonXsecLines,SpecPhotonXSecInterpolated,PhotoIonFirstLine,PhotoIonLastLine
USE MOD_MCC_Vars      ,ONLY: ReacToPhotoReac
USE MOD_DSMC_Vars     ,ONLY: ChemReac
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: ZeroLine,iLine,iReac,iPhotoReac
REAL                  :: LostEnergy(2)
!===================================================================================================================================
IF(.NOT.MPIRoot) RETURN

! 1. Check whether the photon energy is sufficient to trigger the chemical reaction
DO iReac = 1, ChemReac%NumOfReact
  IF(TRIM(ChemReac%ReactModel(iReac)).NE.'phIonXSec') CYCLE
  DO iLine = PhotoIonFirstLine, PhotoIonLastLine
    ! Add photon energy to formation energy (which is negative and should be positive after the addition)
    IF(ChemReac%EForm(iReac)+SpecPhotonXSecInterpolated(iLine,1)*eV2Joule.LE.0.0)THEN
      iPhotoReac = ReacToPhotoReac(iReac)
      ! Check if this photo reaction has a cross-section unequal zero
      ! If it is zero, this reaction would never happen anyway, therefore ignore it
      IF(SpecPhotonXSecInterpolated(iLine,2+iPhotoReac).LE.0.) THEN
        CYCLE
      ELSE
        SWRITE (*,*) iLine,SpecPhotonXSecInterpolated(iLine,1:2+iPhotoReac-1),"      requires fix       ",&
        SpecPhotonXSecInterpolated(iLine,2+iPhotoReac+1:)
      END IF
      ! Abort if photon energy is not sufficient and the photo reaction has a cross-section greater zero
      SWRITE (UNIT_stdOut,*) "      iLine      energy-fraction                              %             iPhotoReac(1)        ........"
      SWRITE (UNIT_stdOut,*) iLine,SpecPhotonXSecInterpolated(iLine,:)
      SWRITE (UNIT_stdOut,*) "-----------------------------------------"
      SWRITE (UNIT_stdOut,*) "iLine                                          = ", iLine
      SWRITE (UNIT_stdOut,*) "iReac                                          = ", iReac
      SWRITE (UNIT_stdOut,*) "iPhotoReac                                     = ", iPhotoReac
      SWRITE (UNIT_stdOut,*) "ChemReac%EForm(iReac)                          = ", ChemReac%EForm(iReac)
      SWRITE (UNIT_stdOut,*) "PhotonEnergy [J]                               = ", SpecPhotonXSecInterpolated(iLine,1)*eV2Joule
      SWRITE (UNIT_stdOut,*) "PhotonEnergy [eV]                              = ", SpecPhotonXSecInterpolated(iLine,1)
      SWRITE (UNIT_stdOut,*) "SpecPhotonXSecInterpolated(iLine,2+iPhotoReac) = ", SpecPhotonXSecInterpolated(iLine,2+iPhotoReac)
      CALL abort(__STAMP__,&
        'Photoionization not possible because the photon energy is too low for this reaction. This is not considered yet.')
    END IF
  END DO ! iLine = , PhotoIonLastLine
END DO

! 2. Sanity check: intermediate lines with zero cross sections are not allowed
ZeroLine = 0
DO iLine = PhotoIonFirstLine, PhotoIonLastLine
  IF(SpecPhotonXSecInterpolated(iLine,NbrOfPhotonXsecReactions+3).LE.0.)THEN
    ZeroLine = iLine
    EXIT
  END IF !
END DO ! iLine = PhotoIonFirstLine, PhotoIonLastLine

IF(ZeroLine.GT.0)THEN
  WRITE(UNIT_StdOut,*) ""
  WRITE(UNIT_StdOut,*) "ERROR: Found intermediate wavelength with no correspoding chemical reaction (cross section > 0)"
  WRITE(UNIT_StdOut,*) "SpecPhotonXSecInterpolated: Line ",ZeroLine," is broken"
  DO iLine = PhotoIonFirstLine, PhotoIonLastLine
    IF(iLine.EQ.ZeroLine) WRITE (*,*) " -------------------------------- HERE -------------------------------- "
    SWRITE (*,*) iLine,SpecPhotonXSecInterpolated(iLine,:)
    IF(iLine.EQ.ZeroLine) WRITE (*,*) " -------------------------------- HERE -------------------------------- "
  END DO ! iLine = 1, dims(2)
  CALL abort(__STAMP__,'Found intermediate wavelength with no correspoding chemical reaction (any cross section > 0)')
END IF ! ZeroLine.GT.0

! 3. Calculate the lost energy (wavelengths for which no reactions are possible)
LostEnergy(:) = 0.
IF(PhotoIonFirstLine.GT.1)THEN
  DO iLine = 1, PhotoIonFirstLine-1
    LostEnergy(1) = LostEnergy(1) + SpecPhotonXSecInterpolated(iLine,2)
  END DO ! iLine = 1, PhotoIonFirstLine-1
END IF ! PhotoIonFirstLine.GT.1
IF(PhotoIonLastLine.LT.NbrOfPhotonXsecLines)THEN
  DO iLine = PhotoIonLastLine+1, NbrOfPhotonXsecLines
    LostEnergy(2) = LostEnergy(2) + SpecPhotonXSecInterpolated(iLine,2)
  END DO ! iLine = PhotoIonLastLine+1, NbrOfPhotonXsecLines
END IF ! PhotoIonLastLine.LT.NbrOfPhotonXsecLines
IF(SUM(LostEnergy).GT.0.)THEN
  WRITE (UNIT_stdOut,'(A,3(F6.2,A))') "Lost energy content in photoionization XSec: ", &
      LostEnergy(1)*100., "% (high photon wavelength),",&
      LostEnergy(2)*100., "% (low photon wavelength),",&
      SUM(LostEnergy)*100., "% (total)"
  WRITE (UNIT_stdOut,'(A,I0,A,I0,A,I0,A)') "Only considering level ",PhotoIonFirstLine," to ",PhotoIonLastLine, " from 1 to ",&
      NbrOfPhotonXsecLines," due to the cut-off"
END IF ! SUM(LostEnergy).GT.0.

END SUBROUTINE CheckPhotoionizationXSec


SUBROUTINE FinalizeMCC()
!----------------------------------------------------------------------------------------------------------------------------------!
! Finalize MCC/XSec variables
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_MCC_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(SpecXSec)
SDEALLOCATE(SpecPhotonXSecInterpolated)
SDEALLOCATE(PhotonDistribution)
SDEALLOCATE(PhotonEnergies)
SDEALLOCATE(PhotoReacToReac)
SDEALLOCATE(ReacToPhotoReac)
END SUBROUTINE FinalizeMCC

END MODULE MOD_MCC_Init
