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

MODULE MOD_Radiation_Init
!===================================================================================================================================
! Initialization of DSMC
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE InitRadiation
  MODULE PROCEDURE InitRadiation
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: InitRadiation, FinalizeRadiation, DefineParametersRadiation
!===================================================================================================================================

CONTAINS




SUBROUTINE DefineParametersRadiation()
!==================================================================================================================================
! Define parameters for Radiation
!==================================================================================================================================
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Radiation")

CALL prms%CreateRealOption(   'Part-Species[$]-RadiationTtrans',       'Translational temperature, K', '0.0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(   'Part-Species[$]-RadiationTelec',        'Electronic excitation temperature, K', '0.0', &
                                                                       numberedmulti=.TRUE.)
CALL prms%CreateRealOption(   'Part-Species[$]-RadiationNumDens',      'Number density, 1/cm3', '0.0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(   'Part-Species[$]-RadiationTvib',         'Vibrational temperature, K', '0.0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(   'Part-Species[$]-RadiationTrot',         'Rotational temperature, K', '0.0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(   'Part-Species[$]-RadiationIonizationEn', 'Ionization Energy, 1/cm', '0.0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(   'Part-Species[$]-RadiationRadius_A',     'Species radius, A', '0.0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(   'Part-Species[$]-Starkex',               'Exponent for the determination of Stark broadening', &
                                                                       '0.0', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(    'Part-Species[$]-NuclCharge',            'Nuclear charge:\n'//&
                                                                       '1: neutral atom\n'//&
                                                                       '2: singly ionized atom', '1', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption('Part-Species[$]-DoRadiation',           'Considering species for radiative emission', '.TRUE.', &
                                                                       numberedmulti=.TRUE.)
CALL prms%CreateRealOption(   'Radiation-MinWaveLen',                  'Lower wavelength limit for radiation calculation', '100.0')
CALL prms%CreateRealOption(   'Radiation-MaxWaveLen',                  'Upper wavelength limit for radiation calculation','1000.0')
CALL prms%CreateIntOption(    'Radiation-WaveLenDiscr',                'Number of discretization points', '10000')
CALL prms%CreateIntOption(    'Radiation-WaveLenReductionFactor',      'Number of discretization points', '1')
CALL prms%CreateIntOption(    'Radiation-WaveLenReductionFactorOutput',      'Number of discretization points', '1')
CALL prms%CreateIntOption(    'Radiation-RadType',                     'Select radiation type:\n'//&
                                                                       '1: particle radiation\n'//&
                                                                       '2: black body radiation\n'//&
                                                                       '3: radiation solver only', '3')
CALL prms%CreateLogicalOption('Radiation-ff',                          'Enable free-free radiation', '.FALSE.')
CALL prms%CreateLogicalOption('Radiation-bf',                          'Enable bound-free radiation (only atomic)', '.FALSE.')
CALL prms%CreateLogicalOption('Radiation-bb-atoms',                    'Enable atomic bound-bound radiation', '.FALSE.')
CALL prms%CreateLogicalOption('Radiation-bb-molecules',                'Enable molecular bound-bound radiation', '.FALSE.')
CALL prms%CreateLogicalOption('Radiation-MacroRadInput',               'Reading in flow field data as radiation input', '.FALSE.')
CALL prms%CreateLogicalOption('Radiation-MacroInput-SortCellsY',       'Sorts Cells in y-direction', '.FALSE.')
CALL prms%CreateLogicalOption('Radiation-UseElectronicExcitation',     'Use el. excitation to populate upper state densitites', '.TRUE.')
CALL prms%CreateRealOption(   'Radiation-NumDensElectrons',            'Electron number density, 1/cm3', '0.0')
CALL prms%CreateRealOption(   'Radiation-TElectrons',                  'Electron temperature, K', '0.0')
CALL prms%CreateStringOption( 'Radiation-Species[$]-SpectraFileName',  'File name of data file', 'none', numberedmulti=.TRUE.)
CALL prms%CreateStringOption( 'Radiation-MacroInput-Filename', &
                              'TO-DO')
END SUBROUTINE DefineParametersRadiation




SUBROUTINE InitRadiation()
!===================================================================================================================================
! Init of DSMC Vars
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,          ONLY : PlanckConst, c
USE MOD_Mesh_Vars,             ONLY : nGlobalElems
USE MOD_Particle_Mesh_Vars,    ONLY : nComputeNodeElems
USE MOD_ReadInTools
USE MOD_PARTICLE_Vars,         ONLY : nSpecies
USE MOD_Radiation_Vars
USE MOD_DSMC_Vars,             ONLY : SpecDSMC
USE MOD_Radiation_ReadIn,      ONLY : Radiation_readin_atoms, Radiation_readin_molecules
USE MOD_Mesh_Tools,            ONLY : GetGlobalElemID
#if USE_MPI
!USE MOD_MPI_Shared_Vars
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#else
USE MOD_Mesh_Vars,             ONLY : nElems
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  CHARACTER(32)         :: hilf
  INTEGER               :: iSpec, iWaveLen, firstElem, lastElem, iElem
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)') ' INIT RADIATION SOLVER...'

RadiationSwitches%RadType       = GETINT('Radiation-RadType',           '3')
ALLOCATE(RadiationInput(nSpecies))
ALLOCATE(SpeciesRadiation(nSpecies))
SpeciesRadiation(:)%nLevels = 0
SpeciesRadiation(:)%nLines = 0

IF (RadiationSwitches%RadType.NE.2) THEN
  DO iSpec = 1, nSpecies
    IF(SpecDSMC(iSpec)%InterID.EQ.4) CYCLE
    WRITE(UNIT=hilf,FMT='(I0)') iSpec
    RadiationInput(iSpec)%Ttrans(4) = GETREAL('Part-Species'//TRIM(hilf)//'-RadiationTtrans')
    RadiationInput(iSpec)%Telec = GETREAL('Part-Species'//TRIM(hilf)//'-RadiationTelec')
    RadiationInput(iSpec)%NumDens = GETREAL('Part-Species'//TRIM(hilf)//'-RadiationNumDens')
    IF((SpecDSMC(iSpec)%InterID.EQ.2) .OR. (SpecDSMC(iSpec)%InterID.EQ.20)) THEN
      RadiationInput(iSpec)%Tvib = GETREAL('Part-Species'//TRIM(hilf)//'-RadiationTvib')
      RadiationInput(iSpec)%Trot = GETREAL('Part-Species'//TRIM(hilf)//'-RadiationTrot')
    END IF
    RadiationInput(iSpec)%IonizationEn = GETREAL('Part-Species'//TRIM(hilf)//'-RadiationIonizationEn')
    RadiationInput(iSpec)%IonizationEn = RadiationInput(iSpec)%IonizationEn *PlanckConst*c*100.

    RadiationInput(iSpec)%DoRadiation = GETLOGICAL('Part-Species'//TRIM(hilf)//'-DoRadiation')

    IF((SpecDSMC(iSpec)%InterID .EQ. 1) .OR. (SpecDSMC(iSpec)%InterID .EQ. 10)) THEN !Only for atoms (1) and atomic ions (10)
      CALL Radiation_readin_atoms(iSpec)
      RadiationInput(iSpec)%Radius = GETREAL('Part-Species'//TRIM(hilf)//'-RadiationRadius_A')
      RadiationInput(iSpec)%Radius = RadiationInput(iSpec)%Radius *1.0E-10
      RadiationInput(iSpec)%Starkex = GETREAL('Part-Species'//TRIM(hilf)//'-Starkex')
      RadiationInput(iSpec)%NuclCharge = GETINT('Part-Species'//TRIM(hilf)//'-NuclCharge')
    END IF

    IF((SpecDSMC(iSpec)%InterID .EQ. 2) .OR. (SpecDSMC(iSpec)%InterID .EQ. 20)) THEN !Only for molecules (2) and molecular ions (20)
      CALL Radiation_readin_molecules(iSpec)
    END IF
  END DO
END IF

RadiationParameter%MinWaveLen   = GETREAL('Radiation-MinWaveLen')
RadiationParameter%MinWaveLen   = RadiationParameter%MinWaveLen*1.E-9
RadiationParameter%MaxWaveLen   = GETREAL('Radiation-MaxWaveLen')
RadiationParameter%MaxWaveLen   = RadiationParameter%MaxWaveLen*1.E-9
RadiationParameter%WaveLenDiscr = GETINT('Radiation-WaveLenDiscr')
RadiationParameter%WaveLenReductionFactor = GETINT('Radiation-WaveLenReductionFactor')
RadiationParameter%WaveLenReductionFactorOutput = GETINT('Radiation-WaveLenReductionFactorOutput')
IF((RadiationSwitches%RadType.EQ.3) .AND. (nGlobalElems.EQ.1)) RadiationParameter%WaveLenReductionFactor = 1
IF(RadiationSwitches%RadType.EQ.2) RadiationParameter%WaveLenReductionFactor = 1
IF (RadiationParameter%WaveLenReductionFactor.NE.1) THEN
  RadiationParameter%WaveLenDiscrCoarse = NINT(REAL(RadiationParameter%WaveLenDiscr)/ REAL(RadiationParameter%WaveLenReductionFactor))
  RadiationParameter%WaveLenReductionFactor = INT(RadiationParameter%WaveLenDiscr/RadiationParameter%WaveLenDiscrCoarse)
  SWRITE(UNIT_stdOut,'(A)') 'Corrected WaveLenReductionFactor is ', RadiationParameter%WaveLenReductionFactor
  ALLOCATE(RadiationParameter%WaveLenCoarse(RadiationParameter%WaveLenDiscrCoarse))
  RadiationParameter%WaveLenCoarse = 0.0
ELSE
  RadiationParameter%WaveLenDiscrCoarse = RadiationParameter%WaveLenDiscr
END IF
IF (RadiationParameter%WaveLenReductionFactorOutput.GT.1) THEN
  RadiationParameter%WaveLenDiscrOutput = NINT(REAL(RadiationParameter%WaveLenDiscr)/ REAL(RadiationParameter%WaveLenReductionFactorOutput))
  RadiationParameter%WaveLenReductionFactorOutput = INT(RadiationParameter%WaveLenDiscr/RadiationParameter%WaveLenDiscrOutput)
  SWRITE(UNIT_stdOut,'(A)') 'Corrected WaveLenReductionFactorOutput is ', RadiationParameter%WaveLenReductionFactorOutput
ELSE
  RadiationParameter%WaveLenDiscrOutput = RadiationParameter%WaveLenDiscr
END IF
IF(RadiationParameter%MinWaveLen.GE.RadiationParameter%MaxWaveLen) THEN
  CALL abort(&
                __STAMP__&
                ,' ERROR: Radiation - maximum wavelength is smaller than minimum wavelength')
END IF
IF(RadiationParameter%WaveLenDiscr.LT.100) THEN
  CALL abort(&
                __STAMP__&
                ,' ERROR: Radiation - number of wavelength discretization points is to small')
END IF
RadiationParameter%WaveLenIncr  = (RadiationParameter%MaxWaveLen - RadiationParameter%MinWaveLen) &
  / (RadiationParameter%WaveLenDiscr-1)

ALLOCATE(RadiationParameter%WaveLen(RadiationParameter%WaveLenDiscr))
RadiationParameter%WaveLen = 0.0

DO iWaveLen = 1, RadiationParameter%WaveLenDiscr
  RadiationParameter%WaveLen(iWaveLen) = RadiationParameter%MinWaveLen + (iWaveLen-1) * RadiationParameter%WaveLenIncr
END DO
IF (RadiationParameter%WaveLenReductionFactor.NE.1) THEN
  RadiationParameter%WaveLenIncrCoarse = (RadiationParameter%MaxWaveLen - RadiationParameter%MinWaveLen) &
  / (RadiationParameter%WaveLenDiscr*RadiationParameter%WaveLenReductionFactor-1)
  DO iWaveLen = 1, NINT(RadiationParameter%WaveLenIncrCoarse)
    RadiationParameter%WaveLenCoarse(iWaveLen) = RadiationParameter%MinWaveLen + (iWaveLen-1) * RadiationParameter%WaveLenIncrCoarse
  END DO
END IF

RadiationSwitches%ff                      = GETLOGICAL('Radiation-ff')
RadiationSwitches%bf                      = GETLOGICAL('Radiation-bf')
RadiationSwitches%bb_at                   = GETLOGICAL('Radiation-bb-atoms')
RadiationSwitches%bb_mol                  = GETLOGICAL('Radiation-bb-molecules')
RadiationSwitches%MacroRadInput           = GETLOGICAL('Radiation-MacroRadInput')
RadiationSwitches%SortCellsY              = GETLOGICAL('Radiation-MacroInput-SortCellsY')
RadiationSwitches%UseElectronicExcitation = GETLOGICAL('Radiation-UseElectronicExcitation')

IF (RadiationSwitches%MacroRadInput) CALL MacroscopicRadiationInput()

NumDensElectrons = GETREAL('Radiation-NumDensElectrons')
TElectrons       = GETREAL('Radiation-TElectrons')

#if USE_MPI
  ! allocate shared array for Radiation_Emission/Absorption_Spec
CALL Allocate_Shared((/RadiationParameter%WaveLenDiscrCoarse,nComputeNodeElems/), Radiation_Emission_Spec_Shared_Win,Radiation_Emission_Spec_Shared)
CALL MPI_WIN_LOCK_ALL(0,Radiation_Emission_Spec_Shared_Win,IERROR)
CALL Allocate_Shared((/INT(RadiationParameter%WaveLenDiscrCoarse,IK)*INT(nGlobalElems,IK)/),Radiation_Absorption_Spec_Shared_Win,Radiation_Absorption_Spec_Shared)
CALL MPI_WIN_LOCK_ALL(0,Radiation_Absorption_Spec_Shared_Win,IERROR)
CALL Allocate_Shared((/INT(RadiationParameter%WaveLenDiscrCoarse,IK)*INT(nGlobalElems,IK)*INT(nSpecies,IK)/),Radiation_Absorption_SpecPercent_Shared_Win,Radiation_Absorption_SpecPercent_Shared)
CALL MPI_WIN_LOCK_ALL(0,Radiation_Absorption_SpecPercent_Shared_Win,IERROR)
CALL Allocate_Shared((/nSpecies,nComputeNodeElems,2/), Radiation_ElemEnergy_Species_Shared_Win,Radiation_ElemEnergy_Species_Shared)
CALL MPI_WIN_LOCK_ALL(0,Radiation_ElemEnergy_Species_Shared_Win,IERROR)

Radiation_Emission_spec => Radiation_Emission_spec_Shared
Radiation_Absorption_Spec(1:RadiationParameter%WaveLenDiscrCoarse ,1:nGlobalElems) => Radiation_Absorption_Spec_Shared
Radiation_Absorption_SpecPercent(1:RadiationParameter%WaveLenDiscrCoarse ,1:nSpecies, 1:nGlobalElems) => Radiation_Absorption_SpecPercent_Shared
Radiation_ElemEnergy_Species => Radiation_ElemEnergy_Species_Shared
#else
! allocate local array for ElemInfo
ALLOCATE(Radiation_Emission_spec(RadiationParameter%WaveLenDiscrCoarse,nElems))
ALLOCATE(Radiation_Absorption_Spec(RadiationParameter%WaveLenDiscrCoarse,nElems))
ALLOCATE(Radiation_Absorption_SpecPercent(RadiationParameter%WaveLenDiscrCoarse,nSpecies,nElems))
ALLOCATE(Radiation_ElemEnergy_Species(nSpecies,nElems,2))
#endif  /*USE_MPI*/

ALLOCATE(Radiation_Absorption_SpeciesWave(RadiationParameter%WaveLenDiscrCoarse,nSpecies))

#if USE_MPI
  firstElem = INT(REAL( myComputeNodeRank   *nComputeNodeElems)/REAL(nComputeNodeProcessors))+1
  lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeElems)/REAL(nComputeNodeProcessors))
#else
  firstElem = 1
  lastElem  = nElems
#endif

DO iElem = firstElem, lastElem
  Radiation_Emission_spec(:,iElem) = 0.0
  Radiation_Absorption_Spec(:,GetGlobalElemID(iElem)) = 0.0
  Radiation_Absorption_SpecPercent(:,:,GetGlobalElemID(iElem)) = 0
  Radiation_ElemEnergy_Species(:,iElem,:) =0.0
END DO
#if USE_MPI
  CALL BARRIER_AND_SYNC(Radiation_Emission_Spec_Shared_Win ,MPI_COMM_SHARED)
  CALL BARRIER_AND_SYNC(Radiation_ElemEnergy_Species_Shared_Win ,MPI_COMM_SHARED)
  CALL BARRIER_AND_SYNC(Radiation_Absorption_Spec_Shared_Win ,MPI_COMM_SHARED)
  CALL BARRIER_AND_SYNC(Radiation_Absorption_SpecPercent_Shared_Win ,MPI_COMM_SHARED)
  IF(nLeaderGroupProcs.GT.1)THEN
    IF(myComputeNodeRank.EQ.0)THEN
      CALL MPI_ALLGATHERV( MPI_IN_PLACE                                     &
                     , 0                                                    &
                     , MPI_DATATYPE_NULL                                    &
                     , Radiation_Absorption_Spec                            &
                     , RadiationParameter%WaveLenDiscrCoarse *recvcountElem &
                     , RadiationParameter%WaveLenDiscrCoarse *displsElem    &
                     , MPI_DOUBLE_PRECISION                                 &
                     , MPI_COMM_LEADERS_SHARED                              &
                     , IERROR)
      CALL MPI_ALLGATHERV( MPI_IN_PLACE                                              &
                     , 0                                                             &
                     , MPI_DATATYPE_NULL                                             &
                     , Radiation_Absorption_SpecPercent                              &
                     , RadiationParameter%WaveLenDiscrCoarse*nSpecies *recvcountElem &
                     , RadiationParameter%WaveLenDiscrCoarse*nSpecies *displsElem    &
                     , MPI_INTEGER2                                                  &
                     , MPI_COMM_LEADERS_SHARED                                       &
                     , IERROR)
    END IF
  END IF
  CALL BARRIER_AND_SYNC(Radiation_Absorption_Spec_Shared_Win ,MPI_COMM_SHARED)
  CALL BARRIER_AND_SYNC(Radiation_Absorption_SpecPercent_Shared_Win ,MPI_COMM_SHARED)
#endif



SWRITE(UNIT_stdOut,'(A)') ' INIT RADIATION SOLVER DONE!'

END SUBROUTINE InitRadiation


SUBROUTINE MacroscopicRadiationInput()
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_io_hdf5
USE MOD_HDF5_Input                ,ONLY: OpenDataFile,CloseDataFile,DatasetExists,ReadArray,GetDataProps
USE MOD_Mesh_Vars                 ,ONLY: offsetElem, nElems
USE MOD_Particle_Vars             ,ONLY: nSpecies
USE MOD_DSMC_Vars                 ,ONLY: SpecDSMC
USE MOD_Radiation_Vars            ,ONLY: RadiationSwitches, MacroRadInputParameters
USE MOD_Mesh_Tools                ,ONLY: GetCNElemID
USE MOD_ReadInTools
USE MOD_Particle_Mesh_Vars        ,ONLY: nComputeNodeElems
#if USE_MPI
USE MOD_Radiation_Vars            ,ONLY: MacroRadInputParameters_Shared,MacroRadInputParameters_Shared_Win
!USE MOD_MPI_Shared_Vars
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#endif
USE MOD_Particle_Mesh_Vars        ,ONLY: ElemMidPoint_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: nVar_HDF5, N_HDF5, nElems_HDF5, iVar, iSpec, iElem, CNElemID, IndexElectronTemp, iSpecElectrons
REAL, ALLOCATABLE                 :: ElemData_HDF5(:,:)
CHARACTER(LEN=300)                :: MacroRadiationInputFile
INTEGER, ALLOCATABLE              :: SortElemInd(:)
REAL, ALLOCATABLE                 :: SortElemYPos(:)
!===================================================================================================================================

MacroRadiationInputFile = GETSTR('Radiation-MacroInput-Filename')
CALL OpenDataFile(MacroRadiationInputFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)

CALL GetDataProps('ElemData',nVar_HDF5,N_HDF5,nElems_HDF5)

#if USE_MPI
! allocate shared array for Radiation_Emission/Absorption_Spec
CALL Allocate_Shared((/nComputeNodeElems,nSpecies,5/), MacroRadInputParameters_Shared_Win,MacroRadInputParameters_Shared)
CALL MPI_WIN_LOCK_ALL(0,MacroRadInputParameters_Shared_Win,IERROR)

MacroRadInputParameters => MacroRadInputParameters_Shared
#else
! allocate local array for ElemInfo
ALLOCATE(MacroRadInputParameters(1:nElems,1:nSpecies,1:5))
#endif  /*USE_MPI*/

ALLOCATE(ElemData_HDF5(1:nVar_HDF5,1:nElems))
! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
  nVar_HDF5  => INT(nVar_HDF5,IK) ,&
  offsetElem => INT(offsetElem,IK),&
  nElems     => INT(nElems,IK)    )
  CALL ReadArray('ElemData',2,(/nVar_HDF5,nElems/),offsetElem,2,RealArray=ElemData_HDF5(:,:))
END ASSOCIATE


IF(RadiationSwitches%SortCellsY) THEN !Sort cells if manually created input is used
  IF(nProcessors.GT.1) THEN
    CALL abort(&
                __STAMP__&
                ,' ERROR: Radiation - Sort cells in y-direction only possible on one processor!')
  END IF
  ALLOCATE(SortElemInd(nElems), SortElemYPos(nElems))
  DO iElem = 1, nElems
    SortElemInd(iElem) = iElem
  END DO
  SortElemYPos(:) = -ElemMidPoint_Shared(2,:)
  CALL BubbleSortID(SortElemYPos, SortElemInd, nElems)

  iVar = 1
  DO iSpec = 1, nSpecies
    DO iElem = 1, nElems
      CNElemID = GetCNElemID(iElem+offsetElem)
      MacroRadInputParameters(SortElemInd(CNElemID),iSpec,1) = MAX(0.,ElemData_HDF5(iVar+ 6,iElem)) !density
      MacroRadInputParameters(SortElemInd(CNElemID),iSpec,2) = MAX(0.,ElemData_HDF5(iVar+ 7,iElem)) !T_vib
      MacroRadInputParameters(SortElemInd(CNElemID),iSpec,3) = MAX(0.,ElemData_HDF5(iVar+ 8,iElem)) !T_rot
      MacroRadInputParameters(SortElemInd(CNElemID),iSpec,4) = MAX(0.,ElemData_HDF5(iVar+ 9,iElem)) !T_elec
      MacroRadInputParameters(SortElemInd(CNElemID),iSpec,5) = MAX(0.,ElemData_HDF5(iVar+11,iElem)) !T_mean
    END DO
    iVar = iVar + DSMC_NVARS
  END DO
ELSE
  iVar = 1
  DO iSpec = 1, nSpecies
    DO iElem = 1, nElems
      CNElemID = GetCNElemID(iElem+offsetElem)
      MacroRadInputParameters(CNElemID,iSpec,1) = MAX(0.,ElemData_HDF5(iVar+ 6,iElem)) !density
      MacroRadInputParameters(CNElemID,iSpec,2) = MAX(0.,ElemData_HDF5(iVar+ 7,iElem)) !T_vib
      MacroRadInputParameters(CNElemID,iSpec,3) = MAX(0.,ElemData_HDF5(iVar+ 8,iElem)) !T_rot
      MacroRadInputParameters(CNElemID,iSpec,4) = MAX(0.,ElemData_HDF5(iVar+ 9,iElem)) !T_elec
      !IF((iSpec.EQ.12) .OR. (iSpec.EQ.13)) MacroRadInputParameters(CNElemID,iSpec,4)=MacroRadInputParameters(CNElemID,iSpec,4)*1.1 !Fe Fe+ +-10percent
      MacroRadInputParameters(CNElemID,iSpec,5) = MAX(0.,ElemData_HDF5(iVar+11,iElem)) !T_mean
    END DO
    iVar = iVar + DSMC_NVARS
  END DO

  IF(.NOT.RadiationSwitches%UseElectronicExcitation) THEN
    iSpecElectrons = 0
    DO iSpec = 1, nSpecies
      IF (SpecDSMC(iSpec)%InterID .EQ. 4) iSpecElectrons = iSpec
    END DO
    IF (iSpecElectrons .EQ. 0) THEN
      PRINT*,  "unknown species number for electrons while reading flow field data"
      STOP
    END IF
    IndexElectronTemp = (iSpecElectrons-1)*DSMC_NVARS+1 + 11 !132 for 11th Species
    DO iElem = 1, nElems
      DO iSpec = 1, nSpecies
        CNElemID = GetCNElemID(iElem+offsetElem)
        IF((SpecDSMC(iSpec)%InterID .EQ. 1) .OR. (SpecDSMC(iSpec)%InterID .EQ. 10) .OR. &
        (SpecDSMC(iSpec)%InterID .EQ. 2) .OR. (SpecDSMC(iSpec)%InterID .EQ. 20)) THEN
          MacroRadInputParameters(CNElemID,iSpec,4) = MAX(0.,ElemData_HDF5(IndexElectronTemp,iElem))
        ELSE IF(SpecDSMC(iSpec)%InterID .EQ. 4) THEN
          CYCLE
        ELSE
          PRINT*, "excitation temperature cannot be matched, unknown InterID for species", iSpec
        END IF
      END DO
    END DO
  END IF
END IF

#if USE_MPI
CALL BARRIER_AND_SYNC(MacroRadInputParameters_Shared_Win ,MPI_COMM_SHARED)
#endif

DEALLOCATE(ElemData_HDF5)

END SUBROUTINE MacroscopicRadiationInput


SUBROUTINE BubbleSortID(a,id,len)
!===================================================================================================================================
! bubble sort, taken from rosetta-wiki and modified for own use
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: len
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                :: a(len)
INTEGER,INTENT(INOUT),OPTIONAL    :: id(len)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: temp
INTEGER                           :: iloop,jloop, temp2
LOGICAL                           :: swapped = .TRUE.
!===================================================================================================================================

IF(PRESENT(id))THEN
  DO jloop=len-1,1,-1
    swapped = .FALSE.
    DO iloop=1,jloop
      IF (a(iloop).GT.a(iloop+1))THEN
        ! switch entries
        temp=a(iloop)
        a(iloop) = a(iloop+1)
        a(iloop+1) = temp
        ! switch ids
        temp2=id(iloop)
        id(iloop) = id(iloop+1)
        id(iloop+1) = temp2
        swapped = .TRUE.
      END IF
    END DO ! iloop
    IF (.NOT. swapped) EXIT
  END DO ! jloop
ELSE
  DO jloop=len-1,1,-1
    swapped = .FALSE.
    DO iloop=1,jloop
      IF (a(iloop).GT.a(iloop+1))THEN
        ! switch entries
        temp=a(iloop)
        a(iloop) = a(iloop+1)
        a(iloop+1) = temp
        swapped = .TRUE.
      END IF
    END DO ! iloop
    IF (.NOT. swapped) EXIT
  END DO ! jloop
END IF
END SUBROUTINE BubbleSortID


SUBROUTINE FinalizeRadiation()
!===================================================================================================================================
!> Deallocating radiation variables
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Radiation_Vars
#if USE_MPI
!USE MOD_MPI_Shared_Vars
USE MOD_MPI_Shared_Vars     ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE(RadiationInput)
SDEALLOCATE(SpeciesRadiation)
SDEALLOCATE(RadiationParameter%WaveLen)

#if USE_MPI
! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
CALL UNLOCK_AND_FREE(Radiation_Emission_Spec_Shared_Win)
CALL UNLOCK_AND_FREE(Radiation_Absorption_Spec_Shared_Win)
CALL UNLOCK_AND_FREE(Radiation_ElemEnergy_Species_Shared_Win)
CALL UNLOCK_AND_FREE(Radiation_Absorption_SpecPercent_Shared_Win)
IF(RadiationSwitches%MacroRadInput) CALL UNLOCK_AND_FREE(MacroRadInputParameters_Shared_Win)
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
#endif /*USE_MPI*/
ADEALLOCATE(MacroRadInputParameters_Shared)
ADEALLOCATE(Radiation_Emission_Spec)
ADEALLOCATE(Radiation_Absorption_Spec)
ADEALLOCATE(Radiation_ElemEnergy_Species)
ADEALLOCATE(Radiation_Absorption_SpecPercent)

END SUBROUTINE FinalizeRadiation

END MODULE MOD_Radiation_Init
