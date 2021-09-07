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

INTERFACE FinalizeRadiation
  MODULE PROCEDURE FinalizeRadiation
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
CALL prms%CreateRealOption(   'Part-Species[$]-RadiationMass_u',       'Molar mass, kg/kmol', '0.0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(   'Part-Species[$]-RadiationRadius_A',     'Species radius, A', '0.0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(   'Part-Species[$]-Starkex',               'Exponent for the determination of Stark broadening', &
                                                                       '0.0', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(    'Part-Species[$]-NuclCharge',            'Nuclear charge:\n'//&
                                                                       '1: neutral atom\n'//&
                                                                       '2: singly ionized atom', '1', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(   'Radiation-MinWaveLen',                  'Lower wavelength limit for radiation calculation', '0.0')
CALL prms%CreateRealOption(   'Radiation-MaxWaveLen',                  'Upper wavelength limit for radiation calculation','1000.0')
CALL prms%CreateIntOption(    'Radiation-WaveLenDiscr',                'Number of discretization points', '10000')
CALL prms%CreateIntOption(    'Radiation-RadType',                     'Select radiation type:\n'//&
                                                                       '1: particle radiation\n'//&
                                                                       '2: black body radiation\n'//&
                                                                       '3: radiation solver only', '3')
CALL prms%CreateLogicalOption('Radiation-ff',                          'Enable free-free radiation', '.FALSE.')
CALL prms%CreateLogicalOption('Radiation-bf',                          'Enable bound-free radiation (only atomic)', '.FALSE.')
CALL prms%CreateLogicalOption('Radiation-bb-atoms',                    'Enable atomic bound-bound radiation', '.FALSE.')
CALL prms%CreateLogicalOption('Radiation-bb-molecules',                'Enable molecular bound-bound radiation', '.FALSE.')
CALL prms%CreateLogicalOption('Radiation-MacroRadInput',               'Reading in flow field data as radiation input', '.FALSE.')
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
USE MOD_Mesh_Vars,             ONLY : nElems, nGlobalElems
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
    WRITE(UNIT=hilf,FMT='(I0)') iSpec
    RadiationInput(iSpec)%Ttrans(4) = GETREAL('Part-Species'//TRIM(hilf)//'-RadiationTtrans','0.0')
    RadiationInput(iSpec)%Telec = GETREAL('Part-Species'//TRIM(hilf)//'-RadiationTelec','0.0')
    RadiationInput(iSpec)%NumDens = GETREAL('Part-Species'//TRIM(hilf)//'-RadiationNumDens','0.0')
    IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
      RadiationInput(iSpec)%Tvib = GETREAL('Part-Species'//TRIM(hilf)//'-RadiationTvib','0.0')
      RadiationInput(iSpec)%Trot = GETREAL('Part-Species'//TRIM(hilf)//'-RadiationTrot','0.0')
    END IF
    RadiationInput(iSpec)%IonizationEn = GETREAL('Part-Species'//TRIM(hilf)//'-RadiationIonizationEn','0.0')
    RadiationInput(iSpec)%IonizationEn = RadiationInput(iSpec)%IonizationEn *PlanckConst*c*100.
    RadiationInput(iSpec)%Mass = GETREAL('Part-Species'//TRIM(hilf)//'-RadiationMass_u','0.0')
    RadiationInput(iSpec)%Mass = RadiationInput(iSpec)%Mass*1.660539040E-27

    IF((SpecDSMC(iSpec)%InterID .EQ. 1) .OR. (SpecDSMC(iSpec)%InterID .EQ. 10)) THEN !Only for atoms (1) and atomic ions (10)
      CALL Radiation_readin_atoms(iSpec)
      RadiationInput(iSpec)%Radius = GETREAL('Part-Species'//TRIM(hilf)//'-RadiationRadius_A','0.0')
      RadiationInput(iSpec)%Radius = RadiationInput(iSpec)%Radius *1.0E-10
      RadiationInput(iSpec)%Starkex = GETREAL('Part-Species'//TRIM(hilf)//'-Starkex','0.0')
      RadiationInput(iSpec)%NuclCharge = GETINT('Part-Species'//TRIM(hilf)//'-NuclCharge','1')
    END IF

    IF((SpecDSMC(iSpec)%InterID .EQ. 2) .OR. (SpecDSMC(iSpec)%InterID .EQ. 20)) THEN !Only for molecules (2) and molecular ions (20)
      CALL Radiation_readin_molecules(iSpec)
    END IF
  END DO
END IF

RadiationParameter%MinWaveLen   = GETREAL('Radiation-MinWaveLen','0.0')
RadiationParameter%MinWaveLen   = RadiationParameter%MinWaveLen*1.E-9
RadiationParameter%MaxWaveLen   = GETREAL('Radiation-MaxWaveLen','1000.0')
RadiationParameter%MaxWaveLen   = RadiationParameter%MaxWaveLen*1.E-9
RadiationParameter%WaveLenDiscr = GETINT('Radiation-WaveLenDiscr','10000')
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

RadiationSwitches%ff            = GETLOGICAL('Radiation-ff',            '.FALSE.')
RadiationSwitches%bf            = GETLOGICAL('Radiation-bf',            '.FALSE.')
RadiationSwitches%bb_at         = GETLOGICAL('Radiation-bb-atoms',      '.FALSE.')
RadiationSwitches%bb_mol        = GETLOGICAL('Radiation-bb-molecules',  '.FALSE.')
RadiationSwitches%MacroRadInput = GETLOGICAL('Radiation-MacroRadInput', '.FALSE.')

IF (RadiationSwitches%MacroRadInput) CALL MacroscopicRadiationInput()

NumDensElectrons = GETREAL('Radiation-NumDensElectrons','0.0')
TElectrons       = GETREAL('Radiation-TElectrons',      '0.0')

#if USE_MPI
  ! allocate shared array for Radiation_Emission/Absorption_Spec  
CALL Allocate_Shared((/RadiationParameter%WaveLenDiscr,nComputeNodeElems/), Radiation_Emission_Spec_Shared_Win,Radiation_Emission_Spec_Shared)
CALL MPI_WIN_LOCK_ALL(0,Radiation_Emission_Spec_Shared_Win,IERROR)
CALL Allocate_Shared_Test((/INT(RadiationParameter%WaveLenDiscr,IK)*INT(nGlobalElems,IK)/),Radiation_Absorption_Spec_Shared_Win,Radiation_Absorption_Spec_Shared)
CALL MPI_WIN_LOCK_ALL(0,Radiation_Absorption_Spec_Shared_Win,IERROR)
CALL Allocate_Shared((/nSpecies,nComputeNodeElems,2/), Radiation_ElemEnergy_Species_Shared_Win,Radiation_ElemEnergy_Species_Shared)
CALL MPI_WIN_LOCK_ALL(0,Radiation_ElemEnergy_Species_Shared_Win,IERROR)

Radiation_Emission_spec => Radiation_Emission_spec_Shared
Radiation_Absorption_spec(1:RadiationParameter%WaveLenDiscr ,1:nGlobalElems) => Radiation_Absorption_spec_Shared
Radiation_ElemEnergy_Species => Radiation_ElemEnergy_Species_Shared
#else
! allocate local array for ElemInfo
ALLOCATE(Radiation_Emission_spec(RadiationParameter%WaveLenDiscr,nElems))
ALLOCATE(Radiation_Absorption_spec(RadiationParameter%WaveLenDiscr,nElems))
ALLOCATE(Radiation_ElemEnergy_Species(nSpecies,nElems,2))
#endif  /*USE_MPI*/

#if USE_MPI
  firstElem = INT(REAL( myComputeNodeRank   *nComputeNodeElems)/REAL(nComputeNodeProcessors))+1
  lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeElems)/REAL(nComputeNodeProcessors))
#else
  firstElem = 1
  lastElem  = nElems
#endif

DO iElem = firstElem, lastElem
  Radiation_Emission_spec(:,iElem) = 0.0
  Radiation_Absorption_spec(:,GetGlobalElemID(iElem)) = 0.0
  Radiation_ElemEnergy_Species(:,iElem,:) =0.0
END DO
#if USE_MPI
  CALL BARRIER_AND_SYNC(Radiation_Emission_Spec_Shared_Win ,MPI_COMM_SHARED)
  CALL BARRIER_AND_SYNC(Radiation_ElemEnergy_Species_Shared_Win ,MPI_COMM_SHARED)
  CALL BARRIER_AND_SYNC(Radiation_Absorption_Spec_Shared_Win ,MPI_COMM_SHARED)
  IF(nLeaderGroupProcs.GT.1)THEN
    IF(myComputeNodeRank.EQ.0)THEN
      CALL MPI_ALLGATHERV( MPI_IN_PLACE                  &
                     , 0                             &
                     , MPI_DATATYPE_NULL             &
                     , Radiation_Absorption_Spec_Shared  &
                     , RadiationParameter%WaveLenDiscr *recvcountElem   &
                     , RadiationParameter%WaveLenDiscr *displsElem      &
                     , MPI_DOUBLE_PRECISION          &
                     , MPI_COMM_LEADERS_SHARED       &
                     , IERROR)
    END IF
  END IF  
  CALL BARRIER_AND_SYNC(Radiation_Absorption_Spec_Shared_Win ,MPI_COMM_SHARED)
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
  USE MOD_Radiation_Vars            ,ONLY: MacroRadInputParameters, MacroRadInputParameters_Shared, MacroRadInputParameters_Shared_Win
  USE MOD_Mesh_Tools                ,ONLY: GetCNElemID
  USE MOD_ReadInTools
  USE MOD_Particle_Mesh_Vars        ,ONLY: nComputeNodeElems
#if USE_MPI
!USE MOD_MPI_Shared_Vars
  USE MOD_MPI_Shared
  USE MOD_MPI_Shared_Vars
#endif
!  USE MOD_Utils                     ,ONLY: BubbleSortID !Laux
!  USE MOD_Particle_Mesh_Vars        ,ONLY: GEO !Laux
  ! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES
  INTEGER                           :: nVar_HDF5, N_HDF5, nElems_HDF5, iVar, iSpec, iElem, CNElemID
  REAL, ALLOCATABLE                 :: ElemData_HDF5(:,:)
  CHARACTER(LEN=300)                :: MacroRadiationInputFile
!  INTEGER, ALLOCATABLE              :: SortElemInd(:)  !Laux
!  REAL, ALLOCATABLE                 :: SortElemYPos(:) !Laux
  !===================================================================================================================================

  MacroRadiationInputFile = GETSTR('Radiation-MacroInput-Filename')
  CALL OpenDataFile(MacroRadiationInputFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)

  CALL GetDataProps('ElemData',nVar_HDF5,N_HDF5,nElems_HDF5)

  ALLOCATE(MacroRadInputParameters(1:nElems,1:nSpecies,1:5))
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


! --- uncomment for Laux Test Case ----------------------------------------------------------------
! reanrrages Elem IDs for 1xnElemsx1 meshes (for 3D and 2D rotationally symmetric meshes)

!  ALLOCATE(SortElemInd(nElems), SortElemYPos(nElems))
!  DO iElem = 1, nElems
!    SortElemInd(iElem) = iElem
!  END DO
!  SortElemYPos(:) = -GEO%ElemMidPoint(2,:)
!  CALL BubbleSortID(SortElemYPos, SortElemInd, nElems)


!  iVar = 1
!  DO iSpec = 1, nSpecies
!    DO iElem = 1, nElems
!      MacroRadInputParameters(SortElemInd(iElem),iSpec,1) = ElemData_HDF5(iVar+ 6,iElem) density
!      MacroRadInputParameters(SortElemInd(iElem),iSpec,2) = ElemData_HDF5(iVar+ 7,iElem) T_vib
!      MacroRadInputParameters(SortElemInd(iElem),iSpec,3) = ElemData_HDF5(iVar+ 8,iElem) T_rot
!      MacroRadInputParameters(SortElemInd(iElem),iSpec,4) = ElemData_HDF5(iVar+ 9,iElem) T_elec
!      MacroRadInputParameters(SortElemInd(iElem),iSpec,5) = ElemData_HDF5(iVar+11,iElem) T_mean
!    END DO
!    iVar = iVar + DSMC_NVARS
!  END DO
! -------------------------------------------------------------------------------------------------

  iVar = 1
  DO iSpec = 1, nSpecies
   DO iElem = 1, nElems
     CNElemID = GetCNElemID(iElem+offsetElem)
     MacroRadInputParameters(CNElemID,iSpec,1) = MAX(0.,ElemData_HDF5(iVar+ 6,iElem)) !density
     MacroRadInputParameters(CNElemID,iSpec,2) = MAX(0.,ElemData_HDF5(iVar+ 7,iElem)) !T_vib
     MacroRadInputParameters(CNElemID,iSpec,3) = MAX(0.,ElemData_HDF5(iVar+ 8,iElem)) !T_rot
     MacroRadInputParameters(CNElemID,iSpec,4) = MAX(0.,ElemData_HDF5(iVar+ 9,iElem)) !T_elec
     MacroRadInputParameters(CNElemID,iSpec,5) = MAX(0.,ElemData_HDF5(iVar+11,iElem)) !T_mean
   END DO
   iVar = iVar + DSMC_NVARS
  END DO

#if USE_MPI
  CALL BARRIER_AND_SYNC(MacroRadInputParameters_Shared_Win ,MPI_COMM_SHARED)
#endif

  DEALLOCATE(ElemData_HDF5)

END SUBROUTINE MacroscopicRadiationInput


SUBROUTINE FinalizeRadiation()
!===================================================================================================================================
!> Deallocating radiation variables
!===================================================================================================================================
! MODULES
USE MOD_Radiation_Vars
#if USE_MPI
!USE MOD_MPI_Shared_Vars    !,ONLY: MPI_COMM_SHARED
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
CALL UNLOCK_AND_FREE(Radiation_Emission_Spec_Shared_Win)
CALL UNLOCK_AND_FREE(Radiation_Absorption_Spec_Shared_Win)
CALL UNLOCK_AND_FREE(Radiation_ElemEnergy_Species_Shared_Win)
#endif /*USE_MPI*/
ADEALLOCATE(Radiation_Emission_Spec)
ADEALLOCATE(Radiation_Absorption_Spec)
ADEALLOCATE(Radiation_ElemEnergy_Species)

END SUBROUTINE FinalizeRadiation




END MODULE MOD_Radiation_Init
