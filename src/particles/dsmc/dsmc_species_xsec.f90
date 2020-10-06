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

MODULE MOD_DSMC_SpecXSec
!===================================================================================================================================
! Contains the Argon Ionization
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE

PUBLIC :: InterpolateCrossSection, XSec_CalcCollisionProb, XSec_CalcReactionProb, XSec_CalcVibRelaxProb
!===================================================================================================================================

CONTAINS

SUBROUTINE MCC_Init()
!===================================================================================================================================
!> Read-in of the collision and vibrational cross-section database and initialization of the null collision method.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Globals_Vars          ,ONLY: ElementaryCharge
USE MOD_PARTICLE_Vars         ,ONLY: nSpecies
USE MOD_DSMC_Vars             ,ONLY: BGGas, SpecDSMC, XSec_Database, SpecXSec, XSec_NullCollision, XSec_Relaxation, CollInf
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER       :: iSpec, jSpec, iCase
REAL          :: TotalProb
INTEGER       :: iVib, nVib, iStep, MaxDim
!===================================================================================================================================

XSec_Database = TRIM(GETSTR('Particles-CollXSec-Database'))
IF(BGGas%NumberOfSpecies.GT.0) THEN
  XSec_NullCollision = GETLOGICAL('Particles-CollXSec-NullCollision')
ELSE
  XSec_NullCollision = .FALSE.
END IF
XSec_Relaxation = .FALSE.

IF(TRIM(XSec_Database).EQ.'none') THEN
  CALL abort(&
  __STAMP__&
  ,'ERROR: No database for the collision cross-section given!')
END IF

ALLOCATE(SpecXSec(CollInf%NumCase))
SpecXSec(:)%UseCollXSec = .FALSE.
SpecXSec(:)%UseVibXSec = .FALSE.

DO iSpec = 1, nSpecies
  TotalProb = 0.
  DO jSpec = iSpec, nSpecies
    iCase = CollInf%Coll_Case(iSpec,jSpec)
    ! Skip species, which shall not be treated with collision cross-sections
    IF(.NOT.SpecDSMC(iSpec)%UseCollXSec.AND..NOT.SpecDSMC(jSpec)%UseCollXSec.AND. &
       .NOT.SpecDSMC(iSpec)%UseVibXSec.AND..NOT.SpecDSMC(jSpec)%UseVibXSec) CYCLE
    ! Skip pairing with itself and pairing with other particle species, if background gas is active
    IF(BGGas%NumberOfSpecies.GT.0) THEN
      IF(iSpec.EQ.jSpec) CYCLE
      IF(.NOT.BGGas%BackgroundSpecies(iSpec).AND..NOT.BGGas%BackgroundSpecies(jSpec)) CYCLE
    END IF
    ! Read-in cross-section data for collisions of particles, allocating CollXSecData within the following routine
    CALL ReadCollXSec(iCase, iSpec, jSpec)
    IF(SpecXSec(iCase)%UseCollXSec) THEN
      IF(SpecDSMC(iSpec)%UseCollXSec.AND.SpecDSMC(jSpec)%UseCollXSec) THEN
        CALL abort(&
          __STAMP__&
          ,'ERROR: Both species defined to use collisional cross-section, define only the source species with UseCollXSec!')
      END IF
      ! Store the energy value in J (read-in was in eV)
      SpecXSec(iCase)%CollXSecData(1,:) = SpecXSec(iCase)%CollXSecData(1,:) * ElementaryCharge
      IF(XSec_NullCollision) THEN
        ! Determine the maximum collision frequency for the null collision method
        CALL DetermineNullCollProb(iCase,iSpec,jSpec)
        TotalProb = TotalProb + SpecXSec(iCase)%ProbNull
        IF(TotalProb.GT.1.0) THEN
          CALL abort(&
          __STAMP__&
          ,'ERROR: Total null collision probability is above unity. Please reduce the time step! Probability is: '&
          ,RealInfoOpt=TotalProb)
        END IF
      END IF
    END IF
    ! Read-in vibrational cross sections
    IF(SpecDSMC(iSpec)%UseVibXSec.OR.SpecDSMC(jSpec)%UseVibXSec) CALL ReadVibXSec(iCase, iSpec, jSpec)
    ! Vibrational relaxation probabilities: Interpolate and store the probability at the effective cross-section levels
    IF(SpecXSec(iCase)%UseVibXSec) THEN
      XSec_Relaxation = .TRUE.
      nVib = SIZE(SpecXSec(iCase)%VibMode)
      DO iVib = 1, nVib
        ! Store the energy value in J (read-in was in eV)
        SpecXSec(iCase)%VibMode(iVib)%XSecData(1,:) = SpecXSec(iCase)%VibMode(iVib)%XSecData(1,:) * ElementaryCharge
      END DO
      IF(SpecXSec(iCase)%UseCollXSec) THEN
        ! Effective collision cross-sections are available
        MaxDim = SIZE(SpecXSec(iCase)%CollXSecData,2)
        ALLOCATE(SpecXSec(iCase)%VibXSecData(1:2,1:MaxDim))
        ! Using the same energy intervals as for the effective cross-sections
        SpecXSec(iCase)%VibXSecData(1,:) = SpecXSec(iCase)%CollXSecData(1,:)
        SpecXSec(iCase)%VibXSecData(2,:) = 0.
        DO iVib = 1, nVib
          ! Interpolate the vibrational cross section at the energy levels of the effective collision cross section and sum-up the
          ! vibrational probability (vibrational cross-section divided by the effective)
          DO iStep = 1, MaxDim
            IF(SpecXSec(iCase)%CollXSecData(2,iStep).GT.0.0) THEN
              SpecXSec(iCase)%VibXSecData(2,iStep) = SpecXSec(iCase)%VibXSecData(2,iStep) + &
                InterpolateCrossSection_Vib(iCase,iVib,SpecXSec(iCase)%CollXSecData(1,iStep)) &
                                                            / SpecXSec(iCase)%CollXSecData(2,iStep)
            END IF
          END DO
        END DO
      END IF    ! SpecXSec(iCase)%UseCollXSec
    END IF      ! SpecXSec(iCase)%UseVibXSec
  END DO        ! jSpec = iSpec, nSpecies
END DO          ! iSpec = 1, nSpecies

END SUBROUTINE MCC_Init


SUBROUTINE ReadCollXSec(iCase,iSpec,jSpec)
!===================================================================================================================================
!> Read-in of collision cross-sections from a given database. Dataset name is composed of SpeciesName-SpeciesName (e.g. Ar-electron)
!> Trying to swap the species indices if dataset not found.
!===================================================================================================================================
! use module
USE MOD_io_hdf5
USE MOD_Globals
USE MOD_DSMC_Vars                 ,ONLY: XSec_Database, SpecXSec, SpecDSMC
USE MOD_HDF5_Input                ,ONLY: DatasetExists
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: iCase, iSpec, jSpec
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=64)                 :: dsetname, dsetname2, spec_pair
INTEGER                           :: err
INTEGER(HSIZE_T), DIMENSION(2)    :: dims,sizeMax
INTEGER(HID_T)                    :: file_id_dsmc                       ! File identifier
INTEGER(HID_T)                    :: dset_id_dsmc                       ! Dataset identifier
INTEGER(HID_T)                    :: filespace                          ! filespace identifier
LOGICAL                           :: DatasetFound
!===================================================================================================================================
spec_pair = TRIM(SpecDSMC(jSpec)%Name)//'-'//TRIM(SpecDSMC(iSpec)%Name)
SWRITE(UNIT_StdOut,'(A)') 'Read collision cross section for '//TRIM(spec_pair)//' from '//TRIM(XSec_Database)

DatasetFound = .FALSE.

! Initialize FORTRAN interface.
CALL H5OPEN_F(err)

! Check if file exists
IF(.NOT.FILEEXISTS(XSec_Database)) THEN
  CALL abort(__STAMP__,'ERROR: Database '//TRIM(XSec_Database)//' does not exist.')
END IF

! Open the file.
CALL H5FOPEN_F (TRIM(XSec_Database), H5F_ACC_RDONLY_F, file_id_dsmc, err)

dsetname = TRIM('/'//TRIM(spec_pair)//'/EFFECTIVE')
CALL DatasetExists(File_ID_DSMC,TRIM(dsetname),DatasetFound)

! Check if the dataset exist
IF(.NOT.DatasetFound) THEN
  ! Try to swap the species names
  dsetname = '/'//TRIM(SpecDSMC(iSpec)%Name)//'-'//TRIM(SpecDSMC(jSpec)%Name)//'/EFFECTIVE'
  CALL DatasetExists(File_ID_DSMC,TRIM(dsetname),DataSetFound)
  IF(DatasetFound) THEN
    spec_pair = TRIM(SpecDSMC(iSpec)%Name)//'-'//TRIM(SpecDSMC(jSpec)%Name)
    SpecXSec(iCase)%UseCollXSec = .TRUE.
  ELSE
    dsetname2 = TRIM(spec_pair)//'/EFFECTIVE'
    SWRITE(UNIT_StdOut,'(A)') 'Dataset not found: ['//TRIM(dsetname2)//']. Using standard collision modelling.'
    SpecXSec(iCase)%UseCollXSec = .FALSE.
  END IF
ELSE
  SpecXSec(iCase)%UseCollXSec = .TRUE.
END IF

IF(SpecXSec(iCase)%UseCollXSec) THEN
  ! Open the dataset.
  CALL H5DOPEN_F(file_id_dsmc, dsetname, dset_id_dsmc, err)
  ! Get the file space of the dataset.
  CALL H5DGET_SPACE_F(dset_id_dsmc, FileSpace, err)
  ! get size
  CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, dims, SizeMax, err)
  ! Read-in the effective cross-sections
  ALLOCATE(SpecXSec(iCase)%CollXSecData(1:2,1:dims(2)))
  SpecXSec(iCase)%CollXSecData = 0.
  ! read data
  CALL H5DREAD_F(dset_id_dsmc, H5T_NATIVE_DOUBLE, SpecXSec(iCase)%CollXSecData(1:2,1:dims(2)), dims, err)
END IF

! Close the file.
CALL H5FCLOSE_F(file_id_dsmc, err)
! Close FORTRAN interface.
CALL H5CLOSE_F(err)

END SUBROUTINE ReadCollXSec


SUBROUTINE ReadVibXSec(iCase,iSpec,jSpec)
!===================================================================================================================================
!> Read-in of vibrational cross-sections from a given database. Dataset name is composed of SpeciesName-SpeciesName (e.g. Ar-electron)
!> Trying to swap the species indices if dataset not found.
!===================================================================================================================================
! use module
USE MOD_io_hdf5
USE MOD_Globals
USE MOD_DSMC_Vars                 ,ONLY: XSec_Database, SpecXSec, SpecDSMC
USE MOD_HDF5_Input                ,ONLY: DatasetExists
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: iCase, iSpec, jSpec
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=64)                 :: dsetname, spec_pair, groupname
INTEGER                           :: err, nVar
INTEGER(HSIZE_T), DIMENSION(2)    :: dims,sizeMax
INTEGER(HID_T)                    :: file_id_dsmc                       ! File identifier
INTEGER(HID_T)                    :: group_id                           ! Group identifier
INTEGER(HID_T)                    :: dset_id_dsmc                       ! Dataset identifier
INTEGER(HID_T)                    :: filespace                          ! filespace identifier
INTEGER(SIZE_T)                   :: size                               ! Size of name
INTEGER(HSIZE_T)                  :: iVib                               ! Index
LOGICAL                           :: GroupFound
INTEGER                           :: storage, nVib, max_corder
!===================================================================================================================================
spec_pair = TRIM(SpecDSMC(jSpec)%Name)//'-'//TRIM(SpecDSMC(iSpec)%Name)
SWRITE(UNIT_StdOut,'(A)') 'Read vibrational cross section for '//TRIM(spec_pair)//' from '//TRIM(XSec_Database)

GroupFound = .FALSE.
SpecXSec(iCase)%UseVibXSec = .FALSE.

! Initialize FORTRAN interface.
CALL H5OPEN_F(err)

! Check if file exists
IF(.NOT.FILEEXISTS(XSec_Database)) THEN
  CALL abort(__STAMP__,'ERROR: Database '//TRIM(XSec_Database)//' does not exist.')
END IF

! Open the file.
CALL H5FOPEN_F (TRIM(XSec_Database), H5F_ACC_RDONLY_F, file_id_dsmc, err)

! Check if the species pair group exists
CALL H5LEXISTS_F(file_id_dsmc, TRIM(spec_pair), GroupFound, err)
IF(.NOT.GroupFound) THEN
  ! Try to swap the species names
  spec_pair = TRIM(SpecDSMC(iSpec)%Name)//'-'//TRIM(SpecDSMC(jSpec)%Name)
  CALL H5LEXISTS_F(file_id_dsmc, TRIM(spec_pair), GroupFound, err)
  IF(.NOT.GroupFound) THEN
    SWRITE(UNIT_StdOut,'(A)') 'No vibrational excitation cross sections found in database, using constant read-in values.'
    RETURN
  END IF
END IF

! Check if the vibrational cross-section group exists
groupname = TRIM(spec_pair)//'/VIBRATION/'
CALL H5LEXISTS_F(file_id_dsmc, TRIM(groupname), GroupFound, err)
IF(.NOT.GroupFound) THEN
  SWRITE(UNIT_StdOut,'(A)') 'No vibrational excitation cross sections found in database, using constant read-in values.'
  RETURN
END IF

IF(GroupFound) THEN
  CALL H5GOPEN_F(file_id_dsmc,TRIM(groupname), group_id, err)
  call H5Gget_info_f(group_id, storage, nVib,max_corder, err)
  ! If cross-section data is found, set the corresponding flag
  IF(nVib.GT.0) THEN
    SWRITE(UNIT_StdOut,'(A,I3,A)') 'Found ', nVib,' vibrational excitation cross section(s) in database.'
    SpecXSec(iCase)%UseVibXSec = .TRUE.
    nVar = 3
  ELSE
    SWRITE(UNIT_StdOut,'(A)') 'No vibrational excitation cross sections found in database, using constant read-in values.'
  END IF
ELSE
  SWRITE(UNIT_StdOut,'(A)') 'No vibrational excitation cross sections found in database, using constant read-in values.'
END IF

IF(SpecXSec(iCase)%UseVibXSec) THEN
  ALLOCATE(SpecXSec(iCase)%VibMode(1:nVib))
  DO iVib = 0, nVib-1
    ! Get name and size of name
    CALL H5Lget_name_by_idx_f(group_id, ".", H5_INDEX_NAME_F, H5_ITER_INC_F, iVib, dsetname, err, size)
    dsetname = TRIM(groupname)//TRIM(dsetname)
    ! Open the dataset.
    CALL H5DOPEN_F(file_id_dsmc, dsetname, dset_id_dsmc, err)
    ! Get the file space of the dataset.
    CALL H5DGET_SPACE_F(dset_id_dsmc, FileSpace, err)
    ! get size
    CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, dims, SizeMax, err)
    ALLOCATE(SpecXSec(iCase)%VibMode(iVib+1)%XSecData(dims(1),dims(2)))
    ! read data
    CALL H5DREAD_F(dset_id_dsmc, H5T_NATIVE_DOUBLE, SpecXSec(iCase)%VibMode(iVib+1)%XSecData, dims, err)
  END DO
  ! Close the group
  CALL H5GCLOSE_F(group_id,err)
END IF

! Close the file.
CALL H5FCLOSE_F(file_id_dsmc, err)
! Close FORTRAN interface.
CALL H5CLOSE_F(err)

END SUBROUTINE ReadVibXSec


SUBROUTINE DetermineNullCollProb(iCase,iSpec,jSpec)
!===================================================================================================================================
!> Routine for the MCC method: calculates the maximal collision frequency for a given species and the collision probability
!===================================================================================================================================
! MODULES
USE MOD_ReadInTools
USE MOD_Globals_Vars          ,ONLY: Pi
USE MOD_Particle_Vars         ,ONLY: Species, ManualTimeStep
USE MOD_DSMC_Vars             ,ONLY: BGGas, SpecXSec
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)            :: iCase                            !< Case index
INTEGER,INTENT(IN)            :: iSpec
INTEGER,INTENT(IN)            :: jSpec
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: MaxDOF, bggSpec
REAL                          :: MaxCollFreq, Mass
REAL,ALLOCATABLE              :: Velocity(:)
!===================================================================================================================================

! Select the background species as the target cloud and use the mass of particle species
IF(BGGas%BackgroundSpecies(iSpec)) THEN
  bggSpec = BGGas%MapSpecToBGSpec(iSpec)
  Mass = Species(jSpec)%MassIC
ELSE
  bggSpec = BGGas%MapSpecToBGSpec(jSpec)
  Mass = Species(iSpec)%MassIC
END IF

MaxDOF = SIZE(SpecXSec(iCase)%CollXSecData,2)
ALLOCATE(Velocity(MaxDOF))

! Determine the mean relative velocity at the given energy level
Velocity(1:MaxDOF) = SQRT(2.) * SQRT(8.*SpecXSec(iCase)%CollXSecData(1,1:MaxDOF)/(Pi*Mass))

! Calculate the maximal collision frequency
MaxCollFreq = MAXVAL(Velocity(1:MaxDOF) * SpecXSec(iCase)%CollXSecData(2,1:MaxDOF) * BGGas%NumberDensity(bggSpec))

! Determine the collision probability
SpecXSec(iCase)%ProbNull = 1. - EXP(-MaxCollFreq*ManualTimeStep)

DEALLOCATE(Velocity)

END SUBROUTINE DetermineNullCollProb


PURE REAL FUNCTION InterpolateCrossSection(iCase,CollisionEnergy)
!===================================================================================================================================
!> Interpolate the collision cross-section [m^2] from the available data at the given collision energy [J]
!> Collision energies below and above the given data will be set at the first and last level of the data set
!> Note: Requires the data to be sorted by ascending energy values
!> Assumption: First species given is the particle species, second species input is the background gas species
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars             ,ONLY: SpecXSec
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)            :: iCase                            !< Case index
REAL,INTENT(IN)               :: CollisionEnergy                  !< Collision energy in [J]
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iDOF, MaxDOF
!===================================================================================================================================

InterpolateCrossSection = 0.
MaxDOF = SIZE(SpecXSec(iCase)%CollXSecData,2)

IF(CollisionEnergy.GT.SpecXSec(iCase)%CollXSecData(1,MaxDOF)) THEN 
  ! If the collision energy is greater than the maximal value, get the cross-section of the last level and leave routine
  InterpolateCrossSection = SpecXSec(iCase)%CollXSecData(2,MaxDOF)
  ! Leave routine
  RETURN
ELSE IF(CollisionEnergy.LE.SpecXSec(iCase)%CollXSecData(1,1)) THEN
  ! If collision energy is below the minimal value, get the cross-section of the first level and leave routine
  InterpolateCrossSection = SpecXSec(iCase)%CollXSecData(2,1)
  ! Leave routine
  RETURN
END IF

DO iDOF = 1, MaxDOF
  ! Check if the stored energy value is above the collision energy
  IF(SpecXSec(iCase)%CollXSecData(1,iDOF).GT.CollisionEnergy) THEN
    ! Interpolate the cross-section from the data set using the current and the energy level below
    InterpolateCrossSection = SpecXSec(iCase)%CollXSecData(2,iDOF-1) &
                                + (CollisionEnergy - SpecXSec(iCase)%CollXSecData(1,iDOF-1)) &
                                / (SpecXSec(iCase)%CollXSecData(1,iDOF) - SpecXSec(iCase)%CollXSecData(1,iDOF-1)) &
                                * (SpecXSec(iCase)%CollXSecData(2,iDOF) - SpecXSec(iCase)%CollXSecData(2,iDOF-1))
    ! Leave routine and do not finish DO loop
    RETURN
  END IF
END DO

END FUNCTION InterpolateCrossSection


PURE REAL FUNCTION InterpolateCrossSection_Vib(iCase,iVib,CollisionEnergy)
!===================================================================================================================================
!> Interpolate the vibrational cross-section data for specific vibrational level at the given collision energy
!> Note: Requires the data to be sorted by ascending energy values
!> Assumption: First species given is the particle species, second species input is the background gas species
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars             ,ONLY: SpecXSec
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)            :: iCase                            !< Case index
INTEGER,INTENT(IN)            :: iVib                             !< 
REAL,INTENT(IN)               :: CollisionEnergy                  !< Collision energy in [J]
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iDOF, MaxDOF
!===================================================================================================================================

InterpolateCrossSection_Vib = 0.
MaxDOF = SIZE(SpecXSec(iCase)%VibMode(iVib)%XSecData,2)

IF(CollisionEnergy.GT.SpecXSec(iCase)%VibMode(iVib)%XSecData(1,MaxDOF)) THEN 
  ! If the collision energy is greater than the maximal value, get the cross-section of the last level and leave routine
  InterpolateCrossSection_Vib = SpecXSec(iCase)%VibMode(iVib)%XSecData(2,MaxDOF)
  ! Leave routine
  RETURN
ELSE IF(CollisionEnergy.LE.SpecXSec(iCase)%VibMode(iVib)%XSecData(1,1)) THEN
  ! If collision energy is below the minimal value, get the cross-section of the first level and leave routine
  InterpolateCrossSection_Vib = SpecXSec(iCase)%VibMode(iVib)%XSecData(2,1)
  ! Leave routine
  RETURN
END IF

DO iDOF = 1, MaxDOF
  ! Check if the stored energy value is above the collision energy
  IF(SpecXSec(iCase)%VibMode(iVib)%XSecData(1,iDOF).GT.CollisionEnergy) THEN
    ! Interpolate the cross-section from the data set using the current and the energy level below
    InterpolateCrossSection_Vib = SpecXSec(iCase)%VibMode(iVib)%XSecData(2,iDOF-1) &
              + (CollisionEnergy - SpecXSec(iCase)%VibMode(iVib)%XSecData(1,iDOF-1)) &
              / (SpecXSec(iCase)%VibMode(iVib)%XSecData(1,iDOF) - SpecXSec(iCase)%VibMode(iVib)%XSecData(1,iDOF-1)) &
              * (SpecXSec(iCase)%VibMode(iVib)%XSecData(2,iDOF) - SpecXSec(iCase)%VibMode(iVib)%XSecData(2,iDOF-1))
    ! Leave routine and do not finish DO loop
    RETURN
  END IF
END DO

END FUNCTION InterpolateCrossSection_Vib


PURE REAL FUNCTION InterpolateVibRelaxProb(iCase,CollisionEnergy)
!===================================================================================================================================
!> Interpolate the vibrational relaxation probability at the same intervals as the effective collision cross-section
!> Note: Requires the data to be sorted by ascending energy values
!> Assumption: First species given is the particle species, second species input is the background gas species
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars             ,ONLY: SpecXSec
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)            :: iCase                            !< Case index
REAL,INTENT(IN)               :: CollisionEnergy                  !< Collision energy in [J]
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iDOF, MaxDOF
!===================================================================================================================================

InterpolateVibRelaxProb = 0.
MaxDOF = SIZE(SpecXSec(iCase)%VibXSecData,2)

IF(CollisionEnergy.GT.SpecXSec(iCase)%VibXSecData(1,MaxDOF)) THEN 
  ! If the collision energy is greater than the maximal value, get the cross-section of the last level and leave routine
  InterpolateVibRelaxProb = SpecXSec(iCase)%VibXSecData(2,MaxDOF)
  ! Leave routine
  RETURN
ELSE IF(CollisionEnergy.LE.SpecXSec(iCase)%VibXSecData(1,1)) THEN
  ! If collision energy is below the minimal value, get the cross-section of the first level and leave routine
  InterpolateVibRelaxProb = SpecXSec(iCase)%VibXSecData(2,1)
  ! Leave routine
  RETURN
END IF

DO iDOF = 1, MaxDOF
  ! Check if the stored energy value is above the collision energy
  IF(SpecXSec(iCase)%VibXSecData(1,iDOF).GT.CollisionEnergy) THEN
    ! Interpolate the cross-section from the data set using the current and the energy level below
    InterpolateVibRelaxProb = SpecXSec(iCase)%VibXSecData(2,iDOF-1) &
              + (CollisionEnergy - SpecXSec(iCase)%VibXSecData(1,iDOF-1)) &
              / (SpecXSec(iCase)%VibXSecData(1,iDOF) - SpecXSec(iCase)%VibXSecData(1,iDOF-1)) &
              * (SpecXSec(iCase)%VibXSecData(2,iDOF) - SpecXSec(iCase)%VibXSecData(2,iDOF-1))
    ! Leave routine and do not finish DO loop
    RETURN
  END IF
END DO

END FUNCTION InterpolateVibRelaxProb


PURE REAL FUNCTION XSec_CalcCollisionProb(iPair,SpecNum1,SpecNum2,CollCaseNum,MacroParticleFactor,Volume,dtCell)
!===================================================================================================================================
!> Calculate the collision probability if collision cross-section data is used. Can be utilized in combination with the regular
!> DSMC collision calculation probability.
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars             ,ONLY: SpecXSec, SpecDSMC, Coll_pData, CollInf, BGGas, XSec_NullCollision
USE MOD_Particle_Vars         ,ONLY: PartSpecies, Species, PartState
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)            :: iPair
REAL,INTENT(IN)               :: SpecNum1, SpecNum2, CollCaseNum, MacroParticleFactor, Volume, dtCell
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: CollEnergy, VeloSquare, SpecNumTarget, SpecNumSource
INTEGER                       :: XSecSpec, XSecPart, targetSpec, iPart_p1, iPart_p2, iSpec_p1, iSpec_p2, iCase
!===================================================================================================================================

iPart_p1 = Coll_pData(iPair)%iPart_p1; iPart_p2 = Coll_pData(iPair)%iPart_p2
iSpec_p1 = PartSpecies(iPart_p1);      iSpec_p2 = PartSpecies(iPart_p2)
iCase = CollInf%Coll_Case(iSpec_p1,iSpec_p2)

IF(SpecXSec(iCase)%UseCollXSec) THEN
  IF(SpecDSMC(iSpec_p1)%UseCollXSec) THEN
    XSecSpec = iSpec_p1; targetSpec = iSpec_p2; XSecPart = iPart_p1; SpecNumTarget = SpecNum2; SpecNumSource = SpecNum1
  ELSE
    XSecSpec = iSpec_p2; targetSpec = iSpec_p1; XSecPart = iPart_p2; SpecNumTarget = SpecNum1; SpecNumSource = SpecNum2
  END IF
  ! Using the kinetic energy of the particle (as is described in Vahedi1995 and Birdsall1991)
  VeloSquare = DOT_PRODUCT(PartState(4:6,XSecPart),PartState(4:6,XSecPart))
  CollEnergy = 0.5 * Species(XSecSpec)%MassIC * VeloSquare
  ! Calculate the collision probability
  IF(BGGas%BackgroundSpecies(targetSpec)) THEN
    ! Correct the collision probability in the case of the second species being a background species as the number of pairs
    ! is either determined based on the null collision probability or in the case of mixture on the species fraction
    IF(XSec_NullCollision) THEN
      XSec_CalcCollisionProb = (1. - EXP(-SQRT(VeloSquare) * InterpolateCrossSection(iCase,CollEnergy) &
            * SpecNumTarget * MacroParticleFactor / Volume * dtCell)) / SpecXSec(iCase)%ProbNull
    ELSE
      XSec_CalcCollisionProb = (1. - EXP(-SQRT(VeloSquare) * InterpolateCrossSection(iCase,CollEnergy) &
            * SpecNumTarget * MacroParticleFactor / Volume * dtCell)) / BGGas%SpeciesFraction(BGGas%MapSpecToBGSpec(targetSpec))
    END IF
  ELSE
    XSec_CalcCollisionProb = (1. - EXP(-SQRT(VeloSquare)*InterpolateCrossSection(iCase,CollEnergy) &
            * SpecNumTarget * MacroParticleFactor/Volume*dtCell)) * SpecNumSource / CollInf%Coll_CaseNum(iCase)
  END IF
ELSE
  XSec_CalcCollisionProb = SpecNum1*SpecNum2/(1 + CollInf%KronDelta(Coll_pData(iPair)%PairType))  &
          * CollInf%Cab(Coll_pData(iPair)%PairType)                           & ! Cab species comb fac
          * MacroParticleFactor / CollCaseNum                                                     &
          * Coll_pData(iPair)%CRela2 ** (0.5-CollInf%omega(iSpec_p1,iSpec_p2)) &
          * dtCell / Volume
END IF

END FUNCTION XSec_CalcCollisionProb


SUBROUTINE XSec_CalcVibRelaxProb(iPair,SpecNum1,SpecNum2,MacroParticleFactor,Volume,dtCell)
!===================================================================================================================================
!> Calculate the collision probability if collision cross-section data is used. Can be utilized in combination with the regular
!> DSMC collision calculation probability.
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars             ,ONLY: SpecXSec, SpecDSMC, Coll_pData, CollInf, BGGas
USE MOD_Particle_Vars         ,ONLY: PartSpecies, Species, PartState
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)            :: iPair
REAL,INTENT(IN)               :: SpecNum1, SpecNum2, MacroParticleFactor, Volume, dtCell
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: CollEnergy, VeloSquare, SpecNumTarget, SpecNumSource
INTEGER                       :: XSecSpec, XSecPart, targetSpec, iPart_p1, iPart_p2, iSpec_p1, iSpec_p2, iCase, iVib, nVib
!===================================================================================================================================

iPart_p1 = Coll_pData(iPair)%iPart_p1; iPart_p2 = Coll_pData(iPair)%iPart_p2
iSpec_p1 = PartSpecies(iPart_p1);      iSpec_p2 = PartSpecies(iPart_p2)
iCase = CollInf%Coll_Case(iSpec_p1,iSpec_p2)
SpecXSec(iCase)%VibProb = 0.

IF(SpecXSec(iCase)%UseVibXSec) THEN
  IF(SpecDSMC(iSpec_p1)%UseVibXSec) THEN
    XSecSpec = iSpec_p1; targetSpec = iSpec_p2; XSecPart = iPart_p1; SpecNumTarget = SpecNum2; SpecNumSource = SpecNum1
  ELSE
    XSecSpec = iSpec_p2; targetSpec = iSpec_p1; XSecPart = iPart_p2; SpecNumTarget = SpecNum1; SpecNumSource = SpecNum2
  END IF
  ! Using the kinetic energy of the particle (as is described in Vahedi1995 and Birdsall1991)
  VeloSquare = DOT_PRODUCT(PartState(4:6,XSecPart),PartState(4:6,XSecPart))
  CollEnergy = 0.5 * Species(XSecSpec)%MassIC * VeloSquare
  nVib = SIZE(SpecXSec(iCase)%VibMode)
  DO iVib = 1, nVib
    ! Calculate the relaxation probability per vibrational mode
    SpecXSec(iCase)%VibMode(iVib)%Prob = (1. - EXP(-SQRT(VeloSquare) * InterpolateCrossSection_Vib(iCase,iVib,CollEnergy) &
                                                    * SpecNumTarget * MacroParticleFactor / Volume * dtCell))
    IF(BGGas%BackgroundSpecies(targetSpec)) THEN
      ! Correct the collision probability in the case of the second species being a background species as the number of pairs
      ! is determined based on the species fraction
      SpecXSec(iCase)%VibMode(iVib)%Prob = SpecXSec(iCase)%VibMode(iVib)%Prob / BGGas%SpeciesFraction(BGGas%MapSpecToBGSpec(targetSpec))
    ELSE
      SpecXSec(iCase)%VibMode(iVib)%Prob = SpecXSec(iCase)%VibMode(iVib)%Prob * SpecNumSource / CollInf%Coll_CaseNum(iCase)
    END IF
    SpecXSec(iCase)%VibProb = SpecXSec(iCase)%VibProb + SpecXSec(iCase)%VibMode(iVib)%Prob
    Coll_pData(iPair)%Prob = Coll_pData(iPair)%Prob + SpecXSec(iCase)%VibMode(iVib)%Prob
  END DO
END IF

END SUBROUTINE XSec_CalcVibRelaxProb


PURE REAL FUNCTION InterpolateCrossSection_Chem(iCase,iPath,CollisionEnergy)
!===================================================================================================================================
!> Interpolate the reaction cross-section data for specific reaction path at the given collision energy
!> Note: Requires the data to be sorted by ascending energy values
!> Assumption: First species given is the particle species, second species input is the background gas species
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars             ,ONLY: SpecXSec
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)            :: iCase                            !< Case index
INTEGER,INTENT(IN)            :: iPath                            !< 
REAL,INTENT(IN)               :: CollisionEnergy                  !< Collision energy in [J]
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iDOF, MaxDOF
!===================================================================================================================================

InterpolateCrossSection_Chem = 0.
MaxDOF = SIZE(SpecXSec(iCase)%ReactionPath(iPath)%XSecData,2)

IF(CollisionEnergy.GT.SpecXSec(iCase)%ReactionPath(iPath)%XSecData(1,MaxDOF)) THEN 
  ! If the collision energy is greater than the maximal value, get the cross-section of the last level and leave routine
  InterpolateCrossSection_Chem = SpecXSec(iCase)%ReactionPath(iPath)%XSecData(2,MaxDOF)
  ! Leave routine
  RETURN
ELSE IF(CollisionEnergy.LE.SpecXSec(iCase)%ReactionPath(iPath)%XSecData(1,1)) THEN
  ! If collision energy is below the minimal value, get the cross-section of the first level and leave routine
  InterpolateCrossSection_Chem = SpecXSec(iCase)%ReactionPath(iPath)%XSecData(2,1)
  ! Leave routine
  RETURN
END IF

DO iDOF = 1, MaxDOF
  ! Check if the stored energy value is above the collision energy
  IF(SpecXSec(iCase)%ReactionPath(iPath)%XSecData(1,iDOF).GT.CollisionEnergy) THEN
    ! Interpolate the cross-section from the data set using the current and the energy level below
    InterpolateCrossSection_Chem = SpecXSec(iCase)%ReactionPath(iPath)%XSecData(2,iDOF-1) &
              + (CollisionEnergy - SpecXSec(iCase)%ReactionPath(iPath)%XSecData(1,iDOF-1)) &
              / (SpecXSec(iCase)%ReactionPath(iPath)%XSecData(1,iDOF) - SpecXSec(iCase)%ReactionPath(iPath)%XSecData(1,iDOF-1)) &
              * (SpecXSec(iCase)%ReactionPath(iPath)%XSecData(2,iDOF) - SpecXSec(iCase)%ReactionPath(iPath)%XSecData(2,iDOF-1))
    ! Leave routine and do not finish DO loop
    RETURN
  END IF
END DO

END FUNCTION InterpolateCrossSection_Chem


SUBROUTINE ReadReacXSec(iCase,iPath)
!===================================================================================================================================
!> Read-in of reaction cross-sections from a given database. Using the effective cross-section database to check whether the
!> group exists. Trying to swap the species indices if dataset not found.
!===================================================================================================================================
! use module
USE MOD_io_hdf5
USE MOD_Globals
USE MOD_Globals_Vars              ,ONLY: ElementaryCharge
USE MOD_DSMC_Vars                 ,ONLY: XSec_Database, SpecXSec, SpecDSMC, ChemReac
USE MOD_HDF5_Input                ,ONLY: DatasetExists
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: iCase, iPath
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=64)                 :: dsetname, groupname, EductPair, dsetname2, ProductPair
INTEGER                           :: err
INTEGER(HSIZE_T), DIMENSION(2)    :: dims,sizeMax
INTEGER(HID_T)                    :: file_id_dsmc                       ! File identifier
INTEGER(HID_T)                    :: group_id                           ! Group identifier
INTEGER(HID_T)                    :: dset_id_dsmc                       ! Dataset identifier
INTEGER(HID_T)                    :: filespace                          ! filespace identifier
INTEGER(SIZE_T)                   :: size                               ! Size of name
INTEGER(HSIZE_T)                  :: iSet                               ! Index
INTEGER                           :: storage, nSets, max_corder
LOGICAL                           :: DataSetFound, GroupFound, ReactionFound
INTEGER                           :: iReac, EductReac(1:3), ProductReac(1:3)
!===================================================================================================================================
iReac = ChemReac%CollCaseInfo(iCase)%ReactionIndex(iPath)
EductReac(1:3) = ChemReac%DefinedReact(iReac,1,1:3)
ProductReac(1:3) = ChemReac%DefinedReact(iReac,2,1:3)

DatasetFound = .FALSE.; GroupFound = .FALSE.

EductPair = TRIM(SpecDSMC(EductReac(1))%Name)//'-'//TRIM(SpecDSMC(EductReac(2))%Name)
SWRITE(UNIT_StdOut,'(A)') 'Read reaction cross section for '//TRIM(EductPair)//' from '//TRIM(XSec_Database)

! Initialize FORTRAN interface.
CALL H5OPEN_F(err)
! Open the file.
CALL H5FOPEN_F (TRIM(XSec_Database), H5F_ACC_RDONLY_F, file_id_dsmc, err)

! Check if the species pair group exists
CALL H5LEXISTS_F(file_id_dsmc, TRIM(EductPair), GroupFound, err)
IF(.NOT.GroupFound) THEN
  ! Try to swap the species names
  EductPair = TRIM(SpecDSMC(EductReac(2))%Name)//'-'//TRIM(SpecDSMC(EductReac(1))%Name)
  CALL H5LEXISTS_F(file_id_dsmc, TRIM(EductPair), GroupFound, err)
  IF(.NOT.GroupFound) THEN
    CALL abort(__STAMP__,&
      'No reaction cross sections found in database for reaction number:',iReac)
  END IF
END IF

groupname = TRIM('/'//TRIM(EductPair)//'/REACTION/')
CALL H5LEXISTS_F(file_id_dsmc, groupname, GroupFound, err)
IF(.NOT.GroupFound) THEN
  CALL abort(__STAMP__,&
    'No reaction cross sections found in database for reaction number:',iReac)
END IF

CALL H5GOPEN_F(file_id_dsmc,TRIM(groupname), group_id, err)
CALL H5Gget_info_f(group_id, storage, nSets,max_corder, err)
IF(nSets.EQ.0) THEN
  CALL abort(__STAMP__,&
    'No reaction cross sections found in database for reaction number:',iReac)
END IF

SELECT CASE(COUNT(ProductReac.GT.0))
CASE(2)
  ProductPair = TRIM(SpecDSMC(ProductReac(1))%Name)//'-'//TRIM(SpecDSMC(ProductReac(2))%Name)
CASE(3)
  ProductPair = TRIM(SpecDSMC(ProductReac(1))%Name)//'-'//TRIM(SpecDSMC(ProductReac(2))%Name)//'-'//TRIM(SpecDSMC(ProductReac(3))%Name)
CASE DEFAULT
  CALL abort(__STAMP__,&
      'Number of products is not supported yet! Reaction number:', iReac)
END SELECT

ReactionFound = .FALSE.

DO iSet = 0, nSets-1
  ! Get name and size of name
  CALL H5Lget_name_by_idx_f(group_id, ".", H5_INDEX_NAME_F, H5_ITER_INC_F, iSet, dsetname, err, size)
  dsetname2 = TRIM(groupname)//TRIM(dsetname)
  ! Look for the correct reaction path
  IF(TRIM(ProductPair).EQ.TRIM(dsetname)) THEN
    ! Open the dataset.
    CALL H5DOPEN_F(file_id_dsmc, dsetname2, dset_id_dsmc, err)
    ! Get the file space of the dataset.
    CALL H5DGET_SPACE_F(dset_id_dsmc, FileSpace, err)
    ! get size
    CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, dims, SizeMax, err)
    ALLOCATE(SpecXSec(iCase)%ReactionPath(iPath)%XSecData(dims(1),dims(2)))
    ! read data
    CALL H5DREAD_F(dset_id_dsmc, H5T_NATIVE_DOUBLE, SpecXSec(iCase)%ReactionPath(iPath)%XSecData, dims, err)
    ReactionFound = .TRUE.
    ! stop looking for other reaction paths
    EXIT
  END IF
END DO

IF(.NOT.ReactionFound) CALL abort(__STAMP__,&
    'No reaction cross-section data found for reaction number:', iReac)

! Store the energy value in J (read-in was in eV)
SpecXSec(iCase)%ReactionPath(iPath)%XSecData(1,:) = SpecXSec(iCase)%ReactionPath(iPath)%XSecData(1,:) * ElementaryCharge

! Close the group
CALL H5GCLOSE_F(group_id,err)
! Close the file.
CALL H5FCLOSE_F(file_id_dsmc, err)
! Close FORTRAN interface.
CALL H5CLOSE_F(err)

END SUBROUTINE ReadReacXSec


SUBROUTINE XSec_CalcReactionProb(iPair,iCase)
!===================================================================================================================================
!> Calculate the collision probability if collision cross-section data is used. Can be utilized in combination with the regular
!> DSMC collision calculation probability.
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars             ,ONLY: SpecDSMC, Coll_pData, CollInf, BGGas, ChemReac, RadialWeighting, DSMC, PartStateIntEn
USE MOD_Particle_Vars         ,ONLY: PartSpecies, Species, VarTimeStep
USE MOD_TimeDisc_Vars         ,ONLY: dt
USE MOD_part_tools            ,ONLY: GetParticleWeight
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)            :: iPair, iCase
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPath, ReacTest, EductReac(1:3), ProductReac(1:3), bggSpec, partSpec, ReactInx(1:3), nPair
INTEGER                       :: NumWeightEduct, NumWeightProd
REAL                          :: EZeroPoint_Prod, dtCell, Weight1, Weight2, Weight3, ReducedMass, ReducedMassUnweighted, WeightProd
REAL                          :: EZeroPoint_Educt
!===================================================================================================================================
WeightProd = 0.; ReactInx = 0
nPair = SIZE(Coll_pData)
NumWeightProd = 2

DO iPath = 1, ChemReac%CollCaseInfo(iCase)%NumOfReactionPaths
  ReacTest = ChemReac%CollCaseInfo(iCase)%ReactionIndex(iPath)
  IF(ChemReac%XSec_Procedure(ReacTest)) THEN
    EductReac(1:3) = ChemReac%DefinedReact(ReacTest,1,1:3)
    ProductReac(1:3) = ChemReac%DefinedReact(ReacTest,2,1:3)
    IF (ChemReac%DefinedReact(ReacTest,1,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
      ReactInx(1) = Coll_pData(iPair)%iPart_p1
      ReactInx(2) = Coll_pData(iPair)%iPart_p2
    ELSE
      ReactInx(1) = Coll_pData(iPair)%iPart_p2
      ReactInx(2) = Coll_pData(iPair)%iPart_p1
    END IF

    IF(EductReac(3).NE.0) THEN
      IF(ChemReac%RecombParticle.EQ.0) THEN
        IF(iPair.LT.(nPair - ChemReac%nPairForRec)) THEN
          ChemReac%LastPairForRec = nPair - ChemReac%nPairForRec
          ReactInx(3) = Coll_pData(ChemReac%LastPairForRec)%iPart_p1
          ChemReac%RecombParticle = ReactInx(3)
        ELSE
          ReactInx(3) = 0
        END IF
      ELSE
        ReactInx(3) = ChemReac%RecombParticle
      END IF
      IF(ReactInx(3).GT.0) THEN
        NumWeightEduct = 3.
        NumWeightProd = 2.
        EductReac(3) = PartSpecies(ReactInx(3))
        ! Save the new reaction index (depending on the third partner) for the case to be used in DSMC_Chemistry
        ChemReac%CollCaseInfo(iCase)%ReactionIndex(iPath) = ChemReac%ReactNumRecomb(EductReac(1), EductReac(2), EductReac(3))
      ELSE
        ! If no third collision partner can be found e.g. last (available) pair, set reaction probability to zero and leave routine
        ChemReac%CollCaseInfo(iCase)%ReactionProb(iPath) = 0.
        RETURN
      END IF
    END IF

    Weight1 = GetParticleWeight(ReactInx(1))
    Weight2 = GetParticleWeight(ReactInx(2))
    IF(EductReac(3).NE.0) Weight3 = GetParticleWeight(ReactInx(3))

    EZeroPoint_Educt = 0.0; EZeroPoint_Prod = 0.0
    ! Testing if the first reacting particle is an atom or molecule, if molecule: is it polyatomic?
    IF((SpecDSMC(EductReac(1))%InterID.EQ.2).OR.(SpecDSMC(EductReac(1))%InterID.EQ.20)) THEN
      EZeroPoint_Educt = EZeroPoint_Educt + SpecDSMC(EductReac(1))%EZeroPoint * Weight1
    END IF
    !---------------------------------------------------------------------------------------------------------------------------------
    ! Testing if the second particle is an atom or molecule, if molecule: is it polyatomic?
    IF((SpecDSMC(EductReac(2))%InterID.EQ.2).OR.(SpecDSMC(EductReac(2))%InterID.EQ.20)) THEN
      EZeroPoint_Educt = EZeroPoint_Educt + SpecDSMC(EductReac(2))%EZeroPoint * Weight2
    END IF
    !---------------------------------------------------------------------------------------------------------------------------------
    IF(EductReac(3).NE.0) THEN
      ! Testing if the third particle is an atom or molecule, if molecule: is it polyatomic?
      IF((SpecDSMC(EductReac(3))%InterID.EQ.2).OR.(SpecDSMC(EductReac(3))%InterID.EQ.20)) THEN
        EZeroPoint_Educt = EZeroPoint_Educt + SpecDSMC(EductReac(3))%EZeroPoint * Weight3
      END IF
    END IF

    IF (RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
      ReducedMass = (Species(EductReac(1))%MassIC *Weight1  * Species(EductReac(2))%MassIC * Weight2) &
        / (Species(EductReac(1))%MassIC * Weight1+ Species(EductReac(2))%MassIC * Weight2)
      ReducedMassUnweighted = ReducedMass * 2./(Weight1+Weight2)
    ELSE
      ReducedMass = CollInf%MassRed(Coll_pData(iPair)%PairType)
      ReducedMassUnweighted = CollInf%MassRed(Coll_pData(iPair)%PairType)
    END IF

    ! Testing if the first produced particle is an atom or molecule, if molecule: is it polyatomic?
    IF((SpecDSMC(ProductReac(1))%InterID.EQ.2).OR.(SpecDSMC(ProductReac(1))%InterID.EQ.20)) THEN
      EZeroPoint_Prod = EZeroPoint_Prod + SpecDSMC(ProductReac(1))%EZeroPoint * Weight1
    END IF
    ! Testing if the second produced particle is an atom or molecule, if molecule: is it polyatomic?
    IF((SpecDSMC(ProductReac(2))%InterID.EQ.2).OR.(SpecDSMC(ProductReac(2))%InterID.EQ.20)) THEN
      EZeroPoint_Prod = EZeroPoint_Prod + SpecDSMC(ProductReac(2))%EZeroPoint * Weight2
    END IF
    IF(ProductReac(3).NE.0) THEN
      NumWeightProd = 3.
      IF(EductReac(3).NE.0) THEN
        WeightProd = Weight3
      ELSE
        WeightProd = Weight1
      END IF
      IF((SpecDSMC(ProductReac(3))%InterID.EQ.2).OR.(SpecDSMC(ProductReac(3))%InterID.EQ.20)) THEN
        EZeroPoint_Prod = EZeroPoint_Prod + SpecDSMC(ProductReac(3))%EZeroPoint*WeightProd
      END IF
    END IF

    ! Select the background species as the target cloud
    IF(BGGas%BackgroundSpecies(EductReac(1))) THEN
      bggSpec = BGGas%MapSpecToBGSpec(EductReac(1))
      partSpec = EductReac(2)
    ELSE
      bggSpec = BGGas%MapSpecToBGSpec(EductReac(2))
      partSpec = EductReac(1)
    END IF
    IF (VarTimeStep%UseVariableTimeStep) THEN
      dtCell = dt * (VarTimeStep%ParticleTimeStep(ReactInx(1)) + VarTimeStep%ParticleTimeStep(ReactInx(2)))*0.5
    ELSE
      dtCell = dt
    END IF
    Coll_pData(iPair)%Ec = 0.5 * ReducedMass * Coll_pData(iPair)%CRela2 &
                + (PartStateIntEn(1,ReactInx(1)) + PartStateIntEn(2,ReactInx(1))) * Weight1 &
                + (PartStateIntEn(1,ReactInx(2)) + PartStateIntEn(2,ReactInx(2))) * Weight2
    IF (DSMC%ElectronicModel) THEN
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(3,ReactInx(1))*Weight1 + PartStateIntEn(3,ReactInx(2))*Weight2
    END IF
    ! Calculate the reaction probability
    IF(((Coll_pData(iPair)%Ec-EZeroPoint_Prod).GE.(-1./NumWeightProd*(Weight1+Weight2+WeightProd)*ChemReac%EForm(ReacTest)))) THEN
      ChemReac%CollCaseInfo(iCase)%ReactionProb(iPath) = (1. - EXP(-SQRT(Coll_pData(iPair)%CRela2) * dtCell &
        * InterpolateCrossSection_Chem(iCase,iPath,Coll_pData(iPair)%Ec-EZeroPoint_Educt) * BGGas%NumberDensity(bggSpec) ))
      Coll_pData(iPair)%Prob = Coll_pData(iPair)%Prob + ChemReac%CollCaseInfo(iCase)%ReactionProb(iPath)
    ELSE
      ChemReac%CollCaseInfo(iCase)%ReactionProb(iPath) = 0.
    END IF
    ! Calculation of reaction rate coefficient
#if (PP_TimeDiscMethod==42)
    IF (.NOT.DSMC%ReservoirRateStatistic) THEN
      ChemReac%NumReac(ReacTest) = ChemReac%NumReac(ReacTest) + ChemReac%CollCaseInfo(iCase)%ReactionProb(iPath)
      ChemReac%ReacCount(ReacTest) = ChemReac%ReacCount(ReacTest) + 1
    END IF
#endif
  END IF
END DO

END SUBROUTINE XSec_CalcReactionProb

END MODULE MOD_DSMC_SpecXSec