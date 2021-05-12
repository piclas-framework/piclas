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
REAL          :: TotalProb, VibCrossSection
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
SpecXSec(:)%CollXSec_Effective = .FALSE.
SpecXSec(:)%SpeciesToRelax = 0

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
    IF(SpecDSMC(iSpec)%UseCollXSec.OR.SpecDSMC(jSpec)%UseCollXSec) CALL ReadCollXSec(iCase, iSpec, jSpec)
    ! Check if both species were given the UseCollXSec flag and store the energy value in Joule
    IF(SpecXSec(iCase)%UseCollXSec) THEN
      IF(SpecDSMC(iSpec)%UseCollXSec.AND.SpecDSMC(jSpec)%UseCollXSec) THEN
        CALL abort(&
          __STAMP__&
          ,'ERROR: Both species defined to use collisional cross-section, define only the source species with UseCollXSec!')
      END IF
      ! Store the energy value in J (read-in was in eV)
      SpecXSec(iCase)%CollXSecData(1,:) = SpecXSec(iCase)%CollXSecData(1,:) * ElementaryCharge
    END IF
    ! Read-in vibrational cross sections
    IF(SpecDSMC(iSpec)%UseVibXSec.OR.SpecDSMC(jSpec)%UseVibXSec) CALL ReadVibXSec(iCase, iSpec, jSpec)
    ! Vibrational relaxation probabilities: Interpolate and store the probability at the collision cross-section levels
    IF(SpecXSec(iCase)%UseVibXSec) THEN
      IF(SpecDSMC(iSpec)%UseVibXSec.AND.SpecDSMC(jSpec)%UseVibXSec) THEN
        CALL abort(&
          __STAMP__&
          ,'ERROR: Both species defined to use vib. cross-section, define only the source species with UseVibXSec!')
      END IF
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
      nVib = SIZE(SpecXSec(iCase)%VibMode)
      DO iVib = 1, nVib
        ! Store the energy value in J (read-in was in eV)
        SpecXSec(iCase)%VibMode(iVib)%XSecData(1,:) = SpecXSec(iCase)%VibMode(iVib)%XSecData(1,:) * ElementaryCharge
      END DO
      IF(SpecXSec(iCase)%UseCollXSec) THEN
        ! Collision cross-sections are available
        MaxDim = SIZE(SpecXSec(iCase)%CollXSecData,2)
        ALLOCATE(SpecXSec(iCase)%VibXSecData(1:2,1:MaxDim))
        ! Using the same energy intervals as for the collision cross-sections
        SpecXSec(iCase)%VibXSecData(1,:) = SpecXSec(iCase)%CollXSecData(1,:)
        SpecXSec(iCase)%VibXSecData(2,:) = 0.
        ! Interpolate the vibrational cross section at the energy levels of the collision collision cross section and sum-up the
        ! vibrational probability (vibrational cross-section divided by the effective)
        DO iStep = 1, MaxDim
          DO iVib = 1, nVib
            VibCrossSection = InterpolateCrossSection_Vib(iCase,iVib,SpecXSec(iCase)%CollXSecData(1,iStep))
            SpecXSec(iCase)%VibXSecData(2,iStep) = SpecXSec(iCase)%VibXSecData(2,iStep) + VibCrossSection
            ! When no effective cross-section is available, the vibrational cross-section has to be added to the collisional
            IF(.NOT.SpecXSec(iCase)%CollXSec_Effective) SpecXSec(iCase)%CollXSecData(2,iStep) &
                                                        = SpecXSec(iCase)%CollXSecData(2,iStep) + VibCrossSection
          END DO
        END DO
      END IF    ! SpecXSec(iCase)%UseCollXSec
    END IF      ! SpecXSec(iCase)%UseVibXSec
    IF(SpecXSec(iCase)%UseCollXSec) THEN
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
CHARACTER(LEN=64)                 :: dsetname, spec_pair
INTEGER                           :: err
INTEGER(HSIZE_T), DIMENSION(2)    :: dims,sizeMax
INTEGER(HID_T)                    :: file_id_dsmc                       ! File identifier
INTEGER(HID_T)                    :: dset_id_dsmc                       ! Dataset identifier
INTEGER(HID_T)                    :: filespace                          ! filespace identifier
LOGICAL                           :: DatasetFound
!===================================================================================================================================
spec_pair = TRIM(SpecDSMC(jSpec)%Name)//'-'//TRIM(SpecDSMC(iSpec)%Name)

DatasetFound = .FALSE.

! Initialize FORTRAN interface.
CALL H5OPEN_F(err)

! Check if file exists
IF(.NOT.FILEEXISTS(XSec_Database)) THEN
  CALL abort(__STAMP__,'ERROR: Database '//TRIM(XSec_Database)//' does not exist.')
END IF

! Open the file.
CALL H5FOPEN_F (TRIM(XSec_Database), H5F_ACC_RDONLY_F, file_id_dsmc, err)

! Check if the species pair group exists
CALL H5LEXISTS_F(file_id_dsmc, TRIM(spec_pair), DatasetFound, err)
IF(.NOT.DatasetFound) THEN
  ! Try to swap the species names
  spec_pair = TRIM(SpecDSMC(iSpec)%Name)//'-'//TRIM(SpecDSMC(jSpec)%Name)
  CALL H5LEXISTS_F(file_id_dsmc, TRIM(spec_pair), DatasetFound, err)
  IF(.NOT.DatasetFound) THEN
    SWRITE(UNIT_StdOut,'(A)') TRIM(spec_pair)//': No data set found. Using standard collision modelling.'
    RETURN
  END IF
END IF

dsetname = TRIM('/'//TRIM(spec_pair)//'/EFFECTIVE')
CALL DatasetExists(File_ID_DSMC,TRIM(dsetname),DatasetFound)

IF(DatasetFound) THEN
  SpecXSec(iCase)%UseCollXSec = .TRUE.
  SpecXSec(iCase)%CollXSec_Effective = .TRUE.
  CALL DatasetExists(File_ID_DSMC,TRIM('/'//TRIM(spec_pair)//'/ELASTIC'),DatasetFound)
  IF(DatasetFound) CALL abort(__STAMP__,'ERROR: Please provide either elastic or effective collision cross-section data '//&
                                             & 'for '//TRIM(spec_pair)//'.')
  SWRITE(UNIT_StdOut,'(A)') TRIM(spec_pair)//': Found EFFECTIVE collision cross section.'
ELSE
  dsetname = TRIM('/'//TRIM(spec_pair)//'/ELASTIC')
  CALL DatasetExists(File_ID_DSMC,TRIM(dsetname),DatasetFound)
  IF(DatasetFound) THEN
    SpecXSec(iCase)%UseCollXSec = .TRUE.
    SpecXSec(iCase)%CollXSec_Effective = .FALSE.
  SWRITE(UNIT_StdOut,'(A)') TRIM(spec_pair)//': Found ELASTIC collision cross section.'
  ELSE
    SWRITE(UNIT_StdOut,'(A)') TRIM(spec_pair)//': No data set found. Using standard collision modelling.'
    RETURN
  END IF
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
!> Read-in of vibrational cross-sections from a database. Dataset name is composed of SpeciesName-SpeciesName (e.g. Ar-electron)
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
    SWRITE(UNIT_StdOut,'(A)') TRIM(spec_pair)//': No vibrational excitation cross sections found, using constant read-in values.'
    RETURN
  END IF
END IF

! Check if the vibrational cross-section group exists
groupname = TRIM(spec_pair)//'/VIBRATION/'
CALL H5LEXISTS_F(file_id_dsmc, TRIM(groupname), GroupFound, err)
IF(.NOT.GroupFound) THEN
  SWRITE(UNIT_StdOut,'(A)') TRIM(spec_pair)//': No vibrational excitation cross sections found, using constant read-in values.'
  RETURN
END IF

IF(GroupFound) THEN
  CALL H5GOPEN_F(file_id_dsmc,TRIM(groupname), group_id, err)
  call H5Gget_info_f(group_id, storage, nVib,max_corder, err)
  ! If cross-section data is found, set the corresponding flag
  IF(nVib.GT.0) THEN
    SWRITE(UNIT_StdOut,'(A,I3,A)') TRIM(spec_pair)//': Found ', nVib,' vibrational excitation cross section(s).'
    SpecXSec(iCase)%UseVibXSec = .TRUE.
    nVar = 3
  ELSE
    SWRITE(UNIT_StdOut,'(A)') TRIM(spec_pair)//': No vibrational excitation cross sections found, using constant read-in values.'
  END IF
ELSE
  SWRITE(UNIT_StdOut,'(A)') TRIM(spec_pair)//': No vibrational excitation cross sections found, using constant read-in values.'
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
USE MOD_Particle_Vars         ,ONLY: Species
USE MOD_TimeDisc_Vars         ,ONLY: ManualTimeStep
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
  IF(SpecXSec(iCase)%VibMode(iVib)%XSecData(1,iDOF).GE.CollisionEnergy) THEN
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


SUBROUTINE XSec_CalcCollisionProb(iPair,SpecNum1,SpecNum2,CollCaseNum,MacroParticleFactor,Volume,dtCell)
!===================================================================================================================================
!> Calculate the collision probability if collision cross-section data is used. Can be utilized in combination with the regular
!> DSMC collision calculation probability.
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars             ,ONLY: SpecXSec, SpecDSMC, Coll_pData, CollInf, BGGas, XSec_NullCollision, RadialWeighting
USE MOD_Particle_Vars         ,ONLY: PartSpecies, Species, VarTimeStep
USE MOD_part_tools            ,ONLY: GetParticleWeight
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)            :: iPair
REAL,INTENT(IN)               :: SpecNum1, SpecNum2, CollCaseNum, MacroParticleFactor, Volume, dtCell
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: CollEnergy, SpecNumTarget, SpecNumSource, Weight1, Weight2, ReducedMass
INTEGER                       :: targetSpec, iPart_p1, iPart_p2, iSpec_p1, iSpec_p2, iCase
!===================================================================================================================================

iPart_p1 = Coll_pData(iPair)%iPart_p1; iPart_p2 = Coll_pData(iPair)%iPart_p2
iSpec_p1 = PartSpecies(iPart_p1);      iSpec_p2 = PartSpecies(iPart_p2)
iCase = CollInf%Coll_Case(iSpec_p1,iSpec_p2)

IF(SpecXSec(iCase)%UseCollXSec) THEN
  IF(SpecDSMC(iSpec_p1)%UseCollXSec) THEN
    targetSpec = iSpec_p2; SpecNumTarget = SpecNum2; SpecNumSource = SpecNum1
  ELSE
    targetSpec = iSpec_p1; SpecNumTarget = SpecNum1; SpecNumSource = SpecNum2
  END IF
  Weight1 = GetParticleWeight(iPart_p1)
  Weight2 = GetParticleWeight(iPart_p2)
  IF (RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
    ReducedMass = (Species(iSpec_p1)%MassIC *Weight1  * Species(iSpec_p2)%MassIC * Weight2) &
      / (Species(iSpec_p1)%MassIC * Weight1+ Species(iSpec_p2)%MassIC * Weight2)
  ELSE
    ReducedMass = CollInf%MassRed(Coll_pData(iPair)%PairType)
  END IF
  ! Using the relative kinetic energy of the particle
  CollEnergy = 0.5 * ReducedMass * Coll_pData(iPair)%CRela2
  ! Calculate the collision probability
  SpecXSec(iCase)%CrossSection = InterpolateCrossSection(iCase,CollEnergy)
  Coll_pData(iPair)%Prob = (1. - EXP(-SQRT(Coll_pData(iPair)%CRela2) * SpecXSec(iCase)%CrossSection * SpecNumTarget * MacroParticleFactor &
                                        / Volume * dtCell))
  ! Correction for conditional probabilities in case of MCC
  IF(BGGas%BackgroundSpecies(targetSpec)) THEN
    ! Correct the collision probability in the case of the second species being a background species as the number of pairs
    ! is either determined based on the null collision probability or on the species fraction
    IF(XSec_NullCollision) THEN
      Coll_pData(iPair)%Prob = Coll_pData(iPair)%Prob / SpecXSec(iCase)%ProbNull
    ELSE
      Coll_pData(iPair)%Prob = Coll_pData(iPair)%Prob / BGGas%SpeciesFraction(BGGas%MapSpecToBGSpec(targetSpec))
    END IF
  ELSE
    ! Using cross-sectional probabilities without background gas
    Coll_pData(iPair)%Prob = Coll_pData(iPair)%Prob * SpecNumSource / CollInf%Coll_CaseNum(iCase)
  END IF
ELSE
  Coll_pData(iPair)%Prob = SpecNum1*SpecNum2/(1 + CollInf%KronDelta(Coll_pData(iPair)%PairType))  &
          * CollInf%Cab(Coll_pData(iPair)%PairType)                           & ! Cab species comb fac
          * MacroParticleFactor / CollCaseNum                                                     &
          * Coll_pData(iPair)%CRela2 ** (0.5-CollInf%omega(iSpec_p1,iSpec_p2)) &
          * dtCell / Volume
END IF

END SUBROUTINE XSec_CalcCollisionProb


SUBROUTINE XSec_CalcVibRelaxProb(iPair,SpecNum1,SpecNum2,MacroParticleFactor,Volume,dtCell)
!===================================================================================================================================
!> Calculate the relaxation probability using cross-section data.
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars             ,ONLY: SpecXSec, SpecDSMC, Coll_pData, CollInf, BGGas, RadialWeighting
USE MOD_Particle_Vars         ,ONLY: PartSpecies, Species, VarTimeStep
USE MOD_part_tools            ,ONLY: GetParticleWeight
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)            :: iPair
REAL,INTENT(IN),OPTIONAL      :: SpecNum1, SpecNum2, MacroParticleFactor, Volume, dtCell
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: CollEnergy, SpecNumTarget, SpecNumSource, Weight1, Weight2, ReducedMass, SumVibCrossSection
INTEGER                       :: targetSpec, iPart_p1, iPart_p2, iSpec_p1, iSpec_p2, iCase, iVib, nVib
!===================================================================================================================================

iPart_p1 = Coll_pData(iPair)%iPart_p1; iPart_p2 = Coll_pData(iPair)%iPart_p2
iSpec_p1 = PartSpecies(iPart_p1);      iSpec_p2 = PartSpecies(iPart_p2)
iCase = CollInf%Coll_Case(iSpec_p1,iSpec_p2)
Weight1 = GetParticleWeight(iPart_p1)
Weight2 = GetParticleWeight(iPart_p2)
SpecXSec(iCase)%VibProb = 0.
SumVibCrossSection = 0.

IF (RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
  ReducedMass = (Species(iSpec_p1)%MassIC *Weight1  * Species(iSpec_p2)%MassIC * Weight2) &
    / (Species(iSpec_p1)%MassIC * Weight1+ Species(iSpec_p2)%MassIC * Weight2)
ELSE
  ReducedMass = CollInf%MassRed(Coll_pData(iPair)%PairType)
END IF
! Using the relative translational energy of the pair
CollEnergy = 0.5 * ReducedMass * Coll_pData(iPair)%CRela2
! Calculate the total vibrational cross-section
nVib = SIZE(SpecXSec(iCase)%VibMode)
DO iVib = 1, nVib
  SumVibCrossSection = SumVibCrossSection + InterpolateCrossSection_Vib(iCase,iVib,CollEnergy)
END DO

IF(SpecXSec(iCase)%UseCollXSec) THEN
  SpecXSec(iCase)%VibProb = SumVibCrossSection
ELSE
  IF(SpecDSMC(iSpec_p1)%UseVibXSec) THEN
    targetSpec = iSpec_p2; SpecNumTarget = SpecNum2; SpecNumSource = SpecNum1
  ELSE
    targetSpec = iSpec_p1; SpecNumTarget = SpecNum1; SpecNumSource = SpecNum2
  END IF
  ! Calculate the total vibrational relaxation probability
  SpecXSec(iCase)%VibProb = (1. - EXP(-SQRT(Coll_pData(iPair)%CRela2) * SpecNumTarget * MacroParticleFactor / Volume &
                                                  * SumVibCrossSection * dtCell))
  IF(BGGas%BackgroundSpecies(targetSpec)) THEN
    ! Correct the collision probability in the case of the second species being a background species as the number of pairs
    ! is determined based on the species fraction
    SpecXSec(iCase)%VibProb = SpecXSec(iCase)%VibProb / BGGas%SpeciesFraction(BGGas%MapSpecToBGSpec(targetSpec))
  ELSE
    SpecXSec(iCase)%VibProb = SpecXSec(iCase)%VibProb * SpecNumSource / CollInf%Coll_CaseNum(iCase)
  END IF
END IF

END SUBROUTINE XSec_CalcVibRelaxProb


SUBROUTINE MCC_Chemistry_Init()
!===================================================================================================================================
!> Read-in of the reaction cross-section database and re-calculation of the null collision probability
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_PARTICLE_Vars ,ONLY: nSpecies
USE MOD_DSMC_Vars     ,ONLY: BGGas, SpecXSec, XSec_NullCollision, CollInf, ChemReac
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSpec, jSpec, iCase, iReac
REAL                  :: TotalProb, ReactionCrossSection
INTEGER               :: iStep, MaxDim
INTEGER               :: iPath, NumPaths
!===================================================================================================================================

IF(BGGas%NumberOfSpecies.LE.0) THEN
  CALL abort(__STAMP__,&
    'Chemistry - Error: Cross-section based chemistry without background gas has not been tested yet!')
END IF

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
    ! When no effective cross-section is available, the total cross-section has to be determined
    IF(.NOT.SpecXSec(iCase)%CollXSec_Effective) THEN
      MaxDim = SIZE(SpecXSec(iCase)%CollXSecData,2)
      NumPaths = ChemReac%CollCaseInfo(iCase)%NumOfReactionPaths
      ! Interpolate the reaction cross section at the energy levels of the collision collision cross section
      DO iPath = 1, NumPaths
        DO iStep = 1, MaxDim
          ReactionCrossSection = InterpolateCrossSection_Chem(iCase,iPath,SpecXSec(iCase)%CollXSecData(1,iStep))
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
        TotalProb = TotalProb + SpecXSec(iCase)%ProbNull
        IF(TotalProb.GT.1.0) THEN
          CALL abort(&
          __STAMP__&
          ,'ERROR: Total null collision probability is above unity. Please reduce the time step! Probability is: '&
          ,RealInfoOpt=TotalProb)
        END IF
      END IF
    END DO
  END DO
END IF


END SUBROUTINE MCC_Chemistry_Init


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
CHARACTER(LEN=128)                :: dsetname, groupname, EductPair, dsetname2, ProductPair
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
INTEGER                           :: iReac, EductReac(1:3), ProductReac(1:4)
!===================================================================================================================================
iReac = ChemReac%CollCaseInfo(iCase)%ReactionIndex(iPath)
EductReac(1:3) = ChemReac%Reactants(iReac,1:3)
ProductReac(1:4) = ChemReac%Products(iReac,1:4)

DatasetFound = .FALSE.; GroupFound = .FALSE.

EductPair = TRIM(SpecDSMC(EductReac(1))%Name)//'-'//TRIM(SpecDSMC(EductReac(2))%Name)

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

SWRITE(UNIT_StdOut,'(A,I0)') 'Read cross section for '//TRIM(EductPair)//' from '//TRIM(XSec_Database)//' for reaction # ', iReac

SELECT CASE(COUNT(ProductReac.GT.0))
CASE(2)
  ProductPair = TRIM(SpecDSMC(ProductReac(1))%Name)//'-'//TRIM(SpecDSMC(ProductReac(2))%Name)
CASE(3)
  ProductPair = TRIM(SpecDSMC(ProductReac(1))%Name)//'-'//TRIM(SpecDSMC(ProductReac(2))%Name)//&
                &'-'//TRIM(SpecDSMC(ProductReac(3))%Name)
CASE(4)
  ProductPair = TRIM(SpecDSMC(ProductReac(1))%Name)//'-'//TRIM(SpecDSMC(ProductReac(2))%Name)//&
                &'-'//TRIM(SpecDSMC(ProductReac(3))%Name)//'-'//TRIM(SpecDSMC(ProductReac(4))%Name)
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


SUBROUTINE XSec_CalcReactionProb(iPair,iCase,SpecNum1,SpecNum2,MacroParticleFactor,Volume)
!===================================================================================================================================
!> Calculate the collision probability if collision cross-section data is used. Can be utilized in combination with the regular
!> DSMC collision calculation probability.
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars             ,ONLY: SpecDSMC, Coll_pData, CollInf, BGGas, ChemReac, RadialWeighting, DSMC, PartStateIntEn, SpecXSec
USE MOD_Particle_Vars         ,ONLY: PartSpecies, Species, VarTimeStep
USE MOD_TimeDisc_Vars         ,ONLY: dt
USE MOD_part_tools            ,ONLY: GetParticleWeight
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)            :: iPair, iCase
REAL,INTENT(IN),OPTIONAL      :: SpecNum1, SpecNum2, MacroParticleFactor, Volume
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPath, ReacTest, EductReac(1:3), ProductReac(1:4), ReactInx(1:2), nPair, iProd
INTEGER                       :: NumWeightProd, targetSpec
REAL                          :: EZeroPoint_Prod, dtCell, Weight(1:4), ReducedMass, ReducedMassUnweighted
REAL                          :: EZeroPoint_Educt, SpecNumTarget, SpecNumSource, CrossSection
!===================================================================================================================================
Weight = 0.; ReactInx = 0
nPair = SIZE(Coll_pData)
NumWeightProd = 2

DO iPath = 1, ChemReac%CollCaseInfo(iCase)%NumOfReactionPaths
  ReacTest = ChemReac%CollCaseInfo(iCase)%ReactionIndex(iPath)
  IF(TRIM(ChemReac%ReactModel(ReacTest)).EQ.'XSec') THEN
    EductReac(1:3) = ChemReac%Reactants(ReacTest,1:3); ProductReac(1:4) = ChemReac%Products(ReacTest,1:4)
    IF (EductReac(1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
      ReactInx(1) = Coll_pData(iPair)%iPart_p1; ReactInx(2) = Coll_pData(iPair)%iPart_p2
    ELSE
      ReactInx(1) = Coll_pData(iPair)%iPart_p2; ReactInx(2) = Coll_pData(iPair)%iPart_p1
    END IF

    Weight(1) = GetParticleWeight(ReactInx(1)); Weight(2) = GetParticleWeight(ReactInx(2))

    EZeroPoint_Educt = 0.0; EZeroPoint_Prod = 0.0
    ! Testing if the first reacting particle is an atom or molecule, if molecule: is it polyatomic?
    IF((SpecDSMC(EductReac(1))%InterID.EQ.2).OR.(SpecDSMC(EductReac(1))%InterID.EQ.20)) THEN
      EZeroPoint_Educt = EZeroPoint_Educt + SpecDSMC(EductReac(1))%EZeroPoint * Weight(1)
    END IF
    !---------------------------------------------------------------------------------------------------------------------------------
    ! Testing if the second particle is an atom or molecule, if molecule: is it polyatomic?
    IF((SpecDSMC(EductReac(2))%InterID.EQ.2).OR.(SpecDSMC(EductReac(2))%InterID.EQ.20)) THEN
      EZeroPoint_Educt = EZeroPoint_Educt + SpecDSMC(EductReac(2))%EZeroPoint * Weight(2)
    END IF

    IF (RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
      ReducedMass = (Species(EductReac(1))%MassIC *Weight(1) * Species(EductReac(2))%MassIC * Weight(2)) &
        / (Species(EductReac(1))%MassIC * Weight(1) + Species(EductReac(2))%MassIC * Weight(2))
      ReducedMassUnweighted = ReducedMass * 2./(Weight(1)+Weight(2))
    ELSE
      ReducedMass = CollInf%MassRed(Coll_pData(iPair)%PairType)
      ReducedMassUnweighted = CollInf%MassRed(Coll_pData(iPair)%PairType)
    END IF

    IF(ProductReac(4).NE.0) THEN
      ! 4 Products
      NumWeightProd = 4
      Weight(3:4) = Weight(1)
    ELSE IF(ProductReac(3).NE.0) THEN
      ! 3 Products
      NumWeightProd = 3
      Weight(3) = Weight(1)
    END IF

    DO iProd = 1, NumWeightProd
      IF((SpecDSMC(ProductReac(iProd))%InterID.EQ.2).OR.(SpecDSMC(ProductReac(iProd))%InterID.EQ.20)) THEN
        EZeroPoint_Prod = EZeroPoint_Prod + SpecDSMC(ProductReac(iProd))%EZeroPoint * Weight(iProd)
      END IF
    END DO

    IF (VarTimeStep%UseVariableTimeStep) THEN
      dtCell = dt * (VarTimeStep%ParticleTimeStep(ReactInx(1)) + VarTimeStep%ParticleTimeStep(ReactInx(2)))*0.5
    ELSE
      dtCell = dt
    END IF
    Coll_pData(iPair)%Ec = 0.5 * ReducedMass * Coll_pData(iPair)%CRela2 &
                + (PartStateIntEn(1,ReactInx(1)) + PartStateIntEn(2,ReactInx(1))) * Weight(1) &
                + (PartStateIntEn(1,ReactInx(2)) + PartStateIntEn(2,ReactInx(2))) * Weight(2)
    IF (DSMC%ElectronicModel.GT.0) THEN
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(3,ReactInx(1))*Weight(1) &
                                                  + PartStateIntEn(3,ReactInx(2))*Weight(2)
    END IF
    ! Check first if sufficient energy is available for the products after the reaction
    IF(((Coll_pData(iPair)%Ec-EZeroPoint_Prod).GE.(-ChemReac%EForm(ReacTest)*SUM(Weight)/NumWeightProd))) THEN
      CrossSection = InterpolateCrossSection_Chem(iCase,iPath,Coll_pData(iPair)%Ec-EZeroPoint_Educt)
      IF(SpecXSec(iCase)%UseCollXSec) THEN
        ! Interpolate the reaction cross-section
        ChemReac%CollCaseInfo(iCase)%ReactionProb(iPath) = CrossSection
      ELSE
        ! Select the background species as the target cloud
        IF(BGGas%BackgroundSpecies(PartSpecies(Coll_pData(iPair)%iPart_p1))) THEN
          targetSpec = PartSpecies(Coll_pData(iPair)%iPart_p1); SpecNumTarget = SpecNum1; SpecNumSource = SpecNum2
        ELSE
          targetSpec = PartSpecies(Coll_pData(iPair)%iPart_p2); SpecNumTarget = SpecNum2; SpecNumSource = SpecNum1
        END IF
        ! Calculate the reaction probability
        ChemReac%CollCaseInfo(iCase)%ReactionProb(iPath) = (1. - EXP(-SQRT(Coll_pData(iPair)%CRela2) * dtCell * SpecNumTarget &
                                                                      * MacroParticleFactor / Volume * CrossSection))
        IF(BGGas%BackgroundSpecies(targetSpec)) THEN
        ! Correct the reaction probability in the case of the second species being a background species as the number of pairs
        ! is based on the species fraction
          ChemReac%CollCaseInfo(iCase)%ReactionProb(iPath) = ChemReac%CollCaseInfo(iCase)%ReactionProb(iPath) &
                                                              / BGGas%SpeciesFraction(BGGas%MapSpecToBGSpec(targetSpec))
        ELSE
          ChemReac%CollCaseInfo(iCase)%ReactionProb(iPath) = ChemReac%CollCaseInfo(iCase)%ReactionProb(iPath) &
                                                              * SpecNumSource / CollInf%Coll_CaseNum(iCase)
        END IF
      END IF
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
