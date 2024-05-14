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

MODULE MOD_MCC_XSec
!===================================================================================================================================
! Contains the Argon Ionization
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

PUBLIC :: ReadCollXSec, ReadVibXSec, ReadElecXSec, ReadReacXSec, ReadReacPhotonXSec, ReadReacPhotonSpectrum
PUBLIC :: InterpolateCrossSection
PUBLIC :: XSec_CalcCollisionProb, XSec_CalcVibRelaxProb, XSec_CalcElecRelaxProb, XSec_CalcReactionProb
PUBLIC :: XSec_ElectronicRelaxation
!===================================================================================================================================

CONTAINS

SUBROUTINE ReadCollXSec(iCase,iSpec,jSpec)
!===================================================================================================================================
!> Read-in of collision cross-sections from a given database. Dataset name is composed of SpeciesName-SpeciesName (e.g. Ar-electron)
!> Trying to swap the species indices if dataset not found.
!===================================================================================================================================
! use module
USE MOD_io_hdf5
USE MOD_Globals
USE MOD_Particle_Vars    ,ONLY: SpeciesDatabase, Species
USE MOD_MCC_Vars         ,ONLY: XSec_Database, SpecXSec
USE MOD_HDF5_Input       ,ONLY: DatasetExists
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
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
CHARACTER(LEN=265)                :: XSecDatabaseName
INTEGER                           :: err
INTEGER(HSIZE_T), DIMENSION(2)    :: dims,sizeMax
INTEGER(HID_T)                    :: file_id_dsmc                       ! File identifier
INTEGER(HID_T)                    :: dset_id_dsmc                       ! Dataset identifier
INTEGER(HID_T)                    :: filespace                          ! filespace identifier
LOGICAL                           :: DatasetFound
!===================================================================================================================================
spec_pair = TRIM(Species(jSpec)%Name)//'-'//TRIM(Species(iSpec)%Name)

DatasetFound = .FALSE.

! Initialize FORTRAN interface.
CALL H5OPEN_F(err)

! Check if file exists
IF((.NOT.FILEEXISTS(XSec_Database)).AND.(SpeciesDatabase.EQ.'none')) THEN
  CALL abort(__STAMP__,'ERROR: Database ['//TRIM(XSec_Database)//'] does not exist.')
END IF

IF(SpeciesDatabase.EQ.'none') THEN
  XSecDatabaseName = TRIM(XSec_Database)
ELSE 
  XSecDatabaseName = TRIM(SpeciesDatabase)
  spec_pair = TRIM('/Cross-Sections/'//Species(jSpec)%Name)//'-'//TRIM(Species(iSpec)%Name)
END IF

! Open the file.
CALL H5FOPEN_F (TRIM(XSecDatabaseName), H5F_ACC_RDONLY_F, file_id_dsmc, err)

! Check if the species pair group exists
CALL H5LEXISTS_F(file_id_dsmc, TRIM(spec_pair), DatasetFound, err)
IF(.NOT.DatasetFound) THEN
  ! Try to swap the species names
  IF (SpeciesDatabase.EQ.'none') THEN
    spec_pair = TRIM(Species(iSpec)%Name)//'-'//TRIM(Species(jSpec)%Name)
  ELSE 
    spec_pair = TRIM('Cross-Sections/'//Species(iSpec)%Name)//'-'//TRIM(Species(jSpec)%Name)
  END IF
  CALL H5LEXISTS_F(file_id_dsmc, TRIM(spec_pair), DatasetFound, err)
  IF(.NOT.DatasetFound) THEN
    LBWRITE(UNIT_StdOut,'(A)') ' | '//TRIM(spec_pair)//': No data set found. Using standard collision modelling.'
    RETURN
  END IF
END IF

IF (SpeciesDatabase.EQ.'none') THEN
  dsetname = TRIM('/'//TRIM(spec_pair)//'/EFFECTIVE')
ELSE 
  dsetname = TRIM(TRIM(spec_pair)//'/EFFECTIVE')
END IF
CALL DatasetExists(File_ID_DSMC,TRIM(dsetname),DatasetFound)

IF(DatasetFound) THEN
  SpecXSec(iCase)%UseCollXSec = .TRUE.
  SpecXSec(iCase)%CollXSec_Effective = .TRUE.
  IF (SpeciesDatabase.EQ.'none') THEN
    CALL DatasetExists(File_ID_DSMC,TRIM('/'//TRIM(spec_pair)//'/ELASTIC'),DatasetFound)
  ELSE 
    CALL DatasetExists(File_ID_DSMC,TRIM(TRIM(spec_pair)//'/ELASTIC'),DatasetFound)
  END IF
  IF(DatasetFound) CALL abort(__STAMP__,'ERROR: Please provide either elastic or effective collision cross-section data '//&
                                             & 'for '//TRIM(spec_pair)//'.')
  LBWRITE(UNIT_StdOut,'(A)') ' | '//TRIM(spec_pair)//': Found EFFECTIVE collision cross section.'
ELSE

  IF (SpeciesDatabase.EQ.'none') THEN
    dsetname = TRIM('/'//TRIM(spec_pair)//'/ELASTIC')
  ELSE
    dsetname = TRIM(TRIM(spec_pair)//'/ELASTIC')
  END IF
  CALL DatasetExists(File_ID_DSMC,TRIM(dsetname),DatasetFound)
  IF(DatasetFound) THEN
    SpecXSec(iCase)%UseCollXSec = .TRUE.
    SpecXSec(iCase)%CollXSec_Effective = .FALSE.
  LBWRITE(UNIT_StdOut,'(A)') ' | '//TRIM(spec_pair)//': Found ELASTIC collision cross section.'
  ELSE
    LBWRITE(UNIT_StdOut,'(A)') ' | '//TRIM(spec_pair)//': No data set found. Using standard collision modelling.'
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
USE MOD_Particle_Vars    ,ONLY: SpeciesDatabase, Species
USE MOD_MCC_Vars         ,ONLY: XSec_Database, SpecXSec
USE MOD_HDF5_Input       ,ONLY: DatasetExists
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
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
CHARACTER(LEN=265)                :: XSecDatabaseName
INTEGER                           :: err
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
spec_pair = TRIM(Species(jSpec)%Name)//'-'//TRIM(Species(iSpec)%Name)

GroupFound = .FALSE.
SpecXSec(iCase)%UseVibXSec = .FALSE.

IF(SpeciesDatabase.EQ.'none') THEN
  XSecDatabaseName = TRIM(XSec_Database)
ELSE 
  XSecDatabaseName = TRIM(SpeciesDatabase)
  spec_pair = TRIM('/Cross-Sections/'//Species(jSpec)%Name)//'-'//TRIM(Species(iSpec)%Name)
END IF

! Initialize FORTRAN interface.
CALL H5OPEN_F(err)

! Check if file exists
IF(.NOT.FILEEXISTS(XSecDatabaseName)) THEN
  CALL abort(__STAMP__,'ERROR: Database '//TRIM(XSec_Database)//' does not exist.')
END IF

! Open the file.
CALL H5FOPEN_F (TRIM(XSecDatabaseName), H5F_ACC_RDONLY_F, file_id_dsmc, err)

! Check if the species pair group exists
CALL H5LEXISTS_F(file_id_dsmc, TRIM(spec_pair), GroupFound, err)
IF(.NOT.GroupFound) THEN
  ! Try to swap the species names
  IF (SpeciesDatabase.EQ.'none') THEN
    spec_pair = TRIM(Species(iSpec)%Name)//'-'//TRIM(Species(jSpec)%Name)
  ELSE 
    spec_pair = TRIM('/Cross-Sections/'//Species(iSpec)%Name)//'-'//TRIM(Species(jSpec)%Name)
  END IF
  CALL H5LEXISTS_F(file_id_dsmc, TRIM(spec_pair), GroupFound, err)
  IF(.NOT.GroupFound) THEN
    LBWRITE(UNIT_StdOut,'(A)') ' | '//TRIM(spec_pair)//': No vibrational excitation cross sections found, using constant read-in values.'
    RETURN
  END IF
END IF

! Check if the vibrational cross-section group exists
groupname = TRIM(spec_pair)//'/VIBRATION/'
CALL H5LEXISTS_F(file_id_dsmc, TRIM(groupname), GroupFound, err)
IF(.NOT.GroupFound) THEN
  LBWRITE(UNIT_StdOut,'(A)') ' | '//TRIM(spec_pair)//': No vibrational excitation cross sections found, using constant read-in values.'
  RETURN
END IF

IF(GroupFound) THEN
  CALL H5GOPEN_F(file_id_dsmc,TRIM(groupname), group_id, err)
  call H5Gget_info_f(group_id, storage, nVib,max_corder, err)
  ! If cross-section data is found, set the corresponding flag
  IF(nVib.GT.0) THEN
    LBWRITE(UNIT_StdOut,'(A,I3,A)') ' | '//TRIM(spec_pair)//': Found ', nVib,' vibrational excitation cross section(s).'
    SpecXSec(iCase)%UseVibXSec = .TRUE.
  ELSE
    LBWRITE(UNIT_StdOut,'(A)') ' | '//TRIM(spec_pair)//': No vibrational excitation cross sections found, using constant read-in values.'
  END IF
ELSE
  LBWRITE(UNIT_StdOut,'(A)') ' | '//TRIM(spec_pair)//': No vibrational excitation cross sections found, using constant read-in values.'
END IF

IF(SpecXSec(iCase)%UseVibXSec) THEN
  ALLOCATE(SpecXSec(iCase)%VibMode(1:nVib))
  DO iVib = 0, nVib-1
    ! Get name and size of name
    CALL H5Lget_name_by_idx_f(group_id, ".", H5_INDEX_NAME_F, H5_ITER_INC_F, iVib, dsetname, err, size)
    READ(dsetname,*) SpecXSec(iCase)%VibMode(iVib+1)%Threshold
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


SUBROUTINE ReadElecXSec(iCase,iSpec,jSpec)
!===================================================================================================================================
!> Read-in of electronic cross-sections from a database. Dataset name is composed of SpeciesName-SpeciesName (e.g. Ar-electron)
!> Trying to swap the species indices if dataset not found.
!===================================================================================================================================
! use module
USE MOD_io_hdf5
USE MOD_Globals
USE MOD_Particle_Vars    ,ONLY: SpeciesDatabase, Species
USE MOD_MCC_Vars         ,ONLY: XSec_Database, SpecXSec
USE MOD_HDF5_Input       ,ONLY: DatasetExists
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
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
CHARACTER(LEN=265)                :: XSecDatabaseName
INTEGER                           :: err
INTEGER(HSIZE_T), DIMENSION(2)    :: dims,sizeMax
INTEGER(HID_T)                    :: file_id_dsmc                       ! File identifier
INTEGER(HID_T)                    :: group_id                           ! Group identifier
INTEGER(HID_T)                    :: dset_id_dsmc                       ! Dataset identifier
INTEGER(HID_T)                    :: filespace                          ! filespace identifier
INTEGER(SIZE_T)                   :: size                               ! Size of name
INTEGER(HSIZE_T)                  :: iElec                              ! Index
LOGICAL                           :: GroupFound
INTEGER                           :: storage, nElec, max_corder
!===================================================================================================================================

GroupFound = .FALSE.
SpecXSec(iCase)%UseElecXSec = .FALSE.
SpecXSec(iCase)%NumElecLevel = 0

IF(SpeciesDatabase.EQ.'none') THEN
  XSecDatabaseName = TRIM(XSec_Database)
  spec_pair = TRIM(Species(jSpec)%Name)//'-'//TRIM(Species(iSpec)%Name)
ELSE 
  XSecDatabaseName = TRIM(SpeciesDatabase)
  spec_pair = TRIM('/Cross-Sections/'//Species(jSpec)%Name)//'-'//TRIM(Species(iSpec)%Name)
END IF

! Initialize FORTRAN interface.
CALL H5OPEN_F(err)

! Check if file exists
IF(.NOT.FILEEXISTS(XSecDatabaseName)) THEN
  CALL abort(__STAMP__,'ERROR: Database '//TRIM(XSecDatabaseName)//' does not exist.')
END IF

! Open the file.
CALL H5FOPEN_F (TRIM(XSecDatabaseName), H5F_ACC_RDONLY_F, file_id_dsmc, err)

! Check if the species pair group exists
CALL H5LEXISTS_F(file_id_dsmc, TRIM(spec_pair), GroupFound, err)
IF(.NOT.GroupFound) THEN
  ! Try to swap the species names
  IF (SpeciesDatabase.EQ.'none') THEN
    spec_pair = TRIM(Species(iSpec)%Name)//'-'//TRIM(Species(jSpec)%Name)
  ELSE 
    spec_pair = TRIM('/Cross-Sections/'//Species(iSpec)%Name)//'-'//TRIM(Species(jSpec)%Name)
  END IF
  CALL H5LEXISTS_F(file_id_dsmc, TRIM(spec_pair), GroupFound, err)
  IF(.NOT.GroupFound) THEN
    LBWRITE(UNIT_StdOut,'(A)') ' | '//TRIM(spec_pair)//': No electronic excitation cross sections found, using constant read-in values.'
    RETURN
  END IF
END IF

! Check if the electronic cross-section group exists
groupname = TRIM(spec_pair)//'/ELECTRONIC/'
CALL H5LEXISTS_F(file_id_dsmc, TRIM(groupname), GroupFound, err)
IF(.NOT.GroupFound) THEN
  LBWRITE(UNIT_StdOut,'(A)') ' | '//TRIM(spec_pair)//': No electronic excitation cross sections found, using constant read-in values.'
  RETURN
END IF

IF(GroupFound) THEN
  CALL H5GOPEN_F(file_id_dsmc,TRIM(groupname), group_id, err)
  call H5Gget_info_f(group_id, storage, nElec,max_corder, err)
  ! If cross-section data is found, set the corresponding flag
  IF(nElec.GT.0) THEN
    LBWRITE(UNIT_StdOut,'(A,I3,A)') TRIM(spec_pair)//': Found ', nElec,' electronic excitation cross section(s).'
    SpecXSec(iCase)%UseElecXSec = .TRUE.
    SpecXSec(iCase)%NumElecLevel = nElec
  ELSE
    LBWRITE(UNIT_StdOut,'(A)') ' | '//TRIM(spec_pair)//': No electronic excitation cross sections found, using constant read-in values.'
  END IF
ELSE
  LBWRITE(UNIT_StdOut,'(A)') ' | '//TRIM(spec_pair)//': No electronic excitation cross sections found, using constant read-in values.'
END IF

IF(SpecXSec(iCase)%UseElecXSec) THEN
  ALLOCATE(SpecXSec(iCase)%ElecLevel(1:nElec))
  DO iElec = 0, nElec-1
    ! Get name and size of name
    CALL H5Lget_name_by_idx_f(group_id, ".", H5_INDEX_NAME_F, H5_ITER_INC_F, iElec, dsetname, err, size)
    READ(dsetname,*) SpecXSec(iCase)%ElecLevel(iElec+1)%Threshold
    SpecXSec(iCase)%ElecLevel(iElec+1)%Counter = 0.
    dsetname = TRIM(groupname)//TRIM(dsetname)
    ! Open the dataset.
    CALL H5DOPEN_F(file_id_dsmc, dsetname, dset_id_dsmc, err)
    ! Get the file space of the dataset.
    CALL H5DGET_SPACE_F(dset_id_dsmc, FileSpace, err)
    ! get size
    CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, dims, SizeMax, err)
    ALLOCATE(SpecXSec(iCase)%ElecLevel(iElec+1)%XSecData(dims(1),dims(2)))
    ! read data
    CALL H5DREAD_F(dset_id_dsmc, H5T_NATIVE_DOUBLE, SpecXSec(iCase)%ElecLevel(iElec+1)%XSecData, dims, err)
  END DO
  ! Close the group
  CALL H5GCLOSE_F(group_id,err)
END IF

! Close the file.
CALL H5FCLOSE_F(file_id_dsmc, err)
! Close FORTRAN interface.
CALL H5CLOSE_F(err)

END SUBROUTINE ReadElecXSec


PPURE REAL FUNCTION InterpolateCrossSection(XSecData,CollisionEnergy)
!===================================================================================================================================
!> Interpolate the collision cross-section [m^2] from the available data at the given collision energy [J]
!> Collision energies below and above the given data will be set at the first and last level of the data set
!> Note: Requires the data to be sorted by ascending energy values
!> Assumption: First species given is the particle species, second species input is the background gas species
!===================================================================================================================================
! MODULES
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,INTENT(IN)               :: XSecData(:,:)                    !< Array for the interpolation [1:2,1:MaxDOF]
REAL,INTENT(IN)               :: CollisionEnergy                  !< Collision energy in [J]
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iDOF, MaxDOF
!===================================================================================================================================

InterpolateCrossSection = 0.
MaxDOF = SIZE(XSecData,2)

IF(CollisionEnergy.GT.XSecData(1,MaxDOF)) THEN
  ! If the collision energy is greater than the maximal value, extrapolate from the last two values
  IF((MaxDOF.LT.2).OR.(XSecData(2,MaxDOF).LE.0.))THEN
    ! If only one value is given or the last cross-section is zero
    InterpolateCrossSection = XSecData(2,MaxDOF)
  ELSE
    ! Extrapolate
    InterpolateCrossSection = XSecData(2,MaxDOF-1) + (CollisionEnergy - XSecData(1,MaxDOF-1)) &
                            / (XSecData(1,MaxDOF) - XSecData(1,MaxDOF-1)) * (XSecData(2,MaxDOF) - XSecData(2,MaxDOF-1))
    ! Check if extrapolation drops under zero
    IF(InterpolateCrossSection.LE.0.) InterpolateCrossSection = 0.
  END IF ! (MaxDOF.LT.2).OR.(XSecData(2,MaxDOF).LE.0.)
  ! Leave routine
  RETURN
ELSE IF(CollisionEnergy.LE.XSecData(1,1)) THEN
  ! If collision energy is below the minimal value, get the cross-section of the first level and leave routine
  InterpolateCrossSection = XSecData(2,1)
  ! Leave routine
  RETURN
END IF

DO iDOF = 1, MaxDOF
  ! Check if the stored energy value is above the collision energy
  IF(XSecData(1,iDOF).GT.CollisionEnergy) THEN
    ! Interpolate the cross-section from the data set using the current and the energy level below
    InterpolateCrossSection = XSecData(2,iDOF-1) + (CollisionEnergy  - XSecData(1,iDOF-1)) &
                              / (XSecData(1,iDOF) - XSecData(1,iDOF-1)) * (XSecData(2,iDOF) - XSecData(2,iDOF-1))
    ! Leave routine and do not finish DO loop
    RETURN
  ELSE IF(CollisionEnergy.EQ.XSecData(1,iDOF)) THEN
    ! In case the collision energy is exactly the stored value
    InterpolateCrossSection = XSecData(2,iDOF)
  END IF
END DO

END FUNCTION InterpolateCrossSection


SUBROUTINE XSec_CalcCollisionProb(iPair,iElem,SpecNum1,SpecNum2,CollCaseNum,MacroParticleFactor,Volume,dtCell)
!===================================================================================================================================
!> Calculate the collision probability if collision cross-section data is used. Can be utilized in combination with the regular
!> DSMC collision calculation probability.
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars             ,ONLY: SpecDSMC, Coll_pData, CollInf, BGGas, RadialWeighting
USE MOD_MCC_Vars              ,ONLY: SpecXSec, XSec_NullCollision
USE MOD_Particle_Vars         ,ONLY: PartSpecies, Species, UseVarTimeStep, usevMPF
USE MOD_part_tools            ,ONLY: GetParticleWeight
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN) :: iPair,iElem
REAL,INTENT(IN)    :: SpecNum1, SpecNum2, CollCaseNum, MacroParticleFactor, Volume, dtCell
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: CollEnergy, SpecNumTarget, SpecNumSource, Weight1, Weight2, ReducedMass, ReducedMassUnweighted
INTEGER            :: targetSpec, iPart_p1, iPart_p2, iSpec_p1, iSpec_p2, iCase, bgSpec
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
  IF (RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep.OR.usevMPF) THEN
    ReducedMass = (Species(iSpec_p1)%MassIC *Weight1  * Species(iSpec_p2)%MassIC * Weight2) &
      / (Species(iSpec_p1)%MassIC * Weight1+ Species(iSpec_p2)%MassIC * Weight2)
    ReducedMassUnweighted = ReducedMass * 2. / (Weight1 + Weight2)
  ELSE
    ReducedMass = CollInf%MassRed(Coll_pData(iPair)%PairType)
    ReducedMassUnweighted = CollInf%MassRed(Coll_pData(iPair)%PairType)
  END IF
  ! Using the relative kinetic energy of the particle pair (real energy value per particle pair, no weighting/scaling factors)
  CollEnergy = 0.5 * ReducedMassUnweighted * Coll_pData(iPair)%CRela2
  ! Calculate the collision probability
  SpecXSec(iCase)%CrossSection = InterpolateCrossSection(SpecXSec(iCase)%CollXSecData,CollEnergy)
  Coll_pData(iPair)%Prob = (1. - EXP(-SQRT(Coll_pData(iPair)%CRela2) * SpecXSec(iCase)%CrossSection * SpecNumTarget * MacroParticleFactor &
                                        / Volume * dtCell))
  ! Correction for conditional probabilities in case of MCC
  IF(BGGas%BackgroundSpecies(targetSpec)) THEN
    ! Correct the collision probability in the case of the second species being a background species as the number of pairs
    ! is either determined based on the null collision probability or on the species fraction
    IF(XSec_NullCollision) THEN
      IF(BGGas%UseDistribution) THEN
        Coll_pData(iPair)%Prob = Coll_pData(iPair)%Prob / SpecXSec(iCase)%ProbNullElem(iElem)
      ELSE
        Coll_pData(iPair)%Prob = Coll_pData(iPair)%Prob / SpecXSec(iCase)%ProbNull
      END IF
    ELSE
      bgSpec = BGGas%MapSpecToBGSpec(targetSpec)
      IF(BGGas%UseDistribution)THEN
        Coll_pData(iPair)%Prob = Coll_pData(iPair)%Prob / BGGas%SpeciesFractionElem(bgSpec,iElem)
      ELSE
        Coll_pData(iPair)%Prob = Coll_pData(iPair)%Prob / BGGas%SpeciesFraction(bgSpec)
      END IF ! BGGas%UseDistribution
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


SUBROUTINE XSec_CalcVibRelaxProb(iPair,iElem,SpecNum1,SpecNum2,MacroParticleFactor,Volume,dtCell)
!===================================================================================================================================
!> Calculate the relaxation probability using cross-section data.
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars             ,ONLY: SpecDSMC, Coll_pData, CollInf, BGGas, RadialWeighting
USE MOD_MCC_Vars              ,ONLY: SpecXSec
USE MOD_Particle_Vars         ,ONLY: PartSpecies, Species, UseVarTimeStep, usevMPF
USE MOD_part_tools            ,ONLY: GetParticleWeight
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)            :: iPair,iElem
REAL,INTENT(IN),OPTIONAL      :: SpecNum1, SpecNum2, MacroParticleFactor, Volume, dtCell
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: CollEnergy, SpecNumTarget, SpecNumSource, Weight1, Weight2, ReducedMass, SumVibCrossSection
REAL                          :: ReducedMassUnweighted
INTEGER                       :: targetSpec, iPart_p1, iPart_p2, iSpec_p1, iSpec_p2, iCase, iVib, nVib, bgSpec
!===================================================================================================================================

iPart_p1 = Coll_pData(iPair)%iPart_p1; iPart_p2 = Coll_pData(iPair)%iPart_p2
iSpec_p1 = PartSpecies(iPart_p1);      iSpec_p2 = PartSpecies(iPart_p2)
iCase = CollInf%Coll_Case(iSpec_p1,iSpec_p2)
Weight1 = GetParticleWeight(iPart_p1)
Weight2 = GetParticleWeight(iPart_p2)
SpecXSec(iCase)%VibProb = 0.
SumVibCrossSection = 0.

IF (RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep.OR.usevMPF) THEN
  ReducedMass = (Species(iSpec_p1)%MassIC *Weight1  * Species(iSpec_p2)%MassIC * Weight2) &
    / (Species(iSpec_p1)%MassIC * Weight1+ Species(iSpec_p2)%MassIC * Weight2)
  ReducedMassUnweighted = ReducedMass * 2./(Weight1 + Weight2)
ELSE
  ReducedMass = CollInf%MassRed(Coll_pData(iPair)%PairType)
  ReducedMassUnweighted = CollInf%MassRed(Coll_pData(iPair)%PairType)
END IF
! Using the relative translational energy of the pair
CollEnergy = 0.5 * ReducedMassUnweighted * Coll_pData(iPair)%CRela2
! Calculate the total vibrational cross-section
nVib = SIZE(SpecXSec(iCase)%VibMode)
DO iVib = 1, nVib
  SumVibCrossSection = SumVibCrossSection + InterpolateCrossSection(SpecXSec(iCase)%VibMode(iVib)%XSecData,CollEnergy)
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
    bgSpec = BGGas%MapSpecToBGSpec(targetSpec)
    IF(BGGas%UseDistribution)THEN
      SpecXSec(iCase)%VibProb = SpecXSec(iCase)%VibProb / BGGas%SpeciesFractionElem(bgSpec,iElem)
    ELSE
      SpecXSec(iCase)%VibProb = SpecXSec(iCase)%VibProb / BGGas%SpeciesFraction(bgSpec)
    END IF ! BGGas%UseDistribution
  ELSE
    SpecXSec(iCase)%VibProb = SpecXSec(iCase)%VibProb * SpecNumSource / CollInf%Coll_CaseNum(iCase)
  END IF
END IF

END SUBROUTINE XSec_CalcVibRelaxProb


SUBROUTINE XSec_ElectronicRelaxation(iPair,iCase,iPart_p1,iPart_p2,DoElec1,DoElec2,ElecLevelRelax)
!===================================================================================================================================
!> Determines whether a relaxation occurs based on the electronic relaxation probability
!> 1. Interpolate the cross-section (MCC) or use the probability (VHS)
!> 2. Determine which electronic level is to be excited
!> 3. Reduce the total collision probability if no electronic excitation occurred
!> 4. Count the number of relaxation process for the relaxation rate (Only with Particles-DSMCReservoirSim = T)
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars             ,ONLY: SpecDSMC, Coll_pData, PartStateIntEn
USE MOD_MCC_Vars              ,ONLY: SpecXSec
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_Particle_Vars         ,ONLY: PartSpecies, PEM, WriteMacroVolumeValues
USE MOD_Particle_Analyze_Vars ,ONLY: CalcRelaxProb
USE MOD_Particle_Vars         ,ONLY: Species, usevMPF, SampleElecExcitation, ExcitationSampleData, ExcitationLevelMapping
USE MOD_DSMC_Vars             ,ONLY: RadialWeighting, SamplingActive
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: iPair, iCase, iPart_p1, iPart_p2
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)           :: DoElec1, DoElec2
INTEGER,INTENT(OUT)           :: ElecLevelRelax
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iSpec_p1, iSpec_p2, iLevel, ElecLevel, ElemID
REAL                          :: ProbSum, ProbElec, iRan
REAL                          :: WeightedParticle
!===================================================================================================================================

iSpec_p1 = PartSpecies(iPart_p1)
iSpec_p2 = PartSpecies(iPart_p2)

ElecLevelRelax = 0

! Excitation only from ground-state
IF(PartStateIntEn(3,iPart_p1).EQ.0.0.AND.PartStateIntEn(3,iPart_p2).EQ.0.0) THEN
  ! 1. Interpolate the cross-section (MCC) or use the probability (VHS)
  IF(SpecXSec(iCase)%UseCollXSec) THEN
    ! Interpolate the electronic cross-section at the current collision energy
    CALL XSec_CalcElecRelaxProb(iPair)
    ProbSum = SpecXSec(iCase)%CrossSection
  ELSE
    ! Probabilities were saved and added to the total collision probability
    ProbSum = Coll_pData(iPair)%Prob
  END IF
  ! Only proceed if any of the electronic excitation probabilities is above zero
  IF(SUM(SpecXSec(iCase)%ElecLevel(:)%Prob).GT.0.) THEN
    ProbElec = 0.
    ! 2. Decide which electronic excitation should occur
    CALL RANDOM_NUMBER(iRan)
    DO iLevel = 1, SpecXSec(iCase)%NumElecLevel
      ProbElec = ProbElec + SpecXSec(iCase)%ElecLevel(iLevel)%Prob
      IF((ProbElec/ProbSum).GT.iRan) THEN
        IF((Species(iSpec_p1)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec_p1)%FullyIonized)) THEN
          DoElec1 = .TRUE.
        ELSE
          DoElec2 = .TRUE.
        END IF
        ElecLevelRelax = iLevel
        EXIT
      END IF
    END DO
    ! 3. Reducing the total collision probability if no electronic excitation occurred
    IF((.NOT.DoElec1).AND.(.NOT.DoElec2)) THEN
      IF(SpecXSec(iCase)%UseCollXSec) THEN
        SpecXSec(iCase)%CrossSection = SpecXSec(iCase)%CrossSection - SUM(SpecXSec(iCase)%ElecLevel(:)%Prob)
      ELSE
        Coll_pData(iPair)%Prob = Coll_pData(iPair)%Prob - SUM(SpecXSec(iCase)%ElecLevel(:)%Prob)
      END IF
    END IF
  END IF  ! SUM(SpecXSec(iCase)%ElecLevel(:)%Prob).GT.0.
END IF    ! Electronic energy = 0, ground-state

! 4. Count the number of relaxation process for the relaxation rate and cell-local sampling
IF(CalcRelaxProb.OR.SamplingActive.OR.WriteMacroVolumeValues) THEN
  IF(ElecLevelRelax.GT.0) THEN
    IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
      ! Weighting factor already included in GetParticleWeight
      WeightedParticle = 1.
    ELSE
      ! Weighting factor should be the same for all species anyway (BGG: first species is the non-BGG particle species)
      WeightedParticle = Species(iSpec_p1)%MacroParticleFactor
    END IF
    IF(DoElec1) THEN
      WeightedParticle = GetParticleWeight(iPart_p1) * WeightedParticle
    ELSE IF(DoElec2) THEN
      WeightedParticle = GetParticleWeight(iPart_p2) * WeightedParticle
    END IF
    SpecXSec(iCase)%ElecLevel(ElecLevelRelax)%Counter = SpecXSec(iCase)%ElecLevel(ElecLevelRelax)%Counter + WeightedParticle
    ! Sampling of the cell-local excitation rate
    IF(SampleElecExcitation) THEN
      ElecLevel = ExcitationLevelMapping(iCase,ElecLevelRelax)
      ElemID = PEM%LocalElemID(iPart_p1)
      ExcitationSampleData(ElecLevel,ElemID) = ExcitationSampleData(ElecLevel,ElemID) + WeightedParticle
    END IF
  END IF
END IF

END SUBROUTINE XSec_ElectronicRelaxation


SUBROUTINE XSec_CalcElecRelaxProb(iPair,SpecNum1,SpecNum2,MacroParticleFactor,Volume,dtCell)
!===================================================================================================================================
!> Calculate the electronic relaxation probability using cross-section data.
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars          ,ONLY: ElementaryCharge
USE MOD_DSMC_Vars             ,ONLY: SpecDSMC, Coll_pData, CollInf, BGGas, RadialWeighting
USE MOD_MCC_Vars              ,ONLY: SpecXSec
USE MOD_Particle_Vars         ,ONLY: PartSpecies, Species, UseVarTimeStep, usevMPF
USE MOD_part_tools            ,ONLY: GetParticleWeight
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: iPair
REAL,INTENT(IN),OPTIONAL      :: SpecNum1, SpecNum2, MacroParticleFactor, Volume, dtCell
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: CollEnergy, SpecNumTarget, SpecNumSource, Weight1, Weight2, ReducedMass
REAL                          :: ReducedMassUnweighted
INTEGER                       :: targetSpec, iPart_p1, iPart_p2, iSpec_p1, iSpec_p2, iCase, iLevel
!===================================================================================================================================

iPart_p1 = Coll_pData(iPair)%iPart_p1; iPart_p2 = Coll_pData(iPair)%iPart_p2
iSpec_p1 = PartSpecies(iPart_p1);      iSpec_p2 = PartSpecies(iPart_p2)
iCase = CollInf%Coll_Case(iSpec_p1,iSpec_p2)
Weight1 = GetParticleWeight(iPart_p1)
Weight2 = GetParticleWeight(iPart_p2)
SpecXSec(iCase)%ElecLevel(:)%Prob = 0.

IF (RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep.OR.usevMPF) THEN
  ReducedMass = (Species(iSpec_p1)%MassIC *Weight1  * Species(iSpec_p2)%MassIC * Weight2) &
    / (Species(iSpec_p1)%MassIC * Weight1+ Species(iSpec_p2)%MassIC * Weight2)
  ReducedMassUnweighted = ReducedMass * 2./(Weight1 + Weight2)
ELSE
  ReducedMass = CollInf%MassRed(Coll_pData(iPair)%PairType)
  ReducedMassUnweighted = CollInf%MassRed(Coll_pData(iPair)%PairType)
END IF
! Using the relative translational energy of the pair
CollEnergy = 0.5 * ReducedMassUnweighted * Coll_pData(iPair)%CRela2
! Calculate the electronic cross-section
DO iLevel = 1, SpecXSec(iCase)%NumElecLevel
  IF(CollEnergy.GT.SpecXSec(iCase)%ElecLevel(iLevel)%Threshold) THEN
    SpecXSec(iCase)%ElecLevel(iLevel)%Prob = InterpolateCrossSection(SpecXSec(iCase)%ElecLevel(iLevel)%XSecData,CollEnergy)
  END IF
END DO

IF(.NOT.SpecXSec(iCase)%UseCollXSec) THEN
  IF(SpecDSMC(iSpec_p1)%UseElecXSec) THEN
    targetSpec = iSpec_p2; SpecNumTarget = SpecNum2; SpecNumSource = SpecNum1
  ELSE
    targetSpec = iSpec_p1; SpecNumTarget = SpecNum1; SpecNumSource = SpecNum2
  END IF
  DO iLevel = 1, SpecXSec(iCase)%NumElecLevel
    IF(CollEnergy.GT.SpecXSec(iCase)%ElecLevel(iLevel)%Threshold) THEN
      ! Calculate the electronic relaxation probability
      SpecXSec(iCase)%ElecLevel(iLevel)%Prob = (1. - EXP(-SQRT(Coll_pData(iPair)%CRela2) * SpecNumTarget * MacroParticleFactor / Volume &
                                                      * SpecXSec(iCase)%ElecLevel(iLevel)%Prob * dtCell))
      IF(BGGas%BackgroundSpecies(targetSpec)) THEN
        ! Correct the collision probability in the case of the second species being a background species as the number of pairs
        ! is determined based on the species fraction
        SpecXSec(iCase)%ElecLevel(iLevel)%Prob = SpecXSec(iCase)%ElecLevel(iLevel)%Prob / BGGas%SpeciesFraction(BGGas%MapSpecToBGSpec(targetSpec))
      ELSE
        SpecXSec(iCase)%ElecLevel(iLevel)%Prob = SpecXSec(iCase)%ElecLevel(iLevel)%Prob * SpecNumSource / CollInf%Coll_CaseNum(iCase)
      END IF
    END IF
  END DO
END IF

END SUBROUTINE XSec_CalcElecRelaxProb


SUBROUTINE ReadReacXSec(iCase,iPath)
!===================================================================================================================================
!> Read-in of reaction cross-sections from a given database. Using the effective cross-section database to check whether the
!> group exists. Trying to swap the species indices if dataset not found.
!===================================================================================================================================
! use module
USE MOD_io_hdf5
USE MOD_Globals
USE MOD_Particle_Vars    ,ONLY: SpeciesDatabase, Species
USE MOD_Globals_Vars     ,ONLY: ElementaryCharge,Joule2eV
USE MOD_DSMC_Vars        ,ONLY: ChemReac
USE MOD_MCC_Vars         ,ONLY: XSec_Database, SpecXSec
USE MOD_HDF5_Input       ,ONLY: DatasetExists
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
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
CHARACTER(LEN=265)                :: XSecDatabaseName
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
REAL, ALLOCATABLE                 :: tempArray(:,:)
CHARACTER(LEN=32)                 :: hilf
REAL                              :: ERatio
!===================================================================================================================================
iReac = ChemReac%CollCaseInfo(iCase)%ReactionIndex(iPath)
EductReac(1:3) = ChemReac%Reactants(iReac,1:3)
ProductReac(1:4) = ChemReac%Products(iReac,1:4)

DatasetFound = .FALSE.; GroupFound = .FALSE.


IF (SpeciesDatabase.EQ.'none') THEN
  EductPair = TRIM(Species(EductReac(1))%Name)//'-'//TRIM(Species(EductReac(2))%Name)
ELSE
  EductPair = TRIM('/Cross-Sections/'//TRIM(Species(EductReac(1))%Name)//'-'//TRIM(Species(EductReac(2))%Name))
END IF
CALL H5LEXISTS_F(file_id_dsmc, TRIM(EductPair), GroupFound, err)

IF(SpeciesDatabase.EQ.'none') THEN
  XSecDatabaseName = TRIM(XSec_Database)
ELSE 
  XSecDatabaseName = TRIM(SpeciesDatabase)
END IF

! Initialize FORTRAN interface.
CALL H5OPEN_F(err)
! Open the file.
CALL H5FOPEN_F (TRIM(XSecDatabaseName), H5F_ACC_RDONLY_F, file_id_dsmc, err)

! Check if the species pair group exists
CALL H5LEXISTS_F(file_id_dsmc, TRIM(EductPair), GroupFound, err)
IF(.NOT.GroupFound) THEN
  ! Try to swap the species names
  IF (SpeciesDatabase.EQ.'none') THEN
    EductPair = TRIM(Species(EductReac(2))%Name)//'-'//TRIM(Species(EductReac(1))%Name)
  ELSE
    EductPair = TRIM('/Cross-Sections/'//TRIM(Species(EductReac(2))%Name)//'-'//TRIM(Species(EductReac(1))%Name))
  END IF
  CALL H5LEXISTS_F(file_id_dsmc, TRIM(EductPair), GroupFound, err)
  IF(.NOT.GroupFound) CALL abort(__STAMP__,'No reaction cross sections found in database for reaction number:',IntInfoOpt=iReac)
END IF

IF (SpeciesDatabase.EQ.'none') THEN
  groupname = TRIM('/'//TRIM(EductPair)//'/REACTION/')
ELSE
  groupname = TRIM(TRIM(EductPair)//'/REACTION/')
END IF
CALL H5LEXISTS_F(file_id_dsmc, groupname, GroupFound, err)
IF(.NOT.GroupFound) CALL abort(__STAMP__,'No reaction cross sections found in database for reaction number:',IntInfoOpt=iReac)

CALL H5GOPEN_F(file_id_dsmc,TRIM(groupname), group_id, err)
CALL H5Gget_info_f(group_id, storage, nSets,max_corder, err)
IF(nSets.EQ.0) CALL abort(__STAMP__,'No reaction cross sections found in database for reaction number:',IntInfoOpt=iReac)

SELECT CASE(COUNT(ProductReac.GT.0))
CASE(2)
  ProductPair = TRIM(Species(ProductReac(1))%Name)//'-'//TRIM(Species(ProductReac(2))%Name)
CASE(3)
  ProductPair = TRIM(Species(ProductReac(1))%Name)//'-'//TRIM(Species(ProductReac(2))%Name)//&
                &'-'//TRIM(Species(ProductReac(3))%Name)
CASE(4)
  ProductPair = TRIM(Species(ProductReac(1))%Name)//'-'//TRIM(Species(ProductReac(2))%Name)//&
                &'-'//TRIM(Species(ProductReac(3))%Name)//'-'//TRIM(Species(ProductReac(4))%Name)
CASE DEFAULT
  CALL abort(__STAMP__,'Number of products is not supported yet! Reaction number:', IntInfoOpt=iReac)
END SELECT

LBWRITE(UNIT_StdOut,'(A,I0)') ' | Reading cross section for '//TRIM(EductPair)//' -> '//TRIM(ProductPair)//' as reaction # ', iReac

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

IF(.NOT.ReactionFound) THEN
  SWRITE(*,*) ' Trying to read the following product pair failed: ', TRIM(ProductPair)
  CALL abort(__STAMP__,'No reaction cross-section data found for reaction number:', IntInfoOpt=iReac)
END IF

! Store the energy value in J (read-in was in eV)
SpecXSec(iCase)%ReactionPath(iPath)%XSecData(1,:) = SpecXSec(iCase)%ReactionPath(iPath)%XSecData(1,:) * ElementaryCharge

!
IF(SpecXSec(iCase)%ReactionPath(iPath)%XSecData(2,1).GT.0.0) THEN
  ! Store old data in temporary array
  ALLOCATE(tempArray(dims(1),dims(2)))
  tempArray = SpecXSec(iCase)%ReactionPath(iPath)%XSecData
  ! Allocate new array with one additional row
  DEALLOCATE(SpecXSec(iCase)%ReactionPath(iPath)%XSecData)
  ALLOCATE(SpecXSec(iCase)%ReactionPath(iPath)%XSecData(dims(1),dims(2)+1))
  ! Create new cross-section at the reaction threshold
  SpecXSec(iCase)%ReactionPath(iPath)%XSecData(1,1) = -ChemReac%EForm(iReac)
  SpecXSec(iCase)%ReactionPath(iPath)%XSecData(2,1) = 0.0
  ! Store the read-in dataset
  SpecXSec(iCase)%ReactionPath(iPath)%XSecData(1:dims(1),2:dims(2)+1) = tempArray(1:dims(1),1:dims(2))
  DEALLOCATE(tempArray)
  ! Sanity check: Is the heat of formation larger than the first energy level of the cross-section data
  ERatio=SpecXSec(iCase)%ReactionPath(iPath)%XSecData(1,1)/SpecXSec(iCase)%ReactionPath(iPath)%XSecData(1,2)
  IF(ERatio.GT.1.0) THEN
    LBWRITE(*,*) '      (Negative) Heat of reaction [J]: ', -ChemReac%EForm(iReac),", [eV]: ",-ChemReac%EForm(iReac)*Joule2eV
    LBWRITE(*,*) ' First energy level from database [J]: ', SpecXSec(iCase)%ReactionPath(iPath)%XSecData(1,2),", [eV]: ",&
    SpecXSec(iCase)%ReactionPath(iPath)%XSecData(1,2)*Joule2eV
    WRITE(UNIT=hilf,FMT='(F8.2)') ERatio
    LBWRITE (*,'(A,I0)') 'Warning: Heat of reaction is factor '//TRIM(ADJUSTL(hilf))//&
        ' greater than the first read-in energy level for reaction number: ',iReac
    IF(ERatio.GT.10.0) CALL abort(__STAMP__,&
        ' Heat of reaction is much greater than the first read-in energy level for reaction number:', iReac)
  END IF
END IF

! Close the group
CALL H5GCLOSE_F(group_id,err)
! Close the file.
CALL H5FCLOSE_F(file_id_dsmc, err)
! Close FORTRAN interface.
CALL H5CLOSE_F(err)

END SUBROUTINE ReadReacXSec


SUBROUTINE ReadReacPhotonXSec(iPhotoReac)
!===================================================================================================================================
!> Read-in of photoionization reaction cross-sections from a given database. Check whether the
!> group exists. Trying to swap the species indices if dataset not found.
!===================================================================================================================================
! use module
USE MOD_io_hdf5
USE MOD_Globals
USE MOD_Globals_Vars     ,ONLY: ElementaryCharge
USE MOD_Particle_Vars    ,ONLY: SpeciesDatabase, Species
USE MOD_DSMC_Vars        ,ONLY: ChemReac
USE MOD_HDF5_Input       ,ONLY: DatasetExists
USE MOD_MCC_Vars         ,ONLY: XSec_Database,SpecPhotonXSec,PhotoReacToReac
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: iPhotoReac
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=128)                :: dsetname, groupname, EductPair, dsetname2, ProductPair
CHARACTER(LEN=265)                :: XSecDatabaseName
INTEGER                           :: err
INTEGER(HSIZE_T), DIMENSION(2)    :: dims,sizeMax
INTEGER(HID_T)                    :: file_id_dsmc                       ! File identifier
INTEGER(HID_T)                    :: group_id                           ! Group identifier
INTEGER(HID_T)                    :: dset_id_dsmc                       ! Dataset identifier
INTEGER(HID_T)                    :: filespace                          ! filespace identifier
INTEGER(SIZE_T)                   :: size                               ! Size of name
INTEGER(HSIZE_T)                  :: iSet                               ! Index
INTEGER                           :: storage, nSets, max_corder
LOGICAL                           :: GroupFound, ReactionFound
INTEGER                           :: iReac, ProductReac(1:4)
!===================================================================================================================================
iReac = PhotoReacToReac(iPhotoReac)

! Sanity check
IF(ChemReac%Reactants(iReac,1).EQ.0) CALL abort(__STAMP__,'phIonXSec reaction: The first educt species ID cannot be zero.')
EductPair = TRIM(Species(ChemReac%Reactants(iReac,1))%Name)//'-photon'
ProductReac(1:4) = ChemReac%Products(iReac,1:4)

GroupFound = .FALSE.

IF(SpeciesDatabase.EQ.'none') THEN
  XSecDatabaseName = TRIM(XSec_Database)
ELSE 
  XSecDatabaseName = TRIM(SpeciesDatabase)
  EductPair = TRIM('Cross-Sections/'//Species(ChemReac%Reactants(iReac,1))%Name)//'-photon'
END IF

! Initialize FORTRAN interface.
CALL H5OPEN_F(err)
! Open the file.
CALL H5FOPEN_F (TRIM(XSecDatabaseName), H5F_ACC_RDONLY_F, file_id_dsmc, err)

! Check if the species pair group exists
CALL H5LEXISTS_F(file_id_dsmc, TRIM(EductPair), GroupFound, err)
IF(.NOT.GroupFound) THEN
  ! Try to swap the species names
  IF (SpeciesDatabase.EQ.'none') THEN
    EductPair = 'photon-'//TRIM(Species(ChemReac%Reactants(iReac,1))%Name)
  ELSE 
  EductPair = '/Cross-Sections/photon-'//TRIM(Species(ChemReac%Reactants(iReac,1))%Name)
  END IF
  CALL H5LEXISTS_F(file_id_dsmc, TRIM(EductPair), GroupFound, err)
  IF(.NOT.GroupFound) CALL abort(__STAMP__,'No reaction cross sections found in database for reaction number:',IntInfoOpt=iReac)
END IF

IF (SpeciesDatabase.EQ.'none') THEN
  groupname = TRIM('/'//TRIM(EductPair)//'/REACTION/')
ELSE 
  groupname = TRIM(TRIM(EductPair)//'/REACTION/')
END IF
CALL H5LEXISTS_F(file_id_dsmc, groupname, GroupFound, err)
IF(.NOT.GroupFound) CALL abort(__STAMP__,'No reaction cross sections found in database for reaction number:',IntInfoOpt=iReac)

CALL H5GOPEN_F(file_id_dsmc,TRIM(groupname), group_id, err)
CALL H5Gget_info_f(group_id, storage, nSets,max_corder, err)
IF(nSets.EQ.0) CALL abort(__STAMP__,'No reaction cross sections found in database for reaction number:',IntInfoOpt=iReac)

LBWRITE(UNIT_StdOut,'(A,I0)') 'Read cross section for '//TRIM(EductPair)//' from '//TRIM(XSecDatabaseName)//' for reaction # ', iReac

SELECT CASE(COUNT(ProductReac.GT.0))
CASE(2)
  ProductPair = TRIM(Species(ProductReac(1))%Name)//'-'//TRIM(Species(ProductReac(2))%Name)
CASE(3)
  ProductPair = TRIM(Species(ProductReac(1))%Name)//'-'//TRIM(Species(ProductReac(2))%Name)//&
                &'-'//TRIM(Species(ProductReac(3))%Name)
CASE(4)
  ProductPair = TRIM(Species(ProductReac(1))%Name)//'-'//TRIM(Species(ProductReac(2))%Name)//&
                &'-'//TRIM(Species(ProductReac(3))%Name)//'-'//TRIM(Species(ProductReac(4))%Name)
CASE DEFAULT
  CALL abort(__STAMP__,'Number of products is not supported yet! Reaction number:', IntInfoOpt=iReac)
END SELECT

! IF (SpeciesDatabase.NE.'none') THEN
!   groupname = TRIM('/Cross-Sections/'//TRIM(ProductPair))
! END IF

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
    ALLOCATE(SpecPhotonXSec(iPhotoReac)%XSecData(dims(1),dims(2)))
    ! read data
    CALL H5DREAD_F(dset_id_dsmc, H5T_NATIVE_DOUBLE, SpecPhotonXSec(iPhotoReac)%XSecData, dims, err)
    ReactionFound = .TRUE.
    ! stop looking for other reaction paths
    EXIT
  END IF
END DO

IF(.NOT.ReactionFound) CALL abort(__STAMP__,'No reaction cross-section data found for reaction number:', IntInfoOpt=iReac)

! Store the energy value in J (read-in was in eV)
!SpecPhotonXSec(iPhotoReac)%XSecData= SpecPhotonXSec(iPhotoReac)%XSecData * ElementaryCharge

! Close the group
CALL H5GCLOSE_F(group_id,err)
! Close the file.
CALL H5FCLOSE_F(file_id_dsmc, err)
! Close FORTRAN interface.
CALL H5CLOSE_F(err)

END SUBROUTINE ReadReacPhotonXSec


SUBROUTINE ReadReacPhotonSpectrum(iPhotoReac)
!===================================================================================================================================
!> TODO:
!===================================================================================================================================
! use module
USE MOD_io_hdf5
USE MOD_Globals
USE MOD_Particle_Vars    ,ONLY: SpeciesDatabase, Species
USE MOD_Globals_Vars     ,ONLY: ElementaryCharge
USE MOD_DSMC_Vars        ,ONLY: ChemReac
USE MOD_MCC_Vars         ,ONLY: XSec_Database, PhotoReacToReac, PhotonSpectrum
USE MOD_HDF5_Input       ,ONLY: DatasetExists
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: iPhotoReac
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=128)                :: dsetname, groupname, EductPair, dsetname2, ProductPair
CHARACTER(LEN=265)                :: XSecDatabaseName
INTEGER                           :: err
INTEGER(HSIZE_T), DIMENSION(2)    :: dims,sizeMax
INTEGER(HID_T)                    :: file_id_dsmc                       ! File identifier
INTEGER(HID_T)                    :: group_id                           ! Group identifier
INTEGER(HID_T)                    :: dset_id_dsmc                       ! Dataset identifier
INTEGER(HID_T)                    :: filespace                          ! filespace identifier
INTEGER(SIZE_T)                   :: size                               ! Size of name
INTEGER(HSIZE_T)                  :: iSet                               ! Index
INTEGER                           :: storage, nSets, max_corder
LOGICAL                           :: GroupFound, ReactionFound
INTEGER                           :: iReac
!===================================================================================================================================
iReac = PhotoReacToReac(iPhotoReac)

! Sanity check
IF(ChemReac%Reactants(iReac,1).EQ.0) CALL abort(__STAMP__,'phIonXSec reaction: The first educt species ID cannot be zero.')
EductPair   = TRIM(Species(ChemReac%Reactants(iReac,1))%Name)//'-photon'
ProductPair = EductPair

GroupFound = .FALSE.

IF(SpeciesDatabase.EQ.'none') THEN
  XSecDatabaseName = TRIM(XSec_Database)
ELSE 
  XSecDatabaseName = TRIM(SpeciesDatabase)
  EductPair   = TRIM('/Cross-Sections/'//Species(ChemReac%Reactants(iReac,1))%Name)//'-photon'
ProductPair = EductPair
END IF

! Initialize FORTRAN interface.
CALL H5OPEN_F(err)
! Open the file.
CALL H5FOPEN_F (TRIM(XSecDatabaseName), H5F_ACC_RDONLY_F, file_id_dsmc, err)

! Check if the species pair group exists
CALL H5LEXISTS_F(file_id_dsmc, TRIM(EductPair), GroupFound, err)
IF(.NOT.GroupFound) THEN
  ! Try to swap the species names
  IF (SpeciesDatabase.NE.'none') THEN
    EductPair = 'photon-'//TRIM(Species(ChemReac%Reactants(iReac,1))%Name)
  ELSE
    EductPair = '/Cross-Sections/photon-'//TRIM(Species(ChemReac%Reactants(iReac,1))%Name)
  END IF
  CALL H5LEXISTS_F(file_id_dsmc, TRIM(EductPair), GroupFound, err)
  IF(.NOT.GroupFound) CALL abort(__STAMP__,'No reaction cross sections found in database for reaction number:',IntInfoOpt=iReac)
END IF

IF (SpeciesDatabase.NE.'none') THEN
    groupname = TRIM('/'//TRIM(EductPair)//'/SPECTRUM/')
  ELSE
    groupname = TRIM(TRIM(EductPair)//'/SPECTRUM/')
  END IF
CALL H5LEXISTS_F(file_id_dsmc, groupname, GroupFound, err)
IF(.NOT.GroupFound) CALL abort(__STAMP__,'No reaction cross sections found in database for reaction number:',IntInfoOpt=iReac)

CALL H5GOPEN_F(file_id_dsmc,TRIM(groupname), group_id, err)
CALL H5Gget_info_f(group_id, storage, nSets,max_corder, err)
IF(nSets.EQ.0) CALL abort(__STAMP__,'No reaction cross sections found in database for reaction number:',IntInfoOpt=iReac)

LBWRITE(UNIT_StdOut,'(A)') 'Read photon energy spectrum for '//TRIM(EductPair)//' from '//TRIM(XSecDatabaseName)

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
    ALLOCATE(PhotonSpectrum(dims(1),dims(2)))
    ! read data
    CALL H5DREAD_F(dset_id_dsmc, H5T_NATIVE_DOUBLE, PhotonSpectrum, dims, err)
    ReactionFound = .TRUE.
    ! stop looking for other reaction paths
    EXIT
  END IF
END DO

IF(.NOT.ReactionFound) CALL abort(__STAMP__,'No reaction cross-section data found for reaction number:', IntInfoOpt=iReac)

! Store the energy value in J (read-in was in eV)
!SpecPhotonXSec(iPhotoReac)%XSecData= SpecPhotonXSec(iPhotoReac)%XSecData * ElementaryCharge

! Close the group
CALL H5GCLOSE_F(group_id,err)
! Close the file.
CALL H5FCLOSE_F(file_id_dsmc, err)
! Close FORTRAN interface.
CALL H5CLOSE_F(err)

END SUBROUTINE ReadReacPhotonSpectrum


SUBROUTINE XSec_CalcReactionProb(iPair,iCase,iElem,SpecNum1,SpecNum2,MacroParticleFactor,Volume)
!===================================================================================================================================
!> Calculate the reaction probability if cross-section data is used. Utilized in combination with the regular DSMC collision
!> calculation probability (non-MCC routines).
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars
USE MOD_Globals               ,ONLY: abort
USE MOD_DSMC_Vars             ,ONLY: SpecDSMC, Coll_pData, CollInf, BGGas, ChemReac, RadialWeighting, DSMC, PartStateIntEn
USE MOD_MCC_Vars              ,ONLY: SpecXSec
USE MOD_Particle_Vars         ,ONLY: PartSpecies, Species, UseVarTimeStep, usevMPF, PartTimeStep
USE MOD_TimeDisc_Vars         ,ONLY: dt
USE MOD_part_tools            ,ONLY: GetParticleWeight
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)            :: iPair, iCase, iElem
REAL,INTENT(IN),OPTIONAL      :: SpecNum1, SpecNum2, MacroParticleFactor, Volume
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPath, ReacTest, EductReac(1:3), ProductReac(1:4), ReactInx(1:2), nPair, iProd
INTEGER                       :: NumWeightProd, targetSpec, bgSpec
REAL                          :: EZeroPoint_Prod, dtCell, Weight(1:4), ReducedMass, ReducedMassUnweighted, CollEnergy, GammaFac
REAL                          :: EZeroPoint_Educt, SpecNumTarget, SpecNumSource, CrossSection, EcRelativistic
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
    IF((Species(EductReac(1))%InterID.EQ.2).OR.(Species(EductReac(1))%InterID.EQ.20)) THEN
      EZeroPoint_Educt = EZeroPoint_Educt + SpecDSMC(EductReac(1))%EZeroPoint * Weight(1)
    END IF
    !---------------------------------------------------------------------------------------------------------------------------------
    ! Testing if the second particle is an atom or molecule, if molecule: is it polyatomic?
    IF((Species(EductReac(2))%InterID.EQ.2).OR.(Species(EductReac(2))%InterID.EQ.20)) THEN
      EZeroPoint_Educt = EZeroPoint_Educt + SpecDSMC(EductReac(2))%EZeroPoint * Weight(2)
    END IF

    IF (RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep.OR.usevMPF) THEN
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
      IF((Species(ProductReac(iProd))%InterID.EQ.2).OR.(Species(ProductReac(iProd))%InterID.EQ.20)) THEN
        EZeroPoint_Prod = EZeroPoint_Prod + SpecDSMC(ProductReac(iProd))%EZeroPoint * Weight(iProd)
      END IF
    END DO

    IF (UseVarTimeStep) THEN
      dtCell = dt * (PartTimeStep(ReactInx(1)) + PartTimeStep(ReactInx(2)))*0.5
    ELSE
      dtCell = dt
    END IF

    ! No relativistic treatment here until the chemistry is done fully considering relativistic effects
    Coll_pData(iPair)%Ec = 0.5 * ReducedMass * Coll_pData(iPair)%cRela2

    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + (PartStateIntEn(1,ReactInx(1)) + PartStateIntEn(2,ReactInx(1))) * Weight(1) &
                                                + (PartStateIntEn(1,ReactInx(2)) + PartStateIntEn(2,ReactInx(2))) * Weight(2)
    IF (DSMC%ElectronicModel.GT.0) THEN
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(3,ReactInx(1))*Weight(1) &
                                                  + PartStateIntEn(3,ReactInx(2))*Weight(2)
    END IF
    ! Check first if sufficient energy is available for the products after the reaction
    ASSOCIATE( ReactionProb => ChemReac%CollCaseInfo(iCase)%ReactionProb(iPath) )
      IF(((Coll_pData(iPair)%Ec-EZeroPoint_Prod).GE.(-ChemReac%EForm(ReacTest)*SUM(Weight)/REAL(NumWeightProd)))) THEN

        ! Check relativistic limit for the interpolation of the XSec data but not for in inquiry above
        IF(Coll_pData(iPair)%cRela2 .LT. RelativisticLimit) THEN
          CollEnergy = (Coll_pData(iPair)%Ec-EZeroPoint_Educt) * 2./(Weight(1)+Weight(2))
        ELSE
          ! Relativistic treatment under the assumption that the velocity of the background species is zero or negligible
          GammaFac = Coll_pData(iPair)%cRela2*c2_inv
          ! Sanity check
          IF(GammaFac.GE.1.0)THEN
            CALL abort(__STAMP__,'X = Coll_pData(iPair)%cRela2*c2_inv >= 1.0, which results in SQRT(1-X) failing.')
          ELSE
            GammaFac = 1./SQRT(1.-GammaFac)
          END IF
          EcRelativistic = (GammaFac-1.) * ReducedMass * c2

          EcRelativistic = EcRelativistic + (PartStateIntEn(1,ReactInx(1)) + PartStateIntEn(2,ReactInx(1))) * Weight(1) &
                                                      + (PartStateIntEn(1,ReactInx(2)) + PartStateIntEn(2,ReactInx(2))) * Weight(2)
          IF (DSMC%ElectronicModel.GT.0) EcRelativistic = EcRelativistic + PartStateIntEn(3,ReactInx(1))*Weight(1) &
                                                                         + PartStateIntEn(3,ReactInx(2))*Weight(2)

          CollEnergy = (EcRelativistic-EZeroPoint_Educt) * 2./(Weight(1)+Weight(2))
        END IF

        CrossSection = InterpolateCrossSection(SpecXSec(iCase)%ReactionPath(iPath)%XSecData,CollEnergy)
        IF(SpecXSec(iCase)%UseCollXSec) THEN
          ! Interpolate the reaction cross-section
          ReactionProb = CrossSection
        ELSE
          ! Select the background species as the target cloud
          IF(BGGas%BackgroundSpecies(PartSpecies(Coll_pData(iPair)%iPart_p1))) THEN
            targetSpec = PartSpecies(Coll_pData(iPair)%iPart_p1); SpecNumTarget = SpecNum1; SpecNumSource = SpecNum2
          ELSE
            targetSpec = PartSpecies(Coll_pData(iPair)%iPart_p2); SpecNumTarget = SpecNum2; SpecNumSource = SpecNum1
          END IF
          ! Calculate the reaction probability
          ReactionProb = (1. - EXP(-SQRT(Coll_pData(iPair)%CRela2) * dtCell * SpecNumTarget &
                                                                        * MacroParticleFactor / Volume * CrossSection))
          IF(BGGas%BackgroundSpecies(targetSpec)) THEN
            ! Correct the reaction probability in the case of the second species being a background species as the number of pairs
            ! is based on the species fraction
            bgSpec = BGGas%MapSpecToBGSpec(targetSpec)
            IF(BGGas%UseDistribution)THEN
              ReactionProb = ReactionProb / BGGas%SpeciesFractionElem(bgSpec,iElem)
            ELSE
              ReactionProb = ReactionProb / BGGas%SpeciesFraction(bgSpec)
            END IF ! BGGas%UseDistribution
          ELSE
            ReactionProb = ReactionProb * SpecNumSource / CollInf%Coll_CaseNum(iCase)
          END IF
        END IF
      ELSE
        ReactionProb = 0.
      END IF
    END ASSOCIATE
    ! Calculation of reaction rate coefficient
    IF (DSMC%ReservoirSimu.AND..NOT.DSMC%ReservoirRateStatistic) THEN
      ChemReac%NumReac(ReacTest) = ChemReac%NumReac(ReacTest) + ChemReac%CollCaseInfo(iCase)%ReactionProb(iPath)
      ChemReac%ReacCount(ReacTest) = ChemReac%ReacCount(ReacTest) + 1
    END IF
  END IF
END DO

END SUBROUTINE XSec_CalcReactionProb

END MODULE MOD_MCC_XSec
