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

INTERFACE InterpolateCrossSection
  MODULE PROCEDURE InterpolateCrossSection
END INTERFACE

PUBLIC :: InterpolateCrossSection
!===================================================================================================================================

CONTAINS

SUBROUTINE MCC_Init()
!===================================================================================================================================
!> Initialization of the MCC algorithm: Read-in of the collision cross-section database and calculation of the maximal collision
!> frequency for the null collision method.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Globals_Vars          ,ONLY: ElementaryCharge
USE MOD_PARTICLE_Vars         ,ONLY: nSpecies
USE MOD_DSMC_Vars             ,ONLY: BGGas, SpecDSMC, XSec_Database, SpecXSec, XSec_NullCollision
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER       :: iSpec, jSpec
REAL          :: TotalProb
!===================================================================================================================================

XSec_Database = TRIM(GETSTR('Particles-CollXSec-Database'))
XSec_NullCollision = GETLOGICAL('Particles-CollXSec-NullCollision')

IF(TRIM(XSec_Database).EQ.'none') THEN
  CALL abort(&
  __STAMP__&
  ,'ERROR: No database for the collision cross-section given!')
END IF

IF (BGGas%NumberOfSpecies.EQ.0) THEN
  CALL abort(&
  __STAMP__&
  ,'ERROR: Usage of read-in collision cross-sections only possible with a background gas!')
END IF

ALLOCATE(SpecXSec(nSpecies,nSpecies))

DO iSpec = 1,nSpecies
  IF(BGGas%BackgroundSpecies(iSpec)) CYCLE
  IF(.NOT.SpecDSMC(iSpec)%UseCollXSec) CYCLE
  TotalProb = 0.
  DO jSpec = 1, nSpecies
    IF(.NOT.BGGas%BackgroundSpecies(jSpec)) CYCLE
    ! Read-in cross-section data for collisions of particles from the background gas and the current species
    CALL ReadCollXSec(iSpec, jSpec)
    ! Store the energy value in J (read-in was in eV)
    SpecXSec(iSpec,jSpec)%CollXSecData(1,:) = SpecXSec(iSpec,jSpec)%CollXSecData(1,:) * ElementaryCharge
    IF(XSec_NullCollision) THEN
      ! Determine the maximum collision frequency for the null collision method
      CALL DetermineNullCollProb(iSpec,jSpec)
      TotalProb = TotalProb + SpecXSec(iSpec,jSpec)%ProbNull
      IF(TotalProb.GT.1.0) THEN
        CALL abort(&
        __STAMP__&
        ,'ERROR: Total null collision probability is above unity. Please reduce the time step! Probability is: '&
        ,RealInfoOpt=TotalProb)
      END IF
    END IF
  END DO
END DO

END SUBROUTINE MCC_Init


SUBROUTINE ReadCollXSec(iSpec,jSpec)
!===================================================================================================================================
!> Read-in of collision cross-sections from a given database. Dataset name is composed of SpeciesName-SpeciesName (e.g. Ar-electron)
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
INTEGER,INTENT(IN)                :: iSpec, jSpec
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=64)                 :: dsetname, dsetname2
INTEGER                           :: err
INTEGER(HSIZE_T), DIMENSION(2)    :: dims,sizeMax
INTEGER(HID_T)                    :: file_id_dsmc                       ! File identifier
INTEGER(HID_T)                    :: dset_id_dsmc                       ! Dataset identifier
INTEGER(HID_T)                    :: filespace                          ! filespace identifier
LOGICAL                           :: DataSetFound
!===================================================================================================================================
dsetname = TRIM(SpecDSMC(jSpec)%Name)//'-'//TRIM(SpecDSMC(iSpec)%Name)
SWRITE(UNIT_StdOut,'(A)') 'Read collision cross section for '//TRIM(dsetname)//' from '//TRIM(XSec_Database)

! Initialize FORTRAN interface.
CALL H5OPEN_F(err)
! Open the file.
CALL H5FOPEN_F (TRIM(XSec_Database), H5F_ACC_RDONLY_F, file_id_dsmc, err)
CALL DatasetExists(File_ID_DSMC,TRIM(dsetname),DataSetFound)

! Check if the dataset exist
IF(.NOT.DataSetFound) THEN
  ! Try to swap the species names
  dsetname = TRIM(SpecDSMC(iSpec)%Name)//'-'//TRIM(SpecDSMC(jSpec)%Name)
  CALL DatasetExists(File_ID_DSMC,TRIM(dsetname),DataSetFound)
  IF(.NOT.DataSetFound) THEN
    dsetname2 = TRIM(SpecDSMC(jSpec)%Name)//'-'//TRIM(SpecDSMC(iSpec)%Name)
    CALL abort(&
    __STAMP__&
    ,'Dataset not found: ['//TRIM(dsetname2)//'] or ['//TRIM(dsetname)//'] not found in ['//TRIM(XSec_Database)//']')
  END IF
END IF

! Open the dataset.
CALL H5DOPEN_F(file_id_dsmc, dsetname, dset_id_dsmc, err)
! Get the file space of the dataset.
CALL H5DGET_SPACE_F(dset_id_dsmc, FileSpace, err)
! get size
CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, dims, SizeMax, err)

ALLOCATE(SpecXSec(iSpec,jSpec)%CollXSecData(dims(1),dims(2)))
! read data
CALL H5DREAD_F(dset_id_dsmc, H5T_NATIVE_DOUBLE, SpecXSec(iSpec,jSpec)%CollXSecData, dims, err)

! Close the file.
CALL H5FCLOSE_F(file_id_dsmc, err)
! Close FORTRAN interface.
CALL H5CLOSE_F(err)

END SUBROUTINE ReadCollXSec


SUBROUTINE DetermineNullCollProb(iSpec,jSpec)
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
INTEGER,INTENT(IN)            :: iSpec                            !< Species index colliding with the background gas
INTEGER,INTENT(IN)            :: jSpec                            !< Species index of the background gas
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: MaxDOF
REAL,ALLOCATABLE              :: Velocity(:)
!===================================================================================================================================

MaxDOF = SIZE(SpecXSec(iSpec,jSpec)%CollXSecData,2)
ALLOCATE(Velocity(MaxDOF))

! Determine the mean relative velocity at the given energy level
Velocity(1:MaxDOF) = SQRT(2.) * SQRT(8.*SpecXSec(iSpec,jSpec)%CollXSecData(1,1:MaxDOF)/(Pi*Species(iSpec)%MassIC))

! Calculate the maximal collision frequency
SpecXSec(iSpec,jSpec)%MaxCollFreq = MAXVAL(Velocity(1:MaxDOF) * SpecXSec(iSpec,jSpec)%CollXSecData(2,1:MaxDOF) &
                                    * BGGas%NumberDensity(BGGas%MapSpecToBGSpec(jSpec)))

! Determine the collision probability
SpecXSec(iSpec,jSpec)%ProbNull = 1. - EXP(-SpecXSec(iSpec,jSpec)%MaxCollFreq*ManualTimeStep)

DEALLOCATE(Velocity)

END SUBROUTINE DetermineNullCollProb


PURE REAL FUNCTION InterpolateCrossSection(iSpec,jSpec,CollisionEnergy)
!===================================================================================================================================
!> Interpolate the collision cross-section [m^2] from the available data at the given collision energy [J]
!> Collision energies below and above the given data will be set at the first and last level of the data set
!> Note: Requires the data to be sorted by ascending energy values
!> Assumption: First species given is the particle species, second species input is the backgroung gas species
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars             ,ONLY: SpecXSec
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)            :: iSpec                            !< Species index colliding with the background gas
INTEGER,INTENT(IN)            :: jSpec                            !< Species index of the background gas
REAL,INTENT(IN)               :: CollisionEnergy                  !< Collision energy in [J]
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iDOF, MaxDOF
!===================================================================================================================================

InterpolateCrossSection = 0.
MaxDOF = SIZE(SpecXSec(iSpec,jSpec)%CollXSecData,2)

IF(CollisionEnergy.GT.SpecXSec(iSpec,jSpec)%CollXSecData(1,MaxDOF)) THEN 
  ! If the collision energy is greater than the maximal value, get the cross-section of the last level and leave routine
  InterpolateCrossSection = SpecXSec(iSpec,jSpec)%CollXSecData(2,MaxDOF)
  ! Leave routine
  RETURN
ELSE IF(CollisionEnergy.LE.SpecXSec(iSpec,jSpec)%CollXSecData(1,1)) THEN
  ! If collision energy is below the minimal value, get the cross-section of the first level and leave routine
  InterpolateCrossSection = SpecXSec(iSpec,jSpec)%CollXSecData(2,1)
  ! Leave routine
  RETURN
END IF

DO iDOF = 1, MaxDOF
  ! Check if the stored energy value is above the collision energy
  IF(SpecXSec(iSpec,jSpec)%CollXSecData(1,iDOF).GT.CollisionEnergy) THEN
    ! Interpolate the cross-section from the data set using the current and the energy level below
    InterpolateCrossSection = SpecXSec(iSpec,jSpec)%CollXSecData(2,iDOF-1) &
                                + (CollisionEnergy - SpecXSec(iSpec,jSpec)%CollXSecData(1,iDOF-1)) &
                                / (SpecXSec(iSpec,jSpec)%CollXSecData(1,iDOF) - SpecXSec(iSpec,jSpec)%CollXSecData(1,iDOF-1)) &
                                * (SpecXSec(iSpec,jSpec)%CollXSecData(2,iDOF) - SpecXSec(iSpec,jSpec)%CollXSecData(2,iDOF-1))
    ! Leave routine and do not finish DO loop
    RETURN
  END IF
END DO

END FUNCTION InterpolateCrossSection

END MODULE MOD_DSMC_SpecXSec