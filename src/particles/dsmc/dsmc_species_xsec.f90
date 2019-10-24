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

INTERFACE XSec_Argon_DravinLotz
  MODULE PROCEDURE XSec_Argon_DravinLotz
END INTERFACE

PUBLIC :: MCC_Init, InterpolateCrossSection
PUBLIC :: XSec_Argon_DravinLotz
!===================================================================================================================================

CONTAINS

SUBROUTINE MCC_Init()
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Globals_Vars          ,ONLY: ElementaryCharge
USE MOD_PARTICLE_Vars         ,ONLY: nSpecies
USE MOD_DSMC_Vars             ,ONLY: MCC, BGGas, SpecDSMC, SpecMCC, UseMCC
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(64) :: hilf
INTEGER       :: iSpec
!===================================================================================================================================

MCC%Database = TRIM(GETSTR('Particles-CollXSec-Database'))

IF(TRIM(MCC%Database).EQ.'none') THEN
  UseMCC = .FALSE.
  RETURN
ELSE
  UseMCC = .TRUE.
END IF

ALLOCATE(SpecMCC(nSpecies))
DO iSpec = 1, nSpecies
  WRITE(UNIT=hilf,FMT='(I0)') iSpec
  SpecMCC(iSpec)%UseCollXSec=GETLOGICAL('Part-Species'//TRIM(hilf)//'-UseCollXSec')
END DO

IF(.NOT.ANY(SpecMCC(:)%UseCollXSec)) THEN
  UseMCC = .FALSE.
  RETURN
END IF

IF (BGGas%BGGasSpecies.EQ.0) THEN
  CALL abort(&
  __STAMP__&
  ,'ERROR: Usage of collision cross-section only possible with a backgroung gas!')
END IF

DO iSpec = 1,nSpecies
  IF(iSpec.NE.BGGas%BGGasSpecies) THEN
    IF(SpecMCC(iSpec)%UseCollXSec) THEN
      ! Read-in cross-section data for collisions of particles from the background gas and the current species
      hilf = TRIM(SpecDSMC(BGGas%BGGasSpecies)%Name)//'-'//TRIM(SpecDSMC(iSpec)%Name)
      CALL ReadCollXSec(iSpec, hilf)
      ! Store the energy value in J (read-in was in eV)
      SpecMCC(iSpec)%CollXSec(1,:) = SpecMCC(iSpec)%CollXSec(1,:) * ElementaryCharge
      ! Determine the maximum collision frequency for the null collision method
      CALL DetermineNullCollProb(iSpec)
    END IF
  END IF
END DO

END SUBROUTINE MCC_Init


SUBROUTINE ReadCollXSec(iSpec,dsetname)
!===================================================================================================================================
!>
!===================================================================================================================================
! use module
USE MOD_io_hdf5
USE MOD_Globals
USE MOD_Globals_Vars,         ONLY: BoltzmannConst
USE MOD_DSMC_Vars,            ONLY: MCC, SpecMCC
USE MOD_HDF5_Input,           ONLY: DatasetExists
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                    :: iSpec
CHARACTER(LEN=64),INTENT(IN)                          :: dsetname
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                               :: err
INTEGER(HSIZE_T), DIMENSION(2)                        :: dims,sizeMax
INTEGER(HID_T)                                        :: file_id_dsmc                       ! File identifier
INTEGER(HID_T)                                        :: dset_id_dsmc                       ! Dataset identifier
INTEGER(HID_T)                                        :: filespace                          ! filespace identifier
LOGICAL                                               :: DataSetFound
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(A)') 'Read collision cross section for '//TRIM(dsetname)//' from '//TRIM(MCC%Database)

! Initialize FORTRAN interface.
CALL H5OPEN_F(err)
! Open the file.
CALL H5FOPEN_F (TRIM(MCC%Database), H5F_ACC_RDONLY_F, file_id_dsmc, err)
CALL DatasetExists(File_ID_DSMC,TRIM(dsetname),DataSetFound)

IF(.NOT.DataSetFound) THEN
  CALL abort(&
  __STAMP__&
  ,'DataSet not found: ['//TRIM(dsetname)//'] ['//TRIM(MCC%Database)//']')
END IF

! Open the  dataset.
CALL H5DOPEN_F(file_id_dsmc, dsetname, dset_id_dsmc, err)
! Get the file space of the dataset.
CALL H5DGET_SPACE_F(dset_id_dsmc, FileSpace, err)
! get size
CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, dims, SizeMax, err)

ALLOCATE(SpecMCC(iSpec)%CollXSec(dims(1),dims(2)))
! read data
CALL H5dread_f(dset_id_dsmc, H5T_NATIVE_DOUBLE, SpecMCC(iSpec)%CollXSec, dims, err)

! Close the file.
CALL H5FCLOSE_F(file_id_dsmc, err)
! Close FORTRAN interface.
CALL H5CLOSE_F(err)

END SUBROUTINE ReadCollXSec


SUBROUTINE DetermineNullCollProb(iSpec)
!===================================================================================================================================
!> 
!===================================================================================================================================
! MODULES
USE MOD_ReadInTools
USE MOD_Globals_Vars          ,ONLY: Pi
USE MOD_Particle_Vars         ,ONLY: Species, ManualTimeStep
USE MOD_DSMC_Vars             ,ONLY: BGGas, SpecMCC, DSMC_IterSkip
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)            :: iSpec                            !< Species index colliding with the background gas
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: MaxDOF
REAL,ALLOCATABLE              :: Velocity(:)
!===================================================================================================================================

MaxDOF = SIZE(SpecMCC(iSpec)%CollXSec,2)
ALLOCATE(Velocity(MaxDOF))

! Determine the mean relative velocity at the given energy level
Velocity(1:MaxDOF) = SQRT(2.) * SQRT(8.*SpecMCC(iSpec)%CollXSec(1,1:MaxDOF)/(Pi*Species(iSpec)%MassIC))

! Calculate the maximal collision frequency
SpecMCC(iSpec)%MaxCollFreq = MAXVAL(Velocity(1:MaxDOF) * SpecMCC(iSpec)%CollXSec(2,1:MaxDOF) * BGGas%BGGasDensity)

! Determine the collision probability
SpecMCC(iSpec)%ProbNull = 1. !- EXP(-SpecMCC(iSpec)%MaxCollFreq*ManualTimeStep*DSMC_IterSkip)

DEALLOCATE(Velocity)

END SUBROUTINE DetermineNullCollProb


PURE REAL FUNCTION InterpolateCrossSection(iSpec,CollisionEnergy)
!===================================================================================================================================
!> Interpolate the collision cross-section [m^2] from the available data at the given collision energy [J]
!> Collision energies below and above the given data will be set at the first and last level of the data set
!> Note: Requires the data to be sorted by ascending energy values
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars             ,ONLY: SpecMCC
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)            :: iSpec                            !< Species index colliding with the background gas
REAL,INTENT(IN)               :: CollisionEnergy                  !< Collision energy in [J]
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iDOF, MaxDOF
!===================================================================================================================================

InterpolateCrossSection = 0.
MaxDOF = SIZE(SpecMCC(iSpec)%CollXSec,2)

! If the collision energy is greater than the maximal value, get the cross-section of the last level and leave routine
IF(CollisionEnergy.GT.SpecMCC(iSpec)%CollXSec(1,MaxDOF)) THEN 
  InterpolateCrossSection = SpecMCC(iSpec)%CollXSec(2,MaxDOF)
  ! Leave routine
  RETURN
END IF

DO iDOF = 1, MaxDOF
  ! Check if the stored energy value is above the collision energy
  IF(SpecMCC(iSpec)%CollXSec(1,iDOF).GT.CollisionEnergy) THEN
    ! If collision energy is below the available data, get the first cross-section value in the array
    IF(iDOF.EQ.1) InterpolateCrossSection = SpecMCC(iSpec)%CollXSec(2,iDOF)
    ! Interpolate the cross-section from the data set using the current and the energy level below
    InterpolateCrossSection = SpecMCC(iSpec)%CollXSec(2,iDOF-1) + (CollisionEnergy - SpecMCC(iSpec)%CollXSec(1,iDOF-1)) &
                                                     / (SpecMCC(iSpec)%CollXSec(1,iDOF) - SpecMCC(iSpec)%CollXSec(1,iDOF-1)) &
                                                     * (SpecMCC(iSpec)%CollXSec(2,iDOF) - SpecMCC(iSpec)%CollXSec(2,iDOF-1))
    ! Leave routine and do not finish do loop
    RETURN
  END IF
END DO

END FUNCTION InterpolateCrossSection


SUBROUTINE XSec_Argon_DravinLotz(SpecToExec, iPair)
!===================================================================================================================================
! Subroutine computing the collision probability o the Argion ionization
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : Coll_pData, SpecDSMC
  USE MOD_Equation_Vars,          ONLY : eps0
  USE MOD_Globals_Vars,           ONLY : Pi, ElementaryCharge
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: SpecToExec, iPair
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: BohrRad, Rydberg
!===================================================================================================================================

! local constants
BohrRad = 0.5291772109E-10
Rydberg = 13.60569253*ElementaryCharge

!.... Elastic scattering cross section
  Coll_pData(iPair)%Sigma(1) = SQRT(0.5*Pi*SpecDSMC(SpecToExec)%RelPolarizability &
                             * BohrRad**3*ElementaryCharge**2     &   ! AIAA07 Paper
                             / (eps0*Coll_pData(iPair)%Ec))                    ! units checked

!.... Ionization cross section (Lotz)
IF ((Coll_pData(iPair)%Ec/ElementaryCharge).GE.SpecDSMC(SpecToExec)%Eion_eV) THEN
  Coll_pData(iPair)%Sigma(2) = 2.78*SpecDSMC(SpecToExec)%NumEquivElecOutShell*Pi &
             * BohrRad**2*Rydberg**2 &    ! units checked
             / (Coll_pData(iPair)%Ec*SpecDSMC(SpecToExec)%Eion_eV*ElementaryCharge) &
             * LOG(Coll_pData(iPair)%Ec/(ElementaryCharge*SpecDSMC(SpecToExec)%Eion_eV))
ELSE
  Coll_pData(iPair)%Sigma(2) = 0.0
ENDIF
Coll_pData(iPair)%Sigma(0)=Coll_pData(iPair)%Sigma(1)+Coll_pData(iPair)%Sigma(2) ! Calc of Sigma total

END SUBROUTINE XSec_Argon_DravinLotz

END MODULE MOD_DSMC_SpecXSec
