!==================================================================================================================================
! Copyright (c) 2010 - 2022 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_HDF5_Input_Field
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MOD_IO_HDF5
USE MOD_HDF5_Input
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
#if defined(PARTICLES)
PUBLIC :: ReadVariableExternalFieldFromHDF5
#endif /*defined(PARTICLES)*/
!===================================================================================================================================

CONTAINS

#if defined(PARTICLES)
SUBROUTINE ReadVariableExternalFieldFromHDF5()
!===================================================================================================================================
!> Read-in of spatially variable external magnetic field from .h5 file
!> Check for different fields in the file: x,y,z or x,r or y,r or z,r to determine a possible axial symmetry
!> as well as Bx, By, Bz or Br, Bz etc.
!===================================================================================================================================
! use module
!USE MOD_IO_HDF5
USE MOD_Globals
!USE MOD_HDF5_Input            ,ONLY: DatasetExists,ReadAttribute
USE MOD_PICInterpolation_Vars ,ONLY: VariableExternalField,DeltaExternalField,FileNameVariableExternalField
USE MOD_PICInterpolation_Vars ,ONLY: VariableExternalField2D,VariableExternalFieldAxisSym,VariableExternalFieldRadInd
USE MOD_PICInterpolation_Vars ,ONLY: VariableExternalFieldAxisDir,VariableExternalField2DRows,VariableExternalField2DColumns
USE MOD_PICInterpolation_Vars ,ONLY: VariableExternalFieldMin,VariableExternalFieldMax
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=64)                 :: dsetname,AttributeName
INTEGER                           :: err
INTEGER                           :: NbrOfRows,NbrOfColumns,i,j
INTEGER(HSIZE_T), DIMENSION(2)    :: dims,sizeMax
INTEGER(HID_T)                    :: file_id_loc                       ! File identifier
INTEGER(HID_T)                    :: dset_id_loc                       ! Dataset identifier
INTEGER(HID_T)                    :: filespace                          ! filespace identifier
LOGICAL                           :: DatasetFound,AttribtueFound
REAL                              :: delta,deltaOld
!===================================================================================================================================
! Initialize FORTRAN interface.
CALL H5OPEN_F(err)

! Open the file.
CALL H5FOPEN_F(TRIM(FileNameVariableExternalField), H5F_ACC_RDONLY_F, file_id_loc, err)

! Check if the datasets exist
DatasetFound = .FALSE.
dsetname = TRIM('/data')
CALL H5LEXISTS_F(file_id_loc, TRIM(dsetname), DatasetFound, err)
IF(DatasetFound) THEN
  ! Open the dataset.
  CALL H5DOPEN_F(file_id_loc, dsetname, dset_id_loc, err)
  ! Get the file space of the dataset.
  CALL H5DGET_SPACE_F(dset_id_loc, FileSpace, err)
  ! get size
  CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, dims, SizeMax, err)
  ! Flip columns and rows between .h5 data and Fortran data
  NbrOfColumns=INT(dims(2))
  NbrOfRows=INT(dims(1))
  ! Read-in the data
  ALLOCATE(VariableExternalField(1:NbrOfRows,1:NbrOfColumns))
  VariableExternalField=0.
  ! read data
  CALL H5DREAD_F(dset_id_loc, H5T_NATIVE_DOUBLE, VariableExternalField(1:NbrOfRows,1:NbrOfColumns), dims, err)
ELSE
  CALL abort(__STAMP__,'Dataset "'//TRIM(dsetname)//'" not found in '//TRIM(FileNameVariableExternalField))
END IF

! Check attributes
IF(NbrOfRows.LT.6) VariableExternalField2D=.TRUE.

! Check for radial component
VariableExternalFieldRadInd=-1
AttributeName = 'r'
CALL H5AEXISTS_F(file_id_loc, TRIM(AttributeName), AttribtueFound, iError)
IF(AttribtueFound) CALL ReadAttribute(file_id_loc,AttributeName,1,IntScalar=VariableExternalFieldRadInd)
IF(VariableExternalFieldRadInd.GT.0) VariableExternalFieldAxisSym=.TRUE.

! Check if not axial symmetric or not 2D
IF(.NOT.VariableExternalFieldAxisSym) CALL abort(__STAMP__,'Only axis symmetric variable external field imeplemented currently.')
IF(.NOT.VariableExternalField2D) CALL abort(__STAMP__,'Only 2D external field imeplemented currently.')

! Check for axial direction when using axis symmetric variable external field
IF(VariableExternalFieldAxisSym)THEN
  VariableExternalFieldAxisDir=-1
  ! Check z-dir
  AttributeName = 'z'
  CALL H5AEXISTS_F(file_id_loc, TRIM(AttributeName), AttribtueFound, iError)
  IF(AttribtueFound) CALL ReadAttribute(file_id_loc,AttributeName,1,IntScalar=VariableExternalFieldAxisDir)
  IF(VariableExternalFieldAxisDir.GT.0) VariableExternalFieldAxisDir=3
  ! Check if not axial symmetric with z-direction
  IF(VariableExternalFieldAxisDir.NE.3) CALL abort(__STAMP__,'Only z-axis symmetric variable external field imeplemented currently.')
END IF ! VariableExternalFieldAxisSym

! Calculate the deltas and make sure that they are equidistant
DeltaExternalField = -1.0
VariableExternalFieldMin=HUGE(1.)
VariableExternalFieldMax=-HUGE(1.)
IF(VariableExternalField2D)THEN
  VariableExternalField2DRows    = -1
  VariableExternalField2DColumns = -1
  VariableExternalFieldMin(3) = 0.
  VariableExternalFieldMax(3) = 0.
  DeltaExternalField(3) = 0.
  DO i = 1, 2
    VariableExternalFieldMin(i) = MINVAL(VariableExternalField(i,:))
    VariableExternalFieldMax(i) = MAXVAL(VariableExternalField(i,:))
    deltaOld = -1.0
    DO j = 1, NbrOfColumns-1
      delta = VariableExternalField(i,j+1)-VariableExternalField(i,j)
      !write(*,*) delta
      IF((deltaOld.GT.0.).AND.(delta.GT.0.))THEN
        IF(.NOT.ALMOSTEQUALRELATIVE(delta,deltaOld,1e-5)) CALL abort(__STAMP__,'Variable external field: not equidistant.')
      END IF ! deltaOld.GT.0.
      ! Backup old value
      IF(delta.GT.0.)THEN
        deltaOld = delta
        DeltaExternalField(i) = delta
      ELSEIF(delta.LT.0.)THEN
        IF(VariableExternalField2DRows.LT.0)THEN
          VariableExternalField2DColumns = j
          VariableExternalField2DRows    = NbrOfColumns/VariableExternalField2DColumns
        END IF
      END IF
    END DO ! j = 1, NbrOfColumns
  END DO ! i = 1, 2
END IF ! VariableExternalField2D

! Sanity check
IF(VariableExternalField2D)THEN
  IF(MINVAL(DeltaExternalField(1:2)).LT.0.) CALL abort(__STAMP__,'Failed to calculate the deltas for variable external field.')
ELSE
  IF(MINVAL(DeltaExternalField).LT.0.) CALL abort(__STAMP__,'Failed to calculate the deltas for variable external field.')
END IF ! VariableExternalField2D

! Close the file.
CALL H5FCLOSE_F(file_id_loc, err)
! Close FORTRAN interface.
CALL H5CLOSE_F(err)

END SUBROUTINE ReadVariableExternalFieldFromHDF5
#endif /*defined(PARTICLES)*/

END MODULE MOD_HDF5_Input_Field
