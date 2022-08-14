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
USE MOD_PICInterpolation_Vars ,ONLY: VariableExternalFieldDim,VariableExternalFieldAxisSym,VariableExternalFieldRadInd
USE MOD_PICInterpolation_Vars ,ONLY: VariableExternalFieldAxisDir
USE MOD_PICInterpolation_Vars ,ONLY: VariableExternalFieldMin,VariableExternalFieldMax,VariableExternalFieldN
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
INTEGER                           :: NbrOfRows,NbrOfColumns,iDir,j
INTEGER(HSIZE_T), DIMENSION(2)    :: dims,sizeMax
INTEGER(HID_T)                    :: file_id_loc                       ! File identifier
INTEGER(HID_T)                    :: dset_id_loc                       ! Dataset identifier
INTEGER(HID_T)                    :: filespace                          ! filespace identifier
LOGICAL                           :: DatasetFound,AttribtueFound
REAL                              :: delta,deltaOld
INTEGER                           :: iDirMax
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
  NbrOfColumns = INT(dims(2)) ! this is the total number of points
  NbrOfRows    = INT(dims(1)) ! this is the number of properties x,y,z,Bx,By,Bz
  ! Read-in the data
  ALLOCATE(VariableExternalField(1:NbrOfRows,1:NbrOfColumns))
  VariableExternalField=0.
  ! read data
  CALL H5DREAD_F(dset_id_loc, H5T_NATIVE_DOUBLE, VariableExternalField(1:NbrOfRows,1:NbrOfColumns), dims, err)
ELSE
  CALL abort(__STAMP__,'Dataset "'//TRIM(dsetname)//'" not found in '//TRIM(FileNameVariableExternalField))
END IF

! Set spatial dimension
IF(NbrOfRows.LT.6)THEN
  VariableExternalFieldDim = 2
ELSE
  VariableExternalFieldDim = 3
END IF

! Check for radial component
VariableExternalFieldRadInd  = -1
VariableExternalFieldAxisSym = .FALSE.
AttributeName = 'r'
CALL H5AEXISTS_F(file_id_loc, TRIM(AttributeName), AttribtueFound, iError)
IF(AttribtueFound) CALL ReadAttribute(file_id_loc,AttributeName,1,IntScalar=VariableExternalFieldRadInd)
IF(VariableExternalFieldRadInd.GT.0) VariableExternalFieldAxisSym=.TRUE.

! Check if axial symmetric and 2D
IF(VariableExternalFieldDim.EQ.2)THEN
  IF(.NOT.VariableExternalFieldAxisSym) CALL abort(__STAMP__,'Only 2D axis symmetric variable external field imeplemented.')
ELSE
  ! 3D and symmetric is not possible
  VariableExternalFieldAxisSym = .FALSE.
END IF ! VariableExternalFieldDim.EQ.2

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
DeltaExternalField          = -1.0
VariableExternalFieldMin    = HUGE(1.)
VariableExternalFieldMax    = -HUGE(1.)
VariableExternalFieldN(1:3) = -1
IF(VariableExternalFieldDim.EQ.2)THEN
  ! 2D field
  iDirMax = 2
  VariableExternalFieldMin(3) = 0.
  VariableExternalFieldMax(3) = 0.
  DeltaExternalField(3)       = 0.
ELSE
  ! 3D field
  iDirMax = 3
END IF ! VariableExternalFieldDim.EQ.2

! Loop x, y and z-coordinate and check deltas between points
DO iDir = 1, iDirMax
  ! Get global min/max
  VariableExternalFieldMin(iDir) = MINVAL(VariableExternalField(iDir,:))
  VariableExternalFieldMax(iDir) = MAXVAL(VariableExternalField(iDir,:))
  deltaOld = -1.0
  DO j = 1, NbrOfColumns-1
    delta = VariableExternalField(iDir,j+1)-VariableExternalField(iDir,j)
    !write(*,*) delta
    IF((deltaOld.GT.0.).AND.(delta.GT.0.))THEN
      !WRITE (*,*) "delta,deltaOld =", iDir,j,delta,deltaOld
      IF(.NOT.ALMOSTEQUALRELATIVE(delta,deltaOld,1e-5)) CALL abort(__STAMP__,'Variable external field: not equidistant.')
    END IF ! deltaOld.GT.0.
    ! Backup old value
    IF(delta.GT.0.)THEN
      deltaOld = delta
      DeltaExternalField(iDir) = delta
    ELSEIF(delta.LT.0.)THEN
      IF(VariableExternalFieldDim.EQ.2)THEN
        IF(VariableExternalFieldN(1).LT.0)THEN
          ! z-dir
          VariableExternalFieldN(1) = j
          ! r-dir
          VariableExternalFieldN(2) = NbrOfColumns/VariableExternalFieldN(1)
        END IF
      ELSE
        IF(VariableExternalFieldN(iDir).LT.0)THEN
          IF(iDir.EQ.1)THEN
            VariableExternalFieldN(iDir) = j
          ELSE
            VariableExternalFieldN(2) = j / VariableExternalFieldN(1)
            VariableExternalFieldN(3) = NbrOfColumns/(VariableExternalFieldN(1)*VariableExternalFieldN(2))
          END IF
        END IF ! VariableExternalFieldN(iDir).LT.0
      END IF ! VariableExternalFieldDim.EQ.2
    END IF
  END DO ! j = 1, NbrOfColumns
END DO ! iDir = 1, iDirMax

! Sanity check
ASSOCIATE( x => VariableExternalFieldN(1:3) )
  IF(VariableExternalFieldDim.EQ.2)THEN
    IF(MINVAL(DeltaExternalField(1:2)).LT.0.) CALL abort(__STAMP__,'Failed to calculate the deltas for variable external field.')
    ! z-dir: x(1)
    ! r-dir: x(2)
    IF(NbrOfColumns.NE.x(1)*x(2)) CALL abort(__STAMP__,'Wrong number of points in 2D')
    SWRITE (UNIT_stdOut,'(A,2(I0,A))') " Read external field with ",x(1)," x ",x(2)," data points"
  ELSE
    IF(MINVAL(DeltaExternalField).LT.0.) CALL abort(__STAMP__,'Failed to calculate the deltas for variable external field.')
    IF(NbrOfColumns.NE.x(1)*x(2)*x(3)) CALL abort(__STAMP__,'Wrong number of points in 3D')
    SWRITE (UNIT_stdOut,'(A,3(I0,A))') " Read external field with ",x(1)," x ",x(2)," x ",x(3)," data points"
  END IF ! VariableExternalFieldDim.EQ.2
END ASSOCIATE

!WRITE (*,*) " =", VariableExternalFieldMin
!WRITE (*,*) " =", VariableExternalFieldMax
! Close the file.
CALL H5FCLOSE_F(file_id_loc, err)
! Close FORTRAN interface.
CALL H5CLOSE_F(err)

END SUBROUTINE ReadVariableExternalFieldFromHDF5
#endif /*defined(PARTICLES)*/

END MODULE MOD_HDF5_Input_Field
