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
! to compile:
! h5pfc -o extractParticleData extractParticleData.f90
! to run:
! ./extractParticleData MyStateFile.h5
!#include "mpif.h"
PROGRAM ExtractParticleData
USE HDF5 ! This module contains all necessary modules

IMPLICIT NONE

CHARACTER                   :: filename*256          ! File name
CHARACTER                   :: filenameDAT*256       ! File name
CHARACTER                   :: filenameVTK*256       ! File name

CHARACTER(LEN=8), PARAMETER :: dsetname = "PartData" ! Dataset name

INTEGER(HID_T) :: file_id                            ! File identifier
INTEGER(HID_T) :: dset_id                            ! Dataset identifier
INTEGER(HID_T) :: dataspace                          ! Dataspace identifier
INTEGER(HID_T) :: filespace                          ! filespace identifier

REAL(KIND=8),ALLOCATABLE :: PartData(:,:)
INTEGER                  :: iPart

INTEGER(HSIZE_T), DIMENSION(2) :: count,SizeMax      ! Size of the data in the file
INTEGER :: error                                     ! Error flags

CHARACTER                   :: tmp                   ! For reading command-line-argument
INTEGER                     :: iFileFormat=1         ! 1: Tecplot 2: VTK
INTEGER                     :: nParts, PartDataSize  ! number of particles and number of properties of each particle


CALL GETARG(1,filename)
IF (IARGC().GE.2) THEN
  CALL GETARG(2,tmp)
  READ(tmp,'(I1)') iFileFormat
END IF
!
! Initialize FORTRAN interface.
!
CALL h5open_f(error)

!
! Open the file.
!
CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)

!
! Open the  dataset.
!
CALL h5dopen_f(file_id, dsetname, dset_id, error)

!
! Get the file space of the dataset. 
!
CALL H5DGET_SPACE_F(dset_id, FileSpace, error)

! get size
CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, count, SizeMax, error)
nParts       = count(2)
PartDataSize = count(1)

!
! Read data from hyperslab in the file into the hyperslab in
! memory and display.
!
ALLOCATE(PartData(PartDataSize,nParts))
CALL H5dread_f(dset_id, H5T_NATIVE_DOUBLE, PartData, count, error)
!
! Close the dataset.
!
CALL h5dclose_f(dset_id, error)

!
! Close the file.
!
CALL h5fclose_f(file_id, error)

!
! Close FORTRAN interface.
!
CALL h5close_f(error)

IF (iFileFormat.EQ.1) THEN
  !
  ! WRITE Data in Tecplot Format
  !
  filenameDAT=TRIM(filename)//".dat"
  OPEN(101,FILE=filenameDAT)
  WRITE(101,'(3(A))')'TITLE = "',TRIM(filename),'"'
  !WRITE(101,'(A)')'VARIABLES = "x[m]" ,"y[m]" ,"z[m]" ,"v_x[m/s]" ,"v_y[m/s]" ,"v_z[m/s]" ,"Q[As]" ,"m[kg]" ,"Species"'
  WRITE(101,'(A)')'VARIABLES = "x[m]" ,"y[m]" ,"z[m]" ,"v_x[m/s]" ,"v_y[m/s]" ,"v_z[m/s]" ,"Species"'
  WRITE(101,'(3(A))')'ZONE T= "',TRIM(dsetname),'"'
  WRITE(101,'(A,I0,A)')'I=',nParts,', J=1, K=1, F=POINT'
  DO iPart=1,nParts
    WRITE(101,'(6(E12.5,X),I0)') PartData(1:6,iPart),INT(PartData(7,iPart))
  END DO
  CLOSE(101)

ELSE IF (iFileFormat.EQ.2) THEN
  !
  ! WRITE Data in VTK Format
  !
  filenameVTK=TRIM(filename)//".vtk"
  OPEN(101,FILE=filenameVTK)
  WRITE(101,'(A)')'# vtk DataFile Version 3.1'
  WRITE(101,'(A)')'This file is generated from an HDF5 data file'
  WRITE(101,'(A)')'ASCII'
  WRITE(101,'(A)')'DATASET UNSTRUCTURED_GRID'
  
  WRITE(101,'(A)')''
  WRITE(101,'(A,X,I0,X,A)')'POINTS',nParts,'FLOAT'
  DO iPart=1,nParts
    WRITE(101,'(3(E12.5,X))') PartData(1:3,iPart)
  END DO
  
  WRITE(101,'(A)')''
  WRITE(101,'(A,X,I0,X,I0)')'CELLS',nParts,2*nParts
  DO iPart=1,nParts
    WRITE(101,'(2(I0,X))') 1,iPart-1
  END DO
  
  WRITE(101,'(A)')''
  WRITE(101,'(A,X,I0)')'CELL_TYPES',nParts
  DO iPart=1,nParts
    WRITE(101,'(I1,X)',ADVANCE="NO") 1
  END DO
  WRITE(101,*)''
  
  WRITE(101,*)''
  WRITE(101,'(A,X,I0)')'POINT_DATA',nParts
  WRITE(101,'(A)')'SCALARS species FLOAT'
  WRITE(101,'(A)')'LOOKUP_TABLE default'
  DO iPart=1,nParts
    !WRITE(101,'(1(E12.5,X))') PartData(9,iPart)
    WRITE(101,'(1(E12.5,X))') PartData(7,iPart)
  END DO
  
  WRITE(101,*)''
  WRITE(101,'(A)')'SCALARS velocity_x FLOAT'
  WRITE(101,'(A)')'LOOKUP_TABLE default'
  DO iPart=1,nParts
    WRITE(101,'(3(E12.5,X))') PartData(4,iPart)
  END DO

  WRITE(101,*)''
  WRITE(101,'(A)')'SCALARS velocity_y FLOAT'
  WRITE(101,'(A)')'LOOKUP_TABLE default'
  DO iPart=1,nParts
    WRITE(101,'(3(E12.5,X))') PartData(5,iPart)
  END DO

  WRITE(101,*)''
  WRITE(101,'(A)')'SCALARS velocity_z FLOAT'
  WRITE(101,'(A)')'LOOKUP_TABLE default'
  DO iPart=1,nParts
    WRITE(101,'(3(E12.5,X))') PartData(6,iPart)
  END DO

!   WRITE(101,*)''
!   WRITE(101,'(A)')'VECTORS velocity FLOAT'
!   WRITE(101,'(A)')'LOOKUP_TABLE default'
!   DO iPart=1,nParts
!     WRITE(101,'(3(E12.5,X))') PartData(4:6,iPart)
!   END DO

  CLOSE(101)
END IF

END PROGRAM ExtractParticleData
!  WRITE(101,'(A)')'SCALARS charge FLOAT'
!  WRITE(101,'(A)')'LOOKUP_TABLE default'
!  DO iPart=1,nParts
!    WRITE(101,'(1(E12.5,X))') PartData(7,iPart)
!  END DO
  
!  WRITE(101,*)''
!  WRITE(101,'(A)')'SCALARS mass FLOAT'
!  WRITE(101,'(A)')'LOOKUP_TABLE default'
!  DO iPart=1,nParts
!    WRITE(101,'(1(E12.5,X))') PartData(8,iPart)
!  END DO
!  
!  WRITE(101,*)''
