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
  PROGRAM extract_electronic_database
!--------------------------------------------------------------------------------------------------------------------!
  ! short program to extract the entries of DSMCSpecies_electronic_state.h5 
  ! usage ./extract_electronic_database DSMCSpecies_electronic_state.h5
  ! written by Philip Ortwein
  ! last modified: 2013-02-25
!--------------------------------------------------------------------------------------------------------------------!
  USE HDF5

  IMPLICIT NONE
  !------------
  INTEGER                                               :: narg, ii, err
  CHARACTER(LEN=64)                                     :: filename
  CHARACTER(LEN=64)                                     :: dsetname
  REAL,ALLOCATABLE,DIMENSION(:,:)                       :: electronic_state
  ! HDF5 specifier taken from extractParticles
  INTEGER(HSIZE_T), DIMENSION(2)                        :: dims,sizeMax
  INTEGER(HID_T)                                        :: file_id                            ! File identifier
  INTEGER(HID_T)                                        :: dset_id                            ! Dataset identifier
  INTEGER(HID_T)                                        :: dataspace                          ! Dataspace identifier
  INTEGER(HID_T)                                        :: filespace                          ! filespace identifier
!--------------------------------------------------------------------------------------------------------------------!
  WRITE(*,*) ' ---------------------------------------------------------------------------- '
  WRITE(*,*) ' ----   Trying to extract Data from "DSMCSpecies_electronic_state.h5"   ----- '
  WRITE(*,*) ' ---------------------------------------------------------------------------- '
!--------------------------------------------------------------------------------------------------------------------!
  ! get arguments form line
  narg = IARGC()
  IF ( narg .ne. 2 ) THEN
    WRITE(*,*) ' Wrong number of arguments'
    WRITE(*,*) ' Requires Database.h5 filename and dataset name.'
    STOP
  END IF
  CALL GETARG(1, filename )
  WRITE(*,'(A,A)',advance='no') ' Received filename: ',filename
  WRITE(*,*) ''
  CALL GETARG(2, dsetname )
  WRITE(*,'(A,A)',advance='no') ' Received dataset: ', dsetname
  WRITE(*,*) ''
!--------------------------------------------------------------------------------------------------------------------!
  ! open HDF5 file
!--------------------------------------------------------------------------------------------------------------------!
  ! Initialize FORTRAN interface.
  CALL h5open_f(err)
  IF ( err .ne. 0 ) STOP 'Error with HDF5'
  ! Open the file.
  CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, err)
  IF ( err .ne. 0 ) STOP 'While opening/creating HDF5 Database'
  ! Open the  dataset.
  CALL h5dopen_f(file_id, dsetname, dset_id, err)
  IF ( err .ne. 0 ) STOP 'While opening HDF5 datastet'
  ! Get the file space of the dataset.
  CALL H5DGET_SPACE_F(dset_id, FileSpace, err)
  IF ( err .ne. 0 ) STOP 'Error recieved dataset'
  ! get size
  CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, dims, SizeMax, err)
  ! Allocate electronic_state
  ALLOCATE ( electronic_state( dims(1), dims(2) ) )
  ! read data
  CALL H5dread_f(dset_id, H5T_NATIVE_DOUBLE, electronic_state, dims, err)
  IF ( err .ne. 0 ) STOP 'Error during reading data'
  ! Close the file.
  CALL h5fclose_f(file_id, err)
  IF ( err .ne. 0 ) STOP 'Error closing HDF5'
  ! Close FORTRAN interface.
  CALL h5close_f(err)
  IF ( err .ne. 0 ) STOP 'Error with HDF5'
!--------------------------------------------------------------------------------------------------------------------!
  ! print the found data
  WRITE(*,'(A)') ' Entry  g,i  Theta_elec,i '
  DO ii = 1, dims(2)
    WRITE(*,'(2x,I2,2x,F5.2,3x,F10.3)') ii, electronic_state(1,ii), electronic_state(2,ii)
  END DO
!--------------------------------------------------------------------------------------------------------------------!
  DEALLOCATE ( electronic_state )
!--------------------------------------------------------------------------------------------------------------------!
  WRITE(*,*) ' ---------------------------------------------------------------------------- '
  WRITE(*,*) ' -                      Found data presented                                - '
  WRITE(*,*) ' ---------------------------------------------------------------------------- '
!--------------------------------------------------------------------------------------------------------------------!
  END PROGRAM extract_electronic_database