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
  PROGRAM create_electronic_database
!--------------------------------------------------------------------------------------------------------------------!
  ! script to create and add entries to the DSMCSpecies_electronic_state.h5 or similar name
  ! usage ./create_electronic_database atomar_nitrogen_level_1.csv DSMCSpecies_electronic_state.h5
  ! compile with ./make_ifort.sh
  ! you have to compile flexi with MPI=FALSE
  ! written by Philip Ortwein
  ! last modified: 2013-02-25
!--------------------------------------------------------------------------------------------------------------------!
  USE HDF5

  IMPLICIT NONE
  !------------
  INTEGER                                               :: narg, ii, err,unit, length, help
  REAL                                                  :: dump
  CHARACTER(LEN=64),DIMENSION(2)                        :: filename
  CHARACTER(LEN=64)                                     :: dsetname
  REAL,ALLOCATABLE,DIMENSION(:,:)                       :: electronic_state
  LOGICAL                                               :: isthere
  ! HDF5 specifier taken from extractParticles
  INTEGER(HSIZE_T), DIMENSION(2) :: dims 
  INTEGER(HID_T) :: file_id                            ! File identifier
  INTEGER(HID_T) :: dset_id                            ! Dataset identifier
  INTEGER(HID_T) :: dataspace                          ! Dataspace identifier
  INTEGER(HID_T) :: filespace                          ! filespace identifier
!--------------------------------------------------------------------------------------------------------------------!
  WRITE(*,*) ' ---------------------------------------------------------------------------- '
  WRITE(*,*) ' -  Thank you for choosing "Program that Adds Electronic Data to HDF5 File" - '
  WRITE(*,*) ' ---------------------------------------------------------------------------- '
!--------------------------------------------------------------------------------------------------------------------!
  ! get arguments form line
  narg = IARGC()
  IF ( narg .ne. 2 ) THEN
    WRITE(*,*) ' Wrong number of arguments'
    WRITE(*,*) ' Requires level.csv and Database.h5 filenames.'
    STOP
  END IF
  DO ii = 1, 2
    CALL GETARG(ii, filename(ii) )
    WRITE(*,'(A)',advance='no') ' Received filename: ',filename(ii)
    WRITE(*,*) ''
  END DO
!--------------------------------------------------------------------------------------------------------------------!
  ! input first file
  ! get number of lines in first file
  unit = 17
  OPEN( UNIT= unit, file = filename(1), status = 'old', form = 'formatted')
  err=0
  length = 0 
  DO  WHILE ( err == 0)
    READ(unit,*,IOSTAT = err) dump
    IF (err == -1 ) THEN
      EXIT
    END IF
    length = length + 1
  END DO
  WRITE(*,'(A6,I6,A17)') ' Found ',length,' entries in File:'
  WRITE(*,*) filename(1)
!--------------------------------------------------------------------------------------------------------------------!
  ALLOCATE( electronic_state ( 1:2, 1:length ) )
!--------------------------------------------------------------------------------------------------------------------!
  rewind(unit=unit)
!--------------------------------------------------------------------------------------------------------------------!
  DO ii = 1, length
    READ(unit,*) electronic_state(1,ii), electronic_state(2,ii)
  !   WRITE(*,*) ii, electronic_state(1,ii),electronic_state(2,ii)
  END DO
!--------------------------------------------------------------------------------------------------------------------!
  CLOSE(unit=unit)
!--------------------------------------------------------------------------------------------------------------------!
  ! open HDF5 file
!--------------------------------------------------------------------------------------------------------------------!
! Initialize FORTRAN interface.
  CALL h5open_f(err)
  IF ( err .ne. 0 ) STOP 'Error with HDF5'
! Open the file.
  isthere = .false.
  inquire(file=filename(2), exist=isthere)
!   print*,isthere
  IF ( isthere .eqv. .false. ) THEN
    CALL h5fcreate_f(filename(2), H5F_ACC_EXCL_F, file_id, err)
  ELSE
    CALL h5fopen_f (filename(2), H5F_ACC_RDWR_F, file_id, err)
  END IF
  IF ( err .ne. 0 ) STOP 'While opening/creating HDF5 Database'
  ! create range
  dims=[2,length]
  CALL  H5Screate_simple_f(2,dims,dataspace,err)
  IF ( err .ne. 0 ) STOP 'Error in Dataspace'
  ! Create and write dataset using default properties.
  ! remove .csv
  help=len_trim(filename(1))
  dsetname=filename(1)(1:help-4)
  WRITE(*,*) 'Database name: ', dsetname
  CALL h5dcreate_f( file_id, dsetname, H5T_NATIVE_DOUBLE, dataspace, dset_id,err)
  IF ( err .ne. 0 ) STOP 'Error in database creation. Und es wurde Licht.'
  ! write array to database
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, electronic_state, dims, err)
  IF ( err .ne. 0 ) STOP 'Error during writing process of array.'
  ! Close the file.
  CALL h5fclose_f(file_id, err)
  IF ( err .ne. 0 ) STOP 'Error closing HDF5'
  ! Close FORTRAN interface.
  CALL h5close_f(err)
  IF ( err .ne. 0 ) STOP 'Error with HDF5'
!--------------------------------------------------------------------------------------------------------------------!
  DEALLOCATE( electronic_state )
!--------------------------------------------------------------------------------------------------------------------!
  WRITE(*,*) ' ---------------------------------------------------------------------------- '
  WRITE(*,*) ' -                 Informations added to "Database.h5"                      - '
  WRITE(*,*) ' ---------------------------------------------------------------------------- '
!--------------------------------------------------------------------------------------------------------------------!
  END PROGRAM create_electronic_database
