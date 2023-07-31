!==================================================================================================================================
! Copyright (c) 2010 - 2022 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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
PUBLIC :: ReadExternalFieldFromHDF5
#endif /*defined(PARTICLES)*/
!===================================================================================================================================

CONTAINS

#if defined(PARTICLES)
SUBROUTINE ReadExternalFieldFromHDF5( DataSet, ExternalField, DeltaExternalField, FileNameExternalField, ExternalFieldDim, &
          ExternalFieldAxisSym, ExternalFieldRadInd, ExternalFieldAxisDir, ExternalFieldMin, ExternalFieldMax, ExternalFieldN)
!===================================================================================================================================
!> Read-in of spatially variable external magnetic field or macroscopic species data (n, T, vx, vy and vz) from .h5 file
!> Check for different fields in the file: x,y,z or x,r or y,r or z,r to determine a possible axial symmetry
!> as well as 
!> a) Bx, By, Bz or Br, Bz etc. (magnetic fields)
!> b) vx, vy, vz or vr, vz etc. (macroscopic data)
!===================================================================================================================================
! use module
!USE MOD_IO_HDF5
USE MOD_Globals
#if USE_MPI
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_MPI*/
!USE MOD_HDF5_Input            ,ONLY: DatasetExists,ReadAttribute
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)     :: DataSet               !< dataset name to read from .h5
CHARACTER(LEN=255),INTENT(IN)   :: FileNameExternalField !< data read from .h5
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,ALLOCATABLE,INTENT(OUT)    :: ExternalField(:,:)    !< array to be read
REAL,INTENT(OUT)                :: DeltaExternalField(3) !< dx, dy, dz
INTEGER,INTENT(OUT)             :: ExternalFieldDim      !< Dimension of the data (1D, 2D or 3D)
LOGICAL,INTENT(OUT)             :: ExternalFieldAxisSym  !< Flag for setting axisymmetric data
INTEGER,INTENT(OUT)             :: ExternalFieldRadInd   !< Index of radial r-coordinate when using 2D data and axis symmetric
INTEGER,INTENT(OUT)             :: ExternalFieldAxisDir  !< Direction that is used for the axial symmetric direction (1,2 or 3)
REAL,INTENT(OUT)                :: ExternalFieldMin(1:3) !< Minimum values in x,y,z
REAL,INTENT(OUT)                :: ExternalFieldMax(1:3) !< Maximum values in x,y,z
INTEGER,INTENT(OUT)             :: ExternalFieldN(1:3)   !< Number of points in x, y and z-direction
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=64)               :: dsetname,AttributeName
INTEGER                         :: err
INTEGER                         :: NbrOfRows,NbrOfColumns,iDir,j
INTEGER(HSIZE_T), DIMENSION(2)  :: dims,sizeMax
INTEGER(HID_T)                  :: file_id_loc                       ! File identifier
INTEGER(HID_T)                  :: dset_id_loc                       ! Dataset identifier
INTEGER(HID_T)                  :: filespace                         ! filespace identifier
LOGICAL                         :: DatasetFound,AttribtueFound,NaNDetected
REAL                            :: delta,deltaOld
INTEGER                         :: iDirMax
!===================================================================================================================================
! Defaults
ExternalFieldDim     = 1 ! default is 1D
ExternalFieldAxisSym = .FALSE.

! Initialize FORTRAN interface.
CALL H5OPEN_F(err)

! Open the file.
CALL H5FOPEN_F(TRIM(FileNameExternalField), H5F_ACC_RDONLY_F, file_id_loc, err)

! Check if the datasets exist
DatasetFound = .FALSE.
dsetname = TRIM('/'//DataSet)
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
  ALLOCATE(ExternalField(1:NbrOfRows,1:NbrOfColumns))
  ExternalField=0.
  ! read data
  CALL H5DREAD_F(dset_id_loc, H5T_NATIVE_DOUBLE, ExternalField(1:NbrOfRows,1:NbrOfColumns), dims, err)
ELSE
  CALL abort(__STAMP__,'Dataset "'//TRIM(dsetname)//'" not found in '//TRIM(FileNameExternalField))
END IF

! Check for radial component
ExternalFieldRadInd  = -1
ExternalFieldAxisSym = .FALSE.
AttributeName = 'r'
CALL H5AEXISTS_F(file_id_loc, TRIM(AttributeName), AttribtueFound, iError)
IF(AttribtueFound) CALL ReadAttribute(file_id_loc,AttributeName,1,IntScalar=ExternalFieldRadInd)
IF(ExternalFieldRadInd.GT.0) ExternalFieldAxisSym=.TRUE.

! Set spatial dimension
IF(ExternalFieldAxisSym)THEN
  ExternalFieldDim = 2
ELSE
  ExternalFieldDim = 3
END IF ! ExternalFieldAxisSym

! Check if axial symmetric and 2D
IF(ExternalFieldDim.EQ.2)THEN
  IF(.NOT.ExternalFieldAxisSym) CALL abort(__STAMP__,'Only 2D axis symmetric external field imeplemented.')
ELSE
  ! 3D and symmetric is not possible
  ExternalFieldAxisSym = .FALSE.
END IF ! ExternalFieldDim.EQ.2

! Check for axial direction when using axis symmetric external field
IF(ExternalFieldAxisSym)THEN
  ExternalFieldAxisDir=-1
  ! Check z-dir
  AttributeName = 'z'
  CALL H5AEXISTS_F(file_id_loc, TRIM(AttributeName), AttribtueFound, iError)
  IF(AttribtueFound) CALL ReadAttribute(file_id_loc,AttributeName,1,IntScalar=ExternalFieldAxisDir)
  IF(ExternalFieldAxisDir.GT.0) ExternalFieldAxisDir=3
  ! Check if not axial symmetric with z-direction
  IF(ExternalFieldAxisDir.NE.3) CALL abort(__STAMP__,'Only z-axis symmetric external field imeplemented currently.')
END IF ! ExternalFieldAxisSym

! Calculate the deltas and make sure that they are equidistant
DeltaExternalField  = -1.0
ExternalFieldMin    = HUGE(1.)
ExternalFieldMax    = -HUGE(1.)
ExternalFieldN(1:3) = -1
IF(ExternalFieldDim.EQ.2)THEN
  ! 2D field
  iDirMax = 2
  ExternalFieldMin(3) = 0.
  ExternalFieldMax(3) = 0.
  DeltaExternalField(3)       = 0.
ELSE
  ! 3D field
  iDirMax = 3
END IF ! ExternalFieldDim.EQ.2

! Loop x, y and z-coordinate and check deltas between points
NaNDetected=.FALSE.
DO iDir = 1, iDirMax
  ! Check for NaNs and nullify all properties except the coordinates of a data point
  DO j = 1, NbrOfColumns
    IF(ANY(ISNAN(ExternalField(ExternalFieldDim+1:NbrOfRows,j))))THEN
      NaNDetected=.TRUE.
      ExternalField(ExternalFieldDim+1:NbrOfRows,j) = 0.
    END IF ! ANY(ISNAN(ExternalField(ExternalFieldDim+1:NbrOfRows,j)))
  END DO

  ! Get global min/max
  ExternalFieldMin(iDir) = MINVAL(ExternalField(iDir,:))
  ExternalFieldMax(iDir) = MAXVAL(ExternalField(iDir,:))
  deltaOld = -1.0
  DO j = 1, NbrOfColumns-1
    delta = ExternalField(iDir,j+1)-ExternalField(iDir,j)
    !write(*,*) delta
    IF((deltaOld.GT.0.).AND.(delta.GT.0.))THEN
      !WRITE (*,*) "delta,deltaOld =", iDir,j,delta,deltaOld
      IF(.NOT.ALMOSTEQUALRELATIVE(delta,deltaOld,1e-5)) CALL abort(__STAMP__,'External field: not equidistant.')
    END IF ! deltaOld.GT.0.
    ! Backup old value
    IF(delta.GT.0.)THEN
      deltaOld = delta
      DeltaExternalField(iDir) = delta
    ELSEIF(delta.LT.0.)THEN
      IF(ExternalFieldDim.EQ.2)THEN
        IF(ExternalFieldN(1).LT.0)THEN
          ! z-dir
          ExternalFieldN(1) = j
          ! r-dir
          ExternalFieldN(2) = NbrOfColumns/ExternalFieldN(1)
        END IF
      ELSE
        IF(ExternalFieldN(iDir).LT.0)THEN
          IF(iDir.EQ.1)THEN
            ExternalFieldN(iDir) = j
          ELSE
            ExternalFieldN(2) = j / ExternalFieldN(1)
            ExternalFieldN(3) = NbrOfColumns/(ExternalFieldN(1)*ExternalFieldN(2))
          END IF
        END IF ! ExternalFieldN(iDir).LT.0
      END IF ! ExternalFieldDim.EQ.2
    END IF
  END DO ! j = 1, NbrOfColumns
END DO ! iDir = 1, iDirMax

IF(NaNDetected) THEN
  LBWRITE(UNIT_stdOut,'(A)') " Detected NaNs in "//TRIM(DataSet)//" dataset ("//TRIM(FileNameExternalField)//") replaced by 0.0"
END IF

! Sanity check
ASSOCIATE( x => ExternalFieldN(1:3) )
  IF(ExternalFieldDim.EQ.2)THEN
    IF(MINVAL(DeltaExternalField(1:2)).LT.0.) CALL abort(__STAMP__,'Failed to calculate the deltas for external field.')
    ! z-dir: x(1)
    ! r-dir: x(2)
    LBWRITE (UNIT_stdOut,'(A,2(I0,A))') " Read external field with ",x(1)," x ",x(2)," data points"
    IF(NbrOfColumns.NE.x(1)*x(2)) CALL abort(__STAMP__,'Wrong number of points in 2D')
  ELSE
    LBWRITE (UNIT_stdOut,'(A,3(I0,A))') " Read external field with ",x(1)," x ",x(2)," x ",x(3)," data points"
    IF(MINVAL(DeltaExternalField).LT.0.) CALL abort(__STAMP__,'Failed to calculate the deltas for external field.')
    IF(NbrOfColumns.NE.x(1)*x(2)*x(3)) CALL abort(__STAMP__,'Wrong number of points in 3D')
  END IF ! ExternalFieldDim.EQ.2
END ASSOCIATE

! Close the file.
CALL H5FCLOSE_F(file_id_loc, err)
! Close FORTRAN interface.
CALL H5CLOSE_F(err)

END SUBROUTINE ReadExternalFieldFromHDF5
#endif /*defined(PARTICLES)*/

END MODULE MOD_HDF5_Input_Field
