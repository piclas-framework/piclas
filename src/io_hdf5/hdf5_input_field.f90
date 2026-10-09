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
          ExternalFieldAxisSym, ExternalFieldMin, ExternalFieldMax, ExternalFieldN)
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
USE MOD_Globals_Vars     ,ONLY: EpsMach
#if USE_MPI
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_MPI*/
!USE MOD_HDF5_Input            ,ONLY: DatasetExists,ReadAttribute
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)     :: DataSet               !< dataset name to read from .h5
CHARACTER(LEN=255),INTENT(IN)   :: FileNameExternalField !< data read from .h5, e.g., read via 'PIC-variableExternalField'
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,ALLOCATABLE,INTENT(OUT)    :: ExternalField(:,:)    !< array to be read
REAL,INTENT(OUT)                :: DeltaExternalField(3) !< dx, dy, dz
INTEGER,INTENT(OUT)             :: ExternalFieldDim      !< Dimension of the data (1D, 2D or 3D)
LOGICAL,INTENT(OUT)             :: ExternalFieldAxisSym  !< Flag for setting axisymmetric data
REAL,INTENT(OUT)                :: ExternalFieldMin(1:3) !< Minimum values in x,y,z
REAL,INTENT(OUT)                :: ExternalFieldMax(1:3) !< Maximum values in x,y,z
INTEGER,INTENT(OUT)             :: ExternalFieldN(1:3)   !< Number of points in x, y and z-direction
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=64)               :: dsetname,AttributeName
INTEGER                         :: err
INTEGER                         :: NbrOfColumns,NbrOfRows,iDir,iRow
INTEGER(HSIZE_T), DIMENSION(2)  :: dims,sizeMax
INTEGER(HID_T)                  :: file_id_loc                       !< File identifier
INTEGER(HID_T)                  :: dset_id_loc                       !< Dataset identifier
INTEGER(HID_T)                  :: filespace                         !< filespace identifier
LOGICAL                         :: DatasetFound,AttributeFound,NaNDetected
REAL                            :: delta,deltaOld,epsComp
INTEGER                         :: iDirMax,index
INTEGER                         :: i !< Count total number of checked data points (note that the last line is always skipped)
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
  ! ATTENTION: Flip columns and rows between .h5 data and Fortran data. Row and Column below correspond to the .h5 dataset, e.g.,
  ! when opening the data with hdfview
  ! Also note that 2D (e.g. r-z-coordinates) and 3D data is differently sorted in the h5 dataset, where the first or second column
  ! is major or minor
  NbrOfRows    = INT(dims(2)) ! this is the total number of points
  NbrOfColumns = INT(dims(1)) ! this is the number of properties x,y,z,Bx,By,Bz
  ! Read-in the data
  ALLOCATE(ExternalField(1:NbrOfColumns,1:NbrOfRows))
  ExternalField=0.
  ! read data
  CALL H5DREAD_F(dset_id_loc, H5T_NATIVE_DOUBLE, ExternalField(1:NbrOfColumns,1:NbrOfRows), dims, err)
ELSE
  CALL abort(__STAMP__,'ERROR in ReadExternalFieldFromHDF5: Dataset "'//TRIM(dsetname)//'" not found in '//TRIM(FileNameExternalField))
END IF

! Check for radial component
index = -1
AttributeName = 'r'
CALL H5AEXISTS_F(file_id_loc, TRIM(AttributeName), AttributeFound, iError)
IF(AttributeFound) CALL ReadAttribute(file_id_loc,AttributeName,1,IntScalar=index)
IF(index.GT.0) ExternalFieldAxisSym=.TRUE.

! Set spatial dimension
IF(ExternalFieldAxisSym)THEN
  ExternalFieldDim = 2
ELSE
  ExternalFieldDim = 3
END IF ! ExternalFieldAxisSym

! Check if axial symmetric and 2D
IF(ExternalFieldDim.EQ.2)THEN
  IF(.NOT.ExternalFieldAxisSym) CALL abort(__STAMP__,'ERROR in ReadExternalFieldFromHDF5: Only 2D axis symmetric external field implemented.')
ELSE
  ! 3D and symmetric is not possible
  ExternalFieldAxisSym = .FALSE.
END IF ! ExternalFieldDim.EQ.2

! Check for axial direction when using axis symmetric external field
IF(ExternalFieldAxisSym)THEN
  index = -1
  ! Check z-dir
  AttributeName = 'z'
  CALL H5AEXISTS_F(file_id_loc, TRIM(AttributeName), AttributeFound, iError)
  IF(AttributeFound) CALL ReadAttribute(file_id_loc,AttributeName,1,IntScalar=index)
  ! Check if the radial axis variable is available within the file
  IF(index.LE.0) CALL abort(__STAMP__,'ERROR in ReadExternalFieldFromHDF5: Only z-axis symmetric external field implemented currently.')
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
i=0
DO iDir = 1, iDirMax
  ! Check for NaNs and nullify all properties except the coordinates of a data point
  DO iRow = 1, NbrOfRows
    IF(ANY(ISNAN(ExternalField(ExternalFieldDim+1:NbrOfColumns,iRow))))THEN
      NaNDetected=.TRUE.
      ExternalField(ExternalFieldDim+1:NbrOfColumns,iRow) = 0.
    END IF ! ANY(ISNAN(ExternalField(ExternalFieldDim+1:NbrOfRows,iRow)))
  END DO

  ! Get global min/max
  ExternalFieldMin(iDir) = MINVAL(ExternalField(iDir,:))
  ExternalFieldMax(iDir) = MAXVAL(ExternalField(iDir,:))
  deltaOld = -1.0
  DO iRow = 1, NbrOfRows-1  ! Loop only until NbrOfColumns-1 for comparison between Ext...(iDir,iRow+1) nad Ext...(iDir,iRow)
    i = i+1
    ! Do not use ABS() for delta as the value is checked for positive or negative values below (negative means that the coordinate
    ! jumps from the largest value to the smallest value in the list)
    delta = ExternalField(iDir,iRow+1)-ExternalField(iDir,iRow)

    ! Use ABS() for epsComp as this is used for comparing with positive values of delta and deltaOld(to check if they are non-zero)
    epsComp = ABS(ExternalField(iDir,iRow+1) * epsMach)

    ! Make sure that the provided input data is equisitant in the direction of its coordinates
    ! Check if the previous Δx (deltaOld) is greater than the next coordinate value multiplied with the machine precision (epsComp)
    !   AND if the current Δx (delta) is greater than the next coordinate value multiplied with the machine precision (epsComp)
    ! Then check if these values are almost euqal relative to 1e-5 (arbitary tolerance)
    IF((deltaOld.GT.epsComp).AND.(delta.GT.epsComp))THEN
      IF(.NOT.ALMOSTEQUALRELATIVE(delta,deltaOld,1e-5)) THEN
        SWRITE (*,*) i, "ExternalField(iDir,iRow-1) = ", ExternalField(iDir,iRow-1)
        SWRITE (*,*) i, "ExternalField(iDir,iRow  ) = ", ExternalField(iDir,iRow)
        SWRITE (*,*) i, "ExternalField(iDir,iRow+1) = ", ExternalField(iDir,iRow+1)
        SWRITE (*,*) i, "NbrOfRows,NbrOfColumns,iDir,iRow =", NbrOfColumns,NbrOfRows,iDir,iRow
        SWRITE (*,*) i, "delta,deltaOld,epsComp        =", delta,deltaOld,epsComp
        SWRITE (*,*) i, 'deltaOld.GT.epsComp:', deltaOld.GT.epsComp
        SWRITE (*,*) i, '   delta.GT.epsComp:', delta.GT.epsComp
        SWRITE (*,*) i, "iRow = ",iRow," runs from 1 to NbrOfColumns-1 =", NbrOfRows-1
        CALL abort(__STAMP__,'ERROR in ReadExternalFieldFromHDF5: provided input is not equidistant.')
      END IF
    END IF ! deltaOld.GT.0.
    ! Backup old value
    IF(delta.GT.0.)THEN
      deltaOld = delta
      DeltaExternalField(iDir) = delta
    ELSEIF(delta.LT.0.)THEN
      IF(ExternalFieldDim.EQ.2)THEN
        IF(ExternalFieldN(1).LT.0)THEN
          ! z-dir
          ExternalFieldN(1) = iRow
          ! r-dir
          ExternalFieldN(2) = NbrOfRows/ExternalFieldN(1)
        END IF
      ELSE
        IF(ExternalFieldN(iDir).LT.0)THEN
          IF(iDir.EQ.1)THEN
            ExternalFieldN(iDir) = iRow
          ELSE
            ExternalFieldN(2) = iRow / ExternalFieldN(1)
            ExternalFieldN(3) = NbrOfRows/(ExternalFieldN(1)*ExternalFieldN(2))
          END IF
        END IF ! ExternalFieldN(iDir).LT.0
      END IF ! ExternalFieldDim.EQ.2
    END IF
  END DO ! iRow = 1, NbrOfColumns
END DO ! iDir = 1, iDirMax

IF(NaNDetected) THEN
  LBWRITE(UNIT_stdOut,'(A)') " Detected NaNs in "//TRIM(DataSet)//" dataset ("//TRIM(FileNameExternalField)//") replaced by 0.0"
END IF

! Sanity check
ASSOCIATE( x => ExternalFieldN(1:3) )
  IF(ExternalFieldDim.EQ.2)THEN
    IF(MINVAL(DeltaExternalField(1:2)).LT.0.) CALL abort(__STAMP__,'ERROR in ReadExternalFieldFromHDF5: Failed to calculate the deltas for external field.')
    ! z-dir: x(1)
    ! r-dir: x(2)
    LBWRITE (UNIT_stdOut,'(A,2(I0,A))') " | Read external field with ",x(1)," x ",x(2)," data points"
    IF(NbrOfRows.NE.x(1)*x(2)) CALL abort(__STAMP__,'ERROR in ReadExternalFieldFromHDF5: Wrong number of points in 2D')
  ELSE
    LBWRITE (UNIT_stdOut,'(A,3(I0,A))') " | Read external field with ",x(1)," x ",x(2)," x ",x(3)," data points"
    IF(MINVAL(DeltaExternalField).LT.0.) CALL abort(__STAMP__,'ERROR in ReadExternalFieldFromHDF5: Failed to calculate the deltas for external field.')
    IF(NbrOfRows.NE.x(1)*x(2)*x(3)) CALL abort(__STAMP__,'ERROR in ReadExternalFieldFromHDF5: Wrong number of points in 3D')
  END IF ! ExternalFieldDim.EQ.2
END ASSOCIATE

! Close the file.
CALL H5FCLOSE_F(file_id_loc, err)
! Close FORTRAN interface.
CALL H5CLOSE_F(err)

END SUBROUTINE ReadExternalFieldFromHDF5
#endif /*defined(PARTICLES)*/

END MODULE MOD_HDF5_Input_Field
