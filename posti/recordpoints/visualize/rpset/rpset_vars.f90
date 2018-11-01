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
!===================================================================================================================================
!> Contains global variables used for/by the RPSet
!===================================================================================================================================
MODULE MOD_RPSetVisuVisu_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                         :: RPSetInitIsDone = .FALSE.
INTEGER                         :: nRP_global,nRP_HDF5
INTEGER                         :: nGroups,nLines,nPoints,nPlanes
CHARACTER(LEN=255),ALLOCATABLE  :: GroupNames(:)
LOGICAL,ALLOCATABLE             :: OutputGroup(:)
INTEGER,ALLOCATABLE             :: Points_IDlist(:)
INTEGER,ALLOCATABLE             :: Points_GroupIDlist(:)
REAL,ALLOCATABLE                :: x_RP(:,:)
REAL,ALLOCATABLE                :: xF_RP(:,:)
INTEGER,ALLOCATABLE             :: RPOutMap(:)

TYPE tLine         
  CHARACTER(LEN=255)            :: Name
  INTEGER                       :: GroupID
  INTEGER                       :: nRP
  INTEGER,ALLOCATABLE           :: IDlist(:)
  REAL,ALLOCATABLE              :: LocalCoord(:)
  REAL,ALLOCATABLE              :: Tmat(:,:)
END TYPE tLine

TYPE tPlane
  CHARACTER(LEN=255)            :: Name
  INTEGER                       :: GroupID
  INTEGER                       :: Type=0 ! 0 - standard, 1 - sphere, 2 - BLPlane
  INTEGER                       :: nRP(2) ! RP resolution in i,j direction i: P1->P2, 
  INTEGER,ALLOCATABLE           :: IDlist(:,:)
  REAL,ALLOCATABLE              :: NormVec(:,:),TangVec(:,:),LocalCoord(:,:,:)
  REAL,ALLOCATABLE              :: BLProps(:,:)
END TYPE tPlane

TYPE(tLine),POINTER             :: Lines(:)
TYPE(tPlane),POINTER            :: Planes(:)

!===================================================================================================================================
END MODULE MOD_RPSetVisuVisu_Vars
