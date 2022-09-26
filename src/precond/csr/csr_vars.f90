!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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
MODULE MOD_CSR_Vars
!===================================================================================================================================
! Contains global variables used by the DG modules.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! sparse test for calculation
INTEGER,ALLOCATABLE,DIMENSION(:,:)    :: nNonZeros
!REAL,ALLOCATABLE,DIMENSION(:,:)       :: Ucsr,Utcsr
TYPE tSparse
 REAL,ALLOCATABLE,DIMENSION(:)        :: Entry
 INTEGER,ALLOCATABLE,DIMENSION(:)     :: IEntry,JEntry
END TYPE
REAL,ALLOCATABLE                      :: DebugMat(:,:)
TYPE(tSparse), ALLOCATABLE            :: SparseMatrix(:,:)
REAL,ALLOCATABLE                      :: L_HatPlusMinus(:,:)
REAL,ALLOCATABLE                      :: GlobalAA(:)
REAL,ALLOCATABLE                      :: GlobalBAA(:,:,:)
INTEGER,ALLOCATABLE                   :: GlobalIA(:),GlobalJA(:), GlobalDiag(:)
INTEGER                               :: nonZerosGlobal
REAL                                  :: epsZero
INTEGER,ALLOCATABLE,DIMENSION(:)      :: nUNonZeros,nLNonZeros
INTEGER                               :: nMTriangle
REAL,ALLOCATABLE,DIMENSION(:,:)       :: DE
TYPE tIL                                                                     ! ILU for each element
 REAL,ALLOCATABLE,DIMENSION(:)        :: Entry
 INTEGER,ALLOCATABLE,DIMENSION(:)     :: IEntry,JEntry
END TYPE
TYPE(tIL), ALLOCATABLE                :: IL(:)
TYPE(tIL), ALLOCATABLE                :: IU(:)
!TYPE tAA                                                                     ! ILU for each element
! REAL,ALLOCATABLE,DIMENSION(:)        :: AA
! INTEGER,ALLOCATABLE,DIMENSION(:)     :: IA,JA,Diag
!END TYPE
!TYPE(taa),ALLOCATABLE,DIMENSION(:)    :: EILU
!===================================================================================================================================
END MODULE MOD_CSR_Vars
