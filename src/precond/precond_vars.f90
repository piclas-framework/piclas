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
#include "piclas.h"

MODULE MOD_Precond_Vars
!===================================================================================================================================
! Contains global variables used by the Timedisc modules.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE      :: invP(:,:,:)  !inverse of block Jacobian for each element (1:nDOF_elem,1:nDOFelem,1:PP_nElems)
INTEGER               :: PrecondType
INTEGER               :: DebugMatrix
LOGICAL               :: doVol,doSurf
LOGICAL               :: UpdatePrecond
#ifdef maxwell
#if USE_LOADBALANCE
LOGICAL               :: UpdatePrecondLB
#endif /*USE_LOADBALANCE*/
#endif /*maxwell*/
LOGICAL               :: PrecondInitIsDone
REAL,ALLOCATABLE      :: invBJ(:,:,:,:,:,:)  !inverse of block Jacobian for each DOF (1:PP_nVar,1:PP_nVar,0:N,0:N,0:N,1:PP_nElems)
REAL,ALLOCATABLE      :: invJ(:,:,:,:,:)  ! inverse of Jacobian (1:PP_nVar,0:N,0:N,0:N,1:PP_nElems)
INTEGER               :: PrecondMethod
INTEGER,ALLOCATABLE   :: NeighborElemID(:,:)  ! 1:6, 1:PP_nElems
INTEGER,ALLOCATABLE   :: MatPatternElemID(:,:) ! number of row, ElemID, row equals diag-elemid
INTEGER,ALLOCATABLE   :: MatPattern(:,:)       ! row, locSideID
INTEGER,ALLOCATABLE   :: Lower(:,:)            ! index range for MatPattern, containing the lower triangular matrix
INTEGER,ALLOCATABLE   :: Upper(:,:)            ! index range for MatPattern, containing the upper triangular matrix
REAL,ALLOCATABLE      :: ProcJacobian(:,:,:,:)
INTEGER               :: BlockSize
INTEGER               :: nBlockSize
REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:) :: invXi,invEta,invZeta
REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:) :: dRdXi,dRdEta,dRdZeta
LOGICAL               :: BuildNvecisDone=.FALSE.
!===================================================================================================================================
END MODULE MOD_Precond_Vars
