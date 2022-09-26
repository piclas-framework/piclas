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
MODULE MOD_DG_Vars
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
! DG basis
REAL,ALLOCATABLE                      :: D(:,:)
REAL,ALLOCATABLE                      :: D_T(:,:)    ! D^T
REAL,ALLOCATABLE                      :: D_Hat(:,:)
REAL,ALLOCATABLE                      :: D_Hat_T(:,:) ! D_Hat^T
REAL,ALLOCATABLE                      :: L_HatMinus(:)
REAL,ALLOCATABLE                      :: L_HatPlus(:)
! DG solution (reference / physical)
REAL,ALLOCATABLE,TARGET               :: U(:,:,:,:,:) ! computed from JU
! DG time derivative
REAL,ALLOCATABLE                      :: Ut(:,:,:,:,:)
! number of array items in U, Ut, gradUx, gradUy, gradUz after allocated
INTEGER                               :: nTotalU
INTEGER                               :: nTotal_vol    !loop i,j,k
INTEGER                               :: nTotal_face   !loop i,j
! interior face values for all elements
REAL,ALLOCATABLE                      :: U_master(:,:,:,:),U_slave(:,:,:,:)
REAL,ALLOCATABLE                      :: Flux_Master(:,:,:,:)
REAL,ALLOCATABLE                      :: Flux_Slave(:,:,:,:)
! face values for Riemann
REAL,ALLOCATABLE                      :: U_Master_loc(:,:,:)
REAL,ALLOCATABLE                      :: U_Slave_loc (:,:,:)
REAL,ALLOCATABLE                      :: Flux_loc(:,:,:)

LOGICAL                               :: DGInitIsDone=.FALSE.
#if defined(IMPA) || defined(ROS)
REAL,ALLOCATABLE                      :: Un(:,:,:,:,:) ! computed from JU
#endif
!===================================================================================================================================
END MODULE MOD_DG_Vars
