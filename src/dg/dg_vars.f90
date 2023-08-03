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
! number of array items in U, Ut, gradUx, gradUy, gradUz after allocated
!INTEGER                               :: nTotal_vol    !loop i,j,k
!INTEGER                               :: nTotal_face   !loop i,j

#if defined(IMPA) || defined(ROS)
REAL,ALLOCATABLE                      :: Un(:,:,:,:,:) ! computed from JU
#endif

! Element local polynomial degrees
INTEGER,ALLOCATABLE                   :: N_DG(:)                !< polynomial degree inside DG element,         size(nElems)
INTEGER,ALLOCATABLE                   :: DG_Elems_master(:)     !< prolongate local polynomial degree to faces, size(nSides)
INTEGER,ALLOCATABLE                   :: DG_Elems_slave(:)      !< prolongate local polynomial degree to faces, size(nSides)

! DG basis, contains the differentiation and interpolation operators
TYPE, PUBLIC :: DG_Basis
  REAL,ALLOCATABLE                      :: D(:,:)               !< Differentiation matrix of size [0..N,0..N], contains the
                                                                !< first derivative of each Lagrange polynomial at each node.

  REAL,ALLOCATABLE                      :: D_T(:,:)             !< Transpose of differentiation matrix, size [0..N,0..N].

  REAL,ALLOCATABLE                      :: D_Hat(:,:)           !< Differentiation matrix premultiplied by
                                                                !< mass matrix, \f$ \hat{D} = M^{-1} D^T M \f$, size [0..N,0..N].

  REAL,ALLOCATABLE                      :: D_Hat_T(:,:)         !< Transpose of differentiation matrix premultiplied by
                                                                !< mass matrix, size [0..N,0..N].

  REAL,ALLOCATABLE                      :: L_HatMinus(:)        !< Lagrange polynomials evaluated at \f$\xi=-1\f$
                                                                !< premultiplied by mass matrix

  REAL,ALLOCATABLE                      :: L_HatPlus(:)         !< Lagrange polynomials evaluated at \f$\xi=+1\f$
                                                                !< premultiplied by mass matrix

  INTEGER                               :: nDOFElem
END TYPE DG_Basis

TYPE(DG_Basis),ALLOCATABLE              :: DGB_N(:)             !< Array of precomputed DG basis tensors for variable N_DG

! DG solution vol
TYPE N_U_Vol
  REAL,ALLOCATABLE  :: U(:,:,:,:)
  REAL,ALLOCATABLE  :: Ut(:,:,:,:)
  REAL,ALLOCATABLE  :: E(:,:,:,:)
  REAL,ALLOCATABLE  :: Et(:,:,:,:)
END TYPE N_U_Vol

! DG solution (JU or U) vectors)
TYPE(N_U_Vol),ALLOCATABLE :: U_N(:)       !< Solution variable for each equation, node and element,

! DG solution face
TYPE N_U_Surf
  REAL,ALLOCATABLE  :: U_master(:,:,:)
  REAL,ALLOCATABLE  :: U_slave(:,:,:)
  REAL,ALLOCATABLE  :: Flux_Master(:,:,:)
  REAL,ALLOCATABLE  :: Flux_Slave(:,:,:)
END TYPE N_U_Surf

! DG solution (JU or U) vectors)
TYPE(N_U_Surf),ALLOCATABLE :: U_Surf_N(:) !< Solution variable for each equation, node and element,

!----------------------------------------------------------------------------------------------------------------------------------
! Auxilliary variables
LOGICAL             :: DGInitIsDone=.FALSE.
!===================================================================================================================================
END MODULE MOD_DG_Vars
