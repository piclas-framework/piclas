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
MODULE MOD_Interpolation_Vars
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
! reserved for Gauss Points with polynomial degree N, all allocated (0:N)
REAL,ALLOCATABLE  :: L_Plus(:), L_Minus(:)       ! L for boundary flux computation at both sides (-1,1)
REAL,ALLOCATABLE  :: L_PlusMinus(:,:)            ! L for boundary flux computation at both sides (-1,1)
REAL,ALLOCATABLE  :: xGP(:)                      ! Gauss point coordinates
REAL,ALLOCATABLE  :: wGP(:)                      ! GP integration weights
REAL,ALLOCATABLE  :: swGP(:)                     ! 1.0/ GP integration weights
REAL,ALLOCATABLE  :: wBary(:)                    ! barycentric weights
REAL,ALLOCATABLE  :: wGPSurf(:,:)                ! wGPSurf(i,j)=wGP(i)*wGP(j)
REAL,ALLOCATABLE  :: NChooseK(:,:)               ! array n over n
REAL,ALLOCATABLE  :: Vdm_Leg(:,:), sVdm_Leg(:,:) !< Legendre Vandermonde matrix
CHARACTER(LEN=255),PARAMETER :: NodeTypeG    = 'GAUSS'                    !< Gauss nodes (-1,1)
CHARACTER(LEN=255),PARAMETER :: NodeTypeGL   = 'GAUSS-LOBATTO'            !< Gauss-Lobatto nodes [-1,1]
CHARACTER(LEN=255),PARAMETER :: NodeTypeCL   = 'CHEBYSHEV-GAUSS-LOBATTO'
CHARACTER(LEN=255),PARAMETER :: NodeTypeVISU = 'VISU'                     !< equidistant nodes [-1,1]
#if (PP_NodeType==1)
  CHARACTER(LEN=255),PARAMETER :: NodeType = 'GAUSS'
#elif (PP_NodeType==2)
  CHARACTER(LEN=255),PARAMETER :: NodeType = 'GAUSS-LOBATTO'
#elif (PP_NodeType==3)
  CHARACTER(LEN=255),PARAMETER :: NodeType = 'CHEBYSHEV-GAUSS-LOBATTO'
#endif
!===================================================================================================================================

LOGICAL           :: InterpolationInitIsDone = .FALSE.
END MODULE MOD_Interpolation_Vars
