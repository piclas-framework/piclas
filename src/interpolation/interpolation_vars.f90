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
REAL,ALLOCATABLE  :: L_Plus(:)                   !< L for boundary flux computation at plus side  (1)
REAL,ALLOCATABLE  :: L_Minus(:)                  !< L for boundary flux computation at minus side (-1)
REAL,ALLOCATABLE  :: L_PlusMinus(:,:)            !< L for boundary flux computation at both sides (-1,1)
REAL,ALLOCATABLE  :: xGP(:)                      !< Gauss point coordinates
REAL,ALLOCATABLE  :: wGP(:)                      !< GP integration weights
REAL,ALLOCATABLE  :: swGP(:)                     !< 1.0/ GP integration weights
REAL,ALLOCATABLE  :: wBary(:)                    !< barycentric weights
REAL,ALLOCATABLE  :: wGPSurf(:,:)                !< wGPSurf(i,j)=wGP(i)*wGP(j)
REAL,ALLOCATABLE  :: NChooseK(:,:)               !< array n over n
REAL,ALLOCATABLE  :: Vdm_Leg(:,:), sVdm_Leg(:,:) !< Legendre Vandermonde matrix
! Analyze variables 
INTEGER           :: NAnalyze                    !< number of analyzation points is NAnalyze+1
REAL,ALLOCATABLE  :: wAnalyze(:)                 !< GL integration weights used for the analyze
REAL,ALLOCATABLE  :: Vdm_GaussN_NAnalyze(:,:)    !< for interpolation to Analyze points (from NodeType nodes to Gauss-Lobatto nodes)
REAL,ALLOCATABLE  :: Uex(:,:,:,:,:)              !< Analytic solution Uex(1:PP_nVar,0:NAnalyze,0:NAnalyze,0:NAnalyze,1:nElems)

!==================================================================================================================================
!@{ Named nodetype parameters
!==================================================================================================================================
CHARACTER(LEN=255),PARAMETER :: NodeTypeG    = 'GAUSS'                    !< Gauss nodes (-1,1)
CHARACTER(LEN=255),PARAMETER :: NodeTypeGL   = 'GAUSS-LOBATTO'            !< Gauss-Lobatto nodes [-1,1]
CHARACTER(LEN=255),PARAMETER :: NodeTypeCL   = 'CHEBYSHEV-GAUSS-LOBATTO'  !< Chebyshev-Gauss-Lobatto nodes [-1,1]
CHARACTER(LEN=255),PARAMETER :: NodeTypeVISU = 'VISU'                     !< equidistant nodes [-1,1]
#if (PP_NodeType==1)
  CHARACTER(LEN=255),PARAMETER :: NodeType = 'GAUSS'
#elif (PP_NodeType==2)
  CHARACTER(LEN=255),PARAMETER :: NodeType = 'GAUSS-LOBATTO'
#elif (PP_NodeType==3)
  CHARACTER(LEN=255),PARAMETER :: NodeType = 'CHEBYSHEV-GAUSS-LOBATTO'
#endif

!==================================================================================================================================
! (Magnetic) background field
!==================================================================================================================================
INTEGER          :: NBG                        !< Polynomial degree of BG-Field
INTEGER          :: BGType                     !< Type of BG-Field (Electric,Magnetic,Both)
INTEGER          :: BGDataSize                 !< Type of BG-Field (Electric,Magnetic,Both)
REAL,ALLOCATABLE :: BGField(:,:,:,:,:)         !< BGField numerical solution (1:x,0:NBG,0:NBG,0:NBG,1:PP_nElems)
REAL,ALLOCATABLE :: BGFieldAnalytic(:,:,:,:,:) !< BGField analytic solution (1:x,0:NBG,0:NBG,0:NBG,1:PP_nElems)
REAL,ALLOCATABLE :: BGField_xGP(:)             !< Gauss point coordinates
REAL,ALLOCATABLE :: BGField_wBary(:)           !< barycentric weights
LOGICAL          :: BGFieldVTKOutput           !< Output the background field in VTK data format
REAL,ALLOCATABLE :: PsiMag(:,:,:,:)            !< Magnetic potential [1:PP_N,1:PP_N,1:PP_N,1:nElems]
!===================================================================================================================================

LOGICAL           :: InterpolationInitIsDone = .FALSE. !< Flag whether the initialization has been completed or not
END MODULE MOD_Interpolation_Vars
