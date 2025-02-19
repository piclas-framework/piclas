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

MODULE MOD_Mesh_Vars_FV

!===================================================================================================================================
!> Contains finite volumes metrics
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE

#if USE_FV
REAL,ALLOCATABLE :: NormVec_FV(:,:,:,:)           !< normal vector for each side       (1:3,0:N,0:N,nSides)
REAL,ALLOCATABLE :: TangVec1_FV(:,:,:,:)          !< tangential vector 1 for each side (1:3,0:N,0:N,nSides)
REAL,ALLOCATABLE :: TangVec2_FV(:,:,:,:)          !< tangential vector 3 for each side (1:3,0:N,0:N,nSides)
REAL,ALLOCATABLE :: SurfElem_FV(:,:,:)            !< surface area for each side        (    0:N,0:N,nSides)
REAL,ALLOCATABLE :: Ja_Face_FV(:,:,:,:,:)         !< surface  metrics for each side
REAL,ALLOCATABLE :: Metrics_fTilde_FV(:,:,:,:,:) !< Metric Terms (first indices 3) on each GaussPoint
REAL,ALLOCATABLE :: Metrics_gTilde_FV(:,:,:,:,:)
REAL,ALLOCATABLE :: Metrics_hTilde_FV(:,:,:,:,:)
REAL,ALLOCATABLE,TARGET :: Elem_xGP_FV(:,:,:,:,:)          !< XYZ positions (first index 1:3) of the volume Gauss Point
REAL,ALLOCATABLE :: Face_xGP_FV(:,:,:,:)            !< XYZ positions (first index 1:3) of the Boundary Face Gauss Point
REAL,ALLOCATABLE :: Face_xGP_PP_1(:,:,:,:)            !< XYZ positions (first index 1:3) of the Boundary Face Gauss Point
REAL,ALLOCATABLE :: NormVec_PP_1(:,:,:,:)           !< normal vector for each side       (1:3,0:N,0:N,nSides)
REAL,ALLOCATABLE :: TangVec1_PP_1(:,:,:,:)          !< tangential vector 1 for each side (1:3,0:N,0:N,nSides)
REAL,ALLOCATABLE :: TangVec2_PP_1(:,:,:,:)          !< tangential vector 3 for each side (1:3,0:N,0:N,nSides)
REAL,ALLOCATABLE :: SurfElem_PP_1(:,:,:)            !< surface area for each side        (    0:N,0:N,nSides)
REAL,ALLOCATABLE :: Ja_Face_PP_1(:,:,:,:,:)         !< surface  metrics for each side
REAL,ALLOCATABLE :: Vdm_CLN_N_PP_1(:,:)
REAL,ALLOCATABLE :: XCL_N_PP_1(:,:,:,:,:)             !< mapping X(xi) P\in N
REAL,ALLOCATABLE :: dXCL_N_PP_1(:,:,:,:,:,:)    !< jacobi matrix of the mapping P\in NGeo
REAL,ALLOCATABLE :: JaCL_N_PP_1(:,:,:,:,:,:)    !< metric terms P\in N
LOGICAL,ALLOCATABLE :: IsPeriodicSide(:)        !< Periodic flag for gradient calculations

#if USE_HDG
INTEGER,ALLOCATABLE :: BoundaryType_FV(:,:) !< BCType    = BoundaryType(BC(SideID),BC_TYPE)
                                            !< BCState   = BoundaryType(BC(SideID),BC_STATE)
#endif

! For surface output
INTEGER          :: nGlobalBCSides=0
INTEGER          :: offsetBCSide=0
INTEGER                                 :: nFVSurfBC                       ! number of surface side BCs
CHARACTER(LEN=255),ALLOCATABLE          :: FVSurfBCName(:)                 ! names of belonging surface BC

#endif /*USE_FV*/
END MODULE MOD_Mesh_Vars_FV