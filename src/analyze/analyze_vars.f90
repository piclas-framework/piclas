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
!> Contains global variables used by the Analyze modules.
!===================================================================================================================================
MODULE MOD_Analyze_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER           :: NAnalyze                    !< number of analyzation points is NAnalyze+1
REAL,ALLOCATABLE  :: wAnalyze(:)                 !< GL integration weights used for the analyze
REAL,ALLOCATABLE  :: Vdm_GaussN_NAnalyze(:,:)    !< for interpolation to Analyze points
REAL              :: Analyze_dt                  !< time difference to trigger analyze output
INTEGER(KIND=8)   :: iAnalyze                    !> count number of next analyze
REAL              :: OutputTimeFixed             !< fixed time for writing state to .h5
LOGICAL           :: CalcPoyntingInt             !< calculate pointing vector integral | only perp to z axis
LOGICAL           :: CalcMeshInfo                !< Output myrank, ElemID and tracking info to ElemData
LOGICAL           :: CalcHaloInfo                !< Output halo element information to ElemData
REAL              :: PoyntingIntCoordErr         !< tolerance in plane searching
REAL              :: PoyntingIntPlaneFactor      !< factor for poyntingintplanes
INTEGER           :: nPoyntingIntPlanes          !< number of planes
REAL,ALLOCATABLE  :: PoyntingIntegral(:)         !< poyntingintegral of each plane
INTEGER           :: AnalyzeCount                !< number of analyzes (for info)
REAL              :: AnalyzeTime                 !< accumulated time of analyzes (for info)
REAL,ALLOCATABLE  :: PosPoyntingInt(:)           !< z-coordinate of plane
REAL,ALLOCATABLE  :: S(:,:,:,:), STEM(:,:,:)     !< vector, abs for TEM waves
LOGICAL           :: DoFieldAnalyze              !< perform analyze
INTEGER           :: FieldAnalyzeStep            !< Analyze is performed each Nth time step
LOGICAL           :: DoCalcErrorNorms            !< perform L2, LInf error calculation
LOGICAL           :: DoSurfModelAnalyze          !< perform analyze for SurfaceModel
LOGICAL           :: CalcEpot                    !< Computation of the energy stored in the electric and
REAL              :: Wel                         !< energy of the electric field
REAL              :: Wmag                        !< energy of the magnetic field
REAL              :: Wphi
REAL              :: Wpsi
! magnetic field
LOGICAL           :: OutputErrorNorms            !< print L2 norms (DG state and particles if present)
#ifdef CODE_ANALYZE
LOGICAL           :: DoCodeAnalyzeOutput         !< print code analyze info to CodeAnalyze.csv (default is TRUE)
#endif /* CODE_ANALYZE */
LOGICAL           :: CalcPointsPerWavelength     !< Flag to compute the points per wavelength in each cell (assume equidistant DOF
!                                                !< distribution within each cell
!                                                !< PPW = (p+1)*lambda / GEO%CharLength
!                                                !<   GEO%CharLength = (V_cell)^(1/3)          characteristic length in the cell
REAL,ALLOCATABLE  :: PPWCell(:)                  !< Points per wavelength for each cell
!===================================================================================================================================
LOGICAL           :: AnalyzeInitIsDone = .FALSE.
END MODULE MOD_Analyze_Vars
