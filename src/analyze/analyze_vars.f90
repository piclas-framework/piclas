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
REAL              :: OutputTimeFixed             !< fixed time for writing state to .h5
LOGICAL           :: CalcPoyntingInt             !< calulate pointing vector integral | only perp to z axis
REAL              :: PoyntingIntCoordErr         !< tolerance in plane searching
INTEGER           :: nPoyntingIntPlanes          !< number of planes
REAL,ALLOCATABLE  :: PosPoyntingInt(:)           !< z-coordinate of plane
REAL,ALLOCATABLE  :: PoyntingIntPlaneFactor(:)   !< plane factor
REAL,ALLOCATABLE  :: S(:,:,:,:), STEM(:,:,:)     !< vector, abs for TEM waves
LOGICAL           :: DoAnalyze                   !< perform analyze
LOGICAL           :: DoSurfModelAnalyze          !< perform analyze for SurfaceModel
LOGICAL           :: DoCalcErrorNorms            !< perform L2, LInf error calculation
LOGICAL           :: CalcEpot                    !< Computation of the energy stored in the electric and
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
