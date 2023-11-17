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
!===================================================================================================================================
!> Contains global variables used by the Analyze modules.
!===================================================================================================================================
MODULE MOD_Analyze_Vars
! MODULES
#if USE_MPI
USE MOD_Globals
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL              :: Analyze_dt                  !< time difference to trigger analyze output
INTEGER(KIND=8)   :: nSkipAnalyze                !< Skip Analyze_dt
REAL              :: SkipAnalyzeWindow           !< Reoccurring time frame when switching between nSkipAnalyze and SkipAnalyzeWindow
REAL              :: SkipAnalyzeSwitchTime       !< Time within the reoccurring time frame, when using nSkipAnalyzeSwitch instead of
                                                 !< nSkipAnalyze
INTEGER(KIND=8)   :: nSkipAnalyzeSwitch          !< Skip Analyze_dt with a different values as nSkipAnalyze
INTEGER(KIND=8)   :: iAnalyze                    !< count number of next analyze
REAL              :: OutputTimeFixed             !< fixed time for writing state to .h5
LOGICAL           :: CalcMeshInfo                !< Output myrank, ElemID and tracking info to ElemData
LOGICAL           :: CalcHaloInfo                !< Output halo element information to ElemData
INTEGER           :: AnalyzeCount                !< number of analyzes (for info)
REAL              :: AnalyzeTime                 !< Accumulated wall time of analyzes (for info)
REAL,ALLOCATABLE  :: S(:,:,:,:), STEM(:,:,:)     !< vector, abs for TEM waves
LOGICAL           :: DoFieldAnalyze              !< perform analyze
LOGICAL           :: DoMeasureAnalyzeTime        !< measure time that is spent in analyze routines and count the number of analysis
                                                 !< calls (to std out stream)
INTEGER(KIND=8)   :: FieldAnalyzeStep            !< Analyze is performed each Nth time step
LOGICAL           :: DoCalcErrorNorms            !< perform L2, LInf error calculation
LOGICAL           :: OutputErrorNormsToH5        !< Set true to write the analytical solution, the L2 and LInf error norms at analyze step to .h5 state file. Default = F
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
#if (PP_nVar>=6)
!-----------------------------------------------------------------------------------------------------------------------------------
!< PoyntingVectorIntegral variables
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL             :: CalcPoyntingInt        !< calculate pointing vector integral | only perp to z axis
INTEGER             :: PoyntingMainDir        !< direction in which the Poynting vector integral is to be computed
LOGICAL,ALLOCATABLE :: isPoyntingIntSide(:)   !< number of all PoyntingInt sides
INTEGER,ALLOCATABLE :: SideIDToPoyntingSide(:)!< plane number used for calculation of Poynting vector
REAL,ALLOCATABLE    :: PosPoyntingInt(:)      !< x- y- or z-coordinate of plane
REAL                :: PoyntingIntCoordErr    !< tolerance in plane searching
INTEGER             :: nPoyntingIntPlanes     !< number of planes
REAL,ALLOCATABLE    :: PoyntingIntegral(:)    !< poyntingintegral of each plane
#endif
#if USE_HDG
LOGICAL             :: CalcAverageElectricPotential!< flag for activating the usage
REAL                :: AverageElectricPotential    !< 2D Landmark: averaged electric field in y-direction (for interpolation)
LOGICAL,ALLOCATABLE :: isAverageElecPotSide(:)     !< Flag for sides that are required for calculating the averaged electric potential
REAL                :: AverageElectricPotentialCoordErr !< tolerance in plane searching
REAL                :: PosAverageElectricPotential      !< x-coordinate of plane
INTEGER             :: AverageElectricPotentialFaces    !< global number of faces
LOGICAL             :: CalcElectricTimeDerivative       !< Calculate the time derivative of E and output to h5
#endif /*USE_HDG*/
!===================================================================================================================================
! --- BoundaryFieldOutput = BFO
LOGICAL                       :: CalcBoundaryFieldOutput !< Flag for activating this output

TYPE tBoundaryFieldOutput
  INTEGER                       :: NFieldBoundaries   !< Total number of boundaries where the field BC is stored to .csv
  INTEGER,ALLOCATABLE           :: FieldBoundaries(:) !< Field-boundary number (iBC)
END TYPE

TYPE(tBoundaryFieldOutput)   :: BFO
!===================================================================================================================================
!-- Electric displacement current

#if USE_MPI
TYPE tMPIGROUP
  INTEGER                     :: ID                     !< MPI communicator ID
  INTEGER                     :: UNICATOR=MPI_COMM_NULL !< MPI communicator for electric displacement current
  INTEGER                     :: nProcs                 !< number of MPI processes for particles
  INTEGER                     :: MyRank                 !< MyRank within communicator
END TYPE
#endif /*USE_MPI*/

TYPE tEDC
  REAL,ALLOCATABLE            :: Current(:)          !< Electric displacement current for each (required) BC index
#if USE_MPI
  TYPE(tMPIGROUP),ALLOCATABLE :: COMM(:)             !< communicator and ID for parallel execution
#endif /*USE_MPI*/
  INTEGER                     :: NBoundaries         !< Total number of boundaries where the electric displacement current is evaluated
  INTEGER,ALLOCATABLE         :: FieldBoundaries(:)  !< Field-boundary number on which the particles are counted
  INTEGER,ALLOCATABLE         :: BCIDToEDCBCID(:)    !< Mapping BCID to EDC BCID (1:nPartBound)
END TYPE

TYPE(tEDC)   :: EDC
!===================================================================================================================================
LOGICAL           :: AnalyzeInitIsDone = .FALSE.
END MODULE MOD_Analyze_Vars
