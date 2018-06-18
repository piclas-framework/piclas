MODULE MOD_ESBGK_Vars
!===================================================================================================================================
! Contains the FP Flow variables
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                                        :: ESBGK_Flow

TYPE tSpeciesESBGK                                                              ! ESBK Species Param
  REAL, ALLOCATABLE                            :: CollFreqPreFactor(:)
END TYPE tSpeciesESBGK

TYPE(tSpeciesESBGK), ALLOCATABLE               :: SpecESBGK(:)                  ! Species DSMC params (nSpec)
LOGICAL                                        :: DoBGKCellAdaptation
REAL                                           :: ESBGKTempCorrectFact

INTEGER                                        :: BGKCollModel                  ! 1 ES-BGK; 2 S-BGK; 3 BGK
INTEGER                                        :: ESBGKModel                    ! 1 Approx Levin; 2 Exact Solution A; 3 Metropolis
REAL                                           :: BGKUnifiedCes
INTEGER                                        :: BGKMinPartPerCell
INTEGER                                        :: BGKAveragingLength
LOGICAL                                        :: BGKDoAveraging
LOGICAL                                        :: BGKDoAveragingCorrect
LOGICAL                                        :: BGKUseQuantVibEn
REAL                                           :: BGKMinCFL = 1.0
INTEGER                                        :: BGKAdaptTimeStep
REAL                                           :: BGKAcceleration
LOGICAL                         	       :: BGKDoVibRelaxation

LOGICAL                                        :: DoBGKCellSplitting
LOGICAL                                        :: BGKSampAdapFac
TYPE tElemSplitCells
  REAL, ALLOCATABLE                            :: AdapFac(:)
  REAL, ALLOCATABLE                            :: SplitCellVolumes(:,:,:)
  INTEGER                                      :: Splitnum(3)
  INTEGER                                      :: CellOrientation(3)
END TYPE

TYPE(tElemSplitCells), ALLOCATABLE             :: ElemSplitCells(:)

REAL                                           :: BGKDiffEn = 0.0
REAL                                           :: BGKDiffEn2 = 0.0
REAL                                           :: BGKDiffEn3 = 0.0
REAL                                           :: BGKDiffEn4 = 0.0
INTEGER                                        :: BGKTest = 0

TYPE tElemNodeAveraging
    TYPE (tNodeAverage), POINTER               :: Root => null()
END TYPE

TYPE tNodeAverage
    TYPE (tNodeAverage), POINTER               :: SubNode1 => null()
    TYPE (tNodeAverage), POINTER               :: SubNode2 => null()
    TYPE (tNodeAverage), POINTER               :: SubNode3 => null()
    TYPE (tNodeAverage), POINTER               :: SubNode4 => null()
    TYPE (tNodeAverage), POINTER               :: SubNode5 => null()
    TYPE (tNodeAverage), POINTER               :: SubNode6 => null()
    TYPE (tNodeAverage), POINTER               :: SubNode7 => null()
    TYPE (tNodeAverage), POINTER               :: SubNode8 => null()
    REAL, ALLOCATABLE                          :: AverageValues(:,:)
    INTEGER                                    :: CorrectStep
END TYPE

TYPE (tElemNodeAveraging), ALLOCATABLE         :: ElemNodeAveraging(:)
!===================================================================================================================================
END MODULE MOD_ESBGK_Vars
