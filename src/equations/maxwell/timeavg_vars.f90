MODULE MOD_Timeaverage_Vars
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
!===================================================================================================================================
INTEGER                        :: nMaxVarAvg              !< max number of variables for averaging
INTEGER                        :: nMaxVarFluc             !< max number of variables for RMS
! Time averaging and fluctuation variables
LOGICAL              :: doCalcTimeAverage   =.FALSE.      !< marks if time averaging should be performed
LOGICAL              :: doCalcFluctuations  =.FALSE.      !< marks if time fluctuations should be computed
REAL   ,ALLOCATABLE  :: UAvg(:,:,:,:,:)                   !< time averaged solution U
REAL   ,ALLOCATABLE  :: UFluc(:,:,:,:,:)                  !< time averaged solution squared (U^2)
LOGICAL,ALLOCATABLE  :: CalcAvg(:)                        !< variables for which time averages should be computed (global indexing)
LOGICAL,ALLOCATABLE  :: CalcFluc(:)                       !< variables for which fluctuations should be computed (global indexing)
INTEGER,ALLOCATABLE  :: iAvg(:)                           !< map from (global) VariableList to index in UAvg array
INTEGER,ALLOCATABLE  :: iFluc(:)                          !< map from (global) VariableList to index in UFluc array
INTEGER,ALLOCATABLE  :: FlucAvgMap(:,:)                   !< map from index in UFluc array to index in UAvg array
                                                          !< (e.g. for mixed term uv: iFluc(1,1) -> u iFluc(2,1) -> v)
INTEGER              :: nVarAvg                           !< number of time averag variables
INTEGER              :: nVarFluc                          !< number of fluctuation variables
INTEGER              :: nVarFlucHasAvg                    !< number of fluctuations depending only on one time average
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesAvgOut(:)       !< time averaged variable names
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesFlucOut(:)      !< fluctuation variable names
REAL                 :: dtAvg                             !< sum of timesteps
REAL                 :: dtOld                             !< dt from previous iteration
LOGICAL              :: DoPoyntingVectorAvg               !< logical if PoyntingVector is sampled
#ifdef PARTICLES
LOGICAL,ALLOCATABLE  :: DoPowerDensity(:)                 !> Sample Power-Density of species
REAL,ALLOCATABLE     :: PowerDensity(:,:,:,:,:,:)         !> Power-Density of species
INTEGER              :: nSpecPowerDensity
#endif /*PARTICLES*/

END MODULE MOD_Timeaverage_Vars
