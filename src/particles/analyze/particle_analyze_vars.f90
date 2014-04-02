MODULE MOD_Particle_Analyze_Vars
!===================================================================================================================================
! Contains global variables used by the Analyze modules.
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
LOGICAL            :: ParticleAnalyzeInitIsDone = .FALSE.
LOGICAL            :: CalcNumSpec
LOGICAL            :: CalcCharge
LOGICAL            :: CalcEpot
LOGICAL            :: CalcEkin
LOGICAL            :: CalcTemp
LOGICAL            :: TrackParticlePosition
INTEGER            :: nEkin
LOGICAL            :: DoAnalyze
LOGICAL            :: IsRestart
LOGICAL            :: ChargeCalcDone
LOGICAL            :: CalcShapeEfficiency
CHARACTER(LEN=256) :: CalcShapeEfficiencyMethod      ! Explanations in particle_analyze.f90
INTEGER            :: ShapeEfficiencyNumber          ! Explanations in particle_analyze.f90
INTEGER            :: PartAnalyzeStep
END MODULE MOD_Particle_Analyze_Vars
