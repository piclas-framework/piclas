!==================================================================================================================================
! Copyright (c) 2015 - 2019 Wladimir Reschke
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
#include "piclas.h"

MODULE MOD_SurfaceModel_Analyze
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
#ifdef PARTICLES
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitSurfModelAnalyze
  MODULE PROCEDURE InitSurfModelAnalyze
END INTERFACE

INTERFACE AnalyzeSurface
  MODULE PROCEDURE AnalyzeSurface
END INTERFACE

INTERFACE WriteDataHeaderInfo
  MODULE PROCEDURE WriteDataHeaderInfo
END INTERFACE

INTERFACE WriteDataInfo
  MODULE PROCEDURE WriteDataInfo
END INTERFACE

#if (PP_TimeDiscMethod==42)
INTERFACE AnalyzeSurfRates
  MODULE PROCEDURE AnalyzeSurfRates
END INTERFACE
PUBLIC:: AnalyzeSurfRates
#endif /*RESERVOIR*/

PUBLIC:: InitSurfModelAnalyze
PUBLIC:: AnalyzeSurface
PUBLIC:: DefineParametersSurfModelAnalyze
!===================================================================================================================================
CONTAINS

!==================================================================================================================================
!> Define parameters for analyze if wallmodel enabled (.csv output)
!==================================================================================================================================
SUBROUTINE DefineParametersSurfModelAnalyze()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
!USE MOD_AnalyzeEquation ,ONLY: DefineParametersAnalyzeEquation
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Surface Analyze")

CALL prms%CreateIntOption(      'Surface-AnalyzeStep'   , 'Analyze is performed each Nth time step for surfaces','1')
CALL prms%CreateLogicalOption(  'Surf-CalcCollCounter'  , 'Analyze the number of surface collision and number of '//&
                                                          'adsorbed particles per species','.FALSE.')
CALL prms%CreateLogicalOption(  'Surf-CalcDesCounter'   , 'Analyze the number of desorbed particles per species','.FALSE.')
CALL prms%CreateLogicalOption(  'Surf-CalcAdsProb'      , 'Analyze the mean probabilty for adsorption per species','.FALSE.')
CALL prms%CreateLogicalOption(  'Surf-CalcDesProb'      , 'Analyze the mean probablity for desorption per species','.FALSE.')
#if (PP_TimeDiscMethod==42) || (PP_TimeDiscMethod==4)
CALL prms%CreateLogicalOption(  'Surf-CalcNumSpec'      , 'Analyze the number of simulated particles per species on surfaces'&
                                                          ,'.FALSE.')
CALL prms%CreateLogicalOption(  'Surf-CalcCoverage'     , 'Analyze the mean surface coverages for each species','.FALSE.')
#if (PP_TimeDiscMethod==42)
CALL prms%CreateLogicalOption(  'Surf-CalcAccomodation' , 'Analyze the mean surface accomodation coefficient','.FALSE.')
CALL prms%CreateLogicalOption(  'Surf-CalcAdsorbRates'  , 'Analyze every refined rate data of gas-surface reactions.\n'//&
                                                          'Enables flags: CalcAdsorbProb / CalcAdsorbE / CalcAdsorbnu.','.FALSE.')
CALL prms%CreateLogicalOption(  'Surf-CalcAdsorbProb'   , 'Analyze eaction probabilities per reaction and species','.FALSE.')
CALL prms%CreateLogicalOption(  'Surf-CalcAdsorbE'      , 'Analyze activation barriers per reaction and species.','.FALSE.')
CALL prms%CreateLogicalOption(  'Surf-CalcAdsorbnu'     , 'Analyze reaction frequencies (nu_r) per reaction and species','.FALSE.')
CALL prms%CreateLogicalOption(  'Surf-CalcSurfRates'    , 'Analyze every refined rate data on the surfaces.\n'//&
                                                          'Enables flags: CalcSurfProb / CalcSurfE / CalcSurfnu.','.FALSE.')
CALL prms%CreateLogicalOption(  'Surf-CalcSurfProb'     , 'Analyze eaction probabilities per reaction and species','.FALSE.')
CALL prms%CreateLogicalOption(  'Surf-CalcSurfE'        , 'Analyze activation barriers per reaction and species.','.FALSE.')
CALL prms%CreateLogicalOption(  'Surf-CalcSurfnu'       , 'Analyze reaction frequencies (nu_r) per reaction and species','.FALSE.')
CALL prms%CreateLogicalOption(  'Surf-CalcHeatFlux'     , 'Analyze the the heat fluxes onto surface and corresponding reaction'//&
                                                          'counters per reaction','.FALSE.')
#endif
#endif

END SUBROUTINE DefineParametersSurfModelAnalyze


SUBROUTINE InitSurfModelAnalyze()
!===================================================================================================================================
!> Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_ReadInTools               ,ONLY: GETLOGICAL, GETINT
USE MOD_Analyze_Vars              ,ONLY: DoSurfModelAnalyze
USE MOD_SurfaceModel_Analyze_Vars
#if (PP_TimeDiscMethod==42)
USE MOD_SurfaceModel_Vars         ,ONLY: Adsorption
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF (SurfModelAnalyzeInitIsDone) THEN
CALL abort(__STAMP__,&
'InitParticleAnalyse already called.',999,999.)
  RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT SURFACE MODEL ANALYZE...'

SurfaceAnalyzeStep = GETINT('Surface-AnalyzeStep')
IF (SurfaceAnalyzeStep.EQ.0) SurfaceAnalyzeStep = HUGE(1)

DoSurfModelAnalyze = .FALSE.

CalcCollCounter = GETLOGICAL('Surf-CalcCollCounter')
CalcDesCounter = GETLOGICAL('Surf-CalcDesCounter')

#if (PP_TimeDiscMethod==42) || (PP_TimeDiscMethod==4)
CalcSurfNumSpec = GETLOGICAL('Surf-CalcNumSpec')
CalcSurfCoverage = GETLOGICAL('Surf-CalcCoverage')
#if (PP_TimeDiscMethod==42)
CalcAccomodation = GETLOGICAL('Surf-CalcAccomodation')
CalcAdsorbRates = GETLOGICAL('Surf-CalcAdsorbRates')
IF (CalcAdsorbRates) THEN
  CalcAdsorbProb  = .TRUE.
  CalcAdsorbE     = .TRUE.
  CalcAdsorbnu    = .TRUE.
ELSE
  CalcAdsorbProb  = GETLOGICAL('Surf-CalcAdsorbProb')
  CalcAdsorbE     = GETLOGICAL('Surf-CalcAdsorbE')
  CalcAdsorbnu    = GETLOGICAL('Surf-CalcAdsorbnu')
  IF (CalcAdsorbProb.OR.CalcAdsorbE.OR.CalcAdsorbnu) CalcAdsorbRates=.TRUE.
END IF
CalcSurfRates = GETLOGICAL('Surf-CalcSurfRates')
IF (CalcSurfRates) THEN
  CalcSurfProb  = .TRUE.
  CalcSurfnu    = .TRUE.
  CalcSurfE     = .TRUE.
ELSE
  CalcSurfProb  = GETLOGICAL('Surf-CalcSurfProb')
  CalcSurfnu    = GETLOGICAL('Surf-CalcSurfnu')
  CalcSurfE     = GETLOGICAL('Surf-CalcSurfE')
  IF (CalcSurfProb.OR.CalcSurfnu.OR.CalcSurfE) CalcSurfRates=.TRUE.
END IF
CalcHeatflux = GETLOGICAL('Surf-CalcHeatFlux')
IF (    CalcSurfNumSpec &
   .OR. CalcSurfRates &
   .OR. CalcSurfCoverage &
   .OR. CalcAccomodation &
   .OR. Adsorption%TPD &
   .OR. CalcAdsorbRates &
   .OR. CalcHeatFlux) &
  DoSurfModelAnalyze = .TRUE.
IF (Adsorption%TPD.AND.((.NOT.CalcSurfRates))) CalcSurfRates = .TRUE.
#else
IF(CalcSurfNumSpec.OR.CalcSurfCoverage.OR.CalcAccomodation) DoSurfModelAnalyze = .TRUE.
#endif
#endif
IF(CalcCollCounter.OR.CalcDesCounter.OR.CalcAdsProb.OR.CalcDesProb) DoSurfModelAnalyze = .TRUE.

SurfModelAnalyzeInitIsDone=.TRUE.

SWRITE(UNIT_stdOut,'(A)')' INIT SURFACE MODEL ANALYZE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitSurfModelAnalyze


SUBROUTINE AnalyzeSurface(Time)
!===================================================================================================================================
!> create/open SurfaceAnalyze.csv and write calculated variables for surface analyze
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars              ,ONLY: DoSurfModelAnalyze
USE MOD_SurfaceModel_Analyze_Vars
USE MOD_Restart_Vars              ,ONLY: DoRestart
USE MOD_Particle_Boundary_Vars    ,ONLY: SurfMesh
USE MOD_Particle_Boundary_Vars    ,ONLY: PartBound, nPartBound
#if USE_MPI
USE MOD_Particle_Boundary_Vars    ,ONLY: SurfCOMM
#endif /*USE_MPI*/
USE MOD_SurfaceModel_Vars         ,ONLY: SurfModel
#if ( PP_TimeDiscMethod ==42)
USE MOD_SurfaceModel_Vars         ,ONLY: Adsorption
#endif /* DSMC*/
USE MOD_Particle_Vars             ,ONLY: nSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: Time
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL             :: isOpen, isRestart, doDistributionData
CHARACTER(LEN=350)  :: outfile
INTEGER             :: unit_index, OutputCounter, iPB
INTEGER             :: SurfCollNum(nSpecies),AdsorptionNum(nSpecies),DesorptionNum(nSpecies)
REAL                :: MeanAdsorptionProb(nSpecies), MeanDesorptionProb(nSpecies)
INTEGER             :: iSpec
#if (PP_TimeDiscMethod ==42)
INTEGER             :: iCase
REAL                :: Accomodation(nSpecies)
REAL,ALLOCATABLE    :: SurfReactRate(:), AdsorptionReactRate(:)
REAL,ALLOCATABLE    :: AdsorptionActE(:), ProperAdsorptionActE(:), Adsorptionnu(:), ProperAdsorptionnu(:)
REAL,ALLOCATABLE    :: SurfaceActE(:), ProperSurfaceActE(:), Surfacenu(:), ProperSurfacenu(:)
REAL,ALLOCATABLE    :: HeatFlux(:,:), AdsReactCount(:), DesReactCount(:)
#endif
#if (PP_TimeDiscMethod ==42) || (PP_TimeDiscMethod ==4)
INTEGER(KIND=8)     :: WallNumSpec(nSpecies), WallNumSpec_SurfDist(nSpecies)
REAL                :: WallCoverage(nSpecies)
#endif
!===================================================================================================================================
IF (.NOT.SurfMesh%SurfOnProc) RETURN
IF (SurfMesh%nMasterSides.EQ.0) RETURN
  isRestart = .FALSE.
  IF ( DoRestart ) THEN
    isRestart = .TRUE.
  END IF
  IF (.NOT.DoSurfModelAnalyze) RETURN
  OutputCounter = 2
  unit_index = 636
  doDistributionData=.FALSE.
  DO iPB=1,nPartBound
    IF (PartBound%SurfaceModel(iPB).EQ.3) THEN
      doDistributionData=.TRUE.
      EXIT
    END IF
  END DO
#if USE_MPI
  IF (SurfCOMM%MPIOutputRoot) THEN
#endif /*USE_MPI*/
    INQUIRE(UNIT   = unit_index , OPENED = isOpen)
    IF (.NOT.isOpen) THEN
      outfile = 'SurfaceAnalyze.csv'
!===================================================================================================================================
! Write Header
!===================================================================================================================================
      IF (isRestart .and. FILEEXISTS(outfile)) THEN
        OPEN(unit_index,file=TRIM(outfile),position="APPEND",status="OLD")
      ELSE
        OPEN(unit_index,file=TRIM(outfile))
        !--- insert header
        WRITE(unit_index,'(A8)',ADVANCE='NO') '001-TIME'
        IF (CalcCollCounter) THEN
          CALL WriteDataHeaderInfo(unit_index,'nSurfColl-Spec',OutputCounter,nSpecies)
          CALL WriteDataHeaderInfo(unit_index,'N_Ads-Spec',OutputCounter,nSpecies)
        END IF
        IF (CalcDesCounter) THEN
          CALL WriteDataHeaderInfo(unit_index,'N_Des-Spec',OutputCounter,nSpecies)
        END IF
        IF (CalcAdsProb) THEN
          CALL WriteDataHeaderInfo(unit_index,'MeanAdsProb-Spec',OutputCounter,nSpecies)
        END IF
        IF (CalcDesProb) THEN
          CALL WriteDataHeaderInfo(unit_index,'MeanDesProb-Spec',OutputCounter,nSpecies)
        END IF
#if (PP_TimeDiscMethod==42) || (PP_TimeDiscMethod==4)
        IF (doDistributionData) THEN
          IF (CalcSurfNumSpec) THEN
            CALL WriteDataHeaderInfo(unit_index,'nSimPart-Wall-Spec',OutputCounter,nSpecies)
            CALL WriteDataHeaderInfo(unit_index,'nSurfPart-Wall-Spec',OutputCounter,nSpecies)
          END IF
          IF (CalcSurfCoverage) THEN
            CALL WriteDataHeaderInfo(unit_index,'Surf-Cov',OutputCounter,nSpecies)
          END IF
#if (PP_TimeDiscMethod==42)
          IF (CalcAccomodation) THEN
            CALL WriteDataHeaderInfo(unit_index,'Alpha-Spec',OutputCounter,nSpecies)
          END IF
          IF (CalcAdsorbProb) THEN
            CALL WriteDataHeaderInfo(unit_index,'Prob_adsorption-Spec',OutputCounter,nSpecies)
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3)',ADVANCE='NO') OutputCounter,'-P_Molec-Adsorb-Spec-',iSpec
              OutputCounter = OutputCounter + 1
              DO iCase = 1,Adsorption%DissNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-P_Dissoc-Spec-',iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
              DO iCase = 1, Adsorption%RecombNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-P_ER-Spec-',iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
            END DO
          END IF
          IF (CalcAdsorbnu) THEN
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3)',ADVANCE='NO') OutputCounter,'-nu-Adsorb-Spec-', iSpec
              OutputCounter = OutputCounter + 1
              DO iCase = 1,Adsorption%DissNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-nu-diss-Spec-', iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
              DO iCase = 1,Adsorption%RecombNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-nu-ER-Spec-', iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
            END DO
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3)',ADVANCE='NO') OutputCounter,'-Proper-nu-Adsorb-Spec-', iSpec
              OutputCounter = OutputCounter + 1
              DO iCase = 1,Adsorption%DissNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-Proper-nu-diss-Spec-', iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
              DO iCase = 1,Adsorption%RecombNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-Proper-nu-ER-Spec-', iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
            END DO
          END IF
          IF (CalcAdsorbE) THEN
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3)',ADVANCE='NO') OutputCounter,'-E-Adsorb-Spec-', iSpec
              OutputCounter = OutputCounter + 1
              DO iCase = 1,Adsorption%DissNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-E-diss-Spec-', iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
              DO iCase = 1,Adsorption%RecombNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-E-ER-Spec-', iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
            END DO
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3)',ADVANCE='NO') OutputCounter,'-Proper-E-Adsorb-Spec-', iSpec
              OutputCounter = OutputCounter + 1
              DO iCase = 1,Adsorption%DissNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-Proper-E-diss-Spec-', iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
              DO iCase = 1,Adsorption%RecombNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-Proper-E-ER-Spec-', iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
            END DO
          END IF
          IF (CalcSurfProb) THEN
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3)',ADVANCE='NO') OutputCounter,'-P-SurfDesorb-Molec-Spec-', iSpec
              OutputCounter = OutputCounter + 1
              DO iCase = 1, Adsorption%DissNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-P-SurfDissoc-Spec-',iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
              DO iCase = 1, Adsorption%RecombNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-P-SurfLH-Spec-',iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
            END DO
            CALL WriteDataHeaderInfo(unit_index,'P-Surfexch-Case',OutputCounter,Adsorption%NumOfExchReact)
          END IF
          IF (CalcSurfnu) THEN
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3)',ADVANCE='NO') OutputCounter,'-nu-Desorb-Spec-', iSpec
              OutputCounter = OutputCounter + 1
              DO iCase = 1,Adsorption%DissNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-nu-Diss-Spec-', iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
              DO iCase = 1,Adsorption%RecombNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-nu-LH-Spec-', iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
            END DO
            CALL WriteDataHeaderInfo(unit_index,'nu-Exch-Reaction',OutputCounter,Adsorption%NumOfExchReact)
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3)',ADVANCE='NO') OutputCounter,'-Proper-nu-Desorb-Spec-', iSpec
              OutputCounter = OutputCounter + 1
              DO iCase = 1,Adsorption%DissNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-Proper-nu-Diss-Spec-', iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
              DO iCase = 1,Adsorption%RecombNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-Proper-nu-LH-Spec-', iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
            END DO
            CALL WriteDataHeaderInfo(unit_index,'Proper-nu-Exch-Reaction',OutputCounter,Adsorption%NumOfExchReact)
          END IF
          IF (CalcSurfE) THEN
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3)',ADVANCE='NO') OutputCounter,'-E-Desorb-Spec-', iSpec
              OutputCounter = OutputCounter + 1
              DO iCase = 1,Adsorption%DissNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-E-Diss-Spec-', iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
              DO iCase = 1,Adsorption%RecombNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-E-LH-Spec-', iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
            END DO
            CALL WriteDataHeaderInfo(unit_index,'E-Exch-Reaction',OutputCounter,Adsorption%NumOfExchReact)
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3)',ADVANCE='NO') OutputCounter,'-Proper-E-Desorb-Spec-', iSpec
              OutputCounter = OutputCounter + 1
              DO iCase = 1,Adsorption%DissNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-Proper-E-Diss-Spec-', iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
              DO iCase = 1,Adsorption%RecombNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-Proper-E-LH-Spec-', iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
            END DO
            CALL WriteDataHeaderInfo(unit_index,'Proper-E-Exch-Reaction',OutputCounter,Adsorption%NumOfExchReact)
          END IF
          IF (CalcHeatFlux) THEN
            CALL WriteDataHeaderInfo(unit_index,'Adsorption-HeatFlux-Spec',OutputCounter,nSpecies)
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3)',ADVANCE='NO') OutputCounter,'-AdsCount-Spec-', iSpec
              OutputCounter = OutputCounter + 1
              DO iCase = 1,Adsorption%DissNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-Count-Diss-Spec-', iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
              DO iCase = 1,Adsorption%RecombNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-Count-ER-Spec-', iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
            END DO
            CALL WriteDataHeaderInfo(unit_index,'Desorption-HeatFlux-Spec',OutputCounter,nSpecies)
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3)',ADVANCE='NO') OutputCounter,'-DesCount-Spec-', iSpec
              OutputCounter = OutputCounter + 1
              DO iCase = 1,Adsorption%DissNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-Count-Diss-Spec-', iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
              DO iCase = 1,Adsorption%RecombNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3)',ADVANCE='NO') &
                    OutputCounter,'-Count-LH-Spec-', iSpec,'-Reaction-', iCase
                OutputCounter = OutputCounter + 1
              END DO
            END DO
            CALL WriteDataHeaderInfo(unit_index,'Count-Exch-Reaction',OutputCounter,Adsorption%NumOfExchReact)
          END IF
          IF (Adsorption%TPD) THEN
            CALL WriteDataHeaderInfo(unit_index,'WallTemp',OutputCounter,1)
          END IF
#endif
        END IF
#endif
        WRITE(unit_index,'(A)') ''
      END IF
    END IF
#if USE_MPI
  END IF
#endif /*USE_MPI*/

!===================================================================================================================================
! Analyze Routines
!===================================================================================================================================
#if (PP_TimeDiscMethod==42)
IF (CalcAccomodation) CALL GetAccCoeff(Accomodation) ! called here because uses wallcollcount that is reset in getcollcounter
#endif
IF (CalcCollCounter) CALL GetCollCounter(SurfCollNum,AdsorptionNum) !collision coutner is reset here
IF (CalcDesCounter) CALL GetDesCounter(DesorptionNum)
IF (CalcAdsProb) CALL GetAdsProb(MeanAdsorptionProb)
IF (CalcDesProb) CALL GetDesProb(MeanDesorptionProb)
#if (PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==42)
IF (doDistributionData) THEN
  IF (CalcSurfNumSpec.OR.CalcSurfCoverage) CALL GetWallNumSpec(WallNumSpec,WallCoverage,WallNumSpec_SurfDist)
#if (PP_TimeDiscMethod==42)
  IF (CalcAdsorbRates) THEN
    SDEALLOCATE(AdsorptionReactRate)
    SDEALLOCATE(AdsorptionActE)
    SDEALLOCATE(ProperAdsorptionActE)
    SDEALLOCATE(Adsorptionnu)
    SDEALLOCATE(ProperAdsorptionnu)
    ALLOCATE(AdsorptionReactRate(1:nSpecies*(Adsorption%ReactNum+1)))
    ALLOCATE(AdsorptionActE(1:nSpecies*(Adsorption%ReactNum+1)))
    ALLOCATE(ProperAdsorptionActE(1:nSpecies*(Adsorption%ReactNum+1)))
    ALLOCATE(Adsorptionnu(1:nSpecies*(Adsorption%ReactNum+1)))
    ALLOCATE(ProperAdsorptionnu(1:nSpecies*(Adsorption%ReactNum+1)))
    CALL GetAdsRates(AdsorptionReactRate,AdsorptionActE,ProperAdsorptionActE,Adsorptionnu,ProperAdsorptionnu)
  END IF
  IF (CalcSurfRates) THEN
    SDEALLOCATE(SurfReactRate)
    SDEALLOCATE(SurfaceActE)
    SDEALLOCATE(ProperSurfaceActE)
    SDEALLOCATE(Surfacenu)
    SDEALLOCATE(ProperSurfacenu)
    ALLOCATE(SurfReactRate(1:nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact))
    ALLOCATE(SurfaceActE(1:nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact))
    ALLOCATE(ProperSurfaceActE(1:nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact))
    ALLOCATE(Surfacenu(1:nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact))
    ALLOCATE(ProperSurfacenu(1:nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact))
    CALL GetSurfRates(SurfReactRate,SurfaceActE,ProperSurfaceActE,Surfacenu,ProperSurfacenu)
  END IF
  IF (CalcHeatFlux) THEN
    SDEALLOCATE(HeatFlux)
    SDEALLOCATE(AdsReactCount)
    SDEALLOCATE(DesReactCount)
    ALLOCATE(HeatFlux(1:2,1:nSpecies))
    ALLOCATE(AdsReactCount(1:nSpecies*(Adsorption%ReactNum+1)))
    ALLOCATE(DesReactCount(1:nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact))
    CALL GetSurfHeatFluxes(HeatFlux,AdsReactCount,DesReactCount)
  END IF
#endif
END IF
#endif
! ResetAllInfo
IF (ALLOCATED(SurfModel%Info)) THEN
  DO iSpec = 1,nSpecies
    SurfModel%Info(iSpec)%WallCollCount = 0
    SurfModel%Info(iSpec)%NumOfAds = 0
    SurfModel%Info(iSpec)%NumOfDes = 0
    SurfModel%Info(iSpec)%MeanProbAds = 0.
    SurfModel%Info(iSpec)%MeanProbAdsCount = 0
    SurfModel%Info(iSpec)%MeanProbDes = 0.
    SurfModel%Info(iSpec)%MeanProbDesCount = 0
#if (PP_TimeDiscMethod==42)
    SurfModel%Info(iSpec)%WallSpecNumCount = 0
    SurfModel%Info(iSpec)%Accomodation = 0.
#endif
  END DO
END IF
!===================================================================================================================================
! Output Analyzed variables
!===================================================================================================================================
#if USE_MPI
IF (SurfCOMM%MPIOutputRoot) THEN
#endif /*USE_MPI*/
  WRITE(unit_index,'(E23.16E3)',ADVANCE='NO') Time
  IF (CalcCollCounter) THEN
    CALL WriteDataInfo(unit_index,nSpecies,IntegerArray=SurfCollNum(:))
    CALL WriteDataInfo(unit_index,nSpecies,IntegerArray=AdsorptionNum(:))
  END IF
  IF (CalcDesCounter) THEN
    CALL WriteDataInfo(unit_index,nSpecies,IntegerArray=DesorptionNum(:))
  END IF
  IF (CalcAdsProb) THEN
    CALL WriteDataInfo(unit_index,nSpecies,RealArray=MeanAdsorptionProb(:))
  END IF
  IF (CalcDesProb) THEN
    CALL WriteDataInfo(unit_index,nSpecies,RealArray=MeanDesorptionProb(:))
  END IF
#if ((PP_TimeDiscMethod==42) || (PP_TimeDiscMethod==4))
! output for adsorption
IF (doDistributiondata) THEN
      IF (CalcSurfNumSpec) THEN
        CALL WriteDataInfo(unit_index,nSpecies,IntegerK8Array=WallNumSpec(:))
        CALL WriteDataInfo(unit_index,nSpecies,IntegerK8Array=WallNumSpec_SurfDist(:))
      END IF
      IF (CalcSurfCoverage) THEN
        CALL WriteDataInfo(unit_index,nSpecies,RealArray=WallCoverage(:))
      END IF
#if (PP_TimeDiscMethod==42)
      IF (CalcAccomodation) THEN
        CALL WriteDataInfo(unit_index,nSpecies,RealArray=Accomodation(:))
      END IF
      IF (CalcAdsorbProb) THEN
        CALL WriteDataInfo(unit_index,nSpecies*(Adsorption%ReactNum+1) ,RealArray=AdsorptionReactRate(:))
      END IF
      IF (CalcAdsorbnu) THEN
        CALL WriteDataInfo(unit_index,nSpecies*(Adsorption%ReactNum+1) ,RealArray=Adsorptionnu(:))
        CALL WriteDataInfo(unit_index,nSpecies*(Adsorption%ReactNum+1) ,RealArray=ProperAdsorptionnu(:))
      END IF
      IF (CalcAdsorbE) THEN
        CALL WriteDataInfo(unit_index,nSpecies*(Adsorption%ReactNum+1) ,RealArray=AdsorptionActE(:))
        CALL WriteDataInfo(unit_index,nSpecies*(Adsorption%ReactNum+1) ,RealArray=ProperAdsorptionActE(:))
      END IF
      IF (CalcSurfProb) THEN
        CALL WriteDataInfo(unit_index,nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact ,RealArray=SurfReactRate(:))
      END IF
      IF (CalcSurfnu) THEN
        CALL WriteDataInfo(unit_index,nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact ,RealArray=Surfacenu(:))
        CALL WriteDataInfo(unit_index,nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact ,RealArray=ProperSurfacenu(:))
      END IF
      IF (CalcSurfE) THEN
        CALL WriteDataInfo(unit_index,nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact ,RealArray=SurfaceActE(:))
        CALL WriteDataInfo(unit_index,nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact ,RealArray=ProperSurfaceActE(:))
      END IF
      IF (CalcHeatFlux) THEN
        CALL WriteDataInfo(unit_index,nSpecies,RealArray=HeatFlux(1,:))
        CALL WriteDataInfo(unit_index,nSpecies*(Adsorption%ReactNum+1) ,RealArray=AdsReactCount(:))
        CALL WriteDataInfo(unit_index,nSpecies,RealArray=HeatFlux(2,:))
        CALL WriteDataInfo(unit_index,nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact ,RealArray=DesReactCount(:))
      END IF
      IF (Adsorption%TPD) THEN
        CALL WriteDataInfo(unit_index,1,RealScalar=Adsorption%TPD_Temp)
      END IF
#endif /*(PP_TimeDiscMethod==42)*/
    END IF
#endif /*(PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==42)*/
    WRITE(unit_index,'(A)') ''
#if USE_MPI
  END IF
#endif /*USE_MPI*/
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE AnalyzeSurface


SUBROUTINE WriteDataHeaderInfo(unit_index,AttribName,OutputCounter,LoopSize)
!===================================================================================================================================
!> writes OutputCounter-AttribNamestring-iLoop into WRITEFORMAT output
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)          :: unit_index
CHARACTER(LEN=*),INTENT(IN) :: AttribName
INTEGER,INTENT(INOUT)       :: OutputCounter
INTEGER,INTENT(IN)          :: LoopSize
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                     :: iLoop
!===================================================================================================================================
DO iLoop = 1, LoopSize
  WRITE(unit_index,'(A1)',ADVANCE='NO') ','
  WRITE(unit_index,'(I3.3,A,A,A,I3.3)',ADVANCE='NO') OutputCounter,'-',TRIM(AttribName),'-',iLoop
  OutputCounter = OutputCounter + 1
END DO
END SUBROUTINE WriteDataHeaderInfo


SUBROUTINE WriteDataInfo(unit_index,nVal,RealScalar,IntegerScalar,StrScalar,LogicalScalar, &
                                  RealArray,IntegerArray,IntegerK8Array,StrArray)
!===================================================================================================================================
!> writes INPUTData into unit_index output
!> only one data input should be given at a time
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER           ,INTENT(IN)          :: unit_index
INTEGER           ,INTENT(IN)          :: nVal
REAL              ,INTENT(IN),OPTIONAL :: RealScalar
INTEGER           ,INTENT(IN),OPTIONAL :: IntegerScalar
CHARACTER(LEN=*)  ,INTENT(IN),OPTIONAL :: StrScalar
REAL              ,INTENT(IN),OPTIONAL :: RealArray(nVal)
INTEGER           ,INTENT(IN),OPTIONAL :: IntegerArray(nVal)
INTEGER(KIND=8)   ,INTENT(IN),OPTIONAL :: IntegerK8Array(nVal)
CHARACTER(LEN=255),INTENT(IN),OPTIONAL :: StrArray(nVal)
LOGICAL           ,INTENT(IN),OPTIONAL :: LogicalScalar
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                     :: iLoop
!===================================================================================================================================
IF(PRESENT(RealArray)) THEN
  DO iLoop = 1, nVal
    WRITE (unit_index, CSVFORMAT, ADVANCE='NO') ',',RealArray(iLoop)
  END DO
END IF
IF(PRESENT(RealScalar)) THEN
  WRITE (unit_index, CSVFORMAT, ADVANCE='NO') ',',RealScalar
END IF

IF(PRESENT(IntegerArray)) THEN
  DO iLoop = 1, nVal
    WRITE (unit_index, CSVFORMAT, ADVANCE='NO') ',',REAL(IntegerArray(iLoop))
  END DO
END IF

IF(PRESENT(IntegerK8Array)) THEN
  DO iLoop = 1, nVal
    WRITE (unit_index, CSVFORMAT, ADVANCE='NO') ',',REAL(IntegerK8Array(iLoop))
  END DO
END IF

IF(PRESENT(IntegerScalar)) THEN
  WRITE (unit_index, CSVFORMAT, ADVANCE='NO') ',',REAL(IntegerScalar)
END IF

IF(PRESENT(StrArray)) THEN
  DO iLoop = 1, nVal
    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
    WRITE(unit_index,'(A)',ADVANCE='NO') TRIM(StrArray(iLoop))
  END DO
END IF

IF(PRESENT(StrScalar)) THEN
  WRITE(unit_index,'(A1)',ADVANCE='NO') ','
  WRITE(unit_index,'(A)',ADVANCE='NO') TRIM(StrScalar)
END IF

IF(PRESENT(LogicalScalar)) THEN
  WRITE(unit_index,'(A1)',ADVANCE='NO') ','
  WRITE(unit_index,'(I1)',ADVANCE='NO') LogicalScalar
END IF
END SUBROUTINE WriteDataInfo


SUBROUTINE GetCollCounter(SurfCollNum,AdsorbNum)
!===================================================================================================================================
!> Calculates species counters for surface collisions
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_SurfaceModel_Vars      ,ONLY: SurfModel
#if USE_MPI
USE MOD_Particle_Boundary_Vars ,ONLY: SurfCOMM
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(OUT) :: SurfCollNum(nSpecies), AdsorbNum(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iSpec
#if USE_MPI
INTEGER            :: ADN(nSpecies)
#endif /*USE_MPI*/
!===================================================================================================================================
DO iSpec = 1,nSpecies
  SurfCollNum(iSpec) = SurfModel%Info(iSpec)%WallCollCount
  AdsorbNum(iSpec) = SurfModel%Info(iSpec)%NumOfAds
END DO
#if USE_MPI
IF (SurfCOMM%MPIOutputRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,SurfCollNum ,nSpecies,MPI_INTEGER,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,AdsorbNum   ,nSpecies,MPI_INTEGER,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
ELSE
  CALL MPI_REDUCE(SurfCollNum ,ADN         ,nSpecies,MPI_INTEGER,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
  CALL MPI_REDUCE(AdsorbNum   ,ADN         ,nSpecies,MPI_INTEGER,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
END IF
#endif /*USE_MPI*/
END SUBROUTINE GetCollCounter


SUBROUTINE GetDesCounter(DesorbNum)
!===================================================================================================================================
!> Calculate counter for desorption for each species
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_SurfaceModel_Vars      ,ONLY: SurfModel
#if USE_MPI
USE MOD_Particle_Boundary_Vars ,ONLY: SurfCOMM
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(OUT) :: DesorbNum(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iSpec
#if USE_MPI
INTEGER            :: DEN(nSpecies)
#endif /*USE_MPI*/
!===================================================================================================================================
DO iSpec = 1,nSpecies
  DesorbNum(iSpec) = SurfModel%Info(iSpec)%NumOfDes
END DO
#if USE_MPI
IF (SurfCOMM%MPIOutputRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,DesorbNum   ,nSpecies,MPI_INTEGER,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
ELSE
  CALL MPI_REDUCE(DesorbNum   ,DEN         ,nSpecies,MPI_INTEGER,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
END IF
#endif /*USE_MPI*/
END SUBROUTINE GetDesCounter


SUBROUTINE GetAdsProb(MeanAdsorptionProb)
!===================================================================================================================================
!> Calculate mean adsorption probability for adsorption for each species
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_SurfaceModel_Vars      ,ONLY: SurfModel
#if USE_MPI
USE MOD_Particle_Boundary_Vars ,ONLY: SurfCOMM
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT) :: MeanAdsorptionProb(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iSpec
#if USE_MPI
REAL               :: MAP(nSpecies)
#endif /*USE_MPI*/
!===================================================================================================================================
DO iSpec = 1,nSpecies
  IF (SurfModel%Info(iSpec)%MeanProbAdsCount.GT.0) THEN
    MeanAdsorptionProb(iSpec) = SurfModel%Info(iSpec)%MeanProbAds / SurfModel%Info(iSpec)%MeanProbAdsCount
  ELSE
    MeanAdsorptionProb(iSpec) = 0.
  END IF
END DO
#if USE_MPI
IF (SurfCOMM%MPIOutputRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,MeanAdsorptionProb   ,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
  MeanAdsorptionProb = MeanAdsorptionProb / REAL(SurfCOMM%nOutputProcs)
ELSE
  CALL MPI_REDUCE(MeanAdsorptionProb   ,MAP         ,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
END IF
#endif /*USE_MPI*/
END SUBROUTINE GetAdsProb


SUBROUTINE GetDesProb(MeanDesorptionProb)
!===================================================================================================================================
!> Calculate mean adsorption probability for adsorption for each species
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_SurfaceModel_Vars      ,ONLY: SurfModel
#if USE_MPI
USE MOD_Particle_Boundary_Vars ,ONLY: SurfCOMM
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT) :: MeanDesorptionProb(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iSpec
#if USE_MPI
REAL               :: MAP(nSpecies)
#endif /*USE_MPI*/
!===================================================================================================================================
DO iSpec = 1,nSpecies
  IF (SurfModel%Info(iSpec)%MeanProbDesCount.GT.0) THEN
    MeanDesorptionProb(iSpec) = SurfModel%Info(iSpec)%MeanProbDes / SurfModel%Info(iSpec)%MeanProbDesCount
  ELSE
    MeanDesorptionProb(iSpec) = 0.
  END IF
END DO
#if USE_MPI
IF (SurfCOMM%MPIOutputRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,MeanDesorptionProb   ,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
  MeanDesorptionProb = MeanDesorptionProb / REAL(SurfCOMM%nOutputProcs)
ELSE
  CALL MPI_REDUCE(MeanDesorptionProb   ,MAP         ,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
END IF
#endif /*USE_MPI*/
END SUBROUTINE GetDesProb


#if (PP_TimeDiscMethod==42) || (PP_TimeDiscMethod==4)
SUBROUTINE GetWallNumSpec(WallNumSpec,WallCoverage,WallNumSpec_SurfDist)
!===================================================================================================================================
! Calculate number of wallparticles for all species
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars                 ,ONLY: BC
USE MOD_Particle_Vars             ,ONLY: Species, PartSpecies, PDM, nSpecies, KeepWallParticles
USE MOD_SurfaceModel_Analyze_Vars
USE MOD_SurfaceModel_Vars         ,ONLY: Adsorption, SurfDistInfo
USE MOD_Particle_Boundary_Vars    ,ONLY: nSurfSample, SurfMesh, PartBound
#if USE_MPI
USE MOD_Particle_Boundary_Vars    ,ONLY: SurfCOMM
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=8), INTENT(OUT)    :: WallNumSpec(nSpecies),WallNumSpec_SurfDist(nSpecies)
REAL           , INTENT(OUT)    :: WallCoverage(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i, iSpec, iSurfSide, p, q, SideID, PartBoundID
REAL                            :: SurfPart
REAL                            :: Coverage(nSpecies)
#if USE_MPI
REAL                            :: RD(nSpecies)
INTEGER(KIND=8)                 :: IDR(nSpecies), ID1(nSpecies), ID2(nSpecies), ID3(nSpecies*2)
#endif /*USE_MPI*/
INTEGER                         :: Coord, AdsorbID, Surfpos, SpecID
INTEGER                         :: adsorbates(nSpecies)
REAL                            :: SubWallNumSpec(nSpecies), WallNumSpec_tmp(2*nSpecies)
!===================================================================================================================================
WallNumSpec = 0
WallNumSpec_SurfDist = 0
SurfPart = 0.
Coverage(:) = 0.
WallCoverage(:) = 0.
WallNumSpec_tmp = 0.
SubWallNumSpec = 0.

DO iSpec=1,nSpecies
DO iSurfSide=1,SurfMesh%nMasterSides
  SideID = SurfMesh%SurfIDToSideID(iSurfSide)
  PartboundID = PartBound%MapToPartBC(BC(SideID))
  IF (PartBound%Reactive(PartboundID)) THEN
  DO q = 1,nSurfSample
    DO p = 1,nSurfSample
      Coverage(iSpec) = Coverage(iSpec) + Adsorption%Coverage(p,q,iSurfSide,iSpec)
      IF ((.NOT.KeepWallParticles) .AND. CalcSurfNumSpec) THEN
        SurfPart = REAL(INT(Adsorption%DensSurfAtoms(iSurfSide) * SurfMesh%SurfaceArea(p,q,iSurfSide),8))
!          WallNumSpec(iSpec) = WallNumSpec(iSpec) + INT( Adsorption%Coverage(p,q,iSurfSide,iSpec) &
!              * SurfPart/Species(iSpec)%MacroParticleFactor)
        ! calculate number of adsorbates for each species
        adsorbates = 0
        DO Coord = 1,3
        DO AdsorbID = 1,SurfDistInfo(p,q,iSurfSide)%nSites(Coord)-SurfDistInfo(p,q,iSurfSide)%SitesRemain(Coord)
          Surfpos = SurfDistInfo(p,q,iSurfSide)%AdsMap(Coord)%UsedSiteMap(SurfDistInfo(p,q,iSurfSide)%SitesRemain(Coord)+AdsorbID)
          SpecID = SurfDistInfo(p,q,iSurfSide)%AdsMap(Coord)%Species(Surfpos)
          adsorbates(SpecID) = adsorbates(SpecID) + 1
        END DO
        END DO
        ! discret simulated particles on surface distribution
        WallNumSpec_SurfDist(iSpec) = WallNumSpec_SurfDist(iSpec) + adsorbates(iSpec)
        ! simulated (gas) particles from discret surface distribution
        SubWallNumSpec(iSpec) = SubWallNumSpec(iSpec) + REAL(adsorbates(iSpec)) / REAL(SurfDistInfo(p,q,iSurfSide)%nSites(3))&
            * SurfPart/Species(iSpec)%MacroParticleFactor
        ! simulated gas particles safed in temporary arrays
        WallNumSpec_tmp(iSpec) = WallNumSpec_tmp(iSpec) + &
            ( SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec) / SurfDistInfo(p,q,iSurfSide)%nSites(3) &
            * SurfPart / Species(iSpec)%MacroParticleFactor )
        WallNumSpec_tmp(iSpec+nSpecies) = WallNumSpec_tmp(iSpec+nSpecies) + SurfDistInfo(p,q,iSurfSide)%desorbnum_tmp(iSpec)&
            - SurfDistInfo(p,q,iSurfSide)%reactnum_tmp(iSpec)
      END IF
    END DO
  END DO
  END IF
END DO
END DO
IF (CalcSurfCoverage .AND. SurfMesh%nMasterSides.GT.0) THEN
  WallCoverage(:) = Coverage(:) / (SurfMesh%nMasterSides*nSurfSample*nSurfSample)
END IF

#if USE_MPI
  IF (SurfCOMM%MPIOutputRoot) THEN
    IF (CalcSurfNumSpec)  THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,SubWallNumSpec      ,nSpecies  ,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
      CALL MPI_REDUCE(MPI_IN_PLACE,WallNumSpec_SurfDist,nSpecies  ,MPI_INTEGER,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
      CALL MPI_REDUCE(MPI_IN_PLACE,WallNumSpec_tmp     ,nSpecies*2,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
    END IF
    IF (CalcSurfCoverage) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,WallCoverage,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
      WallCoverage = WallCoverage / REAL(SurfCOMM%nOutputProcs)
    END IF
  ELSE
    IF (CalcSurfNumSpec) THEN
      CALL MPI_REDUCE(SubWallNumSpec      ,ID1,nSpecies  ,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
      CALL MPI_REDUCE(WallNumSpec_SurfDist,ID2,nSpecies  ,MPI_INTEGER,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
      CALL MPI_REDUCE(WallNumSpec_tmp     ,ID3,nSpecies*2,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
    END IF
    IF (CalcSurfCoverage) CALL MPI_REDUCE(WallCoverage,RD,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
  END IF
#endif /*USE_MPI*/

  IF (KeepWallParticles.AND.CalcSurfNumSpec) THEN
    DO i=1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(i) .AND. PDM%ParticleAtWall(i)) THEN
        WallNumSpec(PartSpecies(i)) = WallNumSpec(PartSpecies(i)) + 1
      END IF
    END DO
#if USE_MPI
  IF (SurfCOMM%MPIOutputRoot) THEN
    IF (CalcSurfNumSpec) CALL MPI_REDUCE(MPI_IN_PLACE,WallNumSpec,nSpecies,MPI_INTEGER,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
  ELSE
    IF (CalcSurfNumSpec) CALL MPI_REDUCE(WallNumSpec ,IDR        ,nSpecies,MPI_INTEGER,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
  END IF
#endif /*USE_MPI*/
  ELSE
    WallNumSpec = INT(SubWallNumSpec)+INT(WallNumSpec_tmp(1:nSpecies))+INT(WallNumSpec_tmp(nSpecies+1:nSpecies*2))
  END IF

END SUBROUTINE GetWallNumSpec

#if (PP_TimeDiscMethod==42)
SUBROUTINE GetAccCoeff(Accomodation)
!===================================================================================================================================
! Calculate accomodation rates for all species
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_SurfaceModel_Vars      ,ONLY: SurfModel
#if USE_MPI
USE MOD_Particle_Boundary_Vars ,ONLY: SurfCOMM
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL   , INTENT(OUT)            :: Accomodation(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSpec
#if USE_MPI
REAL                            :: AC(nSpecies)
#endif /*USE_MPI*/
!===================================================================================================================================

Accomodation(:) = 0.
DO iSpec = 1,nSpecies
  IF (SurfModel%Info(iSpec)%WallCollCount.GT.0) THEN
    Accomodation(iSpec) = SurfModel%Info(iSpec)%Accomodation / REAL(SurfModel%Info(iSpec)%WallCollCount)
  ELSE
    Accomodation(iSpec) = 0.
  END IF
END DO

#if USE_MPI
IF (SurfCOMM%MPIOutputRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,Accomodation,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
  Accomodation= Accomodation/ REAL(SurfCOMM%nOutputProcs)
ELSE
  CALL MPI_REDUCE(Accomodation,AC          ,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
END IF
#endif /*USE_MPI*/

END SUBROUTINE GetAccCoeff


SUBROUTINE GetAdsRates(ReactRate,AdsorbActE,ProperAdsorbActE,Adsorbnu,ProperAdsorbnu)
!===================================================================================================================================
! Calculate adsorption, desorption and accomodation rates for all species
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption, SurfModel
#if USE_MPI
USE MOD_Particle_Boundary_Vars ,ONLY: SurfCOMM
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL   , INTENT(OUT)            :: ReactRate(nSpecies*(Adsorption%ReactNum+1))
REAL   , INTENT(OUT)            :: AdsorbActE(nSpecies*(Adsorption%ReactNum+1))
REAL   , INTENT(OUT)            :: ProperAdsorbActE(nSpecies*(Adsorption%ReactNum+1))
REAL   , INTENT(OUT)            :: Adsorbnu(nSpecies*(Adsorption%ReactNum+1))
REAL   , INTENT(OUT)            :: ProperAdsorbnu(nSpecies*(Adsorption%ReactNum+1))
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSpec, iCase, iReact
#if USE_MPI
REAL                            :: RR(nSpecies*Adsorption%ReactNum)
#endif /*USE_MPI*/
!===================================================================================================================================

iCase = 1
DO iSpec = 1, nSpecies
  DO iReact = 1, Adsorption%ReactNum+1
    IF (SurfModel%ProperInfo(iSpec)%AdsReactCount(iReact).GT.0) THEN
      ReactRate(iCase) = SurfModel%ProperInfo(iSpec)%NumAdsReact(iReact) &
          / REAL(SurfModel%ProperInfo(iSpec)%AdsReactCount(iReact))
    ELSE
      ReactRate(iCase) = 0.
    END IF
    iCase = iCase + 1
  END DO
END DO

iCase = 1
DO iSpec = 1,nSpecies
  DO iReact = 1,Adsorption%ReactNum+1
    IF (SurfModel%ProperInfo(iSpec)%AdsReactCount(iReact).GT.0) THEN
      AdsorbActE(iCase) = SurfModel%ProperInfo(iSpec)%MeanAdsActE(iReact) &
          / REAL(SurfModel%ProperInfo(iSpec)%AdsReactCount(iReact))
      Adsorbnu(iCase) = SurfModel%ProperInfo(iSpec)%MeanAdsnu(iReact) &
          / REAL(SurfModel%ProperInfo(iSpec)%AdsReactCount(iReact))
    ELSE
      AdsorbActE(iCase) = 0.
      Adsorbnu(iCase)   = 0.
    END IF
    IF (SurfModel%ProperInfo(iSpec)%ProperAdsReactCount(iReact).GT.0) THEN
      ProperAdsorbActE(iCase) = SurfModel%ProperInfo(iSpec)%ProperAdsActE(iReact) &
          / REAL(SurfModel%ProperInfo(iSpec)%ProperAdsReactCount(iReact))
      ProperAdsorbnu(iCase) = SurfModel%ProperInfo(iSpec)%ProperAdsnu(iReact) &
          / REAL(SurfModel%ProperInfo(iSpec)%ProperAdsReactCount(iReact))
    ELSE
      ProperAdsorbActE(iCase) = 0.
      ProperAdsorbnu(iCase) = 0.
    END IF
    iCase = iCase + 1
  END DO
END DO

#if USE_MPI
IF (SurfCOMM%MPIOutputRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,ReactRate,nSpecies*(Adsorption%ReactNum+1),MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,AdsorbActE  ,nSpecies*Adsorption%ReactNum,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
  ReactRate  = ReactRate   / REAL(SurfCOMM%nOutputProcs)
  AdsorbActE = AdsorbActE  / REAL(SurfCOMM%nOutputProcs)
ELSE
  CALL MPI_REDUCE(ReactRate   ,RR       ,nSpecies*(Adsorption%ReactNum+1),MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
  CALL MPI_REDUCE(AdsorbActE  ,RR          ,nSpecies*Adsorption%ReactNum,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
END IF
#endif /*USE_MPI*/

DO iSpec = 1,nSpecies
  SurfModel%ProperInfo(iSpec)%NumAdsReact(:) = 0.
  SurfModel%ProperInfo(iSpec)%AdsReactCount(:) = 0
  SurfModel%ProperInfo(iSpec)%MeanAdsActE(:) = 0.
  SurfModel%ProperInfo(iSpec)%ProperAdsActE(:) = 0.
  SurfModel%ProperInfo(iSpec)%MeanAdsnu(:) = 0.
  SurfModel%ProperInfo(iSpec)%ProperAdsnu(:) = 0.
END DO

END SUBROUTINE GetAdsRates


SUBROUTINE GetSurfRates(ReactRate,SurfaceActE,ProperSurfaceActE,Surfacenu,ProperSurfacenu)
!===================================================================================================================================
! Calculate adsorption, desorption and accomodation rates for all species
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption, SurfModel
#if USE_MPI
USE MOD_Particle_Boundary_Vars ,ONLY: SurfCOMM
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL   , INTENT(OUT)            :: ReactRate(nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact)
REAL   , INTENT(OUT)            :: SurfaceActE(nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact)
REAL   , INTENT(OUT)            :: ProperSurfaceActE(nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact)
REAL   , INTENT(OUT)            :: Surfacenu(nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact)
REAL   , INTENT(OUT)            :: ProperSurfacenu(nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSpec, iReact, iCase
#if USE_MPI
INTEGER                         :: commSize
REAL                            :: RR(nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact)
#endif /*USE_MPI*/
!===================================================================================================================================

IF (.NOT.DSMC%ReservoirRateStatistic) THEN
  iCase = 1
  DO iSpec = 1,nSpecies
    DO iReact = 1, Adsorption%ReactNum+1
      IF (SurfModel%ProperInfo(iSpec)%SurfReactCount(iReact).GT.0) THEN
        ReactRate(iCase) = SurfModel%ProperInfo(iSpec)%NumSurfReact(iReact) &
            / REAL(SurfModel%ProperInfo(iSpec)%SurfReactCount(iReact))
      ELSE
        ReactRate(iCase) = 0.
      END IF
      iCase = iCase + 1
    END DO
  END DO
END IF

iCase = 1
DO iSpec = 1,nSpecies
  DO iReact = 1,Adsorption%ReactNum+1
    IF (SurfModel%ProperInfo(iSpec)%SurfReactCount(iReact).GT.0) THEN
      SurfaceActE(iCase) = SurfModel%ProperInfo(iSpec)%MeanSurfActE(iReact) &
          / REAL(SurfModel%ProperInfo(iSpec)%SurfReactCount(iReact))
    ELSE
      SurfaceActE(iCase) = 0.
    END IF
    iCase = iCase + 1
  END DO
END DO

iCase = 1
DO iSpec = 1,nSpecies
  DO iReact = 1,Adsorption%ReactNum+1
    IF (SurfModel%ProperInfo(iSpec)%ProperSurfReactCount(iReact).GT.0) THEN
      ProperSurfaceActE(iCase) = SurfModel%ProperInfo(iSpec)%ProperSurfActE(iReact) &
          / REAL(SurfModel%ProperInfo(iSpec)%ProperSurfReactCount(iReact))
    ELSE
      ProperSurfaceActE(iCase) = 0.
    END IF
    iCase = iCase + 1
  END DO
END DO

iCase = 1
DO iSpec = 1,nSpecies
  DO iReact = 1,Adsorption%ReactNum+1
    IF (SurfModel%ProperInfo(iSpec)%SurfReactCount(iReact).GT.0) THEN
      Surfacenu(iCase) = SurfModel%ProperInfo(iSpec)%MeanSurfnu(iReact) &
          / REAL(SurfModel%ProperInfo(iSpec)%SurfReactCount(iReact))
    ELSE
      Surfacenu(iCase) = 0.
    END IF
    iCase = iCase + 1
  END DO
END DO

iCase = 1
DO iSpec = 1,nSpecies
  DO iReact = 1,Adsorption%ReactNum+1
    IF (SurfModel%ProperInfo(iSpec)%ProperSurfReactCount(iReact).GT.0) THEN
      ProperSurfacenu(iCase) = SurfModel%ProperInfo(iSpec)%ProperSurfnu(iReact) &
          / REAL(SurfModel%ProperInfo(iSpec)%ProperSurfReactCount(iReact))
    ELSE
      ProperSurfacenu(iCase) = 0.
    END IF
    iCase = iCase + 1
  END DO
END DO

#if USE_MPI
commSize = nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact
IF (SurfCOMM%MPIOutputRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE ,ReactRate        ,commSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE ,SurfaceActE      ,commSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE ,ProperSurfaceActE,commSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
  ReactRate   = ReactRate   / REAL(SurfCOMM%nOutputProcs)
  SurfaceActE = SurfaceActE / REAL(SurfCOMM%nOutputProcs)
  ProperSurfaceActE = ProperSurfaceActE / REAL(SurfCOMM%nOutputProcs)
ELSE
  CALL MPI_REDUCE(ReactRate         ,RR ,commSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
  CALL MPI_REDUCE(SurfaceActE       ,RR ,commSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
  CALL MPI_REDUCE(ProperSurfaceActE ,RR ,commSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
END IF
#endif /*USE_MPI*/

DO iSpec = 1,nSpecies
  SurfModel%ProperInfo(iSpec)%MeanSurfActE = 0.
  SurfModel%ProperInfo(iSpec)%ProperSurfActE = 0.
  SurfModel%ProperInfo(iSpec)%MeanSurfnu = 0.
  SurfModel%ProperInfo(iSpec)%ProperSurfnu = 0.
  SurfModel%ProperInfo(iSpec)%NumSurfReact = 0.
  SurfModel%ProperInfo(iSpec)%SurfReactCount = 0
  SurfModel%ProperInfo(iSpec)%ProperSurfReactCount = 0
END DO

END SUBROUTINE GetSurfRates


SUBROUTINE GetSurfHeatFluxes(HeatFlux,AdsReactCount,DesReactCount)
!===================================================================================================================================
!> Calculate heat fluxes on surface resulting from enthalpy of reaction for all species
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption, SurfModel
USE MOD_Particle_Boundary_Vars ,ONLY: SurfMesh
#if USE_MPI
USE MOD_Particle_Boundary_Vars ,ONLY: SurfCOMM
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL   , INTENT(OUT)            :: HeatFlux(1:2,1:nSpecies)
REAL   , INTENT(OUT)            :: AdsReactCount(1:nSpecies*(Adsorption%ReactNum+1))
REAL   , INTENT(OUT)            :: DesReactCount(1:nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSpec, iCase, iReact
#if USE_MPI
INTEGER                         :: commSize1, commSize2
REAL                            :: HE(1:2,1:nSpecies)
REAL                            :: RA(1:nSpecies*(Adsorption%ReactNum+1))
REAL                            :: RD(1:nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact)
#endif /*USE_MPI*/
!===================================================================================================================================

IF(SurfMesh%SurfOnProc)THEN
  ! analyze heatflux to surface for each species
  DO iSpec = 1,nSpecies
    HeatFlux(1,iSpec) = -SurfModel%ProperInfo(iSpec)%HeatFlux(1)
    HeatFlux(2,iSpec) = -SurfModel%ProperInfo(iSpec)%HeatFlux(2)
  END DO
  ! analyze number of reactions for each species and each reaction
  iCase = 1
  DO iSpec = 1,nSpecies
    DO iReact = 1,Adsorption%ReactNum+1
      AdsReactCount(iCase) = SurfModel%ProperInfo(iSpec)%HeatFluxAdsCount(iReact)
      iCase = iCase + 1
    END DO
  END DO
  iCase = 1
  DO iSpec = 1,nSpecies
    DO iReact = 1,Adsorption%ReactNum+1
      DesReactCount(iCase) = SurfModel%ProperInfo(iSpec)%HeatFluxDesCount(iReact)
      iCase = iCase + 1
    END DO
  END DO
ELSE
  HeatFlux(:,:) = 0.
  AdsReactCount(:) = 0.
  DesReactCount(:) = 0.
END IF

#if USE_MPI
commSize1 = nSpecies*(Adsorption%ReactNum+1)
commSize2 = nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact
IF (SurfCOMM%MPIOutputRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,HeatFlux(1,:),nSpecies ,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,HeatFlux(2,:),nSpecies ,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,AdsReactCount,commSize1,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,DesReactCount,commSize2,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
ELSE
  CALL MPI_REDUCE(HeatFlux(1,:),HE(1,:)     ,nSpecies ,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
  CALL MPI_REDUCE(HeatFlux(2,:),HE(2,:)     ,nSpecies ,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
  CALL MPI_REDUCE(AdsReactCount,RA          ,commSize1,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
  CALL MPI_REDUCE(DesReactCount,RD          ,commSize2,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%OutputCOMM,IERROR)
END IF
#endif /*USE_MPI*/

IF(SurfMesh%SurfOnProc)THEN
  DO iSpec = 1,nSpecies
    SurfModel%ProperInfo(iSpec)%HeatFlux(:) = 0.
    SurfModel%ProperInfo(iSpec)%HeatFluxAdsCount(:) = 0.
    SurfModel%ProperInfo(iSpec)%HeatFluxDesCount(:) = 0.
  END DO
END IF

END SUBROUTINE GetSurfHeatFluxes


SUBROUTINE AnalyzeSurfRates(AnalyzeCase,SpecID,ReactionID,EAct,nuReact,Probability)
!===================================================================================================================================
!> Routine analyzing reaction rates at surfaces for SMCR
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars         ,ONLY: DSMC
USE MOD_SurfaceModel_Vars ,ONLY: SurfModel
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)            :: AnalyzeCase      !1: meansurfrate, 2: propersurfrate
INTEGER, INTENT(IN)            :: SpecID
INTEGER, INTENT(IN)            :: ReactionID
REAL, INTENT(IN)               :: EAct
REAL, INTENT(IN)               :: nuReact
REAL, INTENT(IN)               :: Probability
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSampleReact
!===================================================================================================================================
iSampleReact = ReactionID + 1

SELECT CASE(AnalyzeCase)
CASE(1)
  IF (.NOT.DSMC%ReservoirRateStatistic) THEN
    SurfModel%ProperInfo(SpecID)%NumSurfReact(iSampleReact) = &
        SurfModel%ProperInfo(SpecID)%NumSurfReact(iSampleReact) + Probability
  END IF
  SurfModel%ProperInfo(SpecID)%MeanSurfActE(iSampleReact) = &
      SurfModel%ProperInfo(SpecID)%MeanSurfActE(iSampleReact) + EAct
  SurfModel%ProperInfo(SpecID)%MeanSurfnu(iSampleReact) = &
      SurfModel%ProperInfo(SpecID)%MeanSurfnu(iSampleReact) + nuReact
  SurfModel%ProperInfo(SpecID)%SurfReactCount(iSampleReact) = &
      SurfModel%ProperInfo(SpecID)%SurfReactCount(iSampleReact) + 1
CASE(2)
  IF (DSMC%ReservoirRateStatistic) THEN
    SurfModel%ProperInfo(SpecID)%NumSurfReact(iSampleReact) = &
        SurfModel%ProperInfo(SpecID)%NumSurfReact(iSampleReact) + 1
  END IF
  SurfModel%ProperInfo(SpecID)%ProperSurfActE(iSampleReact) = &
      SurfModel%ProperInfo(SpecID)%ProperSurfActE(iSampleReact) + EAct
  SurfModel%ProperInfo(SpecID)%ProperSurfnu(iSampleReact) = &
      SurfModel%ProperInfo(SpecID)%ProperSurfnu(iSampleReact) + nuReact
  SurfModel%ProperInfo(SpecID)%ProperSurfReactCount(iSampleReact) = &
      SurfModel%ProperInfo(SpecID)%ProperSurfReactCount(iSampleReact) + 1
CASE DEFAULT
  CALL abort(&
__STAMP__,&
'ERROR: analyze case in AnalyzeSurfRates not defined!',AnalyzeCase)
END SELECT

END SUBROUTINE AnalyzeSurfRates
#endif /*(PP_TimeDiscMethod==42)*/
#endif /*(PP_TimeDiscMethod==42) || (PP_TimeDiscMethod==4)*/

#endif /*PARTICLES*/

END MODULE MOD_SurfaceModel_Analyze
