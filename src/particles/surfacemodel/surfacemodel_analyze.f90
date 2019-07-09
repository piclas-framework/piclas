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

#if (PP_TimeDiscMethod==42) || (PP_TimeDiscMethod==4)
INTERFACE WriteDataHeaderInfo
  MODULE PROCEDURE WriteDataHeaderInfo
END INTERFACE

INTERFACE WriteDataInfo
  MODULE PROCEDURE WriteDataInfo
END INTERFACE
#endif /* DSMC*/

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

CALL prms%CreateIntOption(      'Surface-AnalyzeStep'     , 'Analyze is performed each Nth time step for surfaces','1')
CALL prms%CreateLogicalOption(  'Surf-CalcNumSpec'        , 'TODO-DEFINE-PARAMETER\n'//&
                                                            'Calculate the number of simulated'//&
                                                            'particles per species on surfaces','.FALSE.')
CALL prms%CreateLogicalOption(  'Surf-CalcCoverage'       , 'TODO-DEFINE-PARAMETER\n'//&
                                                            'Calculate the surface coverages for'//&
                                                            'each species','.FALSE.')
CALL prms%CreateLogicalOption(  'Surf-CalcAccomodation'   , 'TODO-DEFINE-PARAMETER\n'//&
                                                            'Calculate the surface accomodation coefficient'&
                                                          ,'.FALSE.')
CALL prms%CreateLogicalOption(  'Surf-CalcEvaporation'    , 'TODO-DEFINE-PARAMETER\n'//&
                                                            'Calculate rate of evaporation [kg/s]','.FALSE.')
CALL prms%CreateLogicalOption(  'Surf-CalcAdsorbRates'    , 'TODO-DEFINE-PARAMETER\n'//&
                                                            'Calcualte the adsorption probabilities of species'&
                                                          ,'.FALSE.')
CALL prms%CreateLogicalOption(  'Surf-CalcSurfRates'      , 'TODO-DEFINE-PARAMETER\n'//&
                                                            'Calculate the surface reaction rate per reaction'//&
                                                            ' (k_r)','.FALSE.')

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
SWRITE(UNIT_stdOut,'(A)') ' INIT SURFACE ANALYZE...'

SurfaceAnalyzeStep = GETINT('Surface-AnalyzeStep','1')
IF (SurfaceAnalyzeStep.EQ.0) SurfaceAnalyzeStep = 123456789

DoSurfModelAnalyze = .FALSE.

#if (PP_TimeDiscMethod==42) || (PP_TimeDiscMethod==4)
CalcSurfNumSpec = GETLOGICAL('Surf-CalcNumSpec','.FALSE.')
CalcSurfCoverage = GETLOGICAL('Surf-CalcCoverage','.FALSE.')
#if (PP_TimeDiscMethod==42)
CalcAccomodation = GETLOGICAL('Surf-CalcAccomodation','.FALSE.')
CalcEvaporation = GETLOGICAL('Surf-CalcEvaporation','.FALSE.')
IF (CalcEvaporation) DoSurfModelAnalyze = .TRUE.
CalcAdsorbRates = GETLOGICAL('Surf-CalcAdsorbRates','.FALSE.')
CalcSurfRates = GETLOGICAL('Surf-CalcSurfRates','.FALSE.')
IF(CalcSurfNumSpec.OR.CalcSurfRates.OR.CalcSurfCoverage.OR.CalcAccomodation.OR.Adsorption%TPD.OR.CalcAdsorbRates) &
  DoSurfModelAnalyze = .TRUE.
IF (Adsorption%TPD.AND.((.NOT.CalcSurfRates))) CalcSurfRates = .TRUE.
#else
IF(CalcSurfNumSpec.OR.CalcSurfCoverage.OR.CalcAccomodation) DoSurfModelAnalyze = .TRUE.
#endif
#endif

SurfModelAnalyzeInitIsDone=.TRUE.

SWRITE(UNIT_stdOut,'(A)')' INIT SURFACE ANALYZE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitSurfModelAnalyze


SUBROUTINE AnalyzeSurface(Time)
!===================================================================================================================================
!> create/open SurfaceDatabase.csv and write calculated variables for surface analyze
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars              ,ONLY: BoltzmannConst
USE MOD_Preproc
USE MOD_Analyze_Vars              ,ONLY: DoSurfModelAnalyze
USE MOD_SurfaceModel_Analyze_Vars
USE MOD_Restart_Vars              ,ONLY: DoRestart
#ifdef MPI
USE MOD_Particle_MPI_Vars         ,ONLY: PartMPI
#endif /*MPI*/
#if ( PP_TimeDiscMethod ==42)
USE MOD_Particle_Boundary_Vars    ,ONLY: SurfMesh
USE MOD_SurfaceModel_Vars         ,ONLY: Adsorption
#endif /* DSMC*/
#if ( PP_TimeDiscMethod ==42) || (PP_TimeDiscMethod==4)
USE MOD_Particle_Vars             ,ONLY: nSpecies, PartSurfaceModel
#endif /* DSMC*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: Time
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL             :: isOpen, isRestart
CHARACTER(LEN=350)  :: outfile
INTEGER             :: unit_index, OutputCounter
#if (PP_TimeDiscMethod ==42)
INTEGER             :: iCase
INTEGER             :: iCov, iSpec
CHARACTER(LEN=350)  :: hilf
REAL                :: Adsorptionrate(nSpecies), Desorptionrate(nSpecies), Accomodation(nSpecies)
REAL                :: EvaporationRate(nSpecies)
INTEGER             :: SurfCollNum(nSpecies),AdsorptionNum(nSpecies),DesorptionNum(nSpecies)
REAL,ALLOCATABLE    :: SurfReactRate(:), AdsorptionReactRate(:), AdsorptionActE(:), SurfaceActE(:), ProperSurfaceActE(:)
#endif
#if (PP_TimeDiscMethod ==42) || (PP_TimeDiscMethod ==4)
INTEGER(KIND=8)     :: WallNumSpec(nSpecies), WallNumSpec_SurfDist(nSpecies)
REAL                :: WallCoverage(nSpecies)
#endif
!===================================================================================================================================
  isRestart = .FALSE.
  IF ( DoRestart ) THEN
    isRestart = .TRUE.
  END IF
  IF (.NOT.DoSurfModelAnalyze) RETURN
  OutputCounter = 2
  unit_index = 636
#ifdef MPI
  IF (PartMPI%MPIRoot) THEN
#endif /* MPI */
    INQUIRE(UNIT   = unit_index , OPENED = isOpen)
    IF (.NOT.isOpen) THEN
#if (PP_TimeDiscMethod==42)
    ! if only the reaction rate is desired (resevoir) the initial temperature
    ! of the second species is added to the filename
      IF (Adsorption%TPD) THEN
          iCov = INT(Adsorption%Coverage(1,1,1,1)*1000)
          WRITE( hilf, '(I4.4)') iCov
        outfile = 'SurfaceDatabase_Cov_'//TRIM(hilf)//'.csv'
      ELSE
        outfile = 'SurfaceDatabase.csv'
      END IF
#else
      outfile = 'SurfaceDatabase.csv'
#endif

!===================================================================================================================================
! Write Header
!===================================================================================================================================
      IF (isRestart .and. FILEEXISTS(outfile)) THEN
        OPEN(unit_index,file=TRIM(outfile),position="APPEND",status="OLD")
      ELSE
        OPEN(unit_index,file=TRIM(outfile))
        !--- insert header
        WRITE(unit_index,'(A8)',ADVANCE='NO') '001-TIME'
#if (PP_TimeDiscMethod==42) || (PP_TimeDiscMethod==4)
        IF (PartSurfaceModel.EQ.3) THEN
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
          IF (CalcAdsorbRates) THEN
            CALL WriteDataHeaderInfo(unit_index,'nSurfColl-Spec',OutputCounter,nSpecies)
            CALL WriteDataHeaderInfo(unit_index,'N_Ads-Spec',OutputCounter,nSpecies)
            CALL WriteDataHeaderInfo(unit_index,'Prob_adsorption-Spec',OutputCounter,nSpecies)
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-P_Molec-Adsorb-Spec-',iSpec,' '
                OutputCounter = OutputCounter + 1
              DO iCase = 1,Adsorption%DissNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3,A5)',ADVANCE='NO') &
                    OutputCounter,'-P_Dissoc-Spec-',iSpec,'-Reaction-', iCase,' '
                OutputCounter = OutputCounter + 1
              END DO
              DO iCase = 1, Adsorption%RecombNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3,A5)',ADVANCE='NO') &
                    OutputCounter,'-P_ER-Spec-',iSpec,'-Reaction-', iCase,' '
                OutputCounter = OutputCounter + 1
              END DO
            END DO
            DO iSpec = 1, nSpecies
              DO iCase = 1,Adsorption%DissNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3,A5)',ADVANCE='NO') &
                    OutputCounter,'-E-diss-Spec-', iSpec,'-Reaction-', iCase,' '
                OutputCounter = OutputCounter + 1
              END DO
              DO iCase = 1,Adsorption%RecombNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3,A5)',ADVANCE='NO') &
                    OutputCounter,'-E-ER-Spec-', iSpec,'-Reaction-', iCase,' '
                OutputCounter = OutputCounter + 1
              END DO
            END DO
          END IF
          IF (CalcSurfRates) THEN
            CALL WriteDataHeaderInfo(unit_index,'N_Des-Spec',OutputCounter,nSpecies)
            CALL WriteDataHeaderInfo(unit_index,'P_Des-Spec',OutputCounter,nSpecies)
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-P-SurfDesorb-Molec-Spec-', iSpec,' '
              OutputCounter = OutputCounter + 1
              DO iCase = 1, Adsorption%DissNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3,A5)',ADVANCE='NO') &
                    OutputCounter,'-P-SurfDissoc-Spec-',iSpec,'-Reaction-', iCase,' '
                OutputCounter = OutputCounter + 1
              END DO
              DO iCase = 1, Adsorption%RecombNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3,A5)',ADVANCE='NO') &
                    OutputCounter,'-P-SurfLH-Spec-',iSpec,'-Reaction-', iCase,' '
                OutputCounter = OutputCounter + 1
              END DO
            END DO
            CALL WriteDataHeaderInfo(unit_index,'P-Surfexch-Case',OutputCounter,Adsorption%NumOfExchReact)
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-E-Desorb-Spec-', iSpec,' '
              OutputCounter = OutputCounter + 1
              DO iCase = 1,Adsorption%DissNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3,A5)',ADVANCE='NO') &
                    OutputCounter,'-E-Diss-Spec-', iSpec,'-Reaction-', iCase,' '
                OutputCounter = OutputCounter + 1
              END DO
              DO iCase = 1,Adsorption%RecombNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3,A5)',ADVANCE='NO') &
                    OutputCounter,'-E-LH-Spec-', iSpec,'-Reaction-', iCase,' '
                OutputCounter = OutputCounter + 1
              END DO
            END DO
            CALL WriteDataHeaderInfo(unit_index,'E-Exch-Reaction',OutputCounter,Adsorption%NumOfExchReact)
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-Proper-E-Desorb-Spec-', iSpec,' '
              OutputCounter = OutputCounter + 1
              DO iCase = 1,Adsorption%DissNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3,A5)',ADVANCE='NO') &
                    OutputCounter,'-Proper-E-Diss-Spec-', iSpec,'-Reaction-', iCase,' '
                OutputCounter = OutputCounter + 1
              END DO
              DO iCase = 1,Adsorption%RecombNum
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3,A5)',ADVANCE='NO') &
                    OutputCounter,'-Proper-E-LH-Spec-', iSpec,'-Reaction-', iCase,' '
                OutputCounter = OutputCounter + 1
              END DO
            END DO
            CALL WriteDataHeaderInfo(unit_index,'Proper-E-Exch-Reaction',OutputCounter,Adsorption%NumOfExchReact)
          END IF
          IF (Adsorption%TPD) THEN
            CALL WriteDataHeaderInfo(unit_index,'WallTemp',OutputCounter,1)
          END IF
        END IF
        IF (CalcEvaporation) THEN
          CALL WriteDataHeaderInfo(unit_index,'Evap-Mass-Spec',OutputCounter,nSpecies)
#endif
        END IF
#endif
        WRITE(unit_index,'(A1)') ' '
      END IF
    END IF
#ifdef MPI
  END IF
#endif /* MPI */

!===================================================================================================================================
! Analyze Routines
!===================================================================================================================================
#if (PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==42)
IF (PartSurfaceModel.EQ.3) THEN
  IF (CalcSurfNumSpec.OR.CalcSurfCoverage) CALL GetWallNumSpec(WallNumSpec,WallCoverage,WallNumSpec_SurfDist)
#if (PP_TimeDiscMethod==42)
  IF (CalcAccomodation) CALL GetAccCoeff(Accomodation)
  IF (CalcAdsorbRates) THEN
    SDEALLOCATE(AdsorptionReactRate)
    SDEALLOCATE(AdsorptionActE)
    ALLOCATE(AdsorptionReactRate(1:nSpecies*(Adsorption%ReactNum+1)))
    ALLOCATE(AdsorptionActE(1:nSpecies*(Adsorption%ReactNum)))
    CALL GetAdsRates(Adsorptionrate,SurfCollNum,AdsorptionNum,AdsorptionReactRate,AdsorptionActE)
  ELSE
    IF(SurfMesh%SurfOnProc)THEN
      DO iSpec = 1,nSpecies
        Adsorption%AdsorpInfo(iSpec)%WallCollCount = 0
      END DO
    END IF
  END IF
  IF (CalcSurfRates) THEN
    SDEALLOCATE(SurfReactRate)
    SDEALLOCATE(SurfaceActE)
    ALLOCATE(SurfReactRate(1:nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact))
    ALLOCATE(SurfaceActE(1:nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact))
    ALLOCATE(ProperSurfaceActE(1:nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact))
    CALL GetSurfRates(Desorptionrate,DesorptionNum,SurfReactRate,SurfaceActE,ProperSurfaceActE)
  END IF
#endif
END IF
#endif
#if (PP_TimeDiscMethod==42)
IF (CalcEvaporation) CALL GetEvaporationRate(EvaporationRate)
#endif /*PP_TimeDiscMethod==42*/
!===================================================================================================================================
! Output Analyzed variables
!===================================================================================================================================
#ifdef MPI
IF (PartMPI%MPIROOT) THEN
#endif    /* MPI */
  WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') Time
#if ((PP_TimeDiscMethod==42) || (PP_TimeDiscMethod==4))
! output for adsorption
    IF (PartSurfaceModel.EQ.3) THEN
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
      IF (CalcAdsorbRates) THEN
        CALL WriteDataInfo(unit_index,nSpecies                         ,IntegerArray=SurfCollNum(:))
        CALL WriteDataInfo(unit_index,nSpecies                         ,IntegerArray=AdsorptionNum(:))
        CALL WriteDataInfo(unit_index,nSpecies                         ,RealArray=Adsorptionrate(:))
        CALL WriteDataInfo(unit_index,nSpecies*(Adsorption%ReactNum+1) ,RealArray=AdsorptionReactRate(:))
        CALL WriteDataInfo(unit_index,nSpecies*Adsorption%ReactNum     ,RealArray=AdsorptionActE(:))
      END IF
      IF (CalcSurfRates) THEN
        CALL WriteDataInfo(unit_index,nSpecies                                                   ,IntegerArray=DesorptionNum(:))
        CALL WriteDataInfo(unit_index,nSpecies                                                   ,RealArray=Desorptionrate(:))
        CALL WriteDataInfo(unit_index,nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact ,RealArray=SurfReactRate(:))
        CALL WriteDataInfo(unit_index,nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact ,RealArray=SurfaceActE(:))
        CALL WriteDataInfo(unit_index,nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact ,RealArray=ProperSurfaceActE(:))
      END IF
      IF (Adsorption%TPD) THEN
        CALL WriteDataInfo(unit_index,1,RealScalar=Adsorption%TPD_Temp)
      END IF
    END IF
    IF (CalcEvaporation) THEN
      CALL WriteDataInfo(unit_index,nSpecies,RealArray=EvaporationRate(:))
#endif /*(PP_TimeDiscMethod==42)*/
    END IF
#endif /*(PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==42)*/
    WRITE(unit_index,'(A1)') ' '
#ifdef MPI
  END IF
#endif /* MPI */
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE AnalyzeSurface


#if (PP_TimeDiscMethod==42) || (PP_TimeDiscMethod==4)
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
  WRITE(unit_index,'(I3.3,A,A,A,I3.3,A3)',ADVANCE='NO') OutputCounter,'-',AttribName,'-',iLoop,'   '
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
CHARACTER(LEN=*)  ,INTENT(IN),OPTIONAL :: StrScalar(1)
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
    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
    WRITE(unit_index,'(ES25.14E3)',ADVANCE='NO') RealArray(iLoop)
  END DO
END IF
IF(PRESENT(RealScalar)) THEN
  WRITE(unit_index,'(A1)',ADVANCE='NO') ','
  WRITE(unit_index,'(ES25.14E3)',ADVANCE='NO') RealScalar
END IF

IF(PRESENT(IntegerArray)) THEN
  DO iLoop = 1, nVal
    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
    WRITE(unit_index,'(I18.1)',ADVANCE='NO') IntegerArray(iLoop)
  END DO
END IF

IF(PRESENT(IntegerK8Array)) THEN
  DO iLoop = 1, nVal
    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
    WRITE(unit_index,'(I18.1)',ADVANCE='NO') IntegerK8Array(iLoop)
  END DO
END IF

IF(PRESENT(IntegerScalar)) THEN
  WRITE(unit_index,'(A1)',ADVANCE='NO') ','
  WRITE(unit_index,'(I18.1)',ADVANCE='NO') IntegerScalar
END IF

IF(PRESENT(StrArray)) THEN
  DO iLoop = 1, nVal
    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
    WRITE(unit_index,'(A)',ADVANCE='NO') StrArray(iLoop)
  END DO
END IF

IF(PRESENT(StrScalar)) THEN
  WRITE(unit_index,'(A1)',ADVANCE='NO') ','
  WRITE(unit_index,'(A)',ADVANCE='NO') StrScalar
END IF

IF(PRESENT(LogicalScalar)) THEN
  WRITE(unit_index,'(A1)',ADVANCE='NO') ','
  WRITE(unit_index,'(L2)',ADVANCE='NO') LogicalScalar
END IF
END SUBROUTINE WriteDataInfo
#endif /*DSMC*/


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
#ifdef MPI
USE MOD_Particle_Boundary_Vars    ,ONLY: SurfCOMM
USE MOD_Particle_MPI_Vars         ,ONLY: PartMPI
#endif /*MPI*/
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
#ifdef MPI
REAL                            :: RD(nSpecies)
INTEGER(KIND=8)                 :: IDR(nSpecies), ID1(nSpecies), ID2(nSpecies), ID3(nSpecies*2)
#endif /*MPI*/
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

IF(SurfMesh%SurfOnProc)THEN
  DO iSpec=1,nSpecies
  DO iSurfSide=1,SurfMesh%nSides
    SideID = Adsorption%SurfSideToGlobSideMap(iSurfSide)
    PartboundID = PartBound%MapToPartBC(BC(SideID))
    IF (PartBound%SolidReactive(PartboundID)) THEN
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
  IF (CalcSurfCoverage) THEN
    WallCoverage(:) = Coverage(:) / (SurfMesh%nSides*nSurfSample*nSurfSample)
  END IF
END IF

#ifdef MPI
  IF (PartMPI%MPIRoot) THEN
    IF (CalcSurfNumSpec)  THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,SubWallNumSpec      ,nSpecies  ,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(MPI_IN_PLACE,WallNumSpec_SurfDist,nSpecies  ,MPI_LONG,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(MPI_IN_PLACE,WallNumSpec_tmp     ,nSpecies*2,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
    END IF
    IF (CalcSurfCoverage) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,WallCoverage,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
      WallCoverage = WallCoverage / REAL(SurfCOMM%nProcs)
    END IF
  ELSE
    IF (CalcSurfNumSpec) THEN
      CALL MPI_REDUCE(SubWallNumSpec      ,ID1,nSpecies  ,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(WallNumSpec_SurfDist,ID2,nSpecies  ,MPI_LONG,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(WallNumSpec_tmp     ,ID3,nSpecies*2,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
    END IF
    IF (CalcSurfCoverage) CALL MPI_REDUCE(WallCoverage,RD,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  END IF
#endif /*MPI*/

  IF (KeepWallParticles.AND.CalcSurfNumSpec) THEN
    DO i=1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(i) .AND. PDM%ParticleAtWall(i)) THEN
        WallNumSpec(PartSpecies(i)) = WallNumSpec(PartSpecies(i)) + 1
      END IF
    END DO
#ifdef MPI
  IF (PartMPI%MPIRoot) THEN
    IF (CalcSurfNumSpec) CALL MPI_REDUCE(MPI_IN_PLACE,WallNumSpec,nSpecies,MPI_LONG,MPI_SUM,0,PartMPI%COMM,IERROR)
  ELSE
    IF (CalcSurfNumSpec) CALL MPI_REDUCE(WallNumSpec ,IDR        ,nSpecies,MPI_LONG,MPI_SUM,0,PartMPI%COMM,IERROR)
  END IF
#endif /*MPI*/
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
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption
USE MOD_Particle_Boundary_Vars ,ONLY: SurfMesh
#ifdef MPI
USE MOD_Particle_Boundary_Vars ,ONLY: SurfCOMM
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*MPI*/
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
#ifdef MPI
REAL                            :: AC(nSpecies)
#endif /*MPI*/
!===================================================================================================================================

Accomodation(:) = 0.
IF(SurfMesh%SurfOnProc)THEN
  IF (DSMC%ReservoirRateStatistic) THEN
    DO iSpec = 1,nSpecies
      IF (Adsorption%AdsorpInfo(iSpec)%WallCollCount.GT.0) THEN
        Accomodation(iSpec) = Adsorption%AdsorpInfo(iSpec)%Accomodation / REAL(Adsorption%AdsorpInfo(iSpec)%WallCollCount)
      ELSE
        Accomodation(iSpec) = 0.
      END IF
    END DO
  ELSE IF (.NOT.DSMC%ReservoirRateStatistic) THEN
    DO iSpec = 1,nSpecies
      IF (Adsorption%AdsorpInfo(iSpec)%WallCollCount.GT.0) THEN
        Accomodation(iSpec) = Adsorption%AdsorpInfo(iSpec)%Accomodation / REAL(Adsorption%AdsorpInfo(iSpec)%WallCollCount)
      ELSE
        Accomodation(iSpec) = 0.
      END IF
    END DO
  END IF
END IF

#ifdef MPI
IF (PartMPI%MPIRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,Accomodation,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  Accomodation= Accomodation/ REAL(SurfCOMM%nProcs)
ELSE
  CALL MPI_REDUCE(Accomodation,AC          ,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
END IF
#endif /*MPI*/

IF(SurfMesh%SurfOnProc)THEN
  DO iSpec = 1,nSpecies
    Adsorption%AdsorpInfo(iSpec)%Accomodation = 0.
  END DO
END IF

END SUBROUTINE GetAccCoeff


SUBROUTINE GetAdsRates(AdsorbRate,SurfCollNum,AdsorbNum,ReactRate,AdsorbActE)
!===================================================================================================================================
! Calculate adsorption, desorption and accomodation rates for all species
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars          ,ONLY: nSpecies, PartSurfaceModel
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh
#ifdef MPI
USE MOD_Particle_Boundary_Vars ,ONLY: SurfCOMM
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL   , INTENT(OUT)            :: AdsorbRate(nSpecies)
REAL   , INTENT(OUT)            :: ReactRate(nSpecies*(Adsorption%ReactNum+1))
REAL   , INTENT(OUT)            :: AdsorbActE(nSpecies*Adsorption%ReactNum)
INTEGER, INTENT(OUT)            :: SurfCollNum(nSpecies), AdsorbNum(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSpec, iCase, iReact
#ifdef MPI
REAL                            :: AD(nSpecies),RR(nSpecies*Adsorption%ReactNum)
INTEGER                         :: ADN(nSpecies)
#endif /*MPI*/
!===================================================================================================================================

IF(SurfMesh%SurfOnProc)THEN
  IF (DSMC%ReservoirRateStatistic) THEN
    DO iSpec = 1,nSpecies
      IF (Adsorption%AdsorpInfo(iSpec)%WallCollCount.GT.0) THEN
        AdsorbRate(iSpec) = REAL(Adsorption%AdsorpInfo(iSpec)%NumOfAds) / REAL(Adsorption%AdsorpInfo(iSpec)%WallCollCount)
      ELSE
        AdsorbRate(iSpec) = 0.
      END IF
    END DO
  ELSE IF (.NOT.DSMC%ReservoirRateStatistic) THEN
    DO iSpec = 1,nSpecies
      IF (Adsorption%AdsorpInfo(iSpec)%WallCollCount.GT.0) THEN
        IF (PartSurfaceModel.EQ.1) THEN
          AdsorbRate(iSpec) = Adsorption%AdsorpInfo(iSpec)%MeanProbAds / REAL(nSurfSample * nSurfSample * SurfMesh%nSides)
        ELSE IF (PartSurfaceModel.EQ.3) THEN
          AdsorbRate(iSpec) = Adsorption%AdsorpInfo(iSpec)%MeanProbAds / REAL(Adsorption%AdsorpInfo(iSpec)%WallCollCount)
        END IF
      ELSE
        AdsorbRate(iSpec)= 0.
      END IF
    END DO
  END IF

  iCase = 1
  DO iSpec = 1, nSpecies
    DO iReact = 1, Adsorption%ReactNum+1
      IF (Adsorption%AdsorpReactInfo(iSpec)%AdsReactCount(iReact).GT.0) THEN
        ReactRate(iCase) = Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(iReact) &
            / REAL(Adsorption%AdsorpReactInfo(iSpec)%AdsReactCount(iReact)) !* REAL(Adsorption%AdsorpInfo(iSpec)%WallCollCount)
      ELSE
        ReactRate(iCase) = 0.
      END IF
      iCase = iCase + 1
    END DO
  END DO

  iCase = 1
  DO iSpec = 1,nSpecies
    DO iReact = 1,Adsorption%ReactNum
      IF (Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(iReact).GT.0) THEN
        IF (PartSurfaceModel.EQ.1) THEN
          AdsorbActE(iCase) = Adsorption%AdsorpReactInfo(iSpec)%MeanAdsActE(iReact) &
              / REAL(nSurfSample * nSurfSample * SurfMesh%nSides)
        ELSE IF (PartSurfaceModel.EQ.3) THEN
          IF (Adsorption%AdsorpReactInfo(iSpec)%AdsReactCount(iReact).GT.0) THEN
            AdsorbActE(iCase) = Adsorption%AdsorpReactInfo(iSpec)%MeanAdsActE(iReact) &
                / REAL(Adsorption%AdsorpReactInfo(iSpec)%AdsReactCount(iReact))
          END IF
        END IF
      ELSE
        AdsorbActE(iCase)= 0.
      END IF
      iCase = iCase + 1
    END DO
  END DO
  DO iSpec = 1,nSpecies
    SurfCollNum(iSpec) = Adsorption%AdsorpInfo(iSpec)%WallCollCount
    AdsorbNum(iSpec) = Adsorption%AdsorpInfo(iSpec)%NumOfAds
  END DO
ELSE
  SurfCollNum(:) = 0
  AdsorbRate(:) = 0.
  AdsorbNum(:) = 0
  ReactRate(:) = 0.
  AdsorbActE(:) = 0.
END IF

#ifdef MPI
IF (PartMPI%MPIRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,AdsorbRate  ,nSpecies                    ,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,SurfCollNum ,nSpecies                    ,MPI_LONG            ,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,AdsorbNum   ,nSpecies                    ,MPI_LONG            ,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,ReactRate   ,nSpecies*(Adsorption%ReactNum+1),MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,AdsorbActE  ,nSpecies*Adsorption%ReactNum,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  AdsorbRate = AdsorbRate  / REAL(SurfCOMM%nProcs)
  SurfCollNum= INT( REAL(SurfCollNum) / REAL(SurfCOMM%nProcs) )
  AdsorbNum  = INT( REAL(AdsorbNum)   / REAL(SurfCOMM%nProcs) )
  ReactRate  = ReactRate   / REAL(SurfCOMM%nProcs)
  AdsorbActE = AdsorbActE  / REAL(SurfCOMM%nProcs)
ELSE
  CALL MPI_REDUCE(AdsorbRate  ,AD          ,nSpecies                    ,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(SurfCollNum ,ADN         ,nSpecies                    ,MPI_LONG            ,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(AdsorbNum   ,ADN         ,nSpecies                    ,MPI_LONG            ,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(ReactRate   ,RR          ,nSpecies*(Adsorption%ReactNum+1),MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(AdsorbActE  ,RR          ,nSpecies*Adsorption%ReactNum,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
END IF
#endif /*MPI*/

IF(SurfMesh%SurfOnProc)THEN
  DO iSpec = 1,nSpecies
    Adsorption%AdsorpInfo(iSpec)%WallCollCount = 0
    Adsorption%AdsorpInfo(iSpec)%MeanProbAds = 0.
    Adsorption%AdsorpInfo(iSpec)%NumOfAds = 0
    Adsorption%AdsorpReactInfo(iSpec)%NumAdsReact(:) = 0.
    Adsorption%AdsorpReactInfo(iSpec)%AdsReactCount(:) = 0
    Adsorption%AdsorpReactInfo(iSpec)%MeanAdsActE(:) = 0.
  END DO
END IF

END SUBROUTINE GetAdsRates


SUBROUTINE GetSurfRates(DesorbRate,DesorbNum,ReactRate,SurfaceActE,ProperSurfaceActE)
!===================================================================================================================================
! Calculate adsorption, desorption and accomodation rates for all species
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars          ,ONLY: nSpecies, PartSurfaceModel
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh
#ifdef MPI
USE MOD_Particle_Boundary_Vars ,ONLY: SurfCOMM
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL   , INTENT(OUT)            :: DesorbRate(nSpecies)
INTEGER, INTENT(OUT)            :: DesorbNum(nSpecies)
REAL   , INTENT(OUT)            :: ReactRate(nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact)
REAL   , INTENT(OUT)            :: SurfaceActE(nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact)
REAL   , INTENT(OUT)            :: ProperSurfaceActE(nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSpec, iReact, iCase
#ifdef MPI
INTEGER                         :: commSize
REAL                            :: DE(nSpecies)
REAL                            :: RR(nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact)
INTEGER                         :: DEN(nSpecies)
#endif /*MPI*/
!===================================================================================================================================

IF(SurfMesh%SurfOnProc)THEN
  IF (DSMC%ReservoirRateStatistic) THEN
    DO iSpec = 1,nSpecies
      IF (Adsorption%AdsorpInfo(iSpec)%WallSpecNumCount.GT.0) THEN
        DesorbRate(iSpec) = REAL(Adsorption%AdsorpInfo(iSpec)%NumOfDes) / REAL(Adsorption%AdsorpInfo(iSpec)%WallSpecNumCount)
      ELSE
        DesorbRate(iSpec) = 0.
      END IF
    END DO
  ELSE IF (.NOT.DSMC%ReservoirRateStatistic) THEN
    iCase = 1
    DO iSpec = 1,nSpecies
      DesorbRate(iSpec)= 0.
      IF (PartSurfaceModel.EQ.1) THEN
        DO iReact = 1, Adsorption%ReactNum+1
          ReactRate(iCase) = Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(iReact) &
              / REAL(nSurfSample * nSurfSample * SurfMesh%nSides)
          iCase = iCase + 1
        END DO
      ELSE IF (PartSurfaceModel.EQ.3) THEN
        DO iReact = 1, Adsorption%ReactNum+1
          IF (Adsorption%AdsorpReactInfo(iSpec)%SurfReactCount(iReact).GT.0) THEN
            ReactRate(iCase) = Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact(iReact) &
                / REAL(Adsorption%AdsorpReactInfo(iSpec)%SurfReactCount(iReact))
          ELSE
            ReactRate(iCase) = 0.
          END IF
          iCase = iCase + 1
        END DO
      END IF
    END DO
  END IF
  DO iSpec = 1,nSpecies
    DesorbNum(iSpec) = Adsorption%AdsorpInfo(iSpec)%NumOfDes
  END DO
ELSE
  DesorbNum(:)  = 0
  DesorbRate(:) = 0.
  ReactRate(:)  = 0.
  SurfaceActE(:)= 0.
  ProperSurfaceActE(:)= 0.
END IF

#ifdef MPI
IF (PartMPI%MPIRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,DesorbRate  ,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,DesorbNum   ,nSpecies,MPI_LONG            ,MPI_SUM,0,PartMPI%COMM,IERROR)
  DesorbRate  = DesorbRate / REAL(SurfCOMM%nProcs)
  DesorbNum   = INT( REAL(DesorbNum) / REAL(SurfCOMM%nProcs) )
ELSE
  CALL MPI_REDUCE(DesorbRate  ,DE          ,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(DesorbNum   ,DEN         ,nSpecies,MPI_LONG            ,MPI_SUM,0,PartMPI%COMM,IERROR)
END IF
#endif /*MPI*/

IF(SurfMesh%SurfOnProc)THEN
  DO iSpec = 1,nSpecies
    Adsorption%AdsorpInfo(iSpec)%MeanProbDes = 0.
    Adsorption%AdsorpInfo(iSpec)%NumOfDes = 0
  END DO
END IF

IF(SurfMesh%SurfOnProc)THEN
  iCase = 1
  DO iSpec = 1,nSpecies
    DO iReact = 1,Adsorption%ReactNum+1
      IF (Adsorption%AdsorpReactInfo(iSpec)%SurfReactCount(iReact).GT.0) THEN
        IF (PartSurfaceModel.EQ.1) THEN
          SurfaceActE(iCase) = Adsorption%AdsorpReactInfo(iSpec)%MeanSurfActE(iReact) &
              / REAL(nSurfSample * nSurfSample * SurfMesh%nSides)
        ELSE IF (PartSurfaceModel.EQ.3) THEN
          SurfaceActE(iCase) = Adsorption%AdsorpReactInfo(iSpec)%MeanSurfActE(iReact) &
              / REAL(Adsorption%AdsorpReactInfo(iSpec)%SurfReactCount(iReact))
        END IF
      ELSE
        SurfaceActE(iCase) = 0.
      END IF
      iCase = iCase + 1
    END DO
  END DO
ELSE
  SurfaceActE(:)= 0.
END IF

IF(SurfMesh%SurfOnProc)THEN
  iCase = 1
  DO iSpec = 1,nSpecies
    DO iReact = 1,Adsorption%ReactNum+1
      IF (Adsorption%AdsorpReactInfo(iSpec)%ProperSurfReactCount(iReact).GT.0) THEN
        IF (PartSurfaceModel.EQ.3) THEN
          ProperSurfaceActE(iCase) = Adsorption%AdsorpReactInfo(iSpec)%ProperSurfActE(iReact) &
              / REAL(Adsorption%AdsorpReactInfo(iSpec)%ProperSurfReactCount(iReact))
        END IF
      ELSE
        ProperSurfaceActE(iCase) = 0.
      END IF
      iCase = iCase + 1
    END DO
  END DO
ELSE
  ProperSurfaceActE(:)= 0.
END IF

#ifdef MPI
commSize = nSpecies*(Adsorption%ReactNum+1)+Adsorption%NumOfExchReact
IF (PartMPI%MPIRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE ,ReactRate        ,commSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE ,SurfaceActE      ,commSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE ,ProperSurfaceActE,commSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  ReactRate   = ReactRate   / REAL(SurfCOMM%nProcs)
  SurfaceActE = SurfaceActE / REAL(SurfCOMM%nProcs)
  ProperSurfaceActE = ProperSurfaceActE / REAL(SurfCOMM%nProcs)
ELSE
  CALL MPI_REDUCE(ReactRate         ,RR ,commSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(SurfaceActE       ,RR ,commSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(ProperSurfaceActE ,RR ,commSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
END IF
#endif /*MPI*/

IF(SurfMesh%SurfOnProc)THEN
  DO iSpec = 1,nSpecies
    Adsorption%AdsorpInfo(iSpec)%MeanProbDes = 0.
    Adsorption%AdsorpInfo(iSPec)%NumOfDes = 0
    Adsorption%AdsorpReactInfo(iSpec)%MeanSurfActE = 0.
    Adsorption%AdsorpReactInfo(iSpec)%ProperSurfActE = 0.
    Adsorption%AdsorpReactInfo(iSpec)%NumSurfReact = 0.
    Adsorption%AdsorpReactInfo(iSpec)%SurfReactCount = 0
    Adsorption%AdsorpReactInfo(iSpec)%ProperSurfReactCount = 0
  END DO
END IF

END SUBROUTINE GetSurfRates


SUBROUTINE GetEvaporationRate(EvaporationRate)
!===================================================================================================================================
! Calculate evaporation rate from number of particles of a species evaporating from surface in the defined analyze time [kg/s]
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars         ,ONLY: Species, nSpecies
USE MOD_Particle_Analyze_Vars
USE MOD_SurfaceModel_Vars     ,ONLY: Liquid
#ifdef MPI
!USE MOD_Particle_Boundary_Vars, ONLY : SurfCOMM
USE MOD_Particle_MPI_Vars     ,ONLY: PartMPI
#endif /*MPI*/
USE MOD_TimeDisc_Vars         ,ONLY: dt
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)               :: EvaporationRate(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSpec
#ifdef MPI
REAL                            :: RD(nSpecies)
#endif /*MPI*/
!===================================================================================================================================
EvaporationRate = 0.

DO iSpec=1,nSpecies
  EvaporationRate(iSpec) = Species(iSpec)%MassIC * Species(iSpec)%MacroParticleFactor &
                        * REAL(Liquid%Info(iSpec)%NumOfDes - Liquid%Info(iSpec)%NumOfAds) / dt
END DO

Liquid%Info(:)%NumOfAds = 0
Liquid%Info(:)%NumOfDes = 0

#ifdef MPI
  IF (PartMPI%MPIRoot) THEN
    CALL MPI_REDUCE(MPI_IN_PLACE,EvaporationRate,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  ELSE
    CALL MPI_REDUCE(EvaporationRate,RD,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  END IF
#endif /*MPI*/

END SUBROUTINE GetEvaporationRate
#endif /*(PP_TimeDiscMethod==42)*/
#endif /*(PP_TimeDiscMethod==42) || (PP_TimeDiscMethod==4)*/

#endif /*PARTICLES*/

END MODULE MOD_SurfaceModel_Analyze
