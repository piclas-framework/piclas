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
IF(CalcCollCounter.OR.CalcDesCounter) DoSurfModelAnalyze = .TRUE.

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
INTEGER             :: iSpec
!===================================================================================================================================
IF (.NOT.SurfMesh%SurfOnProc) RETURN
IF (SurfMesh%nOutputSides.EQ.0) RETURN
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
        WRITE(unit_index,'(A)') ''
      END IF
    END IF
#if USE_MPI
  END IF
#endif /*USE_MPI*/

!===================================================================================================================================
! Analyze Routines
!===================================================================================================================================
IF (CalcCollCounter) CALL GetCollCounter(SurfCollNum,AdsorptionNum) !collision counter is reset here
IF (CalcDesCounter) CALL GetDesCounter(DesorptionNum)
! ResetAllInfo
IF (ALLOCATED(SurfModel%Info)) THEN
  DO iSpec = 1,nSpecies
    SurfModel%Info(iSpec)%WallCollCount = 0
    SurfModel%Info(iSpec)%NumOfAds = 0
    SurfModel%Info(iSpec)%NumOfDes = 0
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



#endif /*PARTICLES*/

END MODULE MOD_SurfaceModel_Analyze
