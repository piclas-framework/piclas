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
                                                          'adsorbed/desorbed particles per species','.FALSE.')
CALL prms%CreateLogicalOption(  'Surf-CalcPorousBCInfo' , 'Calculate output of porous BCs such pumping speed, removal '//&
                                                             'probability and pressure (normalized with the given pressure). '//&
                                                             'Values are averaged over the whole porous BC.' , '.FALSE.')
!-- BoundaryParticleOutput
CALL prms%CreateLogicalOption(  'CalcBoundaryParticleOutput', 'Count number of particles exiting for species X on boundary X' , '.FALSE.')
CALL prms%CreateIntOption(      'BPO-NPartBoundaries'       , 'Number of boundaries used for CalcBoundaryParticleOutput')
CALL prms%CreateIntArrayOption( 'BPO-PartBoundaries'        , 'Vector (length BPO-NPartBoundaries) with the numbers of each Part-Boundary')
CALL prms%CreateIntOption(      'BPO-NSpecies'              , 'Number of species used for CalcBoundaryParticleOutput')
CALL prms%CreateIntArrayOption( 'BPO-Species'               , 'Vector (length BPO-NSpecies) with the corresponding Species IDs')

END SUBROUTINE DefineParametersSurfModelAnalyze


SUBROUTINE InitSurfModelAnalyze()
!===================================================================================================================================
!> Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_ReadInTools               ,ONLY: GETLOGICAL,GETINT,GETINTARRAY
USE MOD_Particle_Vars             ,ONLY: nSpecies
USE MOD_Analyze_Vars              ,ONLY: DoSurfModelAnalyze
USE MOD_SurfaceModel_Vars         ,ONLY: nPorousBC
USE MOD_SurfaceModel_Analyze_Vars
USE MOD_Particle_Boundary_Vars    ,ONLY: nPartBound,PartBound
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iPartBound,iSpec
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

!-- Surface Collision Counter
CalcSurfCollCounter = GETLOGICAL('Surf-CalcCollCounter')
IF(CalcSurfCollCounter) THEN
  DoSurfModelAnalyze = .TRUE.
  ! allocate info and constants
  ALLOCATE(SurfAnalyzeCount(1:nSpecies),SurfAnalyzeNumOfAds(1:nSpecies),SurfAnalyzeNumOfDes(1:nSpecies))
  SurfAnalyzeCount = 0; SurfAnalyzeNumOfAds = 0;  SurfAnalyzeNumOfDes = 0
END IF

!-- Porous Boundaries
IF(nPorousBC.GT.0) THEN
  ! Output for porous BC: Pump averaged values
  CalcPorousBCInfo = GETLOGICAL('Surf-CalcPorousBCInfo')
  IF(CalcPorousBCInfo) THEN
    DoSurfModelAnalyze = .TRUE.
    ALLOCATE(PorousBCOutput(1:5,1:nPorousBC))
    PorousBCOutput = 0.
  END IF
END IF

!-- BoundaryParticleOutput (after mapping of PartBound on FieldBound and determination of PartBound types = open, reflective etc.)
CalcBoundaryParticleOutput = GETLOGICAL('CalcBoundaryParticleOutput')
IF(CalcBoundaryParticleOutput)THEN
  DoSurfModelAnalyze = .TRUE.
  BPO%NPartBoundaries = GETINT('BPO-NPartBoundaries')
  BPO%PartBoundaries  = GETINTARRAY('BPO-PartBoundaries',BPO%NPartBoundaries)
  BPO%NSpecies        = GETINT('BPO-NSpecies')
  BPO%Species         = GETINTARRAY('BPO-Species',BPO%NSpecies)
  IF(BPO%NPartBoundaries.EQ.0.OR.BPO%NSpecies.EQ.0)THEN
    CALL abort(&
    __STAMP__&
    ,'BPO-NPartBoundaries or BPO-NSpecies is zero, which is not allowed')
  END IF ! BPO%NPartBoundaries.EQ.0.OR.BPO%NSpecies.EQ.0
  ALLOCATE(BPO%RealPartOut(1:BPO%NPartBoundaries,1:BPO%NSpecies))
  BPO%RealPartOut = 0.

  ALLOCATE(BPO%SpecIDToBPOSpecID(1:nSPecies))
  BPO%SpecIDToBPOSpecID = -1
  DO iSpec = 1, BPO%NSpecies
    BPO%SpecIDToBPOSpecID(BPO%Species(iSpec)) = iSpec
  END DO ! iSpec = 1, BPO%NSpecies

  ALLOCATE(BPO%BCIDToBPOBCID(1:nPartBound))
  BPO%BCIDToBPOBCID     = -1
  DO iPartBound = 1, BPO%NPartBoundaries
    BPO%BCIDToBPOBCID(BPO%PartBoundaries(iPartBound)) = iPartBound
    ! Sanity check BC types: BPO%PartBoundaries(iPartBound) = 1 (open)
    ! Add more BCs to the vector if required
    IF(.NOT.ANY(PartBound%TargetBoundCond(BPO%PartBoundaries(iPartBound)).EQ.(/1/)))THEN
      SWRITE(UNIT_stdOut,'(A)')'\nError for CalcBoundaryParticleOutput=T\n'
      SWRITE(UNIT_stdOut,'(A,I0)')'  iPartBound = ',BPO%PartBoundaries(iPartBound)
      SWRITE(UNIT_stdOut,'(A,A)') '  SourceName = ',TRIM(PartBound%SourceBoundName(BPO%PartBoundaries(iPartBound)))
      SWRITE(UNIT_stdOut,'(A,I0)')'   Condition = ',PartBound%TargetBoundCond(BPO%PartBoundaries(iPartBound))
      SWRITE(UNIT_stdOut,'(A)')'\n  Conditions are'//&
          '  OpenBC          = 1  \n'//&
          '                  ReflectiveBC    = 2  \n'//&
          '                  PeriodicBC      = 3  \n'//&
          '                  SimpleAnodeBC   = 4  \n'//&
          '                  SimpleCathodeBC = 5  \n'//&
          '                  RotPeriodicBC   = 6  \n'//&
          '                  SymmetryBC      = 10 \n'//&
          '                  SymmetryAxis    = 11 '
      CALL abort(&
          __STAMP__&
          ,'PartBound%TargetBoundCond(BPO%PartBoundaries(iPartBound)) is not implemented for CalcBoundaryParticleOutput',&
          IntInfoOpt=PartBound%TargetBoundCond(BPO%PartBoundaries(iPartBound)))
    END IF
  END DO

END IF ! CalcBoundaryParticleOutput

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
USE MOD_Particle_Boundary_Vars    ,ONLY: nComputeNodeSurfSides,PartBound
#if USE_MPI
USE MOD_Particle_MPI_Vars         ,ONLY: PartMPI
#endif /*USE_MPI*/
USE MOD_SurfaceModel_Vars         ,ONLY: nPorousBC
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
LOGICAL             :: isOpen, isRestart
CHARACTER(LEN=350)  :: outfile
INTEGER             :: unit_index, OutputCounter
INTEGER             :: SurfCollNum(nSpecies),AdsorptionNum(nSpecies),DesorptionNum(nSpecies)
INTEGER             :: iPartBound,iSpec
!===================================================================================================================================
IF ((nComputeNodeSurfSides.EQ.0).AND.(.NOT.CalcBoundaryParticleOutput)) RETURN
isRestart = .FALSE.
IF ( DoRestart ) THEN
  isRestart = .TRUE.
END IF
IF (.NOT.DoSurfModelAnalyze) RETURN
OutputCounter = 2
unit_index = 636
#if USE_MPI
IF (PartMPI%MPIRoot) THEN
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
      IF (CalcSurfCollCounter) THEN
        CALL WriteDataHeaderInfo(unit_index,'nSurfColl-Spec',OutputCounter,nSpecies)
        CALL WriteDataHeaderInfo(unit_index,'N_Ads-Spec',OutputCounter,nSpecies)
        CALL WriteDataHeaderInfo(unit_index,'N_Des-Spec',OutputCounter,nSpecies)
      END IF
      IF(CalcPorousBCInfo) THEN ! calculate porous boundary condition output (pumping speed, removal probability, pressure
                                ! normalized with given pressure. Values are averaged over the whole porous BC
        CALL WriteDataHeaderInfo(unit_index,'PumpSpeed-Measure-PorousBC',OutputCounter,nPorousBC)
        CALL WriteDataHeaderInfo(unit_index,'PumpSpeed-Control-PorousBC',OutputCounter,nPorousBC)
        CALL WriteDataHeaderInfo(unit_index,'RemovalProbability-PorousBC',OutputCounter,nPorousBC)
        CALL WriteDataHeaderInfo(unit_index,'PressureNorm-PorousBC',OutputCounter,nPorousBC)
      END IF
      IF (CalcBoundaryParticleOutput) THEN
        DO iPartBound = 1, BPO%NPartBoundaries
          DO iSpec = 1, BPO%NSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A19,I3.3,A1,A)',ADVANCE='NO') OutputCounter,'-nRealPartOut-Spec-', BPO%Species(iSpec),'-',&
                TRIM(PartBound%SourceBoundName(BPO%PartBoundaries(iPartBound)))
            OutputCounter = OutputCounter + 1
          END DO
        END DO
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
IF (CalcSurfCollCounter)        CALL GetCollCounter(SurfCollNum,AdsorptionNum,DesorptionNum)
IF (CalcPorousBCInfo)           CALL GetPorousBCInfo()
IF (CalcBoundaryParticleOutput) CALL GetBoundaryParticleOutput()
!===================================================================================================================================
! Output Analyzed variables
!===================================================================================================================================
#if USE_MPI
IF (PartMPI%MPIRoot) THEN
#endif /*USE_MPI*/
  WRITE(unit_index,'(E23.16E3)',ADVANCE='NO') Time
  IF (CalcSurfCollCounter) THEN
    CALL WriteDataInfo(unit_index,nSpecies,IntegerArray=SurfCollNum(:))
    CALL WriteDataInfo(unit_index,nSpecies,IntegerArray=AdsorptionNum(:))
    CALL WriteDataInfo(unit_index,nSpecies,IntegerArray=DesorptionNum(:))
  END IF
  IF(CalcPorousBCInfo) THEN
    CALL WriteDataInfo(unit_index,nPorousBC,RealArray=PorousBCOutput(2,:))
    CALL WriteDataInfo(unit_index,nPorousBC,RealArray=PorousBCOutput(3,:))
    CALL WriteDataInfo(unit_index,nPorousBC,RealArray=PorousBCOutput(4,:))
    CALL WriteDataInfo(unit_index,nPorousBC,RealArray=PorousBCOutput(5,:))
  END IF
  IF (CalcBoundaryParticleOutput) THEN
    DO iPartBound = 1, BPO%NPartBoundaries
      DO iSpec = 1, BPO%NSpecies
        CALL WriteDataInfo(unit_index,1,RealArray=(/BPO%RealPartOut(iPartBound,iSpec)/))
        ! Reset PartMPI%MPIRoot counters after writing the data to the file,
        ! non-PartMPI%MPIRoot are reset in GetBoundaryParticleOutput()
        BPO%RealPartOut(iPartBound,iSpec) = 0.
      END DO
    END DO
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


SUBROUTINE GetCollCounter(SurfCollNum,AdsorbNum, DesorbNum)
!===================================================================================================================================
!> Calculates species counters for surface collisions: total, absorbed and desorbed
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars             ,ONLY: nSpecies
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: SurfAnalyzeCount, SurfAnalyzeNumOfAds, SurfAnalyzeNumOfDes
#if USE_MPI
USE MOD_Particle_MPI_Vars         ,ONLY: PartMPI
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(OUT) :: SurfCollNum(nSpecies), AdsorbNum(nSpecies), DesorbNum(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iSpec
!===================================================================================================================================

DO iSpec = 1,nSpecies
  SurfCollNum(iSpec) = SurfAnalyzeCount(iSpec)
  AdsorbNum(iSpec) = SurfAnalyzeNumOfAds(iSpec)
  DesorbNum(iSpec) = SurfAnalyzeNumOfDes(iSpec)
END DO

#if USE_MPI
IF (PartMPI%MPIRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,SurfCollNum ,nSpecies,MPI_INTEGER,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,AdsorbNum   ,nSpecies,MPI_INTEGER,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,DesorbNum   ,nSpecies,MPI_INTEGER,MPI_SUM,0,PartMPI%COMM,IERROR)
ELSE
  CALL MPI_REDUCE(SurfCollNum ,SurfCollNum ,nSpecies,MPI_INTEGER,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(AdsorbNum   ,AdsorbNum   ,nSpecies,MPI_INTEGER,MPI_SUM,0,PartMPI%COMM,IERROR)
  CALL MPI_REDUCE(DesorbNum   ,DesorbNum   ,nSpecies,MPI_INTEGER,MPI_SUM,0,PartMPI%COMM,IERROR)
END IF
#endif /*USE_MPI*/

! Reset counters
SurfAnalyzeCount = 0
SurfAnalyzeNumOfAds = 0
SurfAnalyzeNumOfDes = 0

END SUBROUTINE GetCollCounter


SUBROUTINE GetPorousBCInfo()
!===================================================================================================================================
!> 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_SurfaceModel_Vars         ,ONLY: nPorousBC
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: PorousBCOutput
#if USE_MPI
USE MOD_Particle_MPI_Vars         ,ONLY: PartMPI
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iPBC
!===================================================================================================================================
#if USE_MPI
IF (PartMPI%MPIRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,PorousBCOutput,5*nPorousBC,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,iError)
ELSE
  CALL MPI_REDUCE(PorousBCOutput,PorousBCOutput,5*nPorousBC,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,iError)
END IF
#endif /*USE_MPI*/

#if USE_MPI
IF (PartMPI%MPIRoot) THEN
#endif /*USE_MPI*/
  DO iPBC = 1, nPorousBC
    IF(PorousBCOutput(1,iPBC).GT.0.0) THEN
      ! Pumping Speed (Output(2)) is the sum of all elements (counter over particles exiting through pump)
      ! Other variales are averaged over the elements
      PorousBCOutput(3:5,iPBC) = PorousBCOutput(3:5,iPBC) / PorousBCOutput(1,iPBC)
    END IF
  END DO
#if USE_MPI
END IF
#endif /*USE_MPI*/

END SUBROUTINE GetPorousBCInfo


SUBROUTINE GetBoundaryParticleOutput()
!===================================================================================================================================
!> Synchronize BoundaryParticleOutput analyze arrays
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: BPO
#if USE_MPI
USE MOD_Particle_MPI_Vars         ,ONLY: PartMPI
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: SendBuf(1:BPO%NPartBoundaries*BPO%NSpecies)
INTEGER :: SendBufSize
!===================================================================================================================================
#if USE_MPI
SendBufSize = BPO%NPartBoundaries*BPO%NSpecies
IF (PartMPI%MPIRoot) THEN
  ! Map 2D array to vector for sending via MPI
  SendBuf = RESHAPE(BPO%RealPartOut,(/SendBufSize/))
  CALL MPI_REDUCE(MPI_IN_PLACE,SendBuf(1:SendBufSize),SendBufSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  ! MAP vector back to 2D array
  BPO%RealPartOut = RESHAPE(SendBuf,(/BPO%NPartBoundaries,BPO%NSpecies/))
ELSE
  ! Map 2D array to vector for sending via MPI
  SendBuf = RESHAPE(BPO%RealPartOut,(/SendBufSize/))
  CALL MPI_REDUCE(SendBuf(1:SendBufSize),0,SendBufSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  ! Reset non PartMPI%MPIRoot counters, PartMPI%MPIRoot counters are reset after writing the data to the file
  BPO%RealPartOut = 0.
END IF
#endif /*USE_MPI*/

END SUBROUTINE GetBoundaryParticleOutput
#endif /*PARTICLES*/

END MODULE MOD_SurfaceModel_Analyze
