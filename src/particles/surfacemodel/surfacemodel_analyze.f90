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
PUBLIC:: FinalizeSurfaceModelAnalyze
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
CALL prms%CreateLogicalOption(  'CalcElectronSEE'         , 'Count the electron emission from BCs where SEE is active','.FALSE.')

END SUBROUTINE DefineParametersSurfModelAnalyze


!===================================================================================================================================
!> Initializes variables necessary for analyse subroutines
!===================================================================================================================================
SUBROUTINE InitSurfModelAnalyze()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_ReadInTools               ,ONLY: GETLOGICAL,GETINT,GETINTARRAY
USE MOD_Particle_Vars             ,ONLY: nSpecies
USE MOD_Analyze_Vars              ,ONLY: DoSurfModelAnalyze
USE MOD_SurfaceModel_Vars         ,ONLY: nPorousBC
USE MOD_SurfaceModel_Analyze_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(SurfModelAnalyzeInitIsDone)THEN
  CALL abort(__STAMP__,'InitParticleAnalyse already called.')
  RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT SURFACE MODEL ANALYZE...'

SurfaceAnalyzeStep = GETINT('Surface-AnalyzeStep')
IF(SurfaceAnalyzeStep.EQ.0) SurfaceAnalyzeStep = HUGE(1)

DoSurfModelAnalyze = .FALSE.

!-- Surface Collision Counter
CalcSurfCollCounter = GETLOGICAL('Surf-CalcCollCounter')
IF(CalcSurfCollCounter)THEN
  DoSurfModelAnalyze = .TRUE.
  ! allocate info and constants
  ALLOCATE(SurfAnalyzeCount(1:nSpecies),SurfAnalyzeNumOfAds(1:nSpecies),SurfAnalyzeNumOfDes(1:nSpecies))
  SurfAnalyzeCount = 0; SurfAnalyzeNumOfAds = 0;  SurfAnalyzeNumOfDes = 0
END IF

!-- Porous Boundaries
IF(nPorousBC.GT.0)THEN
  ! Output for porous BC: Pump averaged values
  CalcPorousBCInfo = GETLOGICAL('Surf-CalcPorousBCInfo')
  IF(CalcPorousBCInfo)THEN
    DoSurfModelAnalyze = .TRUE.
    ALLOCATE(PorousBCOutput(1:5,1:nPorousBC))
    PorousBCOutput = 0.
  END IF
END IF

!-- BoundaryParticleOutput (after mapping of PartBound on FieldBound and determination of PartBound types = open, reflective etc.)
CalcBoundaryParticleOutput = GETLOGICAL('CalcBoundaryParticleOutput')
IF(CalcBoundaryParticleOutput) CALL InitBoundaryParticleOutput()

!-- Electron SEE emission counter
CalcElectronSEE = GETLOGICAL('CalcElectronSEE','.FALSE.')
IF(CalcElectronSEE) CALL InitCalcElectronSEE()

SurfModelAnalyzeInitIsDone=.TRUE.

SWRITE(UNIT_stdOut,'(A)')' INIT SURFACE MODEL ANALYZE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitSurfModelAnalyze


!===================================================================================================================================
!> create/open SurfaceAnalyze.csv and write calculated variables for surface analyze
!===================================================================================================================================
SUBROUTINE AnalyzeSurface(Time)
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
USE MOD_Particle_Vars             ,ONLY: nSpecies,UseNeutralization,NeutralizationBalanceGlobal
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: Time
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL             :: isOpen
CHARACTER(LEN=350)  :: outfile
INTEGER             :: unit_index, OutputCounter
INTEGER             :: SurfCollNum(nSpecies),AdsorptionNum(nSpecies),DesorptionNum(nSpecies)
INTEGER             :: iPartBound,iSpec
!===================================================================================================================================
IF((nComputeNodeSurfSides.EQ.0).AND.(.NOT.CalcBoundaryParticleOutput).AND.(.NOT.UseNeutralization).AND.(.NOT.CalcElectronSEE)) RETURN
IF(.NOT.DoSurfModelAnalyze) RETURN
OutputCounter = 2
unit_index = 636
#if USE_MPI
IF(PartMPI%MPIRoot)THEN
#endif /*USE_MPI*/
  INQUIRE(UNIT   = unit_index , OPENED = isOpen)
  IF(.NOT.isOpen)THEN
    outfile = 'SurfaceAnalyze.csv'
!===================================================================================================================================
! Write Header
!===================================================================================================================================
    IF(DoRestart.AND.FILEEXISTS(outfile))THEN
      OPEN(unit_index,file=TRIM(outfile),position="APPEND",status="OLD")
    ELSE
      OPEN(unit_index,file=TRIM(outfile))
      !--- insert header
      WRITE(unit_index,'(A8)',ADVANCE='NO') '001-TIME'
      IF(CalcSurfCollCounter)THEN
        CALL WriteDataHeaderInfo(unit_index,'nSurfColl-Spec',OutputCounter,nSpecies)
        CALL WriteDataHeaderInfo(unit_index,'N_Ads-Spec',OutputCounter,nSpecies)
        CALL WriteDataHeaderInfo(unit_index,'N_Des-Spec',OutputCounter,nSpecies)
      END IF
      IF(CalcPorousBCInfo)THEN ! calculate porous boundary condition output (pumping speed, removal probability, pressure
                                ! normalized with given pressure. Values are averaged over the whole porous BC
        CALL WriteDataHeaderInfo(unit_index,'PumpSpeed-Measure-PorousBC',OutputCounter,nPorousBC)
        CALL WriteDataHeaderInfo(unit_index,'PumpSpeed-Control-PorousBC',OutputCounter,nPorousBC)
        CALL WriteDataHeaderInfo(unit_index,'RemovalProbability-PorousBC',OutputCounter,nPorousBC)
        CALL WriteDataHeaderInfo(unit_index,'PressureNorm-PorousBC',OutputCounter,nPorousBC)
      END IF
      IF(CalcBoundaryParticleOutput)THEN
        DO iPartBound = 1, BPO%NPartBoundaries
          DO iSpec = 1, BPO%NSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A19,I3.3,A1,A)',ADVANCE='NO') OutputCounter,'-nRealPartOut-Spec-', BPO%Species(iSpec),'-',&
                TRIM(PartBound%SourceBoundName(BPO%PartBoundaries(iPartBound)))
            OutputCounter = OutputCounter + 1
          END DO
        END DO
      END IF
      IF(UseNeutralization)THEN ! Ion thruster neutralization current (virtual cathode electrons)
        WRITE(unit_index,'(A1,I3.3,A)',ADVANCE='NO') ',',OutputCounter,'-NeutralizationParticles'
        OutputCounter = OutputCounter + 1
      END IF ! UseNeutralization
      IF(CalcElectronSEE)THEN
        DO iPartBound = 1, SEE%NPartBoundaries
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A)',ADVANCE='NO') OutputCounter,'-nRealElectronsEmmited-'//&
              TRIM(PartBound%SourceBoundName(SEE%PartBoundaries(iPartBound)))
          OutputCounter = OutputCounter + 1
        END DO ! iPartBound = 1, SEE%NPartBoundaries
      END IF ! CalcElectronSEE
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
IF (CalcBoundaryParticleOutput) CALL SyncBoundaryParticleOutput()
IF (CalcElectronSEE)            CALL SyncElectronSEE()
!===================================================================================================================================
! Output Analyzed variables
!===================================================================================================================================
#if USE_MPI
IF(PartMPI%MPIRoot)THEN
#endif /*USE_MPI*/
  WRITE(unit_index,'(E23.16E3)',ADVANCE='NO') Time
  IF(CalcSurfCollCounter)THEN
    CALL WriteDataInfo(unit_index,nSpecies,IntegerArray=SurfCollNum(:))
    CALL WriteDataInfo(unit_index,nSpecies,IntegerArray=AdsorptionNum(:))
    CALL WriteDataInfo(unit_index,nSpecies,IntegerArray=DesorptionNum(:))
  END IF
  IF(CalcPorousBCInfo)THEN
    CALL WriteDataInfo(unit_index,nPorousBC,RealArray=PorousBCOutput(2,:))
    CALL WriteDataInfo(unit_index,nPorousBC,RealArray=PorousBCOutput(3,:))
    CALL WriteDataInfo(unit_index,nPorousBC,RealArray=PorousBCOutput(4,:))
    CALL WriteDataInfo(unit_index,nPorousBC,RealArray=PorousBCOutput(5,:))
  END IF
  IF(CalcBoundaryParticleOutput)THEN
    DO iPartBound = 1, BPO%NPartBoundaries
      DO iSpec = 1, BPO%NSpecies
        CALL WriteDataInfo(unit_index,1,RealArray=(/BPO%RealPartOut(iPartBound,iSpec)/))
        ! Reset PartMPI%MPIRoot counters after writing the data to the file,
        ! non-PartMPI%MPIRoot are reset in SyncBoundaryParticleOutput()
        BPO%RealPartOut(iPartBound,iSpec) = 0.
      END DO
    END DO
  END IF
  IF(UseNeutralization) CALL WriteDataInfo(unit_index,1,RealArray=(/REAL(NeutralizationBalanceGlobal)/))
  IF(CalcElectronSEE)THEN
    DO iPartBound = 1, SEE%NPartBoundaries
      CALL WriteDataInfo(unit_index,1,RealArray=(/SEE%RealElectronOut(iPartBound)/))
        ! Reset PartMPI%MPIRoot counters after writing the data to the file,
        ! non-PartMPI%MPIRoot are reset in SyncBoundaryParticleOutput()
        SEE%RealElectronOut(iPartBound) = 0.
    END DO ! iPartBound = 1, SEE%NPartBoundaries
  END IF ! CalcElectronSEE
  WRITE(unit_index,'(A)') ''
#if USE_MPI
END IF
#endif /*USE_MPI*/
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE AnalyzeSurface


!===================================================================================================================================
!> writes OutputCounter-AttribNamestring-iLoop into WRITEFORMAT output
!===================================================================================================================================
SUBROUTINE WriteDataHeaderInfo(unit_index,AttribName,OutputCounter,LoopSize)
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


!===================================================================================================================================
!> writes INPUTData into unit_index output
!> only one data input should be given at a time
!===================================================================================================================================
SUBROUTINE WriteDataInfo(unit_index,nVal,RealScalar,IntegerScalar,StrScalar,LogicalScalar, &
                                  RealArray,IntegerArray,IntegerK8Array,StrArray)
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
IF(PRESENT(RealArray))THEN
  DO iLoop = 1, nVal
    WRITE (unit_index, CSVFORMAT, ADVANCE='NO') ',',RealArray(iLoop)
  END DO
END IF
IF(PRESENT(RealScalar))THEN
  WRITE (unit_index, CSVFORMAT, ADVANCE='NO') ',',RealScalar
END IF

IF(PRESENT(IntegerArray))THEN
  DO iLoop = 1, nVal
    WRITE (unit_index, CSVFORMAT, ADVANCE='NO') ',',REAL(IntegerArray(iLoop))
  END DO
END IF

IF(PRESENT(IntegerK8Array))THEN
  DO iLoop = 1, nVal
    WRITE (unit_index, CSVFORMAT, ADVANCE='NO') ',',REAL(IntegerK8Array(iLoop))
  END DO
END IF

IF(PRESENT(IntegerScalar))THEN
  WRITE (unit_index, CSVFORMAT, ADVANCE='NO') ',',REAL(IntegerScalar)
END IF

IF(PRESENT(StrArray))THEN
  DO iLoop = 1, nVal
    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
    WRITE(unit_index,'(A)',ADVANCE='NO') TRIM(StrArray(iLoop))
  END DO
END IF

IF(PRESENT(StrScalar))THEN
  WRITE(unit_index,'(A1)',ADVANCE='NO') ','
  WRITE(unit_index,'(A)',ADVANCE='NO') TRIM(StrScalar)
END IF

IF(PRESENT(LogicalScalar))THEN
  WRITE(unit_index,'(A1)',ADVANCE='NO') ','
  WRITE(unit_index,'(I1)',ADVANCE='NO') LogicalScalar
END IF
END SUBROUTINE WriteDataInfo


!===================================================================================================================================
!> Calculates species counters for surface collisions: total, absorbed and desorbed
!===================================================================================================================================
SUBROUTINE GetCollCounter(SurfCollNum,AdsorbNum, DesorbNum)
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
IF(PartMPI%MPIRoot)THEN
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


!===================================================================================================================================
!> 
!===================================================================================================================================
SUBROUTINE GetPorousBCInfo()
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
IF(PartMPI%MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,PorousBCOutput,5*nPorousBC,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,iError)
ELSE
  CALL MPI_REDUCE(PorousBCOutput,PorousBCOutput,5*nPorousBC,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,iError)
END IF
#endif /*USE_MPI*/

#if USE_MPI
IF(PartMPI%MPIRoot)THEN
#endif /*USE_MPI*/
  DO iPBC = 1, nPorousBC
    IF(PorousBCOutput(1,iPBC).GT.0.0)THEN
      ! Pumping Speed (Output(2)) is the sum of all elements (counter over particles exiting through pump)
      ! Other variables are averaged over the elements
      PorousBCOutput(3:5,iPBC) = PorousBCOutput(3:5,iPBC) / PorousBCOutput(1,iPBC)
    END IF
  END DO
#if USE_MPI
END IF
#endif /*USE_MPI*/

END SUBROUTINE GetPorousBCInfo


!===================================================================================================================================
!> Synchronize BoundaryParticleOutput analyze arrays
!===================================================================================================================================
SUBROUTINE SyncBoundaryParticleOutput()
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
IF(PartMPI%MPIRoot)THEN
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

END SUBROUTINE SyncBoundaryParticleOutput


!===================================================================================================================================
!> Synchronize CalcElectronSEE analyze arrays
!===================================================================================================================================
SUBROUTINE SyncElectronSEE()
! MODULES
USE MOD_Globals
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: SEE
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
!===================================================================================================================================
#if USE_MPI
IF (PartMPI%MPIRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,SEE%RealElectronOut,SEE%NPartBoundaries,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
ELSE
  CALL MPI_REDUCE(SEE%RealElectronOut,0,SEE%NPartBoundaries,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  ! Reset non PartMPI%MPIRoot counters, PartMPI%MPIRoot counters are reset after writing the data to the file
  SEE%RealElectronOut = 0.
END IF
#endif /*USE_MPI*/

END SUBROUTINE SyncElectronSEE


!===================================================================================================================================
!> Allocate the required arrays (mappings and containers) for BoundaryParticleOutput
!===================================================================================================================================
SUBROUTINE InitBoundaryParticleOutput()
! MODULES
USE MOD_Globals                   ,ONLY: abort,UNIT_stdOut,MPIRoot
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: BPO
USE MOD_Particle_Boundary_Vars    ,ONLY: nPartBound,PartBound
USE MOD_ReadInTools               ,ONLY: GETLOGICAL,GETINT,GETINTARRAY
USE MOD_Analyze_Vars              ,ONLY: DoSurfModelAnalyze
USE MOD_Particle_Vars             ,ONLY: nSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iPartBound,iSpec
!===================================================================================================================================
DoSurfModelAnalyze = .TRUE.
BPO%NPartBoundaries = GETINT('BPO-NPartBoundaries')
BPO%PartBoundaries  = GETINTARRAY('BPO-PartBoundaries',BPO%NPartBoundaries)
BPO%NSpecies        = GETINT('BPO-NSpecies')
BPO%Species         = GETINTARRAY('BPO-Species',BPO%NSpecies)
IF(BPO%NPartBoundaries.EQ.0.OR.BPO%NSpecies.EQ.0) CALL abort(__STAMP__,'BPO-NPartBoundaries or BPO-NSpecies is zero.')
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
    CALL abort(__STAMP__&
        ,'PartBound%TargetBoundCond(BPO%PartBoundaries(iPartBound)) is not implemented for CalcBoundaryParticleOutput',&
        IntInfoOpt=PartBound%TargetBoundCond(BPO%PartBoundaries(iPartBound)))
  END IF ! .NOT.ANY(PartBound%TargetBoundCond(BPO%PartBoundaries(iPartBound)).EQ. ...
END DO ! iPartBound = 1, BPO%NPartBoundaries

END SUBROUTINE InitBoundaryParticleOutput


!===================================================================================================================================
!> Allocate the required arrays (mappings and containers) for secondary electron emission analysis
!===================================================================================================================================
SUBROUTINE InitCalcElectronSEE()
! MODULES
USE MOD_Globals                   ,ONLY: abort!,UNIT_stdOut,MPIRoot
USE MOD_Analyze_Vars              ,ONLY: DoSurfModelAnalyze
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: SEE
USE MOD_Particle_Boundary_Vars    ,ONLY: nPartBound,PartBound
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iPartBound
!===================================================================================================================================
DoSurfModelAnalyze = .TRUE.

! Count number of different SEE boundaries
SEE%NPartBoundaries = 0
DO iPartBound=1,nPartBound
  IF(.NOT.PartBound%Reactive(iPartBound)) CYCLE
  SELECT CASE(PartBound%SurfaceModel(iPartBound))
  CASE(SEE_MODELS_ID)
    SEE%NPartBoundaries = SEE%NPartBoundaries +1
  END SELECT
END DO ! iPartBound=1,nPartBound

! Sanity check
IF(SEE%NPartBoundaries.EQ.0) CALL abort(__STAMP__,'No SEE boundaries found for counting the emitted electrons')

! Create Mapping
ALLOCATE(SEE%PartBoundaries(SEE%NPartBoundaries))
SEE%NPartBoundaries = 0
DO iPartBound=1,nPartBound
  IF(.NOT.PartBound%Reactive(iPartBound)) CYCLE
  SELECT CASE(PartBound%SurfaceModel(iPartBound))
  CASE(SEE_MODELS_ID)
    SEE%NPartBoundaries = SEE%NPartBoundaries +1
    SEE%PartBoundaries(SEE%NPartBoundaries) = iPartBound
  END SELECT
END DO ! iPartBound=1,nPartBound

! Allocate the container
ALLOCATE(SEE%RealElectronOut(1:SEE%NPartBoundaries))
SEE%RealElectronOut = 0.

! Create Mapping
ALLOCATE(SEE%BCIDToSEEBCID(1:nPartBound))
SEE%BCIDToSEEBCID     = -1
DO iPartBound = 1, SEE%NPartBoundaries
  SEE%BCIDToSEEBCID(SEE%PartBoundaries(iPartBound)) = iPartBound
END DO ! iPartBound = 1, BPO%NPartBoundaries

END SUBROUTINE InitCalcElectronSEE


!===================================================================================================================================
!> Deallocate surface model vars
!===================================================================================================================================
SUBROUTINE FinalizeSurfaceModelAnalyze()
! MODULES
USE MOD_SurfaceModel_Analyze_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! === Surface Analyze Vars
SurfModelAnalyzeInitIsDone=.FALSE.

SDEALLOCATE(SurfAnalyzeCount)
SDEALLOCATE(SurfAnalyzeNumOfAds)
SDEALLOCATE(SurfAnalyzeNumOfDes)

! Boundary Particle Output
IF(CalcBoundaryParticleOutput)THEN
  SDEALLOCATE(BPO%RealPartOut)
  SDEALLOCATE(BPO%PartBoundaries)
  SDEALLOCATE(BPO%BCIDToBPOBCID)
  SDEALLOCATE(BPO%Species)
  SDEALLOCATE(BPO%SpecIDToBPOSpecID)
END IF ! CalcBoundaryParticleOutput

! Electron emission counter
IF(CalcElectronSEE)THEN
  SDEALLOCATE(SEE%RealElectronOut)
  SDEALLOCATE(SEE%PartBoundaries)
  SDEALLOCATE(SEE%BCIDToSEEBCID)
END IF ! CalcElectronSEE

END SUBROUTINE FinalizeSurfaceModelAnalyze

#endif /*PARTICLES*/


END MODULE MOD_SurfaceModel_Analyze
