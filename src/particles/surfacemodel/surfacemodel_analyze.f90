!==================================================================================================================================
! Copyright (c) 2015 - 2019 Wladimir Reschke
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
                                                          'probability and pressure. Values are averaged over the whole porous BC.'//&
                                                          'Disabled per default, but automatically enabled if a sensor is detected.')
!-- BoundaryParticleOutput
CALL prms%CreateLogicalOption(  'CalcBoundaryParticleOutput', 'Count number of particles exiting for species X on boundary X' , '.FALSE.')
CALL prms%CreateIntOption(      'BPO-NPartBoundaries'       , 'Number of boundaries used for CalcBoundaryParticleOutput')
CALL prms%CreateIntArrayOption( 'BPO-PartBoundaries'        , 'Vector (length BPO-NPartBoundaries) with the numbers of each Part-Boundary', no=0)
CALL prms%CreateIntOption(      'BPO-NSpecies'              , 'Number of species used for CalcBoundaryParticleOutput')
CALL prms%CreateIntArrayOption( 'BPO-Species'               , 'Vector (length BPO-NSpecies) with the corresponding Species IDs', no=0)
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
USE MOD_SurfaceModel_Vars         ,ONLY: nPorousBC, PorousBC
USE MOD_SurfaceModel_Analyze_Vars
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
USE MOD_Restart_Vars              ,ONLY: DoRestart,RestartTime
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
LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT SURFACE MODEL ANALYZE...'

! Initialize with restart time or zero
IF(DoRestart)THEN
  SurfModelAnalyzeSampleTime = RestartTime
ELSE
  SurfModelAnalyzeSampleTime = 0.
END IF ! DoRestart

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
  IF(ANY(PorousBC(:)%Type.EQ.'sensor')) THEN
    ! If a sensor was defined, set the default value to TRUE
    CalcPorousBCInfo = GETLOGICAL('Surf-CalcPorousBCInfo','.TRUE.')
  ELSE
    CalcPorousBCInfo = GETLOGICAL('Surf-CalcPorousBCInfo','.FALSE.')
  END IF
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
CALL InitCalcElectronSEE() ! This routine calls GETLOGICAL('CalcElectronSEE','.FALSE.')

SurfModelAnalyzeInitIsDone=.TRUE.

LBWRITE(UNIT_stdOut,'(A)')' INIT SURFACE MODEL ANALYZE DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')

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
USE MOD_SurfaceModel_Vars         ,ONLY: nPorousBC, PorousBC
USE MOD_Particle_Vars             ,ONLY: nSpecies,UseNeutralization,NeutralizationBalanceGlobal,Species
#if USE_HDG
USE MOD_Analyze_Vars              ,ONLY: EDC
USE MOD_Analyze_Vars              ,ONLY: CalcElectricTimeDerivative
USE MOD_HDG_Vars                  ,ONLY: UseBiasVoltage,BiasVoltage,BVDataLength
#if USE_MPI
USE MOD_HDG                       ,ONLY: SynchronizeBV
#endif /*USE_MPI*/
#endif /*USE_HDG*/
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
INTEGER             :: iBC,iPartBound,iSEE,iBPO,iSpec
REAL                :: charge,TotalElectricCharge
#if USE_HDG
INTEGER             :: iEDCBC,i,iBoundary,iPartBound2
#endif /*USE_HDG*/
!===================================================================================================================================
IF((nComputeNodeSurfSides.EQ.0).AND.(.NOT.CalcBoundaryParticleOutput).AND.(.NOT.UseNeutralization).AND.(.NOT.CalcElectronSEE)) RETURN
IF(.NOT.DoSurfModelAnalyze) RETURN
SurfModelAnalyzeSampleTime = Time - SurfModelAnalyzeSampleTime ! Set SurfModelAnalyzeSampleTime=Time at the end of this routine
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
        CALL WriteDataHeaderInfo(unit_index,'nSurfColl-Spec',OutputCounter,LoopSize=nSpecies)
        CALL WriteDataHeaderInfo(unit_index,'N_Ads-Spec',OutputCounter,LoopSize=nSpecies)
        CALL WriteDataHeaderInfo(unit_index,'N_Des-Spec',OutputCounter,LoopSize=nSpecies)
      END IF
      ! Calculate porous boundary condition output (pumping speed, removal probability, pressure)
      ! Values are averaged over the whole porous BC
      IF(CalcPorousBCInfo)THEN
        DO iBC = 1, nPorousBC
          IF(PorousBC(iBC)%Type.EQ.'pump') THEN
            CALL WriteDataHeaderInfo(unit_index,'PumpSpeed-Measure-Pump',OutputCounter,iLoop_in=iBC)
            CALL WriteDataHeaderInfo(unit_index,'PumpSpeed-Control-Pump',OutputCounter,iLoop_in=iBC)
            CALL WriteDataHeaderInfo(unit_index,'RemovalProbability-Pump',OutputCounter,iLoop_in=iBC)
            CALL WriteDataHeaderInfo(unit_index,'Pressure-Pump',OutputCounter,iLoop_in=iBC)
          ELSE IF(PorousBC(iBC)%Type.EQ.'sensor') THEN
            CALL WriteDataHeaderInfo(unit_index,'Pressure-Sensor',OutputCounter,iLoop_in=iBC)
          END IF
        END DO
      END IF
      IF(CalcBoundaryParticleOutput)THEN
        ! Output of particle fluxes
        DO iBPO = 1, BPO%NPartBoundaries
          DO iSpec = 1, BPO%NSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I3.3,A1,A)',ADVANCE='NO') OutputCounter,'-Flux-Spec-', BPO%Species(iSpec),'-',&
                TRIM(PartBound%SourceBoundName(BPO%PartBoundaries(iBPO)))
            OutputCounter = OutputCounter + 1
          END DO
        END DO
        ! Output of total electric current
        IF(BPO%OutputTotalElectricCurrent)THEN
          DO iBPO = 1, BPO%NPartBoundaries
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A)',ADVANCE='NO') OutputCounter,'-TotalElectricCurrent-'//&
                TRIM(PartBound%SourceBoundName(BPO%PartBoundaries(iBPO)))
            OutputCounter = OutputCounter + 1
          END DO
        END IF ! BPO%OutputTotalElectricCurrent
      END IF
#if USE_HDG
      ! Output of bias voltage variables
      IF(UseBiasVoltage)THEN
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,'(I3.3,A)',ADVANCE='NO') OutputCounter,'-BiasVoltage[V]'
        OutputCounter = OutputCounter + 1
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,'(I3.3,A)',ADVANCE='NO') OutputCounter,'-BiasVoltageIonExcess[C]'
        OutputCounter = OutputCounter + 1
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,'(I3.3,A)',ADVANCE='NO') OutputCounter,'-BiasVoltageUpdateTime[s]'
        OutputCounter = OutputCounter + 1
      END IF ! UseBiasVoltage
#endif /*USE_HDG*/
      IF(UseNeutralization)THEN ! Ion thruster neutralization current (virtual cathode electrons)
        CALL WriteDataHeaderInfo(unit_index,'NeutralizationParticles',OutputCounter)
      END IF ! UseNeutralization
      IF(CalcElectronSEE)THEN
        DO iSEE = 1, SEE%NPartBoundaries
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A)',ADVANCE='NO') OutputCounter,'-ElectricCurrentSEE-'//&
              TRIM(PartBound%SourceBoundName(SEE%PartBoundaries(iSEE)))
          OutputCounter = OutputCounter + 1
        END DO ! iSEE = 1, SEE%NPartBoundaries
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
    DO iBC = 1, nPorousBC
      IF(PorousBC(iBC)%Type.EQ.'pump') THEN
        CALL WriteDataInfo(unit_index,RealScalar=PorousBCOutput(2,iBC))
        CALL WriteDataInfo(unit_index,RealScalar=PorousBCOutput(3,iBC))
        CALL WriteDataInfo(unit_index,RealScalar=PorousBCOutput(4,iBC))
      END IF
      CALL WriteDataInfo(unit_index,RealScalar=PorousBCOutput(5,iBC))
    END DO
  END IF
  IF(CalcBoundaryParticleOutput)THEN
    ! Output of particle fluxes
    DO iBPO = 1, BPO%NPartBoundaries
      DO iSpec = 1, BPO%NSpecies
        IF(ABS(SurfModelAnalyzeSampleTime).LE.0.0)THEN
          CALL WriteDataInfo(unit_index,RealScalar=0.0)
        ELSE
          CALL WriteDataInfo(unit_index,RealScalar=BPO%RealPartOut(iBPO,iSpec)/SurfModelAnalyzeSampleTime)
        END IF ! ABS(SurfModelAnalyzeSampleTime).LE.0.0
      END DO
    END DO
    ! Output of total electric current
    IF(BPO%OutputTotalElectricCurrent)THEN
      DO iBPO = 1, BPO%NPartBoundaries
        TotalElectricCharge = 0. ! Initialize total charge for each boundary
        IF(ABS(SurfModelAnalyzeSampleTime).LE.0.0)THEN
          CALL WriteDataInfo(unit_index,1,RealArray=(/0.0/))
        ELSE
          ! Sum over all fluxes (only if the species has a charge)
          DO iSpec = 1, BPO%NSpecies
            charge = Species(BPO%Species(iSpec))%ChargeIC
            ! Impacting charged particles: positive number for positive ions (+) and negative number for electrons (-)
            IF(ABS(charge).GT.0.0) TotalElectricCharge = TotalElectricCharge + BPO%RealPartOut(iBPO,iSpec)*charge
          END DO
          ! Released secondary electrons (always a positive number). SEE%BCIDToSEEBCID(iPartBound) yields the iSEEBCIndex
          IF(CalcElectronSEE)THEN
            ! Get particle boundary ID
            iPartBound = BPO%PartBoundaries(iBPO)
            ! Get SEE ID
            iSEE = SEE%BCIDToSEEBCID(iPartBound)
            ! Add SEE current if this BC has secondary electron emission
            IF(iSEE.GT.0) TotalElectricCharge = TotalElectricCharge + SEE%RealElectronOut(iSEE)
          END IF ! CalcElectronSEE
          TotalElectricCharge = TotalElectricCharge/SurfModelAnalyzeSampleTime
#if USE_HDG
          ! Add electric displacement current to total electric current
          IF(CalcElectricTimeDerivative)THEN
            ! Get iBC index (field) and EDC index (displacement current)
            iBC    = BPO%FieldBoundaries(iBPO)
            iEDCBC = EDC%BCIDToEDCBCID(iBC)
            TotalElectricCharge = TotalElectricCharge + EDC%Current(iEDCBC)
          END IF ! CalcElectricTimeDerivative
#endif /*USE_HDG*/
          ! Sampling time has already been considered due to the displacement current
          CALL WriteDataInfo(unit_index,RealScalar=TotalElectricCharge)
        END IF ! ABS(SurfModelAnalyzeSampleTime).LE.0.0
      END DO ! iBPO = 1, BPO%NPartBoundaries
    END IF ! BPO%OutputTotalElectricCurrent
#if USE_HDG
    ! Bias voltage
    IF(UseBiasVoltage)THEN
      TotalElectricCharge = 0. ! Initialize sum over all boundaries
      ! Ion excess
      DO iBoundary = 1, BiasVoltage%NPartBoundaries
        iPartBound = BiasVoltage%PartBoundaries(iBoundary)
        iBPO = BPO%BCIDToBPOBCID(iPartBound)
        ! Sum over all fluxes (only if the species has a charge)
        DO iSpec = 1, BPO%NSpecies
          charge = Species(BPO%Species(iSpec))%ChargeIC
          ! Impacting charged particles: positive number for positive ions (+) and negative number for electrons (-)
          IF(ABS(charge).GT.0.0) TotalElectricCharge = TotalElectricCharge + BPO%RealPartOut(iBPO,iSpec)*charge
        END DO
        ! Released secondary electrons (always a positive number). SEE%BCIDToSEEBCID(iPartBound) yields the iSEEBCIndex
        IF(CalcElectronSEE)THEN
          ! Get particle boundary ID
          iPartBound2 = BPO%PartBoundaries(iBPO)
          ! Sanity Check
          IF(iPartBound.NE.iPartBound2) CALL abort(__STAMP__,'AnalyzeSurface(): Wrong particle boundary encountered!')
          ! Get SEE ID
          iSEE = SEE%BCIDToSEEBCID(iPartBound)
          ! Add SEE current if this BC has secondary electron emission
          IF(iSEE.GT.0) TotalElectricCharge = TotalElectricCharge + SEE%RealElectronOut(iSEE)
        END IF ! CalcElectronSEE
      END DO ! iBoundary = 1, BiasVoltage%NPartBoundaries

      ASSOCIATE( V => BiasVoltage%BVData(1), Q => BiasVoltage%BVData(2), tBV => BiasVoltage%BVData(3) )
        ! Add total electric charge (without displacement current!)
        Q = Q + TotalElectricCharge

        ! Simulation time threshold
        IF(time.GE.tBV)THEN
          ! Update time
          IF(BiasVoltage%Frequency.GT.0.0) tBV = tBV + 1.0/BiasVoltage%Frequency

          ! Update Voltage
          IF(Q.GT.0.0)THEN
            ! Increase voltage
            V = V + BiasVoltage%Delta
          ELSEIF(Q.LT.0.0)THEN
            ! Decrease voltage
            V = V - BiasVoltage%Delta
          ELSE
            ! do nothing
          END IF ! Q.LT.0.0

          ! Reset ion excess counter
          Q = 0.
        END IF ! time.GE.tBV
      END ASSOCIATE

      ! Write: Voltage, Ion excess and simulation update time
      DO i = 1, BVDataLength
        CALL WriteDataInfo(unit_index, RealScalar=BiasVoltage%BVData(i))
      END DO ! i = 1, 3
    END IF ! UseBiasVoltage
#endif /*USE_HDG*/

    ! Reset BPO containers
    DO iPartBound = 1, BPO%NPartBoundaries
      DO iSpec = 1, BPO%NSpecies
        ! Reset PartMPI%MPIRoot counters after writing the data to the file,
        ! non-PartMPI%MPIRoot are reset in SyncBoundaryParticleOutput()
        BPO%RealPartOut(iPartBound,iSpec) = 0.
      END DO ! iSpec = 1, BPO%NSpecies
    END DO ! iPartBound = 1, BPO%NPartBoundaries
  END IF ! CalcBoundaryParticleOutput

  IF(UseNeutralization) CALL WriteDataInfo(unit_index,RealScalar=REAL(NeutralizationBalanceGlobal))

  IF(CalcElectronSEE)THEN
    DO iPartBound = 1, SEE%NPartBoundaries
      IF(ABS(SurfModelAnalyzeSampleTime).LE.0.0)THEN
        CALL WriteDataInfo(unit_index,RealScalar=0.0)
      ELSE
        CALL WriteDataInfo(unit_index,RealScalar=SEE%RealElectronOut(iPartBound)/SurfModelAnalyzeSampleTime)
      END IF ! ABS(SurfModelAnalyzeSampleTime).LE.0.0
        ! Reset PartMPI%MPIRoot counters after writing the data to the file,
        ! non-PartMPI%MPIRoot are reset in SyncBoundaryParticleOutput()
        SEE%RealElectronOut(iPartBound) = 0.
    END DO ! iPartBound = 1, SEE%NPartBoundaries
  END IF ! CalcElectronSEE
  WRITE(unit_index,'(A)') ''
#if USE_MPI
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
#if USE_HDG
  ! Bias voltage: Update values of sub-communicator processes (broadcast from MPIRoot to all sub-communicator processes)
  IF(UseBiasVoltage) CALL SynchronizeBV()
#endif /*USE_HDG*/
#endif /*USE_MPI*/
!-----------------------------------------------------------------------------------------------------------------------------------
SurfModelAnalyzeSampleTime = Time ! Backup "old" time value for next output
END SUBROUTINE AnalyzeSurface


!===================================================================================================================================
!> Writes OutputCounter-AttribNameString(-iLoop) into WRITEFORMAT output
!===================================================================================================================================
SUBROUTINE WriteDataHeaderInfo(unit_index,AttribName,OutputCounter,LoopSize,iLoop_in)
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
INTEGER,INTENT(IN),OPTIONAL :: LoopSize
INTEGER,INTENT(IN),OPTIONAL :: iLoop_in
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                     :: iLoop
!===================================================================================================================================
IF(PRESENT(LoopSize))THEN
  DO iLoop = 1, LoopSize
    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
    WRITE(unit_index,'(I3.3,A,A,A,I3.3)',ADVANCE='NO') OutputCounter,'-',TRIM(AttribName),'-',iLoop
    OutputCounter = OutputCounter + 1
  END DO
ELSE IF(PRESENT(iLoop_in)) THEN
  iLoop = iLoop_in
  WRITE(unit_index,'(A1)',ADVANCE='NO') ','
  WRITE(unit_index,'(I3.3,A,A,A,I3.3)',ADVANCE='NO') OutputCounter,'-',TRIM(AttribName),'-',iLoop
  OutputCounter = OutputCounter + 1
ELSE
  WRITE(unit_index,'(A1)',ADVANCE='NO') ','
  WRITE(unit_index,'(I3.3,A,A)',ADVANCE='NO') OutputCounter,'-',TRIM(AttribName)
  OutputCounter = OutputCounter + 1
END IF

END SUBROUTINE WriteDataHeaderInfo


!===================================================================================================================================
!> Writes input data into unit_index output: scalar or array+nVal can be given as input
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
INTEGER           ,INTENT(IN),OPTIONAL :: nVal
REAL              ,INTENT(IN),OPTIONAL :: RealScalar
INTEGER           ,INTENT(IN),OPTIONAL :: IntegerScalar
CHARACTER(LEN=*)  ,INTENT(IN),OPTIONAL :: StrScalar
REAL              ,INTENT(IN),OPTIONAL :: RealArray(:)
INTEGER           ,INTENT(IN),OPTIONAL :: IntegerArray(:)
INTEGER(KIND=8)   ,INTENT(IN),OPTIONAL :: IntegerK8Array(:)
CHARACTER(LEN=255),INTENT(IN),OPTIONAL :: StrArray(:)
LOGICAL           ,INTENT(IN),OPTIONAL :: LogicalScalar
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                     :: iLoop
!===================================================================================================================================

! Real
IF(PRESENT(RealArray))THEN
  IF(PRESENT(nVal)) THEN
    DO iLoop = 1, nVal
      WRITE (unit_index, CSVFORMAT, ADVANCE='NO') ',',RealArray(iLoop)
    END DO
  ELSE
    CALL abort(__STAMP__,'ERROR WriteDataInfo: nVal required for array output!')
  END IF
END IF

IF(PRESENT(RealScalar))THEN
  WRITE (unit_index, CSVFORMAT, ADVANCE='NO') ',',RealScalar
END IF

! Integer
IF(PRESENT(IntegerArray))THEN
  IF(PRESENT(nVal)) THEN
    DO iLoop = 1, nVal
      WRITE (unit_index, CSVFORMAT, ADVANCE='NO') ',',REAL(IntegerArray(iLoop))
    END DO
  ELSE
    CALL abort(__STAMP__,'ERROR WriteDataInfo: nVal required for array output!')
  END IF
END IF

IF(PRESENT(IntegerK8Array))THEN
  IF(PRESENT(nVal)) THEN
    DO iLoop = 1, nVal
      WRITE (unit_index, CSVFORMAT, ADVANCE='NO') ',',REAL(IntegerK8Array(iLoop))
    END DO
  ELSE
    CALL abort(__STAMP__,'ERROR WriteDataInfo: nVal required for array output!')
  END IF
END IF

IF(PRESENT(IntegerScalar))THEN
  WRITE (unit_index, CSVFORMAT, ADVANCE='NO') ',',REAL(IntegerScalar)
END IF

! String
IF(PRESENT(StrArray))THEN
  IF(PRESENT(nVal)) THEN
    DO iLoop = 1, nVal
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,'(A)',ADVANCE='NO') TRIM(StrArray(iLoop))
    END DO
  ELSE
    CALL abort(__STAMP__,'ERROR WriteDataInfo: nVal required for array output!')
  END IF
END IF

IF(PRESENT(StrScalar))THEN
  WRITE(unit_index,'(A1)',ADVANCE='NO') ','
  WRITE(unit_index,'(A)',ADVANCE='NO') TRIM(StrScalar)
END IF

! Logical
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
!> Communicate porous BC data across all ranks
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
  CALL MPI_REDUCE(MPI_IN_PLACE  , PorousBCOutput, 5*nPorousBC, MPI_DOUBLE_PRECISION, MPI_SUM, 0, PartMPI%COMM, iError)
ELSE
  CALL MPI_REDUCE(PorousBCOutput, PorousBCOutput, 5*nPorousBC, MPI_DOUBLE_PRECISION, MPI_SUM, 0, PartMPI%COMM, iError)
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
#if USE_MPI
USE MOD_Globals
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: BPO
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
#if USE_MPI
REAL    :: SendBuf(1:BPO%NPartBoundaries*BPO%NSpecies)
INTEGER :: SendBufSize
#endif /*USE_MPI*/
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
#if USE_MPI
USE MOD_Globals
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: SEE
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
  CALL MPI_REDUCE(MPI_IN_PLACE        , SEE%RealElectronOut, SEE%NPartBoundaries,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
ELSE
  CALL MPI_REDUCE(SEE%RealElectronOut , 0                  , SEE%NPartBoundaries,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
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
USE MOD_Globals                   ,ONLY: CollectiveStop,UNIT_stdOut
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: BPO
USE MOD_Particle_Boundary_Vars    ,ONLY: nPartBound,PartBound
USE MOD_ReadInTools               ,ONLY: GETLOGICAL,GETINT,GETINTARRAY
USE MOD_Analyze_Vars              ,ONLY: DoSurfModelAnalyze
USE MOD_Particle_Vars             ,ONLY: nSpecies,Species
#if USE_MPI
USE MOD_Globals                   ,ONLY: MPIRoot
#endif /*USE_MPI*/
#if USE_HDG
USE MOD_Globals                   ,ONLY: abort
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: SurfaceAnalyzeStep
USE MOD_Analyze_Vars              ,ONLY: CalcElectricTimeDerivative,FieldAnalyzeStep
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
#endif /*USE_HDG*/
USE MOD_Mesh_Vars                 ,ONLY: nBCs
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iPartBound,iSpec,iBPO,iBC
!===================================================================================================================================
DoSurfModelAnalyze = .TRUE.
BPO%NPartBoundaries = GETINT('BPO-NPartBoundaries')
BPO%PartBoundaries  = GETINTARRAY('BPO-PartBoundaries',BPO%NPartBoundaries)
BPO%NSpecies        = GETINT('BPO-NSpecies')
BPO%Species         = GETINTARRAY('BPO-Species',BPO%NSpecies)
IF(BPO%NPartBoundaries.EQ.0.OR.BPO%NSpecies.EQ.0) CALL CollectiveStop(__STAMP__,'BPO-NPartBoundaries or BPO-NSpecies is zero.')
ALLOCATE(BPO%RealPartOut(1:BPO%NPartBoundaries,1:BPO%NSpecies))
BPO%RealPartOut = 0.

! Initialize
BPO%OutputTotalElectricCurrent = .FALSE.

ALLOCATE(BPO%SpecIDToBPOSpecID(1:nSpecies))
BPO%SpecIDToBPOSpecID = -1
DO iSpec = 1, BPO%NSpecies
  ! Sanity check
  IF(BPO%Species(iSpec).GT.nSpecies) CALL CollectiveStop(__STAMP__,&
      'BPO-Species contains a wrong species ID, which is greater than nSpecies')
  BPO%SpecIDToBPOSpecID(BPO%Species(iSpec)) = iSpec
  ! Activate output of total electric current when at least one species has a charge unequal to zero
  IF(ABS(Species(BPO%Species(iSpec))%ChargeIC).GT.0.0) BPO%OutputTotalElectricCurrent = .TRUE.
END DO ! iSpec = 1, BPO%NSpecies

#if USE_HDG
!-- Electric displacement current: Sanity check
IF(BPO%OutputTotalElectricCurrent.AND.CalcElectricTimeDerivative)THEN
  IF(SurfaceAnalyzeStep.LT.FieldAnalyzeStep) CALL abort(__STAMP__,&
      'SurfaceAnalyzeStep < FieldAnalyzeStep is not allowed for CalcElectricTimeDerivative=T and CalcBoundaryParticleOutput=T')
  IF(MOD(SurfaceAnalyzeStep,FieldAnalyzeStep).NE.0) CALL abort(__STAMP__,&
      'MOD(SurfaceAnalyzeStep,FieldAnalyzeStep) must be zero for CalcElectricTimeDerivative=T and CalcBoundaryParticleOutput=T')
  ! Display warning if the two analyze step variables are not the same. They should be to have the same sampling time for both
  ! outputs when they are added together
  IF(SurfaceAnalyzeStep.NE.FieldAnalyzeStep)THEN
    LBWRITE (*,*) "WARNING: SurfaceAnalyzeStep should be equal to FieldAnalyzeStep for CalcElectricTimeDerivative=T and CalcBoundaryParticleOutput=T"
  END IF ! SurfaceAnalyzeStep.NE.FieldAnalyzeStep
END IF ! BPO%OutputTotalElectricCurrent.AND.CalcElectricTimeDerivative
#endif /*USE_HDG*/

ALLOCATE(BPO%BCIDToBPOBCID(1:nPartBound))
BPO%BCIDToBPOBCID = -1
ALLOCATE(BPO%FieldBoundaries(1:BPO%NPartBoundaries))
BPO%FieldBoundaries = -1
DO iBPO = 1, BPO%NPartBoundaries
  ! Get particle boundary ID
  iPartBound = BPO%PartBoundaries(iBPO)
  ! Fill mapping from iBPO to field BC index
  DO iBC = 1, nBCs
    ! Find matching iBC
    IF(PartBound%MapToPartBC(iBC).EQ.iPartBound)THEN
      BPO%FieldBoundaries(iBPO) = iBC
    END IF ! PartBound%MapToPartBC(1:nBCs).EQ.iPartBound
  END DO ! iBC = 1, nBCs
  ! Fill mapping from iPartBound to iBPO
  BPO%BCIDToBPOBCID(iPartBound) = iBPO
  ! Sanity check BC types: iPartBound = 1 (open) or 2 (ReflectiveBC)
  ! Add more BCs to the vector if required
  IF(.NOT.ANY(PartBound%TargetBoundCond(iPartBound).EQ.(/1,2/)))THEN
    ! Check if species swap is used
    IF(PartBound%NbrOfSpeciesSwaps(iPartBound).GT.0)THEN
      ! nothing to do for now
    ELSE
      SWRITE(UNIT_stdOut,'(A)')'\nError for CalcBoundaryParticleOutput=T\n'
      SWRITE(UNIT_stdOut,'(A,I0)')'  iPartBound = ',iPartBound
      SWRITE(UNIT_stdOut,'(A,A)') '  SourceName = ',TRIM(PartBound%SourceBoundName(iPartBound))
      SWRITE(UNIT_stdOut,'(A,I0)')'   Condition = ',PartBound%TargetBoundCond(iPartBound)
      SWRITE(UNIT_stdOut,'(A)')'\n  Conditions are'//&
                          '  OpenBC          = 1  \n'//&
          '                  ReflectiveBC    = 2  \n'//&
          '                  PeriodicBC      = 3  \n'//&
          '                  SimpleAnodeBC   = 4  \n'//&
          '                  SimpleCathodeBC = 5  \n'//&
          '                  RotPeriodicBC   = 6  \n'//&
          '                  SymmetryBC      = 10 \n'//&
          '                  SymmetryAxis    = 11 '
      CALL CollectiveStop(__STAMP__&
          ,'PartBound%TargetBoundCond(iPartBound) is not implemented for CalcBoundaryParticleOutput',&
          IntInfo=PartBound%TargetBoundCond(iPartBound))
    END IF ! PartBound%NbrOfSpeciesSwaps(iPartBound).GT.0
  END IF ! .NOT.ANY(PartBound%TargetBoundCond(iPartBound).EQ. ...

  ! Check reflective BCs with or without species swap activated
  IF(PartBound%TargetBoundCond(iPartBound).EQ.2)THEN
    ! Check if species swap is used
    IF(PartBound%NbrOfSpeciesSwaps(iPartBound).GT.0)THEN
      ! nothing to do for now
    ELSE
      ! Check the surface model
      SELECT CASE(PartBound%SurfaceModel(iPartBound))
      CASE(SEE_MODELS_ID)
        ! all secondary electron models
      CASE DEFAULT
        CALL CollectiveStop(__STAMP__,'CalcBoundaryParticleOutput not implemented for this '//&
        'PartBound%SurfaceModel(iPartBound). Either select different surface model or activate NbrOfSpeciesSwaps',&
            IntInfo=PartBound%SurfaceModel(iPartBound))
      END SELECT
    END IF ! PartBound%NbrOfSpeciesSwaps(iPartBound).GT.0
  END IF ! PartBound%TargetBoundCond(BPO%PartBoundaries(iPartBound).EQ.2)
END DO ! iPartBound = 1, BPO%NPartBoundaries

END SUBROUTINE InitBoundaryParticleOutput


!===================================================================================================================================
!> Allocate the required arrays (mappings and containers) for secondary electron emission analysis, which tracks the number of
!> electrons that are emitted from a surface
!>
!> 1) Check if secondary electron emission occurs
!> 1.1) Count number of different SEE boundaries via reflective particle BC
!> 1.2) Count number of different photon SEE boundaries
!> 2.) Create Mapping from SEE BC index to particle BC index
!> 2.1) Create mapping for reactive surface SEE (non-photon impacts)
!> 2.2) Create mapping for photon-surface SEE
!> 3) Create Mapping from particle BC index to SEE BC index
!===================================================================================================================================
SUBROUTINE InitCalcElectronSEE()
! MODULES
USE MOD_Globals
USE MOD_Globals                   ,ONLY: abort,StringBeginsWith
USE MOD_Analyze_Vars              ,ONLY: DoSurfModelAnalyze
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: SEE,CalcBoundaryParticleOutput,CalcElectronSEE
USE MOD_Particle_Boundary_Vars    ,ONLY: nPartBound,PartBound
USE MOD_ReadInTools               ,ONLY: GETLOGICAL,PrintOption
USE MOD_Particle_Vars             ,ONLY: Species,nSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iSEE,iPartBound,NPartBoundariesReflectiveSEE,NPartBoundariesPhotonSEE,iInit,iSpec
!===================================================================================================================================

! Initialize
CalcElectronSEE = .FALSE.

! 1) Check if secondary electron emission occurs
! 1.1) Count number of different SEE boundaries via reflective particle BC
SEE%NPartBoundaries = 0
DO iPartBound=1,nPartBound
  IF(.NOT.PartBound%Reactive(iPartBound)) CYCLE
  SELECT CASE(PartBound%SurfaceModel(iPartBound))
  CASE(SEE_MODELS_ID)
    SEE%NPartBoundaries = SEE%NPartBoundaries + 1
  END SELECT
END DO ! iPartBound=1,nPartBound
NPartBoundariesReflectiveSEE = SEE%NPartBoundaries

! 1.2) Count number of different photon SEE boundaries
DO iSpec=1,nSpecies
  DO iInit=1, Species(iSpec)%NumberOfInits
    IF(Species(iSpec)%Init(iInit)%PartBCIndex.GT.0)THEN
      IF(StringBeginsWith(Species(iSpec)%Init(iInit)%SpaceIC,'photon_SEE')) SEE%NPartBoundaries = SEE%NPartBoundaries + 1
    END IF
  END DO
END DO
NPartBoundariesPhotonSEE = SEE%NPartBoundaries - NPartBoundariesReflectiveSEE

! If not SEE boundaries exist, no measurement of the current can be performed
IF(SEE%NPartBoundaries.EQ.0) RETURN

! Automatically activate when CalcBoundaryParticleOutput=T
IF(CalcBoundaryParticleOutput)THEN
  CalcElectronSEE = .TRUE.
  CALL PrintOption('SEE current measurement activated (CalcBoundaryParticleOutput=T): CalcElectronSEE','INFO',&
      LogOpt=CalcElectronSEE)
ELSE
  CalcElectronSEE = GETLOGICAL('CalcElectronSEE','.FALSE.')
END IF ! CalcBoundaryParticleOutput

! Check if analysis has been activated
IF(.NOT.CalcElectronSEE) RETURN

! Automatically activate surface model analyze flag
DoSurfModelAnalyze = .TRUE.

! 2.) Create Mapping from SEE BC index to particle BC index
! 2.1) Create mapping for reactive surface SEE (non-photon impacts)
ALLOCATE(SEE%PartBoundaries(SEE%NPartBoundaries))
SEE%NPartBoundaries = 0
DO iPartBound=1,nPartBound
  IF(.NOT.PartBound%Reactive(iPartBound)) CYCLE
  SELECT CASE(PartBound%SurfaceModel(iPartBound))
  CASE(SEE_MODELS_ID)
    SEE%NPartBoundaries = SEE%NPartBoundaries + 1
    SEE%PartBoundaries(SEE%NPartBoundaries) = iPartBound
  END SELECT
END DO ! iPartBound=1,nPartBound
! 2.2) Create mapping for photon-surface SEE
DO iSpec=1,nSpecies
  DO iInit=1, Species(iSpec)%NumberOfInits
    IF(Species(iSpec)%Init(iInit)%PartBCIndex.GT.0)THEN
      IF(StringBeginsWith(Species(iSpec)%Init(iInit)%SpaceIC,'photon_SEE'))THEN
        SEE%NPartBoundaries = SEE%NPartBoundaries + 1
        SEE%PartBoundaries(SEE%NPartBoundaries) = Species(iSpec)%Init(iInit)%PartBCIndex
      END IF
    END IF
  END DO
END DO

! Allocate the container
ALLOCATE(SEE%RealElectronOut(1:SEE%NPartBoundaries))
SEE%RealElectronOut = 0.

! 3) Create Mapping from particle BC index to SEE BC index
ALLOCATE(SEE%BCIDToSEEBCID(1:nPartBound))
SEE%BCIDToSEEBCID     = -1
DO iSEE = 1, SEE%NPartBoundaries
  iPartBound = SEE%PartBoundaries(iSEE)
  SEE%BCIDToSEEBCID(iPartBound) = iSEE
END DO ! iSEE = 1, SEE%NPartBoundaries

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
  SDEALLOCATE(BPO%FieldBoundaries)
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
