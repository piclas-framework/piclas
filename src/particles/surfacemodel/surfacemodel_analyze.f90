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
CALL prms%CreateLogicalOption(  'CalcCurrentSEE'            , 'Count the electron emission from BCs where SEE is active','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcEnergyViolationSEE'    , 'Count the number and amount of energy conservation violations during SEE '//&
                                                              ' using the Chung-Everhart distribution for multiple secondaries','.FALSE.')

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
CALL InitCalcElectronSEE() ! This routine calls GETLOGICAL('CalcCurrentSEE','.FALSE.')

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
USE MOD_Particle_Boundary_Vars    ,ONLY: PartBound
USE MOD_SurfaceModel_Vars         ,ONLY: nPorousBC, PorousBC
USE MOD_Particle_Vars             ,ONLY: nSpecies,UseNeutralization,NeutralizationBalanceGlobal,Species,VarTimeStep
#if USE_MPI
USE MOD_Particle_Boundary_Vars    ,ONLY: SurfCOMM
#endif /*USE_MPI*/
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
IF(.NOT.DoSurfModelAnalyze) RETURN

! Only proceed with processors, which have a surface side (as determined in InitParticleBoundarySurfSides)
#if USE_MPI
IF(SurfCOMM%UNICATOR.EQ.MPI_COMM_NULL) RETURN
#endif /*USE_MPI*/

SurfModelAnalyzeSampleTime = Time - SurfModelAnalyzeSampleTime ! Set SurfModelAnalyzeSampleTime=Time at the end of this routine
OutputCounter = 2
unit_index = 636
IF(MPIRoot)THEN
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
      IF(CalcCurrentSEE)THEN
        DO iSEE = 1, SEE%NPartBoundaries
          ASSOCIATE(BoundaryNameSEE => PartBound%SourceBoundName(SEE%PartBoundaries(iSEE)))
            IF(CalcElectronSEE) THEN
              CALL WriteDataHeaderInfo(unit_index,'ElectricCurrentSEE-'//TRIM(BoundaryNameSEE),OutputCounter)
            END IF
            IF(CalcPhotonSEE) THEN
              CALL WriteDataHeaderInfo(unit_index,'ElectricCurrentPhotonSEE-'//TRIM(BoundaryNameSEE),OutputCounter)
            END IF
            IF(CalcEnergyViolationSEE) THEN
              CALL WriteDataHeaderInfo(unit_index,'EnergyViolationRatioSEE-'//TRIM(BoundaryNameSEE),OutputCounter)
              CALL WriteDataHeaderInfo(unit_index,'EnergyViolationAdditionSEE-'//TRIM(BoundaryNameSEE),OutputCounter)
            END IF
          END ASSOCIATE
        END DO ! iSEE = 1, SEE%NPartBoundaries
      END IF ! CalcCurrentSEE
      WRITE(unit_index,'(A)') ''
    END IF
  END IF
END IF

!===================================================================================================================================
! Analyze Routines
!===================================================================================================================================
IF (CalcSurfCollCounter)        CALL GetCollCounter(SurfCollNum,AdsorptionNum,DesorptionNum)
IF (CalcPorousBCInfo)           CALL GetPorousBCInfo()
#if USE_MPI
IF (CalcBoundaryParticleOutput) CALL SyncBoundaryParticleOutput()
IF (CalcCurrentSEE)             CALL SyncElectronSEE()
#endif /*USE_MPI*/
!===================================================================================================================================
! Output Analyzed variables
!===================================================================================================================================
#if USE_MPI
IF(MPIRoot)THEN
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
          ! Scaling the number of particles depending on the species-specific timestep factor
          IF(VarTimeStep%UseSpeciesSpecific) BPO%RealPartOut(iBPO,iSpec) = BPO%RealPartOut(iBPO,iSpec) &
                                                                            / Species(BPO%Species(iSpec))%TimeStepFactor
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
          IF(CalcCurrentSEE)THEN
            ! Get particle boundary ID
            iPartBound = BPO%PartBoundaries(iBPO)
            ! Get SEE ID
            iSEE = SEE%BCIDToSEEBCID(iPartBound)
            ! Skip boundaries without SEE
            IF(iSEE.GT.0) THEN
              ! Add SEE current if this BC has secondary electron emission
              IF(CalcElectronSEE) TotalElectricCharge = TotalElectricCharge + SEE%RealElectronOut(iSEE)
              ! Add SEE current if this BC has photon-based SEE
              IF(CalcPhotonSEE) TotalElectricCharge = TotalElectricCharge + SEE%RealElectronOutPhoton(iSEE)
            END IF
          END IF ! CalcCurrentSEE
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
        IF(CalcCurrentSEE)THEN
          ! Get particle boundary ID
          iPartBound2 = BPO%PartBoundaries(iBPO)
          ! Sanity Check
          IF(iPartBound.NE.iPartBound2) CALL abort(__STAMP__,'AnalyzeSurface(): Wrong particle boundary encountered!')
          ! Get SEE ID
          iSEE = SEE%BCIDToSEEBCID(iPartBound)
          ! Skip boundaries without SEE
          IF(iSEE.GT.0) THEN
            ! Add SEE current if this BC has secondary electron emission
            IF(CalcElectronSEE) TotalElectricCharge = TotalElectricCharge + SEE%RealElectronOut(iSEE)
            ! Add SEE current if this BC has photon-based SEE
            IF(CalcPhotonSEE) TotalElectricCharge = TotalElectricCharge + SEE%RealElectronOutPhoton(iSEE)
          END IF
        END IF ! CalcCurrentSEE
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
        ! Reset MPIRoot counters after writing the data to the file,
        ! non-MPIRoot are reset in SyncBoundaryParticleOutput()
        BPO%RealPartOut(iPartBound,iSpec) = 0.
      END DO ! iSpec = 1, BPO%NSpecies
    END DO ! iPartBound = 1, BPO%NPartBoundaries
  END IF ! CalcBoundaryParticleOutput

  IF(UseNeutralization) CALL WriteDataInfo(unit_index,RealScalar=REAL(NeutralizationBalanceGlobal))

  IF(CalcCurrentSEE)THEN
    DO iPartBound = 1, SEE%NPartBoundaries
      IF(CalcElectronSEE) THEN
        IF(ABS(SurfModelAnalyzeSampleTime).LE.0.0)THEN
          CALL WriteDataInfo(unit_index,RealScalar=0.0)
        ELSE
          CALL WriteDataInfo(unit_index,RealScalar=SEE%RealElectronOut(iPartBound)/SurfModelAnalyzeSampleTime)
        END IF ! ABS(SurfModelAnalyzeSampleTime).LE.0.0
        ! Reset MPIRoot counters after writing the data to the file, non-MPIRoot are reset in SyncBoundaryParticleOutput()
        SEE%RealElectronOut(iPartBound) = 0.
      END IF
      IF(CalcPhotonSEE) THEN
        IF(ABS(SurfModelAnalyzeSampleTime).LE.0.0)THEN
          CALL WriteDataInfo(unit_index,RealScalar=0.0)
        ELSE
          CALL WriteDataInfo(unit_index,RealScalar=SEE%RealElectronOutPhoton(iPartBound)/SurfModelAnalyzeSampleTime)
        END IF ! ABS(SurfModelAnalyzeSampleTime).LE.0.0
        ! Reset MPIRoot counters after writing the data to the file, non-MPIRoot are reset in SyncBoundaryParticleOutput()
        SEE%RealElectronOutPhoton(iPartBound) = 0.
      END IF
      IF(CalcEnergyViolationSEE) THEN
        IF(ABS(SEE%EventCount(iPartBound)).LE.0.0)THEN
          CALL WriteDataInfo(unit_index,RealScalar=0.0)
          CALL WriteDataInfo(unit_index,RealScalar=0.0)
        ELSE
          CALL WriteDataInfo(unit_index,RealScalar=SEE%EnergyConsViolationCount(iPartBound)/SEE%EventCount(iPartBound))
          CALL WriteDataInfo(unit_index,RealScalar=SEE%EnergyConsViolationSum(iPartBound)/SEE%EventCount(iPartBound))
        END IF ! ABS(SEE%EventCount(iPartBound)).LE.0.0
        ! Reset MPIRoot counters after writing the data to the file, non-MPIRoot are reset in SyncBoundaryParticleOutput()
        SEE%EventCount(iPartBound) = 0.
        SEE%EnergyConsViolationCount(iPartBound) = 0.
        SEE%EnergyConsViolationSum(iPartBound) = 0.
      END IF
    END DO ! iPartBound = 1, SEE%NPartBoundaries
  END IF ! CalcCurrentSEE
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
USE MOD_Particle_Boundary_Vars    ,ONLY: SurfCOMM
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
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,SurfCollNum ,nSpecies,MPI_INTEGER,MPI_SUM,0,SurfCOMM%UNICATOR,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,AdsorbNum   ,nSpecies,MPI_INTEGER,MPI_SUM,0,SurfCOMM%UNICATOR,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,DesorbNum   ,nSpecies,MPI_INTEGER,MPI_SUM,0,SurfCOMM%UNICATOR,IERROR)
ELSE
  CALL MPI_REDUCE(SurfCollNum ,SurfCollNum ,nSpecies,MPI_INTEGER,MPI_SUM,0,SurfCOMM%UNICATOR,IERROR)
  CALL MPI_REDUCE(AdsorbNum   ,AdsorbNum   ,nSpecies,MPI_INTEGER,MPI_SUM,0,SurfCOMM%UNICATOR,IERROR)
  CALL MPI_REDUCE(DesorbNum   ,DesorbNum   ,nSpecies,MPI_INTEGER,MPI_SUM,0,SurfCOMM%UNICATOR,IERROR)
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
USE MOD_Particle_Boundary_Vars    ,ONLY: SurfCOMM
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
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE  , PorousBCOutput, 5*nPorousBC, MPI_DOUBLE_PRECISION, MPI_SUM, 0, SurfCOMM%UNICATOR, iError)
ELSE
  CALL MPI_REDUCE(PorousBCOutput, PorousBCOutput, 5*nPorousBC, MPI_DOUBLE_PRECISION, MPI_SUM, 0, SurfCOMM%UNICATOR, iError)
END IF
#endif /*USE_MPI*/

IF(MPIRoot)THEN
  DO iPBC = 1, nPorousBC
    IF(PorousBCOutput(1,iPBC).GT.0.0)THEN
      ! Pumping Speed (Output(2)) is the sum of all elements (counter over particles exiting through pump)
      ! Other variables are averaged over the elements
      PorousBCOutput(3:5,iPBC) = PorousBCOutput(3:5,iPBC) / PorousBCOutput(1,iPBC)
    END IF
  END DO
END IF

END SUBROUTINE GetPorousBCInfo


!===================================================================================================================================
!> Synchronize BoundaryParticleOutput analyze arrays
!===================================================================================================================================
#if USE_MPI
SUBROUTINE SyncBoundaryParticleOutput()
! MODULES
USE MOD_Globals
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: BPO
USE MOD_Particle_Boundary_Vars    ,ONLY: SurfCOMM
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
SendBufSize = BPO%NPartBoundaries*BPO%NSpecies

! Map 2D array to vector for sending via MPI
SendBuf = RESHAPE(BPO%RealPartOut,(/SendBufSize/))
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,SendBuf(1:SendBufSize),SendBufSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%UNICATOR,IERROR)
  ! MAP vector back to 2D array
  BPO%RealPartOut = RESHAPE(SendBuf,(/BPO%NPartBoundaries,BPO%NSpecies/))
ELSE
  CALL MPI_REDUCE(SendBuf(1:SendBufSize),0,SendBufSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%UNICATOR,IERROR)
  ! Reset non SurfCOMM%UNICATOR counters, SurfCOMM%UNICATOR counters are reset after writing the data to the file
  BPO%RealPartOut = 0.
END IF

END SUBROUTINE SyncBoundaryParticleOutput
#endif /*USE_MPI*/


!===================================================================================================================================
!> Synchronize CalcCurrentSEE analyze arrays
!===================================================================================================================================
#if USE_MPI
SUBROUTINE SyncElectronSEE()
! MODULES
USE MOD_Globals
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: SEE, CalcElectronSEE, CalcPhotonSEE, CalcEnergyViolationSEE
USE MOD_Particle_Boundary_Vars    ,ONLY: SurfCOMM
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF (MPIRoot) THEN
  IF(CalcElectronSEE) CALL MPI_REDUCE(MPI_IN_PLACE, SEE%RealElectronOut         , SEE%NPartBoundaries,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%UNICATOR,IERROR)
  IF(CalcPhotonSEE) CALL MPI_REDUCE(MPI_IN_PLACE  , SEE%RealElectronOutPhoton   , SEE%NPartBoundaries,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%UNICATOR,IERROR)
  IF(CalcEnergyViolationSEE) THEN
    CALL MPI_REDUCE(MPI_IN_PLACE                  , SEE%EventCount              , SEE%NPartBoundaries,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%UNICATOR,IERROR)
    CALL MPI_REDUCE(MPI_IN_PLACE                  , SEE%EnergyConsViolationCount, SEE%NPartBoundaries,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%UNICATOR,IERROR)
    CALL MPI_REDUCE(MPI_IN_PLACE                  , SEE%EnergyConsViolationSum  , SEE%NPartBoundaries,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%UNICATOR,IERROR)
  END IF
ELSE
  ! Reset non MPIRoot counters, MPIRoot counters are reset after writing the data to the file
  IF(CalcElectronSEE) THEN
    CALL MPI_REDUCE(SEE%RealElectronOut     , 0             , SEE%NPartBoundaries,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%UNICATOR,IERROR)
    SEE%RealElectronOut = 0.
  END IF
  IF(CalcPhotonSEE) THEN
    CALL MPI_REDUCE(SEE%RealElectronOutPhoton , 0             , SEE%NPartBoundaries,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%UNICATOR,IERROR)
    SEE%RealElectronOutPhoton = 0.
  END IF
  IF(CalcEnergyViolationSEE) THEN
    CALL MPI_REDUCE(SEE%EventCount                            , 0             , SEE%NPartBoundaries,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%UNICATOR,IERROR)
    CALL MPI_REDUCE(SEE%EnergyConsViolationCount              , 0             , SEE%NPartBoundaries,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%UNICATOR,IERROR)
    CALL MPI_REDUCE(SEE%EnergyConsViolationSum                , 0             , SEE%NPartBoundaries,MPI_DOUBLE_PRECISION,MPI_SUM,0,SurfCOMM%UNICATOR,IERROR)
    ! Reset non MPIRoot counters, MPIRoot counters are reset after writing the data to the file
    SEE%EventCount = 0.
    SEE%EnergyConsViolationCount = 0.
    SEE%EnergyConsViolationSum = 0.
  END IF
END IF

END SUBROUTINE SyncElectronSEE
#endif /*USE_MPI*/


!===================================================================================================================================
!> Allocate the required arrays (mappings and containers) for BoundaryParticleOutput
!===================================================================================================================================
SUBROUTINE InitBoundaryParticleOutput()
! MODULES
USE MOD_Globals
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: BPO
USE MOD_Particle_Boundary_Vars    ,ONLY: nPartBound,PartBound
USE MOD_ReadInTools               ,ONLY: GETLOGICAL,GETINT,GETINTARRAY
USE MOD_Analyze_Vars              ,ONLY: DoSurfModelAnalyze
USE MOD_Particle_Vars             ,ONLY: nSpecies,Species
#if USE_MPI
USE MOD_Globals                   ,ONLY: MPIRoot
USE MOD_Mesh_Vars                 ,ONLY: BoundaryName
#endif /*USE_MPI*/
#if USE_HDG
USE MOD_Globals                   ,ONLY: abort
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: SurfaceAnalyzeStep
USE MOD_Analyze_Vars              ,ONLY: CalcElectricTimeDerivative,FieldAnalyzeStep
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
#endif /*USE_HDG*/
USE MOD_Mesh_Vars                 ,ONLY: nBCs,BC,nSides
USE MOD_TimeDisc_Vars             ,ONLY: iter
USE MOD_Mesh_Vars                 ,ONLY: N_SurfMesh
USE MOD_Interpolation_Vars        ,ONLY: N_Inter
USE MOD_DG_Vars                   ,ONLY: DG_Elems_master,DG_Elems_slave
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iPartBound,iSpec,iBPO,iBC,BCSideID,Nloc
REAL,ALLOCATABLE  :: TotalSurfArea(:)
CHARACTER(LEN=32) :: hilf
!REAL              :: area
!===================================================================================================================================
DoSurfModelAnalyze  = .TRUE.
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
      SWRITE(UNIT_stdOut,'(A)')'\n  Conditions and availability are'//&
                          '  OpenBC                    = 1  (yes)\n'//&
          '                  ReflectiveBC              = 2  (yes)\n'//&
          '                  PeriodicBC                = 3  (no)\n'//&
          '                  RotPeriodicBC             = 6  (no)\n'//&
          '                  RotPeriodicInterPlaneBC   = 6  (no)\n'//&
          '                  SymmetryBC                = 10 (no)\n'//&
          '                  SymmetryAxis              = 11 (no)'
      CALL CollectiveStop(__STAMP__,'PartBound%TargetBoundCond(iPartBound) is not implemented for CalcBoundaryParticleOutput',&
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
      CASE(2)
        ! Event probability model
      CASE(SEE_MODELS_ID)
        ! all secondary electron models
      CASE (VDL_MODEL_ID)!,SEE_VDL_MODEL_ID)
        ! VDL - virtual dielectric layer
      CASE DEFAULT
        WRITE(UNIT=hilf,FMT='(I0)') iPartBound
        CALL CollectiveStop(__STAMP__,'CalcBoundaryParticleOutput not implemented for '//&
        'PartBound%SurfaceModel(iPartBound='//TRIM(hilf)//'). Either select different surface model or activate NbrOfSpeciesSwaps',&
            IntInfo=PartBound%SurfaceModel(iPartBound))
      END SELECT
    END IF ! PartBound%NbrOfSpeciesSwaps(iPartBound).GT.0
  END IF ! PartBound%TargetBoundCond(BPO%PartBoundaries(iPartBound).EQ.2)
END DO ! iPartBound = 1, BPO%NPartBoundaries

! Display the total area of the all BPO%NPartBoundaries. MPI requires all-reduce to root process
IF(iter.EQ.0)THEN ! First iteration: Only output this information once
  ALLOCATE(TotalSurfArea(1:BPO%NPartBoundaries))
  TotalSurfArea = 0.
  ! Loop over all sides (to include inner BCs, which are not part of nBCSides) and get surface area
  DO BCSideID = 1,nSides
    IF(BC(BCSideID).EQ.0) CYCLE
    ! Get particle boundary ID
    iPartBound = PartBound%MapToPartBC(BC(BCSideID))
    ! get BPO boundary ID
    iBPO = BPO%BCIDToBPOBCID(iPartBound)
    ! Check if this boundary is tracked
    IF(iBPO.GT.0)THEN
      Nloc = MAX(DG_Elems_master(BCSideID),DG_Elems_slave(BCSideID))
      TotalSurfArea(iBPO) = TotalSurfArea(iBPO) + SUM(N_SurfMesh(BCSideID)%SurfElem(:,:)*N_Inter(Nloc)%wGPSurf(:,:))
    END IF ! iBPO.GT.0
  END DO ! BCSideID = 1,nSides
#if USE_MPI
  IF(MPIroot)THEN
    CALL MPI_REDUCE(MPI_IN_PLACE , TotalSurfArea , BPO%NPartBoundaries , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)
    WRITE(UNIT_stdOut,'(A)') "Total area used for BoundaryParticleOutput (BPO):"
    DO iBPO = 1, BPO%NPartBoundaries
      IF(iBPO.GT.9)THEN
        WRITE(UNIT_stdOut,'(A,I0,A,ES15.7,A)') "BPO-",iBPO,": ",TotalSurfArea(iBPO)," "//TRIM(BoundaryName(BPO%FieldBoundaries(iBPO)))
      ELSE
        WRITE(UNIT_stdOut,'(A,I0,A,ES16.7,A)') "BPO-",iBPO,": ",TotalSurfArea(iBPO)," "//TRIM(BoundaryName(BPO%FieldBoundaries(iBPO)))
      END IF ! iBPO.GT.9
    END DO ! iBPO = 1, BPO%NPartBoundaries
  ELSE
    CALL MPI_REDUCE(TotalSurfArea  , 0           , BPO%NPartBoundaries , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)
  END IF
#endif /*USE_MPI*/
  DEALLOCATE(TotalSurfArea)
END IF ! iter.EQ.0

END SUBROUTINE InitBoundaryParticleOutput


!===================================================================================================================================
!> Allocate the required arrays (mappings and containers) for secondary electron emission analysis, which tracks the number of
!> electrons that are emitted from a surface
!>
!> 1) Check if secondary electron emission occurs
!> 1.1) Count number of different SEE boundaries via reflective particle BC
!> 1.2) Count number of different photon SEE boundaries
!> 1.3) Count number of ray tracing photon SEE boundaries. WARNING: The combination of Init-based photon SEE and ray tracing photon
!       SEE is not allowed
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
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: SEE,CalcBoundaryParticleOutput,CalcCurrentSEE,CalcElectronSEE,CalcPhotonSEE,CalcEnergyViolationSEE
USE MOD_Particle_Boundary_Vars    ,ONLY: nPartBound,PartBound
USE MOD_ReadInTools               ,ONLY: GETLOGICAL,PrintOption
USE MOD_Particle_Vars             ,ONLY: Species,nSpecies
USE MOD_SurfaceModel_Vars         ,ONLY: SurfModEnergyDistribution
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iSEE,iPartBound,iInit,iSpec, PartBoundID
INTEGER           :: NPartBoundariesReflectiveSEE,NPartBoundariesPhotonSEE
LOGICAL           :: IsMappedBC(1:nPartBound)
!===================================================================================================================================

! Initialize
CalcElectronSEE = .FALSE.
CalcPhotonSEE   = .FALSE.
IsMappedBC      = .FALSE.
NPartBoundariesPhotonSEE = 0

! 1) Check if secondary electron emission occurs
! 1.1) Count number of different SEE boundaries via reflective particle BC
SEE%NPartBoundaries = 0
DO iPartBound=1,nPartBound
  IF(.NOT.PartBound%Reactive(iPartBound)) CYCLE
  SELECT CASE(PartBound%SurfaceModel(iPartBound))
  CASE(SEE_MODELS_ID)
    SEE%NPartBoundaries = SEE%NPartBoundaries + 1
    ! Mark BC as mapped for SEE
    IsMappedBC(iPartBound) = .TRUE.
  END SELECT
END DO ! iPartBound=1,nPartBound
NPartBoundariesReflectiveSEE = SEE%NPartBoundaries

! 1.2) Count number of different photon SEE boundaries
DO iSpec=1,nSpecies
  DO iInit=1, Species(iSpec)%NumberOfInits
    PartBoundID = Species(iSpec)%Init(iInit)%PartBCIndex
    IF(PartBoundID.GT.0)THEN
      IF(StringBeginsWith(Species(iSpec)%Init(iInit)%SpaceIC,'photon_SEE')) THEN
        NPartBoundariesPhotonSEE = NPartBoundariesPhotonSEE + 1
        ! Check if BC has already been flagged for regular SEE to avoid double count
        IF(.NOT.IsMappedBC(PartBoundID)) THEN
          SEE%NPartBoundaries = SEE%NPartBoundaries + 1
          IsMappedBC(PartBoundID) = .TRUE.
        END IF
      END IF
    END IF
  END DO
END DO

! Sanity check: The combination of Init-based photon SEE and ray tracing photon SEE is not allowed
IF(NPartBoundariesPhotonSEE.GT.0) THEN
  IF(ANY(PartBound%PhotonSEEYield(:).GT.0.0)) CALL abort(__STAMP__,'The combination of Init-based + ray tracing photon SEE is not allowed!')
END IF

! 1.3) Count number of ray tracing photon SEE boundaries
DO iPartBound = 1, nPartBound
  IF(PartBound%PhotonSEEYield(iPartBound).GT.0.0)THEN
    ! Remove the following check when the model is implemented (deposited charge holes by SEE)
    IF(PartBound%Dielectric(iPartBound)) CALL abort(__STAMP__,'Dielectric surfaces and ray tracing ist not implemented')
    NPartBoundariesPhotonSEE = NPartBoundariesPhotonSEE + 1
    ! Check if BC has already been flagged for regular SEE to avoid double count
    IF(.NOT.IsMappedBC(iPartBound)) THEN
      SEE%NPartBoundaries = SEE%NPartBoundaries + 1
      IsMappedBC(iPartBound) = .TRUE.
    END IF
  END IF ! PartBound%PhotonSEEYield(iPartBound).GT.0.0
END DO ! iPartBound = 1, nPartBound

! If not SEE boundaries exist, no measurement of the current can be performed
IF(SEE%NPartBoundaries.EQ.0) RETURN

! Automatically activate when CalcBoundaryParticleOutput=T
IF(CalcBoundaryParticleOutput)THEN
  CalcCurrentSEE = .TRUE.
  CALL PrintOption('SEE current activated (CalcBoundaryParticleOutput=T): CalcCurrentSEE','INFO',&
      LogOpt=CalcCurrentSEE)
ELSE
  CalcCurrentSEE = GETLOGICAL('CalcCurrentSEE','.FALSE.')
END IF ! CalcBoundaryParticleOutput

CalcEnergyViolationSEE = GETLOGICAL('CalcEnergyViolationSEE')

IF(CalcEnergyViolationSEE) THEN
  CalcEnergyViolationSEE = .FALSE.
  ! Sanity check: Energy violation output only with Chung-Everhart distribution
  DO iPartBound = 1, nPartBound
    IF(TRIM(SurfModEnergyDistribution(iPartBound)).EQ.'Chung-Everhart-cosine') CalcEnergyViolationSEE = .TRUE.
  END DO
  IF(CalcEnergyViolationSEE) THEN
    ! Automatically activate CalcCurrentSEE when CalcEnergyViolationSEE=T
    IF(.NOT.CalcCurrentSEE) THEN
      CalcCurrentSEE = .TRUE.
      CALL PrintOption('SEE current activated (CalcEnergyViolationSEE=T): CalcCurrentSEE','INFO',LogOpt=CalcCurrentSEE)
    END IF
  ELSE
    ! Deactivate CalcEnergyViolationSEE when no SEE model uses the Chung-Everhart distribution
    CALL PrintOption('SEE energy violation output deactivated, no SEE using the Chung-Everhart distribution: CalcEnergyViolationSEE',&
                      'INFO',LogOpt=CalcEnergyViolationSEE)
  END IF
END IF ! CalcBoundaryParticleOutput

! Check if analysis has been activated
IF(.NOT.CalcCurrentSEE) RETURN

! Automatically activate surface model analyze flag
DoSurfModelAnalyze = .TRUE.

! 2.) Create Mapping from SEE BC index to particle BC index
! 2.1) Create mapping for reactive surface SEE (non-photon impacts)
ALLOCATE(SEE%PartBoundaries(SEE%NPartBoundaries))
SEE%NPartBoundaries = 0
IsMappedBC          = .FALSE.
DO iPartBound=1,nPartBound
  IF(.NOT.PartBound%Reactive(iPartBound)) CYCLE
  SELECT CASE(PartBound%SurfaceModel(iPartBound))
  CASE(SEE_MODELS_ID)
    SEE%NPartBoundaries = SEE%NPartBoundaries + 1
    SEE%PartBoundaries(SEE%NPartBoundaries) = iPartBound
    IsMappedBC(iPartBound) = .TRUE.
  END SELECT
END DO ! iPartBound=1,nPartBound
! 2.2) Create mapping for photon-surface SEE (Init-based)
DO iSpec=1,nSpecies
  DO iInit=1, Species(iSpec)%NumberOfInits
    PartBoundID = Species(iSpec)%Init(iInit)%PartBCIndex
    IF(PartBoundID.GT.0)THEN
      IF(StringBeginsWith(Species(iSpec)%Init(iInit)%SpaceIC,'photon_SEE'))THEN
        IF(.NOT.IsMappedBC(PartBoundID)) THEN
          SEE%NPartBoundaries = SEE%NPartBoundaries + 1
          SEE%PartBoundaries(SEE%NPartBoundaries) = PartBoundID
          IsMappedBC(PartBoundID) = .TRUE.
        END IF
      END IF
    END IF
  END DO
END DO
! 2.3) Create mapping for photon-surface SEE (ray tracing)
DO iPartBound = 1, nPartBound
  IF(PartBound%PhotonSEEYield(iPartBound).GT.0.0)THEN
    IF(.NOT.IsMappedBC(iPartBound)) THEN
      SEE%NPartBoundaries = SEE%NPartBoundaries + 1
      SEE%PartBoundaries(SEE%NPartBoundaries) = iPartBound
      IsMappedBC(iPartBound) = .TRUE.
    END IF
  END IF ! PartBound%PhotonSEEYield(iPartBound).GT.0.0
END DO ! iPartBound = 1, nPartBound

! Allocate the container
IF(NPartBoundariesReflectiveSEE.GT.0) THEN
  CalcElectronSEE = .TRUE.
  ALLOCATE(SEE%RealElectronOut(1:SEE%NPartBoundaries))
  SEE%RealElectronOut = 0.
END IF

IF(NPartBoundariesPhotonSEE.GT.0) THEN
  CalcPhotonSEE = .TRUE.
  ALLOCATE(SEE%RealElectronOutPhoton(1:SEE%NPartBoundaries))
  SEE%RealElectronOutPhoton = 0.
END IF

IF(CalcEnergyViolationSEE) THEN
  ALLOCATE(SEE%EventCount(1:SEE%NPartBoundaries))
  SEE%EventCount = 0.
  ALLOCATE(SEE%EnergyConsViolationCount(1:SEE%NPartBoundaries))
  SEE%EnergyConsViolationCount = 0.
  ALLOCATE(SEE%EnergyConsViolationSum(1:SEE%NPartBoundaries))
  SEE%EnergyConsViolationSum = 0.
END IF

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
IF(CalcCurrentSEE)THEN
  SDEALLOCATE(SEE%RealElectronOut)
  SDEALLOCATE(SEE%RealElectronOutPhoton)
  SDEALLOCATE(SEE%EventCount)
  SDEALLOCATE(SEE%EnergyConsViolationCount)
  SDEALLOCATE(SEE%EnergyConsViolationSum)
  SDEALLOCATE(SEE%PartBoundaries)
  SDEALLOCATE(SEE%BCIDToSEEBCID)
END IF ! CalcCurrentSEE

END SUBROUTINE FinalizeSurfaceModelAnalyze

#endif /*PARTICLES*/


END MODULE MOD_SurfaceModel_Analyze
