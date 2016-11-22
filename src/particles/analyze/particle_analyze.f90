#include "boltzplatz.h"

MODULE MOD_Particle_Analyze
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
#ifdef PARTICLES
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitParticleAnalyze
  MODULE PROCEDURE InitParticleAnalyze
END INTERFACE

INTERFACE FinalizeParticleAnalyze
  MODULE PROCEDURE FinalizeParticleAnalyze
END INTERFACE

INTERFACE AnalyzeParticles
  MODULE PROCEDURE AnalyzeParticles
END INTERFACE

INTERFACE CalcKineticEnergy
  MODULE PROCEDURE CalcKineticEnergy
END INTERFACE

INTERFACE CalcShapeEfficiencyR
  MODULE PROCEDURE CalcShapeEfficiencyR
END INTERFACE

INTERFACE CalcEkinPart
  MODULE PROCEDURE CalcEkinPart
END INTERFACE

PUBLIC:: InitParticleAnalyze, FinalizeParticleAnalyze!, CalcPotentialEnergy
PUBLIC:: CalcKineticEnergy, CalcEkinPart, AnalyzeParticles
#if (PP_TimeDiscMethod == 42)
PUBLIC :: ElectronicTransition, WriteEletronicTransition
#endif
!===================================================================================================================================

CONTAINS

SUBROUTINE InitParticleAnalyze()
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_Preproc
  USE MOD_Analyze_Vars            ,ONLY: DoAnalyze
  USE MOD_Particle_Analyze_Vars  !,ONLY:ParticleAnalyzeInitIsDone, CalcCharge, CalcEkin, CalcEpot 
  USE MOD_ReadInTools             ,ONLY: GETLOGICAL, GETINT, GETSTR, GETINTARRAY, GETREALARRAY, GETREAL
  USE MOD_Particle_Vars           ,ONLY: nSpecies
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER   :: dir, VeloDirs_hilf(4)
!===================================================================================================================================
  IF (ParticleAnalyzeInitIsDone) THEN
CALL abort(__STAMP__,&
'InitParticleAnalyse already called.',999,999.)
    RETURN
  END IF
  SWRITE(UNIT_StdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE ANALYZE...'

  PartAnalyzeStep = GETINT('Part-AnalyzeStep','1')
  IF (PartAnalyzeStep.EQ.0) PartAnalyzeStep = 123456789

  DoAnalyze = .FALSE.
  CalcEpot = GETLOGICAL('CalcPotentialEnergy','.FALSE.')
  IF(CalcEpot) DoAnalyze = .TRUE.
  DoVerifyCharge = GETLOGICAL('PIC-VerifyCharge','.FALSE.')
  CalcCharge = GETLOGICAL('CalcCharge','.FALSE.')
  IF(CalcCharge) DoAnalyze = .TRUE. 
  CalcEkin = GETLOGICAL('CalcKineticEnergy','.FALSE.')
  CalcEint = GETLOGICAL('CalcInternalEnergy','.FALSE.')
  CalcTemp = GETLOGICAL('CalcTemp','.FALSE.')
  IF(CalcTemp.OR.CalcEint) DoAnalyze = .TRUE.
  IF(CalcEkin) THEN
    DoAnalyze = .TRUE.
    IF (nSpecies .GT. 1) THEN
      nEkin = nSpecies + 1
    ELSE
     nEkin = 1
    END IF
  END IF
  CalcPartBalance = GETLOGICAL('CalcPartBalance','.FALSE.')
  IF (CalcPartBalance) THEN
    DoAnalyze = .TRUE.
    SDEALLOCATE(nPartIn)
    SDEALLOCATE(nPartOut)
    SDEALLOCATE(PartEkinIn)
    SDEALLOCATE(PartEkinOut)
    ALLOCATE( nPartIn(nSpecies)     &
            , nPartOut(nSpecies)    &
            , PartEkinOut(nSpecies) &
            , PartEkinIn(nSpecies)  )
    nPartIn=0
    nPartOut=0
    PartEkinOut=0.
    PartEkinIn=0.
#if defined(LSERK) || defined(IMEX) || defined(IMPA) 
    SDEALLOCATE( nPartInTmp)
    SDEALLOCATE( PartEkinInTmp)
    ALLOCATE( nPartInTmp(nSpecies)     &
            , PartEkinInTmp(nSpecies)  )
    PartEkinInTmp=0.
    nPartInTmp=0
#endif
  END IF
  CalcVelos = GETLOGICAL('CalcVelos','.FALSE')
  IF (CalcVelos) THEN
    DoAnalyze=.TRUE.
    VeloDirs_hilf = GetIntArray('VelocityDirections',4,'1,1,1,1') ! x,y,z,abs -> 0/1 = T/F
    VeloDirs(:) = .FALSE.
    DO dir = 1,4
      IF (VeloDirs_hilf(dir) .EQ. 1) THEN
        VeloDirs(dir) = .TRUE.
      END IF
    END DO
    IF ((.NOT. VeloDirs(1)) .AND. (.NOT. VeloDirs(2)) .AND. &
        (.NOT. VeloDirs(3)) .AND. (.NOT. VeloDirs(4))) THEN
      CALL abort(&
        __STAMP__&
        ,'No VelocityDirections set in CalcVelos!')
    END IF
  END IF
  TrackParticlePosition = GETLOGICAL('Part-TrackPosition','.FALSE.')
  IF(TrackParticlePosition)THEN
    printDiff=GETLOGICAL('printDiff','.FALSE.')
    IF(printDiff)THEN
      printDiffTime=GETREAL('printDiffTime','12.')
      printDiffVec=GETREALARRAY('printDiffVec',6,'0.,0.,0.,0.,0.,0.')
    END IF
  END IF
  CalcNumSpec   = GETLOGICAL('CalcNumSpec','.FALSE.')
  CalcCollRates = GETLOGICAL('CalcCollRates','.FALSE.')
  CalcReacRates = GETLOGICAL('CalcReacRates','.FALSE.')
  IF(CalcNumSpec.OR.CalcCollRates.OR.CalcReacRates) DoAnalyze = .TRUE.
  CalcShapeEfficiency = GETLOGICAL('CalcShapeEfficiency','.FALSE.')
  IF (CalcShapeEfficiency) THEN
    DoAnalyze = .TRUE.
    CalcShapeEfficiencyMethod = GETSTR('CalcShapeEfficiencyMethod','AllParts')
    SELECT CASE(CalcShapeEfficiencyMethod)
    CASE('AllParts')  ! All currently available Particles are used
    CASE('SomeParts') ! A certain percentage of currently available Particles is used
      ShapeEfficiencyNumber = GETINT('ShapeEfficiencyNumber','100')  ! in percent
    CASE DEFAULT
      CALL abort(&
          __STAMP__&
          , ' CalcShapeEfficiencyMethod not implemented: ')

    END SELECT
  END IF

  IsRestart = GETLOGICAL('IsRestart','.FALSE.')

  ParticleAnalyzeInitIsDone=.TRUE.

  SWRITE(UNIT_stdOut,'(A)')' INIT PARTCILE ANALYZE DONE!'
  SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitParticleAnalyze


SUBROUTINE AnalyzeParticles(Time)
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
! MPI-INFO:
! important change: all MPI-communication is done directly in the corresponding subroutines
! a) its easier and cleaner
! b) reduces the probability of errors (routines are again fixed for MPI...)
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_Preproc
  USE MOD_Analyze_Vars,          ONLY: DoAnalyze
  USE MOD_Particle_Analyze_Vars!,ONLY: ParticleAnalyzeInitIsDone,CalcCharge,CalcEkin,IsRestart
  USE MOD_PARTICLE_Vars,         ONLY: nSpecies, BoltzmannConst
  USE MOD_DSMC_Vars,             ONLY: CollInf, useDSMC, CollisMode, ChemReac
  USE MOD_Restart_Vars,          ONLY: DoRestart
  USE MOD_AnalyzeField,          ONLY: CalcPotentialEnergy
  USE MOD_DSMC_Vars,             ONLY: DSMC
#if (PP_TimeDiscMethod==2 || PP_TimeDiscMethod==4 || PP_TimeDiscMethod==42 || (PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506))
  USE MOD_TimeDisc_Vars,         ONLY : iter, dt
#endif
  USE MOD_PIC_Analyze,           ONLY: CalcDepositedCharge
#ifdef MPI
  USE MOD_LoadBalance_Vars,      ONLY: tCurrent
  USE MOD_Particle_MPI_Vars,     ONLY: PartMPI
#endif /*MPI*/
#if ( PP_TimeDiscMethod ==42)
  USE MOD_DSMC_Vars,             ONLY: Adsorption,BGGas, SpecDSMC
  USE MOD_Particle_Vars,         ONLY: Species
  USE MOD_Particle_Mesh_Vars,    ONLY : GEO
#endif
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL,INTENT(IN)                 :: Time
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  LOGICAL             :: isOpen, FileExists
  CHARACTER(LEN=350)  :: outfile
  INTEGER             :: unit_index, iSpec, OutputCounter
  INTEGER(KIND=8)     :: SimNumSpec(nSpecies+1)
  REAL                :: WEl, WMag, NumSpec(nSpecies+1)
  REAL                :: Ekin(nSpecies + 1), Temp(nSpecies+1)
  REAL                :: IntEn(nSpecies,3),IntTemp(nSpecies,3),TempTotal(nSpecies+1), Xi_Vib(nSpecies)
  REAL                :: MaxCollProb, MeanCollProb
#ifdef MPI
  REAL                :: sumMeanCollProb
#endif /*MPI*/
#ifdef MPI
  REAL                :: RECBR(nSpecies)
  INTEGER             :: RECBIM(nSpecies)
#endif /*MPI*/
  REAL, ALLOCATABLE   :: CRate(:), RRate(:)
#if (PP_TimeDiscMethod ==42)
  INTEGER             :: ii, iunit, iCase, iTvib,jSpec, WallNumSpec(nSpecies)
  CHARACTER(LEN=64)   :: DebugElectronicStateFilename
  CHARACTER(LEN=350)  :: hilf
  REAL                :: WallCoverage(nSpecies), Accomodation(nSpecies), Adsorptionrate(nSpecies), Desorptionrate(nSpecies)
  INTEGER             :: iCov
  REAL                :: NumSpecTmp(nSpecies+1)
#endif
  REAL                :: PartVtrans(nSpecies,4) ! macroscopic velocity (drift velocity) A. Frohn: kinetische Gastheorie
  REAL                :: PartVtherm(nSpecies,4) ! microscopic velocity (eigen velocity) PartVtrans + PartVtherm = PartVtotal
  INTEGER             :: dir
#ifdef MPI 
! load balance
  REAL                :: tLBStart,tLBEnd
#endif /*MPI*/
!===================================================================================================================================
  IF ( DoRestart ) THEN
    isRestart = .true.
  END IF
  IF (.NOT.DoAnalyze) RETURN
#ifdef MPI 
  tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
  IF (useDSMC) THEN
    IF (CollisMode.NE.0) THEN
      SDEALLOCATE(CRate)
      ALLOCATE(CRate(CollInf%NumCase + 1))
      IF (CollisMode.EQ.3) THEN
        SDEALLOCATE(RRate)
        ALLOCATE(RRate(ChemReac%NumOfReact))
        RRate = 0.0
      END IF
    END IF
  END IF
  OutputCounter = 2
  unit_index = 535
#ifdef MPI
!#ifdef PARTICLES
  IF (PartMPI%MPIRoot) THEN
#endif    /* MPI */
    INQUIRE(UNIT   = unit_index , OPENED = isOpen)
    IF (.NOT.isOpen) THEN
#if (PP_TimeDiscMethod==42)
    ! if only the reaction rate is desired (resevoir) the initial temperature
    ! of the second species is added to the filename
      IF (DSMC%ReservoirSimuRate) THEN
        IF ( SpecDSMC(1)%InterID .EQ. 2 .OR. SpecDSMC(1)%InterID .EQ. 20 ) THEN
          iTvib = INT(SpecDSMC(1)%Init(0)%Tvib)
          WRITE( hilf, '(I5.5)') iTvib
          outfile = 'Database_Tvib_'//TRIM(hilf)//'.csv'
        ELSE
          !iTvib = INT(SpecDSMC(1)%Telec )
          iTvib = INT(Species(1)%Init(0)%MWTemperatureIC) !wrong name, if MWTemp is defined in %Init!!!
          WRITE( hilf, '(I5.5)') iTvib
          outfile = 'Database_Ttrans_'//TRIM(hilf)//'.csv'
        END IF
      ELSE IF (Adsorption%TPD) THEN
          iCov = INT(Adsorption%Coverage(1,1,1,1)*1000)
          WRITE( hilf, '(I4.4)') iCov
        outfile = 'Database_Cov_'//TRIM(hilf)//'.csv'
      ELSE
        outfile = 'Database.csv'
      END IF
#else
      outfile = 'Database.csv'
#endif

      INQUIRE(file=TRIM(outfile),EXIST=FileExists)
      IF (isRestart .and. FileExists) THEN
        OPEN(unit_index,file=TRIM(outfile),position="APPEND",status="OLD")
        !CALL FLUSH (unit_index)
      ELSE
        OPEN(unit_index,file=TRIM(outfile))
        !CALL FLUSH (unit_index)
        !--- insert header
        WRITE(unit_index,'(A6,A5)',ADVANCE='NO') 'TIME', ' '
        IF (CalcNumSpec) THEN
          DO iSpec = 1, nSpecies + 1
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A12,I3.3,A5)',ADVANCE='NO') OutputCounter,'-nPart-Spec-', iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
        END IF
        IF (CalcCharge) THEN
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A12,A5)',ADVANCE='NO') OutputCounter,'-Charge     '
          OutputCounter = OutputCounter + 1
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A16,A5)',ADVANCE='NO') OutputCounter,'-Charge-absError'
          OutputCounter = OutputCounter + 1
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A16,A5)',ADVANCE='NO') OutputCounter,'-Charge-relError'
          OutputCounter = OutputCounter + 1
        END IF
        IF (CalcPartBalance) THEN
          DO iSpec=1, nSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A14,I3.3,A5)',ADVANCE='NO') OutputCounter,'-nPartIn-Spec-',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
          DO iSpec=1, nSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A15,I3.3,A5)',ADVANCE='NO') OutputCounter,'-nPartOut-Spec-',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
        END IF
        IF (CalcEpot) THEN 
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A11)',ADVANCE='NO') OutputCounter,'-W-El      '
          OutputCounter = OutputCounter + 1
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A11)',ADVANCE='NO') OutputCounter,'-W-Mag    '
          OutputCounter = OutputCounter + 1
        END IF
        IF (CalcEkin) THEN
          DO iSpec=1, nEkin
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-Ekin-',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
        END IF
        IF (CalcTemp) THEN
          DO iSpec=1, nSpecies + 1
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' -TempTra-',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
        END IF
        IF (CalcVelos) THEN
          DO iSpec=1, nSpecies
            IF (VeloDirs(1)) THEN
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' Velo_Xtrans',iSpec,' '
              OutputCounter = OutputCounter + 1
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' Velo_Xtherm',iSpec,' '
              OutputCounter = OutputCounter + 1
            END IF
            IF (VeloDirs(2)) THEN
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' Velo_Ytrans',iSpec,' '
              OutputCounter = OutputCounter + 1
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' Velo_Ytherm',iSpec,' '
              OutputCounter = OutputCounter + 1
            END IF
            IF (VeloDirs(3)) THEN
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' Velo_Ztrans',iSpec,' '
              OutputCounter = OutputCounter + 1
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' Velo_Ztherm',iSpec,' '
              OutputCounter = OutputCounter + 1
            END IF
            IF (VeloDirs(4)) THEN
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' AbsVelo_trans',iSpec,' '
              OutputCounter = OutputCounter + 1
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' AbsVelo_therm',iSpec,' '
              OutputCounter = OutputCounter + 1
            END IF
          END DO
        END IF
        IF (CalcPartBalance) THEN
          DO iSpec=1, nSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A8,I3.3,A5)',ADVANCE='NO') OutputCounter,'-EkinIn-',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
          DO iSpec=1, nSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A9,I3.3,A5)',ADVANCE='NO') OutputCounter,'-EkinOut-',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
        END IF
#if (PP_TimeDiscMethod==1000)
        IF (CollisMode.GT.1) THEN
          DO iSpec=1, nSpecies         
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' EVib',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
          DO iSpec=1, nSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' ERot',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
          IF ( DSMC%ElectronicModel ) THEN
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' EElec',iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
          END IF
          DO iSpec=1, nSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' TempVib',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
          DO iSpec=1, nSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' TempRot',iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
        END IF
#endif
        IF (CollisMode.GT.1) THEN
          IF(CalcEint) THEN
            DO iSpec=1, nSpecies         
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-EVib',iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
            DO iSpec=1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-ERot',iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
            IF (DSMC%ElectronicModel) THEN
              DO iSpec = 1, nSpecies
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-EElec',iSpec,' '
                OutputCounter = OutputCounter + 1
              END DO
            END IF
          END IF
          IF(CalcTemp) THEN
            DO iSpec=1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-TempVib',iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
            DO iSpec=1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-XiVibMean',iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
            DO iSpec=1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-TempRot',iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
            DO iSpec=1, nSpecies+1
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-TempTotal',iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
            IF ( DSMC%ElectronicModel ) THEN
              DO iSpec=1, nSpecies
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-TempElec',iSpec,' '
                OutputCounter = OutputCounter + 1
              END DO
            END IF
          END IF
        END IF
        IF(DSMC%CalcQualityFactors) THEN
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,'-Pmean',' '
          OutputCounter = OutputCounter + 1
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,'-Pmax',' '
          OutputCounter = OutputCounter + 1
        END IF
#if (PP_TimeDiscMethod==42)
        IF (DSMC%WallModel.GE.1) THEN
!                 IF (CalcWallNumSpec) THEN
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-nPart-Wall-Spec-', iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
!                 END IF
!                 IF (CalcWallCoverage) THEN
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,'-Surf-Cov-', iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
!                 END IF
          DO iSpec = 1, nSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' alpha', iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
          DO iSpec = 1, nSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' Pads', iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
          DO iSpec = 1, nSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' Pdes', iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
          DO iSpec = 1, nSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' Nads', iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
          DO iSpec = 1, nSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' Ndes', iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
          DO iSpec = 1, nSpecies
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' E-ads', iSpec,' '
            OutputCounter = OutputCounter + 1
          END DO
          IF (Adsorption%TPD) THEN
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' WallTemp', iSpec,' '
          END IF
          OutputCounter = OutputCounter + 1
        END IF
        IF(CalcCollRates) THEN
          DO iSpec = 1, nSpecies
            DO jSpec = iSpec, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' CollRate', iSpec, '+', jSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
          END DO
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,' TotalCollRate',' '
          OutputCounter = OutputCounter + 1
        END IF
        IF(CalcReacRates) THEN
          IF(CollisMode.EQ.3) THEN
            DO iCase=1, ChemReac%NumOfReact
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' Reaction', iCase,' '
              OutputCounter = OutputCounter + 1
            END DO
          END IF
        END IF
#endif
        WRITE(unit_index,'(A1)') ' '
      END IF
    END IF
#ifdef MPI
  END IF
#endif    /* MPI */

!===================================================================================================================================
! Analyze Routines
!===================================================================================================================================

!-----------------------------------------------------------------------------------------------------------------------------------
! Determine the maximal collision probability for whole reservoir and mean collision probability (only for one cell)
  IF((iter.GT.0).AND.(DSMC%CalcQualityFactors).AND.(DSMC%CollProbMeanCount.GT.0)) THEN
    MaxCollProb = DSMC%CollProbMax
    MeanCollProb = DSMC%CollProbMean / DSMC%CollProbMeanCount
  ELSE
    MaxCollProb = 0.0
    MeanCollProb = 0.0
  END IF

  ! computes the real and simulated number of particles
  CALL CalcNumPartsofSpec(NumSpec,SimNumSpec)
!-----------------------------------------------------------------------------------------------------------------------------------
! Calculate total temperature of each molecular species (Laux, p. 109)
  IF(CalcTemp) THEN
    CALL CalcTemperature(NumSpec,Temp,IntTemp,IntEn,TempTotal,Xi_Vib) ! contains MPI Communication
  ELSEIF(CalcEint.AND.(CollisMode.GT.1)) THEN
    CALL CalcIntTempsAndEn(IntTemp, IntEn,NumSpec)
  END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! Other Analyze Routines
  IF(CalcCharge) CALL CalcDepositedCharge() ! mpi communication done in calcdepositedcharge
  IF(CalcEpot) CALL CalcPotentialEnergy(WEl,WMag)
  IF(CalcEkin) CALL CalcKineticEnergy(Ekin)
  IF(TrackParticlePosition) CALL TrackingParticlePosition(time)
  PartVtrans=0.
  PartVtherm=0.
  IF(CalcVelos) CALL CalcVelocities(PartVtrans, PartVtherm)
!===================================================================================================================================
! MPI Communication for values which are not YET communicated
! all routines ABOVE contains the required MPI-Communication
!===================================================================================================================================

! MPI Communication
#ifdef MPI
  IF (PartMPI%MPIRoot) THEN
    tLBStart = LOCALTIME() ! LB Time Start
    IF (CalcPartBalance)THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,nPartIn(:)    ,nSpecies,MPI_INTEGER         ,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(MPI_IN_PLACE,nPartOUt(:)   ,nSpecies,MPI_INTEGER         ,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(MPI_IN_PLACE,PartEkinIn(:) ,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(MPI_IN_PLACE,PartEkinOut(:),nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
    END IF
      ! Determining the maximal (MPI_MAX) and mean (MPI_SUM) collision probabilities
      CALL MPI_REDUCE(MPI_IN_PLACE,MaxCollProb,1, MPI_DOUBLE_PRECISION, MPI_MAX,0, PartMPI%COMM, IERROR)
      CALL MPI_REDUCE(MeanCollProb,sumMeanCollProb,1, MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
      MeanCollProb = sumMeanCollProb / PartMPI%nProcs
  ELSE ! no Root
    tLBStart = LOCALTIME() ! LB Time Start
    IF (CalcPartBalance)THEN
      CALL MPI_REDUCE(nPartIn,RECBIM    ,nSpecies,MPI_INTEGER         ,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(nPartOut,RECBIM   ,nSpecies,MPI_INTEGER         ,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(PartEkinIn,RECBR ,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
      CALL MPI_REDUCE(PartEkinOut,RECBR,nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
    END IF
      CALL MPI_REDUCE(MaxCollProb,RECBR,1,MPI_DOUBLE_PRECISION,MPI_MAX,0, PartMPI%COMM, IERROR)
      CALL MPI_REDUCE(MeanCollProb,sumMeanCollProb,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, PartMPI%COMM, IERROR)
    tLBEnd = LOCALTIME() ! LB Time End
    tCurrent(14)=tCurrent(14)+tLBEnd-tLBStart
  END IF
#endif

#ifdef MPI 
  tLBEnd = LOCALTIME() ! LB Time End
  tCurrent(14)=tCurrent(14)+tLBEnd-tLBStart
#endif /*MPI*/
#ifdef MPI 
tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
!-----------------------------------------------------------------------------------------------------------------------------------
#if (PP_TimeDiscMethod==1000)
  IF (CollisMode.GT.1) CALL CalcIntTempsAndEn(IntTemp, IntEn,NumSpec) 
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! Calculate the collision rates and reaction rate coefficients (Arrhenius-type chemistry)
#if (PP_TimeDiscMethod==42)
  IF(CalcCollRates) CALL CollRates(CRate)
  IF(CalcReacRates) THEN
    IF ((CollisMode.EQ.3).AND.(iter.GT.0)) THEN
      NumSpecTmp = NumSpec
      IF(BGGas%BGGasSpecies.NE.0) THEN
        NumSpecTmp(BGGas%BGGasSpecies) = (BGGas%BGGasDensity * GEO%MeshVolume / Species(BGGas%BGGasSpecies)%MacroParticleFactor)
        NumSpecTmp(nSpecies+1) = NumSpecTmp(nSpecies+1)+NumSpecTmp(BGGas%BGGasSpecies)
      END IF
      CALL ReacRates(RRate, NumSpecTmp)
    END IF
  END IF
IF (DSMC%WallModel.GE.1) THEN
  CALL GetWallNumSpec(WallNumSpec,WallCoverage)
  CALL CalcSurfRates(WallNumSpec,Accomodation,Adsorptionrate,Desorptionrate)
END IF
#endif /*PP_TimeDiscMethod==42*/
!-----------------------------------------------------------------------------------------------------------------------------------
  IF (CalcShapeEfficiency) CALL CalcShapeEfficiencyR()   ! This will NOT be placed in the file but directly in "out"
!===================================================================================================================================
! Output Routines
!===================================================================================================================================
#ifdef MPI
IF (PartMPI%MPIROOT) THEN
#endif    /* MPI */
  WRITE(unit_index,104,ADVANCE='NO') Time
    IF (CalcNumSpec) THEN
      DO iSpec=1, nSpecies + 1
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,104,ADVANCE='NO') REAL(SimNumSpec(iSpec))
      END DO
    END IF
    IF (CalcCharge) THEN 
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,104,ADVANCE='NO') PartCharge(1)
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,104,ADVANCE='NO') PartCharge(2)       
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,104,ADVANCE='NO') PartCharge(3)
    END IF
    IF (CalcPartBalance) THEN
      DO iSpec=1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,104,ADVANCE='NO') REAL(nPartIn(iSpec))
      END DO
      DO iSpec=1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,104,ADVANCE='NO') REAL(nPartOut(iSpec))
      END DO
    END IF
    IF (CalcEpot) THEN 
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,104,ADVANCE='NO') WEl
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,104,ADVANCE='NO') WMag
    END IF
    IF (CalcEkin) THEN 
      DO iSpec=1, nEkin
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,104,ADVANCE='NO') Ekin(iSpec)
      END DO
    END IF
    IF (CalcTemp) THEN
      DO iSpec=1, nSpecies + 1
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,104,ADVANCE='NO') Temp(iSpec)
      END DO
    END IF
    IF (CalcVelos) THEN
      DO iSpec=1, nSpecies
        DO dir = 1,4
          IF (VeloDirs(dir)) THEN
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,104,ADVANCE='NO') PartVtrans(iSpec,dir)
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,104,ADVANCE='NO') PartVtherm(iSpec,dir)
          END IF
        END DO
      END DO
    END IF
    IF (CalcPartBalance) THEN
      DO iSpec=1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,104,ADVANCE='NO') PartEkinIn(iSpec)
      END DO
      DO iSpec=1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,104,ADVANCE='NO') PartEkinOut(iSpec)
      END DO
    END IF

#if (PP_TimeDiscMethod==1000)
    IF (CollisMode.GT.1) THEN
      DO iSpec=1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,104,ADVANCE='NO') IntEn(iSpec,1)
      END DO
      DO iSpec=1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,104,ADVANCE='NO') IntEn(iSpec,2)
      END DO
      IF ( DSMC%ElectronicModel ) THEN
        DO iSpec=1, nSpecies
        ! currently set to one
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,104,ADVANCE='NO') IntEn(iSpec,3)
        END DO
      END IF
      DO iSpec=1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,104,ADVANCE='NO') IntTemp(iSpec,1)
      END DO
      DO iSpec=1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,104,ADVANCE='NO') IntTemp(iSpec,2)
      END DO
    END IF
#endif

    IF (CollisMode.GT.1) THEN
      IF(CalcEint) THEN
        DO iSpec=1, nSpecies
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,104,ADVANCE='NO') IntEn(iSpec,1)
        END DO
        DO iSpec=1, nSpecies
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,104,ADVANCE='NO') IntEn(iSpec,2)
        END DO
        IF (DSMC%ElectronicModel) THEN
          DO iSpec=1, nSpecies
          ! currently set to one
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,104,ADVANCE='NO') IntEn(iSpec,3)
          END DO
        END IF
      END IF
      IF(CalcTemp) THEN
        DO iSpec=1, nSpecies
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,104,ADVANCE='NO') IntTemp(iSpec,1)
        END DO
        DO iSpec=1, nSpecies
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,104,ADVANCE='NO') Xi_Vib(iSpec)
        END DO
        DO iSpec=1, nSpecies
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,104,ADVANCE='NO') IntTemp(iSpec,2)
        END DO
        DO iSpec=1, nSpecies+1
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,104,ADVANCE='NO') TempTotal(iSpec)
        END DO
        IF ( DSMC%ElectronicModel ) THEN
          DO iSpec=1, nSpecies
          ! currently set to one
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,104,ADVANCE='NO') IntTemp(iSpec,3)
          END DO
        END IF
      END IF
    END IF
    IF(DSMC%CalcQualityFactors) THEN
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,104,ADVANCE='NO') MeanCollProb
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,104,ADVANCE='NO') MaxCollProb
    END IF
#if (PP_TimeDiscMethod==42)
! output for adsorption
    IF (DSMC%WallModel.GE.1) THEN
!         IF (CalcWallNumSpec) THEN
      DO iSpec=1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,'(I18.1)',ADVANCE='NO') WallNumSpec(iSpec)
      END DO
!         END IF
!         IF (CalcWallCoverage) THEN
      DO iSpec=1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,104,ADVANCE='NO') WallCoverage(iSpec)
      END DO
      DO iSpec = 1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,104,ADVANCE='NO') Accomodation(iSpec)
      END DO
      DO iSpec = 1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,104,ADVANCE='NO') Adsorptionrate(iSpec)
      END DO
      DO iSpec = 1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,104,ADVANCE='NO') Desorptionrate(iSpec)
      END DO
      DO iSpec = 1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,'(I18.1)',ADVANCE='NO') Adsorption%AdsorpInfo(iSpec)%NumOfAds
      END DO
      DO iSpec = 1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,'(I18.1)',ADVANCE='NO') Adsorption%AdsorpInfo(iSpec)%NumOfDes
      END DO
      DO iSpec = 1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,104,ADVANCE='NO') Adsorption%AdsorpInfo(iSpec)%MeanEads
      END DO
      IF (Adsorption%TPD) THEN
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,104,ADVANCE='NO') Adsorption%TPD_Temp
      END IF
    END IF
    IF(CalcCollRates) THEN
      DO iCase=1, CollInf%NumCase +1 
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,104,ADVANCE='NO') CRate(iCase)
      END DO
    END IF
    IF(CalcReacRates) THEN
      DO iCase=1, ChemReac%NumOfReact 
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,104,ADVANCE='NO') RRate(iCase)
      END DO
    END IF
#endif /*(PP_TimeDiscMethod==42)*/
    WRITE(unit_index,'(A1)') ' ' 
#ifdef MPI
  END IF
#endif    /* MPI */

  104    FORMAT (e25.14)

!-----------------------------------------------------------------------------------------------------------------------------------
  IF( CalcPartBalance) CALL CalcParticleBalance()
!-----------------------------------------------------------------------------------------------------------------------------------
#ifdef MPI 
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(14)=tCurrent(14)+tLBEnd-tLBStart
#endif /*MPI*/

#if ( PP_TimeDiscMethod ==42 )
! hard coded
! array not allocated
  IF ( DSMC%ElectronicModel ) THEN
  IF (DSMC%ReservoirSimuRate) THEN
    IF(Time.GT.0.) CALL ElectronicTransition( Time, NumSpec )
  END IF
END IF
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
#if ( PP_TimeDiscMethod ==42 )
  ! Debug Output for initialized electronic state
  IF ( DSMC%ElectronicModel .AND. DSMC%ReservoirSimuRate) THEN
    DO iSpec = 1, nSpecies
      IF ( SpecDSMC(iSpec)%InterID .ne. 4 ) THEN
        IF (  SpecDSMC(iSpec)%levelcounter(0) .ne. 0) THEN
          WRITE(DebugElectronicStateFilename,'(I2.2)') iSpec
          iunit = 485
          DebugElectronicStateFilename = 'End_Electronic_State_Species_'//trim(DebugElectronicStateFilename)//'.dat'
          OPEN(unit=iunit,file=DebugElectronicStateFilename,form='formatted',status='unknown')
          DO ii = 0, SpecDSMC(iSpec)%MaxElecQuant - 1                         !has to be changed when using %Init definitions!!!
            WRITE(iunit,'(I3.1,3x,F12.7)') ii, REAL( SpecDSMC(iSpec)%levelcounter(ii) ) / &
                                               REAL( Species(iSpec)%Init(0)%initialParticleNumber)
          END DO
          CLOSE(iunit)
        END IF
      END IF
    END DO
  END IF
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE AnalyzeParticles

! all other analysis with particles
SUBROUTINE CalcShapeEfficiencyR()
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Particle_Analyze_Vars,   ONLY : CalcShapeEfficiencyMethod, ShapeEfficiencyNumber
USE MOD_Mesh_Vars,               ONLY : nElems, Elem_xGP
USE MOD_Particle_Mesh_Vars,      ONLY : GEO
USE MOD_PICDepo_Vars
USE MOD_Particle_Vars
USE MOD_PreProc
USE MOD_Particle_MPI_Vars,       ONLY : PartMPI
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                     :: NbrOfComps, NbrWithinRadius, NbrOfElems, NbrOfElemsWithinRadius
REAL                     :: RandVal1
LOGICAL                  :: chargedone(1:nElems), WITHIN
INTEGER                  :: kmin, kmax, lmin, lmax, mmin, mmax                           !
INTEGER                  :: kk, ll, mm, ppp,m,l,k, i                                             !
INTEGER                  :: ElemID
REAL                     :: radius, deltax, deltay, deltaz
!===================================================================================================================================

NbrOfComps = 0.
NbrOfElems = 0.
NbrWithinRadius = 0.
NbrOfElemsWithinRadius = 0.
SELECT CASE(CalcShapeEfficiencyMethod)
CASE('AllParts')
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      chargedone(:) = .FALSE.
      !-- determine which background mesh cells (and interpolation points within) need to be considered
      kmax = INT((PartState(i,1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
      kmax = MIN(kmax,GEO%FIBGMimax)
      kmin = INT((PartState(i,1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
      kmin = MAX(kmin,GEO%FIBGMimin)
      lmax = INT((PartState(i,2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
      lmax = MIN(lmax,GEO%FIBGMjmax)
      lmin = INT((PartState(i,2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
      lmin = MAX(lmin,GEO%FIBGMjmin)
      mmax = INT((PartState(i,3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
      mmax = MIN(mmax,GEO%FIBGMkmax)
      mmin = INT((PartState(i,3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
      mmin = MAX(mmin,GEO%FIBGMkmin)
      !-- go through all these cells
      DO kk = kmin,kmax
        DO ll = lmin, lmax
          DO mm = mmin, mmax
            !--- go through all mapped elements not done yet
            DO ppp = 1,GEO%FIBGM(kk,ll,mm)%nElem
              WITHIN=.FALSE.
              ElemID = GEO%FIBGM(kk,ll,mm)%Element(ppp)
              IF (.NOT.chargedone(ElemID)) THEN
                NbrOfElems = NbrOfElems + 1.
                !--- go through all gauss points
                DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
                  NbrOfComps = NbrOfComps + 1.
                  !-- calculate distance between gauss and particle
                  deltax = PartState(i,1) - Elem_xGP(1,k,l,m,ElemID)
                  deltay = PartState(i,2) - Elem_xGP(2,k,l,m,ElemID)
                  deltaz = PartState(i,3) - Elem_xGP(3,k,l,m,ElemID)
                  radius = deltax * deltax + deltay * deltay + deltaz * deltaz
                  IF (radius .LT. r2_sf) THEN
                    WITHIN=.TRUE.
                    NbrWithinRadius = NbrWithinRadius + 1.
                  END IF
                END DO; END DO; END DO
                chargedone(ElemID) = .TRUE.
              END IF
              IF(WITHIN) NbrOfElemsWithinRadius = NbrOfElemsWithinRadius + 1.
            END DO ! ppp
          END DO ! mm
        END DO ! ll
      END DO ! kk
    END IF ! inside
  END DO ! i
IF(NbrOfComps.GT.0.0)THEN
#ifdef MPI
  WRITE(*,*) 'ShapeEfficiency (Proc,%,%Elems)',PartMPI%MyRank,100*NbrWithinRadius/NbrOfComps,100*NbrOfElemsWithinRadius/NbrOfElems
  WRITE(*,*) 'ShapeEfficiency (Elems) for Proc',PartMPI%MyRank,'is',100*NbrOfElemsWithinRadius/NbrOfElems,'%'
#else
  WRITE(*,*) 'ShapeEfficiency (%,%Elems)',100*NbrWithinRadius/NbrOfComps, 100*NbrOfElemsWithinRadius/NbrOfElems
  WRITE(*,*) 'ShapeEfficiency (Elems) is',100*NbrOfElemsWithinRadius/NbrOfElems,'%'
#endif
END IF
CASE('SomeParts')
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      CALL RANDOM_NUMBER(RandVal1)
      IF(RandVal1.LT.REAL(ShapeEfficiencyNumber)/100)THEN
        chargedone(:) = .FALSE.
        !-- determine which background mesh cells (and interpolation points within) need to be considered
        kmax = INT((PartState(i,1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
        kmax = MIN(kmax,GEO%FIBGMimax)
        kmin = INT((PartState(i,1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
        kmin = MAX(kmin,GEO%FIBGMimin)
        lmax = INT((PartState(i,2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
        lmax = MIN(lmax,GEO%FIBGMjmax)
        lmin = INT((PartState(i,2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
        lmin = MAX(lmin,GEO%FIBGMjmin)
        mmax = INT((PartState(i,3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
        mmax = MIN(mmax,GEO%FIBGMkmax)
        mmin = INT((PartState(i,3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
        mmin = MAX(mmin,GEO%FIBGMkmin)
        !-- go through all these cells
        DO kk = kmin,kmax
          DO ll = lmin, lmax
            DO mm = mmin, mmax
              !--- go through all mapped elements not done yet
              DO ppp = 1,GEO%FIBGM(kk,ll,mm)%nElem
              WITHIN=.FALSE.
                ElemID = GEO%FIBGM(kk,ll,mm)%Element(ppp)
                IF (.NOT.chargedone(ElemID)) THEN
                  NbrOfElems = NbrOfElems + 1
                  !--- go through all gauss points
                  DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
                    NbrOfComps = NbrOfComps + 1
                    !-- calculate distance between gauss and particle
                    deltax = PartState(i,1) - Elem_xGP(1,k,l,m,ElemID)
                    deltay = PartState(i,2) - Elem_xGP(2,k,l,m,ElemID)
                    deltaz = PartState(i,3) - Elem_xGP(3,k,l,m,ElemID)
                    radius = deltax * deltax + deltay * deltay + deltaz * deltaz
                    IF (radius .LT. r2_sf) THEN
                      NbrWithinRadius = NbrWithinRadius + 1
                      WITHIN=.TRUE.
                    END IF
                    END DO; END DO; END DO
                    chargedone(ElemID) = .TRUE.
                  END IF
                  IF(WITHIN) NbrOfElemsWithinRadius = NbrOfElemsWithinRadius + 1
                END DO ! ppp
              END DO ! mm
            END DO ! ll
          END DO ! kk
        END IF  ! RandVal
      END IF ! inside
  END DO ! i
IF(NbrOfComps.GT.0)THEN
#ifdef MPI
  WRITE(*,*) 'ShapeEfficiency (Proc,%,%Elems)',PartMPI%MyRank,100*NbrWithinRadius/NbrOfComps,100*NbrOfElemsWithinRadius/NbrOfElems
  WRITE(*,*) 'ShapeEfficiency (Elems) for Proc',PartMPI%MyRank,'is',100*NbrOfElemsWithinRadius/NbrOfElems,'%'
#else
  WRITE(*,*) 'ShapeEfficiency (%,%Elems)',100*NbrWithinRadius/NbrOfComps,100*NbrOfElemsWithinRadius/NbrOfElems
  WRITE(*,*) 'ShapeEfficiency (Elems) is',100*NbrOfElemsWithinRadius/NbrOfElems,'%'
#endif
END IF
END SELECT
END SUBROUTINE CalcShapeEfficiencyR


#if (PP_TimeDiscMethod==42)
SUBROUTINE GetWallNumSpec(WallNumSpec,WallCoverage)
!===================================================================================================================================
! Calculate number of wallparticles for all species
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars,      ONLY : Species, PartSpecies, PDM, nSpecies, KeepWallParticles
USE MOD_DSMC_Vars,          ONLY : Adsorption, DSMC
USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(OUT)            :: WallNumSpec(nSpecies)
REAL   , INTENT(OUT)            :: WallCoverage(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i, iSurfSide, p, q
REAL              :: Surface, Coverage(nSpecies), SurfPartDens
!===================================================================================================================================
  WallNumSpec = 0
  SurfPartDens = 0
  Surface = 0.
  Coverage(:) = 0.
  DO i=1,nSpecies
  DO iSurfSide=1,SurfMesh%nSides
    DO q = 1,nSurfSample
      DO p = 1,nSurfSample
        Surface = Surface + SurfMesh%SurfaceArea(p,q,iSurfSide)
        Coverage(i) = Coverage(i) + Adsorption%Coverage(p,q,iSurfSide,i)
      END DO
    END DO
    SurfPartDens = SurfPartDens + Adsorption%DensSurfAtoms(iSurfSide)/Species(i)%MacroParticleFactor
  END DO
  END DO
  WallCoverage(:) = Coverage(:) / (SurfMesh%nSides*nSurfSample*nSurfSample)
  SurfPartDens = SurfPartDens / SurfMesh%nSides
  
  IF (KeepWallParticles) THEN
    DO i=1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(i) .AND. PDM%ParticleAtWall(i)) THEN
        WallNumSpec(PartSpecies(i)) = WallNumSpec(PartSpecies(i)) + 1
      END IF
    END DO
  ELSE 
    DO i=1,nSpecies
      IF ( (Adsorption%Coordination(i).EQ.2) .AND. (DSMC%WallModel.EQ.3) ) THEN
        WallNumSpec(i) = INT(WallCoverage(i) * Surface * SurfPartDens * 2.)
      ELSE
        WallNumSpec(i) = INT(WallCoverage(i) * Surface * SurfPartDens)
      END IF
      IF (WallNumSpec(i).EQ.0) THEN
       WallCoverage(i) = 0.
      END IF
    END DO
  END IF
    
END SUBROUTINE GetWallNumSpec

SUBROUTINE CalcSurfRates(WallNumSpec,Accomodation,Adsorbrate,Desorbrate)
!===================================================================================================================================
! Calculate number of wallparticles for all species
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars,      ONLY : nSpecies
USE MOD_DSMC_Vars,          ONLY : Adsorption, DSMC  
USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(IN)             :: WallNumSpec(nSpecies)
REAL   , INTENT(OUT)            :: Adsorbrate(nSpecies), Desorbrate(nSpecies), Accomodation(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iSpec
!===================================================================================================================================

IF (DSMC%ReservoirRateStatistic) THEN
  DO iSpec = 1,nSpecies
    IF (Adsorption%AdsorpInfo(iSpec)%WallCollCount.GT.0) THEN
        Adsorbrate(iSpec) = Adsorption%AdsorpInfo(iSpec)%NumOfAds / Adsorption%AdsorpInfo(iSpec)%WallCollCount
        Accomodation(iSpec) = Adsorption%AdsorpInfo(iSpec)%Accomodation / Adsorption%AdsorpInfo(iSpec)%WallCollCount
        Desorbrate(iSpec) = Adsorption%AdsorpInfo(iSpec)%NumOfDes / WallNumSpec(iSpec)
      ELSE
        Adsorbrate(iSpec) = 0.
        Accomodation(iSpec) = 0.
        Desorbrate(iSpec) = 0.
      END IF
  END DO
ELSE

  IF ((.NOT.DSMC%ReservoirRateStatistic).AND.(DSMC%WallModel.EQ.3)) THEN
    DO iSpec = 1,nSpecies
      IF (Adsorption%AdsorpInfo(iSpec)%WallCollCount.GT.0) THEN
        Adsorption%AdsorpInfo(iSpec)%MeanProbAds = Adsorption%AdsorpInfo(iSpec)%MeanProbAds / (REAL(nSurfSample) * REAL(nSurfSample) &
                                            * REAL(SurfMesh%nSides) * REAL(Adsorption%AdsorpInfo(iSpec)%WallCollCount))
        Accomodation(iSpec) = Adsorption%AdsorpInfo(iSpec)%Accomodation / Adsorption%AdsorpInfo(iSpec)%WallCollCount
      ELSE
        Adsorption%AdsorpInfo(iSpec)%MeanProbAds = 0.
        Accomodation(iSpec) = 0.
      END IF
    END DO
  END IF
  
  DO iSpec = 1,nSpecies
    Adsorbrate(iSpec) = Adsorption%AdsorpInfo(iSpec)%MeanProbAds
    Desorbrate(iSpec) = Adsorption%AdsorpInfo(iSpec)%MeanProbDes
  END DO
END IF

END SUBROUTINE CalcSurfRates
#endif

SUBROUTINE CalcParticleBalance()
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars,      ONLY : nPartIn,nPartOut,PartEkinIn,PartEkinOut
#if defined(LSERK) || defined(IMEX) || defined(IMPA)
!#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506)
USE MOD_Particle_Analyze_Vars,      ONLY : nPartInTmp,PartEkinInTmp
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

#if defined(LSERK) || defined(IMEX) || defined(IMPA)
!#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506)
nPartIn=nPartInTmp
nPartOut=0
PartEkinIn=PartEkinInTmp
PartEkinOut=0.
nPartInTmp=0
PartEkinInTmp=0.
#else
nPartIn=0
nPartOut=0
PartEkinIn=0.
PartEkinOut=0.
#endif

END SUBROUTINE CalcParticleBalance


SUBROUTINE CalcKineticEnergy(Ekin) 
!===================================================================================================================================
! compute the kinetic energy of particles
! for velocity <1e3 non-relativistic formula is used, for larger velocities the relativistic kinetic energy is computed
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Equation_Vars,          ONLY : c2, c2_inv
USE MOD_Particle_Vars,          ONLY : PartState, PartSpecies, Species, PDM
USE MOD_PARTICLE_Vars,          ONLY : nSpecies, PartMPF, usevMPF
USE MOD_Particle_Analyze_Vars,  ONLY : nEkin
#ifndef PP_HDG
USE MOD_PML_Vars,               ONLY : DoPML,xyzPhysicalMinMax
#endif /*PP_HDG*/ 
#ifdef MPI
USE MOD_Particle_MPI_Vars,      ONLY : PartMPI
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Ekin(nEkin) 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i
REAL(KIND=8)              :: partV2, Gamma
#ifdef MPI
REAL                      :: RD(nEKIN)
#endif /*MPI*/
!===================================================================================================================================

Ekin = 0.!d0
IF (nEkin .GT. 1 ) THEN
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
#ifndef PP_HDG
      IF(DoPML)THEN
        IF (PartState(i,1) .GE. xyzPhysicalMinMax(1) .AND. PartState(i,1) .LE. xyzPhysicalMinMax(2) .AND. &
            PartState(i,2) .GE. xyzPhysicalMinMax(3) .AND. PartState(i,2) .LE. xyzPhysicalMinMax(4) .AND. &
            PartState(i,3) .GE. xyzPhysicalMinMax(5) .AND. PartState(i,3) .LE. xyzPhysicalMinMax(6)) THEN        
          CYCLE
        END IF
      ENDIF
#endif /*PP_HDG*/ 
      partV2 = PartState(i,4) * PartState(i,4) &
              + PartState(i,5) * PartState(i,5) &
              + PartState(i,6) * PartState(i,6)
      IF ( partV2 .LT. 1e6) THEN  ! |v| < 1000
  !       Ekin = Ekin + 0.5 *  Species(PartSpecies(i))%MassIC * partV2 * PartMPF(i)            
        IF(usevMPF) THEN
          Ekin(nSpecies+1)     = Ekin(nSpecies+1)     + 0.5 * Species(PartSpecies(i))%MassIC * partV2 * PartMPF(i)
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + 0.5 * Species(PartSpecies(i))%MassIC * partV2 * PartMPF(i)
        ELSE
          Ekin(nSpecies+1) = Ekin(nSpecies+1) + 0.5 *  Species(PartSpecies(i))%MassIC * partV2 &
                                            *  Species(PartSpecies(i))%MacroParticleFactor
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + 0.5 *  Species(PartSpecies(i))%MassIC * partV2 &
                                            *  Species(PartSpecies(i))%MacroParticleFactor
        END IF != usevMPF
      ELSE ! partV2 > 1e6
  !       Ekin = Ekin + (gamma - 1) * mass * MPF *c^2
        Gamma = partV2*c2_inv
        Gamma = 1./SQRT(1.-Gamma)
        IF(usevMPF) THEN
          Ekin(nSpecies+1)     = Ekin(nSpecies+1)     + PartMPF(i) * (Gamma-1.) * Species(PartSpecies(i))%MassIC * c2
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + PartMPF(i) * (Gamma-1.) * Species(PartSpecies(i))%MassIC * c2
        ELSE
          Ekin(nSpecies+1)     = Ekin(nSpecies+1)     + (Gamma-1.) * Species(PartSpecies(i))%MassIC &
                                                                   * Species(PartSpecies(i))%MacroParticleFactor * c2
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + (Gamma-1.) * Species(PartSpecies(i))%MassIC &
                                                                   * Species(PartSpecies(i))%MacroParticleFactor * c2
        END IF !=usevMPF
      END IF ! partV2
    END IF ! (PDM%ParticleInside(i))
  END DO ! i=1,PDM%ParticleVecLength
ELSE ! nEkin = 1 : only 1 species
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
#ifndef PP_HDG
      IF(DoPML)THEN
        IF (PartState(i,1) .GE. xyzPhysicalMinMax(1) .AND. PartState(i,1) .LE. xyzPhysicalMinMax(2) .AND. &
            PartState(i,2) .GE. xyzPhysicalMinMax(3) .AND. PartState(i,2) .LE. xyzPhysicalMinMax(4) .AND. &
            PartState(i,3) .GE. xyzPhysicalMinMax(5) .AND. PartState(i,3) .LE. xyzPhysicalMinMax(6)) THEN        
          CYCLE
        END IF
      ENDIF
#endif /*PP_HDG*/ 
      partV2 = PartState(i,4) * PartState(i,4) &
             + PartState(i,5) * PartState(i,5) &
             + PartState(i,6) * PartState(i,6)
      IF ( partV2 .LT. 1e6) THEN  ! |v| < 1000
        IF(usevMPF) THEN
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + 0.5 *  Species(PartSpecies(i))%MassIC * partV2 &
                                            * PartMPF(i)            
        ELSE
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + 0.5 *  Species(PartSpecies(i))%MassIC * partV2 &
                                            *  Species(PartSpecies(i))%MacroParticleFactor
        END IF ! usevMPF
      ELSE ! partV2 > 1e6
        Gamma = partV2*c2_inv
        Gamma = 1./SQRT(1.-Gamma)
        IF(usevMPF)THEN
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + PartMPF(i) * (Gamma-1.) &
                      * Species(PartSpecies(i))%MassIC * c2
        ELSE
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + (Gamma -1.) &
                      * Species(PartSpecies(i))%MassIC &
                      * Species(PartSpecies(i))%MacroParticleFactor * c2
        END IF ! useuvMPF

      END IF ! par2
    END IF ! particle inside
  END DO ! particleveclength
END IF

#ifdef MPI
IF(PartMPI%MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,Ekin,nEkin,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM, IERROR)
ELSE
  CALL MPI_REDUCE(Ekin  ,RD        ,nEkin,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM, IERROR)
END IF
#endif /*MPI*/

END SUBROUTINE CalcKineticEnergy


SUBROUTINE CalcNumPartsOfSpec(NumSpec,SimNumSpec) 
!===================================================================================================================================
! computes the number of simulated particles AND number of real particles within the domain
! CAUTION: SimNumSpec equals NumSpec only for constant MPF, not vMPF
! Last section of the routine contains the MPI-communication
! Please not, we do not take the paricle number of the BGGas into account
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Particle_Vars,          ONLY: PartMPF, usevMPF, PDM,Species,nSpecies,PartSpecies
USE MOD_Particle_Analyze_Vars,  ONLY: CalcNumSpec
USE MOD_DSMC_Vars,              ONLY: BGGas
USE MOD_Particle_Mesh_Vars,     ONLY: Geo
#ifdef MPI
USE MOD_Particle_MPI_Vars,      ONLY: PartMPI
#endif /*MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(OUT)                   :: NumSpec(nSpecies+1)
INTEGER(KIND=8),INTENT(OUT)        :: SimNumSpec(nSpecies+1)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iPart
#ifdef MPI
REAL                               :: RD(nSpecies+1)
INTEGER(KIND=8)                    :: ID(nSpecies+1)
#endif /*MPI*/
!===================================================================================================================================

IF(usevMPF)THEN ! for MPF differentiate between real particle number and simulated particles
  NumSpec    = 0.
  SimNumSpec = 0
  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
        NumSpec(PartSpecies(iPart))    = NumSpec(PartSpecies(iPart)) + PartMPF(iPart)          ! NumSpec = real particle number
        SimNumSpec(PartSpecies(iPart)) = SimNumSpec(PartSpecies(iPart)) + 1                    ! NumSpec =  particle number
    END IF
  END DO
  IF(BGGas%BGGasSpecies.NE.0) THEN
    NumSpec(BGGas%BGGasSpecies) = BGGas%BGGasDensity * GEO%MeshVolume / Species(BGGas%BGGasSpecies)%MacroParticleFactor
    SimNumSpec(BGGas%BGGasSpecies) = INT(NumSpec(BGGas%BGGasSpecies))
  END IF
  NumSpec(nSpecies + 1)    = SUM(NumSpec(1:nSpecies))
  SimNumSpec(nSpecies + 1) = SUM(SimNumSpec(1:nSpecies))
ELSE ! no mpf-simulated particle number = real particle number
  SimNumSpec = 0
  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
      SimNumSpec(PartSpecies(iPart)) = SimNumSpec(PartSpecies(iPart)) + 1                      ! NumSpec =  particle number
    END IF
  END DO
  NumSpec = REAL(SimNumSpec)
  IF(BGGas%BGGasSpecies.NE.0) THEN
  !  NumSpec(BGGas%BGGasSpecies) = (BGGas%BGGasDensity * GEO%MeshVolume / Species(BGGas%BGGasSpecies)%MacroParticleFactor)
  !  SimNumSpec(BGGas%BGGasSpecies) = 1.
    NumSpec(BGGas%BGGasSpecies) = 0.
    SimNumSpec(BGGas%BGGasSpecies) = 0
  END IF
  SimNumSpec(nSpecies + 1) = SUM(SimNumSpec(1:nSpecies))
  NumSpec(nSpecies + 1) = SUM(NumSpec(1:nSpecies))
END IF

#ifdef MPI
IF (PartMPI%MPIRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,NumSpec    ,nSpecies+1,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  IF(CalcNumSpec) & 
  CALL MPI_REDUCE(MPI_IN_PLACE,SimNumSpec ,nSpecies+1,MPI_LONG            ,MPI_SUM,0,PartMPI%COMM,IERROR)
ELSE
  CALL MPI_REDUCE(NumSpec     ,RD         ,nSpecies+1,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  IF(CalcNumSpec) & 
  CALL MPI_REDUCE(SimNumSpec  ,ID         ,nSpecies+1,MPI_LONG            ,MPI_SUM,0,PartMPI%COMM,IERROR)
END IF
#endif /*MPI*/

END SUBROUTINE CalcNumPartsOfSpec


SUBROUTINE CalcTemperature(NumSpec,Temp,IntTemp,IntEn,TempTotal,Xi_Vib)
!===================================================================================================================================
! computes the temperature, subroutine performs all the MPI communication, do to it at ONE place and CORRECT
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_PARTICLE_Vars,         ONLY: nSpecies
USE MOD_Particle_MPI_Vars,     ONLY: PartMPI
USE MOD_DSMC_Vars,             ONLY: SpecDSMC, PolyatomMolDSMC,CollisMode
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
REAL, INTENT(IN)                   :: NumSpec(nSpecies+1)    ! number of real particles (already GLOBAL number)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL, INTENT(OUT)                  :: Temp(nSpecies+1)
REAL, INTENT(OUT)                  :: IntEn(nSpecies,3)
REAL, INTENT(OUT)                  :: IntTemp(nSpecies,3)
REAL, INTENT(OUT)                  :: TempTotal(nSpecies+1)
REAL, INTENT(OUT)                  :: Xi_Vib(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iPolyatMole,iDOF,iSpec
!===================================================================================================================================


! next, calctranstemp
CALL CalcTransTemp(Temp, NumSpec)

IF (CollisMode.GT.1) THEN
  CALL CalcIntTempsAndEn(IntTemp, IntEn,NumSpec)
  IF(PartMPI%MPIRoot)THEN
    TempTotal = 0.0
    Xi_Vib = 0.0
    DO iSpec = 1, nSpecies
      IF(((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)).AND.(NumSpec(iSpec).GT.0)) THEN
        IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
          IF(IntTemp(iSpec,1).GT.0) THEN
            iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
            DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
              Xi_Vib(iSpec) = Xi_Vib(iSpec) + 2*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/IntTemp(iSpec,1) &
                            /(EXP(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/IntTemp(iSpec,1)) - 1)
            END DO
          ELSE
            Xi_Vib(iSpec) = 0.0
          END IF
          TempTotal(iSpec) = (3*Temp(iSpec)+SpecDSMC(iSpec)%Xi_Rot*IntTemp(iSpec,2) &
                              + Xi_Vib(iSpec)*IntTemp(iSpec,1)) &
                              / (3+SpecDSMC(iSpec)%Xi_Rot+Xi_Vib(iSpec))
        ELSE
          IF(IntTemp(iSpec,1).GT.0) THEN
            Xi_Vib(iSpec) = 2*SpecDSMC(iSpec)%CharaTVib/IntTemp(iSpec,1)/(EXP(SpecDSMC(iSpec)%CharaTVib/IntTemp(iSpec,1)) - 1)
          ELSE
            Xi_Vib(iSpec) = 0.0
          END IF
          TempTotal(iSpec) = (3*Temp(iSpec) + 2*IntTemp(iSpec,2) + Xi_Vib(iSpec)*IntTemp(iSpec,1))/(5+Xi_Vib(iSpec))
        END IF
      ELSE
        TempTotal(iSpec) = Temp(iSpec)
      END IF
      TempTotal(nSpecies+1) = TempTotal(nSpecies+1) + TempTotal(iSpec)*NumSpec(iSpec)
    END DO
    TempTotal(nSpecies+1) = TempTotal(nSpecies+1) / NumSpec(nSpecies+1)
  END IF
END IF

END SUBROUTINE CalcTemperature


SUBROUTINE CalcTransTemp(Temp, NumSpec)
!===================================================================================================================================
! calculate the translational temperature of each species
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars,         ONLY:PartState, PartSpecies, Species, PDM, nSpecies, BoltzmannConst, PartMPF, usevMPF
#if (PP_TimeDiscMethod==1000)
USE MOD_LD_Vars,               ONLY:BulkValues
#endif
USE MOD_Particle_MPI_Vars,     ONLY:PartMPI
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: NumSpec(:)    !< global number of REAL particles in domain
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Temp(:)       !< output value is already the GLOBAL temperature
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i, iSpec
REAL              ::  TempDirec(nSpecies,3)
#if (PP_TimeDiscMethod!=1000)
REAL              :: PartVandV2(nSpecies, 6), Mean_PartV2(nSpecies, 3), MeanPartV_2(nSpecies,3)
#endif
#ifdef MPI
REAL              :: RD(nSpecies*6)
#endif
!===================================================================================================================================

! Compute velocity averages
Temp = 0.0
! Sum up velocity
#if (PP_TimeDiscMethod==1000)
DO iSpec=1, nSpecies
  TempDirec(iSpec,1:3) = BulkValues(1)%BulkTemperature
  Temp(iSpec) = BulkValues(1)%BulkTemperature
END DO
#else 
PartVandV2 = 0.
DO i=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(i)) THEN
    IF (usevMPF) THEN 
      PartVandV2(PartSpecies(i),1:3) = PartVandV2(PartSpecies(i),1:3) + PartState(i,4:6) * PartMPF(i)
      PartVandV2(PartSpecies(i),4:6) = PartVandV2(PartSpecies(i),4:6) + PartState(i,4:6)**2 * PartMPF(i)
    ELSE
      PartVandV2(PartSpecies(i),1:3) = PartVandV2(PartSpecies(i),1:3) + PartState(i,4:6)
      PartVandV2(PartSpecies(i),4:6) = PartVandV2(PartSpecies(i),4:6) + PartState(i,4:6)**2
    END IF
  END IF
END DO

#ifdef MPI
IF(PartMPI%MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,PartVandV2,nSpecies*6,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM, IERROR)
ELSE
  CALL MPI_REDUCE(PartVandV2  ,RD        ,nSpecies*6,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM, IERROR)
END IF
#endif /*MPI*/

IF(PartMPI%MPIRoot)THEN
  DO iSpec=1, nSpecies
    IF(NumSpec(iSpec).NE.0) THEN
      ! Compute velocity averages
      MeanPartV_2(iSpec,1:3)  = (PartVandV2(iSpec,1:3) / NumSpec(iSpec))**2       ! < |v| >**2
      Mean_PartV2(iSpec,1:3)  =  PartVandV2(iSpec,4:6) / NumSpec(iSpec)           ! < |v|**2 >
    ELSE
      MeanPartV_2(iSpec,1:3) = 0.
      Mean_PartV2(iSpec,1:3) = 0.
    END IF
    ! Compute temperatures
    TempDirec(iSpec,1:3) = Species(iSpec)%MassIC * (Mean_PartV2(iSpec,1:3) - MeanPartV_2(iSpec,1:3)) &
         / BoltzmannConst ! Temp calculation is limitedt to one species
    Temp(iSpec) = (TempDirec(iSpec,1) + TempDirec(iSpec,2) + TempDirec(iSpec,3))/3
    Temp(nSpecies + 1) = Temp(nSpecies + 1) + Temp(iSpec)*NumSpec(iSpec)
  END DO
  Temp(nSpecies+1) = Temp(nSpecies + 1) / NumSpec(nSpecies+1)
END IF
#endif
END SUBROUTINE CalcTransTemp


SUBROUTINE CalcVelocities(PartVtrans, PartVtherm)
!===================================================================================================================================
! Calculates the drift and eigen velocity of all particles: PartVtotal = PartVtrans + PartVtherm 
! PartVtrans(nSpecies,4) ! macroscopic velocity (drift velocity) A. Frohn: kinetische Gastheorie
! PartVtherm(nSpecies,4) ! microscopic velocity (eigen velocity)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Analyze_Vars, ONLY: VeloDirs
USE MOD_Particle_Vars,         ONLY: PartState, PartSpecies, PDM, nSpecies, PartMPF, usevMPF
USE MOD_Particle_MPI_Vars,     ONLY: PartMPI
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: PartVtrans(:,:), PartVtherm(:,:) 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iSpec
INTEGER                        :: i
INTEGER                        :: NumSpecloc(nSpecies)
INTEGER                        :: NumSpecglob(nSpecies)
INTEGER                        :: dir
REAL                           :: RealNumSpecloc(nSpecies), RealNumSpecglob(nSpecies)
REAL                           :: PartVglob(nSpecies,4), PartVthermglob(nSpecies,4)
!===================================================================================================================================
! Compute velocity averages
  PartVtrans = 0.
  PartVglob = 0.
  PartVthermglob = 0.
  PartVtherm = 0.
  NumSpecloc = 0
  NumSpecglob = 0
  RealNumSpecloc = 0.
  RealNumSpecglob = 0.

  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      IF (usevMPF) THEN
        RealNumSpecloc(PartSpecies(i)) = RealNumSpecloc(PartSpecies(i)) + PartMPF(i)
      END IF
      NumSpecloc(PartSpecies(i)) = NumSpecloc(PartSpecies(i)) + 1
      DO dir = 1,3
        IF (VeloDirs(dir) .OR. VeloDirs(4)) THEN
          IF (usevMPF) THEN 
            PartVtrans(PartSpecies(i),dir) = PartVtrans(PartSpecies(i),dir) + PartState(i,dir+3) * PartMPF(i)
          ELSE
            PartVtrans(PartSpecies(i),dir) = PartVtrans(PartSpecies(i),dir) + PartState(i,dir+3)
          END IF
        END IF
      END DO
    END IF
  END DO
#ifdef MPI
  CALL MPI_ALLREDUCE(PartVtrans, PartVglob, 4*nSpecies, MPI_DOUBLE_PRECISION, MPI_SUM,  PartMPI%COMM, IERROR)
  PartVtrans = PartVglob
#endif
  IF (usevMPF) THEN
#ifdef MPI
    CALL MPI_ALLREDUCE(RealNumSpecloc, RealNumSpecglob, nSpecies, MPI_DOUBLE_PRECISION, MPI_SUM,  PartMPI%COMM, IERROR)
    RealNumSpecloc = RealNumSpecglob
#endif
    DO dir = 1,3
      IF (VeloDirs(dir) .OR. VeloDirs(4)) THEN
        DO iSpec = 1,nSpecies
          IF(RealNumSpecloc(iSpec).EQ.0)THEN
            PartVtrans(iSpec,dir) = 0.
          ELSE
            PartVtrans(iSpec,dir) = PartVtrans(iSpec,dir)/RealNumSpecloc(iSpec)
          END IF
        END DO ! iSpec = 1,nSpecies
      END IF
    END DO
  ELSE !no vMPF
#ifdef MPI
    CALL MPI_ALLREDUCE(NumSpecloc, NumSpecglob, nSpecies, MPI_INTEGER, MPI_SUM, PartMPI%COMM, IERROR)
    NumSpecloc = NumSpecglob
#endif
    DO dir = 1,3
      IF (VeloDirs(dir) .OR. VeloDirs(4)) THEN
        DO iSpec = 1,nSpecies
          IF(NumSpecloc(iSpec).EQ.0)THEN
            PartVtrans(iSpec,dir) = 0.
          ELSE
            PartVtrans(iSpec,dir) = PartVtrans(iSpec,dir)/NumSpecloc(iSpec)
          END IF
        END DO ! iSpec = 1,nSpecies
      END IF
    END DO
  END IF !usevMPF
  
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      DO dir = 1,3
        IF (VeloDirs(dir) .OR. VeloDirs(4)) THEN
          IF (usevMPF) THEN 
            PartVtherm(PartSpecies(i),dir) = PartVtherm(PartSpecies(i),dir) + PartMPF(i) * &
                (PartState(i,dir+3) - PartVtrans(PartSpecies(i),dir))*(PartState(i,dir+3) - PartVtrans(PartSpecies(i),dir))
          ELSE
            PartVtherm(PartSpecies(i),dir) = PartVtherm(PartSpecies(i),dir) + &
                (PartState(i,dir+3) - PartVtrans(PartSpecies(i),dir))*(PartState(i,dir+3) - PartVtrans(PartSpecies(i),dir))
          END IF
        END IF
      END DO
    END IF
  END DO
#ifdef MPI
  CALL MPI_ALLREDUCE(PartVtherm, PartVthermglob, 4*nSpecies, MPI_DOUBLE_PRECISION, MPI_SUM,  PartMPI%COMM, IERROR)
  PartVtherm = PartVthermglob
#endif
  DO dir = 1,3
    IF (VeloDirs(dir) .OR. VeloDirs(4)) THEN
      DO iSpec = 1,nSpecies
        IF (usevMPF) THEN
          IF(RealNumSpecloc(iSpec).EQ.0)THEN
            PartVtherm(iSpec,dir)=0.
          ELSE
            PartVtherm(iSpec,dir)=PartVtherm(iSpec,dir)/RealNumSpecloc(iSpec)
          END IF
        ELSE
          IF(NumSpecloc(iSpec).EQ.0)THEN
            PartVtherm(iSpec,dir)=0.
          ELSE
            PartVtherm(iSpec,dir)=PartVtherm(iSpec,dir)/NumSpecloc(iSpec)
          END IF
        END IF
      END DO ! iSpec = 1,nSpecies
    END IF
  END DO
 ! calc abolute value
  IF (VeloDirs(4)) THEN
    PartVtrans(:,4) = SQRT(PartVtrans(:,1)*PartVtrans(:,1) + PartVtrans(:,2)*PartVtrans(:,2) + PartVtrans(:,3)*PartVtrans(:,3))
    PartVtherm(:,4) = PartVtherm(:,1) + PartVtherm(:,2) + PartVtherm(:,3)
  END IF
  PartVtherm(:,:) = SQRT(PartVtherm(:,:))
END SUBROUTINE CalcVelocities


SUBROUTINE CalcIntTempsAndEn(IntTemp, IntEn,NumSpec)
!===================================================================================================================================
! Calculation of internal Temps (TVib, TRot, Telec) and gives back the global values
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars,      ONLY: PartSpecies, Species, PDM, nSpecies, BoltzmannConst, PartMPF, usevMPF
USE MOD_DSMC_Vars,          ONLY: PartStateIntEn, SpecDSMC, DSMC
USE MOD_DSMC_Analyze,       ONLY: CalcTVib, CalcTelec, CalcTVibPoly
USE MOD_Particle_MPI_Vars,  ONLY: PartMPI
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)               :: NumSpec(nSpecies+1)    ! number of real particles (already GLOBAL number)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: IntTemp(nSpecies,3) , IntEn(nSpecies,3) 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iPart, iSpec
REAL                           :: EVib(nSpecies), ERot(nSpecies), Eelec(nSpecies), tempVib, NumSpecTemp
#ifdef MPI
REAL                           :: RD(nSpecies)
#endif /*MPI*/
!===================================================================================================================================
EVib    = 0.
ERot    = 0.
Eelec   = 0.
! set electronic state to zero
IntTemp(:,3) = 0.

! Sum up internal energies
DO iPart=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    IF (usevMPF) THEN
      EVib(PartSpecies(iPart)) = EVib(PartSpecies(iPart)) + PartStateIntEn(iPart,1) * PartMPF(iPart)
      ERot(PartSpecies(iPart)) = ERot(PartSpecies(iPart)) + PartStateIntEn(iPart,2) * PartMPF(iPart)
      IF ( DSMC%ElectronicModel ) THEN
        Eelec(PartSpecies(iPart)) = Eelec(PartSpecies(iPart)) + PartStateIntEn(iPart,3) * PartMPF(iPart)
      END IF
    ELSE
      EVib(PartSpecies(iPart)) = EVib(PartSpecies(iPart)) + PartStateIntEn(iPart,1)
      ERot(PartSpecies(iPart)) = ERot(PartSpecies(iPart)) + PartStateIntEn(iPart,2)
      IF ( DSMC%ElectronicModel ) THEN
        Eelec(PartSpecies(iPart)) = Eelec(PartSpecies(iPart)) + PartStateIntEn(iPart,3)
      END IF
    END IF
  END IF
END DO

#ifdef MPI
IF(PartMPI%MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,EVib ,nSpecies, MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,ERot ,nSpecies, MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,Eelec,nSpecies, MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
ELSE
  CALL MPI_REDUCE(EVib        ,RD   ,nSpecies, MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
  CALL MPI_REDUCE(ERot        ,RD   ,nSpecies, MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
  CALL MPI_REDUCE(Eelec       ,RD   ,nSpecies, MPI_DOUBLE_PRECISION, MPI_SUM,0, PartMPI%COMM, IERROR)
END IF
#endif /*MPI*/

! final computation is only done for the root
IF(PartMPI%MPIRoot)THEN
  ! Calc TVib, TRot
  DO iSpec = 1, nSpecies
    NumSpecTemp = NumSpec(iSpec)
    IF(((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)).AND.(NumSpecTemp.GT.0.0)) THEN
      IF (SpecDSMC(iSpec)%PolyatomicMol.AND.(SpecDSMC(iSpec)%Xi_Rot.EQ.3)) THEN 
        IntTemp(iSpec,2) = 2.0*ERot(iSpec)/(3.0*BoltzmannConst*NumSpecTemp)  !Calc TRot          
      ELSE
        IntTemp(iSpec,2) = ERot(iSpec)/(BoltzmannConst*NumSpecTemp)  !Calc TRot
      END IF
      IF (EVib(iSpec)/NumSpecTemp.GT.SpecDSMC(iSpec)%EZeroPoint) THEN
        IF (SpecDSMC(iSpec)%PolyatomicMol) THEN
          IntTemp(iSpec,1) = CalcTVibPoly(EVib(iSpec)/NumSpecTemp, iSpec)
        ELSE
          IF (DSMC%VibEnergyModel.EQ.0) THEN              ! SHO-model
            tempVib = (EVib(iSpec)/(NumSpecTemp*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)-DSMC%GammaQuant)
            IF ((tempVib.GT.0.0) &
              .OR.(EVib(iSpec)/(NumSpecTemp*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib).LE.DSMC%GammaQuant)) THEN
              IntTemp(iSpec,1) = SpecDSMC(iSpec)%CharaTVib/LOG(1 + 1/(EVib(iSpec) & 
                                /(NumSpecTemp*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)-DSMC%GammaQuant))
            END IF
          ELSE                                            ! TSHO-model
            IntTemp(iSpec,1) = CalcTVib(SpecDSMC(iSpec)%CharaTVib, EVib(iSpec)/NumSpecTemp, SpecDSMC(iSpec)%MaxVibQuant)
          END IF
        END IF
      ELSE
        IntTemp(iSpec,1) = 0
      END IF
    ELSE
      IntTemp(iSpec,1) = 0
      IntTemp(iSpec,2) = 0
    END IF
    IF ( DSMC%ElectronicModel ) THEN
      IF ((NumSpecTemp.GT.0).AND.(SpecDSMC(iSpec)%InterID.NE.4) ) THEN
        IntTemp(iSpec,3) = CalcTelec(Eelec(iSpec)/NumSpecTemp,iSpec)
        IntEn(iSpec,3) = Eelec(iSpec)
      ELSE
        IntEn(iSpec,3) = 0.0
      END IF
    END IF
    IntEn(iSpec,1) = EVib(iSpec)
    IntEn(iSpec,2) = ERot(iSpec)
    IF(.NOT.usevMPF) IntEn(iSpec,:) = IntEn(iSpec,:) * Species(iSpec)%MacroParticleFactor
  END DO
END IF

END SUBROUTINE CalcIntTempsAndEn

#if (PP_TimeDiscMethod==42)
SUBROUTINE CollRates(CRate) 
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars,          ONLY: CollInf, DSMC
USE MOD_TimeDisc_Vars,      ONLY: dt
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: CRate(:) 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iCase
!===================================================================================================================================

  DO iCase=1, CollInf%NumCase + 1
    CRate(iCase) =  DSMC%NumColl(iCase) / dt
  END DO
  DSMC%NumColl = 0
END SUBROUTINE CollRates

SUBROUTINE ReacRates(RRate, NumSpec) 
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars,          ONLY: ChemReac, DSMC
USE MOD_TimeDisc_Vars,      ONLY: dt
USE MOD_Particle_Vars,      ONLY: Species, nSpecies
USE MOD_Particle_Mesh_Vars, ONLY: GEO
USE MOD_Particle_MPI_Vars,  ONLY: PartMPI
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: RRate(:) 
REAL,INTENT(IN)                 :: NumSpec(:) ! is the global number of real particles
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iReac
#ifdef MPI
REAL                            :: RD(1:ChemReac%NumOfReact)
#endif /*MPI*/
!===================================================================================================================================

#ifdef MPI
IF(PartMPI%MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE    ,ChemReac%NumReac,ChemReac%NumOfReact,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
ELSE
  CALL MPI_REDUCE(ChemReac%NumReac,RD              ,ChemReac%NumOfReact,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
END IF
#endif /*MPI*/

IF(PartMPI%MPIRoot)THEN
  DO iReac=1, ChemReac%NumOfReact
    IF ((NumSpec(ChemReac%DefinedReact(iReac,1,1)).GT.0).AND.(NumSpec(ChemReac%DefinedReact(iReac,1,2)).GT.0)) THEN
      SELECT CASE(TRIM(ChemReac%ReactType(iReac)))
      CASE('R')
        IF (DSMC%ReservoirRateStatistic) THEN ! Calculation of rate constant through actual number of allowed reactions
          RRate(iReac) = ChemReac%NumReac(iReac) * Species(ChemReac%DefinedReact(iReac,2,1))%MacroParticleFactor &
                     * GEO%MeshVolume**2 / (dt &
                     * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,1)) &
                     * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,2)) &
                     * Species(ChemReac%DefinedReact(iReac,1,3))%MacroParticleFactor * NumSpec(nSpecies+1))
        ! Calculation of rate constant through mean reaction probability (using mean reaction prob and sum of coll prob)
        ELSEIF(ChemReac%ReacCount(iReac).GT.0) THEN
          RRate(iReac) = ChemReac%NumReac(iReac) * ChemReac%ReacCollMean(iReac) * GEO%MeshVolume**2 &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor / (dt * ChemReac%ReacCount(iReac)             &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1))     &
               * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2))    &
               * Species(ChemReac%DefinedReact(iReac,1,3))%MacroParticleFactor*NumSpec(nSpecies+1))
        END IF
      CASE('D')
        IF (DSMC%ReservoirRateStatistic) THEN ! Calculation of rate constant through actual number of allowed reactions
          RRate(iReac) = ChemReac%NumReac(iReac) * Species(ChemReac%DefinedReact(iReac,2,1))%MacroParticleFactor &
                       * GEO%MeshVolume / (dt &
                       * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1)) &
                       * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2)))
        ! Calculation of rate constant through mean reaction probability (using mean reaction prob and sum of coll prob)
        ELSEIF(ChemReac%ReacCount(iReac).GT.0) THEN
          RRate(iReac) = ChemReac%NumReac(iReac) * ChemReac%ReacCollMean(iReac) &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor* GEO%MeshVolume / (dt * ChemReac%ReacCount(iReac) &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1))         &
               * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2)))
        END IF
      CASE('E')
        IF (DSMC%ReservoirRateStatistic) THEN ! Calculation of rate constant through actual number of allowed reactions
          RRate(iReac) = ChemReac%NumReac(iReac) * Species(ChemReac%DefinedReact(iReac,2,1))%MacroParticleFactor &
                       * GEO%MeshVolume / (dt &
                       * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1)) &
                       * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2)))
        ! Calculation of rate constant through mean reaction probability (using mean reaction prob and sum of coll prob)
        ELSEIF(ChemReac%ReacCount(iReac).GT.0) THEN 
          RRate(iReac) = ChemReac%NumReac(iReac) * ChemReac%ReacCollMean(iReac) &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor * GEO%MeshVolume / (dt * ChemReac%ReacCount(iReac) &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1))         &
               * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2)))
        END IF
      CASE('i')
        IF (DSMC%ReservoirRateStatistic) THEN ! Calculation of rate constant through actual number of allowed reactions
          RRate(iReac) = ChemReac%NumReac(iReac) * Species(ChemReac%DefinedReact(iReac,2,1))%MacroParticleFactor &
                       * GEO%MeshVolume / (dt &
                       * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1)) &
                       * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2)))
        ! Calculation of rate constant through mean reaction probability (using mean reaction prob and sum of coll prob)
        ELSEIF(ChemReac%ReacCount(iReac).GT.0) THEN
          RRate(iReac) = ChemReac%NumReac(iReac) * ChemReac%ReacCollMean(iReac) &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor * GEO%MeshVolume / (dt * ChemReac%ReacCount(iReac) &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1))         &
               * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2)))
        END IF
      CASE('r')
        IF (DSMC%ReservoirRateStatistic) THEN ! Calculation of rate constant through actual number of allowed reactions
          RRate(iReac) = ChemReac%NumReac(iReac) * Species(ChemReac%DefinedReact(iReac,2,1))%MacroParticleFactor &
                     * GEO%MeshVolume**2 / (dt &
                     * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,1)) &
                     * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,2)) &
!                     * Species(ChemReac%DefinedReact(iReac,1,3))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,3)) )
                     * Species(ChemReac%DefinedReact(iReac,1,3))%MacroParticleFactor * NumSpec(nSpecies+1))
        ! Calculation of rate constant through mean reaction probability (using mean reaction prob and sum of coll prob)
        ELSEIF(ChemReac%ReacCount(iReac).GT.0) THEN
          RRate(iReac) = ChemReac%NumReac(iReac) * ChemReac%ReacCollMean(iReac) * GEO%MeshVolume**2 &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor / (dt * ChemReac%ReacCount(iReac)             &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1))     &
               * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2))    &
               * Species(ChemReac%DefinedReact(iReac,1,3))%MacroParticleFactor*NumSpec(nSpecies+1))
!               * Species(ChemReac%DefinedReact(iReac,1,3))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,3)) )
        END IF
      CASE('x')
        IF (DSMC%ReservoirRateStatistic) THEN ! Calculation of rate constant through actual number of allowed reactions
          RRate(iReac) = ChemReac%NumReac(iReac) * Species(ChemReac%DefinedReact(iReac,2,1))%MacroParticleFactor &
                       * GEO%MeshVolume / (dt &
                       * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1)) &
                       * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2)))
        ! Calculation of rate constant through mean reaction probability (using mean reaction prob and sum of coll prob)
        ELSEIF(ChemReac%ReacCount(iReac).GT.0) THEN
      RRate(iReac) = ChemReac%NumReac(iReac) * ChemReac%ReacCollMean(iReac) &
           * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor * GEO%MeshVolume / (dt * ChemReac%ReacCount(iReac) &
           * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1))         &
           * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2)))
        END IF
      CASE('iQK')
        IF (DSMC%ReservoirRateStatistic) THEN ! Calculation of rate constant through actual number of allowed reactions
          RRate(iReac) = ChemReac%NumReac(iReac) * Species(ChemReac%DefinedReact(iReac,2,1))%MacroParticleFactor &
                       * GEO%MeshVolume / (dt &
                       * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1)) &
                       * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2)))
        ! Calculation of rate constant through mean reaction probability (using mean reaction prob and sum of coll prob)
        ELSEIF(ChemReac%ReacCount(iReac).GT.0) THEN
      RRate(iReac) = ChemReac%NumReac(iReac) * ChemReac%ReacCollMean(iReac) &
           * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor * GEO%MeshVolume / (dt * ChemReac%ReacCount(iReac) &
           * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1))         &
           * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2)))
        END IF
      CASE('rQK')
        IF (DSMC%ReservoirRateStatistic) THEN ! Calculation of rate constant through actual number of allowed reactions
          RRate(iReac) = ChemReac%NumReac(iReac) * Species(ChemReac%DefinedReact(iReac,2,1))%MacroParticleFactor &
                     * GEO%MeshVolume**2 / (dt &
                     * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,1)) &
                     * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,2)) &
                     * Species(ChemReac%DefinedReact(iReac,1,3))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,3)) )
        ! Calculation of rate constant through mean reaction probability (using mean reaction prob and sum of coll prob)
        ELSEIF(ChemReac%ReacCount(iReac).GT.0) THEN
          RRate(iReac) = ChemReac%NumReac(iReac) * ChemReac%ReacCollMean(iReac) * GEO%MeshVolume**2 &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor / (dt * ChemReac%ReacCount(iReac)             &
               * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,1))     &
               * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,2))    &
               * Species(ChemReac%DefinedReact(iReac,1,3))%MacroParticleFactor*NumSpec(ChemReac%DefinedReact(iReac,1,3)) )
        END IF
      END SELECT
    END IF
  END DO
END IF
ChemReac%NumReac = 0.

END SUBROUTINE ReacRates
#endif 

#if ( PP_TimeDiscMethod == 42)
SUBROUTINE ElectronicTransition (  Time, NumSpec )
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,          ONLY: DSMC, SpecDSMC
  USE MOD_TimeDisc_Vars,      ONLY: dt
  USE MOD_Particle_Vars,      ONLY: nSpecies, Species
  USE MOD_Particle_Mesh_Vars, ONLY: GEO
  ! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  REAL,INTENT(IN)                :: Time
  REAL,INTENT(IN)               :: NumSpec(:)
  INTEGER                        :: iSpec, iSpec2, iQua1, iQua2, MaxElecQua
! accary of kf
!===================================================================================================================================

  IF ( DSMC%ElectronicModel ) THEN
  ! kf = d n_of_N^i / dt / ( n_of_N^i n_of_M) )
    DO iSpec = 1, nSpecies
      IF ( SpecDSMC(iSpec)%InterID .ne. 4 ) THEN
        DO iSpec2 = 1, nSpecies
     ! calculaction of kf for each reaction
!       MaxElecQua = SpecDSMC(iSpec)%MaxElecQuant
          MaxElecQua = 2
        ! for first tests only consider the first 10 transition levels
          DO iQua1 = 0, MaxElecQua
            DO iQua2 = 0, MaxElecQua
          ! calculate kf
          ! kf = ( d n_of_N^i / d t )  / ( n_of_N^i n_of_M )
              IF ( (NumSpec(iSpec2) .ne. 0) .and. (SpecDSMC(iSpec)%levelcounter(iQua1) .ne. 0 ) ) THEN
                SpecDSMC(iSpec)%ElectronicTransition(iSpec2,iQua1,iQua2) = &
                                      SpecDSMC(iSpec)%ElectronicTransition(iSpec2,iQua1,iQua2) * GEO%MeshVolume &
                                    / ( dt * SpecDSMC(iSpec)%levelcounter(iQua1) * NumSpec(iSpec2) *           &
                                        Species(iSpec2)%MacroParticleFactor )
              END IF
!             print*,SpecDSMC(iSpec)%ElectronicTransition(iSpec2,iQua1,iQua2)
            END DO
          END DO
        END DO
      END IF
    END DO
    CALL WriteEletronicTransition( Time )
    ! nullyfy
    DO iSpec = 1, nSpecies
      IF(SpecDSMC(iSpec)%InterID.ne.4)THEN
        SpecDSMC(iSpec)%ElectronicTransition = 0
      END IF
    END DO
  END IF
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE
#endif

#if ( PP_TimeDiscMethod == 42)
SUBROUTINE WriteEletronicTransition ( Time )
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
  USE MOD_Globals,            ONLY: Getfreeunit
  USE MOD_DSMC_Vars,          ONLY: DSMC, SpecDSMC
  USE MOD_Particle_Vars,      ONLY: nSpecies
!   ! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  REAL,INTENT(IN)                :: Time
  INTEGER                        :: iSpec, iSpec2, iunit, iQua1, iQua2, MaxElecQua, ii
  CHARACTER(LEN=128)             :: FileNameTransition
  LOGICAL                        :: bExist
! accary of kf
!-----------------------------------------------------------------------------------------------------------------------------------
  bExist = .false.
  IF ( DSMC%ElectronicModel ) THEN
  ! kf = d n_of_N^i / dt / ( n_of_N^i n_of_M) )
    DO iSpec = 1, nSpecies
      IF (SpecDSMC(iSpec)%InterID .ne. 4 ) THEN
!        MaxElecQua = SpecDSMC(iSpec)%MaxElecQuant
        MaxElecQua = 2
        ! output to transition file
        FileNameTransition = trim(SpecDSMC(iSpec)%Name)//'_Transition.csv'
        INQUIRE( FILE = FileNameTransition, EXIST=bExist)
!-----------------------------------------------------------------------------------------------------------------------------------
        IF ( bExist .EQV. .false. ) THEN 
          iunit=GETFREEUNIT()
          OPEN(UNIT=iunit,FILE=FileNameTransition,FORM='FORMATTED',STATUS='UNKNOWN')
!         ! writing header
          WRITE(iunit,'(A6,A5)',ADVANCE='NO') 'TIME', ' '
          ii = 2
          DO iSpec2 = 1, nSpecies
            DO iQua1 = 0, MaxElecQua
              DO iQua2 = 0, MaxElecQua
                WRITE(iunit,'(I3.3,A,I2.2,A,I2.2,A,I2.2,A)', ADVANCE='NO') ii,'_Species',iSpec2,'_',iQua1,'_to_',iQua2,'  '
                WRITE(iunit,'(A1)',ADVANCE='NO') ','
                ii = ii + 1
              END DO
            END DO
          END DO
        ELSE
!-----------------------------------------------------------------------------------------------------------------------------------
!         ! writing header
          iunit=GETFREEUNIT()
          OPEN(unit=iunit,FILE=FileNameTransition,FORM='Formatted',POSITION='APPEND',STATUS='old')
          WRITE(iunit,104,ADVANCE='NO') TIME
          DO iSpec2 = 1, nSpecies
!       ! calculaction of kf for each reaction
            DO iQua1 = 0, MaxElecQua
              DO iQua2 = 0, MaxElecQua
              ! write values to file
                WRITE(iunit,'(A1)',ADVANCE='NO') ','
                WRITE(iunit,104,ADVANCE='NO') SpecDSMC(iSpec)%ElectronicTransition(iSpec2,iQua1,iQua2)
              END DO
            END DO
          END DO
          WRITE(iunit,'(A)') ' '
        END IF
        CLOSE(unit = iunit)
      END IF
    END DO
!-----------------------------------------------------------------------------------------------------------------------------------
  END IF
!-----------------------------------------------------------------------------------------------------------------------------------
104    FORMAT (e25.14)
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE WriteEletronicTransition
#endif

SUBROUTINE TrackingParticlePosition(time) 
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars,          ONLY:PartState, PDM, PEM
USE MOD_Particle_MPI_Vars,      ONLY:PartMPI
USE MOD_Particle_Analyze_Vars,  ONLY:printDiff,printDiffVec,printDiffTime
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: time
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,iunit,iPartState
CHARACTER(LEN=60) :: TrackingFilename!,hilf
LOGICAL            :: fexist
REAL               :: diffPos,diffVelo
!===================================================================================================================================

iunit=GETFREEUNIT()
!WRITE(UNIT=hilf,FMT='(I6.6)') MyRank
!TrackingFilename = ('MyRank'//TRIM(hilf)//'_ParticlePosition.csv')
TrackingFilename = ('ParticlePosition.csv')

INQUIRE(FILE = TrackingFilename, EXIST=fexist)
IF(.NOT.fexist) THEN 
 IF(PartMPI%MPIRoot)THEN
   iunit=GETFREEUNIT()
   OPEN(UNIT=iunit,FILE=TrackingFilename,FORM='FORMATTED',STATUS='UNKNOWN')
   !CALL FLUSH (iunit)
    ! writing header
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'TIME', ' '
    WRITE(iunit,'(A1)',ADVANCE='NO') ','
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartNum', ' '
    WRITE(iunit,'(A1)',ADVANCE='NO') ','
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartPosX', ' '
    WRITE(iunit,'(A1)',ADVANCE='NO') ','
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartPosY', ' '
    WRITE(iunit,'(A1)',ADVANCE='NO') ','
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartPosZ', ' '
    WRITE(iunit,'(A1)',ADVANCE='NO') ','
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartVelX', ' '
    WRITE(iunit,'(A1)',ADVANCE='NO') ','
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartVelY', ' '
    WRITE(iunit,'(A1)',ADVANCE='NO') ','
    WRITE(iunit,'(A8,A5)',ADVANCE='NO') 'PartVelZ', ' '
    CLOSE(iunit)
  END IF
ELSE
  iunit=GETFREEUNIT()
  OPEN(unit=iunit,FILE=TrackingFileName,FORM='Formatted',POSITION='APPEND',STATUS='old')
  !CALL FLUSH (iunit)
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      WRITE(iunit,104,ADVANCE='NO') TIME
      WRITE(iunit,'(A1)',ADVANCE='NO') ','
      WRITE(iunit,'(I12)',ADVANCE='NO') i
      DO iPartState=1,6
        WRITE(iunit,'(A1)',ADVANCE='NO') ','
        WRITE(iunit,104,ADVANCE='NO') PartState(i,iPartState)
      END DO
      WRITE(iunit,'(A1)',ADVANCE='NO') ','
      WRITE(iunit,'(I12)',ADVANCE='NO') PEM%Element(i)
      WRITE(iunit,'(A)') ' '
     END IF
  END DO
  CLOSE(iunit)
END IF
IF (printDiff) THEN
  diffPos=0.
  diffVelo=0.
  IF (time.GE.printDiffTime) THEN
    printDiff=.FALSE.
    DO iPartState=1,3
      diffPos=diffPos+(printDiffVec(iPartState)-PartState(1,iPartState))**2
      diffVelo=diffVelo+(printDiffVec(iPartState+3)-PartState(1,iPartState+3))**2
    END DO
    WRITE(*,'(A,e24.14,x,e24.14)') 'L2-norm from printDiffVec: ',SQRT(diffPos),SQRT(diffVelo)
  END IF
END IF
104    FORMAT (e25.14)

END SUBROUTINE TrackingParticlePosition

Function CalcEkinPart(iPart) 
!===================================================================================================================================
! computes the kinetic energy of one particle
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Equation_Vars,          ONLY : c2, c2_inv
USE MOD_Particle_Vars,          ONLY : PartState, PartSpecies, Species
USE MOD_PARTICLE_Vars,          ONLY : PartMPF, usevMPF
USE MOD_Particle_Vars,          ONLY : PartLorentzType 
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(IN)                 :: iPart
REAL                               :: CalcEkinPart
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                               :: partV2, gamma1, Ekin
!===================================================================================================================================

IF (PartLorentzType.EQ.5)THEN
  ! gamma v is pushed instead of gamma, therefore, only the relativistic kinetic energy is computed
  ! compute gamma
  gamma1=SQRT(1.0+DOT_PRODUCT(PartState(iPart,4:6),PartState(iPart,4:6))*c2_inv)
  IF(usevMPF)THEN
    Ekin=PartMPF(iPart)*(gamma1-1.0)*Species(PartSpecies(iPart))%MassIC*c2
  ELSE
    Ekin= (gamma1-1.0)* Species(PartSpecies(iPart))%MassIC*Species(PartSpecies(iPart))%MacroParticleFactor*c2
  END IF
ELSE
  partV2 = PartState(iPart,4) * PartState(iPart,4) &
         + PartState(iPart,5) * PartState(iPart,5) &
         + PartState(iPart,6) * PartState(iPart,6)
  
  IF(usevMPF)THEN
    IF (partV2.LT.1e6)THEN
      Ekin= 0.5 * Species(PartSpecies(iPart))%MassIC * partV2 * PartMPF(iPart)            
    ELSE
      gamma1=partV2*c2_inv
      gamma1=1.0/SQRT(1.-gamma1)
      Ekin=PartMPF(iPart)*(gamma1-1.0)*Species(PartSpecies(iPart))%MassIC*c2
    END IF ! ipartV2
  ELSE ! novMPF
    IF (partV2.LT.1e6)THEN
      Ekin= 0.5*Species(PartSpecies(iPart))%MassIC*partV2* Species(PartSpecies(iPart))%MacroParticleFactor
    ELSE
      gamma1=partV2*c2_inv
      gamma1=1.0/SQRT(1.-gamma1)
      Ekin= (gamma1-1.0)* Species(PartSpecies(iPart))%MassIC*Species(PartSpecies(iPart))%MacroParticleFactor*c2
    END IF ! ipartV2
  END IF ! usevMPF
END IF
CalcEkinPart=Ekin
END FUNCTION CalcEkinPart
 
SUBROUTINE FinalizeParticleAnalyze()
!===================================================================================================================================
! Finalizes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Particle_Analyze_Vars,ONLY:ParticleAnalyzeInitIsDone
! IMPLICIT VARIABLE HANDLINGDGInitIsDone
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
ParticleAnalyzeInitIsDone = .FALSE.
END SUBROUTINE FinalizeParticleAnalyze
#endif /*PARTICLES*/


END MODULE MOD_Particle_Analyze
