#include "boltzplatz.h"

MODULE MOD_Particle_Analyze
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
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

INTERFACE CalcPotentialEnergy
  MODULE PROCEDURE CalcPotentialEnergy
END INTERFACE

INTERFACE CalcShapeEfficiencyR
  MODULE PROCEDURE CalcShapeEfficiencyR
END INTERFACE

PUBLIC:: InitParticleAnalyze, FinalizeParticleAnalyze, AnalyzeParticles, CalcKineticEnergy, CalcPotentialEnergy
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
USE MOD_Particle_Analyze_Vars  !,ONLY:ParticleAnalyzeInitIsDone, CalcCharge, CalcEkin, CalcEpot, DoAnalyze
USE MOD_ReadInTools             ,ONLY: GETLOGICAL, GETINT, GETSTR
USE MOD_Particle_Vars           ,ONLY: nSpecies

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!CHARACTER(LEN=40)                :: DefStr
!===================================================================================================================================
IF (ParticleAnalyzeInitIsDone) THEN
  CALL abort(__STAMP__,'InitParticleAnalyse already called.',999,999.)
  RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE ANALYZE...'

 DoAnalyze = .FALSE.
 CalcCharge = GETLOGICAL('CalcCharge','.FALSE.')
 IF(CalcCharge) DoAnalyze = .TRUE. 
 CalcEpot = GETLOGICAL('CalcPotentialEnergy','.FALSE.')
 IF(CalcEpot) DoAnalyze = .TRUE.
 CalcEkin = GETLOGICAL('CalcKineticEnergy','.FALSE.')
 CalcTemp = GETLOGICAL('CalcTransTemp','.FALSE.')
 IF (CalcTemp) CalcEkin = .TRUE.
 IF(CalcEkin) THEN
   DoAnalyze = .TRUE.
   IF (nSpecies .GT. 1) THEN
    nEkin = nSpecies + 1
   ELSE
    nEkin = 1
   END IF
 END IF
 CalcNumSpec = GETLOGICAL('CalcNumSpec','.FALSE.')
 IF(CalcNumSpec) DoAnalyze = .TRUE.
 CalcShapeEfficiency = GETLOGICAL('CalcShapeEfficiency','.FALSE.')
 IF (CalcShapeEfficiency) THEN
   DoAnalyze = .TRUE.
   CalcShapeEfficiencyMethod = GETSTR('CalcShapeEfficiencyMethod','AllParts')
   SELECT CASE(CalcShapeEfficiencyMethod)
   CASE('AllParts')  ! All currently available Particles are used
   CASE('SomeParts') ! A certain percentage of currently available Particles is used
     ShapeEfficiencyNumber = GETINT('ShapeEfficiencyNumber','100')  ! in percent
   CASE DEFAULT
     SWRITE(*,*) 'ERROR: CalcShapeEfficiencyMethod',CalcShapeEfficiencyMethod,'does not exist!'
     STOP
   END SELECT
 END IF

IsRestart = GETLOGICAL('IsRestart','.FALSE.')


PartAnalyzeStep = GETINT('Part-AnalyzeStep','1')
IF (PartAnalyzeStep.EQ.0) PartAnalyzeStep = 123456789


ChargeCalcDone = .FALSE.

ParticleAnalyzeInitIsDone=.TRUE.

SWRITE(UNIT_stdOut,'(A)')' INIT PARTCILE ANALYZE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitParticleAnalyze



SUBROUTINE AnalyzeParticles(Time)
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars!,ONLY:ParticleAnalyzeInitIsDone, CalcCharge, CalcEkin, CalcEpot, DoAnalyze, IsRestart
USE MOD_PARTICLE_Vars,    ONLY: nSpecies
USE MOD_DSMC_Vars,        ONLY: DSMC, CollInf, useDSMC, CollisMode, ChemReac, SpecDSMC
USE MOD_Restart_Vars, ONLY: DoRestart
#ifdef MPI
  USE MOD_part_MPI_Vars,      ONLY : PMPIVAR
#endif
#if ( PP_TimeDiscMethod ==42)
USE MOD_Particle_Vars,         ONLY : Species
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
LOGICAL             :: isOpen, FileExists                                          !
CHARACTER(LEN=350)  :: outfile ,hilf                                               !
INTEGER             :: unit_index, unit_index_2, iSpec,  jSpec, iCase
REAL                :: WEl, WMag
REAL                :: Ekin(nSpecies + 1), Temp(nSpecies), IntTemp(nSpecies,3), IntEn(nSpecies,3)
INTEGER             :: NumSpec(nSpecies), OutputCounter, iTvib
#ifdef MPI
INTEGER             :: sumNumSpec(nSpecies)
REAL                :: sumTemp(nSpecies), sumIntTemp(nSpecies),sumIntEn(nSpecies), sumEkin(nSpecies + 1)
REAL                :: sumWEl, sumWMag
#endif
REAL, ALLOCATABLE   :: CRate(:), RRate(:)
#if (PP_TimeDiscMethod ==42)
INTEGER           :: ii, iunit
CHARACTER(LEN=64) :: DebugElectronicStateFilename
#endif 
!===================================================================================================================================
IF ( DoRestart ) THEN
  isRestart = .true.
END IF
IF (DoAnalyze) THEN
!SWRITE(UNIT_StdOut,'(132("-"))')
!SWRITE(UNIT_stdOut,'(A)') ' PERFORMING PARTICLE ANALYZE...'
IF (useDSMC) THEN
  ALLOCATE(CRate(CollInf%NumCase + 1))
  IF (CollisMode.EQ.3) ALLOCATE(RRate(ChemReac%NumOfReact))
END IF
OutputCounter = 2
unit_index = 535
#ifdef MPI
 IF (PMPIVAR%iProc .EQ. 0) THEN
#endif    /* MPI */
    INQUIRE(UNIT   = unit_index , OPENED = isOpen)
    IF (.NOT.isOpen) THEN
#if (PP_TimeDiscMethod==42)
      ! if only the reaction rate is desired (resevoir) the initial temperature
      ! of the second species is added to the filename
       IF ( DSMC%ReservoirSimuRate .EQV. .true. ) THEN
        IF ( SpecDSMC(1)%InterID .eq. 2 .or. SpecDSMC(1)%InterID .eq. 20 ) THEN
          iTvib = INT(SpecDSMC(1)%TVib)
          WRITE( hilf, '(I5.5)') iTvib
          outfile = 'Database_Tvib_'//TRIM(hilf)//'.csv'
        ELSE
          iTvib = INT(Species(1)%MWTemperatureIC )
          WRITE( hilf, '(I5.5)') iTvib
          outfile = 'Database_Ttr_'//TRIM(hilf)//'.csv'
        END IF
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
            DO iSpec = 1, nSpecies
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(A14,I3.3)',ADVANCE='NO') 'PartNum-Spec-', iSpec
            END DO              
          END IF
          IF (CalcEpot) THEN 
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(A14)',ADVANCE='NO') 'W-El'
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,'(A14)',ADVANCE='NO') 'W-Mag'
          END IF
          IF (CalcEkin) THEN
            DO iSpec=1, nEkin
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' Ekin',iSpec,' '
              OutputCounter = OutputCounter + 1
            END DO
            IF (CalcTemp) THEN
              DO iSpec=1, nSpecies
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' Temp',iSpec,' '
                OutputCounter = OutputCounter + 1
              END DO
            END IF
#if (PP_TimeDiscMethod==42)
              IF (CollisMode.NE.1) THEN
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
                IF ( DSMC%ElectronicState ) THEN
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
                IF ( DSMC%ElectronicState ) THEN
                  DO iSpec=1, nSpecies
                    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                    WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' TempElec',iSpec,' '
                    OutputCounter = OutputCounter + 1
                  END DO
                END IF
              END IF
              WRITE(unit_index,'(A1)',ADVANCE='NO') ','
              WRITE(unit_index,'(I3.3,A,A5)',ADVANCE='NO') OutputCounter,' Pmax' ,' '
              OutputCounter = OutputCounter + 1                       
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
              IF (CollisMode.EQ.3) THEN
                DO iCase=1, ChemReac%NumOfReact 
                  WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                  WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' Reaction', iCase,' '
                  OutputCounter = OutputCounter + 1
                END DO
              END IF
              DO iSpec = 1, nSpecies
                WRITE(unit_index,'(A1)',ADVANCE='NO') ','
                WRITE(unit_index,'(I3.3,A,I3.3,A5)',ADVANCE='NO') OutputCounter,' PartNum', iSpec,' '
                OutputCounter = OutputCounter + 1
              END DO              
#endif
          END IF
          WRITE(unit_index,'(A14)') ' ' 
       END IF
    END IF
#ifdef MPI
 END IF
#endif    /* MPI */


!IF (CalcCharge.AND.(.NOT.ChargeCalcDone)) CALL CalcDepositedCharge()
IF (CalcEpot) THEN
  CALL CalcPotentialEnergy(WEl,WMag)
END IF
IF (CalcEkin) THEN
  CALL CalcKineticEnergy(Ekin)
END IF
IF (CalcTemp) THEN
  CALL CalcTemperature(Temp, NumSpec)
END IF
! MPI Communication
#ifdef MPI
IF (CalcEpot) THEN
  CALL MPI_REDUCE(WEl , sumWEl , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, PMPIVAR%COMM, IERROR)
  WEl = sumWEl
  CALL MPI_REDUCE(WMag, sumWMag , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, PMPIVAR%COMM, IERROR)
  WMag = sumWMag
END IF
IF (CalcEkin) THEN
  CALL MPI_REDUCE(Ekin(:) , sumEkin , nEkin , MPI_DOUBLE_PRECISION, MPI_SUM,0, PMPIVAR%COMM, IERROR)
  Ekin(:) = sumEkin(:)
END IF
IF (CalcTemp) THEN
  sumTemp = 0
  CALL MPI_REDUCE(Temp   , sumTemp    , nSpecies, MPI_DOUBLE_PRECISION, MPI_SUM,0, PMPIVAR%COMM, IERROR)
  Temp = sumTemp / PMPIVAR%nProcs
END IF
#endif

#if (PP_TimeDiscMethod==42)
CALL CollRates(CRate)
IF (CollisMode.EQ.3) CALL ReacRates(RRate, NumSpec)
IF (CollisMode.NE.1) CALL CalcIntTempsAndEn(IntTemp, IntEn)
CALL GetNumSpec(NumSpec)
! currently, calculation of internal electronic energy not implemented !
#ifdef MPI
! average over all cells
  IF (CollisMode.NE.1) THEN
    CALL MPI_REDUCE(IntTemp(:,1), sumIntTemp , nSpecies , MPI_DOUBLE_PRECISION, MPI_SUM,0, PMPIVAR%COMM, IERROR)
    IntTemp(:,1) = sumIntTemp / PMPIVAR%nProcs
    CALL MPI_REDUCE(IntTemp(:,2), sumIntTemp , nSpecies , MPI_DOUBLE_PRECISION, MPI_SUM,0,PMPIVAR%COMM, IERROR)
    IntTemp(:,2) = sumIntTemp / PMPIVAR%nProcs
    CALL MPI_REDUCE(IntEn(:,1) , sumIntEn   , nSpecies , MPI_DOUBLE_PRECISION, MPI_SUM,0, PMPIVAR%COMM, IERROR)
    IntEn(:,1) = sumIntEn 
    CALL MPI_REDUCE(IntEn(:,2) , sumIntEn   , nSpecies , MPI_DOUBLE_PRECISION, MPI_SUM,0, PMPIVAR%COMM, IERROR)
    IntEn(:,2) = sumIntEn 
    IF ( DSMC%ElectronicState ) THEN
      CALL MPI_REDUCE(IntTemp(:,3), sumIntTemp , nSpecies , MPI_DOUBLE_PRECISION, MPI_SUM,0, PMPIVAR%COMM, IERROR)
      IntTemp(:,3) = sumIntTemp / PMPIVAR%nProcs
      CALL MPI_REDUCE(IntEN(:,3), sumIntEn , nSpecies , MPI_DOUBLE_PRECISION, MPI_SUM,0, PMPIVAR%COMM, IERROR)
      IntEn(:,3) = sumIntEn
    END IF
  END IF 
  CALL MPI_REDUCE(NumSpec, sumNumSpec, nSpecies, MPI_INTEGER, MPI_SUM,0, PMPIVAR%COMM, IERROR)
  NumSpec = sumNumSpec
#endif /*MPI*/
#endif /*PP_TimeDiscMethod==42*/

IF (CalcShapeEfficiency) CALL CalcShapeEfficiencyR()   ! This will NOT be placed in the file but directly in "out"

#ifdef MPI
 IF (PMPIVAR%iProc .EQ. 0) THEN
#endif    /* MPI */
   WRITE(unit_index,104,ADVANCE='NO') Time
   IF (CalcNumSpec) THEN
     DO iSpec=1, nSpecies
       WRITE(unit_index,'(A1)',ADVANCE='NO') ','
       WRITE(unit_index,'(I10.1)',ADVANCE='NO') NumSpec(iSpec)
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
     IF (CalcTemp) THEN
       DO iSpec=1, nSpecies
         WRITE(unit_index,'(A1)',ADVANCE='NO') ','
         WRITE(unit_index,104,ADVANCE='NO') Temp(iSpec)
       END DO
     END IF
#if (PP_TimeDiscMethod==42)
      IF (CollisMode.NE.1) THEN
        DO iSpec=1, nSpecies
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,104,ADVANCE='NO') IntEn(iSpec,1)
        END DO
        DO iSpec=1, nSpecies
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,104,ADVANCE='NO') IntEn(iSpec,2)
        END DO
        IF ( DSMC%ElectronicState ) THEN
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
        IF ( DSMC%ElectronicState ) THEN
          DO iSpec=1, nSpecies
          ! currently set to one
            WRITE(unit_index,'(A1)',ADVANCE='NO') ','
            WRITE(unit_index,104,ADVANCE='NO') IntTemp(iSpec,3)
          END DO
        END IF
      END IF
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,104,ADVANCE='NO') DSMC%CollProbMax
      DO iCase=1, CollInf%NumCase +1 
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,104,ADVANCE='NO') CRate(iCase)
      END DO
      IF (CollisMode.EQ.3) THEN
        DO iCase=1, ChemReac%NumOfReact 
          WRITE(unit_index,'(A1)',ADVANCE='NO') ','
          WRITE(unit_index,104,ADVANCE='NO') RRate(iCase)
        END DO
      END IF
      DO iSpec=1, nSpecies
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,'(I10.1)',ADVANCE='NO') NumSpec(iSpec)
      END DO
#endif
   END IF
   WRITE(unit_index,'(A1)') ' ' 
#ifdef MPI
 END IF
#endif    /* MPI */

104    FORMAT (e25.14)

!SWRITE(UNIT_stdOut,'(A)')' PARTCILE ANALYZE DONE!'
!SWRITE(UNIT_StdOut,'(132("-"))')
ELSE
!SWRITE(UNIT_stdOut,'(A)')' NO PARTCILE ANALYZE TO DO!'
!SWRITE(UNIT_StdOut,'(132("-"))')
END IF ! DoAnalyze

#if ( PP_TimeDiscMethod ==42 )
IF ( DSMC%ElectronicState ) THEN
  CALL ElectronicTransition( Time, NumSpec )
END IF
#endif
#if ( PP_TimeDiscMethod ==42 )
  ! Debug Output for initialized electronic state
IF ( DSMC%ElectronicState ) THEN
  DO iSpec = 1, nSpecies
    IF ( SpecDSMC(iSpec)%InterID .ne. 4 ) THEN
      IF (  SpecDSMC(iSpec)%levelcounter(0) .ne. 0) THEN
        WRITE(DebugElectronicStateFilename,'(I2.2)') iSpec
        iunit = 485
        DebugElectronicStateFilename = 'End_Electronic_State_Species_'//trim(DebugElectronicStateFilename)//'.dat'
        OPEN(unit=iunit,file=DebugElectronicStateFilename,form='formatted',status='unknown')
        DO ii = 0, SpecDSMC(iSpec)%MaxElecQuant - 1
          WRITE(iunit,'(I3.1,3x,F12.7)') ii, REAL( SpecDSMC(iSpec)%levelcounter(ii) ) / REAL( Species(iSpec)%initialParticleNumber)
        END DO
        close(iunit)
      END IF
    END IF
  END DO
END IF
#endif
END SUBROUTINE AnalyzeParticles

SUBROUTINE CalcShapeEfficiencyR()
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Particle_Analyze_Vars,   ONLY : CalcShapeEfficiencyMethod, ShapeEfficiencyNumber
USE MOD_Mesh_Vars,               ONLY : nElems, Elem_xGP
USE MOD_PICDepo_Vars
USE MOD_Particle_Vars
USE MOD_PreProc
USE MOD_part_MPI_Vars,            ONLY : PMPIVAR
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                     :: NbrOfComps, NbrWithinRadius, NbrOfElems, NbrOfElemsWithinRadius
REAL                     :: RandVal1, RandVal3(1:3)
REAL, ALLOCATABLE        :: RandomParts(:,:)
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
      lmax = MIN(lmax,GEO%FIBGMkmax)
      lmin = INT((PartState(i,2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
      lmin = MAX(lmin,GEO%FIBGMkmin)
      mmax = INT((PartState(i,3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
      mmax = MIN(mmax,GEO%FIBGMlmax)
      mmin = INT((PartState(i,3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
      mmin = MAX(mmin,GEO%FIBGMlmin)
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
  WRITE(*,*) 'ShapeEfficiency (Proc,%,%Elems)',PMPIVAR%iProc,100*NbrWithinRadius/NbrOfComps,100*NbrOfElemsWithinRadius/NbrOfElems
  WRITE(*,*) 'ShapeEfficiency (Elems) for Proc',PMPIVAR%iProc,'is',100*NbrOfElemsWithinRadius/NbrOfElems,'%'
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
        lmax = MIN(lmax,GEO%FIBGMkmax)
        lmin = INT((PartState(i,2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
        lmin = MAX(lmin,GEO%FIBGMkmin)
        mmax = INT((PartState(i,3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
        mmax = MIN(mmax,GEO%FIBGMlmax)
        mmin = INT((PartState(i,3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
        mmin = MAX(mmin,GEO%FIBGMlmin)
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
  WRITE(*,*) 'ShapeEfficiency (Proc,%,%Elems)',PMPIVAR%iProc,100*NbrWithinRadius/NbrOfComps,100*NbrOfElemsWithinRadius/NbrOfElems
  WRITE(*,*) 'ShapeEfficiency (Elems) for Proc',PMPIVAR%iProc,'is',100*NbrOfElemsWithinRadius/NbrOfElems,'%'
#else
  WRITE(*,*) 'ShapeEfficiency (%,%Elems)',100*NbrWithinRadius/NbrOfComps,100*NbrOfElemsWithinRadius/NbrOfElems
  WRITE(*,*) 'ShapeEfficiency (Elems) is',100*NbrOfElemsWithinRadius/NbrOfElems,'%'
#endif
END IF
END SELECT
END SUBROUTINE CalcShapeEfficiencyR

SUBROUTINE GetNumSpec(NumSpec)
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Equation_Vars,      ONLY : c2 
USE MOD_Particle_Vars,      ONLY : PartSpecies, Species, PDM, nSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(OUT)            :: NumSpec(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i, iSpec
!===================================================================================================================================
  NumSpec = 0
! Sum up velocity 
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      NumSpec(PartSpecies(i)) = NumSpec(PartSpecies(i)) + 1
    END IF
  END DO
END SUBROUTINE GetNumSpec


SUBROUTINE CalcPotentialEnergy(WEl, WMag) 
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,          ONLY : nElems, sJ
USE MOD_Interpolation_Vars, ONLY : wGP
USE MOD_Equation_Vars,      ONLY : smu0, eps0 
USE MOD_DG_Vars,ONLY:U
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: WEl, WMag 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iElem
INTEGER           :: i,j,k
REAL              :: J_N(1,0:PP_N,0:PP_N,0:PP_N)
REAL              :: WEl_tmp, WMag_tmp, E_abs, B_abs 
!===================================================================================================================================

Wel=0.
WMag=0.

DO iElem=1,nElems
  !--- Calculate and save volume of element iElem
  WEl_tmp=0. 
  WMag_tmp=0. 
  J_N(1,0:PP_N,0:PP_N,0:PP_N)=1./sJ(:,:,:,iElem)
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
! in electromagnetische felder by henke 2011 - springer
! WMag = 1/(2mu) * int_V B^2 dV 
    E_abs = U(1,i,j,k,iElem)*U(1,i,j,k,iElem) &
          + U(2,i,j,k,iElem)*U(2,i,j,k,iElem) &
          + U(3,i,j,k,iElem)*U(3,i,j,k,iElem)
    B_abs = U(4,i,j,k,iElem)*U(4,i,j,k,iElem) &
          + U(5,i,j,k,iElem)*U(5,i,j,k,iElem) &
          + U(6,i,j,k,iElem)*U(6,i,j,k,iElem)
    WEl_tmp  = WEl_tmp  + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * E_abs 
    WMag_tmp = WMag_tmp + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * B_abs
  END DO; END DO; END DO
  WEl = WEl + WEl_tmp
  WMag = WMag + WMag_tmp
END DO

WEl = WEl * eps0 * 0.5 
WMag = WMag * smu0 * 0.5

END SUBROUTINE CalcPotentialEnergy



SUBROUTINE CalcKineticEnergy(Ekin) 
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Equation_Vars,          ONLY : c2 
USE MOD_Particle_Vars,          ONLY : PartState, PartSpecies, Species, PDM
USE MOD_PARTICLE_Vars,          ONLY : nSpecies, PartMPF, usevMPF
USE MOD_Particle_Analyze_Vars,  ONLY : nEkin
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Ekin(:) 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i
REAL(KIND=8)              :: partV2, Gamma  
!===================================================================================================================================

Ekin = 0.
IF (nEkin .GT. 1 ) THEN
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      partV2 = PartState(i,4) * PartState(i,4) &
              + PartState(i,5) * PartState(i,5) &
              + PartState(i,6) * PartState(i,6)
      IF ( partV2 .LT. 1e6) THEN  ! |v| < 1000
  !       Ekin = Ekin + 0.5 *  Species(PartSpecies(i))%MassIC * partV2 &
  !                     * PartMPF(i)            
        IF(usevMPF) THEN
          Ekin(nSpecies+1) = Ekin(nSpecies+1) + 0.5 *  Species(PartSpecies(i))%MassIC * partV2 &
                                            * PartMPF(i)            
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + 0.5 *  Species(PartSpecies(i))%MassIC * partV2 &
                                            * PartMPF(i)            
        ELSE
          Ekin(nSpecies+1) = Ekin(nSpecies+1) + 0.5 *  Species(PartSpecies(i))%MassIC * partV2 &
                                            *  Species(PartSpecies(i))%MacroParticleFactor
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + 0.5 *  Species(PartSpecies(i))%MassIC * partV2 &
                                            *  Species(PartSpecies(i))%MacroParticleFactor
        END IF != usevMPF
      ELSE ! partV2 > 1e6
  !       Ekin = Ekin + (gamma - 1) * mass * MPF *c^2
        Gamma = partV2/c2      
        Gamma = 1./SQRT(1.-Gamma)
        IF(usevMPF) THEN
          Ekin(nSpecies+1) = Ekin(nSpecies+1) + PartMPF(i) * (Gamma-1.) &
                        * Species(PartSpecies(i))%MassIC * c2
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + PartMPF(i) * (Gamma-1.) &
                        * Species(PartSpecies(i))%MassIC * c2
        ELSE
          Ekin(nSpecies+1) = Ekin(nSpecies+1) + (Gamma-1.) &
                        * Species(PartSpecies(i))%MassIC &
                        * Species(PartSpecies(i))%MacroParticleFactor * c2
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + (Gamma-1.) &
                        * Species(PartSpecies(i))%MassIC &
                        * Species(PartSpecies(i))%MacroParticleFactor * c2
        END IF !=usevMPF
      END IF ! partV2
    END IF
  END DO
ELSE ! nEkin = 1 : only 1 species
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
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
        Gamma = partV2/c2      
        Gamma = 1./SQRT(1.-Gamma)
        IF(usevMPF) THEN
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + PartMPF(i) * (Gamma-1.) &
                        * Species(PartSpecies(i))%MassIC * c2
        ELSE
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + (Gamma -1.) &
                        * Species(PartSpecies(i))%MassIC &
                        * Species(PartSpecies(i))%MacroParticleFactor * c2
        END IF ! usevMPF
      END IF ! partV2
    END IF ! particle inside
  END DO ! particleveclength
END IF

END SUBROUTINE CalcKineticEnergy



SUBROUTINE CalcTemperature(Temp, NumSpec)
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Equation_Vars,      ONLY : c2 
USE MOD_Particle_Vars,      ONLY : PartState, PartSpecies, Species, PDM, nSpecies, BoltzmannConst, PartMPF, usevMPF
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Temp(:) 
INTEGER, INTENT(OUT)            :: NumSpec(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i, iSpec
REAL              :: Mean_PartV2(nSpecies, 3), MeanPartV_2(nSpecies,3) , TempDirec(nSpecies,3), RealNumSpec(nSpecies) 
REAL              :: PartV(nSpecies, 3), PartV2(nSpecies,3) !PartVx, PartVY, PartVZ, PartVx2, PartVY2, PartVZ2
!===================================================================================================================================
! Compute velocity averages
  PartV = 0
  PartV2 = 0
  NumSpec = 0
  RealNumSpec = 0
! Sum up velocity 
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      IF (usevMPF) THEN 
        PartV(PartSpecies(i),1:3) = PartV(PartSpecies(i),1:3) + PartState(i,4:6) * PartMPF(i)
        PartV2(PartSpecies(i),1:3) = PartV2(PartSpecies(i),1:3) + PartState(i,4:6)**2 * PartMPF(i) 
        RealNumSpec(PartSpecies(i)) = RealNumSpec(PartSpecies(i)) + PartMPF(i)
        NumSpec(PartSpecies(i)) = NumSpec(PartSpecies(i)) + 1
      ELSE
        PartV(PartSpecies(i),1:3) = PartV(PartSpecies(i),1:3) + PartState(i,4:6)
        PartV2(PartSpecies(i),1:3) = PartV2(PartSpecies(i),1:3) + PartState(i,4:6)**2 
        NumSpec(PartSpecies(i)) = NumSpec(PartSpecies(i)) + 1
      END IF
    END IF
  END DO
  DO iSpec=1, nSpecies
    IF(NumSpec(iSpec).NE.0) THEN
      ! Compute velocity averages
      IF (usevMPF) THEN
        MeanPartV_2(iSpec,1:3)  = (PartV(iSpec,1:3) / RealNumSpec(iSpec))**2       ! < |v| >**2
        Mean_PartV2(iSpec,1:3)  = PartV2(iSpec,1:3) / RealNumSpec(iSpec)           ! < |v|**2 >
      ELSE
        MeanPartV_2(iSpec,1:3)  = (PartV(iSpec,1:3) / NumSpec(iSpec))**2       ! < |v| >**2
        Mean_PartV2(iSpec,1:3)  = PartV2(iSpec,1:3) / NumSpec(iSpec)           ! < |v|**2 >
      END IF
    ELSE
      MeanPartV_2(iSpec,1:3) = 0.
      Mean_PartV2(iSpec,1:3) = 0.
    END IF
    ! Compute temperatures
    TempDirec(iSpec,1:3) = Species(iSpec)%MassIC * (Mean_PartV2(iSpec,1:3) - MeanPartV_2(iSpec,1:3)) &
         / BoltzmannConst ! Temp calculation is limitedt to one species
    Temp(iSpec) = (TempDirec(iSpec,1) + TempDirec(iSpec,2) + TempDirec(iSpec,3))/3
  END DO
END SUBROUTINE CalcTemperature



SUBROUTINE CalcIntTempsAndEn(IntTemp, IntEn)
!===================================================================================================================================
! Calculation of internal Temps (TVib, TRot)
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars,      ONLY : PartSpecies, Species, PDM, nSpecies, BoltzmannConst, PartMPF, usevMPF
USE MOD_DSMC_Vars,          ONLY : PartStateIntEn, SpecDSMC, DSMC
USE MOD_DSMC_Analyze,       ONLY : CalcTVib, CalcTelec
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: IntTemp(:,:) , IntEn(:,:) 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i, iSpec
INTEGER           :: NumSpec(nSpecies)
REAL              :: EVib(nSpecies), ERot(nSpecies), Eelec(nSpecies), RealNumSpec(nSpecies)
!REAL              :: CalcTVib
!===================================================================================================================================
NumSpec = 0
EVib = 0
ERot = 0
! set electronic state to zero
IntTemp(:,3) = 0
RealNumSpec  = 0

! Sum up internal energies
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      IF (usevMPF) THEN
        EVib(PartSpecies(i)) = EVib(PartSpecies(i)) + PartStateIntEn(i,1) * PartMPF(i)
        ERot(PartSpecies(i)) = ERot(PartSpecies(i)) + PartStateIntEn(i,2) * PartMPF(i)
        NumSpec(PartSpecies(i)) = NumSpec(PartSpecies(i)) + 1
        RealNumSpec(PartSpecies(i)) = RealNumSpec(PartSpecies(i)) + PartMPF(i)
        IF ( DSMC%ElectronicState ) THEN
          Eelec(PartSpecies(i)) = Eelec(PartSpecies(i)) + PartStateIntEn(i,3) * PartMPF(i)
        END IF
      ELSE
        EVib(PartSpecies(i)) = EVib(PartSpecies(i)) + PartStateIntEn(i,1)
        ERot(PartSpecies(i)) = ERot(PartSpecies(i)) + PartStateIntEn(i,2)
        IF ( DSMC%ElectronicState ) THEN
          Eelec(PartSpecies(i)) = Eelec(PartSpecies(i)) + PartStateIntEn(i,3)
        END IF
       NumSpec(PartSpecies(i)) = NumSpec(PartSpecies(i)) + 1
      END IF
    END IF
  END DO
! Calc TVib, TRot
  DO iSpec = 1, nSpecies
    IF((SpecDSMC(iSpec)%InterID.EQ.2).AND.(NumSpec(iSpec).GT.0)) THEN !NumPart gt 0 and species is molecule
      IF (usevMPF) THEN
        IntTemp(iSpec,2) = ERot(iSpec)/(BoltzmannConst*RealNumSpec(iSpec))  !Calc TRot
      ELSE
        IntTemp(iSpec,2) = ERot(iSpec)/(BoltzmannConst*NumSpec(iSpec))  !Calc TRot
      END IF
      IF (EVib(iSpec).GT.0) THEN
        IF (DSMC%VibEnergyModel.EQ.0) THEN              ! SHO-model
          IF (usevMPF) THEN
            IntTemp(iSpec,1) = SpecDSMC(iSpec)%CharaTVib/LOG(1 + 1/(EVib(iSpec) & 
                            /(RealNumSpec(iSpec)*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)-DSMC%GammaQuant))
          ELSE
            IntTemp(iSpec,1) = SpecDSMC(iSpec)%CharaTVib/LOG(1 + 1/(EVib(iSpec) & 
                            /(NumSpec(iSpec)*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)-DSMC%GammaQuant))
          END IF
        ELSE                                            ! TSHO-model
          IF (usevMPF) THEN
            IntTemp(iSpec,1) = CalcTVib(SpecDSMC(iSpec)%CharaTVib, EVib(iSpec)/RealNumSpec(iSpec), SpecDSMC(iSpec)%MaxVibQuant)
          ELSE
            IntTemp(iSpec,1) = CalcTVib(SpecDSMC(iSpec)%CharaTVib, EVib(iSpec)/NumSpec(iSpec), SpecDSMC(iSpec)%MaxVibQuant)
          END IF
        END IF
      ELSE
        IntTemp(iSpec,1) = 0        
      END IF
    ELSE
      IntTemp(iSpec,1) = 0
      IntTemp(iSpec,2) = 0
    END IF
    IF(usevMPF) THEN
      IF ( DSMC%ElectronicState ) THEN
        IF ( NumSpec(iSpec) .GT. 0 ) THEN
          IntTemp(iSpec,3) = CalcTelec( Eelec(iSpec)/RealNumSpec(iSpec),iSpec )
        END IF
        IntEn(iSpec,3) = Eelec(iSpec)
      END IF
      IntEn(iSpec,1) = EVib(iSpec)
      IntEn(iSpec,2) = ERot(iSpec)   
    ELSE
      IF ( DSMC%ElectronicState ) THEN
        IF ( NumSpec(iSpec) .GT. 0 ) THEN
          IntTemp(iSpec,3) = CalcTelec( Eelec(iSpec)/NumSpec(iSpec),iSpec )
        END IF
        IntEn(iSpec,3) = Eelec(iSpec) * Species(iSpec)%MacroParticleFactor
      END IF
      IntEn(iSpec,1) = EVib(iSpec) * Species(iSpec)%MacroParticleFactor 
      IntEn(iSpec,2) = ERot(iSpec) * Species(iSpec)%MacroParticleFactor 
    END IF
  END DO

END SUBROUTINE CalcIntTempsAndEn



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
USE MOD_DSMC_Vars,          ONLY: CollInf, DSMC, ChemReac
USE MOD_TimeDisc_Vars,      ONLY: dt
USE MOD_Particle_Vars,      ONLY: GEO, Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: RRate(:) 
INTEGER, INTENT(IN)            :: NumSpec(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iReac
!===================================================================================================================================
  DO iReac=1, ChemReac%NumOfReact
    SELECT CASE(ChemReac%ReactType(iReac))
    CASE('R')
        RRate(iReac) = ChemReac%NumReac(iReac) * Species(ChemReac%DefinedReact(iReac,2,1))%MacroParticleFactor &
                     * GEO%Volume(1)**2 / (dt &
                     * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,1)) &
                     * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,2)) &
                     * Species(ChemReac%DefinedReact(iReac,1,3))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,3)) )
    CASE('D')
        RRate(iReac) = ChemReac%NumReac(iReac) * Species(ChemReac%DefinedReact(iReac,2,1))%MacroParticleFactor &
                     * GEO%Volume(1) / (dt &
                     * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,1)) &
                     * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,2)) )
    CASE('E')
        RRate(iReac) = ChemReac%NumReac(iReac) * Species(ChemReac%DefinedReact(iReac,2,1))%MacroParticleFactor &
                     * GEO%Volume(1) / (dt &
                     * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,1)) &
                     * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,2)) )
    CASE('I')
        RRate(iReac) = ChemReac%NumReac(iReac) * Species(ChemReac%DefinedReact(iReac,2,1))%MacroParticleFactor &
                     * GEO%Volume(1) / (dt &
                     * Species(ChemReac%DefinedReact(iReac,1,1))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,1)) &
                     * Species(ChemReac%DefinedReact(iReac,1,2))%MacroParticleFactor * NumSpec(ChemReac%DefinedReact(iReac,1,2)) )
    END SELECT
  END DO
  ChemReac%NumReac = 0
END SUBROUTINE ReacRates

#if ( PP_TimeDiscMethod == 42)
SUBROUTINE ElectronicTransition (  Time, NumSpec )
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,          ONLY: DSMC, SpecDSMC
  USE MOD_TimeDisc_Vars,      ONLY: dt
  USE MOD_Particle_Vars,      ONLY: GEO, nSpecies, Species
  ! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  REAL,INTENT(IN)                :: Time
  INTEGER, INTENT(IN)            :: NumSpec(:)
  INTEGER                        :: iSpec, iSpec2, iunit, iQua1, iQua2, MaxElecQua, ii
  CHARACTER(LEN=128)             :: FileNameTransition
  LOGICAL                        :: bExist
! accary of kf
!===================================================================================================================================

  IF ( DSMC%ElectronicState ) THEN
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
                                      SpecDSMC(iSpec)%ElectronicTransition(iSpec2,iQua1,iQua2) * GEO%Volume(1) &
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
      SpecDSMC(iSpec)%ElectronicTransition = 0
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
  IF ( DSMC%ElectronicState ) THEN
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



END MODULE MOD_Particle_Analyze
