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

!===================================================================================================================================
!> Module contains the routines for load balancing
!===================================================================================================================================
MODULE MOD_LoadBalance
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
#if USE_MPI
INTERFACE InitLoadBalance
  MODULE PROCEDURE InitLoadBalance
END INTERFACE

INTERFACE FinalizeLoadBalance
  MODULE PROCEDURE FinalizeLoadBalance
END INTERFACE

#if USE_LOADBALANCE
INTERFACE ComputeElemLoad
  MODULE PROCEDURE ComputeElemLoad
END INTERFACE
#endif /*USE_LOADBALANCE*/

INTERFACE LoadBalance
  MODULE PROCEDURE LoadBalance
END INTERFACE

PUBLIC::InitLoadBalance,FinalizeLoadBalance,LoadBalance
#if USE_LOADBALANCE
PUBLIC::ComputeElemLoad
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/

PUBLIC::DefineParametersLoadBalance
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersLoadBalance()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("LoadBalance")

#if USE_LOADBALANCE
CALL prms%CreateLogicalOption( 'DoLoadBalance'                ,  "Set flag for doing dynamic LoadBalance.", '.FALSE.')
CALL prms%CreateIntOption(     'LoadBalanceSample'            ,  "Define number of iterations (before analyze_dt)"//&
  " that are used for calculation of elemtime information"    , value='1')
CALL prms%CreateLogicalOption( 'PartWeightLoadBalance'        ,  "Set flag for doing LoadBalance with partMPIWeight instead of"//&
  " elemtimes. Elemtime array in state file is filled with nParts*PartMPIWeight for each Elem. "//&
  " If Flag [TRUE] LoadBalanceSample is set to 0 and vice versa.", '.FALSE.')
CALL prms%CreateRealOption(    'Load-DeviationThreshold'       ,  "Define threshold for dynamic load-balancing.\n"//&
  "Restart performed if (Maxweight-Targetweight)/Targetweight > defined value." , value='0.10')
#endif /*USE_LOADBALANCE*/
CALL prms%CreateRealOption(    'Particles-MPIWeight'          ,  "Define weight of particles for elem loads.\n"//&
  "(only used if ElemTime does not exist or DoLoadBalance=F).", value='0.02')
CALL prms%CreateIntOption(     'WeightDistributionMethod'     ,  "Method for distributing the elem to procs.\n"//&
  "DEFAULT: 1 if Elemtime exits else -1\n"//&
  "-1: elements are equally distributed\n"//&
  " 0: distribute to procs using elemloads\n"//&
  " 1: distribute to procs using elemloads, last proc recieves least\n"//&
  " 2: NOT WORKING\n"//&
  " 3: TODO DEFINE\n"//&
  " 4: TODO DEFINE\n"//&
  " 5/6: iterative smoothing of loads towards last proc\n")

END SUBROUTINE DefineParametersLoadBalance

#if USE_MPI
SUBROUTINE InitLoadBalance()
!===================================================================================================================================
! init load balancing, new initialization of variables for load balancing
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_LoadBalance_Vars
USE MOD_ReadInTools      ,ONLY: GETLOGICAL, GETREAL, GETINT
USE MOD_ReadInTools      ,ONLY: PrintOption
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT LOAD BALANCE ...'

#ifdef PARTICLES
! Read particle MPI weight in order to determine the ElemTime when no time measurement is performed
! Must be read in init (only once) and before the first load balance is determined because if no ElemTimes are used they are
! calculated with ParticleMPIWeight
ParticleMPIWeight = GETREAL('Particles-MPIWeight','0.02')
IF (ParticleMPIWeight.LT.0.0) THEN
  CALL abort(&
      __STAMP__&
      ,' ERROR: Particle weight cannot be negative!')
END IF
#endif /*PARTICLES*/

#if USE_LOADBALANCE
IF(nProcessors.EQ.1)THEN
  DoLoadBalance=.FALSE. ! deactivate loadbalance for single computations
  CALL PrintOption('No LoadBalance (nProcessors=1): DoLoadBalance','INFO',LogOpt=DoLoadBalance)
  DeviationThreshold=HUGE(1.0)
ELSE
  DoLoadBalance = GETLOGICAL('DoLoadBalance','F')
  DeviationThreshold  = GETREAL('Load-DeviationThreshold','0.10')
END IF
LoadBalanceSample   = GETINT('LoadBalanceSample')
PerformPartWeightLB = GETLOGICAL('PartWeightLoadBalance','F')
IF (PerformPartWeightLB) THEN
  LoadBalanceSample = 0 ! deactivate loadbalance sampling of elemtimes if balancing with partweight is enabled
  CALL PrintOption('PartWeightLoadBalance = T : LoadBalanceSample','INFO',IntOpt=LoadBalanceSample)
ELSE IF (LoadBalanceSample.EQ.0) THEN
  PerformPartWeightLB = .TRUE. ! loadbalance (elemtimes) is done with partmpiweight if loadbalancesampling is set to zero
  CALL PrintOption('LoadbalanceSample = 0 : PartWeightLoadBalance','INFO',LogOpt=PerformPartWeightLB)
END IF
#else
DoLoadBalance=.FALSE. ! deactivate loadbalance if no preproc flag is set
DeviationThreshold=HUGE(1.0)
LoadBalanceSample  = 0
PerformPartWeightLB = .FALSE.
#endif /*USE_LOADBALANCE*/
nLoadBalance = 0
nLoadBalanceSteps = 0
PerformLBSample = .FALSE.

#if USE_LOADBALANCE
ALLOCATE( tCurrent(1:LB_NTIMES) )
! Allocation length (1:number of loadbalance times)
! look into piclas.h for more info about time names
tCurrent=0.
#endif /*USE_LOADBALANCE*/

InitLoadBalanceIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT LOAD BALANCE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitLoadBalance


#if USE_LOADBALANCE
SUBROUTINE ComputeElemLoad()
!----------------------------------------------------------------------------------------------------------------------------------!
! compute the element load
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_TimeDisc_Vars          ,ONLY: time
USE MOD_LoadBalance_Vars       ,ONLY: ElemTime,nLoadBalance,tCurrent
#ifndef PP_HDG
USE MOD_PML_Vars               ,ONLY: DoPML,nPMLElems,ElemToPML
#endif /*PP_HDG*/
USE MOD_LoadBalance_Vars       ,ONLY: DeviationThreshold, PerformLoadBalance, LoadBalanceSample
#ifdef PARTICLES
USE MOD_LoadBalance_Vars       ,ONLY: nPartsPerElem,nDeposPerElem,nTracksPerElem
USE MOD_LoadBalance_Vars       ,ONLY: nSurfacefluxPerElem,nPartsPerBCElem,nSurfacePartsPerElem
USE MOD_Particle_Tracking_vars ,ONLY: DoRefMapping
USE MOD_PICDepo_Vars           ,ONLY: DepositionType
USE MOD_LoadBalance_Vars       ,ONLY: nPartsPerElem
USE MOD_LoadBalance_Vars       ,ONLY: ParticleMPIWeight
#endif /*PARTICLES*/
USE MOD_LoadDistribution       ,ONLY: WriteElemTimeStatistics
USE MOD_LoadBalance_Vars       ,ONLY: CurrentImbalance,PerformLBSample
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iElem
#ifdef PARTICLES
INTEGER(KIND=8)       :: HelpSum
REAL                  :: stotalDepos,stotalParts,sTotalTracks,stotalSurfacefluxes,sTotalBCParts,sTotalSurfaceParts
#endif /*PARTICLES*/
!====================================================================================================================================

! If elem times are calculated by time measurement (PerformLBSample) and no Partweight Loadbalance is enabled
IF(PerformLBSample .AND. LoadBalanceSample.GT.0) THEN

  ! number of load balance calls to Compute Elem Load
  nLoadBalance=nLoadBalance+1

#ifdef PARTICLES
  ! ----------------------------------------------
  ! calculate weightings
  stotalDepos=1.0
  sTotalTracks=1.0
  stotalSurfacefluxes=1.0
  stotalBCParts=1.0
  stotalSurfaceParts=1.0
  ! calculate and weight particle number per element
  helpSum=SUM(nPartsPerElem)
  IF(helpSum.GT.0) THEN
    stotalParts=1.0/REAL(helpSum)
  ELSE
    stotalParts=1.0/REAL(PP_nElems)
    nPartsPerElem=1
  END IF
  ! set and weight tracks per element
  IF (DoRefMapping) THEN
    helpSum=SUM(nTracksPerElem)
    IF(SUM(nTracksPerElem).GT.0) THEN
      sTotalTracks=1.0/REAL(helpSum)
    END IF
  END IF
  ! set and weight depositions per element
  helpSum=SUM(nDeposPerElem)
  IF(helpSum.GT.0) THEN
    stotalDepos=1.0/REAL(helpSum)
  END IF
  ! calculate and weight number of surface fluxes (during setvelo,setinterenergies...)
  helpSum=SUM(nSurfacefluxPerElem)
  IF(helpSum.GT.0) THEN
    stotalSurfacefluxes=1.0/REAL(helpSum)
  END IF
  ! calculate weight number of particles used in sampling of bc element for adaptive particle bcs
  helpSum=SUM(nPartsPerBCElem)
  IF(helpSum.GT.0) THEN
    stotalBCParts=1.0/REAL(helpSum)
  ELSE
    stotalBCParts=1.0/REAL(PP_nElems)
    nPartsPerBCElem=1
  END IF
  ! calculate and weight number of surfaces (catalytic, liquid) used for evaporation and surface chemistry
  helpSum=SUM(nSurfacePartsPerElem)
  IF(helpSum.GT.0) THEN
    stotalSurfaceParts=1.0/REAL(helpSum)
  END IF
  ! ----------------------------------------------
#endif /*PARTICLES*/

  ! distribute times of different routines on elements with respective weightings
  DO iElem=1,PP_nElems
    ! time used in dg routines
#ifdef HDG
    ElemTime(iElem) = ElemTime(iElem) + tCurrent(LB_DG)*ElemHDGSides(iElem)/TotalHDGSides + tCurrent(LB_DGANALYZE)/REAL(PP_nElems)
#else
    ElemTime(iElem) = ElemTime(iElem) + (tCurrent(LB_DG) + tCurrent(LB_DGANALYZE))/REAL(PP_nElems)
#endif /*HDG*/
#ifndef PP_HDG
    ! time used in pml routines
    IF(DoPML)THEN
      IF(ElemToPML(iElem).GT.0 ) ElemTime(iElem) = ElemTime(iElem) + tCurrent(LB_PML)/REAL(nPMLElems)
    END IF
#endif /*PP_HDG*/

#ifdef PARTICLES
    ! add particle LB times to elements with respective weightings
    ElemTime(iElem) = ElemTime(iElem)                                             &
        + tCurrent(LB_ADAPTIVE)      * nPartsPerBCElem(iElem)*sTotalBCParts       &
        + tCurrent(LB_INTERPOLATION)   * nPartsPerElem(iElem)*sTotalParts         &
        + tCurrent(LB_PUSH)            * nPartsPerElem(iElem)*sTotalParts         &
        + tCurrent(LB_UNFP)            * nPartsPerElem(iElem)*sTotalParts         &
        + tCurrent(LB_DSMC)            * nPartsPerElem(iElem)*sTotalParts         &
        + tCurrent(LB_CARTMESHDEPO)    * nPartsPerElem(iElem)*sTotalParts         &
        + tCurrent(LB_TRACK)          * nTracksPerElem(iElem)*sTotalTracks        &
        + tCurrent(LB_SURFFLUX)  * nSurfacefluxPerElem(iElem)*stotalSurfacefluxes &
        + tCurrent(LB_SURF)     * nSurfacePartsPerElem(iElem)*stotalSurfaceParts
    ! e.g. 'shape_function', 'shape_function_1d', 'shape_function_cylindrical'
    IF(TRIM(DepositionType(1:MIN(14,LEN(TRIM(ADJUSTL(DepositionType)))))).EQ.'shape_function')THEN
      ElemTime(iElem) = ElemTime(iElem)                              &
          + tCurrent(LB_DEPOSITION) * nDeposPerElem(iElem)*stotalDepos
    END IF
#endif /*PARTICLES*/
  END DO ! iElem=1,PP_nElems
#ifdef PARTICLES
! If no Elem times are calculated but Partweight Loadbalance is enabled
ELSE IF(PerformLBSample .AND. LoadBalanceSample.EQ.0) THEN
  ! number of load balance calls to Compute Elem Load
  nLoadBalance=nLoadBalance+1
  ! no time measurement and particles are present: simply add the ParticleMPIWeight times the number of particles present
  DO iElem=1,PP_nElems
    ElemTime(iElem) = ElemTime(iElem) &
        + nPartsPerElem(iElem)*ParticleMPIWeight + 1.0
  END DO ! iElem=1,PP_nElems
  IF((MAXVAL(nPartsPerElem).GT.0).AND.(MAXVAL(ElemTime).LE.1.0))THEN
    IPWRITE (*,*) "parts, time =", MAXVAL(nPartsPerElem),MAXVAL(ElemTime)
    CALL abort(&
        __STAMP__&
        ,' ERROR: MAXVAL(nPartsPerElem).GT.0 but MAXVAL(ElemTime).LE.1.0 with ParticleMPIWeight=',RealInfoOpt=ParticleMPIWeight)
  END IF
#endif /*PARTICLES*/
END IF

! Determine sum of balance and calculate target balanced weight and communicate via MPI_ALLREDUCE
CALL ComputeImbalance()

! Fill .csv file for performance analysis and load balance: write data line
CALL WriteElemTimeStatistics(WriteHeader=.FALSE.,time=time)

! only check if imbalance is > a given threshold
PerformLoadBalance=.FALSE.
IF(CurrentImbalance.GT.DeviationThreshold) PerformLoadBalance=.TRUE.

#ifdef PARTICLES
nTracksPerElem       = 0
nDeposPerElem        = 0
nPartsPerElem        = 0
nSurfacefluxPerElem  = 0
nPartsPerBCElem      = 0
nSurfacePartsPerElem = 0
#endif /*PARTICLES*/
tCurrent = 0.

END SUBROUTINE ComputeElemLoad
#endif /*USE_LOADBALANCE*/


SUBROUTINE LoadBalance()
!===================================================================================================================================
! routine perfoming the load balancing stuff
!===================================================================================================================================
! USED MODULES
USE MOD_Globals
USE MOD_Globals_vars     ,ONLY: InitializationWallTime
USE MOD_Preproc
USE MOD_Restart          ,ONLY: Restart
USE MOD_Piclas_Init  ,ONLY: InitPiclas,FinalizePiclas
USE MOD_LoadBalance_Vars ,ONLY: ElemTime,nLoadBalanceSteps,NewImbalance,MinWeight,MaxWeight
#ifdef PARTICLES
USE MOD_PICDepo_Vars     ,ONLY: DepositionType
USE MOD_Particle_MPI     ,ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#endif /*PARTICLES*/
USE MOD_LoadBalance_Vars ,ONLY: CurrentImbalance, MaxWeight, MinWeight
USE MOD_LoadBalance_Vars ,ONLY: Currentimbalance, PerformLoadBalance
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!REAL,INTENT(IN)      :: Currentimbalance
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: LB_Time,LB_StartTime
!===================================================================================================================================
! only do load-balance if necessary
IF(.NOT.PerformLoadBalance) THEN
  ElemTime=0.
  InitializationWallTime=0.
  RETURN
END IF

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' PERFORMING LOAD BALANCE ...'
! Measure init duration
LB_StartTime=PICLASTIME()

nLoadBalanceSteps=nLoadBalanceSteps+1
! finialize all arrays
CALL FinalizePiclas(IsLoadBalance=.TRUE.)
! reallocate
CALL InitPiclas(IsLoadBalance=.TRUE.) ! determines new imbalance in InitMesh() -> ReadMesh()

! restart
CALL Restart()

! compute imbalance
!CALL ComputeImbalance(NewImbalance,MaxWeight,MinWeight,ElemTime)  ! --> new imbalance has already been calculated in mesh_readin

! zero ElemTime, the measurement starts again
ElemTime=0.

IF( NewImbalance.GT.CurrentImbalance ) THEN
  SWRITE(UNIT_stdOut,'(A)') ' WARNING: LoadBalance not successful!'
ELSE
  SWRITE(UNIT_stdOut,'(A)') ' LoadBalance successful!'
END IF
SWRITE(UNIT_stdOut,'(A25,ES15.7)') ' OldImbalance: ', CurrentImbalance
SWRITE(UNIT_stdOut,'(A25,ES15.7)') ' NewImbalance: ', NewImbalance
SWRITE(UNIT_stdOut,'(A25,ES15.7)') ' MaxWeight:    ', MaxWeight
SWRITE(UNIT_stdOut,'(A25,ES15.7)') ' MinWeight:    ', MinWeight

#ifdef PARTICLES
! e.g. 'shape_function', 'shape_function_1d', 'shape_function_cylindrical'
IF(TRIM(DepositionType(1:MIN(14,LEN(TRIM(ADJUSTL(DepositionType)))))).EQ.'shape_function')THEN
  ! open receive buffer for number of particles
  CALL IRecvNbofParticles()
  ! send number of particles
  CALL SendNbOfParticles()
  ! finish communication of number of particles and send particles
  CALL MPIParticleSend()
  ! finish communication
  CALL MPIParticleRecv()
END IF
#endif /*PARTICLES*/

! Measure init duration
LB_Time=PICLASTIME()
InitializationWallTime=LB_Time-LB_StartTime
SWRITE(UNIT_stdOut,'(A,F14.2,A)') ' INITIALIZATION DONE! [',InitializationWallTime,' sec ]'
SWRITE(UNIT_stdOut,'(A)')' LOAD BALANCE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE LoadBalance


SUBROUTINE ComputeImbalance()
!----------------------------------------------------------------------------------------------------------------------------------!
! subroutine to compute the imbalance
! Maxweight:        maximum weight of all processes
! Minweight:        minimum weight of all processes
! targetweight:     current target weight
! CurrentImbalance: rel. deviation of current imbalance
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_LoadBalance_Vars,    ONLY:WeightSum, TargetWeight,CurrentImbalance, MaxWeight, MinWeight
USE MOD_LoadBalance_Vars,    ONLY:ElemTime, PerformLBSample, PerformPartWeightLB
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

IF(.NOT.PerformLBSample .AND. .NOT.PerformPartWeightLB) THEN
  WeightSum        = 0.
  TargetWeight     = 0.
  CurrentImbalance = -1.0
ELSE
  WeightSum=SUM(ElemTime)

  IF(ALMOSTZERO(WeightSum))THEN
    IPWRITE(*,*) 'Info: The measured time of all elems is zero. ALMOSTZERO(WeightSum)=.TRUE., WeightSum=',WeightSum
  END IF

  CALL MPI_ALLREDUCE(WeightSum,TargetWeight,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,iError)
  CALL MPI_ALLREDUCE(WeightSum,MaxWeight   ,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,iError)
  CALL MPI_ALLREDUCE(WeightSum,MinWeight   ,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,iError)

  WeightSum    = TargetWeight ! Set total weight for writing to file
  TargetWeight = TargetWeight/nProcessors ! Calculate the average value that is supposed to be the optimally distributed weight

  !IF(ALMOSTZERO(WeightSum))CALL abort( __STAMP__&
      !,' ERROR: after ALLREDUCE, weight sum cannot be zero! WeightSum=',RealInfoOpt=WeightSum)

  ! new computation of current imbalance
  IF(ABS(TargetWeight).EQ.0.)THEN
    CurrentImbalance = 0.
  ELSE IF(ABS(TargetWeight).LT.0.0)THEN
    SWRITE(UNIT_stdOut,'(A,F14.2,A1)')&
        ' ERROR: after ALLREDUCE, WeightSum/TargetWeight cannot be zero! TargetWeight=[',TargetWeight,']'
    CurrentImbalance = HUGE(1.0)
  ELSE
    CurrentImbalance =  (MaxWeight-TargetWeight ) / TargetWeight
  END IF
  SWRITE(UNIT_stdOut,'(A25,ES15.7)') ' MaxWeight:        ', MaxWeight
  SWRITE(UNIT_stdOut,'(A25,ES15.7)') ' MinWeight:        ', MinWeight
  SWRITE(UNIT_stdOut,'(A25,ES15.7)') ' TargetWeight:     ', TargetWeight
  SWRITE(UNIT_stdOut,'(A25,ES15.7)') ' CurrentImbalance: ', CurrentImbalance

END IF

END SUBROUTINE ComputeImbalance


SUBROUTINE FinalizeLoadBalance()
!===================================================================================================================================
! Deallocate arrays
!===================================================================================================================================
! MODULES
USE MOD_LoadBalance_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE( tCurrent  )
InitLoadBalanceIsDone = .FALSE.

END SUBROUTINE FinalizeLoadBalance
#endif /*USE_MPI*/

END MODULE MOD_LoadBalance
