!==================================================================================================================================
! Copyright (c) 2010 - 2019 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

#if defined(PARTICLES) && USE_HDG
MODULE MOD_Part_BR_Elecron_Fluid
!===================================================================================================================================
!> Module for particle insertion via electron density conversion from fluid model to kinetic particles
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: InitSwitchBRElectronModel
PUBLIC :: InitializeVariablesElectronFluidRegions
PUBLIC :: SwitchBRElectronModel
PUBLIC :: CreateElectronsFromBRFluid
PUBLIC :: GetNextBRSwitchTime
!===================================================================================================================================
CONTAINS


SUBROUTINE InitSwitchBRElectronModel()
!----------------------------------------------------------------------------------------------------------------------------------!
! Initialize variables (only once, never during load balance restart) for switching between BR electron fluid model and fully
! kinetic model in HDG simulations
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_HDG_Vars
USE MOD_Globals                ,ONLY: myrank, abort, UNIT_StdOut
USE MOD_ReadInTools            ,ONLY: GETLOGICAL,GETREAL,PrintOption
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

!===================================================================================================================================
! Set possibility of either converting kinetic electrons to BR fluid or vice versa during restart
BRConvertElectronsToFluid     = GETLOGICAL('BRConvertElectronsToFluid') ! default is FALSE
BRConvertFluidToElectrons     = GETLOGICAL('BRConvertFluidToElectrons') ! default is FALSE
BRElectronsRemoved            = .FALSE.
BRConvertElectronsToFluidTime = GETREAL('BRConvertElectronsToFluidTime') ! switch from kinetic to BR electron fluid
BRConvertFluidToElectronsTime = GETREAL('BRConvertFluidToElectronsTime') ! switch from BR electron fluid to kinetic electrons
! ! (create new electron particles)
BRTimeStepMultiplier          = GETREAL('BRTimeStepMultiplier') ! Factor that is multiplied with the ManualTimeStep when using BR model
BRConvertModelRepeatedly      = GETLOGICAL('BRConvertModelRepeatedly') ! Repeat the switch between BR and kinetic multiple times

IF(BRConvertModelRepeatedly.AND.&
  ((BRConvertElectronsToFluidTime.LT.0.).OR.(BRConvertFluidToElectronsTime.LT.0.)))THEN
  IPWRITE(UNIT_StdOut,*) "BRConvertModelRepeatedly      =", BRConvertModelRepeatedly
  IPWRITE(UNIT_StdOut,*) "BRConvertElectronsToFluidTime =", BRConvertElectronsToFluidTime
  IPWRITE(UNIT_StdOut,*) "BRConvertFluidToElectronsTime =", BRConvertFluidToElectronsTime
   CALL abort(__STAMP__,&
     'BRConvertModelRepeatedly=T and either BRConvertElectronsToFluidTime or BRConvertFluidToElectronsTime is not set correctly.')
END IF

IF((BRConvertElectronsToFluid.OR.BRConvertFluidToElectrons).AND.&
  ((BRConvertElectronsToFluidTime.GT.0.).OR.(BRConvertFluidToElectronsTime.GT.0.)))THEN
  CALL abort(__STAMP__,'BR electron model: Use either fixed conversion or times but not both!')
END IF

BRConvertMode = 0 ! Initialize

! Both times are given: Two or more switches
IF((BRConvertElectronsToFluidTime.GE.0.).AND.(BRConvertFluidToElectronsTime.GE.0.))THEN
  IF(BRConvertElectronsToFluidTime.GT.BRConvertFluidToElectronsTime)THEN
    ! Mode=1: BR -> kin -> BR (when BRConvertElectronsToFluidTime > BRConvertFluidToElectronsTime)
    IF(BRConvertModelRepeatedly) THEN
       BRConvertMode = 1
     ELSE
       BRConvertMode = -1
     END IF
  ELSEIF(BRConvertFluidToElectronsTime.GT.BRConvertElectronsToFluidTime)THEN
    ! Mode=2: kin -> BR -> kin (when BRConvertFluidToElectronsTime > BRConvertElectronsToFluidTime)
    IF(BRConvertModelRepeatedly) THEN
      BRConvertMode = 2
    ELSE
      BRConvertMode = -2
    END IF
  ELSE
    CALL abort(__STAMP__,'BRConvertFluidToElectronsTime == BRConvertElectronsToFluidTime is not allowed!')
  END IF ! BRConvertElectronsToFluidTime.GT.BRConvertFluidToElectronsTime
ELSEIF(BRConvertElectronsToFluidTime.GE.0.)THEN
  BRConvertMode = 3 ! Single Switch
ELSEIF(BRConvertFluidToElectronsTime.GE.0.)THEN
  BRConvertMode = 3 ! Single Switch
END IF ! (BRConvertElectronsToFluidTime.GE.0.).AND.(BRConvertFluidToElectronsTime.GE.0.)

CALL PrintOption('Switch BR Electron <-> Kinetic: BRConvertMode (zero means OFF)' , 'INFO' , IntOpt=BRConvertMode)

END SUBROUTINE InitSwitchBRElectronModel


SUBROUTINE InitializeVariablesElectronFluidRegions()
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_HDG_Vars
USE MOD_ReadInTools
USE MOD_Particle_Vars
USE MOD_Restart_Vars        ,ONLY: DoRestart,RestartTime
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(32)         :: hilf, hilf2
INTEGER               :: iRegions
REAL                  :: phimax_tmp
!===================================================================================================================================
!-- Read parameters for region mapping
NbrOfRegions = GETINT('NbrOfRegions','0')
UseBRElectronFluid = .FALSE. ! Initialize
IF (NbrOfRegions .GT. 0) THEN
  UseBRElectronFluid = .TRUE.
  ALLOCATE(RegionBounds(1:6,1:NbrOfRegions))
  DO iRegions=1,NbrOfRegions
    WRITE(UNIT=hilf2,FMT='(I0)') iRegions
    RegionBounds(1:6,iRegions) = GETREALARRAY('RegionBounds'//TRIM(hilf2),6,'0. , 0. , 0. , 0. , 0. , 0.')
  END DO

  CALL MapBRRegionToElem()
  ALLOCATE(RegionElectronRef(1:3,1:NbrOfRegions))
  DO iRegions=1,NbrOfRegions
    WRITE(UNIT=hilf2,FMT='(I0)') iRegions
    ! 1:3 - rho_ref, phi_ref, and Te[eV]
    RegionElectronRef(1:3,iRegions) = GETREALARRAY('Part-RegionElectronRef'//TRIM(hilf2),3,'0. , 0. , 1.')
    WRITE(UNIT=hilf,FMT='(G0)') RegionElectronRef(2,iRegions)
    phimax_tmp = GETREAL('Part-RegionElectronRef'//TRIM(hilf2)//'-PhiMax',TRIM(hilf))
    IF (phimax_tmp.NE.RegionElectronRef(2,iRegions)) THEN !shift reference point (rho_ref, phi_ref) to phi_max:
      RegionElectronRef(1,iRegions) = RegionElectronRef(1,iRegions) &
        * EXP((phimax_tmp-RegionElectronRef(2,iRegions))/RegionElectronRef(3,iRegions))
      RegionElectronRef(2,iRegions) = phimax_tmp
      SWRITE(*,*) 'WARNING: BR-reference point is shifted to:', RegionElectronRef(1:2,iRegions)
    END IF
  END DO
END IF

! Check whether it is a restart or a fresh computation
IF(.NOT.DoRestart)THEN ! When starting at t=0
  
  ! With switch BR -> kin -> BR (both cases: -2 switch twice and +2 switch multiple times)
  IF(ABS(BRConvertMode).EQ.2) UseBRElectronFluid=.FALSE.

  ! Run kinetic simulation and then switch to BR
  IF((BRConvertElectronsToFluidTime.GT.0.).AND.(BRConvertMode.EQ.3)) UseBRElectronFluid=.FALSE.

ELSE ! Restart (Important: also load balance restarts)

  ! --------------
  ! Single switch
  ! --------------
  IF(BRConvertMode.EQ.3)THEN

    ! When restarting a simulation that has already been switched from BR -> kin
    IF( (BRConvertFluidToElectronsTime.GT.0.).AND.&
        (BRConvertFluidToElectronsTime.LT.RestartTime)) UseBRElectronFluid= .FALSE.

    ! When restarting a simulation that has not yet been switched from kin -> BR
    IF( (BRConvertElectronsToFluidTime.GT.0.).AND.&
        (BRConvertElectronsToFluidTime.GE.RestartTime)) UseBRElectronFluid= .FALSE.

  ! --------------
  ! 2 Switches
  ! --------------
  ELSEIF(BRConvertMode.EQ.-1)THEN ! 2 switches: BR -> kin -> BR
    IF(BRConvertFluidToElectronsTime.LT.RestartTime)THEN
      BRConvertFluidToElectronsTime = -1. ! BR -> kin has already happened
      IF(BRConvertElectronsToFluidTime.GE.RestartTime) UseBRElectronFluid= .FALSE. ! not yet switched from kin -> BR
    END IF

  ELSEIF(BRConvertMode.EQ.-2)THEN ! 2 switches: kin -> BR -> kin
    IF(BRConvertElectronsToFluidTime.LT.RestartTime)THEN
      BRConvertElectronsToFluidTime = -1.
      IF(BRConvertFluidToElectronsTime.LT.RestartTime) UseBRElectronFluid = .FALSE. ! BR -> kin has already happened
    ELSE
      UseBRElectronFluid = .FALSE. ! still kinetic: kin -> BR has not happened yet
    END IF

  ! --------------
  ! Multiple Switches
  ! --------------
  ELSEIF(BRConvertMode.EQ.1)THEN ! Multiple switches: BR -> kin -> BR
    ASSOCIATE( t       => MOD(RestartTime,BRConvertElectronsToFluidTime) ,&
               tBR2Kin => BRConvertFluidToElectronsTime                  )
      ! Check if restart falls in kinetic region
      IF(t.GT.0.0)THEN
        ! not at t=0 (restart) or t=BRConvertElectronsToFluidTime
        IF(t.GT.tBR2Kin)THEN
          UseBRElectronFluid = .FALSE. ! not yet switched kin -> BR
          ! fix tolerance issue: restart from HDG-BR file, but comes out as non-BR above -> this is actually false: t.GT.tBR2Kin
          IF(ALMOSTEQUALRELATIVE(t,tBR2Kin,1e-6)) UseBRElectronFluid = .TRUE.
        END IF ! t.GT.tBR2Kin
        ! fix tolerance issue: t = 2.1175823681357508E-022 -> restart from kinetic at t=BRConvertElectronsToFluidTime
        ! catch MOD rest with arbitrary limit 1e-9
        IF(t.LT.1e-9) UseBRElectronFluid = .FALSE. ! restart from kinetic simulation
      ELSE
        IF(RestartTime.GT.0.) UseBRElectronFluid = .FALSE. ! not yet switched kin -> BR
      END IF ! t.GT.0.0
    END ASSOCIATE


  ELSEIF(BRConvertMode.EQ.2)THEN ! Multiple switches: kin -> BR -> kin
    ! Check if restart falls in kinetic region
    IF(MOD(RestartTime,BRConvertFluidToElectronsTime).LE.BRConvertElectronsToFluidTime)THEN
      ! 1.) Restart at 0 must be kinetic
      ! 2.) Restart at exactly BRConvertFluidToElectronsTime must be BR
      IF(RestartTime.GT.0.)THEN
        IF(MOD(RestartTime,BRConvertFluidToElectronsTime).GT.0.) UseBRElectronFluid = .FALSE. ! not yet switched kin -> BR
      ELSE
        UseBRElectronFluid = .FALSE. ! restart from t=0 file
      END IF ! RestartTime.GT.0.
    ELSE
      ! fix tolerance issue: restart from HDG-kinetic file, but comes out as BR because MOD rest is a little bit smaller than 
      ! BRConvertElectronsToFluidTime, but is actually the same value. Example:
      !   MOD(RestartTime,BRConvertFluidToElectronsTime) = 5.000000000000001e-6
      !   BRConvertElectronsToFluidTime                  = 4.999999999999990e-6
      IF(ALMOSTEQUALRELATIVE(MOD(RestartTime,BRConvertFluidToElectronsTime),BRConvertElectronsToFluidTime,1e-5)) UseBRElectronFluid = .FALSE.
    END IF ! MOD(RestartTime,BRConvertFluidToElectronsTime).LE.BRConvertElectronsToFluidTime

  END IF ! BRConvertMode.EQ.3

END IF ! .NOT.DoRestart

! Sanity Check
IF((.NOT.UseBRElectronFluid).AND.BRConvertElectronsToFluid)THEN
  SWRITE(UNIT_StdOut,*) "UseBRElectronFluid        =", UseBRElectronFluid
  SWRITE(UNIT_StdOut,*) "BRConvertElectronsToFluid =", BRConvertElectronsToFluid
  CALL abort(__STAMP__,'UseBRElectronFluid and BRConvertElectronsToFluid MUST be both true. Define BR electron fluid mode!')
END IF ! (.NOT.UseBRElectronFluid).AND.BRConvertElectronsToFluid

! Sanity Check
IF(UseBRElectronFluid.AND.BRConvertFluidToElectrons)THEN
  SWRITE(UNIT_StdOut,*) "UseBRElectronFluid        =", UseBRElectronFluid
  SWRITE(UNIT_StdOut,*) "BRConvertFluidToElectrons =", BRConvertFluidToElectrons
  CALL abort(__STAMP__,'UseBRElectronFluid and BRConvertFluidToElectrons CANNOT both be true. Deactivate BR electron fluid model!')
END IF ! UseBRElectronFluid.AND.BRConvertFluidToElectrons

END SUBROUTINE InitializeVariablesElectronFluidRegions


!===================================================================================================================================
!> For BR Electron / fully kinetic model switch, get the next time a switch is going to be performed
!> Either one/two or multiple switches are possible depending on the user settings
!===================================================================================================================================
SUBROUTINE GetNextBRSwitchTime()
! MODULES
USE MOD_Globals       ,ONLY: abort
USE MOD_TimeDisc_Vars ,ONLY: dt_Min,Time
USE MOD_HDG_Vars      ,ONLY: BRConvertFluidToElectronsTime,BRConvertElectronsToFluidTime,BRConvertMode
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: t             !> time in relative frame
REAL :: tMin,tMax     !> sorted t1 and t2
REAL :: tBRSwitchDiff !> delta to next switch (BR or kinetic)
!===================================================================================================================================
!tBRSwitchDiff = tBRSwitch-Time ! Time to BR<->kinetic switch, use extra variable so number doesn't change due to numerical errors
ASSOCIATE( t1 => BRConvertFluidToElectronsTime,&
           t2 => BRConvertElectronsToFluidTime &
         )
  SELECT CASE (BRConvertMode)
    CASE (3,-1,-2) ! Single switch OR 2 switches: BR -> kin -> BR (t1 < t2) OR kin -> BR -> kin (t1 > t2)
      t = Time
    CASE (1)  ! Multiple switches: BR -> kin -> BR (t1 < t2)
      t = MOD(Time,t2)
    CASE (2)  ! Multiple switches: kin -> BR -> kin (t1 > t2)
      t = MOD(Time,t1)
    CASE DEFAULT
      CALL abort(__STAMP__,'Unknown value for BRConvertMode =',IntInfoOpt=BRConvertMode)
  END SELECT

  ! Calculate delta time
  IF((t1.GT.0.0).AND.(t2.GT.0.0))THEN
    tMin=MIN(t1,t2)
    tMax=MAX(t1,t2)
    IF(t.GT.tMin)THEN
      tBRSwitchDiff = tMax
    ELSE
      tBRSwitchDiff = tMin
    END IF ! t.GT.t1
  ELSEIF(t1.GT.0.0)THEN
    tBRSwitchDiff = t1
  ELSEIF(t2.GT.0.0)THEN
    tBRSwitchDiff = t2
  END IF ! (t1.GT.0.0).AND.(t2.GT.0)

  ! Catch tolerance issue which leads to a timestep of 1.0123123E-23 (basically zero)
  IF(.NOT.ALMOSTEQUALRELATIVE(tBRSwitchDiff,t,1e-5))THEN
    tBRSwitchDiff = tBRSwitchDiff - t
  ELSE
    tBRSwitchDiff = 0.0
  END IF ! .NOT.ALMOSTEQUALRELATIVE(tBRSwitchDiff,t,1e-5)

  ! Set dt_Min(DT_BR_SWITCH)
  IF(tBRSwitchDiff.GT.0.0)THEN
    dt_Min(DT_BR_SWITCH) = tBRSwitchDiff
  ELSE
    dt_Min(DT_BR_SWITCH) = HUGE(1.)
  END IF ! tBRSwitchDiff.GT.0.0

!WRITE (*,*) "dt_Min(DT_BR_SWITCH),t2,t1 =", dt_Min(DT_BR_SWITCH),t2,t1
END ASSOCIATE
!WRITE (*,*) "dt_Min(DT_BR_SWITCH) =", dt_Min(DT_BR_SWITCH)
!read*

END SUBROUTINE GetNextBRSwitchTime


SUBROUTINE SwitchBRElectronModel()
!----------------------------------------------------------------------------------------------------------------------------------!
! Switch between BR electron fluid and kinetic model
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_HDG_Vars
USE MOD_TimeDisc_Vars         ,ONLY: time,iter,dt_Min
USE MOD_Elem_Mat              ,ONLY: Elem_Mat
USE MOD_part_operations       ,ONLY: RemoveAllElectrons
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL        :: debug
!===================================================================================================================================
!debug=.true.
debug=.false.

ASSOCIATE( tBR2Kin => BRConvertFluidToElectronsTime ,&
           tKin2BR => BRConvertElectronsToFluidTime )
  ! BR -> kinetic
  IF(UseBRElectronFluid.AND.tBR2Kin.GT.0.0)THEN
    IF((.NOT.BRConvertModelRepeatedly).AND.(time.GE.tBR2Kin)                                             .OR.&
       ((BRConvertMode.EQ.1)          .AND.(GreaterEqualWithTolerance(MOD(time,tKin2BR),tBR2Kin)))       .OR.&
       ((BRConvertMode.EQ.2)          .AND.(LesserThanWithTolerance(MOD(time,tBR2Kin),tKin2BR,tBR2Kin))) )THEN
     IF(debug)THEN
      write(*,*) ""
      IPWRITE(UNIT_StdOut,*) "SWITCH TO kinetic ?",time
      write(*,*) ""
      read*
     END IF ! debug
      CALL CreateElectronsFromBRFluid(.FALSE.) ! Use BR electron fluid model density to create kinetic electrons in each cell
      CALL Elem_Mat(iter)          ! Recompute elem matrices
      UseBRElectronFluid = .FALSE. ! Deactivate BR fluid
      IF((.NOT.BRConvertModelRepeatedly).AND.(BRConvertMode.EQ.-2)) tKin2BR = -1.0 ! deactivate kin -> BR
    END IF
  ENDIF

  ! kinetic -> BR
  IF(.NOT.UseBRElectronFluid.AND.tKin2BR.GT.0.0)THEN
    IF((.NOT.BRConvertModelRepeatedly).AND.(time.GE.tKin2BR)                                             .OR.&
       ((BRConvertMode.EQ.1)          .AND.(LesserThanWithTolerance(MOD(time,tKin2BR),tBR2Kin,tKin2BR))) .OR.&
       ((BRConvertMode.EQ.2)          .AND.(GreaterEqualWithTolerance(MOD(time,tBR2Kin),tKin2BR)))       )THEN
     IF(debug)THEN
      write(*,*) ""
      IPWRITE(UNIT_StdOut,*) "MOD(time,tKin2BR),tBR2Kin,MOD(time,tKin2BR).LT.tBR2Kin =", MOD(time,tKin2BR),tBR2Kin,MOD(time,tKin2BR).LT.tBR2Kin
      IPWRITE(UNIT_StdOut,*) "SWITCH TO BR ?",time
      write(*,*) ""
      read*
     END IF ! debug
     IF(NbrOfRegions.EQ.0) CALL abort(__STAMP__,'SwitchBRElectronModel(): Cannot switch [kin -> BR] as no BR regions are defined!')
      CALL RemoveAllElectrons()    ! Remove all electron particles from the simulation
      CALL Elem_Mat(iter)          ! Recompute elem matrices
      UseBRElectronFluid = .TRUE.  ! Activate BR fluid
      IF((.NOT.BRConvertModelRepeatedly).AND.(BRConvertMode.EQ.-1))tBR2Kin = -1.0 ! deactivate BR -> kin
    END IF
  END IF ! .NOT.UseBRElectronFluid.AND.BRConvertE.GT.0.0
END ASSOCIATE

! Restore the initial dt_Min from InitTimeStep() as it might have been changed in the previous time step
dt_Min(DT_MIN) = BRTimeStepBackup
! Adjust the time step when BR electron fluid is active. Usually BR electron time step is XX times larger than the fully kinetic
IF(UseBRElectronFluid) dt_Min(DT_MIN) = BRTimeStepMultiplier*dt_Min(DT_MIN)
CALL GetNextBRSwitchTime()
END SUBROUTINE SwitchBRElectronModel


PURE LOGICAL FUNCTION GreaterEqualWithTolerance(a,b)
!===================================================================================================================================
! Check if a >= b or a is almost equal to b via ALMOSTEQUALRELATIVE
! Catch tolerance issues when a is only an epsilon smaller than b but the inquiry should be that they are equal
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: a,b !< Two real numbers for comparison
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,PARAMETER :: tol=1.0e-6 ! fix for tolerance issues
!===================================================================================================================================
IF((a.GE.b).OR.(ALMOSTEQUALRELATIVE(a,b,tol)))THEN
  GreaterEqualWithTolerance = .TRUE.
ELSE
  GreaterEqualWithTolerance = .FALSE.
END IF
END FUNCTION GreaterEqualWithTolerance


PURE LOGICAL FUNCTION LesserThanWithTolerance(a,b,c)
!===================================================================================================================================
! Check if a < b and NOT a is almost equal to b via ALMOSTEQUALRELATIVE
! Catch tolerance issues when a<b returns a false positive, because the numbers are actually the same (with an epsilon difference)
! Example: a = 5.4999999999999995E-006
!          b = 5.4999999999999996E-006
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: a,b,c !< Two real numbers for comparison
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,PARAMETER :: tol=1.0e-6 ! fix for tolerance issues
!===================================================================================================================================
IF((a.LT.b).AND.(.NOT.ALMOSTEQUALRELATIVE(a,b,tol)))THEN
  LesserThanWithTolerance = .TRUE.
ELSE
  ! Catch tolerance issue, when MOD(A,B) should actually be 0, but instead it is B
  IF(ALMOSTEQUALRELATIVE(a,c,tol))THEN
    LesserThanWithTolerance = .TRUE. ! a=MOD(A,B) returns B because e.g. A/B=5
  ELSE
    LesserThanWithTolerance = .FALSE.
  END IF ! 
END IF
END FUNCTION LesserThanWithTolerance


SUBROUTINE CreateElectronsFromBRFluid(CreateFromRestartFile)
!----------------------------------------------------------------------------------------------------------------------------------!
! 1.) reconstruct the electron phase space using the integrated charge density in each cell 
!     a.) from ElectronDensityCell and ElectronTemperatureCell that are read from .h5 state file (only during restart)
!     b.) from BR electron model variables in each cell (only during the simulation)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Globals             ,ONLY: abort,MPIRoot,UNIT_stdOut,IK,MPI_COMM_WORLD
USE MOD_Globals_Vars        ,ONLY: ElementaryCharge,BoltzmannConst
USE MOD_PreProc             ,ONLY: PP_nElems
USE MOD_Particle_Vars       ,ONLY: PDM,PEM,PartState,nSpecies,Species,PartSpecies,usevMPF
USE MOD_PIC_Analyze         ,ONLY: CalculateBRElectronsPerCell
USE MOD_Mesh_Vars           ,ONLY: NGeo,XCL_NGeo,XiCL_NGeo,wBaryCL_NGeo,offsetElem
USE MOD_DSMC_Vars           ,ONLY: CollisMode,DSMC,PartStateIntEn
USE MOD_part_emission_tools ,ONLY: CalcVelocity_maxwell_lpn
USE MOD_DSMC_Vars           ,ONLY: useDSMC
USE MOD_Eval_xyz            ,ONLY: TensorProductInterpolation
USE MOD_HDF5_input          ,ONLY: OpenDataFile,CloseDataFile,ReadArray
USE MOD_HDF5_Input          ,ONLY: File_ID,DatasetExists
USE MOD_Mesh_Vars           ,ONLY: offsetElem
USE MOD_Restart_Vars        ,ONLY: RestartFile
USE MOD_Particle_Mesh_Vars  ,ONLY: ElemVolume_Shared
USE MOD_HDG_Vars            ,ONLY: ElemToBRRegion,RegionElectronRef
USE MOD_Mesh_Tools          ,ONLY: GetCNElemID
USE MOD_Part_Tools          ,ONLY: UpdateNextFreePosition
USE MOD_TimeDisc_Vars       ,ONLY: time
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
LOGICAL,INTENT(IN)             :: CreateFromRestartFile
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL,DIMENSION(1,1:PP_nElems)  :: ElectronDensityCell,ElectronTemperatureCell
REAL,ALLOCATABLE               :: ElectronDensityCell(:,:),ElectronTemperatureCell(:,:)
INTEGER                        :: ElemCharge,ElecSpecIndx,iSpec,iElem,iPart,ParticleIndexNbr,RegionID
REAL                           :: PartPosRef(1:3),ElemTemp
CHARACTER(32)                  :: hilf
LOGICAL                        :: ElectronDensityCellExists,ElectronTemperatureCellExists
REAL                           :: MPF,ElectronNumberCell
INTEGER                        :: BRNbrOfElectronsCreated
!===================================================================================================================================
BRNbrOfElectronsCreated=0
! ---------------------------------------------------------------------------------------------------------------------------------
! 0.) Read the data
! ---------------------------------------------------------------------------------------------------------------------------------
IF(CreateFromRestartFile)THEN
  ALLOCATE(ElectronDensityCell(1,1:PP_nElems))
  ALLOCATE(ElectronTemperatureCell(1,1:PP_nElems))
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
  CALL DatasetExists(File_ID,'ElectronDensityCell',ElectronDensityCellExists)
  IF(ElectronDensityCellExists)THEN
    ! Associate construct for integer KIND=8 possibility
    ASSOCIATE (&
          PP_nElems   => INT(PP_nElems,IK)   ,&
          offsetElem  => INT(offsetElem,IK)   )
      CALL ReadArray('ElectronDensityCell',2,(/1_IK,PP_nElems/),offsetElem,2,RealArray=ElectronDensityCell)
    END ASSOCIATE
  ELSE
    CALL abort(__STAMP__,'ElectronDensityCell container not found in state file.'//&
        ' This is required for CreateElectronsFromBRFluid(). Set BRConvertFluidToElectrons = F to continue the kinetic simulation')
  END IF

  CALL DatasetExists(File_ID,'ElectronTemperatureCell',ElectronTemperatureCellExists)
  IF(ElectronTemperatureCellExists)THEN
    ! Associate construct for integer KIND=8 possibility
    ASSOCIATE (&
          PP_nElems   => INT(PP_nElems,IK)   ,&
          offsetElem  => INT(offsetElem,IK)   )
      CALL ReadArray('ElectronTemperatureCell',2,(/1_IK,PP_nElems/),offsetElem,2,RealArray=ElectronTemperatureCell)
    END ASSOCIATE
  ELSE
    CALL abort(__STAMP__,'ElectronTemperatureCell container not found in state file.'//&
        ' This is required for CreateElectronsFromBRFluid(). Set BRConvertFluidToElectrons = F to continue the kinetic simulation')
  END IF
  CALL CloseDataFile()
END IF ! CreateFromRestartFile

! ---------------------------------------------------------------------------------------------------------------------------------
! 1.) reconstruct electrons
! ---------------------------------------------------------------------------------------------------------------------------------
SWRITE(UNIT_stdOut,'(A,ES25.14E3,A)')' CreateElectronsFromBRFluid(): Reconstructing electrons at t=',time,' from BR electron fluid density in each cell'

! Loop over all species and find the index corresponding to the electron species: take the first electron species that is
! encountered
ElecSpecIndx = -1 ! Initialize the species index for the electron species with -1
DO iSpec = 1, nSpecies
  IF (Species(iSpec)%ChargeIC.GE.0.0) CYCLE
    IF(NINT(Species(iSpec)%ChargeIC/(-ElementaryCharge)).EQ.1)THEN
      ElecSpecIndx = iSpec
    EXIT
  END IF
END DO
IF (ElecSpecIndx.LE.0) CALL abort(&
  __STAMP__&
  ,'Electron species not found. Cannot create electrons without the defined species!')

WRITE(UNIT=hilf,FMT='(I0)') iSpec
SWRITE(UNIT_stdOut,'(A)')'  Using iSpec='//TRIM(hilf)//' as electron species index from BR fluid conversion.'

IF (usevMPF) THEN
  CALL abort(__STAMP__,'vMPF not implemented yet in CreateElectronsFromBRFluid().')
ELSE
  MPF = Species(ElecSpecIndx)%MacroParticleFactor
END IF

! Loop over all elements
DO iElem=1,PP_nElems

  ! Set electron charge number for each cell
  IF(CreateFromRestartFile)THEN
    ElemCharge=NINT(ElectronDensityCell(1,iElem)*ElemVolume_Shared(GetCNElemID(iElem+offSetElem))/MPF)
  ELSE
    RegionID=ElemToBRRegion(iElem)
    CALL CalculateBRElectronsPerCell(iElem,RegionID,ElectronNumberCell)
    ElemCharge=NINT(ElectronNumberCell/MPF)
  END IF ! CreateFromRestartFile

  ! Create electrons, 1 electron for each charge of each element
  DO iPart=1,ElemCharge
    BRNbrOfElectronsCreated = BRNbrOfElectronsCreated + 1

    ! Set the next free position in the particle vector list
    PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + 1
    ParticleIndexNbr            = PDM%nextFreePosition(PDM%CurrentNextFreePosition)
    IF (ParticleIndexNbr.EQ.0) THEN
      CALL Abort(&
      __STAMP__&
      ,'ERROR in CreateElectronsFromBRFluid(): New Particle Number greater max Part Num!')
    END IF
    PDM%ParticleVecLength       = PDM%ParticleVecLength + 1

    !Set new SpeciesID of new particle (electron)
    PDM%ParticleInside(ParticleIndexNbr) = .true.
    PartSpecies(ParticleIndexNbr) = ElecSpecIndx

    ! Place the electron randomly in the reference cell
    CALL RANDOM_NUMBER(PartPosRef(1:3)) ! get random reference space
    PartPosRef(1:3)=PartPosRef(1:3)*2. - 1. ! map (0,1) -> (-1,1)

    ! Get the physical coordinates that correspond to the reference coordinates
    CALL TensorProductInterpolation(PartPosRef(1:3),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,iElem) &
                      ,PartState(1:3,ParticleIndexNbr)) !Map into phys. space

    ! Set the internal energies (vb, rot and electronic) to zero if needed
    IF ((useDSMC).AND.(CollisMode.GT.1)) THEN
      PartStateIntEn(1,ParticleIndexNbr) = 0.
      PartStateIntEn(2,ParticleIndexNbr) = 0.
      IF (DSMC%ElectronicModel.GT.0)  PartStateIntEn(3,ParticleIndexNbr) = 0.
    END IF

    ! Set the element ID of the electron to the current element ID
    PEM%GlobalElemID(ParticleIndexNbr) = iElem + offsetElem

    ! Set the electron velocity using the Maxwellian distribution (use the function that is suitable for small numbers)
    IF(CreateFromRestartFile)THEN
      ElemTemp = ElectronTemperatureCell(1,iElem)
    ELSE
      ElemTemp = RegionElectronRef(3,RegionID)*ElementaryCharge/BoltzmannConst ! convert eV to K
    END IF ! CreateFromRestartFile
    CALL CalcVelocity_maxwell_lpn(ElecSpecIndx, PartState(4:6,ParticleIndexNbr),Temperature=ElemTemp)
  END DO
END DO

IF(CreateFromRestartFile)THEN
  DEALLOCATE(ElectronDensityCell)
  DEALLOCATE(ElectronTemperatureCell)
END IF ! CreateFromRestartFile

#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,BRNbrOfElectronsCreated,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,iError)
#endif /*USE_MPI*/

SWRITE(UNIT_StdOut,'(A,I0,A)') '  Created a total of ',BRNbrOfElectronsCreated,' electrons.'
!read*

! Update
CALL UpdateNextFreePosition()


END SUBROUTINE CreateElectronsFromBRFluid


SUBROUTINE MapBRRegionToElem()
!----------------------------------------------------------------------------------------------------------------------------------!
! map a particle region to element
! check only element barycenter, nothing else
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_HDG_Vars           ,ONLY: NbrOfRegions, RegionBounds,ElemToBRRegion
USE MOD_Mesh_Vars          ,ONLY: ElemBaryNGeo
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
 INTEGER                :: iElem, iRegions
!===================================================================================================================================
SDEALLOCATE(ElemToBRRegion)
ALLOCATE(ElemToBRRegion(1:PP_nElems))
ElemToBRRegion=0

DO iElem=1,PP_nElems
  DO iRegions=1,NbrOfRegions
    IF ((ElemBaryNGeo(1,iElem).LT.RegionBounds(1,iRegions)).OR.(ElemBaryNGEO(1,iElem).GE.RegionBounds(2,iRegions))) CYCLE
    IF ((ElemBaryNGeo(2,iElem).LT.RegionBounds(3,iRegions)).OR.(ElemBaryNGEO(2,iElem).GE.RegionBounds(4,iRegions))) CYCLE
    IF ((ElemBaryNGeo(3,iElem).LT.RegionBounds(5,iRegions)).OR.(ElemBaryNGEO(3,iElem).GE.RegionBounds(6,iRegions))) CYCLE
    IF (ElemToBRRegion(iElem).EQ.0) THEN
      ElemToBRRegion(iElem)=iRegions
    ELSE
      CALL ABORT(__STAMP__,'Defined regions are overlapping')
    END IF
  END DO ! iRegions=1,NbrOfRegions
END DO ! iElem=1,PP_nElems
END SUBROUTINE MapBRRegionToElem


END MODULE MOD_Part_BR_Elecron_Fluid
#endif /*defined(PARTICLES) && USE_HDG*/
