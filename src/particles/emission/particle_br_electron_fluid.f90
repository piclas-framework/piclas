!==================================================================================================================================
! Copyright (c) 2010 - 2019 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_Part_BR_Elecron_Fluid
#if defined(PARTICLES) && USE_HDG
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
PUBLIC :: DefineParametersBR
PUBLIC :: InitSwitchBRElectronModel
PUBLIC :: InitializeVariablesElectronFluidRegions
PUBLIC :: SwitchBRElectronModel
PUBLIC :: CreateElectronsFromBRFluid
PUBLIC :: GetNextBRSwitchTime
PUBLIC :: UpdateVariableRefElectronTemp
PUBLIC :: UpdateNonlinVolumeFac
!===================================================================================================================================
CONTAINS

!===================================================================================================================================
!> description
!===================================================================================================================================
SUBROUTINE DefineParametersBR()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL prms%CreateIntOption(      'BRNbrOfRegions'                   , 'Number of regions to be mapped to Elements', '0')
CALL prms%CreateStringOption(   'BRVariableElectronTemp'           , 'Variable electron reference temperature when using Boltzmann relation electron model (default is using a constant temperature)','constant')
CALL prms%CreateRealArrayOption('BRRegionBounds[$]'                , 'BRRegionBounds ((xmin,xmax,ymin,...)|1:BRNbrOfRegions)', '0. , 0. , 0. , 0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-RegionElectronRef[$]'        , 'rho_ref, phi_ref, and Te[eV] for Region#', '0. , 0. , 1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-RegionElectronRef[$]-PhiMax' , 'max. expected phi for Region#\n (linear approx. above! def.: phi_ref)', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'BRConvertElectronsToFluid'        , 'Remove all electrons when using BR electron fluid', '.FALSE.')
CALL prms%CreateLogicalOption(  'BRConvertFluidToElectrons'        , 'Create electrons from BR electron fluid (requires ElectronDensityCell ElectronTemperatureCell from .h5 state file)', '.FALSE.')
CALL prms%CreateRealOption(     'BRConvertFluidToElectronsTime'    , "Time when BR fluid electrons are to be converted to kinetic particles", '-1.0')
CALL prms%CreateRealOption(     'BRConvertElectronsToFluidTime'    , "Time when kinetic electrons should be converted to BR fluid electrons", '-1.0')
CALL prms%CreateLogicalOption(  'BRConvertModelRepeatedly'         , 'Repeat the switch between BR and kinetic multiple times', '.FALSE.')
CALL prms%CreateRealOption(     'BRTimeStepMultiplier'             , "Factor that is multiplied with the ManualTimeStep when using BR model", '1.0')
CALL prms%CreateLogicalOption(  'BRAutomaticElectronRef'           , 'Automatically obtain the reference parameters (from a fully kinetic simulation)', '.FALSE.')
END SUBROUTINE DefineParametersBR


!===================================================================================================================================
!> Initialize variables (only once, never during load balance restart) for switching between BR electron fluid model and fully
!> kinetic model in HDG simulations
!===================================================================================================================================
SUBROUTINE InitSwitchBRElectronModel()
! MODULES                                                                                                                          !
USE MOD_HDG_Vars
USE MOD_Globals     ,ONLY: abort, UNIT_StdOut
USE MOD_ReadInTools ,ONLY: GETLOGICAL,GETREAL,PrintOption
#if USE_MPI
USE MOD_Globals     ,ONLY: myrank
#endif /*USE_MPI*/
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
  CALL abort(__STAMP__,'BR electron model: Use either a) fixed conversion or b) times, but not both!')
END IF

BRConvertMode = 0 ! Initialize
DeltaTimeBRWindow = -1. ! Initialize

! Both times are given: Two or more switches
IF((BRConvertElectronsToFluidTime.GE.0.).AND.(BRConvertFluidToElectronsTime.GE.0.))THEN
  IF(BRConvertElectronsToFluidTime.GT.BRConvertFluidToElectronsTime)THEN
    ! Mode=1: BR -> kin -> BR (when BRConvertElectronsToFluidTime > BRConvertFluidToElectronsTime)
    IF(BRConvertModelRepeatedly) THEN
       BRConvertMode = 1
     ELSE
       BRConvertMode = -1
     END IF
    DeltaTimeBRWindow = BRConvertFluidToElectronsTime
  ELSEIF(BRConvertFluidToElectronsTime.GT.BRConvertElectronsToFluidTime)THEN
    ! Mode=2: kin -> BR -> kin (when BRConvertFluidToElectronsTime > BRConvertElectronsToFluidTime)
    IF(BRConvertModelRepeatedly) THEN
      BRConvertMode = 2
    ELSE
      BRConvertMode = -2
    END IF
    DeltaTimeBRWindow = BRConvertFluidToElectronsTime - BRConvertElectronsToFluidTime
  ELSE
    CALL abort(__STAMP__,'BRConvertFluidToElectronsTime == BRConvertElectronsToFluidTime is not allowed!')
  END IF ! BRConvertElectronsToFluidTime.GT.BRConvertFluidToElectronsTime
ELSEIF(BRConvertElectronsToFluidTime.GE.0.)THEN
  BRConvertMode = 3 ! Single Switch
  !DeltaTimeBRWindow = HUGE(1.) ! no relaxation
ELSEIF(BRConvertFluidToElectronsTime.GE.0.)THEN
  BRConvertMode = 3 ! Single Switch
  DeltaTimeBRWindow = BRConvertFluidToElectronsTime
END IF ! (BRConvertElectronsToFluidTime.GE.0.).AND.(BRConvertFluidToElectronsTime.GE.0.)

CALL PrintOption('Switch BR Electron <-> Kinetic: BRConvertMode (zero means OFF)' , 'INFO' , IntOpt=BRConvertMode)

! Depending on usage of variable BR ref. electron temperature, do a sanity check
IF(CalcBRVariableElectronTemp.AND.DeltaTimeBRWindow.LT.0.0) CALL abort(__STAMP__,&
    'DeltaTimeBRWindow is negative. Set value for BRConvertFluidToElectronsTime.')

END SUBROUTINE InitSwitchBRElectronModel


!===================================================================================================================================
!> Initialize the variables first
!===================================================================================================================================
SUBROUTINE InitializeVariablesElectronFluidRegions()
! MODULES
USE MOD_Globals
USE MOD_HDG_Vars
USE MOD_ReadInTools   ,ONLY: GETSTR,GETLOGICAL,GETREALARRAY,GETINT,GETREAL,PrintOption
USE MOD_Particle_Vars ,ONLY: nSpecies,Species
USE MOD_Restart_Vars  ,ONLY: DoRestart,RestartTime
USE MOD_TimeDisc_Vars ,ONLY: Time
USE MOD_HDF5_Input    ,ONLY: DatasetExists,ReadArray
USE MOD_IO_HDF5       ,ONLY: OpenDataFile,CloseDataFile,File_ID
USE MOD_Restart_Vars  ,ONLY: RestartFile
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(32) :: hilf, hilf2
INTEGER       :: iRegions,iSpec,iInit
REAL          :: phimax_tmp
LOGICAL       :: RegionElectronRefExists
REAL          :: RegionElectronRefHDF5(3) !< RegionElectronRefHDF5(rho0[C/m^3],phi0[V],Te[eV])|1:BRNbrOfRegions) when using
                                          !< BRAutomaticElectronRef
!===================================================================================================================================
!-- Read parameters for region mapping
BRNbrOfRegions = GETINT('BRNbrOfRegions','0')
UseBRElectronFluid = .FALSE. ! Initialize
CalcBRVariableElectronTemp = .FALSE. ! Initialize
IF(BRNbrOfRegions.GT.0)THEN
  UseBRElectronFluid = .TRUE.

  !--- Set BR electron region(s)
  ALLOCATE(BRRegionBounds(1:6,1:BRNbrOfRegions))
  DO iRegions=1,BRNbrOfRegions
    WRITE(UNIT=hilf2,FMT='(I0)') iRegions
    BRRegionBounds(1:6,iRegions) = GETREALARRAY('BRRegionBounds'//TRIM(hilf2),6,'0. , 0. , 0. , 0. , 0. , 0.')
  END DO

  !--- Create mapping of element ID to BR electron region and set reference variables
  CALL MapBRRegionToElem()
  ALLOCATE(RegionElectronRef(1:3,1:BRNbrOfRegions))
  DO iRegions=1,BRNbrOfRegions
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

  !--- Set variable reference electron temperature
  BRVariableElectronTemp = GETSTR('BRVariableElectronTemp')
  SELECT CASE(TRIM(BRVariableElectronTemp))
  CASE('constant')  ! Default, nothing to do
  CASE('linear','exp','linear-phi','exp-phi') ! Linear/Exponential drop towards background temperature
    ! BGGas temperature
    BRVariableElectronTempValue = -1. ! initialize
    DO iSpec=1,nSpecies
      DO iInit=1, Species(iSpec)%NumberOfInits
        IF(TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'background')THEN
          BRVariableElectronTempValue = Species(iSpec)%Init(iInit)%MWTemperatureIC
          CALL PrintOption('Final value for variable BR reference electron temperature','INFO',RealOpt=BRVariableElectronTempValue)
          EXIT
        END IF ! TRIM(Species(iSpec)%Init(iInit)%SpaceIC).EQ.'background')
      END DO ! iInit=1, Species(iSpec)%NumberOfInits
    END DO ! iSpec=1,nSpecies
    IF(BRVariableElectronTempValue.LT.0.0) CALL abort(__STAMP__,'Variable reference electron temperature: final value is negative.')
  CASE DEFAULT
    CALL abort(__STAMP__,'Unknown method for BRVariableElectronTemp: '//TRIM(BRVariableElectronTemp))
  END SELECT
  IF((TRIM(BRVariableElectronTemp).NE.'').AND.(TRIM(BRVariableElectronTemp).NE.'constant')) CalcBRVariableElectronTemp=.TRUE.
  IF(CalcBRVariableElectronTemp)THEN
    ALLOCATE(RegionElectronRefBackup(1:3,1:BRNbrOfRegions))
    DO iRegions=1,BRNbrOfRegions
      RegionElectronRefBackup(1,iRegions) = 0.
      RegionElectronRefBackup(2,iRegions) = RegionElectronRef(2,iRegions)
      RegionElectronRefBackup(3,iRegions) = RegionElectronRef(3,iRegions)
    END DO
  END IF ! CalcBRVariableElectronTemp

  !--- Set reference parameters from .h5 file and update automatically
  BRAutomaticElectronRef = GETLOGICAL('BRAutomaticElectronRef')
  IF(BRAutomaticElectronRef)THEN
    IF(BRNbrOfRegions.GT.1) CALL abort(__STAMP__,'BRAutomaticElectronRef is only implemented for BRNbrOfRegions = 1')
    IF(CalcBRVariableElectronTemp.AND.BRAutomaticElectronRef) &
        CALL abort(__STAMP__,'BRVariableElectronTemp AND BRAutomaticElectronRef cannot both be used at the same time.')

    ! Initialize for all procs
    RegionElectronRefHDF5 = -1.
    ! Restart: Only root reads state file to prevent access with a large number of processors
    IF(DoRestart)THEN
      IF(MPIRoot)THEN
        CALL OpenDataFile(RestartFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
        CALL DatasetExists(File_ID,'RegionElectronRef',RegionElectronRefExists)
        IF(RegionElectronRefExists)THEN
          CALL ReadArray('RegionElectronRef',2,(/1_IK,3_IK/),0_IK,2,RealArray=RegionElectronRefHDF5)
          WRITE(UNIT_stdOut,'(3(A,ES10.2E3))') " Read RegionElectronRef from restart file ["//TRIM(RestartFile)//"] rho0[C/m^3]: ",&
              RegionElectronRefHDF5(1),", phi0[V]: ",RegionElectronRefHDF5(2),", Te[eV]: ",RegionElectronRefHDF5(3)
        END IF ! RegionElectronRefExists
        CALL CloseDataFile()
      END IF ! MPIRoot
#if USE_MPI
      ! Broadcast from root to other processors
      CALL MPI_BCAST(RegionElectronRefHDF5,3, MPI_DOUBLE_PRECISION,0,MPI_COMM_PICLAS,iERROR)
#endif /*USE_MPI*/
    END IF ! DoRestart

    ! Use reference parameters if they were read from .h5
    IF((RegionElectronRefHDF5(1).GT.0.).AND.(RegionElectronRefHDF5(3).GT.0.))THEN! rho_ref and T_ref > 0
      DO iRegions=1,BRNbrOfRegions
        RegionElectronRef(1,iRegions) = RegionElectronRefHDF5(1)
        RegionElectronRef(2,iRegions) = RegionElectronRefHDF5(2)
        RegionElectronRef(3,iRegions) = RegionElectronRefHDF5(3)
      END DO
    END IF ! (RegionElectronRefHDF5(1).GT.0.).AND.(RegionElectronRefHDF5(3).GT.0.)

    ! Find all elements, where the reference values are to be extracted from
    CALL InitBRAutomaticElectronRefElements()
  END IF ! BRAutomaticElectronRef
END IF ! BRNbrOfRegions.GT.0

!--- Check whether it is a restart or a fresh computation
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

!--- Sanity Check
IF((.NOT.UseBRElectronFluid).AND.BRConvertElectronsToFluid)THEN
  SWRITE(UNIT_StdOut,*) "UseBRElectronFluid        =", UseBRElectronFluid
  SWRITE(UNIT_StdOut,*) "BRConvertElectronsToFluid =", BRConvertElectronsToFluid
  CALL abort(__STAMP__,'UseBRElectronFluid and BRConvertElectronsToFluid MUST be both true. Define BR electron fluid mode!')
END IF ! (.NOT.UseBRElectronFluid).AND.BRConvertElectronsToFluid

!--- Sanity Check
IF(UseBRElectronFluid.AND.BRConvertFluidToElectrons)THEN
  SWRITE(UNIT_StdOut,*) "UseBRElectronFluid        =", UseBRElectronFluid
  SWRITE(UNIT_StdOut,*) "BRConvertFluidToElectrons =", BRConvertFluidToElectrons
  CALL abort(__STAMP__,'UseBRElectronFluid and BRConvertFluidToElectrons CANNOT both be true. Deactivate BR electron fluid model!')
END IF ! UseBRElectronFluid.AND.BRConvertFluidToElectrons

!--- Check variable reference electron temperature
IF(CalcBRVariableElectronTemp) THEN
  ! For BR Electron / fully kinetic model switch, get the next time a switch is going to be performed
  time=RestartTime
  CALL GetNextBRSwitchTime()
  ! Depending on kinetic/BR model, set the reference electron temperature
  CALL CalculateVariableRefElectronTemp(0.)
END IF

END SUBROUTINE InitializeVariablesElectronFluidRegions


!===================================================================================================================================
!> When changing the reference electron temperature, first calculate the new temperature and the update the temperature-dependent
!> matrices for the non-linear HDG solver
!===================================================================================================================================
SUBROUTINE UpdateVariableRefElectronTemp(tShift)
! MODULES
USE MOD_HDG_vars      ,ONLY: UseBRElectronFluid,HDGNonLinSolver
USE MOD_Elem_Mat      ,ONLY: Elem_Mat,BuildPrecond
USE MOD_TimeDisc_Vars ,ONLY: iter
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)  :: tShift ! temporal shift for electron temperature calculation, calculates temperature, e.g., for t^n or t^n+1
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! Set new reference electron temperature for at t^n or t^n+1
CALL CalculateVariableRefElectronTemp(tShift)
! Calculate NonlinVolumeFac(r,iElem)=RegionElectronRef(1,RegionID) / (RegionElectronRef(3,RegionID)*eps0)
IF(UseBRElectronFluid.AND.(HDGNonLinSolver.EQ.1)) CALL UpdateNonlinVolumeFac(.FALSE.)
! Pre-compute HDG local element matrices
CALL Elem_Mat(iter)
! Build a block-diagonal preconditioner for the lambda system
CALL BuildPrecond()
END SUBROUTINE UpdateVariableRefElectronTemp


!===================================================================================================================================
!> When automatically calculating the reference potential, as well as electron density and temperature, also update the
!> temperature- and density-dependent matrices for the non-linear HDG solver
!===================================================================================================================================
SUBROUTINE UpdateBRAutomaticElectronRef()
! MODULES
USE MOD_HDG_vars               ,ONLY: HDGNonLinSolver
USE MOD_Particle_Analyze_Tools ,ONLY: CalculateElectronIonDensityCell,CalculateElectronTemperatureCell
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! ------------
! TODO: calculate n_e and T_e only in the relevant cells (1st: only communicate values between these procs, 2nd: all other procs)
! ------------
! Update electron density in each cell
CALL CalculateElectronIonDensityCell()
! Update electron temperature in each cell
CALL CalculateElectronTemperatureCell()
! Set new reference potential, electron density and temperature
CALL CalculateBRAutomaticElectronRef()
! Calculate NonlinVolumeFac(r,iElem)=RegionElectronRef(1,RegionID) / (RegionElectronRef(3,RegionID)*eps0)
IF(HDGNonLinSolver.EQ.1) CALL UpdateNonlinVolumeFac(.FALSE.)
END SUBROUTINE UpdateBRAutomaticElectronRef


!===================================================================================================================================
!> For BR Electron / fully kinetic model: set current reference electron temperature
!===================================================================================================================================
SUBROUTINE CalculateVariableRefElectronTemp(tAdd)
! MODULES
USE MOD_PreProc
USE MOD_Globals       ,ONLY: abort
USE MOD_HDG_Vars      ,ONLY: BRNbrOfRegions,UseBRElectronFluid,RegionElectronRefBackup,RegionElectronRef,DeltaTimeBRWindow
USE MOD_HDG_Vars      ,ONLY: BRVariableElectronTemp,BRVariableElectronTempValue
USE MOD_Globals_Vars  ,ONLY: ElementaryCharge,BoltzmannConst
USE MOD_TimeDisc_Vars ,ONLY: dt_Min
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN) :: tAdd
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER         :: iRegions
!===================================================================================================================================
IF(UseBRElectronFluid)THEN
  ! Sanity check
  IF(ABS(DeltaTimeBRWindow).LE.0.0) CALL abort(__STAMP__,&
      'DeltaTimeBRWindow is zero. Therefore, 1/DeltaTimeBRWindow is not possible!')

  ASSOCIATE( &
        T2   => BRVariableElectronTempValue*BoltzmannConst/ElementaryCharge ,& ! convert K to eV
        dt   => dt_Min(DT_BR_SWITCH)+tAdd                                   ,&
        Phi2 => 0.2 ) ! Volt
    ! Select scaling function
    SELECT CASE(TRIM(BRVariableElectronTemp))
    CASE('linear') ! Linear drop towards background temperature
      DO iRegions=1,BRNbrOfRegions
        RegionElectronRef(3,iRegions) = RegionElectronRefBackup(3,iRegions)
        ASSOCIATE(T1   => RegionElectronRefBackup(3,iRegions) )
          ! Fallback
          IF(dt.EQ.HUGE(1.))THEN
            RegionElectronRef(3,iRegions) = T2 ! fall back to starting temperature
          ELSE
            RegionElectronRef(3,iRegions) = ((T1-T2)/DeltaTimeBRWindow)*dt + T2
          END IF ! dt_Min(DT_BR_SWITCH).EQ.HUGE(1.)
        END ASSOCIATE
      END DO
    CASE('exp') ! Exponential drop towards background temperature
      DO iRegions=1,BRNbrOfRegions
        RegionElectronRef(3,iRegions) = RegionElectronRefBackup(3,iRegions)
        ASSOCIATE(T1   => RegionElectronRefBackup(3,iRegions) )
          ! Fallback
          IF(dt.EQ.HUGE(1.))THEN
            RegionElectronRef(3,iRegions) = T2 ! fall back to starting temperature
          ELSE
            RegionElectronRef(3,iRegions) = (T1-T2)*EXP((-DeltaTimeBRWindow+dt)*3e6) + T2
          END IF ! dt_Min(DT_BR_SWITCH).EQ.HUGE(1.)
        END ASSOCIATE
      END DO
    CASE('linear-phi') ! Linear drop towards background temperature
      DO iRegions=1,BRNbrOfRegions
        RegionElectronRef(3,iRegions) = RegionElectronRefBackup(3,iRegions)
        RegionElectronRef(2,iRegions) = RegionElectronRefBackup(2,iRegions)
        ASSOCIATE(T1   => RegionElectronRefBackup(3,iRegions) ,&
                  Phi1 => RegionElectronRefBackup(2,iRegions))
          ! Fallback
          IF(dt.EQ.HUGE(1.))THEN
            RegionElectronRef(3,iRegions) = T2 ! fall back to starting temperature
            RegionElectronRef(2,iRegions) = Phi2 ! fall back to starting voltage
          ELSE
            RegionElectronRef(3,iRegions) = ((T1-T2)/DeltaTimeBRWindow)*dt + T2
            RegionElectronRef(2,iRegions) = ((Phi1-Phi2)/DeltaTimeBRWindow)*dt + Phi2
          END IF ! dt_Min(DT_BR_SWITCH).EQ.HUGE(1.)
        END ASSOCIATE
      END DO
    CASE('exp-phi') ! Exponential drop towards background temperature
      DO iRegions=1,BRNbrOfRegions
        RegionElectronRef(3,iRegions) = RegionElectronRefBackup(3,iRegions)
        RegionElectronRef(2,iRegions) = RegionElectronRefBackup(2,iRegions)
        ASSOCIATE(T1   => RegionElectronRefBackup(3,iRegions) ,&
                  Phi1 => RegionElectronRefBackup(2,iRegions))
          ! Fallback
          IF(dt.EQ.HUGE(1.))THEN
            RegionElectronRef(3,iRegions) = T2 ! fall back to starting temperature
            RegionElectronRef(2,iRegions) = Phi2 ! fall back to starting voltage
          ELSE
            RegionElectronRef(3,iRegions) = (T1-T2)*EXP((-DeltaTimeBRWindow+dt)*3e6) + T2
            RegionElectronRef(2,iRegions) = (Phi1-Phi2)*EXP((-DeltaTimeBRWindow+dt)*3e6) + Phi2
          END IF ! dt_Min(DT_BR_SWITCH).EQ.HUGE(1.)
        END ASSOCIATE
      END DO
    END SELECT
  END ASSOCIATE

  ! Sanity check
  DO iRegions=1,BRNbrOfRegions
    IF(RegionElectronRef(3,iRegions).LT.0.0) CALL abort(__STAMP__,'Negative ref. electron temp!')
  END DO

ELSE
  DO iRegions=1,BRNbrOfRegions
    RegionElectronRef(3,iRegions) = 0.
  END DO
END IF ! UseBRElectronFluid

END SUBROUTINE CalculateVariableRefElectronTemp


!===================================================================================================================================
!> For BR Electron: set current reference electron temperature
!===================================================================================================================================
SUBROUTINE CalculateBRAutomaticElectronRef()
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_HDG_Vars              ,ONLY: nBRAverageElems,BRAverageElemToElem,BRNbrOfRegions,nBRAverageElemsGlobal,RegionElectronRef
USE MOD_Globals_Vars          ,ONLY: ElementaryCharge,BoltzmannConst
USE MOD_Particle_Mesh_Vars    ,ONLY: ElemVolume_Shared
USE MOD_Mesh_Vars             ,ONLY: offSetElem
USE MOD_ChangeBasis           ,ONLY: ChangeBasis3D
USE MOD_DG_Vars               ,ONLY: U
USE MOD_Interpolation_Vars    ,ONLY: NAnalyze,Vdm_GaussN_NAnalyze,wAnalyze
USE MOD_Mesh_Vars             ,ONLY: Elem_xGP,sJ
USE MOD_Mesh_Tools            ,ONLY: GetCNElemID
USE MOD_Particle_Analyze_Vars ,ONLY: ElectronDensityCell,ElectronTemperatureCell
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iElem,iBRElem,iRegions,CNElemID
REAL    :: RegionElectronRefNew(3),RegionElectronRefNewGlobal(3),phiElem(1:PP_nVar)
REAL    :: U_NAnalyze(1:PP_nVar,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL    :: Coords_NAnalyze(3,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL    :: J_NAnalyze(1,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL    :: J_N(1,0:PP_N,0:PP_N,0:PP_N)
INTEGER :: k,l,m
!===================================================================================================================================
RegionElectronRefNew = 0. ! Initialize

! Add up all contributions
DO iBRElem = 1, nBRAverageElems
  iElem = BRAverageElemToElem(iBRElem)
  ! n_e
  RegionElectronRefNew(1) = RegionElectronRefNew(1) + ElectronDensityCell(iElem)

  ! phi: integrate over cell volume and divide by cell volume to get the integral average value
  phiElem=0.
  ! Interpolate the physical position Elem_xGP to the analyze position, needed for exact function
  CALL ChangeBasis3D(3,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,Elem_xGP(1:3,:,:,:,iElem),Coords_NAnalyze(1:3,:,:,:))
  ! Interpolate the Jacobian to the analyze grid: be careful we interpolate the inverse of the inverse of the jacobian ;-)
  J_N(1,0:PP_N,0:PP_N,0:PP_N)=1./sJ(:,:,:,iElem)
  CALL ChangeBasis3D(1,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,J_N(1:1,0:PP_N,0:PP_N,0:PP_N),J_NAnalyze(1:1,:,:,:))
  ! Interpolate the solution to the analyze grid
  CALL ChangeBasis3D(PP_nVar,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,U(1:PP_nVar,:,:,:,iElem),U_NAnalyze(1:PP_nVar,:,:,:))
  DO m=0,NAnalyze
    DO l=0,NAnalyze
      DO k=0,NAnalyze
        phiElem(1:PP_nVar) = phiElem(1:PP_nVar) + U_NAnalyze(:,k,l,m)*wAnalyze(k)*wAnalyze(l)*wAnalyze(m)*J_NAnalyze(1,k,l,m)
      END DO ! k
    END DO ! l
  END DO ! m
  CNElemID = GetCNElemID(iElem+offSetElem)
  phiElem = phiElem / ElemVolume_Shared(CNElemID)
  RegionElectronRefNew(2) = RegionElectronRefNew(2) + phiElem(1)

  ! T_e
  RegionElectronRefNew(3) = RegionElectronRefNew(3) + ElectronTemperatureCell(iElem)
END DO ! iBRElem = 1, nBRAverageElems

! Sum the number of elements (later required for averaging the cell-constant values globally on each processor)
#if USE_MPI
CALL MPI_ALLREDUCE(RegionElectronRefNew , RegionElectronRefNewGlobal , 3 , MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_PICLAS , IERROR)
#else
RegionElectronRefNewGlobal = RegionElectronRefNew
#endif /*USE_MPI*/

! Calculate the new values
RegionElectronRefNewGlobal    = RegionElectronRefNewGlobal / nBRAverageElemsGlobal
RegionElectronRefNewGlobal(1) = RegionElectronRefNewGlobal(1) * ElementaryCharge ! convert 1/m^3 to C/m^3
RegionElectronRefNewGlobal(3) = RegionElectronRefNewGlobal(3) * BoltzmannConst/ElementaryCharge ! convert K to eV

! Assign new values
DO iRegions=1,BRNbrOfRegions
  RegionElectronRef(1,iRegions) = RegionElectronRefNewGlobal(1)
  RegionElectronRef(2,iRegions) = RegionElectronRefNewGlobal(2)
  RegionElectronRef(3,iRegions) = RegionElectronRefNewGlobal(3)
END DO

END SUBROUTINE CalculateBRAutomaticElectronRef


!===================================================================================================================================
!> For BR Electron / fully kinetic model switch, get the next time a switch is going to be performed
!> Either one/two or multiple switches are possible depending on the user settings
!===================================================================================================================================
SUBROUTINE GetNextBRSwitchTime()
! MODULES
USE MOD_Globals       ,ONLY: abort,mpiroot
USE MOD_TimeDisc_Vars ,ONLY: dt_Min,Time,tEnd
USE MOD_HDG_Vars      ,ONLY: BRConvertFluidToElectronsTime,BRConvertElectronsToFluidTime,BRConvertMode
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: t             !> time in relative frame
REAL :: tMin,tMax     !> sorted t1 and t2
REAL :: tBRSwitchDiff !> delta to next switch (BR or kinetic), in relative time frame, i.e., MOD
REAL :: tBRSwitch     !> time at next switch (BR or kinetic), in relative time frame, i.e., MOD
LOGICAL::debug
!===================================================================================================================================
debug=.false.
!tBRSwitchDiff = tBRSwitch-Time ! Time to BR<->kinetic switch, use extra variable so number doesn't change due to numerical errors
ASSOCIATE( t1 => BRConvertFluidToElectronsTime,&
           t2 => BRConvertElectronsToFluidTime &
         )
  SELECT CASE (BRConvertMode)
    CASE (3,-1,-2) ! Single switch OR 2 switches: BR -> kin -> BR (t1 < t2) OR kin -> BR -> kin (t1 > t2)
      t = Time
    CASE (1)  ! Multiple switches: BR -> kin -> BR (t1 < t2)
      t = MOD(Time,t2)
      IF(ALMOSTEQUALRELATIVE(t,t2,1e-5)) t=0. ! Prevent that modulus returns t ~= t2 (due to tolerance)
    CASE (2)  ! Multiple switches: kin -> BR -> kin (t1 > t2)
      IF(debug.and.mpiroot) WRITE (*,*) "Time,t1 =", Time,t1
      t = MOD(Time,t1)
      IF(ALMOSTEQUALRELATIVE(t,t1,1e-5)) t=0. ! Prevent that modulus returns t ~= t1 (due to tolerance)
    CASE DEFAULT
      CALL abort(__STAMP__,'Unknown value for BRConvertMode =',IntInfoOpt=BRConvertMode)
  END SELECT

  ! Calculate delta time
  IF((t1.GT.0.0).AND.(t2.GT.0.0))THEN
    tMin=MIN(t1,t2)
    tMax=MAX(t1,t2)
    IF(debug.and.mpiroot) WRITE (*,*) "t,tmin,t.GE.tMin                   =", t,tmin,t.GE.tMin,ALMOSTEQUALRELATIVE(t,tMin,1e-5)
    IF((t.GE.tMin).OR.(ALMOSTEQUALRELATIVE(t,tMin,1e-5)))THEN
      tBRSwitch = tMax
    ELSE
      tBRSwitch = tMin
    END IF ! t.GT.t1
  ELSEIF(t1.GT.0.0)THEN
    tBRSwitch = t1
  ELSEIF(t2.GT.0.0)THEN
    tBRSwitch = t2
  END IF ! (t1.GT.0.0).AND.(t2.GT.0)

  ! Check if the switch time has already been passed, i.e., catch negative deltas
  IF((tBRSwitch.LT.t).AND.(tEnd.GE.t)) tBRSwitch = tEnd

  ! Catch equal times
  IF((ALMOSTEQUALRELATIVE(tBRSwitch,t,1e-5)).AND.(tEnd.GE.t)) tBRSwitch = tEnd

  IF(debug.and.mpiroot) WRITE (*,*) "tBRSwitch,t i                      =", tBRSwitch,t
  ! Catch tolerance issue which leads to a timestep of 1.0123123E-23 (basically zero)
  IF(ALMOSTEQUALRELATIVE(tBRSwitch,t,1e-5))THEN
    tBRSwitchDiff = 0.0
  ELSE
    tBRSwitchDiff = tBRSwitch - t
  END IF

  ! Final sanity check: Set dt_Min(DT_BR_SWITCH) and catch 0.0 or negative time delta
  IF(tBRSwitchDiff.GT.0.0)THEN
    dt_Min(DT_BR_SWITCH) = tBRSwitchDiff
  ELSE
    dt_Min(DT_BR_SWITCH) = HUGE(1.)
  END IF ! tBRSwitchDiff.GT.0.0
END ASSOCIATE

IF(debug.and.mpiroot) WRITE (*,*) "tBRSwitchDiff,dt_Min(DT_BR_SWITCH) =", tBRSwitchDiff,dt_Min(DT_BR_SWITCH)

END SUBROUTINE GetNextBRSwitchTime


!===================================================================================================================================
!> Switch between BR electron fluid and kinetic model
!===================================================================================================================================
SUBROUTINE SwitchBRElectronModel()
! MODULES
USE MOD_Globals
USE MOD_HDG_Vars
USE MOD_TimeDisc_Vars    ,ONLY: time,iter
USE MOD_Elem_Mat         ,ONLY: Elem_Mat,BuildPrecond
USE MOD_part_operations  ,ONLY: RemoveAllElectrons
USE MOD_DSMC_ChemInit    ,ONLY: InitReactionPaths
USE MOD_DSMC_Vars        ,ONLY: ChemReac,CollInf,UseDSMC,CollisMode
USE MOD_MCC_Vars         ,ONLY: XSec_NullCollision, SpecXSec
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: nLoadBalanceSteps
#endif /*USE_LOADBALANCE*/
USE MOD_MCC_Init         ,ONLY: MCC_Chemistry_Init
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL :: debug,SwitchToBR,SwitchToKin
INTEGER :: iCase
!===================================================================================================================================
!debug=.true.
debug=.false.

! check if a switch happens now to update the variable reference electron temperature or activate chemical reactions with electrons
SwitchToBR=.FALSE.
SwitchToKin=.FALSE.

ASSOCIATE( tBR2Kin => BRConvertFluidToElectronsTime ,&
           tKin2BR => BRConvertElectronsToFluidTime )
  ! BR -> kinetic
  IF(UseBRElectronFluid.AND.tBR2Kin.GT.0.0)THEN
    IF((.NOT.BRConvertModelRepeatedly).AND.(time.GE.tBR2Kin)                                             .OR.&
       ((BRConvertMode.EQ.1)          .AND.(GreaterEqualWithTolerance(MOD(time,tKin2BR),tBR2Kin)))       .OR.&
       ((BRConvertMode.EQ.2)          .AND.(LesserThanWithTolerance(MOD(time,tBR2Kin),tKin2BR,tBR2Kin))) )THEN
      IF(debug)THEN
        IPWRITE(UNIT_StdOut,*) "\nSWITCH TO kinetic ?",time,"\n"
        read*
      END IF ! debug
      CALL CreateElectronsFromBRFluid(.FALSE.) ! Use BR electron fluid model density to create kinetic electrons in each cell
      UseBRElectronFluid = .FALSE. ! Deactivate BR fluid
      CALL Elem_Mat(iter)          ! Recompute elem matrices
      SwitchToKin=.TRUE.! check if a switch happens now to update the chemical reactions (activate electron products)
      IF((.NOT.BRConvertModelRepeatedly).AND.(BRConvertMode.EQ.-2)) tKin2BR = -1.0 ! deactivate kin -> BR
      ! (Re-)activate Null-Collision (if the read-in parameter is set true)
      XSec_NullCollision = BRNullCollisionDefault
    END IF
  ENDIF

  ! kinetic -> BR
  IF((.NOT.UseBRElectronFluid).AND.(tKin2BR.GT.0.0).AND.(.NOT.SwitchToKin))THEN
    IF((.NOT.BRConvertModelRepeatedly).AND.(time.GE.tKin2BR)                                             .OR.&
       ((BRConvertMode.EQ.1)          .AND.(LesserThanWithTolerance(MOD(time,tKin2BR),tBR2Kin,tKin2BR))) .OR.&
       ((BRConvertMode.EQ.2)          .AND.(GreaterEqualWithTolerance(MOD(time,tBR2Kin),tKin2BR)))       )THEN
      IF(debug)THEN
        IPWRITE(UNIT_StdOut,*) "\nMOD(time,tKin2BR),tBR2Kin,MOD(time,tKin2BR).LT.tBR2Kin =", MOD(time,tKin2BR),tBR2Kin,MOD(time,tKin2BR).LT.tBR2Kin
        IPWRITE(UNIT_StdOut,*) "SWITCH TO BR ?",time,"\n"
        read*
      END IF ! debug
      IF(BRNbrOfRegions.EQ.0) CALL abort(__STAMP__,'SwitchBRElectronModel(): Cannot switch [kin -> BR] as no BR regions are defined!')
      ! When switching from kin to BR and using automatic ref. value determination, update ref. values here and subsequent matrices
      IF(BRAutomaticElectronRef) CALL UpdateBRAutomaticElectronRef()
      CALL RemoveAllElectrons()  ! Remove all electron particles from the simulation
      UseBRElectronFluid = .TRUE.! Activate BR fluid (must be set after UpdateBRAutomaticElectronRef because of determination of Te)
      ! Recompute the matrices: also consider BRAutomaticElectronRef here
      IF(.NOT.CalcBRVariableElectronTemp)THEN! for variable Te, these matrices are updated in UpdateVariableRefElectronTemp
        CALL Elem_Mat(iter) ! Recompute elem matrices
        CALL BuildPrecond() ! Build a block-diagonal preconditioner for the lambda system
      END IF
      SwitchToBR=.TRUE.! check if a switch happens now to update the variable reference electron temperature
      IF((.NOT.BRConvertModelRepeatedly).AND.(BRConvertMode.EQ.-1))tBR2Kin = -1.0 ! deactivate BR -> kin
      ! Recompute lambda: force iteration
      !CALL  RecomputeLambda(time)
      ! Deactivate Null-Collision
      XSec_NullCollision = .FALSE.
    END IF
  END IF ! .NOT.UseBRElectronFluid.AND.BRConvertE.GT.0.0
END ASSOCIATE

! For BR Electron / fully kinetic model switch, get the next time a switch is going to be performed
CALL GetNextBRSwitchTime()

! Depending on kinetic/BR model, update values and matrices
IF(SwitchToBR.AND.CalcBRVariableElectronTemp) CALL UpdateVariableRefElectronTemp(0.)

! Update reaction paths (specifically the ones that involve electrons, which are deactivated for UseBRElectronFluid = .FALSE.)
IF(UseDSMC)THEN
  IF((SwitchToBR.OR.SwitchToKin).AND.(CollisMode.EQ.3))THEN
    DO iCase = 1, CollInf%NumCase
      SDEALLOCATE(ChemReac%CollCaseInfo(iCase)%ReactionIndex)
      SDEALLOCATE(ChemReac%CollCaseInfo(iCase)%ReactionProb)
    END DO
    ChemReac%CollCaseInfo(:)%NumOfReactionPaths = 0 ! Re-initialize
    CALL InitReactionPaths()

    ! Initialize MCC model: Read-in of the reaction cross-section database and re-calculation of the null collision probability
    IF(ChemReac%AnyXSecReaction) THEN
      DO iCase = 1, CollInf%NumCase
        SDEALLOCATE(SpecXSec(iCase)%ReactionPath)
      END DO
      CALL MCC_Chemistry_Init()
    END IF

  END IF ! (SwitchToBR.OR.SwitchToKin).AND.(CollisMode.EQ.3)
END IF ! UseDSMC


#if USE_LOADBALANCE
! When switching BR <-> kin, reset the number of load balances to 0
IF(SwitchToBR.OR.SwitchToKin)THEN
  nLoadBalanceSteps = 0
  SWRITE (*,*) " Switching BR <-> kin: Setting nLoadBalanceSteps=0"
END IF ! SwitchToBR.OR.SwitchToKin
#endif /*USE_LOADBALANCE*/

END SUBROUTINE SwitchBRElectronModel


!===================================================================================================================================
!> Check if a >= b or a is almost equal to b via ALMOSTEQUALRELATIVE
!> Catch tolerance issues when a is only an epsilon smaller than b but the inquiry should be that they are equal
!===================================================================================================================================
PPURE LOGICAL FUNCTION GreaterEqualWithTolerance(a,b)
! MODULES
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


!===================================================================================================================================
!> Check if a < b and NOT a is almost equal to b via ALMOSTEQUALRELATIVE
!> Catch tolerance issues when a<b returns a false positive, because the numbers are actually the same (with an epsilon difference)
!> Example: a = 5.4999999999999995E-006
!>          b = 5.4999999999999996E-006
!===================================================================================================================================
PPURE LOGICAL FUNCTION LesserThanWithTolerance(a,b,c)
! MODULES
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


!===================================================================================================================================
!> 1.) reconstruct the electron phase space using the integrated charge density in each cell
!>     a.) from ElectronDensityCell and ElectronTemperatureCell that are read from .h5 state file (only during restart)
!>     b.) from BR electron model variables in each cell (only during the simulation)
!===================================================================================================================================
SUBROUTINE CreateElectronsFromBRFluid(CreateFromRestartFile)
! MODULES
USE MOD_Globals
USE MOD_Globals             ,ONLY: abort,MPIRoot,UNIT_stdOut,IK,MPI_COMM_PICLAS
USE MOD_Globals_Vars        ,ONLY: ElementaryCharge,BoltzmannConst
USE MOD_PreProc
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
USE MOD_Part_Tools          ,ONLY: UpdateNextFreePosition, GetNextFreePosition
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
CHARACTER(1)                   :: hilf2
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
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)
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
  hilf2='T'
ELSE
  hilf2='F'
END IF ! CreateFromRestartFile

! ---------------------------------------------------------------------------------------------------------------------------------
! 1.) reconstruct electrons
! ---------------------------------------------------------------------------------------------------------------------------------
SWRITE(UNIT_stdOut,'(A,ES25.14E3,A)')' CreateElectronsFromBRFluid(): Reconstructing electrons at t=',time,&
                                     ' from BR electron fluid density in each cell (CreateFromRestartFile='//hilf2//')'

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
IF (ElecSpecIndx.LE.0) CALL abort(__STAMP__,'Electron species not found. Cannot create electrons without the defined species!')

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

    ParticleIndexNbr            = GetNextFreePosition()

    !Set new SpeciesID of new particle (electron)
    PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
    PDM%isNewPart(ParticleIndexNbr) = .TRUE.
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
    PEM%GlobalElemID(ParticleIndexNbr)     = iElem + offsetElem
    PEM%LastGlobalElemID(ParticleIndexNbr) = PEM%GlobalElemID(ParticleIndexNbr)

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
CALL MPI_ALLREDUCE(MPI_IN_PLACE,BRNbrOfElectronsCreated,1,MPI_INTEGER,MPI_SUM,MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/

SWRITE(UNIT_StdOut,'(A,I0,A)') '  Created a total of ',BRNbrOfElectronsCreated,' electrons.'
!read*

! Update
CALL UpdateNextFreePosition()


END SUBROUTINE CreateElectronsFromBRFluid


!===================================================================================================================================
!> map a particle region to element
!> check only element barycenter, nothing else
!===================================================================================================================================
SUBROUTINE MapBRRegionToElem()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_HDG_Vars           ,ONLY: BRNbrOfRegions, BRRegionBounds,ElemToBRRegion
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
  DO iRegions=1,BRNbrOfRegions
    IF ((ElemBaryNGeo(1,iElem).LT.BRRegionBounds(1,iRegions)).OR.(ElemBaryNGEO(1,iElem).GE.BRRegionBounds(2,iRegions))) CYCLE
    IF ((ElemBaryNGeo(2,iElem).LT.BRRegionBounds(3,iRegions)).OR.(ElemBaryNGEO(2,iElem).GE.BRRegionBounds(4,iRegions))) CYCLE
    IF ((ElemBaryNGeo(3,iElem).LT.BRRegionBounds(5,iRegions)).OR.(ElemBaryNGEO(3,iElem).GE.BRRegionBounds(6,iRegions))) CYCLE
    IF (ElemToBRRegion(iElem).EQ.0) THEN
      ElemToBRRegion(iElem)=iRegions
    ELSE
      CALL ABORT(__STAMP__,'Defined regions are overlapping')
    END IF
  END DO ! iRegions=1,BRNbrOfRegions
END DO ! iElem=1,PP_nElems
END SUBROUTINE MapBRRegionToElem


!===================================================================================================================================
!> Set NonlinVolumeFac for each element depending on
!===================================================================================================================================
SUBROUTINE UpdateNonlinVolumeFac(NullifyField)
! MODULES
USE MOD_PreProc
USE MOD_HDG_Vars     ,ONLY: NonlinVolumeFac,RegionElectronRef,ElemToBRRegion
USE MOD_Globals_Vars ,ONLY: eps0
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
LOGICAL,INTENT(IN)      :: NullifyField !< Set NonlinVolumeFac = 0 if this variable is true
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i,j,k,r,iElem,RegionID
!===================================================================================================================================
IF(NullifyField) THEN
  NonlinVolumeFac = 0.
ELSE
  DO iElem=1,PP_nElems
    RegionID=ElemToBRRegion(iElem)
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      r=k*(PP_N+1)**2+j*(PP_N+1) + i+1
      NonlinVolumeFac(r,iElem)=RegionElectronRef(1,RegionID) / (RegionElectronRef(3,RegionID)*eps0)
    END DO; END DO; END DO !i,j,k
  END DO !iElem
END IF ! NewtonAdaptStartValue

END SUBROUTINE UpdateNonlinVolumeFac


!===================================================================================================================================
!> Find all elements, where the reference values are to be extracted from
!===================================================================================================================================
SUBROUTINE InitBRAutomaticElectronRefElements()
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_HDF5_Output_Fields ,ONLY: WriteBRAverageElemToHDF5
USE MOD_HDG_Vars           ,ONLY: nBRAverageElems,nBRAverageElemsGlobal,BRAverageElemToElem
USE MOD_Mesh_Tools         ,ONLY: GetCNElemID
USE MOD_Mesh_Vars          ,ONLY: offSetElem
USE MOD_Particle_Mesh_Vars ,ONLY: ElemCharLength_Shared,NodeCoords_Shared,ElemInfo_Shared
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL :: isBRAverageElem(1:PP_nElems) ! < Flag every element T/F if it belongs to the averaging region
INTEGER :: iElem,CNElemID,first,iBRElem
REAL    :: tol
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)') ' GET ELEMENTS TO CALCULATE AVERAGE BR ELECTRON REFERENCE PARAMETERS ...'
isBRAverageElem = .FALSE. ! Initialize
nBRAverageElems = 0 ! Initialize
! Loop all elements and find the ones that are in the averaging region
DO iElem=1,PP_nElems
  ! Check min value in x and y for each element (from nodes)
  CNElemID = GetCNElemID(iElem+offSetElem)
  tol      = ElemCharLength_Shared(CNElemID)/1000.
  first    = ElemInfo_Shared(ELEM_FIRSTNODEIND,CNElemID)
  IF(MINVAL(NodeCoords_Shared(1,first+1:first+8)).GT.tol) CYCLE
  IF(MINVAL(NodeCoords_Shared(2,first+1:first+8)).GT.tol) CYCLE
  isBRAverageElem(iElem) = .TRUE.
  nBRAverageElems        = nBRAverageElems + 1
END DO!iElem

! Create mapping
IF(nBRAverageElems.GT.0)THEN
  ALLOCATE(BRAverageElemToElem(nBRAverageElems))
  iBRElem=0
  DO iElem = 1, PP_nElems
    IF(isBRAverageElem(iElem))THEN
      iBRElem = iBRElem + 1
      BRAverageElemToElem(iBRElem) = iElem
    END IF ! isBRAverageElem(iElem)
  END DO ! iElem = 1, PP_nElems
END IF ! nBRAverageElems.GT.0

! Create .h5 file with info, which elements are within the averaging region
CALL WriteBRAverageElemToHDF5(isBRAverageElem)

! Sum the number of elements (later required for averaging the cell-constant values globally on each processor)
#if USE_MPI
CALL MPI_ALLREDUCE(nBRAverageElems , nBRAverageElemsGlobal , 1 , MPI_INTEGER , MPI_SUM , MPI_COMM_PICLAS , IERROR)
#else
nBRAverageElemsGlobal = nBRAverageElems
#endif /*USE_MPI*/
SWRITE(UNIT_stdOut,'(A,I10,A)') ' Found a total of',nBRAverageElemsGlobal,&
    ' elements for the calculation of average BR electron reference parameters.'

! If no elements are found, abort the simulation (sanity check)
IF(MPIRoot)THEN
  IF(nBRAverageElemsGlobal.EQ.0)THEN
    CALL abort(__STAMP__,'InitBRAutomaticElectronRefElements(): Found zero elements for averaging the BR electron reference parameters')
  END IF ! nBRAverageElemsGlobal.EQ.0
END IF ! MPIRoot

SWRITE(UNIT_stdOut,'(A)') ' ... AVERAGE BR ELECTRON REFERENCE PARAMETERS INITIALIZATION DONE.'
END SUBROUTINE InitBRAutomaticElectronRefElements


#endif /*defined(PARTICLES) && USE_HDG*/
END MODULE MOD_Part_BR_Elecron_Fluid
