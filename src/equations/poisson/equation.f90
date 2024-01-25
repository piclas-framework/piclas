!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_Equation
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitEquation
  MODULE PROCEDURE InitEquation
END INTERFACE
INTERFACE ExactFunc
  MODULE PROCEDURE ExactFunc
END INTERFACE
INTERFACE CalcSource
  MODULE PROCEDURE CalcSource
END INTERFACE
INTERFACE CalcSourceHDG
  MODULE PROCEDURE CalcSourceHDG
END INTERFACE
INTERFACE DivCleaningDamping
  MODULE PROCEDURE DivCleaningDamping
END INTERFACE
INTERFACE FinalizeEquation
  MODULE PROCEDURE FinalizeEquation
END INTERFACE

PUBLIC :: InitEquation,ExactFunc,CalcSource,FinalizeEquation, CalcSourceHDG,DivCleaningDamping
PUBLIC :: DefineParametersEquation
!===================================================================================================================================
CONTAINS

!==================================================================================================================================
!> Define parameters for equation
!==================================================================================================================================
SUBROUTINE DefineParametersEquation()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
#if defined(PARTICLES)
USE MOD_HDG_Vars    ,ONLY: CPPDataLength
USE MOD_ReadInTools ,ONLY: addStrListEntry
#endif /*defined(PARTICLES)*/
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Equation")
CALL prms%CreateIntOption(      'IniExactFunc'     , 'Define exact function necessary for linear scalar advection')
CALL prms%CreateRealArrayOption('RefState'         , 'State(s) for electric potential (amplitude, frequency and phase shift).', multiple=.TRUE., no=3)
CALL prms%CreateRealArrayOption('IniWavenumber'    , 'TODO-DEFINE-PARAMETER' , '1. , 1. , 1.')
CALL prms%CreateRealArrayOption('IniCenter'        , 'TODO-DEFINE-PARAMETER' , '0. , 0. , 0.')
CALL prms%CreateRealOption(     'IniAmplitude'     , 'TODO-DEFINE-PARAMETER' , '0.1')
CALL prms%CreateRealOption(     'IniHalfwidth'     , 'TODO-DEFINE-PARAMETER' , '0.1')

CALL prms%CreateIntOption(      'chitensWhichField', 'TODO-DEFINE-PARAMETER', '-1')
CALL prms%CreateRealOption(     'chitensValue'     , 'TODO-DEFINE-PARAMETER', '-1.0')
CALL prms%CreateRealOption(     'chitensRadius'    , 'TODO-DEFINE-PARAMETER', '-1.0')

CALL prms%CreateIntOption(      'AlphaShape'       , 'TODO-DEFINE-PARAMETER', '2')
CALL prms%CreateRealOption(     'r_cutoff'         , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Modified for curved and shape-function influence (c*dt*SafetyFactor+r_cutoff)' , '1.0')

! Special BC with linear potential ramp (constant in time)
CALL prms%CreateRealArrayOption('LinPhiBasePoint'  , 'Origin(s) of coordinate system for linear potential ramp for BoundaryType = (/7,1/)'       , multiple=.TRUE. , no=3 )
CALL prms%CreateRealArrayOption('LinPhiNormal'     , 'Normal(s) vector of coordinate system for linear potential ramp for BoundaryType = (/7,1/)', multiple=.TRUE. , no=3 )
CALL prms%CreateRealOption(     'LinPhiHeight'     , 'Interval(s) for ramping from 0 to LinPhi potential ramp for BoundaryType = (/7,1/)'        , multiple=.TRUE.        )
CALL prms%CreateRealOption(     'LinPhi'           , 'Potential(s) for ramping from 0 for BoundaryType = (/7,1/)'                                , multiple=.TRUE.        )

#if defined(PARTICLES)
! Special BC with floating potential that is defined by the absorbed power of the charged particles
CALL prms%CreateRealArrayOption('CoupledPowerPotential' , 'Controlled power input: Supply vector of form (/min, start, max/) for the minimum, start (t=0) and maximum electric potential that is applied at BoundaryType = (/2,2/).', no=CPPDataLength )
CALL prms%CreateRealOption(     'CoupledPowerTarget'    , 'Controlled power input: Target input power to which the electric potential is adjusted for BoundaryType = (/2,2/)' )
CALL prms%CreateRealOption(     'CoupledPowerRelaxFac'  , 'Relaxation factor for calculation of new electric potential due to defined Target input power. Default = 0.05 (which is 5%)', '0.05' )
CALL prms%CreateRealOption(     'CoupledPowerFrequency' , 'Adaption frequency for coupled power using the integrated power for adjustment', '0.0' )

CALL prms%CreateIntFromStringOption('CoupledPowerMode', 'Define mode used for coupled power adjustment:\n'//&
                                                        ' instantaneous (1): xxx\n'//&
                                                        ' moving-average (2): xxx\n'//&
                                                        ' integrated (3): xxx\n','instantaneous')
CALL addStrListEntry('CoupledPowerMode' , 'instantaneous' , 1)
CALL addStrListEntry('CoupledPowerMode' , 'moving-average', 2)
CALL addStrListEntry('CoupledPowerMode' , 'integrated'    , 3)
#endif /*defined(PARTICLES)*/

END SUBROUTINE DefineParametersEquation


SUBROUTINE InitEquation()
!===================================================================================================================================
! Init Poisson euqation system
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Equation_Vars
USE MOD_HDG_vars
USE MOD_Globals_Vars       ,ONLY: PI
USE MOD_ReadInTools        ,ONLY: GETREALARRAY,GETREAL,GETINT,CountOption
USE MOD_Interpolation_Vars ,ONLY: InterpolationInitIsDone
USE MOD_Mesh_Vars          ,ONLY: BoundaryName,BoundaryType,nBCs
USE MOD_Mesh_Vars          ,ONLY: nSides
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: chitensValue,chitensRadius ! Deprecated variables, remove in future (by the end of 2017)
INTEGER            :: chitensWhichField          ! Deprecated variables, remove in future (by the end of 2017)
INTEGER            :: iState                     ! i-th RefState
INTEGER            :: i,BCType,BCState
CHARACTER(LEN=255) :: BCName
INTEGER            :: nRefStateMax
INTEGER            :: nLinState,nLinStateMax
INTEGER,PARAMETER  :: BCTypeRefstate(1:4)=(/5,51,52,60/)
CHARACTER(LEN=32)  :: hilf
!===================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.EquationInitIsDone)THEN
   LBWRITE(*,*) "InitPoisson not ready to be called or already called."
   RETURN
END IF
LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT POISSON...'

! Read in boundary parameters
IniExactFunc = GETINT('IniExactFunc')

! Sanity checks
SELECT CASE (IniExactFunc)
CASE(800,801,900) ! Dielectric slab on electrode (left) with plasma between slab and other electrode opposite
#if ! (defined(CODE_ANALYZE) && USE_PETSC && PARTICLES)
  CALL abort(__STAMP__,'IniExactFunc=800,801,900 requires PICLAS_CODE_ANALYZE=ON, PICLAS_PETSC=ON and PICLAS_PARTICLES=ON')
#endif /*! (defined(CODE_ANALYZE) && USE_PETSC && PARTICLES)*/
END SELECT

! Sanity Check BCs
nRefStateMax = 0
nLinStateMax = 0
! Defaults
DO i=1,nBCs
  BCType  = BoundaryType(i,BC_TYPE)
  BCState = BoundaryType(i,BC_STATE)
  BCName  = BoundaryName(i)
  IF(ANY(BCType.EQ.BCTypeRefstate).AND.BCState.LE.0) THEN
    SWRITE(*,'(A)') "Error found for the following boundary condition"
    SWRITE(*,'(A,I0)') "   BC: ",i
    SWRITE(*,'(A,I0)') " Type: ",BCType
    SWRITE(*,'(A,I0)') "State: ",BCState
    SWRITE(*,'(A)')    " Name: "//TRIM(BCName)
    WRITE(UNIT=hilf,FMT='(I0)') BCType
    CALL abort(__STAMP__,'BCState is <= 0 for BCType='//TRIM(hilf)//' is not allowed! Set a positive integer for the n-th RefState')
  ELSEIF(ANY(BCType.EQ.BCTypeRefstate).AND.BCState.GT.0)THEN
    nRefStateMax = MAX(nRefStateMax,BCState)
  ELSEIF(BCType.EQ.7.AND.BCState.GT.0)THEN
    nLinStateMax = MAX(nLinStateMax,BCState)
  END IF
END DO

#if defined(PARTICLES)
! Coupled Power Potential: Adjust the electric potential on a BC to meet a specific power absorbed by the charged particles
CALL InitCoupledPowerPotential()
#endif /*defined(PARTICLES)*/

! Read linear potential boundaries
nLinState = CountOption('LinPhi')
IF(nLinStateMax.GT.nLinState)THEN
  SWRITE(*,'(A,I0)') "nLinStateMax: ",nLinStateMax
  SWRITE(*,'(A,I0,A)') "   nLinState: ",nLinState," (number of times LinPhi = X occurrences in the parameter file)"
  CALL abort(__STAMP__&
      ,'nLinStateMax > nLinState: The given LinState number for boundary type 7 is larger than the supplied LinPhi values. '//&
       'Define the correct number of LinPhi via, e.g., \n\n  LinPhi = 120.0 ! LinPhi Nbr 1: Voltage\n  LinPhi = 500.0 '//&
       '! LinPhi Nbr 2: Voltage\n and the corresponding LinPhiBasePoint, LinPhiNormal and LinPhiHeight parameters')
END IF ! nLinStateMax.GT.nLinState
IF(nLinState .GT. 0)THEN
  ALLOCATE(LinPhiBasePoint(3,nLinState))
  ALLOCATE(LinPhiNormal(3,nLinState))
  ALLOCATE(LinPhiHeight(nLinState))
  ALLOCATE(LinPhi(nLinState))
  DO iState=1,nLinState
    ! Read linear potential parameters
    LinPhiBasePoint(1:3,iState) = GETREALARRAY('LinPhiBasePoint',3)
    LinPhiNormal(1:3,iState)    = GETREALARRAY('LinPhiNormal',3)
    LinPhiNormal(1:3,iState)    = UNITVECTOR(LinPhiNormal(1:3,iState))
    LinPhiHeight(iState)        = GETREAL('LinPhiHeight')
    LinPhi(iState)              = GETREAL('LinPhi')
  END DO
END IF

! Read Boundary information / RefStates / perform sanity check
nRefState=CountOption('RefState')
IF(nRefStateMax.GT.nRefState)THEN
  SWRITE(*,'(A,I0)') "nRefStateMax: ",nRefStateMax
  SWRITE(*,'(A,I0,A)') "   nRefState: ",nRefState," (number of times RefState = (/x,x,x/) occurrences in the parameter file)"
  CALL abort(__STAMP__&
      ,'nRefStateMax > nRefState: The given RefState number for boundary type 5 is larger than the supplied RefStates. '//&
       'Define the correct number of RefStates via, e.g., \n\n  RefState = (/100.0 , 13.56E6 , -1.57079632679/) '//&
       '! RefState Nbr 1: Voltage, Frequency and Phase shift\n  RefState = (/ 50.0 , 13.56E6 ,  1.57079632679/) '//&
       '! RefState Nbr 2: Voltage, Frequency and Phase shift\n')
END IF ! nRefStateMax.GT.nRefState

IF(nRefState .GT. 0)THEN
  ALLOCATE(RefState(3,nRefState))
  DO iState=1,nRefState
    RefState(1:3,iState) = GETREALARRAY('RefState',3)
  END DO
END IF


! Read the velocity vector from ini file
IniWavenumber = GETREALARRAY('IniWavenumber',3,'1.,1.,1.')
IniCenter     = GETREALARRAY('IniCenter',3,'0.,0.,0.')
IniAmplitude  = GETREAL('IniAmplitude','0.1')
IniHalfwidth  = GETREAL('IniHalfwidth','0.1')

chitensWhichField = GETINT( 'chitensWhichField','-1')
chitensValue      = GETREAL('chitensValue','-1.0')
chitensRadius     = GETREAL('chitensRadius','-1.0')
IF(chitensWhichField.GT.0.0.OR.&
   chitensValue     .GT.0.0.OR.&
   chitensRadius    .GT.0.0)THEN
  CALL abort(__STAMP__,'chitensWhichField, chitensValue and chitensRadius are no longer supported. Deactivate them!')
END IF
ALLOCATE(chitens(3,3,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
ALLOCATE(chitensInv(3,3,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
ALLOCATE(chitens_face(3,3,0:PP_N,0:PP_N,nSides))

! initialize
chitens=0.
chitens(1,1,:,:,:,:)=1.
chitens(2,2,:,:,:,:)=1.
chitens(3,3,:,:,:,:)=1.
chitensInv=chitens
chitens_face=0.
chitens_face(1,1,:,:,:)=1.
chitens_face(2,2,:,:,:)=1.
chitens_face(3,3,:,:,:)=1.

alpha_shape = GETINT('AlphaShape','2')
rCutoff     = GETREAL('r_cutoff','1.')
! Compute factor for shape function
ShapeFuncPrefix = 1./(2. * beta(1.5, REAL(alpha_shape) + 1.) * REAL(alpha_shape) + 2. * beta(1.5, REAL(alpha_shape) + 1.)) &
                * (REAL(alpha_shape) + 1.)/(PI*(rCutoff**3))

ALLOCATE(E(1:3,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
E=0.


EquationInitIsDone=.TRUE.
LBWRITE(UNIT_stdOut,'(A)')' INIT POISSON DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitEquation


#if defined(PARTICLES)
!===================================================================================================================================
!> Initialize containers for coupled power potential (CPP) for adjusting the electric potential on a BC to meet a specific power
!> absorbed by the charged particles
!>
!> 1.) Get global number of coupled power potential boundaries in [1:nBCs]
!> 2.) Get parameters
!> 3.) Establish sub-communicator (all BCs directly connected to a coupled power potential boundary); MPIRoot is always in the group)
!> 4.) When restarting, load the last known values of CoupledPowerPotential(1:3) from the .h5 state file
!===================================================================================================================================
SUBROUTINE InitCoupledPowerPotential()
! MODULES
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
USE MOD_Globals          ,ONLY: CollectiveStop
USE MOD_HDG_Vars         ,ONLY: UseCoupledPowerPotential,CoupledPowerPotential,CoupledPowerTarget,CoupledPowerRelaxFac
USE MOD_HDG_Vars         ,ONLY: CPPDataLength,CoupledPowerMode,CoupledPowerFrequency
USE MOD_ReadInTools      ,ONLY: GETREALARRAY,GETREAL,GETINTFROMSTR,CountOption
USE MOD_Mesh_Vars        ,ONLY: BoundaryType,nBCs
#if USE_MPI
USE MOD_Globals          ,ONLY: IERROR,MPI_COMM_NULL,MPI_DOUBLE_PRECISION,MPI_COMM_PICLAS,MPI_INFO_NULL,MPI_UNDEFINED,MPIRoot
USE MOD_Globals          ,ONLY: UNIT_StdOut
USE MOD_HDG_Vars         ,ONLY: CPPCOMM
USE MOD_Mesh_Vars        ,ONLY: nBCSides,BC
#endif /*USE_MPI*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, PARAMETER :: BCTypeCPP(1:3) = (/2,52,60/) ! BCType which allows coupled power potential control
!                                                  !  2: DC + coupled power for DC potential adjustment
!                                                  ! 52: AC with cos(wt) + coupled power for AC potential adjustment + bias voltage
!                                                  ! 60: AC with cos(wt) + coupled power for AC potential adjustment
INTEGER            :: iBC,CPPBoundaries,BCType,BCState
#if USE_MPI
LOGICAL            :: BConProc
INTEGER            :: color,SideID
#endif /*USE_MPI*/
!===================================================================================================================================

! 1.) Get global number of coupled power potential boundaries in [1:nBCs]
UseCoupledPowerPotential = .FALSE.
CPPBoundaries = 0
! Loop over global BC list
DO iBC=1,nBCs
  BCType = BoundaryType(iBC,BC_TYPE)
  IF(.NOT.ANY(BCType.EQ.BCTypeCPP)) CYCLE ! Skip non-CPP boundaries
  BCState = BoundaryType(iBC,BC_STATE) ! BCState of the corresponding BCType
  IF((BCType.EQ.2).AND.(BCState.NE.2)) CYCLE ! BCType=2 must be combined with BCState=2
  CPPBoundaries=CPPBoundaries+1
  IF(BCState.LE.0) CALL CollectiveStop(__STAMP__,' BCState for FPC must be >0! BCState=',IntInfo=BCState)
END DO

IF(CPPBoundaries.EQ.0) RETURN ! Already determined in HDG initialization

UseCoupledPowerPotential = .TRUE.

! 2.) Get parameters
#if USE_LOADBALANCE
! Do not set during load balance in order to keep the old value
IF(.NOT.PerformLoadBalance)THEN
#endif /*USE_LOADBALANCE*/
  ! Electric potential (/min, start, max/) Note that start is only required at t=0 and is used for the BC potential
  CoupledPowerPotential = 0.
  CoupledPowerPotential(1:3) = GETREALARRAY('CoupledPowerPotential',3)
  ! Get target power value
  CoupledPowerTarget = GETREAL('CoupledPowerTarget')
  IF(CoupledPowerTarget.LE.0.) CALL CollectiveStop(__STAMP__,'CoupledPowerTarget must be > 0.')
  ! Get relaxation factor
  CoupledPowerRelaxFac = GETREAL('CoupledPowerRelaxFac')
  IF(CoupledPowerRelaxFac.LE.0.) CALL CollectiveStop(__STAMP__,'CoupledPowerRelaxFac must be > 0.')
  ! Get frequency within which the integrated power is calculated
  CoupledPowerFrequency = GETREAL('CoupledPowerFrequency')
  IF(CoupledPowerFrequency.LT.0.) CALL CollectiveStop(__STAMP__,'CoupledPowerFrequency must be >= 0.')
  ! Sanity check
  CoupledPowerMode = GETINTFROMSTR('CoupledPowerMode')
  IF(.NOT.ANY(CoupledPowerMode.EQ.(/1,2,3/))) CALL CollectiveStop(__STAMP__,'CoupledPowerMode with unknown value!')
  IF((CoupledPowerMode.EQ.3).AND.(ABS(CoupledPowerFrequency).LE.0.))&
      CALL CollectiveStop(__STAMP__,'CoupledPowerMode=3 requires a positive value for CoupledPowerFrequency!')
  ! Update time
  IF(CoupledPowerFrequency.GT.0.0) CoupledPowerPotential(5) = 1.0 / CoupledPowerFrequency
#if USE_LOADBALANCE
END IF ! .NOT.PerformLoadBalance
#endif /*USE_LOADBALANCE*/

#if USE_MPI
! 3.) Establish sub-communicator (all BCs directly connected to a coupled power potential boundary); MPIRoot is always in the group
BConProc = .FALSE.
IF(MPIRoot)THEN
  BConProc = .TRUE.
ELSE
  CoupledPowerPotential(1:CPPDataLength) = 0.
  ! Check local sides
  DO SideID=1,nBCSides
    iBC    = BC(SideID)
    BCType = BoundaryType(iBC,BC_TYPE)
    IF(.NOT.ANY(BCType.EQ.BCTypeCPP)) CYCLE ! Skip other boundaries
    BCState = BoundaryType(iBC,BC_STATE) ! BCState of the corresponding BCType
    IF((BCType.EQ.2).AND.(BCState.NE.2)) CYCLE ! BCType=2 must be combined with BCState=2
    BConProc = .TRUE.
  END DO ! SideID=1,nBCSides
END IF ! MPIRoot

! create new communicator
color = MERGE(CPPBoundaries, MPI_UNDEFINED, BConProc)

! set communicator id
CPPCOMM%ID = CPPBoundaries

! create new emission communicator for coupled power potential communication. Pass MPI_INFO_NULL as rank to follow the original ordering
CALL MPI_COMM_SPLIT(MPI_COMM_PICLAS, color, MPI_INFO_NULL, CPPCOMM%UNICATOR, iError)

! Find my rank on the shared communicator, comm size and proc name
IF(BConProc)THEN
  CALL MPI_COMM_RANK(CPPCOMM%UNICATOR, CPPCOMM%MyRank, iError)
  CALL MPI_COMM_SIZE(CPPCOMM%UNICATOR, CPPCOMM%nProcs, iError)

  ! inform about size of emission communicator
  IF (CPPCOMM%MyRank.EQ.0) THEN
#if USE_LOADBALANCE
    IF(.NOT.PerformLoadBalance)&
#endif /*USE_LOADBALANCE*/
        WRITE(UNIT_StdOut,'(A,I0,A)') ' Coupled power potential boundary condition: Emission-Communicator on ',&
            CPPCOMM%nProcs,' procs'
  END IF
END IF ! BConProc
#endif /*USE_MPI*/

! 4.) When restarting, load the last known values of CoupledPowerPotential(1:3) from the .h5 state file
CALL ReadCPPDataFromH5()

END SUBROUTINE InitCoupledPowerPotential


!===================================================================================================================================
!> Read the coupled power potential (CPP) data from a .h5 state file.
!> 1. The MPI root process reads the info and checks data consistency
!> 2. The MPI root process distributes the information among the sub-communicator processes connected to the CPP boundary.
!===================================================================================================================================
SUBROUTINE ReadCPPDataFromH5()
! MODULES
USE MOD_io_hdf5
USE MOD_Globals          ,ONLY: UNIT_stdOut,MPIRoot,IK,abort
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
#endif /*USE_LOADBALANCE*/
USE MOD_IO_HDF5          ,ONLY: OpenDataFile,CloseDataFile,File_ID
USE MOD_Restart_Vars     ,ONLY: DoRestart,RestartFile
USE MOD_HDF5_Input       ,ONLY: DatasetExists,ReadArray,GetDataSize
USE MOD_HDG_Vars         ,ONLY: CoupledPowerPotential,CPPDataLength
#if USE_MPI
USE MOD_Equation_Tools   ,ONLY: SynchronizeCPP
#endif /*USE_MPI*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(255) :: ContainerName
LOGICAL        :: CPPExists
REAL           :: CPPDataHDF5(1:CPPDataLength)
!===================================================================================================================================
! Only required during restart
IF(.NOT.DoRestart) RETURN

#if USE_LOADBALANCE
! Do not try to read the data from .h5 if load balance is performed without creating a .h5 restart file
IF(PerformLoadBalance.AND..NOT.(UseH5IOLoadBalance)) RETURN
#endif /*USE_LOADBALANCE*/

! 1. The MPI root process reads the info and checks data consistency
! Only root reads the values and distributes them via MPI Broadcast
IF(MPIRoot)THEN
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
  ! Check old parameter name
  ContainerName='CoupledPowerPotential'
  CALL DatasetExists(File_ID,TRIM(ContainerName),CPPExists)
  ! Check for new parameter name
  IF(CPPExists)THEN
    CALL ReadArray(TRIM(ContainerName) , 2 , (/1_IK , INT(CPPDataLength,IK)/) , 0_IK , 1 , RealArray=CPPDataHDF5)
    WRITE(UNIT_stdOut,'(6(A,ES10.2E3))') " Read CoupledPowerPotential from restart file ["//TRIM(RestartFile)//&
        "] min[V]: "                 ,CPPDataHDF5(1),&
        ", current[V]: "             ,CPPDataHDF5(2),&
        ", max[V]: "                 ,CPPDataHDF5(3),&
        ", intPower[Ws]: "           ,CPPDataHDF5(4),&
        ", next adjustment time[s]: ",CPPDataHDF5(5),&
        ", last-energy[J]: "         ,CPPDataHDF5(6)
    CoupledPowerPotential = CPPDataHDF5
  END IF ! CPPExists
  CALL CloseDataFile()
END IF ! MPIRoot

#if USE_MPI
! 2. The MPI root process distributes the information among the sub-communicator processes for each EPC
CALL SynchronizeCPP()
#endif /*USE_MPI*/

END SUBROUTINE ReadCPPDataFromH5
#endif /*defined(PARTICLES)*/

SUBROUTINE ExactFunc(ExactFunction,x,resu,t,ElemID,iRefState,iLinState,BCState)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals         ,ONLY: Abort
USE MOD_Globals_Vars    ,ONLY: PI,ElementaryCharge,eps0
USE MOD_Equation_Vars   ,ONLY: IniCenter,IniHalfwidth,IniAmplitude,RefState,LinPhi,LinPhiHeight,LinPhiNormal,LinPhiBasePoint
#if defined(PARTICLES)
USE MOD_HDG_Vars        ,ONLY: CoupledPowerPotential,UseCoupledPowerPotential,BiasVoltage,UseBiasVoltage
USE MOD_Particle_Vars   ,ONLY: Species,nSpecies!,PartState,PDM
#endif /*defined(PARTICLES)*/
USE MOD_Dielectric_Vars ,ONLY: DielectricRatio,Dielectric_E_0,DielectricRadiusValue,DielectricEpsR
USE MOD_Mesh_Vars       ,ONLY: ElemBaryNGeo
USE MOD_HDG_Vars        ,ONLY: FPC,EPC
USE MOD_TimeDisc_Vars   ,ONLY: time
#if USE_MPI
USE MOD_Globals         ,ONLY: mpiroot
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: x(3)
INTEGER,INTENT(IN)              :: ExactFunction    ! determines the exact function
INTEGER,INTENT(IN),OPTIONAL     :: ElemID           ! ElemID
REAL,INTENT(IN),OPTIONAl        :: t                ! time
INTEGER,INTENT(IN),OPTIONAL     :: iRefState        ! i-th reference state
INTEGER,INTENT(IN),OPTIONAL     :: iLinState        ! i-th linear potential state
INTEGER,INTENT(IN),OPTIONAL     :: BCState          ! BCState
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Resu(1:PP_nVar)    ! state in conservative variables
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: Omega,r1,r2,r_2D,r_3D,r_bary,cos_theta,eps1,eps2,xi,a(3),b(3),Q
#if defined(PARTICLES)
INTEGER                         :: i!,iPart
#endif /*defined(PARTICLES)*/
!===================================================================================================================================
SELECT CASE (ExactFunction)
#if defined(PARTICLES)
CASE(-5) ! Bias voltage DC boundary
  Resu(:) = BiasVoltage%BVData(1)
#endif /*defined(PARTICLES)*/
CASE(-4) ! Electric potential condition (EPC) where charge accumulated over one time step is removed and creates a voltage
  Resu(:) = EPC%Voltage(EPC%Group(BCState,2))
CASE(-3) ! linear function with base point, normal vector and heigh: Requires BoundaryType = (/7,X/)
  ASSOCIATE( height => LinPhiHeight(iLinState)    ,&
             phi    => LinPhi(iLinState)          )
    a = x(1:3) - LinPhiBasePoint(1:3,iLinState)
    b = LinPhiNormal(1:3,iLinState)
    xi = DOT_PRODUCT(a,b)
    IF(xi.GT.0.)THEN
      IF(xi.GT.height)THEN
        Resu(:)=phi
      ELSE
        Resu(:)=phi*xi/height
      END IF ! x(3)
    ELSE
      Resu(:)=0.0
    END IF ! xi.GT.0
  END ASSOCIATE
CASE(-2) ! Signal without zero-crossing (always positive or negative), otherwise  like CASE(-1):
         ! Amplitude, Frequency and Phase Shift supplied by RefState
  ! RefState(1,iRefState): amplitude
  ! RefState(2,iRefState): frequency
  ! efState(3,iRefState): phase shift
  Omega   = 2.*PI*RefState(2,iRefState)
  r1      = RefState(1,iRefState) / 2.0
  Resu(:) = r1*(COS(Omega*t+RefState(3,iRefState)) + 1.0)
CASE(-1) ! Signal with zero-crossing: Amplitude, Frequency and Phase Shift supplied by RefState
  ! RefState(1,iRefState): amplitude
  ! RefState(2,iRefState): frequency
  ! RefState(3,iRefState): phase shift
  Omega   = 2.*PI*RefState(2,iRefState)
#if defined(PARTICLES)
  ! Match coupled power by adjusting the coefficient of the cos function
  IF(UseCoupledPowerPotential)THEN
    Resu(:) = CoupledPowerPotential(2)*COS(Omega*t+RefState(3,iRefState))
  ELSE
#endif /*defined(PARTICLES)*/
    Resu(:) = RefState(1,iRefState)*COS(Omega*t+RefState(3,iRefState))
#if defined(PARTICLES)
  END IF ! UseCoupledPowerPotential
  ! Add bias potential (only if bias voltage model is activated, BCType is 51 for DC or 52 for AC)
  IF(UseBiasVoltage) Resu(:) = Resu(:) + BiasVoltage%BVData(1)
#endif /*defined(PARTICLES)*/
CASE(0) ! constant 0.
    Resu(:)=0.
#if defined(PARTICLES)
CASE(2) ! Floating voltage that depends on the coupled power. User-specified power input as target value sets the potential.
  Resu(:) = CoupledPowerPotential(2)
#endif /*defined(PARTICLES)*/
CASE(1001) ! linear in y-z
    Resu(:)=x(2)*2340 + x(3)*2340
CASE(102) !linear: z=-1: 0, z=1, 1000
  resu(:)=(1+x(3))*1000.
CASE(103) ! dipole
  r1=SQRT(SUM((x(:)-(IniCenter(:)-(/IniHalfwidth,0.,0./)))**2)) !+1.0E-3
  r2=SQRT(SUM((x(:)-(IniCenter(:)+(/IniHalfwidth,0.,0./)))**2)) !+1.0E-3
  resu(:)=IniAmplitude*(1/r2-1/r1)
CASE(104) ! solution to Laplace's equation: Phi_xx + Phi_yy + Phi_zz = 0
  resu(1) = ( COS(x(1))+SIN(x(1)) )*( COS(x(2))+SIN(x(2)) )*( COSH(SQRT(2.0)*x(3))+SINH(SQRT(2.0)*x(3)) )
CASE(105) ! 3D periodic test case
  resu(1)= SIN(x(1) + 1) * SIN(x(2) + 2) * SIN(x(3) + 3)
CASE(200) ! Dielectric Sphere of Radius R in constant electric field E_0 from book:
  ! John David Jackson, Classical Electrodynamics, 3rd edition, New York: Wiley, 1999.
  ! E_0       : constant electric field in z-direction far away from sphere
  ! R         : constant radius of the sphere
  ! eps_outer : dielectric constant of surrounding medium
  ! eps_inner : dielectric constant of sphere
  ! DielectricRatio = eps_inner / eps_outer (set in dielectric init)

  ! set radius and angle for DOF position x(1:3)
  r_2D   = SQRT(x(1)**2+x(2)**2)
  r_3D   = SQRT(x(1)**2+x(2)**2+x(3)**2)
  IF(r_3D.EQ.0.0)THEN
    cos_theta = 0.0
  ELSE
    cos_theta = x(3) / r_3D
  END IF
  IF(PRESENT(ElemID))THEN ! if ElemID is present, use for bary center determination versus sphere radius
    r_bary = SQRT(ElemBaryNGeo(1,ElemID)**2+ElemBaryNGeo(2,ElemID)**2+ElemBaryNGeo(3,ElemID)**2)
  ELSE
    r_bary = r_3D
  END IF

  ! depending on the radius the solution for the potential is different for inner/outer parts of the domain
  IF(r_bary.LT.DielectricRadiusValue)THEN ! inside sphere: DOF and element bary center
    ! Phi_inner = - (3 / (2 + eps_inner / eps_outer)) * E_1 * z
    !resu(1:PP_nVar) = -(3./(DielectricRatio+2.))*Dielectric_E_0*x(3)
    resu(1:PP_nVar) = -(3./(DielectricRatio+2.))*r_3D*cos_theta*Dielectric_E_0
  ELSEIF(r_bary.GE.DielectricRadiusValue)THEN ! outside sphere
    ! Phi_outer = ( (eps_inner / eps_outer - 1 )/( eps_inner / eps_outer + 2 ) * ( R^3/r^3 )   - 1 ) * E_0 * z
    resu(1:PP_nVar) =  ( (DielectricRatio-1)        / (DielectricRatio+2)       ) *&
                       Dielectric_E_0*(DielectricRadiusValue**3/r_3D**2)*cos_theta-Dielectric_E_0 * r_3D*cos_theta

                       !( (DielectricRadiusValue**3) / (r_3D**3) ) - 1 )*(Dielectric_E_0 * x(3))
                       !( (DielectricRadiusValue**3) / ((r_2D**2+x(3)**2)**(3./2.)) ) - 1 )*(Dielectric_E_0 * x(3))
  ELSE
    IF(PRESENT(ElemID))THEN
      SWRITE(*,*) "ElemID                ",ElemID
      SWRITE(*,*) "ElemBaryNGeo(1:3)     ",ElemBaryNGeo(1,ElemID),ElemBaryNGeo(2,ElemID),ElemBaryNGeo(3,ElemID),r_bary
    END IF
    SWRITE(*,*) "x(1),x(2),x(3)        ",x(1),x(2),x(3)
    SWRITE(*,*) "DielectricRadiusValue ",DielectricRadiusValue
    CALL abort(__STAMP__,'Dielectric sphere. Invalid radius for exact function!')
  END IF

  ! varphi = ATAN2(x(2),x(1)) ! only needed for the electric field
  !   E_r,inner = 0
  !   E_z,inner = (3 / (2 + eps_inner / eps_outer)) * E_0
  !
  !   E_r,outer = 3 * ( (eps_inner / eps_outer - 1 )/( eps_inner / eps_outer + 2 ) * ( R^3/r^4 ) ) * E_0 * z
  !   E_z,inner =   ( - (eps_inner / eps_outer - 1 )/( eps_inner / eps_outer + 2 ) * ( R^3/r^3 )   + 1 ) * E_0
CASE(300) ! Dielectric Slab in z-direction of half width R in constant electric field E_0: adjusted from CASE(200)
  ! R = DielectricRadiusValue
  ! DielectricRatio = eps/eps0

  ! for BC, not ElemID will be given
  IF(.NOT.PRESENT(ElemID))THEN
    resu(1:PP_nVar) = -Dielectric_E_0*x(3)*(-DielectricRadiusValue   *((DielectricRatio-1.)/(DielectricRatio))/(abs(x(3))) + 1)
    RETURN
  END IF

  ! depending on the radius the solution for the potential is different for inner/outer parts of the domain
  IF(ABS(ElemBaryNGeo(3,ElemID)).LT.DielectricRadiusValue)THEN ! inside box: DOF and element bary center
    ! Phi_inner = ?

    ! marcel
    !resu(1:PP_nVar) = -Dielectric_E_0*x(3)*(1 - (DielectricRatio-1)/(DielectricRatio+2))
    resu(1:PP_nVar) = -Dielectric_E_0*x(3)*(-((DielectricRatio-1.)/(DielectricRatio)) + 1)

    ! linear
    !resu(1:PP_nVar) = -(1./(DielectricRatio))*Dielectric_E_0*x(3)

    ! from sphere
    ! Phi_inner = - (3 / (2 + eps_inner / eps_outer)) * E_1 * z
    !resu(1:PP_nVar) = -(3./(DielectricRatio+2.))*x(3)*Dielectric_E_0



  ELSEIF(ABS(ElemBaryNGeo(3,ElemID)).GT.DielectricRadiusValue)THEN ! outside sphere
    ! Phi_outer = ?
    !resu(1:PP_nVar) =( ( (DielectricRatio-1)        / (DielectricRatio+2)       ) *&
                       !( (DielectricRadiusValue**3) / (x(3)**(3.)) ) - 1 )*(Dielectric_E_0 * x(3))
                       !( (DielectricRadiusValue**3) / ((x(3)**2)**(3./2.)) ) - 1 )*(Dielectric_E_0 * x(3))

    ! marcel
    resu(1:PP_nVar) = -Dielectric_E_0*x(3)*(-DielectricRadiusValue**2*((DielectricRatio-1.)/(DielectricRatio))/(x(3)**2) + 1)

    ! linear
    resu(1:PP_nVar) = -Dielectric_E_0*x(3)*(-DielectricRadiusValue   *((DielectricRatio-1.)/(DielectricRatio))/(abs(x(3))) + 1)

    !resu(1:PP_nVar) = -Dielectric_E_0*(x(3) - sign(1.0,x(3))*((DielectricRatio-1)/(DielectricRatio+2))*((2**3)/(x(3)**2)) )
    !resu(1:PP_nVar) =( ( (DielectricRatio-1)        / (DielectricRatio+2)       ) *&
                       !( (2.0**3) / ((x(3)**3) ) - 1 )*(Dielectric_E_0 * x(3))
  ELSE
    SWRITE(*,*) "ElemID                ",ElemID
    SWRITE(*,*) "x(1),x(2),x(3)        ",x(1),x(2),x(3)
    SWRITE(*,*) "ElemBaryNGeo(1:3)     ",ElemBaryNGeo(1,ElemID),ElemBaryNGeo(2,ElemID),ElemBaryNGeo(3,ElemID)
    SWRITE(*,*) "DielectricRadiusValue ",DielectricRadiusValue
    CALL abort(__STAMP__,'Dielectric sphere. Invalid radius for exact function!')
  END IF
CASE(301) ! like CASE=300, but only in positive z-direction the dielectric region is assumed
  ! R = DielectricRadiusValue
  ! DielectricRatio = eps/eps0

  ! for BC, not ElemID will be given
  IF(.NOT.PRESENT(ElemID))THEN
    IF(x(3).GT.0.0)THEN ! inside dielectric
      resu(1:PP_nVar) = -(Dielectric_E_0/DielectricRatio)*(x(3)-DielectricRadiusValue)
    ELSE
      resu(1:PP_nVar) = -(Dielectric_E_0)*(x(3)-DielectricRadiusValue/DielectricRatio)
    END IF
    RETURN
  END IF

  ! depending on the radius the solution for the potential is different for inner/outer parts of the domain
  IF(     (ABS(ElemBaryNGeo(3,ElemID)).LT.2.0*DielectricRadiusValue).AND.(ElemBaryNGeo(3,ElemID).GT.0.0))THEN ! inside box: DOF and element bary center
    resu(1:PP_nVar) = -(Dielectric_E_0/DielectricRatio)*(x(3)-DielectricRadiusValue)
  ELSEIF( (ABS(ElemBaryNGeo(3,ElemID)).GT.DielectricRadiusValue).OR.(ElemBaryNGeo(3,ElemID).LT.0.0) )THEN ! outside sphere
    resu(1:PP_nVar) = -(Dielectric_E_0)*(x(3)-DielectricRadiusValue/DielectricRatio)
  ELSE
    SWRITE(*,*) "ElemID                ",ElemID
    SWRITE(*,*) "x(1),x(2),x(3)        ",x(1),x(2),x(3)
    SWRITE(*,*) "ElemBaryNGeo(1:3)     ",ElemBaryNGeo(1,ElemID),ElemBaryNGeo(2,ElemID),ElemBaryNGeo(3,ElemID)
    SWRITE(*,*) "DielectricRadiusValue ",DielectricRadiusValue
    CALL abort(__STAMP__,'Dielectric sphere. Invalid radius for exact function!')
  END IF

CASE(400,401) ! Point Source in Dielectric Region with epsR_1  = 1 for x < 0 (vacuum)
  !                                                epsR_2 != 1 for x > 0 (dielectric region)
  ! DielectricRadiusValue is used as distance between dielectric interface and position of charged point particle
  ! set radius and angle for DOF position x(1:3)
  ! Limitations:
  ! only valid for eps_2 = 1*eps0
  ! and q = 1
  r_2D   = SQRT(x(1)**2+x(2)**2)
  r1 = SQRT(r_2D**2 + (DielectricRadiusValue-x(3))**2)
  r2 = SQRT(r_2D**2 + (DielectricRadiusValue+x(3))**2)

  IF(ExactFunction.EQ.400)THEN
    ! Vacuum bottom, dielectric top
    eps1 = DielectricEpsR*eps0
    eps2 = 1.0*eps0
  ELSEIF(ExactFunction.EQ.401)THEN
    ! Dielectric bottom, vacuum top
    eps1 = 1.0*eps0
    eps2 = DielectricEpsR*eps0
  END IF ! ExactFunction.EQ.400
#if defined(PARTICLES)
  IF(ALLOCATED(Species).AND.(nSpecies.GT.0))THEN
    IF(ABS(Species(1)%ChargeIC).GT.0.)THEN
      Q = Species(1)%ChargeIC/(4.0*PI)
    END IF ! Species(1)%ChargeIC.GT.0
  ELSE
    Q = 0.
  END IF ! ALLOCATED(Species
#else
  Q = 1.0/(4.0*PI)
#endif /*defined(PARTICLES)*/
  ASSOCIATE( eps12 => eps1+eps2 )

    IF(x(3).GT.0.0)THEN
      IF(ALL((/ x(1).EQ.0.0,  x(2).EQ.0.0, x(3).EQ.DielectricRadiusValue /)))THEN
        print*, "HERE?!?!?!"
      END IF
      IF((r1.LE.0.0).OR.(r2.LE.0.0))THEN
        SWRITE(*,*) "r1=",r1
        SWRITE(*,*) "r2=",r2
        CALL abort(__STAMP__,'ExactFunc=400: Point source in dielectric region. Cannot evaluate the exact function at the singularity!')
      END IF
      resu(1:PP_nVar) = (Q/eps1)*(1./r1 - ((eps2-eps1)/eps12)*(1./r2) )
    ELSE
      IF(r1.LE.0.0)THEN
        SWRITE(*,*) "r1=",r1
        CALL abort(__STAMP__,'Point source in dielectric region: Cannot evaluate the exact function at the singularity!')
      END IF
      resu(1:PP_nVar) = (2.0*Q/eps12) * 1./r1
    END IF
  END ASSOCIATE
CASE(500) ! Coaxial capacitor with Floating Boundary Condition (FPC) with from
  ! Chen 2020 "A hybridizable discontinuous Galerkin method for simulation of electrostatic problems with floating potential conductors".
  r_2D = SQRT(x(1)**2+x(2)**2)
  Q = 0. ! Initialize
  IF(ALLOCATED(FPC%Charge)) Q = FPC%Charge(1) ! [C] - accumulated charge on iUniqueFPCBC = 1
  ASSOCIATE( &
        V0  => 0                  ,& ! [V]
        V1  => 10.0               ,& ! [V]
        r0  => 0.1e-2             ,& ! [m]
        r1  => 2.0e-2             ,& ! [m]
        r2  => 0.8e-2             ,& ! [m]
        r3  => 1.2e-2             ,& ! [m]
        !eps => ElementaryCharge    & ! [e]
        eps => eps0                &
        )
    ASSOCIATE( C20 => LOG(r2/r0) , C31 => LOG(r3/r1) )
      ASSOCIATE( b1 => (V1 - V0 - C20*Q/(2*PI*eps))/(C20-C31) )
        ASSOCIATE( b0 => b1+Q/(2*PI*eps) )
          ASSOCIATE( a0 => V0-b0*LOG(r0) , a1 => V1-b1*LOG(r1) )
            ! Check if point is located in first or second region
            IF(r_2D.LT.(r2+r3)/2.0)THEN
              resu = a0 + b0 * LOG(r_2D)
            ELSE
              resu = a1 + b1 * LOG(r_2D)
            END IF ! r.LT.(r2+r3)/2.0
          END ASSOCIATE
        END ASSOCIATE
      END ASSOCIATE
    END ASSOCIATE
  END ASSOCIATE
#if !(USE_PETSC)
  CALL abort(__STAMP__,'ExactFunc=500 requires PICLAS_PETSC=ON')
#endif /*!(USE_PETSC)*/
CASE(600) ! 2 cubes with two different charges
  IF(ALLOCATED(FPC%Charge))THEN
    FPC%Charge(1)=5.0
    FPC%Charge(2)=10.0
  END IF ! ALLOCATED(FPC%Charge)
  resu = 0.
#if !(USE_PETSC)
  CALL abort(__STAMP__,'ExactFunc=600 requires PICLAS_PETSC=ON')
#endif /*!(USE_PETSC)*/
CASE(700) ! Analytical solution of a charged particle moving in cylindrical coordinates between two grounded walls
#if defined(PARTICLES)
  eps1 = -ElementaryCharge/(4.0*PI*eps0)
  !eps1 = -ElementaryCharge/(PI*eps0)
  !eps1 = -ElementaryCharge/((4*PI*eps0)**(3./2.))
  !IF(ALLOCATED(Species))THEN
    !eps1 = -Species(1)%ChargeIC/(4.0*PI*eps0)
    !eps1 = Species(1)%ChargeIC/(PI*eps0)
  !ELSE
    !eps1=0
  !END IF ! ALLOCATED(Species)
  r1 = x(1)**2 + x(2)**2
  resu = 0.
  !DO iPart = 1, PDM%ParticleVecLength
    !IF(.NOT.PDM%ParticleInside(iPart)) CYCLE
    !ASSOCIATE( z => x(3), zq => PartState(3,1), H => 80e-3, iMax => 20 )
    ASSOCIATE( z => x(3), zq => 70e-3-5e6*time, H => 80e-3, iMax => 20 )
      resu = resu + 1./SQRT(r1+(z-zq)**2) - 1./SQRT(r1+(z+zq)**2)
      DO i = 1, iMax
        resu = resu + 1./SQRT(r1+(z+REAL(2*i)*H-zq)**2) + 1./SQRT(r1+(z-REAL(2*i)*H-zq)**2)
      END DO ! i = 1, iMax
      DO i = 1, iMax
        resu = resu - 1./SQRT(r1+(z+REAL(2*i)*H+zq)**2) - 1./SQRT(r1+(z-REAL(2*i)*H+zq)**2)
      END DO ! i = 1, iMax
    END ASSOCIATE
    resu = eps1 * resu
  !END DO ! iPart = 1, PDM%ParticleVecLength
#else
  CALL abort(__STAMP__,'ExactFunc=700 requires PARTICLES=ON')
#endif /*defined(PARTICLES)*/
CASE(800,801,900) ! Dielectric slab on electrode (left) with plasma between slab and other electrode opposite
  resu = 0.
  ASSOCIATE( x     => x(1)   , &
             y     => x(2)   , &
             z     => x(3)   , &
             L     =>  1e-3  , &
             d     =>  1e-8  , &
             eps1  =>  1e1   , &
             sigma =>  1e-2  , &
             rho0  => -1e-4  , &
             Phi0  => -1e1    )
    ASSOCIATE( PhiF => ((d/L)/((d/L)+eps1))*( L*(sigma + 0.5*rho0*L)/eps0 + Phi0 ) )
      ASSOCIATE( a => 0.5*rho0*L/eps0 + (Phi0 - PhiF)/L ,&
                 b => PhiF)
        IF(x.GE.0.0)THEN
          resu = -0.5*rho0*x**2/eps0 + a*x + b
        ELSE
          resu = b * (x/d + 1.0)
        END IF ! x.GE.0.0
      END ASSOCIATE
    END ASSOCIATE
  END ASSOCIATE
CASE DEFAULT
  CALL abort(__STAMP__,'Exactfunction not specified!', IntInfoOpt=ExactFunction)
END SELECT ! ExactFunction

END SUBROUTINE ExactFunc



SUBROUTINE CalcSource(Ut)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Globals,ONLY:Abort
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:IniExactFunc
USE MOD_Equation_Vars,ONLY:IniCenter,IniHalfwidth,IniAmplitude
USE MOD_Mesh_Vars,ONLY:Elem_xGP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)              :: Ut(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iElem
REAL                             :: r1,r2
REAL,DIMENSION(3)                :: dx1,dx2,dr1dx,dr2dx,dr1dx2,dr2dx2
!===================================================================================================================================
SELECT CASE (IniExactFunc)
CASE(0) ! Particles
#ifdef PARTICLES
!  DO iElem=1,PP_nElems
!    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
!      !  Get source from Particles
!      Ut(1:3,i,j,k,iElem) = Ut(1:3,i,j,k,iElem) - eps0inv * PartSource(1:3,i,j,k,iElem)
!      Ut(  8,i,j,k,iElem) = Ut(  8,i,j,k,iElem) + eps0inv * PartSource(  4,i,j,k,iElem) * c_corr
!    END DO; END DO; END DO
!  END DO
#endif /*PARTICLES*/
CASE(103)
DO iElem=1,PP_nElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
         dx1=(Elem_xGP(:,i,j,k,iElem)-(IniCenter(:)-(/IniHalfwidth,0.,0./)))
         dx2=(Elem_xGP(:,i,j,k,iElem)-(IniCenter(:)+(/IniHalfwidth,0.,0./)))
         r1=SQRT(SUM(dx1**2))
         r2=SQRT(SUM(dx2**2))
         dr1dx(:)= r1*dx1
         dr2dx(:)= r2*dx2
         dr1dx2(:)= r1+dr1dx(:)*dx1
         dr2dx2(:)= r2+dr2dx(:)*dx2
         Ut(:,i,j,k,iElem)=Ut(:,i,j,k,iElem)- IniAmplitude*( SUM((r1*dr1dx2(:)-2*dr1dx(:)**2)/(r1*r1*r1)) &
                                 -SUM((r2*dr2dx2(:)-2*dr2dx(:)**2)/(r2*r2*r2)) )
      END DO !i
    END DO !j
  END DO !k
END DO ! iElem=1,nElems

CASE DEFAULT
!  CALL abort(__STAMP__,&
             !'Exactfunction not specified!',999,999.)
END SELECT ! ExactFunction
END SUBROUTINE CalcSource


SUBROUTINE DivCleaningDamping()
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,       ONLY : U
USE MOD_Equation_Vars, ONLY : fDamping,DoParabolicDamping
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iElem
!===================================================================================================================================
IF(DoParabolicDamping) RETURN
DO iElem=1,PP_nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    !  Get source from Particles
    U(7:8,i,j,k,iElem) = U(7:8,i,j,k,iElem) * fDamping
  END DO; END DO; END DO
END DO
END SUBROUTINE DivCleaningDamping


#if defined(CODE_ANALYZE)
      SUBROUTINE CalcSourceHDG(i,j,k,iElem,resu, Phi, warning_linear, warning_linear_phi)
#else
PPURE SUBROUTINE CalcSourceHDG(i,j,k,iElem,resu, Phi, warning_linear, warning_linear_phi)
#endif /*defined(CODE_ANALYZE)*/
!===================================================================================================================================
! Determine the right-hand-side of Poisson's equation (either by an analytic function or deposition of charge from particles)
! TODO: currently particles are enforced, which means that they over-write the exact function solution because
! the combination of both has not been specified
! How should this function work???
! for dielectric regions DO NOT apply the scaling factor Eps_R here (which is const. in HDG due to current implementation) because
! it is in the tensor "chitens"
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars          ,ONLY: Elem_xGP
USE MOD_Globals_Vars       ,ONLY: eps0
#ifdef PARTICLES
USE MOD_PICDepo_Vars       ,ONLY: PartSource,DoDeposition
USE MOD_HDG_Vars           ,ONLY: ElemToBRRegion,UseBRElectronFluid,RegionElectronRef
#if IMPA
USE MOD_LinearSolver_Vars  ,ONLY: ExplicitPartSource
#endif
#if defined(CODE_ANALYZE)
USE MOD_Mesh_Vars          ,ONLY: offSetElem
USE MOD_Particle_Mesh_Vars ,ONLY: BoundsOfElem_Shared
#endif /*defined(CODE_ANALYZE)*/
#endif /*PARTICLES*/
USE MOD_Equation_Vars      ,ONLY: IniExactFunc
USE MOD_Equation_Vars      ,ONLY: IniCenter,IniHalfwidth,IniAmplitude
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: i, j, k,iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Resu(PP_nVar)    ! state in conservative variables
LOGICAL,INTENT(OUT),OPTIONAL    :: warning_linear
REAL,INTENT(OUT),OPTIONAL       :: warning_linear_phi
REAL,INTENT(IN),OPTIONAL        :: Phi
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: xvec(3)
REAL                            :: r1,r2
REAL,DIMENSION(3)               :: dx1,dx2,dr1dx,dr2dx,dr1dx2,dr2dx2
#ifdef PARTICLES
REAL                            :: source_e ! Electron charge density for Boltzmann relation (electrons as isothermal fluid!)
INTEGER                         :: RegionID
#if defined(CODE_ANALYZE)
REAL                            :: ElemCharLengthX
#endif /*defined(CODE_ANALYZE)*/
#endif /*PARTICLES*/
!===================================================================================================================================
ASSOCIATE( x => Elem_xGP(1,i,j,k,iElem), y => Elem_xGP(2,i,j,k,iElem), z => Elem_xGP(3,i,j,k,iElem))
IF(PRESENT(warning_linear)) warning_linear=.FALSE. ! Initialize
IF(PRESENT(warning_linear_phi)) warning_linear_phi=0. ! Initialize
! Calculate IniExactFunc before particles are superimposed, because the IniExactFunc might be needed by the CalcError function
SELECT CASE (IniExactFunc)
CASE(0) ! Particles
  resu=0. ! empty
CASE(103)
  dx1=(Elem_xGP(1:3,i,j,k,iElem)-(IniCenter(:)-(/IniHalfwidth,0.,0./)))
  dx2=(Elem_xGP(1:3,i,j,k,iElem)-(IniCenter(:)+(/IniHalfwidth,0.,0./)))
  r1=SQRT(SUM(dx1**2))
  r2=SQRT(SUM(dx2**2))
  dr1dx(:)= r1*dx1
  dr2dx(:)= r2*dx2
  dr1dx2(:)= r1+dr1dx(:)*dx1
  dr2dx2(:)= r2+dr2dx(:)*dx2
  resu(1)=- IniAmplitude*( SUM((r1*dr1dx2(:)-2*dr1dx(:)**2)/(r1*r1*r1)) - SUM((r2*dr2dx2(:)-2*dr2dx(:)**2)/(r2*r2*r2)) )
CASE(105) ! 3D periodic test case
  xvec(1:3) = Elem_xGP(1:3,i,j,k,iElem)
  resu(1)=-3 * SIN(xvec(1) + 1) * SIN(xvec(2) + 2) * SIN(xvec(3) + 3)
CASE DEFAULT
  resu=0.
END SELECT ! ExactFunction

#if defined(PARTICLES)
#if defined(CODE_ANALYZE)
! Specific source terms after particle deposition
SELECT CASE(IniExactFunc)
CASE(800) ! plasma between electrodes + particles
  ! Check Dirichlet elements
  IF(ElemHasDirichletBC(iElem))THEN
    ! Get length on element in 1D
    ASSOCIATE( Bounds => BoundsOfElem_Shared(1:2,1:3,iElem + offsetElem) ) ! 1-2: Min, Max value; 1-3: x,y,z
      ElemCharLengthX = ABS(Bounds(2,1)-Bounds(1,1)) ! ABS(max - min)
      ! Add scaled value in BC elements
      IF((x.GT.0.0).AND.(x.LT.1.0e-3))THEN
        IF(x.GT.0.5e-3)THEN
          ! Negative gradient
          PartSource(4,i,j,k,iElem) = PartSource(4,i,j,k,iElem) - 1e-4*(1e-3-x)/ElemCharLengthX
        ELSE
          ! Positive gradient
          PartSource(4,i,j,k,iElem) = PartSource(4,i,j,k,iElem) - 1e-4*x/ElemCharLengthX
        END IF ! x.GT.0.5e-3
      END IF ! (x.GT.0.0).AND.(x.LT.1.0e-3)
      !WRITE (*,*) "x,source =", x,PartSource(4,i,j,k,iElem),NINT(x*1e9), " nm", " TRUE"
    END ASSOCIATE
  ELSE
    ! Add constant value
    IF((x.GT.0.0).AND.(x.LT.1.0e-3)) PartSource(4,i,j,k,iElem) = PartSource(4,i,j,k,iElem) - 1e-4
  !WRITE (*,*) "x,source =", x,PartSource(4,i,j,k,iElem),NINT(x*1e9), " nm"
  END IF ! ElemHasDirichletBC(iElem)
CASE(801) ! plasma between electrodes + particles: Linear source
  IF(x.GT.0.0)THEN
    PartSource(4,i,j,k,iElem) = PartSource(4,i,j,k,iElem) - 1e-4*(1.0 - Elem_xGP(2,i,j,k,iElem)/1e-3)
  END IF ! x.GT.0.0
END SELECT
#endif /*defined(CODE_ANALYZE)*/

IF(DoDeposition)THEN
  source_e=0.
  IF(UseBRElectronFluid.AND.PRESENT(Phi))THEN
    RegionID=ElemToBRRegion(iElem)
    IF (RegionID .NE. 0) THEN
      source_e = Phi-RegionElectronRef(2,RegionID)
      IF (source_e .LT. 0.) THEN
        source_e = RegionElectronRef(1,RegionID) &         !--- Boltzmann relation (electrons as isothermal fluid!)
            * EXP( (source_e) / RegionElectronRef(3,RegionID) )
      ELSE
        ! Store delta for output
        IF(PRESENT(warning_linear_phi)) warning_linear_phi = MAX(warning_linear_phi,source_e)
        ! Linear approximation from Taylor expansion O(1)
        source_e = RegionElectronRef(1,RegionID) &         !--- linearized boltzmann relation at positive exponent
            * (1. + ((source_e) / RegionElectronRef(3,RegionID)) )
        IF(PRESENT(warning_linear)) warning_linear = .TRUE.
      END IF
      !source_e = RegionElectronRef(1,RegionID) &         !--- boltzmann relation (electrons as isothermal fluid!)
      !* EXP( (Phi-RegionElectronRef(2,RegionID)) / RegionElectronRef(3,RegionID) )
    END IF
  END IF ! UseBRElectronFluid
#if IMPA
  resu(1)= - (PartSource(4,i,j,k,iElem)+ExplicitPartSource(4,i,j,k,iElem)-source_e)/eps0
#else
  resu(1)= - (PartSource(4,i,j,k,iElem)-source_e)/eps0
#endif
END IF
#endif /*defined(PARTICLES)*/
END ASSOCIATE

END SUBROUTINE CalcSourceHDG


#if defined(PARTICLES) && defined(CODE_ANALYZE)
!===================================================================================================================================
!> Check if elements has at least one side that is a Dirichlet BC
!===================================================================================================================================
LOGICAL FUNCTION ElemHasDirichletBC(iElem) RESULT(L)
! MODULES
USE MOD_Particle_Mesh_Tools ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Mesh_Vars  ,ONLY: SideInfo_Shared
USE MOD_Mesh_Vars           ,ONLY: offSetElem
USE MOD_Mesh_Vars           ,ONLY: BoundaryType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: BCType,BCIndex,GlobalSideID,iLocSide
!===================================================================================================================================
L = .FALSE.
! Check 6 local sides for Dirichlet BC
DO iLocSide = 1, 6
  GlobalSideID = GetGlobalNonUniqueSideID(iElem+offSetElem,iLocSide)
  BCIndex = SideInfo_Shared(SIDE_BCID,GlobalSideID)
  ! Only check BC sides with BC index > 0
  IF(BCIndex.GT.0)THEN
    ! Get boundary type
    BCType = BoundaryType(BCIndex,BC_TYPE)
    ! Check if Dirichlet BC has been found
    SELECT CASE(BCType)
    CASE(HDGDIRICHLETBCSIDEIDS) ! Dirichlet
      L = .TRUE.
      RETURN
    END SELECT ! BCType
  END IF ! BCIndex.GT.0
END DO ! iLocSide = 1, 6
END FUNCTION ElemHasDirichletBC
#endif /*defined(PARTICLES) && defined(CODE_ANALYZE)*/


FUNCTION shapefunc(r)
!===================================================================================================================================
! Implementation of (possibly several different) shapefunctions
!===================================================================================================================================
! MODULES
  USE MOD_Equation_Vars, ONLY : shapeFuncPrefix, alpha_shape, rCutoff
! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
    REAL                 :: r         ! radius / distance to center
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
    REAL                 :: shapefunc ! sort of a weight for the source
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
   IF (r.GE.rCutoff) THEN
     shapefunc = 0.0
   ELSE
     shapefunc = ShapeFuncPrefix *(1-(r/rCutoff)**2)**alpha_shape
   END IF
END FUNCTION shapefunc

FUNCTION beta(z,w)
   IMPLICIT NONE
   REAL beta, w, z
   beta = GAMMA(z)*GAMMA(w)/GAMMA(z+w)           ! n   - kind=8
END FUNCTION beta

SUBROUTINE FinalizeEquation()
!===================================================================================================================================
! Deallocate the vars !!!!
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
EquationInitIsDone = .FALSE.
SDEALLOCATE(chitens)
SDEALLOCATE(chitensInv)
SDEALLOCATE(chitens_face)
SDEALLOCATE(E)
SDEALLOCATE(Et)
SDEALLOCATE(RefState)
SDEALLOCATE(LinPhiBasePoint)
SDEALLOCATE(LinPhiNormal)
SDEALLOCATE(LinPhiHeight)
SDEALLOCATE(LinPhi)
END SUBROUTINE FinalizeEquation

END MODULE MOD_Equation

