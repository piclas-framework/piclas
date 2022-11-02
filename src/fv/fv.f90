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

MODULE MOD_FV
!===================================================================================================================================
! Contains the initialization of the DG global variables
! Computes the different DG spatial operators/residuals(Ut) using U
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
INTERFACE FillIni
  MODULE PROCEDURE FillIni
END INTERFACE

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitFV
  MODULE PROCEDURE InitFV
END INTERFACE

INTERFACE FV_main
  MODULE PROCEDURE FV_main
END INTERFACE

INTERFACE FinalizeFV
  MODULE PROCEDURE FinalizeFV
END INTERFACE

PUBLIC::InitFV,FinalizeFV,FV_main,DefineParametersFV
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersFV()
! MODULES
USE MOD_ReadInTools       ,ONLY: prms
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("Finite Volumes")

CALL prms%CreateLogicalOption('FV-Reconstruction' , 'Use reconstruction for finite volumes', '.TRUE.')
CALL prms%CreateIntOption('FV-LimiterType',"Type of slope limiter of second order reconstruction", '1')

END SUBROUTINE DefineParametersFV


SUBROUTINE InitFV()
!===================================================================================================================================
! Allocate global variable U (solution) and Ut (dg time derivative).
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ReadInTools        ,ONLY: GETLOGICAL, GETINT
USE MOD_FV_Vars
USE MOD_FV_Metrics         ,ONLY: InitFV_Metrics
USE MOD_FV_Limiter         ,ONLY: InitFV_Limiter
USE MOD_Restart_Vars       ,ONLY: DoRestart,RestartInitIsDone
USE MOD_Interpolation_Vars ,ONLY: InterpolationInitIsDone
USE MOD_Mesh_Vars          ,ONLY: nSides, Face_xGP, NormVec
USE MOD_Mesh_Vars          ,ONLY: MeshInitIsDone
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#if !(USE_HDG)
USE MOD_LoadBalance_Vars   ,ONLY: UseH5IOLoadBalance
#endif /*!(USE_HDG)*/
#endif /*USE_LOADBALANCE*/
#if USE_MPI
USE MOD_MPI_Vars
USE MOD_MPI               ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
REAL                      :: Geotemp(6,1:nSides), dx_slave_temp(PP_nVar+1,1:nSides)
#endif /*USE_MPI*/
!===================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.(.NOT.MeshInitIsDone).OR.(.NOT.RestartInitIsDone).OR.FVInitIsDone) CALL abort(__STAMP__,&
    'InitFV not ready to be called or already called.')
LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT FV...'

doFVReconstruction=GETLOGICAL('FV-Reconstruction')

#if USE_LOADBALANCE && !(USE_HDG)
IF (.NOT.(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance))) THEN
#endif /*USE_LOADBALANCE && !(USE_HDG)*/
  ! the local DG solution in physical and reference space
  ALLOCATE( U(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
  U=0.
#if !(USE_HDG)
#if USE_LOADBALANCE
END IF
#endif /*USE_LOADBALANCE*/
! the time derivative computed with the DG scheme
ALLOCATE(Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
Ut=0.
#endif /*USE_HDG*/

nTotal_face=1!(PP_N+1)*(PP_N+1)
nTotal_vol=nTotal_face*1!(PP_N+1)
nTotalU=PP_nVar*nTotal_vol*PP_nElems

! U is filled with the ini solution
IF(.NOT.DoRestart) CALL FillIni()

! We store the interior data at the each element face
!ALLOCATE(U_Minus(PP_nVar,0:PP_N,0:PP_N,sideID_minus_lower:sideID_minus_upper))
!ALLOCATE(U_Plus(PP_nVar,0:PP_N,0:PP_N,sideID_plus_lower:sideID_plus_upper))
!U_Minus=0.
!U_Plus=0.

ALLOCATE(U_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(U_slave(PP_nVar,0:PP_N,0:PP_N,1:nSides))
U_master=0.
U_slave=0.

! ! Allocate arrays to hold the face flux to reduce memory churn
! ALLOCATE(U_Master_loc(1:PP_nVar        ,0:PP_N,0:PP_N))
! ALLOCATE(U_Slave_loc (1:PP_nVar        ,0:PP_N,0:PP_N))

! unique flux per side
! additional fluxes for the CFS-PML auxiliary variables (no PML: PMLnVar=0)
! additional fluxes for the CFS-PML auxiliary variables (no PML: PMLnVar=0)
ALLOCATE(Flux_Master(1:PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(Flux_Slave (1:PP_nVar,0:PP_N,0:PP_N,1:nSides))
! ALLOCATE(Flux_loc   (1:PP_nVar,0:PP_N,0:PP_N))
Flux_Master=0.
Flux_Slave=0.

IF (doFVReconstruction) THEN
  LimiterType = GETINT('FV-LimiterType')
  CALL InitFV_Limiter()

  ALLOCATE(FV_dx_slave(1:nSides))
  ALLOCATE(FV_dx_master(1:nSides))
#if (PP_TimeDiscMethod==600) /*DVM*/
  ALLOCATE(DVMtraj_slave(PP_nVar,1:nSides))
  ALLOCATE(DVMtraj_master(PP_nVar,1:nSides))
#endif

  ALLOCATE(FV_gradU(1:PP_nVar,1:nSides))

  !calculate face to center distances for reconstruction
#if USE_MPI
  !send face coordinates because MPI slave sides don't have them
  Geotemp=0.
  Geotemp(1:3,:)=Face_xGP(:,0,0,1:nSides)
  Geotemp(4:6,:)=NormVec(:,0,0,1:nSides)
  CALL StartReceiveMPIData(6,Geotemp,1,nSides,RecRequest_Geo,SendID=1) ! Receive YOUR
  CALL StartSendMPIData(6,Geotemp,1,nSides,SendRequest_Geo,SendID=1) ! Send MINE
  CALL FinishExchangeMPIData(SendRequest_Geo,RecRequest_Geo,SendID=1)
  Face_xGP(:,0,0,1:nSides)=Geotemp(1:3,:)
  NormVec(:,0,0,1:nSides)=Geotemp(4:6,:)

  ! distances for MPI sides - send direction
  CALL StartReceiveMPIData(PP_nVar+1,dx_slave_temp,1,nSides,RecRequest_U,SendID=2) ! Receive MINE
  CALL InitFV_Metrics(dx_slave_temp,doMPISides=.TRUE.)
  CALL StartSendMPIData(PP_nVar+1,dx_slave_temp,1,nSides,SendRequest_U,SendID=2) ! Send YOUR
#endif /*USE_MPI*/
  ! distances for BCSides, InnerSides and MPI sides - receive direction
  CALL InitFV_Metrics(dx_slave_temp,doMPISides=.FALSE.)
#if USE_MPI
  CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=2)
#endif /*USE_MPI*/
  FV_dx_slave(:)=dx_slave_temp(PP_nVar+1,:)
#if (PP_TimeDiscMethod==600)
  DVMtraj_slave(1:PP_nVar,:)=dx_slave_temp(1:PP_nVar,:)
#endif
END IF

FVInitIsDone=.TRUE.
LBWRITE(UNIT_stdOut,'(A)')' INIT DG DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitFV


SUBROUTINE FV_main(t,tStage,doSource)
!===================================================================================================================================
! Computes the FV time derivative
! U and Ut are allocated
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Vector
USE MOD_FV_Vars           ,ONLY: U,Ut,Flux_Master,Flux_Slave, doFVReconstruction
USE MOD_FV_Vars           ,ONLY: U_master, U_slave, FV_gradU
USE MOD_SurfInt            ,ONLY: SurfInt
USE MOD_ProlongToFace     ,ONLY: ProlongToFace, CalcFVGradients!, LimitFVGradients
USE MOD_FillFlux          ,ONLY: FillFlux
USE MOD_Equation          ,ONLY: CalcSource
USE MOD_Interpolation     ,ONLY: ApplyJacobian
! USE MOD_FillMortar        ,ONLY: U_Mortar,Flux_Mortar
#if USE_MPI
USE MOD_Mesh_Vars         ,ONLY: nSides
USE MOD_MPI_Vars
USE MOD_MPI               ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#if defined(PARTICLES) && defined(LSERK)
USE MOD_Particle_Vars     ,ONLY: DelayTime
USE MOD_TimeDisc_Vars     ,ONLY: time
#endif /*defined(PARTICLES) && defined(LSERK)*/
#ifdef PARTICLES
USE MOD_Particle_MPI      ,ONLY: MPIParticleSend,MPIParticleRecv
#endif /*PARTICLES*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers,ONLY: LBStartTime,LBPauseTime,LBSplitTime
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: t,tStage
LOGICAL,INTENT(IN)              :: doSource
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_LOADBALANCE
REAL                            :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================

! prolong the solution to the face integration points for flux computation
#if USE_MPI
! Prolong to face for MPI sides - send direction
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
CALL StartReceiveMPIData(PP_nVar,U_slave,1,nSides,RecRequest_U,SendID=2) ! Receive MINE
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DGCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
CALL ProlongToFace(U,U_master,U_slave,doMPISides=.TRUE.)
! CALL U_Mortar(U_master,U_slave,doMPISides=.TRUE.) !mortars not yet implemented
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/
CALL StartSendMPIData(PP_nVar,U_slave,1,nSides,SendRequest_U,SendID=2) ! Send YOUR
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DGCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/

! Prolong to face for BCSides, InnerSides and MPI sides - receive direction
CALL ProlongToFace(U,U_master,U_slave,doMPISides=.FALSE.)
! CALL U_Mortar(U_master,U_slave,doMPISides=.FALSE.)

#if USE_MPI
#if defined(PARTICLES) && defined(LSERK)
IF (time.GE.DelayTime) THEN
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_INTERPOLATION,tLBStart)
#endif /*USE_LOADBALANCE*/
#if USE_MPI
  CALL MPIParticleSend()
#endif /*USE_MPI*/
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
END IF
#endif /*defined(PARTICLES) && defined(LSERK)*/
#endif /*USE_MPI*/

!Nullify time derivative array
Ut=0.

IF (doFVReconstruction) THEN

#if USE_MPI
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/
  ! Complete send / receive of prolongtoface results
  CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=2)

  ! Gradient calculation for reconstruction
  CALL StartReceiveMPIData(PP_nVar,FV_gradU,1,nSides,RecRequest_gradUx,SendID=1) ! Receive YOUR
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_DGCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
  ! fill the global gradient list
  CALL CalcFVGradients(doMPISides=.TRUE.)
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/

  CALL StartSendMPIData(PP_nVar,FV_gradU,1,nSides,SendRequest_gradUx,SendID=1) ! Send MINE
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_DGCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/

  ! fill the all gradients on this proc
  CALL CalcFVGradients(doMPISides=.FALSE.)

! prolong the solution to the faces by applying reconstruction gradients
#if USE_MPI
! Prolong to face for MPI sides - send direction
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  ! Complete send / receive of gradients
  CALL FinishExchangeMPIData(SendRequest_gradUx,RecRequest_gradUx,SendID=1)

  ! prolong the solution to the faces by applying reconstruction gradients
  CALL StartReceiveMPIData(PP_nVar,U_slave,1,nSides,RecRequest_U,SendID=2) ! Receive MINE
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_DGCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
  CALL ProlongToFace(U,U_master,U_slave,FV_gradU,doMPISides=.TRUE.)
! CALL U_Mortar(U_master,U_slave,doMPISides=.TRUE.)
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/
  CALL StartSendMPIData(PP_nVar,U_slave,1,nSides,SendRequest_U,SendID=2) ! Send YOUR
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_DGCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/

  ! Prolong to face for BCSides, InnerSides and MPI sides - receive direction
  CALL ProlongToFace(U,U_master,U_slave,FV_gradU,doMPISides=.FALSE.)

END IF ! end of reconstruction


#if USE_MPI
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/
! Complete send / receive of prolongtoface results
CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=2)

! Initialization of the time derivative
!Flux=0. !don't nullify the fluxes if not really needed (very expensive)
CALL StartReceiveMPIData(PP_nVar,Flux_Slave,1,nSides,RecRequest_Flux,SendID=1) ! Receive YOUR
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DGCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
! fill the global surface flux list
CALL FillFlux(t,Flux_Master,Flux_Slave,U_master,U_slave,doMPISides=.TRUE.)
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/

CALL StartSendMPIData(PP_nVar,Flux_Slave,1,nSides,SendRequest_Flux,SendID=1) ! Send MINE
!CALL StartExchangeMPIData(PP_nVar,Flux,1,nSides,SendRequest_Flux,RecRequest_Flux,SendID=1) ! Send MINE - receive YOUR
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DGCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/

! fill the all surface fluxes on this proc
CALL FillFlux(t,Flux_Master,Flux_Slave,U_master,U_slave,doMPISides=.FALSE.)
! CALL Flux_Mortar(Flux_Master,Flux_Slave,doMPISides=.FALSE.)
! compute surface integral contribution and add to ut
CALL SurfInt(Flux_Master,Flux_Slave,Ut,doMPISides=.FALSE.)

#if USE_MPI
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/
! Complete send / receive
CALL FinishExchangeMPIData(SendRequest_Flux,RecRequest_Flux,SendID=1)
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DGCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/

!FINALIZE Fluxes for MPI Sides
! CALL Flux_Mortar(Flux_Master,Flux_Slave,doMPISides=.TRUE.)
CALL SurfInt(Flux_Master,Flux_Slave,Ut,doMPISides=.TRUE.)
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif

! Swap to right sign and divide by element volume (1/sJ)
CALL ApplyJacobian(Ut,toPhysical=.TRUE.,toSwap=.TRUE.)

! Add Source Terms
IF(doSource) CALL CalcSource(tStage,1.0,Ut)

#if USE_LOADBALANCE
CALL LBSplitTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/

#if defined(PARTICLES) && defined(LSERK)
#if USE_MPI
IF (time.GE.DelayTime) THEN
  CALL MPIParticleRecv()
END IF
#if USE_LOADBALANCE
CALL LBSplitTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/
#endif /*defined(PARTICLES) && defined(LSERK)*/

END SUBROUTINE FV_main


SUBROUTINE FillIni()
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_FV_Vars,ONLY:U
USE MOD_Mesh_Vars,ONLY:Elem_xGP
USE MOD_Equation_Vars,ONLY:IniExactFunc
USE MOD_Equation,ONLY:ExactFunc
#ifdef maxwell
USE MOD_Equation_Vars,ONLY:DoExactFlux
#endif /*maxwell*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iElem
!===================================================================================================================================
! Determine Size of the Loops, i.e. the number of grid cells in the
! corresponding directions
#ifdef maxwell
IF(DoExactFlux.AND.(IniExactFunc.NE.16)) RETURN ! IniExactFunc=16 is pulsed laser mixed IC+BC
#endif /*maxwell*/
DO iElem=1,PP_nElems
    CALL ExactFunc(IniExactFunc,0.,0,Elem_xGP(1:3,0,0,0,iElem),U(1:PP_nVar,0,0,0,iElem))
END DO ! iElem=1,PP_nElems
END SUBROUTINE FillIni



SUBROUTINE FinalizeFV()
!===================================================================================================================================
! Deallocate global variable U (solution) and Ut (dg time derivative).
!===================================================================================================================================
! MODULES
USE MOD_FV_Vars
#if USE_LOADBALANCE && !(USE_HDG)
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
#endif /*USE_LOADBALANCE && !(USE_HDG)*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(Ut)
SDEALLOCATE(U_master)
SDEALLOCATE(U_slave)
SDEALLOCATE(FLUX_Master)
SDEALLOCATE(FLUX_Slave)
! SDEALLOCATE(U_Master_loc)
! SDEALLOCATE(U_Slave_loc)
! SDEALLOCATE(Flux_loc)
IF (doFVReconstruction) THEN
  SDEALLOCATE(FV_dx_slave)
  SDEALLOCATE(FV_dx_master)
END IF

! Do not deallocate the solution vector during load balance here as it needs to be communicated between the processors
#if USE_LOADBALANCE && !(USE_HDG)
IF(.NOT.(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)))THEN
#endif /*USE_LOADBALANCE && !(USE_HDG)*/
  SDEALLOCATE(U)
#if USE_LOADBALANCE && !(USE_HDG)
END IF
#endif /*USE_LOADBALANCE && !(USE_HDG)*/

FVInitIsDone = .FALSE.
END SUBROUTINE FinalizeFV

END MODULE MOD_FV
