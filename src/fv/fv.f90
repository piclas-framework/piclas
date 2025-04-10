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
! Contains the initialization of the FV global variables
! Computes the different FV spatial operators/residuals(Ut_FV) using U_FV
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

PUBLIC::InitFV,FinalizeFV,FV_main
!===================================================================================================================================

CONTAINS

SUBROUTINE InitFV()
!===================================================================================================================================
! Allocate global variable U_FV (solution) and Ut_FV (FV time derivative).
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_FV_Vars
USE MOD_Restart_Vars       ,ONLY: DoRestart,RestartInitIsDone
USE MOD_Interpolation_Vars ,ONLY: InterpolationInitIsDone
USE MOD_Mesh_Vars          ,ONLY: nSides
USE MOD_Mesh_Vars          ,ONLY: MeshInitIsDone
USE MOD_Gradients          ,ONLY: InitGradients
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
USE MOD_LoadBalance_Vars   ,ONLY: UseH5IOLoadBalance
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
!===================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.(.NOT.MeshInitIsDone).OR.(.NOT.RestartInitIsDone).OR.FVInitIsDone) CALL abort(__STAMP__,&
    'InitFV not ready to be called or already called.')
LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT FV...'

#if USE_LOADBALANCE
IF (.NOT.(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance))) THEN
#endif /*USE_LOADBALANCE*/
  ! the local FV solution
  ALLOCATE( U_FV(PP_nVar_FV,0:0,0:0,0:0,PP_nElems))
  U_FV=0.
#if USE_LOADBALANCE
END IF
#endif /*USE_LOADBALANCE)*/
! the time derivative computed with the FV scheme
ALLOCATE(Ut_FV(PP_nVar_FV,0:0,0:0,0:0,PP_nElems))
Ut_FV=0.

! U_FV is filled with the ini solution
IF(.NOT.DoRestart) CALL FillIni()

#ifdef drift_diffusion
ALLOCATE(U_master_FV(PP_nVar_FV+3,0:0,0:0,1:nSides))
ALLOCATE(U_slave_FV(PP_nVar_FV+3,0:0,0:0,1:nSides))
#else
ALLOCATE(U_master_FV(PP_nVar_FV,0:0,0:0,1:nSides))
ALLOCATE(U_slave_FV(PP_nVar_FV,0:0,0:0,1:nSides))
#endif

U_master_FV=0.
U_slave_FV=0.

! unique flux per side
ALLOCATE(Flux_Master_FV(1:PP_nVar_FV,0:0,0:0,1:nSides))
ALLOCATE(Flux_Slave_FV (1:PP_nVar_FV,0:0,0:0,1:nSides))
Flux_Master_FV=0.
Flux_Slave_FV=0.

CALL InitGradients(PP_nVar_FV)

FVInitIsDone=.TRUE.
LBWRITE(UNIT_stdOut,'(A)')' INIT FV DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitFV


SUBROUTINE FV_main(t,tStage,doSource)
!===================================================================================================================================
! Computes the FV time derivative
! U_FV and Ut_FV are allocated
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Vector
USE MOD_FV_Vars           ,ONLY: U_FV,Ut_FV,Flux_Master_FV,Flux_Slave_FV
USE MOD_FV_Vars           ,ONLY: U_master_FV,U_slave_FV
USE MOD_SurfInt           ,ONLY: SurfInt
USE MOD_Prolong_FV        ,ONLY: ProlongToFace_FV
USE MOD_Gradients         ,ONLY: GetGradients
USE MOD_FillFlux          ,ONLY: FillFlux
USE MOD_Equation_FV       ,ONLY: CalcSource_FV
USE MOD_Interpolation     ,ONLY: ApplyJacobian
USE MOD_FillMortar        ,ONLY: U_Mortar,Flux_Mortar
USE MOD_Particle_Mesh_Vars,ONLY: ElemVolume_Shared
USE MOD_Interpolation_Vars,ONLY: wGP
#if USE_MPI
USE MOD_Mesh_Vars         ,ONLY: nSides
USE MOD_MPI_Vars
USE MOD_MPI               ,ONLY: StartReceiveMPIDataFV,StartSendMPIDataFV,FinishExchangeMPIData
#if defined(PARTICLES) && defined(LSERK)
USE MOD_Particle_Vars     ,ONLY: DelayTime
USE MOD_TimeDisc_Vars     ,ONLY: time
#endif /*defined(PARTICLES) && defined(LSERK)*/
#ifdef PARTICLES
USE MOD_Particle_MPI      ,ONLY: MPIParticleSend,MPIParticleRecv
USE MOD_Mesh_Vars         ,ONLY: offsetElem
USE MOD_Mesh_Tools        ,ONLY: GetCNElemID
#endif /*PARTICLES*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers,ONLY: LBStartTime,LBPauseTime,LBSplitTime
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/
#ifdef drift_diffusion
USE MOD_Equation_Vars     ,ONLY: E
USE MOD_MPI               ,ONLY: StartReceiveMPIData,StartSendMPIData
#endif
REAL,INTENT(IN)                 :: t,tStage
LOGICAL,INTENT(IN)              :: doSource
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: CNElemID, iElem,i,j,k,ElemID
#ifdef drift_diffusion
REAL                            :: U_DD(1:PP_nVar_FV+3,0:0,0:0,0:0,PP_nElems) ! U_FV(1:PP_nVar_FV) + E(1:3)
#endif
#if USE_LOADBALANCE
REAL                            :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================

! Compute the FV solution gradients (LB times are measured inside GetGradients)
CALL GetGradients(U_FV(:,0,0,0,:),output=.FALSE.) ! this might trigger a copy of U_FV -> the useless dimensions should be removed someday

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

#ifdef drift_diffusion
U_DD(:,:,:,:,:) = 0.
DO ElemID = 1, PP_nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    U_DD(PP_nVar_FV+1:PP_nVar_FV+3,0,0,0,ElemID) = U_DD(PP_nVar_FV+1:PP_nVar_FV+3,0,0,0,ElemID) &
                                                 + wGP(i)*wGP(j)*wGP(k)*E(1:3,i,j,k,ElemID)/((PP_N+1.)**3) !need jacobi here for noncartesian
  END DO; END DO; END DO
  U_DD(1:PP_nVar_FV,0,0,0,ElemID) = U_FV(1:PP_nVar_FV,0,0,0,ElemID)
END DO

ASSOCIATE(  PP_nVar_tmp => PP_nVar_FV+3 , &
            U_tmp => U_DD)
#else
ASSOCIATE(  PP_nVar_tmp => PP_nVar_FV, &
            U_tmp => U_FV)
#endif

! prolong the solution to the faces by applying reconstruction gradients
#if USE_MPI
! Prolong to face for MPI sides - send direction
#if USE_LOADBALANCE
CALL LBSplitTime(LB_FV,tLBStart)
#endif /*USE_LOADBALANCE*/
CALL StartReceiveMPIDataFV(PP_nVar_tmp,U_slave_FV(:,0,0,:),1,nSides,RecRequest_U,SendID=2) ! Receive MINE
#if USE_LOADBALANCE
CALL LBSplitTime(LB_FVCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/

CALL ProlongToFace_FV(U_tmp,U_master_FV,U_slave_FV,doMPISides=.TRUE.)
CALL U_Mortar(U_master_FV,U_slave_FV,doMPISides=.TRUE.) !not working

#if USE_LOADBALANCE
CALL LBSplitTime(LB_FV,tLBStart)
#endif /*USE_LOADBALANCE*/
CALL StartSendMPIDataFV(PP_nVar_tmp,U_slave_FV(:,0,0,:),1,nSides,SendRequest_U,SendID=2) ! Send YOUR
#if USE_LOADBALANCE
CALL LBSplitTime(LB_FVCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/

! Prolong to face for BCSides, InnerSides and MPI sides - receive direction
CALL ProlongToFace_FV(U_tmp,U_master_FV,U_slave_FV,doMPISides=.FALSE.)
CALL U_Mortar(U_master_FV,U_slave_FV,doMPISides=.FALSE.) !not working

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
Ut_FV=0.

#if USE_MPI
#if USE_LOADBALANCE
CALL LBSplitTime(LB_FV,tLBStart)
#endif /*USE_LOADBALANCE*/
! Complete send / receive of prolongtoface results
CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=2)

END ASSOCIATE

! Initialization of the time derivative
!Flux=0. !don't nullify the fluxes if not really needed (very expensive)
CALL StartReceiveMPIDataFV(PP_nVar_FV,Flux_Slave_FV(:,0,0,:),1,nSides,RecRequest_Flux,SendID=1) ! Receive YOUR
#if USE_LOADBALANCE
CALL LBSplitTime(LB_FVCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
! fill the global surface flux list
CALL FillFlux(t,Flux_Master_FV,Flux_Slave_FV,U_master_FV,U_slave_FV,doMPISides=.TRUE.)
#if USE_LOADBALANCE
CALL LBSplitTime(LB_FV,tLBStart)
#endif /*USE_LOADBALANCE*/
CALL StartSendMPIDataFV(PP_nVar_FV,Flux_Slave_FV(:,0,0,:),1,nSides,SendRequest_Flux,SendID=1) ! Send MINE

#if USE_LOADBALANCE
CALL LBSplitTime(LB_FVCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/

! fill the all surface fluxes on this proc
CALL FillFlux(t,Flux_Master_FV,Flux_Slave_FV,U_master_FV,U_slave_FV,doMPISides=.FALSE.)
CALL Flux_Mortar(Flux_Master_FV,Flux_Slave_FV,doMPISides=.FALSE.)
! compute surface integral contribution and add to ut
CALL SurfInt(Flux_Master_FV,Flux_Slave_FV,Ut_FV,doMPISides=.FALSE.)

#if USE_MPI
#if USE_LOADBALANCE
CALL LBSplitTime(LB_FV,tLBStart)
#endif /*USE_LOADBALANCE*/
! Complete send / receive
CALL FinishExchangeMPIData(SendRequest_Flux,RecRequest_Flux,SendID=1)
#if USE_LOADBALANCE
CALL LBSplitTime(LB_FVCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/

!FINALIZE Fluxes for MPI Sides
CALL Flux_Mortar(Flux_Master_FV,Flux_Slave_FV,doMPISides=.TRUE.)
CALL SurfInt(Flux_Master_FV,Flux_Slave_FV,Ut_FV,doMPISides=.TRUE.)
#if USE_LOADBALANCE
CALL LBSplitTime(LB_FV,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif

! Swap to right sign and divide by element volume (1/sJ)
! CALL ApplyJacobian(Ut_FV,toPhysical=.TRUE.,toSwap=.TRUE.)
DO iElem=1,PP_nElems
#if USE_MPI && defined(PARTICLES)
  CNElemID=GetCNElemID(iElem+offSetElem)
#else
  CNElemID=iElem
#endif
  Ut_FV(:,:,:,:,iElem)=-Ut_FV(:,:,:,:,iElem)/ElemVolume_Shared(CNElemID)
END DO

! Add Source Terms
IF(doSource) CALL CalcSource_FV(tStage,1.0,Ut_FV)

#if USE_LOADBALANCE
CALL LBSplitTime(LB_FV,tLBStart)
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
USE MOD_FV_Vars,ONLY:U_FV
USE MOD_Mesh_Vars_FV,ONLY:Elem_xGP_FV
USE MOD_Equation_Vars_FV,ONLY:IniExactFunc_FV
USE MOD_Equation_FV,ONLY:ExactFunc_FV
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
DO iElem=1,PP_nElems
    CALL ExactFunc_FV(IniExactFunc_FV,0.,0,Elem_xGP_FV(1:3,0,0,0,iElem),U_FV(1:PP_nVar_FV,0,0,0,iElem))
END DO ! iElem=1,PP_nElems
END SUBROUTINE FillIni



SUBROUTINE FinalizeFV()
!===================================================================================================================================
! Deallocate global variable U_FV (solution) and Ut_FV (dg time derivative).
!===================================================================================================================================
! MODULES
USE MOD_FV_Vars
USE MOD_Gradients          ,ONLY: FinalizeGradients
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(Ut_FV)
SDEALLOCATE(U_master_FV)
SDEALLOCATE(U_slave_FV)
SDEALLOCATE(Flux_Master_FV)
SDEALLOCATE(Flux_Slave_FV)
CALL FinalizeGradients()
! SDEALLOCATE(U_Master_loc)
! SDEALLOCATE(U_Slave_loc)
! SDEALLOCATE(Flux_loc)

! Do not deallocate the solution vector during load balance here as it needs to be communicated between the processors
#if USE_LOADBALANCE
IF(.NOT.(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)))THEN
#endif /*USE_LOADBALANCE*/
  SDEALLOCATE(U_FV)
#if USE_LOADBALANCE
END IF
#endif /*USE_LOADBALANCE*/

FVInitIsDone = .FALSE.
END SUBROUTINE FinalizeFV

END MODULE MOD_FV