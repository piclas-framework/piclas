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

MODULE MOD_Gradients
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ---------------------------------------------------------------------------------------------------------------------
PUBLIC::InitGradients,GetGradients,FinalizeGradients,DefineParametersGradients
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersGradients()
! MODULES
USE MOD_ReadInTools       ,ONLY: prms
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("Gradients")

CALL prms%CreateIntOption('Grad-LimiterType',"Type of slope limiter of second order reconstruction", '9')
CALL prms%CreateRealOption('Grad-VktK',"K parameter for Venkatakrishnan limiter", '1.')

END SUBROUTINE DefineParametersGradients

SUBROUTINE InitGradients(ini_dim)
!===================================================================================================================================
! Initialization of the gradient variables
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools           ,ONLY: GETINT, GETREAL, GETREALARRAY
USE MOD_Mesh_Vars             ,ONLY: nElems, nSides
USE MOD_Mesh_Vars_FV          ,ONLY: Face_xGP_FV
USE MOD_Gradient_Metrics      ,ONLY: InitGradMetrics, BuildGradSideMatrix
USE MOD_Gradient_Vars
#if USE_MPI
USE MOD_MPI_Vars
USE MOD_MPI                   ,ONLY: StartReceiveMPIDataFV,StartSendMPIDataFV,FinishExchangeMPIData
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)    :: ini_dim
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT GRADIENTS...'
ALLOCATE(Grad_dx_slave(1:3,1:nSides))
ALLOCATE(Grad_dx_master(1:3,1:nSides))
ALLOCATE(Grad_SysSol_slave(1:3,1:nSides))
ALLOCATE(Grad_SysSol_master(1:3,1:nSides))
ALLOCATE(Grad_SysSol_BC(1:3,1:nSides))
Grad_dx_slave = 0.; Grad_dx_master = 0.; Grad_SysSol_slave = 0.; Grad_SysSol_master = 0.; Grad_SysSol_BC = 0.

Grad_DIM = ini_dim
ALLOCATE(Var_slave(Grad_DIM,nSides))
ALLOCATE(Var_master(Grad_DIM,nSides))
ALLOCATE(Diff_side(Grad_DIM,nSides))
ALLOCATE(Gradient_elem(3,Grad_DIM,nElems))
Var_slave = 0.; Var_master = 0.; Diff_side = 0.; Gradient_elem = 0.

GradLimiterType=GETINT('Grad-LimiterType')
GradLimVktK=GETREAL('Grad-VktK')
SELECT CASE(GradLimiterType)
CASE(0)
  LBWRITE(UNIT_stdOut,*)'Limiter = 0 -> first order FV'
CASE(1) !minmax
  LBWRITE(UNIT_stdOut,*)'Using Barth-Jespersen Limiter'
CASE(4) !venkatakrishnan
  LBWRITE(UNIT_stdOut,*)'Using Venkatakrishnan limiter with K =', GradLimVktK
CASE(9) ! no limiter (central)
  LBWRITE(UNIT_stdOut,*)'Not using any limiter'
CASE DEFAULT
  CALL abort(__STAMP__,'Limiter type not implemented.')
END SELECT

!calculate face to center distances for reconstruction
#if USE_MPI
!send face coordinates because MPI slave sides don't have them
CALL StartReceiveMPIDataFV(3,Face_xGP_FV(:,0,0,:),1,nSides,RecRequest_Geo,SendID=1) ! Receive YOUR
CALL StartSendMPIDataFV(3,Face_xGP_FV(:,0,0,:),1,nSides,SendRequest_Geo,SendID=1) ! Send MINE
CALL FinishExchangeMPIData(SendRequest_Geo,RecRequest_Geo,SendID=1)

! distances for MPI sides - send direction
CALL StartReceiveMPIDataFV(3,Grad_dx_slave,1,nSides,RecRequest_U,SendID=2) ! Receive MINE
CALL StartReceiveMPIDataFV(3,Grad_dx_master,1,nSides,RecRequest_U2,SendID=1)! Receive YOUR
CALL InitGradMetrics(doMPISides=.TRUE.)
! CALL Dx_Mortar(Grad_dx_master,Grad_dx_slave,doMPISides=.TRUE.)
! CALL Dx_Mortar2(Grad_dx_master,Grad_dx_slave,doMPISides=.TRUE.)
CALL StartSendMPIDataFV(3,Grad_dx_slave,1,nSides,SendRequest_U,SendID=2) ! Send YOUR
CALL StartSendMPIDataFV(3,Grad_dx_master,1,nSides,SendRequest_U2,SendID=1) ! Send MINE
#endif /*USE_MPI*/
! distances for BCSides, InnerSides and MPI sides - receive direction
CALL InitGradMetrics(doMPISides=.FALSE.)
! CALL Dx_Mortar(Grad_dx_master,Grad_dx_slave,doMPISides=.FALSE.)
! CALL Dx_Mortar2(Grad_dx_master,Grad_dx_slave,doMPISides=.FALSE.)
#if USE_MPI
CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=2)
CALL FinishExchangeMPIData(SendRequest_U2,RecRequest_U2,SendID=1)
#endif /*USE_MPI*/
! build neighbour-side matrices and invert them
CALL BuildGradSideMatrix(Grad_dx_master,Grad_dx_slave)

LBWRITE(UNIT_stdOut,'(A)')' INIT GRADIENTS DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitGradients

SUBROUTINE GetGradients(VarForGradients)
!===================================================================================================================================
!> Main routine for the gradient calculation
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Gradient_Vars
USE MOD_Mesh_Vars           ,ONLY: nElems,nSides,ElemToSide
#if USE_MPI
USE MOD_MPI                 ,ONLY: StartReceiveMPIDataFV,StartSendMPIDataFV,FinishExchangeMPIData
USE MOD_MPI_Vars
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers  ,ONLY: LBStartTime,LBSplitTime
#endif /*USE_LOADBALANCE*/
#endif /*MPI*/
#ifdef drift_diffusion
USE MOD_Equation_Vars_FV    ,ONLY: EFluid_GradSide
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)         :: VarForGradients(Grad_DIM,nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: ElemID, SideID, locSideID, flip
REAL                    :: gradWeight, maxDiff(Grad_DIM), minDiff(Grad_DIM)
#if USE_LOADBALANCE
REAL                            :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================

! prolong the solution to the face for grad computation
#if USE_MPI
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
! Prolong to face for MPI sides - send direction
CALL StartReceiveMPIDataFV(Grad_DIM,Var_slave,1,nSides,RecRequest_U,SendID=2) ! Receive MINE
#if USE_LOADBALANCE
CALL LBSplitTime(LB_FVCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
CALL ProlongToFace_ElemCopy(Grad_DIM,VarForGradients,Var_master,Var_slave,doMPISides=.TRUE.)
! CALL U_Mortar(Var_master,Var_slave,doMPISides=.TRUE.)
#if USE_LOADBALANCE
CALL LBSplitTime(LB_FV,tLBStart)
#endif /*USE_LOADBALANCE*/
CALL StartSendMPIDataFV(Grad_DIM,Var_slave,1,nSides,SendRequest_U,SendID=2) ! Send YOUR
#if USE_LOADBALANCE
CALL LBSplitTime(LB_FVCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/

! Prolong to face for BCSides, InnerSides and MPI sides - receive direction
CALL ProlongToFace_ElemCopy(Grad_DIM,VarForGradients,Var_master,Var_slave,doMPISides=.FALSE.)
! CALL U_Mortar(U_master,U_slave,doMPISides=.FALSE.)

#if USE_MPI
#if USE_LOADBALANCE
CALL LBSplitTime(LB_FV,tLBStart)
#endif /*USE_LOADBALANCE*/
! Complete send / receive of prolongtoface results
CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=2)
CALL StartReceiveMPIDataFV(Grad_DIM,Diff_side,1,nSides,RecRequest_gradUx,SendID=1) ! Receive YOUR
#if USE_LOADBALANCE
CALL LBSplitTime(LB_FVCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/

! fill the global neighbour difference list
CALL CalcDiff(doMPISides=.TRUE.)
#if USE_LOADBALANCE
CALL LBSplitTime(LB_FV,tLBStart)
#endif /*USE_LOADBALANCE*/

CALL StartSendMPIDataFV(Grad_DIM,Diff_side,1,nSides,SendRequest_gradUx,SendID=1) ! Send MINE
! Complete send / receive of gradients (before mpiFALSE gradients because bc grads need further grads)
CALL FinishExchangeMPIData(SendRequest_gradUx,RecRequest_gradUx,SendID=1)
#if USE_LOADBALANCE
CALL LBSplitTime(LB_FVCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/

! fill all the neighbour differences on this proc
CALL CalcDiff(doMPISides=.FALSE.)

#ifdef drift_diffusion
EFluid_GradSide(:)=Diff_side(1,:)/(SQRT((Grad_dx_master(1,:)-Grad_dx_slave(1,:))**2 &
                                       +(Grad_dx_master(2,:)-Grad_dx_slave(2,:))**2 &
                                       +(Grad_dx_master(3,:)-Grad_dx_slave(3,:))**2))
#endif

! least square method to get elem gradient from neighbour differences
DO ElemID = 1, nElems
  Gradient_elem(:,:,ElemID) = 0.
  maxDiff = 0.
  minDiff = 0.
#if PP_dim == 3
  DO locSideID=1,6
#else
  DO locSideID=2,5
#endif
    SideID=ElemToSide(E2S_SIDE_ID,locSideID,ElemID)
    flip=ElemToSide(E2S_FLIP,locSideID,ElemID)
    gradWeight = 1/VECNORM(Grad_dx_master(:,SideID)-Grad_dx_slave(:,SideID))**2
    IF (flip.EQ.0) THEN !master
      maxDiff = MAX(maxDiff,-Diff_side(:,SideID))
      minDiff = MIN(minDiff,-Diff_side(:,SideID))
      Gradient_elem(1,:,ElemID) = Gradient_elem(1,:,ElemID) &
                                          - gradWeight*Grad_SysSol_master(1,SideID)*Diff_side(:,SideID)
      Gradient_elem(2,:,ElemID) = Gradient_elem(2,:,ElemID) &
                                          - gradWeight*Grad_SysSol_master(2,SideID)*Diff_side(:,SideID)
#if PP_dim == 3
      Gradient_elem(3,:,ElemID) = Gradient_elem(3,:,ElemID) &
                                          - gradWeight*Grad_SysSol_master(3,SideID)*Diff_side(:,SideID)
#endif
    ELSE !slave
      maxDiff = MAX(maxDiff,Diff_side(:,SideID))
      minDiff = MIN(minDiff,Diff_side(:,SideID))
      Gradient_elem(1,:,ElemID) = Gradient_elem(1,:,ElemID) &
                                          + gradWeight*Grad_SysSol_slave(1,SideID)*Diff_side(:,SideID)
      Gradient_elem(2,:,ElemID) = Gradient_elem(2,:,ElemID) &
                                          + gradWeight*Grad_SysSol_slave(2,SideID)*Diff_side(:,SideID)
#if PP_dim == 3
      Gradient_elem(3,:,ElemID) = Gradient_elem(3,:,ElemID) &
                                          + gradWeight*Grad_SysSol_slave(3,SideID)*Diff_side(:,SideID)
#endif
    END IF
  END DO

  CALL GradLimiter(ElemID,minDiff,maxDiff,Gradient_elem(:,:,ElemID))

END DO

#if USE_LOADBALANCE
CALL LBSplitTime(LB_FV,tLBStart)
#endif /*USE_LOADBALANCE*/

END SUBROUTINE GetGradients

SUBROUTINE CalcDiff(doMPISides)
!===================================================================================================================================
! Compute difference between neighbour elements
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Gradient_Vars            ,ONLY: Grad_dx_master, Grad_dx_slave, Grad_SysSol_BC, Grad_DIM, Var_slave, Var_master, Diff_side
USE MOD_GetBoundaryGrad          ,ONLY: GetBoundaryGrad
USE MOD_Mesh_Vars_FV             ,ONLY: NormVec_FV,Face_xGP_FV
USE MOD_Mesh_Vars                ,ONLY: firstBCSide,lastBCSide,firstInnerSide, lastInnerSide
USE MOD_Mesh_Vars                ,ONLY: firstMPISide_MINE,lastMPISide_MINE, SideToElem, ElemToSide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)      :: doMPISides
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: SideID, SideID2, lastSideID, firstSideID_wo_BC, ElemID, locSideID, locSideID2
REAL                            :: diffUinside(Grad_DIM), gradWeight
!===================================================================================================================================
! Set the side range according to MPI or no MPI
IF(doMPISides)THEN
! fill only flux for MINE MPISides (where the local proc is master)
  firstSideID_wo_BC = firstMPISide_MINE
  lastSideID =  lastMPISide_MINE
ELSE
! fill only InnerSides that do not need communication
  firstSideID_wo_BC = firstInnerSide ! for fluxes
  lastSideID = lastInnerSide
END IF

DO SideID=firstSideID_wo_BC,lastSideID
  Diff_side(:,SideID) = Var_master(:,SideID) - Var_slave(:,SideID)
END DO

! 2. Compute the gradients at the boundary conditions: 1..nBCSides
IF (.NOT.doMPISides) THEN
  DO SideID=firstBCSide,lastBCSide
#ifdef discrete_velocity
    !second-order boundaries need bc defined from inner slopes
    diffUinside = 0.
    ElemID     = SideToElem(S2E_ELEM_ID,SideID)  !element is always master (slave=boundary ghost cell)
    locSideID  = SideToElem(S2E_LOC_SIDE_ID,SideID)
    DO locSideID2=1,6
      SideID2=ElemToSide(E2S_SIDE_ID,locSideID2,ElemID)
      IF (SideID2.LE.lastBCSide) CYCLE
      gradWeight = 1/VECNORM(Grad_dx_master(:,SideID2)-Grad_dx_slave(:,SideID2))**2
      diffUinside(:) = diffUinside(:) - Grad_dx_master(1,SideID)*gradWeight*Grad_SysSol_BC(1,SideID2)*Diff_side(:,SideID2) &
                                      - Grad_dx_master(2,SideID)*gradWeight*Grad_SysSol_BC(2,SideID2)*Diff_side(:,SideID2) &
                                      - Grad_dx_master(3,SideID)*gradWeight*Grad_SysSol_BC(3,SideID2)*Diff_side(:,SideID2)
    END DO
    CALL GetBoundaryGrad(SideID,Diff_side(:,SideID),diffUinside,&
                          Var_master(:,SideID),&
                          NormVec_FV(:,0,0,SideID),&
                          Face_xGP_FV(:,0,0,SideID))
#else
    CALL GetBoundaryGrad(SideID,Diff_side(:,SideID),&
                          Var_master(:,SideID),&
                          NormVec_FV(:,0,0,SideID),&
                          Face_xGP_FV(:,0,0,SideID))
#endif
  END DO
END IF

END SUBROUTINE CalcDiff

SUBROUTINE GradLimiter(ElemID,minDiff,maxDiff,Gradient)
!===================================================================================================================================
! Gradient limiters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Gradient_Vars       ,ONLY: Grad_dx_master, Grad_dx_slave, GradLimiterType, GradLimVktK, Grad_DIM
USE MOD_Mesh_Vars           ,ONLY: ElemToSide
USE MOD_Particle_Mesh_Vars  ,ONLY: ElemVolume_Shared
#if USE_MPI && defined(PARTICLES)
USE MOD_Mesh_Tools          ,ONLY: GetCNElemID
USE MOD_Mesh_Vars           ,ONLY: offsetElem
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)      :: ElemID
REAL,INTENT(IN)         :: minDiff(Grad_DIM)
REAL,INTENT(IN)         :: maxDiff(Grad_DIM)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)      :: Gradient(3,Grad_DIM)

!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iVar,locSideID,SideID,flip, CNElemID
REAL                            :: limiter(Grad_DIM), a, b, w
!===================================================================================================================================
limiter = 1.

#if USE_MPI && defined(PARTICLES)
CNElemID=GetCNElemID(ElemID+offSetElem)
#else
CNElemID=ElemID
#endif
w = (GradLimVktK**3)*ElemVolume_Shared(CNElemID)

#if PP_dim == 3
DO locSideID=1,6
#else
DO locSideID=2,5
#endif
  SideID=ElemToSide(E2S_SIDE_ID,locSideID,ElemID)
  flip=ElemToSide(E2S_FLIP,locSideID,ElemID)
  DO iVar=1,Grad_DIM
    SELECT CASE(GradLimiterType)
    CASE(0)
      limiter(iVar) = 0.
    CASE(1) !minmax
      IF (flip.EQ.0) THEN !master
        b = DOT_PRODUCT(Grad_dx_master(:,SideID),Gradient(:,iVar))
      ELSE !slave
        b = DOT_PRODUCT(Grad_dx_slave(:,SideID),Gradient(:,iVar))
      END IF
      IF (b.GT.0.) THEN
        limiter(iVar) = MIN(limiter(iVar),MIN(1.,maxDiff(iVar)/b))
      ELSE iF (b.LT.0.) THEN
        limiter(iVar) = MIN(limiter(iVar),MIN(1.,minDiff(iVar)/b))
      ELSE
        limiter(iVar) = MIN(limiter(iVar),1.)
      END IF
    CASE(4) !venkatakrishnan
      IF (flip.EQ.0) THEN !master
        b = DOT_PRODUCT(Grad_dx_master(:,SideID),Gradient(:,iVar))
      ELSE !slave
        b = DOT_PRODUCT(Grad_dx_slave(:,SideID),Gradient(:,iVar))
      END IF
      IF (b.GT.0.) THEN
        a = maxDiff(iVar)
        limiter(iVar) = MIN(limiter(iVar),(a**2+2*a*b+w)/(a**2+2*b**2+a*b+w))
      ELSE iF (b.LT.0.) THEN
        a = minDiff(iVar)
        limiter(iVar) = MIN(limiter(iVar),(a**2+2*a*b+w)/(a**2+2*b**2+a*b+w))
      ELSE
        limiter(iVar) = MIN(limiter(iVar),1.)
      END IF
    CASE(9) ! no limiter (central)
      ! already equal to 1
    CASE DEFAULT
      CALL abort(__STAMP__,'Limiter type not implemented for unstructured meshes.')
    END SELECT
  END DO
END DO

Gradient(1,:) = limiter(:)*Gradient(1,:)
Gradient(2,:) = limiter(:)*Gradient(2,:)
Gradient(3,:) = limiter(:)*Gradient(3,:)

END SUBROUTINE GradLimiter

SUBROUTINE ProlongToFace_ElemCopy(VarDim,ElemVar,SideVar_master,SideVar_slave,doMPISides)
!===================================================================================================================================
! Simply copies element-based variable to sides
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,          ONLY: nSides, SideToElem, nElems
USE MOD_Mesh_Vars,          ONLY: firstBCSide,firstInnerSide, lastInnerSide
USE MOD_Mesh_Vars,          ONLY: firstMPISide_YOUR,lastMPISide_YOUR,lastMPISide_MINE,firstMortarMPISide,lastMortarMPISide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)              :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides +InnerSides +MPISides MINE
INTEGER,INTENT(IN)              :: VarDim
REAL,INTENT(IN)                 :: ElemVar(VarDim,nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)              :: SideVar_master(VarDim,1:nSides)
REAL,INTENT(INOUT)              :: SideVar_slave(VarDim,1:nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: ElemID,SideID,firstSideID,lastSideID
!===================================================================================================================================
IF(doMPISides)THEN
! only YOUR MPI Sides are filled
firstSideID = firstMPISide_YOUR
lastSideID  = lastMPISide_YOUR
ELSE
! (Mortar-)InnerSides are filled
firstSideID = firstInnerSide
lastSideID  = lastInnerSide
END IF

DO SideID=firstSideID,lastSideID
! neighbor side !ElemID=-1 if not existing
ElemID     = SideToElem(S2E_NB_ELEM_ID,SideID)
!Copy Bulk Velocity in U for Dense Fokker Planck instead of Finite Volumes Solution
IF (ElemID.LT.0) CYCLE !mpi-mortar whatever
SideVar_slave(:,SideID) = ElemVar(:,ElemID)
END DO

! Second process Minus/Master sides, U_Minus is always MINE
! master side, flip=0
IF(doMPISides)THEN
! only MPI mortars
firstSideID = firstMortarMPISide
lastSideID =  lastMortarMPISide
ELSE
! BCSides, (Mortar-)InnerSides and MINE MPISides are filled
firstSideID = firstBCSide
lastSideID =  lastMPISide_MINE
END IF

DO SideID=firstSideID,lastSideID
ElemID    = SideToElem(S2E_ELEM_ID,SideID)
IF (ElemID.LT.0) CYCLE !for small mortar sides without info on big master element
SideVar_master(:,SideID)=ElemVar(:,ElemID)
END DO !SideID

END SUBROUTINE ProlongToFace_ElemCopy

SUBROUTINE FinalizeGradients()
!===================================================================================================================================
! deallocate gradient variables
!===================================================================================================================================
! MODULES
USE MOD_Gradient_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(Grad_dx_slave)
SDEALLOCATE(Grad_dx_master)
SDEALLOCATE(Grad_SysSol_slave)
SDEALLOCATE(Grad_SysSol_master)
SDEALLOCATE(Grad_SysSol_BC)
SDEALLOCATE(Var_slave)
SDEALLOCATE(Var_master)
SDEALLOCATE(Diff_side)
SDEALLOCATE(Gradient_elem)

END SUBROUTINE FinalizeGradients

END MODULE MOD_Gradients

