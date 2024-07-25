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
#if USE_PETSC
#include "petsc/finclude/petsc.h"
#endif

!===================================================================================================================================
!> Module for the HDG method
!===================================================================================================================================
MODULE MOD_HDG_Linear
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
#if USE_HDG
PUBLIC :: HDGLinear
#endif /*USE_HDG*/
!===================================================================================================================================

CONTAINS

#if USE_HDG
!===================================================================================================================================
!> Linear HDG solver
!===================================================================================================================================
SUBROUTINE HDGLinear(time)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_HDG_Vars
USE MOD_DG_Vars            ,ONLY: DG_Elems_master,N_DG_Mapping,U_N
USE MOD_Equation           ,ONLY: CalcSourceHDG,ExactFunc
USE MOD_Equation_Vars      ,ONLY: IniExactFunc
USE MOD_Mesh_Vars          ,ONLY: BoundaryType,nSides,BC,N_SurfMesh
USE MOD_Mesh_Vars          ,ONLY: ElemToSide, offSetElem
USE MOD_Interpolation_Vars ,ONLY: NMax,PREF_VDM
USE MOD_Elem_Mat           ,ONLY: PostProcessGradientHDG
#if USE_PETSC
USE MOD_Elem_Mat           ,ONLY: ChangeBasisSmat
#endif /*USE_PETSC*/
USE MOD_FillMortar_HDG     ,ONLY: SmallToBigMortar_HDG
#if (PP_nVar==1)
!USE MOD_Equation_Vars      ,ONLY: E
#elif (PP_nVar==3)
USE MOD_Equation_Vars      ,ONLY: B
#else
USE MOD_Equation_Vars      ,ONLY: B, E
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers ,ONLY: LBStartTime,LBPauseTime,LBSplitTime
#endif /*USE_LOADBALANCE*/
#if USE_PETSC
USE PETSc
USE MOD_Mesh_Vars          ,ONLY: SideToElem
#if USE_MPI
USE MOD_MPI                ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
USE MOD_MPI_Vars
#endif
USE MOD_FillMortar_HDG     ,ONLY: BigToSmallMortar_HDG
#endif
#if USE_MPI
USE MOD_MPI                ,ONLY: Mask_MPIsides
#endif
USE MOD_Globals_Vars       ,ONLY: ElementaryCharge,eps0
USE MOD_ChangeBasis        ,ONLY: ChangeBasis2D
USE MOD_HDG_Tools          ,ONLY: CG_solver,DisplayConvergence
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: time !time
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: lambdatmp(1:nGP_face(NMax))
INTEGER :: i,j,k,r,p,q,iElem, iVar,NSideMin,Nloc, jNloc, iNloc
INTEGER :: BCsideID,BCType,BCState,SideID,iLocSide
REAL    :: RHS_facetmp(nGP_face(NMax))
REAL    :: rtmp(nGP_vol(NMax))
INTEGER :: DOFindices(nGP_face(NMax))

!LOGICAL :: converged
#if (PP_nVar!=1)
REAL    :: BTemp(3,3,nGP_vol,PP_nElems)
#endif
#if USE_LOADBALANCE
REAL    :: tLBStart
#endif /*USE_LOADBALANCE*/
#if USE_PETSC
PetscErrorCode       :: ierr
PetscScalar, POINTER :: lambda_pointer(:)
KSPConvergedReason   :: reason
PetscInt             :: iterations
PetscReal            :: petscnorm
INTEGER              :: ElemID,iBCSide,PETScLocalID
INTEGER              :: DOF_start, DOF_stop
REAL                 :: timeStartPiclas,timeEndPiclas
REAL                 :: RHS_conductor(nGP_face(NMax))
INTEGER              :: jLocSide
REAL                 :: Smatloc(nGP_face(NMax),nGP_face(NMax))
#endif
#if USE_PETSC
INTEGER              :: iUniqueFPCBC
#endif /*USE_PETSC*/
!===================================================================================================================================
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
DO iVar = 1, PP_nVar
  !Dirichlet boundary conditions
#if (PP_nVar!=1)
  IF (iVar.EQ.4) THEN
#endif
  DO BCsideID=1,nDirichletBCSides
    SideID=DirichletBC(BCsideID)
    Nloc = N_SurfMesh(SideID)%NSideMin
    BCType =BoundaryType(BC(SideID),BC_TYPE)
    BCState=BoundaryType(BC(SideID),BC_STATE)
    SELECT CASE(BCType)
    CASE(2) ! exact BC = Dirichlet BC !! ExactFunc via BCState (time is optional)
      ! Determine the exact BC state
      DO q=0,Nloc; DO p=0,Nloc
        r=q*(Nloc+1) + p+1
        CALL ExactFunc(BCState,N_SurfMesh(SideID)%Face_xGP(:,p,q),HDG_Surf_N(SideID)%lambda(iVar,r:r),t=time)
      END DO; END DO !p,q
    CASE(4) ! exact BC = Dirichlet BC !! Zero potential
      DO q=0,Nloc; DO p=0,Nloc
        r=q*(Nloc+1) + p+1
       HDG_Surf_N(SideID)%lambda(iVar,r:r)=0.
      END DO; END DO !p,q
    CASE(5,51,52,60) ! exact BC = Dirichlet BC !! ExactFunc via RefState (time is optional) for reference state (with zero crossing)
      DO q=0,Nloc; DO p=0,Nloc
        r=q*(Nloc+1) + p+1
        CALL ExactFunc(-1,N_SurfMesh(SideID)%Face_xGP(:,p,q),HDG_Surf_N(SideID)%lambda(iVar,r:r),t=time,iRefState=BCState)
      END DO; END DO !p,q
    CASE(6) ! exact BC = Dirichlet BC !! ExactFunc via RefState (time is optional) for reference state (without zero crossing)
      DO q=0,Nloc; DO p=0,Nloc
        r=q*(Nloc+1) + p+1
        CALL ExactFunc(-2,N_SurfMesh(SideID)%Face_xGP(:,p,q),HDG_Surf_N(SideID)%lambda(iVar,r:r),t=time,iRefState=BCState)
      END DO; END DO !p,q
    CASE(7) ! exact BC = Dirichlet BC !! ExactFunc via LinState (time is optional)for linear potential function
      DO q=0,Nloc; DO p=0,Nloc
        r=q*(Nloc+1) + p+1
        CALL ExactFunc(-3,N_SurfMesh(SideID)%Face_xGP(:,p,q),HDG_Surf_N(SideID)%lambda(iVar,r:r),t=time,iLinState=BCState)
      END DO; END DO !p,q
    CASE(8) ! exact BC = Dirichlet BC !! ExactFunc via electric potential and decharing of a surface
      DO q=0,Nloc; DO p=0,Nloc
        r=q*(Nloc+1) + p+1
        CALL ExactFunc(-4,N_SurfMesh(SideID)%Face_xGP(:,p,q),HDG_Surf_N(SideID)%lambda(iVar,r:r),t=time,BCState=BCState)
      END DO; END DO !p,q
    CASE(50) ! exact BC = Dirichlet BC !! ExactFunc via DC bias voltage
      DO q=0,Nloc; DO p=0,Nloc
        r=q*(Nloc+1) + p+1
        CALL ExactFunc(-5,N_SurfMesh(SideID)%Face_xGP(:,p,q),HDG_Surf_N(SideID)%lambda(iVar,r:r),t=time,BCState=BCState)
      END DO; END DO !p,q
    END SELECT ! BCType
  END DO !BCsideID=1,nDirichletBCSides
#if (PP_nVar!=1)
  END IF
#endif

  !neumann BC
  DO BCsideID=1,nNeumannBCSides
    SideID  = NeumannBC(BCsideID)
    Nloc    = DG_Elems_master(SideID)
    BCType  = BoundaryType(BC(SideID),BC_TYPE)
    BCState = BoundaryType(BC(SideID),BC_STATE)
    SELECT CASE(BCType)
    CASE(10) !neumann q=0 !! Zero gradient
      DO q=0,Nloc; DO p=0,Nloc
        r=q*(Nloc+1) + p+1
        HDG_Surf_N(SideID)%qn_face(iVar,r)= 0.
      END DO; END DO !p,q
    END SELECT ! BCType
  END DO !BCsideID=1,nNeumannBCSides

!for magnetostatic only neumann
#if (PP_nVar!=1)
  IF (iVar.LT.4) THEN
    DO BCsideID=1,nDirichletBCSides
!      SideID=DirichletBC(BCsideID)
      DO q=0,Nloc; DO p=0,Nloc
        r=q*(Nloc+1) + p+1
        qn_face_MagStat(iVar,r,BCSideID)= 0.
      END DO; END DO !p,q
    END DO !BCsideID=1,nDirichletBCSides
  END IF
#endif

  ! Floating boundary BCs
#if USE_PETSC
  IF(UseFPC) THEN
#if USE_MPI
    ! Communicate the accumulated charged on each BC to all processors on the communicator
    DO iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
      ASSOCIATE( COMM => FPC%COMM(iUniqueFPCBC)%UNICATOR)
        IF(FPC%COMM(iUniqueFPCBC)%UNICATOR.NE.MPI_COMM_NULL)THEN
          CALL MPI_ALLREDUCE(MPI_IN_PLACE, FPC%ChargeProc(iUniqueFPCBC), 1, MPI_DOUBLE_PRECISION, MPI_SUM, COMM, IERROR)
          FPC%Charge(iUniqueFPCBC) = FPC%Charge(iUniqueFPCBC) + FPC%ChargeProc(iUniqueFPCBC)
        END IF ! FPC%COMM(iUniqueFPCBC)%UNICATOR.NE.MPI_COMM_NULL
      END ASSOCIATE
    END DO ! iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
#else
    FPC%Charge = FPC%Charge + FPC%ChargeProc
#endif /*USE_MPI*/
    FPC%ChargeProc = 0.
    ! Apply charge to RHS, which this is done below: RHS_conductor(1)=FPC%Charge(iUniqueFPCBC)/eps0
  END IF ! UseFPC
#endif /*USE_PETSC*/

  ! Set potential to zero (only one process does this)
  IF(SetZeroPotentialDOF) HDG_Surf_N(1)%lambda(iVar,1) = 0.
END DO

!volume source (volume RHS of u system)
DO iElem=1,PP_nElems
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
    r=k*(Nloc+1)**2+j*(Nloc+1) + i+1
    CALL CalcSourceHDG(i,j,k,iElem,HDG_Vol_N(iElem)%RHS_vol(1:PP_nVar,r))
  END DO; END DO; END DO !i,j,k
  DO iVar = 1, PP_nVar
    HDG_Vol_N(iElem)%RHS_Vol(iVar,:)=-HDG_Vol_N(iElem)%JwGP_vol(:)*HDG_Vol_N(iElem)%RHS_vol(iVar,:)
  END DO
END DO !iElem

!replace lambda with exact function (debugging)
IF(ExactLambda)THEN
  DO SideID=1,nSides
    Nloc = N_SurfMesh(SideID)%NSideMin
    DO q=0,Nloc; DO p=0,Nloc
      r=q*(Nloc+1) + p+1
      CALL ExactFunc(IniExactFunc,N_SurfMesh(SideID)%Face_xGP(:,p,q),HDG_Surf_N(SideID)%lambda( 1:PP_nVar,r))
    END DO; END DO !p,q
  END DO
END IF

!prepare RHS_face ( RHS for lambda system.)
DO iVar = 1, PP_nVar
  DO SideID = 1, nSides
    HDG_Surf_N(SideID)%RHS_face(iVar,:)=0.
  END DO ! SideID = 1, nSides
  DO iElem=1,PP_nElems
    Nloc = N_DG_Mapping(2,iElem+offSetElem)
    !rtmp=MATMUL(InvDhat(:,:,iElem),-RHS_loc(:,iElem))
    CALL DSYMV('U',nGP_vol(Nloc),1., HDG_Vol_N(iElem)%InvDhat(:,:),nGP_vol(Nloc), &
                               -HDG_Vol_N(iElem)%RHS_vol(iVar,:),1,0., &
                               rtmp(1:nGP_face(Nloc)),1)

    DO iLocSide=1,6
      SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
      CALL DGEMV('N',nGP_face(Nloc),nGP_vol(Nloc),1., &
                          HDG_Vol_N(iElem)%Ehat(:,:,iLocSide), nGP_face(Nloc), &
                          rtmp(1:nGP_face(Nloc)),1,0.,& ! 1: add to RHS_face, 0: set value
                          RHS_facetmp(1:nGP_face(Nloc)),1)

      NSideMin = N_SurfMesh(SideID)%NSideMin
      IF(Nloc.EQ.NSideMin)THEN
        HDG_Surf_N(SideID)%RHS_face(iVar,:) = HDG_Surf_N(SideID)%RHS_face(iVar,:) + RHS_facetmp(1:nGP_face(Nloc))
      ELSE
        ! From high to low
        CALL ChangeBasis2D(1, Nloc, NSideMin, TRANSPOSE(PREF_VDM(NSideMin,Nloc)%Vdm) , RHS_facetmp(1:nGP_face(Nloc)), RHS_facetmp(1:nGP_face(NSideMin)))
        HDG_Surf_N(SideID)%RHS_face(iVar,:) = HDG_Surf_N(SideID)%RHS_face(iVar,:) + RHS_facetmp(1:nGP_face(NSideMin))
      END IF ! Nloc.NE.NSideMin
    END DO
  END DO !iElem
END DO !ivar

!add Neumann
DO BCsideID=1,nNeumannBCSides
  SideID=NeumannBC(BCsideID)
  HDG_Surf_N(SideID)%RHS_face(:,:) = HDG_Surf_N(SideID)%RHS_face(:,:) + HDG_Surf_N(SideID)%qn_face(:,:)
END DO

#if USE_PETSC
! add Dirichlet contribution
DO iBCSide=1,nDirichletBCSides
  BCSideID=DirichletBC(iBCSide)
  ElemID    = SideToElem(S2E_ELEM_ID,BCSideID)
  jLocSide  = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
  ! At BCSides, jNLoc = Nmax
  jNloc     = N_SurfMesh(BCsideID)%NSideMin
  DO iLocSide=1,6
    SideID = ElemToSide(E2S_SIDE_ID,iLocSide,ElemID)
    iNloc     = N_SurfMesh(SideID)%NSideMin
    IF(MaskedSide(SideID).GT.0) CYCLE

    ! TODO PETSC P-Adaption - Improvement: Store V^T * S * V in Smat
    CALL ChangeBasisSmat(Smatloc, HDG_Vol_N(ElemID)%Smat(:,:,iLocSide,jLocSide), jNloc, iNloc, jNloc)

    CALL DGEMV('N',nGP_face(iNloc),nGP_face(jNloc),-1., &
                          Smatloc(1:nGP_face(iNloc),1:nGP_face(jNloc)), nGP_face(iNloc), &
                          HDG_Surf_N(BCSideID)%lambda(1,:),1,1.,& !add to RHS_face
                          HDG_Surf_N(SideID)%RHS_face(1,:),1)
  END DO
END DO
#endif /*USE_PETSC*/

#if (PP_nVar!=1)
DO iVar = 1, PP_nVar
  IF (iVar.LT.4) THEN
    DO BCsideID=1,nDirichletBCSides
      SideID=DirichletBC(BCsideID)
      RHS_face(iVar,:,SideID)=RHS_face(iVar,:,SideID)+qn_face_MagStat(iVar,:,BCSideID)
    END DO !BCsideID=1,nDirichletBCSides
  END IF
END DO
#endif

DO iVar=1, PP_nVar
#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/

#if USE_MPI
CALL Mask_MPIsides('RHS_face',iVar)
#endif /*USE_MPI*/

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
#if !USE_PETSC
CALL SmallToBigMortar_HDG(iVar,0) ! RHS_face(1:PP_nVar,1:nGP_Face,1:nSides))
#endif /*USE_PETSC*/

#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/

! SOLVE

#if USE_PETSC
! Fill right hand side
PetscCallA(VecZeroEntries(PETScRHS,ierr))
TimeStartPiclas=PICLASTIME()
DO SideID=1,nSides
  IF(MaskedSide(SideID).GT.0) CYCLE

  ! TODO Create a function to map localToGlobalDOFs
  Nloc = N_SurfMesh(SideID)%NSideMin
  DO i=1,nGP_face(Nloc)
    DOFindices(i) = i + OffsetGlobalPETScDOF(SideID) - 1
  END DO

  PetscCallA(VecSetValues(PETScRHS,nGP_face(Nloc),DOFindices(1:nGP_face(Nloc)),HDG_Surf_N(SideID)%RHS_face(1,:),ADD_VALUES,ierr))
END DO

! TODO We had to use ADD_VALUES when filling the RHS. Maybe loop over sideID=1,nSides-nSides_YOUR?
PetscCallA(VecAssemblyBegin(PETScRHS,ierr))
PetscCallA(VecAssemblyEnd(PETScRHS,ierr))

! The MPIRoot process has charge and voltage of all FPCs, there, this process sets all conductor RHS information
IF(UseFPC) THEN
  IF(MPIRoot)THEN
    DO iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
      PetscCallA(VecSetValues(PETScRHS,1,nGlobalPETScDOFs-1-FPC%nUniqueFPCBounds+iUniqueFPCBC,FPC%Charge(iUniqueFPCBC)/eps0,INSERT_VALUES,ierr))
    END DO
  END IF
END IF

! Reset the RHS of the first DOF if ZeroPotential must be set
IF(MPIroot .AND. SetZeroPotentialDOF) THEN
  PetscCallA(VecSetValue(PETScRHS,0,0,INSERT_VALUES,ierr))
END IF

PetscCallA(VecAssemblyBegin(PETScRHS,ierr))
PetscCallA(VecAssemblyEnd(PETScRHS,ierr))

! Calculate lambda
PetscCallA(KSPSolve(PETScSolver,PETScRHS,PETScSolution,ierr))
TimeEndPiclas=PICLASTIME()
PetscCallA(KSPGetIterationNumber(PETScSolver,iterations,ierr))
PetscCallA(KSPGetConvergedReason(PETScSolver,reason,ierr))
PetscCallA(KSPGetResidualNorm(PETScSolver,petscnorm,ierr))
! reason - negative value indicates diverged, positive value converged, see KSPConvergedReason
!  -2: KSP_DIVERGED_NULL
!  -3: KSP_DIVERGED_ITS            -> Ran out of iterations before any convergence criteria was reached
!  -4: KSP_DIVERGED_DTOL           -> norm(r) >= dtol*norm(b)
!  -5: KSP_DIVERGED_BREAKDOWN      -> Breakdown in the Krylov method: the method could not continue to enlarge the Krylov space
!  -6: KSP_DIVERGED_BREAKDOWN_BICG -> Breakdown in the Krylov method: the method could not continue to enlarge the Krylov space
!  -7: KSP_DIVERGED_NONSYMMETRIC   -> It appears the operator or preconditioner is not symmetric and this Krylov method
!  -8: KSP_DIVERGED_INDEFINITE_PC  -> It appears the preconditioner is indefinite
!  -9: KSP_DIVERGED_NANORINF
! -10: KSP_DIVERGED_INDEFINITE_MAT
! -11: KSP_DIVERGED_PC_FAILED      -> It was not possible to build or use the requested preconditioner
! -11: KSP_DIVERGED_PCSETUP_FAILED_DEPRECATED
IF(reason.LT.0)THEN
  SWRITE(*,*) 'Attention: PETSc not converged! Reason: ', reason
END IF
IF(MPIroot) CALL DisplayConvergence(TimeEndPiclas-TimeStartPiclas, iterations, petscnorm)

! Fill element local lambda for post processing
! Get the local DOF subarray
PetscCallA(VecScatterBegin(PETScScatter, PETScSolution, PETScSolutionLocal, INSERT_VALUES, SCATTER_FORWARD,ierr))
PetscCallA(VecScatterEnd(  PETScScatter, PETScSolution, PETScSolutionLocal, INSERT_VALUES, SCATTER_FORWARD,ierr))
PetscCallA(VecGetArrayReadF90(PETScSolutionLocal,lambda_pointer,ierr))
DOF_stop = 0
DO SideID=1,nSides
  IF(MaskedSide(SideID).GT.0) CYCLE
  Nloc = N_SurfMesh(SideID)%NSideMin
  DOF_start = 1 + DOF_stop
  DOF_stop = DOF_start + nGP_face(Nloc) - 1
  ! MORTARS: All Mortar stuff (BigToSmall, ...) is done when solving the system!
  HDG_Surf_N(SideID)%lambda(1,:) = lambda_pointer(DOF_start:DOF_stop)
END DO

! Fill Conductor lambda
IF(UseFPC) THEN
  FPC%VoltageProc = 0. ! nullify just to be safe
  DO BCsideID=1,nConductorBCsides
    SideID       = ConductorBC(BCSideID)
    Nloc         = N_SurfMesh(SideID)%NSideMin
    BCState      = BoundaryType(BC(SideID),BC_STATE)
    iUniqueFPCBC = FPC%Group(BCState,2)
    DO i=1,nGP_face(Nloc)
      HDG_Surf_N(SideID)%lambda(1,i) = lambda_pointer(nLocalPETScDOFs - FPC%nUniqueFPCBounds + iUniqueFPCBC)
    END DO
    ! Copy the value to the FPC container for output to .csv (only mpi root does this)
    FPC%VoltageProc(iUniqueFPCBC) = HDG_Surf_N(SideID)%lambda(1,1)
  END DO
  IF(MPIRoot)THEN
    DO iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
      FPC%Voltage(iUniqueFPCBC) = lambda_pointer(nLocalPETScDOFs - FPC%nUniqueFPCBounds + iUniqueFPCBC)
    END DO ! iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
  END IF ! MPIRoot
END IF ! UseFPC

PetscCallA(VecRestoreArrayReadF90(PETScSolutionLocal,lambda_pointer,ierr))
#else
  ! HDGLinear
  CALL CG_solver(iVar)
#endif /*USE_PETSC*/
  !POST PROCESSING

#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
  !post processing:
  DO iElem=1,PP_nElems
    Nloc = N_DG_Mapping(2,iElem+offSetElem)
    ! for post-proc
    DO iLocSide=1,6
      SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
      NSideMin = N_SurfMesh(SideID)%NSideMin
      IF(Nloc.EQ.NSideMin)THEN
        lambdatmp(1:nGP_face(Nloc)) = HDG_Surf_N(SideID)%lambda(iVar,:)
      ELSE
        ! From low to high
        CALL ChangeBasis2D(1, NSideMin, Nloc, PREF_VDM(NSideMin,Nloc)%Vdm , HDG_Surf_N(SideID)%lambda(iVar,1:nGP_face(NSideMin)), lambdatmp(1:nGP_face(Nloc)))
      END IF ! Nloc.EQ.NSideMin
      CALL DGEMV('T',nGP_face(Nloc),nGP_vol(Nloc),1., &
                          HDG_Vol_N(iElem)%Ehat(:,:,iLocSide), nGP_face(Nloc), &
                          lambdatmp(1:nGP_face(Nloc)),1,1.,& ! 1: add to RHS_face, 0: set value
                          HDG_Vol_N(iElem)%RHS_vol(iVar,:),1)
    END DO
    !U(:,iElem)=MATMUL(InvDhat(:,:,iElem),-RHS_loc(:,iElem))
    CALL DSYMV('U',nGP_vol(Nloc),1., HDG_Vol_N(iElem)%InvDhat(:,:),nGP_vol(Nloc), &
                               -HDG_Vol_N(iElem)%RHS_vol(iVar,:),1,0., &
                               U_N(iElem)%U(iVar,:,:,:),1)
  END DO !iElem
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/
END DO !iVar

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
#if (PP_nVar==1)
  CALL PostProcessGradientHDG()
#elif (PP_nVar==3)
  DO iVar=1, PP_nVar
    CALL PostProcessGradient(U_out(iVar,:,:),lambda(iVar,:,:),BTemp(iVar,:,:,:))
  END DO
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    r=k*(PP_N+1)**2+j*(PP_N+1) + i+1
    B(1,i,j,k,:) = BTemp(3,2,r,:) - BTemp(2,3,r,:)
    B(2,i,j,k,:) = BTemp(1,3,r,:) - BTemp(3,1,r,:)
    B(3,i,j,k,:) = BTemp(2,1,r,:) - BTemp(1,2,r,:)
  END DO; END DO; END DO !i,j,k
#else
  DO iVar=1, 3
    CALL PostProcessGradient(U_out(iVar,:,:),lambda(iVar,:,:),BTemp(iVar,:,:,:))
  END DO
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    r=k*(PP_N+1)**2+j*(PP_N+1) + i+1
    B(1,i,j,k,:) = BTemp(3,2,r,:) - BTemp(2,3,r,:)
    B(2,i,j,k,:) = BTemp(1,3,r,:) - BTemp(3,1,r,:)
    B(3,i,j,k,:) = BTemp(2,1,r,:) - BTemp(1,2,r,:)
  END DO; END DO; END DO !i,j,k
  CALL PostProcessGradient(U_out(4,:,:),lambda(4,:,:),E)
#endif
#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/

END SUBROUTINE HDGLinear
#endif /*USE_HDG*/


END MODULE MOD_HDG_Linear
