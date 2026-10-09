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
USE MOD_Interpolation_Vars ,ONLY: NMax,PREF_VDM,N_Inter
USE MOD_Elem_Mat           ,ONLY: PostProcessGradientHDG
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
USE MOD_Globals_Vars       ,ONLY: eps0
USE PETSc
USE MOD_Mesh_Vars          ,ONLY: SideToElem,nGlobalMortarSides
USE MOD_HDG_Vars_PETSc
#if defined(PARTICLES)
USE MOD_HDG_Readin         ,ONLY: UpdateChargeOnCMBC
#endif /*defined(PARTICLES*/
#if USE_MPI
USE MOD_MPI                ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
USE MOD_MPI_Vars
#endif
USE MOD_FillMortar_HDG     ,ONLY: BigToSmallMortar_HDG
#endif
#if USE_MPI
USE MOD_MPI_HDG            ,ONLY: Mask_MPIsides
#endif
USE MOD_Globals_Vars       ,ONLY: ElementaryCharge
USE MOD_ChangeBasis        ,ONLY: ChangeBasis2D
USE MOD_HDG_Tools          ,ONLY: CG_solver,DisplayConvergence
USE MOD_Interpolation_Vars ,ONLY: N_Inter
#if defined(PARTICLES)
USE MOD_Globals_Vars       ,ONLY: eps0
#if USE_MPI
USE MOD_PICDepo_MPI        ,ONLY: ExchangeSurfNodeSourceMPI
USE MOD_Particle_Boundary_Vars ,ONLY: Do2DSurfaceCharge
#endif /*USE_MPI*/
USE MOD_PICDepo_Vars       ,ONLY: SurfNodeSource,FEMVertexID2DepoSurfNodeID,Vdm_EQ_N,IsDepoSurfSide,pq2iNode,SurfNodeArea
USE MOD_Mesh_Vars          ,ONLY: NonUniqueGlobalSideIDToNonUniqueGlobalNodeID,NonUniqueGlobalNodeIDToFEMVertexID
USE MOD_Mesh_Vars          ,ONLY: SideToNonUniqueGlobalSide
#endif /*defined(PARTICLES)*/
#if defined(PARTICLES)
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound
#endif /*defined(PARTICLES)*/
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: time !time
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: lambdatmp(1:nGP_face(NMax))
INTEGER :: i,j,k,r,p,q,iElem, iVar,NSide,Nloc
INTEGER :: BCsideID,BCType,BCState,SideID,iLocSide
REAL    :: RHS_facetmp(nGP_face(NMax))
REAL    :: rtmp(nGP_vol(NMax))
#if defined(PARTICLES)
REAL                 :: src
INTEGER              :: FEMVertexID
INTEGER :: iPartBound, iNode, NonUniqueNodeID, NonUniqueGlobalSideID, iDepoSurfNodeID
REAL    :: SubArea
REAL,ALLOCATABLE     :: SurfNodeSourceEquiN1(:,:,:),SurfNodeSourceNodeTypeNloc(:,:,:)
#endif /*defined(PARTICLES)*/
#if (PP_nVar!=1)
REAL    :: BTemp(3,3,nGP_vol,PP_nElems)
#endif /*(PP_nVar!=1)*/
#if USE_LOADBALANCE
REAL    :: tLBStart
#endif /*USE_LOADBALANCE*/
#if USE_PETSC
INTEGER              :: DOFindices(nGP_face(NMax))
INTEGER              :: jNloc, iNloc
PetscErrorCode       :: ierr
PetscScalar, POINTER :: lambda_pointer(:)
KSPConvergedReason   :: reason
PetscInt             :: iterations
PetscReal            :: petscnorm
INTEGER              :: ElemID,iBCSide
INTEGER              :: DOF_start, DOF_stop
REAL                 :: timeStartPiclas,timeEndPiclas
INTEGER              :: jLocSide
REAL                 :: Smatloc(nGP_face(NMax),nGP_face(NMax))
INTEGER              :: iUniqueFPCBC
INTEGER              :: iMortar,iType
#endif /*USE_PETSC*/
REAL                 :: chitens_face(3,3),DCBiasVoltage(1:PP_nVar)
!===================================================================================================================================
! Dummy for chitens_face(:,:,p,q,SideID)
chitens_face=0.0
chitens_face(1,1)=1.
chitens_face(2,2)=1.
chitens_face(3,3)=1.
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
    Nloc = N_SurfMesh(SideID)%NSide
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
    CASE(53) ! exact BC = Dirichlet BC !! ExactFunc via RefState, bias voltage instead of phase shift with sine function
      DO q=0,Nloc; DO p=0,Nloc
        r=q*(Nloc+1) + p+1
        CALL ExactFunc(-53,N_SurfMesh(SideID)%Face_xGP(:,p,q),HDG_Surf_N(SideID)%lambda(iVar,r:r),t=time,iRefState=BCState)
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
    CASE(11) !neumann q*n=1 !test
      DO q=0,Nloc; DO p=0,Nloc
        r=q*(Nloc+1) + p+1
        HDG_Surf_N(SideID)%qn_face(iVar,r)=SUM((/1.,1.,1./)  &
                            *MATMUL(chitens_face(:,:),N_SurfMesh(SideID)%NormVec(:,p,q)))*N_SurfMesh(SideID)%SurfElem(p,q)*N_Inter(Nloc)%wGP(p)*N_Inter(Nloc)%wGP(q)
      END DO; END DO !p,q
    CASE(12) !neumann q*n=1 !test
      DO q=0,Nloc; DO p=0,Nloc
        r=q*(Nloc+1) + p+1
        HDG_Surf_N(SideID)%qn_face(iVar,r)=SUM((/-1.45e7,1.,1./)  &
                            *MATMUL(chitens_face(:,:),N_SurfMesh(SideID)%NormVec(:,p,q)))*N_SurfMesh(SideID)%SurfElem(p,q)*N_Inter(Nloc)%wGP(p)*N_Inter(Nloc)%wGP(q)
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
#endif /*PP_nVar!=1*/

  ! Floating boundary BCs
#if USE_PETSC
  IF(UseFPC) THEN
#if USE_MPI
    ! Communicate the accumulated charged on each BC to all processors on the communicator
    DO iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
      ASSOCIATE( COMM => FPC%COMM(iUniqueFPCBC)%UNICATOR)
        IF(COMM.NE.MPI_COMM_NULL)THEN
          CALL MPI_ALLREDUCE(MPI_IN_PLACE, FPC%ChargeProc(iUniqueFPCBC), 1, MPI_DOUBLE_PRECISION, MPI_SUM, COMM, IERROR)
          FPC%Charge(iUniqueFPCBC) = FPC%Charge(iUniqueFPCBC) + FPC%ChargeProc(iUniqueFPCBC)
        END IF ! COMM.NE.MPI_COMM_NULL
      END ASSOCIATE
    END DO ! iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
#else
    FPC%Charge = FPC%Charge + FPC%ChargeProc
#endif /*USE_MPI*/
    FPC%ChargeProc = 0.
    ! Apply charge to RHS, which is done below: RHS_conductor(1)=FPC%Charge(iUniqueFPCBC)/eps0
  END IF ! UseFPC

#if defined(PARTICLES)
  ! Communicate the accumulated charged on each BC to MPIRoot
  IF(UseCircuitModel) CALL UpdateChargeOnCMBC()
#endif /*defined(PARTICLES*/
#endif /*USE_PETSC*/

  ! Set potential to zero (only one process does this)
  IF(ZeroPotentialSide>0) HDG_Surf_N(ZeroPotentialSide)%lambda(iVar,1) = 0.
END DO ! iVar = 1, PP_nVar

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
    Nloc = N_SurfMesh(SideID)%NSide
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

      NSide = N_SurfMesh(SideID)%NSide
      IF(Nloc.NE.NSide)THEN
        CALL ChangeBasis2D(1, Nloc, NSide, TRANSPOSE(PREF_VDM(NSide,Nloc)%Vdm) , RHS_facetmp(1:nGP_face(Nloc)), RHS_facetmp(1:nGP_face(NSide)))
      END IF ! Nloc.NE.NSide
      HDG_Surf_N(SideID)%RHS_face(iVar,:) = HDG_Surf_N(SideID)%RHS_face(iVar,:) + RHS_facetmp(1:nGP_face(NSide))
    END DO
  END DO !iElem
END DO !ivar

!add Neumann
DO BCsideID=1,nNeumannBCSides
  SideID=NeumannBC(BCsideID)
  HDG_Surf_N(SideID)%RHS_face(:,:) = HDG_Surf_N(SideID)%RHS_face(:,:) + HDG_Surf_N(SideID)%qn_face(:,:)
END DO

#if defined(PARTICLES)
! Add Distributed Capacitance BC
#if USE_MPI
IF(Do2DSurfaceCharge) CALL ExchangeSurfNodeSourceMPI()
#endif /*USE_MPI*/
! NOTE: DistriCapBC only includes nBCSides and not the inner BCs
IF(nDistriCapBCsides.GT.0) ALLOCATE(SurfNodeSourceEquiN1(1:1,0:1,0:1),SurfNodeSourceNodeTypeNloc(1:1,0:Nmax,0:Nmax))
! Loop over all local BC sides, where the DCBC model is active
DO BCsideID=1,nDistriCapBCsides
  SideID     = DistriCapBC(BCsideID)             ! Get side index
  BCState    = BoundaryType(BC(SideID),BC_STATE)
  iPartBound = PartBound%MapToPartBC(BC(SideID)) ! Get particle boundary index
  Nloc       = N_SurfMesh(SideID)%NSide          ! Get polynomial degree of side
  NonUniqueGlobalSideID = SideToNonUniqueGlobalSide(1,SideID) ! Get global side index

  ! Map surface charge from vertices to SideID surface with N=1
  ! NOTE: Might not work because the side is not always oriented in the master ordering
  DO q=0,1; DO p=0,1
    ! Get local node index
    ! Use mapping p,q -> iNode
    iNode = pq2iNode(p,q,SideID)
    ! IPWRITE(*,*) 'p,q,iNode,2*q + p + 1:', p,q,iNode,2*q + p + 1
    ! Mapping from non-unique global side index to non-unique global node index
    NonUniqueNodeID = NonUniqueGlobalSideIDToNonUniqueGlobalNodeID(iNode,NonUniqueGlobalSideID)
    ! Sanity check
    IF (NonUniqueNodeID.LE.0) THEN
      IPWRITE(*,*) 'NonUniqueNodeID,NonUniqueGlobalSideID,IsDepoSurfSide(NonUniqueGlobalSideID),SideID:',&
                    NonUniqueNodeID,NonUniqueGlobalSideID,IsDepoSurfSide(NonUniqueGlobalSideID),SideID
      CALL abort(__STAMP__,' NonUniqueNodeID <= 0')
    END IF ! NonUniqueNodeID.LE.0
    ! Mapping from NonUniqueNodeID to FEMVertexID
    FEMVertexID = NonUniqueGlobalNodeIDToFEMVertexID(NonUniqueNodeID)
    ! Get surface deposition node index
    iDepoSurfNodeID = FEMVertexID2DepoSurfNodeID(FEMVertexID)
    ! Sanity check
    IF (iDepoSurfNodeID.LE.0) THEN
      IPWRITE(*,*) 'NonUniqueNodeID,FEMVertexID,iDepoSurfNodeID,NonUniqueGlobalSideID,IsDepoSurfSide(NonUniqueGlobalSideID),SideID:',&
                    NonUniqueNodeID,FEMVertexID,iDepoSurfNodeID,NonUniqueGlobalSideID,IsDepoSurfSide(NonUniqueGlobalSideID),SideID
      CALL abort(__STAMP__,' iDepoSurfNodeID <= 0')
    END IF ! iDepoSurfNodeID.LE.0
    ! Store in 2D temporary array
    SurfNodeSourceEquiN1(1:1,p,q) = SurfNodeSource(iDepoSurfNodeID)/SurfNodeArea(iDepoSurfNodeID)
  END DO; END DO ! q=0,1; DO p=0,1

  ! Map from equidistant side nodes (N=1) to side node type p-q system (Nloc)
  CALL ChangeBasis2D(1, 1, Nloc, Vdm_EQ_N(Nloc)%Vdm, SurfNodeSourceEquiN1(1:1,0:1,0:1), SurfNodeSourceNodeTypeNloc(1:1,0:Nloc,0:Nloc))

  ! Map from N=1 to N=Nloc
  DO q=0,Nloc; DO p=0,Nloc
    ! 2D p,q-index to 1D r-index
    r=q*(Nloc+1) + p+1
    ! Apply linear or constant electric potential
    IF (BCState.GT.0) THEN
      CALL ExactFunc(-3,N_SurfMesh(SideID)%Face_xGP(:,p,q),DCBiasVoltage(:),iLinState=BCState)
#if (PP_nVar!=1)
      CALL abort(__STAMP__,' LinPhi model not implemented for PP_nVar!=1', IERROR)
#endif /*PP_nVar!=1*/
    ELSE
      DCBiasVoltage = PartBound%DCBiasVoltage(iPartBound)
    END IF ! BCState.GT.0
    ! IPWRITE(*,*) 'SurfElem(p,q),SurfNodeSourceNodeTypeNloc(1,p,q),SurfNodeSourceNodeTypeNloc(1,p,q)/(N_SurfMesh(SideID)%SurfElem(p,q)):', &
    !               N_SurfMesh(SideID)%SurfElem(p,q),SurfNodeSourceNodeTypeNloc(1,p,q),SurfNodeSourceNodeTypeNloc(1,p,q)/(N_Inter(Nloc)%wGP(p)*N_Inter(Nloc)%wGP(q)*N_SurfMesh(SideID)%SurfElem(p,q))
    SubArea = N_Inter(Nloc)%wGP(p)*N_Inter(Nloc)%wGP(q)*N_SurfMesh(SideID)%SurfElem(p,q)
    src = SubArea * ( &
          PartBound%DCPermittivity(iPartBound) * DCBiasVoltage(1) / PartBound%DCThickness(iPartBound) + & ! DCBC
          (SurfNodeSourceNodeTypeNloc(1,p,q) +                         & ! Surface charge due to deposited particles
           PartBound%DCSurfaceChargeDensity(iPartBound) & ! Surface charge due to analytical expression
          )/eps0 )
    ! src = ( SubArea*PartBound%DCPermittivity(iPartBound) * PartBound%DCBiasVoltage(iPartBound) / PartBound%DCThickness(iPartBound) + & ! DCBC
    !       ( SurfNodeSourceNodeTypeNloc(1,p,q)        +                         & ! Surface charge due to deposited particles
    !         SubArea*PartBound%DCSurfaceChargeDensity(iPartBound) & ! Surface charge due to analytical expression
    !       )/eps0 )
    HDG_Surf_N(SideID)%RHS_face(1,r) = HDG_Surf_N(SideID)%RHS_face(1,r) + src
  END DO; END DO ! p,q
  ! read*
END DO
#else
IF(nDistriCapBCsides.GT.0) CALL Abort(__STAMP__,'ERROR: Distributed capacitance requires PARTICLES=ON')
#endif /*defined(PARTICLES)*/

#if USE_PETSC
! add Dirichlet contribution
DO iBCSide=1,nDirichletBCSides
  BCSideID=DirichletBC(iBCSide)
  ElemID    = SideToElem(S2E_ELEM_ID,BCSideID)
  jLocSide  = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
  ! At BCSides, jNLoc = Nmax
  jNloc     = N_SurfMesh(BCsideID)%NSide
  DO iLocSide=1,6
    SideID = ElemToSide(E2S_SIDE_ID,iLocSide,ElemID)
    iNloc     = N_SurfMesh(SideID)%NSide
    IF(MaskedSide(SideID).GT.0) CYCLE

    ! TODO PETSC P-Adaption - Improvement: Store V^T * S * V in Smat
    Smatloc(1:nGP_face(jNloc),1:nGP_face(jNloc)) = HDG_Vol_N(ElemID)%Smat(:,:,iLocSide,jLocSide)
    ! 1. S_{Ij} = (V^T)_{Ii} * S_{ij}
    IF(jNloc.NE.iNloc)THEN
      DO j=1,nGP_face(jNloc)
        CALL ChangeBasis2D(1, jNloc, iNloc, TRANSPOSE(PREF_VDM(iNloc,jNloc)%Vdm), Smatloc(1:nGP_face(jNloc),j), Smatloc(1:nGP_face(iNloc),j))
      END DO
    END IF

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

#if USE_MPI && !USE_PETSC
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
TimeStartPiclas=PICLASTIME()
! Fill right hand side
PetscCallA(VecZeroEntries(PETScRHS,ierr))
DO SideID=1,nSides
  IF(MaskedSide(SideID).GT.0) CYCLE

  Nloc = N_SurfMesh(SideID)%NSide
  DO i=1,nGP_face(Nloc)
    DOFindices(i) = i + OffsetGlobalPETScDOF(SideID) - 1
  END DO

  PetscCallA(VecSetValues(PETScRHS,nGP_face(Nloc),DOFindices(1:nGP_face(Nloc)),HDG_Surf_N(SideID)%RHS_face(1,:),ADD_VALUES,ierr))
END DO

! 1st assembly call is required because "Calls to VecSetValues() with the INSERT_VALUES and ADD_VALUES options cannot be mixed
! without intervening calls to the assembly routines."
PetscCallA(VecAssemblyBegin(PETScRHS,ierr))
PetscCallA(VecAssemblyEnd(PETScRHS,ierr))

! The MPIRoot process has charge and voltage of all FPCs, therefore, this process sets all conductor RHS information
IF(UseFPC) THEN
  IF(MPIRoot)THEN
    DO iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
      PetscCallA(VecSetValues(PETScRHS,1,[nGlobalPETScDOFs-1-FPC%nUniqueFPCBounds+iUniqueFPCBC],[FPC%Charge(iUniqueFPCBC)/eps0],INSERT_VALUES,ierr))
    END DO
  END IF
END IF

#if defined(PARTICLES)
! The MPIRoot process has charge and voltage of the CMBC, therefore, this process sets all RHS information
IF(UseCircuitModel) THEN
  IF(MPIRoot)THEN
    CALL ExactFunc(-1,(/0.,0.,0./),CMBC%VoltageRF(1:1),t=time,iRefState=CMBC%RefState)
    src = ( CMBC%Charge + CMBC%Capacitance*CMBC%VoltageRF(1) ) / eps0
    ! IPWRITE(*,*) 'hdg_linear.f90 ',__LINE__,'     C=',CMBC%Capacitance, 'Phi_rf=',CMBC%VoltageRF(1),'Phi_a=',CMBC%Voltage,'( C*Phi_rf + Q ) / eps0=',src,"Q=",CMBC%Charge
    PetscCallA(VecSetValues(PETScRHS,1,[nGlobalPETScDOFs-1],[src],INSERT_VALUES,ierr))
  END IF ! MPIRoot
END IF ! UseCircuitModel
#endif /*defined(PARTICLES*/

! Reset the RHS of the first DOF if ZeroPotential must be set
IF(mpiRoot.AND.ZeroPotentialDOF.GE.0) THEN
  ! NOTE: sice PETSc 3.24, VecSetValue checks the data type of the arguments. The value 0 is not allowed, but 0.0 is correct
  PetscCallA(VecSetValue(PETScRHS,ZeroPotentialDOF,0.0,INSERT_VALUES,ierr))
END IF

! 2nd assembly call is required because "These values may be cached, so VecAssemblyBegin() and VecAssemblyEnd() MUST be called after
! all calls to VecSetValues() have been completed."
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

#if ((PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR > 22)) || (PETSC_VERSION_MAJOR >= 4)
! PETSc 3.26.0 and onwards requires special TYPE() for variable "reason", which cannot be compared with Fortran INT
#if ((PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 24))
! 3.24: Deprecate KSP_CONVERGED_RTOL_NORMAL in favor of KSP_CONVERGED_RTOL_NORMAL_EQUATIONS and
!                 KSP_CONVERGED_ATOL_NORMAL in favor of KSP_CONVERGED_ATOL_NORMAL_EQUATIONS
IF(.NOT.(reason.EQ.KSP_CONVERGED_RTOL_NORMAL_EQUATIONS     .OR.&
         reason.EQ.KSP_CONVERGED_ATOL_NORMAL_EQUATIONS     .OR.&
#else
IF(.NOT.(reason.EQ.KSP_CONVERGED_RTOL_NORMAL               .OR.&
         reason.EQ.KSP_CONVERGED_ATOL_NORMAL               .OR.&
#endif
         reason.EQ.KSP_CONVERGED_RTOL                      .OR.&
         reason.EQ.KSP_CONVERGED_ATOL                      .OR.&
         reason.EQ.KSP_CONVERGED_ITS                       .OR.&
         reason.EQ.KSP_CONVERGED_NEG_CURVE                 .OR.&
         reason.EQ.KSP_CONVERGED_CG_NEG_CURVE_DEPRECATED   .OR.&
         reason.EQ.KSP_CONVERGED_CG_CONSTRAINED_DEPRECATED .OR.&
         reason.EQ.KSP_CONVERGED_STEP_LENGTH               .OR.&
         reason.EQ.KSP_CONVERGED_HAPPY_BREAKDOWN)) THEN
#else
IF(reason.LT.0)THEN
#endif
  ! Output used memory
  CALL WarningMemusage(Mode=1,Threshold=5.0) ! Mode=1: Memory per node (NOT over all nodes), Threshold=5.0% of total node memory
  !  View solver converged reason
  PetscCallA(KSPConvergedReasonView(PETScSolver,PETSC_VIEWER_STDOUT_WORLD,ierr))
  !  View solver info
  PetscCallA(KSPView(PETScSolver,PETSC_VIEWER_STDOUT_WORLD,ierr))
  IPWRITE(*,*) 'PETSc convergence failed with reason:', reason
  CALL Abort(__STAMP__,'ERROR: PETSc not converged')
END IF ! reason.LT.0 (i.e. not converged)

! Only the MPIRoot determines the computation time required for PETSc.
IF(MPIroot) THEN
  ! NOTE: this measurement includes the time spent in the MPI routines of PETSc, which is not included in time measurements of
  ! piclas modules to distinguish between calculation time and MPI time. Only calculation time is relevant for load balancing, but
  ! PETSc does its own balancing method so this information here is just for comparison with the time required for the particles.
  PETScFieldTime = TimeEndPiclas-TimeStartPiclas
  CALL DisplayConvergence(PETScFieldTime, iterations, petscnorm)
END IF ! MPIroot

! Fill element local lambda for post processing
! Get the local DOF subarray
PetscCallA(VecScatterBegin(PETScScatter, PETScSolution, PETScSolutionLocal, INSERT_VALUES, SCATTER_FORWARD,ierr))
PetscCallA(VecScatterEnd(  PETScScatter, PETScSolution, PETScSolutionLocal, INSERT_VALUES, SCATTER_FORWARD,ierr))
PetscCallA(VecGetArrayRead(PETScSolutionLocal,lambda_pointer,ierr))
DOF_stop = 0
DO SideID=1,nSides
  IF(MaskedSide(SideID).GT.0) CYCLE
  Nloc = N_SurfMesh(SideID)%NSide
  DOF_start = 1 + DOF_stop
  DOF_stop = DOF_start + nGP_face(Nloc) - 1
  IF(nGlobalMortarSides.LE.0)THEN ! No Mortars in mesh
    HDG_Surf_N(SideID)%lambda(1,:) = lambda_pointer(DOF_start:DOF_stop)
  ELSE ! Mortars are present
    IF(SmallMortarInfo(SideID).EQ.0) THEN
      HDG_Surf_N(SideID)%lambda(1,:) = lambda_pointer(DOF_start:DOF_stop)
    ELSE
      ! lambda_small = M * lambda_big
      iType = SmallMortarType(1,SideID)
      iMortar = SmallMortarType(2,SideID)
      HDG_Surf_N(SideID)%lambda(1,:) = MATMUL(N_Inter(Nloc)%IntMatMortar(:,:,iMortar,iType),lambda_pointer(DOF_start:DOF_stop))
    END IF
  END IF ! nGlobalMortarSides.GT.0
END DO

! Fill Conductor lambda
IF(UseFPC) THEN
  FPC%VoltageProc = 0. ! nullify just to be safe
  DO BCsideID=1,nConductorBCsides
    SideID       = ConductorBC(BCSideID)
    Nloc         = N_SurfMesh(SideID)%NSide
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

#if defined(PARTICLES)
! Fill circuit model lambda
IF(UseCircuitModel) THEN
  CMBC%Voltage = 0. ! Nullify just to be safe
  DO BCsideID = 1,nCircuitModelBCsides
    SideID  = CircuitModelBC(BCSideID)
    Nloc    = N_SurfMesh(SideID)%NSide
    BCState = BoundaryType(BC(SideID),BC_STATE)
    DO i=1,nGP_face(Nloc)
      HDG_Surf_N(SideID)%lambda(1,i) = lambda_pointer(nLocalPETScDOFs)
    END DO
  END DO
  ! MPIRoot sets global value for BC and I/O
  IF(MPIRoot) CMBC%Voltage = lambda_pointer(nLocalPETScDOFs)
  ! IPWRITE(*,*) 'hdg_linear.f90      CMBC%Voltage:', CMBC%Voltage
END IF ! UseCircuitModel
#endif /*defined(PARTICLES*/

PetscCallA(VecRestoreArrayRead(PETScSolutionLocal,lambda_pointer,ierr))
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
      NSide = N_SurfMesh(SideID)%NSide
      IF(Nloc.EQ.NSide)THEN
        lambdatmp(1:nGP_face(Nloc)) = HDG_Surf_N(SideID)%lambda(iVar,:)
      ELSE
        CALL ChangeBasis2D(1, NSide, Nloc, PREF_VDM(NSide,Nloc)%Vdm , HDG_Surf_N(SideID)%lambda(iVar,1:nGP_face(NSide)), lambdatmp(1:nGP_face(Nloc)))
      END IF ! Nloc.EQ.NSide
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