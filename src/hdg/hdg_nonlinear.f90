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
MODULE MOD_HDG_NonLinear
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
#if USE_HDG
#if defined(PARTICLES)
PUBLIC :: HDGNewton
#endif /*defined(PARTICLES)*/
#endif /*USE_HDG*/
!===================================================================================================================================

CONTAINS

#if USE_HDG
!===================================================================================================================================
!> HDG non-linear solver via Newton's method
!===================================================================================================================================
#if defined(PARTICLES)
!SUBROUTINE HDGNewton(time,U_out,td_iter,ForceCGSolverIteration_opt)
SUBROUTINE HDGNewton(time,td_iter,ForceCGSolverIteration_opt)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_HDG_Vars
USE MOD_Equation               ,ONLY: CalcSourceHDG,ExactFunc
USE MOD_FillMortar_HDG         ,ONLY: SmallToBigMortar_HDG
USE MOD_TimeDisc_Vars          ,ONLY: IterDisplayStep,DoDisplayIter
USE MOD_Globals_Vars           ,ONLY: eps0
USE MOD_Mesh_Vars              ,ONLY: N_SurfMesh,BoundaryType,nSides,BC,offSetElem
USE MOD_Mesh_Vars              ,ONLY: ElemToSide
USE MOD_Interpolation_Vars     ,ONLY: NMax,PREF_VDM
USE MOD_Elem_Mat               ,ONLY: PostProcessGradientHDG, Elem_Mat,BuildPrecond
USE MOD_Restart_Vars           ,ONLY: DoRestart,RestartTime
#if (PP_nVar==1)
!USE MOD_Equation_Vars          ,ONLY: E
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers     ,ONLY: LBStartTime,LBSplitTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
USE MOD_HDG_Vars               ,ONLY: ElemToBRRegion,RegionElectronRef
USE MOD_DG_Vars                ,ONLY: U_N,N_DG_Mapping
#if USE_MPI
USE MOD_MPI                    ,ONLY: Mask_MPIsides
#endif
USE MOD_HDG_Tools              ,ONLY: CG_solver,DisplayConvergence
USE MOD_ChangeBasis            ,ONLY: ChangeBasis2D
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: time !time
INTEGER(KIND=8),INTENT(IN)  :: td_iter
#if defined(PARTICLES)
LOGICAL,INTENT(IN),OPTIONAL :: ForceCGSolverIteration_opt ! set converged=F in first step (only BR electron fluid)
#endif /*defined(PARTICLES)*/
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!REAL,INTENT(INOUT)  :: U_out(PP_nVar,nGP_vol(PP_N),PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER :: iVar=1
REAL    :: lambdatmp(1:nGP_face(NMax))
REAL    :: RHS_facetmp(nGP_face(NMax))
INTEGER :: i,j,k,r,p,q,iElem, iter,RegionID,Nloc,NSide
INTEGER :: BCsideID,BCType,BCState,SideID,iLocSide
!REAL    :: RHS_face(PP_nVar,nGP_face(PP_N),nSides)
REAL    :: rtmp(nGP_vol(PP_N)),Norm_r2!,Norm_r2_old
LOGICAL :: converged, beLinear
LOGICAL :: warning_linear
REAL    :: warning_linear_phi
#if (PP_nVar!=1)
REAL    :: BTemp(3,3,nGP_vol(PP_N),PP_nElems)
#endif
#if USE_LOADBALANCE
REAL    :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
#if (PP_nVar!=1)
  CALL CollectiveStop(__STAMP__,'Nonlinear Newton solver only available with EQ-system Poisson!')
#endif
Norm_r2=0.

!Dirichlet boundary conditions
DO BCsideID=1,nDirichletBCSides
  SideID=DirichletBC(BCsideID)
  BCType =BoundaryType(BC(SideID),BC_TYPE)
  BCState=BoundaryType(BC(SideID),BC_STATE)
  SELECT CASE(BCType)
  CASE(2) ! exact BC = Dirichlet BC !! ExactFunc via BCState (time is optional)
    DO q=0,PP_N; DO p=0,PP_N
      r=q*(PP_N+1) + p+1
      CALL ExactFunc(BCState,N_SurfMesh(SideID)%Face_xGP(:,p,q),HDG_Surf_N(SideID)%lambda(PP_nVar,r:r),time)
    END DO; END DO !p,q
  CASE(4) ! exact BC = Dirichlet BC !! Zero potential
    DO q=0,PP_N; DO p=0,PP_N
      r=q*(PP_N+1) + p+1
      HDG_Surf_N(SideID)%lambda(PP_nVar,r:r)= 0.
    END DO; END DO !p,q
  CASE(5) ! exact BC = Dirichlet BC !! ExactFunc via RefState (time is optional) for reference state (with zero crossing)
    DO q=0,PP_N; DO p=0,PP_N
      r=q*(PP_N+1) + p+1
      CALL ExactFunc(-1,N_SurfMesh(SideID)%Face_xGP(:,p,q),HDG_Surf_N(SideID)%lambda(PP_nVar,r:r),t=time,iRefState=BCState)
    END DO; END DO !p,q
  CASE(6) ! exact BC = Dirichlet BC !! ExactFunc via RefState (time is optional) for reference state (without zero crossing)
    DO q=0,PP_N; DO p=0,PP_N
      r=q*(PP_N+1) + p+1
      CALL ExactFunc(-2,N_SurfMesh(SideID)%Face_xGP(:,p,q),HDG_Surf_N(SideID)%lambda(PP_nVar,r:r),t=time,iRefState=BCState)
    END DO; END DO !p,q
  CASE(7) ! exact BC = Dirichlet BC !! ExactFunc via LinState (time is optional)for linear potential function
    DO q=0,PP_N; DO p=0,PP_N
      r=q*(PP_N+1) + p+1
      CALL ExactFunc(-3,N_SurfMesh(SideID)%Face_xGP(:,p,q),HDG_Surf_N(SideID)%lambda(PP_nVar,r:r),t=time,iLinState=BCState)
    END DO; END DO !p,q
    CASE(8,50,51,52,60) ! exact BC = Dirichlet BC !! ExactFunc via electric potential and decharing of a surface
      CALL abort(__STAMP__,'Dirichlet BC=8,50,51,52,60 model not implemented for HDG Newton!')
  END SELECT ! BCType
END DO !BCsideID=1,nDirichletBCSides

!neumann BC
DO BCsideID=1,nNeumannBCSides
  SideID  = NeumannBC(BCsideID)
  BCType  = BoundaryType(BC(SideID),BC_TYPE)
  BCState = BoundaryType(BC(SideID),BC_STATE)
  SELECT CASE(BCType)
  CASE(10) !neumann q=0 !! Zero gradient
    DO q=0,PP_N; DO p=0,PP_N
      r=q*(PP_N+1) + p+1
      HDG_Surf_N(SideID)%qn_face(PP_nVar,r)= 0.
    END DO; END DO !p,q
  END SELECT ! BCType
END DO !BCsideID=1,nNeumannBCSides

warning_linear=.FALSE.
DO iElem=1,PP_nElems
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
    r=k*(Nloc+1)**2+j*(Nloc+1) + i+1
    !CALL CalcSourceHDG(i,j,k,iElem,HDG_Vol_N(iElem)%RHS_vol(1:PP_nVar,r),U_out(1,r,iElem),warning_linear,warning_linear_phi)
    CALL CalcSourceHDG(i,j,k,iElem,HDG_Vol_N(iElem)%RHS_vol(1:PP_nVar,r),U_N(iElem)%U(1,i,j,k),warning_linear,warning_linear_phi)
  END DO; END DO; END DO !i,j,k
  HDG_Vol_N(iElem)%RHS_Vol(iVar,:)=-HDG_Vol_N(iElem)%JwGP_vol(:)*HDG_Vol_N(iElem)%RHS_vol(iVar,:)
END DO !iElem
IF (warning_linear) THEN
  SWRITE(*,'(A,ES10.2E3)') 'WARNING: during iteration at least one DOF resulted in a phi > phi_max.\n'//&
    '=> Increase Part-RegionElectronRef#-PhiMax if already steady! Phi-Phi_ref=',warning_linear_phi
END IF

!prepare RHS_face ( RHS for lambda system.)
DO iElem=1,PP_nElems
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
    r=k*(Nloc+1)**2+j*(Nloc+1) + i+1
    HDG_Vol_N(iElem)%RHS_vol(PP_nVar,r)=HDG_Vol_N(iElem)%RHS_vol(PP_nVar,r)&
                                       +HDG_Vol_N(iElem)%JwGP_vol(r)*U_N(iElem)%U(1,i,j,k)*HDG_Vol_N(iElem)%NonlinVolumeFac(r)
  END DO; END DO; END DO !i,j,k
END DO

!RHS_face(PP_nVar,:,:) =0.
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
    ! TODO NSideMin - LOW/HIGH
    IF(Nloc.EQ.NSide)THEN
      HDG_Surf_N(SideID)%RHS_face(iVar,:) = HDG_Surf_N(SideID)%RHS_face(iVar,:) + RHS_facetmp(1:nGP_face(Nloc))
    ELSE
      CALL ChangeBasis2D(1, Nloc, NSide, TRANSPOSE(PREF_VDM(NSide,Nloc)%Vdm) , RHS_facetmp(1:nGP_face(Nloc)), RHS_facetmp(1:nGP_face(NSide)))
      HDG_Surf_N(SideID)%RHS_face(iVar,:) = HDG_Surf_N(SideID)%RHS_face(iVar,:) + RHS_facetmp(1:nGP_face(NSide))
    END IF ! Nloc.EQ.NSide
  END DO
END DO !iElem


DO BCsideID=1,nNeumannBCSides
  SideID=NeumannBC(BCsideID)
  HDG_Surf_N(SideID)%RHS_face(:,:)=HDG_Surf_N(SideID)%RHS_face(:,:)+HDG_Surf_N(SideID)%qn_face(:,:)
END DO

#if USE_LOADBALANCE
CALL LBSplitTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/
#if USE_MPI
!CALL Mask_MPIsides(PP_nVar,RHS_Face)
CALL Mask_MPIsides('RHS_face',1)
#endif /*USE_MPI*/
CALL SmallToBigMortar_HDG(PP_nVar,0)!RHS_face(1:PP_nVar,1:nGP_Face,1:nSides))

#if USE_LOADBALANCE
CALL LBSplitTime(LB_DGCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/


! SOLVE
!CALL CheckNonLinRes(RHS_face(1,:,:),lambda(1,:,:),converged,Norm_r2)
CALL CheckNonLinRes(converged,Norm_r2)
IF(PRESENT(ForceCGSolverIteration_opt))THEN
  IF(ForceCGSolverIteration_opt)THEN
    ! Due to the removal of electrons during restart
    converged=.false.
    SWRITE(UNIT_StdOut,*) "Forcing initial CG solver to iterate by setting converged=F (Norm_r2=",Norm_r2,")"
  END IF ! ForceCGSolverIteration_opt
END IF ! ForceCGSolverIteration_opt
IF (converged) THEN
  SWRITE(*,*) 'HDGNewton: Newton Iteration has converged in 0 steps...'
ELSE
  !CALL CG_solver(RHS_face(PP_nVar,:,:),lambda(PP_nVar,:,:))
  !post processing:
  CALL CG_solver(1)

  DO iElem=1,PP_nElems
    Nloc = N_DG_Mapping(2,iElem+offSetElem)
    ! for post-proc
    DO iLocSide=1,6
      SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
      NSide = N_SurfMesh(SideID)%NSide
      ! TODO NSideMin - LOW/HIGH
      IF(Nloc.EQ.NSide)THEN
        lambdatmp(1:nGP_face(Nloc)) = HDG_Surf_N(SideID)%lambda(iVar,:)
      ELSE
        ! From low to high
        CALL ChangeBasis2D(1, NSide, Nloc, PREF_VDM(NSide,Nloc)%Vdm , HDG_Surf_N(SideID)%lambda(iVar,1:nGP_face(NSide)), lambdatmp(1:nGP_face(Nloc)))
      END IF ! Nloc.EQ.NSide
      CALL DGEMV('T',nGP_face(Nloc),nGP_vol(Nloc),1., &
                          HDG_Vol_N(iElem)%Ehat(:,:,iLocSide), nGP_face(Nloc), &
                          lambdatmp(1:nGP_face(Nloc)),1,1.,& !add to RHS_face
                          HDG_Vol_N(iElem)%RHS_vol(PP_nVar,:),1)
    END DO
    CALL DSYMV('U',nGP_vol(Nloc),1., HDG_Vol_N(iElem)%InvDhat(:,:),nGP_vol(Nloc), &
                               -HDG_Vol_N(iElem)%RHS_vol(PP_nVar,:),1,0., &
                               U_N(iElem)%U(1,:,:,:),1)
  END DO !iElem

  IF(NewtonAdaptStartValue) THEN
    IF ((.NOT.DoRestart.AND.ALMOSTEQUAL(time,0.)).OR.(DoRestart.AND.ALMOSTEQUAL(time,RestartTime))) THEN
      DO iElem=1,PP_nElems
        RegionID=ElemToBRRegion(iElem)
        Nloc = N_DG_Mapping(2,iElem+offSetElem)
        DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
          r=k*(Nloc+1)**2+j*(Nloc+1) + i+1
          IF (NewtonExactSourceDeriv) THEN
            !HDG_Vol_N(iElem)%NonlinVolumeFac(r) = RegionElectronRef(1,RegionID)/ (RegionElectronRef(3,RegionID)*eps0) &
                         !* EXP( (U_out(1,r,iElem)-RegionElectronRef(2,RegionID)) / RegionElectronRef(3,RegionID) )
            HDG_Vol_N(iElem)%NonlinVolumeFac(r) = RegionElectronRef(1,RegionID)/ (RegionElectronRef(3,RegionID)*eps0) &
                         * EXP( (U_N(iElem)%U(1,i,j,k)-RegionElectronRef(2,RegionID)) / RegionElectronRef(3,RegionID) )
          ELSE
            HDG_Vol_N(iElem)%NonlinVolumeFac(r)=RegionElectronRef(1,RegionID) / (RegionElectronRef(3,RegionID)*eps0)
          END IF
        END DO; END DO; END DO !i,j,k
      END DO !iElem
      CALL Elem_Mat(td_iter)
      CALL BuildPrecond()
    END IF
  END IF

  converged =.false.
  beLinear=.false.
  AdaptIterNewton = AdaptIterNewtonOld
  DO iter=1,MaxIterNewton

    IF (.NOT.beLinear) THEN
      IF ((iter.EQ.AdaptIterNewtonToLinear)) THEN !.OR.(iter.GT.3*AdaptIterNewtonOld)) THEN
                                                 !removed second cond. to ensure fast convergence with very small AdaptIterNewton
        IF(MPIroot) WRITE(*,*) 'Info: Solver not converging with exact source derivative, switching to linearization (The linear way, baby)'
        DO iElem=1,PP_nElems
          RegionID=ElemToBRRegion(iElem)
          Nloc = N_DG_Mapping(2,iElem+offSetElem)
          DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
            r=k*(Nloc+1)**2+j*(Nloc+1) + i+1
            HDG_Vol_N(iElem)%NonlinVolumeFac(r)=RegionElectronRef(1,RegionID) / (RegionElectronRef(3,RegionID)*eps0)
          END DO; END DO; END DO !i,j,k
        END DO !iElem
        CALL Elem_Mat(td_iter)
        CALL BuildPrecond()
        AdaptIterNewton = 0
        beLinear=.true.
      END IF
    END IF

    IF (AdaptIterNewton.GT.0) THEN
      IF (MOD(iter,AdaptIterNewton).EQ.0) THEN
        DO iElem=1,PP_nElems
          Nloc = N_DG_Mapping(2,iElem+offSetElem)
          RegionID=ElemToBRRegion(iElem)
          DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
            r=k*(Nloc+1)**2+j*(Nloc+1) + i+1
            IF (NewtonExactSourceDeriv) THEN
              !HDG_Vol_N(iElem)%NonlinVolumeFac(r) = RegionElectronRef(1,RegionID)/ (RegionElectronRef(3,RegionID)*eps0) &
                         !* EXP( (U_out(1,r,iElem)-RegionElectronRef(2,RegionID)) / RegionElectronRef(3,RegionID) )
              HDG_Vol_N(iElem)%NonlinVolumeFac(r) = RegionElectronRef(1,RegionID)/ (RegionElectronRef(3,RegionID)*eps0) &
                         * EXP( (U_N(iElem)%U(1,i,j,k)-RegionElectronRef(2,RegionID)) / RegionElectronRef(3,RegionID) )
            ELSE
              HDG_Vol_N(iElem)%NonlinVolumeFac(r)=RegionElectronRef(1,RegionID) / (RegionElectronRef(3,RegionID)*eps0)
            END IF
          END DO; END DO; END DO !i,j,k
        END DO !iElem
        CALL Elem_Mat(td_iter)
        CALL BuildPrecond()
      END IF
    END IF
    !volume source (volume RHS of u system)
    !SWRITE(*,*) '!!!!!!!!!!!!!!!!!', iter

    warning_linear=.FALSE.
    DO iElem=1,PP_nElems
      Nloc = N_DG_Mapping(2,iElem+offSetElem)
      DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
        r=k*(Nloc+1)**2+j*(Nloc+1) + i+1
        !CALL CalcSourceHDG(i,j,k,iElem,HDG_Vol_N(iElem)%RHS_vol(1:PP_nVar,r),U_out(1,r,iElem),warning_linear,warning_linear_phi)
        CALL CalcSourceHDG(i,j,k,iElem,HDG_Vol_N(iElem)%RHS_vol(1:PP_nVar,r),U_N(iElem)%U(1,i,j,k),warning_linear,warning_linear_phi)
      END DO; END DO; END DO !i,j,k
      HDG_Vol_N(iElem)%RHS_Vol(iVar,:)=-HDG_Vol_N(iElem)%JwGP_vol(:)*HDG_Vol_N(iElem)%RHS_vol(iVar,:)
    END DO !iElem
    IF (warning_linear) THEN
      !SWRITE(*,'(A,F5.2,A,F5.2,A)') 'HDGNewton WARNING: during iteration at least one DOF resulted in a phi > phi_max. '//&
        !'=> Increase Part-RegionElectronRef#-PhiMax if already steady! Phi-Phi_ref=',warning_linear_phi,'(Phi_ref=',RegionElectronRef(2,RegionID),')'
      SWRITE(*,'(A,ES10.2E3)') 'HDGNewton WARNING: at least one DOF resulted in phi > phi_max. '//&
        'Increase Part-RegionElectronRef#-PhiMax to shift the ref. point! Phi-Phi_ref=',warning_linear_phi!,' (Phi_ref=',RegionElectronRef(2,RegionID),')'
    END IF

    !prepare RHS_face ( RHS for lamdba system.)
    DO iElem=1,PP_nElems
      Nloc = N_DG_Mapping(2,iElem+offSetElem)
      DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
        r=k*(Nloc+1)**2+j*(Nloc+1) + i+1
        HDG_Vol_N(iElem)%RHS_vol(iVar,r) = HDG_Vol_N(iElem)%RHS_vol(iVar,r)&
                                         + HDG_Vol_N(iElem)%JwGP_vol(r)*U_N(iElem)%U(iVar,i,j,k)*HDG_Vol_N(iElem)%NonlinVolumeFac(r)
      END DO; END DO; END DO !i,j,k
    END DO !iElem

    !RHS_face(PP_nVar,:,:) =0.
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
      ! TODO NSideMin - LOW/HIGH
      IF(Nloc.EQ.NSide)THEN
        HDG_Surf_N(SideID)%RHS_face(iVar,:) = HDG_Surf_N(SideID)%RHS_face(iVar,:) + RHS_facetmp(1:nGP_face(Nloc))
      ELSE
        ! From high to low
        CALL ChangeBasis2D(1, Nloc, NSide, TRANSPOSE(PREF_VDM(NSide,Nloc)%Vdm) , RHS_facetmp(1:nGP_face(Nloc)), RHS_facetmp(1:nGP_face(NSide)))
        HDG_Surf_N(SideID)%RHS_face(iVar,:) = HDG_Surf_N(SideID)%RHS_face(iVar,:) + RHS_facetmp(1:nGP_face(NSide))
      END IF ! Nloc.EQ.NSide
      END DO
    END DO !iElem

    !add Neumann
    DO BCsideID=1,nNeumannBCSides
      SideID=NeumannBC(BCsideID)
      !RHS_face(:,:,SideID) = RHS_face(:,:,SideID) + HDG_Surf_N(SideID)%qn_face(:,:)
      HDG_Surf_N(SideID)%RHS_face(:,:) = HDG_Surf_N(SideID)%RHS_face(:,:) + HDG_Surf_N(SideID)%qn_face(:,:)
    END DO


#if USE_MPI
  CALL Mask_MPIsides('RHS_face',1)
#endif /*USE_MPI*/
  !CALL SmallToBigMortar_HDG(PP_nVar,HDG_Surf_N(1)%RHS_face(1:PP_nVar,1:nGP_Face)) ! 1 -> 1:nSIdes
  CALL SmallToBigMortar_HDG(1,0) ! RHS_face(1:PP_nVar,1:nGP_Face,1:nSides))

    ! SOLVE
    !CALL CheckNonLinRes(HDG_Surf_N(1)%RHS_face(1,:),lambda(1,:,:),converged,Norm_r2)
    CALL CheckNonLinRes(converged,Norm_r2)
    IF (converged) THEN
      IF(DoDisplayIter)THEN
        IF(HDGDisplayConvergence.AND.(MOD(td_iter,IterDisplayStep).EQ.0)) THEN
          SWRITE(*,*) 'HDGNewton: Newton Iteration has converged in ',iter,' steps...'
        END IF
      END IF
      EXIT
    ELSE IF (iter.EQ.MaxIterNewton) THEN
      IPWRITE(UNIT_StdOut,*) "Norm_r2       =", Norm_r2
      IPWRITE(UNIT_StdOut,*) "iter          =", iter
      IPWRITE(UNIT_StdOut,*) "MaxIterNewton =", MaxIterNewton
      CALL abort(__STAMP__,'HDGNewton: Newton Iteration has NOT converged!')
    END IF

    CALL CG_solver(1)!RHS_face(PP_nVar,:,:),lambda(PP_nVar,:,:))

    !post processing:
    DO iElem=1,PP_nElems
      Nloc = N_DG_Mapping(2,iElem+offSetElem)
      ! for post-proc
      DO iLocSide=1,6
        SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
        NSide = N_SurfMesh(SideID)%NSide
        ! TODO NSideMin - LOW/HIGH
        IF(Nloc.EQ.NSide)THEN
          lambdatmp(1:nGP_face(Nloc)) = HDG_Surf_N(SideID)%lambda(iVar,:)
        ELSE
          ! From low to high
          CALL ChangeBasis2D(1, NSide, Nloc, PREF_VDM(NSide,Nloc)%Vdm , HDG_Surf_N(SideID)%lambda(iVar,1:nGP_face(NSide)), lambdatmp(1:nGP_face(Nloc)))
        END IF ! Nloc.EQ.NSide
        CALL DGEMV('T',nGP_face(Nloc),nGP_vol(Nloc),1., &
                            HDG_Vol_N(iElem)%Ehat(:,:,iLocSide), nGP_face(Nloc), &
                            lambdatmp(1:nGP_face(Nloc)),1,1.,& ! 1: add to RHS_face, 0: set value
                            HDG_Vol_N(iElem)%RHS_vol(PP_nVar,:),1)
      END DO ! iLocSide=1,6
      CALL DSYMV('U',nGP_vol(Nloc),1., HDG_Vol_N(iElem)%InvDhat(:,:),nGP_vol(Nloc), &
                                 -HDG_Vol_N(iElem)%RHS_vol(PP_nVar,:),1,0., &
                                 U_N(iElem)%U(PP_nVar,:,:,:),1)
    END DO ! iElem=1,PP_nElems
  END DO ! iter=1,MaxIterNewto
END IF ! converged

#if (PP_nVar==1)
CALL PostProcessGradientHDG()
#endif

#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/
END SUBROUTINE HDGNewton


!===================================================================================================================================
!> Determine the residual of the HDG solution
!===================================================================================================================================
!SUBROUTINE CheckNonLinRes(RHS,lambda,converged,Norm_R2)
SUBROUTINE CheckNonLinRes(converged,Norm_R2)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_HDG_Vars  ,ONLY: nGP_face
USE MOD_HDG_Vars  ,ONLY: EpsNonLinear
USE MOD_Mesh_Vars ,ONLY: nSides
#if USE_MPI
USE MOD_Mesh_Vars ,ONLY: nMPISides_YOUR
#endif /*USE_MPI*/
#if defined(MEASURE_MPI_WAIT)
USE MOD_MPI_Vars  ,ONLY: MPIW8TimeField,MPIW8CountField
#endif /*defined(MEASURE_MPI_WAIT)*/
USE MOD_HDG_Tools ,ONLY: evalresidual,VectorDotProductRR
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!REAL, INTENT(IN)    :: RHS(nGP_face(PP_N)*nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!REAL, INTENT(INOUT)    :: lambda(nGP_face(PP_N)*nSides)
LOGICAL, INTENT(INOUT) :: converged
REAL, INTENT(OUT)      :: Norm_r2
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL,DIMENSION(nGP_face(PP_N)*nSides) :: R
INTEGER                         :: VecSize
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=8)   :: CounterStart,CounterEnd
REAL(KIND=8)      :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================
#if USE_MPI
! not use MPI_YOUR sides for vector_dot_product!!!
VecSize=(nSides-nMPIsides_YOUR)*nGP_face(PP_N)
#else
VecSize=nSides*nGP_face(PP_N)
#endif /*USE_MPI*/
CALL EvalResidual(1) !RHS,lambda,R)

!CALL VectorDotProduct(VecSize,R(1:VecSize),R(1:VecSize),Norm_R2) !Z=V
CALL VectorDotProductRR(Norm_R2) !Z=V
!  print*, Norm_R2

#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/

#if USE_MPI
IF(MPIroot) converged=(Norm_R2.LT.EpsNonLinear**2)
CALL MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_PICLAS,iError)
#else
converged=(Norm_R2.LT.EpsNonLinear**2)
#endif /*USE_MPI*/

#if defined(MEASURE_MPI_WAIT)
CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
MPIW8TimeField(3)  = MPIW8TimeField(3) + REAL(CounterEnd-CounterStart,8)/Rate
MPIW8CountField(3) = MPIW8CountField(3) + 1_8
#endif /*defined(MEASURE_MPI_WAIT)*/
END SUBROUTINE CheckNonLinRes
#endif /*defined(PARTICLES)*/




#endif /*USE_HDG*/


END MODULE MOD_HDG_NonLinear
