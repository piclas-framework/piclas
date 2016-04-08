#include "boltzplatz.h"

MODULE MOD_Eval_xyz
!===================================================================================================================================
! Changes a 3D Tensor Product Lagrange Points of Lagrange Basis of degree N_In to  
! Lagrange points of a Lagrange Basis for one point, using two
! arbitrary point disributions xi_In(0:N_In) and xi_Out 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------

!#ifdef PARTICLES
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE eval_xyz_curved
  MODULE PROCEDURE eval_xyz_curved
END INTERFACE

INTERFACE eval_xyz_elemcheck
  MODULE PROCEDURE eval_xyz_elemcheck
END INTERFACE

INTERFACE eval_xyz_part2
  MODULE PROCEDURE eval_xyz_part2
END INTERFACE

INTERFACE eval_xyz_poly
  MODULE PROCEDURE eval_xyz_poly
END INTERFACE


PUBLIC :: eval_xyz_curved,eval_xyz_elemcheck, eval_xyz_part2,eval_xyz_poly
!#endif /*PARTICLES*/
!===================================================================================================================================

CONTAINS

!#ifdef PARTICLES
SUBROUTINE eval_xyz_curved(x_in,NVar,N_in,U_In,U_Out,ElemID,PartID)
!===================================================================================================================================
! interpolate a 3D tensor product Lagrange basis defined by (N_in+1) 1D interpolation point positions x
! first get xi,eta,zeta from x,y,z...then do tenso product interpolation
! xi is defined in the 1DrefElem xi=[-1,1]
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Basis,                   ONLY:LagrangeInterpolationPolys
USE MOD_Interpolation_Vars,      ONLY:xGP,wBary
USE MOD_Mesh_Vars,               ONLY:dXCL_NGeo,Elem_xGP,XCL_NGeo,NGeo,wBaryCL_NGeo,XiCL_NGeo
USE MOD_Particle_Mesh_Vars,      ONLY:RefMappingGuess
USE MOD_Particle_Mesh_Vars,      ONLY:XiEtaZetaBasis,ElemBaryNGeo,slenXiEtaZetaBasis
USE MOD_PICInterpolation_Vars,   ONLY:NBG,BGField,useBGField,BGDataSize,BGField_wBary, BGField_xGP,BGType
USE MOD_Mesh_Vars,               ONLY:CurvedElem,wBaryCL_NGeo1,XiCL_NGeo1
!USE MOD_Particle_Vars,           ONLY:PartPosRef
!USE MOD_Mesh_Vars,ONLY: X_CP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: NVar                                  ! 6 (Ex, Ey, Ez, Bx, By, Bz) 
INTEGER,INTENT(IN)  :: N_In                                  ! usually PP_N
INTEGER,INTENT(IN)  :: ElemID                                 ! elem index
REAL,INTENT(IN)     :: U_In(1:NVar,0:N_In,0:N_In,0:N_In)   ! elem state
REAL,INTENT(IN)     :: x_in(3)                                  ! physical position of particle 
INTEGER,INTENT(IN),OPTIONAL :: PartID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: U_Out(1:NVar)  ! Interpolated state
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER             :: i,j,k
REAL                :: xi(3)
!REAL                :: X3D_Buf1(1:NVar,0:N_In,0:N_In)  ! first intermediate results from 1D interpolations
!REAL                :: X3D_Buf2(1:NVar,0:N_In) ! second intermediate results from 1D interpolations
REAL                :: Winner_Dist,Dist!,abortcrit
REAL, PARAMETER     :: EPSONE=1.00000001
INTEGER             :: iDir
!REAL                :: Lag(1:3,0:NGeo)
REAL                :: L_xi(3,0:PP_N), L_eta_zeta
REAL                :: Ptild(1:3),XiLinear(1:6)
REAL                :: XCL_NGeo1(1:3,0:1,0:1,0:1)
REAL                :: dXCL_NGeo1(1:3,1:3,0:1,0:1,0:1)
REAL                :: XiA, XiB
LOGICAL             :: Found
! h5-external e,b field
REAL,ALLOCATABLE    :: L_xi_BGField(:,:), U_BGField(:)
!===================================================================================================================================

! get initial guess by nearest GP search ! simple guess
!IF(CurvedElem(ElemID))THEN
!  Xi=0.
!  !XiA=-1.05
!  !XiB=1.05
!  !CALL RefElemBisection(Xi(1),XiA,XiB,X_In(1),wBaryCL_NGeo,XiCL_NGeo,XCL_NGeo(1,:,:,:,ElemID),NGeo,Found)
!  !XiA=-1.05
!  !XiB=1.05
!  !CALL RefElemBisection(Xi(2),XiA,XiB,X_In(2),wBaryCL_NGeo,XiCL_NGeo,XCL_NGeo(2,:,:,:,ElemID),NGeo,Found)
!  !XiA=-1.05
!  !XiB=1.05
!  !CALL RefElemBisection(Xi(3),XiA,XiB,X_In(3),wBaryCL_NGeo,XiCL_NGeo,XCL_NGeo(3,:,:,:,ElemID),NGeo,Found)
!ELSE
  SELECT CASE(RefMappingGuess)
  CASE(1)
    Ptild=X_in - ElemBaryNGeo(:,ElemID)
    ! plus coord system (1-3) and minus coord system (4-6)
    DO iDir=1,6
      XiLinear(iDir)=DOT_PRODUCT(Ptild,XiEtaZetaBasis(:,iDir,ElemID))*slenXiEtaZetaBasis(iDir,ElemID)
    END DO
    ! compute guess as average value
    DO iDir=1,3
      Xi(iDir)=0.5*(XiLinear(iDir)-XiLinear(iDir+3))
    END DO 
    !IF(MAXVAL(ABS(Xi)).GT.epsOne) Xi=0.
    IF(MAXVAL(ABS(Xi)).GT.epsOne) Xi=LimitXi(Xi)
  CASE(2) 
    ! compute distance on Gauss Points
    Winner_Dist=HUGE(1.)
    DO i=0,N_in; DO j=0,N_in; DO k=0,N_in
      Dist=SUM((x_in(:)-Elem_xGP(:,i,j,k,ElemID))*(x_in(:)-Elem_xGP(:,i,j,k,ElemID)))
      IF (Dist.LT.Winner_Dist) THEN
        Winner_Dist=Dist
        Xi(:)=(/xGP(i),xGP(j),xGP(k)/) ! start value
      END IF
    END DO; END DO; END DO
  CASE(3) 
    ! compute distance on XCL Points
    Winner_Dist=HUGE(1.)
    DO i=0,NGeo; DO j=0,NGeo; DO k=0,NGeo
      Dist=SUM((x_in(:)-XCL_NGeo(:,i,j,k,ElemID))*(x_in(:)-XCL_NGeo(:,i,j,k,ElemID)))
      IF (Dist.LT.Winner_Dist) THEN
        Winner_Dist=Dist
        Xi(:)=(/XiCL_NGeo(i),XiCL_NGeo(j),XiCL_NGeo(k)/) ! start value
      END IF
    END DO; END DO; END DO
  CASE(4)
    ! trival guess 
    xi=0.
  END SELECT
!END IF

IF(CurvedElem(ElemID))THEN
  CALL RefElemNewton(Xi,X_In,wBaryCL_NGeo,XiCL_NGeo,XCL_NGeo(:,:,:,:,ElemID),dXCL_NGeo(:,:,:,:,:,ElemID),NGeo,ElemID,Mode=1)
ELSE
  ! fill dummy XCL_NGeo1
  XCL_NGeo1(1:3,0,0,0) = XCL_NGeo(1:3, 0  , 0  , 0  ,ElemID)
  XCL_NGeo1(1:3,1,0,0) = XCL_NGeo(1:3,NGeo, 0  , 0  ,ElemID)
  XCL_NGeo1(1:3,0,1,0) = XCL_NGeo(1:3, 0  ,NGeo, 0  ,ElemID)
  XCL_NGeo1(1:3,1,1,0) = XCL_NGeo(1:3,NGeo,NGeo, 0  ,ElemID)
  XCL_NGeo1(1:3,0,0,1) = XCL_NGeo(1:3, 0  , 0  ,NGeo,ElemID)
  XCL_NGeo1(1:3,1,0,1) = XCL_NGeo(1:3,NGeo, 0  ,NGeo,ElemID)
  XCL_NGeo1(1:3,0,1,1) = XCL_NGeo(1:3, 0  ,NGeo,NGeo,ElemID)
  XCL_NGeo1(1:3,1,1,1) = XCL_NGeo(1:3,NGeo,NGeo,NGeo,ElemID)
  ! fill dummy dXCL_NGeo1
  dXCL_NGeo1(1:3,1:3,0,0,0) = dXCL_NGeo(1:3,1:3, 0  , 0  , 0  ,ElemID)
  dXCL_NGeo1(1:3,1:3,1,0,0) = dXCL_NGeo(1:3,1:3,NGeo, 0  , 0  ,ElemID)
  dXCL_NGeo1(1:3,1:3,0,1,0) = dXCL_NGeo(1:3,1:3, 0  ,NGeo, 0  ,ElemID)
  dXCL_NGeo1(1:3,1:3,1,1,0) = dXCL_NGeo(1:3,1:3,NGeo,NGeo, 0  ,ElemID)
  dXCL_NGeo1(1:3,1:3,0,0,1) = dXCL_NGeo(1:3,1:3, 0  , 0  ,NGeo,ElemID)
  dXCL_NGeo1(1:3,1:3,1,0,1) = dXCL_NGeo(1:3,1:3,NGeo, 0  ,NGeo,ElemID)
  dXCL_NGeo1(1:3,1:3,0,1,1) = dXCL_NGeo(1:3,1:3, 0  ,NGeo,NGeo,ElemID)
  dXCL_NGeo1(1:3,1:3,1,1,1) = dXCL_NGeo(1:3,1:3,NGeo,NGeo,NGeo,ElemID)
  CALL RefElemNewton(Xi,X_In,wBaryCL_NGeo1,XiCL_NGeo1,XCL_NGeo1,dXCL_NGeo1,1,ElemID,Mode=1)
END IF

! IF(ANY(ABS(Xi).GT.1.1)) THEN
! !IF((NewtonIter.GE.4).AND.(ANY(ABS(Xi).GT.1.5)))THEN
!   IPWRITE(UNIT_stdOut,*) ' Particle not inside of element, force!!!'
!   IPWRITE(UNIT_stdOut,*) ' Newton-Iter', NewtonIter
!   IPWRITE(UNIT_stdOut,*) ' xi  ', xi(1)
!   IPWRITE(UNIT_stdOut,*) ' eta ', xi(2)
!   IPWRITE(UNIT_stdOut,*) ' zeta', xi(3)
!   IPWRITE(UNIT_stdOut,*) ' PartPos', X_in
!   IPWRITE(UNIT_stdOut,*) ' ElemID', ElemID+offSetElem
!   IF(PRESENT(PartID)) IPWRITE(UNIT_stdOut,*) ' PartID', PartID
!   CALL abort(__STAMP__, &
!       'Particle Not inSide of Element, ElemID, iPart',ElemID,REAL(PartID))
! END IF

! 2.1) get "Vandermonde" vectors
CALL LagrangeInterpolationPolys(xi(1),N_in,xGP,wBary,L_xi(1,:))
CALL LagrangeInterpolationPolys(xi(2),N_in,xGP,wBary,L_xi(2,:))
CALL LagrangeInterpolationPolys(xi(3),N_in,xGP,wBary,L_xi(3,:))

! "more efficient" - Quote Thomas B.
U_out(:)=0
DO k=0,N_in
  DO j=0,N_in
    L_eta_zeta=L_xi(2,j)*L_xi(3,k)
    DO i=0,N_in
      U_out = U_out + U_IN(:,i,j,k)*L_xi(1,i)*L_Eta_Zeta
    END DO ! i=0,N_In
  END DO ! j=0,N_In
END DO ! k=0,N_In

IF(useBGField)THEN
  ! use of BG-Field with possible different polynomial order and nodetype
  ALLOCATE( L_xi_BGField(3,0:NBG)             &
          , U_BGField(1:BGDataSize)           )
!          , X3D_tmp1(BGDataSize,0:NBG,0:NBG) &
!          , X3D_tmp2(BGDataSize,0:NBG)       &
!          , X3D_tmp3(BGDataSize)             )
  CALL LagrangeInterpolationPolys(xi(1),NBG,BGField_xGP,BGField_wBary,L_xi_BGField(1,:))
  CALL LagrangeInterpolationPolys(xi(2),NBG,BGField_xGP,BGField_wBary,L_xi_BGField(2,:))
  CALL LagrangeInterpolationPolys(xi(3),NBG,BGField_xGP,BGField_wBary,L_xi_BGField(3,:))
  
  U_BGField(:)=0
  DO k=0,NBG
    DO j=0,NBG
      L_eta_zeta=L_xi_BGField(2,j)*L_xi_BGField(3,k)
      DO i=0,NBG
        U_BGField = U_BGField + BGField(:,i,j,k,ElemID)*L_xi_BGField(1,i)*L_Eta_Zeta
      END DO ! i=0,NBG
    END DO ! j=0,NBG
  END DO ! k=0,NBG

  SELECT CASE(BGType)
  CASE(1)
    U_Out(1:3)=U_Out(1:3)+U_BGField
  CASE(2)
    U_Out(4:6)=U_Out(4:6)+U_BGField
  CASE(3)
    U_Out=U_Out+U_BGField
  END SELECT
  DEALLOCATE( L_xi_BGField, U_BGField)! X3d_tmp1, x3d_tmp2, x3d_tmp3)
END IF ! useBGField

END SUBROUTINE eval_xyz_curved


SUBROUTINE eval_xyz_elemcheck(x_in,xi,ElemID,PartID,DoReUseMap)
!===================================================================================================================================
! interpolate a 3D tensor product Lagrange basis defined by (N_in+1) 1D interpolation point positions x
! first get xi,eta,zeta from x,y,z...then do tenso product interpolation
! xi is defined in the 1DrefElem xi=[-1,1]
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Basis,                   ONLY:LagrangeInterpolationPolys
USE MOD_Interpolation_Vars,      ONLY:xGP
USE MOD_Particle_Mesh_Vars,      ONLY:RefMappingGuess,RefMappingEps
USE MOD_Particle_Mesh_Vars,      ONLY:XiEtaZetaBasis,ElemBaryNGeo,slenXiEtaZetaBasis!,ElemRadiusNGeo
USE MOD_Mesh_Vars,               ONLY:dXCL_NGeo,Elem_xGP,XCL_NGeo,NGeo,wBaryCL_NGeo,XiCL_NGeo,NGeo
USE MOD_Mesh_Vars,               ONLY:CurvedElem,wBaryCL_NGeo1,XiCL_NGeo1
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)          :: ElemID                                 ! elem index
REAL,INTENT(IN)             :: x_in(3)                                  ! physical position of particle 
INTEGER,INTENT(IN),OPTIONAL :: PartID
LOGICAL,INTENT(IN),OPTIONAL :: DoReUseMap
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)          :: xi(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                    :: i,j,k
REAL                       :: epsOne
REAL                       :: Winner_Dist,Dist
INTEGER                    :: idir
REAL                       :: Ptild(1:3),XiLinear(1:6)
REAL                       :: XCL_NGeo1(1:3,0:1,0:1,0:1)
REAL                       :: dXCL_NGeo1(1:3,1:3,0:1,0:1,0:1)
REAL                       :: XiA,XiB
!===================================================================================================================================


epsOne=1.0+RefMappingEps
IF(.NOT.PRESENT(DoReUseMap))THEN
  !IF(CurvedElem(ElemID))THEN
  !  Xi=0.
  !  !print*,'partpois',x_in
  !  !XiA=-1.05
  !  !XiB=1.05
  !  !print*,'x'
  !  !CALL RefElemBisection(Xi(1),XiA,XiB,X_In(1),wBaryCL_NGeo,XiCL_NGeo,XCL_NGeo(1,:,:,:,ElemID),NGeo,Found)
  !  !XiA=-1.05
  !  !XiB=1.05
  !  !print*,'y'
  !  !CALL RefElemBisection(Xi(2),XiA,XiB,X_In(2),wBaryCL_NGeo,XiCL_NGeo,XCL_NGeo(2,:,:,:,ElemID),NGeo,Found)
  !  !XiA=-1.05
  !  !XiB=1.05
  !  !print*,'z'
  !  !CALL RefElemBisection(Xi(3),XiA,XiB,X_In(3),wBaryCL_NGeo,XiCL_NGeo,XCL_NGeo(3,:,:,:,ElemID),NGeo,Found)
  !  !print*,'guess',xi
  !ELSE
    SELECT CASE(RefMappingGuess)
    CASE(1)
      Ptild=X_in - ElemBaryNGeo(:,ElemID)
      ! plus coord system (1-3) and minus coord system (4-6)
      DO iDir=1,6
        XiLinear(iDir)=DOT_PRODUCT(Ptild,XiEtaZetaBasis(:,iDir,ElemID))*slenXiEtaZetaBasis(iDir,ElemID)
      END DO
      ! compute guess as average value
      DO iDir=1,3
        Xi(iDir)=0.5*(XiLinear(iDir)-XiLinear(iDir+3))
      END DO 
      IF(MAXVAL(ABS(Xi)).GT.epsOne) Xi=LimitXi(Xi)
    CASE(2)
      Winner_Dist=HUGE(1.)
      DO i=0,PP_N; DO j=0,PP_N; DO k=0,PP_N
        Dist=SUM((x_in(:)-Elem_xGP(:,i,j,k,ElemID))*(x_in(:)-Elem_xGP(:,i,j,k,ElemID)))
        IF (Dist.LT.Winner_Dist) THEN
          Winner_Dist=Dist
          Xi(:)=(/xGP(i),xGP(j),xGP(k)/) ! start value
        END IF
      END DO; END DO; END DO
    CASE(3) 
      ! compute distance on XCL Points
      Winner_Dist=HUGE(1.)
      DO i=0,NGeo; DO j=0,NGeo; DO k=0,NGeo
        Dist=SUM((x_in(:)-XCL_NGeo(:,i,j,k,ElemID))*(x_in(:)-XCL_NGeo(:,i,j,k,ElemID)))
        IF (Dist.LT.Winner_Dist) THEN
          Winner_Dist=Dist
          Xi(:)=(/XiCL_NGeo(i),XiCL_NGeo(j),XiCL_NGeo(k)/) ! start value
        END IF
      END DO; END DO; END DO
    CASE(4)
      ! trival guess, cell mean point
      xi=0.
    END SELECT
  !END IF
END IF

IF(CurvedElem(ElemID))THEN
  CALL RefElemNewton(Xi,X_In,wBaryCL_NGeo,XiCL_NGeo,XCL_NGeo(:,:,:,:,ElemID),dXCL_NGeo(:,:,:,:,:,ElemID),NGeo,ElemID,Mode=2)
ELSE
  ! fill dummy XCL_NGeo1
  IF(NGeo.EQ.1)THEN
    CALL RefElemNewton(Xi,X_In,wBaryCL_NGeo,XiCL_NGeo,XCL_NGeo(:,:,:,:,ElemID),dXCL_NGeo(:,:,:,:,:,ElemID),NGeo,ElemID,Mode=2)
  ELSE
    XCL_NGeo1(1:3,0,0,0) = XCL_NGeo(1:3, 0  , 0  , 0  ,ElemID)
    XCL_NGeo1(1:3,1,0,0) = XCL_NGeo(1:3,NGeo, 0  , 0  ,ElemID)
    XCL_NGeo1(1:3,0,1,0) = XCL_NGeo(1:3, 0  ,NGeo, 0  ,ElemID)
    XCL_NGeo1(1:3,1,1,0) = XCL_NGeo(1:3,NGeo,NGeo, 0  ,ElemID)
    XCL_NGeo1(1:3,0,0,1) = XCL_NGeo(1:3, 0  , 0  ,NGeo,ElemID)
    XCL_NGeo1(1:3,1,0,1) = XCL_NGeo(1:3,NGeo, 0  ,NGeo,ElemID)
    XCL_NGeo1(1:3,0,1,1) = XCL_NGeo(1:3, 0  ,NGeo,NGeo,ElemID)
    XCL_NGeo1(1:3,1,1,1) = XCL_NGeo(1:3,NGeo,NGeo,NGeo,ElemID)
    ! fill dummy dXCL_NGeo1
    dXCL_NGeo1(1:3,1:3,0,0,0) = dXCL_NGeo(1:3,1:3, 0  , 0  , 0  ,ElemID)
    dXCL_NGeo1(1:3,1:3,1,0,0) = dXCL_NGeo(1:3,1:3,NGeo, 0  , 0  ,ElemID)
    dXCL_NGeo1(1:3,1:3,0,1,0) = dXCL_NGeo(1:3,1:3, 0  ,NGeo, 0  ,ElemID)
    dXCL_NGeo1(1:3,1:3,1,1,0) = dXCL_NGeo(1:3,1:3,NGeo,NGeo, 0  ,ElemID)
    dXCL_NGeo1(1:3,1:3,0,0,1) = dXCL_NGeo(1:3,1:3, 0  , 0  ,NGeo,ElemID)
    dXCL_NGeo1(1:3,1:3,1,0,1) = dXCL_NGeo(1:3,1:3,NGeo, 0  ,NGeo,ElemID)
    dXCL_NGeo1(1:3,1:3,0,1,1) = dXCL_NGeo(1:3,1:3, 0  ,NGeo,NGeo,ElemID)
    dXCL_NGeo1(1:3,1:3,1,1,1) = dXCL_NGeo(1:3,1:3,NGeo,NGeo,NGeo,ElemID)
    CALL RefElemNewton(Xi,X_In,wBaryCL_NGeo1,XiCL_NGeo1,XCL_NGeo1,dXCL_NGeo1,1,ElemID,Mode=2)
  END IF
END IF

!print*,'found xi',xi
!read*

END SUBROUTINE eval_xyz_elemcheck


SUBROUTINE Eval_xyz_Poly(Xi_in,NVar,N_in,xGP_in,wBary_In,U_In,U_Out)
!===================================================================================================================================
! interpolate a 3D tensor product Lagrange basis defined by (N_in+1) 1D interpolation point positions x
!===================================================================================================================================
! MODULES
USE MOD_Basis,                 ONLY: LagrangeInterpolationPolys
!USE MOD_Interpolation_Vars,    ONLY: wBary,xGP
!USE MOD_Mesh_Vars,             ONLY: wBaryCL_NGeo,XiCL_NGeo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)        :: NVar, N_in
REAL,INTENT(IN)           :: U_In(1:NVar,0:N_In,0:N_In,0:N_In)   ! solution
REAL,INTENT(IN)           :: xi_in(3)                            ! reference space position of particle 
REAL,INTENT(IN)           :: xGP_In(0:N_in)
REAL,INTENT(IN)           :: wBary_In(0:N_in)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)          :: U_Out(1:NVar)  ! Interpolated state
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                   :: i,j,k
REAL,DIMENSION(3,0:N_in)  :: L_xi        
REAL                      :: L_eta_zeta
!===================================================================================================================================

! 
CALL LagrangeInterpolationPolys(xi_in(1),N_in,xGP_in,wBary_In,L_xi(1,:))
CALL LagrangeInterpolationPolys(xi_in(2),N_in,xGP_in,wBary_In,L_xi(2,:))
CALL LagrangeInterpolationPolys(xi_in(3),N_in,xGP_in,wBary_In,L_xi(3,:))

U_out(:)=0
DO k=0,N_in
  DO j=0,N_in
    L_eta_zeta=L_xi(2,j)*L_xi(3,k)
    DO i=0,N_in
      U_out = U_out + U_IN(:,i,j,k)*L_xi(1,i)*L_eta_zeta
    END DO ! i=0,N_In
  END DO ! j=0,N_In
END DO ! k=0,N_In


END SUBROUTINE Eval_xyz_poly


SUBROUTINE eval_xyz_part2(xi_in,NVar,N_in,U_In,U_Out,ElemID)
!===================================================================================================================================
! interpolate a 3D tensor product Lagrange basis defined by (N_in+1) 1D interpolation point positions x
! hoewver, particle is already mapped to reference space -1|1
! xi instead of x_in
!===================================================================================================================================
! MODULES
USE MOD_Basis,                 ONLY: LagrangeInterpolationPolys
USE MOD_Interpolation_Vars,    ONLY: wBary,xGP
USE MOD_PICInterpolation_Vars, ONLY:NBG,BGField,useBGField,BGDataSize,BGField_xGP,BGField_wBary,BGType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)        :: NVar                                  ! 6 (Ex, Ey, Ez, Bx, By, Bz) 
INTEGER,INTENT(IN)        :: N_In                                  ! usually PP_N
INTEGER,INTENT(IN)        :: ElemID                                 ! elem index
REAL,INTENT(IN)           :: U_In(1:NVar,0:N_In,0:N_In,0:N_In)   ! elem state
REAL,INTENT(IN)           :: xi_in(3)                              ! reference space position of particle 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)          :: U_Out(1:NVar)  ! Interpolated state
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER             :: i,j,k
!REAL                :: X3D_Buf1(1:NVar,0:N_In,0:N_In)  ! first intermediate results from 1D interpolations
!REAL                :: X3D_Buf2(1:NVar,0:N_In) ! second intermediate results from 1D interpolations
REAL                :: L_xi(3,0:N_in), L_eta_zeta
!REAL                :: buff,buff2
! h5-external e,b field
REAL,ALLOCATABLE    :: L_xi_BGField(:,:), U_BGField(:)
!===================================================================================================================================

! 2.1) get "Vandermonde" vectors
CALL LagrangeInterpolationPolys(xi_in(1),N_in,xGP,wBary,L_xi(1,:))
CALL LagrangeInterpolationPolys(xi_in(2),N_in,xGP,wBary,L_xi(2,:))
CALL LagrangeInterpolationPolys(xi_in(3),N_in,xGP,wBary,L_xi(3,:))


! "more efficient" - Quote Thomas B.
U_out(:)=0
DO k=0,N_in
  DO j=0,N_in
    L_eta_zeta=L_xi(2,j)*L_xi(3,k)
    DO i=0,N_in
      U_out = U_out + U_IN(:,i,j,k)*L_xi(1,i)*L_Eta_Zeta
    END DO ! i=0,N_In
  END DO ! j=0,N_In
END DO ! k=0,N_In

!! 2.2) do the tensor product thing 
!X3D_buf1=0.
!! first direction iN_In
!DO k=0,N_In
!  DO j=0,N_In
!    DO i=0,N_In
!      X3D_Buf1(:,j,k)=X3D_Buf1(:,j,k)+Lag2(1,i)*X3D_In(:,i,j,k)
!    END DO
!  END DO
!END DO
!X3D_buf2=0.
!! second direction jN_In
!DO k=0,N_In
!  DO j=0,N_In
!    X3D_Buf2(:,k)=X3D_Buf2(:,k)+Lag2(2,j)*X3D_Buf1(:,j,k)
!  END DO
!END DO
!X3D_Out=0.
!! last direction kN_In
!DO k=0,N_In
!  X3D_Out(:)=X3D_Out(:)+Lag2(3,k)*X3D_Buf2(:,k)
!END DO

IF(useBGField)THEN
  ! use of BG-Field with possible different polynomial order and nodetype
  ALLOCATE( L_xi_BGField(3,0:NBG)             &
          , U_BGField(1:BGDataSize)           )
!          , X3D_tmp1(BGDataSize,0:NBG,0:NBG) &
!          , X3D_tmp2(BGDataSize,0:NBG)       &
!          , X3D_tmp3(BGDataSize)             )
  CALL LagrangeInterpolationPolys(xi_in(1),NBG,BGField_xGP,BGField_wBary,L_xi_BGField(1,:))
  CALL LagrangeInterpolationPolys(xi_in(2),NBG,BGField_xGP,BGField_wBary,L_xi_BGField(2,:))
  CALL LagrangeInterpolationPolys(xi_in(3),NBG,BGField_xGP,BGField_wBary,L_xi_BGField(3,:))
  
  U_BGField(:)=0
  DO k=0,NBG
    DO j=0,NBG
      L_eta_zeta=L_xi_BGField(2,j)*L_xi_BGField(3,k)
      DO i=0,NBG
        U_BGField = U_BGField + BGField(:,i,j,k,ElemID)*L_xi_BGField(1,i)*L_Eta_Zeta
      END DO ! i=0,NBG
    END DO ! j=0,NBG
  END DO ! k=0,NBG


  !! 2.2) do the tensor product thing 
  !X3D_tmp1=0.
  !! first direction iN_In
  !DO k=0,NBG
  !  DO j=0,NBG
  !    DO i=0,NBG
  !      X3D_tmp1(:,j,k)=X3D_tmp1(:,j,k)+L_xi_BGField(1,i)*BGField(:,i,j,k,ElemID)
  !    END DO
  !  END DO
  !END DO
  !X3D_tmp2=0.
  !! second direction jN_In
  !DO k=0,NBG
  !  DO j=0,NBG
  !    X3D_tmp2(:,k)=X3D_tmp2(:,k)+L_xi_BGField(2,j)*X3D_tmp1(:,j,k)
  !  END DO
  !END DO
  !X3D_tmp3=0.
  !! last direction kN_In
  !DO k=0,NBG
  !  X3D_tmp3(:)=X3D_tmp3(:)+L_xi_BGField(3,k)*X3D_tmp2(:,k)
  !END DO
  SELECT CASE(BGType)
  CASE(1)
    U_Out(1:3)=U_Out(1:3)+U_BGField
  CASE(2)
    U_Out(4:6)=U_Out(4:6)+U_BGField
  CASE(3)
    U_Out=U_Out+U_BGField
  END SELECT
  DEALLOCATE( L_xi_BGField, U_BGFIeld)! X3d_tmp1, x3d_tmp2, x3d_tmp3)
END IF ! useBGField


END SUBROUTINE eval_xyz_part2


SUBROUTINE RefElemNewton(Xi,X_In,wBaryCL_N_In,XiCL_N_In,XCL_N_In,dXCL_N_In,N_In,ElemID,Mode)
!=================================================================================================================================
! Netwon for finding the position inside the reference element [-1,1] for an arbitrary physical point
!=================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_Basis,                   ONLY:LagrangeInterpolationPolys
USE MOD_Particle_Mesh_Vars,      ONLY:RefMappingEps!,ElemRadiusNGeo
USE MOD_Mesh_Vars,               ONLY:offsetElem
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
INTEGER,INTENT(IN)               :: N_In,ElemID
INTEGER,INTENT(IN)               :: Mode
REAL,INTENT(IN)                  :: X_in(3) ! position in physical space 
REAL,INTENT(IN)                  :: XiCL_N_in(0:N_In)               ! position of CL points in reference space
REAL,INTENT(IN)                  ::  XCL_N_in(3,0:N_In,0:N_in,0:N_In) ! position of CL points in physical space
REAL,INTENT(IN)                  :: dXCL_N_in(3,3,0:N_In,0:N_in,0:N_In) ! derivation of CL points
REAL,INTENT(IN)                  :: wBaryCL_N_in(0:N_In) ! derivation of CL points
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(INOUT)               :: Xi(3) ! position in reference element
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: Lag(1:3,0:N_In), F(1:3)
INTEGER                          :: NewTonIter,i,j,k
REAL                             :: deltaXi(1:3),deltaXi2
REAL                             :: Jac(1:3,1:3),sdetJac,sJac(1:3,1:3)
REAL                             :: buff,buff2
!===================================================================================================================================


! initial guess
CALL LagrangeInterpolationPolys(Xi(1),N_In,XiCL_N_in,wBaryCL_N_in,Lag(1,:))
CALL LagrangeInterpolationPolys(Xi(2),N_In,XiCL_N_in,wBaryCL_N_in,Lag(2,:))
CALL LagrangeInterpolationPolys(Xi(3),N_In,XiCL_N_in,wBaryCL_N_in,Lag(3,:))
! F(xi) = x(xi) - x_in
F=-x_in ! xRp
DO k=0,N_In
  DO j=0,N_In
    buff=Lag(2,j)*Lag(3,k)
    DO i=0,N_In
      F=F+XCL_N_in(:,i,j,k)*Lag(1,i)*buff !Lag(2,j)*Lag(3,k)
    END DO !l=0,N_In
  END DO !i=0,N_In
END DO !j=0,N_In

IF(ALL(ABS(F).LT.epsMach)) THEN
  deltaXi2=0.
ELSE
  deltaXi2=1. !HUGE(1.0)
END IF

NewtonIter=0
!abortCrit=ElemRadiusN_in(ElemID)*ElemRadiusN_in(ElemID)*RefMappingEps
DO WHILE((deltaXi2.GT.RefMappingEps).AND.(NewtonIter.LT.100))
  NewtonIter=NewtonIter+1

  ! caution, dXCL_NGeo is transposed of required matrix
  Jac=0.
  DO k=0,N_In
    DO j=0,N_In
      buff=Lag(2,j)*Lag(3,k)
      DO i=0,N_In
        buff2=Lag(1,i)*buff
        Jac(1,1:3)=Jac(1,1:3)+dXCL_N_in(1:3,1,i,j,k)*buff2
        Jac(2,1:3)=Jac(2,1:3)+dXCL_N_in(1:3,2,i,j,k)*buff2
        Jac(3,1:3)=Jac(3,1:3)+dXCL_N_in(1:3,3,i,j,k)*buff2
      END DO !i=0,N_In
    END DO !j=0,N_In
  END DO !k=0,N_In
  
  ! Compute inverse of Jacobian
  sdetJac=getDet(Jac)
  IF(sdetJac.GT.0.) THEN
   sdetJac=1./sdetJac
  ELSE !shit
   ! Newton has not converged !?!?
   IF(Mode.EQ.1)THEN
    CALL abort(&
__STAMP__&
, 'Newton in FindXiForPartPos singular. iter,sdetJac',NewtonIter,sDetJac)
   ELSE
     Xi(1)=HUGE(1.0)
     Xi(2)=Xi(1)
     Xi(3)=Xi(1)
     RETURN
   END IF
  ENDIF 
  sJac=getInv(Jac,sdetJac)

  ! Iterate Xi using Newton step
  ! Use FAIL
  !Xi = Xi - MATMUL(sJac,F)
  deltaXi=MATMUL(sJac,F)
  Xi = Xi - deltaXI!MATMUL(sJac,F)
  deltaXi2=DOT_PRODUCT(deltaXi,deltaXi)


  IF(ANY(ABS(Xi).GT.1.8)) THEN
    IF(Mode.EQ.1)THEN
      IPWRITE(UNIT_stdOut,*) ' Particle not inside of element, force!!!'
      IPWRITE(UNIT_stdOut,*) ' Newton-Iter', NewtonIter
      IPWRITE(UNIT_stdOut,*) ' xi  ', xi(1:3)
      IPWRITE(UNIT_stdOut,*) ' PartPos', X_in
      IPWRITE(UNIT_stdOut,*) ' ElemID', ElemID+offSetElem
      !IF(PRESENT(PartID)) IPWRITE(UNIT_stdOut,*) ' PartID', PartID
      CALL abort(&
__STAMP__&
,'Particle Not inSide of Element, ElemID,',ElemID)
    ELSE
      EXIT
    END IF
  END IF
  
  ! Compute function value
  CALL LagrangeInterpolationPolys(Xi(1),N_In,XiCL_N_in,wBaryCL_N_in,Lag(1,:))
  CALL LagrangeInterpolationPolys(Xi(2),N_In,XiCL_N_in,wBaryCL_N_in,Lag(2,:))
  CALL LagrangeInterpolationPolys(Xi(3),N_In,XiCL_N_in,wBaryCL_N_in,Lag(3,:))
  ! F(xi) = x(xi) - x_in
  F=-x_in ! xRp
  DO k=0,N_In
    DO j=0,N_In
      buff=Lag(2,j)*Lag(3,k)
      DO i=0,N_In
        buff2=Lag(1,i)*buff
        F=F+XCL_N_in(:,i,j,k)*buff2
      END DO !l=0,N_In
    END DO !i=0,N_In
  END DO !j=0,N_In
END DO !newton
!print*,'newton iter', newtoniter

END SUBROUTINE RefElemNewton


!RECURSIVE SUBROUTINE RefElemBisection(XiOut,XiA,XiB,X_In,wBaryCL_NGeo,XiCL_NGeo,XCL_NGeo,N_In,Found)
!===================================================================================================================================
! bisection algorithm in reference-element to get initial guess
!===================================================================================================================================
! MODULES                                                                                                                          !
!USE MOD_Basis,                   ONLY:LagrangeInterpolationPolys
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!! INPUT VARIABLES 
!INTEGER,INTENT(IN)               :: N_In
!REAL,INTENT(IN)                  :: X_in ! position in physical space 
!REAL,INTENT(IN)                  :: XiCL_NGeo(0:N_In)               ! position of CL points in reference space
!REAL,INTENT(IN)                  :: XCL_NGeo(0:N_In,0:N_in,0:N_In) ! position of CL points in physical space
!REAL,INTENT(IN)                  :: wBaryCL_NGeo(0:N_In) ! derivation of CL points
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!REAL,INTENT(INOUT)               :: XiOut ! position in reference element
!REAL,INTENT(INOUT)               :: XiA,XiB
!LOGICAL,INTENT(OUT)              :: Found
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL                             :: Lag(0:N_In), F(1:3), Lag2(0:N_In), Lag3(0:N_In)
!REAL                             :: buff,buff2,buff3,XiOutA,XiOutB,dummy
!INTEGER                          :: iter,i,j,k
!LOGICAL                          :: DoXi, DoEta,DoZeta
!LOGICAL                          :: FoundA,FoundB
!===================================================================================================================================

!FoundB=.FALSE.
!FoundA=.FALSE.
!Found =.FALSE.
!
!XiOut=0.5*(XiA+XiB)
!
!! compute f(a) 
!CALL LagrangeInterpolationPolys(XiA,N_In,XiCL_NGeo,wBaryCL_NGeo,Lag(:))
!! f(b)
!CALL LagrangeInterpolationPolys(XiB,N_In,XiCL_NGeo,wBaryCL_NGeo,Lag2(:))
!! f(c)
!CALL LagrangeInterpolationPolys(XiOut,N_In,XiCL_NGeo,wBaryCL_NGeo,Lag3(:))
!
!F(1:3)=-X_in
!DO k=0,N_In
!  DO j=0,N_In
!    buff =Lag (j)*Lag (k)
!    buff2=Lag2(j)*Lag2(k)
!    buff3=Lag3(j)*Lag3(k)
!    DO i=0,N_In
!      F(1)=F(1)+XCL_NGeo(i,j,k)*Lag (i)*buff
!      F(2)=F(2)+XCL_NGeo(i,j,k)*Lag2(i)*buff2
!      F(3)=F(3)+XCL_NGeo(i,j,k)*Lag3(i)*buff2
!    END DO !l=0,N_In
!  END DO !i=0,N_In
!END DO !j=0,N_In
!
!print*,'a , b, c',XiA,XiB,XiOut
!print*,'fa,fb,fc',F(1),F(2),F(3)
!
!IF(F(1)*F(2).GT.0.)THEN
!  dummy=XiB
!  IF(F(1)*F(3).LT.0)THEN
!    XiB=XiOut
!    CALL RefElemBisection(XiOutA,XiA,XiB,X_In,wBaryCL_NGeo,XiCL_NGeo,XCL_NGeo(:,:,:),N_In,FoundA)
!  END IF
!  IF(F(2)*F(3).LT.0)THEN
!    XiA=XiOut
!    XiB=dummy
!    CALL RefElemBisection(XiOutB,XiA,XiB,X_In,wBaryCL_NGeo,XiCL_NGeo,XCL_NGeo(:,:,:),N_In,FoundB)
!  END IF
!  print*,'found double',XiOutA,XiOutB
!  IF(FoundA .AND. FoundB)THEN
!    IF((ABS(XiOutA).GT. 1.0) .AND.(ABS(XiOutB).GT.1.0))THEN
!      XiOut=0.5*(XiOutB+XiOutA)
!      Found=.TRUE.
!      RETURN
!    ELSE IF(ABS(XiOutA).GT. 1.0)THEN
!      XiOut=XiOutB
!      RETURN
!    ELSE
!      XiOut=XiOutA
!      RETURN
!    END IF
!  END IF
!END IF
!
!
!buff=F(3)*F(1)
!IF(buff.LT.0.)THEN
!  XiB = XiOut
!ELSE IF (buff.GT.0.)THEN
!  XiA = XiOut
!ELSE
!  Found=.TRUE.
!  RETURN
!END IF
!
!print*,'iterations'
!DO iter=1,4
!  print*,'iter',iter
!  XiOut=0.5*(XiA+XiB)
!  ! compute f(a) 
!  CALL LagrangeInterpolationPolys(XiA,N_In,XiCL_NGeo,wBaryCL_NGeo,Lag(:))
!  ! f(b)
!  CALL LagrangeInterpolationPolys(XiB,N_In,XiCL_NGeo,wBaryCL_NGeo,Lag2(:))
!  ! f(c)
!  CALL LagrangeInterpolationPolys(XiOut,N_In,XiCL_NGeo,wBaryCL_NGeo,Lag3(:))
!  
!  F(1:3)=-X_in
!  DO k=0,N_In
!    DO j=0,N_In
!      buff =Lag (j)*Lag (k)
!      buff2=Lag2(j)*Lag2(k)
!      buff3=Lag3(j)*Lag3(k)
!      DO i=0,N_In
!        F(1)=F(1)+XCL_NGeo(i,j,k)*Lag (i)*buff
!        F(2)=F(2)+XCL_NGeo(i,j,k)*Lag2(i)*buff2
!        F(3)=F(3)+XCL_NGeo(i,j,k)*Lag3(i)*buff2
!      END DO !l=0,N_In
!    END DO !i=0,N_In
!  END DO !j=0,N_In
!  print*,'a , b, c',XiA,XiB,XiOut
!  print*,'fa,fb,fc',F(1),F(2),F(3)
!  read*
!
!  IF(F(1)*F(2).GT.0.)THEN
!    dummy=XiB
!    IF(F(1)*F(3).LT.0)THEN
!      XiB=XiOut
!      CALL RefElemBisection(XiOut,XiA,XiB,X_In,wBaryCL_NGeo,XiCL_NGeo,XCL_NGeo(:,:,:),N_In,FoundA)
!      IF(FoundA) RETURN
!    END IF
!    IF(F(2)*F(3).LT.0)THEN
!      XiA=XiOut
!      XiB=dummy
!      CALL RefElemBisection(XiOut,XiA,XiB,X_In,wBaryCL_NGeo,XiCL_NGeo,XCL_NGeo(:,:,:),N_In,FoundB)
!      IF(FoundB) RETURN
!    END IF
!  END IF
!  buff=F(3)*F(1)
!  IF(buff.LT.0.)THEN
!    XiB = XiOut
!  ELSE IF (buff.GT.0.)THEN
!    XiA = XiOut
!  ELSE
!    Found=.TRUE.
!    RETURN
!  END IF
!END DO ! iter=1,4
!
!Found=.TRUE.
!
!END SUBROUTINE RefElemBisection

FUNCTION getDet(Mat)
!=================================================================================================================================
! compute determinant of 3x3 matrix
!=================================================================================================================================
  ! MODULES
  ! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: Mat(3,3)
!---------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL             :: getDet
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!=================================================================================================================================
getDet=   ( Mat(1,1) * Mat(2,2) - Mat(1,2) * Mat(2,1) ) * Mat(3,3) &
        + ( Mat(1,2) * Mat(2,3) - Mat(1,3) * Mat(2,2) ) * Mat(3,1) &
        + ( Mat(1,3) * Mat(2,1) - Mat(1,1) * Mat(2,3) ) * Mat(3,2)
END FUNCTION getDet


FUNCTION getInv(Mat,sdet)
!=================================================================================================================================
! compute inverse of 3x3 matrix, needs sDet=1/det(Mat)
!=================================================================================================================================
  ! MODULES
  ! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: Mat(3,3),sDet
!---------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL             :: getInv(3,3)
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!=================================================================================================================================
getInv(1,1) = ( Mat(2,2) * Mat(3,3) - Mat(2,3) * Mat(3,2) ) * sdet
getInv(1,2) = ( Mat(1,3) * Mat(3,2) - Mat(1,2) * Mat(3,3) ) * sdet
getInv(1,3) = ( Mat(1,2) * Mat(2,3) - Mat(1,3) * Mat(2,2) ) * sdet
getInv(2,1) = ( Mat(2,3) * Mat(3,1) - Mat(2,1) * Mat(3,3) ) * sdet
getInv(2,2) = ( Mat(1,1) * Mat(3,3) - Mat(1,3) * Mat(3,1) ) * sdet
getInv(2,3) = ( Mat(1,3) * Mat(2,1) - Mat(1,1) * Mat(2,3) ) * sdet
getInv(3,1) = ( Mat(2,1) * Mat(3,2) - Mat(2,2) * Mat(3,1) ) * sdet
getInv(3,2) = ( Mat(1,2) * Mat(3,1) - Mat(1,1) * Mat(3,2) ) * sdet
getInv(3,3) = ( Mat(1,1) * Mat(2,2) - Mat(1,2) * Mat(2,1) ) * sdet
END FUNCTION getInv 


FUNCTION LimitXi(Xi)
!=================================================================================================================================
! lilmit xi to [-1,1]
!=================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: Xi(3)
!---------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL             :: LimitXi(3)
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!=================================================================================================================================

LimitXi=MAX(MIN(1.0d0,XI),-1.0d0)

END FUNCTION LimitXi 
!#endif /*PARTICLES*/

END MODULE MOD_Eval_xyz
