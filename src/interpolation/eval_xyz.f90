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
SUBROUTINE eval_xyz_curved(x_in,NVar,N_in,U_In,U_Out,iElem,PartID)
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
USE MOD_Particle_Mesh_Vars,      ONLY:MappingGuess,epsMapping,ElemRadiusNGeo
USE MOD_Particle_Mesh_Vars,      ONLY:XiEtaZetaBasis,ElemBaryNGeo,slenXiEtaZetaBasis
USE MOD_PICInterpolation_Vars,   ONLY:NBG,BGField,useBGField,BGDataSize,BGField_wBary, BGField_xGP,BGType
USE MOD_Mesh_Vars,               ONLY:offsetElem
!USE MOD_Particle_Vars,           ONLY:PartPosRef
!USE MOD_Mesh_Vars,ONLY: X_CP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: NVar                                  ! 6 (Ex, Ey, Ez, Bx, By, Bz) 
INTEGER,INTENT(IN)  :: N_In                                  ! usually PP_N
INTEGER,INTENT(IN)  :: iElem                                 ! elem index
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
REAL                :: deltaXi(1:3),deltaXi2
INTEGER             :: NewTonIter
!REAL                :: X3D_Buf1(1:NVar,0:N_In,0:N_In)  ! first intermediate results from 1D interpolations
!REAL                :: X3D_Buf2(1:NVar,0:N_In) ! second intermediate results from 1D interpolations
REAL                :: Winner_Dist,Dist,abortcrit
REAL, PARAMETER     :: EPSONE=1.00000001
INTEGER             :: iDir
REAL                :: F(1:3),Lag(1:3,0:NGeo)
REAL                :: L_xi(3,0:PP_N), L_eta_zeta
REAL                :: Jac(1:3,1:3),sdetJac,sJac(1:3,1:3)
REAL                :: buff,buff2
REAL                :: Ptild(1:3),XiLinear(1:6)
! h5-external e,b field
REAL,ALLOCATABLE    :: L_xi_BGField(:,:), U_BGField(:)
!===================================================================================================================================

! get initial guess by nearest GP search ! simple guess
SELECT CASE(MappingGuess)
CASE(1)
  Ptild=X_in - ElemBaryNGeo(:,iElem)
  ! plus coord system (1-3) and minus coord system (4-6)
  DO iDir=1,6
    XiLinear(iDir)=DOT_PRODUCT(Ptild,XiEtaZetaBasis(:,iDir,iElem))*slenXiEtaZetaBasis(iDir,iElem)
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
    Dist=SUM((x_in(:)-Elem_xGP(:,i,j,k,iElem))*(x_in(:)-Elem_xGP(:,i,j,k,iElem)))
    IF (Dist.LT.Winner_Dist) THEN
      Winner_Dist=Dist
      Xi(:)=(/xGP(i),xGP(j),xGP(k)/) ! start value
    END IF
  END DO; END DO; END DO
CASE(3) 
  ! compute distance on XCL Points
  Winner_Dist=HUGE(1.)
  DO i=0,NGeo; DO j=0,NGeo; DO k=0,NGeo
    Dist=SUM((x_in(:)-XCL_NGeo(:,i,j,k,iElem))*(x_in(:)-XCL_NGeo(:,i,j,k,iElem)))
    IF (Dist.LT.Winner_Dist) THEN
      Winner_Dist=Dist
      Xi(:)=(/XiCL_NGeo(i),XiCL_NGeo(j),XiCL_NGeo(k)/) ! start value
    END IF
  END DO; END DO; END DO
CASE(4)
  ! trival guess 
  xi=0.
END SELECT

! initial guess
CALL LagrangeInterpolationPolys(Xi(1),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
CALL LagrangeInterpolationPolys(Xi(2),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
CALL LagrangeInterpolationPolys(Xi(3),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))
! F(xi) = x(xi) - x_in
F=-x_in ! xRp
DO k=0,NGeo
  DO j=0,NGeo
    DO i=0,NGeo
      F=F+XCL_NGeo(:,i,j,k,iElem)*Lag(1,i)*Lag(2,j)*Lag(3,k)
      !print*,'XCL',i,j,k,XCL_NGeo(:,i,j,k,iElem)
    END DO !l=0,NGeo
  END DO !i=0,NGeo
END DO !j=0,NGeo

NewtonIter=0
!abortCrit=ElemRadiusNGeo(iElem)*ElemRadiusNGeo(iElem)*epsMapping
deltaXi2=HUGE(1.0)
!DO WHILE ((SUM(F*F).GT.abortCrit).AND.(NewtonIter.LT.50))
DO WHILE((deltaXi2.GT.epsMapping).AND.(NewtonIter.LT.100))
!DO WHILE ((DOT_PRODUCT(F,F).GT.abortCrit).AND.(NewtonIter.LT.100))
  NewtonIter=NewtonIter+1

  ! caution, dXCL_NGeo is transposed of required matrix
  Jac=0.
  DO k=0,NGeo
    DO j=0,NGeo
      buff=Lag(2,j)*Lag(3,k)
      DO i=0,NGeo
        buff2=Lag(1,i)*buff
        Jac(1,1:3)=Jac(1,1:3)+dXCL_NGeo(1:3,1,i,j,k,iElem)*buff2
        Jac(2,1:3)=Jac(2,1:3)+dXCL_NGeo(1:3,2,i,j,k,iElem)*buff2
        Jac(3,1:3)=Jac(3,1:3)+dXCL_NGeo(1:3,3,i,j,k,iElem)*buff2
      END DO !i=0,NGeo
    END DO !j=0,NGeo
  END DO !k=0,NGeo
  
  ! Compute inverse of Jacobian
  sdetJac=getDet(Jac)
  IF(sdetJac.GT.0.) THEN
   sdetJac=1./sdetJac
  ELSE !shit
   ! Newton has not converged !?!?
   CALL abort(__STAMP__, &
        'Newton in FindXiForPartPos singular. iter,sdetJac',NewtonIter,sDetJac)
  ENDIF 
  sJac=getInv(Jac,sdetJac)

  ! Iterate Xi using Newton step
  ! Use FAIL
  !Xi = Xi - MATMUL(sJac,F)
  deltaXi=MATMUL(sJac,F)
  Xi = Xi - deltaXI!MATMUL(sJac,F)
  deltaXi2=DOT_PRODUCT(deltaXi,deltaXi)


  IF(ANY(ABS(Xi).GT.1.8)) THEN
  !IF((NewtonIter.GE.4).AND.(ANY(ABS(Xi).GT.1.5)))THEN
    IPWRITE(UNIT_stdOut,*) ' Particle not inside of element, force!!!'
    IPWRITE(UNIT_stdOut,*) ' Newton-Iter', NewtonIter
    IPWRITE(UNIT_stdOut,*) ' xi  ', xi(1:3)
    IPWRITE(UNIT_stdOut,*) ' PartPos', X_in
    IPWRITE(UNIT_stdOut,*) ' ElemID', iElem+offSetElem
    IF(PRESENT(PartID)) IPWRITE(UNIT_stdOut,*) ' PartID', PartID
    CALL abort(__STAMP__, &
        'Particle Not inSide of Element, iElem, iPart',iElem,REAL(PartID))
  END IF
  
  ! Compute function value
  CALL LagrangeInterpolationPolys(Xi(1),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
  CALL LagrangeInterpolationPolys(Xi(2),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
  CALL LagrangeInterpolationPolys(Xi(3),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))
  ! F(xi) = x(xi) - x_in
  F=-x_in ! xRp
  DO k=0,NGeo
    DO j=0,NGeo
      buff=Lag(2,j)*Lag(3,k)
      DO i=0,NGeo
        buff2=Lag(1,i)*buff
        F=F+XCL_NGeo(:,i,j,k,iElem)*buff2
      END DO !l=0,NGeo
    END DO !i=0,NGeo
  END DO !j=0,NGeo
END DO !newton

! IF(ANY(ABS(Xi).GT.1.1)) THEN
! !IF((NewtonIter.GE.4).AND.(ANY(ABS(Xi).GT.1.5)))THEN
!   IPWRITE(UNIT_stdOut,*) ' Particle not inside of element, force!!!'
!   IPWRITE(UNIT_stdOut,*) ' Newton-Iter', NewtonIter
!   IPWRITE(UNIT_stdOut,*) ' xi  ', xi(1)
!   IPWRITE(UNIT_stdOut,*) ' eta ', xi(2)
!   IPWRITE(UNIT_stdOut,*) ' zeta', xi(3)
!   IPWRITE(UNIT_stdOut,*) ' PartPos', X_in
!   IPWRITE(UNIT_stdOut,*) ' ElemID', iElem+offSetElem
!   IF(PRESENT(PartID)) IPWRITE(UNIT_stdOut,*) ' PartID', PartID
!   CALL abort(__STAMP__, &
!       'Particle Not inSide of Element, iElem, iPart',iElem,REAL(PartID))
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
  CALL LagrangeInterpolationPolys(xi(1),NBG,BGField_xGP,BGField_wBary,L_xi_BGField(1,:))
  CALL LagrangeInterpolationPolys(xi(2),NBG,BGField_xGP,BGField_wBary,L_xi_BGField(2,:))
  CALL LagrangeInterpolationPolys(xi(3),NBG,BGField_xGP,BGField_wBary,L_xi_BGField(3,:))
  
  U_BGField(:)=0
  DO k=0,N_in
    DO j=0,N_in
      L_eta_zeta=L_xi_BGField(2,j)*L_xi_BGField(3,k)
      DO i=0,N_in
        U_BGField = U_BGField + BGField(:,i,j,k,iElem)*L_xi_BGField(1,i)*L_Eta_Zeta
      END DO ! i=0,N_In
    END DO ! j=0,N_In
  END DO ! k=0,N_In


  !! 2.2) do the tensor product thing 
  !X3D_tmp1=0.
  !! first direction iN_In
  !DO k=0,NBG
  !  DO j=0,NBG
  !    DO i=0,NBG
  !      X3D_tmp1(:,j,k)=X3D_tmp1(:,j,k)+L_xi_BGField(1,i)*BGField(:,i,j,k,iElem)
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
  DEALLOCATE( L_xi_BGField, U_BGField)! X3d_tmp1, x3d_tmp2, x3d_tmp3)
END IF ! useBGField

END SUBROUTINE eval_xyz_curved


SUBROUTINE eval_xyz_elemcheck(x_in,xi,iElem,PartID,DoReUseMap)
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
USE MOD_Particle_Mesh_Vars,      ONLY:MappingGuess,epsMapping
USE MOD_Particle_Mesh_Vars,      ONLY:XiEtaZetaBasis,ElemBaryNGeo,slenXiEtaZetaBasis,ElemRadiusNGeo
USE MOD_Mesh_Vars,               ONLY:dXCL_NGeo,Elem_xGP,XCL_NGeo,NGeo,wBaryCL_NGeo,XiCL_NGeo,NGeo
USE MOD_Mesh_Vars,               ONLY:offsetElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)          :: iElem                                 ! elem index
REAL,INTENT(IN)             :: x_in(3)                                  ! physical position of particle 
INTEGER,INTENT(IN),OPTIONAL :: PartID
LOGICAL,INTENT(IN),OPTIONAL :: DoReUseMap
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)          :: xi(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                    :: i,j,k
REAL                       :: deltaXi(1:3),deltaXi2
REAL                       :: epsOne
INTEGER                    :: NewTonIter
REAL                       :: Winner_Dist,Dist
INTEGER                    :: idir
REAL                       :: F(1:3),Lag(1:3,0:NGeo)
REAL                       :: Jac(1:3,1:3),sdetJac,sJac(1:3,1:3)
REAL                       :: buff,buff2,abortcrit
REAL                       :: Ptild(1:3),XiLinear(1:6)
!===================================================================================================================================

epsOne=1.0+epsMapping
IF(.NOT.PRESENT(DoReUseMap))THEN
  SELECT CASE(MappingGuess)
  CASE(1)
    Ptild=X_in - ElemBaryNGeo(:,iElem)
    ! plus coord system (1-3) and minus coord system (4-6)
    DO iDir=1,6
      XiLinear(iDir)=DOT_PRODUCT(Ptild,XiEtaZetaBasis(:,iDir,iElem))*slenXiEtaZetaBasis(iDir,iElem)
    END DO
    ! compute guess as average value
    DO iDir=1,3
      Xi(iDir)=0.5*(XiLinear(iDir)-XiLinear(iDir+3))
    END DO 
    IF(MAXVAL(ABS(Xi)).GT.epsOne) Xi=LimitXi(Xi)
  CASE(2)
    Winner_Dist=HUGE(1.)
    DO i=0,PP_N; DO j=0,PP_N; DO k=0,PP_N
      Dist=SUM((x_in(:)-Elem_xGP(:,i,j,k,iElem))*(x_in(:)-Elem_xGP(:,i,j,k,iElem)))
      IF (Dist.LT.Winner_Dist) THEN
        Winner_Dist=Dist
        Xi(:)=(/xGP(i),xGP(j),xGP(k)/) ! start value
      END IF
    END DO; END DO; END DO
  CASE(3) 
    ! compute distance on XCL Points
    Winner_Dist=HUGE(1.)
    DO i=0,NGeo; DO j=0,NGeo; DO k=0,NGeo
      Dist=SUM((x_in(:)-XCL_NGeo(:,i,j,k,iElem))*(x_in(:)-XCL_NGeo(:,i,j,k,iElem)))
      IF (Dist.LT.Winner_Dist) THEN
        Winner_Dist=Dist
        Xi(:)=(/XiCL_NGeo(i),XiCL_NGeo(j),XiCL_NGeo(k)/) ! start value
      END IF
    END DO; END DO; END DO
  CASE(4)
    ! trival guess, cell mean point
    xi=0.
  END SELECT
END IF

! initial guess
CALL LagrangeInterpolationPolys(Xi(1),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
CALL LagrangeInterpolationPolys(Xi(2),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
CALL LagrangeInterpolationPolys(Xi(3),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))
! F(xi) = x(xi) - x_in
F=-x_in ! xRp
DO k=0,NGeo
  DO j=0,NGeo
    buff=Lag(2,j)*Lag(3,k)
    DO i=0,NGeo
      buff2=Lag(1,i)*buff
      F=F+XCL_NGeo(:,i,j,k,iElem)*buff2
      !F=F+XCL_NGeo(:,i,j,k,iElem)*Lag(1,i)*Lag(2,j)*Lag(3,k)
    END DO !l=0,NGeo
  END DO !i=0,NGeo
END DO !j=0,NGeo

NewtonIter=0
!abortCrit=ElemRadiusNGeo(iElem)*ElemRadiusNGeo(iElem)*epsMapping
!abortCrit=ElemRadiusNGeo(iElem)*ElemRadiusNGeo(iElem)*epsMapping
deltaXi2=HUGE(1.0)
!abortCrit=DOT_PRODUCT(F,F)*epsMapping
!DO WHILE ((SUM(F*F).GT.abortCrit).AND.(NewtonIter.LT.50))
!@DO WHILE ((DOT_PRODUCT(F,F).GT.abortCrit).AND.(NewtonIter.LT.100))
DO WHILE((deltaXi2.GT.epsMapping).AND.(NewtonIter.LT.100))
  NewtonIter=NewtonIter+1
  ! 
  ! caution, dXCL_NGeo is transposed of required matrix
  Jac=0.
  DO k=0,NGeo
    DO j=0,NGeo
      buff=Lag(2,j)*Lag(3,k)
      DO i=0,NGeo
        buff2=Lag(1,i)*buff
        Jac(1,1:3)=Jac(1,1:3)+dXCL_NGeo(1:3,1,i,j,k,iElem)*buff2
        Jac(2,1:3)=Jac(2,1:3)+dXCL_NGeo(1:3,2,i,j,k,iElem)*buff2
        Jac(3,1:3)=Jac(3,1:3)+dXCL_NGeo(1:3,3,i,j,k,iElem)*buff2
      END DO !i=0,NGeo
    END DO !j=0,NGeo
  END DO !k=0,NGeo
  

  ! Compute inverse of Jacobian
  sdetJac=getDet(Jac)
  IF(sdetJac.GT.0.) THEN
   sdetJac=1./sdetJac
  ELSE !shit
   ! Newton has not converged !?!?
   Xi(1)=HUGE(1.0)
   Xi(2)=Xi(1)
   Xi(3)=Xi(1)
   EXIT
   !IPWRITE(*,*) ' GlobalElemID ', OffSetElem+iElem
   !IPWRITE(*,*) ' PartPos ', X_in
   !CALL abort(__STAMP__, &
   !     'Newton in FindXiForPartPos singular. iter,sdetJac',NewtonIter,sDetJac)
  ENDIF 
  sJac=getInv(Jac,sdetJac)
  
  ! Iterate Xi using Newton step
  ! Use FAIL
  deltaXi=MATMUL(sJac,F)
  Xi = Xi - deltaXI!MATMUL(sJac,F)
  deltaXi2=DOT_PRODUCT(deltaXi,deltaXi)

  !IF((NewtonIter.GE.4).AND.(ANY(ABS(Xi).GT.1.5)))THEN
  IF(ANY(ABS(Xi).GT.1.8))THEN
!    IF(PRESENT(PartID)) THEN
!      IPWRITE(UNIT_stdOut,*) 'ParticleID', PartID
!      IPWRITE(UNIT_stdOut,*) ' Particle not inside of element!!!'
!      IPWRITE(UNIT_stdOut,*) ' Element', iElem
!      IPWRITE(UNIT_stdOut,*) ' xi  ', xi(1)
!      IPWRITE(UNIT_stdOut,*) ' eta ', xi(2)
!      IPWRITE(UNIT_stdOut,*) ' zeta', xi(3)
!      IPWRITE(UNIT_stdOut,*) ' PartPos', X_in
!      !CALL abort(__STAMP__, &
!      !    'Particle Not inSide of Element')
!    ELSE
      EXIT
!    END IF
  END IF
  
  ! Compute function value
  CALL LagrangeInterpolationPolys(Xi(1),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
  CALL LagrangeInterpolationPolys(Xi(2),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
  CALL LagrangeInterpolationPolys(Xi(3),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))
  ! F(xi) = x(xi) - x_in
  F=-x_in ! xRp
  DO k=0,NGeo
    DO j=0,NGeo
      buff=Lag(2,j)*Lag(3,k)
      DO i=0,NGeo
        buff2=Lag(1,i)*buff
        F=F+XCL_NGeo(:,i,j,k,iElem)*buff2
        !F=F+XCL_NGeo(:,i,j,k,iElem)*Lag(1,i)*Lag(2,j)*Lag(3,k)
      END DO !l=0,NGeo
    END DO !i=0,NGeo
  END DO !j=0,NGeo
END DO !newton


!IF(MAXVAL(ABS(Xi)).GT.1.0) THEN
! ! IF(PRESENT(PartID).AND.PartID.EQ.813) THEN
!  IPWRITE(UNIT_stdOut,*) ' Particle not inside of element!!!'
!  IF(PRESENT(PartID)) IPWRITE(UNIT_stdOut,*) 'ParticleID', PartID
!  IPWRITE(UNIT_stdOut,*) ' Element', iElem
!  IPWRITE(UNIT_stdOut,*) ' xi  ', xi(:)
!  !  IPWRITE(UNIT_stdOut,*) ' PartPos', X_in
!  !END IF
!END IF
!

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


SUBROUTINE eval_xyz_part2(xi_in,NVar,N_in,U_In,U_Out,iElem)
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
INTEGER,INTENT(IN)        :: iElem                                 ! elem index
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
REAL                :: F(1:3)!,Lag(1:3,0:NGeo)
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
  DO k=0,N_in
    DO j=0,N_in
      L_eta_zeta=L_xi_BGField(2,j)*L_xi_BGField(3,k)
      DO i=0,N_in
        U_BGField = U_BGField + BGField(:,i,j,k,iElem)*L_xi_BGField(1,i)*L_Eta_Zeta
      END DO ! i=0,N_In
    END DO ! j=0,N_In
  END DO ! k=0,N_In


  !! 2.2) do the tensor product thing 
  !X3D_tmp1=0.
  !! first direction iN_In
  !DO k=0,NBG
  !  DO j=0,NBG
  !    DO i=0,NBG
  !      X3D_tmp1(:,j,k)=X3D_tmp1(:,j,k)+L_xi_BGField(1,i)*BGField(:,i,j,k,iElem)
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
