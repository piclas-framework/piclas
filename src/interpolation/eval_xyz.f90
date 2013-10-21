MODULE MOD_Eval_xyz
!===================================================================================================================================
! Changes a 3D Tensor Product Lagrange Points of Lagrange Basis of degree N_In to  
! Lagrange points of a Lagrange Basis for one point, using two
! arbitrary point disributions xi_In(0:N_In) and xi_Out 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
  INTEGER :: errorflag
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE eval_xyz
  MODULE PROCEDURE eval_xyz 
END INTERFACE

INTERFACE eval_xyz_fast
  MODULE PROCEDURE eval_xyz_fast
END INTERFACE

PUBLIC :: eval_xyz, eval_xyz_fast, eval_xyz_part2, Calc_F, Calc_dF_inv 
!===================================================================================================================================

CONTAINS

SUBROUTINE eval_xyz_fast(x_in,NVar,N_in,X3D_In,X3D_Out,iElem)
!===================================================================================================================================
! interpolate a 3D tensor product Lagrange basis defined by (N_in+1) 1D interpolation point positions x
! first get xi,eta,zeta from x,y,z...then do tenso product interpolation
! xi is defined in the 1DrefElem xi=[-1,1]
!===================================================================================================================================
! MODULES
USE MOD_Basis,ONLY: LagrangeInterpolationPolys
USE MOD_Interpolation_Vars,ONLY: wBary,xGP
USE MOD_Particle_Vars,ONLY:GEO
USE MOD_PICInterpolation_Vars, ONLY: ElemT_inv
!USE MOD_Mesh_Vars,ONLY: X_CP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: NVar                                  ! 6 (Ex, Ey, Ez, Bx, By, Bz) 
INTEGER,INTENT(IN)  :: N_In                                  ! usually PP_N
INTEGER,INTENT(IN)  :: iElem                                 ! elem index
REAL,INTENT(IN)     :: X3D_In(1:NVar,0:N_In,0:N_In,0:N_In)   ! elem state
REAL,INTENT(IN)     :: x_in(3)                                  ! physical position of particle 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: X3D_Out(1:NVar)  ! Interpolated state
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER             :: iN_In,jN_In,kN_In,iN_Out,jN_Out,kN_Out,i,j,k, iNode
REAL                :: X3D_Buf1(1:NVar,0:N_In,0:N_In)  ! first intermediate results from 1D interpolations
REAL                :: X3D_Buf2(1:NVar,0:N_In) ! second intermediate results from 1D interpolations
REAL                :: xi(3)
REAL,DIMENSION(3,0:N_in)  :: L_xi        
REAL              :: P(3,8), F(3), dF_inv(3,3), s(3) 
REAL, PARAMETER   :: EPS=1E-8
REAL              :: PT(3,8), T_inv(3,3), DP(3)
INTEGER           :: n_Newton, iPoint
!===================================================================================================================================
errorflag = 0
! --------------------------------------------------
! 1.) Mapping: get xi,eta,zeta value from x,y,z
! --------------------------------------------------
! 1.1.) initial guess from linear part:
DO iNode = 1,8
  P(1:3,iNode) = GEO%NodeCoords(1:3,GEO%ElemToNodeID(iNode,iElem))
END DO
!PT(1:3,1:8)    = ElemPT(1:3,1:8,iElem) 
T_inv(1:3,1:3) = ElemT_inv(1:3,1:3,iElem)


! transform also the physical coordinate of the point 
! into the unit element (this is the solution of the 
! linear problem already 
xi = 0.
DP = x_in - P(:,1)
DO i=1,3
  DO j=1,3
    xi(i)= xi(i) + T_inv(i,j) * DP(j) 
  END DO
END DO
xi = xi - (/1.,1.,1./)

! 1.2.) Newton-Method to solve non-linear part
!       If linear elements then F should becom 0 and no 
!       Newton step is required.

F = Calc_F(xi,x_in,P)
n_Newton = 0.
DO WHILE(SUM(ABS(F)).GE.EPS) 
  dF_inv = Calc_dF_inv(xi,P)
  s=0.
  DO j = 1,3
    DO k = 1,3
      s(j) = s(j) + dF_inv(j,k) * F(k) 
    END DO ! k
  END DO ! j
  xi = xi - s
  F = Calc_F(xi,x_in,P)
  n_Newton = n_Newton + 1
END DO ! i
!print*, "n-Steps", n_Newton
!print*,"F", F
!print*, "xi", xi
!read*

! --------------------------------------------------
! 2.) Interpolation
! --------------------------------------------------
! 2.1) get "Vandermonde" vectors
DO i=1,3
  CALL LagrangeInterpolationPolys(xi(i),N_in,xGP,wBary,L_xi(i,:))
END DO

! 2.2) do the tensor product thing 
X3D_buf1=0.
! first direction iN_In
DO kN_In=0,N_In
  DO jN_In=0,N_In
    DO iN_In=0,N_In
      X3D_Buf1(:,jN_In,kN_In)=X3D_Buf1(:,jN_In,kN_In)+L_xi(1,iN_In)*X3D_In(:,iN_In,jN_In,kN_In)
    END DO
  END DO
END DO
X3D_buf2=0.
! second direction jN_In
DO kN_In=0,N_In
  DO jN_In=0,N_In
    X3D_Buf2(:,kN_In)=X3D_Buf2(:,kN_In)+L_xi(2,jN_In)*X3D_Buf1(:,jN_In,kN_In)
  END DO
END DO
X3D_Out=0.
! last direction kN_In
DO kN_In=0,N_In
  X3D_Out(:)=X3D_Out(:)+L_xi(3,kN_In)*X3D_Buf2(:,kN_In)
END DO
END SUBROUTINE eval_xyz_fast

SUBROUTINE eval_xyz_part2(xi,NVar,N_in,X3D_In,X3D_Out,iElem)
!===================================================================================================================================
! Same as eval_xyz_fast, with the position already mapped to -1|1 space (xi instead of x_in)
!===================================================================================================================================
! MODULES
USE MOD_Basis,ONLY: LagrangeInterpolationPolys
USE MOD_Interpolation_Vars,ONLY: wBary,xGP
USE MOD_Particle_Vars,ONLY:GEO
USE MOD_PICInterpolation_Vars, ONLY: ElemT_inv
!USE MOD_Mesh_Vars,ONLY: X_CP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: NVar                                  ! 6 (Ex, Ey, Ez, Bx, By, Bz) 
INTEGER,INTENT(IN)  :: N_In                                  ! usually PP_N
INTEGER,INTENT(IN)  :: iElem                                 ! elem index
REAL,INTENT(IN)     :: X3D_In(1:NVar,0:N_In,0:N_In,0:N_In)   ! elem state
REAL,INTENT(IN)     :: xi(3)                                  ! physical position of particle 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: X3D_Out(1:NVar)  ! Interpolated state
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER             :: iN_In,jN_In,kN_In,iN_Out,jN_Out,kN_Out,i,j,k, iNode
REAL                :: X3D_Buf1(1:NVar,0:N_In,0:N_In)  ! first intermediate results from 1D interpolations
REAL                :: X3D_Buf2(1:NVar,0:N_In) ! second intermediate results from 1D interpolations
REAL,DIMENSION(3,0:N_in)  :: L_xi        
REAL              :: P(3,8), F(3), dF_inv(3,3), s(3) 
REAL, PARAMETER   :: EPS=1E-8
REAL              :: PT(3,8), T_inv(3,3), DP(3)
INTEGER           :: n_Newton, iPoint
!===================================================================================================================================
errorflag = 0

! --------------------------------------------------
! 2.) Interpolation
! --------------------------------------------------
! 2.1) get "Vandermonde" vectors
DO i=1,3
  CALL LagrangeInterpolationPolys(xi(i),N_in,xGP,wBary,L_xi(i,:))
END DO

! 2.2) do the tensor product thing 
X3D_buf1=0.
! first direction iN_In
DO kN_In=0,N_In
  DO jN_In=0,N_In
    DO iN_In=0,N_In
      X3D_Buf1(:,jN_In,kN_In)=X3D_Buf1(:,jN_In,kN_In)+L_xi(1,iN_In)*X3D_In(:,iN_In,jN_In,kN_In)
    END DO
  END DO
END DO
X3D_buf2=0.
! second direction jN_In
DO kN_In=0,N_In
  DO jN_In=0,N_In
    X3D_Buf2(:,kN_In)=X3D_Buf2(:,kN_In)+L_xi(2,jN_In)*X3D_Buf1(:,jN_In,kN_In)
  END DO
END DO
X3D_Out=0.
! last direction kN_In
DO kN_In=0,N_In
  X3D_Out(:)=X3D_Out(:)+L_xi(3,kN_In)*X3D_Buf2(:,kN_In)
END DO
END SUBROUTINE eval_xyz_part2


SUBROUTINE eval_xyz(x_in,NVar,N_in,X3D_In,X3D_Out,iElem)
!===================================================================================================================================
! interpolate a 3D tensor product Lagrange basis defined by (N_in+1) 1D interpolation point positions x
! first get xi,eta,zeta from x,y,z...then do tenso product interpolation
! xi is defined in the 1DrefElem xi=[-1,1]
!===================================================================================================================================
! MODULES
USE MOD_Basis,ONLY: LagrangeInterpolationPolys
USE MOD_Interpolation_Vars,ONLY: wBary,xGP
USE MOD_Particle_Vars,ONLY:GEO
!USE MOD_Mesh_Vars,ONLY: X_CP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: NVar                                  ! 6 (Ex, Ey, Ez, Bx, By, Bz) 
INTEGER,INTENT(IN)  :: N_In                                  ! usually PP_N
INTEGER,INTENT(IN)  :: iElem                                 ! elem index
REAL,INTENT(IN)     :: X3D_In(1:NVar,0:N_In,0:N_In,0:N_In)   ! elem state
REAL,INTENT(IN)     :: x_in(3)                                  ! physical position of particle 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: X3D_Out(1:NVar)  ! Interpolated state
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER             :: iN_In,jN_In,kN_In,iN_Out,jN_Out,kN_Out,i,j,k, iNode
REAL                :: X3D_Buf1(1:NVar,0:N_In,0:N_In)  ! first intermediate results from 1D interpolations
REAL                :: X3D_Buf2(1:NVar,0:N_In) ! second intermediate results from 1D interpolations
REAL                :: xi(3)
REAL,DIMENSION(3,0:N_in)  :: L_xi        
REAL              :: detjb
REAL              :: K1(3)
REAL              :: KM_inv(3,3)
REAL              :: P(3,8), F(3), dF_inv(3,3), s(3) 
REAL, PARAMETER   :: EPS=1E-8
REAL              :: m, x(3)
REAL              :: PT(3,8), T(3,3), T_inv(3,3), DP(3), xT(3)
INTEGER           :: n_Newton
!===================================================================================================================================
! --------------------------------------------------
! 1.) Mapping: get xi,eta,zeta value from x,y,z
! --------------------------------------------------
! 1.1.) initial guess from linear part:
!P = X_CP(:,:,iElem)
!print*,"Pos:",x_in
DO iNode = 1,8
  P(1:3,iNode) = GEO%NodeCoords(1:3,GEO%ElemToNodeID(iNode,iElem))
  !print*, "iNode",iNode, P(:,iNode) 
END DO
!WRITE(*,*)'Not implemented yet!!!'
!STOP
!m = 1000
!x = (/0.3,0.4,0.5/)
!x = x * m
x = x_in
!P(:,1) = (/-1.,-1.,-1./)
!P(:,2) = (/ 1.,-1.,-1./)
!P(:,3) = (/ 1., 1.,-1./)
!P(:,4) = (/-1., 1.,-1./)
!P(:,5) = (/-1.,-1., 1./)
!P(:,6) = (/1.,-1., 1./)
!P(:,7) = (/1.,1.2, 1./)
!P(:,8) = (/-1.2, 1., 1./)
!P = P*m

! 1.2.) trafo into unit element due to better matrix condition. 
!       If physical element has dimension of several magnitudes higher
!       than the unit element the abort criteria for the Newton algorithm
!       has big problems to reach the desired precission. Additionally the
!       matrices becomes bad conditioned which causes problem in the inverting
!       algorithm. 
T(:,1) = 0.5 * (P(:,2)-P(:,1))
T(:,2) = 0.5 * (P(:,4)-P(:,1))
T(:,3) = 0.5 * (P(:,5)-P(:,1))
T_inv = Calc_inv(T)
PT = 0.
PT(:,1) = (/-1.,-1.,-1./)
PT(:,2) = (/ 1.,-1.,-1./)
PT(:,4) = (/-1., 1.,-1./)
PT(:,5) = (/-1.,-1., 1./)

! build points in unit elem
DP = P(:,3) - P(:,1)
DO i=1,3
  DO j=1,3
    PT(i,3)= PT(i,3) + T_inv(i,j) * DP(j) 
  END DO
END DO
PT(:,3) = PT(:,3) - (/1.,1.,1./)

DP = P(:,6) - P(:,1)
DO i=1,3
  DO j=1,3
    PT(i,6)= PT(i,6) + T_inv(i,j) * DP(j) 
  END DO
END DO
PT(:,6) = PT(:,6) - (/1.,1.,1./)

DP = P(:,7) - P(:,1)
DO i=1,3
  DO j=1,3
    PT(i,7)= PT(i,7) + T_inv(i,j) * DP(j) 
  END DO
END DO
PT(:,7) = PT(:,7) - (/1.,1.,1./)

DP = P(:,8) - P(:,1)
DO i=1,3
  DO j=1,3
    PT(i,8)= PT(i,8) + T_inv(i,j) * DP(j) 
  END DO
END DO
PT(:,8) = PT(:,8) - (/1.,1.,1./)

! transform also the physical coordinate of the point 
! into the unit element 
DP = x - P(:,1)
DO i=1,3
  DO j=1,3
    xT(i)= xT(i) + T_inv(i,j) * DP(j) 
  END DO
END DO
xT = xT - (/1.,1.,1./)

! build linear guess
K1 = 0.
DO i = 1,8 
  K1(:) = K1(:) + PT(:,i)
END DO
K1 = 0.125 * K1

KM_inv = Calc_KM_inv(PT)

xi = 0.
DO i=1,3
  DO j=1,3
    xi(i) = xi(i) + KM_inv(i,j) * xT(j)  
  END DO
  xi(i) = xi(i) - K1(i) 
END DO
!print*, "xi",xi
!read*

! 1.2.) Newton-Method to solve non-linear part
!       If linear elements then F should becom 0 and no 
!       Newton step is required.

F = Calc_F(xi,xT,PT)
!n_Newton = 0.
DO WHILE(SUM(ABS(F)).GE.EPS) 
  dF_inv = Calc_dF_inv(xi,PT)
  s=0.
  DO j = 1,3
    DO k = 1,3
      s(j) = s(j) + dF_inv(j,k) * F(k) 
    END DO ! k
  END DO ! j
  xi = xi - s
  F = Calc_F(xi,xT,PT)
  !n_Newton = n_Newton + 1
END DO ! i
!print*, "n-Steps", n_Newton
!print*,"F", F
!print*, "xi", xi
!read*

! --------------------------------------------------
! 2.) Interpolation
! --------------------------------------------------
! 2.1) get "Vandermonde" vectors
DO i=1,3
  CALL LagrangeInterpolationPolys(xi(i),N_in,xGP,wBary,L_xi(i,:))
END DO

! 2.2) do the tensor product thing 
X3D_buf1=0.
! first direction iN_In
DO kN_In=0,N_In
  DO jN_In=0,N_In
    DO iN_In=0,N_In
      X3D_Buf1(:,jN_In,kN_In)=X3D_Buf1(:,jN_In,kN_In)+L_xi(1,iN_In)*X3D_In(:,iN_In,jN_In,kN_In)
    END DO
  END DO
END DO
X3D_buf2=0.
! second direction jN_In
DO kN_In=0,N_In
  DO jN_In=0,N_In
    X3D_Buf2(:,kN_In)=X3D_Buf2(:,kN_In)+L_xi(2,jN_In)*X3D_Buf1(:,jN_In,kN_In)
  END DO
END DO
X3D_Out=0.
! last direction kN_In
DO kN_In=0,N_In
  X3D_Out(:)=X3D_Out(:)+L_xi(3,kN_In)*X3D_Buf2(:,kN_In)
END DO
END SUBROUTINE eval_xyz


FUNCTION Calc_F(xi,x,P)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)          :: xi(3)      ! 
REAL,INTENT(IN)          :: x(3)       ! 
REAL,INTENT(IN)          :: P(3,8)     ! 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                     :: Calc_F(3)  !  
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
Calc_F =0.125 *(P(:,1)*(1-xi(1)) * (1-xi(2)) * (1-xi(3))  &
              + P(:,2)*(1+xi(1)) * (1-xi(2)) * (1-xi(3))  &
              + P(:,3)*(1+xi(1)) * (1+xi(2)) * (1-xi(3))  &
              + P(:,4)*(1-xi(1)) * (1+xi(2)) * (1-xi(3))  &
              + P(:,5)*(1-xi(1)) * (1-xi(2)) * (1+xi(3))  &
              + P(:,6)*(1+xi(1)) * (1-xi(2)) * (1+xi(3))  &
              + P(:,7)*(1+xi(1)) * (1+xi(2)) * (1+xi(3))  &
              + P(:,8)*(1-xi(1)) * (1+xi(2)) * (1+xi(3))) &
              - x;
END FUNCTION Calc_F 


FUNCTION Calc_dF_inv(xi,P)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)          :: xi(3)      ! 
REAL,INTENT(IN)          :: P(3,8)     ! 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                     :: Calc_dF_inv(3,3)  !  
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                     :: dF(3,3)  !  
REAL                     :: dF_inv(3,3)  !  
REAL                     :: detjb
INTEGER                  :: i
!===================================================================================================================================
dF = 0.
dF(:,1)=0.125 * ((P(:,2)-P(:,1))*(1-xi(2))*(1-xi(3)) &
                +(P(:,3)-P(:,4))*(1+xi(2))*(1-xi(3)) &
                +(P(:,6)-P(:,5))*(1-xi(2))*(1+xi(3)) &
                +(P(:,7)-P(:,8))*(1+xi(2))*(1+xi(3)))
          
dF(:,2)=0.125 * ((P(:,4)-P(:,1))*(1-xi(1))*(1-xi(3)) &
                +(P(:,3)-P(:,2))*(1+xi(1))*(1-xi(3)) &
                +(P(:,8)-P(:,5))*(1-xi(1))*(1+xi(3)) &
                +(P(:,7)-P(:,6))*(1+xi(1))*(1+xi(3)))

dF(:,3)=0.125 * ((P(:,5)-P(:,1))*(1-xi(1))*(1-xi(2)) &
                +(P(:,6)-P(:,2))*(1+xi(1))*(1-xi(2)) &
                +(P(:,7)-P(:,3))*(1+xi(1))*(1+xi(2)) &
                +(P(:,8)-P(:,4))*(1-xi(1))*(1+xi(2)))

! Determines the determinant of xj and checks for zero values
!
detjb = dF(1, 1) * dF(2, 2) * dF(3, 3) &
      + dF(1, 2) * dF(2, 3) * dF(3, 1) &
      + dF(1, 3) * dF(2, 1) * dF(3, 2) &
      - dF(1, 3) * dF(2, 2) * dF(3, 1) &
      - dF(1, 2) * dF(2, 1) * dF(3, 3) &
      - dF(1, 1) * dF(2, 3) * dF(3, 2)
IF ( detjb <= 0.d0 ) then
  WRITE(*,*)"Negative determinant of Jacobian in calc_df_inv:"
  WRITE(*,*)"dF",dF
  WRITE(*,*)"Determinant is:",detjb
  errorflag = 1
  !STOP 
END IF
!
! Determines the inverse of xj
!
dF_inv(1, 1) = (dF(2, 2) * dF(3, 3) &
              - dF(2, 3) * dF(3, 2) ) / detjb
dF_inv(1, 2) = (dF(1, 3) * dF(3, 2) &
              - dF(1, 2) * dF(3, 3) ) / detjb
dF_inv(1, 3) = (dF(1, 2) * dF(2, 3) &
              - dF(1, 3) * dF(2, 2) ) / detjb
dF_inv(2, 1) = (dF(2, 3) * dF(3, 1) &
              - dF(2, 1) * dF(3, 3) ) / detjb
dF_inv(2, 2) = (dF(1, 1) * dF(3, 3) &
              - dF(1, 3) * dF(3, 1) ) / detjb
dF_inv(2, 3) = (dF(1, 3) * dF(2, 1) &
              - dF(1, 1) * dF(2, 3) ) / detjb
dF_inv(3, 1) = (dF(2, 1) * dF(3, 2) &
              - dF(2, 2) * dF(3, 1) ) / detjb
dF_inv(3, 2) = (dF(1, 2) * dF(3, 1) &
              - dF(1, 1) * dF(3, 2) ) / detjb
dF_inv(3, 3) = (dF(1, 1) * dF(2, 2) &
               - dF(1, 2) * dF(2, 1) ) / detjb

Calc_dF_inv = dF_inv
END FUNCTION Calc_dF_inv 


FUNCTION Calc_KM_inv(P)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)          :: P(3,8)     ! 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                     :: Calc_KM_inv(3,3)  !  
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                     :: KM(3,3)  !  
REAL                     :: KM_inv(3,3)  !  
REAL                     :: detjb
!===================================================================================================================================
KM = 0.
KM_inv = 0.

KM(:,1) = - P(:,1) + P(:,2) + P(:,3) - P(:,4) &
          - P(:,5) + P(:,6) + P(:,7) - P(:,8) 
KM(:,2) = - P(:,1) - P(:,2) + P(:,3) + P(:,4) &
          - P(:,5) - P(:,6) + P(:,7) + P(:,8) 
KM(:,3) = - P(:,1) - P(:,2) - P(:,3) - P(:,4) &
          + P(:,5) + P(:,6) + P(:,7) + P(:,8) 

KM = 0.125 * KM
! Determines the determinant of xj and checks for zero values
!
detjb = KM (1, 1) * KM (2, 2) * KM (3, 3) &
      + KM (1, 2) * KM (2, 3) * KM (3, 1) &
      + KM (1, 3) * KM (2, 1) * KM (3, 2) &
      - KM (1, 3) * KM (2, 2) * KM (3, 1) &
      - KM (1, 2) * KM (2, 1) * KM (3, 3) &
      - KM (1, 1) * KM (2, 3) * KM (3, 2)
IF ( detjb <= 0.d0 ) then
  WRITE(*,*)"Negative determinant of Jacobian in Calc_KM_inv"
  WRITE(*,*)"Determinant is:",detjb
  WRITE(*,*)"KM:",KM
  STOP 
END IF
!
! Determines the inverse of xj
!
KM_inv (1, 1) = (KM (2, 2) * KM (3, 3) &
               - KM (2, 3) * KM (3, 2) ) / detjb
KM_inv (1, 2) = (KM (1, 3) * KM (3, 2) &
               - KM (1, 2) * KM (3, 3) ) / detjb
KM_inv (1, 3) = (KM (1, 2) * KM (2, 3) &
               - KM (1, 3) * KM (2, 2) ) / detjb
KM_inv (2, 1) = (KM (2, 3) * KM (3, 1) &
               - KM (2, 1) * KM (3, 3) ) / detjb
KM_inv (2, 2) = (KM (1, 1) * KM (3, 3) &
               - KM (1, 3) * KM (3, 1) ) / detjb
KM_inv (2, 3) = (KM (1, 3) * KM (2, 1) &
               - KM (1, 1) * KM (2, 3) ) / detjb
KM_inv (3, 1) = (KM (2, 1) * KM (3, 2) &
               - KM (2, 2) * KM (3, 1) ) / detjb
KM_inv (3, 2) = (KM (1, 2) * KM (3, 1) &
               - KM (1, 1) * KM (3, 2) ) / detjb
KM_inv (3, 3) = (KM (1, 1) * KM (2, 2) &
               - KM (1, 2) * KM (2, 1) ) / detjb

Calc_KM_inv = KM_inv
END FUNCTION Calc_KM_inv 


FUNCTION Calc_inv(M)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)          :: M(3,3)     ! 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                     :: Calc_inv(3,3)  !  
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                     :: M_inv(3,3)  !  
REAL                     :: detjb
!===================================================================================================================================
M_inv = 0.

! Determines the determinant of xj and checks for zero values
!
detjb = M (1, 1) * M (2, 2) * M (3, 3) &
      + M (1, 2) * M (2, 3) * M (3, 1) &
      + M (1, 3) * M (2, 1) * M (3, 2) &
      - M (1, 3) * M (2, 2) * M (3, 1) &
      - M (1, 2) * M (2, 1) * M (3, 3) &
      - M (1, 1) * M (2, 3) * M (3, 2)
IF ( detjb <= 0.d0 ) then
  WRITE(*,*)"Negative determinant of Jacobian in M_inv"
  WRITE(*,*)"Determinant is:",detjb
  WRITE(*,*)"KM:",M_inv
  STOP 
END IF
!
! Determines the inverse of xj
!
M_inv (1, 1) = (M (2, 2) * M (3, 3) &
              - M (2, 3) * M (3, 2) ) / detjb
M_inv (1, 2) = (M (1, 3) * M (3, 2) &
              - M (1, 2) * M (3, 3) ) / detjb
M_inv (1, 3) = (M (1, 2) * M (2, 3) &
              - M (1, 3) * M (2, 2) ) / detjb
M_inv (2, 1) = (M (2, 3) * M (3, 1) &
              - M (2, 1) * M (3, 3) ) / detjb
M_inv (2, 2) = (M (1, 1) * M (3, 3) &
              - M (1, 3) * M (3, 1) ) / detjb
M_inv (2, 3) = (M (1, 3) * M (2, 1) &
              - M (1, 1) * M (2, 3) ) / detjb
M_inv (3, 1) = (M (2, 1) * M (3, 2) &
              - M (2, 2) * M (3, 1) ) / detjb
M_inv (3, 2) = (M (1, 2) * M (3, 1) &
              - M (1, 1) * M (3, 2) ) / detjb
M_inv (3, 3) = (M (1, 1) * M (2, 2) &
              - M (1, 2) * M (2, 1) ) / detjb

Calc_inv = M_inv
END FUNCTION Calc_inv 


END MODULE MOD_Eval_xyz
