#include "boltzplatz.h"

MODULE MOD_Riemann
!===================================================================================================================================
! Contains routines to compute the riemann (Advection, Diffusion) for a given Face
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

INTERFACE Riemann
  MODULE PROCEDURE Riemann
END INTERFACE

INTERFACE RiemannPML                  ! wird in src/dg/fillflux.f90 aufgerufen (zusaetzlich 24 Hilfsgroeßen)
  MODULE PROCEDURE RiemannPML
END INTERFACE

INTERFACE ExactFlux
  MODULE PROCEDURE ExactFlux
END INTERFACE

PUBLIC::Riemann,RiemannPML,Exactflux
!===================================================================================================================================

CONTAINS

SUBROUTINE Riemann(F,U_L,U_R,nv)
!===================================================================================================================================
! Computes the numerical flux
! Conservative States are rotated into normal direction in this routine and are NOT backrotated: don't use it after this routine!!
!===================================================================================================================================
! MODULES
USE MOD_PreProc ! PP_N
USE MOD_Equation_Vars,ONLY:eta_c,c,c2,c_corr,c_corr_c,c_corr_c2
USE MOD_Equation_Vars,ONLY:CentralFlux
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N),INTENT(IN) :: U_L,U_R
REAL,INTENT(IN)                                  :: nv(3,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                                 :: F(PP_nVar,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                                             :: n_loc(3),A_p(8,8),A_n(8,8), A(8,8)
INTEGER                                          :: Count_1,Count_2
!REAL                                             :: D(3,3)                  ! auxiliary matrices used 
!REAL                                             :: E(3,3), E_trans(3,3)    ! auxiliary matrices used
!===================================================================================================================================

IF(CentralFlux)THEN
  ! Gauss point i,j
  DO Count_2=0,PP_N
    DO Count_1=0,PP_N
      n_loc(:)=nv(:,Count_1,Count_2)
      A=0.
      ! check again 7 and 8
      A(5,1)=n_loc(3)
      A(6,1)=-n_loc(2)
      A(8,1)=c_corr*n_loc(1)
      A(4,2)=-n_loc(3)
      A(6,2)=n_loc(1)
      A(8,2)=c_corr*n_loc(2)
      A(4,3)=n_loc(2)
      A(5,3)=-n_loc(1)
      A(8,3)=c_corr*n_loc(3)

      A(2,4)=       -c2*n_loc(3)
      A(3,4)=        c2*n_loc(2)
      A(7,4)= c_corr_c2*n_loc(1)
      A(1,5)=        c2*n_loc(3)
      A(3,5)=       -c2*n_loc(1)
      A(7,5)= c_corr_c2*n_loc(2)
      A(1,6)=       -c2*n_loc(2)
      A(2,6)=        c2*n_loc(1)
      A(7,6)= c_corr_c2*n_loc(3)

      A(4,7)= c_corr*n_loc(1)
      A(5,7)= c_corr*n_loc(2)
      A(6,7)= c_corr*n_loc(3)

      A(1,8)= c_corr_c2*n_loc(1)
      A(2,8)= c_corr_c2*n_loc(2)
      A(3,8)= c_corr_c2*n_loc(3)

      F(:,Count_1,Count_2)=0.5*(MATMUL(A,U_R(:,Count_1,Count_2))+MATMUL(A,U_L(:,Count_1,Count_2)))
    END DO
  END DO
ELSE
  ! Gauss point i,j
  DO Count_2=0,PP_N
    DO Count_1=0,PP_N
      n_loc(:)=nv(:,Count_1,Count_2)
  
  !--- for original version see below (easier to understand)
  
      A_p(7,1:3)=0.
      A_p(1:3,7)=0.
      A_p(8,4:7)=0.
      A_p(4:7,8)=0.
    
      !D-Teilmatrix: Since chi and gamma is equal we
      ! consider D(chi,gamma) = D(gamma,chi)
      ! ATTENTION: if chi .ne. gamma this have to be changed. 
      ! Then we need D_1 and D_2 (see commented section below)
      A_p(1,1) = c + n_loc(1)*n_loc(1)*eta_c   !  D(1,1)=(1.+n_loc(1)*n_loc(1)*(eta-1.))*c
      A_p(1,2) = n_loc(1)*n_loc(2)*eta_c            !  D(1,2)=n_loc(1)*n_loc(2)*(eta-1)*c
      A_p(1,3) = n_loc(1)*n_loc(3)*eta_c            !  D(1,3)=n_loc(1)*n_loc(3)*(eta-1)*c
      A_p(2,1) = A_p(1,2)                          !  D(2,1)=n_loc(1)*n_loc(2)*(eta-1)*c
      A_p(2,2) = c + n_loc(2)*n_loc(2)*eta_c   !  D(2,2)=(1.+n_loc(2)*n_loc(2)*(eta-1.))*c
      A_p(2,3) = n_loc(2)*n_loc(3)*eta_c            !  D(2,3)=n_loc(2)*n_loc(3)*(eta-1)*c
      A_p(3,1) = A_p(1,3)                          !  D(3,1)=n_loc(1)*n_loc(3)*(eta-1)*c
      A_p(3,2) = A_p(2,3)                          !  D(3,2)=n_loc(2)*n_loc(3)*(eta-1)*c     
      A_p(3,3) = c+n_loc(3)*n_loc(3)*eta_c     !  D(3,3)=(1.+n_loc(3)*n_loc(3)*(mu-1.))*c
      ! epsilon-Teilmatrix
      !E_trans=transpose(E)
      A_p(1,4:6)= (/0.,c2*n_loc(3),-c2*n_loc(2)/)
      A_p(2,4:6)= (/-c2*n_loc(3),0.,c2*n_loc(1)/)
      A_p(3,4:6)= (/c2*n_loc(2),-c2*n_loc(1),0./)
      A_p(4,1:3)= (/0.,-n_loc(3),n_loc(2)/)
      A_p(5,1:3)= (/n_loc(3),0.,-n_loc(1)/)
      A_p(6,1:3)= (/-n_loc(2),n_loc(1),0./)
      !composition of the Matrix
      !positive A-Matrx
      A_p(4:6,4:6)=A_p(1:3,1:3)
      !negative A-Matrix
      A_n(1:3,1:3)=-A_p(1:3,1:3)
      A_n(1:3,4:6)= A_p(1:3,4:6)   ! c*c*E(:,:)
      A_n(4:6,1:3)= A_p(4:6,1:3)
      A_n(4:6,4:6)=-A_p(1:3,1:3)
    
     ! !positive A-Matrix-Divergence-Correction-Term
      A_p(1,8) = c_corr_c2*n_loc(1)
      A_p(2,8) = c_corr_c2*n_loc(2)
      A_p(3,8) = c_corr_c2*n_loc(3)
      A_p(4,7) = c_corr*n_loc(1)
      A_p(5,7) = c_corr*n_loc(2)
      A_p(6,7) = c_corr*n_loc(3)
      A_p(7,4) = c_corr_c2*n_loc(1)
      A_p(7,5) = c_corr_c2*n_loc(2)
      A_p(7,6) = c_corr_c2*n_loc(3)
      A_p(7,7) = c_corr_c
      A_p(8,1) = c_corr*n_loc(1)
      A_p(8,2) = c_corr*n_loc(2)
      A_p(8,3) = c_corr*n_loc(3)
      A_p(8,8) = c_corr_c
      !negative A-Matrix-Divergence-Correction-Term
      A_n(1:7,8) = A_p(1:7,8) !c_corr*c*c*n(1)
      A_n(1:6,7) = A_p(1:6,7) !c_corr*n(1)
      A_n(7,1:6)=  A_p(7,1:6)
      A_n(7,7)  = -A_p(7,7)
      A_n(8,1:7)= A_p(8,1:7)
      A_n(8,8)  =-A_p(8,8)
      ! Warum 0.5 -> Antwort im Taube/Dumbser-Paper. Im Munz/Schneider Paper fehlt das 0.5 lustigerweise.
  
  !--- Original Version:
  
  !    A_p=0.
  !    A_n=0.
  !  
  !    !D-Teilmatrix: Since chi and gamma is equal we
  !    ! consider D(chi,gamma) = D(gamma,chi)
  !    ! ATTENTION: if chi .ne. gamma this have to be changed. 
  !    ! Then we need D_1 and D_2
  !    D(1,1) = c + n_loc(1)*n_loc(1)*eta_c   !  D(1,1)=(1.+n_loc(1)*n_loc(1)*(eta-1.))*c
  !    D(1,2) = n_loc(1)*n_loc(2)*eta_c            !  D(1,2)=n_loc(1)*n_loc(2)*(eta-1)*c
  !    D(1,3) = n_loc(1)*n_loc(3)*eta_c            !  D(1,3)=n_loc(1)*n_loc(3)*(eta-1)*c
  !    D(2,1) = D(1,2)                          !  D(2,1)=n_loc(1)*n_loc(2)*(eta-1)*c
  !    D(2,2) = c + n_loc(2)*n_loc(2)*eta_c   !  D(2,2)=(1.+n_loc(2)*n_loc(2)*(eta-1.))*c
  !    D(2,3) = n_loc(2)*n_loc(3)*eta_c            !  D(2,3)=n_loc(2)*n_loc(3)*(eta-1)*c
  !    D(3,1) = D(1,3)                          !  D(3,1)=n_loc(1)*n_loc(3)*(eta-1)*c
  !    D(3,2) = D(2,3)                          !  D(3,2)=n_loc(2)*n_loc(3)*(eta-1)*c     
  !    D(3,3) = c+n_loc(3)*n_loc(3)*eta_c     !  D(3,3)=(1.+n_loc(3)*n_loc(3)*(mu-1.))*c
  !    ! epsilon-Teilmatrix
  !    !E_trans=transpose(E)
  !    E(1,:)= (/0.,n_loc(3),-n_loc(2)/)
  !    E(2,:)= (/-n_loc(3),0.,n_loc(1)/)
  !    E(3,:)= (/n_loc(2),-n_loc(1),0./)
  !    E_trans(1,:)= (/0.,-n_loc(3),n_loc(2)/)
  !    E_trans(2,:)= (/n_loc(3),0.,-n_loc(1)/)
  !    E_trans(3,:)= (/-n_loc(2),n_loc(1),0./)
  !    !composition of the Matrix
  !    !positive A-Matrx
  !    A_p(1:3,1:3)=D(:,:)
  !    A_p(1:3,4:6)=c2*E(:,:) ! c*c*E(:,:)
  !    A_p(4:6,1:3)=E_trans(:,:)
  !    A_p(4:6,4:6)=D(:,:)
  !    !negative A-Matrix
  !    A_n(1:3,1:3)=-D(:,:)
  !    A_n(1:3,4:6)= A_p(1:3,4:6)   ! c*c*E(:,:)
  !    A_n(4:6,1:3)= E_trans(:,:)
  !    A_n(4:6,4:6)=-D(:,:)
  !  
  !   ! !positive A-Matrix-Divergence-Correction-Term
  !    A_p(1,8) = c_corr_c2*n_loc(1)
  !    A_p(2,8) = c_corr_c2*n_loc(2)
  !    A_p(3,8) = c_corr_c2*n_loc(3)
  !    A_p(4,7) = c_corr*n_loc(1)
  !    A_p(5,7) = c_corr*n_loc(2)
  !    A_p(6,7) = c_corr*n_loc(3)
  !    A_p(7,4) = c_corr_c2*n_loc(1)
  !    A_p(7,5) = c_corr_c2*n_loc(2)
  !    A_p(7,6) = c_corr_c2*n_loc(3)
  !    A_p(7,7) = c_corr_c
  !    A_p(8,1) = c_corr*n_loc(1)
  !    A_p(8,2) = c_corr*n_loc(2)
  !    A_p(8,3) = c_corr*n_loc(3)
  !    A_p(8,8) = c_corr_c
  !    !negative A-Matrix-Divergence-Correction-Term
  !    A_n(1,8) = A_p(1,8) !c_corr*c*c*n(1)
  !    A_n(2,8) = A_p(2,8) !c_corr*c*c*n(2)
  !    A_n(3,8) = A_p(3,8) !c_corr*c*c*n(3)
  !    A_n(4,7) = A_p(4,7) !c_corr*n(1)
  !    A_n(5,7) = A_p(5,7) !c_corr*n(2)
  !    A_n(6,7) = A_p(6,7) !c_corr*n(3)
  !    !A_n(7,:)=(/0.,0.,0.,c_corr*c*c*n(1),c_corr*c*c*n(2),c_corr*c*c*n(3),-c_corr*c,0./)
  !    A_n(7,1:6)=  A_p(7,1:6)
  !    A_n(7,7)  = -A_p(7,7)
  !    !A_n(8,:)=(/c_corr*n(1),c_corr*n(2),c_corr*n(3),0.,0.,0.,0.,-c_corr*c/)
  !    A_n(8,1:7)= A_p(8,1:7)
  !    A_n(8,8)  =-A_p(8,8)
  !    ! Warum 0.5 -> Antwort im Taube/Dumbser-Paper. Im Munz/Schneider Paper fehlt das 0.5 lustigerweise.
      
      F(:,Count_1,Count_2)=0.5*(MATMUL(A_n,U_R(:,Count_1,Count_2))+MATMUL(A_p,U_L(:,Count_1,Count_2)))
    END DO
  END DO
END IF

END SUBROUTINE Riemann


SUBROUTINE RiemannPML(F,U_L,U_R,nv)
!===================================================================================================================================
! Computes the numerical flux
! Conservative States are rotated into normal direction in this routine and are NOT backrotatet: don't use it after this routine!!
!===================================================================================================================================
! MODULES
USE MOD_PreProc ! PP_N
USE MOD_Equation_Vars,ONLY:c,c2,c_corr,c_corr_c,c_corr_c2
USE MOD_PML_vars,     ONLY:PMLnVar
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N),INTENT(IN) :: U_L,U_R
REAL,INTENT(IN)                                  :: nv(3,0:PP_N,0:PP_N) !,t1(3,0:PP_N,0:PP_N),t2(3,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                                 :: F(PP_nVar+PMLnVar,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                                             :: n_loc(3)
INTEGER                                          :: Count_1,Count_2
!REAL                                             :: D(3,3)                  ! auxiliary matrices used 
!REAL                                             :: E(3,3), E_trans(3,3)    ! auxiliary matrices used
!===================================================================================================================================
! Gauss point i,j

DO Count_2=0,PP_N
  DO Count_1=0,PP_N
    n_loc(:)=nv(:,Count_1,Count_2)

    ! A^-*U_R + A^+*U_L
    !===============================================================================================================================
    !P1, P2, P3 and Ex
    F(9,Count_1,Count_2)  = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
          (/-c_corr_c*n_loc(1)*n_loc(1),-c_corr_c*n_loc(1)*n_loc(2),-c_corr_c*n_loc(1)*n_loc(3),0.,0.,0.,0.,c_corr_c2*n_loc(1)/) )+&
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
          (/ c_corr_c*n_loc(1)*n_loc(1), c_corr_c*n_loc(1)*n_loc(2), c_corr_c*n_loc(1)*n_loc(3),0.,0.,0.,0.,c_corr_c2*n_loc(1)/) ))
  
    F(10,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
          (/-c*n_loc(2)*n_loc(2), c*n_loc(1)*n_loc(2),0.,0.,0.,-c2*n_loc(2), c*n_loc(2)*n_loc(3),0./) )+ &
                  DOT_PRODUCT(U_L(:,Count_1,Count_2), &
          (/ c*n_loc(2)*n_loc(2),-c*n_loc(1)*n_loc(2),0.,0.,0.,-c2*n_loc(2),-c*n_loc(2)*n_loc(3),0./)))
 
    F(11,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
          (/-c*n_loc(3)*n_loc(3),0., c*n_loc(1)*n_loc(3),0.,c2*n_loc(3),0.,-c*n_loc(2)*n_loc(3),0./) )+ &
                  DOT_PRODUCT(U_L(:,Count_1,Count_2), &
          (/ c*n_loc(3)*n_loc(3),0.,-c*n_loc(1)*n_loc(3),0.,c2*n_loc(3),0., c*n_loc(2)*n_loc(3),0./)))
    !!==============================================================================================================================
    F(1,Count_1,Count_2)  = F(9,Count_1,Count_2) + F(10,Count_1,Count_2) + F(11,Count_1,Count_2)
    !!==============================================================================================================================
    !Q1, Q2, Q3 and Ey
    F(12,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
          (/ c*n_loc(1)*n_loc(2),-c*n_loc(1)*n_loc(1),0.,0.,0.,c2*n_loc(1),-c*n_loc(1)*n_loc(3),0./) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
          (/-c*n_loc(1)*n_loc(2), c*n_loc(1)*n_loc(1),0.,0.,0.,c2*n_loc(1), c*n_loc(1)*n_loc(3),0./) ))
  
    F(13,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
(/-c_corr_c*n_loc(1)*n_loc(2),-c_corr_c*n_loc(2)*n_loc(2),-c_corr_c*n_loc(2)*n_loc(3),0.,0.,0.,0.,c_corr_c2*n_loc(2)/) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
(/ c_corr_c*n_loc(1)*n_loc(2), c_corr_c*n_loc(2)*n_loc(2), c_corr_c*n_loc(2)*n_loc(3),0.,0.,0.,0.,c_corr_c2*n_loc(2)/)))
 
    F(14,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
          (/0.,-c*n_loc(3)*n_loc(3), c*n_loc(2)*n_loc(3),-c2*n_loc(3),0.,0., c*n_loc(1)*n_loc(3),0./) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
          (/0., c*n_loc(3)*n_loc(3),-c*n_loc(2)*n_loc(3),-c2*n_loc(3),0.,0.,-c*n_loc(1)*n_loc(3),0./)))
    !!==============================================================================================================================
    F(2,Count_1,Count_2)  = F(12,Count_1,Count_2) + F(13,Count_1,Count_2) + F(14,Count_1,Count_2)
    !!==============================================================================================================================
    !R1, R2, R3 and Ez
    F(15,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/ c*n_loc(1)*n_loc(3),0.,-c*n_loc(1)*n_loc(1),0.,-c2*n_loc(1),0., c*n_loc(1)*n_loc(2),0./) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/-c*n_loc(1)*n_loc(3),0., c*n_loc(1)*n_loc(1),0.,-c2*n_loc(1),0.,-c*n_loc(1)*n_loc(2),0./) ))
  
    F(16,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/0., c*n_loc(2)*n_loc(3),-c*n_loc(2)*n_loc(2),c2*n_loc(2),0.,0.,-c*n_loc(1)*n_loc(2),0./) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/0.,-c*n_loc(2)*n_loc(3), c*n_loc(2)*n_loc(2),c2*n_loc(2),0.,0., c*n_loc(1)*n_loc(2),0./)))
 
    F(17,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
(/-c_corr_c*n_loc(1)*n_loc(3),-c_corr_c*n_loc(2)*n_loc(3),-c_corr_c*n_loc(3)*n_loc(3),0.,0.,0.,0.,c_corr_c2*n_loc(3)/))+&
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
(/ c_corr_c*n_loc(1)*n_loc(3), c_corr_c*n_loc(2)*n_loc(3), c_corr_c*n_loc(3)*n_loc(3),0.,0.,0.,0.,c_corr_c2*n_loc(3)/)))
    !!==============================================================================================================================
    F(3,Count_1,Count_2)  = F(15,Count_1,Count_2) + F(16,Count_1,Count_2) + F(17,Count_1,Count_2)
    !!==============================================================================================================================
    !L1, L2, L3 and Bx
    F(18,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/0.,0.,0.,-c_corr_c*n_loc(1)*n_loc(1),-c_corr_c*n_loc(1)*n_loc(2),-c_corr_c*n_loc(1)*n_loc(3),c_corr*n_loc(1),0./) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/0.,0.,0., c_corr_c*n_loc(1)*n_loc(1), c_corr_c*n_loc(1)*n_loc(2), c_corr_c*n_loc(1)*n_loc(3),c_corr*n_loc(1),0./) ))
  
    F(19,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/0.,0., n_loc(2),-c*n_loc(2)*n_loc(2), c*n_loc(1)*n_loc(2),0.,0.,-c*n_loc(2)*n_loc(3)/) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/0.,0., n_loc(2), c*n_loc(2)*n_loc(2),-c*n_loc(1)*n_loc(2),0.,0., c*n_loc(2)*n_loc(3)/)))
 
    F(20,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/0.,-n_loc(3),0.,-c*n_loc(3)*n_loc(3),0., c*n_loc(1)*n_loc(3),0., c*n_loc(2)*n_loc(3)/) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/0.,-n_loc(3),0., c*n_loc(3)*n_loc(3),0.,-c*n_loc(1)*n_loc(3),0.,-c*n_loc(2)*n_loc(3)/)))
    !!==============================================================================================================================
    F(4,Count_1,Count_2)  = F(18,Count_1,Count_2) + F(19,Count_1,Count_2) + F(20,Count_1,Count_2)
    !!==============================================================================================================================
    !!M1, M2, M3 and By
    F(21,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/0.,0.,-n_loc(1), c*n_loc(1)*n_loc(2),-c*n_loc(1)*n_loc(1),0.,0., c*n_loc(1)*n_loc(3)/) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/0.,0.,-n_loc(1),-c*n_loc(1)*n_loc(2), c*n_loc(1)*n_loc(1),0.,0.,-c*n_loc(1)*n_loc(3)/) ))
  
    F(22,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/0.,0.,0.,-c_corr_c*n_loc(1)*n_loc(2),-c_corr_c*n_loc(2)*n_loc(2),-c_corr_c*n_loc(2)*n_loc(3),c_corr*n_loc(2),0./) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/0.,0.,0., c_corr_c*n_loc(1)*n_loc(2), c_corr_c*n_loc(2)*n_loc(2), c_corr_c*n_loc(2)*n_loc(3),c_corr*n_loc(2),0./)))
 
    F(23,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/n_loc(3),0.,0.,0.,-c*n_loc(3)*n_loc(3), c*n_loc(2)*n_loc(3),0.,-c*n_loc(1)*n_loc(3)/) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/n_loc(3),0.,0.,0., c*n_loc(3)*n_loc(3),-c*n_loc(2)*n_loc(3),0., c*n_loc(1)*n_loc(3)/)))
    !!==============================================================================================================================
    F(5,Count_1,Count_2)  = F(21,Count_1,Count_2) + F(22,Count_1,Count_2) + F(23,Count_1,Count_2)
    !!==============================================================================================================================
    !N1, N2, N3 and Bz
    F(24,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/0.,n_loc(1),0., c*n_loc(1)*n_loc(3),0.,-c*n_loc(1)*n_loc(1),0.,-c*n_loc(1)*n_loc(2)/) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/0.,n_loc(1),0.,-c*n_loc(1)*n_loc(3),0., c*n_loc(1)*n_loc(1),0., c*n_loc(1)*n_loc(2)/) ))
  
    F(25,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/-n_loc(2),0.,0.,0., c*n_loc(2)*n_loc(3),-c*n_loc(2)*n_loc(2),0., c*n_loc(1)*n_loc(2)/) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/-n_loc(2),0.,0.,0.,-c*n_loc(2)*n_loc(3), c*n_loc(2)*n_loc(2),0.,-c*n_loc(1)*n_loc(2)/)))
 
    F(26,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/0.,0.,0.,-c_corr_c*n_loc(1)*n_loc(3),-c_corr_c*n_loc(2)*n_loc(3),-c_corr_c*n_loc(3)*n_loc(3),c_corr*n_loc(3),0./) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/0.,0.,0., c_corr_c*n_loc(1)*n_loc(3), c_corr_c*n_loc(2)*n_loc(3), c_corr_c*n_loc(3)*n_loc(3),c_corr*n_loc(3),0./)))
    !!==============================================================================================================================
    F(6,Count_1,Count_2)  = F(24,Count_1,Count_2) + F(25,Count_1,Count_2) + F(26,Count_1,Count_2)
    !!==============================================================================================================================
    !S1, S2, S3 and PhiB
    F(27,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
    (/0., c_corr_c*n_loc(1)*n_loc(3),-c_corr_c*n_loc(1)*n_loc(2),c_corr_c2*n_loc(1),0.,0.,-c_corr_c*n_loc(1)*n_loc(1),0./) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
    (/0.,-c_corr_c*n_loc(1)*n_loc(3), c_corr_c*n_loc(1)*n_loc(2),c_corr_c2*n_loc(1),0.,0., c_corr_c*n_loc(1)*n_loc(1),0./) ))
  
    F(28,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
    (/-c_corr_c*n_loc(2)*n_loc(3),0., c_corr_c*n_loc(1)*n_loc(2),0.,c_corr_c2*n_loc(2),0.,-c_corr_c*n_loc(2)*n_loc(2),0./) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
    (/ c_corr_c*n_loc(2)*n_loc(3),0.,-c_corr_c*n_loc(1)*n_loc(2),0.,c_corr_c2*n_loc(2),0., c_corr_c*n_loc(2)*n_loc(2),0./)))
 
    F(29,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
    (/ c_corr_c*n_loc(2)*n_loc(3),-c_corr_c*n_loc(1)*n_loc(3),0.,0.,0.,c_corr_c2*n_loc(3),-c_corr_c*n_loc(3)*n_loc(3),0./) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
    (/-c_corr_c*n_loc(2)*n_loc(3), c_corr_c*n_loc(1)*n_loc(3),0.,0.,0.,c_corr_c2*n_loc(3), c_corr_c*n_loc(3)*n_loc(3),0./)))
    !!==============================================================================================================================
    F(7,Count_1,Count_2)  = F(27,Count_1,Count_2) + F(28,Count_1,Count_2) + F(29,Count_1,Count_2)
    !!==============================================================================================================================
    !T1, T2, T3 and PhiE
    F(30,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/c_corr*n_loc(1),0.,0.,0.,-c_corr_c*n_loc(1)*n_loc(3), c_corr_c*n_loc(1)*n_loc(2),0.,-c_corr_c*n_loc(1)*n_loc(1)/) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/c_corr*n_loc(1),0.,0.,0., c_corr_c*n_loc(1)*n_loc(3),-c_corr_c*n_loc(1)*n_loc(2),0., c_corr_c*n_loc(1)*n_loc(1)/) ))
  
    F(31,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/0.,c_corr*n_loc(2),0., c_corr_c*n_loc(2)*n_loc(3),0.,-c_corr_c*n_loc(1)*n_loc(2),0.,-c_corr_c*n_loc(2)*n_loc(2)/) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/0.,c_corr*n_loc(2),0.,-c_corr_c*n_loc(2)*n_loc(3),0., c_corr_c*n_loc(1)*n_loc(2),0., c_corr_c*n_loc(2)*n_loc(2)/)))
 
    F(32,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/0.,0.,c_corr*n_loc(3),-c_corr_c*n_loc(2)*n_loc(3), c_corr_c*n_loc(1)*n_loc(3),0.,0.,-c_corr_c*n_loc(3)*n_loc(3)/) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/0.,0.,c_corr*n_loc(3), c_corr_c*n_loc(2)*n_loc(3),-c_corr_c*n_loc(1)*n_loc(3),0.,0., c_corr_c*n_loc(3)*n_loc(3)/)))
    !!==============================================================================================================================
    !F(8) passt, durch ausgabe getestet
    F(8,Count_1,Count_2)  = F(30,Count_1,Count_2) + F(31,Count_1,Count_2) + F(32,Count_1,Count_2)
    !!==============================================================================================================================
  END DO
END DO

END SUBROUTINE RiemannPML


SUBROUTINE ExactFlux(t,tDeriv,Flux_Master,Flux_Slave,U_Master, U_slave,NormVec,Face_xGP,SurfElem)
!===================================================================================================================================
! Routine to add an exact function to a Riemann-Problem Face
! used at PML interfaces to emit a ave
! The ExactFlux is a non-conservative flux and is only emitted in ONE direction
! mapping
!               |
!       Master  |  Slave
!               |
!        nvec ----->
!               |
!               |
! Flux_Master   | Flux_Slave
! CAUTION: This routine has to multiply the new flux with the SurfElem
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_PreProc
USE MOD_Equation_Vars, ONLY: IniExactFunc,DoExactFlux,FluxDir,c,c2,c_corr,c_corr_c,c_corr_c2
USE MOD_Equation,      ONLY: ExactFunc
USE MOD_PML_Vars,      ONLY: xyzPhysicalMinMax
USE MOD_PML_vars,      ONLY: PMLnVar
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
REAL,INTENT(IN)       :: t           ! time
INTEGER,INTENT(IN)    :: tDeriv      ! deriv
REAL,INTENT(IN)       :: NormVec(1:3,0:PP_N,0:PP_N)
REAL,INTENT(IN)       :: Face_xGP(1:3,0:PP_N,0:PP_N)
REAL,INTENT(IN)       :: U_master(1:PP_nVar,0:PP_N,0:PP_N)
REAL,INTENT(IN)       :: U_slave (1:PP_nVar,0:PP_N,0:PP_N)
REAL,INTENT(IN)       :: SurfElem (0:PP_N,0:PP_N)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(INOUT)    :: Flux_Master(1:PP_nVar+PMLnVar,0:PP_N,0:PP_N)
REAL,INTENT(INOUT)    :: Flux_Slave (1:PP_nVar+PMLnVar,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: p,q
! assume exact solution is solution of RP
REAL                  :: U_Face_loc(1:PP_nVar)
REAL                  :: A(1:8,1:8)
REAL                  :: Flux_loc(1:PP_nVar)
REAL                  :: n_loc(1:3)
! assume, that the RP has to be solved with U_ex
!LOGICAL               :: UseMaster
!REAL                  :: U_Master_loc(1:PP_nVar,0:PP_N,0:PP_N)
!REAL                  :: U_Slave_loc (1:PP_nVar,0:PP_N,0:PP_N)
!REAL                  :: Flux_loc(1:PP_nVar+PMLnVar,0:PP_N,0:PP_N)
!===================================================================================================================================

IF(.NOT.DoExactFlux) RETURN

!UseMaster=.TRUE.
!! emission over plane, hence, first entry decides orientation of  plane
!IF(NormVec(FluxDir,0,0).GT.0)THEN
!  UseMaster=.FALSE.
!ELSE IF(NormVec(FluxDir,0,0).LT.0)THEN
!  UseMaster=.TRUE.
!ELSE
!  CALL abort(&
!__STAMP__&
!,'wired mesh?')
!END IF
!
!U_Slave_loc=0.
!U_Master_loc=0.
!DO q=0,PP_N
!  DO p=0,PP_N
!    IF(ALMOSTEQUALTOTOLERANCE(Face_xGP(FluxDir,p,q),xyzPhysicalMinMax(FluxDir*2-1),1e-4))THEN
!      ! the second state is always zero and already computed
!      IF(UseMaster)THEN
!        CALL ExactFunc(IniExactFunc,t,tDeriv,Face_xGP(:,p,q),U_Slave_loc(:,p,q))
!        !U_Master_loc(:,p,q)=U_Master(:,p,q)
!      ELSE
!        CALL ExactFunc(IniExactFunc,t,tDeriv,Face_xGP(:,p,q),U_Master_loc(:,p,q))
!        !U_Slave_loc(:,p,q)=U_Slave(:,p,q)
!      END IF
!    END IF
!  END DO ! p
!END DO ! q
!
!!CALL RiemannPML(Flux_loc(1:32,:,:),U_Master_loc(:,:,:),U_Slave_loc(:,:,:), NormVec(:,:,:))
!CALL Riemann(Flux_loc(1:8,:,:),U_Master_loc(:,:,:),U_Slave_loc(:,:,:), NormVec(:,:,:))
!IF(Usemaster)THEN
!  DO q=0,PP_N
!    DO p=0,PP_N
!      Flux_Master(:,p,q)=Flux_Master(:,p,q)+Flux_loc(:,p,q)*SurfElem(p,q)
!    END DO ! p
!  END DO ! q
!ELSE
!  DO q=0,PP_N
!    DO p=0,PP_N
!      Flux_Slave(:,p,q)=Flux_Slave(:,p,q)+Flux_loc(:,p,q)*SurfElem(p,q)
!    END DO ! p
!  END DO ! q
!END IF

! assume that U_ex is solution of RP at interface and add to flux
DO q=0,PP_N
  DO p=0,PP_N
    U_Face_loc=0.
    IF(ALMOSTEQUALTOTOLERANCE(Face_xGP(FluxDir,p,q),xyzPhysicalMinMax(FluxDir*2-1),1e-4))THEN
      CALL ExactFunc(IniExactFunc,t,tDeriv,Face_xGP(:,p,q),U_Face_loc)
      n_loc(1:3)=NormVec(1:3,p,q)
      A(1,1:4)=0.
      A(1, 5 )=c2*n_loc(3)
      A(1, 6 )=-c2*n_loc(2)
      A(1, 7 )=0.
      A(1, 8 )=c_corr_c2*n_loc(1)
      A(2,1:3)=0.
      A(2, 4 )=-c2*n_loc(3)
      A(2, 5 )=0.
      A(2, 6 )=c2*n_loc(1)
      A(2, 7 )=0.
      A(2, 8 )=c_corr_c2*n_loc(2)
      A(3,1:3)=0.
      A(3, 4 )=c2*n_loc(2)
      A(3, 5 )=-c2*n_loc(1)
      A(3, 6 )=0.
      A(3, 7 )=0.
      A(3, 8 )=c_corr_c2*n_loc(3)

      A(4, 1 )=0.
      A(4, 2 )=-n_loc(3)
      A(4, 3 )=n_loc(2)
      A(4,4:6)=0.
      A(4, 7 )=c_corr*n_loc(1)
      A(4, 8 )=0.

      A(5, 1 )=n_loc(3)
      A(5, 2 )=0.
      A(5, 3 )=-n_loc(1)
      A(5,4:6)=0.
      A(5, 7 )=c_corr*n_loc(2)
      A(5, 8 )=0.

      A(6, 1 )=-n_loc(2)
      A(6, 2 )=n_loc(1)
      A(6, 3 )=0.
      A(6,4:6)=0.
      A(6, 7 )=c_corr*n_loc(3)
      A(6, 8 )=0.

      A(7,1:3)=0.
      A(7, 4 )=c_corr_c2*n_loc(1)
      A(7, 5 )=c_corr_c2*n_loc(2)
      A(7, 6 )=c_corr_c2*n_loc(3)
      A(7, 7 )=0.
      A(7, 8 )=0.

      A(8, 1 )=c_corr*n_loc(1)
      A(8, 2 )=c_corr*n_loc(2)
      A(8, 3 )=c_corr*n_loc(3)
      A(8,4:6)=0.
      A(8, 7 )=0.
      A(8, 8 )=0.
      !Flux_loc=0.
      Flux_loc(1:PP_nVar)=MATMUL(A,U_Face_loc)
      ! mapping
      !         |
      ! Master  |  Slave
      !         |
      !  nvec ----->
      !         |
      !         |
      ! FluxM   | FluxSlave
      ! PO: are the signs correct?
      IF(NormVec(FluxDir,p,q).GT.0)THEN
        Flux_Slave(1:PP_nVar,p,q)=Flux_Slave(1:PP_nVar,p,q)+Flux_loc*SurfElem(p,q)
      ELSE
        Flux_Master(1:PP_nVar,p,q)=Flux_Master(1:PP_nVar,p,q)+Flux_loc*SurfElem(p,q) ! or sign change?
      END IF
    END IF
  END DO ! p
END DO ! q

END SUBROUTINE ExactFlux


END MODULE MOD_Riemann
