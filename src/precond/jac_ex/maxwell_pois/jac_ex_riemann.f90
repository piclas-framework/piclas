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

MODULE MOD_JacExRiemann
!===================================================================================================================================
! Contains routines to compute the jacobian of the riemann (Flux-Vector Splitting, Aplus, Aminus) for a given Face
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

INTERFACE ConstructJacRiemann
  MODULE PROCEDURE ConstructJacRiemann
END INTERFACE

INTERFACE ConstructJacRiemannDielectric
  MODULE PROCEDURE ConstructJacRiemannDielectric
END INTERFACE

INTERFACE ConstructJacBCRiemann
  MODULE PROCEDURE ConstructJacBCRiemann
END INTERFACE

INTERFACE ConstructJacNeighborRiemann
  MODULE PROCEDURE ConstructJacNeighborRiemann
END INTERFACE

PUBLIC::ConstructJacRiemann,ConstructJacBCRiemann,ConstructJacNeighborRiemann,ConstructJacRiemannDielectric

!===================================================================================================================================

CONTAINS

SUBROUTINE ConstructJacRiemann(nVec_loc,SurfElem_loc,Aside)
!===================================================================================================================================
! Computes the Aplus and the Aminus in the Flux-Vector Splittung of the numerical flux
! Conservative States are rotated into normal direction in this routine and are NOT backrotatet: don't use it after this routine!!
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals_Vars  ,ONLY: c,c2
USE MOD_Equation_Vars ,ONLY: eta_c,c_corr,c_corr_c,c_corr_c2
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                                  :: nVec_loc(1:3,0:PP_N,0:PP_N)
REAL,INTENT(IN)                                  :: SurfElem_loc(0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!REAL,DIMENSION(PP_nVar,PP_nVar,0:PP_N,0:PP_N),INTENT(OUT)  :: Aside
REAL,DIMENSION(8,8,0:PP_N,0:PP_N),INTENT(OUT)  :: Aside
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                             :: n_loc(3), A_p(8,8)
INTEGER                                          :: i,j
!===================================================================================================================================
! Gauss point i,j

  DO j=0 ,PP_N
    DO i =0,PP_N
      n_loc(:)=nVec_loc(:,i,j)

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

      ! calculate the contribution on the flux jacobian through the surface
      Aside (:,:,i,j)=0.5*A_p(:,:)*SurfElem_loc(i,j)
    END DO
  END DO

END SUBROUTINE ConstructJacRiemann


SUBROUTINE ConstructJacRiemannDielectric(nVec_loc,SurfElem_loc,Aside,locSideID,ElemID)
!===================================================================================================================================
! Computes the Aplus and the Aminus in the Flux-Vector Splittung of the numerical flux
! prolongation of the dielectric values to the side
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals_Vars       ,ONLY: c,c2
USE MOD_Equation_Vars      ,ONLY: eta_c,c_corr,c_corr_c,c_corr_c2
USE MOD_Dielectric_Vars    ,ONLY: DoDielectric,isDielectricElem, ElemToDielectric,DielectricConstant_inv
USE MOD_Interpolation_Vars ,ONLY: L_Minus,L_Plus
USE MOD_ProlongToFace      ,ONLY: ProlongToFace_Elementlocal
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                                  :: nVec_loc(1:3,0:PP_N,0:PP_N)
REAL,INTENT(IN)                                  :: SurfElem_loc(0:PP_N,0:PP_N)
INTEGER,INTENT(IN)                               :: locSideID,ElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(8,8,0:PP_N,0:PP_N),INTENT(OUT)    :: Aside
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                             :: n_loc(3), A_p(1:8,1:8)
INTEGER                                          :: p,q,l
REAL                                             :: DielectricConstant_inv_Face(0:PP_N,0:PP_N)
REAL                                             :: eta_c_dielectric,c_dielectric,c2_dielectric
!===================================================================================================================================


CALL ProlongToFace_Elementlocal(nVar=1                                                      &
                               ,locSideID=locSideID                                         &
                               ,Uvol=DielectricConstant_inv(:,:,:,ElemToDielectric(ElemID)) &
                               ,Uface=DielectricConstant_inv_Face                           )

! Gauss point p,q
DO q=0 ,PP_N
  DO p =0,PP_N
    ! normal vector of element
    n_loc(:)=nVec_loc(:,p,q)

    ! set dielectric values
    c_dielectric     = c*DielectricConstant_inv_Face(p,q)            !         c/sqrt(EpsR*MuR)
    eta_c_dielectric = c_corr_c - c_dielectric                       ! (chi - 1./sqrt(EpsR*MuR))*c
    c2_dielectric    = c_dielectric*c_dielectric                     !          c**2/(EpsR*MuR)

    !--- for original version see below (easier to understand)
    A_p(7,1:3)=0.
    A_p(1:3,7)=0.
    A_p(8,4:7)=0.
    A_p(4:7,8)=0.

    !D-Teilmatrix: Since chi and gamma is equal we
    ! consider D(chi,gamma) = D(gamma,chi)
    ! ATTENTION: if chi .ne. gamma this have to be changed.
    ! Then we need D_1 and D_2 (see commented section below)
    A_p(1,1) = c_dielectric + n_loc(1)*n_loc(1)*eta_c        !  D(1,1)=(1./sqrt(EpsR*MuR)+n_loc(1)*n_loc(1)*(chi-1.))*c
    A_p(1,2) = n_loc(1)*n_loc(2)*eta_c_dielectric            !  D(1,2)=n_loc(1)*n_loc(2)*(chi-1./sqrt(EpsR*MuR))*c
    A_p(1,3) = n_loc(1)*n_loc(3)*eta_c_dielectric            !  D(1,3)=n_loc(1)*n_loc(3)*(chi-1./sqrt(EpsR*MuR))*c
    A_p(2,1) = A_p(1,2)                                      !  D(2,1)=n_loc(1)*n_loc(2)*(chi-1./sqrt(EpsR*MuR))*c
    A_p(2,2) = c_dielectric + n_loc(2)*n_loc(2)*eta_c        !  D(2,2)=(1./sqrt(EpsR*MuR)+n_loc(2)*n_loc(2)*(chi-1.))*c
    A_p(2,3) = n_loc(2)*n_loc(3)*eta_c_dielectric            !  D(2,3)=n_loc(2)*n_loc(3)*(chi-1./sqrt(EpsR*MuR))*c
    A_p(3,1) = A_p(1,3)                                      !  D(3,1)=n_loc(1)*n_loc(3)*(chi-1./sqrt(EpsR*MuR))*c
    A_p(3,2) = A_p(2,3)                                      !  D(3,2)=n_loc(2)*n_loc(3)*(chi-1./sqrt(EpsR*MuR))*c
    A_p(3,3) = c_dielectric+n_loc(3)*n_loc(3)*eta_c          !  D(3,3)=(1./sqrt(EpsR*MuR)+n_loc(3)*n_loc(3)*(mu-1.))*c
    ! epsilon-Teilmatrix
    !E_trans=transpose(E)
    A_p(1,4:6)= (/0.,c2_dielectric*n_loc(3),-c2_dielectric*n_loc(2)/)
    A_p(2,4:6)= (/-c2_dielectric*n_loc(3),0.,c2_dielectric*n_loc(1)/)
    A_p(3,4:6)= (/c2_dielectric*n_loc(2),-c2_dielectric*n_loc(1),0./)
    A_p(4,1:3)= (/0.,-n_loc(3),n_loc(2)/)
    A_p(5,1:3)= (/n_loc(3),0.,-n_loc(1)/)
    A_p(6,1:3)= (/-n_loc(2),n_loc(1),0./)
    !composition of the Matrix
    !positive A-Matrx
    A_p(4:6,4:6)=A_p(1:3,1:3)

    !positive A-Matrix-Divergence-Correction-Term
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


    ! calculate the contribution on the flux jacobian through the surface
    Aside (:,:,p,q)=0.5*A_p(:,:)*SurfElem_loc(p,q)
  END DO
END DO

END SUBROUTINE ConstructJacRiemannDielectric


SUBROUTINE ConstructJacNeighborRiemann(nVec_loc,SurfElem_loc,Aside)
!===================================================================================================================================
! Computes the Aplus and the Aminus in the Flux-Vector Splittung of the numerical flux
! Conservative States are rotated into normal direction in this routine and are NOT backrotatet: don't use it after this routine!!
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals_Vars  ,ONLY: c,c2
USE MOD_Equation_Vars ,ONLY: eta_c,c_corr,c_corr_c,c_corr_c2
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                                  :: nVec_loc(1:3,0:PP_N,0:PP_N)
REAL,INTENT(IN)                                  :: SurfElem_loc(0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(8,8,0:PP_N,0:PP_N),INTENT(OUT)  :: Aside
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                             :: n_loc(3), A_m(8,8)
INTEGER                                          :: i,j
!===================================================================================================================================
! Gauss point i,j

  DO j=0 ,PP_N
    DO i =0,PP_N
      n_loc(:)=nVec_loc(:,i,j)

!--- for original version see below (easier to understand)

      A_m(7,1:3)=0.
      A_m(1:3,7)=0.
      A_m(8,4:7)=0.
      A_m(4:7,8)=0.

     !D-Teilmatrix: Since chi and gamma is equal we
     ! consider D(chi,gamma) = D(gamma,chi)
     ! ATTENTION: if chi .ne. gamma this have to be changed.
     ! Then we need D_1 and D_2 (see commented section below)
      A_m(1,1) =-1.0*(c + n_loc(1)*n_loc(1)*eta_c ) !  D(1,1)=(1.+n_loc(1)*n_loc(1)*(eta-1.))*c
      A_m(1,2) =-1.0*(n_loc(1)*n_loc(2)*eta_c     )      !  D(1,2)=n_loc(1)*n_loc(2)*(eta-1)*c
      A_m(1,3) =-1.0*(n_loc(1)*n_loc(3)*eta_c     )      !  D(1,3)=n_loc(1)*n_loc(3)*(eta-1)*c
      A_m(2,1) =      A_m(1,2)                         !  D(2,1)=n_loc(1)*n_loc(2)*(eta-1)*c
      A_m(2,2) =-1.0*(c + n_loc(2)*n_loc(2)*eta_c ) !  D(2,2)=(1.+n_loc(2)*n_loc(2)*(eta-1.))*c
      A_m(2,3) =-1.0*(n_loc(2)*n_loc(3)*eta_c     )      !  D(2,3)=n_loc(2)*n_loc(3)*(eta-1)*c
      A_m(3,1) =      A_m(1,3)                         !  D(3,1)=n_loc(1)*n_loc(3)*(eta-1)*c
      A_m(3,2) =      A_m(2,3)                         !  D(3,2)=n_loc(2)*n_loc(3)*(eta-1)*c
      A_m(3,3) =-1.0*(c+n_loc(3)*n_loc(3)*eta_c   ) !  D(3,3)=(1.+n_loc(3)*n_loc(3)*(mu-1.))*c
      ! epsilon-Teilmatrix
      !E_trans=transpose(E)
      A_m(1,4:6)= (/0.,c2*n_loc(3),-c2*n_loc(2)/)
      A_m(2,4:6)= (/-c2*n_loc(3),0.,c2*n_loc(1)/)
      A_m(3,4:6)= (/c2*n_loc(2),-c2*n_loc(1),0./)
      A_m(4,1:3)= (/0.,-n_loc(3),n_loc(2)/)
      A_m(5,1:3)= (/n_loc(3),0.,-n_loc(1)/)
      A_m(6,1:3)= (/-n_loc(2),n_loc(1),0./)
      !composition of the Matrix
      !negative A-Matrx
      A_m(4:6,4:6)=A_m(1:3,1:3)

     ! !negative A-Matrix-Divergence-Correction-Term
      A_m(1,8) = c_corr_c2*n_loc(1)
      A_m(2,8) = c_corr_c2*n_loc(2)
      A_m(3,8) = c_corr_c2*n_loc(3)
      A_m(4,7) = c_corr*n_loc(1)
      A_m(5,7) = c_corr*n_loc(2)
      A_m(6,7) = c_corr*n_loc(3)
      A_m(7,4) = c_corr_c2*n_loc(1)
      A_m(7,5) = c_corr_c2*n_loc(2)
      A_m(7,6) = c_corr_c2*n_loc(3)
      A_m(7,7) = -1.0*c_corr_c
      A_m(8,1) = c_corr*n_loc(1)
      A_m(8,2) = c_corr*n_loc(2)
      A_m(8,3) = c_corr*n_loc(3)
      A_m(8,8) = -1.0*c_corr_c

      ! calculate the contribution on the flux jacobian through the surface
      Aside (:,:,i,j)=0.5*A_m(:,:)*SurfElem_loc(i,j)
    END DO
  END DO

END SUBROUTINE ConstructJacNeighborRiemann

SUBROUTINE ConstructJacBCRiemann(BCType,nVec_loc,SurfElem_loc,Aside)
!===================================================================================================================================
! Computes the jacobin of the Am of the numerical flux
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals_Vars  ,ONLY: c,c2
USE MOD_Equation_Vars ,ONLY: c_corr,c_corr_c,c_corr_c2
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                               :: BCType
REAL,INTENT(IN)                                  :: nVec_loc(1:3,0:PP_N,0:PP_N)
REAL,INTENT(IN)                                  :: SurfElem_loc(0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(8,8,0:PP_N,0:PP_N),INTENT(OUT)  :: Aside
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                             :: n_loc(3), A(8,8)
INTEGER                                          :: i,j
REAL                                             :: c_corr_c_pc
!===================================================================================================================================

SELECT CASE(BCType)
  CASE(1) !Periodic already filled!
    ! has to be nullyfied
    Aside(:,:,:,:) = 0.
  CASE(2) ! exact BC = Dirichlet BC !!
  ! as long as it does not depend on interior value
    ! has to be nullyfied as long as it does not depend on inner values
    Aside(:,:,:,:) = 0.
  CASE(3,5) ! 1st order absorbing BC
          ! Silver-Mueller BC - Munz et al. 2000 / Computer Physics Communication 130, 83-117

! Gauss point i,j
    DO j=0 ,PP_N
      DO i =0,PP_N
        n_loc(:)=nVec_loc(:,i,j)

        A(1:3,4:7)=0.
        A(4:7,1:3)=0.
        !A(4:7,4:7)=0.
        A(8,4:7)  =0.
        A(4:7,8)  =0.

        !D-Teilmatrix: Since chi and gamma is equal we
        ! consider D(chi,gamma) = D(gamma,chi)
        ! ATTENTION: if chi .ne. gamma this have to be changed.
        ! Then we need D_1 and D_2 (see commented section below)
        A(1,1) = -n_loc(1)*n_loc(1)*c_corr_c   !  D(1,1)=-n_loc(1)*n_loc(1)*c_corr_c
        A(1,2) = -n_loc(1)*n_loc(2)*c_corr_c   !  D(1,2)=-n_loc(1)*n_loc(2)*c_corr_c
        A(1,3) = -n_loc(1)*n_loc(3)*c_corr_c   !  D(1,3)=-n_loc(1)*n_loc(3)*c_corr_c
        A(2,1) =  A(1,2)                       !  D(2,1)=-n_loc(1)*n_loc(2)*c_corr_c
        A(2,2) = -n_loc(2)*n_loc(2)*c_corr_c   !  D(2,2)=-n_loc(2)*n_loc(2)*c_corr_c
        A(2,3) = -n_loc(2)*n_loc(3)*c_corr_c   !  D(2,3)=-n_loc(2)*n_loc(3)*c_corr_c
        A(3,1) =  A(1,3)                       !  D(3,1)=-n_loc(1)*n_loc(3)*c_corr_c
        A(3,2) =  A(2,3)                       !  D(3,2)=-n_loc(2)*n_loc(3)*c_corr_c
        A(3,3) = -n_loc(3)*n_loc(3)*c_corr_c   !  D(3,3)=-n_loc(3)*n_loc(3)*c_corr_c

        !composition of the Matrix
        !positive A-Matrx
        A(4:6,4:6)=A(1:3,1:3)

        ! A-Matrix-Divergence-Correction-Term
        A(1,8) =-c_corr_c2*n_loc(1)
        A(2,8) =-c_corr_c2*n_loc(2)
        A(3,8) =-c_corr_c2*n_loc(3)
        A(4,7) =-c_corr*n_loc(1)
        A(5,7) =-c_corr*n_loc(2)
        A(6,7) =-c_corr*n_loc(3)
        A(7,4) = c_corr_c2*n_loc(1)
        A(7,5) = c_corr_c2*n_loc(2)
        A(7,6) = c_corr_c2*n_loc(3)
        A(7,7) = c_corr_c
        A(8,1) = c_corr*n_loc(1)
        A(8,2) = c_corr*n_loc(2)
        A(8,3) = c_corr*n_loc(3)
        A(8,8) = c_corr_c

        Aside (:,:,i,j)=0.5*A(:,:)*SurfElem_loc(i,j)
    END DO
  END DO

CASE(4) ! perfectly conducting surface (MunzOmnesSchneider 2000, pp. 97-98)
  ! Determine the exact BC state
  ! jacobin of the influx
  c_corr_c_pc = c_corr_c+c
  DO j=0,PP_N
    DO i=0,PP_N
        n_loc(:)=nVec_loc(:,i,j)
        ! zeros
        A=0.
        A(1:3,7) = 0.
        A(4:7,8) = 0.
        A(7,1:3) = 0.
        A(8,4:7) = 0.
        ! first diagonal 3x3 block
        A(1,1) = c-c_corr_c_pc*n_loc(1)*n_loc(1)
        A(2,2) = c-c_corr_c_pc*n_loc(2)*n_loc(2)
        A(3,3) = c-c_corr_c_pc*n_loc(3)*n_loc(3)
        A(1,2) =  -c_corr_c_pc*n_loc(1)*n_loc(2)
        A(1,3) =  -c_corr_c_pc*n_loc(1)*n_loc(3)
        A(2,1) = A(1,2)
        A(2,3) =  -c_corr_c_pc*n_loc(2)*n_loc(3)
        A(3,1) = A(1,3)
        A(3,2) = A(2,3)
        ! second diagonal 3x3 block
        A(4,4) = c_corr_c_pc*n_loc(1)*n_loc(1)-c
        A(5,5) = c_corr_c_pc*n_loc(2)*n_loc(2)-c
        A(6,6) = c_corr_c_pc*n_loc(3)*n_loc(3)-c
        A(4,5) = c_corr_c_pc*n_loc(1)*n_loc(2)
        A(4,6) = c_corr_c_pc*n_loc(1)*n_loc(3)
        A(5,4) = A(4,5)
        A(5,6) = c_corr_c_pc*n_loc(2)*n_loc(3)
        A(6,4) = A(4,6)
        A(6,5) = A(5,6)
        ! lower 3x3 block
        A(4,1:3) = (/0.,n_loc(3),-n_loc(2)/)
        A(5,1:3) = (/-n_loc(3),0.,n_loc(1)/)
        A(6,1:3) = (/n_loc(2),-n_loc(1),0./)
        ! left 3x3 block
        A(1,4:6) = (/0.,c2*n_loc(3),-c2*n_loc(2)/)
        A(2,4:6) = (/-c2*n_loc(3),0.,c2*n_loc(1)/)
        A(3,4:6) = (/ c2*n_loc(2),-c2*n_loc(1),0./)
         ! last terms | divergence correction
        A(1,8) = -c_corr_c2*n_loc(1)
        A(2,8) = -c_corr_c2*n_loc(2)
        A(3,8) = -c_corr_c2*n_loc(3)
        A(4,7) = c_corr*n_loc(1)
        A(5,7) = c_corr*n_loc(2)
        A(6,7) = c_corr*n_loc(3)
        A(7,7) = -c_corr_c
        A(7,4) = -c_corr_c2*n_loc(1)
        A(7,5) = -c_corr_c2*n_loc(2)
        A(7,6) = -c_corr_c2*n_loc(3)
        A(8,1) = c_corr*n_loc(1)
        A(8,2) = c_corr*n_loc(2)
        A(8,3) = c_corr*n_loc(3)
        A(8,8) = c_corr_c

        Aside (:,:,i,j)=0.5*A(:,:)*SurfElem_loc(i,j)
    END DO ! p
  END DO ! q

  !CASE(5) ! open boundary condition with no force of divergenz cleaning
  !  ! has to be nullyfied
  !  Aside(:,:,:,:) = 0.


END SELECT

END SUBROUTINE ConstructJacBCRiemann


END MODULE MOD_JacExRiemann
