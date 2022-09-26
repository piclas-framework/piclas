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

MODULE MOD_ChangeBasis
!===================================================================================================================================
! Changes a 2D or 3D Tensor Product Lagrange Points of Lagrange Basis of degree N_In to
! Lagrange points of a Lagrange Basis N_Out, using two
! arbitrary point disributions xi_In(0:N_In) and xi_Out(0:N_Out)
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
! no interface because dim on one input argument is one -> Arg(:,:) INPUT: Array(:,:,1)
!INTERFACE ChangeBasis3D
  !MODULE PROCEDURE ChangeBasis3D
!END INTERFACE

INTERFACE ChangeBasis3D_XYZ
  MODULE PROCEDURE ChangeBasis3D_XYZ
END INTERFACE

! no interface because dim on one input argument is one -> Arg(:,:) INPUT: Array(:,:,1)
!INTERFACE ChangeBasis2D
  !MODULE PROCEDURE ChangeBasis2D
!END INTERFACE

PUBLIC :: ChangeBasis3D
PUBLIC :: ChangeBasis2D
PUBLIC :: ChangeBasis3D_XYZ
!===================================================================================================================================

CONTAINS



PPURE SUBROUTINE ChangeBasis3D(Dim1,N_In,N_Out,Vdm,X3D_In,X3D_Out)
!===================================================================================================================================
! interpolate a 3D tensor product Lagrange basis defined by (N_in+1) 1D interpolation point positions xi_In(0:N_In)
! to another 3D tensor product node positions (number of nodes N_out+1)
! defined by (N_out+1) interpolation point  positions xi_Out(0:N_Out)
!  xi is defined in the 1DrefElem xi=[-1,1]
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: Dim1,N_In,N_Out
REAL,INTENT(IN)     :: X3D_In(1:Dim1,0:N_In,0:N_In,0:N_In)
REAL,INTENT(IN)     :: Vdm(0:N_Out,0:N_In)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: X3D_Out(1:Dim1,0:N_Out,0:N_Out,0:N_Out)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iN_In,jN_In,kN_In,iN_Out,jN_Out,kN_Out
REAL                :: X3D_Buf1(1:Dim1,0:N_Out,0:N_In,0:N_In)  ! first intermediate results from 1D interpolations
REAL                :: X3D_Buf2(1:Dim1,0:N_Out,0:N_Out,0:N_In) ! second intermediate results from 1D interpolations
!===================================================================================================================================
X3D_buf1=0.
! first direction iN_In
DO kN_In=0,N_In
  DO jN_In=0,N_In
    DO iN_In=0,N_In
      DO iN_Out=0,N_Out
        X3D_Buf1(:,iN_Out,jN_In,kN_In)=X3D_Buf1(:,iN_Out,jN_In,kN_In)+Vdm(iN_Out,iN_In)*X3D_In(:,iN_In,jN_In,kN_In)
      END DO
    END DO
  END DO
END DO
X3D_buf2=0.
! second direction jN_In
DO kN_In=0,N_In
  DO jN_In=0,N_In
    DO jN_Out=0,N_Out
      DO iN_Out=0,N_Out
        X3D_Buf2(:,iN_Out,jN_Out,kN_In)=X3D_Buf2(:,iN_Out,jN_Out,kN_In)+Vdm(jN_Out,jN_In)*X3D_Buf1(:,iN_Out,jN_In,kN_In)
      END DO
    END DO
  END DO
END DO
X3D_Out=0.
! last direction kN_In
DO kN_In=0,N_In
  DO kN_Out=0,N_Out
    DO jN_Out=0,N_Out
      DO iN_Out=0,N_Out
        X3D_Out(:,iN_Out,jN_Out,kN_Out)=X3D_Out(:,iN_Out,jN_Out,kN_Out)+Vdm(kN_Out,kN_In)*X3D_Buf2(:,iN_Out,jN_Out,kN_In)
      END DO
    END DO
  END DO
END DO
END SUBROUTINE ChangeBasis3D


SUBROUTINE ChangeBasis3D_XYZ(Dim1,NIn,NOut,Vdm_xi,Vdm_eta,Vdm_zeta,X3D_In,X3D_Out)
!==================================================================================================================================
!> Interpolate a 3D tensor product Lagrange polynomial defined by (NIn+1) 1D Lagrange basis functions of order (Nin) and node
!> positions xi_In(0:Nin) to another 3D tensor product Lagrange basis defined by (NOut+1) 1D interpolation points on the node
!> positions xi_out(0:NOut) using DIFFERENT 1D Vdm matrices in the xi,eta and zeta directions.
!> xi is defined in the 1D referent element \f$ \xi \in [-1,1] \f$.
!==================================================================================================================================
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: Dim1                                    !< Number of variables
INTEGER,INTENT(IN)  :: NIn                                     !< Input polynomial degree, no. of points = NIn+1
INTEGER,INTENT(IN)  :: NOut                                    !< Output polynomial degree, no. of points = NOut+1
REAL,INTENT(IN)     :: X3D_In(1:Dim1,0:NIn,0:NIn,0:NIn)        !< Input field, dimensions must match Dim1,NIn
REAL,INTENT(OUT)    :: X3D_Out(1:Dim1,0:NOut,0:NOut,0:NOut)    !< Output field, dimensions must match Dim1,NOut
REAL,INTENT(IN)     :: Vdm_xi(0:NOut,0:NIn)                    !< 1D Vandermonde In -> Out xi direction
REAL,INTENT(IN)     :: Vdm_eta(0:NOut,0:NIn)                   !< 1D Vandermonde In -> Out eta direction
REAL,INTENT(IN)     :: Vdm_zeta(0:NOut,0:NIn)                  !< 1D Vandermonde In -> Out zeta direction

!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iNIn,jNIn,kNIn,iN_Out,jN_Out,kN_Out
REAL                :: X3D_Buf1(1:Dim1,0:NOut,0:NIn,0:NIn)     ! first intermediate results from 1D interpolations
REAL                :: X3D_Buf2(1:Dim1,0:NOut,0:NOut,0:NIn)    ! second intermediate results from 1D interpolations
!==================================================================================================================================
X3D_buf1=0.
! first direction iNIn
DO kNIn=0,NIn
  DO jNIn=0,NIn
    DO iNIn=0,NIn
      DO iN_Out=0,NOut
        X3D_Buf1(:,iN_Out,jNIn,kNIn)=X3D_Buf1(:,iN_Out,jNIn,kNIn)+Vdm_xi(iN_Out,iNIn)*X3D_In(:,iNIn,jNIn,kNIn)
      END DO
    END DO
  END DO
END DO
X3D_buf2=0.
! second direction jNIn
DO kNIn=0,NIn
  DO jNIn=0,NIn
    DO jN_Out=0,NOut
      DO iN_Out=0,NOut
        X3D_Buf2(:,iN_Out,jN_Out,kNIn)=X3D_Buf2(:,iN_Out,jN_Out,kNIn)+Vdm_eta(jN_Out,jNIn)*X3D_Buf1(:,iN_Out,jNIn,kNIn)
      END DO
    END DO
  END DO
END DO
X3D_Out=0.
! last direction kNIn
DO kNIn=0,NIn
  DO kN_Out=0,NOut
    DO jN_Out=0,NOut
      DO iN_Out=0,NOut
        X3D_Out(:,iN_Out,jN_Out,kN_Out)=X3D_Out(:,iN_Out,jN_Out,kN_Out)+Vdm_zeta(kN_Out,kNIn)*X3D_Buf2(:,iN_Out,jN_Out,kNIn)
      END DO
    END DO
  END DO
END DO
END SUBROUTINE ChangeBasis3D_XYZ


SUBROUTINE ChangeBasis2D(Dim1,N_In,N_Out,Vdm,X2D_In,X2D_Out)
!===================================================================================================================================
! interpolate a 2D tensor product Lagrange basis defined by (N_in+1) 1D interpolation point positions xi_In(0:N_In)
! to another 2D tensor product node positions (number of nodes N_out+1)
! defined by (N_out+1) interpolation point  positions xi_Out(0:N_Out)
!  xi is defined in the 1DrefElem xi=[-1,1]
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: Dim1,N_In,N_Out
REAL,INTENT(IN)     :: X2D_In(1:Dim1,0:N_In,0:N_In)
REAL,INTENT(IN)     :: Vdm(0:N_Out,0:N_In)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: X2D_Out(1:Dim1,0:N_Out,0:N_Out)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iN_In,jN_In,iN_Out,jN_Out
REAL                :: X2D_Buf1(1:Dim1,0:N_Out,0:N_In)  ! first intermediate results from 1D interpolations
!===================================================================================================================================
X2D_buf1=0.
! first direction iN_In
DO jN_In=0,N_In
  DO iN_In=0,N_In
    DO iN_Out=0,N_Out
      X2D_Buf1(:,iN_Out,jN_In)=X2D_Buf1(:,iN_Out,jN_In)+Vdm(iN_Out,iN_In)*X2D_In(:,iN_In,jN_In)
    END DO
  END DO
END DO
X2D_Out=0.
! second direction jN_In
DO jN_In=0,N_In
  DO jN_Out=0,N_Out
    DO iN_Out=0,N_Out
      X2D_Out(:,iN_Out,jN_Out)=X2D_Out(:,iN_Out,jN_Out)+Vdm(jN_Out,jN_In)*X2D_Buf1(:,iN_Out,jN_In)
    END DO
  END DO
END DO
END SUBROUTINE ChangeBasis2D

END MODULE MOD_changeBasis
