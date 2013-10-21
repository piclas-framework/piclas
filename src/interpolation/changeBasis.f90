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
INTERFACE ChangeBasis3D 
  MODULE PROCEDURE ChangeBasis3D
END INTERFACE

INTERFACE ChangeBasis2D 
  MODULE PROCEDURE ChangeBasis2D
END INTERFACE

PUBLIC :: ChangeBasis3D
PUBLIC :: ChangeBasis2D
!===================================================================================================================================

CONTAINS



SUBROUTINE ChangeBasis3D(Dim1,N_In,N_Out,Vdm,X3D_In,X3D_Out)
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
