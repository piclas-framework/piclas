!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (gitlab.com/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

MODULE MOD_Eval_xyz
!===================================================================================================================================
!> Contains routines for transformation from reference to physical space and vice versa:
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

INTERFACE GetPositionInRefElem
  MODULE PROCEDURE GetPositionInRefElem
END INTERFACE

INTERFACE TensorProductInterpolation
  MODULE PROCEDURE TensorProductInterpolation
END INTERFACE

INTERFACE EvaluateFieldAtPhysPos
  MODULE PROCEDURE EvaluateFieldAtPhysPos
END INTERFACE

INTERFACE EvaluateFieldAtRefPos
  MODULE PROCEDURE EvaluateFieldAtRefPos
END INTERFACE

PUBLIC :: GetPositionInRefElem
PUBLIC :: TensorProductInterpolation
PUBLIC :: EvaluateFieldAtPhysPos
PUBLIC :: EvaluateFieldAtRefPos
!===================================================================================================================================

CONTAINS

SUBROUTINE EvaluateFieldAtPhysPos(x_in,NVar,N_in,U_In,U_Out,ElemID,PartID)
!===================================================================================================================================
!> 1) Get position within reference element (x_in -> xi=[-1,1]) by inverting the mapping
!> 2) interpolate DG solution to position (U_In -> U_Out(x_in))
!> 3) interpolate backgroundfield to position ( U_Out -> U_Out(x_in)+BG_field(x_in) )
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Basis,                   ONLY:LagrangeInterpolationPolys
USE MOD_Interpolation_Vars,      ONLY:xGP,wBary
USE MOD_Mesh_Vars,               ONLY:dXCL_NGeo,XCL_NGeo,NGeo,wBaryCL_NGeo,XiCL_NGeo
USE MOD_PICInterpolation_Vars,   ONLY:NBG,BGField,useBGField,BGDataSize,BGField_wBary, BGField_xGP,BGType
USE MOD_Mesh_Vars,               ONLY:CurvedElem,wBaryCL_NGeo1,XiCL_NGeo1
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: NVar                                  !< 6 (Ex, Ey, Ez, Bx, By, Bz)
INTEGER,INTENT(IN)  :: N_In                                  !< usually PP_N
INTEGER,INTENT(IN)  :: ElemID                                !< Element index
REAL,INTENT(IN)     :: U_In(1:NVar,0:N_In,0:N_In,0:N_In)     !< State in Element
REAL,INTENT(IN)     :: x_in(3)                               !< position in physical space
INTEGER,INTENT(IN),OPTIONAL :: PartID                        !< particle ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: U_Out(1:NVar)                         !< Interpolated state at physical position x_in
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k
REAL                :: xi(3)
REAL, PARAMETER     :: EPSONE=1.00000001
REAL                :: L_xi(3,0:PP_N), L_eta_zeta
REAL                :: XCL_NGeo1(1:3,0:1,0:1,0:1)
REAL                :: dXCL_NGeo1(1:3,1:3,0:1,0:1,0:1)
! h5-external e,b field
REAL,ALLOCATABLE    :: L_xi_BGField(:,:), U_BGField(:)
!===================================================================================================================================

CALL GetRefNewtonStartValue(X_in,Xi,ElemID)

IF(CurvedElem(ElemID))THEN
  CALL RefElemNewton(Xi,X_In,wBaryCL_NGeo,XiCL_NGeo,XCL_NGeo(:,:,:,:,ElemID),dXCL_NGeo(:,:,:,:,:,ElemID) &
                    ,NGeo,ElemID,Mode=1,PartID=PartID)
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
  CALL RefElemNewton(Xi,X_In,wBaryCL_NGeo1,XiCL_NGeo1,XCL_NGeo1,dXCL_NGeo1,1,ElemID,Mode=1,PartID=PartID)
END IF

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

END SUBROUTINE EvaluateFieldAtPhysPos


SUBROUTINE GetPositionInRefElem(x_in,xi,ElemID,DoReUseMap,ForceMode)
!===================================================================================================================================
!> Get Position within reference element (x_in -> xi=[-1,1])
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Basis,                   ONLY:LagrangeInterpolationPolys
USE MOD_Mesh_Vars,               ONLY:dXCL_NGeo,XCL_NGeo,NGeo,wBaryCL_NGeo,XiCL_NGeo,NGeo
USE MOD_Mesh_Vars,               ONLY:CurvedElem,wBaryCL_NGeo1,XiCL_NGeo1
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)          :: ElemID                                 !< element index
REAL,INTENT(IN)             :: x_in(3)                                !< position in physical space
LOGICAL,INTENT(IN),OPTIONAL :: DoReUseMap                             !< flag if start values for newton elem mapping already exists
LOGICAL,INTENT(IN),OPTIONAL :: ForceMode                              !< flag for mode change in RefElemNewton
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)          :: xi(1:3)                                !< position in reference element
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                    :: iMode
REAL                       :: XCL_NGeo1(1:3,0:1,0:1,0:1)
REAL                       :: dXCL_NGeo1(1:3,1:3,0:1,0:1,0:1)
!===================================================================================================================================

iMode=2
IF(PRESENT(ForceMode)) iMode=1
IF(.NOT.PRESENT(DoReUseMap))THEN
  CALL GetRefNewtonStartValue(X_in,Xi,ElemID)
END IF

IF(CurvedElem(ElemID))THEN
  CALL RefElemNewton(Xi,X_In,wBaryCL_NGeo,XiCL_NGeo,XCL_NGeo(:,:,:,:,ElemID),dXCL_NGeo(:,:,:,:,:,ElemID),NGeo,ElemID,Mode=iMode)
ELSE
  ! fill dummy XCL_NGeo1
  IF(NGeo.EQ.1)THEN
    CALL RefElemNewton(Xi,X_In,wBaryCL_NGeo,XiCL_NGeo,XCL_NGeo(:,:,:,:,ElemID),dXCL_NGeo(:,:,:,:,:,ElemID),NGeo,ElemID,Mode=iMode)
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
    CALL RefElemNewton(Xi,X_In,wBaryCL_NGeo1,XiCL_NGeo1,XCL_NGeo1,dXCL_NGeo1,1,ElemID,Mode=iMode)
  END IF
END IF

END SUBROUTINE GetPositionInRefElem


SUBROUTINE TensorProductInterpolation(Xi_in,NVar,N_in,xGP_in,wBary_In,U_In,U_Out)
!===================================================================================================================================
!> Interpolates a 3D tensor product Lagrange basis defined by (N_in+1) 1D interpolation points to the position Xi
!===================================================================================================================================
! MODULES
USE MOD_Basis,                 ONLY: LagrangeInterpolationPolys
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: NVar                                  !< 6 (Ex, Ey, Ez, Bx, By, Bz)
INTEGER,INTENT(IN)  :: N_In                                  !< usually PP_N
REAL,INTENT(IN)     :: U_In(1:NVar,0:N_In,0:N_In,0:N_In)     !< State in Element
REAL,INTENT(IN)     :: Xi_in(3)                              !< position in reference element
REAL,INTENT(IN)     :: xGP_In(0:N_in)
REAL,INTENT(IN)     :: wBary_In(0:N_in)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: U_Out(1:NVar)                         !< Interpolated state at reference position xi_in
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: i,j,k
REAL,DIMENSION(3,0:N_in)  :: L_xi
REAL                      :: L_eta_zeta
!===================================================================================================================================
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
END SUBROUTINE TensorProductInterpolation


SUBROUTINE EvaluateFieldAtRefPos(xi_in,NVar,N_in,U_In,U_Out,ElemID)
!===================================================================================================================================
!> 1) interpolate DG solution to position (U_In -> U_Out(xi_in))
!> 2) interpolate backgroundfield to position ( U_Out -> U_Out(xi_in)+BG_field(xi_in) )
!===================================================================================================================================
! MODULES
USE MOD_Basis,                 ONLY: LagrangeInterpolationPolys
USE MOD_Interpolation_Vars,    ONLY: wBary,xGP
USE MOD_PICInterpolation_Vars, ONLY:NBG,BGField,useBGField,BGDataSize,BGField_xGP,BGField_wBary,BGType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: NVar                                  !< 6 (Ex, Ey, Ez, Bx, By, Bz)
INTEGER,INTENT(IN)  :: N_In                                  !< usually PP_N
INTEGER,INTENT(IN)  :: ElemID                                !< Element index
REAL,INTENT(IN)     :: U_In(1:NVar,0:N_In,0:N_In,0:N_In)     !< State in Element
REAL,INTENT(IN)     :: Xi_in(3)                              !< position in reference element
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: U_Out(1:NVar)                         !< Interpolated state at reference position xi_in
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


END SUBROUTINE EvaluateFieldAtRefPos


SUBROUTINE RefElemNewton(Xi,X_In,wBaryCL_N_In,XiCL_N_In,XCL_N_In,dXCL_N_In,N_In,ElemID,Mode,PartID)
!=================================================================================================================================
!> Netwon for finding the position inside the reference element [-1,1] for an arbitrary physical point
!=================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_Basis,                   ONLY:LagrangeInterpolationPolys
USE MOD_Particle_Mesh_Vars,      ONLY:RefMappingEps
USE MOD_Mesh_Vars,               ONLY:offsetElem
#if defined(IMPA)
USE MOD_Particle_Vars,           ONLY:PartIsImplicit
#endif
#if defined(IMAP) || defined(ROS)
USE MOD_Particle_Vars,           ONLY:PEM,LastPartPos
#endif /*IMPA or ROS*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: N_In,ElemID
INTEGER,INTENT(IN)               :: Mode
INTEGER,INTENT(IN),OPTIONAL      :: PartID
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
REAL                             :: Lag(1:3,0:N_In), F(1:3),Xi_Old(1:3)
INTEGER                          :: NewTonIter,i,j,k
REAL                             :: deltaXi(1:3),deltaXi2
REAL                             :: Jac(1:3,1:3),sdetJac,sJac(1:3,1:3)
REAL                             :: buff,buff2, Norm_F, Norm_F_old,lambda
INTEGER                          :: iArmijo
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

Norm_F=DOT_PRODUCT(F,F)
Norm_F_old=Norm_F
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

  ! Armijo step size control
  deltaXi=MATMUL(sJac,F)
  deltaXi2=DOT_PRODUCT(deltaXi,deltaXi)
  Xi_Old=Xi

  Norm_F_old=Norm_F
  Norm_F=Norm_F*2.
  lambda=1.0
  iArmijo=1
  DO WHILE(Norm_F.GT.Norm_F_old*(1.-0.0001*lambda) .AND.iArmijo.LE.8)

    Xi = Xi_Old - lambda*deltaXI!MATMUL(sJac,F)

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
    lambda=0.2*lambda
    iArmijo=iArmijo+1
    Norm_F=DOT_PRODUCT(F,F)
  END DO ! Armijo iteration

  ! check xi value for plausibility
  IF(ANY(ABS(Xi).GT.1.5)) THEN
    IF(Mode.EQ.1)THEN
      IPWRITE(UNIT_stdOut,*) ' Particle not inside of element, force!!!'
      IPWRITE(UNIT_stdOut,*) ' Newton-Iter', NewtonIter
      IPWRITE(UNIT_stdOut,*) ' xi  ', xi(1:3)
      IPWRITE(UNIT_stdOut,*) ' PartPos', X_in
      IPWRITE(UNIT_stdOut,*) ' ElemID', ElemID+offSetElem
      IF(PRESENT(PartID)) IPWRITE(UNIT_stdOut,*) ' PartID', PartID
#if defined(IMPA)
      IF(PRESENT(PartID)) IPWRITE(UNIT_stdOut,*) ' implicit?', PartisImplicit(PartID)
#endif
#if defined(IMAP) || defined(ROS)
      IF(PRESENT(PartID)) IPWRITE(UNIT_stdOut,*) ' last?', LastPartPos(1:3,PartID)
      IF(PRESENT(PartID)) IPWRITE(UNIT_stdOut,*) ' ElemID-N', PEM%ElementN(PartID)+offSetElem
#endif /*IMPA or ROS*/
        CALL abort(&
  __STAMP__&
  ,'Particle Not inSide of Element, ElemID,',ElemID)
    ELSE
      EXIT
    END IF
  END IF

END DO !newton


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
!> compute determinant of 3x3 matrix
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
!> compute inverse of 3x3 matrix, needs sDet=1/det(Mat)
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


SUBROUTINE GetRefNewtonStartValue(X_in,Xi,ElemID)
!===================================================================================================================================
!> Returns the initial value/ guess for the Newton's algorithm
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Preproc,                 ONLY:PP_N,PP_nElems
USE MOD_Particle_Mesh_Vars,      ONLY:RefMappingGuess,RefMappingEps
USE MOD_Particle_Mesh_Vars,      ONLY:XiEtaZetaBasis,slenXiEtaZetaBasis
USE MOD_Mesh_Vars,               ONLY:Elem_xGP,XCL_NGeo
USE MOD_Interpolation_Vars,      ONLY:xGP
USE MOD_Mesh_Vars,               ONLY:XCL_NGeo,NGeo,XiCL_NGeo,ElemBaryNGeo
USE MOD_Particle_Tracking_vars,  ONLY:DoRefMapping
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: ElemID
REAL,INTENT(IN)                :: X_in(1:3)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(INOUT)             :: Xi(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: Ptild(1:3),XiLinear(1:6)
REAL                          :: Winner_Dist,Dist
REAL                          :: epsOne
INTEGER                       :: iDir
INTEGER                       :: i,j,k
REAL                          :: dX,dY,dZ
INTEGER                       :: RefMappingGuessLoc
!===================================================================================================================================

epsOne=1.0+RefMappingEps
RefMappingGuessLoc=RefMappingGuess
! the location of the Gauss-points within halo elements is not communicated. Instead of looking for the closest Gauss-point, the
! closest CL-point is used
IF(ElemID.GT.PP_nElems)THEN
  IF(DoRefMapping)THEN
    IF(RefMappingGuess.EQ.2) RefMappingGuessLoc=3
  ELSE
    IF(RefMappingGuess.EQ.2) RefMappingGuessLoc=1
  END IF
END IF
SELECT CASE(RefMappingGuessLoc)
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
  ! limit xi to [-1,1]
  IF(MAXVAL(ABS(Xi)).GT.epsOne) Xi=MAX(MIN(1.0d0,Xi),-1.0d0)
CASE(2)
  ! compute distance on Gauss Points
  Winner_Dist=SQRT(DOT_PRODUCT((x_in(:)-Elem_xGP(:,0,0,0,ElemID)),(x_in(:)-Elem_xGP(:,0,0,0,ElemID))))
  Xi(:)=(/xGP(0),xGP(0),xGP(0)/) ! start value
  DO i=0,PP_N; DO j=0,PP_N; DO k=0,PP_N
    dX=ABS(X_in(1) - Elem_xGP(1,i,j,k,ElemID))
    IF(dX.GT.Winner_Dist) CYCLE
    dY=ABS(X_in(2) - Elem_xGP(2,i,j,k,ElemID))
    IF(dY.GT.Winner_Dist) CYCLE
    dZ=ABS(X_in(3) - Elem_xGP(3,i,j,k,ElemID))
    IF(dZ.GT.Winner_Dist) CYCLE
    Dist=SQRT(dX*dX+dY*dY+dZ*dZ)
    IF (Dist.LT.Winner_Dist) THEN
      Winner_Dist=Dist
      Xi(:)=(/xGP(i),xGP(j),xGP(k)/) ! start value
    END IF
  END DO; END DO; END DO
CASE(3)
  ! compute distance on XCL Points
  Winner_Dist=SQRT(DOT_PRODUCT((x_in(:)-XCL_NGeo(:,0,0,0,ElemID)),(x_in(:)-XCL_NGeo(:,0,0,0,ElemID))))
  Xi(:)=(/XiCL_NGeo(0),XiCL_NGeo(0),XiCL_NGeo(0)/) ! start value
  DO i=0,NGeo; DO j=0,NGeo; DO k=0,NGeo
    dX=ABS(X_in(1) - XCL_NGeo(1,i,j,k,ElemID))
    IF(dX.GT.Winner_Dist) CYCLE
    dY=ABS(X_in(2) - XCL_NGeo(2,i,j,k,ElemID))
    IF(dY.GT.Winner_Dist) CYCLE
    dZ=ABS(X_in(3) - XCL_NGeo(3,i,j,k,ElemID))
    IF(dZ.GT.Winner_Dist) CYCLE
    Dist=SQRT(dX*dX+dY*dY+dZ*dZ)
    IF (Dist.LT.Winner_Dist) THEN
      Winner_Dist=Dist
      Xi(:)=(/XiCL_NGeo(i),XiCL_NGeo(j),XiCL_NGeo(k)/) ! start value
    END IF
  END DO; END DO; END DO
CASE(4)
  ! trivial guess
  xi=0.
END SELECT

END SUBROUTINE GetRefNewtonStartValue

END MODULE MOD_Eval_xyz
