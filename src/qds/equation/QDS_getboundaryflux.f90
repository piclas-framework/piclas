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
#if USE_QDS_DG
#include "piclas.h"

MODULE MOD_QDS_GetBoundaryFlux
!===================================================================================================================================
!> Contains the routines to
!> - initialze the QDS method for DG + equation
!===================================================================================================================================
! MODULES
!USE MOD_io_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE GetBoundaryFluxQDS
  MODULE PROCEDURE GetBoundaryFluxQDS
END INTERFACE

PUBLIC::GetBoundaryFluxQDS
!===================================================================================================================================
CONTAINS

SUBROUTINE GetBoundaryFluxQDS(t,tDeriv, Flux, U_Minus, NormVec, TangVec1, TangVec2, BCFace_xGP)
!===================================================================================================================================
! Computes the boundary values for a given Cartesian mesh face (defined by FaceID)
! BCType: 1...periodic, 2...exact BC
! Attention 1: this is only a tensor of local values U_Face and has to be stored into the right U_Left or U_Right in
!              SUBROUTINE CalcSurfInt
! Attention 2: U_FacePeriodic is only needed in the case of periodic boundary conditions
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals,            ONLY:Abort,CROSS
USE MOD_QDS_Equation,       ONLY:QDS_ExactFunc
USE MOD_Mesh_Vars    ,      ONLY:nBCSides,nBCs,BoundaryType
USE MOD_Equation_Vars,      ONLY:nBCByType,BCSideID,nBCByType
USE MOD_QDS_Equation_Vars,  ONLY:U_Face_old,QDSnVar
USE MOD_QDS_Riemann,        ONLY:RiemannQDS
USE MOD_QDS_Equation,       ONLY:QDS_Q2U
USE MOD_QDS_DG_vars,        ONLY:QDSnVar_macro
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                    :: t
INTEGER,INTENT(IN)                   :: tDeriv
REAL,INTENT(IN)                      :: U_Minus(     QDSnVar,0:PP_N,0:PP_N,1:nBCSides)
REAL,INTENT(IN)                      :: NormVec(           3,0:PP_N,0:PP_N,1:nBCSides)
REAL,INTENT(IN),OPTIONAL             :: TangVec1(          3,0:PP_N,0:PP_N,1:nBCSides)
REAL,INTENT(IN),OPTIONAL             :: TangVec2(          3,0:PP_N,0:PP_N,1:nBCSides)
REAL,INTENT(IN)                      :: BCFace_xGP(        3,0:PP_N,0:PP_N,1:nBCSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: Flux(1:QDSnVar,0:PP_N,0:PP_N,1:nBCSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iBC,iSide,p,q,SideID
INTEGER                              :: BCType,BCState,nBCLoc
REAL                                 :: epsBC
REAL                                 :: U_Face_loc(QDSnVar,0:PP_N,0:PP_N),v_nor(3), V_par(3), v_nor_abs
INTEGER                              :: iPart, iPart1, iPart2, iPart3, locPartNum
REAL                                 :: Temp
REAL                                 :: QDSMacroValues(1:QDSnVar_macro)
!===================================================================================================================================

DO iBC=1,nBCs
  IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType(iBC,BC_TYPE)
  BCState=BoundaryType(iBC,BC_STATE)
  nBCLoc =nBCByType(iBC)
  SELECT CASE(BCType)
  CASE(1) !Periodic already filled!

  CASE(2) ! exact BC = Dirichlet BC !!
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          CALL QDS_ExactFunc(BCState,t,tDeriv,BCFace_xGP(:,p,q,SideID), QDSMacroValues(1:QDSnVar_macro))
          IF (QDSMacroValues(1).GT.0.0) THEN
            IF (QDSMacroValues(6).LT.0.0) QDSMacroValues(6) = 0.0
              CALL QDS_Q2U(QDSMacroValues(1:QDSnVar_macro), U_Face_loc(:,p,q))
          ELSE
            U_Face_loc(:,p,q) = 0.0
          END IF
          END DO ! p
      END DO ! q
      CALL RiemannQDS(Flux(1:QDSnVar,:,:,SideID),U_Minus(:,:,:,SideID),U_Face_loc(:,:,:),NormVec(:,:,:,SideID))
   END DO

  CASE(3) ! 1st order absorbing BC
    ! U_Face_loc=0.
    ! CALL RiemannQDS(Flux(1:QDSnVar,:,:,SideID),U_Minus(:,:,:,SideID),U_Face_loc(:,:,:),NormVec(:,:,:,SideID))
    ! Flux = 0.0
    U_Face_loc               = 0.5*(U_Minus(:,:,:,SideID)+U_Face_old(:,:,:,SideID))
    U_Face_old(:,:,:,SideID) = U_Minus(:,:,:,SideID)
    CALL RiemannQDS(Flux(1:QDSnVar,:,:,SideID),U_Minus(:,:,:,SideID),U_Face_loc(:,:,:),NormVec(:,:,:,SideID))

  CASE(4) ! perfectly conducting surface (MunzOmnesSchneider 2000, pp. 97-98)
!    ! Determine the exact BC state
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      ! Determine the exact BC state
      DO q=0,PP_N
        DO p=0,PP_N
          DO iPart = 0, 7
            v_nor_abs = NormVec(1,p,q,SideID)*U_Minus(2+5*iPart,p,q,SideID) + &
                        NormVec(2,p,q,SideID)*U_Minus(3+5*iPart,p,q,SideID) + &
                        NormVec(3,p,q,SideID)*U_Minus(4+5*iPart,p,q,SideID)
            v_nor(1:3) = v_nor_abs*NormVec(1:3,p,q,SideID)
            v_par(1) = U_Minus(2+5*iPart,p,q,SideID) - v_nor(1)
            v_par(2) = U_Minus(3+5*iPart,p,q,SideID) - v_nor(2)
            v_par(3) = U_Minus(4+5*iPart,p,q,SideID) - v_nor(3)
            U_Face_loc(1+5*iPart,p,q) = U_Minus(1+5*iPart,p,q,SideID)
            U_Face_loc(2+5*iPart,p,q) = v_par(1) - v_nor(1)
            U_Face_loc(3+5*iPart,p,q) = v_par(2) - v_nor(2)
            U_Face_loc(4+5*iPart,p,q) = v_par(3) - v_nor(3)
            U_Face_loc(5+5*iPart,p,q) = U_Minus(5+5*iPart,p,q,SideID)
          END DO
        END DO ! p
      END DO ! q
      ! Dirichlet means that we use the gradients from inside the grid cell
      CALL RiemannQDS(Flux(1:QDSnVar,:,:,SideID),U_Minus(:,:,:,SideID),U_Face_loc(:,:,:),NormVec(:,:,:,SideID))
      Flux = 0.0
    END DO

  CASE(10) ! symmetry BC (perfect MAGNETIC conductor, PMC)

  CASE(20) ! exact BC = Dirichlet BC !!
!    ! SPECIAL BC: BCState uses readin state
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
!      ! Dirichlet means that we use the gradients from inside the grid cell
!      CALL RiemannQDS(Flux(1:QDSnVar,:,:,SideID),U_Minus( :,:,:,SideID),BCData(:,:,:,SideID),NormVec(:,:,:,SideID))
!    END DO

      DO q=0,PP_N; DO p=0,PP_N
        locPartNum = 0
        DO iPart1=1,2; DO iPart2=1,2; DO iPart3=1,2
          Temp = NormVec(1,p,q,SideID)*U_Minus(2+locPartNum*5,p,q,SideID) &
                +NormVec(2,p,q,SideID)*U_Minus(3+locPartNum*5,p,q,SideID) &
                +NormVec(3,p,q,SideID)*U_Minus(4+locPartNum*5,p,q,SideID)
          Flux(1+locPartNum*5,p,q,SideID) =  U_Minus(1+locPartNum*5,p,q,SideID)*Temp
          Flux(2+locPartNum*5,p,q,SideID) =  U_Minus(2+locPartNum*5,p,q,SideID)*Temp
          Flux(3+locPartNum*5,p,q,SideID) =  U_Minus(3+locPartNum*5,p,q,SideID)*Temp
          Flux(4+locPartNum*5,p,q,SideID) =  U_Minus(4+locPartNum*5,p,q,SideID)*Temp
          Flux(5+locPartNum*5,p,q,SideID) =  U_Minus(5+locPartNum*5,p,q,SideID)*Temp


    !      U_Face_loc(1+locPartNum*5,p,q) =  U_Minus(1+locPartNum*5,p,q,SideID)
    !      U_Face_loc(2+locPartNum*5,p,q) =  U_Minus(2+locPartNum*5,p,q,SideID)
    !      U_Face_loc(3+locPartNum*5,p,q) =  U_Minus(3+locPartNum*5,p,q,SideID)
    !      U_Face_loc(4+locPartNum*5,p,q) =  U_Minus(4+locPartNum*5,p,q,SideID)
    !      U_Face_loc(5+locPartNum*5,p,q) =  U_Minus(5+locPartNum*5,p,q,SideID)
          locPartNum = locPartNum + 1
        END DO; END DO; END DO
      END DO; END DO
    END DO


  CASE DEFAULT ! unknown BCType
    CALL abort(&
__STAMP__&
        ,'no BC defined in maxwell/getboundaryflux.f90!')
  END SELECT ! BCType
END DO ! iBC=1,nBC

IF(1.EQ.2)THEN
  epsBC=TangVec1(1,1,1,1)
  epsBC=TangVec2(1,1,1,1)
END IF

END SUBROUTINE GetBoundaryFluxQDS

END MODULE MOD_QDS_GetBoundaryFlux
#endif /*USE_QDS_DG*/
