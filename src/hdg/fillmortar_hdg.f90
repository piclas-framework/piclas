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

MODULE MOD_FillMortar_HDG
!===================================================================================================================================
! Add comments please!
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
! reshape vector nGP_face to 0:PP_N,0:PP_N
!INTERFACE BigToSmallMortar_HDG
!  MODULE PROCEDURE BigToSmallMortar_HDG
!END INTERFACE

! reshape vector nGP_face to 0:PP_N,0:PP_N
!INTERFACE SmallToBigMortar_HDG
!  MODULE PROCEDURE SmallToBigMortar_HDG
!END INTERFACE

PUBLIC::BigToSmallMortar_HDG
PUBLIC::SmallToBigMortar_HDG

!===================================================================================================================================

CONTAINS

SUBROUTINE BigToSmallMortar_HDG(nVar_in,lambda_in)
!===================================================================================================================================
!> fills small non-conforming sides with data for master side with data from the corresponding large side, using 1D interpolation 
!> operators M_0_1,M_0_2
!>
!> REMARK: NO doMPISides, because nMortarMPIsides=0 has to be guaranteed!
!
!     Type 1               Type 2              Type3
!      eta                  eta                 eta
!       ^                    ^                   ^
!       |                    |                   |
!   +---+---+            +---+---+           +---+---+
!   | 3 | 4 |            |   2   |           |   |   |
!   +---+---+ --->  xi   +---+---+ --->  xi  + 1 + 2 + --->  xi
!   | 1 | 2 |            |   1   |           |   |   |
!   +---+---+            +---+---+           +---+---+
!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mortar_Vars, ONLY: M_0_1,M_0_2
USE MOD_Mesh_Vars,   ONLY: MortarType,MortarInfo
USE MOD_Mesh_Vars,   ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_Mesh_Vars,   ONLY: nSides 
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(IN) :: nVar_in
REAL,INTENT(INOUT) :: lambda_in(1:nVar_in,0:PP_N,0:PP_N,1:nSides) !< (INOUT) can be U or Grad_Ux/y/z_master
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER      :: p,q,l
INTEGER      :: iMortar,nMortars
INTEGER      :: MortarSideID,SideID,locSide,flip
REAL         :: U_tmp( 1:nVar_in,0:PP_N,0:PP_N,1:4)
REAL         :: U_tmp2(1:nVar_in,0:PP_N,0:PP_N,1:2)
!===================================================================================================================================


DO MortarSideID=firstMortarInnerSide,lastMortarInnerSide
  !
  SELECT CASE(MortarType(1,MortarSideID))
  CASE(1) !1->4
    !first in eta
    ! The following q- and l-loop are two MATMULs: (ATTENTION M_0_1 and M_0_2 are already transposed in mortar.f90)
    !    U_tmp2(:,p,:,1)  =  M_0_1 * lambda_in(p,:,MortarSideID)
    !    U_tmp2(:,p,:,2)  =  M_0_2 * lambda_in(p,:,MortarSideID)
    DO q=0,PP_N
      DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction 
        U_tmp2(:,p,q,1)=                  M_0_1(0,q)*lambda_in(:,p,0,MortarSideID)
        U_tmp2(:,p,q,2)=                  M_0_2(0,q)*lambda_in(:,p,0,MortarSideID)
        DO l=1,PP_N
          U_tmp2(:,p,q,1)=U_tmp2(:,p,q,1)+M_0_1(l,q)*lambda_in(:,p,l,MortarSideID)
          U_tmp2(:,p,q,2)=U_tmp2(:,p,q,2)+M_0_2(l,q)*lambda_in(:,p,l,MortarSideID)
        END DO
      END DO
    END DO
    ! then in xi
    DO q=0,PP_N ! for every eta-layer perform Mortar operation in xi-direction 
      ! The following p- and l-loop are four MATMULs: (ATTENTION M_0_1 and M_0_2 are already transposed in mortar.f90)
      !    U_tmp(:,:,q,1)  =  M_0_1 * U_tmp2(:,:,q,1)
      !    U_tmp(:,:,q,2)  =  M_0_2 * U_tmp2(:,:,q,1)
      !    U_tmp(:,:,q,3)  =  M_0_1 * U_tmp2(:,:,q,2)
      !    U_tmp(:,:,q,4)  =  M_0_2 * U_tmp2(:,:,q,2)
      DO p=0,PP_N
        U_tmp(:,p,q,1)=                 M_0_1(0,p)*U_tmp2(:,0,q,1)
        U_tmp(:,p,q,2)=                 M_0_2(0,p)*U_tmp2(:,0,q,1)
        U_tmp(:,p,q,3)=                 M_0_1(0,p)*U_tmp2(:,0,q,2)
        U_tmp(:,p,q,4)=                 M_0_2(0,p)*U_tmp2(:,0,q,2)
        DO l=1,PP_N
          U_tmp(:,p,q,1)=U_tmp(:,p,q,1)+M_0_1(l,p)*U_tmp2(:,l,q,1)
          U_tmp(:,p,q,2)=U_tmp(:,p,q,2)+M_0_2(l,p)*U_tmp2(:,l,q,1)
          U_tmp(:,p,q,3)=U_tmp(:,p,q,3)+M_0_1(l,p)*U_tmp2(:,l,q,2)
          U_tmp(:,p,q,4)=U_tmp(:,p,q,4)+M_0_2(l,p)*U_tmp2(:,l,q,2)
        END DO !l=1,PP_N
      END DO
    END DO 

  CASE(2) !1->2 in eta
    ! The following q- and l-loop are two MATMULs: (ATTENTION M_0_1 and M_0_2 are already transposed in mortar.f90)
    !    U_tmp(:,p,:,1)  =  M_0_1 * lambda_in(:,p,:,MortarSideID)
    !    U_tmp(:,p,:,2)  =  M_0_2 * lambda_in(:,p,:,MortarSideID)
    DO q=0,PP_N
      DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction 
        U_tmp(:,p,q,1)=                 M_0_1(0,q)*lambda_in(:,p,0,MortarSideID)
        U_tmp(:,p,q,2)=                 M_0_2(0,q)*lambda_in(:,p,0,MortarSideID)
        DO l=1,PP_N
          U_tmp(:,p,q,1)=U_tmp(:,p,q,1)+M_0_1(l,q)*lambda_in(:,p,l,MortarSideID)
          U_tmp(:,p,q,2)=U_tmp(:,p,q,2)+M_0_2(l,q)*lambda_in(:,p,l,MortarSideID)
        END DO
      END DO
    END DO

  CASE(3) !1->2 in xi
    DO q=0,PP_N ! for every eta-layer perform Mortar operation in xi-direction
      ! The following p- and l-loop are two MATMULs: (ATTENTION M_0_1 and M_0_2 are already transposed in mortar.f90)
      !    U_tmp(:,:,q,1)  =  M_0_1 * lambda_in(:,:,q,MortarSideID)
      !    U_tmp(:,:,q,2)  =  M_0_2 * lambda_in(:,:,q,MortarSideID)
      DO p=0,PP_N
        U_tmp(:,p,q,1)=                 M_0_1(0,p)*lambda_in(:,0,q,MortarSideID)
        U_tmp(:,p,q,2)=                 M_0_2(0,p)*lambda_in(:,0,q,MortarSideID)
        DO l=1,PP_N
          U_tmp(:,p,q,1)=U_tmp(:,p,q,1)+M_0_1(l,p)*lambda_in(:,l,q,MortarSideID)
          U_tmp(:,p,q,2)=U_tmp(:,p,q,2)+M_0_2(l,p)*lambda_in(:,l,q,MortarSideID)
        END DO
      END DO
    END DO
  END SELECT ! mortarType(SideID)

  nMortars=MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
  locSide=MortarType(2,MortarSideID)
  DO iMortar=1,nMortars
    SideID= MortarInfo(MI_SIDEID,iMortar,locSide) !small SideID
    flip  = MortarInfo(MI_FLIP,iMortar,locSide)
    SELECT CASE(flip)
      CASE(0) ! master side
        lambda_in(:,:,:,SideID)=U_tmp(:,:,:,iMortar)
      CASE(1:4) ! slave side
        STOP 'BigToSmallMortar_HDG: small sides should not be slave!!!!'
!        DO q=0,PP_N; DO p=0,PP_N
!          U_in_slave(p,q,SideID)=U_tmp(FS2M(1,p,q,flip), &
!                                       FS2M(2,p,q,flip),iMortar)
!        END DO; END DO ! q, p
    END SELECT !flip(iMortar)
  END DO !iMortar
END DO !MortarSideID
END SUBROUTINE BigToSmallMortar_HDG

SUBROUTINE SmallToBigMortar_HDG(nVar_in,mv_in)
!===================================================================================================================================
! fills master side from small non-conforming sides, Using 1D projection operators M_1_0,M_2_0
!
!> This routine takes the contribution from the matrix-vector product in HDG of the small element sides  
!> and adds to the big sides, using the transpose of the BigToSmall interpolation operator
!> and also sets then the small mortar contribution to zero!
!>
!> REMARK: NO doMPISides, because nMortarMPIsides=0 has to be guaranteed!
!
!     Type 1               Type 2              Type3
!      eta                  eta                 eta
!       ^                    ^                   ^
!       |                    |                   |
!   +---+---+            +---+---+           +---+---+
!   | 3 | 4 |            |   2   |           |   |   |
!   +---+---+ --->  xi   +---+---+ --->  xi  + 1 + 2 + --->  xi
!   | 1 | 2 |            |   1   |           |   |   |
!   +---+---+            +---+---+           +---+---+
!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mortar_Vars, ONLY: M_0_1,M_0_2
USE MOD_Mesh_Vars,   ONLY: MortarType,MortarInfo,nSides
USE MOD_Mesh_Vars,   ONLY: firstMortarInnerSide,lastMortarInnerSide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: nVar_in
REAL,INTENT(INOUT) :: mv_in(nVar_in,0:PP_N,0:PP_N,1:nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER  :: p,q,l
INTEGER  :: iMortar,nMortars
INTEGER  :: MortarSideID,SideID,iSide,flip
REAL     :: mv_tmp( nVar_in,0:PP_N,0:PP_N,1:4)
REAL     :: mv_tmp2(nVar_in,0:PP_N,0:PP_N,1:2)
REAL     :: M_1_0(0:PP_N,0:PP_N),M_2_0(0:PP_N,0:PP_N) 
!===================================================================================================================================
M_1_0 = TRANSPOSE(M_0_1)
M_2_0 = TRANSPOSE(M_0_2)

DO MortarSideID=firstMortarInnerSide,lastMortarInnerSide  !Big SideID

  nMortars=MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
  iSide=MortarType(2,MortarSideID)  !index of Big Side in MortarInfo
  DO iMortar=1,nMortars
    SideID = MortarInfo(MI_SIDEID,iMortar,iSide) !small sideID
    flip   = MortarInfo(MI_FLIP,iMortar,iSide)
    SELECT CASE(flip)
    CASE(0) ! master side
      mv_tmp(:,:,:,iMortar)=mv_in(:,:,:,SideID)
      !!!!! SET small mortar side contribution to zero here!!!
      mv_in(:,:,:,SideID)  = 0.
    CASE(1:4) ! slave sides (should only occur for MPI)
      STOP 'SmallToBigMortar_HDG small side should not be slave!!!'
!      DO q=0,PP_N; DO p=0,PP_N
!        mv_tmp(FS2M(1,p,q,flip),FS2M(2,p,q,flip),iMortar)=-mv_Slave(p,q,SideID)
!      END DO; END DO
    END SELECT
  END DO

  SELECT CASE(MortarType(1,MortarSideID))
  CASE(1) !1->4
    ! first in xi
    DO q=0,PP_N ! for every eta-layer perform Mortar operation in xi-direction 
      ! The following p- and l-loop are four MATMULs: (ATTENTION M_1_0 and M_2_0 are already transposed in mortar.f90)
      !    mv_tmp2(:,:,q,1)  =  M_1_0 * mv_tmp(:,:,q,1) + M_2_0 * mv_tmp(:,:,q,2)
      !    mv_tmp2(:,:,q,2)  =  M_1_0 * mv_tmp(:,:,q,1) + M_2_0 * mv_tmp(:,:,q,2)
      DO p=0,PP_N
        mv_tmp2(:,p,q,1)=                     M_1_0(0,p)*mv_tmp(:,0,q,1)+M_2_0(0,p)*mv_tmp(:,0,q,2)
        mv_tmp2(:,p,q,2)=                     M_1_0(0,p)*mv_tmp(:,0,q,3)+M_2_0(0,p)*mv_tmp(:,0,q,4)
        DO l=1,PP_N
          mv_tmp2(:,p,q,1)=mv_tmp2(:,p,q,1) + M_1_0(l,p)*mv_tmp(:,l,q,1)+M_2_0(l,p)*mv_tmp(:,l,q,2)
          mv_tmp2(:,p,q,2)=mv_tmp2(:,p,q,2) + M_1_0(l,p)*mv_tmp(:,l,q,3)+M_2_0(l,p)*mv_tmp(:,l,q,4)
        END DO
      END DO
    END DO
    !then in eta
    ! The following q- and l-loop are two MATMULs: (ATTENTION M_1_0 and M_2_0 are already transposed in mortar.f90)
    !    mv_in(p,:,MortarSideID)  +=  M_1_0 * mv_tmp2(p,:,1) + M_2_0 * mv_tmp2(p,:,2)
    DO q=0,PP_N
      DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction 
        DO l=0,PP_N
          mv_in(:,p,q,MortarSideID)=mv_in(:,p,q,MortarSideID) + M_1_0(l,q)*mv_tmp2(:,p,l,1)+M_2_0(l,q)*mv_tmp2(:,p,l,2)
        END DO
      END DO
    END DO

  CASE(2) !1->2 in eta
    ! TODO why not q-loop first?
    DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction 
      ! The following q- and l-loop are two MATMULs: (ATTENTION M_1_0 and M_2_0 are already transposed in mortar.f90)
      !    mv_in(p,:,MortarSideID)  +=  M_1_0 * mv_tmp(p,:,1) + M_2_0 * mv_tmp(p,:,2)
      DO q=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction 
        DO l=0,PP_N
          mv_in(:,p,q,MortarSideID)=mv_in(:,p,q,MortarSideID) + M_1_0(l,q)*mv_tmp(:,p,l,1)+M_2_0(l,q)*mv_tmp(:,p,l,2)
        END DO
      END DO
    END DO

  CASE(3) !1->2 in xi
    DO q=0,PP_N ! for every eta-layer perform Mortar operation in xi-direction 
      ! The following p- and l-loop are two MATMULs: (ATTENTION M_1_0 and M_2_0 are already transposed in mortar.f90)
      !    mv_in(:,q,MortarSideID) + =   M_1_0 * mv_tmp(:,q,1) + M_2_0 * mv_tmp(:,q,2)
      DO p=0,PP_N
        DO l=0,PP_N
          mv_in(:,p,q,MortarSideID)=mv_in(:,p,q,MortarSideID) + M_1_0(l,p)*mv_tmp(:,l,q,1)+M_2_0(l,p)*mv_tmp(:,l,q,2)
        END DO
      END DO
    END DO

  END SELECT ! mortarType(MortarSideID)
END DO !MortarSideID
END SUBROUTINE SmallToBigMortar_HDG

END MODULE MOD_FillMortar_HDG
