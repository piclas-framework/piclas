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

INTERFACE InitMortar_HDG
  MODULE PROCEDURE InitMortar_HDG
END INTERFACE

! reshape vector nGP_face to 0:PP_N,0:PP_N
!INTERFACE BigToSmallMortar_HDG
!  MODULE PROCEDURE BigToSmallMortar_HDG
!END INTERFACE

! reshape vector nGP_face to 0:PP_N,0:PP_N
!INTERFACE SmallToBigMortar_HDG
!  MODULE PROCEDURE SmallToBigMortar_HDG
!END INTERFACE

INTERFACE SmallToBigMortarPrecond_HDG
  MODULE PROCEDURE SmallToBigMortarPrecond_HDG
END INTERFACE

PUBLIC::InitMortar_HDG
PUBLIC::BigToSmallMortar_HDG
PUBLIC::SmallToBigMortar_HDG
PUBLIC::SmallToBigMortarPrecond_HDG

!===================================================================================================================================

CONTAINS


SUBROUTINE InitMortar_HDG()
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Mortar_Vars     ,ONLY: M_0_1,M_0_2
USE MOD_HDG_Vars        ,ONLY: MaskedSide,SmallMortarInfo,IntMatMortar,nGP_Face
USE MOD_Mesh_Vars       ,ONLY: nSides,MortarType,MortarInfo
#if USE_PETSC
USE MOD_HDG_Vars        ,ONLY: SmallMortarType
#if USE_MPI
USE MOD_MPI             ,ONLY: StartReceiveMPIDataInt,StartSendMPIDataInt,FinishExchangeMPIData
USE MOD_MPI_Vars
#endif /*USE_MPI*/
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER     :: SideID,iSide,mtype,iMortar,nMortars
INTEGER     :: i,j,p,q,iGP_face,jGP_face
REAL        :: dkron(0:PP_N,0:PP_N)
!===================================================================================================================================
  ALLOCATE(SmallMortarInfo(1:nSides))
  SmallMortarInfo=0
#if USE_PETSC
  ALLOCATE(SmallMortarType(2,1:nSides))
  SmallMortarType=-1
#endif

  !=-1: small mortar neighbor (slave)
  !=0: not a small mortar
  != 1: small mortar at big side

  DO SideID=1,nSides
    mtype=MortarType(1,SideID)
    IF(mtype.EQ.-1) THEN
      SmallMortarInfo(SideID)= 0
    ELSEIF(mtype.EQ.0)THEN
      SmallMortarInfo(SideID)= 1
    ELSEIF(mtype.EQ.-10)THEN
      SmallMortarInfo(SideID)=-1
    ELSEIF(mtype.GT.0) THEN
      !is a big mortar side:
      nMortars=MERGE(4,2,mtype.EQ.1)
      iSide=MortarType(2,SideID)  !index of Big Side in MortarInfo
      DO iMortar=1,nMortars
        SmallMortarInfo(  MortarInfo(MI_SIDEID,iMortar,iSide) ) = 1  !small sideID
#if USE_PETSC
        SmallMortarType(1,MortarInfo(MI_SIDEID,iMortar,iSide)) = mtype
        SmallMortarType(2,MortarInfo(MI_SIDEID,iMortar,iSide)) = iMortar
#endif
        IF(MortarType(1,MortarInfo(MI_SIDEID,iMortar,iSide)).NE.0) THEN
          IPWRITE(*,*)SideID,MortarType(1,MortarInfo(MI_SIDEID,iMortar,iSide))
          CALL abort(__STAMP__&
                    ,'InitMortar_HDG: check failed')
        END IF
      END DO
    ELSE
      CALL abort(__STAMP__&
                ,'InitMortar_HDG: this case should not appear!!')
    END IF
    ! TODO new value for mortars
    IF(SmallMortarInfo(SideID).NE.0) MaskedSide(SideID)=1!.TRUE.
  END DO !SideID=1,nSides

#if USE_PETSC
#if USE_MPI
  ! Send Mortar data to small sides on other processor (needed for mtype=-10 sides)
  CALL StartReceiveMPIDataInt(2,SmallMortarType,1,nSides, RecRequest_U,SendID=1) ! Receive YOUR
  CALL StartSendMPIDataInt(   2,SmallMortarType,1,nSides,SendRequest_U,SendID=1) ! Send MINE
  CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=1)
#endif /*USE_MPI*/
#endif

  ! not efficient but simpler: build full interpolation matrices for 3 mortartypes
  dkron=0.
  DO i=0,PP_N
    dkron(i,i)=1. !kronecker
  END DO

  ALLOCATE(IntMatMortar(nGP_face,nGP_Face,4,3)) ! 4-1 , 2-1 in eta, 2-1 in xi
  IntMatMortar=0.
  jGP_face=0
  DO j=0,PP_N; DO i=0,PP_N
    jGP_Face=jGP_Face+1
    iGP_Face=0
    DO q=0,PP_N; DO p=0,PP_N
      iGP_Face=iGP_Face+1
      !type 1: update
      ! U_tmp(:,p,q,1)=U_tmp(:,p,q,1)+M_0_1(i,p)*M_0_1(j,q)*lambda_in(:,i,j,MortarSideID)
      ! U_tmp(:,p,q,2)=U_tmp(:,p,q,2)+M_0_2(i,p)*M_0_1(j,q)*lambda_in(:,i,j,MortarSideID)
      ! U_tmp(:,p,q,3)=U_tmp(:,p,q,3)+M_0_1(i,p)*M_0_2(j,q)*lambda_in(:,i,j,MortarSideID)
      ! U_tmp(:,p,q,4)=U_tmp(:,p,q,4)+M_0_2(i,p)*M_0_2(j,q)*lambda_in(:,i,j,MortarSideID)
      IntMatMortar(iGP_Face,jGP_Face,1,1) =IntMatMortar(iGP_Face,jGP_Face,1,1) + M_0_1(i,p)*M_0_1(j,q)
      IntMatMortar(iGP_Face,jGP_Face,2,1) =IntMatMortar(iGP_Face,jGP_Face,2,1) + M_0_2(i,p)*M_0_1(j,q)
      IntMatMortar(iGP_Face,jGP_Face,3,1) =IntMatMortar(iGP_Face,jGP_Face,3,1) + M_0_1(i,p)*M_0_2(j,q)
      IntMatMortar(iGP_Face,jGP_Face,4,1) =IntMatMortar(iGP_Face,jGP_Face,4,1) + M_0_2(i,p)*M_0_2(j,q)

      !type 2: update
      ! U_tmp(:,p,q,1)=U_tmp(:,p,q,1)+dkron(i,p)*M_0_1(j,q)*lambda_in(:,i,j,MortarSideID)
      ! U_tmp(:,p,q,2)=U_tmp(:,p,q,2)+dkron(i,p)*M_0_2(j,q)*lambda_in(:,i,j,MortarSideID)
      IntMatMortar(iGP_Face,jGP_Face,1,2) =IntMatMortar(iGP_Face,jGP_Face,1,2) + dkron(i,p)*M_0_1(j,q)
      IntMatMortar(iGP_Face,jGP_Face,2,2) =IntMatMortar(iGP_Face,jGP_Face,2,2) + dkron(i,p)*M_0_2(j,q)

      !type 3: update
      ! U_tmp(:,p,q,1)=U_tmp(:,p,q,1)+M_0_1(i,p)*dkron(j,q)*lambda_in(:,i,j,MortarSideID)
      ! U_tmp(:,p,q,2)=U_tmp(:,p,q,2)+M_0_2(i,p)*dkron(j,q)*lambda_in(:,i,j,MortarSideID)
      IntMatMortar(iGP_Face,jGP_Face,1,3) =IntMatMortar(iGP_Face,jGP_Face,1,3) + M_0_1(i,p)*dkron(j,q)
      IntMatMortar(iGP_Face,jGP_Face,2,3) =IntMatMortar(iGP_Face,jGP_Face,2,3) + M_0_2(i,p)*dkron(j,q)
    END DO; END DO !p,q (iGP_face)
  END DO; END DO !i,j (jGP_face)
END SUBROUTINE InitMortar_HDG


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
        CALL abort(__STAMP__&
                  ,'BigToSmallMortar_HDG: small sides should not be slave!!!!')
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
      CALL abort(__STAMP__&
                ,'SmallToBigMortar_HDG small side should not be slave!!!')
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


SUBROUTINE SmallToBigMortarPrecond_HDG(whichPrecond)
!===================================================================================================================================
!> Mortar matvec action: matrix of small Mortar side (j) SMat_j is multiplied with a small lambda lambda_small_j,
!> which was computed by interpolation with lambda_small_j = IMat_j lambda_big
!> mv_small_j = SMat_j * lambda_small_j = SMat_j *IMat_j * lambda_big
!> then the result is 'interpolated back' and added to big:
!>  mv_big +=  (Imat_j)^T * mv_small_j = (Imat_j)^T SMat_j Imat_j lambda_big
!>
!> Here the "naked" small Side Matrix (1:nGP_face,1:nGP_face) must be multiplied with the Interpolation matrix from the left
!> and its transpose from the right.
!>
!> not optimized implementation, since buildPrecond is not so oftern called!
!>
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,   ONLY: MortarType,MortarInfo
USE MOD_Mesh_Vars,   ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_HDG_Vars,    ONLY: nGP_Face,Precond,InvPrecondDiag,IntMatMortar
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: whichPrecond
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER     :: i,k,iSide,iMortar,nMortars,mtype,MortarSideID,SmallSideID
!===================================================================================================================================
  IF(whichPrecond.EQ.0) RETURN


  DO MortarSideID=firstMortarInnerSide,lastMortarInnerSide  !Big SideID

    nMortars=MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
    mtype=MortarType(1,MortarSideID)
    iSide=MortarType(2,MortarSideID)  !index of Big Side in MortarInfo
    DO iMortar=1,nMortars
      SmallSideID = MortarInfo(MI_SIDEID,iMortar,iSide)
      SELECT CASE(whichPrecond)
      CASE(1) !side-block matrix
        Precond(:,:,MortarSideID) = Precond(:,:,MortarSideID)                           &
                                    + MATMUL(TRANSPOSE(IntMatMortar(:,:,iMortar,mtype)),      &
                                             MATMUL(Precond(:,:,SmallSideID),IntMatMortar(:,:,iMortar,mtype)))
        Precond(:,:,SmallSideID)=0.
      CASE(2) !only diagonal part of the  matrix , M_ij= I_ki D_kk I_kj, only i=j
        DO i=1,nGP_Face
          DO k=1,nGP_Face
            InvPrecondDiag(i,MortarSideID)= InvPrecondDiag(i,MortarSideID)                  &
                                + IntMatMortar(k,i,iMortar,mtype)*InvPrecondDiag(k,SmallSideID)*IntMatMortar(k,i,iMortar,mtype)
          END DO !k
        END DO !i
        InvPrecondDiag(:,SmallSideID)=0.
      END SELECT
    END DO !iMortar
  END DO !MortarSideID

END SUBROUTINE SmallToBigMortarPrecond_HDG

END MODULE MOD_FillMortar_HDG
