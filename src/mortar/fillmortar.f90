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

MODULE MOD_FillMortar
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
INTERFACE U_Mortar
  MODULE PROCEDURE U_Mortar
END INTERFACE

#if !(USE_HDG)
INTERFACE Flux_Mortar
  MODULE PROCEDURE Flux_Mortar
END INTERFACE
PUBLIC:: Flux_Mortar
#endif /*USE_HDG*/

PUBLIC::U_Mortar

#if PP_Lifting==2 /*BR2*/
INTERFACE Flux_Mortar_BR2
  MODULE PROCEDURE Flux_Mortar_BR2
END INTERFACE

PUBLIC::Flux_Mortar_BR2
#endif /*PP_Lifting==2 */

!===================================================================================================================================

CONTAINS

SUBROUTINE U_Mortar(U_in_master,U_in_slave,doMPISides)
!===================================================================================================================================
!> fills small non-conforming sides with data for master side with data from the corresponding large side, using 1D interpolation
!> operators M_0_1,M_0_2
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
USE MOD_Preproc
USE MOD_Mortar_Vars, ONLY: M_0_1,M_0_2
USE MOD_Mesh_Vars,   ONLY: MortarType,MortarInfo
USE MOD_Mesh_Vars,   ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_Mesh_Vars,   ONLY: firstMortarMPISide,lastMortarMPISide
USE MOD_Mesh_Vars,   ONLY: FS2M,nSides
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: U_in_master(1:PP_nVar,0:PP_N,0:PP_N,1:nSides) !< (INOUT) can be U or Grad_Ux/y/z_master
REAL,INTENT(INOUT) :: U_in_slave( 1:PP_nVar,0:PP_N,0:PP_N,1:nSides) !< (INOUT) can be U or Grad_Ux/y/z_master
LOGICAL,INTENT(IN) :: doMPISides                                 !< flag whether MPI sides are processed
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER      :: p,q,l
INTEGER      :: iMortar,nMortars
INTEGER      :: firstMortarSideID,lastMortarSideID
INTEGER      :: MortarSideID,SideID,locSide,flip
REAL     :: U_tmp( PP_nVar,0:PP_N,0:PP_N,1:4)
REAL     :: U_tmp2(PP_nVar,0:PP_N,0:PP_N,1:2)
REAL,POINTER :: M1(:,:),M2(:,:)
!===================================================================================================================================

firstMortarSideID = MERGE(firstMortarMPISide,firstMortarInnerSide,doMPISides)
 lastMortarSideID = MERGE( lastMortarMPISide, lastMortarInnerSide,doMPISides)

M1=>M_0_1; M2=>M_0_2

DO MortarSideID=firstMortarSideID,lastMortarSideID
  !
  SELECT CASE(MortarType(1,MortarSideID))
  CASE(1) !1->4
    !first in eta
    ! The following q- and l-loop are two MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
    !    U_tmp2(iVar,p,:,1)  =  M1 * U_in_master(iVar,p,:,MortarSideID)
    !    U_tmp2(iVar,p,:,2)  =  M2 * U_in_master(iVar,p,:,MortarSideID)
    DO q=0,PP_N
      DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction
        U_tmp2(:,p,q,1)=                  M1(0,q)*U_in_master(:,p,0,MortarSideID)
        U_tmp2(:,p,q,2)=                  M2(0,q)*U_in_master(:,p,0,MortarSideID)
        DO l=1,PP_N
          U_tmp2(:,p,q,1)=U_tmp2(:,p,q,1)+M1(l,q)*U_in_master(:,p,l,MortarSideID)
          U_tmp2(:,p,q,2)=U_tmp2(:,p,q,2)+M2(l,q)*U_in_master(:,p,l,MortarSideID)
        END DO
      END DO
    END DO
    ! then in xi
    DO q=0,PP_N ! for every eta-layer perform Mortar operation in xi-direction
      ! The following p- and l-loop are four MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
      !    U_tmp(iVar,:,q,1)  =  M1 * U_tmp2(iVar,:,q,1)
      !    U_tmp(iVar,:,q,2)  =  M2 * U_tmp2(iVar,:,q,1)
      !    U_tmp(iVar,:,q,3)  =  M1 * U_tmp2(iVar,:,q,2)
      !    U_tmp(iVar,:,q,4)  =  M2 * U_tmp2(iVar,:,q,2)
      DO p=0,PP_N
        U_tmp(:,p,q,1)=                 M1(0,p)*U_tmp2(:,0,q,1)
        U_tmp(:,p,q,2)=                 M2(0,p)*U_tmp2(:,0,q,1)
        U_tmp(:,p,q,3)=                 M1(0,p)*U_tmp2(:,0,q,2)
        U_tmp(:,p,q,4)=                 M2(0,p)*U_tmp2(:,0,q,2)
        DO l=1,PP_N
          U_tmp(:,p,q,1)=U_tmp(:,p,q,1)+M1(l,p)*U_tmp2(:,l,q,1)
          U_tmp(:,p,q,2)=U_tmp(:,p,q,2)+M2(l,p)*U_tmp2(:,l,q,1)
          U_tmp(:,p,q,3)=U_tmp(:,p,q,3)+M1(l,p)*U_tmp2(:,l,q,2)
          U_tmp(:,p,q,4)=U_tmp(:,p,q,4)+M2(l,p)*U_tmp2(:,l,q,2)
        END DO !l=1,PP_N
      END DO
    END DO

  CASE(2) !1->2 in eta
    ! The following q- and l-loop are two MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
    !    U_tmp(iVar,p,:,1)  =  M1 * U_in_master(iVar,p,:,MortarSideID)
    !    U_tmp(iVar,p,:,2)  =  M2 * U_in_master(iVar,p,:,MortarSideID)
    DO q=0,PP_N
      DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction
        U_tmp(:,p,q,1)=                 M1(0,q)*U_in_master(:,p,0,MortarSideID)
        U_tmp(:,p,q,2)=                 M2(0,q)*U_in_master(:,p,0,MortarSideID)
        DO l=1,PP_N
          U_tmp(:,p,q,1)=U_tmp(:,p,q,1)+M1(l,q)*U_in_master(:,p,l,MortarSideID)
          U_tmp(:,p,q,2)=U_tmp(:,p,q,2)+M2(l,q)*U_in_master(:,p,l,MortarSideID)
        END DO
      END DO
    END DO

  CASE(3) !1->2 in xi
    DO q=0,PP_N ! for every eta-layer perform Mortar operation in xi-direction
      ! The following p- and l-loop are two MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
      !    U_tmp(iVar,:,q,1)  =  M1 * U_in_master(iVar,:,q,MortarSideID)
      !    U_tmp(iVar,:,q,2)  =  M2 * U_in_master(iVar,:,q,MortarSideID)
      DO p=0,PP_N
        U_tmp(:,p,q,1)=                 M1(0,p)*U_in_master(:,0,q,MortarSideID)
        U_tmp(:,p,q,2)=                 M2(0,p)*U_in_master(:,0,q,MortarSideID)
        DO l=1,PP_N
          U_tmp(:,p,q,1)=U_tmp(:,p,q,1)+M1(l,p)*U_in_master(:,l,q,MortarSideID)
          U_tmp(:,p,q,2)=U_tmp(:,p,q,2)+M2(l,p)*U_in_master(:,l,q,MortarSideID)
        END DO
      END DO
    END DO
  END SELECT ! mortarType(SideID)

  nMortars=MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
  locSide=MortarType(2,MortarSideID)
  DO iMortar=1,nMortars
    SideID= MortarInfo(MI_SIDEID,iMortar,locSide)
    flip  = MortarInfo(MI_FLIP,iMortar,locSide)
    SELECT CASE(flip)
      CASE(0) ! master side
        U_in_master(:,:,:,SideID)=U_tmp(:,:,:,iMortar)
      CASE(1:4) ! slave side
        DO q=0,PP_N; DO p=0,PP_N
          U_in_slave(:,p,q,SideID)=U_tmp(:,FS2M(1,p,q,flip), &
                                           FS2M(2,p,q,flip),iMortar)
        END DO; END DO ! q, p
    END SELECT !flip(iMortar)
  END DO !iMortar
END DO !MortarSideID
END SUBROUTINE U_Mortar

#if !(USE_HDG)
SUBROUTINE Flux_Mortar(Flux_Master,Flux_Slave,doMPISides)
!===================================================================================================================================
! fills master side from small non-conforming sides, Using 1D projection operators M_1_0,M_2_0
!
!> This routine is used to project the numerical flux at the small sides of the nonconforming interface to the corresponding large
!>  ones.
!>
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
USE MOD_Preproc
USE MOD_Mortar_Vars, ONLY: M_1_0,M_2_0
USE MOD_Mesh_Vars,   ONLY: MortarType,MortarInfo,nSides
USE MOD_Mesh_Vars,   ONLY: firstMortarInnerSide,lastMortarInnerSide,FS2M
USE MOD_Mesh_Vars,   ONLY: firstMortarMPISide,lastMortarMPISide
USE MOD_PML_vars,    ONLY:PMLnVar
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT) :: Flux_Master(1:PP_nVar+PMLnVar,0:PP_N,0:PP_N,1:nSides)
REAL,INTENT(INOUT) :: Flux_Slave(1:PP_nVar+PMLnVar,0:PP_N,0:PP_N,1:nSides)
LOGICAL,INTENT(IN) :: doMPISides
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: p,q,l
INTEGER  :: iMortar,nMortars
INTEGER  :: firstMortarSideID,lastMortarSideID
INTEGER  :: MortarSideID,SideID,iSide,flip
REAL         :: Flux_tmp( PP_nVar+PMLnVar,0:PP_N,0:PP_N,1:4)
REAL         :: Flux_tmp2(PP_nVar+PMLnVar,0:PP_N,0:PP_N,1:2)
REAL,POINTER :: M1(:,:),M2(:,:)
!===================================================================================================================================

firstMortarSideID = MERGE(firstMortarMPISide,firstMortarInnerSide,doMPISides)
 lastMortarSideID = MERGE( lastMortarMPISide, lastMortarInnerSide,doMPISides)

M1=>M_1_0; M2=>M_2_0
DO MortarSideID=firstMortarSideID,lastMortarSideID

  nMortars=MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
  iSide=MortarType(2,MortarSideID)
  DO iMortar=1,nMortars
    SideID = MortarInfo(MI_SIDEID,iMortar,iSide)
    flip   = MortarInfo(MI_FLIP,iMortar,iSide)
    SELECT CASE(flip)
    CASE(0) ! master side
      Flux_tmp(:,:,:,iMortar)=Flux_Master(:,:,:,SideID)
    CASE(1:4) ! slave sides (should only occur for MPI)
      DO q=0,PP_N; DO p=0,PP_N
        Flux_tmp(:,FS2M(1,p,q,flip),FS2M(2,p,q,flip),iMortar)=-Flux_Slave(:,p,q,SideID)
      END DO; END DO
    END SELECT
  END DO

  SELECT CASE(MortarType(1,MortarSideID))
  CASE(1) !1->4
    ! first in xi
    DO q=0,PP_N ! for every eta-layer perform Mortar operation in xi-direction
      ! The following p- and l-loop are four MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
      !    Flux_tmp2(iVar,:,q,1)  =  M1 * Flux_tmp(iVar,:,q,1) + M2 * Flux_tmp(iVar,:,q,2)
      !    Flux_tmp2(iVar,:,q,2)  =  M1 * Flux_tmp(iVar,:,q,1) + M2 * Flux_tmp(iVar,:,q,2)
      DO p=0,PP_N
        Flux_tmp2(:,p,q,1)=                       M1(0,p)*Flux_tmp(:,0,q,1)+M2(0,p)*Flux_tmp(:,0,q,2)
        Flux_tmp2(:,p,q,2)=                       M1(0,p)*Flux_tmp(:,0,q,3)+M2(0,p)*Flux_tmp(:,0,q,4)
        DO l=1,PP_N
          Flux_tmp2(:,p,q,1)=Flux_tmp2(:,p,q,1) + M1(l,p)*Flux_tmp(:,l,q,1)+M2(l,p)*Flux_tmp(:,l,q,2)
          Flux_tmp2(:,p,q,2)=Flux_tmp2(:,p,q,2) + M1(l,p)*Flux_tmp(:,l,q,3)+M2(l,p)*Flux_tmp(:,l,q,4)
        END DO
      END DO
    END DO
    !then in eta
    ! The following q- and l-loop are two MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
    !    Flux(iVar,p,:,MortarSideID)  =  M1 * Flux_tmp2(iVar,p,:,1) + M2 * Flux_tmp2(iVar,p,:,2)
    DO q=0,PP_N
      DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction
        Flux_Master(:,p,q,MortarSideID)=                               M1(0,q)*Flux_tmp2(:,p,0,1)+M2(0,q)*Flux_tmp2(:,p,0,2)
        DO l=1,PP_N
          Flux_Master(:,p,q,MortarSideID)=Flux_Master(:,p,q,MortarSideID) + M1(l,q)*Flux_tmp2(:,p,l,1)+M2(l,q)*Flux_tmp2(:,p,l,2)
        END DO
      END DO
    END DO

  CASE(2) !1->2 in eta
    ! TODO why not q-loop first?
    DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction
      ! The following q- and l-loop are two MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
      !    Flux(iVar,p,:,MortarSideID)  =  M1 * Flux_tmp(iVar,p,:,1) + M2 * Flux_tmp(iVar,p,:,2)
      DO q=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction
        Flux_Master(:,p,q,MortarSideID)=                                  M1(0,q)*Flux_tmp(:,p,0,1)+M2(0,q)*Flux_tmp(:,p,0,2)
        DO l=1,PP_N
          Flux_Master(:,p,q,MortarSideID)=Flux_Master(:,p,q,MortarSideID) + M1(l,q)*Flux_tmp(:,p,l,1)+M2(l,q)*Flux_tmp(:,p,l,2)
        END DO
      END DO
    END DO

  CASE(3) !1->2 in xi
    DO q=0,PP_N ! for every eta-layer perform Mortar operation in xi-direction
      ! The following p- and l-loop are two MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
      !    Flux(iVar,:,q,MortarSideID)  =   M1 * Flux_tmp(iVar,:,q,1) + M2 * Flux_tmp(iVar,:,q,2)
      DO p=0,PP_N
        Flux_Master(:,p,q,MortarSideID)=                                    M1(0,p)*Flux_tmp(:,0,q,1)+M2(0,p)*Flux_tmp(:,0,q,2)
        DO l=1,PP_N
          Flux_Master(:,p,q,MortarSideID)=Flux_Master(:,p,q,MortarSideID) + M1(l,p)*Flux_tmp(:,l,q,1)+M2(l,p)*Flux_tmp(:,l,q,2)
        END DO
      END DO
    END DO

  END SELECT ! mortarType(MortarSideID)
END DO !MortarSideID
END SUBROUTINE Flux_Mortar
#endif /*USE_HDG*/

END MODULE MOD_FillMortar
