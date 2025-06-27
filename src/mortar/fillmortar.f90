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
#if !(USE_FV) || (USE_HDG)
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

#if !(USE_HDG) || (USE_FV)
PUBLIC:: Flux_Mortar
#endif /*USE_HDG*/

PUBLIC::U_Mortar

!===================================================================================================================================

CONTAINS

SUBROUTINE U_Mortar(doDielectricSides, doMPISides)
!===================================================================================================================================
!> fills small non-conforming sides with data for master side with data from the corresponding large side, using 1D interpolation
!> operators M_0_1,M_0_2

!     Type 1               Type 2              Type3
!      eta                  eta                 eta
!       ^                    ^                   ^
!       |                    |                   |
!   +---+---+            +---+---+           +---+---+
!   | 3 | 4 |            |   2   |           |   |   |
!   +---+---+ --->  xi   +---+---+ --->  xi  + 1 + 2 + --->  xi
!   | 1 | 2 |            |   1   |           |   |   |
!   +---+---+            +---+---+           +---+---+

!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Mortar_Vars     ,ONLY: N_Mortar
USE MOD_Mesh_Vars       ,ONLY: MortarType,MortarInfo, offSetElem, SideToElem
USE MOD_Mesh_Vars       ,ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_Mesh_Vars       ,ONLY: firstMortarMPISide,lastMortarMPISide
USE MOD_Mesh_Vars       ,ONLY: N_Mesh
USE MOD_DG_Vars         ,ONLY: N_DG_Mapping, U_Surf_N
USE MOD_Dielectric_Vars ,ONLY: DielectricSurf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: doDielectricSides  !  .TRUE.: use DielectricSurf(:)%Dielectric_dummy_Master and
!                                        !              DielectricSurf(:)%Dielectric_dummy_Slave
!                                        ! .FALSE.: use U_Surf_N(:)%U_master and
!                                        !              U_Surf_N(:)%U_slave
LOGICAL,INTENT(IN) :: doMPISides         !< flag whether MPI sides are processed
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER      :: p,q,l, ElemID, Nloc
INTEGER      :: iMortar,nMortars
INTEGER      :: firstMortarSideID,lastMortarSideID
INTEGER      :: MortarSideID,SideID,locSide,flip
!===================================================================================================================================

firstMortarSideID = MERGE(firstMortarMPISide,firstMortarInnerSide,doMPISides)
lastMortarSideID = MERGE( lastMortarMPISide, lastMortarInnerSide,doMPISides)

DO MortarSideID=firstMortarSideID,lastMortarSideID
  ElemID    = SideToElem(S2E_ELEM_ID,MortarSideID)
  Nloc = N_DG_Mapping(2,ElemID+offSetElem)
  SELECT CASE(MortarType(1,MortarSideID))
  CASE(1) !1->4
   !first in eta
   ! The following q- and l-loop are two MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
   !    U_tmp2(iVar,p,:,1)  =  M1 * U_in_master(iVar,p,:,MortarSideID)
   !    U_tmp2(iVar,p,:,2)  =  M2 * U_in_master(iVar,p,:,MortarSideID)
   IF(doDielectricSides)THEN
     N_Mortar(Nloc)%U_tmp  = 0.
     N_Mortar(Nloc)%U_tmp2 = 0.
     DO q=0,Nloc
       DO p=0,Nloc ! for every xi-layer perform Mortar operation in eta-direction
         N_Mortar(Nloc)%U_tmp2(1,p,q,1)= N_Mortar(Nloc)%M_0_1(0,q)*DielectricSurf(MortarSideID)%Dielectric_dummy_Master(1,p,0)
         N_Mortar(Nloc)%U_tmp2(1,p,q,2)= N_Mortar(Nloc)%M_0_2(0,q)*DielectricSurf(MortarSideID)%Dielectric_dummy_Master(1,p,0)
         DO l=1,Nloc
           N_Mortar(Nloc)%U_tmp2(1,p,q,1)=N_Mortar(Nloc)%U_tmp2(1,p,q,1)+N_Mortar(Nloc)%M_0_1(l,q)*DielectricSurf(MortarSideID)%Dielectric_dummy_Master(1,p,l)
           N_Mortar(Nloc)%U_tmp2(1,p,q,2)=N_Mortar(Nloc)%U_tmp2(1,p,q,2)+N_Mortar(Nloc)%M_0_2(l,q)*DielectricSurf(MortarSideID)%Dielectric_dummy_Master(1,p,l)
         END DO
       END DO
     END DO
   ELSE
     DO q=0,Nloc
       DO p=0,Nloc ! for every xi-layer perform Mortar operation in eta-direction
         N_Mortar(Nloc)%U_tmp2(1:PP_nVar,p,q,1)= N_Mortar(Nloc)%M_0_1(0,q)*U_Surf_N(MortarSideID)%U_master(:,p,0)
         N_Mortar(Nloc)%U_tmp2(1:PP_nVar,p,q,2)= N_Mortar(Nloc)%M_0_2(0,q)*U_Surf_N(MortarSideID)%U_master(:,p,0)
         DO l=1,Nloc
           N_Mortar(Nloc)%U_tmp2(1:PP_nVar,p,q,1)=N_Mortar(Nloc)%U_tmp2(1:PP_nVar,p,q,1)+&
                                                  N_Mortar(Nloc)%M_0_1(l,q)*U_Surf_N(MortarSideID)%U_master(:,p,l)
           N_Mortar(Nloc)%U_tmp2(1:PP_nVar,p,q,2)=N_Mortar(Nloc)%U_tmp2(1:PP_nVar,p,q,2)+&
                                                  N_Mortar(Nloc)%M_0_2(l,q)*U_Surf_N(MortarSideID)%U_master(:,p,l)
         END DO
       END DO
     END DO
   END IF ! doDielectricSides
   ! then in xi
   DO q=0,Nloc ! for every eta-layer perform Mortar operation in xi-direction
     ! The following p- and l-loop are four MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
     !    U_tmp(iVar,:,q,1)  =  M1 * U_tmp2(iVar,:,q,1)
     !    U_tmp(iVar,:,q,2)  =  M2 * U_tmp2(iVar,:,q,1)
     !    U_tmp(iVar,:,q,3)  =  M1 * U_tmp2(iVar,:,q,2)
     !    U_tmp(iVar,:,q,4)  =  M2 * U_tmp2(iVar,:,q,2)
     DO p=0,Nloc
       N_Mortar(Nloc)%U_tmp(:,p,q,1)= N_Mortar(Nloc)%M_0_1(0,p)*N_Mortar(Nloc)%U_tmp2(:,0,q,1)
       N_Mortar(Nloc)%U_tmp(:,p,q,2)= N_Mortar(Nloc)%M_0_2(0,p)*N_Mortar(Nloc)%U_tmp2(:,0,q,1)
       N_Mortar(Nloc)%U_tmp(:,p,q,3)= N_Mortar(Nloc)%M_0_1(0,p)*N_Mortar(Nloc)%U_tmp2(:,0,q,2)
       N_Mortar(Nloc)%U_tmp(:,p,q,4)= N_Mortar(Nloc)%M_0_2(0,p)*N_Mortar(Nloc)%U_tmp2(:,0,q,2)
       DO l=1,Nloc
         N_Mortar(Nloc)%U_tmp(:,p,q,1)=N_Mortar(Nloc)%U_tmp(:,p,q,1)+N_Mortar(Nloc)%M_0_1(l,p)*N_Mortar(Nloc)%U_tmp2(:,l,q,1)
         N_Mortar(Nloc)%U_tmp(:,p,q,2)=N_Mortar(Nloc)%U_tmp(:,p,q,2)+N_Mortar(Nloc)%M_0_2(l,p)*N_Mortar(Nloc)%U_tmp2(:,l,q,1)
         N_Mortar(Nloc)%U_tmp(:,p,q,3)=N_Mortar(Nloc)%U_tmp(:,p,q,3)+N_Mortar(Nloc)%M_0_1(l,p)*N_Mortar(Nloc)%U_tmp2(:,l,q,2)
         N_Mortar(Nloc)%U_tmp(:,p,q,4)=N_Mortar(Nloc)%U_tmp(:,p,q,4)+N_Mortar(Nloc)%M_0_2(l,p)*N_Mortar(Nloc)%U_tmp2(:,l,q,2)
       END DO !l=1,Nloc
     END DO
   END DO

  CASE(2) !1->2 in eta
   ! The following q- and l-loop are two MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
   !    U_tmp(iVar,p,:,1)  =  M1 * U_in_master(iVar,p,:,MortarSideID)
   !    U_tmp(iVar,p,:,2)  =  M2 * U_in_master(iVar,p,:,MortarSideID)
   IF(doDielectricSides)THEN
     DO q=0,Nloc
       DO p=0,Nloc ! for every xi-layer perform Mortar operation in eta-direction
         N_Mortar(Nloc)%U_tmp(1,p,q,1)= N_Mortar(Nloc)%M_0_1(0,q)*DielectricSurf(MortarSideID)%Dielectric_dummy_Master(1,p,0)
         N_Mortar(Nloc)%U_tmp(1,p,q,2)= N_Mortar(Nloc)%M_0_2(0,q)*DielectricSurf(MortarSideID)%Dielectric_dummy_Master(1,p,0)
         DO l=1,Nloc
           N_Mortar(Nloc)%U_tmp(1,p,q,1)=N_Mortar(Nloc)%U_tmp(1,p,q,1)+N_Mortar(Nloc)%M_0_1(l,q)*DielectricSurf(MortarSideID)%Dielectric_dummy_Master(1,p,l)
           N_Mortar(Nloc)%U_tmp(1,p,q,2)=N_Mortar(Nloc)%U_tmp(1,p,q,2)+N_Mortar(Nloc)%M_0_2(l,q)*DielectricSurf(MortarSideID)%Dielectric_dummy_Master(1,p,l)
         END DO
       END DO
     END DO
   ELSE
     DO q=0,Nloc
       DO p=0,Nloc ! for every xi-layer perform Mortar operation in eta-direction
         N_Mortar(Nloc)%U_tmp(1:PP_nVar,p,q,1)= N_Mortar(Nloc)%M_0_1(0,q)*U_Surf_N(MortarSideID)%U_master(:,p,0)
         N_Mortar(Nloc)%U_tmp(1:PP_nVar,p,q,2)= N_Mortar(Nloc)%M_0_2(0,q)*U_Surf_N(MortarSideID)%U_master(:,p,0)
         DO l=1,Nloc
           N_Mortar(Nloc)%U_tmp(1:PP_nVar,p,q,1)=N_Mortar(Nloc)%U_tmp(1:PP_nVar,p,q,1)+&
                                                 N_Mortar(Nloc)%M_0_1(l,q)*U_Surf_N(MortarSideID)%U_master(:,p,l)
           N_Mortar(Nloc)%U_tmp(1:PP_nVar,p,q,2)=N_Mortar(Nloc)%U_tmp(1:PP_nVar,p,q,2)+&
                                                 N_Mortar(Nloc)%M_0_2(l,q)*U_Surf_N(MortarSideID)%U_master(:,p,l)
         END DO
       END DO
     END DO
   END IF ! doDielectricSides

  CASE(3) !1->2 in xi
   DO q=0,Nloc ! for every eta-layer perform Mortar operation in xi-direction
     ! The following p- and l-loop are two MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
     !    U_tmp(iVar,:,q,1)  =  M1 * U_in_master(iVar,:,q,MortarSideID)
     !    U_tmp(iVar,:,q,2)  =  M2 * U_in_master(iVar,:,q,MortarSideID)
     IF(doDielectricSides)THEN
       DO p=0,Nloc
         N_Mortar(Nloc)%U_tmp(1,p,q,1)= N_Mortar(Nloc)%M_0_1(0,p)*DielectricSurf(MortarSideID)%Dielectric_dummy_Master(1,0,q)
         N_Mortar(Nloc)%U_tmp(1,p,q,2)= N_Mortar(Nloc)%M_0_2(0,p)*DielectricSurf(MortarSideID)%Dielectric_dummy_Master(1,0,q)
         DO l=1,Nloc
           N_Mortar(Nloc)%U_tmp(1,p,q,1)=N_Mortar(Nloc)%U_tmp(1,p,q,1)+&
                                         N_Mortar(Nloc)%M_0_1(l,p)*DielectricSurf(MortarSideID)%Dielectric_dummy_Master(1,l,q)
           N_Mortar(Nloc)%U_tmp(1,p,q,2)=N_Mortar(Nloc)%U_tmp(1,p,q,2)+&
                                         N_Mortar(Nloc)%M_0_2(l,p)*DielectricSurf(MortarSideID)%Dielectric_dummy_Master(1,l,q)
         END DO
       END DO
     ELSE
       DO p=0,Nloc
         N_Mortar(Nloc)%U_tmp(1:PP_nVar,p,q,1)= N_Mortar(Nloc)%M_0_1(0,p)*U_Surf_N(MortarSideID)%U_master(:,0,q)
         N_Mortar(Nloc)%U_tmp(1:PP_nVar,p,q,2)= N_Mortar(Nloc)%M_0_2(0,p)*U_Surf_N(MortarSideID)%U_master(:,0,q)
         DO l=1,Nloc
           N_Mortar(Nloc)%U_tmp(1:PP_nVar,p,q,1)=N_Mortar(Nloc)%U_tmp(1:PP_nVar,p,q,1)+&
                                                 N_Mortar(Nloc)%M_0_1(l,p)*U_Surf_N(MortarSideID)%U_master(:,l,q)
           N_Mortar(Nloc)%U_tmp(1:PP_nVar,p,q,2)=N_Mortar(Nloc)%U_tmp(1:PP_nVar,p,q,2)+&
                                                 N_Mortar(Nloc)%M_0_2(l,p)*U_Surf_N(MortarSideID)%U_master(:,l,q)
         END DO
       END DO
     END IF ! doDielectricSides
   END DO
  END SELECT ! mortarType(SideID)

  nMortars=MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
  locSide=MortarType(2,MortarSideID)
  DO iMortar=1,nMortars
   SideID= MortarInfo(MI_SIDEID,iMortar,locSide)
   flip  = MortarInfo(MI_FLIP,iMortar,locSide)
   IF(doDielectricSides)THEN
     SELECT CASE(flip)
       CASE(0) ! master side
         DielectricSurf(SideID)%Dielectric_dummy_Master(1,:,:)=N_Mortar(Nloc)%U_tmp(1,:,:,iMortar)
       CASE(1:4) ! slave side
         DO q=0,Nloc; DO p=0,Nloc
           DielectricSurf(SideID)%Dielectric_dummy_Slave(1,p,q)=N_Mortar(Nloc)%U_tmp(1,N_Mesh(Nloc)%FS2M(1,p,q,flip), &
                                                                                       N_Mesh(Nloc)%FS2M(2,p,q,flip),iMortar)
         END DO; END DO ! q, p
     END SELECT !flip(iMortar)
   ELSE
     SELECT CASE(flip)
       CASE(0) ! master side
         U_Surf_N(SideID)%U_master(:,:,:)=N_Mortar(Nloc)%U_tmp(1:PP_nVar,:,:,iMortar)
       CASE(1:4) ! slave side
         DO q=0,Nloc; DO p=0,Nloc
           U_Surf_N(SideID)%U_slave(:,p,q)=N_Mortar(Nloc)%U_tmp(1:PP_nVar,N_Mesh(Nloc)%FS2M(1,p,q,flip), &
                                                                          N_Mesh(Nloc)%FS2M(2,p,q,flip),iMortar)
         END DO; END DO ! q, p
     END SELECT !flip(iMortar)
   END IF ! doDielectricSides
  END DO !iMortar
END DO !MortarSideID
END SUBROUTINE U_Mortar

#if !(USE_HDG) || (USE_FV)
SUBROUTINE Flux_Mortar(doMPISides)
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
USE MOD_Mortar_Vars, ONLY: N_Mortar
USE MOD_Mesh_Vars,   ONLY: MortarType,MortarInfo,offSetElem,SideToElem
USE MOD_Mesh_Vars,   ONLY: firstMortarInnerSide,lastMortarInnerSide,N_Mesh
USE MOD_Mesh_Vars,   ONLY: firstMortarMPISide,lastMortarMPISide
USE MOD_DG_Vars,     ONLY: N_DG_Mapping, U_Surf_N
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: p,q,l, ElemID, Nloc
INTEGER  :: iMortar,nMortars
INTEGER  :: firstMortarSideID,lastMortarSideID
INTEGER  :: MortarSideID,SideID,iSide,flip
!===================================================================================================================================

firstMortarSideID = MERGE(firstMortarMPISide,firstMortarInnerSide,doMPISides)
lastMortarSideID  = MERGE( lastMortarMPISide, lastMortarInnerSide,doMPISides)

DO MortarSideID=firstMortarSideID,lastMortarSideID
  ElemID   = SideToElem(S2E_ELEM_ID,MortarSideID)
  Nloc     = N_DG_Mapping(2,ElemID+offSetElem)
  nMortars = MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
  iSide    = MortarType(2,MortarSideID)

  ! Loop Mortar faces
  DO iMortar=1,nMortars
    SideID = MortarInfo(MI_SIDEID,iMortar,iSide)
    flip   = MortarInfo(MI_FLIP,iMortar,iSide)
    SELECT CASE(flip)
    CASE(0) ! master side
      N_Mortar(Nloc)%U_tmp(:,:,:,iMortar) = U_Surf_N(SideID)%Flux_Master(:,:,:)
    CASE(1:4) ! slave sides (should only occur for MPI)
      DO q=0,Nloc; DO p=0,Nloc
        N_Mortar(Nloc)%U_tmp(:,N_Mesh(Nloc)%FS2M(1,p,q,flip),N_Mesh(Nloc)%FS2M(2,p,q,flip),iMortar)=-U_Surf_N(SideID)%Flux_Slave(:,p,q)
      END DO; END DO
    END SELECT
  END DO

  SELECT CASE(MortarType(1,MortarSideID))
  CASE(1) !1->4
    ! first in xi
    DO q=0,Nloc ! for every eta-layer perform Mortar operation in xi-direction
      ! The following p- and l-loop are four MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
      !    Flux_tmp2(iVar,:,q,1)  =  M1 * Flux_tmp(iVar,:,q,1) + M2 * Flux_tmp(iVar,:,q,2)
      !    Flux_tmp2(iVar,:,q,2)  =  M1 * Flux_tmp(iVar,:,q,1) + M2 * Flux_tmp(iVar,:,q,2)
      DO p=0,Nloc
        N_Mortar(Nloc)%U_tmp2(:,p,q,1) = N_Mortar(Nloc)%M_1_0(0,p)*N_Mortar(Nloc)%U_tmp(:,0,q,1) + &
                                         N_Mortar(Nloc)%M_2_0(0,p)*N_Mortar(Nloc)%U_tmp(:,0,q,2)
        N_Mortar(Nloc)%U_tmp2(:,p,q,2) = N_Mortar(Nloc)%M_1_0(0,p)*N_Mortar(Nloc)%U_tmp(:,0,q,3) + &
                                         N_Mortar(Nloc)%M_2_0(0,p)*N_Mortar(Nloc)%U_tmp(:,0,q,4)
        DO l=1,Nloc
          N_Mortar(Nloc)%U_tmp2(:,p,q,1)=N_Mortar(Nloc)%U_tmp2(:,p,q,1) + &
               N_Mortar(Nloc)%M_1_0(l,p)*N_Mortar(Nloc)%U_tmp( :,l,q,1) + &
               N_Mortar(Nloc)%M_2_0(l,p)*N_Mortar(Nloc)%U_tmp( :,l,q,2)
          N_Mortar(Nloc)%U_tmp2(:,p,q,2)=N_Mortar(Nloc)%U_tmp2(:,p,q,2) + &
               N_Mortar(Nloc)%M_1_0(l,p)*N_Mortar(Nloc)%U_tmp( :,l,q,3) + &
               N_Mortar(Nloc)%M_2_0(l,p)*N_Mortar(Nloc)%U_tmp( :,l,q,4)
        END DO
      END DO
    END DO
    !then in eta
    ! The following q- and l-loop are two MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
    !    Flux(iVar,p,:,MortarSideID)  =  M1 * Flux_tmp2(iVar,p,:,1) + M2 * Flux_tmp2(iVar,p,:,2)
    DO q=0,Nloc
      DO p=0,Nloc ! for every xi-layer perform Mortar operation in eta-direction
        U_Surf_N(MortarSideID)%Flux_Master(:,p,q) = N_Mortar(Nloc)%M_1_0(0,q)*N_Mortar(Nloc)%U_tmp2(:,p,0,1) + &
                                                    N_Mortar(Nloc)%M_2_0(0,q)*N_Mortar(Nloc)%U_tmp2(:,p,0,2)
        DO l=1,Nloc
          U_Surf_N(MortarSideID)%Flux_Master(:,p,q) = U_Surf_N(MortarSideID)%Flux_Master(:,p,q)                + &
                                                      N_Mortar(Nloc)%M_1_0(l,q)*N_Mortar(Nloc)%U_tmp2(:,p,l,1) + &
                                                      N_Mortar(Nloc)%M_2_0(l,q)*N_Mortar(Nloc)%U_tmp2(:,p,l,2)
        END DO
      END DO
    END DO

  CASE(2) !1->2 in eta
    ! TODO why not q-loop first?
    DO p=0,Nloc ! for every xi-layer perform Mortar operation in eta-direction
      ! The following q- and l-loop are two MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
      !    Flux(iVar,p,:,MortarSideID)  =  M1 * Flux_tmp(iVar,p,:,1) + M2 * Flux_tmp(iVar,p,:,2)
      DO q=0,Nloc ! for every xi-layer perform Mortar operation in eta-direction
        U_Surf_N(MortarSideID)%Flux_Master(:,p,q)= N_Mortar(Nloc)%M_1_0(0,q)*N_Mortar(Nloc)%U_tmp(:,p,0,1) + &
                                                   N_Mortar(Nloc)%M_2_0(0,q)*N_Mortar(Nloc)%U_tmp(:,p,0,2)
        DO l=1,Nloc
          U_Surf_N(MortarSideID)%Flux_Master(:,p,q)=U_Surf_N(MortarSideID)%Flux_Master(:,p,q)               + &
                                                    N_Mortar(Nloc)%M_1_0(l,q)*N_Mortar(Nloc)%U_tmp(:,p,l,1) + &
                                                    N_Mortar(Nloc)%M_2_0(l,q)*N_Mortar(Nloc)%U_tmp(:,p,l,2)
        END DO
      END DO
    END DO

  CASE(3) !1->2 in xi
    DO q=0,Nloc ! for every eta-layer perform Mortar operation in xi-direction
      ! The following p- and l-loop are two MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
      !    Flux(iVar,:,q,MortarSideID)  =   M1 * Flux_tmp(iVar,:,q,1) + M2 * Flux_tmp(iVar,:,q,2)
      DO p=0,Nloc
        U_Surf_N(MortarSideID)%Flux_Master(:,p,q) = N_Mortar(Nloc)%M_1_0(0,p)*N_Mortar(Nloc)%U_tmp(:,0,q,1) + &
                                                    N_Mortar(Nloc)%M_2_0(0,p)*N_Mortar(Nloc)%U_tmp(:,0,q,2)
        DO l=1,Nloc
          U_Surf_N(MortarSideID)%Flux_Master(:,p,q) = U_Surf_N(MortarSideID)%Flux_Master(:,p,q)               + &
                                                      N_Mortar(Nloc)%M_1_0(l,p)*N_Mortar(Nloc)%U_tmp(:,l,q,1) + &
                                                      N_Mortar(Nloc)%M_2_0(l,p)*N_Mortar(Nloc)%U_tmp(:,l,q,2)
        END DO
      END DO
    END DO

  END SELECT ! mortarType(MortarSideID)
END DO !MortarSideID
END SUBROUTINE Flux_Mortar
#endif !(USE_HDG) || (USE_FV)

#endif /*!(USE_FV) || (USE_HDG)*/
END MODULE MOD_FillMortar