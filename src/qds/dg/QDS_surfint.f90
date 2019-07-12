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

MODULE MOD_QDS_SurfInt
!===================================================================================================================================
!> Contains the routines to
!> - determine the surface integral for the QDS-DG method
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
INTERFACE SurfIntQDS
  MODULE PROCEDURE SurfIntQDS
END INTERFACE

PUBLIC::SurfIntQDS
!===================================================================================================================================
CONTAINS
SUBROUTINE SurfIntQDS(Flux_Master,Flux_Slave,Ut,doMPISides)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
#if (PP_NodeType>1)
USE MOD_DG_Vars,            ONLY:L_HatPlus,L_HatMinus
#endif
USE MOD_Mesh_Vars,          ONLY:SideToElem
USE MOD_Mesh_Vars,          ONLY:nSides
USE MOD_Mesh_Vars,          ONLY:firstMPISide_YOUR,lastMPISide_MINE
USE MOD_QDS_DG_Vars,        ONLY:nQDSElems
USE MOD_QDS_Equation_vars,  ONLY:QDSnVar
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides+InnerSides+MPISides MINE
REAL,INTENT(IN)    :: Flux_Master(1:QDSnVar,0:PP_N,0:PP_N,nSides)
REAL,INTENT(IN)    :: Flux_Slave(1:QDSnVar,0:PP_N,0:PP_N,nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)   :: Ut(QDSnVar,0:PP_N,0:PP_N,0:PP_N,1:nQDSElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: ElemID,Flip,SideID,locSideID
INTEGER            :: firstSideID,lastSideID
#if (PP_NodeType>1)
REAL               ::L_HatMinus0,L_HatPlusN
#endif
!===================================================================================================================================

IF(doMPISides)THEN
  ! MPI YOUR
  firstSideID = firstMPISide_YOUR
   lastSideID = nSides
ELSE
  ! inner sides and MPI mine
  firstSideID = 1
   lastSideID = lastMPISide_MINE
END IF

#if (PP_NodeType>1)
L_HatMinus0 = L_HatMinus(0)
L_HatPlusN  = L_HatPlus(PP_N)
#endif
DO SideID=firstSideID,lastSideID
  ! neighbor side
  ElemID    = SideToElem(S2E_NB_ELEM_ID,SideID)
  locSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  flip      = SideToElem(S2E_FLIP,SideID)
  ! ignore MPI-faces and boundary faces
  IF(ElemID.LT.0) CYCLE ! boundary side is BC or MPI side
  CALL CalcSurfInt2(Flux_Slave(1:QDSnVar,0:PP_N,0:PP_N,SideID),Ut,Flip,ElemID,locSideID)
END DO ! SideID=1,nSides


DO SideID=firstSideID,lastSideID
  ! master side, flip=0
  ElemID    = SideToElem(S2E_ELEM_ID,SideID)
  locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
  flip      = 0
  IF(ElemID.LT.0) CYCLE ! if master is MPI side
  CALL CalcSurfInt2(Flux_Master(1:QDSnVar,0:PP_N,0:PP_N,SideID),Ut,Flip,ElemID,locSideID)
END DO ! SideID=1,nSides

END SUBROUTINE SurfIntQDS


SUBROUTINE CalcSurfInt2(Flux,Ut,flip,ElemID,locSideID)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,            ONLY:L_HatPlus,L_HatMinus
USE MOD_QDS_DG_Vars,        ONLY:nQDSElems
USE MOD_QDS_Equation_vars,  ONLY:QDSnVar
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: Flux(1:QDSnVar,0:PP_N,0:PP_N)
INTEGER,INTENT(IN) :: flip,ElemID,locSideID!,SideID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)   :: Ut(QDSnVar,0:PP_N,0:PP_N,0:PP_N,1:nQDSElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q,l
#if (PP_NodeType>1)
REAL            ::L_HatMinus0,L_HatPlusN
#endif
!===================================================================================================================================
#if (PP_NodeType>1)
L_HatMinus0 = L_HatMinus(0)
L_HatPlusN  = L_HatPlus(PP_N)
#endif

#if (PP_NodeType==1)
  SELECT CASE(locSideID)
!===================================================================================================================================
  CASE(XI_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut( :,l,p,q,ElemID)=Ut(:,l,p,q,ElemID)+Flux(:,q,p)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut( :,l,p,q,ElemID)=Ut(:,l,p,q,ElemID)-Flux(:,p,q)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut( :,l,p,q,ElemID)=Ut(:,l,p,q,ElemID)-Flux(:,PP_N-q,p)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut( :,l,p,q,ElemID)=Ut(:,l,p,q,ElemID)-Flux(:,PP_N-p,PP_N-q)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut( :,l,p,q,ElemID)=Ut(:,l,p,q,ElemID)-Flux(:,q,PP_N-p)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    END SELECT
!===================================================================================================================================
  CASE(ETA_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut( :,p,l,q,ElemID)=Ut(:,p,l,q,ElemID)+Flux(:,p,q)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut( :,p,l,q,ElemID)=Ut(:,p,l,q,ElemID)-Flux(:,q,p)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut( :,p,l,q,ElemID)=Ut(:,p,l,q,ElemID)-Flux(:,PP_N-p,q)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut( :,p,l,q,ElemID)=Ut(:,p,l,q,ElemID)-Flux(:,PP_N-q,PP_N-p)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut( :,p,l,q,ElemID)=Ut(:,p,l,q,ElemID)-Flux(:,p,PP_N-q)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    END SELECT
!===================================================================================================================================
  CASE(ZETA_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)= Ut(:,p,q,l,ElemID)+Flux(:,q,p)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(1) ! slave side, SideID=q,jSide=p
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)= Ut(:,p,q,l,ElemID)-Flux(:,p,q)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)= Ut(:,p,q,l,ElemID)-Flux(:,PP_N-q,p)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)= Ut(:,p,q,l,ElemID)-Flux(:,PP_N-p,PP_N-q)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)= Ut(:,p,q,l,ElemID)-Flux(:,q,PP_N-p)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    END SELECT
!===================================================================================================================================
  CASE(XI_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut(:,l,p,q,ElemID)=Ut(:,l,p,q,ElemID)+Flux(:,p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut(:,l,p,q,ElemID)=Ut(:,l,p,q,ElemID)-Flux(:,q,p)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut(:,l,p,q,ElemID)=Ut(:,l,p,q,ElemID)-Flux(:,PP_N-p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut(:,l,p,q,ElemID)=Ut(:,l,p,q,ElemID)-Flux(:,PP_N-q,PP_N-p)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut(:,l,p,q,ElemID)=Ut(:,l,p,q,ElemID)-Flux(:,p,PP_N-q)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    END SELECT
!===================================================================================================================================
  CASE(ETA_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut(:,p,l,q,ElemID)=Ut(:,p,l,q,ElemID)+Flux(:,PP_N-p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut(:,p,l,q,ElemID)=Ut(:,p,l,q,ElemID)-Flux(:,q,PP_N-p)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut(:,p,l,q,ElemID)=Ut(:,p,l,q,ElemID)-Flux(:,p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut(:,p,l,q,ElemID)=Ut(:,p,l,q,ElemID)-Flux(:,PP_N-q,p)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut(:,p,l,q,ElemID)=Ut(:,p,l,q,ElemID)-Flux(:,PP_N-p,PP_N-q)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    END SELECT
!===================================================================================================================================
  CASE(ZETA_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)=Ut(:,p,q,l,ElemID)+Flux(:,p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(1) ! slave side, SideID=q,jSide=p
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)=Ut(:,p,q,l,ElemID)-Flux(:,q,p)*L_hatPlus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)=Ut(:,p,q,l,ElemID)-Flux(:,PP_N-p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)=Ut(:,p,q,l,ElemID)-Flux(:,PP_N-q,PP_N-p)*L_hatPlus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)=Ut(:,p,q,l,ElemID)-Flux(:,p,PP_N-q)*L_hatPlus(l)
      END DO; END DO; END DO ! p,q,l
    END SELECT
  END SELECT !locSideID
#else
  !update local grid cell
  SELECT CASE(locSideID)
!===================================================================================================================================
  CASE(XI_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,0,p,q,ElemID)=Ut(:,0,p,q,ElemID)+Flux(:,q,p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,0,p,q,ElemID)=Ut(:,0,p,q,ElemID)-Flux(:,p,q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,0,p,q,ElemID)=Ut(:,0,p,q,ElemID)-Flux(:,PP_N-q,p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,0,p,q,ElemID)=Ut(:,0,p,q,ElemID)-Flux(:,PP_N-p,PP_N-q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,0,p,q,ElemID)=Ut(:,0,p,q,ElemID)-Flux(:,q,PP_N-p)*L_hatMinus0
      END DO; END DO ! p,q
    END SELECT

  ! switch to right hand system for ETA_PLUS direction
!===================================================================================================================================
  CASE(ETA_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,0,q,ElemID)=Ut(:,p,0,q,ElemID)+Flux(:,p,q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,0,q,ElemID)=Ut(:,p,0,q,ElemID)-Flux(:,q,p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,0,q,ElemID)=Ut(:,p,0,q,ElemID)-Flux(:,PP_N-p,q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,0,q,ElemID)=Ut(:,p,0,q,ElemID)-Flux(:,PP_N-q,PP_N-p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,0,q,ElemID)=Ut(:,p,0,q,ElemID)-Flux(:,p,PP_N-q)*L_hatMinus0
      END DO; END DO ! p,q
    END SELECT

  ! switch to right hand system for ZETA_MINUS direction
!===================================================================================================================================
  CASE(ZETA_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,0,ElemID)=Ut(:,p,q,0,ElemID)+Flux(:,q,p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,0,ElemID)=Ut(:,p,q,0,ElemID)-Flux(:,p,q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,0,ElemID)=Ut(:,p,q,0,ElemID)-Flux(:,PP_N-q,p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,0,ElemID)=Ut(:,p,q,0,ElemID)-Flux(:,PP_N-p,PP_N-q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,0,ElemID)=Ut(:,p,q,0,ElemID)-Flux(:,q,PP_N-p)*L_hatMinus0
      END DO; END DO ! p,q
    END SELECT

!===================================================================================================================================
  CASE(XI_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,PP_N,p,q,ElemID)=Ut(:,PP_N,p,q,ElemID)+Flux(:,p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,PP_N,p,q,ElemID)=Ut(:,PP_N,p,q,ElemID)-Flux(:,q,p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,PP_N,p,q,ElemID)=Ut(:,PP_N,p,q,ElemID)-Flux(:,PP_N-p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,PP_N,p,q,ElemID)=Ut(:,PP_N,p,q,ElemID)-Flux(:,PP_N-q,PP_N-p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,PP_N,p,q,ElemID)=Ut(:,PP_N,p,q,ElemID)-Flux(:,p,PP_N-q)*L_hatPlusN
      END DO; END DO ! p,q
    END SELECT

  ! switch to right hand system for ETA_PLUS direction
!===================================================================================================================================
  CASE(ETA_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,PP_N,q,ElemID)=Ut(:,p,PP_N,q,ElemID)+Flux(:,PP_N-p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,PP_N,q,ElemID)=Ut(:,p,PP_N,q,ElemID)-Flux(:,q,PP_N-p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,PP_N,q,ElemID)=Ut(:,p,PP_N,q,ElemID)-Flux(:,p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,PP_N,q,ElemID)=Ut(:,p,PP_N,q,ElemID)-Flux(:,PP_N-q,p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,PP_N,q,ElemID)=Ut(:,p,PP_N,q,ElemID)-Flux(:,PP_N-p,PP_N-q)*L_hatPlusN
      END DO; END DO ! p,q
    END SELECT

  ! switch to right hand system for ZETA_MINUS direction
!===================================================================================================================================
  CASE(ZETA_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,PP_N,ElemID)=Ut(:,p,q,PP_N,ElemID)+Flux(:,p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,PP_N,ElemID)=Ut(:,p,q,PP_N,ElemID)-Flux(:,q,p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,PP_N,ElemID)=Ut(:,p,q,PP_N,ElemID)-Flux(:,PP_N-p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,PP_N,ElemID)=Ut(:,p,q,PP_N,ElemID)-Flux(:,PP_N-q,PP_N-p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,PP_N,ElemID)=Ut(:,p,q,PP_N,ElemID)-Flux(:,p,PP_N-q)*L_hatPlusN
      END DO; END DO ! p,q
    END SELECT
  END SELECT !locSideID
#endif
END SUBROUTINE CalcSurfInt2


END MODULE MOD_QDS_SurfInt
#endif /*USE_QDS_DG*/
