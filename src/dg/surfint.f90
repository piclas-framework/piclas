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

MODULE MOD_SurfInt
!===================================================================================================================================
! Contains the different Surface integral formulations
! Computes the Surface integral for all faces using U and updates Ut
! Computes only inner surface integrals!
! Surface integrals are separated for each direction
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
INTERFACE SurfInt
  MODULE PROCEDURE SurfInt
END INTERFACE

PUBLIC::SurfInt
!===================================================================================================================================
CONTAINS


SUBROUTINE SurfInt(doMPISides)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars   ,ONLY: U_Surf_N,N_DG_Mapping,DGB_N,U_N
USE MOD_Mesh_Vars ,ONLY: SideToElem, offSetElem
USE MOD_PML_Vars  ,ONLY: DoPML,PMLnVar,isPMLElem
USE MOD_Mesh_Vars ,ONLY: nSides
USE MOD_Mesh_Vars ,ONLY: firstMPISide_YOUR,lastMPISide_MINE
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides+InnerSides+MPISides MINE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: ElemID,Flip,SideID,locSideID,Nloc
INTEGER            :: firstSideID,lastSideID
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

DO SideID=firstSideID,lastSideID
  ! neighbor side
  ElemID    = SideToElem(S2E_NB_ELEM_ID,SideID)
  locSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  flip      = SideToElem(S2E_FLIP,SideID)
  ! ignore MPI-faces and boundary faces
  IF(ElemID.LT.0) CYCLE ! boundary side is BC or MPI side
  Nloc      = N_DG_Mapping(2,ElemID+offSetElem)
  ASSOCIATE( Ut         => U_N(ElemID)%Ut(:,:,:,:) ,&
             L_hatMinus => DGB_N(Nloc)%L_HatMinus  ,&
             L_hatPlus  => DGB_N(Nloc)%L_HatPlus   )
    IF(DoPML)THEN
      IF(isPMLElem(ElemID))THEN
        CALL CalcSurfIntPML(Ut,L_hatMinus,L_hatPlus,Nloc,U_Surf_N(SideID)%Flux_Slave(1:PP_nVar+PMLnVar,0:Nloc,0:Nloc),flip,ElemID,locSideID)
      ELSE
        CALL CalcSurfInt(Ut,L_hatMinus,L_hatPlus,Nloc,U_Surf_N(SideID)%Flux_Slave(1:PP_nVar,0:Nloc,0:Nloc),Flip,ElemID,locSideID)
      END IF
    ELSE
      CALL CalcSurfInt(Ut,L_hatMinus,L_hatPlus,Nloc,U_Surf_N(SideID)%Flux_Slave(1:PP_nVar,0:Nloc,0:Nloc),Flip,ElemID,locSideID)
    END IF
  END ASSOCIATE
END DO ! SideID=1,nSides

DO SideID=firstSideID,lastSideID
  ! master side, flip=0
  ElemID    = SideToElem(S2E_ELEM_ID,SideID)
  locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
  flip      = 0  
  IF(ElemID.LT.0) CYCLE ! if master is MPI side
  Nloc      = N_DG(2,ElemID+offSetElem)
  ASSOCIATE( Ut         => U_N(ElemID)%Ut(:,:,:,:) ,&
             L_hatMinus => DGB_N(Nloc)%L_HatMinus  ,&
             L_hatPlus  => DGB_N(Nloc)%L_HatPlus   )
    IF(DoPML)THEN
      IF(isPMLElem(ElemID))THEN
        CALL CalcSurfIntPML(Ut,L_hatMinus,L_hatPlus,Nloc,U_Surf_N(SideID)%Flux_Master(1:PP_nVar+PMLnVar,0:Nloc,0:Nloc),flip,ElemID,locSideID)
      ELSE
        CALL CalcSurfInt(Ut,L_hatMinus,L_hatPlus,Nloc,U_Surf_N(SideID)%Flux_Master(1:PP_nVar,0:Nloc,0:Nloc),Flip,ElemID,locSideID)
      END IF
    ELSE
      CALL CalcSurfInt(Ut,L_hatMinus,L_hatPlus,Nloc,U_Surf_N(SideID)%Flux_Master(1:PP_nVar,0:Nloc,0:Nloc),Flip,ElemID,locSideID)
    END IF
  END ASSOCIATE
END DO ! SideID=1,nSides

END SUBROUTINE SurfInt

SUBROUTINE CalcSurfInt(Ut,L_hatMinus,L_hatPlus,Nloc,Flux,flip,ElemID,locSideID)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT) :: Ut(1:PP_nVar,0:Nloc,0:Nloc,0:Nloc)
REAL,INTENT(IN)    :: L_hatMinus(0:Nloc)
REAL,INTENT(IN)    :: L_hatPlus(0:Nloc)
INTEGER,INTENT(IN) :: Nloc
REAL,INTENT(IN)    :: Flux(1:PP_nVar,0:Nloc,0:Nloc)
INTEGER,INTENT(IN) :: flip,ElemID,locSideID!,SideID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q,l
!INTEGER            :: firstSideID,lastSideID
#if (PP_NodeType>1)
REAL            ::L_hatMinus0,L_hatPlusN
#endif
!===================================================================================================================================
#if (PP_NodeType>1)
L_hatMinus0 = L_HatMinus(0)
L_hatPlusN  = L_HatPlus(Nloc)
#endif

#if (PP_NodeType==1)
  SELECT CASE(locSideID)
!===================================================================================================================================
  CASE(XI_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO q=0,Nloc; DO p=0,Nloc; DO l=0,Nloc
        Ut( :,l,p,q)=Ut(:,l,p,q)+Flux(:,q,p)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,Nloc; DO p=0,Nloc; DO l=0,Nloc
        Ut( :,l,p,q)=Ut(:,l,p,q)-Flux(:,p,q)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,Nloc; DO p=0,Nloc; DO l=0,Nloc
        Ut( :,l,p,q)=Ut(:,l,p,q)-Flux(:,Nloc-q,p)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,Nloc; DO p=0,Nloc; DO l=0,Nloc
        Ut( :,l,p,q)=Ut(:,l,p,q)-Flux(:,Nloc-p,Nloc-q)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,Nloc; DO p=0,Nloc; DO l=0,Nloc
        Ut( :,l,p,q)=Ut(:,l,p,q)-Flux(:,q,Nloc-p)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    END SELECT
!===================================================================================================================================
  CASE(ETA_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO q=0,Nloc; DO l=0,Nloc; DO p=0,Nloc
        Ut( :,p,l,q)=Ut(:,p,l,q)+Flux(:,p,q)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,Nloc; DO l=0,Nloc; DO p=0,Nloc
        Ut( :,p,l,q)=Ut(:,p,l,q)-Flux(:,q,p)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,Nloc; DO l=0,Nloc; DO p=0,Nloc
        Ut( :,p,l,q)=Ut(:,p,l,q)-Flux(:,Nloc-p,q)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,Nloc; DO l=0,Nloc; DO p=0,Nloc
        Ut( :,p,l,q)=Ut(:,p,l,q)-Flux(:,Nloc-q,Nloc-p)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,Nloc; DO l=0,Nloc; DO p=0,Nloc
        Ut( :,p,l,q)=Ut(:,p,l,q)-Flux(:,p,Nloc-q)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    END SELECT
!===================================================================================================================================
  CASE(ZETA_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO l=0,Nloc; DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,l)= Ut(:,p,q,l)+Flux(:,q,p)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(1) ! slave side, SideID=q,jSide=p
      DO l=0,Nloc; DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,l)= Ut(:,p,q,l)-Flux(:,p,q)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO l=0,Nloc; DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,l)= Ut(:,p,q,l)-Flux(:,Nloc-q,p)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO l=0,Nloc; DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,l)= Ut(:,p,q,l)-Flux(:,Nloc-p,Nloc-q)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO l=0,Nloc; DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,l)= Ut(:,p,q,l)-Flux(:,q,Nloc-p)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    END SELECT
!===================================================================================================================================
  CASE(XI_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO q=0,Nloc; DO p=0,Nloc; DO l=0,Nloc
        Ut(:,l,p,q)=Ut(:,l,p,q)+Flux(:,p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,Nloc; DO p=0,Nloc; DO l=0,Nloc
        Ut(:,l,p,q)=Ut(:,l,p,q)-Flux(:,q,p)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,Nloc; DO p=0,Nloc; DO l=0,Nloc
        Ut(:,l,p,q)=Ut(:,l,p,q)-Flux(:,Nloc-p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,Nloc; DO p=0,Nloc; DO l=0,Nloc
        Ut(:,l,p,q)=Ut(:,l,p,q)-Flux(:,Nloc-q,Nloc-p)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,Nloc; DO p=0,Nloc; DO l=0,Nloc
        Ut(:,l,p,q)=Ut(:,l,p,q)-Flux(:,p,Nloc-q)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    END SELECT
!===================================================================================================================================
  CASE(ETA_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO q=0,Nloc; DO l=0,Nloc; DO p=0,Nloc
        Ut(:,p,l,q)=Ut(:,p,l,q)+Flux(:,Nloc-p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,Nloc; DO l=0,Nloc; DO p=0,Nloc
        Ut(:,p,l,q)=Ut(:,p,l,q)-Flux(:,q,Nloc-p)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,Nloc; DO l=0,Nloc; DO p=0,Nloc
        Ut(:,p,l,q)=Ut(:,p,l,q)-Flux(:,p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,Nloc; DO l=0,Nloc; DO p=0,Nloc
        Ut(:,p,l,q)=Ut(:,p,l,q)-Flux(:,Nloc-q,p)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,Nloc; DO l=0,Nloc; DO p=0,Nloc
        Ut(:,p,l,q)=Ut(:,p,l,q)-Flux(:,Nloc-p,Nloc-q)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    END SELECT
!===================================================================================================================================
  CASE(ZETA_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO l=0,Nloc; DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,l)=Ut(:,p,q,l)+Flux(:,p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(1) ! slave side, SideID=q,jSide=p
      DO l=0,Nloc; DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,l)=Ut(:,p,q,l)-Flux(:,q,p)*L_hatPlus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO l=0,Nloc; DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,l)=Ut(:,p,q,l)-Flux(:,Nloc-p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO l=0,Nloc; DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,l)=Ut(:,p,q,l)-Flux(:,Nloc-q,Nloc-p)*L_hatPlus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO l=0,Nloc; DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,l)=Ut(:,p,q,l)-Flux(:,p,Nloc-q)*L_hatPlus(l)
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
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,0,p,q)=Ut(:,0,p,q)+Flux(:,q,p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,0,p,q)=Ut(:,0,p,q)-Flux(:,p,q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,0,p,q)=Ut(:,0,p,q)-Flux(:,Nloc-q,p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,0,p,q)=Ut(:,0,p,q)-Flux(:,Nloc-p,Nloc-q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,0,p,q)=Ut(:,0,p,q)-Flux(:,q,Nloc-p)*L_hatMinus0
      END DO; END DO ! p,q
    END SELECT

  ! switch to right hand system for ETA_PLUS direction
!===================================================================================================================================
  CASE(ETA_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,0,q)=Ut(:,p,0,q)+Flux(:,p,q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,0,q)=Ut(:,p,0,q)-Flux(:,q,p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,0,q)=Ut(:,p,0,q)-Flux(:,Nloc-p,q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,0,q)=Ut(:,p,0,q)-Flux(:,Nloc-q,Nloc-p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,0,q)=Ut(:,p,0,q)-Flux(:,p,Nloc-q)*L_hatMinus0
      END DO; END DO ! p,q
    END SELECT

  ! switch to right hand system for ZETA_MINUS direction
!===================================================================================================================================
  CASE(ZETA_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,0)=Ut(:,p,q,0)+Flux(:,q,p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,0)=Ut(:,p,q,0)-Flux(:,p,q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,0)=Ut(:,p,q,0)-Flux(:,Nloc-q,p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,0)=Ut(:,p,q,0)-Flux(:,Nloc-p,Nloc-q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,0)=Ut(:,p,q,0)-Flux(:,q,Nloc-p)*L_hatMinus0
      END DO; END DO ! p,q
    END SELECT

!===================================================================================================================================
  CASE(XI_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,Nloc,p,q)=Ut(:,Nloc,p,q)+Flux(:,p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,Nloc,p,q)=Ut(:,Nloc,p,q)-Flux(:,q,p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,Nloc,p,q)=Ut(:,Nloc,p,q)-Flux(:,Nloc-p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,Nloc,p,q)=Ut(:,Nloc,p,q)-Flux(:,Nloc-q,Nloc-p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,Nloc,p,q)=Ut(:,Nloc,p,q)-Flux(:,p,Nloc-q)*L_hatPlusN
      END DO; END DO ! p,q
    END SELECT

  ! switch to right hand system for ETA_PLUS direction
!===================================================================================================================================
  CASE(ETA_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,Nloc,q)=Ut(:,p,Nloc,q)+Flux(:,Nloc-p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,Nloc,q)=Ut(:,p,Nloc,q)-Flux(:,q,Nloc-p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,Nloc,q)=Ut(:,p,Nloc,q)-Flux(:,p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,Nloc,q)=Ut(:,p,Nloc,q)-Flux(:,Nloc-q,p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,Nloc,q)=Ut(:,p,Nloc,q)-Flux(:,Nloc-p,Nloc-q)*L_hatPlusN
      END DO; END DO ! p,q
    END SELECT

  ! switch to right hand system for ZETA_MINUS direction
!===================================================================================================================================
  CASE(ZETA_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,Nloc)=Ut(:,p,q,Nloc)+Flux(:,p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,Nloc)=Ut(:,p,q,Nloc)-Flux(:,q,p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,Nloc)=Ut(:,p,q,Nloc)-Flux(:,Nloc-p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,Nloc)=Ut(:,p,q,Nloc)-Flux(:,Nloc-q,Nloc-p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,Nloc)=Ut(:,p,q,Nloc)-Flux(:,p,Nloc-q)*L_hatPlusN
      END DO; END DO ! p,q
    END SELECT
  END SELECT !locSideID
#endif
END SUBROUTINE CalcSurfInt


SUBROUTINE CalcSurfIntPML(Ut,L_hatMinus,L_hatPlus,Nloc,Flux,flip,ElemID,locSideID)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_PML_Vars ,ONLY: PMLnVar,ElemToPML,U2t
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT) :: Ut(1:PP_nVar,0:Nloc,0:Nloc,0:Nloc)
REAL,INTENT(IN)    :: L_hatMinus(0:Nloc)
REAL,INTENT(IN)    :: L_hatPlus(0:Nloc)
INTEGER,INTENT(IN) :: Nloc
REAL,INTENT(IN)    :: Flux(1:PP_nVar+PMLnVar,0:Nloc,0:Nloc)
INTEGER,INTENT(IN) :: flip,ElemID,locSideID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q,l
#if (PP_NodeType>1)
REAL            ::L_HatMinus0,L_HatPlusN
#endif
!===================================================================================================================================
#if (PP_NodeType>1)
L_hatMinus0 = L_hatMinus(0)
L_hatPlusN  = L_hatPlus(Nloc)
#endif

#if (PP_NodeType==1)
  SELECT CASE(locSideID)
!===================================================================================================================================
  CASE(XI_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO q=0,Nloc; DO p=0,Nloc; DO l=0,Nloc
        Ut( :,l,p,q)           =Ut( :,l,p,q)           +Flux(1: 8,q,p)*L_hatMinus(l)
        U2t(:,l,p,q,ElemToPML(ElemID))=U2t(:,l,p,q,ElemToPML(ElemID))+Flux(9:32,q,p)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,Nloc; DO p=0,Nloc; DO l=0,Nloc
        Ut( :,l,p,q)           =Ut( :,l,p,q)           -Flux(1: 8,p,q)*L_hatMinus(l)
        U2t(:,l,p,q,ElemToPML(ElemID))=U2t(:,l,p,q,ElemToPML(ElemID))-Flux(9:32,p,q)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,Nloc; DO p=0,Nloc; DO l=0,Nloc
        Ut( :,l,p,q)           =Ut( :,l,p,q)           -Flux(1: 8,Nloc-q,p)*L_hatMinus(l)
        U2t(:,l,p,q,ElemToPML(ElemID))=U2t(:,l,p,q,ElemToPML(ElemID))-Flux(9:32,Nloc-q,p)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,Nloc; DO p=0,Nloc; DO l=0,Nloc
        Ut( :,l,p,q)           =Ut( :,l,p,q)           -Flux(1: 8,Nloc-p,Nloc-q)*L_hatMinus(l)
        U2t(:,l,p,q,ElemToPML(ElemID))=U2t(:,l,p,q,ElemToPML(ElemID))-Flux(9:32,Nloc-p,Nloc-q)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,Nloc; DO p=0,Nloc; DO l=0,Nloc
        Ut( :,l,p,q)           =Ut( :,l,p,q)           -Flux(1: 8,q,Nloc-p)*L_hatMinus(l)
        U2t(:,l,p,q,ElemToPML(ElemID))=U2t(:,l,p,q,ElemToPML(ElemID))-Flux(9:32,q,Nloc-p)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    END SELECT
!===================================================================================================================================
  CASE(ETA_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO q=0,Nloc; DO l=0,Nloc; DO p=0,Nloc
        Ut( :,p,l,q)           =Ut( :,p,l,q)           +Flux(1: 8,p,q)*L_hatMinus(l)
        U2t(:,p,l,q,ElemToPML(ElemID))=U2t(:,p,l,q,ElemToPML(ElemID))+Flux(9:32,p,q)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,Nloc; DO l=0,Nloc; DO p=0,Nloc
        Ut( :,p,l,q)           =Ut( :,p,l,q)           -Flux(1: 8,q,p)*L_hatMinus(l)
        U2t(:,p,l,q,ElemToPML(ElemID))=U2t(:,p,l,q,ElemToPML(ElemID))-Flux(9:32,q,p)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,Nloc; DO l=0,Nloc; DO p=0,Nloc
        Ut( :,p,l,q)           =Ut( :,p,l,q)           -Flux(1: 8,Nloc-p,q)*L_hatMinus(l)
        U2t(:,p,l,q,ElemToPML(ElemID))=U2t(:,p,l,q,ElemToPML(ElemID))-Flux(9:32,Nloc-p,q)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,Nloc; DO l=0,Nloc; DO p=0,Nloc
        Ut( :,p,l,q)           =Ut( :,p,l,q)           -Flux(1: 8,Nloc-q,Nloc-p)*L_hatMinus(l)
        U2t(:,p,l,q,ElemToPML(ElemID))=U2t(:,p,l,q,ElemToPML(ElemID))-Flux(9:32,Nloc-q,Nloc-p)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,Nloc; DO l=0,Nloc; DO p=0,Nloc
        Ut( :,p,l,q)           =Ut( :,p,l,q)           -Flux(1: 8,p,Nloc-q)*L_hatMinus(l)
        U2t(:,p,l,q,ElemToPML(ElemID))=U2t(:,p,l,q,ElemToPML(ElemID))-Flux(9:32,p,Nloc-q)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    END SELECT
!===================================================================================================================================
  CASE(ZETA_MINUS)
!===================================================================================================================================

    SELECT CASE(flip)
    CASE(0) ! master side
      DO l=0,Nloc; DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,l)           = Ut( :,p,q,l)           +Flux(1: 8,q,p)*L_hatMinus(l)
        U2t(:,p,q,l,ElemToPML(ElemID))=U2t(:,p,q,l,ElemToPML(ElemID))+Flux(9:32,q,p)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(1) ! slave side, SideID=q,jSide=p
      DO l=0,Nloc; DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,l)           = Ut( :,p,q,l)           -Flux(1: 8,p,q)*L_hatMinus(l)
        U2t(:,p,q,l,ElemToPML(ElemID))=U2t(:,p,q,l,ElemToPML(ElemID))-Flux(9:32,p,q)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO l=0,Nloc; DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,l)           = Ut( :,p,q,l)           -Flux(1: 8,Nloc-q,p)*L_hatMinus(l)
        U2t(:,p,q,l,ElemToPML(ElemID))=U2t(:,p,q,l,ElemToPML(ElemID))-Flux(9:32,Nloc-q,p)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO l=0,Nloc; DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,l)           = Ut( :,p,q,l)           -Flux(1: 8,Nloc-p,Nloc-q)*L_hatMinus(l)
        U2t(:,p,q,l,ElemToPML(ElemID))=U2t(:,p,q,l,ElemToPML(ElemID))-Flux(9:32,Nloc-p,Nloc-q)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO l=0,Nloc; DO q=0,Nloc; DO p=0,Nloc
        Ut(:,p,q,l)           = Ut( :,p,q,l)           -Flux(1: 8,q,Nloc-p)*L_hatMinus(l)
        U2t(:,p,q,l,ElemToPML(ElemID))=U2t(:,p,q,l,ElemToPML(ElemID))-Flux(9:32,q,Nloc-p)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    END SELECT
!===================================================================================================================================
  CASE(XI_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO q=0,Nloc; DO p=0,Nloc; DO l=0,Nloc
        Ut( :,l,p,q)           =Ut( :,l,p,q)           +Flux(1: 8,p,q)*L_hatPlus(l)
        U2t(:,l,p,q,ElemToPML(ElemID))=U2t(:,l,p,q,ElemToPML(ElemID))+Flux(9:32,p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,Nloc; DO p=0,Nloc; DO l=0,Nloc
        Ut( :,l,p,q)           =Ut( :,l,p,q)           -Flux(1: 8,q,p)*L_hatPlus(l)
        U2t(:,l,p,q,ElemToPML(ElemID))=U2t(:,l,p,q,ElemToPML(ElemID))-Flux(9:32,q,p)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,Nloc; DO p=0,Nloc; DO l=0,Nloc
        Ut( :,l,p,q)           =Ut( :,l,p,q)           -Flux(1: 8,Nloc-p,q)*L_hatPlus(l)
        U2t(:,l,p,q,ElemToPML(ElemID))=U2t(:,l,p,q,ElemToPML(ElemID))-Flux(9:32,Nloc-p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,Nloc; DO p=0,Nloc; DO l=0,Nloc
        Ut( :,l,p,q)           =Ut( :,l,p,q)           -Flux(1: 8,Nloc-q,Nloc-p)*L_hatPlus(l)
        U2t(:,l,p,q,ElemToPML(ElemID))=U2t(:,l,p,q,ElemToPML(ElemID))-Flux(9:32,Nloc-q,Nloc-p)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,Nloc; DO p=0,Nloc; DO l=0,Nloc
        Ut( :,l,p,q)           =Ut( :,l,p,q)           -Flux(1: 8,p,Nloc-q)*L_hatPlus(l)
        U2t(:,l,p,q,ElemToPML(ElemID))=U2t(:,l,p,q,ElemToPML(ElemID))-Flux(9:32,p,Nloc-q)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    END SELECT
!===================================================================================================================================
  CASE(ETA_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO q=0,Nloc; DO l=0,Nloc; DO p=0,Nloc
        Ut( :,p,l,q)           =Ut( :,p,l,q)           +Flux(1: 8,Nloc-p,q)*L_hatPlus(l)
        U2t(:,p,l,q,ElemToPML(ElemID))=U2t(:,p,l,q,ElemToPML(ElemID))+Flux(9:32,Nloc-p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,Nloc; DO l=0,Nloc; DO p=0,Nloc
        Ut( :,p,l,q)           =Ut( :,p,l,q)           -Flux(1: 8,q,Nloc-p)*L_hatPlus(l)
        U2t(:,p,l,q,ElemToPML(ElemID))=U2t(:,p,l,q,ElemToPML(ElemID))-Flux(9:32,q,Nloc-p)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,Nloc; DO l=0,Nloc; DO p=0,Nloc
        Ut( :,p,l,q)           =Ut( :,p,l,q)           -Flux(1: 8,p,q)*L_hatPlus(l)
        U2t(:,p,l,q,ElemToPML(ElemID))=U2t(:,p,l,q,ElemToPML(ElemID))-Flux(9:32,p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,Nloc; DO l=0,Nloc; DO p=0,Nloc
        Ut( :,p,l,q)           =Ut( :,p,l,q)           -Flux(1: 8,Nloc-q,p)*L_hatPlus(l)
        U2t(:,p,l,q,ElemToPML(ElemID))=U2t(:,p,l,q,ElemToPML(ElemID))-Flux(9:32,Nloc-q,p)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,Nloc; DO l=0,Nloc; DO p=0,Nloc
        Ut( :,p,l,q)           =Ut( :,p,l,q)           -Flux(1: 8,Nloc-p,Nloc-q)*L_hatPlus(l)
        U2t(:,p,l,q,ElemToPML(ElemID))=U2t(:,p,l,q,ElemToPML(ElemID))-Flux(9:32,Nloc-p,Nloc-q)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    END SELECT
!===================================================================================================================================
  CASE(ZETA_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO l=0,Nloc; DO q=0,Nloc; DO p=0,Nloc
        Ut( :,p,q,l)           =Ut( :,p,q,l)           +Flux(1: 8,p,q)*L_hatPlus(l)
        U2t(:,p,q,l,ElemToPML(ElemID))=U2t(:,p,q,l,ElemToPML(ElemID))+Flux(9:32,p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(1) ! slave side, SideID=q,jSide=p
      DO l=0,Nloc; DO q=0,Nloc; DO p=0,Nloc
        Ut( :,p,q,l)           =Ut( :,p,q,l)           -Flux(1: 8,q,p)*L_hatPlus(l)
        U2t(:,p,q,l,ElemToPML(ElemID))=U2t(:,p,q,l,ElemToPML(ElemID))-Flux(9:32,q,p)*L_hatPlus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO l=0,Nloc; DO q=0,Nloc; DO p=0,Nloc
        Ut( :,p,q,l)           =Ut( :,p,q,l)           -Flux(1: 8,Nloc-p,q)*L_hatPlus(l)
        U2t(:,p,q,l,ElemToPML(ElemID))=U2t(:,p,q,l,ElemToPML(ElemID))-Flux(9:32,Nloc-p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO l=0,Nloc; DO q=0,Nloc; DO p=0,Nloc
        Ut( :,p,q,l)           =Ut( :,p,q,l)           -Flux(1: 8,Nloc-q,Nloc-p)*L_hatPlus(l)
        U2t(:,p,q,l,ElemToPML(ElemID))=U2t(:,p,q,l,ElemToPML(ElemID))-Flux(9:32,Nloc-q,Nloc-p)*L_hatPlus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO l=0,Nloc; DO q=0,Nloc; DO p=0,Nloc
        Ut( :,p,q,l)           =Ut( :,p,q,l)           -Flux(1: 8,p,Nloc-q)*L_hatPlus(l)
        U2t(:,p,q,l,ElemToPML(ElemID))=U2t(:,p,q,l,ElemToPML(ElemID))-Flux(9:32,p,Nloc-q)*L_hatPlus(l)
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
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,0,p,q)           =Ut (:,0,p,q,          ElemID )+Flux(1: 8,q,p)*L_hatMinus0
        U2t(:,0,p,q,ElemToPML(ElemID))=U2t(:,0,p,q,ElemToPML(ElemID))+Flux(9:32,q,p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,0,p,q)           =Ut (:,0,p,q)           -Flux(1:8  ,p,q)*L_hatMinus0
        U2t(:,0,p,q,ElemToPML(ElemID))=U2t(:,0,p,q,ElemToPML(ElemID))-Flux(9:32,p,q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,Nloc; DO p=0,Nloc
        Ut( :,0,p,q)           =Ut (:,0,p,q)           -Flux(1:8 ,Nloc-q,p)*L_hatMinus0
        U2t(:,0,p,q,ElemToPML(ElemID))=U2t(:,0,p,q,ElemToPML(ElemID))-Flux(9:32,Nloc-q,p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,0,p,q)           =Ut (:,0,p,q)           -Flux(1:8 ,Nloc-p,Nloc-q)*L_hatMinus0
        U2t(:,0,p,q,ElemToPML(ElemID))=U2t(:,0,p,q,ElemToPML(ElemID))-Flux(9:32,Nloc-p,Nloc-q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,0,p,q)           =Ut (:,0,p,q)           -Flux(1:8 ,q,Nloc-p)*L_hatMinus0
        U2t(:,0,p,q,ElemToPML(ElemID))=U2t(:,0,p,q,ElemToPML(ElemID))-Flux(9:32,q,Nloc-p)*L_hatMinus0
      END DO; END DO ! p,q
    END SELECT

  ! switch to right hand system for ETA_PLUS direction
!===================================================================================================================================
  CASE(ETA_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,p,0,q)           =Ut (:,p,0,q)           +Flux(1:8 ,p,q)*L_hatMinus0
        U2t(:,p,0,q,ElemToPML(ElemID))=U2t(:,p,0,q,ElemToPML(ElemID))+Flux(9:32,p,q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,p,0,q)           =Ut (:,p,0,q)           -Flux(1:8 ,q,p)*L_hatMinus0
        U2t(:,p,0,q,ElemToPML(ElemID))=U2t(:,p,0,q,ElemToPML(ElemID))-Flux(9:32,q,p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,p,0,q)           =Ut (:,p,0,q)           -Flux(1:8 ,Nloc-p,q)*L_hatMinus0
        U2t(:,p,0,q,ElemToPML(ElemID))=U2t(:,p,0,q,ElemToPML(ElemID))-Flux(9:32,Nloc-p,q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,p,0,q)           =Ut (:,p,0,q)           -Flux(1:8 ,Nloc-q,Nloc-p)*L_hatMinus0
        U2t(:,p,0,q,ElemToPML(ElemID))=U2t(:,p,0,q,ElemToPML(ElemID))-Flux(9:32,Nloc-q,Nloc-p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,p,0,q)           =Ut (:,p,0,q)           -Flux(1:8 ,p,Nloc-q)*L_hatMinus0
        U2t(:,p,0,q,ElemToPML(ElemID))=U2t(:,p,0,q,ElemToPML(ElemID))-Flux(9:32,p,Nloc-q)*L_hatMinus0
      END DO; END DO ! p,q
    END SELECT

  ! switch to right hand system for ZETA_MINUS direction
!===================================================================================================================================
  CASE(ZETA_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,p,q,0)           =Ut (:,p,q,0)           +Flux(1:8 ,q,p)*L_hatMinus0
        U2t(:,p,q,0,ElemToPML(ElemID))=U2t(:,p,q,0,ElemToPML(ElemID))+Flux(9:32,q,p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,p,q,0)           =Ut (:,p,q,0)           -Flux(1:8 ,p,q)*L_hatMinus0
        U2t(:,p,q,0,ElemToPML(ElemID))=U2t(:,p,q,0,ElemToPML(ElemID))-Flux(9:32,p,q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,p,q,0)           =Ut (:,p,q,0)           -Flux(1:8 ,Nloc-q,p)*L_hatMinus0
        U2t(:,p,q,0,ElemToPML(ElemID))=U2t(:,p,q,0,ElemToPML(ElemID))-Flux(9:32,Nloc-q,p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,p,q,0)           =Ut (:,p,q,0)           -Flux(1:8 ,Nloc-p,Nloc-q)*L_hatMinus0
        U2t(:,p,q,0,ElemToPML(ElemID))=U2t(:,p,q,0,ElemToPML(ElemID))-Flux(9:32,Nloc-p,Nloc-q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,p,q,0)           =Ut (:,p,q,0)           -Flux(1:8 ,q,Nloc-p)*L_hatMinus0
        U2t(:,p,q,0,ElemToPML(ElemID))=U2t(:,p,q,0,ElemToPML(ElemID))-Flux(9:32,q,Nloc-p)*L_hatMinus0
      END DO; END DO ! p,q
    END SELECT

!===================================================================================================================================
  CASE(XI_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,Nloc,p,q)           =Ut (:,Nloc,p,q)           +Flux(1:8 ,p,q)*L_hatPlusN
        U2t(:,Nloc,p,q,ElemToPML(ElemID))=U2t(:,Nloc,p,q,ElemToPML(ElemID))+Flux(9:32,p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,Nloc,p,q)           =Ut (:,Nloc,p,q)           -Flux(1:8 ,q,p)*L_hatPlusN
        U2t(:,Nloc,p,q,ElemToPML(ElemID))=U2t(:,Nloc,p,q,ElemToPML(ElemID))-Flux(9:32,q,p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,Nloc,p,q)           =Ut (:,Nloc,p,q)           -Flux(1:8 ,Nloc-p,q)*L_hatPlusN
        U2t(:,Nloc,p,q,ElemToPML(ElemID))=U2t(:,Nloc,p,q,ElemToPML(ElemID))-Flux(9:32,Nloc-p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,Nloc,p,q)           =Ut (:,Nloc,p,q)           -Flux(1:8 ,Nloc-q,Nloc-p)*L_hatPlusN
        U2t(:,Nloc,p,q,ElemToPML(ElemID))=U2t(:,Nloc,p,q,ElemToPML(ElemID))-Flux(9:32,Nloc-q,Nloc-p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,Nloc,p,q)           =Ut (:,Nloc,p,q)           -Flux(1:8 ,p,Nloc-q)*L_hatPlusN
        U2t(:,Nloc,p,q,ElemToPML(ElemID))=U2t(:,Nloc,p,q,ElemToPML(ElemID))-Flux(9:32,p,Nloc-q)*L_hatPlusN
      END DO; END DO ! p,q
    END SELECT

  ! switch to right hand system for ETA_PLUS direction
!===================================================================================================================================
  CASE(ETA_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,p,Nloc,q)           =Ut (:,p,Nloc,q)           +Flux(1:8 ,Nloc-p,q)*L_hatPlusN
        U2t(:,p,Nloc,q,ElemToPML(ElemID))=U2t(:,p,Nloc,q,ElemToPML(ElemID))+Flux(9:32,Nloc-p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,p,Nloc,q)           =Ut (:,p,Nloc,q)           -Flux(1:8 ,q,Nloc-p)*L_hatPlusN
        U2t(:,p,Nloc,q,ElemToPML(ElemID))=U2t(:,p,Nloc,q,ElemToPML(ElemID))-Flux(9:32,q,Nloc-p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,p,Nloc,q)           =Ut (:,p,Nloc,q)           -Flux(1:8 ,p,q)*L_hatPlusN
        U2t(:,p,Nloc,q,ElemToPML(ElemID))=U2t(:,p,Nloc,q,ElemToPML(ElemID))-Flux(9:32,p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,p,Nloc,q)           =Ut (:,p,Nloc,q)           -Flux(1:8 ,Nloc-q,p)*L_hatPlusN
        U2t(:,p,Nloc,q,ElemToPML(ElemID))=U2t(:,p,Nloc,q,ElemToPML(ElemID))-Flux(9:32,Nloc-q,p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,p,Nloc,q)           =Ut (:,p,Nloc,q)           -Flux(1:8 ,Nloc-p,Nloc-q)*L_hatPlusN
        U2t(:,p,Nloc,q,ElemToPML(ElemID))=U2t(:,p,Nloc,q,ElemToPML(ElemID))-Flux(9:32,Nloc-p,Nloc-q)*L_hatPlusN
      END DO; END DO ! p,q
    END SELECT

  ! switch to right hand system for ZETA_MINUS direction
!===================================================================================================================================
  CASE(ZETA_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,p,q,Nloc)           =Ut (:,p,q,Nloc)           +Flux(1:8 ,p,q)*L_hatPlusN
        U2t(:,p,q,Nloc,ElemToPML(ElemID))=U2t(:,p,q,Nloc,ElemToPML(ElemID))+Flux(9:32,p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,p,q,Nloc)           =Ut (:,p,q,Nloc)           -Flux(1:8 ,q,p)*L_hatPlusN
        U2t(:,p,q,Nloc,ElemToPML(ElemID))=U2t(:,p,q,Nloc,ElemToPML(ElemID))-Flux(9:32,q,p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,p,q,Nloc)           =Ut (:,p,q,Nloc)           -Flux(1:8 ,Nloc-p,q)*L_hatPlusN
        U2t(:,p,q,Nloc,ElemToPML(ElemID))=U2t(:,p,q,Nloc,ElemToPML(ElemID))-Flux(9:32,Nloc-p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,p,q,Nloc)           =Ut (:,p,q,Nloc)           -Flux(1:8 ,Nloc-q,Nloc-p)*L_hatPlusN
        U2t(:,p,q,Nloc,ElemToPML(ElemID))=U2t(:,p,q,Nloc,ElemToPML(ElemID))-Flux(9:32,Nloc-q,Nloc-p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,Nloc; DO p=0,Nloc
        Ut (:,p,q,Nloc)           =Ut (:,p,q,Nloc)           -Flux(1:8 ,p,Nloc-q)*L_hatPlusN
        U2t(:,p,q,Nloc,ElemToPML(ElemID))=U2t(:,p,q,Nloc,ElemToPML(ElemID))-Flux(9:32,p,Nloc-q)*L_hatPlusN
      END DO; END DO ! p,q
    END SELECT
  END SELECT !locSideID
#endif
END SUBROUTINE CalcSurfIntPML

END MODULE MOD_SurfInt
