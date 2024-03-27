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

MODULE MOD_ProlongToFace
!===================================================================================================================================
! Contains routines to interpolate the interior solution to the boundary
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
INTERFACE ProlongToFace
  MODULE PROCEDURE ProlongToFace_SideBased
END INTERFACE

INTERFACE ProlongToFace_TypeBased
  MODULE PROCEDURE ProlongToFace_TypeBased
END INTERFACE

! NO interface because of possibility to MAP arrays, dimension reduction, increase
!INTERFACE ProlongToFace_Elementlocal
!  MODULE PROCEDURE ProlongToFace_Elementlocal
!END INTERFACE

PUBLIC::ProlongToFace
PUBLIC::ProlongToFace_Elementlocal
PUBLIC::ProlongToFace_TypeBased, ProlongToFace_Side
!===================================================================================================================================

CONTAINS


SUBROUTINE ProlongToFace_SideBased(Uvol,Uface_master,Uface_slave,doMPISides)
!===================================================================================================================================
! Interpolates the interior volume data (stored at the Gauss or Gauss-Lobatto points) to the surface
! integration points, using fast 1D Interpolation and store in global side structure
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Interpolation_Vars, ONLY: N_Inter
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: SideToElem
USE MOD_Mesh_Vars,          ONLY: firstBCSide,firstInnerSide
USE MOD_Mesh_Vars,          ONLY: firstMPISide_YOUR,lastMPISide_YOUR,lastMPISide_MINE,nSides,firstMortarMPISide,lastMortarMPISide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)              :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides +InnerSides +MPISides MINE
REAL,INTENT(IN)                 :: Uvol(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)              :: Uface_master(PP_nVar,0:PP_N,0:PP_N,1:nSides)
REAL,INTENT(INOUT)              :: Uface_slave(PP_nVar,0:PP_N,0:PP_N,1:nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: l,p,q,ElemID,SideID,flip,LocSideID,firstSideID,lastSideID
REAL                            :: Uface(PP_nVar,0:PP_N,0:PP_N)
!===================================================================================================================================
IF(doMPISides)THEN
  ! only YOUR MPI Sides are filled
  firstSideID = firstMPISide_YOUR
  lastSideID  = lastMPISide_YOUR
ELSE
  ! (Mortar-)InnerSides and MINE MPISides are filled
  firstSideID = firstInnerSide
  lastSideID  = lastMPISide_MINE
END IF

ASSOCIATE( L_Minus => N_Inter(PP_N)%L_Minus ,&
           L_Plus  => N_Inter(PP_N)%L_Plus  )

DO SideID=firstSideID,lastSideID
  ! neighbor side !ElemID,locSideID and flip =-1 if not existing
  ElemID     = SideToElem(S2E_NB_ELEM_ID,SideID)
  locSideID  = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  flip         = SideToElem(S2E_FLIP,SideID)
#if (PP_NodeType==1) /* for Gauss-points*/
  SELECT CASE(locSideID)
  CASE(XI_MINUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface(:,q,p)=Uvol(:,0,p,q,ElemID)*L_Minus(0)
        DO l=1,PP_N
          ! switch to right hand system
          Uface(:,q,p)=Uface(:,q,p)+Uvol(:,l,p,q,ElemID)*L_Minus(l)
        END DO ! l
      END DO ! p
    END DO ! q
  CASE(ETA_MINUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface(:,p,q)=Uvol(:,p,0,q,ElemID)*L_Minus(0)
        DO l=1,PP_N
          Uface(:,p,q)=Uface(:,p,q)+Uvol(:,p,l,q,ElemID)*L_Minus(l)
        END DO ! l
      END DO ! p
    END DO ! q
  CASE(ZETA_MINUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface(:,q,p)=Uvol(:,p,q,0,ElemID)*L_Minus(0)
        DO l=1,PP_N
          ! switch to right hand system
          Uface(:,q,p)=Uface(:,q,p)+Uvol(:,p,q,l,ElemID)*L_Minus(l)
        END DO ! l
      END DO ! p
    END DO ! q
  CASE(XI_PLUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface(:,p,q)=Uvol(:,0,p,q,ElemID)*L_Plus(0)
        DO l=1,PP_N
          Uface(:,p,q)=Uface(:,p,q)+Uvol(:,l,p,q,ElemID)*L_Plus(l)
        END DO ! l
      END DO ! p
    END DO ! q
  CASE(ETA_PLUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface(:,PP_N-p,q)=Uvol(:,p,0,q,ElemID)*L_Plus(0)
        DO l=1,PP_N
          ! switch to right hand system
          Uface(:,PP_N-p,q)=Uface(:,PP_N-p,q)+Uvol(:,p,l,q,ElemID)*L_Plus(l)
        END DO ! l
      END DO ! p
    END DO ! q
  CASE(ZETA_PLUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface(:,p,q)=Uvol(:,p,q,0,ElemID)*L_Plus(0)
        DO l=1,PP_N
          Uface(:,p,q)=Uface(:,p,q)+Uvol(:,p,q,l,ElemID)*L_Plus(l)
        END DO ! l
      END DO ! p
    END DO ! q
  END SELECT
#else /* for Gauss-Lobatto-points*/
  SELECT CASE(locSideID)
  CASE(XI_MINUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface(:,q,p)=Uvol(:,0,p,q,ElemID)
      END DO ! p
    END DO ! q
  CASE(ETA_MINUS)
    Uface(:,:,:)=Uvol(:,:,0,:,ElemID)
  CASE(ZETA_MINUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface(:,q,p)=Uvol(:,p,q,0,ElemID)
      END DO ! p
    END DO ! q
  CASE(XI_PLUS)
    Uface(:,:,:)=Uvol(:,PP_N,:,:,ElemID)
  CASE(ETA_PLUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface(:,PP_N-p,q)=Uvol(:,p,PP_N,q,ElemID)
      END DO ! p
    END DO ! q
  CASE(ZETA_PLUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface(:,p,q)=Uvol(:,p,q,PP_N,ElemID)
      END DO ! p
    END DO ! q
  END SELECT
#endif
  SELECT CASE(flip)
    !CASE(0) ! master side
    !  Uface_Minus(:,:,:,SideID)=Uface(:,:,:)
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,PP_N
        DO p=0,PP_N
          Uface_slave(:,p,q,SideID)=Uface(:,q,p)
        END DO ! p
      END DO ! q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,PP_N
        DO p=0,PP_N
          Uface_slave(:,p,q,SideID)=Uface(:,PP_N-p,q)
        END DO ! p
      END DO ! q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,PP_N
        DO p=0,PP_N
          Uface_slave(:,p,q,SideID)=Uface(:,PP_N-q,PP_N-p)
        END DO ! p
      END DO ! q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,PP_N
        DO p=0,PP_N
          Uface_slave(:,p,q,SideID)=Uface(:,p,PP_N-q)
        END DO ! p
      END DO ! q
  END SELECT
END DO !SideID

! Second process Minus/Master sides, U_Minus is always MINE
! master side, flip=0
IF(doMPISides)THEN
  ! only MPI mortars
  firstSideID = firstMortarMPISide
   lastSideID =  lastMortarMPISide
ELSE
  ! BCSides, (Mortar-)InnerSides and MINE MPISides are filled
  firstSideID = firstBCSide
   lastSideID =  lastMPISide_MINE
END IF
DO SideID=firstSideID,lastSideID
  ElemID    = SideToElem(S2E_ELEM_ID,SideID)
  locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
#if (PP_NodeType==1) /* for Gauss-points*/
  SELECT CASE(locSideID)
  CASE(XI_MINUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface_master(:,q,p,SideID)=Uvol(:,0,p,q,ElemID)*L_Minus(0)
        DO l=1,PP_N
          ! switch to right hand system
          Uface_master(:,q,p,SideID)=Uface_master(:,q,p,SideID)+Uvol(:,l,p,q,ElemID)*L_Minus(l)
        END DO ! l
      END DO ! p
    END DO ! q
  CASE(ETA_MINUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface_master(:,p,q,SideID)=Uvol(:,p,0,q,ElemID)*L_Minus(0)
        DO l=1,PP_N
          Uface_master(:,p,q,SideID)=Uface_master(:,p,q,SideID)+Uvol(:,p,l,q,ElemID)*L_Minus(l)
        END DO ! l
      END DO ! p
    END DO ! q
  CASE(ZETA_MINUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface_master(:,q,p,SideID)=Uvol(:,p,q,0,ElemID)*L_Minus(0)
        DO l=1,PP_N
          ! switch to right hand system
          Uface_master(:,q,p,SideID)=Uface_master(:,q,p,SideID)+Uvol(:,p,q,l,ElemID)*L_Minus(l)
        END DO ! l
      END DO ! p
    END DO ! q
  CASE(XI_PLUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface_master(:,p,q,SideID)=Uvol(:,0,p,q,ElemID)*L_Plus(0)
        DO l=1,PP_N
          Uface_master(:,p,q,SideID)=Uface_master(:,p,q,SideID)+Uvol(:,l,p,q,ElemID)*L_Plus(l)
        END DO ! l
      END DO ! p
    END DO ! q
  CASE(ETA_PLUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface_master(:,PP_N-p,q,SideID)=Uvol(:,p,0,q,ElemID)*L_Plus(0)
        DO l=1,PP_N
          ! switch to right hand system
          Uface_master(:,PP_N-p,q,SideID)=Uface_master(:,PP_N-p,q,SideID)+Uvol(:,p,l,q,ElemID)*L_Plus(l)
        END DO ! l
      END DO ! p
    END DO ! q
  CASE(ZETA_PLUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface_master(:,p,q,SideID)=Uvol(:,p,q,0,ElemID)*L_Plus(0)
        DO l=1,PP_N
          Uface_master(:,p,q,SideID)=Uface_master(:,p,q,SideID)+Uvol(:,p,q,l,ElemID)*L_Plus(l)
        END DO ! l
      END DO ! p
    END DO ! q
  END SELECT
#else /* for Gauss-Lobatto-points*/
  SELECT CASE(locSideID)
  CASE(XI_MINUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface_master(:,q,p,SideID)=Uvol(:,0,p,q,ElemID)
      END DO ! p
    END DO ! q
  CASE(ETA_MINUS)
    Uface_master(:,:,:,SideID)=Uvol(:,:,0,:,ElemID)
  CASE(ZETA_MINUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface_master(:,q,p,SideID)=Uvol(:,p,q,0,ElemID)
      END DO ! p
    END DO ! q
  CASE(XI_PLUS)
    Uface_master(:,:,:,SideID)=Uvol(:,PP_N,:,:,ElemID)
  CASE(ETA_PLUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface_master(:,PP_N-p,q,SideID)=Uvol(:,p,PP_N,q,ElemID)
      END DO ! p
    END DO ! q
  CASE(ZETA_PLUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface_master(:,p,q,SideID)=Uvol(:,p,q,PP_N,ElemID)
      END DO ! p
    END DO ! q
  END SELECT
#endif
END DO !SideID

END ASSOCIATE

END SUBROUTINE ProlongToFace_SideBased



SUBROUTINE ProlongToFace_TypeBased(doMPISides)
!===================================================================================================================================
! Interpolates the interior volume data (stored at the Gauss or Gauss-Lobatto points) to the surface
! integration points, using fast 1D Interpolation and store in global side structure
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Interpolation_Vars ,ONLY: N_Inter
USE MOD_DG_Vars            ,ONLY: U_N,DG_Elems_slave,DG_Elems_master,U_Surf_N, N_DG_Mapping
USE MOD_Mesh_Vars          ,ONLY: SideToElem, offSetElem
USE MOD_Mesh_Vars          ,ONLY: firstBCSide,firstInnerSide
USE MOD_Mesh_Vars          ,ONLY: firstMPISide_YOUR,lastMPISide_YOUR,lastMPISide_MINE,nSides,firstMortarMPISide,lastMortarMPISide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides +InnerSides +MPISides MINE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: l,p,q,nbElemID,ElemID,SideID,flip,LocSideID,firstSideID,lastSideID,Nloc
REAL,ALLOCATABLE   :: Uface(:,:,:)
!===================================================================================================================================
IF(doMPISides)THEN
  ! only YOUR MPI Sides are filled
  firstSideID = firstMPISide_YOUR
  lastSideID  = lastMPISide_YOUR
ELSE
  ! (Mortar-)InnerSides and MINE MPISides are filled
  firstSideID = firstInnerSide
  lastSideID  = lastMPISide_MINE
END IF

DO SideID=firstSideID,lastSideID
  ! neighbor side !nbElemID,locSideID and flip =-1 if not existing
  nbElemID   = SideToElem(S2E_NB_ELEM_ID,SideID)
  IF(nbElemID.LE.0) CYCLE
  locSideID  = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  flip       = SideToElem(S2E_FLIP,SideID)
  Nloc = N_DG_Mapping(2,nbElemID+offSetElem)
  !CALL ProlongToFace_Side(PP_nVar, Nloc, locSideID, flip, nbElemID,  U_N(nbElemID)%U, U_Surf_N(SideID)%U_slave)
  CALL ProlongToFace_Side(PP_nVar, Nloc, locSideID, flip, U_N(nbElemID)%U, U_Surf_N(SideID)%U_slave)
  
END DO !SideID

! Second process Minus/Master sides, U_Minus is always MINE
! master side, flip=0
IF(doMPISides)THEN
  ! only MPI mortars
  firstSideID = firstMortarMPISide
   lastSideID =  lastMortarMPISide
ELSE
  ! BCSides, (Mortar-)InnerSides and MINE MPISides are filled
  firstSideID = firstBCSide
   lastSideID =  lastMPISide_MINE
END IF
DO SideID=firstSideID,lastSideID
  ElemID    = SideToElem(S2E_ELEM_ID,SideID)
  IF(ElemID.LE.0) CYCLE
  locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
  Nloc = N_DG_Mapping(2,ElemID+offSetElem)

  !CALL ProlongToFace_Side(PP_nVar, Nloc, locSideID, 0, ElemID,  U_N(ElemID)%U, U_Surf_N(SideID)%U_master)
  CALL ProlongToFace_Side(PP_nVar, Nloc, locSideID, 0, U_N(ElemID)%U, U_Surf_N(SideID)%U_master)
END DO !SideID

END SUBROUTINE ProlongToFace_TypeBased


!SUBROUTINE ProlongToFace_Side(Nvar, Nloc, locSideID, flip, ElemID,  U, USide)
SUBROUTINE ProlongToFace_Side(Nvar, Nloc, locSideID, flip, U, USide)
!===================================================================================================================================
! Interpolates the interior volume data (stored at the Gauss or Gauss-Lobatto points) to the surface
! integration points, using fast 1D Interpolation and store in global side structure
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Interpolation_Vars ,ONLY: N_Inter
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN):: locSideID, flip, Nloc, nVar !,ElemID
REAL, INTENT(INOUT):: U(1:nVar,0:Nloc,0:Nloc,0:Nloc)
REAL, INTENT(INOUT):: USide(1:nvar,0:Nloc,0:Nloc)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: l,p,q
REAL               :: Uface(Nvar,0:Nloc,0:Nloc)
!===================================================================================================================================
  ! neighbor side !ElemID,locSideID and flip =-1 if not existing

#if (PP_NodeType==1) /* for Gauss-points*/
  SELECT CASE(locSideID)
  CASE(XI_MINUS)
    DO q=0,Nloc
      DO p=0,Nloc
        Uface(:,q,p)=U(:,0,p,q)*N_Inter(Nloc)%L_Minus(0)
        DO l=1,Nloc
          ! switch to right hand system
          Uface(:,q,p)=Uface(:,q,p)+U(:,l,p,q)*N_Inter(Nloc)%L_Minus(l)
        END DO ! l
      END DO ! p
    END DO ! q
  CASE(ETA_MINUS)
    DO q=0,Nloc
      DO p=0,Nloc
        Uface(:,p,q)=U(:,p,0,q)*N_Inter(Nloc)%L_Minus(0)
        DO l=1,Nloc
          Uface(:,p,q)=Uface(:,p,q)+U(:,p,l,q)*N_Inter(Nloc)%L_Minus(l)
        END DO ! l
      END DO ! p
    END DO ! q
  CASE(ZETA_MINUS)
    DO q=0,Nloc
      DO p=0,Nloc
        Uface(:,q,p)=U(:,p,q,0)*N_Inter(Nloc)%L_Minus(0)
        DO l=1,Nloc
          ! switch to right hand system
          Uface(:,q,p)=Uface(:,q,p)+U(:,p,q,l)*N_Inter(Nloc)%L_Minus(l)
        END DO ! l
      END DO ! p
    END DO ! q
  CASE(XI_PLUS)
    DO q=0,Nloc
      DO p=0,Nloc
        Uface(:,p,q)=U(:,0,p,q)*N_Inter(Nloc)%L_Plus(0)
        DO l=1,Nloc
          Uface(:,p,q)=Uface(:,p,q)+U(:,l,p,q)*N_Inter(Nloc)%L_Plus(l)
        END DO ! l
      END DO ! p
    END DO ! q
  CASE(ETA_PLUS)
    DO q=0,Nloc
      DO p=0,Nloc
        Uface(:,Nloc-p,q)=U(:,p,0,q)*N_Inter(Nloc)%L_Plus(0)
        DO l=1,Nloc
          ! switch to right hand system
          Uface(:,Nloc-p,q)=Uface(:,Nloc-p,q)+U(:,p,l,q)*N_Inter(Nloc)%L_Plus(l)
        END DO ! l
      END DO ! p
    END DO ! q
  CASE(ZETA_PLUS)
    DO q=0,Nloc
      DO p=0,Nloc
        Uface(:,p,q)=U(:,p,q,0)*N_Inter(Nloc)%L_Plus(0)
        DO l=1,Nloc
          Uface(:,p,q)=Uface(:,p,q)+U(:,p,q,l)*N_Inter(Nloc)%L_Plus(l)
        END DO ! l
      END DO ! p
    END DO ! q
  END SELECT
#else /* for Gauss-Lobatto-points*/
  SELECT CASE(locSideID)
  CASE(XI_MINUS)
    DO q=0,Nloc
      DO p=0,Nloc
        Uface(:,q,p)=U(:,0,p,q)
      END DO ! p
    END DO ! q
  CASE(ETA_MINUS)
    Uface(:,:,:)=U(:,:,0,:)
  CASE(ZETA_MINUS)
    DO q=0,Nloc
      DO p=0,Nloc
        Uface(:,q,p)=U(:,p,q,0)
      END DO ! p
    END DO ! q
  CASE(XI_PLUS)
    Uface(:,:,:)=U(:,Nloc,:,:)
  CASE(ETA_PLUS)
    DO q=0,Nloc
      DO p=0,Nloc
        Uface(:,Nloc-p,q)=U(:,p,Nloc,q)
      END DO ! p
    END DO ! q
  CASE(ZETA_PLUS)
    DO q=0,Nloc
      DO p=0,Nloc
        Uface(:,p,q)=U(:,p,q,Nloc)
      END DO ! p
    END DO ! q
  END SELECT
#endif
  SELECT CASE(flip)
    CASE(0) ! master side
      USide(:,:,:)=Uface(:,:,:)
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,Nloc
        DO p=0,Nloc
          USide(:,p,q)=Uface(:,q,p)
        END DO ! p
      END DO ! q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,Nloc
        DO p=0,Nloc
          USide(:,p,q)=Uface(:,Nloc-p,q)
        END DO ! p
      END DO ! q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,Nloc
        DO p=0,Nloc
          USide(:,p,q)=Uface(:,Nloc-q,Nloc-p)
        END DO ! p
      END DO ! q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,Nloc
        DO p=0,Nloc
          USide(:,p,q)=Uface(:,p,Nloc-q)
        END DO ! p
      END DO ! q
  END SELECT
END SUBROUTINE ProlongToFace_Side


SUBROUTINE ProlongToFace_Elementlocal(nVar,locSideID,Uvol,Uface, Nloc)
!===================================================================================================================================
! Interpolates the interior volume data (stored at the Gauss or Gauss-Lobatto points) to the surface
! integration points, using fast 1D Interpolation and does NOT rotate into global coordinate system
! Face integration points are still volume-IJK orientated
! Has to be called for each element face separately
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Interpolation_Vars, ONLY: N_Inter
USE MOD_PreProc
USE MOD_Mappings,           ONLY: CGNS_VolToSide_IJK
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)              :: nVar
INTEGER,INTENT(IN)              :: locSideID
REAL,INTENT(IN)                 :: Uvol(1:nVar,0:Nloc,0:Nloc,0:Nloc)
INTEGER,INTENT(IN)               :: Nloc
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Uface(1:nVar,0:Nloc,0:Nloc)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k
INTEGER                         :: pq(1:3)
!===================================================================================================================================

#if (PP_NodeType==1) /* for Gauss-points*/
Uface=0.
DO k=0,Nloc
  DO j=0,Nloc
    DO i=0,Nloc
      pq=CGNS_VolToSide_IJK(i,j,k,locSideID)
      Uface(:,pq(1),pq(2))=Uface(:,pq(1),pq(2))+Uvol(:,i,j,k)*N_Inter(Nloc)%L_PlusMinus(pq(3),locSideID)
    END DO ! i=0,PP_N
  END DO ! j=0,PP_N
END DO ! k=0,PP_N
#else /* for Gauss-Lobatto-points*/
CALL abort(__STAMP__,'this is wrong. see flexi, which has an additional S2V2 flip')
SELECT CASE(locSideID)
CASE(XI_MINUS)
  Uface(:,:,:)=Uvol(:,0,:,:)
CASE(ETA_MINUS)
  Uface(:,:,:)=Uvol(:,:,0,:)
CASE(ZETA_MINUS)
  Uface(:,:,:)=Uvol(:,:,:,0)
CASE(XI_PLUS)
  Uface(:,:,:)=Uvol(:,Nloc,:,:)
CASE(ETA_PLUS)
  Uface(:,:,:)=Uvol(:,:,Nloc,:)
CASE(ZETA_PLUS)
  Uface(:,:,:)=Uvol(:,:,:,Nloc)
END SELECT
#endif

END SUBROUTINE ProlongToFace_Elementlocal

END MODULE MOD_ProlongToFace
