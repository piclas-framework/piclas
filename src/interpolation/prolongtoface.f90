#include "boltzplatz.h"

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
  MODULE PROCEDURE ProlongToFace_sideBased
!  MODULE PROCEDURE ProlongToFace_SideBased2
!  MODULE PROCEDURE ProlongToFace_SideBased4
END INTERFACE

PUBLIC::ProlongToFace
!===================================================================================================================================

CONTAINS

SUBROUTINE ProlongToFace_SideBased(Uvol,Uface_Minus,Uface_Plus,doMPISides)
!===================================================================================================================================
! Interpolates the interior volume data (stored at the Gauss or Gauss-Lobatto points) to the surface
! integration points, using fast 1D Interpolation and store in global side structure
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Interpolation_Vars, ONLY: L_Minus,L_Plus
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: SideToElem
USE MOD_Mesh_Vars,          ONLY: nSides,nBCSides,nInnerSides,nMPISides_MINE,nMPISides_YOUR
USE MOD_Mesh_Vars,          ONLY: SideID_minus_lower,SideID_minus_upper
USE MOD_Mesh_Vars,          ONLY: SideID_plus_lower,SideID_plus_upper
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)              :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides +InnerSides +MPISides MINE 
REAL,INTENT(IN)                 :: Uvol(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)              :: Uface_Minus(PP_nVar,0:PP_N,0:PP_N,sideID_minus_lower:sideID_minus_upper)
REAL,INTENT(INOUT)              :: Uface_Plus(PP_nVar,0:PP_N,0:PP_N,sideID_plus_lower:sideID_plus_upper)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                         :: i,l,p,q,ElemID(2),SideID,flip(2),LocSideID(2),firstSideID,lastSideID
REAL                            :: Uface(PP_nVar,0:PP_N,0:PP_N)
!===================================================================================================================================
IF(doMPISides)THEN
  ! only YOUR MPI Sides are filled
  firstSideID = nBCSides+nInnerSides+nMPISides_MINE+1
  lastSideID  = firstSideID-1+nMPISides_YOUR 
  flip(1)      = -1
ELSE
  ! BCSides, InnerSides and MINE MPISides are filled
  firstSideID = 1
  lastSideID  = nBCSides+nInnerSides+nMPISides_MINE
  flip(1)      = 0
END IF
DO SideID=firstSideID,lastSideID
  ! master side, flip=0
  ElemID(1)     = SideToElem(S2E_ELEM_ID,SideID)  
  locSideID(1) = SideToElem(S2E_LOC_SIDE_ID,SideID)
  ! neighbor side !ElemID,locSideID and flip =-1 if not existing
  ElemID(2)     = SideToElem(S2E_NB_ELEM_ID,SideID)
  locSideID(2) = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  flip(2)      = SideToElem(S2E_FLIP,SideID)
  DO i=1,2 !first maste then slave side
#if (PP_NodeType==1) /* for Gauss-points*/
    SELECT CASE(locSideID(i))
    CASE(XI_MINUS)
      DO q=0,PP_N
        DO p=0,PP_N
          Uface(:,q,p)=Uvol(:,0,p,q,ElemID(i))*L_Minus(0)
          DO l=1,PP_N
            ! switch to right hand system
            Uface(:,q,p)=Uface(:,q,p)+Uvol(:,l,p,q,ElemID(i))*L_Minus(l)
          END DO ! l
        END DO ! p
      END DO ! q
    CASE(ETA_MINUS)
      DO q=0,PP_N
        DO p=0,PP_N
          Uface(:,p,q)=Uvol(:,p,0,q,ElemID(i))*L_Minus(0)
          DO l=1,PP_N
            Uface(:,p,q)=Uface(:,p,q)+Uvol(:,p,l,q,ElemID(i))*L_Minus(l)
          END DO ! l
        END DO ! p
      END DO ! q
    CASE(ZETA_MINUS)
      DO q=0,PP_N
        DO p=0,PP_N
          Uface(:,q,p)=Uvol(:,p,q,0,ElemID(i))*L_Minus(0)
          DO l=1,PP_N
            ! switch to right hand system
            Uface(:,q,p)=Uface(:,q,p)+Uvol(:,p,q,l,ElemID(i))*L_Minus(l)
          END DO ! l
        END DO ! p
      END DO ! q
    CASE(XI_PLUS)
      DO q=0,PP_N
        DO p=0,PP_N
          Uface(:,p,q)=Uvol(:,0,p,q,ElemID(i))*L_Plus(0)
          DO l=1,PP_N
            Uface(:,p,q)=Uface(:,p,q)+Uvol(:,l,p,q,ElemID(i))*L_Plus(l)
          END DO ! l
        END DO ! p
      END DO ! q
    CASE(ETA_PLUS)
      DO q=0,PP_N
        DO p=0,PP_N
          Uface(:,PP_N-p,q)=Uvol(:,p,0,q,ElemID(i))*L_Plus(0)
          DO l=1,PP_N
            ! switch to right hand system
            Uface(:,PP_N-p,q)=Uface(:,PP_N-p,q)+Uvol(:,p,l,q,ElemID(i))*L_Plus(l)
          END DO ! l
        END DO ! p
      END DO ! q
    CASE(ZETA_PLUS)
      DO q=0,PP_N
        DO p=0,PP_N
          Uface(:,p,q)=Uvol(:,p,q,0,ElemID(i))*L_Plus(0)
          DO l=1,PP_N
            Uface(:,p,q)=Uface(:,p,q)+Uvol(:,p,q,l,ElemID(i))*L_Plus(l)
          END DO ! l
        END DO ! p
      END DO ! q
    END SELECT
#else /* for Gauss-Lobatto-points*/
    SELECT CASE(locSideID(i))
    CASE(XI_MINUS)
      DO q=0,PP_N
        DO p=0,PP_N
          Uface(:,q,p)=Uvol(:,0,p,q,ElemID(i))
        END DO ! p
      END DO ! q
    CASE(ETA_MINUS)
      Uface(:,:,:)=Uvol(:,:,0,:,ElemID(i))
    CASE(ZETA_MINUS)
      DO q=0,PP_N
        DO p=0,PP_N
          Uface(:,q,p)=Uvol(:,p,q,0,ElemID(i))
        END DO ! p
      END DO ! q
    CASE(XI_PLUS)
      Uface(:,:,:)=Uvol(:,PP_N,:,:,ElemID(i))
    CASE(ETA_PLUS)
      DO q=0,PP_N
        DO p=0,PP_N
          Uface(:,PP_N-p,q)=Uvol(:,p,PP_N,q,ElemID(i))
        END DO ! p
      END DO ! q
    CASE(ZETA_PLUS)
      DO q=0,PP_N
        DO p=0,PP_N
          Uface(:,p,q)=Uvol(:,p,q,PP_N,ElemID(i))
        END DO ! p
      END DO ! q
    END SELECT
#endif
    SELECT CASE(Flip(i))
      CASE(0) ! master side
        Uface_Minus(:,:,:,SideID)=Uface(:,:,:)
      CASE(1) ! slave side, SideID=q,jSide=p
        DO q=0,PP_N
          DO p=0,PP_N
            Uface_Plus(:,p,q,SideID)=Uface(:,q,p)
          END DO ! p
        END DO ! q
      CASE(2) ! slave side, SideID=N-p,jSide=q
        DO q=0,PP_N
          DO p=0,PP_N
            Uface_Plus(:,p,q,SideID)=Uface(:,PP_N-p,q)
          END DO ! p
        END DO ! q
      CASE(3) ! slave side, SideID=N-q,jSide=N-p
        DO q=0,PP_N
          DO p=0,PP_N
            Uface_Plus(:,p,q,SideID)=Uface(:,PP_N-q,PP_N-p)
          END DO ! p
        END DO ! q
      CASE(4) ! slave side, SideID=p,jSide=N-q
        DO q=0,PP_N
          DO p=0,PP_N
            Uface_Plus(:,p,q,SideID)=Uface(:,p,PP_N-q)
          END DO ! p
        END DO ! q
    END SELECT
  END DO !i=1,2, masterside & slave side 
END DO !SideID
END SUBROUTINE ProlongToFace_SideBased

#ifdef dontcompilethis
SUBROUTINE ProlongToFace_SideBased2(Uvol,Uface_Minus,Uface_Plus,doMPISides)
!===================================================================================================================================
! Interpolates the interior volume data (stored at the Gauss or Gauss-Lobatto points) to the surface
! integration points, using fast 1D Interpolation and store in global side structure
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Interpolation_Vars, ONLY: L_Minus,L_Plus
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: SideToElem2
USE MOD_Mesh_Vars,          ONLY: nSides,nBCSides,nInnerSides,nMPISides_MINE,nMPISides_YOUR
USE MOD_Mesh_Vars,          ONLY: SideID_minus_lower,SideID_minus_upper
USE MOD_Mesh_Vars,          ONLY: SideID_plus_lower,SideID_plus_upper
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)              :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides +InnerSides +MPISides MINE 
REAL,INTENT(IN)                 :: Uvol(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)              :: Uface_Minus(PP_nVar,0:PP_N,0:PP_N,sideID_minus_lower:sideID_minus_upper)
REAL,INTENT(INOUT)              :: Uface_Plus(PP_nVar,0:PP_N,0:PP_N,sideID_plus_lower:sideID_plus_upper)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                         :: l,p,q,ElemID,iSide,SideID,flip,LocSideID,firstSide,lastSide
REAL                            :: Uface(PP_nVar,0:PP_N,0:PP_N)
!===================================================================================================================================
IF(doMPISides)THEN
  ! only YOUR MPI Sides are filled
  firstSide = nBCSides+2*nInnerSides+nMPISides_MINE+1
  lastSide  = firstSide-1+nMPISides_YOUR
ELSE
  ! BCSides, InnerSides and MINE MPISides are filled
  firstSide = 1
  lastSide  = nBCSides+2*nInnerSides+nMPISides_MINE
END IF
DO iSide=firstSide,lastSide
  ! master side, flip=0
  ElemID    = SideToElem2(S2E2_ELEM_ID,iSide)  
  SideID    = SideToElem2(S2E2_SIDE_ID,iSide)  
  locSideID = SideToElem2(S2E2_LOC_SIDE_ID,iSide)
  flip      = SideToElem2(S2E2_FLIP,iSide)
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
  SELECT CASE(Flip)
    CASE(0) ! master side
      Uface_Minus(:,:,:,SideID)=Uface(:,:,:)
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,PP_N
        DO p=0,PP_N
          Uface_Plus(:,p,q,SideID)=Uface(:,q,p)
        END DO ! p
      END DO ! q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,PP_N
        DO p=0,PP_N
          Uface_Plus(:,p,q,SideID)=Uface(:,PP_N-p,q)
        END DO ! p
      END DO ! q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,PP_N
        DO p=0,PP_N
          Uface_Plus(:,p,q,SideID)=Uface(:,PP_N-q,PP_N-p)
        END DO ! p
      END DO ! q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,PP_N
        DO p=0,PP_N
          Uface_Plus(:,p,q,SideID)=Uface(:,p,PP_N-q)
        END DO ! p
      END DO ! q
  END SELECT
END DO !iSide
END SUBROUTINE ProlongToFace_SideBased2


SUBROUTINE ProlongToFace_SideBased4(Uvol,Uface_Minus,Uface_Plus,doMPISides)
!===================================================================================================================================
! Interpolates the interior volume data (stored at the Gauss or Gauss-Lobatto points) to the surface
! integration points, using fast 1D Interpolation and store in global side structure
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Interpolation_Vars, ONLY: L_Minus,L_Plus
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: SideToElem
USE MOD_Mesh_Vars,          ONLY: nSides,nBCSides,nInnerSides,nMPISides_MINE,nMPISides_YOUR
USE MOD_Mesh_Vars,          ONLY: SideID_minus_lower,SideID_minus_upper
USE MOD_Mesh_Vars,          ONLY: SideID_plus_lower,SideID_plus_upper
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)              :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides +InnerSides +MPISides MINE 
REAL,INTENT(IN)                 :: Uvol(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)              :: Uface_Minus(PP_nVar,0:PP_N,0:PP_N,sideID_minus_lower:sideID_minus_upper)
REAL,INTENT(INOUT)              :: Uface_Plus(PP_nVar,0:PP_N,0:PP_N,sideID_plus_lower:sideID_plus_upper)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                         :: i,l,p,q,ElemID(2),SideID,flip(2),LocSideID(2),firstSideID,lastSideID
!===================================================================================================================================
IF(doMPISides)THEN
  ! only YOUR MPI Sides are filled
  firstSideID = nBCSides+nInnerSides+nMPISides_MINE+1
  lastSideID  = firstSideID-1+nMPISides_YOUR 
  flip(1)      = -1 ! keine 
ELSE
  ! BCSides, InnerSides and MINE MPISides are filled
  firstSideID = 1
  lastSideID  = nBCSides+nInnerSides+nMPISides_MINE
  flip(1)      = 0
END IF

DO SideID=firstSideID,lastSideID
  ! master side, flip=0
  ElemID(1)     = SideToElem(S2E_ELEM_ID,SideID)  
  locSideID(1) = SideToElem(S2E_LOC_SIDE_ID,SideID)
  ! neighbor side !ElemID,locSideID and flip =-1 if not existing
  ElemID(2)     = SideToElem(S2E_NB_ELEM_ID,SideID)
  locSideID(2) = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  flip(2)      = SideToElem(S2E_FLIP,SideID)
  DO i=1,2 !first maste then slave side
    SELECT CASE(flip(i))
    CASE(0)
      SELECT CASE(locSideID(i))
      CASE(XI_MINUS)
        DO q=0,PP_N; DO p=0,PP_N
          Uface_minus(:,q,p,SideID)=Uvol(:,0,p,q,ElemID(i))*L_Minus(0)
          DO l=1,PP_N ! switch to right hand system
            Uface_minus(:,q,p,SideID)=Uface_minus(:,q,p,SideID)+Uvol(:,l,p,q,ElemID(i))*L_Minus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(ETA_MINUS)
        DO q=0,PP_N; DO p=0,PP_N
          Uface_minus(:,p,q,SideID)=Uvol(:,p,0,q,ElemID(i))*L_Minus(0)
          DO l=1,PP_N
            Uface_minus(:,p,q,SideID)=Uface_minus(:,p,q,SideID)+Uvol(:,p,l,q,ElemID(i))*L_Minus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(ZETA_MINUS)
        DO q=0,PP_N; DO p=0,PP_N
            Uface_minus(:,q,p,SideID)=Uvol(:,p,q,0,ElemID(i))*L_Minus(0)
            DO l=1,PP_N ! switch to right hand system
              Uface_minus(:,q,p,SideID)=Uface_minus(:,q,p,SideID)+Uvol(:,p,q,l,ElemID(i))*L_Minus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(XI_PLUS)
        DO q=0,PP_N; DO p=0,PP_N
            Uface_minus(:,p,q,SideID)=Uvol(:,0,p,q,ElemID(i))*L_Plus(0)
            DO l=1,PP_N
              Uface_minus(:,p,q,SideID)=Uface_minus(:,p,q,SideID)+Uvol(:,l,p,q,ElemID(i))*L_Plus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(ETA_PLUS)
        DO q=0,PP_N; DO p=0,PP_N
          Uface_minus(:,PP_N-p,q,SideID)=Uvol(:,p,0,q,ElemID(i))*L_Plus(0)
          DO l=1,PP_N ! switch to right hand system
            Uface_minus(:,PP_N-p,q,SideID)=Uface_minus(:,PP_N-p,q,SideID)+Uvol(:,p,l,q,ElemID(i))*L_Plus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(ZETA_PLUS)
        DO q=0,PP_N; DO p=0,PP_N
          Uface_minus(:,p,q,SideID)=Uvol(:,p,q,0,ElemID(i))*L_Plus(0)
          DO l=1,PP_N
            Uface_minus(:,p,q,SideID)=Uface_minus(:,p,q,SideID)+Uvol(:,p,q,l,ElemID(i))*L_Plus(l)
        END DO; END DO; END DO ! l,p,q
      END SELECT
 
    CASE(1)
      SELECT CASE(locSideID(i))
      CASE(XI_MINUS)
        DO q=0,PP_N; DO p=0,PP_N
          Uface_plus(:,p,q,SideID)=Uvol(:,0,p,q,ElemID(i))*L_Minus(0)
          DO l=1,PP_N ! switch to right hand system
            Uface_plus(:,p,q,SideID)=Uface_plus(:,p,q,SideID)+Uvol(:,l,p,q,ElemID(i))*L_Minus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(ETA_MINUS)
        DO q=0,PP_N; DO p=0,PP_N
          Uface_plus(:,q,p,SideID)=Uvol(:,p,0,q,ElemID(i))*L_Minus(0)
          DO l=1,PP_N
            Uface_plus(:,q,p,SideID)=Uface_plus(:,q,p,SideID)+Uvol(:,p,l,q,ElemID(i))*L_Minus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(ZETA_MINUS)
        DO q=0,PP_N; DO p=0,PP_N
            Uface_plus(:,p,q,SideID)=Uvol(:,p,q,0,ElemID(i))*L_Minus(0)
            DO l=1,PP_N ! switch to right hand system
              Uface_plus(:,p,q,SideID)=Uface_plus(:,p,q,SideID)+Uvol(:,p,q,l,ElemID(i))*L_Minus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(XI_PLUS)
        DO q=0,PP_N; DO p=0,PP_N
            Uface_plus(:,q,p,SideID)=Uvol(:,0,p,q,ElemID(i))*L_Plus(0)
            DO l=1,PP_N
              Uface_plus(:,q,p,SideID)=Uface_plus(:,q,p,SideID)+Uvol(:,l,p,q,ElemID(i))*L_Plus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(ETA_PLUS)
        DO q=0,PP_N; DO p=0,PP_N
          Uface_plus(:,q,PP_N-p,SideID)=Uvol(:,p,0,q,ElemID(i))*L_Plus(0)
          DO l=1,PP_N ! switch to right hand system
            Uface_plus(:,q,PP_N-p,SideID)=Uface_plus(:,q,PP_N-p,SideID)+Uvol(:,p,l,q,ElemID(i))*L_Plus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(ZETA_PLUS)
        DO q=0,PP_N; DO p=0,PP_N
          Uface_plus(:,q,p,SideID)=Uvol(:,p,q,0,ElemID(i))*L_Plus(0)
          DO l=1,PP_N
            Uface_plus(:,q,p,SideID)=Uface_plus(:,q,p,SideID)+Uvol(:,p,q,l,ElemID(i))*L_Plus(l)
        END DO; END DO; END DO ! l,p,q
      END SELECT
 
    CASE(2)
      SELECT CASE(locSideID(i))
      CASE(XI_MINUS)
        DO q=0,PP_N; DO p=0,PP_N
          Uface_plus(:,PP_N-q,p,SideID)=Uvol(:,0,p,q,ElemID(i))*L_Minus(0)
          DO l=1,PP_N ! switch to right hand system
            Uface_plus(:,PP_N-p,q,SideID)=Uface_plus(:,PP_N-q,p,SideID)+Uvol(:,l,p,q,ElemID(i))*L_Minus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(ETA_MINUS)
        DO q=0,PP_N; DO p=0,PP_N
          Uface_plus(:,PP_N-p,q,SideID)=Uvol(:,p,0,q,ElemID(i))*L_Minus(0)
          DO l=1,PP_N
            Uface_plus(:,PP_N-p,q,SideID)=Uface_plus(:,PP_N-p,q,SideID)+Uvol(:,p,l,q,ElemID(i))*L_Minus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(ZETA_MINUS)
        DO q=0,PP_N; DO p=0,PP_N
            Uface_plus(:,PP_N-q,p,SideID)=Uvol(:,p,q,0,ElemID(i))*L_Minus(0)
            DO l=1,PP_N ! switch to right hand system
              Uface_plus(:,PP_N-q,p,SideID)=Uface_plus(:,PP_N-q,p,SideID)+Uvol(:,p,q,l,ElemID(i))*L_Minus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(XI_PLUS)
        DO q=0,PP_N; DO p=0,PP_N
            Uface_plus(:,PP_N-p,q,SideID)=Uvol(:,0,p,q,ElemID(i))*L_Plus(0)
            DO l=1,PP_N
              Uface_plus(:,PP_N-p,q,SideID)=Uface_plus(:,PP_N-p,q,SideID)+Uvol(:,l,p,q,ElemID(i))*L_Plus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(ETA_PLUS)
        DO q=0,PP_N; DO p=0,PP_N
          Uface_plus(:,p,q,SideID)=Uvol(:,p,0,q,ElemID(i))*L_Plus(0)
          DO l=1,PP_N ! switch to right hand system
            Uface_plus(:,p,q,SideID)=Uface_plus(:,p,q,SideID)+Uvol(:,p,l,q,ElemID(i))*L_Plus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(ZETA_PLUS)
        DO q=0,PP_N; DO p=0,PP_N
          Uface_plus(:,PP_N-p,q,SideID)=Uvol(:,p,q,0,ElemID(i))*L_Plus(0)
          DO l=1,PP_N
            Uface_plus(:,PP_N-p,q,SideID)=Uface_plus(:,PP_N-p,q,SideID)+Uvol(:,p,q,l,ElemID(i))*L_Plus(l)
        END DO; END DO; END DO ! l,p,q
      END SELECT
 
    CASE(3)
      SELECT CASE(locSideID(i))
      CASE(XI_MINUS)
        DO q=0,PP_N; DO p=0,PP_N
          Uface_plus(:,PP_N-p,PP_N-q,SideID)=Uvol(:,0,p,q,ElemID(i))*L_Minus(0)
          DO l=1,PP_N ! switch to right hand system
            Uface_plus(:,PP_N-p,PP_N-q,SideID)=Uface_plus(:,PP_N-p,PP_N-q,SideID)+Uvol(:,l,p,q,ElemID(i))*L_Minus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(ETA_MINUS)
        DO q=0,PP_N; DO p=0,PP_N
          Uface_plus(:,PP_N-q,PP_N-p,SideID)=Uvol(:,p,0,q,ElemID(i))*L_Minus(0)
          DO l=1,PP_N
            Uface_plus(:,PP_N-q,PP_N-p,SideID)=Uface_plus(:,PP_N-q,PP_N-p,SideID)+Uvol(:,p,l,q,ElemID(i))*L_Minus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(ZETA_MINUS)
        DO q=0,PP_N; DO p=0,PP_N
            Uface_plus(:,PP_N-p,PP_N-q,SideID)=Uvol(:,p,q,0,ElemID(i))*L_Minus(0)
            DO l=1,PP_N ! switch to right hand system
              Uface_plus(:,PP_N-p,PP_N-q,SideID)=Uface_plus(:,PP_N-p,PP_N-q,SideID)+Uvol(:,p,q,l,ElemID(i))*L_Minus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(XI_PLUS)
        DO q=0,PP_N; DO p=0,PP_N
            Uface_plus(:,PP_N-q,PP_N-p,SideID)=Uvol(:,0,p,q,ElemID(i))*L_Plus(0)
            DO l=1,PP_N
              Uface_plus(:,PP_N-q,PP_N-p,SideID)=Uface_plus(:,PP_N-q,PP_N-p,SideID)+Uvol(:,l,p,q,ElemID(i))*L_Plus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(ETA_PLUS)
        DO q=0,PP_N; DO p=0,PP_N
          Uface_plus(:,PP_N-q,p,SideID)=Uvol(:,p,0,q,ElemID(i))*L_Plus(0)
          DO l=1,PP_N ! switch to right hand system
            Uface_plus(:,PP_N-q,p,SideID)=Uface_plus(:,PP_N-q,p,SideID)+Uvol(:,p,l,q,ElemID(i))*L_Plus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(ZETA_PLUS)
        DO q=0,PP_N; DO p=0,PP_N
          Uface_plus(:,PP_N-q,PP_N-p,SideID)=Uvol(:,p,q,0,ElemID(i))*L_Plus(0)
          DO l=1,PP_N
            Uface_plus(:,PP_N-q,PP_N-p,SideID)=Uface_plus(:,PP_N-q,PP_N-p,SideID)+Uvol(:,p,q,l,ElemID(i))*L_Plus(l)
        END DO; END DO; END DO ! l,p,q
      END SELECT
      
    CASE(4)
      SELECT CASE(locSideID(i))
      CASE(XI_MINUS)
        DO q=0,PP_N; DO p=0,PP_N
          Uface_plus(:,q,PP_N-p,SideID)=Uvol(:,0,p,q,ElemID(i))*L_Minus(0)
          DO l=1,PP_N ! switch to right hand system
            Uface_plus(:,q,PP_N-p,SideID)=Uface_plus(:,q,PP_N-p,SideID)+Uvol(:,l,p,q,ElemID(i))*L_Minus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(ETA_MINUS)
        DO q=0,PP_N; DO p=0,PP_N
          Uface_plus(:,p,PP_N-q,SideID)=Uvol(:,p,0,q,ElemID(i))*L_Minus(0)
          DO l=1,PP_N
            Uface_plus(:,p,PP_N-q,SideID)=Uface_plus(:,p,PP_N-q,SideID)+Uvol(:,p,l,q,ElemID(i))*L_Minus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(ZETA_MINUS)
        DO q=0,PP_N; DO p=0,PP_N
            Uface_plus(:,q,PP_N-p,SideID)=Uvol(:,p,q,0,ElemID(i))*L_Minus(0)
            DO l=1,PP_N ! switch to right hand system
              Uface_plus(:,q,PP_N-p,SideID)=Uface_plus(:,q,PP_N-p,SideID)+Uvol(:,p,q,l,ElemID(i))*L_Minus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(XI_PLUS)
        DO q=0,PP_N; DO p=0,PP_N
            Uface_plus(:,p,PP_N-q,SideID)=Uvol(:,0,p,q,ElemID(i))*L_Plus(0)
            DO l=1,PP_N
              Uface_plus(:,p,PP_N-q,SideID)=Uface_plus(:,p,PP_N-q,SideID)+Uvol(:,l,p,q,ElemID(i))*L_Plus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(ETA_PLUS)
        DO q=0,PP_N; DO p=0,PP_N
          Uface_plus(:,PP_N-p,PP_N-q,SideID)=Uvol(:,p,0,q,ElemID(i))*L_Plus(0)
          DO l=1,PP_N ! switch to right hand system
            Uface_plus(:,PP_N-p,PP_N-q,SideID)=Uface_plus(:,PP_N-p,PP_N-q,SideID)+Uvol(:,p,l,q,ElemID(i))*L_Plus(l)
        END DO; END DO; END DO ! l,p,q
      CASE(ZETA_PLUS)
        DO q=0,PP_N; DO p=0,PP_N
          Uface_plus(:,p,PP_N-q,SideID)=Uvol(:,p,q,0,ElemID(i))*L_Plus(0)
          DO l=1,PP_N
            Uface_plus(:,p,PP_N-q,SideID)=Uface_plus(:,p,PP_N-q,SideID)+Uvol(:,p,q,l,ElemID(i))*L_Plus(l)
        END DO; END DO; END DO ! l,p,q
      END SELECT
    END SELECT
  END DO !i=1,2, masterside & slave side 
END DO !SideID
END SUBROUTINE ProlongToFace_SideBased4
#endif dontcompilethis


END MODULE MOD_ProlongToFace
