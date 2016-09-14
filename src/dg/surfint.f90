#include "boltzplatz.h"

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
  !MODULE PROCEDURE SurfInt1
  MODULE PROCEDURE SurfInt2
END INTERFACE

PUBLIC::SurfInt
!===================================================================================================================================
CONTAINS

#ifdef DONTCOMPILETHIS
SUBROUTINE SurfInt1(Flux,Ut,doMPISides)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,            ONLY: L_HatPlus,L_HatMinus
USE MOD_Mesh_Vars,          ONLY: SideToElem
USE MOD_Mesh_Vars,          ONLY: nSides,nBCSides,nInnerSides,nMPISides_MINE,nMPISides_YOUR
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides+InnerSides+MPISides MINE  
REAL,INTENT(IN)    :: Flux(1:PP_nVar,0:PP_N,0:PP_N,nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
REAL,INTENT(INOUT)   :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: Flux_loc(1:PP_nVar,0:PP_N,0:PP_N)
INTEGER            :: i,ElemID(2),p,q,l,Flip(2),SideID,locSideID(2)
INTEGER            :: firstSideID,lastSideID
!===================================================================================================================================
IF(doMPISides)THEN 
  ! surfInt only for YOUR MPISides
  firstSideID = nBCSides+nInnerSides+nMPISides_MINE +1
  lastSideID  = firstSideID-1+nMPISides_YOUR 
ELSE
  ! fill only InnerSides
  firstSideID = 1
  lastSideID  = nBCSides+nInnerSides+nMPISides_MINE
END IF
DO SideID=firstSideID,lastSideID
  ! master side, flip=0
  ElemID(1)     = SideToElem(S2E_ELEM_ID,SideID)  
  locSideID(1) = SideToElem(S2E_LOC_SIDE_ID,SideID)
  flip(1)      = 0
  ! neighbor side
  ElemID(2)     = SideToElem(S2E_NB_ELEM_ID,SideID)
  locSideID(2) = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  flip(2)      = SideToElem(S2E_FLIP,SideID)
  DO i=1,2
    ! prepare flux (with correct inverse flip and sign)
    ! store local face data into global data structure U_Minus and U_Plus
    SELECT CASE(Flip(i))
      CASE(0) ! master side
        Flux_loc(:,:,:)=Flux(:,:,:,SideID)
      CASE(1) ! slave side, SideID=q,jSide=p
        DO q=0,PP_N
          DO p=0,PP_N
            Flux_loc(:,q,p)=-Flux(:,p,q,SideID)
          END DO ! p
        END DO ! q
      CASE(2) ! slave side, SideID=N-p,jSide=q
        DO q=0,PP_N
          DO p=0,PP_N
            Flux_loc(:,PP_N-p,q)=-Flux(:,p,q,SideID)
          END DO ! p
        END DO ! q
      CASE(3) ! slave side, SideID=N-q,jSide=N-p
        DO q=0,PP_N
          DO p=0,PP_N
            Flux_loc(:,PP_N-q,PP_N-p)=-Flux(:,p,q,SideID)
          END DO ! p
        END DO ! q
      CASE(4) ! slave side, SideID=p,jSide=N-q
        DO q=0,PP_N
          DO p=0,PP_N
            Flux_loc(:,p,PP_N-q)=-Flux(:,p,q,SideID)
          END DO ! p
        END DO ! q
    END SELECT

  ! update DG time derivative with corresponding SurfInt contribution
#if (PP_NodeType==1)
    SELECT CASE(locSideID(i))
    CASE(XI_MINUS)
      DO q=0,PP_N
        DO p=0,PP_N
          DO l=0,PP_N
            ! update of local grid cell
            ! switch to right hand system for XI_MINUS direction
            Ut(:,l,p,q,ElemID(i))=Ut(:,l,p,q,ElemID(i))+Flux_loc(:,q,p)*L_hatMinus(l)
          END DO ! l
        END DO ! p
      END DO ! q
    CASE(ETA_MINUS)
      DO q=0,PP_N
        DO l=0,PP_N
          DO p=0,PP_N
            ! update of local grid cell
            ! switch to right hand system for ETA_PLUS direction
            Ut(:,p,l,q,ElemID(i))=Ut(:,p,l,q,ElemID(i))+Flux_loc(:,p,q)*L_hatMinus(l)
          END DO ! l
        END DO ! p
      END DO ! q
    CASE(ZETA_MINUS)
      DO l=0,PP_N
        DO q=0,PP_N
          DO p=0,PP_N
            ! update of local grid cell
            ! switch to right hand system for ZETA_MINUS direction
            Ut(:,p,q,l,ElemID(i))=Ut(:,p,q,l,ElemID(i))+Flux_loc(:,q,p)*L_hatMinus(l)
          END DO ! l
        END DO ! p
      END DO ! q
    CASE(XI_PLUS)
      DO q=0,PP_N
        DO p=0,PP_N
          DO l=0,PP_N
            ! update of local grid cell
            ! switch to right hand system for XI_MINUS direction
            Ut(:,l,p,q,ElemID(i))=Ut(:,l,p,q,ElemID(i))+Flux_loc(:,p,q)*L_hatPlus(l)
          END DO ! l
        END DO ! p
      END DO ! q
    CASE(ETA_PLUS)
      DO q=0,PP_N
        DO l=0,PP_N
          DO p=0,PP_N
            ! update of local grid cell
            ! switch to right hand system for ETA_PLUS direction
            Ut(:,p,l,q,ElemID(i))=Ut(:,p,l,q,ElemID(i))+Flux_loc(:,PP_N-p,q)*L_hatPlus(l)
          END DO ! l
        END DO ! p
      END DO ! q
    CASE(ZETA_PLUS)
      DO l=0,PP_N
        DO q=0,PP_N
          DO p=0,PP_N
            ! update of local grid cell
            ! switch to right hand system for ZETA_MINUS direction
            Ut(:,p,q,l,ElemID(i))=Ut(:,p,q,l,ElemID(i))+Flux_loc(:,p,q)*L_hatPlus(l)
          END DO ! l
        END DO ! p
      END DO ! q
    END SELECT !locSideID
#else
    SELECT CASE(locSideID(i))
    CASE(XI_MINUS)
      DO q=0,PP_N
        DO p=0,PP_N
          Ut(:,0,p,q,ElemID(i))=Ut(:,0,p,q,ElemID(i))+Flux_loc(:,q,p)*L_hatMinus(0)
        END DO ! p
      END DO ! q
    CASE(ETA_MINUS)
      DO q=0,PP_N
          DO p=0,PP_N
            ! update of local grid cell
            ! switch to right hand system for ETA_PLUS direction
            Ut(:,p,0,q,ElemID(i))=Ut(:,p,0,q,ElemID(i))+Flux_loc(:,p,q)*L_hatMinus(0)
        END DO ! p
      END DO ! q
    CASE(ZETA_MINUS)
      DO q=0,PP_N
        DO p=0,PP_N
          ! update of local grid cell
          ! switch to right hand system for ZETA_MINUS direction
          Ut(:,p,q,0,ElemID(i))=Ut(:,p,q,0,ElemID(i))+Flux_loc(:,q,p)*L_hatMinus(0)
        END DO ! p
      END DO ! q
    CASE(XI_PLUS)
      DO q=0,PP_N
        DO p=0,PP_N
            ! update of local grid cell
            ! switch to right hand system for XI_MINUS direction
            Ut(:,PP_N,p,q,ElemID(i))=Ut(:,PP_N,p,q,ElemID(i))+Flux_loc(:,p,q)*L_hatPlus(PP_N)
        END DO ! p
      END DO ! q
    CASE(ETA_PLUS)
      DO q=0,PP_N
          DO p=0,PP_N
            ! update of local grid cell
            ! switch to right hand system for ETA_PLUS direction
            Ut(:,p,PP_N,q,ElemID(i))=Ut(:,p,PP_N,q,ElemID(i))+Flux_loc(:,PP_N-p,q)*L_hatPlus(PP_N)
        END DO ! p
      END DO ! q
    CASE(ZETA_PLUS)
      DO q=0,PP_N
        DO p=0,PP_N
          ! update of local grid cell
          ! switch to right hand system for ZETA_MINUS direction
          Ut(:,p,q,PP_N,ElemID(i))=Ut(:,p,q,PP_N,ElemID(i))+Flux_loc(:,p,q)*L_hatPlus(PP_N)
        END DO ! p
      END DO ! q
    END SELECT !locSideID
#endif /*(PP_NodeType==1)*/
  END DO ! i=1,2 master side, slave side
END DO ! SideID=1,nSides
END SUBROUTINE SurfInt1
#endif /*DONTCOMPILETHIS*/

! surfint optimized for performance (but ugly to read)
SUBROUTINE SurfInt2(Flux,Ut,doMPISides)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,            ONLY: L_HatPlus,L_HatMinus
USE MOD_Mesh_Vars,          ONLY: SideToElem
USE MOD_Mesh_Vars,          ONLY: nSides,nBCSides,nInnerSides,nMPISides_MINE,nMPISides_YOUR
USE MOD_PML_Vars,           ONLY: DoPML,PMLnVar,ElemToPML,U2t,isPMLElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides+InnerSides+MPISides MINE  
REAL,INTENT(IN)    :: Flux(1:PP_nVar+PMLnVar,0:PP_N,0:PP_N,nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)   :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: ElemID,p,q,l,Flip,SideID,locSideID
INTEGER            :: firstSideID,lastSideID
#if (PP_NodeType>1)
REAL               ::L_HatMinus0,L_HatPlusN 
#endif
!===================================================================================================================================
IF(doMPISides)THEN 
  ! surfInt only for YOUR MPISides
  firstSideID = nBCSides+nInnerSides+nMPISides_MINE +1
  lastSideID  = firstSideID-1+nMPISides_YOUR 
ELSE
  ! fill only InnerSides
  firstSideID = 1
  lastSideID  = nBCSides+nInnerSides+nMPISides_MINE
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
  IF(ElemID.LT.1) CYCLE 
  IF(DoPML)THEN 
    IF(isPMLElem(ElemID))THEN
      CALL CalcSurfInt2PML(Flux,Ut,flip,ElemID,locSideID,SideID)
    ELSE
      CALL CalcSurfInt2(Flux(1:PP_nVar,:,:,SideID),Ut,Flip,ElemID,locSideID)
    END IF
  ELSE
    CALL CalcSurfInt2(Flux(1:PP_nVar,:,:,SideID),Ut,Flip,ElemID,locSideID)
  END IF
END DO ! SideID=1,nSides


DO SideID=firstSideID,lastSideID
  ! master side, flip=0
  ElemID    = SideToElem(S2E_ELEM_ID,SideID)  
  locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
  flip      = 0
  IF(DoPML)THEN 
    IF(isPMLElem(ElemID))THEN
      CALL CalcSurfInt2PML(Flux,Ut,flip,ElemID,locSideID,SideID)
    ELSE
      CALL CalcSurfInt2(Flux(1:PP_nVar,:,:,SideID),Ut,Flip,ElemID,locSideID)
    END IF
  ELSE
    CALL CalcSurfInt2(Flux(1:PP_nVar,:,:,SideID),Ut,Flip,ElemID,locSideID)
  END IF
END DO ! SideID=1,nSides

END SUBROUTINE SurfInt2


SUBROUTINE CalcSurfInt2(Flux,Ut,flip,ElemID,locSideID)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,            ONLY: L_HatPlus,L_HatMinus
USE MOD_Mesh_Vars,          ONLY: SideToElem
USE MOD_Mesh_Vars,          ONLY: nSides,nBCSides,nInnerSides,nMPISides_MINE,nMPISides_YOUR
USE MOD_PML_Vars,           ONLY: DoPML,PMLnVar,ElemToPML,U2t,isPMLElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!REAL,INTENT(IN)    :: Flux(1:PP_nVar+PMLnVar,0:PP_N,0:PP_N,nSides)
REAL,INTENT(IN)    :: Flux(1:PP_nVar,0:PP_N,0:PP_N)
INTEGER,INTENT(IN) :: flip,ElemID,locSideID!,SideID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)   :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q,l
!INTEGER            :: firstSideID,lastSideID
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


SUBROUTINE CalcSurfInt2PML(Flux,Ut,flip,ElemID,locSideID,SideID)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,            ONLY: L_HatPlus,L_HatMinus
USE MOD_Mesh_Vars,          ONLY: SideToElem
USE MOD_Mesh_Vars,          ONLY: nSides,nBCSides,nInnerSides,nMPISides_MINE,nMPISides_YOUR
USE MOD_PML_Vars,           ONLY: DoPML,PMLnVar,ElemToPML,U2t,isPMLElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: Flux(1:PP_nVar+PMLnVar,0:PP_N,0:PP_N,nSides)
INTEGER,INTENT(IN) :: flip,ElemID,locSideID,SideID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)   :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,p,q,l
INTEGER            :: firstSideID,lastSideID
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
        Ut( :,l,p,q,ElemID)           =Ut( :,l,p,q,ElemID)           +Flux(1: 8,q,p,SideID)*L_hatMinus(l)
        U2t(:,l,p,q,ElemToPML(ElemID))=U2t(:,l,p,q,ElemToPML(ElemID))+Flux(9:32,q,p,SideID)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut( :,l,p,q,ElemID)           =Ut( :,l,p,q,ElemID)           -Flux(1: 8,p,q,SideID)*L_hatMinus(l)
        U2t(:,l,p,q,ElemToPML(ElemID))=U2t(:,l,p,q,ElemToPML(ElemID))-Flux(9:32,p,q,SideID)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut( :,l,p,q,ElemID)           =Ut( :,l,p,q,ElemID)           -Flux(1: 8,PP_N-q,p,SideID)*L_hatMinus(l)
        U2t(:,l,p,q,ElemToPML(ElemID))=U2t(:,l,p,q,ElemToPML(ElemID))-Flux(9:32,PP_N-q,p,SideID)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut( :,l,p,q,ElemID)           =Ut( :,l,p,q,ElemID)           -Flux(1: 8,PP_N-p,PP_N-q,SideID)*L_hatMinus(l)
        U2t(:,l,p,q,ElemToPML(ElemID))=U2t(:,l,p,q,ElemToPML(ElemID))-Flux(9:32,PP_N-p,PP_N-q,SideID)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut( :,l,p,q,ElemID)           =Ut( :,l,p,q,ElemID)           -Flux(1: 8,q,PP_N-p,SideID)*L_hatMinus(l)
        U2t(:,l,p,q,ElemToPML(ElemID))=U2t(:,l,p,q,ElemToPML(ElemID))-Flux(9:32,q,PP_N-p,SideID)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    END SELECT
!===================================================================================================================================
  CASE(ETA_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut( :,p,l,q,ElemID)           =Ut( :,p,l,q,ElemID)           +Flux(1: 8,p,q,SideID)*L_hatMinus(l)
        U2t(:,p,l,q,ElemToPML(ElemID))=U2t(:,p,l,q,ElemToPML(ElemID))+Flux(9:32,p,q,SideID)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut( :,p,l,q,ElemID)           =Ut( :,p,l,q,ElemID)           -Flux(1: 8,q,p,SideID)*L_hatMinus(l)
        U2t(:,p,l,q,ElemToPML(ElemID))=U2t(:,p,l,q,ElemToPML(ElemID))-Flux(9:32,q,p,SideID)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut( :,p,l,q,ElemID)           =Ut( :,p,l,q,ElemID)           -Flux(1: 8,PP_N-p,q,SideID)*L_hatMinus(l)
        U2t(:,p,l,q,ElemToPML(ElemID))=U2t(:,p,l,q,ElemToPML(ElemID))-Flux(9:32,PP_N-p,q,SideID)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut( :,p,l,q,ElemID)           =Ut( :,p,l,q,ElemID)           -Flux(1: 8,PP_N-q,PP_N-p,SideID)*L_hatMinus(l)
        U2t(:,p,l,q,ElemToPML(ElemID))=U2t(:,p,l,q,ElemToPML(ElemID))-Flux(9:32,PP_N-q,PP_N-p,SideID)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut( :,p,l,q,ElemID)           =Ut( :,p,l,q,ElemID)           -Flux(1: 8,p,PP_N-q,SideID)*L_hatMinus(l)
        U2t(:,p,l,q,ElemToPML(ElemID))=U2t(:,p,l,q,ElemToPML(ElemID))-Flux(9:32,p,PP_N-q,SideID)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    END SELECT
!===================================================================================================================================
  CASE(ZETA_MINUS)
!===================================================================================================================================

    SELECT CASE(flip)
    CASE(0) ! master side
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)           = Ut( :,p,q,l,ElemID)           +Flux(1: 8,q,p,SideID)*L_hatMinus(l)
        U2t(:,p,q,l,ElemToPML(ElemID))=U2t(:,p,q,l,ElemToPML(ElemID))+Flux(9:32,q,p,SideID)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(1) ! slave side, SideID=q,jSide=p
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)           = Ut( :,p,q,l,ElemID)           -Flux(1: 8,p,q,SideID)*L_hatMinus(l)
        U2t(:,p,q,l,ElemToPML(ElemID))=U2t(:,p,q,l,ElemToPML(ElemID))-Flux(9:32,p,q,SideID)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)           = Ut( :,p,q,l,ElemID)           -Flux(1: 8,PP_N-q,p,SideID)*L_hatMinus(l)
        U2t(:,p,q,l,ElemToPML(ElemID))=U2t(:,p,q,l,ElemToPML(ElemID))-Flux(9:32,PP_N-q,p,SideID)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)           = Ut( :,p,q,l,ElemID)           -Flux(1: 8,PP_N-p,PP_N-q,SideID)*L_hatMinus(l)
        U2t(:,p,q,l,ElemToPML(ElemID))=U2t(:,p,q,l,ElemToPML(ElemID))-Flux(9:32,PP_N-p,PP_N-q,SideID)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)           = Ut( :,p,q,l,ElemID)           -Flux(1: 8,q,PP_N-p,SideID)*L_hatMinus(l)
        U2t(:,p,q,l,ElemToPML(ElemID))=U2t(:,p,q,l,ElemToPML(ElemID))-Flux(9:32,q,PP_N-p,SideID)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    END SELECT
!===================================================================================================================================
  CASE(XI_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut( :,l,p,q,ElemID)           =Ut( :,l,p,q,ElemID)           +Flux(1: 8,p,q,SideID)*L_hatPlus(l)
        U2t(:,l,p,q,ElemToPML(ElemID))=U2t(:,l,p,q,ElemToPML(ElemID))+Flux(9:32,p,q,SideID)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut( :,l,p,q,ElemID)           =Ut( :,l,p,q,ElemID)           -Flux(1: 8,q,p,SideID)*L_hatPlus(l)
        U2t(:,l,p,q,ElemToPML(ElemID))=U2t(:,l,p,q,ElemToPML(ElemID))-Flux(9:32,q,p,SideID)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut( :,l,p,q,ElemID)           =Ut( :,l,p,q,ElemID)           -Flux(1: 8,PP_N-p,q,SideID)*L_hatPlus(l)
        U2t(:,l,p,q,ElemToPML(ElemID))=U2t(:,l,p,q,ElemToPML(ElemID))-Flux(9:32,PP_N-p,q,SideID)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut( :,l,p,q,ElemID)           =Ut( :,l,p,q,ElemID)           -Flux(1: 8,PP_N-q,PP_N-p,SideID)*L_hatPlus(l)
        U2t(:,l,p,q,ElemToPML(ElemID))=U2t(:,l,p,q,ElemToPML(ElemID))-Flux(9:32,PP_N-q,PP_N-p,SideID)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut( :,l,p,q,ElemID)           =Ut( :,l,p,q,ElemID)           -Flux(1: 8,p,PP_N-q,SideID)*L_hatPlus(l)
        U2t(:,l,p,q,ElemToPML(ElemID))=U2t(:,l,p,q,ElemToPML(ElemID))-Flux(9:32,p,PP_N-q,SideID)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    END SELECT
!===================================================================================================================================
  CASE(ETA_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut( :,p,l,q,ElemID)           =Ut( :,p,l,q,ElemID)           +Flux(1: 8,PP_N-p,q,SideID)*L_hatPlus(l)
        U2t(:,p,l,q,ElemToPML(ElemID))=U2t(:,p,l,q,ElemToPML(ElemID))+Flux(9:32,PP_N-p,q,SideID)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut( :,p,l,q,ElemID)           =Ut( :,p,l,q,ElemID)           -Flux(1: 8,q,PP_N-p,SideID)*L_hatPlus(l)
        U2t(:,p,l,q,ElemToPML(ElemID))=U2t(:,p,l,q,ElemToPML(ElemID))-Flux(9:32,q,PP_N-p,SideID)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut( :,p,l,q,ElemID)           =Ut( :,p,l,q,ElemID)           -Flux(1: 8,p,q,SideID)*L_hatPlus(l)
        U2t(:,p,l,q,ElemToPML(ElemID))=U2t(:,p,l,q,ElemToPML(ElemID))-Flux(9:32,p,q,SideID)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut( :,p,l,q,ElemID)           =Ut( :,p,l,q,ElemID)           -Flux(1: 8,PP_N-q,p,SideID)*L_hatPlus(l)
        U2t(:,p,l,q,ElemToPML(ElemID))=U2t(:,p,l,q,ElemToPML(ElemID))-Flux(9:32,PP_N-q,p,SideID)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut( :,p,l,q,ElemID)           =Ut( :,p,l,q,ElemID)           -Flux(1: 8,PP_N-p,PP_N-q,SideID)*L_hatPlus(l)
        U2t(:,p,l,q,ElemToPML(ElemID))=U2t(:,p,l,q,ElemToPML(ElemID))-Flux(9:32,PP_N-p,PP_N-q,SideID)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    END SELECT
!===================================================================================================================================
  CASE(ZETA_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut( :,p,q,l,ElemID)           =Ut( :,p,q,l,ElemID)           +Flux(1: 8,p,q,SideID)*L_hatPlus(l)
        U2t(:,p,q,l,ElemToPML(ElemID))=U2t(:,p,q,l,ElemToPML(ElemID))+Flux(9:32,p,q,SideID)*L_hatPlus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(1) ! slave side, SideID=q,jSide=p
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut( :,p,q,l,ElemID)           =Ut( :,p,q,l,ElemID)           -Flux(1: 8,q,p,SideID)*L_hatPlus(l)
        U2t(:,p,q,l,ElemToPML(ElemID))=U2t(:,p,q,l,ElemToPML(ElemID))-Flux(9:32,q,p,SideID)*L_hatPlus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut( :,p,q,l,ElemID)           =Ut( :,p,q,l,ElemID)           -Flux(1: 8,PP_N-p,q,SideID)*L_hatPlus(l)
        U2t(:,p,q,l,ElemToPML(ElemID))=U2t(:,p,q,l,ElemToPML(ElemID))-Flux(9:32,PP_N-p,q,SideID)*L_hatPlus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut( :,p,q,l,ElemID)           =Ut( :,p,q,l,ElemID)           -Flux(1: 8,PP_N-q,PP_N-p,SideID)*L_hatPlus(l)
        U2t(:,p,q,l,ElemToPML(ElemID))=U2t(:,p,q,l,ElemToPML(ElemID))-Flux(9:32,PP_N-q,PP_N-p,SideID)*L_hatPlus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut( :,p,q,l,ElemID)           =Ut( :,p,q,l,ElemID)           -Flux(1: 8,p,PP_N-q,SideID)*L_hatPlus(l)
        U2t(:,p,q,l,ElemToPML(ElemID))=U2t(:,p,q,l,ElemToPML(ElemID))-Flux(9:32,p,PP_N-q,SideID)*L_hatPlus(l)
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
        Ut (:,0,p,q,ElemID)           =Ut (:,0,p,q,          ElemID )+Flux(1: 8,q,p,SideID)*L_hatMinus0
        U2t(:,0,p,q,ElemToPML(ElemID))=U2t(:,0,p,q,ElemToPML(ElemID))+Flux(9:32,q,p,SideID)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,0,p,q,ElemID)           =Ut (:,0,p,q,ElemID)           -Flux(1:8  ,p,q,SideID)*L_hatMinus0
        U2t(:,0,p,q,ElemToPML(ElemID))=U2t(:,0,p,q,ElemToPML(ElemID))-Flux(9:32,p,q,SideID)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        Ut( :,0,p,q,ElemID)           =Ut (:,0,p,q,ElemID)           -Flux(1:8 ,PP_N-q,p,SideID)*L_hatMinus0
        U2t(:,0,p,q,ElemToPML(ElemID))=U2t(:,0,p,q,ElemToPML(ElemID))-Flux(9:32,PP_N-q,p,SideID)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,0,p,q,ElemID)           =Ut (:,0,p,q,ElemID)           -Flux(1:8 ,PP_N-p,PP_N-q,SideID)*L_hatMinus0
        U2t(:,0,p,q,ElemToPML(ElemID))=U2t(:,0,p,q,ElemToPML(ElemID))-Flux(9:32,PP_N-p,PP_N-q,SideID)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,0,p,q,ElemID)           =Ut (:,0,p,q,ElemID)           -Flux(1:8 ,q,PP_N-p,SideID)*L_hatMinus0
        U2t(:,0,p,q,ElemToPML(ElemID))=U2t(:,0,p,q,ElemToPML(ElemID))-Flux(9:32,q,PP_N-p,SideID)*L_hatMinus0
      END DO; END DO ! p,q
    END SELECT
  
  ! switch to right hand system for ETA_PLUS direction
!===================================================================================================================================
  CASE(ETA_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,p,0,q,ElemID)           =Ut (:,p,0,q,ElemID)           +Flux(1:8 ,p,q,SideID)*L_hatMinus0
        U2t(:,p,0,q,ElemToPML(ElemID))=U2t(:,p,0,q,ElemToPML(ElemID))+Flux(9:32,p,q,SideID)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,p,0,q,ElemID)           =Ut (:,p,0,q,ElemID)           -Flux(1:8 ,q,p,SideID)*L_hatMinus0
        U2t(:,p,0,q,ElemToPML(ElemID))=U2t(:,p,0,q,ElemToPML(ElemID))-Flux(9:32,q,p,SideID)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,p,0,q,ElemID)           =Ut (:,p,0,q,ElemID)           -Flux(1:8 ,PP_N-p,q,SideID)*L_hatMinus0
        U2t(:,p,0,q,ElemToPML(ElemID))=U2t(:,p,0,q,ElemToPML(ElemID))-Flux(9:32,PP_N-p,q,SideID)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,p,0,q,ElemID)           =Ut (:,p,0,q,ElemID)           -Flux(1:8 ,PP_N-q,PP_N-p,SideID)*L_hatMinus0
        U2t(:,p,0,q,ElemToPML(ElemID))=U2t(:,p,0,q,ElemToPML(ElemID))-Flux(9:32,PP_N-q,PP_N-p,SideID)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,p,0,q,ElemID)           =Ut (:,p,0,q,ElemID)           -Flux(1:8 ,p,PP_N-q,SideID)*L_hatMinus0
        U2t(:,p,0,q,ElemToPML(ElemID))=U2t(:,p,0,q,ElemToPML(ElemID))-Flux(9:32,p,PP_N-q,SideID)*L_hatMinus0
      END DO; END DO ! p,q
    END SELECT
  
  ! switch to right hand system for ZETA_MINUS direction
!===================================================================================================================================
  CASE(ZETA_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,p,q,0,ElemID)           =Ut (:,p,q,0,ElemID)           +Flux(1:8 ,q,p,SideID)*L_hatMinus0
        U2t(:,p,q,0,ElemToPML(ElemID))=U2t(:,p,q,0,ElemToPML(ElemID))+Flux(9:32,q,p,SideID)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,p,q,0,ElemID)           =Ut (:,p,q,0,ElemID)           -Flux(1:8 ,p,q,SideID)*L_hatMinus0
        U2t(:,p,q,0,ElemToPML(ElemID))=U2t(:,p,q,0,ElemToPML(ElemID))-Flux(9:32,p,q,SideID)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,p,q,0,ElemID)           =Ut (:,p,q,0,ElemID)           -Flux(1:8 ,PP_N-q,p,SideID)*L_hatMinus0
        U2t(:,p,q,0,ElemToPML(ElemID))=U2t(:,p,q,0,ElemToPML(ElemID))-Flux(9:32,PP_N-q,p,SideID)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,p,q,0,ElemID)           =Ut (:,p,q,0,ElemID)           -Flux(1:8 ,PP_N-p,PP_N-q,SideID)*L_hatMinus0
        U2t(:,p,q,0,ElemToPML(ElemID))=U2t(:,p,q,0,ElemToPML(ElemID))-Flux(9:32,PP_N-p,PP_N-q,SideID)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,p,q,0,ElemID)           =Ut (:,p,q,0,ElemID)           -Flux(1:8 ,q,PP_N-p,SideID)*L_hatMinus0
        U2t(:,p,q,0,ElemToPML(ElemID))=U2t(:,p,q,0,ElemToPML(ElemID))-Flux(9:32,q,PP_N-p,SideID)*L_hatMinus0
      END DO; END DO ! p,q
    END SELECT
  
!===================================================================================================================================
  CASE(XI_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,PP_N,p,q,ElemID)           =Ut (:,PP_N,p,q,ElemID)           +Flux(1:8 ,p,q,SideID)*L_hatPlusN
        U2t(:,PP_N,p,q,ElemToPML(ElemID))=U2t(:,PP_N,p,q,ElemToPML(ElemID))+Flux(9:32,p,q,SideID)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,PP_N,p,q,ElemID)           =Ut (:,PP_N,p,q,ElemID)           -Flux(1:8 ,q,p,SideID)*L_hatPlusN
        U2t(:,PP_N,p,q,ElemToPML(ElemID))=U2t(:,PP_N,p,q,ElemToPML(ElemID))-Flux(9:32,q,p,SideID)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,PP_N,p,q,ElemID)           =Ut (:,PP_N,p,q,ElemID)           -Flux(1:8 ,PP_N-p,q,SideID)*L_hatPlusN
        U2t(:,PP_N,p,q,ElemToPML(ElemID))=U2t(:,PP_N,p,q,ElemToPML(ElemID))-Flux(9:32,PP_N-p,q,SideID)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,PP_N,p,q,ElemID)           =Ut (:,PP_N,p,q,ElemID)           -Flux(1:8 ,PP_N-q,PP_N-p,SideID)*L_hatPlusN
        U2t(:,PP_N,p,q,ElemToPML(ElemID))=U2t(:,PP_N,p,q,ElemToPML(ElemID))-Flux(9:32,PP_N-q,PP_N-p,SideID)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,PP_N,p,q,ElemID)           =Ut (:,PP_N,p,q,ElemID)           -Flux(1:8 ,p,PP_N-q,SideID)*L_hatPlusN
        U2t(:,PP_N,p,q,ElemToPML(ElemID))=U2t(:,PP_N,p,q,ElemToPML(ElemID))-Flux(9:32,p,PP_N-q,SideID)*L_hatPlusN
      END DO; END DO ! p,q
    END SELECT
  
  ! switch to right hand system for ETA_PLUS direction
!===================================================================================================================================
  CASE(ETA_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,p,PP_N,q,ElemID)           =Ut (:,p,PP_N,q,ElemID)           +Flux(1:8 ,PP_N-p,q,SideID)*L_hatPlusN
        U2t(:,p,PP_N,q,ElemToPML(ElemID))=U2t(:,p,PP_N,q,ElemToPML(ElemID))+Flux(9:32,PP_N-p,q,SideID)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,p,PP_N,q,ElemID)           =Ut (:,p,PP_N,q,ElemID)           -Flux(1:8 ,q,PP_N-p,SideID)*L_hatPlusN
        U2t(:,p,PP_N,q,ElemToPML(ElemID))=U2t(:,p,PP_N,q,ElemToPML(ElemID))-Flux(9:32,q,PP_N-p,SideID)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,p,PP_N,q,ElemID)           =Ut (:,p,PP_N,q,ElemID)           -Flux(1:8 ,p,q,SideID)*L_hatPlusN
        U2t(:,p,PP_N,q,ElemToPML(ElemID))=U2t(:,p,PP_N,q,ElemToPML(ElemID))-Flux(9:32,p,q,SideID)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,p,PP_N,q,ElemID)           =Ut (:,p,PP_N,q,ElemID)           -Flux(1:8 ,PP_N-q,p,SideID)*L_hatPlusN
        U2t(:,p,PP_N,q,ElemToPML(ElemID))=U2t(:,p,PP_N,q,ElemToPML(ElemID))-Flux(9:32,PP_N-q,p,SideID)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,p,PP_N,q,ElemID)           =Ut (:,p,PP_N,q,ElemID)           -Flux(1:8 ,PP_N-p,PP_N-q,SideID)*L_hatPlusN
        U2t(:,p,PP_N,q,ElemToPML(ElemID))=U2t(:,p,PP_N,q,ElemToPML(ElemID))-Flux(9:32,PP_N-p,PP_N-q,SideID)*L_hatPlusN
      END DO; END DO ! p,q
    END SELECT

  ! switch to right hand system for ZETA_MINUS direction
!===================================================================================================================================
  CASE(ZETA_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,p,q,PP_N,ElemID)           =Ut (:,p,q,PP_N,ElemID)           +Flux(1:8 ,p,q,SideID)*L_hatPlusN
        U2t(:,p,q,PP_N,ElemToPML(ElemID))=U2t(:,p,q,PP_N,ElemToPML(ElemID))+Flux(9:32,p,q,SideID)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,p,q,PP_N,ElemID)           =Ut (:,p,q,PP_N,ElemID)           -Flux(1:8 ,q,p,SideID)*L_hatPlusN
        U2t(:,p,q,PP_N,ElemToPML(ElemID))=U2t(:,p,q,PP_N,ElemToPML(ElemID))-Flux(9:32,q,p,SideID)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,p,q,PP_N,ElemID)           =Ut (:,p,q,PP_N,ElemID)           -Flux(1:8 ,PP_N-p,q,SideID)*L_hatPlusN
        U2t(:,p,q,PP_N,ElemToPML(ElemID))=U2t(:,p,q,PP_N,ElemToPML(ElemID))-Flux(9:32,PP_N-p,q,SideID)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,p,q,PP_N,ElemID)           =Ut (:,p,q,PP_N,ElemID)           -Flux(1:8 ,PP_N-q,PP_N-p,SideID)*L_hatPlusN
        U2t(:,p,q,PP_N,ElemToPML(ElemID))=U2t(:,p,q,PP_N,ElemToPML(ElemID))-Flux(9:32,PP_N-q,PP_N-p,SideID)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        Ut (:,p,q,PP_N,ElemID)           =Ut (:,p,q,PP_N,ElemID)           -Flux(1:8 ,p,PP_N-q,SideID)*L_hatPlusN
        U2t(:,p,q,PP_N,ElemToPML(ElemID))=U2t(:,p,q,PP_N,ElemToPML(ElemID))-Flux(9:32,p,PP_N-q,SideID)*L_hatPlusN
      END DO; END DO ! p,q
    END SELECT
  END SELECT !locSideID
#endif
END SUBROUTINE CalcSurfInt2PML

END MODULE MOD_SurfInt
