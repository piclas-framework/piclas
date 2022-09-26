!==================================================================================================================================
! Copyright (c) 2020 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_Mesh_Tools
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES (PUBLIC)
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE GetGlobalElemID
  PROCEDURE GetGlobalElemID
END INTERFACE

INTERFACE GetCNElemID
  PROCEDURE GetCNElemID
END INTERFACE

INTERFACE GetGlobalSideID
  PROCEDURE GetGlobalSideID
END INTERFACE

INTERFACE GetCNSideID
  PROCEDURE GetCNSideID
END INTERFACE

!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: GetGlobalElemID
PUBLIC :: GetCNElemID
PUBLIC :: GetGlobalSideID
PUBLIC :: GetCNSideID
#if USE_HDG
PUBLIC :: LambdaSideToMaster
#if USE_MPI
PUBLIC :: GetMasteriLocSides
#endif /*USE_MPI*/
#if USE_LOADBALANCE
PUBLIC :: BuildSideToNonUniqueGlobalSide
#endif /*USE_LOADBALANCE*/
#endif /*USE_HDG*/
!----------------------------------------------------------------------------------------------------------------------------------

ABSTRACT INTERFACE
  PPURE INTEGER FUNCTION GetGlobalElemIDInterface(iElem)
    INTEGER,INTENT(IN) :: iElem
  END FUNCTION
END INTERFACE

PROCEDURE(GetGlobalElemIDInterface),POINTER :: GetGlobalElemID    !< pointer defining the mapping: compute-node element ID -> global element ID

ABSTRACT INTERFACE
  PPURE INTEGER FUNCTION GetCNElemIDInterface(iElem)
    INTEGER,INTENT(IN) :: iElem
  END FUNCTION
END INTERFACE

PROCEDURE(GetCNElemIDInterface),POINTER     :: GetCNElemID        !< pointer defining the mapping: global element ID -> compute-node element ID

ABSTRACT INTERFACE
  PPURE INTEGER FUNCTION GetGlobalSideIDInterface(iSide)
    INTEGER,INTENT(IN) :: iSide
  END FUNCTION
END INTERFACE

PROCEDURE(GetGlobalSideIDInterface),POINTER :: GetGlobalSideID    !< pointer defining the mapping: compute-node element ID -> global element ID

ABSTRACT INTERFACE
  PPURE INTEGER FUNCTION GetCNSideIDInterface(iSide)
    INTEGER,INTENT(IN) :: iSide
  END FUNCTION
END INTERFACE

PROCEDURE(GetCNSideIDInterface),POINTER     :: GetCNSideID        !< pointer defining the mapping: global element ID -> compute-node element ID

! Initialization routines
INTERFACE InitGetGlobalElemID
  MODULE PROCEDURE InitGetGlobalElemID
END INTERFACE

INTERFACE InitGetCNElemID
  MODULE PROCEDURE InitGetCNElemID
END INTERFACE

INTERFACE InitGetGlobalSideID
  MODULE PROCEDURE InitGetGlobalSideID
END INTERFACE

INTERFACE InitGetCNSideID
  MODULE PROCEDURE InitGetCNSideID
END INTERFACE

PUBLIC::InitGetGlobalElemID
PUBLIC::InitGetCNElemID
PUBLIC::InitGetGlobalSideID
PUBLIC::InitGetCNSideID
!===================================================================================================================================
CONTAINS

!==================================================================================================================================!
!> Initialize GetGlobalElemID function (mapping of compute-node element ID to global element ID)
!==================================================================================================================================!
SUBROUTINE InitGetGlobalElemID()
! MODULES
#if USE_MPI
USE MOD_MPI_Shared_Vars, ONLY:nComputeNodeProcessors,nProcessors_Global
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: dummy
!==================================================================================================================================
#if USE_MPI
IF (nComputeNodeProcessors.EQ.nProcessors_Global) THEN
  GetGlobalElemID => GetGlobalElemID_iElem
ELSE
  GetGlobalElemID => GetGlobalElemID_fromTotalElem
END IF
#else
GetGlobalElemID => GetGlobalElemID_iElem
#endif

! Suppress compiler warning
RETURN
#if USE_MPI
dummy=GetGlobalElemID_fromTotalElem(1)
#endif
dummy=GetGlobalElemID_iElem(1)
END SUBROUTINE InitGetGlobalElemID


!==================================================================================================================================!
!> Initialize GetGlobalSideID function (mapping of compute-node side ID to global side ID)
!==================================================================================================================================!
SUBROUTINE InitGetGlobalSideID()
! MODULES
#if USE_MPI
USE MOD_MPI_Shared_Vars, ONLY:nComputeNodeProcessors,nProcessors_Global
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: dummy
!==================================================================================================================================
#if USE_MPI
IF (nComputeNodeProcessors.EQ.nProcessors_Global) THEN
  GetGlobalSideID => GetGlobalSideID_iSide
ELSE
  GetGlobalSideID => GetGlobalSideID_fromTotalSide
END IF
#else
GetGlobalSideID => GetGlobalSideID_iSide
#endif

! Suppress compiler warning
RETURN
#if USE_MPI
dummy=GetGlobalSideID_fromTotalSide(1)
#endif
dummy=GetGlobalSideID_iSide(1)
END SUBROUTINE InitGetGlobalSideID


!==================================================================================================================================!
!> Get the compute-node element ID in case of MPI=OFF or single compute node (CN)
!==================================================================================================================================!
PPURE INTEGER FUNCTION GetGlobalElemID_iElem(iElem)
! MODULES
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER :: GetGlobalElemID_iElem
!===================================================================================================================================
GetGlobalElemID_iElem = iElem
END FUNCTION GetGlobalElemID_iElem


!==================================================================================================================================!
!> Get the compute-node element ID in case of MPI=OFF or single compute node (CN)
!==================================================================================================================================!
PPURE FUNCTION GetGlobalSideID_iSide(iSide)
! MODULES
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iSide
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: GetGlobalSideID_iSide
!===================================================================================================================================
GetGlobalSideID_iSide = iSide
END FUNCTION GetGlobalSideID_iSide


#if USE_MPI
!==================================================================================================================================!
!> Get the global element ID in case of MPI=ON for single or multiple compute nodes (CN)
!==================================================================================================================================!
PPURE INTEGER FUNCTION GetGlobalElemID_fromTotalElem(iElem)
! MODULES
USE MOD_MPI_Shared_Vars, ONLY:CNTotalElem2GlobalElem
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!INTEGER :: GetGlobalElemID_fromTotalElem
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
GetGlobalElemID_fromTotalElem = CNTotalElem2GlobalElem(iElem)
END FUNCTION GetGlobalElemID_fromTotalElem
#endif /*USE_MPI*/


#if USE_MPI
!==================================================================================================================================!
!> Get the global element ID in case of MPI=ON for single or multiple compute nodes (CN)
!==================================================================================================================================!
PPURE INTEGER FUNCTION GetGlobalSideID_fromTotalSide(iSide)
! MODULES
USE MOD_MPI_Shared_Vars, ONLY:CNTotalSide2GlobalSide
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iSide
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!INTEGER :: GetGlobalSideID_fromTotalSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
GetGlobalSideID_fromTotalSide = CNTotalSide2GlobalSide(iSide)
END FUNCTION GetGlobalSideID_fromTotalSide
#endif /*USE_MPI*/


!==================================================================================================================================!
!> Initialize GetCNElemID function (mapping of global element ID to compute-node element ID)
!==================================================================================================================================!
SUBROUTINE InitGetCNElemID()
! MODULES
#if USE_MPI
USE MOD_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors,nProcessors_Global
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: dummy
!==================================================================================================================================
#if USE_MPI
IF (nComputeNodeProcessors.EQ.nProcessors_Global) THEN
  GetCNElemID => CNElemID_is_iElem
ELSE
  GetCNElemID => GetGlobalElem2CNTotalElem
END IF
#else
GetCNElemID => CNElemID_is_iElem
#endif

! Suppress compiler warning
RETURN
#if USE_MPI
dummy=GetGlobalElem2CNTotalElem(1)
#endif
dummy=CNElemID_is_iElem(1)
END SUBROUTINE InitGetCNElemID


!==================================================================================================================================!
!> Initialize GetCNSideID function (mapping of global element ID to compute-node element ID)
!==================================================================================================================================!
SUBROUTINE InitGetCNSideID()
! MODULES
#if USE_MPI
USE MOD_MPI_Shared_Vars ,ONLY: nComputeNodeProcessors,nProcessors_Global
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: dummy
!==================================================================================================================================
#if USE_MPI
IF (nComputeNodeProcessors.EQ.nProcessors_Global) THEN
  GetCNSideID => CNSideID_is_iSide
ELSE
  GetCNSideID => GetGlobalSide2CNTotalSide
END IF
#else
GetCNSideID => CNSideID_is_iSide
#endif

! Suppress compiler warning
RETURN
#if USE_MPI
dummy=GetGlobalSide2CNTotalSide(1)
#endif
dummy=CNSideID_is_iSide(1)
END SUBROUTINE InitGetCNSideID


!==================================================================================================================================!
!> Get the CN element ID in case of MPI=OFF or single compute node (CN)
!==================================================================================================================================!
PPURE INTEGER FUNCTION CNElemID_is_iElem(iElem)
! MODULES
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iElem ! Global and local element ID are the same
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER :: CNElemID_is_iElem
!===================================================================================================================================
CNElemID_is_iElem = iElem
END FUNCTION CNElemID_is_iElem


!==================================================================================================================================!
!> Get the CN element ID in case of MPI=OFF or single compute node (CN)
!==================================================================================================================================!
PPURE INTEGER FUNCTION CNSideID_is_iSide(iSide)
! MODULES
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iSide ! Global and local element ID are the same
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER :: CNSideID_is_iSide
!===================================================================================================================================
CNSideID_is_iSide = iSide
END FUNCTION CNSideID_is_iSide


#if USE_MPI
!==================================================================================================================================!
!> Get the CN element ID in case of MPI=ON for single or multiple compute nodes (CN)
!==================================================================================================================================!
PPURE INTEGER FUNCTION GetGlobalElem2CNTotalElem(iElem)
! MODULES
USE MOD_MPI_Shared_Vars, ONLY:GlobalElem2CNTotalElem
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iElem ! Global element ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!INTEGER :: GetGlobalElem2CNTotalElem
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
GetGlobalElem2CNTotalElem = GlobalElem2CNTotalElem(iElem)
END FUNCTION GetGlobalElem2CNTotalElem
#endif /*USE_MPI*/


#if USE_MPI
!==================================================================================================================================!
!> Get the CN element ID in case of MPI=ON for single or multiple compute nodes (CN)
!==================================================================================================================================!
PPURE INTEGER FUNCTION GetGlobalSide2CNTotalSide(iSide)
! MODULES
USE MOD_MPI_Shared_Vars, ONLY:GlobalSide2CNTotalSide
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iSide ! Global element ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!INTEGER :: GetGlobalSide2CNTotalSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
GetGlobalSide2CNTotalSide = GlobalSide2CNTotalSide(iSide)
END FUNCTION GetGlobalSide2CNTotalSide
#endif /*USE_MPI*/

#if USE_HDG
#if USE_MPI
!===================================================================================================================================
!> Setup array that contains the exchanged iLocSides from master to slaves: Send MINE, receive YOUR direction
!> The dimensions are dummy sized (too big, but then the send/recv from lambda can be used)
!===================================================================================================================================
SUBROUTINE GetMasteriLocSides()
! MODULES
USE MOD_PreProc
USE MOD_globals   ,ONLY: abort,MPI_COMM_WORLD
USE MOD_Mesh_Vars ,ONLY: MortarType,SideToElem,MortarInfo
USE MOD_Mesh_Vars ,ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_HDG_Vars  ,ONLY: nGP_face, iLocSides
USE MOD_Mesh_Vars ,ONLY: nSides
#if USE_MPI
USE MOD_MPI_Vars  ,ONLY: RecRequest_U,SendRequest_U
USE MOD_MPI       ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#endif /*USE_MPI*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iSide,SideID,iLocSide,iMortar,nMortars,MortarSideID
!===================================================================================================================================
! Exchange iLocSides from master to slaves: Send MINE, receive YOUR direction
SDEALLOCATE(iLocSides)
ALLOCATE(iLocSides(PP_nVar,nGP_face,nSides))
iLocSides = -100.
DO iSide = 1, nSides
  ! Get local side ID for each side ID
  iLocSides(:,:,iSide) = REAL(SideToElem(S2E_LOC_SIDE_ID,iSide))

  ! Small virtual mortar master side (blue) is encountered with MortarType(1,iSide) = 0
  ! Blue (small mortar master) side writes as yellow (big mortar master) for consistency (you spin me right round baby right round)
  IF(MortarType(1,iSide).EQ.0)THEN
    ! check all my big mortar sides and find the one to which the small virtual is connected
    Check1: DO MortarSideID=firstMortarInnerSide,lastMortarInnerSide
      nMortars=MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
      DO iMortar=1,nMortars
        SideID= MortarInfo(MI_SIDEID,iMortar,MortarType(2,MortarSideID)) !small SideID
        IF(iSide.EQ.SideID)THEN
          iLocSide = SideToElem(S2E_LOC_SIDE_ID,MortarSideID)
          IF(iLocSide.NE.-1)THEN ! MINE side (big mortar)
            iLocSides(:,:,iSide) = REAL(iLocSide)
          ELSE
            CALL abort(__STAMP__,'This big mortar side must be master')
          END IF !iLocSide.NE.-1
          EXIT Check1
        END IF ! iSide.EQ.SideID
      END DO !iMortar
    END DO Check1 !MortarSideID
  END IF ! MortarType(1,iSide).EQ.0
END DO

! At Mortar interfaces: Send my loc side ID (normal master or small mortar master sides) to the slave sides
CALL StartReceiveMPIData( 1 , iLocSides , 1 , nSides , RecRequest_U  , SendID=1 ) ! Receive YOUR
CALL StartSendMPIData(    1 , iLocSides , 1 , nSides , SendRequest_U , SendID=1 ) ! Send MINE
CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=1)
END SUBROUTINE GetMasteriLocSides
#endif /*USE_MPI*/


!===================================================================================================================================
!> Transform lambda solution from local coordinate system into master orientation for iSide and return as array "MasterSide"
!===================================================================================================================================
SUBROUTINE LambdaSideToMaster(iSide,MasterSide)
! MODULES
USE MOD_PreProc
USE MOD_globals   ,ONLY: abort
USE MOD_HDG_Vars  ,ONLY: lambda, nGP_face, iLocSides
USE MOD_Mesh_Vars ,ONLY: MortarType,SideToElem,MortarInfo
USE MOD_Mesh_Vars ,ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_Mappings  ,ONLY: CGNS_SideToVol2
USE MOD_Mesh_Vars ,ONLY: lastMPISide_MINE
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                           :: iSide
REAL,DIMENSION(PP_nVar,nGP_face),INTENT(OUT) :: MasterSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: p,q,r,rr,pq(1:2)
INTEGER          :: SideID,iLocSide,iMortar,nMortars,MortarSideID
!===================================================================================================================================

! Check if side is a master (iSide<=lastMPISide_MINE) or a slave (iSide>lastMPISide_MINE) side
! When slave sides are encountered, get the iLocSide ID from the neighbouring master, because later the orientation of the data
! is assumed to be in the master orientation
IF(iSide.GT.lastMPISide_MINE)THEN
  iLocSide = NINT(iLocSides(1,1,iSide))
ELSE
  iLocSide = SideToElem(S2E_LOC_SIDE_ID,iSide)
END IF ! iSide.GT.lastMPISide_MINE

!1 of 2: Master element with iLocSide = SideToElem(S2E_LOC_SIDE_ID,iSide)
IF(iLocSide.NE.-1)THEN ! MINE side
  DO q=0,PP_N
    DO p=0,PP_N
      pq=CGNS_SideToVol2(PP_N,p,q,iLocSide)
      r  = q    *(PP_N+1)+p    +1
      rr = pq(2)*(PP_N+1)+pq(1)+1
      MasterSide(:,r:r) = lambda(:,rr:rr,iSide)
    END DO
  END DO !p,q
  RETURN
END IF !iLocSide.NE.-1

! 2 of 2: Small virtual mortar master side (blue) is encountered with MortarType(1,iSide) = 0
! Blue (small mortar master) side writes as yellow (big mortar master) for consistency (you spin me right round baby right round)
IF(MortarType(1,iSide).EQ.0)THEN
  ! check all my big mortar sides and find the one to which the small virtual is connected
  Check2: DO MortarSideID=firstMortarInnerSide,lastMortarInnerSide
    nMortars=MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
    !locSide=MortarType(2,MortarSideID)
    DO iMortar=1,nMortars
      SideID= MortarInfo(MI_SIDEID,iMortar,MortarType(2,MortarSideID)) !small SideID
      IF(iSide.EQ.SideID)THEN
        iLocSide = SideToElem(S2E_LOC_SIDE_ID,MortarSideID)
        IF(iLocSide.NE.-1)THEN ! MINE side (big mortar)
          DO q=0,PP_N
            DO p=0,PP_N
              pq=CGNS_SideToVol2(PP_N,p,q,iLocSide)
              r  = q    *(PP_N+1)+p    +1
              rr = pq(2)*(PP_N+1)+pq(1)+1
              MasterSide(:,r:r) = lambda(:,rr:rr,iSide)
            END DO
          END DO !p,q
        ELSE
          CALL abort(__STAMP__,'This big mortar side must be master')
        END IF !iLocSide.NE.-1
        EXIT Check2
      END IF ! iSide.EQ.SideID
    END DO !iMortar
  END DO Check2 !MortarSideID
END IF ! MortarType(1,iSide).EQ.0

END SUBROUTINE LambdaSideToMaster


#if USE_LOADBALANCE
!===================================================================================================================================
!> Build mapping: Side index -> Non-unique global side index
!> Requires: Elem
!===================================================================================================================================
SUBROUTINE BuildSideToNonUniqueGlobalSide()
! MODULES
#if USE_DEBUG
USE MOD_Globals   ,ONLY: myrank,UNIT_StdOut,MPI_COMM_WORLD
#endif /*USE_DEBUG*/
USE MOD_Globals   ,ONLY: iError
USE MOD_Mesh_Vars ,ONLY: MortarType,ElemInfo,SideToNonUniqueGlobalSide,nSides,nElems,ElemToSide,offsetElem,MortarInfo
USE MOD_Mesh_Vars ,ONLY: GlobalUniqueSideID
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iElem
INTEGER :: NonUniqueGlobalSide,LocSideID,SideID,MortarSideID
INTEGER :: iMortar,nMortars
INTEGER :: locMortarSide
#if USE_DEBUG
INTEGER :: checkRank
#endif /*USE_DEBUG*/
!===================================================================================================================================
#if USE_DEBUG
checkRank = -1
#endif /*USE_DEBUG*/
! Mapping from side ID to on-unique global side ID
ALLOCATE(SideToNonUniqueGlobalSide(2,nSides))
SideToNonUniqueGlobalSide = -1
! Loop over all elements and check each of the 6 side, but also consider Mortar sides when a big Mortar is encountered
DO iElem=1,nElems
  ! Loop over all six sides (do not forget Mortars)
  nMortars = 0

  DO LocSideID=1,6
    ! get side ID
    SideID = ElemToSide(E2S_SIDE_ID,LocSideID,iElem)
    ! get non-unique global side ID
    NonUniqueGlobalSide = ElemInfo(3,iElem+offsetElem) + LocSideID + nMortars
    ! Note that Mortar-SideIDs will not be filled here
    IF(SideToNonUniqueGlobalSide(1,SideID).EQ.-1)THEN
      SideToNonUniqueGlobalSide(1,SideID) = NonUniqueGlobalSide
    ELSE
      SideToNonUniqueGlobalSide(2,SideID) = NonUniqueGlobalSide
    END IF ! SideToNonUniqueGlobalSide(SideID).EQ.-1
#if USE_DEBUG
        IF(myrank.eq.checkRank)THEN
        IPWRITE(UNIT_StdOut,*) &
            "LocSideID,SideID,GlobalUniqueSideID(SideID),NonUniqueGlobalSide,MortarType(1,SideID),nMortars                   ="&
           , LocSideID,SideID,GlobalUniqueSideID(SideID),NonUniqueGlobalSide,MortarType(1,SideID),nMortars
        END IF ! myrank.eq.0
#endif /*USE_DEBUG*/
    ! For Mortar sides, find the side IDs
    ! Exclude MortarType(1,SideID)==0, which corresponds to small mortar slave sides (on the same proc as the large))
    ! Note that MortarType(1,MortarSideID) is also 0 because these are the small mortar master sides
    IF(MortarType(1,SideID).GT.0)THEN
      ! Loop over all Mortar sides
      DO iMortar=1,MERGE(4,2,MortarType(1,SideID).EQ.1)
        locMortarSide = MortarType(2,SideID)
        ! get side ID
        MortarSideID = MortarInfo(MI_SIDEID,iMortar,locMortarSide)
        ! get non-unique global side ID
        NonUniqueGlobalSide = ElemInfo(3,iElem+offsetElem) + LocSideID + iMortar + nMortars
#if USE_DEBUG
        IF(myrank.eq.checkRank)THEN
          IPWRITE(UNIT_StdOut,*) &
              "LocSideID,MortarSideID,GlobalUniqueSideID(SideID),NonUniqueGlobalSide,MortarType(1,MortarSideID),iMortar,SideID =",&
          "           ?",MortarSideID,GlobalUniqueSideID(SideID),NonUniqueGlobalSide,MortarType(1,MortarSideID),iMortar,SideID
        END IF ! myrank.eq.0
#endif /*USE_DEBUG*/
        ! Only fill empty non-unique global sides
        IF(SideToNonUniqueGlobalSide(1,MortarSideID).EQ.-1)THEN
          SideToNonUniqueGlobalSide(1,MortarSideID) = NonUniqueGlobalSide
        ELSE
          SideToNonUniqueGlobalSide(2,MortarSideID) = NonUniqueGlobalSide
        END IF ! SideToNonUniqueGlobalSide(MortarSideID).EQ.-1
      END DO
      ! Get the number of Mortar sides, either 2 or 4 and count all encountered sides within one element
      nMortars = nMortars + MERGE(4,2,MortarType(1,SideID).EQ.1)
    END IF ! MortarType(1,SideID).eq.1
  END DO
END DO ! iElem
#if USE_DEBUG
IF(myrank.eq.0.AND.checkRank.GT.-1) read*; CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif /*USE_DEBUG*/
END SUBROUTINE BuildSideToNonUniqueGlobalSide
#endif /*USE_LOADBALANCE*/
#endif /*USE_HDG*/



END MODULE MOD_Mesh_Tools
