!==================================================================================================================================
! Copyright (c) 2020 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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
!----------------------------------------------------------------------------------------------------------------------------------

ABSTRACT INTERFACE
  PURE INTEGER FUNCTION GetGlobalElemIDInterface(iElem)
    INTEGER,INTENT(IN) :: iElem
  END FUNCTION
END INTERFACE

PROCEDURE(GetGlobalElemIDInterface),POINTER :: GetGlobalElemID    !< pointer defining the mapping: compute-node element ID -> global element ID

ABSTRACT INTERFACE
  PURE INTEGER FUNCTION GetCNElemIDInterface(iElem)
    INTEGER,INTENT(IN) :: iElem
  END FUNCTION
END INTERFACE

PROCEDURE(GetCNElemIDInterface),POINTER     :: GetCNElemID        !< pointer defining the mapping: global element ID -> compute-node element ID

ABSTRACT INTERFACE
  PURE INTEGER FUNCTION GetGlobalSideIDInterface(iSide)
    INTEGER,INTENT(IN) :: iSide
  END FUNCTION
END INTERFACE

PROCEDURE(GetGlobalSideIDInterface),POINTER :: GetGlobalSideID    !< pointer defining the mapping: compute-node element ID -> global element ID

ABSTRACT INTERFACE
  PURE INTEGER FUNCTION GetCNSideIDInterface(iSide)
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
PURE FUNCTION GetGlobalElemID_iElem(iElem)
! MODULES
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: GetGlobalElemID_iElem
!===================================================================================================================================
GetGlobalElemID_iElem = iElem
END FUNCTION GetGlobalElemID_iElem


!==================================================================================================================================!
!> Get the compute-node element ID in case of MPI=OFF or single compute node (CN)
!==================================================================================================================================!
PURE FUNCTION GetGlobalSideID_iSide(iSide)
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
PURE FUNCTION GetGlobalElemID_fromTotalElem(iElem)
! MODULES
USE MOD_MPI_Shared_Vars, ONLY:CNTotalElem2GlobalElem
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER :: GetGlobalElemID_fromTotalElem
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
PURE FUNCTION GetGlobalSideID_fromTotalSide(iSide)
! MODULES
USE MOD_MPI_Shared_Vars, ONLY:CNTotalSide2GlobalSide
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iSide
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER :: GetGlobalSideID_fromTotalSide
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
PURE FUNCTION CNElemID_is_iElem(iElem)
! MODULES
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iElem ! Global and local element ID are the same
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: CNElemID_is_iElem
!===================================================================================================================================
CNElemID_is_iElem = iElem
END FUNCTION CNElemID_is_iElem


!==================================================================================================================================!
!> Get the CN element ID in case of MPI=OFF or single compute node (CN)
!==================================================================================================================================!
PURE FUNCTION CNSideID_is_iSide(iSide)
! MODULES
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iSide ! Global and local element ID are the same
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: CNSideID_is_iSide
!===================================================================================================================================
CNSideID_is_iSide = iSide
END FUNCTION CNSideID_is_iSide


#if USE_MPI
!==================================================================================================================================!
!> Get the CN element ID in case of MPI=ON for single or multiple compute nodes (CN)
!==================================================================================================================================!
PURE FUNCTION GetGlobalElem2CNTotalElem(iElem)
! MODULES
USE MOD_MPI_Shared_Vars, ONLY:GlobalElem2CNTotalElem
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iElem ! Global element ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER :: GetGlobalElem2CNTotalElem
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
PURE FUNCTION GetGlobalSide2CNTotalSide(iSide)
! MODULES
USE MOD_MPI_Shared_Vars, ONLY:GlobalSide2CNTotalSide
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: iSide ! Global element ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER :: GetGlobalSide2CNTotalSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
GetGlobalSide2CNTotalSide = GlobalSide2CNTotalSide(iSide)
END FUNCTION GetGlobalSide2CNTotalSide
#endif /*USE_MPI*/

END MODULE MOD_Mesh_Tools
