!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "piclas.h"
!===================================================================================================================================
!> Contains the routines to exchange data using MPI-3 shared memory
!===================================================================================================================================
MODULE MOD_MPI_Shared
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE

#if USE_MPI
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE DefineParametersMPIShared
  MODULE PROCEDURE DefineParametersMPIShared
END INTERFACE

INTERFACE InitMPIShared
  MODULE PROCEDURE InitMPIShared
END INTERFACE

INTERFACE FinalizeMPIShared
  MODULE PROCEDURE FinalizeMPIShared
END INTERFACE

#ifdef DEBUG_MEMORY
INTERFACE Allocate_Shared_DEBUG
#else
INTERFACE Allocate_Shared
#endif /*DEBUG_MEMORY*/
  MODULE PROCEDURE Allocate_Shared_Logical_1
  MODULE PROCEDURE Allocate_Shared_Logical_2
  MODULE PROCEDURE Allocate_Shared_Int_1
  MODULE PROCEDURE Allocate_Shared_Int_2
  MODULE PROCEDURE Allocate_Shared_Int_3
  MODULE PROCEDURE Allocate_Shared_Int_4
  MODULE PROCEDURE Allocate_Shared_Real_1
  MODULE PROCEDURE Allocate_Shared_Real_2
  MODULE PROCEDURE Allocate_Shared_Real_3
  MODULE PROCEDURE Allocate_Shared_Real_4
  MODULE PROCEDURE Allocate_Shared_Real_5
  MODULE PROCEDURE Allocate_Shared_Real_6
END INTERFACE

PUBLIC::DefineParametersMPIShared
PUBLIC::InitMPIShared
PUBLIC::FinalizeMPIShared

#ifdef DEBUG_MEMORY
PUBLIC::Allocate_Shared_DEBUG
#else
PUBLIC::Allocate_Shared
#endif /*DEBUG_MEMORY*/

INTERFACE UNLOCK_AND_FREE_DUMMY
  MODULE PROCEDURE UNLOCK_AND_FREE_1
END INTERFACE

INTERFACE BARRIER_AND_SYNC
  MODULE PROCEDURE BARRIER_AND_SYNC
END INTERFACE

INTERFACE MPI_SIZE
  MODULE PROCEDURE MPI_SIZE
END INTERFACE

PUBLIC::UNLOCK_AND_FREE_DUMMY
PUBLIC::BARRIER_AND_SYNC
PUBLIC::MPI_SIZE
!==================================================================================================================================
CONTAINS

!==================================================================================================================================
!> Define parameters for MPI-3 shared memory
!==================================================================================================================================
SUBROUTINE DefineParametersMPIShared()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection         ("MPI Shared")

END SUBROUTINE DefineParametersMPIShared


!==================================================================================================================================
!> Initialize for MPI-3 shared memory
!==================================================================================================================================
SUBROUTINE InitMPIShared()
! MODULES
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: color
!==================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT MPI SHARED COMMUNICATION ...'

! Save the global number of procs
nProcessors_Global = nProcessors

! Split the node communicator (shared memory) from the global communicator on physical processor or node level
#if USE_CORE_SPLIT
  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,myRank,0,MPI_COMM_SHARED,iError)
#else
  ! Note that using SharedMemoryMethod=OMPI_COMM_TYPE_CORE somehow does not work in every case (intel/amd processors)
  ! Also note that OMPI_COMM_TYPE_CORE is undefined when not using OpenMPI
  CALL MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD,SharedMemoryMethod,0,MPI_INFO_NULL,MPI_COMM_SHARED,IERROR)
#endif

! Find my rank on the shared communicator, comm size and proc name
CALL MPI_COMM_RANK(MPI_COMM_SHARED, myComputeNodeRank,IERROR)
CALL MPI_COMM_SIZE(MPI_COMM_SHARED, nComputeNodeProcessors,IERROR)

! MPI3 shared implementation currently only works with equal procs per node
IF (MOD(nProcessors_Global,nComputeNodeProcessors).NE.0) &
  CALL ABORT(__STAMP__,'MPI shared communication currently only supported with equal procs per node!')

IF (nProcessors_Global/nComputeNodeProcessors.EQ.1) THEN
  SWRITE(UNIT_stdOUt,'(A,I0,A,I0,A)') ' | Starting shared communication with ',nComputeNodeProcessors,' procs on ',1,' node'
ELSE
  SWRITE(UNIT_stdOUt,'(A,I0,A,I0,A)') ' | Starting shared communication with ',nComputeNodeProcessors,' procs on ',         &
                                                            nProcessors_Global/nComputeNodeProcessors,' nodes'
END IF

! Send rank of compute node root to all procs on shared comm
IF (myComputeNodeRank.EQ.0) ComputeNodeRootRank = myRank
CALL MPI_BCAST(ComputeNodeRootRank,1,MPI_INTEGER,0,MPI_COMM_SHARED,IERROR)

! now split global communicator into small group leaders and the others
MPI_COMM_LEADERS_SHARED=MPI_COMM_NULL
myLeaderGroupRank=-1
color = MERGE(101,MPI_UNDEFINED,myComputeNodeRank.EQ.0)
CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,color,0,MPI_COMM_LEADERS_SHARED,IERROR)
IF(myComputeNodeRank.EQ.0)THEN
  CALL MPI_COMM_RANK(MPI_COMM_LEADERS_SHARED,myLeaderGroupRank,IERROR)
  CALL MPI_COMM_SIZE(MPI_COMM_LEADERS_SHARED,nLeaderGroupProcs,IERROR)
END IF

! Leaders inform every processor on their node about the leader group rank and group size
CALL MPI_BCAST(myLeaderGroupRank,1,MPI_INTEGER,0,MPI_COMM_SHARED,IERROR)
CALL MPI_BCAST(nLeaderGroupProcs,1,MPI_INTEGER,0,MPI_COMM_SHARED,IERROR)

! Create MPI_Info for shared memory windows
CALL MPI_INFO_CREATE(MPI_INFO_SHARED_LOOSE,IERROR)
CALL MPI_INFO_SET(   MPI_INFO_SHARED_LOOSE,'accumulate_ordering','none',IERROR)
! Only root allocates, size differs between ranks
!CALL MPI_INFO_SET(   MPI_INFO_SHARED_LOOSE,'same_size'          ,'true',IERROR)
CALL MPI_INFO_SET(   MPI_INFO_SHARED_LOOSE,'same_disp_unit'     ,'true',IERROR)

! synchronize everything or bad things will happen
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)

MPISharedInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')      ' INIT MPI SHARED COMMUNICATION DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitMPIShared


!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Logical_1(nVal,SM_WIN,DataPointer&
#ifdef DEBUG_MEMORY
        ,SM_WIN_NAME&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(1)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
LOGICAL,INTENT(OUT),POINTER               :: DataPointer(:)           !> Pointer to the RMA window
#ifdef DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
!==================================================================================================================================
! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#ifdef DEBUG_MEMORY
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(&
__STAMP__&
,'ERROR: Datapointer (Logical1) already associated')

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN,IERROR)

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR,IERROR)
END IF

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer,nVal)

END SUBROUTINE ALLOCATE_SHARED_LOGICAL_1

!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Logical_2(nVal,SM_WIN,DataPointer&
#ifdef DEBUG_MEMORY
        ,SM_WIN_NAME&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(2)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
LOGICAL,INTENT(OUT),POINTER               :: DataPointer(:,:)         !> Pointer to the RMA window
#ifdef DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
!==================================================================================================================================
! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#ifdef DEBUG_MEMORY
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(&
__STAMP__&
,'ERROR: Datapointer (Logical1) already associated')

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN,IERROR)

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR,IERROR)
END IF

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer,nVal)

END SUBROUTINE ALLOCATE_SHARED_LOGICAL_2

!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Int_1(nVal,SM_WIN,DataPointer&
#ifdef DEBUG_MEMORY
        ,SM_WIN_NAME&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(1)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
INTEGER,INTENT(OUT),POINTER               :: DataPointer(:)         !> Pointer to the RMA window
#ifdef DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
!==================================================================================================================================
! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#ifdef DEBUG_MEMORY
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(&
__STAMP__&
,'ERROR: Datapointer (Int1) already associated')

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN,IERROR)

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR,IERROR)
END IF

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer,nVal)

END SUBROUTINE ALLOCATE_SHARED_INT_1


!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Int_2(nVal,SM_WIN,DataPointer&
#ifdef DEBUG_MEMORY
        ,SM_WIN_NAME&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(2)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
INTEGER,INTENT(OUT),POINTER               :: DataPointer(:,:)         !> Pointer to the RMA window
#ifdef DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
!==================================================================================================================================
! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#ifdef DEBUG_MEMORY
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(&
__STAMP__&
,'ERROR: Datapointer (Int2) already associated')

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN,IERROR)

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR,IERROR)
END IF

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer,nVal)

END SUBROUTINE ALLOCATE_SHARED_INT_2


!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Int_3(nVal,SM_WIN,DataPointer&
#ifdef DEBUG_MEMORY
        ,SM_WIN_NAME&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(3)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
INTEGER,INTENT(OUT),POINTER               :: DataPointer(:,:,:)       !> Pointer to the RMA window
#ifdef DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
!==================================================================================================================================
! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#ifdef DEBUG_MEMORY
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(&
__STAMP__&
,'ERROR: Datapointer (Int3) already associated')

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN,IERROR)

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR,IERROR)
END IF

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer,nVal)

END SUBROUTINE ALLOCATE_SHARED_INT_3


!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Int_4(nVal,SM_WIN,DataPointer&
#ifdef DEBUG_MEMORY
        ,SM_WIN_NAME&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(4)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
INTEGER,INTENT(OUT),POINTER               :: DataPointer(:,:,:,:)     !> Pointer to the RMA window
#ifdef DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
!==================================================================================================================================
! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#ifdef DEBUG_MEMORY
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(&
__STAMP__&
,'ERROR: Datapointer (Int3) already associated')

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN,IERROR)

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR,IERROR)
END IF

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer,nVal)

END SUBROUTINE ALLOCATE_SHARED_INT_4


!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Real_1(nVal,SM_WIN,DataPointer&
#ifdef DEBUG_MEMORY
        ,SM_WIN_NAME&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(1)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
REAL   ,INTENT(OUT),POINTER               :: DataPointer(:)         !> Pointer to the RMA window
#ifdef DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
!==================================================================================================================================
! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#ifdef DEBUG_MEMORY
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(&
__STAMP__&
,'ERROR: Datapointer (Real1) already associated')

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN,IERROR)


! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR,IERROR)
END IF

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer,nVal)

END SUBROUTINE ALLOCATE_SHARED_REAL_1


!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Real_2(nVal,SM_WIN,DataPointer&
#ifdef DEBUG_MEMORY
        ,SM_WIN_NAME&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(2)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
REAL   ,INTENT(OUT),POINTER               :: DataPointer(:,:)         !> Pointer to the RMA window
#ifdef DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
!==================================================================================================================================
! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#ifdef DEBUG_MEMORY
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(&
__STAMP__&
,'ERROR: Datapointer (Real2) already associated')

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN,IERROR)

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR,IERROR)
END IF

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer,nVal)

END SUBROUTINE ALLOCATE_SHARED_REAL_2


!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Real_3(nVal,SM_WIN,DataPointer&
#ifdef DEBUG_MEMORY
        ,SM_WIN_NAME&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(3)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
REAL   ,INTENT(OUT),POINTER               :: DataPointer(:,:,:)       !> Pointer to the RMA window
#ifdef DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
!==================================================================================================================================
! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#ifdef DEBUG_MEMORY
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(&
__STAMP__&
,'ERROR: Datapointer (Real3) already associated')

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN,IERROR)

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR,IERROR)
END IF

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer,nVal)

END SUBROUTINE ALLOCATE_SHARED_REAL_3


!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Real_4(nVal,SM_WIN,DataPointer&
#ifdef DEBUG_MEMORY
        ,SM_WIN_NAME&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(4)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
REAL   ,INTENT(OUT),POINTER               :: DataPointer(:,:,:,:)     !> Pointer to the RMA window
#ifdef DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
!==================================================================================================================================
! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#ifdef DEBUG_MEMORY
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(&
__STAMP__&
,'ERROR: Datapointer (Real4) already associated')

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN,IERROR)

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR,IERROR)
END IF

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer,nVal)

END SUBROUTINE ALLOCATE_SHARED_REAL_4


!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Real_5(nVal,SM_WIN,DataPointer&
#ifdef DEBUG_MEMORY
        ,SM_WIN_NAME&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(5)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
REAL   ,INTENT(OUT),POINTER               :: DataPointer(:,:,:,:,:)   !> Pointer to the RMA window
#ifdef DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
!==================================================================================================================================
! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#ifdef DEBUG_MEMORY
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(&
__STAMP__&
,'ERROR: Datapointer (Real5) already associated')

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN,IERROR)

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR,IERROR)
END IF

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer,nVal)

END SUBROUTINE ALLOCATE_SHARED_REAL_5


!==================================================================================================================================
!> Allocate data with MPI-3 shared memory option
!==================================================================================================================================
SUBROUTINE Allocate_Shared_Real_6(nVal,SM_WIN,DataPointer&
#ifdef DEBUG_MEMORY
        ,SM_WIN_NAME&
#endif /*DEBUG_MEMORY*/
)
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: nVal(6)                  !> Local number of variables in each rank
INTEGER,INTENT(OUT)                       :: SM_WIN                   !> Shared memory window
REAL   ,INTENT(OUT),POINTER               :: DataPointer(:,:,:,:,:,:) !> Pointer to the RMA window
#ifdef DEBUG_MEMORY
CHARACTER(LEN=*),INTENT(IN)               :: SM_WIN_NAME              !> Shared memory window name
#endif /*DEBUG_MEMORY*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                               :: SM_PTR                   !> Base pointer, translated to DataPointer later
INTEGER                                   :: DISP_UNIT                !> Displacement unit
INTEGER(KIND=MPI_ADDRESS_KIND)            :: WIN_SIZE                 !> Size of the allocated memory window on current proc
!==================================================================================================================================
! Only node MPI root actually allocates the memory, all other nodes allocate memory with zero length but use the same displacement
WIN_SIZE  = MERGE(MPI_SIZE(PRODUCT(INT(nVal,KIND=8)),KIND(DataPointer)),INT(0,MPI_ADDRESS_KIND),myComputeNodeRank.EQ.0)
DISP_UNIT = 1

#ifdef DEBUG_MEMORY
LWRITE(UNIT_stdOut,'(A,I7,A65,I20)') "myrank=",myrank," Allocated "//TRIM(SM_WIN_NAME)//" with WIN_SIZE = ",WIN_SIZE
#endif /*DEBUG_MEMORY*/

IF (ASSOCIATED(DataPointer)) CALL abort(&
__STAMP__&
,'ERROR: Datapointer (Real6) already associated')

! Allocate MPI-3 remote memory access (RMA) type memory window
CALL MPI_WIN_ALLOCATE_SHARED(WIN_SIZE, DISP_UNIT, MPI_INFO_SHARED_LOOSE, MPI_COMM_SHARED, SM_PTR, SM_WIN,IERROR)

! Node MPI root already knows the location in virtual memory, all other find it here
IF (myComputeNodeRank.NE.0) THEN
  CALL MPI_WIN_SHARED_QUERY(SM_WIN, 0, WIN_SIZE, DISP_UNIT, SM_PTR,IERROR)
END IF

! SM_PTR can now be associated with a Fortran pointer and thus used to access the shared data
CALL C_F_POINTER(SM_PTR, DataPointer,nVal)

END SUBROUTINE ALLOCATE_SHARED_REAL_6


!==================================================================================================================================
!> Unlock and free shared memory array
!==================================================================================================================================
SUBROUTINE BARRIER_AND_SYNC(SharedWindow,Communicator) !,Barrier_Opt)
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(INOUT)       :: SharedWindow !> Shared memory window
INTEGER,INTENT(INOUT)       :: Communicator !> Shared memory communicator
! LOGICAL,INTENT(IN)          :: Barrier_Opt  !
! LOCAL VARIABLES
! LOGICAL                     :: Barrier
!==================================================================================================================================
! Barrier = MERGE(Barrier_Opt,.TRUE.,PRESENT(Barrier_Opt)

CALL MPI_WIN_SYNC(SharedWindow,iError)
! IF (Barrier) CALL MPI_BARRIER (Communicator,iError)
CALL MPI_BARRIER (Communicator,iError)
CALL MPI_WIN_SYNC(SharedWindow,iError)

! IF(iError.NE.0)THEN
!   CALL abort(__STAMP__,'ERROR in MPI_WIN_SYNC() for '//TRIM(SM_WIN_NAME)//': iError returned non-zero value =',IntInfoOpt=iError)
! END IF ! iError.NE.0
END SUBROUTINE BARRIER_AND_SYNC


!==================================================================================================================================
!> Unlock and free shared memory array
!==================================================================================================================================
SUBROUTINE UNLOCK_AND_FREE_1(SharedWindow,SM_WIN_NAME)
! MODULES
USE MOD_Globals
!USE MOD_ReadInTools
#ifdef DEBUG_MEMORY
USE MOD_MPI_Shared_Vars ,ONLY:myComputeNodeRank
#endif /*DEBUG_MEMORY*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(INOUT)       :: SharedWindow !> Shared memory window
CHARACTER(LEN=*),INTENT(IN) :: SM_WIN_NAME  !> Shared memory window name
! LOCAL VARIABLES
!==================================================================================================================================
#ifdef DEBUG_MEMORY
LWRITE(UNIT_stdOut,'(A,I7,A85)') "myrank=",myrank," Unlocking "//TRIM(SM_WIN_NAME)//" with MPI_WIN_UNLOCK_ALL()"
#endif /*DEBUG_MEMORY*/
CALL MPI_WIN_UNLOCK_ALL(SharedWindow,iError)
#ifdef DEBUG_MEMORY
LWRITE(UNIT_stdOut,'(A,I7,A85)') "myrank=",myrank," Freeing window "//TRIM(SM_WIN_NAME)//" with       MPI_WIN_FREE()"
#endif /*DEBUG_MEMORY*/
CALL MPI_WIN_FREE(      SharedWindow,iError)

IF(iError.NE.0)THEN
  CALL abort(__STAMP__,'ERROR in UNLOCK_AND_FREE_1() for '//TRIM(SM_WIN_NAME)//': iError returned non-zero value =',IntInfoOpt=iError)
END IF ! iError.NE.0
END SUBROUTINE UNLOCK_AND_FREE_1


!==================================================================================================================================
!> Finalize for MPI-3 shared memory
!==================================================================================================================================
SUBROUTINE FinalizeMPIShared()
! MODULES
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_MPI_Shared_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

! Free MPI_INFO objects
CALL MPI_INFO_FREE(MPI_INFO_SHARED_LOOSE,IERROR)

! Free the shared communicator
IF(MPI_COMM_SHARED        .NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(MPI_COMM_SHARED        ,IERROR)
IF(MPI_COMM_LEADERS_SHARED.NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(MPI_COMM_LEADERS_SHARED,IERROR)
MPISharedInitIsDone=.FALSE.

END SUBROUTINE FinalizeMPIShared


!===================================================================================================================================
!
!===================================================================================================================================
FUNCTION MPI_SIZE(nVal,VarSize)
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=8),INTENT(IN) :: nVal
INTEGER,INTENT(IN)         :: VarSize
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPI_SIZE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

IF (INT(nVal*INT(VarSize,KIND=8),KIND=8).LT.INT(HUGE(INT(1,KIND=MPI_ADDRESS_KIND)),KIND=8)) THEN
  MPI_SIZE = INT(nVal,KIND=MPI_ADDRESS_KIND) * INT(VarSize,KIND=MPI_ADDRESS_KIND)
ELSE
  MPI_SIZE = -1
  CALL ABORT(__STAMP__,'MPI_SIZE for shared array too large!')
END IF

END FUNCTION MPI_SIZE
#endif /* USE_MPI */
END MODULE MOD_MPI_Shared
