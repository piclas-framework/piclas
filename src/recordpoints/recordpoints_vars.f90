MODULE MOD_RecordPoints_Vars
!===================================================================================================================================
! Contains the parameters needed for the Navier Stokes calculation
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255) :: RPDefFile
LOGICAL            :: RecordPointsInitIsDone = .FALSE.
INTEGER            :: RP_OutputInterval
INTEGER            :: RP_SamplingOffset
LOGICAL            :: RP_inUse  = .FALSE.
LOGICAL            :: RP_onProc = .FALSE.
INTEGER            :: nRP_loc
INTEGER            :: nRP_global
INTEGER            :: offsetRP_loc
INTEGER(KIND=8)    :: iSample,iSample_lastWrite
INTEGER,ALLOCATABLE:: OffsetRP(:,:)
INTEGER,ALLOCATABLE:: RP_ElemID(:)
REAL,ALLOCATABLE   :: xi_RP(:,:)
REAL,ALLOCATABLE   :: L_xi_RP(:,:)
REAL,ALLOCATABLE   :: L_eta_RP(:,:)
REAL,ALLOCATABLE   :: L_zeta_RP(:,:)
REAL,ALLOCATABLE   :: u_RP(:,:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! MPI Communicator for RPs
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER            :: myRPrank
INTEGER            :: RP_COMM
INTEGER            :: nRP_Procs
!-----------------------------------------------------------------------------------------------------------------------------------
! Output Buffer
!-----------------------------------------------------------------------------------------------------------------------------------
TYPE tRPset
  REAL,ALLOCATABLE     :: data(:,:,:)
  INTEGER(KIND=8)      :: Offset 
  TYPE(tRPset),POINTER :: nextset
END TYPE tRPset

TYPE(tRPset),POINTER   :: firstset, actualset
!===================================================================================================================================

INTERFACE getNewRPset
  MODULE PROCEDURE getNewRPset
END INTERFACE getNewRPset

PUBLIC :: getNewRPset

CONTAINS 

SUBROUTINE getNewRPset(RPset,Offset)
!===================================================================================================================================
! Read RP parameters from ini file and RP definitions from HDF5 
!===================================================================================================================================
! MODULES
USE MOD_Preproc
IMPLICIT NONE
! INPUT VARIABLES
INTEGER(KIND=8),INTENT(IN)    :: Offset
! OUTPUT VARIABLES
TYPE(tRPset),POINTER          :: RPset
!===================================================================================================================================
ALLOCATE(RPset)
ALLOCATE(RPset%data(0:PP_nVar,1:nRP_loc,1:RP_OutputInterval))
RPset%data=0.
NULLIFY(RPset%nextset)
RPset%Offset=Offset
END SUBROUTINE getNewRPset

END MODULE MOD_recordPoints_Vars
