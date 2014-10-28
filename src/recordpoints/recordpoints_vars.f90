MODULE MOD_RecordPoints_Vars
!===================================================================================================================================
! Variables needed for the evaluation of the record points
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255) :: RPDefFile               ! file with elementlocal parametric RP coords
LOGICAL            :: RecordPointsInitIsDone = .FALSE.
LOGICAL            :: RP_inUse  = .FALSE.
LOGICAL            :: RP_onProc = .FALSE.
LOGICAL            :: RP_fileExists = .FALSE. ! flag if RP file for analyze level has been created
INTEGER            :: RP_Buffersize           ! no. of time samples (size of RP_Data)
INTEGER            :: RP_MaxBuffersize        ! max. allowed no. of time samples
INTEGER            :: RP_SamplingOffset       ! sampling rate (each .. iterations)
INTEGER            :: nRP                     ! no. of RP on proc
INTEGER            :: nGlobalRP               ! total no. of RP
INTEGER            :: offsetRP                ! offset for each proc in global RP list
INTEGER            :: iSample=0               ! no of samples in array
INTEGER            :: nSamples=0              ! total no. samples in case of multiple io steps
INTEGER,ALLOCATABLE:: RP_ElemID(:)            ! mapping from RP->Elem (nRP)
REAL,ALLOCATABLE   :: L_xi_RP(:,:)            ! lagrange basis evaluated at RPs coords
REAL,ALLOCATABLE   :: L_eta_RP(:,:)            
REAL,ALLOCATABLE   :: L_zeta_RP(:,:)          
REAL,ALLOCATABLE   :: RP_Data(:,:,:)          ! solution evaluated at RPs (nvar,nRP,nSamples)
REAL,ALLOCATABLE   :: lastSample(:,:)         ! solution evaluated at RPs (nvar,nRP,nSamples)
CHARACTER(LEN=255) :: StrVarNames(PP_nVar)
!-----------------------------------------------------------------------------------------------------------------------------------
! MPI Communicator for RPs
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER            :: myRPrank
INTEGER            :: RP_COMM
INTEGER            :: nRP_Procs

END MODULE MOD_recordPoints_Vars
