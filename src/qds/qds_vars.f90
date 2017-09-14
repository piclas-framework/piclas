MODULE MOD_QDS_Vars
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                   :: QDSInitIsDone = .FALSE.
LOGICAL                   :: DoQDS                    ! true/false switch for QDS calculation procedures
REAL,ALLOCATABLE          :: UQDS(:,:,:,:,:)          ! U2( 1:QDSnVar,i,j,k,nQDSElems)
REAL,ALLOCATABLE          :: UQDSt(:,:,:,:,:)         ! U2t(1:QDSnVar,i,j,k,nQDSElems)
INTEGER                   :: QDSnVar                  ! is zero or 24 depending
INTEGER                   :: nQDSElems                ! number of QDS elements
REAL,ALLOCATABLE          :: UQDS_Master(:,:,:,:)
REAL,ALLOCATABLE          :: UQDS_Slave(:,:,:,:)
REAL,ALLOCATABLE          :: FluxQDS_Master(:,:,:,:)
REAL,ALLOCATABLE          :: FluxQDS_Slave(:,:,:,:)
REAL                      :: QDSSpeciesMass=9.1093826E-31 ! electron
REAL,ALLOCATABLE          :: QDSMacroValues(:,:,:,:,:)
REAL,ALLOCATABLE          :: GaussHermitWeiAbs(:,:)
REAL                      :: QDSSpecDOF
REAL                      :: QDSMaxVelo
CHARACTER(LEN=255),DIMENSION(40),PARAMETER :: StrVarNames(40)=(/ CHARACTER(LEN=255) :: 'PartMass1', &
                                                                                      'VeloX1', &
                                                                                      'VeloY1', &
                                                                                      'VeloZ1', &
                                                                                      'Energy1', &
                                                                                      'PartMass2', &
                                                                                      'VeloX2', &
                                                                                      'VeloY2', &
                                                                                      'VeloZ2', &
                                                                                      'Energy2', &
                                                                                      'PartMass3', &
                                                                                      'VeloX3', &
                                                                                      'VeloY3', &
                                                                                      'VeloZ3', &
                                                                                      'Energy3', &
                                                                                      'PartMass4', &
                                                                                      'VeloX4', &
                                                                                      'VeloY4', &
                                                                                      'VeloZ4', &
                                                                                      'Energy4', &
                                                                                      'PartMass5', &
                                                                                      'VeloX5', &
                                                                                      'VeloY5', &
                                                                                      'VeloZ5', &
                                                                                      'Energy5', &
                                                                                      'PartMass6', &
                                                                                      'VeloX6', &
                                                                                      'VeloY6', &
                                                                                      'VeloZ6', &
                                                                                      'Energy6', &
                                                                                      'PartMass7', &
                                                                                      'VeloX7', &
                                                                                      'VeloY7', &
                                                                                      'VeloZ7', &
                                                                                      'Energy7', &
                                                                                      'PartMass8', &
                                                                                      'VeloX8', &
                                                                                      'VeloY8', &
                                                                                      'VeloZ8', &
                                                                                      'Energy8'/)

END MODULE MOD_QDS_Vars
