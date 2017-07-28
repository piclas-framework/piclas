MODULE MOD_Dielectric_Vars
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
! Dielectric region damping factor
LOGICAL             :: DoDielectric                   ! true/false switch for Dielectric calculation procedures
LOGICAL             :: DielectricInitIsDone           ! initialisation flag
LOGICAL,ALLOCATABLE :: isDielectricElem(:)            ! true if iElem is an element located within the Dielectric region
LOGICAL,ALLOCATABLE :: isDielectricFace(:)            ! true if iFace is a Face located wihtin or on the boarder (interface) of the
!                                                     ! Dielectric region
LOGICAL,ALLOCATABLE :: isDielectricInterFace(:)       ! true if iFace is a Face located on the boarder (interface) of the Dielectric
!                                                     ! region
INTEGER             :: DielectricwriteFields          ! output Eps field for debug
INTEGER             :: Dielectricspread               ! if true Eps_x=Eps_y=Eps_z for all Dielectric cells
INTEGER             :: DielectricprintInfo            ! 0=only root prints Dielectric info, 1=all procs print Dielectric info
INTEGER             :: DielectricprintInfoProcs       ! number of procs taking part in Dielectric info printing
REAL,DIMENSION(6)   :: xyzPhysicalMinMax              ! physical boundary coordinates, outside = Dielectric region
REAL,DIMENSION(6)   :: xyzDielectricMinMax            ! Dielectric boundary coordinates, outside = Dielectric region
LOGICAL             :: useDielectricMinMax            ! switch between 'xyzPhysicalMinMax' and 'xyzDielectricMinMax'
REAL                :: DielectricEps0                 ! damping constant for Dielectric region shift
REAL                :: DielectricMu0                  ! CFS-Dielectric aplha factor for complex frequency shift
! mapping variables
INTEGER             :: nDielectricElems,nDielectricFaces,nDielectricInterFaces          ! number of Dielectric elements and faces
!                                                                                       ! (mapping)
INTEGER,ALLOCATABLE :: DielectricToElem(:),DielectricToFace(:),DielectricInterToFace(:) ! mapping to total element/face list
INTEGER,ALLOCATABLE :: ElemToDielectric(:),FaceToDielectric(:),FaceTODielectricInter(:) ! mapping to Dielectric element/face list
!
REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:) :: DielectricEps
REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:) :: DielectricMu
REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:) :: DielectricEpsGlobal
REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:) :: DielectricMuGlobal
! gradients
!===================================================================================================================================
END MODULE MOD_Dielectric_Vars
