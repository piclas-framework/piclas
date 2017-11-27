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
LOGICAL             :: DielectricCheckRadius          ! instead of a bounding box region for setting a dielectric area, use radius
REAL                :: DielectricRadiusValue          ! radius for setting dielectric element ON/OFF
INTEGER             :: Dielectricspread               ! if true Eps_x=Eps_y=Eps_z for all Dielectric cells
REAL,DIMENSION(6)   :: xyzPhysicalMinMaxDielectric    ! physical   boundary coordinates, outside = Dielectric region
REAL,DIMENSION(6)   :: xyzDielectricMinMax            ! Dielectric boundary coordinates, outside = physical region
LOGICAL             :: useDielectricMinMax            ! switch between 'xyzPhysicalMinMax' and 'xyzDielectricMinMax'
CHARACTER(255)      :: DielectricTestCase             ! special test cases, e.g., "fish eye lens" Maxwell 1860
REAL                :: DielectricEpsR                 ! for Dielectric region shift
REAL                :: DielectricEpsR_inv             ! 1./EpsR
#ifdef PP_HDG
REAL                :: DielectricRatio                ! set dielectric ratio e_io = eps_inner/eps_outer for dielectric sphere
REAL                :: Dielectric_E_0                 ! axial electric field strength in x-direction of the dielectric sphere setup
#endif /*PP_HDG*/
REAL                :: DielectricMuR                  ! MuR
REAL                :: DielectricRmax                 ! maximum radius for dielectric material distribution
REAL                :: DielectricConstant_RootInv     ! 1./sqrt(EpsR*MuR)
REAL                :: eta_c_dielectric               ! ( chi - 1./sqrt(EpsR*MuR) ) * c
REAL                :: c_dielectric                   ! c/sqrt(EpsR*MuR)
REAL                :: c2_dielectric                  ! c**2/(EpsR*MuR)
! mapping variables
INTEGER             :: nDielectricElems,nDielectricFaces,nDielectricInterFaces          ! number of Dielectric elements and faces
!                                                                                       ! (mapping)
INTEGER,ALLOCATABLE :: DielectricToElem(:),DielectricToFace(:),DielectricInterToFace(:) ! mapping to total element/face list
INTEGER,ALLOCATABLE :: ElemToDielectric(:),FaceToDielectric(:),FaceToDielectricInter(:) ! mapping to Dielectric element/face list
!
REAL,ALLOCATABLE,DIMENSION(:,:,:,:)   :: DielectricEps
REAL,ALLOCATABLE,DIMENSION(:,:,:,:)   :: DielectricMu
REAL,ALLOCATABLE,DIMENSION(:,:,:,:)   :: DielectricConstant_inv         ! 1./(EpsR*MuR)
REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:) :: DielectricGlobal               ! contains DielectricEps and DielectricMu for HDF5 output
REAL,ALLOCATABLE,DIMENSION(:,:,:)     :: Dielectric_Master ! face array containing 1./SQRT(EpsR*MuR) for each DOF
REAL,ALLOCATABLE,DIMENSION(:,:,:)     :: Dielectric_Slave
! gradients
!===================================================================================================================================
END MODULE MOD_Dielectric_Vars
