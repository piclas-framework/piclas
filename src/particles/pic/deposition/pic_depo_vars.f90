MODULE MOD_PICDepo_Vars 
!===================================================================================================================================
! Contains the variables for the particle deposition
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE                      :: source(:,:,:,:,:)  ! source(1:4,PP_N,PP_N,PP_N,nElems)
REAL,ALLOCATABLE                      :: GaussBorder(:)     ! 1D coords of gauss points in -1|1 space
INTEGER,ALLOCATABLE                   :: GaussBGMIndex(:,:,:,:,:) ! Background mesh index of gausspoints (1:3,PP_N,PP_N,PP_N,nElems)
REAL,ALLOCATABLE                      :: GaussBGMFactor(:,:,:,:,:) ! BGM factor of gausspoints (1:3,PP_N,PP_N,PP_N,nElems)
REAL,ALLOCATABLE                      :: GPWeight(:,:,:,:,:,:,:) ! Weights for splines deposition (check pic_depo for details)
CHARACTER(LEN=256)                    :: DepositionType     ! Type of Deposition-Method
INTEGER,ALLOCATABLE                   :: PartToFIBGM(:,:)   ! Mapping form Particle to FIBGM
REAL,ALLOCATABLE                      :: ElemRadius2_SF(:)  ! elem radius plus radius_sf
REAL, ALLOCATABLE                     :: BGMSource(:,:,:,:)
REAL                                  :: r_sf               ! cutoff radius of shape function
REAL                                  :: r2_sf              ! cutoff radius of shape function * cutoff radius of shape function
REAL                                  :: r2_sf_inv          ! 1/cutoff radius of shape function * cutoff radius of shape function
REAL                                  :: w_sf               ! shapefuntion weight
REAL                                  :: r_sf0              ! minimal shape function radius
REAL                                  :: r_sf_scale         ! scaling of shape function radius
REAL                                  :: BetaFac            ! betafactor of shape-function || integral =1
INTEGER                               :: sf1d_dir           ! direction of 1D shape function 
INTEGER                               :: NDepo              ! polynomial degree of delta distri
REAL,ALLOCATABLE                      :: tempcharge(:)      ! temp-charge for epo. kernal
REAL,ALLOCATABLE                      :: NDepoChooseK(:,:)               ! array n over n
REAL,ALLOCATABLE                      :: wBaryNDepo(:)      ! barycentric weights for deposition
REAL,ALLOCATABLE                      :: swGPNDepo(:)       ! integration weights for deposition
REAL,ALLOCATABLE                      :: sJNDepo(:,:,:,:)   ! sj on ndepo
REAL,ALLOCATABLE                      :: XiNDepo(:)         ! gauss position of barycenters
REAL,ALLOCATABLE                      :: Vdm_NDepo_GaussN(:,:) ! VdM between different polynomial degrees
LOGICAL                               :: DoChangeBasis      ! Change polynomial degree
LOGICAL                               :: DoSFEqui           ! use equidistant points for SF
INTEGER                               :: SfRadiusInt        ! radius integer for cylindrical and spherical shape function
REAL,ALLOCATABLE                      :: ElemDepo_xGP(:,:,:,:,:)  ! element xGPs for deposition 
REAL,ALLOCATABLE                      :: Vdm_EquiN_GaussN(:,:)  ! Vdm from equidistant points to Gauss Points
INTEGER                               :: alpha_sf           ! shapefuntion exponent 
REAL                                  :: BGMdeltas(3)       ! Backgroundmesh size in x,y,z
REAL                                  :: FactorBGM(3)       ! Divider for BGM (to allow real numbers)
REAL                                  :: BGMVolume          ! Volume of a BGM Cell
INTEGER                               :: BGMminX            ! Local minimum BGM Index in x
INTEGER                               :: BGMminY            ! Local minimum BGM Index in y
INTEGER                               :: BGMminZ            ! Local minimum BGM Index in z
INTEGER                               :: BGMmaxX            ! Local maximum BGM Index in x
INTEGER                               :: BGMmaxY            ! Local maximum BGM Index in y
INTEGER                               :: BGMmaxZ            ! Local maximum BGM Index in z
LOGICAL                               :: Periodic_Depo      ! Flag for periodic treatment for deposition
INTEGER                               :: DeltaType          ! Flag
INTEGER                               :: NKnots
REAL,ALLOCATABLE                      :: Knots(:)
LOGICAL                               :: OutputSource       ! write the source to hdf5
REAL,ALLOCATABLE                      :: CellVolWeightFac(:)
INTEGER                               :: VolIntOrder
REAL,ALLOCATABLE                      :: VolInt_X(:)
REAL,ALLOCATABLE                      :: VolInt_W(:)
REAL,ALLOCATABLE                      :: CellVolWeight_Volumes(:,:,:,:)
INTEGER                               :: NbrOfSFdepoFixes                  !Number of fixes for shape func depo at planar BCs
REAL    , ALLOCATABLE                 :: SFdepoFixesGeo(:,:,:)             !1:nFixes;1:2(base,normal);1:3(x,y,z) normal outwards!!!
REAL    , ALLOCATABLE                 :: SFdepoFixesChargeMult(:)          !multiplier for mirrored charges (wall: -1.0, sym: 1.0)
!REAL,ALLOCATABLE                      :: Vdm_BernSteinN_GaussN(:,:)
!REAL,ALLOCATABLE                      :: sVdm_BernSteinN_GaussN(:,:)
!===================================================================================================================================
END MODULE MOD_PICDepo_Vars
