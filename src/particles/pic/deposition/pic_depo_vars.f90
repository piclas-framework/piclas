MODULE MOD_PICDepo_Vars
!===================================================================================================================================
! Contains the constant Advection Velocity Vector used for the linear scalar advection equation
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
LOGICAL                               :: DoSFEqui           ! use equidistant points for SF
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
!REAL,ALLOCATABLE                      :: Vdm_BernSteinN_GaussN(:,:)
!REAL,ALLOCATABLE                      :: sVdm_BernSteinN_GaussN(:,:)
!===================================================================================================================================
END MODULE MOD_PICDepo_Vars
