MODULE MOD_PICInterpolation_Vars
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
REAL,ALLOCATABLE                    :: FieldAtParticle(:,:) ! (PIC%maxParticleNumber,6) 2nd index: Ex,Ey,Ez,Bx,By,Bz
CHARACTER(LEN=256)                    :: InterpolationType    ! Type of Interpolation-Method
REAL                                  :: externalField(6)     ! ext field is added to the maxwell-solver-field
LOGICAL                               :: DoInterpolation      ! Flag for interpolation
LOGICAL                               :: useBGField           ! Flag for BGField via h5-File
INTEGER                               :: NBG                  ! Polynomial degree of BG-Field
INTEGER                               :: BGType               ! Type of BG-Field (Electric,Magnetic,Both)
INTEGER                               :: BGDataSize           ! Type of BG-Field (Electric,Magnetic,Both)
REAL, ALLOCATABLE                     :: BGField(:,:,:,:,:)   ! BGField data
                                                              ! (1:x,0:NBG,0:NBG,0:NBG,1:PP_nElems)
REAL,ALLOCATABLE                      :: BGField_xGP(:)       ! Gauss point coordinates
REAL,ALLOCATABLE                      :: BGField_wGP(:)       ! GP integration weights
REAL,ALLOCATABLE                      :: BGField_wBary(:)     ! barycentric weights

                                                              

INTEGER                               :: Interpolation_p_IDW  ! power parameter for Inverse distance weighting
! particle position interpolation method variables:

REAL, ALLOCATABLE                     :: ElemT_inv(:,:,:)     ! (1:3,1:3,1:nElems) trafo matrix to transform 
                                                              ! physical coordinate into unit elem 



CHARACTER(LEN=256)                     :: FileNameCurvedExternalField  ! filename containing the externanl field csv tabe
LOGICAL                               :: usecurvedExternalField       ! use given external field. only for Bz variation in z
REAL,ALLOCATABLE                      :: CurvedExternalField(:,:)     ! z - Pos , Bz
REAL                                  :: DeltaExternalField
INTEGER                               :: nIntPoints                   ! number of all interpolation points of curved external field
!===================================================================================================================================
END MODULE MOD_PICInterpolation_Vars
