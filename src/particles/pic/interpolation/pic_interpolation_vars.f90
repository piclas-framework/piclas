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
REAL,ALLOCATABLE                      :: FieldAtParticle(:,:) ! (PIC%maxParticleNumber,6) 2nd index: Ex,Ey,Ez,Bx,By,Bz
CHARACTER(LEN=256)                    :: InterpolationType    ! Type of Interpolation-Method
REAL                                  :: externalField(6)     ! ext field is added to the maxwell-solver-field
LOGICAL                               :: DoInterpolation      ! Flag for interpolation
LOGICAL                               :: useVTKFileEField      ! Flag for BGField via VTK-File 
LOGICAL                               :: useVTKFileBField      ! Flag for BGField via VTK-File 
! particle position interpolation method variables:
!REAL, ALLOCATABLE                     :: ElemPT(:,:,:)        ! (1:3,1:8,1:nElems) node coordinates in
                                                              ! transformed elem
REAL, ALLOCATABLE                     :: ElemT_inv(:,:,:)     ! (1:3,1:3,1:nElems) trafo matrix to transform 
                                                              ! physical coordinate into unit elem 
!REAL, ALLOCATABLE                     :: ElemK_inv(:,:,:)     ! (1:3,1:3,1:nElems) trafo matrix for linear guess 
!REAL, ALLOCATABLE                     :: ElemK1(:,:)        ! (1:3,1:nElems) vector required for linear guess

REAL, ALLOCATABLE                     :: BGEfieldAtNode(:,:)  ! (1:3,1:nNodes) BackGround-E-Field at nodes
REAL, ALLOCATABLE                     :: BGBfieldAtNode(:,:)  ! (1:3,1:nNodes) BackGround-B-Field at nodes
REAL                                  :: eps_distance = 0     ! epsillon distance accuracy of VTK BG Field ReadIns

LOGICAL                               :: usecurvedExternalField       ! use given external field. only for Bz variation in z
CHARACTER(LEN=40)                     :: FileNameCurvedExternalField  ! filename containing the externanl field csv tabe
REAL,ALLOCATABLE                      :: CurvedExternalField(:,:)     ! z - Pos , Bz
REAL                                  :: DeltaExternalField
INTEGER                               :: nIntPoints                   ! number of all interpolation points of curved external field
!===================================================================================================================================
END MODULE MOD_PICInterpolation_Vars
