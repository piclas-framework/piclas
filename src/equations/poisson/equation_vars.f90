MODULE MOD_Equation_Vars
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
REAL              :: Pi
REAL              :: IniWavenumber(3) ! wavenumbers in 3 directions (sinus periodic with exactfunc=6) 
INTEGER           :: IniExactFunc
REAL              :: IniCenter(3)
REAL              :: IniAmplitude 
REAL              :: IniHalfwidth 

! needed for various stuff (compilation)
REAL              :: c_corr
REAL              :: c_corr2    !c_corr^2
REAL              :: c_corr_c   !c_corr*c
REAL              :: c_corr_c2  !c_corr*c^2
REAL              :: eta_c      !(c_corr -1 )*c
REAL              :: fDamping          
LOGICAL           :: DoParabolicDamping



REAL,ALLOCATABLE  :: chitens(:,:,:,:,:,:) !diffusion 3x3 tensor on each gausspoint
REAL,ALLOCATABLE  :: chitensInv(:,:,:,:,:,:) ! inverse of diffusion 3x3 tensor on each gausspoint
REAL,ALLOCATABLE  :: chitens_face(:,:,:,:,:) !diffusion 3x3 tensor on each face gausspoint



LOGICAL           :: EquationInitIsDone=.FALSE.
REAL              :: eps0 
REAL              :: mu0, smu0 
REAL              :: c
REAL              :: c2
REAL              :: c2_inv
REAL              :: c_inv
INTEGER           :: alpha_shape
REAL              :: shapeFuncPrefix
REAL              :: rCutoff
REAL,ALLOCATABLE  :: E(:,:,:,:,:)
! Boundary condition arrays
REAL,ALLOCATABLE     :: BCData(:,:,:,:)
INTEGER,ALLOCATABLE  :: nBCByType(:)
INTEGER,ALLOCATABLE  :: BCSideID(:,:)
! can specify BC state
CHARACTER(LEN=255):: BCStateFile

CHARACTER(LEN=255),DIMENSION(4),PARAMETER :: StrVarNames(4)=(/ CHARACTER(LEN=255) :: 'Phi'           , &
                                                                                     'ElectricFieldX', &
                                                                                     'ElectricFieldY', &
                                                                                     'ElectricFieldZ'/)
!===================================================================================================================================
END MODULE MOD_Equation_Vars
