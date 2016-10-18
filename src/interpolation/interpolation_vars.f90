MODULE MOD_Interpolation_Vars
!===================================================================================================================================
! Contains global variables used by the DG modules.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! reserved for Gauss Points with polynomial degree N, all allocated (0:N)
REAL,ALLOCATABLE  :: L_Plus(:), L_Minus(:)       ! L for boundary flux computation at both sides (-1,1)
REAL,ALLOCATABLE  :: L_PlusMinus(:,:)            ! L for boundary flux computation at both sides (-1,1)
REAL,ALLOCATABLE  :: xGP(:)                      ! Gauss point coordinates
REAL,ALLOCATABLE  :: wGP(:)                      ! GP integration weights
REAL,ALLOCATABLE  :: swGP(:)                     ! 1.0/ GP integration weights
REAL,ALLOCATABLE  :: wBary(:)                    ! barycentric weights
REAL,ALLOCATABLE  :: wGPSurf(:,:)                ! wGPSurf(i,j)=wGP(i)*wGP(j)
REAL,ALLOCATABLE  :: NChooseK(:,:)               ! array n over n
REAL,ALLOCATABLE  :: Vdm_Leg(:,:), sVdm_Leg(:,:) !< Legendre Vandermonde matrix
CHARACTER(LEN=255),PARAMETER :: NodeTypeG    = 'GAUSS'                    !< Gauss nodes (-1,1)
CHARACTER(LEN=255),PARAMETER :: NodeTypeGL   = 'GAUSS-LOBATTO'            !< Gauss-Lobatto nodes [-1,1]
CHARACTER(LEN=255),PARAMETER :: NodeTypeCL   = 'CHEBYSHEV-GAUSS-LOBATTO'  
CHARACTER(LEN=255),PARAMETER :: NodeTypeVISU = 'VISU'                     !< equidistant nodes [-1,1]
#if (PP_NodeType==1)
  CHARACTER(LEN=255),PARAMETER :: NodeType = 'GAUSS'
#elif (PP_NodeType==2)
  CHARACTER(LEN=255),PARAMETER :: NodeType = 'GAUSS-LOBATTO'
#elif (PP_NodeType==3)
  CHARACTER(LEN=255),PARAMETER :: NodeType = 'CHEBYSHEV-GAUSS-LOBATTO'
#endif
!===================================================================================================================================

LOGICAL           :: InterpolationInitIsDone = .FALSE.
END MODULE MOD_Interpolation_Vars
