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
REAL,ALLOCATABLE  :: xGP(:)                      ! Gauss point coordinates
REAL,ALLOCATABLE  :: wGP(:)                      ! GP integration weights
REAL,ALLOCATABLE  :: swGP(:)                     ! 1.0/ GP integration weights
REAL,ALLOCATABLE  :: wBary(:)                    ! barycentric weights
REAL,ALLOCATABLE  :: wGPSurf(:,:)                ! wGPSurf(i,j)=wGP(i)*wGP(j)
REAL,ALLOCATABLE  :: NChooseK(:,:)                ! array n over n
CHARACTER(LEN=255)::StrNodeType
!===================================================================================================================================

LOGICAL           :: InterpolationInitIsDone = .FALSE.
END MODULE MOD_Interpolation_Vars
