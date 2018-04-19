MODULE MOD_Interfaces_Vars
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
LOGICAL             :: InterfacesInitIsDone = .FALSE.
INTEGER,ALLOCATABLE :: InterfaceRiemann(:)            ! face identifier for switching between different Riemann solvers
!-----------------------------------------------------------------------------------------------------------------------------------
! GEOMETRY - PRE-DEFINED COORDINATES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                        :: GeometryIsSet=.FALSE.  ! use pre-defined coordinates, e.g., a gyrotron tube
REAL,ALLOCATABLE               :: Geometry(:,:)          ! Coordinates of the pre-defined geometry
INTEGER                        :: GeometryNPoints        ! Number of Coordinates (first array dimension of Geometry(:,:))
REAL,ALLOCATABLE               :: GeometryMin(:)         ! Minimum value of Geometry coordinates of each column
REAL,ALLOCATABLE               :: GeometryMax(:)         ! Maximum value of Geometry coordinates of each column
END MODULE MOD_Interfaces_Vars
