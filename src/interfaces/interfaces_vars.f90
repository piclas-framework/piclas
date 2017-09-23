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
END MODULE MOD_Interfaces_Vars
