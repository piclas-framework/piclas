#include "boltzplatz.h"

MODULE MOD_Interpolation
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MOD_basis
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part
! ----------------------------------------------------------------------------------------------------------------------------------
! Public Part
! ----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitInterpolation
   MODULE PROCEDURE InitInterpolation
END INTERFACE

INTERFACE FinalizeInterpolation
   MODULE PROCEDURE FinalizeInterpolation
END INTERFACE

INTERFACE ApplyJacobian
   MODULE PROCEDURE ApplyJacobian
END INTERFACE

PUBLIC::InitInterpolation
PUBLIC::ApplyJacobian
PUBLIC::FinalizeInterpolation


!===================================================================================================================================

CONTAINS


SUBROUTINE InitInterpolation()
!============================================================================================================================
! Initialize basis for Gauss-points of order N. 
! Prepares Differentiation matrices D, D_Hat, Basis at the boundaries L(1), L(-1), L_Hat(1), L_Hat(-1)
! Gaussnodes xGP and weights wGP, wBary 
!============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Interpolation_Vars
USE MOD_Restart_Vars,ONLY:DoRestart,N_Restart,BuildNewMesh,WriteNewMesh,InterpolateSolution
USE MOD_ReadInTools,ONLY:GETINT,GETREAL
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
!input parameters
!----------------------------------------------------------------------------------------------------------------------------
!output parameters
!----------------------------------------------------------------------------------------------------------------------------
!local variables
!============================================================================================================================
IF (InterpolationInitIsDone) THEN
  CALL abort(__STAMP__,'InitInterpolation already called.',999,999.)
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT INTERPOLATION...'
 
! Access ini-file
#if PP_N == N
PP_N=GETINT('N','2')   ! N could be set by readin_HDF5 routine -> postproctool
#endif
CALL initInterpolationBasis(PP_N)


InterpolationInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT INTERPOLATION DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitInterpolation


SUBROUTINE InitInterpolationBasis(N_in)
!============================================================================================================================
! Initialize basis for Gauss-points of order N. 
! Prepares Differentiation matrices D, D_Hat, Basis at the boundaries L(1), L(-1), L_Hat(1), L_Hat(-1)
! Gaussnodes xGP and weights wGP, wBary 
!============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Interpolation_Vars,ONLY:xGP,wGP,wBary,L_Minus,L_Plus,StrNodeType
USE MOD_Basis,ONLY:LegendreGaussNodesAndWeights,LegGaussLobNodesAndWeights,ChebyGaussLobNodesAndWeights
USE MOD_Basis,ONLY:BarycentricWeights,LagrangeInterpolationPolys
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                         :: N_in
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!============================================================================================================================
! Allocate global variables, needs to go somewhere else later
ALLOCATE(xGP(0:N_in), wGP(0:N_in), wBary(0:N_in))
ALLOCATE(L_Minus(0:N_in), L_Plus(0:N_in))

! Compute Nodes and weights for Gauss or GaussLobatto-Nodes
#if (PP_NodeType==1)
    StrNodeType='GAUSS'
    CALL LegendreGaussNodesAndWeights(N_in,xGP,wGP)
    SWRITE(UNIT_stdOut,'(A)') ' NodeType: GAUSS'
#elif (PP_NodeType==2)
    StrNodeType='GAUSS-LOBATTO'
    CALL LegGaussLobNodesAndWeights(N_in,xGP,wGP)
    SWRITE(UNIT_stdOut,'(A)') ' NodeType: GAUSS-LOBATTO'
#elif (PP_NodeType==3)
    StrNodeType='CHEBYSHEV-GAUSS-LOBATTO'
    CALL ChebyGaussLobNodesAndWeights(N_in,xGP,wGP)
    SWRITE(UNIT_stdOut,'(A)') ' NodeType: CHEBYSHEV-GAUSS-LOBATTO'
    SWRITE(UNIT_stdOut,'(A)') ' WARNING: CGL-nodes NOT fully implemented, Gauss-Nodes used, Enter to continue!'
    READ(*,*)
#else
    SWRITE(UNIT_stdOut,'(A)') ' NodeType NOT implemented'
    CALL abort(__STAMP__,'Code stopped!',999,999.)
#endif

CALL BarycentricWeights(N_in,xGP,wBary)

!! interpolate to left and right face (1 and -1) and pre-divide by mass matrix
CALL LagrangeInterpolationPolys(1.,N_in,xGP,wBary,L_Plus)
CALL LagrangeInterpolationPolys(-1.,N_in,xGP,wBary,L_Minus)
END SUBROUTINE InitInterpolationBasis

SUBROUTINE ApplyJacobian(U,toPhysical,toSwap)
!===================================================================================================================================
! Convert solution between physical <-> reference space
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,ONLY:sJ
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT) :: U(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems)
LOGICAL,INTENT(IN) :: toPhysical
LOGICAL,INTENT(IN) :: toSwap
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iElem,iVar
!===================================================================================================================================
IF(toPhysical)THEN
  IF(toSwap)THEN
    DO iElem=1,PP_nElems
      DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N; DO iVar=1,PP_nVar
          U(iVar,i,j,k,iElem)=-U(iVar,i,j,k,iElem)*sJ(i,j,k,iElem)
      END DO; END DO; END DO; END DO 
    END DO
  ELSE
    DO iElem=1,PP_nElems
      DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N; DO iVar=1,PP_nVar
        U(iVar,i,j,k,iElem)=U(iVar,i,j,k,iElem)*sJ(i,j,k,iElem)
      END DO; END DO; END DO; END DO
    END DO
  END IF
ELSE
  IF(toSwap)THEN
    DO iElem=1,PP_nElems
      DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N; DO iVar=1,PP_nVar
        U(iVar,i,j,k,iElem)=-U(iVar,i,j,k,iElem)/sJ(i,j,k,iElem)
      END DO; END DO; END DO; END DO
    END DO
  ELSE
    DO iElem=1,PP_nElems
      DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N; DO iVar=1,PP_nVar
        U(iVar,i,j,k,iElem)=U(iVar,i,j,k,iElem)/sJ(i,j,k,iElem)
      END DO; END DO; END DO; END DO
    END DO
  END IF
END IF
END SUBROUTINE ApplyJacobian

SUBROUTINE FinalizeInterpolation()
!============================================================================================================================
! Deallocate all global interpolation variables.
!============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Interpolation_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
!input parameters
!----------------------------------------------------------------------------------------------------------------------------
!output parameters
!----------------------------------------------------------------------------------------------------------------------------
!local variables
!============================================================================================================================
! Deallocate global variables, needs to go somewhere else later
SDEALLOCATE(xGP)
SDEALLOCATE(wGP)
SDEALLOCATE(wBary)
SDEALLOCATE(L_Minus)
SDEALLOCATE(L_Plus)

InterpolationInitIsDone = .FALSE.
END SUBROUTINE FinalizeInterpolation

END MODULE MOD_Interpolation

