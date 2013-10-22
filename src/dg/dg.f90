#include "boltzplatz.h"

MODULE MOD_DG
!===================================================================================================================================
! Contains the initialization of the DG global variables
! Computes the different DG spatial operators/residuals(Ut) using U 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
INTERFACE FillIni
  MODULE PROCEDURE FillIni
END INTERFACE

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitDG
  MODULE PROCEDURE InitDG
END INTERFACE

INTERFACE DGTimeDerivative_weakForm
  MODULE PROCEDURE DGTimeDerivative_weakForm
END INTERFACE

INTERFACE FinalizeDG
  MODULE PROCEDURE FinalizeDG
END INTERFACE


PUBLIC::InitDG,DGTimeDerivative_weakForm,FinalizeDG,DGTimeDerivative_WoSource_weakForm
!===================================================================================================================================

CONTAINS

SUBROUTINE InitDG()
!===================================================================================================================================
! Allocate global variable U (solution) and Ut (dg time derivative).
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars
USE MOD_Restart_Vars,       ONLY: DoRestart,RestartInitIsDone
USE MOD_Interpolation_Vars, ONLY: xGP,wGP,wBary,InterpolationInitIsDone
USE MOD_Mesh_Vars,          ONLY: nSides,nInnerSides,nBCSides
USE MOD_Mesh_Vars,          ONLY: SideID_plus_lower,SideID_plus_upper
USE MOD_Mesh_Vars,          ONLY: SideID_minus_lower,SideID_minus_upper
USE MOD_Mesh_Vars,          ONLY: MeshInitIsDone 
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.(.NOT.MeshInitIsDone).OR.(.NOT.RestartInitIsDone).OR.DGInitIsDone)THEN
   CALL abort(__STAMP__,'InitDG not ready to be called or already called.',999,999.)
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT DG...'

CALL initDGbasis(PP_N,xGP,wGP,wBary)
! the local DG solution
ALLOCATE(U(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
! the time derivative computed with the DG scheme
ALLOCATE(Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
nTotalU=PP_nVar*(PP_N+1)*(PP_N+1)*(PP_N+1)*PP_nElems

IF(.NOT.DoRestart)THEN
  ! U is filled with the ini solution
  CALL FillIni()
END IF
! Ut is set to zero because it is successively updated with DG contributions
Ut=0.

! We store the interior data at the each element face
ALLOCATE(U_Minus(PP_nVar,0:PP_N,0:PP_N,sideID_minus_lower:sideID_minus_upper))
ALLOCATE(U_Plus(PP_nVar,0:PP_N,0:PP_N,sideID_plus_lower:sideID_plus_upper))
U_Minus=0.
U_Plus=0.

! unique flux per side
ALLOCATE(Flux(PP_nVar,0:PP_N,0:PP_N,1:nSides))
Flux=0.

DGInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT DG DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitDG


SUBROUTINE InitDGbasis(N_in,xGP,wGP,wBary)
!===================================================================================================================================
! Allocate global variable U (solution) and Ut (dg time derivative).
!===================================================================================================================================
! MODULES
USE MOD_Basis,ONLY:LegendreGaussNodesAndWeights,LegGaussLobNodesAndWeights,BarycentricWeights
USE MOD_Basis,ONLY:PolynomialDerivativeMatrix,LagrangeInterpolationPolys
USE MOD_DG_Vars,ONLY:D,D_Hat,L_HatMinus,L_HatPlus
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                         :: N_in
REAL,INTENT(IN),DIMENSION(0:N_in)          :: xGP,wGP,wBary
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL,DIMENSION(0:N_in,0:N_in)              :: M,Minv
REAL,DIMENSION(0:N_in)                     :: L_minus,L_plus        
INTEGER                                    :: iMass         
!===================================================================================================================================
ALLOCATE(D_Hat(0:N_in,0:N_in), L_HatMinus(0:N_in), L_HatPlus(0:N_in))
ALLOCATE(D(0:N_in,0:N_in))
! Compute Differentiation matrix D for given Gausspoints
CALL PolynomialDerivativeMatrix(N_in,xGP,D)

! Build D_Hat matrix. (D^ = M^(-1) * D^T * M
M(:,:)=0.
Minv(:,:)=0.
DO iMass=0,N_in
  M(iMass,iMass)=wGP(iMass)
  Minv(iMass,iMass)=1./wGP(iMass)
END DO
D_Hat(:,:) = -MATMUL(Minv,MATMUL(TRANSPOSE(D),M))

! interpolate to left and right face (1 and -1) and pre-divide by mass matrix
CALL LagrangeInterpolationPolys(1.,N_in,xGP,wBary,L_Plus)
L_HatPlus(:) = MATMUL(Minv,L_Plus)
CALL LagrangeInterpolationPolys(-1.,N_in,xGP,wBary,L_Minus)
L_HatMinus(:) = MATMUL(Minv,L_Minus)
END SUBROUTINE InitDGbasis

SUBROUTINE DGTimeDerivative_weakForm(t,tStage,tDeriv)
!===================================================================================================================================
! Computes the DG time derivative consisting of Volume Integral and Surface integral for the whole field
! U and Ut are allocated
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_DG_Vars,       ONLY: U,Ut,nTotalU,U_Plus,U_Minus,Flux
USE MOD_SurfInt,       ONLY: SurfInt
USE MOD_VolInt,        ONLY: VolInt
USE MOD_ProlongToFace, ONLY: ProlongToFace
USE MOD_FillFlux,      ONLY: FillFlux,FillFlux_BC
USE MOD_Mesh_Vars,     ONLY: sJ,Elem_xGP,nSides,nBCSides,nInnerSides
USE MOD_Equation,      ONLY: CalcSource
USE MOD_Equation_Vars, ONLY: IniExactFunc
#ifdef MPI
USE MOD_MPI_Vars
USE MOD_MPI,           ONLY:StartExchangeMPIData,FinishExchangeMPIData
USE MOD_Mesh_Vars,     ONLY:SideID_plus_upper,SideID_plus_lower
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: t,tStage
INTEGER,INTENT(IN)              :: tDeriv
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iElem,i,j,k,iVar
!===================================================================================================================================

! prolong the solution to the face integration points for flux computation
#ifdef MPI
! Prolong to face for MPI sides - send direction
CALL ProlongToFace(U,U_Minus,U_Plus,doMPiSides=.TRUE.)
CALL StartExchangeMPIData(U_Plus,SideID_plus_lower,SideID_plus_upper,SendRequest_U,RecRequest_U,SendID=2) ! Send YOUR - receive MINE
#endif /*MPI*/

! Prolong to face for BCSides, InnerSides and MPI sides - receive direction
CALL ProlongToFace(U,U_Minus,U_Plus,doMPISides=.FALSE.)

#ifdef MPI
! Complete send / receive
CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=2) !Send YOUR - receive MINE
#endif /*MPI*/

Ut=0.
CALL VolInt(Ut)

! Initialization of the time derivative
!Flux=0. !don't nullify the fluxes if not really needed (very expensive)
#ifdef MPI
! fill the global surface flux list
CALL FillFlux(Flux,doMPISides=.TRUE.)
CALL StartExchangeMPIData(Flux,1,nSides,SendRequest_Flux,RecRequest_Flux,SendID=1) ! Send MINE - receive YOUR
#endif /* MPI*/

! fill the all surface fluxes on this proc
CALL FillFlux_BC(t,tDeriv,Flux)
CALL FillFlux(Flux,doMPISides=.FALSE.)
! compute surface integral contribution and add to ut
CALL SurfInt(Flux,Ut,doMPISides=.FALSE.)
!! compute volume integral contribution and add to ut
!CALL VolInt(Ut)

#ifdef MPI
! Complete send / receive
CALL FinishExchangeMPIData(SendRequest_Flux,RecRequest_Flux,SendID=1) !Send MINE -receive YOUR
!FINALIZE Fluxes for MPI Sides
CALL SurfInt(Flux,Ut,doMPISides=.TRUE.)
#endif

! We have to take the inverse of the Jacobians into account
DO iElem=1,PP_nElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        DO iVar=1,PP_nVar
          Ut(iVar,i,j,k,iElem) = - Ut(iVar,i,j,k,iElem) * sJ(i,j,k,iElem)
        END DO ! iVar
      END DO !i
    END DO !j
  END DO !k
END DO ! iElem=1,nElems

! Add Source Terms
CALL CalcSource(tStage)

END SUBROUTINE DGTimeDerivative_weakForm


SUBROUTINE DGTimeDerivative_WoSource_weakForm(t,tStage,tDeriv)
!===================================================================================================================================
! Computes the DG time derivative consisting of Volume Integral and Surface integral for the whole field
! U and Ut are allocated
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_DG_Vars,       ONLY: U,Ut,nTotalU,U_Plus,U_Minus,Flux
USE MOD_SurfInt,       ONLY: SurfInt
USE MOD_VolInt,        ONLY: VolInt
USE MOD_ProlongToFace, ONLY: ProlongToFace
USE MOD_FillFlux,      ONLY: FillFlux,FillFlux_BC
USE MOD_Mesh_Vars,     ONLY: sJ,Elem_xGP,nSides,nBCSides,nInnerSides
!USE MOD_Equation,      ONLY: CalcSource
USE MOD_Equation_Vars, ONLY: IniExactFunc
#ifdef MPI
USE MOD_MPI_Vars
USE MOD_MPI,           ONLY:StartExchangeMPIData,FinishExchangeMPIData
USE MOD_Mesh_Vars,     ONLY:SideID_plus_upper,SideID_plus_lower
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: t,tStage
INTEGER,INTENT(IN)              :: tDeriv
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iElem,i,j,k,iVar
!===================================================================================================================================

! prolong the solution to the face integration points for flux computation
#ifdef MPI
! Prolong to face for MPI sides - send direction
CALL ProlongToFace(U,U_Minus,U_Plus,doMPiSides=.TRUE.)
CALL StartExchangeMPIData(U_Plus,SideID_plus_lower,SideID_plus_upper,SendRequest_U,RecRequest_U,SendID=2) ! Send YOUR - receive MINE
#endif /*MPI*/

! Prolong to face for BCSides, InnerSides and MPI sides - receive direction
CALL ProlongToFace(U,U_Minus,U_Plus,doMPISides=.FALSE.)

#ifdef MPI
! Complete send / receive
CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=2) !Send YOUR - receive MINE
#endif /*MPI*/

Ut=0.
CALL VolInt(Ut)

! Initialization of the time derivative
!Flux=0. !don't nullify the fluxes if not really needed (very expensive)
#ifdef MPI
! fill the global surface flux list
CALL FillFlux(Flux,doMPISides=.TRUE.)
CALL StartExchangeMPIData(Flux,1,nSides,SendRequest_Flux,RecRequest_Flux,SendID=1) ! Send MINE - receive YOUR
#endif /* MPI*/

! fill the all surface fluxes on this proc
CALL FillFlux_BC(t,tDeriv,Flux)
CALL FillFlux(Flux,doMPISides=.FALSE.)
! compute surface integral contribution and add to ut
CALL SurfInt(Flux,Ut,doMPISides=.FALSE.)
!! compute volume integral contribution and add to ut
!CALL VolInt(Ut)

#ifdef MPI
! Complete send / receive
CALL FinishExchangeMPIData(SendRequest_Flux,RecRequest_Flux,SendID=1) !Send MINE -receive YOUR
!FINALIZE Fluxes for MPI Sides
CALL SurfInt(Flux,Ut,doMPISides=.TRUE.)
#endif

! We have to take the inverse of the Jacobians into account
DO iElem=1,PP_nElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        DO iVar=1,PP_nVar
          Ut(iVar,i,j,k,iElem) = - Ut(iVar,i,j,k,iElem) * sJ(i,j,k,iElem)
        END DO ! iVar
      END DO !i
    END DO !j
  END DO !k
END DO ! iElem=1,nElems


END SUBROUTINE DGTimeDerivative_WoSource_weakForm



SUBROUTINE FillIni()
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars,ONLY:U
USE MOD_Mesh_Vars,ONLY:Elem_xGP
USE MOD_Equation_Vars,ONLY:IniExactFunc
USE MOD_Equation,ONLY:ExactFunc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iElem
!===================================================================================================================================
! Determine Size of the Loops, i.e. the number of grid cells in the
! corresponding directions
DO iElem=1,PP_nElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        CALL ExactFunc(IniExactFunc,0.,0,Elem_xGP(1:3,i,j,k,iElem),U(1:PP_nVar,i,j,k,iElem))
      END DO ! i
    END DO ! j
  END DO !k
END DO ! iElem=1,PP_nElems
END SUBROUTINE FillIni



SUBROUTINE FinalizeDG()
!===================================================================================================================================
! Deallocate global variable U (solution) and Ut (dg time derivative).
!===================================================================================================================================
! MODULES
USE MOD_DG_Vars,ONLY:D_Hat,L_HatMinus,L_HatPlus,U,Ut, U_Minus,U_Plus
USE MOD_DG_Vars,ONLY:Flux,DGInitIsDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================

SDEALLOCATE(D_Hat)
SDEALLOCATE(L_HatMinus)
SDEALLOCATE(L_HatPlus)
SDEALLOCATE(Ut)
SDEALLOCATE(U)
SDEALLOCATE(U_Minus)
SDEALLOCATE(U_Plus)
DGInitIsDone = .FALSE.
END SUBROUTINE FinalizeDG

END MODULE MOD_DG
