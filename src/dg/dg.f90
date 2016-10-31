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

#ifndef PP_HDG
INTERFACE DGTimeDerivative_weakForm
  MODULE PROCEDURE DGTimeDerivative_weakForm
END INTERFACE
#endif /*PP_HDG*/

INTERFACE FinalizeDG
  MODULE PROCEDURE FinalizeDG
END INTERFACE

PUBLIC::InitDG,FinalizeDG
#ifndef PP_HDG
PUBLIC::DGTimeDerivative_weakForm
#endif /*PP_HDG*/
#ifdef PP_POIS
PUBLIC::DGTimeDerivative_weakForm_Pois
#endif
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
USE MOD_Mesh_Vars,          ONLY: nSides
USE MOD_Mesh_Vars,          ONLY: MeshInitIsDone
#ifdef OPTIMIZED
USE MOD_Riemann,            ONLY: GetRiemannMatrix
#endif /*OPTIMIZED*/
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
   CALL abort(&
       __STAMP__&
       ,'InitDG not ready to be called or already called.',999,999.)
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT DG...'

CALL initDGbasis(PP_N,xGP,wGP,wBary)
! the local DG solution in physical and reference space
ALLOCATE( U(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
! the time derivative computed with the DG scheme
ALLOCATE(Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
U=0.
Ut=0.

nTotal_face=(PP_N+1)*(PP_N+1)
nTotal_vol=nTotal_face*(PP_N+1)
nTotalU=PP_nVar*nTotal_vol*PP_nElems

! U is filled with the ini solution
IF(.NOT.DoRestart) CALL FillIni()

! We store the interior data at the each element face
!ALLOCATE(U_Minus(PP_nVar,0:PP_N,0:PP_N,sideID_minus_lower:sideID_minus_upper))
!ALLOCATE(U_Plus(PP_nVar,0:PP_N,0:PP_N,sideID_plus_lower:sideID_plus_upper))
!U_Minus=0.
!U_Plus=0.

ALLOCATE(U_master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(U_slave(PP_nVar,0:PP_N,0:PP_N,1:nSides))
U_master=0.
U_slave=0.


#ifdef OPTIMIZED
  CALL GetRiemannMatrix()
#endif /*OPTIMIZED*/

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
USE MOD_Globals
USE MOD_Basis     ,ONLY:LegendreGaussNodesAndWeights,LegGaussLobNodesAndWeights,BarycentricWeights
USE MOD_Basis     ,ONLY:PolynomialDerivativeMatrix,LagrangeInterpolationPolys
USE MOD_DG_Vars   ,ONLY:D,D_T,D_Hat,D_Hat_T,L_HatMinus,L_HatPlus
#ifdef PP_HDG
#ifdef MPI
USE MOD_PreProc
USE MOD_MPI_vars,      ONLY:SendRequest_Geo,RecRequest_Geo
USE MOD_MPI,           ONLY:StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
USE MOD_Mesh_Vars,     ONLY:NormVec,TangVec1,TangVec2,SurfElem!,Face_xGP
#endif /*MPI*/
#endif /*PP_HDG*/
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
#ifdef PP_HDG
#ifdef MPI  
REAL                                       :: Geotemp(10,0:PP_N,0:PP_N,1:nSides)
#endif /*MPI*/
#endif /*PP_HDG*/
!===================================================================================================================================
ALLOCATE(L_HatMinus(0:N_in), L_HatPlus(0:N_in))
ALLOCATE(D(0:N_in,0:N_in), D_T(0:N_in,0:N_in))
ALLOCATE(D_Hat(0:N_in,0:N_in), D_Hat_T(0:N_in,0:N_in))
! Compute Differentiation matrix D for given Gausspoints
CALL PolynomialDerivativeMatrix(N_in,xGP,D)
D_T=TRANSPOSE(D)

! Build D_Hat matrix. (D^ = M^(-1) * D^T * M
M(:,:)=0.
Minv(:,:)=0.
DO iMass=0,N_in
  M(iMass,iMass)=wGP(iMass)
  Minv(iMass,iMass)=1./wGP(iMass)
END DO
D_Hat(:,:) = -MATMUL(Minv,MATMUL(TRANSPOSE(D),M))
D_Hat_T=TRANSPOSE(D_hat)

! interpolate to left and right face (1 and -1) and pre-divide by mass matrix
CALL LagrangeInterpolationPolys(1.,N_in,xGP,wBary,L_Plus)
L_HatPlus(:) = MATMUL(Minv,L_Plus)
CALL LagrangeInterpolationPolys(-1.,N_in,xGP,wBary,L_Minus)
L_HatMinus(:) = MATMUL(Minv,L_Minus)

#ifdef PP_HDG
#ifdef MPI
! exchange is in initDGbasis as InitMesh() and InitMPI() is needed
Geotemp=0.
Geotemp(1,:,:,:)=SurfElem(:,:,1:nSides)
Geotemp(2:4,:,:,:)=NormVec(:,:,:,1:nSides)
Geotemp(5:7,:,:,:)=TangVec1(:,:,:,1:nSides)
Geotemp(8:10,:,:,:)=TangVec2(:,:,:,1:nSides)
!Geotemp(11:13,:,:,:)=Face_xGP(:,:,:,SideID_minus_lower:SideID_minus_upper)

IPWRITE(*,*) size(Geotemp)
CALL StartReceiveMPIData(10,Geotemp,1,nSides,RecRequest_Geo ,SendID=1) ! Receive MINE
CALL StartSendMPIData(   10,Geotemp,1,nSides,SendRequest_Geo,SendID=1) ! Send YOUR
CALL FinishExchangeMPIData(SendRequest_Geo,RecRequest_Geo,SendID=1)                                 ! Send YOUR - receive MINE

SurfElem(:,:,1:nSides)=Geotemp(1,:,:,:)
NormVec(:,:,:,1:nSides)=Geotemp(2:4,:,:,:)
TangVec1(:,:,:,1:nSides)=Geotemp(5:7,:,:,:)
TangVec2(:,:,:,1:nSides)=Geotemp(8:10,:,:,:)
!Face_xGP(:,:,:,SideID_minus_lower:SideID_minus_upper)=Geotemp(11:13,:,:,:)

#endif /*MPI*/
#endif /*PP_HDG*/
END SUBROUTINE InitDGbasis


#ifndef PP_HDG
SUBROUTINE DGTimeDerivative_weakForm(t,tStage,tDeriv,doSource)
!===================================================================================================================================
! Computes the DG time derivative consisting of Volume Integral and Surface integral for the whole field
! U and Ut are allocated
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Vector
USE MOD_DG_Vars,          ONLY:U,Ut,U_master,U_slave,Flux !,nTotalU
USE MOD_SurfInt,          ONLY:SurfInt
USE MOD_VolInt,           ONLY:VolInt
USE MOD_ProlongToFace,    ONLY:ProlongToFace
USE MOD_FillFlux,         ONLY:FillFlux
USE MOD_Mesh_Vars,        ONLY:nSides
USE MOD_Equation,         ONLY:CalcSource
USE MOD_Interpolation,    ONLY:ApplyJacobian
USE MOD_FillMortar,       ONLY:U_Mortar,Flux_Mortar
#ifdef MPI
USE MOD_MPI_Vars
USE MOD_MPI,              ONLY:StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
USE MOD_LoadBalance_Vars, ONLY:tCurrent
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: t,tStage
INTEGER,INTENT(IN)              :: tDeriv
LOGICAL,INTENT(IN)              :: doSource
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#ifdef MPI
! load balance
REAL                            :: tLBStart,tLBEnd
#endif /*MPI*/
!===================================================================================================================================

! prolong the solution to the face integration points for flux computation
#ifdef MPI
! Prolong to face for MPI sides - send direction
tLBStart = LOCALTIME() ! LB Time Start
CALL StartReceiveMPIData(PP_nVar,U_slave,1,nSides,RecRequest_U,SendID=2) ! Receive MINE
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(2)=tCurrent(2)+tLBEnd-tLBStart
tLBStart = LOCALTIME() ! LB Time Start
CALL ProlongToFace(U,U_master,U_slave,doMPISides=.TRUE.)
CALL U_Mortar(U_master,U_slave,doMPISides=.TRUE.)
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(1)=tCurrent(1)+tLBEnd-tLBStart
tLBStart = LOCALTIME() ! LB Time Start
CALL StartSendMPIData(PP_nVar,U_slave,1,nSides,SendRequest_U,SendID=2) ! Send YOUR
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(2)=tCurrent(2)+tLBEnd-tLBStart
tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/

! Prolong to face for BCSides, InnerSides and MPI sides - receive direction
CALL ProlongToFace(U,U_master,U_slave,doMPISides=.FALSE.)
CALL U_Mortar(U_master,U_slave,doMPISides=.FALSE.)
! Nullify arrays
! NOTE: IF NEW DG_VOLINT AND LIFTING_VOLINT ARE USED AND CALLED FIRST,
!       ARRAYS DO NOT NEED TO BE NULLIFIED, OTHERWISE THEY HAVE TO!
!CALL VNullify(nTotalU,Ut)
Ut=0.
! compute volume integral contribution and add to ut, first half of all elements
CALL VolInt(Ut,dofirstElems=.TRUE.)

#ifdef MPI
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(1)=tCurrent(1)+tLBEnd-tLBStart
! Complete send / receive
tLBStart = LOCALTIME() ! LB Time Start
CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=2) !Send YOUR - receive MINE
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(2)=tCurrent(2)+tLBEnd-tLBStart
#endif /*MPI*/

! Initialization of the time derivative
!Flux=0. !don't nullify the fluxes if not really needed (very expensive)
#ifdef MPI
tLBStart = LOCALTIME() ! LB Time Start
CALL StartReceiveMPIData(PP_nVar,Flux,1,nSides,RecRequest_Flux,SendID=1) ! Receive MINE
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(2)=tCurrent(2)+tLBEnd-tLBStart
! fill the global surface flux list
tLBStart = LOCALTIME() ! LB Time Start
CALL FillFlux(t,tDeriv,Flux,U_master,U_slave,doMPISides=.TRUE.)
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(1)=tCurrent(1)+tLBEnd-tLBStart

tLBStart = LOCALTIME() ! LB Time Start
CALL StartSendMPIData(PP_nVar,Flux,1,nSides,SendRequest_Flux,SendID=1) ! Send YOUR
!CALL StartExchangeMPIData(PP_nVar,Flux,1,nSides,SendRequest_Flux,RecRequest_Flux,SendID=1) ! Send MINE - receive YOUR
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(2)=tCurrent(2)+tLBEnd-tLBStart
tLBStart = LOCALTIME() ! LB Time Start
#endif /* MPI*/

! fill the all surface fluxes on this proc
CALL FillFlux(t,tDeriv,Flux,U_master,U_slave,doMPISides=.FALSE.)
CALL Flux_Mortar(Flux,doMPISides=.FALSE.)
! compute surface integral contribution and add to ut
CALL SurfInt(Flux,Ut,doMPISides=.FALSE.)

! compute volume integral contribution and add to ut
CALL VolInt(Ut,dofirstElems=.FALSE.)

#ifdef MPI
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(1)=tCurrent(1)+tLBEnd-tLBStart
! Complete send / receive
tLBStart = LOCALTIME() ! LB Time Start
CALL FinishExchangeMPIData(SendRequest_Flux,RecRequest_Flux,SendID=1) !Send MINE -receive YOUR
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(2)=tCurrent(2)+tLBEnd-tLBStart

!FINALIZE Fluxes for MPI Sides
tLBStart = LOCALTIME() ! LB Time Start
CALL Flux_Mortar(Flux,doMPISides=.TRUE.)
CALL SurfInt(Flux,Ut,doMPISides=.TRUE.)
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(1)=tCurrent(1)+tLBEnd-tLBStart
tLBStart = LOCALTIME() ! LB Time Start
#endif

! swap and map to physical space
CALL ApplyJacobian(Ut,toPhysical=.TRUE.,toSwap=.TRUE.)

! Add Source Terms
IF(doSource) CALL CalcSource(tStage,1.0,Ut)

#ifdef MPI
tLBEnd = LOCALTIME() ! LB Time End
tCurrent(1)=tCurrent(1)+tLBEnd-tLBStart
#endif /*MPI*/

END SUBROUTINE DGTimeDerivative_weakForm
#endif /*PP_HDG*/

#ifdef PP_POIS

SUBROUTINE DGTimeDerivative_weakForm_Pois(t,tStage,tDeriv)
!===================================================================================================================================
! Computes the DG time derivative consisting of Volume Integral and Surface integral for the whole field
! U and Ut are allocated
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Vector
USE MOD_Equation,      ONLY: VolInt_Pois,FillFlux_Pois,ProlongToFace_Pois, SurfInt_Pois
USE MOD_GetBoundaryFlux, ONLY: FillFlux_BC_Pois
USE MOD_Mesh_Vars,     ONLY: sJ,Elem_xGP,nSides
USE MOD_Equation,      ONLY: CalcSource_Pois
USE MOD_Equation_Vars, ONLY: IniExactFunc,Phi,Phit,Phi_master,Phi_slave,FluxPhi,nTotalPhi
USE MOD_Interpolation, ONLY: ApplyJacobian
#ifdef MPI
USE MOD_MPI_Vars
USE MOD_MPI,           ONLY:StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
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

!CALL ProlongToFace(Phi,Phi_Minus,Phi_slave,doMPiSides=.TRUE.)
CALL StartReceiveMPIData(4,Phi_slave,1,nSides,RecRequest_U,SendID=2) ! Receive MINE
CALL ProlongToFace_Pois(Phi,Phi_master,Phi_slave,doMPiSides=.TRUE.)

!CALL StartExchangeMPIData(Phi_slave,SideID_plus_lower,SideID_plus_upper,SendRequest_U,RecRequest_U,SendID=2) 
CALL StartSendMPIData(4,Phi_slave,1,nSides,SendRequest_U,SendID=2) ! Send YOUR
! Send YOUR - receive MINE
#endif /*MPI*/

! Prolong to face for BCSides, InnerSides and MPI sides - receive direction
!CALL ProlongToFace(Phi,Phi_Minus,Phi_slave,doMPISides=.FALSE.)
CALL ProlongToFace_Pois(Phi,Phi_master,Phi_slave,doMPISides=.FALSE.)

Phit=0.
CALL VolInt_Pois(Phit)
!print*,'Phi',Phit(:,1,1,1,4)
!read*


#ifdef MPI
! Complete send / receive
CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=2) !Send YOUR - receive MINE
#endif /*MPI*/


! Initialization of the time derivative
!Flux=0. !don't nullify the fluxes if not really needed (very expensive)
#ifdef MPI
! fill the global surface flux list
CALL StartReceiveMPIData(4,FluxPhi,1,nSides,RecRequest_Flux,SendID=1) ! Receive MINE
CALL FillFlux_Pois(FluxPhi,doMPISides=.TRUE.)

!CALL StartExchangeMPIData(FluxPhi,1,nSides,SendRequest_Flux,RecRequest_Flux,SendID=1) ! Send MINE - receive YOUR
CALL StartSendMPIData(4,FluxPhi,1,nSides,SendRequest_Flux,SendID=1) ! Send YOUR
!CALL StartExchangeMPIData(4,FluxPhi,1,nSides,SendRequest_Flux,RecRequest_Flux,SendID=1)
#endif /* MPI*/

! fill the all surface fluxes on this proc
CALL FillFlux_BC_Pois(t,tDeriv,FluxPhi)
CALL FillFlux_Pois(FluxPhi,doMPISides=.FALSE.)
! compute surface integral contribution and add to ut
CALL SurfInt_Pois(FluxPhi,Phit,doMPISides=.FALSE.)
!! compute volume integral contribution and add to ut
!CALL VolInt(Ut)

#ifdef MPI
! Complete send / receive
CALL FinishExchangeMPIData(SendRequest_Flux,RecRequest_Flux,SendID=1) !Send MINE -receive YOUR
!FINALIZE Fluxes for MPI Sides
CALL SurfInt_Pois(FluxPhi,Phit,doMPISides=.TRUE.)
#endif

! We have to take the inverse of the Jacobians into account
DO iElem=1,PP_nElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        DO iVar=1,4
          Phit(iVar,i,j,k,iElem) = - Phit(iVar,i,j,k,iElem) * sJ(i,j,k,iElem)
        END DO ! iVar
      END DO !i
    END DO !j
  END DO !k
END DO ! iElem=1,nElems

! Add Source Terms
CALL CalcSource_Pois(tStage)

END SUBROUTINE DGTimeDerivative_weakForm_Pois

#endif

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
#ifdef PP_HDG
        CALL ExactFunc(IniExactFunc,Elem_xGP(1:3,i,j,k,iElem),U(1:PP_nVar,i,j,k,iElem))
#else
        CALL ExactFunc(IniExactFunc,0.,0,Elem_xGP(1:3,i,j,k,iElem),U(1:PP_nVar,i,j,k,iElem))
#endif
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
USE MOD_DG_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
SDEALLOCATE(D)
SDEALLOCATE(D_T)
SDEALLOCATE(D_Hat)
SDEALLOCATE(D_Hat_T)
SDEALLOCATE(L_HatMinus)
SDEALLOCATE(L_HatPlus)
SDEALLOCATE(Ut)
SDEALLOCATE(U)
SDEALLOCATE(U_master)
SDEALLOCATE(U_slave)
SDEALLOCATE(FLUX)
DGInitIsDone = .FALSE.
END SUBROUTINE FinalizeDG

END MODULE MOD_DG
