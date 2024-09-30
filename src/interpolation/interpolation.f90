!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (piclas.boltzplatz.eu/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

MODULE MOD_Interpolation
!===================================================================================================================================
!> Contains routines to prepare for interpolation procedures:
!> - Initialize interpolation variables
!> - Calculate node positions and weights
!> - Build Vandermonde matrices
!> - Build derivative matrices
!> Also contains routines to map the solution between physical and reference space.
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

INTERFACE GetNodesAndWeights
   MODULE PROCEDURE GetNodesAndWeights
END INTERFACE

INTERFACE GetDerivativeMatrix
   MODULE PROCEDURE GetDerivativeMatrix
END INTERFACE

INTERFACE ApplyJacobian
   MODULE PROCEDURE ApplyJacobian
END INTERFACE

INTERFACE InitInterpolationBasis
   MODULE PROCEDURE InitInterpolationBasis
END INTERFACE

PUBLIC::InitInterpolation
PUBLIC::ApplyJacobian
PUBLIC::FinalizeInterpolation
PUBLIC::GetNodesAndWeights
PUBLIC::GetDerivativeMatrix
PUBLIC::GetVandermonde
PUBLIC::InitInterpolationBasis

PUBLIC::DefineParametersInterpolation

!===================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersInterpolation()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("Interpolation")
CALL prms%CreateIntOption('N'        , "Polynomial degree of computation to represent to solution")
CALL prms%CreateIntOption('NAnalyze' , 'Polynomial degree at which analysis is performed (e.g. for L2 errors).\n'//&
                                       'Default: 2*(N+1).')
END SUBROUTINE DefineParametersInterpolation


SUBROUTINE InitInterpolation(NIn,NAnalyzeIn)
!============================================================================================================================
! Initialize basis for Gauss-points of order N.
! Prepares Differentiation matrices D, D_Hat, Basis at the boundaries L(1), L(-1), L_Hat(1), L_Hat(-1)
! Gaussnodes xGP and weights wGP, wBary
!============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Interpolation_Vars
USE MOD_ReadInTools        ,ONLY: GETINT,CountOption
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
!input parameters
INTEGER,INTENT(IN),OPTIONAL :: NIn         !< optional polynomial degree
INTEGER,INTENT(IN),OPTIONAL :: NAnalyzeIn  !< optional analyze polynomial degree
!----------------------------------------------------------------------------------------------------------------------------
!output parameters
!----------------------------------------------------------------------------------------------------------------------------
!local variables
CHARACTER(LEN=40)           :: DefStr
#if !(PP_N == N)
INTEGER                     :: Ntmp
#endif /*!(PP_N == N)*/
!============================================================================================================================
IF (InterpolationInitIsDone) CALL CollectiveStop(__STAMP__,'InitInterpolation already called.')
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT INTERPOLATION...'

! Access ini-file
#if PP_N == N
IF(PRESENT(Nin))THEN
  PP_N = NIn
ELSE
  PP_N = GETINT('N')
END IF
#else
IF(PRESENT(Nin))THEN
  Ntmp = NIn
ELSE
  Ntmp=PP_N
  IF(CountOption('N').EQ.1) Ntmp=GETINT('N')
END IF
IF(PP_N.NE.Ntmp) CALL CollectiveStop(__STAMP__,'N in ini-file is different from hard-compiled N. Ini/Compiled:',Ntmp,REAL(PP_N))
#endif

SWRITE(UNIT_stdOut,'(A)') ' NodeType: '//NodeType
!CALL InitInterpolationBasis(PP_N, xGP ,wGP, swGP,wBary ,L_Minus ,L_Plus , L_PlusMinus, wGPSurf, Vdm_Leg ,sVdm_Leg)
CALL InitInterpolationBasis(PP_N, xGP, wGP, wBary, L_Minus , L_Plus, L_PlusMinus, swGP=swGP, wGPSurf=wGPSurf)

! Set the default analyze polynomial degree NAnalyze to 2*(N+1)
IF(PRESENT(NAnalyzeIn))THEN
  NAnalyze = NAnalyzeIn
ELSE
  WRITE(DefStr,'(i4)') 2*(PP_N+1)
  NAnalyze = GETINT('NAnalyze',DefStr)
END IF ! PRESENT(NAnalyzeIn)

! Initialize the basis functions for the analyze polynomial
CALL InitAnalyzeBasis(PP_N,NAnalyze,xGP,wBary)

InterpolationInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT INTERPOLATION DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitInterpolation


SUBROUTINE InitInterpolationBasis(N_in, xGP ,wGP, wBary ,L_Minus ,L_Plus , L_PlusMinus,swGP, wGPSurf, Vdm_Leg ,sVdm_Leg, NodeType_in)
!============================================================================================================================
! Initialize basis for Gauss-points of order N.
! Calculate positions of Gauss-points, integration weights and barycentric weights. Prepare basis evaluation at -1 and +1.
! Prepares Differentiation matrices D, D_Hat, Basis at the boundaries L(1), L(-1), L_Hat(1), L_Hat(-1)
! Gaussnodes xGP and weights wGP, wBary
!============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Interpolation_Vars,  ONLY:NodeType
USE MOD_Basis,               ONLY:LagrangeInterpolationPolys,buildLegendreVdm
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                  :: N_in
REAL,ALLOCATABLE,DIMENSION(:),  INTENT(OUT)         :: xGP                !< Gauss points in [-1,1]
REAL,ALLOCATABLE,DIMENSION(:),  INTENT(OUT)         :: wGP                !< Integration weights
REAL,ALLOCATABLE,DIMENSION(:),  INTENT(OUT)         :: wBary              !< Barycentric weights
REAL,ALLOCATABLE,DIMENSION(:),  INTENT(OUT)         :: L_Minus            !< Lagrange polynomials at -1
REAL,ALLOCATABLE,DIMENSION(:),  INTENT(OUT)         :: L_Plus             !< Lagrange polynomials at +1
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT)         :: L_PlusMinus        !< Vandermonde Nodal->Modal
REAL,ALLOCATABLE,DIMENSION(:),  INTENT(OUT),OPTIONAL:: swGP                !< Integration weights
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT),OPTIONAL:: wGPSurf            !< Vandermonde Nodal->Modal
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT),OPTIONAL:: Vdm_Leg            !< Vandermonde Nodal->Modal
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT),OPTIONAL:: sVdm_Leg           !< Vandermonde Modal->Nodal
CHARACTER(LEN=*),INTENT(IN),OPTIONAL                :: NodeType_in        !< Type of 1D points
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                     :: i,j
!============================================================================================================================
! Allocate global variables, needs to go somewhere else later
ALLOCATE(xGP(0:N_in), wGP(0:N_in), wBary(0:N_in))
ALLOCATE(L_Minus(0:N_in), L_Plus(0:N_in))
ALLOCATE(L_PlusMinus(0:N_in,6))

IF(PRESENT(NodeType_in))THEN
  CALL GetNodesAndWeights(N_in , NodeType_in , xGP , wGP , wBary)
ELSE
  CALL GetNodesAndWeights(N_in , NodeType    , xGP , wGP , wBary)
END IF ! PRESENT(NodeType_in)

IF(PRESENT(wGPSurf))THEN
  ALLOCATE(wGPSurf(0:N_in,0:N_in))
  DO i=0,N_in;DO j=0,N_in;
    wGPSurf(i,j)  = wGP(i)*wGP(j)
  END DO; END DO;
END IF

IF(PRESENT(swGP))THEN
  ALLOCATE(swGP(0:N_in))
  DO i=0,N_in
    swGP(i)=1.0/wGP(i)
  END DO ! i
END IF

!! interpolate to left and right face (1 and -1) and pre-divide by mass matrix
CALL LagrangeInterpolationPolys(1.,N_in,xGP,wBary,L_Plus)
CALL LagrangeInterpolationPolys(-1.,N_in,xGP,wBary,L_Minus)
L_PlusMinus(:,  XI_MINUS) = L_Minus
L_PlusMinus(:, ETA_MINUS) = L_Minus
L_PlusMinus(:,ZETA_MINUS) = L_Minus
L_PlusMinus(:,  XI_PLUS)  = L_Plus
L_PlusMinus(:, ETA_PLUS)  = L_Plus
L_PlusMinus(:,ZETA_PLUS)  = L_Plus

IF(PRESENT(Vdm_Leg).AND.PRESENT(sVdm_Leg))THEN
! Vandermonde from NODAL <--> MODAL
ALLOCATE(Vdm_Leg(0:N_in,0:N_in),sVdm_Leg(0:N_in,0:N_in))
CALL buildLegendreVdm(N_in,xGP,Vdm_Leg,sVdm_Leg)
END IF

END SUBROUTINE InitInterpolationBasis


SUBROUTINE GetNodesAndWeights(N_in,NodeType_in,xIP,wIP,wIPBary)
!==================================================================================================================================
! Compute 1D nodes and weights for several node types in interval [-1,1]
!==================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Basis,       ONLY: LegendreGaussNodesAndWeights,LegGaussLobNodesAndWeights,ChebyGaussLobNodesAndWeights,BarycentricWeights
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                 :: N_in            !< Polynomial degree
CHARACTER(LEN=*),INTENT(IN)        :: NodeType_in     !< Type of 1D points
REAL,INTENT(OUT)                   :: xIP(0:N_in)     !< Position of nodes
REAL,INTENT(OUT),OPTIONAL          :: wIP(0:N_in)     !< Integration weights
REAL,INTENT(OUT),OPTIONAL          :: wIPBary(0:N_in) !< Barycentric weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: i
!==================================================================================================================================
IF(PRESENT(wIP))THEN
  SELECT CASE(TRIM(NodeType_in))
  CASE('GAUSS')
    CALL LegendreGaussNodesAndWeights(N_in,xIP,wIP)
  CASE('GAUSS-LOBATTO')
    CALL LegGaussLobNodesAndWeights(N_in,xIP,wIP)
  CASE('CHEBYSHEV-GAUSS-LOBATTO')
    CALL ChebyGaussLobNodesAndWeights(N_in,xIP,wIP)
  CASE('VISU')
    DO i=0,N_in
      xIP(i) = 2.*REAL(i)/REAL(N_in) - 1.
    END DO
    ! Trapez rule for integration !!!
    wIP(:) = 2./REAL(N_in)
    wIP(0) = 0.5*wIP(0)
    wIP(N_in) = 0.5*wIP(N_in)
  CASE('VISU_INNER')
    DO i=0,N_in
      xIP(i) = 1./REAL(N_in+1)+2.*REAL(i)/REAL(N_in+1) - 1.
    END DO
    ! first order intergration !!!
    wIP=2./REAL(N_in+1)
  CASE DEFAULT
    CALL Abort(__STAMP__,&
      'NodeType "'//TRIM(NodeType_in)//'" in GetNodesAndWeights not found!')
  END SELECT
ELSE
  SELECT CASE(TRIM(NodeType_in))
  CASE('GAUSS')
    CALL LegendreGaussNodesAndWeights(N_in,xIP)
  CASE('GAUSS-LOBATTO')
    CALL LegGaussLobNodesAndWeights(N_in,xIP)
  CASE('CHEBYSHEV-GAUSS-LOBATTO')
    CALL ChebyGaussLobNodesAndWeights(N_in,xIP)
  CASE('VISU')
    DO i=0,N_in
      xIP(i) = 2.*REAL(i)/REAL(N_in) - 1.
    END DO
  CASE('VISU_INNER')
    DO i=0,N_in
      xIP(i) = 1./REAL(N_in+1)+2.*REAL(i)/REAL(N_in+1) - 1.
    END DO
  CASE DEFAULT
    CALL Abort(__STAMP__,&
      'NodeType "'//TRIM(NodeType_in)//'" in GetNodesAndWeights not found!')
  END SELECT
END IF !present wIP
IF(PRESENT(wIPBary)) CALL BarycentricWeights(N_in,xIP,wIPBary)
END SUBROUTINE GetNodesAndWeights

SUBROUTINE GetVandermonde(N_in,NodeType_in,N_out,NodeType_out,Vdm_In_Out,Vdm_Out_In,modal)
!==================================================================================================================================
!> Build a Vandermonde-Matrix from/to different node types and polynomial degrees.
!==================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Basis   ,ONLY: BarycentricWeights,InitializeVandermonde
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                 :: N_in                       !> Input polynomial degree
INTEGER,INTENT(IN)                 :: N_out                      !> Output polynomial degree
CHARACTER(LEN=*),INTENT(IN)        :: NodeType_in                !> Type of 1D input points
CHARACTER(LEN=*),INTENT(IN)        :: NodeType_out               !> Type of 1D output points
LOGICAL,INTENT(IN),OPTIONAL        :: modal                      !> Switch if a modal Vandermonde should be build
REAL,INTENT(OUT)                   :: Vdm_In_out(0:N_out,0:N_in) !> Vandermonde In->Out
REAL,INTENT(OUT),OPTIONAL          :: Vdm_Out_In(0:N_in,0:N_out) !> Vandermonde Out->in
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: i
REAL                               :: xIP_in(0:N_in)
REAL                               :: xIP_out(0:N_out)
REAL                               :: wBary_in(0:N_in)
REAL                               :: wBary_out(0:N_out)
REAL                               ::  Vdm_Leg_in( 0:N_in,0:N_in)
REAL                               :: sVdm_Leg_in( 0:N_in,0:N_in)
REAL                               ::  Vdm_Leg_out(0:N_out,0:N_out)
REAL                               :: sVdm_Leg_out(0:N_out,0:N_out)
LOGICAL                            :: modalLoc
!==================================================================================================================================
modalLoc=.FALSE.
IF(PRESENT(modal)) modalLoc=modal

! Check if change Basis is needed
IF((TRIM(NodeType_out).EQ.TRIM(NodeType_in)).AND.(N_in.EQ.N_out))THEN
  Vdm_In_Out=0.
  DO i=0,N_in
    Vdm_In_out(i,i)=1.
  END DO
  IF(PRESENT(Vdm_Out_In))THEN
    Vdm_Out_In=0.
    DO i=0,N_Out
      Vdm_Out_In(i,i)=1.
    END DO
  END IF
ELSE
  ! Input points
  CALL GetNodesAndWeights(N_in,NodeType_in,xIP_in)
  CALL BarycentricWeights(N_in,xIP_in,wBary_in)
  ! Output points
  CALL GetNodesAndWeights(N_out,NodeType_out,xIP_out)

  IF(modalLoc)THEN
    CALL buildLegendreVdm(N_In, xIP_in, Vdm_Leg_in, sVdm_Leg_in)
    CALL buildLegendreVdm(N_Out,xIP_out,Vdm_Leg_out,sVdm_Leg_out)
  END IF

  IF((N_Out.LT.N_In).AND.modalLoc)THEN
    Vdm_In_Out=MATMUL(Vdm_Leg_Out(0:N_Out,0:N_Out),sVdm_Leg_In(0:N_Out,0:N_In))
  ELSE
    CALL InitializeVandermonde(N_in,N_out,wBary_in,xIP_in,xIP_out,Vdm_In_Out)
  END IF
  IF(PRESENT(Vdm_Out_In))THEN
    IF((N_In.LT.N_Out).AND.modalLoc)THEN
      Vdm_Out_In=MATMUL(Vdm_Leg_In(0:N_In,0:N_In),sVdm_Leg_Out(0:N_In,0:N_Out))
    ELSE
      CALL BarycentricWeights(N_out,xIP_out,wBary_out)
      CALL InitializeVandermonde(N_out,N_in,wBary_out,xIP_out,xIP_in,Vdm_Out_In)
    END IF
  END IF
END IF
END SUBROUTINE GetVandermonde


SUBROUTINE GetDerivativeMatrix(N_in,NodeType_in,D)
!==================================================================================================================================
!> Compute polynomial derivative matrix. D(i,j) = Derivative of basis function j evaluated at node point i.
!==================================================================================================================================
! MODULES
USE MOD_Basis,             ONLY:PolynomialDerivativeMatrix
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                 :: N_in                       !< Polynomial degree
CHARACTER(LEN=255),INTENT(IN)      :: NodeType_in                !< Type of 1D input points
REAL,INTENT(OUT)                   :: D(0:N_in,0:N_in)           !< Derivative matrix
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                               :: xIP(0:N_in)
!==================================================================================================================================
CALL GetNodesAndWeights(N_in,NodeType_in,xIP)
CALL PolynomialDerivativeMatrix(N_in,xIP,D)
END SUBROUTINE GetDerivativeMatrix


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


SUBROUTINE InitAnalyzeBasis(N_in,Nanalyze_in,xGP,wBary)
!===================================================================================================================================
! Build analyze nodes (Gauss-Lobatto) and corresponding Vandermonde matrix
!===================================================================================================================================
! MODULES
USE MOD_Interpolation_Vars ,ONLY: wAnalyze ! GL integration weights used for the analyze
USE MOD_Interpolation_Vars ,ONLY: Vdm_GaussN_NAnalyze
USE MOD_Basis              ,ONLY: LegGaussLobNodesAndWeights,InitializeVandermonde
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)   :: N_in          !< input polynomial degree
INTEGER,INTENT(IN)   :: Nanalyze_in   !< polynomial degree of analysis polynomial
REAL,INTENT(IN)      :: xGP(0:N_in)   !< interpolation points
REAL,INTENT(IN)      :: wBary(0:N_in) !< barycentric weights
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL ,DIMENSION(0:Nanalyze_in) :: XiAnalyze ! GL nodes
!===================================================================================================================================
ALLOCATE(wAnalyze(0:NAnalyze_in),Vdm_GaussN_NAnalyze(0:NAnalyze_in,0:N_in))
! Build analyze nodes (Gauss-Lobatto)
CALL LegGaussLobNodesAndWeights(NAnalyze_in,XiAnalyze,wAnalyze)
! Build analyze Vandermonde matrix which maps from NodeType nodes to Gauss-Lobatto nodes
CALL InitializeVandermonde(N_in,NAnalyze_in,wBary,xGP,XiAnalyze,Vdm_GaussN_NAnalyze)
END SUBROUTINE InitAnalyzeBasis


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
SDEALLOCATE(swGP)
SDEALLOCATE(wGPSurf)
SDEALLOCATE(wBary)
SDEALLOCATE(NChooseK)
SDEALLOCATE(L_Minus)
SDEALLOCATE(L_Plus)
SDEALLOCATE(L_PlusMinus)
SDEALLOCATE(Vdm_GaussN_NAnalyze)
SDEALLOCATE(wAnalyze)

InterpolationInitIsDone = .FALSE.
END SUBROUTINE FinalizeInterpolation

END MODULE MOD_Interpolation

