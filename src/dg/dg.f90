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

#if !(USE_HDG)
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
INTERFACE DGTimeDerivative_weakForm
  MODULE PROCEDURE DGTimeDerivative_weakForm
END INTERFACE
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/
#endif /*USE_HDG*/

INTERFACE FinalizeDG
  MODULE PROCEDURE FinalizeDG
END INTERFACE

PUBLIC::InitDG,FinalizeDG
#if !(USE_HDG)
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
PUBLIC::DGTimeDerivative_weakForm
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/
#endif /*USE_HDG*/
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
USE MOD_Restart_Vars       ,ONLY: DoRestart,RestartInitIsDone
USE MOD_Interpolation_Vars ,ONLY: N_Inter,InterpolationInitIsDone,Nmax,Nmin
USE MOD_Mesh_Vars          ,ONLY: nSides,nElems, offSetElem
USE MOD_Mesh_Vars          ,ONLY: MeshInitIsDone
#if ! (USE_HDG)
USE MOD_PML_Vars           ,ONLY: PMLnVar ! Additional fluxes for the CFS-PML auxiliary variables
#endif /*USE_HDG*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
USE MOD_LoadBalance_Vars   ,ONLY: UseH5IOLoadBalance
#endif /*USE_LOADBALANCE*/
#if USE_MPI
USE MOD_MPI                ,ONLY: StartExchange_DG_Elems,FinishExchangeMPIData
!USE MOD_MPI_Vars           ,ONLY: SendRequest_U,RecRequest_U,SendRequest_U2,RecRequest_U2
#endif /*USE_MPI*/
#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)|| (PP_TimeDiscMethod==6)
USE MOD_TimeDisc_Vars      ,ONLY: Ut_N
USE MOD_PML_Vars           ,ONLY: DoPML,isPMLElem
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: Nloc,iElem,iSide
!===================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.(.NOT.MeshInitIsDone).OR.(.NOT.RestartInitIsDone).OR.DGInitIsDone) CALL abort(__STAMP__,&
    'InitDG not ready to be called or already called.')
LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT DG...'

! allocate arrays of pre-computed dg basis tensors
ALLOCATE(DGB_N(Nmin:Nmax))

DO Nloc=Nmin,Nmax
  ! Pre-compute the dg operator building blocks (differentiation matrices and prolongation operators)
  CALL InitDGBasis(Nloc, N_Inter(Nloc)%xGP, N_Inter(Nloc)%wGP, N_Inter(Nloc)%L_minus, N_Inter(Nloc)%L_plus, &
                   DGB_N(Nloc)%D, DGB_N(Nloc)%D_T, DGB_N(Nloc)%D_Hat, DGB_N(Nloc)%D_Hat_T, DGB_N(Nloc)%L_HatMinus, DGB_N(Nloc)%L_HatPlus )
END DO

#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
#if USE_LOADBALANCE && !(USE_HDG)
! Not "LB via MPI" means during 1st initialisation
IF (.NOT.(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance))) THEN
#endif /*USE_LOADBALANCE && !(USE_HDG)*/
  ! the local DG solution in physical and reference space
  ALLOCATE(U_N(1:nElems))
  DO iElem = 1, nElems
    Nloc = N_DG_Mapping(2,iElem+offSetElem)
    ALLOCATE(U_N(iElem)%U(PP_nVar,0:Nloc,0:Nloc,0:Nloc))
    U_N(iElem)%U = 0.
#if !(USE_HDG)
    ALLOCATE(U_N(iElem)%Ut(PP_nVar,0:Nloc,0:Nloc,0:Nloc))
    U_N(iElem)%Ut = 0.
    IF(DoPML)THEN
      IF(isPMLElem(iElem))THEN
        ALLOCATE(U_N(iElem)%U2(PMLnVar,0:Nloc,0:Nloc,0:Nloc))
        U_N(iElem)%U2 = 0.
        ALLOCATE(U_N(iElem)%U2t(PMLnVar,0:Nloc,0:Nloc,0:Nloc))
        U_N(iElem)%U2t = 0.
      END IF ! isPMLElem(iElem)
    END IF ! DoPML
#endif /*!(USE_HDG)*/
  END DO ! iElem = 1, nElems
#if USE_LOADBALANCE && !(USE_HDG)
END IF
#endif /*USE_LOADBALANCE && !(USE_HDG)*/

! Allocate additional containers
#if USE_HDG
DO iElem = 1, nElems
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  ALLOCATE(U_N(iElem)%E(1:3,0:Nloc,0:Nloc,0:Nloc))
  ALLOCATE(U_N(iElem)%Et(1:3,0:Nloc,0:Nloc,0:Nloc))
  U_N(iElem)%E = 0.
  U_N(iElem)%Et = 0.
END DO ! iElem = 1, nElems
#else
#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)|| (PP_TimeDiscMethod==6)
! the time derivative computed with the DG scheme
ALLOCATE(Ut_N(nElems))
#endif /*(PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)|| (PP_TimeDiscMethod==6)*/

! the time derivative computed with the DG scheme
DO iElem = 1, nElems
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)|| (PP_TimeDiscMethod==6)
  ALLOCATE(Ut_N(iElem)%Ut_temp(PP_nVar,0:Nloc,0:Nloc,0:Nloc))
  Ut_N(iElem)%Ut_temp = 0.
  IF(DoPML)THEN
    IF(isPMLElem(iElem))THEN
      ALLOCATE(Ut_N(iElem)%U2t_temp(PMLnVar,0:Nloc,0:Nloc,0:Nloc))
      Ut_N(iElem)%Ut_temp = 0.
    END IF ! isPMLElem(iElem)
  END IF ! DoPML
#endif
END DO ! iElem = 1, nElems
#endif /*USE_HDG*/

#if IMPA || ROS
ALLOCATE( Un(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
Un=0.
#endif
!nTotal_face = (PP_N+1)*(PP_N+1)
!nTotal_vol  = nTotal_face*(PP_N+1)
!nTotalU     = PP_nVar*nTotal_vol*nElems

! U is filled with the ini solution
IF(.NOT.DoRestart) CALL FillIni()

! We store the interior data at the each element face
!ALLOCATE(U_Minus(PP_nVar,0:PP_N,0:PP_N,sideID_minus_lower:sideID_minus_upper))
!ALLOCATE(U_Plus(PP_nVar,0:PP_N,0:PP_N,sideID_plus_lower:sideID_plus_upper))
!U_Minus=0.
!U_Plus=0.
ALLOCATE(U_Surf_N(1:nSides))

DO iSide = 1, nSides
  Nloc = DG_Elems_master(iSide)
  ALLOCATE(U_Surf_N(iSide)%U_master(1:PP_nVar,0:Nloc,0:Nloc))
  Nloc = DG_Elems_slave(iSide)
  ALLOCATE(U_Surf_N(iSide)%U_slave( 1:PP_nVar,0:Nloc,0:Nloc))
  U_Surf_N(iSide)%U_master = 0.
  U_Surf_N(iSide)%U_slave  = 0.
END DO ! iSide = 1, nSides

#if !(USE_HDG)
! unique flux per side
! additional fluxes for the CFS-PML auxiliary variables (no PML: PMLnVar=0)
! additional fluxes for the CFS-PML auxiliary variables (no PML: PMLnVar=0)

DO iSide = 1, nSides
  Nloc = DG_Elems_master(iSide)
  ALLOCATE(U_Surf_N(iSide)%Flux_Master(1:PP_nVar+PMLnVar,0:Nloc,0:Nloc))
  Nloc = DG_Elems_slave(iSide)
  ALLOCATE(U_Surf_N(iSide)%Flux_Slave( 1:PP_nVar+PMLnVar,0:Nloc,0:Nloc))
  U_Surf_N(iSide)%Flux_Master = 0.
  U_Surf_N(iSide)%Flux_Slave  = 0.
END DO ! iSide = 1, nSides
#endif /*USE_HDG*/

#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/
DGInitIsDone=.TRUE.
LBWRITE(UNIT_stdOut,'(A)')' INIT DG DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitDG


!==================================================================================================================================
!> Allocate and initialize the building blocks for the DG operator: Differentiation matrices and prolongation operators
!==================================================================================================================================
SUBROUTINE InitDGBasis(N_in,xGP,wGP,L_Minus,L_Plus,D,D_T,D_Hat,D_Hat_T,L_HatMinus,L_HatPlus)
! MODULES
USE MOD_Globals
USE MOD_Basis              ,ONLY: LegendreGaussNodesAndWeights,LegGaussLobNodesAndWeights,BarycentricWeights
USE MOD_Basis              ,ONLY: PolynomialDerivativeMatrix,LagrangeInterpolationPolys
#if USE_HDG
USE MOD_Interpolation_Vars ,ONLY: Nmax

#if USE_MPI
USE MOD_PreProc
USE MOD_MPI_Vars           ,ONLY: SurfExchange, nNbProcs, DataSizeSurfRecMax, DataSizeSurfSendMax, DataSizeSurfRecMin, DataSizeSurfSendMin
USE MOD_DG_Vars            ,ONLY: N_DG_Mapping,DG_Elems_master,DG_Elems_slave
USE MOD_MPI                ,ONLY: StartReceiveMPISurfDataType,StartSendMPISurfDataType,FinishExchangeMPISurfDataType
USE MOD_Mesh_Vars          ,ONLY: N_SurfMesh,nSides, offSetElem
USE MOD_Interpolation_Vars ,ONLY: NInfo,PREF_VDM,N_Inter
#endif /*USE_MPI*/
#endif /*USE_HDG*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                             :: N_in                   !< Polynomial degree
REAL,DIMENSION(0:N_in),INTENT(IN)              :: xGP                    !< Gauss/Gauss-Lobatto Nodes
REAL,DIMENSION(0:N_in),INTENT(IN)              :: wGP                    !< Gauss/Gauss-Lobatto Weights
REAL,DIMENSION(0:N_in),INTENT(IN)              :: L_Minus                !< Values of lagrange polynomials at \f$ \xi = -1 \f$
REAL,DIMENSION(0:N_in),INTENT(IN)              :: L_Plus                 !< Values of lagrange polynomials at \f$ \xi = +1 \f$
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT)    :: D                      !< Differentation matrix
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT)    :: D_T                    !< Transpose of differentation matrix
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT)    :: D_Hat                  !< Differentiation matrix premultiplied by mass matrix,
                                                                         !< \f$ \hat{D} = M^{-1} D^T M \f$
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT)    :: D_Hat_T                !< Transpose of D_Hat matrix \f$ \hat{D}^T \f$
REAL,ALLOCATABLE,DIMENSION(:)  ,INTENT(OUT)    :: L_HatMinus             !< Values of lagrange polynomials at \f$ \xi = -1 \f$
                                                                         !< premultiplied with mass matrix
REAL,ALLOCATABLE,DIMENSION(:)  ,INTENT(OUT)    :: L_HatPlus              !< Values of lagrange polynomials at \f$ \xi = +1 \f$
                                                                         !< premultiplied with mass matrix
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(0:N_in,0:N_in)              :: M,Minv
INTEGER                                    :: iMass
#if USE_HDG
#if USE_MPI
!REAL                                       :: Geotemp(10,0:PP_N,0:PP_N,1:nSides)
INTEGER                 :: iSide, i,nRecVal, nSendVal, Nloc, iNbProc, NSideMax, NSideMin, p, q, SideID_end, SideID_start
#endif /*USE_MPI*/
#endif /*USE_HDG*/
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
D_Hat_T    = TRANSPOSE(D_hat)

! interpolate to left and right face (1 and -1) and pre-divide by mass matrix
L_HatPlus(:)  = MATMUL(Minv,L_Plus )
L_HatMinus(:) = MATMUL(Minv,L_Minus)

#if USE_HDG
! exchange is in InitDGBasis as InitMesh() and InitMPI() is needed
!CALL abort(__STAMP__,'not implemented: but is it actually required?')
!Geotemp=0.
!Geotemp(1,:,:,:)    = N_SurfMesh(?)%SurfElem(:,:)
!Geotemp(2:4,:,:,:)  = N_SurfMesh(?)%NormVec(:,:,:)
!Geotemp(5:7,:,:,:)  = N_SurfMesh(?)%TangVec1(:,:,:)
!Geotemp(8:10,:,:,:) = N_SurfMesh(?)%TangVec2(:,:,:)
!!Geotemp(11:13,:,:,:)=Face_xGP(:,:,:,SideID_minus_lower:SideID_minus_upper)
!CALL StartReceiveMPIData(10,Geotemp,1,nSides,RecRequest_Geo ,SendID=1) ! Receive MINE
!CALL StartSendMPIData(   10,Geotemp,1,nSides,SendRequest_Geo,SendID=1) ! Send YOUR
!CALL FinishExchangeMPIData(SendRequest_Geo,RecRequest_Geo,SendID=1)

#endif /*USE_HDG*/
END SUBROUTINE InitDGBasis


#if !(USE_HDG)
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
SUBROUTINE DGTimeDerivative_weakForm(t,tStage,tDeriv,doSource)
!===================================================================================================================================
! Computes the DG time derivative consisting of Volume Integral and Surface integral for the whole field
! U and Ut are allocated
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Vector
USE MOD_DG_Vars           ,ONLY: U_N
USE MOD_SurfInt           ,ONLY: SurfInt
USE MOD_VolInt            ,ONLY: VolInt
USE MOD_ProlongToFace     ,ONLY: ProlongToFace_TypeBased
USE MOD_FillFlux          ,ONLY: FillFlux
USE MOD_Equation          ,ONLY: CalcSource
USE MOD_Interpolation     ,ONLY: ApplyJacobian
USE MOD_PML_Vars          ,ONLY: DoPML
USE MOD_Mesh_Vars         ,ONLY: nElems
!USE MOD_FillMortar        ,ONLY: U_Mortar,Flux_Mortar
#if USE_MPI
!USE MOD_PML_Vars          ,ONLY: PMLnVar
!USE MOD_Mesh_Vars         ,ONLY: nSides
USE MOD_MPI_Vars
USE MOD_MPI                ,ONLY: StartExchange_DG_Elems,StartReceiveMPIDataType,StartSendMPIDataType,FinishExchangeMPIDataType
#if defined(PARTICLES) && defined(LSERK)
USE MOD_Particle_Vars      ,ONLY: DelayTime
USE MOD_TimeDisc_Vars      ,ONLY: time
#endif /*defined(PARTICLES) && defined(LSERK)*/
#ifdef PARTICLES
USE MOD_Particle_MPI       ,ONLY: MPIParticleSend,MPIParticleRecv
#endif /*PARTICLES*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers ,ONLY: LBStartTime,LBPauseTime,LBSplitTime
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/
USE MOD_PML_Vars           ,ONLY: nPMLElems,PMLToElem
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
#if USE_LOADBALANCE
REAL                            :: tLBStart
#endif /*USE_LOADBALANCE*/
INTEGER                         :: iElem,iPML
!===================================================================================================================================

! prolong the solution to the face integration points for flux computation
#if USE_MPI
! Prolong to face for MPI sides - send direction
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
CALL StartReceiveMPIDataType(RecRequest_U,SendID=2) ! Receive MINE
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DGCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
CALL ProlongToFace_TypeBased(doMPISides=.TRUE.)
!CALL U_Mortar(U_master,U_slave,doMPISides=.TRUE.)
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/
CALL StartSendMPIDataType(SendRequest_U,SendID=2) ! Send YOUR
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DGCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/

! Prolong to face for BCSides, InnerSides and MPI sides - receive direction
CALL ProlongToFace_TypeBased(doMPISides=.FALSE.)
!CALL U_Mortar(U_master,U_slave,doMPISides=.FALSE.)

#if USE_MPI
#if defined(PARTICLES) && defined(LSERK)
IF (time.GE.DelayTime) THEN
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_INTERPOLATION,tLBStart)
#endif /*USE_LOADBALANCE*/
#if USE_MPI
  CALL MPIParticleSend()
#endif /*USE_MPI*/
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
END IF
#endif /*defined(PARTICLES) && defined(LSERK)*/
#endif /*USE_MPI*/

! Nullify arrays
! NOTE: IF NEW DG_VOLINT AND LIFTING_VOLINT ARE USED AND CALLED FIRST,
!       ARRAYS DO NOT NEED TO BE NULLIFIED, OTHERWISE THEY HAVE TO!
!CALL VNullify(nTotalU,Ut)
DO iElem = 1, nElems
  U_N(iElem)%Ut = 0.
END DO ! iElem = 1, nElems

IF(DoPML)THEN ! Set U2t for auxiliary variables to zero
  DO iPML=1,nPMLElems
    iElem = PMLToElem(iPML)
    U_N(iElem)%U2t = 0.
  END DO ! iPML=1,nPMLElems
END IF ! DoPML
! compute volume integral contribution and add to ut, first half of all elements
CALL VolInt(dofirstElems=.TRUE.)

#if USE_MPI
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/
! Complete send / receive
CALL FinishExchangeMPIDataType(SendRequest_U,RecRequest_U,SendID=2) !Send YOUR - receive MINE

! Initialization of the time derivative
!Flux=0. !don't nullify the fluxes if not really needed (very expensive)
!CALL StartReceiveMPIDataType(PP_nVar+PMLnVar,Flux_Slave,1,nSides,RecRequest_Flux,SendID=1) ! Receive MINE
CALL StartReceiveMPIDataType(RecRequest_Flux,SendID=1) ! Receive MINE
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DGCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
! fill the global surface flux list
CALL FillFlux(t,tDeriv,doMPISides=.TRUE.)
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/

!CALL StartSendMPIData(PP_nVar+PMLnVar,Flux_Slave,1,nSides,SendRequest_Flux,SendID=1) ! Send YOUR
CALL StartSendMPIDataType(SendRequest_Flux,SendID=1) ! Send YOUR
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DGCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/

! fill the all surface fluxes on this proc
CALL FillFlux(t,tDeriv,doMPISides=.FALSE.)
!CALL Flux_Mortar(Flux_Master,Flux_Slave,doMPISides=.FALSE.)
! compute surface integral contribution and add to ut
CALL SurfInt(doMPISides=.FALSE.)

! compute volume integral contribution and add to ut
CALL VolInt(dofirstElems=.FALSE.)

#if USE_MPI
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/
! Complete send / receive
CALL FinishExchangeMPIDataType(SendRequest_Flux,RecRequest_Flux,SendID=1) !Send MINE -receive YOUR
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DGCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/

!FINALIZE Fluxes for MPI Sides
!CALL Flux_Mortar(Flux_Master,Flux_Slave,doMPISides=.TRUE.)
CALL SurfInt(doMPISides=.TRUE.)
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif

! swap and map to physical space
CALL ApplyJacobian(toPhysical=.TRUE.,toSwap=.TRUE.)

! Add Source Terms
IF(doSource) CALL CalcSource(tStage,1.0)

#if USE_LOADBALANCE
CALL LBSplitTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/

#if defined(PARTICLES) && defined(LSERK)
#if USE_MPI
IF (time.GE.DelayTime) THEN
  CALL MPIParticleRecv()
END IF
#if USE_LOADBALANCE
CALL LBSplitTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/
#endif /*defined(PARTICLES) && defined(LSERK)*/

END SUBROUTINE DGTimeDerivative_weakForm
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/
#endif /*!(USE_HDG)*/


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
USE MOD_Equation           ,ONLY: VolInt_Pois,FillFlux_Pois,ProlongToFace_Pois, SurfInt_Pois
USE MOD_GetBoundaryFlux    ,ONLY: FillFlux_BC_Pois
USE MOD_Mesh_Vars          ,ONLY: sJ,Elem_xGP,nSides
USE MOD_Equation           ,ONLY: CalcSource_Pois
USE MOD_Equation_Vars      ,ONLY: IniExactFunc,Phi,Phit,Phi_master,Phi_slave,FluxPhi,nTotalPhi
USE MOD_Interpolation      ,ONLY: ApplyJacobian
#if USE_MPI
USE MOD_MPI_Vars
USE MOD_MPI                ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers ,ONLY: LBStartTime,LBPauseTime,LBSplitTime
#endif /*USE_LOADBALANCE*/
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
#if USE_MPI
REAL    :: tLBStart
#endif /*USE_MPI*/
!===================================================================================================================================

! prolong the solution to the face integration points for flux computation
#if USE_MPI
! Prolong to face for MPI sides - send direction
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
!CALL ProlongToFace(Phi,Phi_Minus,Phi_slave,doMPiSides=.TRUE.)
CALL StartReceiveMPIData(4,Phi_slave,1,nSides,RecRequest_U,SendID=2) ! Receive MINE
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DGCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
CALL ProlongToFace_Pois(Phi,Phi_master,Phi_slave,doMPiSides=.TRUE.)
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/

!CALL StartExchangeMPIData(Phi_slave,SideID_plus_lower,SideID_plus_upper,SendRequest_U,RecRequest_U,SendID=2)
CALL StartSendMPIData(4,Phi_slave,1,nSides,SendRequest_U,SendID=2) ! Send YOUR
! Send YOUR - receive MINE
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DGCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/

! Prolong to face for BCSides, InnerSides and MPI sides - receive direction
!CALL ProlongToFace(Phi,Phi_Minus,Phi_slave,doMPISides=.FALSE.)
CALL ProlongToFace_Pois(Phi,Phi_master,Phi_slave,doMPISides=.FALSE.)

Phit=0.
CALL VolInt_Pois(Phit)
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/


#if USE_MPI
! Complete send / receive
CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=2) !Send YOUR - receive MINE


! Initialization of the time derivative
!Flux=0. !don't nullify the fluxes if not really needed (very expensive)
! fill the global surface flux list
CALL StartReceiveMPIData(4,FluxPhi,1,nSides,RecRequest_Flux,SendID=1) ! Receive MINE
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DGCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
CALL FillFlux_Pois(FluxPhi,doMPISides=.TRUE.)
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/

!CALL StartExchangeMPIData(FluxPhi,1,nSides,SendRequest_Flux,RecRequest_Flux,SendID=1) ! Send MINE - receive YOUR
CALL StartSendMPIData(4,FluxPhi,1,nSides,SendRequest_Flux,SendID=1) ! Send YOUR
!CALL StartExchangeMPIData(4,FluxPhi,1,nSides,SendRequest_Flux,RecRequest_Flux,SendID=1)
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DGCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/

! fill the all surface fluxes on this proc
CALL FillFlux_BC_Pois(t,tDeriv,FluxPhi)
CALL FillFlux_Pois(FluxPhi,doMPISides=.FALSE.)
! compute surface integral contribution and add to ut
CALL SurfInt_Pois(FluxPhi,Phit,doMPISides=.FALSE.)
!! compute volume integral contribution and add to ut
!CALL VolInt(Ut)

#if USE_MPI
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/
! Complete send / receive
CALL FinishExchangeMPIData(SendRequest_Flux,RecRequest_Flux,SendID=1) !Send MINE -receive YOUR
#if USE_LOADBALANCE
CALL LBSplitTime(LB_DGCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
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
#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart)
#endif /*USE_LOADBALANCE*/

END SUBROUTINE DGTimeDerivative_weakForm_Pois

#endif

SUBROUTINE FillIni()
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars       ,ONLY: U_N,N_DG_Mapping
USE MOD_Mesh_Vars     ,ONLY: N_VolMesh, offSetElem
USE MOD_Equation_Vars ,ONLY: IniExactFunc
USE MOD_Equation      ,ONLY: ExactFunc
#ifdef maxwell
USE MOD_Equation_Vars ,ONLY: DoExactFlux
#endif /*maxwell*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iElem,Nloc
!===================================================================================================================================
! Determine Size of the Loops, i.e. the number of grid cells in the
! corresponding directions
#ifdef maxwell
IF(DoExactFlux.AND.(IniExactFunc.NE.16)) RETURN ! IniExactFunc=16 is pulsed laser mixed IC+BC
#endif /*maxwell*/
DO iElem=1,PP_nElems
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  DO k=0,Nloc
    DO j=0,Nloc
      DO i=0,Nloc
#if USE_HDG
        CALL ExactFunc(IniExactFunc,     N_VolMesh(iElem)%Elem_xGP(1:3,i,j,k),U_N(iElem)%U(1:PP_nVar,i,j,k),ElemID=iElem)
#else
        CALL ExactFunc(IniExactFunc,0.,0,N_VolMesh(iElem)%Elem_xGP(1:3,i,j,k),U_N(iElem)%U(1:PP_nVar,i,j,k))
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
USE MOD_globals          ,ONLY: abort
USE MOD_DG_Vars
#if USE_LOADBALANCE && ! (USE_HDG)
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
#endif /*USE_LOADBALANCE && ! (USE_HDG)*/
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
#if !(USE_HDG)
USE MOD_TimeDisc_Vars    ,ONLY: Ut_N
#endif /*!(USE_HDG)*/
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(DGB_N)
#if IMPA || ROS
SDEALLOCATE(Un)
#endif
SDEALLOCATE(U_Surf_N)

! Do not deallocate the solution vector during load balance here as it needs to be communicated between the processors
#if USE_LOADBALANCE && !(USE_HDG)
IF(.NOT.(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)))THEN
#endif /*USE_LOADBALANCE && !(USE_HDG)*/
  ! Keep for load balance and deallocate/reallocate after communication
  SDEALLOCATE(U_N)
#if USE_LOADBALANCE && !(USE_HDG)
END IF
#endif /*USE_LOADBALANCE && !(USE_HDG)*/

#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
#if !(USE_HDG)
SDEALLOCATE(Ut_N)
#endif /*!(USE_HDG)*/
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/

DGInitIsDone = .FALSE.
END SUBROUTINE FinalizeDG

END MODULE MOD_DG
