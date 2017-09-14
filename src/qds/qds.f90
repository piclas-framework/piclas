#include "boltzplatz.h"

MODULE MOD_QDS
!===================================================================================================================================
!> Contains the routines to
!> - 
!===================================================================================================================================
! MODULES
!USE MOD_io_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitQDS
  MODULE PROCEDURE InitQDS
END INTERFACE
INTERFACE QDSTimeDerivative
  MODULE PROCEDURE QDSTimeDerivative
END INTERFACE
INTERFACE FinalizeQDS
  MODULE PROCEDURE FinalizeQDS
END INTERFACE
INTERFACE QDSReCalculateDGValues
  MODULE PROCEDURE QDSReCalculateDGValues
END INTERFACE
INTERFACE QDSCalculateMacroValues
  MODULE PROCEDURE QDSCalculateMacroValues
END INTERFACE

PUBLIC::InitQDS
PUBLIC::FinalizeQDS
PUBLIC::QDSTimeDerivative
PUBLIC::QDSReCalculateDGValues
PUBLIC::QDSCalculateMacroValues
!===================================================================================================================================
CONTAINS


SUBROUTINE InitQDS
!===================================================================================================================================
!> Allocate all QDS variables, determine
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_QDS_Vars
USE MOD_Particle_Vars,   ONLY:BoltzmannConst
USE MOD_Globals,         ONLY:abort,myrank,UNIT_stdOut,mpiroot,iError
USE MOD_ReadInTools,     ONLY:GETLOGICAL
USE MOD_Mesh_Vars,       ONLY:nSides
USE MOD_Globals_Vars,    ONLY:PI
USE MOD_Restart_Vars,    ONLY:DoRestart
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: tempNorm
INTEGER :: iWeight,i,j,k
REAL    :: Velo(3), Temp, Dens, Mass
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT QDS...' 
DoQDS                      = GETLOGICAL('DoQDS','.FALSE.')
IF(DoQDS)THEN
  QDSnVar=40
  nQDSElems=PP_nElems

  ALLOCATE(UQDS       (1:QDSnVar,0:PP_N,0:PP_N,0:PP_N,1:nQDSElems))        
  ALLOCATE(UQDSt      (1:QDSnVar,0:PP_N,0:PP_N,0:PP_N,1:nQDSElems))
  UQDS=0.
  UQDSt=0.

  ALLOCATE(UQDS_master(1:QDSnVar,0:PP_N,0:PP_N,1:nSides))
  ALLOCATE(UQDS_slave( 1:QDSnVar,0:PP_N,0:PP_N,1:nSides))
  UQDS_master=0.
  UQDS_slave=0.

  ALLOCATE(FluxQDS_Master(1:QDSnVar,0:PP_N,0:PP_N,1:nSides)) 
  ALLOCATE(FluxQDS_Slave( 1:QDSnVar,0:PP_N,0:PP_N,1:nSides)) 
  FluxQDS_Master=0.
  FluxQDS_Slave=0.

  ALLOCATE(QDSMacroValues(6,0:PP_N,0:PP_N,0:PP_N,nQDSElems))
  QDSMacroValues=0.

  ALLOCATE(GaussHermitWeiAbs(2,2))
  GaussHermitWeiAbs=0.

  ! Weights
  GaussHermitWeiAbs(1,1:2) = 0.8862269254527580136491 ! 0.886227
  GaussHermitWeiAbs(1,:) = GaussHermitWeiAbs(1,:)*SQRT(PI)/SUM(GaussHermitWeiAbs(1,:))
  !Absisc
  GaussHermitWeiAbs(2,1) = -0.7071067811865475244008 !-0.707107
  GaussHermitWeiAbs(2,2) =  0.7071067811865475244008 ! 0.707107
  tempNorm = 0.0
  DO iWeight = 1, 2
    tempNorm = tempNorm + GaussHermitWeiAbs(1,iWeight)*GaussHermitWeiAbs(2,iWeight)**2
  END DO
  GaussHermitWeiAbs(2,:)= GaussHermitWeiAbs(2,:)*SQRT(SQRT(PI)/(2.*tempNorm))

  Temp = 278.687
  Velo =(/0.,0.,0./) !(/1020.882,0.,0./)
  Dens= 2.633459376E25
  Mass = 4.651734101E-26
  QDSSpeciesMass=Mass

  !QDS_Species = GETINT('Particles-QDSSpecies','0')

  QDSSpecDOF=3
  QDSMaxVelo=2*(ABS(Velo(1))+SQRT(2.*BoltzmannConst*Temp/Mass)*ABS(GaussHermitWeiAbs(2,1)))



  !IF(.NOT.DoRestart) CALL FillIniQDS()

  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
  !  i=0; j=0; k=0
  !  QDSMacroValues(1,i,j,k,1) = Dens*wGP(i)*wGP(j)*wGP(k)/sJ(i,j,k,1)*Mass
    QDSMacroValues(1,i,j,k,88) = Dens*Mass
    QDSMacroValues(2:4,i,j,k,88) = QDSMacroValues(1,i,j,k,88)*Velo(1:3)
    QDSMacroValues(6,i,j,k,88) = Temp
  END DO; END DO; END DO


ELSE
  QDSnVar=0
END IF
QDSInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT QDS DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitQDS


SUBROUTINE QDSTimeDerivative(t,tStage,tDeriv,doSource)
!===================================================================================================================================
! Computes the DG time derivative consisting of Volume Integral and Surface integral for the whole field
! UQDS and UQDSt are allocated
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_QDS_Vars,      ONLY: UQDS,UQDSt
USE MOD_QDS_Vars,      ONLY: nQDSElems,QDSnVar
USE MOD_Mesh_Vars,     ONLY: sJ
USE MOD_Vector
USE MOD_QDS_Vars,         ONLY:UQDS,UQDSt,UQDS_master,UQDS_Slave,FluxQDS_Master,FluxQDS_Slave !,nTotalU
USE MOD_ProlongToFace,    ONLY:ProlongToFaceQDS
USE MOD_Mesh_Vars,        ONLY:nSides
USE MOD_Interpolation,    ONLY:ApplyJacobianQDS
#ifdef MPI
USE MOD_MPI_Vars
USE MOD_MPI,              ONLY:StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
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
INTEGER                         :: i,j,k,iQDSElem,iQDSVar
!===================================================================================================================================

! prolong the solution to the face integration points for flux computation
#ifdef MPI
! Prolong to face for MPI sides - send direction
CALL StartReceiveMPIData(QDSnVar,UQDS_Slave,1,nSides,RecRequest_U,SendID=2) ! Receive MINE
CALL ProlongToFaceQDS(UQDS,UQDS_Master,UQDS_Slave,doMPISides=.TRUE.)
!CALL U_Mortar(UQDS_Master,UQDS_Slave,doMPISides=.TRUE.)
CALL StartSendMPIData(QDSnVar,UQDS_Slave,1,nSides,SendRequest_U,SendID=2) ! Send YOUR
#endif /*MPI*/

! Prolong to face for BCSides, InnerSides and MPI sides - receive direction
CALL ProlongToFaceQDS(UQDS,UQDS_Master,UQDS_Slave,doMPISides=.FALSE.)
!CALL U_Mortar(UQDS_Master,UQDS_Slave,doMPISides=.FALSE.)
! Nullify arrays
! NOTE: IF NEW DG_VOLINT AND LIFTING_VOLINT ARE USED AND CALLED FIRST,
!       ARRAYS DO NOT NEED TO BE NULLIFIED, OTHERWISE THEY HAVE TO!
!CALL VNullify(nTotalU,UQDSt)
UQDSt=0.
! compute volume integral contribution and add to ut, first half of all elements
CALL VolIntQDS(UQDSt,dofirstElems=.TRUE.)

#ifdef MPI
! Complete send / receive
CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=2) !Send YOUR - receive MINE
#endif /*MPI*/

! Initialization of the time derivative
!Flux=0. !don't nullify the fluxes if not really needed (very expensive)
#ifdef MPI
CALL StartReceiveMPIData(QDSnVar,FluxQDS_Slave,1,nSides,RecRequest_Flux,SendID=1) ! Receive MINE
! fill the global surface flux list
CALL FillFluxQDS(t,tDeriv,FluxQDS_Master,FluxQDS_Slave,UQDS_Master,UQDS_Slave,doMPISides=.TRUE.)

CALL StartSendMPIData(QDSnVar,FluxQDS_Slave,1,nSides,SendRequest_Flux,SendID=1) ! Send YOUR
#endif /* MPI*/

! fill the all surface fluxes on this proc
CALL FillFluxQDS(t,tDeriv,FluxQDS_Master,FluxQDS_Slave,UQDS_Master,UQDS_Slave,doMPISides=.FALSE.)
!CALL Flux_Mortar(FluxQDS_Master,FluxQDS_Slave,doMPISides=.FALSE.)
! compute surface integral contribution and add to ut
CALL SurfIntQDS(FluxQDS_Master,FluxQDS_Slave,UQDSt,doMPISides=.FALSE.)

! compute volume integral contribution and add to ut
CALL VolIntQDS(UQDSt,dofirstElems=.FALSE.)

#ifdef MPI
! Complete send / receive
CALL FinishExchangeMPIData(SendRequest_Flux,RecRequest_Flux,SendID=1) !Send MINE -receive YOUR

!FINALIZE Fluxes for MPI Sides
!CALL Flux_Mortar(FluxQDS_Master,FluxQDS_Slave,doMPISides=.TRUE.)
CALL SurfIntQDS(FluxQDS_Master,FluxQDS_Slave,UQDSt,doMPISides=.TRUE.)
#endif

! swap and map to physical space
CALL ApplyJacobianQDS(UQDSt,toPhysical=.TRUE.,toSwap=.TRUE.)

! Add Source Terms
!IF(doSource) CALL CalcSourceQDS(tStage,1.0,Ut)


END SUBROUTINE QDSTimeDerivative


SUBROUTINE FillIniQDS()
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,      ONLY:Elem_xGP
!USE MOD_Equation_Vars, ONLY:IniExactFunc
USE MOD_QDS_Vars,       ONLY:nQDSElems,UQDS,QDSnVar
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iQDSElems
!===================================================================================================================================
DO iQDSElems=1,nQDSElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
#ifdef PP_HDG
        !CALL ExactFuncQDS(IniExactFunc,     Elem_xGP(1:3,i,j,k,iQDSElems),UQDS(1:QDSnVar,i,j,k,iQDSElems))
        CALL ExactFuncQDS(1,     Elem_xGP(1:3,i,j,k,iQDSElems),UQDS(1:QDSnVar,i,j,k,iQDSElems))
#else
        !CALL ExactFuncQDS(IniExactFunc,0.,0,Elem_xGP(1:3,i,j,k,iQDSElems),UQDS(1:QDSnVar,i,j,k,iQDSElems))
        CALL ExactFuncQDS(1,0.,0,Elem_xGP(1:3,i,j,k,iQDSElems),UQDS(1:QDSnVar,i,j,k,iQDSElems))
#endif
      END DO ! i
    END DO ! j
  END DO !k
END DO ! iQDSElems=1,nQDSElems
END SUBROUTINE FillIniQDS


SUBROUTINE ExactFuncQDS(ExactFunction,t,tDeriv,x,resu) 
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_QDS_Vars,        ONLY:QDSnVar
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: t
INTEGER,INTENT(IN)              :: tDeriv           ! determines the time derivative of the function
REAL,INTENT(IN)                 :: x(3)              
INTEGER,INTENT(IN)              :: ExactFunction    ! determines the exact function
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Resu(1:QDSnVar)    ! state in conservative variables
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                            :: Resu_t(1:QDSnVar),Resu_tt(1:QDSnVar) ! state in conservative variables
REAL                            :: Cent(3)
!===================================================================================================================================
Cent=x
SELECT CASE (ExactFunction)
#ifdef PARTICLES
CASE(0) ! Particles
  Resu=0.
  !resu(1:3)= x(1:3)!*x(1) 
#endif
CASE(1) ! Constant 
  Resu=1.
  Resu_t=0.
  Resu_tt=0.

CASE DEFAULT
  SWRITE(*,*)'Exact function not specified'
END SELECT ! ExactFunction



# if (PP_TimeDiscMethod==1)
! For O3 RK, the boundary condition has to be adjusted
! Works only for O3 RK!!
SELECT CASE(tDeriv)
CASE(0)
  ! resu = g(t)
CASE(1)
  ! resu = g(t) + dt/3*g'(t)
  Resu=Resu + dt/3.*Resu_t
CASE(2)
  ! resu = g(t) + 3/4 dt g'(t) +5/16 dt^2 g''(t)
  Resu=Resu + 0.75*dt*Resu_t+5./16.*dt*dt*Resu_tt
CASE DEFAULT
  ! Stop, works only for 3 Stage O3 LS RK
  CALL abort(&
      __STAMP__&
      ,'Exactfuntion works only for 3 Stage O3 LS RK!',999,999.)
END SELECT
#endif
END SUBROUTINE ExactFuncQDS


SUBROUTINE FillFluxQDS(t,tDeriv,FluxQDS_Master,FluxQDS_Slave,UQDS_Master,UQDS_Slave,doMPISides)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_GLobals
USE MOD_PreProc
USE MOD_Mesh_Vars,       ONLY:NormVec,SurfElem
USE MOD_Mesh_Vars,       ONLY:nSides,nBCSides
USE MOD_Mesh_Vars,       ONLY:NormVec,TangVec1, tangVec2, SurfElem,Face_xGP
USE MOD_Mesh_Vars,       ONLY:firstMPISide_MINE,lastMPISide_MINE,firstInnerSide,firstBCSide,lastInnerSide
USE MOD_Mesh_Vars,       ONLY:SideToElem
USE MOD_QDS_Vars,        ONLY:QDSnVar
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  != .TRUE. only MINE MPISides are filled, =.FALSE. InnerSides  
REAL,INTENT(IN)    :: t           ! time
INTEGER,INTENT(IN) :: tDeriv      ! deriv
REAL,INTENT(IN)    :: UQDS_Master(QDSnVar,0:PP_N,0:PP_N,1:nSides)
REAL,INTENT(IN)    :: UQDS_slave (QDSnVar,0:PP_N,0:PP_N,1:nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: FluxQDS_Master(1:QDSnVar,0:PP_N,0:PP_N,nSides)
REAL,INTENT(OUT)   :: FluxQDS_Slave(1:QDSnVar,0:PP_N,0:PP_N,nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,p,q,firstSideID_wo_BC,firstSideID ,lastSideID,ElemID
!===================================================================================================================================

! fill flux for sides ranging between firstSideID and lastSideID using RiemannQDS solver
IF(doMPISides)THEN 
  ! fill only flux for MINE MPISides
  ! fill only flux for MINE MPISides (where the local proc is master) 
  firstSideID_wo_BC = firstMPISide_MINE
  firstSideID = firstMPISide_MINE
  lastSideID =  lastMPISide_MINE
ELSE
  ! fill only InnerSides
  ! fill only InnerSides that do not need communication
  firstSideID_wo_BC = firstInnerSide ! for fluxes
  firstSideID = firstBCSide    ! include BCs for master sides
  lastSideID = lastInnerSide
END IF

! Compute fluxes on PP_N, no additional interpolation required
DO SideID=firstSideID,lastSideID
  CALL RiemannQDS(FluxQDS_Master(1:QDSnVar,:,:,SideID),UQDS_Master( :,:,:,SideID),UQDS_Slave(  :,:,:,SideID),NormVec(:,:,:,SideID))
END DO ! SideID
  
IF(.NOT.doMPISides)THEN
  CALL GetBoundaryFluxQDS(t,tDeriv,FluxQDS_Master (1:QDSnVar,0:PP_N,0:PP_N,1:nBCSides) &
                                  ,UQDS_Master    (1:QDSnVar,0:PP_N,0:PP_N,1:nBCSides) &
                                  ,NormVec        (1:3      ,0:PP_N,0:PP_N,1:nBCSides) &
                                  ,TangVec1       (1:3      ,0:PP_N,0:PP_N,1:nBCSides) &
                                  ,TangVec2       (1:3      ,0:PP_N,0:PP_N,1:nBCSides) &
                                  ,Face_XGP       (1:3      ,0:PP_N,0:PP_N,1:nBCSides) )
END IF

! Apply surface element size
DO SideID=firstSideID,lastSideID
  DO q=0,PP_N; DO p=0,PP_N
    FluxQDS_Master(:,p,q,SideID)=FluxQDS_Master(:,p,q,SideID)*SurfElem(p,q,SideID)
  END DO; END DO
END DO

! copy flux from Master side to slave side, DO not change sign
FluxQDS_slave(:,:,:,firstSideID:lastSideID) = FluxQDS_master(:,:,:,firstSideID:lastSideID)


END SUBROUTINE FillFluxQDS


SUBROUTINE VolIntQDS(Ut,dofirstElems)
!===================================================================================================================================
! Computes the volume integral of the weak DG form a la Kopriva
! Attention 1: 1/J(i,j,k) is not yet accounted for
! Attention 2: ut is initialized and is updated with the volume flux derivatives
!===================================================================================================================================
! MODULES
USE MOD_DG_Vars,           ONLY:D_hat
USE MOD_Mesh_Vars,         ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Dielectric_Vars,   ONLY: DoDielectric,isDielectricElem
USE MOD_PreProc
USE MOD_Flux,ONLY:EvalFlux3D,EvalFlux3DDielectric                      ! computes volume fluxes in local coordinates
USE MOD_QDS_Vars,        ONLY:QDSnVar,nQDSElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                                  :: Ut(QDSnVar,0:PP_N,0:PP_N,0:PP_N,1:nQDSElems)
LOGICAL,INTENT(IN)                                  :: dofirstElems
! Adds volume contribution to time derivative Ut contained in MOD_DG_Vars (=aufschmutzen!)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(QDSnVar,0:PP_N,0:PP_N,0:PP_N)      :: f,g,h                ! volume fluxes at all Gauss points
REAL,DIMENSION(QDSnVar)                           :: fTilde,gTilde,hTilde ! auxiliary variables needed to store the fluxes at one GP
INTEGER                                           :: i,j,k,iElem
INTEGER                                           :: l                    ! row index for matrix vector product
INTEGER                                           :: firstElemID, lastElemID
!===================================================================================================================================

IF(dofirstElems)THEN
  firstElemID = 1
  lastElemID  = nQDSElems/2+1
ELSE ! second half of elements
  firstElemID = nQDSElems/2+2
  lastElemID  = nQDSElems
END IF

DO iElem=firstElemID,lastElemID
!DO iElem=1,nQDSElems
  ! Cut out the local DG solution for a grid cell iElem and all Gauss points from the global field
  ! Compute for all Gauss point values the Cartesian flux components
  CALL EvalFlux3DQDS(iElem,f,g,h)

  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        fTilde=f(:,i,j,k)
        gTilde=g(:,i,j,k)
        hTilde=h(:,i,j,k)
        ! Compute the transformed fluxes with the metric terms
        ! Attention 1: we store the transformed fluxes in f,g,h again
        f(:,i,j,k) = fTilde(:)*Metrics_fTilde(1,i,j,k,iElem) + &
                     gTilde(:)*Metrics_fTilde(2,i,j,k,iElem) + &
                     hTilde(:)*Metrics_fTilde(3,i,j,k,iElem)
        g(:,i,j,k) = fTilde(:)*Metrics_gTilde(1,i,j,k,iElem) + &
                     gTilde(:)*Metrics_gTilde(2,i,j,k,iElem) + &
                     hTilde(:)*Metrics_gTilde(3,i,j,k,iElem)
        h(:,i,j,k) = fTilde(:)*Metrics_hTilde(1,i,j,k,iElem) + &
                     gTilde(:)*Metrics_hTilde(2,i,j,k,iElem) + &
                     hTilde(:)*Metrics_hTilde(3,i,j,k,iElem)
      END DO ! i
    END DO ! j
  END DO ! k


  DO l=0,PP_N
    DO k=0,PP_N
      DO j=0,PP_N
        DO i=0,PP_N
          ! Update the time derivative with the spatial derivatives of the transformed fluxes
          Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + D_hat(i,l)*f(:,l,j,k) + &
                                                  D_hat(j,l)*g(:,i,l,k) + &
                                                  D_hat(k,l)*h(:,i,j,l)
        END DO !i
      END DO ! j
    END DO ! k
  END DO ! l
END DO ! iElem
END SUBROUTINE VolIntQDS


SUBROUTINE SurfIntQDS(Flux_Master,Flux_Slave,Ut,doMPISides)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
#if (PP_NodeType>1)
USE MOD_DG_Vars,            ONLY: L_HatPlus,L_HatMinus
#endif
USE MOD_Mesh_Vars,          ONLY: SideToElem
USE MOD_Mesh_Vars,          ONLY: nSides
USE MOD_Mesh_Vars,          ONLY: firstMPISide_YOUR,lastMPISide_MINE
USE MOD_QDS_Vars,           ONLY:QDSnVar,nQDSElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides+InnerSides+MPISides MINE  
REAL,INTENT(IN)    :: Flux_Master(1:QDSnVar,0:PP_N,0:PP_N,nSides)
REAL,INTENT(IN)    :: Flux_Slave(1:QDSnVar,0:PP_N,0:PP_N,nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)   :: Ut(QDSnVar,0:PP_N,0:PP_N,0:PP_N,1:nQDSElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: ElemID,Flip,SideID,locSideID
INTEGER            :: firstSideID,lastSideID
#if (PP_NodeType>1)
REAL               ::L_HatMinus0,L_HatPlusN 
#endif
!===================================================================================================================================

IF(doMPISides)THEN
  ! MPI YOUR
  firstSideID = firstMPISide_YOUR
   lastSideID = nSides
ELSE
  ! inner sides and MPI mine
  firstSideID = 1
   lastSideID = lastMPISide_MINE
END IF

#if (PP_NodeType>1)
L_HatMinus0 = L_HatMinus(0)
L_HatPlusN  = L_HatPlus(PP_N)
#endif
DO SideID=firstSideID,lastSideID
  ! neighbor side
  ElemID    = SideToElem(S2E_NB_ELEM_ID,SideID)
  locSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  flip      = SideToElem(S2E_FLIP,SideID)
  ! ignore MPI-faces and boundary faces
  IF(ElemID.LT.0) CYCLE ! boundary side is BC or MPI side
  CALL CalcSurfInt2(Flux_Slave(1:QDSnVar,0:PP_N,0:PP_N,SideID),Ut,Flip,ElemID,locSideID)
END DO ! SideID=1,nSides


DO SideID=firstSideID,lastSideID
  ! master side, flip=0
  ElemID    = SideToElem(S2E_ELEM_ID,SideID)  
  locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
  flip      = 0
  IF(ElemID.LT.0) CYCLE ! if master is MPI side
  CALL CalcSurfInt2(Flux_Master(1:QDSnVar,0:PP_N,0:PP_N,SideID),Ut,Flip,ElemID,locSideID)
END DO ! SideID=1,nSides

END SUBROUTINE SurfIntQDS


SUBROUTINE CalcSurfInt2(Flux,Ut,flip,ElemID,locSideID)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,            ONLY: L_HatPlus,L_HatMinus
USE MOD_QDS_Vars,           ONLY: QDSnVar,nQDSElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: Flux(1:QDSnVar,0:PP_N,0:PP_N)
INTEGER,INTENT(IN) :: flip,ElemID,locSideID!,SideID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)   :: Ut(QDSnVar,0:PP_N,0:PP_N,0:PP_N,1:nQDSElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q,l
#if (PP_NodeType>1)
REAL            ::L_HatMinus0,L_HatPlusN 
#endif
!===================================================================================================================================
#if (PP_NodeType>1)
L_HatMinus0 = L_HatMinus(0)
L_HatPlusN  = L_HatPlus(PP_N)
#endif

#if (PP_NodeType==1)
  SELECT CASE(locSideID)
!===================================================================================================================================
  CASE(XI_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut( :,l,p,q,ElemID)=Ut(:,l,p,q,ElemID)+Flux(:,q,p)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut( :,l,p,q,ElemID)=Ut(:,l,p,q,ElemID)-Flux(:,p,q)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut( :,l,p,q,ElemID)=Ut(:,l,p,q,ElemID)-Flux(:,PP_N-q,p)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut( :,l,p,q,ElemID)=Ut(:,l,p,q,ElemID)-Flux(:,PP_N-p,PP_N-q)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut( :,l,p,q,ElemID)=Ut(:,l,p,q,ElemID)-Flux(:,q,PP_N-p)*L_hatMinus(l)
      END DO; END DO; END DO ! l,p,q
    END SELECT
!===================================================================================================================================
  CASE(ETA_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut( :,p,l,q,ElemID)=Ut(:,p,l,q,ElemID)+Flux(:,p,q)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut( :,p,l,q,ElemID)=Ut(:,p,l,q,ElemID)-Flux(:,q,p)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut( :,p,l,q,ElemID)=Ut(:,p,l,q,ElemID)-Flux(:,PP_N-p,q)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut( :,p,l,q,ElemID)=Ut(:,p,l,q,ElemID)-Flux(:,PP_N-q,PP_N-p)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut( :,p,l,q,ElemID)=Ut(:,p,l,q,ElemID)-Flux(:,p,PP_N-q)*L_hatMinus(l)
      END DO; END DO; END DO ! p,l,q
    END SELECT
!===================================================================================================================================
  CASE(ZETA_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)= Ut(:,p,q,l,ElemID)+Flux(:,q,p)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(1) ! slave side, SideID=q,jSide=p
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)= Ut(:,p,q,l,ElemID)-Flux(:,p,q)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)= Ut(:,p,q,l,ElemID)-Flux(:,PP_N-q,p)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)= Ut(:,p,q,l,ElemID)-Flux(:,PP_N-p,PP_N-q)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)= Ut(:,p,q,l,ElemID)-Flux(:,q,PP_N-p)*L_hatMinus(l)
      END DO; END DO; END DO ! p,q,l
    END SELECT
!===================================================================================================================================
  CASE(XI_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut(:,l,p,q,ElemID)=Ut(:,l,p,q,ElemID)+Flux(:,p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut(:,l,p,q,ElemID)=Ut(:,l,p,q,ElemID)-Flux(:,q,p)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut(:,l,p,q,ElemID)=Ut(:,l,p,q,ElemID)-Flux(:,PP_N-p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut(:,l,p,q,ElemID)=Ut(:,l,p,q,ElemID)-Flux(:,PP_N-q,PP_N-p)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,PP_N; DO p=0,PP_N; DO l=0,PP_N
        Ut(:,l,p,q,ElemID)=Ut(:,l,p,q,ElemID)-Flux(:,p,PP_N-q)*L_hatPlus(l)
      END DO; END DO; END DO ! l,p,q
    END SELECT
!===================================================================================================================================
  CASE(ETA_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut(:,p,l,q,ElemID)=Ut(:,p,l,q,ElemID)+Flux(:,PP_N-p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(1) ! slave side, SideID=q,jSide=p
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut(:,p,l,q,ElemID)=Ut(:,p,l,q,ElemID)-Flux(:,q,PP_N-p)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut(:,p,l,q,ElemID)=Ut(:,p,l,q,ElemID)-Flux(:,p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut(:,p,l,q,ElemID)=Ut(:,p,l,q,ElemID)-Flux(:,PP_N-q,p)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO q=0,PP_N; DO l=0,PP_N; DO p=0,PP_N
        Ut(:,p,l,q,ElemID)=Ut(:,p,l,q,ElemID)-Flux(:,PP_N-p,PP_N-q)*L_hatPlus(l)
      END DO; END DO; END DO ! p,l,q
    END SELECT
!===================================================================================================================================
  CASE(ZETA_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0) ! master side
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)=Ut(:,p,q,l,ElemID)+Flux(:,p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(1) ! slave side, SideID=q,jSide=p
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)=Ut(:,p,q,l,ElemID)-Flux(:,q,p)*L_hatPlus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(2) ! slave side, SideID=N-p,jSide=q
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)=Ut(:,p,q,l,ElemID)-Flux(:,PP_N-p,q)*L_hatPlus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(3) ! slave side, SideID=N-q,jSide=N-p
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)=Ut(:,p,q,l,ElemID)-Flux(:,PP_N-q,PP_N-p)*L_hatPlus(l)
      END DO; END DO; END DO ! p,q,l
    CASE(4) ! slave side, SideID=p,jSide=N-q
      DO l=0,PP_N; DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,l,ElemID)=Ut(:,p,q,l,ElemID)-Flux(:,p,PP_N-q)*L_hatPlus(l)
      END DO; END DO; END DO ! p,q,l
    END SELECT
  END SELECT !locSideID
#else
  !update local grid cell
  SELECT CASE(locSideID)
!===================================================================================================================================
  CASE(XI_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,0,p,q,ElemID)=Ut(:,0,p,q,ElemID)+Flux(:,q,p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,0,p,q,ElemID)=Ut(:,0,p,q,ElemID)-Flux(:,p,q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,0,p,q,ElemID)=Ut(:,0,p,q,ElemID)-Flux(:,PP_N-q,p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,0,p,q,ElemID)=Ut(:,0,p,q,ElemID)-Flux(:,PP_N-p,PP_N-q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,0,p,q,ElemID)=Ut(:,0,p,q,ElemID)-Flux(:,q,PP_N-p)*L_hatMinus0
      END DO; END DO ! p,q
    END SELECT
  
  ! switch to right hand system for ETA_PLUS direction
!===================================================================================================================================
  CASE(ETA_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,0,q,ElemID)=Ut(:,p,0,q,ElemID)+Flux(:,p,q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,0,q,ElemID)=Ut(:,p,0,q,ElemID)-Flux(:,q,p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,0,q,ElemID)=Ut(:,p,0,q,ElemID)-Flux(:,PP_N-p,q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,0,q,ElemID)=Ut(:,p,0,q,ElemID)-Flux(:,PP_N-q,PP_N-p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,0,q,ElemID)=Ut(:,p,0,q,ElemID)-Flux(:,p,PP_N-q)*L_hatMinus0
      END DO; END DO ! p,q
    END SELECT
  
  ! switch to right hand system for ZETA_MINUS direction
!===================================================================================================================================
  CASE(ZETA_MINUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,0,ElemID)=Ut(:,p,q,0,ElemID)+Flux(:,q,p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,0,ElemID)=Ut(:,p,q,0,ElemID)-Flux(:,p,q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,0,ElemID)=Ut(:,p,q,0,ElemID)-Flux(:,PP_N-q,p)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,0,ElemID)=Ut(:,p,q,0,ElemID)-Flux(:,PP_N-p,PP_N-q)*L_hatMinus0
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,0,ElemID)=Ut(:,p,q,0,ElemID)-Flux(:,q,PP_N-p)*L_hatMinus0
      END DO; END DO ! p,q
    END SELECT
  
!===================================================================================================================================
  CASE(XI_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,PP_N,p,q,ElemID)=Ut(:,PP_N,p,q,ElemID)+Flux(:,p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,PP_N,p,q,ElemID)=Ut(:,PP_N,p,q,ElemID)-Flux(:,q,p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,PP_N,p,q,ElemID)=Ut(:,PP_N,p,q,ElemID)-Flux(:,PP_N-p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,PP_N,p,q,ElemID)=Ut(:,PP_N,p,q,ElemID)-Flux(:,PP_N-q,PP_N-p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,PP_N,p,q,ElemID)=Ut(:,PP_N,p,q,ElemID)-Flux(:,p,PP_N-q)*L_hatPlusN
      END DO; END DO ! p,q
    END SELECT
  
  ! switch to right hand system for ETA_PLUS direction
!===================================================================================================================================
  CASE(ETA_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,PP_N,q,ElemID)=Ut(:,p,PP_N,q,ElemID)+Flux(:,PP_N-p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,PP_N,q,ElemID)=Ut(:,p,PP_N,q,ElemID)-Flux(:,q,PP_N-p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,PP_N,q,ElemID)=Ut(:,p,PP_N,q,ElemID)-Flux(:,p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,PP_N,q,ElemID)=Ut(:,p,PP_N,q,ElemID)-Flux(:,PP_N-q,p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,PP_N,q,ElemID)=Ut(:,p,PP_N,q,ElemID)-Flux(:,PP_N-p,PP_N-q)*L_hatPlusN
      END DO; END DO ! p,q
    END SELECT

  ! switch to right hand system for ZETA_MINUS direction
!===================================================================================================================================
  CASE(ZETA_PLUS)
!===================================================================================================================================
    SELECT CASE(flip)
    CASE(0)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,PP_N,ElemID)=Ut(:,p,q,PP_N,ElemID)+Flux(:,p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(1)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,PP_N,ElemID)=Ut(:,p,q,PP_N,ElemID)-Flux(:,q,p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(2)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,PP_N,ElemID)=Ut(:,p,q,PP_N,ElemID)-Flux(:,PP_N-p,q)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(3)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,PP_N,ElemID)=Ut(:,p,q,PP_N,ElemID)-Flux(:,PP_N-q,PP_N-p)*L_hatPlusN
      END DO; END DO ! p,q
    CASE(4)
      DO q=0,PP_N; DO p=0,PP_N
        Ut(:,p,q,PP_N,ElemID)=Ut(:,p,q,PP_N,ElemID)-Flux(:,p,PP_N-q)*L_hatPlusN
      END DO; END DO ! p,q
    END SELECT
  END SELECT !locSideID
#endif
END SUBROUTINE CalcSurfInt2


SUBROUTINE GetBoundaryFluxQDS(t,tDeriv, Flux, U_Minus, NormVec, TangVec1, TangVec2, BCFace_xGP)
!===================================================================================================================================
! Computes the boundary values for a given Cartesian mesh face (defined by FaceID)
! BCType: 1...periodic, 2...exact BC
! Attention 1: this is only a tensor of local values U_Face and has to be stored into the right U_Left or U_Right in
!              SUBROUTINE CalcSurfInt
! Attention 2: U_FacePeriodic is only needed in the case of periodic boundary conditions
!===================================================================================================================================
! MODULES
USE MOD_Globals,        ONLY:Abort,CROSS
USE MOD_PreProc
USE MOD_Equation,       ONLY:ExactFunc
USE MOD_Equation_vars,  ONLY:c,c_inv
USE MOD_Mesh_Vars    ,  ONLY:nBCSides,nBCs,BoundaryType
USE MOD_Equation_Vars,  ONLY:nBCByType,BCSideID
USE MOD_Equation_Vars,  ONLY:BCData,nBCByType
USE MOD_QDS_Vars,       ONLY:QDSnVar
!USE MOD_Equation_Vars,  ONLY:IniExactFunc! richtig with particles???
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                    :: t
INTEGER,INTENT(IN)                   :: tDeriv
REAL,INTENT(IN)                      :: U_Minus(     QDSnVar,0:PP_N,0:PP_N,1:nBCSides)
REAL,INTENT(IN)                      :: NormVec(           3,0:PP_N,0:PP_N,1:nBCSides)
REAL,INTENT(IN),OPTIONAL             :: TangVec1(          3,0:PP_N,0:PP_N,1:nBCSides)
REAL,INTENT(IN),OPTIONAL             :: TangVec2(          3,0:PP_N,0:PP_N,1:nBCSides)
REAL,INTENT(IN)                      :: BCFace_xGP(        3,0:PP_N,0:PP_N,1:nBCSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: Flux( QDSnVar,0:PP_N,0:PP_N,1:nBCSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iBC,iSide,p,q,SideID
INTEGER                              :: BCType,BCState,nBCLoc
REAL                                 :: n_loc(3),resul(QDSnVar),epsBC
REAL                                 :: U_Face_loc(QDSnVar,0:PP_N,0:PP_N),v_nor(3), V_par(3), v_nor_abs
INTEGER                              :: iPart, iPart1, iPart2, iPart3, locPartNum
REAL                                 :: PartDens, Velo(3), Temp, Mass
!===================================================================================================================================

DO iBC=1,nBCs
  IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType(iBC,BC_TYPE)
  BCState=BoundaryType(iBC,BC_STATE)
  nBCLoc =nBCByType(iBC)
  SELECT CASE(BCType)
  CASE(1) !Periodic already filled!

  CASE(2) ! exact BC = Dirichlet BC !!
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          CALL ExactFunc(BCState,t,tDeriv,BCFace_xGP(:,p,q,SideID),U_Face_loc(:,p,q))
        END DO ! p
      END DO ! q
      ! Dirichlet means that we use the gradients from inside the grid cell
      CALL RiemannQDS(Flux(1:QDSnVar,:,:,SideID),U_Minus(:,:,:,SideID),U_Face_loc(  :,:,:), NormVec(:,:,:,SideID))
   END DO

  CASE(3) ! 1st order absorbing BC 
    U_Face_loc=0.
    CALL RiemannQDS(Flux(1:QDSnVar,:,:,SideID),U_Minus(:,:,:,SideID),U_Face_loc(:,:,:),NormVec(:,:,:,SideID))
    Flux = 0.0
  
  CASE(4) ! perfectly conducting surface (MunzOmnesSchneider 2000, pp. 97-98)
!    ! Determine the exact BC state
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      ! Determine the exact BC state
      DO q=0,PP_N
        DO p=0,PP_N
          DO iPart = 0, 7
            v_nor_abs = NormVec(1,p,q,SideID)*U_Minus(2+5*iPart,p,q,SideID) + &
                        NormVec(2,p,q,SideID)*U_Minus(3+5*iPart,p,q,SideID) + &
                        NormVec(3,p,q,SideID)*U_Minus(4+5*iPart,p,q,SideID)
            v_nor(1:3) = v_nor_abs*NormVec(1:3,p,q,SideID)
            v_par(1) = U_Minus(2+5*iPart,p,q,SideID) - v_nor(1)
            v_par(2) = U_Minus(3+5*iPart,p,q,SideID) - v_nor(2)
            v_par(3) = U_Minus(4+5*iPart,p,q,SideID) - v_nor(3)
            U_Face_loc(1+5*iPart,p,q) = U_Minus(1+5*iPart,p,q,SideID)
            U_Face_loc(2+5*iPart,p,q) = v_par(1) - v_nor(1)
            U_Face_loc(3+5*iPart,p,q) = v_par(2) - v_nor(2)
            U_Face_loc(4+5*iPart,p,q) = v_par(3) - v_nor(3)
            U_Face_loc(5+5*iPart,p,q) = U_Minus(5+5*iPart,p,q,SideID)
          END DO
        END DO ! p
      END DO ! q
      ! Dirichlet means that we use the gradients from inside the grid cell
      CALL RiemannQDS(Flux(1:QDSnVar,:,:,SideID),U_Minus(:,:,:,SideID),U_Face_loc(:,:,:),NormVec(:,:,:,SideID))
      Flux = 0.0
    END DO

  CASE(10) ! symmetry BC (perfect MAGNETIC conductor, PMC)

  CASE(20) ! exact BC = Dirichlet BC !!
!    ! SPECIAL BC: BCState uses readin state
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
!      ! Dirichlet means that we use the gradients from inside the grid cell
!      CALL RiemannQDS(Flux(1:QDSnVar,:,:,SideID),U_Minus( :,:,:,SideID),BCData(:,:,:,SideID),NormVec(:,:,:,SideID))
!    END DO

      DO q=0,PP_N; DO p=0,PP_N
        locPartNum = 0    
        DO iPart1=1,2; DO iPart2=1,2; DO iPart3=1,2
          Temp = NormVec(1,p,q,SideID)*U_Minus(2+locPartNum*5,p,q,SideID) &
                +NormVec(2,p,q,SideID)*U_Minus(3+locPartNum*5,p,q,SideID) &
                +NormVec(3,p,q,SideID)*U_Minus(4+locPartNum*5,p,q,SideID)
          Flux(1+locPartNum*5,p,q,SideID) =  U_Minus(1+locPartNum*5,p,q,SideID)*Temp
          Flux(2+locPartNum*5,p,q,SideID) =  U_Minus(2+locPartNum*5,p,q,SideID)*Temp
          Flux(3+locPartNum*5,p,q,SideID) =  U_Minus(3+locPartNum*5,p,q,SideID)*Temp
          Flux(4+locPartNum*5,p,q,SideID) =  U_Minus(4+locPartNum*5,p,q,SideID)*Temp
          Flux(5+locPartNum*5,p,q,SideID) =  U_Minus(5+locPartNum*5,p,q,SideID)*Temp
    
    
    !      U_Face_loc(1+locPartNum*5,p,q) =  U_Minus(1+locPartNum*5,p,q,SideID)
    !      U_Face_loc(2+locPartNum*5,p,q) =  U_Minus(2+locPartNum*5,p,q,SideID)
    !      U_Face_loc(3+locPartNum*5,p,q) =  U_Minus(3+locPartNum*5,p,q,SideID)
    !      U_Face_loc(4+locPartNum*5,p,q) =  U_Minus(4+locPartNum*5,p,q,SideID)
    !      U_Face_loc(5+locPartNum*5,p,q) =  U_Minus(5+locPartNum*5,p,q,SideID)
          locPartNum = locPartNum + 1
        END DO; END DO; END DO
      END DO; END DO
    END DO

  
  CASE DEFAULT ! unknown BCType
    CALL abort(&
__STAMP__&
        ,'no BC defined in maxwell/getboundaryflux.f90!')
  END SELECT ! BCType
END DO ! iBC=1,nBC

IF(1.EQ.2)THEN
  epsBC=TangVec1(1,1,1,1)
  epsBC=TangVec2(1,1,1,1)
END IF

END SUBROUTINE GetBoundaryFluxQDS


SUBROUTINE RiemannQDS(F,U_L,U_R,nv)
!===================================================================================================================================
! Computes the numerical flux
! Conservative States are rotated into normal direction in this routine and are NOT backrotatet: don't use it after this routine!!
!===================================================================================================================================
! MODULES
USE MOD_PreProc ! PP_N
USE MOD_QDS_Vars,     ONLY:QDSnVar,QDSMaxVelo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(QDSnVar,0:PP_N,0:PP_N),INTENT(IN) :: U_L,U_R
REAL,INTENT(IN)                                  :: nv(3,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                                 :: F(QDSnVar,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                                             :: n_loc(3),A_p(4,4),A_n(4,4)
INTEGER                                          :: Count_1,Count_2, iVar
!REAL                                             :: D(3,3)                  ! auxiliary matrices used 
!REAL                                             :: E(3,3), E_trans(3,3)    ! auxiliary matrices used
REAL                                            :: Lambda_L, Lambda_R, velocompL, velocompR,LambdaMax, LambdaAbs
!===================================================================================================================================
!Lax-Friedrich
DO iVar=0,7
  DO Count_2=0,PP_N; DO Count_1=0,PP_N 
    velocompL = U_L(2+iVar*5,Count_1,Count_2)*nv(1,Count_1,Count_2) + U_L(3+iVar*5,Count_1,Count_2)*nv(2,Count_1,Count_2) &
            + U_L(4+iVar*5,Count_1,Count_2)*nv(3,Count_1,Count_2)
    velocompR = U_R(2+iVar*5,Count_1,Count_2)*nv(1,Count_1,Count_2) + U_R(3+iVar*5,Count_1,Count_2)*nv(2,Count_1,Count_2) &
            + U_R(4+iVar*5,Count_1,Count_2)*nv(3,Count_1,Count_2)
    !IF (ABS(velocompL).GT.ABS(velocompR)) THEN
      !LambdaMax = ABS(velocompL)
    !ELSE
      !LambdaMax = ABS(velocompR)
    !END IF
!    LambdaMax = MERGE(velocompL, velocompR, ABS(velocompL).GT.ABS(velocompR))
    LambdaMax=QDSMaxVelo

!    Lambda_L = 0.5 * (LambdaMax + ABS(LambdaMax))
!    Lambda_R = 0.5 * (LambdaMax - ABS(LambdaMax))
!    F(1 + iVar*5,Count_1,Count_2) =  (Lambda_L * U_L(1 + iVar*5,Count_1,Count_2) + Lambda_R * U_R(1 + iVar*5,Count_1,Count_2)) 
!    F(2 + iVar*5,Count_1,Count_2) =  (Lambda_L * U_L(2 + iVar*5,Count_1,Count_2) + Lambda_R * U_R(2 + iVar*5,Count_1,Count_2))
!    F(3 + iVar*5,Count_1,Count_2) =  (Lambda_L * U_L(3 + iVar*5,Count_1,Count_2) + Lambda_R * U_R(3 + iVar*5,Count_1,Count_2))
!    F(4 + iVar*5,Count_1,Count_2) =  (Lambda_L * U_L(4 + iVar*5,Count_1,Count_2) + Lambda_R * U_R(4 + iVar*5,Count_1,Count_2))
!    F(5 + iVar*5,Count_1,Count_2) =  (Lambda_L * U_L(5 + iVar*5,Count_1,Count_2) + Lambda_R * U_R(5 + iVar*5,Count_1,Count_2)) 
    


     F(1 + iVar*5,Count_1,Count_2) = 0.5*(velocompL* U_L(1 + iVar*5,Count_1,Count_2) + velocompR* U_R(1 + iVar*5,Count_1,Count_2)) &
          - 0.5* LambdaMax * (U_R(1 + iVar*5,Count_1,Count_2)- U_L(1 + iVar*5,Count_1,Count_2))
     F(2 + iVar*5,Count_1,Count_2) = 0.5*(velocompL* U_L(2 + iVar*5,Count_1,Count_2) + velocompR* U_R(2 + iVar*5,Count_1,Count_2)) &
          - 0.5* LambdaMax * (U_R(2 + iVar*5,Count_1,Count_2)- U_L(2 + iVar*5,Count_1,Count_2))
     F(3 + iVar*5,Count_1,Count_2) = 0.5*(velocompL* U_L(3 + iVar*5,Count_1,Count_2) + velocompR* U_R(3 + iVar*5,Count_1,Count_2)) &
          - 0.5* LambdaMax * (U_R(3 + iVar*5,Count_1,Count_2)- U_L(3 + iVar*5,Count_1,Count_2))
     F(4 + iVar*5,Count_1,Count_2) = 0.5*(velocompL* U_L(4 + iVar*5,Count_1,Count_2) + velocompR* U_R(4 + iVar*5,Count_1,Count_2)) &
          - 0.5* LambdaMax * (U_R(4 + iVar*5,Count_1,Count_2)- U_L(4 + iVar*5,Count_1,Count_2))
     F(5 + iVar*5,Count_1,Count_2) = 0.5*(velocompL* U_L(5 + iVar*5,Count_1,Count_2) + velocompR* U_R(5 + iVar*5,Count_1,Count_2)) &
          - 0.5* LambdaMax * (U_R(5 + iVar*5,Count_1,Count_2)- U_L(5 + iVar*5,Count_1,Count_2))

!     F(1 + iVar*5,Count_1,Count_2) = 0.5*(velocompL* U_L(1 + iVar*5,Count_1,Count_2) + velocompR* U_R(1 + iVar*5,Count_1,Count_2))
!     F(2 + iVar*5,Count_1,Count_2) = 0.5*(velocompL* U_L(2 + iVar*5,Count_1,Count_2) + velocompR* U_R(2 + iVar*5,Count_1,Count_2))
!     F(3 + iVar*5,Count_1,Count_2) = 0.5*(velocompL* U_L(3 + iVar*5,Count_1,Count_2) + velocompR* U_R(3 + iVar*5,Count_1,Count_2))
!     F(4 + iVar*5,Count_1,Count_2) = 0.5*(velocompL* U_L(4 + iVar*5,Count_1,Count_2) + velocompR* U_R(4 + iVar*5,Count_1,Count_2))
!     F(5 + iVar*5,Count_1,Count_2) = 0.5*(velocompL* U_L(5 + iVar*5,Count_1,Count_2) + velocompR* U_R(5 + iVar*5,Count_1,Count_2))

     
!    LambdaMax = MAX(ABS(velocompL), ABS(velocompR))
!    F(1 + iVar*5,Count_1,Count_2) =  0.5*(velocompL * U_L(1 + iVar*5,Count_1,Count_2) +velocompR * U_R(1 + iVar*5,Count_1,Count_2) &
!          + LambdaMax *(U_L(1 + iVar*5,Count_1,Count_2) -  U_R(1 + iVar*5,Count_1,Count_2)))
!    F(2 + iVar*5,Count_1,Count_2) =  0.5*(velocompL * U_L(2 + iVar*5,Count_1,Count_2) +velocompR * U_R(2 + iVar*5,Count_1,Count_2) &
!          + LambdaMax *(U_L(2 + iVar*5,Count_1,Count_2) -  U_R(2 + iVar*5,Count_1,Count_2)))
!    F(3 + iVar*5,Count_1,Count_2) =  0.5*(velocompL * U_L(3 + iVar*5,Count_1,Count_2) +velocompR * U_R(3 + iVar*5,Count_1,Count_2) &
!          + LambdaMax *(U_L(3 + iVar*5,Count_1,Count_2) -  U_R(3 + iVar*5,Count_1,Count_2)))
!    F(4 + iVar*5,Count_1,Count_2) =  0.5*(velocompL * U_L(4 + iVar*5,Count_1,Count_2) +velocompR * U_R(4 + iVar*5,Count_1,Count_2) &
!          + LambdaMax *(U_L(4 + iVar*5,Count_1,Count_2) -  U_R(4 + iVar*5,Count_1,Count_2)))
!    F(5 + iVar*5,Count_1,Count_2) =  0.5*(velocompL * U_L(5 + iVar*5,Count_1,Count_2) +velocompR * U_R(5 + iVar*5,Count_1,Count_2) &
!          + LambdaMax *(U_L(5 + iVar*5,Count_1,Count_2) -  U_R(5 + iVar*5,Count_1,Count_2)))
  END DO; END DO
END DO


END SUBROUTINE RiemannQDS


SUBROUTINE EvalFlux3DQDS(iElem,f,g,h)
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
USE MOD_PreProc ! PP_N
USE MOD_QDS_Vars, ONLY:QDSnVar,UQDS
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                 :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(QDSnVar,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: f,g,h    ! Cartesian fluxes (iVar,i,j,k)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: Uin(QDSnVar)
INTEGER             :: i,j,k,iVar
!===================================================================================================================================
DO k=0,PP_N
  DO j=0,PP_N
    DO i=0,PP_N
      Uin=UQDS(:,i,j,k,iElem)
      DO iVar=0,7
        ! hier der physikalische Fluss ohne die Divergenzkorrektur!
        !A
        f(1+5*iVar,i,j,k) = UQDS(2+5*iVar,i,j,k,iElem)*UQDS(1+5*iVar,i,j,k,iElem) 
        f(2+5*iVar,i,j,k) = UQDS(2+5*iVar,i,j,k,iElem)*UQDS(2+5*iVar,i,j,k,iElem) 
        f(3+5*iVar,i,j,k) = UQDS(2+5*iVar,i,j,k,iElem)*UQDS(3+5*iVar,i,j,k,iElem) 
        f(4+5*iVar,i,j,k) = UQDS(2+5*iVar,i,j,k,iElem)*UQDS(4+5*iVar,i,j,k,iElem) 
        f(5+5*iVar,i,j,k) = UQDS(2+5*iVar,i,j,k,iElem)*UQDS(5+5*iVar,i,j,k,iElem) 
        !B
        g(1+5*iVar,i,j,k) = UQDS(3+5*iVar,i,j,k,iElem)*UQDS(1+5*iVar,i,j,k,iElem) 
        g(2+5*iVar,i,j,k) = UQDS(3+5*iVar,i,j,k,iElem)*UQDS(2+5*iVar,i,j,k,iElem) 
        g(3+5*iVar,i,j,k) = UQDS(3+5*iVar,i,j,k,iElem)*UQDS(3+5*iVar,i,j,k,iElem) 
        g(4+5*iVar,i,j,k) = UQDS(3+5*iVar,i,j,k,iElem)*UQDS(4+5*iVar,i,j,k,iElem) 
        g(5+5*iVar,i,j,k) = UQDS(3+5*iVar,i,j,k,iElem)*UQDS(5+5*iVar,i,j,k,iElem) 
        !C
        h(1+5*iVar,i,j,k) = UQDS(4+5*iVar,i,j,k,iElem)*UQDS(1+5*iVar,i,j,k,iElem) 
        h(2+5*iVar,i,j,k) = UQDS(4+5*iVar,i,j,k,iElem)*UQDS(2+5*iVar,i,j,k,iElem) 
        h(3+5*iVar,i,j,k) = UQDS(4+5*iVar,i,j,k,iElem)*UQDS(3+5*iVar,i,j,k,iElem) 
        h(4+5*iVar,i,j,k) = UQDS(4+5*iVar,i,j,k,iElem)*UQDS(4+5*iVar,i,j,k,iElem) 
        h(5+5*iVar,i,j,k) = UQDS(4+5*iVar,i,j,k,iElem)*UQDS(5+5*iVar,i,j,k,iElem) 
      END DO
    END DO ! i
  END DO ! j
END DO ! k
END SUBROUTINE EvalFlux3DQDS

SUBROUTINE QDSReCalculateDGValues()
!===================================================================================================================================
! Get the constant advection velocity vector from the ini file 
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars,    ONLY:PI
USE MOD_QDS_Vars 
USE MOD_PARTICLE_Vars,      ONLY : BoltzmannConst
!USE MOD_QDS_Vars,         ONLY :  QDS_Species
USE MOD_PreProc
USE MOD_QDS_Vars,           ONLY : QDSSpeciesMass
USE MOD_Mesh_Vars,          ONLY : sJ
USE MOD_Interpolation_Vars, ONLY : wGP
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iElem, k, j, i, locPartNum, iPart1, iPart2, iPart3
REAL        :: checkSum
!===================================================================================================================================
!QDSSpeciesMass=Species(QDS_Species)%MassIC
! Read the maximum number of time steps MaxIter and the end time TEnd from ini file
DO iElem = 1, nQDSElems


IF(1.EQ.2)THEN
  IF (QDSMacroValues(1,0,0,0,iElem).NE.0.0) then
      print*, 'build parts', iElem
      print*, QDSMacroValues(1,0,0,0,iElem)/(QDSSpeciesMass*wGP(0)*wGP(0)*wGP(0)/sJ(0,0,0,iElem)), &
      QDSMacroValues(2:4,0,0,0,iElem)/QDSMacroValues(1,0,0,0,iElem), QDSMacroValues(6,0,0,0,iElem)
  end if
  IF (QDSMacroValues(1,1,1,1,iElem).NE.0.0) then
      print*, 'zwei', iElem
      print*, QDSMacroValues(1,1,1,1,iElem)/(QDSSpeciesMass*wGP(1)*wGP(1)*wGP(1)/sJ(1,1,1,iElem)), &
      QDSMacroValues(2:4,1,1,1,iElem)/QDSMacroValues(1,1,1,1,iElem), QDSMacroValues(6,1,1,1,iElem)
  end if
END IF


END DO
DO iElem = 1, nQDSElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        locPartNum = 0
        IF (QDSMacroValues(1,i,j,k,iElem).GT.0.0) THEN
          IF (QDSMacroValues(6,i,j,k,iElem).LT.0.0) QDSMacroValues(6,i,j,k,iElem) = 0.0
          DO iPart1=1,2; DO iPart2=1,2; DO iPart3=1,2
            UQDS(1+locPartNum*5,i,j,k,iElem) =  QDSMacroValues(1,i,j,k,iElem)*GaussHermitWeiAbs(1,iPart1) &
                *GaussHermitWeiAbs(1,iPart2)*GaussHermitWeiAbs(1,iPart3)/(PI*SQRT(PI))
            UQDS(2+locPartNum*5,i,j,k,iElem) =  QDSMacroValues(2,i,j,k,iElem) / QDSMacroValues(1,i,j,k,iElem) &
                  + SQRT(2.*BoltzmannConst*QDSMacroValues(6,i,j,k,iElem)/QDSSpeciesMass)*GaussHermitWeiAbs(2,iPart1)
            UQDS(3+locPartNum*5,i,j,k,iElem) =  QDSMacroValues(3,i,j,k,iElem) / QDSMacroValues(1,i,j,k,iElem) &
                  + SQRT(2.*BoltzmannConst*QDSMacroValues(6,i,j,k,iElem)/QDSSpeciesMass)*GaussHermitWeiAbs(2,iPart2)
            UQDS(4+locPartNum*5,i,j,k,iElem) =  QDSMacroValues(4,i,j,k,iElem) / QDSMacroValues(1,i,j,k,iElem) &
                  + SQRT(2.*BoltzmannConst*QDSMacroValues(6,i,j,k,iElem)/QDSSpeciesMass)*GaussHermitWeiAbs(2,iPart3)
            UQDS(5+locPartNum*5,i,j,k,iElem) =(QDSSpecDOF-3.)*BoltzmannConst*QDSMacroValues(6,i,j,k,iElem) &
                  /(QDSSpeciesMass*2.)
            locPartNum = locPartNum + 1
          END DO; END DO; END DO
        ELSE
          UQDS(:,i,j,k,iElem) = 0.0
        END IF
      END DO ! i
    END DO ! j
  END DO
END DO
!read*
END SUBROUTINE QDSReCalculateDGValues

SUBROUTINE QDSCalculateMacroValues()
!===================================================================================================================================
! Get the constant advection velocity vector from the ini file 
!===================================================================================================================================
! MODULES
USE MOD_QDS_Vars 
USE MOD_PARTICLE_Vars,      ONLY : BoltzmannConst
!USE MOD_QDS_Vars,           ONLY : QDS_Species
USE MOD_QDS_Vars,           ONLY : QDSSpeciesMass
USE MOD_PreProc

USE MOD_Mesh_Vars,          ONLY : sJ
USE MOD_Interpolation_Vars, ONLY : wGP
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iElem, k, j, i, iPart
REAL :: Temp, Velo(3), Dens, Mass
!===================================================================================================================================
!QDSSpeciesMass=Species(QDS_Species)%MassIC
! Read the maximum number of time steps MaxIter and the end time TEnd from ini file
DO iElem = 1, nQDSElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        QDSMacroValues(:,i,j,k,iElem) = 0.0
        DO iPart=0,7
          !IF (i.eq.1.and.j.eq.1.and.k.eq.1.and.ielem.eq.1) print*, UQDS(1+iPart*5,i,j,k,iElem)
          QDSMacroValues(1,i,j,k,iElem) = QDSMacroValues(1,i,j,k,iElem) + UQDS(1+iPart*5,i,j,k,iElem)
          QDSMacroValues(2,i,j,k,iElem) = QDSMacroValues(2,i,j,k,iElem) + UQDS(2+iPart*5,i,j,k,iElem)*UQDS(1+iPart*5,i,j,k,iElem)
          QDSMacroValues(3,i,j,k,iElem) = QDSMacroValues(3,i,j,k,iElem) + UQDS(3+iPart*5,i,j,k,iElem)*UQDS(1+iPart*5,i,j,k,iElem)
          QDSMacroValues(4,i,j,k,iElem) = QDSMacroValues(4,i,j,k,iElem) + UQDS(4+iPart*5,i,j,k,iElem)*UQDS(1+iPart*5,i,j,k,iElem)
          QDSMacroValues(5,i,j,k,iElem) = QDSMacroValues(5,i,j,k,iElem) + UQDS(1+iPart*5,i,j,k,iElem) &
            *(0.5*(UQDS(2+iPart*5,i,j,k,iElem)**2.+&
                   UQDS(3+iPart*5,i,j,k,iElem)**2.+&
                   UQDS(4+iPart*5,i,j,k,iElem)**2.)+&
                   UQDS(5+iPart*5,i,j,k,iElem))
!          IF (k.eq.0.and.j.eq.0.and.i.eq.0) print*,1+iPart*5, UQDS(1+iPart*5,i,j,k,iElem),UQDS(2+iPart*5,i,j,k,iElem),& 
!            UQDS(3+iPart*5,i,j,k,iElem)  ,UQDS(4+iPart*5,i,j,k,iElem),UQDS(5+iPart*5,i,j,k,iElem)
        END DO        
        IF (QDSMacroValues(1,i,j,k,iElem).GT.0.0) THEN
          QDSMacroValues(6,i,j,k,iElem) = (2.*QDSMacroValues(5,i,j,k,iElem) &
            -(QDSMacroValues(2,i,j,k,iElem)**2.0+QDSMacroValues(3,i,j,k,iElem)**2.0+QDSMacroValues(4,i,j,k,iElem)**2.0) &
            /QDSMacroValues(1,i,j,k,iElem)) / (QDSMacroValues(1,i,j,k,iElem)*3.) *QDSSpeciesMass /BoltzmannConst
        ELSE
          QDSMacroValues(6,i,j,k,iElem) = 0.0
        END IF
      END DO ! i
    END DO ! j
  END DO
END DO

!  Temp = 278.687
!  Velo =(/0.,0.,0./) !(/1020.882,0.,0./)
!  Dens= 2.633459376E25
!  Mass = 4.651734101E-26
!DO k=0,PP_N
!  DO j=0,PP_N
!    DO i=0,PP_N
!  QDSMacroValues(1,i,j,k,1) = Dens*wGP(i)*wGP(j)*wGP(k)/sJ(i,j,k,1)*Mass
!  QDSMacroValues(2:4,i,j,k,1) = QDSMacroValues(1,i,j,k,1)*Velo(1:3)
!  QDSMacroValues(6,i,j,k,1) = Temp
!END DO; END DO; END DO
!print*,'energie and temp', QDSMacroValues(5,0,0,0,1),  QDSMacroValues(6,0,0,0,1)
!read*
END SUBROUTINE QDSCalculateMacroValues




SUBROUTINE FinalizeQDS()
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_QDS_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
QDSInitIsDone = .FALSE.
SDEALLOCATE(UQDS)
SDEALLOCATE(UQDSt)
SDEALLOCATE(UQDS_Master)
SDEALLOCATE(UQDS_Slave)
SDEALLOCATE(FluxQDS_Master)
SDEALLOCATE(FluxQDS_Slave)
END SUBROUTINE FinalizeQDS


END MODULE MOD_QDS
