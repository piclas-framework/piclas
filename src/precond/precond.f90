#include "boltzplatz.h"

MODULE MOD_Precond
!===================================================================================================================================
! Module for the Block-Jacobi Preconditioner  
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitPrecond
  MODULE PROCEDURE InitPrecond
END INTERFACE

#ifdef maxwell
INTERFACE BuildPrecond
  MODULE PROCEDURE BuildPrecond 
END INTERFACE
#endif /*maxwell*/

INTERFACE FinalizePrecond
  MODULE PROCEDURE FinalizePrecond
END INTERFACE

#ifdef maxwell
INTERFACE BuildnVecSurf
  MODULE PROCEDURE BuildnVecSurf
END INTERFACE
#endif /*maxwell*/

PUBLIC :: InitPrecond
#ifdef maxwell
PUBLIC :: BuildPrecond
PUBLIC :: BuildnVecSurf
#endif /*maxwell*/
PUBLIC :: FinalizePrecond
!===================================================================================================================================

CONTAINS

SUBROUTINE InitPrecond()
!===================================================================================================================================
! Initialize preconditioner and call initialize of type of preconditioner
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Precond_Vars
#ifdef maxwell
!USE MOD_CSR_Vars,          ONLY:EILU
USE MOD_LinearSolver_Vars, ONLY:nDOFelem
USE MOD_Jac_ex,            ONLY:InitJac_ex
USE MOD_Jac_fd,            ONLY:InitJac_Fd
USE MOD_JacDG,             ONLY:InitJacDG
USE MOD_ReadInTools,       ONLY:GETINT
USE MOD_SparseILU,         ONLY:InitSparseILU
USE MOD_ILU,               ONLY:InitILU
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: locSize
!===================================================================================================================================
IF(PrecondInitIsDone)THEN
   SWRITE(*,*) "InitPrecond already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PRECONDITIONER...'

#ifdef maxwell
PrecondType   = GETINT('PrecondType','0')
PrecondMethod = GETINT('PrecondMethod','0')
DebugMatrix   = GETINT('DebugMatrix','0')
IF(PrecondType.EQ.0) PrecondMethod=-1
!IF(PrecondType.EQ.9) PrecondMethod=3

IF(PrecondType.NE.0)THEN
  SWRITE(UNIT_stdOut,'(A)') ' BUILDING LOCAL nVec,Surf ...'
  ! Now the side based normal vector and surface element are rotated in the correct
  ! local system
  CALL BuildnVecSurf()
  SWRITE(UNIT_stdOut,'(A)') ' ... DONE.'
END IF

ALLOCATE( neighborElemID    (1:6,1:PP_nElems) )!&
        !, neighborlocSideID (1:6,1:PP_nElems) )
neighborElemID=-1
!neighborlocSideID=-1

SELECT CASE(PrecondType)
  CASE(0)
    SWRITE(UNIT_stdOut,'(A)') ' No Preconditioner'
  CASE(1)
    SWRITE(UNIT_stdOut,'(A)') ' BJ with FD-Matrix'
    ALLOCATE(invP(1:nDOFelem,1:nDOFelem,PP_nElems))
    CALL InitJac_FD()
    !CALL InitJac_Ex()
    !CALL InitJacDG()
  CASE(2)
    SWRITE(UNIT_stdOut,'(A)') ' BJ'
    ALLOCATE(invP(1:nDOFelem,1:nDOFelem,PP_nElems))
    !CALL InitJac_FD()
    CALL InitJac_Ex()
    !CALL InitJacDG()
  CASE(3)
    SWRITE(UNIT_stdOut,'(A)') ' Element-ILU(0)'
    ! sparse ILU0
    CALL InitJac_Ex()
    !CALL InitJacDG()
    CALL InitSparseILU()
  CASE(4)
    SWRITE(UNIT_stdOut,'(A)') ' Element-BILU(0)'
    CALL InitJac_Ex()
    !CALL InitJacDG()
    CALL InitILU()
  CASE DEFAULT
    CALL abort(__STAMP__,'Preconditioner not implemented!',999,999.)
END SELECT

#else
PrecondType=0
#endif /*maxwell*/

PrecondInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PRECONDITIONER DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitPrecond

#ifdef maxwell
SUBROUTINE BuildPrecond(t,tStage,tDeriv,alpha,dt)
!===================================================================================================================================
! Build preconditioner for each element, calls a type of preconditioner 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Basis             ,ONLY:GetInverse
USE MOD_LinearSolver_Vars ,ONLY:nDOFelem,nDOFLine,mass
USE MOD_Precond_Vars      ,ONLY:invP,PrecondType,DebugMatrix,invBJ,invJ, neighborElemID, ProcJacobian
USE MOD_Jac_ex            ,ONLY:Jac_ex, Jac_Ex_Neighbor,Jac_ex1D
USE MOD_Jac_FD            ,ONLY:Jac_FD_slow!,Jac_FD
USE MOD_Jac_FD_Vars,       ONLY:reps0,XK
USE MOD_JacDG             ,ONLY:BuildJacDG
USE MOD_Mesh_Vars         ,ONLY:nBCSides,ElemToSide,SideToElem
!USE MOD_Prepare_FD    ,ONLY:PrepareFD
USE MOD_Precond_Vars      ,ONLy:invXi,invEta,invZeta,dRdXi,dRdZeta,dRdEta
USE MOD_DG,                ONLY: DGTimeDerivative_WeakForm
USE MOD_DG_Vars,           ONLY: U,Ut
USE MOD_SparseILU,         ONLY: BuildILU0
USE MOD_ILU,               ONLY: BuildBILU0BCSR
USE MOD_CSR,               ONLY: CSR
USE MOD_Mesh_Vars,         ONLY: nInnerSides,nMPISides_MINE,nMPISides_YOUR,nSides
!USE MOD_CSR_Vars,      ONLY: DiagonalEntries
#ifdef MPI
USE MOD_MPI_Vars
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: alpha,dt,t,tStage
INTEGER,INTENT(IN) :: tDeriv
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: r,s,iElem,ilocSide,SideID,iElem2
REAL,DIMENSION(:,:),ALLOCATABLE :: Ploc,Ploc1!,Ploc2
REAL,DIMENSION(:,:,:,:,:),ALLOCATABLE :: BJDOF
!REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: dRdXi,dRdEta,dRdZeta
!REAL,DIMENSION(:,:,:,:),ALLOCATABLE   :: GlobalJacobian
!REAL,DIMENSION(:,:),ALLOCATABLE       :: GlobalJacobian2,GlobalJacobianTild
REAL               :: dummy
REAL               :: TimeStart(3)
REAL               :: TimeEnd(3)
REAL               :: TotalTime(3)
REAL               :: delta(0:PP_N,0:PP_N)
REAL               :: tmpMat(PP_nVar,PP_nVar),tmp2Mat(PP_nVar,PP_nVar)
CHARACTER(LEN=255) :: Filename
CHARACTER(LEN=17)  :: strfmt
INTEGER            :: i,j,k,vn1,vn2,mm,iVar,nn,oo,p,q,vni,vnj,vnk,vnmm,vnoo,vnnn,l,vnl
REAL               :: coeff
REAL               :: tmp(1:nDOFLine,1:nDOFLine)
INTEGER            :: nTotalDOF,NBElemID,nlocSides
!===================================================================================================================================

IF(PrecondType.EQ.0) RETURN !NO PRECONDITIONER
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' BUILD PRECONDITIONER...'

delta=0.
DO i=0,PP_N
  delta(i,i)=1.
END DO

SELECT CASE(PrecondType)
  CASE(1,2)
    ALLOCATE( Ploc (1:nDOFElem,1:nDOFElem), &
              Ploc1(1:nDOFElem,1:nDOFElem) )
  CASE(3,4) 
    ALLOCATE( Ploc (1:nDOFElem,1:nDOFElem))
  CASE(204,215)
    ALLOCATE( Ploc (1:nDOFElem,1:nDOFElem), &
              Ploc1(1:nDOFElem,1:nDOFElem))!&
              !Ploc2(1:nDOFElem,1:nDOFElem)  )
END SELECT

TotalTime=0.

! PrecondType .LT. 99 - element local
DO iElem=1,PP_nElems
  CALL CPU_TIME(TimeStart(1))
  !  WRITE(UNIT_stdOut,'(A,I6,A)')'Build Jacobian, iElem=',iElem,' ... '
  ! nullify
  DO s=1,nDOFelem
    DO r=1,nDOFelem
      Ploc(r,s)=0.
    END DO !r
  END DO !s
  IF(PrecondType.EQ.1) THEN 
    ! obtained by finite difference
    !Prepare Linearisation State 
    CALL DGTimeDerivative_WeakForm(t,tStage,tDeriv,doSource=.FALSE.)
    ! finit differences per Element ! never to use ... 
    CALL Jac_FD_slow(t,tStage,tDeriv,iElem,Ploc)
  ELSE
    ! analytic per Element
    CALL Jac_ex(iElem,Ploc)
  END IF 
  IF(DebugMatrix.NE.0) Ploc1=Ploc

  !IF(PrecondType.NE.60)THEN
  ! add contibution I-alpha*dt*dRdU
  coeff=-alpha*dt
  DO s=1,nDOFelem
    DO r=1,nDOFelem
      Ploc(r,s)=coeff*Ploc(r,s)
    END DO !r
    Ploc(s,s)=Ploc(s,s)+1.
  END DO !s
  DO s=0,nDOFelem-1,PP_nVar
    r=0
    DO k=0,PP_N
      DO j=0,PP_N
        DO i=0,PP_N
          Ploc(r+1:r+PP_nVar,s+1:s+PP_nVar) = Ploc(r+1:r+PP_nVar,s+1:s+PP_nVar)*Mass(1,i,j,k,iElem)
          r=r+PP_nVar
        END DO !i
      END DO !j
    END DO !k
  END DO ! s
  CALL CPU_TIME(TimeEnd(1))

  CALL CPU_TIME(TimeStart(2))
  SELECT CASE(PrecondType)
  CASE(1,2) ! element block jacobi
    !invert Ploc => invP(:,:,iElem)
    invP(:,:,iElem)=getInverse(nDOFelem,Ploc)
  CASE(3) ! ilu(0) of element
    CALL BuildILU0(Ploc,iElem)
  CASE(4)
    CALL BuildBILU0BCSR(Ploc,iElem)
  END SELECT
  CALL CPU_TIME(TimeEnd(2))
  TotalTime=TotalTime+(TimeEnd-TimeStart)
  IF((PrecondType.LE.2).AND.DebugMatrix.NE.0)THEN
    CALL CheckBJPrecond(Ploc1,Ploc,invP(:,:,iElem),iElem,TotalTime(3))
  END IF ! DebugMatrix
END DO ! ! iELEM

WRITE(UNIT_stdOut,'(A,F11.3,A)')' TOTAL DERIVATING TIME =[',TotalTime(1),' ]'
WRITE(UNIT_stdOut,'(A,F11.3,A)')' TOTAL INVERTING  TIME =[',TotalTime(2),' ]'
IF(debugMatrix.GT.2) WRITE(UNIT_stdOut,'(A,F11.3,A)')' TOTAL CHECKING  TIME =[',TotalTime(3),' ]'
SWRITE(UNIT_stdOut,'(A)')' BUILD PRECONDITIONER DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')


SELECT CASE(PrecondType)
  CASE(1,2) 
    DEALLOCATE( Ploc, Ploc1)
  CASE(3,4)
    DEALLOCATE( Ploc )
END SELECT

END SUBROUTINE  BuildPrecond


SUBROUTINE BuildnVecSurf()
!===================================================================================================================================
! used for BR2: normal vectors, outward pointing and sorted in ijk element fashion!!! 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars          ,ONLY: Normvec,SurfElem,ElemToSide,nSides
USE MOD_Precond_Vars       ,ONLY: nVec,Surf,BuildNVecisDone
USE MOD_Mesh_Vars          ,ONLY:SideID_plus_lower,SideID_plus_upper,SideID_minus_upper,SideID_minus_lower
#ifdef MPI
! exchange of normal vector and surface element
USE MOD_MPI_Vars
USE MOD_MPI,              ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER                                      :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                      :: p,q,iSide,SideID,Flip
REAL,ALLOCATABLE                             :: NormVecPlus(:,:,:,:)  
!===================================================================================================================================

IF(BuildNvecisDone) RETURN


! global
! this are the local normal vector and surface element for each interpolation point
ALLOCATE(nVec(1:3,0:PP_N,0:PP_N,1:6,PP_nElems))
ALLOCATE(Surf(0:PP_N,0:PP_N,1:6,PP_nElems))

! usefull only for Precondioner and local sides
ALLOCATE(       NormVecPlus(4,0:PP_N,0:PP_N,1:nSides)) 

!DO SideID=SideID_Plus_Lower,SideID_Minus_Upper
DO SideID=SideID_minus_Lower,SideID_Minus_Upper
  normVecPlus (1:3,:,:,SideID) = -normVec(1:3,:,:,SideID)
  normVecPlus ( 4 ,:,:,SideID) = SurfElem(:,:,SideID)
END DO

! like flux
! communicate sideID_minus_upper : sideID_plus_upper
! performed in buildnvecsurf because then the necessary MPI Comm is build

! communicate sideID_minus_upper : sideID_plus_upper
! communicate sideID_minus_upper : sideID_plus_upper
#ifdef MPI
!CALL StartExchangeMPIData(3,normVecPlus,SideID_Plus_Lower,SideID_Plus_Upper,SendRequest_Flux,RecRequest_Flux,SendID=1)
CALL StartReceiveMPIData(4,normVecPlus,1,nSides,RecRequest_Flux,SendID=1) ! Receive MINE
CALL StartSendMPIData   (4,normVecPlus,1,nSides,SendRequest_Flux,SendID=1) ! Send YOUR
!CALL StartExchangeMPIData(4,normVecPlus,1,nSides,SendRequest_Flux,RecRequest_Flux,SendID=1)
CALL FinishExchangeMPIData(SendRequest_Flux,RecRequest_Flux,SendID=1) !Send MINE -receive YOUR
#endif /*MPI*/

DO iElem=1,PP_nElems
  SideID=ElemToSide(E2S_SIDE_ID,XI_MINUS,iElem)
  Flip  =ElemToSide(E2S_FLIP,XI_MINUS,iElem)
  SELECT CASE(flip)
  CASE(0) !master
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,q,p,XI_MINUS,iElem)            =  NormVec(:,p,q,SideID)
      Surf(q,p,XI_MINUS,iElem)              =  SurfElem(p,q,SideID)
    END DO; END DO 
  CASE(1) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,p,q,XI_MINUS,iElem)            =   NormVecPlus(1:3,p,q,SideID)
      Surf(p,q,XI_MINUS,iElem)              =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  CASE(2) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      !nVecSurf(:,PP_N-q,p,XI_MINUS,iElem)      = -NormVec(:,p,q,SideID)*surfElem(p,q,SideID)
      ! INTERCHANGED WITH FLIP 4!
      nVec(:,q,PP_N-p,XI_MINUS,iElem)        =   NormVecPlus(1:3,p,q,SideID)
      Surf(  q,PP_N-p,XI_MINUS,iElem)        =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  CASE(3) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,PP_N-p,PP_N-q,XI_MINUS,iElem)   =   NormVecPlus(1:3,p,q,SideID)
      Surf(  PP_N-p,PP_N-q,XI_MINUS,iElem)   =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  CASE(4) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      !nVecSurf(:,q,PP_N-p,XI_MINUS,iElem)      = -NormVec(:,p,q,SideID)*surfElem(p,q,SideID)
      ! INTERCHANGED WITH FLIP 2!
      nVec(:,PP_N-q,p,XI_MINUS,iElem)        =   NormVecPlus(1:3,p,q,SideID)
      Surf(  PP_N-q,p,XI_MINUS,iElem)        =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  END SELECT
  SideID=ElemToSide(E2S_SIDE_ID,XI_PLUS,iElem)
  Flip  =ElemToSide(E2S_FLIP,XI_PLUS,iElem)
  SELECT CASE(flip)
  CASE(0) !master
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,p,q,XI_PLUS,iElem)              =   NormVec(:,p,q,SideID)
      Surf(  p,q,XI_PLUS,iElem)              =  SurfElem(p,q,SideID)
    END DO; END DO 
  CASE(1) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,q,p,XI_PLUS,iElem)              =   NormVecPlus(1:3,p,q,SideID)
      Surf(  q,p,XI_PLUS,iElem)              =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  CASE(2) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,PP_N-p,q,XI_PLUS,iElem)         =   NormVecPlus(1:3,p,q,SideID)
      Surf(  PP_N-p,q,XI_PLUS,iElem)         =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  CASE(3) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,PP_N-q,PP_N-p,XI_PLUS,iElem)    =   NormVecPlus(1:3,p,q,SideID)
      Surf(  PP_N-q,PP_N-p,XI_PLUS,iElem)    =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  CASE(4) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,p,PP_N-q,XI_PLUS,iElem)         =   NormVecPlus(1:3,p,q,SideID)
      Surf(  p,PP_N-q,XI_PLUS,iElem)         =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  END SELECT
  SideID=ElemToSide(E2S_SIDE_ID,ETA_MINUS,iElem)
  Flip  =ElemToSide(E2S_FLIP,ETA_MINUS,iElem)
  SELECT CASE(flip)
  CASE(0) !master
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,p,q,ETA_MINUS,iElem)            =   NormVec(:,p,q,SideID)
      Surf(  p,q,ETA_MINUS,iElem)            =  SurfElem(p,q,SideID)
    END DO; END DO 
  CASE(1) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,q,p,ETA_MINUS,iElem)            =   NormVecPlus(1:3,p,q,SideID)
      Surf(  q,p,ETA_MINUS,iElem)            =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  CASE(2) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,PP_N-p,q,ETA_MINUS,iElem)       =   NormVecPlus(1:3,p,q,SideID)
      Surf(  PP_N-q,p,ETA_MINUS,iElem)       =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  CASE(3) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,PP_N-q,PP_N-p,ETA_MINUS,iElem)  =   NormVecPlus(1:3,p,q,SideID)
      Surf(  PP_N-q,PP_N-p,ETA_MINUS,iElem)  =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  CASE(4) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,p,PP_N-q,ETA_MINUS,iElem)       =   NormVecPlus(1:3,p,q,SideID)
      Surf(  p,PP_N-q,ETA_MINUS,iElem)       =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  END SELECT

  SideID=ElemToSide(E2S_SIDE_ID,ETA_PLUS,iElem)
  Flip  =ElemToSide(E2S_FLIP,ETA_PLUS,iElem)
  SELECT CASE(flip)
  CASE(0) !master
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,PP_N-p,q,ETA_PLUS,iElem)        =   NormVec(:,p,q,SideID)
      Surf(  PP_N-p,q,ETA_PLUS,iElem)        =  SurfElem(p,q,SideID)
    END DO; END DO 
  CASE(1) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      !nVecSurf(:,q,PP_N-p,ETA_PLUS,iElem)           = -NormVec(:,p,q,SideID)*surfElem(p,q,SideID)
      ! INTERCHANGED WITH FLIP 3!
      nVec(:,PP_N-q,p,ETA_PLUS,iElem)        =   NormVecPlus(1:3,p,q,SideID)
      Surf(  PP_N-q,p,ETA_PLUS,iElem)        =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  CASE(2) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,p,q,ETA_PLUS,iElem)            =   NormVecPlus(1:3,p,q,SideID)
      Surf(  p,q,ETA_PLUS,iElem)            =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  CASE(3) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      !nVecSurf(:,PP_N-q,p,ETA_PLUS,iElem)           = -NormVec(:,p,q,SideID)*surfElem(p,q,SideID)
      ! INTERCHANGED WITH FLIP 1!
      nVec(:,q,PP_N-p,ETA_PLUS,iElem)       =   NormVecPlus(1:3,p,q,SideID)
      Surf(  q,PP_N-p,ETA_PLUS,iElem)       =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  CASE(4) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,PP_N-p,PP_N-q,ETA_PLUS,iElem)  =   NormVecPlus(1:3,p,q,SideID)
      Surf(  PP_N-p,PP_N-q,ETA_PLUS,iElem)  =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  END SELECT
  SideID=ElemToSide(E2S_SIDE_ID,ZETA_MINUS,iElem)
  Flip  =ElemToSide(E2S_FLIP,ZETA_MINUS,iElem)
  SELECT CASE(flip)
  CASE(0) !master
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,q,p,ZETA_MINUS,iElem)           =   NormVec(:,p,q,SideID)
      Surf(  q,p,ZETA_MINUS,iElem)           =  SurfElem(p,q,SideID)
    END DO; END DO 
  CASE(1) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,p,q,ZETA_MINUS,iElem)           =   NormVecPlus(1:3,p,q,SideID)
      Surf(  p,q,ZETA_MINUS,iElem)           =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  CASE(2) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      !nVecSurf(:,PP_N-q,p,ZETA_MINUS,iElem)      = -NormVec(:,p,q,SideID)*surfElem(p,q,SideID)
      ! INTERCHANGED WITH FLIP 4!
      nVec(:,q,PP_N-p,ZETA_MINUS,iElem)      =   NormVecPlus(1:3,p,q,SideID)
      Surf(  q,PP_N-p,ZETA_MINUS,iElem)      =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  CASE(3) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,PP_N-p,PP_N-q,ZETA_MINUS,iElem) =   NormVecPlus(1:3,p,q,SideID)
      Surf(  PP_N-p,PP_N-q,ZETA_MINUS,iElem) =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  CASE(4) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      !nVecSurf(:,q,PP_N-p,ZETA_MINUS,iElem)      = -NormVec(:,p,q,SideID)*surfElem(p,q,SideID)
      ! INTERCHANGED WITH FLIP 2!
      nVec(:,PP_N-q,p,ZETA_MINUS,iElem)      =   NormVecPlus(1:3,p,q,SideID)
      Surf(  PP_N-q,p,ZETA_MINUS,iElem)      =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  END SELECT
  SideID=ElemToSide(E2S_SIDE_ID,ZETA_PLUS,iElem)
  Flip  =ElemToSide(E2S_FLIP,ZETA_PLUS,iElem)
  SELECT CASE(flip)
  CASE(0) !master
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,p,q,ZETA_PLUS,iElem)           =   NormVec(:,p,q,SideID)
      Surf(  p,q,ZETA_PLUS,iElem)           =  SurfElem(p,q,SideID)
    END DO; END DO 
  CASE(1) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,q,p,ZETA_PLUS,iElem)           =   NormVecPlus(1:3,p,q,SideID)
      Surf(  q,p,ZETA_PLUS,iElem)           =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  CASE(2) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,PP_N-p,q,ZETA_PLUS,iElem)      =   NormVecPlus(1:3,p,q,SideID)
      Surf(  PP_N-p,q,ZETA_PLUS,iElem)      =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  CASE(3) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,PP_N-q,p,ZETA_PLUS,iElem)      =   NormVecPlus(1:3,p,q,SideID)
      Surf(  PP_N-q,p,ZETA_PLUS,iElem)      =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  CASE(4) !slave, flip normal!!
    DO q=0,PP_N; DO p=0,PP_N
      nVec(:,p,PP_N-q,ZETA_PLUS,iElem)      =   NormVecPlus(1:3,p,q,SideID)
      Surf(  p,PP_N-q,ZETA_PLUS,iElem)      =  NormVecPlus(4,p,q,SideID)
    END DO; END DO 
  END SELECT
END DO !iElem

DEALLOCATE(normVecPlus)

BuildNvecisDone=.TRUE.

END SUBROUTINE BuildnVecSurf

SUBROUTINE CheckBJPrecond(Ploc1,Ploc,invPloc,iElem,Time)
!===================================================================================================================================
! Debug routine for checking block jacobian preconditioners
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_LinearSolver_Vars,      ONLY:nDOFelem
USE MOD_Precond_Vars,       ONLY:PrecondType,DebugMatrix
#ifdef MPI
USE MOD_MPI_Vars
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                                      :: Time
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(1:nDOFElem,1:nDOFElem),INTENT(IN)      :: Ploc,invPloc,Ploc1
INTEGER,INTENT(IN)                                    :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
!LOCAL VARIABLES
REAL,DIMENSION(:,:), ALLOCATABLE                      :: diff
REAL                                                  :: dummy
CHARACTER(LEN=255)                                    :: Filename
CHARACTER(LEN=17)                                     :: strfmt
INTEGER                                               :: r,s
REAL                                                  :: Time1,Time2
!===================================================================================================================================

ALLOCATE( diff(1:nDOFElem,1:nDOFElem))

! output of BlockJac  | only derivativ
IF(DebugMatrix.GE.1)THEN
#ifdef MPI
  WRITE(Filename,'(A,I2.2,A,I2.2,A,I4.4,A)')'BlockJac_Rank_',MyRank,'_Type_',PreCondType,'_Mat_', iElem,'.dat'
#else
  WRITE(Filename,'(A,I2.2,A,I4.4,A)')'BlockJac_',PreCondType,'_Mat_', iElem,'.dat'
#endif /*MPI*/
  WRITE(strfmt,'(A1,I4,A12)')'(',nDOFelem,'(1X,E23.16))'
  WRITE(UNIT_stdOut,*)'Debug Block Jacobian to:',TRIM(Filename)
  OPEN (UNIT=103,FILE=TRIM(Filename),STATUS='REPLACE')
  DO r=1,nDOFelem
    WRITE(103,strfmt)Ploc1(r,:)
  END DO
  CLOSE(103)
END IF !DebugMatrix >1

! output of BlockPrecond, no inverse
IF(DebugMatrix.GE.2)THEN

#ifdef MPI
  WRITE(Filename,'(A,I2.2,A,I2.2,A,I4.4,A)')'Precond_Rank_',MyRank,'_Type_',PreCondType,'_Mat_', iElem,'.dat'
#else
  WRITE(Filename,'(A,I2.2,A,I4.4,A)')'Precond_',PreCondType,'_Mat_', iElem,'.dat'
#endif /*MPI*/
  WRITE(strfmt,'(A1,I4,A12)')'(',nDOFelem,'(1X,E23.16))'
  WRITE(UNIT_stdOut,*)'Debug Precond (no Inverse) to:',TRIM(Filename)
  OPEN (UNIT=103,FILE=TRIM(Filename),STATUS='REPLACE')
  DO r=1,nDOFelem
    WRITE(103,strfmt)Ploc(r,:)
  END DO
  CLOSE(103)
END IF !DebugMatrix >1

! output of Inverse
IF(DebugMatrix.GE.3)THEN
#ifdef MPI
  WRITE(Filename,'(A,I2.2,A,I2.2,A,I4.4,A)')'Precond_Rank_',MyRank,'_Type__',PreCondType,'_InvMat_', iElem,'.dat'
#else
  WRITE(Filename,'(A,I2.2,A,I4.4,A)')'Precond_',PreCondType,'_InvMat_', iElem,'.dat'
#endif /*MPI*/
  WRITE(strfmt,'(A1,I4,A12)')'(',nDOFelem,'(1X,E23.16))'
  WRITE(UNIT_stdOut,*)'Debug Precond to:',TRIM(Filename)
  OPEN (UNIT=103,FILE=TRIM(Filename),STATUS='REPLACE')
  DO r=1,nDOFelem
    WRITE(103,strfmt)invPloc(r,:)
  END DO
  CLOSE(103)
END IF !DebugMatrix>2

! sanity check of inverse
IF((DebugMatrix.GE.4).OR.(DebugMatrix.LT.0))THEN
  !CHECK INVERSE
  CALL CPU_TIME(Time1)
  WRITE(UNIT_stdOut,'(A)')'    =>Check invert... '
  dummy=0.
  DO s=1,nDOFelem
    DO r=1,nDOFelem
      dummy=dummy+ABS(SUM(invPloc(r,:)*Ploc(:,s)))
    END DO
  END DO
  dummy=(dummy-REAL(nDOFelem))/REAL(nDOFelem*nDOFelem) !relative error per matrix entry

  IF(dummy.GT. 1.0E-08) THEN
    IPWRITE(UNIT_stdOut,*)'WARNING!!! accuracy problems in with preconditioner inverse..',dummy
  END IF
  IF(dummy.NE. dummy) THEN  !NAN
    IPWRITE(UNIT_stdOut,*)'WARNING!!! NAN problem in with preconditioner inverse..',dummy
    STOP
  END IF
  CALL CPU_TIME(Time2)
  WRITE(UNIT_stdOut,'(A,F11.3,A)')'      ... Check invert done. time=[',Time2-Time1,' ]'
END IF !DebugMatrix>3

Time=Time2-Time1

DEALLOCATE(diff)

END SUBROUTINE CheckBJPrecond
#endif /*maxwell*/

SUBROUTINE FinalizePrecond()
!===================================================================================================================================
! Finalizes variables 
!===================================================================================================================================
! MODULES
USE MOD_Precond_Vars,ONLY:invP,PrecondInitIsDone,PrecondType
USE MOD_Precond_Vars,ONLy:invXi,invEta,invZeta,dRdXi,dRdZeta,dRdEta
USE MOD_SparseILU   ,ONLY:FinalizeSparseILU
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(invP)
!IF(PrecondType.EQ.210)THEN
!  DEALLCOATE(dRdXi)
!  DEALLCOATE(dRdEta)
!  DEALLCOATE(dRdZeta)
!  DEALLCOATE(invXi)
!  DEALLCOATE(invEta)
!  DEALLCOATE(dRdXi)
!!SDEALLCOATE(dRdEta)
!!SDEALLCOATE(dRdZeta)
!!SDEALLCOATE(invXi)
!!SDEALLCOATE(invEta)
!!SDEALLCOATE(invZeta)
!END IF
PrecondInitIsDone = .FALSE.
IF(PrecondType.EQ.22) CALL FinalizeSparseILU

END SUBROUTINE FinalizePrecond


END MODULE MOD_Precond
