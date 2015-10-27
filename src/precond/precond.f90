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
IF(PrecondType.EQ.9) PrecondMethod=3

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
    SWRITE(*,*) ' No Preconditioner'
  CASE(1,2)
    ALLOCATE(invP(1:nDOFelem,1:nDOFelem,PP_nElems))
    CALL InitJac_FD()
    CALL InitJac_Ex()
    CALL InitJacDG()
  CASE(3,22)
    ! sparse ILU0
    CALL InitJac_Ex()
    CALL InitJacDG()
    CALL InitSparseILU()
  CASE(25,26)
    CALL InitJac_Ex()
    CALL InitJacDG()
    CALL InitILU()
  CASE(9)
    ! use GMRES for each element
    ! PP_nVar x PP_nVar Jacobi
    ALLOCATE( invBJ (1:PP_nVar,1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)  )
    CALL InitJac_Ex()
    CALL InitJacDG() ! here we build dR/dU U operator
  CASE(200,201,202,203)
    ! 1D-Tensor BJ preconditioner
    locSize=PP_nVar*(PP_N+1)
    ALLOCATE( invXi  (locSize,locSize,0:PP_N,0:PP_N,1:PP_nElems) &
            , invEta (locSize,locSize,0:PP_N,0:PP_N,1:PP_nElems) &
            , invZeta(locSize,locSize,0:PP_N,0:PP_N,1:PP_nElems) )
    ALLOCATE( dRdXi  (locSize,locSize,0:PP_N,0:PP_N,1:1) &
            , dRdEta (locSize,locSize,0:PP_N,0:PP_N,1:1) &
            , dRdZeta(locSize,locSize,0:PP_N,0:PP_N,1:1) )
    CALL InitJac_Ex()
    CALL InitJacDG()
  CASE(205,210,211,212)
    ! 1D-Tensor BJ preconditioner
    locSize=PP_nVar*(PP_N+1)
    ALLOCATE( invXi  (locSize,locSize,0:PP_N,0:PP_N,1:PP_nElems) &
            , invEta (locSize,locSize,0:PP_N,0:PP_N,1:PP_nElems) &
            , invZeta(locSize,locSize,0:PP_N,0:PP_N,1:PP_nElems) )
    ALLOCATE( dRdXi  (locSize,locSize,0:PP_N,0:PP_N,1:PP_nElems) &
            , dRdEta (locSize,locSize,0:PP_N,0:PP_N,1:PP_nElems) &
            , dRdZeta(locSize,locSize,0:PP_N,0:PP_N,1:PP_nElems) )
    CALL InitJac_Ex()
    CALL InitJacDG()

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
USE MOD_LinearSolver_Vars ,ONLY:nDOFelem,nDOFLine
USE MOD_Precond_Vars      ,ONLY:invP,PrecondType,DebugMatrix,invBJ,invJ, neighborElemID, ProcJacobian
USE MOD_Jac_ex            ,ONLY:Jac_ex, Jac_Ex_Neighbor,Jac_ex1D
USE MOD_Jac_FD            ,ONLY:Jac_FD_slow!,Jac_FD
USE MOD_Jac_FD_Vars,       ONLY:reps0,XK
USE MOD_JacDG             ,ONLY:JacBJDOF
USE MOD_Mesh_Vars         ,ONLY:nBCSides,ElemToSide,SideToElem
!USE MOD_Prepare_FD    ,ONLY:PrepareFD
USE MOD_Precond_Vars      ,ONLy:invXi,invEta,invZeta,dRdXi,dRdZeta,dRdEta
USE MOD_DG,                ONLY: DGTimeDerivative_WeakForm
USE MOD_DG_Vars,           ONLY: U,Ut
USE MOD_SparseILU,         ONLY: BuildILU0
USE MOD_ILU,               ONLY: BuildBILU0,BuildBILU0BCSR
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
  CASE(1,2,3)
    ALLOCATE( Ploc (1:nDOFElem,1:nDOFElem), &
              Ploc1(1:nDOFElem,1:nDOFElem) )
  CASE(9)
    ALLOCATE( Ploc (1:nDOFElem,1:nDOFElem))
  CASE(22,25,26,30,31) 
    ALLOCATE( Ploc (1:nDOFElem,1:nDOFElem))
  CASE(202,205,210,211)
    ALLOCATE( Ploc (1:nDOFElem,1:nDOFElem), &
              Ploc1(1:nDOFElem,1:nDOFElem))!&
              !Ploc2(1:nDOFElem,1:nDOFElem)  )
END SELECT

TotalTime=0.

! selection of precond type
! PrecondType .LT. 99 - element local
IF( PrecondType.LT.99) THEN ! element local preconditoner
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
    CALL CPU_TIME(TimeEnd(1))

    CALL CPU_TIME(TimeStart(2))
    SELECT CASE(PrecondType)
    CASE(1,2) ! element block jacobi
      !invert Ploc => invP(:,:,iElem)
      invP(:,:,iElem)=getInverse(nDOFelem,Ploc)
    CASE(3) ! inverse of PP_nVarxPP_nVar block
      vn1=PP_nVar*(PP_N+1)
      vn2=vn1*(PP_N+1)
      ! alpha*dt and the calcution with the unity matrix is already performed
      !CALL CPU_TIME(TimeStart(2))
      DO k=0,PP_N
        DO j=0,PP_N
          DO i=0,PP_N
            s = vn2 * k + vn1 * j + PP_nVar * i
            invBJ(:,:,i,j,k,iElem)=getInverse(PP_nVar,Ploc(s+1:s+PP_nVar,s+1:s+PP_nVar))
          END DO ! i
        END DO ! j
      END DO ! k
    CASE(22) ! ilu(0) of element
      CALL BuildILU0(Ploc,iElem)
    CASE(25)
      CALL BuildBILU0(Ploc,iElem)
    CASE(26)
      CALL BuildBILU0BCSR(Ploc,iElem)
    END SELECT
    CALL CPU_TIME(TimeEnd(2))
    TotalTime=TotalTime+(TimeEnd-TimeStart)
    IF((PrecondType.LE.2).AND.DebugMatrix.NE.0)THEN
      CALL CheckBJPrecond(Ploc1,Ploc,invP(:,:,iElem),iElem,TotalTime(3))
    END IF ! DebugMatrix
  END DO ! ! iELEM
ELSEIF(PrecondType.GT.199)THEN
  coeff=-alpha*dt
  DO iElem=1,PP_nElems
    CALL CPU_TIME(TimeStart(1))
    SELECT CASE(PrecondType)
    CASE (201,202,203)
      dRdXi  =0.
      dRdEta =0.
      dRdZeta=0.
      CALL Jac_ex1D(dRdXi(:,:,:,:,1),dRdEta(:,:,:,:,1),dRdZeta(:,:,:,:,1),iElem)
      !SELECT CASE(PrecondType)
      !CASE(200,201,203)
        DO q=0,PP_N
          DO p=0,PP_N
            DO s=1,nDOFLine
              DO r=1,nDOFLine
                dRdXi  (r,s,p,q,1)=coeff*dRdXi  (r,s,p,q,1)
                dRdEta (r,s,p,q,1)=coeff*dRdEta (r,s,p,q,1)
                dRdZeta(r,s,p,q,1)=coeff*dRdZeta(r,s,p,q,1)
              END DO ! r
              dRdXi  (s,s,p,q,1)=dRdXi  (s,s,p,q,1) +1.!/3.
              dRdEta (s,s,p,q,1)=dRdEta (s,s,p,q,1) +1.!/3.
              dRdZeta(s,s,p,q,1)=dRdZeta(s,s,p,q,1) +1.!/3.
            END DO ! s
          END DO ! p
        END DO !q
      !CASE(202)
      !  DO q=0,PP_N
      !    DO p=0,PP_N
      !      DO s=1,nDOFLine
      !        DO r=1,nDOFLine
      !          dRdXi  (r,s,p,q,1)=coeff*dRdXi  (r,s,p,q,1)
      !          dRdEta (r,s,p,q,1)=coeff*dRdEta (r,s,p,q,1)
      !          dRdZeta(r,s,p,q,1)=coeff*dRdZeta(r,s,p,q,1)
      !        END DO ! r
      !        dRdXi  (s,s,p,q,1)=dRdXi  (s,s,p,q,1) +1.!/3.
      !        dRdEta (s,s,p,q,1)=dRdEta (s,s,p,q,1) +1.!/3.
      !        dRdZeta(s,s,p,q,1)=dRdZeta(s,s,p,q,1) +1.!/3.
      !      END DO ! s
      !    END DO ! p
      !  END DO !q
      !END SELECT
    CASE(205,210,211)
      ! nullify
      DO s=1,nDOFelem
        DO r=1,nDOFelem
          Ploc(r,s)=0.
        END DO !r
      END DO !s
      ! analytic per Element
      CALL Jac_ex(iElem,Ploc)
      ! add contibution I-alpha*dt*dRdU
      DO s=1,nDOFelem
        DO r=1,nDOFelem
          Ploc(r,s)=coeff*Ploc(r,s)
        END DO !r
        Ploc(s,s)=Ploc(s,s)+1.
      END DO !s
    CASE(212)
      dRdXi  (:,:,:,:,iElem)=0.
      dRdEta (:,:,:,:,iElem)=0.
      dRdZeta(:,:,:,:,iElem)=0.
      CALL Jac_ex1D(dRdXi(:,:,:,:,iElem),dRdEta(:,:,:,:,iElem),dRdZeta(:,:,:,:,iElem),iElem)
      !SELECT CASE(PrecondType)
      !CASE(200,201,203)
        DO q=0,PP_N
          DO p=0,PP_N
            DO s=1,nDOFLine
              DO r=1,nDOFLine
                dRdXi  (r,s,p,q,iElem)=coeff*dRdXi  (r,s,p,q,iElem)
                dRdEta (r,s,p,q,iElem)=coeff*dRdEta (r,s,p,q,iElem)
                dRdZeta(r,s,p,q,iElem)=coeff*dRdZeta(r,s,p,q,iElem)
              END DO ! r
              dRdXi  (s,s,p,q,iElem)=dRdXi  (s,s,p,q,iElem) +1.!/3.
              dRdEta (s,s,p,q,iElem)=dRdEta (s,s,p,q,iElem) +1.!/3.
              dRdZeta(s,s,p,q,iElem)=dRdZeta(s,s,p,q,iElem) +1.!/3.
            END DO ! s
          END DO ! p
        END DO !q
    CASE DEFAULT
      CALL abort(__STAMP__,&
          ' Preconditioner not implemented.',PrecondType,999.)

    END SELECT

    CALL CPU_TIME(TimeEnd(1))

    CALL CPU_TIME(TimeStart(2))
    SELECT CASE(PrecondType)
    CASE(200,201)
      ! build invXi,invEta and invZeta
      ! compute 1D inverses
      DO q=0,PP_N
        DO p=0,PP_N
          invXi  (:,:,p,q,iElem)=getInverse(nDOFLine,dRdXi  (:,:,p,q,1))
          invEta (:,:,p,q,iElem)=getInverse(nDOFLine,dRdEta (:,:,p,q,1))
          invZeta(:,:,p,q,iElem)=getInverse(nDOFLine,dRdZeta(:,:,p,q,1))
        END DO ! j
      END DO ! k
    CASE(203)
      ! build invXi,invEta and invZeta
      ! compute 1D inverses
      DO q=0,PP_N
        DO p=0,PP_N
          tmp=getInverse(nDOFLine,dRdXi  (:,:,p,q,1))
          DO l=1,nDOFLine
            invXi  (:,l,p,q,iElem) = tmp(l,:)
          END DO ! l
          tmp=getInverse(nDOFLine,dRdEta (:,:,p,q,1))
          DO l=1,nDOFLine
            invEta (:,l,p,q,iElem) = tmp(l,:)
          END DO ! l
          tmp=getInverse(nDOFLine,dRdZeta(:,:,p,q,1))
          DO l=1,nDOFLine
            invZeta  (:,l,p,q,iElem) = tmp(l,:)
          END DO ! l
        END DO ! j
      END DO ! k
    CASE(202)
      ! build invXi,invEta and invZeta
      ! compute 1D inverses
      DO q=0,PP_N
        DO p=0,PP_N
          invXi  (:,:,p,q,iElem)=getInverse(nDOFLine,dRdXi  (:,:,p,q,1))
          invZeta(:,:,p,q,iElem)=getInverse(nDOFLine,dRdZeta(:,:,p,q,1))
        END DO ! j
      END DO ! k
      ! map xi derivatives to r,s
      vn1 = PP_nVar * (PP_N + 1)
      vn2 = vn1 * (PP_N +1)
      Ploc=0.
      Ploc1=0.
      DO k=0,PP_N
        DO j=0,PP_N
          DO i=0,PP_N
            vni=PP_nVar*i
            r = vn2 * k + vn1 * j + PP_nVar * i
            DO l=0,PP_N
              vnl=l*PP_nVar
              s = vn2 * k + vn1 * j + PP_nVar * l
              Ploc(r+1:r+PP_nVar,s+1:s+PP_nVar) = dRdXi    (vni+1:vni+PP_nVar,vnl+1:vnl+PP_nVar,j,k,1)
            END DO ! l
          END DO ! i
        END DO ! j
      END DO ! k
      ! and zeta
      DO k=0,PP_N
        vnk=PP_nVar*k
        DO j=0,PP_N
          DO i=0,PP_N
            r = vn2 * k + vn1 * j + PP_nVar * i
            DO l=0,PP_N
              vnl=l*PP_nVar
              s = vn2 * l + vn1 * j + PP_nVar * i
              Ploc1(r+1:r+PP_nVar,s+1:s+PP_nVar) = dRdZeta    (vnk+1:vnk+PP_nVar,vnl+1:vnl+PP_nVar,i,j,1)
            END DO ! l
          END DO ! i
        END DO ! j
      END DO ! k
      ! rotate all into eta system
      !IF(MPIRoot)THEN
      !  IF(iElem.EQ.1)THEN
      !    WRITE(strfmt,'(A1,I4,A12)')'(',nDOFelem,'(1X,E23.16))'
      !    WRITE(UNIT_stdOut,*)'Debug Block Jacobian to:'
      !    OPEN (UNIT=103,FILE='MPIRank0_Elem1_Xi.dat',STATUS='REPLACE')
      !    DO r=1,nDOFelem
      !      WRITE(103,strfmt)Ploc(r,:)
      !    END DO
      !    CLOSE(103)
      !    OPEN (UNIT=103,FILE='MPIRank0_Elem1_Zeta.dat',STATUS='REPLACE')
      !    DO r=1,nDOFelem
      !      WRITE(103,strfmt)Ploc1(r,:)
      !    END DO
      !    CLOSE(103)
      !    Ploc=0.
      !    DO k=0,PP_N
      !      DO j=0,PP_N
      !        vnj=PP_nVar*j
      !        DO i=0,PP_N
      !          r = vn2 * k + vn1 * j + PP_nVar * i
      !          DO l=0,PP_N
      !            vnl=l*PP_nVar
      !            s = vn2 * k + vn1 * l + PP_nVar * i
      !            Ploc(r+1:r+PP_nVar,s+1:s+PP_nVar) = dRdEta    (vnj+1:vnj+PP_nVar,vnl+1:vnl+PP_nVar,i,k)
      !          END DO ! l
      !        END DO ! i
      !      END DO ! j
      !    END DO ! k
      !    OPEN (UNIT=103,FILE='MPIRank0_Elem1_Eta.dat',STATUS='REPLACE')
      !    DO r=1,nDOFelem
      !      WRITE(103,strfmt)Ploc(r,:)
      !    END DO
      !    CLOSE(103)
      !  END IF
      !END IF

      ! order xi and zeta after eta
      dRdXi  =0.
      dRdZeta=0.
      DO k=0,PP_N
        DO j=0,PP_N
          vnj=PP_nVar*j
          DO i=0,PP_N
            r = vn2 * k + vn1 * j + PP_nVar * i
            DO l=0,PP_N
              vnl=l*PP_nVar
              s = vn2 * k + vn1 * l + PP_nVar * i
              dRdXi  (vnj+1:vnj+PP_nVar,vnl+1:vnl+PP_nVar,i,k,1)=Ploc (r+1:r+PP_nVar,s+1:s+PP_nVar) 
              dRdZeta(vnj+1:vnj+PP_nVar,vnl+1:vnl+PP_nVar,i,k,1)=Ploc1(r+1:r+PP_nVar,s+1:s+PP_nVar) 
            END DO ! l
          END DO ! i
        END DO ! j
      END DO ! k

!      ! compute 1D diagonal blocks 8x8
!      DO q=0,PP_N
!        DO p=0,PP_N
!          DO l=0,PP_N
!            vnl=PP_nVar*l
!            dRdXi  (vnl+1:vnl+PP_nVar,vnl+1:vnl+PP_nVar,p,q)=&
!                                                              getInverse(PP_nVar,dRdXi  (vnl+1:vnl+PP_nVar,vnl+1:vnl+PP_nVar,p,q))
!            dRdZeta(vnl+1:vnl+PP_nVar,vnl+1:vnl+PP_nVar,p,q)=&
!                                                              getInverse(PP_nVar,dRdZeta(vnl+1:vnl+PP_nVar,vnl+1:vnl+PP_nVar,p,q))
!          END DO ! l
!        END DO ! j
!      END DO ! k
!      ! compute stuff
!      DO q=0,PP_N
!        DO p=0,PP_N
!          DO l=0,PP_N
!            dRdEta(vnl+1:vnl+PP_nVar,vnl+1:vnl+PP_nVar,p,q)=MATMUL(dRdEta (vnl+1:vnl+PP_nVar,vnl+1:vnl+PP_nVar,p,q)&
!                                                                  ,dRdZeta(vnl+1:vnl+PP_nVar,vnl+1:vnl+PP_nVar,p,q))
!            dRdEta(vnl+1:vnl+PP_nVar,vnl+1:vnl+PP_nVar,p,q)=MATMUL(dRdXi  (vnl+1:vnl+PP_nVar,vnl+1:vnl+PP_nVar,p,q)&
!                                                                  ,dRdEta (vnl+1:vnl+PP_nVar,vnl+1:vnl+PP_nVar,p,q))
!          END DO ! l
!        END DO ! j
!      END DO ! k

      ! compute 1D diagonal blocks 8x8
      DO q=0,PP_N
        DO p=0,PP_N
          dRdXi  (:,:,p,q,1)=getInverse(nDOFLine,dRdXi  (:,:,p,q,1))
          dRdZeta(:,:,p,q,1)=getInverse(nDOFLine,dRdZeta(:,:,p,q,1))
        END DO ! p
      END DO ! q
      ! compute stuff
      DO q=0,PP_N
        DO p=0,PP_N
          dRdEta(:,:,p,q,1)=MATMUL(dRdEta (:,:,p,q,1),dRdZeta(:,:,p,q,1))
          dRdEta(:,:,p,q,1)=MATMUL(dRdXi  (:,:,p,q,1),dRdEta (:,:,p,q,1))
        END DO ! p
      END DO ! q

      ! add rest
      dRdEta=dRdEta+dRdZeta+dRdXi
       ! finally, compte eta-inv
       ! compute 1D inverses
       DO q=0,PP_N
         DO p=0,PP_N
           invEta (:,:,p,q,iElem)=getInverse(nDOFLine,dRdEta (:,:,p,q,1))
         END DO ! j
       END DO ! k
      !IF(MPIRoot)THEN
      !  IF(iElem.EQ.1)THEN
      !    WRITE(strfmt,'(A1,I4,A12)')'(',nDOFelem,'(1X,E23.16))'
      !    WRITE(UNIT_stdOut,*)'Debug Block Jacobian to:'
      !    Ploc=0.
      !    DO k=0,PP_N
      !      DO j=0,PP_N
      !        vnj=PP_nVar*j
      !        DO i=0,PP_N
      !          r = vn2 * k + vn1 * j + PP_nVar * i
      !          DO l=0,PP_N
      !            vnl=l*PP_nVar
      !            s = vn2 * k + vn1 * l + PP_nVar * i
      !            Ploc(r+1:r+PP_nVar,s+1:s+PP_nVar) = invEta    (vnj+1:vnj+PP_nVar,vnl+1:vnl+PP_nVar,i,k,iElem)
      !          END DO ! l
      !        END DO ! i
      !      END DO ! j
      !    END DO ! k
      !    OPEN (UNIT=103,FILE='MPIRank0_Elem1_Etainv.dat',STATUS='REPLACE')
      !    DO r=1,nDOFelem
      !      WRITE(103,strfmt)Ploc(r,:)
      !    END DO
      !    CLOSE(103)
      !  END IF
      !END IF
    CASE(205,210,211)
      vn1 = PP_nVar * (PP_N + 1)
      vn2 = vn1 * (PP_N +1)
      DO k=0,PP_N
        vnk=PP_nVar*k
        DO j=0,PP_N
          vnj=PP_nVar*j
          DO i=0,PP_N
            vni=PP_nVar*i
            r = vn2 * k + vn1 * j + PP_nVar * i
            DO l=0,PP_N
              vnl=l*PP_nVar
              ! xi
              s = vn2 * k + vn1 * j + PP_nVar * l
              dRdXi  (vni+1:vni+PP_nVar,vnl+1:vnl+PP_nVar,j,k,iElem)=Ploc (r+1:r+PP_nVar,s+1:s+PP_nVar) 
              ! eta
              s = vn2 * k + vn1 * l + PP_nVar * i
              dRdEta (vnj+1:vnj+PP_nVar,vnl+1:vnl+PP_nVar,i,k,iElem)=Ploc (r+1:r+PP_nVar,s+1:s+PP_nVar) 
              ! zeta
              s = vn2 * l + vn1 * j + PP_nVar * i
              dRdZeta(vnk+1:vnk+PP_nVar,vnl+1:vnl+PP_nVar,i,j,iElem)=Ploc (r+1:r+PP_nVar,s+1:s+PP_nVar) 
            END DO ! l
          END DO ! i
        END DO ! j
      END DO ! k
      ! compute 1D inverses
      DO q=0,PP_N
        DO p=0,PP_N
          invXi  (:,:,p,q,iElem)=getInverse(nDOFLine,dRdXi  (:,:,p,q,iElem))
          invEta (:,:,p,q,iElem)=getInverse(nDOFLine,dRdEta (:,:,p,q,iElem))
          invZeta(:,:,p,q,iElem)=getInverse(nDOFLine,dRdZeta(:,:,p,q,iElem))
        END DO ! j
      END DO ! k
      IF(PrecondType.GE.210)THEN
        vn1 = PP_nVar * (PP_N + 1)
        vn2 = vn1 * (PP_N +1)
        DO k=0,PP_N
          vnk=PP_nVar*k
          DO j=0,PP_N
            vnj=PP_nVar*j
            DO i=0,PP_N
              vni=PP_nVar*i
              r = vn2 * k + vn1 * j + PP_nVar * i
              DO l=0,PP_N
                vnl=l*PP_nVar
                ! xi
                s = vn2 * k + vn1 * j + PP_nVar * l
                dRdXi  (vni+1:vni+PP_nVar,vnl+1:vnl+PP_nVar,j,k,iElem)=Ploc (r+1:r+PP_nVar,s+1:s+PP_nVar) 
                Ploc (r+1:r+PP_nVar,s+1:s+PP_nVar)=0.
                ! eta
                s = vn2 * k + vn1 * l + PP_nVar * i
                dRdEta (vnj+1:vnj+PP_nVar,vnl+1:vnl+PP_nVar,i,k,iElem)=Ploc (r+1:r+PP_nVar,s+1:s+PP_nVar) 
                Ploc (r+1:r+PP_nVar,s+1:s+PP_nVar)=0.
                ! zeta
                s = vn2 * l + vn1 * j + PP_nVar * i
                dRdZeta(vnk+1:vnk+PP_nVar,vnl+1:vnl+PP_nVar,i,j,iElem)=Ploc (r+1:r+PP_nVar,s+1:s+PP_nVar) 
                Ploc (r+1:r+PP_nVar,s+1:s+PP_nVar)=0.
              END DO ! l
            END DO ! i
          END DO ! j
        END DO ! k
      END IF ! PrecondType .GE. 210
    CASE(212)
      ! build invXi,invEta and invZeta
      ! compute 1D inverses
      DO q=0,PP_N
        DO p=0,PP_N
          invXi  (:,:,p,q,iElem)=getInverse(nDOFLine,dRdXi  (:,:,p,q,iElem))
          invEta (:,:,p,q,iElem)=getInverse(nDOFLine,dRdEta (:,:,p,q,iElem))
          invZeta(:,:,p,q,iElem)=getInverse(nDOFLine,dRdZeta(:,:,p,q,iElem))
        END DO ! j
      END DO ! k
    CASE DEFAULT
      CALL abort(__STAMP__,&
          ' Preconditioner not implemented.',PrecondType,999.)
    END SELECT
    CALL CPU_TIME(TimeEnd(2))
    TotalTime=TotalTime+(TimeEnd-TimeStart)
  END DO ! iElem
END IF

WRITE(UNIT_stdOut,'(A,F11.3,A)')' TOTAL DERIVATING TIME =[',TotalTime(1),' ]'
WRITE(UNIT_stdOut,'(A,F11.3,A)')' TOTAL INVERTING  TIME =[',TotalTime(2),' ]'
IF(debugMatrix.GT.2) WRITE(UNIT_stdOut,'(A,F11.3,A)')' TOTAL CHECKING  TIME =[',TotalTime(3),' ]'
SWRITE(UNIT_stdOut,'(A)')' BUILD PRECONDITIONER DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')


SELECT CASE(PrecondType)
  CASE(1,2) 
    DEALLOCATE( Ploc, Ploc1)
  CASE(3)
    DEALLOCATE( Ploc, Ploc1)
  CASE(22) 
    DEALLOCATE( Ploc )
  CASE(200,201,203,205) 
    DEALLOCATE( dRdXi,dRdEta,dRdZeta)
  CASE(202) 
    DEALLOCATE( dRdXi,dRdEta,dRdZeta,Ploc,Ploc1)!,Ploc2)
  CASE(100)
    DEALLOCATE( ProcJacobian )
!    DEALLOCATE( GlobalJacobian    )   
!    DEALLOCATE( GlobalJacobian2   )  
!    DEALLOCATE( GlobalJacobianTild) 
END SELECT

END SUBROUTINE  BuildPrecond


SUBROUTINE Transform(GlobalJac,uhelp,nDOFElem)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(OUT)             :: GlobalJac(nDOFElem)
REAL,INTENT(IN)              :: uhelp(nDOFElem)
INTEGER,INTENT(IN)           :: nDOFElem!,nTotalDOF
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

GlobalJac(:)=uhelp(:)

END SUBROUTINE Transform

SUBROUTINE SortMatrix()
!===================================================================================================================================
! Sort the matrix to get the element-pattern of the Jacoabian Matrix
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Precond_Vars,      ONLY: MatPattern,MatPatternElemID,neighborElemID,Lower,Upper
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER       :: length=7
INTEGER                 :: ElemList(1:length)
INTEGER                 :: iList,iElem,ilocSide
LOGICAL                 :: first
!===================================================================================================================================

Lower=-1
Upper=-1
DO iElem=1,PP_nElems
  ! fill ElemList
  ElemList(1) = iElem
  DO ilocSide=1,6
    ElemList(ilocSide+1) = neighborElemID(ilocSide,iElem)
  END DO ! ilocSide
  CALL BubbleSort(ElemList)
  MatPatternElemID(:,iElem)=ElemList
  DO iList=1,length
    IF(ElemList(iList).EQ.-1) CYCLE
    IF(iElem.EQ.ElemList(iList)) MatPattern(iList,iElem) = 0
    DO ilocSide=1,6
      IF(neighborElemID(ilocSide,iElem).EQ.ElemList(iList)) MatPattern(iList,iElem)=ilocSide
    END DO ! ilocSide
  END DO ! iList
  ! now, fill lower index range
  first=.TRUE.
  DO iList=1,length
    IF(ElemList(iList).EQ.-1) CYCLE
    IF(first.AND.ElemList(iList).LT.iElem)THEN
      Lower(1,iElem)=iList
      first=.FALSE.
    END IF
    IF(.NOT.first.AND.ElemList(iList).LT.iElem)THEN
      Lower(2,iElem)=iList
    END IF
  END DO ! iList
  ! and upper
  first=.TRUE.
  DO iList=1,length
    IF(ElemList(iList).EQ.-1) CYCLE
    IF(first.AND.ElemList(iList).GT.iElem)THEN
      Upper(1,iElem)=iList
      first=.FALSE.
    END IF
    IF(.NOT.first.AND.ElemList(iList).GT.iElem)THEN
      Upper(2,iElem)=iList
    END IF
  END DO ! iList
  ! correct for minus 1
  IF(Lower(1,iElem).EQ.-1) Lower(2,iElem)=-2
  IF(Upper(1,iElem).EQ.-1) Upper(2,iElem)=-2
!  WRITE(*,'(A,8I6)') ' iElem and  MatPattern', iElem, MatPatternElemID(:,iElem)
!  WRITE(*,'(A,8I6)') ' iElem and  MatPattern', iElem, MatPattern(:,iElem)
!  WRITE(*,'(A,2I6)') ' Lower ',Lower(:,iElem)
!  WRITE(*,'(A,2I6)') ' Upper ',Upper(:,iElem)
END DO ! iElem

END SUBROUTINE SortMatrix

SUBROUTINE BubbleSort(SortList)
!===================================================================================================================================
! code copied & modified from rosettacode.org
! modified for integer array with output of sorted array, unsorted array remains
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)   :: SortList(1:7)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER       :: length=7
INTEGER                 :: temp
INTEGER                 :: i, j
LOGICAL                 :: swapped = .TRUE.
!===================================================================================================================================

DO j = length-1, 1, -1
  swapped = .FALSE.
  DO i = 1, j
    IF (sortlist(i) > sortlist(i+1)) THEN
      temp = sortlist(i)
      sortlist(i) = sortlist(i+1)
      sortlist(i+1) = temp
      swapped = .TRUE.
    END IF
  END DO
  IF (.NOT. swapped) EXIT
END DO

END SUBROUTINE BubbleSort

SUBROUTINE WriteJacobian(GlobalJacobian)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_LinearSolver_Vars,      ONLY:nDOFelem
USE MOD_Precond_Vars,           ONLY:neighborElemID
USE MOD_Mesh_Vars         ,ONLY:nBCSides,ElemToSide,SideToElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)              :: GlobalJacobian (nDOFElem,nDOFElem,0:6,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: ElemList(0:6),SideID
INTEGER(KIND=8)     :: r,s,iList,iElem,NBElemID,ilocSide, rout,sout
CHARACTER(LEN=255)           :: filename
REAL                         :: epsZero
!===================================================================================================================================

  epsZero=1e-12
  WRITE(Filename,'(A,I2.2,A)')'MyRank_',myRank,'_Jac_SparseM.dat'
  WRITE(*,*)'Debug Jacobian to:',TRIM(Filename)
  OPEN (UNIT=103,FILE=TRIM(Filename),STATUS='REPLACE')


  ! output of each elem row
  DO iElem=1,PP_nElems
    ! fill ElemList
    ElemList(0) = iElem
    DO ilocSide=1,6
      ElemList(ilocSide) = neighborElemID(ilocSide,iElem)
    END DO ! ilocSide
    CALL BubbleSort(ElemList)
!    print*,'ElemID',iElem
!    print*,ElemList
!    print*,'----'
    DO iList=0,6
      IF(iElem.EQ.ElemList(iList))THEN
        ! output of ielem
      !  print*,'iElem',iElem
        rout=(iElem-1)*nDOFElem
        sout=rout
      !  print*,'rout,sout',rout,sout
        DO r=1,nDOFElem
          DO s=1, nDOFElem
            IF(ABS(GlobalJacobian(r,s,0,iElem)).GT.epsZero)THEN
              WRITE(103,'(I12,2x,I12,2x,E23.16)') rout+r,sout+s,GlobalJacobian(r,s,0,iElem)
            END IF
          END DO !s 
        END DO !r
      ELSE
        DO ilocSide=1,6
          SideID = ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
          IF(SideID.LE.nBCSides) CYCLE
          NBElemID=neighborElemID(ilocSide,iElem)
          !print*,'NBElemID',NBElemID
          IF(NBElemID.EQ.ElemList(iList))THEN
          !  print*,'valid iElem,NBElemID',iElem,NBElemID
            ! output of NBElemID
            rout=(iElem-1)*nDOFElem
            sout=(NBElemID-1)*nDOFElem
          !  print*,'rout,sout',rout,sout
           DO r=1,nDOFElem
             DO s=1, nDOFElem
               IF(ABS(GlobalJacobian(r,s,ilocSide,iElem)).GT.epsZero)THEN
                 WRITE(103,'(I12,2x,I12,2x,E23.16)') rout+r,sout+s,GlobalJacobian(r,s,ilocSide,iElem)
               END IF
             END DO !s 
           END DO !r
          END IF
        END DO ! ilocSide2
      END IF 
      !read*
    END DO ! ilocSide
  END DO ! iElem

!    NBElemID=neighborElemID(ilocSide,iElem)
!  DO r=1,nTotalDOF
!    DO s=1, nTotalDOF
!      IF(ABS(GlobalJacobian(r,s)).GT.epsZero)THEN
!        WRITE(103,'(I6,2x,I6,2x,E23.16)') r,s,GlobalJacobian(r,s)
!      END IF
!    END DO !s 
!  END DO !r

  CLOSE(103)

END SUBROUTINE WriteJacobian

SUBROUTINE WriteGlobalJacobian(GlobalJacobian,GlobalJacobian2,nTotalDOF)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)              :: GlobalJacobian (nTotalDOF,nTotalDOF)
REAL,INTENT(IN)              :: GlobalJacobian2(nTotalDOF,nTotalDOF)
INTEGER,INTENT(IN)           :: nTotalDOF
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: r,s
CHARACTER(LEN=255)           :: filename
CHARACTER(LEN=17)  :: strfmt
REAL                         :: epsZero
!===================================================================================================================================

  epsZero=1e-12
  WRITE(Filename,'(A,I2.2,A,I4.4,A)')'GlobalJac_SparseM.dat'
  !WRITE(strfmt,'(A1,I4,A12)')'(',nTotalDOF,'(1X,E23.16))'
  WRITE(*,*)'Debug Jacobian to:',TRIM(Filename)
  OPEN (UNIT=103,FILE=TRIM(Filename),STATUS='REPLACE')

  DO r=1,nTotalDOF
    DO s=1, nTotalDOF
      IF(ABS(GlobalJacobian(r,s)).GT.epsZero)THEN
        WRITE(103,'(I6,2x,I6,2x,E23.16)') r,s,GlobalJacobian(r,s)
      END IF
    END DO !s 
  END DO !r

  !DO r=1,nTotalDOF
  !  WRITE(103,strfmt)GlobalJacobian(r,:)
  !END DO
  CLOSE(103)

!  WRITE(Filename,'(A,I2.2,A,I4.4,A)')'GlobalJac_corr_SparseM.dat'
!  !WRITE(strfmt,'(A1,I4,A12)')'(',nTotalDOF,'(1X,E23.16))'
!  WRITE(*,*)'Debug Jacobian to:',TRIM(Filename)
!  OPEN (UNIT=103,FILE=TRIM(Filename),STATUS='REPLACE')
!
!  DO r=1,nTotalDOF
!    DO s=1, nTotalDOF
!      IF(ABS(GlobalJacobian2(r,s)).GT.epsZero)THEN
!        WRITE(103,'(I6,2x,I6,2x,E23.16)') r,s,GlobalJacobian2(r,s)
!      END IF
!    END DO !s 
!  END DO !r

  !DO r=1,nTotalDOF
  !  WRITE(103,strfmt)GlobalJacobian2(r,:)
  !END DO
  CLOSE(103)
END SUBROUTINE WriteGlobalJacobian


!SUBROUTINE ApplyPrecond_Elem(v,z)
!!=================================================================================================================================
!! Apply BJ Preconditioner which is constructed per element
!!=================================================================================================================================
!! MODULES
!USE MOD_PreProc
!USE MOD_LinearSolver_Vars ,ONLY:nDOFelem
!USE MOD_Precond_Vars      ,ONLY:invP,PrecondMethod
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!!---------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!REAL,INTENT(IN)              :: v(nDOFelem,PP_nElems)
!!---------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!REAL,INTENT(OUT)             :: z(nDOFelem,PP_nElems)
!!---------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                      :: r,s,iElem
!!=================================================================================================================================
!
!SELECT CASE(PrecondMethod)
!  CASE(0) ! using a loop
!    ! loop over all elements
!    DO iElem=1,PP_nElems
!      z(:,iElem) = 0.
!      DO s=1,nDOFElem
!        DO r=1,nDOFElem
!          z(r,iElem) = z(r,iElem)+invP(r,s,iElem)*v(s,iElem)
!        END DO ! r
!      END DO ! s
!    END DO ! iElem
!  CASE(1) ! Matmul
!    DO iElem=1,PP_nElems
!      z(:,iElem) = MATMUL(invP(:,:,iElem),v(:,iElem))
!    END DO !iElem
!END SELECT
!
!END SUBROUTINE  ApplyPrecond_Elem
!
!SUBROUTINE ApplyPrecond_DOF(v,z)
!!=================================================================================================================================
!! Apply BJ Preconditioner which is constructed per DOF
!!=================================================================================================================================
!! MODULES
!USE MOD_PreProc
!USE MOD_LinearSolver_Vars ,ONLY:nDOFelem
!USE MOD_Precond_Vars      ,ONLY:invBJ,PrecondMethod
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!!---------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!REAL,INTENT(IN)              :: v(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems)
!!---------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!REAL,INTENT(OUT)             :: z(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems)
!!---------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                      :: iElem,i,j,k,iVar,jVar
!!=================================================================================================================================
!
!SELECT CASE(PrecondMethod)
!  CASE(0) ! using a loop
!    ! loop over all DOFS
!    z=0.
!    DO iElem=1,PP_nElems
!      DO k=0,PP_N
!        DO j=0,PP_N
!          DO i=0,PP_N
!            DO jVar=1,PP_nVar
!              DO iVar=1,PP_nVar
!                z(iVar,i,j,k,iElem) = z(iVar,i,j,k,iElem) + invBJ(iVar,jVar,i,j,k,iElem)*v(jVar,i,j,k,iElem)
!              END DO ! iVar
!            END DO ! jVar
!          END DO ! i
!        END DO ! j
!      END DO ! k
!    END DO ! Elem
!  CASE(1) ! Matmul
!    DO iElem=1,PP_nElems
!      DO k=0,PP_N
!        DO j=0,PP_N
!          DO i=0,PP_N
!            z(:,i,j,k,iElem) = MATMUL(invBJ(:,:,i,j,k,iElem),v(:,i,j,k,iElem))
!          END DO ! i
!        END DO ! j
!      END DO ! k
!    END DO !iElem
!END SELECT
!
!END SUBROUTINE  ApplyPrecond_DOF
!
!SUBROUTINE ApplyPrecond_Jacobi(v,z)
!!=================================================================================================================================
!! Apply Jacobi Preconditioner which is constructed per DOF
!!=================================================================================================================================
!! MODULES
!USE MOD_PreProc
!USE MOD_Precond_Vars  ,ONLY:invJ,PrecondMethod
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!!---------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!REAL,INTENT(IN)              :: v(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems)
!!---------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!REAL,INTENT(OUT)             :: z(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems)
!!---------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                      :: iElem,i,j,k,iVar
!!=================================================================================================================================
!
! ! loop over all DOFS
!DO iElem=1,PP_nElems
!  DO k=0,PP_N
!    DO j=0,PP_N
!      DO i=0,PP_N
!        DO iVar=1,PP_nVar
!          z(iVar,i,j,k,iElem) =  invJ(iVar,i,j,k,iElem)*v(iVar,i,j,k,iElem)
!        END DO ! iVar
!      END DO ! i
!    END DO ! j
!  END DO ! k
!END DO ! Elem
!
!END SUBROUTINE  ApplyPrecond_Jacobi
!

!SUBROUTINE ApplyPrecond_GMRES(coeff,Vin,Vout)
!===================================================================================================================================
! Uses matrix free to solve the linear system
! Attention: We use DeltaX=0 as our initial guess   ! why not Un??
!            X0 is allready stored in U
!===================================================================================================================================
! MODULES
!USE MOD_PreProc
!USE MOD_Globals
!USE MOD_JacDG,                ONLY: JacDG
!USE MOD_LinearSolver_Vars,    ONLY: eps_LinearSolver
!USE MOD_LinearSolver_Vars,    ONLY: nRestarts
!USE MOD_LinearOperator,       ONLY: ElementVectorDotProduct
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!REAL,INTENT(IN)    :: coeff, Vin(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!REAL,INTENT(OUT)   :: Vout(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL,ALLOCATABLE,DIMENSION(:,:,:,:)   :: Un, W,Z,R0
!REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:) :: V
!REAL,ALLOCATABLE                        :: Gam(:),C(:),S(:),H(:,:),Alp(:)
!REAL                                    :: Norm_R0,Resu,Temp,Bet
!REAL                                    :: AbortCrit
!INTEGER                                 :: Restart
!INTEGER                                 :: m,nn,o,IterPrecond
!!REAL                                    :: tS,tE
!INTEGER                                 :: nKDimPrecond
!INTEGER                                 :: iElem
!===================================================================================================================================

!nKDimPrecond=3
!
!ALLOCATE( Gam(1:nKDimPrecond+1)           &
!        , C  (1:nKDimPrecond)             &
!        , S  (1:nKDimPrecond)             &
!        , H  (1:nKDimPrecond+1,1:nKDIMPrecond+1) &
!        , Alp(1:nKDimPrecond)             )
!
!ALLOCATE( Un(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)           &
!        , W (1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)           &
!        , Z (1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)           &
!        , R0(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)           &
!        , V (1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nKDimPrecond)   )
!
!! time measurement
!!CALL CPU_TIME(tS)
!DO iElem=1,PP_nElems
!  ! start GMRES
!  Restart=0
!  Un=Vin(:,:,:,:,iElem)
!  IterPrecond=0
!!  DO WHILE (Restart<nRestarts)
!    ! start residuum berrechnen
!    CALL JacDG(coeff,iElem,Un,R0)
!    R0=Un-R0
!    CALL ElementVectorDotProduct(R0,R0,Norm_R0)
!    Norm_R0=SQRT(Norm_R0)
!    !AbortCrit=Norm_R0*eps_LinearSolver
!    AbortCrit=Norm_R0*1e-3
!    ! GMRES(m)  inner loop
!    V(:,:,:,:,1)=R0/Norm_R0
!    Gam(1)=Norm_R0
!    DO m=1,nKDimPrecond
!      IterPrecond=IterPrecond+1
!      ! matrix vector
!      CALL JacDG(coeff,iElem,V(:,:,:,:,m),W)
!      ! Gram-Schmidt
!      DO nn=1,m
!        CALL ElementVectorDotProduct(V(:,:,:,:,nn),W,H(nn,m))
!        W=W-H(nn,m)*V(:,:,:,:,nn)
!      END DO !nn
!      CALL ElementVectorDotProduct(W,W,Resu)
!      H(m+1,m)=SQRT(Resu)
!      ! Givens Rotation
!      DO nn=1,m-1
!        Temp     =   C(nn)*H(nn,m) + S(nn)*H(nn+1,m)
!        H(nn+1,m) = - S(nn)*H(nn,m) + C(nn)*H(nn+1,m)
!        H(nn,m)   =   Temp
!      END DO !nn
!      Bet=SQRT(H(m,m)*H(m,m)+H(m+1,m)*H(m+1,m))
!      S(m)=H(m+1,m)/Bet
!      C(m)=H(m,m)/Bet 
!      H(m,m)=Bet
!      Gam(m+1)=-S(m)*Gam(m)
!      Gam(m)=C(m)*Gam(m)
!      IF ((ABS(Gam(m+1)).LE.AbortCrit) .OR. (m.EQ.nKDimPrecond)) THEN !converge or max Krylov reached
!        DO nn=m,1,-1
!           Alp(nn)=Gam(nn) 
!           DO o=nn+1,m
!             Alp(nn)=Alp(nn) - H(nn,o)*Alp(o)
!           END DO !o
!           Alp(nn)=Alp(nn)/H(nn,nn)
!        END DO !nn
!        DO nn=1,m
!          Un=Un+Alp(nn)*V(:,:,:,:,nn)
!        END DO !nn
!        IF ((ABS(Gam(m+1)).LE.AbortCrit).OR.(M.EQ.nKDimPrecond)) THEN !converged
!          Vout(:,:,:,:,iElem)=Un
!          
!  !        CALL CPU_TIME(tE)
!  !        ! Debug Ausgabe, Anzahl der Iterationen...
!  !        SWRITE(UNIT_stdOut,'(A22,I5)')      ' Iter Precond       : ',IterPrecond
!  !        SWRITE(UNIT_stdOut,'(A22,F16.9)')   ' Time in P-GMRES    : ',tE-tS
!  !        SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R0            : ',Gam(1)
!  !        SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R             : ',Gam(m+1)
!          !SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R/Norm_R0       : ',Gam(m+1)/Gam(1)
!          EXIT
!        END IF  ! converged
!      ELSE ! no convergence, next iteration   ((ABS(Gam(m+1)).LE.AbortCrit) .OR. (m.EQ.nKDim)) 
!        V(:,:,:,:,m+1)=W/H(m+1,m)
!      END IF ! ((ABS(Gam(m+1)).LE.AbortCrit) .OR. (m.EQ.nKDim))
!    END DO ! m 
!    ! Restart needed
!!    Restart=Restart+1
!!  END DO ! Restart
!
!END DO ! iElem
!DEALLOCATE(Un,W,R0,Z,V)
!DEALLOCATE(Gam,C,S,H,Alp)
!
!RETURN
!
!CALL abort(__STAMP__, &
!     'GMRES_M NOT CONVERGED WITH RESTARTS AND GMRES ITERATIONS:',Restart,REAL(IterPrecond))
!
!END SUBROUTINE ApplyPrecond_GMRES

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

!SUBROUTINE Apply_J(BJ,iElem)
!===================================================================================================================================
!! Deallocate global variables
!===================================================================================================================================
!! MODULES
!USE MOD_PreProc
!USE MOD_Implicit_Vars ,ONLY:nDOFElem
!USE MOD_Mesh_Vars     ,ONLY:sJ
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!INTEGER,INTENT(IN) :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!REAL,INTENT(INOUT) :: BJ(nDOFelem,nDOFelem)
!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES 
!INTEGER :: r,s,i,j,k
!===================================================================================================================================
!DO s=0,nDOFelem-1,PP_nVar
!  r=0
!  DO k=0,PP_N
!    DO j=0,PP_N
!      DO i=0,PP_N
!        !BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = -sJ(i,j,k,iElem)*BJ(r+1:r+PP_nVar,s+1:s+PP_nVar)
!        BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar)/sJ(i,j,k,iElem)
!        r=r+PP_nVar
!      END DO !i
!    END DO !j
!  END DO !k
!END DO ! s
!END SUBROUTINE Apply_J
!
!SUBROUTINE Add_J(BJ,iElem)
!===================================================================================================================================
!! Deallocate global variables
!===================================================================================================================================
!! MODULES
!USE MOD_PreProc
!USE MOD_Implicit_Vars ,ONLY:nDOFElem
!USE MOD_Mesh_Vars     ,ONLY:sJ
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!INTEGER,INTENT(IN) :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!REAL,INTENT(INOUT) :: BJ(nDOFelem,nDOFelem)
!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES 
!INTEGER :: r,s,i,j,k,vn1,vn2
!===================================================================================================================================
!
!vn1 = PP_nVar * (PP_N + 1)
!vn2 = vn1 * (PP_N +1)
!
!
!DO k = 0,PP_N
!   DO j = 0,PP_N
!     DO i = 0,PP_N
!       r = vn2 * k + vn1 * j + PP_nVar * i
!       s = r
!       BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) +1./sJ(i,j,k,iElem)
!   END DO ! i
!  END DO ! j
!END DO ! k
!
!END SUBROUTINE Add_J
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

#ifdef DOLDMIST 
!/*do not compile, backup stored */
SUBROUTINE BuildPrecond(t,tStage,tDeriv,alpha,dt)
!===================================================================================================================================
! Build preconditioner for each element, calls a type of preconditioner 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Basis         ,ONLY:GetInverse
USE MOD_LinearSolver_Vars ,ONLY:nDOFelem
USE MOD_Precond_Vars  ,ONLY:invP,PrecondType,DebugMatrix,invBJ,invJ
USE MOD_Jac_ex        ,ONLY:Jac_ex,Jac_ex1D
USE MOD_Jac_FD        ,ONLY:Jac_FD_slow!,Jac_FD
USE MOD_Jac_FD_Vars,   ONLY:reps0,XK
USE MOD_JacDG         ,ONLY:JacBJDOF
!USE MOD_Prepare_FD    ,ONLY:PrepareFD
USE MOD_DG,            ONLY: DGTimeDerivative_WeakForm
USE MOD_DG_Vars,       ONLY: U,Ut
USE MOD_SparseILU,     ONLY: BuildILU0
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
INTEGER            :: r,s,iElem
REAL,DIMENSION(:,:),ALLOCATABLE :: Ploc,Ploc1
REAL,DIMENSION(:,:,:,:,:),ALLOCATABLE :: BJDOF
REAL               :: dummy
REAL               :: TimeStart(3)
REAL               :: TimeEnd(3)
REAL               :: TotalTime(3)
CHARACTER(LEN=255) :: Filename
CHARACTER(LEN=17)  :: strfmt
INTEGER            :: i,j,k,vn1,vn2,mm
REAL               :: coeff
!===================================================================================================================================
IF(PrecondType.EQ.0) RETURN !NO PRECONDITIONER
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' BUILD PRECONDITIONER...'

SELECT CASE(PrecondType)
  CASE(1,2,5,6,7,8,9,10,11)
    ALLOCATE( Ploc (1:nDOFElem,1:nDOFElem), &
              Ploc1(1:nDOFElem,1:nDOFElem) )
  CASE(3)
    ALLOCATE( BJDOF(1:PP_nVar,1:PP_nVar,0:PP_N,0:PP_N,0:PP_N) )
    ALLOCATE( Ploc (1:nDOFElem,1:nDOFElem))
  CASE(22) 
    ALLOCATE( Ploc (1:nDOFElem,1:nDOFElem))
END SELECT

TotalTime=0.

!Prepare Linearisation State 
CALL DGTimeDerivative_WeakForm(t,tStage,tDeriv,doSource=.FALSE.)

DO iElem=1,PP_nElems
  CALL CPU_TIME(TimeStart(1))
!  WRITE(UNIT_stdOut,'(A,I6,A)')'Build Jacobian, iElem=',iElem,' ... '
 
  !compute derivative of DG residual, dRdU  
  SELECT CASE(PrecondType)
    CASE(1)
      DO s=1,nDOFelem
        DO r=1,nDOFelem
          Ploc(r,s)=0.
        END DO !r
      END DO !s
      ! finit differences per Element ! never to use ... 
      CALL Jac_FD_slow(t,tStage,tDeriv,iElem,Ploc)
    CASE(2,3,5,6,7,8,9,10)
      DO s=1,nDOFelem
        DO r=1,nDOFelem
          Ploc(r,s)=0.
        END DO !r
      END DO !s
      ! analytic per Element
      CALL Jac_ex(iElem,Ploc)
!    CASE(3)
!      ! BJ only per DOF
!      coeff = alpha*dt
!      CALL JacBJDOF(coeff,iElem,BJDOF)
!      CALL CPU_TIME(TimeEnd(1))
    CASE(22)
      Ploc=0.
      CALL Jac_ex(iElem,Ploc)
    CASE(60)
      dRdXi=0.
      dRdEta=0.
      dRdZeta=0.
      Jac_ex1D(dRdXi,dRdEta,dRdZeta,ielem)
  END SELECT

  ! compute stuff for element block-jaboci
  SELECT CASE(PrecondType)
    CASE(1,2)
      ! add contibution I-alpha*dt*dRdU
      DO s=1,nDOFelem
        DO r=1,nDOFelem
          Ploc(r,s)=-alpha*dt*Ploc(r,s)
        END DO !r
        Ploc(s,s)=Ploc(s,s)+1.
      END DO !s
      !CALL Add_J(Ploc(:,:),iElem)
      ! testing 3
      CALL CPU_TIME(TimeEnd(1))

      CALL CPU_TIME(TimeStart(2))

      !invert Ploc => invP(:,:,iElem)
      invP(:,:,iElem)=getInverse(nDOFelem,Ploc)
      !CALL Apply_J(invP(:,:,iElem),iElem)
      CALL CPU_TIME(TimeEnd(2))
      TotalTime=TotalTime+(TimeEnd-TimeStart)
      ! Debug subroutine
      IF(DebugMatrix.NE.0) THEN
        CALL CheckBJPrecond(Ploc1,Ploc,invP(:,:,iElem),iElem,TotalTime(3))
      END IF ! DebugMatrix
    CASE(3)
      DO s=1,nDOFelem
        DO r=1,nDOFelem
          Ploc(r,s)=-alpha*dt*Ploc(r,s)
        END DO !r
        Ploc(s,s)=Ploc(s,s)+1.
      END DO !s
      CALL CPU_TIME(TimeEnd(1))
      CALL CPU_TIME(TimeStart(2))
      vn1=PP_nVar*(PP_N+1)
      vn2=vn1*(PP_N+1)
      ! alpha*dt and the calcution with the unity matrix is already performed
      !CALL CPU_TIME(TimeStart(2))
      DO k=0,PP_N
        DO j=0,PP_N
          DO i=0,PP_N
            !invert Ploc => invP(:,:,iElem)
            !print*,BJDOF(:,:,i,j,k)
            !read*
            s = vn2 * k + vn1 * j + PP_nVar * i
            !CALL BUILDINVERSE(PP_nVar,BJDOF(:,:,i,j,k),invBJ(:,:,i,j,k,iElem))
            !CALL BUILDINVERSE(PP_nVar,Ploc(s+1:s+PP_nVar,s+1:s+PP_nVar),invBJ(:,:,i,j,k,iElem))
            invBJ(:,:,i,j,k,iElem)=getInverse(PP_nVar,Ploc(s+1:s+PP_nVar))
          END DO ! i
        END DO ! j
      END DO ! k
      CALL CPU_TIME(TimeEnd(2))
      TotalTime=TotalTime+(TimeEnd-TimeStart)
      TotalTime(3) = 0. 
      ! test
      !CHECK INVERSE
      DO k=0,PP_N
        DO j=0,PP_N
          DO i=0,PP_N
            dummy=0.
            s = vn2 * k + vn1 * j + PP_nVar * i
            dummy=dummy+SUM(ABS(MATMUL(invBJ(:,:,i,j,k,iElem),Ploc(s+1:s+PP_nVar,s+1:s+PP_nVar))))
            dummy=(dummy-REAL(PP_nVar))/REAL(PP_nVar*PP_nVar) !relative error per matrix entry
            IF(dummy.GT. 1.0E-06) THEN
              IPWRITE(UNIT_stdOut,*)'WARNING!!! accuracy problems in with preconditioner inverse..',dummy
            END IF
            IF(dummy.NE. dummy) THEN  !NAN
              IPWRITE(UNIT_stdOut,*)'WARNING!!! NAN problem in with preconditioner inverse..',dummy
              STOP
            END IF
          END DO
        END DO
      END DO    
    CASE(5)
      ! Jacobi Preconditioner
      DO s=1,nDOFelem
        DO r=1,nDOFelem
          Ploc(r,s)=-alpha*dt*Ploc(r,s)
        END DO !r
        Ploc(s,s)=Ploc(s,s)+1.
      END DO !s
      CALL CPU_TIME(TimeEnd(1))
      CALL CPU_TIME(TimeStart(2))
      vn1=PP_nVar*(PP_N+1)
      vn2=vn1*(PP_N+1)
      ! alpha*dt and the calcution with the unity matrix is already performed
      !CALL CPU_TIME(TimeStart(2))
      DO k=0,PP_N
        DO j=0,PP_N
          DO i=0,PP_N
            DO r = 1,PP_nVar
              !read*
               s = vn2 * k + vn1 * j + PP_nVar * i
              invJ(r,i,j,k,iElem) = 1.0/Ploc(s+r,s+r)
            END DO ! r
          END DO ! i
        END DO ! j
      END DO ! k
      CALL CPU_TIME(TimeEnd(2))
      TotalTime=TotalTime+(TimeEnd-TimeStart)
      TotalTime(3) = 0. 
    CASE(6)
      ! Jacobi Preconditioner
      ! Spaltensumme
      DO s=1,nDOFelem
        DO r=1,nDOFelem
          Ploc(r,s)=-alpha*dt*Ploc(r,s)
        END DO !r
        Ploc(s,s)=Ploc(s,s)+1.
      END DO !s
      CALL CPU_TIME(TimeEnd(1))
      CALL CPU_TIME(TimeStart(2))
      vn1=PP_nVar*(PP_N+1)
      vn2=vn1*(PP_N+1)
      ! alpha*dt and the calcution with the unity matrix is already performed
      !CALL CPU_TIME(TimeStart(2))
      DO k=0,PP_N
        DO j=0,PP_N
          DO i=0,PP_N
            DO r = 1,PP_nVar
              s = vn2 * k + vn1 * j + PP_nVar * i
              dummy=SUM(ABS(Ploc(:,s+r)))
              invJ(r,i,j,k,iElem) = 1.0/dummy
            END DO ! r
          END DO ! i
        END DO ! j
      END DO ! k
      CALL CPU_TIME(TimeEnd(2))
      TotalTime=TotalTime+(TimeEnd-TimeStart)
      TotalTime(3) = 0. 
    CASE(7)
      ! Jacobi Preconditioner
      ! Zeilensumme
      DO s=1,nDOFelem
        DO r=1,nDOFelem
          Ploc(r,s)=-alpha*dt*Ploc(r,s)
        END DO !r
        Ploc(s,s)=Ploc(s,s)+1.
      END DO !s
      CALL CPU_TIME(TimeEnd(1))
      CALL CPU_TIME(TimeStart(2))
      vn1=PP_nVar*(PP_N+1)
      vn2=vn1*(PP_N+1)
      ! alpha*dt and the calcution with the unity matrix is already performed
      !CALL CPU_TIME(TimeStart(2))
      DO k=0,PP_N
        DO j=0,PP_N
          DO i=0,PP_N
            DO r = 1,PP_nVar
              s = vn2 * k + vn1 * j + PP_nVar * i
              dummy=SUM(ABS(Ploc(s+r,:)))
              invJ(r,i,j,k,iElem) = 1.0/dummy
            END DO ! r
          END DO ! i
        END DO ! j
      END DO ! k
      CALL CPU_TIME(TimeEnd(2))
      TotalTime=TotalTime+(TimeEnd-TimeStart)
      TotalTime(3) = 0. 
    CASE(8)
      ! Jacobi Preconditioner
      ! Spalt Euklid
      DO s=1,nDOFelem
        DO r=1,nDOFelem
          Ploc(r,s)=-alpha*dt*Ploc(r,s)
        END DO !r
        Ploc(s,s)=Ploc(s,s)+1.
      END DO !s
      CALL CPU_TIME(TimeEnd(1))
      CALL CPU_TIME(TimeStart(2))
      vn1=PP_nVar*(PP_N+1)
      vn2=vn1*(PP_N+1)
      ! alpha*dt and the calcution with the unity matrix is already performed
      !CALL CPU_TIME(TimeStart(2))
      DO k=0,PP_N
        DO j=0,PP_N
          DO i=0,PP_N
            DO r = 1,PP_nVar
              s = vn2 * k + vn1 * j + PP_nVar * i
              dummy=SUM(Ploc(:,s+r)*Ploc(:,s+r))
              invJ(r,i,j,k,iElem) = 1.0/SQRT(dummy)
            END DO ! r
          END DO ! i
        END DO ! j
      END DO ! k
      CALL CPU_TIME(TimeEnd(2))
      TotalTime=TotalTime+(TimeEnd-TimeStart)
      TotalTime(3) = 0. 
    CASE(9)
      ! Jacobi Preconditioner
      ! Zeil Euklid
      DO s=1,nDOFelem
        DO r=1,nDOFelem
          Ploc(r,s)=-alpha*dt*Ploc(r,s)
        END DO !r
        Ploc(s,s)=Ploc(s,s)+1.
      END DO !s
      CALL CPU_TIME(TimeEnd(1))
      CALL CPU_TIME(TimeStart(2))
      vn1=PP_nVar*(PP_N+1)
      vn2=vn1*(PP_N+1)
      ! alpha*dt and the calcution with the unity matrix is already performed
      !CALL CPU_TIME(TimeStart(2))
      DO k=0,PP_N
        DO j=0,PP_N
          DO i=0,PP_N
            DO r = 1,PP_nVar
              s = vn2 * k + vn1 * j + PP_nVar * i
              dummy=SUM(Ploc(s+r,:)*Ploc(s+r,:))
              invJ(r,i,j,k,iElem) = 1.0/SQRT(dummy)
            END DO ! r
          END DO ! i
        END DO ! j
      END DO ! k
      CALL CPU_TIME(TimeEnd(2))
      TotalTime=TotalTime+(TimeEnd-TimeStart)
      TotalTime(3) = 0. 
    CASE(10)
      ! Jacobi Preconditioner
      ! Spalte Maximum 
      DO s=1,nDOFelem
        DO r=1,nDOFelem
          Ploc(r,s)=-alpha*dt*Ploc(r,s)
        END DO !r
        Ploc(s,s)=Ploc(s,s)+1.
      END DO !s
      CALL CPU_TIME(TimeEnd(1))
      CALL CPU_TIME(TimeStart(2))
      vn1=PP_nVar*(PP_N+1)
      vn2=vn1*(PP_N+1)
      ! alpha*dt and the calcution with the unity matrix is already performed
      !CALL CPU_TIME(TimeStart(2))
      DO k=0,PP_N
        DO j=0,PP_N
          DO i=0,PP_N
            DO r = 1,PP_nVar
              s = vn2 * k + vn1 * j + PP_nVar * i
              dummy=ABS(Ploc(1,s+r))
              DO mm = 2,nDOFElem
                dummy=MAX(ABS(Ploc(mm,s+r)),dummy)
              END DO
              invJ(r,i,j,k,iElem) = 1.0/dummy
            END DO ! r
          END DO ! i
        END DO ! j
      END DO ! k
      CALL CPU_TIME(TimeEnd(2))
      TotalTime=TotalTime+(TimeEnd-TimeStart)
      TotalTime(3) = 0. 
    CASE(11)
      ! Jacobi Preconditioner
      ! Zeile Maximum
      DO s=1,nDOFelem
        DO r=1,nDOFelem
          Ploc(r,s)=-alpha*dt*Ploc(r,s)
        END DO !r
        Ploc(s,s)=Ploc(s,s)+1.
      END DO !s
      CALL CPU_TIME(TimeEnd(1))
      CALL CPU_TIME(TimeStart(2))
      vn1=PP_nVar*(PP_N+1)
      vn2=vn1*(PP_N+1)
      ! alpha*dt and the calcution with the unity matrix is already performed
      !CALL CPU_TIME(TimeStart(2))
      DO k=0,PP_N
        DO j=0,PP_N
          DO i=0,PP_N
            DO r = 1,PP_nVar
              s = vn2 * k + vn1 * j + PP_nVar * i
              dummy=ABS(Ploc(s+r,1))
              DO mm = 2,nDOFElem
                dummy=MAX(ABS(Ploc(s+r,mm)),dummy)
              END DO
              invJ(r,i,j,k,iElem) = 1.0/dummy
            END DO ! r
          END DO ! i
        END DO ! j
      END DO ! k
      CALL CPU_TIME(TimeEnd(2))
      TotalTime=TotalTime+(TimeEnd-TimeStart)
      TotalTime(3) = 0. 
    CASE(22)
      ! add contibution I-alpha*dt*dRdU
      DO s=1,nDOFelem
        DO r=1,nDOFelem
          Ploc(r,s)=-alpha*dt*Ploc(r,s)
        END DO !r
        Ploc(s,s)=Ploc(s,s)+1.
      END DO !s
      ! testing 3
      CALL CPU_TIME(TimeEnd(1))

      CALL CPU_TIME(TimeStart(2))
      CALL BuildILU0(Ploc,iElem)
      CALL CPU_TIME(TimeEnd(2))
      TotalTime=TotalTime+(TimeEnd-TimeStart)
    END SELECT
END DO !iElem

WRITE(UNIT_stdOut,'(A,F11.3,A)')' TOTAL BUILDING  TIME =[',TotalTime(1),' ]'
WRITE(UNIT_stdOut,'(A,F11.3,A)')' TOTAL INVERTING TIME =[',TotalTime(2),' ]'
IF(debugMatrix.GT.2) WRITE(UNIT_stdOut,'(A,F11.3,A)')' TOTAL CHECKING  TIME =[',TotalTime(3),' ]'
SWRITE(UNIT_stdOut,'(A)')' BUILD PRECONDITIONER DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

SELECT CASE(PrecondType)
  CASE(1,2) 
    DEALLOCATE( Ploc, Ploc1)
  CASE(3)
    DEALLOCATE(BJDOF)
  CASE(22) 
    DEALLOCATE( Ploc )
END SELECT

END SUBROUTINE  BuildPrecond
#endif /* do not compile, backup */

END MODULE MOD_Precond
