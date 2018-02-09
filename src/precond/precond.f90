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

PUBLIC :: InitPrecond
#ifdef maxwell
PUBLIC :: BuildPrecond
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
    CALL abort(&
        __STAMP__&
        ,'Preconditioner not implemented!',999,999.)
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
USE MOD_LinearSolver_Vars ,ONLY:nDOFelem,mass
USE MOD_Precond_Vars      ,ONLY:invP,PrecondType,DebugMatrix
USE MOD_Jac_ex            ,ONLY:Jac_ex, Jac_Ex_Neighbor,Jac_ex1D
USE MOD_Jac_FD            ,ONLY:Jac_FD_slow!,Jac_FD
USE MOD_JacDG             ,ONLY:BuildJacDG
USE MOD_DG,                ONLY: DGTimeDerivative_WeakForm
USE MOD_SparseILU,         ONLY: BuildILU0
USE MOD_ILU,               ONLY: BuildBILU0BCSR
USE MOD_CSR,               ONLY: CSR
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
INTEGER            :: r,s,iElem
REAL,DIMENSION(:,:),ALLOCATABLE :: Ploc,Ploc1!,Ploc2
REAL               :: TimeStart(3)
REAL               :: TimeEnd(3)
REAL               :: TotalTime(3)
REAL               :: delta(0:PP_N,0:PP_N)
INTEGER            :: i,j,k
REAL               :: coeff
#ifdef MPI
REAL               :: TotalTimeMPI(3)
#endif /*MPI*/
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
TimeStart=0.
TimeEnd=0.

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

#ifdef MPI
CALL MPI_REDUCE(TotalTime,TotalTimeMPI ,3,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,IERROR)
IF(MPIRoot) THEN
  TotalTime=TotalTimeMPI
END IF

#endif /*MPI*/

SWRITE(UNIT_stdOut,'(A,F11.3,A)')' TOTAL DERIVATING TIME =[',TotalTime(1),' ]'
SWRITE(UNIT_stdOut,'(A,F11.3,A)')' TOTAL INVERTING  TIME =[',TotalTime(2),' ]'
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
    CALL abort(&
        __STAMP__&
        ,' NAN! Problem with preconditioner inverse.')
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
