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
USE MOD_ReadInTools,       ONLY:GETINT,GETLOGICAL
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
UpdatePrecond = GETLOGICAL('UpdatePrecond','.FALSE.')

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

  CASE(201)
    ! 1D-Tensor BJ preconditioner
    locSize=PP_nVar*(PP_N+1)
    ALLOCATE( invXi  (locSize,locSize,0:PP_N,0:PP_N,1:PP_nElems) &
            , invEta (locSize,locSize,0:PP_N,0:PP_N,1:PP_nElems) &
            , invZeta(locSize,locSize,0:PP_N,0:PP_N,1:PP_nElems) )
    ALLOCATE( dRdXi  (locSize,locSize,0:PP_N,0:PP_N,1:PP_nElems) &
            , dRdEta (locSize,locSize,0:PP_N,0:PP_N,1:PP_nElems) &
            , dRdZeta(locSize,locSize,0:PP_N,0:PP_N,1:PP_nElems) )
    CALL InitJac_Ex()

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
USE MOD_Mathtools         ,ONLY: INVERSE
USE MOD_Precond_Vars      ,ONLY: invXi,invEta,invZeta,dRdXi,dRdZeta,dRdEta
USE MOD_LinearSolver_Vars ,ONLY: nDOFelem,mass,nDOFLine
USE MOD_Precond_Vars      ,ONLY: invP,PrecondType,DebugMatrix
USE MOD_Jac_ex            ,ONLY: Jac_ex, Jac_Ex_Neighbor,Jac_ex1D
USE MOD_Jac_FD            ,ONLY: Jac_FD_slow
USE MOD_JacDG             ,ONLY: BuildJacDG
USE MOD_DG                ,ONLY: DGTimeDerivative_WeakForm
USE MOD_SparseILU         ,ONLY: BuildILU0
USE MOD_ILU               ,ONLY: BuildBILU0BCSR
USE MOD_CSR               ,ONLY: CSR
#if USE_MPI
USE MOD_MPI_Vars
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars  ,ONLY: PerformLBSample
USE MOD_LoadBalance_Vars  ,ONLY: ElemTime,ElemTimeField
#endif /*USE_LOADBALANCE*/
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
INTEGER            :: i,j,k,p,q
INTEGER            :: oo,mm,nn,ll,v1,v2
REAL               :: coeff
#if USE_MPI
REAL               :: TotalTimeMPI(3)
#endif /*USE_MPI*/
#if USE_LOADBALANCE
REAL               :: ElemTimePrecond
#endif /*USE_LOADBALANCE*/
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
  SELECT CASE(PrecondType)
  CASE DEFAULT ! normal, block preconditioner
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
#ifdef IMPA
    coeff=-alpha*dt
    DO s=1,nDOFelem
      DO r=1,nDOFelem
        Ploc(r,s)=coeff*Ploc(r,s)
      END DO !r
      Ploc(s,s)=Ploc(s,s)+1.
    END DO !s
#else
    coeff=alpha*dt
    DO s=1,nDOFelem
      DO r=1,nDOFelem
        Ploc(r,s)=-Ploc(r,s)
      END DO !r
      Ploc(s,s)=Ploc(s,s)+coeff
    END DO !s
#endif
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
  CASE(201) ! 1D-Tensor-Product preconditioner
    ! build 1D Tensor line derivative
    dRdXi  (:,:,:,:,iElem)=0.
    dRdEta (:,:,:,:,iElem)=0.
    dRdZeta(:,:,:,:,iElem)=0.
    ! derivative for 1D
    CALL Jac_ex1D(dRdXi  (:,:,:,:,iElem) &
                 ,dRdEta (:,:,:,:,iElem) &
                 ,dRdZeta(:,:,:,:,iElem) ,iElem )
    ! apply coefficient
#ifdef IMPA
    coeff=-alpha*dt
    DO q=0,PP_N
      DO p=0,PP_N
        DO s=1,nDOFLine
          DO r=1,nDOFLine
            dRdXi  (r,s,p,q,iElem)=coeff*dRdXi  (r,s,p,q,iElem)
            dRdEta (r,s,p,q,iElem)=coeff*dRdEta (r,s,p,q,iElem)
            dRdZeta(r,s,p,q,iElem)=coeff*dRdZeta(r,s,p,q,iElem)
          END DO ! r
          dRdXi  (s,s,p,q,iElem)=dRdXi  (s,s,p,q,iElem) +1.
          dRdEta (s,s,p,q,iElem)=dRdEta (s,s,p,q,iElem) +1.
          dRdZeta(s,s,p,q,iElem)=dRdZeta(s,s,p,q,iElem) +1.
        END DO ! s
      END DO ! p
    END DO !q
#else
    ! change sign and add 1/(gamma_ii*dt)
    coeff=alpha*dt
    DO q=0,PP_N
      DO p=0,PP_N
        DO s=1,nDOFLine
          DO r=1,nDOFLine
            dRdXi  (r,s,p,q,iElem)=-dRdXi  (r,s,p,q,iElem)
            dRdEta (r,s,p,q,iElem)=-dRdEta (r,s,p,q,iElem)
            dRdZeta(r,s,p,q,iElem)=-dRdZeta(r,s,p,q,iElem)
          END DO ! r
          dRdXi  (s,s,p,q,iElem)=dRdXi  (s,s,p,q,iElem) +coeff
          dRdEta (s,s,p,q,iElem)=dRdEta (s,s,p,q,iElem) +coeff
          dRdZeta(s,s,p,q,iElem)=dRdZeta(s,s,p,q,iElem) +coeff
        END DO ! s
      END DO ! p
    END DO !q
#endif
    DO oo=0,PP_N; DO nn=0,PP_N; DO mm=0,PP_N
      v1=0
      DO ll=0,PP_N
        v2=mm*PP_nVar
        dRdXi  (v1+1:v1+PP_nVar,v2+1:v2+PP_nVar,nn,oo,iElem)=dRdXi  (v1+1:v1+PP_nVar,v2+1:v2+PP_nVar,nn,oo,iElem) &
                                                              *Mass(1,ll,nn,oo,iElem)
        v2=nn*PP_nVar
        dRdEta (v1+1:v1+PP_nVar,v2+1:v2+PP_nVar,mm,oo,iElem)=dRdEta (v1+1:v1+PP_nVar,v2+1:v2+PP_nVar,mm,oo,iElem) &
                                                              *Mass(1,mm,ll,oo,iElem)
        v2=oo*PP_nVar
        dRdZeta(v1+1:v1+PP_nVar,v2+1:v2+PP_nVar,mm,nn,iElem)=dRdZeta(v1+1:v1+PP_nVar,v2+1:v2+PP_nVar,mm,nn,iElem) &
                                                            *Mass(1,mm,nn,ll,iElem)
        v1=v1+PP_nVar
      END DO !ll
    END DO; END DO; END DO! mm,nn,oo=0,PP_N
  END SELECT
  CALL CPU_TIME(TimeEnd(1))
  CALL CPU_TIME(TimeStart(2))
  SELECT CASE(PrecondType)
  CASE(1,2) ! element block jacobi
    !invert Ploc => invP(:,:,iElem)
    invP(:,:,iElem)=INVERSE(Ploc)
  CASE(3) ! ilu(0) of element
    CALL BuildILU0(Ploc,iElem)
  CASE(4)
    CALL BuildBILU0BCSR(Ploc,iElem)
  CASE(201)
    ! compute the inverse of the 1D preconditioner
    DO q=0,PP_N
      DO p=0,PP_N
        invXi  (:,:,p,q,iElem)=INVERSE(dRdXi  (:,:,p,q,iElem))
        invEta (:,:,p,q,iElem)=INVERSE(dRdEta (:,:,p,q,iElem))
        invZeta(:,:,p,q,iElem)=INVERSE(dRdZeta(:,:,p,q,iElem))
      END DO ! j
    END DO ! k
  END SELECT
  CALL CPU_TIME(TimeEnd(2))
  TotalTime=TotalTime+(TimeEnd-TimeStart)
#if USE_LOADBALANCE
  IF(PerformLBSample)THEN
    ElemTimePrecond = SUM(TimeEnd(1:2)-TimeStart(1:2))
    ElemTimeField   = ElemTimeField  +ElemTimePrecond
    ElemTime(iElem) = ElemTime(iElem)+ElemTimePrecond
  END IF
#endif /*USE_LOADBALANCE*/
  IF((PrecondType.LE.2).AND.DebugMatrix.NE.0)THEN
    CALL CheckBJPrecond(Ploc1,Ploc,invP(:,:,iElem),iElem,TotalTime(3))
  END IF ! DebugMatrix
END DO ! ! iELEM

#if USE_MPI
CALL MPI_REDUCE(TotalTime,TotalTimeMPI ,3,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_PICLAS,IERROR)
IF(MPIRoot) THEN
  TotalTime=TotalTimeMPI
END IF

#endif /*USE_MPI*/

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
#if USE_MPI
USE MOD_MPI_Vars
#endif /*USE_MPI*/
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
#if USE_MPI
  WRITE(Filename,'(A,I2.2,A,I2.2,A,I4.4,A)')'BlockJac_Rank_',MyRank,'_Type_',PreCondType,'_Mat_', iElem,'.dat'
#else
  WRITE(Filename,'(A,I2.2,A,I4.4,A)')'BlockJac_',PreCondType,'_Mat_', iElem,'.dat'
#endif /*USE_MPI*/
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

#if USE_MPI
  WRITE(Filename,'(A,I2.2,A,I2.2,A,I4.4,A)')'Precond_Rank_',MyRank,'_Type_',PreCondType,'_Mat_', iElem,'.dat'
#else
  WRITE(Filename,'(A,I2.2,A,I4.4,A)')'Precond_',PreCondType,'_Mat_', iElem,'.dat'
#endif /*USE_MPI*/
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
#if USE_MPI
  WRITE(Filename,'(A,I2.2,A,I2.2,A,I4.4,A)')'Precond_Rank_',MyRank,'_Type__',PreCondType,'_InvMat_', iElem,'.dat'
#else
  WRITE(Filename,'(A,I2.2,A,I4.4,A)')'Precond_',PreCondType,'_InvMat_', iElem,'.dat'
#endif /*USE_MPI*/
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
USE MOD_Precond_Vars,ONLY:invP,PrecondInitIsDone,PrecondType,neighborElemID
USE MOD_Precond_Vars,ONLY:invXi,invEta,invZeta,dRdXi,dRdZeta,dRdEta
USE MOD_SparseILU   ,ONLY:FinalizeSparseILU
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(invP)
SDEALLOCATE(neighborElemID)
SDEALLOCATE(invXi)
SDEALLOCATE(invEta)
SDEALLOCATE(invZeta)
SDEALLOCATE(dRdXi)
SDEALLOCATE(dRdEta)
SDEALLOCATE(dRdZeta)
PrecondInitIsDone = .FALSE.
IF(PrecondType.EQ.3) CALL FinalizeSparseILU

END SUBROUTINE FinalizePrecond


END MODULE MOD_Precond
