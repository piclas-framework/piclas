#include "boltzplatz.h"


MODULE MOD_PML
!===================================================================================================================================
!  
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitPML
  MODULE PROCEDURE InitPML
END INTERFACE
INTERFACE FinalizePML
  MODULE PROCEDURE FinalizePML
END INTERFACE
INTERFACE CalcPMLSource
  MODULE PROCEDURE CalcPMLSource
END INTERFACE
INTERFACE TransformPMLVars
  MODULE PROCEDURE TransformPMLVars
END INTERFACE
INTERFACE PMLTimeDerivative
  MODULE PROCEDURE PMLTimeDerivative
END INTERFACE
INTERFACE BacktransformPMLVars
  MODULE PROCEDURE BacktransformPMLVars
END INTERFACE
INTERFACE ProbePML
  MODULE PROCEDURE ProbePML
END INTERFACE
PUBLIC::InitPML,FinalizePML,CalcPMLSource,TransformPMLVars,PMLTimeDerivative,BacktransformPMLVars,ProbePML
!===================================================================================================================================
CONTAINS

SUBROUTINE InitPML()
!===================================================================================================================================
!  Initialize perfectly matched layer
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ReadInTools
USE MOD_PML_Vars,            ONLY: PMLzeta,U2,U2t,Probes,DoPML,ntotalPML
USE MOD_PML_Vars,            ONLY: nPMLElems,ElemtoPML,PMLtoElem
USE MOD_PML_Vars,            ONLY: PMLzeta0,xyzPhysicalMinMax,PMLzetaShape,PMLspread,PMLwritezeta, PMLzetaNorm
USE MOD_Mesh_Vars,           ONLY: Elem_xGP,Face_xGP,nBCSides  ! for PML region: xyz position of the Gauss points and Face Gauss points
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iElem,iProbe,iPMLElem
REAL                :: xyzMinMax(6),xi,L,XiN
REAL                :: zetaVec,zetaVecABS
LOGICAL,ALLOCATABLE :: isPMLElem(:)
INTEGER             :: PMLID,nGlobalPMLElems
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PML...'

!===================================================================================================================================
! Readin
!===================================================================================================================================


! get information of PML size
PMLzeta0               = GETREAL('PMLzeta0','0.')
xyzPhysicalMinMax(1:6) = GETREALARRAY('xyzPhysicalMinMax',6,'0.0,0.0,0.0,0.0,0.0,0.0')
PMLzetaShape           = GETINT('PMLzetaShape','0')
PMLspread              = GETINT('PMLspread','0')
PMLwriteZeta           = GETINT('PMLwriteZeta','0')
PMLzetaNorm            = GETLOGICAL('PMLzetaNorm','.FALSE.')
! only for Maxwell, PP_nVar=8
#if PP_nVar == 8
  DoPML                = GETLOGICAL('DoPML','.FALSE.')
#else
  DoPML=.FALSE.
#endif

IF(.NOT.DoPML) THEN
  SWRITE(UNIT_stdOut,'(A)') ' No PML region detected. '
#if PP_nVar == 1
  SWRITE(UNIT_stdOut,'(A)') ' Equation system does not support a PML '
#endif
#if PP_nVar == 4
  SWRITE(UNIT_stdOut,'(A)') ' Equation system electrostatic does not support a PML '
#endif
  nPMLElems=0
  RETURN
END IF

ALLOCATE(Probes%Distance(3,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
Probes%Distance=HUGE(1.)! Set Probes%Distance to a large number (lowest number counts later on)
!print *, "xyzPhysicalMinMax",xyzPhysicalMinMax
!print *, "PMLzeta size & shape",SIZE(PMLzeta),shape(PMLzeta)
!print *, "U2 size & shape  ",SIZE(U2),shape(U2)
!print *, "PMLzeta0",PMLzeta0

!===================================================================================================================================
! check if Element is PMLElem
!===================================================================================================================================

ALLOCATE(isPMLElem(1:PP_nElems))
isPMLElem=.FALSE.

DO iElem=1,PP_nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  ! x-PML region
  IF (Elem_xGP(1,i,j,k,iElem) .LT. xyzPhysicalMinMax(1) .OR. Elem_xGP(1,i,j,k,iElem) .GT. xyzPhysicalMinMax(2)) THEN        
    isPMLElem(iElem) = .TRUE.
  END IF
  ! y-PML region
  IF (Elem_xGP(2,i,j,k,iElem) .LT. xyzPhysicalMinMax(3) .OR. Elem_xGP(2,i,j,k,iElem) .GT. xyzPhysicalMinMax(4)) THEN        
    isPMLElem(iElem) = .TRUE.
  END IF
  ! z-PML region
  IF (Elem_xGP(3,i,j,k,iElem) .LT. xyzPhysicalMinMax(5) .OR. Elem_xGP(3,i,j,k,iElem) .GT. xyzPhysicalMinMax(6)) THEN        
    isPMLElem(iElem) = .TRUE.
  END IF
END DO; END DO; END DO; END DO !iElem,k,i,j

! Get number of PML Elems
nPMLElems = 0

DO iElem=1,PP_nElems
  IF(isPMLElem(iElem))THEN
    nPMLElems=nPMLElems+1
  END IF
END DO ! iElem

#ifdef MPI
  CALL MPI_REDUCE(nPMLElems,nGlobalPMLElems,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
#else
  nGlobalPMLElems=nPMLElems
#endif /*MPI*/
SWRITE(UNIT_stdOut,'(A,I10,A)') ' Found ', nGlobalPMLElems,' Elems inside of PML region.'

ALLOCATE(ElemToPML(PP_nElems)&
        ,PMLtoElem(nPMLElems))

ElemToPML=0
PMLtoElem=0

! Create array with mapping
PMLID=0
DO iElem=1,PP_nElems
  IF(isPMLElem(iElem))THEN
    PMLID=PMLID+1
    ElemToPML(iElem) = PMLID
    PMLToElem(PMLID) = iElem
  END IF
END DO

DEALLOCATE(isPMLElem)
ALLOCATE(PMLzeta(1:3,0:PP_N,0:PP_N,0:PP_N,1:nPMLElems))
ALLOCATE(U2     (1:6,0:PP_N,0:PP_N,0:PP_N,1:nPMLElems))        
ALLOCATE(U2t    (1:6,0:PP_N,0:PP_N,0:PP_N,1:nPMLElems))
nTotalPML=6*(PP_N+1)**3
! zero
PMLzeta=0.
U2 =0.0
U2t=0.0

!===================================================================================================================================
!determine PMLzeta values for each interpolation point
!===================================================================================================================================
xyzMinMax(:) = (/MINVAL(Face_xGP(1,:,:,1:nBCSides)),MAXVAL(Face_xGP(1,:,:,1:nBCSides)),MINVAL(Face_xGP(2,:,:,1:nBCSides)),&
                 MAXVAL(Face_xGP(2,:,:,1:nBCSides)),MINVAL(Face_xGP(3,:,:,1:nBCSides)),MAXVAL(Face_xGP(3,:,:,1:nBCSides))/)
!print *, "xyzMinMax",xyzMinMax
SELECT CASE (PMLzetaShape)
CASE(0) ! Constant Distribution of the Damping Coefficient
  DO iElem=1,PP_nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        ! x-PML region
        IF (Elem_xGP(1,i,j,k,iElem) .LT. xyzPhysicalMinMax(1) .OR. Elem_xGP(1,i,j,k,iElem) .GT. xyzPhysicalMinMax(2)) THEN        
          PMLzeta(1,i,j,k,ElemToPML(iElem)) = PMLzeta0
        END IF
        ! y-PML region
        IF (Elem_xGP(2,i,j,k,iElem) .LT. xyzPhysicalMinMax(3) .OR. Elem_xGP(2,i,j,k,iElem) .GT. xyzPhysicalMinMax(4)) THEN        
          PMLzeta(2,i,j,k,ElemToPML(iElem)) = PMLzeta0
        END IF
        ! z-PML region
        IF (Elem_xGP(3,i,j,k,iElem) .LT. xyzPhysicalMinMax(5) .OR. Elem_xGP(3,i,j,k,iElem) .GT. xyzPhysicalMinMax(6)) THEN        
          PMLzeta(3,i,j,k,ElemToPML(iElem)) = PMLzeta0
        END IF
  END DO; END DO; END DO; END DO !iElem,k,i,j
CASE(1) ! Linear Distribution of the Damping Coefficient
  DO iElem=1,PP_nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        ! x-PML region
        IF     (Elem_xGP(1,i,j,k,iElem) .LT. xyzPhysicalMinMax(1)) THEN        
          PMLzeta(1,i,j,k,ElemToPML(iElem)) = PMLzeta0*(ABS(Elem_xGP(1,i,j,k,iElem))-ABS(xyzPhysicalMinMax(1)))/&
                                      (ABS(xyzMinMax(1))-ABS(xyzPhysicalMinMax(1)))
        ELSEIF (Elem_xGP(1,i,j,k,iElem) .GT. xyzPhysicalMinMax(2)) THEN
          PMLzeta(1,i,j,k,ElemToPML(iElem)) = PMLzeta0*(ABS(Elem_xGP(1,i,j,k,iElem))-ABS(xyzPhysicalMinMax(2)))/&
                                      (ABS(xyzMinMax(2))-ABS(xyzPhysicalMinMax(2)))
        END IF
        ! y-PML region
        IF     (Elem_xGP(2,i,j,k,iElem) .LT. xyzPhysicalMinMax(3)) THEN        
          PMLzeta(2,i,j,k,ElemToPML(iElem)) = PMLzeta0*(ABS(Elem_xGP(2,i,j,k,iElem))-ABS(xyzPhysicalMinMax(3)))/&
                                      (ABS(xyzMinMax(3))-ABS(xyzPhysicalMinMax(3)))
        ELSEIF (Elem_xGP(2,i,j,k,iElem) .GT. xyzPhysicalMinMax(4)) THEN    
          PMLzeta(2,i,j,k,ElemToPML(iElem)) = PMLzeta0*(ABS(Elem_xGP(2,i,j,k,iElem))-ABS(xyzPhysicalMinMax(4)))/&
                                      (ABS(xyzMinMax(4))-ABS(xyzPhysicalMinMax(4)))
        END IF
        ! z-PML region
        IF     (Elem_xGP(3,i,j,k,iElem) .LT. xyzPhysicalMinMax(5)) THEN        
          PMLzeta(3,i,j,k,ElemToPML(iElem)) = PMLzeta0*(ABS(Elem_xGP(3,i,j,k,iElem))-ABS(xyzPhysicalMinMax(5)))/&
                                      (ABS(xyzMinMax(5))-ABS(xyzPhysicalMinMax(5)))
        ELSEIF (Elem_xGP(3,i,j,k,iElem) .GT. xyzPhysicalMinMax(6)) THEN    
          PMLzeta(3,i,j,k,ElemToPML(iElem)) = PMLzeta0*(ABS(Elem_xGP(3,i,j,k,iElem))-ABS(xyzPhysicalMinMax(6)))/&
                                      (ABS(xyzMinMax(6))-ABS(xyzPhysicalMinMax(6)))
        END IF
  END DO; END DO; END DO; END DO !iElem,k,i,j
CASE(2) ! Sinusoidal  Distribution of the Damping Coefficient
  DO iElem=1,PP_nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        ! x-PML region
        IF     (Elem_xGP(1,i,j,k,iElem) .LT. xyzPhysicalMinMax(1)) THEN
          xi                  = ABS(Elem_xGP(1,i,j,k,iElem))-ABS(xyzPhysicalMinMax(1))
          L                   = ABS(xyzMinMax(1))-ABS(xyzPhysicalMinMax(1))
          PMLzeta(1,i,j,k,ElemToPML(iElem)) = PMLzeta0*(xi/L-SIN(2*ACOS(-1.)*xi/L)/(2*ACOS(-1.)))
        ELSEIF (Elem_xGP(1,i,j,k,iElem) .GT. xyzPhysicalMinMax(2)) THEN
          xi                  = ABS(Elem_xGP(1,i,j,k,iElem))-ABS(xyzPhysicalMinMax(2))
          L                   = ABS(xyzMinMax(2))-ABS(xyzPhysicalMinMax(2))
          PMLzeta(1,i,j,k,ElemToPML(iElem)) = PMLzeta0*(xi/L-SIN(2*ACOS(-1.)*xi/L)/(2*ACOS(-1.)))
        END IF
        ! y-PML region
        IF     (Elem_xGP(2,i,j,k,iElem) .LT. xyzPhysicalMinMax(3)) THEN
          xi                  = ABS(Elem_xGP(2,i,j,k,iElem))-ABS(xyzPhysicalMinMax(3))
          L                   = ABS(xyzMinMax(3))-ABS(xyzPhysicalMinMax(3))
          PMLzeta(2,i,j,k,ElemToPML(iElem)) = PMLzeta0*(xi/L-SIN(2*ACOS(-1.)*xi/L)/(2*ACOS(-1.)))
        ELSEIF (Elem_xGP(2,i,j,k,iElem) .GT. xyzPhysicalMinMax(4)) THEN
          xi                  = ABS(Elem_xGP(2,i,j,k,iElem))-ABS(xyzPhysicalMinMax(4))
          L                   = ABS(xyzMinMax(4))-ABS(xyzPhysicalMinMax(4))
          PMLzeta(2,i,j,k,ElemToPML(iElem)) = PMLzeta0*(xi/L-SIN(2*ACOS(-1.)*xi/L)/(2*ACOS(-1.)))
        END IF
        ! z-PML region
        IF     (Elem_xGP(3,i,j,k,iElem) .LT. xyzPhysicalMinMax(5)) THEN
          xi                  = ABS(Elem_xGP(3,i,j,k,iElem))-ABS(xyzPhysicalMinMax(5))
          L                   = ABS(xyzMinMax(5))-ABS(xyzPhysicalMinMax(5))
          PMLzeta(3,i,j,k,ElemToPML(iElem)) = PMLzeta0*(xi/L-SIN(2*ACOS(-1.)*xi/L)/(2*ACOS(-1.)))
        ELSEIF (Elem_xGP(3,i,j,k,iElem) .GT. xyzPhysicalMinMax(6)) THEN
          xi                  = ABS(Elem_xGP(3,i,j,k,iElem))-ABS(xyzPhysicalMinMax(6))
          L                   = ABS(xyzMinMax(6))-ABS(xyzPhysicalMinMax(6))
          PMLzeta(3,i,j,k,ElemToPML(iElem)) = PMLzeta0*(xi/L-SIN(2*ACOS(-1.)*xi/L)/(2*ACOS(-1.)))
        ENDIF 
  END DO; END DO; END DO; END DO !iElem,k,i,j
CASE(3) ! polynomial
        ! f''(0) = 0, f'(1) = 0
  DO iElem=1,PP_nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        ! x-PML region
        IF     (Elem_xGP(1,i,j,k,iElem) .LT. xyzPhysicalMinMax(1)) THEN
          xi                  = ABS(Elem_xGP(1,i,j,k,iElem))-ABS(xyzPhysicalMinMax(1))
          L                   = ABS(xyzMinMax(1))-ABS(xyzPhysicalMinMax(1))
          XiN                 = xi/L
          PMLzeta(1,i,j,k,ElemToPML(iElem)) = PMLzeta0*(-3*XiN**4+4*XiN**3)
        ELSEIF (Elem_xGP(1,i,j,k,iElem) .GT. xyzPhysicalMinMax(2)) THEN
          xi                  = ABS(Elem_xGP(1,i,j,k,iElem))-ABS(xyzPhysicalMinMax(2))
          L                   = ABS(xyzMinMax(2))-ABS(xyzPhysicalMinMax(2))
          XiN                 = xi/L
          PMLzeta(1,i,j,k,ElemToPML(iElem)) = PMLzeta0*(-3*XiN**4+4*XiN**3)
        END IF
        ! y-PML region
        IF     (Elem_xGP(2,i,j,k,iElem) .LT. xyzPhysicalMinMax(3)) THEN
          xi                  = ABS(Elem_xGP(2,i,j,k,iElem))-ABS(xyzPhysicalMinMax(3))
          L                   = ABS(xyzMinMax(3))-ABS(xyzPhysicalMinMax(3))
          XiN                 = xi/L
          PMLzeta(2,i,j,k,ElemToPML(iElem)) = PMLzeta0*(-3*XiN**4+4*XiN**3)
        ELSEIF (Elem_xGP(2,i,j,k,iElem) .GT. xyzPhysicalMinMax(4)) THEN
          xi                  = ABS(Elem_xGP(2,i,j,k,iElem))-ABS(xyzPhysicalMinMax(4))
          L                   = ABS(xyzMinMax(4))-ABS(xyzPhysicalMinMax(4))
          XiN                 = xi/L
          PMLzeta(2,i,j,k,ElemToPML(iElem)) = PMLzeta0*(-3*XiN**4+4*XiN**3)
        END IF
        ! z-PML region
        IF     (Elem_xGP(3,i,j,k,iElem) .LT. xyzPhysicalMinMax(5)) THEN
          xi                  = ABS(Elem_xGP(3,i,j,k,iElem))-ABS(xyzPhysicalMinMax(5))
          L                   = ABS(xyzMinMax(5))-ABS(xyzPhysicalMinMax(5))
          XiN                 = xi/L
          PMLzeta(3,i,j,k,ElemToPML(iElem)) = PMLzeta0*(-3*XiN**4+4*XiN**3)
        ELSEIF (Elem_xGP(3,i,j,k,iElem) .GT. xyzPhysicalMinMax(6)) THEN
          xi                  = ABS(Elem_xGP(3,i,j,k,iElem))-ABS(xyzPhysicalMinMax(6))
          L                   = ABS(xyzMinMax(6))-ABS(xyzPhysicalMinMax(6))
          PMLzeta(3,i,j,k,ElemToPML(iElem)) = PMLzeta0*(xi/L-SIN(2*ACOS(-1.)*xi/L)/(2*ACOS(-1.)))
          XiN                 = xi/L
        PMLzeta(3,i,j,k,ElemToPML(iElem)) = PMLzeta0*(-3*XiN**4+4*XiN**3)
      ENDIF 
END DO; END DO; END DO; END DO !iElem,k,i,j

CASE DEFAULT
!  CALL abort(__STAMP__,'Shape function for damping coefficient in PML region not specified!',999,999.)
END SELECT ! PMLzetaShape

!Test: Set All PMLzeta values for a direction to PMLzeta
IF (PMLspread.EQ.1) THEN
  DO iElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        IF (PMLzeta(1,i,j,k,iElem) .GT. 0.0) PMLzeta(:,i,j,k,iElem)=PMLzeta(1,i,j,k,iElem)
        IF (PMLzeta(2,i,j,k,iElem) .GT. 0.0) PMLzeta(:,i,j,k,iElem)=PMLzeta(2,i,j,k,iElem)
        IF (PMLzeta(3,i,j,k,iElem) .GT. 0.0) PMLzeta(:,i,j,k,iElem)=PMLzeta(3,i,j,k,iElem) 
  END DO; END DO; END DO; END DO !iElem,k,i,j
END IF

!PMLzetaNorm=.TRUE.
! Normalize zeta if multiple direction
IF (PMLzetaNorm) THEN
  DO iElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        zetaVecABS=SQRT(PMLzeta(1,i,j,k,iElem)**2 &
                       +PMLzeta(2,i,j,k,iElem)**2 &
                       +PMLzeta(3,i,j,k,iElem)**2 )
        zetaVec=MAX(PMLzeta(1,i,j,k,iElem),0.)
        zetaVec=MAX(PMLzeta(2,i,j,k,iElem),zetaVec)
        zetaVec=MAX(PMLzeta(3,i,j,k,iElem),zetaVec)
        PMLzeta(:,i,j,k,iElem) = PMLzeta(:,i,j,k,iElem)/zetaVecABS*zetaVec
  END DO; END DO; END DO; END DO !iElem,k,i,j
END IF

!===================================================================================================================================
! write PMLzeta field
!===================================================================================================================================
IF (PMLwritezeta.EQ.1) THEN
  OPEN(unit=110,file='PMLzeta.dat',status='unknown')
  DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        write(110,'(ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7)')&
                   Elem_xGP(1,i,j,k,PMLtoElem(iPMLElem)),'  ', &
                   Elem_xGP(2,i,j,k,PMLtoElem(iPMLElem)),'  ', &
                   Elem_xGP(3,i,j,k,PMLtoElem(iPMLElem)),'  ', &
                   PMLzeta(1,i,j,k,iPMLElem),'  ',             &
                   PMLzeta(2,i,j,k,iPMLElem),'  ',             &
                   PMLzeta(3,i,j,k,iPMLElem)
  END DO; END DO; END DO; END DO !iElem,k,j,i
  CLOSE(110)
END IF
!===================================================================================================================================
! find element index for probe points
!===================================================================================================================================
Probes%Coordinates(1,:)=(/5.9,5.9,5.9/)
Probes%Coordinates(2,:)=(/5.9,5.9,0.0/)
Probes%Coordinates(3,:)=(/5.9,0.0,0.0/)
DO iProbe=1,3;
  DO iElem=1,PP_nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    Probes%Distance(iProbe,i,j,k,iElem)=SQRT((Elem_xGP(1,i,j,k,iElem)-Probes%Coordinates(iProbe,1))**2+&
                                             (Elem_xGP(2,i,j,k,iElem)-Probes%Coordinates(iProbe,2))**2+&
                                             (Elem_xGP(3,i,j,k,iElem)-Probes%Coordinates(iProbe,3))**2)
  END DO; END DO; END DO; END DO;
  !print *, "MINVAL(Probes%Distance(iProbe,:,:,:,:))        ",MINVAL(Probes%Distance(iProbe,:,:,:,:))
  Probes%iElemMinLoc(iProbe,:)=MINLOC(Probes%Distance(iProbe,:,:,:,:))-(/1,1,1,0/) !weil I,J,K bei 0 beginnen
  !print *, "MINLOC(Probes%Distance(iProbe,:,:,:,:))        ",MINLOC(Probes%Distance(iProbe,:,:,:,:))
  Probes%Element(iProbe)=Probes%iElemMinLoc(iProbe,4)
  !print *, "Probes%Element(iProbe)                         ",Probes%Element(iProbe)
  !print *, "Probes%iElemMinLoc(iProbe,:)                 ",Probes%iElemMinLoc(iProbe,:)
END DO !iProbe,iElem,k,i,j
!NodeMap(:,1)=(/1,4,3,2/)
!RESULT = MINVAL(ARRAY, DIM [, MASK]) 

SWRITE(UNIT_stdOut,'(A)') ' INIT PML DONE...'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitPML



SUBROUTINE CalcPMLSource()
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals,       ONLY: abort
USE MOD_DG_Vars,       ONLY: Ut,U
USE MOD_PML_Vars,      ONLY: nPMLElems,PMLtoElem
USE MOD_PML_Vars,      ONLY: PMLzeta,U2
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iPMLElem
!===================================================================================================================================
DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N 
      Ut(1,i,j,k,PMLtoElem(iPMLElem)) = Ut(1,i,j,k,PMLtoElem(iPMLElem)) - &
((PMLzeta(2,i,j,k,iPMLElem)+PMLzeta(3,i,j,k,iPMLElem)-PMLzeta(1,i,j,k,iPMLElem))*U(1,i,j,k,PMLtoElem(iPMLElem))+&
(PMLzeta(1,i,j,k,iPMLElem)-PMLzeta(2,i,j,k,iPMLElem))*(PMLzeta(1,i,j,k,iPMLElem)-PMLzeta(3,i,j,k,iPMLElem))*U2(1,i,j,k,iPMLElem))
                         
      Ut(2,i,j,k,PMLtoElem(iPMLElem)) = Ut(2,i,j,k,PMLtoElem(iPMLElem)) - &
((PMLzeta(1,i,j,k,iPMLElem)+PMLzeta(3,i,j,k,iPMLElem)-PMLzeta(2,i,j,k,iPMLElem))*U(2,i,j,k,PMLtoElem(iPMLElem))+&
(PMLzeta(2,i,j,k,iPMLElem)-PMLzeta(1,i,j,k,iPMLElem))*(PMLzeta(2,i,j,k,iPMLElem)-PMLzeta(3,i,j,k,iPMLElem))*U2(2,i,j,k,iPMLElem))
                         
      Ut(3,i,j,k,PMLtoElem(iPMLElem)) = Ut(3,i,j,k,PMLtoElem(iPMLElem)) - &
((PMLzeta(1,i,j,k,iPMLElem)+PMLzeta(2,i,j,k,iPMLElem)-PMLzeta(3,i,j,k,iPMLElem))*U(3,i,j,k,PMLtoElem(iPMLElem))+&
(PMLzeta(3,i,j,k,iPMLElem)-PMLzeta(1,i,j,k,iPMLElem))*(PMLzeta(3,i,j,k,iPMLElem)-PMLzeta(2,i,j,k,iPMLElem))*U2(3,i,j,k,iPMLElem))
                         
      Ut(4,i,j,k,PMLtoElem(iPMLElem)) = Ut(4,i,j,k,PMLtoElem(iPMLElem)) - &
((PMLzeta(2,i,j,k,iPMLElem)+PMLzeta(3,i,j,k,iPMLElem)-PMLzeta(1,i,j,k,iPMLElem))*U(4,i,j,k,PMLtoElem(iPMLElem))+&
(PMLzeta(1,i,j,k,iPMLElem)-PMLzeta(2,i,j,k,iPMLElem))*(PMLzeta(1,i,j,k,iPMLElem)-PMLzeta(3,i,j,k,iPMLElem))*U2(4,i,j,k,iPMLElem))
                         
      Ut(5,i,j,k,PMLtoElem(iPMLElem)) = Ut(5,i,j,k,PMLtoElem(iPMLElem)) - &
((PMLzeta(1,i,j,k,iPMLElem)+PMLzeta(3,i,j,k,iPMLElem)-PMLzeta(2,i,j,k,iPMLElem))*U(5,i,j,k,PMLtoElem(iPMLElem))+&
(PMLzeta(2,i,j,k,iPMLElem)-PMLzeta(1,i,j,k,iPMLElem))*(PMLzeta(2,i,j,k,iPMLElem)-PMLzeta(3,i,j,k,iPMLElem))*U2(5,i,j,k,iPMLElem))
                         
      Ut(6,i,j,k,PMLtoElem(iPMLElem)) = Ut(6,i,j,k,PMLtoElem(iPMLElem)) - &
((PMLzeta(1,i,j,k,iPMLElem)+PMLzeta(2,i,j,k,iPMLElem)-PMLzeta(3,i,j,k,iPMLElem))*U(6,i,j,k,PMLtoElem(iPMLElem))+&
(PMLzeta(3,i,j,k,iPMLElem)-PMLzeta(1,i,j,k,iPMLElem))*(PMLzeta(3,i,j,k,iPMLElem)-PMLzeta(2,i,j,k,iPMLElem))*U2(6,i,j,k,iPMLElem))
END DO; END DO; END DO; END DO !nPMLElems,k,j,i
END SUBROUTINE CalcPMLSource



SUBROUTINE PMLTimeDerivative()
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
USE MOD_Globals,       ONLY : abort
USE MOD_PreProc
USE MOD_DG_Vars,       ONLY : U
USE MOD_PML_Vars,      ONLY : PMLzeta,U2,U2t
USE MOD_PML_Vars,      ONLY: nPMLElems,PMLtoElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iPMLElem
!===================================================================================================================================
! Add Source Terms
DO iPMLElem=1,nPMLElems; 
U2t(1,:,:,:,iPMLElem)=U(1,:,:,:,PMLtoElem(iPMLElem))-PMLzeta(1,:,:,:,iPMLElem)*U2(1,:,:,:,iPMLElem)
U2t(2,:,:,:,iPMLElem)=U(2,:,:,:,PMLtoElem(iPMLElem))-PMLzeta(2,:,:,:,iPMLElem)*U2(2,:,:,:,iPMLElem)
U2t(3,:,:,:,iPMLElem)=U(3,:,:,:,PMLtoElem(iPMLElem))-PMLzeta(3,:,:,:,iPMLElem)*U2(3,:,:,:,iPMLElem)
U2t(4,:,:,:,iPMLElem)=U(4,:,:,:,PMLtoElem(iPMLElem))-PMLzeta(1,:,:,:,iPMLElem)*U2(4,:,:,:,iPMLElem)
U2t(5,:,:,:,iPMLElem)=U(5,:,:,:,PMLtoElem(iPMLElem))-PMLzeta(2,:,:,:,iPMLElem)*U2(5,:,:,:,iPMLElem)
U2t(6,:,:,:,iPMLElem)=U(6,:,:,:,PMLtoElem(iPMLElem))-PMLzeta(3,:,:,:,iPMLElem)*U2(6,:,:,:,iPMLElem)
END DO ! iElem=1,PP_nElems,k,j,i
END SUBROUTINE PMLTimeDerivative



SUBROUTINE TransformPMLVars()
!===================================================================================================================================
! Transform E to E_tilde by adding the auxiliary PML fields to the State Vector
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars,       ONLY: U
USE MOD_PML_Vars,      ONLY: PMLzeta,U2
USE MOD_PML_Vars,      ONLY: nPMLElems,PMLtoElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iPMLElem
!===================================================================================================================================

DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        U(1,i,j,k,PMLtoElem(iPMLElem))=U(1,i,j,k,PMLtoElem(iPMLElem))+PMLzeta(1,i,j,k,iPMLElem)*U2(1,i,j,k,iPMLElem)
        U(2,i,j,k,PMLtoElem(iPMLElem))=U(2,i,j,k,PMLtoElem(iPMLElem))+PMLzeta(2,i,j,k,iPMLElem)*U2(2,i,j,k,iPMLElem)
        U(3,i,j,k,PMLtoElem(iPMLElem))=U(3,i,j,k,PMLtoElem(iPMLElem))+PMLzeta(3,i,j,k,iPMLElem)*U2(3,i,j,k,iPMLElem)
        U(4,i,j,k,PMLtoElem(iPMLElem))=U(4,i,j,k,PMLtoElem(iPMLElem))+PMLzeta(1,i,j,k,iPMLElem)*U2(4,i,j,k,iPMLElem)
        U(5,i,j,k,PMLtoElem(iPMLElem))=U(5,i,j,k,PMLtoElem(iPMLElem))+PMLzeta(2,i,j,k,iPMLElem)*U2(5,i,j,k,iPMLElem)
        U(6,i,j,k,PMLtoElem(iPMLElem))=U(6,i,j,k,PMLtoElem(iPMLElem))+PMLzeta(3,i,j,k,iPMLElem)*U2(6,i,j,k,iPMLElem)
END DO; END DO; END DO; END DO ! iElem=1,PP_nElems,k,j,i
END SUBROUTINE TransformPMLVars

SUBROUTINE BacktransformPMLVars()
!===================================================================================================================================
! Transform E_tilde to E by subtracting the auxiliary PML fields from the State Vector
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars,       ONLY: U
USE MOD_PML_Vars,      ONLY: PMLzeta,U2
USE MOD_PML_Vars,      ONLY: nPMLElems,PMLtoElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iPMLElem
!===================================================================================================================================

DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        U(1,i,j,k,PMLtoElem(iPMLElem))=U(1,i,j,k,PMLtoElem(iPMLElem))-PMLzeta(1,i,j,k,iPMLElem)*U2(1,i,j,k,iPMLElem)
        U(2,i,j,k,PMLtoElem(iPMLElem))=U(2,i,j,k,PMLtoElem(iPMLElem))-PMLzeta(2,i,j,k,iPMLElem)*U2(2,i,j,k,iPMLElem)
        U(3,i,j,k,PMLtoElem(iPMLElem))=U(3,i,j,k,PMLtoElem(iPMLElem))-PMLzeta(3,i,j,k,iPMLElem)*U2(3,i,j,k,iPMLElem)
        U(4,i,j,k,PMLtoElem(iPMLElem))=U(4,i,j,k,PMLtoElem(iPMLElem))-PMLzeta(1,i,j,k,iPMLElem)*U2(4,i,j,k,iPMLElem)
        U(5,i,j,k,PMLtoElem(iPMLElem))=U(5,i,j,k,PMLtoElem(iPMLElem))-PMLzeta(2,i,j,k,iPMLElem)*U2(5,i,j,k,iPMLElem)
        U(6,i,j,k,PMLtoElem(iPMLElem))=U(6,i,j,k,PMLtoElem(iPMLElem))-PMLzeta(3,i,j,k,iPMLElem)*U2(6,i,j,k,iPMLElem)
END DO; END DO; END DO; END DO ! iElem=1,PP_nElems,k,j,i
END SUBROUTINE BacktransformPMLVars



SUBROUTINE ProbePML(t)
!===================================================================================================================================
!  write E,B-field at probe points
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars,            ONLY : U
USE MOD_PML_Vars,           ONLY : Probes
!#ifdef PARTICLES
!USE MOD_Eval_xyz,           ONLY : eval_xyz
!#endif /*PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: field2(3,6)
INTEGER                         :: iProbe
!===================================================================================================================================
! calc Probe interpolation field
!===================================================================================================================================
!DO iProbe=1,3;
  !!CALL eval_xyz(Probes%Coordinates(iProbe,:),6,PP_N,U(1:6,:,:,:,Probes%Element(iProbe)),field(:),Probes%Element(iProbe))
  !field2(iProbe,:)=U(1:6,Probes%iElemMinLoc(iProbe,1),& ! i
                         !Probes%iElemMinLoc(iProbe,2),& ! j
                         !Probes%iElemMinLoc(iProbe,3),& ! k
                         !Probes%iElemMinLoc(iProbe,4))  ! iElem
!END DO !iProbe
!OPEN(unit=111,FILE='probes.dat',ACCESS='APPEND',STATUS='UNKNOWN')
!write(111,'(ES15.7,A,&
           !& ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A, &
           !& ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A, &
           !& ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7)') &
!t,'  ',&
!field2(1,1),'  ',field2(1,2),'  ',field2(1,3),'  ',field2(1,4),'  ',field2(1,5),'  ',field2(1,6),'  ',&
!field2(2,1),'  ',field2(2,2),'  ',field2(2,3),'  ',field2(2,4),'  ',field2(2,5),'  ',field2(2,6),'  ',&
!field2(3,1),'  ',field2(3,2),'  ',field2(3,3),'  ',field2(3,4),'  ',field2(3,5),'  ',field2(3,6)
!CLOSE(111)
END SUBROUTINE ProbePML



SUBROUTINE FinalizePML()
!===================================================================================================================================
!  
!===================================================================================================================================
! MODULES
USE MOD_PML_Vars,            ONLY: PMLzeta,U2,U2t
USE MOD_PML_Vars,            ONLY: ElemtoPML,PMLtoElem,DoPML
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(.NOT.DoPML) RETURN
SDEALLOCATE(PMLzeta)
SDEALLOCATE(U2)
SDEALLOCATE(U2t)
SDEALLOCATE(PMLtoElem)
SDEALLOCATE(ElemtoPML)
END SUBROUTINE FinalizePML




END MODULE MOD_PML

