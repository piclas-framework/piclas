PROGRAM RayBilinearPatch
!================================================================================================================================
! this program computes the intersection of a ray with an bi-linear patch
! first try is based on a certain paper      
!================================================================================================================================
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! nodes and particle 
REAL,DIMENSION(3,4)           :: xNode
REAL,DIMENSION(3)             :: PartState, lastPartState
!--------------------------------------------------------------------------------------------------------------------------------
! vector equation
REAL,DIMENSION(3)             :: q
!--------------------------------------------------------------------------------------------------------------------------------
! bi-linear equation
REAL,DIMENSION(3,4)           :: BiCoeff
!--------------------------------------------------------------------------------------------------------------------------------
! information for quatric equation solver
REAL                          :: A,B,C 
REAL,DIMENSION(4)             :: a1,a2
INTEGER                       :: nRoot
REAL,DIMENSION(2)             :: u,v,t
REAL,DIMENSION(3)             :: q1
LOGICAL,DIMENSION(2)          :: isIntersection
!--------------------------------------------------------------------------------------------------------------------------------
! pic paper
REAL                          :: beta1, beta2,beta3,beta4, lenq
!--------------------------------------------------------------------------------------------------------------------------------
! functions
REAL                          :: Computet
!================================================================================================================================


!--------------------------------------------------------------------------------------------------------------------------------
! insert nodes and particle position
!--------------------------------------------------------------------------------------------------------------------------------

! bi linear
xNode(1:3,1) = [0.0,0.0,1.8]
xNode(1:3,2) = [1.0,0.0,0.70]
xNode(1:3,3) = [0.0,1.0,0.70]
xNode(1:3,4) = [1.0,1.0,1.8]

!!! planar - required to be a right hand system
!xNode(1:3,1) = [0.0,0.0,1.8]
!xNode(1:3,2) = [1.0,0.0,1.8]
!xNode(1:3,3) = [1.0,1.0,1.8]
!xNode(1:3,4) = [0.0,1.0,1.8]

! wrong index
!xNode(1:3,4) = [0.0,0.0,1.8]
!xNode(1:3,1) = [1.0,0.0,1.8]
!xNode(1:3,2) = [1.0,1.0,1.8]
!xNode(1:3,3) = [0.0,1.0,1.8]

!lastPartState = [0.8,0.7,0.8]
!PartState     = [0.5,0.5,2.5]


!lastPartState = [0.8,0.7,0.8]
!PartState     = [0.8,0.7,2.5]

!PartState     = [0.8,0.7,0.8]
!lastPartState = [0.8,0.7,2.5]

PartState     = [0.4,0.7,0.4]
lastPartState = [0.4,0.7,2.2]

! two intersections
!PartState     = [0.8,0.9,1.4]
!lastPartState = [0.1,0.2,1.2]


!--------------------------------------------------------------------------------------------------------------------------------
! primarily computations
!--------------------------------------------------------------------------------------------------------------------------------

! particle vector
! r+t*(PartState-lastPartState)
! r+t*q
q = PartState - lastPartState

! bi-linear equation

! setting the first node to p0
! can be saved per side
BiCoeff(:,1) = xNode(:,4) - xNode(:,3) - xNode(:,2) + xNode(:,1)
BiCoeff(:,2) = xNode(:,3) - xNode(:,1)
BiCoeff(:,3) = xNode(:,2) - xNode(:,1)
BiCoeff(:,4) = xNode(:,1)

!--------------------------------------------------------------------------------------------------------------------------------
! ray intersection with patch, Ramsey 2004 paper
!--------------------------------------------------------------------------------------------------------------------------------

WRITE(*,*) '--------------------------------------------------------------------------------------------------------------------'
WRITE(*,*) '    Ray alogirthm                                                                                                   '
WRITE(*,*) '--------------------------------------------------------------------------------------------------------------------'
WRITE(*,*) ''

! prepare solver  ! has to be computed for each patch
a1(1)= BiCoeff(1,1)*q(3) - BiCoeff(3,1)*q(1)
a1(2)= BiCoeff(1,2)*q(3) - BiCoeff(3,2)*q(1)
a1(3)= BiCoeff(1,3)*q(3) - BiCoeff(3,3)*q(1)
a1(4)= (BiCoeff(1,4)-lastPartState(1))*q(3) &
     - (BiCoeff(3,4)-lastpartState(3))*q(1)

a2(1)= BiCoeff(2,1)*q(3) - BiCoeff(3,1)*q(2)
a2(2)= BiCoeff(2,2)*q(3) - BiCoeff(3,2)*q(2)
a2(3)= BiCoeff(2,3)*q(3) - BiCoeff(3,3)*q(2)
a2(4)= (BiCoeff(2,4)-lastPartState(2))*q(3) &
     - (BiCoeff(3,4)-lastpartState(3))*q(2)

A = a2(1)*a1(3)-a1(1)*a2(3)
B = a2(1)*a1(4)-a1(1)*a2(4)+a2(2)*a1(3)-a1(2)*a2(3)
C = a1(4)*a2(2)-a1(2)*a2(4)
print*,'A,B,C', A,B,C
CALL QuatricSolver(A,B,C,nRoot,v(1),v(2))
print*,'nRoot,v', nRoot,v

isIntersection=.FALSE.
IF(nRoot.EQ.0)THEN
  WRITE(*,*) ' no intersection'
ELSE IF (nRoot.EQ.1) THEN
  IF((v(1).GE.0.).AND.v(1).LT.1.)THEN
    u(1)=v(1)*(a2(1)-a1(1))+a2(2)-a1(2)
    u(1)=1.0/u(1)
    u(1)=(v(1)*(a1(3)-a2(3))+a1(4)-a2(4))*u(1)
    IF((u(1).GE.0).AND.u(1).LT.1)THEN
      q1=u(1)*v(1)*BiCoeff(:,1)+u(1)*BiCoeff(:,2)+v(1)*BiCoeff(:,3)+BiCoeff(:,4)-lastPartState
      t(1)=Computet(q1,q)
      IF((t(1).GE.0.).AND.(t(1).LE.1.0))THEN
        WRITE(*,*) ' One Intersection'
        WRITE(*,*) ' t ', t(1)
        WRITE(*,*) ' Intersection at ', lastPartState+t(1)*q
        isIntersection(1)=.TRUE.
      END IF 
    END IF
  END IF
ELSE 
  IF((v(1).GE.0.).AND.v(1).LT.1.)THEN
    !u(1)=v(1)*a2(1)+a2(2)
    !u(1)=1.0/u(1)
    !u(1)=(-v(1)*a2(3)-a2(4))*u(1)
    u(1)=v(1)*(a2(1)-a1(1))+a2(2)-a1(2)
    u(1)=1.0/u(1)
    u(1)=(v(1)*(a1(3)-a2(3))+a1(4)-a2(4))*u(1)
    IF((u(1).GE.0.).AND.u(1).LT.1.)THEN
      q1=u(1)*v(1)*BiCoeff(:,1)+u(1)*BiCoeff(:,2)+v(1)*BiCoeff(:,3)+BiCoeff(:,4)-lastPartState
      t(1)=Computet(q1,q)
      IF((t(1).GE.0.).AND.(t(1).LE.1.0))THEN
        WRITE(*,*) ' One Intersection'
        WRITE(*,*) ' t ', t(1)
        WRITE(*,*) ' Intersection at ', lastPartState+t(1)*q
        isIntersection(1)=.TRUE.
      END IF 
    END IF
  END IF
  IF((v(2).GE.0.).AND.v(2).LT.1.)THEN
    u(2)=v(2)*a2(1)-v(2)*a1(1)+a2(2)-a1(2)
    u(2)=1.0/u(2)
    u(2)=(v(2)*a1(3)-v(2)*a2(3)+a1(4)-a2(4))*u(2)
    IF((u(2).GE.0.).AND.u(2).LT.1.)THEN
      q1=u(2)*v(2)*BiCoeff(:,1)+u(2)*BiCoeff(:,2)+v(2)*BiCoeff(:,3)+BiCoeff(:,4)-lastPartState
      t(2)=Computet(q1,q)
      IF((t(2).GE.0.).AND.(t(2).LE.1.0))THEN
        WRITE(*,*) ' Second Intersection'
        WRITE(*,*) ' t ', t(2)
        WRITE(*,*) ' Intersection at ', lastPartState+t(2)*q
        isIntersection(2)=.TRUE.
      END IF 
    END IF
  END IF
END IF
WRITE(*,*) ' Following intersections ', isIntersection


!--------------------------------------------------------------------------------------------------------------------------------
! pic paper from Haselbach 2007 "An efficient and robust particle-localization for unstructured grids"
!--------------------------------------------------------------------------------------------------------------------------------

WRITE(*,*) '--------------------------------------------------------------------------------------------------------------------'
WRITE(*,*) '    PIC alogirthm                                                                                                   '
WRITE(*,*) '--------------------------------------------------------------------------------------------------------------------'
WRITE(*,*) ''

! prepare side coefficients
beta1 = BiCoeff(1,1)*q(1)+BiCoeff(2,1)*q(2)+BiCoeff(3,1)*q(3)
beta2 = BiCoeff(1,2)*q(1)+BiCoeff(2,2)*q(2)+BiCoeff(3,2)*q(3)
beta3 = BiCoeff(1,3)*q(1)+BiCoeff(2,3)*q(2)+BiCoeff(3,3)*q(3)
beta4 = (BiCoeff(1,4)-lastPartState(1))*q(1) &
      + (BiCoeff(2,4)-lastPartState(2))*q(2) &
      + (BiCoeff(3,4)-lastPartState(3))*q(3)
! xz
a1(1)=BiCoeff(3,1)-beta1*q(3)-BiCoeff(1,1)+beta1*q(1)
a1(2)=BiCoeff(3,2)-beta2*q(3)-BiCoeff(1,2)+beta2*q(1)
a1(3)=BiCoeff(3,3)-beta3*q(3)-BiCoeff(1,3)+beta3*q(1)
a1(4)=BiCoeff(3,4)-beta4*q(3)-BiCoeff(1,4)+beta4*q(1)

! yz
a2(1)=BiCoeff(3,1)-beta1*q(3)-BiCoeff(2,1)+beta1*q(2)
a2(2)=BiCoeff(3,2)-beta2*q(3)-BiCoeff(2,2)+beta2*q(2)
a2(3)=BiCoeff(3,3)-beta3*q(3)-BiCoeff(2,3)+beta3*q(2)
a2(4)=BiCoeff(3,4)-beta4*q(3)-BiCoeff(2,4)+beta4*q(2)


A = a1(1)*a2(3)-a2(1)*a1(3)
B = a1(1)*a2(4)-a2(1)*a1(4)+a1(2)*a2(3)-a2(2)*a1(3)
C = a1(2)*a1(4)-a2(2)*a1(4)
!print*,'A,B,C',A,B,C
CALL QuatricSolver(A,B,C,nRoot,v(1),v(2))
!print*,nRoot

isIntersection=.FALSE.
IF(nRoot.EQ.0)THEN
  WRITE(*,*) ' no intersection'
ELSE IF (nRoot.EQ.1) THEN
  IF((v(1).GE.0.).AND.v(1).LT.1.)THEN
    u(1)=v(1)*a1(1)+a1(2)
    u(1)=1.0/u(1)
    u(1)=(-v(1)*a1(3)-a1(4))*u(1)
    IF((u(1).GE.0).AND.u(1).LT.1)THEN
      t(1)=(u(1)*v(1)*BiCoeff(1,1)+u(1)*BiCoeff(1,2)+v(1)*BiCoeff(1,3)+BiCoeff(1,4)-lastPartState(1))*q(1) &
          +(u(1)*v(1)*BiCoeff(2,1)+u(1)*BiCoeff(2,2)+v(1)*BiCoeff(2,3)+BiCoeff(2,4)-lastPartState(2))*q(2) &
          +(u(1)*v(1)*BiCoeff(3,1)+u(1)*BiCoeff(3,2)+v(1)*BiCoeff(3,3)+BiCoeff(3,4)-lastPartState(2))*q(3)
      IF((t(1).GE.0.).AND.(t(1).LE.1.0))THEN
        WRITE(*,*) ' One Intersection'
        WRITE(*,*) ' t ', t(1)
        WRITE(*,*) ' Intersection at ', lastPartState+t(1)*q
        isIntersection(1)=.TRUE.
      END IF 
    END IF
  END IF
ELSE 
  IF((v(1).GE.0.).AND.v(1).LT.1.)THEN
    u(1)=v(1)*a1(1)+a1(2)
    u(1)=1.0/u(1)
    u(1)=(-v(1)*a1(3)-a1(4))*u(1)
    IF((u(1).GE.0.).AND.u(1).LT.1.)THEN
      t(1)=(u(1)*v(1)*BiCoeff(1,1)+u(1)*BiCoeff(1,2)+v(1)*BiCoeff(1,3)+BiCoeff(1,4)-lastPartState(1))*q(1) &
          +(u(1)*v(1)*BiCoeff(2,1)+u(1)*BiCoeff(2,2)+v(1)*BiCoeff(2,3)+BiCoeff(2,4)-lastPartState(2))*q(2) &
          +(u(1)*v(1)*BiCoeff(3,1)+u(1)*BiCoeff(3,2)+v(1)*BiCoeff(3,3)+BiCoeff(3,4)-lastPartState(2))*q(3)

      IF((t(1).GE.0.).AND.(t(1).LE.1.0))THEN
        WRITE(*,*) ' One Intersection'
        WRITE(*,*) ' t ', t(1)
        WRITE(*,*) ' Intersection at ', lastPartState+t(1)*q
        isIntersection(1)=.TRUE.
      END IF 
    END IF
  END IF
  IF((v(2).GE.0.).AND.v(2).LT.1.)THEN
    u(2)=v(2)*a1(1)+a1(2)
    u(2)=1.0/u(2)
    u(2)=(-v(2)*a1(3)-a1(4))*u(2)
    IF((u(2).GE.0.).AND.u(2).LT.1.)THEN
      t(2)=(u(2)*v(2)*BiCoeff(1,1)+u(2)*BiCoeff(1,2)+v(2)*BiCoeff(1,3)+BiCoeff(1,4)-lastPartState(1))*q(1) &
          +(u(2)*v(2)*BiCoeff(2,1)+u(2)*BiCoeff(2,2)+v(2)*BiCoeff(2,3)+BiCoeff(2,4)-lastPartState(2))*q(2) &
          +(u(2)*v(2)*BiCoeff(3,1)+u(2)*BiCoeff(3,2)+v(2)*BiCoeff(3,3)+BiCoeff(3,4)-lastPartState(2))*q(3)
      IF((t(2).GE.0.).AND.(t(2).LE.1.0))THEN
        WRITE(*,*) ' Second Intersection'
        WRITE(*,*) ' t ', t(2)
        WRITE(*,*) ' Intersection at ', lastPartState+t(2)*q
        isIntersection(2)=.TRUE.
      END IF 
    END IF
  END IF
END IF
WRITE(*,*) ' Following intersections ', isIntersection


!================================================================================================================================
END PROGRAM RayBilinearPatch

SUBROUTINE QuatricSolver(A,B,C,nRoot,r1,r2)
!================================================================================================================================
! subroutine to compute the modified a,b,c equation, parameter already mapped in final version
!================================================================================================================================
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)         :: A,B,C
!--------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(OUT)     :: nRoot
REAL,INTENT(OUT)        :: R1,R2
!--------------------------------------------------------------------------------------------------------------------------------
! local variables
REAL                    :: eps=1e-12
!================================================================================================================================

IF(ABS(a).LT.eps)THEN
  nRoot=1
  R1=-c/b
  R2=0.
ELSE
  IF(B.LT.0) THEN
    nRoot = 0
    R1=0.
    R2=0.
  ELSE IF (B.EQ.0)THEN
    nRoot =1
    R1 = -0.5*B/A
    R2 = 0.
  ELSE
    nRoot=2
    R1 = SQRT(B*B-4.0*A*C)
    R2 = -R1
    R1 = -B+R1
    R1 = 0.5*R1/A
    R2 = -B+R2
    R2 = 0.5*R2/A
  END IF
END IF

END SUBROUTINE QuatricSolver

FUNCTION Computet(q1,q)
!================================================================================================================================
! compute the required vector length to intersection
!================================================================================================================================
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(3),INTENT(IN)         :: q1,q
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                                 :: Computet
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: t
REAL                                 :: eps=1e-6
!================================================================================================================================

IF((ABS(q(1)).GE.ABS(q(2))).AND.(ABS(q(1)).GT.ABS(q(3))))THEN
  t = q1(1)/q(1)-eps
ELSE IF(ABS(q(2)).GE.ABS(q(3)))THEN
  t = q1(2)/q(2)-eps
ELSE
  t = q1(3)/q(3)-eps
END IF

Computet=t

END FUNCTION Computet

