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
! own pic apporach
REAL,DIMENSION(2)             :: Xi,Eta
REAL,DIMENSION(3)             :: LocCoeff
!--------------------------------------------------------------------------------------------------------------------------------
! pic paper
REAL                          :: beta1, beta2,beta3,beta4, lenq
REAL,DIMENSION(3)             :: as,bs,cs,ds,dummy
!--------------------------------------------------------------------------------------------------------------------------------
! functions
REAL                          :: Computet
REAL                          :: epsZero,epsOne,eps
LOGICAL                       :: inPlane
!--------------------------------------------------------------------------------------------------------------------------------
! speed test of different implementations
INTEGER                       :: iloop,nloop
REAL                          :: tstart, tend, Ran
REAL,DIMENSION(4)             :: time
INTEGER,DIMENSION(4)          :: nInter
!================================================================================================================================

INTEGER, DIMENSION(:), ALLOCATABLE :: seed
INTEGER::n

CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))
seed=24324210
CALL RANDOM_SEED(PUT = seed)

epsZero=-1e-12
eps    =1e-12
epsOne =1.0+1e-12

time=0.
nInter=3
nloop=1e7
!nloop=100

DO iloop=1,nloop

  !--------------------------------------------------------------------------------------------------------------------------------
  ! insert nodes and particle position
  !--------------------------------------------------------------------------------------------------------------------------------
  
  !jWRITE(*,'(A,10x,I5,13x,I5,13x,I5)') ' Intersections', nInter(1), nInter(2), nInter(3)
  !jread*
  ! bi linear
  xNode(1:3,1) = [0.0,0.0,1.8]
  xNode(1:3,2) = [1.0,0.0,0.70]
  xNode(1:3,3) = [0.0,1.0,0.70]
  xNode(1:3,4) = [1.0,1.0,1.8]
  
  ! with random numbers
  CALL RANDOM_NUMBER(RAN)
  lastPartState(1) = 0.8-RAN
  CALL RANDOM_NUMBER(RAN)
  lastPartState(2) = 1.3*RAN
  CALL RANDOM_NUMBER(RAN)
  lastPartState(3) = RAN
  
  CALL RANDOM_NUMBER(RAN)
  PartState(1) = (1.0-RAN)*2.13
  CALL RANDOM_NUMBER(RAN)
  PartState(2) = RAN
  CALL RANDOM_NUMBER(RAN)
  PartState(3) = (1.0-RAN)*3.212754

  !print*,PartState
  !print*,lastPartState

  !--------------------------------------------------------------------------------------------------------------------------------
  ! primarily computations
  !--------------------------------------------------------------------------------------------------------------------------------
  

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
  
  CALL CPU_TIME(tstart)
  ! particle vector
  ! r+t*(PartState-lastPartState)
  ! r+t*q
  q = PartState - lastPartState

  
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
  !print*,'A,B,C', A,B,C
  CALL QuatricSolver(A,B,C,nRoot,v(1),v(2))
  !print*,'nRoot,v', nRoot,v
  
  inPlane=.FALSE.
  IF(nRoot.EQ.0)THEN
    IF((ABS(A).LT.eps).AND.(ABS(B).LT.eps).AND.(ABS(C).LT.eps))THEN
      inPlane=.TRUE.
    END IF
  END IF
  
  isIntersection=.FALSE.
  IF(nRoot.EQ.0)THEN
    IF(inPlane)THEN
      nInter(1)=nInter(1)+1
    END IF
  ELSE IF (nRoot.EQ.1) THEN
   IF((v(1).GE.epsZero).AND.v(1).LT.epsOne)THEN
     u(1)=v(1)*(a2(1)-a1(1))+a2(2)-a1(2)
     u(1)=1.0/u(1)
     u(1)=(v(1)*(a1(3)-a2(3))+a1(4)-a2(4))*u(1)
     !print*,'u(1)',u(1)
     IF((u(1).GE.epsZero).AND.u(1).LT.epsOne)THEN
       q1=u(1)*v(1)*BiCoeff(:,1)+u(1)*BiCoeff(:,2)+v(1)*BiCoeff(:,3)+BiCoeff(:,4)-lastPartState
       t(1)=Computet(q1,q)
       IF((t(1).GE.0.).AND.(t(1).LE.epsOne))THEN
         !WRITE(*,*) ' One Intersection'
         !WRITE(*,*) ' t ', t(1)
         !WRITE(*,*) ' Intersection at ', lastPartState+t(1)*q
         !isIntersection(1)=.TRUE.
         nInter(1)=nInter(1)+1
       END IF 
     END IF
   END IF
  ELSE 
    IF((v(1).GE.epsZero).AND.v(1).LT.epsOne)THEN
      !u(1)=v(1)*a2(1)+a2(2)
      !u(1)=1.0/u(1)
      !u(1)=(-v(1)*a2(3)-a2(4))*u(1)
      u(1)=v(1)*(a2(1)-a1(1))+a2(2)-a1(2)
      u(1)=1.0/u(1)
      u(1)=(v(1)*(a1(3)-a2(3))+a1(4)-a2(4))*u(1)
      IF((u(1).GE.epsZero).AND.u(1).LT.epsOne)THEN
        q1=u(1)*v(1)*BiCoeff(:,1)+u(1)*BiCoeff(:,2)+v(1)*BiCoeff(:,3)+BiCoeff(:,4)-lastPartState
        t(1)=Computet(q1,q)
        IF((t(1).GE.epsZero).AND.(t(1).LE.epsOne))THEN
          !WRITE(*,*) ' One Intersection'
          !WRITE(*,*) ' t ', t(1)
          !WRITE(*,*) ' Intersection at ', lastPartState+t(1)*q
          !isIntersection(1)=.TRUE.
          nInter(1)=nInter(1)+1
        END IF 
      END IF
    END IF
    IF((v(2).GE.epsZero).AND.v(2).LT.epsOne)THEN
      u(2)=v(2)*a2(1)-v(2)*a1(1)+a2(2)-a1(2)
      u(2)=1.0/u(2)
      u(2)=(v(2)*a1(3)-v(2)*a2(3)+a1(4)-a2(4))*u(2)
      IF((u(2).GE.epsZero).AND.u(2).LT.epsOne)THEN
        q1=u(2)*v(2)*BiCoeff(:,1)+u(2)*BiCoeff(:,2)+v(2)*BiCoeff(:,3)+BiCoeff(:,4)-lastPartState
        t(2)=Computet(q1,q)
        IF((t(2).GE.epsZero).AND.(t(2).LE.epsOne))THEN
          !WRITE(*,*) ' Second Intersection'
          !WRITE(*,*) ' t ', t(2)
          !WRITE(*,*) ' Intersection at ', lastPartState+t(2)*q
          !isIntersection(2)=.TRUE.
          nInter(1)=nInter(1)+1
        END IF 
      END IF
    END IF
  END IF
  CALL CPU_TIME(tend)   
  time(1)=time(1)+tend-tstart
  
  !--------------------------------------------------------------------------------------------------------------------------------
  ! pic paper from Haselbach 2007 "An efficient and robust particle-localization for unstructured grids"
  !--------------------------------------------------------------------------------------------------------------------------------
  
  CALL CPU_TIME(tstart)
  ! particle vector
  ! r+t*(PartState-lastPartState)
  ! r+t*q
  q = PartState - lastPartState

  
  ! prepare side coefficients
  ! caution override D
  lenq=q(1)*q(1)+q(2)*q(2)+q(3)*q(3)
  lenq=SQRT(lenq)
  q=q/lenq
  BiCoeff(:,4) = BiCoeff(:,4) - lastPartState
  beta1 = BiCoeff(1,1)*q(1)+BiCoeff(2,1)*q(2)+BiCoeff(3,1)*q(3)
  beta2 = BiCoeff(1,2)*q(1)+BiCoeff(2,2)*q(2)+BiCoeff(3,2)*q(3)
  beta3 = BiCoeff(1,3)*q(1)+BiCoeff(2,3)*q(2)+BiCoeff(3,3)*q(3)
  beta4 = BiCoeff(1,4)*q(1)+BiCoeff(2,4)*q(2)+BiCoeff(3,4)*q(3)
  
  ! debug
  as = BiCoeff(:,1) - beta1*q
  bs = BiCoeff(:,2) - beta2*q
  cs = BiCoeff(:,3) - beta3*q
  ds = BiCoeff(:,4) - beta4*q
  
  a1(1) = as(3)-as(1)
  a1(2) = bs(3)-bs(1)
  a1(3) = cs(3)-cs(1)
  a1(4) = ds(3)-ds(1)
  
  a2(1) = as(3)-as(2)
  a2(2) = bs(3)-bs(2)
  a2(3) = cs(3)-cs(2)
  a2(4) = ds(3)-ds(2)
  
  !!!!! old stuff with nice error
  !! xz
  ! a1(1)= BiCoeff(3,1)-beta1*q(3)-BiCoeff(1,1)+beta1*q(1)
  ! a1(2)= BiCoeff(3,2)-beta2*q(3)-BiCoeff(1,2)+beta2*q(1)
  ! a1(3)= BiCoeff(3,3)-beta3*q(3)-BiCoeff(1,3)+beta3*q(1)
  ! a1(4)= BiCoeff(3,4)-beta4*q(3)-BiCoeff(1,4)+beta3*q(1)
  
  ! yz
  !a2(1)= BiCoeff(3,1)-beta1*q(3)-BiCoeff(2,1)+beta1*q(2)
  !a2(2)= BiCoeff(3,2)-beta2*q(3)-BiCoeff(2,2)+beta2*q(2)
  !a2(3)= BiCoeff(3,3)-beta3*q(3)-BiCoeff(2,3)+beta3*q(2)
  !a2(4)= BiCoeff(3,4)-beta4*q(3)-BiCoeff(2,4)+beta3*q(2)
  
  A = a1(1)*a2(3)-a2(1)*a1(3)
  B = a1(1)*a2(4)-a2(1)*a1(4)+a1(2)*a2(3)-a2(2)*a1(3)
  C = a1(2)*a2(4)-a2(2)*a1(4)
  !print*,'A,B,C',A,B,C
  !print*,' help of abc-formula',B*B-4.0*A*C
  !print*,' sqrt', SQRT(b*b-4.0*a*c)
  CALL QuatricSolver(A,B,C,nRoot,v(1),v(2))
  !print*,nRoot,v
  
  inPlane=.FALSE.
  IF(nRoot.EQ.0)THEN
    IF((ABS(A).LT.eps).AND.(ABS(B).LT.eps).AND.(ABS(C).LT.eps))THEN
      inPlane=.TRUE.
    END IF
  END IF
  
  lenq=lenq+eps
  isIntersection=.FALSE.
  IF(nRoot.EQ.0)THEN
    IF(inPlane)THEN
      nInter(2)=nInter(2)+1
    END IF
  ELSE IF (nRoot.EQ.1) THEN
    IF((v(1).GE.epsZero).AND.(v(1).LT.epsOne))THEN
      u(1)=v(1)*(a2(1)-a1(1))+a2(2)-a1(2)
      u(1)=1.0/u(1)
      u(1)=(v(1)*(a1(3)-a2(3))+a1(4)-a2(4))*u(1)
      IF((u(1).GE.0).AND.u(1).LT.epsOne)THEN
        t(1)=(u(1)*v(1)*BiCoeff(1,1)+u(1)*BiCoeff(1,2)+v(1)*BiCoeff(1,3)+BiCoeff(1,4))*q(1) &
            +(u(1)*v(1)*BiCoeff(2,1)+u(1)*BiCoeff(2,2)+v(1)*BiCoeff(2,3)+BiCoeff(2,4))*q(2) &
            +(u(1)*v(1)*BiCoeff(3,1)+u(1)*BiCoeff(3,2)+v(1)*BiCoeff(3,3)+BiCoeff(3,4))*q(3)
        IF((t(1).GE.epsZero).AND.(t(1).LE.lenq))THEN
          !WRITE(*,*) ' One Intersection'
          !!WRITE(*,*) ' t ', t(1)
          !WRITE(*,*) ' Intersection at ', lastPartState+t(1)*q
          !isIntersection(1)=.TRUE.
          nInter(2)=nInter(2)+1
        END IF 
      END IF
    END IF
  ELSE 
    IF((v(1).GE.epsZero).AND.v(1).LT.epsOne)THEN
      !u(1)=v(1)*a1(1)+a1(2)
      !u(1)=1.0/u(1)
      !u(1)=(-v(1)*a1(3)-a1(4))*u(1)
      u(1)=v(1)*(a2(1)-a1(1))+a2(2)-a1(2)
      u(1)=1.0/u(1)
      u(1)=(v(1)*(a1(3)-a2(3))+a1(4)-a2(4))*u(1)
      IF((u(1).GE.epsZero).AND.u(1).LT.epsOne)THEN
        t(1)=(u(1)*v(1)*BiCoeff(1,1)+u(1)*BiCoeff(1,2)+v(1)*BiCoeff(1,3)+BiCoeff(1,4))*q(1) &
            +(u(1)*v(1)*BiCoeff(2,1)+u(1)*BiCoeff(2,2)+v(1)*BiCoeff(2,3)+BiCoeff(2,4))*q(2) &
            +(u(1)*v(1)*BiCoeff(3,1)+u(1)*BiCoeff(3,2)+v(1)*BiCoeff(3,3)+BiCoeff(3,4))*q(3)
        IF((t(1).GE.epsZero).AND.(t(1).LE.lenq))THEN
          !WRITE(*,*) ' One Intersection'
          !WRITE(*,*) ' t ', t(1)
          !WRITE(*,*) ' Intersection at ', lastPartState+t(1)*q
          !isIntersection(1)=.TRUE.
          nInter(2)=nInter(2)+1
        END IF 
      END IF
    END IF
    IF((v(2).GE.epsZero).AND.v(2).LT.epsOne)THEN
      u(2)=v(2)*(a2(1)-a1(1))+a2(2)-a1(2)
      u(2)=1.0/u(2)
      u(2)=(v(2)*(a1(3)-a2(3))+a1(4)-a2(4))*u(2)
      !u(2)=v(2)*a1(1)+a1(2)
      !u(2)=1.0/u(2)
      !u(2)=(-v(2)*a1(3)-a1(4))*u(2)
      IF((u(2).GE.epsZero).AND.u(2).LT.epsOne)THEN
        t(2)=(u(2)*v(2)*BiCoeff(1,1)+u(2)*BiCoeff(1,2)+v(2)*BiCoeff(1,3)+BiCoeff(1,4))*q(1) &
            +(u(2)*v(2)*BiCoeff(2,1)+u(2)*BiCoeff(2,2)+v(2)*BiCoeff(2,3)+BiCoeff(2,4))*q(2) &
            +(u(2)*v(2)*BiCoeff(3,1)+u(2)*BiCoeff(3,2)+v(2)*BiCoeff(3,3)+BiCoeff(3,4))*q(3)
        IF((t(2).GE.epsZero).AND.(t(2).LE.lenq))THEN
          !WRITE(*,*) ' Second Intersection'
          !WRITE(*,*) ' t ', t(2)
          !WRITE(*,*) ' Intersection at ', lastPartState+t(2)*q
          !isIntersection(2)=.TRUE.
          nInter(2)=nInter(2)+1
        END IF 
      END IF
    END IF
  END IF
  !WRITE(*,*) ' Following intersections ', isIntersection
  CALL CPU_TIME(tend)   
  time(2)=time(2)+tend-tstart  
 
  !--------------------------------------------------------------------------------------------------------------------------------
  ! own algorithm based on 
  ! pic paper from Haselbach 2007 "An efficient and robust particle-localization for unstructured grids"
  !--------------------------------------------------------------------------------------------------------------------------------
  
  ! preliminary computation
  ! resort the nodes counter-clock-wise || mathematic positive
  dummy        = xNode(1:3,3)
  xNode(1:3,1) = xNode(1:3,1)
  xNode(1:3,2) = xNode(1:3,2)
  xNode(1:3,3) = xNode(1:3,4)
  xNode(1:3,4) = dummy
  
  BiCoeff(:,1) = 0.25*( xNode(:,1) -xNode(:,2) +xNode(:,3) -xNode(:,4))
  BiCoeff(:,2) = 0.25*(-xNode(:,1) +xNode(:,2) +xNode(:,3) -xNode(:,4))
  BiCoeff(:,3) = 0.25*(-xNode(:,1) -xNode(:,2) +xNode(:,3) +xNode(:,4))
  BiCoeff(:,4) = 0.25*( xNode(:,1) +xNode(:,2) +xNode(:,3) +xNode(:,4))
  
  CALL CPU_TIME(tstart)
  ! prepare side coefficients
  q = PartState - lastPartState
  lenq=q(1)*q(1)+q(2)*q(2)+q(3)*q(3)
  lenq=SQRT(lenq)
  q=q/lenq
  LocCoeff = BiCoeff(:,4) - lastPartState
  beta1 = BiCoeff(1,1)*q(1)+BiCoeff(2,1)*q(2)+BiCoeff(3,1)*q(3)
  beta2 = BiCoeff(1,2)*q(1)+BiCoeff(2,2)*q(2)+BiCoeff(3,2)*q(3)
  beta3 = BiCoeff(1,3)*q(1)+BiCoeff(2,3)*q(2)+BiCoeff(3,3)*q(3)
  beta4 =  LocCoeff(1)*q(1)+ LocCoeff(2)*q(2)+ LocCoeff(3)*q(3)
  
  ! debug
  as = BiCoeff(:,1) - beta1*q
  bs = BiCoeff(:,2) - beta2*q
  cs = BiCoeff(:,3) - beta3*q
  ds =  LocCoeff(:) - beta4*q
  
  a1(1) = as(3)-as(1)
  a1(2) = bs(3)-bs(1)
  a1(3) = cs(3)-cs(1)
  a1(4) = ds(3)-ds(1)
  
  a2(1) = as(3)-as(2)
  a2(2) = bs(3)-bs(2)
  a2(3) = cs(3)-cs(2)
  a2(4) = ds(3)-ds(2)
  
  A = a1(1)*a2(3)-a2(1)*a1(3)
  B = a1(1)*a2(4)-a2(1)*a1(4)+a1(2)*a2(3)-a2(2)*a1(3)
  C = a1(2)*a2(4)-a2(2)*a1(4)

  !print*,' help of abc-formula',B*B-4.0*A*C
  ! print*,' sqrt', SQRT(b*b-4.0*a*c)
  CALL QuatricSolver(A,B,C,nRoot,Eta(1),Eta(2))
  !print*,nRoot,Eta
!  IF(iloop.EQ.34)THEN
!    print*,eta
!  END IF
  
  inPlane=.FALSE.
  IF(nRoot.EQ.0)THEN
    IF((ABS(A).LT.eps).AND.(ABS(B).LT.eps).AND.(ABS(C).LT.eps))THEN
      inPlane=.TRUE.
    END IF
  END IF
  
  lenq=lenq+eps
  isIntersection=.FALSE.
  IF(nRoot.EQ.0)THEN
    IF(inPlane)THEN
      nInter(3)=nInter(3)+1
    END IF
  ELSE IF (nRoot.EQ.1) THEN
    IF(ABS(eta(1)).LT.epsOne)THEN
      xi(1)=eta(1)*(a2(1)-a1(1))+a2(2)-a1(2)
      xi(1)=1.0/xi(1)
      xi(1)=(eta(1)*(a1(3)-a2(3))+a1(4)-a2(4))*xi(1)
      IF(ABS(xi(1)).LT.epsOne)THEN
        t(1)=(xi(1)*eta(1)*BiCoeff(1,1)+xi(1)*BiCoeff(1,2)+eta(1)*BiCoeff(1,3)+LocCoeff(1))*q(1) &
            +(xi(1)*eta(1)*BiCoeff(2,1)+xi(1)*BiCoeff(2,2)+eta(1)*BiCoeff(2,3)+LocCoeff(2))*q(2) &
            +(xi(1)*eta(1)*BiCoeff(3,1)+xi(1)*BiCoeff(3,2)+eta(1)*BiCoeff(3,3)+LocCoeff(3))*q(3)
        IF((t(1).GE.epsZero).AND.(t(1).LE.lenq))THEN
          !WRITE(*,*) ' One Intersection'
          !!WRITE(*,*) ' t ', t(1)
          !WRITE(*,*) ' Intersection at ', lastPartState+t(1)*q
          !isIntersection(1)=.TRUE.
          nInter(3)=nInter(3)+1
        END IF 
      END IF
    END IF
  ELSE 
    IF(ABS(eta(1)).LT.epsOne)THEN
      xi(1)=eta(1)*(a2(1)-a1(1))+a2(2)-a1(2)
      xi(1)=1.0/xi(1)
      xi(1)=(eta(1)*(a1(3)-a2(3))+a1(4)-a2(4))*xi(1)
      !print*,'xi1',xi(1)
      IF(ABS(xi(1)).LT.epsOne)THEN
        t(1)=(xi(1)*eta(1)*BiCoeff(1,1)+xi(1)*BiCoeff(1,2)+eta(1)*BiCoeff(1,3)+LocCoeff(1))*q(1) &
            +(xi(1)*eta(1)*BiCoeff(2,1)+xi(1)*BiCoeff(2,2)+eta(1)*BiCoeff(2,3)+LocCoeff(2))*q(2) &
            +(xi(1)*eta(1)*BiCoeff(3,1)+xi(1)*BiCoeff(3,2)+eta(1)*BiCoeff(3,3)+LocCoeff(3))*q(3)
        !print*,t(1),lenq
        IF((t(1).GE.epsZero).AND.(t(1).LE.lenq))THEN
          !WRITE(*,*) ' One Intersection'
          !WRITE(*,*) ' t ', t(1)
          !WRITE(*,*) ' Intersection at ', lastPartState+t(1)*q
          !isIntersection(1)=.TRUE.
          nInter(3)=nInter(3)+1
        END IF 
      END IF
    END IF
    IF(ABS(eta(2)).LT.epsOne)THEN
      xi(2)=eta(2)*(a2(1)-a1(1))+a2(2)-a1(2)
      xi(2)=1.0/xi(2)
      xi(2)=(eta(2)*(a1(3)-a2(3))+a1(4)-a2(4))*xi(2)
      !print*,'xi2',xi(2)
      IF(ABS(xi(2)).LT.epsOne)THEN
        t(2)=(xi(2)*eta(2)*BiCoeff(1,1)+xi(2)*BiCoeff(1,2)+eta(2)*BiCoeff(1,3)+LocCoeff(1))*q(1) &
            +(xi(2)*eta(2)*BiCoeff(2,1)+xi(2)*BiCoeff(2,2)+eta(2)*BiCoeff(2,3)+LocCoeff(2))*q(2) &
            +(xi(2)*eta(2)*BiCoeff(3,1)+xi(2)*BiCoeff(3,2)+eta(2)*BiCoeff(3,3)+LocCoeff(3))*q(3)
        !print*,t(2),lenq
        IF((t(2).GE.epsZero).AND.(t(2).LE.lenq))THEN
          !WRITE(*,*) ' Second Intersection'
          !WRITE(*,*) ' t ', t(2)
          !WRITE(*,*) ' Intersection at ', lastPartState+t(2)*q
          !isIntersection(2)=.TRUE.
          nInter(3)=nInter(3)+1
        END IF 
      END IF
    END IF
  END IF
  CALL CPU_TIME(tend)   
  time(3)=time(3)+tend-tstart
  
  !print*,iloop,nInter(1),nInter(3)
!  IF(nInter(3).NE.nInter(1))THEN
!    print*,"partpos",PartState
!    print*,"lastpartpos",lastPartState
!    print*,'iloop',iloop
!    EXIT
!    !STOP
!  END IF

  !--------------------------------------------------------------------------------------------------------------------------------
  ! modified ray intersection with patch, Ramsey 2004 paper
  ! instead of 0,1 the interval is -1;1
  !--------------------------------------------------------------------------------------------------------------------------------
  
  ! generate coefficients
  BiCoeff(:,1) = 0.25*( xNode(:,1) -xNode(:,2) +xNode(:,3) -xNode(:,4))
  BiCoeff(:,2) = 0.25*(-xNode(:,1) +xNode(:,2) +xNode(:,3) -xNode(:,4))
  BiCoeff(:,3) = 0.25*(-xNode(:,1) -xNode(:,2) +xNode(:,3) +xNode(:,4))
  BiCoeff(:,4) = 0.25*( xNode(:,1) +xNode(:,2) +xNode(:,3) +xNode(:,4))
  
  CALL CPU_TIME(tstart)
  ! particle vector / not normalized
  q = PartState - lastPartState

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
  !print*,'A,B,C', A,B,C
  CALL QuatricSolver(A,B,C,nRoot,Eta(1),Eta(2))
  !print*,nRoot,Eta
!  IF(iloop.EQ.34)THEN
!    print*,eta
!  END IF

  !print*,'nRoot,v', nRoot,v
  
  inPlane=.FALSE.
  IF(nRoot.EQ.0)THEN
    IF((ABS(A).LT.eps).AND.(ABS(B).LT.eps).AND.(ABS(C).LT.eps))THEN
      inPlane=.TRUE.
    END IF
  END IF
  
  isIntersection=.FALSE.
  IF(nRoot.EQ.0)THEN
    IF(inPlane)THEN
      nInter(4)=nInter(4)+1
    END IF
  ELSE IF (nRoot.EQ.1) THEN
    IF(ABS(eta(1)).LT.epsOne)THEN
     xi(1)=eta(1)*(a2(1)-a1(1))+a2(2)-a1(2)
     xi(1)=1.0/xi(1)
     xi(1)=(eta(1)*(a1(3)-a2(3))+a1(4)-a2(4))*xi(1)
     !print*,'xi(1)',xi(1)
     IF(ABS(xi(1)).LT.epsOne)THEN
       q1=xi(1)*eta(1)*BiCoeff(:,1)+xi(1)*BiCoeff(:,2)+eta(1)*BiCoeff(:,3)+BiCoeff(:,4)-lastPartState
       t(1)=Computet(q1,q)
       IF((t(1).GE.epsZero).AND.(t(1).LE.epsOne))THEN
         !WRITE(*,*) ' One Intersection'
         !WRITE(*,*) ' t ', t(1)
         !WRITE(*,*) ' Intersection at ', lastPartState+t(1)*q
         !isIntersection(1)=.TRUE.
         nInter(4)=nInter(4)+1
       END IF 
     END IF
   END IF
  ELSE 
    IF(ABS(eta(1)).LT.epsOne)THEN
      !xi(1)=eta(1)*a2(1)+a2(2)
      !xi(1)=1.0/xi(1)
      !xi(1)=(-eta(1)*a2(3)-a2(4))*xi(1)
      xi(1)=eta(1)*(a2(1)-a1(1))+a2(2)-a1(2)
      xi(1)=1.0/xi(1)
      xi(1)=(eta(1)*(a1(3)-a2(3))+a1(4)-a2(4))*xi(1)
      IF(ABS(xi(1)).LT.epsOne)THEN
        q1=xi(1)*eta(1)*BiCoeff(:,1)+xi(1)*BiCoeff(:,2)+eta(1)*BiCoeff(:,3)+BiCoeff(:,4)-lastPartState
        t(1)=Computet(q1,q)
        IF((t(1).GE.epsZero).AND.(t(1).LE.epsOne))THEN
          !WRITE(*,*) ' One Intersection'
          !WRITE(*,*) ' t ', t(1)
          !WRITE(*,*) ' Intersection at ', lastPartState+t(1)*q
          !isIntersection(1)=.TRxiE.
          nInter(4)=nInter(4)+1
        END IF 
      END IF
    END IF
    IF(ABS(eta(2)).LT.epsOne)THEN
      xi(2)=eta(2)*a2(1)-eta(2)*a1(1)+a2(2)-a1(2)
      xi(2)=1.0/xi(2)
      xi(2)=(eta(2)*a1(3)-eta(2)*a2(3)+a1(4)-a2(4))*xi(2)
      IF(ABS(xi(2)).LT.epsOne)THEN
        q1=xi(2)*eta(2)*BiCoeff(:,1)+xi(2)*BiCoeff(:,2)+eta(2)*BiCoeff(:,3)+BiCoeff(:,4)-lastPartState
        t(2)=Computet(q1,q)
        IF((t(2).GE.epsZero).AND.(t(2).LE.epsOne))THEN
          !WRITE(*,*) ' Second Intersection'
          !WRITE(*,*) ' t ', t(2)
          !WRITE(*,*) ' Intersection at ', lastPartState+t(2)*q
          !isIntersection(2)=.TRxiE.
          nInter(4)=nInter(4)+1
        END IF 
      END IF
    END IF
  END IF
  CALL CPU_TIME(tend)   
  time(4)=time(4)+tend-tstart

!  IF(nInter(4).NE.nInter(1))THEN
!    print*,"partpos",PartState
!    print*,"lastpartpos",lastPartState
!    print*,'iloop',iloop
!    EXIT
!    !STOP
!  END IF

END DO ! iloop

WRITE(*,*) ' '
WRITE(*,'(A)') '                    Ray-Paper         PIC-Paper           Own-PIC           Own-Ray'
WRITE(*,'(A,5x,I10,8x,I10,8x,I10,8x,I10)') ' Intersections', nInter(1), nInter(2), nInter(3),nInter(4)
WRITE(*,'(A,12x,F12.9,6x,F12.9,6x,F12.9,6x,F12.9)') ' time', time(1), time(2), time(3), time(4)
WRITE(*,'(A,9x,F12.9,6x,F12.9,6x,F12.9,6x,F12.9)') ' time,n1', time(1)/time(1), time(2)/time(1), time(3)/time(1), time(4)/time(1)


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
REAL                    :: eps=1e-12, radicant
!================================================================================================================================

radicant = B*B-4.0*A*C
IF(ABS(a).LT.eps)THEN
  IF(ABS(b).LT.eps)THEN
    nRoot=0
    R1=0.
    R2=0.
  ELSE
    nRoot=1
    R1=-c/b
    R2=0.
  END IF
ELSE
  IF(radicant.LT.0) THEN
    nRoot = 0
    R1=0.
    R2=0.
  ELSE IF (ABS(radicant).LT.eps)THEN
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
