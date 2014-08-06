#include "boltzplatz.h"

MODULE MOD_Particle_Tracking
!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC

INTERFACE ParticleTracking
  MODULE PROCEDURE ParticleTracking
END INTERFACE

PUBLIC::ParticleTracking
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================

CONTAINS

SUBROUTINE ParticleTracking()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals,                  ONLY:abort
USE MOD_Mesh_Vars,                ONLY:ElemToSide,nBCSides
USE MOD_Particle_Vars,            ONLY:PEM,PDM
USE MOD_Particle_Vars,            ONLY:PartState,LastPartPos
USE MOD_Particle_Surfaces_Vars,   ONLY:epsilontol,SideIsPlanar,epsilonOne,neighborElemID,neighborlocSideID
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart,ElemID
INTEGER                       :: ilocSide,SideID
INTEGER                       :: iInterSect
LOGICAL                       :: PartisDone,dolocSide(1:6)
!REAL                          :: alpha(1:6),xietaIntersect(1:2,1:6)
REAL                          :: alpha,xi,eta!xietaIntersect(1:2,1:6)
REAL                          :: PartTrajectory(1:3)
!===================================================================================================================================

DO iPart=1,PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart))THEN
    PartisDone=.FALSE.
    ElemID = PEM%lastElement(iPart)
   !Element = PEM%lastElement(i)
    print*,ElemID
    PartTrajectory=PartState(iPart,1:3) - LastPartPos(iPart,1:3)
    ! track particle vector until the final particle position is achieved
    alpha=-1.
    dolocSide=.TRUE.
    DO WHILE (.NOT.PartisDone)
      DO ilocSide=1,6
        IF(.NOT.dolocSide(ilocSide)) CYCLE
        SideID=ElemToSide(E2S_SIDE_ID,ilocSide,ElemID) 
        IF(SideIsPlanar(SideID))THEN
          !CALL ComputePlanarIntersection(PartTrajectory,iPart,SideID,ElemID,alpha,xietaIntersect(1,ilocSide) &
          !                              ,XiEtaIntersect(2,ilocSide))
          CALL ComputePlanarIntersection(PartTrajectory,iPart,SideID,alpha,xi,eta)
        ELSE
          CALL ComputeBiLinearIntersection(PartTrajectory,iPart,SideID,alpha,xi,eta)
          !CALL ComputeBiLinearIntersection(PartTrajectory,iPart,SideID,ElemID,alpha,xietaIntersect(1,ilocSide) &
          !                              ,XiEtaIntersect(2,ilocSide))
        END IF
        !print*,ilocSide,alpha
        print*,'ilocSide,alpha,xi,eta',ilocSide,alpha,xi,eta
        !print*,'neighborElemID',neighborElemID(ilocSide,ElemID)
        ! check after each side if particle went through checked side
        IF(alpha.GT.epsilontol)THEN ! or minus epsilontol
          !IF(alpha+epsilontol.GE.epsilonOne) PartisDone=.TRUE.
          IF(SideID.LE.nBCSides)THEN
            CALL abort(__STAMP__,&
                ' Boundary interaction not implemented for new method.',999,999.)
          END IF
          iInterSect=INT((ABS(xi)-epsilontol)/1.0)+INT((ABS(eta)-epsilontol)/1.0)
          IF(iInterSect.GT.0)THEN
            CALL abort(__STAMP__,&
                ' Particle went through edge or node. Not implemented yet.',999,999.)
          ELSE
            dolocSide=.TRUE.
            dolocSide(neighborlocSideID(ilocSide,ElemID))=.FALSE.
            ElemID=neighborElemID(ilocSide,ElemID)
            !print*,'new particle positon',ElemID
!            CALL abort(__STAMP__,&
!                ' Particle mapping to neighbor elem not verified!',999,999.)
            EXIT
          END IF ! iInteSect
        END IF
      END DO ! ilocSide
      !sop
      !read*
      ! no intersection found
      IF(alpha.EQ.-1.0)THEN
        PEM%Element(iPart) = ElemID
        PartisDone=.TRUE.
      END IF
    END DO ! PartisDone=.FALSE.
  END IF ! Part inside
END DO ! iPart

END SUBROUTINE ParticleTracking

!SUBROUTINE ComputePlanarIntersection(PartTrajectory,iPart,SideID,alpha,xi,eta)
!!==================================================================================================================================
!! Compute the Intersection with planar surface
!!==================================================================================================================================
!! MODULES
!USE MOD_Particle_Vars,           ONLY:LastPartPos
!USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonbilinear,BiLinearCoeff, SideNormVec,epsilontol,SideDistance,epsilonOne
!!USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonOne,SideIsPlanar,BiLinearCoeff,SideNormVec
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!! INPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
!INTEGER,INTENT(IN)                :: iPart,SideID!,ElemID
!!----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!REAL,INTENT(OUT)                  :: alpha,xi,eta
!!----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!REAL                              :: coeffA,coeffB,xInter(3)
!REAL                              :: Axz, Bxz, Cxz
!REAL                              :: Ayz, Byz, Cyz
!!==================================================================================================================================
!
!! set alpha to minus 1, asume no intersection
!alpha=-1.0
!xi=-2.
!eta=-2.
!
!! check if the particle can intersect with the planar plane
!! if the normVec point in the opposite direction, cycle
!coeffA=DOT_PRODUCT(SideNormVec(1:3,SideID),PartTrajectory)
!!print*,'coeffA',coeffA
!IF(ABS(coeffA).LT.+epsilontol)THEN
!  ! particle tangential to surface ==> no interesection, particle remains in element
!  RETURN
!END IF
!
!! distance of plane fromn origion minus trajectory start point times normal vector of side
!coeffB=SideDistance(SideID)-DOT_PRODUCT(SideNormVec(1:3,SideID),LastPartPos(iPart,1:3))
!
!alpha=coeffB/coeffA
!!print*,'coeffB',coeffB
!!print*,'alpha',alpha
!!read*
!
!IF((alpha.GT.epsilonOne).OR.(alpha.LT.epsilontol))THEN
!  alpha=-1.0
!  RETURN
!END IF
!
!! compute intersection
!xInter(1:3) =LastPartPos(iPart,1:3)+alpha*PartTrajectory(1:3)
!
!!! theoretically, can be computed in advance
!Axz = BiLinearCoeff(1,2,SideID) - BiLinearCoeff(3,2,SideID)
!Bxz = BiLinearCoeff(1,3,SideID) - BiLinearCoeff(3,3,SideID)
!Cxz = xInter(1) - BiLinearCoeff(1,4,SideID) - xInter(3) + BiLinearCoeff(3,4,SideID)
!
!Ayz = BiLinearCoeff(2,2,SideID) - BiLinearCoeff(3,2,SideID)
!Byz = BiLinearCoeff(2,3,SideID) - BiLinearCoeff(3,3,SideID)
!Cyz = xInter(2) - BiLinearCoeff(2,4,SideID) - xInter(3) + BiLinearCoeff(3,4,SideID)
!
!print*,'Bxz,Byz',Bxz,Byz
!
!IF(ABS(Bxz).LT.epsilontol)THEN
!  xi = Axz + Bxz*Ayz/Byz
!  ! check denominator
!  xi = (Cxz - Bxz*Ayz/Byz*Cyz)/xi
!ELSE
!  xi = Ayz + Byz*Axz/Bxz
!  ! check denominator
!  xi = (Cyz - Byz*Axz/Bxz*Cxz)/xi
!END IF
!
!IF(ABS(xi).GT.epsilonOne) THEN 
!  ! xi outside of possible range
!  alpha=-1.0
!  RETURN
!END IF
!
!eta = Bxz+Byz
!eta = (Cxz+Cyz - (Axz+Ayz)*xi) / eta
!
!!print*,'xi,eta',xi,eta
!IF(ABS(eta).GT.epsilonOne) THEN 
!  ! eta outside of possible range
!  alpha=-1.0
!  RETURN
!END IF
!
!! here, eta,xi,alpha are computed
!
!END SUBROUTINE ComputePlanarIntersection

SUBROUTINE ComputePlanarIntersection(PartTrajectory,iPart,SideID,alpha,xi,eta)
!===================================================================================================================================
! Compute the Intersection with planar surface
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars,           ONLY:LastPartPos
USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonbilinear,BiLinearCoeff, SideNormVec,epsilontol,SideDistance,epsilonOne
!USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonOne,SideIsPlanar,BiLinearCoeff,SideNormVec
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
INTEGER,INTENT(IN)                :: iPart,SideID!,ElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: alpha,xi,eta
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(2:4)                 :: a1,a2  ! array dimension from 2:4 according to bi-linear surface
REAL                                :: coeffA,coeffB
!===================================================================================================================================

! set alpha to minus 1, asume no intersection
!print*,PartTrajectory
alpha=-1.0
xi=-2.
eta=-2.

! compute distance of lastPartPos with planar plane

coeffA=DOT_PRODUCT(SideNormVec(1:3,SideID),PartTrajectory)
!print*,'coeffA',coeffA
!read*
! corresponding to particle starting in plane
! interaction should be computed in last step
IF(ABS(coeffA).LT.+epsilontol)THEN 
  RETURN
END IF
! distance of plane fromn origion minus trajectory start point times normal vector of side
coeffB=SideDistance(SideID)-DOT_PRODUCT(SideNormVec(1:3,SideID),LastPartPos(iPart,1:3))

alpha=coeffB/coeffA
!!print*,'coeffB',coeffB
!!print*,'alpha',alpha
!!read*

IF((alpha.GT.epsilonOne).OR.(alpha.LT.-epsilontol))THEN
  alpha=-1.0
  RETURN
END IF


a1(2)= BilinearCoeff(1,2,SideID)*PartTrajectory(3) - BilinearCoeff(3,2,SideID)*PartTrajectory(1)
a1(3)= BilinearCoeff(1,3,SideID)*PartTrajectory(3) - BilinearCoeff(3,3,SideID)*PartTrajectory(1)
a1(4)= (BilinearCoeff(1,4,SideID)-LastPartPos(iPart,1))*PartTrajectory(3) &
     - (BilinearCoeff(3,4,SideID)-LastPartPos(iPart,3))*PartTrajectory(1)

a2(2)= BilinearCoeff(2,2,SideID)*PartTrajectory(3) - BilinearCoeff(3,2,SideID)*PartTrajectory(2)
a2(3)= BilinearCoeff(2,3,SideID)*PartTrajectory(3) - BilinearCoeff(3,3,SideID)*PartTrajectory(2)
a2(4)= (BilinearCoeff(2,4,SideID)-LastPartPos(iPart,2))*PartTrajectory(3) &
     - (BilinearCoeff(3,4,SideID)-LastPartPos(iPart,3))*PartTrajectory(2)

!print*,'a23,a13',a2(3),a1(3)
!print*,'a22,a12',a2(2),a1(2)

!! caution with accuracy
IF(ABS(a2(3)).LT.epsilontol)THEN ! term c is close to zero ==> eta is zero
  eta=0.
  IF(ABS(a2(2)).LT.epsilontol)THEN
    xi=0.
  ELSE
    ! compute xi
    xi=a1(2)-a2(2)
    xi=1.0/xi
    xi=(a2(4)-a1(4))*xi
  END IF
!  IF(ABS(xi).GT.epsilonOne)THEN
!    RETURN
!  END IF
ELSE ! a2(3) not zero
  IF(ABS(a2(2)).LT.epsilontol)THEN
    xi=0.
    eta=a1(3)-a2(3)
    eta=1.0/eta
    eta=(a2(4)-a1(4))*eta
  ELSE
    xi = a1(2) - a1(3)*a2(2)/a2(3)
    xi = 1.0/xi
    xi = (-a1(4)-a1(3)*a2(4)/a2(3))*xi
    ! check distance of xi 
  !  IF(ABS(xi).GT.epsilonOne)THEN
  !    RETURN
  !  END IF
    ! compute eta
    eta=a1(3)-a2(3)
    eta=1.0/eta
    eta=((a2(2)-a1(2))*xi+a2(4)-a1(4))*eta
  END IF
END IF

!xi = a1(2) - a1(3)*a2(2)/a2(3)
!xi = 1.0/xi
!xi = (-a1(4)-a1(3)*a2(4)/a2(3))*xi
!! check distance of xi 
IF(ABS(xi).GT.epsilonOne)THEN
  alpha=-1.0
  RETURN
END IF
!! compute eta
!eta=a1(3)-a2(3)
!eta=1.0/eta
!eta=((a2(2)-a1(2))*xi+a2(4)-a1(4))*eta
!
IF(ABS(eta).GT.epsilonOne)THEN
  alpha=-1.0
  RETURN
END IF

!! compute distance with intersection
!IF((ABS(PartTrajectory(1)).GE.ABS(PartTrajectory(2))).AND.(ABS(PartTrajectory(1)).GT.ABS(PartTrajectory(3))))THEN
!  alpha =xi*BilinearCoeff(1,2,SideID)+eta*BilinearCoeff(1,3,SideID)+BilinearCoeff(1,4,SideID) -lastPartPos(iPart,1)
!  alpha = alpha/ PartTrajectory(1)
!ELSE IF(ABS(PartTrajectory(2)).GE.ABS(PartTrajectory(3)))THEN
!  alpha =xi*BilinearCoeff(2,2,SideID)+eta*BilinearCoeff(2,3,SideID)+BilinearCoeff(2,4,SideID) -lastPartPos(iPart,2)
!  alpha = alpha/ PartTrajectory(2)
!ELSE
!  alpha =xi*BilinearCoeff(3,2,SideID)+eta*BilinearCoeff(3,3,SideID)+BilinearCoeff(3,4,SideID) -lastPartPos(iPart,3)
!  alpha = alpha/ PartTrajectory(3)
!END IF
!
!IF((alpha.LT.epsilontol).OR.(alpha.GT.epsilonOne)) alpha=-1.0

END SUBROUTINE ComputePlanarIntersection

SUBROUTINE ComputeBiLinearIntersection(PartTrajectory,iPart,SideID,alpha,xitild,etatild)
!===================================================================================================================================
! Compute the Intersection with planar surface
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars,           ONLY:LastPartPos
USE MOD_Mesh_Vars,               ONLY:nBCSides
USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonbilinear,BiLinearCoeff, epsilontol,epsilonOne
!USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonOne,SideIsPlanar,BiLinearCoeff,SideNormVec
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
INTEGER,INTENT(IN)                :: iPart,SideID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: alpha,xitild,etatild
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(4)                 :: a1,a2
REAL                              :: A,B,C
REAL                              :: xi(2),eta(2),t(2), q1(3)
INTEGER                           :: nInter,nRoot
!===================================================================================================================================

! set alpha to minus one // no interesction
alpha=-1.0
xitild=-2.0
etatild=-2.0

a1(1)= BilinearCoeff(1,1,SideID)*PartTrajectory(3) - BilinearCoeff(3,1,SideID)*PartTrajectory(1)
a1(2)= BilinearCoeff(1,2,SideID)*PartTrajectory(3) - BilinearCoeff(3,2,SideID)*PartTrajectory(1)
a1(3)= BilinearCoeff(1,3,SideID)*PartTrajectory(3) - BilinearCoeff(3,3,SideID)*PartTrajectory(1)
a1(4)= (BilinearCoeff(1,4,SideID)-LastPartPos(iPart,1))*PartTrajectory(3) &
     - (BilinearCoeff(3,4,SideID)-LastPartPos(iPart,3))*PartTrajectory(1)

a2(1)= BilinearCoeff(2,1,SideID)*PartTrajectory(3) - BilinearCoeff(3,1,SideID)*PartTrajectory(2)
a2(2)= BilinearCoeff(2,2,SideID)*PartTrajectory(3) - BilinearCoeff(3,2,SideID)*PartTrajectory(2)
a2(3)= BilinearCoeff(2,3,SideID)*PartTrajectory(3) - BilinearCoeff(3,3,SideID)*PartTrajectory(2)
a2(4)= (BilinearCoeff(2,4,SideID)-LastPartPos(iPart,2))*PartTrajectory(3) &
     - (BilinearCoeff(3,4,SideID)-LastPartPos(iPart,3))*PartTrajectory(2)

A = a2(1)*a1(3)-a1(1)*a2(3)
B = a2(1)*a1(4)-a1(1)*a2(4)+a2(2)*a1(3)-a1(2)*a2(3)
C = a1(4)*a2(2)-a1(2)*a2(4)
!print*,'A,B,C', A,B,C
CALL QuatricSolver(A,B,C,nRoot,Eta(1),Eta(2))
!print*,nRoot,Eta
!  IF(iloop.EQ.34)THEN
!    print*,eta
!  END IF

IF(nRoot.EQ.0)THEN
  RETURN
END IF

IF (nRoot.EQ.1) THEN
  IF(ABS(eta(1)).LT.epsilonOne)THEN
    xi(1)=eta(1)*(a2(1)-a1(1))+a2(2)-a1(2)
    xi(1)=1.0/xi(1)
    xi(1)=(eta(1)*(a1(3)-a2(3))+a1(4)-a2(4))*xi(1)
    IF(ABS(xi(1)).LT.epsilonOne)THEN
      !q1=xi(1)*eta(1)*BilinearCoeff(:,1)+xi(1)*BilinearCoeff(:,2)+eta(1)*BilinearCoeff(:,3)+BilinearCoeff(:,4)-lastPartState
      t(1)=ComputeSurfaceDistance(xi(1),eta(1),PartTrajectory,iPart,SideID)
      IF((t(1).GE.+epsilontol).AND.(t(1).LE.epsilonOne))THEN
        alpha=t(1)
        xitild=xi(1)
        etatild=eta(1)
        RETURN
      ELSE ! t is not in range
        RETURN
      END IF
    ELSE ! xi not in range
      RETURN
    END IF ! xi .lt. epsilonOne
  ELSE ! eta not in reange
    RETURN 
  END IF ! eta .lt. epsilonOne
ELSE 
  nInter=0
  IF(ABS(eta(1)).LT.epsilonOne)THEN
    xi(1)=eta(1)*(a2(1)-a1(1))+a2(2)-a1(2)
    xi(1)=1.0/xi(1)
    xi(1)=(eta(1)*(a1(3)-a2(3))+a1(4)-a2(4))*xi(1)
    IF(ABS(xi(1)).LT.epsilonOne)THEN
      ! q1=xi(1)*eta(1)*BilinearCoeff(:,1)+xi(1)*BilinearCoeff(:,2)+eta(1)*BilinearCoeff(:,3)+BilinearCoeff(:,4)-lastPartState
      !  WRITE(*,*) ' t ', t(2)
      !  WRITE(*,*) ' Intersection at ', lastPartState+t(2)*q
      t(1)=ComputeSurfaceDistance(xi(1),eta(1),PartTrajectory,iPart,SideID)
      IF((t(1).LT.epsilontol).OR.(t(1).GT.epsilonOne))THEN
        t(1)=-2.0
      ELSE
        nInter=nInter+1
      END IF
!      IF((t(1).LT.epsilontol).AND.(t(1).GT.epsilonOne))THEN
!        t(1)=-1
!      END IF
    END IF
  END IF
  IF(ABS(eta(2)).LT.epsilonOne)THEN
    xi(2)=eta(2)*a2(1)-eta(2)*a1(1)+a2(2)-a1(2)
    xi(2)=1.0/xi(2)
    xi(2)=(eta(2)*a1(3)-eta(2)*a2(3)+a1(4)-a2(4))*xi(2)
    IF(ABS(xi(2)).LT.epsilonOne)THEN
      ! q1=xi(2)*eta(2)*BilinearCoeff(:,1)+xi(2)*BilinearCoeff(:,2)+eta(2)*BilinearCoeff(:,3)+BilinearCoeff(:,4)-lastPartState
      t(2)=ComputeSurfaceDistance(xi(2),eta(2),PartTrajectory,iPart,SideID)
      IF((t(2).LT.epsilontol).OR.(t(2).GT.epsilonOne))THEN
        t(2)=-2.0
      ELSE
        nInter=nInter+1
      END IF
!      IF((t(2).LT.epsilontol).AND.(t(2).GT.epsilonOne))THEN
!        t(2)=-1
!      END IF
      !IF((t(2).GE.epsZero).AND.(t(2).LE.epsOne))THEN
      !!  WRITE(*,*) ' Second Intersection'
      !!  WRITE(*,*) ' t ', t(2)
      !!  WRITE(*,*) ' Intersection at ', lastPartState+t(2)*q
      !END IF 
    END IF
  END IF
  ! if no intersection, return
  IF(nInter.EQ.0) RETURN
  IF(SideID.LE.nBCSides)THEN
    IF(ABS(t(1)).LT.ABS(t(2)))THEN
      alpha=t(1)
      xitild=xi(1)
      etatild=eta(1)
    ELSE
      alpha=t(2)
      xitild=xi(2)
      etatild=eta(2)
    END IF
  ELSE ! no BC Side
    ! if two intersections, return, particle re-enters element
    IF(nInter.EQ.2) RETURN
    IF(ABS(t(1)).LT.ABS(t(2)))THEN
      alpha=t(1)
      xitild=xi(1)
      etatild=eta(1)
    ELSE
      alpha=t(2)
      xitild=xi(2)
      etatild=eta(2)
    END IF
  END IF ! SideID.LT.nCBSides
END IF ! nRoot

END SUBROUTINE ComputeBiLinearIntersection

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

FUNCTION ComputeSurfaceDistance(xi,eta,PartTrajectory,iPart,SideID)
!================================================================================================================================
! compute the required vector length to intersection
!================================================================================================================================
USE MOD_Particle_Surfaces_Vars,   ONLY:epsilontol,BiLinearCoeff
USE MOD_Particle_Vars,            ONLY:PartState,LastPartPos
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(3),INTENT(IN)         :: PartTrajectory
REAL,INTENT(IN)                      :: xi,eta
INTEGER,INTENT(IN)                   :: iPart,SideID
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                                 :: ComputeSurfaceDistance
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: t
!================================================================================================================================

IF((ABS(PartTrajectory(1)).GE.ABS(PartTrajectory(2))).AND.(ABS(PartTrajectory(1)).GT.ABS(PartTrajectory(3))))THEN
  t =xi*eta*BiLinearCoeff(1,1,SideID)+xi*BilinearCoeff(1,2,SideID)+eta*BilinearCoeff(1,3,SideID)+BilinearCoeff(1,4,SideID) &
             -lastPartPos(iPart,1)
  t = t/ PartTrajectory(1)-epsilontol 
ELSE IF(ABS(PartTrajectory(2)).GE.ABS(PartTrajectory(3)))THEN
  t =xi*eta*BilinearCoeff(2,1,SideID)+xi*BilinearCoeff(2,2,SideID)+eta*BilinearCoeff(2,3,SideID)+BilinearCoeff(2,4,SideID) &
             -lastPartPos(iPart,2)
  t = t/ PartTrajectory(2)-epsilontol 
ELSE
  t =xi*eta*BilinearCoeff(3,1,SideID)+xi*BilinearCoeff(3,2,SideID)+eta*BilinearCoeff(3,3,SideID)+BilinearCoeff(3,4,SideID) &
             -lastPartPos(iPart,3)
  t = t/ PartTrajectory(3)-epsilontol 
END IF

ComputeSurfaceDistance=t

END FUNCTION ComputeSurfaceDistance

END MODULE MOD_Particle_Tracking
