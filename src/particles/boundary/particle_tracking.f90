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
    PartTrajectory=PartState(1:3,iPart) - LastPartPos(1:3,iPart)
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
      END DO ! ilocSide
      IF(alpha.GT.0.)THEN
        IF(alpha+epsilontol.GE.epsilonOne) PartisDone=.TRUE.
        IF(SideID.LT.nBCSides)THEN
          CALL abort(__STAMP__,&
              ' Boundary interaction not implemented for new method.',999,999.)
        END IF
        iInterSect=INT((xi+epsilontol)/1.0)+INT((eta+epsilontol)/1.0)
        IF(iInterSect.GT.0)THEN
          CALL abort(__STAMP__,&
              ' Particle went through edge or node. Not implemented yet.',999,999.)
        ELSE
          dolocSide=.TRUE.
          dolocSide(neighborlocSideID(ilocSide,ElemID))=.FALSE.
          ElemID=neighborElemID(ilocSide,ElemID)
          CALL abort(__STAMP__,&
              ' Particle mapping to neighbor elem not verified!',999,999.)
          EXIT
        END IF
      ELSE
        PartisDone=.TRUE.
      END IF
    END DO
  ELSE
    CYCLE
  END IF
END DO ! iPart

END SUBROUTINE ParticleTracking

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
INTEGER,INTENT(IN)                :: iPart,SideID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: alpha,xi,eta
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: coeffA,coeffB,xInter(3)
REAL                              :: Axz, Bxz
REAL                              :: Ayz, Byz
!===================================================================================================================================

! check if the particle can intersect with the planar plane
! if the normVec point in the opposite direction, cycle
coeffA=DOT_PRODUCT(SideNormVec(1:3,SideID),PartTrajectory)
IF(coeffA.LT.-epsilontol)THEN
  alpha=-1.0
  RETURN
END IF

coeffB=DOT_PRODUCT(SideNormVec(1:3,SideID),LastPartPos(1:3,iPart))

coeffB=SideDistance(SideID)-coeffB

alpha=coeffB/coeffA

IF(alpha.GT.epsilonOne) THEN
  alpha=-1.0
  RETURN
END IF

! compute intersection
xInter(1:3) =LastPartPos(1:3,iPart)+alpha*PartTrajectory

! theoretically, can be computed in advance
Axz = BiLinearCoeff(1,2,SideID) - BiLinearCoeff(3,2,SideID)
Bxz = BiLinearCoeff(1,3,SideID) - BiLinearCoeff(3,3,SideID)

Ayz = BiLinearCoeff(2,2,SideID) - BiLinearCoeff(3,2,SideID)
Byz = BiLinearCoeff(2,3,SideID) - BiLinearCoeff(3,3,SideID)

xi = Ayz*Bxz - Axz
xi = (BiLinearCoeff(1,4,SideID)-LastPartPos(1,iPart) - BiLinearCoeff(2,4,SideID)+LastPartPos(2,SideID))/xi

eta = (xi*-1.0*Axz+BiLinearCoeff(3,4,SideID)-LastPartPos(3,iPart) - BiLinearCoeff(1,4,SideID)+LastPartPos(1,SideID)) / Bxz

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
INTEGER                           :: whichInter,Inter1,Inter2,nRoot
!===================================================================================================================================

a1(1)= BilinearCoeff(1,1,SideID)*PartTrajectory(3) - BilinearCoeff(3,1,SideID)*PartTrajectory(1)
a1(2)= BilinearCoeff(1,2,SideID)*PartTrajectory(3) - BilinearCoeff(3,2,SideID)*PartTrajectory(1)
a1(3)= BilinearCoeff(1,3,SideID)*PartTrajectory(3) - BilinearCoeff(3,3,SideID)*PartTrajectory(1)
a1(4)= (BilinearCoeff(1,4,SideID)-LastPartPos(1,iPart))*PartTrajectory(3) &
     - (BilinearCoeff(3,4,SideID)-LastPartPos(3,iPart))*PartTrajectory(1)

a2(1)= BilinearCoeff(2,1,SideID)*PartTrajectory(3) - BilinearCoeff(3,1,SideID)*PartTrajectory(2)
a2(2)= BilinearCoeff(2,2,SideID)*PartTrajectory(3) - BilinearCoeff(3,2,SideID)*PartTrajectory(2)
a2(3)= BilinearCoeff(2,3,SideID)*PartTrajectory(3) - BilinearCoeff(3,3,SideID)*PartTrajectory(2)
a2(4)= (BilinearCoeff(2,4,SideID)-LastPartPos(2,iPart))*PartTrajectory(3) &
     - (BilinearCoeff(3,4,SideID)-LastPartPos(3,iPart))*PartTrajectory(2)

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
  alpha=-1.
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
     IF((t(1).GE.-epsilontol).AND.(t(1).LE.epsilonOne))THEN
       alpha=t(1)
       xitild=xi(1)
       etatild=eta(1)
       RETURN
     ELSE
       alpha=-1
       RETURN
     END IF
   END IF
 END IF
ELSE 
  t=-1.0
  IF(ABS(eta(1)).LT.epsilonOne)THEN
    xi(1)=eta(1)*(a2(1)-a1(1))+a2(2)-a1(2)
    xi(1)=1.0/xi(1)
    xi(1)=(eta(1)*(a1(3)-a2(3))+a1(4)-a2(4))*xi(1)
    IF(ABS(xi(1)).LT.epsilonOne)THEN
      ! q1=xi(1)*eta(1)*BilinearCoeff(:,1)+xi(1)*BilinearCoeff(:,2)+eta(1)*BilinearCoeff(:,3)+BilinearCoeff(:,4)-lastPartState
      !  WRITE(*,*) ' t ', t(2)
      !  WRITE(*,*) ' Intersection at ', lastPartState+t(2)*q
      t(1)=ComputeSurfaceDistance(xi(1),eta(1),PartTrajectory,iPart,SideID)
      IF((t(1).GE.epsilontol).AND.(t(1).LT.epsilonOne))THEN
        alpha=t(1)
        xitild=xi(1)
        etatild=eta(1)
      ELSE
        alpha=-1
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
      IF((t(2).GE.epsilontol).AND.(t(2).LT.epsilonOne))THEN
        alpha=t(2)
        xitild=xi(2)
        etatild=eta(2)
      ELSE
        alpha=-1
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
  !IF(SideID.LT.nBCSides)THEN
  !  IF(t(1).LT.t(2))THEN
  !    IF((t(1).GE.epsilontol).AND.(t(1).LT.epsilonOne))THEN
  !      alpha=t(1)
  !      xitild=xi(1)
  !      etatild=eta(1)
  !    ELSE
  !      alpha=-1.0
  !    END IF
  !  ELSE
  !    IF((t(2).GE.epsilontol).AND.(t(2).LT.epsilonOne))THEN
  !      alpha=t(2)
  !      xitild=xi(2)
  !      etatild=eta(2)
  !    ELSE
  !      alpha=-1.0
  !    END IF
  !  END IF
  !ELSE ! no BC side , assume that particle can re-enter cell
  !  IF((t(1).GE.epsilontol).AND.(t(1).LE.epsilonOne).AND.(t(2).GE.epsilontol).AND.(t(2).LE.epsilonOne))THEN
  !    alpha=-1 ! partilce remains in cell
  !  ELSE IF((t(1).GE.epsilontol).AND.(t(1).LE.epsilonOne).AND.(t(2).LT.epsilontol).OR.(t(2).GT.epsilonOne))THEN
  !    alpha=t(1)
  !    xitild=xi(1)
  !    etatild=eta(1)
  !  ELSE IF((t(2).GE.epsilontol).AND.(t(2).LE.epsilonOne).AND.(t(1).LT.epsilontol).OR.(t(1).GT.epsilonOne))THEN
  !    alpha=t(2)
  !    xitild=xi(2)
  !    etatild=eta(2)
  !  ELSE
  !    alpha=-1.0 ! particle move a distance which is too short
  !  END IF
  !END IF 
  ! alternative version
  !Inter1=INT((t(1)+epsilontol)/(1.0+2*epsilontol))
  !Inter2=INT((t(2)+epsilontol)/(1.0+2*epsilontol))
  Inter1=INT((t(1))/(epsilonOne))
  Inter2=INT((t(2))/(epsilonOne))
  whichInter=Inter1+Inter2
  IF(SideID.LT.nBCSides)THEN
    SELECT CASE(whichInter)
    CASE(0)
      IF(t(1).LT.t(2))THEN
        alpha=t(1)
        xitild=xi(1)
        etatild=eta(1)
      ELSE
        alpha=t(2)
        xitild=xi(2)
        etatild=eta(2)
      END IF
    CASE(1)
    IF(Inter1.NE.0)THEN
      alpha=t(1)
      xitild=xi(1)
      etatild=eta(1)
    ELSE
      alpha=t(2)
      xitild=xi(2)
      etatild=eta(2)
    END IF
    CASE DEFAULT
      alpha=-1.0 ! no Intersection
    END SELECT
  ELSE ! no BC Side
    SELECT CASE(whichInter)
    CASE(0)
      alpha=-1.0 ! particle leace and reenter cell
    CASE(1)
      IF(Inter1.NE.0)THEN
        alpha=t(1)
        xitild=xi(1)
        etatild=eta(1)
      ELSE
        alpha=t(2)
        xitild=xi(2)
        etatild=eta(2)
      END IF
    CASE DEFAULT
      alpha=-1.0 ! particle move not a enough long distance
    END SELECT
  END IF
END IF

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
             -lastPartPos(1,iPart)
  t = t/ PartTrajectory(1)-epsilontol 
ELSE IF(ABS(PartTrajectory(2)).GE.ABS(PartTrajectory(3)))THEN
  t =xi*eta*BilinearCoeff(2,1,SideID)+xi*BilinearCoeff(2,2,SideID)+eta*BilinearCoeff(2,3,SideID)+BilinearCoeff(2,4,SideID) &
             -lastPartPos(2,iPart)
  t = t/ PartTrajectory(2)-epsilontol 
ELSE
  t =xi*eta*BilinearCoeff(3,1,SideID)+xi*BilinearCoeff(3,2,SideID)+eta*BilinearCoeff(3,3,SideID)+BilinearCoeff(3,4,SideID) &
             -lastPartPos(3,iPart)
  t = t/ PartTrajectory(3)-epsilontol 
END IF

ComputeSurfaceDistance=t

END FUNCTION ComputeSurfaceDistance

END MODULE MOD_Particle_Tracking
