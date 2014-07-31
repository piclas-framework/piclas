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
IMPLICIT NONE
USE  Particle_Vars,               ONLY: PDM
USE MOD_Particle_Vars,            ONLY: PartState,LastPartPos
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart
INTEGER                       :: ilocSide,SideID
LOGICAL                       :: PartisDone,hitBC
REAL                          :: alpha(1:6),xietaIntersect(1:2,1:6)
REAL                          :: PartTrajectory(1:3)
!===================================================================================================================================

DO iPart=1,PDM%ParticleVecLength
  IF(PDM%ParticleInside)THEn
    PartisDone=.FALSE.
    ElemID = PEM%lastElement(i)
    PartTrajectory=PartState(1:3,iPart) - LastPartPos(:,iPart)
    ! track particle vector until the final particle position is achieved
    alphaIntersect=-1.
    DO WHILE (.NOT.PartisDone)
      DO ilocSide=1,6
        SideID=ElemToSide(E2S_SIDE_ID,ilocSide,iElem) 
        IF(SideIsPlanar(SideID))THEN
          CALL ComputePlanarIntersection(PartTrajectory,iPart,SideID,ElemID,alpha(ilocSide),xietaIntersect(1,ilocSide) &
                                        ,XiEtaIntersect(2,ilocSide))
        ELSE
          CALL ComputeBiLinearIntersection(PartTrajectory,iPart,SideID,ElemID,alpha(ilocSide),xietaIntersect(1,ilocSide),
                                        ,XiEtaIntersect(2,ilocSide))

        END IF
      END DO ! ilocSide

    IF(hitBC) CYCLE
    END DO
    PEM%Element(i) = ElemID
  ELSE
    CYCLE
  END IF
END DO ! iPart

END SUBROUTINE ParticleTracking

SUBROUTINE ComputePlanarIntersection(PartTrajectory,iPart,SideID,ElemID,alpha,xi,eta)
!===================================================================================================================================
! Compute the Intersection with planar surface
!===================================================================================================================================
! MODULES
IMPLICIT NONE
USE MOD_Particle_Vars,           ONLY:LastPartPos
USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonbilinear,BiLinearCoeff, SideNormVec,epsilontol,SideDistance
!USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonOne,SideIsPlanar,BiLinearCoeff,SideNormVec
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
INTEGER,INTENT(IN)                :: iPart,SideID,ElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: alpha,xi,eta
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: coeffA,coeffB,xInter(3)
!===================================================================================================================================

! check if the particle can intersect with the planar plane
! if the normVec point in the opposite direction, cycle
coeffA=DOT_PRODUCT(SideNormVec(1:3,SideID),PartTrajectory)
IF(coeffA.LT.-epsilontol)THEN
  alpha=-1.
  RETURN
END IF

coeffB=DOT_PRODUCT(SideNormVec(1:3,SideID),LastPartPos(1:3,iPart))

coeffB=SideDistance(SideID)-coeffB

alpha=coeffB/coeffA

IF(alpha.GT.epsilonOne) RETURN

! compute intersection
xInter(1:3) =LastPartPos(1:3,iPart)+alpha*PartTrajectory



END SUBROUTINE ComputePlanarIntersection
SUBROUTINE ComputeBiLinearIntersection(PartTrajectory,iPart,SideID,ElemID,alpha,xi,eta)
!===================================================================================================================================
! Compute the Intersection with planar surface
!===================================================================================================================================
! MODULES
IMPLICIT NONE
USE MOD_Particle_Vars,           ONLY:LastPartPos
USE MOD_Mesh_Vars,               ONLY:nBCSides
!USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonOne,SideIsPlanar,BiLinearCoeff,SideNormVec
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
INTEGER,INTENT(IN)                :: iPart,SideID,ElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: alpha,xi,eta
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================


END SUBROUTINE ComputeBiLinearIntersection

END MODULE MOD_Particle_Tracking
