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

MODULE MOD_Particle_Tracking
!===================================================================================================================================
!> General routines concerning particle tracking
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC::PerformTracking
PUBLIC::ParticleInsideCheck
PUBLIC::ParticleSanityCheck
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================

CONTAINS

SUBROUTINE PerformTracking()
!===================================================================================================================================
!> Routine called from the timedisc to call the selected tracking routine
!===================================================================================================================================
! MODULES
USE MOD_Globals                  ,ONLY: abort,CROSS
USE MOD_Particle_Tracking_vars   ,ONLY: TrackingMethod
USE MOD_Particle_Tracing         ,ONLY: ParticleTracing
USE MOD_Particle_RefTracking     ,ONLY: ParticleRefTracking
USE MOD_Particle_TriaTracking    ,ONLY: ParticleTriaTracking
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

#ifdef EXTRAE
CALL extrae_eventandcounters(int(9000001), int8(50))
#endif /*EXTRAE*/
SELECT CASE(TrackingMethod)
CASE(REFMAPPING)
  CALL ParticleRefTracking()
CASE(TRACING)
  CALL ParticleTracing()
CASE(TRIATRACKING)
  CALL ParticleTriaTracking()
CASE DEFAULT
  CALL abort(__STAMP__,'TrackingMethod not implemented! TrackingMethod =',IntInfoOpt=TrackingMethod)
END SELECT
#ifdef EXTRAE
CALL extrae_eventandcounters(int(9000001), int8(0))
#endif /*EXTRAE*/

END SUBROUTINE PerformTracking


LOGICAL FUNCTION ParticleInsideCheck(Position,iPart,GlobalElemID)
!===================================================================================================================================
!> Checks if the position is inside the element with the appropriate routine depending on the TrackingMethod
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
USE MOD_Particle_Localization   ,ONLY: PartInElemCheck
USE MOD_Particle_Mesh_Tools     ,ONLY: ParticleInsideQuad3D
USE MOD_Eval_xyz                ,ONLY: GetPositionInRefElem
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
USE MOD_Particle_Vars           ,ONLY: PartPosRef
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: Position(3)
INTEGER, INTENT(IN)             :: iPart,GlobalElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

ParticleInsideCheck = .FALSE.

SELECT CASE(TrackingMethod)
CASE(REFMAPPING)
  CALL GetPositionInRefElem(Position,PartPosRef(1:3,iPart),GlobalElemID)
  IF (MAXVAL(ABS(PartPosRef(1:3,iPart))).LT.1.0) ParticleInsideCheck=.TRUE.
CASE(TRACING)
  CALL PartInElemCheck(Position,iPart,GlobalElemID,ParticleInsideCheck)
CASE(TRIATRACKING)
  CALL ParticleInsideQuad3D(Position,GlobalElemID,ParticleInsideCheck)
CASE DEFAULT
  CALL abort(__STAMP__,'TrackingMethod not implemented! TrackingMethod =',IntInfoOpt=TrackingMethod)
END SELECT

END FUNCTION ParticleInsideCheck


SUBROUTINE ParticleSanityCheck(PartID)
!===================================================================================================================================
! this routine checks the LastPartPos and PartPosition for sanity
! 1) check if LastPartPos is within globalminmax of proc
! 2) check if ParticlePosition is within globalminmax
! 3) check if PartPosRef is within the element
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Mesh_Vars,              ONLY:offsetElem
USE MOD_Particle_Localization,  ONLY:PartInElemCheck
USE MOD_Particle_Mesh_Vars,     ONLY:GEO
USE MOD_Particle_Mesh_Vars,     ONLY:ElemBaryNGeo_Shared
USE MOD_Particle_Tracking_Vars, ONLY:TrackingMethod
USE MOD_Particle_Vars,          ONLY:PEM,PDM,LastPartPos,PartState
USE MOD_TimeDisc_Vars,          ONLY:iStage
#ifdef IMPA
USE MOD_Particle_Vars,          ONLY:PartIsImplicit,PartDtFrac
#endif /*IMPA*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: PartID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: ElemID
LOGICAL                          :: IsHit
REAL                             :: IntersectionPoint(1:3)
!===================================================================================================================================

IF(   (LastPartPos(1,PartID).GT.GEO%xmaxglob) &
  .OR.(LastPartPos(1,PartID).LT.GEO%xminglob) &
  .OR.(LastPartPos(2,PartID).GT.GEO%ymaxglob) &
  .OR.(LastPartPos(2,PartID).LT.GEO%yminglob) &
  .OR.(LastPartPos(3,PartID).GT.GEO%zmaxglob) &
  .OR.(LastPartPos(3,PartID).LT.GEO%zminglob) ) THEN
  IPWRITE(UNIt_stdOut,'(I0,A18,L1)')                            ' ParticleInside ', PDM%ParticleInside(PartID)
#ifdef IMPA
  IPWRITE(UNIt_stdOut,'(I0,A18,L1)')                            ' PartIsImplicit ', PartIsImplicit(PartID)
  IPWRITE(UNIt_stdOut,'(I0,A18,E27.16)')                       ' PartDtFrac ', PartDtFrac(PartID)
#endif /*IMPA*/
  IPWRITE(UNIt_stdOut,'(I0,A18,L1)')                            ' PDM%IsNewPart ', PDM%IsNewPart(PartID)
  IPWRITE(UNIt_stdOut,'(I0,A18,1X,A18,1X,A18)')                  '    min ', ' value ', ' max '
  IPWRITE(UNIt_stdOut,'(I0,A2,1X,E27.16,1X,E27.16,1X,E27.16)') ' x', GEO%xminglob, LastPartPos(1,PartID), GEO%xmaxglob
  IPWRITE(UNIt_stdOut,'(I0,A2,1X,E27.16,1X,E27.16,1X,E27.16)') ' y', GEO%yminglob, LastPartPos(2,PartID), GEO%ymaxglob
  IPWRITE(UNIt_stdOut,'(I0,A2,1X,E27.16,1X,E27.16,1X,E27.16)') ' z', GEO%zminglob, LastPartPos(3,PartID), GEO%zmaxglob
  CALL abort(&
         __STAMP__ &
         ,' LastPartPos outside of mesh. PartID=, iStage',PartID,REAL(iStage))
END IF
IF(   (PartState(1,PartID).GT.GEO%xmaxglob) &
  .OR.(PartState(1,PartID).LT.GEO%xminglob) &
  .OR.(PartState(2,PartID).GT.GEO%ymaxglob) &
  .OR.(PartState(2,PartID).LT.GEO%yminglob) &
  .OR.(PartState(3,PartID).GT.GEO%zmaxglob) &
  .OR.(PartState(3,PartID).LT.GEO%zminglob) ) THEN
  IPWRITE(UNIt_stdOut,'(I0,A18,L1)')                            ' ParticleInside ', PDM%ParticleInside(PartID)
#ifdef IMPA
      IPWRITE(UNIt_stdOut,'(I0,A18,L1)')                        ' PartIsImplicit ', PartIsImplicit(PartID)
      IPWRITE(UNIt_stdOut,'(I0,A18,E27.16)')                   ' PartDtFrac ', PartDtFrac(PartID)
#endif /*IMPA*/
  IPWRITE(UNIt_stdOut,'(I0,A18,3(1X,E27.16))')                  ' LastPartPos    ', LastPartPos(1:3,PartID)
  IPWRITE(UNIt_stdOut,'(I0,A18,3(1X,E27.16))')                  ' Velocity       ', PartState(4:6,PartID)
  IPWRITE(UNIt_stdOut,'(I0,A18,L1)')                            ' PDM%IsNewPart ', PDM%IsNewPart(PartID)
  IPWRITE(UNIt_stdOut,'(I0,A18,1X,A18,1X,A18)')                  '    min ', ' value ', ' max '
  IPWRITE(UNIt_stdOut,'(I0,A2,1X,E27.16,1X,E27.16,1X,E27.16)') ' x', GEO%xminglob, PartState(1,PartID), GEO%xmaxglob
  IPWRITE(UNIt_stdOut,'(I0,A2,1X,E27.16,1X,E27.16,1X,E27.16)') ' y', GEO%yminglob, PartState(2,PartID), GEO%ymaxglob
  IPWRITE(UNIt_stdOut,'(I0,A2,1X,E27.16,1X,E27.16,1X,E27.16)') ' z', GEO%zminglob, PartState(3,PartID), GEO%zmaxglob
  CALL abort(&
     __STAMP__ &
     ,' PartPos outside of mesh. PartID=, iStage',PartID,REAL(iStage))
END IF
IF(TrackingMethod.NE.REFMAPPING)THEN
  ElemID=PEM%GlobalElemID(PartID)
#ifdef CODE_ANALYZE
  CALL PartInElemCheck(PartState(1:3,PartID),PartID,ElemID,isHit,IntersectionPoint,CodeAnalyze_Opt=.TRUE.)
#else
  CALL PartInElemCheck(PartState(1:3,PartID),PartID,ElemID,isHit,IntersectionPoint)
#endif /*CODE_ANALYZE*/
  IF(.NOT.isHit)THEN  ! particle not inside
    IPWRITE(UNIT_stdOut,'(I0,A)') ' PartPos not inside of element! '
    IF(ElemID.GE.offSetElem+1.AND.ElemID.LE.offSetElem+PP_nElems)THEN
      IPWRITE(UNIT_stdOut,'(I0,A,I0)')       ' ElemID             ', ElemID
    END IF
    IPWRITE(UNIT_stdOut,'(I0,A,3(1X,E15.8))') ' ElemBaryNGeo:      ', ElemBaryNGeo_Shared(1:3,ElemID)
    IPWRITE(UNIT_stdOut,'(I0,A,3(1X,E15.8))') ' IntersectionPoint: ', IntersectionPoint
    IPWRITE(UNIT_stdOut,'(I0,A,3(1X,E15.8))') ' LastPartPos:       ', LastPartPos(1:3,PartID)
    IPWRITE(UNIT_stdOut,'(I0,A,3(1X,E15.8))') ' PartPos:           ', PartState(1:3,PartID)
    CALL abort(&
    __STAMP__ &
    ,'PartID=. ',PartID)
  END IF
END IF

END SUBROUTINE ParticleSanityCheck

END MODULE MOD_Particle_Tracking
