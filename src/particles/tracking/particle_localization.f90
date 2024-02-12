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

MODULE MOD_Particle_Localization
!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! required variables
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
INTERFACE SinglePointToElement
  MODULE PROCEDURE SinglePointToElement
END INTERFACE

INTERFACE LocateParticleInElement
  MODULE PROCEDURE LocateParticleInElement
END INTERFACE

INTERFACE PartInElemCheck
  MODULE PROCEDURE PartInElemCheck
END INTERFACE

INTERFACE CountPartsPerElem
  MODULE PROCEDURE CountPartsPerElem
END INTERFACE

PUBLIC:: SinglePointToElement
PUBLIC:: LocateParticleInElement
PUBLIC:: PartInElemCheck
PUBLIC:: CountPartsPerElem
PUBLIC:: PARTHASMOVED

CONTAINS

SUBROUTINE LocateParticleInElement(PartID,doHALO)
!----------------------------------------------------------------------------------------------------------------------------------!
! Finds a single particle in its host element
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
USE MOD_Particle_Vars          ,ONLY: PDM,PEM,PartState,PartPosRef
USE MOD_part_operations        ,ONLY: RemoveParticle
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: PartID
LOGICAL,INTENT(IN) :: doHalo
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: ElemID
!===================================================================================================================================
ElemID = SinglePointToElement(PartState(1:3,PartID),doHALO=doHALO)
PEM%GlobalElemID(PartID) = ElemID
IF(ElemID.EQ.-1)THEN
  CALL RemoveParticle(PartID)
ELSE
  PDM%ParticleInside(PartID)=.TRUE.
  PDM%isNewPart(PartID)     = .TRUE.
  IF(TrackingMethod.EQ.REFMAPPING) CALL GetPositionInRefElem(PartState(1:3,PartID),PartPosRef(1:3,PartID),ElemID)
END IF ! ElemID.EQ.-1
END SUBROUTINE LocateParticleInElement


!===================================================================================================================================
!> This function maps a 3D point to an element
!> returns elementID or -1 in case no element was found
!===================================================================================================================================
INTEGER FUNCTION SinglePointToElement(Pos3D,doHALO,doEmission_opt)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_Mesh_Vars              ,ONLY: offsetElem
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemRadius2NGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO,ElemEpsOneCell
USE MOD_Particle_Mesh_Vars     ,ONLY: FIBGM_nElems, FIBGM_offsetElem, FIBGM_Element
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID,GetGlobalElemID
USE MOD_Particle_Mesh_Tools    ,ONLY: ParticleInsideQuad3D
USE MOD_Particle_Tracking_Vars ,ONLY: Distance,ListDistance,TrackingMethod
USE MOD_Utils                  ,ONLY: InsertionSort
#if USE_MPI
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemBaryNGeo_Shared
#else
USE MOD_Mesh_Vars              ,ONLY: ElemBaryNGeo
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,INTENT(IN)    :: Pos3D(1:3)
LOGICAL,INTENT(IN) :: doHalo
LOGICAL,INTENT(IN),OPTIONAL :: doEmission_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iBGMElem,nBGMElems,ElemID,CNElemID,iBGM,jBGM,kBGM
REAL    :: Distance2, RefPos(1:3)
REAL    :: Det(6,2)
LOGICAL :: InElementCheck, doEmission
!===================================================================================================================================
#if USE_MPI
ASSOCIATE(ElemBaryNGeo => ElemBaryNGeo_Shared)
#endif

SinglePointToElement = -1

! --- get background mesh cell of point
iBGM = CEILING((Pos3D(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))
iBGM = MAX(MIN(GEO%FIBGMimax,iBGM),GEO%FIBGMimin)
jBGM = CEILING((Pos3D(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
jBGM = MAX(MIN(GEO%FIBGMjmax,jBGM),GEO%FIBGMjmin)
kBGM = CEILING((Pos3D(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
kBGM = MAX(MIN(GEO%FIBGMkmax,kBGM),GEO%FIBGMkmin)

!--- check all cells associated with this background mesh cell
nBGMElems = FIBGM_nElems(iBGM,jBGM,kBGM)

! get closest element barycenter
Distance = -1.

ListDistance=0
DO iBGMElem = 1,nBGMElems
  ElemID   = FIBGM_Element(FIBGM_offsetElem(iBGM,jBGM,kBGM)+iBGMElem)
  CNElemID = GetCNElemID(ElemID)

  Distance2 = SUM((Pos3D(1:3)-ElemBaryNGeo(1:3,CNElemID))**2.)

  ! element in range
  Distance(iBGMElem)     = MERGE(Distance2,-1.,Distance2.LE.ElemRadius2NGeo(CNElemID))
  ListDistance(iBGMElem) = ElemID
END DO ! nBGMElems

IF(ALMOSTEQUAL(MAXVAL(Distance),-1.))THEN
  RETURN
END IF

IF (nBGMElems.GT.1) CALL InsertionSort(Distance(1:nBGMElems),ListDistance(1:nBGMElems),nBGMElems)

! loop through sorted list and start by closest element
InElementCheck = .FALSE.
DO iBGMElem = 1,nBGMElems
  IF (ALMOSTEQUAL(Distance(iBGMElem),-1.)) CYCLE

  ElemID = ListDistance(iBGMElem)

  IF (.NOT.DoHALO) THEN
    IF (ElemID.LT.offsetElem+1 .OR. ElemID.GT.offsetElem+PP_nElems) CYCLE
  END IF

  SELECT CASE(TrackingMethod)
    CASE(TRIATRACKING)
      CALL ParticleInsideQuad3D(Pos3D(1:3),ElemID,InElementCheck,Det)

    CASE(TRACING)
      CALL GetPositionInRefElem(Pos3D(1:3),RefPos,ElemID)
      IF (MAXVAL(ABS(RefPos)).LE.1.0) InElementCheck = .TRUE.

    CASE(REFMAPPING)
      CNElemID = GetCNElemID(ElemID)
      CALL GetPositionInRefElem(Pos3D(1:3),RefPos,ElemID)

      IF (PRESENT(doEmission_opt)) THEN; doEmission = doEmission_opt
      ELSE;                              doEmission = .FALSE.
      END IF

      IF(doEmission) THEN
        IF (MAXVAL(ABS(RefPos)).LE.1.0) InElementCheck = .TRUE.
      ELSE
        IF (MAXVAL(ABS(RefPos)).LE.ElemEpsOneCell(CNElemID)) InElementCheck = .TRUE.
      END IF
  END SELECT

  IF (InElementCheck) THEN
    SinglePointToElement = ElemID
    RETURN
  END IF
END DO ! iBGMElem

#if USE_MPI
END ASSOCIATE
#endif /*USE_MPI*/

END FUNCTION SinglePointToElement


SUBROUTINE PartInElemCheck(PartPos_In,PartID,ElemID,FoundInElem,IntersectPoint_Opt, &
#ifdef CODE_ANALYZE
                           Sanity_Opt,Tol_Opt,CodeAnalyze_Opt)
#else
                           Tol_Opt)
#endif /*CODE_ANALYZE*/
!===================================================================================================================================
! Checks if particle is in Element
!===================================================================================================================================
! MODULES
USE MOD_Globals                ,ONLY: VECNORM
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID,GetCNSideID
USE MOD_Particle_Intersection  ,ONLY: ComputePlanarRectIntersection
USE MOD_Particle_Intersection  ,ONLY: ComputePlanarCurvedIntersection
USE MOD_Particle_Intersection  ,ONLY: ComputeBiLinearIntersection
USE MOD_Particle_Intersection  ,ONLY: ComputeCurvedIntersection
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Mesh_Vars     ,ONLY: SideInfo_Shared
USE MOD_Particle_Surfaces      ,ONLY: CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Surfaces_Vars ,ONLY: SideType,SideNormVec
USE MOD_Particle_Vars          ,ONLY: LastPartPos
#ifdef CODE_ANALYZE
USE MOD_Globals                ,ONLY: MyRank,UNIT_stdout
USE MOD_Particle_Tracking_Vars ,ONLY: PartOut,MPIRankOut
USE MOD_Particle_Surfaces      ,ONLY: OutputBezierControlPoints
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
USE MOD_Particle_Intersection  ,ONLY: OutputTrajectory
#endif /*CODE_ANALYZE*/
#if USE_MPI
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemBaryNGeo_Shared
#else
USE MOD_Mesh_Vars              ,ONLY: ElemBaryNGeo
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                       :: ElemID,PartID
REAL,INTENT(IN)                          :: PartPos_In(1:3)
#ifdef CODE_ANALYZE
LOGICAL,INTENT(IN),OPTIONAL              :: CodeAnalyze_Opt
LOGICAL,INTENT(IN),OPTIONAL              :: Sanity_Opt
#endif /*CODE_ANALYZE*/
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)                      :: FoundInElem
REAL,INTENT(OUT),OPTIONAL                :: IntersectPoint_Opt(1:3)
REAL,INTENT(OUT),OPTIONAL                :: Tol_Opt
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                  :: CNElemID
INTEGER                                  :: ilocSide,flip,SideID,CNSideID
REAL                                     :: PartTrajectory(1:3),NormVec(1:3)
REAL                                     :: lengthPartTrajectory,PartPos(1:3),LastPosTmp(1:3)
LOGICAL                                  :: isHit
REAL                                     :: alpha,eta,xi,IntersectPoint(1:3)
!===================================================================================================================================

#if USE_MPI
ASSOCIATE(ElemBaryNGeo => ElemBaryNGeo_Shared)
#endif

IF(PRESENT(tol_Opt)) tol_Opt=-1.

CNElemID = GetCNElemID(ElemID)

! virtual move to element barycenter
LastPosTmp(1:3)         = LastPartPos(1:3,PartID)
LastPartPos(1:3,PartID) = ElemBaryNGeo(1:3,CNElemID)
PartPos(1:3)            = PartPos_In(1:3)

! get trajectory from element barycenter to current position
PartTrajectory       = PartPos - LastPartPos(1:3,PartID)
lengthPartTrajectory = VECNORM(PartTrajectory(1:3))

! output the part trajectory
#ifdef CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      IPWRITE(UNIT_stdout,*) ' --------------------------------------------- '
      IPWRITE(UNIT_stdout,*) ' PartInElemCheck '
      CALL OutputTrajectory(PartID,PartPos,PartTrajectory,lengthPartTrajectory)
    END IF
  END IF
#endif /*CODE_ANALYZE*/

! we found the particle on the element barycenter
IF (ALMOSTZERO(lengthPartTrajectory)) THEN
  FoundInElem = .TRUE.
  LastPartPos(1:3,PartID) = LastPosTmp(1:3)
  RETURN
END IF

! normalize the part trajectory vector
PartTrajectory = PartTrajectory/lengthPartTrajectory

! reset intersection counter and alpha
isHit = .FALSE.
alpha = -1.

DO ilocSide = 1,6
  SideID   = GetGlobalNonUniqueSideID(ElemID,iLocSide)
  CNSideID = GetCNSideID(SideID)
  flip     = MERGE(0,MOD(SideInfo_Shared(SIDE_FLIP,SideID),10),SideInfo_Shared(SIDE_ID,SideID).GT.0)

  SELECT CASE(SideType(CNSideID))
    CASE(PLANAR_RECT)
      CALL ComputePlanarRectIntersection(  ishit,PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,flip,SideID)
    CASE(PLANAR_CURVED)
      CALL ComputePlanarCurvedIntersection( isHit,PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,flip,SideID)
    CASE(BILINEAR,PLANAR_NONRECT)
        CALL ComputeBiLinearIntersection(     isHit,PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,     SideID,ElemCheck_Opt=.TRUE.)
    CASE(CURVED)
      CALL ComputeCurvedIntersection(       isHit,PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,flip,SideID,ElemCheck_Opt=.TRUE.)
  END SELECT

#ifdef CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,'(15("="))')
      WRITE(UNIT_stdout,'(A)')           '     | Output after compute intersection (PartInElemCheck): '
      WRITE(UNIT_stdout,'(2(A,I0),A,L1)') '     | SideType: ',SideType(CNSideID)  ,' | SideID: ',SideID,'| Hit: ',isHit
      WRITE(UNIT_stdout,'(2(A,G0))')     '     | LengthPT: ',LengthPartTrajectory,' | Alpha: ',Alpha
      WRITE(UNIT_stdout,'(A,2(1X,G0))')   '     | Intersection xi/eta: ',xi,eta
    END IF
  END IF
  IF(PRESENT(Sanity_Opt))THEN
    IF(Sanity_Opt)THEN
      IF(alpha.GT.-1)THEN
        ! alpha is going from barycenter to point
        ! here, the tolerance for the ratio alpha/LengthPartTrajectory for tracing with element-corners is determined.
        IF(PRESENT(tol_Opt)) tol_Opt=MAX(ABS(1.-alpha/LengthPartTrajectory),tol_Opt)
        ! mark element as trouble element if rel. tol from alpha/LengthPartTrajectory to 1 > 1e-4
        ! tolerance 1e-4 is from ANSA_BOX grid (experimental, arbitrary)
        IF(ALMOSTEQUALRELATIVE(alpha/LengthPartTrajectory,1.,0.002)) THEN
          alpha=-1
        ELSE
          print*,'alpha',alpha,LengthPartTrajectory,alpha/LengthPartTrajectory,ABS(1.-alpha/LengthPartTrajectory),tol_Opt
        END IF
      END IF
    END IF
  END IF
  ! Dirty fix for PartInElemCheck if Lastpartpos is almost on side (tolerance issues)
  IF(PRESENT(CodeAnalyze_Opt))THEN
    IF(CodeAnalyze_Opt)THEN
      IF((alpha)/LengthPartTrajectory.GT.0.9)THEN
        alpha = -1.0
      END IF
    END IF
  END IF
#endif /*CODE_ANALYZE*/
  IF(alpha.GT.-1)THEN
    SELECT CASE(SideType(CNSideID))
      CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
        CNSideID = GetCNSideID(SideID)
        NormVec  = SideNormVec(1:3,CNSideID)
      CASE(BILINEAR)
        CALL CalcNormAndTangBilinear(nVec=NormVec,xi=xi,eta=eta,SideID=SideID)
      CASE(CURVED)
        CALL CalcNormAndTangBezier(nVec=NormVec,xi=xi,eta=eta,SideID=SideID)
    END SELECT
    IF(flip.NE.0) NormVec = -NormVec
    IntersectPoint = LastPartPos(1:3,PartID)+alpha*PartTrajectory

#ifdef CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,*) '     | alpha          ',alpha
      WRITE(UNIT_stdout,*) '     | Normal vector  ',NormVec
      WRITE(UNIT_stdout,*) '     | PartTrajectory ',PartTrajectory
      WRITE(UNIT_stdout,*) '     | Dotprod        ',DOT_PRODUCT(NormVec,PartTrajectory)
      WRITE(UNIT_stdout,*) '     | Point 2        ', LastPartPos(1:3,PartID)+alpha*PartTrajectory+NormVec
      WRITE(UNIT_stdout,*) '     | Beziercontrolpoints3d-x'
      CALL OutputBezierControlPoints(BezierControlPoints3D_in=BezierControlPoints3D(1:3,:,:,SideID))
    END IF
  END IF
#endif /*CODE_ANALYZE*/

    IF(DOT_PRODUCT(NormVec,PartTrajectory).LT.0.)THEN
      alpha=-1.0
    ELSE
      EXIT
    END IF
  END IF
END DO ! ilocSide

FoundInElem=.TRUE.

IF(PRESENT(IntersectPoint_Opt)) IntersectPoint_Opt=0.
IF(alpha.GT.-1) THEN
  FoundInElem=.FALSE.
  IF(PRESENT(IntersectPoint_Opt)) IntersectPoint_Opt=IntersectPoint
END IF

! reset the LastParPos to its original value
LastPartPos(1:3,PartID) = LastPosTmp(1:3)

#if USE_MPI
END ASSOCIATE
#endif /*USE_MPI*/

END SUBROUTINE PartInElemCheck


SUBROUTINE CountPartsPerElem(ResetNumberOfParticles)
!===================================================================================================================================
! count number of particles in element
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_LoadBalance_Vars ,ONLY: nPartsPerElem
USE MOD_Particle_Vars    ,ONLY: PDM,PEM
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: ResetNumberOfParticles
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iPart, ElemID
!===================================================================================================================================
! DO NOT NULL this here, if e.g. this routine is called in between RK-stages in which particles are created
IF(ResetNumberOfParticles)THEN
  nPartsPerElem=0
END IF
! loop over all particles and add them up
DO iPart=1,PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart))THEN
    ElemID = PEM%LocalElemID(iPart)
    IF(ElemID.LE.PP_nElems)THEN
      nPartsPerElem(ElemID)=nPartsPerElem(ElemID)+1
    END IF
  END IF
END DO ! iPart=1,PDM%ParticleVecLength

END SUBROUTINE CountPartsPerElem


PPURE FUNCTION PARTHASMOVED(lengthPartTrajectory,ElemRadiusNGeo)
!================================================================================================================================
! check if particle has moved significantly within an element
!================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                      :: lengthPartTrajectory
REAL,INTENT(IN)                      :: ElemRadiusNGeo
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL                              :: PARTHASMOVED
!================================================================================================================================

IF(ALMOSTZERO(lengthPartTrajectory/ElemRadiusNGeo))THEN
  PARTHASMOVED=.FALSE.
ELSE
  PARTHASMOVED=.TRUE.
END IF

END FUNCTION PARTHASMOVED

END MODULE MOD_Particle_Localization