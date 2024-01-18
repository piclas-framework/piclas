!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (gitlab.com/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

MODULE MOD_Photon_Tracking
!===================================================================================================================================
! Routines for photon tracking in radiave transfer solver
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC

INTERFACE PhotonTriaTracking
  MODULE PROCEDURE PhotonTriaTracking
END INTERFACE

PUBLIC :: PhotonTriaTracking, Photon2DSymTracking
PUBLIC :: InitPhotonSurfSample
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> Allocate photon surface sampling containers
!===================================================================================================================================
SUBROUTINE InitPhotonSurfSample()
! MODULES
USE MOD_Globals
USE MOD_Photon_TrackingVars    ,ONLY: PhotonSurfSideArea,PhotonSurfSideSamplingMidPoints
USE MOD_Particle_Boundary_Vars ,ONLY: nComputeNodeSurfTotalSides
USE MOD_Particle_Boundary_Vars ,ONLY: SurfSide2GlobalSide
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared_Vars        ,ONLY: myComputeNodeRank,nComputeNodeProcessors
USE MOD_Photon_TrackingVars    ,ONLY: PhotonSurfSideSamplingMidPoints_Shared,PhotonSurfSideSamplingMidPoints_Shared_Win
USE MOD_Photon_TrackingVars    ,ONLY: PhotonSurfSideArea_Shared,PhotonSurfSideArea_Shared_Win
#else
USE MOD_Particle_Boundary_Vars ,ONLY: nGlobalSurfSides
#endif /*USE_MPI*/
USE MOD_Particle_Vars          ,ONLY: Symmetry
USE MOD_Basis                  ,ONLY: LegendreGaussNodesAndWeights
USE MOD_DSMC_Symmetry          ,ONLY: DSMC_2D_CalcSymmetryArea, DSMC_1D_CalcSymmetryArea
USE MOD_Mesh_Vars              ,ONLY: NGeo
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars     ,ONLY: SideInfo_Shared,NodeCoords_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemSideNodeID_Shared
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
USE MOD_Particle_Surfaces      ,ONLY: EvaluateBezierPolynomialAndGradient
USE MOD_RayTracing_Vars        ,ONLY: Ray
USE MOD_Interpolation          ,ONLY: GetNodesAndWeights
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
! Space-separated list of input and output types. Use: (int|real|logical|...)_(in|out|inout)_dim(n)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                :: iSide,firstSide,lastSide
! surface area
INTEGER                                :: SideID,ElemID,CNElemID,LocSideID
INTEGER                                :: p,q,iSample,jSample
INTEGER                                :: TriNum, Node1, Node2
REAL                                   :: area,nVal
REAL,DIMENSION(2,3)                    :: gradXiEta3D
REAL,DIMENSION(:),ALLOCATABLE          :: Xi_NGeo,wGP_NGeo
REAL                                   :: XiOut(1:2),E,F,G,D,tmp1,tmpI2,tmpJ2
REAL                                   :: xNod(3), Vector1(3), Vector2(3), nx, ny, nz
LOGICAL                                :: UseBezierControlPointsForArea
REAL,ALLOCATABLE                       :: xIP_VISU(:),wIP_VISU(:)
REAL,ALLOCATABLE                       :: RayXiEQ_SurfSample(:)            ! position of RayXiEQ_SurfSample
REAL                                   :: dRayXiEQ_SurfSample              ! deltaXi in [-1,1]
!===================================================================================================================================

#if USE_MPI
CALL Allocate_Shared((/Ray%nSurfSample,Ray%nSurfSample,nComputeNodeSurfTotalSides/),PhotonSurfSideArea_Shared_Win,PhotonSurfSideArea_Shared)
CALL MPI_WIN_LOCK_ALL(0,PhotonSurfSideArea_Shared_Win,IERROR)
CALL Allocate_Shared((/3,Ray%nSurfSample,Ray%nSurfSample,nComputeNodeSurfTotalSides/),PhotonSurfSideSamplingMidPoints_Shared_Win,PhotonSurfSideSamplingMidPoints_Shared)
CALL MPI_WIN_LOCK_ALL(0,PhotonSurfSideSamplingMidPoints_Shared_Win,IERROR)
PhotonSurfSideArea => PhotonSurfSideArea_Shared
PhotonSurfSideSamplingMidPoints => PhotonSurfSideSamplingMidPoints_Shared

firstSide = INT(REAL( myComputeNodeRank   )*REAL(nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1))*REAL(nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))
#else
ALLOCATE(PhotonSurfSideArea(1:Ray%nSurfSample,1:Ray%nSurfSample,1:nComputeNodeSurfTotalSides))
ALLOCATE(PhotonSurfSideSamplingMidPoints(1:3,1:Ray%nSurfSample,1:Ray%nSurfSample,1:nComputeNodeSurfTotalSides))

firstSide = 1
lastSide  = nGlobalSurfSides
#endif /*USE_MPI*/

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
  PhotonSurfSideArea=0.
  PhotonSurfSideSamplingMidPoints=0.
#if USE_MPI
END IF
CALL BARRIER_AND_SYNC(PhotonSurfSideArea_Shared_Win,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(PhotonSurfSideSamplingMidPoints_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/

! Calculate equidistant surface points
ALLOCATE(RayXiEQ_SurfSample(0:Ray%nSurfSample))
dRayXiEQ_SurfSample =2./REAL(Ray%nSurfSample)
DO q=0,Ray%nSurfSample
  RayXiEQ_SurfSample(q) = dRayXiEQ_SurfSample * REAL(q) - 1.
END DO

! get interpolation points and weights
ALLOCATE( Xi_NGeo( 0:NGeo)  &
        , wGP_NGeo(0:NGeo) )
CALL LegendreGaussNodesAndWeights(NGeo,Xi_NGeo,wGP_NGeo)

! compute area of sub-faces
tmp1=dRayXiEQ_SurfSample/2.0 !(b-a)/2

ALLOCATE(xIP_VISU(0:Ray%nSurfSample),wIP_VISU(0:Ray%nSurfSample))
! Build basis for surface sampling on VISU nodes, and not on Ray%NodeType, which default to VISU_INNER
! because VISU nodes are hard-coded in piclas2vtk
CALL GetNodesAndWeights(Ray%nSurfSample, 'VISU', xIP_VISU, wIP=wIP_VISU)

DO iSide = firstSide,LastSide
  ! get global SideID. This contains only nonUniqueSide, no special mortar treatment required
  SideID = SurfSide2GlobalSide(SURF_SIDEID,iSide)

  UseBezierControlPointsForArea = .FALSE.

  IF (TrackingMethod.EQ.TRIATRACKING) THEN
    ElemID    = SideInfo_Shared(SIDE_ELEMID ,SideID)
    CNElemID  = GetCNElemID(ElemID)
    LocSideID = SideInfo_Shared(SIDE_LOCALID,SideID)
    IF((Symmetry%Order.NE.3).AND.Ray%nSurfSample.GT.1) CALL abort(__STAMP__,'Ray%nSurfSample>1 not implemented for this symmetry!')

    IF(Symmetry%Order.EQ.3) THEN
      ! Check if triangles are used for the calculation of the surface area or not
      IF(Ray%nSurfSample.GT.1)THEN
        ! Do not use triangles
        UseBezierControlPointsForArea = .TRUE.
      ELSE
        xNod(1:3) = NodeCoords_Shared(1:3,ElemSideNodeID_Shared(1,LocSideID,CNElemID)+1)
        area = 0.
        DO TriNum = 1,2
          Node1 = TriNum+1     ! normal = cross product of 1-2 and 1-3 for first triangle
          Node2 = TriNum+2     !          and 1-3 and 1-4 for second triangle
          Vector1(1:3) = NodeCoords_Shared(1:3,ElemSideNodeID_Shared(Node1,LocSideID,CNElemID)+1) - xNod(1:3)
          Vector2(1:3) = NodeCoords_Shared(1:3,ElemSideNodeID_Shared(Node2,LocSideID,CNElemID)+1) - xNod(1:3)
          nx = - Vector1(2) * Vector2(3) + Vector1(3) * Vector2(2) !NV (inwards)
          ny = - Vector1(3) * Vector2(1) + Vector1(1) * Vector2(3)
          nz = - Vector1(1) * Vector2(2) + Vector1(2) * Vector2(1)
          nVal = SQRT(nx*nx + ny*ny + nz*nz)
          area = area + nVal/2.
        END DO
        PhotonSurfSideArea(1,1,iSide) = area
      END IF ! Ray%nSurfSample.GT.1
    ELSE IF(Symmetry%Order.EQ.2) THEN
      PhotonSurfSideArea(1,1,iSide) = DSMC_2D_CalcSymmetryArea(LocSideID, CNElemID)
    ELSE IF(Symmetry%Order.EQ.1) THEN
      PhotonSurfSideArea(1,1,iSide) = DSMC_1D_CalcSymmetryArea(LocSideID, CNElemID)
    END IF
  ELSE ! TrackingMethod.NE.TRIATRACKING
    UseBezierControlPointsForArea = .TRUE.
  END IF ! TrackingMethod.EQ.TRIATRACKIN

  ! Instead of triangles use Bezier control points (curved or triangle tracking with Ray%nSurfSample>1)
  IF(UseBezierControlPointsForArea)THEN
    DO jSample=1,Ray%nSurfSample
      DO iSample=1,Ray%nSurfSample
        area=0.
        tmpI2=(RayXiEQ_SurfSample(iSample-1)+RayXiEQ_SurfSample(iSample))/2. ! (a+b)/2
        tmpJ2=(RayXiEQ_SurfSample(jSample-1)+RayXiEQ_SurfSample(jSample))/2. ! (a+b)/2
        ASSOCIATE( xi => 0.5*(xIP_VISU(iSample)+xIP_VISU(iSample-1)), eta => 0.5*(xIP_VISU(jSample)+xIP_VISU(jSample-1)) )
          CALL EvaluateBezierPolynomialAndGradient((/xi,eta/),NGeo,3,BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID) &
              ,Point=PhotonSurfSideSamplingMidPoints(1:3,iSample,jSample,iSide))
        END ASSOCIATE
        DO q=0,NGeo
          DO p=0,NGeo
            XiOut(1)=tmp1*Xi_NGeo(p)+tmpI2
            XiOut(2)=tmp1*Xi_NGeo(q)+tmpJ2
            CALL EvaluateBezierPolynomialAndGradient(XiOut,NGeo,3,BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID) &
                                                    ,Gradient=gradXiEta3D)
            ! calculate first fundamental form
            E=DOT_PRODUCT(gradXiEta3D(1,1:3),gradXiEta3D(1,1:3))
            F=DOT_PRODUCT(gradXiEta3D(1,1:3),gradXiEta3D(2,1:3))
            G=DOT_PRODUCT(gradXiEta3D(2,1:3),gradXiEta3D(2,1:3))
            D=SQRT(E*G-F*F)
            area = area+tmp1*tmp1*D*wGP_NGeo(p)*wGP_NGeo(q)
          END DO
        END DO
        PhotonSurfSideArea(iSample,jSample,iSide) = area
      END DO ! iSample=1,Ray%nSurfSample
    END DO ! jSample=1,Ray%nSurfSample
  END IF ! UseBezierControlPointsForArea

END DO ! iSide = firstSide,lastSide

#if USE_MPI
CALL BARRIER_AND_SYNC(PhotonSurfSideArea_Shared_Win,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(PhotonSurfSideSamplingMidPoints_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/

END SUBROUTINE InitPhotonSurfSample


!===================================================================================================================================
!> Deallocate photon surface sampling containers
!===================================================================================================================================
SUBROUTINE FinalizePhotonSurfSample()
! MODULES
USE MOD_Globals
USE MOD_Photon_TrackingVars ,ONLY: PhotonSurfSideArea,PhotonSurfSideSamplingMidPoints
#if USE_MPI
USE MOD_MPI_Shared_Vars     ,ONLY: MPI_COMM_SHARED
USE MOD_Photon_TrackingVars ,ONLY: PhotonSurfSideSamplingMidPoints_Shared,PhotonSurfSideSamplingMidPoints_Shared_Win
USE MOD_Photon_TrackingVars ,ONLY: PhotonSurfSideArea_Shared,PhotonSurfSideArea_Shared_Win,PhotonSampWall_Shared_Win_allocated
USE MOD_MPI_Shared
#endif /*USE_MPI*/
USE MOD_Photon_TrackingVars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
#if USE_MPI
! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
CALL UNLOCK_AND_FREE(PhotonSurfSideSamplingMidPoints_Shared_Win)
CALL UNLOCK_AND_FREE(PhotonSurfSideArea_Shared_Win)
IF(PhotonSampWall_Shared_Win_allocated) CALL UNLOCK_AND_FREE(PhotonSampWall_Shared_Win)
PhotonSampWall_Shared_Win_allocated = .FALSE.
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
ADEALLOCATE(PhotonSampWall_Shared)
ADEALLOCATE(PhotonSurfSideSamplingMidPoints_Shared)
ADEALLOCATE(PhotonSurfSideArea_Shared)
#endif /*USE_MPI*/

ADEALLOCATE(PhotonSurfSideSamplingMidPoints)
ADEALLOCATE(PhotonSurfSideArea)
END SUBROUTINE FinalizePhotonSurfSample


SUBROUTINE PhotonTriaTracking()
!===================================================================================================================================
! Routine for tracking of moving particles and boundary interaction using triangulated sides.
! 1) Loop over all particles that are still inside
!    2) Perform tracking until the particle is considered "done" (either localized or deleted)
!       2a) Perform a check based on the determinant of (3x3) matrix of the vectors from the particle position to the nodes of each
!           triangle (ParticleInsideQuad3D)
!       2b) If particle is not within the given element in a), the side through which the particle went is determined by checking
!           each side of the element (ParticleThroughSideCheck3DFast)
!       2c) If no sides are found, the particle is deleted (very rare case). If multiple possible sides are found, additional
!           treatment is required, where the particle path is reconstructed (which side was crossed first) by comparing the ratio
!           the determinants
!    3) In case of a boundary, determine the intersection and perform the appropriate boundary interaction (GetBoundaryInteraction)
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound
USE MOD_Photon_TrackingVars     ,ONLY: PhotonProps,PhotonModeBPO,UsePhotonTriaTracking
USE MOD_RadiationTrans_Vars     ,ONLY: RadObservation_Emission, RadObservationPointMethod, RadObservation_EmissionPart
USE MOD_Photon_TrackingTools    ,ONLY: PhotonThroughSideCheck3DFast, PhotonIntersectionWithSide,CalcAbsoprtion,PhotonOnLineOfSight
USE MOD_Photon_TrackingTools    ,ONLY: PerfectPhotonReflection, DiffusePhotonReflection, CalcWallAbsoprtion, PointInObsCone
USE MOD_Photon_TrackingTools    ,ONLY: PeriodicPhotonBC
USE MOD_Photon_TrackingTools    ,ONLY: PhotonIntersectSensor
USE MOD_Particle_Boundary_Tools ,ONLY: StoreBoundaryParticleProperties
USE MOD_part_tools              ,ONLY: StoreLostPhotonProperties
USE MOD_RadiationTrans_Vars     ,ONLY: RadiationAbsorptionModel
USE MOD_RayTracing_Vars         ,ONLY: RayForceAbsorption
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalNonUniqueSideID
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(KIND=8),PARAMETER        :: MaxIterPhoton(1:2)=(/1000,1000000000/) ! Maximum number of cycles in the do while loop for each photon for bilinear and TriaTracking
INTEGER(KIND=8)                  :: IterPhoton(2)
INTEGER                          :: NblocSideID, NbElemID, ind, nbSideID, nMortarElems, BCType, localSideID, iPBC, GlobSideID
INTEGER                          :: ElemID,OldElemID,nlocSides
INTEGER                          :: LocalSide
INTEGER                          :: NrOfThroughSides, ind2
INTEGER                          :: SideID,TempSideID,iLocSide
INTEGER                          :: TriNum, LocSidesTemp(1:6),TriNumTemp(1:6), GlobSideTemp(1:6)
INTEGER                          :: SecondNrOfThroughSides, indSide 
INTEGER                          :: DoneLastElem(1:4,1:6) ! 1:3: 1=Element,2=LocalSide,3=TriNum 1:2: 1=last 2=beforelast
LOGICAL                          :: ThroughSide, Done, TempPointSelect(2), LastInterPointSelect(2),DummyPointSelect(2),FallBack
LOGICAL                          :: oldElemIsMortar, isMortarSideTemp(1:6), doCheckSide, InterPointSelect(2), InterPointSelectTemp(2,6)
REAL                             :: minRatio, intersecDist, intersecDistVec(3)
REAL                             :: IntersectionPos(1:3), IntersectionPosTemp(1:3)
REAL                             :: DistTemp(1:6)
LOGICAL                          :: PhotonLost
!===================================================================================================================================
Done                 = .FALSE.
ElemID               = PhotonProps%ElemID
SideID               = 0
GlobSideID           = 0
DoneLastElem(:,:)    = 0
InterPointSelect     = .FALSE.
LastInterPointSelect = .FALSE.
IterPhoton           = 0_8
TriNum               = 1

! 1) Loop tracking until Photon is considered "done" (either absorbed or deleted)
THEWHILELOOP: DO WHILE (.NOT.Done)
  IterPhoton(2) = IterPhoton(2) + 1_8 ! Stop when MaxIterPhoton(2) is reached with this counter
  IF(IterPhoton(2).GE.MaxIterPhoton(2))THEN
    CALL StoreLostPhotonProperties(ElemID,TRIM(__FILE__),__LINE__,99999)
    Done = .TRUE.
    EXIT THEWHILELOOP
  END IF
  PhotonLost=.FALSE.
  InterPointSelectTemp = .FALSE.
  oldElemIsMortar = .FALSE.
  NrOfThroughSides = 0
  LocSidesTemp(:) = 0
  DistTemp = 0.0  
  TriNumTemp(:) = 0
  GlobSideTemp = 0
  isMortarSideTemp = .FALSE.
  FallBack = .FALSE.
  nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)
  LocSideLoop: DO iLocSide=1,nlocSides
    TempSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide
    localSideID = SideInfo_Shared(SIDE_LOCALID,TempSideID)
    ! Side is not one of the 6 local sides
    IF (localSideID.LE.0) CYCLE
    NbElemID = SideInfo_Shared(SIDE_NBELEMID,TempSideID)
    IF (NbElemID.LT.0) THEN ! Mortar side
      nMortarElems = MERGE(4,2,SideInfo_Shared(SIDE_NBELEMID,TempSideID).EQ.-1)
      DO ind = 1, nMortarElems
        nbSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide + ind
        NbElemID = SideInfo_Shared(SIDE_NBELEMID,nbSideID)
        ! If small mortar element not defined, abort. Every available information on the compute-node is kept in shared memory, so
        ! no way to recover it during runtime
        IF (NbElemID.LT.1) CALL ABORT(__STAMP__,'Small mortar element not defined!',ElemID)
        ! For small mortar sides, SIDE_NBSIDEID contains the SideID of the corresponding big mortar side
        nbSideID = SideInfo_Shared(SIDE_NBSIDEID,nbSideID)
        NblocSideID =  SideInfo_Shared(SIDE_LOCALID,nbSideID)
        DO TriNum = 1,2
          ThroughSide = .FALSE.
          CALL PhotonThroughSideCheck3DFast(NblocSideID,NbElemID,ThroughSide,TriNum, .TRUE.)
          IF (ThroughSide) THEN
            ! Store the information for this side for future checks, if this side was already treated
            oldElemIsMortar = .TRUE.
            NrOfThroughSides = NrOfThroughSides + 1
            LocSidesTemp(NrOfThroughSides) = NblocSideID
            TriNumTemp(NrOfThroughSides) = TriNum
            GlobSideTemp(NrOfThroughSides) = nbSideID
            isMortarSideTemp(NrOfThroughSides)  = .TRUE.
            SideID = nbSideID
            LocalSide = NblocSideID
          END IF
        END DO
      END DO
    ELSE  ! Regular side
      ! Select A) TriaTracking or 2) Tracing on bilinear sides
      IF(UsePhotonTriaTracking)THEN
        ! A) TriaTracking
        DO TriNum = 1,2
          ThroughSide = .FALSE.
          CALL PhotonThroughSideCheck3DFast(localSideID,ElemID,ThroughSide,TriNum)
          IF (ThroughSide) THEN
            NrOfThroughSides = NrOfThroughSides + 1
            LocSidesTemp(NrOfThroughSides) = localSideID
            TriNumTemp(NrOfThroughSides) = TriNum
            GlobSideTemp(NrOfThroughSides) = TempSideID
            SideID = TempSideID
            LocalSide = localSideID
          END IF
        END DO
      ELSE
        ! 2) Tracing on bilinear sides (bilinear algorithm for intersection calculation))
        GlobSideID   = GetGlobalNonUniqueSideID(ElemID,localSideID)
        ThroughSide = .FALSE.
        !CALL ComputeBiLinearIntersection(foundHit,PartTrajectory,lengthPartTrajectory,locAlpha,xi,eta,iPart,SideID,alpha2=currentIntersect%alpha)
        IF (SideInfo_Shared(SIDE_NBELEMID,TempSideID).EQ.DoneLastElem(1,1)) THEN
          DummyPointSelect = LastInterPointSelect 
        ELSE
          DummyPointSelect = (/.FALSE.,.FALSE./)
        END IF
        CALL PhotonComputeBiLinearIntersection(ThroughSide,GlobSideID, intersecDist,TempPointSelect, DummyPointSelect)
!        print*,intersecDist, SideInfo_Shared(SIDE_NBELEMID,TempSideID), DoneLastElem(1,1), DummyPointSelect, 'ThroughSide',ThroughSide,'which', TempPointSelect
        IF(ThroughSide)THEN
          NrOfThroughSides = NrOfThroughSides + 1
          LocSidesTemp(NrOfThroughSides) = localSideID
          GlobSideTemp(NrOfThroughSides) = TempSideID
          DistTemp(NrOfThroughSides) = intersecDist
          SideID = TempSideID
          LocalSide = localSideID
          InterPointSelect = TempPointSelect(1:2)
          InterPointSelectTemp(1:2,NrOfThroughSides) = TempPointSelect(1:2)
          IntersectionPos(1:3) = PhotonProps%PhotonPos(1:3) + intersecDist*PhotonProps%PhotonDirection(1:3)
!          print*, NrOfThroughSides, SideID, IntersectionPos
        END IF ! foundHit
      END IF ! UsePhotonTriaTracking
    END IF  ! Mortar or regular side
  END DO LocSideLoop ! iLocSide=1,6

  IF ((NrOfThroughSides.EQ.0).AND.(.NOT.UsePhotonTriaTracking)) THEN
    ! Check if tracking fails to converge
    IterPhoton(1) = IterPhoton(1) + 1_8 ! Stop when MaxIterPhoton(2) is reached with this counter
    IF(IterPhoton(1).GE.MaxIterPhoton(1))THEN
      CALL StoreLostPhotonProperties(ElemID,TRIM(__FILE__),__LINE__,9999)
      Done = .TRUE.
      EXIT THEWHILELOOP
    END IF
    FallBack = .TRUE.
    LocSideLoop2: DO iLocSide=1,nlocSides
      TempSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide
      localSideID = SideInfo_Shared(SIDE_LOCALID,TempSideID)
      ! Side is not one of the 6 local sides
      IF (localSideID.LE.0) CYCLE
      NbElemID = SideInfo_Shared(SIDE_NBELEMID,TempSideID)
      IF (NbElemID.LT.0) THEN ! Mortar side
        CALL ABORT(__STAMP__,'Mortars not allowed for photon tracing!')
      ELSE  ! Regular side
        ! A) TriaTracking
        DO TriNum = 1,2
          ThroughSide = .FALSE.
          CALL PhotonThroughSideCheck3DFast(localSideID,ElemID,ThroughSide,TriNum)
          IF (ThroughSide) THEN
            CALL PhotonIntersectionWithSide(localSideID,ElemID,TriNum, IntersectionPos, PhotonLost)
            IF (PhotonLost) CYCLE
            NrOfThroughSides = NrOfThroughSides + 1
            LocSidesTemp(NrOfThroughSides) = localSideID
            TriNumTemp(NrOfThroughSides) = TriNum
            GlobSideTemp(NrOfThroughSides) = TempSideID
            ind = MAXLOC(ABS(PhotonProps%PhotonDirection),1)
            intersecDist = (IntersectionPos(ind) - PhotonProps%PhotonPos(ind))/PhotonProps%PhotonDirection(ind)
            InterPointSelect = (/.TRUE.,.TRUE./)
            InterPointSelectTemp(1:2,NrOfThroughSides) = InterPointSelect
            SideID = TempSideID
            LocalSide = localSideID
          END IF
        END DO
      END IF
    END DO LocSideLoop2
  END IF
!  
  TriNum = TriNumTemp(1)
  ! ----------------------------------------------------------------------------
  ! Additional treatment if particle did not cross any sides or it crossed multiple sides
  IF (NrOfThroughSides.NE.1) THEN
    ! 2c) If no sides are found, the particle is deleted (very rare case). If multiple possible sides are found, additional
    ! treatment is required, where the particle path is reconstructed (which side was crossed first) by comparing the ratio
    ! the determinants
    IF (NrOfThroughSides.EQ.0) THEN
      ! Particle appears to have not crossed any of the checked sides (NrOfThroughSides=0). Deleted!
      CALL StoreLostPhotonProperties(ElemID,TRIM(__FILE__),__LINE__,999)
      Done = .TRUE.
      EXIT THEWHILELOOP
    ELSE IF (NrOfThroughSides.GT.1) THEN
      IF(UsePhotonTriaTracking.OR.FallBack)THEN
        ! Use the slower search method if particle appears to have crossed more than one side (possible for irregular hexagons
        ! and in the case of mortar elements)
        SecondNrOfThroughSides = 0
        minRatio = 1E90
        oldElemIsMortar = .FALSE.
        DO ind2 = 1, NrOfThroughSides
          doCheckSide = .TRUE.
          ! Check if this side was already treated
          DO indSide = 2, 6
            IF((DoneLastElem(1,indSide).EQ.ElemID).AND. &
               (DoneLastElem(4,indSide).EQ.GlobSideTemp(ind2)).AND. &
               (DoneLastElem(3,indSide).EQ.TriNumTemp(ind2))) THEN
              doCheckSide = .FALSE.
            END IF
          END DO
          IF (doCheckSide) THEN
            IF (isMortarSideTemp(ind2)) THEN  ! Mortar side
              NbElemID = SideInfo_Shared(SIDE_ELEMID,GlobSideTemp(ind2))
              ! Get the determinant between the old and new particle position and the nodes of the triangle which was crossed
              CALL PhotonIntersectionWithSide(LocSidesTemp(ind2),NbElemID,TriNumTemp(ind2), IntersectionPosTemp, PhotonLost, .TRUE.)
              IF(PhotonLost)THEN
                ! Error in Photon TriaTracking! PhotonIntersectionWithSide() cannot determine intersection because photon is
                ! parallel to side
                CALL StoreLostPhotonProperties(ElemID,TRIM(__FILE__),__LINE__,999)
                Done = .TRUE.
                EXIT THEWHILELOOP
              END IF ! PhotonLost
              intersecDistVec(1:3) = IntersectionPosTemp(1:3) - PhotonProps%PhotonLastPos(1:3)
              intersecDist = DOT_PRODUCT(intersecDistVec, intersecDistVec)
              ! If the particle is inside the neighboring mortar element, it moved through this side
              ! Ratio is always negative since detM(=detLastPartPos) is negative or zero, i.e. maximum abs is wanted
              ! The closer the intersected side is to the last particle position the greater the absolute ratio will be
              IF (intersecDist.LT.minRatio) THEN
                IntersectionPos = IntersectionPosTemp
                minRatio = intersecDist
                SecondNrOfThroughSides = SecondNrOfThroughSides + 1
                SideID = GlobSideTemp(ind2)
                LocalSide = LocSidesTemp(ind2)
                TriNum = TriNumTemp(ind2)
                oldElemIsMortar = .TRUE.
              END IF
            ELSE  ! Regular side
              CALL PhotonIntersectionWithSide(LocSidesTemp(ind2),ElemID,TriNumTemp(ind2), IntersectionPosTemp, PhotonLost)
              IF(PhotonLost)THEN
                ! Error in Photon TriaTracking! PhotonIntersectionWithSide() cannot determine intersection because photon is
                ! parallel to side
                CALL StoreLostPhotonProperties(ElemID,TRIM(__FILE__),__LINE__,999)
                Done = .TRUE.
                EXIT THEWHILELOOP
              END IF ! PhotonLost
              intersecDistVec(1:3) = IntersectionPosTemp(1:3) - PhotonProps%PhotonLastPos(1:3)
              intersecDist = DOT_PRODUCT(intersecDistVec, intersecDistVec)
              IF (intersecDist.LT.minRatio) THEN
                IntersectionPos = IntersectionPosTemp
                minRatio = intersecDist
                SecondNrOfThroughSides = SecondNrOfThroughSides + 1
                SideID = GlobSideTemp(ind2)
                LocalSide = LocSidesTemp(ind2)
                TriNum = TriNumTemp(ind2)
                oldElemIsMortar = .FALSE.
              END IF
            END IF  ! isMortarSideTemp = T/F
          END IF  ! doCheckSide
        END DO  ! ind2 = 1, NrOfThroughSides
        ! Particle that went through multiple sides first, but did not cross any sides during the second check -> Deleted!
        IF (SecondNrOfThroughSides.EQ.0) THEN
          CALL StoreLostPhotonProperties(ElemID,TRIM(__FILE__),__LINE__,999)
          Done = .TRUE.
          EXIT THEWHILELOOP
        END IF
      ELSE
        ind2 = MINLOC(DistTemp(1:NrOfThroughSides),1)
        IF (DistTemp(ind2).LE.0.0) THEN
          IPWRITE(UNIT_StdOut,*) "NrOfThroughSides =", NrOfThroughSides, "DistTemp(1:NrOfThroughSides) =",DistTemp(1:NrOfThroughSides)
          CALL abort(__STAMP__,' ERROR: Side Distance is negative!')
        END IF
        SideID = GlobSideTemp(ind2)
        LocalSide = LocSidesTemp(ind2)
        TriNum = 1
        intersecDist = DistTemp(ind2)
        oldElemIsMortar = .FALSE.
        InterPointSelect = InterPointSelectTemp(1:2,ind2)
        IntersectionPos(1:3) = PhotonProps%PhotonLastPos(1:3) + intersecDist*PhotonProps%PhotonDirection(1:3)
      END IF
    END IF  ! NrOfThroughSides.EQ.0/.GT.1
  END IF  ! NrOfThroughSides.NE.1

  ! Dummy flag
  !IF((.NOT.UsePhotonTriaTracking).AND.(.NOT.FallBack)) TriNum=1

  ! Dummy flag: only if both flags are false, set TriNum=1
  IF(.NOT.(UsePhotonTriaTracking.OR.FallBack)) TriNum=1
  ! ----------------------------------------------------------------------------
  ! 3) In case of a boundary, perform the appropriate boundary interaction
  IF (SideInfo_Shared(SIDE_BCID,SideID).GT.0) THEN
    OldElemID=ElemID
    iPBC = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
    BCType = PartBound%TargetBoundCond(iPBC)
    SELECT CASE(BCType)
    CASE(1) !PartBound%OpenBC)
      IF(UsePhotonTriaTracking.OR.Fallback)THEN
        IF(NrOfThroughSides.LT.2)THEN
          CALL PhotonIntersectionWithSide(LocalSide,ElemID,TriNum, IntersectionPos, PhotonLost)
          IF(PhotonLost)THEN
            ! Error in open BC Photon TriaTracking! PhotonIntersectionWithSide() cannot determine intersection because photon is
            ! parallel to side
            CALL StoreLostPhotonProperties(ElemID,TRIM(__FILE__),__LINE__,999)
            Done = .TRUE.
            EXIT THEWHILELOOP
          END IF ! PhotonLost
        END IF ! NrOfThroughSides.LT.2
      END IF
      CALL CalcAbsoprtion(IntersectionPos(1:3),ElemID, DONE)
      IF (RadObservationPointMethod.EQ.1) THEN
        IF (PointInObsCone(IntersectionPos(1:3))) THEN
          IF (PhotonIntersectSensor(IntersectionPos(1:3), PhotonProps%PhotonDirection(1:3))) THEN
            RadObservation_Emission(PhotonProps%WaveLength) = RadObservation_Emission(PhotonProps%WaveLength) + PhotonProps%PhotonEnergy
            RadObservation_EmissionPart(PhotonProps%WaveLength) = RadObservation_EmissionPart(PhotonProps%WaveLength) + 1
          END IF
        END IF
      ELSE IF (RadObservationPointMethod.EQ.2) THEN
        IF (PhotonOnLineOfSight(PhotonProps%PhotonDirection(1:3))) THEN
          RadObservation_Emission(PhotonProps%WaveLength) = RadObservation_Emission(PhotonProps%WaveLength) + PhotonProps%PhotonEnergy
          RadObservation_EmissionPart(PhotonProps%WaveLength) = RadObservation_EmissionPart(PhotonProps%WaveLength) + 1
        END IF
      END IF
      DONE = .TRUE.

    CASE(2) ! PartBound%ReflectiveBC
      ! Backup photon direction for ray tracing
      IF(RadiationAbsorptionModel.EQ.0) PhotonProps%PhotonDirectionBeforeReflection(1:3) = PhotonProps%PhotonDirection(1:3)

      ! Check if specular of diffuse reflection
      IF (PartBound%PhotonSpecularReflection(iPBC)) THEN
        ! Specular reflection
        IF ((NrOfThroughSides.LT.2).AND.(UsePhotonTriaTracking.OR.Fallback)) THEN
          CALL PerfectPhotonReflection(LocalSide, ElemID, TriNum, IntersectionPos, .FALSE.)
        ELSE
          !TriNum=1
          CALL PerfectPhotonReflection(LocalSide, ElemID, TriNum, IntersectionPos, .TRUE.)
        END IF
      ELSE
        ! Diffuse reflection
        IF ((NrOfThroughSides.LT.2).AND.(UsePhotonTriaTracking.OR.Fallback)) THEN
          CALL DiffusePhotonReflection(LocalSide, ElemID, TriNum, IntersectionPos, .FALSE.)
        ELSE
          CALL DiffusePhotonReflection(LocalSide, ElemID, TriNum, IntersectionPos, .TRUE.)
        END IF
      END IF

      ! Check if ray tracing is active
      IF(RadiationAbsorptionModel.EQ.0)THEN
        CALL CalcAbsoprtion(IntersectionPos(1:3), ElemID, DONE, before = .TRUE.)
        IF (.NOT.DONE) THEN
          CALL CalcWallAbsoprtion(IntersectionPos(1:3),SideID, DONE, RayForceAbsorption)
          CALL CalcAbsoprtion(IntersectionPos(1:3), ElemID, DONE, before = .FALSE.)
        END IF ! .NOT.DONE
      ELSE
        CALL CalcAbsoprtion(IntersectionPos(1:3), ElemID, DONE)
        IF (.NOT.DONE) CALL CalcWallAbsoprtion(IntersectionPos(1:3),SideID, DONE)
      END IF ! RadiationAbsorptionModel.EQ.0

    CASE(3) ! PartBound%PeriodicBC
      IF((NrOfThroughSides.LT.2).AND.(UsePhotonTriaTracking.OR.Fallback))THEN
        CALL PhotonIntersectionWithSide(LocalSide,ElemID,TriNum, IntersectionPos, PhotonLost)      
        IF(PhotonLost)THEN
          ! Error in periodic Photon TriaTracking! PhotonIntersectionWithSide() cannot determine intersection because photon is
          ! parallel to side
          CALL StoreLostPhotonProperties(ElemID,TRIM(__FILE__),__LINE__,999)
          Done = .TRUE.
          EXIT THEWHILELOOP
        END IF ! PhotonLost
      END IF ! NrOfThroughSides.LT.2
      CALL CalcAbsoprtion(IntersectionPos(1:3), ElemID, DONE)
      ! Move photon across periodic BC
      CALL PeriodicPhotonBC(LocalSide,ElemID,TriNum,IntersectionPos,.TRUE.,SideID)
    CASE DEFAULT
      CALL abort(__STAMP__,' ERROR: PartBound not associated!. (unknown case)',BCType,999.)
    END SELECT !PartBound%MapToPartBC(BC(SideID)


    IF ((BCType.EQ.2).OR.(BCType.EQ.10)) THEN
      DoneLastElem(:,:) = 0
      LastInterPointSelect = .FALSE.
    ELSE
      DO ind2= 5, 1, -1
        DoneLastElem(:,ind2+1) = DoneLastElem(:,ind2)
      END DO
      DoneLastElem(1,1) = OldElemID
      DoneLastElem(2,1) = LocalSide
      DoneLastElem(3,1) = TriNum
      DoneLastElem(4,1) = SideID
      LastInterPointSelect= InterPointSelect
    END IF
  ELSE  ! BC(SideID).LE.0
    DO ind2= 5, 1, -1
      DoneLastElem(:,ind2+1) = DoneLastElem(:,ind2)
    END DO
    DoneLastElem(1,1) = ElemID
    DoneLastElem(2,1) = LocalSide
    DoneLastElem(3,1) = TriNum
    DoneLastElem(4,1) = SideID
    LastInterPointSelect = InterPointSelect

    IF (oldElemIsMortar) THEN
      ElemID = SideInfo_Shared(SIDE_ELEMID,SideID)
      CALL PhotonIntersectionWithSide(LocalSide,ElemID,TriNum, IntersectionPos, PhotonLost, IsMortar=.TRUE.)
    ELSE
      ! For bilinear tracing, the intersection point is already calculated by the tracking algorithm itself
      IF(UsePhotonTriaTracking) CALL PhotonIntersectionWithSide(LocalSide,ElemID,TriNum, IntersectionPos, PhotonLost)
    END IF

    ! Check if lost during intersection
    IF(PhotonLost)THEN
      ! Error in Photon TriaTracking! PhotonIntersectionWithSide() cannot determine intersection because photon is parallel to side
      CALL StoreLostPhotonProperties(ElemID,TRIM(__FILE__),__LINE__,999)
      Done = .TRUE.
      EXIT THEWHILELOOP
    ELSE
      ! Absorption
      IF (oldElemIsMortar) THEN
        CALL CalcAbsoprtion(IntersectionPos(1:3),DoneLastElem(1,1), DONE)
      ELSE
        CALL CalcAbsoprtion(IntersectionPos(1:3),ElemID, DONE)
        ElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)
      END IF
    END IF ! PhotonLost
  END IF  ! BC(SideID).GT./.LE. 0

  ! Check if output to PartStateBoundary is activated
  IF(PhotonModeBPO.EQ.2)THEN
    CALL StoreBoundaryParticleProperties(0,&
         999,&
         IntersectionPos,&
         UNITVECTOR(PhotonProps%PhotonDirection(1:3)),(/0.0,0.0,1.0/),&
         iPartBound=0,&
         mode=2,&
         MPF_optIN=0.0,&
         Velo_optIN=PhotonProps%PhotonDirection(1:3))
  END IF ! PhotonModeBPO.EQ.2

  IF (ElemID.LT.1) CALL abort(__STAMP__ ,'ERROR: Element not defined! Please increase the size of the halo region (HaloEpsVelo)!')
END DO THEWHILELOOP ! .NOT.PartisDone


END SUBROUTINE PhotonTriaTracking



SUBROUTINE Photon2DSymTracking()
!===================================================================================================================================
! Routine for tracking of moving particles and boundary interaction using triangulated sides.
! 1) Loop over all particles that are still inside
!    2) Perform tracking until the particle is considered "done" (either localized or deleted)
!       2a) Perform a check based on the determinant of (3x3) matrix of the vectors from the particle position to the nodes of each
!           triangle (ParticleInsideQuad3D)
!       2b) If particle is not within the given element in a), the side through which the particle went is determined by checking
!           each side of the element (ParticleThroughSideCheck3DFast)
!       2c) If no sides are found, the particle is deleted (very rare case). If multiple possible sides are found, additional
!           treatment is required, where the particle path is reconstructed (which side was crossed first) by comparing the ratio
!           the determinants
!    3) In case of a boundary, determine the intersection and perform the appropriate boundary interaction (GetBoundaryInteraction)
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Boundary_Vars,      ONLY:PartBound
USE MOD_Photon_TrackingVars,         ONLY:PhotonProps
USE MOD_RadiationTrans_Vars,         ONLY:RadObservation_Emission, RadObservationPointMethod,RadObservation_EmissionPart
USE MOD_Photon_TrackingTools,        ONLY:CalcAbsoprtion, CalcWallAbsoprtion, DiffusePhotonReflection2D, PointInObsCone
USE MOD_Photon_TrackingTools,        ONLY:PhotonIntersectionWithSide2D, RotatePhotonIn2DPlane, PerfectPhotonReflection2D
USE MOD_Photon_TrackingTools,        ONLY:PhotonIntersectSensor, PhotonOnLineOfSight
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: NbElemID, ind, nbSideID, nMortarElems, BCType, NblocSideID, OrigElem
INTEGER                          :: ElemID, OldElemID, LocalSide, NrOfThroughSides, localSideID, nlocSides
INTEGER                          :: SideID, TempSideID, iLocSide, correctSide, LastSide
INTEGER                          :: LocSidesTemp(1:6), GlobSideTemp(1:6)
LOGICAL                          :: oldElemIsMortar, isMortarSideTemp(1:6), isLastSide, ThroughSide, Done
REAL                             :: IntersectionPos(1:6), IntersectionPosTemp(1:6,4), DistanceTemp(1:6), Distance
!===================================================================================================================================
Done = .FALSE.
ElemID = PhotonProps%ElemID
OrigElem = ElemID
SideID = 0
LastSide = 0
! 1) Loop tracking until Photon is considered "done" (either absorbed or deleted)
DO WHILE (.NOT.Done)
  oldElemIsMortar = .FALSE.
  NrOfThroughSides = 0
  LocSidesTemp = 0
  GlobSideTemp = 0
  DistanceTemp = 0.0
  IntersectionPosTemp = 0.0
  isMortarSideTemp = .FALSE.
  nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)
  DO iLocSide=1,nlocSides
    isLastSide = .FALSE.
    TempSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide
    localSideID = SideInfo_Shared(SIDE_LOCALID,TempSideID)
    ! Side is not one of the 6 local sides
    IF (localSideID.LE.0) CYCLE
    IF (SideIsSymSide(TempSideID)) CYCLE
    IF (SideInfo_Shared(SIDE_BCID,TempSideID).GT.0) THEN
      IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,TempSideID))).EQ.11) CYCLE
    END IF
    IF (LastSide.EQ.TempSideID) isLastSide = .TRUE.
    NbElemID = SideInfo_Shared(SIDE_NBELEMID,TempSideID)
    IF (NbElemID.LT.0) THEN ! Mortar side
      nMortarElems = MERGE(4,2,SideInfo_Shared(SIDE_NBELEMID,TempSideID).EQ.-1)
      DO ind = 1, nMortarElems
        isLastSide = .FALSE.
        nbSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide + ind
        NbElemID = SideInfo_Shared(SIDE_NBELEMID,nbSideID)
        ! If small mortar element not defined, abort. Every available information on the compute-node is kept in shared memory, so
        ! no way to recover it during runtime
        IF (LastSide.EQ.nbSideID) isLastSide = .TRUE.
        IF (NbElemID.LT.1) CALL abort(__STAMP__,'Small mortar element not defined!',ElemID)
        ! For small mortar sides, SIDE_NBSIDEID contains the SideID of the corresponding big mortar side
        nbSideID = SideInfo_Shared(SIDE_NBSIDEID,nbSideID)
        NblocSideID =  SideInfo_Shared(SIDE_LOCALID,nbSideID)
        ThroughSide = .FALSE.
        CALL PhotonIntersectionWithSide2D(NblocSideID,NbElemID,ThroughSide,IntersectionPos, isLastSide,Distance)
        IF (ThroughSide) THEN
          ! Store the information for this side for future checks, if this side was already treated
          oldElemIsMortar = .TRUE.
          NrOfThroughSides = NrOfThroughSides + 1
          LocSidesTemp(NrOfThroughSides) = NblocSideID
          GlobSideTemp(NrOfThroughSides) = nbSideID
          DistanceTemp(NrOfThroughSides) = Distance
          IntersectionPosTemp(1:3,NrOfThroughSides) = IntersectionPos(1:3)
          isMortarSideTemp(NrOfThroughSides)  = .TRUE.
          SideID = nbSideID
          LocalSide = NblocSideID
        END IF
      END DO
    ELSE  ! Regular side
      ThroughSide = .FALSE.
      CALL PhotonIntersectionWithSide2D(localSideID,ElemID,ThroughSide,IntersectionPos, isLastSide, Distance)
      IF (ThroughSide) THEN
        NrOfThroughSides = NrOfThroughSides + 1
        SideID = TempSideID
        LocalSide = localSideID
        LocSidesTemp(NrOfThroughSides) = localSideID
        GlobSideTemp(NrOfThroughSides) = TempSideID
        DistanceTemp(NrOfThroughSides) = Distance
        IntersectionPosTemp(1:3,NrOfThroughSides) = IntersectionPos(1:3)
      END IF
    END IF  ! Mortar or regular side
  END DO  ! iLocSide=1,4
  ! ----------------------------------------------------------------------------
  ! Addition treatment if particle did not cross any sides or it crossed multiple sides
  IF (NrOfThroughSides.NE.1) THEN
    ! 2c) If no sides are found, the particle is deleted (very rare case). If multiple possible sides are found, additional
    ! treatment is required, where the particle path is reconstructed (which side was crossed first) by comparing the ratio
    ! the determinants
    IF (NrOfThroughSides.EQ.0) THEN
      ! Particle appears to have not crossed any of the checked sides. Deleted!
      IPWRITE(*,*) 'Error in Photon 2DAxisTracking! Photon lost. Element:', ElemID
      IPWRITE(*,*) 'LastPos: ', PhotonProps%PhotonLastPos(1:3)
      IPWRITE(*,*) 'Pos:     ', PhotonProps%PhotonPos(1:3)
      IPWRITE(*,*) 'Direction:', PhotonProps%PhotonDirection(1:3)
      IPWRITE(*,*) 'OrigElem:', OrigElem
      IPWRITE(*,*) 'Photon deleted!'
      Done = .TRUE.
      EXIT
    ELSE IF (NrOfThroughSides.GT.1) THEN
      correctSide = MINLOC(DistanceTemp(1:NrOfThroughSides),DIM=1)
      oldElemIsMortar = isMortarSideTemp(correctSide)
      SideID = GlobSideTemp(correctSide)
      LocalSide = LocSidesTemp(correctSide)
      IntersectionPos(1:3) = IntersectionPosTemp(1:3,correctSide)
    END IF
  END IF  ! NrOfThroughSides.NE.1
  ! ----------------------------------------------------------------------------
  ! 3) In case of a boundary, perform the appropriate boundary interaction
  IF (SideInfo_Shared(SIDE_BCID,SideID).GT.0) THEN
    OldElemID=ElemID
    BCType = PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)))
    SELECT CASE(BCType)
    CASE(1) !PartBound%OpenBC)
      CALL CalcAbsoprtion(IntersectionPos(1:3),ElemID, DONE)
      DONE = .TRUE.
      IF (RadObservationPointMethod.EQ.1) THEN
        IF (PointInObsCone(IntersectionPos(1:3))) THEN
          IF (PhotonIntersectSensor(IntersectionPos(1:3), PhotonProps%PhotonDirection(1:3))) THEN
            RadObservation_Emission(PhotonProps%WaveLength) = RadObservation_Emission(PhotonProps%WaveLength) + PhotonProps%PhotonEnergy
            RadObservation_EmissionPart(PhotonProps%WaveLength) = RadObservation_EmissionPart(PhotonProps%WaveLength) + 1
          END IF
        END IF
      ELSE IF (RadObservationPointMethod.EQ.2) THEN
        IF (PhotonOnLineOfSight(PhotonProps%PhotonDirection(1:3))) THEN
          RadObservation_Emission(PhotonProps%WaveLength) = RadObservation_Emission(PhotonProps%WaveLength) + PhotonProps%PhotonEnergy
          RadObservation_EmissionPart(PhotonProps%WaveLength) = RadObservation_EmissionPart(PhotonProps%WaveLength) + 1
        END IF
      END IF
      CYCLE
    CASE(2)
      IF (PartBound%PhotonSpecularReflection(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)))) THEN
        CALL CalcAbsoprtion(IntersectionPos(1:3),ElemID, DONE)
        IF (.NOT.DONE) CALL PerfectPhotonReflection2D(LocalSide,ElemID, IntersectionPos)
      ELSE
        CALL CalcAbsoprtion(IntersectionPos(1:3),ElemID, DONE)
        IF (.NOT.DONE) CALL CalcWallAbsoprtion(IntersectionPos(1:3),SideID, DONE)
        IF (.NOT.DONE) CALL DiffusePhotonReflection2D(LocalSide,ElemID, IntersectionPos)
      END IF
      LastSide = SideID
    CASE DEFAULT
      CALL abort(__STAMP__,' ERROR: PartBound not associated!. (unknown case)',999,999.)
    END SELECT !PartBound%MapToPartBC(BC(SideID)
  ELSE  ! BC(SideID).LE.0
    IF (oldElemIsMortar) THEN
      CALL CalcAbsoprtion(IntersectionPos(1:3),ElemID, DONE)
      ElemID = SideInfo_Shared(SIDE_ELEMID,SideID)
      LastSide = SideID
    ELSE
      CALL CalcAbsoprtion(IntersectionPos(1:3),ElemID, DONE)
      ElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)
      LastSide = SideInfo_Shared(SIDE_NBSIDEID,SideID)
    END IF
  END IF  ! BC(SideID).GT./.LE. 0

  IF (.NOT.DONE) CALL RotatePhotonIn2DPlane(IntersectionPos(1:3))
  IF (ElemID.LT.1) THEN
    CALL abort(__STAMP__ ,'ERROR: Element not defined! Please increase the size of the halo region (HaloEpsVelo)!')
  END IF
END DO  ! .NOT.PartisDone
END SUBROUTINE Photon2DSymTracking

SUBROUTINE PhotonComputeBiLinearIntersection(isHit,SideID, Dist,PointSelect,OldPointSelect,ElemCheck_Opt)
!===================================================================================================================================
! Compute the Intersection with planar surface, improved version by
! Haselbacher, A.; Najjar, F. M. & Ferry, J. P., An efficient and robust particle-localization algorithm for unstructured grids
! Journal of Computational Physics, Elsevier BV, 2007, 225, 2198-2213
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: EpsMach
USE MOD_Utils                  ,ONLY: QuadraticSolver
USE MOD_Mesh_Tools             ,ONLY: GetCNSideID, GetCNElemID
USE MOD_Particle_Mesh_Vars     ,ONLY: SideInfo_Shared
USE MOD_Particle_Surfaces_Vars ,ONLY: BaseVectors0,BaseVectors1,BaseVectors2,BaseVectors3,SideNormVec,epsilonTol!,BaseVectorsScale
USE MOD_Particle_Surfaces      ,ONLY: CalcNormAndTangBilinear
USE MOD_Photon_TrackingVars    ,ONLY: PhotonProps
#ifdef CODE_ANALYZE
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
USE MOD_Particle_Tracking_Vars ,ONLY: MPIRankOut
USE MOD_Mesh_Vars              ,ONLY: NGeo
#endif /*CODE_ANALYZE*/
#if USE_MPI
!USE MOD_Mesh_Vars              ,ONLY: BC
#endif /*USE_MPI*/
USE MOD_Particle_Intersection   ,ONLY: ComputeXi,ComputeSurfaceDistance2
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: SideID
LOGICAL,INTENT(IN),OPTIONAL       :: ElemCheck_Opt
REAL, INTENT(OUT)                 :: Dist
LOGICAL,INTENT(OUT)               :: PointSelect(2)
LOGICAL,INTENT(IN)                :: OldPointSelect(2)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)               :: isHit
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(4)                 :: a1,a2
REAL,DIMENSION(1:3,1:4)           :: BiLinearCoeff,NormalCoeff
REAL                              :: A,B,C,alpha, xitild, etatild
REAL                              :: xi(2),eta(2),t(2),scaleFac
INTEGER                           :: CNSideID,InterType,nRoot
LOGICAL                           :: ElemCheck
!===================================================================================================================================

! set alpha to minus one // no intersection
alpha    = -1.0
xitild   = -2.0
etatild  = -2.0
Dist = 0.0
isHit    = .FALSE.
PointSelect = .FALSE.
CNSideID = GetCNSideID(SideID)
! compute initial vectors
BiLinearCoeff(:,1) = 0.25*BaseVectors3(:,SideID)
BiLinearCoeff(:,2) = 0.25*BaseVectors1(:,SideID)
BiLinearCoeff(:,3) = 0.25*BaseVectors2(:,SideID)
BiLinearCoeff(:,4) = 0.25*BaseVectors0(:,SideID)

#ifdef CODE_ANALYZE
  IF(MPIRANKOUT.EQ.MyRank)THEN
    WRITE(UNIT_stdout,'(110("-"))')
    WRITE(UNIT_stdout,'(A)') '     | Output of bilinear intersection equation constants: '
    WRITE(UNIT_stdout,'(A,3(1X,G0))') '     | SideNormVec  : ',SideNormVec(1:3,CNSideID)
    WRITE(UNIT_stdout,'(A,4(1X,G0))') '     | BilinearCoeff: ',BilinearCoeff(1,1:4)
    WRITE(UNIT_stdout,'(A,4(1X,G0))') '     | BilinearCoeff: ',BilinearCoeff(2,1:4)
    WRITE(UNIT_stdout,'(A,4(1X,G0))') '     | BilinearCoeff: ',BilinearCoeff(3,1:4)
    WRITE(UNIT_stdout,'(A,3(1X,G0))') '     | Beziercontrolpoint1: ',BezierControlPoints3D(:,0,0,SideID)
    WRITE(UNIT_stdout,'(A,3(1X,G0))') '     | Beziercontrolpoint2: ',BezierControlPoints3D(:,NGeo,0,SideID)
    WRITE(UNIT_stdout,'(A,3(1X,G0))') '     | Beziercontrolpoint3: ',BezierControlPoints3D(:,0,NGeo,SideID)
    WRITE(UNIT_stdout,'(A,3(1X,G0))') '     | Beziercontrolpoint4: ',BezierControlPoints3D(:,NGeo,NGeo,SideID)
  END IF
#endif /*CODE_ANALYZE*/

! Check if the site can be encountered. Both vectors are already normalized
scaleFac = DOT_PRODUCT(PhotonProps%PhotonDirection,SideNormVec(1:3,CNSideID))
IF (ABS(scaleFac).LT.epsilontol) RETURN

! Haselbacher et al. define d = d - r_p
BiLinearCoeff(:,4) = BiLinearCoeff(:,4) - PhotonProps%PhotonPos(1:3)

! Calculate component normal to ray
NormalCoeff(:,1) = BiLinearCoeff(:,1) - SUM(BiLinearCoeff(:,1)*PhotonProps%PhotonDirection(:))*PhotonProps%PhotonDirection
NormalCoeff(:,2) = BiLinearCoeff(:,2) - SUM(BiLinearCoeff(:,2)*PhotonProps%PhotonDirection(:))*PhotonProps%PhotonDirection
NormalCoeff(:,3) = BiLinearCoeff(:,3) - SUM(BiLinearCoeff(:,3)*PhotonProps%PhotonDirection(:))*PhotonProps%PhotonDirection
NormalCoeff(:,4) = BiLinearCoeff(:,4) - SUM(BiLinearCoeff(:,4)*PhotonProps%PhotonDirection(:))*PhotonProps%PhotonDirection

! A1 is X_xz = X_z - X_x
A1(:) = NormalCoeff(3,:) - NormalCoeff(1,:)
! A2 is X_yz = X_z - X_y
A2(:) = NormalCoeff(3,:) - NormalCoeff(2,:)

! Bring into quadratic form
A = a1(1)*a2(3) - a2(1)*a1(3)
B = a1(1)*a2(4) - a2(1)*a1(4) + a1(2)*a2(3) - a2(2)*a1(3)
C = a1(2)*a2(4) - a2(2)*a1(4)

! Scale with <PartTraj.,NormVec>^2 and cell-scale (~area) for getting coefficients at least approx. in the order of 1
!scaleFac = scaleFac**2 * BaseVectorsScale(SideID) !<...>^2 * cell-scale
!scaleFac = 1./scaleFac
!A = A * scaleFac
!B = B * scaleFac
!C = C * scaleFac

CALL QuadraticSolver(A,B,C,nRoot,Eta(1),Eta(2))

#ifdef CODE_ANALYZE
  IF(MPIRANKOUT.EQ.MyRank)THEN
    WRITE(UNIT_stdout,'(A)') '     | Output after QuadraticSolver: '
    WRITE(UNIT_stdout,'(A,I0,A,2(1X,G0))') '     | number of root: ',nRoot,' | Eta: ',Eta(1:2)
  END IF
#endif /*CODE_ANALYZE*/

! nRoot equals the number of possible intersections with the bilinear surface. However, only values between [-1,1] are valid
SELECT CASE(nRoot)
  ! No intersection
  CASE(0) ! nRoot = 0
    RETURN

  ! One possible intersection
  CASE(1) ! nRoot = 1
#ifdef CODE_ANALYZE
    IF(MPIRANKOUT.EQ.MyRank)THEN
      WRITE(UNIT_stdout,'(A)') '     | nRoot = 1 '
    END IF
#endif /*CODE_ANALYZE*/
    ! Check if eta is valid
    IF (ABS(eta(1)).LE.1.0) THEN
      ! check for Xi only, if eta is possible
      xi(1) = ComputeXi(eta(1),A1=A1,A2=A2)

      IF (Xi(1).EQ.HUGE(1.)) THEN
        RETURN
        IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' Both denominators zero when calculating Xi in bilinear intersection'
        IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' global SideID:      ', SideID
        IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' global ElemID:      ', SideInfo_Shared(SIDE_ELEMID,SideID)
        IPWRITE(UNIT_stdOut,'(I0,A,3(1X,ES25.17E3))') ' LastPhotonPos:   ', PhotonProps%PhotonPos(1:3)
        IPWRITE(UNIT_stdOut,'(I0,A,3(1X,ES25.17E3))') ' PhotonDirection:       ', PhotonProps%PhotonDirection(1:3)        
        CALL ABORT(__STAMP__,'Invalid intersection with bilinear side!',SideID)
      END IF

      IF( ABS(xi(1)).LE.1.0) THEN
        ! compute alpha only with valid xi and eta
        t(1) = ComputeSurfaceDistance2(BiLinearCoeff,xi(1),eta(1),PhotonProps%PhotonDirection)
        IF ((t(1).GT.ABS(EpsMach)).AND.(.NOT.OldPointSelect(1))) THEN
          alpha   = t(1)
          xitild  = xi(1)
          etatild = eta(1)
          isHit   = .TRUE.
          Dist = t(1)
          PointSelect(1) = .TRUE.
#ifdef CODE_ANALYZE
          IF(MPIRANKOUT.EQ.MyRank)THEN
            WRITE(UNIT_stdout,'(A,G0,A,G0)') '     | t(1): ',t(1),' | epsilonTolerance: ',epsilontol
          END IF
#endif /*CODE_ANALYZE*/
          ! This is the only possible intersection, so we are done
          RETURN
        ELSE ! t is not in range
          RETURN
        END IF
      ELSE ! xi not in range
        RETURN
      END IF ! ABS(xi(1)).LE.1.0
    ELSE ! eta not in range
      RETURN
    END IF ! ABS(eta(1)).LE.1.0

  CASE(2) ! nRoot = 2
#ifdef CODE_ANALYZE
    IF(MPIRANKOUT.EQ.MyRank)THEN
      WRITE(UNIT_stdout,'(A)') '     | nRoot = 2 '
    END IF
#endif /*CODE_ANALYZE*/
    InterType = 0
    t(:)      =-1.

    ! Check if eta(1)) is valid
    IF (ABS(eta(1)).LE.1.0) THEN
      ! check for Xi only, if eta is possible
      xi(1) = ComputeXi(eta(1),A1=A1,A2=A2)

      IF (Xi(1).EQ.HUGE(1.)) THEN
        return
!        IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' Both denominators zero when calculating Xi in bilinear intersection'
!        IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' global SideID:      ', SideID
!        IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' global ElemID:      ', SideInfo_Shared(SIDE_ELEMID,SideID)
!        IPWRITE(UNIT_stdOut,'(I0,A,3(1X,ES25.17E3))') ' LastPhotonPos:   ', PhotonProps%PhotonPos(1:3)
!        IPWRITE(UNIT_stdOut,'(I0,A,3(1X,ES25.17E3))') ' PhotonDirection:       ', PhotonProps%PhotonDirection(1:3)
!        CALL ABORT(__STAMP__,'Invalid intersection with bilinear side!',SideID)
      END IF

      IF( ABS(xi(1)).LE.1.0) THEN
        ! compute alpha only with valid xi and eta
        t(1) = ComputeSurfaceDistance2(BiLinearCoeff,xi(1),eta(1),PhotonProps%PhotonDirection)

#ifdef CODE_ANALYZE
        IF(MPIRANKOUT.EQ.MyRank)THEN
          WRITE(UNIT_stdout,'(A,G0,A,G0)') '     | xi: ',xi(1),' | t: ',t(1)
        END IF
#endif /*CODE_ANALYZE*/

        IF (t(1).GT.ABS(EpsMach)) THEN
          InterType = InterType+1
          isHit     = .TRUE.
#ifdef CODE_ANALYZE
          IF(MPIRANKOUT.EQ.MyRank)THEN
            WRITE(UNIT_stdout,'(A,E15.8,A,E15.8)') '     | t(1): ',t(1),' | epsilonTolerance: ',epsilontol
          END IF
#endif /*CODE_ANALYZE*/
        END IF
      END IF ! ABS(xi(1)).LE.1.0
    END IF ! ABS(eta(1)).LE.1.0

    ! Check if eta(2) is valid
    IF (ABS(eta(2)).LE.1.0) THEN
      ! check for Xi only, if eta is possible
      xi(2) = ComputeXi(eta(2),A1=A1,A2=A2)

      IF (Xi(2).EQ.HUGE(1.)) THEN
        return
!        IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' Both denominators zero when calculating Xi in bilinear intersection'
!        IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' global SideID:      ', SideID
!        IPWRITE(UNIT_stdOut,'(I0,A,I0)') ' global ElemID:      ', SideInfo_Shared(SIDE_ELEMID,SideID)
!        IPWRITE(UNIT_stdOut,'(I0,A,3(1X,ES25.17E3))') ' LastPhotonPos:   ', PhotonProps%PhotonPos(1:3)
!        IPWRITE(UNIT_stdOut,'(I0,A,3(1X,ES25.17E3))') ' PhotonDirection:       ', PhotonProps%PhotonDirection(1:3)
!        CALL ABORT(__STAMP__,'Invalid intersection with bilinear side!',SideID)
      END IF

      IF( ABS(xi(2)).LE.1.0) THEN
        ! compute alpha only with valid xi and eta
        t(2) = ComputeSurfaceDistance2(BiLinearCoeff,xi(2),eta(2),PhotonProps%PhotonDirection)

#ifdef CODE_ANALYZE
        IF(MPIRANKOUT.EQ.MyRank)THEN          
          WRITE(UNIT_stdout,'(A,G0,A,G0)') '     | xi: ',xi(2),' | t: ',t(2)
        END IF
#endif /*CODE_ANALYZE*/

        IF (t(2).GT.ABS(EpsMach)) THEN
          InterType = InterType+2
          isHit     = .TRUE.
#ifdef CODE_ANALYZE
          IF(MPIRANKOUT.EQ.MyRank)THEN
            WRITE(UNIT_stdout,'(A,E15.8,A,E15.8)') '     | t(2): ',t(2),' | epsilonTolerance: ',epsilontol
          END IF
#endif /*CODE_ANALYZE*/
        END IF
      END IF ! ABS(xi(2)).LE.1.0
    END IF ! ABS(eta(2)).LE.1.0

    SELECT CASE(InterType)
      ! No intersection found, return
      CASE(0)
        RETURN

      ! First intersection is only hit
      CASE(1)
        IF (.NOT.OldPointSelect(1)) THEN
          alpha  =t  (1)
          xitild =xi (1)
          etatild=eta(1)
          Dist = t(1)
          PointSelect(1)=.TRUE.
        ELSE
          isHit = .FALSE.
          RETURN
        END IF
      ! Second intersection is only hit
      CASE(2)
        IF (.NOT.OldPointSelect(2)) THEN
          alpha  =t  (2)
          xitild =xi (2)
          etatild=eta(2)
          Dist = t(2)
          PointSelect(2)=.TRUE.
        ELSE
          isHit = .FALSE.
          RETURN
        END IF

      ! Two intersections found, decide on the correct one
      CASE(3)
        ! If side is a BC side, take only the intersection encountered first
        IF (SideInfo_Shared(SIDE_BCID,SideID).GT.0) THEN
          ! Check if the element is supposed to be checked
          ElemCheck = .FALSE.
          IF(PRESENT(ElemCheck_Opt))THEN
            ElemCheck = ElemCheck_Opt
          END IF

          IF(ElemCheck)THEN
            alpha  =-1
            xitild =-2
            etatild=-2
          ELSE
            ! Apparently we don't care about the direction of the PartTrajectory
            IF((ABS(t(1)).LT.ABS(t(2))).AND.(.NOT.OldPointSelect(1)))THEN
              alpha  =t  (1)
              xitild =xi (1)
              etatild=eta(1)
              Dist = t(1)
              PointSelect(1)=.TRUE.
            ELSE IF (.NOT.OldPointSelect(2)) THEN
              alpha  =t  (2)
              xitild =xi (2)
              etatild=eta(2)
              Dist = t(2)
              PointSelect(2)=.TRUE.
            ELSE
              isHit = .FALSE.
            END IF
          END IF
        ! Inner side with double intersection, particle leaves and enters element
        ELSE
          alpha  =-1
          xitild = 0.
          etatild= 0.
          isHit  = .FALSE.
        END IF
    END SELECT ! InterType
END SELECT ! nRoot

END SUBROUTINE PhotonComputeBiLinearIntersection

END MODULE MOD_Photon_Tracking
