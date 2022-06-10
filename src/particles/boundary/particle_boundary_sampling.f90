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

MODULE MOD_Particle_Boundary_Sampling
!===================================================================================================================================
!! Determines how particles interact with a given boundary condition
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE DefineParametersParticleBoundarySampling
  MODULE PROCEDURE DefineParametersParticleBoundarySampling
END INTERFACE

INTERFACE InitParticleBoundarySampling
  MODULE PROCEDURE InitParticleBoundarySampling
END INTERFACE

INTERFACE FinalizeParticleBoundarySampling
  MODULE PROCEDURE FinalizeParticleBoundarySampling
END INTERFACE

INTERFACE WriteSurfSampleToHDF5
  MODULE PROCEDURE WriteSurfSampleToHDF5
END INTERFACE

PUBLIC::DefineParametersParticleBoundarySampling
PUBLIC::InitParticleBoundarySampling
PUBLIC::WriteSurfSampleToHDF5, WriteSurfSampleChemToHDF5
PUBLIC::FinalizeParticleBoundarySampling
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for particle surface sampling
!==================================================================================================================================
SUBROUTINE DefineParametersParticleBoundarySampling()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Particle Boundary Sampling")
CALL prms%CreateIntOption(      'DSMC-nSurfSample'  , 'Define polynomial degree of particle BC sampling. Default: NGeo', '1')
CALL prms%CreateLogicalOption(  'CalcSurfaceImpact' , 'Sample average impact energy of particles for each species (trans, rot, '//&
                                                      'vib), impact vector and angle.','.FALSE.')
END SUBROUTINE DefineParametersParticleBoundarySampling


SUBROUTINE InitParticleBoundarySampling()
!===================================================================================================================================
! Initialization of particle boundary sampling
! default: use for sampling same polynomial degree as NGeo
! 1) all procs identify surfaces for sampling on the node (plus halo region)
! 2) the compute-node leaders communicate the number of sampling surfaces
! 3) every proc holds a copy of the complete sampling region on the node
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Basis                   ,ONLY: LegendreGaussNodesAndWeights
USE MOD_DSMC_Symmetry           ,ONLY: DSMC_2D_CalcSymmetryArea, DSMC_1D_CalcSymmetryArea
USE MOD_Mesh_Vars               ,ONLY: NGeo,nBCs,BoundaryName
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfOnNode
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfSample,dXiEQ_SurfSample,PartBound,XiEQ_SurfSample
USE MOD_Particle_Boundary_Vars  ,ONLY: nComputeNodeSurfSides,nComputeNodeSurfTotalSides,nComputeNodeSurfOutputSides
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfBC,SurfBCName
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfTotalSides
USE MOD_Particle_Boundary_Vars  ,ONLY: GlobalSide2SurfSide,SurfSide2GlobalSide
USE MOD_SurfaceModel_Vars       ,ONLY: nPorousBC
USE MOD_Particle_Boundary_Vars  ,ONLY: CalcSurfaceImpact
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfSideArea,SurfSampSize
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallState
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallPumpCapacity
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallImpactEnergy
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallImpactVector
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallImpactAngle
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallImpactNumber
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared,NodeCoords_Shared
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemSideNodeID_Shared
USE MOD_Particle_Surfaces       ,ONLY: EvaluateBezierPolynomialAndGradient
USE MOD_Particle_Surfaces_Vars  ,ONLY: BezierControlPoints3D
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
USE MOD_Particle_Vars           ,ONLY: nSpecies,VarTimeStep
USE MOD_Particle_Vars           ,ONLY: Symmetry
USE MOD_ReadInTools             ,ONLY: GETINT,GETLOGICAL,GETINTARRAY
#if USE_MPI
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemInfo_Shared
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_LEADERS_SURF,mySurfRank
USE MOD_MPI_Shared_Vars         ,ONLY: myComputeNodeRank,nComputeNodeProcessors
USE MOD_Particle_Mesh_Vars      ,ONLY: nNonUniqueGlobalSides
!USE MOD_Particle_Mesh_Vars      ,ONLY: offsetComputeNodeElem,nComputeNodeElems
USE MOD_MPI_Shared_Vars         ,ONLY: myLeaderGroupRank,nLeaderGroupProcs
USE MOD_Particle_Boundary_Vars  ,ONLY: GlobalSide2SurfSide_Shared,GlobalSide2SurfSide_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfSide2GlobalSide_Shared,SurfSide2GlobalSide_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfSideArea_Shared,SurfSideArea_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallState_Shared,SampWallState_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallPumpCapacity_Shared,SampWallPumpCapacity_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallImpactEnergy_Shared,SampWallImpactEnergy_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallImpactVector_Shared,SampWallImpactVector_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallImpactAngle_Shared,SampWallImpactAngle_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallImpactNumber_Shared,SampWallImpactNumber_Shared_Win
USE MOD_Particle_MPI_Boundary_Sampling,ONLY: InitSurfCommunication
USE MOD_Particle_Boundary_Vars  ,ONLY: nComputeNodeInnerBCs
#else
USE MOD_MPI_Shared_Vars         ,ONLY: mySurfRank
USE MOD_Particle_Mesh_Vars      ,ONLY: nComputeNodeSides
USE MOD_Particle_Boundary_Vars  ,ONLY: nOutputSides
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                :: iBC
INTEGER                                :: iSide,firstSide,lastSide,iSurfSide,GlobalSideID
INTEGER                                :: nSurfSidesProc
INTEGER                                :: offsetSurfTotalSidesProc
INTEGER,ALLOCATABLE                    :: GlobalSide2SurfSideProc(:,:)
!INTEGER,ALLOCATABLE                    :: SurfSide2GlobalSideProc(:,:)
CHARACTER(20)                          :: hilf
CHARACTER(LEN=255),ALLOCATABLE         :: BCName(:)
! surface area
INTEGER                                :: SideID,ElemID,CNElemID,LocSideID
INTEGER                                :: p,q,iSample,jSample
INTEGER                                :: TriNum, Node1, Node2
REAL                                   :: area,nVal
REAL,DIMENSION(2,3)                    :: gradXiEta3D
REAL,DIMENSION(:),ALLOCATABLE          :: Xi_NGeo,wGP_NGeo
REAL                                   :: XiOut(1:2),E,F,G,D,tmp1,tmpI2,tmpJ2
REAL                                   :: xNod, zNod, yNod, Vector1(3), Vector2(3), nx, ny, nz
#if USE_MPI
INTEGER                                :: offsetSurfSidesProc
INTEGER                                :: GlobalElemID,GlobalElemRank
INTEGER                                :: sendbuf,recvbuf
INTEGER                                :: NbGlobalElemID, NbElemRank, NbLeaderID, nSurfSidesTmp
#endif /*USE_MPI*/
INTEGER                                :: NbGlobalSideID
!===================================================================================================================================

! Get input parameters
SWRITE(UNIT_stdOut,'(A)') ' INIT SURFACE SAMPLING ...'

WRITE(UNIT=hilf,FMT='(I0)') NGeo
nSurfSample = GETINT('DSMC-nSurfSample',TRIM(hilf))

IF((nSurfSample.GT.1).AND.(TrackingMethod.EQ.TRIATRACKING)) &
  CALL abort(__STAMP__,'nSurfSample cannot be >1 if TrackingMethod = triatracking')

! Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z) and angle
CalcSurfaceImpact = GETLOGICAL('CalcSurfaceImpact')

! Allocate shared array for surf sides
#if USE_MPI
CALL Allocate_Shared((/3,nNonUniqueGlobalSides/),GlobalSide2SurfSide_Shared_Win,GlobalSide2SurfSide_Shared)
CALL MPI_WIN_LOCK_ALL(0,GlobalSide2SurfSide_Shared_Win,IERROR)
GlobalSide2SurfSide => GlobalSide2SurfSide_Shared
#else
ALLOCATE(GlobalSide2SurfSide(1:3,1:nComputeNodeSides))
#endif /*USE_MPI*/

! only CN root nullifies
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /* USE_MPI*/
  GlobalSide2SurfSide = -1.
#if USE_MPI
END IF

CALL BARRIER_AND_SYNC(GlobalSide2SurfSide_Shared_Win,MPI_COMM_SHARED)
#endif /* USE_MPI*/

! get number of BC-Sides
#if USE_MPI
! NO HALO REGION REDUCTION
firstSide = INT(REAL( myComputeNodeRank   *nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1)*nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))
ALLOCATE(GlobalSide2SurfSideProc(1:3,firstSide:lastSide))
        !,SurfSide2GlobalSideProc(1:3,1         :INT(nNonUniqueGlobalSides/REAL(nComputeNodeProcessors))))
#else
firstSide = 1
lastSide  = nComputeNodeSides
ALLOCATE(GlobalSide2SurfSideProc(1:3,1:nComputeNodeSides))
        !,SurfSide2GlobalSideProc(1:3,1:nComputeNodeSides))
#endif /*USE_MPI*/

GlobalSide2SurfSideProc    = -1
!SurfSide2GlobalSideProc    = -1
nComputeNodeSurfSides      = 0
nSurfSidesProc             = 0

! check every BC side
DO iSide = firstSide,lastSide
  ! ignore non-BC sides
  IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE

#if USE_MPI
  ! ignore sides outside of halo region
  IF (ElemInfo_Shared(ELEM_HALOFLAG,SideInfo_Shared(SIDE_ELEMID,iSide)).EQ.0) CYCLE
#endif /*USE_MPI*/

  ! count number of reflective and inner BC sides
  IF ( (PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,iSide))).EQ.PartBound%ReflectiveBC) .OR. &
       (PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,iSide))).EQ.PartBound%RotPeriodicBC) ) THEN
    nSurfSidesProc = nSurfSidesProc + 1
    ! check if element for this side is on the current compute-node
    ! IF ((SideInfo_Shared(SIDE_ID,iSide).GT.ElemInfo_Shared(ELEM_FIRSTSIDEIND,offsetComputeNodeElem+1))                  .AND. &
    !     (SideInfo_Shared(SIDE_ID,iSide).LE.ElemInfo_Shared(ELEM_LASTSIDEIND ,offsetComputeNodeElem+nComputeNodeElems))) THEN
!    IF ((iSide.GE.(ElemInfo_Shared(ELEM_FIRSTSIDEIND,offsetComputeNodeElem+1)+1))                  .AND. &
!        (iSide.LE.ElemInfo_Shared(ELEM_LASTSIDEIND ,offsetComputeNodeElem+nComputeNodeElems))) THEN
!      nComputeNodeSurfSides  = nComputeNodeSurfSides + 1
!    END IF

    ! TODO: Add another check to determine the surface side in halo_eps from current proc. Node-wide halo can become quite large with
    !       with 128 procs!

    ! Write local mapping from Side to Surf side. The rank is already correct, the offset must be corrected by the proc offset later
    GlobalSide2SurfSideProc(SURF_SIDEID,iSide) = nSurfSidesProc
#if USE_MPI
    GlobalSide2SurfSideProc(SURF_RANK  ,iSide) = ElemInfo_Shared(ELEM_RANK,SideInfo_Shared(SIDE_ELEMID,iSide))
    ! get global Elem ID
    GlobalElemID   = SideInfo_Shared(SIDE_ELEMID,iSide)
    GlobalElemRank = ElemInfo_Shared(ELEM_RANK,GlobalElemID)
    ! running on one node, everything belongs to us
    IF (nLeaderGroupProcs.EQ.1) THEN
      GlobalSide2SurfSideProc(SURF_LEADER,iSide) = myLeaderGroupRank
    ELSE
      ! find the compute node
      GlobalSide2SurfSideProc(SURF_LEADER,iSide) = INT(GlobalElemRank/nComputeNodeProcessors)
        END IF
#else
    GlobalSide2SurfSideProc(SURF_RANK  ,iSide) = 0
    GlobalSide2SurfSideProc(SURF_LEADER,iSide) = GlobalSide2SurfSideProc(SURF_RANK,iSide)
#endif /*USE_MPI*/

#if USE_MPI
    ! check if element for this side is on the current compute-node. Alternative version to the check above
    IF (GlobalSide2SurfSideProc(SURF_LEADER,iSide).EQ.myLeaderGroupRank) THEN
#endif /*USE_MPI*/
      nComputeNodeSurfSides  = nComputeNodeSurfSides + 1
#if USE_MPI
    END IF
#endif /*USE_MPI*/
  END IF ! reflective side
END DO

! Find CN global number of total surf sides and write Side to Surf Side mapping into shared array
#if USE_MPI
sendbuf = nSurfSidesProc - nComputeNodeSurfSides
recvbuf = 0
CALL MPI_EXSCAN(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,iError)
offsetSurfTotalSidesProc   = recvbuf
! last proc knows CN total number of BC elems
sendbuf = offsetSurfTotalSidesProc + nSurfSidesProc - nComputeNodeSurfSides
CALL MPI_BCAST(sendbuf,1,MPI_INTEGER,nComputeNodeProcessors-1,MPI_COMM_SHARED,iError)
nComputeNodeSurfTotalSides = sendbuf

! Find CN global number of local surf sides and write Side to Surf Side mapping into shared array
sendbuf = nComputeNodeSurfSides
recvbuf = 0
CALL MPI_EXSCAN(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,iError)
offsetSurfSidesProc   = recvbuf
! last proc knows CN total number of BC elems
sendbuf = offsetSurfSidesProc + nComputeNodeSurfSides
CALL MPI_BCAST(sendbuf,1,MPI_INTEGER,nComputeNodeProcessors-1,MPI_COMM_SHARED,iError)
nComputeNodeSurfSides = sendbuf
nComputeNodeSurfTotalSides = nComputeNodeSurfTotalSides + nComputeNodeSurfSides

! increment SURF_SIDEID by offset
nSurfSidesTmp = 0
DO iSide = firstSide,lastSide
  IF (GlobalSide2SurfSideProc(SURF_SIDEID,iSide).EQ.-1) CYCLE

  ! sort compute-node local sides first
  IF (GlobalSide2SurfSideProc(SURF_LEADER,iSide).EQ.myLeaderGroupRank) THEN
    nSurfSidesTmp = nSurfSidesTmp + 1

    GlobalSide2SurfSide(:          ,iSide) = GlobalSide2SurfSideProc(:,iSide)
    GlobalSide2SurfSide(SURF_SIDEID,iSide) = nSurfSidesTmp + offsetSurfSidesProc
  END IF
END DO

nSurfSidesTmp = 0
DO iSide = firstSide,lastSide
  IF (GlobalSide2SurfSideProc(SURF_SIDEID,iSide).EQ.-1) CYCLE

  ! sampling sides in halo region follow at the end
  IF (GlobalSide2SurfSideProc(SURF_LEADER,iSide).NE.myLeaderGroupRank) THEN
    nSurfSidesTmp = nSurfSidesTmp + 1

    GlobalSide2SurfSide(:          ,iSide) = GlobalSide2SurfSideProc(:,iSide)
    GlobalSide2SurfSide(SURF_SIDEID,iSide) = nSurfSidesTmp + nComputeNodeSurfSides + offsetSurfTotalSidesProc
  END IF
END DO
#else
offsetSurfTotalSidesProc  = 0
nComputeNodeSurfTotalSides = nSurfSidesProc
GlobalSide2SurfSide(:,firstSide:lastSide) = GlobalSide2SurfSideProc(:,firstSide:lastSide)
#endif /*USE_MPI*/

! Build inverse mapping
#if USE_MPI
CALL Allocate_Shared((/3,nComputeNodeSurfTotalSides/),SurfSide2GlobalSide_Shared_Win,SurfSide2GlobalSide_Shared)
CALL MPI_WIN_LOCK_ALL(0,SurfSide2GlobalSide_Shared_Win,IERROR)
SurfSide2GlobalSide => SurfSide2GlobalSide_Shared

DO iSide = firstSide,lastSide
  IF (GlobalSide2SurfSideProc(SURF_SIDEID,iSide).EQ.-1) CYCLE

  SurfSide2GlobalSide(:          ,GlobalSide2SurfSide(SURF_SIDEID,iSide)) = GlobalSide2SurfSide(:,iSide)
  SurfSide2GlobalSide(SURF_SIDEID,GlobalSide2SurfSide(SURF_SIDEID,iSide)) = iSide
END DO

CALL BARRIER_AND_SYNC(GlobalSide2SurfSide_Shared_Win,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(SurfSide2GlobalSide_Shared_Win,MPI_COMM_SHARED)
#else
ALLOCATE(SurfSide2GlobalSide(1:1,1:nComputeNodeSurfTotalSides))
!SurfSide2GlobalSide = SurfSide2GlobalSideProc(:,1:nComputeNodeSurfTotalSides)
DO iSide = firstSide,lastSide
  IF (GlobalSide2SurfSide(SURF_SIDEID,iSide).EQ.-1) CYCLE
  SurfSide2GlobalSide(SURF_SIDEID,GlobalSide2SurfSide(SURF_SIDEID,iSide)) =iSide
END DO
#endif /*USE_MPI*/

! Determine the number of surface output sides (inner BCs are not counted twice)
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
  nComputeNodeInnerBCs = 0
#endif /*USE_MPI*/
  nComputeNodeSurfOutputSides = 0
  DO iSurfSide = 1,nComputeNodeSurfSides
    GlobalSideID = SurfSide2GlobalSide(SURF_SIDEID,iSurfSide)
    ! Check if the surface side has a neighbor (and is therefore an inner BCs)
    IF(SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID).GT.0) THEN
      ! Abort inner BC + Mortar! (too complex and confusing to implement)
      ! This test catches large Mortar sides, i.e.,  SideInfo_Shared(SIDE_NBELEMID,NonUniqueGlobalSideID) gives the 2 or 4
      ! connecting small Mortar sides. It is assumed that inner BC result in being flagged as a "SurfSide" and therefore are checked
      ! here.
      IF(SideInfo_Shared(SIDE_LOCALID,GlobalSideID).EQ.-1)THEN
        IPWRITE(UNIT_StdOut,'(I12,A,I0)')   " NonUniqueGlobalSideID                               = ",GlobalSideID
        IPWRITE(UNIT_StdOut,'(I12,A,I0)')   " SideInfo_Shared(SIDE_LOCALID,NonUniqueGlobalSideID) = ",&
            SideInfo_Shared(SIDE_LOCALID,GlobalSideID)
        IPWRITE(UNIT_StdOut,'(I12,A,I0,A)') " SideInfo_Shared(SIDE_ELEMID,NonUniqueGlobalSideID)  = ",&
            SideInfo_Shared(SIDE_ELEMID,GlobalSideID)," (GlobalElemID)"
        CALL abort(__STAMP__,'Inner BC + Mortar is not implemented!')
      END IF
      ! Only add the side with the smaller index
      NbGlobalSideID = SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID)
      IF(GlobalSideID.GT.NbGlobalSideID)THEN
#if USE_MPI
        !--- switcheroo check 1 of 2: Non-HALO sides
        ! Only required for sampling on the larger NonUniqueGlobalSideID of the two sides of the inner BC
        ! Count larger inner BCs as these may have to be sent to a different leader processor
        NbGlobalElemID = SideInfo_Shared(SIDE_ELEMID,NbGlobalSideID)
        NbElemRank = ElemInfo_Shared(ELEM_RANK,NbGlobalElemID)
        NbLeaderID = INT(NbElemRank/nComputeNodeProcessors)
        IF(NbLeaderID.NE.INT(myRank/nComputeNodeProcessors))THEN
          nComputeNodeInnerBCs(1) = nComputeNodeInnerBCs(1) + 1
        END IF
#endif
        CYCLE! Skip sides with the larger index
      END IF
    END IF
    nComputeNodeSurfOutputSides = nComputeNodeSurfOutputSides + 1
  END DO
#if USE_MPI
  !--- switcheroo check 2 of 2: HALO sides
  ! Count number of inner BC in halo region
  ! Only required for sampling on the larger NonUniqueGlobalSideID of the two sides of the inner BC
  DO iSurfSide = nComputeNodeSurfSides+1, nComputeNodeSurfTotalSides
    GlobalSideID = SurfSide2GlobalSide(SURF_SIDEID,iSurfSide)
    ! Check if the surface side has a neighbor (and is therefore an inner BCs)
    IF(SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID).GT.0) THEN
      ! Only add the side with the smaller index
      IF(GlobalSideID.GT.SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID))THEN
        ! Count larger inner BCs as these may have to be sent to a different leader processor
        nComputeNodeInnerBCs(2) = nComputeNodeInnerBCs(2) + 1
      END IF
    END IF
  END DO ! iSurfSide = nComputeNodeSurfSides+1, nComputeNodeSurfTotalSides
END IF
#endif

! free temporary arrays
DEALLOCATE(GlobalSide2SurfSideProc)
!DEALLOCATE(SurfSide2GlobalSideProc)

! flag if there is at least one surf side on the node (sides in halo region do also count)
SurfOnNode = MERGE(.TRUE.,.FALSE.,nComputeNodeSurfTotalSides.GT.0)

!> Energy + Force + nSpecies
SurfSampSize = SAMPWALL_NVARS+nSpecies
IF(VarTimeStep%UseVariableTimeStep) SurfSampSize = SurfSampSize + 1

!> Leader communication
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
  CALL InitSurfCommunication()
END IF
! The leaders are synchronized at this point, but behind the other procs. nSurfTotalSides is only required when compiled without
! MPI, so perform latency hiding by postponing synchronization
#else
mySurfRank      = 0
nSurfTotalSides = nComputeNodeSurfTotalSides
nOutputSides    = nComputeNodeSurfOutputSides
#endif /* USE_MPI */

! surface sampling array do not need to be allocated if there are no sides within halo_eps range
IF(.NOT.SurfOnNode) RETURN

! create boundary name mapping for surfaces SurfaceBC number mapping
nSurfBC = 0
ALLOCATE(BCName(1:nBCs))
DO iBC=1,nBCs
  BCName=''
END DO
DO iBC=1,nBCs
  ! inner side (can be just in the name list from preproc although already sorted out)
  IF (PartBound%MapToPartBC(iBC).EQ.-1) CYCLE

  ! count number of reflective BCs
  IF ( (PartBound%TargetBoundCond(PartBound%MapToPartBC(iBC)).EQ.PartBound%ReflectiveBC).OR. &
       (PartBound%TargetBoundCond(PartBound%MapToPartBC(iBC)).EQ.PartBound%RotPeriodicBC)    ) THEN
    nSurfBC         = nSurfBC + 1
    BCName(nSurfBC) = BoundaryName(iBC)
  END IF
END DO

IF (nSurfBC.GE.1) THEN
ALLOCATE(SurfBCName(1:nSurfBC))
  DO iBC=1,nSurfBC
    SurfBCName(iBC) = BCName(iBC)
  END DO
END IF
DEALLOCATE(BCName)

! allocate everything. Caution: SideID is surfSideID, not global ID
!> First local arrays for boundary sampling
!>>> Pointer structure type ruins possibility of ALLREDUCE, need separate arrays for everything
ALLOCATE(SampWallState(1:SurfSampSize,1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides))
SampWallState = 0.
!
IF(nPorousBC.GT.0) THEN
  ALLOCATE(SampWallPumpCapacity(1:nComputeNodeSurfTotalSides))
  SampWallPumpCapacity = 0.
END IF
! Sampling of impact energy for each species (trans, rot, vib, elec), impact vector (x,y,z) and angle
IF (CalcSurfaceImpact) THEN
  ALLOCATE( SampWallImpactEnergy(1:nSpecies,1:4,1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides) &
          , SampWallImpactVector(1:nSpecies,1:3,1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides) &
          , SampWallImpactAngle (1:nSpecies,    1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides) &
          , SampWallImpactNumber(1:nSpecies,    1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides))
  SampWallImpactEnergy = 0.
  SampWallImpactVector = 0.
  SampWallImpactAngle  = 0.
  SampWallImpactNumber = 0.
END IF

#if USE_MPI
!> Then shared arrays for boundary sampling
CALL Allocate_Shared((/SurfSampSize,nSurfSample,nSurfSample,nComputeNodeSurfTotalSides/),SampWallState_Shared_Win,SampWallState_Shared)
CALL MPI_WIN_LOCK_ALL(0,SampWallState_Shared_Win,IERROR)
! #else
! ALLOCATE(SampWallState_Shared(1:SurfSampSize,1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides))
! SampWallState_Shared = 0.
! #endif /*USE_MPI*/

! #if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
! #endif /*USE_MPI*/
  SampWallState_Shared = 0.
! #if USE_MPI
END IF
CALL BARRIER_AND_SYNC(SampWallState_Shared_Win,MPI_COMM_SHARED)
!
IF(nPorousBC.GT.0) THEN
  CALL Allocate_Shared((/nComputeNodeSurfTotalSides/),SampWallPumpCapacity_Shared_Win,SampWallPumpCapacity_Shared)
  CALL MPI_WIN_LOCK_ALL(0,SampWallPumpCapacity_Shared_Win,IERROR)
  IF (myComputeNodeRank.EQ.0) THEN
    SampWallPumpCapacity_Shared = 0.
  END IF
  CALL BARRIER_AND_SYNC(SampWallPumpCapacity_Shared_Win,MPI_COMM_SHARED)
END IF
! Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z) and angle
IF (CalcSurfaceImpact) THEN
  ! nSpecies*4
  CALL Allocate_Shared((/nSpecies,4,nSurfSample,nSurfSample,nComputeNodeSurfTotalSides/),SampWallImpactEnergy_Shared_Win,SampWallImpactEnergy_Shared)
  CALL MPI_WIN_LOCK_ALL(0,SampWallImpactEnergy_Shared_Win,IERROR)
  ! nSpecies*3
  CALL Allocate_Shared((/nSpecies,3,nSurfSample,nSurfSample,nComputeNodeSurfTotalSides/),SampWallImpactVector_Shared_Win,SampWallImpactVector_Shared)
  CALL MPI_WIN_LOCK_ALL(0,SampWallImpactVector_Shared_Win,IERROR)
  ! nSpecies
  CALL Allocate_Shared((/nSpecies,  nSurfSample,nSurfSample,nComputeNodeSurfTotalSides/),SampWallImpactAngle_Shared_Win,SampWallImpactAngle_Shared)
  CALL MPI_WIN_LOCK_ALL(0,SampWallImpactAngle_Shared_Win,IERROR)
  CALL Allocate_Shared((/nSpecies,  nSurfSample,nSurfSample,nComputeNodeSurfTotalSides/),SampWallImpactNumber_Shared_Win,SampWallImpactNumber_Shared)
  CALL MPI_WIN_LOCK_ALL(0,SampWallImpactNumber_Shared_Win,IERROR)
  IF (myComputeNodeRank.EQ.0) THEN
    SampWallImpactEnergy_Shared = 0.
    SampWallImpactVector_Shared = 0.
    SampWallImpactAngle_Shared  = 0.
    SampWallImpactNumber_Shared = 0.
  END IF
  CALL BARRIER_AND_SYNC(SampWallImpactEnergy_Shared_Win,MPI_COMM_SHARED)
  CALL BARRIER_AND_SYNC(SampWallImpactVector_Shared_Win,MPI_COMM_SHARED)
  CALL BARRIER_AND_SYNC(SampWallImpactAngle_Shared_Win ,MPI_COMM_SHARED)
  CALL BARRIER_AND_SYNC(SampWallImpactNumber_Shared_Win,MPI_COMM_SHARED)
END IF
!#else
!CALL ABORT(__STAMP__,'Still needs to be implemented for MPI=OFF')
#endif /*USE_MPI*/

! Surf sides are shared, array calculation can be distributed
#if USE_MPI
CALL Allocate_Shared((/nSurfSample,nSurfSample,nComputeNodeSurfTotalSides/),SurfSideArea_Shared_Win,SurfSideArea_Shared)
CALL MPI_WIN_LOCK_ALL(0,SurfSideArea_Shared_Win,IERROR)
SurfSideArea => SurfSideArea_Shared

firstSide = INT(REAL( myComputeNodeRank   *nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1)*nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))
#else
ALLOCATE(SurfSideArea(1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides))

firstSide = 1
lastSide  = nSurfTotalSides
#endif /*USE_MPI*/

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
  SurfSideArea=0.
#if USE_MPI
END IF
CALL BARRIER_AND_SYNC(SurfSideArea_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/

! Calculate equidistant surface points
ALLOCATE(XiEQ_SurfSample(0:nSurfSample))
dXiEQ_SurfSample =2./REAL(nSurfSample)
DO q=0,nSurfSample
  XiEQ_SurfSample(q) = dXiEQ_SurfSample * REAL(q) - 1.
END DO

! get interpolation points and weights
ALLOCATE( Xi_NGeo( 0:NGeo)  &
        , wGP_NGeo(0:NGeo) )
CALL LegendreGaussNodesAndWeights(NGeo,Xi_NGeo,wGP_NGeo)

! compute area of sub-faces
tmp1=dXiEQ_SurfSample/2.0 !(b-a)/2

DO iSide = firstSide,LastSide
  ! get global SideID. This contains only nonUniqueSide, no special mortar treatment required
  SideID = SurfSide2GlobalSide(SURF_SIDEID,iSide)

  IF (TrackingMethod.EQ.TRIATRACKING) THEN
    ElemID    = SideInfo_Shared(SIDE_ELEMID ,SideID)
    CNElemID  = GetCNElemID(ElemID)
    LocSideID = SideInfo_Shared(SIDE_LOCALID,SideID)
    area = 0.
    xNod = NodeCoords_Shared(1,ElemSideNodeID_Shared(1,LocSideID,CNElemID)+1)
    yNod = NodeCoords_Shared(2,ElemSideNodeID_Shared(1,LocSideID,CNElemID)+1)
    zNod = NodeCoords_Shared(3,ElemSideNodeID_Shared(1,LocSideID,CNElemID)+1)
    IF(Symmetry%Order.EQ.3) THEN
      DO TriNum = 1,2
        Node1 = TriNum+1     ! normal = cross product of 1-2 and 1-3 for first triangle
        Node2 = TriNum+2     !          and 1-3 and 1-4 for second triangle
        Vector1(1) = NodeCoords_Shared(1,ElemSideNodeID_Shared(Node1,LocSideID,CNElemID)+1) - xNod
        Vector1(2) = NodeCoords_Shared(2,ElemSideNodeID_Shared(Node1,LocSideID,CNElemID)+1) - yNod
        Vector1(3) = NodeCoords_Shared(3,ElemSideNodeID_Shared(Node1,LocSideID,CNElemID)+1) - zNod
        Vector2(1) = NodeCoords_Shared(1,ElemSideNodeID_Shared(Node2,LocSideID,CNElemID)+1) - xNod
        Vector2(2) = NodeCoords_Shared(2,ElemSideNodeID_Shared(Node2,LocSideID,CNElemID)+1) - yNod
        Vector2(3) = NodeCoords_Shared(3,ElemSideNodeID_Shared(Node2,LocSideID,CNElemID)+1) - zNod
        nx = - Vector1(2) * Vector2(3) + Vector1(3) * Vector2(2) !NV (inwards)
        ny = - Vector1(3) * Vector2(1) + Vector1(1) * Vector2(3)
        nz = - Vector1(1) * Vector2(2) + Vector1(2) * Vector2(1)
        nVal = SQRT(nx*nx + ny*ny + nz*nz)
        area = area + nVal/2.
      END DO
      SurfSideArea(1,1,iSide) = area
    ELSE IF(Symmetry%Order.EQ.2) THEN
      SurfSideArea(1,1,iSide) = DSMC_2D_CalcSymmetryArea(LocSideID, CNElemID)
    ELSE IF(Symmetry%Order.EQ.1) THEN
      SurfSideArea(1,1,iSide) = DSMC_1D_CalcSymmetryArea(LocSideID, CNElemID)
    END IF
  ELSE ! TrackingMethod.NE.TRIATRACKING
    ! call here stephens algorithm to compute area
    DO jSample=1,nSurfSample
      DO iSample=1,nSurfSample
        area=0.
        tmpI2=(XiEQ_SurfSample(iSample-1)+XiEQ_SurfSample(iSample))/2. ! (a+b)/2
        tmpJ2=(XiEQ_SurfSample(jSample-1)+XiEQ_SurfSample(jSample))/2. ! (a+b)/2
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
        SurfSideArea(iSample,jSample,iSide) = area
      END DO ! iSample=1,nSurfSample
    END DO ! jSample=1,nSurfSample
  END IF
END DO ! iSide = firstSide,lastSide

#if USE_MPI
CALL BARRIER_AND_SYNC(SurfSideArea_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/

! get the full area of all surface sides
area=0.

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
  DO iSide = 1,nComputeNodeSurfSides
    area = area + SUM(SurfSideArea(:,:,iSide))
  END DO ! iSide = 1,nComputeNodeSurfTotalSides
#if USE_MPI
  ! surf leaders communicate total surface area
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,area,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEADERS_SURF,IERROR)
END IF

! surf leaders inform the other procs on their node
CALL MPI_BCAST(area,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/

! de-allocate temporary interpolation points and weights
DEALLOCATE(Xi_NGeo,wGP_NGeo)

#if USE_MPI
! Delayed synchronization
CALL MPI_BCAST(nSurfTotalSides,1,MPI_INTEGER,0,MPI_COMM_SHARED,iError)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

IF (mySurfRank.EQ.0) THEN
#endif
  WRITE(UNIT_StdOut,'(A,I8)')       ' | Number of sampling sides:           '    , nSurfTotalSides
  WRITE(UNIT_StdOut,'(A,ES10.4E2)') ' | Surface-Area:                         ', Area
  WRITE(UNIT_stdOut,'(A)') ' INIT SURFACE SAMPLING DONE'
#if USE_MPI
END IF
#endif

END SUBROUTINE InitParticleBoundarySampling


SUBROUTINE WriteSurfSampleToHDF5(MeshFileName,OutputTime)
!===================================================================================================================================
!> write the final values of the surface sampling to a HDF5 state file
!> additional performs all the final required computations
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: ProjectName
USE MOD_DSMC_Vars               ,ONLY: MacroSurfaceVal,MacroSurfaceSpecVal, CollisMode
USE MOD_HDF5_Output             ,ONLY: WriteAttributeToHDF5,WriteArrayToHDF5,WriteHDF5Header
USE MOD_IO_HDF5
USE MOD_MPI_Shared_Vars         ,ONLY: mySurfRank
USE MOD_SurfaceModel_Vars       ,ONLY: nPorousBC
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfSample,CalcSurfaceImpact
USE MOD_Particle_Boundary_Vars  ,ONLY: nOutputSides
USE MOD_Particle_boundary_Vars  ,ONLY: nComputeNodeSurfOutputSides,offsetComputeNodeSurfOutputSide
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfBC,SurfBCName, PartBound
USE MOD_Particle_Vars           ,ONLY: nSpecies
#if USE_MPI
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfTotalSides
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_LEADERS_SURF
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: MeshFileName
REAL,INTENT(IN)                      :: OutputTime
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)                  :: FileName,FileString,Statedummy
CHARACTER(LEN=255)                  :: H5_Name
CHARACTER(LEN=255)                  :: NodeTypeTemp
CHARACTER(LEN=255)                  :: SpecID, PBCID
CHARACTER(LEN=255),ALLOCATABLE      :: Str2DVarNames(:)
INTEGER                             :: nVar2D, nVar2D_Spec, nVar2D_Total, nVarCount, iSpec, iPBC
REAL                                :: tstart,tend
!===================================================================================================================================

#if USE_MPI
! Return if not a sampling leader
IF (MPI_COMM_LEADERS_SURF.EQ.MPI_COMM_NULL) RETURN
CALL MPI_BARRIER(MPI_COMM_LEADERS_SURF,iERROR)

! Return if no sampling sides
IF (nSurfTotalSides      .EQ.0) RETURN
#endif /*USE_MPI*/

IF (mySurfRank.EQ.0) THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE DSMCSurfSTATE TO HDF5 FILE...'
  tstart=LOCALTIME()
END IF

FileName   = TIMESTAMP(TRIM(ProjectName)//'_DSMCSurfState',OutputTime)
FileString = TRIM(FileName)//'.h5'

! Create dataset attribute "SurfVarNames"
nVar2D      = 5
nVar2D_Spec = 1

! Sampling of impact energy for each species (trans, rot, vib, elec), impact vector (x,y,z), angle, total number, number per second: Add 10 variables
IF (CalcSurfaceImpact) nVar2D_Spec = nVar2D_Spec + 10

IF (nPorousBC.GT.0)    nVar2D = nVar2D + nPorousBC

IF (ANY(PartBound%UseAdaptedWallTemp)) nVar2D = nVar2D + 1

nVar2D_Total = nVar2D + nVar2D_Spec*nSpecies

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
#if USE_MPI
IF (mySurfRank.EQ.0) THEN
#endif
  CALL OpenDataFile(FileString,create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
  Statedummy = 'DSMCSurfState'

  ! Write file header
  CALL WriteHDF5Header(Statedummy,File_ID)
  CALL WriteAttributeToHDF5(File_ID,'DSMC_nSurfSample',1,IntegerScalar=nSurfSample)
  CALL WriteAttributeToHDF5(File_ID,'DSMC_nSpecies'   ,1,IntegerScalar=nSpecies)
  CALL WriteAttributeToHDF5(File_ID,'DSMC_CollisMode' ,1,IntegerScalar=CollisMode)
  CALL WriteAttributeToHDF5(File_ID,'MeshFile'        ,1,StrScalar=(/TRIM(MeshFileName)/))
  CALL WriteAttributeToHDF5(File_ID,'Time'            ,1,RealScalar=OutputTime)
  CALL WriteAttributeToHDF5(File_ID,'BC_Surf'         ,nSurfBC,StrArray=SurfBCName)
  CALL WriteAttributeToHDF5(File_ID,'N',1             ,IntegerScalar=nSurfSample)
  NodeTypeTemp='VISU'
  CALL WriteAttributeToHDF5(File_ID,'NodeType'        ,1,StrScalar=(/NodeTypeTemp/))

  ALLOCATE(Str2DVarNames(1:nVar2D_Total))
  Str2DVarNames(:) = ''
  nVarCount        = 1
  DO iSpec = 1,nSpecies
    WRITE(SpecID,'(I3.3)') iSpec
    CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Spec'//TRIM(SpecID)//'_SimPartPerIter')
    ! Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z) and angle
    IF(CalcSurfaceImpact)THEN
      ! Add average impact energy for each species (trans, rot, vib)
      CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Spec'//TRIM(SpecID)//'_ImpactEnergyTrans')
      CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Spec'//TRIM(SpecID)//'_ImpactEnergyRot')
      CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Spec'//TRIM(SpecID)//'_ImpactEnergyVib')
      CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Spec'//TRIM(SpecID)//'_ImpactEnergyElec')
      ! Add average impact vector for each species (x,y,z)
      CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Spec'//TRIM(SpecID)//'_ImpactVectorX')
      CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Spec'//TRIM(SpecID)//'_ImpactVectorY')
      CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Spec'//TRIM(SpecID)//'_ImpactVectorZ')
      ! Add average impact angle for each species
      CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Spec'//TRIM(SpecID)//'_ImpactAngle')
      ! Add number of impacts
      CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Spec'//TRIM(SpecID)//'_ImpactNumber')
      ! Add number of impacts per second
      CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Spec'//TRIM(SpecID)//'_ImpactFlux')
    END IF ! CalcSurfaceImpact

  END DO ! iSpec=1,nSpecies

  ! fill varnames for total values
  CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Total_ForcePerAreaX')
  CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Total_ForcePerAreaY')
  CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Total_ForcePerAreaZ')
  CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Total_HeatFlux')
  CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Total_SimPartPerIter')

  IF(nPorousBC.GT.0) THEN
    DO iPBC = 1, nPorousBC
      WRITE(PBCID,'(I2.2)') iPBC
      CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'PorousBC'//TRIM(PBCID)//'_PumpCapacity')
    END DO
  END IF

  IF (ANY(PartBound%UseAdaptedWallTemp)) CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Wall_Temperature')

  CALL WriteAttributeToHDF5(File_ID,'VarNamesSurface',nVar2D_Total,StrArray=Str2DVarNames)

  CALL CloseDataFile()
  DEALLOCATE(Str2DVarNames)
#if USE_MPI
END IF

CALL MPI_BARRIER(MPI_COMM_LEADERS_SURF,iERROR)

CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_LEADERS_SURF)
#else
CALL OpenDataFile(FileString,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
#endif

nVarCount=0
WRITE(H5_Name,'(A)') 'SurfaceData'
ASSOCIATE (&
      nVar2D_Total         => INT(nVar2D_Total,IK)                    , &
      nSurfSample          => INT(nSurfSample,IK)                     , &
      nGlobalSides         => INT(nOutputSides,IK)                    , &
      nLocalSides          => INT(nComputeNodeSurfOutputSides,IK)     , &
      offsetSurfSide       => INT(offsetComputeNodeSurfOutputSide,IK) , &
      nVar2D_Spec          => INT(nVar2D_Spec,IK)                     , &
      nVar2D               => INT(nVar2D,IK))
  DO iSpec = 1,nSpecies
    CALL WriteArrayToHDF5(DataSetName=H5_Name             , rank=4                                           , &
                            nValGlobal =(/nVar2D_Total      , nSurfSample , nSurfSample , nGlobalSides   /)  , &
                            nVal       =(/nVar2D_Spec       , nSurfSample , nSurfSample , nLocalSides/)      , &
                            offset     =(/INT(nVarCount,IK) , 0_IK        , 0_IK        , offsetSurfSide/)   , &
                            collective =.FALSE.                                                              , &
                            RealArray  = MacroSurfaceSpecVal(1:nVar2D_Spec,1:nSurfSample,1:nSurfSample,1:nLocalSides,iSpec))
    nVarCount = nVarCount + INT(nVar2D_Spec)
  END DO
  CALL WriteArrayToHDF5(DataSetName=H5_Name            , rank=4                                              , &
                        nValGlobal =(/nVar2D_Total     , nSurfSample, nSurfSample , nGlobalSides/)           , &
                        nVal       =(/nVar2D           , nSurfSample, nSurfSample , nLocalSides/)            , &
                        offset     =(/INT(nVarCount,IK), 0_IK       , 0_IK        , offsetSurfSide/)         , &
                        collective =.FALSE.                                                                  , &
                        RealArray  = MacroSurfaceVal(1:nVar2D,1:nSurfSample,1:nSurfSample,1:nLocalSides))
END ASSOCIATE

CALL CloseDataFile()

IF (mySurfRank.EQ.0) THEN
  tend=LOCALTIME()
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',tend-tstart,'s]'
END IF

END SUBROUTINE WriteSurfSampleToHDF5


SUBROUTINE WriteSurfSampleChemToHDF5(MeshFileName,OutputTime)
!===================================================================================================================================
!> write the final values of the surface sampling to a HDF5 state file
!> additional performs all the final required computations
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: ProjectName
USE MOD_DSMC_Vars               ,ONLY: CollisMode
USE MOD_HDF5_Output             ,ONLY: WriteAttributeToHDF5,WriteArrayToHDF5,WriteHDF5Header
USE MOD_IO_HDF5
USE MOD_MPI_Shared_Vars         ,ONLY: mySurfRank
USE MOD_SurfaceModel_Vars       ,ONLY: ChemWallProp_Shared_Win, ChemWallProp
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfSample,CalcSurfaceImpact
USE MOD_Particle_Boundary_Vars  ,ONLY: nOutputSides, nComputeNodeSurfSides
USE MOD_Particle_boundary_Vars  ,ONLY: nComputeNodeSurfOutputSides,offsetComputeNodeSurfOutputSide
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfBC,SurfBCName, PartBound
USE MOD_Particle_Vars           ,ONLY: nSpecies
#if USE_MPI
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfTotalSides
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_LEADERS_SURF, MPI_COMM_SHARED
USE MOD_MPI_Shared
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: MeshFileName
REAL,INTENT(IN)                      :: OutputTime
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)                  :: FileName,FileString,Statedummy
CHARACTER(LEN=255)                  :: H5_Name
CHARACTER(LEN=255)                  :: NodeTypeTemp
CHARACTER(LEN=255)                  :: SpecID, PBCID
CHARACTER(LEN=255),ALLOCATABLE      :: Str2DVarNames(:)
INTEGER                             :: nVar2D, nVar2D_Spec, nVar2D_Total, nVarCount, iSpec, iPBC, iSurfSide, nVar2D_Heat
INTEGER                             :: p,q,OutputCounter
REAL                                :: tstart,tend
REAL, ALLOCATABLE                   :: MacroSurfaceSpecChemVal(:,:,:,:,:)
REAL, ALLOCATABLE                   :: MacroSurfaceHeatVal(:,:,:,:)
!===================================================================================================================================

#if USE_MPI
! Return if not a sampling leader
IF (MPI_COMM_LEADERS_SURF.EQ.MPI_COMM_NULL) RETURN
CALL MPI_BARRIER(MPI_COMM_LEADERS_SURF,iERROR)

! Return if no sampling sides
IF (nSurfTotalSides      .EQ.0) RETURN
#endif /*USE_MPI*/

IF (mySurfRank.EQ.0) THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE DSMCSurfChemSTATE TO HDF5 FILE...'
  tstart=LOCALTIME()
END IF

FileName   = TIMESTAMP(TRIM(ProjectName)//'_DSMCSurfChemState',OutputTime)
FileString = TRIM(FileName)//'.h5'

! Create dataset attribute "SurfVarNames"
nVar2D      = 0
nVar2D_Spec = 1
nVar2D_Heat = 1


nVar2D_Total = nVar2D + nVar2D_Spec*nSpecies + nVar2D_Heat

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
#if USE_MPI
IF (mySurfRank.EQ.0) THEN
#endif
  CALL OpenDataFile(FileString,create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
  Statedummy = 'DSMCSurfChemState'

  ! Write file header
  CALL WriteHDF5Header(Statedummy,File_ID)
  CALL WriteAttributeToHDF5(File_ID,'DSMC_nSurfSample',1,IntegerScalar=nSurfSample)
  CALL WriteAttributeToHDF5(File_ID,'DSMC_nSpecies'   ,1,IntegerScalar=nSpecies)
  CALL WriteAttributeToHDF5(File_ID,'DSMC_CollisMode' ,1,IntegerScalar=CollisMode)
  CALL WriteAttributeToHDF5(File_ID,'MeshFile'        ,1,StrScalar=(/TRIM(MeshFileName)/))
  CALL WriteAttributeToHDF5(File_ID,'Time'            ,1,RealScalar=OutputTime)
  CALL WriteAttributeToHDF5(File_ID,'BC_Surf'         ,nSurfBC,StrArray=SurfBCName)
  CALL WriteAttributeToHDF5(File_ID,'N',1             ,IntegerScalar=nSurfSample)
  NodeTypeTemp='VISU'
  CALL WriteAttributeToHDF5(File_ID,'NodeType'        ,1,StrScalar=(/NodeTypeTemp/))

  ALLOCATE(Str2DVarNames(1:nVar2D_Total))
  Str2DVarNames(:) = ''
  nVarCount        = 1
  DO iSpec = 1,nSpecies
    WRITE(SpecID,'(I3.3)') iSpec
    CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Spec'//TRIM(SpecID)//'_Coverage')
  END DO ! iSpec=1,nSpecies

!  ! fill varnames for total values
!  CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Total_Coverage')

  CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Heat_Flux')

  CALL WriteAttributeToHDF5(File_ID,'VarNamesSurface',nVar2D_Total,StrArray=Str2DVarNames)
  CALL CloseDataFile()
  DEALLOCATE(Str2DVarNames)
#if USE_MPI
END IF

CALL MPI_BARRIER(MPI_COMM_LEADERS_SURF,iERROR)

CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_LEADERS_SURF)
#else
CALL OpenDataFile(FileString,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
#endif

ALLOCATE(MacroSurfaceSpecChemVal(1:nVar2D_Spec , 1:nSurfSample , 1:nSurfSample , nComputeNodeSurfOutputSides , nSpecies))
MacroSurfaceSpecChemVal = 0.

ALLOCATE(MacroSurfaceHeatVal(1:nVar2D_Heat , 1:nSurfSample , 1:nSurfSample , nComputeNodeSurfOutputSides))
MacroSurfaceHeatVal = 0.


OutputCounter = 0
DO iSurfSide = 1,nComputeNodeSurfSides
  !================== INNER BC CHECK TODO !
!  GlobalSideID = SurfSide2GlobalSide(SURF_SIDEID,iSurfSide)
!  IF(SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID).GT.0) THEN
!    IF(GlobalSideID.LT.SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID)) THEN
!      SurfSideNb = GlobalSide2SurfSide(SURF_SIDEID,SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID))
!      SampWallState(:,:,:,iSurfSide) = SampWallState(:,:,:,iSurfSide) + SampWallState(:,:,:,SurfSideNb)
!    ELSE
!      CYCLE
!    END IF
!  END IF
  !================== INNER BC CHECK
  OutputCounter = OutputCounter + 1
  DO q = 1,nSurfSample
    DO p = 1,nSurfSample
      ! --- Total output (force per area, heat flux, simulation particle impact per iteration)
      ! --- Species-specific output
      MacroSurfaceHeatVal(1,p,q,OutputCounter) = SUM(ChemWallProp(:,2,p, q, iSurfSide))
      DO iSpec=1,nSpecies
        ! Species-specific counter of simulation particle impacts per iteration
        MacroSurfaceSpecChemVal(1,p,q,OutputCounter,iSpec) = ChemWallProp(iSpec,1,p, q, iSurfSide)
      END DO ! iSpec=1,nSpecies
    END DO ! q=1,nSurfSample
  END DO ! p=1,nSurfSample
END DO ! iSurfSide=1,nComputeNodeSurfSides


nVarCount=0
WRITE(H5_Name,'(A)') 'SurfaceData'
ASSOCIATE (&
      nVar2D_Total         => INT(nVar2D_Total,IK)                    , &
      nSurfSample          => INT(nSurfSample,IK)                     , &
      nGlobalSides         => INT(nOutputSides,IK)                    , &
      nLocalSides          => INT(nComputeNodeSurfOutputSides,IK)     , &
      offsetSurfSide       => INT(offsetComputeNodeSurfOutputSide,IK) , &
      nVar2D_Spec          => INT(nVar2D_Spec,IK)                     , &
      nVar2D               => INT(nVar2D,IK))
  DO iSpec = 1,nSpecies
    CALL WriteArrayToHDF5(DataSetName=H5_Name             , rank=4                                           , &
                            nValGlobal =(/nVar2D_Total      , nSurfSample , nSurfSample , nGlobalSides   /)  , &
                            nVal       =(/nVar2D_Spec       , nSurfSample , nSurfSample , nLocalSides/)      , &
                            offset     =(/INT(nVarCount,IK) , 0_IK        , 0_IK        , offsetSurfSide/)   , &
                            collective =.FALSE.                                                              , &
                            RealArray  = MacroSurfaceSpecChemVal(1:nVar2D_Spec,1:nSurfSample,1:nSurfSample,1:nLocalSides,iSpec))
    nVarCount = nVarCount + INT(nVar2D_Spec)
  END DO
!  CALL WriteArrayToHDF5(DataSetName=H5_Name            , rank=4                                              , &
!                        nValGlobal =(/nVar2D_Total     , nSurfSample, nSurfSample , nGlobalSides/)           , &
!                        nVal       =(/nVar2D           , nSurfSample, nSurfSample , nLocalSides/)            , &
!                        offset     =(/INT(nVarCount,IK), 0_IK       , 0_IK        , offsetSurfSide/)         , &
!                        collective =.FALSE.                                                                  , &
!                        RealArray  = MacroSurfaceVal(1:nVar2D,1:nSurfSample,1:nSurfSample,1:nLocalSides))

  CALL WriteArrayToHDF5(DataSetName=H5_Name            , rank=4                                              , &
                        nValGlobal =(/nVar2D_Total     , nSurfSample, nSurfSample , nGlobalSides/)           , &
                        nVal       =(/nVar2D           , nSurfSample, nSurfSample , nLocalSides/)            , &
                        offset     =(/INT(nVarCount,IK), 0_IK       , 0_IK        , offsetSurfSide/)         , &
                        collective =.FALSE.                                                                  , &
                        RealArray  = MacroSurfaceHeatVal(1:nVar2D_Spec,1:nSurfSample,1:nSurfSample,1:nLocalSides))

END ASSOCIATE

CALL CloseDataFile()

IF (mySurfRank.EQ.0) THEN
  tend=LOCALTIME()
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',tend-tstart,'s]'
END IF

END SUBROUTINE WriteSurfSampleChemToHDF5


SUBROUTINE AddVarName(StrArray,ArrayDim,idx,VarName)
!----------------------------------------------------------------------------------------------------------------------------------!
! description
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: ArrayDim
CHARACTER(LEN=*),INTENT(INOUT) :: StrArray(ArrayDim)
INTEGER,INTENT(INOUT)          :: idx
CHARACTER(LEN=*),INTENT(IN)    :: VarName
! Space-separated list of input and output types. Use: (int|real|logical|...)_(in|out|inout)_dim(n)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
StrArray(idx)=TRIM(VarName)
idx=idx+1
END SUBROUTINE AddVarName


SUBROUTINE FinalizeParticleBoundarySampling()
!===================================================================================================================================
! deallocate everything
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_DSMC_Vars                      ,ONLY: DSMC
USE MOD_Particle_Boundary_Vars
USE MOD_Particle_Vars               ,ONLY: WriteMacroSurfaceValues
USE MOD_SurfaceModel_Vars           ,ONLY: nPorousBC
USE MOD_Particle_Mesh_Vars          ,ONLY: GEO
#if USE_MPI
USE MOD_MPI_Shared_Vars                ,ONLY: MPI_COMM_SHARED,MPI_COMM_LEADERS_SURF
USE MOD_MPI_Shared
USE MOD_Particle_MPI_Boundary_Sampling ,ONLY: FinalizeSurfCommunication
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! Return if nothing was allocated
IF (.NOT.WriteMacroSurfaceValues.AND..NOT.DSMC%CalcSurfaceVal.AND..NOT.(ANY(PartBound%Reactive)).AND..NOT.(nPorousBC.GT.0).AND..NOT.GEO%RotPeriodicBC) RETURN

! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
#if USE_MPI
! Mapping arrays are allocated even if the node does not have sampling surfaces
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
CALL UNLOCK_AND_FREE(GlobalSide2SurfSide_Shared_Win)
CALL UNLOCK_AND_FREE(SurfSide2GlobalSide_Shared_Win)
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
ADEALLOCATE(GlobalSide2SurfSide_Shared)
ADEALLOCATE(SurfSide2GlobalSide_Shared)
#endif /*USE_MPI*/

! Return if no sampling surfaces on node
IF (.NOT.SurfOnNode) RETURN

! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

CALL UNLOCK_AND_FREE(SampWallState_Shared_Win)
CALL UNLOCK_AND_FREE(SurfSideArea_Shared_Win)
IF(nPorousBC.GT.0) CALL UNLOCK_AND_FREE(SampWallPumpCapacity_Shared_Win)
IF (CalcSurfaceImpact) THEN
  CALL UNLOCK_AND_FREE(SampWallImpactEnergy_Shared_Win)
  CALL UNLOCK_AND_FREE(SampWallImpactVector_Shared_Win)
  CALL UNLOCK_AND_FREE(SampWallImpactAngle_Shared_Win)
  CALL UNLOCK_AND_FREE(SampWallImpactNumber_Shared_Win)
END IF

CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

! Communication is handled in particle_mpi_boundary_sampling.f90
CALL FinalizeSurfCommunication()

IF(MPI_COMM_LEADERS_SURF.NE.MPI_COMM_NULL) THEN
  CALL MPI_COMM_FREE(MPI_COMM_LEADERS_SURF,iERROR)
END IF

ADEALLOCATE(SampWallState_Shared)
ADEALLOCATE(SampWallPumpCapacity_Shared)
ADEALLOCATE(SampWallImpactEnergy_Shared)
ADEALLOCATE(SampWallImpactVector_Shared)
ADEALLOCATE(SampWallImpactAngle_Shared)
ADEALLOCATE(SampWallImpactNumber_Shared)
ADEALLOCATE(SurfSideArea_Shared)
#endif /*USE_MPI*/

! Then, free the pointers or arrays
SDEALLOCATE(XiEQ_SurfSample)
SDEALLOCATE(SurfBCName)
SDEALLOCATE(SampWallState)
SDEALLOCATE(SampWallPumpCapacity)
SDEALLOCATE(SampWallImpactEnergy)
SDEALLOCATE(SampWallImpactVector)
SDEALLOCATE(SampWallImpactAngle)
SDEALLOCATE(SampWallImpactNumber)
ADEALLOCATE(SurfSideArea)
ADEALLOCATE(GlobalSide2SurfSide)
ADEALLOCATE(SurfSide2GlobalSide)

END SUBROUTINE FinalizeParticleBoundarySampling

END MODULE MOD_Particle_Boundary_Sampling
