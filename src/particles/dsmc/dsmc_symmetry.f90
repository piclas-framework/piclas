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

MODULE MOD_DSMC_Symmetry
!===================================================================================================================================
!> Routines for 2D (planar/axisymmetric) simulations
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
PUBLIC :: DSMC_1D_InitVolumes, DSMC_2D_InitVolumes, DSMC_2D_InitRadialWeighting, DSMC_2D_RadialWeighting, DSMC_2D_SetInClones
PUBLIC :: DSMC_2D_CalcSymmetryArea, DSMC_1D_CalcSymmetryArea, DSMC_2D_CalcSymmetryAreaSubSides
PUBLIC :: DefineParametersParticleSymmetry, Init_Symmetry, DSMC_2D_TreatIdenticalParticles
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for particles
!==================================================================================================================================
SUBROUTINE DefineParametersParticleSymmetry()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE

CALL prms%SetSection("Particle Symmetry")
CALL prms%CreateIntOption(    'Particles-Symmetry-Order',  &
                              'Order of the Simulation 1, 2 or 3 D', '3')
CALL prms%CreateLogicalOption('Particles-Symmetry2D', 'Activating a 2D simulation on a mesh with one cell in z-direction in the '//&
                              'xy-plane (y ranging from 0 to the domain boundaries)', '.FALSE.')
CALL prms%CreateLogicalOption('Particles-Symmetry2DAxisymmetric', 'Activating an axisymmetric simulation with the same mesh '//&
                              'requirements as for the 2D case (y is then the radial direction)', '.FALSE.')
CALL prms%CreateLogicalOption('Particles-RadialWeighting', 'Activates a radial weighting in y for the axisymmetric '//&
                              'simulation based on the particle position.', '.FALSE.')
CALL prms%CreateRealOption(   'Particles-RadialWeighting-PartScaleFactor', 'Axisymmetric radial weighting factor, defining '//&
                              'the linear increase of the weighting factor (e.g. factor 2 means that the weighting factor will '//&
                              'be twice as large at the outer radial domain boundary than at the rotational axis')
CALL prms%CreateLogicalOption('Particles-RadialWeighting-CellLocalWeighting', 'Enables a cell-local radial weighting, '//&
                              'where every particle has the same weighting factor within a cell', '.FALSE.')
CALL prms%CreateIntOption(    'Particles-RadialWeighting-CloneMode',  &
                              'Radial weighting: Select between methods for the delayed insertion of cloned particles:/n'//&
                              '1: Chronological, 2: Random', '2')
CALL prms%CreateIntOption(    'Particles-RadialWeighting-CloneDelay', &
                              'Radial weighting:  Delay (number of iterations) before the stored cloned particles are inserted '//&
                              'at the position they were cloned', '2')
CALL prms%CreateIntOption(    'Particles-RadialWeighting-SurfFluxSubSides', &
                              'Radial weighting: Split the surface flux side into the given number of subsides, reduces the '//&
                              'error in the particle distribution across the cell (visible in the number density)', '20')

END SUBROUTINE DefineParametersParticleSymmetry


SUBROUTINE DSMC_2D_InitVolumes()
!===================================================================================================================================
!> Routine determines the symmetry sides and calculates the 2D (area faces in symmetry plane) and axisymmetric volumes (cells are
!> revolved around the symmetry axis). The symmetry side (SymmetrySide array) will be used later on to determine in which two
!> directions the quadtree shall refine the mesh, skipping the z-dimension to avoid an unnecessary refinement. Additionally,
!> symmetry sides will be skipped during tracking (SideIsSymSide_Shared array).
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: Pi
USE MOD_PreProc
USE MOD_Mesh_Vars               ,ONLY: nElems,offsetElem,nBCSides,SideToElem
USE MOD_Particle_Vars           ,ONLY: Symmetry
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO,LocalVolume,MeshVolume
USE MOD_DSMC_Vars               ,ONLY: SymmetrySide
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemVolume_Shared,ElemCharLength_Shared
USE MOD_Particle_Mesh_Vars      ,ONLY: NodeCoords_Shared,ElemSideNodeID_Shared, SideInfo_Shared, SideIsSymSide_Shared, SideIsSymSide
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces       ,ONLY: CalcNormAndTangTriangle
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars      ,ONLY: nNonUniqueGlobalSides, offsetComputeNodeElem, ElemInfo_Shared
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemVolume_Shared_Win,ElemCharLength_Shared_Win,SideIsSymSide_Shared_Win
USE MOD_MPI_Shared_Vars         ,ONLY: myComputeNodeRank,nComputeNodeProcessors,MPI_COMM_LEADERS_SHARED
#else
USE MOD_Particle_Mesh_Vars      ,ONLY: nComputeNodeSides
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: SideID, iLocSide, iNode, BCSideID, locElemID, CNElemID, iSide
REAL                            :: radius, triarea(2)
#if USE_MPI
REAL                            :: CNVolume
INTEGER                         :: offsetElemCNProc
#endif /*USE_MPI*/
LOGICAL                         :: SymmetryBCExists
INTEGER                         :: firstElem, lastElem, firstSide, lastSide
!===================================================================================================================================

#if USE_MPI
CALL Allocate_Shared((/nNonUniqueGlobalSides/),SideIsSymSide_Shared_Win,SideIsSymSide_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideIsSymSide_Shared_Win,IERROR)
SideIsSymSide => SideIsSymSide_Shared
! only CN root nullifies
IF(myComputeNodeRank.EQ.0) SideIsSymSide = .FALSE.
! This sync/barrier is required as it cannot be guaranteed that the zeros have been written to memory by the time the MPI_REDUCE
! is executed (see MPI specification). Until the Sync is complete, the status is undefined, i.e., old or new value or utter nonsense.
CALL BARRIER_AND_SYNC(SideIsSymSide_Shared_Win,MPI_COMM_SHARED)
#else
ALLOCATE(SideIsSymSide(nComputeNodeSides))
SideIsSymSide = .FALSE.
#endif  /*USE_MPI*/

! Flag of symmetry sides to be skipped during tracking
#if USE_MPI
  firstSide = INT(REAL( myComputeNodeRank   *nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))+1
  lastSide  = INT(REAL((myComputeNodeRank+1)*nNonUniqueGlobalSides)/REAL(nComputeNodeProcessors))
#else
  firstSide = 1
  lastSide  = nComputeNodeSides
#endif

DO iSide = firstSide, lastSide
  SideIsSymSide(iSide) = .FALSE.
  ! ignore non-BC sides
  IF (SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE
#if USE_MPI
  ! ignore sides outside of halo region
  IF (ElemInfo_Shared(ELEM_HALOFLAG,SideInfo_Shared(SIDE_ELEMID,iSide)).EQ.0) CYCLE
#endif /*USE_MPI*/
  IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,iSide))).EQ.PartBound%SymmetryBC) &
    SideIsSymSide(iSide) = .TRUE.
END DO

SymmetryBCExists = .FALSE.
ALLOCATE(SymmetrySide(1:nElems,1:2))                ! 1: GlobalSide, 2: LocalSide
SymmetrySide = -1

! Sanity check: mesh has to be centered at z = 0
IF(.NOT.ALMOSTEQUALRELATIVE(GEO%zmaxglob,ABS(GEO%zminglob),1e-5)) THEN
  SWRITE(*,*) 'Maximum dimension in z:', GEO%zmaxglob
  SWRITE(*,*) 'Minimum dimension in z:', GEO%zminglob
  SWRITE(*,*) 'Deviation', (ABS(GEO%zmaxglob)-ABS(GEO%zminglob))/ABS(GEO%zminglob), ' > 1e-5'
  CALL abort(__STAMP__&
    ,'ERROR: Please orient your mesh with one cell in z-direction around 0, |z_min| = z_max !')
END IF

! Calculation of the correct volume and characteristic length
DO BCSideID=1,nBCSides
  locElemID = SideToElem(S2E_ELEM_ID,BCSideID)
  iLocSide = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
  SideID=GetGlobalNonUniqueSideID(offsetElem+locElemID,iLocSide)
  IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))).EQ.PartBound%SymmetryBC) THEN
    CNElemID = GetCNElemID(SideInfo_Shared(SIDE_ELEMID,SideID))
    iLocSide = SideInfo_Shared(SIDE_LOCALID,SideID)
    ! Exclude the symmetry axis (y=0)
    IF(Symmetry%Axisymmetric) THEN
      IF(MAXVAL(NodeCoords_Shared(2,ElemSideNodeID_Shared(:,iLocSide,CNElemID)+1)).LE.0.0) CYCLE
    END IF
    ! The z-plane with the positive z component is chosen
    IF(MINVAL(NodeCoords_Shared(3,ElemSideNodeID_Shared(:,iLocSide,CNElemID)+1)).GT.(GEO%zmaxglob+GEO%zminglob)/2.) THEN
      IF(SymmetrySide(locElemID,1).GT.0) THEN
        CALL abort(__STAMP__&
          ,'ERROR: PICLas could not determine a unique symmetry surface for 2D/axisymmetric calculation!'//&
          ' Please orient your mesh with x as the symmetry axis and positive y as the second/radial direction!')
      END IF
      SymmetrySide(locElemID,1) = BCSideID
      SymmetrySide(locElemID,2) = iLocSide
      ! The volume calculated at this point (final volume for the 2D case) corresponds to the cell face area (z-dimension=1) in
      ! the xy-plane.
      ElemVolume_Shared(CNElemID) = 0.0
      CALL CalcNormAndTangTriangle(area=triarea(1),TriNum=1, SideID=SideID)
      CALL CalcNormAndTangTriangle(area=triarea(2),TriNum=2, SideID=SideID)
      ElemVolume_Shared(CNElemID) = triarea(1) + triarea(2)
      ! Characteristic length is compared to the mean free path as the condition to refine the mesh. For the 2D/axisymmetric case
      ! the third dimension is not considered as particle interaction occurs in the xy-plane, effectively reducing the refinement
      ! requirement.
      ElemCharLength_Shared(CNElemID) = SQRT(ElemVolume_Shared(CNElemID))
      ! Axisymmetric case: The volume is multiplied by the circumference to get the volume of the ring. The cell face in the
      ! xy-plane is rotated around the x-axis. The radius is the middle point of the cell face.
      IF (Symmetry%Axisymmetric) THEN
        radius = 0.
        DO iNode = 1, 4
          radius = radius + NodeCoords_Shared(2,ElemSideNodeID_Shared(iNode,iLocSide,CNElemID)+1)
        END DO
        radius = radius / 4.
        ElemVolume_Shared(CNElemID) = ElemVolume_Shared(CNElemID) * 2. * Pi * radius
      END IF
      SymmetryBCExists = .TRUE.
    END IF      ! Greater z-coord
  END IF
END DO

IF(.NOT.SymmetryBCExists) THEN
  CALL abort(__STAMP__&
    ,'At least one symmetric BC (in the xy-plane) has to be defined for 2D simulations')
END IF

! LocalVolume & MeshVolume: Recalculate the volume of the mesh of a single process and the total mesh volume
#if USE_MPI
FirstElem = offsetElem - offsetComputeNodeElem + 1
LastElem  = offsetElem - offsetComputeNodeElem + nElems
#else
firstElem = 1
lastElem  = nElems
#endif  /*USE_MPI*/

LocalVolume = SUM(ElemVolume_Shared(FirstElem:LastElem))

#if USE_MPI
CALL BARRIER_AND_SYNC(SideIsSymSide_Shared_Win ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(ElemVolume_Shared_Win    ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(ElemCharLength_Shared_Win,MPI_COMM_SHARED)
! Compute-node mesh volume
offsetElemCNProc = offsetElem - offsetComputeNodeElem
CNVolume = SUM(ElemVolume_Shared(offsetElemCNProc+1:offsetElemCNProc+nElems))
CALL MPI_ALLREDUCE(MPI_IN_PLACE,CNVolume,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_SHARED,iError)
IF (myComputeNodeRank.EQ.0) THEN
  ! All-reduce between node leaders
  CALL MPI_ALLREDUCE(CNVolume,MeshVolume,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEADERS_SHARED,IERROR)
END IF
! Broadcast from node leaders to other processors on the same node
CALL MPI_BCAST(MeshVolume,1, MPI_DOUBLE_PRECISION,0,MPI_COMM_SHARED,iERROR)
#else
MeshVolume = LocalVolume
#endif /*USE_MPI*/

END SUBROUTINE DSMC_2D_InitVolumes


SUBROUTINE DSMC_1D_InitVolumes()
!===================================================================================================================================
!> Routine determines a symmetry side and calculates the 1D (area faces at x axis) and axisymmetric volumes (cells are
!> revolved around the symmetry axis). The symmetry side will be used later on to determine in which two directions the quadtree
!> shall refine the mesh, skipping the z-dimension to avoid an unnecessary refinement.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars               ,ONLY: nElems, offsetElem
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO,LocalVolume,MeshVolume
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemVolume_Shared,ElemCharLength_Shared
USE MOD_Particle_Mesh_Vars      ,ONLY: NodeCoords_Shared,ElemSideNodeID_Shared, SideInfo_Shared
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalNonUniqueSideID
#if USE_MPI
USE MOD_MPI_Shared              ,ONLY: BARRIER_AND_SYNC
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED
USE MOD_Particle_Mesh_Vars      ,ONLY: offsetComputeNodeElem
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemVolume_Shared_Win,ElemCharLength_Shared_Win
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iLocSide, iElem, SideID, iDim, CNElemID
REAL                            :: X(2), MaxCoord, MinCoord
LOGICAL                         :: SideInPlane, X1Occupied
INTEGER                         :: firstElem,lastElem
!===================================================================================================================================

#if USE_MPI
firstElem = offsetElem - offsetComputeNodeElem + 1
lastElem  = offsetElem - offsetComputeNodeElem + nElems
#else
firstElem = 1
lastElem  = nElems
#endif  /*USE_MPI*/

ALLOCATE(GEO%XMinMax(2,nElems))

IF(.NOT.ALMOSTEQUALRELATIVE(GEO%ymaxglob,ABS(GEO%yminglob),1e-5)) THEN
  SWRITE(*,*) 'Maximum dimension in y:', GEO%ymaxglob
  SWRITE(*,*) 'Minimum dimension in y:', GEO%yminglob
  SWRITE(*,*) 'Deviation', (ABS(GEO%ymaxglob)-ABS(GEO%yminglob))/ABS(GEO%yminglob), ' > 1e-5'
  CALL abort(__STAMP__&
    ,'ERROR: Please orient your mesh with one cell in y-direction around 0, |z_min| = z_max !')
END IF

IF(.NOT.ALMOSTEQUALRELATIVE(GEO%zmaxglob,ABS(GEO%zminglob),1e-5)) THEN
  SWRITE(*,*) 'Maximum dimension in z:', GEO%zmaxglob
  SWRITE(*,*) 'Minimum dimension in z:', GEO%zminglob
  SWRITE(*,*) 'Deviation', (ABS(GEO%zmaxglob)-ABS(GEO%zminglob))/ABS(GEO%zminglob), ' > 1e-5'
  CALL abort(__STAMP__&
    ,'ERROR: Please orient your mesh with one cell in z-direction around 0, |z_min| = z_max !')
END IF

DO iElem = 1,nElems
  ! Check if all sides of the element are parallel to xy-, xz-, or yz-plane and Sides parallel to xy-,and xz-plane are symmetric
  ! And determine xmin and xmax of the element
  X = 0
  X1Occupied = .FALSE.
  DO iLocSide = 1,6
    SideInPlane=.FALSE.
    SideID = GetGlobalNonUniqueSideID(offsetElem+iElem,iLocSide)
    CNElemID = GetCNElemID(SideInfo_Shared(SIDE_ELEMID,SideID))
    DO iDim=1,3
      MaxCoord = MAXVAL(NodeCoords_Shared(iDim,ElemSideNodeID_Shared(:,iLocSide,CNElemID)+1))
      MinCoord = MINVAL(NodeCoords_Shared(iDim,ElemSideNodeID_Shared(:,iLocSide,CNElemID)+1))
      IF(ALMOSTALMOSTEQUAL(MaxCoord,MinCoord).OR.(ALMOSTZERO(MaxCoord).AND.ALMOSTZERO(MinCoord))) THEN
        IF(SideInPlane) CALL abort(__STAMP__&
          ,'ERROR: Please orient your mesh with all element sides parallel to xy-,xz-,or yz-plane')
        SideInPlane=.TRUE.
        IF(iDim.GE.2)THEN
          IF(PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))).NE.PartBound%SymmetryBC) &
            CALL abort(__STAMP__,&
              'ERROR: Sides parallel to xy-,and xz-plane has to be the symmetric boundary condition')
        END IF
        IF(iDim.EQ.1) THEN
          IF(X1Occupied) THEN
            X(2) = NodeCoords_Shared(1,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)
          ELSE
            X(1) = NodeCoords_Shared(1,ElemSideNodeID_Shared(1,iLocSide,CNElemID)+1)
            X1Occupied = .TRUE.
          END IF
        END IF !iDim.EQ.1
      END IF
    END DO !iDim=1,3
    IF(.NOT.SideInPlane) THEN
      IPWRITE(*,*) 'ElemID:',iElem,'SideID:',iLocSide
      DO iDim=1,4
        IPWRITE(*,*) 'Node',iDim,'x:',NodeCoords_Shared(1,ElemSideNodeID_Shared(iDim,iLocSide,CNElemID)+1), &
                                 'y:',NodeCoords_Shared(2,ElemSideNodeID_Shared(iDim,iLocSide,CNElemID)+1), &
                                 'z:',NodeCoords_Shared(3,ElemSideNodeID_Shared(iDim,iLocSide,CNElemID)+1)
      END DO
      CALL abort(__STAMP__&
    ,'ERROR: Please orient your mesh with all element sides parallel to xy-,xz-,or yz-plane')
    END IF
  END DO ! iLocSide = 1,6
  GEO%XMinMax(1,iElem) = MINVAL(X)
  GEO%XMinMax(2,iElem) = MAXVAL(X)
  ElemVolume_Shared(CNElemID) = ABS(X(1)-X(2))
  ElemCharLength_Shared(CNElemID) = ElemVolume_Shared(CNElemID)
END DO
! LocalVolume & MeshVolume: Recalculate the volume of the mesh of a single process and the total mesh volume
LocalVolume = SUM(ElemVolume_Shared(firstElem:lastElem))
#if USE_MPI
CALL BARRIER_AND_SYNC(ElemVolume_Shared_Win    ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(ElemCharLength_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/
MeshVolume = SUM(ElemVolume_Shared)

END SUBROUTINE DSMC_1D_InitVolumes


SUBROUTINE DSMC_2D_InitRadialWeighting()
!===================================================================================================================================
!> Read-in and initialize the variables required for the cloning procedures. Two modes with a delayed clone insertion are available:
!> 1: Insert the clones after the delay in the same chronological order as they were created
!> 2: Choose a random list of particles to insert after the delay buffer is full with clones
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Restart_Vars            ,ONLY: DoRestart
USE MOD_DSMC_Vars               ,ONLY: RadialWeighting, ClonedParticles
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars        ,ONLY: DoLoadBalance, UseH5IOLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! Clone read-in during load balance is currently only supported via the HDF5 output
#if USE_LOADBALANCE
IF(DoLoadBalance.AND.(.NOT.UseH5IOLoadBalance)) THEN
  CALL abort(__STAMP__,'ERROR: Radial weighting only supports a load balance using an HDF5 output (UseH5IOLoadBalance = T)!')
END IF
#endif /*USE_LOADBALANCE*/

! Linear increasing weighting factor in the radial direction up to the domain boundary
RadialWeighting%PartScaleFactor = GETREAL('Particles-RadialWeighting-PartScaleFactor')
IF(RadialWeighting%PartScaleFactor.LT.1.) THEN
  CALL Abort(__STAMP__,'ERROR in 2D axisymmetric simulation: PartScaleFactor has to be greater than 1!',RealInfoOpt=RadialWeighting%PartScaleFactor)
END IF
RadialWeighting%CloneMode = GETINT('Particles-RadialWeighting-CloneMode')
RadialWeighting%CloneInputDelay = GETINT('Particles-RadialWeighting-CloneDelay')
! Cell local radial weighting (all particles have the same weighting factor within a cell)
RadialWeighting%CellLocalWeighting = GETLOGICAL('Particles-RadialWeighting-CellLocalWeighting')

! Number of subsides to split the surface flux sides into, otherwise a wrong distribution of particles across large cells will be
! inserted, visible in the number density as an increase in the number density closer the axis (e.g. resulting in a heat flux peak)
! (especially when using mortar meshes)
RadialWeighting%nSubSides=GETINT('Particles-RadialWeighting-SurfFluxSubSides')

RadialWeighting%NextClone = 0
RadialWeighting%CloneVecLengthDelta = 100
RadialWeighting%CloneVecLength = RadialWeighting%CloneVecLengthDelta

SELECT CASE(RadialWeighting%CloneMode)
  CASE(1)
    IF(RadialWeighting%CloneInputDelay.LT.1) THEN
      CALL Abort(__STAMP__,'ERROR in 2D axisymmetric simulation: Clone delay should be greater than 0')
    END IF
    ALLOCATE(RadialWeighting%ClonePartNum(0:(RadialWeighting%CloneInputDelay-1)))
    ALLOCATE(ClonedParticles(1:RadialWeighting%CloneVecLength,0:(RadialWeighting%CloneInputDelay-1)))
    RadialWeighting%ClonePartNum = 0
    IF(.NOT.DoRestart) RadialWeighting%CloneDelayDiff = 1
  CASE(2)
    IF(RadialWeighting%CloneInputDelay.LT.2) THEN
      CALL Abort(__STAMP__,'ERROR in 2D axisymmetric simulation: Clone delay should be greater than 1')
    END IF
    ALLOCATE(RadialWeighting%ClonePartNum(0:RadialWeighting%CloneInputDelay))
    ALLOCATE(ClonedParticles(1:RadialWeighting%CloneVecLength,0:RadialWeighting%CloneInputDelay))
    RadialWeighting%ClonePartNum = 0
    IF(.NOT.DoRestart) RadialWeighting%CloneDelayDiff = 0
  CASE DEFAULT
    CALL Abort(__STAMP__,'ERROR in Radial Weighting of 2D/Axisymmetric: The selected cloning mode is not available! Choose between 1 and 2.'//&
        ' CloneMode=1: Delayed insertion of clones; CloneMode=2: Delayed randomized insertion of clones')
END SELECT

END SUBROUTINE DSMC_2D_InitRadialWeighting


SUBROUTINE DSMC_2D_RadialWeighting(iPart,iElem)
!===================================================================================================================================
!> Routine for the treatment of particles with enabled radial weighting (weighting factor is increasing linearly with increasing y)
!> 1.) Determine the new particle weight and decide whether to clone or to delete the particle
!> 2a.) Particle cloning, if the local weighting factor is smaller than the previous (particle travelling downwards)
!> 2b.) Particle deletion, if the local weighting factor is greater than the previous (particle travelling upwards)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars               ,ONLY: RadialWeighting, DSMC, PartStateIntEn, useDSMC, CollisMode, AmbipolElecVelo
USE MOD_DSMC_Vars               ,ONLY: ClonedParticles, VibQuantsPar, SpecDSMC, PolyatomMolDSMC, ElectronicDistriPart
USE MOD_Particle_Vars           ,ONLY: PartMPF, PartSpecies, PartState, Species, LastPartPos
USE MOD_TimeDisc_Vars           ,ONLY: iter
USE MOD_part_operations         ,ONLY: RemoveParticle
USE MOD_part_tools              ,ONLY: CalcRadWeightMPF
USE Ziggurat
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: iPart, iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: SpecID, iPolyatMole, cloneIndex, DelayCounter
REAL                            :: DeleteProb, iRan, NewMPF, CloneProb, OldMPF
LOGICAL                         :: DoCloning
!===================================================================================================================================
DoCloning = .FALSE.
DeleteProb = 0.
SpecID = PartSpecies(iPart)

IF (.NOT.(PartMPF(iPart).GT.Species(SpecID)%MacroParticleFactor)) RETURN

! 1.) Determine the new particle weight and decide whether to clone or to delete the particle
NewMPF = CalcRadWeightMPF(PartState(2,iPart),SpecID,iPart)
OldMPF = PartMPF(iPart)
CloneProb = (OldMPF/NewMPF)-INT(OldMPF/NewMPF)
CALL RANDOM_NUMBER(iRan)
IF((CloneProb.GT.iRan).AND.(NewMPF.LT.OldMPF)) THEN
  DoCloning = .TRUE.
  IF(INT(OldMPF/NewMPF).GT.1) THEN
    IPWRITE(*,*) 'New weighting factor:', NewMPF, 'Old weighting factor:', OldMPF
    CALL Abort(&
        __STAMP__,&
      'ERROR in 2D axisymmetric simulation: More than one clone per particle is not allowed! Reduce the time step or'//&
        ' the radial weighting factor! Cloning probability is:',RealInfoOpt=CloneProb)
  END IF
END IF
PartMPF(iPart) = NewMPF

IF(DoCloning) THEN
  ! 2a.) Particle cloning, if the local weighting factor is smaller than the previous (particle travelling downwards)
  ! Get the list number to store the clones, depending on the chosen clone mode
  SELECT CASE(RadialWeighting%CloneMode)
  CASE(1)
  ! ######## Clone Delay ###################################################################################################
  ! Insertion of the clones after a defined delay, all clones are collected in a single list and inserted before boundary
  ! treatment in the next time step at their original positions
    DelayCounter = MOD((INT(iter,4)+RadialWeighting%CloneDelayDiff-1),RadialWeighting%CloneInputDelay)
  CASE(2)
  ! ######## Clone Random Delay #############################################################################################
  ! A list, which is RadialWeighting%CloneInputDelay + 1 long, is filled with clones to be inserted. After the list
  ! is full, NextClone gives the empty particle list, whose clones were inserted during the last SetInClones step
    IF((INT(iter,4)+RadialWeighting%CloneDelayDiff).LE.RadialWeighting%CloneInputDelay) THEN
      DelayCounter = INT(iter,4)+RadialWeighting%CloneDelayDiff
    ELSE
      DelayCounter = RadialWeighting%NextClone
    END IF
  END SELECT
  ! Storing the particle information
  IF(RadialWeighting%ClonePartNum(DelayCounter)+1.GT.RadialWeighting%CloneVecLength) CALL IncreaseClonedParticlesType()
  RadialWeighting%ClonePartNum(DelayCounter) = RadialWeighting%ClonePartNum(DelayCounter) + 1
  cloneIndex = RadialWeighting%ClonePartNum(DelayCounter)
  ClonedParticles(cloneIndex,DelayCounter)%PartState(1:6)= PartState(1:6,iPart)
  IF (useDSMC.AND.(CollisMode.GT.1)) THEN
    ClonedParticles(cloneIndex,DelayCounter)%PartStateIntEn(1:2) = PartStateIntEn(1:2,iPart)
    IF(DSMC%ElectronicModel.GT.0) THEN
      ClonedParticles(cloneIndex,DelayCounter)%PartStateIntEn(3) =   PartStateIntEn(3,iPart)
      IF ((DSMC%ElectronicModel.EQ.2).AND.(.NOT.((SpecDSMC(SpecID)%InterID.EQ.4).OR.SpecDSMC(SpecID)%FullyIonized))) THEN
        IF(ALLOCATED(ClonedParticles(cloneIndex,DelayCounter)%DistriFunc)) &
          DEALLOCATE(ClonedParticles(cloneIndex,DelayCounter)%DistriFunc)
        ALLOCATE(ClonedParticles(cloneIndex,DelayCounter)%DistriFunc(1:SpecDSMC(SpecID)%MaxElecQuant))
        ClonedParticles(cloneIndex,DelayCounter)%DistriFunc(:) = ElectronicDistriPart(iPart)%DistriFunc(:)
      END IF
    END IF
    IF ((DSMC%DoAmbipolarDiff).AND.(Species(SpecID)%ChargeIC.GT.0.0)) THEN
      IF(ALLOCATED(ClonedParticles(cloneIndex,DelayCounter)%AmbiPolVelo)) &
        DEALLOCATE(ClonedParticles(cloneIndex,DelayCounter)%AmbiPolVelo)
      ALLOCATE(ClonedParticles(cloneIndex,DelayCounter)%AmbiPolVelo(1:3))
      ClonedParticles(cloneIndex,DelayCounter)%AmbiPolVelo(1:3) = AmbipolElecVelo(iPart)%ElecVelo(1:3)
    END IF
    IF(SpecDSMC(SpecID)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(SpecID)%SpecToPolyArray
      IF(ALLOCATED(ClonedParticles(cloneIndex,DelayCounter)%VibQuants)) &
        DEALLOCATE(ClonedParticles(cloneIndex,DelayCounter)%VibQuants)
      ALLOCATE(ClonedParticles(cloneIndex,DelayCounter)%VibQuants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF))
      ClonedParticles(cloneIndex,DelayCounter)%VibQuants(:) = VibQuantsPar(iPart)%Quants(:)
    END IF
  END IF
  ClonedParticles(cloneIndex,DelayCounter)%Species = SpecID
  ClonedParticles(cloneIndex,DelayCounter)%Element = iElem
  ClonedParticles(cloneIndex,DelayCounter)%LastPartPos(1:3) = LastPartPos(1:3,iPart)
  ClonedParticles(cloneIndex,DelayCounter)%WeightingFactor = PartMPF(iPart)
ELSE
! ######## Particle Delete #######################################################################################################
! 2b.) Particle deletion, if the local weighting factor is greater than the previous (particle travelling upwards)
  IF(NewMPF.GT.OldMPF) THEN
    ! Start deleting particles after the clone delay has passed and particles are also inserted
    IF((INT(iter,4)+RadialWeighting%CloneDelayDiff).LE.RadialWeighting%CloneInputDelay) RETURN
    DeleteProb = 1. - CloneProb
    IF (DeleteProb.GT.0.5) THEN
      IPWRITE(*,*) 'New weighting factor:', NewMPF, 'Old weighting factor:', OldMPF
      CALL abort(__STAMP__,&
        'ERROR in Radial Weighting of 2D/Axisymmetric: The deletion probability is higher than 0.5! Reduce the time step or'//&
        ' the radial weighting factor! Deletion probability is:',RealInfoOpt=DeleteProb)
    END IF
    CALL RANDOM_NUMBER(iRan)
    IF(DeleteProb.GT.iRan) THEN
      CALL RemoveParticle(iPart)
    END IF
  END IF
END IF

END SUBROUTINE DSMC_2D_RadialWeighting


SUBROUTINE DSMC_2D_SetInClones()
!===================================================================================================================================
!> Insertion of cloned particles during the previous time steps. Clones insertion is delayed by at least one time step to avoid the
!> avalanche phenomenon (identical particles travelling on the same path, not colliding due to zero relative velocity).
!> 1.) Chose which list to insert depending on the clone mode
!> 2.) Insert the clones at the position they were created
!> 3.) Reset the list
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars               ,ONLY: ClonedParticles, PartStateIntEn, useDSMC, CollisMode, DSMC, RadialWeighting
USE MOD_DSMC_Vars               ,ONLY: AmbipolElecVelo
USE MOD_DSMC_Vars               ,ONLY: VibQuantsPar, SpecDSMC, PolyatomMolDSMC, SamplingActive, ElectronicDistriPart
USE MOD_Particle_Vars           ,ONLY: PDM, PEM, PartSpecies, PartState, LastPartPos, PartMPF, WriteMacroVolumeValues, Species
USE MOD_Particle_Vars           ,ONLY: UseVarTimeStep, PartTimeStep
USE MOD_Particle_TimeStep       ,ONLY: GetParticleTimeStep
USE MOD_TimeDisc_Vars           ,ONLY: iter
USE MOD_Particle_Analyze_Vars   ,ONLY: CalcPartBalance, nPartIn
USE MOD_Part_Tools              ,ONLY: GetNextFreePosition
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iPart, PositionNbr, iPolyatMole, DelayCounter, locElemID
REAL                            :: iRan
!===================================================================================================================================

! 1.) Chose which list to insert depending on the clone mode
SELECT CASE(RadialWeighting%CloneMode)
CASE(1)
  ! During the first iterations the delay counter refers to the empty clone array (which is filled during the following tracking)
  ! Afterwards, the MODULUS counts up from zero to CloneInputDelay-1
  DelayCounter = MOD((INT(iter,4)+RadialWeighting%CloneDelayDiff-1),RadialWeighting%CloneInputDelay)
CASE(2)
  ! During the first iterations, check if number of iterations is less than the input delay and leave routine. Afterwards, a
  ! random clone list from the previous time steps is chosen.
  IF((INT(iter,4)+RadialWeighting%CloneDelayDiff).GT.RadialWeighting%CloneInputDelay) THEN
    CALL RANDOM_NUMBER(iRan)
    ! Choosing random clone between 0 and CloneInputDelay
    DelayCounter = INT((RadialWeighting%CloneInputDelay+1)*iRan)
    DO WHILE (DelayCounter.EQ.RadialWeighting%NextClone)
      CALL RANDOM_NUMBER(iRan)
      DelayCounter = INT((RadialWeighting%CloneInputDelay+1)*iRan)
    END DO
    ! Save the chosen list as the next available list to store clones in the next time step
    RadialWeighting%NextClone = DelayCounter
  ELSE
    RETURN
  END IF
END SELECT

IF(RadialWeighting%ClonePartNum(DelayCounter).EQ.0) RETURN

! 2.) Insert the clones at the position they were created
DO iPart = 1, RadialWeighting%ClonePartNum(DelayCounter)
  PositionNbr = GetNextFreePosition()
  ! Copy particle parameters
  PDM%ParticleInside(PositionNbr) = .TRUE.
  PDM%IsNewPart(PositionNbr) = .TRUE.
  PDM%dtFracPush(PositionNbr) = .FALSE.
  PartState(1:5,PositionNbr) = ClonedParticles(iPart,DelayCounter)%PartState(1:5)
  ! Creating a relative velocity in the z-direction
  PartState(6,PositionNbr) = - ClonedParticles(iPart,DelayCounter)%PartState(6)
  IF (useDSMC.AND.(CollisMode.GT.1)) THEN
    PartStateIntEn(1:2,PositionNbr) = ClonedParticles(iPart,DelayCounter)%PartStateIntEn(1:2)
    IF(DSMC%ElectronicModel.GT.0) THEN
      PartStateIntEn(3,PositionNbr) = ClonedParticles(iPart,DelayCounter)%PartStateIntEn(3)
      IF ((DSMC%ElectronicModel.EQ.2).AND.(.NOT.((SpecDSMC(ClonedParticles(iPart,DelayCounter)%Species)%InterID.EQ.4) &
          .OR.SpecDSMC(ClonedParticles(iPart,DelayCounter)%Species)%FullyIonized))) THEN
        IF(ALLOCATED(ElectronicDistriPart(PositionNbr)%DistriFunc)) DEALLOCATE(ElectronicDistriPart(PositionNbr)%DistriFunc)
        ALLOCATE(ElectronicDistriPart(PositionNbr)%DistriFunc(1:SpecDSMC(ClonedParticles(iPart,DelayCounter)%Species)%MaxElecQuant))
        ElectronicDistriPart(PositionNbr)%DistriFunc(:) = ClonedParticles(iPart,DelayCounter)%DistriFunc(:)
      END IF
    END IF
    IF ((DSMC%DoAmbipolarDiff).AND.(Species(ClonedParticles(iPart,DelayCounter)%Species)%ChargeIC.GT.0.0)) THEN
      IF(ALLOCATED(AmbipolElecVelo(PositionNbr)%ElecVelo)) DEALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo)
      ALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo(1:3))
      AmbipolElecVelo(PositionNbr)%ElecVelo(1:2) = ClonedParticles(iPart,DelayCounter)%AmbiPolVelo(1:2)
      AmbipolElecVelo(PositionNbr)%ElecVelo(3) = -ClonedParticles(iPart,DelayCounter)%AmbiPolVelo(3)
    END IF
    IF(SpecDSMC(ClonedParticles(iPart,DelayCounter)%Species)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(ClonedParticles(iPart,DelayCounter)%Species)%SpecToPolyArray
      IF(ALLOCATED(VibQuantsPar(PositionNbr)%Quants)) DEALLOCATE(VibQuantsPar(PositionNbr)%Quants)
      ALLOCATE(VibQuantsPar(PositionNbr)%Quants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF))
      VibQuantsPar(PositionNbr)%Quants(:) = ClonedParticles(iPart,DelayCounter)%VibQuants(:)
    END IF
  END IF
  PartSpecies(PositionNbr) = ClonedParticles(iPart,DelayCounter)%Species
  ! Set the global element number with the offset
  PEM%GlobalElemID(PositionNbr) = ClonedParticles(iPart,DelayCounter)%Element
  PEM%LastGlobalElemID(PositionNbr) = PEM%GlobalElemID(PositionNbr)
  locElemID = PEM%LocalElemID(PositionNbr)
  LastPartPos(1:3,PositionNbr) = ClonedParticles(iPart,DelayCounter)%LastPartPos(1:3)
  PartMPF(PositionNbr) =  ClonedParticles(iPart,DelayCounter)%WeightingFactor
  IF (UseVarTimeStep) THEN
    PartTimeStep(PositionNbr) = GetParticleTimeStep(PartState(1,PositionNbr),PartState(2,PositionNbr),locElemID)
  END IF
  ! Counting the number of clones per cell
  IF(SamplingActive.OR.WriteMacroVolumeValues) THEN
    IF(DSMC%CalcQualityFactors) DSMC%QualityFacSamp(locElemID,5) = DSMC%QualityFacSamp(locElemID,5) + 1
  END IF
  IF(CalcPartBalance) THEN
    nPartIn(PartSpecies(PositionNbr))=nPartIn(PartSpecies(PositionNbr)) + 1
  END IF ! CalcPartBalance
END DO

! 3.) Reset the list
RadialWeighting%ClonePartNum(DelayCounter) = 0

! 3.1) Reduce ClonedParticles if necessary
CALL ReduceClonedParticlesType()

END SUBROUTINE DSMC_2D_SetInClones


REAL FUNCTION DSMC_2D_CalcSymmetryArea(iLocSide,iElem, ymin, ymax)
!===================================================================================================================================
!> Calculates the actual area of an element for 2D simulations (plane/axisymmetric) regardless of the mesh dimension in z
!> Utilized in the particle emission (surface flux) and boundary sampling
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars          ,ONLY: Pi
USE MOD_Particle_Vars         ,ONLY: Symmetry
USE MOD_Particle_Mesh_Vars    ,ONLY: NodeCoords_Shared, ElemSideNodeID_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: iElem,iLocSide           !> iElem is the compute-node element ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, OPTIONAL, INTENT(OUT)   :: ymax,ymin
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iNode
REAL                          :: P(1:2,1:4), Pmin(2), Pmax(2), Length, MidPoint
!===================================================================================================================================

Pmin = HUGE(Pmin)
Pmax = -HUGE(Pmax)

DO iNode = 1,4
  P(1:2,iNode) = NodeCoords_Shared(1:2,ElemSideNodeID_Shared(iNode,iLocSide,iElem)+1)
END DO

Pmax(1) = MAXVAL(P(1,:))
Pmax(2) = MAXVAL(P(2,:))
Pmin(1) = MINVAL(P(1,:))
Pmin(2) = MINVAL(P(2,:))

IF (PRESENT(ymax).AND.PRESENT(ymin)) THEN
  ymin = Pmin(2)
  ymax = Pmax(2)
END IF

Length = SQRT((Pmax(1)-Pmin(1))**2 + (Pmax(2)-Pmin(2))**2)

MidPoint = (Pmax(2)+Pmin(2)) / 2.
IF(Symmetry%Axisymmetric) THEN
  DSMC_2D_CalcSymmetryArea = Length * MidPoint * Pi * 2.
  ! Area of the cells on the rotational symmetry axis is set to one
  IF(.NOT.(DSMC_2D_CalcSymmetryArea.GT.0.0)) DSMC_2D_CalcSymmetryArea = 1.
ELSE
  DSMC_2D_CalcSymmetryArea = Length
END IF
RETURN

END FUNCTION DSMC_2D_CalcSymmetryArea


REAL FUNCTION DSMC_1D_CalcSymmetryArea(iLocSide,iElem)
!===================================================================================================================================
!> Calculates the actual area of an element for 1D simulations regardless of the mesh dimension in z and y
!> Utilized in the particle emission (surface flux) and boundary sampling
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Mesh_Vars    ,ONLY: NodeCoords_Shared, ElemSideNodeID_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: iElem,iLocSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: Pmin, Pmax, Length
!===================================================================================================================================

Pmax = MAXVAL(NodeCoords_Shared(1,ElemSideNodeID_Shared(:,iLocSide,iElem)+1))
Pmin = MINVAL(NodeCoords_Shared(1,ElemSideNodeID_Shared(:,iLocSide,iElem)+1))

Length = ABS(Pmax-Pmin)

! IF(Symmetry%Axisymmetric) THEN
!   MidPoint = (Pmax(2)+Pmin(2)) / 2.
!   DSMC_1D_CalcSymmetryArea = Length * MidPoint * Pi * 2.
!   ! Area of the cells on the rotational symmetry axis is set to one
!   IF(.NOT.(DSMC_1D_CalcSymmetryArea.GT.0.0)) DSMC_1D_CalcSymmetryArea = 1.
! ELSE
  DSMC_1D_CalcSymmetryArea = Length
  IF (DSMC_1D_CalcSymmetryArea.EQ.0.) DSMC_1D_CalcSymmetryArea = 1.
! END IF
RETURN

END FUNCTION DSMC_1D_CalcSymmetryArea


FUNCTION DSMC_2D_CalcSymmetryAreaSubSides(iLocSide,iElem)
!===================================================================================================================================
!> Calculates the area of the subsides for the insertion with the surface flux
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars              ,ONLY: Pi
USE MOD_DSMC_Vars                 ,ONLY: RadialWeighting
USE MOD_Particle_Mesh_Vars        ,ONLY: NodeCoords_Shared,ElemSideNodeID_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: iLocSide,iElem           !> iElem is the compute-node element ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                              :: DSMC_2D_CalcSymmetryAreaSubSides(RadialWeighting%nSubSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: iNode
REAL                              :: P(1:2,1:4), Pmin(2), Pmax(2), MidPoint, PminTemp, PmaxTemp, Length
!===================================================================================================================================

Pmin = HUGE(Pmin)
Pmax = -HUGE(Pmax)

DO iNode = 1,4
  P(1:2,iNode) = NodeCoords_Shared(1:2,ElemSideNodeID_Shared(iNode,iLocSide,iElem)+1)
END DO

Pmax(1) = MAXVAL(P(1,:))
Pmax(2) = MAXVAL(P(2,:))
Pmin(1) = MINVAL(P(1,:))
Pmin(2) = MINVAL(P(2,:))
Length = SQRT((Pmax(1)-Pmin(1))**2 + (Pmax(2)-Pmin(2))**2)

DO iNode = 1, RadialWeighting%nSubSides
  PminTemp = Pmin(2) + (Pmax(2) - Pmin(2))/RadialWeighting%nSubSides*(iNode-1.)
  PmaxTemp = Pmin(2) + (Pmax(2) - Pmin(2))/RadialWeighting%nSubSides*iNode
  MidPoint = (PmaxTemp+PminTemp) / 2.
  DSMC_2D_CalcSymmetryAreaSubSides(iNode) = Length/RadialWeighting%nSubSides * MidPoint * Pi * 2.
END DO

RETURN

END FUNCTION DSMC_2D_CalcSymmetryAreaSubSides


SUBROUTINE Init_Symmetry()
!===================================================================================================================================
!> Initialize if a 2D/1D Simulation is performed and which type
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars    ,ONLY: Symmetry
USE MOD_DSMC_Vars        ,ONLY: RadialWeighting
USE MOD_ReadInTools      ,ONLY: GETLOGICAL,GETINT
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                 :: Symmetry2D
!===================================================================================================================================
Symmetry%Order = GETINT('Particles-Symmetry-Order')
Symmetry2D = GETLOGICAL('Particles-Symmetry2D')
IF(Symmetry2D.AND.(Symmetry%Order.EQ.3)) THEN
  Symmetry%Order = 2
  LBWRITE(*,*) 'WARNING: Particles-Symmetry-Order is set to 2 because of Particles-Symmetry2D=.TRUE. .'
  LBWRITE(*,*) 'Set Particles-Symmetry-Order=2 and remove Particles-Symmetry2D to avoid this warning'
ELSE IF(Symmetry2D) THEN
  CALL ABORT(__STAMP__&
    ,'ERROR: 2D Simulations either with Particles-Symmetry-Order=2 or (but not recommended) with Symmetry2D=.TRUE.')
END IF

IF((Symmetry%Order.LE.0).OR.(Symmetry%Order.GE.4)) CALL ABORT(__STAMP__&
,'Particles-Symmetry-Order (space dimension) has to be in the range of 1 to 3')

Symmetry%Axisymmetric = GETLOGICAL('Particles-Symmetry2DAxisymmetric')
IF(Symmetry%Axisymmetric.AND.(Symmetry%Order.EQ.3)) CALL ABORT(__STAMP__&
  ,'ERROR: Axisymmetric simulations only for 1D or 2D')
IF(Symmetry%Axisymmetric.AND.(Symmetry%Order.EQ.1))CALL ABORT(__STAMP__&
  ,'ERROR: Axisymmetric simulations are only implemented for Particles-Symmetry-Order=2 !')
IF(Symmetry%Axisymmetric) THEN
  RadialWeighting%DoRadialWeighting = GETLOGICAL('Particles-RadialWeighting')
ELSE
  RadialWeighting%DoRadialWeighting = .FALSE.
  RadialWeighting%PerformCloning = .FALSE.
END IF

END SUBROUTINE Init_Symmetry


SUBROUTINE DSMC_2D_TreatIdenticalParticles(iPair, nPair, nPart, iElem, iPartIndx_Node)
!===================================================================================================================================
!> Check if particle pairs have a zero relative velocity (and thus a collision probability of zero), if they do, break up the pair
!> and use either a left-over particle (uneven number of particles in a cell) or swap the collision partners with the next pair in
!> the list.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars             ,ONLY: Coll_pData, DSMC, ChemReac, SamplingActive, CollInf, CollisMode
USE MOD_Particle_Vars         ,ONLY: PartSpecies, PartState, WriteMacroVolumeValues
USE MOD_part_tools            ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPair
INTEGER, INTENT(IN)           :: nPair
INTEGER, INTENT(IN)           :: nPart
INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(INOUT)        :: iPartIndx_Node(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart_p1, iPart_p2, tempPart, cSpec1, cSpec2, iCase
!===================================================================================================================================
! Two particles with the exact same velocities at the same positions -> clones that did not interact with other particles/walls
IF (Coll_pData(iPair)%CRela2.EQ.0.0) THEN
  IF ((CollisMode.LT.3).AND.(nPart.EQ.1)) THEN
    ! Uneven number of particles in the cell, a single particle is left without a pair
    ! Removing the pairs from the weighting factor and the case num sums
    CollInf%SumPairMPF(Coll_pData(iPair)%PairType) = CollInf%SumPairMPF(Coll_pData(iPair)%PairType) &
      -(GetParticleWeight(Coll_pData(iPair)%iPart_p1) + GetParticleWeight(Coll_pData(iPair)%iPart_p2))*0.5
    CollInf%Coll_CaseNum(Coll_pData(iPair)%PairType) = CollInf%Coll_CaseNum(Coll_pData(iPair)%PairType) - 1
    ! Swapping particle without a pair with the first particle of the current pair
    tempPart = Coll_pData(iPair)%iPart_p1
    Coll_pData(iPair)%iPart_p1 = iPartIndx_Node(1)
    iPartIndx_Node(1) = tempPart
    IF (CollisMode.EQ.3) ChemReac%RecombParticle = iPartIndx_Node(1)
    IF (CollInf%ProhibitDoubleColl)  CollInf%OldCollPartner(iPartIndx_Node(1)) = 0
    ! Increase the appropriate case number and set the right pair type
    iPart_p1 = Coll_pData(iPair)%iPart_p1; iPart_p2 = Coll_pData(iPair)%iPart_p2
    cSpec1 = PartSpecies(iPart_p1); cSpec2 = PartSpecies(iPart_p2)
    iCase = CollInf%Coll_Case(cSpec1, cSpec2)
    CollInf%Coll_CaseNum(iCase) = CollInf%Coll_CaseNum(iCase) + 1
    Coll_pData(iPair)%PairType = iCase
    ! Adding the pair to the sums of the number of collisions (with and without weighting factor)
    CollInf%SumPairMPF(iCase) = CollInf%SumPairMPF(iCase) + (GetParticleWeight(iPart_p1) + GetParticleWeight(iPart_p2))*0.5
    ! Calculation of the relative velocity for the new first pair
    Coll_pData(iPair)%CRela2 = (PartState(4,iPart_p1) - PartState(4,iPart_p2))**2 &
                             + (PartState(5,iPart_p1) - PartState(5,iPart_p2))**2 &
                             + (PartState(6,iPart_p1) - PartState(6,iPart_p2))**2
  ELSE IF (iPair.LT.nPair) THEN
    IF (.NOT.Coll_pData(iPair+1)%NeedForRec) THEN
    ! "Partner-Tausch": if there are pairs ahead in the pairing list, the next is pair is broken up and collision partners
    ! are swapped
      CollInf%Coll_CaseNum(Coll_pData(iPair)%PairType) = CollInf%Coll_CaseNum(Coll_pData(iPair)%PairType) - 1
      CollInf%SumPairMPF(Coll_pData(iPair)%PairType) = CollInf%SumPairMPF(Coll_pData(iPair)%PairType) &
        - 0.5 * (GetParticleWeight(Coll_pData(iPair)%iPart_p1) + GetParticleWeight(Coll_pData(iPair)%iPart_p2))
      CollInf%Coll_CaseNum(Coll_pData(iPair+1)%PairType) = CollInf%Coll_CaseNum(Coll_pData(iPair+1)%PairType) - 1
      CollInf%SumPairMPF(Coll_pData(iPair+1)%PairType) = CollInf%SumPairMPF(Coll_pData(iPair+1)%PairType) &
        - 0.5 * (GetParticleWeight(Coll_pData(iPair+1)%iPart_p1) + GetParticleWeight(Coll_pData(iPair+1)%iPart_p2))
      ! Breaking up the next pair and swapping partners
      tempPart = Coll_pData(iPair)%iPart_p1
      Coll_pData(iPair)%iPart_p1 = Coll_pData(iPair + 1)%iPart_p1
      Coll_pData(iPair + 1)%iPart_p1 = tempPart
      ! Increase the appropriate case number and set the right pair type
      iPart_p1 = Coll_pData(iPair)%iPart_p1; iPart_p2 = Coll_pData(iPair)%iPart_p2
      cSpec1 = PartSpecies(iPart_p1); cSpec2 = PartSpecies(iPart_p2)
      iCase = CollInf%Coll_Case(cSpec1, cSpec2)
      CollInf%Coll_CaseNum(iCase) = CollInf%Coll_CaseNum(iCase) + 1
      CollInf%SumPairMPF(iCase) = CollInf%SumPairMPF(iCase) + 0.5 * (GetParticleWeight(iPart_p1) + GetParticleWeight(iPart_p2))
      Coll_pData(iPair)%PairType = iCase
      ! Calculation of the relative velocity for the new first pair
      Coll_pData(iPair)%CRela2 = (PartState(4,iPart_p1) - PartState(4,iPart_p2))**2 &
                               + (PartState(5,iPart_p1) - PartState(5,iPart_p2))**2 &
                               + (PartState(6,iPart_p1) - PartState(6,iPart_p2))**2
      CollInf%SumPairMPF(iCase) = CollInf%SumPairMPF(iCase) + 0.5 * (GetParticleWeight(iPart_p1) + GetParticleWeight(iPart_p2))
      IF(Coll_pData(iPair)%CRela2.EQ.0.0) THEN
        ! If the relative velocity is still zero, add the pair to the identical particles count
        IF(SamplingActive.OR.WriteMacroVolumeValues) THEN
          IF(DSMC%CalcQualityFactors) DSMC%QualityFacSamp(iElem,6) = DSMC%QualityFacSamp(iElem,6) + 1
        END IF
      END IF
      ! Increase the appropriate case number and set the right pair type
      iPart_p1 = Coll_pData(iPair+1)%iPart_p1; iPart_p2 = Coll_pData(iPair+1)%iPart_p2
      cSpec1 = PartSpecies(iPart_p1); cSpec2 = PartSpecies(iPart_p2)
      iCase = CollInf%Coll_Case(cSpec1, cSpec2)
      CollInf%Coll_CaseNum(iCase) = CollInf%Coll_CaseNum(iCase) + 1
      CollInf%SumPairMPF(iCase) = CollInf%SumPairMPF(iCase) + 0.5 * (GetParticleWeight(iPart_p1) + GetParticleWeight(iPart_p2))
      ! Calculation of the relative velocity for the new follow-up pair
      Coll_pData(iPair+1)%CRela2 = (PartState(4,iPart_p1) - PartState(4,iPart_p2))**2 &
                                 + (PartState(5,iPart_p1) - PartState(5,iPart_p2))**2 &
                                 + (PartState(6,iPart_p1) - PartState(6,iPart_p2))**2
      Coll_pData(iPair+1)%PairType = iCase
      CollInf%SumPairMPF(iCase) = CollInf%SumPairMPF(iCase) + 0.5 * (GetParticleWeight(iPart_p1) + GetParticleWeight(iPart_p2))
    ELSE
      ! For the last pair, only invert the velocity in z and calculate new relative velocity
      iPart_p1 = Coll_pData(iPair)%iPart_p1; iPart_p2 = Coll_pData(iPair)%iPart_p2
      PartState(6,iPart_p1) = - PartState(6,iPart_p1)
      Coll_pData(iPair)%CRela2 = (PartState(4,iPart_p1) - PartState(4,iPart_p2))**2 &
                               + (PartState(5,iPart_p1) - PartState(5,iPart_p2))**2 &
                               + (PartState(6,iPart_p1) - PartState(6,iPart_p2))**2
      ! Add the pair to the number of identical particles output
      IF(SamplingActive.OR.WriteMacroVolumeValues) THEN
        IF(DSMC%CalcQualityFactors) DSMC%QualityFacSamp(iElem,6) = DSMC%QualityFacSamp(iElem,6) + 1
      END IF
    END IF  ! nPart.EQ.1/iPair.LT.nPair
  ELSE
    ! For the last pair, only invert the velocity in z and calculate new relative velocity
    iPart_p1 = Coll_pData(iPair)%iPart_p1; iPart_p2 = Coll_pData(iPair)%iPart_p2
    PartState(6,iPart_p1) = - PartState(6,iPart_p1)
    Coll_pData(iPair)%CRela2 = (PartState(4,iPart_p1) - PartState(4,iPart_p2))**2 &
                             + (PartState(5,iPart_p1) - PartState(5,iPart_p2))**2 &
                             + (PartState(6,iPart_p1) - PartState(6,iPart_p2))**2
    ! Add the pair to the number of identical particles output
    IF(SamplingActive.OR.WriteMacroVolumeValues) THEN
      IF(DSMC%CalcQualityFactors) DSMC%QualityFacSamp(iElem,6) = DSMC%QualityFacSamp(iElem,6) + 1
    END IF
  END IF  ! nPart.EQ.1/iPair.LT.nPair
END IF    ! Coll_pData(iPair)%CRela2.EQ.0.0

END SUBROUTINE DSMC_2D_TreatIdenticalParticles


SUBROUTINE IncreaseClonedParticlesType()
!===================================================================================================================================
!> Increases RadialWeighting%CloneVecLength and the ClonedParticles(iPart,iDelay) type
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars             ,ONLY: RadialWeighting, ClonedParticles, tClonedParticles
USE MOD_Particle_Vars         ,ONLY: PDM
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: NewSize,i,ii,ALLOCSTAT
TYPE (tClonedParticles), ALLOCATABLE :: ClonedParticles_new(:,:)
!===================================================================================================================================

NewSize = MAX(CEILING(RadialWeighting%CloneVecLength * (1+PDM%MaxPartNumIncrease)),RadialWeighting%CloneVecLength+RadialWeighting%CloneVecLengthDelta)

SELECT CASE(RadialWeighting%CloneMode)
CASE(1)
  ALLOCATE(ClonedParticles_new(1:NewSize,0:(RadialWeighting%CloneInputDelay-1)),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL ABORT(&
__STAMP__&
,'Cannot allocate increased new ClonedParticles Array')
  DO ii=0,RadialWeighting%CloneInputDelay-1
    DO i=1,RadialWeighting%ClonePartNum(ii)
      ClonedParticles_new(i,ii)%Species=ClonedParticles(i,ii)%Species
      ClonedParticles_new(i,ii)%PartState(1:6)=ClonedParticles(i,ii)%PartState(1:6)
      ClonedParticles_new(i,ii)%PartStateIntEn(1:3)=ClonedParticles(i,ii)%PartStateIntEn(1:3)
      ClonedParticles_new(i,ii)%Element=ClonedParticles(i,ii)%Element
      ClonedParticles_new(i,ii)%LastPartPos(1:3)=ClonedParticles(i,ii)%LastPartPos(1:3)
      ClonedParticles_new(i,ii)%WeightingFactor=ClonedParticles(i,ii)%WeightingFactor
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%VibQuants,ClonedParticles_new(i,ii)%VibQuants)
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%DistriFunc,ClonedParticles_new(i,ii)%DistriFunc)
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%AmbiPolVelo,ClonedParticles_new(i,ii)%AmbiPolVelo)
    END DO
  END DO
  DEALLOCATE(ClonedParticles)
  CALL MOVE_ALLOC(ClonedParticles_New,ClonedParticles)
CASE(2)
  ALLOCATE(ClonedParticles_new(1:NewSize,0:RadialWeighting%CloneInputDelay),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL ABORT(&
__STAMP__&
,'Cannot allocate increased new ClonedParticles Array')
  DO ii=0,RadialWeighting%CloneInputDelay
    DO i=1,RadialWeighting%ClonePartNum(ii)
      ClonedParticles_new(i,ii)%Species=ClonedParticles(i,ii)%Species
      ClonedParticles_new(i,ii)%PartState(1:6)=ClonedParticles(i,ii)%PartState(1:6)
      ClonedParticles_new(i,ii)%PartStateIntEn(1:3)=ClonedParticles(i,ii)%PartStateIntEn(1:3)
      ClonedParticles_new(i,ii)%Element=ClonedParticles(i,ii)%Element
      ClonedParticles_new(i,ii)%LastPartPos(1:3)=ClonedParticles(i,ii)%LastPartPos(1:3)
      ClonedParticles_new(i,ii)%WeightingFactor=ClonedParticles(i,ii)%WeightingFactor
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%VibQuants,ClonedParticles_new(i,ii)%VibQuants)
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%DistriFunc,ClonedParticles_new(i,ii)%DistriFunc)
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%AmbiPolVelo,ClonedParticles_new(i,ii)%AmbiPolVelo)
    END DO
  END DO
  DEALLOCATE(ClonedParticles)
  CALL MOVE_ALLOC(ClonedParticles_New,ClonedParticles)
CASE DEFAULT
  CALL Abort(&
      __STAMP__,&
    'ERROR in Radial Weighting of 2D/Axisymmetric: The selected cloning mode is not available! Choose between 1 and 2.'//&
      ' CloneMode=1: Delayed insertion of clones; CloneMode=2: Delayed randomized insertion of clones')
END SELECT

RadialWeighting%CloneVecLength = NewSize

END SUBROUTINE IncreaseClonedParticlesType


SUBROUTINE ReduceClonedParticlesType()
!===================================================================================================================================
!> Reduces RadialWeighting%CloneVecLength and the ClonedParticles(iPart,iDelay) type
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars             ,ONLY: RadialWeighting, ClonedParticles, tClonedParticles
USE MOD_Particle_Vars         ,ONLY: PDM
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: NewSize,i,ii,ALLOCSTAT
TYPE (tClonedParticles), ALLOCATABLE :: ClonedParticles_new(:,:)
!===================================================================================================================================
IF (MAXVAL(RadialWeighting%ClonePartNum(:)).GE.PDM%maxParticleNumber/(1.+PDM%MaxPartNumIncrease)**2) RETURN

NewSize = MAX(CEILING(MAXVAL(RadialWeighting%ClonePartNum(:))*(1.+PDM%MaxPartNumIncrease)),1)

IF (NewSize.GT.RadialWeighting%CloneVecLength-RadialWeighting%CloneVecLengthDelta) RETURN

SELECT CASE(RadialWeighting%CloneMode)
CASE(1)
  ALLOCATE(ClonedParticles_new(1:NewSize,0:(RadialWeighting%CloneInputDelay-1)),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL ABORT(&
__STAMP__&
,'Cannot allocate increased new ClonedParticles Array')
  DO ii=0,RadialWeighting%CloneInputDelay-1
    DO i=1,RadialWeighting%ClonePartNum(ii)
      ClonedParticles_new(i,ii)%Species=ClonedParticles(i,ii)%Species
      ClonedParticles_new(i,ii)%PartState(1:6)=ClonedParticles(i,ii)%PartState(1:6)
      ClonedParticles_new(i,ii)%PartStateIntEn(1:3)=ClonedParticles(i,ii)%PartStateIntEn(1:3)
      ClonedParticles_new(i,ii)%Element=ClonedParticles(i,ii)%Element
      ClonedParticles_new(i,ii)%LastPartPos(1:3)=ClonedParticles(i,ii)%LastPartPos(1:3)
      ClonedParticles_new(i,ii)%WeightingFactor=ClonedParticles(i,ii)%WeightingFactor
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%VibQuants,ClonedParticles_new(i,ii)%VibQuants)
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%DistriFunc,ClonedParticles_new(i,ii)%DistriFunc)
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%AmbiPolVelo,ClonedParticles_new(i,ii)%AmbiPolVelo)
    END DO
  END DO
  DEALLOCATE(ClonedParticles)
  CALL MOVE_ALLOC(ClonedParticles_New,ClonedParticles)
CASE(2)
  ALLOCATE(ClonedParticles_new(1:NewSize,0:RadialWeighting%CloneInputDelay),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL ABORT(&
__STAMP__&
,'Cannot allocate increased new ClonedParticles Array')
  DO ii=0,RadialWeighting%CloneInputDelay
    DO i=1,RadialWeighting%ClonePartNum(ii)
      ClonedParticles_new(i,ii)%Species=ClonedParticles(i,ii)%Species
      ClonedParticles_new(i,ii)%PartState(1:6)=ClonedParticles(i,ii)%PartState(1:6)
      ClonedParticles_new(i,ii)%PartStateIntEn(1:3)=ClonedParticles(i,ii)%PartStateIntEn(1:3)
      ClonedParticles_new(i,ii)%Element=ClonedParticles(i,ii)%Element
      ClonedParticles_new(i,ii)%LastPartPos(1:3)=ClonedParticles(i,ii)%LastPartPos(1:3)
      ClonedParticles_new(i,ii)%WeightingFactor=ClonedParticles(i,ii)%WeightingFactor
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%VibQuants,ClonedParticles_new(i,ii)%VibQuants)
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%DistriFunc,ClonedParticles_new(i,ii)%DistriFunc)
      CALL MOVE_ALLOC(ClonedParticles(i,ii)%AmbiPolVelo,ClonedParticles_new(i,ii)%AmbiPolVelo)
    END DO
  END DO
  DEALLOCATE(ClonedParticles)
  CALL MOVE_ALLOC(ClonedParticles_New,ClonedParticles)
CASE DEFAULT
  CALL Abort(&
      __STAMP__,&
    'ERROR in Radial Weighting of 2D/Axisymmetric: The selected cloning mode is not available! Choose between 1 and 2.'//&
      ' CloneMode=1: Delayed insertion of clones; CloneMode=2: Delayed randomized insertion of clones')
END SELECT

RadialWeighting%CloneVecLength = NewSize

END SUBROUTINE ReduceClonedParticlesType

END MODULE MOD_DSMC_Symmetry
