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

MODULE MOD_Particle_Boundary_Sampling
!===================================================================================================================================
!> Particle boundary sampling: calculation and output of heat flux, forces, and impact properties at boundaries
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
PUBLIC::CalcSurfaceValues
PUBLIC::WriteSurfSampleToHDF5
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
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfTotalSideOnNode
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfSample,dXiEQ_SurfSample,PartBound,XiEQ_SurfSample
USE MOD_Particle_Boundary_Vars  ,ONLY: nComputeNodeSurfSides,nComputeNodeSurfTotalSides,nComputeNodeSurfOutputSides
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfBC,SurfBCName
USE MOD_Particle_Boundary_Vars  ,ONLY: nGlobalSurfSides
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfSide2GlobalSide
USE MOD_SurfaceModel_Vars       ,ONLY: nPorousBC
USE MOD_Particle_Boundary_Vars  ,ONLY: CalcSurfaceImpact
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfSideArea,SurfSampSize,SurfOutputSize,SurfSpecOutputSize
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallState, SWIVarTimeStep, SWIStickingCoefficient
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallPumpCapacity
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallImpactEnergy
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallImpactVector
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallImpactAngle
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallImpactNumber
USE MOD_Particle_Boundary_Vars  ,ONLY: MacroSurfaceVal,MacroSurfaceSpecVal
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared,NodeCoords_Shared
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemSideNodeID_Shared
USE MOD_Particle_Surfaces       ,ONLY: EvaluateBezierPolynomialAndGradient
USE MOD_Particle_Surfaces_Vars  ,ONLY: BezierControlPoints3D
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
USE MOD_Particle_Vars           ,ONLY: nSpecies,UseVarTimeStep,VarTimeStep
USE MOD_Particle_Vars           ,ONLY: Symmetry
USE MOD_ReadInTools             ,ONLY: GETINT,GETLOGICAL,GETINTARRAY
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_LEADERS_SURF,mySurfRank
USE MOD_MPI_Shared_Vars         ,ONLY: myComputeNodeRank,nComputeNodeProcessors
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfSideArea_Shared,SurfSideArea_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallState_Shared,SampWallState_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallPumpCapacity_Shared,SampWallPumpCapacity_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallImpactEnergy_Shared,SampWallImpactEnergy_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallImpactVector_Shared,SampWallImpactVector_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallImpactAngle_Shared,SampWallImpactAngle_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallImpactNumber_Shared,SampWallImpactNumber_Shared_Win
USE MOD_Particle_MPI_Boundary_Sampling,ONLY: InitSurfCommunication
#else
USE MOD_MPI_Shared_Vars         ,ONLY: mySurfRank
USE MOD_Particle_Mesh_Vars      ,ONLY: nComputeNodeSides
USE MOD_Particle_Boundary_Vars  ,ONLY: nGlobalOutputSides
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars        ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                :: iBC
INTEGER                                :: iSide,firstSide,lastSide
CHARACTER(LEN=255),ALLOCATABLE         :: BCName(:)
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
!===================================================================================================================================

! Get input parameters
LBWRITE(UNIT_stdOut,'(A)') ' INIT SURFACE SAMPLING ...'

! Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z) and angle
CalcSurfaceImpact = GETLOGICAL('CalcSurfaceImpact')

! Flag if there is at least one surf side on the node (sides in halo region do also count)
SurfTotalSideOnNode = MERGE(.TRUE.,.FALSE.,nComputeNodeSurfTotalSides.GT.0)

!> Setting the number of sampling (SurfSampSize -> SampWallState) and output (SurfOutputSize -> MacroSurfaceVal) variables
!> Optional sampling variables require an additional SampWallIndex (SWI)
! Default: Energy + Force + nSpecies
SurfSampSize = SAMPWALL_NVARS+nSpecies
! Default: Heatflux + Force + Total impact counter + iBC
SurfOutputSize = MACROSURF_NVARS
! Optional variables (number of sampling and output variables can differ)
! Variable time step (required for correct heat flux calculation)
IF(UseVarTimeStep.OR.VarTimeStep%UseSpeciesSpecific) THEN
  SurfSampSize = SurfSampSize + 1
  SWIVarTimeStep = SurfSampSize
END IF
! Sticking coefficient (empirical model)
IF(ANY(PartBound%SurfaceModel.EQ.1)) THEN
  SurfSampSize = SurfSampSize + 1
  SWIStickingCoefficient = SurfSampSize
  SurfOutputSize = SurfOutputSize + 1
END IF
! Pump capacity (Porous BC)
IF (nPorousBC.GT.0) SurfOutputSize = SurfOutputSize + nPorousBC
! Wall temperature (Adaptive radiative-equilibrium BC)
IF (PartBound%OutputWallTemp) SurfOutputSize = SurfOutputSize + 1
!> Additional species-specific output (SurfSpecOutputSize -> MacroSurfaceSpecVal)
! Species-specific counter of impacts
SurfSpecOutputSize = 1
! Sampling of impact energy for each species (trans, rot, vib, elec), impact vector (x,y,z), angle, total number, number per second: Add 10 variables
IF (CalcSurfaceImpact) SurfSpecOutputSize = SurfSpecOutputSize + 10

!> Leader communication
#if USE_MPI
IF (myComputeNodeRank.EQ.0) CALL InitSurfCommunication()
! The leaders are synchronized at this point, but behind the other procs. nGlobalSurfSides is only required when compiled without
! MPI, so perform latency hiding by postponing synchronization
#else
mySurfRank         = 0
nGlobalSurfSides   = nComputeNodeSurfTotalSides
nGlobalOutputSides = nComputeNodeSurfOutputSides
#endif /* USE_MPI */

! surface sampling array do not need to be allocated if there are no sides within halo_eps range
IF(.NOT.SurfTotalSideOnNode) RETURN

!> Allocate the output container
ALLOCATE(MacroSurfaceVal(1:SurfOutputSize         , 1:nSurfSample , 1:nSurfSample , nComputeNodeSurfOutputSides))
MacroSurfaceVal     = 0.
ALLOCATE(MacroSurfaceSpecVal(1:SurfSpecOutputSize , 1:nSurfSample , 1:nSurfSample , nComputeNodeSurfOutputSides , nSpecies))
MacroSurfaceSpecVal = 0.

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
  IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(iBC)).EQ.PartBound%ReflectiveBC) THEN
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

firstSide = INT(REAL( myComputeNodeRank   )*REAL(nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1))*REAL(nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))
#else
ALLOCATE(SurfSideArea(1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides))

firstSide = 1
lastSide  = nGlobalSurfSides
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

  UseBezierControlPointsForArea = .FALSE.

  IF (TrackingMethod.EQ.TRIATRACKING) THEN
    ElemID    = SideInfo_Shared(SIDE_ELEMID ,SideID)
    CNElemID  = GetCNElemID(ElemID)
    LocSideID = SideInfo_Shared(SIDE_LOCALID,SideID)
    IF((Symmetry%Order.NE.3).AND.nSurfSample.GT.1) CALL abort(__STAMP__,'nSurfSample>1 not implemented for this symmetry!')

    IF(Symmetry%Order.EQ.3) THEN
      ! Check if triangles are used for the calculation of the surface area or not
      IF(nSurfSample.GT.1)THEN
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
      SurfSideArea(1,1,iSide) = area
      END IF ! nSurfSample.GT.1
    ELSE IF(Symmetry%Order.EQ.2) THEN
      SurfSideArea(1,1,iSide) = DSMC_2D_CalcSymmetryArea(LocSideID, CNElemID)
    ELSE IF(Symmetry%Order.EQ.1) THEN
      SurfSideArea(1,1,iSide) = DSMC_1D_CalcSymmetryArea(LocSideID, CNElemID)
    END IF
  ELSE ! TrackingMethod.NE.TRIATRACKING
    UseBezierControlPointsForArea = .TRUE.
  END IF ! TrackingMethod.EQ.TRIATRACKIN

  ! Instead of triangles use Bezier control points (curved or triangle tracking with nSurfSample>1)
  IF(UseBezierControlPointsForArea)THEN
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
  END IF ! UseBezierControlPointsForArea

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
CALL MPI_BCAST(nGlobalSurfSides,1,MPI_INTEGER,0,MPI_COMM_SHARED,iError)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

IF (mySurfRank.EQ.0) THEN
#endif
  LBWRITE(UNIT_StdOut,'(A,I8)')       ' | Number of sampling sides:           '    , nGlobalSurfSides
  LBWRITE(UNIT_StdOut,'(A,ES10.4E2)') ' | Surface-Area:                         ', Area
  LBWRITE(UNIT_stdOut,'(A)') ' INIT SURFACE SAMPLING DONE'
#if USE_MPI
END IF
#endif

END SUBROUTINE InitParticleBoundarySampling


SUBROUTINE CalcSurfaceValues(during_dt_opt)
!===================================================================================================================================
!> Calculates macroscopic surface values from samples
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars               ,ONLY: StefanBoltzmannConst
USE MOD_DSMC_Vars                  ,ONLY: DSMC
USE MOD_Mesh_Vars                  ,ONLY: MeshFile
USE MOD_Particle_Boundary_Vars     ,ONLY: SurfTotalSideOnNode
USE MOD_SurfaceModel_Vars          ,ONLY: nPorousBC
USE MOD_Particle_Boundary_Vars     ,ONLY: nSurfSample,CalcSurfaceImpact
USE MOD_Particle_Boundary_Vars     ,ONLY: SurfSide2GlobalSide, GlobalSide2SurfSide, PartBound
USE MOD_Particle_Boundary_Vars     ,ONLY: nComputeNodeSurfSides, BoundaryWallTemp
USE MOD_Particle_Boundary_Vars     ,ONLY: PorousBCInfo_Shared,MapSurfSideToPorousSide_Shared
USE MOD_Particle_Boundary_vars     ,ONLY: SurfOutputSize, SWIVarTimeStep, SWIStickingCoefficient
USE MOD_Particle_Boundary_Vars     ,ONLY: MacroSurfaceVal, MacroSurfaceSpecVal
USE MOD_Particle_Mesh_Vars         ,ONLY: SideInfo_Shared
USE MOD_Particle_Vars              ,ONLY: WriteMacroSurfaceValues,nSpecies,MacroValSampTime,UseVarTimeStep,Symmetry,VarTimeStep
USE MOD_Particle_Vars              ,ONLY: Species
USE MOD_Restart_Vars               ,ONLY: RestartTime
USE MOD_TimeDisc_Vars              ,ONLY: TEnd
USE MOD_Timedisc_Vars              ,ONLY: time,dt
#if USE_MPI
USE MOD_MPI_Shared                 ,ONLY: BARRIER_AND_SYNC
USE MOD_MPI_Shared_Vars            ,ONLY: MPI_COMM_LEADERS_SURF, MPI_COMM_SHARED
USE MOD_Particle_Boundary_Vars     ,ONLY: SampWallPumpCapacity_Shared
USE MOD_Particle_Boundary_vars     ,ONLY: SampWallState_Shared,SampWallImpactNumber_Shared,SampWallImpactEnergy_Shared
USE MOD_Particle_Boundary_vars     ,ONLY: SampWallImpactVector_Shared,SampWallImpactAngle_Shared
USE MOD_Particle_Boundary_vars     ,ONLY: SurfSideArea_Shared
USE MOD_Particle_MPI_Boundary_Sampling,ONLY: ExchangeSurfData
USE MOD_Particle_Boundary_Vars    ,ONLY: BoundaryWallTemp_Shared_Win
#else
USE MOD_Particle_Boundary_Vars     ,ONLY: SampWallPumpCapacity
USE MOD_Particle_Boundary_vars     ,ONLY: SampWallState,SampWallImpactNumber,SampWallImpactEnergy
USE MOD_Particle_Boundary_vars     ,ONLY: SampWallImpactVector,SampWallImpactAngle
USE MOD_Particle_Boundary_vars     ,ONLY: SurfSideArea
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL, INTENT(IN), OPTIONAL      :: during_dt_opt !routine was called during timestep (i.e. before iter=iter+1, time=time+dt...)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iSpec,iSurfSide,p,q, iPBC, nVarCount, OutputCounter
REAL                               :: TimeSample, ActualTime, TimeSampleTemp, CounterSum, nImpacts, IterNum
LOGICAL                            :: during_dt
INTEGER                            :: idx, GlobalSideID, SurfSideNb, iBC
!===================================================================================================================================

IF (PRESENT(during_dt_opt)) THEN
  during_dt=during_dt_opt
ELSE
  during_dt=.FALSE.
END IF
IF (during_dt) THEN
  ActualTime=time+dt
ELSE
  ActualTime=time
END IF

! Determine the sampling time for the calculation of fluxes
IF (WriteMacroSurfaceValues) THEN
  ! Elapsed time since last sampling (variable dt's possible!)
  TimeSample = Time - MacroValSampTime
  MacroValSampTime = Time
ELSE IF (RestartTime.GT.(1-DSMC%TimeFracSamp)*TEnd) THEN
  ! Sampling at the end of the simulation: When a restart is performed and the sampling starts immediately, determine the correct sampling time
  ! (e.g. sampling is set to 20% of tend = 1s, and restart is performed at 0.9s, sample time = 0.1s)
  TimeSample = Time - RestartTime
ELSE
  ! Sampling at the end of the simulation: calculated from the user given input
  TimeSample = (Time-(1-DSMC%TimeFracSamp)*TEnd)
END IF

IF(ALMOSTZERO(TimeSample)) RETURN

IF(.NOT.SurfTotalSideOnNode) RETURN

#if USE_MPI
CALL ExchangeSurfData()

! Only surface sampling leaders take part in the remainder of this routine
IF (MPI_COMM_LEADERS_SURF.EQ.MPI_COMM_NULL) THEN
  IF (PartBound%OutputWallTemp) THEN
    CALL BARRIER_AND_SYNC(BoundaryWallTemp_Shared_Win,MPI_COMM_SHARED)
  END IF
  RETURN
END IF
#endif /*USE_MPI*/

#if USE_MPI
ASSOCIATE(SampWallState        => SampWallState_Shared           ,&
          SampWallImpactNumber => SampWallImpactNumber_Shared    ,&
          SampWallImpactEnergy => SampWallImpactEnergy_Shared    ,&
          SampWallImpactVector => SampWallImpactVector_Shared    ,&
          SampWallImpactAngle  => SampWallImpactAngle_Shared     ,&
          SampWallPumpCapacity => SampWallPumpCapacity_Shared    ,&
          SurfSideArea         => SurfSideArea_Shared)
#endif

OutputCounter = 0
! Determine the total number of iterations
IterNum = REAL(NINT(TimeSample / dt))

DO iSurfSide = 1,nComputeNodeSurfSides
  !================== INNER BC CHECK
  GlobalSideID = SurfSide2GlobalSide(SURF_SIDEID,iSurfSide)
  IF(SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID).GT.0) THEN
    IF(GlobalSideID.LT.SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID)) THEN
      SurfSideNb = GlobalSide2SurfSide(SURF_SIDEID,SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID))
      ! Add your contribution to my inner BC
      SampWallState(:,:,:,iSurfSide) = SampWallState(:,:,:,iSurfSide) + SampWallState(:,:,:,SurfSideNb)
    ELSE
      CYCLE
    END IF
  END IF
  !================== ROTATIONALLY PERIODIC BC CHECK
  IF(PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,GlobalSideID))).EQ.PartBound%RotPeriodicBC) CYCLE
  !================== INTER PLANE BC CHECK
  IF(PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,GlobalSideID))).EQ.PartBound%RotPeriodicInterPlaneBC) CYCLE

  OutputCounter = OutputCounter + 1

  DO q = 1,nSurfSample
    DO p = 1,nSurfSample
      ! --- Default output (force per area, heat flux, simulation particle impact per iteration, boundary index)
      CounterSum = SUM(SampWallState(SAMPWALL_NVARS+1:SAMPWALL_NVARS+nSpecies,p,q,iSurfSide))

      IF(CounterSum.GT.0.0) THEN
        ! Correct the sample time in the case of a cell local time step with the average time step factor for each side
        IF(UseVarTimeStep .OR. VarTimeStep%UseSpeciesSpecific) THEN
          TimeSampleTemp = TimeSample * SampWallState(SWIVarTimeStep,p,q,iSurfSide) / CounterSum
        ELSE
          TimeSampleTemp = TimeSample
        END IF

        ! Force per area in x,y,z-direction
        MacroSurfaceVal(1:3,p,q,OutputCounter) = SampWallState(SAMPWALL_DELTA_MOMENTUMX:SAMPWALL_DELTA_MOMENTUMZ,p,q,iSurfSide) &
                                              / (SurfSideArea(p,q,iSurfSide)*TimeSampleTemp)
        ! Deleting the y/z-component for 1D/2D/axisymmetric simulations
        IF(Symmetry%Order.LT.3) MacroSurfaceVal(Symmetry%Order+1:3,p,q,iSurfSide) = 0.
        ! Heat flux (energy difference per second per area -> W/m2)
        MacroSurfaceVal(4,p,q,OutputCounter) = (SampWallState(SAMPWALL_ETRANSOLD,p,q,iSurfSide)  &
                                          + SampWallState(SAMPWALL_EROTOLD  ,p,q,iSurfSide)  &
                                          + SampWallState(SAMPWALL_EVIBOLD  ,p,q,iSurfSide)  &
                                          + SampWallState(SAMPWALL_EELECOLD ,p,q,iSurfSide)  &
                                          - SampWallState(SAMPWALL_ETRANSNEW,p,q,iSurfSide)  &
                                          - SampWallState(SAMPWALL_EROTNEW  ,p,q,iSurfSide)  &
                                          - SampWallState(SAMPWALL_EVIBNEW  ,p,q,iSurfSide)  &
                                          - SampWallState(SAMPWALL_EELECNEW ,p,q,iSurfSide)) &
                                            / (SurfSideArea(p,q,iSurfSide) * TimeSampleTemp)
      END IF

      ! Number of simulation particle impacts per iteration
      MacroSurfaceVal(5,p,q,OutputCounter) = CounterSum / IterNum

      ! Boundary index
      MacroSurfaceVal(6,p,q,OutputCounter) = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,GlobalSideID))

      ! --- Output of optional variables
      nVarCount = MACROSURF_NVARS
      ! Output of the pump capacity
      IF(nPorousBC.GT.0) THEN
        DO iPBC=1, nPorousBC
          IF(MapSurfSideToPorousSide_Shared(iSurfSide).EQ.0) CYCLE
          IF(PorousBCInfo_Shared(1,MapSurfSideToPorousSide_Shared(iSurfSide)).EQ.iPBC) THEN
            nVarCount = nVarCount + 1
            ! Pump capacity is already in cubic meter per second (diving by the number of iterations)
            MacroSurfaceVal(nVarCount,p,q,OutputCounter) = SampWallPumpCapacity(iSurfSide) / IterNum
          END IF
        END DO
      END IF
      ! Output of the wall temperature
      IF (PartBound%OutputWallTemp) THEN
        IF ((MacroSurfaceVal(4,p,q,OutputCounter).GT.0.0).AND.PartBound%AdaptWallTemp) THEN
          iBC = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,GlobalSideID))
          BoundaryWallTemp(p,q,iSurfSide) = (MacroSurfaceVal(4,p,q,OutputCounter) &
              /(StefanBoltzmannConst*PartBound%RadiativeEmissivity(iBC)))**(1./4.)
        END IF
        nVarCount = nVarCount + 1
        MacroSurfaceVal(nVarCount,p,q,OutputCounter) = BoundaryWallTemp(p,q,iSurfSide)
      END IF
      ! Output of the sticking probability
      IF(ANY(PartBound%SurfaceModel.EQ.1)) THEN
        nVarCount = nVarCount + 1
        IF(nVarCount.GT.SurfOutputSize) CALL Abort(__STAMP__,'ERROR in CalcSurfaceValues: Number of output variables greater than the allocated array!')
        IF(CounterSum.GT.0) MacroSurfaceVal(nVarCount,p,q,OutputCounter) = SampWallState(SWIStickingCoefficient,p,q,iSurfSide) / CounterSum
      END IF
      ! --- Species-specific output (in a separate array)
      DO iSpec=1,nSpecies
        idx = 1
        ! Species-specific counter of simulation particle impacts per iteration
        MacroSurfaceSpecVal(idx,p,q,OutputCounter,iSpec) = SampWallState(SAMPWALL_NVARS+iSpec,p,q,iSurfSide) / IterNum
        ! Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z) and angle
        IF(CalcSurfaceImpact)THEN
          nImpacts = SampWallImpactNumber(iSpec,p,q,iSurfSide)
          IF(nImpacts.GT.0.)THEN
            ! Add average impact energy for each species (trans, rot, vib)
            idx = idx + 1
            MacroSurfaceSpecVal(idx,p,q,OutputCounter,iSpec) = SampWallImpactEnergy(iSpec,1,p,q,iSurfSide) / nImpacts
            idx = idx + 1
            MacroSurfaceSpecVal(idx,p,q,OutputCounter,iSpec) = SampWallImpactEnergy(iSpec,2,p,q,iSurfSide) / nImpacts
            idx = idx + 1
            MacroSurfaceSpecVal(idx,p,q,OutputCounter,iSpec) = SampWallImpactEnergy(iSpec,3,p,q,iSurfSide) / nImpacts
            idx = idx + 1
            MacroSurfaceSpecVal(idx,p,q,OutputCounter,iSpec) = SampWallImpactEnergy(iSpec,4,p,q,iSurfSide) / nImpacts

            ! Add average impact vector (x,y,z) for each species
            idx = idx + 1
            MacroSurfaceSpecVal(idx,p,q,OutputCounter,iSpec) = SampWallImpactVector(iSpec,1,p,q,iSurfSide) / nImpacts
            idx = idx + 1
            MacroSurfaceSpecVal(idx,p,q,OutputCounter,iSpec) = SampWallImpactVector(iSpec,2,p,q,iSurfSide) / nImpacts
            idx = idx + 1
            MacroSurfaceSpecVal(idx,p,q,OutputCounter,iSpec) = SampWallImpactVector(iSpec,3,p,q,iSurfSide) / nImpacts

            ! Add average impact angle for each species
            idx = idx + 1
            MacroSurfaceSpecVal(idx,p,q,OutputCounter,iSpec) = SampWallImpactAngle(iSpec,p,q,iSurfSide) / nImpacts

            ! Add number of impacts (real particles with weighting factor)
            idx = idx + 1
            MacroSurfaceSpecVal(idx,p,q,OutputCounter,iSpec) = nImpacts

            ! Add number of impacts per second
            idx = idx + 1
            IF(VarTimeStep%UseSpeciesSpecific) THEN
              TimeSampleTemp = TimeSample * Species(iSpec)%TimeStepFactor
            ELSE
              TimeSampleTemp = TimeSample
            END IF
            MacroSurfaceSpecVal(idx,p,q,OutputCounter,iSpec) = nImpacts / (SurfSideArea(p,q,iSurfSide) * TimeSampleTemp)
          END IF ! nImpacts.GT.0.
        END IF ! CalcSurfaceImpact
      END DO ! iSpec=1,nSpecies
    END DO ! q=1,nSurfSample
  END DO ! p=1,nSurfSample
END DO ! iSurfSide=1,nComputeNodeSurfSides

#if USE_MPI
END ASSOCIATE
IF (PartBound%OutputWallTemp) THEN
  CALL BARRIER_AND_SYNC(BoundaryWallTemp_Shared_Win,MPI_COMM_SHARED)
END IF
#endif /*USE_MPI*/

CALL WriteSurfSampleToHDF5(TRIM(MeshFile),ActualTime)

MacroSurfaceVal = 0.
MacroSurfaceSpecVal = 0.

END SUBROUTINE CalcSurfaceValues


SUBROUTINE WriteSurfSampleToHDF5(MeshFileName,OutputTime)
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
USE MOD_SurfaceModel_Vars       ,ONLY: nPorousBC
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfSample,CalcSurfaceImpact
USE MOD_Particle_Boundary_Vars  ,ONLY: nGlobalOutputSides
USE MOD_Particle_boundary_Vars  ,ONLY: nComputeNodeSurfOutputSides,offsetComputeNodeSurfOutputSide
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfBC,SurfBCName, PartBound
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfOutputSize,SurfSpecOutputSize
USE MOD_Particle_Boundary_Vars  ,ONLY: MacroSurfaceVal,MacroSurfaceSpecVal
USE MOD_Particle_Vars           ,ONLY: nSpecies
#if USE_MPI
USE MOD_Particle_Boundary_Vars  ,ONLY: nGlobalSurfSides
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
INTEGER                             :: nVar2D_Total, nVarCount, iSpec, iPBC
REAL                                :: StartT,EndT
!===================================================================================================================================

#if USE_MPI
! Return if not a sampling leader
IF (MPI_COMM_LEADERS_SURF.EQ.MPI_COMM_NULL) RETURN
CALL MPI_BARRIER(MPI_COMM_LEADERS_SURF,iERROR)

! Return if no sampling sides
IF (nGlobalSurfSides      .EQ.0) RETURN
#endif /*USE_MPI*/

IF (mySurfRank.EQ.0) THEN
  WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' WRITE DSMCSurfSTATE TO HDF5 FILE...'
  StartT=LOCALTIME()
END IF

FileName   = TIMESTAMP(TRIM(ProjectName)//'_DSMCSurfState',OutputTime)
FileString = TRIM(FileName)//'.h5'

nVar2D_Total = SurfOutputSize + SurfSpecOutputSize*nSpecies

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

  ! fill varnames for total values (add new variables to MACROSURF_NVARS)
  CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Total_ForcePerAreaX')
  CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Total_ForcePerAreaY')
  CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Total_ForcePerAreaZ')
  CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Total_HeatFlux')
  CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Total_SimPartPerIter')
  CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'iBC')

  IF(nPorousBC.GT.0) THEN
    DO iPBC = 1, nPorousBC
      WRITE(PBCID,'(I2.2)') iPBC
      CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'PorousBC'//TRIM(PBCID)//'_PumpCapacity')
    END DO
  END IF

  IF (PartBound%OutputWallTemp) CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Wall_Temperature')
  IF (ANY(PartBound%SurfaceModel.EQ.1)) CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Sticking_Coefficient')

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
! WARNING: Only the sampling leaders write the data to .h5
ASSOCIATE (&
      nVar2D_Total         => INT(nVar2D_Total,IK)                    , &
      nSurfSample          => INT(nSurfSample,IK)                     , &
      nGlobalSides         => INT(nGlobalOutputSides,IK)                    , &
      nLocalSides          => INT(nComputeNodeSurfOutputSides,IK)     , &
      offsetSurfSide       => INT(offsetComputeNodeSurfOutputSide,IK) , &
      SurfOutputSize       => INT(SurfOutputSize,IK)                  , &
      SurfSpecOutputSize   => INT(SurfSpecOutputSize,IK))
  DO iSpec = 1,nSpecies
    CALL WriteArrayToHDF5(DataSetName=H5_Name             , rank=4                                           , &
                            nValGlobal =(/nVar2D_Total      , nSurfSample , nSurfSample , nGlobalSides   /)  , &
                            nVal       =(/SurfSpecOutputSize       , nSurfSample , nSurfSample , nLocalSides/)      , &
                            offset     =(/INT(nVarCount,IK) , 0_IK        , 0_IK        , offsetSurfSide/)   , &
                            collective =.FALSE.                                                              , &
                            RealArray  = MacroSurfaceSpecVal(1:SurfSpecOutputSize,1:nSurfSample,1:nSurfSample,1:nLocalSides,iSpec))
    nVarCount = nVarCount + INT(SurfSpecOutputSize)
  END DO
  CALL WriteArrayToHDF5(DataSetName=H5_Name            , rank=4                                              , &
                        nValGlobal =(/nVar2D_Total     , nSurfSample, nSurfSample , nGlobalSides/)           , &
                        nVal       =(/SurfOutputSize           , nSurfSample, nSurfSample , nLocalSides/)            , &
                        offset     =(/INT(nVarCount,IK), 0_IK       , 0_IK        , offsetSurfSide/)         , &
                        collective =.FALSE.                                                                  , &
                        RealArray  = MacroSurfaceVal(1:SurfOutputSize,1:nSurfSample,1:nSurfSample,1:nLocalSides))
END ASSOCIATE

CALL CloseDataFile()

IF (mySurfRank.EQ.0) THEN
  EndT=LOCALTIME()
  CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE., rank=mySurfRank)
END IF

END SUBROUTINE WriteSurfSampleToHDF5


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
!USE MOD_DSMC_Vars                      ,ONLY: DSMC
USE MOD_Particle_Boundary_Vars
!USE MOD_Particle_Vars                  ,ONLY: WriteMacroSurfaceValues
USE MOD_SurfaceModel_Vars              ,ONLY: nPorousBC
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
!IF (.NOT.WriteMacroSurfaceValues.AND..NOT.DSMC%CalcSurfaceVal.AND..NOT.(ANY(PartBound%Reactive))) RETURN

! Return if no sampling surfaces on node
IF (.NOT.SurfTotalSideOnNode) RETURN

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
SDEALLOCATE(MacroSurfaceVal)
SDEALLOCATE(MacroSurfaceSpecVal)

END SUBROUTINE FinalizeParticleBoundarySampling

END MODULE MOD_Particle_Boundary_Sampling
