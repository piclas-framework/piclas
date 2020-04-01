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

MODULE MOD_Particle_BGM
!===================================================================================================================================
!> Contains
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

PUBLIC::DefineParametersParticleBGM
PUBLIC::BuildBGMAndIdentifyHaloRegion

CONTAINS

!==================================================================================================================================
!> Define parameters for particle backgroundmesh
!==================================================================================================================================
SUBROUTINE DefineParametersParticleBGM()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection('BGM')

! Background mesh init variables
CALL prms%CreateRealArrayOption('Part-FIBGMdeltas'&
  , 'Define the deltas for the cartesian Fast-Init-Background-Mesh.'//&
  ' They should be of the similar size as the smallest cells of the used mesh for simulation.'&
  , '1. , 1. , 1.')
CALL prms%CreateRealArrayOption('Part-FactorFIBGM'&
  , 'Factor with which the background mesh will be scaled.'&
  , '1. , 1. , 1.')

END SUBROUTINE DefineParametersParticleBGM


SUBROUTINE BuildBGMAndIdentifyHaloRegion()
!===================================================================================================================================
!> computes the BGM-indices of an element and maps the number of element and which element to each BGM cell
!> BGM is only saved for compute-node-mesh + halo-region on shared memory
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars              ,ONLY: nElems,offsetElem,nGlobalElems
USE MOD_Partilce_Periodic_BC   ,ONLY: InitPeriodicBC
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO,FIBGM_nElems,FIBGM_Element,FIBGM_offsetElem
USE MOD_Particle_Tracking_Vars ,ONLY: Distance,ListDistance
USE MOD_Equation_Vars          ,ONLY: c
USE MOD_ReadInTools            ,ONLY: GETREAL, GetRealArray, PrintOption
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
USE MOD_TimeDisc_Vars          ,ONLY: iStage,nRKStages,RK_c
#endif
#if ! (USE_HDG)
USE MOD_DG                     ,ONLY: DGTimeDerivative_weakForm
USE MOD_CalcTimeStep           ,ONLY: CalcTimeStep
USE MOD_TimeDisc_Vars          ,ONLY: time
#endif /*USE_HDG*/
#if USE_MPI
USE MOD_MPI_Shared_Vars
USE MOD_MPI_Shared             ,ONLY: Allocate_Shared
USE MOD_PICDepo_Vars           ,ONLY: DepositionType, r_sf
USE MOD_Particle_MPI_Vars      ,ONLY: SafetyFactor,halo_eps_velo,halo_eps,halo_eps2
USE MOD_Particle_Vars          ,ONLY: manualtimestep
#else
USE MOD_Mesh_Vars              ,ONLY: NodeCoords, ElemToBGM_Shared,BoundsOfElem_Shared,NodeCoords_Shared
USE MOD_Mesh_Vars              ,ONLY: ElemInfo_Shared,SideInfo_Shared, NodeInfo_Shared
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem
REAL                           :: xmin, xmax, ymin, ymax, zmin, zmax
INTEGER                        :: iBGM, jBGM, kBGM
INTEGER                        :: BGMimax, BGMimin, BGMjmax, BGMjmin, BGMkmax, BGMkmin
INTEGER                        :: BGMiDelta,BGMjDelta,BGMkDelta
INTEGER                        :: BGMCellXmax, BGMCellXmin, BGMCellYmax, BGMCellYmin, BGMCellZmax, BGMCellZmin
#if USE_MPI
INTEGER                        :: iSide, SideID
INTEGER                        :: ElemID
REAL                           :: deltaT
REAL                           :: globalDiag
INTEGER,ALLOCATABLE            :: sendbuf(:,:,:), recvbuf(:,:,:)
INTEGER,ALLOCATABLE            :: offsetElemsInBGMCell(:,:,:)
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
INTEGER                        :: nHaloElems, nMPISidesShared, currentOffset, moveBGMindex
INTEGER,ALLOCATABLE            :: offsetCNHalo2GlobalElem(:), offsetMPISideShared(:)
REAL,ALLOCATABLE               :: BoundsOfElemCenter(:), MPISideBoundsOfElemCenter(:,:)
LOGICAL                        :: ElemInsideHalo
INTEGER                        :: FirstElem, LastElem, firstHaloElem, lastHaloElem
INTEGER                        :: offsetNodeID, nNodeIDs, firstNodeID, lastNodeID
#else
INTEGER,ALLOCATABLE            :: ElemToBGM(:,:)
REAL,POINTER                   :: NodeCoordsPointer(:,:,:,:,:)
#endif
!===================================================================================================================================

! Read parameter for FastInitBackgroundMesh (FIBGM)
GEO%FIBGMdeltas(1:3) = GETREALARRAY('Part-FIBGMdeltas',3,'1. , 1. , 1.')
GEO%FactorFIBGM(1:3) = GETREALARRAY('Part-FactorFIBGM',3,'1. , 1. , 1.')
GEO%FIBGMdeltas(1:3) = 1./GEO%FactorFIBGM(1:3) * GEO%FIBGMdeltas(1:3)

#if USE_MPI
MPISharedSize = INT(6*nGlobalElems,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/6,nGlobalElems/),ElemToBGM_Shared_Win,ElemToBGM_Shared)
CALL Allocate_Shared(MPISharedSize,(/2,3,nGlobalElems/),BoundsOfElem_Shared_Win,BoundsOfElem_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemToBGM_Shared_Win,IERROR)
CALL MPI_WIN_LOCK_ALL(0,BoundsOfElem_Shared_Win,IERROR)
firstElem = INT(REAL(myComputeNodeRank*nGlobalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nGlobalElems)/REAL(nComputeNodeProcessors))
#else
! In order to use only one type of variables VarName_Shared in code structure such as tracking etc. for NON_MPI
! the same variables are allocated on the single proc and used from mesh_vars instead of mpi_shared_vars
ALLOCATE(ElemToBGM_Shared(1:6,1:nElems))
ALLOCATE(BoundsOfElem_Shared(1:2,1:3,1:nElems)) ! 1-2: Min, Max value; 1-3: x,y,z
firstElem = 1
lastElem  = nElems
#endif  /*USE_MPI*/

moveBGMindex = 1 ! BGM indeces must be >1 --> move by 1
DO iElem = firstElem, lastElem
  offSetNodeID=ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem)
  nNodeIDs=ElemInfo_Shared(ELEM_LASTNODEIND,iElem)-ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem)
  firstNodeID = offsetNodeID+1
  lastNodeID = offsetNodeID+nNodeIDs

  xmin=MINVAL(NodeCoords_Shared(1,firstNodeID:lastNodeID))
  xmax=MAXVAL(NodeCoords_Shared(1,firstNodeID:lastNodeID))
  ymin=MINVAL(NodeCoords_Shared(2,firstNodeID:lastNodeID))
  ymax=MAXVAL(NodeCoords_Shared(2,firstNodeID:lastNodeID))
  zmin=MINVAL(NodeCoords_Shared(3,firstNodeID:lastNodeID))
  zmax=MAXVAL(NodeCoords_Shared(3,firstNodeID:lastNodeID))

  BoundsOfElem_Shared(1,1,iElem) = xmin
  BoundsOfElem_Shared(2,1,iElem) = xmax
  BoundsOfElem_Shared(1,2,iElem) = ymin
  BoundsOfElem_Shared(2,2,iElem) = ymax
  BoundsOfElem_Shared(1,3,iElem) = zmin
  BoundsOfElem_Shared(2,3,iElem) = zmax

  ! BGM indeces must be >0 --> move by 1
  ElemToBGM_Shared(1,iElem) = FLOOR((xmin-GEO%xminglob)/GEO%FIBGMdeltas(1)) + moveBGMindex
  ElemToBGM_Shared(2,iElem) = FLOOR((xmax-GEO%xminglob)/GEO%FIBGMdeltas(1)) + moveBGMindex
  ElemToBGM_Shared(3,iElem) = FLOOR((ymin-GEO%yminglob)/GEO%FIBGMdeltas(2)) + moveBGMindex
  ElemToBGM_Shared(4,iElem) = FLOOR((ymax-GEO%yminglob)/GEO%FIBGMdeltas(2)) + moveBGMindex
  ElemToBGM_Shared(5,iElem) = FLOOR((zmin-GEO%zminglob)/GEO%FIBGMdeltas(3)) + moveBGMindex
  ElemToBGM_Shared(6,iElem) = FLOOR((zmax-GEO%zminglob)/GEO%FIBGMdeltas(3)) + moveBGMindex
END DO ! iElem = firstElem, lastElem

#if USE_MPI
CALL MPI_WIN_SYNC(ElemToBGM_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(BoundsOfElem_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
#endif  /*USE_MPI*/

!CALL InitPeriodicBC()

! deallocate stuff // required for dynamic load balance
#if USE_LOADBALANCE
IF (ALLOCATED(GEO%FIBGM)) THEN
  DO iBGM=GEO%FIBGMimin,GEO%FIBGMimax
    DO jBGM=GEO%FIBGMjmin,GEO%FIBGMjmax
      DO kBGM=GEO%FIBGMkmin,GEO%FIBGMkmax
        SDEALLOCATE(GEO%FIBGM(iBGM,jBGM,kBGM)%Element)
        !SDEALLOCATE(GEO%FIBGM(iBGM,jBGM,kBGM)%ShapeProcs)
        !SDEALLOCATE(GEO%FIBGM(iBGM,jBGM,kBGM)%PaddingProcs)
        !SDEALLOCATE(GEO%FIBGM(i,k,l)%SharedProcs)
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
  DEALLOCATE(GEO%FIBGM)
END IF
#endif /*USE_LOADBALANCE*/

SafetyFactor  =GETREAL('Part-SafetyFactor','1.0')
halo_eps_velo =GETREAL('Particles-HaloEpsVelo','0')

IF (nComputeNodeProcessors.EQ.nProcessors_Global) THEN
  halo_eps  = 0.
  halo_eps2 = 0.
ELSE
  IF (ManualTimeStep.EQ.0.0) THEN
#if !(USE_HDG)
    CALL DGTimeDerivative_weakForm(time,time,0,doSource=.TRUE.)
    deltaT = CalcTimeStep()
#else
     CALL abort(&
  __STAMP__&
  , 'ManualTimeStep.EQ.0.0 -> ManualTimeStep is not defined correctly! Particles-ManualTimeStep = ',RealInfoOpt=ManualTimeStep)
#endif /*USE_HDG*/
  ELSE
    deltaT=ManualTimeStep
  END IF
  IF (halo_eps_velo.EQ.0) halo_eps_velo = c
#if (PP_TimeDiscMethod==4 || PP_TimeDiscMethod==200 || PP_TimeDiscMethod==42 || PP_TimeDiscMethod==43)
  IF (halo_eps_velo.EQ.c) THEN
     CALL abort(&
  __STAMP__&
  , 'halo_eps_velo.EQ.c -> Halo Eps Velocity for MPI not defined')
  END IF
#endif
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
  halo_eps = RK_c(2)
  DO iStage=2,nRKStages-1
    halo_eps = MAX(halo_eps,RK_c(iStage+1)-RK_c(iStage))
  END DO
  halo_eps = MAX(halo_eps,1.-RK_c(nRKStages))
  CALL PrintOption('max. RKdtFrac','CALCUL.',RealOpt=halo_eps)
  halo_eps = halo_eps*halo_eps_velo*deltaT*SafetyFactor !dt multiplied with maximum RKdtFrac
#else
  halo_eps = halo_eps_velo*deltaT*SafetyFactor ! for RK too large
#endif

  ! Check whether halo_eps is smaller than shape function radius
  ! e.g. 'shape_function', 'shape_function_1d', 'shape_function_cylindrical', 'shape_function_spherical', 'shape_function_simple'
  IF(TRIM(DepositionType(1:MIN(14,LEN(TRIM(ADJUSTL(DepositionType)))))).EQ.'shape_function')THEN
    IF(halo_eps.LT.r_sf)THEN
      SWRITE(UNIT_stdOut,'(A)') ' halo_eps is smaller than shape function radius. Setting halo_eps=r_sf'
      halo_eps = halo_eps + r_sf
      CALL PrintOption('max. RKdtFrac','CALCUL.',RealOpt=halo_eps)
    END IF
  END IF

  ! limit halo_eps to diagonal of bounding box
  globalDiag = SQRT( (GEO%xmaxglob-GEO%xminglob)**2 &
                   + (GEO%ymaxglob-GEO%yminglob)**2 &
                   + (GEO%zmaxglob-GEO%zminglob)**2 )
  IF(halo_eps.GT.globalDiag)THEN
    CALL PrintOption('unlimited halo distance','CALCUL.',RealOpt=halo_eps)
    SWRITE(UNIT_stdOut,'(A38)') ' |   limitation of halo distance  |    '
    halo_eps=globalDiag
  END IF

  halo_eps2=halo_eps*halo_eps
  CALL PrintOption('halo distance','CALCUL.',RealOpt=halo_eps)
END IF

moveBGMindex = 1 ! BGM indices must be >0 --> move by 1
! enlarge BGM with halo region (all element outside of this region will be cut off)
BGMimin = FLOOR((MAX(GEO%CNxmin-halo_eps,GEO%xminglob)-GEO%xminglob)/GEO%FIBGMdeltas(1)) + moveBGMindex
BGMimax = FLOOR((MIN(GEO%CNxmax+halo_eps,GEO%xmaxglob)-GEO%xminglob)/GEO%FIBGMdeltas(1)) + moveBGMindex
BGMjmin = FLOOR((MAX(GEO%CNymin-halo_eps,GEO%yminglob)-GEO%yminglob)/GEO%FIBGMdeltas(2)) + moveBGMindex
BGMjmax = FLOOR((MIN(GEO%CNymax+halo_eps,GEO%ymaxglob)-GEO%yminglob)/GEO%FIBGMdeltas(2)) + moveBGMindex
BGMkmin = FLOOR((MAX(GEO%CNzmin-halo_eps,GEO%zminglob)-GEO%zminglob)/GEO%FIBGMdeltas(3)) + moveBGMindex
BGMkmax = FLOOR((MIN(GEO%CNzmax+halo_eps,GEO%zmaxglob)-GEO%zminglob)/GEO%FIBGMdeltas(3)) + moveBGMindex
! write function-local BGM indices into global variables
GEO%FIBGMimin=BGMimin
GEO%FIBGMimax=BGMimax
GEO%FIBGMjmin=BGMjmin
GEO%FIBGMjmax=BGMjmax
GEO%FIBGMkmin=BGMkmin
GEO%FIBGMkmax=BGMkmax
! initialize BGM min/max indices using GEO min/max distances
!GEO%FIBGMimax = INT((GEO%xmax-GEO%xminglob)/GEO%FIBGMdeltas(1))+1  + moveBGMindex
!GEO%FIBGMimin = INT((GEO%xmin-GEO%xminglob)/GEO%FIBGMdeltas(1))-1  + moveBGMindex
!GEO%FIBGMjmax = INT((GEO%ymax-GEO%yminglob)/GEO%FIBGMdeltas(2))+1  + moveBGMindex
!GEO%FIBGMjmin = INT((GEO%ymin-GEO%yminglob)/GEO%FIBGMdeltas(2))-1  + moveBGMindex
!GEO%FIBGMkmax = INT((GEO%zmax-GEO%zminglob)/GEO%FIBGMdeltas(3))+1  + moveBGMindex
!GEO%FIBGMkmin = INT((GEO%zmin-GEO%zminglob)/GEO%FIBGMdeltas(3))-1  + moveBGMindex

ALLOCATE(GEO%FIBGM(BGMimin:BGMimax,BGMjmin:BGMjmax,BGMkmin:BGMkmax))

! null number of element per BGM cell
DO kBGM = BGMkmin,BGMkmax
  DO jBGM = BGMjmin,BGMjmax
    DO iBGM = BGMimin,BGMimax
      GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = 0
    END DO ! kBGM
  END DO ! jBGM
END DO ! iBGM

#if USE_MPI
! check which element is inside of compute-node domain (1),
! check which element is inside of compute-node halo (2)
! and which element is outside of compute-node domain (0)
! first do coarse check with BGM
IF (nComputeNodeProcessors.EQ.nProcessors_Global) THEN
  ElemInfo_Shared(ELEM_HALOFLAG,firstElem:lastElem)=1
ELSE
  ElemInfo_Shared(ELEM_HALOFLAG,firstElem:lastElem)=0
  DO iElem = firstElem, lastElem
    BGMCellXmin = ElemToBGM_Shared(1,iElem)
    BGMCellXmax = ElemToBGM_Shared(2,iElem)
    BGMCellYmin = ElemToBGM_Shared(3,iElem)
    BGMCellYmax = ElemToBGM_Shared(4,iElem)
    BGMCellZmin = ElemToBGM_Shared(5,iElem)
    BGMCellZmax = ElemToBGM_Shared(6,iElem)
    ! add current element to number of BGM-elems
    ! ATTENTION: THIS ONLY ADDS THE ELEMENT TO THE BGM CELLS ON THE NODE WHILE
    ! SKIPPING BGM CELLS OUTSIDE. WE END UP WITH PARTIALLY ADDED ELEMENTS
    DO iBGM = BGMCellXmin,BGMCellXmax
      IF(iBGM.LT.BGMimin) CYCLE
      IF(iBGM.GT.BGMimax) CYCLE
      DO jBGM = BGMCellYmin,BGMCellYmax
        IF(jBGM.LT.BGMjmin) CYCLE
        IF(jBGM.GT.BGMjmax) CYCLE
        DO kBGM = BGMCellZmin,BGMCellZmax
          IF(kBGM.LT.BGMkmin) CYCLE
          IF(kBGM.GT.BGMkmax) CYCLE
          !GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1
          IF(iElem.GE.offsetComputeNodeElem+1 .AND. iElem.LE.offsetComputeNodeElem+nComputeNodeElems) THEN
            ElemInfo_Shared(ELEM_HALOFLAG,iElem)=1 ! compute-node element
          ELSE
            ElemInfo_Shared(ELEM_HALOFLAG,iElem)=2 ! halo element
          END IF
        END DO ! kBGM
      END DO ! jBGM
    END DO ! iBGM
  END DO ! iElem
  CALL MPI_WIN_SYNC(ElemInfo_Shared_Win,IERROR)
  CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

  ! sum up potential halo elements and create correct offset mapping via ElemInfo_Shared
  nHaloElems = 0
  DO iElem = 1, nGlobalElems
    IF (ElemInfo_Shared(ELEM_HALOFLAG,iElem).EQ.2) THEN
      nHaloElems = nHaloElems + 1
    END IF
  END DO
  ALLOCATE(offsetCNHalo2GlobalElem(1:nHaloElems))
  nHaloElems = 0
  DO iElem = 1, nGlobalElems
    IF (ElemInfo_Shared(ELEM_HALOFLAG,iElem).EQ.2) THEN
      nHaloElems = nHaloElems + 1
      offsetCNHalo2GlobalElem(nHaloElems) = iElem
    END IF
  END DO

  ! sum all MPI-side of compute-node and create correct offset mapping in SideInfo_Shared
  nMPISidesShared = 0
  DO iSide = 1, nNonUniqueGlobalSides
    IF (SideInfo_Shared(SIDEINFOSIZE+1,iSide).EQ.2) THEN
      nMPISidesShared = nMPISidesShared + 1
    END IF
  END DO
  ALLOCATE(offsetMPISideShared(nMPISidesShared))

  nMPISidesShared = 0
  DO iSide = 1, nNonUniqueGlobalSides
    IF (SideInfo_Shared(SIDEINFOSIZE+1,iSide).EQ.2) THEN
      nMPISidesShared = nMPISidesShared + 1
      offsetMPISideShared(nMPISidesShared) = iSide
    END IF
  END DO

  ! Distribute nHaloElements evenly on compute-node procs
  IF (nHaloElems.GT.nComputeNodeProcessors) THEN
    firstHaloElem = INT(REAL( myComputeNodeRank   *nHaloElems)/REAL(nComputeNodeProcessors))+1
    lastHaloElem  = INT(REAL((myComputeNodeRank+1)*nHaloElems)/REAL(nComputeNodeProcessors))
  ELSE
    firstHaloElem = myComputeNodeRank + 1
    IF (myComputeNodeRank.LT.nHaloElems) THEN
      lastHaloElem = myComputeNodeRank + 1
    ELSE
      lastHaloElem = 0
    END IF
  END IF

  ALLOCATE(MPISideBoundsOfElemCenter(1:4,1:nMPISidesShared))
  DO iSide = 1, nMPISidesShared
    SideID = offsetMPISideShared(iSide)
    ElemID = SideInfo_Shared(SIDE_ELEMID,SideID)
    MPISideBoundsOfElemCenter(1:3,iSide) = (/ SUM(BoundsOfElem_Shared(1:2,1,ElemID)), &
                                              SUM(BoundsOfElem_Shared(1:2,2,ElemID)), &
                                              SUM(BoundsOfElem_Shared(1:2,3,ElemID)) /) / 2.
    MPISideBoundsOfElemCenter(4,iSide) = VECNORM ((/ BoundsOfElem_Shared(2,1,ElemID)-BoundsOfElem_Shared(1,1,ElemID), &
                                                      BoundsOfElem_Shared(2,2,ElemID)-BoundsOfElem_Shared(1,2,ElemID), &
                                                      BoundsOfElem_Shared(2,3,ElemID)-BoundsOfElem_Shared(1,3,ElemID) /) / 2.)
  END DO

  ! do refined check: (refined halo region reduction)
  ! check the bounding box of each element in compute-nodes' halo domain
  ! against the bounding boxes of the elements of the MPI-surface (inter compute-node MPI sides)
  ALLOCATE(BoundsOfElemCenter(1:4))
  DO iElem = firstHaloElem, lastHaloElem
    ElemID = offsetCNHalo2GlobalElem(iElem)
    ElemInsideHalo = .FALSE.
    BoundsOfElemCenter(1:3) = (/ SUM(BoundsOfElem_Shared(1:2,1,ElemID)), &
                                 SUM(BoundsOfElem_Shared(1:2,2,ElemID)), &
                                 SUM(BoundsOfElem_Shared(1:2,3,ElemID)) /) / 2.
    BoundsOfElemCenter(4) = VECNORM ((/ BoundsOfElem_Shared(2,1,ElemID)-BoundsOfElem_Shared(1,1,ElemID), &
                                        BoundsOfElem_Shared(2,2,ElemID)-BoundsOfElem_Shared(1,2,ElemID), &
                                        BoundsOfElem_Shared(2,3,ElemID)-BoundsOfElem_Shared(1,3,ElemID) /) / 2.)
    DO iSide = 1, nMPISidesShared
      ! compare distance of centers with sum of element outer radii+halo_eps
      IF (VECNORM(BoundsOfElemCenter(1:3)-MPISideBoundsOfElemCenter(1:3,iSide)) &
          .GT. halo_eps+BoundsOfElemCenter(4)+MPISideBoundsOfElemCenter(4,iSide) ) CYCLE
      ElemInsideHalo = .TRUE.
      EXIT
    END DO ! iSide = 1, nMPISidesShared
    IF (.NOT.ElemInsideHalo) THEN
      ElemInfo_Shared(ELEM_HALOFLAG,ElemID)=0
    ELSE
      ! Only add element to BGM if inside halo region on node.
      ! THIS IS WRONG. WE ARE WORKING ON THE CN HALO REGION. IF WE OMIT THE
      ! ELEMENT HERE, WE LOOSE IT. IF WE KEEP IT, WE BREAK AT 589. YOUR CALL.
!      print *, ElemID,BoundsOfElem_Shared(1,1,iElem),GEO%xminglob
      BGMCellXmin = MAX(ElemToBGM_Shared(1,ElemID),BGMimin)
      BGMCellXmax = MIN(ElemToBGM_Shared(2,ElemID),BGMimax)
      BGMCellYmin = MAX(ElemToBGM_Shared(3,ElemID),BGMjmin)
      BGMCellYmax = MIN(ElemToBGM_Shared(4,ElemID),BGMjmax)
      BGMCellZmin = MAX(ElemToBGM_Shared(5,ElemID),BGMkmin)
      BGMCellZmax = MIN(ElemToBGM_Shared(6,ElemID),BGMkmax)
      ! add current element to number of BGM-elems
      DO iBGM = BGMCellXmin,BGMCellXmax
        DO jBGM = BGMCellYmin,BGMCellYmax
          DO kBGM = BGMCellZmin,BGMCellZmax
            GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1
          END DO ! kBGM
        END DO ! jBGM
      END DO ! iBGM
    END IF
  END DO ! iElem = firstHaloElem, lastHaloElem
END IF ! nComputeNodeProcessors.EQ.nProcessors_Global
CALL MPI_WIN_SYNC(ElemInfo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif  /*USE_MPI*/

!--- compute number of elements in each background cell
DO iElem = offsetElem+1, offsetElem+nElems
  BGMCellXmin = ElemToBGM_Shared(1,iElem)
  BGMCellXmax = ElemToBGM_Shared(2,iElem)
  BGMCellYmin = ElemToBGM_Shared(3,iElem)
  BGMCellYmax = ElemToBGM_Shared(4,iElem)
  BGMCellZmin = ElemToBGM_Shared(5,iElem)
  BGMCellZmax = ElemToBGM_Shared(6,iElem)
  ! add current element to number of BGM-elems
  DO iBGM = BGMCellXmin,BGMCellXmax
    DO jBGM = BGMCellYmin,BGMCellYmax
      DO kBGM = BGMCellZmin,BGMCellZmax
        GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
END DO ! iElem

! alternative nElem count with cycles
!DO iElem = firstElem, lastElem
!  IF (ElemInfo_Shared(ELEM_HALOFLAG,iElem).EQ.0) CYCLE
!  BGMCellXmin = ElemToBGM_Shared(1,iElem)
!  BGMCellXmax = ElemToBGM_Shared(2,iElem)
!  BGMCellYmin = ElemToBGM_Shared(3,iElem)
!  BGMCellYmax = ElemToBGM_Shared(4,iElem)
!  BGMCellZmin = ElemToBGM_Shared(5,iElem)
!  BGMCellZmax = ElemToBGM_Shared(6,iElem)
!  ! add current element to number of BGM-elems
!  DO iBGM = BGMCellXmin,BGMCellXmax
!    DO jBGM = BGMCellYmin,BGMCellYmax
!      DO kBGM = BGMCellZmin,BGMCellZmax
!        GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1
!      END DO ! kBGM
!    END DO ! jBGM
!  END DO ! iBGM
!END DO ! iElem

#if USE_MPI
ALLOCATE(sendbuf(BGMimin:BGMimax,BGMjmin:BGMjmax,BGMkmin:BGMkmax))
ALLOCATE(recvbuf(BGMimin:BGMimax,BGMjmin:BGMjmax,BGMkmin:BGMkmax))
! find max nelems and offset in each BGM cell
DO iBGM = BGMimin,BGMimax
  DO jBGM = BGMjmin,BGMjmax
    DO kBGM = BGMkmin,BGMkmax
      sendbuf(iBGM,jBGM,kBGM)=GEO%FIBGM(iBGM,jBGM,kBGM)%nElem
      recvbuf(iBGM,jBGM,kBGM)=0
    END DO ! kBGM
  END DO ! jBGM
END DO ! iBGM

ALLOCATE(offsetElemsInBGMCell(BGMimin:BGMimax,BGMjmin:BGMjmax,BGMkmin:BGMkmax))
CALL MPI_EXSCAN(sendbuf(:,:,:),recvbuf(:,:,:),((BGMimax-BGMimin)+1)*((BGMjmax-BGMjmin)+1)*((BGMkmax-BGMkmin)+1) &
                ,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,iError)
offsetElemsInBGMCell=recvbuf
DEALLOCATE(recvbuf)

! last proc of compute-node calculates total number of elements in each BGM-cell
! after this loop sendbuf of last proc contains nElems per BGM cell
IF(myComputeNodeRank.EQ.nComputeNodeProcessors-1)THEN
  DO iBGM = BGMimin,BGMimax
    DO jBGM = BGMjmin,BGMjmax
      DO kBGM = BGMkmin,BGMkmax
        sendbuf(iBGM,jBGM,kBGM)=offsetElemsInBGMCell(iBGM,jBGM,kBGM)+GEO%FIBGM(iBGM,jBGM,kBGM)%nElem
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
END IF

! allocated shared memory for nElems per BGM cell
! MPI shared memory is continouos, beginning from 1. All shared arrays have to
! be shifted from BGM[i]min to 1
BGMiDelta=BGMimax-BGMimin
BGMjDelta=BGMjmax-BGMjmin
BGMkDelta=BGMkmax-BGMkmin
MPISharedSize = INT((BGMiDelta+1)*(BGMjDelta+1)*(BGMkDelta+1),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/BGMiDelta+1,BGMjDelta+1,BGMkDelta+1/) &
                    ,FIBGM_nElems_Shared_Win,FIBGM_nElems_Shared)
CALL MPI_WIN_LOCK_ALL(0,FIBGM_nElems_Shared_Win,IERROR)
FIBGM_nElems => FIBGM_nElems_Shared
! allocated shared memory for BGM cell offset in 1D array of BGM to element mapping
MPISharedSize = INT((BGMiDelta+1)*(BGMjDelta+1)*(BGMkDelta+1),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/BGMiDelta+1,BGMjDelta+1,BGMkDelta+1/) &
                    ,FIBGM_offsetElem_Shared_Win,FIBGM_offsetElem_Shared)
CALL MPI_WIN_LOCK_ALL(0,FIBGM_offsetElem_Shared_Win,IERROR)
FIBGM_offsetElem => FIBGM_offsetElem_Shared

! last proc of compute-node writes into shared memory to make nElems per BGM accessible for every proc
IF(myComputeNodeRank.EQ.nComputeNodeProcessors-1)THEN
  currentOffset = 0
  DO iBGM = 1,BGMiDelta+1
    DO jBGM = 1,BGMjDelta+1
      DO kBGM = 1,BGMkDelta+1
        ! senfbuf and recvbuf have to stay on original position. Shift 1 --> BGMimin
        FIBGM_nElems(iBGM,jBGM,kBGM)     = sendbuf(BGMimin+iBGM-1,BGMjmin+jBGM-1,BGMkmin+kBGM-1)
        FIBGM_offsetElem(iBGM,jBGM,kBGM) = currentOffset
        currentOffset = currentoffset    + sendbuf(BGMimin+iBGM-1,BGMjmin+jBGM-1,BGMkmin+kBGM-1)
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
END IF
DEALLOCATE(sendbuf)
CALL MPI_WIN_SYNC(FIBGM_nElems_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(FIBGM_offsetElem_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#else /*NOT USE_MPI*/
ALLOCATE( FIBGM_nElems    (BGMiDelta+1,BGMjDelta+1,BGMkDelta+1))
ALLOCATE( FIBGM_offsetElem(BGMiDelta+1,BGMjDelta+1,BGMkDelta+1))
currentOffset = 0
DO iBGM = 1,BGMiDelta+1
  DO jBGM = 1,BGMjDelta+1
    DO kBGM = 1,BGMkDelta+1
      FIBGM_nElems(iBGM,jBGM,kBGM)     = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem
      FIBGM_offsetElem(iBGM,jBGM,kBGM) = currentOffset
      currentOffset = currentoffset    + GEO%FIBGM(iBGM,jBGM,kBGM)%nElem
    END DO ! kBGM
  END DO ! jBGM
END DO ! iBGM
#endif  /*USE_MPI*/

#if USE_MPI
! allocate 1D array for mapping of BGM cell to Element indeces
MPISharedSize = INT((FIBGM_offsetElem(BGMiDelta+1,BGMjDelta+1,BGMkDelta+1) + &
                     FIBGM_nElems    (BGMiDelta+1,BGMjDelta+1,BGMkDelta+1))  &
                     ,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/FIBGM_offsetElem(BGMiDelta+1,BGMjDelta+1,BGMkDelta+1)   + &
                                     FIBGM_nElems    (BGMiDelta+1,BGMjDelta+1,BGMkDelta+1)/)   &
                     ,FIBGM_Element_Shared_Win,FIBGM_Element_Shared)
CALL MPI_WIN_LOCK_ALL(0,FIBGM_Element_Shared_Win,IERROR)
FIBGM_Element => FIBGM_Element_Shared
#else
ALLOCATE( FIBGM_Element(1:FIBGM_offsetElem(BGMiDelta+1,BGMjDelta+1,BGMkDelta+1) + &
                          FIBGM_nElems    (BGMiDelta+1,BGMjDelta+1,BGMkDelta+1)))
#endif  /*USE_MPI*/

DO iBGM = BGMimin,BGMimax
  DO jBGM = BGMjmin,BGMjmax
    DO kBGM = BGMkmin,BGMkmax
      GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = 0
    END DO ! kBGM
  END DO ! jBGM
END DO ! iBGM

#if USE_MPI
! We might need to expand the halo BGM region


IF (nComputeNodeProcessors.NE.nProcessors_Global) THEN
  DO iElem = firstHaloElem, lastHaloElem
    ElemID = offsetCNHalo2GlobalElem(iElem)
    IF (ElemInfo_Shared(ELEM_HALOFLAG,ElemID).EQ.0) CYCLE
    BGMCellXmin = ElemToBGM_Shared(1,ElemID)
    BGMCellXmax = ElemToBGM_Shared(2,ElemID)
    BGMCellYmin = ElemToBGM_Shared(3,ElemID)
    BGMCellYmax = ElemToBGM_Shared(4,ElemID)
    BGMCellZmin = ElemToBGM_Shared(5,ElemID)
    BGMCellZmax = ElemToBGM_Shared(6,ElemID)
    print *,'halo ',GEO%FIBGMimin,GEO%FIBGMimax,GEO%FIBGMjmin,GEO%FIBGMjmax,GEO%FIBGMkmin,GEO%FIBGMkmax
    print *,'shape', SHAPE(GEO%FIBGM),'rank',myRank
!    print *,BGMCellXmin,BGMCellXmax,BGMCellYmin,BGMCellYmax,BGMCellZmin,BGMCellZmax
    print *,ElemInfo_Shared(ELEM_HALOFLAG,ElemID)
    ! add current Element to BGM-Elem
    DO kBGM = BGMCellZmin,BGMCellZmax
      DO jBGM = BGMCellYmin,BGMCellYmax
        DO iBGM = BGMCellXmin,BGMCellXmax
          GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1
          FIBGM_Element( FIBGM_offsetElem(iBGM,jBGM,kBGM)            & ! offset of BGM cell in 1D array
                              + offsetElemsInBGMCell(iBGM,jBGM,kBGM) & ! offset of BGM nElems in local proc
                              + GEO%FIBGM(iBGM,jBGM,kBGM)%nElem) = ElemID
        END DO ! kBGM
      END DO ! jBGM
    END DO ! iBGM
  END DO ! iElem = firstHaloElem, lastHaloElem
END IF
#endif  /*USE_MPI*/

DO iElem = offsetElem+1, offsetElem+nElems
  BGMCellXmin = ElemToBGM_Shared(1,iElem)
  BGMCellXmax = ElemToBGM_Shared(2,iElem)
  BGMCellYmin = ElemToBGM_Shared(3,iElem)
  BGMCellYmax = ElemToBGM_Shared(4,iElem)
  BGMCellZmin = ElemToBGM_Shared(5,iElem)
  BGMCellZmax = ElemToBGM_Shared(6,iElem)
  ! add current Element to BGM-Elem
  DO kBGM = BGMCellZmin,BGMCellZmax
    DO jBGM = BGMCellYmin,BGMCellYmax
      DO iBGM = BGMCellXmin,BGMCellXmax
        GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1
        FIBGM_Element( FIBGM_offsetElem(iBGM,jBGM,kBGM) & ! offset of BGM cell in 1D array
#if USE_MPI
                            + offsetElemsInBGMCell(iBGM,jBGM,kBGM)    & ! offset of BGM nElems in local proc
#endif  /*USE_MPI*/
                            + GEO%FIBGM(iBGM,jBGM,kBGM)%nElem         ) = iElem
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
END DO ! iElem

!--- map elements to background cells
! alternative if nElem is counted with cycles
!DO iElem = firstElem, lastElem
!  IF (ElemInfo_Shared(ELEM_HALOFLAG,iElem).EQ.0) CYCLE
!  BGMCellXmin = ElemToBGM_Shared(1,iElem)
!  BGMCellXmax = ElemToBGM_Shared(2,iElem)
!  BGMCellYmin = ElemToBGM_Shared(3,iElem)
!  BGMCellYmax = ElemToBGM_Shared(4,iElem)
!  BGMCellZmin = ElemToBGM_Shared(5,iElem)
!  BGMCellZmax = ElemToBGM_Shared(6,iElem)
!  ! add current Element to BGM-Elem
!  DO kBGM = BGMCellZmin,BGMCellZmax
!    DO jBGM = BGMCellYmin,BGMCellYmax
!      DO iBGM = BGMCellXmin,BGMCellXmax
!        GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1
!        FIBGM_Element( FIBGM_offsetElem(iBGM,jBGM,kBGM) & ! offset of BGM cell in 1D array
!#if USE_MPI
!                       + offsetElemsInBGMCell(iBGM,jBGM,kBGM)    & ! offset of BGM nElems in local proc
!#endif  /*USE_MPI*/
!                       + GEO%FIBGM(iBGM,jBGM,kBGM)%nElem         ) = iElem
!      END DO ! kBGM
!    END DO ! jBGM
!  END DO ! iBGM
!END DO ! iElem
#if USE_MPI
DEALLOCATE(offsetElemsInBGMCell)

CALL MPI_WIN_SYNC(FIBGM_Element_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

! sum up Number of all elements on current compute-node (including halo region)
IF (nComputeNodeProcessors.EQ.nProcessors_Global) THEN
  nComputeNodeTotalElems = nGlobalElems
  nComputeNodeTotalSides = nNonUniqueGlobalSides
  nComputeNodeTotalNodes = nNonUniqueGlobalNodes
ELSE
  nComputeNodeTotalElems = 0
  nComputeNodeTotalSides = 0
  nComputeNodeTotalNodes = 0
  DO iElem = 1, nGlobalElems
    IF (ElemInfo_Shared(ELEM_HALOFLAG,iElem).EQ.2 .OR. ElemInfo_Shared(ELEM_HALOFLAG,iElem).EQ.1) THEN
      nComputeNodeTotalElems = nComputeNodeTotalElems + 1
    END IF
  END DO
  ALLOCATE(CNTotalElem2GlobalElem(1:nComputeNodeTotalElems))
  ALLOCATE(GlobalElem2CNTotalElem(1:nGlobalElems))
  nComputeNodeTotalElems = 0
  nComputeNodeTotalElems = 0
  GlobalElem2CNTotalElem(1:nGlobalElems) = -1
  DO iElem = 1,nGlobalElems
    IF (ElemInfo_Shared(ELEM_HALOFLAG,iElem).EQ.1) THEN
      nComputeNodeTotalElems = nComputeNodeTotalElems + 1
      CNTotalElem2GlobalElem(nComputeNodeTotalElems) = iElem
      GlobalElem2CNTotalElem(iElem) = nComputeNodeTotalElems
      nComputeNodeTotalSides = nComputeNodeTotalElems &
                             + (ElemInfo_Shared(ELEM_LASTSIDEIND,iElem) - ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem))
      nComputeNodeTotalNodes = nComputeNodeTotalNodes &
                             + (ElemInfo_Shared(ELEM_LASTNODEIND,iElem) - ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem))
    END IF
  END DO
  DO iElem = 1,nGlobalElems
    IF (ElemInfo_Shared(ELEM_HALOFLAG,iElem).EQ.2) THEN
      nComputeNodeTotalElems = nComputeNodeTotalElems + 1
      CNTotalElem2GlobalElem(nComputeNodeTotalElems) = iElem
      GlobalElem2CNTotalElem(iElem) = nComputeNodeTotalElems
      nComputeNodeTotalSides = nComputeNodeTotalElems &
                             + (ElemInfo_Shared(ELEM_LASTSIDEIND,iElem) - ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem))
      nComputeNodeTotalNodes = nComputeNodeTotalNodes &
                             + (ElemInfo_Shared(ELEM_LASTNODEIND,iElem) - ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem))
    END IF
  END DO
END IF
#endif  /*USE_MPI*/

! and get max number of bgm-elems
ALLOCATE(Distance    (1:MAXVAL(FIBGM_nElems)) &
        ,ListDistance(1:MAXVAL(FIBGM_nElems)) )

END SUBROUTINE BuildBGMAndIdentifyHaloRegion


END MODULE MOD_Particle_BGM
