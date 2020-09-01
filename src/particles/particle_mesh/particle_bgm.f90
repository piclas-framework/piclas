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

INTERFACE DefineParametersParticleBGM
    MODULE PROCEDURE DefineParametersParticleBGM
END INTERFACE

INTERFACE BuildBGMAndIdentifyHaloRegion
    MODULE PROCEDURE BuildBGMAndIdentifyHaloRegion
END INTERFACE

INTERFACE FinalizeBGM
    MODULE PROCEDURE FinalizeBGM
END INTERFACE

PUBLIC::DefineParametersParticleBGM
PUBLIC::BuildBGMAndIdentifyHaloRegion
PUBLIC :: FinalizeBGM

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
USE MOD_Globals_Vars           ,ONLY: c
USE MOD_Preproc
USE MOD_Mesh_Vars              ,ONLY: nElems,offsetElem,nGlobalElems
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Periodic_BC   ,ONLY: InitPeriodicBC
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod,Distance,ListDistance
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
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,iLocSide,SideID
INTEGER                        :: FirstElem,LastElem
INTEGER                        :: firstNodeID,lastNodeID
INTEGER                        :: offsetNodeID,nNodeIDs,currentOffset,moveBGMindex
REAL                           :: xmin, xmax, ymin, ymax, zmin, zmax
INTEGER                        :: iBGM, jBGM, kBGM
INTEGER                        :: BGMimax, BGMimin, BGMjmax, BGMjmin, BGMkmax, BGMkmin
INTEGER                        :: BGMiDelta,BGMjDelta,BGMkDelta
INTEGER                        :: BGMCellXmax, BGMCellXmin, BGMCellYmax, BGMCellYmin, BGMCellZmax, BGMCellZmin
#if USE_MPI
INTEGER                        :: iSide
INTEGER                        :: ElemID
REAL                           :: deltaT
REAL                           :: globalDiag
INTEGER,ALLOCATABLE            :: sendbuf(:,:,:), recvbuf(:,:,:)
INTEGER,ALLOCATABLE            :: offsetElemsInBGMCell(:,:,:)
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
INTEGER                        :: nHaloElems,nMPISidesShared
INTEGER,ALLOCATABLE            :: offsetCNHalo2GlobalElem(:), offsetMPISideShared(:)
REAL,ALLOCATABLE               :: BoundsOfElemCenter(:), MPISideBoundsOfElemCenter(:,:)
LOGICAL                        :: ElemInsideHalo
INTEGER                        :: firstHaloElem,lastHaloElem
! FIBGMToProc
INTEGER                        :: iProc,ProcRank,nFIBGMToProc,MessageSize
INTEGER                        :: BGMiminglob,BGMimaxglob,BGMjminglob,BGMjmaxglob,BGMkminglob,BGMkmaxglob
LOGICAL,ALLOCATABLE            :: FIBGMToProcTmp(:,:,:,:)
INTEGER,ALLOCATABLE            :: FIBGM_nTotalElemsTmp(:,:,:)
#else
REAL                           :: halo_eps
#endif /*USE_MPI*/
#ifdef CODE_ANALYZE
INTEGER,ALLOCATABLE            :: NumberOfElements(:)
#endif /*CODE_ANALYZE*/
!===================================================================================================================================

! Read parameter for FastInitBackgroundMesh (FIBGM)
GEO%FIBGMdeltas(1:3) = GETREALARRAY('Part-FIBGMdeltas',3,'1. , 1. , 1.')
GEO%FactorFIBGM(1:3) = GETREALARRAY('Part-FactorFIBGM',3,'1. , 1. , 1.')
GEO%FIBGMdeltas(1:3) = 1./GEO%FactorFIBGM(1:3) * GEO%FIBGMdeltas(1:3)

! Read periodic vectors from parameter file
CALL InitPeriodicBC()

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

moveBGMindex = 1 ! BGM indices must be >1 --> move by 1

! Use NodeCoords only for TriaTracking since Tracing and RefMapping have potentially curved elements, only BezierControlPoints form
! convec hull
SELECT CASE(TrackingMethod)
  CASE(TRIATRACKING)
    DO iElem = firstElem, lastElem
      offSetNodeID = ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem)
      nNodeIDs     = ElemInfo_Shared(ELEM_LASTNODEIND ,iElem)-ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem)
      firstNodeID  = offsetNodeID+1
      lastNodeID   = offsetNodeID+nNodeIDs

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

      ! BGM indices must be >0 --> move by 1
      ElemToBGM_Shared(1,iElem) = FLOOR((xmin-GEO%xminglob)/GEO%FIBGMdeltas(1)) + moveBGMindex
      ElemToBGM_Shared(2,iElem) = FLOOR((xmax-GEO%xminglob)/GEO%FIBGMdeltas(1)) + moveBGMindex
      ElemToBGM_Shared(3,iElem) = FLOOR((ymin-GEO%yminglob)/GEO%FIBGMdeltas(2)) + moveBGMindex
      ElemToBGM_Shared(4,iElem) = FLOOR((ymax-GEO%yminglob)/GEO%FIBGMdeltas(2)) + moveBGMindex
      ElemToBGM_Shared(5,iElem) = FLOOR((zmin-GEO%zminglob)/GEO%FIBGMdeltas(3)) + moveBGMindex
      ElemToBGM_Shared(6,iElem) = FLOOR((zmax-GEO%zminglob)/GEO%FIBGMdeltas(3)) + moveBGMindex
    END DO ! iElem = firstElem, lastElem

  CASE(TRACING,REFMAPPING)
    DO iElem = firstElem, lastElem
      xmin= HUGE(1.)
      xmax=-HUGE(1.)
      ymin= HUGE(1.)
      ymax=-HUGE(1.)
      zmin= HUGE(1.)
      zmax=-HUGE(1.)

      DO iLocSide = 1,6
        SideID = GetGlobalNonUniqueSideID(iElem,iLocSide)
        xmin = MIN(xmin,MINVAL(BezierControlPoints3D(1,:,:,SideID)))
        xmax = MAX(xmax,MAXVAL(BezierControlPoints3D(1,:,:,SideID)))
        ymin = MIN(ymin,MINVAL(BezierControlPoints3D(2,:,:,SideID)))
        ymax = MAX(ymax,MAXVAL(BezierControlPoints3D(2,:,:,SideID)))
        zmin = MIN(zmin,MINVAL(BezierControlPoints3D(3,:,:,SideID)))
        zmax = MAX(zmax,MAXVAL(BezierControlPoints3D(3,:,:,SideID)))
      END DO

      ! Restrict to domain extent
      xmin = MAX(xmin,GEO%xminglob)
      xmax = MIN(xmax,GEO%xmaxglob)
      ymin = MAX(ymin,GEO%yminglob)
      ymax = MIN(ymax,GEO%ymaxglob)
      zmin = MAX(zmin,GEO%zminglob)
      zmax = MIN(zmax,GEO%zmaxglob)

      BoundsOfElem_Shared(1,1,iElem) = xmin
      BoundsOfElem_Shared(2,1,iElem) = xmax
      BoundsOfElem_Shared(1,2,iElem) = ymin
      BoundsOfElem_Shared(2,2,iElem) = ymax
      BoundsOfElem_Shared(1,3,iElem) = zmin
      BoundsOfElem_Shared(2,3,iElem) = zmax

      ! BGM indices must be >0 --> move by 1
      ElemToBGM_Shared(1,iElem) = FLOOR((xmin-GEO%xminglob)/GEO%FIBGMdeltas(1)) + moveBGMindex
      ElemToBGM_Shared(2,iElem) = FLOOR((xmax-GEO%xminglob)/GEO%FIBGMdeltas(1)) + moveBGMindex
      ElemToBGM_Shared(3,iElem) = FLOOR((ymin-GEO%yminglob)/GEO%FIBGMdeltas(2)) + moveBGMindex
      ElemToBGM_Shared(4,iElem) = FLOOR((ymax-GEO%yminglob)/GEO%FIBGMdeltas(2)) + moveBGMindex
      ElemToBGM_Shared(5,iElem) = FLOOR((zmin-GEO%zminglob)/GEO%FIBGMdeltas(3)) + moveBGMindex
      ElemToBGM_Shared(6,iElem) = FLOOR((zmax-GEO%zminglob)/GEO%FIBGMdeltas(3)) + moveBGMindex
    END DO ! iElem = firstElem, lastElem
END SELECT

#if USE_MPI
CALL MPI_WIN_SYNC(ElemToBGM_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(BoundsOfElem_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
#endif  /*USE_MPI*/

! deallocate stuff // required for dynamic load balance
#if USE_LOADBALANCE
IF (ALLOCATED(GEO%FIBGM)) THEN
  DO iBGM=GEO%FIBGMimin,GEO%FIBGMimax
    DO jBGM=GEO%FIBGMjmin,GEO%FIBGMjmax
      DO kBGM=GEO%FIBGMkmin,GEO%FIBGMkmax
        SDEALLOCATE(GEO%FIBGM(iBGM,jBGM,kBGM)%Element)
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
  DEALLOCATE(GEO%FIBGM)
END IF
#endif /*USE_LOADBALANCE*/

#if USE_MPI
SafetyFactor  =GETREAL('Part-SafetyFactor','1.0')
halo_eps_velo =GETREAL('Particles-HaloEpsVelo','0')

IF (nComputeNodeProcessors.EQ.nProcessors_Global) THEN
#endif /*USE_MPI*/
  halo_eps  = 0.
#if USE_MPI
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
#endif /*USE_MPI*/

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
    BoundsOfElemCenter(1:3) = (/ SUM(   BoundsOfElem_Shared(1:2,1,ElemID)), &
                                 SUM(   BoundsOfElem_Shared(1:2,2,ElemID)), &
                                 SUM(   BoundsOfElem_Shared(1:2,3,ElemID)) /) / 2.
    BoundsOfElemCenter(4) = VECNORM ((/ BoundsOfElem_Shared(2  ,1,ElemID)-BoundsOfElem_Shared(1,1,ElemID), &
                                        BoundsOfElem_Shared(2  ,2,ElemID)-BoundsOfElem_Shared(1,2,ElemID), &
                                        BoundsOfElem_Shared(2  ,3,ElemID)-BoundsOfElem_Shared(1,3,ElemID) /) / 2.)
    DO iSide = 1, nMPISidesShared
      ! compare distance of centers with sum of element outer radii+halo_eps
      IF (VECNORM(BoundsOfElemCenter(1:3)-MPISideBoundsOfElemCenter(1:3,iSide)) &
          .GT. halo_eps+BoundsOfElemCenter(4)+MPISideBoundsOfElemCenter(4,iSide) ) CYCLE
      ElemInsideHalo = .TRUE.
      EXIT
    END DO ! iSide = 1, nMPISidesShared
    IF (.NOT.ElemInsideHalo) THEN
      ElemInfo_Shared(ELEM_HALOFLAG,ElemID) = 0
    ELSE
      ! Only add element to BGM if inside halo region on node.
      ! THIS IS WRONG. WE ARE WORKING ON THE CN HALO REGION. IF WE OMIT THE
      ! ELEMENT HERE, WE LOOSE IT. IF WE KEEP IT, WE BREAK AT 589. YOUR CALL.
      CALL AddElementToFIBGM(ElemID)
    END IF
  END DO ! iElem = firstHaloElem, lastHaloElem
END IF ! nComputeNodeProcessors.EQ.nProcessors_Global
CALL MPI_WIN_SYNC(ElemInfo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

IF (GEO%nPeriodicVectors.GT.0) CALL CheckPeriodicSides()
#else
ElemInfo_Shared(ELEM_HALOFLAG,:) = 1
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

BGMiDelta = BGMimax - BGMimin
BGMjDelta = BGMjmax - BGMjmin
BGMkDelta = BGMkmax - BGMkmin
! allocated shared memory for nElems per BGM cell
! MPI shared memory is continuous, beginning from 1. All shared arrays have to
! be shifted to BGM[i]min with pointers
ALLOCATE(offsetElemsInBGMCell(BGMimin:BGMimax,BGMjmin:BGMjmax,BGMkmin:BGMkmax))
CALL MPI_EXSCAN(sendbuf(:,:,:),recvbuf(:,:,:),(BGMiDelta+1)*(BGMjDelta+1)*(BGMkDelta+1),MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,iError)
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
! MPI shared memory is continuous, beginning from 1. All shared arrays have to
! be shifted to BGM[i]min with pointers
MPISharedSize = INT((BGMiDelta+1)*(BGMjDelta+1)*(BGMkDelta+1),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/(BGMiDelta+1)*(BGMjDelta+1)*(BGMkDelta+1)/) &
                    ,FIBGM_nElems_Shared_Win,FIBGM_nElems_Shared)
CALL MPI_WIN_LOCK_ALL(0,FIBGM_nElems_Shared_Win,IERROR)
! allocated shared memory for BGM cell offset in 1D array of BGM to element mapping
MPISharedSize = INT((BGMiDelta+1)*(BGMjDelta+1)*(BGMkDelta+1),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/(BGMiDelta+1)*(BGMjDelta+1)*(BGMkDelta+1)/) &
                    ,FIBGM_offsetElem_Shared_Win,FIBGM_offsetElem_Shared)
CALL MPI_WIN_LOCK_ALL(0,FIBGM_offsetElem_Shared_Win,IERROR)
FIBGM_nElems     (BGMimin:BGMimax,BGMjmin:BGMjmax,BGMkmin:BGMkmax) => FIBGM_nElems_Shared
FIBGM_offsetElem (BGMimin:BGMimax,BGMjmin:BGMjmax,BGMkmin:BGMkmax) => FIBGM_offsetElem_Shared

! last proc of compute-node writes into shared memory to make nElems per BGM accessible for every proc
IF(myComputeNodeRank.EQ.nComputeNodeProcessors-1)THEN
  currentOffset = 0
  DO iBGM = BGMimin,BGMimax
    DO jBGM = BGMjmin,BGMjmax
      DO kBGM = BGMkmin,BGMkmax
        ! senfbuf and recvbuf have to stay on original position. Shift 1 --> BGMimin
        FIBGM_nElems(iBGM,jBGM,kBGM)     = sendbuf(iBGM,jBGM,kBGM)
        FIBGM_offsetElem(iBGM,jBGM,kBGM) = currentOffset
        currentOffset = currentoffset    + sendbuf(iBGM,jBGM,kBGM)
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
END IF
DEALLOCATE(sendbuf)
CALL MPI_WIN_SYNC(FIBGM_nElems_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(FIBGM_offsetElem_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#else /*NOT USE_MPI*/
ALLOCATE( FIBGM_nElems    (BGMimin:BGMimax,BGMjmin:BGMjmax,BGMjmin:BGMjmax))
ALLOCATE( FIBGM_offsetElem(BGMimin:BGMimax,BGMjmin:BGMjmax,BGMjmin:BGMjmax))
currentOffset = 0
  DO iBGM = BGMimin,BGMimax
    DO jBGM = BGMjmin,BGMjmax
      DO kBGM = BGMkmin,BGMkmax
      FIBGM_nElems(iBGM,jBGM,kBGM)     = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem
      FIBGM_offsetElem(iBGM,jBGM,kBGM) = currentOffset
      currentOffset = currentoffset    + GEO%FIBGM(iBGM,jBGM,kBGM)%nElem
    END DO ! kBGM
  END DO ! jBGM
END DO ! iBGM
#endif  /*USE_MPI*/

#if USE_MPI
! allocate 1D array for mapping of BGM cell to Element indeces
MPISharedSize = INT((FIBGM_offsetElem(BGMimax,BGMjmax,BGMkmax) + FIBGM_nElems    (BGMimax,BGMjmax,BGMkmax))  &
                     ,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/FIBGM_offsetElem(BGMimax,BGMjmax,BGMkmax)   + &
                                     FIBGM_nElems    (BGMimax,BGMjmax,BGMkmax)/)   &
                     ,FIBGM_Element_Shared_Win,FIBGM_Element_Shared)
CALL MPI_WIN_LOCK_ALL(0,FIBGM_Element_Shared_Win,IERROR)
FIBGM_Element => FIBGM_Element_Shared
#else
ALLOCATE( FIBGM_Element(1:FIBGM_offsetElem(BGMimax,BGMjmax,BGMkmax) + &
                          FIBGM_nElems    (BGMimax,BGMjmax,BGMkmax)))
#endif  /*USE_MPI*/

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
  FIBGM_Element = -1
#if USE_MPI
END IF
#endif /*USE_MPI*/
CALL MPI_WIN_SYNC(FIBGM_Element_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

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
  ! Only add non-peri halo elems
  DO iElem = firstHaloElem, lastHaloElem
    ElemID = offsetCNHalo2GlobalElem(iElem)
    IF (ElemInfo_Shared(ELEM_HALOFLAG,ElemID).NE.2) CYCLE

    BGMCellXmin = MAX(ElemToBGM_Shared(1,ElemID),BGMimin)
    BGMCellXmax = MIN(ElemToBGM_Shared(2,ElemID),BGMimax)
    BGMCellYmin = MAX(ElemToBGM_Shared(3,ElemID),BGMjmin)
    BGMCellYmax = MIN(ElemToBGM_Shared(4,ElemID),BGMjmax)
    BGMCellZmin = MAX(ElemToBGM_Shared(5,ElemID),BGMkmin)
    BGMCellZmax = MIN(ElemToBGM_Shared(6,ElemID),BGMkmax)

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

! Add local elements
DO iElem = offsetElem+1, offsetElem+nElems
  ! find element extent on BGM
  BGMCellXmin = MAX(ElemToBGM_Shared(1,iElem),BGMimin)
  BGMCellXmax = MIN(ElemToBGM_Shared(2,iElem),BGMimax)
  BGMCellYmin = MAX(ElemToBGM_Shared(3,iElem),BGMjmin)
  BGMCellYmax = MIN(ElemToBGM_Shared(4,iElem),BGMjmax)
  BGMCellZmin = MAX(ElemToBGM_Shared(5,iElem),BGMkmin)
  BGMCellZmax = MIN(ElemToBGM_Shared(6,iElem),BGMkmax)

  ! add current element to BGM-Elem
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

! Abort if FIBGM_Element still contains unfilled entries
IF (ANY(FIBGM_Element.EQ.-1)) CALL ABORT(__STAMP__,'Error while filling FIBGM element array')

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
    IF (ElemInfo_Shared(ELEM_HALOFLAG,iElem).NE.0) THEN
      nComputeNodeTotalElems = nComputeNodeTotalElems + 1
    END IF
  END DO
  ALLOCATE(CNTotalElem2GlobalElem(1:nComputeNodeTotalElems))
  ALLOCATE(GlobalElem2CNTotalElem(1:nGlobalElems))
  nComputeNodeTotalElems = 0
  nComputeNodeTotalElems = 0
  GlobalElem2CNTotalElem(1:nGlobalElems) = -1
  ! CN-local elements
  DO iElem = 1,nGlobalElems
    IF (ElemInfo_Shared(ELEM_HALOFLAG,iElem).EQ.1) THEN
      nComputeNodeTotalElems = nComputeNodeTotalElems + 1
      CNTotalElem2GlobalElem(nComputeNodeTotalElems) = iElem
      GlobalElem2CNTotalElem(iElem) = nComputeNodeTotalElems
      nComputeNodeTotalSides = nComputeNodeTotalSides &
                             + (ElemInfo_Shared(ELEM_LASTSIDEIND,iElem) - ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem))
      nComputeNodeTotalNodes = nComputeNodeTotalNodes &
                             + (ElemInfo_Shared(ELEM_LASTNODEIND,iElem) - ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem))
    END IF
  END DO
  ! CN-halo elements (non-periodic)
  DO iElem = 1,nGlobalElems
    IF (ElemInfo_Shared(ELEM_HALOFLAG,iElem).EQ.2) THEN
      nComputeNodeTotalElems = nComputeNodeTotalElems + 1
      CNTotalElem2GlobalElem(nComputeNodeTotalElems) = iElem
      GlobalElem2CNTotalElem(iElem) = nComputeNodeTotalElems
      nComputeNodeTotalSides = nComputeNodeTotalSides &
                             + (ElemInfo_Shared(ELEM_LASTSIDEIND,iElem) - ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem))
      nComputeNodeTotalNodes = nComputeNodeTotalNodes &
                             + (ElemInfo_Shared(ELEM_LASTNODEIND,iElem) - ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem))
    END IF
  END DO
  ! CN-halo elements (periodic)
  DO iElem = 1,nGlobalElems
    IF (ElemInfo_Shared(ELEM_HALOFLAG,iElem).EQ.3) THEN
      nComputeNodeTotalElems = nComputeNodeTotalElems + 1
      CNTotalElem2GlobalElem(nComputeNodeTotalElems) = iElem
      GlobalElem2CNTotalElem(iElem) = nComputeNodeTotalElems
      nComputeNodeTotalSides = nComputeNodeTotalSides &
                             + (ElemInfo_Shared(ELEM_LASTSIDEIND,iElem) - ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem))
      nComputeNodeTotalNodes = nComputeNodeTotalNodes &
                             + (ElemInfo_Shared(ELEM_LASTNODEIND,iElem) - ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem))
    END IF
  END DO
END IF

#ifdef CODE_ANALYZE
! Sanity checks
IF (  SUM(ElemInfo_Shared(ELEM_HALOFLAG,:)  ,MASK=ElemInfo_Shared(ELEM_HALOFLAG,:).EQ.1).NE.nComputeNodeElems) &
  CALL ABORT(__STAMP__,'Error with number of local elements on compute node')

IF ((SUM(ElemInfo_Shared(ELEM_HALOFLAG,:)  ,MASK=ElemInfo_Shared(ELEM_HALOFLAG,:).EQ.1) &
    +SUM(ElemInfo_Shared(ELEM_HALOFLAG,:)/2,MASK=ElemInfo_Shared(ELEM_HALOFLAG,:).EQ.2) &
    +SUM(ElemInfo_Shared(ELEM_HALOFLAG,:)/3,MASK=ElemInfo_Shared(ELEM_HALOFLAG,:).EQ.3)).NE.nComputeNodeTotalElems) &
  CALL ABORT(__STAMP__,'Error with number of halo elements on compute node')

! Debug output
IF (myRank.EQ.0) THEN
  SWRITE(Unit_StdOut,'(A)') ' DETERMINED compute-node (CN) halo region ...'
  SWRITE(Unit_StdOut,'(A)') ' | CN Rank | Local Elements | Halo Elements (non-peri) | Halo Elements (peri) |'
  CALL FLUSH(UNIT_stdOut)
  ALLOCATE(NumberOfElements(3*nLeaderGroupProcs))
END IF

IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_GATHER((/ SUM(ElemInfo_Shared(ELEM_HALOFLAG,:)  ,MASK=ElemInfo_Shared(ELEM_HALOFLAG,:).EQ.1),  &
                     SUM(ElemInfo_Shared(ELEM_HALOFLAG,:)/2,MASK=ElemInfo_Shared(ELEM_HALOFLAG,:).EQ.2),  &
                     SUM(ElemInfo_Shared(ELEM_HALOFLAG,:)/3,MASK=ElemInfo_Shared(ELEM_HALOFLAG,:).EQ.3)/) &
                 , 3 , MPI_INTEGER                                                                        &
                 , NumberOfElements                                                                       &
                 , 3 , MPI_INTEGER                                                                        &
                 , 0 , MPI_COMM_LEADERS_SHARED ,iError)
END IF

IF (myRank.EQ.0) THEN
  DO iProc = 0,nLeaderGroupProcs-1
    WRITE(Unit_StdOut,'(A,I7,A,I15,A,I25,A,I21,A)')  &
                                      ' |>',iProc, &
                                      ' |'  ,NumberOfElements(iProc*3+1), &
                                      ' |'  ,NumberOfElements(iProc*3+2), &
                                      ' |'  ,NumberOfElements(iProc*3+3), ' |'
  END DO
END IF
CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif /*CODE_ANALYZE*/

! Loop over all elements and build a global FIBGM to processor mapping. This is required to identify potential emission procs
BGMiminglob = 0 + moveBGMindex
BGMimaxglob = FLOOR((GEO%xmaxglob-GEO%xminglob)/GEO%FIBGMdeltas(1)) + moveBGMindex
BGMjminglob = 0 + moveBGMindex
BGMjmaxglob = FLOOR((GEO%ymaxglob-GEO%yminglob)/GEO%FIBGMdeltas(2)) + moveBGMindex
BGMkminglob = 0 + moveBGMindex
BGMkmaxglob = FLOOR((GEO%zmaxglob-GEO%zminglob)/GEO%FIBGMdeltas(3)) + moveBGMindex

GEO%FIBGMiminglob = BGMiminglob
GEO%FIBGMimaxglob = BGMimaxglob
GEO%FIBGMjminglob = BGMjminglob
GEO%FIBGMjmaxglob = BGMjmaxglob
GEO%FIBGMkminglob = BGMkminglob
GEO%FIBGMkmaxglob = BGMkmaxglob

firstElem = INT(REAL( myComputeNodeRank   *nGlobalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nGlobalElems)/REAL(nComputeNodeProcessors))

! Flag each FIBGM element proc positive
MPISharedSize = INT((BGMimaxglob-BGMiminglob+1)*(BGMjmaxglob-BGMjminglob+1)*(BGMkmaxglob-BGMkminglob+1)                        &
                    ,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/(BGMimaxglob-BGMiminglob+1)*(BGMjmaxglob-BGMjminglob+1)*(BGMkmaxglob-BGMkminglob+1)/)     &
                    ,FIBGM_nTotalElems_Shared_Win,FIBGM_nTotalElems_Shared)
CALL MPI_WIN_LOCK_ALL(0,FIBGM_nTotalElems_Shared_Win,IERROR)
FIBGM_nTotalElems(BGMiminglob:BGMimaxglob,BGMjminglob:BGMjmaxglob,BGMkminglob:BGMkmaxglob) => FIBGM_nTotalElems_Shared

ALLOCATE(FIBGMToProcTmp      (BGMiminglob:BGMimaxglob,BGMjminglob:BGMjmaxglob,BGMkminglob:BGMkmaxglob,0:nProcessors_Global-1) ,&
         FIBGM_nTotalElemsTmp(BGMiminglob:BGMimaxglob,BGMjminglob:BGMjmaxglob,BGMkminglob:BGMkmaxglob))
FIBGMToProcTmp = .FALSE.
FIBGM_nTotalElemsTmp = 0

DO iElem = firstElem,lastElem
  ProcRank = ElemInfo_Shared(ELEM_RANK,iElem)

  DO kBGM = ElemToBGM_Shared(5,iElem),ElemToBGM_Shared(6,iElem)
    DO jBGM = ElemToBGM_Shared(3,iElem),ElemToBGM_Shared(4,iElem)
      DO iBGM = ElemToBGM_Shared(1,iElem),ElemToBGM_Shared(2,iElem)
        FIBGMToProcTmp      (iBGM,jBGM,kBGM,ProcRank) = .TRUE.
        FIBGM_nTotalElemsTmp(iBGM,jBGM,kBGM)          = FIBGM_nTotalElemsTmp(iBGM,jBGM,kBGM) + 1
      END DO
    END DO
  END DO
END DO

! Sum global number of elements per FIBGM cell
MessageSize = (BGMimaxglob-BGMiminglob+1)*(BGMjmaxglob-BGMjminglob+1)*(BGMkmaxglob-BGMkminglob+1)
CALL MPI_REDUCE(FIBGM_nTotalElemsTmp,FIBGM_nTotalElems,MessageSize,MPI_INTEGER,MPI_SUM,0,MPI_COMM_SHARED,iERROR)
CALL MPI_WIN_SYNC(FIBGM_nTotalElems_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
DEALLOCATE(FIBGM_nTotalElemsTmp)

! Perform logical OR and place data on CN root
MessageSize = (BGMimaxglob-BGMiminglob+1)*(BGMjmaxglob-BGMjminglob+1)*(BGMkmaxglob-BGMkminglob+1)*nProcessors_Global
IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE  ,FIBGMToProcTmp,MessageSize,MPI_LOGICAL,MPI_LOR,0,MPI_COMM_SHARED,iError)
ELSE
  CALL MPI_REDUCE(FIBGMToProcTmp,FIBGMToProcTmp,MessageSize,MPI_LOGICAL,MPI_LOR,0,MPI_COMM_SHARED,iError)
  DEALLOCATE(FIBGMToProcTmp)
END IF
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

! Allocate shared array to hold the mapping
MPISharedSize = INT(3*(BGMimaxglob-BGMiminglob+1)*(BGMjmaxglob-BGMjminglob+1)*(BGMkmaxglob-BGMkminglob+1),MPI_ADDRESS_KIND)&
                                                                                                         *MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/2,BGMimaxglob-BGMiminglob+1,BGMjmaxglob-BGMjminglob+1,BGMkmaxglob-BGMkminglob+1/), &
                                    FIBGMToProc_Shared_Win,FIBGMToProc_Shared)
CALL MPI_WIN_LOCK_ALL(0,FIBGMToProc_Shared_Win,IERROR)
FIBGMToProc => FIBGMToProc_Shared

IF (myComputeNodeRank.EQ.0) THEN
  FIBGMToProc = 0
END IF
CALL MPI_WIN_SYNC(FIBGMToProc_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

! CN root build the mapping to avoid further communication
IF (myComputeNodeRank.EQ.0) THEN
  nFIBGMToProc = 0

  DO kBGM = BGMkminglob,BGMkmaxglob
    DO jBGM = BGMjminglob,BGMjmaxglob
      DO iBGM = BGMiminglob,BGMimaxglob
        ! Save current offset
        FIBGMToProc(FIBGM_FIRSTPROCIND,iBGM,jBGM,kBGM) = nFIBGMToProc
        ! Save number of procs per FIBGM element
        DO iProc = 0,nProcessors_Global-1
          ! Proc belongs to current FIBGM cell
          IF (FIBGMToProcTmp(iBGM,jBGM,kBGM,iProc)) THEN
            nFIBGMToProc = nFIBGMToProc + 1
            FIBGMToProc(FIBGM_NPROCS,iBGM,jBGM,kBGM) = FIBGMToProc(FIBGM_NPROCS,iBGM,jBGM,kBGM) + 1
          END IF
        END DO
      END DO
    END DO
  END DO
END IF

! Synchronize array and communicate the information to other procs on CN node
CALL MPI_WIN_SYNC(FIBGMToProc_Shared_Win,IERROR)
CALL MPI_BCAST(nFIBGMToProc,1,MPI_INTEGER,0,MPI_COMM_SHARED,iError)

! Allocate shared array to hold the proc information
MPISharedSize = INT(nFIBGMToProc,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nFIBGMToProc/),FIBGMProcs_Shared_Win,FIBGMProcs_Shared)
CALL MPI_WIN_LOCK_ALL(0,FIBGMProcs_Shared_Win,IERROR)
FIBGMProcs => FIBGMProcs_Shared

IF (myComputeNodeRank.EQ.0) THEN
  FIBGMProcs= -1
END IF
CALL MPI_WIN_SYNC(FIBGMProcs_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

! CN root fills the information
IF (myComputeNodeRank.EQ.0) THEN
  nFIBGMToProc = 0

  DO kBGM = BGMkminglob,BGMkmaxglob
    DO jBGM = BGMjminglob,BGMjmaxglob
      DO iBGM = BGMiminglob,BGMimaxglob
        ! Save proc ID
        DO iProc = 0,nProcessors_Global-1
          ! Proc belongs to current FIBGM cell
          IF (FIBGMToProcTmp(iBGM,jBGM,kBGM,iProc)) THEN
            nFIBGMToProc = nFIBGMToProc + 1
            FIBGMProcs(nFIBGMToProc) = iProc
          END IF
        END DO
      END DO
    END DO
  END DO

  DEALLOCATE(FIBGMToProcTmp)
END IF

CALL MPI_WIN_SYNC(FIBGMProcs_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/

! and get max number of bgm-elems
ALLOCATE(Distance    (1:MAXVAL(FIBGM_nElems)) &
        ,ListDistance(1:MAXVAL(FIBGM_nElems)) )

! ElemToBGM is only used during init. First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

CALL MPI_WIN_UNLOCK_ALL(ElemToBGM_Shared_Win,iError)
CALL MPI_WIN_FREE(ElemToBGM_Shared_Win,iError)

CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

! Then, free the pointers or arrays
ADEALLOCATE(ElemToBGM_Shared)
#endif /*USE_MPI*/

END SUBROUTINE BuildBGMAndIdentifyHaloRegion


SUBROUTINE FinalizeBGM()
!===================================================================================================================================
! Deallocates variables for the particle background mesh
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_MPI_Shared_Vars
USE MOD_Particle_Mesh_Vars
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

!CALL MPI_WIN_UNLOCK_ALL(ElemToBGM_Shared_Win,iError)
!CALL MPI_WIN_FREE(ElemToBGM_Shared_Win,iError)
CALL MPI_WIN_UNLOCK_ALL(BoundsOfElem_Shared_Win,iError)
CALL MPI_WIN_FREE(BoundsOfElem_Shared_Win,iError)
CALL MPI_WIN_UNLOCK_ALL(FIBGM_nTotalElems_Shared_Win,iError)
CALL MPI_WIN_FREE(FIBGM_nTotalElems_Shared_Win,iError)
CALL MPI_WIN_UNLOCK_ALL(FIBGM_nElems_Shared_Win,iError)
CALL MPI_WIN_FREE(FIBGM_nElems_Shared_Win,iError)
CALL MPI_WIN_UNLOCK_ALL(FIBGM_offsetElem_Shared_Win,iError)
CALL MPI_WIN_FREE(FIBGM_offsetElem_Shared_Win,iError)
CALL MPI_WIN_UNLOCK_ALL(FIBGM_Element_Shared_Win,iError)
CALL MPI_WIN_FREE(FIBGM_Element_Shared_Win,iError)
CALL MPI_WIN_UNLOCK_ALL(FIBGMToProc_Shared_Win,iError)
CALL MPI_WIN_FREE(FIBGMToProc_Shared_Win,iError)
CALL MPI_WIN_UNLOCK_ALL(FIBGMProcs_Shared_Win,iError)
CALL MPI_WIN_FREE(FIBGMProcs_Shared_Win,iError)

CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

! Then, free the pointers or arrays
SDEALLOCATE(CNTotalElem2GlobalElem)
SDEALLOCATE(GlobalElem2CNTotalElem)
#endif /*USE_MPI*/

!ADEALLOCATE(ElemToBGM_Shared)
ADEALLOCATE(BoundsOfElem_Shared)
ADEALLOCATE(FIBGM_nTotalElems_Shared)
ADEALLOCATE(FIBGM_nElems_Shared)
ADEALLOCATE(FIBGM_offsetElem_Shared)
ADEALLOCATE(FIBGM_Element_Shared)
ADEALLOCATE(FIBGMToProc)
ADEALLOCATE(FIBGMToProc_Shared)
ADEALLOCATE(FIBGMProcs)
ADEALLOCATE(FIBGMProcs_Shared)

END SUBROUTINE FinalizeBGM


#if USE_MPI
SUBROUTINE CheckPeriodicSides()
!===================================================================================================================================
!> checks the elements against periodic distance
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars              ,ONLY: BoundaryType,nGlobalElems
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound
USE MOD_MPI_Shared_Vars
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemInfo_Shared,SideInfo_Shared,BoundsOfElem_Shared
USE MOD_Particle_MPI_Vars      ,ONLY: halo_eps
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                        :: hasPeriodic
INTEGER                        :: iElem,firstElem,lastElem
INTEGER                        :: iSide,firstSide,lastSide
INTEGER                        :: iPeriodicElem
INTEGER                        :: iPeriodicVector,jPeriodicVector
INTEGER                        :: nPeriodicElems
REAL                           :: BoundsOfElemCenter(1:4)
REAL,ALLOCATABLE               :: PeriodicSideBoundsOfElemCenter(:,:)
INTEGER,ALLOCATABLE            :: nPeriodicVectorsPerElem(:,:)
!===================================================================================================================================

firstElem = INT(REAL( myComputeNodeRank   *nGlobalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nGlobalElems)/REAL(nComputeNodeProcessors))

! count number of elements with periodic sides
nPeriodicElems = 0
DO iElem = 1,nGlobalElems
  ! only consider elements within the DG region
  IF (ElemInfo_Shared(ELEM_HALOFLAG,iElem).EQ.0) CYCLE

  firstSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem) + 1
  lastSide  = ElemInfo_Shared(ELEM_LASTSIDEIND ,iElem)

  DO iSide = firstSide,lastSide
    ! Side is a boundary side
    IF (SideInfo_Shared(SIDE_BCID,iSide).GT.0) THEN
      ! Boundary is a periodic boundary
      IF (PartBound%TargetBoundCond(SideInfo_Shared(SIDE_BCID,iSide)).EQ.3) THEN
        nPeriodicElems = nPeriodicElems + 1
        EXIT
      END IF
    END IF
  END DO
END DO

! return if there are no periodic elements on the compute node or inside the halo region
IF (nPeriodicElems.EQ.0) RETURN

ALLOCATE(PeriodicSideBoundsOfElemCenter(1:4,1:nPeriodicElems))
ALLOCATE(nPeriodicVectorsPerElem(1:3,1:nPeriodicElems))

nPeriodicElems = 0
nPeriodicVectorsPerElem = 0

! Every proc checks every periodic element. It is assumed that there are not "too many" periodic elements. In order to parallelize
! this loop, three communications steps would be required. First the offset of the periodic elements must be communicated, then a
! shared array holding the metrics of the periodic elements must be allocated and filled. Finally, the shared array must be
! synchronized. The communication overhead would most likely exceed the calculation effort.
DO iElem = 1,nGlobalElems
  ! only consider elements within the DG or halo region
  IF (ElemInfo_Shared(ELEM_HALOFLAG,iElem).EQ.0) CYCLE

  firstSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem) + 1
  lastSide  = ElemInfo_Shared(ELEM_LASTSIDEIND ,iElem)

  hasPeriodic = .FALSE.
  DO iSide = firstSide,lastSide
    ! Side is a boundary side
    IF (SideInfo_Shared(SIDE_BCID,iSide).GT.0) THEN
    ! Boundary is a periodic boundary
      IF (PartBound%TargetBoundCond(SideInfo_Shared(SIDE_BCID,iSide)).EQ.3) THEN
        IF (.NOT.hasPeriodic) THEN
          nPeriodicElems = nPeriodicElems + 1
          hasPeriodic    = .TRUE.
        END IF

        ! check the orientation of the periodic side and flip vector to point in the other direction
        IF     (BoundaryType(SideInfo_Shared(SIDE_BCID,iSide),BC_ALPHA).GT.0) THEN
          nPeriodicVectorsPerElem(    BoundaryType(SideInfo_Shared(SIDE_BCID,iSide),BC_ALPHA) ,nPeriodicElems) = -1
        ELSEIF (BoundaryType(SideInfo_Shared(SIDE_BCID,iSide),BC_ALPHA).LT.0) THEN
          nPeriodicVectorsPerElem(ABS(BoundaryType(SideInfo_Shared(SIDE_BCID,iSide),BC_ALPHA)),nPeriodicElems) = 1
        END IF
      END IF
    END IF
  END DO

  IF (hasPeriodic) THEN
    PeriodicSideBoundsOfElemCenter(1:3,nPeriodicElems) = (/ SUM(BoundsOfElem_Shared(1:2,1,iElem)),                         &
                                                            SUM(BoundsOfElem_Shared(1:2,2,iElem)),                         &
                                                            SUM(BoundsOfElem_Shared(1:2,3,iElem)) /) / 2.
    PeriodicSideBoundsOfElemCenter(4,nPeriodicElems) = VECNORM ((/ BoundsOfElem_Shared(2,1,iElem)-BoundsOfElem_Shared(1,1,iElem), &
                                                                   BoundsOfElem_Shared(2,2,iElem)-BoundsOfElem_Shared(1,2,iElem), &
                                                                   BoundsOfElem_Shared(2,3,iElem)-BoundsOfElem_Shared(1,3,iElem) /) / 2.)
  END IF
END DO

! The code below changes ElemInfo_Shared, identification of periodic elements must complete before
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)

! This is a distributed loop. Nonetheless, the load will be unbalanced due to the location of the space-filling curve. Still,
! this approach is again preferred compared to the communication overhead.
DO iElem = firstElem,lastElem
  ! only consider elements that are not already flagged
  IF (ElemInfo_Shared(ELEM_HALOFLAG,iElem).NE.0) CYCLE

  BoundsOfElemCenter(1:3) = (/ SUM(   BoundsOfElem_Shared(1:2,1,iElem)),                                                   &
                               SUM(   BoundsOfElem_Shared(1:2,2,iElem)),                                                   &
                               SUM(   BoundsOfElem_Shared(1:2,3,iElem)) /) / 2.
  BoundsOfElemCenter(4) = VECNORM ((/ BoundsOfElem_Shared(2  ,1,iElem)-BoundsOfElem_Shared(1,1,iElem),                     &
                                      BoundsOfElem_Shared(2  ,2,iElem)-BoundsOfElem_Shared(1,2,iElem),                     &
                                      BoundsOfElem_Shared(2  ,3,iElem)-BoundsOfElem_Shared(1,3,iElem) /) / 2.)

  DO iPeriodicElem = 1,nPeriodicElems
    ! element might be already added back
    IF (ElemInfo_Shared(ELEM_HALOFLAG,iElem).NE.0) EXIT

    SELECT CASE(SUM(ABS(nPeriodicVectorsPerElem(:,iPeriodicElem))))

      CASE(1)
        ! check the only possible periodic vector
        iPeriodicVector = FINDLOC(ABS(nPeriodicVectorsPerElem(:,iPeriodicElem)),1,1)
        IF (iPeriodicVector.EQ.-1) &
          CALL ABORT(__STAMP__,'Error determining periodic vector!')
        ! check if element is within halo_eps of periodically displaced element
        IF (VECNORM( BoundsOfElemCenter(1:3)                                                                               &
                   + GEO%PeriodicVectors(1:3,iPeriodicVector) * nPeriodicVectorsPerElem(iPeriodicVector,iPeriodicElem)     &
                   - PeriodicSideBoundsOfElemCenter(1:3,iPeriodicElem))                                                    &
                .LE. halo_eps+BoundsOfElemCenter(4)+PeriodicSideBoundsOfElemCenter(4,iPeriodicElem) ) THEN
          ! add element back to halo region
          ElemInfo_Shared(ELEM_HALOFLAG,iElem) = 3
!          CALL AddElementToFIBGM(iElem)
        END IF

      CASE(2)
        ! check the two possible periodic vectors. Begin with checking the first periodic vector, followed by the combination of
        ! the first periodic vector with the other. Finally check the second periodic vector, i.e. 1, 1+2, 2
        DO iPeriodicVector = 1,3
          ! element might be already added back
          IF (ElemInfo_Shared(ELEM_HALOFLAG,iElem).NE.0) EXIT

          IF (nPeriodicVectorsPerElem(iPeriodicVector,iPeriodicElem).EQ.0) CYCLE

          ! check if element is within halo_eps of periodically displaced element
          IF (VECNORM( BoundsOfElemCenter(1:3)                                                                             &
                     + GEO%PeriodicVectors(1:3,iPeriodicVector) * nPeriodicVectorsPerElem(iPeriodicVector,iPeriodicElem)   &
                     - PeriodicSideBoundsOfElemCenter(1:3,iPeriodicElem))                                                  &
                    .LE. halo_eps+BoundsOfElemCenter(4)+PeriodicSideBoundsOfElemCenter(4,iPeriodicElem) ) THEN
            ! add element back to halo region
            ElemInfo_Shared(ELEM_HALOFLAG,iElem) = 3
!            CALL AddElementToFIBGM(iElem)
            EXIT
          END IF

          DO jPeriodicVector = 1,3
            IF (nPeriodicVectorsPerElem(jPeriodicVector,iPeriodicElem).EQ.0) CYCLE
            IF (iPeriodicVector.GE.jPeriodicVector) CYCLE

            ! check if element is within halo_eps of periodically displaced element
            IF (VECNORM( BoundsOfElemCenter(1:3)                                                                           &
                       + GEO%PeriodicVectors(1:3,iPeriodicVector) * nPeriodicVectorsPerElem(iPeriodicVector,iPeriodicElem) &
                       + GEO%PeriodicVectors(1:3,jPeriodicVector) * nPeriodicVectorsPerElem(jPeriodicVector,iPeriodicElem) &
                       - PeriodicSideBoundsOfElemCenter(1:3,iPeriodicElem))                                                &
                    .LE. halo_eps+BoundsOfElemCenter(4)+PeriodicSideBoundsOfElemCenter(4,iPeriodicElem) ) THEN
              ! add element back to halo region
              ElemInfo_Shared(ELEM_HALOFLAG,iElem) = 3
!              CALL AddElementToFIBGM(iElem)
              EXIT
            END IF
          END DO
        END DO

      CASE(3)
        ! check the three periodic vectors. Begin with checking the first periodic vector, followed by the combination of
        ! the first periodic vector with the others. Then check the other combinations, i.e. 1, 1+2, 1+3, 2, 2+3, 3, 1+2+3
        DO iPeriodicVector = 1,3
          ! element might be already added back
          IF (ElemInfo_Shared(ELEM_HALOFLAG,iElem).NE.0) EXIT

          ! check if element is within halo_eps of periodically displaced element
          IF (VECNORM( BoundsOfElemCenter(1:3)                                                                             &
                     + GEO%PeriodicVectors(1:3,iPeriodicVector) * nPeriodicVectorsPerElem(iPeriodicVector,iPeriodicElem)   &
                     - PeriodicSideBoundsOfElemCenter(1:3,iPeriodicElem))                                                  &
                    .LE. halo_eps+BoundsOfElemCenter(4)+PeriodicSideBoundsOfElemCenter(4,iPeriodicElem) ) THEN
            ! add element back to halo region
            ElemInfo_Shared(ELEM_HALOFLAG,iElem) = 3
!            CALL AddElementToFIBGM(iElem)
            EXIT
          END IF

          DO jPeriodicVector = 1,3
            IF (iPeriodicVector.GE.jPeriodicVector) CYCLE

            ! check if element is within halo_eps of periodically displaced element
            IF (VECNORM( BoundsOfElemCenter(1:3)                                                                           &
                       + GEO%PeriodicVectors(1:3,iPeriodicVector) * nPeriodicVectorsPerElem(iPeriodicVector,iPeriodicElem) &
                       + GEO%PeriodicVectors(1:3,jPeriodicVector) * nPeriodicVectorsPerElem(jPeriodicVector,iPeriodicElem) &
                       - PeriodicSideBoundsOfElemCenter(1:3,iPeriodicElem))                                                &
                    .LE. halo_eps+BoundsOfElemCenter(4)+PeriodicSideBoundsOfElemCenter(4,iPeriodicElem) ) THEN
              ! add element back to halo region
              ElemInfo_Shared(ELEM_HALOFLAG,iElem) = 3
!              CALL AddElementToFIBGM(iElem)
              EXIT
            END IF

          END DO
        END DO

        ! check if element is within halo_eps of periodically displaced element
        IF (VECNORM( BoundsOfElemCenter(1:3)                                                                               &
                   + GEO%PeriodicVectors(1:3,1) * nPeriodicVectorsPerElem(1,iPeriodicElem)                                 &
                   + GEO%PeriodicVectors(1:3,2) * nPeriodicVectorsPerElem(2,iPeriodicElem)                                 &
                   + GEO%PeriodicVectors(1:3,3) * nPeriodicVectorsPerElem(3,iPeriodicElem)                                 &
                   - PeriodicSideBoundsOfElemCenter(1:3,iPeriodicElem))                                                    &
                .LE. halo_eps+BoundsOfElemCenter(4)+PeriodicSideBoundsOfElemCenter(4,iPeriodicElem) ) THEN
          ! add element back to halo region
          ElemInfo_Shared(ELEM_HALOFLAG,iElem) = 3
!          CALL AddElementToFIBGM(iElem)
        END IF

      CASE DEFAULT
        CALL ABORT(__STAMP__,'Periodic vectors .LT.1 or .GT.3!')
      END SELECT
  END DO
END DO

END SUBROUTINE CheckPeriodicSides


SUBROUTINE AddElementToFIBGM(ElemID)
!===================================================================================================================================
!> adds an element to all corresponding FIBGM cells and ensures correct bounds
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemToBGM_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: ElemID
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iBGM,jBGM,kBGM
INTEGER                        :: BGMCellXmax,BGMCellXmin,BGMCellYmax,BGMCellYmin,BGMCellZmax,BGMCellZmin
!===================================================================================================================================

BGMCellXmin = MAX(ElemToBGM_Shared(1,ElemID),GEO%FIBGMimin)
BGMCellXmax = MIN(ElemToBGM_Shared(2,ElemID),GEO%FIBGMimax)
BGMCellYmin = MAX(ElemToBGM_Shared(3,ElemID),GEO%FIBGMjmin)
BGMCellYmax = MIN(ElemToBGM_Shared(4,ElemID),GEO%FIBGMjmax)
BGMCellZmin = MAX(ElemToBGM_Shared(5,ElemID),GEO%FIBGMkmin)
BGMCellZmax = MIN(ElemToBGM_Shared(6,ElemID),GEO%FIBGMkmax)

! add current element to number of BGM-elems
DO iBGM = BGMCellXmin,BGMCellXmax
  DO jBGM = BGMCellYmin,BGMCellYmax
    DO kBGM = BGMCellZmin,BGMCellZmax
      GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1
    END DO ! kBGM
  END DO ! jBGM
END DO ! iBGM

END SUBROUTINE


#if GCC_VERSION < 90000
PURE FUNCTION FINDLOC(Array,Value,Dim)
!===================================================================================================================================
!> Implements a subset of the intrinsic FINDLOC function for Fortran < 2008
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: Array(:)
INTEGER,INTENT(IN)             :: Value
INTEGER,INTENT(IN)             :: Dim
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER                        :: FINDLOC
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iVar
!===================================================================================================================================
DO iVar = 1,SIZE(ARRAY,1)
  IF (Array(iVar).EQ.Value) THEN
    FINDLOC = iVar
    RETURN
  END IF
END DO

! Return error code -1 if the value was not found
FINDLOC = -1

END FUNCTION FINDLOC
#endif /*GCC_VERSION < 90000*/
#endif /*USE_MPI*/


END MODULE MOD_Particle_BGM
