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

MODULE MOD_Particle_MPI_Halo
!===================================================================================================================================
! Contains routines to build the halo exchange
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

#if USE_MPI
INTERFACE IdentifyPartExchangeProcs
  MODULE PROCEDURE IdentifyPartExchangeProcs
END INTERFACE

INTERFACE FinalizePartExchangeProcs
  MODULE PROCEDURE FinalizePartExchangeProcs
END INTERFACE

PUBLIC :: IdentifyPartExchangeProcs
PUBLIC :: FinalizePartExchangeProcs
!===================================================================================================================================

CONTAINS

SUBROUTINE IdentifyPartExchangeProcs
!===================================================================================================================================
! Identifies processors in physical range for particle exchange communication. This communication has to occur at every RK step and
! would be too costly if done as an all-to-all communication
!> Procs on the same compute node are always assumed to be in range since communication is handled on proc
!> Procs in the compute node halo region are only considered in range if they lie within halo_eps of the mesh on the current proc
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: c,ProjectName
USE MOD_Preproc
USE MOD_HDF5_Output_ElemData    ,ONLY: WriteMyInvisibleRankToHDF5
USE MOD_IO_HDF5                 ,ONLY: AddToElemData,ElementOut
USE MOD_Mesh_Vars               ,ONLY: nElems,offsetElem,myInvisibleRank
USE MOD_Mesh_Tools              ,ONLY: GetGlobalElemID,GetCNElemID
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
USE MOD_MPI_Vars                ,ONLY: offsetElemMPI
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO,MeshHasRotPeriodic,nComputeNodeElems
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemInfo_Shared,SideInfo_Shared,BoundsOfElem_Shared,NodeCoords_Shared
USE MOD_Particle_MPI_Vars       ,ONLY: SafetyFactor,halo_eps,halo_eps_velo,MPI_halo_eps,halo_eps_woshape
USE MOD_Particle_MPI_Vars       ,ONLY: nExchangeProcessors,ExchangeProcToGlobalProc,GlobalProcToExchangeProc,CheckExchangeProcs
USE MOD_Particle_MPI_Vars       ,ONLY: AbortExchangeProcs
USE MOD_Particle_Surfaces_Vars  ,ONLY: BezierControlPoints3D
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
USE MOD_PICDepo_Vars            ,ONLY: DepositionType, ShapeElemProcSend_Shared, ShapeElemProcSend_Shared_Win
USE MOD_PICDepo_Vars            ,ONLY: SendElemShapeID, CNRankToSendRank, nShapeExchangeProcs
USE MOD_PICDepo_Vars            ,ONLY: ShapeMapping,CNShapeMapping,r_sf
USE MOD_ReadInTools             ,ONLY: PrintOption
USE MOD_TimeDisc_Vars           ,ONLY: ManualTimeStep
#if ! (USE_HDG)
USE MOD_CalcTimeStep            ,ONLY: CalcTimeStep
#endif
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
USE MOD_TimeDisc_Vars           ,ONLY: nRKStages,RK_c
#elif (PP_TimeDiscMethod==400)
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Partner identification
LOGICAL                        :: ProcInRange
INTEGER                        :: iPeriodicVector,jPeriodicVector,iPeriodicDir,jPeriodicDir,kPeriodicDir
INTEGER,DIMENSION(2),PARAMETER :: DirPeriodicVector = (/-1,1/)
REAL,DIMENSION(6)              :: xCoordsProc,xCoordsOrigin
INTEGER                        :: iElem,ElemID,firstElem,lastElem,NbElemID,localElem
INTEGER                        :: iSide,SideID,firstSide,lastSide,iLocSide
INTEGER                        :: iMortar,nMortarElems,NbSideID
INTEGER                        :: iProc,HaloProc
INTEGER                        :: nExchangeSides
INTEGER,ALLOCATABLE            :: ExchangeSides(:)
REAL,ALLOCATABLE               :: BoundsOfElemCenter(:),MPISideBoundsOfElemCenter(:,:)
INTEGER                        :: ExchangeProcLeader
LOGICAL,ALLOCATABLE            :: MPISideElem(:)
LOGICAL                        :: ProcHasExchangeElem
! halo_eps reconstruction
REAL                           :: MPI_halo_eps_velo,MPI_halo_diag,vec(1:3),deltaT, MPI_halo_eps_woshape
LOGICAL                        :: fullMesh
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
INTEGER                        :: iStage
#elif (PP_TimeDiscMethod==400)
REAL                           :: NbElemBounds(2,3)
REAL,ALLOCATABLE               :: MPISideBoundsOfNbElemCenter(:,:)
#endif
! shape function
INTEGER                        :: GlobalElemID,GlobalElemRank,GlobalLeaderRank
LOGICAL,ALLOCATABLE            :: FlagShapeElem(:)
! Non-symmetric particle exchange
INTEGER,ALLOCATABLE            :: SendRequest(:),RecvRequest(:),SendShapeElemID(:),RecvProcsElems(:)
LOGICAL,ALLOCATABLE            :: GlobalProcToRecvProc(:), RecvProcs(:)
LOGICAL,ALLOCATABLE            :: CommFlag(:)
INTEGER                        :: nNonSymmetricExchangeProcs,nNonSymmetricExchangeProcsGlob
INTEGER                        :: nExchangeProcessorsGlobal, nSendShapeElems, CNElemID, exElem, exProc, jProc, ProcID
REAL,ALLOCATABLE               :: RotBoundsOfElemCenter(:)
#if (PP_TimeDiscMethod==400)
LOGICAL                        :: SideIsRotPeriodic
INTEGER                        :: BCindex
#endif
!=================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' IDENTIFYING Particle Exchange Processors  ...'

! Allocate arrays
ALLOCATE(GlobalProcToExchangeProc(EXCHANGE_PROC_SIZE,0:nProcessors_Global-1))
GlobalProcToExchangeProc(:,:) = -1

! Identify all procs on same node
nExchangeProcessors = 0

! Communicate non-symmetric part exchange partners to catch non-symmetric proc identification due to inverse distance calculation
ALLOCATE(GlobalProcToRecvProc(0:nProcessors_Global-1), &
         SendRequest         (0:nProcessors_Global-1), &
         RecvRequest         (0:nProcessors_Global-1), &
         CommFlag            (0:nProcessors_Global-1))

GlobalProcToRecvProc = .FALSE.

! Identify all procs with elements in range. This includes checking the procs on the compute-node as they might lie far apart
IF (nProcessors.EQ.1) THEN
  SWRITE(UNIT_stdOut,'(A)') ' | Running on one processor. Particle exchange communication disabled.'
  SWRITE(UNIT_stdOut,'(A)') ' IDENTIFYING Particle Exchange Processors DONE!'
  SWRITE(UNIT_StdOut,'(132("-"))')
  RETURN
END IF

! Open receive buffer for non-symmetric exchange identification
IF (CheckExchangeProcs) THEN
  DO iProc = 0,nProcessors_Global-1
    IF (iProc.EQ.myRank) CYCLE

    CALL MPI_IRECV( GlobalProcToRecvProc(iProc)  &
                  , 1                            &
                  , MPI_LOGICAL                  &
                  , iProc                        &
                  , 1999                         &
                  , MPI_COMM_WORLD               &
                  , RecvRequest(iProc)           &
                  , IERROR)
  END DO
END IF
!> Count all MPI sides on current proc.
firstElem = offsetElem+1
lastElem  = offsetElem+nElems

IF(StringBeginsWith(DepositionType,'shape_function'))THEN
  ALLOCATE(FlagShapeElem(1:nComputeNodeTotalElems))
  FlagShapeElem = .FALSE.
END IF

!>>> For all elements, loop over the six sides and check if the neighbor element is on the current proc
!>>> Special care for big mortar sides, here the SIDE_ELEMID must be used
nExchangeSides = 0

#if !(PP_TimeDiscMethod==400)
ALLOCATE(MPISideElem(offsetElem+1:offsetElem+nElems))
MPISideElem = .FALSE.
#endif /*!(PP_TimeDiscMethod==400)*/

DO iElem = firstElem,lastElem
#if !(PP_TimeDiscMethod==400)
  ! Element already flagged
  IF (MPISideElem(iElem)) CYCLE
#endif /*!(PP_TimeDiscMethod==400)*/

  DO iLocSide = 1,6
    SideID   = GetGlobalNonUniqueSideID(iElem,iLocSide)
    NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)

    ! Mortar side
    IF (NbElemID.LT.0) THEN
      nMortarElems = MERGE(4,2,SideInfo_Shared(SIDE_NBELEMID,SideID).EQ.-1)

      DO iMortar = 1,nMortarElems
        NbSideID = -SideInfo_Shared(SIDE_LOCALID,SideID + iMortar)
        ! If small mortar side not defined, skip it for now, likely not inside the halo region
        IF (NbSideID.LT.1) CYCLE

        NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideID + iMortar)
        ! If small mortar element not defined, skip it for now, likely not inside the halo region
        IF (NbElemID.LT.1) CYCLE

        ! If any of the small mortar sides is not on the local proc, the side is a MPI side
        IF (NbElemID.LT.firstElem .OR. NbElemID.GT.lastElem) THEN
          nExchangeSides = nExchangeSides + 1
#if !(PP_TimeDiscMethod==400)
          MPISideElem(iElem) = .TRUE.
#endif /*!(PP_TimeDiscMethod==400)*/
          EXIT
        END IF
      END DO

    ! regular side or small mortar side
    ELSE
      ! Only check inner (MPI interfaces) and boundary sides (NbElemID.EQ.0)
      ! Boundary sides cannot be discarded because one proc might have MPI interfaces and the other might not
      ! NbElemID.LT.firstElem is always true for NbElemID.EQ.0 because firstElem.GE.1
      IF (NbElemID.LT.firstElem .OR. NbElemID.GT.lastElem) THEN
        nExchangeSides = nExchangeSides + 1
#if !(PP_TimeDiscMethod==400)
        MPISideElem(iElem) = .TRUE.
#endif /*!(PP_TimeDiscMethod==400)*/
      END IF
    END IF
  END DO
END DO

IF (nComputeNodeProcessors.GT.1.AND.nExchangeSides.EQ.0) &
  CALL ABORT(__STAMP__,'Found no side connectivity between processor domains')

!> Build mapping for all MPI sides on current proc
ALLOCATE(ExchangeSides(1:nExchangeSides))

nExchangeSides   = 0
#if !(PP_TimeDiscMethod==400)
MPISideElem = .FALSE.
#endif /*!(PP_TimeDiscMethod==400)*/

DO iElem = firstElem,lastElem
#if !(PP_TimeDiscMethod==400)
  ! Element already flagged
  IF (MPISideElem(iElem)) CYCLE
#endif /*!(PP_TimeDiscMethod==400)*/

  DO iLocSide = 1,6
    SideID   = GetGlobalNonUniqueSideID(iElem,iLocSide)
    NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)

    ! Mortar side
    IF (NbElemID.LT.0) THEN
      nMortarElems = MERGE(4,2,SideInfo_Shared(SIDE_NBELEMID,SideID).EQ.-1)

      DO iMortar = 1,nMortarElems
        NbSideID = -SideInfo_Shared(SIDE_LOCALID,SideID + iMortar)
        ! If small mortar side not defined, skip it for now, likely not inside the halo region
        IF (NbSideID.LT.1) CYCLE

        NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideID + iMortar)
        ! If small mortar element not defined, skip it for now, likely not inside the halo region
        IF (NbElemID.LT.1) CYCLE

        ! If any of the small mortar sides is not on the local proc, the side is a MPI side
        IF (NbElemID.LT.firstElem .OR. NbElemID.GT.lastElem) THEN
          nExchangeSides = nExchangeSides + 1
          ExchangeSides(nExchangeSides) = SideID
#if !(PP_TimeDiscMethod==400)
          MPISideElem(iElem) = .TRUE.
#endif /*!(PP_TimeDiscMethod==400)*/
          EXIT
        END IF
      END DO

    ! regular side or small mortar side
    ELSE
      ! Only check inner (MPI interfaces) and boundary sides (NbElemID.EQ.0)
      ! Boundary sides cannot be discarded because one proc might have MPI interfaces and the other might not
      ! NbElemID.LT.firstElem is always true for NbElemID.EQ.0 because firstElem.GE.1
      IF (NbElemID.LT.firstElem .OR. NbElemID.GT.lastElem) THEN
        nExchangeSides = nExchangeSides + 1
        ExchangeSides(nExchangeSides) = SideID
#if !(PP_TimeDiscMethod==400)
        MPISideElem(iElem) = .TRUE.
#endif /*!(PP_TimeDiscMethod==400)*/
      END IF
    END IF
  END DO
END DO

#if !(PP_TimeDiscMethod==400)
DEALLOCATE(MPISideElem)
#endif /*!(PP_TimeDiscMethod==400)*/

!> Build metrics for all MPI sides on current proc
ALLOCATE(BoundsOfElemCenter(   1:5))
ALLOCATE(RotBoundsOfElemCenter(1:4))
ALLOCATE(MPISideBoundsOfElemCenter(  1:4,1:nExchangeSides))
#if (PP_TimeDiscMethod==400)
ALLOCATE(MPISideBoundsOfNbElemCenter(1:4,1:nExchangeSides))
#endif /*(PP_TimeDiscMethod==400)*/

DO iSide = 1, nExchangeSides
  SideID = ExchangeSides(iSide)
  ElemID = SideInfo_Shared(SIDE_ELEMID,SideID)
  MPISideBoundsOfElemCenter(1:3,iSide) = (/ SUM(  BoundsOfElem_Shared(1:2,1,ElemID)), &
                                            SUM(  BoundsOfElem_Shared(1:2,2,ElemID)), &
                                            SUM(  BoundsOfElem_Shared(1:2,3,ElemID)) /) / 2.
  MPISideBoundsOfElemCenter(4,iSide) = VECNORM ((/BoundsOfElem_Shared(2  ,1,ElemID)-BoundsOfElem_Shared(1,1,ElemID), &
                                                  BoundsOfElem_Shared(2  ,2,ElemID)-BoundsOfElem_Shared(1,2,ElemID), &
                                                  BoundsOfElem_Shared(2  ,3,ElemID)-BoundsOfElem_Shared(1,3,ElemID) /) / 2.)

#if (PP_TimeDiscMethod==400)
  NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)

  ! Mortar elements or normal elements (inkl. rot periodic sides)
  ! Check if SideID is a rot periodic BC
  BCindex = SideInfo_Shared(SIDE_BCID,SideID)
  IF(BCindex.GT.0)THEN
    SideIsRotPeriodic = PartBound%TargetBoundCond(PartBound%MapToPartBC(BCindex)).EQ.PartBound%RotPeriodicBC
    ! NbElemID for BCSide neighbour is zero. Therefore, change it to the ID of the considered element
    NbElemID = SideInfo_Shared(SIDE_ELEMID,SideID)
  ELSE
    SideIsRotPeriodic = .FALSE.
  END IF ! BCindex.GT.0

  ! Only mortar (MPI interfaces) sides: large mortar side
  IF ((NbElemID.LT.1).AND.(.NOT.SideIsRotPeriodic)) THEN
    MPISideBoundsOfNbElemCenter(1:4,iSide) = 0.0
    nMortarElems = MERGE(4,2,SideInfo_Shared(SIDE_NBELEMID,SideID).EQ.-1)
    NbElemBounds(1,:) = HUGE(1.)
    NbElemBounds(2,:) = -HUGE(1.)
    DO iMortar = 1, nMortarElems
      NbSideID = -SideInfo_Shared(SIDE_LOCALID,SideID + iMortar)
      ! If small mortar side not defined, skip it for now, likely not inside the halo region
      IF (NbSideID.LT.1) CALL abort(__STAMP__,'Neighbour side for exchange side missing!')

      NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideID + iMortar)
      ! If small mortar element not defined, skip it for now, likely not inside the halo region
      IF (NbElemID.LT.1) CALL abort(__STAMP__,'Neighbour element for exchange side missing!')

      ! If any of the small mortar sides is not on the local proc, the side is a MPI side
      IF (NbElemID.LT.firstElem .OR. NbElemID.GT.lastElem) THEN
        NbElemBounds(1,:) = MIN(NbElemBounds(1,:),BoundsOfElem_Shared(1,:,NbElemID))
        NbElemBounds(2,:) = MAX(NbElemBounds(2,:),BoundsOfElem_Shared(2,:,NbElemID))
      END IF
    END DO ! iMortar = 1, nMortarElems
    MPISideBoundsOfNbElemCenter(1:3,iSide) = (/ SUM(  NbElemBounds(1:2,1)), &
                                                SUM(  NbElemBounds(1:2,2)), &
                                                SUM(  NbElemBounds(1:2,3)) /) / 2.
    MPISideBoundsOfNbElemCenter(4,iSide) = VECNORM ((/NbElemBounds(2  ,1)-NbElemBounds(1,1), &
                                                      NbElemBounds(2  ,2)-NbElemBounds(1,2), &
                                                      NbElemBounds(2  ,3)-NbElemBounds(1,3) /) / 2.)
  ELSE
    ! Non-mortar (MPI interfaces) sides and rotperiodic sides (BC sides)
    MPISideBoundsOfNbElemCenter(1:3,iSide) = (/ SUM(  BoundsOfElem_Shared(1:2,1,NbElemID)), &
                                                SUM(  BoundsOfElem_Shared(1:2,2,NbElemID)), &
                                                SUM(  BoundsOfElem_Shared(1:2,3,NbElemID)) /) / 2.
    MPISideBoundsOfNbElemCenter(4,iSide) = VECNORM ((/BoundsOfElem_Shared(2  ,1,NbElemID)-BoundsOfElem_Shared(1,1,NbElemID), &
                                                      BoundsOfElem_Shared(2  ,2,NbElemID)-BoundsOfElem_Shared(1,2,NbElemID), &
                                                      BoundsOfElem_Shared(2  ,3,NbElemID)-BoundsOfElem_Shared(1,3,NbElemID) /) / 2.)
  END IF
#endif /*(PP_TimeDiscMethod==400)*/
END DO

!> Check all elements in the CN halo region against local MPI sides. Check is identical to particle_bgm.f90
!>>> Check the bounding box of each element in compute-nodes' halo domain against the bounding boxes of the elements of the
!>>> MPI-surface (local proc MPI sides)

! if running on one node, halo_eps is meaningless. Get a representative MPI_halo_eps for MPI proc identification
fullMesh = .FALSE.
IF (halo_eps.LE.0.) THEN
  ! reconstruct halo_eps_velo
  IF (halo_eps_velo.EQ.0.) THEN
    MPI_halo_eps_velo = c
  ELSE
    MPI_halo_eps_velo = halo_eps_velo
  END IF

  ! reconstruct deltaT
  deltaT = 0.
  IF (ManualTimeStep.GT.0.) THEN
    deltaT    = ManualTimeStep
#if ! (USE_HDG)
  ELSE
    deltaT    = CalcTimeStep()
#endif
  END IF

  ! calculate halo_eps
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
  MPI_halo_eps = RK_c(2)
  DO iStage=2,nRKStages-1
    MPI_halo_eps = MAX(MPI_halo_eps,RK_c(iStage+1)-RK_c(iStage))
  END DO
  MPI_halo_eps = MAX(MPI_halo_eps,1.-RK_c(nRKStages))
  MPI_halo_eps = MPI_halo_eps*MPI_halo_eps_velo*deltaT*SafetyFactor
#else
  MPI_halo_eps = MPI_halo_eps_velo*deltaT*SafetyFactor
#endif

  vec(1)   = GEO%xmaxglob-GEO%xminglob
  vec(2)   = GEO%ymaxglob-GEO%yminglob
  vec(3)   = GEO%zmaxglob-GEO%zminglob
  MPI_halo_diag = VECNORM(vec)

  MPI_halo_eps_woshape = MPI_halo_eps
  ! Check whether MPI_halo_eps is smaller than shape function radius e.g. 'shape_function'
  IF(StringBeginsWith(DepositionType,'shape_function'))THEN
    IF(r_sf.LT.0.) CALL abort(__STAMP__,'Shape function radius is below zero; not correctly set yet? r_sf=',RealInfoOpt=r_sf)
    MPI_halo_eps = MPI_halo_eps + r_sf
    CALL PrintOption('MPI_halo_eps from shape function radius','CALCUL.',RealOpt=MPI_halo_eps)
  END IF

  ! compare halo_eps against global diagonal and reduce if necessary
  IF (.NOT.ALMOSTZERO(MPI_halo_eps).AND.(MPI_halo_diag.GE.MPI_halo_eps)) THEN
    SWRITE(UNIT_stdOUt,'(A,E11.3)') ' | No halo_eps given. Reconstructed to ',MPI_halo_eps
  ELSEIF (.NOT.ALMOSTZERO(MPI_halo_eps).AND.(MPI_halo_diag.LT.MPI_halo_eps)) THEN
    fullMesh = .TRUE.
    MPI_halo_eps = MPI_halo_diag
    SWRITE(UNIT_stdOUt,'(A,E11.3)') ' | No halo_eps given. Reconstructed to global diag with ',MPI_halo_eps
  ! halo_eps still at zero. Set it to global diagonal
  ELSE
    fullMesh = .TRUE.
    MPI_halo_eps = MPI_halo_diag
    SWRITE(UNIT_stdOUt,'(A,F11.3)') ' | No halo_eps given and could not be reconstructed. Using global diag with ',MPI_halo_eps
  END IF
ELSE
  vec(1)   = GEO%xmaxglob-GEO%xminglob
  vec(2)   = GEO%ymaxglob-GEO%yminglob
  vec(3)   = GEO%zmaxglob-GEO%zminglob
  MPI_halo_diag = VECNORM(vec)

  IF (MPI_halo_diag.LE.halo_eps) fullMesh = .TRUE.
  MPI_halo_eps = halo_eps
  MPI_halo_eps_woshape = halo_eps_woshape
END IF

! Identify all procs with elements in range
xCoordsProc(1) = GEO%xmin
xCoordsProc(2) = GEO%xmax
xCoordsProc(3) = GEO%ymin
xCoordsProc(4) = GEO%ymax
xCoordsProc(5) = GEO%zmin
xCoordsProc(6) = GEO%zmax

! Use a named loop so the entire element can be cycled
ElemLoop:  DO iElem = 1,nComputeNodeTotalElems
  ElemID    = GetGlobalElemID(iElem)
  localElem = ElemID - offSetElem
  HaloProc  = ElemInfo_Shared(ELEM_RANK,ElemID)

  BoundsOfElemCenter(1:3) = (/SUM(      BoundsOfElem_Shared(1:2,1,ElemID)), &
                              SUM(      BoundsOfElem_Shared(1:2,2,ElemID)), &
                              SUM(      BoundsOfElem_Shared(1:2,3,ElemID)) /) / 2.
  BoundsOfElemCenter(4)   = VECNORM ((/ BoundsOfElem_Shared(2  ,1,ElemID)-BoundsOfElem_Shared(1,1,ElemID), &
                                        BoundsOfElem_Shared(2  ,2,ElemID)-BoundsOfElem_Shared(1,2,ElemID), &
                                        BoundsOfElem_Shared(2  ,3,ElemID)-BoundsOfElem_Shared(1,3,ElemID) /) / 2.)

#if (PP_TimeDiscMethod==400)
  IF (HaloProc.EQ.myRank) THEN
    DO iSide = 1, nExchangeSides
      ! compare distance of centers with sum of element outer radii+halo_eps
      IF (VECNORM(BoundsOfElemCenter(1:3)-MPISideBoundsOfNbElemCenter(1:3,iSide)) &
          .GT. MPI_halo_eps_woshape+BoundsOfElemCenter(4)+MPISideBoundsOfNbElemCenter(4,iSide)) THEN

        ! Also check periodic directions. Only MPI sides of the local proc are
        ! taken into account, so do not perform additional case distinction
        SELECT CASE(GEO%nPeriodicVectors)
          ! One periodic vector
          CASE(1)
            DO iPeriodicDir = 1,2
              IF (VECNORM( BoundsOfElemCenter(1:3)                                                       &
                         + GEO%PeriodicVectors(1:3,1) * DirPeriodicVector(iPeriodicDir)                  &
                         - MPISideBoundsOfNbElemCenter(1:3,iSide))                                       &
                      .LE. MPI_halo_eps_woshape+BoundsOfElemCenter(4)                                    & !-BoundsOfElemCenter(5) &
                         + MPISideBoundsOfNbElemCenter(4,iSide) ) THEN
                ! flag the proc as exchange proc (in halo region)
                IsExchangeElem(localElem) = .TRUE.
                CYCLE ElemLoop
              END IF
            END DO

          ! Two periodic vectors. Also check linear combination, see particle_bgm.f90
          CASE(2)
            DO iPeriodicVector = 1,2
              DO iPeriodicDir = 1,2
                ! check if element is within halo_eps of periodically displaced element
                IF (VECNORM( BoundsOfElemCenter(1:3)                                                     &
                           + GEO%PeriodicVectors(1:3,iPeriodicVector) * DirPeriodicVector(iPeriodicDir)  &
                           - MPISideBoundsOfNbElemCenter(1:3,iSide))                                     &
                        .LE. MPI_halo_eps_woshape+BoundsOfElemCenter(4)                                  & !-BoundsOfElemCenter(5) &
                           + MPISideBoundsOfNbElemCenter(4,iSide)) THEN
                  ! flag the proc as exchange proc (in halo region)
                  IsExchangeElem(localElem) = .TRUE.
                  CYCLE ElemLoop
                END IF
              END DO
            END DO

            DO iPeriodicDir = 1,2
              DO jPeriodicDir = 1,2
                  ! check if element is within halo_eps of periodically displaced element
                  IF (VECNORM( BoundsOfElemCenter(1:3)                                                   &
                             + GEO%PeriodicVectors(1:3,1) * DirPeriodicVector(iPeriodicDir)              &
                             + GEO%PeriodicVectors(1:3,2) * DirPeriodicVector(jPeriodicDir)              &
                             - MPISideBoundsOfNbElemCenter(1:3,iSide))                                   &
                          .LE. MPI_halo_eps_woshape+BoundsOfElemCenter(4)                                & !-BoundsOfElemCenter(5) &
                             + MPISideBoundsOfNbElemCenter(4,iSide)) THEN
                    ! flag the proc as exchange proc (in halo region)
                    IsExchangeElem(localElem) = .TRUE.
                    CYCLE ElemLoop
                  END IF
                END DO
              END DO

          ! Three periodic vectors. Also check linear combination, see particle_bgm.f90
          CASE(3)
            ! check the three periodic vectors. Begin with checking the first periodic vector, followed by the combination of
            ! the first periodic vector with the others. Then check the other combinations, i.e. 1, 1+2, 1+3, 2, 2+3, 3, 1+2+3
            DO iPeriodicVector = 1,3
              DO iPeriodicDir = 1,2
                ! element might be already added back
                ! check if element is within halo_eps of periodically displaced element
                IF (VECNORM( BoundsOfElemCenter(1:3)                                                     &
                           + GEO%PeriodicVectors(1:3,iPeriodicVector) * DirPeriodicVector(iPeriodicDir)  &
                           - MPISideBoundsOfNbElemCenter(1:3,iSide))                                     &
                        .LE. MPI_halo_eps_woshape+BoundsOfElemCenter(4)                                  & !-BoundsOfElemCenter(5) &
                           + MPISideBoundsOfNbElemCenter(4,iSide)) THEN
                  ! flag the proc as exchange proc (in halo region)
                  IsExchangeElem(localElem) = .TRUE.
                  CYCLE ElemLoop
                END IF

                DO jPeriodicVector = 1,3
                  DO jPeriodicDir = 1,2
                    IF (iPeriodicVector.GE.jPeriodicVector) CYCLE

                    ! check if element is within halo_eps of periodically displaced element
                    IF (VECNORM( BoundsOfElemCenter(1:3)                                                      &
                               + GEO%PeriodicVectors(1:3,iPeriodicVector) * DirPeriodicVector(iPeriodicDir)   &
                               + GEO%PeriodicVectors(1:3,jPeriodicVector) * DirPeriodicVector(jPeriodicDir)   &
                               - MPISideBoundsOfNbElemCenter(1:3,iSide))                                      &
                            .LE. MPI_halo_eps_woshape+BoundsOfElemCenter(4)                                   & !-BoundsOfElemCenter(5) &
                               + MPISideBoundsOfNbElemCenter(4,iSide)) THEN
                      ! flag the proc as exchange proc (in halo region)
                      IsExchangeElem(localElem) = .TRUE.
                      CYCLE ElemLoop
                    END IF
                  END DO
                END DO
              END DO
            END DO

            ! check if element is within halo_eps of periodically displaced element
            DO iPeriodicDir = 1,2
              DO jPeriodicDir = 1,2
                DO kPeriodicDir = 1,2
                  IF (VECNORM( BoundsOfElemCenter(1:3)                                                   &
                             + GEO%PeriodicVectors(1:3,1) * DirPeriodicVector(iPeriodicDir)              &
                             + GEO%PeriodicVectors(1:3,2) * DirPeriodicVector(jPeriodicDir)              &
                             + GEO%PeriodicVectors(1:3,3) * DirPeriodicVector(kPeriodicDir)              &
                             - MPISideBoundsOfNbElemCenter(1:3,iSide))                                   &
                          .LE. MPI_halo_eps_woshape+BoundsOfElemCenter(4)                                & !-BoundsOfElemCenter(5) &
                             + MPISideBoundsOfNbElemCenter(4,iSide) ) THEN
                    ! flag the proc as exchange proc (in halo region)
                    IsExchangeElem(localElem) = .TRUE.
                    CYCLE ElemLoop
                  END IF
                END DO
              END DO
            END DO
          ! No periodic vectors, element out of range
          CASE(0)
              ! Do nothing
          CASE DEFAULT
            CALL ABORT(__STAMP__,'Invalid number of periodic vectors in particle_mpi_halo.f90')
        END SELECT

        ! Check rot periodic Elems and if iSide is on rot periodic BC
        IF(MeshHasRotPeriodic) THEN
          DO iPeriodicDir = 1,2
            ASSOCIATE( alpha => GEO%RotPeriodicAngle * DirPeriodicVector(iPeriodicDir) )
              SELECT CASE(GEO%RotPeriodicAxi)
                CASE(1) ! x-rotation axis
                  RotBoundsOfElemCenter(1) = BoundsOfElemCenter(1)
                  RotBoundsOfElemCenter(2) = COS(alpha)*BoundsOfElemCenter(2) - SIN(alpha)*BoundsOfElemCenter(3)
                  RotBoundsOfElemCenter(3) = SIN(alpha)*BoundsOfElemCenter(2) + COS(alpha)*BoundsOfElemCenter(3)
                CASE(2) ! y-rotation axis
                  RotBoundsOfElemCenter(1) = COS(alpha)*BoundsOfElemCenter(1) + SIN(alpha)*BoundsOfElemCenter(3)
                  RotBoundsOfElemCenter(2) = BoundsOfElemCenter(2)
                  RotBoundsOfElemCenter(3) =-SIN(alpha)*BoundsOfElemCenter(1) + COS(alpha)*BoundsOfElemCenter(3)
                CASE(3) ! z-rotation axis
                  RotBoundsOfElemCenter(1) = COS(alpha)*BoundsOfElemCenter(1) - SIN(alpha)*BoundsOfElemCenter(2)
                  RotBoundsOfElemCenter(2) = SIN(alpha)*BoundsOfElemCenter(1) + COS(alpha)*BoundsOfElemCenter(2)
                  RotBoundsOfElemCenter(3) = BoundsOfElemCenter(3)
              END SELECT
            END ASSOCIATE
            ! check if element is within halo_eps of rotationally displaced element
            IF (VECNORM( RotBoundsOfElemCenter(1:3)                                           &
                       - MPISideBoundsOfNbElemCenter(1:3,iSide))                              &
                    .LE. MPI_halo_eps_woshape+BoundsOfElemCenter(4)                           & !-BoundsOfElemCenter(5) &
                       + MPISideBoundsOfNbElemCenter(4,iSide) ) THEN
              ! flag the proc as exchange proc (in halo region)
              IsExchangeElem(localElem) = .TRUE.
              CYCLE ElemLoop
            END IF
          END DO
          ! End check rot periodic Elems and if iSide is on rot periodic BC
        END IF ! GEO%RotPeriodicBC

      ! Element is in range of not-periodically displaced MPI side
      ELSE
        IsExchangeElem(localElem) = .TRUE.
        CYCLE ElemLoop
      END IF
    END DO ! iSide = 1, nExchangeSides
    CYCLE ElemLoop
  END IF ! HaloProc.EQ.myRank
#else
  ! Skip own myrank
  IF (HaloProc.EQ.myRank) CYCLE
#endif /*(PP_TimeDiscMethod==400)*/

  ! Skip if the proc is already flagged, only if the exact elements are not required (.NOT.shape_function)
  IF(.NOT.StringBeginsWith(DepositionType,'shape_function'))THEN
    SELECT CASE(GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc))
      ! Proc not previously encountered, check if possibly in range
      CASE(-1)
        IF(GEO%RotPeriodicBC) THEN
          GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 0
        ELSE
          firstElem = offsetElemMPI(HaloProc)+1
          lastElem  = offsetElemMPI(HaloProc +1)

          SELECT CASE(TrackingMethod)
            ! Build mesh min/max on BezierControlPoints for possibly curved elements
            CASE(REFMAPPING,TRACING)
              firstSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,firstElem)+1
              lastSide  = ElemInfo_Shared(ELEM_LASTSIDEIND ,lastElem)

              xCoordsOrigin(1) = MINVAL(BezierControlPoints3D(1,:,:,firstSide:lastSide))
              xCoordsOrigin(2) = MAXVAL(BezierControlPoints3D(1,:,:,firstSide:lastSide))
              xCoordsOrigin(3) = MINVAL(BezierControlPoints3D(2,:,:,firstSide:lastSide))
              xCoordsOrigin(4) = MAXVAL(BezierControlPoints3D(2,:,:,firstSide:lastSide))
              xCoordsOrigin(5) = MINVAL(BezierControlPoints3D(3,:,:,firstSide:lastSide))
              xCoordsOrigin(6) = MAXVAL(BezierControlPoints3D(3,:,:,firstSide:lastSide))

            ! TriaTracking does not have curved elements, nodeCoords are sufficient
            CASE(TRIATRACKING)
              xCoordsOrigin(1) = MINVAL(NodeCoords_Shared(1,ElemInfo_Shared(ELEM_FIRSTNODEIND,firstElem) + 1 &
                                                           :ElemInfo_Shared(ELEM_LASTNODEIND ,lastElem)))
              xCoordsOrigin(2) = MAXVAL(NodeCoords_Shared(1,ElemInfo_Shared(ELEM_FIRSTNODEIND,firstElem) + 1 &
                                                           :ElemInfo_Shared(ELEM_LASTNODEIND ,lastElem)))
              xCoordsOrigin(3) = MINVAL(NodeCoords_Shared(2,ElemInfo_Shared(ELEM_FIRSTNODEIND,firstElem) + 1 &
                                                           :ElemInfo_Shared(ELEM_LASTNODEIND ,lastElem)))
              xCoordsOrigin(4) = MAXVAL(NodeCoords_Shared(2,ElemInfo_Shared(ELEM_FIRSTNODEIND,firstElem) + 1 &
                                                           :ElemInfo_Shared(ELEM_LASTNODEIND ,lastElem)))
              xCoordsOrigin(5) = MINVAL(NodeCoords_Shared(3,ElemInfo_Shared(ELEM_FIRSTNODEIND,firstElem) + 1 &
                                                           :ElemInfo_Shared(ELEM_LASTNODEIND ,lastElem)))
              xCoordsOrigin(6) = MAXVAL(NodeCoords_Shared(3,ElemInfo_Shared(ELEM_FIRSTNODEIND,firstElem) + 1 &
                                                           :ElemInfo_Shared(ELEM_LASTNODEIND ,lastElem)))
          END SELECT

          ! Keep direction to account for accuracy issues
          IF (myRank.LT.HaloProc) THEN
            ProcInRange = HaloBoxInProc(xCoordsOrigin,xCoordsProc,MPI_halo_eps,GEO%nPeriodicVectors,GEO%PeriodicVectors)
          ELSE
            ProcInRange = HaloBoxInProc(xCoordsProc,xCoordsOrigin,MPI_halo_eps,GEO%nPeriodicVectors,GEO%PeriodicVectors)
          END IF

          ! Check if proc is in range
          IF (.NOT.ProcInRange) THEN
            ! Proc definitely not in range
            GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = -2
            CYCLE
          ELSE
            ! Proc possible in range
            GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 0
          END IF
        END IF
      ! Proc definitely not in range or already flagged
      CASE(-2,1,2)
        CYCLE
    END SELECT
  END IF

  DO iSide = 1, nExchangeSides
    ! compare distance of centers with sum of element outer radii+halo_eps
    IF (VECNORM(BoundsOfElemCenter(1:3)-MPISideBoundsOfElemCenter(1:3,iSide)) &
        .GT. MPI_halo_eps+BoundsOfElemCenter(4)+MPISideBoundsOfElemCenter(4,iSide)) THEN

      ! Also check periodic directions. Only MPI sides of the local proc are
      ! taken into account, so do not perform additional case distinction
      SELECT CASE(GEO%nPeriodicVectors)
        ! One periodic vector
        CASE(1)
          DO iPeriodicDir = 1,2
            IF (VECNORM( BoundsOfElemCenter(1:3)                                                       &
                       + GEO%PeriodicVectors(1:3,1) * DirPeriodicVector(iPeriodicDir)                  &
                       - MPISideBoundsOfElemCenter(1:3,iSide))                                         &
                    .LE. MPI_halo_eps+BoundsOfElemCenter(4)                                            & !-BoundsOfElemCenter(5) &
                       + MPISideBoundsOfElemCenter(4,iSide) ) THEN
              ! flag the proc as exchange proc (in halo region)
              IF(StringBeginsWith(DepositionType,'shape_function'))THEN
                IF (ElemInfo_Shared(ELEM_HALOFLAG,ElemID).NE.4) FlagShapeElem(iElem) = .TRUE.
                IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc).EQ.2) CYCLE ElemLoop
              END IF
              GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
              GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
              nExchangeProcessors = nExchangeProcessors + 1
              CYCLE ElemLoop
            END IF
          END DO

        ! Two periodic vectors. Also check linear combination, see particle_bgm.f90
        CASE(2)
          DO iPeriodicVector = 1,2
            DO iPeriodicDir = 1,2
              ! check if element is within halo_eps of periodically displaced element
              IF (VECNORM( BoundsOfElemCenter(1:3)                                                    &
                         + GEO%PeriodicVectors(1:3,iPeriodicVector) * DirPeriodicVector(iPeriodicDir) &
                         - MPISideBoundsOfElemCenter(1:3,iSide))                                      &
                      .LE. MPI_halo_eps+BoundsOfElemCenter(4)                                         & !-BoundsOfElemCenter(5) &
                         + MPISideBoundsOfElemCenter(4,iSide)) THEN
                ! flag the proc as exchange proc (in halo region)
                IF(StringBeginsWith(DepositionType,'shape_function'))THEN
                  IF (ElemInfo_Shared(ELEM_HALOFLAG,ElemID).NE.4) FlagShapeElem(iElem) = .TRUE.
                  IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc).EQ.2) CYCLE ElemLoop
                END IF
                GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
                GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
                nExchangeProcessors = nExchangeProcessors + 1
                CYCLE ElemLoop
              END IF
            END DO
          END DO

          DO iPeriodicDir = 1,2
            DO jPeriodicDir = 1,2
                ! check if element is within halo_eps of periodically displaced element
                IF (VECNORM( BoundsOfElemCenter(1:3)                                                  &
                           + GEO%PeriodicVectors(1:3,1) * DirPeriodicVector(iPeriodicDir)             &
                           + GEO%PeriodicVectors(1:3,2) * DirPeriodicVector(jPeriodicDir)             &
                           - MPISideBoundsOfElemCenter(1:3,iSide))                                    &
                        .LE. MPI_halo_eps+BoundsOfElemCenter(4)                                       & !-BoundsOfElemCenter(5) &
                           + MPISideBoundsOfElemCenter(4,iSide)) THEN
                  ! flag the proc as exchange proc (in halo region)
                  IF(StringBeginsWith(DepositionType,'shape_function'))THEN
                    IF (ElemInfo_Shared(ELEM_HALOFLAG,ElemID).NE.4) FlagShapeElem(iElem) = .TRUE.
                    IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc).EQ.2) CYCLE ElemLoop
                  END IF
                  GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
                  GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
                  nExchangeProcessors = nExchangeProcessors + 1
                  CYCLE ElemLoop
                END IF
              END DO
            END DO

        ! Three periodic vectors. Also check linear combination, see particle_bgm.f90
        CASE(3)
          ! check the three periodic vectors. Begin with checking the first periodic vector, followed by the combination of
          ! the first periodic vector with the others. Then check the other combinations, i.e. 1, 1+2, 1+3, 2, 2+3, 3, 1+2+3
          DO iPeriodicVector = 1,3
            DO iPeriodicDir = 1,2
              ! element might be already added back
              ! check if element is within halo_eps of periodically displaced element
              IF (VECNORM( BoundsOfElemCenter(1:3)                                                      &
                         + GEO%PeriodicVectors(1:3,iPeriodicVector) * DirPeriodicVector(iPeriodicDir)   &
                         - MPISideBoundsOfElemCenter(1:3,iSide))                                        &
                      .LE. MPI_halo_eps+BoundsOfElemCenter(4)                                           & !-BoundsOfElemCenter(5) &
                         + MPISideBoundsOfElemCenter(4,iSide)) THEN
                ! flag the proc as exchange proc (in halo region)
                IF(StringBeginsWith(DepositionType,'shape_function'))THEN
                  IF (ElemInfo_Shared(ELEM_HALOFLAG,ElemID).NE.4) FlagShapeElem(iElem) = .TRUE.
                  IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc).EQ.2) CYCLE ElemLoop
                END IF
                GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
                GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
                nExchangeProcessors = nExchangeProcessors + 1
                CYCLE ElemLoop
              END IF

              DO jPeriodicVector = 1,3
                DO jPeriodicDir = 1,2
                  IF (iPeriodicVector.GE.jPeriodicVector) CYCLE

                  ! check if element is within halo_eps of periodically displaced element
                  IF (VECNORM( BoundsOfElemCenter(1:3)                                                      &
                             + GEO%PeriodicVectors(1:3,iPeriodicVector) * DirPeriodicVector(iPeriodicDir)   &
                             + GEO%PeriodicVectors(1:3,jPeriodicVector) * DirPeriodicVector(jPeriodicDir)   &
                             - MPISideBoundsOfElemCenter(1:3,iSide))                                        &
                          .LE. MPI_halo_eps+BoundsOfElemCenter(4)                                           & !-BoundsOfElemCenter(5) &
                             + MPISideBoundsOfElemCenter(4,iSide)) THEN
                    ! flag the proc as exchange proc (in halo region)
                    IF(StringBeginsWith(DepositionType,'shape_function'))THEN
                      IF (ElemInfo_Shared(ELEM_HALOFLAG,ElemID).NE.4) FlagShapeElem(iElem) = .TRUE.
                      IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc).EQ.2) CYCLE ElemLoop
                    END IF
                    GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
                    GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
                    nExchangeProcessors = nExchangeProcessors + 1
                    CYCLE ElemLoop
                  END IF
                END DO
              END DO
            END DO
          END DO

          ! check if element is within halo_eps of periodically displaced element
          DO iPeriodicDir = 1,2
            DO jPeriodicDir = 1,2
              DO kPeriodicDir = 1,2
                IF (VECNORM( BoundsOfElemCenter(1:3)                                                        &
                           + GEO%PeriodicVectors(1:3,1) * DirPeriodicVector(iPeriodicDir)                   &
                           + GEO%PeriodicVectors(1:3,2) * DirPeriodicVector(jPeriodicDir)                   &
                           + GEO%PeriodicVectors(1:3,3) * DirPeriodicVector(kPeriodicDir)                   &
                           - MPISideBoundsOfElemCenter(1:3,iSide))                                          &
                        .LE. MPI_halo_eps+BoundsOfElemCenter(4)                                             & !-BoundsOfElemCenter(5) &
                           + MPISideBoundsOfElemCenter(4,iSide) ) THEN
                  ! flag the proc as exchange proc (in halo region)
                  IF(StringBeginsWith(DepositionType,'shape_function'))THEN
                    IF (ElemInfo_Shared(ELEM_HALOFLAG,ElemID).NE.4) FlagShapeElem(iElem) = .TRUE.
                    IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc).EQ.2) CYCLE ElemLoop
                  END IF
                  GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
                  GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
                  nExchangeProcessors = nExchangeProcessors + 1
                  CYCLE ElemLoop
                END IF
              END DO
            END DO
          END DO
        ! No periodic vectors, element out of range
        CASE(0)
            ! Do nothing
        CASE DEFAULT
          CALL ABORT(__STAMP__,'Invalid number of periodic vectors in particle_mpi_halo.f90')
      END SELECT

      ! Check rot periodic Elems and if iSide is on rot periodic BC
      IF(MeshHasRotPeriodic) THEN
        DO iPeriodicDir = 1,2
          ASSOCIATE( alpha => GEO%RotPeriodicAngle * DirPeriodicVector(iPeriodicDir) )
            SELECT CASE(GEO%RotPeriodicAxi)
              CASE(1) ! x-rotation axis
                RotBoundsOfElemCenter(1) = BoundsOfElemCenter(1)
                RotBoundsOfElemCenter(2) = COS(alpha)*BoundsOfElemCenter(2) - SIN(alpha)*BoundsOfElemCenter(3)
                RotBoundsOfElemCenter(3) = SIN(alpha)*BoundsOfElemCenter(2) + COS(alpha)*BoundsOfElemCenter(3)
              CASE(2) ! y-rotation axis
                RotBoundsOfElemCenter(1) = COS(alpha)*BoundsOfElemCenter(1) + SIN(alpha)*BoundsOfElemCenter(3)
                RotBoundsOfElemCenter(2) = BoundsOfElemCenter(2)
                RotBoundsOfElemCenter(3) =-SIN(alpha)*BoundsOfElemCenter(1) + COS(alpha)*BoundsOfElemCenter(3)
              CASE(3) ! z-rotation axis
                RotBoundsOfElemCenter(1) = COS(alpha)*BoundsOfElemCenter(1) - SIN(alpha)*BoundsOfElemCenter(2)
                RotBoundsOfElemCenter(2) = SIN(alpha)*BoundsOfElemCenter(1) + COS(alpha)*BoundsOfElemCenter(2)
                RotBoundsOfElemCenter(3) = BoundsOfElemCenter(3)
            END SELECT
          END ASSOCIATE
          ! check if element is within halo_eps of rotationally displaced element
          IF (VECNORM( RotBoundsOfElemCenter(1:3)                                           &
                     - MPISideBoundsOfElemCenter(1:3,iSide))                                &
                  .LE. MPI_halo_eps+BoundsOfElemCenter(4)                                   & !-BoundsOfElemCenter(5) &
                     + MPISideBoundsOfElemCenter(4,iSide) ) THEN
            ! flag the proc as exchange proc (in halo region)
            IF(StringBeginsWith(DepositionType,'shape_function'))THEN
              IF (ElemInfo_Shared(ELEM_HALOFLAG,ElemID).NE.4) FlagShapeElem(iElem) = .TRUE.
              IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc).EQ.2) CYCLE ElemLoop
            END IF
            GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
            GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
            nExchangeProcessors = nExchangeProcessors + 1
            CYCLE ElemLoop
          END IF
        END DO
        ! End check rot periodic Elems and if iSide is on rot periodic BC
      END IF ! GEO%RotPeriodicBC

    ! Element is in range of not-periodically displaced MPI side
    ELSE
      IF(StringBeginsWith(DepositionType,'shape_function'))THEN
        IF (ElemInfo_Shared(ELEM_HALOFLAG,ElemID).NE.4) FlagShapeElem(iElem) = .TRUE.
        IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc).EQ.2) CYCLE ElemLoop
      END IF
      GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
      GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
      nExchangeProcessors = nExchangeProcessors + 1
      CYCLE ElemLoop
    END IF
  END DO ! iSide = 1, nExchangeSides
END DO ElemLoop

#if (PP_TimeDiscMethod==400)
SDEALLOCATE(MPISideBoundsOfNbElemCenter)
#endif /*(PP_TimeDiscMethod==400)*/

! Remove elements if the halo proc contains only internal elements, i.e. we cannot possibly reach the halo element
!
!   CN1     CN2    > If a processor contains large changes in element size, internal elements might intersect with
!  _ _ _    _ _ _  > the MPI sides. Since a processor checks all potential halo elements against only its own exchange
! |_|_|_|  |_|_|_| > sides, the large elements will only take effect for the proc not containing it. However, if the
! |_|_|_|  |_| | | > proc flags only large internal elements without flagging a single exchange element, there is no
! |_|_|_|  |_|_|_| > way for a particle to actually reach.
! |_|_|_|  |_| | |
! |_|_|_|  |_|_|_| > This routine therefore checks for the presence of exchange sides on the procs and unflags the
!                  > proc if none is found.
!
DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE

  IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc).LE.0) CYCLE

  ProcHasExchangeElem = .FALSE.
  ! Use a named loop so the entire element can be cycled
ExchangeLoop: DO iElem = offsetElemMPI(iProc)+1,offsetElemMPI(iProc+1)
    ! Ignore elements outside nComputeNodeTotalElems
    IF (ElemInfo_Shared(ELEM_HALOFLAG,iElem).LT.0) CYCLE

    DO iLocSide = 1,6
      SideID   = GetGlobalNonUniqueSideID(iElem,iLocSide)
      NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)

      ! Mortar side
      IF (NbElemID.LT.0) THEN
        nMortarElems = MERGE(4,2,SideInfo_Shared(SIDE_NBELEMID,SideID).EQ.-1)

        DO iMortar = 1,nMortarElems
          NbSideID = -SideInfo_Shared(SIDE_LOCALID,SideID + iMortar)
          ! If small mortar side not defined, skip it for now, likely not inside the halo region
          IF (NbSideID.LT.1) CYCLE

          NbElemID = SideInfo_Shared(SIDE_NBELEMID,SideID + iMortar)
          ! If small mortar element not defined, skip it for now, likely not inside the halo region
          IF (NbElemID.LT.1) CYCLE

          ! If any of the small mortar sides is not on the halo proc, the side is a MPI side
          IF (NbElemID.LT.offsetElemMPI(iProc)+1 .OR. NbElemID.GT.offsetElemMPI(iProc+1)) THEN
            ProcHasExchangeElem = .TRUE.
            EXIT ExchangeLoop
          END IF
        END DO

      ! regular side or small mortar side
      ELSE
        ! Only check inner (MPI interfaces) and boundary sides (NbElemID.EQ.0)
        ! Boundary sides cannot be discarded because one proc might have MPI interfaces and the other might not
        IF (NbElemID.LT.offsetElemMPI(iProc)+1 .OR. NbElemID.GT.offsetElemMPI(iProc+1)) THEN
          ProcHasExchangeElem = .TRUE.
          EXIT ExchangeLoop
        END IF
      END IF
    END DO
  END DO ExchangeLoop ! iElem = offsetElemMPI(iProc)+1,offsetElemMPI(iProc+1)

  ! Processor has halo elements but no MPI sides, remove the exchange processor
  IF (.NOT.ProcHasExchangeElem) THEN
    GlobalProcToExchangeProc(:,iProc) = 0
    nExchangeProcessors = nExchangeProcessors - 1
  END IF
END DO ! iProc = 1,nExchangeProcessors

! Notify every proc if it was identified by the local proc
IF(CheckExchangeProcs)THEN
  DO iProc = 0,nProcessors_Global-1
    IF (iProc.EQ.myRank) CYCLE

    ! CommFlag holds the information if the local proc wants to communicate with iProc. Cannot be a logical because ISEND might not
    ! return before the next value is written
    CommFlag(iProc) = MERGE(.TRUE.,.FALSE.,GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc).GT.0)
    CALL MPI_ISEND( CommFlag(iProc)              &
                  , 1                            &
                  , MPI_LOGICAL                  &
                  , iProc                        &
                  , 1999                         &
                  , MPI_COMM_WORLD               &
                  , SendRequest(iProc)           &
                  , IERROR)
  END DO

  ! Finish communication
  DO iProc = 0,nProcessors_Global-1
    IF (iProc.EQ.myRank) CYCLE

    CALL MPI_WAIT(RecvRequest(iProc),MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
    CALL MPI_WAIT(SendRequest(iProc),MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END DO

  ! Append previously not found procs to list of exchange processors
  nNonSymmetricExchangeProcs = 0
  ! Flag elements with processor IDs, which we cannot reach (but they reach us)
  ALLOCATE(myInvisibleRank(1:nElems))
  myInvisibleRank=-1

  DO iProc = 0,nProcessors_Global-1
    IF (iProc.EQ.myRank) CYCLE

    ! Ignore procs that are already flagged or not requesting communication
    IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc).GT.0) CYCLE
    ! GlobalProcToRecvProc(iProc) is true if iProc has flagged myrank
    IF (.NOT.GlobalProcToRecvProc(iProc)) CYCLE

    ! Found a previously missing proc
    GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc) = 2
    GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc) = nExchangeProcessors
    nNonSymmetricExchangeProcs = nNonSymmetricExchangeProcs + 1
    nExchangeProcessors        = nExchangeProcessors + 1
    DO iElem=1,nElems
      myInvisibleRank(iElem)=iProc
    END DO ! iElem=1,nElems
  END DO

  ! On smooth grids, nNonSymmetricExchangeProcs should be zero. Only output if previously missing particle exchange procs are found
  CALL MPI_ALLREDUCE(nNonSymmetricExchangeProcs, nNonSymmetricExchangeProcsGlob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERROR)
  ! Check sum of nNonSymmetricExchangeProcs over all processors
  IF(nNonSymmetricExchangeProcsGlob.GT.0)THEN
    SWRITE(UNIT_StdOut,'(X,131("~"))')
    SWRITE(Unit_StdOut,'(A,I0,A)') ' | Found ',nNonSymmetricExchangeProcsGlob, ' previously missing non-symmetric particle exchange procs'
    SWRITE(Unit_StdOut,'(A)')" | See ElemData container 'myInvisibleRank' for more information on which MPI ranks are non-symmetric"
    SWRITE(Unit_StdOut,'(A)')" | This information is written to "//TRIM(ProjectName)//"_MyInvisibleRank.h5 (only when CheckExchangeProcs=T)"
    SWRITE(Unit_StdOut,'(A)')" | This check is optional. You can disable it via CheckExchangeProcs = F"
    SWRITE(UNIT_StdOut,'(X,131("~"))')
    CALL AddToElemData(ElementOut,'myInvisibleRank',LongIntArray=myInvisibleRank)
    CALL WriteMyInvisibleRankToHDF5()
    ! Only root aborts
    IF(AbortExchangeProcs) CALL CollectiveStop(__STAMP__," Non-symmetric particle exchange procs > 0. This abort is optional."//&
        " You can disable it via AbortExchangeProcs = F.\n          See message above regarding 'myInvisibleRank' output for details.")
  ELSE
    SDEALLOCATE(myInvisibleRank)
  END IF ! nNonSymmetricExchangeProcsGlob.GT.0
END IF

SDEALLOCATE(GlobalProcToRecvProc)
SDEALLOCATE(RecvRequest)
SDEALLOCATE(SendRequest)
SDEALLOCATE(CommFlag)

! Build reverse mapping
!-- EXCHANGE_PROC_TYPE information is currently unused and either -1 (no communication) or 2 (communication). Can be used to
!-- implement check if exchange partner is on the same compute node, so build it here
ALLOCATE(ExchangeProcToGlobalProc(2,0:nExchangeProcessors-1))

DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE

  IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc).GT.0) THEN
    ! Find it the other proc is on the same compute node
    ExchangeProcLeader = INT(GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc)/nComputeNodeProcessors)
    IF (ExchangeProcLeader.EQ.myLeaderGroupRank) THEN
      GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc) = 1
      ExchangeProcToGlobalProc(EXCHANGE_PROC_TYPE,GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc)) = 1
      ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc)) = iProc
    ELSE
      ExchangeProcToGlobalProc(EXCHANGE_PROC_TYPE,GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc)) = 2
      ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc)) = iProc
    END IF
  END IF
END DO

! -- Average number of exchange processors
CALL MPI_REDUCE(nExchangeProcessors,nExchangeProcessorsGlobal,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
SWRITE(UNIT_stdOut,'(A,I0,A)') ' | Started particle exchange communication with average ', &
                                 nExchangeProcessorsGlobal/nProcessors_Global            , &
                                 ' partners per proc'

! Build shape function mapping
IF(StringBeginsWith(DepositionType,'shape_function'))THEN
  CALL Allocate_Shared((/nComputeNodeTotalElems,nProcessors/),ShapeElemProcSend_Shared_Win,ShapeElemProcSend_Shared)
  CALL MPI_WIN_LOCK_ALL(0,ShapeElemProcSend_Shared_Win,IERROR)
  IF (myComputeNodeRank.EQ.0)  ShapeElemProcSend_Shared = .FALSE.
  CALL BARRIER_AND_SYNC(ShapeElemProcSend_Shared_Win,MPI_COMM_SHARED)

  ! 1st loop to determine the number of nSendShapeElems
  nSendShapeElems = 0
  DO iELem = 1,nComputeNodeTotalElems
    IF (FlagShapeElem(iElem)) nSendShapeElems = nSendShapeElems + 1
  END DO

  ALLOCATE(SendShapeElemID(1:nSendShapeElems), SendElemShapeID(1:nComputeNodeTotalElems), CNRankToSendRank(0:nProcessors-1))
  ! Initialize
  SendShapeElemID = -1
  SendElemShapeID = -1
  CNRankToSendRank = -1

  ! 2nd loop to fill the array of size nSendShapeElems
  nSendShapeElems = 0
  DO iELem = 1,nComputeNodeTotalElems
    IF (FlagShapeElem(iElem)) THEN
      nSendShapeElems                  = nSendShapeElems + 1
      SendShapeElemID(nSendShapeElems) = iElem
!      SendElemShapeID(iElem)           = nSendShapeElems
    END IF
  END DO

  ! 1 of 2: Inner-Node Communication
  IF (myComputeNodeRank.EQ.0) THEN
    ALLOCATE(ShapeMapping(1:nComputeNodeProcessors-1), &
             RecvRequest (1:nComputeNodeProcessors-1))

    DO iProc = 1,nComputeNodeProcessors-1
        CALL MPI_IRECV( ShapeMapping(iProc)%nRecvShapeElems &
                      , 1                                   &
                      , MPI_INTEGER                         &
                      , iProc                               &
                      , 2001                                &
                      , MPI_COMM_SHARED                     &
                      , RecvRequest(iProc)                  &
                      , IERROR)
    END DO

    DO iProc = 1,nComputeNodeProcessors-1
      CALL MPI_WAIT(RecvRequest(iProc),MPIStatus,IERROR)
      IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
    END DO

    DO iProc = 1,nComputeNodeProcessors-1
      IF (ShapeMapping(iProc)%nRecvShapeElems.EQ.0) CYCLE
      ALLOCATE(ShapeMapping(iProc)%RecvShapeElemID(1:ShapeMapping(iProc)%nRecvShapeElems))

      CALL MPI_IRECV( ShapeMapping(iProc)%RecvShapeElemID   &
                    , ShapeMapping(iProc)%nRecvShapeElems   &
                    , MPI_INTEGER                           &
                    , iProc                                 &
                    , 2001                                  &
                    , MPI_COMM_SHARED                       &
                    , RecvRequest(iProc)                    &
                    , IERROR)
    END DO

    DO iProc = 1,nComputeNodeProcessors-1
      IF (ShapeMapping(iProc)%nRecvShapeElems.EQ.0) CYCLE
      CALL MPI_WAIT(RecvRequest(iProc),MPIStatus,IERROR)
      IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
    END DO

    DO iProc = 1, nComputeNodeProcessors -1
      DO iElem = 1, ShapeMapping(iProc)%nRecvShapeElems
        ShapeElemProcSend_Shared(ShapeMapping(iProc)%RecvShapeElemID(iElem), iProc+1+ComputeNodeRootRank) = .TRUE.
      END DO
    END DO
    ! Own contribution
    DO iELem = 1,nSendShapeElems
      ShapeElemProcSend_Shared(SendShapeElemID(iElem), 1+ComputeNodeRootRank) = .TRUE.
    END DO

    SDEALLOCATE(ShapeMapping)
    SDEALLOCATE(RecvRequest)
    ! 2 of 2: Multi-node communication
    ! Second stage of communication, identify and send inter-compute-node information
    IF (nLeaderGroupProcs.GT.1) THEN
      CALL MPI_WIN_SYNC(ShapeElemProcSend_Shared_Win,iError)
      ALLOCATE(CNShapeMapping(0:nLeaderGroupProcs-1), &
               SendRequest   (0:nLeaderGroupProcs-1), &
               RecvRequest   (0:nLeaderGroupProcs-1))

      DO iProc = 0,nLeaderGroupProcs-1
        IF (iProc.EQ.myLeaderGroupRank) CYCLE

        CALL MPI_IRECV( CNShapeMapping(iProc)%nRecvShapeElems(1:2)   &
                      , 2                                       &
                      , MPI_INTEGER                             &
                      , iProc                                   &
                      , 2002                                    &
                      , MPI_COMM_LEADERS_SHARED                 &
                      , RecvRequest(iProc)                      &
                      , IERROR)
      END DO

      CNShapeMapping%nSendShapeElems(1) = 0
      CNShapeMapping%nSendShapeElems(2) = 0
      ! Count number of elems per CN
      DO iElem = nComputeNodeElems+1,nComputeNodeTotalElems
        GlobalElemID = GetGlobalElemID(iElem)
        IF ((ElemInfo_Shared(ELEM_HALOFLAG,GlobalElemID).GE.2).AND.(ElemInfo_Shared(ELEM_HALOFLAG,GlobalElemID).NE.4)) THEN
          GlobalElemRank = ElemInfo_Shared(ELEM_RANK,GlobalElemID)
          GlobalLeaderRank = INT(GlobalElemRank/nComputeNodeProcessors)
          CNShapeMapping(GlobalLeaderRank)%nSendShapeElems(1) = CNShapeMapping(GlobalLeaderRank)%nSendShapeElems(1) + 1
          CNShapeMapping(GlobalLeaderRank)%nSendShapeElems(2) = MAX(CNShapeMapping(GlobalLeaderRank)%nSendShapeElems(2), COUNT(ShapeElemProcSend_Shared(iElem, :)))
        END IF
      END DO

      DO iProc = 0,nLeaderGroupProcs-1
        IF (iProc.EQ.myLeaderGroupRank) CYCLE

        CALL MPI_ISEND( CNShapeMapping(iProc)%nSendShapeElems(1:2)   &
                      , 2                                       &
                      , MPI_INTEGER                             &
                      , iProc                                   &
                      , 2002                                    &
                      , MPI_COMM_LEADERS_SHARED                 &
                      , SendRequest(iProc)                     &
                      , IERROR)
      END DO
      DO iProc = 0,nLeaderGroupProcs-1
        IF (iProc.EQ.myLeaderGroupRank) CYCLE

        CALL MPI_WAIT(RecvRequest(iProc),MPIStatus,IERROR)
        IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
        CALL MPI_WAIT(SendRequest(iProc),MPIStatus,IERROR)
        IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
      END DO

      DO iProc = 0,nLeaderGroupProcs-1
        IF (iProc.EQ.myLeaderGroupRank) CYCLE

        IF (CNShapeMapping(iProc)%nRecvShapeElems(1).EQ.0) CYCLE

        ALLOCATE(CNShapeMapping(iProc)%RecvShapeElemID(CNShapeMapping(iProc)%nRecvShapeElems(1)))
        ALLOCATE(CNShapeMapping(iProc)%RecvShapeProcElemID(CNShapeMapping(iProc)%nRecvShapeElems(1),CNShapeMapping(iProc)%nRecvShapeElems(2)))
!        ALLOCATE(CNShapeMapping(iProc)%RecvBuffer(1:4,0:PP_N,0:PP_N,0:PP_N, 1:CNShapeMapping(iProc)%nRecvShapeElems))

        CALL MPI_IRECV( CNShapeMapping(iProc)%RecvShapeElemID   &
                      , CNShapeMapping(iProc)%nRecvShapeElems(1)   &
                      , MPI_INTEGER                             &
                      , iProc                                   &
                      , 2002                                    &
                      , MPI_COMM_LEADERS_SHARED                 &
                      , RecvRequest(iProc)                      &
                      , IERROR)

      END DO

      DO iProc = 0,nLeaderGroupProcs-1
        IF (iProc.EQ.myLeaderGroupRank) CYCLE
        IF (CNShapeMapping(iProc)%nSendShapeElems(1).EQ.0) CYCLE
        ALLOCATE(CNShapeMapping(iProc)%SendShapeElemID(CNShapeMapping(iProc)%nSendShapeElems(1)))
        ALLOCATE(CNShapeMapping(iProc)%SendShapeProcElemID(CNShapeMapping(iProc)%nSendShapeElems(1),CNShapeMapping(iProc)%nSendShapeElems(2)))
        CNShapeMapping(iProc)%nSendShapeElems(1) = 0
      END DO


      DO iElem = nComputeNodeElems+1,nComputeNodeTotalElems
        GlobalElemID = GetGlobalElemID(iElem)
        IF ((ElemInfo_Shared(ELEM_HALOFLAG,GlobalElemID).GE.2).AND.(ElemInfo_Shared(ELEM_HALOFLAG,GlobalElemID).NE.4)) THEN
          GlobalElemRank = ElemInfo_Shared(ELEM_RANK,GlobalElemID)
          GlobalLeaderRank = INT(GlobalElemRank/nComputeNodeProcessors)

          CNShapeMapping(GlobalLeaderRank)%nSendShapeElems(1) = CNShapeMapping(GlobalLeaderRank)%nSendShapeElems(1) + 1
          CNShapeMapping(GlobalLeaderRank)%SendShapeElemID(CNShapeMapping(GlobalLeaderRank)%nSendShapeElems(1)) = GlobalElemID
          CNShapeMapping(GlobalLeaderRank)%SendShapeProcElemID(CNShapeMapping(GlobalLeaderRank)%nSendShapeElems(1),:) = 0
          jProc = 0
          DO iProc = 0, nComputeNodeProcessors-1
            IF (ShapeElemProcSend_Shared(iElem, iProc+1+ComputeNodeRootRank)) THEN
              jProc = jProc + 1
              CNShapeMapping(GlobalLeaderRank)%SendShapeProcElemID(CNShapeMapping(GlobalLeaderRank)%nSendShapeElems(1),jProc) = iProc+ComputeNodeRootRank+1
            END IF
          END DO
        END IF
      END DO

      DO iProc = 0,nLeaderGroupProcs-1
        IF (iProc.EQ.myLeaderGroupRank) CYCLE

        IF (CNShapeMapping(iProc)%nSendShapeElems(1).EQ.0) CYCLE

        CALL MPI_ISEND( CNShapeMapping(iProc)%SendShapeElemID   &
                      , CNShapeMapping(iProc)%nSendShapeElems(1)   &
                      , MPI_INTEGER                             &
                      , iProc                                   &
                      , 2002                                    &
                      , MPI_COMM_LEADERS_SHARED                 &
                      , SendRequest(iProc)                      &
                      , IERROR)
      END DO

      DO iProc = 0,nLeaderGroupProcs-1
        IF (iProc.EQ.myLeaderGroupRank) CYCLE

        IF (CNShapeMapping(iProc)%nRecvShapeElems(1).NE.0) THEN
          CALL MPI_WAIT(RecvRequest(iProc),MPIStatus,IERROR)
          IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
        END IF

        IF (CNShapeMapping(iProc)%nSendShapeElems(1).NE.0) THEN
          CALL MPI_WAIT(SendRequest(iProc),MPIStatus,IERROR)
          IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
        END IF

        IF (CNShapeMapping(iProc)%nRecvShapeElems(1).GT.0) THEN
          CALL MPI_IRECV( CNShapeMapping(iProc)%RecvShapeProcElemID   &
                , CNShapeMapping(iProc)%nRecvShapeElems(1)*CNShapeMapping(iProc)%nRecvShapeElems(2)   &
                , MPI_INTEGER                             &
                , iProc                                   &
                , 2012                                    &
                , MPI_COMM_LEADERS_SHARED                 &
                , RecvRequest(iProc)                      &
                , IERROR)
        END IF
        IF (CNShapeMapping(iProc)%nSendShapeElems(1).GT.0) THEN
          CALL MPI_ISEND( CNShapeMapping(iProc)%SendShapeProcElemID   &
                , CNShapeMapping(iProc)%nSendShapeElems(1)*CNShapeMapping(iProc)%nSendShapeElems(2)   &
                , MPI_INTEGER                             &
                , iProc                                   &
                , 2012                                    &
                , MPI_COMM_LEADERS_SHARED                 &
                , SendRequest(iProc)                      &
                , IERROR)
        END IF
      END DO
      DO iProc = 0,nLeaderGroupProcs-1
        IF (iProc.EQ.myLeaderGroupRank) CYCLE
        IF (CNShapeMapping(iProc)%nRecvShapeElems(1).NE.0) THEN
          CALL MPI_WAIT(RecvRequest(iProc),MPIStatus,IERROR)
          IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
        END IF

        IF (CNShapeMapping(iProc)%nSendShapeElems(1).NE.0) THEN
          CALL MPI_WAIT(SendRequest(iProc),MPIStatus,IERROR)
          IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
        END IF
      END DO

      DO iProc = 0,nLeaderGroupProcs-1
        IF (iProc.EQ.myLeaderGroupRank) CYCLE
        DO iElem = 1, CNShapeMapping(iProc)%nRecvShapeElems(1)
          CNElemID = GetCNElemID(CNShapeMapping(iProc)%RecvShapeElemID(iElem))
          DO jProc = 1, CNShapeMapping(iProc)%nRecvShapeElems(2)
            ProcID = CNShapeMapping(iProc)%RecvShapeProcElemID(iElem, jProc)
            IF (ProcID.EQ.0) EXIT
            ShapeElemProcSend_Shared(CNElemID, ProcID) = .TRUE.
          END DO
        END DO
      END DO

      SDEALLOCATE(SendRequest)
      SDEALLOCATE(RecvRequest)
      SDEALLOCATE(CNShapeMapping)
    END IF ! nLeaderGroupProcs.GT.1

  ! .NOT. ComputeNodeRoot
  ELSE
    ALLOCATE(SendRequest(1))
    CALL MPI_ISEND( nSendShapeElems                          &
                  , 1                                        &
                  , MPI_INTEGER                              &
                  , 0                                        &
                  , 2001                                     &
                  , MPI_COMM_SHARED                          &
                  , SendRequest(1)                           &
                  , IERROR)

    CALL MPI_WAIT(SendRequest(1),MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)

    IF (nSendShapeElems.GE.1) THEN
      CALL MPI_ISEND( SendShapeElemID                        &
                    , nSendShapeElems                        &
                    , MPI_INTEGER                            &
                    , 0                                      &
                    , 2001                                   &
                    , MPI_COMM_SHARED                        &
                    , SendRequest(1)                         &
                    , IERROR)

      CALL MPI_WAIT(SendRequest(1),MPIStatus,IERROR)
      IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
    END IF

    DEALLOCATE(SendRequest)
  END IF

  CALL BARRIER_AND_SYNC(ShapeElemProcSend_Shared_Win,MPI_COMM_SHARED)

  !REcv part1
  ALLOCATE(RecvProcs(nProcessors), RecvProcsElems(nProcessors))
  RecvProcs = .FALSE.
  RecvProcsElems = 0
  DO iElem = 1, nElems
    CNElemID = GetCNElemID(iELem+offSetElem)
    DO iProc = 1, nProcessors
      IF (myRank.EQ.iProc-1) CYCLE
      IF (ShapeElemProcSend_Shared(CNElemID,iProc)) THEN
        RecvProcs(iProc) = .TRUE.
        RecvProcsElems(iProc) = RecvProcsElems(iProc) + 1
      END IF
    END DO
  END DO

  nShapeExchangeProcs = COUNT(RecvProcs)
  ALLOCATE(ShapeMapping(1:nShapeExchangeProcs), &
           RecvRequest (1:nShapeExchangeProcs), SendRequest(1:nShapeExchangeProcs))

  !REcv Part2
  exProc = 0
  DO iProc = 1, nProcessors
    IF (myRank.EQ.iProc-1) CYCLE
    IF (RecvProcs(iProc)) THEN
      exProc = exProc + 1
      ShapeMapping(exProc)%Rank = iProc -1
      CNRankToSendRank(ShapeMapping(exProc)%Rank) = exProc
      ShapeMapping(exProc)%nRecvShapeElems = RecvProcsElems(iProc)
      ShapeMapping(exProc)%nSendShapeElems = 0
      ALLOCATE(ShapeMapping(exProc)%RecvShapeElemID(1:RecvProcsElems(iProc)), &
            ShapeMapping(exProc)%RecvBuffer(4,0:PP_N,0:PP_N,0:PP_N,1:RecvProcsElems(iProc)))
      exElem = 0
      DO iElem = 1, nElems
        CNElemID = GetCNElemID(iELem+offSetElem)
        IF (ShapeElemProcSend_Shared(CNElemID,iProc)) THEN
          exElem = exElem + 1
          ShapeMapping(exProc)%RecvShapeElemID(exElem) = iElem+offSetElem
        END IF
      END DO
    END IF
  END DO

  DO iProc = 1,nShapeExchangeProcs
    CALL MPI_IRECV( ShapeMapping(iProc)%nSendShapeElems   &
                  , 1                                       &
                  , MPI_INTEGER                             &
                  , ShapeMapping(iProc)%Rank                                    &
                  , 2003                                    &
                  , MPI_COMM_WORLD                 &
                  , RecvRequest(iProc)                      &
                  , IERROR)
  END DO

  DO iProc = 1,nShapeExchangeProcs
    CALL MPI_ISEND( ShapeMapping(iProc)%nRecvShapeElems   &
                  , 1                                       &
                  , MPI_INTEGER                             &
                  , ShapeMapping(iProc)%Rank                                   &
                  , 2003                                    &
                  , MPI_COMM_WORLD                 &
                  , SendRequest(iProc)                     &
                  , IERROR)
  END DO

  DO iProc = 1,nShapeExchangeProcs
    CALL MPI_WAIT(RecvRequest(iProc),MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
    CALL MPI_WAIT(SendRequest(iProc),MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
    ALLOCATE(ShapeMapping(iProc)%SendShapeElemID(1:ShapeMapping(iProc)%nSendShapeElems), &
        ShapeMapping(iProc)%SendBuffer(4,0:PP_N,0:PP_N,0:PP_N,1:ShapeMapping(iProc)%nSendShapeElems))
  END DO

  DO iProc = 1, nShapeExchangeProcs
    CALL MPI_IRECV( ShapeMapping(iProc)%SendShapeElemID   &
              , ShapeMapping(iProc)%nSendShapeElems   &
              , MPI_INTEGER                             &
              , ShapeMapping(iProc)%Rank                &
              , 2003                                    &
              , MPI_COMM_WORLD                 &
              , RecvRequest(iProc)                      &
              , IERROR)
  END DO
  DO iProc = 1, nShapeExchangeProcs
    CALL MPI_ISEND( ShapeMapping(iProc)%RecvShapeElemID   &
                  , ShapeMapping(iProc)%nRecvShapeElems   &
                  , MPI_INTEGER                             &
                  , ShapeMapping(iProc)%Rank                &
                  , 2003                                    &
                  , MPI_COMM_WORLD                         &
                  , SendRequest(iProc)                      &
                  , IERROR)
  END DO

  DO iProc = 1,nShapeExchangeProcs
    CALL MPI_WAIT(RecvRequest(iProc),MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
    CALL MPI_WAIT(SendRequest(iProc),MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
    DO iELem = 1, ShapeMapping(iProc)%nSendShapeElems
      SendElemShapeID(GetCNElemID(ShapeMapping(iProc)%SendShapeElemID(iElem))) = iElem
    END DO
  END DO

  CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
  CALL UNLOCK_AND_FREE(ShapeElemProcSend_Shared_Win)
  CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

  ADEALLOCATE(ShapeElemProcSend_Shared)
END IF

SWRITE(UNIT_stdOut,'(A)') ' IDENTIFYING Particle Exchange Processors DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE IdentifyPartExchangeProcs


!===================================================================================================================================
! Deallocates arrays for halo exchange
!===================================================================================================================================
SUBROUTINE FinalizePartExchangeProcs()
! MODULES
USE MOD_Particle_MPI_Vars       ,ONLY: ExchangeProcToGlobalProc,GlobalProcToExchangeProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE(ExchangeProcToGlobalProc)
SDEALLOCATE(GlobalProcToExchangeProc)

END SUBROUTINE FinalizePartExchangeProcs


PPURE FUNCTION HaloBoxInProc(CartNodes,CartProc,halo_eps,nPeriodicVectors,PeriodicVectors)
!===================================================================================================================================
! Check if bounding box is on proc by comparing against the other bounding box extended by halo_eps
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)             :: CartNodes(6)
REAL,INTENT(IN)             :: CartProc( 6)
REAL,INTENT(IN)             :: halo_eps
INTEGER,INTENT(IN)          :: nPeriodicVectors
REAL,INTENT(IN),OPTIONAL    :: PeriodicVectors(3,nPeriodicVectors)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL                     :: HaloBoxInProc
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,DIMENSION(2),PARAMETER :: DirPeriodicVector = [-1,1]
INTEGER                        :: iPeriodicVector,jPeriodicVector,iPeriodicDir,jPeriodicDir,kPeriodicDir
REAL,DIMENSION(1:6)            :: xCordsPeri
!===================================================================================================================================

HaloBoxInProc = .FALSE.

! Check whether the bounding boxes intersect
IF (   ((CartNodes(1).LE.CartProc(2)+halo_eps).AND.(CartNodes(2).GE.CartProc(1)-halo_eps))  &
  .AND.((CartNodes(3).LE.CartProc(4)+halo_eps).AND.(CartNodes(4).GE.CartProc(3)-halo_eps))  &
  .AND.((CartNodes(5).LE.CartProc(6)+halo_eps).AND.(CartNodes(6).GE.CartProc(5)-halo_eps))) HaloBoxInProc = .TRUE.

! Also check periodic directions. Only MPI sides of the local proc are taken into account, so do not perform additional case distinction
SELECT CASE(nPeriodicVectors)
  ! One periodic vector
  CASE(1)
    DO iPeriodicDir = 1,2
      ! check if element is within halo_eps of periodically displaced element
      xCordsPeri(1:2) = CartNodes(1:2) + PeriodicVectors(1,1) * DirPeriodicVector(iPeriodicDir)
      xCordsPeri(3:4) = CartNodes(3:4) + PeriodicVectors(2,1) * DirPeriodicVector(iPeriodicDir)
      xCordsPeri(5:6) = CartNodes(5:6) + PeriodicVectors(3,1) * DirPeriodicVector(iPeriodicDir)

      ! Check whether the bounding boxes intersect
      IF (   ((xCordsPeri(1).LE.CartProc(2)+halo_eps).AND.(xCordsPeri(2).GE.CartProc(1)-halo_eps))  &
        .AND.((xCordsPeri(3).LE.CartProc(4)+halo_eps).AND.(xCordsPeri(4).GE.CartProc(3)-halo_eps))  &
        .AND.((xCordsPeri(5).LE.CartProc(6)+halo_eps).AND.(xCordsPeri(6).GE.CartProc(5)-halo_eps))) THEN
        HaloBoxInProc = .TRUE.
        RETURN
      END IF
    END DO

  ! Two periodic vectors. Also check linear combination, see particle_bgm.f90
  CASE(2)
    DO iPeriodicVector = 1,2
      DO iPeriodicDir = 1,2
        ! check if element is within halo_eps of periodically displaced element
        xCordsPeri(1:2) = CartNodes(1:2) + PeriodicVectors(1,iPeriodicVector) * DirPeriodicVector(iPeriodicDir)
        xCordsPeri(3:4) = CartNodes(3:4) + PeriodicVectors(2,iPeriodicVector) * DirPeriodicVector(iPeriodicDir)
        xCordsPeri(5:6) = CartNodes(5:6) + PeriodicVectors(3,iPeriodicVector) * DirPeriodicVector(iPeriodicDir)

        ! Check whether the bounding boxes intersect
        IF (   ((xCordsPeri(1).LE.CartProc(2)+halo_eps).AND.(xCordsPeri(2).GE.CartProc(1)-halo_eps))  &
          .AND.((xCordsPeri(3).LE.CartProc(4)+halo_eps).AND.(xCordsPeri(4).GE.CartProc(3)-halo_eps))  &
          .AND.((xCordsPeri(5).LE.CartProc(6)+halo_eps).AND.(xCordsPeri(6).GE.CartProc(5)-halo_eps))) THEN
          HaloBoxInProc = .TRUE.
          RETURN
        END IF
      END DO
    END DO

    DO iPeriodicDir = 1,2
      DO jPeriodicDir = 1,2
        ! check if element is within halo_eps of periodically displaced element
        xCordsPeri(1:2) = CartNodes(1:2) + PeriodicVectors(1,1) * DirPeriodicVector(iPeriodicDir) &
                                         + PeriodicVectors(1,2) * DirPeriodicVector(jPeriodicDir)
        xCordsPeri(3:4) = CartNodes(3:4) + PeriodicVectors(2,1) * DirPeriodicVector(iPeriodicDir) &
                                         + PeriodicVectors(2,2) * DirPeriodicVector(jPeriodicDir)
        xCordsPeri(5:6) = CartNodes(5:6) + PeriodicVectors(3,1) * DirPeriodicVector(iPeriodicDir) &
                                         + PeriodicVectors(3,2) * DirPeriodicVector(jPeriodicDir)

        ! Check whether the bounding boxes intersect
        IF (   ((xCordsPeri(1).LE.CartProc(2)+halo_eps).AND.(xCordsPeri(2).GE.CartProc(1)-halo_eps))  &
          .AND.((xCordsPeri(3).LE.CartProc(4)+halo_eps).AND.(xCordsPeri(4).GE.CartProc(3)-halo_eps))  &
          .AND.((xCordsPeri(5).LE.CartProc(6)+halo_eps).AND.(xCordsPeri(6).GE.CartProc(5)-halo_eps))) THEN
          HaloBoxInProc = .TRUE.
          RETURN
        END IF
      END DO
    END DO

  ! Three periodic vectors. Also check linear combination, see particle_bgm.f90
  CASE(3)
    ! check the three periodic vectors. Begin with checking the first periodic vector, followed by the combination of
    ! the first periodic vector with the others. Then check the other combinations, i.e. 1, 1+2, 1+3, 2, 2+3, 3, 1+2+3
    DO iPeriodicVector = 1,3
      DO iPeriodicDir = 1,2
        ! check if element is within halo_eps of periodically displaced element
        xCordsPeri(1:2) = CartNodes(1:2) + PeriodicVectors(1,iPeriodicVector) * DirPeriodicVector(iPeriodicDir)
        xCordsPeri(3:4) = CartNodes(3:4) + PeriodicVectors(2,iPeriodicVector) * DirPeriodicVector(iPeriodicDir)
        xCordsPeri(5:6) = CartNodes(5:6) + PeriodicVectors(3,iPeriodicVector) * DirPeriodicVector(iPeriodicDir)

        ! Check whether the bounding boxes intersect
        IF (   ((xCordsPeri(1).LE.CartProc(2)+halo_eps).AND.(xCordsPeri(2).GE.CartProc(1)-halo_eps))  &
          .AND.((xCordsPeri(3).LE.CartProc(4)+halo_eps).AND.(xCordsPeri(4).GE.CartProc(3)-halo_eps))  &
          .AND.((xCordsPeri(5).LE.CartProc(6)+halo_eps).AND.(xCordsPeri(6).GE.CartProc(5)-halo_eps))) THEN
          HaloBoxInProc = .TRUE.
          RETURN
        END IF

        DO jPeriodicVector = 1,3
          DO jPeriodicDir = 1,2
            IF (iPeriodicVector.GE.jPeriodicVector) CYCLE
            ! check if element is within halo_eps of periodically displaced element
            xCordsPeri(1:2) = CartNodes(1:2) + PeriodicVectors(1,iPeriodicVector) * DirPeriodicVector(iPeriodicDir) &
                                             + PeriodicVectors(1,jPeriodicVector) * DirPeriodicVector(jPeriodicDir)
            xCordsPeri(3:4) = CartNodes(3:4) + PeriodicVectors(2,iPeriodicVector) * DirPeriodicVector(iPeriodicDir) &
                                             + PeriodicVectors(2,jPeriodicVector) * DirPeriodicVector(jPeriodicDir)
            xCordsPeri(5:6) = CartNodes(5:6) + PeriodicVectors(3,iPeriodicVector) * DirPeriodicVector(iPeriodicDir) &
                                             + PeriodicVectors(3,jPeriodicVector) * DirPeriodicVector(jPeriodicDir)

            ! Check whether the bounding boxes intersect
            IF (   ((xCordsPeri(1).LE.CartProc(2)+halo_eps).AND.(xCordsPeri(2).GE.CartProc(1)-halo_eps))  &
              .AND.((xCordsPeri(3).LE.CartProc(4)+halo_eps).AND.(xCordsPeri(4).GE.CartProc(3)-halo_eps))  &
              .AND.((xCordsPeri(5).LE.CartProc(6)+halo_eps).AND.(xCordsPeri(6).GE.CartProc(5)-halo_eps))) THEN
              HaloBoxInProc = .TRUE.
              RETURN
            END IF
          END DO
        END DO
      END DO
    END DO

    ! check if element is within halo_eps of periodically displaced element
    DO iPeriodicDir = 1,2
      DO jPeriodicDir = 1,2
        DO kPeriodicDir = 1,2
          ! check if element is within halo_eps of periodically displaced element
          xCordsPeri(1:2) = CartNodes(1:2) + PeriodicVectors(1,1) * DirPeriodicVector(iPeriodicDir) &
                                           + PeriodicVectors(1,2) * DirPeriodicVector(jPeriodicDir) &
                                           + PeriodicVectors(1,3) * DirPeriodicVector(kPeriodicDir)
          xCordsPeri(3:4) = CartNodes(3:4) + PeriodicVectors(2,1) * DirPeriodicVector(iPeriodicDir) &
                                           + PeriodicVectors(2,2) * DirPeriodicVector(jPeriodicDir) &
                                           + PeriodicVectors(2,3) * DirPeriodicVector(kPeriodicDir)
          xCordsPeri(5:6) = CartNodes(5:6) + PeriodicVectors(3,1) * DirPeriodicVector(iPeriodicDir) &
                                           + PeriodicVectors(3,2) * DirPeriodicVector(jPeriodicDir) &
                                           + PeriodicVectors(3,3) * DirPeriodicVector(kPeriodicDir)

          ! Check whether the bounding boxes intersect
          IF (   ((xCordsPeri(1).LE.CartProc(2)+halo_eps).AND.(xCordsPeri(2).GE.CartProc(1)-halo_eps))  &
            .AND.((xCordsPeri(3).LE.CartProc(4)+halo_eps).AND.(xCordsPeri(4).GE.CartProc(3)-halo_eps))  &
            .AND.((xCordsPeri(5).LE.CartProc(6)+halo_eps).AND.(xCordsPeri(6).GE.CartProc(5)-halo_eps))) THEN
            HaloBoxInProc = .TRUE.
            RETURN
          END IF
        END DO
      END DO
    END DO

END SELECT

END FUNCTION HaloBoxInProc
#endif /*USE_MPI*/

END MODULE MOD_Particle_MPI_Halo
