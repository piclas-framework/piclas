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
USE MOD_Globals_Vars            ,ONLY: c
USE MOD_Preproc
USE MOD_Mesh_Vars               ,ONLY: nElems,offsetElem
USE MOD_MPI_Shared_Vars
USE MOD_Mesh_Tools              ,ONLY: GetGlobalElemID
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_MPI_Vars       ,ONLY: halo_eps,halo_eps_velo
USE MOD_Particle_MPI_Vars       ,ONLY: nExchangeProcessors,ExchangeProcToGlobalProc,GlobalProcToExchangeProc
USE MOD_Particle_Vars           ,ONLY: ManualTimeStep
USE MOD_PICDepo_Vars            ,ONLY: DepositionType
USE MOD_PICDepo_Vars            ,ONLY: nSendShapeElems,SendShapeElemID
USE MOD_PICDepo_Vars            ,ONLY: ShapeMapping,CNShapeMapping
#if ! (USE_HDG)
USE MOD_CalcTimeStep            ,ONLY: CalcTimeStep
#endif
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
USE MOD_TimeDisc_Vars           ,ONLY: nRKStages,RK_c
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iPeriodicVector,jPeriodicVector,iPeriodicDir
INTEGER,DIMENSION(2)           :: DirPeriodicVector = [-1,1]
INTEGER                        :: iElem,ElemID,firstElem,lastElem,NbElemID
INTEGER                        :: iSide,SideID,iLocSide
!INTEGER                        :: firstSide,lastSide
INTEGER                        :: iMortar,nMortarElems,NbSideID
INTEGER                        :: iProc,HaloProc
!INTEGER                        :: GlobalProcID
INTEGER                        :: nExchangeSides
!INTEGER                        :: nExchangeProcs
INTEGER,ALLOCATABLE            :: ExchangeSides(:)
REAL,ALLOCATABLE               :: BoundsOfElemCenter(:),MPISideBoundsOfElemCenter(:,:)
INTEGER                        :: ExchangeProcLeader
! halo_eps reconstruction
REAL                           :: MPI_halo_eps,MPI_halo_eps_velo,MPI_halo_diag,vec(1:3),deltaT
LOGICAL                        :: fullMesh
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
INTEGER                        :: iStage
#endif
! shape function
INTEGER                        :: GlobalElemID,GlobalElemRank,GlobalLeaderRank
LOGICAL,ALLOCATABLE            :: FlagShapeElem(:)
! Non-symmetric particle exchange
INTEGER,ALLOCATABLE            :: SendRequest(:),RecvRequest(:)
LOGICAL,ALLOCATABLE            :: GlobalProcToRecvProc(:)
LOGICAL                        :: CommFlag
INTEGER                        :: nNonSymmetricExchangeProcs,nNonSymmetricExchangeProcsGlob
!=================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' IDENTIFYING Particle Exchange Processors  ...'

! Allocate arrays
ALLOCATE(GlobalProcToExchangeProc(EXCHANGE_PROC_SIZE,0:nProcessors_Global-1))
GlobalProcToExchangeProc(:,:) = -1

! Identify all procs on same node
nExchangeProcessors = 0
!DO iProc = ComputeNodeRootRank,ComputeNodeRootRank+nComputeNodeProcessors-1
!  ! Build mapping global to compute-node
!  GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc) = 1
!  GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc) = nExchangeProcessors
!  nExchangeProcessors = nExchangeProcessors + 1
!END DO

! Communicate non-symmetric part exchange partners to catch non-symmetric proc identification due to inverse distance calculation
ALLOCATE(GlobalProcToRecvProc(0:nProcessors_Global-1), &
         SendRequest         (0:nProcessors_Global-1), &
         RecvRequest         (0:nProcessors_Global-1))

GlobalProcToRecvProc = .FALSE.

! Identify all procs with elements in range. This includes checking the procs on the compute-node as they might lie far apart
IF (nProcessors.EQ.1) THEN
  SWRITE(UNIT_stdOut,'(A)') ' | Running on one processor. Particle exchange communication disabled.'
  SWRITE(UNIT_stdOut,'(A)') ' IDENTIFYING Particle Exchange Processors DONE!'
  SWRITE(UNIT_StdOut,'(132("-"))')
  RETURN
END IF

! Open receive buffer for non-symmetric exchange identification
DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE

  CALL MPI_IRECV( GlobalProcToRecvProc(iProc)  &
                , 1                            &
                , MPI_INTEGER                  &
                , iProc                        &
                , 1999                         &
                , MPI_COMM_WORLD               &
                , RecvRequest(iProc)           &
                , IERROR)
END DO

!> Count all MPI sides on current proc.
firstElem = offsetElem+1
lastElem  = offsetElem+nElems

IF(TRIM(DepositionType(1:MIN(14,LEN(TRIM(ADJUSTL(DepositionType)))))).EQ.'shape_function')THEN
  ALLOCATE(FlagShapeElem(1:nComputeNodeElems))
  FlagShapeElem = .FALSE.
END IF

! This approach does not work, we get only the MPI sides pointing into the compute-node halo region
!DO iSide = firstSide,lastSide
!  IF (SideInfo_Shared(SIDE_NBELEMTYPE,iSide).EQ.2) THEN
!    nExchangeSides = nExchangeSides + 1
!  END IF
!END DO

!>>> For all element, loop over the six sides and check if the neighbor element is on the current proc
!>>> Special care for big mortar sides, here the SIDE_ELEMID must be used
nExchangeSides = 0

DO iElem = firstElem,lastElem
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

        NbElemID = SideInfo_Shared(SIDE_ELEMID,NbSideID)
        ! If small mortar element not defined, skip it for now, likely not inside the halo region
        IF (NbElemID.LT.1) CYCLE

        ! If any of the small mortar sides is not on the local proc, the side is a MPI side
        IF (NbElemID.LT.firstElem .OR. NbElemID.GT.lastElem) THEN
          nExchangeSides = nExchangeSides + 1
          EXIT
        END IF
      END DO

    ! regular side or small mortar side
    ELSE
      IF (NbElemID.LT.firstElem .OR. NbElemID.GT.lastElem) THEN
        nExchangeSides = nExchangeSides + 1
      END IF
    END IF
  END DO
END DO

IF (nComputeNodeProcessors.GT.1.AND.nExchangeSides.EQ.0) &
  CALL ABORT(__STAMP__,'Found no side connectivity between processor domains')

!> Build mapping for all MPI sides on current proc
ALLOCATE(ExchangeSides(1:nExchangeSides))

nExchangeSides = 0

! This approach does not work, we get only the MPI sides pointing into the compute-node halo region
!DO iSide = firstSide,lastSide
!  IF (SideInfo_Shared(SIDE_NBELEMTYPE,iSide).EQ.2) THEN
!    nExchangeSides = nExchangeSides + 1
!    ExchangeSides(nExchangeSides) = SideInfo_Shared(SIDE_ID,iSide)
!  END IF
!END DO

DO iElem = firstElem,lastElem
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

        NbElemID = SideInfo_Shared(SIDE_ELEMID,NbSideID)
        ! If small mortar element not defined, skip it for now, likely not inside the halo region
        IF (NbElemID.LT.1) CYCLE

        ! If any of the small mortar sides is not on the local proc, the side is a MPI side
        IF (NbElemID.LT.firstElem .OR. NbElemID.GT.lastElem) THEN
          nExchangeSides = nExchangeSides + 1
          ExchangeSides(nExchangeSides) = SideID
          EXIT
        END IF
      END DO

    ! regular side or small mortar side
    ELSE
      IF (NbElemID.LT.firstElem .OR. NbElemID.GT.lastElem) THEN
        nExchangeSides = nExchangeSides + 1
        ExchangeSides(nExchangeSides) = SideID
      END IF
    END IF
  END DO
END DO

!> Build metrics for all MPI sides on current proc
ALLOCATE(BoundsOfElemCenter(1:4))
ALLOCATE(MPISideBoundsOfElemCenter(1:4,1:nExchangeSides))

DO iSide = 1, nExchangeSides
  SideID = ExchangeSides(iSide)
  ElemID = SideInfo_Shared(SIDE_ELEMID,SideID)
  MPISideBoundsOfElemCenter(1:3,iSide) = (/ SUM(BoundsOfElem_Shared(1:2,1,ElemID)), &
                                            SUM(BoundsOfElem_Shared(1:2,2,ElemID)), &
                                            SUM(BoundsOfElem_Shared(1:2,3,ElemID)) /) / 2.
  MPISideBoundsOfElemCenter(4,iSide) = VECNORM ((/BoundsOfElem_Shared(2,1,ElemID)-BoundsOfElem_Shared(1,1,ElemID), &
                                                  BoundsOfElem_Shared(2,2,ElemID)-BoundsOfElem_Shared(1,2,ElemID), &
                                                  BoundsOfElem_Shared(2,3,ElemID)-BoundsOfElem_Shared(1,3,ElemID) /) / 2.)
END DO

!> Check all elements in the CN halo region against local MPI sides. Check is identical to particle_bgm.f90
!>>> Check the bounding box of each element in compute-nodes' halo domain against the bounding boxes of the
!>>> of the elements of the MPI-surface (local proc MPI sides)

! if running on one node, halo_eps is meaningless. Get a representative MPI_halo_eps for BC side identification
fullMesh = .FALSE.
IF (halo_eps.EQ.0) THEN
  ! reconstruct halo_eps_velo
  IF (halo_eps_velo.EQ.0) THEN
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
  MPI_halo_eps = MPI_halo_eps*MPI_halo_eps_velo*deltaT
#else
  MPI_halo_eps = MPI_halo_eps_velo*deltaT
#endif

  vec(1)   = GEO%xmaxglob-GEO%xminglob
  vec(2)   = GEO%ymaxglob-GEO%yminglob
  vec(3)   = GEO%zmaxglob-GEO%zminglob
  MPI_halo_diag = VECNORM(vec)

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
END IF

! Use a named loop so the entire element can be cycled
ElemLoop:  DO iElem = 1,nComputeNodeTotalElems
  ElemID   = GetGlobalElemID(iElem)
  HaloProc = ElemInfo_Shared(ELEM_RANK,ElemID)

  IF (HaloProc.EQ.myRank) CYCLE

!#if CODE_ANALYZE
!    ! Sanity checks. Elems in halo region must have ELEM_HALOFLAG=2 and the proc must not be flagged yet
!    IF (ElemInfo_Shared(ELEM_HALOFLAG,ElemID).NE.2) THEN
!      IPWRITE(UNIT_stdOut,*) 'Element ID:',ElemID,'Halo Flag: ',ElemInfo_Shared(ELEM_HALOFLAG,ElemID)
!      CALL ABORT(__STAMP__,  'Element found in range of halo elements while not flagged as such!')
!    END IF
!
!    IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc).EQ.1) THEN
!      IPWRITE(UNIT_stdOut,*) 'Element ID:',ElemID,'Halo Proc: ',HaloProc
!      CALL ABORT(__STAMP__, 'Proc claimed to have elements both on compute node and in halo region!')
!    END IF
!#endif

  ! Skip if the proc is already flagged, only if the exact elements are not required (.NOT.shape_function)
  IF(.NOT.TRIM(DepositionType(1:MIN(14,LEN(TRIM(ADJUSTL(DepositionType)))))).EQ.'shape_function')THEN
    IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc).EQ.2) CYCLE
  END IF

  BoundsOfElemCenter(1:3) = (/SUM(BoundsOfElem_Shared(1:2,1,ElemID)), &
                              SUM(BoundsOfElem_Shared(1:2,2,ElemID)), &
                              SUM(BoundsOfElem_Shared(1:2,3,ElemID)) /) / 2.
  BoundsOfElemCenter(4)   = VECNORM ((/ BoundsOfElem_Shared(2,1,ElemID)-BoundsOfElem_Shared(1,1,ElemID), &
                                        BoundsOfElem_Shared(2,2,ElemID)-BoundsOfElem_Shared(1,2,ElemID), &
                                        BoundsOfElem_Shared(2,3,ElemID)-BoundsOfElem_Shared(1,3,ElemID) /) / 2.)
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
              .LE. MPI_halo_eps+BoundsOfElemCenter(4)+MPISideBoundsOfElemCenter(4,iSide) ) THEN
                ! flag the proc as exchange proc (in halo region)
                GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
                GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
                nExchangeProcessors = nExchangeProcessors + 1
                IF(TRIM(DepositionType(1:MIN(14,LEN(TRIM(ADJUSTL(DepositionType)))))).EQ.'shape_function')THEN
                  FlagShapeElem(iElem) = .TRUE.
                END IF
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
                        .LE. MPI_halo_eps+BoundsOfElemCenter(4)+MPISideBoundsOfElemCenter(4,iSide) ) THEN
                ! flag the proc as exchange proc (in halo region)
                GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
                GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
                nExchangeProcessors = nExchangeProcessors + 1
                IF(TRIM(DepositionType(1:MIN(14,LEN(TRIM(ADJUSTL(DepositionType)))))).EQ.'shape_function')THEN
                  FlagShapeElem(iElem) = .TRUE.
                END IF
                CYCLE ElemLoop
              END IF

              DO jPeriodicVector = 1,2

                ! check if element is within halo_eps of periodically displaced element
                IF (VECNORM( BoundsOfElemCenter(1:3)                                                    &
                           + GEO%PeriodicVectors(1:3,iPeriodicVector) * DirPeriodicVector(iPeriodicDir) &
                           + GEO%PeriodicVectors(1:3,jPeriodicVector) * DirPeriodicVector(iPeriodicDir) &
                           - MPISideBoundsOfElemCenter(1:3,iSide))                                      &
                        .LE. MPI_halo_eps+BoundsOfElemCenter(4)+MPISideBoundsOfElemCenter(4,iSide) ) THEN
                  ! flag the proc as exchange proc (in halo region)
                  GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
                  GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
                  nExchangeProcessors = nExchangeProcessors + 1
                  IF(TRIM(DepositionType(1:MIN(14,LEN(TRIM(ADJUSTL(DepositionType)))))).EQ.'shape_function')THEN
                    FlagShapeElem(iElem) = .TRUE.
                  END IF
                  CYCLE ElemLoop
                END IF
              END DO
            END DO
          END DO

        ! Two periodic vectors. Also check linear combination, see particle_bgm.f90
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
                        .LE. MPI_halo_eps+BoundsOfElemCenter(4)+MPISideBoundsOfElemCenter(4,iSide) ) THEN
                ! flag the proc as exchange proc (in halo region)
                GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
                GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
                nExchangeProcessors = nExchangeProcessors + 1
                IF(TRIM(DepositionType(1:MIN(14,LEN(TRIM(ADJUSTL(DepositionType)))))).EQ.'shape_function')THEN
                  FlagShapeElem(iElem) = .TRUE.
                END IF
                CYCLE ElemLoop
              END IF

              DO jPeriodicVector = 1,3
                IF (iPeriodicVector.GE.jPeriodicVector) CYCLE

                ! check if element is within halo_eps of periodically displaced element
                IF (VECNORM( BoundsOfElemCenter(1:3)                                                    &
                           + GEO%PeriodicVectors(1:3,iPeriodicVector) * DirPeriodicVector(iPeriodicDir) &
                           + GEO%PeriodicVectors(1:3,jPeriodicVector) * DirPeriodicVector(iPeriodicDir) &
                           - MPISideBoundsOfElemCenter(1:3,iSide))                                      &
                        .LE. MPI_halo_eps+BoundsOfElemCenter(4)+MPISideBoundsOfElemCenter(4,iSide) ) THEN
                  ! flag the proc as exchange proc (in halo region)
                  GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
                  GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
                  nExchangeProcessors = nExchangeProcessors + 1
                  IF(TRIM(DepositionType(1:MIN(14,LEN(TRIM(ADJUSTL(DepositionType)))))).EQ.'shape_function')THEN
                    FlagShapeElem(iElem) = .TRUE.
                  END IF
                  CYCLE ElemLoop
                END IF

              END DO
            END DO
          END DO

          ! check if element is within halo_eps of periodically displaced element
          DO iPeriodicDir = 1,2
            IF (VECNORM( BoundsOfElemCenter(1:3)                                                        &
                       + GEO%PeriodicVectors(1:3,1) * DirPeriodicVector(iPeriodicDir)                   &
                       + GEO%PeriodicVectors(1:3,2) * DirPeriodicVector(iPeriodicDir)                   &
                       + GEO%PeriodicVectors(1:3,3) * DirPeriodicVector(iPeriodicDir)                   &
                       - MPISideBoundsOfElemCenter(1:3,iSide))                                          &
                    .LE. MPI_halo_eps+BoundsOfElemCenter(4)+MPISideBoundsOfElemCenter(4,iSide) ) THEN
              ! flag the proc as exchange proc (in halo region)
              GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
              GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
              nExchangeProcessors = nExchangeProcessors + 1
              IF(TRIM(DepositionType(1:MIN(14,LEN(TRIM(ADJUSTL(DepositionType)))))).EQ.'shape_function')THEN
                FlagShapeElem(iElem) = .TRUE.
              END IF
              CYCLE ElemLoop
            END IF
          END DO

        ! No periodic vectors, element out of range
        CASE(0)
          ! Do nothing

        CASE DEFAULT
          CALL ABORT(__STAMP__,'Invalid number of periodic vectors in particle_mpi_halo.f90')

      END SELECT

    ! Element is in range of not-periodically displaced MPI side
    ELSE
      GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,HaloProc) = 2
      GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,HaloProc) = nExchangeProcessors
      nExchangeProcessors = nExchangeProcessors + 1
      IF(TRIM(DepositionType(1:MIN(14,LEN(TRIM(ADJUSTL(DepositionType)))))).EQ.'shape_function')THEN
        FlagShapeElem(iElem) = .TRUE.
      END IF
      CYCLE ElemLoop
    END IF
  END DO ! iSide = 1, nExchangeSides
END DO ElemLoop

! Notify every proc which was identified by the local proc
DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE

  ! CommFlag holds the information if the local proc wants to communicate with iProc
  CommFlag = MERGE(.TRUE.,.FALSE.,GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc).NE.-1)
  CALL MPI_ISEND( CommFlag                     &
                , 1                            &
                , MPI_INTEGER                  &
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
DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE

  ! Ignore procs that are already flagged or not requesting communication
  IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc) .NE.-1) CYCLE
  IF (.NOT.GlobalProcToRecvProc(iProc)) CYCLE

  ! Found a previously missing proc
  nNonSymmetricExchangeProcs = nNonSymmetricExchangeProcs + 1
  nExchangeProcessors        = nExchangeProcessors + 1

  GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc) = 2
  GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc) = nExchangeProcessors
END DO

DEALLOCATE(GlobalProcToRecvProc,RecvRequest,SendRequest)
CALL MPI_REDUCE(nNonSymmetricExchangeProcs,nNonSymmetricExchangeProcsGlob,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
IF (nNonSymmetricExchangeProcsGlob.GT.1) THEN
  SWRITE(Unit_StdOut,'(A,I0,A)') ' | Found ',nNonSymmetricExchangeProcsGlob, &
                                 ' previously missing non-symmetric particle exchange procs'
END IF

! Build reverse mapping
!-- EXCHANGE_PROC_TYPE information is currently unused and either -1 (no communication) or 2 (communication). Can be used to
!-- implement check if exchange partner is on the same compute node, so build it here
ALLOCATE(ExchangeProcToGlobalProc(2,0:nExchangeProcessors-1))

! Loop through all procs and build reverse mapping
!DO iProc = 0,nComputeNodeProcessors-1
!  ExchangeProcToGlobalProc(EXCHANGE_PROC_TYPE,iProc) = 1
!  ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,iProc) = iProc + ComputeNodeRootRank
!END DO
!nExchangeProcessors = nComputeNodeProcessors

nExchangeProcessors = 0
DO iProc = 0,nProcessors_Global-1
  IF (GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc).NE.-1) THEN
    ! Find it the other proc is on the same compute node
    ExchangeProcLeader = INT(GlobalProcToExchangeProc(EXCHANGE_PROC_RANK,iProc)/nComputeNodeProcessors)
    IF (ExchangeProcLeader.EQ.myLeaderGroupRank) THEN
      GlobalProcToExchangeProc(EXCHANGE_PROC_TYPE,iProc) = 1
      ExchangeProcToGlobalProc(EXCHANGE_PROC_TYPE,nExchangeProcessors) = 1
      ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,nExchangeProcessors) = iProc
    ELSE
      ExchangeProcToGlobalProc(EXCHANGE_PROC_TYPE,nExchangeProcessors) = 2
      ExchangeProcToGlobalProc(EXCHANGE_PROC_RANK,nExchangeProcessors) = iProc
    END IF
    nExchangeProcessors = nExchangeProcessors +1
  END IF
END DO

! Build shapeFunction mapping
IF(TRIM(DepositionType(1:MIN(14,LEN(TRIM(ADJUSTL(DepositionType)))))).EQ.'shape_function')THEN
  nSendShapeElems = 0
  DO iELem = 1,nComputeNodeTotalElems
    IF (FlagShapeElem(iElem)) nSendShapeElems = nSendShapeElems + 1
  END DO
!  nSendShapeElems = SUM(FlagShapeElem)

  ALLOCATE(SendShapeElemID(1:nSendShapeElems))

  SendShapeElemID = -1

  nSendShapeElems = 0
  DO iELem = 1,nComputeNodeTotalElems
    IF (FlagShapeElem(iElem)) THEN
      nSendShapeElems = nSendShapeElems + 1
      SendShapeElemID(nSendShapeElems) = iElem
    END IF
  END DO

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

    ! Second stage of communication, identify and send inter-compute-node information
    IF (nLeaderGroupProcs.EQ.1) THEN
      DEALLOCATE(RecvRequest)
      ALLOCATE(CNShapeMapping(0:nLeaderGroupProcs-1), &
               SendRequest   (0:nLeaderGroupProcs-1), &
               RecvRequest   (0:nLeaderGroupProcs-1))

      DO iProc = 0,nLeaderGroupProcs-1
        IF (iProc.EQ.myLeaderGroupRank) CYCLE

        CALL MPI_IRECV( CNShapeMapping(iProc)%nRecvShapeElems   &
                      , 1                                       &
                      , MPI_INTEGER                             &
                      , iProc                                   &
                      , 2002                                    &
                      , MPI_COMM_LEADERS_SHARED                 &
                      , RecvRequest(iProc)                      &
                      , IERROR)
      END DO

      CNShapeMapping%nSendShapeElems = 0
      ! Count number of elems per CN
      DO iElem = nComputeNodeElems+1,nComputeNodeTotalElems
        GlobalElemID = GetGlobalElemID(iElem)
        IF (ElemInfo_Shared(ELEM_HALOFLAG,GlobalElemID).EQ.2) THEN
          GlobalElemRank = ElemInfo_Shared(ELEM_RANK,GlobalElemID)
          GlobalLeaderRank = INT(GlobalElemRank/nComputeNodeProcessors)
          CNShapeMapping(GlobalLeaderRank)%nSendShapeElems = CNShapeMapping(GlobalLeaderRank)%nSendShapeElems + 1
        END IF
      END DO

      DO iProc = 0,nLeaderGroupProcs-1
        IF (iProc.EQ.myLeaderGroupRank) CYCLE

        CALL MPI_ISEND( CNShapeMapping(iProc)%nSendShapeElems   &
                      , 1                                       &
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

        IF (CNShapeMapping(iProc)%nRecvShapeElems.EQ.0) CYCLE

        ALLOCATE(CNShapeMapping(iProc)%RecvShapeElemID(CNShapeMapping(iProc)%nRecvShapeElems))

        CALL MPI_IRECV( CNShapeMapping(iProc)%RecvShapeElemID   &
                      , CNShapeMapping(iProc)%nRecvShapeElems   &
                      , MPI_INTEGER                             &
                      , iProc                                   &
                      , 2002                                    &
                      , MPI_COMM_LEADERS_SHARED                 &
                      , RecvRequest(iProc)                      &
                      , IERROR)
      END DO

      DO iProc = 0,nLeaderGroupProcs-1
        IF (iProc.EQ.myLeaderGroupRank) CYCLE
        IF (CNShapeMapping(iProc)%nSendShapeElems.EQ.0) CYCLE
        ALLOCATE(CNShapeMapping(iProc)%SendShapeElemID(CNShapeMapping(iProc)%nSendShapeElems))
        CNShapeMapping(iProc)%nSendShapeElems = 0
      END DO


      DO iElem = nComputeNodeElems+1,nComputeNodeTotalElems
        GlobalElemID = GetGlobalElemID(iElem)
        IF (ElemInfo_Shared(ELEM_HALOFLAG,GlobalElemID).EQ.2) THEN
          GlobalElemRank = ElemInfo_Shared(ELEM_RANK,GlobalElemID)
          GlobalLeaderRank = INT(GlobalElemRank/nComputeNodeProcessors)

          CNShapeMapping(GlobalLeaderRank)%nSendShapeElems = CNShapeMapping(GlobalLeaderRank)%nSendShapeElems + 1
          CNShapeMapping(GlobalLeaderRank)%SendShapeElemID(CNShapeMapping(GlobalLeaderRank)%nSendShapeElems) = GlobalElemID
        END IF
      END DO

      DO iProc = 0,nLeaderGroupProcs-1
        IF (iProc.EQ.myLeaderGroupRank) CYCLE

        IF (CNShapeMapping(iProc)%nSendShapeElems.EQ.0) CYCLE

        ALLOCATE(CNShapeMapping(iProc)%SendShapeElemID(CNShapeMapping(iProc)%nSendShapeElems))

        CALL MPI_ISEND( CNShapeMapping(iProc)%SendShapeElemID   &
                      , CNShapeMapping(iProc)%nSendShapeElems   &
                      , MPI_INTEGER                             &
                      , iProc                                   &
                      , 2002                                    &
                      , MPI_COMM_LEADERS_SHARED                 &
                      , SendRequest(iProc)                      &
                      , IERROR)
      END DO

      DO iProc = 0,nLeaderGroupProcs-1
        IF (iProc.EQ.myLeaderGroupRank) CYCLE

        IF (CNShapeMapping(iProc)%nRecvShapeElems.NE.0) THEN
          CALL MPI_WAIT(RecvRequest(iProc),MPIStatus,IERROR)
          IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
        END IF

        IF (CNShapeMapping(iProc)%nSendShapeElems.NE.0) THEN
          CALL MPI_WAIT(SendRequest(iProc),MPIStatus,IERROR)
          IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
        END IF
      END DO
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

    IF (nSendShapeElems.GT.1) THEN
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
#endif /*USE_MPI*/

END MODULE MOD_Particle_MPI_Halo
