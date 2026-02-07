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

!===================================================================================================================================
!> Module for the HDG method
!===================================================================================================================================
MODULE MOD_HDG_Init
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
#if USE_HDG
PUBLIC :: InitFPC
PUBLIC :: InitEPC
#if defined(PARTICLES)
PUBLIC :: InitBV
#endif /*defined(PARTICLES)*/
#endif /*USE_HDG*/
!===================================================================================================================================

CONTAINS

#if USE_HDG
!===================================================================================================================================
!> Create containers and communicators for each floating boundary condition where impacting charges are accumulated.
!>
!> 1.) Loop over all field BCs and check if the current processor is either the MPI root or has at least one of the BCs that
!>     contribute to the total floating boundary condition. If yes, then this processor is part of the communicator
!> 2.) Create Mapping from floating boundary condition BC index to field BC index
!> 3.) Create Mapping from field BC index to floating boundary condition BC index
!> 4.0) Check if field BC is on current proc (or MPI root)
!> 4.1.) Each processor loops over all of his elements
!> 4.2.) Loop over all compute-node elements (every processor loops over all of these elements)
!> 5.) Create MPI sub-communicators
!===================================================================================================================================
SUBROUTINE InitFPC()
! MODULES
USE MOD_Globals  ! ,ONLY: MPIRoot,iError,myrank,UNIT_stdOut,MPI_COMM_PICLAS
USE MOD_Preproc
USE MOD_Mesh_Vars          ,ONLY: nBCs,BoundaryType
USE MOD_Analyze_Vars       ,ONLY: DoFieldAnalyze
USE MOD_HDG_Vars           ,ONLY: UseFPC,FPC
USE MOD_Mesh_Vars          ,ONLY: nBCSides,BC
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
#if USE_MPI && defined(PARTICLES)
USE MOD_Mesh_Tools         ,ONLY: GetGlobalElemID
USE MOD_Globals            ,ONLY: ElementOnProc
USE MOD_Particle_Mesh_Vars ,ONLY: ElemInfo_Shared,BoundsOfElem_Shared,SideInfo_Shared
USE MOD_MPI_Shared_Vars    ,ONLY: nComputeNodeTotalElems
USE MOD_Mesh_Vars          ,ONLY: nElems, offsetElem
USE MOD_Particle_MPI_Vars  ,ONLY: halo_eps
#endif /*USE_MPI && defined(PARTICLES)*/
#if USE_MPI
USE MOD_MPI_Shared_Vars    ,ONLY: nComputeNodeProcessors,nProcessors_Global
#endif /*USE_MPI*/
USE MOD_Equation_Vars      ,ONLY: IniExactFunc
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
USE MOD_HDG_Readin         ,ONLY: ReadFPCDataFromH5
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, PARAMETER  :: BCTypeFPC = 20
INTEGER             :: BCType,BCState,iUniqueFPCBC
INTEGER             :: SideID,iBC
#if USE_MPI
INTEGER             :: color,WithSides
#if defined(PARTICLES)
INTEGER             :: iElem,iCNElem
REAL                :: iElemCenter(1:3),iGlobElemCenter(1:3)
REAL                :: iElemRadius,iGlobElemRadius
INTEGER             :: iGlobElem,BCIndex,iSide
#endif /*defined(PARTICLES)*/
#endif /*USE_MPI*/
CHARACTER(5)        :: hilf,hilf2
REAL                :: StartT,EndT
!===================================================================================================================================

! Get global number of FPC boundaries in [1:nBCs], they might belong to the same group (will be reduced to "nUniqueFPCBounds" below)
! FPC boundaries with the same BCState will be in the same group (electrically connected)
UseFPC = .FALSE.
FPC%nFPCBounds = 0
FPC%nUniqueFPCBounds = 0
DO iBC=1,nBCs
  BCType = BoundaryType(iBC,BC_TYPE)
  IF(BCType.NE.BCTypeFPC) CYCLE ! Skip non-FPC boundaries
  BCState = BoundaryType(iBC,BC_STATE) ! State is iFPC
  FPC%nFPCBounds=FPC%nFPCBounds+1
  IF(BCState.LE.0) CALL CollectiveStop(__STAMP__,' BCState for FPC must be >0! BCState=',IntInfo=BCState)
END DO

IF(FPC%nFPCBounds.EQ.0) RETURN ! Already determined in HDG initialization

UseFPC = .TRUE.

#if !(USE_PETSC)
CALL CollectiveStop(__STAMP__,'FPC model requires compilation with LIBS_USE_PETSC=ON')
#endif /*!(USE_PETSC)*/

GETTIME(StartT)
LBWRITE(UNIT_stdOut,'(A)')' | INIT FPC ...'

ALLOCATE(FPC%Group(1:FPC%nFPCBounds,3))
FPC%Group = 0 ! Initialize

! 1.) Loop over all field BCs and check if the current processor is either the MPI root or has at least one of the BCs that
! contribute to the total floating boundary condition. If yes, then this processor is part of the communicator
DO iBC=1,nBCs
  BCType = BoundaryType(iBC,BC_TYPE)
  IF(BCType.NE.BCTypeFPC) CYCLE ! Skip non-FPC boundaries
  BCState = BoundaryType(iBC,BC_STATE) ! State is iFPC
  WRITE(UNIT=hilf,FMT='(I0)') BCState
  WRITE(UNIT=hilf2,FMT='(I0)') FPC%nFPCBounds
  IF(BCState.GT.FPC%nFPCBounds) CALL CollectiveStop(__STAMP__,&
      'BCState='//TRIM(hilf)//' must be smaller or equal than the total number of '//TRIM(hilf2)//' FPC boundaries!')
  FPC%Group(BCState,1) = FPC%Group(BCState,1) + 1
  IF(FPC%Group(BCState,1).EQ.1) THEN
    FPC%nUniqueFPCBounds = FPC%nUniqueFPCBounds +1 ! Only count once
    FPC%Group(BCState,2) = FPC%nUniqueFPCBounds
  END IF
END DO

! Automatically activate surface model analyze flag
DoFieldAnalyze = .TRUE.

! 2.) Create Mapping from floating boundary condition BC index to BCState
ALLOCATE(FPC%BCState(FPC%nUniqueFPCBounds))
FPC%BCState = 0
DO iBC=1,nBCs
  BCType = BoundaryType(iBC,BC_TYPE)
  IF(BCType.NE.BCTypeFPC) CYCLE
  BCState = BoundaryType(iBC,BC_STATE) ! State is iFPC
  iUniqueFPCBC = FPC%Group(BCState,2)
  FPC%BCState(iUniqueFPCBC) = BCState
END DO

! Allocate the containers
! This container is not deallocated for MPIRoot when performing load balance (only root needs this info to write it to .csv)
IF(.NOT.ALLOCATED(FPC%Voltage))THEN
  ALLOCATE(FPC%Voltage(1:FPC%nUniqueFPCBounds))
  FPC%Voltage = 0.
END IF ! .NOT.ALLOCATED(FPC%Voltage)
ALLOCATE(FPC%VoltageProc(1:FPC%nUniqueFPCBounds))
FPC%VoltageProc = 0.
! This container is not deallocated for MPIRoot when performing load balance as this process updates the information on the new
! sub-communicator processes during load balance
IF(.NOT.ALLOCATED(FPC%Charge))THEN
  ALLOCATE(FPC%Charge(1:FPC%nUniqueFPCBounds))
  FPC%Charge = 0.
END IF ! .NOT.ALLOCATED(FPC%Charge)
ALLOCATE(FPC%ChargeProc(1:FPC%nUniqueFPCBounds))
FPC%ChargeProc = 0.

! Set initial value depending on IniExactFunc
SELECT CASE (IniExactFunc)
CASE(800,900,901,1000,1100) ! Dielectric slab on electrode (left) with plasma between slab and other electrode opposite
  ! Set initial value
  FPC%Charge(1)  = 1.0e-2*(GEO%ymaxglob-GEO%yminglob)*(GEO%zmaxglob-GEO%zminglob) ! C/m2
  FPC%Voltage(1) = 1.1293922903231239 ! V
CASE(801) ! Dielectric slab on electrode (left) with plasma between slab and other electrode opposite: 2D case
  ! Set initial value
  FPC%Charge(1)  = 5e-11*(1.0 - 0.05) ! C/m2
  FPC%Charge(2)  = 5e-11*(1.0 - 0.15) ! C/m2
  FPC%Charge(3)  = 5e-11*(1.0 - 0.25) ! C/m2
  FPC%Charge(4)  = 5e-11*(1.0 - 0.35) ! C/m2
  FPC%Charge(5)  = 5e-11*(1.0 - 0.45) ! C/m2
  FPC%Charge(6)  = 5e-11*(1.0 - 0.55) ! C/m2
  FPC%Charge(7)  = 5e-11*(1.0 - 0.65) ! C/m2
  FPC%Charge(8)  = 5e-11*(1.0 - 0.75) ! C/m2
  FPC%Charge(9)  = 5e-11*(1.0 - 0.85) ! C/m2
  FPC%Charge(10) = 5e-11*(1.0 - 0.95) ! C/m2
END SELECT

!! 3.) Create Mapping from field BC index to floating boundary condition BC index
!ALLOCATE(FPC%BCIDToFPCBCID(nBCs))
!FPC%BCIDToFPCBCID = -1
!DO iFPCBC = 1, FPC%NBoundaries
!  iBC = EDC%FieldBoundaries(iEDCBC)
!  EDC%BCIDToEDCBCID(iBC) = iEDCBC
!END DO ! iEDCBC = 1, EDC%NBoundaries

! Get processor-local number of FPC sides associated with each i-th FPC boundary
! Check local sides
DO SideID=1,nBCSides
  iBC    = BC(SideID)
  BCType = BoundaryType(iBC,BC_TYPE)
  IF(BCType.NE.BCTypeFPC) CYCLE ! Skip non-FPC boundaries
  BCState = BoundaryType(iBC,BC_STATE) ! BCState corresponds to iFPC
  FPC%Group(BCState,3) = FPC%Group(BCState,3) + 1
END DO ! SideID=1,nBCSides

#if USE_MPI
! 4.0) Check if field BC is on current proc (or MPI root)
ALLOCATE(FPC%BConProc(FPC%nUniqueFPCBounds))
FPC%BConProc = .FALSE.

! Check if single-node or multi-node run
IF (nComputeNodeProcessors.EQ.nProcessors_Global) THEN
  ! For single-node execution, simply add all processes to the communicators and do not bother measuring the distance as the gain
  ! in performance is negligible here
  FPC%BConProc = .TRUE.
ELSE
  ! Check local sides
  DO SideID=1,nBCSides
    iBC    = BC(SideID)
    BCType = BoundaryType(iBC,BC_TYPE)
    IF(BCType.NE.BCTypeFPC) CYCLE ! Skip non-FPC boundaries
    BCState = BoundaryType(iBC,BC_STATE) ! BCState corresponds to iFPC
    iUniqueFPCBC = FPC%Group(BCState,2)
    FPC%BConProc(iUniqueFPCBC) = .TRUE.
  END DO ! SideID=1,nBCSides
END IF ! nComputeNodeProcessors.EQ.nProcessors_Global

#if defined(PARTICLES)
  ! Check if all FPCs have already been found
IF(.NOT.(ALL(FPC%BConProc)))THEN

  ! Check whether this information has already been created before to skip the costly search below
  !CALL ReadFPCCommunicationFromH5()

    ! Particles might impact the FPC on another proc/node. Therefore check if a particle can travel from a local element to an
    ! element that has at least one side, which is an FPC
    ! 4.1.) Each processor loops over all of his elements
    iElemLoop: DO iElem = 1+offsetElem, nElems+offsetElem

      iElemCenter(1:3) = (/ SUM(BoundsOfElem_Shared(1:2,1,iElem)),&
                            SUM(BoundsOfElem_Shared(1:2,2,iElem)),&
                            SUM(BoundsOfElem_Shared(1:2,3,iElem)) /) / 2.
      iElemRadius = VECNORM3D ((/ BoundsOfElem_Shared(2,1,iElem)-BoundsOfElem_Shared(1,1,iElem),&
                                  BoundsOfElem_Shared(2,2,iElem)-BoundsOfElem_Shared(1,2,iElem),&
                                  BoundsOfElem_Shared(2,3,iElem)-BoundsOfElem_Shared(1,3,iElem) /) / 2.)

      ! 4.2.) Loop over all compute-node elements (every processor loops over all of these elements)
      ! Loop ALL compute-node elements (use global element index)
      iCNElemLoop: DO iCNElem = 1,nComputeNodeTotalElems
        iGlobElem = GetGlobalElemID(iCNElem)

        ! Skip my own elements as they have already been tested when the local sides are checked
        IF(ElementOnProc(iGlobElem)) CYCLE iCNElemLoop

        ! Check if one of the six sides of the compute-node element is a FPC
        ! Note that iSide is in the range of 1:nNonUniqueGlobalSides
        DO iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iGlobElem)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,iGlobElem)
          ! Get BC index of the global side index
          BCIndex = SideInfo_Shared(SIDE_BCID,iSide)
          ! Only check BC sides with BC index > 0
          IF(BCIndex.GT.0)THEN
            ! Get boundary type
            BCType = BoundaryType(BCIndex,BC_TYPE)
            ! Check if FPC has been found
            IF(BCType.EQ.BCTypeFPC)THEN

              ! Check if the BC can be reached
              iGlobElemCenter(1:3) = (/ SUM(BoundsOfElem_Shared(1:2,1,iGlobElem)),&
                                        SUM(BoundsOfElem_Shared(1:2,2,iGlobElem)),&
                                        SUM(BoundsOfElem_Shared(1:2,3,iGlobElem)) /) / 2.
              iGlobElemRadius = VECNORM3D ((/ BoundsOfElem_Shared(2,1,iGlobElem)-BoundsOfElem_Shared(1,1,iGlobElem),&
                                              BoundsOfElem_Shared(2,2,iGlobElem)-BoundsOfElem_Shared(1,2,iGlobElem),&
                                              BoundsOfElem_Shared(2,3,iGlobElem)-BoundsOfElem_Shared(1,3,iGlobElem) /) / 2.)

              ! check if compute-node element "iGlobElem" is within halo_eps of processor-local element "iElem"
            ! TODO: what about periodic vectors?
              IF (VECNORM3D( iElemCenter(1:3) - iGlobElemCenter(1:3) ) .LE. ( halo_eps + iElemRadius + iGlobElemRadius ) )THEN
                BCState = BoundaryType(BCIndex,BC_STATE) ! BCState corresponds to iFPC
                IF(BCState.LT.1) CALL abort(__STAMP__,'BCState cannot be <1',IntInfoOpt=BCState)
                iUniqueFPCBC = FPC%Group(BCState,2)
                ! Flag the i-th FPC
              FPC%BConProc(iUniqueFPCBC) = .TRUE.
                ! Check if all FPCs have been found -> exit complete loop
              IF(ALL(FPC%BConProc)) EXIT iElemLoop
                ! Go to next element
                CYCLE iCNElemLoop
              END IF ! VECNORM3D( ...
            END IF ! BCType.EQ.BCTypeFPC
          END IF ! BCIndex.GT.0
        END DO ! iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iGlobElem)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,iGlobElem)
      END DO iCNElemLoop ! iCNElem = 1,nComputeNodeTotalElems
    END DO iElemLoop ! iElem = 1, nElems
END IF ! .NOT.(ALL(FPC%BConProc))
#endif /*defined(PARTICLES)*/

#if USE_MPI
! Store FPC%BConProc information to skip the costly search above for the next run of the same setup. Store the info before adding
! the MPIRoot by force in the next step.
! CALL WriteFPCCommunicationToH5()

! Store similar to FPCDataHDF5
! 1.) nProcessors
! 2.) nFPCs
! 3.) halo_eps as integer value
! and check when reading that the value are unchanged or less procs, smaller halo_eps is still okay
#endif /*USE_MPI*/

! Storing the FPC info requires that the MPIRoot also takes part in the search even though it is not required as the MPIRoot is
! always part of the communicator anyway as the MPIRoot writes the FPC charge/potential information to the .csv file.
! The MPIRoot is part of all communicators
IF(MPIRoot) FPC%BConProc = .TRUE.

! 5.) Create MPI sub-communicators
ALLOCATE(FPC%COMM(FPC%nUniqueFPCBounds))
DO iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
  ! create new communicator
  color = MERGE(iUniqueFPCBC, MPI_UNDEFINED, FPC%BConProc(iUniqueFPCBC))

  ! set communicator id
  FPC%COMM(iUniqueFPCBC)%ID=iUniqueFPCBC

  ! create new emission communicator for floating boundary condition communication. Pass MPI_INFO_NULL as rank to follow the original ordering
  CALL MPI_COMM_SPLIT(MPI_COMM_PICLAS, color, 0, FPC%COMM(iUniqueFPCBC)%UNICATOR, iError)

  ! Find my rank on the shared communicator, comm size and proc name
  IF(FPC%BConProc(iUniqueFPCBC))THEN
    CALL MPI_COMM_RANK(FPC%COMM(iUniqueFPCBC)%UNICATOR, FPC%COMM(iUniqueFPCBC)%MyRank, iError)
    CALL MPI_COMM_SIZE(FPC%COMM(iUniqueFPCBC)%UNICATOR, FPC%COMM(iUniqueFPCBC)%nProcs, iError)

    ! inform about size of emission communicator
    IF (FPC%COMM(iUniqueFPCBC)%MyRank.EQ.0) THEN
#if USE_LOADBALANCE
      IF(.NOT.PerformLoadBalance)&
#endif /*USE_LOADBALANCE*/
          WRITE(UNIT_StdOut,'(A,I0,A,I0,A,I0)') ' Floating boundary condition: Emission-Communicator ',iUniqueFPCBC,' on ',&
              FPC%COMM(iUniqueFPCBC)%nProcs,' procs for BCState ',FPC%BCState(iUniqueFPCBC)
    END IF
  END IF ! FPC%BConProc(iUniqueFPCBC)
END DO ! iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
DEALLOCATE(FPC%BConProc)

! Get the number of procs that actually have a local BC side that is an FPC (required for voltage output to .csv)
! Procs might have zero FPC sides but are in the group because 1.) MPIRoot or 2.) the FPC is in the halo region
! Because only the MPI root process writes the .csv data, the information regarding the voltage on each FPC must be
! communicated with this process even though it might not be connected to each FPC boundary
DO iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
  ASSOCIATE( COMM => FPC%COMM(iUniqueFPCBC)%UNICATOR, nProcsWithSides => FPC%COMM(iUniqueFPCBC)%nProcsWithSides )
    IF(FPC%COMM(iUniqueFPCBC)%UNICATOR.NE.MPI_COMM_NULL)THEN
      ! Check if the current processor is actually connected to the FPC via a BC side
      IF(FPC%Group(FPC%BCState(iUniqueFPCBC),3).EQ.0)THEN
        WithSides = 0
      ELSE
        WithSides = 1
      END IF ! FPC%Group(FPC%BCState(iUniqueFPCBC),3).EQ.0
      ! Calculate the sum across the sub-communicator. Only the MPI root process needs this information
      IF(MPIRoot)THEN
        CALL MPI_REDUCE(WithSides, nProcsWithSides, 1 ,MPI_INTEGER, MPI_SUM, 0, COMM, iError)
        ! Sanity check
        IF(nProcsWithSides.EQ.0) CALL abort(__STAMP__,'Found FPC with no processors connected to it')
      ELSE
        CALL MPI_REDUCE(WithSides, 0              , 1 ,MPI_INTEGER, MPI_SUM, 0, COMM, IError)
      END IF ! MPIRoot
    END IF ! FPC%COMM(iUniqueFPCBC)%UNICATOR.NE.MPI_COMM_NULL
  END ASSOCIATE
END DO ! iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
#endif /*USE_MPI*/

! When restarting, load the deposited charge on each FPC from the .h5 state file
CALL ReadFPCDataFromH5()

GETTIME(EndT)
LBWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' | INIT FPC'
CALL DisplayMessageAndTime(EndT-StartT, 'DONE!')

END SUBROUTINE InitFPC


!===================================================================================================================================
!> Create containers and communicators for each electric potential boundary condition where impacting charges are removed and
!> subsequently an electric potential is created.
!>
!> 1.) Loop over all field BCs and check if the current processor is either the MPI root or has at least one of the BCs that
!>     contribute to the total electric potential boundary condition. If yes, then this processor is part of the communicator
!> 2.) Create Mapping from electric potential boundary condition BC index to field BC index
!> 3.) Check if field BC is on current proc (or MPI root)
!> 3.1) Each processor loops over all of his elements
!> 3.2) Loop over all compute-node elements (every processor loops over all of these elements)
!> 4.) Create MPI sub-communicators
!===================================================================================================================================
SUBROUTINE InitEPC()
! MODULES
USE MOD_Globals  ! ,ONLY: MPIRoot,iError,myrank,UNIT_stdOut,MPI_COMM_PICLAS
USE MOD_Preproc
USE MOD_Mesh_Vars          ,ONLY: nBCs,BoundaryType
USE MOD_Analyze_Vars       ,ONLY: DoFieldAnalyze
USE MOD_HDG_Vars           ,ONLY: UseEPC,EPC
USE MOD_Mesh_Vars          ,ONLY: nBCSides,BC
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
#if USE_MPI && defined(PARTICLES)
USE MOD_Mesh_Tools         ,ONLY: GetGlobalElemID
USE MOD_Globals            ,ONLY: ElementOnProc
USE MOD_Particle_Mesh_Vars ,ONLY: ElemInfo_Shared,BoundsOfElem_Shared,SideInfo_Shared
USE MOD_MPI_Shared_Vars    ,ONLY: nComputeNodeTotalElems
USE MOD_Particle_MPI_Vars  ,ONLY: halo_eps
USE MOD_Mesh_Vars          ,ONLY: nElems, offsetElem
#endif /*USE_MPI && defined(PARTICLES)*/
USE MOD_ReadInTools        ,ONLY: GETREALARRAY
USE MOD_HDG_Readin         ,ONLY: ReadEPCDataFromH5
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, PARAMETER  :: BCTypeEPC = 8
INTEGER             :: BCType,BCState,iUniqueEPCBC
INTEGER             :: SideID,iBC
#if USE_MPI
INTEGER             :: color,WithSides
LOGICAL,ALLOCATABLE :: BConProc(:)
#if defined(PARTICLES)
INTEGER             :: iElem,iCNElem
REAL                :: iElemCenter(1:3),iGlobElemCenter(1:3)
REAL                :: iElemRadius,iGlobElemRadius
INTEGER             :: iGlobElem,BCIndex,iSide
#endif /*defined(PARTICLES)*/
#endif /*USE_MPI*/
CHARACTER(5)        :: hilf,hilf2
!===================================================================================================================================

! Get global number of EPC boundaries in [1:nBCs], they might belong to the same group (will be reduced to "nUniqueEPCBounds" below)
! EPC boundaries with the same BCState will be in the same group (electrically connected)
UseEPC = .FALSE.
EPC%nEPCBounds = 0
EPC%nUniqueEPCBounds = 0
DO iBC=1,nBCs
  BCType = BoundaryType(iBC,BC_TYPE)
  IF(BCType.NE.BCTypeEPC) CYCLE ! Skip non-EPC boundaries
  BCState = BoundaryType(iBC,BC_STATE) ! State is iEPC
  EPC%nEPCBounds=EPC%nEPCBounds+1
  IF(BCState.LE.0) CALL CollectiveStop(__STAMP__,' BCState for EPC must be >0! BCState=',IntInfo=BCState)
END DO

IF(EPC%nEPCBounds.EQ.0) RETURN ! Already determined in HDG initialization

UseEPC = .TRUE.

ALLOCATE(EPC%Group(1:EPC%nEPCBounds,3))
EPC%Group = 0 ! Initialize

! 1.) Loop over all field BCs and check if the current processor is either the MPI root or has at least one of the BCs that
! contribute to the total electric potential boundary condition. If yes, then this processor is part of the communicator
DO iBC=1,nBCs
  BCType = BoundaryType(iBC,BC_TYPE)
  IF(BCType.NE.BCTypeEPC) CYCLE ! Skip non-EPC boundaries
  BCState = BoundaryType(iBC,BC_STATE) ! State is iEPC
  WRITE(UNIT=hilf,FMT='(I0)') BCState
  WRITE(UNIT=hilf2,FMT='(I0)') EPC%nEPCBounds
  IF(BCState.GT.EPC%nEPCBounds) CALL CollectiveStop(__STAMP__,&
      'BCState='//TRIM(hilf)//' must be smaller or equal than the total number of '//TRIM(hilf2)//' EPC boundaries!')
  EPC%Group(BCState,1) = EPC%Group(BCState,1) + 1
  IF(EPC%Group(BCState,1).EQ.1) THEN
    EPC%nUniqueEPCBounds = EPC%nUniqueEPCBounds +1 ! Only count once
    EPC%Group(BCState,2) = EPC%nUniqueEPCBounds
  END IF
END DO

! Automatically activate surface model analyze flag
DoFieldAnalyze = .TRUE.

! Read resistances for each unique EPC
EPC%Resistance  = GETREALARRAY('EPC-Resistance',EPC%nUniqueEPCBounds)

! 2.) Create Mapping from electric potential boundary condition BC index to BCState
ALLOCATE(EPC%BCState(EPC%nUniqueEPCBounds))
EPC%BCState = 0
DO iBC=1,nBCs
  BCType = BoundaryType(iBC,BC_TYPE)
  IF(BCType.NE.BCTypeEPC) CYCLE
  BCState = BoundaryType(iBC,BC_STATE) ! State is iEPC
  iUniqueEPCBC = EPC%Group(BCState,2)
  EPC%BCState(iUniqueEPCBC) = BCState
END DO

! Allocate the containers
! This container is not deallocated for MPIRoot when performing load balance (only root needs this info to write it to .csv)
IF(.NOT.ALLOCATED(EPC%Voltage))THEN
  ALLOCATE(EPC%Voltage(1:EPC%nUniqueEPCBounds))
  EPC%Voltage = 0.
END IF ! .NOT.ALLOCATED(EPC%Voltage)
ALLOCATE(EPC%VoltageProc(1:EPC%nUniqueEPCBounds))
EPC%VoltageProc = 0.
! This container is not deallocated for MPIRoot when performing load balance as this process updates the information on the new
! sub-communicator processes during load balance
IF(.NOT.ALLOCATED(EPC%Charge))THEN
  ALLOCATE(EPC%Charge(1:EPC%nUniqueEPCBounds))
  EPC%Charge = 0.
END IF ! .NOT.ALLOCATED(EPC%Charge)
ALLOCATE(EPC%ChargeProc(1:EPC%nUniqueEPCBounds))
EPC%ChargeProc = 0.

! Get processor-local number of EPC sides associated with each i-th EPC boundary
! Check local sides
DO SideID=1,nBCSides
  iBC    = BC(SideID)
  BCType = BoundaryType(iBC,BC_TYPE)
  IF(BCType.NE.BCTypeEPC) CYCLE ! Skip non-EPC boundaries
  BCState = BoundaryType(iBC,BC_STATE) ! BCState corresponds to iEPC
  EPC%Group(BCState,3) = EPC%Group(BCState,3) + 1
END DO ! SideID=1,nBCSides

#if USE_MPI
! 3) Check if field BC is on current proc (or MPI root)
ALLOCATE(BConProc(EPC%nUniqueEPCBounds))
BConProc = .FALSE.
IF(MPIRoot)THEN
  BConProc = .TRUE.
ELSE

  ! Check local sides
  DO SideID=1,nBCSides
    iBC    = BC(SideID)
    BCType = BoundaryType(iBC,BC_TYPE)
    IF(BCType.NE.BCTypeEPC) CYCLE ! Skip non-EPC boundaries
    BCState = BoundaryType(iBC,BC_STATE) ! BCState corresponds to iEPC
    iUniqueEPCBC = EPC%Group(BCState,2)
    BConProc(iUniqueEPCBC) = .TRUE.
  END DO ! SideID=1,nBCSides

#if defined(PARTICLES)
  ! Check if all FPCs have already been found
  IF(.NOT.(ALL(BConProc)))THEN
    ! Particles might impact the EPC on another proc/node. Therefore check if a particle can travel from a local element to an
    ! element that has at least one side, which is an EPC
    ! 3.1) Each processor loops over all of his elements
    iElemLoop: DO iElem = 1+offsetElem, nElems+offsetElem

      iElemCenter(1:3) = (/ SUM(BoundsOfElem_Shared(1:2,1,iElem)),&
                            SUM(BoundsOfElem_Shared(1:2,2,iElem)),&
                            SUM(BoundsOfElem_Shared(1:2,3,iElem)) /) / 2.
      iElemRadius = VECNORM3D ((/ BoundsOfElem_Shared(2,1,iElem)-BoundsOfElem_Shared(1,1,iElem),&
                                  BoundsOfElem_Shared(2,2,iElem)-BoundsOfElem_Shared(1,2,iElem),&
                                  BoundsOfElem_Shared(2,3,iElem)-BoundsOfElem_Shared(1,3,iElem) /) / 2.)

      ! 3.2) Loop over all compute-node elements (every processor loops over all of these elements)
      ! Loop ALL compute-node elements (use global element index)
      iCNElemLoop: DO iCNElem = 1,nComputeNodeTotalElems
        iGlobElem = GetGlobalElemID(iCNElem)

        ! Skip my own elements as they have already been tested when the local sides are checked
        IF(ElementOnProc(iGlobElem)) CYCLE iCNElemLoop

        ! Check if one of the six sides of the compute-node element is a EPC
        ! Note that iSide is in the range of 1:nNonUniqueGlobalSides
        DO iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iGlobElem)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,iGlobElem)
          ! Get BC index of the global side index
          BCIndex = SideInfo_Shared(SIDE_BCID,iSide)
          ! Only check BC sides with BC index > 0
          IF(BCIndex.GT.0)THEN
            ! Get boundary type
            BCType = BoundaryType(BCIndex,BC_TYPE)
            ! Check if EPC has been found
            IF(BCType.EQ.BCTypeEPC)THEN

              ! Check if the BC can be reached
              iGlobElemCenter(1:3) = (/ SUM(BoundsOfElem_Shared(1:2,1,iGlobElem)),&
                                        SUM(BoundsOfElem_Shared(1:2,2,iGlobElem)),&
                                        SUM(BoundsOfElem_Shared(1:2,3,iGlobElem)) /) / 2.
              iGlobElemRadius = VECNORM3D ((/ BoundsOfElem_Shared(2,1,iGlobElem)-BoundsOfElem_Shared(1,1,iGlobElem),&
                                              BoundsOfElem_Shared(2,2,iGlobElem)-BoundsOfElem_Shared(1,2,iGlobElem),&
                                              BoundsOfElem_Shared(2,3,iGlobElem)-BoundsOfElem_Shared(1,3,iGlobElem) /) / 2.)

              ! check if compute-node element "iGlobElem" is within halo_eps of processor-local element "iElem"
              IF (VECNORM3D( iElemCenter(1:3) - iGlobElemCenter(1:3) ) .LE. ( halo_eps + iElemRadius + iGlobElemRadius ) )THEN
                BCState = BoundaryType(BCIndex,BC_STATE) ! BCState corresponds to iEPC
                IF(BCState.LT.1) CALL abort(__STAMP__,'BCState cannot be <1',IntInfoOpt=BCState)
                iUniqueEPCBC = EPC%Group(BCState,2)
                ! Flag the i-th EPC
                BConProc(iUniqueEPCBC) = .TRUE.
                ! Check if all FPCs have been found -> exit complete loop
                IF(ALL(BConProc)) EXIT iElemLoop
                ! Go to next element
                CYCLE iCNElemLoop
              END IF ! VECNORM3D( ...
            END IF ! BCType.EQ.BCTypeEPC
          END IF ! BCIndex.GT.0
        END DO ! iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iGlobElem)+1,ElemInfo_Shared(ELEM_LASTSIDEIND,iGlobElem)
      END DO iCNElemLoop ! iCNElem = 1,nComputeNodeTotalElems
    END DO iElemLoop ! iElem = 1, nElems
  END IF ! .NOT.(ALL(BConProc))
#endif /*defined(PARTICLES)*/

END IF ! MPIRoot

! 4.) Create MPI sub-communicators
ALLOCATE(EPC%COMM(EPC%nUniqueEPCBounds))
DO iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
  ! Create new communicator
  color = MERGE(iUniqueEPCBC, MPI_UNDEFINED, BConProc(iUniqueEPCBC))

  ! Set communicator id
  EPC%COMM(iUniqueEPCBC)%ID=iUniqueEPCBC

  ! Create new emission communicator for Electric potential boundary condition communication.
  ! Pass MPI_INFO_NULL as rank to follow the original ordering
  CALL MPI_COMM_SPLIT(MPI_COMM_PICLAS, color, 0, EPC%COMM(iUniqueEPCBC)%UNICATOR, iError)

  ! Find my rank on the shared communicator, comm size and proc name
  IF(BConProc(iUniqueEPCBC))THEN
    CALL MPI_COMM_RANK(EPC%COMM(iUniqueEPCBC)%UNICATOR, EPC%COMM(iUniqueEPCBC)%MyRank, iError)
    CALL MPI_COMM_SIZE(EPC%COMM(iUniqueEPCBC)%UNICATOR, EPC%COMM(iUniqueEPCBC)%nProcs, iError)

    ! inform about size of emission communicator
    IF (EPC%COMM(iUniqueEPCBC)%MyRank.EQ.0) THEN
#if USE_LOADBALANCE
      IF(.NOT.PerformLoadBalance)&
#endif /*USE_LOADBALANCE*/
          WRITE(UNIT_StdOut,'(A,I0,A,I0,A,I0)') ' Electric potential boundary condition: Emission-Communicator ',iUniqueEPCBC,' on ',&
              EPC%COMM(iUniqueEPCBC)%nProcs,' procs for BCState ',EPC%BCState(iUniqueEPCBC)
    END IF
  END IF ! BConProc(iUniqueEPCBC)
END DO ! iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
DEALLOCATE(BConProc)

! Get the number of procs that actually have a local BC side that is an EPC (required for voltage output to .csv)
! Procs might have zero EPC sides but are in the group because 1.) MPIRoot or 2.) the EPC is in the halo region
! Because only the MPI root process writes the .csv data, the information regarding the voltage on each EPC must be
! communicated with this process even though it might not be connected to each EPC boundary
DO iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
  ASSOCIATE( COMM => EPC%COMM(iUniqueEPCBC)%UNICATOR, nProcsWithSides => EPC%COMM(iUniqueEPCBC)%nProcsWithSides )
    IF(EPC%COMM(iUniqueEPCBC)%UNICATOR.NE.MPI_COMM_NULL)THEN
      ! Check if the current processor is actually connected to the EPC via a BC side
      IF(EPC%Group(EPC%BCState(iUniqueEPCBC),3).EQ.0)THEN
        WithSides = 0
      ELSE
        WithSides = 1
      END IF ! EPC%Group(EPC%BCState(iUniqueEPCBC),3).EQ.0
      ! Calculate the sum across the sub-communicator. Only the MPI root process needs this information
      IF(MPIRoot)THEN
        CALL MPI_REDUCE(WithSides, nProcsWithSides, 1 ,MPI_INTEGER, MPI_SUM, 0, COMM, iError)
        ! Sanity check
        IF(nProcsWithSides.EQ.0) CALL abort(__STAMP__,'Found EPC with no processors connected to it')
      ELSE
        CALL MPI_REDUCE(WithSides, 0              , 1 ,MPI_INTEGER, MPI_SUM, 0, COMM, IError)
      END IF ! MPIRoot
    END IF ! EPC%COMM(iUniqueEPCBC)%UNICATOR.NE.MPI_COMM_NULL
  END ASSOCIATE
END DO ! iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
#endif /*USE_MPI*/

! When restarting, load the deposited charge on each EPC from the .h5 state file
CALL ReadEPCDataFromH5()

END SUBROUTINE InitEPC


#if defined(PARTICLES)
!===================================================================================================================================
!> Create containers and communicators for each bias-voltage boundary condition where impacting charges are removed and
!> subsequently an electric potential is created (the particle communication is part of the BPO analysis and required here).
!>
!> 1.) Activate bias voltage and check number of boundaries
!> 2.) Get bias voltage parameters
!> 3.) Check if actual bias voltage BC is on current process (or MPI root)
!> 4.) Create MPI sub-communicators
!===================================================================================================================================
SUBROUTINE InitBV()
! MODULES
USE MOD_Globals                   ,ONLY: CollectiveStop,UNIT_StdOut
USE MOD_ReadInTools               ,ONLY: GETLOGICAL,GETREAL,GETINT,GETINTARRAY
USE MOD_Particle_Boundary_Vars    ,ONLY: PartBound
USE MOD_Mesh_Vars                 ,ONLY: nBCs,BoundaryType,BoundaryName
USE MOD_HDG_Vars                  ,ONLY: UseBiasVoltage,BiasVoltage
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: CalcBoundaryParticleOutput,BPO
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
#if USE_MPI
USE MOD_Globals                   ,ONLY: IERROR,MPI_COMM_NULL,MPI_DOUBLE_PRECISION,MPI_COMM_PICLAS,MPI_INFO_NULL,MPI_UNDEFINED,MPIRoot
USE MOD_Mesh_Vars                 ,ONLY: nBCSides,BC
#endif /*USE_MPI*/
USE MOD_HDG_Readin                ,ONLY: ReadBVDataFromH5
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER, PARAMETER :: BCTypeBV(1:3) = (/50,51,52/) ! BCType which allows bias voltage control
!                                                  ! 50: pure DC potential
!                                                  ! 51: cos(wt) function with DC bias
!                                                  ! 52: cos(wt) function with DC bias + coupled power for AC potential adjustment
INTEGER             :: BCType,BVBoundaries,BCState,iBoundary
INTEGER             :: iBC,iPBC
#if USE_MPI
INTEGER             :: color,SideID
LOGICAL             :: BConProc
#endif /*USE_MPI*/
!===================================================================================================================================

!> 1.) Activate bias voltage and check number of boundaries
! Activate model
UseBiasVoltage = GETLOGICAL('UseBiasVoltage')

! Count the number of boundaries that allow bias voltage
BVBoundaries = 0
DO iBC=1,nBCs
  BCType = BoundaryType(iBC,BC_TYPE)
  IF(.NOT.ANY(BCType.EQ.BCTypeBV)) CYCLE ! Skip other boundaries
  BVBoundaries = BVBoundaries + 1
END DO

! Skip the following if bias voltage is not active
IF(.NOT.UseBiasVoltage)THEN
  ! Sanity check before returning: bias voltage BCs cannot be used without activating the bias voltage model
  IF(BVBoundaries.GT.0) CALL CollectiveStop(__STAMP__,' Bias voltage BCs require UseBiasVoltage=T!')

  ! Exit this subroutine
  RETURN
END IF

! CalcBoundaryParticleOutput=T and boundaries must be set correctly
IF(.NOT.CalcBoundaryParticleOutput) CALL CollectiveStop(__STAMP__,' UseBiasVoltage=T requires CalcBoundaryParticleOutput=T!')

! Check the number of boundaries that allow bias voltage: Must be exactly 1
IF(BVBoundaries.NE.1) CALL CollectiveStop(__STAMP__,' UseBiasVoltage=T requires exactly one boundary with this feature!')

!> 2.) Get bias voltage parameters
BiasVoltage%NPartBoundaries = GETINT('BiasVoltage-NPartBoundaries')
BiasVoltage%PartBoundaries  = GETINTARRAY('Biasvoltage-PartBoundaries',biasvoltage%npartboundaries)
BiasVoltage%Frequency       = GETReal('BiasVoltage-Frequency')
BiasVoltage%Delta           = GETReal('BiasVoltage-Delta')
#if USE_LOADBALANCE
! Do not nullify during load balance in order to keep the old value on the MPIRoot
IF((.NOT.PerformLoadBalance).OR.(.NOT.MPIRoot))THEN
#endif /*USE_LOADBALANCE*/
  BiasVoltage%BVData = 0.
  ! Update time
  IF(BiasVoltage%Frequency.GT.0.0) BiasVoltage%BVData(3) = 1.0/BiasVoltage%Frequency
#if USE_LOADBALANCE
END IF
#endif /*USE_LOADBALANCE*/

IF(BiasVoltage%NPartBoundaries.LT.1) CALL CollectiveStop(__STAMP__,' UseBiasVoltage=T requires one or more particle boundaries!')

DO iBoundary=1,BiasVoltage%NPartBoundaries
  iPBC = BiasVoltage%PartBoundaries(iBoundary)
  IF(.NOT.ANY(iPBC.EQ.BPO%PartBoundaries(:))) CALL CollectiveStop(__STAMP__,&
      'One of Biasvoltage-PartBoundaries not defined in any BPO-PartBoundaries')
  iBC = PartBound%MapToFieldBC(iPBC)
  IF(iBC.GT.SIZE(BoundaryName)) CALL CollectiveStop(__STAMP__,'BiasVoltage-PartBoundaries BC index maps to wrong field BCID= ',&
    IntInfo=iBC)
  BCType  = BoundaryType(iBC,BC_TYPE)
  BCState = BoundaryType(iBC,BC_STATE)
  LBWRITE(UNIT_stdOut,'(A,I0,A,I0,A)') ' Activated bias voltage by collecting currents from ['//TRIM(BoundaryName(iBC))&
      //'] with BCType [',BCType,'] and BCState [',BCState,']'
END DO

#if USE_MPI
!> 3.) Check if actual bias voltage BC is on current process (or MPI root)
BConProc = .FALSE.
IF(MPIRoot)THEN
  BConProc = .TRUE.
ELSE
  ! Check local sides
  DO SideID=1,nBCSides
    iBC    = BC(SideID)
    BCType = BoundaryType(iBC,BC_TYPE)
    IF(.NOT.ANY(BCType.EQ.BCTypeBV)) CYCLE ! Skip other boundaries
    BConProc = .TRUE.
  END DO ! SideID=1,nBCSides
END IF ! MPIRoot

! 4.) Create MPI sub-communicators
! Create new communicator
color = MERGE(BVBoundaries, MPI_UNDEFINED, BConProc)

! Set communicator id
BiasVoltage%COMM%ID = BVBoundaries

! Create new emission communicator for electric potential boundary condition communication. Pass MPI_INFO_NULL as rank to follow the original ordering
CALL MPI_COMM_SPLIT(MPI_COMM_PICLAS, color, 0, BiasVoltage%COMM%UNICATOR, iError)

! Find my rank on the shared communicator, comm size and process name
IF(BConProc)THEN
  CALL MPI_COMM_RANK(BiasVoltage%COMM%UNICATOR, BiasVoltage%COMM%MyRank, iError)
  CALL MPI_COMM_SIZE(BiasVoltage%COMM%UNICATOR, BiasVoltage%COMM%nProcs, iError)

  ! Inform about size of emission communicator
  IF (BiasVoltage%COMM%MyRank.EQ.0) THEN
#if USE_LOADBALANCE
    IF(.NOT.PerformLoadBalance)&
#endif /*USE_LOADBALANCE*/
        WRITE(UNIT_StdOut,'(A,I0,A,I0)') ' Bias voltage communicator on ',BiasVoltage%COMM%nProcs,' procs for BCState ',BCState
  END IF
END IF ! BConProc
#endif /*USE_MPI*/

! When restarting, load the deposited charge on each EPC from the .h5 state file
CALL ReadBVDataFromH5()

END SUBROUTINE InitBV
#endif /*defined(PARTICLES)*/
#endif /*USE_HDG*/

END MODULE MOD_HDG_Init
