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

MODULE MOD_Mesh_ReadIn
!===================================================================================================================================
!> \brief Module containing routines to read the mesh and BCs from a HDF5 file
!>
!> This module contains the following routines related to mesh IO
!> - parallel HDF5-based mesh IO
!> - readin of mesh coordinates and connectivity
!> - readin of boundary conditions
!===================================================================================================================================
! MODULES
USE MOD_HDF5_Input
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
INTEGER              :: nNodeIDs

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE ReadMesh
  MODULE PROCEDURE ReadMesh
END INTERFACE

INTERFACE Qsort1Int
  MODULE PROCEDURE Qsort1Int
END INTERFACE

INTERFACE INVMAP
  MODULE PROCEDURE INVMAP
END INTERFACE

INTERFACE FinalizeMeshReadin
  MODULE PROCEDURE FinalizeMeshReadin
END INTERFACE

PUBLIC :: FinalizeMeshReadin
PUBLIC::ReadMesh,Qsort1Int,INVMAP
!===================================================================================================================================

CONTAINS

SUBROUTINE ReadBCs()
!===================================================================================================================================
!> This module will read boundary conditions from the HDF5 mesh file and from the parameter file.
!> The parameters defined in the mesh file can be overridden by those defined in the parameter file, by specifying the parameters
!> name and a new boundary condition set: a user-defined boundary condition consists of a type and a state.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars        ,ONLY: BoundaryName,BoundaryType,nBCs,nUserBCs
#if USE_HDG
USE MOD_Mesh_Vars        ,ONLY: ChangedPeriodicBC
#endif /*USE_HDG*/
USE MOD_ReadInTools      ,ONLY: GETINTARRAY,CountOption,GETSTR
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
LOGICAL,ALLOCATABLE            :: UserBCFound(:)
LOGICAL                        :: NameCheck,LengthCheck
CHARACTER(LEN=255),ALLOCATABLE :: BCNames(:)
INTEGER, ALLOCATABLE           :: BCMapping(:),BCType(:,:)
INTEGER                        :: iBC,iUserBC,OriginalBC,NewBC
INTEGER                        :: Offset=0 ! Every process reads all BCs
!===================================================================================================================================
! read in boundary conditions from ini file, will overwrite BCs from meshfile!
nUserBCs = CountOption('BoundaryName')
IF(nUserBCs.GT.0)THEN
  ALLOCATE(BoundaryName(1:nUserBCs))
  ALLOCATE(BoundaryType(1:nUserBCs,2))
  DO iBC=1,nUserBCs
    BoundaryName(iBC)   = GETSTR('BoundaryName')
    BoundaryType(iBC,:) = GETINTARRAY('BoundaryType',2) !(/Type,State/)
  END DO
END IF !nUserBCs>0

! Read boundary names from data file
CALL GetDataSize(File_ID,'BCNames',nDims,HSize)
CHECKSAFEINT(HSize(1),4)
nBCs=INT(HSize(1),4)
DEALLOCATE(HSize)
ALLOCATE(BCNames(nBCs))
ALLOCATE(BCMapping(nBCs))
ALLOCATE(UserBCFound(nUserBCs))

! Associate construct for integer KIND=8 possibility
ASSOCIATE ( nBCs   => INT(nBCs,IK)   ,&
            Offset => INT(Offset,IK) )
  CALL ReadArray('BCNames',1,(/nBCs/),Offset,1,StrArray=BCNames)  ! Type is a dummy type only
END ASSOCIATE
! User may have redefined boundaries in the ini file. So we have to create mappings for the boundaries.
BCMapping=0
UserBCFound=.FALSE.
IF(nUserBCs.GT.0)THEN
  DO iBC=1,nBCs
    DO iUserBC=1,nUserBCs
      ! Check if BoundaryName(iUserBC) is a substring of BCNames(iBC)
      NameCheck = INDEX(TRIM(BCNames(iBC)),TRIM(BoundaryName(iUserBC))).NE.0
      ! Check if both strings have equal length
      LengthCheck = LEN(TRIM(BCNames(iBC))).EQ.LEN(TRIM(BoundaryName(iUserBC)))
      ! Check if both strings are equal (length has to be checked because index checks for substrings!)
      IF(NameCheck.AND.LengthCheck)THEN
        BCMapping(iBC)=iUserBC
        UserBCFound(iUserBC)=.TRUE.
      END IF ! NameCheck.AND.LengthCheck
    END DO
  END DO
END IF

! Check if all BCs were found
DO iUserBC=1,nUserBCs
  IF(.NOT.UserBCFound(iUserBC)) CALL Abort(__STAMP__,'Boundary condition in parameter file not found: '//TRIM(BoundaryName(iUserBC)))
END DO
DEALLOCATE(UserBCFound)

! Read boundary types from data file
CALL GetDataSize(File_ID,'BCType',nDims,HSize)
IF((HSize(1).NE.4).OR.(HSize(2).NE.nBCs)) STOP 'Problem in readBC'
DEALLOCATE(HSize)
ALLOCATE(BCType(4,nBCs))
offset=0

! Associate construct for integer KIND=8 possibility
ASSOCIATE ( nBCs   => INT(nBCs,IK)   ,&
            Offset => INT(Offset,IK) )
  CALL ReadArray('BCType',2,(/4_IK,nBCs/),Offset,1,IntegerArray_i4=BCType)
END ASSOCIATE
! Now apply boundary mappings
#if USE_HDG
ChangedPeriodicBC=.FALSE. ! set true if BCs are changed from periodic to non-periodic
#endif /*USE_HDG*/
IF(nUserBCs .GT. 0)THEN
  LBWRITE(Unit_StdOut,'(A)')' REMAPPING BOUNDARY CONDITIONS...'
  DO iBC=1,nBCs
    IF(BCMapping(iBC).NE.0)THEN
      ! Compare new and original BC type (from mesh file)
      OriginalBC = BoundaryType(BCMapping(iBC),1)
      NewBC      = BCType(1,iBC)
      ! non-periodic to periodic
      IF((OriginalBC.EQ.1).AND.(NewBC.NE.1)) CALL abort(__STAMP__,'Remapping non-periodic to periodic BCs is not possible!')
#if USE_HDG
      ! periodic to non-periodic
      IF((NewBC.EQ.1).AND.(OriginalBC.NE.1))THEN
        ChangedPeriodicBC=.TRUE.
        ! Currently, remapping periodic to non-periodic BCs is not allowed. TODO: implement nGlobalUniqueSides determination.
        CALL abort(__STAMP__,'Remapping periodic to non-periodic BCs is currently not possible for HDG (changes nGlobalUniqueSides)')
      END IF
#endif /*USE_HDG*/
      ! Output
      LBWRITE(Unit_StdOut,'(A,A50,A,I4,I4,A,I4,I4)') ' |     Boundary in HDF file found |  ',TRIM(BCNames(iBC)), &
                                      ' was ', NewBC,BCType(3,iBC), ' is set to ',BoundaryType(BCMapping(iBC),1:2)
      BCType(1,iBC) = BoundaryType(BCMapping(iBC),BC_TYPE)
      BCType(3,iBC) = BoundaryType(BCMapping(iBC),BC_STATE)
    END IF
  END DO
END IF
IF(ALLOCATED(BoundaryName)) DEALLOCATE(BoundaryName)
IF(ALLOCATED(BoundaryType)) DEALLOCATE(BoundaryType)
ALLOCATE(BoundaryName(nBCs))
ALLOCATE(BoundaryType(nBCs,3))
BoundaryName = BCNames
BoundaryType(:,BC_TYPE)  = BCType(1,:)
BoundaryType(:,BC_STATE) = BCType(3,:)
BoundaryType(:,BC_ALPHA) = BCType(4,:)
LBWRITE(UNIT_StdOut,'(132("."))')
LBWRITE(Unit_StdOut,'(A,A15,A20,A10,A10,A10)')' BOUNDARY CONDITIONS','|','Name','Type','State','Alpha'
DO iBC=1,nBCs
  LBWRITE(*,'(A,A33,A20,I10,I10,I10)')' |','|',TRIM(BoundaryName(iBC)),BoundaryType(iBC,:)
END DO
LBWRITE(UNIT_StdOut,'(132("."))')
DEALLOCATE(BCNames,BCType,BCMapping)
END SUBROUTINE ReadBCs


SUBROUTINE ReadMesh(FileString,ReadNodes)
!===================================================================================================================================
!> This subroutine reads the mesh from the HDF5 mesh file. The connectivity and further relevant information as flips
!> (i.e. the orientation of sides towards each other) is already contained in the mesh file.
!> For parallel computations the number of elements will be distributed equally onto all processors and each processor only reads
!> its own subset of the mesh.
!> For a documentation of the mesh format see the documentation provided with HOPR (hopr-project.org)
!> The arrays ElemInfo, SideInfo and NodeCoords are read, alongside with the boundary condition data.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars         ,ONLY: ReadMeshWallTime
USE MOD_IO_HDF5
USE MOD_Mesh_Vars            ,ONLY: tElem,tSide
USE MOD_Mesh_Vars            ,ONLY: NGeo
USE MOD_Mesh_Vars            ,ONLY: NodeCoords
USE MOD_Mesh_Vars            ,ONLY: offsetElem,nElems,nGlobalElems,nNodes
USE MOD_Mesh_Vars            ,ONLY: nSides,nInnerSides,nBCSides,nMPISides,nAnalyzeSides,nGlobalMortarSides
USE MOD_Mesh_Vars            ,ONLY: nMortarSides
USE MOD_Mesh_Vars            ,ONLY: nGlobalUniqueSidesFromMesh
USE MOD_Mesh_Vars            ,ONLY: useCurveds
USE MOD_Mesh_Vars            ,ONLY: BoundaryType
USE MOD_Mesh_Vars            ,ONLY: MeshInitIsDone
USE MOD_Mesh_Vars            ,ONLY: Elems!,Nodes
USE MOD_Mesh_Vars            ,ONLY: GETNEWELEM,GETNEWSIDE
USE MOD_Mesh_Vars            ,ONLY: ElemInfo,SideInfo
USE MOD_Particle_Mesh_Vars   ,ONLY: nComputeNodeElems,nNonUniqueGlobalSides,nNonUniqueGlobalNodes
#if USE_MPI
USE MOD_MPI_Vars             ,ONLY: nMPISides_Proc,nNbProcs,NbProc,offsetElemMPI
USE MOD_LoadBalance_Tools    ,ONLY: DomainDecomposition
USE MOD_MPI_Shared_Vars      ,ONLY: ComputeNodeRootRank,nComputeNodeProcessors
#endif /*USE_MPI*/
#ifdef PARTICLES
USE MOD_Particle_Mesh_Readin, ONLY: ReadMeshBasics
USE MOD_Particle_Mesh_Readin, ONLY: ReadMeshSideNeighbors
USE MOD_Particle_Mesh_Readin, ONLY: StartCommunicateMeshReadin,FinishCommunicateMeshReadin
USE MOD_Particle_Vars        ,ONLY: VarTimeStep
USE MOD_LoadBalance_Vars     ,ONLY: nPartsPerElem
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars     ,ONLY: nDeposPerElem,nSurfacePartsPerElem,nTracksPerElem,nPartsPerBCElem,nSurfacefluxPerElem
! Restart without HDF5
USE MOD_Particle_Mesh_Vars   ,ONLY: ElemInfo_Shared,SideInfo_Shared,NodeCoords_Shared
#endif /*USE_LOADBALANCE*/
#endif /*PARTICLES*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars     ,ONLY: PerformLoadBalance,UseH5IOLoadBalance,offsetElemMPIOld
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: FileString
LOGICAL,INTENT(IN)           :: ReadNodes  !< calls ReadMeshElems() and ReadMeshNodes() if true (always true for PARTICLES=ON)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER            :: aElem
TYPE(tSide),POINTER            :: aSide,bSide
REAL,ALLOCATABLE               :: NodeCoordsTmp(:,:,:,:,:)
INTEGER                        :: BCindex
INTEGER                        :: iElem,ElemID
INTEGER                        :: iNode
INTEGER                        :: iLocSide,nbLocSide
INTEGER                        :: iSide
INTEGER                        :: FirstSideInd,LastSideInd,FirstElemInd,LastElemInd
INTEGER                        :: nPeriodicSides,nMPIPeriodics
INTEGER                        :: ReduceData(11)
INTEGER                        :: nSideIDs,offsetSideID
INTEGER                        :: iMortar,jMortar,nMortars
#if USE_MPI
INTEGER                        :: iNbProc
INTEGER                        :: iProc
INTEGER,ALLOCATABLE            :: MPISideCount(:)
#endif /*USE_MPI*/
LOGICAL                        :: doConnection
LOGICAL                        :: oriented
#ifdef PARTICLES
REAL, ALLOCATABLE              :: GlobVarTimeStep(:)
#endif
REAL                           :: StartT,EndT
INTEGER                        :: NGeoOld
!===================================================================================================================================
IF(MESHInitIsDone) RETURN

#if defined(PARTICLES) && USE_LOADBALANCE
IF (.NOT.PerformLoadBalance) THEN
#endif /*defined(PARTICLES) && USE_LOADBALANCE*/
  IF(MPIRoot) THEN
    IF(.NOT.FILEEXISTS(FileString))  CALL Abort(__STAMP__,'readMesh from data file "'//TRIM(FileString)//'" does not exist')
  END IF
  SWRITE(UNIT_stdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A)',ADVANCE="NO")' READ MESH FROM DATA FILE "'//TRIM(FileString)//'" ...'
  GETTIME(StartT)

  ! Get ElemInfo from Mesh file
  CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)
  CALL GetDataSize(File_ID,'ElemInfo',nDims,HSize)
  CALL ReadAttribute(File_ID,'nUniqueSides',1,IntScalar=nGlobalUniqueSidesFromMesh)
  CALL ReadAttribute(File_ID,'nSides',1,IntScalar=nNonUniqueGlobalSides)
  CALL ReadAttribute(File_ID,'nNodes',1,IntScalar=nNonUniqueGlobalNodes)
  CALL CloseDataFile()
  ! INTEGER KIND=4 check for number of elements
  CHECKSAFEINT(HSize(2),4)
  nGlobalElems=INT(HSize(2),4) !global number of elements
  ! INTEGER KIND=4 check for number of nodes
  CHECKSAFEINT(8_8*INT(nGlobalElems,8),4)
  DEALLOCATE(HSize)
  IF(MPIRoot.AND.(nGlobalElems.LT.nProcessors))CALL abort(__STAMP__&
      ,' Number of elements < number of processors',nGlobalElems,REAL(nProcessors))
  GETTIME(EndT)
  ReadMeshWallTime=EndT-StartT
  CALL DisplayMessageAndTime(ReadMeshWallTime, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
#if defined(PARTICLES) && USE_LOADBALANCE
END IF
#endif /*defined(PARTICLES) && USE_LOADBALANCE*/

#if USE_LOADBALANCE
IF (PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)) THEN
  SDEALLOCATE(offsetElemMPIOld)
  ALLOCATE(   offsetElemMPIOld(0:nProcessors))
  offsetElemMPIOld = offsetElemMPI
END IF
#endif /*USE_LOADBALANCE*/

#if USE_MPI
SDEALLOCATE(offsetElemMPI)
ALLOCATE(   offsetElemMPI(0:nProcessors))
offsetElemMPI = 0
#endif /*USE_MPI*/

!----------------------------------------------------------------------------------------------------------------------------
!                              DOMAIN DECOMPOSITION
!----------------------------------------------------------------------------------------------------------------------------
#if USE_MPI
CALL DomainDecomposition()
#else /*USE_MPI*/
nElems=nGlobalElems   ! Local number of Elements
offsetElem=0          ! Offset is the index of first entry, hdf5 array starts at 0-.GT. -1
#endif /*USE_MPI*/

!----------------------------------------------------------------------------------------------------------------------------
!                              ALLOCATE element counters
!----------------------------------------------------------------------------------------------------------------------------
#ifdef PARTICLES
! Re-allocate nPartsPerElem depending on new number of elements
IF(.NOT.ALLOCATED(nPartsPerElem))THEN
  ALLOCATE(nPartsPerElem(1:nElems))
END IF
nPartsPerElem=0
CALL AddToElemData(ElementOut,'nPartsPerElem',LongIntArray=nPartsPerElem(:))
#if USE_LOADBALANCE
ALLOCATE(nDeposPerElem(1:nElems))
nDeposPerElem=0
ALLOCATE(nTracksPerElem(1:nElems))
nTracksPerElem=0
ALLOCATE(nSurfacefluxPerElem(1:nElems))
nSurfacefluxPerElem=0
ALLOCATE(nPartsPerBCElem(1:nElems))
nPartsPerBCElem=0
ALLOCATE(nSurfacePartsPerElem(1:nElems))
nSurfacePartsPerElem=0
#endif /*USE_LOADBALANCE*/
#endif /*PARTICLES*/

!----------------------------------------------------------------------------------------------------------------------------
!                              VARIABLE TIME STEP
!----------------------------------------------------------------------------------------------------------------------------
#ifdef PARTICLES
IF(VarTimeStep%UseDistribution) THEN
  IF(ALLOCATED(VarTimeStep%ElemFac)) THEN
    ALLOCATE(GlobVarTimeStep(nGlobalElems))
    GlobVarTimeStep(1:nGlobalelems) = VarTimeStep%ElemFac(1:nGlobalelems)
    ! Allocate new array for the time step distribution (going from global time step distribution to proc local with nElems)
    DEALLOCATE(VarTimeStep%ElemFac)
    ALLOCATE(VarTimeStep%ElemFac(nElems))
    ! And now get the new variable time steps for the respective elements
    VarTimeStep%ElemFac(1:nElems) = GlobVarTimeStep(offsetElem+1:offsetElem+nElems)
    ! Global distribution is not required anymore
    DEALLOCATE(GlobVarTimeStep)
  ELSE
    ! Allocate the array for the element-wise time step factor
    ALLOCATE(VarTimeStep%ElemFac(nElems))
    VarTimeStep%ElemFac = 1.0
#if USE_MPI
    ! Allocate the array for the element-wise weighting factor
    ALLOCATE(VarTimeStep%ElemWeight(nElems))
    VarTimeStep%ElemWeight = 1.0
#endif
  END IF
END IF
#endif

#if defined(PARTICLES) && USE_LOADBALANCE
IF (.NOT.PerformLoadBalance) THEN
#endif /*defined(PARTICLES) && USE_LOADBALANCE*/
  CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)
  CALL ReadBCs()
#if defined(PARTICLES) && USE_LOADBALANCE
END IF
#endif /*defined(PARTICLES) && USE_LOADBALANCE*/
!----------------------------------------------------------------------------------------------------------------------------
!                              ELEMENTS
!----------------------------------------------------------------------------------------------------------------------------

!read local ElemInfo from data file
FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+nElems
ALLOCATE(Elems(                   FirstElemInd:LastElemInd))
ALLOCATE(ElemInfo(ELEMINFOSIZE_H5,FirstElemInd:LastElemInd))

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      ElemInfoSize => INT(ELEMINFOSIZE_H5,IK) ,&
      nElems       => INT(nElems,IK)       ,&
      offsetElem   => INT(offsetElem,IK)  )
#if defined(PARTICLES) && USE_LOADBALANCE
  IF (PerformLoadBalance) THEN
    ElemInfo(1:ElemInfoSize,:) = ElemInfo_Shared(1:ElemInfoSize,offsetElem+1:offsetElem+nElems)
  ELSE
#endif /*defined(PARTICLES) && USE_LOADBALANCE*/
    CALL ReadArray('ElemInfo',2,(/ElemInfoSize,nElems/),offsetElem,2,IntegerArray_i4=ElemInfo(1:ElemInfoSize,:))
#if defined(PARTICLES) && USE_LOADBALANCE
  END IF
#endif /*defined(PARTICLES) && USE_LOADBALANCE*/
END ASSOCIATE

DO iElem=FirstElemInd,LastElemInd
  iSide=ElemInfo(ELEM_FIRSTSIDEIND,iElem) !first index -1 in Sideinfo
  iNode=ElemInfo(ELEM_FIRSTNODEIND,iElem) !first index -1 in NodeInfo
  Elems(iElem)%ep=>GETNEWELEM()
  aElem=>Elems(iElem)%ep
  aElem%Ind    = iElem
  aElem%Type   = ElemInfo(ELEM_TYPE,iElem)
  ! Sanity check: Allow only specific element types
  SELECT CASE(aElem%Type)
  CASE(108,118,208)
    ! linear hex (108), non-linear hex (118), spline hex (208)
  CASE DEFAULT
    ! Abort if non-hexahedral meshes are read
    CALL abort(__STAMP__,'aElem%Type is NOT allowed: ',IntInfoOpt=aElem%Type)
  END SELECT
  aElem%Zone   = ElemInfo(ELEM_ZONE,iElem)
END DO

! Get number of compute-node elements (required for simulations with PARTICLES=ON/OFF)
#if USE_MPI
nComputeNodeElems = offsetElemMPI(ComputeNodeRootRank+nComputeNodeProcessors) - offsetElemMPI(ComputeNodeRootRank)
#else
nComputeNodeElems = nElems
#endif /*USE_MPI*/

!#ifdef PARTICLES
! Get ElemInfo_Shared(1:ELEMINFOSIZE,1:nGlobalElems)
IF(ReadNodes) CALL ReadMeshElems()
!#endif

!----------------------------------------------------------------------------------------------------------------------------
!                              SIDES
!----------------------------------------------------------------------------------------------------------------------------

! get offset of local side indices in the global sides
offsetSideID = ElemInfo(ELEM_FIRSTSIDEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nSideIDs     = ElemInfo(ELEM_LASTSIDEIND,LastElemInd)-ElemInfo(ELEM_FIRSTSIDEIND,FirstElemInd)
!read local SideInfo from data file
FirstSideInd = offsetSideID+1
LastSideInd  = offsetSideID+nSideIDs
ALLOCATE(SideInfo(SIDEINFOSIZE,FirstSideInd:LastSideInd))

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      SideInfoSize   => INT(SIDEINFOSIZE_H5,IK)   ,&
      nSideIDs       => INT(nSideIDs,IK)       ,&
      offsetSideID   => INT(offsetSideID,IK)  )
#if defined(PARTICLES) && USE_LOADBALANCE
  IF (PerformLoadBalance) THEN
    SideInfo(1:SideInfoSize,:) = SideInfo_Shared(1:SideInfoSize,offsetSideID+1:offsetSideID+nSideIDs)
  ELSE
#endif /*defined(PARTICLES) && USE_LOADBALANCE*/
    CALL ReadArray('SideInfo',2,(/SideInfoSize,nSideIDs/),offsetSideID,2,IntegerArray_i4=SideInfo(1:SideInfoSize,:))
    CALL CloseDataFile()
#if defined(PARTICLES) && USE_LOADBALANCE
  END IF
#endif /*defined(PARTICLES) && USE_LOADBALANCE*/
END ASSOCIATE

#ifdef PARTICLES
! Get SideInfo_Shared(1:SIDEINFOSIZE+1,1:nNonUniqueGlobalSides)
CALL ReadMeshSides()
#endif

! iterate over all local elements and within each element over all sides
DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  iSide=ElemInfo(ELEM_FIRSTSIDEIND,iElem) !first index -1 in Sideinfo
  !build up sides of the element according to CGNS standard
  ! assign flip
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    iSide=iSide+1
    ! ALLOCATE MORTAR
    ElemID=SideInfo(SIDE_NBELEMID,iSide) !IF nbElemID <0, this marks a mortar master side.
                                         ! The number (-1,-2,-3) is the Type of mortar
    IF(ElemID.LT.0)THEN ! mortar Sides attached!
      aSide%MortarType=ABS(ElemID)
      SELECT CASE(aSide%MortarType)
      CASE(1)
        aSide%nMortars=4
      CASE(2,3)
        aSide%nMortars=2
      END SELECT
      ALLOCATE(aSide%MortarSide(aSide%nMortars))
      DO iMortar=1,aSide%nMortars
        aSide%MortarSide(iMortar)%sp=>GETNEWSIDE()  !mortarType=0
      END DO
    ELSE
      aSide%nMortars=0
    END IF
    IF(SideInfo(SIDE_TYPE,iSide).LT.0) aSide%MortarType=-10 !marks small neighbor  side as belonging to a mortar

    IF(aSide%MortarType.LE.0)THEN
      aSide%Elem=>aElem
      oriented=(SideInfo(SIDE_ID,iSide).GT.0)
      aSide%Ind=ABS(SideInfo(SIDE_ID,iSide))
      IF(oriented)THEN !oriented side
        aSide%flip=0
#ifdef PARTICLES
        aSide%BC_Alpha=99
#endif /*PARTICLES*/
      ELSE !not oriented
        aSide%flip=MOD(Sideinfo(SIDE_FLIP,ISIde),10)
        IF((aSide%flip.LT.0).OR.(aSide%flip.GT.4)) STOP 'NodeID doesnt belong to side'
#ifdef PARTICLES
        aSide%BC_Alpha=-99
#endif /*PARTICLES*/
      END IF
    ELSE !mortartype>0
#if defined(PARTICLES) || USE_HDG
      ! Store global unique side index
      aSide%Ind=-ABS(SideInfo(SIDE_ID,iSide))
#endif /*defined(PARTICLES) || USE_HDG*/
      DO iMortar=1,aSide%nMortars
        iSide=iSide+1
        aSide%mortarSide(iMortar)%sp%Elem=>aElem
        IF(SideInfo(SIDE_ID,iSide).LT.0) STOP 'Problem in Mortar readin,should be flip=0'
        aSide%mortarSide(iMortar)%sp%flip=0
        aSide%mortarSide(iMortar)%sp%Ind=ABS(SideInfo(SIDE_ID,iSide))
#ifdef PARTICLES
        aSide%BC_Alpha=99
#endif /*PARTICLES*/
      END DO !iMortar
    END IF
  END DO !i=1,locnSides
END DO !iElem

! build up side connection
DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  iSide=ElemInfo(ELEM_FIRSTSIDEIND,iElem) !first index -1 in Sideinfo
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    iSide=iSide+1
    ! LOOP over mortars, if no mortar, then LOOP is executed once
    nMortars=aSide%nMortars
    DO iMortar=0,nMortars
      IF(iMortar.GT.0)THEN
        iSide=iSide+1
        aSide=>aElem%Side(iLocSide)%sp%mortarSide(iMortar)%sp
      END IF
      elemID  = SideInfo(SIDE_NBELEMID,iSide)
      BCindex = SideInfo(SIDE_BCID,iSide)

      doConnection=.TRUE. ! for periodic sides if BC is reassigned as non periodic
      IF(BCindex.NE.0)THEN !BC
        aSide%BCindex = BCindex
        IF((BoundaryType(aSide%BCindex,BC_TYPE).NE.1).AND.&
           (BoundaryType(aSide%BCindex,BC_TYPE).NE.100))THEN ! Reassignment from periodic to non-periodic
          doConnection=.FALSE.
          aSide%flip  =0
#ifdef PARTICLES
          aSide%BC_Alpha=0
#endif /*PARTICLES*/
          IF(iMortar.EQ.0) aSide%mortarType  = 0
          IF(iMortar.EQ.0) aSide%nMortars    = 0
          elemID            = 0
        END IF
      ELSE
        aSide%BCindex = 0
      END IF

      !no connection for mortar master
      IF(aSide%mortarType.GT.0) CYCLE
      IF(.NOT.doConnection) CYCLE
      IF(ASSOCIATED(aSide%connection)) CYCLE

      ! check if neighbor on local proc or MPI connection
      IF(elemID.NE.0)THEN !connection
        IF((elemID.LE.LastElemInd).AND.(elemID.GE.FirstElemInd))THEN !local
          !TODO: Check if this is still ok
          DO nbLocSide=1,6
            bSide=>Elems(elemID)%ep%Side(nbLocSide)%sp
            ! LOOP over mortars, if no mortar, then LOOP is executed once
            nMortars=bSide%nMortars
            DO jMortar=0,nMortars
              IF(jMortar.GT.0) bSide=>Elems(elemID)%ep%Side(nbLocSide)%sp%mortarSide(jMortar)%sp

              IF(bSide%ind.EQ.aSide%ind)THEN
                aSide%connection=>bSide
                bSide%connection=>aSide
                EXIT
              END IF !bSide%ind.EQ.aSide%ind
            END DO !jMortar
          END DO !nbLocSide
        ELSE !MPI connection
#if USE_MPI
          aSide%connection=>GETNEWSIDE()
          aSide%connection%flip=aSide%flip
          aSide%connection%Elem=>GETNEWELEM()
          aSide%NbProc = ELEMIPROC(elemID)
#else
          CALL abort(__STAMP__, ' ElemID of neighbor not in global Elem list ')
#endif
        END IF
      END IF
#ifdef PARTICLES
      CALL ReadMeshSideNeighbors(elemID,iSide)
#endif
    END DO !iMortar
  END DO !iLocSide
END DO !iElem


!----------------------------------------------------------------------------------------------------------------------------
!                              NODES
!----------------------------------------------------------------------------------------------------------------------------
!#ifdef PARTICLES
! Particles want to node coordinates in the old 2D format, hence this read-in happens twice
! Get NodeInfo_Shared(1:8*nGlobalElems)
! Get NodeCoords_Shared(1:3,1:8*nGlobalElems)
IF(ReadNodes) CALL ReadMeshNodes()
!#endif

! get physical coordinates
#if defined(PARTICLES) && USE_LOADBALANCE
IF (.NOT.PerformLoadBalance) &
#endif /*defined(PARTICLES) && USE_LOADBALANCE*/
  CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)

! Backup required if useCurveds=F
NGeoOld = NGeo
! Linear mesh: set polynomial degree of geometry to 1 and rebuild NodeMap
IF(.NOT.useCurveds) NGeo = 1

ALLOCATE(NodeCoords(3,0:NGeo,0:NGeo,0:NGeo,nElems))

#if defined(PARTICLES) && USE_LOADBALANCE
IF (PerformLoadBalance) THEN
  NodeCoords = RESHAPE(NodeCoords_Shared(1:3,(NGeo+1)**3*offsetElem+1:(NGeo+1)**3*(offsetElem+nElems)),(/3,NGeo+1,NGeo+1,NGeo+1,nElems/))
ELSE
#endif /*defined(PARTICLES) && USE_LOADBALANCE*/
  IF(useCurveds)THEN
    CALL ReadArray('NodeCoords',2,(/3_IK,INT(nElems*(NGeo+1)**3,IK)/),INT(offsetElem*(NGeo+1)**3,IK),2,RealArray=NodeCoords)
  ELSE
    ALLOCATE(NodeCoordsTmp(3,0:NGeoOld,0:NGeoOld,0:NGeoOld,nElems))
    CALL ReadArray('NodeCoords',2,(/3_IK,INT(nElems*(NGeoOld+1)**3,IK)/),INT(offsetElem*(NGeoOld+1)**3,IK),2,RealArray=NodeCoordsTmp)
    ! throw away all nodes except the 8 corner nodes of each hexa
    NodeCoords(:,0,0,0,:) = NodeCoordsTmp(: , 0       , 0       , 0       , :)
    NodeCoords(:,1,0,0,:) = NodeCoordsTmp(: , NGeoOld , 0       , 0       , :)
    NodeCoords(:,0,1,0,:) = NodeCoordsTmp(: , 0       , NGeoOld , 0       , :)
    NodeCoords(:,1,1,0,:) = NodeCoordsTmp(: , NGeoOld , NGeoOld , 0       , :)
    NodeCoords(:,0,0,1,:) = NodeCoordsTmp(: , 0       , 0       , NGeoOld , :)
    NodeCoords(:,1,0,1,:) = NodeCoordsTmp(: , NGeoOld , 0       , NGeoOld , :)
    NodeCoords(:,0,1,1,:) = NodeCoordsTmp(: , 0       , NGeoOld , NGeoOld , :)
    NodeCoords(:,1,1,1,:) = NodeCoordsTmp(: , NGeoOld , NGeoOld , NGeoOld , :)
    DEALLOCATE(NodeCoordsTmp)
  END IF ! useCurveds
#if defined(PARTICLES) && USE_LOADBALANCE
END IF
#endif /*defined(PARTICLES) && USE_LOADBALANCE*/

#if defined(PARTICLES) && USE_LOADBALANCE
IF (.NOT.PerformLoadBalance) &
#endif /*defined(PARTICLES) && USE_LOADBALANCE*/
  CALL CloseDataFile()

#ifdef PARTICLES
! Start non-blocking communication of mesh information
CALL StartCommunicateMeshReadin()
#endif

! Readin of mesh is now finished

!----------------------------------------------------------------------------------------------------------------------------
!                              COUNT SIDES
!----------------------------------------------------------------------------------------------------------------------------
! Readin is now finished
nBCSides=0
nAnalyzeSides=0
nMortarSides=0
nSides=0
nPeriodicSides=0
nMPIPeriodics=0
nMPISides=0
#if USE_MPI
ALLOCATE(MPISideCount(0:nProcessors-1))
MPISideCount=0
#endif
DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    ! LOOP over mortars, if no mortar, then LOOP is executed once
    nMortars=aSide%nMortars
    DO iMortar=0,nMortars
      IF(iMortar.GT.0) aSide=>aElem%Side(iLocSide)%sp%mortarSide(iMortar)%sp
      aSide%tmp=0
    END DO !iMortar
  END DO !iLocSide
END DO !iElem
DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    nMortars=aSide%nMortars
    DO iMortar=0,nMortars
      IF(iMortar.GT.0) aSide=>aElem%Side(iLocSide)%sp%mortarSide(iMortar)%sp

      IF(aSide%tmp.EQ.0)THEN
        nSides=nSides+1
        aSide%tmp=-1 !used as marker
        IF(ASSOCIATED(aSide%connection)) aSide%connection%tmp=-1
        IF(aSide%BCindex.NE.0)THEN !side is BC or periodic side
          nAnalyzeSides=nAnalyzeSides+1
          IF(ASSOCIATED(aSide%connection))THEN
            IF(BoundaryType(aSide%BCindex,BC_TYPE).EQ.1)THEN
              nPeriodicSides=nPeriodicSides+1
#if USE_MPI
              IF(aSide%NbProc.NE.-1) nMPIPeriodics=nMPIPeriodics+1
#endif
            END IF
          ELSE
            IF(aSide%MortarType.EQ.0)THEN !really a BC side
              nBCSides=nBCSides+1
            END IF
          END IF
        END IF
        IF(aSide%MortarType.GT.0) nMortarSides=nMortarSides+1
#if USE_MPI
        IF(aSide%NbProc.NE.-1) THEN
          nMPISides=nMPISides+1
          MPISideCount(aSide%NbProc)=MPISideCount(aSide%NbProc)+1
        END IF
#endif
      END IF
    END DO !iMortar
  END DO !iLocSide
END DO !iElem
nInnerSides=nSides-nBCSides-nMPISides-nMortarSides !periodic side count to inner side!!!

LOGWRITE(*,*)'-------------------------------------------------------'
LOGWRITE(*,'(A22,I8)')'nSides:',nSides
LOGWRITE(*,'(A22,I8)')'nBCSides:',nBCSides
LOGWRITE(*,'(A22,I8)')'nMortarSides:',nMortarSides
LOGWRITE(*,'(A22,I8)')'nInnerSides:',nInnerSides
LOGWRITE(*,'(A22,I8)')'nMPISides:',nMPISides
LOGWRITE(*,*)'-------------------------------------------------------'
 !now MPI sides
#if USE_MPI
nNBProcs=0
DO iProc=0,nProcessors-1
  IF(iProc.EQ.myRank) CYCLE
  IF(MPISideCount(iProc).GT.0) nNBProcs=nNbProcs+1
END DO
IF(nNbProcs.EQ.0)THEN !MPI + 1Proc case !
  ALLOCATE(NbProc(1),nMPISides_Proc(1))
  nNbProcs=1
  NbProc=0
  nMPISides_Proc=0
ELSE
  ALLOCATE(NbProc(nNbProcs),nMPISides_Proc(1:nNbProcs))
  iNbProc=0
  DO iProc=0,nProcessors-1
    IF(iProc.EQ.myRank) CYCLE
    IF(MPISideCount(iProc).GT.0) THEN
      iNbProc=iNbProc+1
      NbProc(iNbProc)=iProc
      ! compute number of MPISides per neighbor proc and divide by two
      nMPISides_Proc(iNBProc)=MPISideCount(iProc)
    END IF
  END DO
END IF
DEALLOCATE(MPISideCount)
#endif /*USE_MPI*/

ReduceData(1)=nElems
ReduceData(2)=nSides
ReduceData(3)=nNodes
ReduceData(4)=nInnerSides
ReduceData(5)=nPeriodicSides
ReduceData(6)=nBCSides
ReduceData(7)=nMPISides
ReduceData(8)=nAnalyzeSides
ReduceData(9)=nMortarSides
ReduceData(10)=nMPIPeriodics
ReduceData(11)=nNodeIDs

#ifdef PARTICLES
! Finish non-blocking communication of mesh information
CALL FinishCommunicateMeshReadin()
#endif

#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,ReduceData,11,MPI_INTEGER,MPI_SUM,MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/

nGlobalMortarSides=ReduceData(9)

IF(MPIRoot)THEN
#if USE_LOADBALANCE
  IF(.NOT.PerformLoadBalance)THEN
#endif /*USE_LOADBALANCE*/
    WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nElems | ',ReduceData(1) !nElems
    WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nNodes, unique | ',ReduceData(3) !nNodes
    WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nNodes, total  | ',ReduceData(11) !nNodes
    WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nSides         | ',ReduceData(2)-ReduceData(7)/2
    WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nSides,    BC  | ',ReduceData(6) !nBCSides
    WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nSides,   MPI  | ',ReduceData(7)/2 !nMPISides
    WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nSides, Inner  | ',ReduceData(4) !nInnerSides
    WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nSides,Mortar  | ',nGlobalMortarSides
    WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nPeriodicSides,Total | ',ReduceData(5)-ReduceData(10)/2
    WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nPeriodicSides,Inner | ',ReduceData(5)-ReduceData(10)
    WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nPeriodicSides,  MPI | ',ReduceData(10)/2 !nPeriodicSides
    WRITE(UNIT_stdOut,'(A,A34,I0)')' |','nAnalyzeSides | ',ReduceData(8) !nAnalyzeSides
    WRITE(UNIT_stdOut,'(A,A34,L1)')' |','useCurveds | ',useCurveds
    WRITE(UNIT_stdOut,'(A,A34,I0)')' |','Ngeo | ',Ngeo
    WRITE(UNIT_stdOut,'(132("."))')
#if USE_LOADBALANCE
  END IF ! .NOT.PerformLoadBalance
#endif /*USE_LOADBALANCE*/
END IF

LOGWRITE_BARRIER

LBWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE ReadMesh


SUBROUTINE ReadMeshElems()
!===================================================================================================================================
!> Create particle mesh arrays for elems:
!> - ElemInfo_Shared(1:ELEMINFOSIZE,1:nGlobalElems)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
USE MOD_Particle_Mesh_Vars
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

#if USE_MPI
#if USE_LOADBALANCE
IF (.NOT.PerformLoadBalance) THEN
#endif /*USE_LOADBALANCE*/
  ! allocate shared array for ElemInfo
  CALL Allocate_Shared((/ELEMINFOSIZE,nGlobalElems/),ElemInfo_Shared_Win,ElemInfo_Shared)
  CALL MPI_WIN_LOCK_ALL(0,ElemInfo_Shared_Win,IERROR)

  ElemInfo_Shared(1:ELEMINFOSIZE_H5,offsetElem+1:offsetElem+nElems) = ElemInfo(:,:)
  ElemInfo_Shared(ELEM_RANK        ,offsetElem+1:offsetElem+nElems) = myRank
  CALL BARRIER_AND_SYNC(ElemInfo_Shared_Win,MPI_COMM_SHARED)
#if USE_LOADBALANCE
ELSEIF(UseH5IOLoadBalance)THEN
  ElemInfo_Shared(ELEM_RANK        ,offsetElem+1:offsetElem+nElems) = myRank
  CALL BARRIER_AND_SYNC(ElemInfo_Shared_Win,MPI_COMM_SHARED)
END IF
#endif /*USE_LOADBALANCE*/
#endif  /*USE_MPI*/

#if USE_MPI
! broadcast elem offset of compute-node root
offsetComputeNodeElem=offsetElem
CALL MPI_BCAST(offsetComputeNodeElem,1, MPI_INTEGER,0,MPI_COMM_SHARED,iERROR)
#else
! allocate local array for ElemInfo
ALLOCATE(ElemInfo_Shared(1:ELEMINFOSIZE,1:nElems))
ElemInfo_Shared(1:ELEMINFOSIZE_H5,1:nElems) = ElemInfo(:,:)
#endif  /*USE_MPI*/

END SUBROUTINE ReadMeshElems


#if defined(PARTICLES)
SUBROUTINE ReadMeshSides()
!===================================================================================================================================
!> Create particle mesh arrays for sides:
!> - SideInfo_Shared(1:SIDEINFOSIZE+1,1:nNonUniqueGlobalSides)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
USE MOD_Particle_Mesh_Vars
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars     ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: FirstElemInd,LastElemInd
INTEGER                        :: nSideIDs,offsetSideID
!===================================================================================================================================

FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+nElems
offsetSideID = ElemInfo(ELEM_FIRSTSIDEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nSideIDs     = ElemInfo(ELEM_LASTSIDEIND ,LastElemInd)-ElemInfo(ELEM_FIRSTSIDEIND,FirstElemInd)

ALLOCATE(SideInfo_Shared_tmp(offsetSideID+1:offsetSideID+nSideIDs))
SideInfo_Shared_tmp = 0

#if USE_MPI
! all procs on my compute-node communicate the number of non-unique sides
CALL MPI_ALLREDUCE(nSideIDs,nComputeNodeSides,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,IERROR)

#if USE_LOADBALANCE
IF (PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/

CALL Allocate_Shared((/SIDEINFOSIZE+1,nNonUniqueGlobalSides/),SideInfo_Shared_Win,SideInfo_Shared)
CALL MPI_WIN_LOCK_ALL(0,SideInfo_Shared_Win,IERROR)
SideInfo_Shared(1                :SIDEINFOSIZE  ,offsetSideID+1:offsetSideID+nSideIDs) = SideInfo(:,:)
SideInfo_Shared(SIDEINFOSIZE_H5+1:SIDEINFOSIZE+1,offsetSideID+1:offsetSideID+nSideIDs) = 0
CALL BARRIER_AND_SYNC(SideInfo_Shared_Win,MPI_COMM_SHARED)
#else
nComputeNodeSides = nSideIDs
ALLOCATE(SideInfo_Shared(1:SIDEINFOSIZE+1         , 1:nSideIDs))
SideInfo_Shared(1                :SIDEINFOSIZE    , 1:nSideIDs) = SideInfo(:,:)
SideInfo_Shared(SIDEINFOSIZE_H5+1:SIDEINFOSIZE+1  , 1:nSideIDs) = 0
#endif /*USE_MPI*/

END SUBROUTINE ReadMeshSides
#endif /*defined(PARTICLES)*/


SUBROUTINE ReadMeshNodes()
!===================================================================================================================================
!> Create particle mesh arrays for nodes:
!> - NodeInfo_Shared(1:8*nGlobalElems)
!> - NodeCoords_Shared(1:3,1:8*nGlobalElems)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_HDF5_Input         ,ONLY: ReadArray,OpenDataFile
USE MOD_IO_HDF5            ,ONLY: CloseDataFile
USE MOD_Mesh_Vars
USE MOD_Particle_Mesh_Vars
USE MOD_Mesh_Tools         ,ONLY: GetCornerNodes
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,iNode
INTEGER                        :: FirstElemInd,LastElemInd
INTEGER                        :: FirstNodeInd,LastNodeInd
INTEGER                        :: nNodeIDs,offsetNodeID
INTEGER,ALLOCATABLE            :: NodeInfo(:),NodeInfoTmp(:,:)
REAL,ALLOCATABLE               :: NodeCoords_indx(:,:)
INTEGER                        :: nNodeInfoIDs,NodeID,NodeCounter
INTEGER                        :: CNS(8)
!===================================================================================================================================

! calculate all offsets
FirstElemInd = offsetElem+1
LastElemInd  = offsetElem+nElems
offsetNodeID = ElemInfo_Shared(ELEM_FIRSTNODEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nNodeIDs     = ElemInfo_Shared(ELEM_LASTNODEIND ,LastElemInd)-ElemInfo_Shared(ELEM_FIRSTNODEIND,FirstElemind)

#if USE_LOADBALANCE
IF (PerformLoadBalance) THEN
  CALL MPI_ALLREDUCE(nNodeIDs,nComputeNodeNodes,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,IERROR)
  RETURN
END IF
#endif /*USE_LOADBALANCE*/

FirstNodeInd = offsetNodeID+1
LastNodeInd  = offsetNodeID+nNodeIDs

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nNodeIDs     => INT(nNodeIDs,IK)     ,&
      offsetNodeID => INT(offsetNodeID,IK) )
  ALLOCATE(NodeCoords_indx(3,nNodeIDs))
  ! read all nodes
  CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)
  CALL ReadArray('NodeCoords',2,(/3_IK,nNodeIDs/),offsetNodeID,2,RealArray=NodeCoords_indx)
  CALL CloseDataFile()
END ASSOCIATE

! Keep all nodes if elements are curved
IF (useCurveds.OR.NGeo.EQ.1) THEN
  MeshWasCurved = .TRUE.
  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (&
        nNodeIDs     => INT(nNodeIDs,IK)     ,&
        offsetNodeID => INT(offsetNodeID,IK) )
    ALLOCATE(NodeInfo(FirstNodeInd:LastNodeInd))
    CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)
    CALL ReadArray('GlobalNodeIDs',1,(/nNodeIDs/),offsetNodeID,1,IntegerArray_i4=NodeInfo)
    CALL CloseDataFile()
  END ASSOCIATE

#if USE_MPI
  ! allocate shared array for NodeInfo
  CALL MPI_ALLREDUCE(nNodeIDs,nComputeNodeNodes,1,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,IERROR)
  CALL Allocate_Shared((/nNonUniqueGlobalNodes/),NodeInfo_Shared_Win,NodeInfo_Shared)
  CALL MPI_WIN_LOCK_ALL(0,NodeInfo_Shared_Win,IERROR)
  NodeInfo_Shared(offsetNodeID+1:offsetNodeID+nNodeIDs) = NodeInfo(:)
  CALL BARRIER_AND_SYNC(NodeInfo_Shared_Win,MPI_COMM_SHARED)

  CALL Allocate_Shared((/3,nNonUniqueGlobalNodes/),NodeCoords_Shared_Win,NodeCoords_Shared)
  CALL MPI_WIN_LOCK_ALL(0,NodeCoords_Shared_Win,IERROR)
  NodeCoords_Shared(:,offsetNodeID+1:offsetNodeID+nNodeIDs) = NodeCoords_indx(:,:)
#else
  nComputeNodeNodes = nNodeIDs
  ALLOCATE(NodeInfo_Shared(1:nNodeIDs))
  NodeInfo_Shared(1:nNodeIDs) = NodeInfo(:)
  ALLOCATE(NodeCoords_Shared(3,nNodeIDs))
  NodeCoords_Shared(:,:) = NodeCoords_indx(:,:)
#endif  /*USE_MPI*/

! Reduce NodeCoords if no curved elements are to be used
ELSE ! .NOT. (useCurveds.OR.NGeo.EQ.1)
  ! every proc needs to allocate the array
  ALLOCATE(NodeInfo(1:nNonUniqueGlobalNodes))

#if USE_MPI
  ! root reads NodeInfo for new mapping
  IF (myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
    ! Associate construct for integer KIND=8 possibility
    ASSOCIATE (nNonUniqueGlobalNodes     => INT(nNonUniqueGlobalNodes,IK))
      CALL OpenDataFile(MeshFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
      CALL ReadArray('GlobalNodeIDs',1,(/nNonUniqueGlobalNodes/),0_IK,1,IntegerArray_i4=NodeInfo)
      CALL CloseDataFile()
    END ASSOCIATE
#if USE_MPI
  END IF

  ! root broadcasts NodeInfo to all procs on compute node
  CALL MPI_BCAST(NodeInfo,nNonUniqueGlobalNodes,MPI_INTEGER,0,MPI_COMM_SHARED,IERROR)
#endif /*USE_MPI*/

  ! Every proc builds new mapping. This step is required for consistency reasons since ElemInfo is not yet communicated
  nNodeInfoIDs = MAXVAL(NodeInfo)
  ALLOCATE(NodeInfoTmp(2,nNodeInfoIDs))
  NodeInfoTmp = 0

  ! Flag unique node IDs we will keep
  DO iNode = 1,nNonUniqueGlobalNodes
    NodeID = NodeInfo(iNode)
    NodeInfoTmp(1,NodeID) = 1
  END DO

  ! Build new NodeInfo IDs
  NodeCounter = 0
  DO iNode = 1,nNodeInfoIDs
    IF (NodeInfoTmp(1,iNode).EQ.0) CYCLE

    NodeCounter = NodeCounter + 1
    NodeInfoTmp(2,iNode) = NodeCounter
  END DO

#if USE_MPI
  CALL Allocate_Shared((/8*nGlobalElems/),NodeInfo_Shared_Win,NodeInfo_Shared)
  CALL MPI_WIN_LOCK_ALL(0,NodeInfo_Shared_Win,IERROR)
#else
  ALLOCATE(NodeInfo_Shared(8*nGlobalElems))
#endif /*USE_MPI*/

  ! Only the 8 corner nodes count for nodes. (NGeo+1)**2 = 8
  nComputeNodeNodes = 8*nComputeNodeElems

#if USE_MPI
  CALL Allocate_Shared((/3,8*nGlobalElems/),NodeCoords_Shared_Win,NodeCoords_Shared)
  CALL MPI_WIN_LOCK_ALL(0,NodeCoords_Shared_Win,IERROR)
#else
  ALLOCATE(NodeCoords_Shared(3,8*nGlobalElems))
#endif  /*USE_MPI*/

  ! the cornernodes are not the first 8 entries (for Ngeo>1) of nodeinfo array so mapping is built
  CNS(1:8) = GetCornerNodes(NGeo)

  ! throw away all nodes except the 8 corner nodes of each hexa
  nNonUniqueGlobalNodes = 8*nGlobalElems

  DO iElem = FirstElemInd,LastElemInd
    FirstNodeInd = ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem) - offsetNodeID
    ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem) = 8*(iElem-1)
    ElemInfo_Shared(ELEM_LASTNODEIND ,iElem) = 8* iElem
    DO iNode = 1,8
      NodeCoords_Shared(:,8*(iElem-1) + iNode) = NodeCoords_indx(:,FirstNodeInd+CNS(iNode))
      NodeInfo_Shared  (  8*(iElem-1) + iNode) = NodeInfoTmp(2,NodeInfo(FirstNodeInd+offsetNodeID+CNS(iNode)))
    END DO
  END DO

  DEALLOCATE(NodeInfoTmp)

END IF ! useCurveds.OR.NGeo.EQ.1

! Update node counters
offsetNodeID = ElemInfo_Shared(ELEM_FIRSTNODEIND,FirstElemInd) ! hdf5 array starts at 0-> -1
nNodeIDs     = ElemInfo_Shared(ELEM_LASTNODEIND ,LastElemInd)-ElemInfo_Shared(ELEM_FIRSTNODEIND,FirstElemind)
FirstNodeInd = offsetNodeID+1
LastNodeInd  = offsetNodeID+nNodeIDs

! Sanity check
IF(ABS(meshScale).LE.0.) CALL abort(__STAMP__,'meshScale is zero')
! scale mesh if desired. Mesh deformation currently not supported!
IF (ABS(meshScale-1.).GT.1e-14) THEN
  NodeCoords_Shared(:,FirstNodeInd:LastNodeInd) = NodeCoords_Shared(:,FirstNodeInd:LastNodeInd) * meshScale
END IF

#if USE_MPI
CALL BARRIER_AND_SYNC(ElemInfo_Shared_Win  ,MPI_COMM_SHARED) ! Only changed here, created in ReadMeshElems()
CALL BARRIER_AND_SYNC(NodeCoords_Shared_Win,MPI_COMM_SHARED) ! Created here
CALL BARRIER_AND_SYNC(NodeInfo_Shared_Win  ,MPI_COMM_SHARED) ! Created here
#endif  /*USE_MPI*/

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
  DEALLOCATE(NodeInfo)
#if USE_MPI
END IF
#endif /*USE_MPI*/
DEALLOCATE(NodeCoords_indx)

END SUBROUTINE ReadMeshNodes


FUNCTION INVMAP(ID,nIDs,ArrID)
!===================================================================================================================================
! find the inverse Mapping p.e. NodeID-> entry in NodeMap (a sorted array of unique NodeIDs), using bisection
! if Index is not in the range, -1 will be returned, if it is in the range, but is not found, 0 will be returned!!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                :: ID            ! ID to search for
INTEGER, INTENT(IN)                :: nIDs          ! size of ArrID
INTEGER, INTENT(IN)                :: ArrID(nIDs)   ! 1D array of IDs
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                            :: INVMAP               ! index of ID in NodeMap array
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: i,maxSteps,low,up,mid
!===================================================================================================================================
INVMAP=0
maxSteps=INT(LOG(REAL(nIDs))*1.4426950408889634556)+1    !1/LOG(2.)=1.4426950408889634556
low=1
up=nIDs
IF((ID.LT.ArrID(low)).OR.(ID.GT.ArrID(up))) THEN
  !WRITE(*,*)'WARNING, Node Index Not in local range -> set to -1'
  INVMAP=-1  ! not in the range!
  RETURN
END IF
IF(ID.EQ.ArrID(low))THEN
  INVMAP=low
ELSEIF(ID.EQ.ArrID(up))THEN
  INVMAP=up
ELSE
  !bisection
  DO i=1,maxSteps
    mid=(up-low)/2+low
    IF(ID .EQ. ArrID(mid))THEN
      INVMAP=mid                     !index found!
      EXIT
    ELSEIF(ID .GT. ArrID(mid))THEN ! seek in upper half
      low=mid
    ELSE
      up=mid
    END IF
  END DO
END IF
END FUNCTION INVMAP


#if USE_MPI
PPURE FUNCTION ELEMIPROC(ElemID)
!===================================================================================================================================
!> Find the id of a processor on which an element with a given ElemID lies, based on the MPI element offsets defined earlier.
!> Use a bisection algorithm for faster search.
!===================================================================================================================================
! MODULES
USE MOD_Globals,   ONLY:nProcessors
USE MOD_MPI_vars,  ONLY:offsetElemMPI
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                :: ElemID     !< (IN)  NodeID to search for
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                            :: ELEMIPROC  !< (OUT) processor id
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: i,maxSteps,low,up,mid
!===================================================================================================================================
ELEMIPROC=0
maxSteps=INT(LOG(REAL(nProcessors))*1.4426950408889634556)+1    !1/LOG(2.)=1.4426950408889634556
low=0
up=nProcessors-1
IF((ElemID.GT.offsetElemMPI(low)).AND.(ElemID.LE.offsetElemMPI(low+1)))THEN
  ELEMIPROC=low
ELSEIF((ElemID.GT.offsetElemMPI(up)).AND.(ElemID.LE.offsetElemMPI(up+1)))THEN
  ELEMIPROC=up
ELSE
  !bisection
  DO i=1,maxSteps
    mid=(up-low)/2+low
    IF((ElemID.GT.offsetElemMPI(mid)).AND.(ElemID.LE.offsetElemMPI(mid+1)))THEN
      ELEMIPROC=mid                     !index found!
      EXIT
    ELSEIF(ElemID .GT. offsetElemMPI(mid+1))THEN ! seek in upper half
      low=mid+1
    ELSE
      up=mid
    END IF
  END DO
END IF
END FUNCTION ELEMIPROC
#endif /*USE_MPI*/

RECURSIVE SUBROUTINE Qsort1Int(A)
!===================================================================================================================================
! QuickSort for integer array A
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(INOUT)            :: A(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: marker
!===================================================================================================================================
IF(SIZE(A).GT.1) THEN
  CALL Partition1Int(A,marker)
  CALL Qsort1Int(A(:marker-1))
  CALL Qsort1Int(A(marker:))
END IF
RETURN
END SUBROUTINE Qsort1Int



SUBROUTINE Partition1Int(A,marker)
!===================================================================================================================================
! Neeeded by QuickSort
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(INOUT)            :: A(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)              :: marker
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: i,j
INTEGER                          :: temp,x
!===================================================================================================================================
x= A(1)
i= 0
j= SIZE(A)+1
DO
  j=j-1
  DO
    IF(A(j).LE.x) EXIT
    j=j-1
  END DO
  i=i+1
  DO
    IF(A(i).GE.x) EXIT
    i=i+1
  END DO
  IF(i.LT.j)THEN
    ! exchange A(i) and A(j)
    temp=A(i)
    A(i)=A(j)
    A(j)=temp
  ELSEIF(i.EQ.j)THEN
    marker=i+1
    RETURN
  ELSE
    marker=i
    RETURN
  ENDIF
END DO
RETURN
END SUBROUTINE Partition1Int


SUBROUTINE FinalizeMeshReadin(meshMode)
!===================================================================================================================================
! Finalizes the shared mesh readin
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
USE MOD_Particle_Mesh_Vars
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
USE MOD_Particle_Vars             ,ONLY: Symmetry
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars          ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: meshMode !<  0: only read and build Elem_xGP,
                               !< -1: as 0 + build connectivity and read node info (automatically read for PARTICLES=ON)
                               !<  1: as 0 + build connectivity
                               !<  2: as 1 + calc metrics
                               !<  3: as 2 but skip InitParticleMesh
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
#if USE_MPI && defined(PARTICLES)
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)

! symmetry sides and elem volumes/characteristic lengths
IF(ABS(meshMode).GT.1)THEN
  IF(Symmetry%Order.EQ.2) CALL UNLOCK_AND_FREE(SideIsSymSide_Shared_Win)
  CALL UNLOCK_AND_FREE(ElemVolume_Shared_Win)
  CALL UNLOCK_AND_FREE(ElemCharLength_Shared_Win)
END IF ! ABS(meshMode).GT.1
#endif /*USE_MPI && defined(PARTICLES)*/

! Then, free the pointers or arrays
ADEALLOCATE(SideIsSymSide_Shared)
ADEALLOCATE(ElemVolume_Shared)
ADEALLOCATE(ElemCharLength_Shared)

#if USE_MPI
! Free communication arrays
SDEALLOCATE(displsElem)
SDEALLOCATE(recvcountElem)
SDEALLOCATE(displsSide)
SDEALLOCATE(recvcountSide)

#if USE_LOADBALANCE
IF (PerformLoadBalance) THEN
  CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
  RETURN
END IF
#endif /*USE_LOADBALANCE*/

! elems
CALL UNLOCK_AND_FREE(ElemInfo_Shared_Win)

! sides
#if defined(PARTICLES)
CALL UNLOCK_AND_FREE(SideInfo_Shared_Win)
#endif /*defined(PARTICLES)*/

! nodes
CALL UNLOCK_AND_FREE(NodeInfo_Shared_Win)
CALL UNLOCK_AND_FREE(NodeCoords_Shared_Win)

CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
#endif /*USE_MPI*/

ADEALLOCATE(ElemInfo_Shared)
ADEALLOCATE(SideInfo_Shared)
ADEALLOCATE(NodeInfo_Shared)
ADEALLOCATE(NodeCoords_Shared)

! Free communication arrays
#if USE_MPI
SDEALLOCATE(displsNode)
SDEALLOCATE(recvcountNode)
#endif /*USE_MPI*/

END SUBROUTINE FinalizeMeshReadin



END MODULE MOD_Mesh_ReadIn
