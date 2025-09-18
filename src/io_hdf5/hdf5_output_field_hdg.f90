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

MODULE MOD_HDF5_Output_Fields_HDG
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
#if USE_HDG
#if defined(PARTICLES)
USE MOD_io_HDF5
USE MOD_HDF5_output
#endif /*defined(PARTICLES)*/
#endif /*USE_HDG*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
#if USE_HDG
#if defined(PARTICLES)
PUBLIC :: WriteBRAverageElemToHDF5
PUBLIC :: WriteSurfVDLToHDF5
#endif /*defined(PARTICLES)*/
#endif /*USE_HDG*/
!===================================================================================================================================

CONTAINS

#if USE_HDG
#if defined(PARTICLES)
SUBROUTINE WriteBRAverageElemToHDF5(isBRAverageElem)
!===================================================================================================================================
! write BRAverageElem field to HDF5 file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars        ,ONLY: MeshFile,nGlobalElems,offsetElem
USE MOD_io_HDF5
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: isBRAverageElem(1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER  :: N_variables=1
REAL               :: BRAverageElem(1:N_variables,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
CHARACTER(LEN=255) :: StrVarNames(1:N_variables)
CHARACTER(LEN=255) :: FileName
REAL               :: StartT,EndT
REAL               :: OutputTime
INTEGER            :: iElem
!===================================================================================================================================
#if USE_LOADBALANCE
IF(PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/
! create global zeta field for parallel output of zeta distribution
StrVarNames(1)='BRAverageElem'
BRAverageElem=0.
DO iElem=1,PP_nElems
  IF(isBRAverageElem(iElem))THEN
    BRAverageElem(:,:,:,:,iElem) = 1.0
  END IF
END DO!iElem
SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE BRAverageElem TO HDF5 FILE...'
GETTIME(StartT)
OutputTime=0.0
! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
CALL GenerateFileSkeleton('BRAverageElem',N_variables,StrVarNames,TRIM(MeshFile),OutputTime,FileNameOut=FileName)
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif
  CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)
CALL WriteAttributeToHDF5(File_ID,'VarNamesBRAverageElem',N_variables,StrArray=StrVarNames)
CALL CloseDataFile()

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
        nGlobalElems    => INT(nGlobalElems,IK)    ,&
        PP_nElems       => INT(PP_nElems,IK)       ,&
        N_variables     => INT(N_variables,IK)     ,&
        N               => INT(PP_N,IK)            ,&
        offsetElem      => INT(offsetElem,IK)      )
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
                          DataSetName='DG_Solution' , rank=5 , &
                          nValGlobal =(/N_variables , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
                          nVal       =(/N_variables , N+1_IK , N+1_IK , N+1_IK , PP_nElems   /) , &
                          offset     =(/       0_IK , 0_IK   , 0_IK   , 0_IK   , offsetElem  /) , &
                          collective =.TRUE.        , RealArray=BRAverageElem)
END ASSOCIATE
GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
END SUBROUTINE WriteBRAverageElemToHDF5


SUBROUTINE WriteSurfVDLToHDF5(OutputTime)
!===================================================================================================================================
!> write the final values of the surface sampling to a HDF5 state file
!> additional performs all the final required computations
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_Particle_Boundary_Vars ,ONLY: nComputeNodeSurfTotalSides
USE MOD_Particle_Boundary_Vars ,ONLY: nComputeNodeSurfOutputSides,nGlobalOutputSides, nSurfBC,N_SurfVDL,nVarSurfData
USE MOD_Particle_Boundary_Vars ,ONLY: offsetComputeNodeSurfOutputSide, SurfBCName, nComputeNodeSurfSides
USE MOD_Particle_Boundary_Vars ,ONLY: SurfSide2GlobalSide, GlobalSide2SurfSide
USE MOD_Particle_Boundary_Vars ,ONLY: SurfTotalSideOnNode
USE MOD_HDF5_Output            ,ONLY: WriteAttributeToHDF5,WriteArrayToHDF5,WriteHDF5Header
USE MOD_Mesh_Vars              ,ONLY: MeshFile,nBCSides,offSetElem,SideToElem,BC,Boundarytype
USE MOD_Particle_Mesh_Vars     ,ONLY: SideInfo_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: mySurfRank
#if USE_MPI
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_SHARED,MPI_COMM_LEADERS_SURF
USE MOD_MPI_Shared_Vars        ,ONLY: myComputeNodeRank
USE MOD_Particle_Boundary_Vars ,ONLY: nGlobalSurfSides
#endif /*USE_MPI*/
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound
USE MOD_ChangeBasis            ,ONLY: ChangeBasis2D
USE MOD_Interpolation_Vars     ,ONLY: NMax,PREF_VDM,NodeType,NMin,NodeTypeVISU
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID
USE MOD_DG_Vars                ,ONLY: N_DG_Mapping
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Globals_Vars           ,ONLY: ProjectName
#if defined(PARTICLES) && USE_HDG
USE MOD_AnalyzeField_HDG       ,ONLY: CalculateElectricPotentialAndFieldBoundaryVDL
#endif /*defined(PARTICLES) && USE_HDG*/
USE MOD_Interpolation          ,ONLY: GetVandermonde,GetNodesAndWeights
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
REAL,INTENT(IN)                :: OutputTime
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)                  :: FileName
CHARACTER(LEN=255)                  :: Statedummy
CHARACTER(LEN=255)                  :: H5_Name
CHARACTER(LEN=255),ALLOCATABLE      :: Str2DVarNames(:)
INTEGER                             :: GlobalSideID, iSurfSide, OutputCounter, Nloc,BCSideID,iLocSide,iElem,BCType,iBC
INTEGER                             :: GlobalElemID,iPartBound,GlobalNonUniqueSideID
INTEGER                             :: nSurfSample
#if USE_MPI
INTEGER                             :: MessageSize
#endif /*USE_MPI*/
REAL                                :: tstart,tend
REAL, ALLOCATABLE                   :: helpArray(:,:,:,:), helpArray2(:,:,:,:)

TYPE InterpolateMe
  REAL,ALLOCATABLE  :: densityVISU(:,:,:,:)
  REAL,ALLOCATABLE  :: J_NAnalyze(:,:,:,:)
  REAL,ALLOCATABLE  :: Coords_NAnalyze(:,:,:,:)
  REAL,ALLOCATABLE  :: xIP_VISU(:),wIP_VISU(:)
  REAL,ALLOCATABLE  :: Vdm_N_EQ_SurfaceData(:,:) ! < Vandermonde mapping from NodeType to equidistant (visu) node set
END TYPE InterpolateMe

TYPE(InterpolateMe),ALLOCATABLE :: NtonSurfSample(:)
!===================================================================================================================================

! nodes without surfaces do not take part in this routine
IF (.NOT.SurfTotalSideOnNode) RETURN

#if defined(PARTICLES) && USE_HDG
! Prolong-to-face electric potential and electric field
CALL CalculateElectricPotentialAndFieldBoundaryVDL()
#endif /*defined(PARTICLES) && USE_HDG*/

! Build mapping from NodeType to equidistant (visu) node set
#if !(PP_NodeType==1 || PP_NodeType==2)
CALL abort(__STAMP__,'Output of VDL surface data only implemented for Gauss or Gauss-Lobatto nodes!')
#endif /*!(PP_NodeType==1 || PP_NodeType==2)*/
! Allocate type for all polynomial degrees from NMin to NMax
ALLOCATE(NtonSurfSample(NMin:NMax))

DO Nloc = NMin, NMax
  ! Allocate all arrays for Nloc
  ALLOCATE(NtonSurfSample(Nloc)%densityVISU(1,0:Nloc,0:Nloc,0:Nloc))
  ALLOCATE(NtonSurfSample(Nloc)%J_NAnalyze(1,0:Nloc,0:Nloc,0:Nloc))
  ALLOCATE(NtonSurfSample(Nloc)%Coords_NAnalyze(3,0:Nloc,0:Nloc,0:Nloc))
  ALLOCATE(NtonSurfSample(Nloc)%xIP_VISU(0:Nloc))
  ALLOCATE(NtonSurfSample(Nloc)%wIP_VISU(0:Nloc))
  ALLOCATE(NtonSurfSample(Nloc)%Vdm_N_EQ_SurfaceData(0:Nloc,0:Nloc))

  ! Allocate and determine Vandermonde mapping from NodeType to equidistant (visu) node set
  CALL GetVandermonde(Nloc, NodeType, Nloc, NodeTypeVISU, NtonSurfSample(Nloc)%Vdm_N_EQ_SurfaceData, modal=.FALSE.)
  ! Required only for integration
  !CALL GetNodesAndWeights(Nloc, NodeTypeVISU, NtonSurfSample(Nloc)%xIP_VISU, wIP=NtonSurfSample(Nloc)%wIP_VISU)
END DO ! Nloc = 1, NMax

nSurfSample = NMax+1

ALLOCATE(helpArray(nVarSurfData,1:nSurfSample,1:nSurfSample,nComputeNodeSurfTotalSides))
helpArray = 0. ! Must be zero because it will be summed up over MPI

! Loop over all local boundary sides
DO BCSideID=1,nBCSides

  ! Exclude periodic sides
  iBC = BC(BCSideID)
  BCType = Boundarytype(iBC,BC_TYPE)
  IF(BCType.EQ.1) CYCLE ! Skip periodic sides

  ! Exclude non-VDL boundaries
  iPartBound = PartBound%MapToPartBC(iBC)
  IF(ABS(PartBound%PermittivityVDL(iPartBound)).GT.0.0)THEN
    iElem                 = SideToElem(S2E_ELEM_ID,BCSideID)
    GlobalElemID          = iElem + offsetElem
    Nloc                  = N_DG_Mapping(2,GlobalElemID)
    iLocSide              = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
    GlobalNonUniqueSideID = GetGlobalNonUniqueSideID(GlobalElemID,iLocSide)
    iSurfSide             = GlobalSide2SurfSide(SURF_SIDEID,GlobalNonUniqueSideID)
    ! Map from Nloc to NMax
    IF(Nloc.EQ.NMax)THEN
      helpArray(1:nVarSurfData-2,1:nSurfSample,1:nSurfSample,iSurfSide) = N_SurfVDL(BCSideID)%U(1:nVarSurfData-2,0:Nloc,0:Nloc)
    ELSE
      ! From low to high
      CALL ChangeBasis2D(nVarSurfData-2, Nloc, NMax, PREF_VDM(Nloc,NMax)%Vdm ,&
             N_SurfVDL(BCSideID)%U(1:nVarSurfData-2,0:Nloc,0:Nloc),&
                         helpArray(1:nVarSurfData-2,1:nSurfSample,1:nSurfSample,iSurfSide))
    END IF ! Nloc.EQ.NMax

    ! Map from Gauss or Gauss-Lobatto nodes to equidistant VISU for .h5 storage
    CALL ChangeBasis2D(nVarSurfData-2, NMax, NMax, NtonSurfSample(NMax)%Vdm_N_EQ_SurfaceData ,&
                                      helpArray(1:nVarSurfData-2,1:nSurfSample,1:nSurfSample,iSurfSide) ,&
                                      helpArray(1:nVarSurfData-2,1:nSurfSample,1:nSurfSample,iSurfSide) )

    ! Set iBC and iPartBound
    helpArray(nVarSurfData-1,1:nSurfSample,1:nSurfSample,iSurfSide) = REAL(iBC)        ! iBC
    helpArray(nVarSurfData  ,1:nSurfSample,1:nSurfSample,iSurfSide) = REAL(iPartBound) ! iPartBound
  END IF ! ABS(PartBound%PermittivityVDL(iPartBound)).GT.0.0
END DO ! BCSideID=1,nBCSides

! Deallocate NtonSurfSample before some processors leave
DO Nloc = NMin, NMax
  DEALLOCATE(NtonSurfSample(Nloc)%densityVISU)
  DEALLOCATE(NtonSurfSample(Nloc)%J_NAnalyze)
  DEALLOCATE(NtonSurfSample(Nloc)%Coords_NAnalyze)
  DEALLOCATE(NtonSurfSample(Nloc)%xIP_VISU)
  DEALLOCATE(NtonSurfSample(Nloc)%wIP_VISU)
  DEALLOCATE(NtonSurfSample(Nloc)%Vdm_N_EQ_SurfaceData)
END DO ! Nloc = 1, NMax
DEALLOCATE(NtonSurfSample)

#if USE_MPI
! collect the information from the proc-local arrays in the compute-node array
MessageSize = nVarSurfData*nSurfSample*nSurfSample*nComputeNodeSurfTotalSides

IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,helpArray,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
ELSE
  CALL MPI_REDUCE(helpArray   ,0        ,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
ENDIF

! Return if not a sampling leader
IF (MPI_COMM_LEADERS_SURF.EQ.MPI_COMM_NULL) RETURN
CALL MPI_BARRIER(MPI_COMM_LEADERS_SURF,iERROR)

! Return if no sampling sides
IF (nGlobalSurfSides.EQ.0) RETURN
#endif /*USE_MPI*/

IF (mySurfRank.EQ.0) THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE VDL SurfState TO HDF5 FILE...'
  tstart=LOCALTIME()
END IF

ALLOCATE(helpArray2(nVarSurfData,1:nSurfSample,1:nSurfSample,nComputeNodeSurfOutputSides))
helpArray2 = 0.

OutputCounter = 0

DO iSurfSide = 1,nComputeNodeSurfSides
  !================== INNER BC CHECK
  GlobalSideID = SurfSide2GlobalSide(SURF_SIDEID,iSurfSide)
  IF(SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID).GT.0) THEN
    IF(GlobalSideID.LT.SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID)) THEN
      ! TODO: INNER BC requires special treatment for the calculation of the VDL potential and corresponding electric field
      ! SurfSideNb = GlobalSide2SurfSide(SURF_SIDEID,SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID))
    ELSE
      CYCLE
    END IF
  END IF
  !================== ROTATIONALLY PERIODIC BC CHECK
  IF(PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,GlobalSideID))).EQ.PartBound%RotPeriodicBC) CYCLE
  !================== INTER PLANE BC CHECK
  IF(PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,GlobalSideID))).EQ.PartBound%RotPeriodicInterPlaneBC) CYCLE
  OutputCounter = OutputCounter + 1
  helpArray2(:,:,:,OutputCounter) = helpArray(:,:,:,iSurfSide)
END DO

FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_VDLSurfState',OutputTime))//'.h5'
SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO') '['//TRIM(FileName)//'] ...'

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
#if USE_MPI
IF (mySurfRank.EQ.0) THEN
#endif /*USE_MPI*/
  CALL OpenDataFile(FileName,create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
  Statedummy = 'VDLSurfState'
  ! Write file header
  CALL WriteHDF5Header(Statedummy,File_ID)
  CALL WriteAttributeToHDF5(File_ID , 'DSMC_nSurfSample' , 1       , IntegerScalar = NMax+1       )
  CALL WriteAttributeToHDF5(File_ID , 'MeshFile'         , 1       , StrScalar     = (/TRIM(MeshFile)/) )
  CALL WriteAttributeToHDF5(File_ID , 'BC_Surf'          , nSurfBC , StrArray      = SurfBCName         )
  CALL WriteAttributeToHDF5(File_ID , 'N'                , 1       , IntegerScalar = NMax        )
  CALL WriteAttributeToHDF5(File_ID , 'NodeType'         , 1       , StrScalar     = (/NodeType/)   )
  CALL WriteAttributeToHDF5(File_ID , 'Time'             , 1       , RealScalar    = OutputTime                 )

  ALLOCATE(Str2DVarNames(1:nVarSurfData))
  ! fill varnames for total values
  Str2DVarNames(1) ='PhiF_From_E'
  Str2DVarNames(2) ='EFieldx'
  Str2DVarNames(3) ='EFieldy'
  Str2DVarNames(4) ='EFieldz'
  Str2DVarNames(5) ='PhiF_Max'
  Str2DVarNames(6) ='E_From_PhiF_Maxx'
  Str2DVarNames(7) ='E_From_PhiF_Maxy'
  Str2DVarNames(8) ='E_From_PhiF_Maxz'
  Str2DVarNames(9) ='PhiF_From_Currents'
  Str2DVarNames(10) ='E_From_PhiF_From_Currentsx'
  Str2DVarNames(11) ='E_From_PhiF_From_Currentsy'
  Str2DVarNames(12) ='E_From_PhiF_From_Currentsz'
  Str2DVarNames(nVarSurfData-1) ='iBC'
  Str2DVarNames(nVarSurfData) ='iPartBound'
  CALL WriteAttributeToHDF5(File_ID,'VarNamesSurface',nVarSurfData,StrArray=Str2DVarNames)

  CALL CloseDataFile()
  DEALLOCATE(Str2DVarNames)
#if USE_MPI
END IF
CALL MPI_BARRIER(MPI_COMM_LEADERS_SURF,iERROR)
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_LEADERS_SURF)
#else
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
#endif /*USE_MPI*/

WRITE(H5_Name,'(A)') 'SurfaceData'

! Note that the output happens on equidistant VISU nodes as piclas2vtk assumes that surface data has been equidistantly sampled
! when using the "SurfaceData" .h5 container, which is used here for convenience
ASSOCIATE (&
      nSurfSample    => INT(nSurfSample                     , IK)  , &
      nGlobalSides   => INT(nGlobalOutputSides              , IK)  , &
      nLocalSides    => INT(nComputeNodeSurfOutputSides     , IK)  , &
      offsetSurfSide => INT(offsetComputeNodeSurfOutputSide , IK)  , &
      nVarSurfData8  => INT(nVarSurfData                    , IK))

  ! WARNING: Only the sampling leaders write the data to .h5
  CALL WriteArrayToHDF5(DataSetName=H5_Name        , rank=4      , &
                        nValGlobal =(/nVarSurfData8, nSurfSample , nSurfSample , nGlobalSides/)   , &
                        nVal       =(/nVarSurfData8, nSurfSample , nSurfSample , nLocalSides/)  , &
                        offset     =(/0_IK         , 0_IK        , 0_IK        , offsetSurfSide/) , &
                        collective =.FALSE.        , &
                        RealArray=helpArray2(1:nVarSurfData,1:nSurfSample,1:nSurfSample,1:nLocalSides))
END ASSOCIATE

DEALLOCATE(helpArray)
DEALLOCATE(helpArray2)

CALL CloseDataFile()

IF (mySurfRank.EQ.0) THEN
  tend=LOCALTIME()
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',tend-tstart,'s]'
END IF

END SUBROUTINE WriteSurfVDLToHDF5
#endif /*defined(PARTICLES)*/
#endif /*USE_HDG*/

END MODULE MOD_HDF5_Output_Fields_HDG