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

MODULE MOD_HDF5_Output_Fields
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MOD_io_HDF5
USE MOD_HDF5_output
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
#if USE_HDG
#if defined(PARTICLES)
PUBLIC :: WriteBRAverageElemToHDF5
PUBLIC :: WriteSurfVDLToHDF5
#endif /*defined(PARTICLES)*/
#else
PUBLIC :: WritePMLzetaGlobalToHDF5
#endif /*USE_HDG*/

PUBLIC :: WriteDielectricGlobalToHDF5
#if (PP_nVar==8)
PUBLIC :: WritePMLDataToHDF5
#endif
PUBLIC :: WriteErrorNormsToHDF5
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/
PUBLIC :: WriteBGFieldToHDF5,WriteBGFieldAnalyticToHDF5
!===================================================================================================================================

CONTAINS


#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
SUBROUTINE WriteDielectricGlobalToHDF5()
!===================================================================================================================================
! write DielectricGlobal field to HDF5 file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Dielectric_Vars ,ONLY: DielectricGlobal
USE MOD_Dielectric_Vars ,ONLY: DielectricVol,isDielectricElem,ElemToDielectric
USE MOD_Mesh_Vars       ,ONLY: MeshFile,nGlobalElems,offsetElem
USE MOD_Globals_Vars    ,ONLY: ProjectName
USE MOD_io_HDF5
USE MOD_ChangeBasis     ,ONLY: ChangeBasis3D
USE MOD_DG_vars                ,ONLY: N_DG_Mapping
USE MOD_Interpolation_Vars     ,ONLY: PREF_VDM,Nmax
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER   :: N_variables=2
CHARACTER(LEN=255),ALLOCATABLE  :: StrVarNames(:)
CHARACTER(LEN=255)  :: FileName
REAL                :: StartT,EndT
REAL                :: OutputTime
INTEGER             :: iElem, Nloc
REAL, ALLOCATABLE   :: DielectricVolTmp(:,:,:,:)
!===================================================================================================================================
#if USE_LOADBALANCE
IF(PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/
! create global Eps field for parallel output of Eps distribution
ALLOCATE(DielectricGlobal(1:N_variables,0:Nmax,0:Nmax,0:Nmax,1:PP_nElems))
ALLOCATE(StrVarNames(1:N_variables))
StrVarNames(1)='DielectricEpsGlobal'
StrVarNames(2)='DielectricMuGlobal'
DielectricGlobal=0.

DO iElem=1,PP_nElems
  IF(isDielectricElem(iElem))THEN
    Nloc = N_DG_Mapping(2,iElem+offSetElem)
    IF (Nloc.EQ.Nmax) THEN
      DielectricGlobal(1,:,:,:,iElem)=DielectricVol(ElemToDielectric(iElem))%DielectricEps(:,:,:)
      DielectricGlobal(2,:,:,:,iElem)=DielectricVol(ElemToDielectric(iElem))%DielectricMu( :,:,:)
    ELSE
      ALLOCATE(DielectricVolTmp(2,0:Nloc,0:Nloc,0:Nloc))
      DielectricVolTmp(1,0:Nloc,0:Nloc,0:Nloc)= DielectricVol(ElemToDielectric(iElem))%DielectricEps(:,:,:)
      DielectricVolTmp(2,0:Nloc,0:Nloc,0:Nloc)= DielectricVol(ElemToDielectric(iElem))%DielectricMu(:,:,:)
      CALL ChangeBasis3D(2,Nloc,NMax,PREF_VDM(Nloc,NMax)%Vdm, &
          DielectricVolTmp(1:2 , 0:Nloc , 0:Nloc , 0:Nloc) , &
          DielectricGlobal(1:2 , 0:NMax , 0:NMax , 0:NMax  , iElem))
      DEALLOCATE(DielectricVolTmp)
    END IF
  END IF
END DO!iElem
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' WRITE DielectricGlobal TO HDF5 FILE...'
GETTIME(StartT)
OutputTime=0.0
! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_DielectricGlobal',OutputTime))//'.h5'
IF(MPIRoot) CALL GenerateFileSkeleton('DielectricGlobal',N_variables,StrVarNames,TRIM(MeshFile),OutputTime)
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif
IF(MPIRoot)THEN
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)
  CALL WriteAttributeToHDF5(File_ID,'VarNamesDielectricGlobal',N_variables,StrArray=StrVarNames)
  CALL CloseDataFile()
END IF ! MPIRoot

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
        nGlobalElems    => INT(nGlobalElems,IK)    ,&
        PP_nElems       => INT(PP_nElems,IK)       ,&
        N_variables     => INT(N_variables,IK)     ,&
        N               => INT(Nmax,IK)            ,&
        offsetElem      => INT(offsetElem,IK)      )
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
                          DataSetName='DG_Solution' , rank=5 , &
                          nValGlobal =(/N_variables , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
                          nVal       =(/N_variables , N+1_IK , N+1_IK , N+1_IK , PP_nElems   /) , &
                          offset     =(/       0_IK , 0_IK   , 0_IK   , 0_IK   , offsetElem  /) , &
                          collective =.TRUE.        , RealArray=DielectricGlobal)
END ASSOCIATE
GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
SDEALLOCATE(DielectricGlobal)
SDEALLOCATE(StrVarNames)
END SUBROUTINE WriteDielectricGlobalToHDF5


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
USE MOD_Globals_Vars     ,ONLY: ProjectName
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
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_BRAverageElem',OutputTime))//'.h5'
IF(MPIRoot) CALL GenerateFileSkeleton('BRAverageElem',N_variables,StrVarNames,TRIM(MeshFile),OutputTime)
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
USE MOD_Photon_TrackingVars    ,ONLY: PhotonSampWall
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound
USE MOD_RayTracing_Vars        ,ONLY: Ray
USE MOD_ChangeBasis            ,ONLY: ChangeBasis2D
USE MOD_Interpolation_Vars     ,ONLY: NMax,PREF_VDM,NodeType,NMin,NodeTypeVISU
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID
USE MOD_DG_Vars                ,ONLY: N_DG_Mapping
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Globals_Vars           ,ONLY: ProjectName
USE MOD_AnalyzeField           ,ONLY: CalculateElectricPotentialAndFieldBoundaryVDL
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
CHARACTER(LEN=255)             :: FileName
CHARACTER(LEN=255)                  :: Statedummy
CHARACTER(LEN=255)                  :: H5_Name
CHARACTER(LEN=255),ALLOCATABLE      :: Str2DVarNames(:)
INTEGER                             :: GlobalSideID, iSurfSide, OutputCounter, SurfSideNb, p, q,Nloc,BCSideID,iLocSide,iElem,BCType
INTEGER                             :: GlobalElemID,iPartBound,GlobalNonUniqueSideID
INTEGER                             :: nSurfSample, MessageSize
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

! Prolong-to-face electric potential and electric field
CALL CalculateElectricPotentialAndFieldBoundaryVDL()

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
helpArray = 0.

! Loop over all local boundary sides
DO BCSideID=1,nBCSides

  ! Exclude periodic sides
  BCType = Boundarytype(BC(BCSideID),BC_TYPE)
  IF(BCType.EQ.1) CYCLE ! Skip periodic side

  ! Exclude non-VDL boundaries
  iPartBound = PartBound%MapToPartBC(BC(BCSideID))
  IF(ABS(PartBound%PermittivityVDL(iPartBound)).GT.0.0)THEN
    iElem                 = SideToElem(S2E_ELEM_ID,BCSideID)
    GlobalElemID          = iElem + offsetElem
    Nloc                  = N_DG_Mapping(2,GlobalElemID)
    iLocSide              = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
    GlobalNonUniqueSideID = GetGlobalNonUniqueSideID(GlobalElemID,iLocSide)
    iSurfSide             = GlobalSide2SurfSide(SURF_SIDEID,GlobalNonUniqueSideID)
    ! Map from Nloc to NMax
    IF(Nloc.EQ.NMax)THEN
      helpArray(1:nVarSurfData,1:nSurfSample,1:nSurfSample,iSurfSide) = N_SurfVDL(BCSideID)%U(1:nVarSurfData,0:Nloc,0:Nloc)
    ELSE
      ! From low to high
      CALL ChangeBasis2D(nVarSurfData, Nloc, NMax, PREF_VDM(Nloc,NMax)%Vdm ,&
             N_SurfVDL(BCSideID)%U(1:nVarSurfData,0:Nloc,0:Nloc),&
                         helpArray(1:nVarSurfData,1:nSurfSample,1:nSurfSample,iSurfSide))
    END IF ! Nloc.EQ.NMax

    ! Map from Gauss or Gauss-Lobatto nodes to equidistant VISU for .h5 storage
    CALL ChangeBasis2D(nVarSurfData, NMax, NMax, NtonSurfSample(NMax)%Vdm_N_EQ_SurfaceData ,&
                                      helpArray(1:nVarSurfData,1:nSurfSample,1:nSurfSample,iSurfSide) ,&
                                      helpArray(1:nVarSurfData,1:nSurfSample,1:nSurfSample,iSurfSide) )
  END IF ! ABS(PartBound%PermittivityVDL(iPartBound)).GT.0.0
END DO ! BCSideID=1,nBCSides

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


#else
SUBROUTINE WritePMLzetaGlobalToHDF5()
!===================================================================================================================================
! write PMLzetaGlobal field to HDF5 file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_PML_Vars           ,ONLY: PMLzetaGlobal,PMLzeta0,isPMLElem,ElemToPML,PML
USE MOD_Mesh_Vars          ,ONLY: MeshFile,nGlobalElems,offsetElem
USE MOD_Globals_Vars       ,ONLY: ProjectName
USE MOD_io_HDF5
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
USE MOD_PML_Vars           ,ONLY: nPMLElems,PMLToElem
USE MOD_Interpolation_Vars ,ONLY: PREF_VDM,Nmax
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_DG_Vars            ,ONLY: N_DG_Mapping
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER   :: N_variables=3
CHARACTER(LEN=255),ALLOCATABLE  :: StrVarNames(:)
CHARACTER(LEN=255)  :: FileName
REAL                :: StartT,EndT
REAL                :: OutputTime
INTEGER             :: iElem,iPMLElem,Nloc
!===================================================================================================================================
#if USE_LOADBALANCE
IF(PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/
! create global zeta field for parallel output of zeta distribution
ALLOCATE(PMLzetaGlobal(1:N_variables,0:NMax,0:NMax,0:NMax,1:PP_nElems))
ALLOCATE(StrVarNames(1:N_variables))
StrVarNames(1)='PMLzetaGlobalX'
StrVarNames(2)='PMLzetaGlobalY'
StrVarNames(3)='PMLzetaGlobalZ'
PMLzetaGlobal=0.
IF(.NOT.ALMOSTZERO(PMLzeta0))THEN
  DO iPMLElem=1,nPMLElems
    iElem = PMLToElem(iPMLElem)
    Nloc  = N_DG_Mapping(2,iElem+offSetElem)
    IF(Nloc.EQ.Nmax)THEN
      PMLzetaGlobal(:,:,:,:,iElem) = PML(iPMLElem)%zeta(:,:,:,:)/PMLzeta0
    ELSE
      CALL ChangeBasis3D(3,Nloc,NMax,PREF_VDM(Nloc,NMax)%Vdm, &
          PML(iPMLElem)%zeta(1:3 , 0:Nloc , 0:Nloc , 0:Nloc)/PMLzeta0 , &
               PMLzetaGlobal(1:3 , 0:NMax , 0:NMax , 0:NMax           , iElem))
    END IF ! Nloc.Eq.Nmax
  END DO
END IF

SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE PMLZetaGlobal TO HDF5 FILE...'
GETTIME(StartT)
OutputTime=0.0
! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_PMLZetaGlobal',OutputTime))//'.h5'
IF(MPIRoot) CALL GenerateFileSkeleton('PMLZetaGlobal',N_variables,StrVarNames,TRIM(MeshFile),OutputTime)
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif
  CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)
CALL WriteAttributeToHDF5(File_ID,'VarNamesPMLzetaGlobal',N_variables,StrArray=StrVarNames)
CALL CloseDataFile()

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
        nGlobalElems    => INT(nGlobalElems,IK)    ,&
        PP_nElems       => INT(PP_nElems,IK)       ,&
        N_variables     => INT(N_variables,IK)     ,&
        N               => INT(NMax,IK)            ,&
        offsetElem      => INT(offsetElem,IK)      )
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
                          DataSetName='DG_Solution' , rank=5 , &
                          nValGlobal =(/N_variables , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
                          nVal       =(/N_variables , N+1_IK , N+1_IK , N+1_IK , PP_nElems   /) , &
                          offset     =(/       0_IK , 0_IK   , 0_IK   , 0_IK   , offsetElem  /) , &
                          collective =.TRUE.        , RealArray=PMLzetaGlobal)
END ASSOCIATE
GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
SDEALLOCATE(PMLzetaGlobal)
SDEALLOCATE(StrVarNames)
END SUBROUTINE WritePMLzetaGlobalToHDF5
#endif /*USE_HDG*/


#if (PP_nVar==8)
SUBROUTINE WritePMLDataToHDF5(FileName)
!===================================================================================================================================
! Write additional (elementwise scalar) data to HDF5
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars          ,ONLY: offsetElem,nGlobalElems,nElems
USE MOD_PML_Vars           ,ONLY: DoPML,PMLToElem,nPMLElems,PMLnVar
USE MOD_DG_Vars            ,ONLY: U_N,N_DG_Mapping
USE MOD_Interpolation_Vars ,ONLY: PREF_VDM,Nmax
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)  :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
REAL,ALLOCATABLE               :: UPML(:,:,:,:,:)
INTEGER                        :: iPML,iElem,Nloc
!===================================================================================================================================

IF(DoPML)THEN
  ALLOCATE(StrVarNames(PMLnVar))
  StrVarNames( 1)='PML-ElectricFieldX-P1'
  StrVarNames( 2)='PML-ElectricFieldX-P2'
  StrVarNames( 3)='PML-ElectricFieldX-P3'
  StrVarNames( 4)='PML-ElectricFieldY-P4'
  StrVarNames( 5)='PML-ElectricFieldY-P5'
  StrVarNames( 6)='PML-ElectricFieldY-P6'
  StrVarNames( 7)='PML-ElectricFieldZ-P7'
  StrVarNames( 8)='PML-ElectricFieldZ-P8'
  StrVarNames( 9)='PML-ElectricFieldZ-P9'
  StrVarNames(10)='PML-MagneticFieldX-P10'
  StrVarNames(11)='PML-MagneticFieldX-P11'
  StrVarNames(12)='PML-MagneticFieldX-P12'
  StrVarNames(13)='PML-MagneticFieldY-P13'
  StrVarNames(14)='PML-MagneticFieldY-P14'
  StrVarNames(15)='PML-MagneticFieldY-P15'
  StrVarNames(16)='PML-MagneticFieldZ-P16'
  StrVarNames(17)='PML-MagneticFieldZ-P17'
  StrVarNames(18)='PML-MagneticFieldZ-P18'
  StrVarNames(19)='PML-PhiB-P19'
  StrVarNames(20)='PML-PhiB-P20'
  StrVarNames(21)='PML-PhiB-P21'
  StrVarNames(22)='PML-PsiE-P22'
  StrVarNames(23)='PML-PsiE-P23'
  StrVarNames(24)='PML-PsiE-P24'

  ALLOCATE(UPML(PMLnVar,0:Nmax,0:Nmax,0:Nmax,nElems))
  UPML=0.0
  DO iPML=1,nPMLElems
    iElem = PMLToElem(iPML)
    Nloc  = N_DG_Mapping(2,iElem+offSetElem)
    IF(Nloc.EQ.Nmax)THEN
      UPML(:,:,:,:,iElem) = U_N(iElem)%U2(:,:,:,:)
    ELSE
      CALL ChangeBasis3D(PMLnVar,Nloc,NMax,PREF_VDM(Nloc,NMax)%Vdm, &
          U_N(iElem)%U2(1:PMLnVar , 0:Nloc , 0:Nloc , 0:Nloc) , &
                   UPML(1:PMLnVar , 0:NMax , 0:NMax , 0:NMax  , iElem))
    END IF ! Nloc.Eq.Nmax
  END DO ! iElem = 1, nElems

  IF(MPIRoot)THEN
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttributeToHDF5(File_ID,'VarNamesPML',PMLnVar,StrArray=StrVarNames)
    CALL CloseDataFile()
  END IF

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nGlobalElems    => INT(nGlobalElems,IK)    ,&
      N               => INT(PP_N,IK)            ,&
      PMLnVar         => INT(PMLnVar,IK)         ,&
      PP_nElems       => INT(PP_nElems,IK)       ,&
      offsetElem      => INT(offsetElem,IK)      )
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
                          DataSetName = 'PML_Solution', rank = 5,&
                          nValGlobal  = (/PMLnVar , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
                          nVal        = (/PMLnVar , N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
                          offset      = (/0_IK    , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
                          collective  = .TRUE.,RealArray = UPML)
END ASSOCIATE

!  CALL WriteArrayToHDF5(DataSetName='PML_Solution', rank=5,&
!                      nValGlobal=(/5,N+1,N+1,N+1,nGlobalElems/),&
!                      nVal=      (/5,N+1,N+1,N+1,PP_nElems/),&
!                      offset=    (/0,      0,     0,     0,     offsetElem/),&
!                      collective=.TRUE., existing=.FALSE., RealArray=UPML)
!
!  CALL CloseDataFile()
  DEALLOCATE(UPML)
  DEALLOCATE(StrVarNames)
END IF ! DoPML


END SUBROUTINE WritePMLDataToHDF5
#endif


SUBROUTINE WriteErrorNormsToHDF5(OutputTime)
!===================================================================================================================================
! Output the exact solution, the L2 error and LInf error to (in NodeTypeGL = 'GAUSS-LOBATTO') to HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: ProjectName
USE MOD_Mesh_Vars          ,ONLY: offsetElem,nGlobalElems,MeshFile
USE MOD_io_HDF5
USE MOD_HDF5_output        ,ONLY: WriteArrayToHDF5,copy_userblock
USE MOD_Output_Vars        ,ONLY: UserBlockTmpFile,userblock_total_len
USE MOD_Interpolation_Vars ,ONLY: Uex,NAnalyze,NodeTypeGL
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
!USE MOD_Equation_Vars      ,ONLY: StrVarNames
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: OutputTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
INTEGER                        :: nVal
INTEGER,PARAMETER              :: UexDataSize=1
REAL                           :: StartT,EndT
!===================================================================================================================================
#if USE_LOADBALANCE
IF(PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' WRITE Analytic solution and L2/LInf norms TO HDF5 FILE...'
GETTIME(StartT)

! Create dataset attribute "VarNames"
ALLOCATE(StrVarNames(1:UexDataSize))
StrVarNames(1)='Phi'
!StrVarNames(2)='BG-MagneticFieldY'
!StrVarNames(3)='BG-MagneticFieldZ'

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
!FileName=TRIM(ProjectName)//'_ErrorNorms.h5'
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_ErrorNorms',OutputTime))//'.h5'
IF(MPIRoot) THEN
  CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.,userblockSize=userblock_total_len)
  ! Write file header
  CALL WriteHDF5Header('DG_Solution',File_ID)
  ! Write dataset properties "Time","MeshFile","NextFile","NodeType","VarNames"
  CALL WriteAttributeToHDF5(File_ID,'N',1,IntegerScalar=PP_N)
  CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=(/TRIM(MeshFile)/))
  CALL WriteAttributeToHDF5(File_ID,'NodeType',1,StrScalar=(/NodeTypeGL/))
  CALL WriteAttributeToHDF5(File_ID,'VarNames',UexDataSize,StrArray=StrVarNames)
  CALL WriteAttributeToHDF5(File_ID,'Time'    ,1,RealScalar=OutputTime)
  CALL CloseDataFile()
  ! Add userblock to hdf5-file
  CALL copy_userblock(TRIM(FileName)//C_NULL_CHAR,TRIM(UserblockTmpFile)//C_NULL_CHAR)
END IF
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)

nVal=nGlobalElems  ! For the MPI case this must be replaced by the global number of elements (sum over all procs)

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
  UexDataSize  => INT(UexDataSize,IK)   ,&
  N            => INT(NAnalyze,IK)         ,&
  PP_nElems    => INT(PP_nElems,IK)    ,&
  offsetElem   => INT(offsetElem,IK)   ,&
  nGlobalElems => INT(nGlobalElems,IK) )
CALL WriteArrayToHDF5(DataSetName='DG_Solution', rank=5 , &
                      nValGlobal=(/UexDataSize , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
                      nVal      =(/UexDataSize , N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
                      offset    =(/0_IK        , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
                      collective=.false., RealArray=Uex)
END ASSOCIATE

CALL CloseDataFile()

DEALLOCATE(StrVarNames)

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
END SUBROUTINE WriteErrorNormsToHDF5
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/


SUBROUTINE WriteBGFieldToHDF5(OutputTime)
!===================================================================================================================================
! Subroutine to write the BField numerical solution to HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: ProjectName
USE MOD_Mesh_Vars          ,ONLY: offsetElem,nGlobalElems, nElems,MeshFile
USE MOD_io_HDF5
USE MOD_HDF5_output        ,ONLY: copy_userblock
USE MOD_HDF5_Output        ,ONLY: WriteArrayToHDF5
USE MOD_Output_Vars        ,ONLY: UserBlockTmpFile,userblock_total_len
USE MOD_Interpolation_Vars ,ONLY: NodeType, N_BG, BGDataSize, BGType, NMax,PREF_VDM
USE MOD_SuperB_Vars        ,ONLY: UseTimeDepCoil,nTimePoints,BGFieldFrequency,BGFieldCurrent
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_DG_Vars            ,ONLY: N_DG_Mapping
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),OPTIONAL         :: OutputTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
INTEGER                        :: nVal, iTimePoint, iElem, Nloc
REAL                           :: StartT,EndT
REAL                           :: BGTemp(BGDataSize,0:NMax,0:NMax,0:NMax, PP_nElems)
REAL, ALLOCATABLE              :: BGTDepTemp(:,:,:,:,:,:)
!===================================================================================================================================
#if USE_LOADBALANCE
IF(PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/
! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
IF(PRESENT(OutputTime))THEN
  FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_BGField',OutputTime))//'.h5'
ELSE
  FileName=TRIM(ProjectName)//'_BGField.h5'
END IF ! PRESENT(OutputTime)

SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' WRITE BG-FIELD ['//TRIM(FileName)//'] TO HDF5 FILE...'
GETTIME(StartT)

IF(UseTimeDepCoil) ALLOCATE(BGTDepTemp(BGDataSize,0:NMax,0:NMax,0:NMax, PP_nElems,nTimePoints))

! Create dataset attribute "VarNames"
ALLOCATE(StrVarNames(1:BGDataSize))
IF(BGType.EQ.1) THEN
  StrVarNames(1)='BG-ElectricFieldX'
  StrVarNames(2)='BG-ElectricFieldY'
  StrVarNames(3)='BG-ElectricFieldZ'
ELSE IF(BGType.EQ.2) THEN
  StrVarNames(1)='BG-MagneticFieldX'
  StrVarNames(2)='BG-MagneticFieldY'
  StrVarNames(3)='BG-MagneticFieldZ'
ELSE IF(BGType.EQ.3) THEN
  StrVarNames(1)='BG-ElectricFieldX'
  StrVarNames(2)='BG-ElectricFieldY'
  StrVarNames(3)='BG-ElectricFieldZ'
  StrVarNames(4)='BG-MagneticFieldX'
  StrVarNames(5)='BG-MagneticFieldY'
  StrVarNames(6)='BG-MagneticFieldZ'
END IF

IF(MPIRoot) THEN
  CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.,userblockSize=userblock_total_len)
  ! Write file header
  CALL WriteHDF5Header('BGField',File_ID) ! File_Type='BGField'
  ! Write dataset properties "Time","MeshFile","NextFile","NodeType","VarNames"
  CALL WriteAttributeToHDF5(File_ID   , 'N'                    , 1          , IntegerScalar=NMax)
  CALL WriteAttributeToHDF5(File_ID   , 'MeshFile'             , 1          , StrScalar=(/TRIM(MeshFile)/))
  CALL WriteAttributeToHDF5(File_ID   , 'NodeType'             , 1          , StrScalar=(/NodeType/))
  CALL WriteAttributeToHDF5(File_ID   , 'VarNames'             , BGDataSize , StrArray=StrVarNames)
  CALL WriteAttributeToHDF5(File_ID   , 'Time'                 , 1          , RealScalar=0.)
  IF(UseTimeDepCoil)THEN
    CALL WriteAttributeToHDF5(File_ID , 'BGFieldFrequency'     , 1          , RealScalar=BGFieldFrequency)
    CALL WriteAttributeToHDF5(File_ID , 'BGFieldCurrent'       , 1          , RealScalar=BGFieldCurrent)
    CALL WriteAttributeToHDF5(File_ID , 'BGFieldTimeDependent' , 1          , LogicalScalar=.TRUE.)
  ELSE
    CALL WriteAttributeToHDF5(File_ID , 'BGFieldTimeDependent' , 1          , LogicalScalar=.FALSE.)
  END IF
  CALL CloseDataFile()
  ! Add userblock to hdf5-file
  CALL copy_userblock(TRIM(FileName)//C_NULL_CHAR,TRIM(UserblockTmpFile)//C_NULL_CHAR)
END IF
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)

nVal=nGlobalElems  ! For the MPI case this must be replaced by the global number of elements (sum over all procs)

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
  BGDataSize8  => INT(BGDataSize,IK)    ,&
  N            => INT(NMax,IK)          ,&
  PP_nElems    => INT(PP_nElems,IK)     ,&
  offsetElem   => INT(offsetElem,IK)    ,&
  nGlobalElems => INT(nGlobalElems,IK)  ,&
  nTimePoints8 => INT(nTimePoints,IK)    )

  DO iElem = 1, INT(PP_nElems)
    Nloc = N_DG_Mapping(2,iElem+offSetElem)
    IF(Nloc.EQ.Nmax)THEN
      BGTemp(:,:,:,:,iElem)   = N_BG(iElem)%BGField(:,:,:,:)
      IF(UseTimeDepCoil) THEN
        DO iTimePoint = 1, nTimePoints
          BGTDepTemp(:,:,:,:,iElem, iTimePoint)   = N_BG(iElem)%BGFieldTDep(:,:,:,:,iTimePoint)
        END DO
      END IF
    ELSE
      CALL ChangeBasis3D(BGDataSize,Nloc,NMax,PREF_VDM(Nloc,NMax)%Vdm, &
          N_BG(iElem)%BGField(1:BGDataSize , 0:Nloc , 0:Nloc , 0:Nloc) , &
                       BGTemp(1:BGDataSize , 0:NMax , 0:NMax , 0:NMax  , iElem))
      IF(UseTimeDepCoil) THEN
        DO iTimePoint = 1, nTimePoints
          CALL ChangeBasis3D(BGDataSize,Nloc,NMax,PREF_VDM(Nloc,NMax)%Vdm, &
              N_BG(iElem)%BGFieldTDep(1:BGDataSize , 0:Nloc , 0:Nloc , 0:Nloc ,         iTimePoint) , &
                           BGTDepTemp(1:BGDataSize , 0:NMax , 0:NMax , 0:NMax , iElem , iTimePoint))
        END DO
      END IF
    END IF ! Nloc.Eq.Nmax
  END DO ! iElem = 1, nElems

  IF(UseTimeDepCoil)THEN
    CALL WriteArrayToHDF5(DataSetName='DG_Solution'   , rank=6 , &
                          nValGlobal=(/BGDataSize8, N+1_IK , N+1_IK , N+1_IK , nGlobalElems, nTimePoints8/) , &
                          nVal      =(/BGDataSize8, N+1_IK , N+1_IK , N+1_IK , PP_nElems   , nTimePoints8/) , &
                          offset    =(/0_IK       , 0_IK   , 0_IK   , 0_IK   , offsetElem  , 0_IK       /) , &
                          collective=.false., RealArray=BGTDepTemp(1:BGDataSize,0:N,0:N,0:N,1:nElems,1:nTimePoints))
  ELSE
    CALL WriteArrayToHDF5(DataSetName='DG_Solution'   , rank=5 , &
                          nValGlobal=(/BGDataSize8, N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
                          nVal      =(/BGDataSize8, N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
                          offset    =(/0_IK       , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
                          collective=.false., RealArray=BGTemp(1:BGDataSize,0:N,0:N,0:N,1:nElems))
  END IF ! UseTimeDepCoil

END ASSOCIATE

CALL CloseDataFile()

DEALLOCATE(StrVarNames)

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
END SUBROUTINE WriteBGFieldToHDF5


SUBROUTINE WriteBGFieldAnalyticToHDF5()
!===================================================================================================================================
! Subroutine to write the BField analytical solution to HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: ProjectName
USE MOD_Mesh_Vars          ,ONLY: offsetElem,nGlobalElems, nElems,MeshFile
USE MOD_io_HDF5
USE MOD_HDF5_output        ,ONLY: WriteArrayToHDF5, copy_userblock
USE MOD_Output_Vars        ,ONLY: UserBlockTmpFile,userblock_total_len
USE MOD_Interpolation_Vars ,ONLY: N_BG, NodeType, BGDataSize, NMax, PREF_VDM
USE MOD_Restart_Vars       ,ONLY: RestartTime
USE MOD_DG_Vars            ,ONLY: N_DG_Mapping
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
INTEGER                        :: nVal, iElem, Nloc
REAL                           :: outputArray(1:BGDataSize,0:NMax,0:NMax,0:NMax,1:nElems)
REAL                           :: StartT,EndT
!===================================================================================================================================
#if USE_LOADBALANCE
IF(PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' WRITE BG-FIELD Analytic solution TO HDF5 FILE...'
GETTIME(StartT)

DO iElem = 1, INT(PP_nElems)
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  IF(Nloc.EQ.Nmax)THEN
    outputArray(:,:,:,:,iElem) = N_BG(iElem)%BGFieldAnalytic(:,:,:,:)
  ELSE
    CALL ChangeBasis3D(BGDataSize,Nloc,NMax,PREF_VDM(Nloc,NMax)%Vdm, &
        N_BG(iElem)%BGFieldAnalytic(1:BGDataSize , 0:Nloc , 0:Nloc , 0:Nloc ) , &
                        outputArray(1:BGDataSize , 0:NMax , 0:NMax , 0:NMax   , iElem))
  END IF ! Nloc.Eq.Nmax
END DO ! iElem = 1, nElems

! Create dataset attribute "VarNames"
ALLOCATE(StrVarNames(1:BGDataSize))
StrVarNames(1)='BG-MagneticFieldX'
StrVarNames(2)='BG-MagneticFieldY'
StrVarNames(3)='BG-MagneticFieldZ'

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(ProjectName)//'_BFieldAnalytic.h5'
IF(MPIRoot) THEN
  CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.,userblockSize=userblock_total_len)
  ! Write file header
  CALL WriteHDF5Header('BField',File_ID) ! File_Type='BField'
  ! Write dataset properties "Time","MeshFile","NextFile","NodeType","VarNames"
  CALL WriteAttributeToHDF5(File_ID,'N',1,IntegerScalar=Nmax)
  CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=(/TRIM(MeshFile)/))
  CALL WriteAttributeToHDF5(File_ID,'NodeType',1,StrScalar=(/NodeType/))
  CALL WriteAttributeToHDF5(File_ID,'VarNames',BGDataSize,StrArray=StrVarNames)
  CALL WriteAttributeToHDF5(File_ID,'Time'    ,1,RealScalar=RestartTime)
  CALL CloseDataFile()
  ! Add userblock to hdf5-file
  CALL copy_userblock(TRIM(FileName)//C_NULL_CHAR,TRIM(UserblockTmpFile)//C_NULL_CHAR)
END IF
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)

nVal=nGlobalElems  ! For the MPI case this must be replaced by the global number of elements (sum over all procs)

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
  BGDataSize   => INT(BGDataSize,IK)   ,&
  N            => INT(NMax,IK)         ,&
  PP_nElems    => INT(PP_nElems,IK)    ,&
  offsetElem   => INT(offsetElem,IK)   ,&
  nGlobalElems => INT(nGlobalElems,IK) )
CALL WriteArrayToHDF5(DataSetName='DG_Solution'    , rank=5 , &
                      nValGlobal=(/BGDataSize , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
                      nVal      =(/BGDataSize , N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
                      offset    =(/0_IK       , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
                      collective=.false., RealArray=outputArray)
END ASSOCIATE

CALL CloseDataFile()

DEALLOCATE(StrVarNames)

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
END SUBROUTINE WriteBGFieldAnalyticToHDF5


END MODULE MOD_HDF5_Output_Fields
