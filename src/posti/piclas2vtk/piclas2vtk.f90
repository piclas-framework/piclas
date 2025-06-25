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

!==================================================================================================================================
!> The PICLAS2VTK tool takes state files written during runtime by PICLAS in the .h5 format and converts them to .vtu files,
!> readable by ParaView. Supports parallel readin.
!> The state files can come from different calculations with different mesh files, equation systems, polynomial degrees and so on.
!> Two modes of usage: command line mode and parameter file mode.
!> In parameter file mode the usage is: piclas2vtk parameter.ini State1.h5 State2.h5 State3.h5 ...
!> In the parameter file the following can be specified:
!> - NVisu: Integer, polynomial degree of visualization basis
!> - NodeTypeVisu: String, node type of visualization basis
!> In command line mode, only the degree of the visualization basis can be directly specified, no parameter file is needed:
!> piclas2vtk --NVisu=INTEGER State1.h5 State2.h5 State3.h5 ...
!> All other options are set to their standard values.
!==================================================================================================================================
PROGRAM piclas2vtk
! MODULES
USE MOD_piclas2vtk_Vars
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_StringTools
USE MOD_Commandline_Arguments
USE MOD_IO_HDF5               ,ONLY: InitIOHDF5,DefineParametersIO
USE MOD_MPI                   ,ONLY: InitMPI
USE MOD_ReadInTools           ,ONLY: prms,PrintDefaultParameterFile
USE MOD_ReadInTools           ,ONLY: GETINT,GETSTR,GETLOGICAL
USE MOD_HDF5_Input            ,ONLY: OpenDataFile,CloseDataFile,ReadAttribute,File_ID,DatasetExists
USE MOD_HDF5_Input            ,ONLY: ISVALIDHDF5FILE
USE MOD_Mesh_Vars             ,ONLY: nGlobalElems,NGeo
USE MOD_Interpolation         ,ONLY: DefineParametersInterpolation,InitInterpolation
USE MOD_Mesh                  ,ONLY: DefineParametersMesh,InitMesh
USE MOD_Mesh_Tools            ,ONLY: InitElemNodeIDs
#ifdef PARTICLES
USE MOD_Particle_Mesh_Tools   ,ONLY: InitParticleGeometry
USE MOD_Particle_Mesh_Vars    ,ONLY: ConcaveElemSide_Shared,ElemSideNodeID_Shared,ElemMidPoint_Shared
USE MOD_Symmetry_Vars         ,ONLY: Symmetry
#endif /*PARTICLES*/
USE MOD_Particle_Mesh_Vars    ,ONLY: NodeCoords_Shared, ElemNodeID_Shared, NodeInfo_Shared
USE MOD_Mesh_Tools            ,ONLY: InitGetCNElemID, InitGetGlobalElemID
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars       ,ONLY: nComputeNodeTotalElems
#endif /*USE_MPI*/
USE MOD_Preproc
USE MOD_Mesh_ReadIn           ,ONLY: FinalizeMeshReadin
USE MOD_Mesh                  ,ONLY: FinalizeMesh
#if USE_MPI
#if defined(PARTICLES)
USE MOD_Particle_Mesh_Vars    ,ONLY: ConcaveElemSide_Shared_Win,ElemSideNodeID_Shared_Win,ElemMidPoint_Shared_Win
#endif /*defined(PARTICLES)*/
USE MOD_Particle_Mesh_Vars    ,ONLY: ElemNodeID_Shared_Win
USE MOD_MPI_Shared_vars       ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared
#endif /*USE_MPI*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Time                              ! Used to track computation time
CHARACTER(LEN=255)             :: NodeTypeVisuOut, InputStateFile, MeshFile, File_Type, DGSolutionDataset
INTEGER                        :: NVisu, iArgs, iArgsStart, TimeStampLength, iExt
LOGICAL                        :: CmdLineMode, NVisuDefault         ! In command line mode only NVisu is specified directly,
                                                                    ! otherwise a parameter file is needed
CHARACTER(LEN=2)               :: NVisuString                       ! String containing NVisu from command line option
CHARACTER(LEN=20)              :: fmtString                         ! String containing options for formatted write
LOGICAL                        :: DGSolutionExists, ElemDataExists, SurfaceDataExists, VisuParticles, PartDataExists, DMDDataExists
LOGICAL                        :: BGFieldExists, ExcitationDataExists, DVMSolutionExists
LOGICAL                        :: VisuAdaptiveInfo, AdaptiveInfoExists
LOGICAL                        :: ReadMeshFinished, ElemMeshInit, SurfMeshInit
LOGICAL                        :: ConvertPointToCellData
INTEGER                        :: iElem, iNode
INTEGER                        :: iMode,iErrorReturn,dmdSingleModeOutput,dmdMaximumModeOutput
CHARACTER(LEN=16)              :: hilf
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
CALL ParseCommandlineArguments()
!CALL DefineParametersMPI()
!CALL DefineParametersIO_HDF5()
! Define parameters for piclas2vtk
CALL prms%SetSection('piclas2vtk')
CALL prms%CreateStringOption( 'NodeTypeVisu'           , 'Node type of the visualization basis: VISU,GAUSS,GAUSS-LOBATTO,CHEBYSHEV-GAUSS-LOBATTO', 'VISU')
CALL prms%CreateIntOption(    'NVisu'                  , 'Number of points at which solution is sampled for visualization.')
CALL prms%CreateIntOption(    'NVisuAdd'               , 'p-adaption: increase local polynomial degree by NVisuAdd','0')
CALL prms%CreateLogicalOption('VisuParticles'          , 'Visualize particles (velocity, species, internal energy).', '.FALSE.')
CALL prms%CreateLogicalOption('VisuAdaptiveInfo'       , 'Visualize the sampled values utilized for the adaptive surface flux and porous BC (velocity, density, pumping speed)', '.FALSE.')
CALL prms%CreateLogicalOption('ConvertPointToCellData' , 'Visualize the DG solution as constant cell data using the VISU INNER nodes', '.FALSE.')
CALL prms%CreateIntOption(    'TimeStampLength'        , 'Length of the floating number time stamp', '21')
CALL prms%CreateIntOption(    'dmdSingleModeOutput'    , 'Convert only a single specific DMD mode.', '0')
CALL prms%CreateIntOption(    'dmdMaximumModeOutput'   , 'Convert all output DMD modes up to this number.', '1000')
CALL prms%CreateLogicalOption('meshCheckWeirdElements' , 'Abort when weird elements are found: it means that part of the element is turned inside-out. ','.TRUE.')
CALL DefineParametersIO()
CALL DefineParametersMesh()
CALL DefineParametersInterpolation()
#if USE_MPI
CALL DefineParametersMPIShared()
#endif /*USE_MPI*/

NVisuDefault = .FALSE.

! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF

IF (nArgs.LT.1) THEN
  ! If not, print out error message containing valid syntax
  CALL CollectiveStop(__STAMP__,&
  'ERROR - Please supply at least one .h5 files or a parameter file (or simply --NVisu=INTEGER) followed by h5 files!')
END IF

! Measure init duration
StartTime=PICLASTIME()

ParameterFile = Args(1)
! Check if first argument is a parameter file or that the NVisu argument has been specified
iExt=INDEX(ParameterFile,'.',BACK = .TRUE.) ! Position of file extension

IF ((nArgs.EQ.1).AND.(ParameterFile(iExt+1:iExt+2) .NE. 'h5')) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Please supply at least one .h5 file!')
END IF

CmdLineMode = .FALSE. ! Initialize
IF(ParameterFile(iExt+1:iExt+3) .EQ. 'ini') THEN
  ! Parameter file has been supplied
  CmdLineMode = .FALSE.
ELSE IF (ParameterFile(iExt+1:iExt+2) .EQ. 'h5') THEN
  NVisuDefault = .TRUE.
  CmdLineMode = .FALSE.
ELSE
  ! Check if the command line argument specifies NVisu
  IF (STRICMP(Args(1)(1:8),'--nvisu=')) THEN
    ! No Paramter file, but NVisu specified per command line option
    CmdLineMode = .TRUE.
  ELSE
  ! Neither parameter file nor NVisu have been specified
  CALL CollectiveStop(__STAMP__,&
    'ERROR - First argument must be a parameter file or NVisu must be specified per --NVisu=INTEGER.')
  END IF
END IF

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(A)')"                                   _..._"
SWRITE(UNIT_stdOut,'(A)')"                                .-'_..._''. .---.                         .-''-."
SWRITE(UNIT_stdOut,'(A)')"_________   _...._      .--.  .' .'      '.\|   |                       .' .-.  ).----.     .----.          ."
SWRITE(UNIT_stdOut,'(A)')"\        |.'      '-.   |__| / .'           |   |                      / .'  / /  \    \   /    /         .'|"
SWRITE(UNIT_stdOut,'(A)')" \        .'```'.    '. .--.. '             |   |                     (_/   / /    '   '. /'   /    .|  .'  |"
SWRITE(UNIT_stdOut,'(A)')"  \      |       \     \|  || |             |   |    __                    / /     |    |'    /   .' |_<    |"
SWRITE(UNIT_stdOut,'(A)')"   |     |        |    ||  || |             |   | .:--.'.         _       / /      |    ||    | .'     ||   | ____"
SWRITE(UNIT_stdOut,'(A)')"   |      \      /    . |  |. '             |   |/ |   \ |      .' |     . '       '.   `'   .''--.  .-'|   | \ .'"
SWRITE(UNIT_stdOut,'(A)')"   |     |\`'-.-'   .'  |  | \ '.          .|   |`' __ | |     .   | /  / /    _.-')\        /    |  |  |   |/  ."
SWRITE(UNIT_stdOut,'(A)')"   |     | '-....-'`    |__|  '. `._____.-'/|   | .'.''| |   .'.'| |//.' '  _.'.-''  \      /     |  |  |    /\  \"
SWRITE(UNIT_stdOut,'(A)')"  .'     '.                     `-.______ / '---'/ /   | |_.'.'.-'  //  /.-'_.'       '----'      |  '.'|   |  \  \"
SWRITE(UNIT_stdOut,'(A)')"'-----------'                            `       \ \._,\ '/.'   \_.'/    _.'                      |   / '    \  \  \"
SWRITE(UNIT_stdOut,'(A)')"                                                  `--'  `'         ( _.-'                         `'-' '------'  '---'"
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(132("="))')

! Set and read in parameters differently depending if piclas2vtk is invoked with a parameter file or not
IF(NVisuDefault.OR.CmdLineMode) THEN
  IF(NVisuDefault) THEN
    NVisu = 1
      ! Since we are not reading a parameter file, some properties of the prms object need to be set
    prms%maxNameLen  = 25 ! meshCheckWeirdElements is the longest option
    prms%maxValueLen = 23 ! CHEBYSHEV-GAUSS-LOBATTO is the longest possible value
  ELSE IF(CmdLineMode) THEN
    ! Read NVisu from the first command line argument
    NVisuString = TRIM(ParameterFile(9:LEN(TRIM(ParameterFile))))
    READ(NVisuString,'(I2.1)') NVisu
  END IF
  ! Formatted output
  WRITE(fmtString,*) prms%maxNameLen
  SWRITE(UNIT_stdOut,'(a3)', ADVANCE='NO')  " | "
  CALL set_formatting("blue")
  SWRITE(UNIT_stdOut,"(a"//fmtString//")", ADVANCE='NO') TRIM('NVisu')
  CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(a3)', ADVANCE='NO')  " | "
  WRITE(fmtString,*) prms%maxValueLen
  SWRITE(UNIT_stdOut,'(I'//fmtString//',A3)',ADVANCE='NO') NVisu,' | '
  CALL set_formatting("green")
  SWRITE(UNIT_stdOut,'(a7)', ADVANCE='NO')  "*CUSTOM"
  CALL clear_formatting()
  SWRITE(UNIT_stdOut,"(a3)") ' | '
ELSE
  ! Parse parameter file
  CALL prms%read_options(ParameterFile)
  ! Read in NVisu from parameter file
  NVisu            = GETINT('NVisu')                  ! Degree of visualization basis
END IF

! p-adaption: increase the cell-local polynomial degree by NVisuAdd
NVisuAdd            = GETINT('NVisuAdd')

! Set necessary parameters for piclas2vtk tool
! If no parameter file has been set, the standard values will be used
NodeTypeVisuOut  = GETSTR('NodeTypeVisu','VISU')    ! Node type of visualization basis
VisuParticles    = GETLOGICAL('VisuParticles','.FALSE.')
VisuAdaptiveInfo    = GETLOGICAL('VisuAdaptiveInfo','.FALSE.')
ConvertPointToCellData    = GETLOGICAL('ConvertPointToCellData','.FALSE.')
IF(ConvertPointToCellData) THEN
  ! Output of DG solution as cell data on the VISU INNER nodes
  PointToCellSwitch = 1
ELSE
  ! Output of DG solution as point data
  PointToCellSwitch = 0
END IF
! Initialization of I/O routines
CALL InitIOHDF5()
! Get length of the floating number time stamp
TimeStampLength = GETINT('TimeStampLength')
IF((TimeStampLength.LT.4).OR.(TimeStampLength.GT.30)) CALL abort(__STAMP__&
    ,'TimeStampLength cannot be smaller than 4 and not larger than 30')
WRITE(UNIT=TimeStampLenStr ,FMT='(I0)') TimeStampLength
WRITE(UNIT=TimeStampLenStr2,FMT='(I0)') TimeStampLength-4
#ifdef PARTICLES
Symmetry%Order=3
#endif
! DMD Stuff
dmdSingleModeOutput  = GETINT('dmdSingleModeOutput')
dmdMaximumModeOutput = GETINT('dmdMaximumModeOutput')

! Measure init duration
GETTIME(Time)
SWRITE(UNIT_stdOut,'(132("="))')
CALL DisplayMessageAndTime(Time-StartTime, ' INITIALIZATION DONE!', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
SWRITE(UNIT_stdOut,'(132("="))')

IF(NVisuDefault) THEN
  iArgsStart = 1
ELSE
  iArgsStart = 2
END IF

ReadMeshFinished = .FALSE.
ElemMeshInit = .FALSE.
SurfMeshInit = .FALSE.

#if USE_MPI
  CALL InitMPIShared()
#endif /*USE_MPI*/

CALL InitInterpolation(NVisu)

! Loop over remaining supplied .h5 files
DO iArgs = iArgsStart,nArgs
  InputStateFile = Args(iArgs)
  ! Check if the argument is a valid .h5 file
  IF(.NOT.ISVALIDHDF5FILE(InputStateFile)) CALL Abort(__STAMP__,'ERROR - Please supply only .h5 files after parameter file.')

  SWRITE(UNIT_stdOut,'(132("="))')
  SWRITE(UNIT_stdOut,'(A,I3,A,I3,A)') 'Processing state ',iArgs-iArgsStart+1,' of ',nArgs-iArgsStart+1,'...'

  ! Open .h5 file
  DGSolutionExists   = .FALSE.
  ElemDataExists     = .FALSE.
  SurfaceDataExists  = .FALSE.
  PartDataExists     = .FALSE.
  AdaptiveInfoExists = .FALSE.
  CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)
  ! Get the type of the .h5 file
  CALL ReadAttribute(File_ID,'File_Type',1,StrScalar=File_Type)
  ! Check which containers are present
  CALL DatasetExists(File_ID , 'DG_Solution'  , DGSolutionExists)
  IF(TRIM(File_Type).EQ.'PartStateBoundary'.OR.TRIM(File_Type).EQ.'PartStateLost') THEN
  ! Enable particle visualization for lost particle container (otherwise you get an empty file)
    DGSolutionExists = .FALSE.
    VisuParticles = .TRUE.
  END IF
  CALL DatasetExists(File_ID , 'ElemData'     , ElemDataExists)
  CALL DatasetExists(File_ID , 'ExcitationData', ExcitationDataExists)
  CALL DatasetExists(File_ID , 'AdaptiveInfo' , AdaptiveInfoExists)
  CALL DatasetExists(File_ID , 'SurfaceData'  , SurfaceDataExists)
  CALL DatasetExists(File_ID , 'PartData'     , PartDataExists)
  CALL DatasetExists(File_ID , 'BGField'      , BGFieldExists) ! deprecated , but allow for backward compatibility
  CALL DatasetExists(File_ID , 'DVM_Solution'  , DVMSolutionExists)
  CALL DatasetExists(File_ID , 'Mode_001_ElectricFieldX_Img'     , DMDDataExists)
  IF(BGFieldExists)THEN
    DGSolutionExists  = .TRUE.
    DGSolutionDataset = 'BGField'
  ELSE
    DGSolutionDataset = 'DG_Solution'
  END IF ! BGFieldExists
  ! Get the name of the mesh file
  CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile)
  CALL CloseDataFile()

  ! Read-in of the mesh
  IF(.NOT.ReadMeshFinished) THEN
    CALL InitMesh(-2,MeshFile_IN=MeshFile)
    CALL InitGetGlobalElemID()
    CALL InitGetCNElemID()
#if USE_MPI
    nComputeNodeTotalElems = nGlobalElems
#endif /*USE_MPI*/
#ifdef PARTICLES
    ! ElemSideNodeID_Shared is required for the connectivity for the conversion of SurfaceData
    CALL InitParticleGeometry()
#endif /*PARTICLES*/
    ! ElemNodeID_Shared is required for the connectivity for the conversion of ElemData
    CALL InitElemNodeIDs()
    ReadMeshFinished = .TRUE.
  END IF
  ! Build connectivity for element/volume output
  IF((ElemDataExists.OR.DVMSolutionExists).AND..NOT.ElemMeshInit) THEN
    CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)
    CALL ReadAttribute(File_ID,'nUniqueNodes',1,IntScalar=nUniqueNodes)
    CALL CloseDataFile()
    ALLOCATE(ElemUniqueNodeID(1:8,1:nGlobalElems))
    ALLOCATE(NodeCoords_Connect(1:3,1:nUniqueNodes))
    DO iElem = 1, nGlobalElems
      DO iNode = 1,8
        ElemUniqueNodeID(iNode,iElem)=ABS(NodeInfo_Shared(ElemNodeID_Shared(iNode,iElem)))
        NodeCoords_Connect(1:3,ElemUniqueNodeID(iNode,iElem)) = NodeCoords_Shared(1:3,ElemNodeID_Shared(iNode,iElem))
      END DO
    END DO
    IF(Ngeo.GT.NVisu) THEN
      CALL abort(__STAMP__,'ERROR: Ngeo as read-in from mesh is greater than the chosen NVisu! Ngeo: ',Ngeo)
    END IF
    ElemMeshInit = .TRUE.
  END IF
  ! Build connectivity for surface output
  IF(SurfaceDataExists.AND..NOT.SurfMeshInit) THEN
#ifdef PARTICLES
    CALL BuildSurfMeshConnectivity(InputStateFile)
    SurfMeshInit = .TRUE.
#else
    SWRITE(*,*) 'WARNING - Conversion of SurfaceData requires the compilation of piclas2vtk with PARTICLES set to ON!'
    SWRITE(*,*) 'WARNING - File contains SurfaceData but output is skipped.'
    SurfaceDataExists = .FALSE.
#endif
  END IF
  ! === DG_Solution (incl. BField etc.) ============================================================================================
  IF(DGSolutionExists) THEN
    IF(DMDDataExists)THEN
      IF(dmdSingleModeOutput.EQ.0)THEN
        DO iMode = 1, dmdMaximumModeOutput
          WRITE(UNIT=hilf,FMT='(I3.3)') iMode
          DGSolutionDataset = 'Mode_'//TRIM(hilf)//'_ElectricFieldX_Img'
          CALL ConvertDGSolution(InputStateFile,NVisu,NodeTypeVisuOut,File_Type,DGSolutionDataset,iErrorReturn)
          IF(iErrorReturn.NE.0) CONTINUE
        END DO ! iMode = 1, 1000
      ELSE
        WRITE(UNIT=hilf,FMT='(I3.3)') dmdSingleModeOutput
        DGSolutionDataset = 'Mode_'//TRIM(hilf)//'_ElectricFieldX_Img'
        CALL ConvertDGSolution(InputStateFile,NVisu,NodeTypeVisuOut,File_Type,DGSolutionDataset,iErrorReturn)
      END IF ! dmdSingleModeOutput.EQ.0
    ELSE
      CALL ConvertDGSolution(InputStateFile,NVisu,NodeTypeVisuOut,File_Type,DGSolutionDataset,iErrorReturn)
    END IF ! DMDDataExists
  END IF
  ! === ElemData ===================================================================================================================
  IF(ElemDataExists) THEN
    CALL ConvertElemData(InputStateFile,'ElemData','VarNamesAdd',iArgs)
  END IF
  ! === ElemData ===================================================================================================================
  IF(ExcitationDataExists) THEN
    CALL ConvertElemData(InputStateFile,'ExcitationData','VarNamesExci',iArgs)
  END IF
  ! === SurfaceData ================================================================================================================
  IF(SurfaceDataExists) THEN
    CALL ConvertSurfaceData(InputStateFile)
  END IF
  ! === PartData ===================================================================================================================
  IF(VisuParticles) THEN
    IF(PartDataExists) THEN
      CALL ConvertPartData(InputStateFile)
    END IF
  END IF
  ! === AdaptiveInfo
  IF(VisuAdaptiveInfo) THEN
    IF(AdaptiveInfoExists) THEN
      CALL ConvertElemData(InputStateFile,'AdaptiveInfo','VarNamesAdaptive',iArgs)
    END IF
  END IF
  ! === DVM_Solution ===============================================================================================================
  IF(DVMSolutionExists) THEN
    CALL ConvertDVMSolution(InputStateFile,NVisu,NodeTypeVisuOut,'DVM_Solution')
  END IF
END DO ! iArgs = 2, nArgs

! Finalize
IF(ReadMeshFinished)THEN
  CALL FinalizeMeshReadin(-2)
  CALL FinalizeMesh() ! Unlock and free N_DG_Mapping_Shared
  ! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
#endif /*USE_MPI*/

#if defined(PARTICLES)
#if USE_MPI
  ! InitParticleGeometry()
  CALL UNLOCK_AND_FREE(ConcaveElemSide_Shared_Win)
  CALL UNLOCK_AND_FREE(ElemSideNodeID_Shared_Win)
  CALL UNLOCK_AND_FREE(ElemMidPoint_Shared_Win)
#endif /*USE_MPI*/
  ! InitParticleGeometry
  ADEALLOCATE(ConcaveElemSide_Shared)
  ADEALLOCATE(ElemSideNodeID_Shared)
  ADEALLOCATE(ElemMidPoint_Shared)
#endif /*defined(PARTICLES)*/

#if USE_MPI
  ! From InitElemNodeIDs
  CALL UNLOCK_AND_FREE(ElemNodeID_Shared_Win)
  CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
#endif /*USE_MPI*/
  ADEALLOCATE(ElemNodeID_Shared)
END IF ! ReadMeshFinished


! Measure processing duration
GETTIME(Time)
#if USE_MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) CALL abort(__STAMP__,'MPI finalize error',iError)
#endif
SWRITE(UNIT_stdOut,'(132("="))')
CALL DisplayMessageAndTime(Time-StartTime, ' piclas2vtk FINISHED!', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
SWRITE(UNIT_stdOut,'(132("="))')

END PROGRAM piclas2vtk


!===================================================================================================================================
!> Subroutine to write 3D point data to VTK format
!===================================================================================================================================
SUBROUTINE WriteDataToVTK_PICLas(dim,data_size,FileString,nVar,VarNameVisu,nNodes,Coords,nElems,Array,ConnectInfo)
! MODULES
USE MOD_Globals
USE MOD_Particle_Boundary_Vars    ,ONLY: nSurfSample
USE MOD_VTK                       ,ONLY: CreateConnectivity
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: dim
INTEGER,INTENT(IN)            :: nVar,nElems,nNodes,data_size           !> Number of nodal output variables
REAL,INTENT(IN)               :: Coords(1:3,0:nSurfSample*(MERGE(1,0,nSurfSample.GT.1)),0:nSurfSample*(MERGE(1,0,nSurfSample.GT.1.AND.dim.GT.1)),0:nSurfSample*(MERGE(1,0,nSurfSample.GT.1.AND.dim.GT.2)),nNodes)                     !> Coordinates x, y and z
REAL,INTENT(IN)               :: Array(1:nVar,0:(nSurfSample-1)*(MERGE(1,0,nSurfSample.GT.1)),0:(nSurfSample-1)*(MERGE(1,0,nSurfSample.GT.1.AND.dim.GT.1)),0:(nSurfSample-1)*(MERGE(1,0,nSurfSample.GT.1.AND.dim.GT.2)),1:nElems)     !< Array with nVar properties for each coordinate
CHARACTER(LEN=*),INTENT(IN)   :: FileString                             !> Output file name
CHARACTER(LEN=*),INTENT(IN)   :: VarNameVisu(nVar)                      !> Variable names
INTEGER,INTENT(IN)            :: ConnectInfo(data_size,nElems)          !> Node connection information
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iVal,iElem,Offset,nBytes,nVTKPoints,nVTKCells,ivtk=44,iVar,iNode, int_size, ElemType
INTEGER                       :: iLen, str_len, NVisu, NVisuCoords,NPlot_p1_2,NVisuValue
CHARACTER(LEN=35)             :: StrOffset,TempStr1,TempStr2
CHARACTER(LEN=200)            :: Buffer, tmp, tmp2, VarNameString
CHARACTER(LEN=1)              :: lf, components_string
REAL(KIND=4)                  :: float
INTEGER,ALLOCATABLE           :: VarNameCombine(:), VarNameCombineLen(:), Vertex(:,:)
REAL                          :: StartT,EndT ! Timer
INTEGER,ALLOCATABLE,TARGET    :: nodeids(:)
!===================================================================================================================================
GETTIME(StartT)

SELECT CASE(data_size)
CASE(8,1) ! VTK_HEXAHEDRON, VTK_VERTEX
  nVTKPoints=nNodes
  nVTKCells=nElems
  NVisu = 1
  NVisuValue = 0
  NVisuCoords = 0
CASE(4)   ! VTK_QUAD
  IF(nSurfSample.GT.1) THEN
    NPlot_p1_2 = (nSurfSample+1)**2
    nVTKPoints = nElems * NPlot_p1_2
    nVTKCells = nElems * nSurfSample**2
    NVisu = nSurfSample
    NVisuValue = nSurfSample - 1
    NVisuCoords = nSurfSample
  ELSE
    nVTKPoints=nNodes
    nVTKCells=nElems
    NVisu = 1
    NVisuValue = 0
    NVisuCoords = 0
  END IF
CASE DEFAULT
  CALL abort(__STAMP__,'Wrong data size given to WriteDataToVTK_PICLas routine!')
END SELECT

ALLOCATE(Vertex(data_size,nVTKCells))
Vertex = 0

! Check if no output is present and print info that no file will be created
IF(nElems.LT.1)THEN
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')"   NOT CREATING OUTPUT FILE "
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') '['//TRIM(FileString)//'] because the dataset is empty...'
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"RETURN"
  RETURN
ELSE
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')"   WRITE 3D DATA TO VTX XML BINARY (VTU) FILE "
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') '['//TRIM(FileString)//'] ...'
END IF

! Prepare output of vector variables as a vector variable suitable for VisIt and Paraview
IF(.NOT.ALLOCATED(VarNameCombine))    ALLOCATE(VarNameCombine   (nVar))
IF(.NOT.ALLOCATED(VarNameCombineLen)) ALLOCATE(VarNameCombineLen(nVar))
VarNameCombine = 0
DO iVar=2,nVar
  ! Get the length of the variable name
  iLen = LEN(TRIM(VarNameVisu(iVar)))
  ! Save the strings in temporary variables
  tmp = VarNameVisu(iVar)
  tmp2 = VarNameVisu(iVar-1)
  ! Compare the strings, while omitting the last two characters to find identify VeloX/Y/Z and so on as vectors
  IF (TRIM(tmp(:iLen-2)) .EQ. TRIM(tmp2(:iLen-2))) THEN
    ! Although the translational temperature and the pressure tensor are given in X/Y/Z or XY/XZ/YZ they are not vectors
    ! (VisIt/Paraview would produce a magnitude variable)
    IF(INDEX(tmp(:iLen-1),'TempTrans').EQ.0 .AND. INDEX(tmp(:iLen-2),'PressTens').EQ.0) THEN
      ! If it is the first occurrence, start counting
      IF (VarNameCombine(iVar-1) .EQ. 0) VarNameCombine(iVar-1) = 1
      VarNameCombine(iVar) = VarNameCombine(iVar-1) + 1
    END IF
  END IF
END DO
VarNameCombineLen = 0
VarNameCombineLen(nVar) = VarNameCombine(nVar)
DO iVar=nVar-1,1,-1
  IF (VarNameCombine(iVar).GT.0) THEN
    VarNameCombineLen(iVar) = MAX(VarNameCombine(iVar), VarNameCombineLen(iVar+1))
  END IF
END DO

! Line feed character
lf = char(10)

! Write file
OPEN(UNIT=ivtk,FILE=TRIM(FileString),ACCESS='STREAM')
! Write header
Buffer='<?xml version="1.0"?>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf;WRITE(ivtk) TRIM(Buffer)

Buffer='  <UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
WRITE(TempStr1,'(I16)')nVTKPoints
WRITE(TempStr2,'(I16)')nVTKCells
Buffer='    <Piece NumberOfPoints="'//TRIM(ADJUSTL(TempStr1))//&
'" NumberOfCells="'//TRIM(ADJUSTL(TempStr2))//'">'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify point data
Buffer='      <PointData> </PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify cell data
Buffer='      <CellData>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=0
WRITE(StrOffset,'(I16)')Offset
IF (nVar .GT.0)THEN
  DO iVar=1,nVar
    IF(VarNameCombine(iVar).EQ.0) THEN
      Buffer='        <DataArray type="Float32" Name="'//TRIM(VarNameVisu(iVar))//&
      '" NumberOfComponents="1" format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
      Offset=Offset+INT(SIZEOF(int_size),4)+nVTKCells*INT(SIZEOF(float),4)
      WRITE(StrOffset,'(I16)')Offset
    ELSE IF (VarNameCombine(iVar).EQ.1) THEN
      str_len = LEN_TRIM(VarNameVisu(iVar))
      WRITE(components_string,'(I1)') VarNameCombineLen(iVar)
      VarNameString = VarNameVisu(iVar)(1:str_len-1)
      Buffer='        <DataArray type="Float32" Name="'//TRIM(VarNameString)//&
      '" NumberOfComponents="'//components_string//'" format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf
      WRITE(ivtk) TRIM(Buffer)
      Offset=Offset+INT(SIZEOF(int_size),4)+nVTKCells*INT(SIZEOF(float),4)*VarNameCombineLen(iVar)
      WRITE(StrOffset,'(I16)')Offset
    END IF
  END DO
END IF
Buffer='      </CellData>'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify coordinate data
Buffer='      <Points>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='        <DataArray type="Float32" Name="Coordinates" NumberOfComponents="3" format="appended"'// &
       ' offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=Offset+INT(SIZEOF(int_size),4)+3*nVTKPoints*INT(SIZEOF(float),4)
WRITE(StrOffset,'(I16)')Offset
Buffer='      </Points>'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify necessary cell data
Buffer='      <Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
! Connectivity
Buffer='        <DataArray type="Int32" Name="connectivity" format="appended"'// &
       ' offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=Offset+INT(SIZEOF(int_size),4)+data_size*nVTKCells*INT(SIZEOF(int_size),4)
WRITE(StrOffset,'(I16)')Offset
! Offsets
Buffer='        <DataArray type="Int32" Name="offsets" format="appended"'// &
       ' offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=Offset+INT(SIZEOF(int_size),4)+nVTKCells*INT(SIZEOF(int_size),4)
WRITE(StrOffset,'(I16)')Offset
! Elem types
Buffer='        <DataArray type="Int32" Name="types" format="appended"'// &
       ' offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='      </Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='    </Piece>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='  </UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
! Prepare append section
Buffer='  <AppendedData encoding="raw">'//lf;WRITE(ivtk) TRIM(Buffer)
! Write leading data underscore
Buffer='_';WRITE(ivtk) TRIM(Buffer)

! Write binary raw data into append section
! cell data
nBytes = nVTKCells*INT(SIZEOF(FLOAT),4)
DO iVal=1,nVar
  IF (VarNameCombine(iVal).EQ.0) THEN
    WRITE(ivtk) nBytes,REAL(Array(iVal,0:NVisuValue,0:NVisuValue*(MERGE(1,0,dim.GT.1)),0:NVisuValue*(MERGE(1,0,dim.GT.2)),1:nElems),4)
  ELSEIF(VarNameCombine(iVal).EQ.1) THEN
    WRITE(ivtk) nBytes*VarNameCombineLen(iVal),REAL(Array(iVal:iVal+VarNameCombineLen(iVal)-1,0:NVisuValue,0:NVisuValue*(MERGE(1,0,dim.GT.1)),0:NVisuValue*(MERGE(1,0,dim.GT.2)),1:nElems),4)
  ENDIF
END DO
! Points
nBytes = 3*nVTKPoints*INT(SIZEOF(FLOAT),4)
WRITE(ivtk) nBytes
WRITE(ivtk) REAL(Coords(1:3,0:NVisuCoords,0:NVisuCoords*(MERGE(1,0,dim.GT.1)),0:NVisuCoords*(MERGE(1,0,dim.GT.2)),1:nNodes),4)
! Connectivity
IF(nSurfSample.GT.1) THEN
  CALL CreateConnectivity(NVisu,nElems,nodeids,dim,DGFV=0)
ELSE
  DO iElem=1,nVTKCells
    DO iNode=1,data_size
      Vertex(iNode,iElem) = ConnectInfo(iNode,iElem)-1
    END DO
  END DO
END IF
nBytes = data_size*nVTKCells*INT(SIZEOF(int_size),4)
WRITE(ivtk) nBytes
IF(nSurfSample.GT.1) THEN
  WRITE(ivtk) nodeids
ELSE
  WRITE(ivtk) Vertex(:,:)
END IF
! Offset
nBytes = nVTKCells*INT(SIZEOF(int_size),4)
WRITE(ivtk) nBytes
WRITE(ivtk) (Offset,Offset=data_size,data_size*nVTKCells,data_size)
! Elem type
SELECT CASE(data_size)
CASE(8)
  ElemType = 12 ! VTK_HEXAHEDRON
CASE(4)
  ElemType = 9  ! VTK_QUAD
CASE(1)
  ElemType = 2  ! VTK_VERTEX
CASE DEFAULT
  CALL abort(__STAMP__,'Wrong data size given to WriteDataToVTK_PICLas routine!')
END SELECT
WRITE(ivtk) nBytes
WRITE(ivtk) (ElemType,iElem=1,nVTKCells)
! Write footer
Buffer=lf//'  </AppendedData>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
CLOSE(ivtk)
GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, ' DONE!', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)

SDEALLOCATE(VarNameCombine)
SDEALLOCATE(VarNameCombineLen)
SDEALLOCATE(Vertex)
SDEALLOCATE(nodeids)

END SUBROUTINE WriteDataToVTK_PICLas


!===================================================================================================================================
!> Convert DVM solution to a cell-based VTK format
!===================================================================================================================================
SUBROUTINE ConvertDVMSolution(InputStateFile,NVisu,NodeTypeVisuOut,ArrayName)
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars          ,ONLY: ProjectName
USE MOD_IO_HDF5               ,ONLY: HSize
USE MOD_HDF5_Input            ,ONLY: OpenDataFile,CloseDataFile,ReadAttribute,File_ID,ReadArray,GetDataSize,GetDataProps
USE MOD_Mesh_Vars             ,ONLY: NodeCoords,nElems,offsetElem,NGeo
USE MOD_Interpolation_Vars    ,ONLY: NodeTypeVisu
USE MOD_Interpolation         ,ONLY: GetVandermonde
USE MOD_ChangeBasis           ,ONLY: ChangeBasis3D
USE MOD_VTK                   ,ONLY: WriteDataToVTK
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)   :: InputStateFile,NodeTypeVisuOut
CHARACTER(LEN=*),INTENT(IN)     :: ArrayName
INTEGER,INTENT(IN)              :: NVisu
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)              :: FileString,NodeType_State
REAL                            :: OutputTime
INTEGER                         :: nDims,nVar,N_State,nElems_State,iElem
CHARACTER(LEN=255),ALLOCATABLE  :: VarNames(:)
REAL,ALLOCATABLE,TARGET         :: U(:,:,:,:,:)
REAL,ALLOCATABLE,TARGET         :: U_Visu(:,:,:,:,:)                 !< Solution on visualization nodes
REAL,ALLOCATABLE                :: Coords_NVisu(:,:,:,:,:)           !< Coordinates of visualization nodes
REAL,ALLOCATABLE                :: Vdm_EQNgeo_NVisu(:,:)             !< Vandermonde from equidistant mesh to visualization nodes
REAL,ALLOCATABLE                :: Vdm_N_NVisu(:,:)                  !< Vandermonde from state to visualization nodes
!===================================================================================================================================
SDEALLOCATE(Vdm_EQNgeo_NVisu)
ALLOCATE(Vdm_EQNgeo_NVisu(0:Ngeo,0:NVisu))
CALL GetVandermonde(Ngeo,NodeTypeVisu,NVisu,NodeTypeVisuOut,Vdm_EQNgeo_NVisu,modal=.FALSE.)
SDEALLOCATE(Coords_NVisu)
ALLOCATE(Coords_NVisu(3,0:NVisu,0:NVisu,0:NVisu,nElems))

! Convert coordinates to visu grid
DO iElem = 1,nElems
  CALL ChangeBasis3D(3,NGeo,NVisu,Vdm_EQNgeo_NVisu,NodeCoords(:,:,:,:,iElem),Coords_NVisu(:,:,:,:,iElem))
END DO

! Read in solution
CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)
OutputTime = 0. ! default
CALL ReadAttribute(File_ID,'Time',1,RealScalar=OutputTime)
CALL GetDataSize(File_ID,TRIM(ArrayName),nDims,HSize)
IF (nDims.NE.5) CALL abort(__STAMP__,'Wrong number of dimensions in state file!')
CALL GetDataProps(TRIM(ArrayName),nVar,N_State,nElems_State,NodeType_State)
IF (nElems.NE.nElems_State) CALL abort(__STAMP__,'Number of elements in state file and mesh file do not match!')

IF (nVar.GT.0) THEN
  ALLOCATE(VarNames(1:nVar))
  CALL ReadAttribute(File_ID,'VarNames',nVar,StrArray=VarNames(1:nVar))

  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (&
        nVar       => INT(nVar,IK)      ,&
        N_State    => INT(N_State, IK)  ,&
        offsetElem => INT(offsetElem,IK),&
        nElems     => INT(nElems,IK)    )
    SDEALLOCATE(U)
    ALLOCATE(U(1:nVar,0:N_State,0:N_State,0:N_State,1:nElems))
    CALL ReadArray(TRIM(ArrayName),5,(/nVar,N_State+1_IK,N_State+1_IK,N_State+1_IK, nElems/),offsetElem,5, &
    RealArray=U(1:nVar,0:N_State,0:N_State,0:N_State,1:nElems))
  END ASSOCIATE

  SDEALLOCATE(Vdm_N_NVisu)
  ALLOCATE(Vdm_N_NVisu(0:N_State,0:NVisu))
  CALL GetVandermonde(N_State,NodeType_State,NVisu,NodeTypeVisuOut,Vdm_N_NVisu,modal=.FALSE.)

  SDEALLOCATE(U_Visu)
  ALLOCATE(U_Visu(nVar,0:NVisu,0:NVisu,0:NVisu,nElems))

  ! Write solution to vtk
  FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution_DVM',OutputTime))//'.vtu'
  ! Interpolate solution to visu grid
  DO iElem = 1,nElems
    CALL ChangeBasis3D(nVar,N_State,NVisu,Vdm_N_NVisu,U(:,:,:,:,iElem),U_Visu(:,:,:,:,iElem))
  END DO
  ! Output to VTK
  CALL WriteDataToVTK(nVar,NVisu,nElems,VarNames,Coords_NVisu,U_Visu,TRIM(FileString),dim=3,DGFV=0)
END IF

SDEALLOCATE(VarNames)
SDEALLOCATE(Vdm_EQNgeo_NVisu)
SDEALLOCATE(Coords_NVisu)
SDEALLOCATE(U)
SDEALLOCATE(Vdm_N_NVisu)
SDEALLOCATE(U_Visu)

CALL CloseDataFile()

END SUBROUTINE ConvertDVMSolution


!===================================================================================================================================
!> Convert the output of the field solver to a VTK output format
!===================================================================================================================================
SUBROUTINE ConvertDGSolution(InputStateFile,NVisu,NodeTypeVisuOut,OutputName,DGSolutionDataset,iErrorReturn)
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars          ,ONLY: ProjectName
USE MOD_HDF5_Input            ,ONLY: OpenDataFile,ReadAttribute,File_ID,ReadArray,GetDataSize,GetDataProps,CloseDataFile
USE MOD_HDF5_Input            ,ONLY: DatasetExists
USE MOD_Mesh_Vars             ,ONLY: NodeCoords,nElems,offsetElem,NGeo
USE MOD_Interpolation_Vars    ,ONLY: NodeTypeVisu
USE MOD_Interpolation         ,ONLY: GetVandermonde
USE MOD_ChangeBasis           ,ONLY: ChangeBasis3D
USE MOD_VTK                   ,ONLY: WriteDataToVTK
USE MOD_IO_HDF5               ,ONLY: HSize
USE MOD_piclas2vtk_Vars       ,ONLY: ElemLocal, Nloc_Visu, PointToCellSwitch, NVisuAdd
USE MOD_ReadInTools           ,ONLY: PrintOption
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN) :: InputStateFile,NodeTypeVisuOut,DGSolutionDataset
CHARACTER(LEN=*),INTENT(IN)   :: OutputName
INTEGER,INTENT(IN)            :: NVisu
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)           :: iErrorReturn
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iElem, iDG, nVar_State, N_State, nElems_State, nVar_Solution, nDims, iField, nFields, Suffix
INTEGER                         :: nDimsOffset, nVar_Source,nVar_TD
INTEGER                         :: iVar, nVarAdd, Nloc, NlocMax, NlocMin, NlocOut, nDOF, iDOF, k, l, m
CHARACTER(LEN=255)              :: MeshFile, NodeType_State, FileString_DG, StrVarNamesTemp(4),StrVarNamesTemp3(3),StrVarNamesTemp4
CHARACTER(LEN=255),ALLOCATABLE  :: StrVarNames(:), StrVarNamesTemp2(:)
REAL                            :: OutputTime
REAL,ALLOCATABLE                :: U2(:,:,:,:,:,:)                   !< Solution from state file with additional dimension, rank=6
REAL,ALLOCATABLE                :: U(:,:,:,:,:)                      !< Solution from state file, rank=5
REAL,ALLOCATABLE                :: U_N_2D(:,:)                       !< Solution from state file, rank=2
REAL,ALLOCATABLE                :: U_N_3D(:,:,:)                     !< Solution from state file, rank=3
REAL,ALLOCATABLE,TARGET         :: U_Visu(:,:,:,:,:)                 !< Solution on visualization nodes
REAL,POINTER                    :: U_Visu_p(:,:,:,:,:)               !< Solution on visualization nodes
REAL,ALLOCATABLE                :: Coords_NVisu(:,:,:,:,:)           !< Coordinates of visualization nodes
REAL,ALLOCATABLE,TARGET         :: Coords_DG(:,:,:,:,:)
REAL,POINTER                    :: Coords_DG_p(:,:,:,:,:)
REAL,ALLOCATABLE                :: Vdm_EQNgeo_NVisu(:,:)             !< Vandermonde from equidistant mesh to visualization nodes
REAL,ALLOCATABLE                :: Vdm_N_NVisu(:,:)                  !< Vandermonde from state to visualization nodes
REAL,ALLOCATABLE                :: ElemData(:,:)                     !< Array for temporary read-in of ElemData container
INTEGER,ALLOCATABLE             :: Nloc_HDF5(:)                      !< Array for temporary read-in of Nloc container
LOGICAL                         :: DGSourceExists,DGTimeDerivativeExists,TimeExists,DGSourceExtExists,DMDMode,DGSolutionDatasetExists
LOGICAL                         :: ElemDataExists, NlocFound
CHARACTER(LEN=16)               :: hilf
CHARACTER(LEN=255)              :: DMDFields(1:16), Dataset, NodeType
CHARACTER(LEN=255),ALLOCATABLE  :: VarNamesAdd(:)
! p-Adaption
TYPE tNGeo
  REAL,ALLOCATABLE              :: Vdm_EQNgeo_NVisu(:,:)        !< Vandermonde from equidistant mesh to visualization nodes
END TYPE tNGeo

TYPE tNVisu
  REAL,ALLOCATABLE              :: Vdm_N_NVisu(:,:)             !< Vandermonde from state to visualization nodes
END TYPE tNVisu

TYPE(tNGeo),ALLOCATABLE         :: NVisuGeo(:)                  !< Container for polynomial degree specific variables [1:NlocMax]
TYPE(tNVisu),ALLOCATABLE        :: NVisuLocal(:,:)              !< Container for polynomial degree specific variables [1:NlocMax,1:NlocMaxVisu]
!===================================================================================================================================
! 1.) Open given file to get the number of elements, the order and the name of the mesh file
CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)

! Check if data set exists, if not return
CALL DatasetExists(File_ID,TRIM(DGSolutionDataset),DGSolutionDatasetExists)
IF(.NOT.DGSolutionDatasetExists)THEN
  CALL CloseDataFile()
  iErrorReturn=1
  RETURN
ELSE
  iErrorReturn=0
END IF ! .NOT.DGSolutionDatasetExists

! Retrieve the element-local N from the ElemData array
NlocFound = .FALSE.
CALL DatasetExists(File_ID,'ElemData',ElemDataExists)
IF(ElemDataExists) THEN
  ! Get size of the ElemData array
  CALL GetDataSize(File_ID,'ElemData',nDims,HSize)
  nVarAdd=INT(HSize(1),4)
  DEALLOCATE(HSize)
  ! Read-in the variable names
  ALLOCATE(VarNamesAdd(1:nVarAdd))
  CALL ReadAttribute(File_ID,'VarNamesAdd',nVarAdd,StrArray=VarNamesAdd(1:nVarAdd))
  ! Loop over the number of variables and find Nloc (or NlocRay, in case of a RadiationVolState), exit the loop and use the last iVar
  DO iVar=1,nVarAdd
    IF(TRIM(VarNamesAdd(iVar)).EQ.'Nloc'.OR.TRIM(VarNamesAdd(iVar)).EQ.'NlocRay') THEN
      NlocFound = .TRUE.
      EXIT
    END IF
  END DO
  IF(NlocFound) THEN
    SDEALLOCATE(Nloc_HDF5)
    SDEALLOCATE(Nloc_Visu)
    ALLOCATE(Nloc_HDF5(1:nElems))
    ALLOCATE(Nloc_Visu(1:nElems))
    ALLOCATE(ElemData(1:nVarAdd,1:nElems))
    ! Associate construct for integer KIND=8 possibility
    ASSOCIATE (&
        nVarAdd     => INT(nVarAdd,IK)    ,&
        offsetElem  => INT(offsetElem,IK) ,&
        nElems      => INT(nElems,IK)     )
      CALL ReadArray('ElemData',2,(/nVarAdd, nElems/),offsetElem,2,RealArray=ElemData(1:nVarAdd,1:nElems))
    END ASSOCIATE
    Nloc_HDF5(1:nElems) = NINT(ElemData(iVar,1:nElems))
    NlocMax = MAXVAL(Nloc_HDF5)
    NlocMin = MINVAL(Nloc_HDF5)
    CALL PrintOption('NlocMin','HDF5',IntOpt=NlocMin)
    CALL PrintOption('NlocMax','HDF5',IntOpt=NlocMax)
    IF(Ngeo.GT.NlocMin) THEN
      CALL abort(__STAMP__,'ERROR: Ngeo as read-in from mesh is greater than the smallest local polynomial degree! Ngeo: ',Ngeo)
    END IF
    DEALLOCATE(ElemData)
    ! Checking whether Nloc varies
    IF(NlocMin.EQ.NlocMax) THEN
      SWRITE(*,*) '| Found element-local polynomial degree ('//TRIM(VarNamesAdd(iVar))//'), but using NVisu since Nloc is constant.'
      Nloc_Visu(1:nElems) = NVisu
      ! Required for Nloc = 1,NlocMax loops
      NlocMax = MAX(NVisu,NlocMax)
      NlocMin = MIN(NVisu,NlocMin)
    ELSE
      SWRITE(*,*) '| Found element-local polynomial degree ('//TRIM(VarNamesAdd(iVar))//'), considering it for the output of each element and adding NVisuAdd.'
      ! Increasing the cell-local polynomial degree by NVisuAdd
      Nloc_Visu(1:nElems) = Nloc_HDF5(1:nElems) + NVisuAdd
      ! NlocMax is increasing accordingly, NlocMin remains the same as it is used to create the Vandermonde
      NlocMax = NlocMax + NVisuAdd
    END IF
  END IF
  DEALLOCATE(VarNamesAdd)
END IF

! Read-in of dimensions of the field array (might have an additional dimension, i.e., rank is 6 instead of 5)
CALL GetDataSize(File_ID,TRIM(DGSolutionDataset),nDims,HSize)
! Check the number of fields in the file, if more than 5 dimensions, the 6th dimensions carries the number of fields
IF(nDims.GE.5) THEN
nFields     = MERGE(1 , INT(HSize(nDims)) , nDims.EQ.5)
nDimsOffset = MERGE(0 , 1                 , nDims.EQ.5)
ELSE
  ! p-adaption data format
  IF (nDims.EQ.2) THEN
    ! For nDims=2, assume shape [1:nVar, 1:nDOF]
    nFields     = 1
    nDimsOffset = 0
  ELSEIF (nDims.EQ.3) THEN
    ! For nDims=3, assume shape [1:nVar, 1:nDOF, 1:nFields] - the nFields e.g. could be 1:nTimePoints for time-dependent data
    nFields = INT(HSize(nDims))
    nDimsOffset = 1
  ELSE
    ! Unknown - abort
    IPWRITE(*,*) 'HSize:', HSize
    CALL abort(__STAMP__,'ERROR: Unknown array shape (HSize) for '//TRIM(DGSolutionDataset)//' encountered.')
  END IF ! nDims.EQ.2
END IF
DEALLOCATE(HSize)

! Check compatibility of ConvertPointToCellData with other features
IF(PointToCellSwitch.EQ.1) THEN
  IF(nFields.GT.1) CALL abort(__STAMP__,'ERROR: ConvertPointToCellData and nFields > 1 is not supported!')
  CALL ReadAttribute(File_ID,'NodeType',1,StrScalar=NodeType)
  IF(TRIM(NodeType).NE.'VISU_INNER') CALL abort(__STAMP__,'ERROR: ConvertPointToCellData is only implemented for point data using NodeType = VISU_INNER!')
END IF

! Default: DGSolutionDataset = 'DG_Solution'
CALL GetDataProps(TRIM(DGSolutionDataset),nVar_Solution,N_State,nElems_State,NodeType_State,nDimsOffset)
CALL PrintOption('N_State','HDF5',IntOpt=N_State)
CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)

! Check if the DG_Source container exists, and if it does save the variable names and increase the nVar_State variable
CALL DatasetExists(File_ID,'DG_Source',DGSourceExists)
IF (DGSourceExists) THEN
  CALL ReadAttribute(File_ID,'VarNamesSource',4,StrArray=StrVarNamesTemp)
  nVar_State  = nVar_Solution + 4
ELSE
  nVar_State = nVar_Solution
END IF
nVar_Source = nVar_State

! Check if the DG_TimeDerivative container exists, and if it does save the variable names and increase the nVar_State variable
CALL DatasetExists(File_ID,'DG_TimeDerivative',DGTimeDerivativeExists)
IF(DGTimeDerivativeExists)THEN
  CALL ReadAttribute(File_ID,'VarNamesTimeDerivative',3,StrArray=StrVarNamesTemp3)
  nVar_State = nVar_State + 3
END IF ! DGTimeDerivativeExists
nVar_TD = nVar_State

! Check if the DG_SourceExt container exists, and if it does save the variable names and increase the nVar_State variable
CALL DatasetExists(File_ID,'DG_SourceExt',DGSourceExtExists)
IF(DGSourceExtExists)THEN
  StrVarNamesTemp4='SurfaceChargeDensity[C/m3]'
  nVar_State = nVar_State + 1
END IF ! DGSourceExtExists

! Save the variable names for the regular DG_Solution in a temporary array
SDEALLOCATE(StrVarNamesTemp2)
ALLOCATE(StrVarNamesTemp2(nVar_Solution))
CALL ReadAttribute(File_ID,'VarNames',nVar_Solution,StrArray=StrVarNamesTemp2)

! Allocate the variable names array used for the output and copy the names from the DG_Solution and DG_Source (if it exists)
SDEALLOCATE(StrVarNames)
! Check if DMD modes are converted
IF(TRIM(DGSolutionDataset(1:MIN(LEN(TRIM(DGSolutionDataset)),5))).EQ.'Mode_')THEN
  DMDMode = .TRUE.
  nVar_State = 16
ELSE
  DMDMode = .FALSE.
END IF ! TRIM(DGSolutionDataset(1:MIN(LEN(TRIM(DGSolutionDataset)),5))).EQ.'Mode_'
ALLOCATE(StrVarNames(nVar_State))
StrVarNames(1:nVar_Solution) = StrVarNamesTemp2
IF(DGSourceExists)         StrVarNames(nVar_Solution+1:nVar_Source) = StrVarNamesTemp(1:4)
IF(DGTimeDerivativeExists) StrVarNames(nVar_Source  +1:nVar_TD)     = StrVarNamesTemp3(1:3)
IF(DGSourceExtExists)      StrVarNames(nVar_TD      +1:nVar_State)  = StrVarNamesTemp4

! Check if a file with time stamp is converted. Read the time stamp from .h5
CALL DatasetExists(File_ID,'Time',TimeExists,attrib=.TRUE.)
IF(TimeExists) CALL ReadAttribute(File_ID,'Time',1,RealScalar=OutputTime)
IF(TRIM(OutputName).EQ.'RadiationVolState') TimeExists = .FALSE.
CALL CloseDataFile()

! Check for 2D array (nVar,nDOF)
IF((nDims.GE.2).AND.(nDims.LT.5)) THEN
  IF(.NOT.NlocFound) CALL abort(__STAMP__,'ERROR: Missing Nloc array for read-in of 2D DG_Solution!')
  nDOF = SUM((Nloc_HDF5(1:nElems)+1)**3)
  IF(nDims.EQ.2) THEN
    ! Allocate local 2D array
    ALLOCATE(U_N_2D(1:nVar_State,1:nDOF))
  ELSEIF(nDims.EQ.3) THEN
    ! Allocate local 3D array
    ALLOCATE(U_N_3D(1:nVar_State,1:nDOF,1:nFields))
  ELSE
    CALL abort(__STAMP__,'ERROR: Unknown nDims >=2 and <5 encountered.')
  END IF
END IF

IF(NlocFound) THEN
  SDEALLOCATE(NVisuGeo)
  NlocMax = NlocMax + PointToCellSwitch
  ALLOCATE(NVisuGeo(1:NlocMax))
  DO Nloc = 1, NlocMax
    ALLOCATE(NVisuGeo(Nloc)%Vdm_EQNgeo_NVisu(0:Nloc,0:NGeo))
    CALL GetVandermonde(Ngeo,NodeTypeVisu,Nloc,NodeTypeVisuOut,NVisuGeo(Nloc)%Vdm_EQNgeo_NVisu,modal=.FALSE.)
  END DO
  SDEALLOCATE(ElemLocal)
  ALLOCATE(ElemLocal(1:nElems))
  ! Convert coordinates to visu grid
  DO iElem = 1 , nElems
    Nloc = Nloc_Visu(iElem) + PointToCellSwitch
    ALLOCATE(ElemLocal(iElem)%Coords_NVisu(1:3,0:Nloc,0:Nloc,0:Nloc))
    ! ALLOCATE(tempArray(3,0:Nloc,0:Nloc,0:Nloc))
    CALL ChangeBasis3D(3, NGeo, Nloc, NVisuGeo(Nloc)%Vdm_EQNgeo_NVisu(0:Nloc,0:NGeo), NodeCoords(1:3,0:NGeo,0:NGeo,0:NGeo,iElem), &
                                                                     ElemLocal(iElem)%Coords_NVisu(1:3,0:Nloc,0:Nloc,0:Nloc))
  END DO
ELSE
SDEALLOCATE(Vdm_EQNgeo_NVisu)
  ALLOCATE(Vdm_EQNgeo_NVisu(0:NVisu,0:NGeo))
CALL GetVandermonde(Ngeo,NodeTypeVisu,NVisu,NodeTypeVisuOut,Vdm_EQNgeo_NVisu,modal=.FALSE.)
SDEALLOCATE(Coords_NVisu)
ALLOCATE(Coords_NVisu(3,0:NVisu,0:NVisu,0:NVisu,nElems))
SDEALLOCATE(Coords_DG)
ALLOCATE(Coords_DG(3,0:NVisu,0:NVisu,0:NVisu,nElems))
! Convert coordinates to visu grid
DO iElem = 1,nElems
    CALL ChangeBasis3D(3, NGeo, NVisu, Vdm_EQNgeo_NVisu(0:NVisu,0:NGeo), NodeCoords(1:3,0:NGeo ,0:NGeo ,0:NGeo ,iElem), &
                                                                       Coords_NVisu(1:3,0:NVisu,0:NVisu,0:NVisu,iElem))
END DO
END IF

! Read in solution
CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nVar_State    => INT(nVar_State    , IK)      , &
      nVar_Solution => INT(nVar_Solution , IK)      , &
      offsetElem    => INT(offsetElem    , IK)      , &
      N_State       => INT(N_State       , IK)      , &
      nElems        => INT(nElems        , IK)      , &
      nDOF          => INT(nDOF          , IK)      , &
      nFields       => INT(nFields       , IK)      )

  ! Check if DMD data is to be converted or normal data
  IF(DMDMode)THEN
    !Sanity check
    IF(LEN(TRIM(DGSolutionDataset)).LT.9) CALL abort(__STAMP__,'DGSolutionDataset name is too short: '//TRIM(DGSolutionDataset))
      DMDFields=(/'ElectricFieldX_Img ',&
                  'ElectricFieldX_Real',&
                  'ElectricFieldY_Img ',&
                  'ElectricFieldY_Real',&
                  'ElectricFieldZ_Img ',&
                  'ElectricFieldZ_Real',&
                  'MagneticFieldX_Img ',&
                  'MagneticFieldX_Real',&
                  'MagneticFieldY_Img ',&
                  'MagneticFieldY_Real',&
                  'MagneticFieldZ_Img ',&
                  'MagneticFieldZ_Real',&
                  'Phi_Img            ',&
                  'Phi_Real           ',&
                  'Psi_Img            ',&
                  'Psi_Real           '/)

      SDEALLOCATE(U)
      ALLOCATE(U(nVar_State,0:N_State,0:N_State,0:N_State,nElems))
      DO iField = 1, 16
        !Dataset='Mode_'//TRIM(hilf)//'_'//DMDFields(iField)
        Dataset=TRIM(DGSolutionDataset(1:9))//DMDFields(iField)
        StrVarNames(iField) = TRIM(Dataset)
        WRITE (*,*) "Converting ... ", TRIM(Dataset)
        CALL ReadArray(TRIM(Dataset),5,(/nVar_Solution,N_State+1_IK,N_State+1_IK,N_State+1_IK,nElems/),offsetElem,5, &
            RealArray=U(iField:iField,0:N_State,0:N_State,0:N_State,1:nElems))
      END DO ! iField = 1, 16

  ELSE
    IF(nFields.EQ.1)THEN
      ! Default: DGSolutionDataset = 'DG_Solution'
      ! Check whether old (nDims=5) or new p-adaption data shape (nDim=2) is used
      IF(nDims.EQ.2) THEN
        IF(offsetElem.NE.0) CALL abort(__STAMP__,'offsetElem must be zero for the following to work (it should also be offsetDOF)')
        ! TODO: Implement offsetDOF
        CALL ReadArray(TRIM(DGSolutionDataset),2,(/nVar_Solution,nDOF/),offsetElem,2, RealArray=U_N_2D(1:nVar_Solution,1:nDOF))
        ! Read nVar_Solution+1:nVar_Source
        IF(DGSourceExists) &
          CALL ReadArray('DG_Source',2,(/nVar_Solution,nDOF/),offsetElem,2, RealArray=U_N_2D(nVar_Solution+1:nVar_Source,1:nDOF))
        ! Read nVar_Source+1:nVar_TD
        IF(DGTimeDerivativeExists) &
          CALL ReadArray('DG_TimeDerivative',2,(/nVar_Solution,nDOF/),offsetElem,2, RealArray=U_N_2D(nVar_Source+1:nVar_TD,1:nDOF))
        ! Read nVar_TD+1:nVar_State
        IF(DGSourceExtExists) &
          CALL ReadArray('DG_SourceExt',2,(/nVar_Solution,nDOF/),offsetElem,2, RealArray=U_N_2D(nVar_TD+1:nVar_State,1:nDOF))
      ELSEIF(nDims.EQ.5) THEN
      SDEALLOCATE(U)
      ALLOCATE(U(nVar_State,0:N_State,0:N_State,0:N_State,nElems))
      CALL ReadArray(TRIM(DGSolutionDataset),5,(/nVar_Solution,N_State+1_IK,N_State+1_IK,N_State+1_IK,nElems/),offsetElem,5, &
          RealArray=U(1:nVar_Solution,0:N_State,0:N_State,0:N_State,1:nElems))
      ! Read nVar_Solution+1:nVar_Source
      IF(DGSourceExists) CALL ReadArray('DG_Source',5,(/4_IK,N_State+1_IK,N_State+1_IK,N_State+1_IK,nElems/),offsetElem,5, &
          RealArray=U(nVar_Solution+1:nVar_Source,0:N_State,0:N_State,0:N_State,1:nElems))
      ! Read nVar_Source+1:nVar_TD
      IF(DGTimeDerivativeExists) CALL ReadArray('DG_TimeDerivative',5,(/3_IK,N_State+1_IK,N_State+1_IK,N_State+1_IK,nElems/),&
          offsetElem,5,RealArray=U(nVar_Source+1:nVar_TD,0:N_State,0:N_State,0:N_State,1:nElems))
      ! Read nVar_TD+1:nVar_State
      IF(DGSourceExtExists) CALL ReadArray('DG_SourceExt',5,(/1_IK,N_State+1_IK,N_State+1_IK,N_State+1_IK,nElems/),&
          offsetElem,5,RealArray=U(nVar_TD+1:nVar_State,0:N_State,0:N_State,0:N_State,1:nElems))
      END IF
    ELSE ! more than one field
      ! Check whether old (nDims=5) or new p-adaption data shape (nDim=3, 2 + 1 for nFields) is used
      IF(nDims.EQ.3) THEN
        IF(offsetElem.NE.0) CALL abort(__STAMP__,'offsetElem must be zero for the following to work (it should also be offsetDOF)')
        ! TODO: Implement offsetDOF
        CALL ReadArray(TRIM(DGSolutionDataset),3,(/nVar_Solution,nDOF,nFields/),offsetElem,2,&
                      RealArray=U_N_3D(1:nVar_Solution,1:nDOF,1:nFields))
      ELSEIF(nDims.EQ.6) THEN
      SDEALLOCATE(U2)
      ALLOCATE(U2(nVar_State,0:N_State,0:N_State,0:N_State,nElems,nFields))
      ! Default: DGSolutionDataset = 'DG_Solution'
      CALL ReadArray(TRIM(DGSolutionDataset),6,(/nVar_Solution,N_State+1_IK,N_State+1_IK,N_State+1_IK,nElems,nFields/),offsetElem,5, &
          RealArray=U2(1:nVar_Solution,0:N_State,0:N_State,0:N_State,1:nElems,1:nFields))
      ELSE
        IPWRITE(*,*) 'nDims:', nDims
        CALL abort(__STAMP__,'Unknown nDims for nFields>1 encountered')
      END IF
    END IF ! nFields.GT.1
  END IF ! TRIM(DGSolutionDataset(1:MIN(LEN(TRIM(DGSolutionDataset)),5))).EQ.'Mode_'
END ASSOCIATE

CALL CloseDataFile()

IF(NlocFound) THEN
  ! New format: p-adaption
  IF((nDims.EQ.2).OR.(nDims.EQ.3)) THEN
    SDEALLOCATE(NVisuLocal)
    ALLOCATE(NVisuLocal(NlocMin:NlocMax,NlocMin:NlocMax))
    DO Nloc = NlocMin, NlocMax
      DO NlocOut = NlocMin, NlocMax
        ALLOCATE(NVisuLocal(Nloc,NlocOut)%Vdm_N_NVisu(0:NlocOut,0:Nloc))
        IF(Nloc.GT.NlocOut) THEN
          CALL GetVandermonde(Nloc,NodeType_State,NlocOut,NodeTypeVisuOut,NVisuLocal(Nloc,NlocOut)%Vdm_N_NVisu,modal=.TRUE.)
        ELSE
          CALL GetVandermonde(Nloc,NodeType_State,NlocOut,NodeTypeVisuOut,NVisuLocal(Nloc,NlocOut)%Vdm_N_NVisu,modal=.FALSE.)
        END IF
      END DO
    END DO
  ELSE IF(nDims.EQ.5) THEN
    ! TODO: Is this even possible? Aren't all outputs with the Nloc container utilizing the new nDims = 2 format?
    SDEALLOCATE(NVisuLocal)
    ALLOCATE(NVisuLocal(N_State:NlocMax,N_State:NlocMax))
    DO Nloc = N_State, NlocMax
      ALLOCATE(NVisuLocal(Nloc,N_State)%Vdm_N_NVisu(0:Nloc,0:N_State))
      CALL GetVandermonde(N_State,NodeType_State,Nloc,NodeTypeVisuOut,NVisuLocal(Nloc,N_State)%Vdm_N_NVisu,modal=.FALSE.)
    END DO
  END IF
  ! Allocate element-local solution vector for visualization
  DO iElem = 1,nElems
    Nloc = Nloc_Visu(iElem)
    ALLOCATE(ElemLocal(iElem)%U_Visu(nVar_State,0:Nloc,0:Nloc,0:Nloc))
  END DO
  ! Write solution into element-local array
  IF(nDims.EQ.2) THEN
    iDOF = 0
    DO iElem = 1,nElems
      Nloc = Nloc_HDF5(iElem)
      ALLOCATE(ElemLocal(iElem)%U(nVar_State,0:Nloc,0:Nloc,0:Nloc))
      DO m=0,Nloc
        DO l=0,Nloc
          DO k=0,Nloc
            iDOF = iDOF + 1
            ! TODO: Treatment required when NodeTypeVisuOut differs from NodeType_State
            ElemLocal(iElem)%U(1:nVar_State,k,l,m) = U_N_2D(1:nVar_State,iDOF)
          END DO ! k
        END DO ! l
      END DO ! m
    END DO
  ELSEIF(nDims.EQ.3) THEN
    iDOF = 0
    DO iElem = 1,nElems
      Nloc = Nloc_HDF5(iElem)
      ALLOCATE(ElemLocal(iElem)%U2(nVar_State,0:Nloc,0:Nloc,0:Nloc,1:nFields))
      DO m=0,Nloc
        DO l=0,Nloc
          DO k=0,Nloc
            iDOF = iDOF + 1
            ! TODO: Treatment required when NodeTypeVisuOut differs from NodeType_State
            ElemLocal(iElem)%U2(1:nVar_State,k,l,m,1:nFields) = U_N_3D(1:nVar_State,iDOF,1:nFields)
          END DO ! k
        END DO ! l
      END DO ! m
    END DO
  END IF
ELSE
  ! Old format: constant polynomial degree
SDEALLOCATE(Vdm_N_NVisu)
ALLOCATE(Vdm_N_NVisu(0:N_State,0:NVisu))
CALL GetVandermonde(N_State,NodeType_State,NVisu,NodeTypeVisuOut,Vdm_N_NVisu,modal=.FALSE.)

SDEALLOCATE(U_Visu)
ALLOCATE(U_Visu(nVar_State,0:NVisu,0:NVisu,0:NVisu,nElems))

! Set DG coords
iDG = 0
DO iElem = 1,nElems
  iDG = iDG + 1
  Coords_DG(:,:,:,:,iDG) = Coords_NVisu(:,:,:,:,iElem)
END DO
Coords_DG_p => Coords_DG(:,:,:,:,1:iDG)
U_Visu_p    => U_Visu(:,:,:,:,1:iDG)
END IF

! Output multiples files, if the DGSolution file contains more than one field (e.g. multiple different times for the same field)
DO iField = 1, nFields

  ! Write solution to vtk
  IF(DMDMode)THEN
    FileString_DG=TRIM(ProjectName)//'_'//TRIM(OutputName)//'_'//TRIM(DGSolutionDataset(1:8))//'.vtu'
  ELSEIF(OutputName.EQ.'State')THEN
    FileString_DG=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution',OutputTime))//'.vtu'
  ELSE
    IF(TimeExists)THEN
      FileString_DG=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(OutputName),OutputTime))//'.vtu'
    ELSE
      FileString_DG=TRIM(ProjectName)//'_'//TRIM(OutputName)//'.vtu'
    END IF ! TimeExists
  END IF

  IF(nFields.GT.1)THEN
    ! Append 001, 002, etc. to file name (before the .vtk suffix)
    Suffix = INDEX(FileString_DG,'.vtu')
    WRITE(UNIT=hilf,FMT='(I3.3)') iField
    hilf = TRIM(hilf)//'.vtu'
    FileString_DG(Suffix:Suffix-1+LEN(TRIM(hilf)))=TRIM(hilf)
  END IF ! nFields.GT.1

  ! Interpolate solution to visu grid
  IF(NlocFound) THEN
    ! New format: p-adaption
    DO iElem = 1,nElems
      Nloc = Nloc_HDF5(iElem)
      NlocOut = Nloc_Visu(iElem)
      IF(nFields.EQ.1)THEN
        IF(PointToCellSwitch.EQ.0) THEN
          IF(nDims.EQ.2) THEN
            ! New output format: might require interpolation from Nloc to NlocOut
            CALL ChangeBasis3D(nVar_State, Nloc, NlocOut, NVisuLocal(Nloc,NlocOut)%Vdm_N_NVisu(0:NlocOut,0:Nloc), &
                               ElemLocal(iElem)%U(1:nVar_State,0:Nloc,0:Nloc,0:Nloc)                            , &
                               ElemLocal(iElem)%U_Visu(1:nVar_State,0:NlocOut,0:NlocOut,0:NlocOut))
          ELSE IF(nDims.EQ.5) THEN
            ! Old output format on Nmax requires interpolation to Nloc  TODO: is this case even relevant?
            CALL ChangeBasis3D(nVar_State, N_State, NlocOut, NVisuLocal(NlocOut,N_State)%Vdm_N_NVisu(0:NlocOut,0:N_State), &
                               U( 1:nVar_State,0:N_State,0:N_State,0:N_State,iElem)                                      , &
                               ElemLocal(iElem)%U_Visu(1:nVar_State,0:NlocOut,0:NlocOut,0:NlocOut))
          END IF
        ELSE IF(PointToCellSwitch.EQ.1) THEN
          IF(nDims.EQ.2) THEN
            ! New output format using the cell-local polynomial degree
            ElemLocal(iElem)%U_Visu(1:nVar_State,0:NlocOut,0:NlocOut,0:NlocOut) = ElemLocal(iElem)%U(1:nVar_State,0:NlocOut,0:NlocOut,0:NlocOut)
          ELSEIF(nDims.EQ.5) THEN
            ! Old output format on Nmax requires interpolation to Nloc TODO: is this case even relevant?
            IF(NlocOut.NE.N_State) THEN
              CALL ChangeBasis3D(nVar_State, N_State, NlocOut, NVisuLocal(NlocOut,N_State)%Vdm_N_NVisu(0:NlocOut,0:N_State), &
                                 U( 1:nVar_State,0:N_State,0:N_State,0:N_State,iElem)                                      , &
                                 ElemLocal(iElem)%U_Visu(1:nVar_State,0:NlocOut,0:NlocOut,0:NlocOut))
            ELSE
              ElemLocal(iElem)%U_Visu(1:nVar_State,0:Nloc,0:Nloc,0:Nloc) = U( 1:nVar_State,0:N_State,0:N_State,0:N_State,iElem)
            END IF
          END IF
        END IF
      ELSE ! more than one field
        IF(PointToCellSwitch.EQ.0) THEN
          IF(nDims.EQ.3) THEN
            ! New output format: might require interpolation from Nloc to NlocOut
            CALL ChangeBasis3D(nVar_State, Nloc, NlocOut, NVisuLocal(Nloc,NlocOut)%Vdm_N_NVisu(0:NlocOut,0:Nloc), &
                               ElemLocal(iElem)%U2(1:nVar_State,0:Nloc,0:Nloc,0:Nloc,iField)                            , &
                               ElemLocal(iElem)%U_Visu(1:nVar_State,0:NlocOut,0:NlocOut,0:NlocOut))
          ELSE IF(nDims.EQ.6) THEN
            ! Old output format on Nmax requires interpolation to Nloc  TODO: is this case even relevant?
            CALL ChangeBasis3D(nVar_State, N_State, NlocOut, NVisuLocal(NlocOut,N_State)%Vdm_N_NVisu(0:NlocOut,0:N_State), &
                               U2( 1:nVar_State,0:N_State,0:N_State,0:N_State,iElem,iField)                                      , &
                               ElemLocal(iElem)%U_Visu(1:nVar_State,0:NlocOut,0:NlocOut,0:NlocOut))
          END IF
        ELSE IF(PointToCellSwitch.EQ.1) THEN
          IF(nDims.EQ.3) THEN
            ! New output format using the cell-local polynomial degree
            ElemLocal(iElem)%U_Visu(1:nVar_State,0:NlocOut,0:NlocOut,0:NlocOut) = ElemLocal(iElem)%U2(1:nVar_State,0:NlocOut,0:NlocOut,0:NlocOut,iField)
          ELSEIF(nDims.EQ.6) THEN
            ! Old output format on Nmax requires interpolation to Nloc TODO: is this case even relevant?
            IF(NlocOut.NE.N_State) THEN
              CALL ChangeBasis3D(nVar_State, N_State, NlocOut, NVisuLocal(NlocOut,N_State)%Vdm_N_NVisu(0:NlocOut,0:N_State), &
                                 U2( 1:nVar_State,0:N_State,0:N_State,0:N_State,iElem,iField)                                      , &
                                 ElemLocal(iElem)%U_Visu(1:nVar_State,0:NlocOut,0:NlocOut,0:NlocOut))
            ELSE
              ElemLocal(iElem)%U_Visu(1:nVar_State,0:Nloc,0:Nloc,0:Nloc) = U2( 1:nVar_State,0:N_State,0:N_State,0:N_State,iElem,iField)
            END IF
          END IF
        END IF
      END IF ! nFields.GT.1
    END DO
  ELSE
    ! Old format: constant polynomial degree
  iDG = 0
  DO iElem = 1,nElems
    iDG = iDG + 1
    IF(nFields.EQ.1)THEN
        CALL ChangeBasis3D(nVar_State, N_State, NVisu, Vdm_N_NVisu             , &
                           U( 1:nVar_State,0:N_State,0:N_State,0:N_State,iElem), &
                           U_Visu(1:nVar_State,0:NVisu,0:NVisu,0:NVisu,iDG))
    ELSE ! more than one field
        CALL ChangeBasis3D(nVar_State, N_State, NVisu, Vdm_N_NVisu                    , &
                           U2(1:nVar_State,0:N_State,0:N_State,0:N_State,iElem,iField), &
                           U_Visu(1:nVar_State,0:NVisu,0:NVisu,0:NVisu,iDG))
    END IF ! nFields.GT.1
  END DO
  END IF

  ! Output to VTK
  IF(NlocFound) THEN
    ! New format: p-adaption
    CALL WriteDataToVTKpAdaption(nVar_State,StrVarNames,TRIM(FileString_DG))
  ELSE
    ! Old format: constant polynomial degree
  CALL WriteDataToVTK(nVar_State,NVisu,iDG,StrVarNames,Coords_DG_p,U_Visu_p,TRIM(FileString_DG),dim=3,DGFV=0)
  END IF

END DO ! iField = 1, nFields

SDEALLOCATE(StrVarNamesTemp2)
SDEALLOCATE(StrVarNames)
SDEALLOCATE(Vdm_EQNgeo_NVisu)
SDEALLOCATE(Coords_NVisu)
SDEALLOCATE(Coords_DG)
SDEALLOCATE(U)
SDEALLOCATE(U_N_2D)
SDEALLOCATE(U2)
SDEALLOCATE(Vdm_N_NVisu)
SDEALLOCATE(U_Visu)
IF(NlocFound) THEN
  DO Nloc = NlocMin, NlocMax
    DO NlocOut = NlocMin, NlocMax
      SDEALLOCATE(NVisuLocal(Nloc,NlocOut)%Vdm_N_NVisu)
    END DO
  END DO
  DO Nloc = 1, NlocMax
    SDEALLOCATE(NVisuGeo(Nloc)%Vdm_EQNgeo_NVisu)
  END DO
  SDEALLOCATE(NVisuGeo)
  SDEALLOCATE(NVisuLocal)
  DO iElem = 1,nElems
    SDEALLOCATE(ElemLocal(iElem)%Coords_NVisu)
    SDEALLOCATE(ElemLocal(iElem)%U)
    SDEALLOCATE(ElemLocal(iElem)%U_Visu)
  END DO
  SDEALLOCATE(ElemLocal)
  SDEALLOCATE(Nloc_HDF5)
  SDEALLOCATE(Nloc_Visu)
END IF

END SUBROUTINE ConvertDGSolution


!===================================================================================================================================
!> Convert particle data (output of each particle) to the VTK format
!===================================================================================================================================
SUBROUTINE ConvertPartData(InputStateFile)
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars ,ONLY: ProjectName
USE MOD_IO_HDF5      ,ONLY: HSize
USE MOD_HDF5_Input   ,ONLY: OpenDataFile,ReadAttribute,File_ID,ReadArray,GetDataSize,CloseDataFile
USE MOD_HDF5_Input   ,ONLY: DatasetExists
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)   :: InputStateFile
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: nParts, nPartsVar, iPart, nDims
INTEGER,ALLOCATABLE             :: ConnectInfo(:,:)
CHARACTER(LEN=255),ALLOCATABLE  :: VarNamesParticle(:), tmpArray(:)
CHARACTER(LEN=255)              :: FileString
REAL, ALLOCATABLE               :: PartData(:,:)
REAL                            :: OutputTime, FileVersionHDFReal
INTEGER                         :: FileVersionHDF5Int
LOGICAL                         :: FileVersionExists
!===================================================================================================================================

CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)
CALL ReadAttribute(File_ID,'Time',1,RealScalar=OutputTime)

! check file version
CALL DatasetExists(File_ID,'File_Version',FileVersionExists,attrib=.TRUE.)
IF (FileVersionExists) THEN
  CALL ReadAttribute(File_ID,'File_Version',1,RealScalar=FileVersionHDFReal)

  IF(FileVersionHDFReal.LT.1.5)THEN
    SWRITE(UNIT_StdOut,'(A)')' '
    SWRITE(UNIT_StdOut,'(A)')' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
    SWRITE(UNIT_StdOut,'(A)')' '
    SWRITE(UNIT_StdOut,'(A)')' Restart file is too old! "File_Version" in restart file < 1.5!'
    SWRITE(UNIT_StdOut,'(A)')' The format used in the restart file is not compatible with this version of PICLas.'
    SWRITE(UNIT_StdOut,'(A)')' Among others, the particle format (PartData) has changed.'
    SWRITE(UNIT_StdOut,'(A)')' Run python script '
    SWRITE(UNIT_StdOut,'(A)')' '
    SWRITE(UNIT_StdOut,'(A)')'     python  ./tools/flip_PartState/flip_PartState.py  --help'
    SWRITE(UNIT_StdOut,'(A)')' '
    SWRITE(UNIT_StdOut,'(A)')' for info regarding the usage and run the script against the restart file, e.g., '
    SWRITE(UNIT_StdOut,'(A)')' '
    SWRITE(UNIT_StdOut,'(A)')'     python  ./tools/flip_PartState/flip_PartState.py  ProjectName_State_000.0000xxxxxx.h5'
    SWRITE(UNIT_StdOut,'(A)')' '
    SWRITE(UNIT_StdOut,'(A)')' to update the format and file version number.'
    SWRITE(UNIT_StdOut,'(A)')' Note that the format can be changed back to the old one by running the script a second time.'
    SWRITE(UNIT_StdOut,'(A)')' '
    SWRITE(UNIT_StdOut,'(A)')' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
    CALL abort(&
    __STAMP__&
    ,'Error in InitRestart(): "File_Version" in restart file < 1.5. See error message above to fix. File version in restart file =',&
    RealInfoOpt=FileVersionHDFReal)
  END IF ! FileVersionHDFReal.LT.1.5
ELSE
  CALL DatasetExists(File_ID,'Piclas_VersionInt',FileVersionExists,attrib=.TRUE.)
  IF (FileVersionExists) THEN
    CALL ReadAttribute(File_ID,'Piclas_VersionInt',1,IntScalar=FileVersionHDF5Int)
  ELSE
    CALL abort(__STAMP__,'Error in InitRestart(): Attribute "Piclas_VersionInt" does not exist!')
  END IF
END IF

! Read-in of dimensions of the particle array (1: Number of particles, 2: Number of variables)
CALL GetDataSize(File_ID,'PartData',nDims,HSize)
! First 3 entries are the particle positions, which are used as the coordinates for the output and not included as a variable
nPartsVar=INT(HSize(1),4)-3
nParts=INT(HSize(2),4)
DEALLOCATE(HSize)
! Allocating the array for the variables and a temporary array since ParticlePositionX,Y,Z are included in the read-in
ALLOCATE(VarNamesParticle(nPartsVar),tmpArray(nPartsVar+3))
CALL ReadAttribute(File_ID,'VarNamesParticles',nPartsVar+3,StrArray=tmpArray)
VarNamesParticle(1:nPartsVar)=tmpArray(4:nPartsVar+3)

! Allocate array (also when no particles exist to prevent to tool from crashing)
ALLOCATE(PartData(1:nPartsVar+3,1:nParts))
PartData = 0.
SDEALLOCATE(ConnectInfo)
ALLOCATE(ConnectInfo(1,1:nParts))
ConnectInfo = 0

ASSOCIATE(nParts    => INT(nParts,IK),  &
          nPartsVar => INT(nPartsVar,IK))
CALL ReadArray('PartData',2,(/nPartsVar+3_IK,nParts/),0_IK,1,RealArray=PartData)
END ASSOCIATE

DO iPart=1,nParts
  ConnectInfo(1,iPart)=iPart
END DO

FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_visuPart',OutputTime))//'.vtu'
 CALL WriteDataToVTK_PICLas(1,1,FileString,nPartsVar,VarNamesParticle,nParts,PartData(1:3,1:nParts),nParts,&
                             PartData(4:nPartsVar+3,1:nParts),ConnectInfo(1,1:nParts))

SDEALLOCATE(VarNamesParticle)
SDEALLOCATE(tmpArray)
SDEALLOCATE(PartData)
SDEALLOCATE(ConnectInfo)

CALL CloseDataFile()

END SUBROUTINE ConvertPartData


!===================================================================================================================================
!> Convert element/volume data (single value per cell, e.g. DSMC/BGK results) to a VTK format
!===================================================================================================================================
SUBROUTINE ConvertElemData(InputStateFile, ArrayName, VarName, iArgs)
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars    ,ONLY: ProjectName
USE MOD_IO_HDF5         ,ONLY: HSize
USE MOD_HDF5_Input      ,ONLY: OpenDataFile,CloseDataFile,ReadAttribute,File_ID,ReadArray,GetDataSize
USE MOD_Mesh_Vars       ,ONLY: nElems, offsetElem
USE MOD_piclas2vtk_Vars ,ONLY: nUniqueNodes, NodeCoords_Connect, ElemUniqueNodeID
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)   :: InputStateFile
CHARACTER(LEN=*),INTENT(IN)     :: ArrayName
CHARACTER(LEN=*),INTENT(IN)     :: VarName
INTEGER, INTENT(IN)             :: iArgs
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)              :: FileString, File_Type, File_Num
REAL                            :: OutputTime
INTEGER                         :: nDims,nVarAdd
CHARACTER(LEN=255),ALLOCATABLE  :: VarNamesAdd(:)
REAL,ALLOCATABLE                :: ElemData(:,:)
!===================================================================================================================================

! Read in solution
CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)
CALL ReadAttribute(File_ID,'File_Type',1,StrScalar=File_Type)
OutputTime = 0. ! default
IF(TRIM(File_Type).NE.'RadiationState') CALL ReadAttribute(File_ID,'Time',1,RealScalar=OutputTime)
CALL GetDataSize(File_ID,TRIM(ArrayName),nDims,HSize)
nVarAdd=INT(HSize(1),4)
DEALLOCATE(HSize)

IF (nVarAdd.GT.0) THEN
  ALLOCATE(VarNamesAdd(1:nVarAdd))
  CALL ReadAttribute(File_ID,TRIM(VarName),nVarAdd,StrArray=VarNamesAdd(1:nVarAdd))
  ALLOCATE(ElemData(1:nVarAdd,1:nElems))

  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (&
        nVarAdd => INT(nVarAdd,IK) ,&
        offsetElem => INT(offsetElem,IK),&
        nElems     => INT(nElems,IK)    )
    CALL ReadArray(TRIM(ArrayName),2,(/nVarAdd, nElems/),offsetElem,2,RealArray=ElemData(1:nVarAdd,1:nElems))
  END ASSOCIATE
  ! Default
  IF(TRIM(File_Type).NE.'RadiationState') FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution_'//TRIM(ArrayName),OutputTime))//'.vtu'
  ! Special file types
  SELECT CASE(TRIM(File_Type))
    CASE('DSMCState','DSMCHOState')
      FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_visuDSMC',OutputTime))//'.vtu'
      IF(TRIM(ArrayName).EQ.'ExcitationData') FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_visuExcitationData',OutputTime))//'.vtu'
    CASE('RadiationState')
      WRITE(File_Num, "(I3.3)") iArgs
      FileString=TRIM(TRIM(ProjectName)//'_RadVisu_'//TRIM(File_Num))//'.vtu'
  END SELECT
  ! TODO: This is probably borked for NGeo>1 because then NodeCoords are not the corner nodes
  ! For the case of high-order output, utilize the CreateConnectivity routine analogous to the ConvertSurfaceData
  CALL WriteDataToVTK_PICLas(3,8,FileString,nVarAdd,VarNamesAdd(1:nVarAdd),nUniqueNodes,NodeCoords_Connect(1:3,1:nUniqueNodes),nElems,&
                              ElemData(1:nVarAdd,1:nElems),ElemUniqueNodeID(1:8,1:nElems))
END IF

SDEALLOCATE(VarNamesAdd)
SDEALLOCATE(ElemData)

CALL CloseDataFile()

END SUBROUTINE ConvertElemData


!===================================================================================================================================
!> Convert surface results to a VTK format
!===================================================================================================================================
SUBROUTINE ConvertSurfaceData(InputStateFile)
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: ProjectName
USE MOD_IO_HDF5                 ,ONLY: HSize
USE MOD_HDF5_Input              ,ONLY: OpenDataFile,CloseDataFile,ReadAttribute,GetDataSize,File_ID,ReadArray
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfSample
USE MOD_piclas2vtk_Vars         ,ONLY: SurfConnect, SurfOutputSideToUniqueSide
USE MOD_Interpolation           ,ONLY: GetVandermonde
USE MOD_Interpolation_Vars      ,ONLY: NMax,NMin
USE MOD_ChangeBasis             ,ONLY: ChangeBasis2D
USE MOD_Interpolation_Vars      ,ONLY: NodeTypeVISU
USE MOD_Mesh_Vars               ,ONLY: N_SurfMesh
USE MOD_ReadInTools             ,ONLY: PrintOption
#if !(PP_TimeDiscMethod==700)
USE MOD_DG_Vars                 ,ONLY: DG_Elems_master,DG_Elems_slave
#endif /*!(PP_TimeDiscMethod==700)*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)   :: InputStateFile
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)              :: FileString, File_Type
CHARACTER(LEN=255),ALLOCATABLE  :: VarNamesSurf_HDF5(:)
INTEGER                         :: nDims, nVarSurf, nSurfaceSidesReadin, SideID, iSurfOutputSide, Nloc
REAL                            :: OutputTime
REAL, ALLOCATABLE               :: tempSurfData(:,:,:,:,:)
REAL,ALLOCATABLE                :: NodeCoords_visu(:,:,:,:,:)     !< Coordinates of visualization nodes
!-----------------------------------------------------------------------------------------------------------------------------------
! Interpolation variables
!-----------------------------------------------------------------------------------------------------------------------------------
TYPE :: VisuInterpolation
REAL,ALLOCATABLE                :: Vdm_N_NVisu(:,:)                    !< Vandermonde from equidistant mesh to visualization nodes
END TYPE VisuInterpolation

TYPE(VisuInterpolation),ALLOCATABLE :: N_Inter_Visu(:)      !< Array of prebuild interpolation matrices
!===================================================================================================================================

! Read in solution
CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)
CALL ReadAttribute(File_ID,'File_Type',1,StrScalar=File_Type)
IF(TRIM(File_Type).NE.'RadiationSurfState') THEN
  CALL ReadAttribute(File_ID,'Time',1,RealScalar=OutputTime)
ELSE
  nSurfSample = 1
END IF
!CALL ReadAttribute(File_ID,'Time',1,RealScalar=OutputTime)
!CALL ReadAttribute(File_ID,'DSMC_nSurfSample',1,IntScalar=nSurfSample)
!IF(nSurfSample.NE.1) THEN
  !CALL abort(&
      !__STAMP__&
      !,'Error in piclas2vtk: Conversion to VTK only possible for DSMC_nSurfSample=1!')
!END IF
CALL ReadAttribute(File_ID,'DSMC_nSurfSample',1,IntScalar=nSurfSample)
CALL PrintOption('DSMC_nSurfSample','HDF5',IntOpt=nSurfSample) ! 'HDF5.'

CALL GetDataSize(File_ID,'SurfaceData',nDims,HSize)
nVarSurf = INT(HSize(1),4)
nSurfaceSidesReadin = INT(HSize(4),4)
ALLOCATE(VarNamesSurf_HDF5(nVarSurf))
CALL ReadAttribute(File_ID,'VarNamesSurface',nVarSurf,StrArray=VarNamesSurf_HDF5(1:nVarSurf))

IF((nVarSurf.GT.0).AND.(SurfConnect%nSurfaceOutputSides.GT.0))THEN
  ALLOCATE(tempSurfData(1:nVarSurf,nSurfSample,nSurfSample,0:0,1:nSurfaceSidesReadin))
  tempSurfData = 0.
  ASSOCIATE(nVarSurf            => INT(nVarSurf,IK),    &
            nSurfSample         => INT(nSurfSample,IK), &
            nSurfaceSidesReadin => INT(nSurfaceSidesReadin,IK))
    CALL ReadArray('SurfaceData',4,(/nVarSurf, nSurfSample, nSurfSample, nSurfaceSidesReadin/), &
                    0_IK,4,RealArray=tempSurfData(:,:,:,0,:))
  END ASSOCIATE

  ! Sanity check
  IF(SurfConnect%nSurfaceOutputSides.NE.nSurfaceSidesReadin)THEN
    WRITE (UNIT_stdOut,*) "SurfConnect%nSurfaceOutputSides =", SurfConnect%nSurfaceOutputSides
    WRITE (UNIT_stdOut,*) "nSurfaceSidesReadin         =", nSurfaceSidesReadin
    CALL abort(__STAMP__,'Error: SurfConnect%nSurfaceOutputSides.NE.nSurfaceSidesReadin')
  END IF ! SurfConnect%nSurfaceOutputSides.NE.nSurfaceSidesReadin
END IF

IF(nSurfSample.GT.1) THEN

  ! Build interpolation basis from NodeType=GAUSS (0:Nloc) to NodeType=VISU (0:nSurfSample)
  ALLOCATE(N_Inter_Visu(Nmin:Nmax))
  ! Loop over all polynomial degrees
  DO Nloc=Nmin,Nmax
    ! Use Nloc here because Face_xGP is built in mesh.f90 using Nloc when piclas2vtk loads the mesh
    ALLOCATE(N_Inter_Visu(Nloc)%Vdm_N_NVisu(0:nSurfSample,0:Nloc))
    ! Note that this uses NodeTypeVISU, which is hard-coded to NodeTypeVISU='VISU'
    CALL GetVandermonde(Nloc , 'GAUSS' , nSurfSample , NodeTypeVISU , N_Inter_Visu(Nloc)%Vdm_N_NVisu , modal=.FALSE.)
  END DO

  ALLOCATE(NodeCoords_visu(1:3,0:nSurfSample,0:nSurfSample,0:0,SurfConnect%nSurfaceOutputSides))
  NodeCoords_visu = 0.
  ! Interpolate mesh onto visu grid
  DO iSurfOutputSide = 1, SurfConnect%nSurfaceOutputSides
    ! Mapping from nSurfaceOutputSides to nSides
    SideID = SurfOutputSideToUniqueSide(iSurfOutputSide)
#if !(PP_TimeDiscMethod==700)
    Nloc   = MAX(DG_Elems_master(SideID),DG_Elems_slave(SideID))
#else
    Nloc = PP_N
#endif /*!(PP_TimeDiscMethod==700)*/
    CALL ChangeBasis2D(3, Nloc, nSurfSample, N_Inter_Visu(Nloc)%Vdm_N_NVisu, &
        N_SurfMesh(SideID)%Face_xGP(1:3 , 0:Nloc        , 0:Nloc)       , &
                    NodeCoords_visu(1:3 , 0:nSurfSample , 0:nSurfSample , 0 , iSurfOutputSide))
  END DO
ELSE
  ALLOCATE(NodeCoords_visu(1:3,0:0,0:0,0:0,1:SurfConnect%nSurfaceNode))
  NodeCoords_visu = 0.
  NodeCoords_visu(1:3,0,0,0,1:SurfConnect%nSurfaceNode) = SurfConnect%NodeCoords(1:3,1:SurfConnect%nSurfaceNode)
END IF

IF(TRIM(File_Type).NE.'RadiationSurfState') THEN
  IF(TRIM(File_Type).NE.'DSMCSurfChemState') THEN
    IF(TRIM(File_Type).NE.'DVMSurfState') THEN
    FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_visuSurf',OutputTime))//'.vtu'
  ELSE
      FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_visuSurfDVM',OutputTime))//'.vtu'
    END IF
  ELSE
    FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_visuSurfChem',OutputTime))//'.vtu'
  END IF
ELSE
  FileString=TRIM(TRIM(ProjectName)//'_RadSurfVisu')//'.vtu'
END IF
!CALL WriteDataToVTK_PICLas(4,FileString,nVarSurf,VarNamesSurf_HDF5,SurfConnect%nSurfaceNode,SurfConnect%NodeCoords(1:3,1:SurfConnect%nSurfaceNode),&
    !SurfConnect%nSurfaceOutputSides,SurfData,SurfConnect%SideSurfNodeMap(1:4,1:SurfConnect%nSurfaceOutputSides))

IF(nSurfSample.GT.1) THEN
  ! Number of nodes is calculated inside the routine in case of super-sampling: nSurfaceOutputSides = nSurfaceNode
  ! Connectivity is build inside, thus no SideSurfNodeMap is required
  CALL WriteDataToVTK_PICLas(2,4,FileString,nVarSurf,VarNamesSurf_HDF5,SurfConnect%nSurfaceOutputSides,NodeCoords_visu,&
      SurfConnect%nSurfaceOutputSides,tempSurfData,SurfConnect%SideSurfNodeMap)
ELSE
  CALL WriteDataToVTK_PICLas(2,4,FileString,nVarSurf,VarNamesSurf_HDF5,SurfConnect%nSurfaceNode,NodeCoords_visu,&
      SurfConnect%nSurfaceOutputSides,tempSurfData,SurfConnect%SideSurfNodeMap)
END IF

SDEALLOCATE(VarNamesSurf_HDF5)
SDEALLOCATE(tempSurfData)

CALL CloseDataFile()

END SUBROUTINE ConvertSurfaceData


!===================================================================================================================================
!> Create connectivities for matching sides
!===================================================================================================================================
SUBROUTINE BuildSurfMeshConnectivity(InputStateFile)
! MODULES
USE MOD_Globals
USE MOD_IO_HDF5            ,ONLY: HSize
USE MOD_HDF5_Input         ,ONLY: OpenDataFile,CloseDataFile,ReadAttribute,GetDataSize,File_ID,ReadArray,GetDataSize
USE MOD_Mesh_Tools         ,ONLY: GetCNElemID
USE MOD_Mesh_Vars          ,ONLY: BoundaryName, offsetElem, ElemToSide
USE MOD_piclas2vtk_Vars    ,ONLY: SurfConnect, SurfOutputSideToUniqueSide
USE MOD_Particle_Mesh_Vars ,ONLY: ElemSideNodeID_Shared,SideInfo_Shared,NodeCoords_Shared,NodeInfo_Shared,nNonUniqueGlobalSides
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)   :: InputStateFile
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255),ALLOCATABLE  :: SurfBCName_HDF5(:)
INTEGER                         :: nDims, iNode, nSurfBC_HDF5, iName
INTEGER                         :: ElemID, CNElemID, iLocSide, iSide, SideID, iNode2, iBC, SurfOutputSideID
INTEGER, ALLOCATABLE            :: TempBCSurfNodes(:), NonUniqueGlobalSideToSurfOutputSide(:)
REAL, ALLOCATABLE               :: TempNodeCoords(:,:)
LOGICAL                         :: IsSortedSurfNode
!===================================================================================================================================
CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)

SWRITE(UNIT_stdOut,'(A,A)')' GET NUMBER AND NAMES OF SURFACE-BCSIDES IN HDF5 FILE... '
CALL GetDataSize(File_ID,'BC_Surf',nDims,HSize,attrib=.true.)
nSurfBC_HDF5 = INT(HSize(1),4)

SWRITE(UNIT_stdOut,'(A3,A45,A3,I33,A13)')' | ','Number of Surface BCs',' | ',nSurfBC_HDF5,' | HDF5    | '
ALLOCATE(SurfBCName_HDF5(nSurfBC_HDF5))
CALL ReadAttribute(File_ID,'BC_Surf',nSurfBC_HDF5,StrArray=SurfBCName_HDF5)
DO iName = 1,nSurfBC_HDF5
  SWRITE(UNIT_stdOut,'(A3,A38,I2.1,A5,A3,A33,A13)')' | ','BC',iName,'Name',' | ',TRIM(SurfBCName_HDF5(iName)),' | HDF5    | '
END DO

! Create mapping from nNonUniqueGlobalSides (as read-in from mesh) to nSurfaceOutputSides (as read-in from SurfaceData)
ALLOCATE(NonUniqueGlobalSideToSurfOutputSide(1:nNonUniqueGlobalSides))
NonUniqueGlobalSideToSurfOutputSide(1:nNonUniqueGlobalSides) = -1

! Loop over all non-unique sides (not using nBCSides here to capture inner BCs as well)
SurfConnect%nSurfaceOutputSides=0
DO iSide = 1,nNonUniqueGlobalSides
  ! Cycle non-BC sides
  IF(SideInfo_Shared(SIDE_BCID,iSide).LE.0) CYCLE
  ! Loop over all read-in surface BCs
  DO iBC=1,nSurfBC_HDF5
    IF((TRIM(BoundaryName(SideInfo_Shared(SIDE_BCID,iSide))) .EQ. TRIM(SurfBCName_HDF5(iBC)))) THEN
      IF(SideInfo_Shared(SIDE_NBSIDEID,iSide).GT.0) THEN
        ! Cycling over non-unique surface sides with a neighbor side (in the output only showing the side with smaller index)
        IF(iSide.GT.SideInfo_Shared(SIDE_NBSIDEID,iSide)) CYCLE
      END IF
      ! Count the number of unique surface sides with a BC and create mapping
      SurfConnect%nSurfaceOutputSides = SurfConnect%nSurfaceOutputSides + 1
      NonUniqueGlobalSideToSurfOutputSide(iSide) = SurfConnect%nSurfaceOutputSides
    END IF
  END DO
END DO

! Build connectivity for the surface mesh
ALLOCATE(TempBCSurfNodes(4*SurfConnect%nSurfaceOutputSides))
TempBCSurfNodes = 0
ALLOCATE(TempNodeCoords(1:3,4*SurfConnect%nSurfaceOutputSides))
TempNodeCoords = 0.0
SDEALLOCATE(SurfConnect%SideSurfNodeMap)
ALLOCATE(SurfConnect%SideSurfNodeMap(1:4,1:SurfConnect%nSurfaceOutputSides))
SurfConnect%SideSurfNodeMap = 0
SurfConnect%nSurfaceNode=0
ALLOCATE(SurfOutputSideToUniqueSide(1:SurfConnect%nSurfaceOutputSides))
SurfOutputSideToUniqueSide(1:SurfConnect%nSurfaceOutputSides) = -1

DO iSide=1, nNonUniqueGlobalSides
  ! Cycling over non-reflective sides and non-unique surf sides (in the output only showing the side with smaller index)
  IF (NonUniqueGlobalSideToSurfOutputSide(iSide).EQ.-1) CYCLE
  ! Create mapping from unique surface output side to regular BC side
  SurfOutputSideID = NonUniqueGlobalSideToSurfOutputSide(iSide)
  ! Create connectivity between surface and nodes
  CNElemID = GetCNElemID(SideInfo_Shared(SIDE_ELEMID,iSide))
  iLocSide = SideInfo_Shared(SIDE_LOCALID,iSide)
  ! Saving SideID (1:nSides) for N_SurfMesh(SideID)%Face_xGP
  ElemID   = SideInfo_Shared(SIDE_ELEMID,iSide) - offsetElem
  SideID   = ElemToSide(E2S_SIDE_ID,iLocSide,ElemID)
  SurfOutputSideToUniqueSide(SurfOutputSideID) = SideID
  ! Loop over the nodes of the surface side
  DO iNode2 = 1, 4
    IsSortedSurfNode = .FALSE.
    DO iNode = 1, SurfConnect%nSurfaceNode
      IF (ABS(NodeInfo_Shared(ElemSideNodeID_Shared(iNode2, iLocSide, CNElemID)+1)).EQ.TempBCSurfNodes(iNode)) THEN
        SurfConnect%SideSurfNodeMap(iNode2,SurfOutputSideID) = iNode
        IsSortedSurfNode = .TRUE.
        EXIT
      END IF
    END DO
    IF(.NOT.IsSortedSurfNode) THEN
      SurfConnect%nSurfaceNode = SurfConnect%nSurfaceNode + 1
      TempBCSurfNodes(SurfConnect%nSurfaceNode) = ABS(NodeInfo_Shared(ElemSideNodeID_Shared(iNode2, iLocSide, CNElemID)+1))
      SurfConnect%SideSurfNodeMap(iNode2,SurfOutputSideID) = SurfConnect%nSurfaceNode
      TempNodeCoords(1:3,SurfConnect%nSurfaceNode) = NodeCoords_Shared(1:3,ElemSideNodeID_Shared(iNode2, iLocSide, CNElemID)+1)
    END IF
  END DO
END DO

SDEALLOCATE(SurfConnect%NodeCoords)
ALLOCATE(SurfConnect%NodeCoords(1:3,1:SurfConnect%nSurfaceNode))
SurfConnect%NodeCoords(1:3,1:SurfConnect%nSurfaceNode) = TempNodeCoords(1:3,1:SurfConnect%nSurfaceNode)
SDEALLOCATE(TempBCSurfNodes)
SDEALLOCATE(TempNodeCoords)
SDEALLOCATE(SurfBCName_HDF5)
CALL CloseDataFile()

END SUBROUTINE BuildSurfMeshConnectivity

!===================================================================================================================================
!> Subroutine to write 3D point data to VTK format using element-local polynomial degree
!===================================================================================================================================
SUBROUTINE WriteDataToVTKpAdaption(nVar,VarNames,FileString)
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars           ,ONLY: nElems
USE MOD_piclas2vtk_Vars     ,ONLY: ElemLocal, Nloc_Visu, PointToCellSwitch
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)          :: nVar                 !< Number of nodal output variables
CHARACTER(LEN=*),INTENT(IN) :: VarNames(nVar)       !< Names of all variables that will be written out
CHARACTER(LEN=*),INTENT(IN) :: FileString           !< Output file name
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iVar,ivtk
INTEGER                     :: nVTKPoints,nVTKCells
INTEGER                     :: nBytes,Offset
INTEGER                     :: INTdummy
INTEGER,PARAMETER           :: SizeINTdummy=STORAGE_SIZE(INTdummy)/8
REAL(KIND=4)                :: FLOATdummy
CHARACTER(LEN=35)           :: StrOffset,TempStr1,TempStr2
CHARACTER(LEN=200)          :: Buffer
CHARACTER(LEN=1)            :: lf
INTEGER                     :: ElemType,iElem
INTEGER,ALLOCATABLE         :: nodeids(:)
INTEGER                     :: Nloc,PointsPerVTKCell
REAL                        :: StartT,EndT ! Timer
!===================================================================================================================================
GETTIME(StartT)

PointsPerVTKCell = 8

SWRITE(UNIT_stdOut,'(A,I1,A)',ADVANCE='NO')"   WRITE 3D DATA TO VTX XML BINARY (VTU) FILE "
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') '['//TRIM(FileString)//'] ...'

nVTKPoints = 0
nVTKCells = 0

DO iElem = 1, nElems
  Nloc = Nloc_Visu(iElem) + PointToCellSwitch
  nVTKPoints = nVTKPoints + (Nloc+1)**3
  nVTKCells  = nVTKCells  +  Nloc**3
END DO

! Line feed character
lf = char(10)

! Write file
OPEN(NEWUNIT=ivtk,FILE=TRIM(FileString),ACCESS='STREAM')
! Write header
Buffer='<?xml version="1.0"?>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify file type
Buffer='  <UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
WRITE(TempStr1,'(I16)')nVTKPoints
WRITE(TempStr2,'(I16)')nVTKCells
Buffer='    <Piece NumberOfPoints="'//TRIM(ADJUSTL(TempStr1))//'" &
        &NumberOfCells="'//TRIM(ADJUSTL(TempStr2))//'">'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify point data
IF(PointToCellSwitch.EQ.0) THEN
  Buffer='      <PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=0
  WRITE(StrOffset,'(I16)')Offset
  DO iVar=1,nVar
    Buffer='        <DataArray type="Float32" Name="'//TRIM(VarNames(iVar))//'" '// &
                      'format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
    ! INTEGER KIND=4 check
    CHECKSAFEINT(INT(Offset,8)+INT(SizeINTdummy,8)+INT(nVTKPoints,8)*INT(SIZEOF_F(FLOATdummy),8),4)
    Offset=          Offset   +    SizeINTdummy   +    nVTKPoints   *    SIZEOF_F(FLOATdummy)
    WRITE(StrOffset,'(I16)')Offset
  END DO
  Buffer='      </PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify cell data
  Buffer='      <CellData> </CellData>'//lf;WRITE(ivtk) TRIM(Buffer)
ELSE IF(PointToCellSwitch.EQ.1) THEN
 ! Specify point data
  Buffer='      <PointData> </PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
  ! Specify cell data
  Buffer='      <CellData>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=0
  WRITE(StrOffset,'(I16)')Offset
  DO iVar=1,nVar
    Buffer='        <DataArray type="Float32" Name="'//TRIM(VarNames(iVar))//'" '// &
                      'format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
    ! INTEGER KIND=4 check
    CHECKSAFEINT(INT(Offset,8)+INT(SizeINTdummy,8)+INT(nVTKCells,8)*INT(SIZEOF_F(FLOATdummy),8),4)
    Offset=          Offset   +    SizeINTdummy   +    nVTKCells   *    SIZEOF_F(FLOATdummy)
    WRITE(StrOffset,'(I16)')Offset
  END DO
  Buffer='      </CellData>'//lf;WRITE(ivtk) TRIM(Buffer)
END IF
! Specify coordinate data
Buffer='      <Points>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='        <DataArray type="Float32" Name="Coordinates" NumberOfComponents="3" format="appended" '// &
                  'offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
! INTEGER KIND=4 check
CHECKSAFEINT(INT(Offset,8)+INT(SizeINTdummy,8)+3_8*INT(nVTKPoints,8)*INT(SIZEOF_F(FLOATdummy),8),4)
Offset=          Offset   +    SizeINTdummy   +3  *    nVTKPoints   *    SIZEOF_F(FLOATdummy)
WRITE(StrOffset,'(I16)')Offset
Buffer='      </Points>'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify necessary cell data
Buffer='      <Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
! Connectivity
Buffer='        <DataArray type="Int32" Name="connectivity" format="appended" '// &
                  'offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
! INTEGER KIND=4 check
CHECKSAFEINT(INT(Offset,8)+INT(SizeINTdummy,8)+INT(PointsPerVTKCell,8)*INT(nVTKCells,8)*INT(SIZEOF_F(FLOATdummy),8),4)
Offset=          Offset   +    SizeINTdummy   +    PointsPerVTKCell   *    nVTKCells   *    SIZEOF_F(FLOATdummy)
WRITE(StrOffset,'(I16)')Offset
! Offsets
Buffer='        <DataArray type="Int32" Name="offsets" format="appended" ' // &
                  'offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
! INTEGER KIND=4 check
CHECKSAFEINT(INT(Offset,8)+INT(SizeINTdummy,8)+INT(nVTKCells,8)*INT(SIZEOF_F(FLOATdummy),8),4)
Offset=          Offset   +    SizeINTdummy   +    nVTKCells   *    SIZEOF_F(FLOATdummy)
WRITE(StrOffset,'(I16)')Offset
! Elem types
Buffer='        <DataArray type="Int32" Name="types" format="appended" '// &
                  'offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='      </Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='    </Piece>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='  </UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
! Prepare append section
Buffer='  <AppendedData encoding="raw">'//lf;WRITE(ivtk) TRIM(Buffer)
! Write leading data underscore
Buffer='_';WRITE(ivtk) TRIM(Buffer)

! Write binary raw data into append section
! Solution data: write size per variable, and follow-up with complete data set per variable
IF(PointToCellSwitch.EQ.0) THEN
  ! Write data as PointData
  nBytes = nVTKPoints*SIZEOF_F(FLOATdummy)
ELSEIF(PointToCellSwitch.EQ.1) THEN
  ! Write data as CellData
  nBytes = nVTKCells*SIZEOF_F(FLOATdummy)
END IF
DO iVar=1,nVar
  WRITE(ivtk) nBytes
  DO iElem = 1, nElems
    Nloc = Nloc_Visu(iElem)
    WRITE(ivtk) REAL(ElemLocal(iElem)%U_Visu(iVar,0:Nloc,0:Nloc,0:Nloc),4)
  END DO
END DO       ! iVar
! Coordinates
nBytes = nVTKPoints*SIZEOF_F(FLOATdummy)*3
WRITE(ivtk) nBytes
DO iElem = 1, nElems
  Nloc = Nloc_Visu(iElem) + PointToCellSwitch
  WRITE(ivtk) REAL(ElemLocal(iElem)%Coords_NVisu(1:3,0:Nloc,0:Nloc,0:Nloc),4)
END DO

! Connectivity and footer
SDEALLOCATE(nodeids)
ALLOCATE(nodeids(PointsPerVTKCell*nVTKCells))
CALL CreateConnectivitypAdaption(PointsPerVTKCell*nVTKCells,nodeids)

nBytes = PointsPerVTKCell*nVTKCells*SizeINTdummy
WRITE(ivtk) nBytes
WRITE(ivtk) nodeids
! Offset
nBytes = nVTKCells*SizeINTdummy
WRITE(ivtk) nBytes
WRITE(ivtk) (Offset,Offset=PointsPerVTKCell,PointsPerVTKCell*nVTKCells,PointsPerVTKCell)
! Elem type
ElemType = 12 ! VTK_HEXAHEDRON
WRITE(ivtk) nBytes
WRITE(ivtk) (ElemType,iElem=1,nVTKCells)

DEALLOCATE(nodeids)

! Footer
lf = char(10)
Buffer=lf//'  </AppendedData>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
CLOSE(ivtk)
GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, ' DONE!', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
END SUBROUTINE WriteDataToVTKpAdaption


!===================================================================================================================================
!> Create connectivity for output of element local polynomial degrees (p-Adaption)
!===================================================================================================================================
SUBROUTINE CreateConnectivitypAdaption(nNodes,nodeids)
USE ISO_C_BINDING
USE MOD_Globals
USE MOD_Mesh_Vars           ,ONLY: nElems
USE MOD_piclas2vtk_Vars     ,ONLY: Nloc_Visu, PointToCellSwitch
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: nNodes             !< number of nodes
INTEGER,INTENT(OUT) :: nodeids(1:nNodes)  !< stores the connectivity
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i,j,k,iElem
INTEGER           :: NodeID,NodeIDElem
INTEGER           :: Nloc, Nloc_elem, Nloc_p1_2
!===================================================================================================================================

! create connectivity
NodeID = 0
NodeIDElem = 0
DO iElem=1,nElems
  Nloc = Nloc_Visu(iElem) + PointToCellSwitch
  Nloc_elem = (Nloc+1)**3
  Nloc_p1_2 = (Nloc+1)**2
  DO k=1,Nloc
    DO j=1,Nloc
      DO i=1,Nloc
        NodeID=NodeID+1
        nodeids(NodeID) = NodeIDElem+i+   j   *(Nloc+1)+(k-1)*Nloc_p1_2-1 !P4(CGNS=tecVisu standard)
        NodeID=NodeID+1
        nodeids(NodeID) = NodeIDElem+i+  (j-1)*(Nloc+1)+(k-1)*Nloc_p1_2-1 !P1
        NodeID=NodeID+1
        nodeids(NodeID) = NodeIDElem+i+1+(j-1)*(Nloc+1)+(k-1)*Nloc_p1_2-1 !P2
        NodeID=NodeID+1
        nodeids(NodeID) = NodeIDElem+i+1+ j   *(Nloc+1)+(k-1)*Nloc_p1_2-1 !P3
        NodeID=NodeID+1
        nodeids(NodeID) = NodeIDElem+i+   j   *(Nloc+1)+ k   *Nloc_p1_2-1 !P8
        NodeID=NodeID+1
        nodeids(NodeID) = NodeIDElem+i+  (j-1)*(Nloc+1)+ k   *Nloc_p1_2-1 !P5
        NodeID=NodeID+1
        nodeids(NodeID) = NodeIDElem+i+1+(j-1)*(Nloc+1)+ k   *Nloc_p1_2-1 !P6
        NodeID=NodeID+1
        nodeids(NodeID) = NodeIDElem+i+1+ j   *(Nloc+1)+ k   *Nloc_p1_2-1 !P7
      END DO
    END DO
  END DO
  NodeIDElem=NodeIDElem+Nloc_elem
END DO
END SUBROUTINE CreateConnectivitypAdaption