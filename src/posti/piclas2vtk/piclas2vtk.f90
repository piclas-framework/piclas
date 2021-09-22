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
MODULE MOD_piclas2vtk_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE

REAL,ALLOCATABLE            :: NodeCoords_Connect(:,:)
INTEGER,ALLOCATABLE         :: ElemUniqueNodeID(:,:)
INTEGER                     :: nUniqueNodes
!----------------------------------------------------------------------------------------------------------------------------------
! Mapping of nodes and surface sides, required for connectivity of elements
!----------------------------------------------------------------------------------------------------------------------------------
TYPE tSurfaceConnect
  INTEGER                         :: nSurfaceNode                 !< Number of Nodes on Surface (reflective)
  INTEGER                         :: nSurfaceBCSides              !< Number of Sides on Surface (reflective)
  REAL, ALLOCATABLE               :: NodeCoords(:,:)
  INTEGER, ALLOCATABLE            :: SideSurfNodeMap(:,:)         !< Mapping from glob Side to SurfaceNodeNum (1:4, nSurfaceBCSides)
END TYPE

TYPE (tSurfaceConnect)               :: SurfConnect

END MODULE MOD_piclas2vtk_Vars

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
USE MOD_Mesh_Vars             ,ONLY: nGlobalElems
USE MOD_Interpolation         ,ONLY: DefineParametersInterpolation,InitInterpolation
USE MOD_Mesh                  ,ONLY: DefineParametersMesh,InitMesh
#ifdef PARTICLES
USE MOD_Particle_Mesh_Tools   ,ONLY: InitParticleGeometry, InitElemNodeIDs
USE MOD_Particle_Mesh_Vars    ,ONLY: NodeCoords_Shared, ElemNodeID_Shared, NodeInfo_Shared
#endif /*PARTICLES*/
USE MOD_Mesh_Tools            ,ONLY: InitGetCNElemID, InitGetGlobalElemID
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars       ,ONLY: nComputeNodeTotalElems
#endif /*USE_MPI*/
USE MOD_Preproc
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
LOGICAL                        :: DGSolutionExists, ElemDataExists, SurfaceDataExists, VisuParticles, PartDataExists
LOGICAL                        :: BGFieldExists
LOGICAL                        :: VisuAdaptiveInfo, AdaptiveInfoExists
LOGICAL                        :: ReadMeshFinished, ElemMeshInit, SurfMeshInit
#ifdef PARTICLES
INTEGER                        :: iElem, iNode
#endif /*PARTICLES*/
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
CALL ParseCommandlineArguments()
!CALL DefineParametersMPI()
!CALL DefineParametersIO_HDF5()
! Define parameters for piclas2vtk
CALL prms%SetSection("piclas2vtk")
CALL prms%CreateStringOption( 'NodeTypeVisu',"Node type of the visualization basis: "//&
                                             "VISU,GAUSS,GAUSS-LOBATTO,CHEBYSHEV-GAUSS-LOBATTO", 'VISU')
CALL prms%CreateIntOption(    'NVisu',       "Number of points at which solution is sampled for visualization.")
CALL prms%CreateLogicalOption('VisuParticles',  "Visualize particles (velocity, species, internal energy).", '.FALSE.')
CALL prms%CreateLogicalOption('VisuAdaptiveInfo', "Visualize the sampled values utilized for the adaptive surface flux and porous BC (velocity, density, pumping speed)", '.FALSE.')
CALL prms%CreateIntOption(    'TimeStampLength', 'Length of the floating number time stamp', '21')
CALL prms%CreateLogicalOption( 'meshCheckWeirdElements'&
  , 'Abort when weird elements are found: it means that part of the element is turned inside-out. ','.TRUE.')
CALL DefineParametersIO()
CALL DefineParametersMesh()
CALL DefineParametersInterpolation()

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

! Set necessary parameters for piclas2vtk tool
! If no parameter file has been set, the standard values will be used
NodeTypeVisuOut  = GETSTR('NodeTypeVisu','VISU')    ! Node type of visualization basis
VisuParticles    = GETLOGICAL('VisuParticles','.FALSE.')
VisuAdaptiveInfo    = GETLOGICAL('VisuAdaptiveInfo','.FALSE.')
! Initialization of I/O routines
CALL InitIOHDF5()
! Get length of the floating number time stamp
TimeStampLength = GETINT('TimeStampLength')
IF((TimeStampLength.LT.4).OR.(TimeStampLength.GT.30)) CALL abort(&
    __STAMP__&
    ,'TimeStampLength cannot be smaller than 4 and not larger than 30')
WRITE(UNIT=TimeStampLenStr ,FMT='(I0)') TimeStampLength
WRITE(UNIT=TimeStampLenStr2,FMT='(I0)') TimeStampLength-4

! Measure init duration
Time=PICLASTIME()
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F14.2,A)') ' INITIALIZATION DONE! [',Time-StartTime,' sec ]'
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
  IF(.NOT.ISVALIDHDF5FILE(InputStateFile)) THEN
    CALL Abort(__STAMP__,'ERROR - Please supply only .h5 files after parameter file.')
  END IF

  SWRITE(UNIT_stdOut,'(132("="))')
  SWRITE(UNIT_stdOut,'(A,I3,A,I3,A)') 'Processing state ',iArgs-iArgsStart+1,' of ',nArgs-iArgsStart+1,'...'

  ! Open .h5 file
  DGSolutionExists   = .FALSE.
  ElemDataExists     = .FALSE.
  SurfaceDataExists  = .FALSE.
  PartDataExists     = .FALSE.
  AdaptiveInfoExists = .FALSE.
  CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
  ! Get the type of the .h5 file
  CALL ReadAttribute(File_ID,'File_Type',1,StrScalar=File_Type)
  ! Check which containers are present
  CALL DatasetExists(File_ID , 'DG_Solution'  , DGSolutionExists)
  IF(TRIM(File_Type).EQ.'PartStateBoundary') DGSolutionExists = .FALSE.
  CALL DatasetExists(File_ID , 'ElemData'     , ElemDataExists)
  CALL DatasetExists(File_ID , 'AdaptiveInfo' , AdaptiveInfoExists)
  CALL DatasetExists(File_ID , 'SurfaceData'  , SurfaceDataExists)
  CALL DatasetExists(File_ID , 'PartData'     , PartDataExists)
  CALL DatasetExists(File_ID , 'BGField'      , BGFieldExists) ! deprecated , but allow for backward compatibility
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
    CALL InitMesh(1,MeshFile_IN=MeshFile)
    CALL InitGetGlobalElemID()
    CALL InitGetCNElemID()
#if USE_MPI
    nComputeNodeTotalElems = nGlobalElems
#endif /*USE_MPI*/
#ifdef PARTICLES
    ! ElemSideNodeID_Shared is required for the connectivity
    CALL InitParticleGeometry()
    ! ElemNodeID_Shared is required for the connectivity
    CALL InitElemNodeIDs()
#endif /*PARTICLES*/
    ReadMeshFinished = .TRUE.
  END IF
  ! Build connectivity for element/volume output
  IF(ElemDataExists.AND..NOT.ElemMeshInit) THEN
#ifdef PARTICLES
    CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
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
    ElemMeshInit = .TRUE.
#else
    SWRITE(*,*) 'WARNING - Conversion of ElemData requires the compilation of piclas2vtk with PARTICLES set to ON!'
    SWRITE(*,*) 'WARNING - File contains ElemData but output is skipped.'
    ElemDataExists = .FALSE.
#endif
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
    CALL ConvertDGSolution(InputStateFile,NVisu,NodeTypeVisuOut,File_Type,DGSolutionDataset)
  END IF
  ! === ElemData ===================================================================================================================
  IF(ElemDataExists) THEN
    CALL ConvertElemData(InputStateFile,'ElemData','VarNamesAdd')
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
      CALL ConvertElemData(InputStateFile,'AdaptiveInfo','VarNamesAdaptive')
    END IF
  END IF
END DO ! iArgs = 2, nArgs

! Measure processing duration
Time=PICLASTIME()
#if USE_MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) THEN
  CALL abort(__STAMP__,&
    'MPI finalize error',iError)
END IF
#endif
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F14.2,A)') ' piclas2vtk FINISHED! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')

END PROGRAM piclas2vtk


!===================================================================================================================================
!> Subroutine to write 3D point data to VTK format
!===================================================================================================================================
SUBROUTINE WriteDataToVTK_PICLas(data_size,FileString,nVar,VarNameVisu,nNodes,Coords,nElems,Value,ConnectInfo)
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: nVar,nElems,nNodes,data_size           !> Number of nodal output variables
REAL,INTENT(IN)               :: Coords(1:3,nNodes), Value(nVar,nElems)
CHARACTER(LEN=*),INTENT(IN)   :: FileString, VarNameVisu(nVar)          !> Output file name
INTEGER,INTENT(IN)            :: ConnectInfo(data_size,nElems)          !> Statevector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iVal,iElem,Offset,nBytes,nVTKPoints,nVTKCells,ivtk=44,iVar,iNode, int_size, ElemType
INTEGER                       :: Vertex(data_size,nElems), iLen, str_len
CHARACTER(LEN=35)             :: StrOffset,TempStr1,TempStr2
CHARACTER(LEN=200)            :: Buffer, tmp, tmp2, VarNameString
CHARACTER(LEN=1)              :: lf, components_string
REAL(KIND=4)                  :: float
INTEGER,ALLOCATABLE           :: VarNameCombine(:), VarNameCombineLen(:)
!===================================================================================================================================

nVTKPoints=nNodes
nVTKCells=nElems

SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')"   WRITE 3D DATA TO VTX XML BINARY (VTU) FILE "
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') '['//TRIM(FileString)//'] ...'
IF(nElems.LT.1)THEN
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
  RETURN
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
  ! Compare the strings, while omitting the last character to find identify VeloX/Y/Z as a vector
  IF (TRIM(tmp(:iLen-1)) .EQ. TRIM(tmp2(:iLen-1))) THEN
    ! Although the translational temperature is given in X/Y/Z its not a vector (VisIt/Paraview would produce a magnitude variable)
    IF(INDEX(tmp(:iLen-1),'TempTrans').EQ.0) THEN
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
    WRITE(ivtk) nBytes,REAL(Value(iVal,1:nVTKCells),4)
  ELSEIF(VarNameCombine(iVal).EQ.1) THEN
    WRITE(ivtk) nBytes*VarNameCombineLen(iVal),REAL(Value(iVal:iVal+VarNameCombineLen(iVal)-1,1:nVTKCells),4)
  ENDIF
END DO
! Points
nBytes = 3*nVTKPoints*INT(SIZEOF(FLOAT),4)
WRITE(ivtk) nBytes
WRITE(ivtk) REAL(Coords(1:3,1:nVTKPoints),4)
! Connectivity
DO iElem=1,nVTKCells
  DO iNode=1,data_size
    Vertex(iNode,iElem) = ConnectInfo(iNode,iElem)-1
  END DO
END DO
nBytes = data_size*nVTKCells*INT(SIZEOF(int_size),4)
WRITE(ivtk) nBytes
WRITE(ivtk) Vertex(:,:)
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
  CALL abort(__STAMP__,&
      'Wrong data size given to WritaDataToVTK_PICLas routine!')
END SELECT
WRITE(ivtk) nBytes
WRITE(ivtk) (ElemType,iElem=1,nVTKCells)
! Write footer
Buffer=lf//'  </AppendedData>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
CLOSE(ivtk)
SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"

SDEALLOCATE(VarNameCombine)
SDEALLOCATE(VarNameCombineLen)

END SUBROUTINE WriteDataToVTK_PICLas


!===================================================================================================================================
!> Convert the output of the field solver to a VTK output format
!===================================================================================================================================
SUBROUTINE ConvertDGSolution(InputStateFile,NVisu,NodeTypeVisuOut,OutputName,DGSolutionDataset)
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
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN) :: InputStateFile,NodeTypeVisuOut,DGSolutionDataset
CHARACTER(LEN=*),INTENT(IN)   :: OutputName
INTEGER,INTENT(IN)            :: NVisu
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iElem, iDG, nVar_State, N_State, nElems_State, nVar_Solution, nDims, iField, nFields, Suffix
INTEGER                         :: nDimsOffset
CHARACTER(LEN=255)              :: MeshFile, NodeType_State, FileString_DG, StrVarNamesTemp(4)
CHARACTER(LEN=255),ALLOCATABLE  :: StrVarNames(:), StrVarNamesTemp2(:)
REAL                            :: OutputTime
REAL,ALLOCATABLE                :: U2(:,:,:,:,:,:)                   !< Solution from state file with additional dimension, rank=6
REAL,ALLOCATABLE                :: U(:,:,:,:,:)                      !< Solution from state file, rank=5
REAL,ALLOCATABLE,TARGET         :: U_Visu(:,:,:,:,:)                 !< Solution on visualization nodes
REAL,POINTER                    :: U_Visu_p(:,:,:,:,:)               !< Solution on visualization nodes
REAL,ALLOCATABLE                :: Coords_NVisu(:,:,:,:,:)           !< Coordinates of visualization nodes
REAL,ALLOCATABLE,TARGET         :: Coords_DG(:,:,:,:,:)
REAL,POINTER                    :: Coords_DG_p(:,:,:,:,:)
REAL,ALLOCATABLE                :: Vdm_EQNgeo_NVisu(:,:)             !< Vandermonde from equidistant mesh to visualization nodes
REAL,ALLOCATABLE                :: Vdm_N_NVisu(:,:)                  !< Vandermonde from state to visualization nodes
LOGICAL                         :: DGSourceExists
CHARACTER(LEN=16)               :: hilf
!===================================================================================================================================
! 1.) Open given file to get the number of elements, the order and the name of the mesh file
CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)

! Read-in of dimensions of the field array (might have an additional dimension, i.e., rank is 6 instead of 5)
CALL GetDataSize(File_ID,TRIM(DGSolutionDataset),nDims,HSize)
! Check the number of fields in the file, if more than 5 dimensions, the 6th dimensions carries the number of fields
nFields     = MERGE(1 , INT(HSize(nDims)) , nDims.EQ.5)
nDimsOffset = MERGE(0 , 1                 , nDims.EQ.5)

! Default: DGSolutionDataset = 'DG_Solution'
CALL GetDataProps(TRIM(DGSolutionDataset),nVar_Solution,N_State,nElems_State,NodeType_State,nDimsOffset)
CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)

! Check if the DG_Source container exists, and if it does save the variable names and increase the nVar_State variable
CALL DatasetExists(File_ID,'DG_Source',DGSourceExists)
IF (DGSourceExists) THEN
  CALL ReadAttribute(File_ID,'VarNamesSource',4,StrArray=StrVarNamesTemp)
  nVar_State = nVar_Solution + 4
ELSE
  nVar_State = nVar_Solution
END IF
! Save the variable names for the regular DG_Solution in a temporary array
SDEALLOCATE(StrVarNamesTemp2)
ALLOCATE(StrVarNamesTemp2(nVar_Solution))
CALL ReadAttribute(File_ID,'VarNames',nVar_Solution,StrArray=StrVarNamesTemp2)

! Allocate the variable names array used for the output and copy the names from the DG_Solution and DG_Source (if it exists)
SDEALLOCATE(StrVarNames)
ALLOCATE(StrVarNames(nVar_State))
StrVarNames(1:nVar_Solution) = StrVarNamesTemp2
IF(DGSourceExists) StrVarNames(nVar_Solution+1:nVar_State) = StrVarNamesTemp(1:4)

! Check if a state file is converted. Read the time stamp from .h5
IF(OutputName.EQ.'State') CALL ReadAttribute(File_ID,'Time',1,RealScalar=OutputTime)
CALL CloseDataFile()

SDEALLOCATE(Vdm_EQNgeo_NVisu)
ALLOCATE(Vdm_EQNgeo_NVisu(0:Ngeo,0:NVisu))
CALL GetVandermonde(Ngeo,NodeTypeVisu,NVisu,NodeTypeVisuOut,Vdm_EQNgeo_NVisu,modal=.FALSE.)
SDEALLOCATE(Coords_NVisu)
ALLOCATE(Coords_NVisu(3,0:NVisu,0:NVisu,0:NVisu,nElems))
SDEALLOCATE(Coords_DG)
ALLOCATE(Coords_DG(3,0:NVisu,0:NVisu,0:NVisu,nElems))

! Convert coordinates to visu grid
DO iElem = 1,nElems
  CALL ChangeBasis3D(3,NGeo,NVisu,Vdm_EQNgeo_NVisu,NodeCoords(:,:,:,:,iElem),Coords_NVisu(:,:,:,:,iElem))
END DO

! Read in solution
CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nVar_State    => INT(nVar_State    , IK)      , &
      nVar_Solution => INT(nVar_Solution , IK)      , &
      offsetElem    => INT(offsetElem    , IK)      , &
      N_State       => INT(N_State       , IK)      , &
      nElems        => INT(nElems        , IK)      , &
      nFields       => INT(nFields       , IK)      )
  IF(nFields.EQ.1)THEN
    SDEALLOCATE(U)
    ALLOCATE(U(nVar_State,0:N_State,0:N_State,0:N_State,nElems))
    ! Default: DGSolutionDataset = 'DG_Solution'
    CALL ReadArray(TRIM(DGSolutionDataset),5,(/nVar_Solution,N_State+1_IK,N_State+1_IK,N_State+1_IK,nElems/),offsetElem,5, &
                    RealArray=U(1:nVar_Solution,0:N_State,0:N_State,0:N_State,1:nElems))
    IF(DGSourceExists) CALL ReadArray('DG_Source',5,(/4_IK,N_State+1_IK,N_State+1_IK,N_State+1_IK,nElems/),offsetElem,5, &
                                       RealArray=U(nVar_Solution+1:nVar_State,0:N_State,0:N_State,0:N_State,1:nElems))
  ELSE ! more than one field
    SDEALLOCATE(U2)
    ALLOCATE(U2(nVar_State,0:N_State,0:N_State,0:N_State,nElems,nFields))
    ! Default: DGSolutionDataset = 'DG_Solution'
    CALL ReadArray(TRIM(DGSolutionDataset),6,(/nVar_Solution,N_State+1_IK,N_State+1_IK,N_State+1_IK,nElems,nFields/),offsetElem,5, &
                    RealArray=U2(1:nVar_Solution,0:N_State,0:N_State,0:N_State,1:nElems,1:nFields))
  END IF ! nFields.GT.1
END ASSOCIATE

CALL CloseDataFile()

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

! Output multiples files, if the DGSolution file contains more than one field (e.g. multiple different times for the same field)
DO iField = 1, nFields

  ! Write solution to vtk
  IF(OutputName.EQ.'State')THEN
    FileString_DG=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution',OutputTime))//'.vtu'
  ELSE
    FileString_DG=TRIM(ProjectName)//'_'//TRIM(OutputName)//'.vtu'
  END IF

  IF(nFields.GT.1)THEN
    ! Append 001, 002, etc. to file name (before the .vtk suffix)
    Suffix = INDEX(FileString_DG,'.vtu')
    WRITE(UNIT=hilf,FMT='(I3.3)') iField
    hilf = TRIM(hilf)//'.vtu'
    FileString_DG(Suffix:Suffix-1+LEN(TRIM(hilf)))=TRIM(hilf)
  END IF ! nFields.GT.1

  ! Interpolate solution to visu grid
  iDG = 0
  DO iElem = 1,nElems
    iDG = iDG + 1
    IF(nFields.EQ.1)THEN
      CALL ChangeBasis3D(nVar_State,N_State,NVisu,Vdm_N_NVisu,U(:,:,:,:,iElem),U_Visu(:,:,:,:,iDG))
    ELSE ! more than one field
      CALL ChangeBasis3D(nVar_State,N_State,NVisu,Vdm_N_NVisu,U2(:,:,:,:,iElem,iField),U_Visu(:,:,:,:,iDG))
    END IF ! nFields.GT.1
  END DO

  ! Output to VTK
  CALL WriteDataToVTK(nVar_State,NVisu,iDG,StrVarNames,Coords_DG_p,U_Visu_p,TRIM(FileString_DG),dim=3,DGFV=0)

END DO ! iField = 1, nFields

SDEALLOCATE(StrVarNamesTemp2)
SDEALLOCATE(StrVarNames)
SDEALLOCATE(Vdm_EQNgeo_NVisu)
SDEALLOCATE(Coords_NVisu)
SDEALLOCATE(Coords_DG)
SDEALLOCATE(U)
SDEALLOCATE(U2)
SDEALLOCATE(Vdm_N_NVisu)
SDEALLOCATE(U_Visu)

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
REAL                            :: OutputTime, FileVersionHDF5
LOGICAL                         :: FileVersionExists
!===================================================================================================================================

CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)
CALL ReadAttribute(File_ID,'Time',1,RealScalar=OutputTime)

! check file version
CALL DatasetExists(File_ID,'File_Version',FileVersionExists,attrib=.TRUE.)
IF (FileVersionExists) THEN
  CALL ReadAttribute(File_ID,'File_Version',1,RealScalar=FileVersionHDF5)
ELSE
  CALL abort(&
      __STAMP__&
      ,'Error in InitRestart(): Attribute "File_Version" does not exist!')
END IF
IF(FileVersionHDF5.LT.1.5)THEN
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
  RealInfoOpt=FileVersionHDF5)
END IF ! FileVersionHDF5.LT.1.5

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
CALL WriteDataToVTK_PICLas(1,FileString,nPartsVar,VarNamesParticle,nParts,PartData(1:3,1:nParts),nParts,&
                            PartData(4:nPartsVar+3,1:nParts),ConnectInfo(1,1:nParts))

SDEALLOCATE(VarNamesParticle)
SDEALLOCATE(tmpArray)
SDEALLOCATE(PartData)
SDEALLOCATE(ConnectInfo)

CALL CloseDataFile()

END SUBROUTINE ConvertPartData


SUBROUTINE ConvertElemData(InputStateFile,ArrayName,VarName)
!===================================================================================================================================
!> Convert element/volume data (single value per cell, e.g. DSMC/BGK results) to a VTK format
!===================================================================================================================================
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
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)              :: FileString, File_Type
REAL                            :: OutputTime
INTEGER                         :: nDims,nVarAdd
CHARACTER(LEN=255),ALLOCATABLE  :: VarNamesAdd(:)
REAL,ALLOCATABLE                :: ElemData(:,:)
!===================================================================================================================================

! Read in solution
CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)
CALL ReadAttribute(File_ID,'Time',1,RealScalar=OutputTime)
CALL ReadAttribute(File_ID,'File_Type',1,StrScalar=File_Type)
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
  SELECT CASE(TRIM(File_Type))
    CASE('State')
      FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution_'//TRIM(ArrayName),OutputTime))//'.vtu'
    CASE('DSMCState','DSMCHOState')
      FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_visuDSMC',OutputTime))//'.vtu'
  END SELECT
  ! TODO: This is probably borked for NGeo>1 because then NodeCoords are not the corner nodes
  CALL WriteDataToVTK_PICLas(8,FileString,nVarAdd,VarNamesAdd(1:nVarAdd),nUniqueNodes,NodeCoords_Connect(1:3,1:nUniqueNodes),nElems,&
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
USE MOD_Globals
USE MOD_Globals_Vars    ,ONLY: ProjectName
USE MOD_IO_HDF5         ,ONLY: HSize
USE MOD_HDF5_Input      ,ONLY: OpenDataFile,CloseDataFile,ReadAttribute,GetDataSize,File_ID,ReadArray
USE MOD_piclas2vtk_Vars ,ONLY: SurfConnect
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)   :: InputStateFile
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)              :: FileString
CHARACTER(LEN=255),ALLOCATABLE  :: VarNamesSurf_HDF5(:)
INTEGER                         :: nDims, nVarSurf, nSurfSample, nSurfaceSidesReadin
REAL                            :: OutputTime
REAL, ALLOCATABLE               :: tempSurfData(:,:,:,:), SurfData(:,:), Coords(:,:)
INTEGER                         :: iSide
!===================================================================================================================================

! Read in solution
CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)
CALL ReadAttribute(File_ID,'Time',1,RealScalar=OutputTime)
CALL ReadAttribute(File_ID,'DSMC_nSurfSample',1,IntScalar=nSurfSample)
IF(nSurfSample.NE.1) THEN
  CALL abort(&
      __STAMP__&
      ,'Error in piclas2vtk: Conversion to VTK only possible for DSMC_nSurfSample=1!')
END IF

CALL GetDataSize(File_ID,'SurfaceData',nDims,HSize)
nVarSurf = INT(HSize(1),4)
nSurfaceSidesReadin = INT(HSize(4),4)
ALLOCATE(VarNamesSurf_HDF5(nVarSurf))
CALL ReadAttribute(File_ID,'VarNamesSurface',nVarSurf,StrArray=VarNamesSurf_HDF5(1:nVarSurf))

IF((nVarSurf.GT.0).AND.(SurfConnect%nSurfaceBCSides.GT.0))THEN
  ALLOCATE(SurfData(1:nVarSurf,1:SurfConnect%nSurfaceBCSides))
  ALLOCATE(tempSurfData(1:nVarSurf,nSurfSample,nSurfSample,1:nSurfaceSidesReadin))
  SurfData=0.
  tempSurfData = 0.
  ASSOCIATE(nVarSurf        => INT(nVarSurf,IK),  &
            nSurfaceSidesReadin => INT(nSurfaceSidesReadin,IK))
    CALL ReadArray('SurfaceData',4,(/nVarSurf, 1_IK, 1_IK, nSurfaceSidesReadin/), &
                    0_IK,4,RealArray=tempSurfData(:,:,:,:))
  END ASSOCIATE

  ! Sanity check
  IF(SurfConnect%nSurfaceBCSides.NE.nSurfaceSidesReadin)THEN
    WRITE (UNIT_stdOut,*) "SurfConnect%nSurfaceBCSides =", SurfConnect%nSurfaceBCSides
    WRITE (UNIT_stdOut,*) "nSurfaceSidesReadin         =", nSurfaceSidesReadin
    CALL abort(&
    __STAMP__&
    ,'Error: SurfConnect%nSurfaceBCSides.NE.nSurfaceSidesReadin')
  END IF ! SurfConnect%nSurfaceBCSides.NE.nSurfaceSidesReadin

  ! Copy data from tmp array
  DO iSide = 1, SurfConnect%nSurfaceBCSides
    SurfData(1:nVarSurf,iSide) = tempSurfData(1:nVarSurf,1,1,iSide)
  END DO ! iSide = 1, SurfConnect%nSurfaceBCSides
END IF

FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_visuSurf',OutputTime))//'.vtu'
CALL WriteDataToVTK_PICLas(4,FileString,nVarSurf,VarNamesSurf_HDF5,SurfConnect%nSurfaceNode,SurfConnect%NodeCoords(1:3,1:SurfConnect%nSurfaceNode),&
    SurfConnect%nSurfaceBCSides,SurfData,SurfConnect%SideSurfNodeMap(1:4,1:SurfConnect%nSurfaceBCSides))

SDEALLOCATE(VarNamesSurf_HDF5)
SDEALLOCATE(SurfData)
SDEALLOCATE(tempSurfData)
SDEALLOCATE(Coords)

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
USE MOD_Mesh_ReadIn        ,ONLY: readMesh
USE MOD_Mesh_Tools         ,ONLY: GetCNElemID
USE MOD_Mesh_Vars          ,ONLY: BoundaryName
USE MOD_piclas2vtk_Vars    ,ONLY: SurfConnect
USE MOD_Particle_Mesh_Vars ,ONLY: ElemSideNodeID_Shared,SideInfo_Shared,NodeCoords_Shared,NodeInfo_Shared
USE MOD_Particle_Mesh_Vars ,ONLY: nNonUniqueGlobalSides
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
INTEGER                         :: nDims, iNode, nSurfBC_HDF5, iName, nSides
INTEGER                         :: CNElemID, iLocSide, iSide, iNode2, iBC, nSurfSides, nUniqueSurfSide
INTEGER, ALLOCATABLE            :: TempBCSurfNodes(:), TempSideSurfNodeMap(:,:), SideToSurfSide(:), NonUnique2UniqueSide(:)
REAL, ALLOCATABLE               :: TempNodeCoords(:,:)
LOGICAL                         :: IsSortedSurfNode
!===================================================================================================================================
CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)

SWRITE(UNIT_stdOut,'(A,A)')' GET NUMBER AND NAMES OF SURFACE-BCSIDES IN HDF5 FILE... '
CALL GetDataSize(File_ID,'BC_Surf',nDims,HSize,attrib=.true.)
nSurfBC_HDF5 = INT(HSize(1),4)

SWRITE(UNIT_stdOut,'(A3,A45,A3,I33,A13)')' | ','Number of Surface BCs',' | ',nSurfBC_HDF5,' | HDF5    | '
ALLOCATE(SurfBCName_HDF5(nSurfBC_HDF5))
CALL ReadAttribute(File_ID,'BC_Surf',nSurfBC_HDF5,StrArray=SurfBCName_HDF5)
DO iName = 1,nSurfBC_HDF5
  SWRITE(UNIT_stdOut,'(A3,A38,I2.1,A5,A3,A33,A13)')' | ','BC',iName,'Name',' | ',TRIM(SurfBCName_HDF5(iName)),' | HDF5    | '
END DO

nSides = nNonUniqueGlobalSides

! create sideid to surfaceid map for all surface-sides contained in statefile
ALLOCATE(SideToSurfSide(1:nSides))
SideToSurfSide(1:nSides) = -1
ALLOCATE(NonUnique2UniqueSide(1:nSides))
NonUnique2UniqueSide(1:nSides) = -1
! if side is surface (BC_reflective defined in State file) then map respective surface side number
nSurfSides = 0
nUniqueSurfSide = 0
DO iSide = 1,nSides
  IF(SideInfo_Shared(SIDE_BCID,iSide).EQ.0) CYCLE
  DO iBC=1,nSurfBC_HDF5
    IF((TRIM(BoundaryName(SideInfo_Shared(SIDE_BCID,iSide))) .EQ. TRIM(SurfBCName_HDF5(iBC)))) THEN
      nSurfSides = nSurfSides + 1
      SideToSurfSide(iSide) = nSurfSides
      IF(SideInfo_Shared(SIDE_NBSIDEID,iSide).GT.0) THEN
        ! Cycling over non-unique sides
        IF(iSide.GT.SideInfo_Shared(SIDE_NBSIDEID,iSide)) CYCLE
      END IF
      nUniqueSurfSide = nUniqueSurfSide + 1
      NonUnique2UniqueSide(nSurfSides) = nUniqueSurfSide
    END IF
  END DO
END DO

! Build connectivity for the surface mesh
ALLOCATE(TempBCSurfNodes(4*nSides))
ALLOCATE(TempNodeCoords(1:3,4*nSides))
ALLOCATE(TempSideSurfNodeMap(1:4,1:nSides))
SurfConnect%nSurfaceNode=0
SurfConnect%nSurfaceBCSides=0

DO iSide=1, nSides
  ! Cycling over non-reflective sides
  IF (SideToSurfSide(iSide).EQ.-1) CYCLE
  ! Cycling over non-unique sides
  IF (NonUnique2UniqueSide(SideToSurfSide(iSide)).EQ.-1) CYCLE
  SurfConnect%nSurfaceBCSides = SurfConnect%nSurfaceBCSides + 1
  CNElemID = GetCNElemID(SideInfo_Shared(SIDE_ELEMID,iSide))
  iLocSide = SideInfo_Shared(SIDE_LOCALID,iSide)
  DO iNode2 = 1, 4
    IsSortedSurfNode = .FALSE.
    DO iNode = 1, SurfConnect%nSurfaceNode
      IF (ABS(NodeInfo_Shared(ElemSideNodeID_Shared(iNode2, iLocSide, CNElemID)+1)).EQ.TempBCSurfNodes(iNode)) THEN
        TempSideSurfNodeMap(iNode2,SurfConnect%nSurfaceBCSides) = iNode
        IsSortedSurfNode = .TRUE.
        EXIT
      END IF
    END DO
    IF(.NOT.IsSortedSurfNode) THEN
      SurfConnect%nSurfaceNode = SurfConnect%nSurfaceNode + 1
      TempBCSurfNodes(SurfConnect%nSurfaceNode) = ABS(NodeInfo_Shared(ElemSideNodeID_Shared(iNode2, iLocSide, CNElemID)+1))
      TempSideSurfNodeMap(iNode2,SurfConnect%nSurfaceBCSides) = SurfConnect%nSurfaceNode
      TempNodeCoords(1:3,SurfConnect%nSurfaceNode) = NodeCoords_Shared(1:3,ElemSideNodeID_Shared(iNode2, iLocSide, CNElemID)+1)
    END IF
  END DO
END DO

SDEALLOCATE(SurfConnect%SideSurfNodeMap)
ALLOCATE(SurfConnect%SideSurfNodeMap(1:4,1:SurfConnect%nSurfaceBCSides))
SurfConnect%SideSurfNodeMap(1:4,1:SurfConnect%nSurfaceBCSides) = TempSideSurfNodeMap(1:4,1:SurfConnect%nSurfaceBCSides)
SDEALLOCATE(SurfConnect%NodeCoords)
ALLOCATE(SurfConnect%NodeCoords(1:3,1:SurfConnect%nSurfaceNode))
SurfConnect%NodeCoords(1:3,1:SurfConnect%nSurfaceNode) = TempNodeCoords(1:3,1:SurfConnect%nSurfaceNode)
SDEALLOCATE(TempBCSurfNodes)
SDEALLOCATE(TempSideSurfNodeMap)
SDEALLOCATE(TempNodeCoords)
SDEALLOCATE(NonUnique2UniqueSide)
SDEALLOCATE(SurfBCName_HDF5)
SDEALLOCATE(SideToSurfSide)
CALL CloseDataFile()

END SUBROUTINE BuildSurfMeshConnectivity


