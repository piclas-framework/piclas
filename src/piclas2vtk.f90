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
!> - useCurveds: Logical, should the mesh be curved or not (if the mesh itself is curved)
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
  INTEGER                         :: nSurfaceNode                 ! Number of Nodes on Surface (reflective)
  INTEGER                         :: nSurfaceBCSides              ! Number of Sides on Surface (reflective)
  INTEGER, ALLOCATABLE            :: NonUnique2UniqueSide(:)      ! Mapping required in the case of inner boundaries
  INTEGER, ALLOCATABLE            :: BCSurfNodes(:)               ! Nodes on Surface (reflective) (nSurfaceNode)
  REAL, ALLOCATABLE               :: NodeCoords(:,:)
  INTEGER, ALLOCATABLE            :: SideSurfNodeMap(:,:)         ! Mapping from glob Side to SurfaceNodeNum (1:4, nSurfaceBCSides)
END TYPE

TYPE (tSurfaceConnect)               :: SurfConnect

END MODULE MOD_piclas2vtk_Vars

PROGRAM piclas2vtk
! MODULES
USE MOD_Globals_Vars
USE MOD_StringTools
USE MOD_Commandline_Arguments
USE MOD_Globals               ,ONLY: doPrintHelp,CollectiveStop,MPIRoot,iError,MPI_COMM_WORLD,PiclasTime,UNIT_stdOut,StartTime,abort
USE MOD_Globals               ,ONLY: SetStackSizeUnlimited
USE MOD_IO_HDF5               ,ONLY: InitIOHDF5,DefineParametersIO
USE MOD_MPI                   ,ONLY: InitMPI
USE MOD_ReadInTools           ,ONLY: prms,PrintDefaultParameterFile
USE MOD_ReadInTools           ,ONLY: GETINT,GETSTR,GETLOGICAL
USE MOD_HDF5_Input            ,ONLY: OpenDataFile,CloseDataFile,GetDataProps,ReadAttribute,File_ID,ReadArray,GetDataSize,DatasetExists
USE MOD_HDF5_Input            ,ONLY: ISVALIDHDF5FILE,ISVALIDMESHFILE
USE MOD_Mesh_ReadIn           ,ONLY: readMesh
USE MOD_Mesh                  ,ONLY: FinalizeMesh
!USE MOD_Mesh_Vars             ,ONLY: NGeo,nElems,offsetElem
USE MOD_Mesh_Vars             ,ONLY: useCurveds,NodeCoords
!USE MOD_Interpolation_Vars    ,ONLY: NodeTypeVisu
!USE MOD_Interpolation         ,ONLY: GetVandermonde
!USE MOD_ChangeBasis           ,ONLY: ChangeBasis3D
!USE MOD_VTK                   ,ONLY: WriteDataToVTK,WriteVTKMultiBlockDataSet
#if USE_MPI
!USE MOD_MPI_Vars              ,ONLY: NbProc,nMPISides_Proc
#endif /*USE_MPI*/
!USE MOD_Analyze               ,ONLY: CalcErrorStateFiles, CalcErrorStateFileSigma
!USE MOD_Interpolation_Vars    ,ONLY: NAnalyze
!USE MOD_Mesh_Vars             ,ONLY: sJ,NGeoRef
USE MOD_Preproc
!USE MOD_Metrics               ,ONLY: CalcMetricsErrorDiff
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Time                              ! Used to track computation time
CHARACTER(LEN=255)             :: NodeTypeVisuOut                   ! Stores user selected type of visualization nodes
INTEGER                        :: NVisu                             ! Polynomial degree of visualization
INTEGER                        :: iArgs                             ! Loop counters
!INTEGER                        :: iElem                             ! Loop counters
INTEGER                        :: iExt                              ! Stores position where the filename extension begins
CHARACTER(LEN=255)             :: InputStateFile
!CHARACTER(LEN=255)             :: MeshFile
!INTEGER                        :: nVar_State,N_State,nElems_State   ! Properties read from state file
!INTEGER                        :: nVar_Solution
!CHARACTER(LEN=255)             :: NodeType_State                    !     "
REAL,ALLOCATABLE               :: U(:,:,:,:,:)                      ! Solution from state file
!REAL,ALLOCATABLE               :: U_first(:,:,:,:,:)                ! Solution from state file
!REAL,ALLOCATABLE               :: U_average(:,:,:,:,:)              ! Solution from state file
!INTEGER                        :: N_average
REAL,ALLOCATABLE,TARGET        :: U_Visu(:,:,:,:,:)                 ! Solution on visualiation nodes
!REAL,POINTER                   :: U_Visu_p(:,:,:,:,:)               ! Solution on visualiation nodes
REAL,ALLOCATABLE               :: Coords_NVisu(:,:,:,:,:)           ! Coordinates of visualisation nodes
!REAL,ALLOCATABLE,TARGET        :: Coords_DG(:,:,:,:,:)
!REAL,POINTER                   :: Coords_DG_p(:,:,:,:,:)
REAL,ALLOCATABLE               :: Vdm_EQNgeo_NVisu(:,:)             ! Vandermonde from equidistand mesh to visualisation nodes
REAL,ALLOCATABLE               :: Vdm_N_NVisu(:,:)                  ! Vandermonde from state to visualisation nodes
INTEGER                        :: nGeo_old,nVar_State_old           ! Variables used to check if we need to reinitialize
INTEGER                        :: N_State_old,nElems_old            !     "
!INTEGER                        :: N_State_first                     ! first state file
CHARACTER(LEN=255)             :: MeshFile_old                      !     "
CHARACTER(LEN=255)             :: NodeType_State_old                !     "
!CHARACTER(LEN=255)             :: FileString_DG
!CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:), StrVarNamesTemp2(:)
!CHARACTER(LEN=255)             :: StrVarNamesTemp(4)
!REAL                           :: OutputTime
!INTEGER                        :: iDG
!CHARACTER(LEN=255)             :: FileString_multiblock
LOGICAL                        :: CmdLineMode, NVisuDefault         ! In command line mode only NVisu is specified directly,
                                                                    ! otherwise a parameter file is needed
CHARACTER(LEN=2)               :: NVisuString                       ! String containing NVisu from command line option
CHARACTER(LEN=20)              :: fmtString                         ! String containing options for formatted write
LOGICAL                        :: CalcDiffError                     ! Use first state file as reference state for L2 error
                                                                    ! calculation with the following state files
LOGICAL                        :: AllowChangedMesh
LOGICAL                        :: CalcDiffSigma                     ! Use last state file as state for L2 sigma calculation
LOGICAL                        :: CalcAverage                       ! Calculate and write arithmetic mean of alle StateFile
!LOGICAL                        :: DGSourceExists, skip
LOGICAL                        :: DGSolutionExists, ElemDataExists, SurfaceDataExists
!CHARACTER(LEN=40)              :: DefStr
INTEGER                        :: iArgsStart
LOGICAL                        :: ReadMeshFinished, ElemMeshInit, SurfMeshInit
LOGICAL                        :: VisuParticles, PartDataExists, BGFieldExists
INTEGER                        :: TimeStampLength
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
CALL prms%CreateLogicalOption('useCurveds',  "Controls usage of high-order information in mesh. Turn off to discard "//&
                                             "high-order data and treat curved meshes as linear meshes.", '.TRUE.')
CALL prms%CreateLogicalOption('CalcDiffError',  "Use first state file as reference state for L2 error calculation "//&
                                                "with the following state files.", '.TRUE.')
CALL prms%CreateLogicalOption('AllowChangedMesh',"Neglect changes mesh, use inits of first mesh (ElemID must match!).", '.FALSE.')
CALL prms%CreateLogicalOption('CalcDiffSigma',  "Use last state file as state for L2 sigma calculation.", '.FALSE.')
CALL prms%CreateLogicalOption('CalcAverage',    "Calculate and write arithmetic mean of alle StateFile.", '.FALSE.')
CALL prms%CreateLogicalOption('VisuSource',     "use DG_Source instead of DG_Solution.", '.FALSE.')
CALL prms%CreateIntOption(    'NAnalyze'         , 'Polynomial degree at which analysis is performed (e.g. for L2 errors).\n'//&
                                                   'Default: 2*N. (needed for CalcDiffError)')
CALL prms%CreateLogicalOption('VisuParticles',  "Visualize particles (velocity, species, internal energy).", '.FALSE.')
CALL prms%CreateLogicalOption('writePartitionInfo',  "Write information about MPI partitions into a file.",'.FALSE.')
CALL prms%CreateIntOption(    'TimeStampLength', 'Length of the floating number time stamp', '21')
CALL DefineParametersIO()

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
    prms%maxNameLen  = 13 ! gatheredWrite is the longest option
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
useCurveds       = GETLOGICAL('useCurveds','.FALSE.')  ! Allow curved mesh or not
CalcDiffError    = GETLOGICAL('CalcDiffError','.FALSE.')
IF(CalcDiffError)THEN
  IF(nArgs.LT.3)CALL abort(__STAMP__,&
      'CalcDiffError needs a minimum of two state files!',iError)
  CalcDiffSigma  = GETLOGICAL('CalcDiffSigma','.FALSE.')
  AllowChangedMesh  = GETLOGICAL('AllowChangedMesh','.FALSE.')
  CalcAverage    = .FALSE.
ELSE
  CalcDiffSigma  = .FALSE.
  AllowChangedMesh = .FALSE. !dummy(?)
  CalcAverage    = GETLOGICAL('CalcAverage','.FALSE.')
END IF
VisuParticles    = GETLOGICAL('VisuParticles','.FALSE.')
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

! Initialize an "old" state to check against - used to determine if we need to reinitialize some variables
nVar_State_old     = 0
N_State_old        = 0
MeshFile_old       = ''
nGeo_old           = 0
nElems_old         = 0
NodeType_State_old = ''

IF(NVisuDefault) THEN
  iArgsStart = 1
ELSE
  iArgsStart = 2
END IF

ReadMeshFinished = .FALSE.
ElemMeshInit = .FALSE.
SurfMeshInit = .FALSE.

! Loop over remaining supplied .h5 files
DO iArgs = iArgsStart,nArgs
  InputStateFile = Args(iArgs)
  ! Check if the argument is a valid .h5 file
  ! IF(.NOT.ISVALIDHDF5FILE(InputStateFile)) THEN
  !   CALL CollectiveStop(__STAMP__,&
  !     'ERROR - Please supply only .h5 files after parameter file.')
  ! END IF

  SWRITE(UNIT_stdOut,'(132("="))')
  SWRITE(UNIT_stdOut,'(A,I3,A,I3,A)') 'Processing state ',iArgs-iArgsStart+1,' of ',nArgs-iArgsStart+1,'...'

  ! Open .h5 file
  DGSolutionExists = .FALSE.; ElemDataExists = .FALSE.; SurfaceDataExists = .FALSE.; BGFieldExists = .FALSE.
  CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
  CALL DatasetExists(File_ID,'DG_Solution',DGSolutionExists)
  CALL DatasetExists(File_ID,'ElemData',ElemDataExists)
  CALL DatasetExists(File_ID,'SurfaceData',SurfaceDataExists)
  CALL DatasetExists(File_ID,'PartData',PartDataExists)
  CALL DatasetExists(File_ID,'BGField',BGFieldExists)

  IF(ElemDataExists.OR.SurfaceDataExists) THEN
    CALL ReadMesh_piclas2vtk(InputStateFile,SurfaceDataExists,ElemMeshInit,SurfMeshInit)
  END IF
  ! === DG_Solution ================================================================================================================
  ! Read in parameters from the State file
!   IF(DGSolutionExists) THEN
!     CALL GetDataProps('DG_Solution',nVar_Solution,N_State,nElems_State,NodeType_State)
!     CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile)
!     CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)

!     ! Check if the DG_Source container exists, and if it does save the variable names and increase the nVar_State variable
!     CALL DatasetExists(File_ID,'DG_Source',DGSourceExists)
!     IF (DGSourceExists) THEN
!       CALL ReadAttribute(File_ID,'VarNamesSource',4,StrArray=StrVarNamesTemp)
!       nVar_State = nVar_Solution + 4
!     ELSE
!       nVar_State = nVar_Solution
!     END IF
!     ! Save the variable names for the regular DG_Solution in a temporary array
!     SDEALLOCATE(StrVarNamesTemp2)
!     ALLOCATE(StrVarNamesTemp2(nVar_Solution))
!     CALL ReadAttribute(File_ID,'VarNames',nVar_Solution,StrArray=StrVarNamesTemp2)

!     ! Allocate the variable names array used for the output and copy the names from the DG_Solution and DG_Source (if it exists)
!     IF (nVar_State.NE.nVar_State_old) THEN
!       SDEALLOCATE(StrVarNames)
!       ALLOCATE(StrVarNames(nVar_State))
!     END IF
!     StrVarNames(1:nVar_Solution) = StrVarNamesTemp2
!     IF (DGSourceExists) THEN
!       StrVarNames(nVar_Solution+1:nVar_State) = StrVarNamesTemp(1:4)
!     END IF

!     CALL ReadAttribute(File_ID,'Time',1,RealScalar=OutputTime)
!     CALL CloseDataFile()
!     ! Check if the mesh has changed
!     IF(CalcDiffError.AND.(iArgs.GT.2).AND.AllowChangedMesh)THEN
!       skip=.TRUE.
!     ELSE
!       skip=.FALSE.
!     END IF
!     IF (.NOT.skip .AND. (TRIM(MeshFile).NE.TRIM(MeshFile_old))) THEN
!       IF(CalcDiffError.AND.(iArgs.GT.2))CALL abort(__STAMP__,&
!       'CalcDiffError needs identical meshes!',iError)
!       IF(CalcAverage.AND.(iArgs.GT.2))CALL abort(__STAMP__,&
!       'CalcAverage needs identical meshes!',iError)
!       ! Check if the file is a valid mesh
!       IF(.NOT.ISVALIDMESHFILE(MeshFile)) THEN
!         CALL CollectiveStop(__STAMP__,&
!           'ERROR - Not a valid mesh file.')
!       END IF
!       ! Deallocate and finalize mesh vars
!       SDEALLOCATE(NodeCoords)
! #if USE_MPI
!       SDEALLOCATE(NbProc)
!       SDEALLOCATE(nMPISides_Proc)
! #endif /*USE_MPI*/
!       CALL FinalizeMesh()

!       ! Read in parameters from mesh file
!       CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
!       CALL ReadAttribute(File_ID,'Ngeo',1,IntegerScalar=NGeo)
!       CALL CloseDataFile()

!       ! Read the mesh itself
!       CALL readMesh(MeshFile)
!       ReadMeshFinished = .TRUE.

!       IF(CalcAverage .AND. iArgs.GT.2)THEN
!         IF (nVar_State_old     .NE. nVar_State .OR.&
!             N_State_old        .NE. N_State .OR.&
!             nGeo_old           .NE. nGeo .OR.&
!             nElems_old         .NE. nElems .OR.&
!             NodeType_State_old .NE. NodeType_State) CALL abort(__STAMP__,&
!         'CalcAverage: different N,nGeo,etc. not implemented yet!')
!       END IF

!       ! Check if ne need to realloacte the Vandermonde from mesh to visualization
!       IF (NGeo.NE.nGeo_old) THEN
!         SDEALLOCATE(Vdm_EQNgeo_NVisu)
!         ALLOCATE(Vdm_EQNgeo_NVisu(0:Ngeo,0:NVisu))
!         CALL GetVandermonde(Ngeo,NodeTypeVisu,NVisu,NodeTypeVisuOut,Vdm_EQNgeo_NVisu,modal=.FALSE.)
!       END IF

!       ! Check if we need to reallocate the coordinate array
!       IF (nElems.NE.nElems_old) THEN
!         SDEALLOCATE(Coords_NVisu)
!         ALLOCATE(Coords_NVisu(3,0:NVisu,0:NVisu,0:NVisu,nElems))
!         SDEALLOCATE(Coords_DG)
!         ALLOCATE(Coords_DG(3,0:NVisu,0:NVisu,0:NVisu,nElems))
!       END IF

!       ! Convert coordinates to visu grid
!       DO iElem = 1,nElems
!         CALL ChangeBasis3D(3,NGeo,NVisu,Vdm_EQNgeo_NVisu,NodeCoords(:,:,:,:,iElem),Coords_NVisu(:,:,:,:,iElem))
!       END DO
!     END IF ! New mesh

!     ! Check if we need to reallocate the solution array
!     IF ((N_State.NE.N_State_old).OR.(nVar_State.NE.nVar_State_old).OR.(nElems.NE.nElems_old)) THEN
!       SDEALLOCATE(U)
!       ALLOCATE(U(nVar_State,0:N_State,0:N_State,0:N_State,nElems))
!     END IF

!     ! Read in solution
!     CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)

!     ! Associate construct for integer KIND=8 possibility
!     ASSOCIATE (&
!           nVar_State => INT(nVar_State,IK) ,&
!           nVar_Solution => INT(nVar_Solution,IK) ,&
!           offsetElem => INT(offsetElem,IK),&
!           N_State    => INT(N_State,IK),&
!           nElems     => INT(nElems,IK)    )
!       CALL ReadArray('DG_Solution',5,(/nVar_Solution,N_State+1_IK,N_State+1_IK,N_State+1_IK,nElems/),offsetElem,5, &
!                       RealArray=U(1:nVar_Solution,0:N_State,0:N_State,0:N_State,1:nElems))
!       IF(DGSourceExists) THEN
!         CALL ReadArray('DG_Source',5,(/4_IK,N_State+1_IK,N_State+1_IK,N_State+1_IK,nElems/),offsetElem,5, &
!                         RealArray=U(nVar_Solution+1:nVar_State,0:N_State,0:N_State,0:N_State,1:nElems))
!       END IF
!     END ASSOCIATE

!     IF(CalcDiffError)THEN
!       IF(iArgs .EQ. 2)THEN
!         ! Set the default analyze polynomial degree NAnalyze to 2*(N+1)
!         WRITE(DefStr,'(i4)') 2*(N_State+1)
!         NAnalyze=GETINT('NAnalyze',DefStr)
!         ! Copy state to 'first'
!         ALLOCATE(U_first(nVar_State,0:N_State,0:N_State,0:N_State,nElems))
!         U_first       = U
!         N_State_first = N_State

!         PP_N=N_State
!         NGeoRef=3*NGeo ! build jacobian at higher degree
!         ALLOCATE(sJ            (  0:N_State,0:N_State,0:N_State,nElems))
!         CALL CalcMetricsErrorDiff()
!       ELSE IF(iArgs .EQ. nArgs .AND. CalcDiffSigma) THEN
!         IF(NAnalyze.LT.2*(N_State+1))CALL abort(__STAMP__,&
!       'CalcDiffError: NAnalyze.LT.2*(N_State+1)! The polynomial degree is too small!',iError)
!         CALL CalcErrorStateFileSigma(nVar_State,N_State,U)

!       ELSE
!         IF(NAnalyze.LT.2*(N_State+1))CALL abort(__STAMP__,&
!       'CalcDiffError: NAnalyze.LT.2*(N_State+1)! The polynomial degree is too small!',iError)
!         CALL CalcErrorStateFiles(nVar_State,N_State_first,N_State,U_first,U)
!       END IF
!     END IF
!     IF(CalcAverage)THEN
!       IF(iArgs .EQ. 2)THEN
!         ALLOCATE(U_average(nVar_State,0:N_State,0:N_State,0:N_State,nElems))
!         U_average = U
!         N_average = 1
!       ELSE
!         U_average = U_average+U
!         N_average = N_average+1
!       END IF
!       IF (iArgs .EQ. nArgs) THEN
!         U = U_average/REAL(N_average)
!       END IF
!     END IF

!     CALL CloseDataFile()

!     IF(CalcAverage .AND. (iArgs.GT.2 .AND. iArgs.LT.nArgs)) THEN
!       CYCLE !go to next file (only output the averaged U, but allocate and set _old-stuff for iArg=2)
!     END IF

!     ! Check if we need to reallocate the Vandermonde from state to visualisation
!     IF ((N_State.NE.N_State_old).OR.(TRIM(NodeType_State).NE.TRIM(NodeType_State_old))) THEN
!       SDEALLOCATE(Vdm_N_NVisu)
!       ALLOCATE(Vdm_N_NVisu(0:N_State,0:NVisu))
!       CALL GetVandermonde(N_State,NodeType_State,NVisu,NodeTypeVisuOut,Vdm_N_NVisu,modal=.FALSE.)
!     END IF

!     ! Check if we need to reallocate the visualisation array
!     IF ((nVar_State.NE.nVar_State_old).OR.(nElems.NE.nElems_old)) THEN
!       SDEALLOCATE(U_Visu)
!       ALLOCATE(U_Visu(nVar_State,0:NVisu,0:NVisu,0:NVisu,nElems))
!     END IF

!     IF(iArgs.EQ.nArgs .OR. .NOT.CalcAverage) THEN
!       ! Interpolate solution to visu grid
!       iDG = 0
!       DO iElem = 1,nElems
!         iDG = iDG + 1
!         CALL ChangeBasis3D(nVar_State,N_State,NVisu,Vdm_N_NVisu,U(:,:,:,:,iElem),U_Visu(:,:,:,:,iDG))
!         Coords_DG(:,:,:,:,iDG) = Coords_NVisu(:,:,:,:,iElem)
!       END DO

!       ! Write solution to vtk
!       FileString_DG=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution',OutputTime))//'.vtu'
!       Coords_DG_p => Coords_DG(:,:,:,:,1:iDG)
!       U_Visu_p => U_Visu(:,:,:,:,1:iDG)
!       CALL WriteDataToVTK(nVar_State,NVisu,iDG,StrVarNames,Coords_DG_p,U_Visu_p,TRIM(FileString_DG),dim=3,DGFV=0)
!       IF (MPIRoot) THEN
!         ! write multiblock file
!         FileString_multiblock=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution',OutputTime))//'.vtm'
!         CALL WriteVTKMultiBlockDataSet(FileString_multiblock,FileString_DG)
!       END IF
!     END IF !iArgs.EQ.nArgs .OR. .NOT.CalcAverage

!     ! Save parameters of this state to later check if we need to reinitialize variables
!     nVar_State_old     = nVar_State
!     N_State_old        = N_State
!     MeshFile_old       = MeshFile
!     nGeo_old           = nGeo
!     nElems_old         = nElems
!     NodeType_State_old = NodeType_State
!   END IF
  ! === ElemData ===================================================================================================================
  IF(ElemDataExists) THEN
    CALL ConvertElemData(InputStateFile)
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
  ! === BField =====================================================================================================================
  IF(BGFieldExists) THEN
    CALL ConvertBGField(InputStateFile,ReadMeshFinished,NVisu,NodeTypeVisuOut)
  END IF
END DO ! iArgs = 2, nArgs

! Finalize
SDEALLOCATE(Vdm_N_NVisu)
SDEALLOCATE(Vdm_EQNgeo_NVisu)
SDEALLOCATE(U)
SDEALLOCATE(U_Visu)
SDEALLOCATE(Coords_NVisu)
SDEALLOCATE(NodeCoords)

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

SUBROUTINE WriteDataToVTK_PICLas(data_size,FileString,nVar,VarNameVisu,nNodes,Coords,nElems,Value,ConnectInfo)
!===================================================================================================================================
! Subroutine to write 3D point data to VTK format
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: nVar,nElems,nNodes,data_size                                 ! Number of nodal output variables
REAL,INTENT(IN)               :: Coords(1:3,nNodes), Value(nVar,nElems)
CHARACTER(LEN=*),INTENT(IN)   :: FileString, VarNameVisu(nVar)   ! Output file name
INTEGER,INTENT(IN)            :: ConnectInfo(data_size,nElems)      ! Statevector
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


SUBROUTINE ReadMesh_piclas2vtk(InputStateFile,SurfaceDataExists,ElemMeshInit,SurfMeshInit)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_IO_HDF5             ,ONLY: HSize
USE MOD_HDF5_Input          ,ONLY: OpenDataFile,CloseDataFile,ReadAttribute,GetDataSize,File_ID,ReadArray,GetDataSize
USE MOD_Mesh_ReadIn         ,ONLY: readMesh
USE MOD_Mesh_Vars           ,ONLY: nElems, nGlobalElems, NGeo
USE MOD_Particle_Mesh_Vars  ,ONLY: nNonUniqueGlobalSides,nNonUniqueGlobalNodes
USE MOD_Particle_Mesh_Vars  ,ONLY: ElemNodeID_Shared,ElemSideNodeID_Shared,NodeCoords_Shared,SideInfo_Shared,NodeInfo_Shared
USE MOD_Particle_Mesh_Vars  ,ONLY: ElemInfo_Shared
USE MOD_piclas2vtk_Vars     ,ONLY: NodeCoords_Connect, ElemUniqueNodeID, nUniqueNodes
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN) :: InputStateFile
LOGICAL,INTENT(IN)            :: SurfaceDataExists
LOGICAL,INTENT(INOUT)         :: ElemMeshInit, SurfMeshInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)            :: MeshFile
INTEGER                       :: nDims, nSideIDs, CornerNodeIDswitch(8), NodeMap(4,6), iElem ,iNode
INTEGER                       :: FirstElemInd, LastElemInd, GlobalSideID, nlocSides, localSideID, nStart, iLocSide
INTEGER                       :: iSide, jlocSide, NbElemID, NbSideID, nlocSidesNb, sideCount
!===================================================================================================================================
IF(.NOT.ElemMeshInit) THEN
  CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
  CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile)
  SWRITE(UNIT_stdOut,'(A)')'READ MESH FROM DATA FILE "'//TRIM(MeshFile)//'" ...'
  SWRITE(UNIT_StdOut,'(132("-"))')
  CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
  CALL ReadAttribute(File_ID,'Ngeo',1,IntegerScalar=NGeo)
  CALL ReadAttribute(File_ID,'nSides',1,IntegerScalar=nNonUniqueGlobalSides)
  CALL ReadAttribute(File_ID,'nNodes',1,IntegerScalar=nNonUniqueGlobalNodes)
  CALL ReadAttribute(File_ID,'nUniqueNodes',1,IntegerScalar=nUniqueNodes)
  CALL GetDataSize(File_ID,'ElemInfo',nDims,HSize)
  CHECKSAFEINT(HSize(2),4)
  nGlobalElems=INT(HSize(2),4)

  FirstElemInd = 1
  LastElemInd = nGlobalElems
  nElems = nGlobalElems

  ASSOCIATE(ELEMINFOSIZE_H5_IK => INT(ELEMINFOSIZE_H5,IK) ,&
            SIDEINFOSIZE_H5_IK => INT(SIDEINFOSIZE_H5,IK) ,&
            nElems             => INT(nElems,IK))
    ALLOCATE(ElemInfo_Shared(ELEMINFOSIZE_H5,FirstElemInd:LastElemInd))
    CALL ReadArray('ElemInfo',2,(/ELEMINFOSIZE_H5_IK,nElems/),0_IK,2,IntegerArray_i4=ElemInfo_Shared(1:ELEMINFOSIZE_H5,:))

    nSideIDs     = ElemInfo_Shared(ELEM_LASTSIDEIND,LastElemInd)-ElemInfo_Shared(ELEM_FIRSTSIDEIND,FirstElemInd)
    ALLOCATE(SideInfo_Shared(SIDEINFOSIZE,1:nSideIDs))
    ASSOCIATE(nSideIDs           => INT(nSideIDs,IK))
      SideInfo_Shared = 0
      CALL ReadArray('SideInfo',2,(/SIDEINFOSIZE_H5_IK,nSideIDs/),0_IK,2,IntegerArray_i4=SideInfo_Shared(1:SIDEINFOSIZE_H5,:))
    END ASSOCIATE
  END ASSOCIATE

  ! Filling SIDE_ELEMID, SIDE_LOCALID, SIDE_NBSIDEID of the SideInfo_Shared-array
  DO iElem = FirstElemInd,LastElemInd
    iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem)
    SideInfo_Shared(SIDE_ELEMID,iSide+1:ElemInfo_Shared(ELEM_LASTSIDEIND,iElem)) = iElem
    sideCount = 0
    nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,iElem) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem)
    DO iLocSide = 1,nlocSides
      iSide = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem) + iLocSide
      ! Big mortar side
      IF (SideInfo_Shared(SIDE_TYPE,iSide).LE.100) THEN
        sideCount = sideCount + 1
        SideInfo_Shared(SIDE_LOCALID,iSide) = sideCount
      ELSE
        ! Mortar case
        SideInfo_Shared(SIDE_LOCALID,iSide) = -1
      END IF
      ! Check all sides on the small element side to find the small mortar side pointing back
      NbElemID    = SideInfo_Shared(SIDE_NBELEMID,iSide)
      IF(NbElemID.EQ.0) THEN
        SideInfo_Shared(SIDE_NBSIDEID,iSide) = 0
      ELSE IF (NbElemID.LT.-1) THEN
        SideInfo_Shared(SIDE_NBSIDEID,iSide) = -1
      ELSE
        nlocSidesNb = ElemInfo_Shared(ELEM_LASTSIDEIND,NbElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,NbElemID)
        DO jLocSide = 1,nlocSidesNb
          NbSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,NbElemID) + jLocSide
          IF (ABS(SideInfo_Shared(SIDE_ID,iSide)).EQ.ABS(SideInfo_Shared(SIDE_ID,NbSideID))) THEN
            SideInfo_Shared(SIDE_NBSIDEID,iSide) = NbSideID
            EXIT
          END IF
        END DO
      END IF
    END DO
  END DO

  ASSOCIATE( nNonUniqueGlobalNodes => INT(nNonUniqueGlobalNodes,IK))
    ALLOCATE(NodeInfo_Shared(1:nNonUniqueGlobalNodes))
    CALL ReadArray('GlobalNodeIDs',1,(/nNonUniqueGlobalNodes/),0_IK,1,IntegerArray_i4=NodeInfo_Shared)
    ALLOCATE(NodeCoords_Shared(3,nNonUniqueGlobalNodes))
    CALL ReadArray('NodeCoords',2,(/3_IK,nNonUniqueGlobalNodes/),0_IK,2,RealArray=NodeCoords_Shared)
  END ASSOCIATE

    ! CGNS Mapping
  CornerNodeIDswitch(1)=1
  CornerNodeIDswitch(2)=(Ngeo+1)
  CornerNodeIDswitch(3)=(Ngeo+1)**2
  CornerNodeIDswitch(4)=(Ngeo+1)*Ngeo+1
  CornerNodeIDswitch(5)=(Ngeo+1)**2*Ngeo+1
  CornerNodeIDswitch(6)=(Ngeo+1)**2*Ngeo+(Ngeo+1)
  CornerNodeIDswitch(7)=(Ngeo+1)**2*Ngeo+(Ngeo+1)**2
  CornerNodeIDswitch(8)=(Ngeo+1)**2*Ngeo+(Ngeo+1)*Ngeo+1

  ALLOCATE(ElemNodeID_Shared(        1:8,1:nElems))
  ALLOCATE(ElemSideNodeID_Shared(1:4,1:6,1:nElems))
  ALLOCATE(ElemUniqueNodeID(1:8,1:nElems))
  ALLOCATE(NodeCoords_Connect(1:3,1:nUniqueNodes))

  ASSOCIATE(CNS => CornerNodeIDswitch)
    NodeMap(:,1)=(/CNS(1),CNS(4),CNS(3),CNS(2)/)
    NodeMap(:,2)=(/CNS(1),CNS(2),CNS(6),CNS(5)/)
    NodeMap(:,3)=(/CNS(2),CNS(3),CNS(7),CNS(6)/)
    NodeMap(:,4)=(/CNS(3),CNS(4),CNS(8),CNS(7)/)
    NodeMap(:,5)=(/CNS(1),CNS(5),CNS(8),CNS(4)/)
    NodeMap(:,6)=(/CNS(5),CNS(6),CNS(7),CNS(8)/)
    DO iElem = FirstElemInd,LastElemInd
      DO iNode = 1,8
        ElemNodeID_Shared(iNode,iElem) = ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem) + CNS(iNode)
        ElemUniqueNodeID(iNode,iElem)=ABS(NodeInfo_Shared(ElemNodeID_Shared(iNode,iElem)))
        NodeCoords_Connect(1:3,ElemUniqueNodeID(iNode,iElem)) = NodeCoords_Shared(1:3,ElemNodeID_Shared(iNode,iElem))
      END DO
      nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,iElem) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem)
      DO iLocSide = 1,nlocSides
        ! Get global SideID
        GlobalSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,iElem) + iLocSide
        IF (SideInfo_Shared(SIDE_LOCALID,GlobalSideID).LE.0) CYCLE
        localSideID = SideInfo_Shared(SIDE_LOCALID,GlobalSideID)
        ! Find start of CGNS mapping from flip
        nStart = MAX(0,MOD(SideInfo_Shared(SIDE_FLIP,GlobalSideID),10)-1)
        ! Shared memory array starts at 1, but NodeID at 0
        ElemSideNodeID_Shared(1:4,localSideID,iElem) = (/ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem)+NodeMap(MOD(nStart,4)+1,localSideID)-1, &
                                                        ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem)+NodeMap(MOD(nStart+1,4)+1,localSideID)-1, &
                                                        ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem)+NodeMap(MOD(nStart+2,4)+1,localSideID)-1, &
                                                        ElemInfo_Shared(ELEM_FIRSTNODEIND,iElem)+NodeMap(MOD(nStart+3,4)+1,localSideID)-1/)
      END DO
    END DO
  END ASSOCIATE
  ElemMeshInit = .TRUE.
  CALL CloseDataFile()
END IF

IF(SurfaceDataExists.AND..NOT.SurfMeshInit) THEN
  CALL BuildSurfMeshConnectivity(InputStateFile)
  SurfMeshInit = .TRUE.
END IF

END SUBROUTINE ReadMesh_piclas2vtk


SUBROUTINE ConvertPartData(InputStateFile)
!===================================================================================================================================
! Subroutine to write 3D point data to VTK format
!===================================================================================================================================
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


SUBROUTINE ConvertElemData(InputStateFile)
!===================================================================================================================================
! Subroutine to write 3D point data to VTK format
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars    ,ONLY: ProjectName
USE MOD_IO_HDF5         ,ONLY: HSize
USE MOD_HDF5_Input      ,ONLY: OpenDataFile,CloseDataFile,ReadAttribute,File_ID,ReadArray,GetDataSize
USE MOD_Mesh_ReadIn     ,ONLY: readMesh
USE MOD_Mesh_Vars       ,ONLY: nElems, offsetElem
USE MOD_piclas2vtk_Vars ,ONLY: nUniqueNodes, NodeCoords_Connect, ElemUniqueNodeID
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
CALL GetDataSize(File_ID,'ElemData',nDims,HSize)
nVarAdd=INT(HSize(1),4)

IF (nVarAdd.GT.0) THEN
  ALLOCATE(VarNamesAdd(1:nVarAdd))
  CALL ReadAttribute(File_ID,'VarNamesAdd',nVarAdd,StrArray=VarNamesAdd(1:nVarAdd))
  ALLOCATE(ElemData(1:nVarAdd,1:nElems))

  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (&
        nVarAdd => INT(nVarAdd,IK) ,&
        offsetElem => INT(offsetElem,IK),&
        nElems     => INT(nElems,IK)    )
    CALL ReadArray('ElemData',2,(/nVarAdd, nElems/),offsetElem,2,RealArray=ElemData(1:nVarAdd,1:nElems))
  END ASSOCIATE
  SELECT CASE(TRIM(File_Type))
    CASE('State')
      FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution_ElemData',OutputTime))//'.vtu'
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


SUBROUTINE ConvertSurfaceData(InputStateFile)
!===================================================================================================================================
!
!===================================================================================================================================
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
CALL ReadAttribute(File_ID,'DSMC_nSurfSample',1,IntegerScalar=nSurfSample)
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


SUBROUTINE BuildSurfMeshConnectivity(InputStateFile)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_IO_HDF5            ,ONLY: HSize
USE MOD_HDF5_Input         ,ONLY: OpenDataFile,CloseDataFile,ReadAttribute,GetDataSize,File_ID,ReadArray,GetDataSize
USE MOD_Mesh_ReadIn        ,ONLY: readMesh
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
CHARACTER(LEN=255)              :: MeshFile
CHARACTER(LEN=255),ALLOCATABLE  :: SurfBCName_HDF5(:)
INTEGER                         :: nDims, iNode, nSurfBC_HDF5, iName, nSides, nBCs
INTEGER                         :: iElem, iLocSide, iSide, iNode2, iBC, nSurfSides, nUniqueSurfSide
INTEGER, ALLOCATABLE            :: TempBCSurfNodes(:), TempSideSurfNodeMap(:,:), SideToSurfSide(:)
REAL, ALLOCATABLE               :: TempNodeCoords(:,:)
LOGICAL                         :: IsSortedSurfNode
!===================================================================================================================================
CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile)
CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
! Read boundary names from data file
CALL GetDataSize(File_ID,'BCNames',nDims,HSize)
CHECKSAFEINT(HSize(1),4)
nBCs=INT(HSize(1),4)
ALLOCATE(BoundaryName(nBCs))

ASSOCIATE(nBCs => INT(nBCs,IK))
  CALL ReadArray('BCNames',1,(/nBCs/),0_IK,1,StrArray=BoundaryName)
END ASSOCIATE
CALL CloseDataFile()

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
ALLOCATE(SurfConnect%NonUnique2UniqueSide(1:nSides))
SurfConnect%NonUnique2UniqueSide(1:nSides) = -1
! if side is surface (BC_reflective defined in State file) then map respective surface side number
nSurfSides = 0
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
      SurfConnect%NonUnique2UniqueSide(nSurfSides) = nUniqueSurfSide
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
  IF (SurfConnect%NonUnique2UniqueSide(SideToSurfSide(iSide)).EQ.-1) CYCLE
  SurfConnect%nSurfaceBCSides = SurfConnect%nSurfaceBCSides + 1
  iElem = SideInfo_Shared(SIDE_ELEMID,iSide)
  iLocSide = SideInfo_Shared(SIDE_LOCALID,iSide)
  DO iNode2 = 1, 4
    IsSortedSurfNode = .FALSE.
    DO iNode = 1, SurfConnect%nSurfaceNode
      IF (ABS(NodeInfo_Shared(ElemSideNodeID_Shared(iNode2, iLocSide, iElem)+1)).EQ.TempBCSurfNodes(iNode)) THEN
        TempSideSurfNodeMap(iNode2,SurfConnect%nSurfaceBCSides) = iNode
        IsSortedSurfNode = .TRUE.
        EXIT
      END IF
    END DO
    IF(.NOT.IsSortedSurfNode) THEN
      SurfConnect%nSurfaceNode = SurfConnect%nSurfaceNode + 1
      TempBCSurfNodes(SurfConnect%nSurfaceNode) = ABS(NodeInfo_Shared(ElemSideNodeID_Shared(iNode2, iLocSide, iElem)+1))
      TempSideSurfNodeMap(iNode2,SurfConnect%nSurfaceBCSides) = SurfConnect%nSurfaceNode
      TempNodeCoords(1:3,SurfConnect%nSurfaceNode) = NodeCoords_Shared(1:3,ElemSideNodeID_Shared(iNode2, iLocSide, iElem)+1)
    END IF
  END DO
END DO

SDEALLOCATE(SurfConnect%BCSurfNodes)
ALLOCATE(SurfConnect%BCSurfNodes(1:SurfConnect%nSurfaceNode))
SurfConnect%BCSurfNodes(1:SurfConnect%nSurfaceNode) = TempBCSurfNodes(1:SurfConnect%nSurfaceNode)
SDEALLOCATE(SurfConnect%SideSurfNodeMap)
ALLOCATE(SurfConnect%SideSurfNodeMap(1:4,1:SurfConnect%nSurfaceBCSides))
SurfConnect%SideSurfNodeMap(1:4,1:SurfConnect%nSurfaceBCSides) = TempSideSurfNodeMap(1:4,1:SurfConnect%nSurfaceBCSides)
SDEALLOCATE(SurfConnect%NodeCoords)
ALLOCATE(SurfConnect%NodeCoords(1:3,1:SurfConnect%nSurfaceNode))
SurfConnect%NodeCoords(1:3,1:SurfConnect%nSurfaceNode) = TempNodeCoords(1:3,1:SurfConnect%nSurfaceNode)
SDEALLOCATE(TempBCSurfNodes)
SDEALLOCATE(TempSideSurfNodeMap)
SDEALLOCATE(TempNodeCoords)
SDEALLOCATE(SurfBCName_HDF5)
SDEALLOCATE(SideToSurfSide)
CALL CloseDataFile()

END SUBROUTINE BuildSurfMeshConnectivity


SUBROUTINE ConvertBGField(InputStateFile,ReadMeshFinished,NVisu,NodeTypeVisuOut)
!===================================================================================================================================
! Subroutine to write 3D point data to VTK format
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: ProjectName
USE MOD_HDF5_Input              ,ONLY: OpenDataFile,GetDataProps,CloseDataFile,ReadAttribute,File_ID,ReadArray,GetDataSize
USE MOD_Mesh_ReadIn             ,ONLY: readMesh
USE MOD_Mesh_Vars               ,ONLY: NGeo, nElems, offsetElem, NodeCoords
USE MOD_Interpolation_Vars      ,ONLY: NodeTypeVisu
USE MOD_Interpolation           ,ONLY: GetVandermonde
USE MOD_ChangeBasis             ,ONLY: ChangeBasis3D
USE MOD_VTK                     ,ONLY: WriteDataToVTK
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)              :: NVisu                             ! Polynomial degree of visualization
CHARACTER(LEN=255),INTENT(IN)   :: InputStateFile, NodeTypeVisuOut
LOGICAL,INTENT(INOUT)           :: ReadMeshFinished
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iDG, iElem
INTEGER                         :: nVar_State,N_State,nElems_State   ! Properties read from state file
CHARACTER(LEN=255)              :: NodeType_State
CHARACTER(LEN=255),ALLOCATABLE  :: StrVarNames(:)
REAL,ALLOCATABLE                :: U(:,:,:,:,:)                      ! Solution from state file
REAL,ALLOCATABLE,TARGET         :: U_Visu(:,:,:,:,:)                 ! Solution on visualiation nodes
REAL,POINTER                    :: U_Visu_p(:,:,:,:,:)               ! Solution on visualiation nodes
REAL,ALLOCATABLE                :: Coords_NVisu(:,:,:,:,:)           ! Coordinates of visualisation nodes
REAL,ALLOCATABLE,TARGET         :: Coords_BField(:,:,:,:,:)
REAL,POINTER                    :: Coords_BField_p(:,:,:,:,:)
REAL,ALLOCATABLE                :: Vdm_EQNgeo_NVisu(:,:)             ! Vandermonde from equidistand mesh to visualisation nodes
REAL,ALLOCATABLE                :: Vdm_N_NVisu(:,:)                  ! Vandermonde from state to visualisation nodes
CHARACTER(LEN=255)              :: FileString_BField, MeshFile
!===================================================================================================================================
! 1.) Open given file to get the number of elements, the order and the name of the mesh file
CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
CALL GetDataProps('BGField',nVar_State,N_State,nElems_State,NodeType_State)
CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)

SDEALLOCATE(StrVarNames)
ALLOCATE(StrVarNames(nVar_State))
CALL ReadAttribute(File_ID,'VarNames',nVar_State,StrArray=StrVarNames)

CALL CloseDataFile()

IF(.NOT.ReadMeshFinished) THEN
! Read in parameters from mesh file
  CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
  CALL ReadAttribute(File_ID,'Ngeo',1,IntegerScalar=NGeo)
  CALL CloseDataFile()
  CALL readMesh(MeshFile)
  ReadMeshFinished = .TRUE.
END IF

SDEALLOCATE(Vdm_EQNgeo_NVisu)
ALLOCATE(Vdm_EQNgeo_NVisu(0:Ngeo,0:NVisu))
CALL GetVandermonde(Ngeo,NodeTypeVisu,NVisu,NodeTypeVisuOut,Vdm_EQNgeo_NVisu,modal=.FALSE.)

SDEALLOCATE(Coords_NVisu)
ALLOCATE(Coords_NVisu(3,0:NVisu,0:NVisu,0:NVisu,nElems))
SDEALLOCATE(Coords_BField)
ALLOCATE(Coords_BField(3,0:NVisu,0:NVisu,0:NVisu,nElems))

! Convert coordinates to visu grid
DO iElem = 1,nElems
  CALL ChangeBasis3D(3,NGeo,NVisu,Vdm_EQNgeo_NVisu,NodeCoords(:,:,:,:,iElem),Coords_NVisu(:,:,:,:,iElem))
END DO

SDEALLOCATE(U)
ALLOCATE(U(nVar_State,0:N_State,0:N_State,0:N_State,nElems))

! Read in solution
CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nVar_State => INT(nVar_State,IK) ,&
      offsetElem => INT(offsetElem,IK),&
      N_State    => INT(N_State,IK),&
      nElems     => INT(nElems,IK)    )
  CALL ReadArray('BGField',5,(/nVar_State,N_State+1_IK,N_State+1_IK,N_State+1_IK,nElems/),offsetElem,5,RealArray=U)
END ASSOCIATE

CALL CloseDataFile()

SDEALLOCATE(Vdm_N_NVisu)
ALLOCATE(Vdm_N_NVisu(0:N_State,0:NVisu))
CALL GetVandermonde(N_State,NodeType_State,NVisu,NodeTypeVisuOut,Vdm_N_NVisu,modal=.FALSE.)

SDEALLOCATE(U_Visu)
ALLOCATE(U_Visu(nVar_State,0:NVisu,0:NVisu,0:NVisu,nElems))

! Interpolate solution to visu grid
iDG = 0
DO iElem = 1,nElems
  iDG = iDG + 1
  CALL ChangeBasis3D(nVar_State,N_State,NVisu,Vdm_N_NVisu,U(:,:,:,:,iElem),U_Visu(:,:,:,:,iDG))
  Coords_BField(:,:,:,:,iDG) = Coords_NVisu(:,:,:,:,iElem)
END DO

! Write solution to vtk
FileString_BField=TRIM(ProjectName)//'_BGField.vtu'
Coords_BField_p => Coords_BField(:,:,:,:,1:iDG)
U_Visu_p => U_Visu(:,:,:,:,1:iDG)
CALL WriteDataToVTK(nVar_State,NVisu,iDG,StrVarNames,Coords_BField_p,U_Visu_p,TRIM(FileString_BField),dim=3,DGFV=0)

END SUBROUTINE ConvertBGField
