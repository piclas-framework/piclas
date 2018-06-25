#include "boltzplatz.h"

!==================================================================================================================================
!> The BOLTZPLATZ2VTK tool takes state files written during runtime by BOLTZPLATZ in the .h5 format and converts them to .vtu files,
!> readable by ParaView. Supports parallel readin. 
!> The state files can come from different calculations with different mesh files, equation systems, polynomial degrees and so on.
!> Two modes of usage: command line mode and parameter file mode.
!> In parameter file mode the usage is: h5boltz2vtk parameter.ini State1.h5 State2.h5 State3.h5 ...
!> In the parameter file the following can be specified:
!> - NVisu: Integer, polynomial degree of visualization basis
!> - NodeTypeVisu: String, node type of visualization basis
!> - useCurveds: Logical, should the mesh be curved or not (if the mesh itself is curved)
!> In command line mode, only the degree of the visualization basis can be directly specified, no parameter file is needed:
!> h5boltz2vtk --NVisu=INTEGER State1.h5 State2.h5 State3.h5 ...
!> All other options are set to their standard values.
!==================================================================================================================================
PROGRAM H5BOLTZ2VTK
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_StringTools
USE MOD_Commandline_Arguments
USE MOD_IO_HDF5,             ONLY: InitIO,DefineParametersIO
USE MOD_IO_HDF5,             ONLY: HSize
USE MOD_MPI,                 ONLY: InitMPI!,DefineParametersMPI
USE MOD_ReadInTools ,        ONLY: prms,PrintDefaultParameterFile
USE MOD_ReadInTools,         ONLY: GETINT,GETSTR,GETLOGICAL
USE MOD_HDF5_Input,          ONLY: OpenDataFile,CloseDataFile,GetDataProps,ReadAttribute,File_ID,ReadArray,GetDataSize,DatasetExists
USE MOD_HDF5_Input,          ONLY: ISVALIDHDF5FILE,ISVALIDMESHFILE
USE MOD_Mesh_ReadIn,         ONLY: readMesh
USE MOD_Mesh,                ONLY: FinalizeMesh
USE MOD_Mesh_Vars,           ONLY: useCurveds,NGeo,nElems,NodeCoords,offsetElem
USE MOD_Interpolation_Vars,  ONLY: NodeTypeCL,NodeTypeVisu
USE MOD_Interpolation,       ONLY: GetVandermonde
USE MOD_ChangeBasis,         ONLY: ChangeBasis3D
USE MOD_VTK,                 ONLY: WriteDataToVTK,WriteVTKMultiBlockDataSet
#ifdef MPI
USE MOD_MPI_Vars,            ONLY: NbProc,nMPISides_Proc
#endif /*MPI*/
USE MOD_Analyze,             ONLY: CalcErrorStateFiles
USE MOD_Analyze_Vars,        ONLY: NAnalyze
USE MOD_Metrics,             ONLY: CalcMetricsErrorDiff
USE MOD_Mesh_Vars,           ONLY: Elem_xGP,dXCL_N,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,sJ,NGeoRef,DetJac_Ref
USE MOD_Mesh_Vars,           ONLY: interpolateFromTree,NormVec,nSides
USE MOD_Mesh_Vars,           ONLY: Face_xGP,NormVec,TangVec1,TangVec2,SurfElem,Ja_Face
USE MOD_PreProc,             ONLY: PP_N
USE MOD_Mesh_Vars,           ONLY: ElemToSide,SideToElem,BC,AnalyzeSide
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Time                              ! Used to track computation time  
CHARACTER(LEN=255)             :: NodeTypeVisuOut                   ! Stores user selected type of visualization nodes
INTEGER                        :: NVisu                             ! Polynomial degree of visualization
INTEGER                        :: iArgs,iElem                       ! Loop counters 
INTEGER                        :: iExt                              ! Stores position where the filename extension begins
CHARACTER(LEN=255)             :: InputStateFile,MeshFile
INTEGER                        :: nVar_State,N_State,nElems_State   ! Properties read from state file
CHARACTER(LEN=255)             :: NodeType_State                    !     "
REAL,ALLOCATABLE               :: U(:,:,:,:,:)                      ! Solution from state file
REAL,ALLOCATABLE               :: U_first(:,:,:,:,:)                ! Solution from state file
REAL,ALLOCATABLE,TARGET        :: U_Visu(:,:,:,:,:)                 ! Solution on visualiation nodes
REAL,POINTER                   :: U_Visu_p(:,:,:,:,:)               ! Solution on visualiation nodes
REAL,ALLOCATABLE               :: Coords_NVisu(:,:,:,:,:)           ! Coordinates of visualisation nodes 
REAL,ALLOCATABLE,TARGET        :: Coords_DG(:,:,:,:,:)
REAL,POINTER                   :: Coords_DG_p(:,:,:,:,:)
REAL,ALLOCATABLE               :: Vdm_EQNgeo_NVisu(:,:)             ! Vandermonde from equidistand mesh to visualisation nodes
REAL,ALLOCATABLE               :: Vdm_N_NVisu(:,:)                  ! Vandermonde from state to visualisation nodes
INTEGER                        :: nGeo_old,nVar_State_old           ! Variables used to check if we need to reinitialize
INTEGER                        :: N_State_old,nElems_old            !     "
INTEGER                        :: N_State_first                     ! first state file
CHARACTER(LEN=255)             :: MeshFile_old                      !     "
CHARACTER(LEN=255)             :: NodeType_State_old                !     "
CHARACTER(LEN=255)             :: FileString_DG
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
REAL                           :: OutputTime
INTEGER                        :: iDG
CHARACTER(LEN=255)             :: FileString_multiblock
INTEGER                        :: nDims,nVarAdd_HDF5,iVarAdd,nLocalElems_DG
LOGICAL                        :: elemDataFound
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesAdd_HDF5(:)
REAL,ALLOCATABLE               :: ElemData_HDF5(:,:)
LOGICAL                        :: CmdLineMode                       ! In command line mode only NVisu is specified directly,
                                                                    ! otherwise a parameter file is needed
CHARACTER(LEN=2)               :: NVisuString                       ! String containing NVisu from command line option
CHARACTER(LEN=20)              :: fmtString                         ! String containing options for formatted write
LOGICAL                        :: CalcDiffError                     ! Use first state file as reference state for L2 error 
                                                                    ! calculation with the following state files
CHARACTER(LEN=40)              :: DefStr
!==================================================================================================================================
CALL InitMPI()
CALL ParseCommandlineArguments()
!CALL DefineParametersMPI()
! Define parameters for H5BOLTZ2VTK
CALL prms%SetSection("H5BOLTZ2VTK")
CALL prms%CreateStringOption( 'NodeTypeVisu',"Node type of the visualization basis: "//&
                                             "VISU,GAUSS,GAUSS-LOBATTO,CHEBYSHEV-GAUSS-LOBATTO", 'VISU')
CALL prms%CreateIntOption(    'NVisu',       "Number of points at which solution is sampled for visualization.")
CALL prms%CreateLogicalOption('useCurveds',  "Controls usage of high-order information in mesh. Turn off to discard "//&
                                             "high-order data and treat curved meshes as linear meshes.", '.TRUE.')
CALL prms%CreateLogicalOption('CalcDiffError',  "Use first state file as reference state for L2 error calculation "//&
                                                "with the following state files.", '.TRUE.')
CALL prms%CreateIntOption(    'NAnalyze'         , 'Polynomial degree at which analysis is performed (e.g. for L2 errors).\n'//&
                                                   'Default: 2*N. (needed for CalcDiffError)')
CALL DefineParametersIO()

! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF

IF (nArgs.LT.2) THEN
  ! If not, print out error message containing valid syntax
  CALL CollectiveStop(__STAMP__,&
  'ERROR - Please supply a parameter file (or directly specify the NVisu option with --NVisu=INTEGER) and at least one .h5 files.')
END IF

! Measure init duration
StartTime=BOLTZPLATZTIME()

ParameterFile = Args(1)
! Check if first argument is a parameter file or that the NVisu argument has been specified
iExt=INDEX(ParameterFile,'.',BACK = .TRUE.) ! Position of file extension
IF(ParameterFile(iExt+1:iExt+3) .EQ. 'ini') THEN
  ! Parameter file has been supplied
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

SWRITE(UNIT_stdOut,'(A)') &
" ________ ___       _______      ___    ___ ___    _______  ___      ___ _________  ___  __        "
SWRITE(UNIT_stdOut,'(A)') &
"|\\  _____\\\\  \\     |\\   ___\\    |\\  \\  /  /|\\  \\  /  ___  \\|\\  \\    /  /|\\___   ___\\\\  \\|\\  \\      "
SWRITE(UNIT_stdOut,'(A)') &
"\\ \\  \\__/\\ \\  \\    \\ \\  \\__/    \\ \\  \\/  / | \\  \\/__/|_/  /\\ \\  \\  /  / ||___ \\  \\_\\ \\  \\/  /|_    "
SWRITE(UNIT_stdOut,'(A)') &
" \\ \\   __\\\\ \\  \\    \\ \\   __\\    \\ \\    / / \\ \\  \\__|//  / /\\ \\  \\/  / /     \\ \\  \\ \\ \\   ___  \\   "
SWRITE(UNIT_stdOut,'(A)') &
"  \\ \\  \\_| \\ \\  \\____\\ \\  \\_/__   /     \\/   \\ \\  \\  /  /_/__\\ \\    / /       \\ \\  \\ \\ \\  \\\\ \\  \\  "
SWRITE(UNIT_stdOut,'(A)') &
"   \\ \\__\\   \\ \\_______\\ \\______\\ /  /\\   \\    \\ \\__\\|\\________\\ \\__/ /         \\ \\__\\ \\ \\__\\\\ \\__\\ "
SWRITE(UNIT_stdOut,'(A)') &
"    \\|__|    \\|_______|\\|______|/__/ /\\ __\\    \\|__| \\|_______|\\|__|/           \\|__|  \\|__| \\|__| "
SWRITE(UNIT_stdOut,'(A)') &
"                                |__|/ \\|__|                                                        "
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(132("="))')
                                                                                                  
! Set and read in parameters differently depending if H5BOLTZ2VTK is invoked with a parameter file or not
IF (CmdLineMode) THEN
  ! Read NVisu from the first command line argument
  NVisuString = TRIM(ParameterFile(9:LEN(TRIM(ParameterFile))))
  READ(NVisuString,'(I2.1)') NVisu
  ! Since we are not reading a parameter file, some properties of the prms object need to be set
  prms%maxNameLen  = 13 ! gatheredWrite is the longest option
  prms%maxValueLen = 23 ! CHEBYSHEV-GAUSS-LOBATTO is the longest possible value
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

! Set necessary parameters for H5BOLTZ2VTK tool
! If no parameter file has been set, the standard values will be used
NodeTypeVisuOut  = GETSTR('NodeTypeVisu','VISU')    ! Node type of visualization basis
useCurveds       = GETLOGICAL('useCurveds','.FALSE.')  ! Allow curved mesh or not
CalcDiffError    = GETLOGICAL('CalcDiffError','.FALSE.')  ! Allow curved mesh or not
IF(CalcDiffError.AND.(nArgs.LT.3))THEN
  ! ! Set the default analyze polynomial degree NAnalyze to 2*(N+1) 
  ! WRITE(DefStr,'(i4)') 2*(PP_N+1)
  ! NAnalyze=GETINT('NAnalyze',DefStr) 
  ! CALL InitAnalyzeBasis(PP_N,NAnalyze,xGP,wBary)
  IF(nArgs.LT.3)CALL abort(__STAMP__,&
      'CalcDiffError needs a minimum of two state files!',iError)
END IF

! Initialization of I/O routines
CALL InitIO()

! Measure init duration
Time=BOLTZPLATZTIME()
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' INITIALIZATION DONE! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')

! Initialize an "old" state to check against - used to determine if we need to reinitialize some variables
nVar_State_old     = 0
N_State_old        = 0
MeshFile_old       = ''
nGeo_old           = 0
nElems_old         = 0
NodeType_State_old = ''

! Loop over remaining supplied .h5 files
DO iArgs = 2,nArgs
  InputStateFile = Args(iArgs)
  ! Check if the argument is a valid .h5 file
  IF(.NOT.ISVALIDHDF5FILE(InputStateFile)) THEN
    CALL CollectiveStop(__STAMP__,&
      'ERROR - Please supply only .h5 files after parameter file.')
  END IF
  
  SWRITE(UNIT_stdOut,'(132("="))')
  SWRITE(UNIT_stdOut,'(A,I3,A,I3,A)') 'Processing state ',iArgs-1,' of ',nArgs-1,'...'

  ! Open .h5 file
  CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
  ! Read in parameters from the State file
  CALL GetDataProps('DG_Solution',nVar_State,N_State,nElems_State,NodeType_State)
  CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile)
  CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)
  ! Check if we need to reallocate the var names array
  IF (nVar_State.NE.nVar_State_old) THEN
    SDEALLOCATE(StrVarNames)
    ALLOCATE(StrVarNames(nVar_State))
  END IF
  CALL ReadAttribute(File_ID,'VarNames',nVar_State,StrArray=StrVarNames)
  CALL ReadAttribute(File_ID,'Time',1,RealScalar=OutputTime)
  CALL CloseDataFile()

  ! Check if the mesh has changed
  IF (TRIM(MeshFile).NE.TRIM(MeshFile_old)) THEN
    IF(CalcDiffError.AND.(iArgs.GT.2))CALL abort(__STAMP__,&
    'CalcDiffError needs identical meshes!',iError)
    ! Check if the file is a valid mesh
    IF(.NOT.ISVALIDMESHFILE(MeshFile)) THEN
      CALL CollectiveStop(__STAMP__,&
        'ERROR - Not a valid mesh file.')
    END IF
    ! Deallocate and finalize mesh vars
    SDEALLOCATE(NodeCoords)
#ifdef MPI
    SDEALLOCATE(NbProc)
    SDEALLOCATE(nMPISides_Proc)
#endif /*MPI*/
    CALL FinalizeMesh()

    ! Read in parameters from mesh file
    CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
    CALL ReadAttribute(File_ID,'Ngeo',1,IntegerScalar=NGeo)
    CALL CloseDataFile()
  
    ! Read the mesh itself
    CALL readMesh(MeshFile) 

    ! Check if ne need to realloacte the Vandermonde from mesh to visualization
    IF (NGeo.NE.nGeo_old) THEN
      SDEALLOCATE(Vdm_EQNgeo_NVisu)
      ALLOCATE(Vdm_EQNgeo_NVisu(0:Ngeo,0:NVisu))
      CALL GetVandermonde(Ngeo,NodeTypeVisu,NVisu,NodeTypeVisuOut,Vdm_EQNgeo_NVisu,modal=.FALSE.)
    END IF

    ! Check if we need to reallocate the coordinate array
    IF (nElems.NE.nElems_old) THEN
      SDEALLOCATE(Coords_NVisu)
      ALLOCATE(Coords_NVisu(3,0:NVisu,0:NVisu,0:NVisu,nElems))
      SDEALLOCATE(Coords_DG)
      ALLOCATE(Coords_DG(3,0:NVisu,0:NVisu,0:NVisu,nElems))
    END IF

    ! Convert coordinates to visu grid
    DO iElem = 1,nElems
      CALL ChangeBasis3D(3,NGeo,NVisu,   Vdm_EQNgeo_NVisu,  NodeCoords(:,:,:,:,iElem),Coords_NVisu  (:,:,:,:,iElem))
    END DO
  END IF ! New mesh

  ! Check if we need to reallocate the solution array
  IF ((N_State.NE.N_State_old).OR.(nVar_State.NE.nVar_State_old).OR.(nElems.NE.nElems_old)) THEN
    SDEALLOCATE(U)
    ALLOCATE(U(nVar_State,0:N_State,0:N_State,0:N_State,nElems))
  END IF

  ! Read in solution
  CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
  CALL ReadArray('DG_Solution',5,(/nVar_State,N_State+1,N_State+1,N_State+1,nElems/),offsetElem,5,RealArray=U)  

  IF(CalcDiffError)THEN
    IF(iArgs <= 2)THEN
      ! Set the default analyze polynomial degree NAnalyze to 2*(N+1) 
      WRITE(DefStr,'(i4)') 2*(N_State+1)
      NAnalyze=GETINT('NAnalyze',DefStr) 
      !CALL InitAnalyzeBasis(N_State,NAnalyze,xGP,wBary)
      ! Copy state to 'first'
      ALLOCATE(U_first(nVar_State,0:N_State,0:N_State,0:N_State,nElems))
      U_first       = U
      N_State_first = N_State




      ! fill ElemToSide, SideToElem,BC
      ALLOCATE(ElemToSide(2,6,nElems))
      ALLOCATE(SideToElem(5,nSides))
      ALLOCATE(BC(1:nSides))
      ALLOCATE(AnalyzeSide(1:nSides))
      ElemToSide  = 0
      SideToElem  = -1   !mapping side to elem, sorted by side ID (for surfint)
      BC          = 0
      AnalyzeSide = 0

      interpolateFromTree=.FALSE.
      PP_N=N_State
      ! allocate type Mesh
      ! volume data
      ALLOCATE(Elem_xGP      (3,0:N_State,0:N_State,0:N_State,nElems))
      ALLOCATE(      dXCL_N(3,3,0:N_State,0:N_State,0:N_State,nElems)) ! temp
      ALLOCATE(Metrics_fTilde(3,0:N_State,0:N_State,0:N_State,nElems))
      ALLOCATE(Metrics_gTilde(3,0:N_State,0:N_State,0:N_State,nElems))
      ALLOCATE(Metrics_hTilde(3,0:N_State,0:N_State,0:N_State,nElems))
      ALLOCATE(sJ            (  0:N_State,0:N_State,0:N_State,nElems))
      NGeoRef=3*NGeo ! build jacobian at higher degree
      ALLOCATE(    DetJac_Ref(1,0:NgeoRef,0:NgeoRef,0:NgeoRef,nElems))

      ! surface data
      ALLOCATE(Face_xGP      (3,0:PP_N,0:PP_N,1:nSides))
      ALLOCATE(NormVec       (3,0:PP_N,0:PP_N,1:nSides)) 
      ALLOCATE(TangVec1      (3,0:PP_N,0:PP_N,1:nSides)) 
      ALLOCATE(TangVec2      (3,0:PP_N,0:PP_N,1:nSides))  
      ALLOCATE(SurfElem      (  0:PP_N,0:PP_N,1:nSides))  
      ALLOCATE(     Ja_Face(3,3,0:PP_N,0:PP_N,             1:nSides)) ! temp
      Face_xGP=0.
      NormVec=0.
      TangVec1=0.
      TangVec2=0.
      SurfElem=0.

      CALL CalcMetricsErrorDiff()


    ELSE
      IF(NAnalyze.LT.2*(N_State+1))CALL abort(__STAMP__,&
    'CalcDiffError: NAnalyze.LT.2*(N_State+1)! The polynomial degree is too small!',iError)
      CALL CalcErrorStateFiles(nVar_State,N_State_first,N_State,U_first,U)
    END IF 
  END IF

  ! check if additional elem data exists
  CALL DatasetExists(File_ID,'ElemData',elemDataFound)
  IF(elemDataFound)THEN
    CALL GetDataSize(File_ID,'ElemData',nDims,HSize)
    nVarAdd_HDF5=INT(HSize(1),4)
    IF (nVarAdd_HDF5.GT.0) THEN
      SDEALLOCATE(VarNamesAdd_HDF5)
      ALLOCATE(VarNamesAdd_HDF5(nVarAdd_HDF5))
      CALL ReadAttribute(File_ID,'VarNamesAdd',nVarAdd_HDF5,StrArray=VarNamesAdd_HDF5)
      SDEALLOCATE(ElemData_HDF5)
      ALLOCATE(ElemData_HDF5(1:nVarAdd_HDF5, nElems))
      CALL ReadArray('ElemData',2,(/nVarAdd_HDF5, nElems/),offsetElem,2,RealArray=ElemData_HDF5)
    END IF
  ELSE
    nVarAdd_HDF5=0
  END IF
  CALL CloseDataFile()

  ! Check if we need to reallocate the Vandermonde from state to visualisation
  IF ((N_State.NE.N_State_old).OR.(TRIM(NodeType_State).NE.TRIM(NodeType_State_old))) THEN
    SDEALLOCATE(Vdm_N_NVisu)
    ALLOCATE(Vdm_N_NVisu(0:N_State,0:NVisu))
    CALL GetVandermonde(N_State,NodeType_State,NVisu,NodeTypeVisuOut,Vdm_N_NVisu,modal=.FALSE.)
  END IF

  ! Check if we need to reallocate the visualisation array
  IF ((nVar_State.NE.nVar_State_old).OR.(nElems.NE.nElems_old)) THEN
    SDEALLOCATE(U_Visu)
    ALLOCATE(U_Visu(nVar_State,0:NVisu,0:NVisu,0:NVisu,nElems))
  END IF

  ! Interpolate solution to visu grid
  iDG = 0
  DO iElem = 1,nElems
    iDG = iDG + 1
    CALL ChangeBasis3D(nVar_State,N_State,NVisu,Vdm_N_NVisu,U(:,:,:,:,iElem),U_Visu(:,:,:,:,iDG))
    Coords_DG(:,:,:,:,iDG) = Coords_NVisu(:,:,:,:,iElem)
  END DO

  ! Write solution to vtk
  FileString_DG=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution',OutputTime))//'.vtu'
  Coords_DG_p => Coords_DG(:,:,:,:,1:iDG)
  U_Visu_p => U_Visu(:,:,:,:,1:iDG)
  CALL WriteDataToVTK(nVar_State,NVisu,iDG,StrVarNames,Coords_DG_p,U_Visu_p,TRIM(FileString_DG),dim=3,DGFV=0)
  IF (MPIRoot) THEN
    ! write multiblock file
    FileString_multiblock=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution',OutputTime))//'.vtm'
    CALL WriteVTKMultiBlockDataSet(FileString_multiblock,FileString_DG)
  END IF

  ! Save parameters of this state to later check if we need to reinitialize variables
  nVar_State_old     = nVar_State
  N_State_old        = N_State
  MeshFile_old       = MeshFile
  nGeo_old           = nGeo
  nElems_old         = nElems
  NodeType_State_old = NodeType_State
END DO ! iArgs = 3, nArgs

! Finalize
SDEALLOCATE(Vdm_N_NVisu)
SDEALLOCATE(Vdm_EQNgeo_NVisu)
SDEALLOCATE(U)
SDEALLOCATE(U_Visu)
SDEALLOCATE(Coords_NVisu)
SDEALLOCATE(NodeCoords)
CALL FinalizeMesh()

! Measure processing duration
Time=BOLTZPLATZTIME()
#ifdef MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) THEN
  CALL abort(__STAMP__,&
    'MPI finalize error',iError)
END IF
#endif
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' H5BOLTZ2VTK FINISHED! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')

END PROGRAM H5BOLTZ2VTK

