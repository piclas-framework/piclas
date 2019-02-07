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
!> In parameter file mode the usage is: h5piclas2vtk parameter.ini State1.h5 State2.h5 State3.h5 ...
!> In the parameter file the following can be specified:
!> - NVisu: Integer, polynomial degree of visualization basis
!> - NodeTypeVisu: String, node type of visualization basis
!> - useCurveds: Logical, should the mesh be curved or not (if the mesh itself is curved)
!> In command line mode, only the degree of the visualization basis can be directly specified, no parameter file is needed:
!> h5piclas2vtk --NVisu=INTEGER State1.h5 State2.h5 State3.h5 ...
!> All other options are set to their standard values.
!==================================================================================================================================
PROGRAM H5PICLAS2VTK
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
USE MOD_Mesh_Vars,           ONLY: useCurveds,NGeo,nElems,NodeCoords,offsetElem, Elems, nNodes
USE MOD_Particle_Mesh_Vars,  ONLY: GEO
USE MOD_Interpolation_Vars,  ONLY: NodeTypeCL,NodeTypeVisu
USE MOD_Interpolation,       ONLY: GetVandermonde
USE MOD_ChangeBasis,         ONLY: ChangeBasis3D
USE MOD_VTK,                 ONLY: WriteDataToVTK,WriteVTKMultiBlockDataSet
USE MOD_Prepare_Mesh,        ONLY: fillMeshInfo
#ifdef MPI
USE MOD_MPI_Vars,            ONLY: NbProc,nMPISides_Proc
#endif /*MPI*/
USE MOD_Analyze,             ONLY: CalcErrorStateFiles, CalcErrorStateFileSigma
USE MOD_Analyze_Vars,        ONLY: NAnalyze
USE MOD_Mesh_Vars,           ONLY: sJ,NGeoRef
USE MOD_PreProc,             ONLY: PP_N
USE MOD_Metrics,             ONLY: CalcMetricsErrorDiff
USE MOD_Particle_Boundary_Vars, ONLY:SurfMeshPosti
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
REAL,ALLOCATABLE               :: U_average(:,:,:,:,:)              ! Solution from state file
INTEGER                        :: N_average
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
CHARACTER(LEN=255)             :: FileString_DG, FileString, File_Type
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
REAL                           :: OutputTime
INTEGER                        :: iDG
CHARACTER(LEN=255)             :: FileString_multiblock
INTEGER                        :: nDims,nVarAdd_HDF5
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesAdd_HDF5(:)
REAL,ALLOCATABLE               :: ElemData_HDF5(:,:)
LOGICAL                        :: CmdLineMode, NVisuDefault         ! In command line mode only NVisu is specified directly,
                                                                    ! otherwise a parameter file is needed
CHARACTER(LEN=2)               :: NVisuString                       ! String containing NVisu from command line option
CHARACTER(LEN=20)              :: fmtString                         ! String containing options for formatted write
LOGICAL                        :: CalcDiffError                     ! Use first state file as reference state for L2 error 
                                                                    ! calculation with the following state files
LOGICAL                        :: AllowChangedMesh
LOGICAL                        :: CalcDiffSigma                     ! Use last state file as state for L2 sigma calculation
LOGICAL                        :: CalcAverage                       ! Calculate and write arithmetic mean of alle StateFile
LOGICAL                        :: VisuSource, DGSourceExists, skip, DGSolutionExists, ElemDataExists, SurfaceDataExists
CHARACTER(LEN=40)              :: DefStr
INTEGER                        :: iArgsStart, iNode, jNode
! Surface
INTEGER                        :: nVarSurf_HDF5, nSurfBC_HDF5, iName, nSurfSample
CHARACTER(LEN=255),ALLOCATABLE :: SurfBCName_HDF5(:), VarNamesSurf_HDF5(:)
REAL, ALLOCATABLE              :: SurfData_HDF5(:,:,:,:)
LOGICAL                        :: MeshInitFinished
!==================================================================================================================================
CALL InitMPI()
CALL ParseCommandlineArguments()
!CALL DefineParametersMPI()
!CALL DefineParametersIO_HDF5()
! Define parameters for H5PICLAS2VTK
CALL prms%SetSection("H5PICLAS2VTK")
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
                                                                                                  
! Set and read in parameters differently depending if H5PICLAS2VTK is invoked with a parameter file or not
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

! Set necessary parameters for H5PICLAS2VTK tool
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
VisuSource    = GETLOGICAL('VisuSource','.FALSE.')

! Initialization of I/O routines
CALL InitIO()

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

MeshInitFinished = .FALSE.

! Loop over remaining supplied .h5 files
DO iArgs = iArgsStart,nArgs
  InputStateFile = Args(iArgs)
  ! Check if the argument is a valid .h5 file
  IF(.NOT.ISVALIDHDF5FILE(InputStateFile)) THEN
    CALL CollectiveStop(__STAMP__,&
      'ERROR - Please supply only .h5 files after parameter file.')
  END IF
  
  SWRITE(UNIT_stdOut,'(132("="))')
  SWRITE(UNIT_stdOut,'(A,I3,A,I3,A)') 'Processing state ',iArgs-1,' of ',nArgs-1,'...'

  ! Open .h5 file
  DGSolutionExists = .FALSE.; ElemDataExists = .FALSE.; SurfaceDataExists = .FALSE.
  CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
  CALL DatasetExists(File_ID,'DG_Solution',DGSolutionExists)
  CALL DatasetExists(File_ID,'ElemData',ElemDataExists)
  CALL DatasetExists(File_ID,'SurfaceData',SurfaceDataExists)

  ! === DG_Solution ================================================================================================================
  ! Read in parameters from the State file
  IF(DGSolutionExists) THEN
    CALL GetDataProps('DG_Solution',nVar_State,N_State,nElems_State,NodeType_State)
    CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile)
    CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)

    IF (VisuSource) THEN
      CALL DatasetExists(File_ID,'DG_Source',DGSourceExists)
    ELSE
      DGSourceExists=.FALSE.
    END IF
    IF (DGSourceExists) THEN
      nVar_State=4
      ! Check if we need to reallocate the var names array
      IF (nVar_State.NE.nVar_State_old) THEN
        SDEALLOCATE(StrVarNames)
        ALLOCATE(StrVarNames(nVar_State))
      END IF
      CALL ReadAttribute(File_ID,'VarNamesSource',nVar_State,StrArray=StrVarNames)
    ELSE
      VisuSource=.FALSE.
      ! Check if we need to reallocate the var names array
      IF (nVar_State.NE.nVar_State_old) THEN
        SDEALLOCATE(StrVarNames)
        ALLOCATE(StrVarNames(nVar_State))
      END IF
      CALL ReadAttribute(File_ID,'VarNames',nVar_State,StrArray=StrVarNames)
    END IF
    CALL ReadAttribute(File_ID,'Time',1,RealScalar=OutputTime)
    CALL CloseDataFile()

    ! Check if the mesh has changed
    IF(CalcDiffError.AND.(iArgs.GT.2).AND.AllowChangedMesh)THEN
      skip=.TRUE.
    ELSE
      skip=.FALSE.
    END IF
    IF (.NOT.skip .AND. (TRIM(MeshFile).NE.TRIM(MeshFile_old))) THEN
      IF(CalcDiffError.AND.(iArgs.GT.2))CALL abort(__STAMP__,&
      'CalcDiffError needs identical meshes!',iError)
      IF(CalcAverage.AND.(iArgs.GT.2))CALL abort(__STAMP__,&
      'CalcAverage needs identical meshes!',iError)
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

      IF(CalcAverage .AND. iArgs.GT.2)THEN
        IF (nVar_State_old     .NE. nVar_State .OR.&
            N_State_old        .NE. N_State .OR.&
            nGeo_old           .NE. nGeo .OR.&
            nElems_old         .NE. nElems .OR.&
            NodeType_State_old .NE. NodeType_State) CALL abort(__STAMP__,&
        'CalcAverage: different N,nGeo,etc. not implemented yet!')
      END IF

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
        CALL ChangeBasis3D(3,NGeo,NVisu,Vdm_EQNgeo_NVisu,NodeCoords(:,:,:,:,iElem),Coords_NVisu(:,:,:,:,iElem))
      END DO
    END IF ! New mesh

    ! Check if we need to reallocate the solution array
    IF ((N_State.NE.N_State_old).OR.(nVar_State.NE.nVar_State_old).OR.(nElems.NE.nElems_old)) THEN
      SDEALLOCATE(U)
      ALLOCATE(U(nVar_State,0:N_State,0:N_State,0:N_State,nElems))
    END IF

    ! Read in solution
    CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)

    ! Associate construct for integer KIND=8 possibility
    ASSOCIATE (&
          nVar_State => INT(nVar_State,IK) ,&
          offsetElem => INT(offsetElem,IK),&
          N_State    => INT(N_State,IK),&
          nElems     => INT(nElems,IK)    )
      CALL ReadArray('DG_Solution',5,(/nVar_State,N_State+1_IK,N_State+1_IK,N_State+1_IK,nElems/),offsetElem,5,RealArray=U)  
    END ASSOCIATE

    IF(VisuSource)THEN
      CALL ReadArray('DG_Source',5,(/nVar_State,N_State+1,N_State+1,N_State+1,nElems/),offsetElem,5,RealArray=U)  
    ELSE
      CALL ReadArray('DG_Solution',5,(/nVar_State,N_State+1,N_State+1,N_State+1,nElems/),offsetElem,5,RealArray=U)  
    END IF

    IF(CalcDiffError)THEN
      IF(iArgs .EQ. 2)THEN
        ! Set the default analyze polynomial degree NAnalyze to 2*(N+1) 
        WRITE(DefStr,'(i4)') 2*(N_State+1)
        NAnalyze=GETINT('NAnalyze',DefStr) 
        !CALL InitAnalyzeBasis(N_State,NAnalyze,xGP,wBary)
        ! Copy state to 'first'
        ALLOCATE(U_first(nVar_State,0:N_State,0:N_State,0:N_State,nElems))
        U_first       = U
        N_State_first = N_State

        PP_N=N_State
        NGeoRef=3*NGeo ! build jacobian at higher degree
        ALLOCATE(sJ            (  0:N_State,0:N_State,0:N_State,nElems))
        CALL CalcMetricsErrorDiff()
      ELSE IF(iArgs .EQ. nArgs .AND. CalcDiffSigma) THEN
        IF(NAnalyze.LT.2*(N_State+1))CALL abort(__STAMP__,&
      'CalcDiffError: NAnalyze.LT.2*(N_State+1)! The polynomial degree is too small!',iError)
        CALL CalcErrorStateFileSigma(nVar_State,N_State,U)

      ELSE
        IF(NAnalyze.LT.2*(N_State+1))CALL abort(__STAMP__,&
      'CalcDiffError: NAnalyze.LT.2*(N_State+1)! The polynomial degree is too small!',iError)
        CALL CalcErrorStateFiles(nVar_State,N_State_first,N_State,U_first,U)
      END IF 
    END IF
    IF(CalcAverage)THEN
      IF(iArgs .EQ. 2)THEN
        ALLOCATE(U_average(nVar_State,0:N_State,0:N_State,0:N_State,nElems))
        U_average = U
        N_average = 1
      ELSE
        U_average = U_average+U
        N_average = N_average+1
      END IF
      IF (iArgs .EQ. nArgs) THEN
        U = U_average/REAL(N_average)
      END IF
    END IF

    CALL CloseDataFile()

    IF(CalcAverage .AND. (iArgs.GT.2 .AND. iArgs.LT.nArgs)) THEN
      CYCLE !go to next file (only output the averaged U, but allocate and set _old-stuff for iArg=2)
    END IF

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

    IF(iArgs.EQ.nArgs .OR. .NOT.CalcAverage) THEN
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
    END IF !iArgs.EQ.nArgs .OR. .NOT.CalcAverage

    ! Save parameters of this state to later check if we need to reinitialize variables
    nVar_State_old     = nVar_State
    N_State_old        = N_State
    MeshFile_old       = MeshFile
    nGeo_old           = nGeo
    nElems_old         = nElems
    NodeType_State_old = NodeType_State
  END IF
  
  ! === ElemData ===================================================================================================================
  IF(ElemDataExists) THEN
    ! DSMC/FV Solution
    CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile)
    CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)
    CALL ReadAttribute(File_ID,'Time',1,RealScalar=OutputTime)
    CALL ReadAttribute(File_ID,'File_Type',1,StrScalar=File_Type)

    IF(.NOT.MeshInitFinished) THEN
      CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
      CALL ReadAttribute(File_ID,'Ngeo',1,IntegerScalar=NGeo)
      CALL CloseDataFile()
      CALL readMesh(MeshFile)
      CALL InitMesh_Connected()
      MeshInitFinished = .TRUE.
    END IF

    ! Read in solution
    CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)

    CALL GetDataSize(File_ID,'ElemData',nDims,HSize)
    nVarAdd_HDF5=INT(HSize(1),4)
    IF (nVarAdd_HDF5.GT.0) THEN
      SDEALLOCATE(VarNamesAdd_HDF5)
      ALLOCATE(VarNamesAdd_HDF5(nVarAdd_HDF5))
      CALL ReadAttribute(File_ID,'VarNamesAdd',nVarAdd_HDF5,StrArray=VarNamesAdd_HDF5)
      SDEALLOCATE(ElemData_HDF5)
      ALLOCATE(ElemData_HDF5(1:nVarAdd_HDF5, nElems))

      ! Associate construct for integer KIND=8 possibility
      ASSOCIATE (&
            nVarAdd_HDF5 => INT(nVarAdd_HDF5,IK) ,&
            offsetElem => INT(offsetElem,IK),&
            nElems     => INT(nElems,IK)    )
        CALL ReadArray('ElemData',2,(/nVarAdd_HDF5, nElems/),offsetElem,2,RealArray=ElemData_HDF5)
      END ASSOCIATE
      SELECT CASE(TRIM(File_Type))
        CASE('State')
          FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution_ElemData',OutputTime))//'.vtu'
        CASE('DSMCHOState')
          FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_DSMCState',OutputTime))//'.vtu'
      END SELECT
      CALL WriteDataToVTK3DElems(nVarAdd_HDF5,ElemData_HDF5,FileString,VarNamesAdd_HDF5)
    END IF
    CALL CloseDataFile()
  END IF
  ! === SurfaceData ================================================================================================================
  IF(SurfaceDataExists) THEN
    ! DSMC/FV Solution
    CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile)
    CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)
    CALL ReadAttribute(File_ID,'Time',1,RealScalar=OutputTime)
    CALL ReadAttribute(File_ID,'File_Type',1,StrScalar=File_Type)
    CALL GetDataProps('SurfaceData',nVarSurf_HDF5,N_State,nElems_State,NodeType_State)
    CALL ReadAttribute(File_ID,'DSMC_nSurfSample',1,IntegerScalar=nSurfSample)
    IF(nSurfSample.EQ.1) THEN
      ! Read-in of mesh information (if not already done for DG solution -> DSMCState case)
      IF(.NOT.MeshInitFinished) THEN
        ! Read in parameters from mesh file
        CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
        CALL ReadAttribute(File_ID,'Ngeo',1,IntegerScalar=NGeo)
        CALL CloseDataFile()
        CALL readMesh(MeshFile)
        CALL InitMesh_Connected()
        MeshInitFinished = .TRUE.
      END IF
    ELSE IF(nSurfSample.GT.1) THEN
      CALL abort(__STAMP__,&
        'DSMC_nSurfSample greater one is not supported yet!')
    ELSE
      CALL abort(__STAMP__,&
        'DSMC_nSurfSample is zero, something might be wrong with your .h5 file!')
    END IF
    ! Read in solution
    CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)

    SWRITE(UNIT_stdOut,'(A,A)')' GET NUMBER AND NAMES OF SURFACE-BCSIDES IN HDF5 FILE... '
    CALL GetDataSize(File_ID,'BC_Surf',nDims,HSize,attrib=.true.)
    nSurfBC_HDF5 = INT(HSize(1),4)

    SWRITE(UNIT_stdOut,'(A3,A45,A3,I33,A13)')' | ','Number of Surface BCs',' | ',nSurfBC_HDF5,' | HDF5    | '
    SDEALLOCATE(SurfBCName_HDF5)
    ALLOCATE(SurfBCName_HDF5(nSurfBC_HDF5))
    CALL ReadAttribute(File_ID,'BC_Surf',nSurfBC_HDF5,StrArray=SurfBCName_HDF5)
    DO iName = 1,nSurfBC_HDF5
      SWRITE(UNIT_stdOut,'(A3,A38,I2.1,A5,A3,A33,A13)')' | ','BC',iName,'Name',' | ',TRIM(SurfBCName_HDF5(iName)),' | HDF5    | '
    END DO
    SWRITE(UNIT_stdOut,'(A)')' DONE!'

    SDEALLOCATE(VarNamesSurf_HDF5)
    ALLOCATE(VarNamesSurf_HDF5(nVarSurf_HDF5))
    CALL ReadAttribute(File_ID,'VarNamesSurface',nVarSurf_HDF5,StrArray=VarNamesSurf_HDF5(1:nVarSurf_HDF5))

    IF((nVarSurf_HDF5.GT.0).AND.(SurfMeshPosti%nSurfaceBCSides.GT.0))THEN
      SDEALLOCATE(SurfData_HDF5)
      ALLOCATE(SurfData_HDF5(1:nVarSurf_HDF5,nSurfSample,nSurfSample,1:SurfMeshPosti%nSurfaceBCSides))
      SurfData_HDF5=0.
      CALL ReadArray('SurfaceData',4,(/nVarSurf_HDF5, nSurfSample, nSurfSample, SurfMeshPosti%nSurfaceBCSides/), &
                      0,4,RealArray=SurfData_HDF5(:,nSurfSample,nSurfSample,:))
    END IF

    FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_SurfaceOutput',OutputTime))//'.vtu'

    CALL WriteDataToVTK2D(nVarSurf_HDF5,SurfData_HDF5,FileString,VarNamesSurf_HDF5,nSurfSample)
    CALL CloseDataFile()
  END IF
END DO ! iArgs = 2, nArgs

! Finalize
SDEALLOCATE(Vdm_N_NVisu)
SDEALLOCATE(Vdm_EQNgeo_NVisu)
SDEALLOCATE(U)
SDEALLOCATE(U_Visu)
SDEALLOCATE(Coords_NVisu)
SDEALLOCATE(NodeCoords)
CALL FinalizeMesh()

! Measure processing duration
Time=PICLASTIME()
#ifdef MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) THEN
  CALL abort(__STAMP__,&
    'MPI finalize error',iError)
END IF
#endif
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F14.2,A)') ' H5PICLAS2VTK FINISHED! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')

END PROGRAM H5PICLAS2VTK

SUBROUTINE WriteDataToVTK3DElems(nVal,ValueCells,FileString,VarNameVisu)
!===================================================================================================================================
! Subroutine to write 3D point data to VTK format
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Mesh_Vars,     ONLY:GEO
USE MOD_Mesh_Vars,              ONLY:nElems,nNodes
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: nVal                                       ! Number of nodal output variables
REAL,INTENT(IN)               :: ValueCells(nVal,nElems)                    ! Statevector 
CHARACTER(LEN=*),INTENT(IN)   :: FileString, VarNameVisu(nVal)              ! Output file name
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iVal,iElem,Offset,nBytes,nVTKElems,nVTKCells,ivtk=44,iVar,iNode
INTEGER            :: int_size
INTEGER            :: Vertex(8,nElems)
INTEGER            :: ElemType
CHARACTER(LEN=35)  :: StrOffset,TempStr1,TempStr2
CHARACTER(LEN=200) :: Buffer
CHARACTER(LEN=1)   :: lf
REAL(KIND=4)       :: float
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')"   WRITE 3D DATA TO VTX XML BINARY (VTU) FILE..."
IF(nElems.LT.1)THEN
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
  RETURN
END IF

! Line feed character
lf = char(10)

! Write file
OPEN(UNIT=ivtk,FILE=TRIM(FileString),ACCESS='STREAM')
! Write header
Buffer='<?xml version="1.0"?>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify file type
nVTKElems=nNodes
nVTKCells=nElems

Buffer='  <UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
WRITE(TempStr1,'(I16)')nVTKElems
WRITE(TempStr2,'(I16)')nVTKCells
Buffer='    <Piece NumberOfPoints="'//TRIM(ADJUSTL(TempStr1))//&
'" NumberOfCells="'//TRIM(ADJUSTL(TempStr2))//'">'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify point data
Buffer='      <PointData> </PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify cell data
Buffer='      <CellData>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=0
WRITE(StrOffset,'(I16)')Offset
IF (nVal .GT.0)THEN
  DO iVar=1,nVal
      Buffer='        <DataArray type="Float32" Name="'//TRIM(VarNameVisu(iVar))//&
      '" NumberOfComponents="1" format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
      Offset=Offset+INT(SIZEOF(int_size),4)+nVTKCells*INT(SIZEOF(float),4)
      WRITE(StrOffset,'(I16)')Offset
  END DO
END IF
Buffer='      </CellData>'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify coordinate data
Buffer='      <Points>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='        <DataArray type="Float32" Name="Coordinates" NumberOfComponents="3" format="appended"'// &
       ' offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=Offset+INT(SIZEOF(int_size),4)+3*nVTKElems*INT(SIZEOF(float),4)
WRITE(StrOffset,'(I16)')Offset
Buffer='      </Points>'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify necessary cell data
Buffer='      <Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
! Connectivity
Buffer='        <DataArray type="Int32" Name="connectivity" format="appended"'// &
       ' offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=Offset+INT(SIZEOF(int_size),4)+8*nVTKCells*INT(SIZEOF(int_size),4)
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
DO iVal=1,nVal
    WRITE(ivtk) nBytes,REAL(ValueCells(iVal,1:nElems),4)
END DO
! Points
nBytes = 3*nVTKElems*INT(SIZEOF(FLOAT),4)
WRITE(ivtk) nBytes
WRITE(ivtk) REAL(GEO%NodeCoords(1:3,1:nNodes),4)
! Connectivity
DO iElem=1,nElems
  DO iNode=1,8
    Vertex(iNode,iElem)= GEO%ElemToNodeID(iNode,iElem)-1
  END DO
END DO
nBytes = 8*nVTKCells*INT(SIZEOF(int_size),4)
WRITE(ivtk) nBytes
WRITE(ivtk) Vertex(:,:)
! Offset
nBytes = nVTKCells*INT(SIZEOF(int_size),4)
WRITE(ivtk) nBytes
WRITE(ivtk) (Offset,Offset=8,8*nVTKCells,8)
! Elem type
ElemType = 12 ! VTK_HEXAHEDRON
WRITE(ivtk) nBytes
WRITE(ivtk) (ElemType,iElem=1,nVTKCells)
! Write footer
Buffer=lf//'  </AppendedData>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
CLOSE(ivtk)
SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END SUBROUTINE WriteDataToVTK3DElems


SUBROUTINE WriteDataToVTK2D(nVal,Value,FileString,VarNameVisu,nSurfSample)
!===================================================================================================================================
! Subroutine to write 2D point data to VTK format
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Mesh_Vars,     ONLY:GEO
USE MOD_Particle_Boundary_Vars, ONLY:SurfMeshPosti
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: nVal,nSurfSample                 ! Number of nodal output variables
REAL,INTENT(IN)               :: Value(nVal,nSurfSample,nSurfSample,SurfMeshPosti%nSurfaceBCSides) ! Statevector 
CHARACTER(LEN=*),INTENT(IN)   :: FileString, VarNameVisu(nVal)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElem,iVal,Offset,nBytes,nVTKElems,nVTKCells,ivtk=45,iVar,iNode
INTEGER            :: int_size
INTEGER            :: Vertex(4,SurfMeshPosti%nSurfaceBCSides)
INTEGER            :: ElemType
CHARACTER(LEN=35)  :: StrOffset,TempStr1,TempStr2
CHARACTER(LEN=200) :: Buffer
CHARACTER(LEN=1)   :: lf
REAL(KIND=4)       :: Float
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')"   WRITE 2D DATA TO VTK XML BINARY (VTU) FILE..."

! Line feed character
lf = char(10)

! Write file
OPEN(UNIT=ivtk,FILE=TRIM(FileString),ACCESS='STREAM')
! Write header
Buffer='<?xml version="1.0"?>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify file type
nVTKElems=SurfMeshPosti%nSurfaceNode
nVTKCells=SurfMeshPosti%nSurfaceBCSides

Buffer='  <UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
WRITE(TempStr1,'(I16)')nVTKElems
WRITE(TempStr2,'(I16)')nVTKCells
Buffer='    <Piece NumberOfPoints="'//TRIM(ADJUSTL(TempStr1))//&
'" NumberOfCells="'//TRIM(ADJUSTL(TempStr2))//'">'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify point data
Buffer='      <PointData> </PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify cell data
Buffer='      <CellData>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=0
WRITE(StrOffset,'(I16)')Offset
DO iVar=1,nVal
    Buffer='        <DataArray type="Float32" Name="'//TRIM(VarNameVisu(iVar))//&
    '" NumberOfComponents="1" format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
      Offset=Offset+INT(SIZEOF(int_size),4)+nVTKCells*INT(SIZEOF(Float),4)
    WRITE(StrOffset,'(I16)')Offset
END DO
Buffer='      </CellData>'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify coordinate data
Buffer='      <Points>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='        <DataArray type="Float32" Name="Coordinates" NumberOfComponents="3" format="appended"'// &
       ' offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
   Offset=Offset+INT(SIZEOF(int_size),4)+3*nVTKElems*INT(SIZEOF(Float),4)
WRITE(StrOffset,'(I16)')Offset
Buffer='      </Points>'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify necessary cell data
Buffer='      <Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
! Connectivity
Buffer='        <DataArray type="Int32" Name="connectivity" format="appended"'// &
       ' offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=Offset+INT(SIZEOF(int_size),4)+4*nVTKCells*INT(SIZEOF(int_size),4)
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
! Point data
nBytes = nVTKCells*INT(SIZEOF(Float),4)
DO iVal=1,nVal
  WRITE(ivtk) nBytes,REAL(Value(iVal,1:nSurfSample,1:nSurfSample,1:SurfMeshPosti%nSurfaceBCSides),4)
END DO
! Points
nBytes = 3*nVTKElems*INT(SIZEOF(Float),4)
WRITE(ivtk) nBytes
DO iNode=1, SurfMeshPosti%nSurfaceNode
  WRITE(ivtk) REAL(GEO%NodeCoords(1:3, SurfMeshPosti%BCSurfNodes(iNode)),4)
END DO
! Connectivity
DO iElem=1,SurfMeshPosti%nSurfaceBCSides
  DO iNode=1,4
    Vertex(iNode,iElem) = SurfMeshPosti%SideSurfNodeMap(iNode,iElem) - 1
  END DO
END DO
nBytes = 4*nVTKCells*INT(SIZEOF(int_size),4)
WRITE(ivtk) nBytes
WRITE(ivtk) Vertex(:,:)
! Offset
nBytes = nVTKCells*INT(SIZEOF(int_size),4)
WRITE(ivtk) nBytes
WRITE(ivtk) (Offset,Offset=4,4*nVTKCells,4)
! Elem type
ElemType = 9 ! VTK_QUAD
WRITE(ivtk) nBytes
WRITE(ivtk) (ElemType,iElem=1,nVTKCells)
! Write footer
Buffer=lf//'  </AppendedData>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
CLOSE(ivtk)
SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END SUBROUTINE WriteDataToVTK2D


SUBROUTINE InitMesh_Connected()
!===================================================================================================================================
! Allocate and generate surface mesh on CL_NGeo points
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
USE MOD_Particle_Boundary_Vars, ONLY:SurfMeshPosti
USE MOD_Particle_Mesh_Vars,     ONLY:GEO
!--------------------------------------------------------------------------------------------------!
! perform Mapping for Surface Output
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE 
! LOCAL VARIABLES
INTEGER                 :: iElem, iLocSide, iSide, iNode, iNode2, iBCSide, iInnerSide, nStart, jNode
INTEGER, ALLOCATABLE    :: TempBCSurfNodes(:), TempSideSurfNodeMap(:,:)
LOGICAL                 :: IsSortedSurfNode
INTEGER                 :: NodeMap(4,6)
! LOCAL VARIABLES
TYPE(tElem),POINTER :: aElem
TYPE(tSide),POINTER :: aSide
INTEGER             :: LocSideID,nSides_flip(0:4)
!===================================================================================================================================

! set side ID, so that BC Sides come first, then InnerSides
DO iElem=1,nElems
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
    aElem%Side(iLocSide)%sp%sideID=-1
  END DO
END DO
iSide=0
iBCSide=0
iInnerSide=nBCSides
DO iElem=1,nElems
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    IF(aSide%nMortars.GT.0) CALL abort(__STAMP__,&
      'Surface results with connectivity on meshes with mortars are not supported!') 
    IF(aSide%sideID.EQ.-1)THEN
      IF(aSide%NbProc.EQ.-1)THEN ! no MPI Sides
        IF(ASSOCIATED(aSide%connection))THEN
          iInnerSide=iInnerSide+1
          iSide=iSide+1
          aSide%SideID=iInnerSide
          aSide%connection%SideID=iInnerSide
        ELSE
          iBCSide=iBCSide+1
          iSide=iSide+1
          aSide%SideID=iBCSide
        END IF !associated connection
      END IF ! .NOT. MPISide
    END IF !sideID NE -1
  END DO ! iLocSide=1,6
END DO !iElem
IF(iSide.NE.(nInnerSides+nBCSides)) CALL abort(__STAMP__,'not all SideIDs are set!')

SDEALLOCATE(ElemToSide)
SDEALLOCATE(SideToElem)
SDEALLOCATE(BC)
ALLOCATE(ElemToSide(2,6,nElems))
ALLOCATE(SideToElem(5,nSides))
ALLOCATE(BC(1:nSides))
ElemToSide  = 0
SideToElem  = -1   !mapping side to elem, sorted by side ID (for surfint)
BC          = 0
! ELement to Side mapping
nSides_flip=0
DO iElem=1,nElems
  aElem=>Elems(iElem)%ep
  DO LocSideID=1,6
    aSide=>aElem%Side(LocSideID)%sp
    ElemToSide(E2S_SIDE_ID,LocSideID,iElem)=aSide%SideID
    ElemToSide(E2S_FLIP,LocSideID,iElem)   =aSide%Flip
    nSides_flip(aSide%flip)=nSides_flip(aSide%flip)+1
  END DO ! LocSideID
END DO ! iElem

! Side to Element mapping, sorted by SideID
DO iElem=1,nElems
  aElem=>Elems(iElem)%ep
  DO LocSideID=1,6
    aSide=>aElem%Side(LocSideID)%sp
    IF(aSide%Flip.EQ.0)THEN !root side
      SideToElem(S2E_ELEM_ID,aSide%SideID)         = iElem !root Element
      SideToElem(S2E_LOC_SIDE_ID,aSide%SideID)     = LocSideID
    ELSE
      SideToElem(S2E_NB_ELEM_ID,aSide%SideID)      = iElem ! element with flipped side
      SideToElem(S2E_NB_LOC_SIDE_ID,aSide%SideID)  = LocSideID
      SideToElem(S2E_FLIP,aSide%SideID)            = aSide%Flip
    END IF
    IF(aSide%sideID .LE. nBCSides) BC(aSide%sideID)=aSide%BCIndex
  END DO ! LocSideID
END DO ! iElem

! Building element to node id and node coords arrays (copied from the first part of the InitParticleGeometry)

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE GEOMETRY INFORMATION...'
NodeMap(:,1)=(/1,4,3,2/)
NodeMap(:,2)=(/1,2,6,5/)
NodeMap(:,3)=(/2,3,7,6/)
NodeMap(:,4)=(/3,4,8,7/)
NodeMap(:,5)=(/1,5,8,4/)
NodeMap(:,6)=(/5,6,7,8/)
SDEALLOCATE(GEO%ElemToNodeID)
SDEALLOCATE(GEO%ElemSideNodeID)
SDEALLOCATE(GEO%NodeCoords)
ALLOCATE(GEO%ElemToNodeID(1:8,1:nElems),GEO%ElemSideNodeID(1:4,1:6,1:nElems),GEO%NodeCoords(1:3,1:nNodes))
GEO%ElemToNodeID(:,:)=0
GEO%ElemSideNodeID(:,:,:)=0
GEO%NodeCoords(:,:)=0.
iNode=0
DO iElem=1,nElems
  DO jNode=1,8
    Elems(iElem+offsetElem)%ep%node(jNode)%np%NodeID=0
  END DO
END DO
DO iElem=1,nElems
  !--- Save corners of sides
  DO jNode=1,8
    IF (Elems(iElem+offsetElem)%ep%node(jNode)%np%NodeID.EQ.0) THEN
      iNode=iNode+1
      Elems(iElem+offsetElem)%ep%node(jNode)%np%NodeID=iNode
      GEO%NodeCoords(1:3,iNode)=Elems(iElem+offsetElem)%ep%node(jNode)%np%x(1:3)
    END IF
    GEO%ElemToNodeID(jNode,iElem)=Elems(iElem+offsetElem)%ep%node(jNode)%np%NodeID
    !GEO%ElemToNodeIDGlobal(jNode,iElem) = Elems(iElem+offsetElem)%ep%node(jNode)%np%ind
  END DO
END DO

DO iElem=1,nElems
  DO iLocSide=1,6
    nStart=MAX(0,ElemToSide(E2S_FLIP,iLocSide,iElem)-1)
    GEO%ElemSideNodeID(1:4,iLocSide,iElem)=(/Elems(iElem+offsetElem)%ep%node(NodeMap(MOD(nStart  ,4)+1,iLocSide))%np%NodeID,&
                                             Elems(iElem+offsetElem)%ep%node(NodeMap(MOD(nStart+1,4)+1,iLocSide))%np%NodeID,&
                                             Elems(iElem+offsetElem)%ep%node(NodeMap(MOD(nStart+2,4)+1,iLocSide))%np%NodeID,&
                                             Elems(iElem+offsetElem)%ep%node(NodeMap(MOD(nStart+3,4)+1,iLocSide))%np%NodeID/)
  END DO
END DO

! Build connectivity for the surface mesh

SDEALLOCATE(TempBCSurfNodes)
SDEALLOCATE(TempSideSurfNodeMap)
SDEALLOCATE(SurfMeshPosti%GlobSideToSurfSideMap)
ALLOCATE(TempBCSurfNodes(4*nSides))
ALLOCATE(TempSideSurfNodeMap(1:4,1:nSides))
ALLOCATE(SurfMeshPosti%GlobSideToSurfSideMap(nSides))
SurfMeshPosti%nSurfaceNode=0
SurfMeshPosti%nSurfaceBCSides=0
SurfMeshPosti%GlobSideToSurfSideMap(1:nSides)=0

DO iSide=1, nSides
  IF(BC(iSide).EQ.0) CYCLE
  IF (BoundaryType(BC(iSide),BC_TYPE).EQ.4) THEN                      ! only reflective boundaries
    SurfMeshPosti%nSurfaceBCSides = SurfMeshPosti%nSurfaceBCSides + 1
    SurfMeshPosti%GlobSideToSurfSideMap(iSide) = SurfMeshPosti%nSurfaceBCSides
    iElem = SideToElem(1,iSide)
    IF (iElem.LT.1) THEN
      iElem = SideToElem(2,iSide)
      iLocSide = SideToElem(4,iSide)
    ELSE
      iLocSide = SideToElem(3,iSide)
    END IF
    DO iNode2 = 1, 4
    IsSortedSurfNode = .false.
      DO iNode = 1, SurfMeshPosti%nSurfaceNode 
        IF (GEO%ElemSideNodeID(iNode2, iLocSide, iElem).EQ.TempBCSurfNodes(iNode)) THEN
        TempSideSurfNodeMap(iNode2,SurfMeshPosti%nSurfaceBCSides) = iNode
        IsSortedSurfNode = .true.
        EXIT
        END IF
      END DO
      IF(.NOT.IsSortedSurfNode) THEN
        SurfMeshPosti%nSurfaceNode = SurfMeshPosti%nSurfaceNode + 1
        TempBCSurfNodes(SurfMeshPosti%nSurfaceNode) = GEO%ElemSideNodeID(iNode2, iLocSide, iElem)
        TempSideSurfNodeMap(iNode2,SurfMeshPosti%nSurfaceBCSides) = SurfMeshPosti%nSurfaceNode
      END IF
    END DO  
  END IF
END DO

SDEALLOCATE(SurfMeshPosti%BCSurfNodes)
ALLOCATE(SurfMeshPosti%BCSurfNodes(1:SurfMeshPosti%nSurfaceNode))
SurfMeshPosti%BCSurfNodes(1:SurfMeshPosti%nSurfaceNode) = TempBCSurfNodes(1:SurfMeshPosti%nSurfaceNode)
SDEALLOCATE(SurfMeshPosti%SideSurfNodeMap)
ALLOCATE(SurfMeshPosti%SideSurfNodeMap(1:4,1:SurfMeshPosti%nSurfaceBCSides))
SurfMeshPosti%SideSurfNodeMap(1:4,1:SurfMeshPosti%nSurfaceBCSides) = TempSideSurfNodeMap(1:4,1:SurfMeshPosti%nSurfaceBCSides)
DEALLOCATE(TempBCSurfNodes)
DEALLOCATE(TempSideSurfNodeMap)

END SUBROUTINE InitMesh_Connected