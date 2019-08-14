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
USE MOD_MPI,                 ONLY: InitMPI!,DefineParametersMPI
USE MOD_ReadInTools ,        ONLY: prms,PrintDefaultParameterFile
USE MOD_ReadInTools,         ONLY: GETINT,GETSTR,GETLOGICAL
USE MOD_HDF5_Input,          ONLY: OpenDataFile,CloseDataFile,GetDataProps,ReadAttribute,File_ID,ReadArray,GetDataSize,DatasetExists
USE MOD_HDF5_Input,          ONLY: ISVALIDHDF5FILE,ISVALIDMESHFILE
USE MOD_Mesh_ReadIn,         ONLY: readMesh
USE MOD_Mesh,                ONLY: FinalizeMesh
#ifdef PARTICLES
USE MOD_Particle_Mesh       ,ONLY: FinalizeParticleMesh
#endif
USE MOD_Mesh_Vars,           ONLY: useCurveds,NGeo,nElems,NodeCoords,offsetElem
USE MOD_Interpolation_Vars,  ONLY: NodeTypeCL,NodeTypeVisu
USE MOD_Interpolation,       ONLY: GetVandermonde
USE MOD_ChangeBasis,         ONLY: ChangeBasis3D
USE MOD_VTK,                 ONLY: WriteDataToVTK,WriteVTKMultiBlockDataSet
USE MOD_Prepare_Mesh,        ONLY: fillMeshInfo
#if USE_MPI
USE MOD_MPI_Vars,            ONLY: NbProc,nMPISides_Proc
#endif /*USE_MPI*/
USE MOD_Analyze,             ONLY: CalcErrorStateFiles, CalcErrorStateFileSigma
USE MOD_Analyze_Vars,        ONLY: NAnalyze
USE MOD_Mesh_Vars,           ONLY: sJ,NGeoRef
USE MOD_PreProc,             ONLY: PP_N
USE MOD_Metrics,             ONLY: CalcMetricsErrorDiff
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
CHARACTER(LEN=255)             :: FileString_DG
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
REAL                           :: OutputTime
INTEGER                        :: iDG
CHARACTER(LEN=255)             :: FileString_multiblock
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
INTEGER                        :: iArgsStart
LOGICAL                        :: MeshInitFinished, ReadMeshFinished
! PartData
LOGICAL                        :: VisuParticles, PartDataExists
INTEGER                        :: TimeStampLength
!===================================================================================================================================
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
CALL prms%CreateLogicalOption('VisuParticles',  "Visualize particles (velocity, species, internal energy).", '.FALSE.')
CALL prms%CreateLogicalOption('writePartitionInfo',  "Write information about MPI partitions into a file.",'.FALSE.')
CALL prms%CreateIntOption(    'TimeStampLength', 'Length of the floating number time stamp', '14')
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
VisuParticles    = GETLOGICAL('VisuParticles','.FALSE.')
! Initialization of I/O routines
CALL InitIO()
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
  SWRITE(UNIT_stdOut,'(A,I3,A,I3,A)') 'Processing state ',iArgs-iArgsStart+1,' of ',nArgs-iArgsStart+1,'...'

  ! Open .h5 file
  DGSolutionExists = .FALSE.; ElemDataExists = .FALSE.; SurfaceDataExists = .FALSE.
  CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
  CALL DatasetExists(File_ID,'DG_Solution',DGSolutionExists)
  CALL DatasetExists(File_ID,'ElemData',ElemDataExists)
  CALL DatasetExists(File_ID,'SurfaceData',SurfaceDataExists)
  CALL DatasetExists(File_ID,'PartData',PartDataExists)

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
#if USE_MPI
      SDEALLOCATE(NbProc)
      SDEALLOCATE(nMPISides_Proc)
#endif /*USE_MPI*/
      CALL FinalizeMesh()

      ! Read in parameters from mesh file
      CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
      CALL ReadAttribute(File_ID,'Ngeo',1,IntegerScalar=NGeo)
      CALL CloseDataFile()

      ! Read the mesh itself
      CALL readMesh(MeshFile)
      ReadMeshFinished = .TRUE.

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
      IF(VisuSource)THEN
        CALL ReadArray('DG_Source',5,(/nVar_State,N_State+1_IK,N_State+1_IK,N_State+1_IK,nElems/),offsetElem,5,RealArray=U)
      ELSE
        CALL ReadArray('DG_Solution',5,(/nVar_State,N_State+1_IK,N_State+1_IK,N_State+1_IK,nElems/),offsetElem,5,RealArray=U)
      END IF
    END ASSOCIATE

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
    CALL ConvertElemData(InputStateFile,ReadMeshFinished,MeshInitFinished)
  END IF
  ! === SurfaceData ================================================================================================================
  IF(SurfaceDataExists) THEN
    CALL ConvertSurfaceData(InputStateFile,ReadMeshFinished,MeshInitFinished)
  END IF
  ! === PartData ===================================================================================================================
  IF(VisuParticles) THEN
    IF(PartDataExists) THEN
      CALL ConvertPartData(InputStateFile)
    END IF
  END IF
  CALL CloseDataFile()
END DO ! iArgs = 2, nArgs

! Finalize
SDEALLOCATE(Vdm_N_NVisu)
SDEALLOCATE(Vdm_EQNgeo_NVisu)
SDEALLOCATE(U)
SDEALLOCATE(U_Visu)
SDEALLOCATE(Coords_NVisu)
SDEALLOCATE(NodeCoords)
! visuSurf
CALL FinalizeMesh()
#ifdef PARTICLES
CALL FinalizeParticleMesh()
#endif
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
SWRITE(UNIT_stdOut,'(A,F14.2,A)') ' H5PICLAS2VTK FINISHED! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')

END PROGRAM H5PICLAS2VTK

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
INTEGER,INTENT(IN)            :: ConnectInfo(data_size,nNodes)      ! Statevector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iVal,iElem,Offset,nBytes,nVTKElems,nVTKCells,ivtk=44,iVar,iNode, int_size, ElemType
INTEGER                       :: Vertex(data_size,nElems), iLen, str_len
CHARACTER(LEN=35)             :: StrOffset,TempStr1,TempStr2
CHARACTER(LEN=200)            :: Buffer, tmp, tmp2, VarNameString
CHARACTER(LEN=1)              :: lf, components_string
REAL(KIND=4)                  :: float
INTEGER,ALLOCATABLE           :: VarNameCombine(:), VarNameCombineLen(:)
!===================================================================================================================================

nVTKElems=nNodes
nVTKCells=nElems

SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')"   WRITE 3D DATA TO VTX XML BINARY (VTU) FILE..."
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
      ! If it is the first occurence, start counting
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
Offset=Offset+INT(SIZEOF(int_size),4)+3*nVTKElems*INT(SIZEOF(float),4)
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
nBytes = 3*nVTKElems*INT(SIZEOF(FLOAT),4)
WRITE(ivtk) nBytes
WRITE(ivtk) REAL(Coords(1:3,1:nVTKElems),4)
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


SUBROUTINE InitMesh_Connected()
!===================================================================================================================================
!> Reusing the routines from Prepare_Mesh and Particle_Mesh to build-up the mesh and its connectivity (for conforming meshes)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
USE MOD_PreProc
USE MOD_Prepare_Mesh            ,ONLY: setLocalSideIDs,fillMeshInfo
#if USE_MPI
USE MOD_Prepare_Mesh            ,ONLY: exchangeFlip
#endif
#ifdef PARTICLES
USE MOD_Particle_Mesh           ,ONLY: InitParticleGeometry
#endif
!--------------------------------------------------------------------------------------------------!
! perform Mapping for Surface Output
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE
! LOCAL VARIABLES
!===================================================================================================================================
PP_nElems=nElems
SWRITE(UNIT_stdOut,'(A)') "NOW CALLING setLocalSideIDs..."
CALL setLocalSideIDs()

#if USE_MPI
SWRITE(UNIT_stdOut,'(A)') "NOW CALLING exchangeFlip..."
CALL exchangeFlip()
#endif

firstBCSide          = 1
firstMortarInnerSide = firstBCSide         +nBCSides
firstInnerSide       = firstMortarInnerSide+nMortarInnerSides
firstMPISide_MINE    = firstInnerSide      +nInnerSides
firstMPISide_YOUR    = firstMPISide_MINE   +nMPISides_MINE
firstMortarMPISide   = firstMPISide_YOUR   +nMPISides_YOUR

lastBCSide           = firstMortarInnerSide-1
lastMortarInnerSide  = firstInnerSide    -1
lastInnerSide        = firstMPISide_MINE -1
lastMPISide_MINE     = firstMPISide_YOUR -1
lastMPISide_YOUR     = firstMortarMPISide-1
lastMortarMPISide    = nSides

SDEALLOCATE(ElemToSide)
SDEALLOCATE(SideToElem)
SDEALLOCATE(BC)
SDEALLOCATE(AnalyzeSide)
ALLOCATE(ElemToSide(2,6,nElems))
ALLOCATE(SideToElem(5,nSides))
ALLOCATE(BC(1:nSides))
ALLOCATE(AnalyzeSide(1:nSides))
ElemToSide  = 0
SideToElem  = -1   !mapping side to elem, sorted by side ID (for surfint)
BC          = 0
AnalyzeSide = 0

SDEALLOCATE(MortarType)
SDEALLOCATE(MortarInfo)
SDEALLOCATE(MortarSlave2MasterInfo)
ALLOCATE(MortarType(2,1:nSides))              ! 1: Type, 2: Index in MortarInfo
ALLOCATE(MortarInfo(MI_FLIP,4,nMortarSides)) ! [1]: 1: Neighbour sides, 2: Flip, [2]: small sides
ALLOCATE(MortarSlave2MasterInfo(1:nSides))
MortarType=-1
MortarInfo=-1

SWRITE(UNIT_stdOut,'(A)') "NOW CALLING fillMeshInfo..."
CALL fillMeshInfo()

#ifdef PARTICLES
CALL InitParticleGeometry()
#else
CALL abort(__STAMP__,&
      'ERROR: Post-processing tool h5piclas2vtk was compiled with PARTICLES=OFF! No DSMC/ElemData output supported!')
#endif

END SUBROUTINE InitMesh_Connected


SUBROUTINE ConvertPartData(InputStateFile)
!===================================================================================================================================
! Subroutine to write 3D point data to VTK format
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,           ONLY: ProjectName
USE MOD_IO_HDF5,                ONLY: HSize
USE MOD_HDF5_Input,             ONLY: OpenDataFile,ReadAttribute,File_ID,ReadArray,GetDataSize
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
REAL, ALLOCATABLE               :: PartData(:,:), tmpPartData(:,:)
REAL                            :: OutputTime
!===================================================================================================================================

CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)
CALL ReadAttribute(File_ID,'Time',1,RealScalar=OutputTime)
! Read-in of dimensions of the particle array (1: Number of particles, 2: Number of variables)
CALL GetDataSize(File_ID,'PartData',nDims,HSize)
! First 3 entries are the particle positions, which are used as the coordinates for the output and not included as a variable
nPartsVar=INT(HSize(2),4)-3
nParts=INT(HSize(1),4)
! Allocating the array for the variables and a temporary array since ParticlePositionX,Y,Z are included in the read-in
ALLOCATE(VarNamesParticle(nPartsVar),tmpArray(nPartsVar+3))
CALL ReadAttribute(File_ID,'VarNamesParticles',nPartsVar+3,StrArray=tmpArray)
VarNamesParticle(1:nPartsVar)=tmpArray(4:nPartsVar+3)

IF(nParts.GT.0) THEN
  ALLOCATE(PartData(1:nPartsVar+3,1:nParts),tmpPartData(1:nParts,1:nPartsVar+3))
  PartData = 0.
  tmpPartData = 0.
  SDEALLOCATE(ConnectInfo)
  ALLOCATE(ConnectInfo(1,1:nParts))
  ConnectInfo = 0
END IF

ASSOCIATE(nParts    => INT(nParts,IK),  &
          nPartsVar => INT(nPartsVar,IK))
CALL ReadArray('PartData',2,(/nParts, nPartsVar+3_IK/),0_IK,1,RealArray=tmpPartData)
END ASSOCIATE

DO iPart=1,nParts
  PartData(1:nPartsVar+3,iPart) = tmpPartData(iPart,1:nPartsVar+3)
  ConnectInfo(1,iPart)=iPart-1
END DO

FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_visuPart',OutputTime))//'.vtu'
CALL WriteDataToVTK_PICLas(1,FileString,nPartsVar,VarNamesParticle,nParts,PartData(1:3,1:nParts),nParts,&
                            PartData(4:nPartsVar+3,1:nParts),ConnectInfo(1,1:nParts))

SDEALLOCATE(VarNamesParticle)
SDEALLOCATE(tmpArray)
SDEALLOCATE(PartData)
SDEALLOCATE(ConnectInfo)
SDEALLOCATE(tmpPartData)

END SUBROUTINE ConvertPartData


SUBROUTINE ConvertElemData(InputStateFile,ReadMeshFinished,MeshInitFinished)
!===================================================================================================================================
! Subroutine to write 3D point data to VTK format
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,           ONLY: ProjectName
USE MOD_IO_HDF5,                ONLY: HSize
USE MOD_HDF5_Input,             ONLY: OpenDataFile,CloseDataFile,ReadAttribute,File_ID,ReadArray,GetDataSize
USE MOD_Mesh_ReadIn,            ONLY: readMesh
USE MOD_Mesh_Vars,              ONLY: NGeo, nElems, nNodes, offsetElem
USE MOD_Particle_Mesh_Vars,     ONLY: GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)   :: InputStateFile
LOGICAL,INTENT(INOUT)           :: ReadMeshFinished,MeshInitFinished
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)              :: FileString, MeshFile, File_Type
REAL                            :: OutputTime
INTEGER                         :: nDims,nVarAdd
CHARACTER(LEN=255),ALLOCATABLE  :: VarNamesAdd(:)
REAL,ALLOCATABLE                :: ElemData(:,:)
!===================================================================================================================================

CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
! DSMC/FV Solution
CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)
CALL ReadAttribute(File_ID,'Time',1,RealScalar=OutputTime)
CALL ReadAttribute(File_ID,'File_Type',1,StrScalar=File_Type)

IF(.NOT.ReadMeshFinished) THEN
  CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
  CALL ReadAttribute(File_ID,'Ngeo',1,IntegerScalar=NGeo)
  CALL CloseDataFile()
  CALL readMesh(MeshFile)
  ReadMeshFinished = .TRUE.
END IF
IF(.NOT.MeshInitFinished) THEN
  CALL InitMesh_Connected()
  MeshInitFinished = .TRUE.
END IF

! Read in solution
CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)

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
    CASE('DSMCHOState')
      FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_visuDSMC',OutputTime))//'.vtu'
  END SELECT
  CALL WriteDataToVTK_PICLas(8,FileString,nVarAdd,VarNamesAdd(1:nVarAdd),nNodes,GEO%NodeCoords(1:3,1:nNodes),nElems,&
                              ElemData(1:nVarAdd,1:nElems),GEO%ElemToNodeID(1:8,1:nElems))
END IF

SDEALLOCATE(VarNamesAdd)
SDEALLOCATE(ElemData)

END SUBROUTINE ConvertElemData


SUBROUTINE ConvertSurfaceData(InputStateFile,ReadMeshFinished,MeshInitFinished)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,           ONLY: ProjectName
USE MOD_IO_HDF5,                ONLY: HSize
USE MOD_HDF5_Input,             ONLY: OpenDataFile,CloseDataFile,ReadAttribute,GetDataSize,File_ID,ReadArray,GetDataSize
USE MOD_Mesh_ReadIn,            ONLY: readMesh
USE MOD_Mesh_Vars,              ONLY: NGeo, SurfConnect
USE MOD_Particle_Mesh_Vars,     ONLY: GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)   :: InputStateFile
LOGICAL,INTENT(INOUT)           :: ReadMeshFinished,MeshInitFinished
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)              :: FileString, MeshFile
REAL                            :: OutputTime
INTEGER                         :: nDims, nVarSurf, nSurfSample, iNode
CHARACTER(LEN=255),ALLOCATABLE  :: VarNamesSurf_HDF5(:)
REAL, ALLOCATABLE               :: tempSurfData(:,:,:,:), SurfData(:,:), Coords(:,:)
!===================================================================================================================================

CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)
CALL ReadAttribute(File_ID,'Time',1,RealScalar=OutputTime)
CALL ReadAttribute(File_ID,'DSMC_nSurfSample',1,IntegerScalar=nSurfSample)
IF(nSurfSample.EQ.1) THEN
  ! Read-in of mesh information (if not already done for DG solution -> DSMCState case)
  IF(.NOT.ReadMeshFinished) THEN
    ! Read in parameters from mesh file
    CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
    CALL ReadAttribute(File_ID,'Ngeo',1,IntegerScalar=NGeo)
    CALL CloseDataFile()
    CALL readMesh(MeshFile)
    ReadMeshFinished = .TRUE.
  END IF
  IF(.NOT.MeshInitFinished) THEN
    CALL InitMesh_Connected()
    MeshInitFinished = .TRUE.
  END IF
  ! Read in solution
  CALL OpenDataFile(InputStateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)

  CALL BuildSurfMeshConnectivity()

  CALL GetDataSize(File_ID,'SurfaceData',nDims,HSize)
  nVarSurf = INT(HSize(1),4)
  ALLOCATE(VarNamesSurf_HDF5(nVarSurf))
  CALL ReadAttribute(File_ID,'VarNamesSurface',nVarSurf,StrArray=VarNamesSurf_HDF5(1:nVarSurf))

  IF((nVarSurf.GT.0).AND.(SurfConnect%nSurfaceBCSides.GT.0))THEN
    ALLOCATE(SurfData(1:nVarSurf,1:SurfConnect%nSurfaceBCSides))
    ALLOCATE(tempSurfData(1:nVarSurf,nSurfSample,nSurfSample,1:SurfConnect%nSurfaceBCSides))
    SurfData=0.
    tempSurfData = 0.
    ASSOCIATE(nVarSurf        => INT(nVarSurf,IK),  &
              nSurfaceBCSides => INT(SurfConnect%nSurfaceBCSides,IK))
      CALL ReadArray('SurfaceData',4,(/nVarSurf, 1_IK, 1_IK, nSurfaceBCSides/), &
                      0_IK,4,RealArray=tempSurfData(:,:,:,:))
    END ASSOCIATE
    SurfData(1:nVarSurf,1:SurfConnect%nSurfaceBCSides) = tempSurfData(1:nVarSurf,1,1,1:SurfConnect%nSurfaceBCSides)
    ALLOCATE(Coords(1:3,SurfConnect%nSurfaceNode))
    Coords = 0.
  END IF

  DO iNode=1, SurfConnect%nSurfaceNode
    Coords(1:3,iNode) = GEO%NodeCoords(1:3, SurfConnect%BCSurfNodes(iNode))
  END DO

  FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_visuSurf',OutputTime))//'.vtu'
  CALL WriteDataToVTK_PICLas(4,FileString,nVarSurf,VarNamesSurf_HDF5,SurfConnect%nSurfaceNode,Coords(1:3,:),&
      SurfConnect%nSurfaceBCSides,SurfData,SurfConnect%SideSurfNodeMap(1:4,1:SurfConnect%nSurfaceBCSides))
ELSE IF(nSurfSample.GT.1) THEN
  CALL abort(__STAMP__,&
    'DSMC_nSurfSample greater one is not supported yet!')
ELSE
  CALL abort(__STAMP__,&
    'DSMC_nSurfSample is zero, something might be wrong with your .h5 file!')
END IF

SDEALLOCATE(VarNamesSurf_HDF5)
SDEALLOCATE(SurfData)
SDEALLOCATE(tempSurfData)
SDEALLOCATE(Coords)

END SUBROUTINE ConvertSurfaceData


SUBROUTINE BuildSurfMeshConnectivity()
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_IO_HDF5,                ONLY: HSize
USE MOD_HDF5_Input,             ONLY: OpenDataFile,CloseDataFile,ReadAttribute,GetDataSize,File_ID,ReadArray,GetDataSize
USE MOD_Mesh_ReadIn,            ONLY: readMesh
USE MOD_Mesh_Vars,              ONLY: SurfConnect, nSides, SideToElem, BC, BoundaryName
USE MOD_Particle_Mesh_Vars,     ONLY: GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: nDims, iNode, nSurfBC_HDF5, iName
CHARACTER(LEN=255),ALLOCATABLE  :: SurfBCName_HDF5(:)
INTEGER                         :: iElem, iLocSide, iSide, iNode2, iBC, nSurfSides
INTEGER, ALLOCATABLE            :: TempBCSurfNodes(:), TempSideSurfNodeMap(:,:), SideToSurfSide(:)
LOGICAL                         :: IsSortedSurfNode
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A,A)')' GET NUMBER AND NAMES OF SURFACE-BCSIDES IN HDF5 FILE... '
CALL GetDataSize(File_ID,'BC_Surf',nDims,HSize,attrib=.true.)
nSurfBC_HDF5 = INT(HSize(1),4)

SWRITE(UNIT_stdOut,'(A3,A45,A3,I33,A13)')' | ','Number of Surface BCs',' | ',nSurfBC_HDF5,' | HDF5    | '
ALLOCATE(SurfBCName_HDF5(nSurfBC_HDF5))
CALL ReadAttribute(File_ID,'BC_Surf',nSurfBC_HDF5,StrArray=SurfBCName_HDF5)
DO iName = 1,nSurfBC_HDF5
  SWRITE(UNIT_stdOut,'(A3,A38,I2.1,A5,A3,A33,A13)')' | ','BC',iName,'Name',' | ',TRIM(SurfBCName_HDF5(iName)),' | HDF5    | '
END DO

! create sideid to surfaceid map for all surface-sides contained in statefile
ALLOCATE(SideToSurfSide(1:nSides))
SideToSurfSide(1:nSides) = -1
! if side is surface (BC_reflective defined in State file) then map respective surface side number
nSurfSides = 0
DO iSide = 1,nSides
  IF(BC(iSide).EQ.0) CYCLE
  DO iBC=1,nSurfBC_HDF5
    IF((TRIM(BoundaryName(BC(iSide))) .EQ. TRIM(SurfBCName_HDF5(iBC)))) THEN
      nSurfSides = nSurfSides + 1
      SideToSurfSide(iSide) = nSurfSides
    END IF
  END DO
END DO

! Build connectivity for the surface mesh
ALLOCATE(TempBCSurfNodes(4*nSides))
ALLOCATE(TempSideSurfNodeMap(1:4,1:nSides))
SurfConnect%nSurfaceNode=0
SurfConnect%nSurfaceBCSides=0

DO iSide=1, nSides
  IF(BC(iSide).EQ.0) CYCLE
  IF (SideToSurfSide(iSide).NE.-1) THEN
    SurfConnect%nSurfaceBCSides = SurfConnect%nSurfaceBCSides + 1
    iElem = SideToElem(1,iSide)
    IF (iElem.LT.1) THEN
      iElem = SideToElem(2,iSide)
      iLocSide = SideToElem(4,iSide)
    ELSE
      iLocSide = SideToElem(3,iSide)
    END IF
    DO iNode2 = 1, 4
    IsSortedSurfNode = .false.
      DO iNode = 1, SurfConnect%nSurfaceNode
        IF (GEO%ElemSideNodeID(iNode2, iLocSide, iElem).EQ.TempBCSurfNodes(iNode)) THEN
          TempSideSurfNodeMap(iNode2,SurfConnect%nSurfaceBCSides) = iNode
          IsSortedSurfNode = .true.
          EXIT
        END IF
      END DO
      IF(.NOT.IsSortedSurfNode) THEN
        SurfConnect%nSurfaceNode = SurfConnect%nSurfaceNode + 1
        TempBCSurfNodes(SurfConnect%nSurfaceNode) = GEO%ElemSideNodeID(iNode2, iLocSide, iElem)
        TempSideSurfNodeMap(iNode2,SurfConnect%nSurfaceBCSides) = SurfConnect%nSurfaceNode
      END IF
    END DO
  END IF
END DO

SDEALLOCATE(SurfConnect%BCSurfNodes)
ALLOCATE(SurfConnect%BCSurfNodes(1:SurfConnect%nSurfaceNode))
SurfConnect%BCSurfNodes(1:SurfConnect%nSurfaceNode) = TempBCSurfNodes(1:SurfConnect%nSurfaceNode)
SDEALLOCATE(SurfConnect%SideSurfNodeMap)
ALLOCATE(SurfConnect%SideSurfNodeMap(1:4,1:SurfConnect%nSurfaceBCSides))
SurfConnect%SideSurfNodeMap(1:4,1:SurfConnect%nSurfaceBCSides) = TempSideSurfNodeMap(1:4,1:SurfConnect%nSurfaceBCSides)
SDEALLOCATE(TempBCSurfNodes)
SDEALLOCATE(TempSideSurfNodeMap)
SDEALLOCATE(SurfBCName_HDF5)
SDEALLOCATE(SideToSurfSide)

END SUBROUTINE BuildSurfMeshConnectivity
