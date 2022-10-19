!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "piclas.h"

!===================================================================================================================================
!> Module containing the main procedures for the visu tool: visu_requestInformation is called by ParaView to create a
!> list of available variables and visu is the main routine which is either called by ParaView to get the data it visualizes
!> or by the standalone tool.
!===================================================================================================================================
MODULE MOD_Visu
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE visu_getVarNamesAndFileType
  MODULE PROCEDURE visu_getVarNamesAndFileType
END INTERFACE

INTERFACE Visu_InitFile
  MODULE PROCEDURE Visu_InitFile
END INTERFACE

INTERFACE visu
  MODULE PROCEDURE visu
END INTERFACE

INTERFACE FinalizeVisu
  MODULE PROCEDURE FinalizeVisu
END INTERFACE

PUBLIC:: visu_getVarNamesAndFileType
PUBLIC:: visu_InitFile
PUBLIC:: visu
PUBLIC:: FinalizeVisu

CONTAINS

!===================================================================================================================================
!> Create a list of available variables for ParaView. This list contains the conservative, primitve and derived quantities
!> that are available in the current equation system as well as the additional variables read from the state file.
!> The additional variables are stored in the datasets 'ElemData' (elementwise data) and 'FieldData' (pointwise data).
!> Also a list of all available boundary names is created for surface visualization.
!===================================================================================================================================
SUBROUTINE visu_getVarNamesAndFileType(mpi_comm_IN,statefile,meshfile,varnames_loc, bcnames_loc)
USE MOD_Globals
USE MOD_Visu_Vars      ,ONLY: FileType,VarNamesHDF5,nBCNamesAll
USE MOD_HDF5_Input     ,ONLY: OpenDataFile,CloseDataFile,GetDataSize,GetVarNames,ISVALIDMESHFILE,ISVALIDHDF5FILE,ReadAttribute
USE MOD_HDF5_Input     ,ONLY: DatasetExists,HSize,nDims,ReadArray
USE MOD_IO_HDF5        ,ONLY: GetDatasetNamesInGroup,File_ID
USE MOD_StringTools    ,ONLY: STRICMP
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                                  :: mpi_comm_IN
CHARACTER(LEN=255),INTENT(IN)                       :: statefile
CHARACTER(LEN=*)  ,INTENT(IN)                       :: meshfile
CHARACTER(LEN=255),INTENT(INOUT),ALLOCATABLE,TARGET :: varnames_loc(:)
CHARACTER(LEN=255),INTENT(INOUT),ALLOCATABLE,TARGET :: bcnames_loc(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                             :: i,j,nVar,dims
LOGICAL                                             :: varnames_found,readDGsolutionVars,sameVars,VarNamesExist, file_exists
CHARACTER(LEN=255),ALLOCATABLE                      :: datasetNames(:)
CHARACTER(LEN=255),ALLOCATABLE                      :: varnames_tmp(:)
CHARACTER(LEN=255),ALLOCATABLE                      :: tmp(:)
CHARACTER(LEN=255)                                  :: MeshFile_loc
INTEGER                                             :: Offset=0 ! Every process reads all BCs
!===================================================================================================================================
sameVars=.FALSE.

IF (ISVALIDMESHFILE(statefile)) THEN      ! MESH
  SDEALLOCATE(varnames_loc)
  FileType='Mesh'
ELSE IF (ISVALIDHDF5FILE(statefile)) THEN ! other file
  SDEALLOCATE(varnames_loc)
  CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=mpi_comm_IN)
  CALL ReadAttribute(File_ID,'File_Type',   1,StrScalar =FileType)

  ! check if variables in state file are the same as in the EQNSYS
  ! and set FileType to 'Generic' if not
  IF (STRICMP(FileType,'State')) THEN
    SDEALLOCATE(VarNamesHDF5)
    CALL GetVarNames("VarNames",VarNamesHDF5,VarNamesExist)
    IF (.NOT.sameVars) FileType='Generic'
  END IF

  !IF (STRICMP(FileType,'State')) THEN
  !  nVar = nVarDepEOS
  !  ALLOCATE(varnames_loc(nVar))
  !  varnames_loc(1:nVar)=DepNames
  !  readDGsolutionVars = .FALSE.
  !ELSE
    nVar=0
    readDGsolutionVars = .TRUE.
  !END IF

  CALL GetDatasetNamesInGroup("/",datasetNames)

  DO i=1,SIZE(datasetNames)
    SDEALLOCATE(varnames_tmp)
    VarNamesExist=.FALSE.
    CALL DatasetExists(File_ID,"VarNames_"//TRIM(datasetNames(i)),varnames_found,attrib=.TRUE.)
    IF (varnames_found) THEN
          CALL GetVarNames("VarNames_"//TRIM(datasetNames(i)),varnames_tmp,VarNamesExist)
    ELSE
      IF (STRICMP(datasetNames(i), "DG_Solution")) THEN
        IF (readDGsolutionVars) THEN
          CALL GetVarNames("VarNames",varnames_tmp,VarNamesExist)
        END IF
      ELSE IF(STRICMP(datasetNames(i), "ElemData")) THEN
        CALL GetVarNames("VarNamesAdd",varnames_tmp,VarNamesExist)
      ELSE IF(STRICMP(datasetNames(i), "PartData")) THEN
        CYCLE ! Do not add particle data
      ELSE IF(STRICMP(datasetNames(i), "PartInt")) THEN
        CYCLE ! Do not add particle location in cells info
      ELSE
        CALL GetDataSize(File_ID,TRIM(datasetNames(i)),dims,HSize)
        IF ((dims.NE.5).AND.(dims.NE.2)) CYCLE ! Do not add datasets to the list that can not contain elementwise or field data
        ALLOCATE(varnames_tmp(INT(HSize(1))))
        DO j=1,INT(HSize(1))
          WRITE(varnames_tmp(j),'(I0)') j
        END DO
        DEALLOCATE(HSize)
        VarNamesExist=.TRUE.
      END IF
    END IF
    IF (.NOT.VarNamesExist) CYCLE

    ! increase array 'varnames_loc'
    IF (nVar.GT.0) THEN
      ALLOCATE(tmp(nVar))
      tmp = varnames_loc
    END IF
    SDEALLOCATE(varnames_loc)
    ALLOCATE(varnames_loc(nVar+SIZE(varnames_tmp)))
    IF (nVar.GT.0) varnames_loc(1:nVar) = tmp(1:nVar)
    SDEALLOCATE(tmp)

    ! copy new varnames from varnames_tmp to varnames_loc
    DO j=1,SIZE(varnames_tmp)
      varnames_loc(nVar+j) = TRIM(datasetNames(i))//":"//TRIM(varnames_tmp(j))
    END DO
    nVar = nVar + SIZE(varnames_tmp)
  END DO

  IF (LEN_TRIM(meshfile).EQ.0) THEN
    ! Save mesh file to get boundary names later
    CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar =MeshFile_loc)
  ELSE
    MeshFile_loc = meshfile
  END IF

  CALL CloseDataFile()

  INQUIRE(FILE=TRIM(MeshFile_loc), EXIST=file_exists)
  IF (file_exists) THEN
    ! Open the mesh file and read all boundary names for surface visualization
    CALL OpenDataFile(MeshFile_loc,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=mpi_comm_IN)
    CALL GetDataSize(File_ID,'BCNames',nDims,HSize)
    CHECKSAFEINT(HSize(1),4)
    nBCNamesAll=INT(HSize(1),4)
    DEALLOCATE(HSize)
    SDEALLOCATE(bcnames_loc)
    ALLOCATE(bcnames_loc(nBCNamesAll))
    ASSOCIATE (&
          Offset      => INT(Offset,IK)      ,&
          nBCNamesAll => INT(nBCNamesAll,IK)  &
          )
      CALL ReadArray('BCNames',1,(/nBCNamesAll/),Offset,1,StrArray=bcnames_loc)
    END ASSOCIATE
    CALL CloseDataFile()
  END IF

  SDEALLOCATE(datasetNames)
END IF
END SUBROUTINE visu_getVarNamesAndFileType

!===================================================================================================================================
!> This routine is used to prepare everything we need to visualize data from a statefile.
!> This includes:
!> * Get the mesh file
!> * Read the desired visualization polynomial degree, the visualization dimennsion, the node type we want to visualize on and the
!>   Dg only option
!> * Decide whether the state file, the mesh file, the visualization polynomial degree or the dg only option changed. This is
!>   needed to decide what parts of the visualization routines should be called.
!> * Call routines that build the distribution between DG elements and the mappings needed to calculate and visualize the
!>   desired variables.
!===================================================================================================================================
SUBROUTINE visu_InitFile(mpi_comm_IN,statefile,postifile)
! MODULES
USE HDF5
USE MOD_Preproc
USE MOD_Globals
USE MOD_Visu_Vars
USE MOD_MPI                ,ONLY: InitMPI
USE MOD_HDF5_Input         ,ONLY: ISVALIDMESHFILE,ISVALIDHDF5FILE,GetArrayAndName
USE MOD_HDF5_Input         ,ONLY: ReadAttribute,File_ID,OpenDataFile,GetDataProps,CloseDataFile,ReadArray,DatasetExists
USE MOD_Interpolation_Vars ,ONLY: NodeType
USE MOD_Globals_Vars       ,ONLY: ProjectName
USE MOD_StringTools        ,ONLY: STRICMP
USE MOD_ReadInTools        ,ONLY: prms,GETINT,GETLOGICAL,addStrListEntry,GETSTR
USE MOD_Posti_Mappings     ,ONLY: Build_mapDepToCalc_mapAllVarsToVisuVars
IMPLICIT NONE
INTEGER,INTENT(IN)               :: mpi_comm_IN
CHARACTER(LEN=255),INTENT(IN)    :: statefile
CHARACTER(LEN=255),INTENT(INOUT) :: postifile
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL
INTEGER                          :: nElems_State,N_HDF5
CHARACTER(LEN=255)               :: NodeType_State, cwd
!===================================================================================================================================
IF (STRICMP(fileType,'Mesh')) THEN
    CALL CollectiveStop(__STAMP__, &
        "FileType==Mesh, but we try to initialize a state file!")
END IF

! open state file to be able to read attributes
CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=mpi_comm_IN)

! read the meshfile attribute from statefile
CALL ReadAttribute(File_ID,'MeshFile',    1,StrScalar =MeshFile_state)

! read options from posti parameter file
CALL prms%read_options(postifile)
NVisu             = GETINT("NVisu")
! again read MeshFile from posti prm file (this overwrites the MeshFile read from the state file)

Meshfile          =  GETSTR("MeshFile",MeshFile_state)
IF (.NOT.FILEEXISTS(MeshFile)) THEN
  !!!!!!
  ! WARNING: GETCWD is a GNU extension to the Fortran standard and will probably not work on other compilers
  CALL GETCWD(cwd)
  !!!!!!
  Meshfile          =  TRIM(cwd) // "/" // TRIM(Meshfile)
END IF
NodeTypeVisuPosti = GETSTR('NodeTypeVisu')
CALL CloseDataFile()

CALL visu_getVarNamesAndFileType(mpi_comm_IN,statefile,"",VarnamesAll,BCNamesAll)

CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=mpi_comm_IN)

! check if state, mesh or NVisu changed
changedStateFile = .NOT.STRICMP(statefile,statefile_old)
changedMeshFile  = .NOT.(STRICMP(MeshFile,MeshFile_old))

SWRITE(*,*) "state file old -> new: ", TRIM(statefile_old), " -> ",TRIM(statefile)
SWRITE(*,*) " mesh file old -> new: ", TRIM(MeshFile_old) , " -> ",TRIM(MeshFile)

! if Mesh or State changed readin some more attributes/parameters
IF (changedStateFile.OR.changedMeshFile) THEN
  CALL GetDataProps('DG_Solution',nVar_State,N_HDF5,nElems_State,NodeType_State)
  IF (.NOT.STRICMP(NodeType_State, NodeType)) THEN
    CALL CollectiveStop(__STAMP__, &
        "NodeType of state does not match with NodeType the visu-posti is compiled with!")
  END IF
  CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar =ProjectName)
  CALL ReadAttribute(File_ID,'Time',        1,RealScalar=OutputTime)
END IF

CALL CloseDataFile()

! Check for changed visualization basis here to take change done for average output into account
IF(NVisu.NE.NVisu_old)THEN
  SWRITE(*,*) "Node type Visu has changed"
END IF ! NVisu.NE.NVisu_old
IF(NodeTypeVisuPosti.NE.NodeTypeVisuPosti_old)THEN
  SWRITE(*,*) "Node type VisuPosti has changed"
END IF ! NodeTypeVisuPosti.NE.NodeTypeVisuPosti_old
changedNVisu     = ((NVisu.NE.NVisu_old) .OR. (NodeTypeVisuPosti.NE.NodeTypeVisuPosti_old))

! set number of dependent and raw variables
!SDEALLOCATE(DepTable)
!SDEALLOCATE(DepSurfaceOnly)
!SDEALLOCATE(DepVolumeOnly)
nVarAll=SIZE(VarnamesAll)
IF (STRICMP(FileType,'State')) THEN
  StateFileMode = .TRUE.
  !nVarDep = nVarDepEOS
  !ALLOCATE(DepTable(nVarDep,0:nVarDep))
  !ALLOCATE(DepSurfaceOnly(nVarDep))
  !ALLOCATE(DepVolumeOnly(nVarDep))
  !DepTable = DepTableEOS
  !DepSurfaceOnly = DepSurfaceOnlyEOS
  !DepVolumeOnly  = DepVolumeOnlyEOS
ELSE
  StateFileMode = .FALSE.
  !nVarDep = 0
  !ALLOCATE(DepTable(nVarDep,0:nVarDep))
  !ALLOCATE(DepSurfaceOnly(nVarDep))
  !ALLOCATE(DepVolumeOnly(nVarDep))
  !DepTable = 0
  !DepSurfaceOnly = 0
  !DepVolumeOnly  = 0
END IF

! build mappings of variables which must be calculated/visualized
! also set withDGOperator flag if a dependent variable requires the evaluation of the DG operator
CALL Build_mapDepToCalc_mapAllVarsToVisuVars()

END SUBROUTINE visu_InitFile

!===================================================================================================================================
!> Main routine of the visualization tool visu. Called either by the ParaView plugin or by the standalone program version.
!===================================================================================================================================
SUBROUTINE visu(mpi_comm_IN, prmfile, postifile, statefile)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_MPI                 ,ONLY: InitMPI
USE MOD_HDF5_Input          ,ONLY: ISVALIDMESHFILE,ISVALIDHDF5FILE,CloseDataFile
USE MOD_Posti_ReadState     ,ONLY: ReadState
USE MOD_Posti_VisuMesh      ,ONLY: VisualizeMesh
!USE MOD_Posti_Calc          ,ONLY: CalcQuantities_DG
USE MOD_Posti_ConvertToVisu ,ONLY: ConvertToVisu_DG,ConvertToVisu_GenericData
USE MOD_ReadInTools         ,ONLY: prms,ExtractParameterFile,FinalizeParameters
USE MOD_StringTools         ,ONLY: STRICMP
USE MOD_Posti_VisuMesh      ,ONLY: BuildVisuCoords
USE MOD_Posti_Mappings      ,ONLY: Build_mapBCSides
USE MOD_IO_HDF5             ,ONLY: InitMPIInfo
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)               :: mpi_comm_IN
CHARACTER(LEN=255),INTENT(INOUT) :: prmfile
CHARACTER(LEN=255),INTENT(INOUT) :: postifile
CHARACTER(LEN=255),INTENT(IN)    :: statefile
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                          :: changedPrmFile
!===================================================================================================================================

!**********************************************************************************************
! General workflow / principles of the visu ParaView-plugin
!
! * all arrays are SDEALLOCATEd just before they are allocated. This is done to keep there
!   content during successive calls of the visu during a ParaView session. They are only
!   deallocated and reallocated, if there content should change. For example the coords of the
!   mesh file only change if the mesh, NVisu or the distribution of DG/FV elements changes.
!
! VISUALIZE MESH: Call 'VisualizeMesh' routine, which reads the mesh, interpolates it to
!   the visu grid and writes it to VTK 3D arrays.
!
! VISUALIZE STATE:
! * There are two different modes:
!   - without gradients: All quantity that should be visualized can be computed without any
!                        gradient computation (mainly the conservative and primitive quantities)
!                        If FV_RECONSTRUCT is enabled and there are FV elements present in
!                        the state, then this mode is not available.
!                        U is read from the state file directly and the EOS is initialized to
!                        perform ConsToPrim conversions based only on the conservative state U.
!
!   - with gradiens: There are quantities that require the computation of gradients or there
!                    are FV elements with FV_RECONSTRUCT enabled. In this case the DG operator
!                    'DGTimeDerivative_weakForm' is called once to fill the gradients and the
!                    reconstruction of the FV subcell method.
!                    This requires the initialization of several modules of the FLEXI.
!                    U is read via a call of 'Restart'. In the DGTimeDerivative_weakForm the
!                    primitive quantities U_Prim and gradUx/y/z as well as gradUxi/eta/zeta are
!                    filled. These are used to calculate the visu-quantities.
!
! * The whole calculation of derived quantities is performed on PP_N and afterwards
!   interpolated to NVisu. This is not the case for FV elements with FV_RECONSTRUCT enabled.
!   These require to reconstruct the solution first to the visu grid and afterwards can
!   calculate the derived quantities on the NVisu_FV grid.
!
! * The dependencies of the visu-quantities on the state-quantities is stored in a dependency
!   integer table 'DepTable' and it corresponding row/col names 'DepNames' (see eos_vars.f90).
!   The string vector 'DepNames' contains all available quantities available for visualization
!   and the 'DepTable' contains in a row the dependencies of a quantity on other quantities.
!
! * The calculation of the visu-quantities is done in two steps:
!   1. calculate all quantities needed for the visu-quantities (stored in UCalc)
!   2. pick the visu-quantities from UCalc and interpolate them to NVisu (stored in UVisu)
!
! * Therefore two mappings from all available quantities to the calc-quantities and the
!   visu-quantities exist:
!   - 'mapAllVarsToVisuVars' is a integer array of size (1:nVarAll), where nVarAll is the total amount
!     of available quantities. This map contains a zero for all not-to-visu-quantities and for
!     all quantities the index where it is stored in 'UVisu'.
!     This mapping is filled from the 'VarName' entries in the parameter file.
!   - 'mapDepToCalc' is the same as mapAllVarsToVisuVars, but for all (intermediate) quantities stored in 'UCalc'.
!     This mapping is filled from the DepTable.
!
! CHANGED system:
! * There are different logical changedXXX variables, which indicated if XXX changed during
!   successive calls of the visu. These variables control the general workflow of the visu.
!   - changedStateFile:     new state file
!   - changedMeshFile:      new mesh file (only possible if changedStateFile==TRUE)
!   - changedVarNames:      new set of variables to visualize
!   - changedNVisu:         new NVisu, new Nodetype
!   - changedFV_Elems:      new distribution of FV/DG elements (only if changedStateFile==TRUE)
!   - changedWithDGOperator: different mode, with/without gradients
!   - changedDGonly:        the visualization of FV elements as DG elements was set or unset
!   - changedNCalc:         the polynomial degree used for calculations changed
!
! WORKFLOW:
! * The main steps are:
!   1. call InitFile for FV/DG distribution and mappings
!   2. read solution          (if changedStateFile or changedWithDGOperator or changedDGonly)
!   3. build mapping for BC sides that should be visualized (done after read solution since some
!      mesh infos are needed)
!   4. read Mesh              (if changedMeshFile)
!   5. compute UCalc          (if changedStateFile or changedVarNames or changedDGonly or changedNCalc)
!   6. convert to UVisu       (if changedStateFile or changedVarNames or changedNVisu or changedDGonly or changedNCalc)
!   7. build visu mesh        (if changedMeshFile  or changedNVisu or changedFV_Elems or changedDGonly)
!   5. - 7. are done seperately for surface variables if surface visualization is turned on
!   8. write VTK arrays       (always!)
!
!**********************************************************************************************
CALL SetStackSizeUnlimited()
#if USE_MPI
CALL InitMPI(mpi_comm_IN)
#else
CALL InitMPI()
#endif /*USE_MPI*/
CALL InitMPIInfo()

CALL FinalizeParameters()
! Read Varnames to visualize and build calc and visu dependencies
CALL prms%SetSection("posti")
CALL prms%CreateStringOption( "MeshFile"        , "Custom mesh file ")
CALL prms%CreateStringOption( "VarName"         , "Names of variables, which should be visualized.", multiple=.TRUE.)
CALL prms%CreateIntOption(    "NVisu"           , "Polynomial degree at which solution is sampled for visualization.")
CALL prms%CreateStringOption( "NodeTypeVisu"    , "NodeType for visualization. Visu, Gauss,Gauss-Lobatto,Visu_inner"    ,"VISU")
CALL prms%CreateStringOption( "BoundaryName"    , "Names of boundaries for surfaces, which should be visualized.", multiple=.TRUE.)

SWRITE (*,*) "READING FROM: ", TRIM(statefile)

changedStateFile      = .FALSE.
changedMeshFile       = .FALSE.
changedNVisu          = .FALSE.
changedVarNames       = .FALSE.
changedWithDGOperator = .FALSE.


IF (ISVALIDMESHFILE(statefile)) THEN ! visualize mesh
  SWRITE(*,*) "MeshFile Mode"
  MeshFileMode = .TRUE.
  MeshFile      = statefile
  nVar_State    = 0
  CALL VisualizeMesh(mpi_comm_IN,postifile,MeshFile)
ELSE IF (ISVALIDHDF5FILE(statefile)) THEN ! visualize state file
  SWRITE(*,*) "State Mode"
  MeshFileMode = .FALSE.
  ! initialize state file
  CALL visu_InitFile(mpi_comm_IN,statefile,postifile)

  ! read solution from state file (either direct or including a evaluation of the DG operator)
  IF (LEN_TRIM(prmfile).EQ.0) THEN
    changedPrmFile = .NOT.STRICMP(prmfile_old, ".piclas_posti.ini")
  ELSE
    changedPrmFile = (prmfile .NE. prmfile_old)
  END IF
  SWRITE (*,*) "changedStateFile     ", changedStateFile
  SWRITE (*,*) "changedMeshFile      ", changedMeshFile
  SWRITE (*,*) "changedNVisu         ", changedNVisu
  SWRITE (*,*) "changedVarNames      ", changedVarNames
  SWRITE (*,*) "changedPrmFile       ", changedPrmFile, TRIM(prmfile_old), " -> ", TRIM(prmfile)
  SWRITE (*,*) "changedBCnames       ", changedBCnames
  IF (changedStateFile.OR.changedWithDGOperator.OR.changedPrmFile) THEN
      CALL ReadState(mpi_comm_IN,prmfile,statefile)
  END IF

  ! build mappings of BC sides for surface visualization
  CALL Build_mapBCSides()

  !! ===== calc solution =====
  !IF (changedStateFile.OR.changedVarNames) THEN
  !  CALL CalcQuantities_DG()
  !END IF

  ! ===== convert solution to visu grid =====
  IF (changedStateFile.OR.changedVarNames.OR.changedNVisu) THEN
    CALL ConvertToVisu_DG()
  END IF

  ! convert generic data to visu grid
  IF (changedStateFile.OR.changedVarNames.OR.changedNVisu.OR.changedBCnames) THEN
    CALL ConvertToVisu_GenericData(mpi_comm_IN,statefile)
  END IF

  ! Convert coordinates to visu grid
  IF (changedMeshFile.OR.changedNVisu) THEN
    CALL BuildVisuCoords()
  END IF

END IF

MeshFile_old          = MeshFile
prmfile_old           = prmfile
statefile_old         = statefile
NVisu_old             = NVisu
!NCalc_old             = NCalc
nVar_State_old        = nVar_State
!withDGOperator_old    = withDGOperator
!DGonly_old            = DGonly
!Avg2D_old             = Avg2D
NodeTypeVisuPosti_old = NodeTypeVisuPosti
!NState_old            = PP_N

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(*,*) "Visu finished for state file: ", TRIM(statefile)
SWRITE(UNIT_StdOut,'(132("="))')
END SUBROUTINE visu

!===================================================================================================================================
!> Deallocate arrays used by visu.
!===================================================================================================================================
SUBROUTINE FinalizeVisu()
USE MOD_Globals
USE MOD_Visu_Vars
#if USE_FV
USE MOD_FV_Vars
#else
USE MOD_DG_Vars
#endif /*FV*/
USE MOD_Mesh_Vars, ONLY: Elem_xGP
IMPLICIT NONE
!===================================================================================================================================
SWRITE (*,*) "VISU FINALIZE"
prmfile_old = ""
statefile_old = ""
MeshFile = ""
MeshFile_old = ""
NodeTypeVisuPosti = "VISU"
NodeTypeVisuPosti_old = ""
NVisu     = -1
NVisu_old = -1
nVar_State_old = -1

SDEALLOCATE(mapDepToCalc)
SDEALLOCATE(mapAllVarsToVisuVars)
SDEALLOCATE(UCalc_DG)

SDEALLOCATE(mapDGElemsToAllElems)

SDEALLOCATE(CoordsVisu_DG)
SDEALLOCATE(UVisu_DG)
SDEALLOCATE(U)
SDEALLOCATE(Elem_xGP)
END SUBROUTINE FinalizeVisu

END MODULE MOD_Visu
