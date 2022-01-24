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
!> Contains routines that convert the calculated FV or DG quantities to the visualization grid.
!> The routines are split into surface and volume data. Also there is a routine that handels generic data like additional arrays
!> or data from non-state files.
!===================================================================================================================================
MODULE MOD_Posti_ConvertToVisu
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE ConvertToVisu_DG
  MODULE PROCEDURE ConvertToVisu_DG
END INTERFACE
PUBLIC:: ConvertToVisu_DG

INTERFACE ConvertToVisu_GenericData
  MODULE PROCEDURE ConvertToVisu_GenericData
END INTERFACE
PUBLIC:: ConvertToVisu_GenericData

CONTAINS

!===================================================================================================================================
!> Perform a ChangeBasis of the calculated volume DG quantities to the visualization grid.
!===================================================================================================================================
SUBROUTINE ConvertToVisu_DG()
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars          ,ONLY: nVarVisu,NodeTypeVisuPosti,nVarDep,NVisu
USE MOD_Visu_Vars          ,ONLY: mapAllVarsToVisuVars,mapDepToCalc
USE MOD_Visu_Vars          ,ONLY: nElems_DG,UCalc_DG,UVisu_DG
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_Interpolation_Vars ,ONLY: NodeType
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElem,iVar,iVarVisu,iVarCalc
REAL,ALLOCATABLE   :: Vdm_N_NVisu(:,:)                  ! Vandermonde from state to visualisation nodes
!===================================================================================================================================
SWRITE(*,*) "[DG] convert to visu grid"

! compute UVisu_DG
ALLOCATE(Vdm_N_NVisu(0:NVisu,0:PP_N))
CALL GetVandermonde(PP_N,NodeType,NVisu,NodeTypeVisuPosti,Vdm_N_NVisu,modal=.FALSE.)

! convert DG solution to UVisu_DG
SDEALLOCATE(UVisu_DG)
ALLOCATE(UVisu_DG(0:NVisu,0:NVisu,0:NVisu,nElems_DG,nVarVisu))
DO iVar=1,nVarDep
  IF (mapAllVarsToVisuVars(iVar).GT.0) THEN
    iVarCalc = mapDepToCalc(iVar)
    iVarVisu = mapAllVarsToVisuVars(iVar)
    DO iElem = 1,nElems_DG
      CALL ChangeBasis3D(1,PP_N,NVisu,Vdm_N_NVisu,UCalc_DG(:,:,:,iElem,iVarCalc),UVisu_DG(:,:,:,iElem,iVarVisu))
    END DO
  END IF
END DO

SDEALLOCATE(Vdm_N_NVisu)
END SUBROUTINE ConvertToVisu_DG

!===================================================================================================================================
!> This routine will read all variables that are not conservative or derived quantities and convert the ones that should be
!> visualized to the visu grid.
!> These variables include the additional data from the ElemData and FieldData datasetes as well as other datasets that are
!> present in the HDF5 file. The variables will be named DATASETNAME:VARIABLENAME if a attribute VarNames_DATASETNAME exist
!> where we can read the variable names from. If this  attribute does not exist, the name will be a generic DATASETNAME:1,2... .
!> For each dataset a new Vandermonde matrix is built to convert from the specific polynomial degree to the visu grid,
!> so the datasets are not limited to one polynomial degree. Either elementwise (2 dimensions) or pointwise (5 dimensions) datasets
!> are allowed.
!> The addtional variables will always be sorted AFTER the conservative or derived quantities.
!> If surface visualization is needed, the quantities will simply be prolonged to the surfaces.
!===================================================================================================================================
SUBROUTINE ConvertToVisu_GenericData(mpi_comm_IN,statefile)
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_IO_HDF5            ,ONLY: HSize
USE MOD_HDF5_Input         ,ONLY: File_ID,GetVarNames
USE MOD_HDF5_Input         ,ONLY: OpenDataFile,ReadArray,CloseDataFile,DatasetExists,ReadAttribute,GetDataSize
USE MOD_Mesh_Vars          ,ONLY: nElems,offsetElem
USE MOD_StringTools        ,ONLY: STRICMP,split_string
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D,ChangeBasis2D
USE MOD_Interpolation_Vars ,ONLY: NodeType
!USE MOD_ProlongToFace      ,ONLY: EvalElemFace
!USE MOD_Mappings           ,ONLY: buildMappings
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)           :: mpi_comm_IN
CHARACTER(LEN=255),INTENT(IN)  :: statefile   !< HDF5 state file
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iVarVisu,iElem_DG,iElem,iVarDataset,iVar,iVar2
INTEGER                        :: substring_count,nDims,nVal,nSize,nSizeZ
CHARACTER(LEN=255)             :: substrings(2),DatasetName,VariableName,DataSetOld
LOGICAL                        :: datasetFound,varnamesExist,datasetChanged
REAL,ALLOCATABLE               :: ElemData(:,:),FieldData(:,:,:,:,:)
REAL,ALLOCATABLE               :: Vdm_DG_Visu(:,:)
REAL,ALLOCATABLE               :: Uface_tmp(:,:,:),Uface(:,:)
CHARACTER(LEN=255),ALLOCATABLE :: DatasetVarNames(:)
INTEGER,ALLOCATABLE            :: S2V2(:,:,:,:,:)
!===================================================================================================================================
SWRITE(*,*) "Convert generic datasets to Visu grid"
! Open HDF5 file
CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=mpi_comm_IN)

DataSetOld = ''  ! Used to decide if arrays and Vandermonde matrix should be re-allocated

! Loop over all generic variables that should be visualized - sorted after the dependant variables
DO iVar=nVarDep+1,nVarAll
  ! Check if this variable should be visualized
  IF ((mapAllVarsToVisuVars(iVar)).GT.0) THEN
    ! The format of the generic data varnames is DATASETNAME:VARIABLENAME - split into DATASETNAME and VARIABLENAME
    CALL split_string(TRIM(VarnamesAll(iVar)),':',substrings,substring_count)
    ! If we find more than one substring, the variable is additional data
    IF (substring_count.GT.1) THEN
      ! Store dataset and variable name
      DatasetName =  TRIM(substrings(1))
      VariableName = TRIM(substrings(2))

      ! Check if we have a new dataset
      datasetChanged = .NOT.STRICMP(TRIM(DataSetName),TRIM(DataSetOld))

      SWRITE(*,*) "Convert variable ",TRIM(VariableName)," from dataset ", TRIM(DatasetName)

      ! Get metadata if dataset changed
      IF (datasetChanged) THEN
        ! Try to open the dataset
        CALL DatasetExists(File_ID,TRIM(DatasetName),datasetFound)
        ! Abort if the dataset was not found
        IF (.NOT.datasetFound)  THEN
          CALL CloseDataFile()
          CALL Abort(__STAMP__,'Dataset '//TRIM(DatasetName)//' does not exist.')
        END IF
        ! Get dimensions of the dataset and store number of variables as well as size of array
        CALL GetDataSize(File_ID,TRIM(DatasetName),nDims,HSize)
        nVal   = INT(HSize(1))
        nSize  = INT(HSize(2))
        DEALLOCATE(HSize)
        nSizeZ = nSize
        SDEALLOCATE(DataSetVarNames)
        CALL GetVarNames("VarNames_"//TRIM(DatasetName),DatasetVarNames,varnamesExist)
      END IF
      SWRITE (*,*) "varnamesExist", varnamesExist

      iVarDataset = 0
      ! loop over all varnames
      DO iVar2=nVarDep+1,nVarAll
        CALL split_string(TRIM(VarnamesAll(iVar2)),':',substrings,substring_count)
        IF (substring_count.GT.1) THEN
          ! if dataset is the same increase variable index inside dataset
          IF (STRICMP(substrings(1),DatasetName)) THEN
            iVarDataset = iVarDataset+1
            ! exit if variablename found
            IF (STRICMP(substrings(2),VariableName)) EXIT
          END IF
        END IF
      END DO

      ! Read in the data if we have a new dataset. Also allocate Vandermonde matrix used in conversion to visu grid.
      IF (datasetChanged.AND.(iVarDataset.GT.0)) THEN
        ASSOCIATE (&
                   nVal       => INT(nVal,IK)       ,&
                   nSizeZ     => INT(nSizeZ,IK)     ,&
                   nSize      => INT(nSize,IK)      ,&
                   nElems     => INT(nElems,IK)     ,&
                   offsetElem => INT(offsetElem,IK)  &
                   )
          SELECT CASE(nDims)
          CASE(2) ! Elementwise data
            ! Allocate array and read dataset
            SDEALLOCATE(ElemData)
            ALLOCATE(ElemData(nVal,nElems))
            CALL ReadArray(TRIM(DatasetName),2,(/nVal,nElems/),offsetElem,2,RealArray=ElemData)
          CASE(5) ! Pointwise data
            ! Allocate array and read dataset
            SDEALLOCATE(FieldData)
            ALLOCATE(FieldData(nVal,nSize,nSize,nSizeZ,nElems))
            CALL ReadArray(TRIM(DatasetName),5,(/nVal,nSize,nSize,nSizeZ,nElems/),offsetElem,5,RealArray=FieldData)
            ! Get Vandermonde matrix used to convert to the visu grid
            SDEALLOCATE(Vdm_DG_Visu)
            ALLOCATE(Vdm_DG_Visu(0:NVisu,0:nSize-1))
            CALL GetVandermonde(INT(nSize,4)-1,NodeType,NVisu,NodeTypeVisuPosti,Vdm_DG_Visu,modal=.FALSE.)
          CASE DEFAULT
            CALL Abort(__STAMP__,'Dataset '//TRIM(DatasetName)//' does not have 2 or 5 dimensions.')
          END SELECT
        END ASSOCIATE

        ! Store current name of dataset
        DataSetOld = TRIM(DatasetName)
      END IF ! New dataset

      ! Get index of visu array that we should write to
      iVarVisu= mapAllVarsToVisuVars(iVar)
      ! Convert the generic data to visu grid
      SELECT CASE(nDims)
      CASE(2) ! Elementwise data
          ! Simply write the elementwise data to all visu points
          DO iElem_DG=1,nElems_DG
            iElem = mapDGElemsToAllElems(iElem_DG)
            UVisu_DG(:,:,:,iElem_DG,iVarVisu) = ElemData(iVarDataset,iElem)
          END DO
      CASE(5) ! Pointwise data
          ! Perform changebasis to visu grid
          DO iElem_DG=1,nElems_DG
            iElem = mapDGElemsToAllElems(iElem_DG)
            CALL ChangeBasis3D(1,nSize-1,NVisu,Vdm_DG_Visu,FieldData(iVarDataset,:,:,:,iElem),UVisu_DG(:,:,:,iElem_DG,iVarVisu))
          END DO
      END SELECT

    END IF ! substring_count.GT.1
  END IF ! mapAllVarsToVisuVars(iVar).GT.0

END DO !iVar=1,

! Close HDF5 file
CALL CloseDataFile()

! Cleanup of allocatable arrays
SDEALLOCATE(ElemData)
SDEALLOCATE(FieldData)
SDEALLOCATE(Vdm_DG_Visu)
SDEALLOCATE(DataSetVarNames)

SDEALLOCATE(S2V2)
SDEALLOCATE(Uface)
SDEALLOCATE(Uface_tmp)

END SUBROUTINE ConvertToVisu_GenericData

END MODULE MOD_Posti_ConvertToVisu
