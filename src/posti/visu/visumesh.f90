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

!=================================================================================================================================
!> Routines to build the mesh for visualization.
!=================================================================================================================================
MODULE MOD_Posti_VisuMesh
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE BuildVisuCoords
  MODULE PROCEDURE BuildVisuCoords
END INTERFACE

INTERFACE VisualizeMesh
  MODULE PROCEDURE VisualizeMesh
END INTERFACE

PUBLIC:: BuildVisuCoords
PUBLIC:: VisualizeMesh

CONTAINS

!=================================================================================================================================
!> Converts the coordinates of the mesh to the visu-mesh.
!=================================================================================================================================
SUBROUTINE BuildVisuCoords()
! MODULES
USE ISO_C_BINDING
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars          ,ONLY: CoordsVisu_DG
USE MOD_Visu_Vars          ,ONLY: NodeTypeVisuPosti
USE MOD_Visu_Vars          ,ONLY: NVisu,nElems_DG,mapDGElemsToAllElems
USE MOD_Interpolation_Vars ,ONLY: NodeType
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_ChangeBasis        ,ONLY: ChangeBasis2D
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_Mesh_Vars          ,ONLY: Elem_xGP
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElem, iElem_DG
REAL,ALLOCATABLE   :: Vdm_N_NVisu(:,:)
!===================================================================================================================================
! Convert coordinates to visu grid
SWRITE (*,*) "[MESH] Convert coordinates to visu grid (DG)"
ALLOCATE(Vdm_N_NVisu(0:NVisu,0:PP_N))
CALL GetVandermonde(PP_N,NodeType,NVisu   ,NodeTypeVisuPosti  ,Vdm_N_NVisu   ,modal=.FALSE.)
! convert coords of DG elements
SDEALLOCATE(CoordsVisu_DG)
ALLOCATE(CoordsVisu_DG(3,0:NVisu,0:NVisu,0:NVisu,nElems_DG))
DO iElem_DG = 1,nElems_DG
  iElem = mapDGElemsToAllElems(iElem_DG)
  CALL ChangeBasis3D(3,PP_N,NVisu,Vdm_N_NVisu,Elem_xGP(:,:,:,:,iElem),CoordsVisu_DG(:,:,:,:,iElem_DG))
END DO
SDEALLOCATE(Vdm_N_NVisu)

END SUBROUTINE BuildVisuCoords


!=================================================================================================================================
!> Visualize mesh only
!> 1. read mesh
!> 2. BuildVisuCoords
!> 3. write mesh to VTK array
!> 4. set length of all other output arrays to zero
!=================================================================================================================================
SUBROUTINE VisualizeMesh(mpi_comm_IN,postifile,meshfile_in)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_ReadInTools   ,ONLY: prms,GETINT
USE MOD_ReadInTools   ,ONLY: FinalizeParameters
#if USE_MPI
USE MOD_MPI           ,ONLY: FinalizeMPI
#endif
USE MOD_Interpolation ,ONLY: DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_Mesh_Vars     ,ONLY: nElems,Ngeo
USE MOD_Mesh          ,ONLY: DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_VTK           ,ONLY: WriteCoordsToVTK_array
USE MOD_HDF5_Input    ,ONLY: ReadAttribute,File_ID,OpenDataFile,CloseDataFile
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)           :: mpi_comm_IN    
CHARACTER(LEN=255),INTENT(IN):: postifile
CHARACTER(LEN=255),INTENT(IN):: meshfile_in
! LOCAL VARIABLES
INTEGER             :: iElem
!===================================================================================================================================
#if USE_MPI
CALL FinalizeMPI()
#endif
CALL FinalizeMesh()
CALL FinalizeInterpolation()

CALL OpenDataFile(meshfile_in,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=mpi_comm_IN)
CALL ReadAttribute(File_ID,'Ngeo',1,IntScalar=Ngeo)
CALL CloseDataFile()

IF (LEN_TRIM(postifile).GT.0) THEN
  ! read options from parameter file
  CALL DefineParametersInterpolation()
  CALL DefineParametersMesh()
  CALL prms%SetSection("posti")
  CALL prms%CreateIntOption('NVisu', "Number of points at which solution is sampled for visualization.")
  CALL prms%read_options(postifile)
  NVisu = GETINT('NVisu') ! Degree of visualization basis
ELSE
  NVisu = 2*NGeo ! TODO: correct?
END IF


! read mesh
CALL InitInterpolation(Ngeo)
CALL InitMesh(meshMode=0, MeshFile_IN=meshfile_in)

! convert to visu grid
nElems_DG = nElems

SDEALLOCATE(mapDGElemsToAllElems)
ALLOCATE(mapDGElemsToAllElems(nElems))
DO iElem=1,nElems
  mapDGElemsToAllElems(iElem) = iElem
END DO
CALL BuildVisuCoords()
DEALLOCATE(mapDGElemsToAllElems)

CALL FinalizeInterpolation()
CALL FinalizeParameters()

END SUBROUTINE VisualizeMesh

END MODULE MOD_Posti_VisuMesh
