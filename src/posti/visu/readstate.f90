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
!> Contains routines to read in the different states
!===================================================================================================================================
MODULE MOD_Posti_ReadState
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE ReadState
  MODULE PROCEDURE ReadState
END INTERFACE

PUBLIC:: ReadState

CONTAINS

!===================================================================================================================================
!> This routine will read in the current state from the statefile. Will call one of two routines: ReadStateWithoutGradients if no
!> gradients have to be visualized or calculated and so no DG operator call is necessary, or ReadStateAndGradients if gradients
!> are needed and we need to calculate the DG operator once.
!> If the DG operator has to be called, we need some parameters. Either a seperate parameter file is passed,  then this one will 
!> be used, or we try to extract the parameter file from the userblock.
!> If both fails and we need to compute the DG operator, the program will abort.
!> If the DG operator should not be called and no parameter file (seperate or from userblock) can be found,
!> we specify PP_N from the state file as our polynomial degree (later needed by InitInterpolation).
!===================================================================================================================================
SUBROUTINE ReadState(mpi_comm_IN,prmfile,statefile)
USE MOD_Globals
USE MOD_PreProc
USE MOD_ReadInTools ,ONLY:ExtractParameterFile
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)               :: mpi_comm_IN
CHARACTER(LEN=255),INTENT(INOUT) :: prmfile      !< FLEXI parameter file, used if DG operator is called
CHARACTER(LEN=255),INTENT(IN)    :: statefile    !< HDF5 state file
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                          :: userblockFound
!===================================================================================================================================
userblockFound = .TRUE. ! Set to true to later test for existing parameters either form userblock or from seperate file
IF (LEN_TRIM(prmfile).EQ.0) THEN ! No seperate parameter file has been given
  ! Try to extract parameter file 
  prmfile = ".piclas_posti.ini"
  CALL ExtractParameterFile(statefile,prmfile,userblockFound)
  ! Only abort if we need some parameters to call the DG operator
  IF (.NOT.userblockFound) THEN
    CALL CollectiveStop(__STAMP__, "No userblock found in state file '"//TRIM(statefile)//"' and no parameter file specified.")
  END IF
END IF
CALL ReadStateWithoutGradients(mpi_comm_IN,prmfile,statefile)
END SUBROUTINE ReadState

!===================================================================================================================================
!> Read 'U' directly from a state file.
!===================================================================================================================================
SUBROUTINE ReadStateWithoutGradients(mpi_comm_IN,prmfile,statefile)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_MPI,                 ONLY: DefineParametersMPI
#if USE_MPI
USE MOD_MPI,                 ONLY: FinalizeMPI
#endif
USE MOD_IO_HDF5,             ONLY: DefineParametersIO,InitIOHDF5
USE MOD_Mesh                ,ONLY: DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_ReadInTools         ,ONLY: prms
USE MOD_ReadInTools         ,ONLY: FinalizeParameters
USE MOD_Mesh_Vars           ,ONLY: nElems,offsetElem
USE MOD_HDF5_Input,          ONLY: OpenDataFile,ReadArray,CloseDataFile
USE MOD_DG_Vars             ,ONLY: U
USE MOD_Interpolation       ,ONLY: DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)           :: mpi_comm_IN
CHARACTER(LEN=255),INTENT(IN):: prmfile       !< FLEXI parameter file, used if DG operator is called
CHARACTER(LEN=255),INTENT(IN):: statefile     !< HDF5 state file
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: meshMode_loc
LOGICAL           :: changedMeshMode
!===================================================================================================================================
CALL FinalizeInterpolation()

! Some features require normal vectors or metrics
meshMode_loc = 0 ! Minimal mesh init
IF (doSurfVisu)  meshMode_loc = MAX(meshMode_loc,2)

! check if the mesh mode has changed from the last time
changedMeshMode = (meshMode_loc.NE.meshMode_old)


#if USE_MPI
IF ((changedMeshFile).OR.(changedMeshMode)) THEN
  CALL FinalizeMPI()
END IF
#endif

! read options from parameter file
CALL FinalizeParameters()
CALL DefineParametersMPI()
CALL DefineParametersIO()
CALL DefineParametersInterpolation()
CALL DefineParametersMesh()
CALL prms%read_options(prmfile)

! Initialization of I/O routines
CALL InitIOHDF5()

CALL InitInterpolation()

! Call mesh init if the mesh file changed or we need a different mesh mode
IF ((changedMeshFile).OR.(changedMeshMode)) THEN
  CALL FinalizeMesh()
  CALL InitMesh(meshMode=meshMode_loc,MeshFile_IN=MeshFile)
END IF

! save old mesh mode for future comparison
meshMode_old = meshMode_loc

SDEALLOCATE(U)
ALLOCATE(U(1:nVar_State,0:PP_N,0:PP_N,0:PP_N,nElems))
CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=mpi_comm_IN)
ASSOCIATE (&
      nVar_State => INT(nVar_State,IK) ,&
      N          => INT(PP_N,IK)       ,&
      nElems     => INT(nElems,IK)     ,&
      offsetElem => INT(offsetElem,IK)  &
      )
  CALL ReadArray('DG_Solution',5,(/nVar_State,N+1_IK,N+1_IK,N+1_IK,nElems/),offsetElem,5,RealArray=U)
END ASSOCIATE
CALL CloseDataFile()

CALL FinalizeParameters()
END SUBROUTINE ReadStateWithoutGradients

END MODULE MOD_Posti_ReadState
