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

MODULE MOD_HDF5_Output_ElemData
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MOD_io_HDF5
USE MOD_HDF5_output
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

#if USE_MPI
PUBLIC :: WriteMyInvisibleRankToHDF5
#endif /*USE_MPI*/
#if defined(PARTICLES)
PUBLIC :: WriteLostRotPeriodicSidesToHDF5
#endif /*defined(PARTICLES)*/

PUBLIC :: WriteAdditionalElemData
!===================================================================================================================================

CONTAINS


SUBROUTINE WriteAdditionalElemData(FileName,ElemList)
!===================================================================================================================================
!> Write additional data for analyze purpose to HDF5.
!> The data is taken from a lists, containing either pointers to data arrays or pointers
!> to functions to generate the data, along with the respective varnames.
!>
!> Two options are available:
!>    1. WriteAdditionalElemData:
!>       Element-wise scalar data, e.g. the timestep or indicators.
!>       The data is collected in a single array and written out in one step.
!>       DO NOT MISUSE NODAL DATA FOR THIS! IT WILL DRASTICALLY INCREASE FILE SIZE AND SLOW DOWN IO!
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars        ,ONLY: offsetElem,nGlobalElems,nElems
USE MOD_LoadBalance_Vars ,ONLY: ElemTime,NullifyElemTime
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)        :: FileName
TYPE(tElementOut),POINTER,INTENT(IN) :: ElemList !< Linked list of arrays to write to file
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
REAL,ALLOCATABLE               :: ElemData(:,:)
INTEGER                        :: nVar,iElem
TYPE(tElementOut),POINTER      :: e
!===================================================================================================================================

IF(.NOT. ASSOCIATED(ElemList)) RETURN

! Count the additional variables
nVar = 0
e=>ElemList
DO WHILE(ASSOCIATED(e))
  nVar=nVar+1
  e=>e%next
END DO

! Allocate variable names and data array
ALLOCATE(StrVarNames(nVar))
ALLOCATE(ElemData(nVar,nElems))

! Fill the arrays
nVar = 0
e=>ElemList
DO WHILE(ASSOCIATED(e))
  nVar=nVar+1
  StrVarNames(nVar)=e%VarName
  IF(ASSOCIATED(e%RealArray))    ElemData(nVar,:)=e%RealArray(1:nElems)
  IF(ASSOCIATED(e%RealScalar))   ElemData(nVar,:)=e%RealScalar
  IF(ASSOCIATED(e%IntArray))     ElemData(nVar,:)=REAL(e%IntArray(1:nElems))
  IF(ASSOCIATED(e%IntScalar))    ElemData(nVar,:)=REAL(e%IntScalar)
  IF(ASSOCIATED(e%LongIntArray)) ElemData(nVar,:)=REAL(e%LongIntArray(1:nElems))
  IF(ASSOCIATED(e%LogArray)) THEN
    DO iElem=1,nElems
      IF(e%LogArray(iElem))THEN
        ElemData(nVar,iElem)=1.
      ELSE
        ElemData(nVar,iElem)=0.
      END IF
    END DO ! iElem=1,PP_nElems
  END IF
  IF(ASSOCIATED(e%eval))       CALL e%eval(ElemData(nVar,:)) ! function fills elemdata
  e=>e%next
END DO

IF(MPIRoot)THEN
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteAttributeToHDF5(File_ID,'VarNamesAdd',nVar,StrArray=StrVarNames)
  CALL CloseDataFile()
END IF

ASSOCIATE (&
      nVar         => INT(nVar,IK)         ,&
      nGlobalElems => INT(nGlobalElems,IK) ,&
      PP_nElems    => INT(PP_nElems,IK)    ,&
      offsetElem   => INT(offsetElem,IK)   )
  CALL GatheredWriteArray(FileName,create = .FALSE.,&
                          DataSetName     = 'ElemData', rank = 2,  &
                          nValGlobal      = (/nVar,nGlobalElems/),&
                          nVal            = (/nVar,PP_nElems   /),&
                          offset          = (/0_IK,offsetElem  /),&
                          collective      = .TRUE.,RealArray = ElemData)
END ASSOCIATE
DEALLOCATE(ElemData,StrVarNames)

! Check if ElemTime is to be nullified (required after user-restart)
! After writing the old ElemTime values to disk, the array must be nullified (because they correspond to the restart file, which
! might have been created with a totally different processor number and distribution)
IF(NullifyElemTime) ElemTime=0.

END SUBROUTINE WriteAdditionalElemData


#if USE_MPI
!===================================================================================================================================
!> write LostRotPeriodicSides field to HDF5 file (as ElemData container)
!===================================================================================================================================
SUBROUTINE WriteMyInvisibleRankToHDF5()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars    ,ONLY: MeshFile
USE MOD_Globals_Vars ,ONLY: ProjectName
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER   :: N_variables=1
CHARACTER(LEN=255),ALLOCATABLE  :: StrVarNames(:)
CHARACTER(LEN=255)  :: FileName
REAL                :: StartT,EndT
REAL                :: OutputTime
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE MyInvisibleRank TO HDF5 FILE...'
GETTIME(StartT)
OutputTime=0.0
ALLOCATE(StrVarNames(1:N_variables))
StrVarNames(1)='dummy'
! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(ProjectName)//'_MyInvisibleRank.h5'
IF(MPIRoot) CALL GenerateFileSkeleton('MyInvisibleRank',N_variables,StrVarNames,TRIM(MeshFile),OutputTime,FileNameIn=FileName)
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif

! Write all 'ElemData' arrays to a single container in the state.h5 file
CALL WriteAdditionalElemData(FileName,ElementOut)

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
END SUBROUTINE WriteMyInvisibleRankToHDF5
#endif /*USE_MPI*/


#if defined(PARTICLES)
!===================================================================================================================================
!> write LostRotPeriodicSides field to HDF5 file (as ElemData container)
!===================================================================================================================================
SUBROUTINE WriteLostRotPeriodicSidesToHDF5()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars    ,ONLY: MeshFile
USE MOD_Globals_Vars ,ONLY: ProjectName
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER   :: N_variables=1
CHARACTER(LEN=255),ALLOCATABLE  :: StrVarNames(:)
CHARACTER(LEN=255)  :: FileName
REAL                :: StartT,EndT
REAL                :: OutputTime
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE LostRotPeriodicSides TO HDF5 FILE...'
GETTIME(StartT)
OutputTime=0.0
ALLOCATE(StrVarNames(1:N_variables))
StrVarNames(1)='dummy'
! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(ProjectName)//'_LostRotPeriodicSides.h5'
IF(MPIRoot) CALL GenerateFileSkeleton('LostRotPeriodicSides',N_variables,StrVarNames,TRIM(MeshFile),OutputTime,FileNameIn=FileName)
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif

! Write all 'ElemData' arrays to a single container in the state.h5 file
CALL WriteAdditionalElemData(FileName,ElementOut)

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
END SUBROUTINE WriteLostRotPeriodicSidesToHDF5
#endif /*defined(PARTICLES)*/


END MODULE MOD_HDF5_Output_ElemData
