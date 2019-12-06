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

MODULE MOD_Tecplot
!===================================================================================================================================
! Module for generic data output in Tecplot fromat
!
! WARNING: WriteDataToTecplot works only for POSTPROCESSING
!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE WriteDataToTecplot
  MODULE PROCEDURE WriteDataToTecplot
END INTERFACE

INTERFACE WriteDataToTecplotBinary
  MODULE PROCEDURE WriteDataToTecplotBinary
END INTERFACE

PUBLIC::WriteDataToTecplotBinary
PUBLIC::WriteDataToTecplot
!===================================================================================================================================

CONTAINS

SUBROUTINE WriteDataToTecplot(NPlot,nElems,nVal,nValMean,VarNames,Coord,Value,FileString,MeanValue)
!===================================================================================================================================
! Subroutine to write point data to Tecplot format
!===================================================================================================================================
! MODULES
!USE MOD_PreProc
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: NPlot                   ! Number of output points .EQ. NAnalyze
INTEGER,INTENT(IN)            :: nElems                  ! Number of Elements
INTEGER,INTENT(IN)            :: nVal                    ! Number of nodal output variables
INTEGER,INTENT(IN)            :: nValMean                ! Number of element mean output variables
CHARACTER(LEN=32),INTENT(IN)  :: VarNames(nVal+nValMean) ! Names of all variables that will be written out
REAL,INTENT(IN)               :: Coord(3,0:NPlot,0:NPlot,0:NPlot,nElems)      ! CoordsVector
REAL,INTENT(IN)               :: Value(1:nVal,0:NPlot,0:NPlot,0:NPlot,nElems) ! Statevector
REAL,INTENT(IN),OPTIONAL      :: MeanValue(1:nValMean,nElems) ! Vector of indicator values
CHARACTER(LEN=*),INTENT(IN)   :: FileString              ! Output file name
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,iVal,Offset
CHARACTER(LEN=255) :: Format_nVal
CHARACTER(LEN=255) :: Format_Title
CHARACTER(LEN=35)  :: VarString
INTEGER            :: iElem
INTEGER            :: NodeIDElem,nPlot_p1_2,nPlot_p1_3
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')" WRITE DATA TO TECPLOT ASCII FILE..."
IF((nValMean.GT.0) .AND. .NOT. PRESENT(MeanValue)) &
  CALL abort(&
  __STAMP__&
  ,'Mean Value >0, but no mean value specified ',999,999.)
!assemble format strings
WRITE(Format_nVal,'(A1,I2,A17)')'(',nVal+nValMean,'(1X,E21.10))'
Format_Title(1:51)='VARIABLES="CoordinateX","CoordinateY","CoordinateZ"'
Offset = 52
DO iVal=1,nVal+nValMean
  WRITE(VarString,'(A2,A,A1)')',"',TRIM(VarNames(iVal)),'"'
  Format_Title(Offset:Offset+LEN(TRIM(VarString)))=TRIM(VarString)
  Offset=Offset+LEN(TRIM(VarString))
END DO

!connected 3D FEM data
NPlot_p1_2=(NPlot+1)**2
NPlot_p1_3=(NPlot+1)**3

OPEN(44,FILE=TRIM(FileString),Status="REPLACE")
WRITE(44,'(A)')Format_Title(1:Offset-1)
WRITE(44,'(A,I8,A,I8,A)')'ZONE T="",DATAPACKING=POINT, NODES=',nElems*NPlot_p1_3,', ELEMENTS=',nElems*(NPlot)**3,',ZONETYPE=FEBRICK'
IF(nValMean.GT.0)THEN
  DO iElem=1,nElems
    DO k=0,NPlot
      DO j=0,NPlot
        DO i=0,NPlot
          WRITE(44,Format_nVal)&
          Coord(:,i,j,k,iElem),Value(:,i,j,k,iElem),MeanValue(:,iElem)
        END DO
      END DO
    END DO
  END DO
ELSE
  DO iElem=1,nElems
    DO k=0,NPlot
      DO j=0,NPlot
        DO i=0,NPlot
          WRITE(44,Format_nVal)&
          Coord(:,i,j,k,iElem),Value(:,i,j,k,iElem)
        END DO
      END DO
    END DO
  END DO
END IF !meanVal>0
!element connectivity
NodeIDElem=0
DO iElem=1,nElems
  DO k=1,NPlot
    DO j=1,NPlot
      DO i=1,NPlot
        !visuHexaElem
        WRITE(44,'(8(I8,1X))')&
          NodeIDElem+i+  (j-1)*(NPlot+1)+(k-1)*NPlot_p1_2,      & !P1(CGNS=tecplot standard)
          NodeIDElem+i+1+(j-1)*(NPlot+1)+(k-1)*NPlot_p1_2,      & !P2
          NodeIDElem+i+1+ j   *(NPlot+1)+(k-1)*NPlot_p1_2,      & !P3
          NodeIDElem+i+   j   *(NPlot+1)+(k-1)*NPlot_p1_2,      & !P4
          NodeIDElem+i+  (j-1)*(NPlot+1)+ k   *NPlot_p1_2,      & !P5
          NodeIDElem+i+1+(j-1)*(NPlot+1)+ k   *NPlot_p1_2,      & !P6
          NodeIDElem+i+1+ j   *(NPlot+1)+ k   *NPlot_p1_2,      & !P7
          NodeIDElem+i+   j   *(NPlot+1)+ k   *NPlot_p1_2         !P8
      END DO
    END DO
  END DO
  NodeIDElem=NodeIDElem+NPlot_p1_3
END DO

CLOSE(44)
SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END SUBROUTINE WriteDataToTecplot



SUBROUTINE WriteDataToTecplotBinary(NPlot,nElems,nVal,nValMean,VarNames,Coord,Value,FileString,MeanValue)
!===================================================================================================================================
! Subroutine to write point data to Tecplot format
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: NPlot                   ! Number of output points .EQ. NAnalyze
INTEGER,INTENT(IN)            :: nElems                  ! Number of Elements
INTEGER,INTENT(IN)            :: nVal                    ! Number of nodal output variables
INTEGER,INTENT(IN)            :: nValMean                ! Number of element mean output variables
CHARACTER(LEN=32),INTENT(IN)  :: VarNames(nVal+nValMean) ! Names of all variables that will be written out
REAL,INTENT(IN)               :: Coord(3,0:NPlot,0:NPlot,0:NPlot,nElems)      ! CoordsVector
REAL,INTENT(IN)               :: Value(1:nVal,0:NPlot,0:NPlot,0:NPlot,nElems) ! Statevector
REAL,INTENT(IN),OPTIONAL      :: MeanValue(1:nValMean,nElems) ! Vector of indicator values
CHARACTER(LEN=*),INTENT(IN)   :: FileString              ! Output file name
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,iVal,Offset,iVar
INTEGER            :: TECINI112,TECZNE112,TECNOD112,TECDAT112,TECEND112
INTEGER            :: array_0(3+nVal+nValMean),array_1(3+nVal+nValMean)
INTEGER            :: Vertex(8,(NPlot+1)**3*nElems)
INTEGER            :: iStat,NodeID
INTEGER            :: iElem
INTEGER            :: NodeIDElem,nPlot_p1_2,nPlot_p1_3
REAL,ALLOCATABLE   :: MeanValueExtended(:,:,:,:,:)
LOGICAL            :: nValCheck
CHARACTER(LEN=255) :: Format_Title
CHARACTER(LEN=35)  :: VarString
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')" WRITE DATA TO TECPLOT BINARY FILE..."

IF(nValMean.GT.0 .AND. .NOT. PRESENT(MeanValue)) &
  CALL abort(&
  __STAMP__&
  ,'Mean Value >0, but no mean value specified',999,999.)
!assemble format strings
Format_Title(1:35)='CoordinateX,CoordinateY,CoordinateZ'
Offset = 36
DO iVal=1,nVal+nValMean
  WRITE(VarString,'(A,A)')',',TRIM(VarNames(iVal))
  Format_Title(Offset:Offset+LEN(TRIM(VarString)))=TRIM(VarString)
  Offset=Offset+LEN(TRIM(VarString))
END DO

!connected 3D FEM data
NPlot_p1_2=(NPlot+1)**2
NPlot_p1_3=(NPlot+1)**3
IF(nValMean.GT.0)THEN
  nValCheck = .true.
END IF
IF(nValCheck) THEN
  ALLOCATE(MeanValueExtended(1:nValMean,0:NPlot,0:NPlot,0:NPlot,1:nElems))
  MeanValueExtended = 0.
  DO iElem=1,nElems
    DO iVar=1,nValMean
      MeanValueExtended(iVar,:,:,:,iElem) = MeanValue(iVar,iElem)
    END DO     ! iVar
  END DO ! iElem
END IF

iStat = TECINI112(                  &
        ""//char(0),                    &       ! Title
        TRIM(Format_Title)//char(0),    &       ! Variables
        TRIM(FileString)//char(0),      &       ! FName
        "."//char(0),                   &       ! ScratchDir
        0,                              &       ! FileType
        0,                              &       ! Debug
        1                               )       ! VIsDouble

array_0 = 0
array_1 = 1

iStat = TECZNE112(      &
        ""//char(0),        &       ! ZoneTitle
        5,                  &       ! ZoneType
        nElems*NPlot_p1_3,  &       ! IMxOrNumPts
        nElems*(NPlot)**3,  &       ! JMxOrNumElements
        0,                  &       ! KMxOrNumFaces
        0,                  &       ! ICellMax
        0,                  &       ! JCellMax
        0,                  &       ! KCellMax
        0.,                 &       ! SolutionTime
        0,                  &       ! StrandID
        0,                  &       ! ParentZone
        1,                  &       ! IsBlock
        0,                  &       ! NumFaceConnections
        0,                  &       ! FaceNeighborMode
        0,                  &       ! TotalNumFaceNodes
        0,                  &       ! NumConnectedBoundaryFaces
        0,                  &       ! TotalNumBoundaryConnections
        array_0,            &       ! PassiveVarList
        array_1,            &       ! ValueLocation
        array_0,            &       ! ShareVarFromZone
        0                   )       ! ShareConnectivityFromZone

DO iVar=1,3
  iStat = TECDAT112(                  &
          nElems*NPlot_p1_3,          &       ! N
          Coord(iVar,:,:,:,:),        &       ! Data
          1                           )       ! IsDouble
END DO       ! iVar

DO iVar=1,nVal
  iStat = TECDAT112(                  &
          nElems*NPlot_p1_3,          &       ! N
          Value(iVar,:,:,:,:),        &       ! Data
          1                           )       ! IsDouble
END DO       ! iVar

IF(nValCheck)THEN
  DO iVar=1,nValMean
    iStat = TECDAT112(                              &
            nElems*NPlot_p1_3,                      &       ! N
            MeanValueExtended(iVar,:,:,:,:),        &       ! Data
            1                                       )       ! IsDouble
  END DO       ! iVar
  DEALLOCATE(MeanValueExtended)
END IF

!element connectivity
NodeID = 0
NodeIDElem = 0
DO iElem=1,nElems
  DO k=1,NPlot
    DO j=1,NPlot
      DO i=1,NPlot
        NodeID = NodeID+1
        !
        Vertex(:,NodeID) = (/                                   &
          NodeIDElem+i+  (j-1)*(NPlot+1)+(k-1)*NPlot_p1_2,      & !P1(CGNS=tecplot standard)
          NodeIDElem+i+1+(j-1)*(NPlot+1)+(k-1)*NPlot_p1_2,      & !P2
          NodeIDElem+i+1+ j   *(NPlot+1)+(k-1)*NPlot_p1_2,      & !P3
          NodeIDElem+i+   j   *(NPlot+1)+(k-1)*NPlot_p1_2,      & !P4
          NodeIDElem+i+  (j-1)*(NPlot+1)+ k   *NPlot_p1_2,      & !P5
          NodeIDElem+i+1+(j-1)*(NPlot+1)+ k   *NPlot_p1_2,      & !P6
          NodeIDElem+i+1+ j   *(NPlot+1)+ k   *NPlot_p1_2,      & !P7
          NodeIDElem+i+   j   *(NPlot+1)+ k   *NPlot_p1_2       /) !P8
      END DO
    END DO
  END DO
  !
  NodeIDElem=NodeIDElem+NPlot_p1_3
END DO

iStat = TECNOD112(Vertex)
iStat = TECEND112()

SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END SUBROUTINE WriteDataToTecplotBinary

END MODULE MOD_Tecplot
