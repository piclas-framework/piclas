#include "boltzplatz.h"
#ifdef PARTICLES

MODULE MOD_OutPutVTK
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

PUBLIC::WriteDataToVTK, WriteDataToVTKBin
!===================================================================================================================================

CONTAINS

SUBROUTINE WriteDataToVTK(NPlot,nElems,nVal,VarNames,Coord,Value,FileString)
!===================================================================================================================================
! Subroutine to write point data to Tecplot format
!===================================================================================================================================
! MODULES
!USE MOD_PreProc
USE MOD_Globals
#ifdef MPI
   USE MOD_part_MPI_Vars, ONLY : PMPIVAR
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: NPlot                   ! Number of output points .EQ. NAnalyze
INTEGER,INTENT(IN)            :: nElems                  ! Number of Elements
INTEGER,INTENT(IN)            :: nVal                    ! Number of nodal output variables
CHARACTER(LEN=32),INTENT(IN)  :: VarNames(nVal) ! Names of all variables that will be written out
REAL,INTENT(IN)               :: Coord(3,0:NPlot,0:NPlot,0:NPlot,nElems)      ! CoordsVector 
REAL,INTENT(IN)               :: Value(1:nVal,0:NPlot,0:NPlot,0:NPlot,nElems) ! Statevector 
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
CHARACTER(LEN=255)            :: myOutputFile
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')" WRITE DATA TO VTK FILE..."  
  OPEN(1112,FILE=FileString,STATUS='replace')
  WRITE(1112,'(A)')'# vtk DataFile Version 2.0'
  WRITE(1112,'(A)')'Debug Mesh '
  WRITE(1112,'(A)')'ASCII'
  WRITE(1112,'(A)')'DATASET UNSTRUCTURED_GRID'
  WRITE(1112,'(A)')''
  WRITE(1112,'(A,I0,A)')'POINTS ',(NPlot+1)**3*nElems,' FLOAT'
  DO iElem=1,nElems
    DO k=0,NPlot
      DO j=0,NPlot
        DO i=0,NPlot
          WRITE(1112,*) Coord(:,i,j,k,iElem)
        END DO 
      END DO 
    END DO 
  END DO 
  WRITE(1112,*)''
  WRITE(1112,'(A,I0,1X,I0)')'CELLS ',nElems*NPlot**3,9*(nElems*NPlot**3)
  NodeIDElem = 0
  NPlot_p1_2=(NPlot+1)**2
  NPlot_p1_3=(NPlot+1)**3
  DO iElem=1,nElems
    DO k=1,NPlot
      DO j=1,NPlot
        DO i=1,NPlot
          WRITE(1112,'(I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0)')8 &
          ,NodeIDElem - 1 +i+  (j-1)*(NPlot+1)+(k-1)*NPlot_p1_2 &  !P1(CGNS=tecplot standard)
          ,NodeIDElem - 1 +i+1+(j-1)*(NPlot+1)+(k-1)*NPlot_p1_2 &  !P2
          ,NodeIDElem - 1 +i+1+ j   *(NPlot+1)+(k-1)*NPlot_p1_2 &  !P3
          ,NodeIDElem - 1 +i+   j   *(NPlot+1)+(k-1)*NPlot_p1_2 & !P4
          ,NodeIDElem - 1 +i+  (j-1)*(NPlot+1)+ k   *NPlot_p1_2 &  !P5
          ,NodeIDElem - 1 +i+1+(j-1)*(NPlot+1)+ k   *NPlot_p1_2 &  !P6
          ,NodeIDElem - 1 +i+1+ j   *(NPlot+1)+ k   *NPlot_p1_2 &  !P7
          ,NodeIDElem - 1 +i+   j   *(NPlot+1)+ k   *NPlot_p1_2   !P8         
        END DO 
      END DO 
    END DO 
    NodeIDElem=NodeIDElem+NPlot_p1_3
  END DO 
  WRITE(1112,*)''
  WRITE(1112,'(A,I0)')'CELL_TYPES ',nElems*NPlot**3
  DO iElem=1,nElems*NPlot**3
    WRITE(1112,'(1X,I0)',ADVANCE="NO")12
  END DO  
  WRITE(1112,*)''
  WRITE(1112,*)''
  WRITE(1112,'(A,I0)')'POINT_DATA ',nElems*(NPlot+1)**3

  DO iVal = 1 , nVal
    WRITE(1112,'(A)') 'SCALARS ', TRIM(VarNames(iVal)), ' FLOAT'
    WRITE(1112,'(A)')'LOOKUP_TABLE default'
    DO iElem=1,nElems
      DO k=0,NPlot
        DO j=0,NPlot
          DO i=0,NPlot
            WRITE(1112,'(E22.15)') Value(iVal,i,j,k,iElem)
          END DO 
        END DO 
      END DO 
    END DO 
  END DO
  CLOSE(1112)

END SUBROUTINE WriteDataToVTK

SUBROUTINE WriteDataToVTKBin(NPlot,nElems,nVal,VarNames,Coord,Value,FileString)
!===================================================================================================================================
! Subroutine to write point data to Tecplot format
!===================================================================================================================================
! MODULES
!USE MOD_PreProc
USE MOD_Globals
#ifdef MPI
   USE MOD_part_MPI_Vars, ONLY : PMPIVAR
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: NPlot                   ! Number of output points .EQ. NAnalyze
INTEGER,INTENT(IN)            :: nElems                  ! Number of Elements
INTEGER,INTENT(IN)            :: nVal                    ! Number of nodal output variables
CHARACTER(LEN=32),INTENT(IN)  :: VarNames(nVal) ! Names of all variables that will be written out
REAL,INTENT(IN)               :: Coord(3,0:NPlot,0:NPlot,0:NPlot,nElems)      ! CoordsVector 
REAL,INTENT(IN)               :: Value(1:nVal,0:NPlot,0:NPlot,0:NPlot,nElems) ! Statevector 
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
CHARACTER(LEN=255)            :: myOutputFile
CHARACTER(LEN=255) :: cbuffer
CHARACTER(LEN=3)   :: dbuffer
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')" WRITE DATA TO VTK FILE..."   
  OPEN(1112,FILE=FileString,&
!        FORM='unformatted',&
!        RECORDTYPE='STREAM_LF',&
!        ACTION='write'                          ,&
        convert='BIG_ENDIAN'                    ,&
        ACCESS='stream'         ,&
        STATUS='replace')
  WRITE(1112) '# vtk DataFile Version 2.0'//char(10)
  WRITE(1112)'Debug Mesh '//char(10)
  WRITE(1112)'BINARY'//char(10)
  WRITE(1112)'DATASET UNSTRUCTURED_GRID'//char(10)
  WRITE(1112)''//char(10)
  WRITE(cbuffer,fmt='(A7,I0,A6)')'POINTS ',(NPlot+1)**3*nElems,' FLOAT'
  WRITE(1112)TRIM(cbuffer)//char(10)
  DO iElem=1,nElems
    DO k=0,NPlot
      DO j=0,NPlot
        DO i=0,NPlot
          WRITE(cbuffer,fmt='(ES22.15,1X,ES22.15,1X,ES22.15)')Coord(:,i,j,k,iElem)
          WRITE(1112) TRIM(cbuffer)//char(10)
        END DO 
      END DO 
    END DO 
  END DO 
  WRITE(1112)''//char(10)
  WRITE(cbuffer,fmt='(A,I0,1X,I0)')'CELLS ',nElems*NPlot**3,9*(nElems*NPlot**3)
  WRITE(1112)TRIM(cbuffer)//char(10)
  NodeIDElem = 0
  NPlot_p1_2=(NPlot+1)**2
  NPlot_p1_3=(NPlot+1)**3
  DO iElem=1,nElems
    DO k=1,NPlot
      DO j=1,NPlot
        DO i=1,NPlot
          WRITE(cbuffer,fmt='(I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0)')8 &
          ,NodeIDElem - 1 +i+  (j-1)*(NPlot+1)+(k-1)*NPlot_p1_2 &  !P1(CGNS=tecplot standard)
          ,NodeIDElem - 1 +i+1+(j-1)*(NPlot+1)+(k-1)*NPlot_p1_2 &  !P2
          ,NodeIDElem - 1 +i+1+ j   *(NPlot+1)+(k-1)*NPlot_p1_2 &  !P3
          ,NodeIDElem - 1 +i+   j   *(NPlot+1)+(k-1)*NPlot_p1_2 & !P4
          ,NodeIDElem - 1 +i+  (j-1)*(NPlot+1)+ k   *NPlot_p1_2 &  !P5
          ,NodeIDElem - 1 +i+1+(j-1)*(NPlot+1)+ k   *NPlot_p1_2 &  !P6
          ,NodeIDElem - 1 +i+1+ j   *(NPlot+1)+ k   *NPlot_p1_2 &  !P7
          ,NodeIDElem - 1 +i+   j   *(NPlot+1)+ k   *NPlot_p1_2   !P8  
          WRITE(1112)TRIM(cbuffer)//char(10)
        END DO
      END DO 
    END DO 
    NodeIDElem=NodeIDElem+NPlot_p1_3
  END DO 
  WRITE(1112)''//char(10)
  WRITE(cbuffer,fmt='(A,I0)')'CELL_TYPES ',nElems*NPlot**3
  WRITE(1112)TRIM(cbuffer)//char(10)
  DO iElem=1,nElems*NPlot**3
    WRITE(dbuffer,fmt='(I0,1X)') 12
    WRITE(1112) dbuffer
  END DO  
  WRITE(1112)''//char(10)
  CLOSE(1112)
END SUBROUTINE WriteDataToVTKBin

END MODULE MOD_OutPutVTK
#endif /*PARTICLES*/
