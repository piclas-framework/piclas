#include "boltzplatz.h"

MODULE  MOD_Particle_Debug
!===================================================================================================================================
!
!===================================================================================================================================
   IMPLICIT NONE                                                                                   !
   PRIVATE                                                                                         !
                                                                                                   !
!----------------------------------------------------------------------------------------------------------------------------------

   INTERFACE            WriteDebugMesh                                                             !
      MODULE PROCEDURE  WriteDebugMesh                                                             !
   END INTERFACE                                                                                   !

   PUBLIC WriteDebugMesh
!===================================================================================================================================

CONTAINS                                                                                           !
                                                                                                   !

SUBROUTINE WriteDebugMesh(n,iFIBG,kFIBG,lFIBG)                                                             !
!===================================================================================================================================
!===================================================================================================================================
   USE MOD_Particle_Vars
   USE MOD_Mesh_Vars,     ONLY : nElems
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
   IMPLICIT NONE                                                                                   !
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION                                                                        !
  INTEGER                            :: n                                                          !
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION                                                                       !
  CHARACTER(LEN=26)                  :: myFileName
  INTEGER                            :: cellType,i,j, iElem                                        !
  INTEGER                            :: nNodes,globSideID,iLocSide, iCorner                        !
  REAL                               :: xNode(1:3,1:8)                                             !
  INTEGER                            :: iFIBG, kFIBG, lFIBG, iNode
!===================================================================================================================================

  nNodes=nElems*8
  WRITE(myFileName,'(A8,I2.2,A4)')'PICDebug_',n,'.vtk'
  OPEN(1112,FILE=myFileName)
  WRITE(1112,'(A)')'# vtk DataFile Version 2.0'
  WRITE(1112,'(A)')'Debug Mesh '
  WRITE(1112,'(A)')'ASCII'
  WRITE(1112,'(A)')'DATASET UNSTRUCTURED_GRID'
  WRITE(1112,'(A)')''
  WRITE(1112,'(A,I0,A)')'POINTS ',nNodes,' FLOAT'
  DO iElem=1,nElems
    DO iCorner=1,8
      WRITE(1112,*)GEO%NodeCoords(1:3,GEO%ElemToNodeID(iCorner,iElem))
    END DO
  END DO
  WRITE(1112,*)''
  WRITE(1112,'(A,I0,1X,I0)')'CELLS ',nElems,nNodes+nElems
  nNodes=0
  DO iElem=1,nElems
    WRITE(1112,'(I0)',ADVANCE="NO")8
    DO j=1,8
      WRITE(1112,'(1X,I0)',ADVANCE="NO")nNodes
      nNodes=nNodes+1
    END DO
    WRITE(1112,*)''
  END DO
  WRITE(1112,*)''
  WRITE(1112,'(A,I0)')'CELL_TYPES ',nElems
  DO iElem=1,nElems
    WRITE(1112,'(1X,I0)',ADVANCE="NO")12
  END DO
  WRITE(1112,*)''
  WRITE(1112,*)''
  WRITE(1112,'(A,I0)')'CELL_DATA ',nElems
  WRITE(1112,'(A)')'SCALARS FIBGMCELL FLOAT'
  WRITE(1112,'(A)')'LOOKUP_TABLE default'
  DO iElem=1,nElems
    IF (ANY(GEO%FIBGM(iFIBG,kFIBG,lFIBG)%Element(:).EQ.iElem)) THEN
      WRITE(1112,'(ES21.15)')1.0
    ELSE
      WRITE(1112,'(ES21.15)')0.0
    END IF
  END DO
  WRITE(1112,*)''
  CLOSE(1112)
  RETURN
END SUBROUTINE WriteDebugMesh


END MODULE MOD_Particle_Debug
