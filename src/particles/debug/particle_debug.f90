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

MODULE  MOD_Particle_Debug
!===================================================================================================================================
! Module containing routines for debug
!===================================================================================================================================
   IMPLICIT NONE                                                                                   !
   PRIVATE                                                                                         !

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
   USE MOD_Mesh_Vars,          ONLY:nElems,XCL_NGeo,NGeo
   USE MOD_Particle_Mesh_Vars, ONLY:GEO
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
   IMPLICIT NONE                                                                                   !
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION                                                                        !
  INTEGER                            :: n                                                          !
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION                                                                       !
  CHARACTER(LEN=26)                  :: myFileName
  INTEGER                            :: j, iElem, nNodes, iFIBG, kFIBG, lFIBG
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
    WRITE(1112,*)XCL_NGeo(1:3,0   ,0   ,0   ,iElem)
    WRITE(1112,*)XCL_NGeo(1:3,NGeo,0   ,0   ,iElem)
    WRITE(1112,*)XCL_NGeo(1:3,0   ,NGeo,0   ,iElem)
    WRITE(1112,*)XCL_NGeo(1:3,NGeo,NGeo,0   ,iElem)
    WRITE(1112,*)XCL_NGeo(1:3,0   ,0   ,NGeo,iElem)
    WRITE(1112,*)XCL_NGeo(1:3,NGeo,0   ,NGeo,iElem)
    WRITE(1112,*)XCL_NGeo(1:3,0   ,NGeo,NGeo,iElem)
    WRITE(1112,*)XCL_NGeo(1:3,NGeo,NGeo,NGeo,iElem)
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
      WRITE(1112,'(ES21.14)')1.0
    ELSE
      WRITE(1112,'(ES21.14)')0.0
    END IF
  END DO
  WRITE(1112,*)''
  CLOSE(1112)
  RETURN
END SUBROUTINE WriteDebugMesh


END MODULE MOD_Particle_Debug
