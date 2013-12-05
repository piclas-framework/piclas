#include "boltzplatz.h"

MODULE  MOD_InitializeBackgroundField
!===================================================================================================================================
! 
!===================================================================================================================================
   IMPLICIT NONE                                                                                   !
   PRIVATE                                                                                         !
!----------------------------------------------------------------------------------------------------------------------------------
   PUBLIC :: InitializeBackgroundField, InitializeBackgroundBField  
!===================================================================================================================================

CONTAINS                                                                                           !
                                                                                                   !
SUBROUTINE InitializeBackgroundField
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Globals!,                ONLY : UNIT_errOut, UNIT_StdOut
USE MOD_ReadInTools
USE MOD_Mesh_Vars,              ONLY : nNodes
USE MOD_PICInterpolation_Vars,  ONLY : BGEfieldAtNode, eps_distance
USE MOD_Particle_Vars,          ONLY : GEO
USE MOD_Part_MPI_Vars,          ONLY : PMPIVAR
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                :: unit_in
INTEGER                :: os !openStatus
CHARACTER(255)         :: VTKfile
CHARACTER(LEN=255)     :: cdummy
INTEGER                :: npoints, iNode, ncells, icell
REAL, ALLOCATABLE      :: VTKNodes(:,:), VTK_Efield_Nodes(:,:)
LOGICAL, ALLOCATABLE   :: IsAssociated1(:), IsAssociated2(:)
REAL                   :: x(3), dist, tempdist
INTEGER                :: VTK2LocNodes(10), iVTKnode, CellX, CellY, CellZ, nFoundNodes, ppp, ElemID
INTEGER                :: iFoundNode, tempnode
!===================================================================================================================================

IF(MPIroot)THEN
  WRITE(UNIT_stdOut,'(132("~"))')
  WRITE(UNIT_stdOut,'(A)')'Reading VTK file for BG E-field...'
END IF

VTKfile = GETSTR('BGEField-VTK-File','blubb')
IF(TRIM(VTKfile).EQ.'blubb')THEN 
  CALL abort(__STAMP__,&
  'ERROR: no VTK-Filename for Background-Field E-Field defined!',999,999.)
END IF 

unit_in = 1123
OPEN(UNIT   = unit_in,              &
     FILE   = VTKfile,              &
     IOSTAT = os,                   &
     STATUS = 'OLD',                &
     ACTION = 'READ',               &
     ACCESS = 'SEQUENTIAL'          )

IF(os.NE.0) THEN  ! File Error
  CALL abort(__STAMP__,&
  'ERROR: cannot open VTK file: '//trim(VTKfile),999,999.)
END IF
  
!  First two lines contain the version number and a title
DO iNode=1,7
  READ(unit_in, '(A)') cdummy                             ! # vtk DataFile Version 2.0
END DO
!READ(unit_in, '(A)') cdummy                             ! cube
!READ(unit_in, '(A)') cdummy                             ! ASCII
!READ(unit_in, '(A)') cdummy                             ! DATASET UNSTRUCTURED_GRID
!READ(unit_in, '(A6,1X,I,1X,A5)') cdummy,npoints,cdummy  ! POINTS ???? float
READ(unit_in,*) cdummy,npoints,cdummy  ! POINTS ???? float
!print*, "check", npoints
!STOP
ALLOCATE(VTKNodes(1:3,npoints))
DO iNode = 0,INT(npoints/3)-1
  READ(unit_in,*) VTKNodes(:,3*iNode+1),VTKNodes(:,3*iNode+2),VTKNodes(:,3*iNode+3)
END DO
IF (MOD(npoints,3).EQ.1) THEN
  READ(unit_in,*) VTKNodes(:,3*iNode+1)
ELSE IF (MOD(npoints,3).EQ.2) THEN
  READ(unit_in,*) VTKNodes(:,3*iNode+1),VTKNodes(:,3*iNode+2)
END IF
READ(unit_in,*) cdummy,ncells,cdummy  ! CELLS ???? ????
DO icell = 1,ncells
  READ(unit_in,*) cdummy !skip cells
END DO
READ(unit_in, '(A)') cdummy  ! blank line
READ(unit_in, '(A)') cdummy  ! blank line
DO icell = 1,ncells
  READ(unit_in,*) cdummy !skip cells
END DO
DO iNode=1,3
  READ(unit_in, '(A)') cdummy  
END DO
ALLOCATE(VTK_Efield_Nodes(1:3,npoints))
DO iNode = 0,INT(npoints/3)-1
  !print*, iNode
  READ(unit_in,*) VTK_EField_Nodes(:,3*iNode+1),VTK_EField_Nodes(:,3*iNode+2),VTK_EField_Nodes(:,3*iNode+3)
END DO
IF (MOD(npoints,3).EQ.1) THEN
  READ(unit_in,*) VTK_EField_Nodes(:,3*iNode+1)
ELSE IF (MOD(npoints,3).EQ.2) THEN
  READ(unit_in,*) VTK_EField_Nodes(:,3*iNode+1),VTK_EField_Nodes(:,3*iNode+2)
END IF
!print*, cdummy
!STOP
  CLOSE(1123)

ALLOCATE(BGEfieldAtNode(1:3,nNodes),   &
         IsAssociated1(1:nNodes)   ,   &   ! prime array
         IsAssociated2(1:nNodes))          ! dummy array
IsAssociated1 = .FALSE.
BGEfieldAtNode = 0.

eps_distance = GETREAL('BGField-VTK-eps','1.E-10')

#ifndef MPI
IF (npoints.NE.nNodes) THEN
  CALL abort(__STAMP__,&
  'ERROR: wrong number of points in VTK-File',999,999.)
END IF
#endif /*MPI*/
DO iVTKnode = 1,npoints
  nFoundNodes = 0
  IsAssociated2 = .FALSE.
  x = VTKNodes(1:3,iVTKnode)
  CellX = INT((x(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
  CellY = INT((x(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
  CellZ = INT((x(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
#ifdef MPI
  IF ((GEO%FIBGMimax.GE.CellX).AND.(GEO%FIBGMimin.LE.CellX)) THEN  
  IF ((GEO%FIBGMkmax.GE.CellY).AND.(GEO%FIBGMkmin.LE.CellY)) THEN  
  IF ((GEO%FIBGMlmax.GE.CellZ).AND.(GEO%FIBGMlmin.LE.CellZ)) THEN  
#endif /* MPI */
  DO ppp = 1,GEO%FIBGM(CellX,CellY,CellZ)%nElem    
    ElemID = GEO%FIBGM(CellX,CellY,CellZ)%Element(ppp)
    DO iNode = 1,8 
      IF (.NOT.IsAssociated1(GEO%ElemToNodeID(iNode,ElemID))) THEN
        IF (ABS(x(1)-GEO%NodeCoords(1,GEO%ElemToNodeID(iNode,ElemID))).LE.eps_distance) THEN
          IF (ABS(x(2)-GEO%NodeCoords(2,GEO%ElemToNodeID(iNode,ElemID))).LE.eps_distance) THEN
            IF (ABS(x(3)-GEO%NodeCoords(3,GEO%ElemToNodeID(iNode,ElemID))).LE.eps_distance) THEN
              IF (.NOT.IsAssociated2(GEO%ElemToNodeID(iNode,ElemID))) THEN 
                IsAssociated2(GEO%ElemToNodeID(iNode,ElemID)) = .TRUE.
                nFoundNodes = nFoundNodes + 1
                IF (nFoundNodes.GT.10) THEN
                  CALL abort(__STAMP__,&
			'ERROR: Found to many nodes in VTK2LocNodes-Mapping. Decrease eps-region',999,999.)
                END IF
                VTK2LocNodes(nFoundNodes) = GEO%ElemToNodeID(iNode,ElemID)
              END IF
            END IF
          END IF
        END IF
      END IF ! IsAssociated1
    END DO
  END DO !ppp
  IF (nFoundNodes.EQ.1) THEN
    BGEfieldAtNode(:,VTK2LocNodes(1)) = VTK_EField_Nodes(:,iVTKnode)
    IsAssociated1(VTK2LocNodes(1)) = .TRUE.
  ELSE IF (nFoundNodes.GT.1) THEN
    dist = HUGE(dist) 
    DO iFoundNode = 1, nFoundNodes
      tempdist = SQRT((x(1)-GEO%NodeCoords(1,VTK2LocNodes(iFoundNode)))**2 + &
                 (x(2)-GEO%NodeCoords(2,VTK2LocNodes(iFoundNode)))**2 + &
                 (x(3)-GEO%NodeCoords(3,VTK2LocNodes(iFoundNode)))**2)
      IF(tempdist.LE.dist) THEN
         dist = tempdist 
         tempnode = iFoundNode
      END IF
    END DO
    BGEfieldAtNode(:,VTK2LocNodes(tempnode)) = VTK_EField_Nodes(:,iVTKnode)
    IsAssociated1(VTK2LocNodes(tempnode)) = .TRUE.
  END IF 
#ifdef MPI
  END IF
  END IF 
  END IF  
#endif /* MPI */
! Muss noch MPI-Fit gemacht werden.  
!#ifdef MPI
!  InsideMyBGM=.TRUE.
!  IF ((CellX.GT.GEO%FIBGMimax).OR.(CellX.LT.GEO%FIBGMimin) .OR. &
!      (CellY.GT.GEO%FIBGMkmax).OR.(CellY.LT.GEO%FIBGMkmin) .OR. &
!      (CellZ.GT.GEO%FIBGMlmax).OR.(CellZ.LT.GEO%FIBGMlmin)) THEN
!    InsideMyBGM=.FALSE.
!  END If
!  IF (InsideMyBGM) THEN
!  END IF
!#endif
END DO ! iNode
IF (.NOT. ALL(IsAssociated1)) THEN
  CALL abort(__STAMP__,&
	'ERROR: Not all nodes mapped for BGEfield!',999,999.)
END IF
END SUBROUTINE InitializeBackgroundField

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!


SUBROUTINE InitializeBackgroundBField
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Globals!,                ONLY : UNIT_errOut, UNIT_StdOut
USE MOD_ReadInTools
USE MOD_Mesh_Vars,              ONLY : nNodes
USE MOD_PICInterpolation_Vars,  ONLY : BGBfieldAtNode, eps_distance
USE MOD_Particle_Vars,          ONLY : GEO
USE MOD_Part_MPI_Vars,          ONLY : PMPIVAR
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                :: unit_in
INTEGER                :: os !openStatus
CHARACTER(255)         :: VTKfile
CHARACTER(LEN=255)     :: cdummy
INTEGER                :: npoints, iNode, ncells, icell
REAL, ALLOCATABLE      :: VTKNodes(:,:), VTK_Bfield_Nodes(:,:)
LOGICAL, ALLOCATABLE   :: IsAssociated1(:), IsAssociated2(:)
REAL                   :: x(3), dist, tempdist
INTEGER                :: VTK2LocNodes(10), iVTKnode, CellX, CellY, CellZ, nFoundNodes, ppp, ElemID
INTEGER                :: iFoundNode, tempnode
!===================================================================================================================================

IF(MPIroot)THEN
  WRITE(UNIT_stdOut,'(132("~"))')
  WRITE(UNIT_stdOut,'(A)')'Reading VTK file for BG B-field...'
END IF

VTKfile = GETSTR('BGBField-VTK-File','blubb')
IF(TRIM(VTKfile).EQ.'blubb')THEN 
  CALL abort(__STAMP__,&
  'ERROR: no VTK-Filename for Background-Field B-Field defined!',999,999.)
END IF 

unit_in = 1123
OPEN(UNIT   = unit_in,              &
     FILE   = VTKfile,              &
     IOSTAT = os,                   &
     STATUS = 'OLD',                &
     ACTION = 'READ',               &
     ACCESS = 'SEQUENTIAL'          )

IF(os.NE.0) THEN  ! File Error
  CALL abort(__STAMP__,&
  'ERROR: cannot open VTK file: '//trim(VTKfile),999,999.)
END IF
  
!  First two lines contain the version number and a title
DO iNode=1,7
  READ(unit_in, '(A)') cdummy                             ! # vtk DataFile Version 2.0
END DO
!READ(unit_in, '(A)') cdummy                             ! cube
!READ(unit_in, '(A)') cdummy                             ! ASCII
!READ(unit_in, '(A)') cdummy                             ! DATASET UNSTRUCTURED_GRID
!READ(unit_in, '(A6,1X,I,1X,A5)') cdummy,npoints,cdummy  ! POINTS ???? float
READ(unit_in,*) cdummy,npoints,cdummy  ! POINTS ???? float
!print*, "check", npoints
!STOP
ALLOCATE(VTKNodes(1:3,npoints))
DO iNode = 0,INT(npoints/3)-1
  READ(unit_in,*) VTKNodes(:,3*iNode+1),VTKNodes(:,3*iNode+2),VTKNodes(:,3*iNode+3)
END DO
IF (MOD(npoints,3).EQ.1) THEN
  READ(unit_in,*) VTKNodes(:,3*iNode+1)
ELSE IF (MOD(npoints,3).EQ.2) THEN
  READ(unit_in,*) VTKNodes(:,3*iNode+1),VTKNodes(:,3*iNode+2)
END IF
READ(unit_in,*) cdummy,ncells,cdummy  ! CELLS ???? ????
DO icell = 1,ncells
  READ(unit_in,*) cdummy !skip cells
END DO
READ(unit_in, '(A)') cdummy  ! blank line
READ(unit_in, '(A)') cdummy  ! blank line
DO icell = 1,ncells
  READ(unit_in,*) cdummy !skip cells
END DO
DO iNode=1,3
  READ(unit_in, '(A)') cdummy  
END DO
ALLOCATE(VTK_Bfield_Nodes(1:3,npoints))
DO iNode = 0,INT(npoints/3)-1
  !print*, iNode
  READ(unit_in,*) VTK_Bfield_Nodes(:,3*iNode+1),VTK_Bfield_Nodes(:,3*iNode+2),VTK_Bfield_Nodes(:,3*iNode+3)
END DO
IF (MOD(npoints,3).EQ.1) THEN
  READ(unit_in,*) VTK_Bfield_Nodes(:,3*iNode+1)
ELSE IF (MOD(npoints,3).EQ.2) THEN
  READ(unit_in,*) VTK_Bfield_Nodes(:,3*iNode+1),VTK_Bfield_Nodes(:,3*iNode+2)
END IF
!print*, cdummy
!STOP

  CLOSE(1123)

ALLOCATE(BGBfieldAtNode(1:3,nNodes),   &
         IsAssociated1(1:nNodes)   ,   &   ! prime array
         IsAssociated2(1:nNodes))          ! dummy array
IsAssociated1 = .FALSE.
BGBfieldAtNode = 0.

IF (eps_distance.EQ.0) THEN
  eps_distance = GETREAL('BGField-VTK-eps','1.E-10')
END IF

#ifndef MPI
IF (npoints.NE.nNodes) THEN
  CALL abort(__STAMP__,&
  'ERROR: wrong number of points in VTK-File',999,999.)
END IF
#endif /*MPI*/
DO iVTKnode = 1,npoints
  nFoundNodes = 0
  IsAssociated2 = .FALSE.
  x = VTKNodes(1:3,iVTKnode)
  CellX = INT((x(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
  CellY = INT((x(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
  CellZ = INT((x(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
#ifdef MPI
  IF ((GEO%FIBGMimax.GE.CellX).AND.(GEO%FIBGMimin.LE.CellX)) THEN  
  IF ((GEO%FIBGMkmax.GE.CellY).AND.(GEO%FIBGMkmin.LE.CellY)) THEN  
  IF ((GEO%FIBGMlmax.GE.CellZ).AND.(GEO%FIBGMlmin.LE.CellZ)) THEN  
#endif /* MPI */
  DO ppp = 1,GEO%FIBGM(CellX,CellY,CellZ)%nElem    
    ElemID = GEO%FIBGM(CellX,CellY,CellZ)%Element(ppp)
    DO iNode = 1,8 
      IF (.NOT.IsAssociated1(GEO%ElemToNodeID(iNode,ElemID))) THEN
        IF (ABS(x(1)-GEO%NodeCoords(1,GEO%ElemToNodeID(iNode,ElemID))).LE.eps_distance) THEN
          IF (ABS(x(2)-GEO%NodeCoords(2,GEO%ElemToNodeID(iNode,ElemID))).LE.eps_distance) THEN
            IF (ABS(x(3)-GEO%NodeCoords(3,GEO%ElemToNodeID(iNode,ElemID))).LE.eps_distance) THEN
              IF (.NOT.IsAssociated2(GEO%ElemToNodeID(iNode,ElemID))) THEN 
                IsAssociated2(GEO%ElemToNodeID(iNode,ElemID)) = .TRUE.
                nFoundNodes = nFoundNodes + 1
                IF (nFoundNodes.GT.10) THEN
                  CALL abort(__STAMP__,&
                  'ERROR: Found to many nodes in VTK2LocNodes-Mapping. Decrease eps-region',999,999.)
                END IF
                VTK2LocNodes(nFoundNodes) = GEO%ElemToNodeID(iNode,ElemID)
              END IF
            END IF
          END IF
        END IF
      END IF ! IsAssociated1
    END DO
  END DO !ppp
  IF (nFoundNodes.EQ.1) THEN
    BGBfieldAtNode(:,VTK2LocNodes(1)) = VTK_Bfield_Nodes(:,iVTKnode)
    IsAssociated1(VTK2LocNodes(1)) = .TRUE.
  ELSE IF (nFoundNodes.GT.1) THEN
    dist = HUGE(dist) 
    DO iFoundNode = 1, nFoundNodes
      tempdist = SQRT((x(1)-GEO%NodeCoords(1,VTK2LocNodes(iFoundNode)))**2 + &
                 (x(2)-GEO%NodeCoords(2,VTK2LocNodes(iFoundNode)))**2 + &
                 (x(3)-GEO%NodeCoords(3,VTK2LocNodes(iFoundNode)))**2)
      IF(tempdist.LE.dist) THEN
         dist = tempdist 
         tempnode = iFoundNode
      END IF
    END DO
    BGBfieldAtNode(:,VTK2LocNodes(tempnode)) = VTK_Bfield_Nodes(:,iVTKnode)
    IsAssociated1(VTK2LocNodes(tempnode)) = .TRUE.
  END IF 
#ifdef MPI
  END IF
  END IF 
  END IF  
#endif /* MPI */
! Muss noch MPI-Fit gemacht werden.  
!#ifdef MPI
!  InsideMyBGM=.TRUE.
!  IF ((CellX.GT.GEO%FIBGMimax).OR.(CellX.LT.GEO%FIBGMimin) .OR. &
!      (CellY.GT.GEO%FIBGMkmax).OR.(CellY.LT.GEO%FIBGMkmin) .OR. &
!      (CellZ.GT.GEO%FIBGMlmax).OR.(CellZ.LT.GEO%FIBGMlmin)) THEN
!    InsideMyBGM=.FALSE.
!  END If
!  IF (InsideMyBGM) THEN
!  END IF
!#endif
END DO ! iNode
IF (.NOT. ALL(IsAssociated1)) THEN
  CALL abort(__STAMP__,&
  'ERROR: Not all nodes mapped for BGBfield!',999,999.)
END IF
END SUBROUTINE InitializeBackgroundBField


END MODULE
