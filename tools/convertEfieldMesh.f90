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
PROGRAM convertEFieldMesh

IMPLICIT NONE

  CHARACTER(LEN=255)          :: tmp1
  CHARACTER(LEN=255)          :: tmp2
  INTEGER                     :: unit_in, os !OpenStatus
  INTEGER                     :: numOfProcs, iLoop, iTemp, iNode, iNode2, same_nodes, iNodeCell, counter
  INTEGER                     :: IWNode, IWNode2, iElem, nSpecies, iSpec, iCell
  CHARACTER(LEN=255)          :: myFileName
  CHARACTER(255)              :: VTKfileEField
  CHARACTER(255)              :: VTKfileNew
  CHARACTER(LEN=255)          :: cdummy
  INTEGER                     :: npoints_old, npoints_new, ncells_old, ncells_new 
  REAL, ALLOCATABLE           :: VTKNodes_old(:,:), VTK_Efield_Nodes(:,:),VTKNodes_new(:,:), New_Efield_Nodes(:,:)
  INTEGER, ALLOCATABLE        :: VTKCells_old(:,:),VTKCells_new(:,:)
  LOGICAL                     :: InElementCheck
  REAL                        :: dist_node(8), sum_dist

  counter = 0
  IF (IARGC().LT.2) THEN
    WRITE(*,*) 'Too less input arguments!'
    WRITE(*,*) 'Input argument 1: original EField VTK File'
    WRITE(*,*) 'Input argument 2: new VTK EField Mesh'
    STOP
  END IF

  CALL GETARG(1,tmp1)
  CALL GETARG(2,tmp2)  
  READ(tmp1,*) VTKfileEField 
  READ(tmp2,*) VTKfileNew 
    
  unit_in = 1112
  WRITE(*,*)"Reading File 1"
  OPEN(UNIT   = unit_in,              &
       FILE   = TRIM(VTKfileEField),  &
       IOSTAT = os,                   &
       STATUS = 'OLD',                &
       ACTION = 'READ',               &
       ACCESS = 'SEQUENTIAL'          )

  DO iNode=1,7
    READ(unit_in, '(A)') cdummy                             ! # vtk DataFile Version 2.0
  END DO
  READ(unit_in,*) cdummy,npoints_old,cdummy  ! POINTS ???? float
  !print*, "check", npoints
  !STOP
  ALLOCATE(VTKNodes_old(1:3,npoints_old))
  DO iNode = 0,INT(npoints_old/3)-1
    READ(unit_in,*) VTKNodes_old(:,3*iNode+1),VTKNodes_old(:,3*iNode+2),VTKNodes_old(:,3*iNode+3)
  END DO
  IF (MOD(npoints_old,3).EQ.1) THEN
    READ(unit_in,*) VTKNodes_old(:,3*iNode+1)
  ELSE IF (MOD(npoints_old,3).EQ.2) THEN
    READ(unit_in,*) VTKNodes_old(:,3*iNode+1),VTKNodes_old(:,3*iNode+2)
  END IF
  READ(unit_in,*) cdummy,ncells_old,cdummy  ! CELLS ???? ????
  ALLOCATE (VTKCells_old(1:8, ncells_old))
  DO icell = 1,ncells_old
    READ(unit_in,*), cdummy, VTKCells_old(:,icell)
  END DO
  VTKCells_old(:,:) = VTKCells_old(:,:) + 1
  READ(unit_in, '(A)') cdummy  ! blank line
  READ(unit_in, '(A)') cdummy  ! blank line
  DO icell = 1,ncells_old
    READ(unit_in,*) cdummy !skip cells
  END DO
  DO iNode=1,3
    READ(unit_in, '(A)') cdummy  
  END DO
  ALLOCATE(VTK_Efield_Nodes(1:3,npoints_old))
  DO iNode = 0,INT(npoints_old/3)-1
    !print*, iNode
    READ(unit_in,*) VTK_EField_Nodes(:,3*iNode+1),VTK_EField_Nodes(:,3*iNode+2),VTK_EField_Nodes(:,3*iNode+3)
  END DO
  IF (MOD(npoints_old,3).EQ.1) THEN
    READ(unit_in,*) VTK_EField_Nodes(:,3*iNode+1)
  ELSE IF (MOD(npoints_old,3).EQ.2) THEN
    READ(unit_in,*) VTK_EField_Nodes(:,3*iNode+1),VTK_EField_Nodes(:,3*iNode+2)
  END IF
  CLOSE(1112)
  WRITE(*,*)"DONE!"
  
  WRITE(*,*)"Reading File 2"
  OPEN(UNIT   = unit_in,              &
       FILE   = TRIM(VTKfileNew),  &
       IOSTAT = os,                   &
       STATUS = 'OLD',                &
       ACTION = 'READ',               &
       ACCESS = 'SEQUENTIAL'          )

  DO iTemp = 1 , 5 !Read Header Data
      READ(unit_in, '(A)') cdummy   
  END DO             
  READ(unit_in, *) cdummy, npoints_new, cdummy
  ALLOCATE (VTKNodes_new(1:3, npoints_new))
  ALLOCATE(New_Efield_Nodes(1:3,npoints_new))
  DO iNode = 1, npoints_new
    READ(unit_in,*) VTKNodes_new(1,iNode),VTKNodes_new(2,iNode),VTKNodes_new(3,iNode)
  END DO   
  READ(unit_in, '(A)') cdummy  
  READ(unit_in,*) cdummy,ncells_new,cdummy
  ALLOCATE (VTKCells_new(1:8, ncells_new))
  DO iNode = 1, ncells_new
    READ(unit_in,*), cdummy, VTKCells_new(:,iNode)    
  END DO
  CLOSE(1112)
  WRITE(*,*)"DONE!"

  WRITE(*,*)"Searching Field of MeshPoints!"
  New_Efield_Nodes = 0.0  
  DO inode=1, npoints_new
    DO icell=1, ncells_old
      CALL ParticleInsideQuad3D(VTKNodes_old(:,VTKCells_old(1,icell)), &
            VTKNodes_old(:,VTKCells_old(4,icell)), &
            VTKNodes_old(:,VTKCells_old(3,icell)), &
            VTKNodes_old(:,VTKCells_old(2,icell)), &
            VTKNodes_new(:,inode), &
            InElementCheck)
      IF(InElementCheck) THEN
        CALL ParticleInsideQuad3D(VTKNodes_old(:,VTKCells_old(3,icell)), &
              VTKNodes_old(:,VTKCells_old(7,icell)), &
              VTKNodes_old(:,VTKCells_old(6,icell)), &
              VTKNodes_old(:,VTKCells_old(2,icell)), &
              VTKNodes_new(:,inode), &
              InElementCheck)
      END IF
      IF(InElementCheck) THEN
        CALL ParticleInsideQuad3D(VTKNodes_old(:,VTKCells_old(6,icell)), &
              VTKNodes_old(:,VTKCells_old(5,icell)), &
              VTKNodes_old(:,VTKCells_old(1,icell)), &
              VTKNodes_old(:,VTKCells_old(2,icell)), &
              VTKNodes_new(:,inode), &
              InElementCheck)
      END IF
      IF(InElementCheck) THEN
        CALL ParticleInsideQuad3D(VTKNodes_old(:,VTKCells_old(5,icell)), &
              VTKNodes_old(:,VTKCells_old(8,icell)), &
              VTKNodes_old(:,VTKCells_old(4,icell)), &
              VTKNodes_old(:,VTKCells_old(1,icell)), &
              VTKNodes_new(:,inode), &
              InElementCheck)
      END IF
      IF(InElementCheck) THEN
        CALL ParticleInsideQuad3D(VTKNodes_old(:,VTKCells_old(8,icell)), &
              VTKNodes_old(:,VTKCells_old(7,icell)), &
              VTKNodes_old(:,VTKCells_old(3,icell)), &
              VTKNodes_old(:,VTKCells_old(4,icell)), &
              VTKNodes_new(:,inode), &
              InElementCheck)
      END IF
      IF(InElementCheck) THEN
        CALL ParticleInsideQuad3D(VTKNodes_old(:,VTKCells_old(5,icell)), &
              VTKNodes_old(:,VTKCells_old(6,icell)), &
              VTKNodes_old(:,VTKCells_old(7,icell)), &
              VTKNodes_old(:,VTKCells_old(8,icell)), &
              VTKNodes_new(:,inode), &
              InElementCheck)
      END IF
      IF(InElementCheck) THEN
        DO iNode2 = 1,8
          dist_node(iNode2) = (VTKNodes_old(1,VTKCells_old(iNode2,icell))-VTKNodes_new(1,inode))**2 + &
                            (VTKNodes_old(2,VTKCells_old(iNode2,icell))-VTKNodes_new(2,inode))**2 + &
                            (VTKNodes_old(3,VTKCells_old(iNode2,icell))-VTKNodes_new(3,inode))**2
        END DO
        dist_node = SQRT(dist_node)
        sum_dist  = SUM(dist_node)
        dist_node = dist_node/sum_dist
        DO iNode2 = 1,8
          New_Efield_Nodes(1,inode) = New_Efield_Nodes(1,inode) + dist_node(iNode2) * VTK_EField_Nodes(1,VTKCells_old(iNode2,icell))
          New_Efield_Nodes(2,inode) = New_Efield_Nodes(2,inode) + dist_node(iNode2) * VTK_EField_Nodes(2,VTKCells_old(iNode2,icell))
          New_Efield_Nodes(3,inode) = New_Efield_Nodes(3,inode) + dist_node(iNode2) * VTK_EField_Nodes(3,VTKCells_old(iNode2,icell))
        END DO
        counter = counter + 1
        EXIT
      END IF
    END DO
  END DO  
  IF (counter.NE.npoints_new) THEN
    WRITE(*,*) "Not all points found!"
    STOP
  END IF
  WRITE(*,*)"DONE!"

  WRITE(*,*)"Write New Output File."

  WRITE(myFileName,'(A)')'CutEField.vtk'
  OPEN(1112,FILE=myFileName,STATUS='replace')
  WRITE(1112,'(A)')'# vtk DataFile Version 2.0 Species: ' 
  WRITE(1112,'(A)')'Debug Mesh '
  WRITE(1112,'(A)')'ASCII'
  WRITE(1112,'(A)')'DATASET UNSTRUCTURED_GRID'
  WRITE(1112,'(A)')''
  WRITE(1112,'(A)')''
  WRITE(1112,'(A)')''
  WRITE(1112,'(A,I0,A)')'POINTS ',npoints_new,' FLOAT'
  DO iNode = 0,INT(npoints_new/3)-1
    WRITE(1112,*) VTKNodes_new(:,3*iNode+1),VTKNodes_new(:,3*iNode+2),VTKNodes_new(:,3*iNode+3)
  END DO
  IF (MOD(npoints_new,3).EQ.1) THEN
    WRITE(1112,*) VTKNodes_new(:,3*iNode+1)
  ELSE IF (MOD(npoints_new,3).EQ.2) THEN
    WRITE(1112,*) VTKNodes_new(:,3*iNode+1),VTKNodes_new(:,3*iNode+2)
  END IF

  WRITE(1112,'(A,I0,1X,I0)')'CELLS ',ncells_new,9*ncells_new
  DO icell=1, ncells_new
    WRITE(1112,'(I0)',ADVANCE="NO")8
    WRITE(1112,'(1X,I0)',ADVANCE="NO") VTKCells_new(:,icell)
    WRITE(1112,*)''
  END DO

  WRITE(1112,*)''
  WRITE(1112,'(A,I0)')'CELL_TYPES ',ncells_new
  DO icell=1,ncells_new
    WRITE(1112,'(I0)')12
  END DO  
  WRITE(1112,*)''
  WRITE(1112,'(A,I0)')'POINT_DATA ',npoints_new
  WRITE(1112,'(A)')'VECTORS E float'
  
  DO iNode = 0,INT(npoints_new/3)-1
    !print*, iNode
     WRITE(1112,*) New_Efield_Nodes(:,3*iNode+1),New_Efield_Nodes(:,3*iNode+2),New_Efield_Nodes(:,3*iNode+3)
  END DO
  IF (MOD(npoints_new,3).EQ.1) THEN
    WRITE(1112,*) New_Efield_Nodes(:,3*iNode+1)
  ELSE IF (MOD(npoints_new,3).EQ.2) THEN
    WRITE(1112,*) New_Efield_Nodes(:,3*iNode+1),New_Efield_Nodes(:,3*iNode+2)
  END IF

  CLOSE(1112)
  WRITE(*,*)"DONE!"
 END PROGRAM convertEFieldMesh

SUBROUTINE ParticleInsideQuad3D(Node1,Node2,Node3,Node4,Point,InElementCheck)                    
!DEC$ ATTRIBUTES FORCEINLINE :: ParticleInsideQuad3D
!--------------------------------------------------------------------------------------------------!
  IMPLICIT NONE                                                                                    !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
  INTEGER                          :: i, Element  ,n                                                 !
  REAL                             :: det(2)                                                       !
! Local variable declaration                                                                       !
  INTEGER                          :: NodeNum                                          !
  REAL                             :: A(1:3,1:4), cross(3)
  REAL, INTENT(IN)                 :: Node1(3),Node2(3),Node3(3),Node4(3),Point(3)
  LOGICAL                          :: ConcaveElemSide
!--------------------------------------------------------------------------------------------------!
  LOGICAL, INTENT(OUT)             :: InElementCheck                                         !
!--------------------------------------------------------------------------------------------------!

!WRITE(*,*) 'Element 51633'
!DO NodeNum = 1,4
!DO iLocSide = 1,6
!WRITE(*,*) GEO%NodeCoords(1,GEO%ElemSideNodeID(NodeNum,iLocSide,51633)), &
!GEO%NodeCoords(2,GEO%ElemSideNodeID(NodeNum,iLocSide,51633)), &
!GEO%NodeCoords(3,GEO%ElemSideNodeID(NodeNum,iLocSide,51633))
!END DO
!END DO
!WRITE(*,*) 'Element 55300'
!DO NodeNum = 1,4
!DO iLocSide = 1,6
!WRITE(*,*) GEO%NodeCoords(1,GEO%ElemSideNodeID(NodeNum,iLocSide,55300)), &
!GEO%NodeCoords(2,GEO%ElemSideNodeID(NodeNum,iLocSide,55300)), &
!GEO%NodeCoords(3,GEO%ElemSideNodeID(NodeNum,iLocSide,55300))
!END DO
!END DO
!STOP

    
   !CHECK WHETHER ElemSide is concave
   A(:,1) = Node1(1:3) - Node4(1:3)
   A(:,2) = Node2(1:3) - Node4(1:3)
   A(:,3) = Node3(1:3) - Node4(1:3)

   !--- compute cross product for vector 1 and 3
   cross(1) = A(2,1) * A(3,3) - A(3,1) * A(2,3)
   cross(2) = A(3,1) * A(1,3) - A(1,1) * A(3,3)
   cross(3) = A(1,1) * A(2,3) - A(2,1) * A(1,3)
   !--- negative determinant of triangle 1 (points 1,3,2):
   det(1) = cross(1) * A(1,2) + &
                     cross(2) * A(2,2) + &
                     cross(3) * A(3,2)
   det(1) = -det(1)
   IF (det(1).lt.0) THEN
     ConcaveElemSide = .TRUE.
   ELSE
     ConcaveElemSide = .FALSE.
   END IF

   ! CHECK POINT
   InElementCheck = .TRUE.

   !--- A = vector from particle to node coords 
   IF(ConcaveElemSide) THEN
     A(:,1) = Node2(1:3) - Point(1:3)
     A(:,2) = Node3(1:3) - Point(1:3)
     A(:,3) = Node4(1:3) - Point(1:3)
     A(:,4) = Node1(1:3) - Point(1:3)
   ELSE
     A(:,1) = Node1(1:3) - Point(1:3)
     A(:,2) = Node2(1:3) - Point(1:3)
     A(:,3) = Node3(1:3) - Point(1:3)
     A(:,4) = Node4(1:3) - Point(1:3)
   END IF


   !--- compute cross product for vector 1 and 3
   cross(1) = A(2,1) * A(3,3) - A(3,1) * A(2,3)
   cross(2) = A(3,1) * A(1,3) - A(1,1) * A(3,3)
   cross(3) = A(1,1) * A(2,3) - A(2,1) * A(1,3)

   !--- negative determinant of triangle 1 (points 1,3,2):
   det(1) = cross(1) * A(1,2) + &
                     cross(2) * A(2,2) + &
                     cross(3) * A(3,2)
   det(1) = -det(1)
   !--- determinant of triangle 2 (points 1,3,4):
   det(2) = cross(1) * A(1,4) + &
                     cross(2) * A(2,4) + &
                     cross(3) * A(3,4)
   IF (det(1).lt.0) THEN
     InElementCheck = .FALSE.
   END IF
   IF (det(2).lt.0) THEN  
     InElementCheck = .FALSE.
   END IF


 RETURN
END SUBROUTINE ParticleInsideQuad3D


