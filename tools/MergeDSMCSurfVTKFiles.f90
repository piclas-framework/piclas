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
PROGRAM mergeDSMCVTKFiles

IMPLICIT NONE

  CHARACTER(LEN=255)          :: tmp1, tmp2
  INTEGER                     :: unit_in, os !OpenStatus
  INTEGER                     :: numOfProcs, DSMCOut_num, iLoop, iTemp, iNode, iNode2, same_nodes, iNodeCell
  INTEGER                     :: IWNode, IWNode2, iElem, withMolecules, nSpecies, iSpec
  INTEGER                     :: SkipNodeNum, SkipCellNodeNum
  INTEGER, ALLOCATABLE        :: SkipNodes(:), SkipCellNodes(:,:)
  CHARACTER(LEN=26)           :: myFileName
  CHARACTER(255)              :: VTKfile
  CHARACTER(LEN=255)          :: cdummy
  INTEGER                     :: npoints_ges, ncells_ges
  INTEGER, ALLOCATABLE        :: npoints(:), ncells(:)
  CHARACTER(LEN=24), ALLOCATABLE      :: VTKNodes_ges(:,:)
  INTEGER, ALLOCATABLE        :: VTKCells_ges(:,:)
  LOGICAL                     :: NodeChanged, NodeSkipped
  TYPE tProc_Inf
    CHARACTER(LEN=24), ALLOCATABLE         :: VTKNodes(:,:)
    INTEGER, ALLOCATABLE      :: VTKCells(:,:)
    REAL, ALLOCATABLE           :: Heatflux(:), Force(:,:), Counter(:)                
  END TYPE
  TYPE(tProc_Inf), ALLOCATABLE     :: Proc_Inf(:)
  
  npoints_ges = 0
  IF (IARGC().LT.2) THEN
    WRITE(*,*) 'Too less input arguments!'
    WRITE(*,*) 'Input argument 1: number of used proc in simulation'
    WRITE(*,*) 'Input argument 2: number of dsmc output'
    STOP
  END IF
  
  CALL SYSTEM('rm DSMCSurfOut_Merged.vtk')
  WRITE(myFileName,'(A22)')'DSMCSurfOut_Merged.vtk'
  OPEN(1112,FILE=myFileName)
  CALL GETARG(1,tmp1)
  CALL GETARG(2,tmp2)  
  READ(tmp1,*) numOfProcs 
  READ(tmp2,*) DSMCOut_num 
  ALLOCATE(Proc_Inf(numOfProcs))
  ALLOCATE(npoints(numOfProcs))
  ALLOCATE(ncells(numOfProcs))
  
  DO iLoop = 1, numOfProcs
    WRITE(*,*)'Reading DSMC VTKFile ',iLoop
    WRITE(VTKfile,'(A17,I5.5,A1,I4.4,A4)')'DSMCSurfaceValues',iLoop - 1,'_',DSMCOut_num,'.vtk'
    unit_in = 1123
    OPEN(UNIT   = unit_in,              &
         FILE   = VTKfile,              &
         IOSTAT = os,                   &
         STATUS = 'OLD',                &
         ACTION = 'READ',               &
         ACCESS = 'SEQUENTIAL'          )
    IF (iLoop.EQ.1) THEN
      READ(unit_in,*) cdummy
      DO iTemp = 1 , 4 !Read Header Data
          READ(unit_in, '(A)') cdummy   
      END DO  
    ELSE
      DO iTemp = 1 , 5 !Read Header Data
          READ(unit_in, '(A)') cdummy   
      END DO             
    END IF
    READ(unit_in, *) cdummy, npoints(iLoop), cdummy
    ALLOCATE (Proc_Inf(iLoop)%VTKNodes(1:3, npoints(iLoop)))
    DO iNode = 1, npoints(iLoop)
      READ(unit_in,*) Proc_Inf(iLoop)%VTKNodes(1,iNode),Proc_Inf(iLoop)%VTKNodes(2,iNode),Proc_Inf(iLoop)%VTKNodes(3,iNode)
    END DO    
    READ(unit_in, '(A)') cdummy  
    READ(unit_in,*) cdummy,ncells(iLoop),cdummy
    ALLOCATE (Proc_Inf(iLoop)%VTKCells(1:4, ncells(iLoop)))
    ALLOCATE(Proc_Inf(iLoop)%Force(3,ncells(iLoop)),&
             Proc_Inf(iLoop)%Heatflux(ncells(iLoop)),&
             Proc_Inf(iLoop)%Counter(ncells(iLoop)))
    DO iNode = 1, ncells(iLoop)
      READ(unit_in,*), cdummy, Proc_Inf(iLoop)%VTKCells(:,iNode)
    END DO
    DO iTemp = 1 , 6 !read dummy stuff
        READ(unit_in, '(A)') cdummy   
    END DO     
    READ(unit_in, '(A)') cdummy
    READ(unit_in, '(A)') cdummy
    DO iNode = 1, ncells(iLoop)
     READ(unit_in,*) Proc_Inf(iLoop)%Heatflux(iNode)
    END DO
    READ(unit_in, '(A)') cdummy
    READ(unit_in, '(A)') cdummy
    DO iNode = 1, ncells(iLoop)
     READ(unit_in,*) Proc_Inf(iLoop)%Force(1,iNode)
    END DO
    READ(unit_in, '(A)') cdummy
    READ(unit_in, '(A)') cdummy
    DO iNode = 1, ncells(iLoop)
     READ(unit_in,*) Proc_Inf(iLoop)%Force(2,iNode)
    END DO
    READ(unit_in, '(A)') cdummy
    READ(unit_in, '(A)') cdummy
    DO iNode = 1, ncells(iLoop)
     READ(unit_in,*) Proc_Inf(iLoop)%Force(3,iNode)
    END DO
    READ(unit_in, '(A)') cdummy
    READ(unit_in, '(A)') cdummy
    DO iNode = 1, ncells(iLoop)
     READ(unit_in,*) Proc_Inf(iLoop)%Counter(iNode)
    END DO     
  END DO
  WRITE(*,*)'DONE!'

  ALLOCATE(VTKNodes_ges(1:3, SUM(npoints)))
  ALLOCATE(VTKCells_ges(1:4, SUM(ncells)))
  WRITE(*,*) 'Start Compare VTKFiles'
  VTKNodes_ges(1:3,1:npoints(1)) =  Proc_Inf(1)%VTKNodes(1:3, 1:npoints(1))
  VTKCells_ges(1:4,1:ncells(1)) = Proc_Inf(1)%VTKCells(1:4, 1:ncells(1))
  npoints_ges = npoints(1)
  ncells_ges = ncells(1)



  DO iLoop = 2, numOfProcs
    WRITE(*,*) 'Compare File ', iLoop
    SkipCellNodeNum = 0
    ALLOCATE (SkipCellNodes(npoints(iLoop),2))
    DO iNode = 0 , npoints(iLoop) - 1 !Search equal points
      DO iNode2 = 0, npoints_ges - 1 
        IF (TRIM(VTKNodes_ges(1,iNode2+1)).EQ.TRIM(Proc_Inf(iLoop)%VTKNodes(1, iNode + 1))) THEN
          IF (TRIM(VTKNodes_ges(2,iNode2+1)).EQ.TRIM(Proc_Inf(iLoop)%VTKNodes(2, iNode + 1))) THEN
            IF (TRIM(VTKNodes_ges(3,iNode2+1)).EQ.TRIM(Proc_Inf(iLoop)%VTKNodes(3, iNode + 1))) THEN
               SkipCellNodeNum = SkipCellNodeNum + 1
               SkipCellNodes(SkipCellNodeNum,1) = iNode2
               SkipCellNodes(SkipCellNodeNum,2) = iNode 
            END IF
          END IF
        END IF  
      END DO
    END DO

    DO iNodeCell = 1, ncells(iLoop) !Sorting the nodes
      DO IWNode = 1, 4
        NodeChanged = .FALSE.
        SkipNodeNum = 0
        DO IWNode2 = 1 , SkipCellNodeNum         
          IF (((Proc_Inf(iLoop)%VTKCells(IWNode, iNodeCell)).EQ.SkipCellNodes(IWNode2,2)).AND.(.NOT.NodeChanged)) THEN
              Proc_Inf(iLoop)%VTKCells(IWNode, iNodeCell) = SkipCellNodes(IWNode2,1)
              NodeChanged = .TRUE.
          ELSE IF ((Proc_Inf(iLoop)%VTKCells(IWNode, iNodeCell)).GT.SkipCellNodes(IWNode2,2)) THEN
              SkipNodeNum = SkipNodeNum + 1
          END IF          
        END DO
        IF(.NOT.NodeChanged) THEN
          Proc_Inf(iLoop)%VTKCells(IWNode, iNodeCell) = Proc_Inf(iLoop)%VTKCells(IWNode, iNodeCell) + npoints_ges -SkipNodeNum
        END IF
      END DO
    END DO
        

    SkipNodeNum = 0
    DO iNode = 1, npoints(iLoop) !Sorting Nodes in total array
      NodeSkipped = .FALSE.
      DO iNode2 = 1, SkipCellNodeNum
        IF ((iNode - 1).EQ.SkipCellNodes(iNode2,2)) NodeSkipped = .TRUE.
      END DO
      IF(.NOT.NodeSkipped) THEN
        SkipNodeNum = SkipNodeNum + 1
        VTKNodes_ges(1:3,npoints_ges + SkipNodeNum) = Proc_Inf(iLoop)%VTKNodes(1:3, iNode)
      END IF
    END DO
    npoints_ges = npoints_ges + SkipNodeNum


    DO iNode = 1, ncells(iLoop) !Sorting cell in total arrays 
      VTKCells_ges(1:4,ncells_ges + iNode) = Proc_Inf(iLoop)%VTKCells(1:4, iNode)
    END DO
    ncells_ges = ncells_ges + ncells(iLoop)


    DEALLOCATE (SkipCellNodes)
  END DO
  
  WRITE(*,*) 'Comparing Finished'
  WRITE(*,*) 'Writing new Output File'

  WRITE(1112,'(A)')'# vtk DataFile Version 2.0'
  WRITE(1112,'(A)')'Debug Mesh '
  WRITE(1112,'(A)')'ASCII'
  WRITE(1112,'(A)')'DATASET UNSTRUCTURED_GRID'
  WRITE(1112,'(A)')''
  WRITE(1112,'(A,I0,A)')'POINTS ',npoints_ges,' FLOAT'
  DO iNode=1, npoints_ges
    WRITE(1112,*) VTKNodes_ges(1:3, iNode)
  END DO
  WRITE(1112,*)''
  WRITE(1112,'(A,I0,1X,I0)')'CELLS ',ncells_ges,5*ncells_ges
  DO iElem=1, ncells_ges
    WRITE(1112,'(I0)',ADVANCE="NO")4
    DO iNode=1, 4
     WRITE(1112,'(1X,I0)',ADVANCE="NO") VTKCells_ges(iNode,iElem)
    END DO
    WRITE(1112,*)''
  END DO
  WRITE(1112,*)''
  WRITE(1112,'(A,I0)')'CELL_TYPES ',ncells_ges
  DO iElem=1,ncells_ges
    WRITE(1112,'(1X,I0)',ADVANCE="NO")9
  END DO  
  WRITE(1112,*)''
  WRITE(1112,*)''
  WRITE(1112,'(A,I0)')'CELL_DATA ',ncells_ges
  WRITE(1112,'(A)')'SCALARS Heatflux FLOAT'
  WRITE(1112,'(A)')'LOOKUP_TABLE default'
  DO iLoop = 1, numOfProcs
    DO iNode = 1, ncells(iLoop)
      WRITE(1112,'(E22.15)') Proc_Inf(iLoop)%Heatflux(iNode)
    END DO
  END DO
  WRITE(1112,'(A)')'SCALARS ForceX FLOAT'
  WRITE(1112,'(A)')'LOOKUP_TABLE default'
  DO iLoop = 1, numOfProcs
    DO iNode = 1, ncells(iLoop)
      WRITE(1112,'(E22.15)') Proc_Inf(iLoop)%Force(1,iNode)
    END DO
  END DO
  WRITE(1112,'(A)')'SCALARS ForceY FLOAT'
  WRITE(1112,'(A)')'LOOKUP_TABLE default'
  DO iLoop = 1, numOfProcs
    DO iNode = 1, ncells(iLoop)
      WRITE(1112,'(E22.15)') Proc_Inf(iLoop)%Force(2,iNode)
    END DO
  END DO
  WRITE(1112,'(A)')'SCALARS ForceZ FLOAT'
  WRITE(1112,'(A)')'LOOKUP_TABLE default'
  DO iLoop = 1, numOfProcs
    DO iNode = 1, ncells(iLoop)
      WRITE(1112,'(E22.15)') Proc_Inf(iLoop)%Force(3,iNode)
    END DO
  END DO
  WRITE(1112,'(A)')'SCALARS Counter FLOAT'
  WRITE(1112,'(A)')'LOOKUP_TABLE default'
  DO iLoop = 1, numOfProcs
    DO iNode = 1, ncells(iLoop)
      WRITE(1112,'(E22.15)') Proc_Inf(iLoop)%Counter(iNode)
    END DO
  END DO 
  CLOSE(1112)

  WRITE(*,*)'DONE! OUTPUT: DSMCOut_Merged.vtk'

END PROGRAM mergeDSMCVTKFiles
