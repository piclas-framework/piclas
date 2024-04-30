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

!==================================================================================================================================
!> The PICLAS2VTK tool takes state files written during runtime by PICLAS in the .h5 format and converts them to .vtu files,
!> readable by ParaView. Supports parallel readin.
!> The state files can come from different calculations with different mesh files, equation systems, polynomial degrees and so on.
!> Two modes of usage: command line mode and parameter file mode.
!> In parameter file mode the usage is: piclas2vtk parameter.ini State1.h5 State2.h5 State3.h5 ...
!> In the parameter file the following can be specified:
!> - NVisu: Integer, polynomial degree of visualization basis
!> - NodeTypeVisu: String, node type of visualization basis
!> In command line mode, only the degree of the visualization basis can be directly specified, no parameter file is needed:
!> piclas2vtk --NVisu=INTEGER State1.h5 State2.h5 State3.h5 ...
!> All other options are set to their standard values.
!==================================================================================================================================
MODULE MOD_piclas2vtk_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE

REAL,ALLOCATABLE            :: NodeCoords_Connect(:,:)
INTEGER,ALLOCATABLE         :: ElemUniqueNodeID(:,:)
INTEGER                     :: nUniqueNodes
!----------------------------------------------------------------------------------------------------------------------------------
! Mapping of nodes and surface sides, required for connectivity of elements
!----------------------------------------------------------------------------------------------------------------------------------
TYPE tSurfaceConnect
  INTEGER                         :: nSurfaceNode                 !< Number of Nodes on Surface (reflective)
  INTEGER                         :: nSurfaceBCSides              !< Number of Sides on Surface (reflective)
  REAL, ALLOCATABLE               :: NodeCoords(:,:)
  INTEGER, ALLOCATABLE            :: SideSurfNodeMap(:,:)         !< Mapping from glob Side to SurfaceNodeNum (1:4, nSurfaceBCSides)
  INTEGER, ALLOCATABLE            :: SurfSideToSide(:)
END TYPE

TYPE (tSurfaceConnect)               :: SurfConnect

END MODULE MOD_piclas2vtk_Vars
