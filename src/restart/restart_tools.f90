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

MODULE MOD_Restart_Tools
!===================================================================================================================================
! Module containing tools/procedures for handling PICLas restarts
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

PUBLIC :: ReadNodeSourceExtFromHDF5
!===================================================================================================================================

CONTAINS


SUBROUTINE ReadNodeSourceExtFromHDF5() 
!----------------------------------------------------------------------------------------------------------------------------------!
! Read NodeSourceExt from h5 file, which is stored as DG solution type field. Map this solution to equidistant-node polynomial and
! then to the global nodes
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_PreProc
USE MOD_Restart_Vars           ,ONLY: N_Restart
USE MOD_Mesh_Vars              ,ONLY: Vdm_N_EQ,OffsetElem
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_PICDepo_Vars           ,ONLY: NodeSourceExt,CellLocNodes_Volumes
USE MOD_Interpolation_Vars     ,ONLY: NodeTypeVISU,NodeType
USE MOD_Interpolation          ,ONLY: GetVandermonde
USE MOD_Dielectric_Vars        ,ONLY: DoDielectric
USE MOD_ChangeBasis            ,ONLY: ChangeBasis3D
USE MOD_HDF5_input             ,ONLY: ReadArray,GetDataSize
USE MOD_HDF5_Input             ,ONLY: File_ID,DatasetExists
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
! Space-separated list of input and output types. Use: (int|real|logical|...)_(in|out|inout)_dim(n)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE                   :: U_local(:,:,:,:,:)
LOGICAL                            :: DG_SourceExtExists
REAL                               :: NodeSourceExtEqui(1,0:1,0:1,0:1)
INTEGER(KIND=IK)                   :: OffsetElemTmp,PP_nElemsTmp,N_RestartTmp
INTEGER                            :: iElem
!===================================================================================================================================
IF(.NOT.DoDielectric) RETURN

! Temp. vars for integer KIND=8 possibility
OffsetElemTmp = INT(OffsetElem,IK)
PP_nElemsTmp  = INT(PP_nElems,IK)
N_RestartTmp  = INT(N_Restart,IK)

CALL DatasetExists(File_ID,'DG_SourceExt',DG_SourceExtExists)

IF(DG_SourceExtExists)THEN

! We need to interpolate the solution to the new computational grid
ALLOCATE(U_local(1,0:N_Restart,0:N_Restart,0:N_Restart,PP_nElems))
CALL ReadArray('DG_SourceExt',5,(/1_IK,N_RestartTmp+1_IK,N_RestartTmp+1_IK,N_RestartTmp+1_IK,PP_nElemsTmp/),&
OffsetElemTmp,5,RealArray=U_local)

! Allocate and determine Vandermonde mapping from NodeType to equidistant (visu) node set
ALLOCATE(Vdm_N_EQ(0:1,N_Restart))
CALL GetVandermonde(N_Restart, NodeType, 1, NodeTypeVISU, Vdm_N_EQ, modal=.FALSE.)

DO iElem =1, PP_nElems
  ! Map G/GL (current node type) to equidistant distribution
  CALL ChangeBasis3D(1, N_Restart, 1, Vdm_N_EQ, U_local(:,:,:,:,iElem),NodeSourceExtEqui(:,:,:,:))

  ASSOCIATE( NodeID => GEO%ElemToNodeID(:,iElem) )
    ! Copy values from equidistant distribution to Nodees
    NodeSourceExt(NodeID(1)) = NodeSourceExtEqui(1,0,0,0) * CellLocNodes_Volumes(NodeID(1))
    NodeSourceExt(NodeID(2)) = NodeSourceExtEqui(1,1,0,0) * CellLocNodes_Volumes(NodeID(2))
    NodeSourceExt(NodeID(3)) = NodeSourceExtEqui(1,1,1,0) * CellLocNodes_Volumes(NodeID(3))
    NodeSourceExt(NodeID(4)) = NodeSourceExtEqui(1,0,1,0) * CellLocNodes_Volumes(NodeID(4))
    NodeSourceExt(NodeID(5)) = NodeSourceExtEqui(1,0,0,1) * CellLocNodes_Volumes(NodeID(5))
    NodeSourceExt(NodeID(6)) = NodeSourceExtEqui(1,1,0,1) * CellLocNodes_Volumes(NodeID(6))
    NodeSourceExt(NodeID(7)) = NodeSourceExtEqui(1,1,1,1) * CellLocNodes_Volumes(NodeID(7))
    NodeSourceExt(NodeID(8)) = NodeSourceExtEqui(1,0,1,1) * CellLocNodes_Volumes(NodeID(8)) 
  END ASSOCIATE
END DO
DEALLOCATE(U_local)
END IF ! DG_SourceExtExists

END SUBROUTINE ReadNodeSourceExtFromHDF5


END MODULE MOD_Restart_Tools
