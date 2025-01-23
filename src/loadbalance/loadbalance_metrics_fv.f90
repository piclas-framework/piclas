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

MODULE MOD_LoadBalance_Metrics_FV
!===================================================================================================================================
!> \brief This module contains routines for computing the geometries volume and surface metric terms.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
#if USE_LOADBALANCE && USE_FV
INTERFACE MoveCoords_FV
  MODULE PROCEDURE MoveCoords_FV
END INTERFACE

INTERFACE MoveMetrics_FV
  MODULE PROCEDURE MoveMetrics_FV
END INTERFACE

PUBLIC::MoveCoords_FV
PUBLIC::MoveMetrics_FV
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> This routine rearranges the coordinates along the space-filling curve
!==================================================================================================================================
SUBROUTINE MoveCoords_FV()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
USE MOD_LoadBalance_Vars   ,ONLY: MPInElemSend,MPInElemRecv,MPIoffsetElemSend,MPIoffsetElemRecv
USE MOD_Mesh_Vars_FV       ,ONLY: Elem_xGP_FV
USE MOD_Mesh_Vars          ,ONLY: nElems
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE                   :: Elem_xGP_LB(:,:,:,:,:)
! Custom data type
INTEGER                            :: MPI_LENGTH(1)
TYPE(MPI_Datatype)                 :: MPI_STRUCT,MPI_TYPE(1)
INTEGER(KIND=MPI_ADDRESS_KIND)     :: MPI_DISPLACEMENT(1)
! Timer
REAL                               :: StartT,EndT,WallTime
!==================================================================================================================================

IF (PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)) THEN
  ! SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') ' Shift volume coordinates during loadbalance...'
  GETTIME(StartT)

  ALLOCATE(Elem_xGP_LB   (3,0:0,0:0,0:0,nElems))
  ASSOCIATE (&
          counts_send  => INT(MPInElemSend     ) ,&
          disp_send    => INT(MPIoffsetElemSend) ,&
          counts_recv  => INT(MPInElemRecv     ) ,&
          disp_recv    => INT(MPIoffsetElemRecv))
    ! Communicate Elem_xGP over MPI
    MPI_LENGTH       = 3
    MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
    MPI_TYPE         = MPI_DOUBLE_PRECISION
    CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
    CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)

    CALL MPI_ALLTOALLV(Elem_xGP_FV,counts_send,disp_send,MPI_STRUCT,Elem_xGP_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
  END ASSOCIATE
  DEALLOCATE(Elem_xGP_FV)
  CALL MOVE_ALLOC(Elem_xGP_LB,Elem_xGP_FV)

  GETTIME(EndT)
  WallTime = EndT-StartT
  ! CALL DisplayMessageAndTime(WallTime,'DONE',DisplayDespiteLB=.TRUE.,DisplayLine=.FALSE.)
END IF

END SUBROUTINE MoveCoords_FV


SUBROUTINE MoveMetrics_FV()
!===================================================================================================================================
!> This routine rearranges the metrics along the space-filling curve
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
USE MOD_LoadBalance_Vars   ,ONLY: MPInElemSend,MPInElemRecv,MPIoffsetElemSend,MPIoffsetElemRecv
USE MOD_Mesh_Vars          ,ONLY: nElems
USE MOD_Mesh_Vars_FV       ,ONLY: JaCL_N_PP_1,Metrics_fTilde_FV,Metrics_gTilde_FV,Metrics_hTilde_FV,XCL_N_PP_1
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Elements
REAL,ALLOCATABLE                   :: XCL_N_LB(         :,:,:,:,:)
REAL,ALLOCATABLE                   :: JaCL_N_LB(      :,:,:,:,:,:)
REAL,ALLOCATABLE                   :: Metrics_fTilde_LB(:,:,:,:,:)
REAL,ALLOCATABLE                   :: Metrics_gTilde_LB(:,:,:,:,:)
REAL,ALLOCATABLE                   :: Metrics_hTilde_LB(:,:,:,:,:)
! Custom data type
INTEGER                            :: MPI_LENGTH(1)
TYPE(MPI_Datatype)                 :: MPI_TYPE(1),MPI_STRUCT
INTEGER(KIND=MPI_ADDRESS_KIND)     :: MPI_DISPLACEMENT(1)
! Timer
REAL                               :: StartT,EndT,WallTime
! !===================================================================================================================================
!
IF (PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)) THEN
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') ' Shift FV mesh metrics during loadbalance...'
  GETTIME(StartT)

  ! volume data
  ALLOCATE(       XCL_N_LB(  3,0:PP_1   ,0:PP_1   ,0:PP_1   ,nElems))
  ASSOCIATE (&
          counts_send  => INT(MPInElemSend     ) ,&
          disp_send    => INT(MPIoffsetElemSend) ,&
          counts_recv  => INT(MPInElemRecv     ) ,&
          disp_recv    => INT(MPIoffsetElemRecv))
    ! Communicate dXCL_N over MPI
    MPI_LENGTH       = 3*(PP_1+1)**3
    MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
    MPI_TYPE         = MPI_DOUBLE_PRECISION
    CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
    CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)

    CALL MPI_ALLTOALLV(XCL_N_PP_1,counts_send,disp_send,MPI_STRUCT,XCL_N_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
  END ASSOCIATE
  DEALLOCATE(XCL_N_PP_1)
  CALL MOVE_ALLOC(XCL_N_LB,XCL_N_PP_1)


  ALLOCATE(      JaCL_N_LB(3,3,0:PP_1   ,0:PP_1   ,0:PP_1   ,nElems))
  ASSOCIATE (&
          counts_send  => INT(MPInElemSend     ) ,&
          disp_send    => INT(MPIoffsetElemSend) ,&
          counts_recv  => INT(MPInElemRecv     ) ,&
          disp_recv    => INT(MPIoffsetElemRecv))
    ! Communicate JaCL_N over MPI
    MPI_LENGTH       = 3*3*(PP_1+1)**3
    MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
    MPI_TYPE         = MPI_DOUBLE_PRECISION
    CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
    CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)

    CALL MPI_ALLTOALLV(JaCL_N_PP_1,counts_send,disp_send,MPI_STRUCT,JaCL_N_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
  END ASSOCIATE
  DEALLOCATE(JaCL_N_PP_1)
  CALL MOVE_ALLOC(JaCL_N_LB,JaCL_N_PP_1)

  ASSOCIATE (&
          counts_send  => INT(MPInElemSend     ) ,&
          disp_send    => INT(MPIoffsetElemSend) ,&
          counts_recv  => INT(MPInElemRecv     ) ,&
          disp_recv    => INT(MPIoffsetElemRecv))
    ! Communicate metrics over MPI
    MPI_LENGTH       = 3
    MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
    MPI_TYPE         = MPI_DOUBLE_PRECISION
    CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
    CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)

    ALLOCATE(Metrics_fTilde_LB(3,0:0   ,0:0   ,0:0   ,nElems))
    CALL MPI_ALLTOALLV(Metrics_fTilde_FV,counts_send,disp_send,MPI_STRUCT,Metrics_fTilde_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
    DEALLOCATE(Metrics_fTilde_FV)
    CALL MOVE_ALLOC(Metrics_fTilde_LB,Metrics_fTilde_FV)

    ALLOCATE(Metrics_gTilde_LB(3,0:0   ,0:0   ,0:0   ,nElems))
    CALL MPI_ALLTOALLV(Metrics_gTilde_FV,counts_send,disp_send,MPI_STRUCT,Metrics_gTilde_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
    DEALLOCATE(Metrics_gTilde_FV)
    CALL MOVE_ALLOC(Metrics_gTilde_LB,Metrics_gTilde_FV)

    ALLOCATE(Metrics_hTilde_LB(3,0:0   ,0:0   ,0:0   ,nElems))
    CALL MPI_ALLTOALLV(Metrics_hTilde_FV,counts_send,disp_send,MPI_STRUCT,Metrics_hTilde_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
    DEALLOCATE(Metrics_hTilde_FV)
    CALL MOVE_ALLOC(Metrics_hTilde_LB,Metrics_hTilde_FV)
  END ASSOCIATE

  GETTIME(EndT)
  WallTime = EndT-StartT
  CALL DisplayMessageAndTime(WallTime,'DONE',DisplayDespiteLB=.TRUE.,DisplayLine=.FALSE.)
END IF

END SUBROUTINE MoveMetrics_FV
#endif /*USE_LOADBALANCE && USE_FV*/

END MODULE MOD_LoadBalance_Metrics_FV
