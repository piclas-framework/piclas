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

MODULE MOD_LoadBalance_Metrics
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
#if USE_LOADBALANCE
INTERFACE ExchangeVolMesh
  MODULE PROCEDURE ExchangeVolMesh
END INTERFACE

INTERFACE ExchangeMetrics
  MODULE PROCEDURE ExchangeMetrics
END INTERFACE

PUBLIC::ExchangeVolMesh
PUBLIC::ExchangeMetrics
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> This routine rearranges the coordinates along the space-filling curve
!==================================================================================================================================
SUBROUTINE ExchangeVolMesh()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
USE MOD_LoadBalance_Vars   ,ONLY: MPInElemSend,MPInElemRecv,MPIoffsetElemSend,MPIoffsetElemRecv
USE MOD_Mesh_Vars          ,ONLY: nElems,N_VolMesh,offSetElem
USE MOD_LoadBalance_Vars   ,ONLY: nElemsOld,offsetElemOld
USE MOD_DG_Vars            ,ONLY: N_DG_Mapping
USE MOD_Interpolation_Vars ,ONLY: Nmax
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE                   :: VolMesh(:,:,:,:,:),VolMesh_LB(:,:,:,:,:)
! Custom data type
INTEGER                            :: Nloc,i,j,k,iElem
INTEGER                            :: MPI_LENGTH(1)
TYPE(MPI_Datatype)                 :: MPI_STRUCT,MPI_TYPE(1)
INTEGER(KIND=MPI_ADDRESS_KIND)     :: MPI_DISPLACEMENT(1)
! Timer
REAL                               :: StartT,EndT,WallTime
!==================================================================================================================================

IF (PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)) THEN
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') ' Exchange N_VolMesh between processes during loadbalance...'
  GETTIME(StartT)

  ! Allocate old distribution of data with size nElemsOld
  ALLOCATE(VolMesh(16,0:NMax,0:NMax,0:NMax,nElemsOld))
  DO iElem = 1, nElemsOld
    Nloc = N_DG_Mapping(2,iElem+offSetElemOld)
    IF(Nloc.EQ.Nmax)THEN
      VolMesh( 1:3 ,:,:,:,iElem) = N_VolMesh(iElem)%Elem_xGP(      :,:,:,:)
      VolMesh( 4:6 ,:,:,:,iElem) = N_VolMesh(iElem)%XCL_N(         :,:,:,:)
      VolMesh( 7:9 ,:,:,:,iElem) = N_VolMesh(iElem)%Metrics_fTilde(:,:,:,:)
      VolMesh(10:12,:,:,:,iElem) = N_VolMesh(iElem)%Metrics_gTilde(:,:,:,:)
      VolMesh(13:15,:,:,:,iElem) = N_VolMesh(iElem)%Metrics_hTilde(:,:,:,:)
      VolMesh(16   ,:,:,:,iElem) = N_VolMesh(iElem)%sJ(              :,:,:)
    ELSE
      VolMesh(:,:,:,:,iElem) = 0.
      DO k=0,Nloc
        DO i=0,Nloc
          DO j=0,Nloc
            VolMesh( 1:3 ,i,j,k,iElem) = N_VolMesh(iElem)%Elem_xGP(      :,i,j,k)
            VolMesh( 4:6 ,i,j,k,iElem) = N_VolMesh(iElem)%XCL_N(         :,i,j,k)
            VolMesh( 7:9 ,i,j,k,iElem) = N_VolMesh(iElem)%Metrics_fTilde(:,i,j,k)
            VolMesh(10:12,i,j,k,iElem) = N_VolMesh(iElem)%Metrics_gTilde(:,i,j,k)
            VolMesh(13:15,i,j,k,iElem) = N_VolMesh(iElem)%Metrics_hTilde(:,i,j,k)
            VolMesh(16   ,i,j,k,iElem) = N_VolMesh(iElem)%sJ(              i,j,k)
          END DO
        END DO
      END DO
    END IF ! Nloc.Eq.Nmax
  END DO ! iElem = 1, nElems

  ! Allocate new distribution of data with size nElems
  ALLOCATE(VolMesh_LB(16,0:NMax,0:NMax,0:NMax,nElems))
  ASSOCIATE (&
          counts_send  => INT(MPInElemSend     ) ,&
          disp_send    => INT(MPIoffsetElemSend) ,&
          counts_recv  => INT(MPInElemRecv     ) ,&
          disp_recv    => INT(MPIoffsetElemRecv))
    ! Communicate VolMesh over MPI
    MPI_LENGTH       = 16*(Nmax+1)**3
    MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
    MPI_TYPE         = MPI_DOUBLE_PRECISION
    CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
    CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)

    CALL MPI_ALLTOALLV(VolMesh,counts_send,disp_send,MPI_STRUCT,VolMesh_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
  END ASSOCIATE
  CALL MOVE_ALLOC(VolMesh_LB,VolMesh)

  DEALLOCATE(N_VolMesh) ! N_VolMesh(1:nElemsOld)
  ! Extract data
  ALLOCATE(N_VolMesh(1:nElems))
  DO iElem = 1, nElems
    Nloc = N_DG_Mapping(2,iElem+offSetElem)
    ALLOCATE(N_VolMesh(iElem)%Elem_xGP(      3,0:Nloc,0:Nloc,0:Nloc))
    ALLOCATE(N_VolMesh(iElem)%XCL_N(         3,0:Nloc,0:Nloc,0:Nloc))
    ALLOCATE(N_VolMesh(iElem)%Metrics_fTilde(3,0:Nloc,0:Nloc,0:Nloc))
    ALLOCATE(N_VolMesh(iElem)%Metrics_gTilde(3,0:Nloc,0:Nloc,0:Nloc))
    ALLOCATE(N_VolMesh(iElem)%Metrics_hTilde(3,0:Nloc,0:Nloc,0:Nloc))
    ALLOCATE(N_VolMesh(iElem)%sJ            (  0:Nloc,0:Nloc,0:Nloc))
    DO k=0,Nloc
      DO j=0,Nloc
        DO i=0,Nloc
          N_VolMesh(iElem)%Elem_xGP(      :,i,j,k) = VolMesh( 1:3 ,i,j,k,iElem)
          N_VolMesh(iElem)%XCL_N(         :,i,j,k) = VolMesh( 4:6 ,i,j,k,iElem)
          N_VolMesh(iElem)%Metrics_fTilde(:,i,j,k) = VolMesh( 7:9 ,i,j,k,iElem)
          N_VolMesh(iElem)%Metrics_gTilde(:,i,j,k) = VolMesh(10:12,i,j,k,iElem)
          N_VolMesh(iElem)%Metrics_hTilde(:,i,j,k) = VolMesh(13:15,i,j,k,iElem)
          N_VolMesh(iElem)%sJ(              i,j,k) = VolMesh(16   ,i,j,k,iElem)
        END DO
      END DO
    END DO
  END DO ! iElem = 1, nElems
  DEALLOCATE(VolMesh)

  GETTIME(EndT)
  WallTime = EndT-StartT
  CALL DisplayMessageAndTime(WallTime,'DONE',DisplayDespiteLB=.TRUE.,DisplayLine=.FALSE.)
END IF

END SUBROUTINE ExchangeVolMesh


SUBROUTINE ExchangeMetrics()
!===================================================================================================================================
!> This routine rearranges the metrics along the space-filling curve
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
USE MOD_LoadBalance_Vars   ,ONLY: MPInElemSend,MPInElemRecv,MPIoffsetElemSend,MPIoffsetElemRecv
! USE MOD_LoadBalance_Vars   ,ONLY: MPInSideSend,MPInSideRecv,MPIoffsetSideSend,MPIoffsetSideRecv
USE MOD_Mesh_Vars          ,ONLY: nElems,NGeo,offSetElem
USE MOD_LoadBalance_Vars   ,ONLY: nElemsOld,offsetElemOld
!USE MOD_Mesh_Vars          ,ONLY: JaCL_N,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,XCL_N!,dXCL_N
! USE MOD_Mesh_Vars          ,ONLY: Face_xGP,NormVec,TangVec1,TangVec2,SurfElem,Ja_Face
!USE MOD_Mesh_Vars          ,ONLY: sJ!,detJac_Ref
USE MOD_Mesh_Vars          ,ONLY: NGeo,XCL_NGeo
#ifdef PARTICLES
USE MOD_Mesh_Vars          ,ONLY: dXCL_NGeo
#endif /*PARTICLES*/
USE MOD_Mesh_Vars          ,ONLY: nElems,N_VolMesh2
USE MOD_DG_Vars            ,ONLY: N_DG_Mapping
USE MOD_Interpolation_Vars ,ONLY: Nmax
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Elements
!REAL,ALLOCATABLE                   :: XCL_N_LB(         :,:,:,:,:)
!REAL,ALLOCATABLE                   :: XCL_N   (         :,:,:,:,:)
! REAL,ALLOCATABLE                   :: dXCL_N_LB(      :,:,:,:,:,:)
REAL,ALLOCATABLE                   :: JaCL_N_LB(      :,:,:,:,:,:)
REAL,ALLOCATABLE                   :: JaCL_N   (      :,:,:,:,:,:)
!REAL,ALLOCATABLE                   :: Metrics_Tilde_LB(:,:,:,:,:)
!REAL,ALLOCATABLE                   :: Metrics_Tilde   (:,:,:,:,:)
!REAL,ALLOCATABLE                   :: sJ_LB            (  :,:,:,:)
!REAL,ALLOCATABLE                   :: sJ               (  :,:,:,:)
REAL,ALLOCATABLE                   :: XCL_NGeo_LB(      :,:,:,:,:)
#ifdef PARTICLES
REAL,ALLOCATABLE                   :: dXCL_NGeo_LB(   :,:,:,:,:,:)
#endif /*PARTICLES*/
! REAL,ALLOCATABLE                   :: DetJac_Ref_LB    (:,:,:,:,:)
! Sides
! REAL,ALLOCATABLE                   :: Face_xGP_LB    (  :,:,:,:)
! REAL,ALLOCATABLE                   :: NormVec_LB     (  :,:,:,:)
! REAL,ALLOCATABLE                   :: TangVec1_LB    (  :,:,:,:)
! REAL,ALLOCATABLE                   :: TangVec2_LB    (  :,:,:,:)
! REAL,ALLOCATABLE                   :: SurfElem_LB    (    :,:,:)
! REAL,ALLOCATABLE                   ::      Ja_Face_LB(:,:,:,:,:)
! Custom data type
INTEGER                            :: Nloc,i,j,k,iElem
INTEGER                            :: MPI_LENGTH(1)
TYPE(MPI_Datatype)                 :: MPI_TYPE(1),MPI_STRUCT
INTEGER(KIND=MPI_ADDRESS_KIND)     :: MPI_DISPLACEMENT(1)
! Timer
REAL                               :: StartT,EndT,WallTime
! !===================================================================================================================================
!
IF (PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance)) THEN
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') ' Shift Ngeo mesh metrics and N_VolMesh2 during loadbalance...'
  GETTIME(StartT)

!  ! Allocate old distribution of data with size nElemsOld
!  ALLOCATE(XCL_N(3,0:NMax,0:NMax,0:NMax,nElemsOld))
!  DO iElem = 1, nElemsOld
!    Nloc = N_DG_Mapping(2,iElem+offSetElemOld)
!    IF(Nloc.EQ.Nmax)THEN
!      XCL_N(:,:,:,:,iElem) = N_VolMesh(iElem)%XCL_N(:,:,:,:)
!    ELSE
!      XCL_N(:,:,:,:,iElem) = 0.
!      DO k=0,Nloc
!        DO i=0,Nloc
!          DO j=0,Nloc
!            XCL_N(:,i,j,k,iElem) = N_VolMesh(iElem)%XCL_N(:,i,j,k)
!          END DO
!        END DO
!      END DO
!    END IF ! Nloc.Eq.Nmax
!  END DO ! iElem = 1, nElems
!
!  ! Allocate new distribution of data with size nElems
!  ALLOCATE(XCL_N_LB(3,0:NMax,0:NMax,0:NMax,nElems))
!  ASSOCIATE (&
!          counts_send  => INT(MPInElemSend     ) ,&
!          disp_send    => INT(MPIoffsetElemSend) ,&
!          counts_recv  => INT(MPInElemRecv     ) ,&
!          disp_recv    => INT(MPIoffsetElemRecv))
!    ! Communicate XCL_N over MPI
!    MPI_LENGTH       = 3*(Nmax+1)**3
!    MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
!    MPI_TYPE         = MPI_DOUBLE_PRECISION
!    CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
!    CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)
!
!    CALL MPI_ALLTOALLV(XCL_N,counts_send,disp_send,MPI_STRUCT,XCL_N_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
!  END ASSOCIATE
!  CALL MOVE_ALLOC(XCL_N_LB,XCL_N)
!
!  ! Extract data
!  DO iElem = 1, nElems
!    Nloc = N_DG_Mapping(2,iElem+offSetElem)
!    ALLOCATE(N_VolMesh(iElem)%XCL_N(3,0:Nloc,0:Nloc,0:Nloc))
!    DO k=0,Nloc
!      DO j=0,Nloc
!        DO i=0,Nloc
!          N_VolMesh(iElem)%XCL_N(:,i,j,k) = XCL_N(1:3,i,j,k,iElem)
!        END DO
!      END DO
!    END DO
!  END DO ! iElem = 1, nElems
!  DEALLOCATE(XCL_N)

  ! ALLOCATE(      dXCL_N_LB(3,3,0:PP_N   ,0:PP_N   ,0:PP_N   ,nElems))
  ! ASSOCIATE (&
  !         counts_send  => INT(MPInElemSend     ) ,&
  !         disp_send    => INT(MPIoffsetElemSend) ,&
  !         counts_recv  => INT(MPInElemRecv     ) ,&
  !         disp_recv    => INT(MPIoffsetElemRecv))
  !   MPI_LENGTH       = 3*3*(PP_N+1)**3
  !   ! Communicate dXCL_N over MPI
  !   MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
  !   MPI_TYPE         = MPI_DOUBLE_PRECISION
  !   CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
  !   CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)
  !
  !   CALL MPI_ALLTOALLV(dXCL_N,counts_send,disp_send,MPI_STRUCT,dXCL_N_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
  ! END ASSOCIATE
  ! DEALLOCATE(dXCL_N)
  ! CALL MOVE_ALLOC(dXCL_N_LB,dXCL_N)

  ALLOCATE(XCL_NGeo_LB(1:3,0:NGeo,0:NGeo,0:NGeo,1:nElems))
  ASSOCIATE (&
          counts_send  => INT(MPInElemSend     ) ,&
          disp_send    => INT(MPIoffsetElemSend) ,&
          counts_recv  => INT(MPInElemRecv     ) ,&
          disp_recv    => INT(MPIoffsetElemRecv))
    ! Communicate XCL_NGeo over MPI
    MPI_LENGTH       = 3*(NGeo+1)**3
    MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
    MPI_TYPE         = MPI_DOUBLE_PRECISION
    CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
    CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)

    CALL MPI_ALLTOALLV(XCL_NGeo,counts_send,disp_send,MPI_STRUCT,XCL_NGeo_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
  END ASSOCIATE
  DEALLOCATE(XCL_NGeo)
  CALL MOVE_ALLOC(XCL_NGeo_LB,XCL_NGeo)

#ifdef PARTICLES
  ALLOCATE(dXCL_NGeo_LB(1:3,1:3,0:NGeo,0:NGeo,0:NGeo,1:nElems))
  ASSOCIATE (&
          counts_send  => INT(MPInElemSend     ) ,&
          disp_send    => INT(MPIoffsetElemSend) ,&
          counts_recv  => INT(MPInElemRecv     ) ,&
          disp_recv    => INT(MPIoffsetElemRecv))
    ! Communicate dXCL_NGeo over MPI
    MPI_LENGTH       = 3*3*(NGeo+1)**3
    MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
    MPI_TYPE         = MPI_DOUBLE_PRECISION
    CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
    CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)

    CALL MPI_ALLTOALLV(dXCL_NGeo,counts_send,disp_send,MPI_STRUCT,dXCL_NGeo_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
  END ASSOCIATE
  DEALLOCATE(dXCL_NGeo)
  CALL MOVE_ALLOC(dXCL_NGeo_LB,dXCL_NGeo)
#endif /*PARTICLES*/


  ! Allocate old distribution of data with size nElemsOld
  ALLOCATE(JaCL_N(3,3,0:NMax,0:NMax,0:NMax,nElemsOld))
  DO iElem = 1, nElemsOld
    Nloc = N_DG_Mapping(2,iElem+offSetElemOld)
    IF(Nloc.EQ.Nmax)THEN
      JaCL_N(:,:,:,:,:,iElem) = N_VolMesh2(iElem)%JaCL_N(:,:,:,:,:)
    ELSE
      JaCL_N(:,:,:,:,:,iElem) = 0.
      DO k=0,Nloc
        DO i=0,Nloc
          DO j=0,Nloc
            JaCL_N(:,:,i,j,k,iElem) = N_VolMesh2(iElem)%JaCL_N(:,:,i,j,k)
          END DO
        END DO
      END DO
    END IF ! Nloc.Eq.Nmax
  END DO ! iElem = 1, nElems

  ! Allocate new distribution of data with size nElems
  ALLOCATE(JaCL_N_LB(3,3,0:NMax,0:NMax,0:NMax,nElems))
  ASSOCIATE (&
          counts_send  => INT(MPInElemSend     ) ,&
          disp_send    => INT(MPIoffsetElemSend) ,&
          counts_recv  => INT(MPInElemRecv     ) ,&
          disp_recv    => INT(MPIoffsetElemRecv))
    ! Communicate JaCL_N over MPI
    MPI_LENGTH       = 3*3*(Nmax+1)**3
    MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
    MPI_TYPE         = MPI_DOUBLE_PRECISION
    CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
    CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)

    CALL MPI_ALLTOALLV(JaCL_N,counts_send,disp_send,MPI_STRUCT,JaCL_N_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
  END ASSOCIATE
  CALL MOVE_ALLOC(JaCL_N_LB,JaCL_N)

  DEALLOCATE(N_VolMesh2) ! N_VolMesh2(1:nElemsOld)
  ! Extract data
  ALLOCATE(N_VolMesh2(1:nElems))
  DO iElem = 1, nElems
    Nloc = N_DG_Mapping(2,iElem+offSetElem)
    ALLOCATE(N_VolMesh2(iElem)%JaCL_N(3,3,0:Nloc,0:Nloc,0:Nloc))
    DO k=0,Nloc
      DO j=0,Nloc
        DO i=0,Nloc
          N_VolMesh2(iElem)%JaCL_N(:,:,i,j,k) = JaCL_N(1:3,1:3,i,j,k,iElem)
        END DO
      END DO
    END DO
  END DO ! iElem = 1, nElems
  DEALLOCATE(JaCL_N)

!  ! Allocate old distribution of data with size nElemsOld
!  ALLOCATE(Metrics_Tilde(9,0:NMax,0:NMax,0:NMax,nElemsOld))
!  DO iElem = 1, nElemsOld
!    Nloc = N_DG_Mapping(2,iElem+offSetElemOld)
!    IF(Nloc.EQ.Nmax)THEN
!      Metrics_Tilde(1:3,:,:,:,iElem) = N_VolMesh(iElem)%Metrics_fTilde(:,:,:,:)
!      Metrics_Tilde(4:6,:,:,:,iElem) = N_VolMesh(iElem)%Metrics_gTilde(:,:,:,:)
!      Metrics_Tilde(7:9,:,:,:,iElem) = N_VolMesh(iElem)%Metrics_hTilde(:,:,:,:)
!    ELSE
!      Metrics_Tilde(:,:,:,:,iElem) = 0.
!      DO k=0,Nloc
!        DO i=0,Nloc
!          DO j=0,Nloc
!            Metrics_Tilde(1:3,i,j,k,iElem) = N_VolMesh(iElem)%Metrics_fTilde(:,i,j,k)
!            Metrics_Tilde(4:6,i,j,k,iElem) = N_VolMesh(iElem)%Metrics_gTilde(:,i,j,k)
!            Metrics_Tilde(7:9,i,j,k,iElem) = N_VolMesh(iElem)%Metrics_hTilde(:,i,j,k)
!          END DO
!        END DO
!      END DO
!    END IF ! Nloc.Eq.Nmax
!  END DO ! iElem = 1, nElems
!
!  ! Allocate new distribution of data with size nElems
!  ALLOCATE(Metrics_Tilde_LB(9,0:NMax,0:NMax,0:NMax,nElems))
!  ASSOCIATE (&
!          counts_send  => INT(MPInElemSend     ) ,&
!          disp_send    => INT(MPIoffsetElemSend) ,&
!          counts_recv  => INT(MPInElemRecv     ) ,&
!          disp_recv    => INT(MPIoffsetElemRecv))
!    ! Communicate Metrics_Tilde over MPI
!    MPI_LENGTH       = 9*(Nmax+1)**3
!    MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
!    MPI_TYPE         = MPI_DOUBLE_PRECISION
!    CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
!    CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)
!
!    CALL MPI_ALLTOALLV(Metrics_Tilde,counts_send,disp_send,MPI_STRUCT,Metrics_Tilde_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
!  END ASSOCIATE
!  CALL MOVE_ALLOC(Metrics_Tilde_LB,Metrics_Tilde)
!
!  ! Extract data
!  DO iElem = 1, nElems
!    Nloc = N_DG_Mapping(2,iElem+offSetElem)
!    ALLOCATE(N_VolMesh(iElem)%Metrics_fTilde(3,0:Nloc,0:Nloc,0:Nloc))
!    ALLOCATE(N_VolMesh(iElem)%Metrics_gTilde(3,0:Nloc,0:Nloc,0:Nloc))
!    ALLOCATE(N_VolMesh(iElem)%Metrics_hTilde(3,0:Nloc,0:Nloc,0:Nloc))
!    DO k=0,Nloc
!      DO j=0,Nloc
!        DO i=0,Nloc
!          N_VolMesh(iElem)%Metrics_fTilde(:,i,j,k) = Metrics_Tilde(1:3,i,j,k,iElem)
!          N_VolMesh(iElem)%Metrics_gTilde(:,i,j,k) = Metrics_Tilde(4:6,i,j,k,iElem)
!          N_VolMesh(iElem)%Metrics_hTilde(:,i,j,k) = Metrics_Tilde(7:9,i,j,k,iElem)
!        END DO
!      END DO
!    END DO
!  END DO ! iElem = 1, nElems
!  DEALLOCATE(Metrics_Tilde)

!  ! Allocate old distribution of data with size nElemsOld
!  ALLOCATE(sJ(0:NMax,0:NMax,0:NMax,nElemsOld))
!  DO iElem = 1, nElemsOld
!    Nloc = N_DG_Mapping(2,iElem+offSetElemOld)
!    IF(Nloc.EQ.Nmax)THEN
!      sJ(:,:,:,iElem) = N_VolMesh(iElem)%sJ(:,:,:)
!    ELSE
!      sJ(:,:,:,iElem) = 0.
!      DO k=0,Nloc
!        DO i=0,Nloc
!          DO j=0,Nloc
!            sJ(i,j,k,iElem) = N_VolMesh(iElem)%sJ(i,j,k)
!          END DO
!        END DO
!      END DO
!    END IF ! Nloc.Eq.Nmax
!  END DO ! iElem = 1, nElems
!
!  ! Allocate new distribution of data with size nElems
!  ALLOCATE(sJ_LB(0:NMax,0:NMax,0:NMax,nElems))
!  ASSOCIATE (&
!          counts_send  => INT(MPInElemSend     ) ,&
!          disp_send    => INT(MPIoffsetElemSend) ,&
!          counts_recv  => INT(MPInElemRecv     ) ,&
!          disp_recv    => INT(MPIoffsetElemRecv))
!    ! Communicate sJ over MPI
!    MPI_LENGTH       = (Nmax+1)**3
!    MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
!    MPI_TYPE         = MPI_DOUBLE_PRECISION
!    CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
!    CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)
!
!    CALL MPI_ALLTOALLV(sJ,counts_send,disp_send,MPI_STRUCT,sJ_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
!  END ASSOCIATE
!  CALL MOVE_ALLOC(sJ_LB,sJ)
!
!  ! Extract data
!  DO iElem = 1, nElems
!    Nloc = N_DG_Mapping(2,iElem+offSetElem)
!    ALLOCATE(N_VolMesh(iElem)%sJ(0:Nloc,0:Nloc,0:Nloc))
!    DO k=0,Nloc
!      DO j=0,Nloc
!        DO i=0,Nloc
!          N_VolMesh(iElem)%sJ(i,j,k) = sJ(i,j,k,iElem)
!        END DO
!      END DO
!    END DO
!  END DO ! iElem = 1, nElems
!  DEALLOCATE(sJ)

  ! ! surface data
  ! ASSOCIATE (&
  !         counts_send  => INT(MPInSideSend     ) ,&
  !         disp_send    => INT(MPIoffsetSideSend) ,&
  !         counts_recv  => INT(MPInSideRecv     ) ,&
  !         disp_recv    => INT(MPIoffsetSideRecv))
  !   ! Communicate scaled Jacobians over MPI
  !   IPWRITE(*,*) 'MPInSideSend,MPIoffsetSideSend:', MPInSideSend,MPIoffsetSideSend
  !   IPWRITE(*,*) 'MPInSideRecv,MPIoffsetSideRecv:', MPInSideSend,MPIoffsetSideSend
  !   MPI_LENGTH       = 3*(PP_N+1)**2
  !   MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
  !   MPI_TYPE         = MPI_DOUBLE_PRECISION
  !   CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
  !   CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)
  !
  !   ! Communicate side metrics over MPI
  !   ALLOCATE(Face_xGP_LB      (3,0:PP_N   ,0:PP_N   ,1:nSides))
  !   CALL MPI_ALLTOALLV(Face_xGP,counts_send,disp_send,MPI_STRUCT,Face_xGP_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
  !   DEALLOCATE(Face_xGP)
  !   CALL MOVE_ALLOC(Face_xGP_LB,Face_xGP)
  !
  !   ALLOCATE(NormVec_LB       (3,0:PP_N   ,0:PP_N   ,1:nSides))
  !   CALL MPI_ALLTOALLV(NormVec,counts_send,disp_send,MPI_STRUCT,NormVec_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
  !   DEALLOCATE(NormVec)
  !   CALL MOVE_ALLOC(NormVec_LB,NormVec)
  !
  !   ALLOCATE(TangVec1_LB      (3,0:PP_N   ,0:PP_N   ,1:nSides))
  !   CALL MPI_ALLTOALLV(TangVec1,counts_send,disp_send,MPI_STRUCT,TangVec1_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
  !   DEALLOCATE(TangVec1)
  !   CALL MOVE_ALLOC(TangVec1_LB,TangVec1)
  !
  !   ALLOCATE(TangVec2_LB      (3,0:PP_N   ,0:PP_N   ,1:nSides))
  !   CALL MPI_ALLTOALLV(TangVec2,counts_send,disp_send,MPI_STRUCT,TangVec2_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
  !   DEALLOCATE(TangVec2)
  !   CALL MOVE_ALLOC(TangVec2_LB,TangVec2)
  ! END ASSOCIATE
  !
  ! ASSOCIATE (&
  !         counts_send  => INT(MPInSideSend     ) ,&
  !         disp_send    => INT(MPIoffsetSideSend) ,&
  !         counts_recv  => INT(MPInSideRecv     ) ,&
  !         disp_recv    => INT(MPIoffsetSideRecv))
  !   ! Communicate scaled Jacobians over MPI
  !   MPI_LENGTH       = (PP_N+1)**2
  !   MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
  !   MPI_TYPE         = MPI_DOUBLE_PRECISION
  !   CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
  !   CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)
  !   ALLOCATE(SurfElem_LB      (  0:PP_N   ,0:PP_N   ,1:nSides))
  !   CALL MPI_ALLTOALLV(SurfElem,counts_send,disp_send,MPI_STRUCT,SurfElem_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
  !   DEALLOCATE(SurfElem)
  !   CALL MOVE_ALLOC(SurfElem_LB,SurfElem)
  ! END ASSOCIATE
  !
  ! ASSOCIATE (&
  !         counts_send  => INT(MPInSideSend     ) ,&
  !         disp_send    => INT(MPIoffsetSideSend) ,&
  !         counts_recv  => INT(MPInSideRecv     ) ,&
  !         disp_recv    => INT(MPIoffsetSideRecv))
  !   ! Communicate scaled Jacobians over MPI
  !   MPI_LENGTH       = 3*3*(PP_N+1)**2
  !   MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
  !   MPI_TYPE         = MPI_DOUBLE_PRECISION
  !   CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
  !   CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)
  !   ALLOCATE(     Ja_Face_LB(3,3,0:PP_N   ,0:PP_N   ,1:nSides))
  !   CALL MPI_ALLTOALLV(Ja_Face,counts_send,disp_send,MPI_STRUCT,Ja_Face_LB,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
  !   DEALLOCATE(Ja_Face)
  !   CALL MOVE_ALLOC(Ja_Face_LB,Ja_Face)
  ! END ASSOCIATE

  GETTIME(EndT)
  WallTime = EndT-StartT
  CALL DisplayMessageAndTime(WallTime,'DONE',DisplayDespiteLB=.TRUE.,DisplayLine=.FALSE.)
END IF

END SUBROUTINE ExchangeMetrics
#endif /*USE_LOADBALANCE*/

END MODULE MOD_LoadBalance_Metrics