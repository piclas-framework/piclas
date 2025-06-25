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

MODULE MOD_MPI_HDG
!===================================================================================================================================
! Module containing subroutines and functions that are required only for Hybrid Discontinuous Galerkin (HDG) methods
!===================================================================================================================================
! MODULES
#if USE_MPI
USE mpi_f08
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
#if USE_MPI && USE_HDG

#if USE_HDG
! !no interface for reshape inout vector
! INTERFACE Mask_MPIsides
!   MODULE PROCEDURE Mask_MPIsides
! END INTERFACE
#endif /*USE_HDG*/
PUBLIC :: Mask_MPIsides
PUBLIC :: StartReceiveMPISurfDataType
PUBLIC :: StartSendMPISurfDataType
PUBLIC :: FinishExchangeMPISurfDataType
#endif /*USE_MPI && USE_HDG*/
!===================================================================================================================================

CONTAINS

#if USE_MPI && USE_HDG
!===================================================================================================================================
!> communicate contribution from MPI slave sides to MPI master sides  and set slaves them to zero afterwards.
!> Select between different modes:
!> 1: mv (MatVec in hdg.f90)
!> 2: RHS_face  (HDGLinear and HDGNewton in hdg.f90)
!> 3: InvPrecondDiag (BuildPrecond in elem_mat.f90)
!> 4: Precond (BuildPrecond in elem_mat.f90)
!===================================================================================================================================
SUBROUTINE Mask_MPIsides(mode, iVar)
! MODULES
USE MOD_GLobals
USE MOD_MPI_Vars
USE MOD_Mesh_Vars ,ONLY: nSides
USE MOD_Mesh_Vars ,ONLY: nMPIsides_YOUR,nMPIsides,nMPIsides_MINE
USE MOD_HDG_Vars  ,ONLY: HDG_Surf_N
! USE MOD_MPI       ,ONLY: StartReceiveMPISurfDataType,StartSendMPISurfDataType,FinishExchangeMPISurfDataType

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: mode     !< Select mv, RHS_face, InvPrecondDiag or Precond
INTEGER, INTENT(IN), OPTIONAL :: iVar
!INTEGER,INTENT(IN   )       :: firstdim !< size of first dimention in array to be sent
!REAL   ,INTENT(INOUT) :: v(firstdim,nGP_face, nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL    :: vbuf(firstdim,nGP_Face,nMPISides_MINE)
INTEGER :: startbuf,endbuf, sendmode, iSide
!===================================================================================================================================

startbuf = nSides-nMPISides+1
endbuf   = nSides-nMPISides+nMPISides_MINE

SELECT CASE(TRIM(mode))
CASE('mv') ! Originally: CALL Mask_MPIsides(1,mv)

  ! HDG_Surf_N(SideID)%mv(iVar,:)
  IF(nMPIsides_MINE.GT.0) THEN
    DO iSide = startbuf, endbuf
      HDG_Surf_N(iSide)%buf(iVar,:) = HDG_Surf_N(iSide)%mv(iVar,:)
    END DO
  END IF
  sendmode = 3
  CALL StartReceiveMPISurfDataType(RecRequest_U , 2 , sendmode)
  CALL StartSendMPISurfDataType(  SendRequest_U , 2 , sendmode , iVar=iVar)
  CALL FinishExchangeMPISurfDataType(SendRequest_U, RecRequest_U, 2, sendmode, iVar=iVar)
  IF(nMPIsides_MINE.GT.0) THEN
    DO iSide = startbuf, endbuf
      HDG_Surf_N(iSide)%mv(iVar,:) = HDG_Surf_N(iSide)%mv(iVar,:) + HDG_Surf_N(iSide)%buf(iVar,:)
    END DO
  END IF
  IF(nMPIsides_YOUR.GT.0) THEN
    DO iSide = nSides-nMPIsides_YOUR+1, nSides
      HDG_Surf_N(iSide)%mv(iVar,:) = 0.
    END DO
  END IF

CASE('Z')

  IF(nMPIsides_MINE.GT.0) THEN
    DO iSide = startbuf, endbuf
      HDG_Surf_N(iSide)%buf(iVar,:) = HDG_Surf_N(iSide)%Z(iVar,:)
    END DO
  END IF
  sendmode = 6
  CALL StartReceiveMPISurfDataType(RecRequest_U , 2 , sendmode)
  CALL StartSendMPISurfDataType(  SendRequest_U , 2 , sendmode , iVar=iVar)
  CALL FinishExchangeMPISurfDataType(SendRequest_U, RecRequest_U, 2, sendmode, iVar=iVar)
  IF(nMPIsides_MINE.GT.0) THEN
    DO iSide = startbuf, endbuf
      HDG_Surf_N(iSide)%Z(iVar,:) = HDG_Surf_N(iSide)%Z(iVar,:) + HDG_Surf_N(iSide)%buf(iVar,:)
    END DO
  END IF
  IF(nMPIsides_YOUR.GT.0) THEN
    DO iSide = nSides-nMPIsides_YOUR+1, nSides
      HDG_Surf_N(iSide)%Z(iVar,:) = 0.
    END DO
  END IF

CASE('RHS_face') ! Originally: CALL Mask_MPIsides(PP_nVar,RHS_face)

  ! HDG_Surf_N(SideID)%RHS_face(iVar,:)
  IF(nMPIsides_MINE.GT.0) THEN
    DO iSide = startbuf, endbuf
      HDG_Surf_N(iSide)%buf(iVar,:) = HDG_Surf_N(iSide)%RHS_face(iVar,:)
    END DO
  END IF
  sendmode = 4
  CALL StartReceiveMPISurfDataType(RecRequest_U , 2 , sendmode)
  CALL StartSendMPISurfDataType(  SendRequest_U , 2 , sendmode , iVar=iVar)
  CALL FinishExchangeMPISurfDataType(SendRequest_U, RecRequest_U, 2, sendmode, iVar=iVar)
  IF(nMPIsides_MINE.GT.0) THEN
    DO iSide = startbuf, endbuf
      HDG_Surf_N(iSide)%RHS_face(iVar,:) = HDG_Surf_N(iSide)%RHS_face(iVar,:) + HDG_Surf_N(iSide)%buf(iVar,:)
    END DO
  END IF
  IF(nMPIsides_YOUR.GT.0) THEN
    DO iSide = nSides-nMPIsides_YOUR+1, nSides
      HDG_Surf_N(iSide)%RHS_face(iVar,:) = 0.
    END DO
  END IF

CASE('Precond') ! CALL Mask_MPIsides(nGP_face,Precond) for PrecondType=1

  ! HDG_Surf_N(SideID)%Precond
  ! NSideMin = N_SurfMesh(SideID)%NSideMin
  ! ALLOCATE(HDG_Surf_N(SideID)%Precond(nGP_face(NSideMin),nGP_face(NSideMin)))
  IF(nMPIsides_MINE.GT.0) THEN
    DO iSide = startbuf, endbuf
      HDG_Surf_N(iSide)%buf2(:,:) = HDG_Surf_N(iSide)%Precond(:,:)
    END DO
  END IF
  sendmode = 8
  CALL StartReceiveMPISurfDataType(RecRequest_U , 2 , sendmode)
  CALL StartSendMPISurfDataType(  SendRequest_U , 2 , sendmode , iVar=1)
  CALL FinishExchangeMPISurfDataType(SendRequest_U, RecRequest_U, 2, sendmode, iVar=1)
  IF(nMPIsides_MINE.GT.0) THEN
    DO iSide = startbuf, endbuf
      HDG_Surf_N(iSide)%Precond(:,:) = HDG_Surf_N(iSide)%Precond(:,:) + HDG_Surf_N(iSide)%buf2(:,:)
    END DO
  END IF
  IF(nMPIsides_YOUR.GT.0) THEN
    DO iSide = nSides-nMPIsides_YOUR+1, nSides
      HDG_Surf_N(iSide)%Precond(:,:) = 0.
    END DO
  END IF

CASE('InvPrecondDiag') ! Originally: CALL Mask_MPIsides(1,InvPrecondDiag) for PrecondType=2

  ! HDG_Surf_N(SideID)%InvPrecondDiag
  ! NSideMin = N_SurfMesh(SideID)%NSideMin
  ! ALLOCATE(HDG_Surf_N(SideID)%InvPrecondDiag(nGP_face(NSideMin)))
  IF(nMPIsides_MINE.GT.0) THEN
    DO iSide = startbuf, endbuf
      HDG_Surf_N(iSide)%buf(1,:) = HDG_Surf_N(iSide)%InvPrecondDiag(:)
    END DO
  END IF
  sendmode = 7
  CALL StartReceiveMPISurfDataType(RecRequest_U , 2 , sendmode)
  CALL StartSendMPISurfDataType(  SendRequest_U , 2 , sendmode , iVar=1)
  CALL FinishExchangeMPISurfDataType(SendRequest_U, RecRequest_U, 2, sendmode, iVar=1)
  IF(nMPIsides_MINE.GT.0) THEN
    DO iSide = startbuf, endbuf
      HDG_Surf_N(iSide)%InvPrecondDiag(:) = HDG_Surf_N(iSide)%InvPrecondDiag(:) + HDG_Surf_N(iSide)%buf(1,:)
    END DO
  END IF
  IF(nMPIsides_YOUR.GT.0) THEN
    DO iSide = nSides-nMPIsides_YOUR+1, nSides
      HDG_Surf_N(iSide)%InvPrecondDiag(:) = 0.
    END DO
  END IF

END SELECT

END SUBROUTINE Mask_MPIsides


!===================================================================================================================================
!> Subroutine does the receive operations for the face data that has to be exchanged between processors (type-based p-adaption).
!===================================================================================================================================
SUBROUTINE StartReceiveMPISurfDataType(MPIRequest, SendID, mode)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
USE MOD_HDG_Vars, ONLY: UseNSideMin
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: SendID                                                 !< defines the send / receive direction -> 1=send MINE
                                                                              !< / receive YOUR, 3=send YOUR / receive MINE
INTEGER,INTENT(IN)  :: mode
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(MPI_Request),INTENT(OUT) :: MPIRequest(nNbProcs)                                   !< communication handles
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
DO iNbProc=1,nNbProcs
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    IF (mode.EQ.1) THEN
      nRecVal = DataSizeSurfRecMax(iNbProc,SendID)
      CALL MPI_IRECV(SurfExchange(iNbProc)%SurfDataRecv(1:nRecVal),nRecVal,MPI_DOUBLE_PRECISION,  &
                      nbProc(iNbProc),0,MPI_COMM_PICLAS,MPIRequest(iNbProc),iError)
    ELSEIF (mode.EQ.8) THEN
      IF(UseNSideMin)THEN
        nRecVal = DataSizeSurfRecMin(iNbProc,SendID)**2
      ELSE
        nRecVal = DataSizeSurfRecMax(iNbProc,SendID)**2
      END IF
      CALL MPI_IRECV(SurfExchange(iNbProc)%SurfDataRecv2(1:nRecVal),nRecVal,MPI_DOUBLE_PRECISION,  &
                      nbProc(iNbProc),0,MPI_COMM_PICLAS,MPIRequest(iNbProc),iError)
    ELSE
      IF(UseNSideMin)THEN
        nRecVal = DataSizeSurfRecMin(iNbProc,SendID)
      ELSE
        nRecVal = DataSizeSurfRecMax(iNbProc,SendID)
      END IF
      CALL MPI_IRECV(SurfExchange(iNbProc)%SurfDataRecv(1:nRecVal),nRecVal,MPI_DOUBLE_PRECISION,  &
                      nbProc(iNbProc),0,MPI_COMM_PICLAS,MPIRequest(iNbProc),iError)
    END IF
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_IRECV',iError)
  ELSE
    MPIRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
END DO !iProc=1,nNBProcs
END SUBROUTINE StartReceiveMPISurfDataType


!===================================================================================================================================
!> See above, but for for send direction (type-based p-adaption).
!===================================================================================================================================
SUBROUTINE StartSendMPISurfDataType(MPIRequest,SendID, mode, iVar)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
USE MOD_DG_Vars   ,ONLY: U_Surf_N,DG_Elems_slave
USE MOD_Mesh_Vars ,ONLY: N_SurfMesh
USE MOD_HDG_Vars  ,ONLY: HDG_Surf_N,UseNSideMin
USE MOD_DG_Vars   ,ONLY: DG_Elems_master,DG_Elems_slave
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)          :: SendID
INTEGER, INTENT(IN)          :: mode
INTEGER, INTENT(IN), OPTIONAL:: ivar
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(MPI_Request), INTENT(OUT)         :: MPIRequest(nNbProcs)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: i,p,q,r,iSide,Nloc
!==================================================================================================================================
DO iNbProc=1,nNbProcs


  IF(nMPISides_send(iNbProc,SendID).GT.0)THEN
    ! Exchange surface elments on Nmax
    IF (mode.EQ.1) THEN
      nSendVal     = DataSizeSurfSendMax(iNbProc,SendID)
      SideID_start = OffsetMPISides_send(iNbProc-1,SendID)+1
      SideID_end   = OffsetMPISides_send(iNbProc,SendID)
      i = 1
      DO iSide = SideID_start, SideID_end
        Nloc = MAX(DG_Elems_master(iSide),DG_Elems_slave(iSide))
        DO p = 0, Nloc
          DO q = 0, Nloc
            SurfExchange(iNbProc)%SurfDataSend(i) = N_SurfMesh(iSide)%SurfElem(p,q)
            i = i + 1
          END DO ! q = 0, N_slave
        END DO ! p = 0, N_slave
      END DO ! iSide = SideID_start, SideID_end
      CALL MPI_ISEND(SurfExchange(iNbProc)%SurfDataSend(1:nSendVal),nSendVal,MPI_DOUBLE_PRECISION,  &
                      nbProc(iNbProc),0,MPI_COMM_PICLAS,MPIRequest(iNbProc),iError)
    ELSEIF(mode.EQ.8)THEN
      ! Exchange HDG solver variables on Nmin
      IF(UseNSideMin)THEN
        nSendVal     = DataSizeSurfSendMin(iNbProc,SendID)**2
      ELSE
        nSendVal     = DataSizeSurfSendMax(iNbProc,SendID)**2
      END IF
      SideID_start = OffsetMPISides_send(iNbProc-1,SendID)+1
      SideID_end   = OffsetMPISides_send(iNbProc,SendID)
      i = 1
      DO iSide = SideID_start, SideID_end
        IF(UseNSideMin)THEN
          Nloc = MIN(DG_Elems_master(iSide),DG_Elems_slave(iSide))
        ELSE
          Nloc = MAX(DG_Elems_master(iSide),DG_Elems_slave(iSide))
        END IF
        DO p = 1, (Nloc+1)**2
          DO q = 1, (Nloc+1)**2
            SurfExchange(iNbProc)%SurfDataSend2(i) = HDG_Surf_N(iSide)%Precond(p,q)
            i = i + 1
          END DO ! q = 0, N_slave
        END DO ! p = 0, N_slave
      END DO ! iSide = SideID_start, SideID_end
      CALL MPI_ISEND(SurfExchange(iNbProc)%SurfDataSend2(1:nSendVal),nSendVal,MPI_DOUBLE_PRECISION,  &
                      nbProc(iNbProc),0,MPI_COMM_PICLAS,MPIRequest(iNbProc),iError)
    ELSE
      ! Exchange HDG solver variables on Nmin
      IF(UseNSideMin)THEN
        nSendVal     = DataSizeSurfSendMin(iNbProc,SendID)
      ELSE
        nSendVal     = DataSizeSurfSendMax(iNbProc,SendID)
      END IF
      SideID_start = OffsetMPISides_send(iNbProc-1,SendID)+1
      SideID_end   = OffsetMPISides_send(iNbProc,SendID)
      i = 1
      DO iSide = SideID_start, SideID_end
        IF(UseNSideMin)THEN
          Nloc = MIN(DG_Elems_master(iSide),DG_Elems_slave(iSide))
        ELSE
          Nloc = MAX(DG_Elems_master(iSide),DG_Elems_slave(iSide))
        END IF
        DO p = 0, Nloc
          DO q = 0, Nloc
            r=p*(Nloc+1) + q+1
            IF (mode.EQ.2) THEN
              SurfExchange(iNbProc)%SurfDataSend(i) = HDG_Surf_N(iSide)%lambda(iVar,r)
            ELSEIF (mode.EQ.3) THEN
              SurfExchange(iNbProc)%SurfDataSend(i) = HDG_Surf_N(iSide)%mv(iVar,r)
            ELSEIF (mode.EQ.4) THEN
              SurfExchange(iNbProc)%SurfDataSend(i) = HDG_Surf_N(iSide)%RHS_face(iVar,r)
            ELSEIF (mode.EQ.5) THEN
              SurfExchange(iNbProc)%SurfDataSend(i) = HDG_Surf_N(iSide)%V(iVar,r)
            ELSEIF (mode.EQ.6) THEN
              SurfExchange(iNbProc)%SurfDataSend(i) = HDG_Surf_N(iSide)%Z(iVar,r)
            ELSEIF (mode.EQ.7) THEN
              SurfExchange(iNbProc)%SurfDataSend(i) = HDG_Surf_N(iSide)%InvPrecondDiag(r)
            ENDIF
            i = i + 1
          END DO ! q = 0, N_slave
        END DO ! p = 0, N_slave
      END DO ! iSide = SideID_start, SideID_end
      CALL MPI_ISEND(SurfExchange(iNbProc)%SurfDataSend(1:nSendVal),nSendVal,MPI_DOUBLE_PRECISION,  &
                      nbProc(iNbProc),0,MPI_COMM_PICLAS,MPIRequest(iNbProc),iError)
    END IF
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error in MPI_ISEND',iError)
  ELSE
    MPIRequest(iNbProc)=MPI_REQUEST_NULL
  END IF
END DO !iProc=1,nNBProcs
END SUBROUTINE StartSendMPISurfDataType


!===================================================================================================================================
!> We have to complete our non-blocking communication operations before we can (re)use the send / receive buffers
!> SendRequest, RecRequest: communication handles
!> SendID: defines the send / receive direction -> 1=send MINE / receive YOUR  2=send YOUR / receive MINE
!===================================================================================================================================
SUBROUTINE FinishExchangeMPISurfDataType(SendRequest,RecRequest,SendID, mode, iVar)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
USE MOD_Mesh_Vars ,ONLY: N_SurfMesh
USE MOD_HDG_Vars  ,ONLY: HDG_Surf_N,UseNSideMin
USE MOD_DG_Vars   ,ONLY: DG_Elems_master,DG_Elems_slave
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)          :: SendID, mode
INTEGER, INTENT(IN), OPTIONAL:: iVar
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(MPI_Request), INTENT(INOUT)       :: SendRequest(nNbProcs),RecRequest(nNbProcs)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=8)               :: CounterStart,CounterEnd
REAL(KIND=8)                  :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
INTEGER                       :: i,p,q,r,iSide,Nloc
!===================================================================================================================================
! Check receive operations first
DO iNbProc=1,nNbProcs
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    CALL MPI_WAIT(RecRequest(iNbProc) ,MPI_STATUS_IGNORE,iError)
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error iyyn MPI_WAIT',iError)
  END IF
END DO !iProc=1,nNBProcs

! Check send operations
DO iNbProc=1,nNbProcs
  IF(nMPISides_send(iNbProc,SendID).GT.0)THEN
    CALL MPI_WAIT(SendRequest(iNbProc),MPI_STATUS_IGNORE,iError)
    IF(iError.NE.MPI_SUCCESS) CALL Abort(__STAMP__,'Error iyyn MPI_WAIT',iError)
  END IF
END DO !iProc=1,nNBProcs

! Unroll data
DO iNbProc=1,nNbProcs
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    SideID_start = OffsetMPISides_rec(iNbProc-1,SendID)+1
    SideID_end   = OffsetMPISides_rec(iNbProc,SendID)

    i = 1
    IF(mode.EQ.1)THEN
      DO iSide = SideID_start, SideID_end
        Nloc = MAX(DG_Elems_master(iSide),DG_Elems_slave(iSide))
        DO p = 0, Nloc
          DO q = 0, Nloc
            N_SurfMesh(iSide)%SurfElem(p,q) = SurfExchange(iNbProc)%SurfDataRecv(i)
            i = i + 1
          END DO ! q = 0, N_slave
        END DO ! p = 0, N_slave
      END DO ! iSide = SideID_start, SideID_end
    ELSEIF(mode.EQ.8)THEN
      DO iSide = SideID_start, SideID_end
        IF(UseNSideMin)THEN
          Nloc = MIN(DG_Elems_master(iSide),DG_Elems_slave(iSide))
        ELSE
          Nloc = MAX(DG_Elems_master(iSide),DG_Elems_slave(iSide))
        END IF
        DO p = 1, (Nloc+1)**2
          DO q = 1, (Nloc+1)**2
            HDG_Surf_N(iSide)%Precond(p,q) = SurfExchange(iNbProc)%SurfDataRecv2(i)
            i = i + 1
          END DO ! q = 0, N_slave
        END DO ! p = 0, N_slave
      END DO ! iSide = SideID_start, SideID_end
    ELSE
      DO iSide = SideID_start, SideID_end
        IF(UseNSideMin)THEN
          Nloc = MIN(DG_Elems_master(iSide),DG_Elems_slave(iSide))
        ELSE
          Nloc = MAX(DG_Elems_master(iSide),DG_Elems_slave(iSide))
        END IF
        DO p = 0, Nloc
          DO q = 0, Nloc
            r=p*(Nloc+1) + q+1
            IF (mode.EQ.2) THEN
              HDG_Surf_N(iSide)%lambda(iVar,r)    = SurfExchange(iNbProc)%SurfDataRecv(i)
            ELSEIF (mode.EQ.3) THEN
              HDG_Surf_N(iSide)%mv(iVar,r)        = SurfExchange(iNbProc)%SurfDataRecv(i)
            ELSEIF (mode.EQ.4) THEN
              HDG_Surf_N(iSide)%RHS_face(iVar,r)  = SurfExchange(iNbProc)%SurfDataRecv(i)
            ELSEIF (mode.EQ.5) THEN
              HDG_Surf_N(iSide)%V(iVar,r)         = SurfExchange(iNbProc)%SurfDataRecv(i)
            ELSEIF (mode.EQ.6) THEN
              HDG_Surf_N(iSide)%Z(iVar,r)         = SurfExchange(iNbProc)%SurfDataRecv(i)
            ELSEIF (mode.EQ.7) THEN
              HDG_Surf_N(iSide)%InvPrecondDiag(r) = SurfExchange(iNbProc)%SurfDataRecv(i)
            ENDIF
            i = i + 1
          END DO ! q = 0, N_slave
        END DO ! p = 0, N_slave
      END DO ! iSide = SideID_start, SideID_end
    END IF ! SendID.EQ.2

  END IF
END DO !iProc=1,nNBProcs

END SUBROUTINE FinishExchangeMPISurfDataType
#endif /*USE_MPI && USE_HDG*/
END MODULE MOD_MPI_HDG
