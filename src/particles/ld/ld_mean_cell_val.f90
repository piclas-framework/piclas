MODULE MOD_LD_mean_cell

!===================================================================================================================================
! module for determination of cell quantities
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

PUBLIC :: CalcMacCellLDValues
!===================================================================================================================================

CONTAINS

!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE CalcMacCellLDValues()

USE MOD_LD_Vars
USE MOD_Mesh_Vars,              ONLY : nElems, ElemToSide, SideToElem
USE MOD_Particle_Vars,          ONLY : GEO, PEM, usevMPF, PartMPF, BoltzmannConst, Species, PartSpecies, PDM

!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                  !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
  INTEGER             :: iElem, iPart, nPart, iPartIndx, LowPartCount
  REAL                :: MPFSum, WeightFak
  REAL                :: CellMass, CellTemp, CellTempMean, CellVelo2, StanVol
  INTEGER             :: NeiElem, nNei
  INTEGER             :: iLocSide, SideID
!--------------------------------------------------------------------------------------------------!

LowPartCount = 0
DO iElem = 1, nElems ! element=cell main loop
  nPart = PEM%pNumber(iElem)
  IF (nPart.GT. 1) THEN ! Are there more than one particle
    BulkValues(iElem)%CellV           = 0.0
    BulkValues(iElem)%MassDens        = 0.0
    BulkValues(iElem)%DegreeOfFreedom = 0.0
    CellTempMean                      = 0.0
    CellMass                          = 0.0
    MPFSum                            = 0.0
    CellVelo2                         = 0.0
    iPartIndx = PEM%pStart(iElem)
    DO ipart = 1, nPart
      IF (usevMPF) THEN
         WeightFak = PartMPF(iPartIndx)
      ELSE
         WeightFak = Species(PartSpecies(iPartIndx))%MacroParticleFactor
      END IF
      BulkValues(iElem)%CellV(1)        = BulkValues(iElem)%CellV(1) + PartStateBulkValues(iPartIndx,1) &
                                        * WeightFak * Species(PartSpecies(iPartIndx))%MassIC
      BulkValues(iElem)%CellV(2)        = BulkValues(iElem)%CellV(2) + PartStateBulkValues(iPartIndx,2) &
                                        * WeightFak * Species(PartSpecies(iPartIndx))%MassIC
      BulkValues(iElem)%CellV(3)        = BulkValues(iElem)%CellV(3) + PartStateBulkValues(iPartIndx,3) &
                                        * WeightFak * Species(PartSpecies(iPartIndx))%MassIC
      CellVelo2                         = CellVelo2 + ( (PartStateBulkValues(iPartIndx,1))**2 &
                                                      + (PartStateBulkValues(iPartIndx,2))**2 &
                                                      + (PartStateBulkValues(iPartIndx,3))**2 ) &
                                                      * WeightFak * Species(PartSpecies(iPartIndx))%MassIC
      CellTempMean                      = CellTempMean + PartStateBulkValues(iPartIndx,4) * WeightFak
      IF (CellTempMean.lt. 0) THEN
        PRINT*,iElem, iPartIndx, PDM%ParticleVecLength, PartStateBulkValues(iPartIndx,4)
        read*
      END IF
      BulkValues(iElem)%DegreeOfFreedom = BulkValues(iElem)%DegreeOfFreedom + PartStateBulkValues(iPartIndx,5) * WeightFak
      CellMass                          = CellMass + WeightFak * Species(PartSpecies(iPartIndx))%MassIC
      MPFSum = MPFSum + WeightFak
      iPartIndx = PEM%pNext(iPartIndx)
    END DO

    BulkValues(iElem)%CellV(1) = BulkValues(iElem)%CellV(1) / CellMass
    BulkValues(iElem)%CellV(2) = BulkValues(iElem)%CellV(2) / CellMass
    BulkValues(iElem)%CellV(3) = BulkValues(iElem)%CellV(3) / CellMass
    CellVelo2                  = CellVelo2 / CellMass
    CellTempMean               = CellTempMean / MPFSum
    BulkValues(iElem)%DegreeOfFreedom = BulkValues(iElem)%DegreeOfFreedom / MPFSum
    BulkValues(iElem)%MassDens = CellMass / GEO%Volume(iElem)
    CellTemp = CellTempMean &
             + CellMass / MPFSum /(BulkValues(iElem)%DegreeOfFreedom * BoltzmannConst) & 
             * (nPart/(nPart-1)) &
             * (CellVelo2 - ( BulkValues(iElem)%CellV(1)**2 &
                            + BulkValues(iElem)%CellV(2)**2 &
                            + BulkValues(iElem)%CellV(3)**2 ))
    BulkValues(iElem)%BulkTemperature = CellTemp
    BulkValues(iElem)%Beta = SQRT(CellMass / MPFSum / (2 * CellTemp * BoltzmannConst))
  ELSE ! if there is no particle => keep old LD-state
    LowPartCount = LowPartCount + 1
  END IF
END DO
! Abfangen von Nulldivisionen
IF (LowPartCount.GE. 1) THEN
  DO iElem = 1, nElems
    nPart = PEM%pNumber(iElem)
    IF (nPart.LE. 1) THEN
      BulkValues(iElem)%CellV           = 0.0
      BulkValues(iElem)%DegreeOfFreedom = 0.0
      BulkValues(iElem)%Beta            = 0.0
      BulkValues(iElem)%MassDens        = 0.0
      nNei = 0
      DO iLocSide=1, 6
        SideID = ElemToSide(1,iLocSide,iElem)
        IF (SideToElem(1,SideID) .EQ. iElem) THEN
          IF (SideToElem(2,SideID).LT. 0) CYCLE
          IF (PEM%pNumber(SideToElem(2,SideID)) .LE. 1) CYCLE
          NeiElem = SideToElem(2,SideID)
          BulkValues(iElem)%CellV(1:3) = BulkValues(iElem)%CellV(1:3) &
                                       + BulkValues(NeiElem)%CellV(1:3)
          BulkValues(iElem)%DegreeOfFreedom = BulkValues(iElem)%DegreeOfFreedom &
                                            + BulkValues(NeiElem)%DegreeOfFreedom
          BulkValues(iElem)%Beta = BulkValues(iElem)%Beta &
                                 + BulkValues(NeiElem)%Beta
          BulkValues(iElem)%MassDens = BulkValues(iElem)%MassDens &
                                     + BulkValues(NeiElem)%MassDens
          nNei = nNei + 1
        ELSE
          IF (SideToElem(1,SideID).LT. 0) CYCLE
          IF (PEM%pNumber(SideToElem(1,SideID)) .LE. 1) CYCLE
          NeiElem = SideToElem(1,SideID)
          BulkValues(iElem)%CellV(1:3) = BulkValues(iElem)%CellV(1:3) &
                                       + BulkValues(NeiElem)%CellV(1:3)
          BulkValues(iElem)%DegreeOfFreedom = BulkValues(iElem)%DegreeOfFreedom &
                                            + BulkValues(NeiElem)%DegreeOfFreedom
          BulkValues(iElem)%Beta = BulkValues(iElem)%Beta &
                                 + BulkValues(NeiElem)%Beta
          BulkValues(iElem)%MassDens = BulkValues(iElem)%MassDens &
                                     + BulkValues(NeiElem)%MassDens
          nNei = nNei + 1
        END IF
      END DO
      IF (nNei.EQ.0) THEN
        PRINT*,'YOU NEED MORE PARTCLES FOR LD!!!'
        PRINT*,iElem, nNei, BulkValues(iElem)%Beta
        STOP
      END IF
      BulkValues(iElem)%CellV(1:3) = BulkValues(iElem)%CellV(1:3) / nNei
      BulkValues(iElem)%DegreeOfFreedom = BulkValues(iElem)%DegreeOfFreedom / nNei
      BulkValues(iElem)%Beta = BulkValues(iElem)%Beta / nNei   
    END IF     
  END DO
END IF

#ifdef MPI
  CALL LD_MPI_Communication
#endif

END SUBROUTINE CalcMacCellLDValues

!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
#ifdef MPI
SUBROUTINE LD_MPI_Communication

USE MOD_LD_Vars
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
USE MOD_part_MPI_Vars, ONLY : PMPIVAR
USE MOD_Mesh_Vars,             ONLY : SideToElem

IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
  REAL  , ALLOCATABLE :: MPI_LDSendBuffer(:)  ! Buffer for MIPSend (nMPISide*5)
  REAL  , ALLOCATABLE :: MPI_LDRecvBuffer(:)  ! Buffer for MIPRecv (nMPISide*5)
  INTEGER              :: iProc, LD_MPISideID_start, LD_MPISideID_end, iSideID, NbrOfSendLDVal
  INTEGER              :: RecvBufferCount, SendBufferCount, SideToElem_ID
!-----------------------------------------------------------------------------------------------------------------------------------

  DO iProc=1,nNbProcs
    NbrOfSendLDVal = 5*(nMPISides_MINE_Proc(iProc) + nMPISides_YOUR_Proc(iProc))
    ALLOCATE(MPI_LDSendBuffer(NbrOfSendLDVal))
    ALLOCATE(MPI_LDRecvBuffer(NbrOfSendLDVal))
    SendBufferCount = 0
    IF(nMPISides_MINE_Proc(iProc).GT.0)THEN
      LD_MPISideID_start = OffsetMPISides_MINE(iProc-1) + 1 !- OffsetMPISides_MINE(0) ! OffsetMPISides_MINE(0) == inner + BC sides
      LD_MPISideID_end   = OffsetMPISides_MINE(iProc) !- OffsetMPISides_MINE(0)
      DO iSideID = LD_MPISideID_start, LD_MPISideID_end
        IF (SideToElem(1,iSideID).LT. 0) THEN
          SideToElem_ID = 2
        ELSE
          SideToElem_ID = 1
        END IF
        MPI_LDSendBuffer(5*SendBufferCount + 1) = BulkValues(SideToElem(SideToElem_ID,iSideID))%CellV(1)
        MPI_LDSendBuffer(5*SendBufferCount + 2) = BulkValues(SideToElem(SideToElem_ID,iSideID))%CellV(2)
        MPI_LDSendBuffer(5*SendBufferCount + 3) = BulkValues(SideToElem(SideToElem_ID,iSideID))%CellV(3)
        MPI_LDSendBuffer(5*SendBufferCount + 4) = BulkValues(SideToElem(SideToElem_ID,iSideID))%Beta
        MPI_LDSendBuffer(5*SendBufferCount + 5) = BulkValues(SideToElem(SideToElem_ID,iSideID))%MassDens
        SendBufferCount = SendBufferCount + 1
      END DO
    END IF
    IF(nMPISides_YOUR_Proc(iProc).GT.0)THEN
      LD_MPISideID_start = OffsetMPISides_YOUR(iProc-1) + 1 !- OffsetMPISides_MINE(0) ! OffsetMPISides_MINE(0) == inner + BC sides
      LD_MPISideID_end   = OffsetMPISides_YOUR(iProc) !- OffsetMPISides_MINE(0)
      DO iSideID = LD_MPISideID_start, LD_MPISideID_end
        IF (SideToElem(1,iSideID).LT. 0) THEN
          SideToElem_ID = 2
        ELSE
          SideToElem_ID = 1
        END IF
        MPI_LDSendBuffer(5*SendBufferCount + 1) = BulkValues(SideToElem(SideToElem_ID,iSideID))%CellV(1)
        MPI_LDSendBuffer(5*SendBufferCount + 2) = BulkValues(SideToElem(SideToElem_ID,iSideID))%CellV(2)
        MPI_LDSendBuffer(5*SendBufferCount + 3) = BulkValues(SideToElem(SideToElem_ID,iSideID))%CellV(3)
        MPI_LDSendBuffer(5*SendBufferCount + 4) = BulkValues(SideToElem(SideToElem_ID,iSideID))%Beta
        MPI_LDSendBuffer(5*SendBufferCount + 5) = BulkValues(SideToElem(SideToElem_ID,iSideID))%MassDens
        SendBufferCount = SendBufferCount + 1
      END DO
    END IF
    IF (PMPIVAR%iProc.LT.nbProc(iProc)) THEN
      CALL MPI_SEND(MPI_LDSendBuffer,NbrOfSendLDVal,MPI_DOUBLE_PRECISION,nbProc(iProc),1101,PMPIVAR%COMM,IERROR)
      CALL MPI_RECV(MPI_LDRecvBuffer,NbrOfSendLDVal,MPI_DOUBLE_PRECISION,nbProc(iProc),1101,PMPIVAR%COMM,MPISTATUS,IERROR)
    ELSE IF (PMPIVAR%iProc.GT.nbProc(iProc)) THEN
      CALL MPI_RECV(MPI_LDRecvBuffer,NbrOfSendLDVal,MPI_DOUBLE_PRECISION,nbProc(iProc),1101,PMPIVAR%COMM,MPISTATUS,IERROR)
      CALL MPI_SEND(MPI_LDSendBuffer,NbrOfSendLDVal,MPI_DOUBLE_PRECISION,nbProc(iProc),1101,PMPIVAR%COMM,IERROR)
    END IF

    RecvBufferCount = 0
    IF(nMPISides_YOUR_Proc(iProc).GT.0)THEN
      LD_MPISideID_start = OffsetMPISides_YOUR(iProc-1) + 1 !- OffsetMPISides_MINE(0) ! OffsetMPISides_MINE(0) == inner + BC sides
      LD_MPISideID_end   = OffsetMPISides_YOUR(iProc) !- OffsetMPISides_MINE(0)
      DO iSideID = LD_MPISideID_start, LD_MPISideID_end
        MPINeighborBulkVal(iSideID,1) = MPI_LDRecvBuffer(5*RecvBufferCount + 1)
        MPINeighborBulkVal(iSideID,2) = MPI_LDRecvBuffer(5*RecvBufferCount + 2)
        MPINeighborBulkVal(iSideID,3) = MPI_LDRecvBuffer(5*RecvBufferCount + 3)
        MPINeighborBulkVal(iSideID,4) = MPI_LDRecvBuffer(5*RecvBufferCount + 4)
        MPINeighborBulkVal(iSideID,5) = MPI_LDRecvBuffer(5*RecvBufferCount + 5)
        RecvBufferCount = RecvBufferCount + 1
      END DO
    END IF
    IF(nMPISides_MINE_Proc(iProc).GT.0)THEN
      LD_MPISideID_start = OffsetMPISides_MINE(iProc-1) + 1 !- OffsetMPISides_MINE(0) ! OffsetMPISides_MINE(0) == inner + BC sides
      LD_MPISideID_end   = OffsetMPISides_MINE(iProc) !- OffsetMPISides_MINE(0)
      DO iSideID = LD_MPISideID_start, LD_MPISideID_end
        MPINeighborBulkVal(iSideID,1) = MPI_LDRecvBuffer(5*RecvBufferCount + 1)
        MPINeighborBulkVal(iSideID,2) = MPI_LDRecvBuffer(5*RecvBufferCount + 2)
        MPINeighborBulkVal(iSideID,3) = MPI_LDRecvBuffer(5*RecvBufferCount + 3)
        MPINeighborBulkVal(iSideID,4) = MPI_LDRecvBuffer(5*RecvBufferCount + 4)
        MPINeighborBulkVal(iSideID,5) = MPI_LDRecvBuffer(5*RecvBufferCount + 5)
        RecvBufferCount = RecvBufferCount + 1 
      END DO
    END IF

    DEALLOCATE(MPI_LDSendBuffer)
    DEALLOCATE(MPI_LDRecvBuffer)

  END DO !iProc=1,nNBProcs



END SUBROUTINE LD_MPI_Communication
#endif
!-----------------------------------------------------------------------------------------------------------------------------------

END MODULE MOD_LD_mean_cell
