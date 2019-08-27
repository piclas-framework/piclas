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

MODULE MOD_LD_mean_cell
!===================================================================================================================================
! module for determination of cell quantities
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE CalcMacCellLDValues
  MODULE PROCEDURE CalcMacCellLDValues
END INTERFACE

PUBLIC :: CalcMacCellLDValues
!===================================================================================================================================

CONTAINS

!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE CalcMacCellLDValues()
!===================================================================================================================================
! module for determination of cell quantities
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_Globals_Vars,           ONLY: BoltzmannConst
  USE MOD_LD_Vars
  USE MOD_Mesh_Vars,              ONLY : nElems
  USE MOD_Particle_Vars,          ONLY : PEM, usevMPF, PartMPF, Species, PartSpecies
  USE MOD_DSMC_Vars,              ONLY : SpecDSMC, CollisMode, LD_MultiTemperaturMod
  USE MOD_LD_internal_Temp
  USE MOD_Particle_Mesh_Vars,     ONLY : GEO
!#if (PP_TimeDiscMethod==1001)
!  USE MOD_LD_DSMC_TOOLS,          ONLY : LD_DSMC_Mean_Bufferzone_A_Val
!#endif
#if (PP_TimeDiscMethod!=1001)
  USE MOD_Mesh_Vars,              ONLY : ElemToSide, SideToElem
#endif
!--------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
   IMPLICIT NONE                                                                                      !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
  INTEGER             :: iElem, iPart, nPart, iPartIndx, LowPartCount
  REAL                :: MPFSum, WeightFak
  REAL                :: CellMass, CellTemp, CellTempMean, CellVelo2
  REAL                :: MeanRefTemp, MeanRefDiameter, MeanOmega
  REAL, PARAMETER     :: PI=3.14159265358979323846_8
#if (PP_TimeDiscMethod!=1001)
  INTEGER             :: NeiElem, nNei
  INTEGER             :: iLocSide, SideID
#endif
!--------------------------------------------------------------------------------------------------!

LowPartCount = 0
DO iElem = 1, nElems ! element=cell main loop
#if (PP_TimeDiscMethod==1001)
IF((BulkValues(iElem)%CellType.EQ.3).OR.(BulkValues(iElem)%CellType.EQ.4)) THEN  ! --- LD Cell ?
#endif
  nPart = PEM%pNumber(iElem)
  IF (nPart.GT. 1) THEN ! Are there more than one particle
    BulkValues(iElem)%CellV           = 0.0
    BulkValues(iElem)%MassDens        = 0.0
    BulkValues(iElem)%DegreeOfFreedom = 0.0
    CellTempMean                      = 0.0
    CellMass                          = 0.0
    MPFSum                            = 0.0
    CellVelo2                         = 0.0
    MeanRefTemp                       = 0.0
    MeanRefDiameter                   = 0.0
    MeanOmega                         = 0.0
    IF (CollisMode.GT.1) THEN
      IF(LD_MultiTemperaturMod.EQ. 1) THEN
        CALL CalcInternalTemp_LD_first(iElem)
      END IF
    END IF
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
!      IF (CellTempMean.lt. 0) THEN
!          SWRITE(UNIT_stdOut,'(A)') 'Element, Temperatur:',iElem, CellTempMean
!          CALL abort(__STAMP__&
!               'ERROR: Temperature is lt zero')
!      END IF
      BulkValues(iElem)%DegreeOfFreedom = BulkValues(iElem)%DegreeOfFreedom + PartStateBulkValues(iPartIndx,5) * WeightFak
      CellMass                          = CellMass + WeightFak * Species(PartSpecies(iPartIndx))%MassIC
!--- for viscousity terms...
      MeanRefTemp = MeanRefTemp + SpecDSMC(PartSpecies(iPartIndx))%TrefVHS * WeightFak
      MeanRefDiameter = MeanRefDiameter + SpecDSMC(PartSpecies(iPartIndx))%DrefVHS * WeightFak
      MeanOmega = MeanOmega + SpecDSMC(PartSpecies(iPartIndx))%omegaVHS * WeightFak
!--- end for viscousity terms
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
!--- for viscousity terms...
    MeanRefTemp = MeanRefTemp / MPFSum
    MeanRefDiameter = MeanRefDiameter / MPFSum
    MeanOmega = MeanOmega / MPFSum
    BulkValues(iElem)%DynamicVisc = 15 * SQRT(PI * (CellMass / MPFSum) * BoltzmannConst * MeanRefTemp) &
                                  / ( 8 * PI * MeanRefDiameter**2 * (2 - MeanOmega) * (3 - MeanOmega) ) &
                                  * ( CellTemp / MeanRefTemp)**(MeanOmega + 0.5)
    BulkValues(iElem)%ThermalCond = 0.25 * (9 + 2 * BulkValues(iElem)%DegreeOfFreedom) &
                                  * BulkValues(iElem)%DynamicVisc &
                                  * BoltzmannConst / (CellMass / MPFSum)
    BulkValues(iElem)%SpezGasConst = BoltzmannConst / (CellMass / MPFSum)
!--- end for viscousity terms
  ELSE ! if there is no particle => keep old LD-state
    LowPartCount = LowPartCount + 1
  END IF
#if (PP_TimeDiscMethod==1001)
ELSE IF(BulkValues(iElem)%CellType.EQ.2) THEN  ! --- Bufferzone_A Cell ?
  CALL LD_DSMC_Mean_Bufferzone_A_Val(iElem)
END IF  ! --- END LD Cell or Bufferzone_A?
#endif
END DO
#if (PP_TimeDiscMethod!=1001)
! Abfangen von Nulldivisionen
IF (LowPartCount.GE. 1) THEN
    SWRITE(UNIT_StdOut,'(A)') 'ATTENTION: YOU NEED MORE PARTCLES FOR LD!!!'
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
          CALL abort(&
          __STAMP__&
          ,'ERROR: YOU NEED MORE PARTCLES FOR LD!!!')
      END IF
      BulkValues(iElem)%CellV(1:3) = BulkValues(iElem)%CellV(1:3) / nNei
      BulkValues(iElem)%DegreeOfFreedom = BulkValues(iElem)%DegreeOfFreedom / nNei
      BulkValues(iElem)%Beta = BulkValues(iElem)%Beta / nNei
    END IF
  END DO
END IF
#endif

#if USE_MPI
  CALL LD_MPI_Communication
#endif

END SUBROUTINE CalcMacCellLDValues

!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
#if USE_MPI
SUBROUTINE LD_MPI_Communication
!===================================================================================================================================
! MPI communication of LD values
!===================================================================================================================================
! MODULES
USE MOD_LD_Vars
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Vars
USE MOD_Particle_MPI_Vars,      ONLY : PartMPI
USE MOD_Mesh_Vars,              ONLY : SideToElem
!--------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
  REAL  , ALLOCATABLE :: MPI_LDSendBuffer(:)  ! Buffer for MIPSend (nMPISide*5)
  REAL  , ALLOCATABLE :: MPI_LDRecvBuffer(:)  ! Buffer for MIPRecv (nMPISide*5)
  INTEGER              :: iProc, LD_MPISideID_start, LD_MPISideID_end, iSideID, NbrOfSendLDVal
  INTEGER              :: RecvBufferCount, SendBufferCount, SideToElem_ID
!-----------------------------------------------------------------------------------------------------------------------------------

  DO iProc=1,nNbProcs
    NbrOfSendLDVal = 8*(nMPISides_MINE_Proc(iProc) + nMPISides_YOUR_Proc(iProc))
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
        MPI_LDSendBuffer(8*SendBufferCount + 1) = BulkValues(SideToElem(SideToElem_ID,iSideID))%CellV(1)
        MPI_LDSendBuffer(8*SendBufferCount + 2) = BulkValues(SideToElem(SideToElem_ID,iSideID))%CellV(2)
        MPI_LDSendBuffer(8*SendBufferCount + 3) = BulkValues(SideToElem(SideToElem_ID,iSideID))%CellV(3)
        MPI_LDSendBuffer(8*SendBufferCount + 4) = BulkValues(SideToElem(SideToElem_ID,iSideID))%Beta
        MPI_LDSendBuffer(8*SendBufferCount + 5) = BulkValues(SideToElem(SideToElem_ID,iSideID))%MassDens
        MPI_LDSendBuffer(8*SendBufferCount + 6) = BulkValues(SideToElem(SideToElem_ID,iSideID))%DynamicVisc
        MPI_LDSendBuffer(8*SendBufferCount + 7) = BulkValues(SideToElem(SideToElem_ID,iSideID))%ThermalCond
        MPI_LDSendBuffer(8*SendBufferCount + 8) = BulkValues(SideToElem(SideToElem_ID,iSideID))%BulkTemperature
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
        MPI_LDSendBuffer(8*SendBufferCount + 1) = BulkValues(SideToElem(SideToElem_ID,iSideID))%CellV(1)
        MPI_LDSendBuffer(8*SendBufferCount + 2) = BulkValues(SideToElem(SideToElem_ID,iSideID))%CellV(2)
        MPI_LDSendBuffer(8*SendBufferCount + 3) = BulkValues(SideToElem(SideToElem_ID,iSideID))%CellV(3)
        MPI_LDSendBuffer(8*SendBufferCount + 4) = BulkValues(SideToElem(SideToElem_ID,iSideID))%Beta
        MPI_LDSendBuffer(8*SendBufferCount + 5) = BulkValues(SideToElem(SideToElem_ID,iSideID))%MassDens
        MPI_LDSendBuffer(8*SendBufferCount + 6) = BulkValues(SideToElem(SideToElem_ID,iSideID))%DynamicVisc
        MPI_LDSendBuffer(8*SendBufferCount + 7) = BulkValues(SideToElem(SideToElem_ID,iSideID))%ThermalCond
        MPI_LDSendBuffer(8*SendBufferCount + 8) = BulkValues(SideToElem(SideToElem_ID,iSideID))%BulkTemperature
        SendBufferCount = SendBufferCount + 1
      END DO
    END IF
    IF (PartMPI%MyRank.LT.nbProc(iProc)) THEN
      CALL MPI_SEND(MPI_LDSendBuffer,NbrOfSendLDVal,MPI_DOUBLE_PRECISION,nbProc(iProc),1101,PartMPI%COMM,IERROR)
      CALL MPI_RECV(MPI_LDRecvBuffer,NbrOfSendLDVal,MPI_DOUBLE_PRECISION,nbProc(iProc),1101,PartMPI%COMM,MPISTATUS,IERROR)
    ELSE IF (PartMPI%MyRank.GT.nbProc(iProc)) THEN
      CALL MPI_RECV(MPI_LDRecvBuffer,NbrOfSendLDVal,MPI_DOUBLE_PRECISION,nbProc(iProc),1101,PartMPI%COMM,MPISTATUS,IERROR)
      CALL MPI_SEND(MPI_LDSendBuffer,NbrOfSendLDVal,MPI_DOUBLE_PRECISION,nbProc(iProc),1101,PartMPI%COMM,IERROR)
    END IF

    RecvBufferCount = 0
    IF(nMPISides_YOUR_Proc(iProc).GT.0)THEN
      LD_MPISideID_start = OffsetMPISides_YOUR(iProc-1) + 1 !- OffsetMPISides_MINE(0) ! OffsetMPISides_MINE(0) == inner + BC sides
      LD_MPISideID_end   = OffsetMPISides_YOUR(iProc) !- OffsetMPISides_MINE(0)
      DO iSideID = LD_MPISideID_start, LD_MPISideID_end
        MPINeighborBulkVal(iSideID,1) = MPI_LDRecvBuffer(8*RecvBufferCount + 1)
        MPINeighborBulkVal(iSideID,2) = MPI_LDRecvBuffer(8*RecvBufferCount + 2)
        MPINeighborBulkVal(iSideID,3) = MPI_LDRecvBuffer(8*RecvBufferCount + 3)
        MPINeighborBulkVal(iSideID,4) = MPI_LDRecvBuffer(8*RecvBufferCount + 4)
        MPINeighborBulkVal(iSideID,5) = MPI_LDRecvBuffer(8*RecvBufferCount + 5)
        MPINeighborBulkVal(iSideID,6) = MPI_LDRecvBuffer(8*RecvBufferCount + 6)
        MPINeighborBulkVal(iSideID,7) = MPI_LDRecvBuffer(8*RecvBufferCount + 7)
        MPINeighborBulkVal(iSideID,8) = MPI_LDRecvBuffer(8*RecvBufferCount + 8)
        RecvBufferCount = RecvBufferCount + 1
      END DO
    END IF
    IF(nMPISides_MINE_Proc(iProc).GT.0)THEN
      LD_MPISideID_start = OffsetMPISides_MINE(iProc-1) + 1 !- OffsetMPISides_MINE(0) ! OffsetMPISides_MINE(0) == inner + BC sides
      LD_MPISideID_end   = OffsetMPISides_MINE(iProc) !- OffsetMPISides_MINE(0)
      DO iSideID = LD_MPISideID_start, LD_MPISideID_end
        MPINeighborBulkVal(iSideID,1) = MPI_LDRecvBuffer(8*RecvBufferCount + 1)
        MPINeighborBulkVal(iSideID,2) = MPI_LDRecvBuffer(8*RecvBufferCount + 2)
        MPINeighborBulkVal(iSideID,3) = MPI_LDRecvBuffer(8*RecvBufferCount + 3)
        MPINeighborBulkVal(iSideID,4) = MPI_LDRecvBuffer(8*RecvBufferCount + 4)
        MPINeighborBulkVal(iSideID,5) = MPI_LDRecvBuffer(8*RecvBufferCount + 5)
        MPINeighborBulkVal(iSideID,6) = MPI_LDRecvBuffer(8*RecvBufferCount + 6)
        MPINeighborBulkVal(iSideID,7) = MPI_LDRecvBuffer(8*RecvBufferCount + 7)
        MPINeighborBulkVal(iSideID,8) = MPI_LDRecvBuffer(8*RecvBufferCount + 8)
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
