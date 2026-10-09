
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

!===================================================================================================================================
!> Module for the HDG method
!===================================================================================================================================
MODULE MOD_HDG_Readin
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
#if USE_HDG
PUBLIC :: ReadFPCDataFromH5
PUBLIC :: ReadEPCDataFromH5
#if USE_MPI
PUBLIC :: SynchronizeChargeOnFPC
PUBLIC :: SynchronizeVoltageOnEPC
#if defined(PARTICLES)
PUBLIC :: SynchronizeBV
#endif /*defined(PARTICLES)*/
#endif /*USE_MPI */
#if defined(PARTICLES)
PUBLIC :: ReadBVDataFromH5
#if USE_PETSC
PUBLIC :: ReadCMBCDataFromH5
PUBLIC :: UpdateChargeOnCMBC
#if USE_MPI && defined(CODE_ANALYZE)
PUBLIC :: BroadcastChargeOnCMBC
#endif /*USE_MPI && defined(CODE_ANALYZE)*/
#endif /*USE_PETSC*/
#endif /*defined(PARTICLES)*/
#endif /*USE_HDG*/
!===================================================================================================================================

CONTAINS

#if USE_HDG
#if defined(PARTICLES)
!===================================================================================================================================
!> Read the bias voltage (BV) data from a .h5 state file.
!> 1. The MPI root process reads the info and checks data consistency
!> 2. The MPI root process distributes the information among the sub-communicator processes connected to the BV boundary.
!===================================================================================================================================
SUBROUTINE ReadBVDataFromH5()
! MODULES
USE MOD_io_hdf5
USE MOD_Globals          ,ONLY: UNIT_stdOut,MPIRoot,IK,abort
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
#endif /*USE_LOADBALANCE*/
USE MOD_IO_HDF5          ,ONLY: OpenDataFile,CloseDataFile,File_ID
USE MOD_Restart_Vars     ,ONLY: DoRestart,RestartFile
USE MOD_HDF5_Input       ,ONLY: DatasetExists,ReadArray,GetDataSize
USE MOD_HDG_Vars         ,ONLY: BiasVoltage,BVDataLength
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(255) :: ContainerName
LOGICAL        :: BVExists
REAL           :: BVDataHDF5(1:BVDataLength)
!===================================================================================================================================
! Only required during restart
IF(.NOT.DoRestart) RETURN

#if USE_LOADBALANCE
! Do not try to read the data from .h5 if load balance is performed without creating a .h5 restart file
IF(PerformLoadBalance.AND..NOT.(UseH5IOLoadBalance)) RETURN
#endif /*USE_LOADBALANCE*/

! 1. The MPI root process reads the info and checks data consistency
! Only root reads the values and distributes them via MPI Broadcast
IF(MPIRoot)THEN
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
  ! Check old parameter name
  ContainerName='BiasVoltage'
  CALL DatasetExists(File_ID,TRIM(ContainerName),BVExists)
  ! Check for new parameter name
  IF(BVExists)THEN
    CALL ReadArray(TRIM(ContainerName) , 2 , (/1_IK , INT(BVDataLength,IK)/) , 0_IK , 1 , RealArray=BVDataHDF5)
    WRITE(UNIT_stdOut,'(3(A,ES10.2E3))') " Read bias voltage from restart file ["//TRIM(RestartFile)//&
        "] Bias voltage [V]: ",BVDataHDF5(1),", Ion excess [C]: ",BVDataHDF5(2),", next adjustment time [s]: ",BVDataHDF5(3)
    BiasVoltage%BVData = BVDataHDF5
  END IF ! BVExists
  CALL CloseDataFile()
END IF ! MPIRoot

#if USE_MPI
! 2. The MPI root process distributes the information among the sub-communicator processes for each EPC
CALL SynchronizeBV()
#endif /*USE_MPI*/
END SUBROUTINE ReadBVDataFromH5


#if USE_MPI
!===================================================================================================================================
!> Communicate the bias voltage values from MPIRoot to sub-communicator processes
!===================================================================================================================================
SUBROUTINE SynchronizeBV()
! MODULES
USE mpi_f08
USE MOD_Globals  ,ONLY: IERROR,MPI_COMM_NULL,MPI_DOUBLE_PRECISION
USE MOD_HDG_Vars ,ONLY: BiasVoltage,BVDataLength
! insert modules here
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(BiasVoltage%COMM%UNICATOR.NE.MPI_COMM_NULL)THEN
  ! Broadcast from root to other processors on the sub-communicator
  CALL MPI_BCAST(BiasVoltage%BVData, BVDataLength, MPI_DOUBLE_PRECISION, 0, BiasVoltage%COMM%UNICATOR, IERROR)
END IF
END SUBROUTINE SynchronizeBV
#endif /*USE_MPI*/


#if USE_PETSC
!===================================================================================================================================
!> Read the Circuit Model BC (CMBC) data from a .h5 state file.
!> 1. The MPI root process reads the info and checks data consistency
!> 2. The MPI root process distributes the information among the sub-communicator processes connected to the CM boundary.
!===================================================================================================================================
SUBROUTINE ReadCMBCDataFromH5()
! MODULES
USE MOD_io_hdf5
USE MOD_Globals          ,ONLY: UNIT_stdOut,MPIRoot,IK,abort
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
#endif /*USE_LOADBALANCE*/
USE MOD_IO_HDF5          ,ONLY: OpenDataFile,CloseDataFile,File_ID
USE MOD_Restart_Vars     ,ONLY: DoRestart,RestartFile
USE MOD_HDF5_Input       ,ONLY: DatasetExists,ReadArray,GetDataSize
USE MOD_HDG_Vars         ,ONLY: CMBC,CMBCDataLength
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(255) :: ContainerName
LOGICAL        :: CMBCExists
REAL           :: CMBCDataHDF5(1:CMBCDataLength)
!===================================================================================================================================
! Only required during restart
IF(.NOT.DoRestart) RETURN

#if USE_LOADBALANCE
! Do not try to read the data from .h5 if load balance is performed without creating a .h5 restart file
IF(PerformLoadBalance.AND..NOT.(UseH5IOLoadBalance)) RETURN
#endif /*USE_LOADBALANCE*/

! 1. The MPI root process reads the info and checks data consistency
! Only root reads the values and distributes them via MPI Broadcast
IF(MPIRoot)THEN
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
  ! Check old parameter name
  ContainerName='CMBC'
  CALL DatasetExists(File_ID,TRIM(ContainerName),CMBCExists)
  ! Check for new parameter name
  IF(CMBCExists)THEN
    CALL ReadArray(TRIM(ContainerName) , 2 , (/1_IK , INT(CMBCDataLength,IK)/) , 0_IK , 1 , RealArray=CMBCDataHDF5)
    WRITE(UNIT_stdOut,'(2(A,ES10.2E3))') " Read circuit model anode voltage and charge from restart file ["//TRIM(RestartFile)//&
        "] Anode voltage [V]: ",CMBCDataHDF5(1),", Anode charge [C]: ",CMBCDataHDF5(2)
    CMBC%Voltage = CMBCDataHDF5(1)
    CMBC%Charge  = CMBCDataHDF5(2)
  END IF ! CMBCExists
  CALL CloseDataFile()
END IF ! MPIRoot

END SUBROUTINE ReadCMBCDataFromH5


!===================================================================================================================================
!> Communicate the Circuit Model accumulated charge values to MPIRoot: Updates CMBC%Charge and nullifies CMBC%ChargeProc
!===================================================================================================================================
SUBROUTINE UpdateChargeOnCMBC()
! MODULES
#if USE_MPI
USE mpi_f08
USE MOD_Globals  ,ONLY: IERROR,MPI_COMM_NULL,MPI_DOUBLE_PRECISION,MPIRoot
#endif /*USE_MPI*/
USE MOD_HDG_Vars ,ONLY: CMBC,CMBCDataLength
! insert modules here
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
#if USE_MPI
! Communicate the accumulated charged on each BC to MPIRoot
ASSOCIATE( COMM => CMBC%COMM%UNICATOR)
  IF(COMM.NE.MPI_COMM_NULL)THEN
    IF(MPIRoot)THEN
      CALL MPI_REDUCE(MPI_IN_PLACE, CMBC%ChargeProc, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM, IERROR)
      ! Update CMBC%Charge on MPIRoot
      CMBC%Charge = CMBC%Charge + CMBC%ChargeProc
    ELSE
      CALL MPI_REDUCE(CMBC%ChargeProc, 0           , 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM, IERROR)
      CMBC%Charge = 0. ! Non-Root processes always have zero
    END IF ! MPIRoot
  END IF ! COMM.NE.MPI_COMM_NULL
END ASSOCIATE
#else
! Update CMBC%Charge on MPIRoot
CMBC%Charge = CMBC%Charge + CMBC%ChargeProc
#endif /*USE_MPI*/
! Reset the coutner
CMBC%ChargeProc = 0.
END SUBROUTINE UpdateChargeOnCMBC


#if USE_MPI && defined(CODE_ANALYZE)
!===================================================================================================================================
!> Communicate the Circuit Model accumulated charge values from MPIRoot to ALL the other processes for analytical potential
!> calculation. This is only required for a subsequent L2 error calculation.
!===================================================================================================================================
SUBROUTINE BroadcastChargeOnCMBC()
! MODULES
USE mpi_f08
USE MOD_Globals       ,ONLY: IERROR,MPI_COMM_PICLAS,MPI_DOUBLE_PRECISION
USE MOD_HDG_Vars      ,ONLY: CMBC
USE MOD_Equation_Vars ,ONLY: IniExactFunc
! insert modules here
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(IniExactFunc.NE.10001) RETURN
! Communicate the accumulated charged on each BC to MPIRoot
CALL MPI_BCAST(CMBC%Charge,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_PICLAS,iError)
! Communicate the RF voltage
CALL MPI_BCAST(CMBC%VoltageRF(1),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_PICLAS,iError)
END SUBROUTINE BroadcastChargeOnCMBC
#endif /*USE_MPI && defined(CODE_ANALYZE)*/
#endif /*USE_PETSC*/
#endif /*defined(PARTICLES)*/


!===================================================================================================================================
!> Read the electric charge that resides on each FPC boundary from .h5 state file.
!> 1. The MPI root process reads the info and checks data consistency
!> 2. The MPI root process distributes the information among the sub-communicator processes for each FPC
!===================================================================================================================================
SUBROUTINE ReadFPCDataFromH5()
! MODULES
USE MOD_io_hdf5
USE MOD_Globals            ,ONLY: UNIT_stdOut,MPIRoot,IERROR,IK,abort
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
#endif /*USE_LOADBALANCE*/
USE MOD_IO_HDF5            ,ONLY: OpenDataFile,CloseDataFile,File_ID
USE MOD_Restart_Vars       ,ONLY: DoRestart,RestartFile
USE MOD_HDF5_Input         ,ONLY: DatasetExists,ReadArray,GetDataSize
USE MOD_HDG_Vars           ,ONLY: FPC
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
! Space-separated list of input and output types. Use: (int|real|logical|...)_(in|out|inout)_dim(n)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iUniqueFPCBC,nFPCHDF5,nDimsFPC
CHARACTER(255)     :: ContainerName
CHARACTER(1000)    :: TmpStr
LOGICAL            :: FPCExists
REAL,ALLOCATABLE   :: FPCDataHDF5(:)
INTEGER(HSIZE_T), POINTER :: HSizeFPC(:)
!===================================================================================================================================
! Only required during restart
IF(.NOT.DoRestart) RETURN

#if USE_LOADBALANCE
! Do not try to read the data from .h5 if load balance is performed without creating a .h5 restart file
IF(PerformLoadBalance.AND..NOT.(UseH5IOLoadBalance)) RETURN
#endif /*USE_LOADBALANCE*/

! 1. The MPI root process reads the info and checks data consistency
! Only root reads the values and distributes them via MPI Broadcast
IF(MPIRoot)THEN
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
  ! Check old parameter name
  ContainerName='FloatingPotentialCharge'
  CALL DatasetExists(File_ID,TRIM(ContainerName),FPCExists)
  ! Check for new parameter name
  IF(FPCExists)THEN
    CALL GetDataSize(File_ID,TRIM(ContainerName),nDimsFPC,HSizeFPC)
    CHECKSAFEINT(HSizeFPC(2),4)
    nFPCHDF5=INT(HSizeFPC(2),4)
    DEALLOCATE(HSizeFPC)
    IF(nFPCHDF5.NE.FPC%nUniqueFPCBounds)THEN
      WRITE(UNIT_StdOut,*) "nFPCHDF5 (restart file) must be equal FPC%nUniqueFPCBounds"
      WRITE(UNIT_StdOut,*) "nFPCHDF5             =", nFPCHDF5
      WRITE(UNIT_StdOut,*) "FPC%nUniqueFPCBounds =", FPC%nUniqueFPCBounds
      CALL abort(__STAMP__,'Restarting with a different number of FPC boundaries, which is not implemented!')
    END IF ! nFPCHDF5.NE.FPC%nUniqueFPCBounds

    ! Allocate the containers
    ALLOCATE(FPCDataHDF5(FPC%nUniqueFPCBounds))
    CALL ReadArray(TRIM(ContainerName) , 2 , (/1_IK , INT(FPC%nUniqueFPCBounds,IK)/) , 0_IK , 1 , RealArray=FPCDataHDF5)
    TmpStr=''
    DO iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
      ! Output in this format: ", 1: 1.312e2 [C]" + ", 2: 3.352e3 [C]" + ...
      WRITE(TmpStr,'(A,I0,A,ES10.3,A)') TRIM(TmpStr)//', ',iUniqueFPCBC,': ',FPCDataHDF5(iUniqueFPCBC),' [C]'
    END DO ! iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
    TmpStr(1:1) = ' ' ! Remove first comma
    TmpStr      = ADJUSTL(TmpStr) ! Remove leading white spaces
    WRITE(UNIT_stdOut,'(A)') " Read floating boundary condition charges from restart file ["//TRIM(RestartFile)//"]: "//TRIM(TmpStr)
    FPC%Charge(:) = FPCDataHDF5
    DEALLOCATE(FPCDataHDF5)
  END IF ! FPCExists
  CALL CloseDataFile()
END IF ! MPIRoot

#if USE_MPI
! 2. The MPI root process distributes the information among the sub-communicator processes for each FPC
CALL SynchronizeChargeOnFPC()
#endif /*USE_MPI*/

END SUBROUTINE ReadFPCDataFromH5


!===================================================================================================================================
!> Read the electric potential (voltage) that was applied on each EPC boundary from .h5 state file.
!> 1. The MPI root process reads the info and checks data consistency
!> 2. The MPI root process distributes the information among the sub-communicator processes for each EPC
!===================================================================================================================================
SUBROUTINE ReadEPCDataFromH5()
! MODULES
USE MOD_io_hdf5
USE MOD_Globals            ,ONLY: UNIT_stdOut,MPIRoot,IERROR,IK,abort
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
#endif /*USE_LOADBALANCE*/
USE MOD_IO_HDF5            ,ONLY: OpenDataFile,CloseDataFile,File_ID
USE MOD_Restart_Vars       ,ONLY: DoRestart,RestartFile
USE MOD_HDF5_Input         ,ONLY: DatasetExists,ReadArray,GetDataSize
USE MOD_HDG_Vars           ,ONLY: EPC
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
! Space-separated list of input and output types. Use: (int|real|logical|...)_(in|out|inout)_dim(n)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iUniqueEPCBC,nEPCHDF5,nDimsEPC
CHARACTER(255)     :: ContainerName
CHARACTER(1000)    :: TmpStr
LOGICAL            :: EPCExists
REAL,ALLOCATABLE   :: EPCDataHDF5(:)
INTEGER(HSIZE_T), POINTER :: HSizeEPC(:)
!===================================================================================================================================
! Only required during restart
IF(.NOT.DoRestart) RETURN

#if USE_LOADBALANCE
! Do not try to read the data from .h5 if load balance is performed without creating a .h5 restart file
IF(PerformLoadBalance.AND..NOT.(UseH5IOLoadBalance)) RETURN
#endif /*USE_LOADBALANCE*/

! 1. The MPI root process reads the info and checks data consistency
! Only root reads the values and distributes them via MPI Broadcast
IF(MPIRoot)THEN
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
  ! Check old parameter name
  ContainerName='ElectricPotentialCondition'
  CALL DatasetExists(File_ID,TRIM(ContainerName),EPCExists)
  ! Check for new parameter name
  IF(EPCExists)THEN
    CALL GetDataSize(File_ID,TRIM(ContainerName),nDimsEPC,HSizeEPC)
    CHECKSAFEINT(HSizeEPC(2),4)
    nEPCHDF5=INT(HSizeEPC(2),4)
    DEALLOCATE(HSizeEPC)
    IF(nEPCHDF5.NE.EPC%nUniqueEPCBounds)THEN
      WRITE(UNIT_StdOut,*) "nEPCHDF5 (restart file) must be equal EPC%nUniqueEPCBounds"
      WRITE(UNIT_StdOut,*) "nEPCHDF5             =", nEPCHDF5
      WRITE(UNIT_StdOut,*) "EPC%nUniqueEPCBounds =", EPC%nUniqueEPCBounds
      CALL abort(__STAMP__,'Restarting with a different number of EPC boundaries, which is not implemented!')
    END IF ! nEPCHDF5.NE.EPC%nUniqueEPCBounds

    ! Allocate the containers
    ALLOCATE(EPCDataHDF5(EPC%nUniqueEPCBounds))
    CALL ReadArray(TRIM(ContainerName) , 2 , (/1_IK , INT(EPC%nUniqueEPCBounds,IK)/) , 0_IK , 1 , RealArray=EPCDataHDF5)
    TmpStr=''
    DO iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
      ! Output in this format: ", 1: 1.312e2 [V]" + ", 2: 3.352e3 [V]" + ...
      WRITE(TmpStr,'(A,I0,A,ES10.3,A)') TRIM(TmpStr)//', ',iUniqueEPCBC,': ',EPCDataHDF5(iUniqueEPCBC),' [V]'
    END DO ! iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
    TmpStr(1:1) = ' ' ! Remove first comma
    TmpStr      = ADJUSTL(TmpStr) ! Remove leading white spaces
    WRITE(UNIT_stdOut,'(A)') " Read electric potential condition from restart file ["//TRIM(RestartFile)//"]: "//TRIM(TmpStr)
    EPC%Voltage(:) = EPCDataHDF5
    DEALLOCATE(EPCDataHDF5)
  END IF ! EPCExists
  CALL CloseDataFile()
END IF ! MPIRoot

#if USE_MPI
! 2. The MPI root process distributes the information among the sub-communicator processes for each EPC
CALL SynchronizeVoltageOnEPC()
#endif /*USE_MPI*/

END SUBROUTINE ReadEPCDataFromH5


#if USE_MPI
!===================================================================================================================================
!> The MPI root process distributes the information among the sub-communicator processes for each FPC
!===================================================================================================================================
SUBROUTINE SynchronizeChargeOnFPC()
! MODULES
USE mpi_f08
USE MOD_HDG_Vars ,ONLY: FPC
USE MOD_Globals  ,ONLY: IERROR,MPI_COMM_NULL,MPI_DOUBLE_PRECISION
! insert modules here
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iUniqueFPCBC
!===================================================================================================================================
DO iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
  IF(FPC%COMM(iUniqueFPCBC)%UNICATOR.NE.MPI_COMM_NULL)THEN
    ! Broadcast from root to other processors on the sub-communicator
    CALL MPI_BCAST(FPC%Charge(iUniqueFPCBC), 1, MPI_DOUBLE_PRECISION, 0, FPC%COMM(iUniqueFPCBC)%UNICATOR, IERROR)
  END IF ! FPC%COMM(iUniqueFPCBC)%UNICATOR.NE.MPI_COMM_NULL
END DO ! iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
END SUBROUTINE SynchronizeChargeOnFPC


!===================================================================================================================================
!> The MPI root process distributes the information among the sub-communicator processes for each EPC
!===================================================================================================================================
SUBROUTINE SynchronizeVoltageOnEPC()
! MODULES
USE mpi_f08
USE MOD_HDG_Vars ,ONLY: EPC
USE MOD_Globals  ,ONLY: IERROR,MPI_COMM_NULL,MPI_DOUBLE_PRECISION
! insert modules here
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iUniqueEPCBC
!===================================================================================================================================
DO iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
  IF(EPC%COMM(iUniqueEPCBC)%UNICATOR.NE.MPI_COMM_NULL)THEN
    ! Broadcast from root to other processors on the sub-communicator
    CALL MPI_BCAST(EPC%Voltage(iUniqueEPCBC), 1, MPI_DOUBLE_PRECISION, 0, EPC%COMM(iUniqueEPCBC)%UNICATOR, IERROR)
  END IF ! EPC%COMM(iUniqueEPCBC)%UNICATOR.NE.MPI_COMM_NULL
END DO ! iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
END SUBROUTINE SynchronizeVoltageOnEPC
#endif /*USE_MPI*/
#endif /*USE_HDG*/

END MODULE MOD_HDG_Readin