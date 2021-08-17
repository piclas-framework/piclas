!==================================================================================================================================
! Copyright (c) 2018 - 2019 Marcel Pfeiffer
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

MODULE MOD_RadTrans_Output
!===================================================================================================================================
! Module for DSMC Sampling and Output
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE WriteRadiationToHDF5
  MODULE PROCEDURE WriteRadiationToHDF5
END INTERFACE


!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: WriteRadiationToHDF5 !, WriteSurfSampleToHDF5
!===================================================================================================================================

CONTAINS

SUBROUTINE WriteRadiationToHDF5()
!===================================================================================================================================
! Writes Radiation values to HDF5
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_PreProc
  USE MOD_io_HDF5
  USE MOD_HDF5_output         ,ONLY: WriteArrayToHDF5,WriteAttributeToHDF5,WriteHDF5Header
  USE MOD_Mesh_Vars           ,ONLY: offsetElem,nGlobalElems, MeshFile
  USE MOD_RadiationTrans_Vars ,ONLY: RadiationElemAbsEnergy_Shared
  USE MOD_RadiationTrans_Vars ,ONLY: Radiation_Emission_Spec_Total, RadTransPhotPerCell
  USE MOD_Globals_Vars        ,ONLY: ProjectName
  USE MOD_Particle_Mesh_Vars  ,ONLY: ElemVolume_Shared
  USE MOD_Radiation_Vars      ,ONLY: RadiationSwitches, Radiation_ElemEnergy_Species, RadiationParameter, Radiation_Absorption_Spec
  USE MOD_Particle_Vars       ,ONLY: nSpecies
  USE MOD_Mesh_Tools         ,ONLY: GetCNElemID
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
  CHARACTER(LEN=255)                  :: FileString,Statedummy
  CHARACTER(LEN=255)                  :: SpecID
  INTEGER                             :: nVal, iElem, nVar, iSpec, nVarCount, nVarSpec, CNElemID, iWave
  REAL, ALLOCATABLE                   :: TempOutput(:,:)
  CHARACTER(LEN=255), ALLOCATABLE     :: StrVarNames(:)
  REAL                                :: AbsTotal
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO') ' WRITE Radiation TO HDF5 FILE...'
  FileString=TRIM(ProjectName)//'_RadiationState.h5'
  Statedummy = 'RadiationState'
  IF (RadiationSwitches%RadType.EQ.1) THEN
    nVarSpec=2               ! _Emission, _Absorption
    nVar=nVarSpec*nSpecies+4 ! nVarSpec + Total_Emission, Total_Absorption, Total_Heatflux, and Total_PhotonNum
  ELSE
    nVar=4
  END IF

  ALLOCATE(StrVarNames(nVar))
  ALLOCATE(TempOutput(nVar, PP_nElems))

  IF (RadiationSwitches%RadType.EQ.1) THEN
    nVarCount=0
    DO iSpec=1, nSpecies
      WRITE(SpecID,'(I3.3)') iSpec
      StrVarNames(nVarCount+1)='Spec'//TRIM(SpecID)//'_Emission'
      StrVarNames(nVarCount+2)='Spec'//TRIM(SpecID)//'_Absorption'
      nVarCount=nVarCount+nVarSpec

    END DO
    StrVarNames(nVarCount+1)='Total_Emission'
    StrVarNames(nVarCount+2)='Total_Absorption'
    StrVarNames(nVarCount+3)='Total_Heatflux'
    StrVarNames(nVarCount+4)='Total_PhotonNum'
  ELSE
    StrVarNames(1)='Total_Emission'
    StrVarNames(2)='Total_Absorption'
    StrVarNames(3)='Total_Heatflux'
    StrVarNames(4)='Total_PhotonNum'
  END IF

  IF(MPIRoot) THEN
    CALL OpenDataFile(FileString,create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteHDF5Header(Statedummy,File_ID)
    CALL WriteAttributeToHDF5(File_ID,'VarNamesAdd',nVar,StrArray=StrVarNames)
    CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=(/TRIM(MeshFile)/))
    CALL CloseDataFile()
  END IF
#if USE_MPI
  CALL MPI_ExchangeRadiationInfo()
#endif /*USE_MPI*/

  CALL OpenDataFile(FileString,create=.false.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_WORLD)
 
  IF (RadiationSwitches%RadType.EQ.1) THEN
    DO iElem=1,PP_nElems
      CNElemID = GetCNElemID(iElem+offSetElem)
      nVarCount=0
      AbsTotal=0.
      DO iSpec=1, nSpecies ! Sum over absorbtion coefficient to determine absorbed energy portion per species
        AbsTotal = AbsTotal + Radiation_ElemEnergy_Species(iSpec,CNElemID,2)
      END DO
      DO iSpec=1, nSpecies
        TempOutput(nVarCount+1, iElem) = Radiation_ElemEnergy_Species(iSpec,CNElemID,1)
!        TempOutput(nVarCount+2, iElem) = Radiation_ElemEnergy_Species(iSpec,iElem,2) !abs coefficient
        TempOutput(nVarCount+2, iElem) = &
             MAX(Radiation_ElemEnergy_Species(iSpec,CNElemID,2)/AbsTotal * RadiationElemAbsEnergy_Shared(iElem+offSetElem)/ ElemVolume_Shared(CNElemID),0.) !lost energy
        nVarCount=nVarCount+nVarSpec
      END DO
      TempOutput((nVarSpec*nSpecies+1), iElem)  = SUM(Radiation_ElemEnergy_Species(:,CNElemID,1))
      TempOutput((nVarSpec*nSpecies+2), iElem)  = RadiationElemAbsEnergy_Shared(iElem+offSetElem)/ ElemVolume_Shared(CNElemID)
      TempOutput(nVarSpec*nSpecies+3, iElem) = SUM(Radiation_ElemEnergy_Species(:,CNElemID,1))- RadiationElemAbsEnergy_Shared(iElem+offSetElem)/ ElemVolume_Shared(CNElemID)
      TempOutput(nVarSpec*nSpecies+4, iElem) = RadTransPhotPerCell(CNElemID)
    END DO
  ELSE IF (RadiationSwitches%RadType.EQ.2) THEN
    DO iElem=1, PP_nElems
      CNElemID = GetCNElemID(iElem+offSetElem)
      TempOutput(1, iElem) = Radiation_Emission_Spec_Total(CNElemID)
      TempOutput(2, iElem) = RadiationElemAbsEnergy_Shared(iElem+offSetElem)/ElemVolume_Shared(CNElemID)
      TempOutput(3, iElem) = Radiation_Emission_Spec_Total(CNElemID)- RadiationElemAbsEnergy_Shared(iElem+offSetElem)/ElemVolume_Shared(CNElemID)
      TempOutput(4, iElem)  = RadTransPhotPerCell(CNElemID)
    END DO
  ELSE
    DO iElem=1, PP_nElems
      CNElemID = GetCNElemID(iElem+offSetElem)
      TempOutput(1, iElem) = Radiation_Emission_Spec_Total(CNElemID)
      TempOutput(2, iElem) = 0.0
      DO iWave = 1, RadiationParameter%WaveLenDiscr
        TempOutput(2, iElem) = TempOutput(2, iElem) + Radiation_Absorption_Spec(iWave, iElem+offSetElem) * RadiationParameter%WaveLenIncr
      END DO
      TempOutput(3, iElem) = Radiation_Emission_Spec_Total(CNElemID) - TempOutput(2, iElem)
      TempOutput(4, iElem) = RadTransPhotPerCell(CNElemID)
    END DO
  END IF

  nVal=nGlobalElems  ! For the MPI case this must be replaced by the global number of elements (sum over all procs)
  CALL WriteArrayToHDF5(DataSetName='ElemData', rank=2,&
                        nValGlobal=(/nVar, nGlobalElems/),&
                        nVal=      (/nVar,   PP_nElems/),&
                        offset=    (/0, offsetElem /),&
                        collective=.TRUE., RealArray=TempOutput(:,:))
  CALL CloseDataFile()
  SWRITE(*,*) 'DONE'

END SUBROUTINE WriteRadiationToHDF5


#if USE_MPI

SUBROUTINE MPI_ExchangeRadiationInfo()
!===================================================================================================================================
! Writes DSMC state values to HDF5
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_PreProc
  USE MOD_RadiationTrans_Vars,   ONLY : RadiationElemAbsEnergy, RadiationElemAbsEnergy_Shared, RadiationElemAbsEnergy_Shared_Win
  USE MOD_Mesh_Vars,              ONLY : nGlobalElems
  USE MOD_MPI_Shared_Vars
  USE MOD_MPI_Shared
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER       :: MessageSize, iELem
!===================================================================================================================================
! collect the information from the proc-local shadow arrays in the compute-node shared array
MessageSize = nGlobalElems

IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_REDUCE(RadiationElemAbsEnergy,RadiationElemAbsEnergy_Shared,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
ELSE
  CALL MPI_REDUCE(RadiationElemAbsEnergy,0                   ,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
ENDIF
CALL BARRIER_AND_SYNC(RadiationElemAbsEnergy_Shared_Win    ,MPI_COMM_SHARED)

IF(nLeaderGroupProcs.GT.1)THEN
  IF(myComputeNodeRank.EQ.0)THEN
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,RadiationElemAbsEnergy_Shared,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEADERS_SHARED,iError)
  END IF
  
  CALL BARRIER_AND_SYNC(RadiationElemAbsEnergy_Shared_Win    ,MPI_COMM_SHARED)
END IF

END SUBROUTINE MPI_ExchangeRadiationInfo


!SUBROUTINE MPI_ExchangeRadiationSurfData() 
!!===================================================================================================================================
!! exchange the surface data
!! only processes with samling sides in their halo region and the original process participate on the communication
!! structure is similar to particle communication
!! each process sends his halo-information directly to the origin process by use of a list, containing the surfsideids for sending
!! the receiving process adds the new data to his own sides
!!===================================================================================================================================
!! MODULES                                                                                                                          !
!!----------------------------------------------------------------------------------------------------------------------------------!
!USE MOD_Globals
!USE MOD_Particle_Boundary_Vars      ,ONLY:SurfComm
!USE MOD_Particle_MPI_Vars           ,ONLY:SurfSendBuf,SurfRecvBuf,SurfExchange
!USE MOD_RadiationTrans_Vars,        ONLY:PhotonSampWall
!!----------------------------------------------------------------------------------------------------------------------------------!
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!! INPUT VARIABLES 
!!----------------------------------------------------------------------------------------------------------------------------------!
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                         :: MessageSize,nValues,iSurfSide,SurfSideID
!INTEGER                         :: iPos,iProc
!INTEGER                         :: recv_status_list(1:MPI_STATUS_SIZE,1:SurfCOMM%nMPINeighbors)
!!===================================================================================================================================

!nValues = 2
!!
!! open receive buffer
!DO iProc=1,SurfCOMM%nMPINeighbors
!  IF(SurfExchange%nSidesRecv(iProc).EQ.0) CYCLE
!  MessageSize=SurfExchange%nSidesRecv(iProc)*nValues
!  CALL MPI_IRECV( SurfRecvBuf(iProc)%content                   &
!                , MessageSize                                  &
!                , MPI_DOUBLE_PRECISION                         &
!                , SurfCOMM%MPINeighbor(iProc)%NativeProcID     &
!                , 1009                                         &
!                , SurfCOMM%COMM                                &
!                , SurfExchange%RecvRequest(iProc)              &
!                , IERROR )
!END DO ! iProc

!! build message
!DO iProc=1,SurfCOMM%nMPINeighbors
!  IF(SurfExchange%nSidesSend(iProc).EQ.0) CYCLE
!  iPos=0
!  SurfSendBuf(iProc)%content = 0.
!  DO iSurfSide=1,SurfExchange%nSidesSend(iProc)
!    SurfSideID=SurfCOMM%MPINeighbor(iProc)%SendList(iSurfSide)
!    SurfSendBuf(iProc)%content(iPos+1:iPos+nValues)= PhotonSampWall(1:2,SurfSideID)
!    iPos=iPos+nValues
!  END DO ! iSurfSide=1,nSurfExchange%nSidesSend(iProc)
!END DO

!! send message
!DO iProc=1,SurfCOMM%nMPINeighbors
!  IF(SurfExchange%nSidesSend(iProc).EQ.0) CYCLE
!  MessageSize=SurfExchange%nSidesSend(iProc)*nValues
!  CALL MPI_ISEND( SurfSendBuf(iProc)%content               &
!                , MessageSize                              &
!                , MPI_DOUBLE_PRECISION                     &
!                , SurfCOMM%MPINeighbor(iProc)%NativeProcID &
!                , 1009                                     &
!                , SurfCOMM%COMM                            &
!                , SurfExchange%SendRequest(iProc)          &
!                , IERROR )
!END DO ! iProc                                                

!! 4) Finish Received number of particles
!DO iProc=1,SurfCOMM%nMPINeighbors
!  IF(SurfExchange%nSidesSend(iProc).NE.0) THEN
!    CALL MPI_WAIT(SurfExchange%SendRequest(iProc),MPIStatus,IERROR)
!    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
!__STAMP__&
!          ,' MPI Communication error', IERROR)
!  END IF
!  IF(SurfExchange%nSidesRecv(iProc).NE.0) THEN
!    CALL MPI_WAIT(SurfExchange%RecvRequest(iProc),recv_status_list(:,iProc),IERROR)
!    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
!__STAMP__&
!          ,' MPI Communication error', IERROR)
!  END IF
!END DO ! iProc

!! add data do my list
!DO iProc=1,SurfCOMM%nMPINeighbors
!  IF(SurfExchange%nSidesRecv(iProc).EQ.0) CYCLE
!  iPos=0
!  DO iSurfSide=1,SurfExchange%nSidesRecv(iProc)
!    SurfSideID=SurfCOMM%MPINeighbor(iProc)%RecvList(iSurfSide)   
!    PhotonSampWall(1:2,SurfSideID)=PhotonSampWall(1:2,SurfSideID) &
!                                     +SurfRecvBuf(iProc)%content(iPos+1:iPos+nValues)
!    iPos=iPos+nValues
!  END DO ! iSurfSide=1,nSurfExchange%nSidesSend(iProc)
!END DO ! iProc

!END SUBROUTINE MPI_ExchangeRadiationSurfData
#endif /*USE_MPI*/


!SUBROUTINE WriteSurfSampleToHDF5()
!!===================================================================================================================================
!!> write the final values of the surface sampling to a HDF5 state file
!!> additional performs all the final required computations
!!===================================================================================================================================
!! MODULES                                                                                                                          !
!!----------------------------------------------------------------------------------------------------------------------------------!
!USE MOD_Globals
!USE MOD_IO_HDF5
!USE MOD_Globals_Vars,               ONLY:ProjectName
!USE MOD_Particle_Boundary_Vars,     ONLY:SurfMesh,offSetSurfSide!,nSurfSample
!USE MOD_HDF5_Output,                ONLY:WriteAttributeToHDF5,WriteArrayToHDF5,WriteHDF5Header
!USE MOD_Particle_Boundary_Vars,     ONLY:SurfCOMM,nSurfBC,SurfBCName
!USE MOD_Mesh_Vars,                  ONLY:MeshFile
!USE MOD_RadiationTrans_Vars,        ONLY:PhotonSampWall
!!----------------------------------------------------------------------------------------------------------------------------------!
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!! INPUT VARIABLES 
!!----------------------------------------------------------------------------------------------------------------------------------!
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!CHARACTER(LEN=255)                  :: FileString,Statedummy
!CHARACTER(LEN=255)                  :: H5_Name
!CHARACTER(LEN=255)                  :: NodeTypeTemp
!CHARACTER(LEN=255),ALLOCATABLE      :: Str2DVarNames(:)
!INTEGER                             :: nVar2D
!REAL                                :: tstart,tend
!REAL, ALLOCATABLE                   :: helpArray(:,:)
!!===================================================================================================================================
!#if USE_MPI
!! Return if not a sampling leader
!IF (MPI_COMM_LEADERS_SURF.EQ.MPI_COMM_NULL) RETURN
!CALL MPI_BARRIER(MPI_COMM_LEADERS_SURF,iERROR)

!! Return if no sampling sides
!IF (nSurfTotalSides      .EQ.0) RETURN
!#endif /*USE_MPI*/
!IF (mySurfRank.EQ.0) THEN
!  WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE Radiation SurfSTATE TO HDF5 FILE...'
!  tstart=LOCALTIME()
!END IF

!FileString=TRIM(ProjectName)//'_RadiationSurfState.h5'
!nVar2D = 2

!! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
!#if USE_MPI
!IF (mySurfRank.EQ.0) THEN
!#endif /*USE_MPI*/
!  CALL OpenDataFile(FileString,create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
!  Statedummy = 'RadiationSurfState'
!  ! Write file header
!  CALL WriteHDF5Header(Statedummy,File_ID)
!  CALL WriteAttributeToHDF5(File_ID,'RadiationnSurfSample',1,IntegerScalar=1)
!  CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=(/TRIM(MeshFile)/))
!  CALL WriteAttributeToHDF5(File_ID,'BC_Surf',nSurfBC,StrArray=SurfBCName)
!  CALL WriteAttributeToHDF5(File_ID,'N',1,IntegerScalar=1)
!  NodeTypeTemp='VISU'
!  CALL WriteAttributeToHDF5(File_ID,'NodeType',1,StrScalar=(/NodeTypeTemp/))

!  ALLOCATE(Str2DVarNames(1:nVar2D))
!  ! fill varnames for total values
!  Str2DVarNames(1) ='PhotonCount'
!  Str2DVarNames(2) ='HeatFlux'

!  CALL WriteAttributeToHDF5(File_ID,'VarNamesSurface',nVar2D,StrArray=Str2DVarNames)

!   CALL CloseDataFile()
!  DEALLOCATE(Str2DVarNames)
!#if USE_MPI
!END IF
!CALL MPI_BARRIER(MPI_COMM_LEADERS_SURF,iERROR)
!CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_LEADERS_SURF)
!#else
!CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
!#endif /*USE_MPI*/


!WRITE(H5_Name,'(A)') 'SurfaceData'
!ASSOCIATE (&
!      nGlobalSides         => INT(nOutputSides,IK) ,&
!      LocalnBCSides        => INT(nComputeNodeSurfOutputSides,IK)     ,&
!      offsetSurfSide       => INT(offsetComputeNodeSurfOutputSide,IK)        ,&
!      nVar2D               => INT(nVar2D,IK))

!  ALLOCATE(helpArray(nVar2D,LocalnBCSides))

!  helpArray(1,1:LocalnBCSides)= PhotonSampWall(1,1:LocalnBCSides)
!  !  SurfaceArea should be changed to 1:SurfMesh%nSides if inner sampling sides exist...
!  helpArray(2,1:LocalnBCSides)= PhotonSampWall(2,1:LocalnBCSides)/SurfSideArea(1,1,1:LocalnBCSides)
!  CALL WriteArrayToHDF5(DataSetName=H5_Name            , rank=4                                     , &
!                        nValGlobal =(/nVar2D     , 1, 1 , nGlobalSides/)  , &
!                        nVal       =(/nVar2D           , 1, 1 , LocalnBCSides/)        , &
!                        offset     =(/0_IK, 0_IK       , 0_IK        , offsetSurfSide/), &
!                        collective =.TRUE.         ,&
!                        RealArray=helpArray(1:nVar2D,1:LocalnBCSides))
!  DEALLOCATE(helpArray)
!END ASSOCIATE

!CALL CloseDataFile()

!IF (mySurfRank.EQ.0) THEN
!  tend=LOCALTIME()
!  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',tend-tstart,'s]'
!END IF

!END SUBROUTINE WriteSurfSampleToHDF5

END MODULE MOD_RadTrans_Output
