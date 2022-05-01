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
PUBLIC :: WriteRadiationToHDF5 , WriteSurfSampleToHDF5
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
  USE MOD_RadiationTrans_Vars ,ONLY: RadiationElemAbsEnergy_Shared, RadObservationPointMethod, RadObservation_Emission, RadObservationPoint
  USE MOD_RadiationTrans_Vars ,ONLY: Radiation_Emission_Spec_Total, RadTransPhotPerCell, RadObservation_EmissionPart
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
  REAL                                :: AbsTotal,tempSpecAbs
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
        IF (AbsTotal.GT.0) THEN
          tempSpecAbs = Radiation_ElemEnergy_Species(iSpec,CNElemID,2)/AbsTotal * RadiationElemAbsEnergy_Shared(iElem+offSetElem)/ ElemVolume_Shared(CNElemID)
        ELSE
          tempSpecAbs = 0.0
        END IF
        TempOutput(nVarCount+2, iElem) = MAX(tempSpecAbs,0.) !lost energy
        nVarCount=nVarCount+nVarSpec
      END DO
      TempOutput((nVarSpec*nSpecies+1), iElem)  = Radiation_Emission_Spec_Total(CNElemID) ! SUM(Radiation_ElemEnergy_Species(:,CNElemID,1))
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
  ASSOCIATE (&
      nVar         => INT(nVar,IK) ,&
      nGlobalElems => INT(nGlobalElems,IK)     ,&
      offsetElem   => INT(offsetElem,IK)        ,&
      PP_nElems    => INT(PP_nElems,IK))
    CALL WriteArrayToHDF5(DataSetName='ElemData', rank=2,&
                          nValGlobal=(/nVar, nGlobalElems/),&
                          nVal=      (/nVar,   PP_nElems/),&
                          offset=    (/0_IK, offsetElem /),&
                          collective=.TRUE., RealArray=TempOutput(:,:))
  END ASSOCIATE
  CALL CloseDataFile()
  SWRITE(*,*) 'DONE'
  
  CALL WriteSurfSampleToHDF5()
  
  IF (RadObservationPointMethod.GT.0) THEN    
    IF (myRank.EQ.0) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,RadObservation_Emission,RadiationParameter%WaveLenDiscrCoarse,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
    ELSE
      CALL MPI_REDUCE(RadObservation_Emission,0                   ,RadiationParameter%WaveLenDiscrCoarse,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
    ENDIF
    IF (myRank.EQ.0) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,RadObservation_EmissionPart,RadiationParameter%WaveLenDiscrCoarse,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
    ELSE
      CALL MPI_REDUCE(RadObservation_EmissionPart,0                   ,RadiationParameter%WaveLenDiscrCoarse,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
    ENDIF
    IF (myRank.EQ.0) THEN
      OPEN(unit=20,file='Radiation_ObservationPoint.csv', status='replace',action='write')
      WRITE(20,*) 'x,y1,y2'
      IF (RadObservationPointMethod.EQ.1) THEN    
        IF (RadiationParameter%WaveLenReductionFactor.NE.1) THEN
          DO iWave=1, RadiationParameter%WaveLenDiscrCoarse
            WRITE(20,*) RadiationParameter%WaveLenCoarse(iWave)*1.E10,',',RadObservation_Emission(iWave)/RadObservationPoint%Area,',',RadObservation_EmissionPart(iWave)
          END DO
        ELSE
          DO iWave=1, RadiationParameter%WaveLenDiscrCoarse
            WRITE(20,*) RadiationParameter%WaveLen(iWave)*1.E10,',',RadObservation_Emission(iWave)/RadObservationPoint%Area,',',RadObservation_EmissionPart(iWave)
          END DO
        END IF
      ELSEIF (RadObservationPointMethod.EQ.2) THEN
        IF (RadiationParameter%WaveLenReductionFactor.NE.1) THEN
          DO iWave=1, RadiationParameter%WaveLenDiscrCoarse
            WRITE(20,*) RadiationParameter%WaveLenCoarse(iWave)*1.E10,',',RadObservation_Emission(iWave),',',RadObservation_EmissionPart(iWave)
          END DO
        ELSE
          DO iWave=1, RadiationParameter%WaveLenDiscrCoarse
            WRITE(20,*) RadiationParameter%WaveLen(iWave)*1.E10,',',RadObservation_Emission(iWave),',',RadObservation_EmissionPart(iWave)
          END DO
        END IF
      END IF
      CLOSE(unit=20)
    END IF
  END IF

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


SUBROUTINE MPI_ExchangeRadiationSurfData() 
!===================================================================================================================================
! exchange the surface data
! only processes with samling sides in their halo region and the original process participate on the communication
! structure is similar to particle communication
! each process sends his halo-information directly to the origin process by use of a list, containing the surfsideids for sending
! the receiving process adds the new data to his own sides
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Boundary_Vars      ,ONLY:SurfOnNode, SurfMapping, nComputeNodeSurfTotalSides, GlobalSide2SurfSide
USE MOD_Particle_MPI_Vars           ,ONLY:SurfSendBuf,SurfRecvBuf
USE MOD_RadiationTrans_Vars         ,ONLY:PhotonSampWall, PhotonSampWall_Shared, PhotonSampWall_Shared_Win
USE MOD_MPI_Shared_Vars             ,ONLY:MPI_COMM_LEADERS_SURF, MPI_COMM_SHARED, nSurfLeaders,myComputeNodeRank,mySurfRank
USE MOD_MPI_Shared              
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: MessageSize,nValues,iSurfSide,SurfSideID, SideID
INTEGER                         :: iPos,iProc
INTEGER                         :: RecvRequest(0:nSurfLeaders-1),SendRequest(0:nSurfLeaders-1)
!===================================================================================================================================
! nodes without sampling surfaces do not take part in this routine
IF (.NOT.SurfOnNode) RETURN

MessageSize = 2*nComputeNodeSurfTotalSides
IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_REDUCE(PhotonSampWall,PhotonSampWall_Shared,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
ELSE
  CALL MPI_REDUCE(PhotonSampWall,0                   ,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
ENDIF

CALL BARRIER_AND_SYNC(PhotonSampWall_Shared_Win         ,MPI_COMM_SHARED)

! prepare buffers for surf leader communication
IF (myComputeNodeRank.EQ.0) THEN
  nValues = 2
  
  ! open receive buffer
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    ! Only open recv buffer if we are expecting sides from this leader node
    IF (SurfMapping(iProc)%nRecvSurfSides.EQ.0) CYCLE

    ! Message is sent on MPI_COMM_LEADERS_SURF, so rank is indeed iProc
    MessageSize = SurfMapping(iProc)%nRecvSurfSides * nValues
    CALL MPI_IRECV( SurfRecvBuf(iProc)%content                   &
                  , MessageSize                                  &
                  , MPI_DOUBLE_PRECISION                         &
                  , iProc                                        &
                  , 1209                                         &
                  , MPI_COMM_LEADERS_SURF                        &
                  , RecvRequest(iProc)                           &
                  , IERROR)
  END DO ! iProc

  ! build message
  DO iProc = 0,nSurfLeaders-1
    ! Ignore myself
    IF (iProc .EQ. mySurfRank) CYCLE
    ! Only assemble message if we are expecting sides to send to this leader node
    IF (SurfMapping(iProc)%nSendSurfSides.EQ.0) CYCLE

    ! Nullify everything
    iPos = 0
    SurfSendBuf(iProc)%content = 0.
    DO iSurfSide = 1,SurfMapping(iProc)%nSendSurfSides
      SideID     = SurfMapping(iProc)%SendSurfGlobalID(iSurfSide)
      SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)
      ! Assemble message
      SurfSendBuf(iProc)%content(iPos+1:iPos+2) = PhotonSampWall_Shared(:,SurfSideID)
      iPos = iPos + 2
      PhotonSampWall_Shared(:,SurfSideID)=0.
    END DO ! iSurfSide = 1,SurfMapping(iProc)%nSendSurfSides
  END DO

  ! send message
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE
    ! Only open recv buffer if we are expecting sides from this leader node
    IF (SurfMapping(iProc)%nSendSurfSides.EQ.0) CYCLE

    ! Message is sent on MPI_COMM_LEADERS_SURF, so rank is indeed iProc
    MessageSize = SurfMapping(iProc)%nSendSurfSides * nValues
    CALL MPI_ISEND( SurfSendBuf(iProc)%content                   &
                  , MessageSize                                  &
                  , MPI_DOUBLE_PRECISION                         &
                  , iProc                                        &
                  , 1209                                         &
                  , MPI_COMM_LEADERS_SURF                        &
                  , SendRequest(iProc)                           &
                  , IERROR)
  END DO ! iProc

  ! Finish received number of sampling surfaces
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    IF (SurfMapping(iProc)%nSendSurfSides.NE.0) THEN
      CALL MPI_WAIT(SendRequest(iProc),MPIStatus,IERROR)
      IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error',IERROR)
    END IF

    IF (SurfMapping(iProc)%nRecvSurfSides.NE.0) THEN
      CALL MPI_WAIT(RecvRequest(iProc),MPIStatus,IERROR)
      IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error',IERROR)
    END IF
  END DO ! iProc

  ! add data do my list
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE
    ! Only open recv buffer if we are expecting sides from this leader node
    IF (SurfMapping(iProc)%nRecvSurfSides.EQ.0) CYCLE

    iPos=0
    DO iSurfSide = 1,SurfMapping(iProc)%nRecvSurfSides
      SideID     = SurfMapping(iProc)%RecvSurfGlobalID(iSurfSide)
      SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)
      PhotonSampWall_Shared(:,SurfSideID) = PhotonSampWall_Shared(:,SurfSideID) &
                                             + SurfRecvBuf(iProc)%content(iPos+1:iPos+2)
      iPos = iPos + 2
    END DO ! iSurfSide = 1,SurfMapping(iProc)%nRecvSurfSides
     ! Nullify buffer
    SurfRecvBuf(iProc)%content = 0.
  END DO ! iProc
END IF

CALL BARRIER_AND_SYNC(PhotonSampWall_Shared_Win         ,MPI_COMM_SHARED)

END SUBROUTINE MPI_ExchangeRadiationSurfData
#endif /*USE_MPI*/


SUBROUTINE WriteSurfSampleToHDF5()
!===================================================================================================================================
!> write the final values of the surface sampling to a HDF5 state file
!> additional performs all the final required computations
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_Globals_Vars,               ONLY:ProjectName
USE MOD_Particle_Boundary_Vars,     ONLY:SurfSideArea, nComputeNodeSurfOutputSides,noutputsides, nSurfTotalSides, nSurfBC
USE MOD_Particle_Boundary_Vars,     ONLY:offsetComputeNodeSurfOutputSide, SurfBCName, SurfSideArea_Shared, nComputeNodeSurfSides
USE MOD_Particle_Boundary_Vars,     ONLY:SurfSide2GlobalSide, GlobalSide2SurfSide
USE MOD_HDF5_Output,                ONLY:WriteAttributeToHDF5,WriteArrayToHDF5,WriteHDF5Header 
USE MOD_Mesh_Vars,                  ONLY:MeshFile
USE MOD_Particle_Mesh_Vars,         ONLY:SideInfo_Shared
USE MOD_RadiationTrans_Vars,        ONLY:PhotonSampWall, PhotonSampWall_Shared
USE MOD_MPI_Shared_Vars,            ONLY:MPI_COMM_LEADERS_SURF,mySurfRank
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)                  :: FileString,Statedummy
CHARACTER(LEN=255)                  :: H5_Name
CHARACTER(LEN=255)                  :: NodeTypeTemp
CHARACTER(LEN=255),ALLOCATABLE      :: Str2DVarNames(:)
INTEGER                             :: nVar2D, GlobalSideID, iSurfSide, OutputCounter, SurfSideNb
REAL                                :: tstart,tend
REAL, ALLOCATABLE                   :: helpArray(:,:)
!===================================================================================================================================
#if USE_MPI
CALL MPI_ExchangeRadiationSurfData()
! Return if not a sampling leader
IF (MPI_COMM_LEADERS_SURF.EQ.MPI_COMM_NULL) RETURN
CALL MPI_BARRIER(MPI_COMM_LEADERS_SURF,iERROR)

! Return if no sampling sides
IF (nSurfTotalSides      .EQ.0) RETURN
#endif /*USE_MPI*/
IF (mySurfRank.EQ.0) THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE Radiation SurfSTATE TO HDF5 FILE...'
  tstart=LOCALTIME()
END IF

FileString=TRIM(ProjectName)//'_RadiationSurfState.h5'
nVar2D = 2

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
#if USE_MPI
IF (mySurfRank.EQ.0) THEN
#endif /*USE_MPI*/
  CALL OpenDataFile(FileString,create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
  Statedummy = 'RadiationSurfState'
  ! Write file header
  CALL WriteHDF5Header(Statedummy,File_ID)
  CALL WriteAttributeToHDF5(File_ID,'RadiationnSurfSample',1,IntegerScalar=1)
  CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=(/TRIM(MeshFile)/))
  CALL WriteAttributeToHDF5(File_ID,'BC_Surf',nSurfBC,StrArray=SurfBCName)
  CALL WriteAttributeToHDF5(File_ID,'N',1,IntegerScalar=1)
  NodeTypeTemp='VISU'
  CALL WriteAttributeToHDF5(File_ID,'NodeType',1,StrScalar=(/NodeTypeTemp/))

  ALLOCATE(Str2DVarNames(1:nVar2D))
  ! fill varnames for total values
  Str2DVarNames(1) ='PhotonCount'
  Str2DVarNames(2) ='HeatFlux'

  CALL WriteAttributeToHDF5(File_ID,'VarNamesSurface',nVar2D,StrArray=Str2DVarNames)

   CALL CloseDataFile()
  DEALLOCATE(Str2DVarNames)
#if USE_MPI
END IF
CALL MPI_BARRIER(MPI_COMM_LEADERS_SURF,iERROR)
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_LEADERS_SURF)
#else
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
#endif /*USE_MPI*/


WRITE(H5_Name,'(A)') 'SurfaceData'
#if USE_MPI
ASSOCIATE(PhotonSampWall        => PhotonSampWall_Shared           ,&
          SurfSideArea         => SurfSideArea_Shared)
#endif

ASSOCIATE (&
      nGlobalSides         => INT(nOutputSides,IK) ,&
      LocalnBCSides        => INT(nComputeNodeSurfOutputSides,IK)     ,&
      offsetSurfSide       => INT(offsetComputeNodeSurfOutputSide,IK)        ,&
      nVar2D               => INT(nVar2D,IK))

  ALLOCATE(helpArray(nVar2D,LocalnBCSides))
  OutputCounter = 0
  DO iSurfSide = 1,nComputeNodeSurfSides
    GlobalSideID = SurfSide2GlobalSide(SURF_SIDEID,iSurfSide)
    IF(SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID).GT.0) THEN
      IF(GlobalSideID.LT.SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID)) THEN
        SurfSideNb = GlobalSide2SurfSide(SURF_SIDEID,SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID))
        PhotonSampWall(:,iSurfSide) = PhotonSampWall(:,iSurfSide) + PhotonSampWall(:,SurfSideNb)
      ELSE
        CYCLE
      END IF
    END IF
    OutputCounter = OutputCounter + 1
    helpArray(1,OutputCounter)= PhotonSampWall(1,iSurfSide)
    !  SurfaceArea should be changed to 1:SurfMesh%nSides if inner sampling sides exist...
    helpArray(2,OutputCounter)= PhotonSampWall(2,iSurfSide)/SurfSideArea(1,1,iSurfSide)
  END DO
  CALL WriteArrayToHDF5(DataSetName=H5_Name            , rank=4                                     , &
                        nValGlobal =(/nVar2D     , 1_IK, 1_IK , nGlobalSides/)  , &
                        nVal       =(/nVar2D           , 1_IK, 1_IK , LocalnBCSides/)        , &
                        offset     =(/0_IK, 0_IK       , 0_IK        , offsetSurfSide/), &
                        collective =.TRUE.         ,&
                        RealArray=helpArray(1:nVar2D,1:LocalnBCSides))
  DEALLOCATE(helpArray)
END ASSOCIATE

#if USE_MPI
END ASSOCIATE
#endif /*USE_MPI*/

CALL CloseDataFile()

IF (mySurfRank.EQ.0) THEN
  tend=LOCALTIME()
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',tend-tstart,'s]'
END IF

END SUBROUTINE WriteSurfSampleToHDF5

END MODULE MOD_RadTrans_Output
