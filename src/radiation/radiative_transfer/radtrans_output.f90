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
! Module for output of radiative transfer solver
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
PUBLIC :: WriteRadiationToHDF5
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
USE MOD_RadiationTrans_Vars ,ONLY: RadObservationPointMethod, RadObservation_Emission, RadObservationPoint
USE MOD_RadiationTrans_Vars ,ONLY: Radiation_Emission_Spec_Total, RadTransPhotPerCell, RadObservation_EmissionPart
USE MOD_RadiationTrans_Vars ,ONLY: ObservationDoConvolution, RadObservation_Emission_Conv
USE MOD_Globals_Vars        ,ONLY: ProjectName
USE MOD_Particle_Mesh_Vars  ,ONLY: ElemVolume_Shared
USE MOD_Radiation_Vars      ,ONLY: RadiationSwitches, Radiation_ElemEnergy_Species, RadiationParameter, Radiation_Absorption_Spec
USE MOD_Particle_Vars       ,ONLY: nSpecies
USE MOD_Mesh_Tools          ,ONLY: GetCNElemID
USE MOD_Photon_TrackingOutput,ONLY:WritePhotonSurfSampleToHDF5
#if USE_MPI
USE MOD_RadiationTrans_Vars ,ONLY: RadiationElemAbsEnergySpec_Shared, RadiationElemAbsEnergy_Shared
#else
USE MOD_RadiationTrans_Vars ,ONLY: RadiationElemAbsEnergySpec, RadiationElemAbsEnergy
#endif /*USE_MPI*/
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
REAL                                :: tmpPartNum, tmpEmission(2)
INTEGER                             :: iWavetmp(2)
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO') ' WRITE Radiation TO HDF5 FILE...'
FileString=TRIM(ProjectName)//'_RadiationState.h5'
Statedummy = 'RadiationState'
IF (RadiationSwitches%RadType.EQ.1) THEN
  nVarSpec=2               ! _Emission, _Absorption
  nVar=nVarSpec*nSpecies+5 ! nVarSpec + Total_Emission, Total_Absorption, Total_Heatflux, and Total_PhotonNum
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
  StrVarNames(nVarCount+5)='Mean_OpticalDepth'
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

CALL OpenDataFile(FileString,create=.false.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)

#if USE_MPI
ASSOCIATE( RadiationElemAbsEnergySpec => RadiationElemAbsEnergySpec_Shared,&
           RadiationElemAbsEnergy     => RadiationElemAbsEnergy_Shared    )
#endif /*USE_MPI*/

IF (RadiationSwitches%RadType.EQ.1) THEN
  DO iElem=1,PP_nElems
    CNElemID = GetCNElemID(iElem+offSetElem)
    nVarCount=0
    DO iSpec=1, nSpecies
      TempOutput(nVarCount+1, iElem) = Radiation_ElemEnergy_Species(iSpec,CNElemID,1)
      TempOutput(nVarCount+2, iElem) = RadiationElemAbsEnergySpec(iSpec, iElem+offSetElem)/ ElemVolume_Shared(CNElemID)
      nVarCount=nVarCount+nVarSpec
    END DO
    TempOutput((nVarSpec*nSpecies+1), iElem)  = Radiation_Emission_Spec_Total(CNElemID)
    TempOutput((nVarSpec*nSpecies+2), iElem)  = SUM(RadiationElemAbsEnergySpec(:, iElem+offSetElem))/ ElemVolume_Shared(CNElemID)
    TempOutput(nVarSpec*nSpecies+3, iElem) = SUM(Radiation_ElemEnergy_Species(:,CNElemID,1))- SUM(RadiationElemAbsEnergySpec(:, iElem+offSetElem))/ ElemVolume_Shared(CNElemID)
    TempOutput(nVarSpec*nSpecies+4, iElem) = RadTransPhotPerCell(CNElemID)
    IF (RadiationElemAbsEnergy(2,iElem+offSetElem).GT.0) THEN
      TempOutput(nVarSpec*nSpecies+5, iElem) = RadiationElemAbsEnergy(1,iElem+offSetElem)/RadiationElemAbsEnergy(2,iElem+offSetElem)
    ELSE
      TempOutput(nVarSpec*nSpecies+5, iElem) = 0.0
    END IF
  END DO
ELSE IF (RadiationSwitches%RadType.EQ.2) THEN
  DO iElem=1, PP_nElems
    CNElemID = GetCNElemID(iElem+offSetElem)
    TempOutput(1, iElem) = Radiation_Emission_Spec_Total(CNElemID)
    TempOutput(2, iElem) = RadiationElemAbsEnergySpec(1,iElem+offSetElem)/ElemVolume_Shared(CNElemID)
    TempOutput(3, iElem) = Radiation_Emission_Spec_Total(CNElemID)- RadiationElemAbsEnergySpec(1,iElem+offSetElem)/ElemVolume_Shared(CNElemID)
    TempOutput(4, iElem)  = RadTransPhotPerCell(CNElemID)
  END DO
ELSE IF (RadiationSwitches%RadType.EQ.3) THEN
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
ELSE IF (RadiationSwitches%RadType.EQ.4) THEN
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
ELSE
  CALL abort(__STAMP__,' ERROR: Radiation type is not implemented! (unknown case)')
END IF

#if USE_MPI
END ASSOCIATE
#endif /*USE_MPI*/

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

CALL WritePhotonSurfSampleToHDF5()

IF (RadObservationPointMethod.GT.0) THEN
#if USE_MPI
  IF (myRank.EQ.0) THEN
    CALL MPI_REDUCE(MPI_IN_PLACE,RadObservation_Emission,RadiationParameter%WaveLenDiscrCoarse,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
  ELSE
    CALL MPI_REDUCE(RadObservation_Emission,0                   ,RadiationParameter%WaveLenDiscrCoarse,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
  ENDIF
  IF (myRank.EQ.0) THEN
    CALL MPI_REDUCE(MPI_IN_PLACE,RadObservation_EmissionPart,RadiationParameter%WaveLenDiscrCoarse,MPI_INTEGER,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
  ELSE
    CALL MPI_REDUCE(RadObservation_EmissionPart,0                   ,RadiationParameter%WaveLenDiscrCoarse,MPI_INTEGER,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
  ENDIF
#endif /*USE_MPI*/
  IF (myRank.EQ.0) THEN
    IF(ObservationDoConvolution) THEN

      CALL SpectralConvolution(RadObservation_Emission,RadObservation_Emission_Conv)
      OPEN(unit=20,file='Radiation_ObservationPoint.csv', status='replace',action='write')
      WRITE(20,*) 'x,y1,y2,y3'
      IF (RadObservationPointMethod.EQ.1) THEN
        IF (RadiationParameter%WaveLenReductionFactor.NE.1) THEN
          DO iWave=1, RadiationParameter%WaveLenDiscrCoarse
            WRITE(20,*) RadiationParameter%WaveLenCoarse(iWave)*1.E10,',',RadObservation_Emission(iWave)/RadObservationPoint%Area,',',RadObservation_EmissionPart(iWave),',',RadObservation_Emission_Conv(iWave)/RadObservationPoint%Area
          END DO
        ELSE
          IF (RadiationParameter%WaveLenReductionFactorOutput.GT.1) THEN
            tmpPartNum=0.; tmpEmission=0.; iWavetmp(1)=0; iWavetmp(2)=1
            DO iWave=1, RadiationParameter%WaveLenDiscrCoarse
              IF (MOD(iWave,RadiationParameter%WaveLenReductionFactorOutput).EQ.0) THEN
                iWavetmp(1) = iWavetmp(1) + 1
                tmpPartNum = tmpPartNum + RadObservation_EmissionPart(iWave)
                tmpEmission(1) = tmpEmission(1) +  RadObservation_Emission(iWave)
                tmpEmission(2) = tmpEmission(2) +  RadObservation_Emission_Conv(iWave)
                WRITE(20,*) (RadiationParameter%WaveLen(iWavetmp(2))+RadiationParameter%WaveLen(iWavetmp(2)+iWavetmp(1)-1))/2.*1.E10,',',tmpEmission(1)/RadObservationPoint%Area,',',tmpPartNum,',',tmpEmission(2)/RadObservationPoint%Area
                tmpPartNum = 0.; tmpEmission= 0.
                iWavetmp(1)=0; iWavetmp(2)=iWave
              ELSE
                iWavetmp(1) = iWavetmp(1) + 1
                tmpPartNum = tmpPartNum + RadObservation_EmissionPart(iWave)
                tmpEmission(1) = tmpEmission(1) +  RadObservation_Emission(iWave)
                tmpEmission(2) = tmpEmission(2) +  RadObservation_Emission_Conv(iWave)
              END IF
            END DO
          ELSE
            DO iWave=1, RadiationParameter%WaveLenDiscrCoarse
              WRITE(20,*) RadiationParameter%WaveLen(iWave)*1.E10,',',RadObservation_Emission(iWave)/RadObservationPoint%Area,',',RadObservation_EmissionPart(iWave),',',RadObservation_Emission_Conv(iWave)/RadObservationPoint%Area
            END DO
          END IF
        END IF
      ELSEIF (RadObservationPointMethod.EQ.2) THEN
        IF (RadiationParameter%WaveLenReductionFactor.NE.1) THEN
          DO iWave=1, RadiationParameter%WaveLenDiscrCoarse
            WRITE(20,*) RadiationParameter%WaveLenCoarse(iWave)*1.E10,',',RadObservation_Emission(iWave),',',RadObservation_EmissionPart(iWave),',',RadObservation_Emission_Conv(iWave)
          END DO
        ELSE
          IF (RadiationParameter%WaveLenReductionFactorOutput.GT.1) THEN
            tmpPartNum=0.; tmpEmission=0.; iWavetmp(1)=0; iWavetmp(2)=1
            DO iWave=1, RadiationParameter%WaveLenDiscrCoarse
              IF (MOD(iWave,RadiationParameter%WaveLenReductionFactorOutput).EQ.0) THEN
                iWavetmp(1) = iWavetmp(1) + 1
                tmpPartNum = tmpPartNum + RadObservation_EmissionPart(iWave)
                tmpEmission(1) = tmpEmission(1) +  RadObservation_Emission(iWave)
                tmpEmission(2) = tmpEmission(2) +  RadObservation_Emission_Conv(iWave)
                WRITE(20,*) (RadiationParameter%WaveLen(iWavetmp(2))+RadiationParameter%WaveLen(iWavetmp(2)+iWavetmp(1)-1))/2.*1.E10,',',tmpEmission(1),',',tmpPartNum,',',tmpEmission(2)
                tmpPartNum = 0.; tmpEmission= 0.
                iWavetmp(1)=0; iWavetmp(2)=iWave
              ELSE
                iWavetmp(1) = iWavetmp(1) + 1
                tmpPartNum = tmpPartNum + RadObservation_EmissionPart(iWave)
                tmpEmission(1) = tmpEmission(1) +  RadObservation_Emission(iWave)
                tmpEmission(2) = tmpEmission(2) +  RadObservation_Emission_Conv(iWave)
              END IF
            END DO
          ELSE
            DO iWave=1, RadiationParameter%WaveLenDiscrCoarse
              WRITE(20,*) RadiationParameter%WaveLen(iWave)*1.E10,',',RadObservation_Emission(iWave),',',RadObservation_EmissionPart(iWave),',',RadObservation_Emission_Conv(iWave)
            END DO
          END IF
        END IF
      END IF
      CLOSE(unit=20)
    ELSE
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
END IF

END SUBROUTINE WriteRadiationToHDF5


#if USE_MPI
SUBROUTINE MPI_ExchangeRadiationInfo()
!===================================================================================================================================
! MPI routine for output of radiative transfer solver
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_RadiationTrans_Vars ,ONLY: RadiationElemAbsEnergy, RadiationElemAbsEnergy_Shared, RadiationElemAbsEnergy_Shared_Win
USE MOD_RadiationTrans_Vars ,ONLY: RadiationElemAbsEnergySpec, RadiationElemAbsEnergySpec_Shared, RadiationElemAbsEnergySpec_Shared_Win
USE MOD_Mesh_Vars           ,ONLY: nGlobalElems
USE MOD_MPI_Shared_Vars
USE MOD_MPI_Shared
USE MOD_Particle_Vars       ,ONLY: nSpecies
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER       :: MessageSize
!===================================================================================================================================
! collect the information from the proc-local shadow arrays in the compute-node shared array
MessageSize = 2*nGlobalElems

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

MessageSize = nSpecies*nGlobalElems

IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_REDUCE(RadiationElemAbsEnergySpec,RadiationElemAbsEnergySpec_Shared,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
ELSE
  CALL MPI_REDUCE(RadiationElemAbsEnergySpec,0                   ,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
ENDIF
CALL BARRIER_AND_SYNC(RadiationElemAbsEnergySpec_Shared_Win    ,MPI_COMM_SHARED)

IF(nLeaderGroupProcs.GT.1)THEN
  IF(myComputeNodeRank.EQ.0)THEN
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,RadiationElemAbsEnergySpec_Shared,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEADERS_SHARED,iError)
  END IF

  CALL BARRIER_AND_SYNC(RadiationElemAbsEnergySpec_Shared_Win    ,MPI_COMM_SHARED)
END IF

END SUBROUTINE MPI_ExchangeRadiationInfo
#endif /*USE_MPI*/


SUBROUTINE SpectralConvolution(RadObservation_Emission, RadObservation_Emission_Conv)
!===================================================================================================================================
! calculates spectral concolution with slit function/instrumental broadening profile/spectral resolution function
!===================================================================================================================================
! MODULES
! USE MOD_Globals
USE MOD_RadiationTrans_Vars ,ONLY: RadObservationPoint
USE MOD_Radiation_Vars      ,ONLY: RadiationParameter
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)    :: RadObservation_Emission(:)
REAL, INTENT(INOUT) :: RadObservation_Emission_Conv(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: topwidth, basewidth
INTEGER :: iWave_min, iWave, i
REAL    :: topwidth_half, basewidth_half, slope
REAL    :: wavelength_min_base, wavelength_max_base, wavelength_min_top, wavelength_max_top
REAL    :: fractionl, fractionr, delta_base, delta_top
!===================================================================================================================================
topwidth = RadObservationPoint%SlitFunction(1)*1.E-10
basewidth = RadObservationPoint%SlitFunction(2)*1.E-10
iWave_min = 1

basewidth_half = 0.5 * basewidth
topwidth_half  = 0.5 * topwidth
slope          = 1. / (basewidth_half-topwidth_half)
RadObservation_Emission_Conv=0.0
DO iWave=1, RadiationParameter%WaveLenDiscr
  wavelength_min_base = RadiationParameter%WaveLen(iWave) - basewidth_half
  wavelength_max_base = RadiationParameter%WaveLen(iWave) + basewidth_half
  wavelength_min_top  = RadiationParameter%WaveLen(iWave) - topwidth_half
  wavelength_max_top  = RadiationParameter%WaveLen(iWave) + topwidth_half
  ! --- start index determination
  DO WHILE(RadiationParameter%WaveLen(iWave_min+1) .LT. wavelength_min_base)
    iWave_min = iWave_min + 1
  END DO
  ! --- slit function
  DO i = iWave_min, RadiationParameter%WaveLenDiscr-1
    IF(RadiationParameter%WaveLen(i) .LT. wavelength_min_base) THEN
      fractionl = 0.
      IF(RadiationParameter%WaveLen(i+1) .GT. wavelength_min_top) THEN
        STOP 'slit function: step width is too big!'
      END IF
      fractionr  = slope * (RadiationParameter%WaveLen(i+1) - RadiationParameter%WaveLen(iWave) + basewidth_half)
      delta_base = RadiationParameter%WaveLen(i+1) - wavelength_min_base
      delta_top  = 0.
    ELSEIF(RadiationParameter%WaveLen(i+1) .LT. wavelength_min_top) THEN
      fractionl  = slope * (RadiationParameter%WaveLen(i  ) - RadiationParameter%WaveLen(iWave) + basewidth_half)
      fractionr  = slope * (RadiationParameter%WaveLen(i+1) - RadiationParameter%WaveLen(iWave) + basewidth_half)
      delta_base = RadiationParameter%WaveLenIncr
      delta_top  = 0.
    ELSEIF(RadiationParameter%WaveLen(i  ) .LT. wavelength_min_top) THEN
      fractionl  = slope * (RadiationParameter%WaveLen(i  ) - RadiationParameter%WaveLen(iWave) + basewidth_half)
      fractionr  = 1.
      delta_base = wavelength_min_top - RadiationParameter%WaveLen(i)
      delta_top  = RadiationParameter%WaveLen(i+1) - wavelength_min_top
    ELSEIF(RadiationParameter%WaveLen(i+1) .LT. wavelength_max_top) THEN
      fractionl  = 0.
      fractionr  = 0.
      delta_base = 0.
      delta_top  = RadiationParameter%WaveLenIncr
    ELSEIF(RadiationParameter%WaveLen(i  ) .LT. wavelength_max_top) THEN
      fractionl  = 1.
      fractionr  = - slope * (RadiationParameter%WaveLen(i+1) - RadiationParameter%WaveLen(iWave) - basewidth_half)
      delta_base = RadiationParameter%WaveLen(i+1) - wavelength_max_top
      delta_top  = wavelength_max_top - RadiationParameter%WaveLen(i)
    ELSEIF(RadiationParameter%WaveLen(i+1) .LT. wavelength_max_base) THEN
      fractionl  = - slope * (RadiationParameter%WaveLen(i  ) - RadiationParameter%WaveLen(iWave) - basewidth_half)
      fractionr  = - slope * (RadiationParameter%WaveLen(i+1) - RadiationParameter%WaveLen(iWave) - basewidth_half)
      delta_base = RadiationParameter%WaveLenIncr
      delta_top  = 0.
    ELSEIF(RadiationParameter%WaveLen(i) .LT. wavelength_max_base) THEN
      fractionl  = - slope * (RadiationParameter%WaveLen(i  ) - RadiationParameter%WaveLen(iWave) - basewidth_half)
      fractionr  = 0.
      delta_base = wavelength_max_base - RadiationParameter%WaveLen(i)
      delta_top  = 0.
    ELSE
      exit
    END IF

    RadObservation_Emission_Conv(iWave) = RadObservation_Emission_Conv(iWave) &
        + ((fractionl+fractionr)*.5*delta_base+delta_top) &
        * RadObservation_Emission(i+1)*1.E10

  END DO
END DO

END SUBROUTINE SpectralConvolution

END MODULE MOD_RadTrans_Output
