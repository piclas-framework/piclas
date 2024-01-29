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

MODULE MOD_RadiationTrans_Init
!===================================================================================================================================
! Initialization of Radiative Transfer
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE InitRadiationTransport
  MODULE PROCEDURE InitRadiationTransport
END INTERFACE

PUBLIC::InitRadiationTransport, DefineParametersRadiationTrans, HALTON, FinalizeRadiationTransport
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for radiative transfer
!==================================================================================================================================
SUBROUTINE DefineParametersRadiationTrans()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Radiation Transport")

CALL prms%CreateLogicalOption('Radiation-AdaptivePhotonNumEmission', 'HM','.FALSE.')
CALL prms%CreateIntOption('Radiation-RadObservationPointMethod', 'HM','0')
CALL prms%CreateIntOption(    'Radiation-DirectionModel', 'HM','1')
CALL prms%CreateIntOption(    'Radiation-NumPhotonsPerCell', 'HM','1')
CALL prms%CreateIntOption(    'Radiation-AbsorptionModel', 'HM','1')
CALL prms%CreateIntOption(    'Radiation-PhotonPosModel', 'HM','1')
CALL prms%CreateIntOption(    'Radiation-PhotonWaveLengthModel', 'HM','1')
CALL prms%CreateRealArrayOption('Radiation-ObservationMidPoint', 'HM', '0.,0.,0.')
CALL prms%CreateRealArrayOption('Radiation-ObservationSlitFunction', 'Slit function for convolution, trapezoid, 1:topwidth[A] 2:basewidth[A]', '0.,0.')
CALL prms%CreateRealOption('Radiation-ObservationDiameter', 'HM')
CALL prms%CreateRealArrayOption('Radiation-ObservationViewDirection', 'HM', '0.,0.,0.')
CALL prms%CreateRealOption('Radiation-ObservationAngularAperture', 'HM')
CALL prms%CreateLogicalOption('Radiation-ObservationCalcFullSpectra','.FALSE.')
CALL prms%CreateLogicalOption('Radiation-ObservationDoConvolution','Consider instrumental broadening?','.FALSE.')
CALL prms%CreateRealOption('Radiation-ShockTubeDiameter', 'Diameter of shock tube in m', '0.0')
CALL prms%CreateIntOption( 'Radiation-nSurfSample'     , 'Define polynomial degree of radiation BC sampling. Default: nSurfSample (which itself defaults to NGeo)')

END SUBROUTINE DefineParametersRadiationTrans

SUBROUTINE InitRadiationTransport()
!===================================================================================================================================
! Initialization of the radiative transfer solver
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars               ,ONLY: Pi, c
USE MOD_ReadInTools
USE MOD_RadiationTrans_Vars
USE MOD_Mesh_Vars                  ,ONLY: nGlobalElems
USE MOD_Particle_Mesh_Vars         ,ONLY: ElemVolume_Shared,ElemMidPoint_Shared, GEO, nComputeNodeElems
USE MOD_Globals_Vars               ,ONLY: BoltzmannConst, PlanckConst
USE MOD_Particle_Boundary_Sampling ,ONLY: InitParticleBoundarySampling
USE MOD_Particle_Boundary_Vars     ,ONLY: nComputeNodeSurfTotalSides,nSurfSample
USE MOD_Radiation_Vars             ,ONLY: RadiationParameter, Radiation_Emission_Spec, Radiation_Absorption_Spec, RadiationSwitches
USE MOD_Radiation_Vars             ,ONLY: Radiation_Absorption_SpecPercent
USE MOD_RadiationTrans_Vars        ,ONLY: RadObservation_Emission
USE MOD_Radiation                  ,ONLY: radiation_main
USE MOD_DSMC_Vars                  ,ONLY: RadialWeighting
USE MOD_Output                     ,ONLY: PrintStatusLineRadiation
USE MOD_Mesh_Tools                 ,ONLY: GetGlobalElemID
USE MOD_Particle_Vars              ,ONLY: Symmetry, nSpecies
USE MOD_MPI_Shared_Vars
USE MOD_MPI_Shared
USE MOD_Particle_Mesh_Build        ,ONLY: BuildMesh2DInfo
USE MOD_SuperB_Tools               ,ONLY: FindLinIndependentVectors, GramSchmidtAlgo
#if USE_MPI
USE MOD_RadiationTrans_Vars        ,ONLY: RadTransObsVolumeFrac_Shared_Win, RadTransObsVolumeFrac_Shared
USE MOD_Radiation_Vars             ,ONLY: Radiation_Absorption_Spec_Shared_Win, RadiationInput
USE MOD_Radiation_Vars             ,ONLY: Radiation_Emission_Spec_Shared_Win, MacroRadInputParameters
USE MOD_Radiation_Vars             ,ONLY: Radiation_Absorption_SpecPercent_Shared_Win
USE MOD_Photon_TrackingVars        ,ONLY: PhotonSampWallProc,PhotonSampWall_Shared_Win_allocated
USE MOD_Photon_TrackingVars        ,ONLY: PhotonSampWall_Shared,PhotonSampWall_Shared_Win,PhotonSampWallProc
#else
USE MOD_Mesh_Vars                  ,ONLY: nElems
#endif
USE MOD_RayTracing_Vars            ,ONLY: Ray
USE MOD_Photon_Tracking            ,ONLY: InitPhotonSurfSample
USE MOD_Photon_TrackingVars        ,ONLY: PhotonSampWall
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iWave, iElem, firstElem, lastElem, ElemDisp, DisplRank, iSpec, currentRank
REAL                  :: LocTemp, ObsLengt, MaxSumTemp(2), GlobalMaxTemp(2), tmp
LOGICAL               :: ElemInCone
REAL,ALLOCATABLE      :: Radiation_ShockTube_Spec(:,:)
INTEGER               :: io_error
CHARACTER(LEN=3)      :: hilf
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT RADIATION TRANSPORT SOLVER ...'

ALLOCATE(RadiationElemAbsEnergy(2,1:nGlobalElems))
RadiationElemAbsEnergy=0.0
ALLOCATE(RadiationElemAbsEnergySpec(nSpecies,1:nGlobalElems))
RadiationElemAbsEnergySpec=0.0

RadiationDirectionModel = GETINT('Radiation-DirectionModel')
RadTrans%NumPhotonsPerCell = GETINT('Radiation-NumPhotonsPerCell')
RadiationAbsorptionModel = GETINT('Radiation-AbsorptionModel')
RadiationPhotonPosModel = GETINT('Radiation-PhotonPosModel')
RadiationPhotonWaveLengthModel = GETINT('Radiation-PhotonWaveLengthModel')
RadEmiAdaptPhotonNum = GETLOGICAL('Radiation-AdaptivePhotonNumEmission')
RadObservationPointMethod = GETINT('Radiation-RadObservationPointMethod')
ObservationDoConvolution = GETLOGICAL('Radiation-ObservationDoConvolution')
RadObservationPoint%ShockTubeDiameter = GETREAL('Radiation-ShockTubeDiameter')
WRITE(UNIT=hilf,FMT='(I0)') nSurfSample
Ray%nSurfSample    = GETINT('Radiation-nSurfSample',hilf)

! Build surface containers
CALL InitPhotonSurfSample()

IF (RadObservationPointMethod.GT.0) THEN
  RadObservationPoint%AngularAperture = GETREAL('Radiation-ObservationAngularAperture')
  RadObservationPoint%Diameter = GETREAL('Radiation-ObservationDiameter')
  RadObservationPoint%MidPoint = GETREALARRAY('Radiation-ObservationMidPoint',3)
  RadObservationPoint%ViewDirection = GETREALARRAY('Radiation-ObservationViewDirection',3)
  IF(.NOT.ALL(RadObservationPoint%ViewDirection(:).EQ.0.)) THEN
    RadObservationPoint%ViewDirection = RadObservationPoint%ViewDirection / VECNORM(RadObservationPoint%ViewDirection)
  END IF
  RadObservationPoint%SlitFunction = GETREALARRAY('Radiation-ObservationSlitFunction',2)
  IF(RadObservationPoint%SlitFunction(1).GT.RadObservationPoint%SlitFunction(2)) THEN
    tmp = RadObservationPoint%SlitFunction(1)
    RadObservationPoint%SlitFunction(1) = RadObservationPoint%SlitFunction(2)
    RadObservationPoint%SlitFunction(2) = tmp
  END IF
  RadObservationPoint%OrthoNormBasis(1:3,1) = RadObservationPoint%ViewDirection(1:3)
  CALL FindLinIndependentVectors(RadObservationPoint%OrthoNormBasis(1:3,1), RadObservationPoint%OrthoNormBasis(1:3,2), RadObservationPoint%OrthoNormBasis(1:3,3))
  CALL GramSchmidtAlgo(RadObservationPoint%OrthoNormBasis(1:3,1), RadObservationPoint%OrthoNormBasis(1:3,2), RadObservationPoint%OrthoNormBasis(1:3,3))
  IF (RadObservationPointMethod.EQ.2) RadObservationPoint%Diameter = 0.0
  ObsLengt = RadObservationPoint%Diameter/(2.*TAN(RadObservationPoint%AngularAperture/2.))
  RadObservationPoint%StartPoint(1:3) = RadObservationPoint%MidPoint(1:3) - ObsLengt*RadObservationPoint%ViewDirection(1:3)
  RadObservationPoint%Area = Pi*RadObservationPoint%Diameter*RadObservationPoint%Diameter/4.
  IF (RadObservationPointMethod.EQ.2) THEN
    RadObservationPoint%CalcFullSpectra = GETLOGICAL('Radiation-ObservationCalcFullSpectra')
    IF (RadObservationPoint%CalcFullSpectra) THEN
      RadEmiAdaptPhotonNum = .TRUE.
      RadTrans%NumPhotonsPerCell = RadiationParameter%WaveLenDiscrCoarse
    END IF
  END IF
END IF

IF(Symmetry%Order.EQ.2) CALL BuildMesh2DInfo()
IF (RadObservationPointMethod.GT.0) THEN
  ALLOCATE(RadObservation_Emission(RadiationParameter%WaveLenDiscrCoarse),RadObservation_Emission_Conv(RadiationParameter%WaveLenDiscrCoarse), &
    RadObservation_EmissionPart(RadiationParameter%WaveLenDiscrCoarse))
  RadObservation_Emission = 0.0
  RadObservation_Emission_Conv = 0.0
  RadObservation_EmissionPart = 0
END IF

#if USE_MPI
  ! allocate shared array for Radiation_Emission/Absorption_Spec
CALL Allocate_Shared((/2,nGlobalElems/),RadiationElemAbsEnergy_Shared_Win,RadiationElemAbsEnergy_Shared)
CALL MPI_WIN_LOCK_ALL(0,RadiationElemAbsEnergy_Shared_Win,IERROR)

IF (myComputeNodeRank.EQ.0) RadiationElemAbsEnergy_Shared = 0.
CALL BARRIER_AND_SYNC(RadiationElemAbsEnergy_Shared_Win,MPI_COMM_SHARED)

  ! allocate shared array for Radiation_Emission/Absorption_Spec
CALL Allocate_Shared((/nSpecies,nGlobalElems/),RadiationElemAbsEnergySpec_Shared_Win,RadiationElemAbsEnergySpec_Shared)
CALL MPI_WIN_LOCK_ALL(0,RadiationElemAbsEnergySpec_Shared_Win,IERROR)

IF (myComputeNodeRank.EQ.0) RadiationElemAbsEnergySpec_Shared = 0.
CALL BARRIER_AND_SYNC(RadiationElemAbsEnergySpec_Shared_Win,MPI_COMM_SHARED)

CALL Allocate_Shared((/nComputeNodeElems/), RadTransPhotPerCell_Shared_Win,RadTransPhotPerCell_Shared)
CALL MPI_WIN_LOCK_ALL(0,RadTransPhotPerCell_Shared_Win,IERROR)
CALL Allocate_Shared((/nComputeNodeElems/), Radiation_Emission_Spec_Total_Shared_Win,Radiation_Emission_Spec_Total_Shared)
CALL MPI_WIN_LOCK_ALL(0,Radiation_Emission_Spec_Total_Shared_Win,IERROR)

CALL Allocate_Shared((/nComputeNodeElems/), RadTransObsVolumeFrac_Shared_Win,RadTransObsVolumeFrac_Shared)
CALL MPI_WIN_LOCK_ALL(0,RadTransObsVolumeFrac_Shared_Win,IERROR)

RadTransPhotPerCell => RadTransPhotPerCell_Shared
Radiation_Emission_Spec_Total => Radiation_Emission_Spec_Total_Shared
RadTransObsVolumeFrac => RadTransObsVolumeFrac_Shared
IF (myComputeNodeRank.EQ.0) THEN
  RadTransPhotPerCell         = 0
  Radiation_Emission_Spec_Total = 0.0
  RadTransObsVolumeFrac     = 1.
END IF
CALL BARRIER_AND_SYNC(RadTransPhotPerCell_Shared_Win ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(RadTransObsVolumeFrac_Shared_Win ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(Radiation_Emission_Spec_Total_Shared_Win ,MPI_COMM_SHARED)

IF (RadiationPhotonWaveLengthModel.EQ.1) THEN
  CALL Allocate_Shared((/nComputeNodeElems/), Radiation_Emission_Spec_Max_Shared_Win,Radiation_Emission_Spec_Max_Shared)
  CALL MPI_WIN_LOCK_ALL(0,Radiation_Emission_Spec_Max_Shared_Win,IERROR)
  Radiation_Emission_Spec_Max => Radiation_Emission_Spec_Max_Shared
  IF (myComputeNodeRank.EQ.0) THEN
    Radiation_Emission_Spec_Max = 0.0
  END IF
  CALL BARRIER_AND_SYNC(Radiation_Emission_Spec_Max_Shared_Win ,MPI_COMM_SHARED)
END IF

ALLOCATE(RadTransPhotPerCellLoc(nComputeNodeElems))
RadTransPhotPerCellLoc = 0

IF (RadObservationPointMethod.EQ.2) THEN
  CALL Allocate_Shared((/7,nComputeNodeElems/), RadObservationPOI_Shared_Win,RadObservationPOI_Shared)
  CALL MPI_WIN_LOCK_ALL(0,RadObservationPOI_Shared_Win,IERROR)
  RadObservationPOI => RadObservationPOI_Shared
  IF (myComputeNodeRank.EQ.0) RadObservationPOI = 0.
  CALL BARRIER_AND_SYNC(RadObservationPOI_Shared_Win ,MPI_COMM_SHARED)
END IF
#else
! allocate local array for ElemInfo
ALLOCATE(RadTransPhotPerCell(nElems),Radiation_Emission_Spec_Total(nElems),RadTransPhotPerCellLoc(nELems), RadTransObsVolumeFrac(nElems))
RadTransPhotPerCell = 0
RadTransPhotPerCellLoc = 0
RadTransObsVolumeFrac = 1.0
Radiation_Emission_Spec_Total=0.0
IF (RadiationPhotonWaveLengthModel.EQ.1) THEN
  ALLOCATE(Radiation_Emission_Spec_Max(nElems))
  Radiation_Emission_Spec_Max=0.0
END IF
IF (RadObservationPointMethod.EQ.2) THEN
  ALLOCATE(RadObservationPOI(7,nElems))
  RadObservationPOI = 0.
END IF
#endif  /*USE_MPI*/


#if USE_MPI
IF (RadiationSwitches%MacroRadInput) THEN
  firstElem = INT(REAL( myComputeNodeRank   *nComputeNodeElems)/REAL(nComputeNodeProcessors))+1
  lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeElems)/REAL(nComputeNodeProcessors))
  IF (nComputeNodeElems.NE.nComputeNodeProcessors) THEN
    MaxSumTemp(1) = 0.0
    DO iSpec = 1, nSpecies
      IF(.NOT.RadiationInput(iSpec)%DoRadiation) CYCLE
      MaxSumTemp(1) = MaxSumTemp(1) + SUM(MacroRadInputParameters(firstElem:lastElem,iSpec,4))
    END DO
    CALL MPI_ALLREDUCE(MaxSumTemp(1), GlobalMaxTemp(1), 1, MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_SHARED,iError)
    GlobalMaxTemp(1) = GlobalMaxTemp(1) / REAL(nComputeNodeProcessors)
    MaxSumTemp(1) = 0.0
    currentRank = 0
    firstElem = 1
    ElemLoop: DO iElem = 1, nComputeNodeElems
      DO iSpec = 1, nSpecies
        IF(.NOT.RadiationInput(iSpec)%DoRadiation) CYCLE
        MaxSumTemp(1) = MaxSumTemp(1) + MacroRadInputParameters(iElem,iSpec,4)
      END DO
      IF ((nComputeNodeElems - iElem).EQ.(nComputeNodeProcessors - currentRank - 1)) THEN
        currentRank = currentRank + 1
        IF (currentRank.EQ.myComputeNodeRank) THEN
          firstElem = iElem + 1
          lastElem = iElem + 1
          EXIT ElemLoop
        ELSE
          CYCLE ElemLoop
        END IF
      END IF
      IF (MaxSumTemp(1).GE.GlobalMaxTemp(1)) THEN
        currentRank = currentRank + 1
        IF (currentRank.GT.myComputeNodeRank) THEN
          lastElem = iElem
          EXIT ElemLoop
        END IF
        IF (currentRank.EQ.myComputeNodeRank) firstElem = MIN(iElem+1, nComputeNodeElems)
        MaxSumTemp(1) = 0.0
      END IF
    END DO ElemLoop
    IF (myRank+1.EQ.nComputeNodeProcessors) lastElem = nComputeNodeElems
  END IF

  MaxSumTemp(2) = REAL(myRank)
  MaxSumTemp(1) = 0.0
  DO iSpec = 1, nSpecies
    IF(.NOT.RadiationInput(iSpec)%DoRadiation) CYCLE
    MaxSumTemp(1) = MaxSumTemp(1) + SUM(MacroRadInputParameters(firstElem:lastElem,iSpec,4))
  END DO
  CALL MPI_ALLREDUCE(MaxSumTemp, GlobalMaxTemp, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC,MPI_COMM_PICLAS,iError)
  DisplRank = NINT(GlobalMaxTemp(2))
ELSE
  firstElem = INT(REAL( myComputeNodeRank   *nComputeNodeElems)/REAL(nComputeNodeProcessors))+1
  lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeElems)/REAL(nComputeNodeProcessors))
  DisplRank = 0
END IF
#else
  firstElem = 1
  lastElem  = nElems
  DisplRank = 0
#endif
SELECT CASE(RadiationSwitches%RadType)
CASE(1) !calls radition solver module
  SWRITE(UNIT_stdOut,'(A)') ' Calculate Radiation Data per Cell ...'
  ElemDisp = INT((lastElem-firstElem+1)/100)
  ElemDisp = MAX(1,ElemDisp)

  DO iElem = firstElem, lastElem
    IF((myRank.EQ.DisplRank).AND.(MOD(iElem-firstElem,ElemDisp).EQ.0)) CALL PrintStatusLineRadiation(REAL(iElem),REAL(firstElem),REAL(lastElem),.FALSE.,DisplRank)
    IF (RadObservationPointMethod.EQ.1) THEN
      CALL ElemInObsCone(iElem, ElemInCone)
      IF (.NOT.ElemInCone) CYCLE
    ELSE IF (RadObservationPointMethod.EQ.2) THEN
      CALL ElemOnLineOfSight(iELem, ElemInCone)
      IF (.NOT.ElemInCone) CYCLE
    END IF
    CALL radiation_main(iElem)
    DO iWave = 1, RadiationParameter%WaveLenDiscrCoarse
      Radiation_Emission_Spec_Total(iElem) = Radiation_Emission_Spec_Total(iElem) &
          + 4.*Pi*Radiation_Emission_Spec(iWave, iElem) * RadiationParameter%WaveLenIncr*RadiationParameter%WaveLenReductionFactor
      IF (RadiationPhotonWaveLengthModel.EQ.1) Radiation_Emission_Spec_Max(iElem) = MAX(Radiation_Emission_Spec_Max(iElem),  &
        4.*Pi*Radiation_Emission_Spec(iWave, iElem) * RadiationParameter%WaveLenIncr*RadiationParameter%WaveLenReductionFactor)
    END DO
    IF (RadiationParameter%WaveLenReductionFactor.GT.1) THEN
      IF (MOD(RadiationParameter%WaveLenDiscr,RadiationParameter%WaveLenDiscrCoarse).NE.0) THEN
        Radiation_Emission_Spec_Total(iElem) = Radiation_Emission_Spec_Total(iElem) &
          + 4.*Pi*Radiation_Emission_Spec(RadiationParameter%WaveLenDiscrCoarse, iElem) * RadiationParameter%WaveLenIncr
        IF (RadiationPhotonWaveLengthModel.EQ.1) Radiation_Emission_Spec_Max(iElem) = MAX(Radiation_Emission_Spec_Max(iElem),  &
          4.*Pi*Radiation_Emission_Spec(RadiationParameter%WaveLenDiscrCoarse, iElem) * RadiationParameter%WaveLenIncr*(RadiationParameter%WaveLenReductionFactor+1.))
      END IF
    END IF
  END DO
CASE(2) ! Black body radiation

  DO iElem = firstElem, lastElem
    IF (ElemMidPoint_Shared(1,iElem).LT.1) THEN
      LocTemp = 10000. !GEO%ElemMidPoint(2,iElem)/GEO%ymaxglob*10000.
    ELSE
      LocTemp = 10000. !GEO%ElemMidPoint(2,iElem)/GEO%ymaxglob*10000
    END IF
    DO iWave = 1, RadiationParameter%WaveLenDiscr
      IF (LocTemp.GT.0.0) Radiation_Emission_Spec(iWave, iElem) = 2.*PlanckConst*c*c/(RadiationParameter%WaveLen(iWave)**5. &
          *(EXP(PlanckConst*c/(RadiationParameter%WaveLen(iWave)*BoltzmannConst*LocTemp))-1.) )
      Radiation_Emission_Spec_Total(iElem) = Radiation_Emission_Spec_Total(iElem) &
          + 4.*Pi*Radiation_Emission_Spec(iWave, iElem) * RadiationParameter%WaveLenIncr
    END DO
  END DO
  DO iElem = firstElem, lastElem
    DO iWave = 1, RadiationParameter%WaveLenDiscr
      Radiation_Absorption_Spec(iWave, GetGlobalElemID(iElem)) = 1.
      Radiation_Absorption_SpecPercent(iWave,:,GetGlobalElemID(iElem)) = 10000
    END DO
  END DO
CASE(3) !only radiation
  SWRITE(UNIT_stdOut,'(A)') ' Calculate Radiation Data per Cell ...'
  DO iElem = firstElem, lastElem
    IF(MPIroot.AND.(MOD(iElem,10).EQ.0)) CALL PrintStatusLineRadiation(REAL(iElem),REAL(firstElem),REAL(lastElem),.FALSE.)
    CALL radiation_main(iElem)
    DO iWave = 1, RadiationParameter%WaveLenDiscr
      Radiation_Emission_Spec_Total(iElem) = Radiation_Emission_Spec_Total(iElem) &
          + 4.*Pi*Radiation_Emission_Spec(iWave, iElem) * RadiationParameter%WaveLenIncr
    END DO
  END DO
CASE(4) !Shocktube mode
  ALLOCATE(Radiation_ShockTube_Spec(RadiationParameter%WaveLenDiscr,nGlobalElems))
  SWRITE(UNIT_stdOut,'(A)') ' Calculate Radiation Data per Cell ...'
  DO iElem = firstElem, lastElem
    IF(MPIroot.AND.(MOD(iElem,10).EQ.0)) CALL PrintStatusLineRadiation(REAL(iElem),REAL(firstElem),REAL(lastElem),.FALSE.)
    CALL radiation_main(iElem)
    DO iWave = 1, RadiationParameter%WaveLenDiscr
      Radiation_Emission_Spec_Total(iElem) = Radiation_Emission_Spec_Total(iElem) &
          + 4.*Pi*Radiation_Emission_Spec(iWave, iElem) * RadiationParameter%WaveLenIncr
      IF(Radiation_Absorption_Spec(iWave, iElem).EQ.0.0) THEN
        Radiation_ShockTube_Spec(iWave,iElem) = Radiation_Emission_Spec(iWave, iElem)*RadObservationPoint%ShockTubeDiameter
      ELSE
        Radiation_ShockTube_Spec(iWave,iElem) = Radiation_Emission_Spec(iWave, iElem)/Radiation_Absorption_Spec(iWave, iElem) * &
          (1.-EXP(-Radiation_Absorption_Spec(iWave, iElem)*RadObservationPoint%ShockTubeDiameter))
      END IF
    END DO
  END DO

  OPEN(unit=40,file='Radiation_Shocktube.csv',status='replace',action='write', iostat=io_error)
    DO iElem=1,nGlobalElems
      WRITE(40,CSVFORMAT,ADVANCE="NO") ',', ElemMidPoint_Shared(1,iElem)
    END DO
    WRITE(40,*)
  DO iWave =1,RadiationParameter%WaveLenDiscr
    WRITE(40,'(E23.16E3)',ADVANCE="NO") RadiationParameter%WaveLen(iWave)*1.E9
    DO iElem = 1,nGlobalElems
      WRITE(40,CSVFORMAT,ADVANCE="NO") ',', Radiation_ShockTube_Spec(iWave,iElem)
    END DO
    WRITE(40,*)
  END DO

CASE DEFAULT
  CALL abort(__STAMP__,' ERROR: Radiation type is not implemented! (unknown case)')
END SELECT


#if USE_MPI
  CALL BARRIER_AND_SYNC(RadTransObsVolumeFrac_Shared_Win ,MPI_COMM_SHARED)
  IF (RadiationPhotonWaveLengthModel.EQ.1) THEN
    CALL BARRIER_AND_SYNC(Radiation_Emission_Spec_Max_Shared_Win,MPI_COMM_SHARED)
  END IF
  CALL BARRIER_AND_SYNC(Radiation_Emission_Spec_Total_Shared_Win,MPI_COMM_SHARED)
  CALL BARRIER_AND_SYNC(Radiation_Emission_Spec_Shared_Win ,MPI_COMM_SHARED)
  CALL BARRIER_AND_SYNC(Radiation_Absorption_Spec_Shared_Win ,MPI_COMM_SHARED)
  IF(nLeaderGroupProcs.GT.1)THEN
    IF(myComputeNodeRank.EQ.0)THEN
      CALL MPI_ALLGATHERV( MPI_IN_PLACE                                     &
                     , 0                                                    &
                     , MPI_DATATYPE_NULL                                    &
                     , Radiation_Absorption_Spec                            &
                     , RadiationParameter%WaveLenDiscrCoarse *recvcountElem &
                     , RadiationParameter%WaveLenDiscrCoarse *displsElem    &
                     , MPI_DOUBLE_PRECISION                                 &
                     , MPI_COMM_LEADERS_SHARED                              &
                     , IERROR)
    END IF
  END IF
  CALL BARRIER_AND_SYNC(Radiation_Absorption_Spec_Shared_Win ,MPI_COMM_SHARED)
  CALL BARRIER_AND_SYNC(Radiation_Absorption_SpecPercent_Shared_Win ,MPI_COMM_SHARED)
  IF(nLeaderGroupProcs.GT.1)THEN
    IF(myComputeNodeRank.EQ.0)THEN
      CALL MPI_ALLGATHERV( MPI_IN_PLACE                                              &
                     , 0                                                             &
                     , MPI_DATATYPE_NULL                                             &
                     , Radiation_Absorption_SpecPercent                              &
                     , RadiationParameter%WaveLenDiscrCoarse *nSpecies*recvcountElem &
                     , RadiationParameter%WaveLenDiscrCoarse *nSpecies*displsElem    &
                     , MPI_INTEGER2                                                  &
                     , MPI_COMM_LEADERS_SHARED                                       &
                     , IERROR)
    END IF
  END IF
  CALL BARRIER_AND_SYNC(Radiation_Absorption_SpecPercent_Shared_Win ,MPI_COMM_SHARED)
  IF (RadObservationPointMethod.EQ.2) CALL BARRIER_AND_SYNC(RadObservationPOI_Shared_Win ,MPI_COMM_SHARED)
  !print*, 'AHAAAA', SUM(RadObservationPOI(7,:))
  !read*
#endif
  RadTrans%GlobalRadiationPower = 0.0
  RadTrans%ScaledGlobalRadiationPower = 0.0
  DO iElem = firstElem, lastElem
    RadTrans%GlobalRadiationPower = RadTrans%GlobalRadiationPower + Radiation_Emission_Spec_Total(iElem)*ElemVolume_Shared(iElem)*RadTransObsVolumeFrac(iElem)
    IF (RadialWeighting%DoRadialWeighting) THEN
      RadTrans%ScaledGlobalRadiationPower = RadTrans%ScaledGlobalRadiationPower  &
        + Radiation_Emission_Spec_Total(iElem)*ElemVolume_Shared(iElem)*RadTransObsVolumeFrac(iElem) &
        /(1. + ElemMidPoint_Shared(2,iElem)/GEO%ymaxglob*(RadialWeighting%PartScaleFactor-1.))
    END IF
  END DO
#if USE_MPI
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,RadTrans%GlobalRadiationPower,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_PICLAS,iError)
  IF (RadialWeighting%DoRadialWeighting) THEN
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,RadTrans%ScaledGlobalRadiationPower,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_PICLAS,iError)
  END IF
#endif /*USE_MPI*/
  RadTrans%GlobalPhotonNum = RadTrans%NumPhotonsPerCell * nGlobalElems




#if USE_MPI
ALLOCATE(PhotonSampWallProc(2,1:Ray%nSurfSample,1:Ray%nSurfSample,1:nComputeNodeSurfTotalSides))
PhotonSampWallProc=0.0
!> Then shared arrays for boundary sampling
CALL Allocate_Shared((/2,Ray%nSurfSample,Ray%nSurfSample,nComputeNodeSurfTotalSides/),PhotonSampWall_Shared_Win,PhotonSampWall_Shared)
PhotonSampWall_Shared_Win_allocated = .TRUE.
CALL MPI_WIN_LOCK_ALL(0,PhotonSampWall_Shared_Win,IERROR)
PhotonSampWall => PhotonSampWall_Shared

IF (myComputeNodeRank.EQ.0) PhotonSampWall = 0.
CALL BARRIER_AND_SYNC(PhotonSampWall_Shared_Win,MPI_COMM_SHARED)
#else
ALLOCATE(PhotonSampWall(2,1:Ray%nSurfSample,1:Ray%nSurfSample,1:nComputeNodeSurfTotalSides))
PhotonSampWall=0.0
#endif

SWRITE(UNIT_stdOut,'(A)')' INIT RADIATION TRANSPORT SOLVER DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitRadiationTransport


SUBROUTINE HALTON( ind, dims, rand )
!===================================================================================================================================
! Halton sequence for reducing stochastical noise in radiative transfer
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INOUTPUT VARIABLES
INTEGER, INTENT(IN)             :: ind, dims
REAL, INTENT(OUT)               :: rand(dims)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: t(dims),j, i1, d
REAL                          :: primeinv(dims)
!===================================================================================================================================
  t(1:dims) = ABS(ind)
  DO i1 = 1, dims
    primeinv(i1) = 1.0/REAL(PRIME(i1))
  END DO
  rand  = 0.0

  DO WHILE (ANY(t(1:dims).NE.0))
    do j = 1, dims
      d = MOD(t(j), PRIME(j))
      rand(j) = rand(j) + REAL(d) * primeinv(j)
      primeinv(j) = primeinv(j) / REAL(PRIME( j ))
      t(j) = ( t(j) / PRIME ( j ) )
    END DO
  END DO

  RETURN
END SUBROUTINE HALTON

SUBROUTINE ElemInObsCone(ElemID, ElemInCone)
!===================================================================================================================================
! Routine to check if element is in opening cone of observation angle of radiative transfer
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Tools                ,ONLY: GetGlobalElemID
USE MOD_Particle_Mesh_Vars        ,ONLY: NodeCoords_Shared,ElemInfo_Shared, BoundsOfElem_Shared, SideIsSymSide, SideInfo_Shared
USE MOD_RadiationTrans_Vars       ,ONLY: RadObservationPoint, RadTransObsVolumeFrac
USE MOD_Particle_Vars             ,ONLY: Symmetry
USE MOD_Particle_Mesh_Tools       ,ONLY: ParticleInsideQuad3D
USE MOD_Photon_TrackingTools      ,ONLY: PhotonIntersectionWithSide2DDir, PhotonThroughSideCheck3DDir
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INOUTPUT VARIABLES
INTEGER, INTENT(IN)             :: ElemID
LOGICAL, INTENT(OUT)            :: ElemInCone
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iNode, MCVar, iGlobalElem, iPoint, iLocSide, nlocSides, TempSideID, localSideID, TriNum
LOGICAL                       :: NodeInCone(8), InsideFlag, ThroughSide
REAL                          :: NodePoint(3), ConeDist, ConeRadius, orthoDist, RandomPos(3)
!===================================================================================================================================
ElemInCone = .FALSE.
NodeInCone = .FALSE.
MCVar = 1000000
DO iNode = 1, 8
  NodePoint(1:3) = NodeCoords_Shared(1:3,ElemInfo_Shared(ELEM_FIRSTNODEIND,ElemID)+iNode)
  ConeDist = DOT_PRODUCT(NodePoint(1:3) - RadObservationPoint%StartPoint(1:3), RadObservationPoint%ViewDirection(1:3))
  ConeRadius = TAN(RadObservationPoint%AngularAperture/2.) * ConeDist
  orthoDist = VECNORM(NodePoint(1:3) - RadObservationPoint%StartPoint(1:3) - ConeDist*RadObservationPoint%ViewDirection(1:3))
  IF (orthoDist.LE.ConeRadius) THEN
    NodeInCone(iNode) = .TRUE.
  END IF
END DO

IF (ALL(NodeInCone)) THEN
  ElemInCone = .TRUE.
ELSE IF (ANY(NodeInCone)) THEN
  iGlobalElem = GetGlobalElemID(ElemID)
  RadTransObsVolumeFrac(ElemID) = 0.0
  ElemInCone = .TRUE.
  ASSOCIATE( Bounds => BoundsOfElem_Shared(1:2,1:3,iGlobalElem) )
    DO iPoint = 1, MCVar
      InsideFlag=.FALSE.
      DO WHILE(.NOT.InsideFlag)
        CALL RANDOM_NUMBER(RandomPos)
        RandomPos = Bounds(1,:) + RandomPos*(Bounds(2,:)-Bounds(1,:))
        IF(Symmetry%Order.LE.2) RandomPos(3) = 0.
        IF(Symmetry%Order.LE.1) RandomPos(2) = 0.
        CALL ParticleInsideQuad3D(RandomPos,iGlobalElem,InsideFlag)
      END DO
      ConeDist = DOT_PRODUCT(RandomPos(1:3) - RadObservationPoint%StartPoint(1:3), RadObservationPoint%ViewDirection(1:3))
      ConeRadius = TAN(RadObservationPoint%AngularAperture/2.) * ConeDist
      orthoDist = VECNORM(RandomPos(1:3) - RadObservationPoint%StartPoint(1:3) - ConeDist*RadObservationPoint%ViewDirection(1:3))
      IF (orthoDist.LE.ConeRadius)  RadTransObsVolumeFrac(ElemID) = RadTransObsVolumeFrac(ElemID) + 1./REAL(MCVar)
    END DO
  END ASSOCIATE
ELSE
  nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)
  SideLoop: DO iLocSide=1,nlocSides
    TempSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide
    localSideID = SideInfo_Shared(SIDE_LOCALID,TempSideID)
    ! Side is not one of the 6 local sides
    IF (localSideID.LE.0) CYCLE
    IF(Symmetry%Axisymmetric) THEN
      IF (SideIsSymSide(TempSideID)) CYCLE
      ThroughSide = .FALSE.
      CALL PhotonIntersectionWithSide2DDir(localSideID,ElemID,ThroughSide, RadObservationPoint%StartPoint(1:3),RadObservationPoint%ViewDirection(1:3))
      IF (ThroughSide) THEN
        ElemInCone = .TRUE.
        ASSOCIATE( Bounds => BoundsOfElem_Shared(1:2,1:3,iGlobalElem) )
          DO iPoint = 1, MCVar
            InsideFlag=.FALSE.
            DO WHILE(.NOT.InsideFlag)
              CALL RANDOM_NUMBER(RandomPos)
              RandomPos = Bounds(1,:) + RandomPos*(Bounds(2,:)-Bounds(1,:))
              IF(Symmetry%Order.LE.2) RandomPos(3) = 0.
              IF(Symmetry%Order.LE.1) RandomPos(2) = 0.
              CALL ParticleInsideQuad3D(RandomPos,iGlobalElem,InsideFlag)
            END DO
            ConeDist = DOT_PRODUCT(RandomPos(1:3) - RadObservationPoint%StartPoint(1:3), RadObservationPoint%ViewDirection(1:3))
            ConeRadius = TAN(RadObservationPoint%AngularAperture/2.) * ConeDist
            orthoDist = VECNORM(RandomPos(1:3) - RadObservationPoint%StartPoint(1:3) - ConeDist*RadObservationPoint%ViewDirection(1:3))
            IF (orthoDist.LE.ConeRadius)  RadTransObsVolumeFrac(ElemID) = RadTransObsVolumeFrac(ElemID) + 1./REAL(MCVar)
          END DO
        END ASSOCIATE
        EXIT SideLoop
      END IF
    ELSE
      DO TriNum = 1,2
        ThroughSide = .FALSE.
        CALL PhotonThroughSideCheck3DDir(localSideID,ElemID,ThroughSide,TriNum, RadObservationPoint%StartPoint(1:3),RadObservationPoint%ViewDirection(1:3))
        IF (ThroughSide) THEN
          ElemInCone = .TRUE.
          ASSOCIATE( Bounds => BoundsOfElem_Shared(1:2,1:3,iGlobalElem) )
            DO iPoint = 1, MCVar
              InsideFlag=.FALSE.
              DO WHILE(.NOT.InsideFlag)
                CALL RANDOM_NUMBER(RandomPos)
                RandomPos = Bounds(1,:) + RandomPos*(Bounds(2,:)-Bounds(1,:))
                IF(Symmetry%Order.LE.2) RandomPos(3) = 0.
                IF(Symmetry%Order.LE.1) RandomPos(2) = 0.
                CALL ParticleInsideQuad3D(RandomPos,iGlobalElem,InsideFlag)
              END DO
              ConeDist = DOT_PRODUCT(RandomPos(1:3) - RadObservationPoint%StartPoint(1:3), RadObservationPoint%ViewDirection(1:3))
              ConeRadius = TAN(RadObservationPoint%AngularAperture/2.) * ConeDist
              orthoDist = VECNORM(RandomPos(1:3) - RadObservationPoint%StartPoint(1:3) - ConeDist*RadObservationPoint%ViewDirection(1:3))
              IF (orthoDist.LE.ConeRadius)  RadTransObsVolumeFrac(ElemID) = RadTransObsVolumeFrac(ElemID) + 1./REAL(MCVar)
            END DO
          END ASSOCIATE
          EXIT SideLoop
        END IF
      END DO
    END IF
  END DO SideLoop
END IF

END SUBROUTINE ElemInObsCone


SUBROUTINE ElemOnLineOfSight(ElemID, ElemInCone)
!===================================================================================================================================
! Routine to check if element is on line of sight of radiative transfer
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Tools                ,ONLY: GetGlobalElemID
USE MOD_Particle_Mesh_Vars        ,ONLY: ElemInfo_Shared, SideInfo_Shared, SideIsSymSide
USE MOD_RadiationTrans_Vars       ,ONLY: RadObservationPoint, RadObservationPOI
USE MOD_Particle_Vars             ,ONLY: Symmetry
USE MOD_Particle_Mesh_Tools       ,ONLY: ParticleInsideQuad3D
USE MOD_Photon_TrackingTools      ,ONLY: PhotonIntersectionWithSide2DDir, PhotonThroughSideCheck3DDir
USE MOD_Particle_Boundary_Vars    ,ONLY: PartBound
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INOUTPUT VARIABLES
INTEGER, INTENT(IN)             :: ElemID
LOGICAL, INTENT(OUT)            :: ElemInCone
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iLocSide, nlocSides, TempSideID, localSideID, TriNum
INTEGER                       :: nThroughSide, BCType
LOGICAL                       :: ThroughSide, IsSymElem
REAL                          :: IntersectionPos(1:3), Distance(2)
REAL                          :: length
!===================================================================================================================================
ElemInCone = .FALSE.
IsSymElem = .FALSE.
nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,ElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID)
nThroughSide = 0
SideLoop: DO iLocSide=1,nlocSides
  TempSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,ElemID) + iLocSide
  localSideID = SideInfo_Shared(SIDE_LOCALID,TempSideID)
  ! Side is not one of the 6 local sides
  IF (localSideID.LE.0) CYCLE
  IF(Symmetry%Axisymmetric) THEN
    IF (SideIsSymSide(TempSideID)) CYCLE
    IF (SideInfo_Shared(SIDE_BCID,TempSideID).GT.0) THEN
      BCType = PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,TempSideID)))
      IF (BCType.EQ.PartBound%SymmetryAxis) IsSymElem = .TRUE.
    END IF
    ThroughSide = .FALSE.
    CALL PhotonIntersectionWithSide2DDir(localSideID,ElemID,ThroughSide, RadObservationPoint%StartPoint(1:3),&
        RadObservationPoint%ViewDirection(1:3), IntersectionPos(1:3), Distance(nThroughSide+1))
    IF (ThroughSide) THEN
      ElemInCone = .TRUE.
      RadObservationPOI(1+nThroughSide*3:3+nThroughSide*3, ElemID) = IntersectionPos(1:3)
      nThroughSide = nThroughSide + 1
      IF (nThroughSide.EQ.2) THEN
        IF (Distance(2).LT.Distance(1)) THEN
          IntersectionPos(1:3) = RadObservationPOI(1:3, ElemID)
          RadObservationPOI(1:3, ElemID) = RadObservationPOI(4:6, ElemID)
          RadObservationPOI(4:6, ElemID) = IntersectionPos(1:3)
        END IF
        RadObservationPOI(7, ElemID) = VECNORM(RadObservationPOI(4:6, ElemID)-RadObservationPOI(1:3, ElemID))
        EXIT SideLoop
      END IF
    END IF
  ELSE
    DO TriNum = 1,2
      ThroughSide = .FALSE.
      CALL PhotonThroughSideCheck3DDir(localSideID,ElemID,ThroughSide,TriNum, RadObservationPoint%StartPoint(1:3),RadObservationPoint%ViewDirection(1:3))
      IF (ThroughSide) THEN
        ElemInCone = .TRUE.
        EXIT SideLoop
      END IF
    END DO
  END IF
END DO SideLoop
IF (ElemInCone.AND.(nThroughSide.NE.2)) THEN
  IF (IsSymElem) THEN
    IF (nThroughSide.NE.1) THEN
      CALL abort(&
      __STAMP__&
      ,' Cannot find 1 POI of LOS in Elem', ElemID)
    END IF
    RadObservationPOI(4:6, ElemID) = RadObservationPOI(1:3, ElemID)
    length = -RadObservationPoint%StartPoint(2)/RadObservationPoint%ViewDirection(2)
    RadObservationPOI(2:3, ElemID) = 0.0
    RadObservationPOI(1, ElemID) = RadObservationPoint%StartPoint(1) + length*RadObservationPoint%ViewDirection(1)
    RadObservationPOI(7, ElemID) = VECNORM(RadObservationPOI(4:6, ElemID)-RadObservationPOI(1:3, ElemID))
    RETURN
  END IF
  CALL abort(&
      __STAMP__&
      ,' Cannot find POI of LOS in Elem', ElemID)
END IF

END SUBROUTINE ElemOnLineOfSight

INTEGER FUNCTION PRIME(n)
!===================================================================================================================================
! contains prime numbers for initialization of Halton sequence
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INOUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, PARAMETER              :: prime_max=50
INTEGER, SAVE                   :: icall = 0, npvec(prime_max)
INTEGER, INTENT(IN)             :: n
!===================================================================================================================================
  IF (icall.EQ.0) THEN
    icall = 1
    npvec(1:50) = (/ &
        2,    3,    5,    7,   11,   13,   17,   19,   23,   29, &
       31,   37,   41,   43,   47,   53,   59,   61,   67,   71, &
       73,   79,   83,   89,   97,  101,  103,  107,  109,  113, &
      127,  131,  137,  139,  149,  151,  157,  163,  167,  173, &
      179,  181,  191,  193,  197,  199,  211,  223,  227,  229 /)
  END IF
  IF (n.EQ.-1) THEN
    PRIME = prime_max
  ELSE IF (n.EQ.0) THEN
    PRIME = 1
  ELSE IF (n.LE.prime_max) THEN
    PRIME = npvec(n)
  ELSE
    PRIME = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PRIME - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal prime index N = ', n
    write ( *, '(a,i8)' ) '  N should be between 1 and PRIME_MAX =', prime_max
    stop 1
  END IF

  RETURN
END FUNCTION PRIME


SUBROUTINE FinalizeRadiationTransport()
!===================================================================================================================================
!> Deallocating radiation variables
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_RadiationTrans_Vars
#if USE_MPI
USE MOD_MPI_Shared_Vars    ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared
USE MOD_Particle_Mesh_Vars ,ONLY: ElemSideNodeID2D_Shared_Win,SideNormalEdge2D_Shared_Win
#endif
USE MOD_Particle_Vars      ,ONLY: Symmetry
USE MOD_Particle_Mesh_Vars ,ONLY: ElemSideNodeID2D_Shared,SideNormalEdge2D_Shared
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

#if USE_MPI
! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
CALL UNLOCK_AND_FREE(RadiationElemAbsEnergy_Shared_Win)
CALL UNLOCK_AND_FREE(RadiationElemAbsEnergySpec_Shared_Win)
CALL UNLOCK_AND_FREE(RadTransPhotPerCell_Shared_Win)
CALL UNLOCK_AND_FREE(Radiation_Emission_Spec_Total_Shared_Win)
CALL UNLOCK_AND_FREE(RadTransObsVolumeFrac_Shared_Win)
IF(RadiationPhotonWaveLengthModel.EQ.1) CALL UNLOCK_AND_FREE(Radiation_Emission_Spec_Max_Shared_Win)
IF(Symmetry%Order.EQ.2)THEN
  CALL UNLOCK_AND_FREE(ElemSideNodeID2D_Shared_Win)
  CALL UNLOCK_AND_FREE(SideNormalEdge2D_Shared_Win)
END IF
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
#endif /*USE_MPI*/
ADEALLOCATE(RadiationElemAbsEnergy_Shared)
ADEALLOCATE(RadiationElemAbsEnergySpec_Shared)
ADEALLOCATE(RadTransPhotPerCell_Shared)
ADEALLOCATE(Radiation_Emission_Spec_Total_Shared)
ADEALLOCATE(RadTransObsVolumeFrac_Shared)
ADEALLOCATE(Radiation_Emission_Spec_Max_Shared)
ADEALLOCATE(ElemSideNodeID2D_Shared)
ADEALLOCATE(SideNormalEdge2D_Shared)

END SUBROUTINE FinalizeRadiationTransport

END MODULE MOD_RadiationTrans_Init
