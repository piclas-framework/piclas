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

MODULE MOD_DSMC_Analyze
!===================================================================================================================================
! Module for DSMC Sampling and Output
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE WriteDSMCToHDF5
  MODULE PROCEDURE WriteDSMCToHDF5
END INTERFACE

INTERFACE CalcTVib
  MODULE PROCEDURE CalcTVib
END INTERFACE

INTERFACE CalcSurfaceValues
  MODULE PROCEDURE CalcSurfaceValues
END INTERFACE

INTERFACE CalcTelec
  MODULE PROCEDURE CalcTelec
END INTERFACE

INTERFACE CalcTVibPoly
  MODULE PROCEDURE CalcTVibPoly
END INTERFACE

INTERFACE CalcMeanFreePath
  MODULE PROCEDURE CalcMeanFreePath
END INTERFACE

INTERFACE CalcGammaVib
  MODULE PROCEDURE CalcGammaVib
END INTERFACE

INTERFACE CalcInstantTransTemp
  MODULE PROCEDURE CalcInstantTransTemp
END INTERFACE

INTERFACE SummarizeQualityFactors
  MODULE PROCEDURE SummarizeQualityFactors
END INTERFACE

INTERFACE DSMCMacroSampling
  MODULE PROCEDURE DSMCMacroSampling
END INTERFACE

INTERFACE SamplingRotVibRelaxProb
  MODULE PROCEDURE SamplingRotVibRelaxProb
END INTERFACE

INTERFACE CalcInstantElecTempXi
  MODULE PROCEDURE CalcInstantElecTempXi
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: DSMC_data_sampling, CalcMeanFreePath,WriteDSMCToHDF5
PUBLIC :: CalcTVib, CalcSurfaceValues, CalcTelec, CalcTVibPoly, CalcGammaVib
PUBLIC :: CalcInstantTransTemp, SummarizeQualityFactors, DSMCMacroSampling
PUBLIC :: SamplingRotVibRelaxProb, CalcInstantElecTempXi
!===================================================================================================================================

CONTAINS

SUBROUTINE CalcSurfaceValues(during_dt_opt)
!===================================================================================================================================
!> Calculates macroscopic surface values from samples
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars                  ,ONLY: MacroSurfaceVal,DSMC,MacroSurfaceSpecVal
USE MOD_Mesh_Vars                  ,ONLY: MeshFile,BC
USE MOD_Particle_Boundary_Sampling ,ONLY: WriteSurfSampleToHDF5
USE MOD_Particle_Boundary_Vars     ,ONLY: SurfMesh,nSurfSample,SampWall,CalcSurfCollis,nPorousBC,PartBound,CalcSurfaceImpact
USE MOD_Particle_Boundary_Vars     ,ONLY: SurfSide2GlobalSide, GlobalSide2SurfSide
USE MOD_Particle_Boundary_Vars     ,ONLY: nComputeNodeSurfSides
USE MOD_Particle_Boundary_vars     ,ONLY: SampWallState_Shared,SampWallImpactNumber_Shared,SampWallImpactEnergy_Shared
USE MOD_Particle_Boundary_vars     ,ONLY: SampWallImpactVector_Shared,SampWallImpactAngle_Shared
USE MOD_Particle_Boundary_vars     ,ONLY: SurfSideArea_Shared
USE MOD_Particle_Boundary_Vars     ,ONLY: PorousBCInfo_Shared,MapSurfSideToPorousSide_Shared,SampWallPumpCapacity_Shared
USE MOD_Particle_Mesh_Vars         ,ONLY: SideInfo_Shared
USE MOD_Particle_Vars              ,ONLY: WriteMacroSurfaceValues,nSpecies,MacroValSampTime,VarTimeStep,Symmetry2D
USE MOD_Restart_Vars               ,ONLY: RestartTime
USE MOD_SurfaceModel_Vars          ,ONLY: Adsorption
USE MOD_TimeDisc_Vars              ,ONLY: TEnd
USE MOD_Timedisc_Vars              ,ONLY: time,dt
#if USE_MPI
USE MOD_MPI_Shared_Vars            ,ONLY: MPI_COMM_LEADERS_SURF
USE MOD_Particle_Boundary_Vars     ,ONLY: SurfOnNode
USE MOD_Particle_MPI_Boundary_Sampling,ONLY: ExchangeSurfData
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL, INTENT(IN), OPTIONAL      :: during_dt_opt !routine was called during timestep (i.e. before iter=iter+1, time=time+dt...)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: BCID
INTEGER                            :: iSpec,iSurfSide,p,q, iReact, nVar, nVarSpec, iPBC, nVarCount
INTEGER                            :: nAdsSamples, iAdsSampl
REAL                               :: TimeSample, ActualTime, TimeSampleTemp, CounterSum
INTEGER, ALLOCATABLE               :: CounterTotal(:), SumCounterTotal(:)              ! Total Wall-Collision counter
LOGICAL                            :: during_dt
INTEGER                            :: idx,OutputCounter
INTEGER                            :: GlobalSideID, SurfSideNb
!===================================================================================================================================

IF (PRESENT(during_dt_opt)) THEN
  during_dt=during_dt_opt
ELSE
  during_dt=.FALSE.
END IF
IF (during_dt) THEN
  ActualTime=time+dt
ELSE
  ActualTime=time
END IF

IF (WriteMacroSurfaceValues) THEN
  TimeSample = Time - MacroValSampTime !elapsed time since last sampling (variable dt's possible!)
  MacroValSampTime = Time
ELSE IF (RestartTime.GT.(1-DSMC%TimeFracSamp)*TEnd) THEN
  TimeSample = Time - RestartTime
ELSE
  TimeSample = (Time-(1-DSMC%TimeFracSamp)*TEnd)
END IF

IF(ALMOSTZERO(TimeSample)) RETURN

IF (CalcSurfCollis%AnalyzeSurfCollis) THEN
  CALL WriteAnalyzeSurfCollisToHDF5(ActualTime,TimeSample)
END IF

IF(.NOT.SurfOnNode) RETURN

#if USE_MPI
CALL ExchangeSurfData()
#endif

! Only surface sampling leaders take part in the remainder of this routine
IF (MPI_COMM_LEADERS_SURF.EQ.MPI_COMM_NULL) RETURN

! Determine the number of variables
nVar = 5
nVarSpec = 1
IF(ANY(PartBound%Reactive)) THEN
  nAdsSamples = 5
  nVar        = nVar + nAdsSamples
  nVarSpec    = nVarSpec + 2 + 2*Adsorption%ReactNum
END IF

! Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z), angle and number: Add 8 to the buffer length
nVarSpec = nVarSpec + 8

IF(nPorousBC.GT.0) THEN
  nVar = nVar + nPorousBC
END IF
! Allocate the output container
ALLOCATE(MacroSurfaceVal(1:nVar         , 1:nSurfSample , 1:nSurfSample , nComputeNodeSurfSides))
MacroSurfaceVal     = 0.
ALLOCATE(MacroSurfaceSpecVal(1:nVarSpec , 1:nSurfSample , 1:nSurfSample , nComputeNodeSurfSides , nSpecies))
MacroSurfaceSpecVal = 0.

IF (CalcSurfCollis%Output) THEN
  ALLOCATE(CounterTotal(1:nSpecies))
  ALLOCATE(SumCounterTotal(1:nSpecies+1))
  CounterTotal(1:nSpecies)      = 0
  SumCounterTotal(1:nSpecies+1) = 0
END IF

OutputCounter = 0

DO iSurfSide = 1,nComputeNodeSurfSides
  OutputCounter = OutputCounter + 1
  !================== INNER BC CHECK
  GlobalSideID = SurfSide2GlobalSide(SURF_SIDEID,iSurfSide)
  IF(SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID).GT.0) THEN
    IF(GlobalSideID.LT.SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID)) THEN
      SurfSideNb = GlobalSide2SurfSide(SURF_SIDEID,SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID))
      SampWallState_Shared(:,:,:,iSurfSide) = SampWallState_Shared(:,:,:,iSurfSide) + SampWallState_Shared(:,:,:,SurfSideNb)
    ELSE
      CYCLE
    END IF
  END IF
  !================== INNER BC CHECK
  DO q = 1,nSurfSample
    DO p = 1,nSurfSample
      CounterSum = SUM(SampWallState_Shared(SAMPWALL_NVARS+1:SAMPWALL_NVARS+nSpecies,p,q,iSurfSide))

      ! even if no impacts happened, sampling is necessary due to surfacemodels -> DO NOT CYCLE
      IF(VarTimeStep%UseVariableTimeStep .AND. CounterSum.GT.0.0) THEN
        TimeSampleTemp = TimeSample * SampWallState_Shared(SAMPWALL_NVARS+nSpecies+1,p,q,iSurfSide) / CounterSum
      ELSE
        TimeSampleTemp = TimeSample
      END IF

      ! Force per area in x,y,z-direction
      MacroSurfaceVal(1:3,p,q,OutputCounter) = SampWallState_Shared(SAMPWALL_DELTA_MOMENTUMX:SAMPWALL_DELTA_MOMENTUMZ,p,q,iSurfSide) &
                                             / (SurfSideArea_Shared(p,q,iSurfSide)*TimeSampleTemp)
      ! Deleting the z-component for 2D/axisymmetric simulations
      IF(Symmetry2D) MacroSurfaceVal(3,p,q,OutputCounter) = 0.
      MacroSurfaceVal(4,p,q,OutputCounter) = (SampWallState_Shared(SAMPWALL_ETRANSOLD,p,q,iSurfSide)  &
                                           +  SampWallState_Shared(SAMPWALL_EROTOLD  ,p,q,iSurfSide)  &
                                           +  SampWallState_Shared(SAMPWALL_EVIBOLD  ,p,q,iSurfSide)  &
                                           -  SampWallState_Shared(SAMPWALL_ETRANSNEW,p,q,iSurfSide)  &
                                           -  SampWallState_Shared(SAMPWALL_EROTNEW  ,p,q,iSurfSide)  &
                                           -  SampWallState_Shared(SAMPWALL_EVIBNEW  ,p,q,iSurfSide)) &
                                           / (SurfSideArea_Shared(p,q,iSurfSide) * TimeSampleTemp)

      nVarCount = 5
      BCID = SideInfo_Shared(SIDE_BCID,SurfSide2GlobalSide(SURF_SIDEID,iSurfSide))
      IF (PartBound%Reactive(PartBound%MapToPartBC(BC(BCID)))) THEN
        ! Currently not working with new halo region
        RETURN

        MacroSurfaceVal(4,p,q,OutputCounter) = MacroSurfaceVal(4,p,q,OutputCounter) + &
                                          (-  SampWall(iSurfSide)%SurfModelState(1,p,q)  &
                                           -  SampWall(iSurfSide)%SurfModelState(2,p,q)  &
                                           -  SampWall(iSurfSide)%SurfModelState(3,p,q)  &
                                           -  SampWall(iSurfSide)%SurfModelState(4,p,q)  &
                                           -  SampWall(iSurfSide)%SurfModelState(5,p,q) )&
                                           / (SurfSideArea_Shared(p,q,iSurfSide) * TimeSampleTemp)

        DO iAdsSampl=1, nAdsSamples
          ! Note: the if-statement is required to prevent the output of "-0" in the .h5 file!
          IF(ABS(SampWall(iSurfSide)%SurfModelState(iAdssampl,p,q)).GT.0.0)THEN
            MacroSurfaceVal(nVarCount+iAdsSampl,p,q,OutputCounter) = (-SampWall(iSurfSide)%SurfModelState(iAdssampl,p,q))&
                                             /(SurfMesh%SurfaceArea(p,q,iSurfSide) * TimeSampleTemp)
          END IF ! ABS(SampWall(iSurfSide)%SurfModelState(iAdssampl,p,q)).GT.0.0
        END DO
        nVarCount = nVarCount + nAdsSamples
      END IF

      IF(nPorousBC.GT.0) THEN
        DO iPBC=1, nPorousBC
          IF(PorousBCInfo_Shared(1,MapSurfSideToPorousSide_Shared(iSurfSide)).EQ.iPBC) THEN
            ! Pump capacity is already in cubic meter per second (diving by the number of iterations)
            MacroSurfaceVal(nVarCount+iPBC,p,q,OutputCounter) = SampWallPumpCapacity_Shared(iSurfSide) * dt / TimeSample
          END IF
        END DO
      END IF

      DO iSpec=1,nSpecies
        ASSOCIATE( nColl    => SampWallState_Shared(SAMPWALL_NVARS+iSpec,p,q,iSurfSide) )
          IF (CalcSurfCollis%Output) CounterTotal(iSpec) = CounterTotal(iSpec) + INT(nColl)

          ! Sum up all Collisions with SpeciesFlags for output
          IF (CalcSurfCollis%SpeciesFlags(iSpec)) THEN
            MacroSurfaceVal(5,p,q,OutputCounter) = MacroSurfaceVal(5,p,q,OutputCounter) + nColl/TimeSample
          END IF

          idx = 1
          MacroSurfaceSpecVal(idx,p,q,OutputCounter,iSpec) = nColl / TimeSample
          IF (PartBound%Reactive(PartBound%MapToPartBC(BC(BCID)))) THEN
            ! calculate accomodation coefficient
            idx = idx + 1
            IF (nColl.EQ.0) THEN
              MacroSurfaceSpecVal(idx,p,q,OutputCounter,iSpec) = 0.
            ELSE
              MacroSurfaceSpecVal(idx,p,q,OutputCounter,iSpec) = SampWall(iSurfSide)%Accomodation(iSpec,p,q) / nColl
            END IF
            ! calculate coverage
            idx = idx + 1
            IF(Adsorption%NumCovSamples.GT.0)THEN
              MacroSurfaceSpecVal(idx,p,q,OutputCounter,iSpec) = SampWall(iSurfSide)%SurfModelState(5+iSpec,p,q) / &
                                                                                                      Adsorption%NumCovSamples
            END IF ! Adsorption%NumCovSamples.GT.0
            ! calculate reaction counters
            DO iReact=1,Adsorption%ReactNum
              ! first part are surface collision processes
              MacroSurfaceSpecVal(idx+iReact,p,q,OutputCounter,iSpec) = &
                  SampWall(iSurfSide)%SurfModelReactCount(iReact,iSpec,p,q) / TimeSample
              ! second part are adsorbate processes
              MacroSurfaceSpecVal(idx+Adsorption%ReactNum+iReact,p,q,OutputCounter,iSpec) = &
                  SampWall(iSurfSide)%SurfModelReactCount(Adsorption%ReactNum+iReact,iSpec,p,q) / TimeSample
            END DO
          END IF

          ! Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z) and angle
          IF(CalcSurfaceImpact)THEN
            ASSOCIATE( nImpacts => SampWallImpactNumber_Shared(iSpec,p,q,iSurfSide))
              IF(nImpacts.GT.0.)THEN
                ! Add average impact energy for each species (trans, rot, vib)
                idx = idx + 1
                MacroSurfaceSpecVal(idx,p,q,OutputCounter,iSpec) = SampWallImpactEnergy_Shared(iSpec,1,p,q,iSurfSide) / nImpacts
                idx = idx + 1
                MacroSurfaceSpecVal(idx,p,q,OutputCounter,iSpec) = SampWallImpactEnergy_Shared(iSpec,2,p,q,iSurfSide) / nImpacts
                idx = idx + 1
                MacroSurfaceSpecVal(idx,p,q,OutputCounter,iSpec) = SampWallImpactEnergy_Shared(iSpec,3,p,q,iSurfSide) / nImpacts

                ! Add average impact vector (x,y,z) for each species
                idx = idx + 1
                MacroSurfaceSpecVal(idx,p,q,OutputCounter,iSpec) = SampWallImpactVector_Shared(iSpec,1,p,q,iSurfSide) / nImpacts
                idx = idx + 1
                MacroSurfaceSpecVal(idx,p,q,OutputCounter,iSpec) = SampWallImpactVector_Shared(iSpec,2,p,q,iSurfSide) / nImpacts
                idx = idx + 1
                MacroSurfaceSpecVal(idx,p,q,OutputCounter,iSpec) = SampWallImpactVector_Shared(iSpec,3,p,q,iSurfSide) / nImpacts

                ! Add average impact angle for each species
                idx = idx + 1
                MacroSurfaceSpecVal(idx,p,q,OutputCounter,iSpec) = SampWallImpactAngle_Shared(iSpec,p,q,iSurfSide) / nImpacts

                ! Add number of impacts
                idx = idx + 1
                MacroSurfaceSpecVal(idx,p,q,OutputCounter,iSpec) = nImpacts
              ELSE
                idx=idx+8
              END IF !
            END ASSOCIATE
          END IF ! CalcSurfaceImpact
        END ASSOCIATE
      END DO ! iSpec=1,nSpecies
    END DO ! q=1,nSurfSample
  END DO ! p=1,nSurfSample
END DO ! iSurfSide=1,SurfMesh%nOutputSides

IF (CalcSurfCollis%Output) THEN
  SumCounterTotal(1:nSpecies) = CounterTotal
  DO iSpec = 1,nSpecies
    ! Sum up all Collisions with SpeciesFlags for output
    IF (CalcSurfCollis%SpeciesFlags(iSpec)) THEN
      SumCounterTotal(nSpecies+1) = SumCounterTotal(nSpecies+1) + SumCounterTotal(iSpec)
    END IF
  END DO

  SWRITE(UNIT_stdOut,'(A)') ' The following species swaps at walls have been sampled:'
  DO iSpec=1,nSpecies
    SWRITE(*,'(A9,I2,A2,E16.9,A6)') ' Species ',iSpec,': '   ,REAL(SumCounterTotal(iSpec))      / TimeSample,' MP/s;'
  END DO
  SWRITE(*,'(A23,E16.9,A6)')        ' All with SpeciesFlag: ',REAL(SumCounterTotal(nSpecies+1)) / TimeSample,' MP/s.'

  DEALLOCATE(CounterTotal)
  DEALLOCATE(SumCounterTotal)
END IF

CALL WriteSurfSampleToHDF5(TRIM(MeshFile),ActualTime)

DEALLOCATE(MacroSurfaceVal,MacroSurfaceSpecVal)

END SUBROUTINE CalcSurfaceValues


REAL FUNCTION CalcTVib(ChaTVib,MeanEVib,nMax)
!===================================================================================================================================
!> Calculation of the vibrational temperature (zero-point search) for the TSHO (Truncated Simple Harmonic Oscillator)
!===================================================================================================================================
! MODULES
USE MOD_Globals       ,ONLY: abort
USE MOD_Globals_Vars  ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars     ,ONLY: DSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL, INTENT(IN)                :: ChaTVib,MeanEVib  ! Charak TVib, mean vibrational Energy of all molecules
INTEGER, INTENT(IN)             :: nMax              ! INT(CharaTDisss/CharaTVib) + 1
REAL(KIND=8)                    :: LowerVal, UpperVal, MiddleVal, MaxPosiVal  ! upper and lower value of zero point search
REAl(KIND=8)                    :: eps_prec=0.1   ! precision of zero point search
REAL(KIND=8)                    :: ZeroVal1, ZeroVal2 ! both fuction values to compare
!===================================================================================================================================

IF (MeanEVib.GT.0) THEN
  !.... Initial limits for a: lower limit = very small value
  !                           upper limit = max. value allowed by system
  !     zero point = CharaTVib / TVib
  LowerVal  = 1.0/(2.0*nMax)                                    ! Tvib is max for nMax => lower limit = 1.0/nMax
  UpperVal  = LOG(HUGE(MiddleVal*nMax))/nMax-1.0/(2.0 * nMax)   ! upper limit = for max possible EXP(nMax*MiddleVal)-value
  MaxPosiVal = LOG(HUGE(MaxPosiVal))  ! maximum value possible in system
  DO WHILE (ABS(LowerVal-UpperVal).GT.eps_prec)                      !  Let's search the zero point by bisection
    MiddleVal = 0.5*(LowerVal+UpperVal)

    IF ((LowerVal.GT.MaxPosiVal).OR.(MiddleVal.GT.MaxPosiVal)) THEN
       CALL Abort(&
__STAMP__&
,'Cannot find zero point in TVib Calculation Function! CharTVib:',RealInfoOpt=ChaTVib)
    END IF

    ! Calc of actual function values
    ZeroVal1 = DSMC%GammaQuant + 1/(EXP(LowerVal)-1) - nMax/(EXP(nMax*LowerVal)-1) - MeanEVib/(ChaTVib*BoltzmannConst)
    ZeroVal2 = DSMC%GammaQuant + 1/(EXP(MiddleVal)-1) - nMax/(EXP(nMax*MiddleVal)-1) - MeanEVib/(ChaTVib*BoltzmannConst)
    ! decision of direction of bisection
    IF (ZeroVal1*ZeroVal2.LT.0) THEN
      UpperVal = MiddleVal
    ELSE
      LowerVal = MiddleVal
    END IF
  END DO
  CalcTVib = ChaTVib/LowerVal ! LowerVal = CharaTVib / TVib
ELSE
  CalcTVib = 0
END IF

RETURN

END FUNCTION CalcTVib

!-----------------------------------------------------------------------------------------------------------------------------------

REAL FUNCTION CalcTelec(MeanEelec, iSpec)
!===================================================================================================================================
!> Calculation of the electronic temperature (zero-point search)
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars  ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars     ,ONLY: SpecDSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)      :: MeanEelec  !< Mean electronic energy
INTEGER, INTENT(IN)   :: iSpec      !< Species index
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER               :: ii
REAL                  :: LowerTemp, UpperTemp, MiddleTemp !< Upper, lower and final value of modified zero point search
REAL,PARAMETER        :: eps_prec=1E-3           !< Relative precision of root-finding algorithm
REAL                  :: TempRatio, SumOne, SumTwo        !< Sums of the electronic partition function
!===================================================================================================================================

IF (MeanEelec.GT.0) THEN
  ! Lower limit: very small value or lowest temperature if ionized
  IF (SpecDSMC(iSpec)%ElectronicState(2,0).EQ.0.0) THEN
    LowerTemp = 1.0
  ELSE
    LowerTemp = SpecDSMC(iSpec)%ElectronicState(2,0)
  END IF
  ! Upper limit: Last excitation level (ionization limit)
  UpperTemp = SpecDSMC(iSpec)%ElectronicState(2,SpecDSMC(iSpec)%MaxElecQuant-1)
  MiddleTemp = LowerTemp
  DO WHILE (.NOT.ALMOSTEQUALRELATIVE(0.5*(LowerTemp + UpperTemp),MiddleTemp,eps_prec))
    MiddleTemp = 0.5*( LowerTemp + UpperTemp)
    SumOne = 0.0
    SumTwo = 0.0
    DO ii = 0, SpecDSMC(iSpec)%MaxElecQuant-1
      TempRatio = SpecDSMC(iSpec)%ElectronicState(2,ii) / MiddleTemp
      IF(CHECKEXP(TempRatio)) THEN
        SumOne = SumOne + SpecDSMC(iSpec)%ElectronicState(1,ii) * EXP(-TempRatio)
        SumTwo = SumTwo + SpecDSMC(iSpec)%ElectronicState(1,ii) * SpecDSMC(iSpec)%ElectronicState(2,ii) * EXP(-TempRatio)
      END IF
    END DO
    IF ( SumTwo / SumOne .GT. MeanEelec / BoltzmannConst ) THEN
      UpperTemp = MiddleTemp
    ELSE
      LowerTemp = MiddleTemp
    END IF
  END DO
  CalcTelec = MiddleTemp
ELSE
  CalcTelec = 0. ! sup
END IF

RETURN

END FUNCTION CalcTelec


REAL FUNCTION CalcTVibPoly(MeanEVib, iSpec)
!===================================================================================================================================
!> Calculation of the vibrational temperature (zero-point search) for polyatomic molecules
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars  ,ONLY: BoltzmannConst, ElementaryCharge
USE MOD_DSMC_Vars     ,ONLY: SpecDSMC, PolyatomMolDSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: MeanEVib  ! Charak TVib, mean vibrational Energy of all molecules
INTEGER, INTENT(IN)             :: iSpec      ! Number of Species
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                 :: iDOF, iPolyatMole
REAL                    :: LowerTemp, UpperTemp, MiddleTemp !< Upper, lower and final value of modified zero point search
REAL                    :: EGuess                           !< Energy value at the current MiddleTemp
REAL,PARAMETER          :: eps_prec=5E-3                    !< Relative precision of root-finding algorithm
!===================================================================================================================================

! lower limit: very small value or lowest temperature if ionized
! upper limit: highest possible temperature
iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
IF (MeanEVib.GT.SpecDSMC(iSpec)%EZeroPoint) THEN
  LowerTemp = 1.0
  UpperTemp = 5.0*SpecDSMC(iSpec)%Ediss_eV*ElementaryCharge/BoltzmannConst
  MiddleTemp = LowerTemp
  DO WHILE (.NOT.ALMOSTEQUALRELATIVE(0.5*(LowerTemp + UpperTemp),MiddleTemp,eps_prec))
    MiddleTemp = 0.5*(LowerTemp + UpperTemp)
    EGuess = SpecDSMC(iSpec)%EZeroPoint
    DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
      ASSOCIATE(CharTVib => PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF))
        IF(CHECKEXP(CharTVib/MiddleTemp)) THEN
          EGuess = EGuess + BoltzmannConst * CharTVib / (EXP(CharTVib/MiddleTemp) - 1.0)
        END IF
      END ASSOCIATE
    END DO
    IF (EGuess.GT.MeanEVib) THEN
      UpperTemp = MiddleTemp
    ELSE
      LowerTemp = MiddleTemp
    END IF
  END DO
  CalcTVibPoly = MiddleTemp
ELSE
  CalcTVibPoly = 0. ! sup
END IF
RETURN

END FUNCTION CalcTVibPoly


REAL FUNCTION CalcMeanFreePath(SpecPartNum, nPart, Volume, opt_omega, opt_temp)
!===================================================================================================================================
!> Calculation of the mean free path for the hard sphere and variable hard sphere (if omega and temperature are given)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars  ,ONLY: Pi
USE MOD_Particle_Vars ,ONLY: Species, nSpecies, usevMPF
USE MOD_DSMC_Vars     ,ONLY: SpecDSMC, RadialWeighting
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: Volume,SpecPartNum(:),nPart
REAL, OPTIONAL, INTENT(IN)      :: opt_omega, opt_temp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                         :: iSpec, jSpec
REAL                            :: DrefMixture, omega, Temp, MFP_Tmp, MacroParticleFactor
!===================================================================================================================================
DrefMixture = 0.0
CalcMeanFreePath = 0.0

IF (nPart.LE.1 .OR. ALL(SpecPartNum.EQ.0.) .OR.Volume.EQ.0) RETURN
! Calculation of mixture reference diameter
IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
  MacroParticleFactor = 1.
ELSE
  MacroParticleFactor = Species(1)%MacroParticleFactor
END IF

! Calculation of mixture reference diameter
DO iSpec = 1, nSpecies
  DrefMixture = DrefMixture + SpecPartNum(iSpec)*SpecDSMC(iSpec)%DrefVHS / nPart
END DO
! Calculation of mean free path for a gas mixture (Bird 1986, p. 96, Eq. 4.77)
! (only defined for a single weighting factor, if omega is present calculation of the mean free path with the VHS model)
IF(PRESENT(opt_omega).AND.PRESENT(opt_temp)) THEN
  omega = opt_omega
  Temp = opt_temp
  IF(Temp.LE.0.0) RETURN
    DO iSpec = 1, nSpecies
      MFP_Tmp = 0.0
      IF(SpecPartNum(iSpec).GT.0.0) THEN ! skipping species not present in the cell
        DO jSpec = 1, nSpecies
          IF(SpecPartNum(jSpec).GT.0.0) THEN ! skipping species not present in the cell
            MFP_Tmp = MFP_Tmp + (Pi*DrefMixture**2.*SpecPartNum(jSpec)*MacroParticleFactor / Volume &
                                  * (SpecDSMC(iSpec)%TrefVHS/Temp)**(omega) &
                                  * SQRT(1+Species(iSpec)%MassIC/Species(jSpec)%MassIC))
          END IF
        END DO
        CalcMeanFreePath = CalcMeanFreePath + (SpecPartNum(iSpec) / nPart) / MFP_Tmp
      END IF
    END DO
ELSE
  DO iSpec = 1, nSpecies
    MFP_Tmp = 0.0
    IF(SpecPartNum(iSpec).GT.0.0) THEN ! skipping species not present in the cell
      DO jSpec = 1, nSpecies
        IF(SpecPartNum(jSpec).GT.0.0) THEN ! skipping species not present in the cell
          MFP_Tmp = MFP_Tmp + (Pi*DrefMixture**2.*SpecPartNum(jSpec)*MacroParticleFactor / Volume &
                                * SQRT(1+Species(iSpec)%MassIC/Species(jSpec)%MassIC))
        END IF
      END DO
      CalcMeanFreePath = CalcMeanFreePath + (SpecPartNum(iSpec) / nPart) / MFP_Tmp
    END IF
  END DO
END IF
RETURN

END FUNCTION CalcMeanFreePath


SUBROUTINE CalcGammaVib()
!===================================================================================================================================
!> calculate Gamma_vib factor necessary for correction of vibrational relaxation according to Gimelshein et al.
!> -> 'Vibrational Relaxation Rates in the DSMC Method', Physics of Fluids V14 No12, 2002
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars ,ONLY: nSpecies
USE MOD_DSMC_Vars     ,ONLY: SpecDSMC, PolyatomMolDSMC, DSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSpec, iDOF, iPolyatMole
REAL                  :: CharaTVib, TempTrans, GammaVib
!===================================================================================================================================

! Calculate GammaVib Factor  = Xi_Vib² * exp(CharaTVib/T_trans) / 2
DO iSpec = 1, nSpecies
  IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    ! First, reset the GammaVib array/value
    IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
      IF (DSMC%PolySingleMode) THEN
        PolyatomMolDSMC(iPolyatMole)%GammaVib = 0.
      ELSE
        SpecDSMC(iSpec)%GammaVib = 0.
      END IF
    ELSE
      SpecDSMC(iSpec)%GammaVib = 0.
    END IF
    TempTrans = DSMC%InstantTransTemp(iSpec)
    IF(TempTrans.GT.0.0) THEN
      IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
        CharaTVib = MAXVAL(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
      ELSE
        CharaTVib = SpecDSMC(iSpec)%CharaTVib
      END IF
      IF(CharaTVib/TempTrans.LT.80) THEN
        ! If CharaTVib/TempTrans is too high the exp function can produce NAN
        ! CharaTVib/TempTrans=80 results in a of GammaVib=2.31020977644213E-31
        IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
          DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
            CharaTVib = PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
            GammaVib = (2.*CharaTVib / (TempTrans *(EXP(CharaTVib/TempTrans)-1.)))**2. * EXP(CharaTVib/TempTrans) / 2.
            IF (DSMC%PolySingleMode) THEN
              PolyatomMolDSMC(iPolyatMole)%GammaVib(iDOF) = GammaVib
            ELSE
              SpecDSMC(iSpec)%GammaVib = SpecDSMC(iSpec)%GammaVib + GammaVib
            END IF
          END DO
        ELSE
          CharaTVib = SpecDSMC(iSpec)%CharaTVib
          SpecDSMC(iSpec)%GammaVib = (2.*CharaTVib / (TempTrans *(EXP(CharaTVib/TempTrans)-1.)))**2. * EXP(CharaTVib/TempTrans) / 2.
        END IF
      END IF  ! CharaTVib/TempTrans.LT.80
    END IF    ! TempTrans.GT.0.0
  END IF      ! (SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)
END DO        ! iSpec = 1, nSpecies

END SUBROUTINE CalcGammaVib


SUBROUTINE CalcInstantElecTempXi(iPartIndx,PartNum)
!===================================================================================================================================
!> Calculation of the instantaneous translational temperature for the cell
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars  ,ONLY: BoltzmannConst
USE MOD_Preproc
USE MOD_DSMC_Vars     ,ONLY: DSMC, CollInf, PartStateIntEn
USE MOD_Particle_Vars ,ONLY: PartSpecies, Species, nSpecies
USE MOD_part_tools    ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)   :: PartNum
INTEGER, INTENT(IN)   :: iPartIndx(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSpec, iPart, SpecPartNum_Simu(nSpecies), PartID, SpecID
REAL                  :: ElecEnergy(nSpecies), partWeight
!===================================================================================================================================
! Actual number of particles, required to avoid calculation of temperature from one particle of the species
SpecPartNum_Simu = 0
! Sum of particle number, might be weighted/multiplied with PartMPF and/or VariableTimeStep
! Setting temperature to zero
DSMC%InstantTXiElec = 0.
ElecEnergy=0.0

DO iPart=1,PartNum
  PartID = iPartIndx(iPart)
  SpecID = PartSpecies(PartID)
  partWeight = GetParticleWeight(PartID)
  ElecEnergy(SpecID) = ElecEnergy(SpecID) + PartStateIntEn(3,PartID) * partWeight
  SpecPartNum_Simu(SpecID) = SpecPartNum_Simu(SpecID) + 1
END DO

DO iSpec=1, nSpecies
  IF(SpecPartNum_Simu(iSpec).GT.1) THEN
    ElecEnergy(iSpec) = ElecEnergy(iSpec) / CollInf%Coll_SpecPartNum(iSpec)
    ! Compute temperatures
    DSMC%InstantTXiElec(1,iSpec) = CalcTelec(ElecEnergy(iSpec), iSpec)
    DSMC%InstantTXiElec(2,iSpec) = 2.*ElecEnergy(iSpec) /(BoltzmannConst*DSMC%InstantTXiElec(1,iSpec))
  ELSE
    DSMC%InstantTXiElec(1:2,iSpec) = 0.0
  END IF
END DO
END SUBROUTINE CalcInstantElecTempXi


SUBROUTINE CalcInstantTransTempOld(iPartIndx,PartNum)
!===================================================================================================================================
!> Calculation of the instantaneous translational temperature for the cell
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars  ,ONLY: BoltzmannConst
USE MOD_Preproc
USE MOD_DSMC_Vars     ,ONLY: DSMC, CollInf
USE MOD_Particle_Vars ,ONLY: PartState, PartSpecies, Species, nSpecies
USE MOD_part_tools    ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)   :: PartNum
INTEGER, INTENT(IN)   :: iPartIndx(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSpec, iPart, SpecPartNum_Simu(nSpecies), PartID, SpecID
REAL                  :: PartV(nSpecies,3), PartV2(nSpecies,3), SumSpecPartNum
REAL                  :: MeanPartV_2(nSpecies,3), Mean_PartV2(nSpecies,3), TempDirec(nSpecies,3)
!===================================================================================================================================

PartV = 0.
PartV2 = 0.
! Actual number of particles, required to avoid calculation of temperature from one particle of the species
SpecPartNum_Simu = 0
! Sum of particle number, might be weighted/multiplied with PartMPF and/or VariableTimeStep
SumSpecPartNum = 0.
! Setting temperature to zero
DSMC%InstantTransTemp = 0.

DO iPart=1,PartNum
  PartID = iPartIndx(iPart)
  SpecID = PartSpecies(PartID)
  PartV(SpecID,1:3) = PartV(SpecID,1:3) + PartState(4:6,PartID) * GetParticleWeight(PartID)
  PartV2(SpecID,1:3) = PartV2(SpecID,1:3) + PartState(4:6,PartID)**2 * GetParticleWeight(PartID)
  SpecPartNum_Simu(SpecID) = SpecPartNum_Simu(SpecID) + 1
END DO

DO iSpec=1, nSpecies
  IF(SpecPartNum_Simu(iSpec).GT.1) THEN
    ! Compute velocity averages
    MeanPartV_2(iSpec,1:3)  = (PartV(iSpec,1:3) / CollInf%Coll_SpecPartNum(iSpec))**2       ! < |v| >**2
    Mean_PartV2(iSpec,1:3)  = PartV2(iSpec,1:3) / CollInf%Coll_SpecPartNum(iSpec)           ! < |v|**2 >
    ! Compute temperatures
    TempDirec(iSpec,1:3) = Species(iSpec)%MassIC * (Mean_PartV2(iSpec,1:3) - MeanPartV_2(iSpec,1:3)) &
                          / BoltzmannConst ! Temp calculation is limitedt to one species
    DSMC%InstantTransTemp(iSpec) = (TempDirec(iSpec,1) + TempDirec(iSpec,2) + TempDirec(iSpec,3)) / 3.
    DSMC%InstantTransTemp(nSpecies + 1) = DSMC%InstantTransTemp(nSpecies + 1)   &
                                          + DSMC%InstantTransTemp(iSpec)*CollInf%Coll_SpecPartNum(iSpec)
    ! Summing up the weights to avoid adding single particles of a species, which do not have a temperature
    SumSpecPartNum = SumSpecPartNum + CollInf%Coll_SpecPartNum(iSpec)
  ELSE
    MeanPartV_2(iSpec,1:3) = 0.
    Mean_PartV2(iSpec,1:3) = 0.
  END IF
END DO

IF(SumSpecPartNum.GT.0) DSMC%InstantTransTemp(nSpecies+1) = DSMC%InstantTransTemp(nSpecies + 1) / SumSpecPartNum

END SUBROUTINE CalcInstantTransTempOld

SUBROUTINE CalcInstantTransTemp(iPartIndx,PartNum)
!===================================================================================================================================
!> Calculation of the instantaneous translational temperature for the cell
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars  ,ONLY: BoltzmannConst
USE MOD_Preproc
USE MOD_DSMC_Vars     ,ONLY: DSMC, CollInf
USE MOD_Particle_Vars ,ONLY: PartState, PartSpecies, Species, nSpecies
USE MOD_part_tools    ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)   :: PartNum
INTEGER, INTENT(IN)   :: iPartIndx(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSpec, iPart, SpecPartNum_Simu(nSpecies), PartID, SpecID
REAL                  :: PartV(nSpecies,3), PartV2(nSpecies), SumSpecPartNum, partWeight, vmag2, V_rel(3), tempweighttotal
REAL                  :: MeanPartV(nSpecies,3), totalweight2(nSpecies),totalweight(nSpecies), Ener(nSpecies), EnerTotal, tempmass
REAL                  :: vBulkTemp(3)
!===================================================================================================================================

PartV = 0.
PartV2 = 0.
! Actual number of particles, required to avoid calculation of temperature from one particle of the species
SpecPartNum_Simu = 0
! Sum of particle number, might be weighted/multiplied with PartMPF and/or VariableTimeStep
SumSpecPartNum = 0.
! Setting temperature to zero
DSMC%InstantTransTemp = 0.
totalweight2 = 0.; totalweight = 0.
DO iPart=1,PartNum
  PartID = iPartIndx(iPart)
  SpecID = PartSpecies(PartID)
  partWeight = GetParticleWeight(PartID)
  PartV(SpecID,1:3) = PartV(SpecID,1:3) + PartState(4:6,PartID) * partWeight
  totalweight(SpecID) = totalweight(SpecID) + partWeight
  totalweight2(SpecID) = totalweight2(SpecID) + partWeight*partWeight
  SpecPartNum_Simu(SpecID) = SpecPartNum_Simu(SpecID) + 1
END DO

DO iSpec=1, nSpecies
  IF(SpecPartNum_Simu(iSpec).GT.1) THEN
    ! Compute velocity averages
    MeanPartV(iSpec,1:3)  = PartV(iSpec,1:3) / totalWeight(iSpec)      
  ELSE
    MeanPartV(iSpec,1:3) = 0.
  END IF
END DO

DO iPart=1,PartNum
  PartID = iPartIndx(iPart)
  SpecID = PartSpecies(PartID)
  partWeight = GetParticleWeight(PartID)
  V_rel(1:3) = PartState(4:6,PartID)-MeanPartV(SpecID,1:3) 
  vmag2 = V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2
  PartV2(SpecID) = PartV2(SpecID) + vmag2 * partWeight
END DO

Ener=0.0; EnerTotal=0.0
tempweighttotal=0.0; tempmass = 0.0; vBulkTemp=0.0
DO iSpec = 1, nSpecies
  IF ((SpecPartNum_Simu(iSpec).GE.2).AND.(.NOT.ALMOSTZERO(PartV2(iSpec)))) THEN
    DSMC%InstantTransTemp(iSpec) = Species(iSpec)%MassIC * PartV2(iSpec) &
        /(3.0*BoltzmannConst*(totalweight(iSpec) - totalweight2(iSpec)/totalweight(iSpec)))
    Ener(iSpec) =  3./2.*BoltzmannConst*DSMC%InstantTransTemp(iSpec) * totalweight(iSpec)
    vmag2 = MeanPartV(iSpec,1)**(2.) + MeanPartV(iSpec,2)**(2.) + MeanPartV(iSpec,3)**(2.)
    Ener(iSpec) = Ener(iSpec) + totalweight(iSpec) * Species(iSpec)%MassIC / 2. * vmag2
    EnerTotal = EnerTotal + Ener(iSpec)
    tempweighttotal = tempweighttotal + totalweight(iSpec)
    tempmass = tempmass +  totalweight(iSpec) * Species(iSpec)%MassIC 
    vBulkTemp(1:3) = vBulkTemp(1:3) + MeanPartV(iSpec,1:3)*totalweight(iSpec) * Species(iSpec)%MassIC 
  END IF
END DO
vBulkTemp(1:3) = vBulkTemp(1:3) / tempmass
vmag2 = vBulkTemp(1)*vBulkTemp(1) + vBulkTemp(2)*vBulkTemp(2) + vBulkTemp(3)*vBulkTemp(3)
EnerTotal = EnerTotal -  tempmass / 2. * vmag2
DSMC%InstantTransTemp(nSpecies+1) = 2. * EnerTotal / (3.*tempweighttotal*BoltzmannConst)

END SUBROUTINE CalcInstantTransTemp

SUBROUTINE DSMC_data_sampling()
!===================================================================================================================================
!> Sampling of variables velocity and energy for DSMC
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars              ,ONLY: PartStateIntEn, DSMC, CollisMode, SpecDSMC, DSMC_Solution
USE MOD_DSMC_Vars              ,ONLY: useDSMC, DSMC_VolumeSample, ConsiderVolumePortions
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_Part_tools             ,ONLY: GetParticleWeight
USE MOD_Particle_Vars          ,ONLY: PartState, PDM, PartSpecies, PEM
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers     ,ONLY: LBStartTime, LBPauseTime
#endif /*USE_LOADBALANCE*/
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemVolume_Shared
USE MOD_Mesh_Vars              ,ONLY: offSetElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart, iElem, iSpec
REAL                          :: partWeight
#if USE_LOADBALANCE
REAL                          :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
DSMC%SampNum = DSMC%SampNum + 1
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
DO iPart=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    iSpec = PartSpecies(iPart)
    iElem = PEM%LocalElemID(iPart)
    partWeight = GetParticleWeight(iPart)
    DSMC_Solution(1:3,iElem,iSpec) = DSMC_Solution(1:3,iElem,iSpec) + PartState(4:6,iPart)*partWeight
    DSMC_Solution(4:6,iElem,iSpec) = DSMC_Solution(4:6,iElem,iSpec) + PartState(4:6,iPart)**2*partWeight
    DSMC_Solution(7,iElem,iSpec) = DSMC_Solution(7,iElem, iSpec) + partWeight  !density number
    IF(useDSMC)THEN
      IF ((CollisMode.EQ.2).OR.(CollisMode.EQ.3)) THEN
        IF ((SpecDSMC(PartSpecies(iPart))%InterID.EQ.2).OR.(SpecDSMC(PartSpecies(iPart))%InterID.EQ.20)) THEN
          DSMC_Solution(8,iElem, iSpec) = DSMC_Solution(8,iElem, iSpec) &
            + (PartStateIntEn(1,iPart) - SpecDSMC(iSpec)%EZeroPoint)*partWeight
          DSMC_Solution(9,iElem, iSpec) = DSMC_Solution(9,iElem, iSpec)+PartStateIntEn(2,iPart)*partWeight
        END IF
        IF (DSMC%ElectronicModel) THEN
          IF ((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
            DSMC_Solution(10,iElem,iSpec)=DSMC_Solution(10,iElem,iSpec)+PartStateIntEn(3,iPart)*partWeight
          END IF
        END IF
      END IF
    END IF
    DSMC_Solution(11,iElem, iSpec) = DSMC_Solution(11,iElem, iSpec) + 1.0 !simpartnum
  END IF
END DO
!IF(ConsiderVolumePortions) THEN
!  ! DO iElem=1,nElems
!  !   DSMC_VolumeSample(iElem) = DSMC_VolumeSample(iElem) + ElemVolume_Shared(GetCNElemID(iElem+offSetElem))*(1.-GEO%MPVolumePortion(iElem))
!  ! END DO
!  CALL abort(&
!__STAMP__&
!  ,' OUTPUT OF MACROBUDDIES NOT IMPLEMENTED YET!')
!ELSE
!  DSMC_VolumeSample(1:nElems) = ElemVolume_Shared(1+offSetElem:nElems+offSetElem)
!END IF
#if USE_LOADBALANCE
CALL LBPauseTime(LB_DSMC,tLBStart)
#endif /*USE_LOADBALANCE*/
END SUBROUTINE DSMC_data_sampling


SUBROUTINE DSMC_output_calc(nVar,nVar_quality,nVarloc,DSMC_MacroVal)
!===================================================================================================================================
!> Subroutine to calculate the solution U for writing into HDF5 format DSMC_output
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_PreProc
USE MOD_BGK_Vars              ,ONLY: BGKInitDone, BGK_QualityFacSamp
USE MOD_DSMC_Vars             ,ONLY: DSMC_Solution, DSMC_VolumeSample, CollisMode, SpecDSMC, DSMC, useDSMC, RadialWeighting
USE MOD_FPFlow_Vars           ,ONLY: FPInitDone, FP_QualityFacSamp
USE MOD_Mesh_Vars             ,ONLY: nElems
USE MOD_Particle_Vars         ,ONLY: Species, nSpecies, WriteMacroVolumeValues, usevMPF, VarTimeStep, Symmetry2D
USE MOD_Particle_VarTimeStep  ,ONLY: CalcVarTimeStep
USE MOD_Restart_Vars          ,ONLY: RestartTime
USE MOD_TimeDisc_Vars         ,ONLY: time,TEnd,iter,dt
USE MOD_Particle_Mesh_Vars    ,ONLY: ElemMidPoint_Shared, ElemVolume_Shared
USE MOD_Mesh_Vars             ,ONLY: offSetElem
USE MOD_Mesh_Tools            ,ONLY: GetCNElemID
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)      :: nVar,nVar_quality,nVarloc
REAL,INTENT(INOUT)      :: DSMC_MacroVal(1:nVar+nVar_quality,nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iElem, iSpec, nVarCount, nSpecTemp, nVarCountRelax
REAL                    :: TVib_TempFac, iter_loc
REAL                    :: MolecPartNum, HeavyPartNum
!===================================================================================================================================
! nullify
DSMC_MacroVal = 0.0

nVarCount=0
DO iElem = 1, nElems ! element/cell main loop
  MolecPartNum = 0.0
  HeavyPartNum = 0.0
  ! Avoid the output and calculation of total values for a single species, associate construct for Total_ points to the same
  ! array as the Macro_ link for a single species and total values are only calculated for nSpecies > 1. Workaround required
  ! since associate construct is around the DO iSpec=1,nSpecies loop.
  IF(nSpecies.EQ.1) THEN
    nSpecTemp = 0
  ELSE
    nSpecTemp = nSpecies
  END IF
  ASSOCIATE ( Total_Velo     => DSMC_MacroVal(nVarLoc*nSpecTemp+1:nVarLoc*nSpecTemp+3,iElem) ,&
              Total_Temp     => DSMC_MacroVal(nVarLoc*nSpecTemp+4:nVarLoc*nSpecTemp+6,iElem) ,&
              Total_TempMean => DSMC_MacroVal(nVarLoc*nSpecTemp+12,iElem)            ,&
              Total_Density  => DSMC_MacroVal(nVarLoc*nSpecTemp+7,iElem)             ,&
              Total_TempVib  => DSMC_MacroVal(nVarLoc*nSpecTemp+8,iElem)             ,&
              Total_TempRot  => DSMC_MacroVal(nVarLoc*nSpecTemp+9,iElem)             ,&
              Total_Tempelec => DSMC_MacroVal(nVarLoc*nSpecTemp+10,iElem)            ,&
              Total_PartNum  => DSMC_MacroVal(nVarLoc*nSpecTemp+11,iElem)            ,&
              SimVolume      => ElemVolume_Shared(GetCNElemID(iElem+offSetElem)) &
              )
    ! compute simulation cell volume
    DO iSpec = 1, nSpecies
      ASSOCIATE ( PartVelo   => DSMC_Solution(1:3,iElem,iSpec) ,&
                  PartVelo2  => DSMC_Solution(4:6,iElem,iSpec) ,&
                  PartNum    => DSMC_Solution(7,iElem,iSpec)   ,&
                  PartEvib   => DSMC_Solution(8,iElem,iSpec)   ,&
                  PartErot   => DSMC_Solution(9,iElem,iSpec)   ,&
                  PartEelec  => DSMC_Solution(10,iElem,iSpec)  ,&
                  SimPartNum => DSMC_Solution(11,iElem,iSpec)  ,&
                  Macro_Velo     => DSMC_MacroVal(nVarLoc*(iSpec-1)+1:nVarLoc*(iSpec-1)+3,iElem) ,&
                  Macro_Temp     => DSMC_MacroVal(nVarLoc*(iSpec-1)+4:nVarLoc*(iSpec-1)+6,iElem) ,&
                  Macro_TempMean => DSMC_MacroVal(nVarLoc*(iSpec-1)+12,iElem)                    ,&
                  Macro_Density  => DSMC_MacroVal(nVarLoc*(iSpec-1)+7,iElem)                     ,&
                  Macro_TempVib  => DSMC_MacroVal(nVarLoc*(iSpec-1)+8,iElem)                     ,&
                  Macro_TempRot  => DSMC_MacroVal(nVarLoc*(iSpec-1)+9,iElem)                     ,&
                  Macro_Tempelec => DSMC_MacroVal(nVarLoc*(iSpec-1)+10,iElem)                    ,&
                  Macro_PartNum  => DSMC_MacroVal(nVarLoc*(iSpec-1)+11,iElem) &
                  )
        IF (PartNum.GT.0.0) THEN
          ! simulation particle number
          Macro_PartNum = PartNum / REAL(DSMC%SampNum)
          ! compute flow velocity
          Macro_Velo = PartVelo / PartNum
          ! compute flow Temperature
          Macro_Temp = Species(iSpec)%MassIC/BoltzmannConst * ( (PartVelo2/PartNum)- (PartVelo/PartNum)**2 )
          ! mean flow Temperature
          Macro_TempMean = (Macro_Temp(1) + Macro_Temp(2) + Macro_Temp(3)) / 3.
          ! compute number density
          IF (SimVolume.GT.0) THEN
            IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
              ! PartNum contains the weighted particle number
              Macro_Density = Macro_PartNum / SimVolume
            ELSE
              Macro_Density = Macro_PartNum*Species(iSpec)%MacroParticleFactor /SimVolume
            END IF
          ELSE
            Macro_Density = 0.
          END IF
          ! Compute total values for a gas mixture (nSpecies > 1)
          IF(nSpecies.GT.1) THEN
            Total_PartNum   = Total_PartNum + Macro_PartNum
            Total_Velo      = Total_Velo + Macro_Velo*Macro_PartNum
            Total_Temp      = Total_Temp + Macro_Temp*Macro_PartNum
            Total_TempMean  = Total_TempMean + Macro_TempMean*Macro_PartNum
            Total_Density   = Total_Density + Macro_Density
          END IF
          ! compute internal energies / has to be changed for vfd
          IF(useDSMC)THEN
            IF ((CollisMode.EQ.2).OR.(CollisMode.EQ.3))THEN
              IF ((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
                IF (DSMC%VibEnergyModel.EQ.0) THEN              ! SHO-model
                  IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
                    IF( (PartEvib/PartNum) .GT. 0.0 ) THEN
                      Macro_TempVib = CalcTVibPoly(PartEvib/PartNum + SpecDSMC(iSpec)%EZeroPoint, iSpec)
                    ELSE
                      Macro_TempVib = 0.0
                    END IF
                  ELSE
                    TVib_TempFac = PartEvib / (PartNum * BoltzmannConst * SpecDSMC(iSpec)%CharaTVib)
                    IF ((PartEvib /PartNum).LE.0.0) THEN
                      Macro_TempVib = 0.0
                    ELSE
                      Macro_TempVib = SpecDSMC(iSpec)%CharaTVib / LOG(1. + 1./(TVib_TempFac))
                    END IF
                  END IF
                ELSE                                            ! TSHO-model
                  Macro_TempVib = CalcTVib(SpecDSMC(iSpec)%CharaTVib, &
                    PartEvib/PartNum + SpecDSMC(iSpec)%EZeroPoint, SpecDSMC(iSpec)%MaxVibQuant)
                END IF
                Macro_TempRot = 2. * PartERot / (PartNum*BoltzmannConst*REAL(SpecDSMC(iSpec)%Xi_Rot))
                MolecPartNum = MolecPartNum + Macro_PartNum
                IF(nSpecies.GT.1) THEN
                  Total_TempVib  = Total_TempVib  + Macro_TempVib*Macro_PartNum
                  Total_TempRot  = Total_TempRot  + Macro_TempRot*Macro_PartNum
                END IF
              END IF
              IF (DSMC%ElectronicModel) THEN
                IF ((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
                  Macro_TempElec = CalcTelec(PartEelec/PartNum, iSpec)
                  HeavyPartNum = HeavyPartNum + Macro_PartNum
                END IF
                IF(nSpecies.GT.1) Total_TempElec = Total_TempElec + Macro_TempElec*Macro_PartNum
              END IF
            END IF
          END IF
        END IF
      END ASSOCIATE
    END DO
    ! Compute total values for a gas mixture (nSpecies > 1)
    IF(nSpecies.GT.1) THEN
      IF (Total_PartNum.GT.0.0) THEN
        Total_Velo = Total_Velo / Total_PartNum
        Total_Temp = Total_Temp / Total_PartNum
        Total_TempMean = Total_TempMean / Total_PartNum
        IF(useDSMC)THEN
          IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.(MolecpartNum.GT.0))THEN
            Total_TempVib = Total_TempVib / MolecPartNum
            Total_TempRot = Total_TempRot / MolecPartNum
          END IF
          IF ( DSMC%ElectronicModel .AND.(HeavyPartNum.GT. 0)) THEN
            Total_TempElec = Total_TempElec / HeavyPartNum
          END IF
        END IF
      END IF
    END IF
    ! Radial weighting, vMPF, variable timestep: Getting the actual number of simulation particles without weighting factors
    IF (usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
      Total_PartNum = 0.0
      DO iSpec = 1, nSpecies
        IF(DSMC%SampNum.GT.0) DSMC_MacroVal(nVarLoc*(iSpec-1)+11,iElem) = DSMC_Solution(11,iElem, iSpec) / REAL(DSMC%SampNum)
        IF(nSpecies.GT.1) Total_PartNum = Total_PartNum + DSMC_MacroVal(nVarLoc*(iSpec-1)+11,iElem)
      END DO
    END IF
  END ASSOCIATE
END DO
! write dsmc quality values
IF (DSMC%CalcQualityFactors) THEN
  IF(WriteMacroVolumeValues) THEN
    iter_loc = REAL(DSMC%SampNum)
  ELSE
    IF (RestartTime.GT.(1-DSMC%TimeFracSamp)*TEnd) THEN
      iter_loc = REAL(iter)
    ELSE
      iter_loc = (Time-(1-DSMC%TimeFracSamp)*TEnd) / dt
    END IF
  END IF
  DO iElem=1,nElems
    nVarCount = nVar
    IF(DSMC%QualityFacSamp(iElem,4).GT.0.0) THEN
      DSMC_MacroVal(nVarCount+1:nVarCount+3,iElem) = DSMC%QualityFacSamp(iElem,1:3) / DSMC%QualityFacSamp(iElem,4)
    END IF
    nVarCount = nVar + 3
    IF(VarTimeStep%UseVariableTimeStep) THEN
      IF(VarTimeStep%UseLinearScaling.AND.Symmetry2D) THEN
        ! 2D/Axisymmetric uses a scaling of the time step per particle, no element values are used. For the output simply the cell
        ! midpoint is used to calculate the time step
        VarTimeStep%ElemFac(iElem) = CalcVarTimeStep(ElemMidPoint_Shared(1,GetCNElemID(iElem + offsetElem)), &
                                                     ElemMidPoint_Shared(2,GetCNElemID(iElem + offsetElem)))
      END IF
      DSMC_MacroVal(nVarCount+1,iElem) = VarTimeStep%ElemFac(iElem)
      nVarCount = nVarCount + 1
    END IF
    IF(RadialWeighting%PerformCloning) THEN
      IF(DSMC%QualityFacSamp(iElem,4).GT.0.0) THEN
        DSMC_MacroVal(nVarCount+1:nVarCount+2,iElem)=DSMC%QualityFacSamp(iElem,5:6) / DSMC%QualityFacSamp(iElem,4)
      END IF
      nVarCount = nVarCount + 2
    END IF
    IF(FPInitDone) THEN
      IF(FP_QualityFacSamp(2,iElem).GT.0) THEN
        ! Mean relaxation factor (mean over all octree subcells)
        DSMC_MacroVal(nVarCount+1,iElem) = FP_QualityFacSamp(1,iElem) / FP_QualityFacSamp(2,iElem)
        ! Mean Prandtl number
        DSMC_MacroVal(nVarCount+2,iElem) = FP_QualityFacSamp(6,iElem) / FP_QualityFacSamp(2,iElem)
      END IF
      IF(FP_QualityFacSamp(4,iElem).GT.0) THEN
        ! Max relaxation factor (maximal value of all octree subcells)
        DSMC_MacroVal(nVarCount+3,iElem) = FP_QualityFacSamp(3,iElem) / FP_QualityFacSamp(4,iElem)
        ! Max rotational relaxation factor
        DSMC_MacroVal(nVarCount+4,iElem) = FP_QualityFacSamp(5,iElem) / FP_QualityFacSamp(4,iElem)
      END IF
      ! Ratio between FP and DSMC usage per cell
      DSMC_MacroVal(nVarCount+5,iElem) = FP_QualityFacSamp(4,iElem) / iter_loc
      nVarCount = nVarCount + 5
    END IF
    IF(BGKInitDone) THEN
      IF(BGK_QualityFacSamp(2,iElem).GT.0) THEN
        ! Mean relaxation factor (mean over all octree subcells)
        DSMC_MacroVal(nVarCount+1,iElem) = BGK_QualityFacSamp(1,iElem) / BGK_QualityFacSamp(2,iElem)
      END IF
      IF(BGK_QualityFacSamp(4,iElem).GT.0) THEN
        ! Max relaxation factor (maximal value of all octree subcells)
        DSMC_MacroVal(nVarCount+2,iElem) = BGK_QualityFacSamp(3,iElem) / BGK_QualityFacSamp(4,iElem)
        ! Max rotational relaxation factor
        DSMC_MacroVal(nVarCount+3,iElem) = BGK_QualityFacSamp(5,iElem) / BGK_QualityFacSamp(4,iElem)
      END IF
      ! Ratio between BGK and DSMC usage per cell
      DSMC_MacroVal(nVarCount+4,iElem) = BGK_QualityFacSamp(4,iElem) / iter_loc
      nVarCount = nVarCount + 4
    END IF
    ! variable rotation and vibration relaxation
    IF(Collismode.GT.1) THEN
      IF((DSMC%RotRelaxProb.GE.2).OR.(DSMC%VibRelaxProb.EQ.2)) THEN
        IF(nSpecies.EQ.1) THEN
          nSpecTemp = 0
        ELSE
          nSpecTemp = nSpecies
        END IF
        DO iSpec=0,nSpecTemp
          nVarCountRelax = 13
          IF(DSMC%RotRelaxProb.GE.2) THEN
            IF(DSMC%QualityFacSampRotSamp(iElem,iSpec+1).GT.0) THEN
              DSMC_MacroVal(nVarLoc*(iSpec)+nVarCountRelax,iElem) = DSMC%QualityFacSampRot(iElem,iSpec+1,2) &
                                                                / REAL(DSMC%QualityFacSampRotSamp(iElem,iSpec+1))
              DSMC_MacroVal(nVarLoc*(iSpec)+nVarCountRelax+1,iElem) = DSMC%QualityFacSampRot(iElem,iSpec+1,1) &
                                                                / REAL(DSMC%QualityFacSampRotSamp(iElem,iSpec+1))
              nVarCountRelax = nVarCountRelax + 2
            END IF
          END IF
          IF((DSMC%VibRelaxProb.EQ.2)) THEN
            IF(DSMC%QualityFacSampVibSamp(iElem,iSpec+1,2).GT.0) DSMC_MacroVal(nVarLoc*(iSpec)+nVarCountRelax,iElem) = &
                                        DSMC%QualityFacSampVib(iElem,iSpec+1,2) / REAL(DSMC%QualityFacSampVibSamp(iElem,iSpec+1,2))
            IF(DSMC%QualityFacSampVibSamp(iElem,iSpec+1,1).GT.0) DSMC_MacroVal(nVarLoc*(iSpec)+nVarCountRelax+1,iElem)=&
                                        DSMC%QualityFacSampVib(iElem,iSpec+1,1) / REAL(DSMC%QualityFacSampVibSamp(iElem,iSpec+1,1))
          END IF
        END DO
      END IF
    END IF
  END DO
  IF (ALLOCATED(DSMC%QualityFacSampRot)) DSMC%QualityFacSampRot = 0.
  IF (ALLOCATED(DSMC%QualityFacSampRotSamp)) DSMC%QualityFacSampRotSamp = 0
  IF (ALLOCATED(DSMC%QualityFacSampVib)) DSMC%QualityFacSampVib = 0.
  IF (ALLOCATED(DSMC%QualityFacSampVibSamp)) DSMC%QualityFacSampVibSamp = 0
END IF

END SUBROUTINE DSMC_output_calc


SUBROUTINE WriteDSMCToHDF5(MeshFileName,OutputTime,FutureTime)
!===================================================================================================================================
!> Subroutine to write the solution U to HDF5 format
!> Is used for postprocessing and for restart
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars     ,ONLY: DSMC, RadialWeighting, CollisMode
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars  ,ONLY: ProjectName
USE MOD_Mesh_Vars     ,ONLY: offsetElem,nGlobalElems, nElems
USE MOD_io_HDF5
USE MOD_HDF5_output   ,ONLY: WriteArrayToHDF5
USE MOD_Particle_Vars ,ONLY: nSpecies, VarTimeStep
USE MOD_BGK_Vars      ,ONLY: BGKInitDone
USE MOD_FPFlow_Vars   ,ONLY: FPInitDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName
REAL,INTENT(IN)                :: OutputTime
REAL,INTENT(IN),OPTIONAL       :: FutureTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName
CHARACTER(LEN=255)             :: SpecID
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
INTEGER                        :: nVar,nVar_quality,nVarloc,nVarCount,ALLOCSTAT, iSpec, nVarRelax
REAL,ALLOCATABLE               :: DSMC_MacroVal(:,:)
REAL                           :: StartT,EndT
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE DSMC TO HDF5 FILE...'
#if USE_MPI
  StartT=MPI_WTIME()
#else
  StartT=LOCALTIME()
#endif

! Create dataset attribute "VarNames"
nVarloc=DSMC_NVARS
nVarRelax=0
IF(DSMC%CalcQualityFactors.AND.(CollisMode.GE.2)) THEN
  IF(DSMC%RotRelaxProb.GE.2) nVarRelax = nVarRelax + 2
  IF(DSMC%VibRelaxProb.EQ.2) nVarRelax = nVarRelax + 2
END IF
IF(nSpecies.EQ.1) THEN
  nVar=nVarloc+nVarRelax
ELSE
  nVar=(nVarloc+nVarRelax)*(nSpecies+1)
END IF

IF (DSMC%CalcQualityFactors) THEN
  nVar_quality=3
  IF(VarTimeStep%UseVariableTimeStep) nVar_quality = nVar_quality + 1
  IF(RadialWeighting%PerformCloning) nVar_quality = nVar_quality + 2
  IF(BGKInitDone) nVar_quality = nVar_quality + 4
  IF(FPInitDone) nVar_quality = nVar_quality + 5
ELSE
  nVar_quality=0
END IF
ALLOCATE(StrVarNames(1:nVar+nVar_quality))
nVarCount=0
IF(nSpecies.GT.1) THEN
  DO iSpec=1,nSpecies
    WRITE(SpecID,'(I3.3)') iSpec
    StrVarNames(nVarCount+DSMC_VELOX      )='Spec'//TRIM(SpecID)//'_VeloX'
    StrVarNames(nVarCount+DSMC_VELOY      )='Spec'//TRIM(SpecID)//'_VeloY'
    StrVarNames(nVarCount+DSMC_VELOZ      )='Spec'//TRIM(SpecID)//'_VeloZ'
    StrVarNames(nVarCount+DSMC_TEMPX      )='Spec'//TRIM(SpecID)//'_TempTransX'
    StrVarNames(nVarCount+DSMC_TEMPY      )='Spec'//TRIM(SpecID)//'_TempTransY'
    StrVarNames(nVarCount+DSMC_TEMPZ      )='Spec'//TRIM(SpecID)//'_TempTransZ'
    StrVarNames(nVarCount+DSMC_NUMDENS    )='Spec'//TRIM(SpecID)//'_NumberDensity'
    StrVarNames(nVarCount+DSMC_TVIB       )='Spec'//TRIM(SpecID)//'_TempVib'
    StrVarNames(nVarCount+DSMC_TROT       )='Spec'//TRIM(SpecID)//'_TempRot'
    StrVarNames(nVarCount+DSMC_TELEC      )='Spec'//TRIM(SpecID)//'_TempElec'
    StrVarNames(nVarCount+DSMC_SIMPARTNUM )='Spec'//TRIM(SpecID)//'_SimPartNum'
    StrVarNames(nVarCount+DSMC_TEMPMEAN   )='Spec'//TRIM(SpecID)//'_TempTransMean'
    nVarCount=nVarCount+nVarloc
    IF(DSMC%CalcQualityFactors.AND.(CollisMode.GE.2)) THEN
      IF(DSMC%RotRelaxProb.GE.2) THEN
        StrVarNames(nVarCount+1              )='Spec'//TRIM(SpecID)//'_DSMC_MaxRotRelaxProb'
        StrVarNames(nVarCount+2              )='Spec'//TRIM(SpecID)//'_DSMC_MeanRotRelaxProb'
        nvarcount=nvarcount+2
      END IF
      IF((DSMC%VibRelaxProb.EQ.2)) THEN
        StrVarNames(nVarCount+1              )='Spec'//TRIM(SpecID)//'_DSMC_MaxVibRelaxProb'
        StrVarNames(nVarCount+2              )='Spec'//TRIM(SpecID)//'_DSMC_MeanVibRelaxProb'
        nvarcount=nvarcount+2
      END IF
    END IF
  END DO ! iSpec=1,nSpecies
END IF
! fill varnames for total values
StrVarNames(nVarCount+DSMC_VELOX      )='Total_VeloX'
StrVarNames(nVarCount+DSMC_VELOY      )='Total_VeloY'
StrVarNames(nVarCount+DSMC_VELOZ      )='Total_VeloZ'
StrVarNames(nVarCount+DSMC_TEMPX      )='Total_TempTransX'
StrVarNames(nVarCount+DSMC_TEMPY      )='Total_TempTransY'
StrVarNames(nVarCount+DSMC_TEMPZ      )='Total_TempTransZ'
StrVarNames(nVarCount+DSMC_NUMDENS    )='Total_NumberDensity'
StrVarNames(nVarCount+DSMC_TVIB       )='Total_TempVib'
StrVarNames(nVarCount+DSMC_TROT       )='Total_TempRot'
StrVarNames(nVarCount+DSMC_TELEC      )='Total_TempElec'
StrVarNames(nVarCount+DSMC_SIMPARTNUM )='Total_SimPartNum'
StrVarNames(nVarCount+DSMC_TEMPMEAN   )='Total_TempTransMean'
nVarCount=nVarCount+nVarloc
IF(DSMC%CalcQualityFactors.AND.(CollisMode.GE.2)) THEN
  IF(DSMC%RotRelaxProb.GE.2) THEN
    StrVarNames(nVarCount+1              )='Total_DSMC_MaxRotRelaxProb'
    StrVarNames(nVarCount+2              )='Total_DSMC_MeanRotRelaxProb'
    nvarcount=nvarcount+2
  END IF
  IF((DSMC%VibRelaxProb.EQ.2)) THEN
    StrVarNames(nVarCount+1              )='Total_DSMC_MaxVibRelaxProb'
    StrVarNames(nVarCount+2              )='Total_DSMC_MeanVibRelaxProb'
    nvarcount=nvarcount+2
  END IF
END IF

IF (DSMC%CalcQualityFactors) THEN
  StrVarNames(nVarCount+1) ='DSMC_MaxCollProb'
  StrVarNames(nVarCount+2) ='DSMC_MeanCollProb'
  StrVarNames(nVarCount+3) ='DSMC_MCS_over_MFP'
  nVarCount=nVarCount+3
  IF(VarTimeStep%UseVariableTimeStep) THEN
    StrVarNames(nVarCount+1) ='VariableTimeStep'
    nVarCount = nVarCount + 1
  END IF
  IF(RadialWeighting%PerformCloning) THEN
    StrVarNames(nVarCount+1) = '2D_ClonesInCell'
    StrVarNames(nVarCount+2) = '2D_IdenticalParticles'
    nVarCount=nVarCount+2
  END IF
  IF(BGKInitDone) THEN
    StrVarNames(nVarCount+1) ='BGK_MeanRelaxationFactor'
    StrVarNames(nVarCount+2) ='BGK_MaxRelaxationFactor'
    StrVarNames(nVarCount+3) ='BGK_MaxRotationRelaxFactor'
    StrVarNames(nVarCount+4) ='BGK_DSMC_Ratio'
    nVarCount=nVarCount+4
  END IF
  IF(FPInitDone) THEN
    StrVarNames(nVarCount+1) ='FP_MeanRelaxationFactor'
    StrVarNames(nVarCount+2) ='FP_MeanPrandtlNumber'
    StrVarNames(nVarCount+3) ='FP_MaxRelaxationFactor'
    StrVarNames(nVarCount+4) ='FP_MaxRotationRelaxFactor'
    StrVarNames(nVarCount+5) ='FP_DSMC_Ratio'
    nVarCount=nVarCount+5
  END IF
END IF

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_DSMCState',OutputTime))//'.h5'
IF(MPIRoot) CALL GenerateDSMCFileSkeleton('DSMCState',nVar+nVar_quality,StrVarNames,MeshFileName,OutputTime,FutureTime)

#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif

CALL OpenDataFile(FileName,create=.false.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_WORLD)

ALLOCATE(DSMC_MacroVal(1:nVar+nVar_quality,nElems), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
  ,' Cannot allocate output array DSMC_MacroVal array!')
END IF
CALL DSMC_output_calc(nVar,nVar_quality,nVarloc+nVarRelax,DSMC_MacroVal)

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
  nVarX        => INT(nVar+nVar_quality,IK)    ,&
  PP_nElems    => INT(PP_nElems,IK)            ,&
  offsetElem   => INT(offsetElem,IK)           ,&
  nGlobalElems => INT(nGlobalElems,IK)         )
  CALL WriteArrayToHDF5(DataSetName='ElemData' , rank=2         , &
                        nValGlobal =(/nVarX    , nGlobalElems/) , &
                        nVal       =(/nVarX    , PP_nElems/)    , &
                        offset     =(/0_IK     , offsetElem/)   , &
                        collective =.false.,  RealArray=DSMC_MacroVal(:,:))
END ASSOCIATE

CALL CloseDataFile()

DEALLOCATE(StrVarNames)
DEALLOCATE(DSMC_MacroVal)
#if USE_MPI
IF(MPIROOT) EndT=MPI_WTIME()
#else
EndT=LOCALTIME()
#endif

SWRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',EndT-StartT,'s]'

END SUBROUTINE WriteDSMCToHDF5


SUBROUTINE GenerateDSMCFileSkeleton(TypeString,nVar,StrVarNames,MeshFileName,OutputTime,FutureTime)
!===================================================================================================================================
!> Subroutine that generates the output file on a single processor and writes all the necessary attributes (better MPI performance)
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Globals_Vars  ,ONLY: ProjectName
USE MOD_io_HDF5
USE MOD_HDF5_Output   ,ONLY: WriteAttributeToHDF5, WriteHDF5Header
USE MOD_Particle_Vars ,ONLY: nSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: TypeString
INTEGER,INTENT(IN)             :: nVar
CHARACTER(LEN=255)             :: StrVarNames(nVar)
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName
REAL,INTENT(IN)                :: OutputTime
REAL,INTENT(IN),OPTIONAL       :: FutureTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName,MeshFile255
!===================================================================================================================================
! Create file
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(TypeString),OutputTime))//'.h5'
CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)

CALL WriteHDF5Header(TRIM('DSMCState'),File_ID)

! Write dataset properties "Time","MeshFile","NextFile","NodeType","VarNames"
CALL WriteAttributeToHDF5(File_ID,'Time',1,RealScalar=OutputTime)
CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=(/TRIM(MeshFileName)/))
IF(PRESENT(FutureTime))THEN
  MeshFile255=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(TypeString),FutureTime))//'.h5'
  CALL WriteAttributeToHDF5(File_ID,'NextFile',1,StrScalar=(/MeshFile255/))
END IF
CALL WriteAttributeToHDF5(File_ID,'VarNamesAdd',nVar,StrArray=StrVarNames)

CALL WriteAttributeToHDF5(File_ID,'NSpecies',1,IntegerScalar=nSpecies)

CALL CloseDataFile()

END SUBROUTINE GenerateDSMCFileSkeleton


SUBROUTINE WriteAnalyzeSurfCollisToHDF5(OutputTime,TimeSample)
!===================================================================================================================================
!> Wrinting AnalyzeSurfCollis-Data to hdf5 file (based on WriteParticleToHDF5 and WriteDSMCToHDF5)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_Globals_Vars           ,ONLY: ProjectName
USE MOD_io_HDF5
USE MOD_HDF5_Output            ,ONLY: WriteAttributeToHDF5, WriteHDF5Header, WriteArrayToHDF5
USE MOD_PICDepo_Vars           ,ONLY: SFResampleAnalyzeSurfCollis, LastAnalyzeSurfCollis, r_SF
USE MOD_Particle_Boundary_Vars ,ONLY: nPartBound, AnalyzeSurfCollis
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: OutputTime, TimeSample
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: Filename, TypeString, H5_Name
INTEGER,ALLOCATABLE            :: SpeciesPositions(:,:)
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)!,params(:)
#if USE_MPI
INTEGER,ALLOCATABLE            :: sendbuf(:),recvbuf(:)
REAL,ALLOCATABLE               :: sendbuf2(:),recvbuf2(:)
INTEGER                        :: iProc
INTEGER                        :: globalNum(0:nProcessors-1), Displace(0:nProcessors-1), RecCount(0:nProcessors-1)
#endif
INTEGER                        :: TotalNumberMPF, counter2, BCTotalNumberMPF
INTEGER,ALLOCATABLE            :: locnPart(:),offsetnPart(:),nPart_glob(:),minnParts(:), iPartCount(:)
INTEGER                        :: iPart, iSpec, counter
REAL,ALLOCATABLE               :: PartData(:,:)
INTEGER                        :: PartDataSize       !number of entries in each line of PartData
REAL                           :: TotalFlowrateMPF, RandVal, BCTotalFlowrateMPF
LOGICAL,ALLOCATABLE            :: PartDone(:)
!===================================================================================================================================
SWRITE(*,*) ' WRITE DSMCSurfCollis TO FILE...'

TypeString='DSMCSurfCollis'
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(TypeString),OutputTime))//'.h5'
PartDataSize=10
ALLOCATE(StrVarNames(PartDataSize))
StrVarNames(1)='ParticlePositionX'
StrVarNames(2)='ParticlePositionY'
StrVarNames(3)='ParticlePositionZ'
StrVarNames(4)='VelocityX'
StrVarNames(5)='VelocityY'
StrVarNames(6)='VelocityZ'
StrVarNames(7)='OldParticlePositionX'
StrVarNames(8)='OldParticlePositionY'
StrVarNames(9)='OldParticlePositionZ'
StrVarNames(10)='BCid'
ALLOCATE(locnPart(1:nSpecies) &
        ,offsetnPart(1:nSpecies) &
        ,nPart_glob(1:nSpecies) &
        ,minnParts(1:nSpecies) &
        ,iPartCount(1:nSpecies) )
#if USE_MPI
ALLOCATE(sendbuf(1:nSpecies) &
        ,recvbuf(1:nSpecies) )
#endif
ALLOCATE(SpeciesPositions( 1:nSpecies,1:MAXVAL(AnalyzeSurfCollis%Number(1:nSpecies)) ))

iPartCount(:)=0
DO iPart=1,AnalyzeSurfCollis%Number(nSpecies+1)
  IF (AnalyzeSurfCollis%Spec(iPart).LT.1 .OR. AnalyzeSurfCollis%Spec(iPart).GT.nSpecies) THEN
    CALL Abort(&
      __STAMP__,&
      'Error 1 in AnalyzeSurfCollis!')
  ELSE
    iPartCount(AnalyzeSurfCollis%Spec(iPart))=iPartCount(AnalyzeSurfCollis%Spec(iPart))+1
    SpeciesPositions(AnalyzeSurfCollis%Spec(iPart),iPartCount(AnalyzeSurfCollis%Spec(iPart)))=iPart
  END IF
END DO
DO iSpec=1,nSpecies
  locnPart(iSpec) = AnalyzeSurfCollis%Number(iSpec)
  IF (iPartCount(iSpec).NE.locnPart(iSpec)) CALL Abort(&
    __STAMP__,&
    'Error 2 in AnalyzeSurfCollis!')
END DO

#if USE_MPI
sendbuf(:)=locnPart(:)
recvbuf(:)=0
CALL MPI_EXSCAN(sendbuf,recvbuf,nSpecies,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,iError)
offsetnPart(:)=recvbuf(:)
sendbuf(:)=recvbuf(:)+locnPart(:)
CALL MPI_BCAST(sendbuf(:),nSpecies,MPI_INTEGER,nProcessors-1,MPI_COMM_WORLD,iError) !last proc knows global number
!global numbers
nPart_glob(:)=sendbuf(:)
DEALLOCATE(sendbuf &
          ,recvbuf )
!LOGWRITE(*,*)'offsetnPart,locnPart,nPart_glob',offsetnPart,locnPart,nPart_glob
CALL MPI_ALLREDUCE(locnPart(:),minnParts(:),nSpecies,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,IERROR)
IF (SFResampleAnalyzeSurfCollis) THEN
  CALL MPI_ALLGATHER(AnalyzeSurfCollis%Number(nSpecies+1), 1, MPI_INTEGER, globalNum, 1, MPI_INTEGER, MPI_COMM_WORLD, IERROR)
  TotalNumberMPF = SUM(globalNum)
ELSE
  CALL MPI_ALLREDUCE(AnalyzeSurfCollis%Number(nSpecies+1),TotalNumberMPF,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERROR)
END IF
#else
offsetnPart(:)=0
nPart_glob(:)=locnPart(:)
minnParts(:)=locnPart(:)
TotalNumberMPF=AnalyzeSurfCollis%Number(nSpecies+1)
#endif
! determine number of parts at BC of interest
BCTotalNumberMPF=0
IF (SFResampleAnalyzeSurfCollis) THEN
  DO iPart=1,AnalyzeSurfCollis%Number(nSpecies+1)
    IF (AnalyzeSurfCollis%BCid(iPart).LT.1 .OR. AnalyzeSurfCollis%BCid(iPart).GT.nPartBound) THEN
      CALL Abort(&
        __STAMP__,&
        'Error 3 in AnalyzeSurfCollis!')
    ELSE IF ( ANY(LastAnalyzeSurfCollis%BCs.EQ.0) .OR. ANY(LastAnalyzeSurfCollis%BCs.EQ.AnalyzeSurfCollis%BCid(iPart)) ) THEN
      BCTotalNumberMPF = BCTotalNumberMPF + 1
    END IF
  END DO
#if USE_MPI
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,BCTotalNumberMPF,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,iError)
#endif
  BCTotalFlowrateMPF=REAL(BCTotalNumberMPF)/TimeSample
END IF
TotalFlowrateMPF=REAL(TotalNumberMPF)/TimeSample

IF(MPIRoot) THEN !create File-Skeleton
  ! Create file
  CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)

  ! Write file header
  CALL WriteHDF5Header(TRIM(TypeString),File_ID)

  ! Write dataset properties "Time","VarNames","nSpecies","TotalFlowrateMPF"
  CALL WriteAttributeToHDF5(File_ID,'Time',1,RealScalar=OutputTime)
  CALL WriteAttributeToHDF5(File_ID,'VarNames',PartDataSize,StrArray=StrVarNames)
  CALL WriteAttributeToHDF5(File_ID,'NSpecies',1,IntegerScalar=nSpecies)
  CALL WriteAttributeToHDF5(File_ID,'TotalFlowrateMPF',1,RealScalar=TotalFlowrateMPF)

  CALL CloseDataFile()
END IF

#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif
CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_WORLD)

IF (SFResampleAnalyzeSurfCollis) THEN
  IF (LastAnalyzeSurfCollis%ReducePartNumber) THEN !reduce saved number of parts to MaxPartNumber
    LastAnalyzeSurfCollis%PartNumberSamp=MIN(BCTotalNumberMPF,LastAnalyzeSurfCollis%PartNumberReduced)
    ALLOCATE(PartDone(1:TotalNumberMPF))
    PartDone(:)=.FALSE.
  ELSE
    LastAnalyzeSurfCollis%PartNumberSamp=BCTotalNumberMPF
  END IF
  SWRITE(*,*) 'Number of saved particles for SFResampleAnalyzeSurfCollis: ',LastAnalyzeSurfCollis%PartNumberSamp
  SDEALLOCATE(LastAnalyzeSurfCollis%WallState)
  SDEALLOCATE(LastAnalyzeSurfCollis%Species)
  ALLOCATE(LastAnalyzeSurfCollis%WallState(6,LastAnalyzeSurfCollis%PartNumberSamp))
  ALLOCATE(LastAnalyzeSurfCollis%Species(LastAnalyzeSurfCollis%PartNumberSamp))
  LastAnalyzeSurfCollis%pushTimeStep = HUGE(LastAnalyzeSurfCollis%pushTimeStep)
#if USE_MPI
  IF (BCTotalNumberMPF.GT.0) THEN
    ALLOCATE(sendbuf2(1:AnalyzeSurfCollis%Number(nSpecies+1)*8))
    ALLOCATE(recvbuf2(1:TotalNumberMPF*8))
    ! Fill sendbufer
    counter2 = 0
    DO iPart=1,AnalyzeSurfCollis%Number(nSpecies+1)
      sendbuf2(counter2+1:counter2+6) = AnalyzeSurfCollis%Data(iPart,1:6)
      sendbuf2(counter2+7)           = REAL(AnalyzeSurfCollis%Spec(iPart))
      sendbuf2(counter2+8)           = REAL(AnalyzeSurfCollis%BCid(iPart))
      counter2 = counter2 + 8
    END DO
    ! Distribute particles to all procs
    counter2 = 0
    DO iProc = 0, nProcessors-1
      RecCount(iProc) = globalNum(iProc) * 8
      Displace(iProc) = counter2
      counter2 = counter2 + globalNum(iProc)*8
    END DO
    CALL MPI_ALLGATHERV(sendbuf2, 8*globalNum(myRank), MPI_DOUBLE_PRECISION, &
      recvbuf2, RecCount, Displace, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERROR)
    ! Add them to particle list
    counter2 = -8 !moved increment before usage, thus: -8 instead of 0
    DO counter = 1, LastAnalyzeSurfCollis%PartNumberSamp
      IF (LastAnalyzeSurfCollis%ReducePartNumber) THEN !reduce saved number of parts (differently in each proc. Could be changed)
        DO !get random (equal!) position between 8*[0,TotalNumberMPF-1] and accept if .NOT.PartDone and with right BC
          CALL RANDOM_NUMBER(RandVal)
          counter2 = MIN(1+INT(RandVal*REAL(TotalNumberMPF)),TotalNumberMPF) !( MIN(1+INT(RandVal*REAL(TotalNumberMPF)),TotalNumberMPF) - 1) *8
          IF (.NOT.PartDone(counter2) .AND. &
            ( ANY(LastAnalyzeSurfCollis%BCs.EQ.0) .OR. ANY(LastAnalyzeSurfCollis%BCs.EQ.INT(recvbuf2(8*counter2))) )) THEN
            PartDone(counter2)=.TRUE.
            counter2 = 8*(counter2-1)
            EXIT
          END IF
        END DO
      ELSE
        counter2 = counter2 + 8
      END IF
      LastAnalyzeSurfCollis%WallState(:,counter) = recvbuf2(counter2+1:counter2+6)
      LastAnalyzeSurfCollis%Species(counter) = INT(recvbuf2(counter2+7))
      IF (ANY(LastAnalyzeSurfCollis%SpeciesForDtCalc.EQ.0) .OR. &
          ANY(LastAnalyzeSurfCollis%SpeciesForDtCalc.EQ.LastAnalyzeSurfCollis%Species(counter))) &
        LastAnalyzeSurfCollis%pushTimeStep = MIN( LastAnalyzeSurfCollis%pushTimeStep &
        , DOT_PRODUCT(LastAnalyzeSurfCollis%NormVecOfWall,LastAnalyzeSurfCollis%WallState(4:6,counter)) )
    END DO
    DEALLOCATE(sendbuf2 &
              ,recvbuf2 )
  END IF
#else
  ! Add particle to list
  counter2 = 0
  DO counter = 1, LastAnalyzeSurfCollis%PartNumberSamp
    IF (LastAnalyzeSurfCollis%ReducePartNumber) THEN !reduce saved number of parts (differently for each proc. Could be changed)
      DO !get random (equal!) position between [1,TotalNumberMPF] and accept if .NOT.PartDone and with right BC
        CALL RANDOM_NUMBER(RandVal)
        counter2 = MIN(1+INT(RandVal*REAL(TotalNumberMPF)),TotalNumberMPF)
        IF (.NOT.PartDone(counter2) .AND. &
          ( ANY(LastAnalyzeSurfCollis%BCs.EQ.0) .OR. ANY(LastAnalyzeSurfCollis%BCs.EQ.AnalyzeSurfCollis%BCid(counter2)) )) THEN
          PartDone(counter2)=.TRUE.
          EXIT
        END IF
      END DO
    ELSE
      counter2 = counter2 + 1
    END IF
    LastAnalyzeSurfCollis%WallState(:,counter) = AnalyzeSurfCollis%Data(counter2,1:6)
    LastAnalyzeSurfCollis%Species(counter) = AnalyzeSurfCollis%Spec(counter2)
    IF (ANY(LastAnalyzeSurfCollis%SpeciesForDtCalc.EQ.0) .OR. &
        ANY(LastAnalyzeSurfCollis%SpeciesForDtCalc.EQ.LastAnalyzeSurfCollis%Species(counter))) &
      LastAnalyzeSurfCollis%pushTimeStep = MIN( LastAnalyzeSurfCollis%pushTimeStep &
      , DOT_PRODUCT(LastAnalyzeSurfCollis%NormVecOfWall,LastAnalyzeSurfCollis%WallState(4:6,counter)) )
  END DO
#endif
  IF (LastAnalyzeSurfCollis%pushTimeStep .LE. 0.) THEN
    CALL Abort(&
      __STAMP__,&
      'Error with SFResampleAnalyzeSurfCollis. Something is wrong with velocities or NormVecOfWall!',&
      999,LastAnalyzeSurfCollis%pushTimeStep)
  ELSE
    LastAnalyzeSurfCollis%pushTimeStep = r_SF / LastAnalyzeSurfCollis%pushTimeStep !dt required for smallest projected velo to cross r_SF
    LastAnalyzeSurfCollis%PartNumberDepo = NINT(BCTotalFlowrateMPF * LastAnalyzeSurfCollis%pushTimeStep)
    SWRITE(*,'(A,E12.5,x,I0)') 'Total Flowrate and to be inserted number of MP for SFResampleAnalyzeSurfCollis: ' &
      ,BCTotalFlowrateMPF, LastAnalyzeSurfCollis%PartNumberDepo
    IF (LastAnalyzeSurfCollis%PartNumberDepo .GT. LastAnalyzeSurfCollis%PartNumberSamp) THEN
      SWRITE(*,*) 'WARNING: PartNumberDepo .GT. PartNumberSamp!'
    END IF
    IF (LastAnalyzeSurfCollis%PartNumberDepo .GT. LastAnalyzeSurfCollis%PartNumThreshold) THEN
      CALL Abort(&
        __STAMP__,&
        'Error with SFResampleAnalyzeSurfCollis: PartNumberDepo .gt. PartNumThreshold',&
        LastAnalyzeSurfCollis%PartNumberDepo,r_SF/LastAnalyzeSurfCollis%pushTimeStep)
    END IF
  END IF
END IF !SFResampleAnalyzeSurfCollis

DO iSpec=1,nSpecies
  ALLOCATE(PartData(PartDataSize,offsetnPart(iSpec)+1:offsetnPart(iSpec)+locnPart(iSpec)))
  DO iPart=1,locnPart(iSpec)
    PartData(1,offsetnPart(iSpec)+iPart)=AnalyzeSurfCollis%Data(SpeciesPositions(iSpec,iPart),1)
    PartData(2,offsetnPart(iSpec)+iPart)=AnalyzeSurfCollis%Data(SpeciesPositions(iSpec,iPart),2)
    PartData(3,offsetnPart(iSpec)+iPart)=AnalyzeSurfCollis%Data(SpeciesPositions(iSpec,iPart),3)
    PartData(4,offsetnPart(iSpec)+iPart)=AnalyzeSurfCollis%Data(SpeciesPositions(iSpec,iPart),4)
    PartData(5,offsetnPart(iSpec)+iPart)=AnalyzeSurfCollis%Data(SpeciesPositions(iSpec,iPart),5)
    PartData(6,offsetnPart(iSpec)+iPart)=AnalyzeSurfCollis%Data(SpeciesPositions(iSpec,iPart),6)
    PartData(7,offsetnPart(iSpec)+iPart)=AnalyzeSurfCollis%Data(SpeciesPositions(iSpec,iPart),7)
    PartData(8,offsetnPart(iSpec)+iPart)=AnalyzeSurfCollis%Data(SpeciesPositions(iSpec,iPart),8)
    PartData(9,offsetnPart(iSpec)+iPart)=AnalyzeSurfCollis%Data(SpeciesPositions(iSpec,iPart),9)
    PartData(10,offsetnPart(iSpec)+iPart)=REAL(AnalyzeSurfCollis%BCid(SpeciesPositions(iSpec,iPart)))
  END DO
  WRITE(H5_Name,'(A,I3.3)') 'SurfCollisData_Spec',iSpec

  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (&
        nPart_glob   => INT(nPart_glob(iSpec),IK)  ,&
        PartDataSize => INT(PartDataSize,IK)       ,&
        locnPart     => INT(locnPart(iSpec),IK)           ,&
        offsetnPart  => INT(offsetnPart(iSpec),IK) )
    IF(minnParts(iSpec).EQ.0)THEN
      CALL WriteArrayToHDF5(DataSetName=TRIM(H5_Name)  , rank=2        , &
                            nValGlobal =(/PartDataSize , nPart_glob /) , &
                            nVal       =(/PartDataSize , locnPart   /) , &
                            offset     =(/0_IK         , offsetnPart/) , &
                            collective =.FALSE.        , RealArray=PartData)
    ELSE
      CALL WriteArrayToHDF5(DataSetName=TRIM(H5_Name)  , rank=2        , &
                            nValGlobal =(/PartDataSize , nPart_glob /) , &
                            nVal       =(/PartDataSize , locnPart   /) , &
                            offset     =(/0_IK         , offsetnPart/) , &
                            collective =.TRUE.         , RealArray=PartData)
    END IF
  END ASSOCIATE
  DEALLOCATE(PartData)
END DO !iSpec

CALL CloseDataFile()
DEALLOCATE(locnPart &
          ,offsetnPart &
          ,nPart_glob &
          ,minnParts &
          ,iPartCount )
DEALLOCATE(SpeciesPositions)
DEALLOCATE(StrVarNames)

END SUBROUTINE WriteAnalyzeSurfCollisToHDF5

SUBROUTINE SamplingRotVibRelaxProb(iElem)
!===================================================================================================================================
!> Update sampling arrays for rotational and vibrational relaxation probability for DSMC Output
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_Particle_Vars          ,ONLY: nSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: MeanProb, MaxProb, PartNum
INTEGER                       :: iSpec
!===================================================================================================================================
  IF(DSMC%RotRelaxProb.EQ.2) THEN
    MeanProb = 0
    MaxProb = 0
    PartNum = 0.
    DO iSpec=1,nSpecies
      IF(DSMC%CalcRotProb(iSpec,3).GT.0) THEN
        DSMC%QualityFacSampRot(iElem,iSpec,2) = DSMC%QualityFacSampRot(iElem,iSpec,2) + DSMC%CalcRotProb(iSpec,2)
        DSMC%QualityFacSampRot(iElem,iSpec,1) = DSMC%QualityFacSampRot(iElem,iSpec,1) &
                                              + DSMC%CalcRotProb(iSpec,1) / DSMC%CalcRotProb(iSpec,3)
        DSMC%QualityFacSampRotSamp(iElem,iSpec) = DSMC%QualityFacSampRotSamp(iElem,iSpec) + 1
        MaxProb = MAX(MaxProb,DSMC%CalcRotProb(iSpec,2))
        MeanProb = MeanProb + DSMC%CalcRotProb(iSpec,1) * DSMC%CalcRotProb(iSpec,3)
        PartNum = PartNum + DSMC%CalcRotProb(iSpec,3)
      END IF
    END DO
    IF((nSpecies.GT.1).AND.(PartNum.GT.0)) THEN
      DSMC%QualityFacSampRot(iElem,nSpecies+1,2) = DSMC%QualityFacSampRot(iElem,nSpecies+1,2) + MaxProb
      DSMC%QualityFacSampRot(iElem,nSpecies+1,1) = DSMC%QualityFacSampRot(iElem,nSpecies+1,1) + MeanProb / PartNum
      DSMC%QualityFacSampRotSamp(iElem,nSpecies+1) = DSMC%QualityFacSampRotSamp(iElem,nSpecies+1) + 1
    END IF
  END IF
  ! Sample vibration relaxation probability
  IF(DSMC%VibRelaxProb.EQ.2) THEN
    MeanProb = 0
    MaxProb = 0
    PartNum = 0
    DO iSpec=1,nSpecies
      IF(DSMC%CalcVibProb(iSpec,2).GT.0) THEN
        DSMC%QualityFacSampVib(iElem,iSpec,2) = DSMC%QualityFacSampVib(iElem,iSpec,2) + DSMC%CalcVibProb(iSpec,2)
        DSMC%QualityFacSampVibSamp(iElem,iSpec,2) = DSMC%QualityFacSampVibSamp(iElem,iSpec,2) + 1
        MaxProb = MAX(MaxProb,DSMC%CalcVibProb(iSpec,2))
      END IF
      IF(DSMC%CalcVibProb(iSpec,3).GT.0) THEN
        DSMC%QualityFacSampVib(iElem,iSpec,1) = DSMC%QualityFacSampVib(iElem,iSpec,1) &
                                              + DSMC%CalcVibProb(iSpec,1) / DSMC%CalcVibProb(iSpec,3)
        DSMC%QualityFacSampVibSamp(iElem,iSpec,1) = DSMC%QualityFacSampVibSamp(iElem,iSpec,1) + 1
        MeanProb = MeanProb + DSMC%CalcVibProb(iSpec,1) * DSMC%CalcVibProb(iSpec,3)
        PartNum = PartNum + DSMC%CalcVibProb(iSpec,3)
      END IF
    END DO
    IF(MaxProb.GT.0.) THEN
      DSMC%QualityFacSampVib(iElem,nSpecies+1,2) = DSMC%QualityFacSampVib(iElem,nSpecies+1,2) + MaxProb
      DSMC%QualityFacSampVibSamp(iElem,nSpecies+1,2) = DSMC%QualityFacSampVibSamp(iElem,nSpecies+1,2) + 1
    END IF
    IF((nSpecies.GT.1).AND.(PartNum.GT.0)) THEN
      DSMC%QualityFacSampVib(iElem,nSpecies+1,1) = DSMC%QualityFacSampVib(iElem,nSpecies+1,1) + MeanProb / PartNum
      DSMC%QualityFacSampVibSamp(iElem,nSpecies+1,1) = DSMC%QualityFacSampVibSamp(iElem,nSpecies+1,1) + 1
    END IF
  END IF
END SUBROUTINE SamplingRotVibRelaxProb


SUBROUTINE SummarizeQualityFactors(iElem)
!===================================================================================================================================
!> Sample quality factors. MCS over MFP, max and mean collision probability
!===================================================================================================================================
! MODULES
USE MOD_TimeDisc_Vars         ,ONLY: time, TEnd
USE MOD_Globals
USE MOD_DSMC_Vars             ,ONLY: DSMC
USE MOD_Particle_Vars         ,ONLY: WriteMacroVolumeValues
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)    :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
  ! mean collision probability of all collision pairs
  IF(DSMC%CollProbMeanCount.GT.0) THEN
    DSMC%QualityFacSamp(iElem,1) = DSMC%QualityFacSamp(iElem,1) + DSMC%CollProbMax
    DSMC%QualityFacSamp(iElem,2) = DSMC%QualityFacSamp(iElem,2) + DSMC%CollProbMean / REAL(DSMC%CollProbMeanCount)
  END IF
  ! mean collision separation distance of actual collisions
  IF(DSMC%CollSepCount.GT.0) DSMC%QualityFacSamp(iElem,3) = DSMC%QualityFacSamp(iElem,3) + DSMC%MCSoverMFP
  ! Counting sample size
  DSMC%QualityFacSamp(iElem,4) = DSMC%QualityFacSamp(iElem,4) + 1.
  ! Sample rotation relaxation probability
  IF((DSMC%RotRelaxProb.EQ.2).OR.(DSMC%VibRelaxProb.EQ.2)) CALL SamplingRotVibRelaxProb(iElem)
END IF
! mean collision separation distance of actual collisions
IF(DSMC%CollSepCount.GT.0) DSMC%QualityFacSamp(iElem,3) = DSMC%QualityFacSamp(iElem,3) + DSMC%MCSoverMFP
! Counting sample size
DSMC%QualityFacSamp(iElem,4) = DSMC%QualityFacSamp(iElem,4) + 1.

END SUBROUTINE SummarizeQualityFactors


SUBROUTINE DSMCMacroSampling()
!===================================================================================================================================
!> Check if sampling should be activated and perform sampling
!===================================================================================================================================
! MODULES
USE MOD_TimeDisc_Vars         ,ONLY: time, TEnd
USE MOD_Globals
USE MOD_DSMC_Vars             ,ONLY: DSMC
USE MOD_DSMC_SteadyState      ,ONLY: QCrit_evaluation, SteadyStateDetection_main
USE MOD_Restart_Vars          ,ONLY: RestartTime
USE MOD_Mesh_Vars             ,ONLY: MeshFile
USE MOD_TimeDisc_Vars         ,ONLY: iter
USE MOD_DSMC_Vars             ,ONLY: UseQCrit, SamplingActive, QCritTestStep, QCritLastTest, UseSSD
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: nOutput
!-----------------------------------------------------------------------------------------------------------------------------------

#if (PP_TimeDiscMethod==42)
! Do not perform sampling in the case of a reservoir simulation
IF (DSMC%ReservoirSimu) RETURN
#endif

IF(UseQCrit) THEN
  ! Use QCriterion (Burt,Boyd) for steady - state detection
  IF((.NOT.SamplingActive).AND.(iter-QCritLastTest.EQ.QCritTestStep)) THEN
    CALL QCrit_evaluation()
    QCritLastTest=iter
    IF(SamplingActive) THEN
      SWRITE(*,*)'Sampling active'
      ! Set TimeFracSamp and DeltaTimeOutput -> correct number of outputs
      DSMC%TimeFracSamp = (TEnd-Time)/TEnd
      DSMC%DeltaTimeOutput = (DSMC%TimeFracSamp * TEnd) / REAL(DSMC%NumOutput)
    ENDIF
  ENDIF
ELSEIF(UseSSD) THEN
  ! Use SSD for steady - state detection
  IF((.NOT.SamplingActive)) THEN
    CALL SteadyStateDetection_main()
    IF(SamplingActive) THEN
      SWRITE(*,*)'Sampling active'
      ! Set TimeFracSamp and DeltaTimeOutput -> correct number of outputs
      DSMC%TimeFracSamp = (TEnd-Time)/TEnd
      DSMC%DeltaTimeOutput = (DSMC%TimeFracSamp * TEnd) / REAL(DSMC%NumOutput)
    ENDIF
  ENDIF
ELSE
  ! Use user given TimeFracSamp
  IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).AND.(.NOT.SamplingActive))  THEN
    SamplingActive=.TRUE.
    SWRITE(*,*)'Sampling active'
  ENDIF
ENDIF
!
! Calculate Entropy using Theorem of Boltzmann
!CALL EntropyCalculation()
!

IF(SamplingActive) THEN
  CALL DSMC_data_sampling()
  IF(DSMC%NumOutput.NE.0) THEN
    nOutput = INT((DSMC%TimeFracSamp * TEnd)/DSMC%DeltaTimeOutput)-DSMC%NumOutput + 1
    IF(Time.GE.((1-DSMC%TimeFracSamp)*TEnd + DSMC%DeltaTimeOutput * nOutput)) THEN
      DSMC%NumOutput = DSMC%NumOutput - 1
      ! Skipping outputs immediately after the first few iterations
      IF(RestartTime.LT.((1-DSMC%TimeFracSamp)*TEnd + DSMC%DeltaTimeOutput * REAL(nOutput))) THEN
        CALL WriteDSMCToHDF5(TRIM(MeshFile),time)
        IF(DSMC%CalcSurfaceVal) CALL CalcSurfaceValues(during_dt_opt=.TRUE.)
      END IF
    END IF
  END IF
END IF

END SUBROUTINE DSMCMacroSampling

END MODULE MOD_DSMC_Analyze
