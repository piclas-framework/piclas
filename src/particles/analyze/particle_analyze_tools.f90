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

MODULE MOD_Particle_Analyze_Tools
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
#ifdef PARTICLES
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE CalcEkinPart
  MODULE PROCEDURE CalcEkinPart
END INTERFACE

INTERFACE CalcEkinPart2
  MODULE PROCEDURE CalcEkinPart2
END INTERFACE

INTERFACE CalcNumPartsOfSpec
  MODULE PROCEDURE CalcNumPartsOfSpec
END INTERFACE

INTERFACE CalcMixtureTemp
  MODULE PROCEDURE CalcMixtureTemp
END INTERFACE

PUBLIC :: CalcEkinPart,CalcEkinPart2
PUBLIC :: CalcNumPartsOfSpec
PUBLIC :: AllocateElectronIonDensityCell,AllocateElectronTemperatureCell,AllocateCalcElectronEnergy
PUBLIC :: CalculateElectronIonDensityCell,CalculateElectronTemperatureCell
PUBLIC :: CalcShapeEfficiencyR
PUBLIC :: CalcKineticEnergy
PUBLIC :: CalcKineticEnergyAndMaximum
PUBLIC :: CalcNumberDensity
PUBLIC :: CalcSurfaceFluxInfo
PUBLIC :: CalcTransTemp
PUBLIC :: CalcMixtureTemp
PUBLIC :: CalcTelec,CalcTVibPoly
#if (PP_TimeDiscMethod==2 || PP_TimeDiscMethod==4 || PP_TimeDiscMethod==300 || PP_TimeDiscMethod==400 || (PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=509) || PP_TimeDiscMethod==120)
PUBLIC :: CalcRelaxProbRotVib
#endif
PUBLIC :: CalcVelocities
#if (PP_TimeDiscMethod==4)
PUBLIC :: CollRates,CalcRelaxRates,CalcRelaxRatesElec,ReacRates
#endif /*(PP_TimeDiscMethod==4)*/
PUBLIC :: CalcPowerDensity
PUBLIC :: CalculatePartElemData
PUBLIC :: CalcCoupledPowerPart, CalcEelec
PUBLIC :: CalcNumberDensityBGGasDistri
#if USE_HDG
PUBLIC :: CalculatePCouplElectricPotential
#endif /*USE_HDG*/
!===================================================================================================================================

CONTAINS


PPURE REAL FUNCTION CalcEkinPart(iPart)
!===================================================================================================================================
! computes the kinetic energy of one particle
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Globals_Vars  ,ONLY: c2, c2_inv, RelativisticLimit
USE MOD_Particle_Vars ,ONLY: PartState, PartSpecies, Species
USE MOD_PARTICLE_Vars ,ONLY: usevMPF
USE MOD_Particle_Vars ,ONLY: PartLorentzType
USE MOD_DSMC_Vars     ,ONLY: RadialWeighting
USE MOD_part_tools    ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(IN)                 :: iPart
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                               :: partV2, gamma1, WeightingFactor
!===================================================================================================================================

IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
  WeightingFactor = GetParticleWeight(iPart)
ELSE
  WeightingFactor = GetParticleWeight(iPart) * Species(PartSpecies(iPart))%MacroParticleFactor
END IF

IF (PartLorentzType.EQ.5)THEN
  ! gamma v is pushed instead of gamma, therefore, only the relativistic kinetic energy is computed
  ! compute gamma
  gamma1=SQRT(1.0+DOTPRODUCT(PartState(4:6,iPart))*c2_inv)
  CalcEkinPart=(gamma1-1.0)*Species(PartSpecies(iPart))%MassIC*c2 * WeightingFactor
ELSE
  partV2 = DOTPRODUCT(PartState(4:6,iPart))
  IF (partV2.LT.RelativisticLimit)THEN ! |v| < 1000000 when speed of light is 299792458
    CalcEkinPart= 0.5 * Species(PartSpecies(iPart))%MassIC * partV2 * WeightingFactor
  ELSE
    gamma1=partV2*c2_inv
    ! Sanity check: Lorentz factor must be below 1.0
    IF(gamma1.GE.1.0)THEN
      CalcEkinPart=-1.0
    ELSE
      gamma1=1.0/SQRT(1.-gamma1)
      CalcEkinPart=(gamma1-1.0)*Species(PartSpecies(iPart))%MassIC*c2 * WeightingFactor
    END IF ! gamma1.GE.1.0
  END IF ! ipartV2
END IF
END FUNCTION CalcEkinPart


PPURE REAL FUNCTION CalcEkinPart2(velocity,Species_IN,WeightingFactor)
!===================================================================================================================================
! computes the kinetic energy of one particle given its velocity, species and weighting factor
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Globals_Vars  ,ONLY: c2,c2_inv,RelativisticLimit
USE MOD_Particle_Vars ,ONLY: Species
USE MOD_Particle_Vars ,ONLY: PartLorentzType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: velocity(1:3)
INTEGER,INTENT(IN)              :: Species_IN
REAL,INTENT(IN)                 :: WeightingFactor
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: partV2, gamma1
!===================================================================================================================================
partV2 = DOT_PRODUCT(velocity,velocity)

IF (PartLorentzType.EQ.5)THEN
  ! gamma v is pushed instead of gamma, therefore, only the relativistic kinetic energy is computed
  ! compute gamma
  gamma1=SQRT(1.0+partV2*c2_inv)
  CalcEkinPart2=(gamma1-1.0)*Species(Species_IN)%MassIC*c2 * WeightingFactor
ELSE
  IF (partV2.LT.RelativisticLimit)THEN ! |v| < 1000000 when speed of light is 299792458
    CalcEkinPart2= 0.5 * Species(Species_IN)%MassIC * partV2 * WeightingFactor
  ELSE
    gamma1=partV2*c2_inv
    ! Sanity check: Lorentz factor must be below 1.0
    IF(gamma1.GE.1.0)THEN
      CalcEkinPart2=-1.0
    ELSE
      gamma1=1.0/SQRT(1.-gamma1)
      CalcEkinPart2=(gamma1-1.0)*Species(Species_IN)%MassIC*c2 * WeightingFactor
    END IF ! gamma1.GE.1.0
  END IF ! ipartV2
END IF
END FUNCTION CalcEkinPart2


SUBROUTINE CalcNumPartsOfSpec(NumSpec,SimNumSpec,CalcNumSpec_IN,CalcSimNumSpec_IN)
!===================================================================================================================================
! Computes the number of simulated particles AND number of real particles within the domain
! Last section of the routine contains the MPI-communication
! CAUTION: SimNumSpec equals NumSpec only for constant weighting factor
! NOTE: Background gas particles are not considered
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Particle_Vars         ,ONLY: PDM,PartSpecies
USE MOD_Particle_Analyze_Vars ,ONLY: nSpecAnalyze
USE MOD_DSMC_Vars             ,ONLY: BGGas, DSMC
USE MOD_Particle_Vars         ,ONLY: nSpecies, Species
USE MOD_part_tools            ,ONLY: GetParticleWeight
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(OUT)                   :: NumSpec(nSpecAnalyze)
INTEGER(KIND=IK)                   :: SimNumSpecMin,SimNumSpecMax
INTEGER(KIND=IK),INTENT(OUT)       :: SimNumSpec(nSpecAnalyze)
LOGICAL,INTENT(IN)                 :: CalcNumSpec_IN,CalcSimNumSpec_IN ! Flags for performing MPI reduce
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iPart, iSpec
!===================================================================================================================================

NumSpec    = 0.
SimNumSpec = 0
DO iPart=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    ! NumSpec = real particle number
    NumSpec(PartSpecies(iPart))    = NumSpec(PartSpecies(iPart)) + GetParticleWeight(iPart)
    ! SimNumSpec = simulated particle number
    SimNumSpec(PartSpecies(iPart)) = SimNumSpec(PartSpecies(iPart)) + 1
  END IF
END DO
IF(BGGas%NumberOfSpecies.GT.0) THEN
  DO iSpec = 1, nSpecies
    IF(BGGas%BackgroundSpecies(iSpec)) THEN
      NumSpec(iSpec) = 0.
      SimNumSpec(iSpec) = 0
    END IF
  END DO
END IF
IF (DSMC%DoAmbipolarDiff) THEN
  NumSpec(DSMC%AmbiDiffElecSpec) = 0.0
  SimNumSpec(DSMC%AmbiDiffElecSpec) = 0.0
  DO iSpec = 1, nSpecies
    IF (iSpec.EQ.DSMC%AmbiDiffElecSpec) CYCLE
    IF(Species(iSpec)%ChargeIC.GT.0.0) THEN
      NumSpec(DSMC%AmbiDiffElecSpec) = NumSpec(DSMC%AmbiDiffElecSpec) + NumSpec(iSpec)
      SimNumSpec(DSMC%AmbiDiffElecSpec) = SimNumSpec(DSMC%AmbiDiffElecSpec) + SimNumSpec(iSpec)
    END IF
  END DO
END IF

IF(nSpecAnalyze.GT.1)THEN
  NumSpec(nSpecAnalyze)    = SUM(NumSpec(1:nSpecies))
  SimNumSpec(nSpecAnalyze) = SUM(SimNumSpec(1:nSpecies))
END IF

#if USE_MPI
IF (MPIRoot) THEN
  IF(CalcNumSpec_IN)    CALL MPI_REDUCE(MPI_IN_PLACE    , NumSpec       , nSpecAnalyze , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
  IF(CalcSimNumSpec_IN) CALL MPI_REDUCE(SUM(SimNumSpec) , SimNumSpecMin , 1            , MPI_INTEGER_INT_KIND , MPI_MIN , 0 , MPI_COMM_PICLAS , IERROR)
  IF(CalcSimNumSpec_IN) CALL MPI_REDUCE(SUM(SimNumSpec) , SimNumSpecMax , 1            , MPI_INTEGER_INT_KIND , MPI_MAX , 0 , MPI_COMM_PICLAS , IERROR)
  IF(CalcSimNumSpec_IN) CALL MPI_REDUCE(MPI_IN_PLACE    , SimNumSpec    , nSpecAnalyze , MPI_INTEGER_INT_KIND , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
ELSE
  IF(CalcNumSpec_IN)    CALL MPI_REDUCE(NumSpec         , 0 , nSpecAnalyze , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
  IF(CalcSimNumSpec_IN) CALL MPI_REDUCE(SUM(SimNumSpec) , 0 , 1            , MPI_INTEGER_INT_KIND , MPI_MIN , 0 , MPI_COMM_PICLAS , IERROR)
  IF(CalcSimNumSpec_IN) CALL MPI_REDUCE(SUM(SimNumSpec) , 0 , 1            , MPI_INTEGER_INT_KIND , MPI_MAX , 0 , MPI_COMM_PICLAS , IERROR)
  IF(CalcSimNumSpec_IN) CALL MPI_REDUCE(SimNumSpec      , 0 , nSpecAnalyze , MPI_INTEGER_INT_KIND , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
END IF
#endif /*USE_MPI*/

! Set global number of particles (info for std out)
IF(CalcSimNumSpec_IN)THEN
  GlobalNbrOfParticlesUpdated = .TRUE.
#if USE_MPI
  IF(MPIRoot)THEN
#endif /*USE_MPI*/
    nGlobalNbrOfParticles(1) = INT(SimNumSpecMin,KIND=IK)
    nGlobalNbrOfParticles(2) = INT(SimNumSpecMax,KIND=IK)
    nGlobalNbrOfParticles(3) = INT(SimNumSpec(nSpecAnalyze),KIND=IK)
    nGlobalNbrOfParticles(4) = MIN(nGlobalNbrOfParticles(1),nGlobalNbrOfParticles(4))
    nGlobalNbrOfParticles(5) = MAX(nGlobalNbrOfParticles(2),nGlobalNbrOfParticles(5))
    nGlobalNbrOfParticles(6) = MAX(nGlobalNbrOfParticles(3),nGlobalNbrOfParticles(6))
#if USE_MPI
  END IF ! MPIRoot
#endif /*USE_MPI*/
END IF ! CalcSimNumSpec_IN

END SUBROUTINE CalcNumPartsOfSpec


!===================================================================================================================================
!> Allocate the required ElemData arrays for electron/ion/neutral density for each element (singly data point in each element),
!> which is required for BR electron model (conversion from BR to fully kinetic)
!===================================================================================================================================
SUBROUTINE AllocateElectronIonDensityCell()
! MODULES
USE MOD_IO_HDF5               ,ONLY: AddToElemData,ElementOut
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars ,ONLY: ElectronDensityCell,IonDensityCell,NeutralDensityCell,ChargeNumberCell
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

!===================================================================================================================================
! Check if ElectronDensityCell has already been allocated (CalcElectronIonDensity=F and UseBRElectronFluid=T)
IF(ALLOCATED(ElectronDensityCell)) RETURN

! electrons
ALLOCATE( ElectronDensityCell(1:PP_nElems) )
ElectronDensityCell=0.0
CALL AddToElemData(ElementOut,'ElectronDensityCell',RealArray=ElectronDensityCell(1:PP_nElems))

! ions
ALLOCATE( IonDensityCell(1:PP_nElems) )
IonDensityCell=0.0
CALL AddToElemData(ElementOut,'IonDensityCell',RealArray=IonDensityCell(1:PP_nElems))

! neutrals
ALLOCATE( NeutralDensityCell(1:PP_nElems) )
NeutralDensityCell=0.0
CALL AddToElemData(ElementOut,'NeutralDensityCell',RealArray=NeutralDensityCell(1:PP_nElems))

! charge number
ALLOCATE( ChargeNumberCell(1:PP_nElems) )
ChargeNumberCell=0.0
CALL AddToElemData(ElementOut,'ChargeNumberCell',RealArray=ChargeNumberCell(1:PP_nElems))

END SUBROUTINE AllocateElectronIonDensityCell


!===================================================================================================================================
!> Allocate the required ElemData arrays for electron/ion/neutral density for each element (singly data point in each element),
!> which is required for BR electron model (conversion from BR to fully kinetic)
!===================================================================================================================================
SUBROUTINE AllocateElectronTemperatureCell()
! MODULES
USE MOD_IO_HDF5               ,ONLY: AddToElemData,ElementOut
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars ,ONLY: ElectronTemperatureCell
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

!===================================================================================================================================
! Check if ElectronTemperatureCell has already been allocated (CalcElectronIonDensity=F and UseBRElectronFluid=T)
IF(ALLOCATED(ElectronTemperatureCell)) RETURN

ALLOCATE( ElectronTemperatureCell(1:PP_nElems) )
ElectronTemperatureCell=0.0
CALL AddToElemData(ElementOut,'ElectronTemperatureCell',RealArray=ElectronTemperatureCell(1:PP_nElems))

END SUBROUTINE AllocateElectronTemperatureCell


!===================================================================================================================================
!> Allocate the required ElemData arrays for electron min/max/average energy for each element (singly data point in each element)
!> The energy is given in electron volt (eV)
!===================================================================================================================================
SUBROUTINE AllocateCalcElectronEnergy()
! MODULES
USE MOD_IO_HDF5               ,ONLY: AddToElemData,ElementOut
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars ,ONLY: ElectronMinEnergyCell,ElectronMaxEnergyCell,ElectronAverageEnergyCell
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

!===================================================================================================================================
! Check if ElectronMinEnergyCell has already been allocated (CalcElectronEnergy=T)
IF(ALLOCATED(ElectronMinEnergyCell)) RETURN

ALLOCATE(ElectronMinEnergyCell(1:PP_nElems))
ElectronMinEnergyCell=0.
CALL AddToElemData(ElementOut,'ElectronMinEnergyCell[eV]',RealArray=ElectronMinEnergyCell(1:PP_nElems))
ALLOCATE(ElectronMaxEnergyCell(1:PP_nElems))
ElectronMaxEnergyCell=0.
CALL AddToElemData(ElementOut,'ElectronMaxEnergyCell[eV]',RealArray=ElectronMaxEnergyCell(1:PP_nElems))
ALLOCATE(ElectronAverageEnergyCell(1:PP_nElems))
ElectronAverageEnergyCell=0.
CALL AddToElemData(ElementOut,'ElectronAverageEnergyCell[eV]',RealArray=ElectronAverageEnergyCell(1:PP_nElems))

END SUBROUTINE AllocateCalcElectronEnergy


SUBROUTINE CalculatePartElemData()
!===================================================================================================================================
! use the plasma frequency per cell to estimate the pic time step
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Analyze_Vars  ,ONLY: CalcPlasmaFrequency,CalcPICTimeStep,CalcElectronIonDensity,CalcPICTimeStepCyclotron
USE MOD_Particle_Analyze_Vars  ,ONLY: CalcElectronTemperature,CalcDebyeLength,CalcIonizationDegree,CalcPointsPerDebyeLength
USE MOD_Particle_Analyze_Vars  ,ONLY: CalcPlasmaParameter,CalcPICCFLCondition,CalcMaxPartDisplacement,CalcElectronEnergy
USE MOD_Particle_Analyze_Vars  ,ONLY: CalcCyclotronFrequency
#if USE_MPI
USE MOD_Globals
USE MOD_Particle_Analyze_Vars ,ONLY: PPDCellResolved,PICTimeCellResolved,PICValidPlasmaCellSum
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
INTEGER :: tmpArray(1:6)
#endif /*USE_MPI*/
!===================================================================================================================================

! electron density
IF(CalcElectronIonDensity) CALL CalculateElectronIonDensityCell()

! Ionization degree: n_i / (n_i + n_n)
! ion density versus sum of ion and neutral density
IF(CalcIonizationDegree) CALL CalculateIonizationCell()

! electron temperature
IF(CalcElectronTemperature) CALL CalculateElectronTemperatureCell()

! electron energies
IF(CalcElectronEnergy) CALL CalculateElectronEnergyCell()

! plasma frequency
IF(CalcPlasmaFrequency) CALL CalculatePlasmaFrequencyCell()

! Cyclotron frequency
IF(CalcCyclotronFrequency) CALL CalculateCyclotronFrequencyAndRadiusCell()

! PIC time step for gyro motion (cyclotron frequency)
IF(CalcPICTimeStepCyclotron) CALL CalculatePICTimeStepCyclotron()

! Debye length
IF(CalcDebyeLength) CALL CalculateDebyeLengthCell()

! Plasma parameter: 4/3 * pi * n_e * lambda_D^3
IF(CalcPlasmaParameter) CALL CalculatePlasmaParameter()

! PIC time step (plasma frequency)
IF(CalcPICTimeStep) CALL CalculatePICTimeStepCell()

! PointsPerDebyeLength: PPD = (p+1)*lambda_D/L_cell
IF(CalcPointsPerDebyeLength) CALL CalculatePPDCell()

! PICCFL Condition: PICCFL = dt/L_cell * SQRT( kB*Te/me )
IF(CalcPICCFLCondition) CALL CalculatePICCFL()

! Compute the maximum displacement of the fastest particle
! MaxPartDisplacement = max(v_iPart)*dT/L_cell <  1.0
IF(CalcMaxPartDisplacement) CALL CalculateMaxPartDisplacement()

! Communicate data
#if USE_MPI
IF(CalcPointsPerDebyeLength.OR.CalcPICTimeStep)THEN
  tmpArray = 0
  IF(CalcPointsPerDebyeLength)THEN
    tmpArray(1) = PPDCellResolved(1)
    tmpArray(2) = PPDCellResolved(2)
    tmpArray(3) = PPDCellResolved(3)
    tmpArray(4) = PPDCellResolved(4)
  END IF ! CalcPointsPerDebyeLength
  IF(CalcPICTimeStep)THEN
    tmpArray(5) = PICTimeCellResolved
  END IF ! CalcPICTimeStep
  tmpArray(6) = PICValidPlasmaCellSum

  ! Collect sum on MPIRoot
  IF(MPIRoot)THEN
    CALL MPI_REDUCE(MPI_IN_PLACE , tmpArray , 6 , MPI_INTEGER , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
    IF(CalcPointsPerDebyeLength)THEN
       PPDCellResolved(1) = tmpArray(1)
       PPDCellResolved(2) = tmpArray(2)
       PPDCellResolved(3) = tmpArray(3)
       PPDCellResolved(4) = tmpArray(4)
    END IF ! CalcPointsPerDebyeLength
    IF(CalcPICTimeStep)THEN
      PICTimeCellResolved = tmpArray(5)
    END IF ! CalcPICTimeStep
    PICValidPlasmaCellSum = tmpArray(6)
  ELSE
    CALL MPI_REDUCE(tmpArray     , 0        , 6 , MPI_INTEGER , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
  END IF ! MPIRoot
END IF ! CalcPointsPerDebyeLength.OR.CalcPICTimeStep

#endif /*USE_MPI*/

END SUBROUTINE CalculatePartElemData


!===================================================================================================================================
! Count the number of electrons per DG cell and divide it by element-volume
!===================================================================================================================================
SUBROUTINE CalculateElectronIonDensityCell()
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Globals_Vars          ,ONLY: ElementaryCharge
USE MOD_Particle_Analyze_Vars ,ONLY: ElectronDensityCell,IonDensityCell,NeutralDensityCell,ChargeNumberCell
USE MOD_Particle_Vars         ,ONLY: Species,PartSpecies,PDM,PEM,usevMPF
USE MOD_Preproc
USE MOD_DSMC_Vars             ,ONLY: RadialWeighting
USE MOD_Particle_Mesh_Vars    ,ONLY: ElemVolume_Shared
USE MOD_Mesh_Vars             ,ONLY: offSetElem
USE MOD_Mesh_Tools            ,ONLY: GetCNElemID
#if USE_HDG
USE MOD_HDG_Vars              ,ONLY: ElemToBRRegion,UseBRElectronFluid
USE MOD_PIC_Analyze           ,ONLY: CalculateBRElectronsPerCell
#endif /*USE_HDG*/
USE MOD_part_tools            ,ONLY: GetParticleWeight
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iPart,iElem
REAL    :: charge, MPF
#if USE_HDG
INTEGER :: RegionID
#endif /*USE_HDG*/
!===================================================================================================================================
! nullify
ElectronDensityCell=0.
     IonDensityCell=0.
 NeutralDensityCell=0.
   ChargeNumberCell=0.

! loop over all particles and count the number of electrons per cell
! CAUTION: we need the number of all real particle instead of simulated particles
DO iPart=1,PDM%ParticleVecLength
  IF(.NOT.PDM%ParticleInside(iPart)) CYCLE
  IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
    MPF = GetParticleWeight(iPart)
  ELSE
    MPF = GetParticleWeight(iPart) * Species(PartSpecies(iPart))%MacroParticleFactor
  END IF
  ASSOCIATE ( &
    ElemID  => PEM%LocalElemID(iPart)                              )  ! Element ID
    ASSOCIATE ( &
      n_e    => ElectronDensityCell(ElemID),& ! Electron density (cell average)
      n_i    => IonDensityCell(ElemID)     ,& ! Ion density (cell average)
      n_n    => NeutralDensityCell(ElemID) ,& ! Neutral density (cell average)
      Z      => ChargeNumberCell(ElemID)   )  ! Charge number (cell average)
      charge = Species(PartSpecies(iPart))%ChargeIC/ElementaryCharge
      IF(PARTISELECTRON(iPart))THEN ! electrons
        n_e = n_e + MPF
      ELSEIF(ABS(charge).GT.0.0)THEN ! ions (positive or negative)
        n_i = n_i + MPF
        Z   = Z   + charge*MPF
      ELSE ! neutrals
        n_n  = n_n + MPF
      END IF
    END ASSOCIATE
  END ASSOCIATE
END DO ! iPart

#if USE_HDG
IF (UseBRElectronFluid) THEN !check for BR electrons
  DO iElem=1,PP_nElems
    RegionID=ElemToBRRegion(iElem)
    IF (RegionID.GT.0) THEN
      IF (ABS(ElectronDensityCell(iElem)).GT.0.0)THEN
        IPWRITE(UNIT_StdOut,*) "iElem =", iElem
        CALL abort( __STAMP__,'Mixed BR and kinetic electrons are not implemented in CalculateElectronIonDensityCell yet!')
      END IF
      CALL CalculateBRElectronsPerCell(iElem,RegionID,ElectronDensityCell(iElem))
    END IF
  END DO ! iElem=1,PP_nElems
END IF
#endif /*USE_HDG*/

! loop over all elements and divide by volume
DO iElem=1,PP_nElems
  ElectronDensityCell(iElem)=ElectronDensityCell(iElem)/ElemVolume_Shared(GetCNElemID(iElem+offSetElem))
       IonDensityCell(iElem)=IonDensityCell(iElem)     /ElemVolume_Shared(GetCNElemID(iElem+offSetElem))
   NeutralDensityCell(iElem)=NeutralDensityCell(iElem) /ElemVolume_Shared(GetCNElemID(iElem+offSetElem))
END DO ! iElem=1,PP_nElems

END SUBROUTINE CalculateElectronIonDensityCell


!===================================================================================================================================
! Count the number of electrons per DG cell and divide it by element-volume
!===================================================================================================================================
SUBROUTINE CalculateElectronTemperatureCell()
! MODULES                                                                                                                          !
USE MOD_Globals               ,ONLY: PARTISELECTRON
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst,ElectronMass
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars ,ONLY: ElectronTemperatureCell
USE MOD_Particle_Vars         ,ONLY: PDM,PEM,usevMPF,Species,PartSpecies,PartState
USE MOD_DSMC_Vars             ,ONLY: RadialWeighting
#if USE_HDG
USE MOD_HDG_Vars              ,ONLY: ElemToBRRegion,UseBRElectronFluid,RegionElectronRef
USE MOD_Globals_Vars          ,ONLY: ElementaryCharge
#endif /*USE_HDG*/
USE MOD_part_tools            ,ONLY: GetParticleWeight
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iPart,iElem,ElemID,Method
REAL    :: nElectronsPerCell(1:PP_nElems)
REAL    :: PartVandV2(1:PP_nElems,1:6)
REAL    :: Mean_PartV2(1:3)
REAL    :: MeanPartV_2(1:3)
REAL    :: TempDirec(1:3)
REAL    :: WeightingFactor
#if USE_HDG
INTEGER :: RegionID
#endif /*USE_HDG*/
!===================================================================================================================================
#if USE_HDG
IF (UseBRElectronFluid) THEN ! check for BR electrons
  DO iElem=1,PP_nElems
    RegionID=ElemToBRRegion(iElem)
    IF (RegionID.GT.0) THEN
      ElectronTemperatureCell(iElem) = RegionElectronRef(3,RegionID)*ElementaryCharge/BoltzmannConst ! convert eV to K
    END IF
  END DO ! iElem=1,PP_nElems
  RETURN ! Mixed BR and kinetic electrons are not implemented yet!
END IF
#endif /*USE_HDG*/

! nullify
ElectronTemperatureCell=0.
nElectronsPerCell      =0.

! hard-coded
Method=1 ! 0: <E> = (3/2)*<k_B*T> (keeps drift velocity, hence, over-estimates the temperature when drift becomes important)
!        ! 1: remove drift from temperature calculation

PartVandV2 = 0.
! 1.   loop over all particles and sum-up the electron energy per cell and count the number of electrons per cell
DO iPart=1,PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart))THEN
    IF(.NOT.PARTISELECTRON(iPart)) CYCLE  ! ignore anything that is not an electron
    ElemID            = PEM%LocalElemID(iPart)
    IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
      WeightingFactor = GetParticleWeight(iPart)
    ELSE
      WeightingFactor = GetParticleWeight(iPart) * Species(PartSpecies(iPart))%MacroParticleFactor
    END IF
    nElectronsPerCell(ElemID) = nElectronsPerCell(ElemID) + WeightingFactor
    ! Determine velocity or kinetic energy
    SELECT CASE(Method)
    CASE(0) ! 1.0   for distributions where the drift is negligible
      ElectronTemperatureCell(ElemID) = ElectronTemperatureCell(ElemID)+CalcEkinPart(iPart)
    CASE(1) ! 1.1   remove drift from distribution
      PartVandV2(ElemID,1:3) = PartVandV2(ElemID,1:3) + PartState(4:6,iPart)    * WeightingFactor
      PartVandV2(ElemID,4:6) = PartVandV2(ElemID,4:6) + PartState(4:6,iPart)**2 * WeightingFactor
    END SELECT
  END IF ! ParticleInside
END DO ! iPart

! 2.   loop over all elements and divide by electrons per cell to get average kinetic energy
SELECT CASE(Method)
CASE(0) ! 2.0   for distributions where the drift is negligible
  DO iElem=1,PP_nElems
    IF(nElectronsPerCell(iElem).GT.0.) THEN
      ! <E> = (3/2)*<k_B*T>
      ElectronTemperatureCell(iElem)  = 2.*ElectronTemperatureCell(iElem)/(3.*nElectronsPerCell(iElem)*BoltzmannConst)
    END IF
  END DO ! iElem=1,PP_nElems
CASE(1) ! 2.1   remove drift from distribution
  DO iElem=1,PP_nElems
    IF(nElectronsPerCell(iElem).LT.2.) THEN ! only calculate the temperature when more than one electron are present
      ElectronTemperatureCell(iElem) = 0.0
    ELSE
      ! Compute velocity averages
      MeanPartV_2(1:3)  = (PartVandV2(iElem,1:3) / nElectronsPerCell(iElem))**2 ! < |v| >**2
      Mean_PartV2(1:3)  =  PartVandV2(iElem,4:6) / nElectronsPerCell(iElem)     ! < |v|**2 >
      ! Compute temperatures
      TempDirec(1:3) = ElectronMass * (Mean_PartV2(1:3) - MeanPartV_2(1:3)) / BoltzmannConst
      ElectronTemperatureCell(iElem) = (TempDirec(1) + TempDirec(2) + TempDirec(3))/3.0
      IF(ElectronTemperatureCell(iElem).LT.0.0)THEN
        ElectronTemperatureCell(iElem)=0.0
      END IF
    END IF
  END DO
END SELECT

END SUBROUTINE CalculateElectronTemperatureCell


!===================================================================================================================================
! Count the number of electrons per DG cell and divide it by element-volume
!===================================================================================================================================
SUBROUTINE CalculateElectronEnergyCell()
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Globals               ,ONLY: PARTISELECTRON
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst,ElectronMass,Joule2eV
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars ,ONLY: ElectronMinEnergyCell,ElectronMaxEnergyCell,ElectronAverageEnergyCell
USE MOD_Particle_Vars         ,ONLY: PDM,PEM,usevMPF,Species,PartSpecies
USE MOD_DSMC_Vars             ,ONLY: RadialWeighting
#if USE_HDG
USE MOD_HDG_Vars              ,ONLY: ElemToBRRegion,UseBRElectronFluid,RegionElectronRef
USE MOD_Globals_Vars          ,ONLY: ElementaryCharge
#endif /*USE_HDG*/
USE MOD_part_tools            ,ONLY: GetParticleWeight
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iPart,iElem,ElemID
REAL    :: nElectronsPerCell(1:PP_nElems),Ekin
REAL    :: WeightingFactor
#if USE_HDG
INTEGER :: RegionID
#endif /*USE_HDG*/
!===================================================================================================================================
ElectronMinEnergyCell     = HUGE(1.) ! Set zero before output to .h5 if unchanged (check if maximum is <= 0.)
ElectronMaxEnergyCell     = 0.
ElectronAverageEnergyCell = 0.
#if USE_HDG
IF (UseBRElectronFluid) THEN ! check for BR electrons
  DO iElem=1,PP_nElems
    RegionID=ElemToBRRegion(iElem)
    IF (RegionID.GT.0) THEN
      ! Assume <E> = (3/2)*<k_B*T>
      ! with T = RegionElectronRef(3,RegionID)*ElementaryCharge/BoltzmannConst ! convert eV to K
      ! T*Joule2eV gives energy in [eV]
      ElectronMinEnergyCell     = 1.5*RegionElectronRef(3,RegionID) ! convert to [eV] -> factors cancel out
      ElectronMaxEnergyCell     = ElectronMinEnergyCell
      ElectronAverageEnergyCell = ElectronMinEnergyCell
    END IF
  END DO ! iElem=1,PP_nElems
  RETURN ! Mixed BR and kinetic electrons are not implemented yet!
ELSE
#endif /*USE_HDG*/
  nElectronsPerCell = 0.
  ! 1. loop over all particles and sum-up the electron energy per cell and count the number of electrons per cell
  DO iPart=1,PDM%ParticleVecLength
    IF(PDM%ParticleInside(iPart))THEN
      IF(.NOT.PARTISELECTRON(iPart)) CYCLE  ! ignore anything that is not an electron
      ElemID            = PEM%LocalElemID(iPart)
      IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
        WeightingFactor = GetParticleWeight(iPart)
      ELSE
        WeightingFactor = GetParticleWeight(iPart) * Species(PartSpecies(iPart))%MacroParticleFactor
      END IF
      nElectronsPerCell(ElemID) = nElectronsPerCell(ElemID) + WeightingFactor
      ! Determine kinetic energy
      Ekin = CalcEkinPart(iPart)
      ElectronAverageEnergyCell(ElemID) = ElectronAverageEnergyCell(ElemID) + Ekin
      ElectronMinEnergyCell(ElemID)     = MIN(ElectronMinEnergyCell(ElemID),Ekin/WeightingFactor) ! Consider energy of single e-
      ElectronMaxEnergyCell(ElemID)     = MAX(ElectronMaxEnergyCell(ElemID),Ekin/WeightingFactor) ! Consider energy of single e-
    END IF ! ParticleInside
  END DO ! iPart

  ! 2. Loop over all elements and calculate the average, also check if minimum is HUGE(1.)
  DO iElem=1,PP_nElems
    IF(nElectronsPerCell(iElem).GT.0.)THEN
      ElectronAverageEnergyCell(iElem) = ElectronAverageEnergyCell(iElem)/nElectronsPerCell(iElem)*Joule2eV
      ElectronMinEnergyCell(iElem)     = ElectronMinEnergyCell(iElem)*Joule2eV
      ElectronMaxEnergyCell(iElem)     = ElectronMaxEnergyCell(iElem)*Joule2eV
    ELSE
      ElectronMinEnergyCell(iElem) = 0. ! Set from HUGE(1.) to zero for output to .h5
    END IF ! nElectronsPerCell(iElem).GT.0.
  END DO ! iElem=1,PP_nElems

#if USE_HDG
END IF
#endif /*USE_HDG*/

END SUBROUTINE CalculateElectronEnergyCell


SUBROUTINE CalcShapeEfficiencyR()
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Particle_Analyze_Vars ,ONLY: CalcShapeEfficiencyMethod, ShapeEfficiencyNumber
USE MOD_Mesh_Vars             ,ONLY: nElems, Elem_xGP
USE MOD_Particle_Mesh_Vars    ,ONLY: GEO
USE MOD_PICDepo_Vars
USE MOD_Particle_Vars
USE MOD_Preproc
#if USE_MPI
USE MOD_Globals
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                     :: NbrOfComps, NbrWithinRadius, NbrOfElems, NbrOfElemsWithinRadius
REAL                     :: RandVal1
LOGICAL                  :: chargedone(1:nElems), WITHIN
INTEGER                  :: kmin, kmax, lmin, lmax, mmin, mmax
INTEGER                  :: kk, ll, mm, ppp,m,l,k, i
INTEGER                  :: ElemID
REAL                     :: radius, deltax, deltay, deltaz
!===================================================================================================================================

NbrOfComps = 0.
NbrOfElems = 0.
NbrWithinRadius = 0.
NbrOfElemsWithinRadius = 0.
SELECT CASE(CalcShapeEfficiencyMethod)
CASE('AllParts')
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      chargedone(:) = .FALSE.
      !-- determine which background mesh cells (and interpolation points within) need to be considered
      kmax = INT((PartState(1,i)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
      kmax = MIN(kmax,GEO%FIBGMimax)
      kmin = INT((PartState(1,i)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
      kmin = MAX(kmin,GEO%FIBGMimin)
      lmax = INT((PartState(2,i)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
      lmax = MIN(lmax,GEO%FIBGMjmax)
      lmin = INT((PartState(2,i)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
      lmin = MAX(lmin,GEO%FIBGMjmin)
      mmax = INT((PartState(3,i)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
      mmax = MIN(mmax,GEO%FIBGMkmax)
      mmin = INT((PartState(3,i)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
      mmin = MAX(mmin,GEO%FIBGMkmin)
      !-- go through all these cells
      DO kk = kmin,kmax
        DO ll = lmin, lmax
          DO mm = mmin, mmax
            !--- go through all mapped elements not done yet
            DO ppp = 1,GEO%FIBGM(kk,ll,mm)%nElem
              WITHIN=.FALSE.
              ElemID = GEO%FIBGM(kk,ll,mm)%Element(ppp)
              IF (.NOT.chargedone(ElemID)) THEN
                NbrOfElems = NbrOfElems + 1.
                !--- go through all gauss points
                DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
                  NbrOfComps = NbrOfComps + 1.
                  !-- calculate distance between gauss and particle
                  deltax = PartState(1,i) - Elem_xGP(1,k,l,m,ElemID)
                  deltay = PartState(2,i) - Elem_xGP(2,k,l,m,ElemID)
                  deltaz = PartState(3,i) - Elem_xGP(3,k,l,m,ElemID)
                  radius = deltax * deltax + deltay * deltay + deltaz * deltaz
                  IF (radius .LT. r2_sf) THEN
                    WITHIN=.TRUE.
                    NbrWithinRadius = NbrWithinRadius + 1.
                  END IF
                END DO; END DO; END DO
                chargedone(ElemID) = .TRUE.
              END IF
              IF(WITHIN) NbrOfElemsWithinRadius = NbrOfElemsWithinRadius + 1.
            END DO ! ppp
          END DO ! mm
        END DO ! ll
      END DO ! kk
    END IF ! inside
  END DO ! i
IF(NbrOfComps.GT.0.0)THEN
#if USE_MPI
  WRITE(*,*) 'ShapeEfficiency (Proc,%,%Elems)',myRank,100*NbrWithinRadius/NbrOfComps,100*NbrOfElemsWithinRadius/NbrOfElems
  WRITE(*,*) 'ShapeEfficiency (Elems) for Proc',myRank,'is',100*NbrOfElemsWithinRadius/NbrOfElems,'%'
#else
  WRITE(*,*) 'ShapeEfficiency (%,%Elems)',100*NbrWithinRadius/NbrOfComps, 100*NbrOfElemsWithinRadius/NbrOfElems
  WRITE(*,*) 'ShapeEfficiency (Elems) is',100*NbrOfElemsWithinRadius/NbrOfElems,'%'
#endif
END IF
CASE('SomeParts')
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      CALL RANDOM_NUMBER(RandVal1)
      IF(RandVal1.LT.REAL(ShapeEfficiencyNumber)/100)THEN
        chargedone(:) = .FALSE.
        !-- determine which background mesh cells (and interpolation points within) need to be considered
        kmax = INT((PartState(1,i)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
        kmax = MIN(kmax,GEO%FIBGMimax)
        kmin = INT((PartState(1,i)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
        kmin = MAX(kmin,GEO%FIBGMimin)
        lmax = INT((PartState(2,i)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
        lmax = MIN(lmax,GEO%FIBGMjmax)
        lmin = INT((PartState(2,i)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
        lmin = MAX(lmin,GEO%FIBGMjmin)
        mmax = INT((PartState(3,i)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
        mmax = MIN(mmax,GEO%FIBGMkmax)
        mmin = INT((PartState(3,i)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
        mmin = MAX(mmin,GEO%FIBGMkmin)
        !-- go through all these cells
        DO kk = kmin,kmax
          DO ll = lmin, lmax
            DO mm = mmin, mmax
              !--- go through all mapped elements not done yet
              DO ppp = 1,GEO%FIBGM(kk,ll,mm)%nElem
              WITHIN=.FALSE.
                ElemID = GEO%FIBGM(kk,ll,mm)%Element(ppp)
                IF (.NOT.chargedone(ElemID)) THEN
                  NbrOfElems = NbrOfElems + 1
                  !--- go through all gauss points
                  DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
                    NbrOfComps = NbrOfComps + 1
                    !-- calculate distance between gauss and particle
                    deltax = PartState(1,i) - Elem_xGP(1,k,l,m,ElemID)
                    deltay = PartState(2,i) - Elem_xGP(2,k,l,m,ElemID)
                    deltaz = PartState(3,i) - Elem_xGP(3,k,l,m,ElemID)
                    radius = deltax * deltax + deltay * deltay + deltaz * deltaz
                    IF (radius .LT. r2_sf) THEN
                      NbrWithinRadius = NbrWithinRadius + 1
                      WITHIN=.TRUE.
                    END IF
                    END DO; END DO; END DO
                    chargedone(ElemID) = .TRUE.
                  END IF
                  IF(WITHIN) NbrOfElemsWithinRadius = NbrOfElemsWithinRadius + 1
                END DO ! ppp
              END DO ! mm
            END DO ! ll
          END DO ! kk
        END IF  ! RandVal
      END IF ! inside
  END DO ! i
IF(NbrOfComps.GT.0)THEN
#if USE_MPI
  WRITE(*,*) 'ShapeEfficiency (Proc,%,%Elems)',myRank,100*NbrWithinRadius/NbrOfComps,100*NbrOfElemsWithinRadius/NbrOfElems
  WRITE(*,*) 'ShapeEfficiency (Elems) for Proc',myRank,'is',100*NbrOfElemsWithinRadius/NbrOfElems,'%'
#else
  WRITE(*,*) 'ShapeEfficiency (%,%Elems)',100*NbrWithinRadius/NbrOfComps,100*NbrOfElemsWithinRadius/NbrOfElems
  WRITE(*,*) 'ShapeEfficiency (Elems) is',100*NbrOfElemsWithinRadius/NbrOfElems,'%'
#endif
END IF
END SELECT
END SUBROUTINE CalcShapeEfficiencyR


PPURE SUBROUTINE CalcKineticEnergy(Ekin)
!===================================================================================================================================
! compute the kinetic energy of particles
! for velocity <1e6 non-relativistic formula is used, for larger velocities the relativistic kinetic energy is computed
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Globals_Vars          ,ONLY: c2, c2_inv, RelativisticLimit
USE MOD_Particle_Vars         ,ONLY: PartState, PartSpecies, Species, PDM, PEM
USE MOD_PARTICLE_Vars         ,ONLY: usevMPF
USE MOD_Particle_Analyze_Vars ,ONLY: nSpecAnalyze
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_DSMC_Vars             ,ONLY: RadialWeighting
#if !(USE_HDG)
USE MOD_PML_Vars              ,ONLY: DoPML,isPMLElem
#endif /*USE_HDG*/
USE MOD_Dielectric_Vars       ,ONLY: DoDielectric,isDielectricElem,DielectricNoParticles
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Ekin(nSpecAnalyze)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,ElemID
REAL(KIND=8)                    :: partV2, GammaFac
REAL                            :: Ekin_loc
!===================================================================================================================================
Ekin    = 0.!d0
IF (nSpecAnalyze.GT.1) THEN
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      ElemID = PEM%LocalElemID(i)
#if !(USE_HDG)
      IF(DoPML)THEN
        IF(isPMLElem(ElemID)) CYCLE
      ENDIF
#endif /*USE_HDG*/
      IF(DoDielectric)THEN
        IF(DielectricNoParticles)THEN
          IF(isDielectricElem(ElemID)) CYCLE
        END IF ! DielectricNoParticles
      ENDIF
      partV2 = DOTPRODUCT(PartState(4:6,i))
      IF ( partV2 .LT. RelativisticLimit) THEN  ! |v| < 1000000 when speed of light is 299792458
        Ekin_loc = 0.5 * Species(PartSpecies(i))%MassIC * partV2
        IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
          ! %MacroParticleFactor is included in the case of vMPF (also in combination with variable time step)
          Ekin(nSpecAnalyze)   = Ekin(nSpecAnalyze)   + Ekin_loc * GetParticleWeight(i)
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * GetParticleWeight(i)
        ELSE
          ! Case for variable time step without vMPF (no species weighting factor applied in GetParticleWeight)
          Ekin(nSpecAnalyze)   = Ekin(nSpecAnalyze)   + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
        END IF != usevMPF
      ELSE ! partV2 > RelativisticLimit
        GammaFac = partV2*c2_inv
        GammaFac = 1./SQRT(1.-GammaFac)
        Ekin_loc = (GammaFac-1.) * Species(PartSpecies(i))%MassIC * c2
        IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
          ! %MacroParticleFactor is included in the case of vMPF (also in combination with variable time step)
          Ekin(nSpecAnalyze)   = Ekin(nSpecAnalyze)   + Ekin_loc * GetParticleWeight(i)
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * GetParticleWeight(i)
        ELSE
          ! Case for variable time step without vMPF (no species weighting factor applied in GetParticleWeight)
          Ekin(nSpecAnalyze)   = Ekin(nSpecAnalyze)   + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
        END IF != usevMPF
      END IF ! partV2
    END IF ! (PDM%ParticleInside(i))
  END DO ! i=1,PDM%ParticleVecLength
ELSE ! nSpecAnalyze = 1 : only 1 species
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      ElemID = PEM%LocalElemID(i)
#if !(USE_HDG)
      IF(DoPML)THEN
        IF(isPMLElem(ElemID)) CYCLE
      ENDIF
#endif /*USE_HDG*/
      IF(DoDielectric)THEN
        IF(DielectricNoParticles)THEN
          IF(isDielectricElem(ElemID)) CYCLE
        END IF ! DielectricNoParticles
      ENDIF
      partV2 = DOTPRODUCT(PartState(4:6,i))
      IF ( partV2 .LT. RelativisticLimit) THEN  ! |v| < 1000000 when speed of light is 299792458
        Ekin_loc = 0.5 *  Species(PartSpecies(i))%MassIC * partV2
        IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * GetParticleWeight(i)
        ELSE
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
        END IF ! usevMPF
      ELSE ! partV2 > RelativisticLimit
        GammaFac = partV2*c2_inv
        GammaFac = 1./SQRT(1.-GammaFac)
        Ekin_loc = (GammaFac-1.) * Species(PartSpecies(i))%MassIC * c2
        IF(usevMPF.OR.RadialWeighting%DoRadialWeighting)THEN
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * GetParticleWeight(i)
        ELSE
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
        END IF ! useuvMPF
      END IF ! partV2 .LT. RelativisticLimit
    END IF ! particle inside
  END DO ! particleveclength
END IF

END SUBROUTINE CalcKineticEnergy


PPURE SUBROUTINE CalcKineticEnergyAndMaximum(Ekin,EkinMax)
!===================================================================================================================================
! compute the kinetic energy of particles
! for velocity <1e6 non-relativistic formula is used, for larger velocities the relativistic kinetic energy is computed
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Globals_Vars          ,ONLY: c2, c2_inv, RelativisticLimit
USE MOD_Particle_Vars         ,ONLY: PartState, PartSpecies, Species, PDM, nSpecies, PEM
USE MOD_PARTICLE_Vars         ,ONLY: usevMPF
USE MOD_Particle_Analyze_Vars ,ONLY: nSpecAnalyze,LaserInteractionEkinMaxRadius,LaserInteractionEkinMaxZPosMin
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_DSMC_Vars             ,ONLY: RadialWeighting
#if !(USE_HDG)
USE MOD_PML_Vars              ,ONLY: DoPML,isPMLElem
#endif /*USE_HDG*/
USE MOD_Dielectric_Vars       ,ONLY: DoDielectric,isDielectricElem,DielectricNoParticles
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Ekin(nSpecAnalyze)
REAL,INTENT(OUT)                :: EkinMax(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,ElemID
REAL(KIND=8)                    :: partV2, GammaFac
REAL                            :: Ekin_loc
!===================================================================================================================================
! default values
Ekin    =  0.
EkinMax = -1.
! set boundaries in order to exclude particles near the boundary (nonphysical velocities)
IF (nSpecAnalyze.GT.1) THEN
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      ElemID = PEM%LocalElemID(i)
#if !(USE_HDG)
      IF(DoPML)THEN
        IF(isPMLElem(ElemID)) CYCLE
      ENDIF
#endif /*USE_HDG*/
      IF(DoDielectric)THEN
        IF(DielectricNoParticles)THEN
          IF(isDielectricElem(ElemID)) CYCLE
        END IF ! DielectricNoParticles
      ENDIF
      partV2 = DOTPRODUCT(PartState(4:6,i))
      IF ( partV2 .LT. RelativisticLimit) THEN  ! |v| < 1000000 when speed of light is 299792458
        Ekin_loc = 0.5 * Species(PartSpecies(i))%MassIC * partV2
        IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
          ! %MacroParticleFactor is included in the case of RadialWeighting (also in combination with variable time step)
          Ekin(nSpecAnalyze)   = Ekin(nSpecAnalyze)   + Ekin_loc * GetParticleWeight(i)
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * GetParticleWeight(i)
        ELSE
          ! Case for variable time step without radial weighting (no regular weighting factor applied in GetParticleWeight)
          Ekin(nSpecAnalyze)   = Ekin(nSpecAnalyze)   + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
        END IF != usevMPF
      ELSE ! partV2 > RelativisticLimit
        GammaFac = partV2*c2_inv
        GammaFac = 1./SQRT(1.-GammaFac)
        Ekin_loc = (GammaFac-1.) * Species(PartSpecies(i))%MassIC * c2
        IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
          Ekin(nSpecAnalyze)   = Ekin(nSpecAnalyze)   + Ekin_loc * GetParticleWeight(i)
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * GetParticleWeight(i)
        ELSE
          Ekin(nSpecAnalyze)   = Ekin(nSpecAnalyze)   + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
        END IF !=usevMPF
      END IF ! partV2
      ! Determine energy of the most energetic particle in [eV]
      IF((SQRT(PartState(1,i)**2 + PartState(2,i)**2).LE.LaserInteractionEkinMaxRadius).OR.&
                                      (PartState(3,i).GE.LaserInteractionEkinMaxZPosMin))THEN
        EkinMax(PartSpecies(i)) = MAX(EkinMax(PartSpecies(i)),Ekin_loc*6.241509e18) ! 6.241509e18 is [J] -> [eV]
      END IF
    END IF ! (PDM%ParticleInside(i))
  END DO ! i=1,PDM%ParticleVecLength
ELSE ! nSpecAnalyze = 1 : only 1 species
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      ElemID = PEM%LocalElemID(i)
#if !(USE_HDG)
      IF(DoPML)THEN
        IF(isPMLElem(ElemID)) CYCLE
      ENDIF
#endif /*USE_HDG*/
      IF(DoDielectric)THEN
        IF(DielectricNoParticles)THEN
          IF(isDielectricElem(ElemID)) CYCLE
        END IF ! DielectricNoParticles
      ENDIF
      partV2 = DOTPRODUCT(PartState(4:6,i))
      IF ( partV2 .LT. RelativisticLimit) THEN ! |v| < 1000000 when speed of light is 299792458
        Ekin_loc = 0.5 *  Species(PartSpecies(i))%MassIC * partV2
        IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * GetParticleWeight(i)
        ELSE
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
        END IF ! usevMPF
      ELSE ! partV2 > RelativisticLimit
        GammaFac = partV2*c2_inv
        GammaFac = 1./SQRT(1.-GammaFac)
        Ekin_loc = (GammaFac-1.) * Species(PartSpecies(i))%MassIC * c2
        IF(usevMPF.OR.RadialWeighting%DoRadialWeighting)THEN
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * GetParticleWeight(i)
        ELSE
          Ekin(PartSpecies(i)) = Ekin(PartSpecies(i)) + Ekin_loc * Species(PartSpecies(i))%MacroParticleFactor*GetParticleWeight(i)
        END IF ! useuvMPF

      END IF ! par2
      ! Determine energy of the most energetic particle in [eV]
      IF((SQRT(PartState(1,i)**2 + PartState(2,i)**2).LE.LaserInteractionEkinMaxRadius).OR.&
                                      (PartState(3,i).GE.LaserInteractionEkinMaxZPosMin))THEN
        EkinMax(PartSpecies(i)) = MAX(EkinMax(PartSpecies(i)),Ekin_loc*6.241509e18) ! 6.241509e18 is [J] -> [eV]
      END IF
    END IF ! particle inside
  END DO ! particleveclength
END IF

END SUBROUTINE CalcKineticEnergyAndMaximum


!===================================================================================================================================
!> Computes the number density per species using the total mesh volume and if necessary particle weights
!> Background gas density is saved as given in the input
!===================================================================================================================================
PPURE SUBROUTINE CalcNumberDensity(NumSpec,NumDens)
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_DSMC_Vars             ,ONLY: BGGas, RadialWeighting
USE MOD_Particle_Analyze_Vars ,ONLY: nSpecAnalyze
USE MOD_Particle_Vars         ,ONLY: Species,nSpecies,usevMPF
USE MOD_Particle_Mesh_Vars    ,ONLY: MeshVolume
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
REAL,INTENT(IN)                   :: NumSpec(nSpecAnalyze)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: NumDens(nSpecAnalyze)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: iSpec,bgSpec
!===================================================================================================================================

! Only root does calculation
IF(.NOT.MPIRoot) RETURN

IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
  NumDens(1:nSpecies) = NumSpec(1:nSpecies) / MeshVolume
ELSE
  NumDens(1:nSpecies) = NumSpec(1:nSpecies) * Species(1:nSpecies)%MacroParticleFactor / MeshVolume
END IF

IF(BGGas%NumberOfSpecies.GT.0) THEN
  DO iSpec = 1, nSpecies
    IF(BGGas%BackgroundSpecies(iSpec)) THEN
      bgSpec = BGGas%MapSpecToBGSpec(iSpec)
      IF(BGGas%UseDistribution) THEN
        ! Use pre-calculate value
        NumDens(iSpec) = BGGas%DistributionNumDens(bgSpec)
      ELSE
        NumDens(iSpec) = BGGas%NumberDensity(bgSpec)
      END IF
    END IF
  END DO
END IF

IF(nSpecAnalyze.GT.1) NumDens(nSpecAnalyze) = SUM(NumDens(1:nSpecies))

END SUBROUTINE CalcNumberDensity


!===================================================================================================================================
!> Computes the number density per species when using BGGas distribution
!===================================================================================================================================
SUBROUTINE CalcNumberDensityBGGasDistri()
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_DSMC_Vars             ,ONLY: BGGas
USE MOD_Particle_Vars         ,ONLY: nSpecies
USE MOD_Particle_Mesh_Vars    ,ONLY: MeshVolume,ElemVolume_Shared
USE MOD_Mesh_Tools            ,ONLY: GetCNElemID
USE MOD_Mesh_Vars             ,ONLY: nElems,offSetElem
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: iSpec,bgSpec,iElem
REAL                              :: DistriNumDens(1:BGGas%NumberOfSpecies)
!===================================================================================================================================
! Initialize
DistriNumDens = 0.
ALLOCATE(BGGas%DistributionNumDens(1:BGGas%NumberOfSpecies))
BGGas%DistributionNumDens = 0.

! Loop over all species and elements and add up the mass within each element
DO iSpec = 1, nSpecies
  IF(BGGas%BackgroundSpecies(iSpec)) THEN
    bgSpec = BGGas%MapSpecToBGSpec(iSpec)
    DistriNumDens(bgSpec) = 0.
    DO iElem = 1, nElems
      ! Calculate mass per element (divide by total mesh volume later on)
      DistriNumDens(bgSpec) = DistriNumDens(bgSpec) &
                            + BGGas%Distribution(bgSpec,7,iElem) * ElemVolume_Shared(GetCNElemID(iElem+offSetElem))
    END DO ! iElem = 1, nElems
  END IF
END DO

! Communicate
#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE , DistriNumDens, BGGas%NumberOfSpecies, MPI_DOUBLE_PRECISION, MPI_SUM, 0,MPI_COMM_PICLAS ,IERROR)
ELSE
  CALL MPI_REDUCE(DistriNumDens, 0.           , BGGas%NumberOfSpecies, MPI_DOUBLE_PRECISION, MPI_SUM, 0,MPI_COMM_PICLAS ,IERROR)
END IF
#endif /*USE_MPI*/
BGGas%DistributionNumDens = DistriNumDens

! Average over total mesh
BGGas%DistributionNumDens = BGGas%DistributionNumDens / MeshVolume

END SUBROUTINE CalcNumberDensityBGGasDistri


SUBROUTINE CalcSurfaceFluxInfo()
!===================================================================================================================================
!> Output for surface flux BCs: calculate the mass flow rate, current or the pressure in adjacent cells per species & surface flux
!> 1) Calculate processor-local values
!> 2) MPI communication
!> 3) Determine final output values
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_TimeDisc_Vars           ,ONLY: dt, iter
USE MOD_Particle_Analyze_Vars   ,ONLY: FlowRateSurfFlux,PressureAdaptiveBC,PartAnalyzeStep
USE MOD_DSMC_Vars               ,ONLY: RadialWeighting
USE MOD_Particle_Vars           ,ONLY: Species,nSpecies,usevMPF,VarTimeStep
USE MOD_Particle_Surfaces_Vars  ,ONLY: BCdata_auxSF, SurfFluxSideSize, SurfMeshSubSideData
USE MOD_Particle_Sampling_Vars  ,ONLY: UseAdaptiveBC, AdaptBCMacroVal, AdaptBCMapElemToSample, AdaptBCAreaSurfaceFlux
USE MOD_Particle_Sampling_Vars  ,ONLY: AdaptBCAverageValBC, AdaptBCAverageMacroVal
USE MOD_Mesh_Vars               ,ONLY: SideToElem
#if USE_MPI
USE MOD_Particle_Analyze_Vars   ,ONLY: nSpecAnalyze
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iSpec, iSF, ElemID, SampleElemID, SurfSideID, iSide, iSample, jSample, currentBC
REAL                :: MacroParticleFactor, dtVar
#if USE_MPI
INTEGER             :: MaxSurfaceFluxBCs
#endif /*USE_MPI*/
!===================================================================================================================================

IF(iter.EQ.0) RETURN

IF(UseAdaptiveBC) PressureAdaptiveBC = 0.
! 1) Calculate the processor-local mass flow rate and sum-up the area weighted pressure
DO iSpec = 1, nSpecies
  ! If usevMPF or DoRadialWeighting then the MacroParticleFactor is already included in the GetParticleWeight
  IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
    MacroParticleFactor = 1.
  ELSE
    MacroParticleFactor = Species(iSpec)%MacroParticleFactor
  END IF
  IF(VarTimeStep%UseSpeciesSpecific) THEN
    dtVar = dt * Species(iSpec)%TimeStepFactor
  ELSE
    dtVar = dt
  END IF
  DO iSF = 1, Species(iSpec)%nSurfacefluxBCs
    ! SampledMassFlow contains the weighted particle number balance (in - out)
    FlowRateSurfFlux(iSpec,iSF) = Species(iSpec)%Surfaceflux(iSF)%SampledMassflow * MacroParticleFactor / dtVar
    IF(Species(iSpec)%Surfaceflux(iSF)%UseEmissionCurrent) THEN
      FlowRateSurfFlux(iSpec,iSF) = FlowRateSurfFlux(iSpec,iSF) * ABS(Species(iSpec)%ChargeIC)
    ELSE
      FlowRateSurfFlux(iSpec,iSF) = FlowRateSurfFlux(iSpec,iSF) * Species(iSpec)%MassIC
    END IF
    ! Reset the mass flow rate
    Species(iSpec)%Surfaceflux(iSF)%SampledMassflow = 0.
    ! Calculate the average pressure
    IF(UseAdaptiveBC) THEN
      currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
      ! Average of the BC for the output
      IF(.NOT.AdaptBCAverageValBC) THEN
        ! Skip processors without a surface flux
        IF (BCdata_auxSF(currentBC)%SideNumber.EQ.0) CYCLE
        ! Average the value of the BC
        ! Loop over sides on the surface flux
        DO iSide=1,BCdata_auxSF(currentBC)%SideNumber
          SurfSideID = BCdata_auxSF(currentBC)%SideList(iSide)
          ElemID = SideToElem(S2E_ELEM_ID,SurfSideID)
          SampleElemID = AdaptBCMapElemToSample(ElemID)
          IF(SampleElemID.GT.0) THEN
            ! Sum up the area weighted pressure in each adjacent element
            DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
              PressureAdaptiveBC(iSpec,iSF) = PressureAdaptiveBC(iSpec,iSF) &
                + AdaptBCMacroVal(6,SampleElemID,iSpec) * SurfMeshSubSideData(iSample,jSample,SurfSideID)%area
            END DO; END DO
          END IF
        END DO
      END IF    ! AdaptBCAverageValBC
    END IF
  END DO
END DO

! 2) Get the sum of the mass flow rate and the sum of the area-weighted area pressures (only in case of a cell-local pressure
!    distribution at the BC, ie. NOT if AdaptBCAverageValBC = T)
#if USE_MPI
MaxSurfaceFluxBCs = MAXVAL(Species(:)%nSurfacefluxBCs)
IF (MPIRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,FlowRateSurfFlux(1:nSpecAnalyze,1:MaxSurfaceFluxBCs),nSpecAnalyze*MaxSurfaceFluxBCs,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
  IF(UseAdaptiveBC.AND.(.NOT.AdaptBCAverageValBC)) THEN
    CALL MPI_REDUCE(MPI_IN_PLACE,PressureAdaptiveBC(1:nSpecAnalyze,1:MaxSurfaceFluxBCs),nSpecAnalyze*MaxSurfaceFluxBCs,&
                    MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
  END IF
ELSE ! no Root
  CALL MPI_REDUCE(FlowRateSurfFlux,FlowRateSurfFlux,nSpecAnalyze*MaxSurfaceFluxBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
  IF(UseAdaptiveBC.AND.(.NOT.AdaptBCAverageValBC)) CALL MPI_REDUCE(PressureAdaptiveBC,PressureAdaptiveBC,nSpecAnalyze*MaxSurfaceFluxBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
END IF
#endif /*USE_MPI*/

! 3) Consider Part-AnalyzeStep for FlowRateSurfFlux and determine the average pressure (value does not depend on the Part-AnalyzeStep)
IF (MPIRoot) THEN
  IF(PartAnalyzeStep.GT.1)THEN
    IF(PartAnalyzeStep.EQ.HUGE(PartAnalyzeStep))THEN
      FlowRateSurfFlux = FlowRateSurfFlux / iter
    ELSE
      FlowRateSurfFlux = FlowRateSurfFlux / MIN(PartAnalyzeStep,iter)
    END IF
  END IF
  IF(UseAdaptiveBC) THEN
    IF(AdaptBCAverageValBC) THEN
      PressureAdaptiveBC(:,:) = AdaptBCAverageMacroVal(3,:,:)
    ELSE
      DO iSpec = 1, nSpecies
        DO iSF = 1, Species(iSpec)%nSurfacefluxBCs
          IF(Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
            ! Use the area sum of the elements included in the sampling, instead of the ideal circular area
            PressureAdaptiveBC(iSpec,iSF) = PressureAdaptiveBC(iSpec,iSF) / AdaptBCAreaSurfaceFlux(iSpec,iSF)
          ELSE
            PressureAdaptiveBC(iSpec,iSF) = PressureAdaptiveBC(iSpec,iSF) / Species(iSpec)%Surfaceflux(iSF)%totalAreaSF
          END IF
        END DO
      END DO
    END IF
  END IF
END IF

END SUBROUTINE CalcSurfaceFluxInfo


SUBROUTINE CalcMixtureTemp(NumSpec,Temp,IntTemp,IntEn,TempTotal,Xi_Vib,Xi_Elec)
!===================================================================================================================================
!> Computes the species-specific and mixture temperature (MPI communication is in the respective subroutines)
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_PARTICLE_Vars             ,ONLY: nSpecies
USE MOD_Particle_Analyze_Vars     ,ONLY: nSpecAnalyze
USE MOD_DSMC_Vars                 ,ONLY: SpecDSMC, CollisMode, DSMC
USE MOD_part_tools                ,ONLY: CalcXiElec
USE MOD_DSMC_Relaxation           ,ONLY: CalcXiVib
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
REAL, INTENT(IN)                  :: NumSpec(nSpecAnalyze)    ! number of real particles (already GLOBAL number)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL, INTENT(IN)                  :: Temp(nSpecAnalyze)
REAL, INTENT(OUT)                 :: IntEn(nSpecAnalyze,3)
REAL, INTENT(OUT)                 :: IntTemp(nSpecies,3)
REAL, INTENT(OUT)                 :: TempTotal(nSpecAnalyze)
REAL, INTENT(OUT)                 :: Xi_Vib(nSpecies), Xi_Elec(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: iSpec
REAL                              :: TempTotalDOF, XiTotal
!===================================================================================================================================

IF (CollisMode.GT.1) THEN
  CALL CalcIntTempsAndEn(NumSpec,IntTemp,IntEn)
  IF(MPIRoot)THEN
    TempTotal = 0.0
    Xi_Vib = 0.0
    Xi_Elec = 0.0
    DO iSpec = 1, nSpecies
      TempTotalDOF = 3.*Temp(iSpec)
      XiTotal = 3.
      ! If the species is molecular, add the vibrational energy to the temperature calculation
      IF(((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)).AND.(NumSpec(iSpec).GT.0)) THEN
        XiTotal = XiTotal + SpecDSMC(iSpec)%Xi_Rot
        TempTotalDOF = TempTotalDOF + SpecDSMC(iSpec)%Xi_Rot*IntTemp(iSpec,2)
        IF(IntTemp(iSpec,1).GT.0) THEN
          CALL CalcXiVib(IntTemp(iSpec,1), iSpec, XiVibTotal=Xi_Vib(iSpec))
          XiTotal = XiTotal + Xi_Vib(iSpec)
          TempTotalDOF = TempTotalDOF + Xi_Vib(iSpec)*IntTemp(iSpec,1)
        END IF
      END IF
      ! If electronic energy is greater zero, added it to the temperature calculation
      IF(IntTemp(iSpec,3).GT.0..AND.(DSMC%ElectronicModel.NE.3)) THEN
        Xi_Elec(iSpec) = CalcXiElec(IntTemp(iSpec,3), iSpec)
        XiTotal = XiTotal + Xi_Elec(iSpec)
        TempTotalDOF = TempTotalDOF + Xi_Elec(iSpec)*IntTemp(iSpec,3)
      END IF
      ! Calculate the species-specific total temperature
      TempTotal(iSpec) = TempTotalDOF / XiTotal
      ! Calculate the total temperature of the mixture (weighted with the particle number)
      IF(nSpecAnalyze.GT.1)THEN
        TempTotal(nSpecAnalyze) = TempTotal(nSpecAnalyze) + TempTotal(iSpec)*NumSpec(iSpec)
      END IF
    END DO
    IF(nSpecAnalyze.GT.1)THEN
      IF(NumSpec(iSpec).NE.0) THEN
        TempTotal(nSpecAnalyze) = TempTotal(nSpecAnalyze) / NumSpec(nSpecAnalyze)
      ELSE
        TempTotal(nSpecAnalyze)= 0.
      END IF
    END IF
  END IF
ELSE
  IF(MPIRoot)THEN
    TempTotal = Temp
    IntTemp   = 0.0
    IntEn     = 0.0
    Xi_Vib    = 0.0
    Xi_Elec   = 0.0
  END IF
END IF

END SUBROUTINE CalcMixtureTemp


SUBROUTINE CalcIntTempsAndEn(NumSpec,IntTemp,IntEn)
!===================================================================================================================================
! Calculation of internal Temps (TVib, TRot, Telec) and gives back the global values
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_Particle_Vars         ,ONLY: PartSpecies, Species, PDM, nSpecies, usevMPF
USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, SpecDSMC, DSMC
USE MOD_Particle_Analyze_Vars ,ONLY: nSpecAnalyze
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_DSMC_Vars             ,ONLY: RadialWeighting
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)               :: NumSpec(nSpecAnalyze)    ! number of real particles (already GLOBAL number)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)               :: IntTemp(nSpecies,3) , IntEn(nSpecAnalyze,3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iPart, iSpec
REAL                           :: EVib(nSpecies), ERot(nSpecies), Eelec(nSpecies), tempVib, NumSpecTemp
#if USE_MPI
REAL                           :: RD(nSpecies)
#endif /*USE_MPI*/
!===================================================================================================================================
EVib    = 0.
ERot    = 0.
Eelec   = 0.
! set electronic state to zero
IntEn(:,:) = 0.
IntTemp(:,:) = 0.

! Sum up internal energies
DO iPart=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    iSpec = PartSpecies(iPart)
    EVib(iSpec) = EVib(iSpec) + PartStateIntEn(1,iPart) * GetParticleWeight(iPart)
    ERot(iSpec) = ERot(iSpec) + PartStateIntEn(2,iPart) * GetParticleWeight(iPart)
    IF (DSMC%ElectronicModel.GT.0) THEN
      IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
        Eelec(iSpec) = Eelec(iSpec) + PartStateIntEn(3,iPart) * GetParticleWeight(iPart)
      END IF
    END IF
  END IF
END DO

#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,EVib ,nSpecies, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,ERot ,nSpecies, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
  IF(DSMC%ElectronicModel.GT.0) CALL MPI_REDUCE(MPI_IN_PLACE,Eelec,nSpecies, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
ELSE
  CALL MPI_REDUCE(EVib        ,RD   ,nSpecies, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
  CALL MPI_REDUCE(ERot        ,RD   ,nSpecies, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
  IF(DSMC%ElectronicModel.GT.0) CALL MPI_REDUCE(Eelec       ,RD   ,nSpecies, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
END IF
#endif /*USE_MPI*/

! final computation is only done for the root
IF(MPIRoot)THEN
  ! Calc TVib, TRot
  DO iSpec = 1, nSpecies
    NumSpecTemp = NumSpec(iSpec)
    IF(((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)).AND.(NumSpecTemp.GT.0.0)) THEN
      IF (SpecDSMC(iSpec)%PolyatomicMol.AND.(SpecDSMC(iSpec)%Xi_Rot.EQ.3)) THEN
        IntTemp(iSpec,2) = 2.0*ERot(iSpec)/(3.0*BoltzmannConst*NumSpecTemp)  !Calc TRot
      ELSE
        IntTemp(iSpec,2) = ERot(iSpec)/(BoltzmannConst*NumSpecTemp)  !Calc TRot
      END IF
      IF (EVib(iSpec)/NumSpecTemp.GT.SpecDSMC(iSpec)%EZeroPoint) THEN
        IF (SpecDSMC(iSpec)%PolyatomicMol) THEN
          IntTemp(iSpec,1) = CalcTVibPoly(EVib(iSpec)/NumSpecTemp, iSpec)
        ELSE
          tempVib = (EVib(iSpec)/(NumSpecTemp*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)-DSMC%GammaQuant)
          IF ((tempVib.GT.0.0) &
            .OR.(EVib(iSpec)/(NumSpecTemp*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib).GT.DSMC%GammaQuant)) THEN
            IntTemp(iSpec,1) = SpecDSMC(iSpec)%CharaTVib/LOG(1 + 1/(EVib(iSpec) &
                              /(NumSpecTemp*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)-DSMC%GammaQuant))
          END IF
        END IF
      ELSE
        IntTemp(iSpec,1) = 0
      END IF
    ELSE
      IntTemp(iSpec,1) = 0
      IntTemp(iSpec,2) = 0
    END IF
    IF(DSMC%ElectronicModel.GT.0) THEN
      IF(NumSpecTemp.GT.0) THEN
        IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
          IntTemp(iSpec,3) = CalcTelec(Eelec(iSpec)/NumSpecTemp,iSpec)
        END IF
      ELSE
        IntEn(iSpec,3) = 0.0
      END IF
    END IF
    IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
      ! MacroParticleFactor is included in the case of RadialWeighting (also in combination with variable time step)
      IntEn(iSpec,1) = EVib(iSpec)
      IntEn(iSpec,2) = ERot(iSpec)
      IF(DSMC%ElectronicModel.GT.0) IntEn(iSpec,3) = Eelec(iSpec)
    ELSE
      IntEn(iSpec,1) = EVib(iSpec) * Species(iSpec)%MacroParticleFactor
      IntEn(iSpec,2) = ERot(iSpec) * Species(iSpec)%MacroParticleFactor
      IF(DSMC%ElectronicModel.GT.0) IntEn(iSpec,3) = Eelec(iSpec) * Species(iSpec)%MacroParticleFactor
    END IF
  END DO
  ! Sums of the energy values
  IF(nSpecAnalyze.GT.1) THEN
    IntEn(nSpecAnalyze,1) = SUM(IntEn(:,1))
    IntEn(nSpecAnalyze,2) = SUM(IntEn(:,2))
    IF(DSMC%ElectronicModel.GT.0) IntEn(nSpecAnalyze,3) = SUM(IntEn(:,3))
  END IF
END IF

END SUBROUTINE CalcIntTempsAndEn


SUBROUTINE CalcTransTemp(NumSpec, Temp)
!===================================================================================================================================
! calculate the translational temperature of each species
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst,Kelvin2eV
USE MOD_Particle_Vars         ,ONLY: nSpecies
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_Particle_Vars         ,ONLY: PartSpecies, PartState, Species, PDM
USE MOD_Particle_Analyze_Vars ,ONLY: nSpecAnalyze
USE MOD_DSMC_Vars             ,ONLY: DSMC, AmbipolElecVelo
USE MOD_Particle_Vars         ,ONLY: CalcBulkElectronTemp,BulkElectronTemp,BulkElectronTempSpecID
#if USE_MPI
USE MOD_SurfaceModel_Vars     ,ONLY: BulkElectronTempSEE,SurfModSEEelectronTempAutoamtic
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: NumSpec(:)    !< global number of REAL particles in domain
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Temp(:)       !< output value is already the GLOBAL temperature
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iSpec
REAL              :: TempDirec(nSpecies,3)
REAL              :: PartVandV2(nSpecies, 6), Mean_PartV2(nSpecies, 3), MeanPartV_2(nSpecies,3)
INTEGER           :: i
!===================================================================================================================================

! Compute velocity averages
Temp = 0.0
! Sum up velocity
PartVandV2 = 0.
DO i=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(i)) THEN
    PartVandV2(PartSpecies(i),1:3) = PartVandV2(PartSpecies(i),1:3) + PartState(4:6,i) * GetParticleWeight(i)
    PartVandV2(PartSpecies(i),4:6) = PartVandV2(PartSpecies(i),4:6) + PartState(4:6,i)**2 * GetParticleWeight(i)
    IF (DSMC%DoAmbipolarDiff) THEN
      IF(Species(PartSpecies(i))%ChargeIC.GT.0.0) THEN
        PartVandV2(DSMC%AmbiDiffElecSpec,1:3) = PartVandV2(DSMC%AmbiDiffElecSpec,1:3) + AmbipolElecVelo(i)%ElecVelo(1:3) * GetParticleWeight(i)
        PartVandV2(DSMC%AmbiDiffElecSpec,4:6) = PartVandV2(DSMC%AmbiDiffElecSpec,4:6) + AmbipolElecVelo(i)%ElecVelo(1:3)**2 * GetParticleWeight(i)
      END IF
    END IF
  END IF
END DO

#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,PartVandV2,nSpecies*6,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS, IERROR)
ELSE
  CALL MPI_REDUCE(PartVandV2  ,PartVandV2,nSpecies*6,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS, IERROR)
END IF
#endif /*USE_MPI*/

IF(MPIRoot)THEN
  DO iSpec=1, nSpecies
    IF(NumSpec(iSpec).NE.0) THEN
      ! Compute velocity averages
      MeanPartV_2(iSpec,1:3)  = (PartVandV2(iSpec,1:3) / NumSpec(iSpec))**2       ! < |v| >**2
      Mean_PartV2(iSpec,1:3)  =  PartVandV2(iSpec,4:6) / NumSpec(iSpec)           ! < |v|**2 >
    ELSE
      MeanPartV_2(iSpec,1:3) = 0.
      Mean_PartV2(iSpec,1:3) = 0.
    END IF
    ! Compute temperatures
    TempDirec(iSpec,1:3) = Species(iSpec)%MassIC * (Mean_PartV2(iSpec,1:3) - MeanPartV_2(iSpec,1:3)) &
         / BoltzmannConst ! Trans Temp calculation is limited to one species
    Temp(iSpec) = (TempDirec(iSpec,1) + TempDirec(iSpec,2) + TempDirec(iSpec,3))/3
    IF(nSpecAnalyze.GT.1)THEN
      Temp(nSpecAnalyze) = Temp(nSpecAnalyze) + Temp(iSpec)*NumSpec(iSpec)
    END IF
  END DO
  IF(nSpecAnalyze.GT.1)THEN
    IF(NumSpec(iSpec).NE.0) THEN
      Temp(nSpecAnalyze)= Temp(nSpecAnalyze) / NumSpec(nSpecAnalyze)
    ELSE
      Temp(nSpecAnalyze)= 0.
    END IF
  END IF
END IF

! For SEE model that used the bulk electron temperature, use the global electron temperature
IF(MPIRoot.AND.CalcBulkElectronTemp)THEN
  ASSOCIATE( Te => BulkElectronTemp ,&
             T  => Temp(BulkElectronTempSpecID)*Kelvin2eV)
    ! Smooth the bulk electron temperature by adjusting the difference by 50%
    Te = Te + 0.5 * (T - Te)
  END ASSOCIATE
END IF
#if USE_MPI
IF(CalcBulkElectronTemp) CALL MPI_BCAST(BulkElectronTemp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_PICLAS,iError)
IF(SurfModSEEelectronTempAutoamtic) BulkElectronTempSEE = BulkElectronTemp
#endif /*USE_MPI*/

END SUBROUTINE CalcTransTemp


REAL FUNCTION CalcTelec(MeanEelec, iSpec)
!===================================================================================================================================
!> Calculation of the electronic temperature (zero-point search)
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars  ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars     ,ONLY: SpecDSMC, DSMC
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

CalcTelec = 0.

SELECT CASE(DSMC%ElectronicModel)
CASE(1,2,4)
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
CASE(3)
  CalcTelec = MeanEelec / BoltzmannConst
END SELECT

RETURN

END FUNCTION CalcTelec


REAL FUNCTION CalcEelec(TElec, iSpec)
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
REAL, INTENT(IN)      :: TElec  !< Mean electronic energy
INTEGER, INTENT(IN)   :: iSpec      !< Species index
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER               :: ii
REAL                  :: TempRatio, SumOne, SumTwo        !< Sums of the electronic partition function
!===================================================================================================================================

IF (TElec.GT.0) THEN
  SumOne = 0.0
  SumTwo = 0.0
  DO ii = 0, SpecDSMC(iSpec)%MaxElecQuant-1
    TempRatio = SpecDSMC(iSpec)%ElectronicState(2,ii) / TElec
    IF(CHECKEXP(TempRatio)) THEN
      SumOne = SumOne + SpecDSMC(iSpec)%ElectronicState(1,ii) * EXP(-TempRatio)
      SumTwo = SumTwo + SpecDSMC(iSpec)%ElectronicState(1,ii) * SpecDSMC(iSpec)%ElectronicState(2,ii) * EXP(-TempRatio)
    END IF
  END DO
  CalcEelec = SumTwo / SumOne * BoltzmannConst
ELSE
  CalcEelec = 0. ! sup
END IF

RETURN

END FUNCTION CalcEelec


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


#if (PP_TimeDiscMethod==2 || PP_TimeDiscMethod==4 || PP_TimeDiscMethod==300 || PP_TimeDiscMethod==400 || (PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=509) || PP_TimeDiscMethod==120)
SUBROUTINE CalcRelaxProbRotVib(RotRelaxProb,VibRelaxProb)
!===================================================================================================================================
! Calculates global rotational and vibrational relaxation probability for PartAnalyse.csv
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars         ,ONLY: nSpecies
USE MOD_DSMC_Vars             ,ONLY: DSMC, VarVibRelaxProb, CollisMode
USE MOD_Mesh_Vars             ,ONLY: nElems, nGlobalElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: RotRelaxProb(2),VibRelaxProb(2)       !< output value is already the GLOBAL RelaxProbs
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iElem, iSpec
REAL                            :: PartNum
!===================================================================================================================================
IF(CollisMode.LT.2) RETURN

! Rot Relax Prob
IF(DSMC%RotRelaxProb.GE.2) THEN
  RotRelaxProb = 0.
  PartNum = 0.
  DO iSpec=1,nSpecies
    RotRelaxProb(1) = MAX(DSMC%CalcRotProb(iSpec,2),RotRelaxProb(1))
    RotRelaxProb(2) = RotRelaxProb(2) + DSMC%CalcRotProb(iSpec,1)
    PartNum  = PartNum + DSMC%CalcRotProb(iSpec,3)
  END DO
  IF(PartNum.GT.1) THEN
    RotRelaxProb(2) = RotRelaxProb(2)/PartNum
  ELSE
    RotRelaxProb(2) = 0
  END IF
#if USE_MPI
  IF(MPIRoot)THEN
    CALL MPI_REDUCE(MPI_IN_PLACE,RotRelaxProb(1),1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_PICLAS, IERROR)
    CALL MPI_REDUCE(MPI_IN_PLACE,RotRelaxProb(2),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS, IERROR)
    RotRelaxProb(2) = RotRelaxProb(2) / REAL(nProcessors)
  ELSE
    CALL MPI_REDUCE(RotRelaxProb(1),RotRelaxProb(1),1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_PICLAS, IERROR)
    CALL MPI_REDUCE(RotRelaxProb(2),RotRelaxProb(2),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS, IERROR)
  END IF
#endif /*USE_MPI*/
ELSE
  RotRelaxProb = DSMC%RotRelaxProb
END IF
! Vib Relax Prob
IF(DSMC%VibRelaxProb.EQ.2) THEN
  VibRelaxProb = 0.
  PartNum = 0.
  DO iSpec=1,nSpecies
    VibRelaxProb(1) = MAX(DSMC%CalcVibProb(iSpec,2),VibRelaxProb(1))
    DO iElem=1,nElems
      VibRelaxProb(2)=VibRelaxProb(2)+VarVibRelaxProb%ProbVibAv(iElem,iSpec)
    END DO
  END DO
  VibRelaxProb(2)=VibRelaxProb(2)/(REAL(nSpecies)*REAL(nGlobalElems))
#if USE_MPI
  IF(MPIRoot)THEN
    CALL MPI_REDUCE(MPI_IN_PLACE,VibRelaxProb(1),1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_PICLAS, IERROR)
    CALL MPI_REDUCE(MPI_IN_PLACE,VibRelaxProb(2),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS, IERROR)
  ELSE
    CALL MPI_REDUCE(VibRelaxProb(1),VibRelaxProb(1),1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_PICLAS, IERROR)
    CALL MPI_REDUCE(VibRelaxProb(2),VibRelaxProb(2),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS, IERROR)
  END IF
#endif /*USE_MPI*/
ELSE
  VibRelaxProb = 0.
#if USE_MPI
  IF(MPIRoot)THEN
    CALL MPI_REDUCE(MPI_IN_PLACE,DSMC%CalcVibProb(1:nSpecies,1),nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS, IERROR)
    CALL MPI_REDUCE(MPI_IN_PLACE,DSMC%CalcVibProb(1:nSpecies,2),nSpecies,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_PICLAS, IERROR)
    CALL MPI_REDUCE(MPI_IN_PLACE,DSMC%CalcVibProb(1:nSpecies,3),nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS, IERROR)
  ELSE
    CALL MPI_REDUCE(DSMC%CalcVibProb(1:nSpecies,1),DSMC%CalcVibProb(1:nSpecies,1),nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS, IERROR)
    CALL MPI_REDUCE(DSMC%CalcVibProb(1:nSpecies,2),DSMC%CalcVibProb(1:nSpecies,2),nSpecies,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_PICLAS, IERROR)
    CALL MPI_REDUCE(DSMC%CalcVibProb(1:nSpecies,3),DSMC%CalcVibProb(1:nSpecies,3),nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS, IERROR)
  END IF
#endif /*USE_MPI*/
  IF(MPIRoot)THEN
    VibRelaxProb(1) = MAXVAL(DSMC%CalcVibProb(1:nSpecies,2))
    IF(SUM(DSMC%CalcVibProb(1:nSpecies,3)).GT.0) THEN
      VibRelaxProb(2) = SUM(DSMC%CalcVibProb(1:nSpecies,1))/SUM(DSMC%CalcVibProb(1:nSpecies,3))
    END IF
  END IF
END IF

END SUBROUTINE CalcRelaxProbRotVib
#endif


SUBROUTINE CalcVelocities(PartVtrans, PartVtherm,NumSpec,SimNumSpec)
!===================================================================================================================================
! Calculates the drift and eigen velocity of all particles: PartVtotal = PartVtrans + PartVtherm
! PartVtrans(nSpecies,4) ! macroscopic velocity (drift velocity) A. Frohn: kinetische Gastheorie
! PartVtherm(nSpecies,4) ! microscopic velocity (eigen velocity)
!
! Note that the thermal velocity corresponds to the root mean square of the total velocity (in three dimensions), which is given by
!
!      v_th = SQRT(3 * kB * T / m)
!
! with kB : Boltzmann's constant
!      T  : temperature
!      m  : mass of the particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Analyze_Vars ,ONLY: VeloDirs
USE MOD_Particle_Vars         ,ONLY: PartState, PartSpecies, PDM, nSpecies, PartMPF, usevMPF
USE MOD_Particle_Analyze_Vars ,ONLY: nSpecAnalyze
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: NumSpec(nSpecAnalyze)
INTEGER(KIND=IK),INTENT(IN)    :: SimNumSpec(nSpecAnalyze)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)               :: PartVtrans(nSpecies,4), PartVtherm(nSpecies,4)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iSpec
INTEGER                        :: i
INTEGER                        :: dir
#if USE_MPI
REAL                           :: RD(nSpecies*4)
#endif /*USE_MPI*/
!===================================================================================================================================
! Compute velocity averages
  PartVtrans = 0.
  PartVtherm = 0.

  ! compute trans. velocity
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      DO dir = 1,3
        IF (VeloDirs(dir) .OR. VeloDirs(4)) THEN
          IF (usevMPF) THEN
            PartVtrans(PartSpecies(i),dir) = PartVtrans(PartSpecies(i),dir) + PartState(dir+3,i) * PartMPF(i)
          ELSE
            PartVtrans(PartSpecies(i),dir) = PartVtrans(PartSpecies(i),dir) + PartState(dir+3,i)
          END IF
        END IF
      END DO
    END IF
  END DO

#if USE_MPI
  IF(MPIRoot)THEN
    CALL MPI_REDUCE(MPI_IN_PLACE,PartVtrans ,4*nSpecies,MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
  ELSE
    CALL MPI_REDUCE(PartVtrans  ,RD         ,4*nSpecies,MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
  END IF
#endif /*USE_MPI*/

  IF(MPIRoot)THEN
    IF (usevMPF) THEN
      DO dir = 1,3
        IF (VeloDirs(dir) .OR. VeloDirs(4)) THEN
          DO iSpec = 1,nSpecies
            IF(NumSpec(iSpec).EQ.0)THEN
              PartVtrans(iSpec,dir) = 0.
            ELSE
              PartVtrans(iSpec,dir) = PartVtrans(iSpec,dir)/NumSpec(iSpec)
            END IF
          END DO ! iSpec = 1,nSpecies
        END IF
      END DO
    ELSE !no vMPF
      DO dir = 1,3
        IF (VeloDirs(dir) .OR. VeloDirs(4)) THEN
          DO iSpec = 1,nSpecies
            IF(SimNumSpec(iSpec).EQ.0)THEN
              PartVtrans(iSpec,dir) = 0.
            ELSE
              PartVtrans(iSpec,dir) = PartVtrans(iSpec,dir)/REAL(SimNumSpec(iSpec),8)
            END IF
          END DO ! iSpec = 1,nSpecies
        END IF
      END DO
    END IF !usevMPF
  END IF

#if USE_MPI
  CALL MPI_BCAST(PartVtrans,4*nSpecies, MPI_DOUBLE_PRECISION,0,MPI_COMM_PICLAS,iERROR)
#endif /*USE_MPI*/

  ! calculate thermal velocity
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      DO dir = 1,3
        IF (VeloDirs(dir) .OR. VeloDirs(4)) THEN
          IF (usevMPF) THEN
            PartVtherm(PartSpecies(i),dir) = PartVtherm(PartSpecies(i),dir) + PartMPF(i) * &
                (PartState(dir+3,i) - PartVtrans(PartSpecies(i),dir))*(PartState(dir+3,i) - PartVtrans(PartSpecies(i),dir))
          ELSE
            PartVtherm(PartSpecies(i),dir) = PartVtherm(PartSpecies(i),dir) + &
                (PartState(dir+3,i) - PartVtrans(PartSpecies(i),dir))*(PartState(dir+3,i) - PartVtrans(PartSpecies(i),dir))
          END IF
        END IF
      END DO
    END IF
  END DO

#if USE_MPI
  IF(MPIRoot)THEN
    CALL MPI_REDUCE(MPI_IN_PLACE,PartVtherm,4*nSpecies,MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
  ELSE
    CALL MPI_REDUCE(PartVtherm  ,RD        ,4*nSpecies,MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
  END IF
#endif /*USE_MPI*/

  IF(MPIRoot)THEN
    DO dir = 1,3
      IF (VeloDirs(dir) .OR. VeloDirs(4)) THEN
        DO iSpec = 1,nSpecies
          IF (usevMPF) THEN
            IF(NumSpec(iSpec).EQ.0)THEN
              PartVtherm(iSpec,dir)=0.
            ELSE
              PartVtherm(iSpec,dir)=PartVtherm(iSpec,dir)/NumSpec(iSpec)
            END IF
          ELSE
            IF(SimNumSpec(iSpec).EQ.0)THEN
              PartVtherm(iSpec,dir)=0.
            ELSE
              PartVtherm(iSpec,dir)=PartVtherm(iSpec,dir)/REAL(SimNumSpec(iSpec),8)
            END IF
          END IF
        END DO ! iSpec = 1,nSpecies
      END IF
    END DO
 !   calc absolute value
    IF (VeloDirs(4)) THEN
      PartVtrans(:,4) = SQRT(PartVtrans(:,1)*PartVtrans(:,1) + PartVtrans(:,2)*PartVtrans(:,2) + PartVtrans(:,3)*PartVtrans(:,3))
      PartVtherm(:,4) = PartVtherm(:,1) + PartVtherm(:,2) + PartVtherm(:,3)
    END IF
    PartVtherm(:,:) = SQRT(PartVtherm(:,:))
  END IF
END SUBROUTINE CalcVelocities


#if (PP_TimeDiscMethod==4)
SUBROUTINE CollRates(CRate)
!===================================================================================================================================
!> Calculate the collision rate per species pairing by diving the summed up variables by the current timestep
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars             ,ONLY: CollInf, DSMC
USE MOD_TimeDisc_Vars         ,ONLY: dt, iter
USE MOD_Particle_Vars         ,ONLY: VarTimeStep
USE MOD_Particle_Analyze_Vars ,ONLY: PartAnalyzeStep
USE MOD_Particle_TimeStep     ,ONLY: GetSpeciesTimeStep
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)              :: CRate(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iCase
REAL                          :: dtVar
!===================================================================================================================================

#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,DSMC%NumColl,CollInf%NumCase + 1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
ELSE
  CALL MPI_REDUCE(DSMC%NumColl,DSMC%NumColl,CollInf%NumCase + 1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
END IF
#endif /*USE_MPI*/

IF(MPIRoot)THEN
  DSMC%NumColl(CollInf%NumCase + 1) = SUM(DSMC%NumColl(1:CollInf%NumCase))
  DO iCase=1, CollInf%NumCase
    ! Species-specific time step
    IF(VarTimeStep%UseSpeciesSpecific.AND..NOT.VarTimeStep%DisableForMCC) THEN
      dtVar = dt * GetSpeciesTimeStep(iCase)
    ELSE
      dtVar = dt
    END IF
    CRate(iCase) =  DSMC%NumColl(iCase) / dtVar
  END DO
  ! Total collision rate is the sum of the case-specific rates
  CRate(CollInf%NumCase + 1) = SUM(CRate(1:CollInf%NumCase))
  ! Consider Part-AnalyzeStep
  IF(PartAnalyzeStep.GT.1)THEN
    IF(PartAnalyzeStep.EQ.HUGE(PartAnalyzeStep))THEN
      DO iCase=1, CollInf%NumCase + 1
        CRate(iCase) = CRate(iCase) / iter
      END DO ! iCase=1, CollInf%NumCase + 1
    ELSE
      DO iCase=1, CollInf%NumCase + 1
        CRate(iCase) = CRate(iCase) / MIN(PartAnalyzeStep,iter)
      END DO ! iCase=1, CollInf%NumCase + 1
    END IF
  END IF
END IF

DSMC%NumColl = 0.

END SUBROUTINE CollRates


SUBROUTINE CalcRelaxRates(NumSpec,VibRelaxProbCase)
!===================================================================================================================================
! Calculates global rotational and vibrational relaxation probability for PartAnalyse.csv
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars         ,ONLY: nSpecies, Species, VarTimeStep
USE MOD_DSMC_Vars             ,ONLY: CollisMode, CollInf
USE MOD_MCC_Vars              ,ONLY: SpecXSec, XSec_Relaxation
USE MOD_Particle_Mesh_Vars    ,ONLY: MeshVolume
USE MOD_TimeDisc_Vars         ,ONLY: dt, iter
USE MOD_Particle_Analyze_Vars ,ONLY: PartAnalyzeStep
USE MOD_Particle_TimeStep     ,ONLY: GetSpeciesTimeStep
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(IN)                 :: NumSpec(:)
REAL,INTENT(OUT)                :: VibRelaxProbCase(CollInf%NumCase)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSpec, iCase, jSpec
REAL                            :: MPF_1, MPF_2, dtVar
!===================================================================================================================================
IF(CollisMode.LT.2) RETURN

#if USE_MPI
IF(XSec_Relaxation) THEN
  IF(MPIRoot)THEN
    CALL MPI_REDUCE(MPI_IN_PLACE,SpecXSec(:)%VibCount,CollInf%NumCase,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS, IERROR)
  ELSE
    CALL MPI_REDUCE(SpecXSec(:)%VibCount,SpecXSec(:)%VibCount,CollInf%NumCase,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS, IERROR)
  END IF
END IF
#endif /*USE_MPI*/
VibRelaxProbCase = 0.
IF(MPIRoot)THEN
  IF(XSec_Relaxation) THEN
    DO iSpec=1,nSpecies
      IF(NumSpec(iSpec).LE.0.0) CYCLE
      MPF_1 = Species(iSpec)%MacroParticleFactor
      DO jSpec = iSpec, nSpecies
        IF(NumSpec(jSpec).LE.0.0) CYCLE
        MPF_2 = Species(jSpec)%MacroParticleFactor
        iCase = CollInf%Coll_Case(iSpec,jSpec)
        ! Species-specific time step
        IF(VarTimeStep%UseSpeciesSpecific.AND..NOT.VarTimeStep%DisableForMCC) THEN
          dtVar = dt * GetSpeciesTimeStep(iCase)
        ELSE
          dtVar = dt
        END IF
        VibRelaxProbCase(iCase) = SpecXSec(iCase)%VibCount * MPF_1 * MeshVolume &
                                  / (dtVar * MPF_1*NumSpec(iSpec) * MPF_2*NumSpec(jSpec))
      END DO
    END DO
  END IF
  ! Consider Part-AnalyzeStep
  IF(PartAnalyzeStep.GT.1) THEN
    IF(PartAnalyzeStep.EQ.HUGE(PartAnalyzeStep))THEN
      DO iCase = 1, CollInf%NumCase
        VibRelaxProbCase(iCase) = VibRelaxProbCase(iCase) / iter
      END DO
    ELSE
      DO iCase = 1, CollInf%NumCase
        VibRelaxProbCase(iCase) = VibRelaxProbCase(iCase) / MIN(PartAnalyzeStep,iter)
      END DO
    END IF
  END IF
END IF
SpecXSec(:)%VibCount = 0.

END SUBROUTINE CalcRelaxRates


SUBROUTINE CalcRelaxRatesElec(ElecRelaxRate)
!===================================================================================================================================
!> Calculates the global electronic relaxation rate per case per level for PartAnalyse.csv
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars         ,ONLY: VarTimeStep
USE MOD_DSMC_Vars             ,ONLY: CollInf
USE MOD_MCC_Vars              ,ONLY: SpecXSec
USE MOD_TimeDisc_Vars         ,ONLY: dt, iter
USE MOD_Particle_Analyze_Vars ,ONLY: PartAnalyzeStep
USE MOD_Particle_TimeStep     ,ONLY: GetSpeciesTimeStep
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT),ALLOCATABLE  :: ElecRelaxRate(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iCase, iLevel, MaxLevel
REAL                          :: dtVar
!===================================================================================================================================

MaxLevel = MAXVAL(SpecXSec(:)%NumElecLevel)
ALLOCATE(ElecRelaxRate(CollInf%NumCase,MaxLevel))
ElecRelaxRate = 0.

DO iCase=1, CollInf%NumCase
  IF(SpecXSec(iCase)%UseElecXSec) THEN
    DO iLevel = 1, SpecXSec(iCase)%NumElecLevel
      ElecRelaxRate(iCase,iLevel) =  SpecXSec(iCase)%ElecLevel(iLevel)%Counter
    END DO
  END IF
END DO

#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,ElecRelaxRate,CollInf%NumCase*MaxLevel,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
ELSE
  CALL MPI_REDUCE(ElecRelaxRate,ElecRelaxRate,CollInf%NumCase*MaxLevel,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
END IF
#endif /*USE_MPI*/

IF(MPIRoot)THEN
  DO iCase=1, CollInf%NumCase
    IF(SpecXSec(iCase)%UseElecXSec) THEN
      ! Species-specific time step
      IF(VarTimeStep%UseSpeciesSpecific.AND..NOT.VarTimeStep%DisableForMCC) THEN
        dtVar = dt * GetSpeciesTimeStep(iCase)
      ELSE
        dtVar = dt
      END IF
      ElecRelaxRate(iCase,:) =  ElecRelaxRate(iCase,:) / dtVar
    END IF
  END DO
  ! Consider Part-AnalyzeStep
  IF(PartAnalyzeStep.GT.1)THEN
    IF(PartAnalyzeStep.EQ.HUGE(PartAnalyzeStep))THEN
      ElecRelaxRate = ElecRelaxRate / iter
    ELSE
      ElecRelaxRate = ElecRelaxRate / MIN(PartAnalyzeStep,iter)
    END IF
  END IF
END IF

DO iCase=1, CollInf%NumCase
  IF(SpecXSec(iCase)%UseElecXSec) THEN
    DO iLevel = 1, SpecXSec(iCase)%NumElecLevel
      SpecXSec(iCase)%ElecLevel(iLevel)%Counter = 0.
    END DO
  END IF
END DO

END SUBROUTINE CalcRelaxRatesElec


SUBROUTINE ReacRates(NumSpec, RRate)
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars             ,ONLY: ChemReac, DSMC
USE MOD_TimeDisc_Vars         ,ONLY: dt, iter
USE MOD_Particle_Vars         ,ONLY: Species, nSpecies, VarTimeStep
USE MOD_Particle_Mesh_Vars    ,ONLY: MeshVolume
USE MOD_Particle_Analyze_Vars ,ONLY: PartAnalyzeStep
USE MOD_Particle_TimeStep     ,ONLY: GetSpeciesTimeStep
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(IN)                 :: NumSpec(:)
REAL,INTENT(OUT)                :: RRate(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iReac, iCase
REAL                            :: dtVar
#if USE_MPI
REAL                            :: RD(1:ChemReac%NumOfReact)
#endif /*USE_MPI*/
!===================================================================================================================================

#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE    ,ChemReac%NumReac,ChemReac%NumOfReact,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
ELSE
  CALL MPI_REDUCE(ChemReac%NumReac,RD              ,ChemReac%NumOfReact,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
END IF
#endif /*USE_MPI*/

IF(MPIRoot)THEN
  DO iReac=1, ChemReac%NumOfReact
    iCase = ChemReac%ReactCase(iReac)
    ! Species-specific time step
    IF(VarTimeStep%UseSpeciesSpecific.AND..NOT.VarTimeStep%DisableForMCC) THEN
      dtVar = dt * GetSpeciesTimeStep(iCase)
    ELSE
      dtVar = dt
    END IF
    IF ((NumSpec(ChemReac%Reactants(iReac,1)).GT.0).AND.(NumSpec(ChemReac%Reactants(iReac,2)).GT.0)) THEN
      IF(ChemReac%Reactants(iReac,3).NE.0) THEN
        ! Recombination reactions with 3 reactants
        IF (DSMC%ReservoirRateStatistic) THEN ! Calculation of rate constant through actual number of allowed reactions
          RRate(iReac) = ChemReac%NumReac(iReac) * Species(ChemReac%Products(iReac,1))%MacroParticleFactor &
                     * MeshVolume**2 / (dtVar &
                     * Species(ChemReac%Reactants(iReac,1))%MacroParticleFactor * NumSpec(ChemReac%Reactants(iReac,1)) &
                     * Species(ChemReac%Reactants(iReac,2))%MacroParticleFactor * NumSpec(ChemReac%Reactants(iReac,2)) &
                     * Species(ChemReac%Reactants(iReac,3))%MacroParticleFactor * NumSpec(nSpecies+1))
        ! Calculation of rate constant through mean reaction probability (using mean reaction prob and sum of coll prob)
        ELSEIF(ChemReac%ReacCount(iReac).GT.0) THEN
          RRate(iReac) = ChemReac%NumReac(iReac) * ChemReac%ReacCollMean(iCase) * MeshVolume**2 &
               * Species(ChemReac%Reactants(iReac,1))%MacroParticleFactor / (dtVar * ChemReac%ReacCount(iReac)             &
               * Species(ChemReac%Reactants(iReac,1))%MacroParticleFactor*NumSpec(ChemReac%Reactants(iReac,1))     &
               * Species(ChemReac%Reactants(iReac,2))%MacroParticleFactor*NumSpec(ChemReac%Reactants(iReac,2))    &
               * Species(ChemReac%Reactants(iReac,3))%MacroParticleFactor*NumSpec(nSpecies+1))
        END IF
      ELSE
        ! Regular reactions with 2 reactants (dissociation, ionization, exchange)
        IF (DSMC%ReservoirRateStatistic) THEN ! Calculation of rate constant through actual number of allowed reactions
          RRate(iReac) = ChemReac%NumReac(iReac) * Species(ChemReac%Products(iReac,1))%MacroParticleFactor &
                       * MeshVolume / (dtVar &
                       * Species(ChemReac%Reactants(iReac,1))%MacroParticleFactor*NumSpec(ChemReac%Reactants(iReac,1)) &
                       * Species(ChemReac%Reactants(iReac,2))%MacroParticleFactor*NumSpec(ChemReac%Reactants(iReac,2)))
        ! Calculation of rate constant through mean reaction probability (using mean reaction prob and sum of coll prob)
        ELSEIF(ChemReac%ReacCount(iReac).GT.0) THEN
          RRate(iReac) = ChemReac%NumReac(iReac) * ChemReac%ReacCollMean(iCase) &
               * Species(ChemReac%Reactants(iReac,1))%MacroParticleFactor* MeshVolume / (dtVar * ChemReac%ReacCount(iReac) &
               * Species(ChemReac%Reactants(iReac,1))%MacroParticleFactor*NumSpec(ChemReac%Reactants(iReac,1))         &
               * Species(ChemReac%Reactants(iReac,2))%MacroParticleFactor*NumSpec(ChemReac%Reactants(iReac,2)))
        END IF
      END IF
    END IF
  END DO
END IF
ChemReac%NumReac = 0.
ChemReac%ReacCount = 0
ChemReac%ReacCollMean = 0.0
! Consider Part-AnalyzeStep
IF(PartAnalyzeStep.GT.1)THEN
  IF(PartAnalyzeStep.EQ.HUGE(PartAnalyzeStep))THEN
    DO iReac=1, ChemReac%NumOfReact
      RRate(iReac) = RRate(iReac) / iter
    END DO ! iReac=1, ChemReac%NumOfReact
  ELSE
    DO iReac=1, ChemReac%NumOfReact
      RRate(iReac) = RRate(iReac) / REAL(MIN(PartAnalyzeStep,iter))
    END DO ! iReac=1, ChemReac%NumOfReact
  END IF
END IF

END SUBROUTINE ReacRates
#endif


SUBROUTINE CalcPowerDensity()
!===================================================================================================================================
! Used to average the source terms per species
!   * compute the power density of the considered species
!     scalar product of < j, E >, with j- current density and E - electric field
!   * compute the charge density of the considered species
!   * compute the current density of the considered species
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Timeaverage_Vars ,ONLY: DoPowerDensity,PowerDensity
USE MOD_Particle_Vars    ,ONLY: nSpecies,PartSpecies,PDM
USE MOD_PICDepo_Vars     ,ONLY: PartSource
USE MOD_Part_RHS         ,ONLY: PartVeloToImp
USE MOD_Preproc
USE MOD_PICDepo          ,ONLY: Deposition
USE MOD_Mesh_Tools       ,ONLY: GetCNElemID
#if ! (USE_HDG)
USE MOD_DG_Vars          ,ONLY: U
#else
#if PP_nVar==1
USE MOD_Equation_Vars    ,ONLY: E
#else
#endif
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iSpec,iSpec2
INTEGER              :: iElem,i,j,k,iPart
LOGICAL              :: doParticle(1:PDM%MaxParticleNumber)
!===================================================================================================================================

iSpec2=0
PowerDensity=0.
DO iSpec=1,nSpecies
  IF(.NOT.DoPowerDensity(iSpec)) CYCLE
  iSpec2=iSpec2+1
  ! mark particle
  DoParticle(:)=.FALSE.
  DO iPart=1,PDM%ParticleVecLength
    IF(PDM%ParticleInside(iPart))THEN
      IF(PartSpecies(iPart).EQ.iSpec)THEN
        DoParticle(iPart)=.TRUE.
      END IF
    END IF ! ParticleInside
  END DO ! iPart


  ! map particle from gamma v to v
  CALL PartVeloToImp(VeloToImp=.FALSE.,doParticle_In=DoParticle(1:PDM%ParticleVecLength))
  ! compute source terms
  ! compute particle source terms on field solver of considered species
  CALL Deposition(doParticle_In=DoParticle(1:PDM%ParticleVecLength))
  ! map particle from v to gamma v
  CALL PartVeloToImp(VeloToImp=.TRUE.,doParticle_In=DoParticle(1:PDM%ParticleVecLength))

  ! compute power density
  DO iElem=1,PP_nElems
    DO k=0,PP_N
      DO j=0,PP_N
        DO i=0,PP_N
          ! 1:3 PowerDensity, 4 charge density
#if !(USE_HDG)
          PowerDensity(1,i,j,k,iElem,iSpec2)=PartSource(1,i,j,k,iElem)*U(1,i,j,k,iElem)
          PowerDensity(2,i,j,k,iElem,iSpec2)=PartSource(2,i,j,k,iElem)*U(2,i,j,k,iElem)
          PowerDensity(3,i,j,k,iElem,iSpec2)=PartSource(3,i,j,k,iElem)*U(3,i,j,k,iElem)
          PowerDensity(4,i,j,k,iElem,iSpec2)=PartSource(4,i,j,k,iElem)
#else
#if PP_nVar==1
          PowerDensity(1,i,j,k,iElem,iSpec2)=PartSource(1,i,j,k,iElem)*E(1,i,j,k,iElem)
          PowerDensity(2,i,j,k,iElem,iSpec2)=PartSource(2,i,j,k,iElem)*E(2,i,j,k,iElem)
          PowerDensity(3,i,j,k,iElem,iSpec2)=PartSource(3,i,j,k,iElem)*E(3,i,j,k,iElem)
#else
          PowerDensity(1:3,i,j,k,iElem,iSpec2)=0.
#endif
          PowerDensity(4,i,j,k,iElem,iSpec2)=PartSource(4,i,j,k,iElem)
#endif
          ! 5:7 current density
          PowerDensity(5,i,j,k,iElem,iSpec2)=PartSource(1,i,j,k,iElem)
          PowerDensity(6,i,j,k,iElem,iSpec2)=PartSource(2,i,j,k,iElem)
          PowerDensity(7,i,j,k,iElem,iSpec2)=PartSource(3,i,j,k,iElem)
        END DO ! i=0,PP_N
      END DO ! j=0,PP_N
    END DO ! k=0,PP_N
  END DO ! iElem=1,PP_nElems
END DO

END SUBROUTINE CalcPowerDensity


SUBROUTINE CalculatePlasmaFrequencyCell()
!===================================================================================================================================
! use the number of electron density to compute the plasma frequency per cell using fixed and global values for the
! electron charge, electronmass and eps0
! CAUTION: if c!=3e8 m/s the computed frequency may be wrong
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars  ,ONLY:ElectronDensityCell,PlasmaFrequencyCell
USE MOD_Globals_Vars           ,ONLY:ElementaryCharge,ElectronMass
USE MOD_Globals_Vars           ,ONLY:eps0
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iElem
!===================================================================================================================================

! nullify
PlasmaFrequencyCell=0.

! loop over all elements and compute the plasma frequency with the use of the electron density
DO iElem=1,PP_nElems
  PlasmaFrequencyCell(iElem) = SQRT((ElectronDensityCell(iElem)*ElementaryCharge*ElementaryCharge)/(ElectronMass*Eps0))
END DO ! iElem=1,PP_nElems

END SUBROUTINE CalculatePlasmaFrequencyCell


SUBROUTINE CalculateCyclotronFrequencyAndRadiusCell()
!===================================================================================================================================
! Determine the (relativistic) electron cyclotron frequency in each cell, which can be calculate without electrons present in the
! cell in the classical (non-relativistic) limit. If electrons are present, their velocity is used to calculate the Lorentz factor
! gamma. From the cyclotron frequency and the electron velocity, the Larmor or gyroradius is calculated
!
! ------------------------------------------------
! omega_c = e*B / m_e     (non-relativistic)
!
!   omega_c: cyclotron frequency
!         e: elementary charge (of an electron, absolute value)
!         B: magnitude of the magnetic flux density at the electron's position
!       m_e: electron rest mass
!
! ------------------------------------------------
! omega_c = e*B / (gamma*m_e) = e*B / (sqrt(1-v_e^2/c^2)*m_e)
!
!   omega_c: cyclotron frequency     (relativistic)
!         e: elementary charge (of an electron, absolute value)
!         B: magnitude of the magnetic flux density at the electron's position
!     gamma: Lorentz factor
!       m_e: electron rest mass
!       v_e: magnitude of velocity
!         c: speed of light
!
! ------------------------------------------------
! r = v_e / omega_c
!   omega_c: cyclotron frequency
!       v_e: magnitude of velocity
!
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Preproc
USE MOD_Globals                ,ONLY: PARTISELECTRON,VECNORM,DOTPRODUCT
USE MOD_Globals_Vars           ,ONLY: c2_inv,RelativisticLimit
USE MOD_Particle_Vars          ,ONLY: PartState
USE MOD_Particle_Analyze_Vars  ,ONLY: CyclotronFrequencyMaxCell,CyclotronFrequencyMinCell,GyroradiusMinCell,GyroradiusMaxCell
USE MOD_Globals_Vars           ,ONLY: ElementaryCharge,ElectronMass
USE MOD_Particle_Vars          ,ONLY: PDM, PEM
USE MOD_PICInterpolation_tools ,ONLY: GetExternalFieldAtParticle,GetInterpolatedFieldPartPos,GetEMField
USE MOD_Interpolation_Vars     ,ONLY: xGP
USE MOD_Mesh_Vars              ,ONLY: Elem_xGP
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iElem,iGlobElem,iPart,i,j,k
REAL    :: PartV,partV2,field(6),B,gamma1,omega_c,gyroradius
LOGICAL :: SetFrequency,SetRadius
!===================================================================================================================================

! nullify
CyclotronFrequencyMaxCell=0.
GyroradiusMinCell=HUGE(1.)
CyclotronFrequencyMinCell=HUGE(1.)
GyroradiusMaxCell=0.

ASSOCIATE( e   => ElementaryCharge,&
           m_e => ElectronMass )
  ! 1.) Loop all particles and check only if electrons
  DO iPart=1,PDM%ParticleVecLength
    IF(.NOT.PDM%ParticleInside(iPart)) CYCLE
    IF(PARTISELECTRON(iPart))THEN
      SetFrequency = .FALSE. ! Initialize
      SetRadius    = .FALSE. ! Initialize
      ! Get magnitude of the electron's velocity and the magnetic field at its location
      PartV  = VECNORM(PartState(4:6,iPart)) ! velocity magnitude
      partV2 = PartV*PartV
      iGlobElem  = PEM%GlobalElemID(iPart)
      iElem  = PEM%LocalElemID(iPart)
      IF (partV2.LT.RelativisticLimit)THEN ! |v| < 1000000 when speed of light is 299792458
        field(1:6)   = GetExternalFieldAtParticle(PartState(1:3,iPart)) + GetInterpolatedFieldPartPos(iGlobElem,iPart)
        B            = VECNORM(field(4:6))
        omega_c      = e*B/m_e
        SetFrequency = .TRUE.
        IF(omega_c.GT.0.) SetRadius = .TRUE.
      ELSE
        gamma1=partV2*c2_inv
        ! Sanity check: Lorentz factor must be below 1.0: gamma in [0,1)
        IF(gamma1.GE.1.0)THEN
          ! don't store this value as cyclotron frequency, keep the zero or an already correctly set value
        ELSE
          field(1:6)   = GetExternalFieldAtParticle(PartState(1:3,iPart)) + GetInterpolatedFieldPartPos(iGlobElem,iPart)
          gamma1       = 1.0/SQRT(1.-gamma1)
          B            = VECNORM(field(4:6))
          omega_c      = e*B/(gamma1*m_e)
          SetFrequency = .TRUE.
          IF(omega_c.GT.0.) SetRadius = .TRUE.
        END IF ! gamma1.GE.1.0
      END IF ! partV2.LT.RelativisticLimit

      ! Check if values were calculated for this particle
      IF(SetFrequency)THEN
        CyclotronFrequencyMaxCell(iElem) = MAX(CyclotronFrequencyMaxCell(iElem), omega_c)
        CyclotronFrequencyMinCell(iElem) = MIN(CyclotronFrequencyMinCell(iElem), omega_c)
        IF(SetRadius)THEN
          gyroradius = PartV / omega_c
          GyroradiusMaxCell(iElem) = MAX(GyroradiusMaxCell(iElem), gyroradius)
          GyroradiusMinCell(iElem) = MIN(GyroradiusMinCell(iElem), gyroradius)
        END IF ! SetRadius
      END IF ! SetFrequency
    END IF ! PARTISELECTRON(iPart)
  END DO ! iPart=1,PDM%ParticleVecLength

  ! 2.) Loop over all elements and compute the Cyclotron frequency for elements, which have not been done yet (get magnetic field
  !     at all DOF and use the maximum)
  DO iElem=1,PP_nElems

    ! Max
    IF(ABS(CyclotronFrequencyMaxCell(iElem)).LE.0.)THEN
      DO k=0,PP_N
        DO j=0,PP_N
          DO i=0,PP_N
            ASSOCIATE( x => Elem_xGP(1,i,j,k,iElem), y => Elem_xGP(2,i,j,k,iElem), z => Elem_xGP(3,i,j,k,iElem))
              field(1:6) = GetExternalFieldAtParticle((/x,y,z/)) + GetEMField(iElem,(/xGP(i),xGP(j),xGP(k)/))
              B = VECNORM(field(4:6))
              CyclotronFrequencyMaxCell(iElem) = MAX(CyclotronFrequencyMaxCell(iElem), e*B/(m_e) )
            END ASSOCIATE
          END DO ! i
        END DO ! j
      END DO ! k
    END IF ! ABS(CyclotronFrequencyMaxCell(iElem)).LE.0.

    ! Min
    IF(ABS(CyclotronFrequencyMinCell(iElem)).EQ.HUGE(1.))THEN
      DO k=0,PP_N
        DO j=0,PP_N
          DO i=0,PP_N
            ASSOCIATE( x => Elem_xGP(1,i,j,k,iElem), y => Elem_xGP(2,i,j,k,iElem), z => Elem_xGP(3,i,j,k,iElem))
              field(1:6) = GetExternalFieldAtParticle((/x,y,z/)) + GetEMField(iElem,(/xGP(i),xGP(j),xGP(k)/))
              B = VECNORM(field(4:6))
              CyclotronFrequencyMinCell(iElem) = MIN(CyclotronFrequencyMinCell(iElem), e*B/(m_e) )
            END ASSOCIATE
          END DO ! i
        END DO ! j
      END DO ! k
    END IF ! ABS(CyclotronFrequencyMinCell(iElem)).LE.0.

    ! Sanity check
    IF(GyroradiusMinCell(iElem).EQ.HUGE(1.)) GyroradiusMinCell(iElem)=0.
  END DO ! iElem=1,PP_nElems
END ASSOCIATE

END SUBROUTINE CalculateCyclotronFrequencyAndRadiusCell


SUBROUTINE CalculatePICTimeStepCyclotron()
!===================================================================================================================================
! use the gyro frequency per cell to estimate the pic time step
! Factor 0.05 = 1/20 from: Qin "Why is Boris algorithm so good?" (2013), PHYSICS OF PLASMAS 20, 084503 (2013)
!
! dt >= 0.05 / omega_c
!   omega_c: electron cyclotron frequency
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars ,ONLY: CyclotronFrequencyMaxCell,PICTimeStepCyclotronCell
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iElem
!===================================================================================================================================

! nullify
PICTimeStepCyclotronCell=0.

! loop over all elements and compute the PIC-timestep with the plasma frequency
DO iElem=1,PP_nElems
  IF(CyclotronFrequencyMaxCell(iElem).LE.0.) CYCLE
  PICTimeStepCyclotronCell(iElem) = 0.05 / CyclotronFrequencyMaxCell(iElem)
END DO ! iElem=1,PP_nElems

END SUBROUTINE CalculatePICTimeStepCyclotron


SUBROUTINE CalculatePICTimeStepCell()
!===================================================================================================================================
! use the plasma frequency per cell to estimate the pic time step
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars ,ONLY: PlasmaFrequencyCell,PICTimeStepCell,PICValidPlasmaCell,PICTimeCellResolved
USE MOD_TimeDisc_Vars         ,ONLY: dt
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iElem
!===================================================================================================================================

! nullify
PICTimeStepCell=0.
PICTimeCellResolved=0

! loop over all elements and compute the PIC-timestep with the plasma frequency
DO iElem=1,PP_nElems
  IF(PlasmaFrequencyCell(iElem).LE.0.) CYCLE
  PICTimeStepCell(iElem) = 0.2 / PlasmaFrequencyCell(iElem)
  IF((PICValidPlasmaCell(iElem).EQ.1).AND.(PICTimeStepCell(iElem).GE.dt)) PICTimeCellResolved = PICTimeCellResolved + 1
END DO ! iElem=1,PP_nElems

END SUBROUTINE CalculatePICTimeStepCell


SUBROUTINE CalculateDebyeLengthCell()
!===================================================================================================================================
! use the number of electron density and electron temperature to compute the cold Debye-length per cell
! CAUTION: use SI-units
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars  ,ONLY:ElectronDensityCell,ElectronTemperatureCell,DebyeLengthCell,QuasiNeutralityCell
USE MOD_Globals_Vars           ,ONLY:ElementaryCharge, BoltzmannConst
USE MOD_Globals_Vars           ,ONLY:eps0
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iElem
!===================================================================================================================================

! nullify
DebyeLengthCell=0.

! loop over all elements and compute the plasma frequency with the use of the electron density
DO iElem=1,PP_nElems
  IF(ElectronDensityCell(iElem).LE.0.0) CYCLE ! ignore cells in which no electrons are present
  IF(QuasiNeutralityCell(iElem).LE.0.0) CYCLE ! ignore cells in which quasi neutrality is not possible
  DebyeLengthCell(iElem) = SQRT( (eps0*BoltzmannConst*ElectronTemperatureCell(iElem))/&
                                 (ElectronDensityCell(iElem)*(ElementaryCharge**2))       )
END DO ! iElem=1,PP_nElems

END SUBROUTINE CalculateDebyeLengthCell


SUBROUTINE CalculatePPDCell()
!===================================================================================================================================
! Calculate the points per Debye length for each cell
! PointsPerDebyeLength: PPD = (p+1)*lambda_D/L_cell
! where L_cell=V^(1/3) is the characteristic cell length determined from the cell volume
! PPDCellX, PPDCellY and PPDCellZ are determined by the average distance in X, Y and Z of each cell
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars ,ONLY: DebyeLengthCell,PPDCell,PPDCellX,PPDCellY,PPDCellZ,PICValidPlasmaCell,PPDCellResolved
USE MOD_Particle_Mesh_Vars    ,ONLY: ElemCharLength_Shared,ElemCharLengthX_Shared,ElemCharLengthY_Shared,ElemCharLengthZ_Shared
USE MOD_Mesh_Vars             ,ONLY: offSetElem
USE MOD_Mesh_Tools            ,ONLY: GetCNElemID
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iElem
!===================================================================================================================================
! Initialize number of elements with PPD>=1
PPDCellResolved = 0
! loop over all elements
DO iElem=1,PP_nElems
  ASSOCIATE( a        => (REAL(PP_N)+1.0)*DebyeLengthCell(iElem) ,&
             CNElemID => GetCNElemID(iElem+offSetElem)            )
    PPDCell(iElem)  = a/ElemCharLength_Shared( CNElemID) ! Determined with characteristic cell length
    PPDCellX(iElem) = a/ElemCharLengthX_Shared(CNElemID) ! Determined from average distance in X
    PPDCellY(iElem) = a/ElemCharLengthY_Shared(CNElemID) ! Determined from average distance in Y
    PPDCellZ(iElem) = a/ElemCharLengthZ_Shared(CNElemID) ! Determined from average distance in Z
    ! Check PIC quality factors
    IF((PICValidPlasmaCell(iElem).EQ.1))THEN
      IF(PPDCell(iElem) .GE.1.0) PPDCellResolved(1) = PPDCellResolved(1) + 1
      IF(PPDCellX(iElem).GE.1.0) PPDCellResolved(2) = PPDCellResolved(2) + 1
      IF(PPDCellY(iElem).GE.1.0) PPDCellResolved(3) = PPDCellResolved(3) + 1
      IF(PPDCellZ(iElem).GE.1.0) PPDCellResolved(4) = PPDCellResolved(4) + 1
    END IF ! PICValidPlasmaCell(iElem)
  END ASSOCIATE
END DO ! iElem=1,PP_nElems

END SUBROUTINE CalculatePPDCell


SUBROUTINE CalculatePICCFL()
!===================================================================================================================================
! Calculate the particle displacement CFL condition for PIC schemes
! PICCFL Condition: PICCFLCell = (p+1)*dt/L_cell * SQRT( kB*Te/me ) / 0.4  <  1.0
! where L_cell=V^(1/3) is the characteristic cell length determined from the cell volume
! PICCFLCellX, PICCFLCellY, PICCFLCellZ are determined by the average distance in X, Y and Z of each cell
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Preproc
USE MOD_Mesh_Vars             ,ONLY: offSetElem
USE MOD_Mesh_Tools            ,ONLY: GetCNElemID
USE MOD_Particle_Analyze_Vars ,ONLY: ElectronTemperatureCell,PICCFLCell,PICCFLCellX,PICCFLCellY,PICCFLCellZ
USE MOD_TimeDisc_Vars         ,ONLY: dt
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst,ElectronMass
USE MOD_Particle_Mesh_Vars    ,ONLY: ElemCharLength_Shared,ElemCharLengthX_Shared,ElemCharLengthY_Shared,ElemCharLengthZ_Shared
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iElem
!===================================================================================================================================
! loop over all elements
DO iElem=1,PP_nElems
  ! Divide the parameter by 0.4 -> the resulting value must always be below 1.0
  ASSOCIATE( a => (REAL(PP_N)+1.0) * 2.5 * dt * SQRT( BoltzmannConst * ElectronTemperatureCell(iElem) / ElectronMass ), &
             CNElemID => GetCNElemID(iElem+offSetElem) &
             )
    PICCFLCell(iElem)  = a/ElemCharLength_Shared(CNElemID) ! determined with characteristic cell length
    PICCFLCellX(iElem) = a/ElemCharLengthX_Shared(CNElemID) ! determined from average distance in X
    PICCFLCellY(iElem) = a/ElemCharLengthY_Shared(CNElemID) ! determined from average distance in Y
    PICCFLCellZ(iElem) = a/ElemCharLengthZ_Shared(CNElemID) ! determined from average distance in Z
  END ASSOCIATE
END DO ! iElem=1,PP_nElems

END SUBROUTINE CalculatePICCFL


SUBROUTINE CalculateMaxPartDisplacement()
!===================================================================================================================================
! Compute the maximum displacement of the fastest particle in each cell
! MaxPartDisplacement = max(v_iPart)*dT/L_cell <  1.0
! where L_cell=V^(1/3) is the characteristic cell length determined from the cell volume
! MaxPartDisplacementX, MaxPartDisplacementY, MaxPartDisplacementZ are determined by the average distance in X, Y and Z of each cell
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals               ,ONLY: VECNORM
USE MOD_Preproc
USE MOD_Mesh_Vars             ,ONLY: nElems, offSetElem
USE MOD_Mesh_Tools            ,ONLY: GetCNElemID
USE MOD_Particle_Analyze_Vars ,ONLY: MaxPartDisplacementCell
USE MOD_Particle_Analyze_Vars ,ONLY: MaxPartDisplacementCellX,MaxPartDisplacementCellY,MaxPartDisplacementCellZ
USE MOD_Particle_Vars         ,ONLY: PDM,PEM,PartState
USE MOD_TimeDisc_Vars         ,ONLY: dt
USE MOD_Particle_Mesh_Vars    ,ONLY: ElemCharLength_Shared,ElemCharLengthX_Shared,ElemCharLengthY_Shared,ElemCharLengthZ_Shared
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iElem,iPart
REAL                 :: MaxVelo(1:nElems,1:3)    ! each component
REAL                 :: MaxVeloAbs(1:nElems,1:3) ! fastest particle in 3D
!===================================================================================================================================
MaxVelo(1:nElems,1:3) = 0.0
MaxVeloAbs(1:nElems,1:3) = 0.0
! loop over all particles
DO iPart = 1, PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart)) THEN
    iElem = PEM%LocalElemID(iPart)
    ! Check velocity of each particle in each direction at get the highest value
    MaxVelo(iElem,1) = MAX(MaxVelo(iElem,1),PartState(4,iPart))
    MaxVelo(iElem,2) = MAX(MaxVelo(iElem,2),PartState(5,iPart))
    MaxVelo(iElem,3) = MAX(MaxVelo(iElem,3),PartState(6,iPart))
    ! Check for fastest particle in cell
    IF(VECNORM(PartState(4:6,iPart)).GT.VECNORM(MaxVeloAbs(iElem,1:3)))THEN
      MaxVeloAbs(iElem,1:3) = PartState(4:6,iPart)
    END IF
  END IF
END DO ! iPart = 1, PDM%ParticleVecLength

! loop over all elements
DO iElem=1,PP_nElems
  ! The resulting value must always be below 1.0
  ASSOCIATE( vAbs => VECNORM(MaxVeloAbs(iElem,1:3)) ,&
             vX   => MaxVelo(iElem,1)               ,&
             vY   => MaxVelo(iElem,2)               ,&
             vZ   => MaxVelo(iElem,3)               ,&
             a    => dt*(REAL(PP_N)+1.0)            ,&
             CNElemID => GetCNElemID(iElem+offSetElem)&
             )
    MaxPartDisplacementCell(iElem)  = a*vAbs/ElemCharLength_Shared(CNElemID)  ! determined with characteristic cell length
    MaxPartDisplacementCellX(iElem) = a*vX  /ElemCharLengthX_Shared(CNElemID)  ! determined from average distance in X
    MaxPartDisplacementCellY(iElem) = a*vY  /ElemCharLengthY_Shared(CNElemID)  ! determined from average distance in Y
    MaxPartDisplacementCellZ(iElem) = a*vZ  /ElemCharLengthZ_Shared(CNElemID)  ! determined from average distance in Z
  END ASSOCIATE
END DO ! iElem=1,PP_nElems

END SUBROUTINE CalculateMaxPartDisplacement


SUBROUTINE CalculateIonizationCell()
!===================================================================================================================================
! 1.) Count the number of ions per DG cell and divide it by element-volume -> ion density n_i
! 2.) Count the number of neutrals per DG cell and divide it by element-volume -> neutral density n_n
! 3.) Calculate the ionization degree: alpha = n_i/(n_i + n_n)
! 4.) Set indicator for valid plasma element: quasi-neutrality is above 0.5 and at least 20 particles are inside the element
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars ,ONLY: IonizationCell,QuasiNeutralityCell,NeutralDensityCell,ElectronDensityCell,IonDensityCell
USE MOD_Particle_Analyze_Vars ,ONLY: ChargeNumberCell,PICValidPlasmaCell,PICValidPlasmaCellSum
USE MOD_Particle_Mesh_Vars    ,ONLY: ElemVolume_Shared
USE MOD_Mesh_Vars             ,ONLY: offSetElem
USE MOD_Mesh_Tools            ,ONLY: GetCNElemID
USE MOD_LoadBalance_Vars      ,ONLY: nPartsPerElem
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iElem
!===================================================================================================================================
! Nullify
IonizationCell        = 0.
QuasiNeutralityCell   = 0.
PICValidPlasmaCellSum = 0

! Loop over all elements
DO iElem=1,PP_nElems
  ASSOCIATE(&
    Q   => QuasiNeutralityCell(iElem) ,& ! Quasi neutral condition approximation
    n_e => ElectronDensityCell(iElem) ,& ! Electron number density (cell average)
    X   => IonizationCell(iElem)      ,& ! Ionization degree (cell average)
    n_i => IonDensityCell(iElem)      ,& ! Ion number density (cell average)
    n_n => NeutralDensityCell(iElem)   ) ! Neutral number density (cell average)

    IF(ABS(n_i + n_n).LE.0.0)THEN ! no particles in cell
      X = 0.0
      Q = 0.0
    ELSE
      ! 0.  Set degree of ionization: X = n_i / (n_i + n_n)
      X  = n_i / (n_i + n_n)

      ! Set quasi neutrality between zero and unity depending on which density is larger
      ! Quasi neutrality holds, when n_e ~ Z_i*n_i (electron density approximately equal to ion density multiplied with charge number)
      ! 1.  Calculate Z_i*n_i (Charge density cell average)
      Q = ChargeNumberCell(iElem) / ElemVolume_Shared(GetCNElemID(iElem+offSetElem))

      ! 2.  Calculate the quasi neutrality parameter: should be near to 1 for quasi-neutrality
      IF(Q.GT.n_e)THEN
        ! 2.1  if Z_i*n_i > n_e -> calculate n_e/Z_i*n_i
        Q = n_e / Q
      ELSE
        ! 2.2  if Z_i*n_i < n_e -> calculate Z_i*n_i/n_e
        IF(ABS(n_e).GT.0.0)THEN
          Q = Q / n_e
        ELSE
          Q = 0.0
        END IF
      END IF
    END IF

    ! 4.) Set indicator for valid plasma element: quasi-neutrality is above 0.5 and at least 20 particles are inside the element
    ! Set PPD and time step validity flag for elements (otherwise the element is considered not part of a plasma state)
    ! Check that quasi-neutrality is above 0.5 and at least 20 particles are inside the element
    IF((Q.GE.0.5).AND.(nPartsPerElem(iElem).GE.20))THEN
      PICValidPlasmaCell(iElem) = 1
      PICValidPlasmaCellSum = PICValidPlasmaCellSum + 1
    ELSE
      PICValidPlasmaCell(iElem) = 0
    END IF ! Q.GE.0.5
  END ASSOCIATE
END DO ! iElem=1,PP_nElems

END SUBROUTINE CalculateIonizationCell


SUBROUTINE CalculatePlasmaParameter()
!===================================================================================================================================
! Calculate the Plasma parameter (here Debye number) for each cell
! Debye number: N_D = 4.0/3.0 * pi * n_e * lambda_D**3
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals_Vars          ,ONLY: PI
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars ,ONLY: DebyeLengthCell,ElectronDensityCell,PlasmaParameterCell
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iElem
!===================================================================================================================================
! loop over all elements
DO iElem=1,PP_nElems
  IF((DebyeLengthCell(iElem).GT.0.0).AND.(ElectronDensityCell(iElem).GT.0.0))THEN
    PlasmaParameterCell(iElem) = (4.0/3.0) * PI * ElectronDensityCell(iElem) * (DebyeLengthCell(iElem)**3)
  ELSE
    PlasmaParameterCell(iElem) = 0.0
  END IF
END DO ! iElem=1,PP_nElems

END SUBROUTINE CalculatePlasmaParameter


!===================================================================================================================================
!> Determines the kinetic energy of a (charged) particle before and after the push, the difference is stored as the coupled power
!===================================================================================================================================
SUBROUTINE CalcCoupledPowerPart(iPart,mode)
! MODULES
USE MOD_Particle_Vars           ,ONLY: PartSpecies, PEM
USE MOD_Particle_Analyze_Vars   ,ONLY: PCoupl, PCouplAverage, PCouplSpec, EDiff
USE MOD_Part_Tools              ,ONLY: isChargedParticle
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemVolume_Shared
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
USE MOD_Mesh_Vars               ,ONLY: offSetElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)              :: iPart                        !< Particle index
CHARACTER(LEN=*),INTENT(IN)     :: mode                         !< Mode: 'before' or 'after' the particle push
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iElem, iSpec
!===================================================================================================================================

IF(.NOT.isChargedParticle(iPart)) RETURN

SELECT CASE(TRIM(mode))
CASE('before')
  ! Kinetic energy before particle push (negative)
  EDiff = (-1.) * CalcEkinPart(iPart)
CASE('after')
  ! Kinetic energy after particle push (positive)
  EDiff         = EDiff + CalcEkinPart(iPart)
  PCoupl        = PCoupl + EDiff
  PCouplAverage = PCouplAverage + EDiff
  iElem         = PEM%LocalElemID(iPart)
  iSpec         = PartSpecies(iPart)
  PCouplSpec(iSpec)%DensityAvgElem(iElem) = PCouplSpec(iSpec)%DensityAvgElem(iElem) &
    + EDiff/ElemVolume_Shared(GetCNElemID(iElem+offSetElem))
END SELECT

END SUBROUTINE CalcCoupledPowerPart


#if USE_HDG
!===================================================================================================================================
!> Determine the electric potential for Dirichlet BCs in combination with a specific power input through these BCs
!===================================================================================================================================
SUBROUTINE CalculatePCouplElectricPotential()
! MODULES
USE MOD_Globals               ,ONLY: MPIRoot
USE MOD_HDG_Vars              ,ONLY: CoupledPowerPotential,CoupledPowerTarget,CoupledPowerRelaxFac,CoupledPowerFrequency
USE MOD_HDG_Vars              ,ONLY: CoupledPowerMode
USE MOD_Particle_Analyze_Vars ,ONLY: PCoupl,PCouplAverage,PCouplIntAverage
USE MOD_Restart_Vars          ,ONLY: RestartTime
USE MOD_TimeDisc_Vars         ,ONLY: Time
#if USE_MPI
USE MOD_Equation_Tools        ,ONLY: SynchronizeCPP
#endif /*USE_MPI*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: PowerRatio
!===================================================================================================================================
IF(MPIRoot)THEN
  ASSOCIATE(&
        Vmin => CoupledPowerPotential(1),&
           V => CoupledPowerPotential(2),&
        Vmax => CoupledPowerPotential(3),&
        VInt => CoupledPowerPotential(4),&
        tCPP => CoupledPowerPotential(5) &
        )

    PowerRatio = 1.0 ! Default
    SELECT CASE(CoupledPowerMode)
    CASE(1)
      ! Adjust electric potential depending on the instantaneous coupled power
      PowerRatio = CoupledPowerTarget / PCoupl
    CASE(2)
      ! Use moving average power
      IF(ABS(Time-RestartTime).GT.0.0) PowerRatio = CoupledPowerTarget / PCouplAverage
    CASE(3)
      ! Use integrated power (via user-defined frequency)
      IF(time.GE.tCPP)THEN ! Simulation time threshold reached
        ! Update time
        IF(CoupledPowerFrequency.GT.0.0) tCPP = tCPP + 1.0 / CoupledPowerFrequency
        ! Update Voltage
        PowerRatio = CoupledPowerTarget / PCouplIntAverage ! PCouplIntAverage is the running integrated average
        ! Reset integrated power value
        VInt = 0.
      END IF ! time.GE.tBV
    END SELECT

    ! Relaxation factor
    V = V * (1.0 + CoupledPowerRelaxFac * (PowerRatio - 1.0))

    ! Keep boundaries
    IF(V.GT.Vmax) V = Vmax
    IF(V.LT.Vmin) V = Vmin
  END ASSOCIATE
END IF ! MPIRoot

#if USE_MPI
! Synchronize CoupledPowerPotential from MPIRoot to all processors on the sub-communicator
CALL SynchronizeCPP()
#endif /*USE_MPI*/

END SUBROUTINE CalculatePCouplElectricPotential
#endif /*USE_HDG*/

#endif /*PARTICLES*/
END MODULE MOD_Particle_Analyze_Tools
