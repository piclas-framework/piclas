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

PUBLIC :: CalcEkinPart,CalcEkinPart2
PUBLIC :: CalcNumPartsOfSpec
PUBLIC :: AllocateElectronIonDensityCell,AllocateElectronTemperatureCell
PUBLIC :: CalculateElectronIonDensityCell,CalculateElectronTemperatureCell
!===================================================================================================================================

CONTAINS


PURE REAL FUNCTION CalcEkinPart(iPart)
!===================================================================================================================================
! computes the kinetic energy of one particle
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Globals_Vars  ,ONLY: c2, c2_inv
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
  IF (partV2.LT.1E12)THEN
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


PURE REAL FUNCTION CalcEkinPart2(velocity,Species_IN,WeightingFactor)
!===================================================================================================================================
! computes the kinetic energy of one particle given its velocity, species and weighting factor
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Globals_Vars  ,ONLY: c2,c2_inv
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
  IF (partV2.LT.1E12)THEN
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
#if USE_MPI
USE MOD_Particle_MPI_Vars     ,ONLY: PartMPI
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(OUT)                   :: NumSpec(nSpecAnalyze)
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
IF (PartMPI%MPIRoot) THEN
  IF(CalcNumSpec_IN)    CALL MPI_REDUCE(MPI_IN_PLACE,NumSpec    ,nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  IF(CalcSimNumSpec_IN) CALL MPI_REDUCE(MPI_IN_PLACE,SimNumSpec ,nSpecAnalyze,MPI_INTEGER_INT_KIND,MPI_SUM,0,PartMPI%COMM,IERROR)
ELSE
  IF(CalcNumSpec_IN)    CALL MPI_REDUCE(NumSpec     ,NumSpec    ,nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  IF(CalcSimNumSpec_IN) CALL MPI_REDUCE(SimNumSpec  ,SimNumSpec ,nSpecAnalyze,MPI_INTEGER_INT_KIND,MPI_SUM,0,PartMPI%COMM,IERROR)
END IF
#endif /*USE_MPI*/

! Set global number of particles (info for std out)
IF(CalcSimNumSpec_IN)THEN
  GlobalNbrOfParticlesUpdated = .TRUE.
#if USE_MPI
  IF(PartMPI%MPIRoot)THEN
#endif /*USE_MPI*/
    nGlobalNbrOfParticles = INT(SimNumSpec(nSpecAnalyze),KIND=IK)
#if USE_MPI
  END IF ! PartMPI%MPIRoot
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
USE MOD_Preproc               ,ONLY: PP_nElems
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
USE MOD_Preproc               ,ONLY: PP_nElems
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


SUBROUTINE CalculateElectronIonDensityCell()
!===================================================================================================================================
! Count the number of electrons per DG cell and divide it by element-volume
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
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
INTEGER              :: iPart,iElem,RegionID
REAL                 :: charge, MPF
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


SUBROUTINE CalculateElectronTemperatureCell()
!===================================================================================================================================
! Count the number of electrons per DG cell and divide it by element-volume
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals               ,ONLY: PARTISELECTRON
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst,ElectronMass,ElementaryCharge
USE MOD_Preproc
USE MOD_Particle_Analyze_Vars ,ONLY: ElectronTemperatureCell
USE MOD_Particle_Vars         ,ONLY: PDM,PEM,usevMPF,Species,PartSpecies,PartState
USE MOD_DSMC_Vars             ,ONLY: RadialWeighting
#if USE_HDG
USE MOD_HDG_Vars              ,ONLY: ElemToBRRegion,UseBRElectronFluid,RegionElectronRef
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
INTEGER :: iPart,iElem,ElemID,Method,RegionID
REAL    :: nElectronsPerCell(1:PP_nElems)
REAL    ::  PartVandV2(1:PP_nElems,1:6)
REAL    :: Mean_PartV2(1:3)
REAL    :: MeanPartV_2(1:3)
REAL    ::   TempDirec(1:3)
REAL    :: WeightingFactor
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
    ElemID                      = PEM%LocalElemID(iPart)
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
#endif /*PARTICLES*/


END MODULE MOD_Particle_Analyze_Tools
