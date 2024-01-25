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

MODULE MOD_DSMC_BGGas
!===================================================================================================================================
! Module for use of a background gas for the simulation of trace species (if number density of bg gas is multiple orders of
! magnitude larger than the trace species)
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
PUBLIC :: DefineParametersBGG, BGGas_Initialize, BGGas_InsertParticles, DSMC_pairing_bggas, BGGas_DeleteParticles
PUBLIC :: BGGas_AssignParticleProperties, BGGas_PhotoIonization, BGGas_InitRegions, BGGas_RegionsSetInternalTemp
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for background gas
!==================================================================================================================================
SUBROUTINE DefineParametersBGG()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Background gas")
CALL prms%CreateLogicalOption(  'Particles-BGGas-UseDistribution', &
                                      'Utilization of a cell-local background gas distribution as read-in from a previous '//&
                                      'DSMC/BGK result using Particles-MacroscopicRestart', '.FALSE.')
! Backgroun gas regions
CALL prms%CreateIntOption(      'Particles-BGGas-nRegions'                    ,'Number of background gas regions', '0')
CALL prms%CreateStringOption(   'Particles-BGGas-Region[$]-Type'              ,'Keyword for particle space condition of species [$] in case of multiple inits' , 'cylinder', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Particles-BGGas-Region[$]-RadiusIC'          ,'Outer radius'                 , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Particles-BGGas-Region[$]-Radius2IC'         ,'Inner radius (e.g. for a ring)' , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Particles-BGGas-Region[$]-BasePointIC'       , 'Base point'         , numberedmulti=.TRUE., no=3)
CALL prms%CreateRealArrayOption('Particles-BGGas-Region[$]-BaseVector1IC'     , 'First base vector'  , numberedmulti=.TRUE., no=3)
CALL prms%CreateRealArrayOption('Particles-BGGas-Region[$]-BaseVector2IC'     , 'Second base vector' , numberedmulti=.TRUE., no=3)
CALL prms%CreateRealOption(     'Particles-BGGas-Region[$]-CylinderHeightIC'  ,'Third measure of cylinder', numberedmulti=.TRUE.)
END SUBROUTINE DefineParametersBGG


SUBROUTINE BGGas_Initialize()
!===================================================================================================================================
!> Initialization of the background gas: compatibility check, array allocation, background species to species mapping and
!> calculation of the molar fraction
!===================================================================================================================================
! MODULES
USE MOD_ReadInTools
USE MOD_Globals               ,ONLY: abort
USE MOD_DSMC_Vars             ,ONLY: BGGas
USE MOD_Mesh_Vars             ,ONLY: nElems
USE MOD_Particle_Vars         ,ONLY: PDM, Species, nSpecies, UseVarTimeStep, VarTimeStep
USE MOD_Restart_Vars          ,ONLY: DoMacroscopicRestart, MacroRestartFileName
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iSpec, bgSpec, iElem
REAL              :: SpeciesDensTmp(1:nSpecies)
!===================================================================================================================================

! 0.) Variable read-in
IF(BGGas%UseDistribution) MacroRestartFileName = GETSTR('Particles-MacroscopicRestart-Filename')

! 1.) Check compatibility with other features and whether required parameters have been read-in
IF(UseVarTimeStep) THEN
  IF(.NOT.VarTimeStep%UseSpeciesSpecific) CALL abort(__STAMP__, &
    'ERROR: 2D/Axisymmetric and variable timestep (except species-specific) are not implemented with a background gas yet!')
END IF

DO iSpec = 1, nSpecies
  IF(BGGas%BackgroundSpecies(iSpec)) THEN
    IF(Species(iSpec)%NumberOfInits.NE.1.AND..NOT.BGGas%UseRegions) CALL abort(__STAMP__, 'BGG species can be used ONLY for BGG!')
    IF(.NOT.BGGas%UseDistribution.AND..NOT.BGGas%UseRegions)THEN
      IF(BGGas%NumberDensity(iSpec).EQ.0.) CALL abort(__STAMP__, 'NumberDensity is zero but must be defined for a background gas!')
    END IF ! .NOT.BGGas%UseDistribution
  END IF
END DO

IF(DoMacroscopicRestart) CALL abort(__STAMP__, 'Constant background gas and macroscopic restart are not compatible!')

! 2.) Allocation
IF(BGGas%UseDistribution.OR.BGGas%UseRegions) THEN
  ALLOCATE(BGGas%Distribution(1:BGGas%NumberOfSpecies,1:10,1:nElems))
  BGGas%Distribution = 0.
  ALLOCATE(BGGas%SpeciesFractionElem(1:BGGas%NumberOfSpecies,1:nElems))
  BGGas%SpeciesFractionElem = 0.
ELSE
  ! Backup densities of all background gas species
  SpeciesDensTmp(1:nSpecies) = BGGas%NumberDensity(1:nSpecies)
  DEALLOCATE(BGGas%NumberDensity)
  ! Re-allocate with the correct size
  ALLOCATE(BGGas%NumberDensity(BGGas%NumberOfSpecies))
  BGGas%NumberDensity = 0.
  ALLOCATE(BGGas%SpeciesFraction(BGGas%NumberOfSpecies))
  BGGas%SpeciesFraction = 0.
END IF

ALLOCATE(BGGas%PairingPartner(PDM%maxParticleNumber))
BGGas%PairingPartner = 0
ALLOCATE(BGGas%MapSpecToBGSpec(nSpecies))
BGGas%MapSpecToBGSpec = 0

ALLOCATE(BGGas%MapBGSpecToSpec(BGGas%NumberOfSpecies))
BGGas%MapBGSpecToSpec = 0
BGGas%MaxMPF = 0.

! 3.) Create a mapping of background species to regular species and vice versa, calculate the molar fraction
bgSpec = 0
DO iSpec = 1, nSpecies
  IF(BGGas%BackgroundSpecies(iSpec)) THEN
    bgSpec = bgSpec + 1
    IF(bgSpec.GT.BGGas%NumberOfSpecies) CALL Abort(__STAMP__,'More background species detected than previously defined!')
    IF((.NOT.BGGas%UseDistribution).AND.(.NOT.BGGas%UseRegions)) BGGas%NumberDensity(bgSpec) = SpeciesDensTmp(iSpec)
    BGGas%MapSpecToBGSpec(iSpec)  = bgSpec
    BGGas%MapBGSpecToSpec(bgSpec) = iSpec
    BGGas%MaxMPF = MAX(BGGas%MaxMPF,Species(iSpec)%MacroParticleFactor)
  END IF
END DO

! 4.) Read-in a background gas distribution
IF(BGGas%UseDistribution) CALL BGGas_ReadInDistribution()

! 5.) Determine species fraction (not yet for background gas regions)
DO bgSpec = 1, BGGas%NumberOfSpecies
  IF(BGGas%UseDistribution) THEN
    DO iElem = 1, nElems
      IF(SUM(BGGas%Distribution(:,7,iElem)).GT.0.)THEN
        BGGas%SpeciesFractionElem(bgSpec,iElem) = BGGas%Distribution(bgSpec,7,iElem) / SUM(BGGas%Distribution(:,7,iElem))
      END IF ! SUM(BGGas%Distribution(:,7,iElem)).GT.0.
    END DO ! iElem = 1, nElems
  ELSEIF(.NOT.BGGas%UseRegions) THEN
    IF(SUM(SpeciesDensTmp).GT.0.)THEN
      BGGas%SpeciesFraction(bgSpec) = BGGas%NumberDensity(bgSpec) / SUM(SpeciesDensTmp)
    END IF ! SUM(SpeciesDensTmp).GT.0.
  END IF
END DO ! bgSpec = 1, BGGas%NumberOfSpecies

END SUBROUTINE BGGas_Initialize


INTEGER FUNCTION BGGas_GetSpecies(iElem)
!===================================================================================================================================
!> Get a species index of the background gas by randomly choosing a species based on the molar fraction
!===================================================================================================================================
! MODULES
USE MOD_Globals               ,ONLY: Abort
USE MOD_DSMC_Vars             ,ONLY: BGGas
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: iRan,FractionSum
INTEGER            :: iSpec
!===================================================================================================================================

BGGas_GetSpecies = 0

IF(BGGas%NumberOfSpecies.GT.1) THEN
  CALL RANDOM_NUMBER(iRan)
  DO iSpec = 1, BGGas%NumberOfSpecies
    ! Add up fractions from 1 to iSpec
    IF(BGGas%UseDistribution)THEN
      FractionSum = SUM(BGGas%SpeciesFractionElem(1:iSpec,iElem))
    ELSE
      FractionSum = SUM(BGGas%SpeciesFraction(1:iSpec))
    END IF ! BGGas%UseDistribution
    ! Check if sum of fractions is met
    IF(FractionSum.GT.iRan) THEN
      BGGas_GetSpecies = BGGas%MapBGSpecToSpec(iSpec)
      RETURN
    END IF
  END DO
ELSE
  BGGas_GetSpecies = BGGas%MapBGSpecToSpec(1)
END IF

IF(BGGas_GetSpecies.EQ.0) CALL Abort(__STAMP__,'ERROR in BGGas: Background gas species is not set correctly!')

END FUNCTION BGGas_GetSpecies


SUBROUTINE BGGas_InsertParticles()
!===================================================================================================================================
!> Creating particles of the background species for each actual simulation particle
!> 1. Initialize particles (loop over the current ParticleVecLength):
!>    a) Get the new index from the nextFreePosition array
!>    b) Same position as the non-BGG particle, set species and initialize the internal energy (if required)
!>    c) Include BGG-particles in PEM%p-Lists (counting towards the amount of particles per cell stored in PEM%pNumber)
!>    d) Map BGG-particle to the non-BGG particle for particle pairing
!> 2. Call SetParticleVelocity: loop over the newly created particles to set the thermal velocity, using nextFreePosition to get
!>    the same indices, possible since the CurrentNextFreePosition was not updated yet
!> 3. Adjust ParticleVecLength and currentNextFreePosition
!===================================================================================================================================
! MODULES
USE MOD_Globals                ,ONLY: Abort
USE MOD_DSMC_Vars              ,ONLY: BGGas
USE MOD_PARTICLE_Vars          ,ONLY: PDM, PartSpecies, PEM
USE MOD_Part_Tools             ,ONLY: GetNextFreePosition
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers      ,ONLY: LBStartTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iNewPart, iPart, PositionNbr, iSpec, LocalElemID
#if USE_LOADBALANCE
REAL              :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

iNewPart=0
PositionNbr = 0
DO iPart = 1, PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    ! Skip background particles that have been created within this loop
    IF(BGGas%BackgroundSpecies(PartSpecies(iPart))) CYCLE
    ! Skip particles outside of any regions
    IF(BGGas%UseRegions) THEN
      IF(BGGas%RegionElemType(PEM%LocalElemID(iPart)).EQ.0) CYCLE
    END IF
    ! Get a free particle index
    iNewPart = iNewPart + 1
    PositionNbr = GetNextFreePosition()
    ! Get the background gas species
    iSpec = BGGas_GetSpecies(PEM%LocalElemID(iPart))
    ! Assign particle properties
    CALL BGGas_AssignParticleProperties(iSpec,iPart,PositionNbr)
    ! Set the pairing index
    BGGas%PairingPartner(iPart) = PositionNbr
    ! Update cell-local particle list
    LocalElemID = PEM%LocalElemID(PositionNbr)
    PEM%pNext(PEM%pEnd(LocalElemID)) = PositionNbr
    PEM%pEnd(LocalElemID) = PositionNbr
    PEM%pNumber(LocalElemID) = PEM%pNumber(LocalElemID) + 1
  END IF
END DO

#if USE_LOADBALANCE
CALL LBPauseTime(LB_DSMC,tLBStart)
#endif /*USE_LOADBALANCE*/
END SUBROUTINE BGGas_InsertParticles


!===================================================================================================================================
!> Assign properties for a single particle from the background gas using the species ID of the background gas, particle index and
!> an index for the background particle. Optional flags can be used to disable the calculation of the velocity and internal energy
!> state (e.g. if an old state is reused)
!===================================================================================================================================
SUBROUTINE BGGas_AssignParticleProperties(SpecID,PartIndex,bggPartIndex,GetVelocity_opt,GetInternalEnergy_opt)
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars           ,ONLY: PDM, PEM, PartState,PartSpecies,PartPosRef, usevMPF, PartMPF
USE MOD_Particle_Vars           ,ONLY: UseVarTimeStep, PartTimeStep
USE MOD_DSMC_Vars               ,ONLY: CollisMode, SpecDSMC, BGGas
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
USE MOD_part_emission_tools     ,ONLY: CalcVelocity_maxwell_lpn, DSMC_SetInternalEnr_LauxVFD
USE MOD_DSMC_PolyAtomicModel    ,ONLY: DSMC_SetInternalEnr_Poly
USE MOD_part_tools              ,ONLY: CalcVelocity_maxwell_particle
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER, INTENT(IN)             :: SpecID             !< Species ID
INTEGER, INTENT(IN)             :: PartIndex          !< ID of simulation particle
INTEGER, INTENT(IN)             :: bggPartIndex       !< ID of the newly created background gas particle
LOGICAL, INTENT(IN), OPTIONAL   :: GetVelocity_opt        !< Default: T, get a new velocity vector from the background gas properties
LOGICAL, INTENT(IN), OPTIONAL   :: GetInternalEnergy_opt  !< Default: T, get a new energy values from the background gas properties
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                         :: bggSpec, LocalElemID
LOGICAL                         :: GetVelocity, GetInternalEnergy
!===================================================================================================================================

IF(PRESENT(GetVelocity_opt)) THEN
  GetVelocity = GetVelocity_opt
ELSE
  GetVelocity = .TRUE.
END IF

IF(PRESENT(GetInternalEnergy_opt)) THEN
  GetInternalEnergy = GetInternalEnergy_opt
ELSE
  GetInternalEnergy = .TRUE.
END IF

! Global element index (Must be before internal energy: BGGas distribution requires the local element ID and uses the background gas particle index to get it)
PEM%GlobalElemID(bggPartIndex) = PEM%GlobalElemID(PartIndex)
PEM%LastGlobalElemID(bggPartIndex) = PEM%GlobalElemID(PartIndex)
LocalElemID = PEM%LocalElemID(PartIndex)
! Position
PartState(1:3,bggPartIndex) = PartState(1:3,PartIndex)
IF(TrackingMethod.EQ.REFMAPPING) PartPosRef(1:3,bggPartIndex)=PartPosRef(1:3,PartIndex)
! Species
PartSpecies(bggPartIndex) = SpecID
! Velocity
IF(GetVelocity) THEN
  IF(BGGas%UseDistribution) THEN
    bggSpec = BGGas%MapSpecToBGSpec(SpecID)
    PartState(4:6,bggPartIndex) = CalcVelocity_maxwell_particle(SpecID,BGGas%Distribution(bggSpec,4:6,LocalElemID)) &
                                  + BGGas%Distribution(bggSpec,1:3,LocalElemID)
  ELSE
    CALL CalcVelocity_maxwell_lpn(FractNbr=SpecID, Vec3D=PartState(4:6,bggPartIndex), iInit=1)
  END IF
END IF
! Internal energy
IF(CollisMode.GT.1) THEN
  IF(GetInternalEnergy) THEN
    IF(SpecDSMC(SpecID)%PolyatomicMol) THEN
      CALL DSMC_SetInternalEnr_Poly(SpecID,1,bggPartIndex,1)
    ELSE
      CALL DSMC_SetInternalEnr_LauxVFD(SpecID,1,bggPartIndex,1)
    END IF
  END IF
END IF
! Simulation flags
PDM%ParticleInside(bggPartIndex) = .TRUE.
PDM%IsNewPart(bggPartIndex)       = .TRUE.
PDM%dtFracPush(bggPartIndex)      = .FALSE.
! Weighting factors
IF(usevMPF) PartMPF(bggPartIndex) = PartMPF(PartIndex)
IF(UseVarTimeStep) PartTimeStep(bggPartIndex) = PartTimeStep(PartIndex)

END SUBROUTINE BGGas_AssignParticleProperties


SUBROUTINE DSMC_pairing_bggas(iElem)
!===================================================================================================================================
! Building of pairs for the background gas
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Analyze          ,ONLY: CalcGammaVib, CalcMeanFreePath
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_DSMC_Vars             ,ONLY: Coll_pData, CollInf, BGGas, CollisMode, ChemReac, PartStateIntEn, DSMC, SelectionProc
USE MOD_Particle_Vars         ,ONLY: PEM,PartSpecies,nSpecies,PartState,Species,usevMPF,Species, WriteMacroVolumeValues
USE MOD_Particle_Mesh_Vars    ,ONLY: ElemVolume_Shared
USE MOD_Mesh_Vars             ,ONLY: offsetElem
USE MOD_DSMC_Collis           ,ONLY: DSMC_perform_collision
USE MOD_DSMC_Relaxation       ,ONLY: FinalizeCalcVibRelaxProb, SumVibRelaxProb, InitCalcVibRelaxProb
USE MOD_TimeDisc_Vars         ,ONLY: TEnd, time
USE MOD_DSMC_CollisionProb    ,ONLY: DSMC_prob_calc
USE MOD_DSMC_Relaxation       ,ONLY: CalcMeanVibQuaDiatomic
USE MOD_Mesh_Tools            ,ONLY: GetCNElemID
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: nPair, iPair, iPart, iLoop, nPart, CNElemID
INTEGER                       :: cSpec1, cSpec2, iCase, iSpec, bggSpec
REAL                          :: iRan, MPF
!===================================================================================================================================

! Skip elements outside of any background gas regions
IF(BGGas%UseRegions) THEN
  IF(BGGas%RegionElemType(iElem).EQ.0) RETURN
END IF

nPart = PEM%pNumber(iElem)
nPair = INT(nPart/2.)

! Routine to increase the sample size for trace background gas species by splitting the simulation particle
IF(usevMPF.AND.ANY(BGGas%TraceSpecies(:))) CALL BGGas_TraceSpeciesSplit(iElem, nPart, nPair)

! Initialize variables
CALL InitCalcVibRelaxProb()

CollInf%Coll_SpecPartNum = 0.
CollInf%Coll_CaseNum = 0
ALLOCATE(Coll_pData(nPair))
Coll_pData%Ec=0

IF (CollisMode.EQ.3) ChemReac%MeanEVib_PerIter(1:nSpecies) = 0.0

! Pairing of particles
iPair = 1
iPart = PEM%pStart(iElem)
DO iLoop = 1, nPart
  iSpec = PartSpecies(iPart)
  ! Counting the number of particles per species
  MPF = GetParticleWeight(iPart)
  CollInf%Coll_SpecPartNum(iSpec) = CollInf%Coll_SpecPartNum(iSpec) + MPF
  ! Calculation of mean vibrational energy per cell and iter, necessary for dissociation probability
  IF (CollisMode.EQ.3) ChemReac%MeanEVib_PerIter(iSpec) = ChemReac%MeanEVib_PerIter(iSpec) + PartStateIntEn(1,iPart) * MPF
  ! Creating pairs for species, which are not the background species
  IF(.NOT.BGGas%BackgroundSpecies(iSpec)) THEN
    Coll_pData(iPair)%iPart_p1 = iPart
    Coll_pData(iPair)%iPart_p2 = BGGas%PairingPartner(iPart)
    iPair = iPair + 1
  END IF
  iPart = PEM%pNext(iPart)
END DO

IF(((CollisMode.GT.1).AND.(SelectionProc.EQ.2)).OR.DSMC%BackwardReacRate.OR.DSMC%CalcQualityFactors) THEN
  ! 1. Case: Inelastic collisions and chemical reactions with the Gimelshein relaxation procedure and variable vibrational
  !           relaxation probability (CalcGammaVib)
  ! 2. Case: Chemical reactions and backward rate require cell temperature for the partition function and equilibrium constant
  ! 3. Case: Temperature required for the mean free path with the VHS model
  ! Instead of calculating the translation temperature, simply the input value of the BG gas is taken. If the other species have
  ! an impact on the temperature, a background gas should not be utilized in the first place.
  DSMC%InstantTransTemp(nSpecies+1) = 0.
  DO iSpec = 1, nSpecies
    IF(BGGas%BackgroundSpecies(iSpec)) THEN
      bggSpec = BGGas%MapSpecToBGSpec(iSpec)
      IF(BGGas%UseDistribution) THEN
        DSMC%InstantTransTemp(nSpecies+1) = DSMC%InstantTransTemp(nSpecies+1) &
                                        + BGGas%SpeciesFractionElem(bggSpec,iElem) * SUM(BGGas%Distribution(bggSpec,4:6,iElem)) / 3.
      ELSE
        DSMC%InstantTransTemp(nSpecies+1) = DSMC%InstantTransTemp(nSpecies+1) &
                                        + BGGas%SpeciesFraction(bggSpec) * Species(iSpec)%Init(1)%MWTemperatureIC
      END IF
    END IF
  END DO
  IF(SelectionProc.EQ.2) CALL CalcGammaVib()
END IF

DO iSpec = 1, nSpecies
  IF(BGGas%BackgroundSpecies(iSpec)) THEN
    CNElemID = GetCNElemID(iElem+offSetElem)
    bggSpec   = BGGas%MapSpecToBGSpec(iSpec)
    IF(BGGas%UseDistribution)THEN
      CollInf%Coll_SpecPartNum(iSpec) = BGGas%Distribution(bggSpec,7,iElem)*ElemVolume_Shared(CNElemID)
    ELSE
      CollInf%Coll_SpecPartNum(iSpec) = BGGas%NumberDensity(bggSpec)       *ElemVolume_Shared(CNElemID)
    END IF ! BGGas%UseDistribution
    ! MPF is multiplied again in ReactionDecision()
    IF(.NOT.usevMPF) CollInf%Coll_SpecPartNum(iSpec) = CollInf%Coll_SpecPartNum(iSpec) / Species(iSpec)%MacroParticleFactor
  END IF
END DO

CollInf%SumPairMPF = 0.

DO iPair = 1, nPair
  cSpec1 = PartSpecies(Coll_pData(iPair)%iPart_p1) !spec of particle 1
  cSpec2 = PartSpecies(Coll_pData(iPair)%iPart_p2) !spec of particle 2
  iCase = CollInf%Coll_Case(cSpec1, cSpec2)
  CollInf%Coll_CaseNum(iCase) = CollInf%Coll_CaseNum(iCase) + 1 !sum of coll case (Sab)
  CollInf%SumPairMPF(iCase) = CollInf%SumPairMPF(iCase) + GetParticleWeight(Coll_pData(iPair)%iPart_p1)
  Coll_pData(iPair)%CRela2 = (PartState(4,Coll_pData(iPair)%iPart_p1) &
                            -  PartState(4,Coll_pData(iPair)%iPart_p2))**2 &
                            + (PartState(5,Coll_pData(iPair)%iPart_p1) &
                            -  PartState(5,Coll_pData(iPair)%iPart_p2))**2 &
                            + (PartState(6,Coll_pData(iPair)%iPart_p1) &
                            -  PartState(6,Coll_pData(iPair)%iPart_p2))**2
  Coll_pData(iPair)%PairType = iCase
  Coll_pData(iPair)%NeedForRec = .FALSE.
END DO

IF(CollisMode.EQ.3) CALL CalcMeanVibQuaDiatomic()

! 6.) Calculate the collision probability and perform the collision if necessary

DO iPair = 1, nPair
  CALL SumVibRelaxProb(iPair)
  IF(.NOT.Coll_pData(iPair)%NeedForRec) THEN
    CALL DSMC_prob_calc(iElem, iPair)
    CALL RANDOM_NUMBER(iRan)
    IF (Coll_pData(iPair)%Prob.ge.iRan) THEN
      CALL DSMC_perform_collision(iPair,iElem)
    END IF
  END IF
END DO
IF(DSMC%CalcQualityFactors) THEN
  IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
    ! Calculation of the mean free path
    DSMC%MeanFreePath = CalcMeanFreePath(REAL(CollInf%Coll_SpecPartNum),SUM(CollInf%Coll_SpecPartNum), &
                          ElemVolume_Shared(GetCNElemID(iElem+offSetElem)), DSMC%InstantTransTemp(nSpecies+1))
    ! Determination of the MCS/MFP for the case without octree
    IF((DSMC%CollSepCount.GT.0.0).AND.(DSMC%MeanFreePath.GT.0.0)) DSMC%MCSoverMFP = (DSMC%CollSepDist/DSMC%CollSepCount) &
                                                                                    / DSMC%MeanFreePath
  END IF
END IF
DEALLOCATE(Coll_pData)
CALL FinalizeCalcVibRelaxProb(iElem)

END SUBROUTINE DSMC_pairing_bggas


SUBROUTINE BGGas_DeleteParticles()
!===================================================================================================================================
! Deletes all background gas particles and updates the particle index list
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars           ,ONLY: BGGas
USE MOD_PARTICLE_Vars       ,ONLY: PDM, PartSpecies
USE MOD_part_tools          ,ONLY: UpdateNextFreePosition
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers  ,ONLY: LBStartTime, LBPauseTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iPart
#if USE_LOADBALANCE
REAL                        :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

DO iPart = 1, PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    ! No Pt_temp=0 necessary, because it is a ghost particle
    IF(BGGas%BackgroundSpecies(PartSpecies(iPart))) PDM%ParticleInside(iPart) = .FALSE.
  END IF
END DO
BGGas%PairingPartner = 0

#if USE_LOADBALANCE
CALL LBPauseTime(LB_DSMC,tLBStart)
#endif /*USE_LOADBALANCE*/

CALL UpdateNextFreePosition()

END SUBROUTINE BGGas_DeleteParticles


SUBROUTINE BGGas_PhotoIonization(iSpec,iInit,TotalNbrOfReactions)
!===================================================================================================================================
!> Particle emission through photo ionization, defined by chemical reactions of the type "phIon" and an emission with SpaceIC
!> "photon_cylinder" for the
!> 0.) Determine the total ionization cross-section
!> 1.) Compute the number of photoionization events in the local domain of each proc
!> 2.) Delete left-over inserted particles
!> 3.) Insert the products of the photoionization reactions
!> 4.) Perform the reaction, distribute the collision energy (including photon energy) and emit electrons perpendicular
!>     to the photon's path
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Analyze           ,ONLY: CalcGammaVib, CalcMeanFreePath
USE MOD_DSMC_Vars              ,ONLY: Coll_pData, CollisMode, ChemReac, PartStateIntEn, DSMC
USE MOD_DSMC_Vars              ,ONLY: SpecDSMC, DSMCSumOfFormedParticles
USE MOD_DSMC_Vars              ,ONLY: newAmbiParts, iPartIndx_NodeNewAmbi, BGGas
USE MOD_Particle_Vars          ,ONLY: PEM, PDM, PartSpecies, PartState, Species, usevMPF, PartMPF, Species, PartPosRef
USE MOD_part_emission_tools    ,ONLY: DSMC_SetInternalEnr_LauxVFD
USE MOD_DSMC_PolyAtomicModel   ,ONLY: DSMC_SetInternalEnr_Poly
USE MOD_part_pos_and_velo      ,ONLY: SetParticleVelocity
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_part_emission_tools    ,ONLY: CalcVelocity_maxwell_lpn
USE MOD_DSMC_ChemReact         ,ONLY: PhotoIonization_InsertProducts
USE MOD_DSMC_AmbipolarDiffusion,ONLY: AD_DeleteParticles
USE MOD_part_tools             ,ONLY: CalcVelocity_maxwell_particle
USE MOD_MCC_Vars               ,ONLY: PhotoIonFirstLine,PhotoIonLastLine,PhotoReacToReac,PhotonEnergies
USE MOD_MCC_Vars               ,ONLY: NbrOfPhotonXsecReactions,SpecPhotonXSecInterpolated,MaxPhotonXSec
USE MOD_Part_Tools             ,ONLY: GetNextFreePosition, IncreaseMaxParticleNumber
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec,iInit
INTEGER, INTENT(INOUT)        :: TotalNbrOfReactions
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart, iPair, iNewPart, iReac, ParticleIndex, NewParticleIndex, bgSpec, NbrOfParticle, LocalElemID
REAL                          :: RandVal,NumTmp,ProbRest
INTEGER                       :: TotalNbrOfReactionsTmp,iCrossSection,NbrCrossSections,iLine,iPhotoReac,iLineStart,iLineEnd,iBGGSpec
INTEGER                       :: NumPhotoIonization(ChemReac%NumOfReact)
REAL                          :: SumCrossSections,CrossSection,MaxCrossSection
!===================================================================================================================================
NumPhotoIonization = 0

IF(TotalNbrOfReactions.LE.0) RETURN
TotalNbrOfReactionsTmp = TotalNbrOfReactions

! Ambipolar diffusion (up to four new particles can be produced during a single photo-ionization event)
newAmbiParts = 0
IF (ALLOCATED(iPartIndx_NodeNewAmbi)) DEALLOCATE(iPartIndx_NodeNewAmbi)
ALLOCATE(iPartIndx_NodeNewAmbi(4*TotalNbrOfReactions))

!> 0.) Calculate sum of cross-sections for photoionization
SumCrossSections = 0.
NbrCrossSections = 0
DO iReac = 1, ChemReac%NumOfReact
  ! Only treat photoionization reactions
  SELECT CASE(ChemReac%ReactModel(iReac))
  CASE('phIon')
    SumCrossSections = SumCrossSections + ChemReac%CrossSection(iReac)
    NbrCrossSections = NbrCrossSections + 1
  CASE('phIonXsec')
    ! already calculated
  CASE DEFAULT
    CYCLE
  END SELECT
END DO ! iReac = 1, ChemReac%NumOfReact

!> 1.) Compute the number of photoionization events in the local domain of each proc
IF(NbrOfPhotonXsecReactions.GT.0)THEN
  iLineStart     = HUGE(1) ! Hugeify
  iLineEnd       = 0       ! Nullify
  PhotonEnergies = 0       ! Nullify
  ! Loop until all reactions have been matched with a specific wavelength
  DO WHILE(TotalNbrOfReactionsTmp.GT.0)
    ! 1.1) Select a wave length (or corresponding photon energy)
    ! Get 1st random number
    CALL RANDOM_NUMBER(RandVal)
    ! Get random wavelength (from line spectrum)
    iLine = INT(RandVal*REAL(PhotoIonLastLine-PhotoIonFirstLine+1)) + PhotoIonFirstLine

    ! Get 2nd random number
    CALL RANDOM_NUMBER(RandVal)
    ! Probe if the line is accepted by comparing against the energy fraction (maximum is 1.)
    IF(RandVal.GT.SpecPhotonXSecInterpolated(iLine,2)/MaxPhotonXSec) CYCLE
    ! Store photon energy for later chemical reaction
    PhotonEnergies(iLine,1) = PhotonEnergies(iLine,1) + 1

    ! 1.2) Select a cross-section
    PDF: DO
      ! Get 3rd random number
      CALL RANDOM_NUMBER(RandVal)
      ! Get random cross-section
      iPhotoReac = INT(RandVal*REAL(NbrOfPhotonXsecReactions) + 1.0)
      ! Get 4th random number
      CALL RANDOM_NUMBER(RandVal)
      ! Check if cross-section is > 0.
      CrossSection = SpecPhotonXSecInterpolated(iLine,2+iPhotoReac)
      IF(SpecPhotonXSecInterpolated(iLine,2+iPhotoReac).LE.0.) CYCLE PDF
      ! Probe if the line is accepted by comparing against the energy fraction (maximum is 1.)
      MaxCrossSection = MAXVAL(SpecPhotonXSecInterpolated(:,3:))
      IF(RandVal.LE.CrossSection/MaxCrossSection) EXIT PDF
    END DO PDF
    ! Store photon reaction for later chemical reaction
    PhotonEnergies(iLine,1+iPhotoReac) = PhotonEnergies(iLine,1+iPhotoReac) + 1
    iLineStart = MIN(iLineStart,iLine)
    iLineEnd   = MAX(iLineEnd,iLine)
    !WRITE (*,*) "iLineStart,iLineEnd =", iLineStart,iLineEnd

    ! 1.3) Reaction and line have been selected
    iReac = PhotoReacToReac(iPhotoReac)
    !IPWRITE(UNIT_StdOut,'(I6,10X,3(A,I3),A,E24.12)') " iLine =",iLine," iPhotoReac =",iPhotoReac," iReac =",iReac," CrossSection =",CrossSection
    NumPhotoIonization(iReac) = NumPhotoIonization(iReac) + 1
    TotalNbrOfReactionsTmp    = TotalNbrOfReactionsTmp - 1
  END DO ! WHILE(TotalNbrOfReactionsTmp.GT.0)
ELSE
  ! Photoionization with const. cross-section data
  iCrossSection = 0
  DO iReac = 1, ChemReac%NumOfReact
    ! Only treat photoionization reactions
    IF(.NOT.StringBeginsWith(ChemReac%ReactModel(iReac),'phIon')) CYCLE
    iCrossSection  = iCrossSection + 1
    IF(iCrossSection.EQ.NbrCrossSections)THEN
      NumPhotoIonization(iReac) = TotalNbrOfReactionsTmp
      EXIT
    END IF ! iCrossSection.EQ.NbrCrossSections
    NumTmp = TotalNbrOfReactionsTmp*ChemReac%CrossSection(iReac)/SumCrossSections
    SumCrossSections = SumCrossSections - ChemReac%CrossSection(iReac)
    NumPhotoIonization(iReac) = INT(NumTmp)
    ProbRest = NumTmp - REAL(NumPhotoIonization(iReac))
    CALL RANDOM_NUMBER(RandVal)
    IF (ProbRest.GT.RandVal) NumPhotoIonization(iReac) = NumPhotoIonization(iReac) + 1
    TotalNbrOfReactionsTmp = TotalNbrOfReactionsTmp - NumPhotoIonization(iReac)
  END DO
END IF ! NbrOfPhotonXsecReactions.GT.0

NbrOfParticle = SUM(NumPhotoIonization)

!> 2.) Delete left-over inserted particles
IF(TotalNbrOfReactions.GT.NbrOfParticle) THEN
  DO iPart = NbrOfParticle+1,TotalNbrOfReactions
    PDM%ParticleInside(GetNextFreePosition(iPart)) = .FALSE.
  END DO
  TotalNbrOfReactions = NbrOfParticle
ELSE IF(TotalNbrOfReactions.LT.NbrOfParticle) THEN
  CALL Abort(__STAMP__,'PhotoIonization: Something is wrong, trying to perform more reactions than anticipated!')
END IF

!> 3.) Insert the products of the photoionization reactions
IF(NbrOfParticle.EQ.0) RETURN

ALLOCATE(Coll_pData(NbrOfParticle))
Coll_pData%Ec=0.
DSMCSumOfFormedParticles = 0

iNewPart = 0; iPair = 0

DO iPart = 1, NbrOfParticle
  ! Loop over the particles with a set position (from SetParticlePosition)
  ParticleIndex = GetNextFreePosition(iPart)
  IF (DSMC%DoAmbipolarDiff) THEN
    newAmbiParts = newAmbiParts + 1
    iPartIndx_NodeNewAmbi(newAmbiParts) = ParticleIndex
  END IF
  iNewPart = iNewPart + 1
  ! Get a new index for the second product
  NewParticleIndex = GetNextFreePosition(iNewPart+NbrOfParticle)
  IF (NewParticleIndex.EQ.0) THEN
    CALL Abort(__STAMP__,'ERROR in PhotoIonization: MaxParticleNumber should be increased!')
  END IF
  IF (DSMC%DoAmbipolarDiff) THEN
    newAmbiParts = newAmbiParts + 1
    iPartIndx_NodeNewAmbi(newAmbiParts) = NewParticleIndex
  END IF
  PartState(1:3,NewParticleIndex) = PartState(1:3,ParticleIndex)
  IF(TrackingMethod.EQ.REFMAPPING)THEN ! here Nearst-GP is missing
    PartPosRef(1:3,NewParticleIndex)=PartPosRef(1:3,ParticleIndex)
  END IF
  ! Species index given from the initialization
  PartSpecies(ParticleIndex) = iSpec
  ! Get the species index of the background gas
  LocalElemID = PEM%LocalElemID(NewParticleIndex)
  bgSpec = BGGas_GetSpecies(LocalElemID)
  PartSpecies(NewParticleIndex) = bgSpec
  ! Particle element
  PEM%GlobalElemID(NewParticleIndex) = PEM%GlobalElemID(ParticleIndex)
  IF(CollisMode.GT.1) THEN
    IF(SpecDSMC(bgSpec)%PolyatomicMol) THEN
      CALL DSMC_SetInternalEnr_Poly(bgSpec,1,NewParticleIndex,1)
    ELSE
      CALL DSMC_SetInternalEnr_LauxVFD(bgSpec,1,NewParticleIndex,1)
    END IF
  END IF
  IF(BGGas%UseDistribution) THEN
    iBGGSpec = BGGas%MapSpecToBGSpec(bgSpec)
    PartState(4:6,NewParticleIndex) = CalcVelocity_maxwell_particle(bgSpec,BGGas%Distribution(iBGGSpec,4:6,LocalElemID)) &
                                  + BGGas%Distribution(iBGGSpec,1:3,LocalElemID)
  ELSE
    CALL CalcVelocity_maxwell_lpn(FractNbr=bgSpec, Vec3D=PartState(4:6,NewParticleIndex), iInit=1)
  END IF
  ! Particle flags
  PDM%ParticleInside(NewParticleIndex)  = .TRUE.
  PDM%IsNewPart(NewParticleIndex)       = .TRUE.
  PDM%dtFracPush(NewParticleIndex)      = .FALSE.
  ! Last element ID
  PEM%LastGlobalElemID(NewParticleIndex) = PEM%GlobalElemID(NewParticleIndex)
  PEM%LastGlobalElemID(ParticleIndex) = PEM%GlobalElemID(ParticleIndex)
  ! Pairing (first particle is the background gas species)
  Coll_pData(iNewPart)%iPart_p1 = NewParticleIndex
  Coll_pData(iNewPart)%iPart_p2 = ParticleIndex
  ! Relative velocity is not required as the relative translational energy will not be considered
  Coll_pData(iNewPart)%CRela2 = 0.
  ! Weighting factor: use the weighting factor of the emission init
  IF(usevMPF) THEN
    PartMPF(ParticleIndex)    = Species(iSpec)%MacroParticleFactor
    PartMPF(NewParticleIndex) = PartMPF(ParticleIndex)
  END IF
  ! Velocity (set it to zero, as it will be subtracted in the chemistry module)
  PartState(4:6,ParticleIndex) = 0.
  ! Internal energies (set it to zero)
  PartStateIntEn(1:2,ParticleIndex) = 0.
  IF(DSMC%ElectronicModel.GT.0) PartStateIntEn(3,ParticleIndex) = 0.
END DO


! Add the particles initialized through the emission and the background particles
! Update the current next free position
IF(iNewPart+NbrOfParticle.GT.0) THEN
  PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + NbrOfParticle + iNewPart
  PDM%ParticleVecLength = MAX(PDM%ParticleVecLength,GetNextFreePosition(0))
END IF
#ifdef CODE_ANALYZE
IF(PDM%ParticleVecLength.GT.PDM%maxParticleNumber) CALL Abort(__STAMP__,'PDM%ParticleVeclength exceeds PDM%maxParticleNumber, Difference:',IntInfoOpt=PDM%ParticleVeclength-PDM%maxParticleNumber)
DO iPart=PDM%ParticleVecLength+1,PDM%maxParticleNumber
  IF (PDM%ParticleInside(iPart)) THEN
    IPWRITE(*,*) iPart,PDM%ParticleVecLength,PDM%maxParticleNumber
    CALL Abort(__STAMP__,'Particle outside PDM%ParticleVeclength',IntInfoOpt=iPart)
  END IF
END DO
#endif


!> 4.) Perform the reaction, distribute the collision energy (including photon energy) and emit electrons perpendicular
!>     to the photon's path
ASSOCIATE(b1          => Species(iSpec)%Init(iInit)%NormalVector1IC(1:3) ,&
          b2          => Species(iSpec)%Init(iInit)%NormalVector2IC(1:3) ,&
          normal      => Species(iSpec)%Init(iInit)%NormalIC             ,&
          PartBCIndex => Species(iSpec)%Init(iInit)%PartBCIndex)
IF(NbrOfPhotonXsecReactions.GT.0)THEN
  DO iPart = 1, SUM(NumPhotoIonization(:))
    ! Loop over all randomized lines (found above)
    DO iLine = iLineStart, iLineEnd
      ! Check if level is occupied
      DO WHILE(PhotonEnergies(iLine,1).GT.0)
        ! Reduce the level counter by one
        PhotonEnergies(iLine,1) = PhotonEnergies(iLine,1) - 1
        ! Check if cross-section is occupied
        DO iPhotoReac = 1, NbrOfPhotonXsecReactions
          IF(PhotonEnergies(iLine,1+iPhotoReac).GT.0)THEN
            ! Reduce cross-section by one
            PhotonEnergies(iLine,1+iPhotoReac) = PhotonEnergies(iLine,1+iPhotoReac) - 1
            iPair = iPair + 1
            CALL PhotoIonization_InsertProducts(iPair, PhotoReacToReac(iPhotoReac), b1, b2, normal, iLineOpt=iLine, PartBCIndex=PartBCIndex)
          END IF ! PhotonEnergies(iLine,1+iPhotoReac).GT.0
        END DO ! iPhotoReac = 1, NbrOfPhotonXsecReactions
      END DO
    END DO ! iLine = iLineStart, iLineEnd
  END DO
ELSE
  DO iReac = 1, ChemReac%NumOfReact
    IF(TRIM(ChemReac%ReactModel(iReac)).NE.'phIon') CYCLE
    DO iPart = 1, NumPhotoIonization(iReac)
      iPair = iPair + 1
      CALL PhotoIonization_InsertProducts(iPair, iReac, b1, b2, normal, PartBCIndex=PartBCIndex)
    END DO
  END DO
END IF ! NbrOfPhotonXsecReactions.GT.0
END ASSOCIATE

DSMCSumOfFormedParticles = 0

DEALLOCATE(Coll_pData)

IF (DSMC%DoAmbipolarDiff) THEN
  ! Every particle created during photo-ionization is stored in the iPartIndx_NodeNewAmbi array, which is treated in the routine
  CALL AD_DeleteParticles()
END IF

END SUBROUTINE BGGas_PhotoIonization


!===================================================================================================================================
!> Read-in of the element data from a DSMC state and utilization as a cell-local background gas distribution
!===================================================================================================================================
SUBROUTINE BGGas_ReadInDistribution()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_io_hdf5
USE MOD_HDF5_Input       ,ONLY: OpenDataFile,CloseDataFile,ReadArray,GetDataSize,ReadAttribute
USE MOD_HDF5_Input       ,ONLY: nDims,HSize,File_ID
USE MOD_Restart_Vars     ,ONLY: MacroRestartFileName
USE MOD_Mesh_Vars        ,ONLY: offsetElem, nElems,nGlobalElems
USE MOD_DSMC_Vars        ,ONLY: BGGas
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: nVarHDF5, nElems_HDF5, iVar, iSpec, iElem, iBGGSpec, nSpecReadin
REAL, ALLOCATABLE                 :: ElemDataHDF5(:,:)
!===================================================================================================================================

LBWRITE(UNIT_stdOut,*) 'BGGas distribution - Using macroscopic values from file: ',TRIM(MacroRestartFileName)

CALL OpenDataFile(MacroRestartFileName,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)

CALL GetDataSize(File_ID,'ElemData',nDims,HSize,attrib=.FALSE.)
nVarHDF5  = INT(HSize(1),4)
IF(nVarHDF5.LT.10) CALL abort(__STAMP__,'Number of variables .h5 file is less than 10')

nElems_HDF5 = INT(HSize(2),4)
IF(nElems_HDF5.NE.nGlobalElems) CALL abort(__STAMP__,'Number of global elements does not match number of elements in .h5 file')

DEALLOCATE(HSize)

ALLOCATE(ElemDataHDF5(1:nVarHDF5,1:nElems))
! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
  nVarHDF5   => INT(nVarHDF5,IK) ,&
  offsetElem => INT(offsetElem,IK),&
  nElems     => INT(nElems,IK)    )
  CALL ReadArray('ElemData',2,(/nVarHDF5,nElems/),offsetElem,2,RealArray=ElemDataHDF5(:,:))
END ASSOCIATE

CALL ReadAttribute(File_ID,'NSpecies',1,IntScalar=nSpecReadin)

! Loop over all the read-in species and map them to the background gas species
iVar = 1
DO iSpec = 1, nSpecReadin
  DO iBGGSpec = 1, BGGas%NumberOfSpecies
    IF(BGGas%DistributionSpeciesIndex(BGGas%MapBGSpecToSpec(iBGGSpec)).EQ.iSpec) THEN
      DO iElem = 1, nElems
        BGGas%Distribution(iBGGSpec,1:10,iElem) = ElemDataHDF5(iVar:iVar-1+10,iElem)
      END DO
      LBWRITE(UNIT_stdOut,*) 'BGGas distribution: Mapped read-in values of species ', iSpec, ' to current species ', &
          BGGas%MapBGSpecToSpec(iBGGSpec)
    END IF
  END DO
  iVar = iVar + DSMC_NVARS
END DO

DEALLOCATE(ElemDataHDF5)

CALL CloseDataFile()

END SUBROUTINE BGGas_ReadInDistribution


!===================================================================================================================================
!> Multi-species background with a defined trace species: Routine splits (clones) the particle species based on the difference
!> between the maximum weighting factor of the background gas and the current trace background species
!===================================================================================================================================
SUBROUTINE BGGas_TraceSpeciesSplit(iElem, nPart, nPair)
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars             ,ONLY: BGGas, CollisMode, PartStateIntEn, DSMC
USE MOD_DSMC_Vars             ,ONLY: DSMC, SpecDSMC, VibQuantsPar, PolyatomMolDSMC
USE MOD_Particle_Vars         ,ONLY: PEM,PartSpecies,PartState,PartMPF,Species
USE MOD_Part_Tools            ,ONLY: GetNextFreePosition
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iElem
INTEGER, INTENT(INOUT)        :: nPart, nPair
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart, iLoop, iNewPart, iSplit, LocalElemID, iSpec, bggSpec, iSplitPart, SplitPartNum, PartIndex
INTEGER                       :: bggPartIndex
REAL                          :: MPFRatio
!===================================================================================================================================
iNewPart = 0

iPart = PEM%pStart(iElem)
DO iLoop = 1, nPart
  iSpec  = PartSpecies(iPart)
  ! Skip background particles that have been created within this loop
  IF(BGGas%BackgroundSpecies(iSpec)) THEN
    iPart = PEM%pNext(iPart)
    CYCLE
  END IF
  bggSpec = PartSpecies(BGGas%PairingPartner(iPart))
  ! Skip pairs with regular background species
  IF(.NOT.BGGas%TraceSpecies(bggSpec)) THEN
    iPart = PEM%pNext(iPart)
    CYCLE
  END IF
  ! Split required if particle MPF is larger than BGGas MPF
  IF(PartMPF(iPart).GT.Species(bggSpec)%MacroParticleFactor) THEN
    MPFRatio = BGGas%MaxMPF/Species(bggSpec)%MacroParticleFactor
    ! Avoid artificial increase of split number due to increasing MPF of the particle species (due to an additional merge process)
    IF(PartMPF(iPart) / Species(bggSpec)%MacroParticleFactor.GT.MPFRatio) THEN
      SplitPartNum = NINT(MPFRatio)
    ELSE
      SplitPartNum = NINT(PartMPF(iPart) / Species(bggSpec)%MacroParticleFactor)
    END IF
    nPair   = nPair + (SplitPartNum - 1)          ! New Number of Pairs
    PartIndex = 0
    PartMPF(iPart) = PartMPF(iPart) / REAL(SplitPartNum) ! set new MPF for test particle
    PartMPF(BGGas%PairingPartner(iPart)) = PartMPF(iPart)
    iSplitPart = 0
    DO iSplit = 2, SplitPartNum
      ! --- Create clone of test particle
      iNewPart = iNewPart + 1
      iSplitPart = iSplitPart + 1
      PartIndex = GetNextFreePosition()
      ! Assign properties but do not use the velocity and energy of the background gas
      CALL BGGas_AssignParticleProperties(iSpec,iPart,PartIndex,GetVelocity_opt=.FALSE.,GetInternalEnergy_opt=.TRUE.)
      ! Copy properties from the particle species
      PartState(4:6,PartIndex) = PartState(4:6,iPart)
      IF(CollisMode.GT.1) THEN
        PartStateIntEn(1:2,PartIndex) = PartStateIntEn(1:2,iPart)
        IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
          IF(ALLOCATED(VibQuantsPar(PartIndex)%Quants)) DEALLOCATE(VibQuantsPar(PartIndex)%Quants)
          ALLOCATE(VibQuantsPar(PartIndex)%Quants(PolyatomMolDSMC(SpecDSMC(iSpec)%SpecToPolyArray)%VibDOF))
          VibQuantsPar(PartIndex)%Quants(:) = VibQuantsPar(iPart)%Quants(:)
        END IF
        IF(DSMC%ElectronicModel.GT.0) PartStateIntEn(3,PartIndex) = PartStateIntEn(3,iPart)
      END IF
      ! Update cell-local particle list
      LocalElemID = PEM%LocalElemID(PartIndex)
      PEM%pNext(PEM%pEnd(LocalElemID)) = PartIndex     ! Next Particle of same Elem (Linked List)
      PEM%pEnd(LocalElemID) = PartIndex
      PEM%pNumber(LocalElemID) = PEM%pNumber(LocalElemID) + 1
      ! --- Create background particles
      iNewPart = iNewPart + 1
      iSplitPart = iSplitPart + 1
      ! Get a free particle index
      bggPartIndex = GetNextFreePosition()
      ! Set the pairing partner
      BGGas%PairingPartner(PartIndex) = bggPartIndex
      ! Assign properties of the background gas
      CALL BGGas_AssignParticleProperties(bggSpec,iPart,bggPartIndex)
      ! Update cell-local particle list
      LocalElemID = PEM%LocalElemID(bggPartIndex)
      PEM%pNext(PEM%pEnd(LocalElemID)) = bggPartIndex
      PEM%pEnd(LocalElemID) = bggPartIndex
      PEM%pNumber(LocalElemID) = PEM%pNumber(LocalElemID) + 1
    END DO
  END IF
  iPart = PEM%pNext(iPart)
END DO
! Set the new number of particles
nPart = PEM%pNumber(iElem)

END SUBROUTINE BGGas_TraceSpeciesSplit


!===================================================================================================================================
!> Initialization of the background gas regions (e.g. geometrical volumes with different background gas properties)
!> 1. Read-in geometry information based on selected type (e.g. cylinder)
!> 2. Determine which elements are the defined regions by comparing the element midpoint
!> 3. Write the corresponding region properties into the BGGas%Distribution array
!> 4. Calculate the element local species fraction in case of a multi-species background gas
!> 5. Activate BGGas%UseDistribution to utilize the same arrays and routines during computation (but do not use it for read-in)
!===================================================================================================================================
SUBROUTINE BGGas_InitRegions()
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_ReadInTools
USE MOD_Mesh_Tools         ,ONLY: GetCNElemID
USE MOD_Mesh_Vars          ,ONLY: nElems, offSetElem
USE MOD_DSMC_Vars          ,ONLY: BGGas
USE MOD_Particle_Vars      ,ONLY: Species,nSpecies
USE MOD_Particle_Mesh_Vars ,ONLY: ElemMidPoint_Shared
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(32)                 :: hilf2
INTEGER                       :: iElem, iSpec, bgSpec, iInit, iReg, CNElemID
REAL                          :: lineVector(3), nodeVec(3), nodeRadius, nodeHeight
!===================================================================================================================================
LBWRITE(UNIT_stdOut,'(A)') ' INIT BACKGROUND GAS REGIONS ...'

ALLOCATE(BGGas%Region(BGGas%nRegions))
ALLOCATE(BGGas%RegionElemType(nElems))
BGGas%RegionElemType = 0

! 1) Read-in of the background gas regions
DO iReg = 1, BGGas%nRegions
  WRITE(UNIT=hilf2,FMT='(I0)') iReg
  BGGas%Region(iReg)%Type = TRIM(GETSTR('Particles-BGGas-Region'//TRIM(hilf2)//'-Type'))
  SELECT CASE(TRIM(BGGas%Region(iReg)%Type))
    CASE('cylinder')
      BGGas%Region(iReg)%RadiusIC               = GETREAL('Particles-BGGas-Region'//TRIM(hilf2)//'-RadiusIC')
      BGGas%Region(iReg)%Radius2IC              = GETREAL('Particles-BGGas-Region'//TRIM(hilf2)//'-Radius2IC')
      BGGas%Region(iReg)%CylinderHeightIC       = GETREAL('Particles-BGGas-Region'//TRIM(hilf2)//'-CylinderHeightIC')
      IF(BGGas%Region(iReg)%Radius2IC.GE.BGGas%Region(iReg)%RadiusIC) CALL abort(__STAMP__,&
          'For this emission type RadiusIC must be greater than Radius2IC!')
      !--- Get BasePointIC
      BGGas%Region(iReg)%BasePointIC   = GETREALARRAY('Particles-BGGas-Region'//TRIM(hilf2)//'-BasePointIC',3)
      !--- Get BaseVector1IC
      BGGas%Region(iReg)%BaseVector1IC = GETREALARRAY('Particles-BGGas-Region'//TRIM(hilf2)//'-BaseVector1IC',3)
      !--- Get BaseVector2IC
      BGGas%Region(iReg)%BaseVector2IC = GETREALARRAY('Particles-BGGas-Region'//TRIM(hilf2)//'-BaseVector2IC',3)
      ! Determine the normal vector of the cylinder from base vectors
      BGGas%Region(iReg)%NormalVector = CROSS(BGGas%Region(iReg)%BaseVector1IC,BGGas%Region(iReg)%BaseVector2IC)
      IF (VECNORM(BGGas%Region(iReg)%NormalVector).EQ.0) THEN
        CALL abort(__STAMP__,'BaseVectors are parallel!')
      ELSE
        BGGas%Region(iReg)%NormalVector = UNITVECTOR(BGGas%Region(iReg)%NormalVector)
      END IF
    CASE DEFAULT
      CALL abort(__STAMP__,'ERROR Background gas regions: Selected region type is not implemented!')
  END SELECT
END DO

DO iElem = 1, nElems
  ! 2) Map elements to regions
  CNElemID = GetCNElemID(iElem+offSetElem)
  DO iReg = 1, BGGas%nRegions
    SELECT CASE(TRIM(BGGas%Region(iReg)%Type))
    CASE('cylinder')
      lineVector = BGGas%Region(iReg)%NormalVector
      ! Node vector relative to cylinder basepoint
      nodeVec(1:3) = ElemMidPoint_Shared(1:3,CNElemID) - BGGas%Region(iReg)%BasePointIC
      ! Node vector projected onto cylinder normal
      nodeHeight = DOT_PRODUCT(nodeVec(1:3),lineVector)
      ! Node radius as the remainder of the node vector length and the projected node height
      nodeRadius = SQRT(DOTPRODUCT(nodeVec(1:3)) - nodeHeight**2)
      ! Check if node is outside of the region
      IF((nodeHeight.GE.0.).AND.(nodeHeight.LE.BGGas%Region(iReg)%CylinderHeightIC) &
        .AND.(nodeRadius.GE.BGGas%Region(iReg)%Radius2IC).AND.(nodeRadius.LE.BGGas%Region(iReg)%RadiusIC)) THEN
        ! Element mid point is inside (positive region number)
        IF(BGGas%RegionElemType(iElem).NE.0) THEN
          CALL abort(__STAMP__,'ERROR Background gas regions: Overlapping regions are not supported!')
        END IF
        BGGas%RegionElemType(iElem) = iReg
      END IF
    END SELECT
  END DO    ! iReg = 1, BGGas%nRegions
  ! 3) Assign the region values to the distribution array
  IF(BGGas%RegionElemType(iElem).NE.0) THEN
    DO iSpec = 1, nSpecies
      ! Skip non-background gas species
      IF(.NOT.BGGas%BackgroundSpecies(iSpec)) CYCLE
      ! Loop over all the inits for the species (different inits for different regions)
      DO iInit = 1, Species(iSpec)%NumberOfInits
        ASSOCIATE(  SpecInit => Species(iSpec)%Init(iInit) )
          IF(BGGas%RegionElemType(iElem).EQ.Species(iSpec)%Init(iInit)%BGGRegion) THEN
            ! Velocity
            BGGas%Distribution(BGGas%MapSpecToBGSpec(iSpec),1:3,iElem) = SpecInit%VeloIC * SpecInit%VeloVecIC
            ! Translational temperature
            BGGas%Distribution(BGGas%MapSpecToBGSpec(iSpec),4:6,iElem) = SpecInit%MWTemperatureIC
            ! Number density
            BGGas%Distribution(BGGas%MapSpecToBGSpec(iSpec),7,iElem) = SpecInit%PartDensity
          END IF        ! BGGas%RegionElemType(iElem).EQ.Species(iSpec)%Init(iInit)%BGGRegion
        END ASSOCIATE
      END DO            ! iInit = 1, Species(iSpec)%NumberOfInits
    END DO              ! iSpec = 1, nSpecies
  END IF                ! BGGas%RegionElemType(iElem).NE.0
  ! 4) Calculate species fraction per element
  DO bgSpec = 1, BGGas%NumberOfSpecies
    IF(SUM(BGGas%Distribution(:,7,iElem)).GT.0.)THEN
      BGGas%SpeciesFractionElem(bgSpec,iElem) = BGGas%Distribution(bgSpec,7,iElem) / SUM(BGGas%Distribution(:,7,iElem))
    END If
  END DO
END DO                  ! iElem = 1, nElems

! 5) Utilizing the same routines after the initialization as the read-in distribution
BGGas%UseDistribution = .TRUE.

LBWRITE(UNIT_stdOut,'(A)') ' BACKGROUND GAS REGIONS DONE!'

END SUBROUTINE BGGas_InitRegions


!===================================================================================================================================
!> Background gas regions: Set the internal temperatures in case of DSMC and CollisMode = 2/3 (not yet available during
!> BGGas_InitRegions). Loop over all elements, species and inits per species to set values for molecules and/or atoms.
!===================================================================================================================================
SUBROUTINE BGGas_RegionsSetInternalTemp()
! MODULES
USE MOD_Mesh_Vars             ,ONLY: nElems
USE MOD_DSMC_Vars             ,ONLY: BGGas, DSMC, SpecDSMC
USE MOD_Particle_Vars         ,ONLY: Species,nSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iElem, iSpec, iInit
!===================================================================================================================================

DO iElem = 1, nElems
  ! Assign the region values to the distribution array
  IF(BGGas%RegionElemType(iElem).NE.0) THEN
    DO iSpec = 1, nSpecies
      ! Skip non-background gas species
      IF(.NOT.BGGas%BackgroundSpecies(iSpec)) CYCLE
      ! Loop over all the inits for the species (different inits for different regions)
      DO iInit = 1, Species(iSpec)%NumberOfInits
        IF(BGGas%RegionElemType(iElem).EQ.Species(iSpec)%Init(iInit)%BGGRegion) THEN
          IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
            ! Vibrational temperature
            BGGas%Distribution(BGGas%MapSpecToBGSpec(iSpec),8,iElem) = SpecDSMC(iSpec)%Init(iInit)%TVib
            ! Rotational temperature
            BGGas%Distribution(BGGas%MapSpecToBGSpec(iSpec),9,iElem) = SpecDSMC(iSpec)%Init(iInit)%TRot
          END IF
          IF (DSMC%ElectronicModel.GT.0) THEN
            ! Electronic temperature
            BGGas%Distribution(BGGas%MapSpecToBGSpec(iSpec),10,iElem) = SpecDSMC(iSpec)%Init(iInit)%Telec
          END IF
        END IF        ! BGGas%RegionElemType(iElem).EQ.Species(iSpec)%Init(iInit)%BGGRegion
      END DO            ! iInit = 1, Species(iSpec)%NumberOfInits
    END DO              ! iSpec = 1, nSpecies
  END IF                ! BGGas%RegionElemType(iElem).NE.0
END DO                  ! iElem = 1, nElems

END SUBROUTINE BGGas_RegionsSetInternalTemp

END MODULE MOD_DSMC_BGGas
