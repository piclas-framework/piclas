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
PUBLIC :: BGGas_Initialize, BGGas_InsertParticles, DSMC_pairing_bggas, BGGas_DeleteParticles, BGGas_AssignParticleProperties
PUBLIC :: BGGas_PhotoIonization
!===================================================================================================================================

CONTAINS

SUBROUTINE BGGas_Initialize()
!===================================================================================================================================
!> Initialization of the background gas: compatibility check, array allocation, background species to species mapping and
!> calculation of the molar fraction
!===================================================================================================================================
! MODULES
USE MOD_Globals


USE MOD_ReadInTools
USE MOD_Globals       ,ONLY: abort
USE MOD_DSMC_Vars     ,ONLY: BGGas
USE MOD_Mesh_Vars     ,ONLY: nElems
USE MOD_Particle_Vars ,ONLY: PDM, Symmetry, Species, nSpecies, VarTimeStep
USE MOD_Restart_Vars  ,ONLY: DoMacroscopicRestart, MacroRestartFileName
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
IF((Symmetry%Order.EQ.2).OR.VarTimeStep%UseVariableTimeStep) THEN
  CALL abort(__STAMP__,'ERROR: 2D/Axisymmetric and variable timestep are not implemented with a background gas yet!')
END IF

DO iSpec = 1, nSpecies
  IF(BGGas%BackgroundSpecies(iSpec)) THEN
    IF(Species(iSpec)%NumberOfInits.NE.1) CALL abort(__STAMP__, 'BGG species can be used ONLY for BGG!')
    IF(.NOT.BGGas%UseDistribution)THEN
      IF(BGGas%NumberDensity(iSpec).EQ.0.) CALL abort(__STAMP__, 'NumberDensity is zero but must be defined for a background gas!')
    END IF ! .NOT.BGGas%UseDistribution
  END IF
END DO

IF(DoMacroscopicRestart) CALL abort(__STAMP__, 'Constant background gas and macroscopic restart are not compatible!')

! 2.) Allocation
IF(BGGas%UseDistribution) THEN
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
    BGGas%MapSpecToBGSpec(iSpec)  = bgSpec
    BGGas%MapBGSpecToSpec(bgSpec) = iSpec
    BGGas%MaxMPF = MAX(BGGas%MaxMPF,Species(iSpec)%MacroParticleFactor)
  END IF
END DO

! 4.) Read-in a background gas distribution
IF(BGGas%UseDistribution) CALL BGGas_ReadInDistribution()

! 5.) Determine species fraction
DO bgSpec = 1, BGGas%NumberOfSpecies
  IF(BGGas%UseDistribution) THEN
    DO iElem = 1, nElems
      BGGas%SpeciesFractionElem(bgSpec,iElem) = BGGas%Distribution(bgSpec,7,iElem) / SUM(BGGas%Distribution(:,7,iElem))
    END DO ! iElem = 1, nElems
  ELSE
    BGGas%NumberDensity(bgSpec)   = SpeciesDensTmp(iSpec)
    BGGas%SpeciesFraction(bgSpec) = BGGas%NumberDensity(bgSpec) / SUM(SpeciesDensTmp)
  END IF
END DO ! bgSpec = 1, BGGas%NumberOfSpecies

END SUBROUTINE BGGas_Initialize


INTEGER FUNCTION BGGas_GetSpecies()
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
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: iRan
INTEGER           :: iSpec
!===================================================================================================================================

BGGas_GetSpecies = 0

IF(BGGas%NumberOfSpecies.GT.1) THEN
  CALL RANDOM_NUMBER(iRan)
  DO iSpec = 1, BGGas%NumberOfSpecies
    IF(SUM(BGGas%SpeciesFraction(1:iSpec)).GT.iRan) THEN
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
    ! Get a free particle index
    iNewPart = iNewPart + 1
    PositionNbr = PDM%nextFreePosition(iNewPart+PDM%CurrentNextFreePosition)
    IF (PositionNbr.EQ.0) THEN
      CALL Abort(__STAMP__,'ERROR in BGGas: MaxParticleNumber should be increased to account for the BGG particles!')
    END IF
    ! Get the background gas species
    iSpec = BGGas_GetSpecies()
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
! Increase the particle vector length and update the linked list
PDM%ParticleVecLength = MAX(PDM%ParticleVecLength,PositionNbr)
PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + iNewPart

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
USE MOD_Particle_Vars           ,ONLY: PDM, PEM, PartState,PartSpecies,PartPosRef, VarTimeStep, usevMPF, PartMPF
USE MOD_DSMC_Vars               ,ONLY: CollisMode, SpecDSMC, BGGas
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
USE MOD_part_emission_tools     ,ONLY: CalcVelocity_maxwell_lpn, DSMC_SetInternalEnr_LauxVFD
USE MOD_DSMC_PolyAtomicModel    ,ONLY: DSMC_SetInternalEnr_Poly
USE MOD_Macro_Restart           ,ONLY: CalcVelocity_maxwell_particle
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
IF(VarTimeStep%UseVariableTimeStep) VarTimeStep%ParticleTimeStep(bggPartIndex) = VarTimeStep%ParticleTimeStep(PartIndex)

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
                                            + BGGas%SpeciesFraction(bggSpec) * SUM(BGGas%Distribution(bggSpec,4:6,iElem)) / 3.
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
    IF(usevMPF) THEN
      CollInf%Coll_SpecPartNum(iSpec) = BGGas%NumberDensity(bggSpec)*ElemVolume_Shared(CNElemID)
    ELSEIF(BGGas%UseDistribution)THEN
      CollInf%Coll_SpecPartNum(iSpec) = BGGas%Distribution(bggSpec,7,iElem)&
                                        * ElemVolume_Shared(CNElemID) / Species(iSpec)%MacroParticleFactor
    ELSE
      CollInf%Coll_SpecPartNum(iSpec) = BGGas%NumberDensity(bggSpec)*ElemVolume_Shared(CNElemID)/Species(iSpec)%MacroParticleFactor
    END IF
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
USE MOD_DSMC_Vars,          ONLY : BGGas
USE MOD_PARTICLE_Vars,      ONLY : PDM, PartSpecies
USE MOD_part_tools,         ONLY : UpdateNextFreePosition
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iPart
!===================================================================================================================================
DO iPart = 1, PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    IF(BGGas%BackgroundSpecies(PartSpecies(iPart))) PDM%ParticleInside(iPart) = .FALSE.
  END IF
END DO
BGGas%PairingPartner = 0

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
USE MOD_Macro_Restart          ,ONLY: CalcVelocity_maxwell_particle
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec,iInit,TotalNbrOfReactions
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart, iPair, iNewPart, iReac, ParticleIndex, NewParticleIndex, bgSpec, NbrOfParticle, LocalElemID
REAL                          :: RandVal,NumTmp,ProbRest
INTEGER                       :: TotalNbrOfReactionsTmp,iCrossSection,NbrCrossSections,iBGGSpec
INTEGER                       :: NumPhotoIonization(ChemReac%NumOfReact)
REAL                          :: SumCrossSections
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
  IF(TRIM(ChemReac%ReactModel(iReac)).NE.'phIon') CYCLE
  SumCrossSections = SumCrossSections + ChemReac%CrossSection(iReac)
  NbrCrossSections = NbrCrossSections + 1
END DO ! iReac = 1, ChemReac%NumOfReact

!> 1.) Compute the number of photoionization events in the local domain of each proc
iCrossSection = 0
DO iReac = 1, ChemReac%NumOfReact
  ! Only treat photoionization reactions
  IF(TRIM(ChemReac%ReactModel(iReac)).NE.'phIon') CYCLE
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

!> 2.) Delete left-over inserted particles
IF(TotalNbrOfReactions.GT.SUM(NumPhotoIonization)) THEN
  DO iPart = SUM(NumPhotoIonization)+1,TotalNbrOfReactions
    PDM%ParticleInside(PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)) = .FALSE.
  END DO
ELSE IF(TotalNbrOfReactions.LT.SUM(NumPhotoIonization)) THEN
  CALL Abort(&
    __STAMP__&
    ,'ERROR in PhotoIonization: Something is wrong, trying to perform more reactions than anticipated!')
END IF

IF(SUM(NumPhotoIonization).EQ.0) RETURN

!> 3.) Insert the products of the photoionization rections
NbrOfParticle = SUM(NumPhotoIonization)

ALLOCATE(Coll_pData(NbrOfParticle))
Coll_pData%Ec=0.
DSMCSumOfFormedParticles = 0

iNewPart = 0; iPair = 0

DO iPart = 1, NbrOfParticle
  ! Loop over the particles with a set position (from SetParticlePosition)
  ParticleIndex = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
  IF (DSMC%DoAmbipolarDiff) THEN
    newAmbiParts = newAmbiParts + 1
    iPartIndx_NodeNewAmbi(newAmbiParts) = ParticleIndex
  END IF
  iNewPart = iNewPart + 1
  ! Get a new index for the second product
  NewParticleIndex = PDM%nextFreePosition(iNewPart+PDM%CurrentNextFreePosition+NbrOfParticle)
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
  bgSpec = BGGas_GetSpecies()
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
    LocalElemID = PEM%LocalElemID(NewParticleIndex)
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
  ! Velocity (set it to zero, as it will be substracted in the chemistry module)
  PartState(4:6,ParticleIndex) = 0.
  ! Internal energies (set it to zero)
  PartStateIntEn(1:2,ParticleIndex) = 0.
  IF(DSMC%ElectronicModel.GT.0) PartStateIntEn(3,ParticleIndex) = 0.
END DO

! Add the particles initialized through the emission and the background particles
PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfParticle + iNewPart
! Update the current next free position
PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + NbrOfParticle + iNewPart

IF(PDM%ParticleVecLength.GT.PDM%MaxParticleNumber) THEN
  CALL Abort(__STAMP__&
    ,'ERROR in PhotoIonization: ParticleVecLength greater than MaxParticleNumber! Increase the MaxParticleNumber to at least: ' &
    , IntInfoOpt=PDM%ParticleVecLength)
END IF

!> 4.) Perform the reaction, distribute the collision energy (including photon energy) and emit electrons perpendicular
!>     to the photon's path
DO iReac = 1, ChemReac%NumOfReact
  ! Only treat photoionization reactions
  IF(TRIM(ChemReac%ReactModel(iReac)).NE.'phIon') CYCLE
  DO iPart = 1, NumPhotoIonization(iReac)
    iPair = iPair + 1
    CALL PhotoIonization_InsertProducts(iPair, iReac, iInit, iSpec)
  END DO
END DO

! Advance particle vector length and the current next free position with newly created particles
PDM%ParticleVecLength = PDM%ParticleVecLength + DSMCSumOfFormedParticles
PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + DSMCSumOfFormedParticles

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
USE MOD_HDF5_Input    ,ONLY: OpenDataFile,CloseDataFile,ReadArray,GetDataSize,ReadAttribute
USE MOD_HDF5_Input    ,ONLY: nDims,HSize,File_ID
USE MOD_Restart_Vars  ,ONLY: MacroRestartFileName
USE MOD_Mesh_Vars     ,ONLY: offsetElem, nElems,nGlobalElems
USE MOD_DSMC_Vars     ,ONLY: BGGas
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

SWRITE(UNIT_stdOut,*) 'BGGas distribution - Using macroscopic values from file: ',TRIM(MacroRestartFileName)

CALL OpenDataFile(MacroRestartFileName,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)

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
      SWRITE(UNIT_stdOut,*) 'BGGas distribution: Mapped read-in values of species ', iSpec, ' to current species ', &
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
USE MOD_Particle_Vars         ,ONLY: PDM,PEM,PartSpecies,PartState,PartMPF,Species
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
      PartIndex = PDM%nextFreePosition(iNewPart+PDM%CurrentNextFreePosition)
      IF (PartIndex.EQ.0) THEN
        CALL Abort(__STAMP__,'ERROR in BGGas: MaxParticleNumber should be increased to account for the BGG particles!')
      END IF
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
      bggPartIndex = PDM%nextFreePosition(iNewPart+PDM%CurrentNextFreePosition)
      IF (bggPartIndex.EQ.0) THEN
        CALL Abort(__STAMP__,'ERROR in BGGas: MaxParticleNumber should be increased to account for the BGG particles!')
      END IF
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
! Increase the particle vector length and the position in the linked list
PDM%ParticleVecLength = MAX(PDM%ParticleVecLength,bggPartIndex)
PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + iNewPart
! Set the new number of particles
nPart = PEM%pNumber(iElem)

END SUBROUTINE BGGas_TraceSpeciesSplit

END MODULE MOD_DSMC_BGGas
