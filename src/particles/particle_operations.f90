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

MODULE MOD_part_operations
!===================================================================================================================================
! Contains tools for particle related operations. This routine is uses MOD_Particle_Boundary_Tools, but not vice versa!
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
PUBLIC :: CreateParticle, RemoveParticle
PUBLIC :: RemoveAllElectrons
!===================================================================================================================================

CONTAINS

SUBROUTINE CreateParticle(SpecID,Pos,GlobElemID,Velocity,RotEnergy,VibEnergy,ElecEnergy,NewPartID,NewMPF)
!===================================================================================================================================
!> creates a single particle at correct array position and assign properties
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars           ,ONLY: PDM, PEM, PartState, LastPartPos, PartSpecies,PartPosRef, Species, usevMPF, PartMPF
USE MOD_Particle_Vars           ,ONLY: UseVarTimeStep, PartTimeStep
USE MOD_Particle_Vars           ,ONLY: UseRotRefFrame, RotRefFrameOmega, PartVeloRotRef
USE MOD_DSMC_Vars               ,ONLY: useDSMC, CollisMode, DSMC, PartStateIntEn, RadialWeighting
USE MOD_DSMC_Vars               ,ONLY: newAmbiParts, iPartIndx_NodeNewAmbi
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
USE MOD_Eval_xyz                ,ONLY: GetPositionInRefElem
USE MOD_part_tools              ,ONLY: CalcRadWeightMPF
USE MOD_Particle_TimeStep       ,ONLY: GetParticleTimeStep
USE MOD_Part_Tools              ,ONLY: InRotRefFrameCheck, GetNextFreePosition
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER, INTENT(IN)           :: SpecID           !< Species ID
REAL, INTENT(IN)              :: Pos(1:3)         !< Position (x,y,z)
INTEGER, INTENT(IN)           :: GlobElemID       !< global element ID
REAL, INTENT(IN)              :: Velocity(1:3)    !< Velocity (vx,vy,vz)
REAL, INTENT(IN)              :: RotEnergy        !< Rotational energy
REAL, INTENT(IN)              :: VibEnergy        !< Vibrational energy
REAL, INTENT(IN)              :: ElecEnergy       !< Electronic energy
INTEGER, INTENT(OUT),OPTIONAL :: NewPartID        !< ID of newly created particle
REAL, INTENT(IN),OPTIONAL     :: NewMPF           !< MPF of newly created particle
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER :: newParticleID
!===================================================================================================================================

newParticleID = GetNextFreePosition()

PartSpecies(newParticleID)     = SpecID
LastPartPos(1:3,newParticleID) = Pos(1:3)
PartState(1:3,newParticleID)   = Pos(1:3)
PartState(4:6,newParticleID)   = Velocity(1:3)

! Set the new reference position here
IF(TrackingMethod.EQ.REFMAPPING)THEN
  CALL GetPositionInRefElem(PartState(1:3,newParticleID),PartPosRef(1:3,newParticleID),GlobElemID)
END IF ! TrackingMethod.EQ.REFMAPPING

IF (useDSMC.AND.(CollisMode.GT.1)) THEN
  PartStateIntEn(1,newParticleID) = VibEnergy
  PartStateIntEn(2,newParticleID) = RotEnergy
  IF (DSMC%ElectronicModel.GT.0) THEN
    PartStateIntEn(3,newParticleID) = ElecEnergy
  ENDIF
  IF (DSMC%DoAmbipolarDiff) THEN
    newAmbiParts = newAmbiParts + 1
    iPartIndx_NodeNewAmbi(newAmbiParts) = newParticleID
  END IF
END IF

PDM%ParticleInside(newParticleID)   = .TRUE.
PDM%dtFracPush(newParticleID)       = .FALSE.
PDM%IsNewPart(newParticleID)        = .TRUE.
PEM%GlobalElemID(newParticleID)     = GlobElemID
PEM%LastGlobalElemID(newParticleID) = GlobElemID

! Set particle time step and weight (if required)
IF (UseVarTimeStep) THEN
  PartTimeStep(newParticleID) = GetParticleTimeStep(PartState(1,newParticleID),PartState(2,newParticleID),PEM%LocalElemID(newParticleID))
END IF

! Set new particle MPF
IF (usevMPF) THEN
  IF(PRESENT(NewMPF))THEN
    ! MPF is already defined via input
    PartMPF(newParticleID) = NewMPF
  ELSE
    ! Check if vMPF (and radial weighting is used) to determine the MPF of the new particle
    IF (RadialWeighting%DoRadialWeighting) THEN
      PartMPF(newParticleID) = CalcRadWeightMPF(PartState(2,newParticleID),SpecID,newParticleID)
    ELSE
      PartMPF(newParticleID) = Species(SpecID)%MacroParticleFactor
    END IF
  END IF ! PRESENT(NewMPF)
END IF ! usevMPF

IF(UseRotRefFrame) THEN
  PDM%InRotRefFrame(newParticleID) = InRotRefFrameCheck(newParticleID)
  IF(PDM%InRotRefFrame(newParticleID)) THEN
    ! Initialize the velocity in the RotRefFrame by transforming the regular velocity
    PartVeloRotRef(1:3,newParticleID) = PartState(4:6,newParticleID) - CROSS(RotRefFrameOmega(1:3),PartState(1:3,newParticleID))
  ELSE
    PartVeloRotRef(1:3,newParticleID) = 0.
  END IF
END IF

IF (PRESENT(NewPartID)) NewPartID=newParticleID

END SUBROUTINE CreateParticle


SUBROUTINE RemoveParticle(PartID,BCID,alpha,crossedBC)
!===================================================================================================================================
!> Removes a single particle "PartID" by setting the required variables.
!> If CalcPartBalance/UseAdaptiveBC/CalcSurfFluxInfo = T: adds/substracts the particle to/from the respective counter
!>  !!!NOTE!!! This routine is inside particle analyze because of circular definition of modules (CalcEkinPart)
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars              ,ONLY: ElementaryCharge
USE MOD_Particle_Vars             ,ONLY: PDM, PartSpecies, Species, PartMPF, usevMPF, PartState, PartPosRef, Pt
USE MOD_Particle_Sampling_Vars    ,ONLY: UseAdaptiveBC, AdaptBCPartNumOut
USE MOD_Particle_Vars             ,ONLY: UseNeutralization, NeutralizationSource, NeutralizationBalance,nNeutralizationElems
USE MOD_Particle_Boundary_Vars    ,ONLY: PartBound
USE MOD_Particle_Analyze_Vars     ,ONLY: CalcPartBalance,nPartOut,PartEkinOut,CalcSurfFluxInfo
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: CalcBoundaryParticleOutput,BPO
USE MOD_Particle_Tracking_Vars    ,ONLY: TrackingMethod
#if defined(IMPA)
USE MOD_Particle_Vars             ,ONLY: PartIsImplicit,DoPartInNewton, PEM, PartLambdaAccept
#endif /*IMPA*/
#if defined(LSERK)
USE MOD_Particle_Vars             ,ONLY: Pt_temp
#endif
USE MOD_Particle_Analyze_Tools    ,ONLY: CalcEkinPart
USE MOD_part_tools                ,ONLY: GetParticleWeight
USE MOD_DSMC_Vars                 ,ONLY: CollInf, AmbipolElecVelo, ElectronicDistriPart, VibQuantsPar
USE MOD_Mesh_Vars                 ,ONLY: BoundaryName
#if USE_HDG
USE MOD_Globals                   ,ONLY: abort
USE MOD_HDG_Vars                  ,ONLY: UseFPC,FPC,UseEPC,EPC
USE MOD_Mesh_Vars                 ,ONLY: BoundaryType
#endif /*USE_HDG*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER, INTENT(IN)           :: PartID
INTEGER, INTENT(IN),OPTIONAL  :: BCID                    !< ID of the boundary the particle crossed
REAL, INTENT(OUT),OPTIONAL    :: alpha                   !< if removed during tracking optional alpha can be set to -1
LOGICAL, INTENT(OUT),OPTIONAL :: crossedBC               !< optional flag is needed if particle removed on BC interaction
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                       :: iSpec, iSF
REAL                          :: MPF
#if USE_HDG
INTEGER                       :: iBC,iUniqueFPCBC,iUniqueEPCBC,BCState
#endif /*USE_HDG*/
!===================================================================================================================================

! Set default values of part arrays

iSpec = PartSpecies(PartID)
! Count the number of particles per species and the kinetic energy per species
IF(CalcPartBalance) THEN
  IF(PRESENT(BCID)) THEN
    IF(PartBound%TargetBoundCond(BCID).NE.7) THEN  !skip crossing InterPlanes
      nPartOut(iSpec)=nPartOut(iSpec) + 1
      PartEkinOut(iSpec)=PartEkinOut(iSpec)+CalcEkinPart(PartID)
    END IF
  ! ELSE
  !   nPartOut(iSpec)=nPartOut(iSpec) + 1
  !   PartEkinOut(iSpec)=PartEkinOut(iSpec)+CalcEkinPart(PartID)
  END IF
END IF ! CalcPartBalance

PDM%ParticleInside(PartID) = .FALSE.
PDM%IsNewPart(PartID)      = .FALSE.
PDM%dtFracPush(PartID)     = .FALSE.
PDM%InRotRefFrame(PartID)  = .FALSE.
PartState(1:6,PartID)      = 0.
IF(TrackingMethod.EQ.REFMAPPING) PartPosRef(1:3,PartID) = -888.
PartSpecies(PartID)        = 0
Pt(1:3,PartID)             = 0.

IF(ALLOCATED(AmbipolElecVelo)) THEN
  SDEALLOCATE(AmbipolElecVelo(PartID)%ElecVelo)
END IF
IF(ALLOCATED(ElectronicDistriPart)) THEN
  SDEALLOCATE(ElectronicDistriPart(PartID)%DistriFunc)
END IF
IF(ALLOCATED(VibQuantsPar)) THEN
  SDEALLOCATE(VibQuantsPar(PartID)%Quants)
END IF

#ifdef IMPA
PartIsImplicit(PartID)   = .FALSE.
DoPartInNewton(PartID)   = .FALSE.
PartLambdaAccept(PartID) = .TRUE.
PEM%PeriodicMoved(PartID)  = .FALSE.
#endif /*IMPA*/

#if defined(LSERK)
Pt_temp(1:6,PartID)   = 0.
#endif


! If a BCID is given (e.g. when a particle is removed at a boundary), check if it is
!   - an adaptive surface flux BC or
!   - the mass flow through the boundary shall be calculated or
!   - the charges impinging on the boundary are to be summed (thruster neutralization)
IF(PRESENT(BCID)) THEN

  ! Check if adaptive BC or surface flux info
  IF(UseAdaptiveBC.OR.CalcSurfFluxInfo) THEN
    DO iSF=1,Species(iSpec)%nSurfacefluxBCs
      IF(Species(iSpec)%Surfaceflux(iSF)%BC.EQ.BCID) THEN
        Species(iSpec)%Surfaceflux(iSF)%SampledMassflow = Species(iSpec)%Surfaceflux(iSF)%SampledMassflow &
                                                          - GetParticleWeight(PartID)
        IF(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType.EQ.4)  AdaptBCPartNumOut(iSpec,iSF) = AdaptBCPartNumOut(iSpec,iSF) + 1
      END IF
    END DO
  END IF ! UseAdaptiveBC.OR.CalcSurfFluxInfo

  ! Ion thruster simulations: Landmark and Liu2010 (SPT-100) if neutralization current is determined from the particle flux over the
  ! neutralization boundary condition instead of looking into the first row of elements along that BC
  IF(UseNeutralization.AND.(nNeutralizationElems.EQ.-1))THEN
    IF(TRIM(BoundaryName(PartBound%MapToFieldBC(BCID))).EQ.TRIM(NeutralizationSource))THEN
      ! Add +1 for electrons and -X for ions: This is opposite to the summation in CountNeutralizationParticles() where the surplus
      ! of ions is calculated and compensated with an equal amount of electrons to force quasi-neutrality in the neutralization
      ! elements.
      NeutralizationBalance = NeutralizationBalance - NINT(Species(iSpec)%ChargeIC/ElementaryCharge)
    END IF
  END IF ! UseNeutralization

  ! Check if BPO boundary is encountered
  IF(CalcBoundaryParticleOutput)THEN
    IF(usevMPF)THEN
      MPF = PartMPF(PartID)
    ELSE
      MPF = Species(iSpec)%MacroParticleFactor
    END IF
    ASSOCIATE( iBPOBC   => BPO%BCIDToBPOBCID(BCID),&
               iBPOSpec => BPO%SpecIDToBPOSpecID(iSpec))
      IF(iBPOBC.GT.0.AND.iBPOSpec.GT.0)THEN! count this species on this BC
        BPO%RealPartOut(iBPOBC,iBPOSpec) = BPO%RealPartOut(iBPOBC,iBPOSpec) + MPF
      END IF ! iBPOBC.GT.0.AND.iBPOSpec.GT.0
    END ASSOCIATE
  END IF ! CalcBoundaryParticleOutput

#if USE_HDG
  ! Check if floating boundary conditions (FPC) are used
  IF(UseFPC)THEN
    iBC = PartBound%MapToFieldBC(BCID)
    IF(iBC.LE.0) CALL abort(__STAMP__,'iBC = PartBound%MapToFieldBC(BCID) must be >0',IntInfoOpt=iBC)
    IF(BoundaryType(iBC,BC_TYPE).EQ.20)THEN ! BCType = BoundaryType(iBC,BC_TYPE)
      IF(usevMPF)THEN
        MPF = PartMPF(PartID)
      ELSE
        MPF = Species(iSpec)%MacroParticleFactor
      END IF
      BCState = BoundaryType(iBC,BC_STATE) ! State is iFPC
      iUniqueFPCBC = FPC%Group(BCState,2)
      FPC%ChargeProc(iUniqueFPCBC) = FPC%ChargeProc(iUniqueFPCBC) + Species(iSpec)%ChargeIC * MPF
    END IF ! BCType.EQ.20
  END IF ! UseFPC

  ! Check if electric potential condition (EPC) is used
  IF(UseEPC)THEN
    iBC = PartBound%MapToFieldBC(BCID)
    IF(iBC.LE.0) CALL abort(__STAMP__,'iBC = PartBound%MapToFieldBC(BCID) must be >0',IntInfoOpt=iBC)
    IF(BoundaryType(iBC,BC_TYPE).EQ.8)THEN ! BCType = BoundaryType(iBC,BC_TYPE)
      IF(usevMPF)THEN
        MPF = PartMPF(PartID)
      ELSE
        MPF = Species(iSpec)%MacroParticleFactor
      END IF
      BCState = BoundaryType(iBC,BC_STATE) ! State is iEPC
      iUniqueEPCBC = EPC%Group(BCState,2)
      EPC%ChargeProc(iUniqueEPCBC) = EPC%ChargeProc(iUniqueEPCBC) + Species(iSpec)%ChargeIC * MPF
    END IF ! BCType.EQ.8
  END IF ! UseEPC
#endif /*USE_HDG*/
END IF ! PRESENT(BCID)

! Tracking-relevant variables (not required if a particle is removed within the domain, e.g. removal due to radial weighting)
IF (PRESENT(alpha)) alpha=-1.
IF (PRESENT(crossedBC)) crossedBC=.TRUE.

! DSMC: Delete the old collision partner index
IF (CollInf%ProhibitDoubleColl) CollInf%OldCollPartner(PartID) = 0

END SUBROUTINE RemoveParticle


SUBROUTINE RemoveAllElectrons()
!===================================================================================================================================
!> Removes all particles for which PARTISELECTRON(iPart) is true, i.e., all electron species
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars ,ONLY: PDM
#if USE_HDG
USE MOD_HDG_Vars      ,ONLY: BRElectronsRemoved
#endif /*USE_HDG*/
USE MOD_TimeDisc_Vars ,ONLY: time
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iPart,NbrOfElectronsRemoved
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A,ES25.14E3,A)')' RemoveAllElectrons(): Using BR electron fluid at t=',time,', removing all electrons.'

NbrOfElectronsRemoved = 0
DO iPart = 1,PDM%ParticleVecLength
  IF(.NOT.PDM%ParticleInside(iPart)) CYCLE
  IF(PARTISELECTRON(iPart))THEN
    CALL RemoveParticle(iPart)
    NbrOfElectronsRemoved = NbrOfElectronsRemoved + 1
  END IF
END DO

#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,NbrOfElectronsRemoved,1,MPI_INTEGER,MPI_SUM,MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/

IF(NbrOfElectronsRemoved.GT.0.AND.MPIRoot) WRITE(UNIT_StdOut,'(A,I0,A)') '  Removed a total of ',NbrOfElectronsRemoved,' electrons.'

#if USE_HDG
IF(NbrOfElectronsRemoved.GT.0)THEN
  BRElectronsRemoved=.TRUE.
ELSE
  BRElectronsRemoved=.FALSE.
END IF ! BRNbrOfElectronsRemove.GT.0
#endif /*USE_HDG*/

END SUBROUTINE RemoveAllElectrons

END MODULE MOD_part_operations
