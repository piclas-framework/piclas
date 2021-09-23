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

MODULE MOD_part_operations
!===================================================================================================================================
! Contains tools for particle related operations. This routine is uses MOD_Particle_Boundary_Tools, but not vice versa!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE CreateParticle
  MODULE PROCEDURE CreateParticle
END INTERFACE

INTERFACE RemoveParticle
  MODULE PROCEDURE RemoveParticle
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: CreateParticle, RemoveParticle
PUBLIC :: RemoveAllElectrons
!===================================================================================================================================

CONTAINS

SUBROUTINE CreateParticle(SpecID,Pos,ElemID,Velocity,RotEnergy,VibEnergy,ElecEnergy,NewPartID)
!===================================================================================================================================
!> creates a single particle at correct array position and assign properties
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars           ,ONLY: PDM, PEM, PartState, LastPartPos, PartSpecies,PartPosRef, VarTimeStep, usevMPF, PartMPF
USE MOD_Particle_Vars           ,ONLY: Species
USE MOD_DSMC_Vars               ,ONLY: useDSMC, CollisMode, DSMC, PartStateIntEn, RadialWeighting
USE MOD_DSMC_Vars               ,ONLY: newAmbiParts, iPartIndx_NodeNewAmbi
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
USE MOD_Eval_xyz                ,ONLY: GetPositionInRefElem
USE MOD_part_tools              ,ONLY: CalcRadWeightMPF
USE MOD_Particle_VarTimeStep    ,ONLY: CalcVarTimeStep
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER, INTENT(IN)           :: SpecID        !< Species ID
REAL, INTENT(IN)              :: Pos(1:3)      !< Position (x,y,z)
INTEGER, INTENT(IN)           :: ElemID        !< global element ID
REAL, INTENT(IN)              :: Velocity(1:3) !< Velocity (vx,vy,vz)
REAL, INTENT(IN)              :: RotEnergy     !< Rotational energy
REAL, INTENT(IN)              :: VibEnergy     !< Vibrational energy
REAL, INTENT(IN)              :: ElecEnergy    !< Electronic energy
INTEGER, INTENT(OUT),OPTIONAL :: NewPartID     !< ID of newly created particle
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER :: newParticleID
!===================================================================================================================================

! Do not increase the ParticleVecLength for Phantom particles!
PDM%ParticleVecLength = PDM%ParticleVecLength + 1 ! Increase particle vector length
newParticleID = PDM%ParticleVecLength
IF(newParticleID.GT.PDM%MaxParticleNumber)THEN
  CALL abort(&
      __STAMP__&
      ,'CreateParticle: newParticleID.GT.PDM%MaxParticleNumber. newParticleID=',IntInfoOpt=newParticleID)
END IF

PartSpecies(newParticleID)     = SpecID
LastPartPos(1:3,newParticleID) = Pos(1:3)
PartState(1:3,newParticleID)   = Pos(1:3)
PartState(4:6,newParticleID)   = Velocity(1:3)

! Set the new reference position here
IF(TrackingMethod.EQ.REFMAPPING)THEN
  CALL GetPositionInRefElem(PartState(1:3,newParticleID),PartPosRef(1:3,newParticleID),ElemID)
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
PEM%GlobalElemID(newParticleID)     = ElemID
PEM%LastGlobalElemID(newParticleID) = ElemID

! Set particle time step and weight (if required)
IF (VarTimeStep%UseVariableTimeStep) THEN
  VarTimeStep%ParticleTimeStep(newParticleID) = &
    CalcVarTimeStep(PartState(1,newParticleID),PartState(2,newParticleID),PEM%LocalElemID(newParticleID))
END IF
IF (usevMPF) THEN
  IF (RadialWeighting%DoRadialWeighting) THEN
    PartMPF(newParticleID) = CalcRadWeightMPF(PartState(2,newParticleID),SpecID,newParticleID)
  ELSE
    PartMPF(newParticleID) = Species(SpecID)%MacroParticleFactor
  END IF
END IF

IF (PRESENT(NewPartID)) NewPartID=newParticleID

END SUBROUTINE CreateParticle


SUBROUTINE RemoveParticle(PartID,BCID,alpha,crossedBC)
!===================================================================================================================================
!> Removes a single particle "PartID" by setting the required variables.
!> If CalcPartBalance/UseAdaptive/CalcAdaptiveBCInfo = T: adds/substracts the particle to/from the respective counter
!>  !!!NOTE!!! This routine is inside particle analyze because of circular definition of modules (CalcEkinPart)
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars             ,ONLY: PDM, PartSpecies, Species, PartMPF, usevMPF
USE MOD_Particle_Sampling_Vars    ,ONLY: UseAdaptive, AdaptBCPartNumOut
USE MOD_Particle_Vars             ,ONLY: UseNeutralization, NeutralizationSource, NeutralizationBalance
USE MOD_Particle_Analyze_Vars     ,ONLY: CalcPartBalance,nPartOut,PartEkinOut,CalcAdaptiveBCInfo
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: CalcBoundaryParticleOutput,BPO
#if defined(IMPA)
USE MOD_Particle_Vars             ,ONLY: PartIsImplicit
USE MOD_Particle_Vars             ,ONLY: DoPartInNewton
#endif /*IMPA*/
USE MOD_Particle_Analyze_Tools    ,ONLY: CalcEkinPart
USE MOD_part_tools                ,ONLY: GetParticleWeight
USE MOD_DSMC_Vars                 ,ONLY: CollInf
USE MOD_Mesh_Vars                 ,ONLY: BoundaryName
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
!===================================================================================================================================

PDM%ParticleInside(PartID) = .FALSE.
#ifdef IMPA
PartIsImplicit(PartID) = .FALSE.
DoPartInNewton(PartID) = .FALSE.
#endif /*IMPA*/

iSpec = PartSpecies(PartID)
! Count the number of particles per species and the kinetic energy per species
IF(CalcPartBalance) THEN
  nPartOut(iSpec)=nPartOut(iSpec) + 1
  PartEkinOut(iSpec)=PartEkinOut(iSpec)+CalcEkinPart(PartID)
END IF ! CalcPartBalance

! If a BCID is given (e.g. when a particle is removed at a boundary), check if it is
!   - an adaptive surface flux BC or
!   - the mass flow through the boundary shall be calculated or
!   - the charges impinging on the boundary are to be summed (thruster neutralization)
IF(PRESENT(BCID)) THEN
  IF(UseAdaptive.OR.CalcAdaptiveBCInfo) THEN
    DO iSF=1,Species(iSpec)%nSurfacefluxBCs
      IF(Species(iSpec)%Surfaceflux(iSF)%BC.EQ.BCID) THEN
        Species(iSpec)%Surfaceflux(iSF)%SampledMassflow = Species(iSpec)%Surfaceflux(iSF)%SampledMassflow &
                                                          - GetParticleWeight(PartID)
        IF(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType.EQ.4)  AdaptBCPartNumOut(iSpec,iSF) = AdaptBCPartNumOut(iSpec,iSF) + 1
      END IF
    END DO
  END IF ! UseAdaptive.OR.CalcAdaptiveBCInfo
  IF(UseNeutralization)THEN
    IF(TRIM(BoundaryName(BCID)).EQ.TRIM(NeutralizationSource))THEN
      ! Add +1 for electrons and -1 for ions
      NeutralizationBalance = NeutralizationBalance - INT(SIGN(1.0, Species(iSpec)%ChargeIC))
    END IF
  END IF ! UseNeutralization
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
CALL MPI_ALLREDUCE(MPI_IN_PLACE,NbrOfElectronsRemoved,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,iError)
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
