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
!===================================================================================================================================

CONTAINS

SUBROUTINE CreateParticle(Species,Pos,ElemID,Velocity,RotEnergy,VibEnergy,ElecEnergy,NewPartID)
!===================================================================================================================================
!> creates a single particle at correct array position and assign properties
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: PDM, PEM, PartState, LastPartPos, PartSpecies,PartPosRef
USE MOD_DSMC_Vars              ,ONLY: useDSMC, CollisMode, DSMC, PartStateIntEn
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER, INTENT(IN)           :: Species
REAL, INTENT(IN)              :: Pos(1:3)
INTEGER, INTENT(IN)           :: ElemID
REAL, INTENT(IN)              :: Velocity(1:3)
REAL, INTENT(IN)              :: RotEnergy
REAL, INTENT(IN)              :: VibEnergy
REAL, INTENT(IN)              :: ElecEnergy
INTEGER, INTENT(OUT),OPTIONAL :: NewPartID
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

PartSpecies(newParticleID)     = Species
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
  IF (DSMC%ElectronicModel) THEN
    PartStateIntEn(3,newParticleID) = ElecEnergy
  ENDIF
END IF

PDM%ParticleInside(newParticleID) = .TRUE.
PDM%dtFracPush(newParticleID)     = .FALSE.
PDM%IsNewPart(newParticleID)      = .TRUE.
PEM%Element(newParticleID)        = ElemID
PEM%lastElement(newParticleID)    = ElemID

! ?????? necessary?
! IF (VarTimeStep%UseVariableTimeStep) THEN
!   VarTimeStep%ParticleTimeStep(newParticleID) &
!     = CalcVarTimeStep(PartState(1,newParticleID),PartState(2,newParticleID),PEM%Element(newParticleID))
! END IF
! IF (RadialWeighting%DoRadialWeighting) THEN
!   PartMPF(newParticleID) = CalcRadWeightMPF(PartState(2,newParticleID), 1,newParticleID)
! END IF
IF (PRESENT(NewPartID)) NewPartID=newParticleID

END SUBROUTINE CreateParticle


SUBROUTINE RemoveParticle(PartID,BCID,alpha,crossedBC)
!===================================================================================================================================
!> Removes a single particle "PartID" by setting the required variables.
!> If CalcPartBalance/UseAdaptive/CalcMassflowRate = T: adds/substracts the particle to/from the respective counter
!>  !!!NOTE!!! This routine is inside particle analyze because of circular definition of modules (CalcEkinPart)
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars           ,ONLY: PDM, PartSpecies, Species, UseAdaptive
USE MOD_Particle_Analyze_Vars   ,ONLY: CalcPartBalance,nPartOut,PartEkinOut,CalcMassflowRate
#if defined(IMPA)
USE MOD_Particle_Vars           ,ONLY: PartIsImplicit
USE MOD_Particle_Vars           ,ONLY: DoPartInNewton
#endif /*IMPA*/
USE MOD_Particle_Analyze_Tools  ,ONLY: CalcEkinPart
USE MOD_part_tools              ,ONLY: GetParticleWeight
USE MOD_DSMC_Vars               ,ONLY: CollInf
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

! If a BCID is given (e.g. when a particle is removed at a boundary), check if its an adaptive surface flux BC or the mass flow
! through the boundary shall be calculated
IF(PRESENT(BCID)) THEN
  IF(UseAdaptive.OR.CalcMassflowRate) THEN
    DO iSF=1,Species(iSpec)%nSurfacefluxBCs
      IF(Species(iSpec)%Surfaceflux(iSF)%BC.EQ.BCID) THEN
        Species(iSpec)%Surfaceflux(iSF)%SampledMassflow = Species(iSpec)%Surfaceflux(iSF)%SampledMassflow &
                                                          - GetParticleWeight(PartID)
        IF(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType.EQ.4)  Species(iSpec)%Surfaceflux(iSF)%AdaptivePartNumOut = &
            Species(iSpec)%Surfaceflux(iSF)%AdaptivePartNumOut + 1
      END IF
    END DO
  END IF ! UseAdaptive.OR.CalcMassflowRate
END IF ! PRESENT(BCID)

! Tracking-relevant variables (not required if a particle is removed within the domain, e.g. removal due to radial weighting)
IF (PRESENT(alpha)) alpha=-1.
IF (PRESENT(crossedBC)) crossedBC=.TRUE.

! DSMC: Delete the old collision partner index
IF (CollInf%ProhibitDoubleColl) CollInf%OldCollPartner(PartID) = 0

END SUBROUTINE RemoveParticle

END MODULE MOD_part_operations
