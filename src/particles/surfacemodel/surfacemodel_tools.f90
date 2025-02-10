!==================================================================================================================================
! Copyright (c) 2015 - 2019 Wladimir Reschke
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

MODULE MOD_SurfaceModel_Tools
!===================================================================================================================================
!> Routines with different application areas in the gas-surface interaction modelling
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
PUBLIC :: MaxwellScattering, PerfectReflection, DiffuseReflection
PUBLIC :: SurfaceModelParticleEmission, SurfaceModelEnergyAccommodation, GetWallTemperature, CalcPostWallCollVelo, CalcRotWallVelo
PUBLIC :: CalcWallTempGradient
!===================================================================================================================================

CONTAINS

SUBROUTINE MaxwellScattering(PartID,SideID,n_loc,SpecularReflectionOnly_opt)
!===================================================================================================================================
!> SurfaceModel = 0, classic DSMC gas-surface interaction model choosing between a perfect specular and a complete diffuse
!> reflection by comparing the given momentum accommodation coefficient (MomentumACC) with a random number
!===================================================================================================================================
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)             :: n_loc(1:3)
INTEGER,INTENT(IN)          :: PartID, SideID
LOGICAL,INTENT(IN),OPTIONAL :: SpecularReflectionOnly_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: RanNum,ACC
INTEGER :: iBC
LOGICAL :: SpecularReflectionOnly
!===================================================================================================================================
iBC = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
ACC = PartBound%MomentumACC(iBC)

! Check if optional parameter was supplied
IF (PRESENT(SpecularReflectionOnly_opt)) THEN; SpecularReflectionOnly = SpecularReflectionOnly_opt
ELSE;                                          SpecularReflectionOnly = PartBound%OnlySpecular(iBC)
END IF

IF (SpecularReflectionOnly) THEN
  CALL PerfectReflection(PartID,SideID,n_loc)
ELSE IF(PartBound%OnlyDiffuse(iBC)) THEN
  CALL DiffuseReflection(PartID,SideID,n_loc)
ELSE
  CALL RANDOM_NUMBER(RanNum)
  IF(RanNum.GE.ACC) THEN
    CALL PerfectReflection(PartID,SideID,n_loc)
  ELSE
    CALL DiffuseReflection(PartID,SideID,n_loc)
  END IF
END IF

END SUBROUTINE MaxwellScattering


SUBROUTINE PerfectReflection(PartID,SideID,n_Loc,opt_Symmetry)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the perfect reflection in 3D
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound
USE MOD_Particle_Vars           ,ONLY: PartState,LastPartPos,PartSpecies,Species,PartLorentzType,UseRotRefSubCycling,nSubCyclingSteps
USE MOD_DSMC_Vars               ,ONLY: DSMC, AmbipolElecVelo
USE MOD_Globals_Vars            ,ONLY: c2_inv
#if defined(LSERK)
USE MOD_Particle_Vars           ,ONLY: Pt_temp
#endif
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackInfo
USE MOD_Particle_Vars           ,ONLY: UseVarTimeStep, PartTimeStep, VarTimeStep
USE MOD_TimeDisc_Vars           ,ONLY: dt,RKdtFrac
USE MOD_Particle_Vars           ,ONLY: PDM, UseRotRefFrame, InRotRefFrame, PartVeloRotRef, RotRefFrameOmega
USE MOD_part_RHS                ,ONLY: CalcPartRHSRotRefFrame
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(IN)                   :: n_loc(1:3)
INTEGER,INTENT(IN)                :: PartID, SideID !,ElemID
LOGICAL,INTENT(IN),OPTIONAL       :: opt_Symmetry
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: WallVelo(3), v_old_Ambi(1:3), NewVeloPush(1:3), OldVelo(1:3)
REAL                                 :: LorentzFac, LorentzFacInv, POI_fak
INTEGER                              :: locBCID, SpecID
LOGICAL                              :: Symmetry
REAL                                 :: POI_vec(1:3)
REAL                                 :: NormNewVeloPush(1:3)
REAL                                 :: dtVar
!===================================================================================================================================
! Initialize
Symmetry = .FALSE.

locBCID   = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
SpecID    = PartSpecies(PartID)
WallVelo  = PartBound%WallVelo(1:3,locBCID)
IF(PRESENT(opt_Symmetry)) Symmetry = opt_Symmetry

! Get Point Of Intersection
POI_vec(1:3) = LastPartPos(1:3,PartID) + TrackInfo%PartTrajectory(1:3)*TrackInfo%alpha

IF(PartBound%RotVelo(locBCID)) THEN
  WallVelo(1:3) = CalcRotWallVelo(locBCID,POI_vec)
END IF

! Set the time step, considering whether a variable particle time step or Runge-Kutta time discretization is used
IF (UseVarTimeStep) THEN
  dtVar = dt*RKdtFrac*PartTimeStep(PartID)
ELSE
  dtVar = dt*RKdtFrac
END IF
! Species-specific time step
IF(VarTimeStep%UseSpeciesSpecific) dtVar = dtVar * Species(SpecID)%TimeStepFactor
IF(UseRotRefSubCycling) dtVar = dtVar / REAL(nSubCyclingSteps)

OldVelo = PartState(4:6,PartID)
IF(UseRotRefFrame) THEN
  ! In case of RotRefFrame utilize the respective velocity
  IF(InRotRefFrame(PartID)) OldVelo = PartVeloRotRef(1:3,PartID)
END IF

IF(SUM(ABS(WallVelo)).GT.0.)THEN
  SELECT CASE(PartLorentzType)
  CASE(3)
    PartState(4:6,PartID) = PartState(4:6,PartID) - 2.*DOT_PRODUCT(PartState(4:6,PartID),n_loc)*n_loc + WallVelo
    ! sanity check of new particle velocity
    LorentzFac=1.0-DOT_PRODUCT(PartState(4:6,PartID),PartState(4:6,PartID))*c2_inv
    IF(LorentzFac.LT.0.) CALL Abort(__STAMP__,'Particle exceeds speed of light! PartID ',PartID)
  CASE(5)
    ! map relativistic momentum to velocity
    LorentzFacInv         = 1.0+DOT_PRODUCT(PartState(4:6,PartID),PartState(4:6,PartID))*c2_inv
    LorentzFacInv         = 1.0/SQRT(LorentzFacInv)
    PartState(4:6,PartID) = LorentzFacInv*PartState(4:6,PartID)
    ! update velocity
    PartState(4:6,PartID) = PartState(4:6,PartID) - 2.*DOT_PRODUCT(PartState(4:6,PartID),n_loc)*n_loc + WallVelo
    ! map back from velocity to relativistic momentum
    LorentzFac=1.0-DOT_PRODUCT(PartState(4:6,PartID),PartState(4:6,PartID))*c2_inv
    IF(LorentzFac.LT.0.) CALL Abort(__STAMP__,'Particle exceeds speed of light! PartID ',PartID)
    LorentzFac=1.0/SQRT(LorentzFac)
    PartState(4:6,PartID) = LorentzFac*PartState(4:6,PartID)
  CASE DEFAULT
      PartState(4:6,PartID) = PartState(4:6,PartID) - 2.*DOT_PRODUCT(PartState(4:6,PartID),n_loc)*n_loc + WallVelo
  END SELECT
ELSE
  PartState(4:6,PartID) = PartState(4:6,PartID) - 2.*DOT_PRODUCT(PartState(4:6,PartID),n_loc)*n_loc
  IF (DSMC%DoAmbipolarDiff) THEN
    IF(Species(PartSpecies(PartID))%ChargeIC.GT.0.0) THEN
      v_old_Ambi = AmbipolElecVelo(PartID)%ElecVelo(1:3)
      AmbipolElecVelo(PartID)%ElecVelo(1:3) = AmbipolElecVelo(PartID)%ElecVelo(1:3) &
                     - 2.*DOT_PRODUCT(AmbipolElecVelo(PartID)%ElecVelo(1:3),n_loc)*n_loc
    END IF
  END IF
END IF

! Considering energy loss due to deformation in granular particles
IF(Species(PartSpecies(PartID))%InterID.EQ.100) PartState(4:6,PartID) = PartState(4:6,PartID) * (1.0 - PartBound%DeformEnergyLoss(locBCID))

! Set particle position on face
LastPartPos(1:3,PartID) = POI_vec(1:3)
TrackInfo%PartTrajectory(1:3)     = TrackInfo%PartTrajectory(1:3)-2.*DOT_PRODUCT(TrackInfo%PartTrajectory(1:3),n_loc)*n_loc
! Mirror the LastPartPos for new particle position
PartState(1:3,PartID) = LastPartPos(1:3,PartID) + TrackInfo%PartTrajectory(1:3)*(TrackInfo%lengthPartTrajectory - TrackInfo%alpha)

IF(UseRotRefFrame) THEN
  ! Check if rotational frame of reference is used, otherwise mirror the LastPartPos for new particle position
  IF(InRotRefFrame(PartID)) THEN
    POI_fak = 1.- (TrackInfo%lengthPartTrajectory-TrackInfo%alpha) / VECNORM(OldVelo*dtVar)
    ! Determine the correct velocity for the subsequent push in case of a rotational frame of reference at POI
    NewVeloPush(1:3) = PartState(4:6,PartID)
    NewVeloPush(1:3) = NewVeloPush(1:3) - CROSS(RotRefFrameOmega(1:3),LastPartPos(1:3,PartID))
    NewVeloPush(1:3) = NewVeloPush(1:3) + CalcPartRHSRotRefFrame(LastPartPos(1:3,PartID),NewVeloPush(1:3)) * (1.0 - POI_fak) * dtVar
      ! Make sure the NewVeloPush is pointing away from the wall
    IF(DOT_PRODUCT(n_loc,NewVeloPush(1:3)).GT.0.) THEN
      ! Normal component of new velo push v = (v dot n / |n|^2) * n, |n| = 1
      NormNewVeloPush(1:3) = DOT_PRODUCT(n_loc,NewVeloPush(1:3)) * n_loc
      ! Nullyfy normal component and keeping rest of NewVeloPush
      NewVeloPush(1:3) = NewVeloPush(1:3) - NormNewVeloPush(1:3)
      ! Move particle a little bit into the domain to avoid losing particles
      NewVeloPush(1:3) = NewVeloPush(1:3) - 1E-6 * n_loc
    END IF
    ! Store the new rotational reference frame velocity
    PartVeloRotRef(1:3,PartID) = NewVeloPush(1:3)
    ! Calc new particle position with NewVeloPush
    PartState(1:3,PartID)   = LastPartPos(1:3,PartID) + (1.0 - POI_fak) * dtVar * NewVeloPush(1:3)
  END IF
END IF

! compute moved particle || rest of movement
TrackInfo%PartTrajectory=PartState(1:3,PartID) - LastPartPos(1:3,PartID)

TrackInfo%lengthPartTrajectory = VECNORM(TrackInfo%PartTrajectory)
IF(ALMOSTZERO(TrackInfo%lengthPartTrajectory)) THEN
  TrackInfo%lengthPartTrajectory= 0.0
ELSE
  TrackInfo%PartTrajectory=TrackInfo%PartTrajectory/TrackInfo%lengthPartTrajectory
END IF
! #endif

#if defined(LSERK) || (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)
!#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506)
   ! correction for Runge-Kutta (correct position!!)
!---------- old ----------
!  absPt_temp=SQRT(Pt_temp(1,PartID)*Pt_temp(1,PartID)+Pt_temp(2,PartID)*Pt_temp(2,PartID)+Pt_temp(3,PartID)*Pt_temp(3,PartID))
!  ! scale PartTrajectory to new Pt_temp
!  Pt_temp(1:3,PartID)=absPt_temp*PartTrajectory(1:3)
!  ! deleate force history
!  Pt_temp(4:6,PartID)=0.
!  ! what happens with force term || acceleration?
!-------------------------
IF (.NOT.ALMOSTZERO(DOT_PRODUCT(WallVelo,WallVelo))) THEN
  PDM%IsNewPart(PartID)=.TRUE. !reconstruction in timedisc during push
#if defined(LSERK)
ELSE
  Pt_temp(1:3,PartID)=Pt_temp(1:3,PartID)-2.*DOT_PRODUCT(Pt_temp(1:3,PartID),n_loc)*n_loc
  IF (Symmetry) THEN !reflect also force history for symmetry
    Pt_temp(4:6,PartID)=Pt_temp(4:6,PartID)-2.*DOT_PRODUCT(Pt_temp(4:6,PartID),n_loc)*n_loc
  ELSE
    Pt_temp(4:6,PartID)=0. !produces best result compared to analytical solution in plate capacitor...
  END IF
#endif  /*LSERK*/
END IF
#endif  /*LSERK || (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)*/

END SUBROUTINE PerfectReflection


SUBROUTINE DiffuseReflection(PartID,SideID,n_loc)
!----------------------------------------------------------------------------------------------------------------------------------!
!> Computes the new particle state (position, velocity, and energy) after a diffuse reflection
!> 1.) Get the wall velocity, temperature and accommodation coefficients
!> 2.) Get the tangential vectors
!> 3.) Calculate new velocity vector (Extended Maxwellian Model)+
!> 4.) Perform vector transformation from the local to the global coordinate system and add wall velocity
!> 5.) Perform internal energy accommodation at the wall
!> 6.) Determine the new particle position after the reflection
!> 7.) Axisymmetric simulation: Rotate the vector back into the symmetry plane
!> 8.) Saving new particle velocity and recompute the trajectory based on new and old particle position
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: TwoepsMach
USE MOD_Particle_Mesh_Vars
USE MOD_DSMC_Vars               ,ONLY: DSMC, AmbipolElecVelo
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound
USE MOD_Particle_Vars           ,ONLY: PartState,LastPartPos,Species,PartSpecies
USE MOD_Particle_Vars           ,ONLY: UseRotRefFrame,InRotRefFrame,PartVeloRotRef
USE MOD_Particle_Vars           ,ONLY: UseVarTimeStep, PartTimeStep, VarTimeStep
USE MOD_TimeDisc_Vars           ,ONLY: dt,RKdtFrac
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
USE MOD_Particle_Vars           ,ONLY: PDM, RotRefFrameOmega,UseRotRefSubCycling,nSubCyclingSteps
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackInfo
USE MOD_part_RHS                ,ONLY: CalcPartRHSRotRefFrame
USE MOD_Symmetry_Vars           ,ONLY: Symmetry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(IN)                   :: n_loc(1:3)
INTEGER,INTENT(IN)                :: PartID, SideID
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: LocSideID, CNElemID, locBCID, SpecID
REAL                              :: WallVelo(1:3), WallTemp, TransACC, VibACC, RotACC, ElecACC
REAL                              :: tang1(1:3), tang2(1:3), NewVelo(3), POI_vec(1:3), NewVeloAmbi(3), VeloC(1:3), VeloCAmbi(1:3)
REAL                              :: POI_fak, TildTrajectory(3), dtVar
! Symmetry
REAL                              :: rotVelY, rotVelZ, rotPosY
REAL                              :: nx, ny, nVal, VelX, VelY, VecX, VecY, Vector1(1:2), OldVelo(1:3)
REAL                              :: NewVeloPush(1:3)
!===================================================================================================================================
! 1.) Get the wall velocity, temperature and accommodation coefficients
locBCID=PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
! get BC values
WallVelo   = PartBound%WallVelo(1:3,locBCID)
WallTemp   = GetWallTemperature(PartID,locBCID,SideID)
TransACC   = PartBound%TransACC(locBCID)
VibACC     = PartBound%VibACC(locBCID)
RotACC     = PartBound%RotACC(locBCID)
ElecACC    = PartBound%ElecACC(locBCID)

SpecID = PartSpecies(PartID)

POI_vec(1:3) = LastPartPos(1:3,PartID) + TrackInfo%PartTrajectory(1:3)*TrackInfo%alpha

IF(PartBound%RotVelo(locBCID)) THEN
  WallVelo(1:3) = CalcRotWallVelo(locBCID,POI_vec)
END IF

! Set the time step, considering whether a variable particle time step or Runge-Kutta time discretization is used
IF (UseVarTimeStep) THEN
  dtVar = dt*RKdtFrac*PartTimeStep(PartID)
ELSE
  dtVar = dt*RKdtFrac
END IF

! Species-specific time step
IF(VarTimeStep%UseSpeciesSpecific) dtVar = dtVar * Species(SpecID)%TimeStepFactor
IF(UseRotRefSubCycling) dtVar = dtVar / REAL(nSubCyclingSteps)

OldVelo = PartState(4:6,PartID)
IF(UseRotRefFrame) THEN
  ! In case of RotRefFrame utilize the respective velocity
  IF(InRotRefFrame(PartID)) OldVelo = PartVeloRotRef(1:3,PartID)
END IF

! 2.) Get the tangential vectors
IF(Symmetry%Axisymmetric) THEN
  ! Storing the old and the new particle position (which is outside the domain), at this point the position is only in the xy-plane
  VelX = PartState(1,PartID) - LastPartPos(1,PartID)
  VelY = PartState(2,PartID) - LastPartPos(2,PartID)

  CNElemID = GetCNElemID(SideInfo_Shared(SIDE_ELEMID,SideID))
  LocSideID = SideInfo_Shared(SIDE_LOCALID,SideID)

  ! Getting the vectors, which span the cell
   Vector1(1:2) = NodeCoords_Shared(1:2,ElemSideNodeID2D_Shared(1,LocSideID, CNElemID))-NodeCoords_Shared(1:2,ElemSideNodeID2D_Shared(2,LocSideID, CNElemID))

  ! Cross product of the two vectors is simplified as Vector1(3) is zero
  nx = Vector1(2)
  ny = -Vector1(1)
  ! Check for the correct orientation of the normal vectors (should be inwards)
  IF ((VelX*nx+VelY*ny).GT.0) THEN
    nx = -Vector1(2)
    ny = Vector1(1)
  END IF

  nVal = SQRT(nx*nx + ny*ny)
  nx = nx/nVal
  ny = ny/nVal
ELSE
  CALL OrthoNormVec(n_loc,tang1,tang2)
END IF

! 3.) Calculate new velocity vector (Extended Maxwellian Model)
VeloC(1:3) = CalcPostWallCollVelo(SpecID,DOTPRODUCT(PartState(4:6,PartID)),WallTemp,TransACC)
IF (DSMC%DoAmbipolarDiff) THEN
  IF(Species(SpecID)%ChargeIC.GT.0.0) THEN
    VeloCAmbi(1:3) = CalcPostWallCollVelo(DSMC%AmbiDiffElecSpec,DOTPRODUCT(AmbipolElecVelo(PartID)%ElecVelo(1:3)),WallTemp,TransACC)
  END IF
END IF

! 4.) Perform vector transformation from the local to the global coordinate system and add wall velocity
!     NewVelo = VeloCx*tang1+CROSS(-n_loc,tang1)*VeloCy-VeloCz*n_loc
IF(Symmetry%Axisymmetric) THEN
  VecX = Vector1(1) / SQRT( Vector1(1)**2 + Vector1(2)**2)
  VecY = Vector1(2) / SQRT( Vector1(1)**2 + Vector1(2)**2)
  NewVelo(1) = VecX*VeloC(1) + nx*VeloC(3)
  NewVelo(2) = VecY*VeloC(1) + ny*VeloC(3)
  NewVelo(3) = VeloC(2)
ELSE
  NewVelo(1:3) = VeloC(1)*tang1(1:3)-tang2(1:3)*VeloC(2)-VeloC(3)*n_loc(1:3)
END IF

NewVelo(1:3) = NewVelo(1:3) + WallVelo(1:3)

IF (DSMC%DoAmbipolarDiff) THEN
  IF(Species(SpecID)%ChargeIC.GT.0.0) THEN
    IF(Symmetry%Axisymmetric) THEN
      NewVeloAmbi(1) = VecX*VeloCAmbi(1) + nx*VeloCAmbi(3)
      NewVeloAmbi(2) = VecY*VeloCAmbi(1) + ny*VeloCAmbi(3)
      NewVeloAmbi(3) = VeloCAmbi(2)
    ELSE
      NewVeloAmbi(1:3) = VeloCAmbi(1)*tang1(1:3)-tang2(1:3)*VeloCAmbi(2)-VeloCAmbi(3)*n_loc(1:3)
    END IF
    NewVeloAmbi(1:3) = NewVeloAmbi(1:3) + WallVelo(1:3)
  END IF
END IF

! 5.) Perform internal energy accommodation at the wall
CALL SurfaceModelEnergyAccommodation(PartID,locBCID,WallTemp)

! 6.) Determine the new particle position after the reflection
LastPartPos(1:3,PartID) = POI_vec(1:3)

! recompute initial position and ignoring preceding reflections and trajectory between current position and recomputed position
!TildPos       =PartState(1:3,PartID)-dt*RKdtFrac*PartState(4:6,PartID)
TildTrajectory = OldVelo * dtVar
! Nullify the components in 1D and 2D to calculate the correct magnitude (2D axisymmetric is not affected)
IF(Symmetry%Order.EQ.3.OR.Symmetry%Axisymmetric) THEN
  POI_fak = VECNORM(TildTrajectory)
ELSE IF(Symmetry%Order.EQ.2.AND..NOT.Symmetry%Axisymmetric) THEN
  POI_fak = VECNORM2D(TildTrajectory(1:2))
ELSE IF(Symmetry%Order.EQ.1) THEN
  POI_fak = ABS(TildTrajectory(1))
ELSE
  CALL Abort(__STAMP__,'Error in DiffuseReflection: Symmetry%Order is not properly defined!')
END IF
POI_fak = 1. - (TrackInfo%lengthPartTrajectory-TrackInfo%alpha)/POI_fak
! travel rest of particle vector
!PartState(1:3,PartID)   = LastPartPos(1:3,PartID) + (1.0 - alpha/lengthPartTrajectory) * dt*RKdtFrac * NewVelo(1:3)
IF (PartBound%Resample(locBCID)) CALL RANDOM_NUMBER(POI_fak) !Resample Equilibirum Distribution

! ! 6a.) Determine the correct velocity for the subsequent push in case of a rotational frame of reference
NewVeloPush(1:3) = NewVelo(1:3)
IF(UseRotRefFrame) THEN
  IF(InRotRefFrame(PartID)) THEN
    NewVeloPush(1:3) = NewVeloPush(1:3) - CROSS(RotRefFrameOmega(1:3),LastPartPos(1:3,PartID))
    NewVeloPush(1:3) = NewVeloPush(1:3) + CalcPartRHSRotRefFrame(LastPartPos(1:3,PartID),NewVeloPush(1:3)) * (1.0 - POI_fak) * dtVar
    ! Store the new rotational reference frame velocity
    PartVeloRotRef(1:3,PartID) = NewVeloPush(1:3)
  END IF
END IF

PartState(1:3,PartID)   = LastPartPos(1:3,PartID) + (1.0 - POI_fak) * dtVar * NewVeloPush(1:3)

! 7.) Axisymmetric simulation: Rotate the vector back into the symmetry plane
IF(Symmetry%Axisymmetric) THEN
  ! Symmetry considerations --------------------------------------------------------
  rotPosY = SQRT(PartState(2,PartID)**2 + (PartState(3,PartID))**2)
  ! Rotation: Vy' =   Vy * cos(alpha) + Vz * sin(alpha) =   Vy * y/y' + Vz * z/y'
  !           Vz' = - Vy * sin(alpha) + Vz * cos(alpha) = - Vy * z/y' + Vz * y/y'
  ! Right-hand system, using new y and z positions after tracking, position vector and velocity vector DO NOT have to
  ! coincide (as opposed to Bird 1994, p. 391, where new positions are calculated with the velocity vector)
  IF (DSMC%DoAmbipolarDiff) THEN
    IF(Species(SpecID)%ChargeIC.GT.0.0) THEN
      rotVelY = (NewVeloAmbi(2)*(PartState(2,PartID))+NewVeloAmbi(3)*PartState(3,PartID))/rotPosY
      rotVelZ = (-NewVeloAmbi(2)*PartState(3,PartID)+NewVeloAmbi(3)*(PartState(2,PartID)))/rotPosY

      NewVeloAmbi(2) = rotVelY
      NewVeloAmbi(3) = rotVelZ
    END IF
  END IF
  rotVelY = (NewVelo(2)*(PartState(2,PartID))+NewVelo(3)*PartState(3,PartID))/rotPosY
  rotVelZ = (-NewVelo(2)*PartState(3,PartID)+NewVelo(3)*(PartState(2,PartID)))/rotPosY

  PartState(2,PartID) = rotPosY
  PartState(3,PartID) = 0.0
  NewVelo(2) = rotVelY
  NewVelo(3) = rotVelZ

  ! First check ensures that the normal vector of the side is not almost entirely in the x-direction (thus avoiding ny = 0)
  IF (ABS(nx).LT.0.999) THEN
    ! If particle has been rotated in the opposite direction of the normal vector, shift last particle position slightly
    ! away from the wall and reset LastSide to check side again.
    IF (NINT(SIGN(1.,rotPosY-POI_vec(2))).NE.NINT(SIGN(1.,ny))) THEN
      LastPartPos(2, PartID) = LastPartPos(2, PartID) + SIGN(1.,ny)*TwoepsMach
      TrackInfo%LastSide = 0
    END IF
  END IF
END IF ! Symmetry%Axisymmetric

IF(Symmetry%Order.LT.3) THEN
  ! y/z-variable is set to zero for the different symmetry cases
  LastPartPos(Symmetry%Order+1:3,PartID) = 0.0
  PartState(Symmetry%Order+1:3,PartID) = 0.0
END IF

! 8.) Saving new particle velocity and recompute the trajectory based on new and old particle position
PartState(4:6,PartID)   = NewVelo(1:3)
IF (DSMC%DoAmbipolarDiff) THEN
  IF(Species(SpecID)%ChargeIC.GT.0.0) AmbipolElecVelo(PartID)%ElecVelo(1:3) = NewVeloAmbi(1:3)
END IF

! Recompute trajectory etc
IF(Symmetry%Axisymmetric) THEN
  TrackInfo%PartTrajectory(1:2)=PartState(1:2,PartID) - LastPartPos(1:2,PartID)
  TrackInfo%PartTrajectory(3) = 0.
  TrackInfo%lengthPartTrajectory=SQRT(TrackInfo%PartTrajectory(1)**2 + TrackInfo%PartTrajectory(2)**2)
ELSE
  TrackInfo%PartTrajectory=PartState(1:3,PartID) - LastPartPos(1:3,PartID)
  TrackInfo%lengthPartTrajectory=VECNORM(TrackInfo%PartTrajectory(1:3))
END IF

IF(ABS(TrackInfo%lengthPartTrajectory).GT.0.) TrackInfo%PartTrajectory=TrackInfo%PartTrajectory/TrackInfo%lengthPartTrajectory

#if defined(LSERK) || (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)
PDM%IsNewPart(PartID)=.TRUE. !reconstruction in timedisc during push
#endif

END SUBROUTINE DiffuseReflection


!===================================================================================================================================
!> Routine for the particle emission at a surface due to a single particle-wall interaction
!===================================================================================================================================
SUBROUTINE SurfaceModelParticleEmission(n_loc, PartID, SideID, ProductSpec, ProductSpecNbr, TempErgy, GlobalElemID, POI_vec, EnergyDistribution)
! MODULES
USE MOD_Globals!                   ,ONLY: OrthoNormVec
USE MOD_Globals_Vars              ,ONLY: eV2Joule
USE MOD_Part_Tools                ,ONLY: VeloFromDistribution
USE MOD_part_operations           ,ONLY: CreateParticle
USE MOD_Particle_Vars             ,ONLY: WriteMacroSurfaceValues,Species,usevMPF,PartMPF,PartState,PartSpecies
USE MOD_Particle_Boundary_Tools   ,ONLY: CalcWallSample, StoreBoundaryParticleProperties
USE MOD_Particle_Boundary_Vars    ,ONLY: Partbound, GlobalSide2SurfSide, DoBoundaryParticleOutputHDF5
USE MOD_Particle_Mesh_Vars        ,ONLY: SideInfo_Shared
USE MOD_DSMC_Vars                 ,ONLY: DSMC, SamplingActive
USE MOD_Particle_Mesh_Vars        ,ONLY: BoundsOfElem_Shared
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: SEE,CalcElectronSEE,CalcEnergyViolationSEE
USE MOD_SurfaceModel_Vars         ,ONLY: SurfModSEEFitCoeff
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)       :: n_loc(1:3)         !< normal vector of the surface
REAL,INTENT(IN)       :: POI_vec(1:3)       !< Point Of Intersection
INTEGER,INTENT(IN)    :: PartID, SideID     !< Particle index and side index
INTEGER,INTENT(IN)    :: GlobalElemID       !< global element ID of the impacting particle (used for creating a new particle)
INTEGER,INTENT(IN)    :: ProductSpec        !< emitted species index
INTEGER,INTENT(IN)    :: ProductSpecNbr     !< number of emitted particles
REAL,INTENT(IN)       :: TempErgy           !< temperature, energy or velocity used for VeloFromDistribution
CHARACTER(LEN=*),INTENT(IN)  :: EnergyDistribution !< energy distribution model used for VeloFromDistribution
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iNewPart, NewPartID, locBCID, SurfSideID, SEEBCID
REAL               :: tang1(1:3), tang2(1:3), WallVelo(1:3), WallTemp, NewVelo(3), BoundsOfElemCenter(1:3),NewPos(1:3),MPF
REAL               :: ImpactEnergy, EnergySumSEE
REAL,PARAMETER     :: eps=1e-6
REAL,PARAMETER     :: eps2=1.0-eps
!===================================================================================================================================
locBCID    = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)
WallTemp   = PartBound%WallTemp(locBCID)
WallVelo   = PartBound%WallVelo(1:3,locBCID)
EnergySumSEE = 0.

IF(PartBound%RotVelo(locBCID)) THEN
  WallVelo(1:3) = CalcRotWallVelo(locBCID,POI_vec)
END IF

CALL OrthoNormVec(n_loc,tang1,tang2)

! Get Elem Center
BoundsOfElemCenter(1:3) = (/SUM(BoundsOfElem_Shared(1:2,1,GlobalElemID)), &
                            SUM(BoundsOfElem_Shared(1:2,2,GlobalElemID)), &
                            SUM(BoundsOfElem_Shared(1:2,3,GlobalElemID)) /) / 2.

! Create new particles
DO iNewPart = 1, ProductSpecNbr
  ! create new particle and assign correct energies
  ! sample newly created velocity
  NewVelo(1:3) = VeloFromDistribution(EnergyDistribution,TempErgy,iNewPart,ProductSpecNbr,locBCID)
  ! Rotate velocity vector from global coordinate system into the surface local coordinates (important: n_loc points outwards)
  NewVelo(1:3) = tang1(1:3)*NewVelo(1) + tang2(1:3)*NewVelo(2) - n_Loc(1:3)*NewVelo(3) + WallVelo(1:3)
  ! Create new position by using POI and moving the particle by eps in the direction of the element center
  NewPos(1:3) = eps*BoundsOfElemCenter(1:3) + eps2*POI_vec(1:3)
  ! Create new particle: in case of vMPF or VarTimeStep, new particle inherits the values of the old particle
  CALL CreateParticle(ProductSpec,NewPos(1:3),GlobalElemID,GlobalElemID,NewVelo(1:3),0.,0.,0.,OldPartID=PartID,NewPartID=NewPartID)
  ! Adding the energy that is transferred from the surface onto the internal energies of the particle
  CALL SurfaceModelEnergyAccommodation(NewPartID,locBCID,WallTemp)
  ! Sampling of newly created particles
  IF((DSMC%CalcSurfaceVal.AND.SamplingActive).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) &
    CALL CalcWallSample(NewPartID,SurfSideID,'new',SurfaceNormal_opt=n_loc)
  ! Store the particle information in PartStateBoundary.h5
  IF(DoBoundaryParticleOutputHDF5) THEN
    IF(usevMPF)THEN
      MPF = PartMPF(NewPartID)
    ELSE
      MPF = Species(ProductSpec)%MacroParticleFactor
    END IF ! usevMPF
    CALL StoreBoundaryParticleProperties(NewPartID,ProductSpec,PartState(1:3,NewPartID),&
          UNITVECTOR(PartState(4:6,NewPartID)),n_loc,iPartBound=locBCID,mode=2,MPF_optIN=MPF)
  END IF ! DoBoundaryParticleOutputHDF5
  ! If more than one secondary electron, sum-up the energy to track energy conservation violations
  IF(CalcEnergyViolationSEE) THEN
    SEEBCID = SEE%BCIDToSEEBCID(locBCID)
    ! Sum-up the energy of all secondaries
    EnergySumSEE = EnergySumSEE + 0.5*Species(ProductSpec)%MassIC*DOTPRODUCT(NewVelo(1:3))
    ! Treatment at the end of the secondaries loop, energy conservation violation is only counted once per impact
    IF(iNewPart.EQ.ProductSpecNbr) THEN
      SEE%EventCount(SEEBCID) = SEE%EventCount(SEEBCID) + 1
      ! Calculated the resulting energy, which should have been distributed (impact energy minus work function)
      ImpactEnergy = 0.5*Species(PartSpecies(PartID))%MassIC*DOTPRODUCT(PartState(4:6,PartID)) - SurfModSEEFitCoeff(4,locBCID)*eV2Joule
      IF(EnergySumSEE.GT.ImpactEnergy) THEN
        ! Count the violation
        SEE%RealElectronEnergyViolationCount(SEEBCID) = SEE%RealElectronEnergyViolationCount(SEEBCID) + 1
        ! Sum-up energy addition as percentage of the impact energy
        SEE%RealElectronEnergyViolationSum(SEEBCID) = SEE%RealElectronEnergyViolationSum(SEEBCID) + (EnergySumSEE - ImpactEnergy) / ImpactEnergy
      END IF
    END IF
  END IF
END DO ! iNewPart = 1, ProductSpecNbr

END SUBROUTINE SurfaceModelParticleEmission


SUBROUTINE SurfaceModelEnergyAccommodation(PartID,locBCID,WallTemp)
!===================================================================================================================================
!> Energy accommodation at the surface: Particle internal energies PartStateIntEn() are sampled at surface temperature
!===================================================================================================================================
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_Particle_Vars         ,ONLY: PartSpecies, Species
USE MOD_Particle_Boundary_Vars,ONLY: PartBound
USE MOD_DSMC_Vars             ,ONLY: CollisMode, PolyatomMolDSMC, useDSMC
USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, SpecDSMC, DSMC, VibQuantsPar
USE MOD_DSMC_ElectronicModel  ,ONLY: RelaxElectronicShellWall
#if (PP_TimeDiscMethod==400)
USE MOD_BGK_Vars              ,ONLY: BGKDoVibRelaxation
#elif (PP_TimeDiscMethod==300)
USE MOD_FPFlow_Vars           ,ONLY: FPDoVibRelaxation
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)    :: PartID, locBCID
REAL,INTENT(IN)       :: WallTemp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: SpecID, vibQuant, vibQuantNew, VibQuantWall
REAL                  :: RanNum
REAL                  :: VibACC, RotACC, ElecACC
REAL                  :: ErotNew, ErotWall, EVibNew
! Polyatomic Molecules
REAL                  :: NormProb, VibQuantNewR
REAL, ALLOCATABLE     :: RanNumPoly(:), VibQuantNewRPoly(:)
INTEGER               :: iPolyatMole, iDOF, VibDOF
INTEGER, ALLOCATABLE  :: VibQuantNewPoly(:), VibQuantWallPoly(:), VibQuantTemp(:)
!-----------------------------------------------------------------------------------------------------------------------------------

! Requires the usage of DSMC and CollisMode = 2/3
IF(.NOT.useDSMC) RETURN
IF(.NOT.(CollisMode.GT.1)) RETURN

SpecID    = PartSpecies(PartID)
VibACC    = PartBound%VibACC(locBCID)
RotACC    = PartBound%RotACC(locBCID)
ElecACC   = PartBound%ElecACC(locBCID)

IF ((Species(SpecID)%InterID.EQ.2).OR.(Species(SpecID)%InterID.EQ.20)) THEN
  !---- Rotational energy accommodation
  IF (SpecDSMC(SpecID)%Xi_Rot.EQ.2) THEN
    CALL RANDOM_NUMBER(RanNum)
    ErotWall = - BoltzmannConst * WallTemp * LOG(RanNum)
  ELSE IF (SpecDSMC(SpecID)%Xi_Rot.EQ.3) THEN
    CALL RANDOM_NUMBER(RanNum)
    ErotWall = RanNum*10. !the distribution function has only non-negligible  values betwenn 0 and 10
    NormProb = SQRT(ErotWall)*EXP(-ErotWall)/(SQRT(0.5)*EXP(-0.5))
    CALL RANDOM_NUMBER(RanNum)
    DO WHILE (RanNum.GE.NormProb)
      CALL RANDOM_NUMBER(RanNum)
      ErotWall = RanNum*10. !the distribution function has only non-negligible  values betwenn 0 and 10
      NormProb = SQRT(ErotWall)*EXP(-ErotWall)/(SQRT(0.5)*EXP(-0.5))
      CALL RANDOM_NUMBER(RanNum)
    END DO
    ErotWall = ErotWall*BoltzmannConst*WallTemp
  END IF
  ErotNew  = PartStateIntEn(2,PartID) + RotACC *(ErotWall - PartStateIntEn(2,PartID))

  PartStateIntEn(2,PartID) = ErotNew

#if (PP_TimeDiscMethod==400)
  IF (BGKDoVibRelaxation) THEN
#elif (PP_TimeDiscMethod==300)
  IF (FPDoVibRelaxation) THEN
#endif
    !---- Vibrational energy accommodation
    IF(SpecDSMC(SpecID)%PolyatomicMol) THEN
      EvibNew = 0.0
      iPolyatMole = SpecDSMC(SpecID)%SpecToPolyArray
      VibDOF = PolyatomMolDSMC(iPolyatMole)%VibDOF
      ALLOCATE(RanNumPoly(VibDOF),VibQuantWallPoly(VibDOF),VibQuantNewRPoly(VibDOF),VibQuantNewPoly(VibDOF),VibQuantTemp(VibDOF))
      CALL RANDOM_NUMBER(RanNumPoly)
      VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
      DO WHILE (ALL(VibQuantWallPoly.GE.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF))
        CALL RANDOM_NUMBER(RanNumPoly)
        VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
      END DO
      VibQuantNewRPoly(:) = VibQuantsPar(PartID)%Quants(:) + VibACC*(VibQuantWallPoly(:) - VibQuantsPar(PartID)%Quants(:))
      VibQuantNewPoly = INT(VibQuantNewRPoly)
      DO iDOF = 1, VibDOF
        EvibNew = EvibNew + (VibQuantNewPoly(iDOF) + DSMC%GammaQuant)*BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
        CALL RANDOM_NUMBER(RanNum)
        IF (RanNum.LT.(VibQuantNewRPoly(iDOF) - VibQuantNewPoly(iDOF))) THEN
          EvibNew = EvibNew + BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
          VibQuantTemp(iDOF) = VibQuantNewPoly(iDOF) + 1
        ELSE
          VibQuantTemp(iDOF) = VibQuantNewPoly(iDOF)
        END IF
      END DO
    ELSE
      VibQuant     = NINT(PartStateIntEn(1,PartID)/(BoltzmannConst*SpecDSMC(SpecID)%CharaTVib) - DSMC%GammaQuant)
      CALL RANDOM_NUMBER(RanNum)
      VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(SpecID)%CharaTVib)
      DO WHILE (VibQuantWall.GE.SpecDSMC(SpecID)%MaxVibQuant)
        CALL RANDOM_NUMBER(RanNum)
        VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(SpecID)%CharaTVib)
      END DO
      VibQuantNewR = VibQuant + VibACC*(VibQuantWall - VibQuant)
      VibQuantNew = INT(VibQuantNewR)
      CALL RANDOM_NUMBER(RanNum)
      IF (RanNum.LT.(VibQuantNewR - VibQuantNew)) THEN
        EvibNew = (VibQuantNew + DSMC%GammaQuant + 1.0d0)*BoltzmannConst*SpecDSMC(SpecID)%CharaTVib
      ELSE
        EvibNew = (VibQuantNew + DSMC%GammaQuant)*BoltzmannConst*SpecDSMC(SpecID)%CharaTVib
      END IF
    END IF

    IF(SpecDSMC(SpecID)%PolyatomicMol) VibQuantsPar(PartID)%Quants(:) = VibQuantTemp(:)
    PartStateIntEn(1,PartID) = EvibNew
#if (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==300)
  END IF ! FPDoVibRelaxation || BGKDoVibRelaxation
#endif
END IF

IF (DSMC%ElectronicModel.GT.0) THEN
  IF((Species(SpecID)%InterID.NE.4).AND.(.NOT.SpecDSMC(SpecID)%FullyIonized).AND.(Species(SpecID)%InterID.NE.100)) THEN
    CALL RANDOM_NUMBER(RanNum)
    IF (RanNum.LT.ElecACC) THEN
      PartStateIntEn(3,PartID) = RelaxElectronicShellWall(PartID, WallTemp)
    END IF
  END IF
END IF

END SUBROUTINE SurfaceModelEnergyAccommodation


REAL FUNCTION GetWallTemperature(PartID,locBCID, SideID)
!===================================================================================================================================
!> Determine the wall temperature, current options: determine a temperature based on an imposed gradient or use a fixed temperature
!===================================================================================================================================
USE MOD_Globals                 ,ONLY: DOTPRODUCT, VECNORM
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound, BoundaryWallTemp, GlobalSide2SurfSide
USE MOD_Particle_Vars           ,ONLY: LastPartPos
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackInfo
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: locBCID, PartID, SideID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: POI(3)
!-----------------------------------------------------------------------------------------------------------------------------------
IF(PartBound%WallTemp2(locBCID).GT.0.0) THEN
  POI(1:3) = LastPartPos(1:3,PartID) + TrackInfo%PartTrajectory(1:3)*TrackInfo%alpha
  GetWallTemperature = CalcWallTempGradient(POI,locBCID)
ELSE IF (PartBound%UseAdaptedWallTemp(locBCID)) THEN
  GetWallTemperature = BoundaryWallTemp(TrackInfo%p,TrackInfo%q,GlobalSide2SurfSide(SURF_SIDEID,SideID))
ELSE
  GetWallTemperature = PartBound%WallTemp(locBCID)
END IF

END FUNCTION GetWallTemperature


PPURE REAL FUNCTION CalcWallTempGradient(PointVec,locBCID)
!===================================================================================================================================
!> Calculation of the wall temperature at a specific position due to the imposed temperature gradient (WallTemp2.GT.0)
!===================================================================================================================================
USE MOD_Globals                 ,ONLY: DOTPRODUCT, VECNORM
USE MOD_Globals_Vars            ,ONLY: EpsMach
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: PointVec(3)
INTEGER, INTENT(IN)             :: locBCID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: Bounds(1:3), TempGradLength, PointVec_projected(1:3), WallTemp2
!-----------------------------------------------------------------------------------------------------------------------------------
ASSOCIATE(PB => PartBound)
PointVec_projected(1:3) = PB%TempGradStart(1:3,locBCID) + DOT_PRODUCT((PointVec(1:3) - PB%TempGradStart(1:3,locBCID)), &
                          PB%TempGradVec(1:3,locBCID)) / DOTPRODUCT(PB%TempGradVec(1:3,locBCID)) * PB%TempGradVec(1:3,locBCID)
TempGradLength = VECNORM(PointVec_projected(1:3)-PB%TempGradStart(1:3,locBCID)) / VECNORM(PB%TempGradVec(1:3,locBCID))

SELECT CASE(PB%TempGradDir(locBCID))
CASE(0)
  ! Position is projected onto the gradient vector
  Bounds(1:3) = PointVec_projected(1:3)
  ! Wall temperature is set to the end value
  WallTemp2   = PB%WallTemp2(locBCID)
CASE(1,2,3)
  ! Simply using the actual position as bounds
  Bounds(1:3) = PointVec(1:3)
  ! Wall temperature is set to the end value
  WallTemp2   = PB%WallTemp2(locBCID)
END SELECT

IF(MINVAL(Bounds(1:3)-PB%TempGradStart(1:3,locBCID)).LT.-EpsMach) THEN
  CalcWallTempGradient = PB%WallTemp(locBCID)
ELSEIF(MINVAL(PB%TempGradEnd(1:3,locBCID)-Bounds(1:3)).LT.-EpsMach) THEN
  CalcWallTempGradient = WallTemp2
ELSE
  IF(TempGradLength.LT.0.0) THEN
    CalcWallTempGradient = PB%WallTemp(locBCID)
  ELSE IF(TempGradLength.GT.1.0) THEN
    CalcWallTempGradient = WallTemp2
  ELSE
    CalcWallTempGradient = PB%WallTemp(locBCID) + TempGradLength * PB%WallTempDelta(locBCID)
  END IF
END IF
END ASSOCIATE

END FUNCTION CalcWallTempGradient


FUNCTION CalcPostWallCollVelo(SpecID,VeloSquare,WallTemp,TransACC)
!===================================================================================================================================
!> Calculate the new velocity vector for the reflected particle based on the wall temperature and accommodation coefficient
!===================================================================================================================================
USE MOD_Globals                 ,ONLY: DOTPRODUCT
USE MOD_Globals_Vars            ,ONLY: BoltzmannConst, PI
USE MOD_Particle_Vars           ,ONLY: Species
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: SpecID
REAL, INTENT(IN)                :: VeloSquare, WallTemp, TransACC
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                            :: CalcPostWallCollVelo(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: RanNum, VeloCrad, VeloCz, Fak_D, EtraOld, EtraNew, Cmr, Phi
!-----------------------------------------------------------------------------------------------------------------------------------

EtraOld     = 0.5 * Species(SpecID)%MassIC * VeloSquare
CALL RANDOM_NUMBER(RanNum)
VeloCrad    = SQRT(-LOG(RanNum))
CALL RANDOM_NUMBER(RanNum)
VeloCz      = SQRT(-LOG(RanNum))
Fak_D       = VeloCrad**2 + VeloCz**2

EtraNew     = EtraOld + TransACC * (BoltzmannConst * WallTemp * Fak_D - EtraOld)
Cmr         = SQRT(2.0 * EtraNew / (Species(SpecID)%MassIC * Fak_D))
CALL RANDOM_NUMBER(RanNum)
Phi     = 2.0 * PI * RanNum

CalcPostWallCollVelo(1)  = Cmr * VeloCrad * COS(Phi) ! tang1
CalcPostWallCollVelo(2)  = Cmr * VeloCrad * SIN(Phi) ! tang2
CalcPostWallCollVelo(3)  = Cmr * VeloCz

END FUNCTION CalcPostWallCollVelo


PPURE FUNCTION CalcRotWallVelo(locBCID,POI)
!===================================================================================================================================
!> Calculation of additional velocity through the rotating wall. The velocity is equal to circumferential speed at
!> the point of intersection (POI):
!> The direction is perpendicular to the rotational axis (vec_axi) AND the distance vector (vec_axi -> POI).
!> Rotation direction based on Right-hand rule. The magnitude of the velocity depends on radius and rotation frequency.
!> Currently implemented: simplified version assuming that the rotational axis is one of the major axis x,y or z.
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals                 ,ONLY: CROSS
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN)    :: locBCID
REAL,INTENT(IN)       :: POI(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                  :: CalcRotWallVelo(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! Case: rotational axis is NOT one of the major axis (x,y,z)
! vec_OrgPOI(1:3) = POI(1:3) - PartBound%RotOrg(1:3,locBCID)
! vec_axi_norm = PartBound%RotAxis(1:3,locBCID) / VECNORM(PartBound%RotAxis(1:3,locBCID))
! vec_a(1:3) = DOT_PRODUCT(vec_axi_norm,vec_OrgPOI) * vec_axi_norm(1:3)
! vec_r(1:3) = vec_OrgPOI(1:3) - vec_a(1:3)
! radius = SQRT( vec_r(1)*vec_r(1) + vec_r(2)*vec_r(2) + vec_r(3)*vec_r(3) )
! circ_speed = 2.0 * PI * radius * PartBound%RotFreq(locBCID)
! vec_t = CROSSNORM(vec_axi_norm,vec_r)
! WallVelo(1:3) = circ_speed * vec_t(1:3)

! Case: rotational is one of the major axis (x,y,z)
CalcRotWallVelo(1:3) = CROSS(PartBound%RotOmega(1:3,locBCID),POI(1:3))

END FUNCTION CalcRotWallVelo


END MODULE MOD_SurfaceModel_Tools
