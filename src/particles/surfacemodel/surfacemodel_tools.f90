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
PUBLIC :: SurfaceModel_ParticleEmission, SurfaceModel_EnergyAccommodation, GetWallTemperature, CalcPostWallCollVelo, CalcRotWallVelo
PUBLIC :: CalcWallTempGradient
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Routine for the particle emission at a surface
!===================================================================================================================================
SUBROUTINE SurfaceModel_ParticleEmission(n_loc, PartID, SideID, ProductSpec, ProductSpecNbr, TempErgy, GlobElemID,POI_vec)
! MODULES
USE MOD_Globals                   ,ONLY: OrthoNormVec
USE MOD_Part_Tools                ,ONLY: VeloFromDistribution
USE MOD_part_operations           ,ONLY: CreateParticle
USE MOD_Particle_Vars             ,ONLY: WriteMacroSurfaceValues
USE MOD_Particle_Boundary_Tools   ,ONLY: CalcWallSample
USE MOD_Particle_Boundary_Vars    ,ONLY: Partbound, GlobalSide2SurfSide
USE MOD_Particle_Mesh_Vars        ,ONLY: SideInfo_Shared
USE MOD_SurfaceModel_Vars         ,ONLY: SurfModEnergyDistribution
USE MOD_DSMC_Vars                 ,ONLY: DSMC, SamplingActive
USE MOD_Particle_Vars             ,ONLY: usevMPF,PartMPF
USE MOD_part_tools                ,ONLY: CalcRadWeightMPF
USE MOD_Particle_Mesh_Vars        ,ONLY: BoundsOfElem_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: n_loc(1:3)       !< normal vector of the surface
REAL,INTENT(IN)    :: POI_vec(1:3)     !< Point Of Intersection
INTEGER,INTENT(IN) :: PartID, SideID   !< Particle index and side index
INTEGER,INTENT(IN) :: GlobElemID       !< global element ID of the impacting particle (used for creating a new particle)
INTEGER,INTENT(IN) :: ProductSpec(2)   !< 1: product species of incident particle (also used for simple reflection)
                                       !< 2: additional species added or removed from surface
                                       !< If productSpec is negative, then the respective particles are adsorbed
                                       !< If productSpec is positive the particle is reflected/emitted
                                       !< with respective species
INTEGER,INTENT(IN) :: ProductSpecNbr   !< number of emitted particles for ProductSpec(1)
REAL,INTENT(IN)    :: TempErgy         !< temperature, energy or velocity used for VeloFromDistribution
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iNewPart, NewPartID, locBCID, SurfSideID
REAL               :: tang1(1:3), tang2(1:3), WallVelo(1:3), WallTemp, NewVelo(3), OldMPF, BoundsOfElemCenter(1:3),NewPos(1:3)
REAL,PARAMETER     :: eps=1e-6
REAL,PARAMETER     :: eps2=1.0-eps
!===================================================================================================================================
locBCID    = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)
WallTemp   = PartBound%WallTemp(locBCID)
WallVelo   = PartBound%WallVelo(1:3,locBCID)

IF(PartBound%RotVelo(locBCID)) THEN
  WallVelo(1:3) = CalcRotWallVelo(locBCID,POI_vec)
END IF

CALL OrthoNormVec(n_loc,tang1,tang2)

! Get Elem Center
BoundsOfElemCenter(1:3) = (/SUM(BoundsOfElem_Shared(1:2,1,GlobElemID)), &
                            SUM(BoundsOfElem_Shared(1:2,2,GlobElemID)), &
                            SUM(BoundsOfElem_Shared(1:2,3,GlobElemID)) /) / 2.

! Create new particles
DO iNewPart = 1, ProductSpecNbr
  ! create new particle and assign correct energies
  ! sample newly created velocity
  NewVelo(1:3) = VeloFromDistribution(SurfModEnergyDistribution(locBCID),TempErgy,iNewPart,ProductSpecNbr)
  ! Rotate velocity vector from global coordinate system into the surface local coordinates (important: n_loc points outwards)
  NewVelo(1:3) = tang1(1:3)*NewVelo(1) + tang2(1:3)*NewVelo(2) - n_Loc(1:3)*NewVelo(3) + WallVelo(1:3)
  ! Create new position by using POI and moving the particle by eps in the direction of the element center
  NewPos(1:3) = eps*BoundsOfElemCenter(1:3) + eps2*POI_vec(1:3)
  IF(usevMPF)THEN
    ! Get MPF of old particle
    OldMPF = PartMPF(PartID)
    ! New particle acquires the MPF of the impacting particle (not necessarily the MPF of the newly created particle species)
    CALL CreateParticle(ProductSpec(2),NewPos(1:3),GlobElemID,NewVelo(1:3),0.,0.,0.,NewPartID=NewPartID, NewMPF=OldMPF)
  ELSE
    ! New particle acquires the MPF of the new particle species
    CALL CreateParticle(ProductSpec(2),NewPos(1:3),GlobElemID,NewVelo(1:3),0.,0.,0.,NewPartID=NewPartID)
  END IF ! usevMPF
  ! Adding the energy that is transferred from the surface onto the internal energies of the particle
  CALL SurfaceModel_EnergyAccommodation(NewPartID,locBCID,WallTemp)
  ! Sampling of newly created particles
  IF((DSMC%CalcSurfaceVal.AND.SamplingActive).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) &
    CALL CalcWallSample(NewPartID,SurfSideID,'new',SurfaceNormal_opt=n_loc)
END DO ! iNewPart = 1, ProductSpecNbr

END SUBROUTINE SurfaceModel_ParticleEmission


SUBROUTINE SurfaceModel_EnergyAccommodation(PartID,locBCID,WallTemp)
!===================================================================================================================================
!> Energy accommodation at the surface: Particle internal energies PartStateIntEn() are sampled at surface temperature
!===================================================================================================================================
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_Particle_Vars         ,ONLY: PartSpecies
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
INTEGER               :: iPolyatMole, iDOF
INTEGER, ALLOCATABLE  :: VibQuantNewPoly(:), VibQuantWallPoly(:), VibQuantTemp(:)
!-----------------------------------------------------------------------------------------------------------------------------------
SpecID    = PartSpecies(PartID)
VibACC    = PartBound%VibACC(locBCID)
RotACC    = PartBound%RotACC(locBCID)
ElecACC   = PartBound%ElecACC(locBCID)

IF (useDSMC) THEN
  IF (CollisMode.GT.1) THEN
    IF ((SpecDSMC(SpecID)%InterID.EQ.2).OR.(SpecDSMC(SpecID)%InterID.EQ.20)) THEN
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
          ALLOCATE(RanNumPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF),VibQuantWallPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF), &
              VibQuantNewRPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF), VibQuantNewPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF), &
              VibQuantTemp(PolyatomMolDSMC(iPolyatMole)%VibDOF))
          CALL RANDOM_NUMBER(RanNumPoly)
          VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
          DO WHILE (ALL(VibQuantWallPoly.GE.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF))
            CALL RANDOM_NUMBER(RanNumPoly)
            VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
          END DO
          VibQuantNewRPoly(:) = VibQuantsPar(PartID)%Quants(:) + VibACC*(VibQuantWallPoly(:) - VibQuantsPar(PartID)%Quants(:))
          VibQuantNewPoly = INT(VibQuantNewRPoly)
          DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
            CALL RANDOM_NUMBER(RanNum)
            IF (RanNum.LT.(VibQuantNewRPoly(iDOF) - VibQuantNewPoly(iDOF))) THEN
              EvibNew = EvibNew + (VibQuantNewPoly(iDOF) + DSMC%GammaQuant + 1.0d0) &
                  * BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
              VibQuantTemp(iDOF) = VibQuantNewPoly(iDOF) + 1
            ELSE
              EvibNew = EvibNew + (VibQuantNewPoly(iDOF) + DSMC%GammaQuant) &
                  * BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
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
      IF((SpecDSMC(SpecID)%InterID.NE.4).AND.(.NOT.SpecDSMC(SpecID)%FullyIonized)) THEN
        CALL RANDOM_NUMBER(RanNum)
        IF (RanNum.LT.ElecACC) THEN
          PartStateIntEn(3,PartID) = RelaxElectronicShellWall(PartID, WallTemp)
        END IF
      END IF
    END IF
  END IF ! CollisMode > 1
END IF ! useDSMC

END SUBROUTINE SurfaceModel_EnergyAccommodation


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
