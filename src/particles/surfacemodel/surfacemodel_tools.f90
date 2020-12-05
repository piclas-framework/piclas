!==================================================================================================================================
! Copyright (c) 2015 - 2019 Wladimir Reschke
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
PUBLIC :: SurfaceModel_ParticleEmission, SurfaceModel_EnergyAccommodation
!===================================================================================================================================

CONTAINS

SUBROUTINE SurfaceModel_ParticleEmission(PartTrajectory, LengthPartTrajectory, alpha, xi, eta, n_loc, PartID, SideID, &
            ProductSpec, ProductSpecNbr, velocityDistribution, TempErgy, &
            IsSpeciesSwap)
!===================================================================================================================================
!> Routine for Selection of Surface interaction
!===================================================================================================================================
USE MOD_Globals                 ,ONLY: abort,UNITVECTOR,OrthoNormVec
USE MOD_Globals_Vars            ,ONLY: PI
USE MOD_Particle_Tracking_Vars  ,ONLY: TriaTracking
USE MOD_Part_Tools              ,ONLY: VeloFromDistribution
USE MOD_part_operations         ,ONLY: CreateParticle, RemoveParticle
USE MOD_Particle_Vars           ,ONLY: PartState,Species,PartSpecies
USE MOD_Globals_Vars            ,ONLY: BoltzmannConst
USE MOD_Particle_Vars           ,ONLY: LastPartPos, PEM
USE MOD_Particle_Boundary_Tools ,ONLY: CalcWallSample, AnalyzeSurfaceCollisions
USE MOD_Particle_Boundary_Tools ,ONLY: CalcRotWallVelo
USE MOD_Particle_Boundary_Vars  ,ONLY: dXiEQ_SurfSample, Partbound, GlobalSide2SurfSide
USE MOD_TimeDisc_Vars           ,ONLY: dt, RKdtFrac
USE MOD_Particle_Surfaces       ,ONLY: CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_SurfaceModel_Vars       ,ONLY: SurfModel
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT-OUTPUT VARIABLES
REAL,INTENT(INOUT)          :: PartTrajectory(1:3), LengthPartTrajectory, alpha
LOGICAL,INTENT(INOUT)       :: IsSpeciesSwap
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)             :: xi, eta
REAL,INTENT(IN)             :: n_loc(1:3)
INTEGER,INTENT(IN)          :: PartID, SideID
INTEGER,INTENT(IN)          :: ProductSpec(2)   !< 1: product species of incident particle (also used for simple reflection)
                                                !< 2: additional species added or removed from surface
                                                !< If productSpec is negative, then the respective particles are adsorbed
                                                !< If productSpec is positive the particle is reflected/emitted
                                                !< with respective species
INTEGER,INTENT(IN)          :: ProductSpecNbr   !< number of emitted particles for ProductSpec(1)
CHARACTER(30),INTENT(IN)    :: velocityDistribution(2)   !< specifying keyword for velocity distribution
REAL,INTENT(IN)             :: TempErgy(2)               !< temperature, energy or velocity used for VeloFromDistribution
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: NewPartID
REAL                             :: RanNum
REAL                             :: Xitild,EtaTild
INTEGER                          :: p,q
REAL                             :: tang1(1:3),tang2(1:3)
INTEGER                          :: SurfSideID, SpecID
! variables for Energy sampling
REAL                             :: oldVelo(1:3)
INTEGER                          :: locBCID
REAL                             :: VeloReal, EtraOld
REAL                             :: EtraWall, EtraNew
REAL                             :: WallVelo(1:3), WallTemp
REAL                             :: TransACC!, VibACC, RotACC
! Polyatomic Molecules
REAL                             :: VeloCrad, Fak_D, NewVelo(3)
REAL                             :: Phi, Cmr, VeloCx, VeloCy, VeloCz
REAL                             :: POI_fak, TildTrajectory(3),POI_vec(3)
INTEGER                          :: iNewPart ! particle counter for newly created particles
!===================================================================================================================================
! compute p and q
! correction of xi and eta, can only be applied if xi & eta are not used later!

POI_vec(1:3) = LastPartPos(1:3,PartID) + PartTrajectory(1:3)*alpha
locBCID=PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)
SpecID = PartSpecies(PartID)
WallTemp = PartBound%WallTemp(locBCID)

IF(PartBound%RotVelo(locBCID)) THEN
  CALL CalcRotWallVelo(locBCID,PartID,POI_vec,WallVelo)
END IF

IF (TriaTracking) THEN
  p=1 ; q=1
ELSE
  Xitild =MIN(MAX(-1.,xi ),0.99)
  Etatild=MIN(MAX(-1.,eta),0.99)
  p=INT((Xitild +1.0)/dXiEQ_SurfSample)+1
  q=INT((Etatild+1.0)/dXiEQ_SurfSample)+1
END IF
! Old particle
IF (ProductSpec(1).LT.0) THEN
  SurfModel%Info(SpecID)%NumOfAds = SurfModel%Info(SpecID)%NumOfAds + 1
END IF

! New particle
IF (ProductSpec(2).LT.0) THEN
  SurfModel%Info(SpecID)%NumOfAds = SurfModel%Info(SpecID)%NumOfAds + 1
END IF
!-----------------------------------------------------------
CALL AnalyzeSurfaceCollisions(PartID,PartTrajectory,alpha,IsSpeciesSwap,locBCID)

IF (ProductSpec(1).LE.0) THEN
  CALL RemoveParticle(PartID,alpha=alpha)
ELSE
  oldVelo(1:3) = PartState(4:6,PartID)
  IF(TRIM(velocityDistribution(1)).NE.'') THEN
    ! sample new velocity for reflected particle
    NewVelo(1:3) = VeloFromDistribution(velocityDistribution(1),ProductSpec(1),TempErgy(1))
    ! important: n_loc points outwards
    PartState(4:6,PartID) = tang1(1:3)*NewVelo(1) + tang2(1:3)*NewVelo(2) - n_Loc(1:3)*NewVelo(3) + WallVelo(1:3)

    ! intersection point with surface
    LastPartPos(1:3,PartID) = POI_vec(1:3)
    ! recompute initial position and ignoring preceding reflections and trajectory between current position and recomputed position
    TildTrajectory=dt*RKdtFrac*oldVelo(1:3)
    POI_fak=1.- (lengthPartTrajectory-alpha)/SQRT(DOT_PRODUCT(TildTrajectory,TildTrajectory))
    ! travel rest of particle vector
    IF (PartBound%Resample(locBCID)) CALL RANDOM_NUMBER(POI_fak) !Resample Equilibirum Distribution
    ! recompute trajectory etc
    PartState(1:3,PartID)   = LastPartPos(1:3,PartID) + (1.0 - POI_fak) * dt*RKdtFrac * PartState(4:6,PartID)
  ELSE
    IF (PartBound%MomentumACC(locBCID).GT.0.0) THEN
      ! diffuse reflection
      TransACC   = PartBound%TransACC(locBCID)
      !VibACC     = PartBound%VibACC(locBCID)
      !RotACC     = PartBound%RotACC(locBCID)
      CALL RANDOM_NUMBER(RanNum)
      VeloCrad    = SQRT(-LOG(RanNum))
      CALL RANDOM_NUMBER(RanNum)
      VeloCz      = SQRT(-LOG(RanNum))
      Fak_D       = VeloCrad**2 + VeloCz**2
      EtraWall    = BoltzmannConst * WallTemp * Fak_D
      VeloReal    = SQRT(DOT_PRODUCT(oldVelo,oldVelo))
      EtraOld     = 0.5 * Species(PartSpecies(PartID))%MassIC * VeloReal**2
      EtraNew     = EtraOld + TransACC * (EtraWall - EtraOld)
      Cmr         = SQRT(2.0 * EtraNew / (Species(ProductSpec(1))%MassIC * Fak_D))
      CALL RANDOM_NUMBER(RanNum)
      Phi     = 2.0 * PI * RanNum
      VeloCx  = Cmr * VeloCrad * COS(Phi) ! tang1
      VeloCy  = Cmr * VeloCrad * SIN(Phi) ! tang2
      VeloCz  = Cmr * VeloCz
      NewVelo = VeloCx*tang1-tang2*VeloCy-VeloCz*n_loc
    ELSE
      ! perfect velocity reflection
      NewVelo(1:3) = oldVelo(1:3) - 2.*DOT_PRODUCT(oldVelo(1:3),n_loc)*n_loc
      ! mass changes, therefore velocity is scaled because momentum remains the same
      NewVelo(1:3) = NewVelo(1:3) * (Species(ProductSpec(1))%MassIC/Species(SpecID)%MassIC)
    END IF
    ! intersection point with surface
    LastPartPos(1:3,PartID) = POI_vec(1:3)
    ! recompute initial position and ignoring preceding reflections and trajectory between current position and recomputed position
    TildTrajectory=dt*RKdtFrac*oldVelo(1:3)
    POI_fak=1.- (lengthPartTrajectory-alpha)/SQRT(DOT_PRODUCT(TildTrajectory,TildTrajectory))
    ! travel rest of particle vector
    IF (PartBound%Resample(locBCID)) CALL RANDOM_NUMBER(POI_fak) !Resample Equilibirum Distribution
    PartState(1:3,PartID)   = LastPartPos(1:3,PartID) + (1.0 - POI_fak) * dt*RKdtFrac * NewVelo(1:3)
    !----  saving new particle velocity
    PartState(4:6,PartID)   = NewVelo(1:3) + WallVelo(1:3)
  END IF
  PartTrajectory=PartState(1:3,PartID) - LastPartPos(1:3,PartID)
  lengthPartTrajectory=SQRT(DOT_PRODUCT(PartTrajectory,PartTrajectory))
  PartTrajectory=PartTrajectory/lengthPartTrajectory
  ! set new species of reflected particle
  PartSpecies(PartID) = ProductSpec(1)
  ! Adding the energy that is transferred from the surface onto the internal energies of the particle
  CALL SurfaceModel_EnergyAccommodation(PartID,locBCID,WallTemp)
END IF
!-----------------------------------------------------------
! Create new particles
IF (ProductSpec(2).GT.0) THEN
  DO iNewPart = 1, ProductSpecNbr
    SurfModel%Info(ProductSpec(2))%NumOfDes = SurfModel%Info(ProductSpec(2))%NumOfDes + 1
    ! create new particle and assign correct energies
    ! sample newly created velocity
    NewVelo(1:3) = VeloFromDistribution(velocityDistribution(2),ProductSpec(2),TempErgy(2))
    ! Rotate velocity vector from global coordinate system into the surface local coordinates (important: n_loc points outwards)
    NewVelo(1:3) = tang1(1:3)*NewVelo(1) + tang2(1:3)*NewVelo(2) - n_Loc(1:3)*NewVelo(3) + WallVelo(1:3)
    ! Create new particle and get a free particle index
    CALL CreateParticle(ProductSpec(2),LastPartPos(1:3,PartID),PEM%GlobalElemID(PartID),NewVelo(1:3),0.,0.,0.,NewPartID=NewPartID)
    ! Adding the energy that is transferred from the surface onto the internal energies of the particle
    CALL SurfaceModel_EnergyAccommodation(NewPartID,locBCID,WallTemp)
  END DO ! iNewPart = 1, ProductSpecNbr
END IF

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
REAL                  :: VibQuantNewR
REAL                  :: VeloReal, RanNum, EtraOld, VeloCrad, Fak_D
REAL                  :: EtraWall, EtraNew
REAL                  :: WallVelo(1:3), TransACC, VibACC, RotACC, ElecACC
REAL                  :: tang1(1:3),tang2(1:3), NewVelo(3), POI_vec(1:3)
REAL                  :: ErotNew, ErotWall, EVibNew, Phi, Cmr, VeloCx, VeloCy, VeloCz
! Polyatomic Molecules
REAL, ALLOCATABLE     :: RanNumPoly(:), VibQuantNewRPoly(:)
REAL                  :: NormProb
INTEGER               :: iPolyatMole, iDOF
INTEGER, ALLOCATABLE  :: VibQuantNewPoly(:), VibQuantWallPoly(:), VibQuantTemp(:)
!-----------------------------------------------------------------------------------------------------------------------------------
SpecID    = PartSpecies(PartID)
TransACC  = PartBound%TransACC(locBCID)
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
          VibQuant     = NINT(PartStateIntEn(1,PartID)/(BoltzmannConst*SpecDSMC(SpecID)%CharaTVib) &
              - DSMC%GammaQuant)
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

    IF (DSMC%ElectronicModel) THEN
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

END MODULE MOD_SurfaceModel_Tools