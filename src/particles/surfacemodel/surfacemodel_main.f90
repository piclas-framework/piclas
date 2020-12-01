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

MODULE MOD_SurfaceModel
!===================================================================================================================================
!> Main Module for surface model
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
INTERFACE ReactiveSurfaceTreatment
  MODULE PROCEDURE ReactiveSurfaceTreatment
END INTERFACE

PUBLIC :: ReactiveSurfaceTreatment
!===================================================================================================================================

CONTAINS


SUBROUTINE ReactiveSurfaceTreatment(PartTrajectory,LengthPartTrajectory,alpha,xi,eta,PartID,SideID,n_Loc,IsSpeciesSwap,&
                              ReflectionIndex)
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
USE MOD_Particle_Boundary_Tools ,ONLY: SurfaceToPartEnergyInternal, CalcWallSample, AnalyzeSurfaceCollisions
USE MOD_Particle_Boundary_Tools ,ONLY: AddPartInfoToSample,CalcRotWallVelo
USE MOD_Particle_Boundary_Vars  ,ONLY: dXiEQ_SurfSample, Partbound, CalcSurfaceImpact, GlobalSide2SurfSide
USE MOD_TimeDisc_Vars           ,ONLY: dt, RKdtFrac
USE MOD_Particle_Surfaces       ,ONLY: CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_SurfaceModel_Vars       ,ONLY: SurfModel
USE MOD_SEE                     ,ONLY: SecondaryElectronEmission
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(INOUT)       :: ReflectionIndex !< has to be set to 1: diffuse , 2: perfect reflection, or 3: reaction
REAL,INTENT(INOUT)          :: PartTrajectory(1:3), LengthPartTrajectory, alpha
REAL,INTENT(IN)             :: xi, eta
REAL,INTENT(IN)             :: n_loc(1:3)
INTEGER,INTENT(IN)          :: PartID
INTEGER,INTENT(IN)          :: SideID
LOGICAL,INTENT(IN)          :: IsSpeciesSwap
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: ProductSpec(2)   !< 1: product species of incident particle (also used for simple reflection)
                                                     !< 2: additional species added or removed from surface
                                                          !< If productSpec is negative, then the respective particles are adsorbed
                                                          !< If productSpec is positive the particle is reflected/emitted
                                                          !< with respective species
INTEGER                          :: ProductSpecNbr   !< number of emitted particles for ProductSpec(1)
CHARACTER(30)                    :: velocityDistribution(2)   !< specifying keyword for velocity distribution
REAL                             :: TempErgy(2)               !< temperature, energy or velocity used for VeloFromDistribution
REAL                             :: reactionEnthalpy     !< negative: transferred to surface / positive: transferred from surface
LOGICAL                          :: SampledEnthalpy
REAL                             :: PartTrajectory2(1:3)
INTEGER                          :: NewPartID
REAL                             :: RanNum
REAL                             :: Xitild,EtaTild
INTEGER                          :: p,q
REAL                             :: tang1(1:3),tang2(1:3)
INTEGER                          :: SurfSideID, SpecID
! variables for Energy sampling
REAL                             :: TransArray(1:6),IntArray(1:6)
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

! =============================
! Workflow:
!
!  1.  Initial surface checks:  Check incident velocity vector and surface normal
!  2.  Select surface model:    Determine what happens at the surface
!  3.  (New) Particle handling: Perform reflection/removal of incident particle, create (multiple) possible new particles
!==============================

!===================================================================================================================================
! 1.) Initial surface checks
! find normal vector two perpendicular tangential vectors (normal_vector points outwards !!!)
!===================================================================================================================================
CALL OrthoNormVec(n_loc,tang1,tang2)

! additional states
locBCID=PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
! get BC values
WallVelo     = PartBound%WallVelo(1:3,locBCID)
WallTemp     = PartBound%WallTemp(locBCID)

POI_vec(1:3) = LastPartPos(1:3,PartID) + PartTrajectory(1:3)*alpha

IF(PartBound%RotVelo(locBCID)) THEN
  CALL CalcRotWallVelo(locBCID,PartID,POI_vec,WallVelo)
END IF

! initialize sampling arrays
TransArray(:) = 0.0
IntArray(:) = 0.0

! compute p and q
! correction of xi and eta, can only be applied if xi & eta are not used later!
IF (TriaTracking) THEN
  p=1 ; q=1
ELSE
  Xitild =MIN(MAX(-1.,xi ),0.99)
  Etatild=MIN(MAX(-1.,eta),0.99)
  p=INT((Xitild +1.0)/dXiEQ_SurfSample)+1
  q=INT((Etatild+1.0)/dXiEQ_SurfSample)+1
END IF

SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID) !SurfMesh%SideIDToSurfID(SideID)
SpecID = PartSpecies(PartID)
! Update wallcollision counter
SurfModel%Info(SpecID)%WallCollCount = SurfModel%Info(SpecID)%WallCollCount + 1

ReflectionIndex = -1 ! has to be reset in SurfaceModel, otherwise abort() will be called
reactionEnthalpy = 0.
SampledEnthalpy = .FALSE.
ProductSpec(1) = SpecID
ProductSpec(2) = 0
ProductSpecNbr = 0
velocityDistribution(1:2)=''
TempErgy(1:2)=WallTemp

!===================================================================================================================================
! 2.) Select surface model
! Here, the surfacemodel decides how the particle is treated on the surface
!===================================================================================================================================
SELECT CASE(PartBound%SurfaceModel(locBCID))
!-----------------------------------------------------------------------------------------------------------------------------------
CASE (5,6,7) ! Copied from CASE(1) and adjusted for secondary e- emission (SEE)
             ! 5: SEE by Levko2015
             ! 6: SEE by Pagonakis2016 (originally from Harrower1956)
             ! 7: SEE-I (bombarding electrons are removed, Ar+ on different materials is considered for SEE)
!-----------------------------------------------------------------------------------------------------------------------------------
  ! Get electron emission probability
  CALL SecondaryElectronEmission(PartBound%SurfaceModel(locBCID),PartID,locBCID,ReflectionIndex,ProductSpec,&
  ProductSpecNbr,TempErgy(2),velocityDistribution)
CASE DEFAULT
    CALL abort(&
      __STAMP__&
      ,'Not implemented anymore!')
END SELECT



!===================================================================================================================================
! 3.) (New) Particle handling
! Here, the incident particle is reflected/adsorbed and an additional product is emitted/adsorbed
!===================================================================================================================================
SELECT CASE(ReflectionIndex)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(1,2) ! (particle is treated in boundary condition)
!-----------------------------------------------------------------------------------------------------------------------------------
  ! 1: Diffuse reflection
  ! 2: Perfect elastic scattering
  RETURN
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(3) ! reactive interaction case
!-----------------------------------------------------------------------------------------------------------------------------------
  ! Old particle
  IF (ProductSpec(1).LT.0) THEN
    SurfModel%Info(SpecID)%NumOfAds = SurfModel%Info(SpecID)%NumOfAds + 1
  END IF

  ! New particle
  IF (ProductSpec(2).LT.0) THEN
    SurfModel%Info(SpecID)%NumOfAds = SurfModel%Info(SpecID)%NumOfAds + 1
  END IF

  !-----------------------------------------------------------
  ! Treat incident particle
  CALL AddPartInfoToSample(PartID,TransArray,IntArray,'old')
  ! Sample momentum, heatflux and collision counter on surface
  CALL CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,IsSpeciesSwap,&
                      impact_opt=CalcSurfaceImpact,PartTrajectory_opt=PartTrajectory,SurfaceNormal_opt=n_loc)
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
    CALL SurfaceToPartEnergyInternal(PartID,WallTemp)
    CALL AddPartInfoToSample(PartID,TransArray,IntArray,'new')
    ! Sample momentum, heatflux and collision counter on surface
    CALL CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,IsSpeciesSwap,emission_opt=.TRUE.)
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

      PartTrajectory2=UNITVECTOR(NewVelo(1:3))

      CALL CreateParticle(ProductSpec(2),LastPartPos(1:3,PartID),PEM%GlobalElemID(PartID),NewVelo(1:3),0.,0.,0.,NewPartID=NewPartID)
      ! Adding the energy that is transferred from the surface onto the internal energies of the particle
      CALL SurfaceToPartEnergyInternal(NewPartID,WallTemp)

      CALL AddPartInfoToSample(NewPartID,TransArray,IntArray,'new')
      CALL CalcWallSample(NewPartID,SurfSideID,p,q,Transarray,IntArray,IsSpeciesSwap,emission_opt=.TRUE.)
    END DO ! iNewPart = 1, ProductSpecNbr
  END IF
!-----------------------------------------------------------------------------------------------------------------------------------
CASE DEFAULT
!-----------------------------------------------------------------------------------------------------------------------------------
  ReflectionIndex = -1
END SELECT

END SUBROUTINE ReactiveSurfaceTreatment


END MODULE MOD_SurfaceModel
