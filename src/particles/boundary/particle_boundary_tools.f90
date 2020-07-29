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

MODULE MOD_Particle_Boundary_Tools
!===================================================================================================================================
! Tools used for boundary interactions
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
INTERFACE AddPartInfoToSample
  MODULE PROCEDURE AddPartInfoToSample
END INTERFACE

INTERFACE CalcWallSample
  MODULE PROCEDURE CalcWallSample
END INTERFACE

INTERFACE AnalyzeSurfaceCollisions
  MODULE PROCEDURE AnalyzeSurfaceCollisions
END INTERFACE

INTERFACE SurfaceToPartEnergyInternal
  MODULE PROCEDURE SurfaceToPartEnergyInternal
END INTERFACE

INTERFACE CountSurfaceImpact
  MODULE PROCEDURE CountSurfaceImpact
END INTERFACE

INTERFACE StoreBoundaryParticleProperties
  MODULE PROCEDURE StoreBoundaryParticleProperties
END INTERFACE

INTERFACE DielectricSurfaceCharge
  MODULE PROCEDURE DielectricSurfaceCharge
END INTERFACE

PUBLIC :: AddPartInfoToSample
PUBLIC :: CalcWallSample
PUBLIC :: AnalyzeSurfaceCollisions
PUBLIC :: SurfaceToPartEnergyInternal
PUBLIC :: CountSurfaceImpact
PUBLIC :: StoreBoundaryParticleProperties
PUBLIC :: SortArray,SortArray2
PUBLIC :: DielectricSurfaceCharge
PUBLIC :: GetWallTemperature
!===================================================================================================================================

CONTAINS


SUBROUTINE AddPartInfoToSample(PartID,Transarray,IntArray,SampleType)
!===================================================================================================================================
!> Adds the velocities and particle energy of a particle to the correct position of transarray and intarray
!>   only performed if sampling is enabled
!===================================================================================================================================
USE MOD_Globals       ,ONLY: abort,VECNORM
USE MOD_Particle_Vars ,ONLY: WriteMacroSurfaceValues
USE MOD_Particle_Vars ,ONLY: PartState, Species, PartSpecies
USE MOD_DSMC_Vars     ,ONLY: CollisMode, useDSMC, SpecDSMC
USE MOD_DSMC_Vars     ,ONLY: PartStateIntEn, DSMC
USE MOD_TimeDisc_Vars ,ONLY: TEnd, time
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)      :: PartID
REAL,INTENT(OUT)        :: TransArray(1:6),IntArray(1:6)
CHARACTER(*),INTENT(IN) :: SampleType
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: VeloReal, ETrans
INTEGER :: ETransID, ERotID, EVibID
!-----------------------------------------------------------------------------------------------------------------------------------
IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
  TransArray(:)=0.
  IntArray(:)=0.
  VeloReal = VECNORM(PartState(4:6,PartID))
  ETrans = 0.5 * Species(PartSpecies(PartID))%MassIC * VeloReal**2
  SELECT CASE (TRIM(SampleType))
  CASE ('old')
    ! must be old_velocity-new_velocity
    TransArray(4:6) = PartState(4:6,PartID)
    ETransID = 1
    ERotID = 1
    EVibID = 4
  CASE ('new')
    ! must be old_velocity-new_velocity
    TransArray(4:6) = -PartState(4:6,PartID)
    ETransID = 3
    ERotID = 3
    EVibID = 6
  CASE DEFAULT
    CALL abort(&
__STAMP__&
,'ERROR in AddPartInfoToSample: wrong SampleType specified. Possible types -> ( old , new )')
  END SELECT
  TransArray(ETransID) = ETrans
  IF (useDSMC .AND. CollisMode.GT.1) THEN
    IF ((SpecDSMC(PartSpecies(PartID))%InterID.EQ.2).OR.SpecDSMC(PartSpecies(PartID))%InterID.EQ.20) THEN
      IntArray(ERotID) = PartStateIntEn(2,PartID)
      IntArray(EVibID) = PartStateIntEn(1,PartID)
    END IF
  END IF
END IF
END SUBROUTINE AddPartInfoToSample


SUBROUTINE CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,IsSpeciesSwap,&
                          emission_opt,impact_opt,PartTrajectory_opt,SurfaceNormal_opt)
!===================================================================================================================================
!> Sample Wall values from Particle collisions
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars
USE MOD_Globals                ,ONLY: abort
USE MOD_DSMC_Vars              ,ONLY: SpecDSMC, useDSMC
USE MOD_DSMC_Vars              ,ONLY: CollisMode, DSMC
USE MOD_Particle_Boundary_Vars ,ONLY: SampWall, CalcSurfCollis
USE MOD_TimeDisc_Vars          ,ONLY: TEnd, time
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                 :: PartID,SurfSideID,p,q
REAL,INTENT(IN)                    :: TransArray(1:6) !1-3 trans energies (old,wall,new), 4-6 diff. trans vel. [Delta_v (x,y,z)]
REAL,INTENT(IN)                    :: IntArray(1:6) ! 1-6 internal energies (rot-old,rot-wall,rot-new,vib-old,vib-wall,vib-new)
LOGICAL,INTENT(IN)                 :: IsSpeciesSwap
LOGICAL,INTENT(IN),OPTIONAL        :: emission_opt
LOGICAL,INTENT(IN),OPTIONAL        :: impact_opt
REAL,INTENT(IN),OPTIONAL           :: PartTrajectory_opt(1:3)
REAL,INTENT(IN),OPTIONAL           :: SurfaceNormal_opt(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL        :: emission_opt_loc
LOGICAL        :: impact_opt_loc
!===================================================================================================================================
! return if sampling is not enabled
IF (.NOT.(&
         (DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)&
         )) RETURN

ASSOCIATE ( MPF  => Species(PartSpecies(PartID))%MacroParticleFactor ,&
            mass => Species(PartSpecies(PartID))%MassIC              ,&
            id   => PartSpecies(PartID)                              )
  !----  Sampling for energy (translation) accommodation at walls
  SampWall(SurfSideID)%State(SAMPWALL_ETRANSOLD,p,q)  = SampWall(SurfSideID)%State(SAMPWALL_ETRANSOLD,p,q)  + TransArray(1) * MPF
  SampWall(SurfSideID)%State(SAMPWALL_ETRANSWALL,p,q) = SampWall(SurfSideID)%State(SAMPWALL_ETRANSWALL,p,q) + TransArray(2) * MPF
  SampWall(SurfSideID)%State(SAMPWALL_ETRANSNEW,p,q)  = SampWall(SurfSideID)%State(SAMPWALL_ETRANSNEW,p,q)  + TransArray(3) * MPF

  !----  Sampling force at walls
  SampWall(SurfSideID)%State(SAMPWALL_DELTA_MOMENTUMX,p,q)= SampWall(SurfSideID)%State(SAMPWALL_DELTA_MOMENTUMX,p,q) &
                                                            + mass * TransArray(4) * MPF
  SampWall(SurfSideID)%State(SAMPWALL_DELTA_MOMENTUMY,p,q)= SampWall(SurfSideID)%State(SAMPWALL_DELTA_MOMENTUMY,p,q) &
                                                            + mass * TransArray(5) * MPF
  SampWall(SurfSideID)%State(SAMPWALL_DELTA_MOMENTUMZ,p,q)= SampWall(SurfSideID)%State(SAMPWALL_DELTA_MOMENTUMZ,p,q) &
                                                            + mass * TransArray(6) * MPF

  IF (useDSMC) THEN
    IF (CollisMode.GT.1) THEN
      IF ((SpecDSMC(PartSpecies(PartID))%InterID.EQ.2).OR.SpecDSMC(PartSpecies(PartID))%InterID.EQ.20) THEN
        !----  Sampling for internal (rotational) energy accommodation at walls
        SampWall(SurfSideID)%State(SAMPWALL_EROTOLD,p,q)  = SampWall(SurfSideID)%State(SAMPWALL_EROTOLD,p,q)  + IntArray(1) * MPF
        SampWall(SurfSideID)%State(SAMPWALL_EROTWALL,p,q) = SampWall(SurfSideID)%State(SAMPWALL_EROTWALL,p,q) + IntArray(2) * MPF
        SampWall(SurfSideID)%State(SAMPWALL_EROTNEW,p,q)  = SampWall(SurfSideID)%State(SAMPWALL_EROTNEW,p,q)  + IntArray(3) * MPF

        !----  Sampling for internal (vibrational) energy accommodation at walls
        SampWall(SurfSideID)%State(SAMPWALL_EVIBOLD,p,q)  = SampWall(SurfSideID)%State(SAMPWALL_EVIBOLD,p,q)  + IntArray(4) * MPF
        SampWall(SurfSideID)%State(SAMPWALL_EVIBWALL,p,q) = SampWall(SurfSideID)%State(SAMPWALL_EVIBWALL,p,q) + IntArray(5) * MPF
        SampWall(SurfSideID)%State(SAMPWALL_EVIBNEW,p,q)  = SampWall(SurfSideID)%State(SAMPWALL_EVIBNEW,p,q)  + IntArray(6) * MPF
      END IF
    END IF
  END IF

  ! if calcwalsample is called with emission_opt (e.g. from particle emission for evaporation) than the collision counters are not
  ! added to sampwall
  IF (PRESENT(emission_opt)) THEN
    emission_opt_loc=emission_opt
  ELSE
    emission_opt_loc=.FALSE.
  END IF
  IF (.NOT.emission_opt_loc) THEN
    !---- Counter for collisions (normal wall collisions - not to count if only SpeciesSwaps to be counted)
    IF (.NOT.CalcSurfCollis%OnlySwaps .AND. .NOT.IsSpeciesSwap) THEN
      SampWall(SurfSideID)%State(SAMPWALL_NVARS+PartSpecies(PartID),p,q) = &
          SampWall(SurfSideID)%State(SAMPWALL_NVARS+PartSpecies(PartID),p,q) + 1
    END IF
  END IF

  ! if calcwalsample is called with impact_opt=.TRUE. (e.g. particles impacting on surfaces)
  IF (PRESENT(impact_opt)) THEN
    impact_opt_loc=impact_opt
  ELSE
    impact_opt_loc=.FALSE.
  END IF

  ! Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z) and angle
  ! only works if impact_opt_loc=CalcSurfaceImpact=T
  IF(impact_opt_loc) CALL CountSurfaceImpact(SurfSideID,PartSpecies(PartID),MPF,TransArray(1),IntArray(4),IntArray(1),&
                                             PartTrajectory_opt,SurfaceNormal_opt,p,q)

END ASSOCIATE

END SUBROUTINE CalcWallSample


SUBROUTINE AnalyzeSurfaceCollisions(PartID,PartTrajectory,alpha,IsSpeciesSwap,locBCID)
!===================================================================================================================================
!> Analyzes values for particle-surface collisions
!===================================================================================================================================
! MODULES
USE MOD_Globals                ,ONLY: abort
USE MOD_Particle_Vars
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_Particle_Boundary_Vars ,ONLY: CalcSurfCollis, AnalyzeSurfCollis
USE MOD_TimeDisc_Vars          ,ONLY: TEnd, time
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                 :: PartID,locBCID
REAL,INTENT(IN)                    :: PartTrajectory(1:3), alpha
LOGICAL,INTENT(IN)                 :: IsSpeciesSwap
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! return if sampling is not enabled
IF (.NOT.(&
         (DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)&
         )) RETURN

!---- Counter for collisions (normal wall collisions - not to count if only SpeciesSwaps to be counted)
IF (.NOT.CalcSurfCollis%OnlySwaps .AND. .NOT.IsSpeciesSwap) THEN
  IF (CalcSurfCollis%AnalyzeSurfCollis .AND. (ANY(AnalyzeSurfCollis%BCs.EQ.0) .OR. ANY(AnalyzeSurfCollis%BCs.EQ.locBCID))) THEN
    AnalyzeSurfCollis%Number(PartSpecies(PartID)) = AnalyzeSurfCollis%Number(PartSpecies(PartID)) + 1
    AnalyzeSurfCollis%Number(nSpecies+1) = AnalyzeSurfCollis%Number(nSpecies+1) + 1
    IF (AnalyzeSurfCollis%Number(nSpecies+1) .GT. AnalyzeSurfCollis%maxPartNumber) THEN
      CALL abort(&
__STAMP__&
,'maxSurfCollisNumber reached!')
    END IF
    AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),1:3) = LastPartPos(1:3,PartID) + alpha * PartTrajectory(1:3)
    AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),4) = PartState(4,PartID)
    AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),5) = PartState(5,PartID)
    AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),6) = PartState(6,PartID)
    AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),7) = LastPartPos(1,PartID)
    AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),8) = LastPartPos(2,PartID)
    AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),9) = LastPartPos(3,PartID)
    AnalyzeSurfCollis%Spec(AnalyzeSurfCollis%Number(nSpecies+1))   = PartSpecies(PartID)
    AnalyzeSurfCollis%BCid(AnalyzeSurfCollis%Number(nSpecies+1))   = locBCID
  END IF
END IF
END SUBROUTINE AnalyzeSurfaceCollisions


SUBROUTINE SurfaceToPartEnergyInternal(PartID,WallTemp)
!===================================================================================================================================
!> Particle internal energies PartStateIntEn() are sampled at surface temperature
!===================================================================================================================================
USE MOD_Particle_Vars ,ONLY: PartSpecies
USE MOD_Globals_Vars  ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars     ,ONLY: CollisMode, PolyatomMolDSMC, useDSMC
USE MOD_DSMC_Vars     ,ONLY: PartStateIntEn, SpecDSMC, DSMC, VibQuantsPar
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: PartID
REAL,INTENT(IN)    :: WallTemp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: RanNum
REAL    :: NormProb
INTEGER :: VibQuant, iDOF, iPolyatMole
!-----------------------------------------------------------------------------------------------------------------------------------
IF (useDSMC .AND. CollisMode.GT.1) THEN
  IF (SpecDSMC(PartSpecies(PartID))%InterID.EQ.2.OR.SpecDSMC(PartSpecies(PartID))%InterID.EQ.20) THEN
    ! Insert new particle with internal energies sampled from surface temperature
    IF(SpecDSMC(PartSpecies(PartID))%PolyatomicMol) THEN
      ! set vibrational energy
      iPolyatMole = SpecDSMC(PartSpecies(PartID))%SpecToPolyArray
      IF(ALLOCATED(VibQuantsPar(PartID)%Quants)) DEALLOCATE(VibQuantsPar(PartID)%Quants)
      ALLOCATE(VibQuantsPar(PartID)%Quants(PolyatomMolDSMC(iPolyatMole)%VibDOF))
      PartStateIntEn( 1,PartID) = 0.0
      DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
        CALL RANDOM_NUMBER(RanNum)
        VibQuant = INT(-LOG(RanNum)*WallTemp/PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF))
        DO WHILE (VibQuant.GE.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(iDOF))
          CALL RANDOM_NUMBER(RanNum)
          VibQuant = INT(-LOG(RanNum)*WallTemp/PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF))
        END DO
        PartStateIntEn( 1,PartID) = PartStateIntEn( 1,PartID) &
                                   + (VibQuant + DSMC%GammaQuant)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
        VibQuantsPar(PartID)%Quants(iDOF)=VibQuant
      END DO
      IF (SpecDSMC(PartSpecies(PartID))%Xi_Rot.EQ.2) THEN
        CALL RANDOM_NUMBER(RanNum)
        PartStateIntEn( 2,PartID) = -BoltzmannConst*WallTemp*LOG(RanNum)
      ELSE IF (SpecDSMC(PartSpecies(PartID))%Xi_Rot.EQ.3) THEN
        CALL RANDOM_NUMBER(RanNum)
        PartStateIntEn( 2,PartID) = RanNum*10 !the distribution function has only non-negligible  values betwenn 0 and 10
        NormProb = SQRT(PartStateIntEn( 2,PartID))*EXP(-PartStateIntEn( 2,PartID))/(SQRT(0.5)*EXP(-0.5))
        CALL RANDOM_NUMBER(RanNum)
        DO WHILE (RanNum.GE.NormProb)
          CALL RANDOM_NUMBER(RanNum)
          PartStateIntEn( 2,PartID) = RanNum*10 !the distribution function has only non-negligible  values betwenn 0 and 10
          NormProb = SQRT(PartStateIntEn( 2,PartID))*EXP(-PartStateIntEn( 2,PartID))/(SQRT(0.5)*EXP(-0.5))
          CALL RANDOM_NUMBER(RanNum)
        END DO
        PartStateIntEn( 2,PartID) = PartStateIntEn( 2,PartID)*BoltzmannConst*WallTemp
      END IF
    ELSE
      ! Set vibrational energy
      CALL RANDOM_NUMBER(RanNum)
      VibQuant = INT(-LOG(RanNum)*WallTemp/SpecDSMC(PartSpecies(PartID))%CharaTVib)
      DO WHILE (VibQuant.GE.SpecDSMC(PartSpecies(PartID))%MaxVibQuant)
        CALL RANDOM_NUMBER(RanNum)
        VibQuant = INT(-LOG(RanNum)*WallTemp/SpecDSMC(PartSpecies(PartID))%CharaTVib)
      END DO
      PartStateIntEn( 1,PartID) = (VibQuant + DSMC%GammaQuant)*SpecDSMC(PartSpecies(PartID))%CharaTVib*BoltzmannConst
      ! Set rotational energy
      CALL RANDOM_NUMBER(RanNum)
      PartStateIntEn( 2,PartID) = -BoltzmannConst*WallTemp*LOG(RanNum)
    END IF
  ELSE
    ! Nullify energy for atomic species
    PartStateIntEn( 1,PartID) = 0.0
    PartStateIntEn( 2,PartID) = 0.0
  END IF
END IF
!End internal energy accomodation
END SUBROUTINE SurfaceToPartEnergyInternal



SUBROUTINE CountSurfaceImpact(SurfSideID,SpecID,MPF,ETrans,EVib,ERot,PartTrajectory,SurfaceNormal,p,q)
!===================================================================================================================================
!> Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z), angle and number of impacts
!>
!===================================================================================================================================
!USE MOD_DSMC_Vars              ,ONLY: SpecDSMC
USE MOD_Particle_Boundary_Vars ,ONLY: SampWall
USE MOD_Globals_Vars           ,ONLY: PI
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: SurfSideID          !< Surface ID
INTEGER,INTENT(IN) :: SpecID              !< Species ID
REAL,INTENT(IN)    :: MPF                 !< Particle macro particle factor
REAL,INTENT(IN)    :: ETrans              !< Translational energy of impacting particle
REAL,INTENT(IN)    :: ERot                !< Rotational energy of impacting particle
REAL,INTENT(IN)    :: EVib                !< Vibrational energy of impacting particle
REAL,INTENT(IN)    :: PartTrajectory(1:3) !< Particle trajectory vector (normalized)
REAL,INTENT(IN)    :: SurfaceNormal(1:3)  !< Surface normal vector (normalized)
INTEGER,INTENT(IN) :: p                 !< Surface sub-faces
INTEGER,INTENT(IN) :: q                 !< Surface sub-faces
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

!----- Sampling of impact energy for each species (trans, rot, vib)
SampWall(SurfSideID)%ImpactEnergy(SpecID,1,p,q)   = SampWall(SurfSideID)%ImpactEnergy(SpecID,1,p,q) + ETrans * MPF
SampWall(SurfSideID)%ImpactEnergy(SpecID,2,p,q) = SampWall(SurfSideID)%ImpactEnergy(SpecID,2,p,q)   + ERot   * MPF
SampWall(SurfSideID)%ImpactEnergy(SpecID,3,p,q) = SampWall(SurfSideID)%ImpactEnergy(SpecID,3,p,q)   + EVib   * MPF

!----- Sampling of impact vector for each species (x,y,z)
SampWall(SurfSideID)%ImpactVector(SpecID,1,p,q)   = SampWall(SurfSideID)%ImpactVector(SpecID,1,p,q) + PartTrajectory(1) * MPF
SampWall(SurfSideID)%ImpactVector(SpecID,2,p,q)   = SampWall(SurfSideID)%ImpactVector(SpecID,2,p,q) + PartTrajectory(2) * MPF
SampWall(SurfSideID)%ImpactVector(SpecID,3,p,q)   = SampWall(SurfSideID)%ImpactVector(SpecID,3,p,q) + PartTrajectory(3) * MPF

!----- Sampling of impact angle for each species
SampWall(SurfSideID)%ImpactAngle(SpecID,p,q) = SampWall(SurfSideID)%ImpactAngle(SpecID,p,q) + &
    (90.-ABS(90.-(180./PI)*ACOS(DOT_PRODUCT(PartTrajectory,SurfaceNormal)))) * MPF

!----- Sampling of impact number for each species
SampWall(SurfSideID)%ImpactNumber(SpecID,p,q) = SampWall(SurfSideID)%ImpactNumber(SpecID,p,q) + MPF

END SUBROUTINE CountSurfaceImpact


SUBROUTINE StoreBoundaryParticleProperties(iPart,PartPos,PartTrajectory,SurfaceNormal)
!----------------------------------------------------------------------------------------------------------------------------------!
! Save particle position, velocity and species to PartDataBoundary container for writing to .h5 later
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals                ,ONLY: abort
USE MOD_Particle_Vars          ,ONLY: usevMPF,PartMPF,PartSpecies,Species,PartState
USE MOD_Particle_Boundary_Vars ,ONLY: PartStateBoundary,PartStateBoundaryVecLength
USE MOD_TimeDisc_Vars          ,ONLY: time
USE MOD_Globals_Vars           ,ONLY: PI
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)   :: iPart
REAL,INTENT(IN)      :: PartPos(1:3)
REAL,INTENT(IN)      :: PartTrajectory(1:3)
REAL,INTENT(IN)      :: SurfaceNormal(1:3)
INTEGER              :: dims(2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                 :: MPF
! Temporary arrays
REAL, ALLOCATABLE    :: PartStateBoundary_tmp(:,:) ! (1:10,1:NParts) 1st index: x,y,z,vx,vy,vz,SpecID,Ekin,MPF,time,impact angle
!                                                  !                 2nd index: 1 to number of boundary-crossed particles
INTEGER              :: ALLOCSTAT
!===================================================================================================================================
IF (usevMPF) THEN
  MPF = PartMPF(iPart)
ELSE
  MPF = Species(PartSpecies(iPart))%MacroParticleFactor
END IF

dims = SHAPE(PartStateBoundary)

ASSOCIATE( iMax => PartStateBoundaryVecLength )
  ! Increase maximum number of boundary-impact particles
  iMax = iMax + 1

  ! Check if array maximum is reached. 
  ! If this happens, re-allocate the arrays and increase their size (every time this barrier is reached, double the size)
  IF(iMax.GT.dims(2))THEN

    ! --- PartStateBoundary ---
    ALLOCATE(PartStateBoundary_tmp(1:10,1:dims(2)), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(&
          __STAMP__&
          ,'ERROR in particle_boundary_tools.f90: Cannot allocate PartStateBoundary_tmp temporary array!')
    ! Save old data
    PartStateBoundary_tmp(1:10,1:dims(2)) = PartStateBoundary(1:10,1:dims(2))

    ! Re-allocate PartStateBoundary to twice the size
    DEALLOCATE(PartStateBoundary)
    ALLOCATE(PartStateBoundary(1:10,1:2*dims(2)), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(&
          __STAMP__&
          ,'ERROR in particle_boundary_tools.f90: Cannot allocate PartStateBoundary array!')
    PartStateBoundary(1:10,        1:  dims(2)) = PartStateBoundary_tmp(1:10,1:dims(2))
    PartStateBoundary(1:10,dims(2)+1:2*dims(2)) = 0.

  END IF

  PartStateBoundary(1:3,iMax) = PartPos
  PartStateBoundary(4:6,iMax) = PartState(4:6,iPart)
  PartStateBoundary(7  ,iMax) = REAL(PartSpecies(iPart))
  PartStateBoundary(8  ,iMax) = MPF
  PartStateBoundary(9  ,iMax) = time
  PartStateBoundary(10 ,iMax) = (90.-ABS(90.-(180./PI)*ACOS(DOT_PRODUCT(PartTrajectory,SurfaceNormal))))
END ASSOCIATE

END SUBROUTINE StoreBoundaryParticleProperties


SUBROUTINE SortArray(EndID,ArrayA,ArrayB)
!----------------------------------------------------------------------------------------------------------------------------------!
! sort ArrayA (e.g. SideIDToSurfID) in ascending order of ArrayB (e.g. GlobalUniqueSideID)
! based on running MINLOC for each element .NE. -1
!
! In the following example sorting from iSide = nBCSides+1 to nSides is performed
!
! BEFORE:
!  
!      iSide         SideIDToSurfID(iSide)      GlobalUniqueSideID(iSide)
!        10                   -1                           3
!        11               --- 10 <--                      10
!        12               |   -1   |                      16
!        13               |   -1   |                      11
!        14               |   -1   |                      15
!        15               --> 11 ---                       4
!        16                   -1                           6
!    
! AFTER:
!    
!      iSide         SideIDToSurfID(iSide)      GlobalUniqueSideID(iSide)
!        10                   -1                           3
!        11                   11                          10
!        12                   -1                          16
!        13                   -1                          11
!        14                   -1                          15
!        15                   10                           4
!        16                   -1                           6
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)    :: EndID
INTEGER,INTENT(INOUT) :: ArrayA(EndID)
INTEGER,INTENT(IN)    :: ArrayB(EndID)
! insert IO variables here
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,idx
LOGICAL :: Mask(EndID) ! only the elements for which Mask is .TRUE. are considered
LOGICAL :: Mask_tmp(EndID)
INTEGER :: ArrayA_temp(EndID)
!===================================================================================================================================
ArrayA_temp = ArrayA
Mask = .TRUE.
DO i = 1, EndID
  IF(ArrayA(i).EQ.-1) THEN
    Mask(i) = .FALSE.
  END IF
END DO
Mask_tmp = Mask
DO i = 1, EndID
  IF(.NOT.Mask_tmp(i)) CYCLE
   idx         = MINLOC(ArrayB,1,Mask)
   ArrayA(idx) = ArrayA_temp(i)
   Mask(idx)   = .FALSE.
END DO

END SUBROUTINE SortArray


RECURSIVE SUBROUTINE SortArray2(StartID,EndID,ArrayA,ArrayB)
!----------------------------------------------------------------------------------------------------------------------------------!
! sort ArrayA (e.g. SideIDToSurfID) in ascending order of ArrayB (e.g. GlobalUniqueSideID)
! based on QuickSort
!
! In the following example sorting from iSide = nBCSides+1 to nSides is performed
!
! BEFORE:
!  
!      iSide         SideIDToSurfID(iSide)      GlobalUniqueSideID(iSide)
!        10                   -1                           3
!        11               --- 10 <--                      10
!        12               |   -1   |                      16
!        13               |   -1   |                      11
!        14               |   -1   |                      15
!        15               --> 11 ---                       4
!        16                   -1                           6
!    
! AFTER:
!    
!      iSide         SideIDToSurfID(iSide)      GlobalUniqueSideID(iSide)
!        10                   -1                           3
!        11                   11                          10
!        12                   -1                          16
!        13                   -1                          11
!        14                   -1                          15
!        15                   10                           4
!        16                   -1                           6
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)    :: StartID,EndID
INTEGER,INTENT(INOUT) :: ArrayA(*)
INTEGER,INTENT(INOUT) :: ArrayB(*)
! insert IO variables here
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,idx
!LOGICAL :: Mask(EndID) ! only the elements for which Mask is .TRUE. are considered
!LOGICAL :: Mask_tmp(EndID)
!INTEGER :: ArrayA_temp(EndID)
INTEGER :: Center,temp
!===================================================================================================================================

Center = ArrayA((StartID+EndID)/2)
i      = StartID
j      = EndID

DO
  DO WHILE (ArrayA(i).LT.Center)
    i=i+1
  END DO
  DO WHILE (Center.LT.ArrayA(j))
    j=j-1
  END DO
  IF(i.GE.j) EXIT
  ! Sort actual array: ArrayA
  temp      = ArrayA(i)
  ArrayA(i) = ArrayA(j)
  ArrayA(j) = temp
  ! Sort corresponding array: ArrayB
  temp      = ArrayB(i)
  ArrayB(i) = ArrayB(j)
  ArrayB(j) = temp

  i=i+1
  j=j-1
END DO

!IF(StartID.LT.i-1) CALL SortArray2(StartID,EndID,ArrayA,ArrayB)

IF(StartID.LT.i-1  ) CALL SortArray2(StartID , i-1   , ArrayA , ArrayB)
IF(j+1    .LT.EndID) CALL SortArray2(j+1     , EndID , ArrayA , ArrayB)


!   RETURN
!   ArrayA_temp = ArrayA
!   Mask = .TRUE.
!   DO i = 1, EndID
!     IF(ArrayA(i).EQ.-1) THEN
!       Mask(i) = .FALSE.
!     END IF
!   END DO
!   Mask_tmp = Mask
!   DO i = 1, EndID
!     IF(.NOT.Mask_tmp(i)) CYCLE
!      idx         = MINLOC(ArrayB,1,Mask)
!      ArrayA(idx) = ArrayA_temp(i)
!      Mask(idx)   = .FALSE.
!   END DO

END SUBROUTINE SortArray2


SUBROUTINE DielectricSurfaceCharge(iPart,ElemID,PartTrajectory,alpha)
!----------------------------------------------------------------------------------------------------------------------------------!
! description
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
USE MOD_Globals         ,ONLY: abort,myrank
USE MOD_Mesh_Vars       ,ONLY: nElems
USE MOD_part_operations ,ONLY: CreateParticle
USE MOD_part_tools      ,ONLY: isChargedParticle
USE MOD_Particle_Vars   ,ONLY: PDM,PartSpecies,LastPartPos
USE MOD_PICDepo_Tools   ,ONLY: DepositParticleOnNodes
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)    :: iPart
INTEGER,INTENT(IN)    :: ElemID
REAL,INTENT(IN)       :: PartTrajectory(1:3), alpha
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: NewPartID
!===================================================================================================================================
! Sanity checks
IF(.NOT.PDM%ParticleInside(iPart))THEN
  IPWRITE (*,*) "iPart  :", iPart
  IPWRITE (*,*) "ElemID :", ElemID
  CALL abort(&
      __STAMP__&
      ,'Dielectric particle-surface interaction: Particle not inside element.')
ELSEIF(PartSpecies(iPart).LT.0)THEN
  IF(myrank.eq.1)THEN
    IPWRITE (*,*) "iPart =", iPart
    IPWRITE (*,*) "PDM%ParticleVecLength =", PDM%ParticleVecLength
  END IF ! myrank.eq.1
  CALL abort(&
      __STAMP__&
      ,'Negative speciesID')
END IF ! PartSpecies(iPart)

IF(isChargedParticle(iPart))THEN
  IF(ElemID.GT.nElems)THEN
    ! Particle is now located in halo element: Create phantom particle, which is sent to new host Processor and removed there (set
    ! negative SpeciesID in order to remove particle in host Processor)
    CALL CreateParticle(-PartSpecies(iPart),LastPartPos(1:3,iPart)+PartTrajectory(1:3)*alpha,ElemID,(/0.,0.,0./),0.,0.,0.,NewPartID)
    ! Set inside to F (it is set to T in SendNbOfParticles if species ID is negative)
    PDM%ParticleInside(NewPartID)=.FALSE.
  ELSE ! Deposit single particle charge on surface here and 
    CALL DepositParticleOnNodes(iPart,LastPartPos(1:3,iPart)+PartTrajectory(1:3)*alpha,ElemID)
  END IF ! ElemID.GT.nElems
END IF ! isChargedParticle(iPart)

END SUBROUTINE DielectricSurfaceCharge


REAL FUNCTION GetWallTemperature(PartID,locBCID,PartTrajectory,alpha)
!===================================================================================================================================
!> Determine the wall temperature, current options: determine a temperature based on an imposed gradient or use a fixed temperature
!===================================================================================================================================
USE MOD_Globals                 ,ONLY: DOTPRODUCT, VECNORM
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound
USE MOD_Particle_Vars           ,ONLY: LastPartPos
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: locBCID, PartID
REAL, INTENT(IN)                :: PartTrajectory(3), alpha
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: TempGradLength, POI(3), POI_projected(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
IF(PartBound%WallTemp2(locBCID).GT.0.0) THEN
  POI(1:3) = LastPartPos(1:3,PartID) + PartTrajectory(1:3)*alpha
  POI_projected(1:3) = PartBound%TempGradStart(1:3,locBCID) &
                      + DOT_PRODUCT((POI(1:3) - PartBound%TempGradStart(1:3,locBCID)),PartBound%TempGradVec(1:3,locBCID)) &
                        / DOTPRODUCT(PartBound%TempGradVec(1:3,locBCID)) * PartBound%TempGradVec(1:3,locBCID)
  TempGradLength = VECNORM(POI_projected(1:3))/VECNORM(PartBound%TempGradVec(1:3,locBCID))
  IF(TempGradLength.LT.0.0) THEN
    GetWallTemperature = PartBound%WallTemp(locBCID)
  ELSE IF(TempGradLength.GT.1.0) THEN
    GetWallTemperature = PartBound%WallTemp2(locBCID)
  ELSE
    GetWallTemperature = PartBound%WallTemp(locBCID) + TempGradLength * PartBound%WallTempDelta(locBCID)
  END IF
ELSE
  GetWallTemperature = PartBound%WallTemp(locBCID)
END IF

END FUNCTION GetWallTemperature

END MODULE MOD_Particle_Boundary_Tools
