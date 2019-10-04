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

INTERFACE LIQUIDEVAP
  MODULE PROCEDURE LIQUIDEVAP
END INTERFACE

INTERFACE LIQUIDREFL
  MODULE PROCEDURE LIQUIDREFL
END INTERFACE

INTERFACE ALPHALIQUID
  MODULE PROCEDURE ALPHALIQUID
END INTERFACE

INTERFACE BETALIQUID
  MODULE PROCEDURE BETALIQUID
END INTERFACE

INTERFACE TSURUTACONDENSCOEFF
  MODULE PROCEDURE TSURUTACONDENSCOEFF
END INTERFACE

INTERFACE CountSurfaceImpact
  MODULE PROCEDURE CountSurfaceImpact
END INTERFACE

INTERFACE BoundaryParticleOutput
  MODULE PROCEDURE BoundaryParticleOutput
END INTERFACE

PUBLIC :: AddPartInfoToSample
PUBLIC :: CalcWallSample
PUBLIC :: AnalyzeSurfaceCollisions
PUBLIC :: SurfaceToPartEnergyInternal
PUBLIC :: LIQUIDEVAP
PUBLIC :: LIQUIDREFL
PUBLIC :: ALPHALIQUID
PUBLIC :: BETALIQUID
PUBLIC :: TSURUTACONDENSCOEFF
PUBLIC :: CountSurfaceImpact
PUBLIC :: BoundaryParticleOutput
PUBLIC :: SortArray
!===================================================================================================================================

CONTAINS


SUBROUTINE AddPartInfoToSample(PartID,Transarray,IntArray,SampleType)
!===================================================================================================================================
!> Adds the velocities and particle energy of a particle to the correct position of transarray and intarray
!>   only performed if sampling is enabled
!===================================================================================================================================
USE MOD_Globals       ,ONLY: abort
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
  VeloReal = SQRT(DOT_PRODUCT(PartState(PartID,4:6),PartState(PartID,4:6)))
  ETrans = 0.5 * Species(PartSpecies(PartID))%MassIC * VeloReal**2
  SELECT CASE (TRIM(SampleType))
  CASE ('old')
    ! must be old_velocity-new_velocity
    TransArray(4:6) = PartState(PartID,4:6)
    ETransID = 1
    ERotID = 1
    EVibID = 4
  CASE ('new')
    ! must be old_velocity-new_velocity
    TransArray(4:6) = -PartState(PartID,4:6)
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
      IntArray(ERotID) = PartStateIntEn(PartID,2)
      IntArray(EVibID) = PartStateIntEn(PartID,1)
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
    AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),1:3) = LastPartPos(PartID,1:3) + alpha * PartTrajectory(1:3)
    AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),4) = PartState(PartID,4)
    AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),5) = PartState(PartID,5)
    AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),6) = PartState(PartID,6)
    AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),7) = LastPartPos(PartID,1)
    AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),8) = LastPartPos(PartID,2)
    AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),9) = LastPartPos(PartID,3)
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
      PartStateIntEn(PartID, 1) = 0.0
      DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
        CALL RANDOM_NUMBER(RanNum)
        VibQuant = INT(-LOG(RanNum)*WallTemp/PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF))
        DO WHILE (VibQuant.GE.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(iDOF))
          CALL RANDOM_NUMBER(RanNum)
          VibQuant = INT(-LOG(RanNum)*WallTemp/PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF))
        END DO
        PartStateIntEn(PartID, 1) = PartStateIntEn(PartID, 1) &
                                   + (VibQuant + DSMC%GammaQuant)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
        VibQuantsPar(PartID)%Quants(iDOF)=VibQuant
      END DO
      IF (SpecDSMC(PartSpecies(PartID))%Xi_Rot.EQ.2) THEN
        CALL RANDOM_NUMBER(RanNum)
        PartStateIntEn(PartID, 2) = -BoltzmannConst*WallTemp*LOG(RanNum)
      ELSE IF (SpecDSMC(PartSpecies(PartID))%Xi_Rot.EQ.3) THEN
        CALL RANDOM_NUMBER(RanNum)
        PartStateIntEn(PartID, 2) = RanNum*10 !the distribution function has only non-negligible  values betwenn 0 and 10
        NormProb = SQRT(PartStateIntEn(PartID, 2))*EXP(-PartStateIntEn(PartID, 2))/(SQRT(0.5)*EXP(-0.5))
        CALL RANDOM_NUMBER(RanNum)
        DO WHILE (RanNum.GE.NormProb)
          CALL RANDOM_NUMBER(RanNum)
          PartStateIntEn(PartID, 2) = RanNum*10 !the distribution function has only non-negligible  values betwenn 0 and 10
          NormProb = SQRT(PartStateIntEn(PartID, 2))*EXP(-PartStateIntEn(PartID, 2))/(SQRT(0.5)*EXP(-0.5))
          CALL RANDOM_NUMBER(RanNum)
        END DO
        PartStateIntEn(PartID, 2) = PartStateIntEn(PartID, 2)*BoltzmannConst*WallTemp
      END IF
    ELSE
      ! Set vibrational energy
      CALL RANDOM_NUMBER(RanNum)
      VibQuant = INT(-LOG(RanNum)*WallTemp/SpecDSMC(PartSpecies(PartID))%CharaTVib)
      DO WHILE (VibQuant.GE.SpecDSMC(PartSpecies(PartID))%MaxVibQuant)
        CALL RANDOM_NUMBER(RanNum)
        VibQuant = INT(-LOG(RanNum)*WallTemp/SpecDSMC(PartSpecies(PartID))%CharaTVib)
      END DO
      PartStateIntEn(PartID, 1) = (VibQuant + DSMC%GammaQuant)*SpecDSMC(PartSpecies(PartID))%CharaTVib*BoltzmannConst
      ! Set rotational energy
      CALL RANDOM_NUMBER(RanNum)
      PartStateIntEn(PartID, 2) = -BoltzmannConst*WallTemp*LOG(RanNum)
    END IF
  ELSE
    ! Nullify energy for atomic species
    PartStateIntEn(PartID, 1) = 0.0
    PartStateIntEn(PartID, 2) = 0.0
  END IF
END IF
!End internal energy accomodation
END SUBROUTINE SurfaceToPartEnergyInternal


PURE REAL FUNCTION LIQUIDEVAP(beta,x,sigma)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: beta,x,sigma
REAL            :: betaLoc
!===================================================================================================================================
betaLoc = beta
IF (betaLoc.GE.2.) betaLoc = 2. - 1e-10
IF (betaLoc.LT.0.) betaLoc = 0.

liquidEvap=(1-betaLoc*exp(-0.5*(x/sigma)**2))/(1-betaLoc/2)  *   x/sigma**2  *  exp(-0.5*(x/sigma)**2)
IF (liquidEvap.LT.0.) liquidEvap = 0.
END FUNCTION

REAL FUNCTION LIQUIDREFL(alpha,beta,x,sigma)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: alpha,beta,x,sigma
REAL            :: betaLoc, alphaLoc
!===================================================================================================================================
betaLoc = beta
IF (betaLoc.GE.2.) betaLoc = 2. - 1e-10
IF (betaLoc.LT.0.) betaLoc = 0.
alphaLoc = alpha
IF (alphaLoc.GT.1.) alphaLoc = 1.
IF (alphaLoc.LT.0.) alphaLoc = 0.

if (alphaLoc.GE.1.) then
  if (betaLoc.LE.0) then
    liquidRefl = x/sigma**2  *  exp(-0.5*(x/sigma)**2)
  else
    liquidRefl = (betaLoc*exp(-0.5*(x/sigma)**2))/(1.-(1.-betaLoc/2.))  *   x/sigma**2  *  exp(-0.5*(x/sigma)**2)
  end if
else
  liquidRefl = (1.-alphaLoc+alphaLoc*betaLoc*exp(-0.5*(x/sigma)**2))/(1.-alphaLoc*(1.-betaLoc/2.)) &
             * x/sigma**2 * exp(-0.5*(x/sigma)**2)
end if

IF (liquidRefl.LT.0.) liquidRefl = 0.
END FUNCTION

PURE FUNCTION ALPHALIQUID(specID,temp) RESULT(alpha)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_SurfaceModel_Vars ,ONLY: SpecSurf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: specID
REAL,INTENT(IN)    :: temp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: alpha
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SELECT CASE (SpecSurf(specID)%condensCase)
CASE (1)
  alpha = SpecSurf(specID)%liquidAlpha
CASE (2)
  alpha = exp(-((4-BETALIQUID(specID,temp))/(2*(2-BETALIQUID(specID,temp)))-1))
END SELECT
IF (alpha.GT.1.) alpha = 1.
IF (alpha.LT.0.) alpha = 0.
END FUNCTION

PURE FUNCTION BETALIQUID(specID,temp) RESULT(beta)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_SurfaceModel_Vars ,ONLY: SpecSurf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: specID
REAL,INTENT(IN)    :: temp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: beta
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SELECT CASE (SpecSurf(specID)%condensCase)
CASE (1)
  beta = SpecSurf(specID)%liquidBeta
CASE (2)
  beta = SpecSurf(specID)%liquidBetaCoeff(1)*temp**5 &
       + SpecSurf(specID)%liquidBetaCoeff(2)*temp**4 &
       + SpecSurf(specID)%liquidBetaCoeff(3)*temp**3 &
       + SpecSurf(specID)%liquidBetaCoeff(4)*temp**2 &
       + SpecSurf(specID)%liquidBetaCoeff(5)*temp    &
       + SpecSurf(specID)%liquidBetaCoeff(6)
END SELECT
IF (beta.GE.2.) beta = 2. - 1e-10
IF (beta.LT.0.) beta=0.
END FUNCTION

FUNCTION TSURUTACONDENSCOEFF(SpecID,normalVelo,temp) RESULT(sigma)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars  ,ONLY: BoltzmannConst
USE MOD_Particle_Vars ,ONLY: Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: specID
REAL,INTENT(IN)    :: normalVelo,temp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: sigma
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
sigma = ALPHALIQUID(specID,temp)*(1-BETALIQUID(specID,temp)*exp(-normalVelo**2*Species(specID)%MassIC/(2*Boltzmannconst*temp)))
IF (sigma.LT.0.) sigma = 0.
IF (sigma.GT.1.) sigma = 1.
END FUNCTION


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
    ABS(0.5*PI - ACOS(DOT_PRODUCT(PartTrajectory,SurfaceNormal))) * (180. / PI) * MPF

!----- Sampling of impact number for each species
SampWall(SurfSideID)%ImpactNumber(SpecID,p,q) = SampWall(SurfSideID)%ImpactNumber(SpecID,p,q) + MPF

END SUBROUTINE CountSurfaceImpact


SUBROUTINE BoundaryParticleOutput(iPart,PartPos) 
!----------------------------------------------------------------------------------------------------------------------------------!
! Save particle position, velocity and species to PartDataBoundary container for writing to .h5 later
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals                ,ONLY: abort
USE MOD_Particle_Vars          ,ONLY: usevMPF,PartMPF,PartSpecies,Species,PartState,PDM
USE MOD_Particle_Boundary_Vars ,ONLY: PartStateBoundary,PartStateBoundaryVecLength,PartStateBoundarySpec
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)  :: iPart
REAL,INTENT(IN)     :: PartPos(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: MPF
!===================================================================================================================================
IF (usevMPF) THEN
  MPF = PartMPF(iPart)
ELSE
  MPF = Species(PartSpecies(iPart))%MacroParticleFactor
END IF

ASSOCIATE( iMax => PartStateBoundaryVecLength )
  iMax = iMax + 1
  IF(iMax.GT.PDM%MaxParticleNumber)THEN
    CALL abort(&
        __STAMP__&
        ,'BoundaryParticleOutput: PartStateBoundaryVecLength.GT.PDM%MaxParticleNumber. iMax=', IntInfoOpt=iMax)
  END IF
  PartStateBoundary(1:3,iMax) = PartPos
  PartStateBoundary(4:6,iMax) = PartState(iPart,4:6)
  PartStateBoundarySpec(iMax) = PartSpecies(iPart)
END ASSOCIATE

END SUBROUTINE BoundaryParticleOutput

SUBROUTINE SortArray(EndID,ArrayA,ArrayB)
!----------------------------------------------------------------------------------------------------------------------------------!
! sort arryA in ascending order of arrayB
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
LOGICAL :: unsorted(EndID)
LOGICAL :: unsorted_tmp(EndID)
INTEGER :: ArrayA_temp(EndID)
!===================================================================================================================================
ArrayA_temp=ArrayA
unsorted = .TRUE.
DO i=1, EndID
  IF(ArrayA(i).EQ.-1) THEN
    unsorted(i)=.FALSE.
  END IF
END DO
unsorted_tmp=unsorted
DO i = 1, EndID
  IF(.NOT.unsorted_tmp(i)) CYCLE
   idx=MINLOC(ArrayB,1,unsorted)
   ArrayA(i) = ArrayA_temp(idx)
   unsorted(idx) = .FALSE.
END DO

END SUBROUTINE SortArray


END MODULE MOD_Particle_Boundary_Tools
