!==================================================================================================================================
! Copyright (c) 2010 - 2019 Wladimir Reschke
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
INTERFACE PartEnergyToSurface
  MODULE PROCEDURE PartEnergyToSurface
END INTERFACE

INTERFACE SurfaceToPartEnergy
  MODULE PROCEDURE SurfaceToPartEnergy
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

PUBLIC :: PartEnergyToSurface
PUBLIC :: SurfaceToPartEnergy
PUBLIC :: LIQUIDEVAP
PUBLIC :: LIQUIDREFL
PUBLIC :: ALPHALIQUID
PUBLIC :: BETALIQUID
PUBLIC :: TSURUTACONDENSCOEFF
!===================================================================================================================================

CONTAINS


SUBROUTINE PartEnergyToSurface(PartID,SpecID,Transarray,IntArray)
!===================================================================================================================================
!> particle energy is zeroed and sampling arrays are written
!===================================================================================================================================
USE MOD_Particle_Vars ,ONLY: WriteMacroSurfaceValues
USE MOD_Particle_Vars ,ONLY: PartState,Species
USE MOD_Globals_Vars  ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars     ,ONLY: CollisMode, PolyatomMolDSMC
USE MOD_DSMC_Vars     ,ONLY: PartStateIntEn, SpecDSMC, DSMC, VibQuantsPar
USE MOD_TimeDisc_Vars ,ONLY: TEnd, time
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: PartID
INTEGER,INTENT(IN) :: SpecID
REAL,INTENT(OUT)   :: TransArray(1:6),IntArray(1:6)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: oldVelo(1:3)
REAL    :: VeloReal
REAL    :: EtraOld, EtraWall, EtraNew
REAL    :: ErotWall, ErotNew
REAL    :: EvibWall, EvibNew
REAL    :: RanNum
INTEGER :: VibQuant, iDOF, iPolyatMole
!-----------------------------------------------------------------------------------------------------------------------------------
IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
  oldVelo(1:3) = PartState(PartID,4:6)
  PartState(PartID,4:6)  = 0.

  VeloReal = SQRT(DOT_PRODUCT(oldVelo,oldVelo))
  EtraOld = 0.5 * Species(SpecID)%MassIC * VeloReal**2
  EtraWall = 0.0
  EtraNew = EtraWall

  TransArray(1) = EtraOld
  TransArray(2) = EtraWall
  TransArray(3) = EtraNew
  ! must be old_velocity-new_velocity
  TransArray(4:6) = oldVelo(1:3)-PartState(PartID,4:6)

  !---- Internal energy accommodation
  IF (CollisMode.GT.1) THEN
    IF (SpecDSMC(SpecID)%InterID.EQ.2) THEN
      !---- Rotational energy accommodation
      CALL RANDOM_NUMBER(RanNum)
      ErotWall = 0
      ErotNew  = 0
      IntArray(1) = PartStateIntEn(PartID,2)
      IntArray(2) = ErotWall
      IntArray(3) = ErotNew
      PartStateIntEn(PartID,2) = ErotNew
      !---- Vibrational energy accommodation
      EvibWall = 0.0
      EvibNew  = 0.0
      IF(SpecDSMC(SpecID)%PolyatomicMol) THEN
        iPolyatMole = SpecDSMC(SpecID)%SpecToPolyArray
        DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
          IntArray(4) = IntArray(4) + (VibQuantsPar(PartID)%Quants(iDOF) + DSMC%GammaQuant) * BoltzmannConst &
                      * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) * Species(SpecID)%MacroParticleFactor
        END DO
      ELSE
        VibQuant     = NINT(PartStateIntEn(PartID,1)/(BoltzmannConst*SpecDSMC(SpecID)%CharaTVib) &
                    - DSMC%GammaQuant)
        IntArray(4) = (VibQuant + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(SpecID)%CharaTVib
      END IF
      IntArray(5) = EvibWall
      IntArray(6) = EvibNew
      PartStateIntEn(PartID,1) = EvibNew
    END IF
  END IF
  !End internal energy accomodation
END IF
END SUBROUTINE PartEnergyToSurface


SUBROUTINE SurfaceToPartEnergy(PartID,SpecID,WallTemp,Transarray,IntArray)
!===================================================================================================================================
!> Particle internal energies are set sampled from surface temperature maxwell distributed and sampling arrays are written
!===================================================================================================================================
USE MOD_Particle_Vars ,ONLY: WriteMacroSurfaceValues
USE MOD_Particle_Vars ,ONLY: PartState,Species
USE MOD_Globals_Vars  ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars     ,ONLY: CollisMode, PolyatomMolDSMC
USE MOD_DSMC_Vars     ,ONLY: PartStateIntEn, SpecDSMC, DSMC, VibQuantsPar
USE MOD_TimeDisc_Vars ,ONLY: TEnd, time
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: PartID
INTEGER,INTENT(IN) :: SpecID
REAL,INTENT(IN)    :: WallTemp
REAL,INTENT(OUT)   :: TransArray(1:6),IntArray(1:6)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL    :: oldVelo(1:3)
REAL    :: VeloReal
!REAL    :: EtraOld, EtraWall, EtraNew
REAL    :: ErotWall, ErotNew
REAL    :: EvibWall, EvibNew
REAL    :: RanNum
REAL    :: NormProb
INTEGER :: VibQuant, iDOF, iPolyatMole
!-----------------------------------------------------------------------------------------------------------------------------------
IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
  ! sample
  TransArray(1) = 0. ! EtraOld
  TransArray(2) = 0. ! EtraWall
  VeloReal = SQRT(DOT_PRODUCT(PartState(PartID,4:6),PartState(PartID,4:6))) !oldVelo,oldVelo))
  TransArray(3) = 0.5 * Species(SpecID)%MassIC * VeloReal**2 ! EtraNew
  ! must be old_velocity-new_velocity
  TransArray(4:6) = -PartState(PartID,4:6)
END IF

! set internal energies of new particle
IF (CollisMode.GT.1) THEN
  IF (SpecDSMC(SpecID)%InterID.EQ.2) THEN
    ErotWall = 0.
    EvibWall = 0.
    ! Insert new particle with internal energies sampled from surface temperature
    IF(SpecDSMC(SpecID)%PolyatomicMol) THEN
      ! set vibrational energy
      iPolyatMole = SpecDSMC(SpecID)%SpecToPolyArray
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
      IF (SpecDSMC(SpecID)%Xi_Rot.EQ.2) THEN
        CALL RANDOM_NUMBER(RanNum)
        PartStateIntEn(PartID, 2) = -BoltzmannConst*WallTemp*LOG(RanNum)
      ELSE IF (SpecDSMC(SpecID)%Xi_Rot.EQ.3) THEN
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
      VibQuant = INT(-LOG(RanNum)*WallTemp/SpecDSMC(SpecID)%CharaTVib)
      DO WHILE (VibQuant.GE.SpecDSMC(SpecID)%MaxVibQuant)
        CALL RANDOM_NUMBER(RanNum)
        VibQuant = INT(-LOG(RanNum)*WallTemp/SpecDSMC(SpecID)%CharaTVib)
      END DO
      !evtl muss partstateinten nochmal geändert werden, mpi, resize etc..
      PartStateIntEn(PartID, 1) = (VibQuant + DSMC%GammaQuant)*SpecDSMC(SpecID)%CharaTVib*BoltzmannConst
      ! Set rotational energy
      CALL RANDOM_NUMBER(RanNum)
      PartStateIntEn(PartID, 2) = -BoltzmannConst*WallTemp*LOG(RanNum)
    END IF
    ErotNew = PartStateIntEn(PartID,2)
    EvibNew = PartStateIntEn(PartID,1)
    !---- Rotational energy
    IntArray(1) = 0.0
    IntArray(2) = ErotWall
    IntArray(3) = ErotNew
    !---- Vibrational energy
    IntArray(4) = 0.0
    IntArray(5) = EvibWall
    IntArray(6) = EvibNew
  ELSE
    ! Nullify energy for atomic species
    PartStateIntEn(PartID, 1) = 0.0
    PartStateIntEn(PartID, 2) = 0.0
    IntArray(1:6)=0.0
  END IF
END IF
!End internal energy accomodation
END SUBROUTINE SurfaceToPartEnergy


REAL FUNCTION LIQUIDEVAP(beta,x,sigma)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: beta,x,sigma
!===================================================================================================================================
liquidEvap=(1-beta*exp(-0.5*(x/sigma)**2))/(1-beta/2)  *   x/sigma**2  *  exp(-0.5*(x/sigma)**2)
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
!===================================================================================================================================
liquidRefl= (1-alpha+alpha*beta*exp(-0.5*(x/sigma)**2))/(1-alpha*(1-beta/2))  *   x/sigma**2  *  exp(-0.5*(x/sigma)**2)
END FUNCTION

FUNCTION ALPHALIQUID(specID,temp) RESULT(alpha)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
!USE MOD_SurfModel_Vars, ONLY:
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
alpha = 0.9
END FUNCTION

FUNCTION BETALIQUID(specID,temp) RESULT(beta)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
!USE MOD_SurfModel_Vars, ONLY:
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
beta = 0.280
END FUNCTION

FUNCTION TSURUTACONDENSCOEFF(SpecID,normalVelo,temp) RESULT(sigma)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
USE MOD_Particle_Vars          ,ONLY: Species
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
END FUNCTION


END MODULE MOD_Particle_Boundary_Tools
