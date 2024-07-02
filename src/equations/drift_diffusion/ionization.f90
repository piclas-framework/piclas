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


MODULE MOD_Ionization
! MODULES
IMPLICIT NONE
PRIVATE

!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC::InsertNewIons
!==================================================================================================================================
CONTAINS

SUBROUTINE InsertNewIons()
!==================================================================================================================================
!> Insert new ions due to background gas ionization
!==================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc ! PP_N
USE MOD_FV_Vars           ,ONLY: U_FV
USE MOD_Equation_Vars     ,ONLY: E
USE MOD_TimeDisc_Vars     ,ONLY: dt
USE MOD_part_operations   ,ONLY: CreateParticle
USE MOD_Particle_Tracking ,ONLY: ParticleInsideCheck
USE MOD_DSMC_Vars         ,ONLY: BGGas
USE MOD_Particle_Vars     ,ONLY: nSpecies, Species
USE MOD_Transport_Data    ,ONLY: CalcDriftDiffusionCoeffAr,CalcDriftDiffusionCoeffH2,CalcIonizationRate
USE MOD_Particle_Mesh_Vars,ONLY: BoundsOfElem_Shared, ElemVolume_Shared
USE MOD_part_tools        ,ONLY: CalcVelocity_maxwell_particle
USE MOD_Mesh_Vars         ,ONLY: offsetElem
USE MOD_Mesh_Tools        ,ONLY: GetCNElemID
USE MOD_part_emission_tools,ONLY: CalcVelocity_maxwell_lpn
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: nPart, iPart, iSpecBG, iSpecIon, iSpec, ElemID, GlobalElemID
REAL                                 :: iRan , RandomPos(1:3), DeltaPartDens, realPartNumber, Velocity(1:3)
REAL                                 :: RotEnergy, VibEnergy, ElecEnergy, ionRate, mu ,D
LOGICAL                              :: InsideFlag
!===================================================================================================================================
DO iSpec = 1, nSpecies
  IF(BGGas%BackgroundSpecies(iSpec)) THEN
    iSpecBG = iSpec
  ELSE
    iSpecIon = iSpec
  END IF
END DO

DO ElemID=1,PP_nElems
  CALL CalcIonizationRate(VECNORM(E(1:3,0,0,0,ElemID)),BGGas%NumberDensity(iSpecBG),ionRate)
  CALL CalcDriftDiffusionCoeffAr(VECNORM(E(1:3,0,0,0,ElemID)),BGGas%NumberDensity(iSpecBG),mu,D) ! do this outside and save coeff

  DeltaPartDens = ionRate*mu*VECNORM(E(1:3,0,0,0,ElemID))*U_FV(1,0,0,0,ElemID)*dt !TODO: interpolate E
  GlobalElemID = ElemID+offsetElem

  ASSOCIATE( Bounds => BoundsOfElem_Shared(1:2,1:3,GlobalElemID) ) ! 1-2: Min, Max value; 1-3: x,y,z
    CALL RANDOM_NUMBER(iRan)
    realPartNumber = DeltaPartDens * ElemVolume_Shared(GetCNElemID(GlobalElemID))
    nPart = INT(realPartNumber / Species(iSpecIon)%MacroParticleFactor + iRan)
    DO iPart = 1, nPart
      InsideFlag=.FALSE.
      DO WHILE(.NOT.InsideFlag)
        CALL RANDOM_NUMBER(RandomPos)
        RandomPos(1:3) = Bounds(1,1:3) + RandomPos(1:3)*(Bounds(2,1:3)-Bounds(1,1:3))
        InsideFlag = ParticleInsideCheck(RandomPos,iPart,GlobalElemID)
      END DO

      IF(BGGas%UseDistribution) THEN
        Velocity = CalcVelocity_maxwell_particle(iSpecIon,BGGas%Distribution(iSpecBG,4:6,ElemID)) &
                                  + BGGas%Distribution(iSpecBG,1:3,ElemID)
      ELSE
        CALL CalcVelocity_maxwell_lpn(FractNbr=iSpecIon, Vec3D=Velocity, iInit=1)
      END IF

      RotEnergy = 0.
      VibEnergy = 0.
      ElecEnergy = 0.

      CALL CreateParticle(iSpecIon,RandomPos,GlobalElemID,Velocity,RotEnergy,VibEnergy,ElecEnergy)

    END DO ! nPart
  END ASSOCIATE
END DO

END SUBROUTINE InsertNewIons

END MODULE MOD_Ionization