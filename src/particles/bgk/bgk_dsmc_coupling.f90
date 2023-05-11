!==================================================================================================================================
! Copyright (c) 2018 - 2019 Marcel Pfeiffer
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

MODULE MOD_BGK_DSMC_Coupling
!===================================================================================================================================
!> Distinction between DSMC and BGK based on a predefined continuum-breakdown-criterion
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
PUBLIC :: CBC_DoDSMC
!===================================================================================================================================

CONTAINS

LOGICAL FUNCTION CBC_DoDSMC(iElem)
!===================================================================================================================================
!> Test of different continuum-breakdown criteria
!> TO-DO: create separate functions for every breakdown criterium
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_BGK_Vars
USE MOD_BGK_CollOperator       ,ONLY: CalcViscosityThermalCondColIntVHS
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
USE MOD_DSMC_Analyze           ,ONLY: CalcMeanFreePath
USE MOD_Particle_Vars          ,ONLY: PEM, nSpecies, PartSpecies, Species, usevMPF, PartState
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemVolume_Shared, ElemToElemMapping, ElemToElemInfo, ElemMidPoint_Shared
USE MOD_Mesh_Vars              ,ONLY: offsetElem, nElems
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID, GetGlobalElemID
USE MOD_part_tools             ,ONLY: GetParticleWeight
USE MOD_DSMC_Vars              ,ONLY: DSMC, SpecDSMC, CollisMode, PartStateIntEn, RadialWeighting, VarWeighting, PolyatomMolDSMC
USE MOD_Particle_Analyze_Tools ,ONLY: CalcTelec,CalcTVibPoly
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iPart, iLoop, iSpec, iVal, jVal, iCoord, iPolyatMole
INTEGER                      :: CNElemID, nNbElems, LocNbElem, GlobNbElem, CnNbElem
REAL                         :: SpecPartNum(nSpecies), MolPartNum(nSpecies), MFP, totalWeight, NbVolume, NbDistance, NbDens
REAL                         :: RefVelo(nSpecies,3), RefVelo2(nSpecies,3), NbVelo(nSpecies,3), NbVelo2(nSpecies,3)
REAL                         :: TempRot, TempVib, TempTotal, RefTemp(3), NbTemp(3), DoF_Rot(nSpecies), DoF_Vib(nSpecies)
REAL                         :: TempTrans(3), SpeciesTemp(nSpecies), DoF_RotTotal, DoF_VibTotal, RefTempTotal, NbTempTotal 
REAL                         :: ERot(nSpecies), EVib(nSpecies), TVib_TempFac, RefDens, RefVeloTotal, NbVeloTotal 
REAL                         :: DensGradient, TempGradient, VeloGradient, Knudsen_Dens, Knudsen_Temp, Knudsen_Velo, Knudsen_NonEq
REAL                         :: NbtotalWeight, NbSpecPartNum(nSpecies), StaticPressure
REAL                         :: BulkVelo(1:3), RelativeVelo(1:3), RelVeloTotal, HeatVector(3), StressTensor(3,3), ChapmanEnskog
!===================================================================================================================================
! Set all needed variables to zero for the element in the loop
DensGradient = 0.0; TempGradient = 0.0; VeloGradient = 0.0; NbVolume = 0.0; RefVelo = 0.0; RefVelo2 = 0.0; SpecPartNum = 0.0
RefVeloTotal = 0.0; RefTempTotal = 0.0; TempTrans = 0.0;  SpeciesTemp = 0.0; ERot = 0.0; EVib  = 0.0
MolPartNum = 0.0; TempVib = 0.0; TempRot = 0.0; DoF_Rot = 0.0; DoF_RotTotal = 0.0; DoF_Vib = 0.0; DoF_VibTotal = 0.0; RelVeloTotal = 0.0
StressTensor = 0.0; HeatVector = 0.0

! Reference element
CNElemID = GetCNElemID(iElem+offSetElem)

! Determine the number of different species and the total particle number in the element by a loop over all particles 
iPart = PEM%pStart(iElem)
DO iLoop = 1, PEM%pNumber(iElem)
  SpecPartNum(PartSpecies(iPart)) = SpecPartNum(PartSpecies(iPart)) + GetParticleWeight(iPart)

  ! Particle velocities (squared) weighted by the particle number per element
  RefVelo(PartSpecies(iPart),1:3) = RefVelo(PartSpecies(iPart),1:3) + PartState(4:6,iPart) * GetParticleWeight(iPart)
  RefVelo2(PartSpecies(iPart),1:3) = RefVelo2(PartSpecies(iPart),1:3) + PartState(4:6,iPart)**2 * GetParticleWeight(iPart)

  ! Rotational and vibrational energy of the reference element
  IF ((SpecDSMC(PartSpecies(iPart))%InterID.EQ.2).OR.(SpecDSMC(PartSpecies(iPart))%InterID.EQ.20)) THEN
    MolPartNum(PartSpecies(iPart)) = MolPartNum(PartSpecies(iPart)) + GetParticleWeight(iPart)
    EVib(PartSpecies(iPart)) = EVib(PartSpecies(iPart)) + (PartStateIntEn(1,iPart) - SpecDSMC(PartSpecies(iPart))%EZeroPoint) * GetParticleWeight(iPart)
    ERot(PartSpecies(iPart)) = ERot(PartSpecies(iPart)) + PartStateIntEn(2,iPart) * GetParticleWeight(iPart)
  END IF
  iPart = PEM%pNext(iPart)
END DO

! Total particle number in the element
totalWeight = SUM(SpecPartNum)

! Number density in the reference element = Particle number/volume
IF(usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.VarWeighting%DoVariableWeighting) THEN
  RefDens = totalWeight / ElemVolume_Shared(CNElemID)
ELSE
  RefDens = totalWeight * Species(1)%MacroParticleFactor / ElemVolume_Shared(CNElemID)
END IF

! Calculation of the total mean velocity under consideration of all directions
DO iLoop = 1, 3
  RefVeloTotal = RefVeloTotal + (SUM(RefVelo(:,iLoop))/totalWeight)**2
  ! Bulk velocity
  BulkVelo(iLoop) = SUM(RefVelo(:,iLoop))/totalWeight
END DO
RefVeloTotal = SQRT(RefVeloTotal)

! Calculation of the total mean translational temperature
DO iSpec = 1, nSpecies
  IF (SpecPartNum(iSpec).GT.1.0) THEN 
    RefTemp  = Species(iSpec)%MassIC/BoltzmannConst * ((RefVelo2(iSpec,:)/SpecPartNum(iSpec)) - (RefVelo(iSpec,:)/SpecPartNum(iSpec))**2)
    RefTempTotal = RefTempTotal + ((RefTemp(1) + RefTemp(2) + RefTemp(3)) / 3.)*SpecPartNum(iSpec)
    ! Translational temperature in x/y/z direction
    TempTrans = TempTrans + RefTemp*SpecPartNum(iSpec)
    ! Translational temperature for each species
    SpeciesTemp(iSpec) = (RefTemp(1) + RefTemp(2) + RefTemp(3)) / 3.

    ! Internal energies
    IF ((CollisMode.EQ.2).OR.(CollisMode.EQ.3))THEN
      IF ((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
        ! Vibrational temperature
        IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
          IF( (EVib(iSpec)/MolPartNum(iSpec)).GT.0.0 ) THEN
            TempVib = TempVib + &
                        CalcTVibPoly(EVib(iSpec)/MolPartNum(iSpec)+SpecDSMC(iSpec)%EZeroPoint, iSpec) * MolPartNum(iSpec)
          END IF
          iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
          DoF_Vib(iSpec) = PolyatomMolDSMC(iPolyatMole)%VibDOF 
        ELSE
          TVib_TempFac = EVib(iSpec) / (MolPartNum(iSpec) * BoltzmannConst * SpecDSMC(iSpec)%CharaTVib)
          IF ((EVib(iSpec) /MolPartNum(iSpec)).GT.0.0) THEN
            TempVib = TempVib + SpecDSMC(iSpec)%CharaTVib / LOG(1. + 1./(TVib_TempFac)) * MolPartNum(iSpec)
          END IF
          DoF_Vib(iSpec) = 1. 
        END IF !PolyatomicMol

        ! Rotational temperature
        TempRot = TempRot + 2. * ERot(iSpec) / (MolPartNum(iSpec)*BoltzmannConst*REAL(SpecDSMC(iSpec)%Xi_Rot)) * MolPartNum(iSpec)
        DoF_Rot(iSpec) = REAL(SpecDSMC(iSpec)%Xi_Rot) 
        DoF_RotTotal = DoF_RotTotal + DoF_Rot(iSpec) * MolPartNum(iSpec)
        DoF_VibTotal = DoF_VibTotal + DoF_Vib(iSpec) * MolPartNum(iSpec)
      END IF ! InterID
    END IF ! CollisMode
  END IF
END DO  
IF (RefTempTotal.GT.0.) THEN
  RefTempTotal = RefTempTotal/totalweight
ELSE 
  RefTempTotal = 0.
END IF
TempTrans = TempTrans/totalweight

! Vibrational and rotational temperature weighted by the particle numbers
IF (SUM(MolPartNum(:)).GT.0) THEN
  TempVib = TempVib / SUM(MolPartNum(:))
  TempRot = TempRot / SUM(MolPartNum(:))
  DoF_VibTotal = DoF_VibTotal / SUM(MolPartNum(:))
  DoF_RotTotal = DoF_RotTotal / SUM(MolPartNum(:))
END IF

! Total temperature and comparison to the mean translational temperature
TempTotal = (RefTempTotal*3. + DoF_VibTotal*TempVib + DoF_RotTotal*TempRot) / (3. + DoF_VibTotal + DoF_RotTotal)

iPart = PEM%pStart(iElem)
DO iLoop = 1, PEM%pNumber(iElem)
  RelativeVelo(1:3) = PartState(4:6,iPart) - BulkVelo(1:3)
  RelVeloTotal = RelativeVelo(1)**2 + RelativeVelo(2)**2 + RelativeVelo(3)**2

  ! Shear stress tensor
  DO iVal = 1,3
    DO jVal = 1,3
      StressTensor(iVal,jVal) = StressTensor(iVal,jVal) + RelativeVelo(iVal)*RelativeVelo(jVal) * GetParticleWeight(iPart)
    END DO
  END DO

  ! Static pressure
  StaticPressure = (StressTensor(1,1) + StressTensor(2,2) + StressTensor(3,3))/3.
  StressTensor(1,1) = StressTensor(1,1) - StaticPressure
  StressTensor(2,2) = StressTensor(2,2) - StaticPressure
  StressTensor(3,3) = StressTensor(3,3) - StaticPressure
 
  ! Heat vector
  HeatVector(1:3) = HeatVector(1:3) + RelativeVelo(1:3) * RelVeloTotal * GetParticleWeight(iPart)

  iPart = PEM%pNext(iPart)
END DO

IF(PEM%pNumber(iElem).GT.2) THEN
  HeatVector = HeatVector*PEM%pNumber(iElem)*PEM%pNumber(iElem)/((PEM%pNumber(iElem)-1.)*(PEM%pNumber(iElem)-2.)*totalWeight)
ELSE
  HeatVector = 0.0
END IF

StressTensor = -(StressTensor/totalWeight)

! 1) Global Knudsen number: 2.55 mm for the Rothe nozzle expansion 
MFP = CalcMeanFreePath(SpecPartNum,SUM(SpecPartNum),ElemVolume_Shared(CNElemID))

! 2) Local Knudsen number: relative gradient between the cell and all neighboring cells for different flow properties
! Considered properties: particle density, velocity, translational temperature
nNbElems = ElemToElemMapping(2,CNElemID)
! Loop over all neighbouring elements
DO iVal = 1, nNbElems 
  CnNbElem = ElemToElemInfo(ElemToElemMapping(1,CNElemID)+iVal)
  GlobNbElem = GetGlobalElemID(CnNbElem)
  LocNbElem = GlobNbElem-offSetElem

  IF ((LocNBElem.LT.1).OR.(LocNBElem.GT.nElems)) CYCLE

  ! Distance between the two midpoints of the cell
  NbDistance = 0.0
  DO iCoord=1, 3
    NbDistance = NbDistance + (ElemMidPoint_Shared(iCoord,CNElemID)-ElemMidPoint_Shared(iCoord,CnNbElem))**2
  END DO
  NbDistance = SQRT(NbDistance)

  ! Add element volume to the total volume of all neighbour elements
  NbVolume = NbVolume + ElemVolume_Shared(CnNbElem)

  NbSpecPartNum = 0.0; NbVelo = 0.0; NbVelo2 = 0.0; NbVeloTotal = 0.0; NbTempTotal = 0.0 

  iPart = PEM%pStart(LocNbElem)
  ! Loop over all particles in the neighbour cells and add them up for the total particle number per species
  DO iLoop = 1, PEM%pNumber(LocNbElem)
    NbSpecPartNum(PartSpecies(iPart)) = NbSpecPartNum(PartSpecies(iPart)) + GetParticleWeight(iPart)

    ! Particle velocity and bulk velocity per species
    NbVelo(PartSpecies(iPart),1:3) = NbVelo(PartSpecies(iPart),1:3) + PartState(4:6,iPart) * GetParticleWeight(iPart)
    NbVelo2(PartSpecies(iPart),1:3) = NbVelo2(PartSpecies(iPart),1:3) + PartState(4:6,iPart)**2 * GetParticleWeight(iPart)

    iPart = PEM%pNext(iPart)
  END DO

  NbtotalWeight = SUM(NbSpecPartNum)

  ! Total number density in the neighbour element
  IF(usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.VarWeighting%DoVariableWeighting) THEN
    NbDens = NbtotalWeight / ElemVolume_Shared(CnNbElem)
  ELSE
    NbDens = NbtotalWeight * Species(1)%MacroParticleFactor / ElemVolume_Shared(CnNbElem)
  END IF

  ! Sum of the density gradient of all neighbour elements, weighted by the volume of the neighbour element
  DensGradient = DensGradient + ABS(RefDens-NbDens/NbDistance)*ElemVolume_Shared(CnNbElem)

  ! Velocity in the neighbour element
  IF (NbDens.GT.0.) THEN
    DO iLoop = 1, 3
      NbVeloTotal = NbVeloTotal + (SUM(NbVelo(:,iLoop))/NbtotalWeight)**2
    END DO
    NbVeloTotal = SQRT(NbVeloTotal)
  END IF
  
  ! Sum of the velocity gradient of all neighbour elements, weighted by the volume of the neighbour element
  VeloGradient = VeloGradient + ABS(RefVeloTotal-NbVeloTotal/NbDistance)*ElemVolume_Shared(CnNbElem)

  ! Total translational temperature in the neighbouring element
  IF (NbDens.GT.0.) THEN
    DO iSpec = 1, nSpecies
      IF (NbSpecPartNum(iSpec).GE.1.) THEN
        NbTemp  = Species(iSpec)%MassIC/BoltzmannConst * ((NbVelo2(iSpec,:)/NbSpecPartNum(iSpec)) - (NbVelo(iSpec,:)/NbSpecPartNum(iSpec))**2)
        NbTempTotal = NbTempTotal + ((NbTemp(1) + NbTemp(2) + NbTemp(3)) / 3.)*NbSpecPartNum(iSpec)
      END IF
    END DO  

    NbTempTotal = NbTempTotal/totalweight
  END IF

  ! Sum of the temperature gradient of all neighbour elements, weighted by the volume of the neighbour element
  TempGradient = TempGradient+ ABS(RefTempTotal-NbTempTotal/NbDistance)*ElemVolume_Shared(CnNbElem)
END DO ! iNbElem

! 2a) Local Knudsen number of the density
DensGradient = DensGradient/NbVolume
Knudsen_Dens = MFP*DensGradient/RefDens

! 2b) Local Knudsen number of the velocity
VeloGradient = VeloGradient/NbVolume
Knudsen_Velo = MFP*VeloGradient/RefVeloTotal

! 2c) Local Knudsen number of the temperature
TempGradient = TempGradient/NbVolume
Knudsen_Temp = MFP*TempGradient/RefTempTotal

! 3) Continuum-breakdown based on a local thermal non-equilibrium
Knudsen_NonEq = SQRT(((TempTrans(1)-TempTotal)**2 + (TempTrans(2)-TempTotal)**2 + (TempTrans(3)-TempTotal)**2 &
              + DoF_RotTotal*(TempRot-TempTotal)**2 + DoF_VibTotal*(TempVib-TempTotal)**2 )/((3. + DoF_VibTotal + DoF_RotTotal)*TempTotal**2))

! 4) Simplififed Chapman-Enskog parameter
ChapmanEnskog = 0.1

! Write out the non-equilibrium parameters
CBC%OutputKnudsen(1,iElem) = MFP/CBC%CharLength
CBC%OutputKnudsen(2,iElem) = Knudsen_Dens
CBC%OutputKnudsen(3,iElem) = Knudsen_Velo
CBC%OutputKnudsen(4,iElem) = Knudsen_Temp
CBC%OutputKnudsen(5,iElem) = Knudsen_NonEq
CBC%OutputKnudsen(6,iElem) = MAXVAL(StressTensor)
CBC%OutputKnudsen(7,iElem) = MAXVAL(HeatVector)

! Definition of continuum-breakdown in the cell and swtich between BGK and DSMC
SELECT CASE (TRIM(CBC%SwitchCriterium))
CASE('GlobalKnudsen')
   IF ((MFP/CBC%CharLength).GE.CBC%MaxGlobalKnudsen) THEN
    CBC_DoDSMC = .TRUE.
  ELSE
    CBC_DoDSMC = .FALSE.
  END IF
CASE('LocalKnudsen')
  IF (Knudsen_Dens.GT.CBC%MaxLocalKnudsen) THEN
    CBC_DoDSMC = .TRUE.
  ELSE IF (Knudsen_Velo.GT.CBC%MaxLocalKnudsen) THEN
    CBC_DoDSMC = .TRUE.
  ELSE IF (Knudsen_Temp.GT.CBC%MaxLocalKnudsen) THEN
    CBC_DoDSMC = .TRUE. 
  ELSE
    CBC_DoDSMC = .FALSE.
  END IF
CASE('ThermNonEq')
  IF (Knudsen_NonEq.GT.CBC%MaxThermNonEq) THEN
    CBC_DoDSMC = .TRUE. 
  ELSE
    CBC_DoDSMC = .FALSE.
  END IF
CASE('Combination')
  IF (Knudsen_Dens.GT.CBC%MaxLocalKnudsen) THEN
    CBC_DoDSMC = .TRUE.
  ELSE IF (Knudsen_Velo.GT.CBC%MaxLocalKnudsen) THEN
    CBC_DoDSMC = .TRUE.
  ELSE IF (Knudsen_Temp.GT.CBC%MaxLocalKnudsen) THEN
    CBC_DoDSMC = .TRUE. 
  ELSE IF (Knudsen_NonEq.GT.CBC%MaxThermNonEq) THEN
    CBC_DoDSMC = .TRUE. 
  ELSE
    CBC_DoDSMC = .FALSE.
  END IF
CASE('ChapmanEnskog')
  IF (ChapmanEnskog.GT.0.2) THEN
    CBC_DoDSMC = .TRUE.
  ELSE
    CBC_DoDSMC = .FALSE.
  END IF
END SELECT

RETURN

END FUNCTION CBC_DoDSMC

END MODULE MOD_BGK_DSMC_Coupling
