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
INTEGER                      :: iPart, iLoop, iSpec, iVal, iCoord, iPolyatMole
INTEGER                      :: CNElemID, nNbElems, LocNbElem, GlobNbElem, CnNbElem
REAL                         :: SpecPartNum(nSpecies), MolPartNum(nSpecies), RefTemp(3), NbTemp(3)
REAL                         :: RefVelo(nSpecies,3), RefVelo2(nSpecies,3), NbVelo(nSpecies,3), NbVelo2(nSpecies,3)
REAL                         :: TempRot, TempVib, TempTotal, DoF_Rot, DoF_Vib
REAL                         :: TempTrans(3)
REAL                         :: ERot(nSpecies), EVib(nSpecies), Knudsen_NonEq, TVib_TempFac
REAL                         :: MFP, totalWeight, NbVolume, NbDistance
REAL                         :: RefDens, RefVeloTotal, RefTempTotal, NbDens, NbVeloTotal, NbTempTotal 
REAL                         :: DensGradient, TempGradient, VeloGradient, Knudsen_Dens, Knudsen_Temp, Knudsen_Velo
!===================================================================================================================================
! Set all needed variables to zero for the element in the loop
DensGradient = 0.
TempGradient = 0.
VeloGradient = 0.
NbVolume     = 0.
RefVelo      = 0.
RefVelo2     = 0.

! Reference element
CNElemID = GetCNElemID(iElem+offSetElem)

! Determine the number of different species and the total particle number in the element by a loop over all particles 
iPart = PEM%pStart(iElem)
SpecPartNum  = 0.
DO iLoop = 1, PEM%pNumber(iElem)
  SpecPartNum(PartSpecies(iPart)) = SpecPartNum(PartSpecies(iPart)) + GetParticleWeight(iPart)

  ! Particle velocities weighted by the particle number per element
  RefVelo(PartSpecies(iPart),1:3) = RefVelo(PartSpecies(iPart),1:3) + PartState(4:6,iPart) * GetParticleWeight(iPart)
  RefVelo2(PartSpecies(iPart),1:3) = RefVelo2(PartSpecies(iPart),1:3) + PartState(4:6,iPart)**2 * GetParticleWeight(iPart)
  iPart = PEM%pNext(iPart)
END DO

totalWeight = SUM(SpecPartNum)

! Number density in the reference element = Particle number/volume
IF(usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.VarWeighting%DoVariableWeighting) THEN
  RefDens = totalWeight / ElemVolume_Shared(CNElemID)
ELSE
  RefDens = totalWeight * Species(1)%MacroParticleFactor / ElemVolume_Shared(CNElemID)
END IF

! Calculation of the total mean velocity under consideration of all directions
RefVeloTotal = 0.
DO iLoop = 1, 3
  RefVeloTotal = RefVeloTotal + (SUM(RefVelo(:,iLoop))/totalWeight)**2
END DO
RefVeloTotal = SQRT(RefVeloTotal)

! Calculation of the total mean translational temperature
RefTempTotal = 0.
TempTrans = 0.
DO iSpec = 1, nSpecies
  RefTemp  = Species(iSpec)%MassIC/BoltzmannConst * ((RefVelo2(iSpec,:)/SpecPartNum(iSpec)) - (RefVelo(iSpec,:)/SpecPartNum(iSpec))**2)
  TempTrans = TempTrans + RefTemp*SpecPartNum(iSpec)
  RefTempTotal = RefTempTotal + ((RefTemp(1) + RefTemp(2) + RefTemp(3)) / 3.)*SpecPartNum(iSpec)
END DO  
RefTempTotal = RefTempTotal/totalweight
TempTrans = TempTrans/totalweight

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
  NbDistance = 0.
  DO iCoord=1, 3
    NbDistance = NbDistance + (ElemMidPoint_Shared(iCoord,CNElemID)-ElemMidPoint_Shared(iCoord,CnNbElem))**2
  END DO
  NbDistance = SQRT(NbDistance)

  ! Add element volume to the total volume of all neighbour elements
  NbVolume = NbVolume + ElemVolume_Shared(CnNbElem)

  SpecPartNum = 0.
  NbVelo      = 0.
  NbVelo2     = 0.
  iPart = PEM%pStart(LocNbElem)
  ! Loop over all particles in the neighbour cells and add them up for the total particle number per species
  DO iLoop = 1, PEM%pNumber(LocNbElem)
    SpecPartNum(PartSpecies(iPart)) = SpecPartNum(PartSpecies(iPart)) + GetParticleWeight(iPart)

    ! Particle velocity and bulk velocity per species
    NbVelo(PartSpecies(iPart),1:3) = NbVelo(PartSpecies(iPart),1:3) + PartState(4:6,iPart) * GetParticleWeight(iPart)
    NbVelo2(PartSpecies(iPart),1:3) = NbVelo2(PartSpecies(iPart),1:3) + PartState(4:6,iPart)**2 * GetParticleWeight(iPart)

    iPart = PEM%pNext(iPart)
  END DO

  totalWeight = SUM(SpecPartNum)

  ! Total number density in the neighbour element
  IF(usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.VarWeighting%DoVariableWeighting) THEN
    NbDens = totalWeight / ElemVolume_Shared(CnNbElem)
  ELSE
    NbDens = totalWeight * Species(1)%MacroParticleFactor / ElemVolume_Shared(CnNbElem)
  END IF

  ! Sum of the density gradient of all neighbour elements, weighted by the volume of the neighbour element
  DensGradient = DensGradient + ABS(RefDens-NbDens/NbDistance)*ElemVolume_Shared(CnNbElem)

  ! Velocity in the neighbour element
  NbVeloTotal = 0.
  IF (NbDens.EQ.0.) THEN
    NbVeloTotal = 0.
  ELSE 
    DO iLoop = 1, 3
      NbVeloTotal = NbVeloTotal + (SUM(NbVelo(:,iLoop))/totalWeight)**2
    END DO
    NbVeloTotal = SQRT(NbVeloTotal)
  END IF
  
  ! Sum of the velocity gradient of all neighbour elements, weighted by the volume of the neighbour element
  VeloGradient = VeloGradient + ABS(RefVeloTotal-NbVeloTotal/NbDistance)*ElemVolume_Shared(CnNbElem)

  ! Total translational temperature in the neighbouring element
  NbTempTotal = 0.
  IF (NbDens.EQ.0.) THEN
    NbTempTotal = 0.
  ELSE
    DO iSpec = 1, nSpecies
      NbTemp  = Species(iSpec)%MassIC/BoltzmannConst * ((NbVelo2(iSpec,:)/SpecPartNum(iSpec)) - (NbVelo(iSpec,:)/SpecPartNum(iSpec))**2)
      NbTempTotal = NbTempTotal + ((NbTemp(1) + NbTemp(2) + NbTemp(3)) / 3.)*SpecPartNum(iSpec)
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
ERot  = 0.
EVib  = 0.

! Loop over all particles in the cell
iPart = PEM%pStart(iElem)
SpecPartNum = 0.
MolPartNum = 0.
DO iLoop = 1, PEM%pNumber(iElem)
  SpecPartNum(PartSpecies(iPart)) = SpecPartNum(PartSpecies(iPart)) + GetParticleWeight(iPart)

  ! Rotational and vibrational energy of the reference element
  IF ((SpecDSMC(PartSpecies(iPart))%InterID.EQ.2).OR.(SpecDSMC(PartSpecies(iPart))%InterID.EQ.20)) THEN
    MolPartNum(PartSpecies(iPart)) = MolPartNum(PartSpecies(iPart)) + GetParticleWeight(iPart)
    EVib(PartSpecies(iPart)) = EVib(PartSpecies(iPart)) + (PartStateIntEn(1,iPart) - SpecDSMC(PartSpecies(iPart))%EZeroPoint) * GetParticleWeight(iPart)
    ERot(PartSpecies(iPart)) = ERot(PartSpecies(iPart)) + PartStateIntEn(2,iPart) * GetParticleWeight(iPart)
  END IF
  iPart = PEM%pNext(iPart)
END DO ! iLoop

totalWeight = SUM(SpecPartNum)

TempVib = 0.
TempRot = 0.
DoF_Rot = 0. 
DoF_Vib = 0.

! Loop over the species to determine the total rotational and vibrational energy and temperature
DO iSpec=1, nSpecies
  IF ((CollisMode.EQ.2).OR.(CollisMode.EQ.3))THEN
    IF ((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
      ! Vibrational temperature
      IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
        IF( (EVib(iSpec)/MolPartNum(iSpec)) .GT. 0.0 ) THEN
          TempVib = TempVib + &
                      CalcTVibPoly(EVib(iSpec)/MolPartNum(iSpec)+SpecDSMC(iSpec)%EZeroPoint, iSpec) * MolPartNum(iSpec)
        END IF
        iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
        DoF_Vib = DoF_Vib + PolyatomMolDSMC(iPolyatMole)%VibDOF * MolPartNum(iSpec)
      ELSE
        TVib_TempFac = EVib(iSpec) / (MolPartNum(iSpec) * BoltzmannConst * SpecDSMC(iSpec)%CharaTVib)
        IF ((EVib(iSpec) /MolPartNum(iSpec)).GT.0.0) THEN
          TempVib = TempVib + SpecDSMC(iSpec)%CharaTVib / LOG(1. + 1./(TVib_TempFac)) * MolPartNum(iSpec)
        END IF
        DoF_Vib = 1. * MolPartNum(iSpec)
      END IF !PolyatomicMol
      ! Rotational temperature
      TempRot = TempRot + 2. * ERot(iSpec) / (MolPartNum(iSpec)*BoltzmannConst*REAL(SpecDSMC(iSpec)%Xi_Rot)) * MolPartNum(iSpec)
      DoF_Rot = DoF_Rot +  REAL(SpecDSMC(iSpec)%Xi_Rot) * MolPartNum(iSpec)
    END IF ! InterID
  END IF ! CollisMode
END DO ! iSpec

! Vibrational and rotational temperature weighted by the particle numbers
IF (SUM(MolPartNum(:)).GT.0) THEN
  TempVib = TempVib / SUM(MolPartNum(:))
  TempRot = TempRot / SUM(MolPartNum(:))
  DoF_Vib = DoF_Vib / SUM(MolPartNum(:))
  DoF_Rot = DoF_Rot / SUM(MolPartNum(:))
END IF

! Total temperature and comparison to the mean translational temperature
TempTotal = (RefTempTotal*3. + DoF_Vib*TempVib + DoF_Rot*TempRot) / (3. + DoF_Vib + DoF_Rot)

Knudsen_NonEq = SQRT(((TempTrans(1)-TempTotal)**2 + (TempTrans(2)-TempTotal)**2 + (TempTrans(3)-TempTotal)**2 &
              + DoF_Rot*(TempRot-TempTotal)**2 + DoF_Vib*(TempVib-TempTotal)**2 )/((3. + DoF_Vib + DoF_Rot)*TempTotal**2))

! Write out the non-equilibrium parameters
BGK_OutputKnudsen(1,iElem) = MFP/0.00255
BGK_OutputKnudsen(2,iElem) = Knudsen_Dens
BGK_OutputKnudsen(3,iElem) = Knudsen_Velo
BGK_OutputKnudsen(4,iElem) = Knudsen_Temp
BGK_OutputKnudsen(5,iElem) = Knudsen_NonEq

! Definition of continuum-breakdown in the cell and swtich between BGK and DSMC
SELECT CASE (TRIM(BGKDSMC_SwitchCriterium))
CASE('GlobalKnudsen')
   IF ((MFP/0.00255).GE.0.1) THEN
    CBC_DoDSMC = .TRUE.
  ELSE
    CBC_DoDSMC = .FALSE.
  END IF
CASE('LocalKnudsen')
  IF (Knudsen_Dens.GT.0.1) THEN
    CBC_DoDSMC = .TRUE.
  ELSE IF (Knudsen_Velo.GT.0.1) THEN
    CBC_DoDSMC = .TRUE.
  ELSE IF (Knudsen_Temp.GT.0.1) THEN
    CBC_DoDSMC = .TRUE. 
  ELSE
    CBC_DoDSMC = .FALSE.
  END IF
CASE('ThermNonEq')
  IF (Knudsen_NonEq.GT.0.05) THEN
    CBC_DoDSMC = .TRUE. 
  ELSE
    CBC_DoDSMC = .FALSE.
  END IF
CASE('Combination')
  IF (Knudsen_Dens.GT.0.1) THEN
    CBC_DoDSMC = .TRUE.
  ELSE IF (Knudsen_Velo.GT.0.1) THEN
    CBC_DoDSMC = .TRUE.
  ELSE IF (Knudsen_Temp.GT.0.1) THEN
    CBC_DoDSMC = .TRUE. 
  ELSE IF (Knudsen_NonEq.GT.0.05) THEN
    CBC_DoDSMC = .TRUE. 
  ELSE
    CBC_DoDSMC = .FALSE.
  END IF
END SELECT

RETURN

END FUNCTION CBC_DoDSMC

END MODULE MOD_BGK_DSMC_Coupling
