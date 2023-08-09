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
PUBLIC :: CBC_DoDSMC, CalcCellProp, CalcDensVeloTemp
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
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iPart, CNElemID, nNbElems, LocNbElem, GlobNbElem, CnNbElem, iVal, iCoord, iLoop
REAL                         :: SpecPartNum(nSpecies), MFP, totalWeight, NbVolume, NbDistance
REAL                         :: RefVelo, NbVelo, RefDens, NbDens, RefTemp, NbTemp, NbtotalWeight, NbSpecPartNum(nSpecies)
REAL                         :: DensGradient, TempGradient, VeloGradient, Knudsen_Dens, Knudsen_Temp, Knudsen_Velo, Knudsen_NonEq
REAL                         :: HeatVector(3), StressTensor(3,3), ChapmanEnskog
!===================================================================================================================================
CBC_DoDSMC = .FALSE.

! Reference element
CNElemID = GetCNElemID(iElem+offSetElem)

! Definition of continuum-breakdown in the cell and swtich between BGK and DSMC
SELECT CASE (TRIM(CBC%SwitchCriterium))
CASE('GlobalKnudsen')
  SpecPartNum = 0.0

  iPart = PEM%pStart(iElem)
  DO iLoop = 1, PEM%pNumber(iElem)
    SpecPartNum(PartSpecies(iPart)) = SpecPartNum(PartSpecies(iPart)) + GetParticleWeight(iPart)
    iPart = PEM%pNext(iPart)
  END DO

  MFP = CalcMeanFreePath(SpecPartNum,SUM(SpecPartNum),ElemVolume_Shared(CNElemID))

  ! Write out the non-equilibrium parameters
  CBC%OutputKnudsen(1,iElem) = MFP/CBC%CharLength

   IF ((MFP/CBC%CharLength).GE.CBC%MaxGlobalKnudsen) THEN
    CBC_DoDSMC = .TRUE.
  END IF

CASE('LocalKnudsen')
  ! Set all needed variables to zero for the element in the loop
  DensGradient = 0.0; TempGradient = 0.0; VeloGradient = 0.0; NbVolume = 0.0

  ! Calculate the density, velocity and temperature of the reference element
  CALL CalcDensVeloTemp(iElem, RefDens, RefVelo, RefTemp, SpecPartNum)

  MFP = CalcMeanFreePath(SpecPartNum,SUM(SpecPartNum),ElemVolume_Shared(CNElemID))

  ! Total particle number in the element
  totalWeight = SUM(SpecPartNum)

  ! Relative gradient between the cell and all neighboring cells for the particle density, velocity, translational temperature
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

    ! Calculate the density, velocity and temperature of the neighbour element
    CALL CalcDensVeloTemp(iElem, NbDens, NbVelo, NbTemp, NbSpecPartNum)
    NbtotalWeight = SUM(NbSpecPartNum)

    ! Sum of the density gradient of all neighbour elements, weighted by the volume of the neighbour element
    DensGradient = DensGradient + ABS((RefDens-NbDens)/NbDistance)*ElemVolume_Shared(CnNbElem)
    ! Sum of the velocity gradient of all neighbour elements, weighted by the volume of the neighbour element
    VeloGradient = VeloGradient + ABS((RefVelo-NbVelo)/NbDistance)*ElemVolume_Shared(CnNbElem)
    ! Sum of the temperature gradient of all neighbour elements, weighted by the volume of the neighbour element
    TempGradient = TempGradient+ ABS((RefTemp-NbTemp)/NbDistance)*ElemVolume_Shared(CnNbElem)
  END DO ! iNbElem

  ! 2a) Local Knudsen number of the density
  DensGradient = DensGradient/NbVolume
  Knudsen_Dens = MFP*DensGradient/RefDens

  ! 2b) Local Knudsen number of the velocity
  VeloGradient = VeloGradient/NbVolume
  Knudsen_Velo = MFP*VeloGradient/RefVelo

  ! 2c) Local Knudsen number of the temperature
  TempGradient = TempGradient/NbVolume
  IF (RefTemp.GT.0.) THEN
    Knudsen_Temp = MFP*TempGradient/RefTemp
  ELSE 
    Knudsen_Temp = 0.
  END IF

  ! Write out the non-equilibrium parameters
  CBC%OutputKnudsen(2,iElem) = Knudsen_Dens
  CBC%OutputKnudsen(3,iElem) = Knudsen_Velo
  CBC%OutputKnudsen(4,iElem) = Knudsen_Temp
  CBC%OutputKnudsen(5,iElem) = MAX(Knudsen_Dens, Knudsen_Velo, Knudsen_Temp)

  IF (MAX(Knudsen_Dens, Knudsen_Velo, Knudsen_Temp).GT.CBC%MaxLocalKnudsen) THEN
    CBC_DoDSMC = .TRUE.
  END IF

CASE('ThermNonEq')
  ! Calculate the density, velocity and temperature of the reference element
  CALL CalcCellProp(iElem, RefDens, RefVelo, RefTemp, SpecPartNum,Knudsen_NonEq,StressTensor,HeatVector)

  ! Total particle number in the element
  totalWeight = SUM(SpecPartNum)

  ! Write out the non-equilibrium parameters
  CBC%OutputKnudsen(6,iElem) = Knudsen_NonEq

  IF (Knudsen_NonEq.GT.CBC%MaxThermNonEq) THEN
    CBC_DoDSMC = .TRUE. 
  END IF

CASE('Combination')
  ! Set all needed variables to zero for the element in the loop
  DensGradient = 0.0; TempGradient = 0.0; VeloGradient = 0.0; NbVolume = 0.0

  ! Calculate the density, velocity and temperature of the reference element
  CALL CalcCellProp(iElem, RefDens, RefVelo, RefTemp, SpecPartNum,Knudsen_NonEq,StressTensor,HeatVector)
  ! Total particle number in the element
  totalWeight = SUM(SpecPartNum)

  MFP = CalcMeanFreePath(SpecPartNum,SUM(SpecPartNum),ElemVolume_Shared(CNElemID))

  ! Relative gradient between the cell and all neighboring cells for the particle density, velocity, translational temperature
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

    ! Calculate the density, velocity and temperature of the neighbour element
    CALL CalcDensVeloTemp(iElem, NbDens, NbVelo, NbTemp, NbSpecPartNum)
    NbtotalWeight = SUM(NbSpecPartNum)

    ! Sum of the density gradient of all neighbour elements, weighted by the volume of the neighbour element
    DensGradient = DensGradient + ABS((RefDens-NbDens)/NbDistance)*ElemVolume_Shared(CnNbElem)
    ! Sum of the velocity gradient of all neighbour elements, weighted by the volume of the neighbour element
    VeloGradient = VeloGradient + ABS((RefVelo-NbVelo)/NbDistance)*ElemVolume_Shared(CnNbElem)
    ! Sum of the temperature gradient of all neighbour elements, weighted by the volume of the neighbour element
    TempGradient = TempGradient+ ABS((RefTemp-NbTemp)/NbDistance)*ElemVolume_Shared(CnNbElem)
  END DO ! iNbElem

  DensGradient = DensGradient/NbVolume
  Knudsen_Dens = MFP*DensGradient/RefDens

  VeloGradient = VeloGradient/NbVolume
  Knudsen_Velo = MFP*VeloGradient/RefVelo

  TempGradient = TempGradient/NbVolume
  IF (RefTemp.GT.0.) THEN
    Knudsen_Temp = MFP*TempGradient/RefTemp
  ELSE 
    Knudsen_Temp = 0.
  END IF

  ! Write out the non-equilibrium parameters
  CBC%OutputKnudsen(2,iElem) = Knudsen_Dens
  CBC%OutputKnudsen(3,iElem) = Knudsen_Velo
  CBC%OutputKnudsen(4,iElem) = Knudsen_Temp
  CBC%OutputKnudsen(5,iElem) = MAX(Knudsen_Dens, Knudsen_Velo, Knudsen_Temp)
  CBC%OutputKnudsen(6,iElem) = Knudsen_NonEq

  IF (MAX(Knudsen_Dens, Knudsen_Velo, Knudsen_Temp).GT.CBC%MaxLocalKnudsen) THEN
    CBC_DoDSMC = .TRUE.
  ELSE IF (Knudsen_NonEq.GT.CBC%MaxThermNonEq) THEN
    CBC_DoDSMC = .TRUE. 
  END IF

CASE('ChapmanEnskog')
  CALL CalcCellProp(iElem, RefDens, RefVelo, RefTemp, SpecPartNum,Knudsen_NonEq,StressTensor,HeatVector)
  ! Write out the non-equilibrium parameters
  CBC%OutputKnudsen(7,iElem) = MAXVAL(StressTensor)
  CBC%OutputKnudsen(8,iElem) = MAXVAL(HeatVector)
  CBC%OutputKnudsen(9,iElem) = MAX(MAXVAL(HeatVector),MAXVAL(StressTensor))

  IF (MAXVAL(StressTensor).GT.(2.8*10.**4)) THEN
    CBC_DoDSMC = .TRUE.
  ELSE IF (MAXVAL(HeatVector).GT.(4.1*10.**7)) THEN
    CBC_DoDSMC = .TRUE.
  END IF

CASE('Output')
  ! Output Global Knudsen
  DensGradient = 0.0; TempGradient = 0.0; VeloGradient = 0.0; NbVolume = 0.0

  ! Calculate the density, velocity and temperature of the reference element
  CALL CalcCellProp(iElem, RefDens, RefVelo, RefTemp, SpecPartNum,Knudsen_NonEq,StressTensor,HeatVector)

  MFP = CalcMeanFreePath(SpecPartNum,SUM(SpecPartNum),ElemVolume_Shared(CNElemID))

  ! Write out the non-equilibrium parameters
  CBC%OutputKnudsen(1,iElem) = MFP/CBC%CharLength

  ! Total particle number in the element
  totalWeight = SUM(SpecPartNum)

  ! Relative gradient between the cell and all neighboring cells for the particle density, velocity, translational temperature
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

    ! Calculate the density, velocity and temperature of the reference element
    CALL CalcDensVeloTemp(iElem, NbDens, NbVelo, NbTemp, NbSpecPartNum)
    NbtotalWeight = SUM(NbSpecPartNum)

    ! Sum of the density gradient of all neighbour elements, weighted by the volume of the neighbour element
    DensGradient = DensGradient + ABS((RefDens-NbDens)/NbDistance)*ElemVolume_Shared(CnNbElem)
    ! Sum of the velocity gradient of all neighbour elements, weighted by the volume of the neighbour element
    VeloGradient = VeloGradient + ABS((RefVelo-NbVelo)/NbDistance)*ElemVolume_Shared(CnNbElem)
    ! Sum of the temperature gradient of all neighbour elements, weighted by the volume of the neighbour element
    TempGradient = TempGradient+ ABS((RefTemp-NbTemp)/NbDistance)*ElemVolume_Shared(CnNbElem)
  END DO ! iNbElem

  DensGradient = DensGradient/NbVolume
  Knudsen_Dens = MFP*DensGradient/RefDens

  VeloGradient = VeloGradient/NbVolume
  Knudsen_Velo = MFP*VeloGradient/RefVelo

  TempGradient = TempGradient/NbVolume
  IF (RefTemp.GT.0.) THEN
    Knudsen_Temp = MFP*TempGradient/RefTemp
  ELSE 
    Knudsen_Temp = 0.
  END IF

  ! Write out the non-equilibrium parameters
  CBC%OutputKnudsen(2,iElem) = Knudsen_Dens
  CBC%OutputKnudsen(3,iElem) = Knudsen_Velo
  CBC%OutputKnudsen(4,iElem) = Knudsen_Temp
  CBC%OutputKnudsen(5,iElem) = MAX(Knudsen_Dens, Knudsen_Velo, Knudsen_Temp)
  CBC%OutputKnudsen(6,iElem) = Knudsen_NonEq
  ! Write out the non-equilibrium parameters
  CBC%OutputKnudsen(7,iElem) = MAXVAL(StressTensor)
  CBC%OutputKnudsen(8,iElem) = MAXVAL(HeatVector)
  CBC%OutputKnudsen(9,iElem) = MAX(MAXVAL(HeatVector),MAXVAL(StressTensor))

END SELECT

RETURN

END FUNCTION CBC_DoDSMC

SUBROUTINE CalcCellProp(ElemID,Density,Velocity,Temperature,SpecPartNum,ThermNonEq,StressTensor,HeatVector)
!===================================================================================================================================
!> Calculate the breakdown properties for the reference cell
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
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: ElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)              :: Density, Velocity, Temperature, SpecPartNum(nSpecies),ThermNonEq, StressTensor(3,3), HeatVector(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart, iLoop, iSpec, iVal, jVal, iCoord, iPolyatMole, CNElemID
REAL                          :: Velo(nSpecies,3), Velo2(nSpecies,3), totalWeight, Temp(3), TempVib, TempTrans(3), TempRot
REAL                          :: EVib(nSpecies), ERot(nSpecies), SpecTemp(nSpecies), MolPartNum(nSpecies), BulkVelo(3)
REAL                          :: DoF_Rot(nSpecies), DoF_Vib(nSpecies), DoF_RotTotal, DoF_VibTotal, TVib_TempFac
REAL                          :: RelativeVelo(3), RelVeloTotal, StaticPressure, TempTotal
!-----------------------------------------------------------------------------------------------------------------------------------
SpecPartNum = 0.0; Velo = 0.0; Velo2 = 0.0; Velocity = 0.0; Temp = 0.0; Temperature = 0.0; TempTrans = 0.0; MolPartNum = 0.0
EVib = 0.0; ERot = 0.0; SpecTemp = 0.0; TempVib = 0.0; TempRot = 0.0; DoF_Rot = 0.0; DoF_RotTotal = 0.0; DoF_Vib = 0.0
DoF_VibTotal = 0.0; StressTensor = 0.0; HeatVector = 0.0; RelativeVelo = 0.0

! Reference element
CNElemID = GetCNElemID(ElemID+offSetElem)

! Determine the number of different species and the total particle number in the element by a loop over all particles 
iPart = PEM%pStart(ElemID)
DO iLoop = 1, PEM%pNumber(ElemID)
  SpecPartNum(PartSpecies(iPart)) = SpecPartNum(PartSpecies(iPart)) + GetParticleWeight(iPart)

  ! Particle velocities (squared) weighted by the particle number per element
  Velo(PartSpecies(iPart),1:3) = Velo(PartSpecies(iPart),1:3) + PartState(4:6,iPart) * GetParticleWeight(iPart)
  Velo2(PartSpecies(iPart),1:3) = Velo2(PartSpecies(iPart),1:3) + PartState(4:6,iPart)**2 * GetParticleWeight(iPart)

  ! Consideration of thermal non-equilibrium
  IF ((SpecDSMC(PartSpecies(iPart))%InterID.EQ.2).OR.(SpecDSMC(PartSpecies(iPart))%InterID.EQ.20)) THEN
    MolPartNum(PartSpecies(iPart)) = MolPartNum(PartSpecies(iPart)) + GetParticleWeight(iPart)
    EVib(PartSpecies(iPart)) = EVib(PartSpecies(iPart)) + (PartStateIntEn(1,iPart) - SpecDSMC(PartSpecies(iPart))%EZeroPoint) * GetParticleWeight(iPart)
    ERot(PartSpecies(iPart)) = ERot(PartSpecies(iPart)) + PartStateIntEn(2,iPart) * GetParticleWeight(iPart)
  END IF
  iPart = PEM%pNext(iPart)
END DO

! Total particle number in the element
totalWeight = SUM(SpecPartNum)

! Number density in the element = Particle number/volume
IF(usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.VarWeighting%DoVariableWeighting) THEN
  Density = totalWeight / ElemVolume_Shared(CNElemID)
ELSE
  Density = totalWeight * Species(1)%MacroParticleFactor / ElemVolume_Shared(CNElemID)
END IF

! Calculation of the total mean velocity under consideration of all directions
DO iLoop = 1, 3
  Velocity = Velocity + (SUM(Velo(:,iLoop))/totalWeight)**2
  BulkVelo(iLoop) = SUM(Velo(:,iLoop))/totalWeight
END DO
Velocity = SQRT(Velocity)

! Calculation of the total mean translational temperature
DO iSpec = 1, nSpecies
  IF (SpecPartNum(iSpec).GE.1.0) THEN 
    Temp  = Species(iSpec)%MassIC/BoltzmannConst * ((Velo2(iSpec,:)/SpecPartNum(iSpec)) - (Velo(iSpec,:)/SpecPartNum(iSpec))**2)
    Temperature = Temperature + ((Temp(1) + Temp(2) + Temp(3)) / 3.)*SpecPartNum(iSpec)

    ! Translational temperature in x/y/z direction
    TempTrans = TempTrans + Temp*SpecPartNum(iSpec)
    ! Translational temperature for each species
    SpecTemp(iSpec) = (Temp(1) + Temp(2) + Temp(3)) / 3.

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
IF (Temperature.GT.0.) THEN
  Temperature = Temperature/totalweight
ELSE 
  Temperature = 0.
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
TempTotal = (Temperature*3. + DoF_VibTotal*TempVib + DoF_RotTotal*TempRot) / (3. + DoF_VibTotal + DoF_RotTotal)

ThermNonEq = SQRT(((TempTrans(1)-TempTotal)**2 + (TempTrans(2)-TempTotal)**2 + (TempTrans(3)-TempTotal)**2 &
            + DoF_RotTotal*(TempRot-TempTotal)**2 + DoF_VibTotal*(TempVib-TempTotal)**2 )/((3. + DoF_VibTotal + DoF_RotTotal)*TempTotal**2))

iPart = PEM%pStart(ElemID)
DO iLoop = 1, PEM%pNumber(ElemID)
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

IF(PEM%pNumber(ElemID).GT.2) THEN
  HeatVector = HeatVector*PEM%pNumber(ELemID)*PEM%pNumber(ElemID)/((PEM%pNumber(ElemID)-1.)*(PEM%pNumber(ElemID)-2.)*totalWeight)
ELSE
  HeatVector = 0.0
END IF

StressTensor = -(StressTensor/totalWeight)

END SUBROUTINE CalcCellProp
  
SUBROUTINE CalcDensVeloTemp(ElemID,Density,Velocity,Temperature,SpecPartNum)
!===================================================================================================================================
!> Calculate the density, velocity and translational temperature for the neighbour cells
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
USE MOD_DSMC_Vars              ,ONLY: DSMC, SpecDSMC, CollisMode, PartStateIntEn, RadialWeighting, VarWeighting
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: ElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)              :: Density, Velocity, Temperature, SpecPartNum(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart, iLoop, iSpec, iCoord, CNElemID
REAL                          :: Velo(nSpecies,3), Velo2(nSpecies,3), totalWeight, Temp(3)
!-----------------------------------------------------------------------------------------------------------------------------------
SpecPartNum = 0.0; Velo = 0.0; Velo2 = 0.0; Velocity = 0.0; Temp = 0.0; Temperature = 0.0

! Reference element
CNElemID = GetCNElemID(ElemID+offSetElem)

! Determine the number of different species and the total particle number in the element by a loop over all particles 
iPart = PEM%pStart(ElemID)
DO iLoop = 1, PEM%pNumber(ElemID)
  SpecPartNum(PartSpecies(iPart)) = SpecPartNum(PartSpecies(iPart)) + GetParticleWeight(iPart)

  ! Particle velocities (squared) weighted by the particle number per element
  Velo(PartSpecies(iPart),1:3) = Velo(PartSpecies(iPart),1:3) + PartState(4:6,iPart) * GetParticleWeight(iPart)
  Velo2(PartSpecies(iPart),1:3) = Velo2(PartSpecies(iPart),1:3) + PartState(4:6,iPart)**2 * GetParticleWeight(iPart)
  iPart = PEM%pNext(iPart)
END DO

! Total particle number in the element
totalWeight = SUM(SpecPartNum)

! Number density in the element = Particle number/volume
IF(usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.VarWeighting%DoVariableWeighting) THEN
  Density = totalWeight / ElemVolume_Shared(CNElemID)
ELSE
  Density = totalWeight * Species(1)%MacroParticleFactor / ElemVolume_Shared(CNElemID)
END IF

! Calculation of the total mean velocity under consideration of all directions
DO iLoop = 1, 3
  Velocity = Velocity + (SUM(Velo(:,iLoop))/totalWeight)**2
END DO
Velocity = SQRT(Velocity)

! Calculation of the total mean translational temperature
DO iSpec = 1, nSpecies
  IF (SpecPartNum(iSpec).GE.1.0) THEN 
    Temp  = Species(iSpec)%MassIC/BoltzmannConst * ((Velo2(iSpec,:)/SpecPartNum(iSpec)) - (Velo(iSpec,:)/SpecPartNum(iSpec))**2)
    Temperature = Temperature + ((Temp(1) + Temp(2) + Temp(3)) / 3.)*SpecPartNum(iSpec)
  END IF
END DO  
IF (Temperature.GT.0.) THEN
  Temperature = Temperature/totalweight
ELSE 
  Temperature = 0.
END IF
END SUBROUTINE CalcDensVeloTemp

END MODULE MOD_BGK_DSMC_Coupling
