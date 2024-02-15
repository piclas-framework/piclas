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

MODULE MOD_FP_DSMC_Coupling
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
PUBLIC :: FP_CBC_DoDSMC
!===================================================================================================================================

CONTAINS

LOGICAL FUNCTION FP_CBC_DoDSMC(iElem)
!===================================================================================================================================
!> Test of different continuum-breakdown criteria for the coupling between FP and DSMC
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_FPFlow_Vars
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
USE MOD_DSMC_Analyze           ,ONLY: CalcMeanFreePath
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_Particle_Vars          ,ONLY: PEM, nSpecies, PartSpecies
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemVolume_Shared, ElemToElemMapping, ElemToElemInfo, ElemMidPoint_Shared
USE MOD_Mesh_Vars              ,ONLY: offsetElem, nElems
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID, GetGlobalElemID
USE MOD_part_tools             ,ONLY: GetParticleWeight
USE MOD_BGK_DSMC_Coupling      ,ONLY: CalcCellProp, CalcDensVeloTemp, CalcCellProp_Samp, CalcDensVeloTemp_Samp
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
REAL                         :: HeatVector(3), StressTensor(3,3) !, ChapmanEnskog
!===================================================================================================================================
! Default case: use of FP
FP_CBC_DoDSMC = .FALSE.

! Reference element
CNElemID = GetCNElemID(iElem+offSetElem)

! Definition of continuum-breakdown in the cell and swtich between BGK and DSMC
SELECT CASE (TRIM(FP_CBC%SwitchCriterium))
! Global Knudsen number: Mean free path / characteristic length
CASE('GlobalKnudsen')
  SpecPartNum = 0.0

  iPart = PEM%pStart(iElem)
  DO iLoop = 1, PEM%pNumber(iElem)
    SpecPartNum(PartSpecies(iPart)) = SpecPartNum(PartSpecies(iPart)) + GetParticleWeight(iPart)
    iPart = PEM%pNext(iPart)
  END DO

  MFP = CalcMeanFreePath(SpecPartNum,SUM(SpecPartNum),ElemVolume_Shared(CNElemID))

  ! Write out the non-equilibrium parameters
  FP_CBC%OutputKnudsen(1,iElem) = MFP/FP_CBC%CharLength

   IF ((MFP/FP_CBC%CharLength).GE.FP_CBC%MaxGlobalKnudsen) THEN
   FP_CBC_DoDSMC = .TRUE.
  END IF

! Local Knudsen number: gradient of the flow-field values
CASE('LocalKnudsen')
  ! Set all needed variables to zero for the element in the loop
  DensGradient = 0.0; TempGradient = 0.0; VeloGradient = 0.0; NbVolume = 0.0

  ! Calculate the density, velocity and temperature of the reference element
  IF ((DSMC%SampNum.GT.0.).AND.FP_CBC%AverageSamp) THEN
    CALL CalcDensVeloTemp_Samp(iElem, RefDens, RefVelo, RefTemp, SpecPartNum)
  ELSE
    CALL CalcDensVeloTemp(iElem, RefDens, RefVelo, RefTemp, SpecPartNum)
  END IF

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
    IF ((DSMC%SampNum.GT.0.).AND.FP_CBC%AverageSamp) THEN
      CALL CalcDensVeloTemp_Samp(iElem, NbDens, NbVelo, NbTemp, NbSpecPartNum)
    ELSE
      CALL CalcDensVeloTemp(iElem, NbDens, NbVelo, NbTemp, NbSpecPartNum)
    END IF
    NbtotalWeight = SUM(NbSpecPartNum)

    ! Sum of the density gradient of all neighbour elements, weighted by the volume of the neighbour element
    DensGradient = DensGradient + ABS(RefDens-NbDens/NbDistance)*ElemVolume_Shared(CnNbElem)
    VeloGradient = VeloGradient + ABS(RefVelo-NbVelo/NbDistance)*ElemVolume_Shared(CnNbElem)
    TempGradient = TempGradient+ ABS(RefTemp-NbTemp/NbDistance)*ElemVolume_Shared(CnNbElem)
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
  FP_CBC%OutputKnudsen(2,iElem) = Knudsen_Dens
  FP_CBC%OutputKnudsen(3,iElem) = Knudsen_Velo
  FP_CBC%OutputKnudsen(4,iElem) = Knudsen_Temp
  FP_CBC%OutputKnudsen(5,iElem) = MAX(Knudsen_Dens, Knudsen_Velo, Knudsen_Temp)

  ! Test if the volume weighted gradient is larger than a predefined value
  IF (MAX(Knudsen_Dens, Knudsen_Velo, Knudsen_Temp).GT.FP_CBC%MaxLocalKnudsen) THEN
   FP_CBC_DoDSMC = .TRUE.
  END IF

! Thermal non-equilibrium: deviation of the individual temperatures from the equilibrium value in the cell
CASE('ThermNonEq')
  ! Calculate the density, velocity and temperature of the reference element
  IF ((DSMC%SampNum.GT.0.).AND.FP_CBC%AverageSamp) THEN
    CALL CalcCellProp_Samp(iElem, RefDens, RefVelo, RefTemp, SpecPartNum,Knudsen_NonEq,StressTensor,HeatVector)
  ELSE
    CALL CalcCellProp(iElem, RefDens, RefVelo, RefTemp, SpecPartNum,Knudsen_NonEq,StressTensor,HeatVector)
  END IF

  ! Total particle number in the element
  totalWeight = SUM(SpecPartNum)

  ! Write out the non-equilibrium parameters
  FP_CBC%OutputKnudsen(6,iElem) = Knudsen_NonEq

  IF (Knudsen_NonEq.GT.FP_CBC%MaxThermNonEq) THEN
   FP_CBC_DoDSMC = .TRUE.
  END IF

! Combination of the local Knudsen number and thermal non-equilibrium
CASE('Combination')
  ! Set all needed variables to zero for the element in the loop
  DensGradient = 0.0; TempGradient = 0.0; VeloGradient = 0.0; NbVolume = 0.0

  ! Calculate the density, velocity and temperature of the reference element
  IF ((DSMC%SampNum.GT.0.).AND.FP_CBC%AverageSamp) THEN
    CALL CalcCellProp_Samp(iElem, RefDens, RefVelo, RefTemp, SpecPartNum,Knudsen_NonEq,StressTensor,HeatVector)
  ELSE
    CALL CalcCellProp(iElem, RefDens, RefVelo, RefTemp, SpecPartNum,Knudsen_NonEq,StressTensor,HeatVector)
  END IF
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
    IF ((DSMC%SampNum.GT.0.).AND.FP_CBC%AverageSamp) THEN
      CALL CalcDensVeloTemp_Samp(iElem, NbDens, NbVelo, NbTemp, NbSpecPartNum)
    ELSE
      CALL CalcDensVeloTemp(iElem, NbDens, NbVelo, NbTemp, NbSpecPartNum)
    END IF
    NbtotalWeight = SUM(NbSpecPartNum)

    ! Sum of the density gradient of all neighbour elements, weighted by the volume of the neighbour element
    DensGradient = DensGradient + ABS(RefDens-NbDens/NbDistance)*ElemVolume_Shared(CnNbElem)
    VeloGradient = VeloGradient + ABS(RefVelo-NbVelo/NbDistance)*ElemVolume_Shared(CnNbElem)
    TempGradient = TempGradient+ ABS(RefTemp-NbTemp/NbDistance)*ElemVolume_Shared(CnNbElem)
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
  FP_CBC%OutputKnudsen(2,iElem) = Knudsen_Dens
  FP_CBC%OutputKnudsen(3,iElem) = Knudsen_Velo
  FP_CBC%OutputKnudsen(4,iElem) = Knudsen_Temp
  FP_CBC%OutputKnudsen(5,iElem) = MAX(Knudsen_Dens, Knudsen_Velo, Knudsen_Temp)
  FP_CBC%OutputKnudsen(6,iElem) = Knudsen_NonEq

  IF (MAX(Knudsen_Dens, Knudsen_Velo, Knudsen_Temp).GT.FP_CBC%MaxLocalKnudsen) THEN
   FP_CBC_DoDSMC = .TRUE.
  ELSE IF (Knudsen_NonEq.GT.FP_CBC%MaxThermNonEq) THEN
   FP_CBC_DoDSMC = .TRUE.
  END IF

! Chapman-Enskog parameter: maximum value of the heat flux vector and the shear stress tensor
CASE('ChapmanEnskog')
  IF ((DSMC%SampNum.GT.0.).AND.FP_CBC%AverageSamp) THEN
    CALL CalcCellProp_Samp(iElem, RefDens, RefVelo, RefTemp, SpecPartNum,Knudsen_NonEq,StressTensor,HeatVector)
  ELSE
    CALL CalcCellProp(iElem, RefDens, RefVelo, RefTemp, SpecPartNum,Knudsen_NonEq,StressTensor,HeatVector)
  END IF

  ! Store the maximum heat flux values for the next iteration
  FP_CBC%Max_HeatVec(iElem) = MAXVAL(HeatVector)
  FP_CBC%Max_StressTens(iElem) = MAXVAL(StressTensor)

  IF (MAXVAL(StressTensor).GT.(FP_CBC%MaxChapmanEnskog*MAXVAL(FP_CBC%Max_StressTens))) THEN
   FP_CBC_DoDSMC = .TRUE.
  ELSE IF (MAXVAL(HeatVector).GT.(FP_CBC%MaxChapmanEnskog*MAXVAL(FP_CBC%Max_HeatVec))) THEN
   FP_CBC_DoDSMC = .TRUE.
  END IF

  ! Write out the non-equilibrium parameters
  FP_CBC%OutputKnudsen(7,iElem) = MAXVAL(StressTensor)/MAXVAL(FP_CBC%Max_StressTens)
  FP_CBC%OutputKnudsen(8,iElem) = MAXVAL(HeatVector)/MAXVAL(FP_CBC%Max_HeatVec)

! Calculation and output of all possible coupling-criteria, use of only BGK or FP in the calculation
CASE('Output')
  ! Output Global Knudsen
  DensGradient = 0.0; TempGradient = 0.0; VeloGradient = 0.0; NbVolume = 0.0

  ! Calculate the density, velocity and temperature of the reference element
  IF ((DSMC%SampNum.GT.0.).AND.FP_CBC%AverageSamp) THEN
    CALL CalcCellProp_Samp(iElem, RefDens, RefVelo, RefTemp, SpecPartNum,Knudsen_NonEq,StressTensor,HeatVector)
  ELSE
    CALL CalcCellProp(iElem, RefDens, RefVelo, RefTemp, SpecPartNum,Knudsen_NonEq,StressTensor,HeatVector)
  END IF

  MFP = CalcMeanFreePath(SpecPartNum,SUM(SpecPartNum),ElemVolume_Shared(CNElemID))

  ! Global Knudsen number
  FP_CBC%OutputKnudsen(1,iElem) = MFP/FP_CBC%CharLength

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
    IF ((DSMC%SampNum.GT.0.).AND.FP_CBC%AverageSamp) THEN
      CALL CalcDensVeloTemp_Samp(iElem, NbDens, NbVelo, NbTemp, NbSpecPartNum)
    ELSE
      CALL CalcDensVeloTemp(iElem, NbDens, NbVelo, NbTemp, NbSpecPartNum)
    END IF
    NbtotalWeight = SUM(NbSpecPartNum)

    ! Sum of the density gradient of all neighbour elements, weighted by the volume of the neighbour element
    DensGradient = DensGradient + ABS(RefDens-NbDens/NbDistance)*ElemVolume_Shared(CnNbElem)
    VeloGradient = VeloGradient + ABS(RefVelo-NbVelo/NbDistance)*ElemVolume_Shared(CnNbElem)
    TempGradient = TempGradient+ ABS(RefTemp-NbTemp/NbDistance)*ElemVolume_Shared(CnNbElem)
  END DO ! iNbElem

  ! Local Knudsen number of the density, velocity and temperature
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

  ! Store the maximum heat flux values for the next iteration
  FP_CBC%Max_HeatVec(iElem) = MAXVAL(HeatVector)
  FP_CBC%Max_StressTens(iElem) = MAXVAL(StressTensor)

  ! Write out the non-equilibrium parameters
  FP_CBC%OutputKnudsen(2,iElem) = Knudsen_Dens
  FP_CBC%OutputKnudsen(3,iElem) = Knudsen_Velo
  FP_CBC%OutputKnudsen(4,iElem) = Knudsen_Temp
  FP_CBC%OutputKnudsen(5,iElem) = MAX(Knudsen_Dens, Knudsen_Velo, Knudsen_Temp)
  FP_CBC%OutputKnudsen(6,iElem) = Knudsen_NonEq
  FP_CBC%OutputKnudsen(7,iElem) = MAXVAL(StressTensor)/MAXVAL(FP_CBC%Max_StressTens)
  FP_CBC%OutputKnudsen(8,iElem) = MAXVAL(HeatVector)/MAXVAL(FP_CBC%Max_HeatVec)

END SELECT

RETURN

END FUNCTION FP_CBC_DoDSMC

END MODULE MOD_FP_DSMC_Coupling
