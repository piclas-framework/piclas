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

MODULE MOD_DSMC_ElectronicModel
!===================================================================================================================================
! module including °qk procedures
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE ElectronicEnergyExchange
  MODULE PROCEDURE ElectronicEnergyExchange
END INTERFACE

INTERFACE InitElectronShell
  MODULE PROCEDURE InitElectronShell
END INTERFACE

INTERFACE ReadSpeciesLevel
  MODULE PROCEDURE ReadSpeciesLevel
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: ElectronicEnergyExchange, InitElectronShell, TVEEnergyExchange, ReadSpeciesLevel, CalcProbCorrFactorElec
PUBLIC :: RelaxElectronicShellWall, LT_ElectronicEnergyExchange, LT_ElectronicExc_ConstructPartList, LT_ElectronicEnergyExchangeChem
!===================================================================================================================================
CONTAINS

SUBROUTINE InitElectronShell(iSpec,iPart,iInit,init_or_sf)
!===================================================================================================================================
!> Set the initial electronic energy of a particle depending on its species and temperature input (through init or surface flux)
!===================================================================================================================================
USE MOD_Globals         ,ONLY: abort
USE MOD_Globals_Vars    ,ONLY: BoltzmannConst
USE MOD_Particle_Vars   ,ONLY: PEM
USE MOD_DSMC_Vars       ,ONLY: SpecDSMC, PartStateIntEn, ElectronicDistriPart, DSMC, BGGas
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: iPart, iSpec, iInit, init_or_sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iQua, LocalElemID
REAL                    :: iRan, ElectronicPartition, ElectronicPartitionTemp, iRan2, tmpExp
REAL                    :: Telec                      ! electronic temperature
!===================================================================================================================================
SELECT CASE (init_or_sf)
CASE(1) ! Initial
  TElec=SpecDSMC(iSpec)%Init(iInit)%TElec
CASE(2) ! SurfaceFlux
  TElec=SpecDSMC(iSpec)%Surfaceflux(iInit)%Telec
CASE DEFAULT
  CALL abort(__STAMP__,'neither iInit nor Surfaceflux defined as reference!')
END SELECT

! Background gas distribution
IF(BGGas%NumberOfSpecies.GT.0) THEN
  IF(BGGas%BackgroundSpecies(iSpec).AND.BGGas%UseDistribution) THEN
    LocalElemID = PEM%LocalElemID(iPart)
    TElec = BGGas%Distribution(BGGas%MapSpecToBGSpec(iSpec),DSMC_TELEC,LocalElemID)
  END IF
END IF

! Check if the temperature is zero
IF(TElec.LE.0.)THEN
  PartStateIntEn(3,iPart) = 0.0
  RETURN
END IF ! TElec.LE.0.

ElectronicPartition  = 0.
ElectronicPartitionTemp = 0.
! calculate sum over all energy levels == partition function for temperature Telec
SELECT CASE(DSMC%ElectronicModel)
CASE(1,4)
  DO iQua = 0, SpecDSMC(iSpec)%MaxElecQuant - 1
    tmpExp = SpecDSMC(iSpec)%ElectronicState(2,iQua) / Telec
    IF (CHECKEXP(tmpExp)) THEN
      ElectronicPartitionTemp = SpecDSMC(iSpec)%ElectronicState(1,iQua) *EXP(-tmpExp)
    ELSE
      ElectronicPartitionTemp = 0.0
    END IF
    IF ( ElectronicPartitionTemp .GT. ElectronicPartition ) THEN
      ElectronicPartition = ElectronicPartitionTemp
    END IF
  END DO
  ElectronicPartitionTemp = 0.
  ! select level
  CALL RANDOM_NUMBER(iRan2)
  DO WHILE ( iRan2 .GE. ElectronicPartitionTemp / ElectronicPartition )
    CALL RANDOM_NUMBER(iRan)
    iQua = int( ( SpecDSMC(iSpec)%MaxElecQuant ) * iRan)
    tmpExp = SpecDSMC(iSpec)%ElectronicState(2,iQua) / Telec
    IF (CHECKEXP(tmpExp)) THEN
      ElectronicPartitionTemp = SpecDSMC(iSpec)%ElectronicState(1,iQua) *EXP(-tmpExp)
    ELSE
      ElectronicPartitionTemp = 0.0
    END IF
    CALL RANDOM_NUMBER(iRan2)
  END DO
  PartStateIntEn(3,iPart) = BoltzmannConst * SpecDSMC(iSpec)%ElectronicState(2,iQua)
CASE(2)
  IF(ALLOCATED(ElectronicDistriPart(iPart)%DistriFunc)) DEALLOCATE(ElectronicDistriPart(iPart)%DistriFunc)
  ALLOCATE(ElectronicDistriPart(iPart)%DistriFunc(1:SpecDSMC(iSpec)%MaxElecQuant))
  PartStateIntEn(3,iPart) = 0.0
  DO iQua = 0, SpecDSMC(iSpec)%MaxElecQuant - 1
    tmpExp = SpecDSMC(iSpec)%ElectronicState(2,iQua) / Telec
    IF (CHECKEXP(tmpExp)) &
      ElectronicPartition = ElectronicPartition + SpecDSMC(iSpec)%ElectronicState(1,iQua) * EXP(-tmpExp)
  END DO
  DO iQua = 0, SpecDSMC(iSpec)%MaxElecQuant - 1
    tmpExp = SpecDSMC(iSpec)%ElectronicState(2,iQua) / Telec
    IF (CHECKEXP(tmpExp)) THEN
      ElectronicDistriPart(iPart)%DistriFunc(iQua+1) = SpecDSMC(iSpec)%ElectronicState(1,iQua)*EXP(-tmpExp)/ElectronicPartition
    ELSE
      ElectronicDistriPart(iPart)%DistriFunc(iQua+1) = 0.0
    END IF
    PartStateIntEn(3,iPart) = PartStateIntEn(3,iPart) + &
        ElectronicDistriPart(iPart)%DistriFunc(iQua+1) * BoltzmannConst * SpecDSMC(iSpec)%ElectronicState(2,iQua)
  END DO
CASE(3)
  ! Initialize in ground state
  PartStateIntEn(3,iPart) = 0.
CASE DEFAULT
  CALL abort(__STAMP__,'ERROR: Unknown electronic relaxation model: ',IntInfoOpt=DSMC%ElectronicModel)
END SELECT

END SUBROUTINE InitElectronShell


FUNCTION RelaxElectronicShellWall(iPart,TWall)
!===================================================================================================================================
!> Function to determine the new electronic state of a particle at the wall temperature
!===================================================================================================================================
USE MOD_Globals
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars             ,ONLY: SpecDSMC, DSMC, ElectronicDistriPart
USE MOD_Particle_Vars         ,ONLY: PartSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPart
REAL, INTENT(IN)              :: TWall
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                          :: RelaxElectronicShellWall
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iQua, iSpec
REAL                          :: iRan, ElectronicPartition, ElectronicPartitionTemp, iRan2, TempRatio, tmpExp
!===================================================================================================================================
ElectronicPartition  = 0.
ElectronicPartitionTemp = 0.
RelaxElectronicShellWall = 0.
iSpec = PartSpecies(iPart)

SELECT CASE(DSMC%ElectronicModel)
CASE(1,4)
  ! calculate sum over all energy levels == partition function for temperature Telec
  DO iQua = 0, SpecDSMC(iSpec)%MaxElecQuant - 1
    TempRatio = SpecDSMC(iSpec)%ElectronicState(2,iQua) / TWall
    IF(CHECKEXP(TempRatio)) THEN
      ElectronicPartitionTemp = SpecDSMC(iSpec)%ElectronicState(1,iQua) * EXP (-TempRatio)
      IF ( ElectronicPartitionTemp .GT. ElectronicPartition ) THEN
        ElectronicPartition = ElectronicPartitionTemp
      END IF
    END IF
  END DO
  CALL RANDOM_NUMBER(iRan)
  iQua = INT(SpecDSMC(iSpec)%MaxElecQuant*iRan)
  TempRatio = SpecDSMC(iSpec)%ElectronicState(2,iQua) / TWall
  IF(CHECKEXP(TempRatio)) THEN
    ElectronicPartitionTemp = SpecDSMC(iSpec)%ElectronicState(1,iQua) * EXP (-TempRatio)
  ELSE
    ElectronicPartitionTemp = 0.
  END IF
  ! select level
  CALL RANDOM_NUMBER(iRan2)
  DO WHILE ( iRan2 .GE. ElectronicPartitionTemp / ElectronicPartition )
    CALL RANDOM_NUMBER(iRan)
    iQua = int( ( SpecDSMC(iSpec)%MaxElecQuant ) * iRan)
    TempRatio = SpecDSMC(iSpec)%ElectronicState(2,iQua) / TWall
    IF(CHECKEXP(TempRatio)) THEN
      ElectronicPartitionTemp  = SpecDSMC(iSpec)%ElectronicState(1,iQua) * EXP(-TempRatio)
      CALL RANDOM_NUMBER(iRan2)
    END IF
  END DO
  RelaxElectronicShellWall = BoltzmannConst * SpecDSMC(iSpec)%ElectronicState(2,iQua)
CASE(2)
  RelaxElectronicShellWall = 0.0
  DO iQua = 0, SpecDSMC(iSpec)%MaxElecQuant - 1
    tmpExp = SpecDSMC(iSpec)%ElectronicState(2,iQua) / TWall
    IF (CHECKEXP(tmpExp)) &
      ElectronicPartition = ElectronicPartition + SpecDSMC(iSpec)%ElectronicState(1,iQua) * EXP(-tmpExp)
  END DO
  DO iQua = 0, SpecDSMC(iSpec)%MaxElecQuant - 1
    tmpExp = SpecDSMC(iSpec)%ElectronicState(2,iQua) / TWall
    IF (CHECKEXP(tmpExp)) THEN
      ElectronicDistriPart(iPart)%DistriFunc(iQua+1) = SpecDSMC(iSpec)%ElectronicState(1,iQua) * EXP(-tmpExp)/ElectronicPartition
    ELSE
      ElectronicDistriPart(iPart)%DistriFunc(iQua+1) = 0.0
    END IF
    RelaxElectronicShellWall = RelaxElectronicShellWall+ &
        ElectronicDistriPart(iPart)%DistriFunc(iQua+1) * BoltzmannConst * SpecDSMC(iSpec)%ElectronicState(2,iQua)
  END DO
CASE(3)
  ! Reset particle to ground state at wall
  RelaxElectronicShellWall = 0.0
CASE DEFAULT
  CALL abort(__STAMP__,'ERROR: Unknown electronic relaxation model: ',IntInfoOpt=DSMC%ElectronicModel)
END SELECT

END FUNCTION RelaxElectronicShellWall


SUBROUTINE ElectronicEnergyExchange(iPair,iPart1,FakXi, NewPart, XSec_Level)
!===================================================================================================================================
!> Electronic energy exchange:
!> Model 1 (Liechty): Simulation particle has a specific electronic energy level
!> Model 2 (Burt): Simulation particle has an electronic energy distribution function attached to it
!> Model 3: Cross-section based excitation
!===================================================================================================================================
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst, ElementaryCharge
USE MOD_DSMC_Vars              ,ONLY: SpecDSMC, PartStateIntEn, RadialWeighting, Coll_pData, DSMC, ElectronicDistriPart, CollInf
USE MOD_Particle_Vars          ,ONLY: PartSpecies, UseVarTimeStep, usevMPF, nSpecies
USE MOD_part_tools             ,ONLY: GetParticleWeight
USE MOD_Particle_Analyze_Tools ,ONLY: CalcTelec
USE MOD_MCC_Vars                ,ONLY: SpecXSec
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPair, iPart1
REAL, INTENT(IN)              :: FakXi
LOGICAL, INTENT(IN),OPTIONAL  :: NewPart
INTEGER, INTENT(IN),OPTIONAL  :: XSec_Level
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iQuaMax, MaxElecQuant, iQua, iSpec, iCase
REAL                          :: iRan, iRan2, gmax, gtemp, PartStateTemp, CollisionEnergy, ETraRel, TransElec, ElectronicPartition
REAL                          :: Eold, Etmp, tmpExp, LocRelaxProb
REAL,ALLOCATABLE              :: DistriOld(:)
!===================================================================================================================================
iSpec = PartSpecies(iPart1)
SELECT CASE(DSMC%ElectronicModel)
CASE(1)
  IF (usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep) THEN
    CollisionEnergy = Coll_pData(iPair)%Ec / GetParticleWeight(iPart1)
  ELSE
    CollisionEnergy = Coll_pData(iPair)%Ec
  END IF

  iQuaMax  = 0
  ! Determine max electronic quant
  MaxElecQuant = SpecDSMC(iSpec)%MaxElecQuant - 1
  ! determine maximal Quant and term according to Eq (7) of Liechty
  gmax = 0.
  PartStateTemp = CollisionEnergy / BoltzmannConst
  DO iQua = 0, MaxElecQuant
    IF (PartStateTemp - SpecDSMC(iSpec)%ElectronicState(2,iQua).GT.0.) THEN
      gtemp = SpecDSMC(iSpec)%ElectronicState(1,iQua) * &
              ( CollisionEnergy - BoltzmannConst * SpecDSMC(iSpec)%ElectronicState(2,iQua))**FakXi
      ! maximal possible Quant before term goes negative
      iQuaMax = iQua
      IF ( gtemp .GT. gmax ) THEN
      ! Quant of largest value of Eq (7)
        gmax = gtemp
      END IF
    ELSE
      EXIT
    END IF
  END DO
  IF(gmax.LE.0.0) THEN
    PartStateIntEn(3,iPart1) = 0.0
    RETURN
  END IF
  CALL RANDOM_NUMBER(iRan)
  iQua = int( ( iQuaMax +1 ) * iRan)
  gtemp = SpecDSMC(iSpec)%ElectronicState(1,iQua) * &
          ( CollisionEnergy - BoltzmannConst * SpecDSMC(iSpec)%ElectronicState(2,iQua))**FakXi
  CALL RANDOM_NUMBER(iRan2)
  ! acceptance-rejection for iQuaElec
  DO WHILE ( iRan2 .GE. gtemp / gmax )
    CALL RANDOM_NUMBER(iRan)
    iQua = int( ( iQuaMax +1 ) * iRan)
    gtemp = SpecDSMC(iSpec)%ElectronicState(1,iQua) * &
            ( CollisionEnergy - BoltzmannConst * SpecDSMC(iSpec)%ElectronicState(2,iQua))**FakXi
    CALL RANDOM_NUMBER(iRan2)
  END DO
  PartStateIntEn(3,iPart1) = BoltzmannConst * SpecDSMC(iSpec)%ElectronicState(2,iQua)
CASE(2)
  IF (PRESENT(NewPart)) THEN
    LocRelaxProb = 1.0
  ELSE
    LocRelaxProb = SpecDSMC(iSpec)%ElecRelaxProb
  END IF
  Eold=  PartStateIntEn(3,iPart1)
  DistriOld = ElectronicDistriPart(iPart1)%DistriFunc
  ETraRel = Coll_pData(iPair)%Ec
  IF (usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep) THEN
    ETraRel = ETraRel / GetParticleWeight(iPart1)
  END IF
  TransElec = DSMC%InstantTransTemp(nSpecies + 1)
  IF (TransElec.LE.0.0) TransElec = 1./(BoltzmannConst*(FakXi+1.))*ETraRel
  ElectronicPartition = 0.0
  DO iQua = 0, SpecDSMC(iSpec)%MaxElecQuant - 1
    tmpExp = SpecDSMC(iSpec)%ElectronicState(2,iQua) / TransElec
    IF (CHECKEXP(tmpExp)) &
      ElectronicPartition = ElectronicPartition + SpecDSMC(iSpec)%ElectronicState(1,iQua) * EXP (-tmpExp)
  END DO
  PartStateIntEn(3,iPart1) = 0.0
  DO iQua = 0, SpecDSMC(iSpec)%MaxElecQuant - 1
    tmpExp = SpecDSMC(iSpec)%ElectronicState(2,iQua) / TransElec
    IF (CHECKEXP(tmpExp)) THEN
      ElectronicDistriPart(iPart1)%DistriFunc(iQua+1) = &
        (1.-LocRelaxProb)*ElectronicDistriPart(iPart1)%DistriFunc(iQua+1) + &
        LocRelaxProb * SpecDSMC(iSpec)%ElectronicState(1,iQua) *EXP (-tmpExp)/ElectronicPartition
    ELSE
      ElectronicDistriPart(iPart1)%DistriFunc(iQua+1) = (1.-LocRelaxProb)*ElectronicDistriPart(iPart1)%DistriFunc(iQua+1)
    END IF
!      ElectronicDistriPart(iPart1)%DistriFunc(iQua+1) =  SpecDSMC(iSpec)%ElectronicState(1,iQua) * &
!              EXP ( - SpecDSMC(iSpec)%ElectronicState(2,iQua) / TransElec)/ElectronicPartition
    PartStateIntEn(3,iPart1) = PartStateIntEn(3,iPart1) + &
        ElectronicDistriPart(iPart1)%DistriFunc(iQua+1) * BoltzmannConst * SpecDSMC(iSpec)%ElectronicState(2,iQua)
  END DO
  IF ((Coll_pData(iPair)%Ec-PartStateIntEn(3,iPart1)*GetParticleWeight(iPart1)).LT.0.0) then
    Etmp = (Coll_pData(iPair)%Ec - (1.-LocRelaxProb)*Eold*GetParticleWeight(iPart1))/(GetParticleWeight(iPart1)*LocRelaxProb)
    TransElec = CalcTelec(Etmp, iSpec)*0.98
    ElectronicPartition = 0.0
    DO iQua = 0, SpecDSMC(iSpec)%MaxElecQuant - 1
      tmpExp = SpecDSMC(iSpec)%ElectronicState(2,iQua) / TransElec
      IF (CHECKEXP(tmpExp)) &
        ElectronicPartition = ElectronicPartition + SpecDSMC(iSpec)%ElectronicState(1,iQua) * EXP (-tmpExp)
    END DO
    PartStateIntEn(3,iPart1) = 0.0
    DO iQua = 0, SpecDSMC(iSpec)%MaxElecQuant - 1
      tmpExp = SpecDSMC(iSpec)%ElectronicState(2,iQua) / TransElec
      IF (CHECKEXP(tmpExp)) THEN
        ElectronicDistriPart(iPart1)%DistriFunc(iQua+1) = (1.-LocRelaxProb)*DistriOld(iQua+1) &
            + LocRelaxProb * SpecDSMC(iSpec)%ElectronicState(1,iQua) * EXP (-tmpExp)/ElectronicPartition
      ELSE
        ElectronicDistriPart(iPart1)%DistriFunc(iQua+1) = (1.-LocRelaxProb)*DistriOld(iQua+1)
      END IF
      PartStateIntEn(3,iPart1) = PartStateIntEn(3,iPart1) + &
          ElectronicDistriPart(iPart1)%DistriFunc(iQua+1) * BoltzmannConst * SpecDSMC(iSpec)%ElectronicState(2,iQua)
    END DO
    IF ((Coll_pData(iPair)%Ec-PartStateIntEn(3,iPart1)*GetParticleWeight(iPart1)).LT.0.0) CALL abort(__STAMP__,&
        'Negative collision energy after electronic excitation relaxation!')
  END IF
CASE(3)

  ! Set optional variable
  IF(PRESENT(XSec_Level))THEN
    IF(XSec_Level.GT.0) THEN
      iCase = CollInf%Coll_Case(PartSpecies(Coll_pData(iPair)%iPart_p1),PartSpecies(Coll_pData(iPair)%iPart_p2))
      ! Set energy according the selected level
      PartStateIntEn(3,iPart1) = SpecXSec(iCase)%ElecLevel(XSec_Level)%Threshold
    ELSE
      ! Reset electronic energy to ground state after a chemical reaction (XSec_Level = 0)
      PartStateIntEn(3,iPart1) = 0.0
    END IF
  ELSE
    CALL abort(__STAMP__,'XSec_Level argument not given to ElectronicEnergyExchange()')
  END IF ! PRESENT(XSec_Level)

CASE DEFAULT
  CALL abort(__STAMP__,'ERROR: Unknown electronic relaxation model: ',IntInfoOpt=DSMC%ElectronicModel)
END SELECT

END SUBROUTINE ElectronicEnergyExchange


SUBROUTINE LT_ElectronicEnergyExchange(iPartIndx_Node, nPart, NodeVolume)
!===================================================================================================================================
!> Subroutine for the cell-local Landau-Teller based electronic excitation relaxation using BGK approach:
!> 1.) Moment calculation: Summing up the relative velocities and their squares
!> 2.) Calculation of the relaxation frequency of the distribution function towards the target electronic temperature
!> 3.) determination of the electronic relaxation frequency
!> 4.) Determine the number of particles undergoing a relaxation
!> 5.) Determine the new electronic state of molecules undergoing a relaxation
!> 6.) Treatment of the electronic energy of molecules
!> 7.) Determine the new PartState
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars           ,ONLY: PartState, Species, PartSpecies, nSpecies, usevMPF, UseVarTimeStep
USE MOD_DSMC_Vars               ,ONLY: SpecDSMC, PartStateIntEn, RadialWeighting, CollInf, ElecRelaxPart
USE MOD_TimeDisc_Vars           ,ONLY: dt
USE MOD_part_tools              ,ONLY: GetParticleWeight, CalcXiElec
USE MOD_Globals_Vars            ,ONLY: BoltzmannConst
USE MOD_Particle_Analyze_Tools  ,ONLY: CalcTelec
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                        :: NodeVolume
INTEGER, INTENT(INOUT)                  :: nPart
INTEGER, INTENT(INOUT)                  :: iPartIndx_Node(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: alpha, CellTemp, dens, NewEn, OldEn, TEqui, dtCell, NewEnElec(nSpecies), iRan, ElectronicPartitionTemp
INTEGER, ALLOCATABLE  :: iPartIndx_NodeRelaxElec(:)
INTEGER               :: iSpec, nSpec(nSpecies), jSpec, nElecRelax, iLoop, iPart, nElecRelaxSpec(nSpecies), iQua
REAL                  :: vBulkAll(3), SpecTemp(nSpecies), TElecSpec(nSpecies), ElecExpSpec(nSpecies)
REAL                  :: totalWeightSpec(nSpecies), totalWeight, partWeight, CellTemptmp
REAL                  :: EElecSpec(nSpecies), Xi_ElecSpec(nSpecies), Xi_Elec_oldSpec(nSpecies), EElecMean(nSpecies)
REAL                  :: collisionfreqSpec(nSpecies),elecrelaxfreqSpec(nSpecies), ElectronicPartition(nSpecies)
LOGICAL               :: SkipEnergyCons(nSpecies)
!===================================================================================================================================
IF(nPart.LT.2) RETURN

! 1.) Moment calculation: Summing up the relative velocities and their squares
CALL CalcMoments_ElectronicExchange(nPart, iPartIndx_Node, nSpec, vBulkAll, totalWeight, totalWeightSpec, &
                      OldEn, EElecSpec, CellTemp, SpecTemp, dtCell)
NewEn = OldEn
IF((CellTemp.LE.0).OR.(MAXVAL(nSpec(:)).EQ.1).OR.(totalWeight.LE.0.0)) RETURN

IF(UseVarTimeStep) THEN
  dtCell = dt * dtCell / totalWeight
ELSE
  dtCell = dt
END IF

IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
  ! totalWeight contains the weighted particle number
  dens = totalWeight / NodeVolume
ELSE
  dens = totalWeight * Species(1)%MacroParticleFactor / NodeVolume
END IF

! Calculation of the electronic degrees of freedom for molecules
Xi_ElecSpec=0.; Xi_Elec_oldSpec=0.; TElecSpec=0.
DO iSpec = 1, nSpecies
  IF (nSpec(iSpec).EQ.0) CYCLE
  IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
    TElecSpec(iSpec)=CalcTelec( EElecSpec(iSpec)/totalWeightSpec(iSpec), iSpec)
    Xi_ElecSpec(iSpec)=CalcXiElec(TElecSpec(iSpec),iSpec)
    Xi_Elec_oldSpec(iSpec) = Xi_ElecSpec(iSpec)
  END IF
END DO

! 3.) Determination of electronic relaxation frequency using the collision frequency
collisionfreqSpec = 0.0
DO iSpec = 1, nSpecies
  DO jSpec = 1, nSpecies
    IF (iSpec.EQ.jSpec) THEN
      CellTemptmp = CellTemp !SpecTemp(iSpec)
    ELSE
      CellTemptmp = CellTemp
    END IF
    collisionfreqSpec(iSpec) = collisionfreqSpec(iSpec) + SpecDSMC(iSpec)%CollFreqPreFactor(jSpec) * totalWeightSpec(jSpec) &
            * (dens / totalWeight) *CellTemptmp**(-CollInf%omega(iSpec,jSpec) +0.5)
  END DO
END DO
elecrelaxfreqSpec(:) = collisionfreqSpec(:) * SpecDSMC(:)%ElecRelaxProb
ElecExpSpec=0.
DO iSpec = 1, nSpecies
  IF (totalWeightSpec(iSpec).GT.0.0) THEN
    EElecMean(iSpec) =  EElecSpec(iSpec)/totalWeightSpec(iSpec)
  ELSE
    EElecMean(iSpec) = 0.0
  END IF
END DO
CALL CalcTEquiMultiElec(nPart, nSpec, CellTemp, TElecSpec, Xi_ElecSpec, Xi_Elec_oldSpec, ElecExpSpec,   &
      TEqui, elecrelaxfreqSpec, dtCell, EElecMean(:))

! 4.) Determine the number of particles undergoing a relaxation (including vibration and rotation)
ALLOCATE(iPartIndx_NodeRelaxElec(nPart))
iPartIndx_NodeRelaxElec = 0

nElecRelaxSpec =0; nElecRelax=0
DO iLoop = 1, nPart
  iPart = iPartIndx_Node(iLoop)
  iSpec = PartSpecies(iPart)
  IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
    IF (.NOT.ElecRelaxPart(iPart)) CYCLE
    partWeight = GetParticleWeight(iPart)
    CALL RANDOM_NUMBER(iRan)
    IF ((1.-ElecExpSpec(iSpec)).GT.iRan) THEN
      nElecRelax = nElecRelax + 1
      nElecRelaxSpec(iSpec) = nElecRelaxSpec(iSpec) + 1
      iPartIndx_NodeRelaxElec(nElecRelax) = iPart
      OldEn = OldEn + PartStateIntEn(3,iPartIndx_NodeRelaxElec(nElecRelax))* partWeight
    END IF
  END IF
END DO
IF ((nElecRelax.EQ.0)) RETURN

ElectronicPartition = 0.
DO iSpec =1, nSpecies
  IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
    ElectronicPartitionTemp = 0.
    ! calculate sum over all energy levels == partition function for temperature Telec
    DO iQua = 0, SpecDSMC(iSpec)%MaxElecQuant - 1
      ElectronicPartitionTemp = SpecDSMC(iSpec)%ElectronicState(1,iQua) * EXP(-SpecDSMC(iSpec)%ElectronicState(2,iQua)/TEqui)
      IF ( ElectronicPartitionTemp .GT. ElectronicPartition(iSpec) ) THEN
        ElectronicPartition(iSpec) = ElectronicPartitionTemp
      END IF
    END DO
  END IF
END DO

! 5.) Determine the new electronical state of molecules undergoing a relaxation
NewEnElec = 0.0
SkipEnergyCons = .FALSE.
DO iLoop = 1, nElecRelax
  iPart = iPartIndx_NodeRelaxElec(iLoop)
  iSpec = PartSpecies(iPart)
  partWeight = GetParticleWeight(iPart)
  CALL RANDOM_NUMBER(iRan)
  iQua = INT( ( SpecDSMC(iSpec)%MaxElecQuant ) * iRan)
  ElectronicPartitionTemp = SpecDSMC(iSpec)%ElectronicState(1,iQua) * EXP(-SpecDSMC(iSpec)%ElectronicState(2,iQua)/TEqui)
  ! select level
  CALL RANDOM_NUMBER(iRan)
  DO WHILE ( iRan .GE. ElectronicPartitionTemp / ElectronicPartition(iSpec) )
    CALL RANDOM_NUMBER(iRan)
    iQua = INT( ( SpecDSMC(iSpec)%MaxElecQuant ) * iRan)
    ElectronicPartitionTemp = SpecDSMC(iSpec)%ElectronicState(1,iQua) * EXP(-SpecDSMC(iSpec)%ElectronicState(2,iQua)/TEqui)
    CALL RANDOM_NUMBER(iRan)
  END DO
  iQua = MAX(iQua,1)
  PartStateIntEn( 3,iPart) = BoltzmannConst*SpecDSMC(iSpec)%ElectronicState(2,iQua)
  IF (TEqui.GT.SpecDSMC(iSpec)%MaxMeanXiElec(2)) THEN
    SkipEnergyCons(iSpec)= .TRUE.
    IF ((OldEn - (PartStateIntEn(3,iPart)*partWeight)).LT.0.0) THEN
      DO WHILE ((OldEn - (PartStateIntEn(3,iPart)*partWeight)).LT.0.0)
        iQua = iQua - 1
        PartStateIntEn(3,iPart) = BoltzmannConst*SpecDSMC(iSpec)%ElectronicState(2,iQua)
        IF (iQua.EQ.0) EXIT
      END DO
    END IF
    OldEn = OldEn - PartStateIntEn(3,iPart) * partWeight
  ELSE
    NewEnElec(iSpec) = NewEnElec(iSpec) + PartStateIntEn(3,iPart) * partWeight
  END IF
END DO
! 6.) Vibrational energy of the molecules: Ensure energy conservation by scaling the new electronic states with the factor alpha
CALL EnergyConsElec(nPart, nElecRelax, nElecRelaxSpec, iPartIndx_NodeRelaxElec, NewEnElec, OldEn, Xi_ElecSpec, SkipEnergyCons)

! 8.) Determine the new particle state and ensure energy conservation by scaling the new velocities with the factor alpha.
alpha = SQRT(OldEn/NewEn)
DO iLoop = 1, nPart
  iPart = iPartIndx_Node(iLoop)
  PartState(4:6,iPart) = vBulkAll(1:3) + alpha*(PartState(4:6,iPart)-vBulkAll(1:3))
END DO

END SUBROUTINE LT_ElectronicEnergyExchange


SUBROUTINE LT_ElectronicEnergyExchangeChem(iPartIndx_Node, nPart)
!===================================================================================================================================
!> Subroutine for the electronic excitation of chemical reacting particles, same steps as LT_ElectronicEnergyExchange but only for
!> chemical reacting particles assuming an thermal equilibrium between translation and electronical excitation after the chemical
!> reaction
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars         ,ONLY: PartState, Species, PartSpecies, nSpecies
USE MOD_DSMC_Vars             ,ONLY: SpecDSMC, DSMC, PartStateIntEn
USE MOD_part_tools            ,ONLY: GetParticleWeight, CalcXiElec
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(INOUT)                  :: nPart
INTEGER, INTENT(INOUT)                  :: iPartIndx_Node(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: alpha, NewEn, OldEn, TEqui, NewEnElec(nSpecies), ElectronicPartitionTemp, iRan
INTEGER, ALLOCATABLE  :: iPartIndx_NodeRelaxElec(:)
INTEGER               :: iSpec, nSpec(nSpecies), nElecRelax, iLoop, iPart, nElecRelaxSpec(nSpecies), iQua
REAL                  :: vBulkAll(3), totalWeightSpec(nSpecies), totalWeight, partWeight, Xi_ElecSpec(nSpecies)
REAL                  :: MaxTemp, MinTemp, TempEn, Xi_ElecTotal, V_rel(3), TotalMass, vmag2, ElectronicPartition(nSpecies)
LOGICAL               :: SkipEnergyCons(nSPecies)
!===================================================================================================================================
IF(nPart.LT.2) RETURN

! Calculate necessary moments
totalWeightSpec = 0.0; vBulkAll=0.0; TotalMass=0.0; nSpec=0
DO iLoop = 1, nPart
  iPart = iPartIndx_Node(iLoop)
  partWeight = GetParticleWeight(iPart)
  iSpec = PartSpecies(iPart)
  totalWeightSpec(iSpec) = totalWeightSpec(iSpec) + partWeight
  vBulkAll(1:3)  =  vBulkAll(1:3) + PartState(4:6,iPart)*Species(iSpec)%MassIC*partWeight
  TotalMass = TotalMass + Species(iSpec)%MassIC*partWeight
  nSpec(iSpec) = nSpec(iSpec) + 1
END DO
totalWeight = SUM(totalWeightSpec)
vBulkAll(1:3) = vBulkAll(1:3) / TotalMass

OldEn=0.0
DO iLoop = 1, nPart
  iPart = iPartIndx_Node(iLoop)
  partWeight = GetParticleWeight(iPart)
  iSpec = PartSpecies(iPart)
  V_rel(1:3)=PartState(4:6,iPart)-vBulkAll(1:3)
  vmag2 = V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2
  OldEn = OldEn + 0.5*Species(iSpec)%MassIC * vmag2*partWeight
END DO

NewEn = OldEn
DO iLoop = 1, nPart
  iPart = iPartIndx_Node(iLoop)
  iSpec = PartSpecies(iPart)
  IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
    partWeight = GetParticleWeight(iPart)
    OldEn = OldEn + PartStateIntEn(3,iPart)* partWeight
  END IF
END DO

! Calculation of equilibrium temperature by using DOF's with a bisection
Xi_ElecSpec=0.
MaxTemp=5.*DSMC%InstantTransTemp(nSpecies + 1)
MinTemp=1E-6
TempEn=0.0
iloop = 0
DO WHILE(.NOT.ALMOSTEQUAL(TempEn, OldEn))
  iloop = iloop + 1
  IF (iLoop.EQ.100) THEN
    DO iSpec = 1, nSpecies
      IF (nSpec(iSpec).EQ.0) THEN
        Xi_ElecSpec(iSpec)= 0.0
        CYCLE
      END IF
      IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
        Xi_ElecSpec(iSpec)=CalcXiElec(DSMC%InstantTransTemp(nSpecies + 1),iSpec)
        TEqui = DSMC%InstantTransTemp(nSpecies + 1)
      END IF
    END DO
    RETURN
  END IF
  Tequi= 0.5*(MaxTemp+MinTemp)
  DO iSpec = 1, nSpecies
    IF (nSpec(iSpec).EQ.0) THEN
      Xi_ElecSpec(iSpec)= 0.0
      CYCLE
    END IF
    IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
      Xi_ElecSpec(iSpec)=CalcXiElec(Tequi,iSpec)
    END IF
  END DO
  Xi_ElecTotal = 0.0
  DO iSpec = 1, nSpecies
    Xi_ElecTotal = Xi_ElecTotal + Xi_ElecSpec(iSpec)*totalWeightSpec(iSpec)
  END DO

  TempEn = (3.*(nPart-1.)/nPart*totalWeight+Xi_ElecTotal)/2.*BoltzmannConst*TEqui
  IF (TempEn.GT.OldEn) THEN
   MaxTemp = TEqui
  ELSE
   MinTemp = TEqui
  END IF
END DO

ElectronicPartition = 0.
DO iSpec =1, nSpecies
  IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
    ElectronicPartitionTemp = 0.
    ! calculate sum over all energy levels == partition function for temperature Telec
    DO iQua = 0, SpecDSMC(iSpec)%MaxElecQuant - 1
      ElectronicPartitionTemp = SpecDSMC(iSpec)%ElectronicState(1,iQua) * EXP(-SpecDSMC(iSpec)%ElectronicState(2,iQua)/TEqui)
      IF ( ElectronicPartitionTemp .GT. ElectronicPartition(iSpec) ) THEN
        ElectronicPartition(iSpec) = ElectronicPartitionTemp
      END IF
    END DO
  END IF
END DO

ALLOCATE(iPartIndx_NodeRelaxElec(nPart))
iPartIndx_NodeRelaxElec = 0
! Determine the new electronic state of molecules undergoing a relaxation
NewEnElec = 0.0
nElecRelax = 0
SkipEnergyCons = .FALSE.
DO iLoop = 1, nPart
  iPart = iPartIndx_Node(iLoop)
  iSpec = PartSpecies(iPart)
  IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
    nElecRelax  = nElecRelax + 1
    nElecRelaxSpec(iSpec) = nElecRelaxSpec(iSpec) + 1
    iPartIndx_NodeRelaxElec(nElecRelax) = iPart
    partWeight = GetParticleWeight(iPart)
    CALL RANDOM_NUMBER(iRan)
    iQua = INT( ( SpecDSMC(iSpec)%MaxElecQuant ) * iRan)
    ElectronicPartitionTemp = SpecDSMC(iSpec)%ElectronicState(1,iQua) * EXP(-SpecDSMC(iSpec)%ElectronicState(2,iQua)/TEqui)
    ! select level
    CALL RANDOM_NUMBER(iRan)
    DO WHILE ( iRan .GE. ElectronicPartitionTemp / ElectronicPartition(iSpec) )
      CALL RANDOM_NUMBER(iRan)
      iQua = INT( ( SpecDSMC(iSpec)%MaxElecQuant ) * iRan)
      ElectronicPartitionTemp = SpecDSMC(iSpec)%ElectronicState(1,iQua) * EXP(-SpecDSMC(iSpec)%ElectronicState(2,iQua)/TEqui)
      CALL RANDOM_NUMBER(iRan)
    END DO
    iQua = MAX(iQua,1)
    PartStateIntEn( 3,iPart) = BoltzmannConst*SpecDSMC(iSpec)%ElectronicState(2,iQua)
    IF (TEqui.GT.SpecDSMC(iSpec)%MaxMeanXiElec(2)) THEN
      SkipEnergyCons(iSpec)= .TRUE.
      IF ((OldEn - (PartStateIntEn(3,iPart)*partWeight)).LT.0.0) THEN
        DO WHILE ((OldEn - (PartStateIntEn(3,iPart)*partWeight)).LT.0.0)
          iQua = iQua - 1
          PartStateIntEn(3,iPart) = BoltzmannConst*SpecDSMC(iSpec)%ElectronicState(2,iQua)
          IF (iQua.EQ.0) EXIT
        END DO
      END IF
      OldEn = OldEn - PartStateIntEn(3,iPart) * partWeight
    ELSE
      NewEnElec(iSpec) = NewEnElec(iSpec) + PartStateIntEn(3,iPart) * partWeight
    END IF
  END IF
END DO

! 7.) Electronical energy of the molecules: Ensure energy conservation by scaling the new electroncal states with the factor alpha
CALL EnergyConsElec(nPart, nElecRelax, nElecRelaxSpec, iPartIndx_NodeRelaxElec, NewEnElec, OldEn, Xi_ElecSpec, SkipEnergyCons)
! 8.) Determine the new particle state and ensure energy conservation by scaling the new velocities with the factor alpha.

alpha = SQRT(OldEn/NewEn)
DO iLoop = 1, nPart
  iPart = iPartIndx_Node(iLoop)
  PartState(4:6,iPart) = vBulkAll(1:3) + alpha*(PartState(4:6,iPart)-vBulkAll(1:3))
END DO

END SUBROUTINE LT_ElectronicEnergyExchangeChem


SUBROUTINE LT_ElectronicExc_ConstructPartList(iPartIndx_NodeTotal, iPartIndx_NodeTotalElecExc,  nPart, nPartRelax)
!===================================================================================================================================
!> Construct List of Particles that must be checked for electronic excitation.
!> Sort out the particles that have already been relaxed by a chemical reaction.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PARTICLE_Vars       ,ONLY: PDM
USE MOD_DSMC_Vars           ,ONLY: newElecRelaxParts, iPartIndx_NodeNewElecRelax
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                  :: iPartIndx_NodeTotal(:),nPart
INTEGER, INTENT(OUT),ALLOCATABLE    :: iPartIndx_NodeTotalElecExc(:)
INTEGER, INTENT(OUT)                :: nPartRelax
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iPart, iLoop
!===================================================================================================================================
nPartRelax = 0
! 1) Treating the particles that already existed within the cell (but may have changed their species)
DO iLoop = 1, nPart
  iPart = iPartIndx_NodeTotal(iLoop)
  IF (PDM%ParticleInside(iPart)) THEN
    nPartRelax = nPartRelax + 1
  END IF
END DO
! 2) Treating particles that were newly created during the time step
DO iLoop = 1, newElecRelaxParts
  iPart = iPartIndx_NodeNewElecRelax(iLoop)
  IF (PDM%ParticleInside(iPart)) THEN
    nPartRelax = nPartRelax + 1
  END IF
END DO

ALLOCATE(iPartIndx_NodeTotalElecExc(nPartRelax))
nPartRelax = 0
! 1) Treating the particles that already existed within the cell (but may have changed their species)
DO iLoop = 1, nPart
  iPart = iPartIndx_NodeTotal(iLoop)
  IF (PDM%ParticleInside(iPart)) THEN
    nPartRelax = nPartRelax + 1
    iPartIndx_NodeTotalElecExc(nPartRelax) = iPart
  END IF
END DO
! 2) Treating particles that were newly created during the time step
DO iLoop = 1, newElecRelaxParts
  iPart = iPartIndx_NodeNewElecRelax(iLoop)
  IF (PDM%ParticleInside(iPart)) THEN
    nPartRelax = nPartRelax + 1
    iPartIndx_NodeTotalElecExc(nPartRelax) = iPart
  END IF
END DO

END SUBROUTINE LT_ElectronicExc_ConstructPartList

SUBROUTINE CalcMoments_ElectronicExchange(nPart, iPartIndx_Node, nSpec, vBulkAll, totalWeight, totalWeightSpec, &
                      OldEn, EElecSpec, CellTemp, SpecTemp, dtCell)
!===================================================================================================================================
!> Moment calculation: Summing up the relative velocities and their squares
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars         ,ONLY: PartState, Species, PartSpecies, nSpecies, UseVarTimeStep, PartTimeStep
USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, SpecDSMC
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nPart, iPartIndx_Node(nPart)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(OUT)          :: nSpec(nSpecies)
REAL, INTENT(OUT)             :: OldEn, EElecSpec(nSpecies)
REAL, INTENT(OUT)             :: CellTemp, SpecTemp(nSpecies), totalWeightSpec(nSpecies)
REAL, INTENT(OUT)             :: vBulkAll(3), totalWeight, dtCell
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iLoop, iPart, iSpec
REAL                          :: V_rel(1:3), vmag2, partWeight, EnerTotal, totalWeightSpec2(nSpecies), vBulkSpec(3,nSpecies),u2
REAL                          :: tempweight, tempweight2, tempmass, vBulkTemp(3), totalWeight2, u2Spec(nSpecies), TotalMass
!===================================================================================================================================

totalWeightSpec = 0.0; totalWeightSpec2=0.0; vBulkAll=0.0; TotalMass=0.0; vBulkSpec=0.0; nSpec=0; dtCell=0.0
DO iLoop = 1, nPart
  iPart = iPartIndx_Node(iLoop)
  partWeight = GetParticleWeight(iPart)
  iSpec = PartSpecies(iPart)
  totalWeightSpec(iSpec) = totalWeightSpec(iSpec) + partWeight
  totalWeightSpec2(iSpec) =   totalWeightSpec2(iSpec) + partWeight*partWeight
  vBulkAll(1:3)  =  vBulkAll(1:3) + PartState(4:6,iPart)*Species(iSpec)%MassIC*partWeight
  TotalMass = TotalMass + Species(iSpec)%MassIC*partWeight
  vBulkSpec(1:3,iSpec) = vBulkSpec(1:3,iSpec) + PartState(4:6,iPart)*partWeight
  nSpec(iSpec) = nSpec(iSpec) + 1
  IF(UseVarTimeStep) THEN
    dtCell = dtCell + PartTimeStep(iPart)*partWeight
  END IF
END DO
totalWeight = SUM(totalWeightSpec)
totalWeight2 = SUM(totalWeightSpec2)
IF ((MAXVAL(nSpec(:)).EQ.1).OR.(totalWeight.LE.0.0)) RETURN
vBulkAll(1:3) = vBulkAll(1:3) / TotalMass
DO iSpec = 1, nSpecies
  IF (nSpec(iSpec).GT.0) vBulkSpec(:,iSpec) = vBulkSpec(:,iSpec) /totalWeightSpec(iSpec)
END DO

u2Spec=0.0; OldEn=0.0; EElecSpec=0.0
DO iLoop = 1, nPart
  iPart = iPartIndx_Node(iLoop)
  partWeight = GetParticleWeight(iPart)
  iSpec = PartSpecies(iPart)
  V_rel(1:3)=PartState(4:6,iPart)-vBulkSpec(1:3,iSpec)
  vmag2 = V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2
  u2Spec(iSpec) = u2Spec(iSpec) + vmag2*partWeight

  V_rel(1:3)=PartState(4:6,iPart)-vBulkAll(1:3)
  vmag2 = V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2
  OldEn = OldEn + 0.5*Species(iSpec)%MassIC * vmag2*partWeight
  IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
    EElecSpec(iSpec) = EElecSpec(iSpec) + PartStateIntEn(3,iPart) * partWeight
  END IF
END DO

IF (nSpecies.GT.1) THEN
  SpecTemp = 0.0
  EnerTotal = 0.0
  tempweight = 0.0; tempweight2 = 0.0; tempmass = 0.0; vBulkTemp = 0.0
  DO iSpec = 1, nSpecies
    IF ((nSpec(iSpec).GE.2).AND.(.NOT.ALMOSTZERO(u2Spec(iSpec)))) THEN
      SpecTemp(iSpec) = Species(iSpec)%MassIC * u2Spec(iSpec) &
          /(3.0*BoltzmannConst*(totalWeightSpec(iSpec) - totalWeightSpec2(iSpec)/totalWeightSpec(iSpec)))
      EnerTotal =  EnerTotal + 3./2.*BoltzmannConst*SpecTemp(iSpec) * totalWeightSpec(iSpec)
      vmag2 = vBulkSpec(1,iSpec)**(2.) + vBulkSpec(2,iSpec)**(2.) + vBulkSpec(3,iSpec)**(2.)
      EnerTotal = EnerTotal + totalWeightSpec(iSpec) * Species(iSpec)%MassIC / 2. * vmag2
      tempweight = tempweight + totalWeightSpec(iSpec)
      tempweight2 = tempweight2 + totalWeightSpec2(iSpec)
      tempmass = tempmass +  totalWeightSpec(iSpec) * Species(iSpec)%MassIC
      vBulkTemp(1:3) = vBulkTemp(1:3) + vBulkSpec(1:3,iSpec)*totalWeightSpec(iSpec) * Species(iSpec)%MassIC
    END IF
  END DO

  vBulkTemp(1:3) = vBulkTemp(1:3) / tempmass
  vmag2 = vBulkTemp(1)*vBulkTemp(1) + vBulkTemp(2)*vBulkTemp(2) + vBulkTemp(3)*vBulkTemp(3)
  EnerTotal = EnerTotal -  tempmass / 2. * vmag2
  CellTemp = 2. * EnerTotal / (3.*tempweight*BoltzmannConst)
ELSE
  u2 = u2Spec(1) / (totalWeight - totalWeight2/totalWeight)
  CellTemp = Species(1)%MassIC * u2 / (3.0*BoltzmannConst)
END IF

END SUBROUTINE CalcMoments_ElectronicExchange


SUBROUTINE CalcTEquiMultiElec(nPart, nSpec, CellTemp, TElecSpec, Xi_ElecSpec, Xi_Elec_oldSpec, ElecExpSpec,   &
      TEqui, elecrelaxfreqSpec, dtCell, meanEelecSpec)
!===================================================================================================================================
! Calculation of the equilibrium temperature for electronic excitation
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars,           ONLY: BoltzmannConst
USE MOD_DSMC_Vars,              ONLY: SpecDSMC
USE MOD_Particle_Vars,          ONLY: nSpecies
USE MOD_part_tools,             ONLY: CalcXiElec
USE MOD_Particle_Analyze_Tools, ONLY: CalcEelec
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: CellTemp, TElecSpec(nSpecies), Xi_Elec_oldSpec(nSpecies), meanEelecSpec(nSpecies)
REAL, INTENT(IN)                :: elecrelaxfreqSpec(nSpecies), dtCell
INTEGER, INTENT(IN)             :: nPart, nSpec(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)               :: Xi_ElecSpec(nSpecies), TEqui, ElecExpSpec(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL                            :: TEqui_Old, betaElec, ElecFracSpec(nSpecies), TEqui_Old2, EelecTtrans(nSpecies), EElecTequi
REAL                            :: eps_prec=1.0E-1
REAL                            :: correctFac,  maxexp, TEquiNumDof
INTEGER                         :: iSpec
!===================================================================================================================================
maxexp = LOG(HUGE(maxexp))
!  Xi_rel = 2.*(2. - CollInf%omega(1,1))
!  correctFac = 1. + (2.*SpecDSMC(1)%CharaTVib / (CellTemp*(EXP(SpecDSMC(1)%CharaTVib / CellTemp)-1.)))**(2.) &
!        * EXP(SpecDSMC(1)%CharaTVib /CellTemp) / (2.*Xi_rel)

correctFac = 1.
ElecFracSpec = 0.0
DO iSpec=1, nSpecies
  IF ((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
    ElecExpSpec(iSpec) = exp(-elecrelaxfreqSpec(iSpec)*dtCell/correctFac)
    ElecFracSpec(iSpec) = nSpec(iSpec)*(1.-ElecExpSpec(iSpec))
    EelecTtrans(iSpec) = CalcEelec(CellTemp, iSpec)
  END IF
END DO
TEqui_Old = 0.0
TEqui = 3.*(nPart-1.)*CellTemp
TEquiNumDof = 3.*(nPart-1.)
DO iSpec=1, nSpecies
  IF ((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
    TEqui = TEqui + Xi_Elec_oldSpec(iSpec)*ElecFracSpec(iSpec)*TElecSpec(iSpec)
    TEquiNumDof = TEquiNumDof + Xi_Elec_oldSpec(iSpec)*ElecFracSpec(iSpec)
  END IF
END DO
TEqui = TEqui / TEquiNumDof
DO WHILE ( ABS( TEqui - TEqui_Old ) .GT. eps_prec )
  DO iSpec = 1, nSpecies
    IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
      Xi_ElecSpec(iSpec) = CalcXiElec(TEqui,iSpec)
      EElecTequi= CalcEelec(TEqui, iSpec)
      IF (ABS(meanEelecSpec(iSpec)-EElecTequi).LT.1E-3) THEN
        ElecExpSpec(iSpec) = exp(-elecrelaxfreqSpec(iSpec)*dtCell/correctFac)
      ELSE
        betaElec = ((meanEelecSpec(iSpec)-EelecTtrans(iSpec))/(meanEelecSpec(iSpec)-EElecTequi))*elecrelaxfreqSpec(iSpec)*dtCell/correctFac
        IF (-betaElec.GT.0.0) THEN
          ElecExpSpec(iSpec) = 0.
        ELSE IF (betaElec.GT.maxexp) THEN
          ElecExpSpec(iSpec) = 0.
        ELSE
          ElecExpSpec(iSpec) = exp(-betaElec)
        END IF
      END IF
      ElecFracSpec(iSpec) = nSpec(iSpec)*(1.-ElecExpSpec(iSpec))
    END IF
  END DO
  TEqui_Old = TEqui
  TEqui_Old2 = TEqui

  TEqui = 3.*(nPart-1.)*CellTemp
  TEquiNumDof = 3.*(nPart-1.)
  DO iSpec=1, nSpecies
    IF ((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
      TEqui = TEqui + Xi_Elec_oldSpec(iSpec)*ElecFracSpec(iSpec)*TElecSpec(iSpec)
      TEquiNumDof = TEquiNumDof + Xi_ElecSpec(iSpec)*ElecFracSpec(iSpec)
    END IF
  END DO
  TEqui = TEqui / TEquiNumDof
END DO

END SUBROUTINE CalcTEquiMultiElec


SUBROUTINE EnergyConsElec(nPart, nElecRelax, nElecRelaxSpec, iPartIndx_NodeRelaxElec, NewEnElec, OldEn, Xi_ElecSpec, SkipEnergyCons)
!===================================================================================================================================
!> Routine to ensure energy conservation electroncal degrees of freedom
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars         ,ONLY: PartSpecies, nSpecies
USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, SpecDSMC
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nPart, nElecRelax, iPartIndx_NodeRelaxElec(:), nElecRelaxSpec(nSpecies)
LOGICAL, INTENT(IN)           :: SkipEnergyCons(nSpecies)
REAL, INTENT(IN)              :: NewEnElec(nSpecies), Xi_ElecSpec(nSpecies)
REAL, INTENT(INOUT)           :: OldEn
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart, iLoop, iSpec, iQuant, iQuaMax
REAL                          :: alpha(nSpecies), partWeight, betaV, iRan, Xi_ElecTotal, Prob
!===================================================================================================================================
IF (ANY(NewEnElec.GT.0.0)) THEN
  Xi_ElecTotal = 0.0
  DO iSpec = 1, nSpecies
    IF (SkipEnergyCons(iSpec)) CYCLE
    Xi_ElecTotal = Xi_ElecTotal + Xi_ElecSpec(iSpec)*nElecRelaxSpec(iSpec)
  END DO

  DO iSpec = 1, nSpecies
    IF (SkipEnergyCons(iSpec)) CYCLE
    IF (NewEnElec(iSpec).GT.0.0) THEN
      alpha(iSpec) = OldEn/NewEnElec(iSpec)*(Xi_ElecSpec(iSpec)*nElecRelaxSpec(iSpec)/(3.*(nPart-1.)+Xi_ElecTotal))
    ELSE
      alpha(iSpec) = 0.0
    END IF
  END DO
  DO iLoop = 1, nElecRelax
    iPart = iPartIndx_NodeRelaxElec(iLoop)
    partWeight = GetParticleWeight(iPart)
    iSpec = PartSpecies(iPart)
    IF (SkipEnergyCons(iSpec)) CYCLE
    betaV = alpha(iSpec)*PartStateIntEn( 3,iPart)
    DO iQuant = 1,  SpecDSMC(iSpec)%MaxElecQuant - 1
      IF (betaV.LT.BoltzmannConst * SpecDSMC(iSpec)%ElectronicState(2,iQuant)) THEN
        iQuaMax = iQuant
        EXIT
      END IF
      IF(iQuant.EQ.SpecDSMC(iSpec)%MaxElecQuant - 1) iQuaMax = iQuant
    END DO
    Prob = (betaV-BoltzmannConst*SpecDSMC(iSpec)%ElectronicState(2,iQuaMax-1)) &
        / (BoltzmannConst*(SpecDSMC(iSpec)%ElectronicState(2,iQuaMax) - SpecDSMC(iSpec)%ElectronicState(2,iQuaMax-1)))
    CALL RANDOM_NUMBER(iRan)
    IF (iRan.GT.Prob) THEN
      iQuant = iQuaMax -1
      PartStateIntEn(3,iPart) = BoltzmannConst*SpecDSMC(iSpec)%ElectronicState(2,iQuaMax-1)
    ELSE
      iQuant = iQuaMax
      PartStateIntEn(3,iPart) = BoltzmannConst*SpecDSMC(iSpec)%ElectronicState(2,iQuaMax)
    END IF
    IF ((OldEn - (PartStateIntEn(3,iPart)*partWeight)).LT.0.0) THEN
      DO WHILE ((OldEn - (PartStateIntEn(3,iPart)*partWeight)).LT.0.0)
        iQuant = iQuant - 1
        PartStateIntEn(3,iPart) = BoltzmannConst*SpecDSMC(iSpec)%ElectronicState(2,iQuant)
        IF (iQuant.EQ.0) EXIT
      END DO
    END IF
    OldEn = OldEn - PartStateIntEn(3,iPart)*partWeight
  END DO
END IF

END SUBROUTINE EnergyConsElec


SUBROUTINE TVEEnergyExchange(CollisionEnergy,iPart1,FakXi)
!===================================================================================================================================
! Translational-vibrational-electronic energy exchange
!===================================================================================================================================
USE MOD_DSMC_Vars             ,ONLY: DSMC, SpecDSMC, PartStateIntEn
USE MOD_Particle_Vars         ,ONLY: PartSpecies
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPart1
REAL, INTENT(IN)              :: FakXi
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(INOUT)           :: CollisionEnergy
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iQuaMax, MaxElecQuant, iQua
INTEGER                       :: jQVib, QMaxVib
REAL                          :: iRan, iRan2, gmax, gtemp, PartStateTemp, iRanVib
!===================================================================================================================================
iQuaMax  = 0
! Determine max electronic quant
MaxElecQuant = SpecDSMC(PartSpecies(iPart1))%MaxElecQuant - 1
! determine maximal Quant and term according to Eq (7) of Liechty
gmax = 0
PartStateTemp = CollisionEnergy / BoltzmannConst
DO iQua = 0, MaxElecQuant
  IF ( (PartStateTemp  &
            - SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua) &
            - DSMC%GammaQuant * SpecDSMC(PartSpecies(iPart1))%CharaTVib) &
      .ge. 0 ) THEN
    gtemp = SpecDSMC(PartSpecies(iPart1))%ElectronicState(1,iQua) * &
            ( CollisionEnergy - BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua) &
            -DSMC%GammaQuant * SpecDSMC(PartSpecies(iPart1))%CharaTVib * BoltzmannConst)**FakXi
    ! maximal possible Quant before term goes negative
    iQuaMax = iQua
    IF ( gtemp .gt. gmax ) THEN
    ! Quant of largest value of Eq (7)
      gmax = gtemp
    END IF
  END IF
END DO
! max iQuant for dicing
QMaxVib = INT(CollisionEnergy/(BoltzmannConst*SpecDSMC(PartSpecies(iPart1))%CharaTVib)  &
            - DSMC%GammaQuant)
QMaxVib = MIN(QMaxVib + 1, SpecDSMC(PartSpecies(iPart1))%MaxVibQuant)
CALL RANDOM_NUMBER(iRan)
CALL RANDOM_NUMBER(iRanVib)
iQua = INT( ( iQuaMax +1 ) * iRan)
jQVib =  INT(iRanVib * QMaxVib)
gtemp =( CollisionEnergy - BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua) &
      -(DSMC%GammaQuant + jQVib) * SpecDSMC(PartSpecies(iPart1))%CharaTVib * BoltzmannConst)
IF (gtemp.LE.0.0) THEN
  gtemp = 0.0
ELSE
  gtemp = SpecDSMC(PartSpecies(iPart1))%ElectronicState(1,iQua) *(gtemp)**FakXi
END IF
CALL RANDOM_NUMBER(iRan2)
! acceptance-rejection for iQuaElec
DO WHILE ( iRan2 .ge. gtemp / gmax )
  CALL RANDOM_NUMBER(iRan)
  CALL RANDOM_NUMBER(iRanVib)
  iQua = int( ( iQuaMax +1 ) * iRan)
  jQVib =  INT(iRanVib * QMaxVib)
  gtemp =( CollisionEnergy - BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua) &
        -(DSMC%GammaQuant + jQVib) * SpecDSMC(PartSpecies(iPart1))%CharaTVib * BoltzmannConst)
  IF (gtemp.LE.0.0) THEN
    gtemp = 0.0
  ELSE
    gtemp = SpecDSMC(PartSpecies(iPart1))%ElectronicState(1,iQua) *(gtemp)**FakXi
  END IF
  CALL RANDOM_NUMBER(iRan2)
END DO

PartStateIntEn(3,iPart1) = BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua)
PartStateIntEn(1,iPart1) = (jQVib + DSMC%GammaQuant) * BoltzmannConst &
              * SpecDSMC(PartSpecies(iPart1))%CharaTVib

END SUBROUTINE TVEEnergyExchange


SUBROUTINE ReadSpeciesLevel ( Dsetname, iSpec )
!===================================================================================================================================
! Subroutine to read the electronic levels from DSMCSpeciesElectronicState.h5
!===================================================================================================================================
! use module
USE MOD_io_hdf5
USE MOD_Globals
USE MOD_DSMC_Vars        ,ONLY: DSMC, SpecDSMC
USE MOD_HDF5_Input       ,ONLY: DatasetExists
USE MOD_part_tools       ,ONLY: CalcXiElec
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                    :: iSpec
CHARACTER(LEN=64),INTENT(IN)                          :: dsetname
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                               :: err
! HDF5 specifier taken from extractParticles
INTEGER(HSIZE_T), DIMENSION(2)                        :: dims,sizeMax
INTEGER(HID_T)                                        :: file_id_dsmc                       ! File identifier
INTEGER(HID_T)                                        :: dset_id_dsmc                       ! Dataset identifier
INTEGER(HID_T)                                        :: filespace                          ! filespace identifier
REAL,ALLOCATABLE                                      :: ElectronicState(:,:)
INTEGER, ALLOCATABLE                                  :: SortElectronicState(:)
INTEGER                                               :: iQua, nQuants, iQuaTemp
REAL                                                  :: tempEnergyDiff, tempEnergy, MaxTemp, MeanTemp, MeanMax, MinTemp
REAL                                                  :: XiMean, XiMeanMax
LOGICAL                                               :: DataSetFound
!===================================================================================================================================
LBWRITE(UNIT_StdOut,'(A)') 'Read electronic level entries '//TRIM(dsetname)//' from '//TRIM(DSMC%ElectronicModelDatabase)
! Initialize FORTRAN interface.
CALL H5OPEN_F(err)
! Open the file.
CALL H5FOPEN_F (TRIM(DSMC%ElectronicModelDatabase), H5F_ACC_RDONLY_F, file_id_dsmc, err)
CALL DatasetExists(File_ID_DSMC,TRIM(dsetname),DataSetFound)
IF(.NOT.DataSetFound)THEN
  CALL abort(&
  __STAMP__&
  ,'DataSet not found: ['//TRIM(dsetname)//'] ['//TRIM(DSMC%ElectronicModelDatabase)//']')
END IF
! Open the  dataset.
CALL H5DOPEN_F(file_id_dsmc, dsetname, dset_id_dsmc, err)
! Get the file space of the dataset.
CALL H5DGET_SPACE_F(dset_id_dsmc, FileSpace, err)
! get size
CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, dims, SizeMax, err)
! Allocate electronic_state
ALLOCATE (ElectronicState( 1:dims(1), 0:dims(2)-1 ) )
! read data
CALL H5dread_f(dset_id_dsmc, H5T_NATIVE_DOUBLE, ElectronicState, dims, err)
CALL SortEnergies(ElectronicState, INT(dims(2)))
IF (ALMOSTEQUAL(DSMC%EpsElecBin, 0.0).OR.(dims(2).EQ.2)) THEN
  ALLOCATE ( SpecDSMC(iSpec)%ElectronicState( 1:dims(1), 0:dims(2)-1 ) )
  SpecDSMC(iSpec)%ElectronicState = ElectronicState
  SpecDSMC(iSpec)%MaxElecQuant  = SIZE( SpecDSMC(iSpec)%ElectronicState,2)
ELSE
  ALLOCATE (SortElectronicState(0:dims(2)-1 ))
  SortElectronicState(0) = 0
  nQuants = 1
  tempEnergy =  ElectronicState(2,1)
  SortElectronicState(1) = nQuants
  DO iQua = 2, INT(dims(2),4)-2
    tempEnergyDiff = DiffElecEnergy(tempEnergy, ElectronicState(2,iQua))
    IF (tempEnergyDiff.LE.DSMC%EpsElecBin) THEN
      SortElectronicState(iQua) = nQuants
    ELSE
      nQuants = nQuants + 1
      SortElectronicState(iQua) = nQuants
      tempEnergy =  ElectronicState(2,iQua)
    END IF
  END DO
  nQuants = nQuants + 1
  SortElectronicState(dims(2)-1) = nQuants

  ALLOCATE ( SpecDSMC(iSpec)%ElectronicState( 1:dims(1), 0:nQuants) )
  SpecDSMC(iSpec)%ElectronicState = 0.0
  DO iQua = 1, INT(dims(2),4)-2
    iQuaTemp = SortElectronicState(iQua)
    SpecDSMC(iSpec)%ElectronicState( 1, iQuaTemp) = SpecDSMC(iSpec)%ElectronicState( 1, iQuaTemp) &
        + ElectronicState(1, iQua)
    SpecDSMC(iSpec)%ElectronicState( 2, iQuaTemp) = SpecDSMC(iSpec)%ElectronicState( 2, iQuaTemp) &
        + ElectronicState(1, iQua)*ElectronicState(2, iQua)
  END DO
  DO iQua = 1, nQuants -1
    SpecDSMC(iSpec)%ElectronicState( 2, iQua) = SpecDSMC(iSpec)%ElectronicState( 2, iQua) &
            / SpecDSMC(iSpec)%ElectronicState( 1, iQua)
  END DO
  SpecDSMC(iSpec)%ElectronicState( 1:2, 0) = ElectronicState(1:2,0)
  SpecDSMC(iSpec)%ElectronicState( 1:2, nQuants) = ElectronicState(1:2,dims(2)-1)
  SpecDSMC(iSpec)%MaxElecQuant  = SIZE( SpecDSMC(iSpec)%ElectronicState,2)
  LBWRITE(UNIT_StdOut,'(A,I5,A,I5,A,A,A)') 'Merged ',dims(2),' Electronic States to ',nQuants, ' for ',TRIM(dsetname),&
      ' (+1 for the ground state)'
END IF
! Close the file.
CALL H5FCLOSE_F(file_id_dsmc, err)
! Close FORTRAN interface.
CALL H5CLOSE_F(err)

! Check if the ground state is defined at 0K
IF(SpecDSMC(iSpec)%ElectronicState(2,0).NE.0.0) THEN
  CALL Abort(&
  __STAMP__,&
'ERROR in electronic energy levels: given ground state is not zero! Species: ', IntInfoOpt=iSpec)
END IF

IF (DSMC%ElectronicModel.EQ.4) THEN
  SpecDSMC(iSpec)%MaxMeanXiElec = 0.
  IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
    MaxTemp = SpecDSMC(iSpec)%ElectronicState(2,SpecDSMC(iSpec)%MaxElecQuant - 1)
    MinTemp = 0.
    MeanTemp = 0.5*(MaxTemp+MinTemp)
    MeanMax = 0.5*(MaxTemp+MeanTemp)
    DO WHILE(MaxTemp-MeanTemp.GT.0.1)
      XiMean=CalcXiElec(MeanTemp, iSpec)
      XiMeanMax=CalcXiElec(MeanMax, iSpec)
      IF (XiMean.GT.XiMeanMax) THEN
        MaxTemp = MeanTemp
      ELSE
        MinTemp = MeanTemp
      END IF
      MeanTemp = 0.5*(MaxTemp+MinTemp)
      MeanMax = 0.5*(MaxTemp+MeanTemp)
    END DO
    SpecDSMC(iSpec)%MaxMeanXiElec(1) = XiMean
    SpecDSMC(iSpec)%MaxMeanXiElec(2) = MeanTemp
  END IF
END IF

END SUBROUTINE ReadSpeciesLevel


SUBROUTINE SortEnergies(ElectronicState, nQuants)
!===================================================================================================================================
!>
!===================================================================================================================================
! use module
  USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)                                   :: nQuants
  REAL,INTENT(INOUT)                                    :: ElectronicState(1:2,0:nQuants-1)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                                               :: iStep, iLoop
  REAL                                                  :: TempState(2)
  LOGICAL                                               :: changed
!===================================================================================================================================
  iStep = nQuants
  changed = .true.
  DO WHILE (changed.OR.iStep.GT.1)
    changed = .false.
    IF (iStep.GT.1) THEN
      iStep = INT(iStep/1.3)
    END IF
    DO iLoop = 0, nQuants - 1 - iStep
      IF (ElectronicState(2,iLoop).GT.ElectronicState(2,iLoop+iStep)) THEN
        TempState(1:2) = ElectronicState(1:2,iLoop+iStep)
        ElectronicState(1:2,iLoop+iStep) = ElectronicState(1:2,iLoop)
        ElectronicState(1:2,iLoop) = TempState(1:2)
        changed = .true.
      END IF
    END DO
  END DO
END SUBROUTINE SortEnergies


SUBROUTINE CalcProbCorrFactorElec()
!===================================================================================================================================
!> Calculates the correction factor for electronical relaxation to achieve a landau-teller relaxation.
!===================================================================================================================================
! use module
  USE MOD_Globals_Vars,       ONLY: BoltzmannConst
  USE MOD_DSMC_Vars,          ONLY: DSMC, SpecDSMC, CollInf
  USE MOD_Particle_Vars,      ONLY: nSpecies
  USE MOD_Particle_Analyze_Tools, ONLY: CalcEelec
  USE MOD_part_tools              ,ONLY: CalcXiElec
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                                               :: iLoop, iSpec, jSpec
  REAL                                                  :: MinTemp, MaxTemp, TEqui, OldEn, NewEn, Xi_rel
  LOGICAL                                               :: doConverge
!===================================================================================================================================
DO iSpec = 1, nSpecies
  SpecDSMC(iSpec)%ElecRelaxCorrectFac = 0.
  IF((SpecDSMC(iSpec)%InterID.EQ.4).OR.(SpecDSMC(iSpec)%FullyIonized)) CYCLE
  DO jSpec = 1, nSpecies
    doConverge = .TRUE.
    Xi_rel = 2.*(2. - CollInf%omega(iSpec,jSpec))
    !Find Equilibrium Temperature
    MinTemp=MIN(DSMC%InstantTXiElec(1,iSpec),DSMC%InstantTransTemp(nSpecies+1))
    MaxTemp=MAX(DSMC%InstantTXiElec(1,iSpec),DSMC%InstantTransTemp(nSpecies+1))
    Tequi= 0.5*(MaxTemp+MinTemp)
    iLoop = 1
    OldEn =DSMC%InstantTXiElec(1,iSpec)*DSMC%InstantTXiElec(2,iSpec)+Xi_rel*DSMC%InstantTransTemp(nSpecies+1)
    NewEn=(CalcXiElec(TEqui, iSpec)+Xi_rel)*Tequi
    DO WHILE(ABS(OldEn-NewEn).GT.1E-2)
      IF (OldEn.LT.NewEn) THEN
       MaxTemp = TEqui
      ELSE
       MinTemp = TEqui
      END IF
      Tequi= 0.5*(MaxTemp+MinTemp)
      NewEn=(CalcXiElec(TEqui, iSpec)+Xi_rel)*Tequi
      iLoop = iLoop + 1
      IF (iLoop.EQ.100) THEN
        doConverge = .FALSE.
        EXIT
      END IF
    END DO
    IF (doConverge) THEN
      SpecDSMC(iSpec)%ElecRelaxCorrectFac(jSpec) = (CalcXiElec(DSMC%InstantTransTemp(nSpecies+1), iSpec)*DSMC%InstantTransTemp(nSpecies+1) &
        -DSMC%InstantTXiElec(1,iSpec)*DSMC%InstantTXiElec(2,iSpec))/(Xi_rel*(DSMC%InstantTransTemp(nSpecies+1)-TEqui))
      SpecDSMC(iSpec)%ElecRelaxCorrectFac(jSpec) = MAX(SpecDSMC(iSpec)%ElecRelaxCorrectFac(jSpec),1.)
    ELSE
      SpecDSMC(iSpec)%ElecRelaxCorrectFac(jSpec) = 1.
    END IF
  END DO
END DO
END SUBROUTINE CalcProbCorrFactorElec

REAL FUNCTION DiffElecEnergy(En1, En2)
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
  USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT(IN)         :: En1, En2
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  IF (ALMOSTEQUAL(En1,0.0)) CALL Abort(&
     __STAMP__,&
    'Energy of Electronic Shell is 0!!!')
  IF (ALMOSTEQUAL(En2,0.0)) CALL Abort(&
     __STAMP__,&
    'Energy of Electronic Shell is 0!!!')

  IF (ALMOSTEQUAL(En1,En2)) THEN
    DiffElecEnergy = 0.0
  ELSE
    DiffElecEnergy = (En2-En1)/En1
  END IF

  RETURN

END FUNCTION DiffElecEnergy

END MODULE MOD_DSMC_ElectronicModel
