!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

INTERFACE TVEEnergyExchange
  MODULE PROCEDURE TVEEnergyExchange
END INTERFACE

INTERFACE ReadSpeciesLevel
  MODULE PROCEDURE ReadSpeciesLevel
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: ElectronicEnergyExchange, InitElectronShell, TVEEnergyExchange, ReadSpeciesLevel, CalcXiElec
PUBLIC :: RelaxElectronicShellWall
!===================================================================================================================================
CONTAINS

SUBROUTINE InitElectronShell(iSpecies,iPart,iInit,init_or_sf)
!===================================================================================================================================
! init electronic shell
!===================================================================================================================================
  USE MOD_Globals,                ONLY : abort, myRank
  USE MOD_Globals_Vars,           ONLY : BoltzmannConst
  USE MOD_DSMC_Vars,              ONLY : SpecDSMC, PartStateIntEn, ElectronicDistriPart, DSMC
  USE MOD_Particle_Vars,          ONLY : Species, PEM
  USE MOD_Mesh_Vars               ,ONLY: offSetElem
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iPart, iSpecies, iInit, init_or_sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: iQua
  REAL                          :: iRan, ElectronicPartition, ElectronicPartitionTemp, iRan2, tmpExp
  REAL                        :: Telec                      ! electronic temperature
!===================================================================================================================================
  SELECT CASE (init_or_sf)
  CASE(1) !iInit
    IF (Species(iSpecies)%Init(iInit)%ElemTElecFileID.EQ.0) THEN
      TElec=SpecDSMC(iSpecies)%Init(iInit)%TElec
    ELSE
      TElec=Species(iSpecies)%Init(iInit)%ElemTElec(PEM%LocalElemID(iPart))
    END IF
  CASE(2) !SurfaceFlux
    Telec=SpecDSMC(iSpecies)%Surfaceflux(iInit)%Telec
  CASE DEFAULT
    CALL abort(&
    __STAMP__&
    ,'neither iInit nor Surfaceflux defined as reference!')
  END SELECT

  ElectronicPartition  = 0.
  ElectronicPartitionTemp = 0.
  ! calculate sum over all energy levels == partition function for temperature Telec
  IF (DSMC%ElectronicDistrModel) THEN
    IF(ALLOCATED(ElectronicDistriPart(iPart)%DistriFunc)) DEALLOCATE(ElectronicDistriPart(iPart)%DistriFunc)
    ALLOCATE(ElectronicDistriPart(iPart)%DistriFunc(1:SpecDSMC(iSpecies)%MaxElecQuant))
    PartStateIntEn(3,iPart) = 0.0
    DO iQua = 0, SpecDSMC(iSpecies)%MaxElecQuant - 1
      tmpExp = SpecDSMC(iSpecies)%ElectronicState(2,iQua) / Telec
      IF (CHECKEXP(tmpExp)) &
        ElectronicPartition = ElectronicPartition + SpecDSMC(iSpecies)%ElectronicState(1,iQua) * EXP(-tmpExp)
    END DO
    DO iQua = 0, SpecDSMC(iSpecies)%MaxElecQuant - 1
      tmpExp = SpecDSMC(iSpecies)%ElectronicState(2,iQua) / Telec
      IF (CHECKEXP(tmpExp)) THEN
        ElectronicDistriPart(iPart)%DistriFunc(iQua+1) = SpecDSMC(iSpecies)%ElectronicState(1,iQua)*EXP(-tmpExp)/ElectronicPartition
      ELSE
        ElectronicDistriPart(iPart)%DistriFunc(iQua+1) = 0.0
      END IF
      PartStateIntEn(3,iPart) = PartStateIntEn(3,iPart) + &
          ElectronicDistriPart(iPart)%DistriFunc(iQua+1) * BoltzmannConst * SpecDSMC(iSpecies)%ElectronicState(2,iQua)
    END DO
  ELSE
    DO iQua = 0, SpecDSMC(iSpecies)%MaxElecQuant - 1
      tmpExp = SpecDSMC(iSpecies)%ElectronicState(2,iQua) / Telec
      IF (CHECKEXP(tmpExp)) THEN
        ElectronicPartitionTemp = SpecDSMC(iSpecies)%ElectronicState(1,iQua) *EXP(-tmpExp)
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
      iQua = int( ( SpecDSMC(iSpecies)%MaxElecQuant ) * iRan)
      tmpExp = SpecDSMC(iSpecies)%ElectronicState(2,iQua) / Telec
      IF (CHECKEXP(tmpExp)) THEN
        ElectronicPartitionTemp = SpecDSMC(iSpecies)%ElectronicState(1,iQua) *EXP(-tmpExp)
      ELSE
        ElectronicPartitionTemp = 0.0
      END IF
      CALL RANDOM_NUMBER(iRan2)
    END DO
#if ( PP_TimeDiscMethod == 42 )
#ifdef CODE_ANALYZE
    SpecDSMC(iSpecies)%levelcounter(iQua) = SpecDSMC(iSpecies)%levelcounter(iQua) + 1
#endif
#endif
    PartStateIntEn(3,iPart) = BoltzmannConst * SpecDSMC(iSpecies)%ElectronicState(2,iQua)
  END IF
END SUBROUTINE InitElectronShell


FUNCTION RelaxElectronicShellWall(iPart,TWall)
!===================================================================================================================================
!> Function to determine the new electronic state of a particle at the wall temperature
!===================================================================================================================================
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
iSpec = PartSpecies(iPart)
IF (DSMC%ElectronicDistrModel) THEN
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
ELSE
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
END IF
END FUNCTION RelaxElectronicShellWall


SUBROUTINE ElectronicEnergyExchange(iPair,iPart1,FakXi, NewPart)
!===================================================================================================================================
! Electronic energy exchange
!===================================================================================================================================
  USE MOD_Globals
  USE MOD_DSMC_Vars,              ONLY : SpecDSMC, PartStateIntEn, RadialWeighting, Coll_pData, DSMC, ElectronicDistriPart  
  USE MOD_Particle_Vars,          ONLY : PartSpecies, VarTimeStep, usevMPF, nSpecies, PEM
  USE MOD_Globals_Vars,           ONLY : BoltzmannConst, Pi
  USE MOD_part_tools              ,ONLY: GetParticleWeight
  USE MOD_DSMC_Analyze,           ONLY: CalcTelec 
  USE MOD_Mesh_Vars              ,ONLY: OffSetElem
#if (PP_TimeDiscMethod==42)
  USE MOD_DSMC_Vars,              ONLY : DSMC
#endif
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iPair, iPart1
  REAL, INTENT(IN)              :: FakXi
  LOGICAL, INTENT(IN),OPTIONAL  :: NewPart
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: iQuaMax, MaxElecQuant, iQua, iSpecies
  REAL                          :: iRan, iRan2, gmax, gtemp, PartStateTemp, CollisionEnergy, ETraRel, TransElec, ElectronicPartition
  REAL                          :: Erel,  Eold, DistriOld(SpecDSMC(PartSpecies(iPart1))%MaxElecQuant), Etmp, tmpExp, LocRelaxProb
!===================================================================================================================================
  IF (DSMC%ElectronicDistrModel) THEN 
    iSpecies = PartSpecies(iPart1)   
    IF (PRESENT(NewPart)) THEN
      LocRelaxProb = 1.0
    ELSE
      LocRelaxProb = SpecDSMC(iSpecies)%ElecRelaxProb
    END IF
    Eold=  PartStateIntEn(3,iPart1)
    DistriOld(:) = ElectronicDistriPart(iPart1)%DistriFunc(:)
    ETraRel = Coll_pData(iPair)%Ec - PartStateIntEn(3,iPart1) * GetParticleWeight(iPart1)
    IF (usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
      ETraRel = ETraRel / GetParticleWeight(iPart1)
    END IF    
    TransElec = DSMC%InstantTransTemp(nSpecies + 1) !CalcTelec(ETraRel, iSpecies)
    IF (TransElec.LE.0.0) TransElec = 1./(BoltzmannConst*(FakXi+1.))*ETraRel
    ElectronicPartition = 0.0
    DO iQua = 0, SpecDSMC(iSpecies)%MaxElecQuant - 1
      tmpExp = SpecDSMC(iSpecies)%ElectronicState(2,iQua) / TransElec
      IF (CHECKEXP(tmpExp)) &
        ElectronicPartition = ElectronicPartition + SpecDSMC(iSpecies)%ElectronicState(1,iQua) * EXP (-tmpExp)
    END DO
    PartStateIntEn(3,iPart1) = 0.0
    DO iQua = 0, SpecDSMC(iSpecies)%MaxElecQuant - 1
      tmpExp = SpecDSMC(iSpecies)%ElectronicState(2,iQua) / TransElec
      IF (CHECKEXP(tmpExp)) THEN
        ElectronicDistriPart(iPart1)%DistriFunc(iQua+1) = &
          (1.-LocRelaxProb)*ElectronicDistriPart(iPart1)%DistriFunc(iQua+1) + &
          LocRelaxProb * SpecDSMC(iSpecies)%ElectronicState(1,iQua) *EXP (-tmpExp)/ElectronicPartition
      ELSE
        ElectronicDistriPart(iPart1)%DistriFunc(iQua+1) = (1.-LocRelaxProb)*ElectronicDistriPart(iPart1)%DistriFunc(iQua+1) 
      END IF
!      ElectronicDistriPart(iPart1)%DistriFunc(iQua+1) =  SpecDSMC(iSpecies)%ElectronicState(1,iQua) * &
!              EXP ( - SpecDSMC(iSpecies)%ElectronicState(2,iQua) / TransElec)/ElectronicPartition
      PartStateIntEn(3,iPart1) = PartStateIntEn(3,iPart1) + &
          ElectronicDistriPart(iPart1)%DistriFunc(iQua+1) * BoltzmannConst * SpecDSMC(iSpecies)%ElectronicState(2,iQua)
    END DO
    IF ((Coll_pData(iPair)%Ec-PartStateIntEn(3,iPart1)*GetParticleWeight(iPart1)).LT.0.0) then
      Etmp = (Coll_pData(iPair)%Ec - (1.-LocRelaxProb)*Eold*GetParticleWeight(iPart1))/(GetParticleWeight(iPart1)*LocRelaxProb)
      TransElec = CalcTelec(Etmp, iSpecies)*0.98
      ElectronicPartition = 0.0
      DO iQua = 0, SpecDSMC(iSpecies)%MaxElecQuant - 1
        tmpExp = SpecDSMC(iSpecies)%ElectronicState(2,iQua) / TransElec
        IF (CHECKEXP(tmpExp)) &
          ElectronicPartition = ElectronicPartition + SpecDSMC(iSpecies)%ElectronicState(1,iQua) * EXP (-tmpExp)
      END DO
      PartStateIntEn(3,iPart1) = 0.0
      DO iQua = 0, SpecDSMC(iSpecies)%MaxElecQuant - 1
        tmpExp = SpecDSMC(iSpecies)%ElectronicState(2,iQua) / TransElec
        IF (CHECKEXP(tmpExp)) THEN
          ElectronicDistriPart(iPart1)%DistriFunc(iQua+1) = (1.-LocRelaxProb)*DistriOld(iQua+1) &
              + LocRelaxProb * SpecDSMC(iSpecies)%ElectronicState(1,iQua) * EXP (-tmpExp)/ElectronicPartition
        ELSE
          ElectronicDistriPart(iPart1)%DistriFunc(iQua+1) = (1.-LocRelaxProb)*DistriOld(iQua+1)
        END IF
        PartStateIntEn(3,iPart1) = PartStateIntEn(3,iPart1) + &
            ElectronicDistriPart(iPart1)%DistriFunc(iQua+1) * BoltzmannConst * SpecDSMC(iSpecies)%ElectronicState(2,iQua)
      END DO
      IF ((Coll_pData(iPair)%Ec-PartStateIntEn(3,iPart1)*GetParticleWeight(iPart1)).LT.0.0) THEN
      CALL abort(&
      __STAMP__&
      ,'Negative collision energy after electronic excitation relaxation!')
      END IF
    END IF



!    IF (usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
!      CollisionEnergy = Coll_pData(iPair)%Ec / GetParticleWeight(iPart1)
!    ELSE
!      CollisionEnergy = Coll_pData(iPair)%Ec
!    END IF
!    iSpecies = PartSpecies(iPart1)
!    CALL RANDOM_NUMBER(iRan)  
!    PartStateIntEn(3,iPart1) = CollisionEnergy * (1.0 - iRan**(1.0/(FakXi+0.5*DSMC%InstantTXiElec(2,iSpecies))))



    
!    iSpecies = PartSpecies(iPart1)

!    TransElec = DSMC%InstantTransTemp(nSpecies + 1) 
!    Erel = (2. - SpecDSMC(iSpecies)%omegaVHS)*BoltzmannConst*TransElec
!    CollisionEnergy = Erel + PartStateIntEn(3,iPart1)

!!    ETraRel = Coll_pData(iPair)%Ec - PartStateIntEn(3,iPart1) * GetParticleWeight(iPart1)
!!    IF (usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
!!      ETraRel = ETraRel / GetParticleWeight(iPart1)
!!    END IF
!!    print*, ETraRel, Erel
!    CollisionEnergy = Erel + DSMC%InstantTXiElec(2,iSpecies)/2.*BoltzmannConst*DSMC%InstantTXiElec(1,iSpecies)

!    iQuaMax  = 0
!    ! Determine max electronic quant
!    MaxElecQuant = SpecDSMC(PartSpecies(iPart1))%MaxElecQuant - 1
!    PartStateTemp = CollisionEnergy / BoltzmannConst
!    DO iQua = 0, MaxElecQuant
!      IF (PartStateTemp - SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua).GT.0.) THEN
!        gmax = gmax + SpecDSMC(PartSpecies(iPart1))%ElectronicState(1,iQua) * &
!                ( CollisionEnergy - BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua))**FakXi 
!        ! maximal possible Quant before term goes negative
!      ELSE
!        EXIT
!      END IF
!    END DO
!    IF(gmax.LE.0.0) THEN
!      PartStateIntEn(3,iPart1) = 0.0
!      RETURN
!    END IF
!    PartStateIntEn(3,iPart1) = 0.0
!    DO iQua = 0, SpecDSMC(iSpecies)%MaxElecQuant - 1
!      IF (PartStateTemp - SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua).GT.0.) THEN
!        ElectronicDistriPart(iPart1)%DistriFunc(iQua+1) = SpecDSMC(PartSpecies(iPart1))%ElectronicState(1,iQua) * &
!                ( CollisionEnergy - BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua))**FakXi/gmax 
!        PartStateIntEn(3,iPart1) = PartStateIntEn(3,iPart1) + &
!            ElectronicDistriPart(iPart1)%DistriFunc(iQua+1) * BoltzmannConst * SpecDSMC(iSpecies)%ElectronicState(2,iQua)
!      ELSE
!        ElectronicDistriPart(iPart1)%DistriFunc(iQua+1) = 0.0
!      END IF
!    END DO
  ELSE
    IF (usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
      CollisionEnergy = Coll_pData(iPair)%Ec / GetParticleWeight(iPart1)
    ELSE
      CollisionEnergy = Coll_pData(iPair)%Ec
    END IF

    iQuaMax  = 0
    ! Determine max electronic quant
    MaxElecQuant = SpecDSMC(PartSpecies(iPart1))%MaxElecQuant - 1
    ! determine maximal Quant and term according to Eq (7) of Liechty
    gmax = 0.
    PartStateTemp = CollisionEnergy / BoltzmannConst
    DO iQua = 0, MaxElecQuant
      IF (PartStateTemp - SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua).GT.0.) THEN
        gtemp = SpecDSMC(PartSpecies(iPart1))%ElectronicState(1,iQua) * &
                ( CollisionEnergy - BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua))**FakXi
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
    gtemp = SpecDSMC(PartSpecies(iPart1))%ElectronicState(1,iQua) * &
            ( CollisionEnergy - BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua))**FakXi
    CALL RANDOM_NUMBER(iRan2)
    ! acceptance-rejection for iQuaElec
    DO WHILE ( iRan2 .GE. gtemp / gmax )
      CALL RANDOM_NUMBER(iRan)
      iQua = int( ( iQuaMax +1 ) * iRan)
      gtemp = SpecDSMC(PartSpecies(iPart1))%ElectronicState(1,iQua) * &
              ( CollisionEnergy - BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua))**FakXi
      CALL RANDOM_NUMBER(iRan2)
    END DO
#if (PP_TimeDiscMethod==42)
! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
  IF (.NOT.DSMC%ReservoirSimuRate) THEN
#endif
   PartStateIntEn(3,iPart1) = BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua)
#if (PP_TimeDiscMethod==42)
  END IF
#endif
  END IF
END SUBROUTINE ElectronicEnergyExchange


SUBROUTINE TVEEnergyExchange(CollisionEnergy,iPart1,FakXi)
!===================================================================================================================================
! Electronic energy exchange
!===================================================================================================================================
  USE MOD_DSMC_Vars,              ONLY : DSMC, SpecDSMC, PartStateIntEn
  USE MOD_Particle_Vars,          ONLY : PartSpecies
  USE MOD_Globals_Vars,           ONLY : BoltzmannConst
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iPart1
  REAL, INTENT(IN)              :: FakXi
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL, INTENT(INOUT)           :: CollisionEnergy                                                !
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: iQuaMax, MaxElecQuant, iQua   ! , iQuaMax3
  INTEGER                       :: jQVib, QMaxVib
  REAL                          :: iRan, iRan2, gmax, gtemp, PartStateTemp, iRanVib
!#if ( PP_TimeDiscMethod==42 )
!  INTEGER                       :: iQuaold
!#endif
!===================================================================================================================================
  iQuaMax  = 0
  ! Determine max electronic quant
  MaxElecQuant = SpecDSMC(PartSpecies(iPart1))%MaxElecQuant - 1
!#if ( PP_TimeDiscMethod==42 )
!  iQuaold=0
!  ! determine old Quant
!  DO iQua = 0, MaxElecQuant
!    IF ( PartStateIntEn(3,iPart1) / BoltzmannConst .ge. &
!      SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua) ) THEN
!      iQuaold = iQua
!    ELSE
!    ! exit loop
!      EXIT
!    END IF
!  END DO
!#endif
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
#if (PP_TimeDiscMethod==42)
! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
  IF (.NOT.DSMC%ReservoirSimuRate) THEN
#endif
   PartStateIntEn(3,iPart1) = BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua)
   PartStateIntEn(1,iPart1) = (jQVib + DSMC%GammaQuant) * BoltzmannConst &
                  * SpecDSMC(PartSpecies(iPart1))%CharaTVib
#if (PP_TimeDiscMethod==42)
  END IF
#endif

!#if (PP_TimeDiscMethod==42)
!    ! list of number of particles in each energy level
!  IF (.NOT.DSMC%ReservoirSimuRate) THEN
!    SpecDSMC(PartSpecies(iPart1))%levelcounter(iQuaold) = SpecDSMC(PartSpecies(iPart1))%levelcounter(iQuaold) - 1
!    SpecDSMC(PartSpecies(iPart1))%levelcounter(iQua)    = SpecDSMC(PartSpecies(iPart1))%levelcounter(iQua)    + 1
!    SpecDSMC(PartSpecies(iPart1))%dtlevelcounter(iQua)  = SpecDSMC(PartSpecies(iPart1))%dtlevelcounter(iQua)  + 1
!  END IF
!  ! collision with X resulting in a transition from i to j
!  IF ( present(iPart2) .AND. (.NOT.usevMPF) ) THEN
!  SpecDSMC(PartSpecies(iPart1))%ElectronicTransition(PartSpecies(iPart2),iQuaold,iQua) = &
!                                SpecDSMC(PartSpecies(iPart1))%ElectronicTransition(PartSpecies(iPart2),iQuaold,iQua) + 1
!  END IF
!#endif
END SUBROUTINE TVEEnergyExchange


SUBROUTINE ReadSpeciesLevel ( Dsetname, iSpec )
!===================================================================================================================================
! Subroutine to read the electronic levels from DSMCSpeciesElectronicState.h5
!===================================================================================================================================
! use module
  USE MOD_io_hdf5
  USE MOD_Globals
  USE MOD_DSMC_Vars,            ONLY: DSMC, SpecDSMC
  USE MOD_HDF5_Input,           ONLY: DatasetExists
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
  REAL                                                  :: tempEnergyDiff, tempEnergy
  LOGICAL                                               :: DataSetFound
!===================================================================================================================================
  SWRITE(UNIT_StdOut,'(A)') 'Read electronic level entries '//TRIM(dsetname)//' from '//TRIM(DSMC%ElectronicModelDatabase)
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
    SWRITE(UNIT_StdOut,'(A,I5,A,I5,A,A,A)') 'Merged ',dims(2),' Electronic States to ',nQuants, ' for ',TRIM(dsetname),&
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

END SUBROUTINE ReadSpeciesLevel


SUBROUTINE SortEnergies(ElectronicState, nQuants)
!===================================================================================================================================
! Subroutine to read the electronic levels from DSMCSpeciesElectronicState.h5
!===================================================================================================================================
! use module
  USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL,INTENT(INOUT)                                    :: ElectronicState(1:2,0:nQuants-1)
  INTEGER, INTENT(IN)                                   :: nQuants
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

REAL FUNCTION DiffElecEnergy(En1, En2)
!===================================================================================================================================
! Calculates area of mesh element
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


PURE REAL FUNCTION CalcXiElec(Telec, iSpec)
!===================================================================================================================================
!> Calculation of the electronic degree of freedom for a given temperature and species
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars               ,ONLY: SpecDSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: Telec  !
INTEGER, INTENT(IN)             :: iSpec      ! Number of Species
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                         :: iQua
REAL                            :: TempRatio, SumOne, SumTwo
!===================================================================================================================================

SumOne = 0.0
SumTwo = 0.0

DO iQua = 0, SpecDSMC(iSpec)%MaxElecQuant-1
  TempRatio = SpecDSMC(iSpec)%ElectronicState(2,iQua)/Telec
  IF(CHECKEXP(TempRatio)) THEN
    SumOne = SumOne + SpecDSMC(iSpec)%ElectronicState(1,iQua)*SpecDSMC(iSpec)%ElectronicState(2,iQua)*EXP(-TempRatio)
    SumTwo = SumTwo + SpecDSMC(iSpec)%ElectronicState(1,iQua)*EXP(-TempRatio)
  END IF
END DO

IF((SumOne.GT.0.0).AND.(SumTwo.GT.0.0)) THEN
  CalcXiElec = 2. * SumOne / (SumTwo * Telec)
ELSE
  CalcXiElec = 0.0
END IF

RETURN

END FUNCTION CalcXiElec

END MODULE MOD_DSMC_ElectronicModel
