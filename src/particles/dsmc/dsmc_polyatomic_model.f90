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

MODULE MOD_DSMC_PolyAtomicModel
!===================================================================================================================================
! Routines for the treatment of polyatomic molecules
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE

INTERFACE DSMC_VibRelaxPoly
  MODULE PROCEDURE DSMC_VibRelaxPoly_ARM_MH
END INTERFACE

INTERFACE DSMC_SetInternalEnr_Poly
  MODULE PROCEDURE DSMC_SetInternalEnr_Poly_ARM_SingleMode
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: InitPolyAtomicMolecs, DSMC_SetInternalEnr_Poly_ARM, DSMC_SetInternalEnr_Poly_MH, DSMC_SetInternalEnr_Poly_MH_FirstPick
PUBLIC :: DSMC_RotRelaxPoly, DSMC_VibRelaxPoly_ARM, DSMC_VibRelaxPoly_MH, DSMC_VibRelaxPoly_ARM_MH
PUBLIC :: DSMC_FindFirstVibPick, DSMC_RelaxVibPolyProduct
!===================================================================================================================================

CONTAINS

SUBROUTINE InitPolyAtomicMolecs(iSpec)
!===================================================================================================================================
! Initialization of variables for polyatomic molecules
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,           ONLY : BoltzmannConst
USE MOD_DSMC_Vars,              ONLY : DSMC, SpecDSMC, PolyatomMolDSMC
USE MOD_ReadInTools
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           ::  iSpec
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  CHARACTER(32)                  :: hilf, hilf2
  INTEGER                        :: iPolyatMole, iVibDOF
!===================================================================================================================================
WRITE(UNIT=hilf,FMT='(I0)') iSpec
iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
PolyatomMolDSMC(iPolyatMole)%LinearMolec = GETLOGICAL('Part-Species'//TRIM(hilf)//'-LinearMolec')
IF (PolyatomMolDSMC(iPolyatMole)%LinearMolec) THEN
  SpecDSMC(iSpec)%Xi_Rot = 2
ELSE
  SpecDSMC(iSpec)%Xi_Rot = 3
END IF
PolyatomMolDSMC(iPolyatMole)%NumOfAtoms = GETINT('Part-Species'//TRIM(hilf)//'-NumOfAtoms')
! TSHO not implemented with polyatomic molecules, but Ediss_eV required for the calculation of polyatomic temp. (upper bound)
SpecDSMC(iSpec)%Ediss_eV   = GETREAL('Part-Species'//TRIM(hilf)//'-Ediss_eV','0.')
IF(SpecDSMC(iSpec)%Ediss_eV.EQ.0.) THEN
  CALL abort(__STAMP__,'ERROR in Polyatomic Species-Ini: Missing dissociation energy, Species: ',iSpec)
END IF
PolyatomMolDSMC(iPolyatMole)%VibDOF = 3*PolyatomMolDSMC(iPolyatMole)%NumOfAtoms - 3 - SpecDSMC(iSpec)%Xi_Rot
ALLOCATE(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(PolyatomMolDSMC(iPolyatMole)%VibDOF))
IF(DSMC%PolySingleMode) THEN
  ALLOCATE(PolyatomMolDSMC(iPolyatMole)%GammaVib(PolyatomMolDSMC(iPolyatMole)%VibDOF))
  PolyatomMolDSMC(iPolyatMole)%GammaVib(:) = 0.0
  ALLOCATE(PolyatomMolDSMC(iPolyatMole)%VibRelaxProb(PolyatomMolDSMC(iPolyatMole)%VibDOF))
  PolyatomMolDSMC(iPolyatMole)%VibRelaxProb(:) = 0.0
END IF
! Read-in of characteristic rotational temperature
ALLOCATE(PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(3))
IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec) THEN
  PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1) = GETREAL('Part-Species'//TRIM(hilf)//'-CharaTempRot','0')
  PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(2:3) = 1
ELSE
  DO iVibDOF = 1,3
    WRITE(UNIT=hilf2,FMT='(I0)') iVibDOF
    PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(iVibDOF) = &
      GETREAL('Part-Species'//TRIM(hilf)//'-CharaTempRot'//TRIM(hilf2),'0')
  END DO
END IF
! Read-in of characteristic vibrational temperature and calculation of zero-point energy
DO iVibDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
  WRITE(UNIT=hilf2,FMT='(I0)') iVibDOF
  PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iVibDOF) = GETREAL('Part-Species'//TRIM(hilf)//'-CharaTempVib'//TRIM(hilf2))
  SpecDSMC(iSpec)%EZeroPoint = SpecDSMC(iSpec)%EZeroPoint &
    + DSMC%GammaQuant*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iVibDOF)*BoltzmannConst
END DO
ALLOCATE(PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(PolyatomMolDSMC(iPolyatMole)%VibDOF))
! Maximum number of quantum number per DOF cut at 80 to reduce computational effort
PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(1:PolyatomMolDSMC(iPolyatMole)%VibDOF) = 80

END SUBROUTINE InitPolyAtomicMolecs


SUBROUTINE DSMC_FindFirstVibPick(iInitTmp, iSpec, init_or_sf)
!===================================================================================================================================
! Burn-in phase for the modified Metropolis-Hasting method for the particle generation of polyatomic molecules
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars,            ONLY : SpecDSMC, PolyatomMolDSMC
USE MOD_Particle_Vars,        ONLY : Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iInitTmp, iSpec, init_or_sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL, ALLOCATABLE            :: iRan(:)
REAL                          :: iRan2, NormProb
INTEGER,ALLOCATABLE          :: iQuant_old(:)
INTEGER                       :: iDOF, iWalk, iPolyatMole, iInit
REAL                          :: TVib                       ! vibrational temperature
!===================================================================================================================================

SELECT CASE (init_or_sf)
  CASE(1) !iInit
    TVib=SpecDSMC(iSpec)%Init(iInitTmp)%TVib
  CASE(2) !SurfaceFlux
    TVib=SpecDSMC(iSpec)%SurfaceFlux(iInitTmp)%TVib
    iInit = iInitTmp + Species(iSpec)%NumberOfInits
  CASE DEFAULT
    CALL abort(&
    __STAMP__&
    ,'Neither iInit nor SurfaceFlux defined as reference!')
END SELECT

iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
ALLOCATE(iRan(PolyatomMolDSMC(iPolyatMole)%VibDOF),iQuant_old(PolyatomMolDSMC(iPolyatMole)%VibDOF))

CALL RANDOM_NUMBER(iRan)
PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(:,iInit) = INT(iRan(:)*PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(:))

DO iWalk=1, 5000
  iQuant_old(:)=PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(:,iInit)
  CALL RANDOM_NUMBER(iRan)
  PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(:,iInit) = iQuant_old(:)+FLOOR(3*iRan(:)-1)
  NormProb = 0.0
  DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
    IF(PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(iDOF,iInit).LT.0) THEN
      PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(iDOF,iInit) = -1*PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(iDOF,iInit) -1
    END IF
    NormProb = NormProb + (iQuant_old(iDOF) &
        -PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(iDOF,iInit))*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) &
        /TVib
  END DO
  NormProb = MIN(1.0,EXP(NormProb))
  CALL RANDOM_NUMBER(iRan2)
  IF (NormProb.LT.iRan2) PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(:,iInit)=iQuant_old(:)
END DO

DEALLOCATE(iRan, iQuant_old)

END SUBROUTINE DSMC_FindFirstVibPick


SUBROUTINE DSMC_SetInternalEnr_Poly_ARM_SingleMode(iSpecies, iInit, iPart, init_or_sf)
!===================================================================================================================================
! Initialization of polyatomic molecules by treating every mode separately in a loop
!===================================================================================================================================
! MODULES
USE MOD_Globals               ,ONLY: Abort
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, SpecDSMC, DSMC,PolyatomMolDSMC,VibQuantsPar,BGGas
USE MOD_Particle_Vars         ,ONLY: PEM, Species
USE MOD_Particle_Sampling_Vars,ONLY: AdaptBCMacroVal, AdaptBCMapElemToSample
USE MOD_DSMC_ElectronicModel  ,ONLY: InitElectronShell
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpecies, iInit, iPart, init_or_sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: iRan, iRan2, NormProb
INTEGER                       :: iQuant, iDOF, iPolyatMole
REAL                          :: TVib                       ! vibrational temperature
REAL                          :: TRot                       ! rotational temperature
INTEGER                       :: ElemID
!===================================================================================================================================
ElemID = PEM%LocalElemID(iPart)
SELECT CASE (init_or_sf)
  CASE(1) !iInit
    TVib=SpecDSMC(iSpecies)%Init(iInit)%TVib
    TRot=SpecDSMC(iSpecies)%Init(iInit)%TRot
  CASE(2) !SurfaceFlux
    IF(Species(iSpecies)%Surfaceflux(iInit)%Adaptive) THEN
      SELECT CASE(Species(iSpecies)%Surfaceflux(iInit)%AdaptiveType)
        CASE(1,3,4) ! Pressure and massflow inlet (pressure/massflow, temperature const)
          TVib=SpecDSMC(iSpecies)%Surfaceflux(iInit)%TVib
          TRot=SpecDSMC(iSpecies)%Surfaceflux(iInit)%TRot
        CASE(2) ! adaptive Outlet/freestream
          TVib = Species(iSpecies)%Surfaceflux(iInit)%AdaptivePressure &
                  / (BoltzmannConst * AdaptBCMacroVal(4,AdaptBCMapElemToSample(ElemID),iSpecies))
          TRot = TVib
        CASE DEFAULT
          CALL abort(&
          __STAMP__&
          ,'Wrong adaptive type for Surfaceflux in vib/rot poly!')
      END SELECT
    ELSE
      TVib=SpecDSMC(iSpecies)%Surfaceflux(iInit)%TVib
      TRot=SpecDSMC(iSpecies)%Surfaceflux(iInit)%TRot
    END IF
  CASE DEFAULT
    CALL abort(&
    __STAMP__&
    ,'Neither iInit nor SurfaceFlux defined as reference!')
END SELECT

! Background gas distribution
IF(BGGas%NumberOfSpecies.GT.0) THEN
  IF(BGGas%BackgroundSpecies(iSpecies).AND.BGGas%UseDistribution) THEN
    TVib = BGGas%Distribution(BGGas%MapSpecToBGSpec(iSpecies),DSMC_TVIB,ElemID)
    TRot = BGGas%Distribution(BGGas%MapSpecToBGSpec(iSpecies),DSMC_TROT,ElemID)
  END IF
END IF

IF (DSMC%ElectronicModel.GT.0) THEN
  CALL InitElectronShell(iSpecies,iPart,iInit,init_or_sf)
ENDIF

! set vibrational energy
iPolyatMole = SpecDSMC(iSpecies)%SpecToPolyArray
IF(ALLOCATED(VibQuantsPar(iPart)%Quants)) DEALLOCATE(VibQuantsPar(iPart)%Quants)
ALLOCATE(VibQuantsPar(iPart)%Quants(PolyatomMolDSMC(iPolyatMole)%VibDOF))
PartStateIntEn( 1,iPart) = 0.0
DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
  CALL RANDOM_NUMBER(iRan)
  iQuant = INT(-LOG(iRan)*TVib/PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF))
  DO WHILE (iQuant.GE.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(iDOF))
    CALL RANDOM_NUMBER(iRan)
    iQuant = INT(-LOG(iRan)*TVib/PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF))
  END DO
  PartStateIntEn( 1,iPart) = PartStateIntEn( 1,iPart) &
                              + (iQuant + DSMC%GammaQuant)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
  VibQuantsPar(iPart)%Quants(iDOF)=iQuant
END DO
IF (SpecDSMC(iSpecies)%Xi_Rot.EQ.2) THEN
  CALL RANDOM_NUMBER(iRan2)
  PartStateIntEn( 2,iPart) = -BoltzmannConst*TRot*LOG(iRan2)
ELSE IF (SpecDSMC(iSpecies)%Xi_Rot.EQ.3) THEN
  CALL RANDOM_NUMBER(iRan2)
  PartStateIntEn( 2,iPart) = iRan2*10 !the distribution function has only non-negligible  values betwenn 0 and 10
  NormProb = SQRT(PartStateIntEn( 2,iPart))*EXP(-PartStateIntEn( 2,iPart))/(SQRT(0.5)*EXP(-0.5))
  CALL RANDOM_NUMBER(iRan2)
  DO WHILE (iRan2.GE.NormProb)
    CALL RANDOM_NUMBER(iRan2)
    PartStateIntEn( 2,iPart) = iRan2*10 !the distribution function has only non-negligible  values betwenn 0 and 10
    NormProb = SQRT(PartStateIntEn( 2,iPart))*EXP(-PartStateIntEn( 2,iPart))/(SQRT(0.5)*EXP(-0.5))
    CALL RANDOM_NUMBER(iRan2)
  END DO
  PartStateIntEn( 2,iPart) = PartStateIntEn( 2,iPart)*BoltzmannConst*TRot
END IF

END SUBROUTINE DSMC_SetInternalEnr_Poly_ARM_SingleMode


SUBROUTINE DSMC_SetInternalEnr_Poly_ARM(iSpec, iInit, iPart, init_or_sf)
!===================================================================================================================================
! Initialization/particle generation of polyatomic molecules with the acceptance-rejection method (extremely slow for molecules with
! more than 3 atoms due to low acceptance probability, only for comparison with Metropolis-Hastings)
!===================================================================================================================================
! MODULES
USE MOD_Globals               ,ONLY: Abort
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, SpecDSMC, DSMC,PolyatomMolDSMC,VibQuantsPar
USE MOD_Particle_Vars         ,ONLY: PEM
USE MOD_DSMC_ElectronicModel  ,ONLY: InitElectronShell
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, INTENT(IN)           :: iSpec, iInit, iPart, init_or_sf
REAL, ALLOCATABLE             :: iRan(:), tempEng(:)
REAL                          :: iRan2, NormProb
INTEGER,ALLOCATABLE           :: iQuant(:)
INTEGER                       :: iDOF,iPolyatMole
REAL                          :: TVib                       ! vibrational temperature
REAL                          :: TRot                       ! rotational temperature
INTEGER                       :: ElemID
!===================================================================================================================================
ElemID = PEM%LocalElemID(iPart)
SELECT CASE (init_or_sf)
  CASE(1) !iInit
    TVib=SpecDSMC(iSpec)%Init(iInit)%TVib
    TRot=SpecDSMC(iSpec)%Init(iInit)%TRot
  CASE(2) !SurfaceFlux
    TVib=SpecDSMC(iSpec)%SurfaceFlux(iInit)%TVib
    TRot=SpecDSMC(iSpec)%SurfaceFlux(iInit)%TRot
  CASE DEFAULT
    CALL abort(&
    __STAMP__&
    ,'Neither iInit nor SurfaceFlux defined as reference!')
END SELECT

IF (DSMC%ElectronicModel.GT.0) THEN
  CALL InitElectronShell(iSpec,iPart,iInit,init_or_sf)
ENDIF

! Set vibrational energy of new molecule
IF (SpecDSMC(iSpec)%PolyatomicMol) THEN
  iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
  IF(ALLOCATED(VibQuantsPar(iPart)%Quants)) DEALLOCATE(VibQuantsPar(iPart)%Quants)
  ALLOCATE(VibQuantsPar(iPart)%Quants(PolyatomMolDSMC(iPolyatMole)%VibDOF))
  ALLOCATE(iRan(PolyatomMolDSMC(iPolyatMole)%VibDOF) &
          ,tempEng(PolyatomMolDSMC(iPolyatMole)%VibDOF) &
          ,iQuant(PolyatomMolDSMC(iPolyatMole)%VibDOF))
  CALL RANDOM_NUMBER(iRan)
  iQuant(:) = INT(iRan(:)*PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(:))
  tempEng(:)=iQuant(:)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:)/TVib
  NormProb = 1.0
  DO iDOF = 1 , PolyatomMolDSMC(iPolyatMole)%VibDOF
    NormProb = NormProb*EXP(-tempEng(iDOF))
  END DO
  CALL RANDOM_NUMBER(iRan2)
  DO WHILE (iRan2.GE.NormProb)
    CALL RANDOM_NUMBER(iRan)
    iQuant(:) = INT(iRan(:)*PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(:))
    tempEng(:)=iQuant(:)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:)/TVib
    NormProb = 1.0
    DO iDOF = 1 , PolyatomMolDSMC(iPolyatMole)%VibDOF
      NormProb = NormProb*EXP(-tempEng(iDOF))
    END DO
    CALL RANDOM_NUMBER(iRan2)
  END DO
  PartStateIntEn( 1,iPart) = 0.0
  VibQuantsPar(iPart)%Quants(:)=iQuant(:)
  DO iDOF = 1 , PolyatomMolDSMC(iPolyatMole)%VibDOF
    PartStateIntEn( 1,iPart)= PartStateIntEn( 1,iPart) &
      +(iQuant(iDOF) + DSMC%GammaQuant)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
  END DO
! Set rotational energy of new molecule
  IF (SpecDSMC(iSpec)%Xi_Rot.EQ.2) THEN
    CALL RANDOM_NUMBER(iRan2)
    PartStateIntEn( 2,iPart) = -BoltzmannConst*TRot*LOG(iRan2)
  ELSE IF (SpecDSMC(iSpec)%Xi_Rot.EQ.3) THEN
    CALL RANDOM_NUMBER(iRan2)
    PartStateIntEn( 2,iPart) = iRan2*10 !the distribution function has only non-negligible  values betwenn 0 and 10
    NormProb = SQRT(PartStateIntEn( 2,iPart))*EXP(-PartStateIntEn( 2,iPart))/(SQRT(0.5)*EXP(-0.5))
    CALL RANDOM_NUMBER(iRan2)
    DO WHILE (iRan2.GE.NormProb)
      CALL RANDOM_NUMBER(iRan2)
      PartStateIntEn( 2,iPart) = iRan2*10 !the distribution function has only non-negligible  values betwenn 0 and 10
      NormProb = SQRT(PartStateIntEn( 2,iPart))*EXP(-PartStateIntEn( 2,iPart))/(SQRT(0.5)*EXP(-0.5))
      CALL RANDOM_NUMBER(iRan2)
    END DO
    PartStateIntEn( 2,iPart) = PartStateIntEn( 2,iPart)*BoltzmannConst*TRot
  END IF
  DEALLOCATE(iRan, tempEng, iQuant)
ELSE
  PartStateIntEn( 1,iPart) = 0
  PartStateIntEn( 2,iPart) = 0
END IF

END SUBROUTINE DSMC_SetInternalEnr_Poly_ARM


SUBROUTINE DSMC_SetInternalEnr_Poly_MH_FirstPick(iSpec, iInit, iPart, init_or_sf)
!===================================================================================================================================
! Initialization/particle generation of polyatomic molecules with modified Metropolis-Hasting method
! Burn-in phase is included for each particle, can be utilized for setting the internal energy regardless of the previous state
!===================================================================================================================================
! MODULES
USE MOD_Globals               ,ONLY: Abort
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, SpecDSMC, DSMC,PolyatomMolDSMC,VibQuantsPar
USE MOD_Particle_Vars         ,ONLY: PEM
USE MOD_DSMC_ElectronicModel  ,ONLY: InitElectronShell
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec, iInit, iPart, init_or_sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL, ALLOCATABLE             :: iRan(:)
REAL                          :: iRan2, NormProb
INTEGER,ALLOCATABLE           :: iQuant(:), iQuant_old(:)
INTEGER                       :: iDOF,iPolyatMole, iWalk
REAL                          :: TVib                       ! vibrational temperature
REAL                          :: TRot                       ! rotational temperature
INTEGER                       :: ElemID
!===================================================================================================================================
ElemID = PEM%LocalElemID(iPart)
  SELECT CASE (init_or_sf)
    CASE(1) !iInit
      TVib=SpecDSMC(iSpec)%Init(iInit)%TVib
      TRot=SpecDSMC(iSpec)%Init(iInit)%TRot
    CASE(2) !SurfaceFlux
      TVib=SpecDSMC(iSpec)%SurfaceFlux(iInit)%TVib
      TRot=SpecDSMC(iSpec)%SurfaceFlux(iInit)%TRot
    CASE DEFAULT
      CALL abort(&
      __STAMP__&
      ,'Neither iInit nor SurfaceFlux defined as reference!')
  END SELECT

  IF (DSMC%ElectronicModel.GT.0) THEN
    CALL InitElectronShell(iSpec,iPart,iInit,init_or_sf)
  ENDIF

  IF (SpecDSMC(iSpec)%PolyatomicMol) THEN
    iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
    ALLOCATE(iRan(PolyatomMolDSMC(iPolyatMole)%VibDOF) &
            ,iQuant(PolyatomMolDSMC(iPolyatMole)%VibDOF) &
            ,iQuant_old(PolyatomMolDSMC(iPolyatMole)%VibDOF))
    IF(ALLOCATED(VibQuantsPar(iPart)%Quants)) DEALLOCATE(VibQuantsPar(iPart)%Quants)
    ALLOCATE(VibQuantsPar(iPart)%Quants(PolyatomMolDSMC(iPolyatMole)%VibDOF))

    CALL RANDOM_NUMBER(iRan)
    iQuant(:) = INT(iRan(:)*PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(:))

    DO iWalk=1, 4000
      iQuant_old(:)=iQuant(:)
      CALL RANDOM_NUMBER(iRan)
      iQuant(:) = iQuant_old(:)+FLOOR(3*iRan(:)-1)
      NormProb = 0.0
      DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
        IF(iQuant(iDOF).LT.0) iQuant(iDOF) = -1*iQuant(iDOF) -1
        NormProb = NormProb + (iQuant_old(iDOF)-iQuant(iDOF))*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / TVib
      END DO
      NormProb = MIN(1.0,EXP(NormProb))
      CALL RANDOM_NUMBER(iRan2)
      IF (NormProb.LT.iRan2) iQuant(:)=iQuant_old(:)
    END DO

    PartStateIntEn( 1,iPart) = 0.0
    DO iDOF = 1 , PolyatomMolDSMC(iPolyatMole)%VibDOF
      PartStateIntEn( 1,iPart)= PartStateIntEn( 1,iPart) &
        +(iQuant(iDOF) + DSMC%GammaQuant)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
    END DO
    VibQuantsPar(iPart)%Quants(:)=iQuant(:)
    DEALLOCATE(iRan, iQuant,iQuant_old)

   !set rotational energy
    IF (SpecDSMC(iSpec)%Xi_Rot.EQ.2) THEN
      CALL RANDOM_NUMBER(iRan2)
      PartStateIntEn( 2,iPart) = -BoltzmannConst*TRot*LOG(iRan2)
    ELSE IF (SpecDSMC(iSpec)%Xi_Rot.EQ.3) THEN
      CALL RANDOM_NUMBER(iRan2)
      PartStateIntEn( 2,iPart) = iRan2*10 !the distribution function has only non-negligible  values betwenn 0 and 10
      NormProb = SQRT(PartStateIntEn( 2,iPart))*EXP(-PartStateIntEn( 2,iPart))/(SQRT(0.5)*EXP(-0.5))
      CALL RANDOM_NUMBER(iRan2)
      DO WHILE (iRan2.GE.NormProb)
        CALL RANDOM_NUMBER(iRan2)
        PartStateIntEn( 2,iPart) = iRan2*10 !the distribution function has only non-negligible  values betwenn 0 and 10
        NormProb = SQRT(PartStateIntEn( 2,iPart))*EXP(-PartStateIntEn( 2,iPart))/(SQRT(0.5)*EXP(-0.5))
        CALL RANDOM_NUMBER(iRan2)
      END DO
      PartStateIntEn( 2,iPart) = PartStateIntEn( 2,iPart)*BoltzmannConst*TRot
    END IF
  ELSE
    PartStateIntEn( 1,iPart) = 0
    PartStateIntEn( 2,iPart) = 0
  END IF

END SUBROUTINE DSMC_SetInternalEnr_Poly_MH_FirstPick


SUBROUTINE DSMC_SetInternalEnr_Poly_MH(iSpec, iInitTmp, iPart, init_or_sf)
!===================================================================================================================================
! Initialization/particle generation of polyatomic molecules with modified Metropolis-Hasting method
! Burn-in phase is NOT included, utilizes LastVibQuantNums as first initial value of the Markov chain
!===================================================================================================================================
! MODULES
USE MOD_Globals               ,ONLY: Abort
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, SpecDSMC, DSMC,PolyatomMolDSMC,VibQuantsPar
USE MOD_Particle_Vars         ,ONLY: Species, PEM
USE MOD_DSMC_ElectronicModel  ,ONLY: InitElectronShell
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec, iInitTmp, iPart, init_or_sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL, ALLOCATABLE             :: iRan(:)
REAL                          :: iRan2, NormProb
INTEGER,ALLOCATABLE           :: iQuant_old(:)
INTEGER                       :: iDOF,iPolyatMole, iWalk, iInit
REAL                          :: TVib                       ! vibrational temperature
REAL                          :: TRot                       ! rotational temperature
INTEGER                       :: ElemID
!===================================================================================================================================
ElemID = PEM%LocalElemID(iPart)
SELECT CASE (init_or_sf)
  CASE(1) !iInit
    TVib=SpecDSMC(iSpec)%Init(iInitTmp)%TVib
    TRot=SpecDSMC(iSpec)%Init(iInitTmp)%TRot
  CASE(2) !SurfaceFlux
    TVib=SpecDSMC(iSpec)%SurfaceFlux(iInitTmp)%TVib
    TRot=SpecDSMC(iSpec)%SurfaceFlux(iInitTmp)%TRot
    iInit = iInitTmp + Species(iSpec)%NumberOfInits
  CASE DEFAULT
    CALL abort(&
    __STAMP__&
    ,'Neither iInit nor SurfaceFlux defined as reference!')
END SELECT

IF (DSMC%ElectronicModel.GT.0) THEN
  CALL InitElectronShell(iSpec,iPart,iInit,init_or_sf)
ENDIF

IF (SpecDSMC(iSpec)%PolyatomicMol) THEN
  iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
  ALLOCATE(iRan(PolyatomMolDSMC(iPolyatMole)%VibDOF) &
          ,iQuant_old(PolyatomMolDSMC(iPolyatMole)%VibDOF))
  IF(ALLOCATED(VibQuantsPar(iPart)%Quants)) DEALLOCATE(VibQuantsPar(iPart)%Quants)
  ALLOCATE(VibQuantsPar(iPart)%Quants(PolyatomMolDSMC(iPolyatMole)%VibDOF))
! Set vibrational energy
  DO iWalk = 1, 150
    iQuant_old(:)=PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(:, iInit)
    CALL RANDOM_NUMBER(iRan)
    PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(:, iInit) = iQuant_old(:)+FLOOR(3*iRan(:)-1)
    NormProb = 0.0
    DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
      IF(PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(iDOF, iInit).LT.0) THEN
          PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(iDOF, iInit) = &
            -1*PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(iDOF, iInit) -1
      END IF
      NormProb = NormProb + (iQuant_old(iDOF) &
        -PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(iDOF, iInit))*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / TVib
    END DO
    NormProb = MIN(1.0,EXP(NormProb))
    CALL RANDOM_NUMBER(iRan2)
    IF (NormProb.LT.iRan2) PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(:, iInit)=iQuant_old(:)
  END DO
  PartStateIntEn( 1,iPart) = 0.0
  DO iDOF = 1 , PolyatomMolDSMC(iPolyatMole)%VibDOF
    PartStateIntEn( 1,iPart)= PartStateIntEn( 1,iPart) &
      +(PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(iDOF, iInit) &
      + DSMC%GammaQuant)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
  END DO
  VibQuantsPar(iPart)%Quants(:)=PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(:, iInit)
  DEALLOCATE(iRan, iQuant_old)
! Set rotational energy
  IF (SpecDSMC(iSpec)%Xi_Rot.EQ.2) THEN
    CALL RANDOM_NUMBER(iRan2)
    PartStateIntEn( 2,iPart) = -BoltzmannConst*TRot*LOG(iRan2)
  ELSE IF (SpecDSMC(iSpec)%Xi_Rot.EQ.3) THEN
    CALL RANDOM_NUMBER(iRan2)
    PartStateIntEn( 2,iPart) = iRan2*10 !the distribution function has only non-negligible  values betwenn 0 and 10
    NormProb = SQRT(PartStateIntEn( 2,iPart))*EXP(-PartStateIntEn( 2,iPart))/(SQRT(0.5)*EXP(-0.5))
    CALL RANDOM_NUMBER(iRan2)
    DO WHILE (iRan2.GE.NormProb)
      CALL RANDOM_NUMBER(iRan2)
      PartStateIntEn( 2,iPart) = iRan2*10 !the distribution function has only non-negligible  values betwenn 0 and 10
      NormProb = SQRT(PartStateIntEn( 2,iPart))*EXP(-PartStateIntEn( 2,iPart))/(SQRT(0.5)*EXP(-0.5))
      CALL RANDOM_NUMBER(iRan2)
    END DO
    PartStateIntEn( 2,iPart) = PartStateIntEn( 2,iPart)*BoltzmannConst*TRot
  END IF
ELSE
  PartStateIntEn( 1,iPart) = 0
  PartStateIntEn( 2,iPart) = 0
END IF

END SUBROUTINE DSMC_SetInternalEnr_Poly_MH


SUBROUTINE DSMC_RelaxVibPolyProduct(iPair, iPart, FakXi, Xi_Vib, WeightProd)
!===================================================================================================================================
! Initialization of the vibrational state of polyatomic molecules created during chemical reactions
! Single mode initialization analagous to DSMC_SetInternalEnr_Poly_ARM_SingleMode
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_Globals_Vars,         ONLY : BoltzmannConst
  USE MOD_DSMC_Vars,            ONLY : PartStateIntEn, SpecDSMC, DSMC, PolyatomMolDSMC, Coll_pData, VibQuantsPar
  USE MOD_Particle_Vars,        ONLY : PartSpecies
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iPart, iPair
  REAL, INTENT(IN)              :: Xi_Vib(:)
  REAL, INTENT(INOUT)           :: FakXi
  REAL, INTENT(IN), OPTIONAL    :: WeightProd
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                          :: iRan, MaxColQua, Weight
  INTEGER                       :: iQua, iQuaMax, iDOF, iPolyatMole
!===================================================================================================================================
  iPolyatMole = SpecDSMC(PartSpecies(iPart))%SpecToPolyArray
  IF (PRESENT(WeightProd)) THEN
    Weight = WeightProd
  ELSE
    Weight = 1.
  END IF
  IF(ALLOCATED(VibQuantsPar(iPart)%Quants)) DEALLOCATE(VibQuantsPar(iPart)%Quants)
  ALLOCATE(VibQuantsPar(iPart)%Quants(PolyatomMolDSMC(iPolyatMole)%VibDOF))
  PartStateIntEn( 1,iPart) = 0.0
  DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
    ! Addition of the zero-point energy part for the respective dofs (avoiding the redistribution of too much vibrational energy)
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec  &
        + DSMC%GammaQuant * BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*Weight
    ! Maximum quantum number calculated with the collision energy
    MaxColQua = Coll_pData(iPair)%Ec/(Weight*BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF))  &
              - DSMC%GammaQuant
    iQuaMax = MIN(INT(MaxColQua) + 1, PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(iDOF))
    CALL RANDOM_NUMBER(iRan)
    iQua = INT(iRan * iQuaMax)
    CALL RANDOM_NUMBER(iRan)
    DO WHILE (iRan.GT.(1 - iQua/MaxColQua)**FakXi)
     !laux diss page 31
     CALL RANDOM_NUMBER(iRan)
     iQua = INT(iRan * iQuaMax)
     CALL RANDOM_NUMBER(iRan)
    END DO
    PartStateIntEn(1,iPart) = PartStateIntEn(1,iPart)     &
      + (iQua + DSMC%GammaQuant) * BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec &
        - (iQua + DSMC%GammaQuant) * BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*Weight
    VibQuantsPar(iPart)%Quants(iDOF) = iQua
    IF (iDOF.LT.PolyatomMolDSMC(iPolyatMole)%VibDOF) FakXi = FakXi - 0.5*Xi_vib(iDOF + 1)
  END DO
END SUBROUTINE DSMC_RelaxVibPolyProduct


SUBROUTINE DSMC_VibRelaxPoly_ARM(iPair, iPart, FakXi)
!===================================================================================================================================
! Vibrational relaxation routine with the acceptance rejection method (slower than Metropolis-Hasting for molecules with more than
! three atoms, use only for comparison)
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, SpecDSMC, PolyatomMolDSMC,VibQuantsPar, Coll_pData, RadialWeighting
USE MOD_Particle_Vars         ,ONLY: PartSpecies, UseVarTimeStep, usevMPF
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_part_tools            ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPair, iPart
REAL, INTENT(IN)              :: FakXi
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL, ALLOCATABLE             :: iRan(:), tempEng(:)
REAL                          :: iRan2, NormProb, tempProb, Ec
INTEGER,ALLOCATABLE           :: iQuant(:), iMaxQuant(:)
INTEGER                       :: iDOF,iPolyatMole, iSpec
!===================================================================================================================================
iSpec = PartSpecies(iPart)
IF (usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep) THEN
  Ec = Coll_pData(iPair)%Ec / GetParticleWeight(iPart)
ELSE
  Ec = Coll_pData(iPair)%Ec
END IF
iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
ALLOCATE(iRan(PolyatomMolDSMC(iPolyatMole)%VibDOF) &
        ,tempEng(PolyatomMolDSMC(iPolyatMole)%VibDOF) &
        ,iQuant(PolyatomMolDSMC(iPolyatMole)%VibDOF) &
        ,iMaxQuant(PolyatomMolDSMC(iPolyatMole)%VibDOF))
NormProb = Ec
DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
  NormProb = NormProb - 0.5*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
END DO
NormProb = NormProb**FakXi
iMaxQuant(:) = INT(Ec/(BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))-0.5) + 1
iMaxQuant(:) = MIN(iMaxQuant(:), PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(:))

DO
  CALL RANDOM_NUMBER(iRan)
  iQuant(:)=INT(iRan(:)*iMaxQuant(:))
  tempEng(:)=(iQuant(:) + 0.5)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:)*BoltzmannConst
  tempProb = 0.0
  DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
    tempProb = tempProb + tempEng(iDOF)
  END DO
  IF (Ec-tempProb.GE.0.0) THEN
    CALL RANDOM_NUMBER(iRan2)
    IF (iRan2.LE.((Ec-tempProb)**FakXi/NormProb)) EXIT
  END IF
END DO
PartStateIntEn(1,iPart)=tempProb
VibQuantsPar(iPart)%Quants(:) = iQuant(:)

DEALLOCATE(iRan ,tempEng ,iQuant ,iMaxQuant)

END SUBROUTINE DSMC_VibRelaxPoly_ARM


SUBROUTINE DSMC_VibRelaxPoly_MH(iPair, iPart,FakXi)
!===================================================================================================================================
! Vibrational relaxation routine with the Metropolis-Hastings method (no burn-in phase)
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, SpecDSMC, PolyatomMolDSMC,VibQuantsPar, Coll_pData, RadialWeighting
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_Particle_Vars         ,ONLY: PartSpecies, UseVarTimeStep, usevMPF
USE MOD_part_tools            ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPair, iPart
REAL, INTENT(IN)              :: FakXi
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL, ALLOCATABLE             :: iRan(:), tempEng(:)
REAL                          :: iRan2, NormProb, tempProb, Ec
INTEGER,ALLOCATABLE           :: iQuant(:), iMaxQuant(:)
INTEGER                       :: iDOF,iPolyatMole, iWalk, iSpec
!===================================================================================================================================
iSpec = PartSpecies(iPart)
IF (usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep) THEN
  Ec = Coll_pData(iPair)%Ec / GetParticleWeight(iPart)
ELSE
  Ec = Coll_pData(iPair)%Ec
END IF
iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
ALLOCATE(iRan(PolyatomMolDSMC(iPolyatMole)%VibDOF) &
        ,tempEng(PolyatomMolDSMC(iPolyatMole)%VibDOF) &
        ,iQuant(PolyatomMolDSMC(iPolyatMole)%VibDOF) &
        ,iMaxQuant(PolyatomMolDSMC(iPolyatMole)%VibDOF))
DO iWalk=1,750
  NormProb = Ec - PartStateIntEn(1,iPart)
  ! Proper modelling of energy transfer between old and new state in chemistry
  NormProb = NormProb**FakXi

  CALL RANDOM_NUMBER(iRan)
  iQuant(:) = VibQuantsPar(iPart)%Quants(:)+FLOOR(3*iRan(:)-1)
  DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
    IF(iQuant(iDOF).LT.0) iQuant(iDOF) = -1*iQuant(iDOF) -1
  END DO

  tempEng(:)=(iQuant(:) + 0.5)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:)*BoltzmannConst
  tempProb = Ec
  DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
    tempProb = tempProb - tempEng(iDOF)
  END DO
  IF(tempProb.GT.0) THEN
    NormProb = MIN(1.0,tempProb**FakXi/NormProb)
    CALL RANDOM_NUMBER(iRan2)
    IF(NormProb.GE.iRan2) THEN
      PartStateIntEn(1,iPart) = 0.0
      DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
        PartStateIntEn(1,iPart) = PartStateIntEn(1,iPart) + tempEng(iDOF)
      END DO
      VibQuantsPar(iPart)%Quants(:) = iQuant(:)
    END IF
  END IF
END DO

END SUBROUTINE DSMC_VibRelaxPoly_MH


SUBROUTINE DSMC_VibRelaxPoly_GibbsSampling(iPair, iPart, FakXi)
!===================================================================================================================================
! Vibrational relaxation (multi-mode) using Gibbs sampling
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, SpecDSMC, PolyatomMolDSMC,VibQuantsPar, Coll_pData, RadialWeighting
USE MOD_Particle_Vars         ,ONLY: PartSpecies, UseVarTimeStep, usevMPF
USE MOD_part_tools            ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPair, iPart
REAL, INTENT(IN)              :: FakXi
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: iRan, iRan2, NormProb, tempProb, Ec, NormProbZero
INTEGER,ALLOCATABLE           :: iQuant(:), iMaxQuant(:)
INTEGER                       :: iDOF, iDOF2, iPolyatMole, iSpec, iLoop
!===================================================================================================================================
iSpec = PartSpecies(iPart)
IF (usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep) THEN
  Ec = Coll_pData(iPair)%Ec / GetParticleWeight(iPart)
ELSE
  Ec = Coll_pData(iPair)%Ec
END IF
iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
ALLOCATE(iQuant(PolyatomMolDSMC(iPolyatMole)%VibDOF),iMaxQuant(PolyatomMolDSMC(iPolyatMole)%VibDOF))
iQuant(:) = VibQuantsPar(iPart)%Quants(:)

NormProbZero = Ec - SpecDSMC(iSpec)%EZeroPoint

iMaxQuant(:) = INT(Ec/(BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))-0.5) + 1
iMaxQuant(:) = MIN(iMaxQuant(:), PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(:))

DO iLoop = 1,4
  DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
    DO
      CALL RANDOM_NUMBER(iRan)
      iQuant(iDOF) = INT(iRan*iMaxQuant(iDOF))
      tempProb = SpecDSMC(iSpec)%EZeroPoint
      DO iDOF2 = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
        tempProb = tempProb + iQuant(iDOF2)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF2)*BoltzmannConst
        IF(iDOF2.NE.iDOF) NormProb = NormProbZero - iQuant(iDOF2)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF2)*BoltzmannConst
      END DO
      CALL RANDOM_NUMBER(iRan2)
      IF ((Ec-tempProb).GE.0.0) THEN
        IF ((iRan2.GT.((Ec-tempProb)/NormProb)**FakXi)) THEN
          DO
            tempProb = tempProb - iQuant(iDOF)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
            CALL RANDOM_NUMBER(iRan)
            iQuant(iDOF) = INT(iRan*iMaxQuant(iDOF))
            tempProb = tempProb  + iQuant(iDOF)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
            CALL RANDOM_NUMBER(iRan2)
            IF ((Ec-tempProb).GE.0.0) THEN
              IF ((iRan2.LE.((Ec-tempProb)/NormProb)**FakXi)) EXIT
            END IF
          END DO
        END IF
        IF ((iRan2.LE.((Ec-tempProb)/NormProb)**FakXi)) EXIT
      END IF
    END DO
  END DO
END DO

PartStateIntEn(1,iPart) = tempProb
VibQuantsPar(iPart)%Quants(:) = iQuant(:)

DEALLOCATE(iQuant ,iMaxQuant)

END SUBROUTINE DSMC_VibRelaxPoly_GibbsSampling


SUBROUTINE DSMC_VibRelaxPoly_ARM_MH(iPair, iPart,FakXi)
!===================================================================================================================================
! Switch between ARM and MH/Gibbs depending on the number of vibrational modes:
! Acceptance Rejection: up to 4 modes (molecules with 3 atoms, linear and non-linear)
! Metropolis-Hastings: from 6 modes (molecules with 4 or more atoms)
! Gibbs sampling: strong dependence on the number of iterations for duration and accuracy, has to be tested more thoroughly for
!                 different molecules
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,            ONLY : SpecDSMC, PolyatomMolDSMC
  USE MOD_Particle_Vars,        ONLY : PartSpecies
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iPair, iPart
  REAL, INTENT(IN)              :: FakXi
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: iPolyatMole
!===================================================================================================================================

  iPolyatMole = SpecDSMC(PartSpecies(iPart))%SpecToPolyArray
  IF(PolyatomMolDSMC(iPolyatMole)%VibDOF.GT.5) THEN
    CALL DSMC_VibRelaxPoly_MH(iPair,iPart,FakXi)
    ! CALL DSMC_VibRelaxPoly_GibbsSampling(iPair,iPart,FakXi)
  ELSE
    CALL DSMC_VibRelaxPoly_ARM(iPair,iPart,FakXi)
  END IF

END SUBROUTINE DSMC_VibRelaxPoly_ARM_MH


SUBROUTINE DSMC_VibRelaxPolySingle(iPair, iPart, FakXi, DOFRelax)
!===================================================================================================================================
! Vibrational relaxation routine for polyatomic molecules, only treating a single given vibrational mode with ARM
! NOTE: Not compatible for radial weighting yet.
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, SpecDSMC, PolyatomMolDSMC, VibQuantsPar, Coll_pData, DSMC, RadialWeighting
USE MOD_Particle_Vars         ,ONLY: PartSpecies, usevMPF, UseVarTimeStep
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_part_tools            ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPair, iPart, DOFRelax
REAL, INTENT(IN)              :: FakXi
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: iRan, MaxColQua, Ec
INTEGER                       :: iPolyatMole, iQua, iQuaMax
!===================================================================================================================================
! Not all vibrational energy is redistributed but only the energy of the selected vibrational degree of freedom
Ec = Coll_pData(iPair)%Ec - PartStateIntEn(1,iPart)*GetParticleWeight(iPart)

iPolyatMole = SpecDSMC(PartSpecies(iPart))%SpecToPolyArray
! Adding the vibrational energy of the selected vibrational mode DOFRelax
Ec = Ec + (VibQuantsPar(iPart)%Quants(DOFRelax) + DSMC%GammaQuant) * BoltzmannConst  &
                                                * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(DOFRelax)*GetParticleWeight(iPart)
IF (usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep) THEN
  Ec = Ec / GetParticleWeight(iPart)
END IF
! Determining the maximal quantum number with the available collision energy
MaxColQua = Ec/(BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(DOFRelax)) - DSMC%GammaQuant
iQuaMax = MIN(INT(MaxColQua) + 1, PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(DOFRelax))
! Get the new vibrational quantum number
CALL RANDOM_NUMBER(iRan)
iQua = INT(iRan * iQuaMax)
CALL RANDOM_NUMBER(iRan)
DO WHILE (iRan.GT.(1 - iQua/MaxColQua)**FakXi)
  !laux diss page 31
  CALL RANDOM_NUMBER(iRan)
  iQua = INT(iRan * iQuaMax)
  CALL RANDOM_NUMBER(iRan)
END DO
! Setting the new vibrational state
PartStateIntEn(1,iPart) = PartStateIntEn(1,iPart) &
  ! Substracting the old energy of the specific mode
  - VibQuantsPar(iPart)%Quants(DOFRelax) * BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(DOFRelax) &
  ! Adding the new energy of the specific mode
  + iQua * BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(DOFRelax)

! Saving the vibrational quantum number
VibQuantsPar(iPart)%Quants(DOFRelax) = iQua

END SUBROUTINE DSMC_VibRelaxPolySingle


SUBROUTINE DSMC_RotRelaxPoly(iPair, iPart,FakXi)
!===================================================================================================================================
! Rotational relaxation routine
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,            ONLY : PartStateIntEn, Coll_pData
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iPart, iPair
  REAL, INTENT(IN)              :: FakXi
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                          :: iRan2, NormProb, tempProb, fak1, fak2, Ec
!===================================================================================================================================
  Ec = Coll_pData(iPair)%Ec
  fak1 = (3.0/2.0+FakXi-1.0)/(3.0/2.0-1.0)
  fak2 = (3.0/2.0+FakXi-1.0)/(FakXi)

  CALL RANDOM_NUMBER(iRan2)
  tempProb = Ec*iRan2
  NormProb = ((fak1*tempProb/Ec)**(3.0/2.0-1.0))*((fak2*(1.0-tempProb/Ec))**(FakXi))
  CALL RANDOM_NUMBER(iRan2)
  DO WHILE (iRan2.GE.NormProb)
    CALL RANDOM_NUMBER(iRan2)
    tempProb = Ec*iRan2
    NormProb = (fak1*tempProb/Ec)**(3.0/2.0-1.0)*(fak2*(1.0-tempProb/Ec))**(FakXi)
    CALL RANDOM_NUMBER(iRan2)
  END DO
  PartStateIntEn(2,iPart)=tempProb

END SUBROUTINE DSMC_RotRelaxPoly


END MODULE MOD_DSMC_PolyAtomicModel
