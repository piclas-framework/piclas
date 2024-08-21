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
PUBLIC :: DSMC_FindFirstVibPick, DSMC_RelaxVibPolyProduct, DSMC_SetInternalEnr
!===================================================================================================================================

CONTAINS

SUBROUTINE InitPolyAtomicMolecs(iSpec)
!===================================================================================================================================
!> Initialization of variables for polyatomic molecules
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars      ,ONLY: BoltzmannConst, PlanckConst, PI
USE MOD_DSMC_Vars         ,ONLY: DSMC, SpecDSMC, PolyatomMolDSMC
USE MOD_ReadInTools
USE MOD_PARTICLE_Vars     ,ONLY: Species, SpeciesDatabase
USE MOD_io_hdf5
USE MOD_HDF5_input        ,ONLY: ReadAttribute, DatasetExists, AttributeExists
#if USE_MPI
USE MOD_LoadBalance_Vars  ,ONLY: PerformLoadBalance
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)       :: iSpec
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(32)             :: hilf, hilf2
CHARACTER(LEN=64)         :: dsetname
INTEGER                   :: iPolyatMole, iVibDOF, IntToLog, err
INTEGER(HID_T)            :: file_id_specdb                       ! File identifier
LOGICAL                   :: AttrExists
!===================================================================================================================================

LBWRITE (UNIT_stdOut,'(68(". "))')
WRITE(UNIT=hilf,FMT='(I0)') iSpec
iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray

IF(SpeciesDatabase.NE.'none') THEN
  ! Initialize FORTRAN interface.
  CALL H5OPEN_F(err)
  CALL H5FOPEN_F (TRIM(SpeciesDatabase), H5F_ACC_RDONLY_F, file_id_specdb, err)
  dsetname = TRIM('/Species/'//TRIM(Species(iSpec)%Name))
  ! Linear molecule
  CALL ReadAttribute(file_id_specdb,'LinearMolec',1,DatasetName = dsetname,IntScalar=IntToLog)
  IF(IntToLog.EQ.1) THEN
    PolyatomMolDSMC(iPolyatMole)%LinearMolec = .TRUE.
  ELSE
    PolyatomMolDSMC(iPolyatMole)%LinearMolec = .FALSE.
  END IF
  CALL PrintOption('LinearMolec, '//TRIM(Species(iSpec)%Name),'DB',LogOpt=PolyatomMolDSMC(iPolyatMole)%LinearMolec)
  ! Number of atoms
  CALL ReadAttribute(file_id_specdb,'NumOfAtoms',1,DatasetName = dsetname,IntScalar=PolyatomMolDSMC(iPolyatMole)%NumOfAtoms)
  CALL PrintOption('NumOfAtoms, '//TRIM(Species(iSpec)%Name),'DB',IntOpt=PolyatomMolDSMC(iPolyatMole)%NumOfAtoms)
  ! Dissociation energy
  ! TSHO not implemented with polyatomic molecules, but Ediss_eV required for the calculation of polyatomic temp. (upper bound)
  CALL ReadAttribute(file_id_specdb,'Ediss_eV',1,DatasetName = dsetname,RealScalar=SpecDSMC(iSpec)%Ediss_eV)
  CALL PrintOption('Ediss_eV, '//TRIM(Species(iSpec)%Name),'DB',RealOpt=SpecDSMC(iSpec)%Ediss_eV)
  ! Close the file.
  CALL H5FCLOSE_F(file_id_specdb, err)
  ! Close FORTRAN interface.
  CALL H5CLOSE_F(err)
END IF

IF(Species(iSpec)%DoOverwriteParameters) THEN
  PolyatomMolDSMC(iPolyatMole)%LinearMolec = GETLOGICAL('Part-Species'//TRIM(hilf)//'-LinearMolec')
  PolyatomMolDSMC(iPolyatMole)%NumOfAtoms = GETINT('Part-Species'//TRIM(hilf)//'-NumOfAtoms')
  SpecDSMC(iSpec)%Ediss_eV   = GETREAL('Part-Species'//TRIM(hilf)//'-Ediss_eV')
END IF

IF(SpecDSMC(iSpec)%Ediss_eV.EQ.0.) THEN
  CALL abort(__STAMP__,'ERROR in Polyatomic Species-Ini: Missing dissociation energy, Species: ',iSpec)
END IF

IF (PolyatomMolDSMC(iPolyatMole)%LinearMolec) THEN
  SpecDSMC(iSpec)%Xi_Rot = 2
ELSE
  SpecDSMC(iSpec)%Xi_Rot = 3
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
ALLOCATE(PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3))

IF(SpeciesDatabase.NE.'none') THEN
  ! Initialize FORTRAN interface.
  CALL H5OPEN_F(err)
  CALL H5FOPEN_F (TRIM(SpeciesDatabase), H5F_ACC_RDONLY_F, file_id_specdb, err)
  dsetname = TRIM('/Species/'//TRIM(Species(iSpec)%Name))
  IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec) THEN
    CALL AttributeExists(file_id_specdb,'MomentOfInertia',TRIM(dsetname), AttrExists=AttrExists)
    IF (AttrExists) THEN
      CALL ReadAttribute(file_id_specdb,'MomentOfInertia',3,DatasetName = dsetname,  &
        RealArray=PolyatomMolDSMC(iPolyatMole)%MomentOfInertia)
      CALL AttributeExists(file_id_specdb,'CharaTempRot',TRIM(dsetname), AttrExists=AttrExists)
      IF(AttrExists)THEN
        CALL ReadAttribute(file_id_specdb,'CharaTempRot',1,DatasetName = dsetname,RealScalar=SpecDSMC(iSpec)%CharaTRot)
      ELSE
        PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1) = PlanckConst**2 / &
        (8 * PI**2 * PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) * BoltzmannConst)
      END IF
    ELSE  ! moment of inertia not found
      CALL AttributeExists(file_id_specdb,'CharaTempRot',TRIM(dsetname), AttrExists=AttrExists)
      IF(AttrExists)THEN
        CALL ReadAttribute(file_id_specdb,'CharaTempRot',1,DatasetName = dsetname,RealScalar=SpecDSMC(iSpec)%CharaTRot)
      ELSE  ! CharaTempRot not found
        PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1) = 0
      END IF
      IF(DSMC%DoRotRelaxQuantized)THEN
        CALL abort(&
        __STAMP__&
        ,'Moment of inertia necessary for quantized rotational energy and is not set for species', iSpec)
      ELSE  ! not needed if DoRotRelaxQuantized false
        PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) = 0
      END IF
    END IF
    CALL PrintOption('CharaTempRot, '//TRIM(Species(iSpec)%Name),'DB',RealOpt=PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))
    CALL PrintOption('MomentOfInertia, '//TRIM(Species(iSpec)%Name),'DB',RealOpt=PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1))
    ! set dummy values
    PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(2:3) = 0
    PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(2:3) = 0
  ELSE
    CALL AttributeExists(file_id_specdb,'MomentOfInertia',TRIM(dsetname), AttrExists=AttrExists)
    IF (AttrExists) THEN
      CALL ReadAttribute(file_id_specdb,'MomentOfInertia',3,DatasetName = dsetname,  &
        RealArray=PolyatomMolDSMC(iPolyatMole)%MomentOfInertia)
      DO iVibDOF = 1,3
        WRITE(UNIT=hilf2,FMT='(I0)') iVibDOF
        CALL AttributeExists(file_id_specdb,TRIM('CharaTempRot'//TRIM(hilf2)),TRIM(dsetname), AttrExists=AttrExists)
        IF (AttrExists) THEN  !read in CharaTempRot
          CALL ReadAttribute(file_id_specdb,TRIM('CharaTempRot'//TRIM(hilf2)),1,DatasetName = dsetname,RealScalar=PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(iVibDOF))
        ELSE  !calculate CharaTempRot
          PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(iVibDOF) = PlanckConst**2 / &
            (8 * PI**2 * PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(iVibDOF) * BoltzmannConst)
        END IF
      END DO
    ELSE  ! moment of inertia not found
      DO iVibDOF = 1,3
        WRITE(UNIT=hilf2,FMT='(I0)') iVibDOF
        CALL AttributeExists(file_id_specdb,TRIM('CharaTempRot'//TRIM(hilf2)),TRIM(dsetname), AttrExists=AttrExists)
        IF (AttrExists) THEN  ! read in CharaTempRot
          CALL ReadAttribute(file_id_specdb,TRIM('CharaTempRot'//TRIM(hilf2)),1,DatasetName = dsetname,RealScalar=PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(iVibDOF))
        ELSE  ! set dummy value to zero
          PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(iVibDOF) = 0
        END IF
      END DO
      IF(DSMC%DoRotRelaxQuantized)THEN
        CALL abort(&
        __STAMP__&
        ,'Moment of inertia necessary for quantized rotational energy and is not set for species', iSpec)
      ELSE  ! not needed if DoRotRelaxQuantized false
        DO iVibDOF = 1,3
          PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(iVibDOF) = 0
        END DO
      END IF
    END IF
    ! Print out variables
    DO iVibDOF = 1,3
      WRITE(UNIT=hilf2,FMT='(I0)') iVibDOF
      CALL PrintOption('MomentOfInertia'//TRIM(hilf2)//' '//TRIM(Species(iSpec)%Name),'DB', &
      RealOpt=PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(iVibDOF))
    END DO
    DO iVibDOF = 1,3
      WRITE(UNIT=hilf2,FMT='(I0)') iVibDOF
      CALL PrintOption('CharaTempRot'//TRIM(hilf2)//' '//TRIM(Species(iSpec)%Name),'DB', &
        RealOpt=PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(iVibDOF))
    END DO
    !//TODO: more checks?
    ! sanity checks for order of moments of inertia
    IF((PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3).EQ.PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(2).AND.  &
      PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3).NE.PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)))THEN
      ! C should be different one (for symmetric tops) - for sperical tops all should be the same
      CALL abort(&
      __STAMP__&
      ,'Moments of inertia in wrong order in database for species', iSpec)
    END IF
  END IF
  ! Read-in of characteristic vibrational temperature
  DO iVibDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
    WRITE(UNIT=hilf2,FMT='(I0)') iVibDOF
    CALL ReadAttribute(file_id_specdb,TRIM('CharaTempVib'//TRIM(hilf2)),1,DatasetName = dsetname, &
      RealScalar=PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iVibDOF))
    CALL PrintOption('CharaTempVib'//TRIM(hilf2)//' '//TRIM(Species(iSpec)%Name),'DB',  &
      RealOpt=PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iVibDOF))
  END DO
  ! Close the file.
  CALL H5FCLOSE_F(file_id_specdb, err)
  ! Close FORTRAN interface.
  CALL H5CLOSE_F(err)
END IF

IF(Species(iSpec)%DoOverwriteParameters) THEN
  IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec) THEN
    PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)   = GETREAL('Part-Species'//TRIM(hilf)//'-MomentOfInertia')
    PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(2:3) = 0
    PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)                   = GETREAL('Part-Species'//TRIM(hilf)//'-CharaTempRot','-1')
    IF(PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1).EQ.-1)THEN
      PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)                 = PlanckConst**2 / &
        (8 * PI**2 * PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) * BoltzmannConst)
    END IF
    PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(2:3)    = 0
  ELSE
    DO iVibDOF = 1,3
      WRITE(UNIT=hilf2,FMT='(I0)') iVibDOF
      PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(iVibDOF) = GETREAL('Part-Species'//TRIM(hilf)//'-MomentOfInertia'//TRIM(hilf2))
      PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(iVibDOF)                 = GETREAL('Part-Species'//TRIM(hilf)//'-CharaTempRot','-1')
      IF(PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(iVibDOF).EQ.-1)THEN
        PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(iVibDOF)               = PlanckConst**2 / &
          (8 * PI**2 * PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(iVibDOF) * BoltzmannConst)
      END IF
    END DO
  END IF
  ! Read-in of characteristic vibrational temperature
  DO iVibDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
    WRITE(UNIT=hilf2,FMT='(I0)') iVibDOF
    PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iVibDOF) = GETREAL('Part-Species'//TRIM(hilf)//'-CharaTempVib'//TRIM(hilf2))
  END DO
END IF

! Calculation of zero-point energy
DO iVibDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
  SpecDSMC(iSpec)%EZeroPoint = SpecDSMC(iSpec)%EZeroPoint + DSMC%GammaQuant*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iVibDOF)*BoltzmannConst
END DO

! save Rotational Group of molecule
IF( PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1).EQ.PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(2).AND. &
    PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(2).EQ.PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3).AND. &
    .NOT.PolyatomMolDSMC(iPolyatMole)%LinearMolec)THEN
  ! sphrical top with A=B=C
  PolyatomMolDSMC(iPolyatMole)%RotationalGroup = 1
ELSE IF(PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1).EQ.PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(2).AND. &
        PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3).NE.PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1).AND. &
        PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3).NE.PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(2))THEN
  ! symmetric top with A=B,
  PolyatomMolDSMC(iPolyatMole)%RotationalGroup = 2
ELSE IF(PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1).NE.PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(2).AND. &
        PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(2).NE.PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3).AND. &
        PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1).NE.PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3))THEN
  ! asymmetric top with all moments of inertia different
  PolyatomMolDSMC(iPolyatMole)%RotationalGroup = 3
ELSE ! set dummy to catch false cases
  PolyatomMolDSMC(iPolyatMole)%RotationalGroup = -1
END IF

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


SUBROUTINE DSMC_SetInternalEnr(iSpec, iInit, iPart, init_or_sf)
!===================================================================================================================================
!> Energy distribution according to dissertation of Laux (diatomic)
!===================================================================================================================================
! MODULES
USE MOD_Globals                 ,ONLY: abort
USE MOD_Globals_Vars            ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars               ,ONLY: PartStateIntEn, SpecDSMC, DSMC, BGGas
USE MOD_Particle_Vars           ,ONLY: Species, PEM
USE MOD_Particle_Sampling_Vars  ,ONLY: AdaptBCMacroVal, AdaptBCMapElemToSample
USE MOD_DSMC_ElectronicModel    ,ONLY: InitElectronShell
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound
USE MOD_DSMC_Relaxation         ,ONLY: DSMC_SetInternalEnr_Diatomic
! USE MOD_DSMC_PolyAtomicModel    ,ONLY: DSMC_SetInternalEnr_Poly
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: iSpec, iInit, iPart, init_or_sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: TVib                       ! vibrational temperature
REAL                            :: TRot                       ! rotational temperature
INTEGER                         :: ElemID
!===================================================================================================================================
! Nullify energy for atomic species
!-----------------------------------------------------------------------------------------------------------------------------------
PartStateIntEn( 1,iPart) = 0
PartStateIntEn( 2,iPart) = 0
!-----------------------------------------------------------------------------------------------------------------------------------
! Set vibrational and rotational energies for molecules
!-----------------------------------------------------------------------------------------------------------------------------------
IF ((Species(iSpec)%InterID.EQ.2).OR.(Species(iSpec)%InterID.EQ.20)) THEN
  ElemID = PEM%LocalElemID(iPart)
  SELECT CASE (init_or_sf)
  CASE(1) !iInit
    TVib=SpecDSMC(iSpec)%Init(iInit)%TVib
    TRot=SpecDSMC(iSpec)%Init(iInit)%TRot
  CASE(2) !SurfaceFlux
    IF(Species(iSpec)%Surfaceflux(iInit)%Adaptive) THEN
      SELECT CASE(Species(iSpec)%Surfaceflux(iInit)%AdaptiveType)
        CASE(1,3,4) ! Pressure and massflow inlet (pressure/massflow, temperature const)
          TVib=SpecDSMC(iSpec)%Surfaceflux(iInit)%TVib
          TRot=SpecDSMC(iSpec)%Surfaceflux(iInit)%TRot
        CASE(2) ! adaptive Outlet/freestream
          TVib = Species(iSpec)%Surfaceflux(iInit)%AdaptivePressure &
                  / (BoltzmannConst * AdaptBCMacroVal(4,AdaptBCMapElemToSample(ElemID),iSpec))
          TRot = TVib
        CASE DEFAULT
          CALL abort(__STAMP__,'ERROR: Wrong adaptive type for Surfaceflux in DSMC_SetInternalEnr!')
      END SELECT
    ELSE
      TVib=SpecDSMC(iSpec)%Surfaceflux(iInit)%TVib
      TRot=SpecDSMC(iSpec)%Surfaceflux(iInit)%TRot
    END IF
  CASE(3) !reactive surface
    TVib=PartBound%WallTemp(iInit)
    TRot=PartBound%WallTemp(iInit)
  CASE(4) !reactive surface
    TVib=PartBound%WallTemp(iInit)
    TRot=PartBound%WallTemp(iInit)
  CASE DEFAULT
    CALL abort(__STAMP__,'ERROR: Neither iInit nor Surfaceflux defined as reference in DSMC_SetInternalEnr!')
  END SELECT
  ! Background gas distribution
  IF(BGGas%NumberOfSpecies.GT.0) THEN
    IF(BGGas%BackgroundSpecies(iSpec).AND.BGGas%UseDistribution) THEN
      TVib = BGGas%Distribution(BGGas%MapSpecToBGSpec(iSpec),DSMC_TVIB,ElemID)
      TRot = BGGas%Distribution(BGGas%MapSpecToBGSpec(iSpec),DSMC_TROT,ElemID)
    END IF
  END IF
  IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
    CALL DSMC_SetInternalEnr_Poly(iSpec, iInit, iPart, init_or_sf)
  ELSE
    CALL DSMC_SetInternalEnr_Diatomic(iSpec, iPart, TRot, TVib)
  END IF
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! Set electronic energy
!-----------------------------------------------------------------------------------------------------------------------------------
IF (DSMC%ElectronicModel.GT.0) THEN
  IF((Species(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
    CALL InitElectronShell(iSpec,iPart,iInit,init_or_sf)
  ELSE
    PartStateIntEn( 3,iPart) = 0.
  END IF
END IF

END SUBROUTINE DSMC_SetInternalEnr

SUBROUTINE DSMC_SetInternalRotEnr_Poly(iSpec, iPart, TRot)
!===================================================================================================================================
! Initialization of internal rotatinal energy of polyatomic molecules
!> Quantized treatment according to LIECHTY, Derek S. State-to-state internal energy relaxation following the quantum-kinetic
!> model in DSMC. In: 44th AIAA Thermophysics Conference. 2013. S. 2901.; doi: 10.2514/6.2013-2901
!===================================================================================================================================
! MODULES
USE MOD_Globals               ,ONLY: Abort
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst, PlanckConst, PI
USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, SpecDSMC, DSMC,PolyatomMolDSMC,VibQuantsPar,BGGas
USE MOD_Particle_Vars         ,ONLY: PEM, Species
USE MOD_Particle_Sampling_Vars,ONLY: AdaptBCMacroVal, AdaptBCMapElemToSample
USE MOD_DSMC_ElectronicModel  ,ONLY: InitElectronShell
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec,iPart
REAL                          :: TRot
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: iRan, iRan2, NormProb, delta, fNorm!, MomentA, MomentB, MomentC
INTEGER                       :: iQuant, J, kQuant, jMax, iPolyatMole
LOGICAL                       :: ARM

INTEGER                       :: JJMAX,KMAX,delta_max,jIter,kIter
REAL                          :: CurrentValue,MaxValue

!===================================================================================================================================
IF(DSMC%DoRotRelaxQuantized) THEN       ! quantized treatment of rotational energy
  iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
  !//TODO make jmax dependent on energy -> quantum number where energy change is less than rel tol
  ARM = .TRUE.
  IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec)THEN        ! check if molecule is linear
    J = NINT(0.5 * (SQRT(2.*TRot/PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)) - 1.))
    jMax = 4.0 * J    ! thumb rule
    CALL RANDOM_NUMBER(iRan)
    iQuant = INT((1+jMax)*iRan)
    DO WHILE (ARM)
      fNorm = (2.*REAL(iQuant) + 1.)*EXP(-REAL(iQuant)*(REAL(iQuant) + 1.)*PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)/TRot) &
      / ((2.*J + 1.)*EXP(-J*(J + 1.)*PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)/TRot))
      CALL RANDOM_NUMBER(iRan)
      IF(fNorm .LT. iRan) THEN
        CALL RANDOM_NUMBER(iRan)
        iQuant = INT((1+jMax)*iRan)
      ELSE
        ARM = .FALSE.
      END IF
    END DO
    PartStateIntEn( 2,iPart) = REAL(iQuant) * (REAL(iQuant) + 1.) * BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)

  ELSE IF(PolyatomMolDSMC(iPolyatMole)%RotationalGroup.EQ.1)THEN       ! molecule is non-linear -> spherical top molecule with all moments of inertia are the same
    J = NINT(0.5 * (SQRT(4.*TRot/PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)) - 1.))
    jMax = 3.0 * J    ! thumb rule
    CALL RANDOM_NUMBER(iRan)
    iQuant = INT((1+jMax)*iRan)
    DO WHILE (ARM)
      fNorm = (2.*REAL(iQuant) + 1.)**2 *EXP(-REAL(iQuant)*(REAL(iQuant) + 1.)*PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)/TRot) &
      / ((2.*J + 1.)**2 *EXP(-J*(J + 1.)*PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)/TRot))
      CALL RANDOM_NUMBER(iRan)
      IF(fNorm .LT. iRan) THEN
        CALL RANDOM_NUMBER(iRan)
        iQuant = INT((1+jMax)*iRan)
      ELSE
        ARM = .FALSE.
      END IF
    END DO
    PartStateIntEn( 2,iPart) = REAL(iQuant) * (REAL(iQuant) + 1.) * BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)

  ELSE IF(PolyatomMolDSMC(iPolyatMole)%RotationalGroup.EQ.2)THEN      ! molecule is non-linear -> symmetric top molecule where two are equal and third is different
    !//TODO MH algorithm
    ! CALL DSMC_RotRelaxSymTopMH(iPart, iPolyatMole)

    ! PRINT *,'sym top start',TRot
    J = NINT(0.5 * (SQRT(16*BoltzmannConst*TRot*PI**2*PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)/PlanckConst**2) - 1.)) ! for max is K = 0
    ! jMax = 3.0 * J    ! thumb rule
    jMax = 250

    !=================================================================================================
    ! Find max value numerically
    ! jIter = 0
    ! JJMAX = 0
    ! KMAX = 0
    ! MaxValue = 0.
    ! DO WHILE (jIter .LE. jMax)
    !     kIter = 0
    !     DO WHILE (kIter .LE. jIter)
    !         IF (kIter.EQ.0) THEN
    !             delta = 1
    !         ELSE
    !             delta = 0
    !         END IF

    !         CurrentValue = (2. - REAL(delta)) * (2. * REAL(jIter) + 1.) * EXP(-PlanckConst**2 / (8 * PI**2 * BoltzmannConst * TRot) * (REAL(jIter) * (REAL(jIter) + 1) / PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) + (1. / PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) - 1. / PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * REAL(kIter)**2.))

    !         IF (CurrentValue .GT. MaxValue) THEN
    !             MaxValue = CurrentValue
    !             JJMAX = jIter
    !             KMAX = kIter
    !         END IF
    !         kIter = kIter + 1
    !     END DO
    !     jIter = jIter + 1
    ! END DO

    ! print out values for Comparison
    ! PRINT *, "JMAX","Kmax",JJMAX, KMAX
    ! PRINT *, "analytic J", J
    ! PRINT *, "JJ", (2.*(2.*J + 1.)*EXP(-PlanckConst**2 / (8*PI**2*BoltzmannConst*TRot) * (J*(J+1) / PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) + (1/PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) - 1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * J**2)))
    ! PRINT *, "num max", ((2.-REAL(delta_max))*(2.*REAL(JJMAX) + 1.)*EXP(-PlanckConst**2 / (8*PI**2*BoltzmannConst*TRot) * (REAL(JJMAX)*(REAL(JJMAX)+1.) / PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) + (1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) - 1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * REAL(KMAX)**2.)))

    !=================================================================================================
    CALL RANDOM_NUMBER(iRan)
    iQuant = INT((1+jMax)*iRan)
    CALL RANDOM_NUMBER(iRan)
    kQuant = INT((2. * jMax + 1.) * iRan) - jMax

    DO WHILE(kQuant**2.GT.iQuant**2)
      CALL RANDOM_NUMBER(iRan)
      iQuant = INT((1+jMax)*iRan)
      CALL RANDOM_NUMBER(iRan)
      kQuant = INT((2. * jMax + 1.) * iRan) - jMax
    END DO

    !=================================================================================================
    DO WHILE (ARM)
      IF(kQuant.EQ.0)THEN   ! set delta for degeneracy in fNorm
        delta=1.
      ELSE
        delta=0.
      END IF

      fNorm = ((2-delta)*(2.*REAL(iQuant) + 1.)*EXP(-PlanckConst**2 / (8*PI**2*BoltzmannConst*TRot) * (REAL(iQuant)*(REAL(iQuant)+1) / PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) + (1/PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) - 1/PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * REAL(kQuant)**2))) &
      / (2.*(2.*J + 1.)*EXP(-PlanckConst**2 / (8*PI**2*BoltzmannConst*TRot) * (J*(J+1) / PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) + (1/PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) - 1/PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * J**2)))

      ! fNorm = ((2.-REAL(delta))*(2.*REAL(iQuant) + 1.)*EXP(-PlanckConst**2 / (8*PI**2*BoltzmannConst*TRot) * (REAL(iQuant)*(REAL(iQuant)+1) / PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) + (1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) - 1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * REAL(kQuant)**2.))) &
      ! / ((2.*J + 1.)*EXP(-PlanckConst**2 / (8*PI**2*BoltzmannConst*TRot) * (J*(J+1.) / PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1))))

      ! fNorm = ((2.-delta)*(2.*REAL(iQuant) + 1.)*EXP(-PlanckConst**2 / (8*PI**2*BoltzmannConst*TRot) * (REAL(iQuant)*(REAL(iQuant)+1.) / PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) + (1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) - 1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * REAL(kQuant)**2.))) &
      ! / MaxValue
      ! PRINT *,'fNorm',fNorm

      ! sanity checks
      IF(fNorm.GT.1.)THEN
        ! PRINT *, "--------------------------------------------------------------------------------------"
      END IF
      IF(kQuant.GT.iQuant)THEN
        CALL abort(&
          __STAMP__&
          ,'Quantum number k (quantum number of the rotational angular momentum component along the unique axis)&
          is bigger than j!', iSpec)
      END IF

      CALL RANDOM_NUMBER(iRan)
      IF(fNorm .LT. iRan) THEN
        CALL RANDOM_NUMBER(iRan)
        iQuant = INT((1+jMax)*iRan)
        CALL RANDOM_NUMBER(iRan)
        kQuant = INT((2. * jMax + 1.) * iRan) - jMax

        DO WHILE(kQuant**2.GT.iQuant**2)
          CALL RANDOM_NUMBER(iRan)
          iQuant = INT((1+jMax)*iRan)
          CALL RANDOM_NUMBER(iRan)
          kQuant = INT((2. * jMax + 1.) * iRan) - jMax
        END DO
      ELSE
        ARM = .FALSE.
      END IF
    END DO

    PartStateIntEn( 2,iPart) = PlanckConst**2. / (8.*PI**2.) * (REAL(iQuant)*(REAL(iQuant)+1.) / PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) + (1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) - 1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * REAL(kQuant)**2.)

  ELSE IF(PolyatomMolDSMC(iPolyatMole)%RotationalGroup.EQ.3)THEN      ! -> asymmetric top molecule, no analytic formula for energy levels
    CALL abort(&
    __STAMP__&
    ,'Quantized treatment of asymmetric top molecules not possible yet!')
    ! //TODO:maybe special case for sightly asymmetrical molecules, e.g. H2O
  ELSE
    CALL abort(&
    __STAMP__&
    ,'Unexpected dimensions of moments of inertia of species iSpec!')
  END IF
ELSE
  ! continous treatment of rotational energy
  IF (SpecDSMC(iSpec)%Xi_Rot.EQ.2) THEN
    CALL RANDOM_NUMBER(iRan)
    PartStateIntEn( 2,iPart) = -BoltzmannConst*TRot*LOG(iRan)
  ELSE IF (SpecDSMC(iSpec)%Xi_Rot.EQ.3) THEN
    CALL RANDOM_NUMBER(iRan)
    PartStateIntEn( 2,iPart) = iRan*10 !the distribution function has only non-negligible  values betwenn 0 and 10
    NormProb = SQRT(PartStateIntEn( 2,iPart))*EXP(-PartStateIntEn( 2,iPart))/(SQRT(0.5)*EXP(-0.5))
    CALL RANDOM_NUMBER(iRan)
    DO WHILE (iRan.GE.NormProb)
      CALL RANDOM_NUMBER(iRan)
      PartStateIntEn( 2,iPart) = iRan*10 !the distribution function has only non-negligible  values betwenn 0 and 10
      NormProb = SQRT(PartStateIntEn( 2,iPart))*EXP(-PartStateIntEn( 2,iPart))/(SQRT(0.5)*EXP(-0.5))
      CALL RANDOM_NUMBER(iRan)
    END DO
    PartStateIntEn( 2,iPart) = PartStateIntEn( 2,iPart)*BoltzmannConst*TRot
  END IF
END IF

END SUBROUTINE DSMC_SetInternalRotEnr_Poly

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
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound
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
    CASE(3) !reactive surface
      TVib=PartBound%WallTemp(iInit)
      TRot=PartBound%WallTemp(iInit)
    CASE(4) !reactive surface
      TVib=PartBound%WallTemp(iInit)
      TRot=PartBound%WallTemp(iInit)
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

! Set vibrational energy
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
! set initial rotational internal energy
CALL DSMC_SetInternalRotEnr_Poly(iSpecies, iPart, TRot)

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
  ! set initial rotational internal energy
  CALL DSMC_SetInternalRotEnr_Poly(iSpec, iPart, TRot)
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

    ! set initial rotational internal energy
    CALL DSMC_SetInternalRotEnr_Poly(iSpec, iPart, TRot)
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

  ! set initial rotational internal energy
  CALL DSMC_SetInternalRotEnr_Poly(iSpec, iPart, TRot)
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
!//TODO also add to other routine????
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
  USE MOD_Globals               ,ONLY: Abort
  USE MOD_Globals_Vars          ,ONLY: BoltzmannConst, PlanckConst, PI
  USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, Coll_pData, SpecDSMC, DSMC,PolyatomMolDSMC
  USE MOD_Particle_Vars         ,ONLY: PartSpecies
  USE, INTRINSIC                :: ieee_arithmetic

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
  REAL                          :: iRan, NormProb, fNorm, tempProb, fak1, fak2, Ec, delta, k!, MomentA, MomentB, MomentC
  INTEGER                       :: iQuant, J1, J2, JStar, kQuant, iPolyatMole, iSpec
  LOGICAL                       :: ARM
!===================================================================================================================================
!//TODO: Function pointer for continous and quant energy
!//TODO: check for other collision routines -> e.g. prohibitted double relaxation
  IF(DSMC%DoRotRelaxQuantized) THEN       ! quantized treatment of rotational energy
    iSpec = PartSpecies(iPart)
    iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
    ARM = .TRUE.
    IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec)THEN        ! check if molecule is linear, same as diatomic
      J2 = INT((-1.+SQRT(1.+(4.*Coll_pData(iPair)%Ec)/(BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))))/2.)
      !//TODO check J1 here and for diatomic again!!
      ! my analytic solution
      J1 = NINT(0.5 * (-1. + SQRT((1+ 4 * Coll_pData(iPair)%Ec / (BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)))/(3. - 2. * SpecDSMC(iSpec)%omega))))
      ! PRINT *,"analytic", J1
      ! wolfram alpha solution (same as analytic)
      ! J1 = NINT(0.5 * (-1. + ((4.*Coll_pData(iPair)%Ec + BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)))/SQRT((3. - 2. * SpecDSMC(iSpec)%omega)*BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)*(4.*Coll_pData(iPair)%Ec + BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)))))
      ! PRINT *,"wolfram", J1
      ! paper solution
      ! J1 = NINT(0.5 * (-1. + SQRT((2. + 3. * SpecDSMC(iSpec)%omega + (8. * Coll_pData(iPair)%Ec) / &
          ! (BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)))/(6. - SpecDSMC(iSpec)%omega))))
      ! J2 much bigger than J1 often
      ! PRINT *, "J1",J1
      ! PRINT *, "J2",J2
      JStar = MIN(J1,J2)

      CALL RANDOM_NUMBER(iRan)
      iQuant = INT((1+J2)*iRan)
      DO WHILE (ARM)
        fNorm = (2.*REAL(iQuant) + 1.)*(Coll_pData(iPair)%Ec - REAL(iQuant)*(REAL(iQuant) + 1.)*BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))**FakXi &
                / ((2.*REAL(JStar) + 1.)*(Coll_pData(iPair)%Ec - REAL(JStar)*(REAL(JStar) + 1.)*BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))**FakXi)
        CALL RANDOM_NUMBER(iRan)
        IF(fNorm .LT. iRan) THEN
          CALL RANDOM_NUMBER(iRan)
          iQuant = INT((1+J2)*iRan)
        ELSE
          ARM = .FALSE.
        END IF
      END DO
      PartStateIntEn( 2,iPart) = REAL(iQuant) * (REAL(iQuant) + 1.) * BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)

    ELSE IF(PolyatomMolDSMC(iPolyatMole)%RotationalGroup.EQ.1)THEN       ! molecule is non-linear -> spherical top molecule with all moments of inertia are the same
      J2 = INT((-1.+SQRT(1.+(4.*Coll_pData(iPair)%Ec)/(BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))))/2.) ! same as linear - only different degeneracy
      J1 = NINT(0.5 * (-1. + SQRT((3.*SpecDSMC(iSpec)%omega-2.+(4.*Coll_pData(iPair)%Ec)/(BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)))/(2. - SpecDSMC(iSpec)%omega))))
      JStar = MIN(J1,J2)

      CALL RANDOM_NUMBER(iRan)
      iQuant = INT((1+J2)*iRan)
      DO WHILE (ARM)
        fNorm = (REAL(iQuant)**2. + REAL(iQuant) + 1.) *(Coll_pData(iPair)%Ec - REAL(iQuant)*(REAL(iQuant) + 1.)*BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))**FakXi &
          / ( (REAL(JStar)**2. + REAL(JStar) + 1.) *(Coll_pData(iPair)%Ec - REAL(JStar)*(REAL(JStar) + 1.)*BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))**FakXi )
        CALL RANDOM_NUMBER(iRan)
        IF(fNorm .LT. iRan) THEN
          CALL RANDOM_NUMBER(iRan)
          iQuant = INT((1+J2)*iRan)
        ELSE
          ARM = .FALSE.
        END IF
      END DO
      PartStateIntEn( 2,iPart) = REAL(iQuant) * (REAL(iQuant) + 1.) * BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)

    ELSE IF(PolyatomMolDSMC(iPolyatMole)%RotationalGroup.EQ.2)THEN      ! molecule is non-linear -> symmetric top molecule where two are equal and third is different
      CALL DSMC_RotRelaxSymTopARM(iPair, iPart,FakXi, iPolyatMole)
      ! CALL DSMC_RotRelaxSymTopMH(iPair, iPart,FakXi, iPolyatMole)
      !CALL DSMC_RotRelaxSymTopGibbsSampling(iPair, iPart,FakXi, iPolyatMole)

    ELSE IF(PolyatomMolDSMC(iPolyatMole)%RotationalGroup.EQ.3)THEN      ! -> asymmetric top molecule, no analytic formula for energy levels
      ! //TODO:maybe special case for sightly asymmetrical molecules, e.g. H2O 10.1139/p57-096
    ELSE
      CALL abort(&
        __STAMP__&
        ,'Unexpected relations of moments of inertia of species!', iSpec)
    END IF
  ELSE  ! continous treatment of rotational energy
    Ec = Coll_pData(iPair)%Ec
    fak1 = (3.0/2.0+FakXi-1.0)/(3.0/2.0-1.0)
    fak2 = (3.0/2.0+FakXi-1.0)/(FakXi)

    CALL RANDOM_NUMBER(iRan)
    tempProb = Ec*iRan
    NormProb = ((fak1*tempProb/Ec)**(3.0/2.0-1.0))*((fak2*(1.0-tempProb/Ec))**(FakXi))
    CALL RANDOM_NUMBER(iRan)
    DO WHILE (iRan.GE.NormProb)
      CALL RANDOM_NUMBER(iRan)
      tempProb = Ec*iRan
      NormProb = (fak1*tempProb/Ec)**(3.0/2.0-1.0)*(fak2*(1.0-tempProb/Ec))**(FakXi)
      CALL RANDOM_NUMBER(iRan)
    END DO
    PartStateIntEn(2,iPart)=tempProb
  END IF

END SUBROUTINE DSMC_RotRelaxPoly

SUBROUTINE DSMC_RotRelaxSymTopARM(iPair, iPart,FakXi, iPolyatMole)
!===================================================================================================================================
! Rotational relaxation routine for symmetric top molecules using the acceptance rejection sampling
!===================================================================================================================================
! MODULES
  USE MOD_Globals               ,ONLY: Abort
  USE MOD_Globals_Vars          ,ONLY: BoltzmannConst, PlanckConst, PI
  USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, Coll_pData, SpecDSMC, DSMC,PolyatomMolDSMC
  USE MOD_Particle_Vars         ,ONLY: PartSpecies

! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iPart, iPair, iPolyatMole
  REAL, INTENT(IN)              :: FakXi
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                          :: iRan, NormProb, fNorm, tempProb, fak1, fak2, Ec, delta, k
  INTEGER                       :: iQuant, J1, J2, JStar, kQuant, iSpec
  LOGICAL                       :: ARM
!===================================================================================================================================
! PRINT *,'sym top routine start'
  ARM = .TRUE.
  !//TODO 16.05 maybe arm for two quantum numbers
  ! k = 0
  ! J2 = INT((-1.+SQRT(1.+(32.*Coll_pData(iPair)%Ec*PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)*PI**2.)/(PlanckConst**2.)))/2.)
  ! PRINT *,"k0", J2

  ! PRINT *, PlanckConst**2. / (8.*PI**2.) * (REAL(J2)*(REAL(J2)+1.) / PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) + (1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) - 1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * REAL(0.)**2)

  ! k = J
  J2 = INT((-PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3)/PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) +SQRT((PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3)/PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1))**2. +(32.*Coll_pData(iPair)%Ec*PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3)*PI**2.)/(PlanckConst**2.)))/2.)
  ! PRINT *,"kJ", J2

  ! PRINT *, PlanckConst**2. / (8.*PI**2.) * (REAL(J2)*(REAL(J2)+1.) / PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) + (1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) - 1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * REAL(J2)**2)
  ! PRINT *,"Ecoll", Coll_pData(iPair)%Ec

  ! wolfram solution
  ! J1 = NINT((-1.+SQRT((FakXi-1)*(32.*PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)*Coll_pData(iPair)%Ec*PI**2.+PlanckConst**2.))/(FakXi-1.)*PlanckConst)/2.)
  ! PRINT *,"wolfram", J1

  ! my analytic solution, k=0
  J1 = NINT(0.5 * (-1. + SQRT((1.+ 32. * Coll_pData(iPair)%Ec * PI**2. * PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) / (PlanckConst**2.))/(3. - 2. * SpecDSMC(iSpec)%omega))))
  ! PRINT *,"analytic J",J1

  JStar = MIN(J1,J2)
  ! PRINT *,'J1,J2,JStar',J1,J2,JStar

  CALL RANDOM_NUMBER(iRan)
  iQuant = INT((1+J2)*iRan)
  IF(iQuant.EQ.0)THEN
    kQuant = 0
    delta = 1.
  ELSE
    CALL RANDOM_NUMBER(iRan)
    kQuant = INT((2. * J2 + 1.) * iRan) - J2
    delta = 0.
  END IF

  DO WHILE(kQuant**2.GT.iQuant**2)
    CALL RANDOM_NUMBER(iRan)
    iQuant = INT((1+J2)*iRan)
    IF(iQuant.EQ.0)THEN
      kQuant = 0
      delta = 1.
    ELSE
      CALL RANDOM_NUMBER(iRan)
      kQuant = INT((2. * J2 + 1.) * iRan) - J2
      delta = 0.
    END IF
  END DO
  DO WHILE (ARM)
    IF((Coll_pData(iPair)%Ec - (PlanckConst**2.)/(PI**2.*8.) * &
      ((REAL(iQuant)*(REAL(iQuant)+1))/(PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) + &
      (1/PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) - 1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * REAL(kQuant)**2.)).LT.0)THEN
      fNorm = -1! set fNorm to less than than 1 to redo loop since the quantum numbers were not possible due to the collision energy
    ELSE
      fNorm = ((2-delta) * (2.*REAL(iQuant)+1.) * (Coll_pData(iPair)%Ec - (PlanckConst**2.)/(PI**2.*8.) * &
      ((REAL(iQuant)*(REAL(iQuant)+1.))/(PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) + &
      (1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) - 1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * REAL(kQuant)**2.))**FakXi / &
      (2. * (2.*REAL(JStar)+1.) * (Coll_pData(iPair)%Ec - (PlanckConst**2.)/(PI**2.*8.) * ((REAL(JStar)*(REAL(JStar)+1.))/(PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) + &
      (1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) - 1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * REAL(JStar)**2.))**FakXi))
      ! PRINT *,'fNorm',fNorm
    END IF

    CALL RANDOM_NUMBER(iRan)
    IF(fNorm .LT. iRan) THEN
      CALL RANDOM_NUMBER(iRan)
      iQuant = INT((1+J2)*iRan)
      IF(iQuant.EQ.0)THEN
        kQuant = 0
        delta = 1.
      ELSE
        CALL RANDOM_NUMBER(iRan)
        kQuant = INT((2. * J2 + 1.) * iRan) - J2
        delta = 0.
      END IF

      DO WHILE(kQuant**2.GT.iQuant**2)
        CALL RANDOM_NUMBER(iRan)
        iQuant = INT((1+J2)*iRan)
        IF(iQuant.EQ.0)THEN
          kQuant = 0
          delta = 1.
        ELSE
          CALL RANDOM_NUMBER(iRan)
          kQuant = INT((2. * J2 + 1.) * iRan) - J2
          delta = 0.
        END IF
      END DO
    ELSE
      ARM = .FALSE.
    END IF
  END DO

  ! sanity checks -> //TODO: J2 vs JStar
  IF(fNorm.GT.1.)THEN
    ! PRINT *, "--------------------------------------------------------------------------------------"
  END IF
  ! sanity check
  IF(kQuant.GT.iQuant)THEN
    CALL abort(&
      __STAMP__&
      ,'Quantum number k (quantum number of the rotational angular momentum component along the unique axis)&
      is bigger than j!', iSpec)
  END IF

  PartStateIntEn( 2,iPart) = PlanckConst**2. / (8.*PI**2.) * (REAL(iQuant)*(REAL(iQuant)+1.) / PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) + (1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) - 1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * REAL(kQuant)**2)

END SUBROUTINE DSMC_RotRelaxSymTopARM

SUBROUTINE DSMC_RotRelaxSymTopMH(iPart, iPolyatMole)
!===================================================================================================================================
! //TODO
!===================================================================================================================================
! MODULES
  USE MOD_Globals               ,ONLY: Abort
  USE MOD_Globals_Vars          ,ONLY: BoltzmannConst, PlanckConst, PI
  USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, Coll_pData, SpecDSMC, DSMC,PolyatomMolDSMC
  USE MOD_Particle_Vars         ,ONLY: PartSpecies

! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iPart, iPolyatMole
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                          :: iRan2, NormProb, fNorm, tempProb, fak1, fak2, Ec, delta, k!, MomentA, MomentB, MomentC
  INTEGER                       :: J1, J2, JStar, kQuant, iSpec, jMax
  LOGICAL                       :: ARM
  INTEGER,ALLOCATABLE           :: iQuant(:), iQuant_old(:)
  REAL,ALLOCATABLE              :: iRan(:)


! LOCAL VARIABLES
  INTEGER                       :: iDOF, iWalk
  REAL                          :: TVib                       ! vibrational temperature
  REAL                          :: TRot                       ! rotational temperature
  INTEGER                       :: ElemID
!===================================================================================================================================

  ALLOCATE(iRan(2), iQuant(2), iQuant_old(2))
  jMax = 250
  DO iWalk = 1, 4000
    iQuant_old(:)=iQuant(:)
    CALL RANDOM_NUMBER(iRan)
    iQuant(:) = INT(iRan(:)*jMax)
    DO WHILE (iQuant(2)**2.GT.iQuant(1)**2)
      CALL RANDOM_NUMBER(iRan)
      iQuant(:) = INT(iRan(:)*jMax)
    END DO

    NormProb = NormProb + PlanckConst**2. / (8.*PI**2.) * (REAL((iQuant_old(1)-iQuant(1)))*(REAL((iQuant_old(1)-iQuant(1)))+1.) / PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) + (1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) - 1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * REAL((iQuant_old(2)-iQuant(2)))**2)

    NormProb = MIN(1.0,EXP(NormProb))

    CALL RANDOM_NUMBER(iRan2)
    IF (NormProb.LT.iRan2) iQuant(:)=iQuant_old(:)
  END DO

  PartStateIntEn( 2,iPart) = PlanckConst**2. / (8.*PI**2.) * (REAL(iQuant(1))*(REAL(iQuant(1))+1.) / PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) + (1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) - 1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * REAL(iQuant(2))**2)
  DEALLOCATE(iRan, iQuant,iQuant_old)

END SUBROUTINE DSMC_RotRelaxSymTopMH

END MODULE MOD_DSMC_PolyAtomicModel