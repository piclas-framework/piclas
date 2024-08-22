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

ABSTRACT INTERFACE
  SUBROUTINE RotRelaxPolyRoutine(iPair,iPart,FakXi)
    INTEGER,INTENT(IN)          :: iPair, iPart               ! index of collision pair
    REAL,INTENT(IN)             :: FakXi
  END SUBROUTINE
END INTERFACE

PROCEDURE(RotRelaxPolyRoutine),POINTER :: RotRelaxPolyRoutineFuncPTR !< pointer defining the function called for rotational relaxation
                                                                !  depending on the RotRelaxModel (continous or quantized)
                                                                !  for polyatomic molecules

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: InitPolyAtomicMolecs, DSMC_SetInternalEnr_Poly_ARM, DSMC_SetInternalEnr_Poly_MH, DSMC_SetInternalEnr_Poly_MH_FirstPick
PUBLIC :: DSMC_RotRelaxPoly, DSMC_VibRelaxPoly_ARM, DSMC_VibRelaxPoly_MH, DSMC_VibRelaxPoly_ARM_MH
PUBLIC :: DSMC_FindFirstVibPick, DSMC_RelaxVibPolyProduct, RotRelaxPolyRoutineFuncPTR, DSMC_SetInternalEnr
!===================================================================================================================================

CONTAINS

SUBROUTINE InitPolyAtomicMolecs(iSpec)
!===================================================================================================================================
!> Initialization of variables for polyatomic molecules
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst, PlanckConst, PI
USE MOD_DSMC_Vars             ,ONLY: DSMC, SpecDSMC, PolyatomMolDSMC
USE MOD_DSMC_ElectronicModel  ,ONLY: ReadRotationalSpeciesLevel
USE MOD_ReadInTools
USE MOD_PARTICLE_Vars         ,ONLY: Species, SpeciesDatabase
USE MOD_io_hdf5
USE MOD_HDF5_input            ,ONLY: ReadAttribute, DatasetExists, AttributeExists
#if USE_MPI
USE MOD_LoadBalance_Vars      ,ONLY: PerformLoadBalance
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
  CALL ReadAttribute(file_id_specdb,'LinearMolec',1,DatasetName = dsetname,IntScalar=IntToLog,ChangeToGroup=.True.)
  IF(IntToLog.EQ.1) THEN
    PolyatomMolDSMC(iPolyatMole)%LinearMolec = .TRUE.
  ELSE
    PolyatomMolDSMC(iPolyatMole)%LinearMolec = .FALSE.
  END IF
  CALL PrintOption('LinearMolec, '//TRIM(Species(iSpec)%Name),'DB',LogOpt=PolyatomMolDSMC(iPolyatMole)%LinearMolec)
  ! Number of atoms
  CALL ReadAttribute(file_id_specdb,'NumOfAtoms',1,DatasetName = dsetname,  &
    IntScalar=PolyatomMolDSMC(iPolyatMole)%NumOfAtoms,ChangeToGroup=.True.)
  CALL PrintOption('NumOfAtoms, '//TRIM(Species(iSpec)%Name),'DB',IntOpt=PolyatomMolDSMC(iPolyatMole)%NumOfAtoms)
  ! Dissociation energy
  ! TSHO not implemented with polyatomic molecules, but Ediss_eV required for the calculation of polyatomic temp. (upper bound)
  CALL ReadAttribute(file_id_specdb,'Ediss_eV',1,DatasetName = dsetname,RealScalar=SpecDSMC(iSpec)%Ediss_eV,ChangeToGroup=.True.)
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
    IF(DSMC%RotRelaxModel.EQ.1)THEN
      CALL AttributeExists(file_id_specdb,'MomentOfInertia',TRIM(dsetname), AttrExists=AttrExists,ChangeToGroup=.True.)
      IF (AttrExists) THEN
        CALL ReadAttribute(file_id_specdb,'MomentOfInertia',1,DatasetName = dsetname, &
          RealScalar=PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1),ChangeToGroup=.True.)
        PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1) = PlanckConst**2 / &
          (8 * PI**2 * PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) * BoltzmannConst)
        CALL PrintOption('MomentOfInertia, '//TRIM(Species(iSpec)%Name),'DB', &
          RealOpt=PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1))
      END IF
    ELSE  ! DSMC RotRelaxModel NE 1
      CALL AttributeExists(file_id_specdb,'CharaTempRot',TRIM(dsetname), AttrExists=AttrExists,ChangeToGroup=.True.)
      IF(AttrExists)THEN
        CALL ReadAttribute(file_id_specdb,'CharaTempRot',1,DatasetName = dsetname,  &
          RealScalar=PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1),ChangeToGroup=.True.)
      ELSE  ! CharaTempRot not found
        CALL AttributeExists(file_id_specdb,'MomentOfInertia',TRIM(dsetname), AttrExists=AttrExists,ChangeToGroup=.True.)
        IF (AttrExists) THEN
          CALL ReadAttribute(file_id_specdb,'MomentOfInertia',1,DatasetName = dsetname, &
            RealScalar=PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1),ChangeToGroup=.True.)
          PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1) = PlanckConst**2 / &
            (8 * PI**2 * PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) * BoltzmannConst)
        ELSE
          PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1) = 0
        END IF
      END IF
    END IF
    CALL PrintOption('CharaTempRot, '//TRIM(Species(iSpec)%Name),'DB',RealOpt=PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))
    ! set dummy values
    PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(2:3) = 0
    PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(2:3) = 0
  ELSE
    IF(DSMC%RotRelaxModel.EQ.1)THEN
      CALL AttributeExists(file_id_specdb,'MomentOfInertia',TRIM(dsetname), AttrExists=AttrExists,ChangeToGroup=.True.)
      IF (AttrExists) THEN
        CALL ReadAttribute(file_id_specdb,'MomentOfInertia',3,DatasetName = dsetname,  &
          RealArray=PolyatomMolDSMC(iPolyatMole)%MomentOfInertia,ChangeToGroup=.True.)
        DO iVibDOF = 1,3
          WRITE(UNIT=hilf2,FMT='(I0)') iVibDOF
          CALL PrintOption('MomentOfInertia'//TRIM(hilf2)//' '//TRIM(Species(iSpec)%Name),'DB', &
            RealOpt=PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(iVibDOF))
          PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(iVibDOF) = PlanckConst**2 / &
            (8 * PI**2 * PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(iVibDOF) * BoltzmannConst)
        END DO
      ELSE
        CALL abort(&
        __STAMP__&
        ,'Moment of inertia necessary for quantized rotational energy and is not set for species', iSpec)
      END IF
      ! sanity checks for order of moments of inertia
      IF((PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3).EQ.PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(2).AND.  &
        PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3).NE.PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)))THEN
        ! C should be different one (for symmetric tops) - for sperical tops all should be the same
        CALL abort(&
        __STAMP__&
        ,'Moments of inertia in wrong order in database for species', iSpec)
      END IF
    ELSE  ! DSMC RotRelaxModel NE 1
      DO iVibDOF = 1,3
        WRITE(UNIT=hilf2,FMT='(I0)') iVibDOF
        CALL AttributeExists(file_id_specdb,TRIM('CharaTempRot'//TRIM(hilf2)),TRIM(dsetname), &
          AttrExists=AttrExists,ChangeToGroup=.True.)
        IF (AttrExists) THEN  ! read in CharaTempRot
          CALL ReadAttribute(file_id_specdb,TRIM('CharaTempRot'//TRIM(hilf2)),1,DatasetName = dsetname, &
            RealScalar=PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(iVibDOF),ChangeToGroup=.True.)
          CALL PrintOption('CharaTempRot'//TRIM(hilf2)//' '//TRIM(Species(iSpec)%Name),'DB', &
            RealOpt=PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(iVibDOF))
        ELSE  ! check for MomentOfInertia
          CALL AttributeExists(file_id_specdb,'MomentOfInertia',TRIM(dsetname), AttrExists=AttrExists,ChangeToGroup=.True.)
          IF (AttrExists) THEN
            CALL ReadAttribute(file_id_specdb,'MomentOfInertia',3,DatasetName = dsetname,  &
              RealArray=PolyatomMolDSMC(iPolyatMole)%MomentOfInertia,ChangeToGroup=.True.)
            WRITE(UNIT=hilf2,FMT='(I0)') iVibDOF
            PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(iVibDOF) = PlanckConst**2 / &
              (8 * PI**2 * PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(iVibDOF) * BoltzmannConst)
          ELSE ! set dummy value to zero
            PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(iVibDOF) = 0
          END IF
        END IF
      END DO
    END IF
    DO iVibDOF = 1,3
      WRITE(UNIT=hilf2,FMT='(I0)') iVibDOF
      CALL PrintOption('CharaTempRot'//TRIM(hilf2)//' '//TRIM(Species(iSpec)%Name),'DB', &
        RealOpt=PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(iVibDOF))
    END DO
  END IF
  ! Read-in of characteristic vibrational temperature
  DO iVibDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
    WRITE(UNIT=hilf2,FMT='(I0)') iVibDOF
    CALL ReadAttribute(file_id_specdb,TRIM('CharaTempVib'//TRIM(hilf2)),1,DatasetName = dsetname, &
      RealScalar=PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iVibDOF),ChangeToGroup=.True.)
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
    IF(DSMC%RotRelaxModel.EQ.1)THEN
      PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)   = GETREAL('Part-Species'//TRIM(hilf)//'-MomentOfInertia')
      PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(2:3) = 0
      PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)      = PlanckConst**2 / (8 * PI**2 *   &
        PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) * BoltzmannConst)
    ELSE  ! DSMC RotRelaxModel NE 1
      PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)                   = GETREAL('Part-Species'//TRIM(hilf)//'-CharaTempRot')
    END IF
    PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(2:3)    = 0
  ELSE
    IF(DSMC%RotRelaxModel.EQ.1)THEN
      DO iVibDOF = 1,3
        WRITE(UNIT=hilf2,FMT='(I0)') iVibDOF
        PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(iVibDOF) = GETREAL('Part-Species'//TRIM(hilf)//'-MomentOfInertia'//TRIM(hilf2))
        PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(iVibDOF)    = PlanckConst**2 / (8 * PI**2 *   &
          PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(iVibDOF) * BoltzmannConst)
      END DO
      ! sanity checks for order of moments of inertia
      IF((PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3).EQ.PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(2).AND.  &
        PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3).NE.PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)))THEN
        ! C should be different one (for symmetric tops) - for sperical tops all should be the same
        CALL abort(&
        __STAMP__&
        ,'Moments of inertia in wrong order in database for species', iSpec)
      END IF
    ELSE
      DO iVibDOF = 1,3
        WRITE(UNIT=hilf2,FMT='(I0)') iVibDOF
        PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(iVibDOF) = GETREAL('Part-Species'//TRIM(hilf)//'-CharaTempRot'//TRIM(hilf2))
      END DO
    END IF
  END IF
  ! Read-in of characteristic vibrational temperature
  DO iVibDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
    WRITE(UNIT=hilf2,FMT='(I0)') iVibDOF
    PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iVibDOF) = GETREAL('Part-Species'//TRIM(hilf)//'-CharaTempVib'//TRIM(hilf2))
  END DO
END IF

! Calculation of zero-point energy
DO iVibDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
  SpecDSMC(iSpec)%EZeroPoint = SpecDSMC(iSpec)%EZeroPoint + DSMC%GammaQuant*  &
    PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iVibDOF)*BoltzmannConst
END DO

! save Rotational Group of molecule with convention I_A .LE. I_B .LE. I_C
IF(.NOT.PolyatomMolDSMC(iPolyatMole)%LinearMolec)THEN
  IF(PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1).EQ.PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(2).AND. &
    PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(2).EQ.PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3))THEN
    ! sphrical top with A=B=C
    PolyatomMolDSMC(iPolyatMole)%RotationalGroup = 1
  ELSE IF(PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1).EQ.PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(2).AND. &
    PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3).GT.PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1))THEN
    ! oblate symmetric top with A=B,
    PolyatomMolDSMC(iPolyatMole)%RotationalGroup = 10
  ELSE IF(PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1).LT.PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(2).AND. &
    PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3).EQ.PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1))THEN
    ! prolate symmetric top with A=B,
    PolyatomMolDSMC(iPolyatMole)%RotationalGroup = 11
  ELSE IF(PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1).NE.PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(2).AND. &
    PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(2).NE.PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3).AND. &
    PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1).NE.PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3))THEN
    ! asymmetric top with all moments of inertia different
    PolyatomMolDSMC(iPolyatMole)%RotationalGroup = 3
  ELSE ! set dummy to catch false cases
    PolyatomMolDSMC(iPolyatMole)%RotationalGroup = -1
  END IF
END IF
! read in rotational levels for asymmetric tops for RotRelaxModel 1
! read in for all species happens only for RotRelaxModel 2
IF(PolyatomMolDSMC(iPolyatMole)%RotationalGroup.EQ.3.AND.DSMC%RotRelaxModel.EQ.1) CALL ReadRotationalSpeciesLevel(iSpec)
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


SUBROUTINE DSMC_SetInternalEnr_Poly_ARM_SingleMode(iSpecies, iInit, iPart, init_or_sf)
!===================================================================================================================================
! Initialization of polyatomic molecules by treating every mode separately in a loop
!===================================================================================================================================
! MODULES
USE MOD_Globals                ,ONLY: Abort
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars              ,ONLY: PartStateIntEn, SpecDSMC, DSMC,PolyatomMolDSMC,VibQuantsPar,BGGas
USE MOD_Particle_Vars          ,ONLY: PEM, Species
USE MOD_Particle_Sampling_Vars ,ONLY: AdaptBCMacroVal, AdaptBCMapElemToSample
USE MOD_DSMC_ElectronicModel   ,ONLY: InitElectronShell
USE MOD_part_tools             ,ONLY: RotInitPolyRoutineFuncPTR
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
REAL                          :: iRan
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
!//TODO
! set initial rotational internal energy
PartStateIntEn( 2,iPart) = RotInitPolyRoutineFuncPTR(iSpecies,TRot,iPart)

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
USE MOD_part_tools            ,ONLY: RotInitPolyRoutineFuncPTR
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
  PartStateIntEn( 2,iPart) = RotInitPolyRoutineFuncPTR(iSpec,TRot,iPart)
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
USE MOD_part_tools            ,ONLY: RotInitPolyRoutineFuncPTR
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
    PartStateIntEn( 2,iPart) = RotInitPolyRoutineFuncPTR(iSpec,TRot,iPart)
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
USE MOD_part_tools            ,ONLY: RotInitPolyRoutineFuncPTR
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
  PartStateIntEn( 2,iPart) = RotInitPolyRoutineFuncPTR(iSpec,TRot,iPart)
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
! Continous rotational relaxation routine
!===================================================================================================================================
! MODULES
  USE MOD_Globals               ,ONLY: Abort
  USE MOD_Globals_Vars          ,ONLY: BoltzmannConst, PlanckConst, PI
  USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, Coll_pData, SpecDSMC
  USE MOD_Particle_Vars         ,ONLY: PartSpecies

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
  REAL                          :: iRan, NormProb, tempProb, fak1, fak2, Ec, LocalFakXi
  INTEGER                       :: iSpec
!===================================================================================================================================
iSpec = PartSpecies(iPart)
IF(SpecDSMC(iSpec)%Xi_Rot.EQ.3) THEN
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
ELSE    ! continous treatment of linear molecules
  ! fix for changed FakXi for polyatomic
  LocalFakXi = FakXi + 0.5*SpecDSMC(iSpec)%Xi_Rot
  CALL RANDOM_NUMBER(iRan)
  PartStateIntEn(2, iPart) = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/LocalFakXi))
END IF

END SUBROUTINE DSMC_RotRelaxPoly


SUBROUTINE DSMC_RotRelaxQuantPoly(iPair, iPart,FakXi)
!===================================================================================================================================
! Quantized rotational relaxation routine for given collision energy
! Different rotational groups are sampled differently:
! - linear molecules, spherical molecules and symmetric top molecules  - Acceptance-Rejection sampling
! - asymmetric top molecules                                           - only possible with database of rotational levels
!===================================================================================================================================
! MODULES
  USE MOD_Globals               ,ONLY: Abort
  USE MOD_Globals_Vars          ,ONLY: BoltzmannConst, PlanckConst, PI
  USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, Coll_pData, SpecDSMC, PolyatomMolDSMC, RadialWeighting
  USE MOD_Particle_Vars         ,ONLY: PartSpecies, UseVarTimeStep, usevMPF
  USE MOD_part_tools            ,ONLY: GetParticleWeight

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
  REAL                          :: iRan, fNorm, Ec, MaxValue, CurrentValue
  INTEGER                       :: iQuant, kQuant, J2, iPolyatMole, iSpec, jIter, kIter, delta
  LOGICAL                       :: ARM
!===================================================================================================================================
IF (usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep) THEN
  Ec = Coll_pData(iPair)%Ec / GetParticleWeight(iPart)
ELSE
  Ec = Coll_pData(iPair)%Ec
END IF
iSpec = PartSpecies(iPart)
iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
ARM = .TRUE.
IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec)THEN        ! check if molecule is linear, same routine as diatomic
  ! calculate maximum allowed energy (all of collision energy in rotatinal energy)
  J2 = INT((-1.+SQRT(1.+(4.*Ec)/(BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))))/2.)
  ! reduce J2 if too big which would correspond to unphysical quantum numbers, necessary for high velocities since otherwise
  ! almost no samples will be accepted
  IF(J2.GT.500) J2 = 500
  ! Find max value of distribution for ARM numerically
  MaxValue = 0.
  DO jIter=0, J2
    CurrentValue = (2.*REAL(jIter) + 1.)*(Ec - REAL(jIter)*(REAL(jIter) + 1.)*BoltzmannConst*SpecDSMC(iSpec)%CharaTRot)**FakXi
    IF (CurrentValue .GT. MaxValue) MaxValue = CurrentValue
  END DO
  CALL RANDOM_NUMBER(iRan)
  iQuant = INT((1+J2)*iRan)
  DO WHILE (ARM)
    fNorm = (2.*REAL(iQuant)+1.)* &
            (Ec-REAL(iQuant)*(REAL(iQuant)+1.)*BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))**FakXi / MaxValue
    CALL RANDOM_NUMBER(iRan)
    IF(fNorm .LT. iRan) THEN
      CALL RANDOM_NUMBER(iRan)
      iQuant = INT((1+J2)*iRan)
    ELSE
      ARM = .FALSE.
    END IF
  END DO
  PartStateIntEn( 2,iPart) = REAL(iQuant) * (REAL(iQuant) + 1.) * BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)

ELSE IF(PolyatomMolDSMC(iPolyatMole)%RotationalGroup.EQ.1)THEN
  ! molecule is non-linear -> spherical top molecule with all moments of inertia are the same
  ! calculate maximum allowed energy (all of collision energy in rotatinal energy)
  J2 = INT((-1.+SQRT(1.+(4.*Ec)/(BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))))/2.)
  ! reduce J2 if too big which would correspond to unphysical quantum numbers, necessary for high velocities since otherwise
  ! almost no samples will be accepted
  IF(J2.GT.500) J2 = 500
  ! Find max value of distribution for ARM numerically
  MaxValue = 0.
  DO jIter=0, J2
    CurrentValue = (2.*REAL(jIter) + 1.)**2 *(Ec - REAL(jIter)*(REAL(jIter) + 1.)*BoltzmannConst* &
                    PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))**FakXi
    IF (CurrentValue .GT. MaxValue) MaxValue = CurrentValue
  END DO
  CALL RANDOM_NUMBER(iRan)
  iQuant = INT((1+J2)*iRan)
  DO WHILE (ARM)
    fNorm = (2.*REAL(iQuant) + 1.)**2 * &
            (Ec - REAL(iQuant)*(REAL(iQuant) + 1.)*BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))**FakXi / MaxValue
    CALL RANDOM_NUMBER(iRan)
    IF(fNorm .LT. iRan) THEN
      CALL RANDOM_NUMBER(iRan)
      iQuant = INT((1+J2)*iRan)
    ELSE
      ARM = .FALSE.
    END IF
  END DO
  PartStateIntEn( 2,iPart) = REAL(iQuant) * (REAL(iQuant) + 1.) * BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)

ELSE IF(PolyatomMolDSMC(iPolyatMole)%RotationalGroup.EQ.10)THEN
  ! molecule is non-linear -> symmetric top molecule where two are equal and third is different -> oblate top
  ! calculate maximum allowed energy (all of collision energy in rotatinal energy) with k = J
    J2 = INT((-PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3)/PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) + &
    SQRT((PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3)/PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1))**2. + &
    (32.*Ec*PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3)*PI**2.)/(PlanckConst**2.)))/2.)
  ! reduce J2 if too big which would correspond to unphysical quantum numbers, necessary for high velocities since otherwise
  ! almost no samples will be accepted
  IF(J2.GT.500) J2 = 500
  ! Find max value of distribution for ARM numerically
  MaxValue = 0.
  DO jIter=0, J2
    ! for k=0 not double degenerate so first term is added outside of loop
    CurrentValue = (2.*REAL(jIter)+1.) * (Ec - (PlanckConst**2.)/(PI**2.*8.) * &
        (REAL(jIter)*(REAL(jIter)+1.))/(PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)))**FakXi
    IF (CurrentValue .GT. MaxValue) MaxValue = CurrentValue
    DO kIter=1, jIter
      CurrentValue = (2.) * (2.*REAL(jIter)+1.) * (Ec - (PlanckConst**2.)/(PI**2.*8.) * &
                    ((REAL(jIter)*(REAL(jIter)+1.))/(PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) + &
                    (1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) - 1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * &
                    REAL(kIter)**2.))**FakXi
      IF (CurrentValue .GT. MaxValue) MaxValue = CurrentValue
    END DO
  END DO
  ! roll quantum numbers
  CALL RANDOM_NUMBER(iRan)
  iQuant = INT((1+J2)*iRan)
  CALL RANDOM_NUMBER(iRan)
  kQuant = INT((2. * J2 + 1.) * iRan) - J2
  ! check if condition |k| <= j is true and reroll unitll it is
  DO WHILE(kQuant**2.GT.iQuant**2)
  CALL RANDOM_NUMBER(iRan)
  iQuant = INT((1+J2)*iRan)
  CALL RANDOM_NUMBER(iRan)
  kQuant = INT((2. * J2 + 1.) * iRan) - J2
  END DO
  ! acceptance rejection sampling
  DO WHILE (ARM)
  ! set delta for degeneracy
  delta=0
  IF(kQuant.EQ.0) delta=1
  ! check if rotational energy is possible
  IF((Ec - (PlanckConst**2.)/(PI**2.*8.) * &
    ((REAL(iQuant)*(REAL(iQuant)+1))/(PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) + &
    (1/PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) - 1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * &
    REAL(kQuant)**2.)).LT.0)THEN
    fNorm = -1! set fNorm to less than than 1 to redo loop since the quantum numbers were not possible due to the collision energy
  ELSE
    fNorm = ((2-REAL(delta)) * (2.*REAL(iQuant)+1.) * (Ec - (PlanckConst**2.)/(PI**2.*8.) * &
            ((REAL(iQuant)*(REAL(iQuant)+1.))/(PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) + &
            (1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) - 1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * &
            REAL(kQuant)**2.))**FakXi) &
            / MaxValue
  END IF
  CALL RANDOM_NUMBER(iRan)
  IF(fNorm .LT. iRan) THEN
    ! roll new quantum numbers
    CALL RANDOM_NUMBER(iRan)
    iQuant = INT((1+J2)*iRan)
    CALL RANDOM_NUMBER(iRan)
    kQuant = INT((2. * J2 + 1.) * iRan) - J2
    ! check if condition |k| <= j is true and reroll unitll it is
    DO WHILE(kQuant**2.GT.iQuant**2)
      CALL RANDOM_NUMBER(iRan)
      iQuant = INT((1+J2)*iRan)
      CALL RANDOM_NUMBER(iRan)
      kQuant = INT((2. * J2 + 1.) * iRan) - J2
    END DO
  ELSE
    ARM = .FALSE.
  END IF
  END DO

  PartStateIntEn( 2,iPart) = PlanckConst**2. / (8.*PI**2.) * (REAL(iQuant)*(REAL(iQuant)+1.) / &
      PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) + (1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) &
      - 1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * REAL(kQuant)**2)

ELSE IF(PolyatomMolDSMC(iPolyatMole)%RotationalGroup.EQ.11)THEN
  ! molecule is non-linear -> symmetric top molecule where two are equal and third is different -> prolate top
  ! function identical to oblate tops other than J2 calculation but for clarity and less if statements separate
  ! calculate maximum allowed energy (all of collision energy in rotatinal energy) with k = 0
    J2 = INT((-1. + SQRT(1. + (32.*Ec*PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)*PI**2.)/(PlanckConst**2.)))/2.)
  ! reduce J2 if too big which would correspond to unphysical quantum numbers, necessary for high velocities since otherwise
  ! almost no samples will be accepted
  IF(J2.GT.500) J2 = 500
  ! Find max value of distribution for ARM numerically
  MaxValue = 0.
  DO jIter=0, J2
    ! for k=0 not double degenerate so first term is added outside of loop
    CurrentValue = (2.*REAL(jIter)+1.) * (Ec - (PlanckConst**2.)/(PI**2.*8.) * &
        (REAL(jIter)*(REAL(jIter)+1.))/(PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)))**FakXi
    IF (CurrentValue .GT. MaxValue) MaxValue = CurrentValue
    DO kIter=1, jIter
      CurrentValue = (2.) * (2.*REAL(jIter)+1.) * (Ec - (PlanckConst**2.)/(PI**2.*8.) * &
                    ((REAL(jIter)*(REAL(jIter)+1.))/(PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) + &
                    (1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) - 1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * &
                    REAL(kIter)**2.))**FakXi
      IF (CurrentValue .GT. MaxValue) MaxValue = CurrentValue
    END DO
  END DO
  ! roll quantum numbers
  CALL RANDOM_NUMBER(iRan)
  iQuant = INT((1+J2)*iRan)
  CALL RANDOM_NUMBER(iRan)
  kQuant = INT((2. * J2 + 1.) * iRan) - J2
  ! check if condition |k| <= j is true and reroll unitll it is
  DO WHILE(kQuant**2.GT.iQuant**2)
  CALL RANDOM_NUMBER(iRan)
  iQuant = INT((1+J2)*iRan)
  CALL RANDOM_NUMBER(iRan)
  kQuant = INT((2. * J2 + 1.) * iRan) - J2
  END DO
  ! acceptance rejection sampling
  DO WHILE (ARM)
  ! set delta for degeneracy
  delta=0
  IF(kQuant.EQ.0) delta=1
  ! check if rotational energy is possible
  IF((Ec - (PlanckConst**2.)/(PI**2.*8.) * &
    ((REAL(iQuant)*(REAL(iQuant)+1))/(PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) + &
    (1/PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) - 1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * &
    REAL(kQuant)**2.)).LT.0)THEN
    fNorm = -1! set fNorm to less than than 1 to redo loop since the quantum numbers were not possible due to the collision energy
  ELSE
    fNorm = ((2-REAL(delta)) * (2.*REAL(iQuant)+1.) * (Ec - (PlanckConst**2.)/(PI**2.*8.) * &
            ((REAL(iQuant)*(REAL(iQuant)+1.))/(PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) + &
            (1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) - 1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * &
            REAL(kQuant)**2.))**FakXi) &
            / MaxValue
  END IF
  CALL RANDOM_NUMBER(iRan)
  IF(fNorm .LT. iRan) THEN
    ! roll new quantum numbers
    CALL RANDOM_NUMBER(iRan)
    iQuant = INT((1+J2)*iRan)
    CALL RANDOM_NUMBER(iRan)
    kQuant = INT((2. * J2 + 1.) * iRan) - J2
    ! check if condition |k| <= j is true and reroll unitll it is
    DO WHILE(kQuant**2.GT.iQuant**2)
      CALL RANDOM_NUMBER(iRan)
      iQuant = INT((1+J2)*iRan)
      CALL RANDOM_NUMBER(iRan)
      kQuant = INT((2. * J2 + 1.) * iRan) - J2
    END DO
  ELSE
    ARM = .FALSE.
  END IF
  END DO

  PartStateIntEn( 2,iPart) = PlanckConst**2. / (8.*PI**2.) * (REAL(iQuant)*(REAL(iQuant)+1.) / &
      PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) + (1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) &
      - 1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * REAL(kQuant)**2)
ELSE IF(PolyatomMolDSMC(iPolyatMole)%RotationalGroup.EQ.3)THEN
  ! -> asymmetric top molecule, no analytic formula for energy levels -> always use rotational levels from database
  CALL DSMC_RotRelaxDatabasePoly(iPair,iPart,FakXi)
ELSE
  CALL abort(__STAMP__,'Unexpected relations of moments of inertia of species!', iSpec)
END IF

END SUBROUTINE DSMC_RotRelaxQuantPoly


SUBROUTINE DSMC_RotRelaxQuantPolyMH(iPair, iPart, FakXi)
!===================================================================================================================================
! rotational relaxation routine with the Metropolis-Hastings method (no burn-in phase) but uses last quantum numbers from inital
! particle insertion - ONLY for comparison with acceptance rejection sapmling
!===================================================================================================================================
! MODULES
  USE MOD_Globals               ,ONLY: Abort
  USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, SpecDSMC, PolyatomMolDSMC, Coll_pData, RadialWeighting, RotQuantsPar
  USE MOD_Globals_Vars          ,ONLY: BoltzmannConst, PlanckConst, PI
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
  REAL                          :: iRan, tempEng, tempProb, NormProb, Ec
  INTEGER                       :: iPolyatMole, iSpec, iWalk
  INTEGER                       :: jMax, iQuant, kQuant, delta, delta_old
!===================================================================================================================================
iSpec = PartSpecies(iPart)
IF (usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep) THEN
  Ec = Coll_pData(iPair)%Ec / GetParticleWeight(iPart)
ELSE
  Ec = Coll_pData(iPair)%Ec
END IF
iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray

IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec)THEN        ! check if molecule is linear, same as diatomic
  jMax = 440
  DO iWalk=1,1000
    NormProb = Ec - PartStateIntEn(2,iPart)
    ! Proper modelling of energy transfer between old and new state in chemistry
    NormProb = NormProb**FakXi
    CALL RANDOM_NUMBER(iRan)
    iQuant = RotQuantsPar(1,iPart)+NINT((2.0 * iRan - 1.0) * jMax)
    IF(iQuant.LT.0) iQuant = -1*iQuant -1
    tempEng= REAL(iQuant) * (REAL(iQuant) + 1.) * BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)
    tempProb = Ec - tempEng
    IF(tempProb.GT.0) THEN
      NormProb = MIN(1.0,(2.*REAL(iQuant) + 1.)*tempProb**FakXi/((2.*REAL(RotQuantsPar(1,iPart)) + 1.)*NormProb))
      CALL RANDOM_NUMBER(iRan)
      IF(NormProb.GE.iRan) THEN
        PartStateIntEn(2,iPart) = tempEng
        RotQuantsPar(1,iPart)=iQuant
      END IF
    END IF
  END DO
ELSE IF(PolyatomMolDSMC(iPolyatMole)%RotationalGroup.EQ.1)THEN
  jMax = 330
  DO iWalk=1,750
    NormProb = Ec - PartStateIntEn(2,iPart)
    ! Proper modelling of energy transfer between old and new state in chemistry
    NormProb = NormProb**FakXi
    CALL RANDOM_NUMBER(iRan)
    iQuant = RotQuantsPar(1,iPart)+NINT((2.0 * iRan - 1.0) * jMax)
    IF(iQuant.LT.0) iQuant = -1*iQuant -1
    tempEng= REAL(iQuant) * (REAL(iQuant) + 1.) * BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)
    tempProb = Ec - tempEng
    IF(tempProb.GT.0) THEN
      NormProb = MIN(1.0,(2.*REAL(iQuant) + 1.)**2.*tempProb**FakXi/((2.*REAL(RotQuantsPar(1,iPart)) + 1.)**2.*NormProb))
      CALL RANDOM_NUMBER(iRan)
      IF(NormProb.GE.iRan) THEN
        PartStateIntEn(2,iPart) = tempEng
        RotQuantsPar(1,iPart)=iQuant
      END IF
    END IF
  END DO
ELSE IF((PolyatomMolDSMC(iPolyatMole)%RotationalGroup.EQ.10).OR.(PolyatomMolDSMC(iPolyatMole)%RotationalGroup.EQ.11))THEN
  ! molecule is non-linear -> symmetric top molecule where two are equal and third is different
  jMax = 75
  DO iWalk=1,1000
    !//TODO what if NormProb is negative? -> leads to NaN if not value is found in iWalk steps
    NormProb = Ec - PartStateIntEn(2,iPart)
    ! Proper modelling of energy transfer between old and new state in chemistry
    NormProb = NormProb**FakXi
    CALL RANDOM_NUMBER(iRan)
    iQuant = RotQuantsPar(1,iPart)+NINT((2.0 * iRan - 1.0) * 100)
    CALL RANDOM_NUMBER(iRan)
    kQuant = RotQuantsPar(2,iPart)+NINT((2.0 * iRan - 1.0) * 100)
    ! set delta for degeneracy
    delta=0
    IF(kQuant.EQ.0) delta=1
    delta_old=0
    IF(RotQuantsPar(2,iPart).EQ.0) delta_old=1
    IF(iQuant.LT.0) iQuant = -1*iQuant -1
    tempEng = PlanckConst**2. / (8.*PI**2.) * (REAL(iQuant)*(REAL(iQuant)+1.) / &
      PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1) + (1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(3) - &
      1./PolyatomMolDSMC(iPolyatMole)%MomentOfInertia(1)) * REAL(kQuant)**2.)
    tempProb = Ec - tempEng
    IF(tempProb.GT.0) THEN
      NormProb = MIN(1.0,(2.-REAL(delta))*(2.*REAL(iQuant) + 1.)*tempProb**FakXi / &
                ((2.-REAL(delta_old))*(2.*REAL(RotQuantsPar(1,iPart)) + 1.)*NormProb))
      CALL RANDOM_NUMBER(iRan)
      IF((NormProb.GE.iRan).AND.(kQuant**2.LE.iQuant**2)) THEN
        PartStateIntEn(2,iPart) = tempEng
        RotQuantsPar(1,iPart)=iQuant
        RotQuantsPar(2,iPart)=kQuant
      END IF
    END IF
  END DO

ELSE IF(PolyatomMolDSMC(iPolyatMole)%RotationalGroup.EQ.3)THEN
  ! -> asymmetric top molecule, no analytic formula for energy levels -> always use rotational levels
  CALL DSMC_RotRelaxDatabasePoly(iPair,iPart,FakXi)
ELSE
  CALL abort(__STAMP__,'Unexpected relations of moments of inertia of species!')
END IF
END SUBROUTINE DSMC_RotRelaxQuantPolyMH


SUBROUTINE DSMC_RotRelaxDatabasePoly(iPair,iPart,FakXi)
!===================================================================================================================================
!> Rotational relaxation with database energy levels
!===================================================================================================================================
  USE MOD_Globals
  USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
  USE MOD_DSMC_Vars              ,ONLY: SpecDSMC, PartStateIntEn, RadialWeighting, Coll_pData
  USE MOD_Particle_Vars          ,ONLY: PartSpecies, UseVarTimeStep, usevMPF
  USE MOD_part_tools             ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iPair, iPart
  REAL, INTENT(IN)              :: FakXi
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: iQuaMax, MaxRotQuant, iQua, iSpec
  REAL                          :: iRan, iRan2, gmax, gtemp, PartStateTemp, CollisionEnergy
!===================================================================================================================================
iSpec = PartSpecies(iPart)
IF (usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep) THEN
  CollisionEnergy = Coll_pData(iPair)%Ec / GetParticleWeight(iPart)
ELSE
  CollisionEnergy = Coll_pData(iPair)%Ec
END IF

iQuaMax  = 0
! Determine max rotational quant
MaxRotQuant = SpecDSMC(iSpec)%MaxRotQuant - 1
gmax = 0.
PartStateTemp = CollisionEnergy / BoltzmannConst
DO iQua = 0, MaxRotQuant
  IF (PartStateTemp - SpecDSMC(iSpec)%RotationalState(2,iQua).GT.0.) THEN
    gtemp = SpecDSMC(iSpec)%RotationalState(1,iQua) * &
            ( CollisionEnergy - BoltzmannConst * SpecDSMC(iSpec)%RotationalState(2,iQua))**FakXi
    ! maximal possible Quant before term goes negative
    iQuaMax = iQua
    IF ( gtemp .GT. gmax ) THEN
      gmax = gtemp
    END IF
  ELSE
    EXIT
  END IF
END DO
CALL RANDOM_NUMBER(iRan)
iQua = int( ( iQuaMax +1 ) * iRan)
gtemp = SpecDSMC(iSpec)%RotationalState(1,iQua) * &
        ( CollisionEnergy - BoltzmannConst * SpecDSMC(iSpec)%RotationalState(2,iQua))**FakXi
CALL RANDOM_NUMBER(iRan2)
! acceptance-rejection for iQuaRot
DO WHILE ( iRan2 .GE. gtemp / gmax )
  CALL RANDOM_NUMBER(iRan)
  iQua = int( ( iQuaMax +1 ) * iRan)
  gtemp = SpecDSMC(iSpec)%RotationalState(1,iQua) * &
          ( CollisionEnergy - BoltzmannConst * SpecDSMC(iSpec)%RotationalState(2,iQua))**FakXi
  CALL RANDOM_NUMBER(iRan2)
END DO
PartStateIntEn(2,iPart) = BoltzmannConst * SpecDSMC(iSpec)%RotationalState(2,iQua)

END SUBROUTINE DSMC_RotRelaxDatabasePoly


END MODULE MOD_DSMC_PolyAtomicModel