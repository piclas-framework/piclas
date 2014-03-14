#include "boltzplatz.h"

MODULE MOD_DSMC_PolyAtomicModel
!===================================================================================================================================
! module including qk procedures
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
  PUBLIC :: InitPolyAtomicMolecs, DSMC_SetInternalEnr_Poly, DSMC_VibRelaxPoly, DSMC_RotRelaxPoly, DSMC_SetInternalEnr_PolyFast
  PUBLIC :: DSMC_VibRelaxPolyFast, DSMC_SetInternalEnr_PolyFastPart2
!-----------------------------------------------------------------------------------------------------------------------------------
  CONTAINS

!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE InitPolyAtomicMolecs(iSpec)
!--------------------------------------------------------------------------------------------------!
! init electronic shell
!--------------------------------------------------------------------------------------------------!
  USE MOD_DSMC_Vars,              ONLY : DSMC, SpecDSMC, PartStateIntEn,PolyatomMolDSMC
  USE MOD_Particle_Vars,          ONLY : BoltzmannConst 
  USE MOD_ReadInTools
!--------------------------------------------------------------------------------------------------!
  IMPLICIT NONE                                                                                    !
!  argument list declaration                                                                       !
!--------------------------------------------------------------------------------------------------!
!  input variable declaration                                                                      !
  INTEGER, INTENT(IN)           ::  iSpec                                                     !
!--------------------------------------------------------------------------------------------------!
!  Local variable declaration                                                                      !
  CHARACTER(32)         :: hilf, hilf2
  INTEGER               :: iPolyatMole, iVibDOF                     
  REAL                  :: JToEv
!--------------------------------------------------------------------------------------------------!
JToEv = 1.602176565E-19  
WRITE(UNIT=hilf,FMT='(I2)') iSpec
iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
SpecDSMC(iSpec)%TVib       = GETREAL('Part-Species'//TRIM(hilf)//'-TempVib','0.')
SpecDSMC(iSpec)%TRot       = GETREAL('Part-Species'//TRIM(hilf)//'-TempRot','0.')  
PolyatomMolDSMC(iPolyatMole)%LinearMolec  =GETLOGICAL('Part-Species'//TRIM(hilf)//'-LinearMolec','.TRUE.')
IF (PolyatomMolDSMC(iPolyatMole)%LinearMolec) THEN
  SpecDSMC(iSpec)%Xi_Rot = 2
ELSE
  SpecDSMC(iSpec)%Xi_Rot = 3
END IF
PolyatomMolDSMC(iPolyatMole)%NumOfAtoms = GETINT('Part-Species'//TRIM(hilf)//'-NumOfAtoms','0')
IF(PolyatomMolDSMC(iPolyatMole)%NumOfAtoms.EQ.0)THEN
  WRITE(*,*) 'ERROR in Number of Atoms in polyatomic Molec Spec: ', iSpec
  STOP
END IF
SpecDSMC(iSpec)%Ediss_eV   = GETREAL('Part-Species'//TRIM(hilf)//'-Ediss_eV','0.')

PolyatomMolDSMC(iPolyatMole)%VibDOF = 3*PolyatomMolDSMC(iPolyatMole)%NumOfAtoms - 3 - SpecDSMC(iSpec)%Xi_Rot
ALLOCATE(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(PolyatomMolDSMC(iPolyatMole)%VibDOF))
ALLOCATE(PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(PolyatomMolDSMC(iPolyatMole)%VibDOF))
DO iVibDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF 
  WRITE(UNIT=hilf2,FMT='(I2)') iVibDOF
  PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iVibDOF) =  &
    GETREAL('Part-Species'//TRIM(hilf)//'-CharaTempVib'//TRIM(hilf2),'0.')
END DO
ALLOCATE(PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(PolyatomMolDSMC(iPolyatMole)%VibDOF))
PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(1:PolyatomMolDSMC(iPolyatMole)%VibDOF) = 80
!PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(1:PolyatomMolDSMC(iPolyatMole)%VibDOF) = &
!        INT(SpecDSMC(iSpec)%Ediss_eV*JToEv &
!        /(BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(1:PolyatomMolDSMC(iPolyatMole)%VibDOF))) + 1

SpecDSMC(iSpec)%VFD_Phi3_Factor = GETREAL('Part-Species'//TRIM(hilf)//'-VFDPhi3','0.')
! Setting the values of Rot-/Vib-RelaxProb to a fix value!
! This should be changed to a calculated value for every coll pair/situation!!!1
SpecDSMC(iSpec)%RotRelaxProb  = 0.2!0.2
SpecDSMC(iSpec)%VibRelaxProb  = 0.02!0.02
SpecDSMC(iSpec)%ElecRelaxProb = 0.01!or 0.02 | Bird: somewhere in range 0.01 .. 0.02
CALL DSMC_FindFirstVibPick(iPolyatMole, iSpec)
END SUBROUTINE InitPolyAtomicMolecs

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE DSMC_SetInternalEnr_Poly(iSpec, iPart)

  USE MOD_DSMC_Vars,            ONLY : PartStateIntEn, SpecDSMC, DSMC,PolyatomMolDSMC,VibQuantsPar
  USE MOD_Particle_Vars,        ONLY : BoltzmannConst
  USE MOD_DSMC_ElectronicModel, ONLY : InitElectronShell
!--------------------------------------------------------------------------------------------------!
! perform chemical init
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE 
! LOCAL VARIABLES
  INTEGER, INTENT(IN)           :: iSpec, iPart
  REAL, ALLOCATABLE             :: iRan(:), tempEng(:)
  REAL                          :: iRan2, NormProb
  INTEGER,ALLOCATABLE           :: iQuant(:)
  INTEGER                       :: iDOF,iPolyatMole
!#ifdef MPI
!#endif
!===================================================================================================================================

!set vibrational energy
IF (SpecDSMC(iSpec)%PolyatomicMol) THEN
  iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
  ALLOCATE(VibQuantsPar(iPart)%Quants(PolyatomMolDSMC(iPolyatMole)%VibDOF))  
  ALLOCATE(iRan(PolyatomMolDSMC(iPolyatMole)%VibDOF) &
          ,tempEng(PolyatomMolDSMC(iPolyatMole)%VibDOF) &
          ,iQuant(PolyatomMolDSMC(iPolyatMole)%VibDOF))
  CALL RANDOM_NUMBER(iRan)
  iQuant(:) = INT(iRan(:)*PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(:))
  tempEng(:)=iQuant(:)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:)/SpecDSMC(iSpec)%TVib
  NormProb = 1.0
  DO iDOF = 1 , PolyatomMolDSMC(iPolyatMole)%VibDOF
    NormProb = NormProb*EXP(-tempEng(iDOF))
  END DO
  CALL RANDOM_NUMBER(iRan2)
  DO WHILE (iRan2.GE.NormProb)
    CALL RANDOM_NUMBER(iRan)
    iQuant(:) = INT(iRan(:)*PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(:))
    tempEng(:)=iQuant(:)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:)/SpecDSMC(iSpec)%TVib
    NormProb = 1.0
    DO iDOF = 1 , PolyatomMolDSMC(iPolyatMole)%VibDOF
      NormProb = NormProb*EXP(-tempEng(iDOF))
    END DO
    CALL RANDOM_NUMBER(iRan2)
  END DO
  !evtl muß partstateinten nochmal geändert werden, mpi, resize etc..
  PartStateIntEn(iPart, 1) = 0.0
  VibQuantsPar(iPart)%Quants(:)=iQuant(:)
  DO iDOF = 1 , PolyatomMolDSMC(iPolyatMole)%VibDOF
    PartStateIntEn(iPart, 1)= PartStateIntEn(iPart, 1) &
      +(iQuant(iDOF) + DSMC%GammaQuant)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
  END DO

!set rotational energy
  IF (SpecDSMC(iSpec)%Xi_Rot.EQ.2) THEN
    CALL RANDOM_NUMBER(iRan2)
    PartStateIntEn(iPart, 2) = -BoltzmannConst*SpecDSMC(iSpec)%TRot*LOG(iRan2)
  ELSE IF (SpecDSMC(iSpec)%Xi_Rot.EQ.3) THEN
    CALL RANDOM_NUMBER(iRan2)
    PartStateIntEn(iPart, 2) = iRan2*10 !the distribution function has only non-negligible  values betwenn 0 and 10
    NormProb = SQRT(PartStateIntEn(iPart, 2))*EXP(-PartStateIntEn(iPart, 2))/(SQRT(0.5)*EXP(-0.5))
    CALL RANDOM_NUMBER(iRan2)
    DO WHILE (iRan2.GE.NormProb)
      CALL RANDOM_NUMBER(iRan2)
      PartStateIntEn(iPart, 2) = iRan2*10 !the distribution function has only non-negligible  values betwenn 0 and 10
      NormProb = SQRT(PartStateIntEn(iPart, 2))*EXP(-PartStateIntEn(iPart, 2))/(SQRT(0.5)*EXP(-0.5))
      CALL RANDOM_NUMBER(iRan2)
    END DO
    PartStateIntEn(iPart, 2) = PartStateIntEn(iPart, 2)*BoltzmannConst*SpecDSMC(iSpec)%TRot
  END IF
  DEALLOCATE(iRan, tempEng, iQuant)
ELSE
  PartStateIntEn(iPart, 1) = 0
  PartStateIntEn(iPart, 2) = 0
END IF


END SUBROUTINE DSMC_SetInternalEnr_Poly

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE DSMC_SetInternalEnr_PolyFast(iSpec, NumParts)

  USE MOD_DSMC_Vars,            ONLY : PartStateIntEn, SpecDSMC, DSMC,PolyatomMolDSMC,VibQuantsPar
  USE MOD_Particle_Vars,        ONLY : BoltzmannConst
  USE MOD_DSMC_ElectronicModel, ONLY : InitElectronShell
!--------------------------------------------------------------------------------------------------!
! perform chemical init
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE 
! LOCAL VARIABLES
  INTEGER, INTENT(IN)           :: iSpec, NumParts
  REAL, ALLOCATABLE             :: iRan(:), tempEng(:)
  REAL                          :: iRan2, NormProb
  INTEGER,ALLOCATABLE           :: iQuant(:), iQuant_old(:)
  INTEGER                       :: iDOF,iPolyatMole,iPart, maxdelta
!#ifdef MPI
!#endif
!===================================================================================================================================
maxdelta=4
!set vibrational energy
!IF (SpecDSMC(iSpec)%PolyatomicMol) THEN
iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray

ALLOCATE(iRan(PolyatomMolDSMC(iPolyatMole)%VibDOF) &
        ,tempEng(PolyatomMolDSMC(iPolyatMole)%VibDOF) &
        ,iQuant(PolyatomMolDSMC(iPolyatMole)%VibDOF) &
        ,iQuant_old(PolyatomMolDSMC(iPolyatMole)%VibDOF))

DO iDOF = 1, NumParts
  ALLOCATE(VibQuantsPar(iDOF)%Quants(PolyatomMolDSMC(iPolyatMole)%VibDOF))
END DO
CALL RANDOM_NUMBER(iRan)
iQuant(:) = INT(iRan(:)*PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(:))
PartStateIntEn(1, 1)=0.0
DO iDOF = 1 , PolyatomMolDSMC(iPolyatMole)%VibDOF
  PartStateIntEn(1, 1)= PartStateIntEn(1, 1) &
    +(iQuant(iDOF) + DSMC%GammaQuant)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
END DO
VibQuantsPar(1)%Quants(:)=iQuant(:)

IF (SpecDSMC(iSpec)%Xi_Rot.EQ.2) THEN
  CALL RANDOM_NUMBER(iRan2)
  PartStateIntEn(1, 2) = -BoltzmannConst*SpecDSMC(iSpec)%TRot*LOG(iRan2)
ELSE IF (SpecDSMC(iSpec)%Xi_Rot.EQ.3) THEN
  CALL RANDOM_NUMBER(iRan2)
  PartStateIntEn(1, 2) = iRan2*10 !the distribution function has only non-negligible  values betwenn 0 and 10
  NormProb = SQRT(PartStateIntEn(1, 2))*EXP(-PartStateIntEn(1, 2))/(SQRT(0.5)*EXP(-0.5))
  CALL RANDOM_NUMBER(iRan2)
  DO WHILE (iRan2.GE.NormProb)
    CALL RANDOM_NUMBER(iRan2)
    PartStateIntEn(1, 2) = iRan2*10 !the distribution function has only non-negligible  values betwenn 0 and 10
    NormProb = SQRT(PartStateIntEn(1, 2))*EXP(-PartStateIntEn(1, 2))/(SQRT(0.5)*EXP(-0.5))
    CALL RANDOM_NUMBER(iRan2)
  END DO
  PartStateIntEn(1, 2) = PartStateIntEn(1, 2)*BoltzmannConst*SpecDSMC(iSpec)%TRot
END IF

DO iPart=2, NumParts
  iQuant_old(:)=iQuant(:)
  CALL RANDOM_NUMBER(iRan)
!  iQuant(:) = iQuant_old(:)+INT(maxdelta*(iRan(:)*2-1))
  iQuant(:) = iQuant_old(:)+FLOOR(3*iRan(:)-1)
  NormProb = 0.0
  DO iDOF = 1 , PolyatomMolDSMC(iPolyatMole)%VibDOF
    IF(iQuant(iDOF).LT.0) iQuant(iDOF) = -1*iQuant(iDOF) - 1
    NormProb = NormProb + (iQuant_old(iDOF)-iQuant(iDOF))*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/SpecDSMC(iSpec)%TVib
  END DO
  NormProb = MIN(1.0,EXP(NormProb))
  CALL RANDOM_NUMBER(iRan2)
  IF (NormProb.GE.iRan2)  THEN
    PartStateIntEn(iPart, 1) = 0.0
    DO iDOF = 1 , PolyatomMolDSMC(iPolyatMole)%VibDOF
      PartStateIntEn(iPart, 1)= PartStateIntEn(iPart, 1) &
        +(iQuant(iDOF) + DSMC%GammaQuant)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
    END DO
  ELSE
    iQuant(:)=iQuant_old(:)
    PartStateIntEn(iPart, 1) = 0.0
    DO iDOF = 1 , PolyatomMolDSMC(iPolyatMole)%VibDOF
      PartStateIntEn(iPart, 1)= PartStateIntEn(iPart, 1) &
        +(iQuant(iDOF) + DSMC%GammaQuant)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
    END DO
  END IF
  VibQuantsPar(iPart)%Quants(:)=iQuant(:)

!set rotational energy
  IF (SpecDSMC(iSpec)%Xi_Rot.EQ.2) THEN
    CALL RANDOM_NUMBER(iRan2)
    PartStateIntEn(iPart, 2) = -BoltzmannConst*SpecDSMC(iSpec)%TRot*LOG(iRan2)
  ELSE IF (SpecDSMC(iSpec)%Xi_Rot.EQ.3) THEN
    CALL RANDOM_NUMBER(iRan2)
    PartStateIntEn(iPart, 2) = iRan2*10 !the distribution function has only non-negligible  values betwenn 0 and 10
    NormProb = SQRT(PartStateIntEn(iPart, 2))*EXP(-PartStateIntEn(iPart, 2))/(SQRT(0.5)*EXP(-0.5))
    CALL RANDOM_NUMBER(iRan2)
    DO WHILE (iRan2.GE.NormProb)
      CALL RANDOM_NUMBER(iRan2)
      PartStateIntEn(iPart, 2) = iRan2*10 !the distribution function has only non-negligible  values betwenn 0 and 10
      NormProb = SQRT(PartStateIntEn(iPart, 2))*EXP(-PartStateIntEn(iPart, 2))/(SQRT(0.5)*EXP(-0.5))
      CALL RANDOM_NUMBER(iRan2)
    END DO
    PartStateIntEn(iPart, 2) = PartStateIntEn(iPart, 2)*BoltzmannConst*SpecDSMC(iSpec)%TRot
  END IF

END DO

DEALLOCATE(iRan, tempEng, iQuant,iQuant_old)
!ELSE
!  PartStateIntEn(iPart, 1) = 0
!  PartStateIntEn(iPart, 2) = 0
!END IF


END SUBROUTINE DSMC_SetInternalEnr_PolyFast

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------


SUBROUTINE DSMC_SetInternalEnr_PolyFastPart(iSpec, iPart)

  USE MOD_DSMC_Vars,            ONLY : PartStateIntEn, SpecDSMC, DSMC,PolyatomMolDSMC,VibQuantsPar
  USE MOD_Particle_Vars,        ONLY : BoltzmannConst
  USE MOD_DSMC_ElectronicModel, ONLY : InitElectronShell
!--------------------------------------------------------------------------------------------------!
! perform chemical init
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE 
! LOCAL VARIABLES
  INTEGER, INTENT(IN)           :: iSpec, iPart
  REAL, ALLOCATABLE             :: iRan(:)
  REAL                          :: iRan2, NormProb
  INTEGER,ALLOCATABLE           :: iQuant(:), iQuant_old(:)
  INTEGER                       :: iDOF,iPolyatMole, iWalk
!#ifdef MPI
!#endif
!===================================================================================================================================

IF (SpecDSMC(iSpec)%PolyatomicMol) THEN
  iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
  ALLOCATE(iRan(PolyatomMolDSMC(iPolyatMole)%VibDOF) &          
          ,iQuant(PolyatomMolDSMC(iPolyatMole)%VibDOF) &
          ,iQuant_old(PolyatomMolDSMC(iPolyatMole)%VibDOF))
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
      NormProb = NormProb + (iQuant_old(iDOF)-iQuant(iDOF))*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/SpecDSMC(iSpec)%TVib
    END DO    
    NormProb = MIN(1.0,EXP(NormProb))
    CALL RANDOM_NUMBER(iRan2)
    IF (NormProb.LT.iRan2) iQuant(:)=iQuant_old(:)   
  END DO

  PartStateIntEn(iPart, 1) = 0.0
  DO iDOF = 1 , PolyatomMolDSMC(iPolyatMole)%VibDOF
    PartStateIntEn(iPart, 1)= PartStateIntEn(iPart, 1) &
      +(iQuant(iDOF) + DSMC%GammaQuant)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
  END DO
  VibQuantsPar(iPart)%Quants(:)=iQuant(:)
  DEALLOCATE(iRan, iQuant,iQuant_old)

 !set rotational energy
  IF (SpecDSMC(iSpec)%Xi_Rot.EQ.2) THEN
    CALL RANDOM_NUMBER(iRan2)
    PartStateIntEn(iPart, 2) = -BoltzmannConst*SpecDSMC(iSpec)%TRot*LOG(iRan2)
  ELSE IF (SpecDSMC(iSpec)%Xi_Rot.EQ.3) THEN
    CALL RANDOM_NUMBER(iRan2)
    PartStateIntEn(iPart, 2) = iRan2*10 !the distribution function has only non-negligible  values betwenn 0 and 10
    NormProb = SQRT(PartStateIntEn(iPart, 2))*EXP(-PartStateIntEn(iPart, 2))/(SQRT(0.5)*EXP(-0.5))
    CALL RANDOM_NUMBER(iRan2)
    DO WHILE (iRan2.GE.NormProb)
      CALL RANDOM_NUMBER(iRan2)
      PartStateIntEn(iPart, 2) = iRan2*10 !the distribution function has only non-negligible  values betwenn 0 and 10
      NormProb = SQRT(PartStateIntEn(iPart, 2))*EXP(-PartStateIntEn(iPart, 2))/(SQRT(0.5)*EXP(-0.5))
      CALL RANDOM_NUMBER(iRan2)
    END DO
    PartStateIntEn(iPart, 2) = PartStateIntEn(iPart, 2)*BoltzmannConst*SpecDSMC(iSpec)%TRot
  END IF
ELSE
  PartStateIntEn(iPart, 1) = 0
  PartStateIntEn(iPart, 2) = 0
END IF

END SUBROUTINE DSMC_SetInternalEnr_PolyFastPart

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE DSMC_SetInternalEnr_PolyFastPart2(iSpec, iPart)

  USE MOD_DSMC_Vars,            ONLY : PartStateIntEn, SpecDSMC, DSMC,PolyatomMolDSMC,VibQuantsPar
  USE MOD_Particle_Vars,        ONLY : BoltzmannConst
  USE MOD_DSMC_ElectronicModel, ONLY : InitElectronShell
!--------------------------------------------------------------------------------------------------!
! perform chemical init
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE 
! LOCAL VARIABLES
  INTEGER, INTENT(IN)           :: iSpec, iPart
  REAL, ALLOCATABLE             :: iRan(:)
  REAL                          :: iRan2, NormProb
  INTEGER,ALLOCATABLE           :: iQuant_old(:)
  INTEGER                       :: iDOF,iPolyatMole, iWalk
!#ifdef MPI
!#endif
!===================================================================================================================================

IF (SpecDSMC(iSpec)%PolyatomicMol) THEN
  iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
  ALLOCATE(iRan(PolyatomMolDSMC(iPolyatMole)%VibDOF) &          
          ,iQuant_old(PolyatomMolDSMC(iPolyatMole)%VibDOF))
  ALLOCATE(VibQuantsPar(iPart)%Quants(PolyatomMolDSMC(iPolyatMole)%VibDOF))

  DO iWalk = 1, 150
    iQuant_old(:)=PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(:)
    CALL RANDOM_NUMBER(iRan)
    PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(:) = iQuant_old(:)+FLOOR(3*iRan(:)-1)
    NormProb = 0.0
    DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
      IF(PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(iDOF).LT.0) &
          PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(iDOF) = -1*PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(iDOF) -1
      NormProb = NormProb + (iQuant_old(iDOF) &
        -PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(iDOF))*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/SpecDSMC(iSpec)%TVib
    END DO    
    NormProb = MIN(1.0,EXP(NormProb))
    CALL RANDOM_NUMBER(iRan2)
    IF (NormProb.LT.iRan2) PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(:)=iQuant_old(:)   
  END DO

  PartStateIntEn(iPart, 1) = 0.0
  DO iDOF = 1 , PolyatomMolDSMC(iPolyatMole)%VibDOF
    PartStateIntEn(iPart, 1)= PartStateIntEn(iPart, 1) &
      +(PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(iDOF) &
      + DSMC%GammaQuant)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
  END DO
  VibQuantsPar(iPart)%Quants(:)=PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(:)
  DEALLOCATE(iRan, iQuant_old)

 !set rotational energy
  IF (SpecDSMC(iSpec)%Xi_Rot.EQ.2) THEN
    CALL RANDOM_NUMBER(iRan2)
    PartStateIntEn(iPart, 2) = -BoltzmannConst*SpecDSMC(iSpec)%TRot*LOG(iRan2)
  ELSE IF (SpecDSMC(iSpec)%Xi_Rot.EQ.3) THEN
    CALL RANDOM_NUMBER(iRan2)
    PartStateIntEn(iPart, 2) = iRan2*10 !the distribution function has only non-negligible  values betwenn 0 and 10
    NormProb = SQRT(PartStateIntEn(iPart, 2))*EXP(-PartStateIntEn(iPart, 2))/(SQRT(0.5)*EXP(-0.5))
    CALL RANDOM_NUMBER(iRan2)
    DO WHILE (iRan2.GE.NormProb)
      CALL RANDOM_NUMBER(iRan2)
      PartStateIntEn(iPart, 2) = iRan2*10 !the distribution function has only non-negligible  values betwenn 0 and 10
      NormProb = SQRT(PartStateIntEn(iPart, 2))*EXP(-PartStateIntEn(iPart, 2))/(SQRT(0.5)*EXP(-0.5))
      CALL RANDOM_NUMBER(iRan2)
    END DO
    PartStateIntEn(iPart, 2) = PartStateIntEn(iPart, 2)*BoltzmannConst*SpecDSMC(iSpec)%TRot
  END IF
ELSE
  PartStateIntEn(iPart, 1) = 0
  PartStateIntEn(iPart, 2) = 0
END IF

END SUBROUTINE DSMC_SetInternalEnr_PolyFastPart2

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE DSMC_FindFirstVibPick(iPolyatMole, iSpec)

  USE MOD_DSMC_Vars,            ONLY : PartStateIntEn, SpecDSMC, DSMC,PolyatomMolDSMC,VibQuantsPar
  USE MOD_Particle_Vars,        ONLY : BoltzmannConst
  USE MOD_DSMC_ElectronicModel, ONLY : InitElectronShell
!--------------------------------------------------------------------------------------------------!
! perform chemical init
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE 
! LOCAL VARIABLES
  INTEGER, INTENT(IN)           :: iPolyatMole, iSpec
  REAL, ALLOCATABLE             :: iRan(:)
  REAL                          :: iRan2, NormProb
  INTEGER,ALLOCATABLE           :: iQuant_old(:)
  INTEGER                       :: iDOF, iWalk
!#ifdef MPI
!#endif
!===================================================================================================================================
 
  ALLOCATE(iRan(PolyatomMolDSMC(iPolyatMole)%VibDOF) &
          ,iQuant_old(PolyatomMolDSMC(iPolyatMole)%VibDOF))

  CALL RANDOM_NUMBER(iRan)
  PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(:) = INT(iRan(:)*PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(:))

  DO iWalk=1, 5000
    iQuant_old(:)=PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(:)
    CALL RANDOM_NUMBER(iRan)
    PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(:) = iQuant_old(:)+FLOOR(3*iRan(:)-1)
    NormProb = 0.0
    DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
      IF(PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(iDOF).LT.0) &
          PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(iDOF) = -1*PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(iDOF) -1
      NormProb = NormProb + (iQuant_old(iDOF) &
          -PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(iDOF))*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/SpecDSMC(iSpec)%TVib
    END DO    
    NormProb = MIN(1.0,EXP(NormProb))
    CALL RANDOM_NUMBER(iRan2)
    IF (NormProb.LT.iRan2) PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(:)=iQuant_old(:)   
  END DO

  DEALLOCATE(iRan, iQuant_old)

END SUBROUTINE DSMC_FindFirstVibPick

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE DSMC_VibRelaxPoly(Ec,iSpec, iPart,FakXi)

  USE MOD_DSMC_Vars,            ONLY : PartStateIntEn, SpecDSMC, DSMC,PolyatomMolDSMC
  USE MOD_Particle_Vars,        ONLY : BoltzmannConst
  USE MOD_DSMC_ElectronicModel, ONLY : InitElectronShell
!--------------------------------------------------------------------------------------------------!
! perform chemical init
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE 
! LOCAL VARIABLES
  INTEGER, INTENT(IN)           :: iSpec, iPart
  REAL, INTENT(IN)              :: FakXi, Ec
  REAL, ALLOCATABLE             :: iRan(:), tempEng(:)
  REAL                          :: iRan2, NormProb, tempProb
  INTEGER,ALLOCATABLE           :: iQuant(:), iMaxQuant(:)
  INTEGER                       :: iDOF,iPolyatMole
!#ifdef MPI
!#endif
!===================================================================================================================================
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
CALL RANDOM_NUMBER(iRan)
iMaxQuant(:) = INT(Ec/(BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))) + 1
iMaxQuant(:) = MIN(iMaxQuant(:), PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(:))

iQuant(:)=INT(iRan(:)*iMaxQuant(:))
tempEng(:)=(iQuant(:) + 0.5)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:)*BoltzmannConst
tempProb = 0.0
DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
  tempProb = tempProb + tempEng(iDOF)
END DO
CALL RANDOM_NUMBER(iRan2)
DO WHILE ((iRan2.GE.((Ec-tempProb)**FakXi/NormProb)).OR.(Ec-tempProb.LT.0.0))
  CALL RANDOM_NUMBER(iRan)
  iQuant(:)=INT(iRan(:)*iMaxQuant(:))
  tempEng(:)=(iQuant(:) + 0.5)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:)*BoltzmannConst
  tempProb = 0.0
  DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
    tempProb = tempProb + tempEng(iDOF)
  END DO
  CALL RANDOM_NUMBER(iRan2)
END DO
PartStateIntEn(iPart,1)=tempProb
DEALLOCATE(iRan ,tempEng ,iQuant ,iMaxQuant)

END SUBROUTINE DSMC_VibRelaxPoly

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------


SUBROUTINE DSMC_VibRelaxPolyFast(Ec,iSpec, iPart,FakXi)

  USE MOD_DSMC_Vars,            ONLY : PartStateIntEn, SpecDSMC, DSMC,PolyatomMolDSMC,VibQuantsPar
  USE MOD_Particle_Vars,        ONLY : BoltzmannConst
  USE MOD_DSMC_ElectronicModel, ONLY : InitElectronShell
!--------------------------------------------------------------------------------------------------!
! perform chemical init
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE 
! LOCAL VARIABLES
  INTEGER, INTENT(IN)           :: iSpec, iPart
  REAL, INTENT(IN)              :: FakXi, Ec
  REAL, ALLOCATABLE             :: iRan(:), tempEng(:)
  REAL                          :: iRan2, NormProb, tempProb, finalProb
  INTEGER,ALLOCATABLE           :: iQuant(:), iMaxQuant(:)
  INTEGER                       :: iDOF,iPolyatMole, iWalk
!#ifdef MPI
!#endif
!===================================================================================================================================
iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
ALLOCATE(iRan(PolyatomMolDSMC(iPolyatMole)%VibDOF) &
        ,tempEng(PolyatomMolDSMC(iPolyatMole)%VibDOF) &
        ,iQuant(PolyatomMolDSMC(iPolyatMole)%VibDOF) &
        ,iMaxQuant(PolyatomMolDSMC(iPolyatMole)%VibDOF))

DO iWalk=1,150
  NormProb = Ec - PartStateIntEn(iPart,1)
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

  NormProb = MIN(1.0,tempProb**FakXi/NormProb)
  CALL RANDOM_NUMBER(iRan2)
  IF ((NormProb.GE.iRan2).AND.(tempProb.GT.0.0))  THEN
    PartStateIntEn(iPart,1) = 0.0
    DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
      PartStateIntEn(iPart,1) = PartStateIntEn(iPart,1) + tempEng(iDOF)
    END DO
    VibQuantsPar(iPart)%Quants(:) = iQuant(:)
  END IF

END DO

END SUBROUTINE DSMC_VibRelaxPolyFast

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE DSMC_RotRelaxPoly(Ec, iPart,FakXi)

  USE MOD_DSMC_Vars,            ONLY : PartStateIntEn, SpecDSMC, DSMC,PolyatomMolDSMC
  USE MOD_Particle_Vars,        ONLY : BoltzmannConst
  USE MOD_DSMC_ElectronicModel, ONLY : InitElectronShell
!--------------------------------------------------------------------------------------------------!
! perform chemical init
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE 
! LOCAL VARIABLES
  INTEGER, INTENT(IN)           :: iPart
  REAL, INTENT(IN)              :: FakXi, Ec
  REAL                          :: iRan2, NormProb, tempProb, fak1, fak2
  INTEGER                       :: iDOF,iPolyatMole
!#ifdef MPI
!#endif
!===================================================================================================================================
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
PartStateIntEn(iPart,2)=tempProb


END SUBROUTINE DSMC_RotRelaxPoly

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

END MODULE MOD_DSMC_PolyAtomicModel

