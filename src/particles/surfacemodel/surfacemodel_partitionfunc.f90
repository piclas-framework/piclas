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

MODULE MOD_SurfaceModel_PartFunc
!===================================================================================================================================
! Module for surface model tools
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
INTERFACE PartitionFuncGas
  MODULE PROCEDURE PartitionFuncGas
END INTERFACE

INTERFACE PartitionFuncSurf
  MODULE PROCEDURE PartitionFuncSurf
END INTERFACE

INTERFACE PartitionFuncActAdsorb
  MODULE PROCEDURE PartitionFuncActAdsorb
END INTERFACE

INTERFACE PartitionFuncActDissAdsorb
  MODULE PROCEDURE PartitionFuncActDissAdsorb
END INTERFACE

INTERFACE PartitionFuncActER
  MODULE PROCEDURE PartitionFuncActER
END INTERFACE

INTERFACE PartitionFuncActDesorb
  MODULE PROCEDURE PartitionFuncActDesorb
END INTERFACE

INTERFACE PartitionFuncActDissSurf
  MODULE PROCEDURE PartitionFuncActDissSurf
END INTERFACE

INTERFACE PartitionFuncActLH
  MODULE PROCEDURE PartitionFuncActLH
END INTERFACE

INTERFACE PartitionFuncActExchSurf
  MODULE PROCEDURE PartitionFuncActExchSurf
END INTERFACE

PUBLIC :: PartitionFuncGas
PUBLIC :: PartitionFuncSurf
PUBLIC :: PartitionFuncActAdsorb
PUBLIC :: PartitionFuncActDissAdsorb
PUBLIC :: PartitionFuncActER
PUBLIC :: PartitionFuncActDesorb
PUBLIC :: PartitionFuncActDissSurf
PUBLIC :: PartitionFuncActLH
PUBLIC :: PartitionFuncActExchSurf
!===================================================================================================================================

CONTAINS


REAL FUNCTION QPartTrans(iSpec, Temp, InDim)
!===================================================================================================================================
!> Calculation of translational partition function for particle species at given temperature
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars  ,ONLY: Pi, PlanckConst, BoltzmannConst
USE MOD_Particle_Vars ,ONLY: Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
REAL, INTENT(IN)              :: Temp
INTEGER, INTENT(IN)           :: InDim
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
!===================================================================================================================================
QPartTrans = 2. * Pi * Species(iSpec)%MassIC * BoltzmannConst * Temp / (PlanckConst**2)
SELECT CASE (InDim)
CASE (1)
  QPartTrans = SQRT(QPartTrans)
CASE (2)
  QPartTrans = QPartTrans**1
CASE DEFAULT
  QPartTrans = QPartTrans**1.5
END SELECT
END FUNCTION QPartTrans


REAL FUNCTION QPartRot(iSpec, Temp)
!===================================================================================================================================
!> Calculation of rotational partition function for particle species at given temperature
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars ,ONLY: Pi, PlanckConst, BoltzmannConst
USE MOD_DSMC_Vars    ,ONLY: SpecDSMC, PolyatomMolDSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
REAL, INTENT(IN)              :: Temp
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                       :: iPolyatMole
!===================================================================================================================================
IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
  IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
    iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
    IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec) THEN
      QPartRot = Temp / (SpecDSMC(iSpec)%SymmetryFactor * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))
    ELSE
      QPartRot = SQRT(Pi) / SpecDSMC(iSpec)%SymmetryFactor * SQRT(Temp**3/( PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)    &
                                                                          * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(2)    &
                                                                          * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(3)))
    END IF
  ELSE
    QPartRot = Temp / (SpecDSMC(iSpec)%SymmetryFactor * SpecDSMC(iSpec)%CharaTRot)
  END IF
ELSE
  QPartRot = 1.
END IF
END FUNCTION QPartRot


REAL FUNCTION QPartVib(iSpec, Temp)
!===================================================================================================================================
!> Calculation of vibrational partition function for particle species at given temperature
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars ,ONLY: SpecDSMC, PolyatomMolDSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
REAL, INTENT(IN)              :: Temp
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                       :: iPolyatMole, iDOF
!===================================================================================================================================
IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
  IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
    iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
    QPartVib = 1.
    DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
      QPartVib = QPartVib / (1. - EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / Temp))
    END DO
  ELSE
    QPartVib = 1. / (1. - EXP(-SpecDSMC(iSpec)%CharaTVib / Temp))
  END IF
ELSE
  QPartVib = 1.
END IF
END FUNCTION QPartVib


REAL FUNCTION QPartElec(iSpec, Temp)
!===================================================================================================================================
!> Calculation of electronic partition function for particle species at given temperature
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars ,ONLY: SpecDSMC, DSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
REAL, INTENT(IN)              :: Temp
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                       :: iDOF
!===================================================================================================================================
IF ( DSMC%ElectronicModel ) THEN
  IF(SpecDSMC(iSpec)%InterID.NE.4) THEN
    QPartElec = 0.
    DO iDOF=0, SpecDSMC(iSpec)%MaxElecQuant - 1
      QPartElec = QPartElec + SpecDSMC(iSpec)%ElectronicState(1,iDOF) * EXP(-SpecDSMC(iSpec)%ElectronicState(2,iDOF) / Temp)
    END DO
  ELSE
    QPartElec = 1.
  END IF
ELSE
  QPartElec = 1.
END IF
END FUNCTION QPartElec


REAL FUNCTION PartitionFuncGas(iSpec, Temp)
!===================================================================================================================================
!> Calculation of partition function for gaseous state for particle species at given temperature
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
REAL, INTENT(IN)              :: Temp
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
REAL                          :: Qtra, Qrot, Qvib, Qelec
!===================================================================================================================================
Qtra = QPartTrans(iSpec,Temp,3)
Qrot = QPartRot(iSpec,Temp)
Qvib = QPartVib(iSpec,Temp)
Qelec = QPartElec(iSpec,Temp)
PartitionFuncGas = Qtra * Qrot * Qvib * Qelec
END FUNCTION PartitionFuncGas


REAL FUNCTION PartitionFuncSurf(iSpec, Temp, CharaSurfTemp)
!===================================================================================================================================
!> Calculate partition function of adsorbates on surface for certain species at given temp
!===================================================================================================================================
! MODULES
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
REAL, INTENT(IN)              :: Temp
REAL, INTENT(IN)              :: CharaSurfTemp
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
REAL                          :: Qtra, Qrot, Qvib, Qelec
!===================================================================================================================================
Qtra = 1.
Qrot = QPartRot(iSpec,Temp)
Qvib = QPartVib(iSpec,Temp) * 1./(1.-EXP(-CharaSurfTemp/Temp))
Qelec = QPartElec(iSpec,Temp)
PartitionFuncSurf = Qtra * Qrot * Qvib * Qelec
END FUNCTION PartitionFuncSurf


REAL FUNCTION PartitionFuncActAdsorb(iSpec, Temp)
!===================================================================================================================================
!> Calculation of Partitionfunction of activated complex (molecular adsorption)
!===================================================================================================================================
! MODULES
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
REAL, INTENT(IN)              :: Temp
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
REAL                          :: Qtra, Qrot, Qvib, Qelec
!===================================================================================================================================
Qtra = QPartTrans(iSpec,Temp,2)
Qrot = QPartRot(iSpec,Temp)
Qvib = QPartVib(iSpec,Temp)
Qelec = QPartElec(iSpec,Temp)
PartitionFuncActAdsorb = Qtra * Qrot * Qvib * Qelec
END FUNCTION PartitionFuncActAdsorb


REAL FUNCTION PartitionFuncActDissAdsorb(EductSpec,ResultSpec1,ResultSpec2,Temp)
!===================================================================================================================================
!> Calculation of Partitionfunction of activated complex (dissociative adsorption)
!===================================================================================================================================
! MODULES
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: EductSpec
INTEGER, INTENT(IN)           :: ResultSpec1
INTEGER, INTENT(IN)           :: ResultSpec2
REAL, INTENT(IN)              :: Temp
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
REAL                          :: Qtra, Qrot, Qvib, Qelec
!===================================================================================================================================
Qtra = QPartTrans(EductSpec,Temp,2)
Qrot = QPartRot(EductSpec,Temp)
Qvib = QPartVib(ResultSpec1,Temp) * QPartVib(ResultSpec2,Temp)
Qelec = QPartElec(EductSpec,Temp)
PartitionFuncActDissAdsorb = Qtra * Qrot * Qvib * Qelec
END FUNCTION PartitionFuncActDissAdsorb


REAL FUNCTION PartitionFuncActER(EductSpec1,EductSpec2,ResultSpec,Temp)
!===================================================================================================================================
!> Calculation of Partitionfunction of activated complex (ER reaction)
!===================================================================================================================================
! MODULES
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: EductSpec1
INTEGER, INTENT(IN)           :: EductSpec2
INTEGER, INTENT(IN)           :: ResultSpec
REAL, INTENT(IN)              :: Temp
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
REAL                          :: Qtra, Qrot, Qvib, Qelec
!===================================================================================================================================
Qtra = QPartTrans(ResultSpec,Temp,2)
Qrot = QPartRot(ResultSpec,Temp)
Qvib = QPartVib(EductSpec1,Temp) * QPartVib(EductSpec2,Temp)
Qelec = QPartElec(ResultSpec,Temp)
PartitionFuncActER = Qtra * Qrot * Qvib * Qelec
END FUNCTION PartitionFuncActER


REAL FUNCTION PartitionFuncActDesorb(iSpec,Temp,SurfDensity)
!===================================================================================================================================
!> Calculation of Partitionfunction of activated complex (desorption)
!===================================================================================================================================
! MODULES
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
REAL, INTENT(IN)              :: Temp
REAL, INTENT(IN)              :: SurfDensity
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
REAL                          :: Qtra, Qrot, Qvib, Qelec
!===================================================================================================================================
Qtra = QPartTrans(iSpec,Temp,2)/SurfDensity
!Qtra = 1.
Qrot = 1.
Qvib = QPartVib(iSpec,Temp)
Qelec = QPartElec(iSpec,Temp)
PartitionFuncActDesorb = Qtra * Qrot * Qvib * Qelec
END FUNCTION PartitionFuncActDesorb


REAL FUNCTION PartitionFuncActLH(EductSpec1,EductSpec2,ResultSpec,Temp,SurfDensity)
!===================================================================================================================================
!> Calculation of Partitionfunction of activated complex (LH reaction at surface)
!===================================================================================================================================
! MODULES
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: EductSpec1
INTEGER, INTENT(IN)           :: EductSpec2
INTEGER, INTENT(IN)           :: ResultSpec
REAL, INTENT(IN)              :: Temp
REAL, INTENT(IN)              :: SurfDensity
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
REAL                          :: Qtra, Qrot, Qvib, Qelec
!===================================================================================================================================
Qtra = QPartTrans(ResultSpec,Temp,2)/SurfDensity
Qrot = 1.
Qvib = QPartVib(EductSpec1,Temp) * QPartVib(EductSpec2,Temp)
Qelec = QPartElec(ResultSpec,Temp)
PartitionFuncActLH = Qtra * Qrot * Qvib * Qelec
END FUNCTION PartitionFuncActLH


REAL FUNCTION PartitionFuncActDissSurf(EductSpec,ResultSpec1,ResultSpec2,Temp,SurfDensity)
!===================================================================================================================================
!> Calculation of Partitionfunction of activated complex (dissociation reaction at surface)
!===================================================================================================================================
! MODULES
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: EductSpec
INTEGER, INTENT(IN)           :: ResultSpec1
INTEGER, INTENT(IN)           :: ResultSpec2
REAL, INTENT(IN)              :: Temp
REAL, INTENT(IN)              :: SurfDensity
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
REAL                          :: Qtra, Qrot, Qvib, Qelec
!===================================================================================================================================
Qtra = QPartTrans(EductSpec,Temp,2)/SurfDensity
Qrot = 1.
Qvib = QPartVib(ResultSpec1,Temp) * QPartVib(ResultSpec2,Temp)
Qelec = QPartElec(EductSpec,Temp)
PartitionFuncActDissSurf = Qtra * Qrot * Qvib * Qelec
END FUNCTION PartitionFuncActDissSurf


REAL FUNCTION PartitionFuncActExchSurf(EductSpec1,EductSpec2,ResultSpec1,ResultSpec2,Temp,SurfDensity)
!===================================================================================================================================
!> Calculation of Partitionfunction of activated complex (dissociation reaction at surface)
!===================================================================================================================================
! MODULES
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: EductSpec1
INTEGER, INTENT(IN)           :: EductSpec2
INTEGER, INTENT(IN)           :: ResultSpec1
INTEGER, INTENT(IN)           :: ResultSpec2
REAL, INTENT(IN)              :: Temp
REAL, INTENT(IN)              :: SurfDensity
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
REAL                          :: Qtra, Qrot, Qvib, Qelec
!===================================================================================================================================
Qtra = QPartTrans(EductSpec1,Temp,2)/SurfDensity * QPartTrans(EductSpec2,Temp,2)/SurfDensity
Qrot = 1.
Qvib = QPartVib(ResultSpec1,Temp) * QPartVib(ResultSpec2,Temp)
Qelec = QPartElec(ResultSpec1,Temp) * QPartElec(ResultSpec2,Temp)
PartitionFuncActExchSurf = Qtra * Qrot * Qvib * Qelec
END FUNCTION PartitionFuncActExchSurf


END MODULE MOD_SurfaceModel_PartFunc


