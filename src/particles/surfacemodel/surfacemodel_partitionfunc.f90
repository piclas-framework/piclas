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

INTERFACE PartitionFuncActAdsorb
  MODULE PROCEDURE PartitionFuncActAdsorb
END INTERFACE

INTERFACE PartitionFuncActER
  MODULE PROCEDURE PartitionFuncActER
END INTERFACE

INTERFACE PartitionFuncAct
  MODULE PROCEDURE PartitionFuncAct
END INTERFACE

INTERFACE PartitionFuncAct_dissoc
  MODULE PROCEDURE PartitionFuncAct_dissoc
END INTERFACE

INTERFACE PartitionFuncAct_recomb
  MODULE PROCEDURE PartitionFuncAct_recomb
END INTERFACE

INTERFACE PartitionFuncAct_exch
  MODULE PROCEDURE PartitionFuncAct_exch
END INTERFACE

INTERFACE PartitionFuncSurf
  MODULE PROCEDURE PartitionFuncSurf
END INTERFACE

PUBLIC :: PartitionFuncGas
PUBLIC :: PartitionFuncActAdsorb
PUBLIC :: PartitionFuncActER
PUBLIC :: PartitionFuncAct
PUBLIC :: PartitionFuncAct_dissoc
PUBLIC :: PartitionFuncAct_recomb
PUBLIC :: PartitionFuncAct_exch
PUBLIC :: PartitionFuncSurf
!===================================================================================================================================

CONTAINS

SUBROUTINE PartitionFuncGas(iSpec, Temp, VarPartitionFuncGas)
!===================================================================================================================================
!> Calculation of partition function for gaseous state particle species at given temperature
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars  ,ONLY: Pi, PlanckConst, BoltzmannConst
USE MOD_DSMC_Vars     ,ONLY: SpecDSMC, PolyatomMolDSMC, DSMC
USE MOD_Particle_Vars ,ONLY: Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
REAL, INTENT(IN)              :: Temp
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: VarPartitionFuncGas
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                       :: iPolyatMole, iDOF
REAL                          :: Qtra, Qrot, Qvib, Qelec
!===================================================================================================================================

Qtra = (2. * Pi * Species(iSpec)%MassIC * BoltzmannConst * Temp / (PlanckConst**2))**(1.5)
IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
  IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
    iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
    IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec) THEN
      Qrot = Temp / (SpecDSMC(iSpec)%SymmetryFactor * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))
    ELSE
      Qrot = SQRT(Pi) / SpecDSMC(iSpec)%SymmetryFactor * SQRT(Temp**3/( PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)    &
                                                                      * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(2)    &
                                                                      * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(3)))
    END IF
    Qvib = 1.
    DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
      Qvib = Qvib / (1. - EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / Temp))
    END DO
  ELSE
    Qrot = Temp / (SpecDSMC(iSpec)%SymmetryFactor * SpecDSMC(iSpec)%CharaTRot)
    Qvib = 1. / (1. - EXP(-SpecDSMC(iSpec)%CharaTVib / Temp))
  END IF
ELSE
  Qrot = 1.
  Qvib = 1.
END IF
IF ( DSMC%ElectronicModel ) THEN
  IF(SpecDSMC(iSpec)%InterID.NE.4) THEN
    Qelec = 0.
    DO iDOF=0, SpecDSMC(iSpec)%MaxElecQuant - 1
      Qelec = Qelec + SpecDSMC(iSpec)%ElectronicState(1,iDOF) * EXP(-SpecDSMC(iSpec)%ElectronicState(2,iDOF) / Temp)
    END DO
  ELSE
    Qelec = 1.
  END IF
  VarPartitionFuncGas = Qtra * Qrot * Qvib * Qelec
ELSE
  VarPartitionFuncGas = Qtra * Qrot * Qvib
END IF

END SUBROUTINE PartitionFuncGas


SUBROUTINE PartitionFuncActAdsorb(iSpec, Temp, VarPartitionFuncAct, Surfdensity)
!===================================================================================================================================
!> Calculation of Partitionfunction of activated complex (molecular adsorption)
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars  ,ONLY: Pi, PlanckConst, BoltzmannConst
USE MOD_DSMC_Vars     ,ONLY: SpecDSMC, PolyatomMolDSMC, DSMC
USE MOD_Particle_Vars ,ONLY: Species
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
REAL, INTENT(IN)              :: Temp
REAL, INTENT(IN)              :: Surfdensity
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: VarPartitionFuncAct
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                       :: iPolyatMole, iDOF
REAL                          :: Qtra, Qrot, Qvib, Qelec
!===================================================================================================================================
Qtra = ((2. * Pi * Species(iSpec)%MassIC * BoltzmannConst * Temp / (PlanckConst**2))/(Surfdensity))
IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
  IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
    iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
       IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec) THEN
         Qrot = Temp / (SpecDSMC(iSpec)%SymmetryFactor * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))
       ELSE
         Qrot = SQRT(Pi) / SpecDSMC(iSpec)%SymmetryFactor * SQRT(Temp**3/( PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)    &
                                                                         * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(2)    &
                                                                         * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(3)))
       END IF
    Qrot = 1.
    Qvib = 1.
    DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
      Qvib = Qvib * 1. / (1. - EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / Temp))
    END DO
  ELSE
       Qrot = Temp / (SpecDSMC(iSpec)%SymmetryFactor * SpecDSMC(iSpec)%CharaTRot)
    Qrot = 1.
    Qvib = 1. / (1. - EXP(-SpecDSMC(iSpec)%CharaTVib / Temp))
  END IF
ELSE
  Qrot = 1.
  Qvib = 1.
END IF
IF ( DSMC%ElectronicModel ) THEN
  IF(SpecDSMC(iSpec)%InterID.NE.4) THEN
    Qelec = 0.
    DO iDOF=0, SpecDSMC(iSpec)%MaxElecQuant - 1
      Qelec = Qelec + SpecDSMC(iSpec)%ElectronicState(1,iDOF) * EXP(-SpecDSMC(iSpec)%ElectronicState(2,iDOF) / Temp)
    END DO
  ELSE
    Qelec = 1.
  END IF
  VarPartitionFuncAct = Qtra * Qrot * Qvib * Qelec
ELSE
  VarPartitionFuncAct = Qtra * Qrot * Qvib
END IF

END SUBROUTINE PartitionFuncActAdsorb


SUBROUTINE PartitionFuncActER(iSpec, Temp, VarPartitionFuncAct, Surfdensity)
!===================================================================================================================================
!> Calculation of Partitionfunction of activated complex (ER reaction)
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars  ,ONLY: Pi, PlanckConst, BoltzmannConst
USE MOD_DSMC_Vars     ,ONLY: SpecDSMC, PolyatomMolDSMC, DSMC
USE MOD_Particle_Vars ,ONLY: Species
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
REAL, INTENT(IN)              :: Temp
REAL, INTENT(IN)              :: Surfdensity
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: VarPartitionFuncAct
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                       :: iPolyatMole, iDOF
REAL                          :: Qtra, Qrot, Qvib, Qelec
!===================================================================================================================================
Qtra = ((2. * Pi * Species(iSpec)%MassIC * BoltzmannConst * Temp / (PlanckConst**2))/(Surfdensity))
IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
  IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
    iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
!       IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec) THEN
!         Qrot = Temp / (SpecDSMC(iSpec)%SymmetryFactor * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))
!       ELSE
!         Qrot = SQRT(Pi) / SpecDSMC(iSpec)%SymmetryFactor * SQRT(Temp**3/( PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)    &
!                                                                         * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(2)    &
!                                                                         * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(3)))
!       END IF
    Qrot = 1.
    Qvib = 1.
    DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
      Qvib = Qvib * 1. / (1. - EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / Temp))
    END DO
  ELSE
!       Qrot = Temp / (SpecDSMC(iSpec)%SymmetryFactor * SpecDSMC(iSpec)%CharaTRot)
    Qrot = 1.
    Qvib = 1. / (1. - EXP(-SpecDSMC(iSpec)%CharaTVib / Temp))
  END IF
ELSE
  Qrot = 1.
  Qvib = 1.
END IF
IF ( DSMC%ElectronicModel ) THEN
  IF(SpecDSMC(iSpec)%InterID.NE.4) THEN
    Qelec = 0.
    DO iDOF=0, SpecDSMC(iSpec)%MaxElecQuant - 1
      Qelec = Qelec + SpecDSMC(iSpec)%ElectronicState(1,iDOF) * EXP(-SpecDSMC(iSpec)%ElectronicState(2,iDOF) / Temp)
    END DO
  ELSE
    Qelec = 1.
  END IF
  VarPartitionFuncAct = Qtra * Qrot * Qvib * Qelec
ELSE
  VarPartitionFuncAct = Qtra * Qrot * Qvib
END IF

END SUBROUTINE PartitionFuncActER


SUBROUTINE PartitionFuncAct(iSpec, Temp, VarPartitionFuncAct)!, Surfdensity)
!===================================================================================================================================
!> Calculation of Partitionfunction of activated complex (molecular desorption)
!===================================================================================================================================
! MODULES
!USE MOD_Globals_Vars ,ONLY: Pi, PlanckConst, BoltzmannConst
USE MOD_DSMC_Vars    ,ONLY: SpecDSMC, PolyatomMolDSMC, DSMC
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
REAL, INTENT(IN)              :: Temp
!REAL, INTENT(IN)              :: Surfdensity
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: VarPartitionFuncAct
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                       :: iPolyatMole, iDOF
REAL                          :: Qtra, Qrot, Qvib, Qelec
!===================================================================================================================================
Qtra = 1.!((2. * Pi * Species(iSpec)%MassIC * BoltzmannConst * Temp / (PlanckConst**2))/(Surfdensity))
IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
  IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
    iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
!       IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec) THEN
!         Qrot = Temp / (SpecDSMC(iSpec)%SymmetryFactor * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))
!       ELSE
!         Qrot = SQRT(Pi) / SpecDSMC(iSpec)%SymmetryFactor * SQRT(Temp**3/( PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)    &
!                                                                         * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(2)    &
!                                                                         * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(3)))
!       END IF
    Qrot = 1.
    Qvib = 1.
    DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
      Qvib = Qvib * 1. / (1. - EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / Temp))
    END DO
  ELSE
!       Qrot = Temp / (SpecDSMC(iSpec)%SymmetryFactor * SpecDSMC(iSpec)%CharaTRot)
    Qrot = 1.
    Qvib = 1. / (1. - EXP(-SpecDSMC(iSpec)%CharaTVib / Temp))
  END IF
ELSE
  Qrot = 1.
  Qvib = 1.
END IF
IF ( DSMC%ElectronicModel ) THEN
  IF(SpecDSMC(iSpec)%InterID.NE.4) THEN
    Qelec = 0.
    DO iDOF=0, SpecDSMC(iSpec)%MaxElecQuant - 1
      Qelec = Qelec + SpecDSMC(iSpec)%ElectronicState(1,iDOF) * EXP(-SpecDSMC(iSpec)%ElectronicState(2,iDOF) / Temp)
    END DO
  ELSE
    Qelec = 1.
  END IF
  VarPartitionFuncAct = Qtra * Qrot * Qvib * Qelec
ELSE
  VarPartitionFuncAct = Qtra * Qrot * Qvib
END IF

END SUBROUTINE PartitionFuncAct


SUBROUTINE PartitionFuncAct_dissoc(iSpec,Prod_Spec1,Prod_Spec2, Temp, VarPartitionFuncAct, Surfdensity)
!SUBROUTINE PartitionFuncAct_dissoc(Prod_Spec1,Prod_Spec2, Temp, VarPartitionFuncAct)
!===================================================================================================================================
!> Calculation of Partitionfunction of activated complex (dissociation at surface)
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars ,ONLY: Pi, PlanckConst, BoltzmannConst
USE MOD_DSMC_Vars    ,ONLY: SpecDSMC, PolyatomMolDSMC, DSMC
USE MOD_Particle_Vars,ONLY: Species
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
INTEGER, INTENT(IN)           :: Prod_Spec1
INTEGER, INTENT(IN)           :: Prod_Spec2
REAL, INTENT(IN)              :: Temp
REAL, INTENT(IN)              :: Surfdensity
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: VarPartitionFuncAct
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                       :: iPolyatMole, iDOF
REAL                          :: Qtra, Qrot, Qvib, Qelec1, Qelec2, Qelec
!===================================================================================================================================
Qtra = ((2. * Pi * Species(iSpec)%MassIC * BoltzmannConst * Temp / (PlanckConst**2))/(Surfdensity))
Qrot = 1.
Qvib = 1.
Qelec1 = 1.
Qelec2 = 1.
Qelec = 1.
IF(SpecDSMC(Prod_Spec1)%InterID.EQ.2) THEN
  IF(SpecDSMC(Prod_Spec1)%PolyatomicMol) THEN
    iPolyatMole = SpecDSMC(Prod_Spec1)%SpecToPolyArray
    DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
      Qvib = Qvib * 1. / (1. - EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / Temp))
    END DO
  ELSE
    Qvib = Qvib * 1. / (1. - EXP(-SpecDSMC(Prod_Spec1)%CharaTVib / Temp))
  END IF
END IF
IF(SpecDSMC(Prod_Spec2)%InterID.EQ.2) THEN
  IF(SpecDSMC(Prod_Spec2)%PolyatomicMol) THEN
    iPolyatMole = SpecDSMC(Prod_Spec2)%SpecToPolyArray
    DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
      Qvib = Qvib * 1. / (1. - EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / Temp))
    END DO
  ELSE
    Qvib = Qvib * 1. / (1. - EXP(-SpecDSMC(Prod_Spec2)%CharaTVib / Temp))
  END IF
END IF
IF ( DSMC%ElectronicModel ) THEN
  IF(SpecDSMC(Prod_Spec1)%InterID.NE.4) THEN
    Qelec1 = 0.
    DO iDOF=0, SpecDSMC(Prod_Spec1)%MaxElecQuant - 1
      Qelec1 = Qelec1 + SpecDSMC(Prod_Spec1)%ElectronicState(1,iDOF) * EXP(-SpecDSMC(Prod_Spec1)%ElectronicState(2,iDOF) / Temp)
    END DO
  END IF
  IF(SpecDSMC(Prod_Spec2)%InterID.NE.4) THEN
    Qelec2 = 0.
    DO iDOF=0, SpecDSMC(Prod_Spec2)%MaxElecQuant - 1
      Qelec2 = Qelec2 + SpecDSMC(Prod_Spec2)%ElectronicState(1,iDOF) * EXP(-SpecDSMC(Prod_Spec2)%ElectronicState(2,iDOF) / Temp)
    END DO
  END IF
  Qelec=Qelec1+Qelec2
END IF
VarPartitionFuncAct = Qtra * Qrot * Qvib * Qelec

END SUBROUTINE PartitionFuncAct_dissoc


SUBROUTINE PartitionFuncAct_recomb(Educt_Spec1, Educt_Spec2,Product_Spec, Temp, VarPartitionFuncAct, Surfdensity)
!===================================================================================================================================
!> Calculation of Partitionfunction of activated complex (recombination for desorption)
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars  ,ONLY: Pi, PlanckConst, BoltzmannConst
USE MOD_DSMC_Vars     ,ONLY: SpecDSMC, PolyatomMolDSMC, DSMC
USE MOD_Particle_Vars ,ONLY: Species
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: Educt_Spec1
INTEGER, INTENT(IN)           :: Educt_Spec2
INTEGER, INTENT(IN)           :: Product_Spec
REAL, INTENT(IN)              :: Temp
REAL, INTENT(IN)              :: Surfdensity
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: VarPartitionFuncAct
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                       :: iPolyatMole, iDOF
REAL                          :: Qtra, Qrot, Qvib, Qelec, Qelec1, Qelec2
!===================================================================================================================================
Qtra = ((2. * Pi * Species(Product_Spec)%MassIC * BoltzmannConst * Temp / (PlanckConst**2))/(Surfdensity))
!Qtra = ((2. * Pi * Species(Educt_Spec1)%MassIC * BoltzmannConst * Temp / (PlanckConst**2))/(Surfdensity))**0.5 &
!     * ((2. * Pi * Species(Educt_Spec2)%MassIC * BoltzmannConst * Temp / (PlanckConst**2))/(Surfdensity))**0.5
Qrot = 1.
Qvib = 1.
Qelec1 = 1.
Qelec2 = 1.
Qelec = 1.
IF(SpecDSMC(Educt_Spec1)%InterID.EQ.2) THEN
  IF(SpecDSMC(Educt_Spec1)%PolyatomicMol) THEN
    iPolyatMole = SpecDSMC(Educt_Spec1)%SpecToPolyArray
    DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
      Qvib = Qvib * 1. / (1. - EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / Temp))
    END DO
  ELSE
    Qvib = Qvib * 1. / (1. - EXP(-SpecDSMC(Educt_Spec1)%CharaTVib / Temp))
  END IF
END IF
IF(SpecDSMC(Educt_Spec2)%InterID.EQ.2) THEN
  IF(SpecDSMC(Educt_Spec2)%PolyatomicMol) THEN
    iPolyatMole = SpecDSMC(Educt_Spec2)%SpecToPolyArray
    DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
      Qvib = Qvib * 1. / (1. - EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / Temp))
    END DO
  ELSE
    Qvib = Qvib * 1. / (1. - EXP(-SpecDSMC(Educt_Spec2)%CharaTVib / Temp))
  END IF
END IF
IF ( DSMC%ElectronicModel ) THEN
  IF(SpecDSMC(Educt_Spec1)%InterID.NE.4) THEN
    Qelec1 = 0.
    DO iDOF=0, SpecDSMC(Educt_Spec1)%MaxElecQuant - 1
      Qelec1 = Qelec1 + SpecDSMC(Educt_Spec1)%ElectronicState(1,iDOF) * EXP(-SpecDSMC(Educt_Spec1)%ElectronicState(2,iDOF) / Temp)
    END DO
  END IF
  IF(SpecDSMC(Educt_Spec2)%InterID.NE.4) THEN
    Qelec2 = 0.
    DO iDOF=0, SpecDSMC(Educt_Spec2)%MaxElecQuant - 1
      Qelec2 = Qelec2 + SpecDSMC(Educt_Spec2)%ElectronicState(1,iDOF) * EXP(-SpecDSMC(Educt_Spec2)%ElectronicState(2,iDOF) / Temp)
    END DO
  END IF
  Qelec=Qelec1+Qelec2
END IF
VarPartitionFuncAct = Qtra * Qrot * Qvib * Qelec

END SUBROUTINE PartitionFuncAct_recomb


SUBROUTINE PartitionFuncAct_exch(Educt_Spec1, Educt_Spec2, Temp, VarPartitionFuncAct, Surfdensity)
!===================================================================================================================================
!> Calculate Partitionfunction of activated complex (exchange reactions at surface)
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars  ,ONLY: PlanckConst, BoltzmannConst, Pi
USE MOD_DSMC_Vars     ,ONLY: SpecDSMC, PolyatomMolDSMC
USE MOD_Particle_Vars ,ONLY: Species
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: Educt_Spec1
INTEGER, INTENT(IN)           :: Educt_Spec2
REAL, INTENT(IN)              :: Temp
REAL, INTENT(IN)              :: Surfdensity
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: VarPartitionFuncAct
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                       :: iPolyatMole, iDOF
REAL                          :: Qtra, Qrot, Qvib, Qelec
!===================================================================================================================================
!   Qtra = ((2. * Pi * Species(Result_Spec)%MassIC * BoltzmannConst * Temp / (PlanckConst**2))/(Surfdensity))**0.5
Qtra = ((2. * Pi * Species(Educt_Spec1)%MassIC * BoltzmannConst * Temp / (PlanckConst**2))/(Surfdensity))**0.5 &
     * ((2. * Pi * Species(Educt_Spec2)%MassIC * BoltzmannConst * Temp / (PlanckConst**2))/(Surfdensity))**0.5
Qrot = 1.
Qvib = 1.
IF(SpecDSMC(Educt_Spec1)%InterID.EQ.2) THEN
  IF(SpecDSMC(Educt_Spec1)%PolyatomicMol) THEN
    iPolyatMole = SpecDSMC(Educt_Spec1)%SpecToPolyArray
    Qrot = 1.
    DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
      Qvib = Qvib * 1. / (1. - EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / Temp))
    END DO
  ELSE
    Qrot = 1.
    Qvib = Qvib * 1. / (1. - EXP(-SpecDSMC(Educt_Spec1)%CharaTVib / Temp))
  END IF
END IF
IF(SpecDSMC(Educt_Spec2)%InterID.EQ.2) THEN
  IF(SpecDSMC(Educt_Spec2)%PolyatomicMol) THEN
    iPolyatMole = SpecDSMC(Educt_Spec2)%SpecToPolyArray
    Qrot = 1.
    DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
      Qvib = Qvib * 1. / (1. - EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / Temp))
    END DO
  ELSE
    Qrot = 1.
    Qvib = Qvib * 1. / (1. - EXP(-SpecDSMC(Educt_Spec2)%CharaTVib / Temp))
  END IF
END IF
!IF(SpecDSMC(iSpec)%InterID.NE.4) THEN
!  Qelec = 0.
!  DO iDOF=0, SpecDSMC(iSpec)%MaxElecQuant - 1
!    Qelec = Qelec + SpecDSMC(iSpec)%ElectronicState(1,iDOF) * EXP(-SpecDSMC(iSpec)%ElectronicState(2,iDOF) / Temp)
!  END DO
!ELSE
!  Qelec = 1.
!END IF
VarPartitionFuncAct = Qtra * Qrot * Qvib! * Qelec

END SUBROUTINE PartitionFuncAct_exch


SUBROUTINE PartitionFuncSurf(iSpec, Temp, VarPartitionFuncSurf, CharaTemp, PartBoundID)
!===================================================================================================================================
!> Calculate partition function of adsorbates on surface for certain species at given partboundID
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars      ,ONLY: PlanckConst, BoltzmannConst
USE MOD_DSMC_Vars         ,ONLY: SpecDSMC, PolyatomMolDSMC, DSMC
!USE MOD_SurfaceModel_Vars ,ONLY: Adsorption
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
REAL, INTENT(IN)              :: Temp
REAL, INTENT(IN)              :: CharaTemp
INTEGER, INTENT(IN)           :: PartBoundID
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: VarPartitionFuncSurf
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                       :: iPolyatMole, iDOF
REAL                          :: Qtra, Qrot, Qvib, Qelec
!===================================================================================================================================
Qtra = 1.
IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
  IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
    iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
    Qvib = 1.
    DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
      Qvib = Qvib * 1. / (1. - EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / Temp))
    END DO
  ELSE
    Qvib = 1. / (1. - EXP(-SpecDSMC(iSpec)%CharaTVib / Temp))
  END IF
ELSE
  Qvib = 1.
END IF
!IF (Adsorption%Coordination(PartBoundID,iSpec).EQ.1) THEN
!Qvib = Qvib * 1. / (1. - EXP(-CharaTemp / Temp))
!ELSE IF (Adsorption%Coordination(PartBoundID,iSpec).EQ.2) THEN
!  IF ((Adsorption%DiCoord(PartBoundID,iSpec).EQ.1) .OR. (Adsorption%DiCoord(PartBoundID,iSpec).EQ.2)) THEN
!    Qvib = Qvib * 1. / (1. - EXP(-CharaTemp / Temp))!**(2. - 1./2.)
!  ELSE
!    Qvib = Qvib * 1. / (1. - EXP(-CharaTemp / Temp))
!  END IF
!ELSE
!Qvib = Qvib * 1. / (1. - EXP(-CharaTemp / Temp))  ! (2-1./REAL(Adsorption%CrystalIndx(1)))
!END IF
Qvib = Qvib * 1. / (1. - EXP(-CharaTemp / Temp))
Qrot = 1.
IF ( DSMC%ElectronicModel ) THEN
  IF(SpecDSMC(iSpec)%InterID.NE.4) THEN
    Qelec = 0.
    DO iDOF=0, SpecDSMC(iSpec)%MaxElecQuant - 1
      Qelec = Qelec + SpecDSMC(iSpec)%ElectronicState(1,iDOF) * EXP(-SpecDSMC(iSpec)%ElectronicState(2,iDOF) / Temp)
    END DO
  ELSE
    Qelec = 1.
  END IF
  VarPartitionFuncSurf = Qtra * Qrot * Qvib * Qelec
ELSE
  VarPartitionFuncSurf = Qtra * Qrot * Qvib
END IF

END SUBROUTINE PartitionFuncSurf


END MODULE MOD_SurfaceModel_PartFunc


