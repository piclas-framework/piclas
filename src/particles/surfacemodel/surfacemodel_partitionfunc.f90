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
!INTERFACE PartitionFuncGas
!  MODULE PROCEDURE PartitionFuncGas
!END INTERFACE

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

INTERFACE AnalyzePartitionTemp
  MODULE PROCEDURE AnalyzePartitionTemp
END INTERFACE

!PUBLIC :: PartitionFuncGas
PUBLIC :: PartitionFuncAct
PUBLIC :: PartitionFuncAct_dissoc
PUBLIC :: PartitionFuncAct_recomb
PUBLIC :: PartitionFuncAct_exch
PUBLIC :: PartitionFuncSurf
PUBLIC :: AnalyzePartitionTemp
!===================================================================================================================================

CONTAINS

!SUBROUTINE PartitionFuncGas(iSpec, Temp, VarPartitionFuncGas)
!!===================================================================================================================================
!!> Calculation of partition function for gaseous state particle species at certain temperature
!!===================================================================================================================================
!! MODULES
!USE MOD_Globals
!USE MOD_Globals_Vars  ,ONLY: PlanckConst, BoltzmannConst
!USE MOD_DSMC_Vars     ,ONLY: SpecDSMC, PolyatomMolDSMC
!USE MOD_Particle_Vars ,ONLY: Species
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!!----------------------------------------------------------------------------------------------------------------------------------!
!! INPUT VARIABLES
!INTEGER, INTENT(IN)           :: iSpec
!REAL, INTENT(IN)              :: Temp
!!----------------------------------------------------------------------------------------------------------------------------------!
!! OUTPUT VARIABLES
!REAL, INTENT(OUT)             :: VarPartitionFuncGas
!!----------------------------------------------------------------------------------------------------------------------------------!
!! OUTPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------!
!! LOCAL VARIABLES
!REAL, PARAMETER               :: Pi=3.14159265358979323846_8
!INTEGER                       :: iPolyatMole, iDOF
!REAL                          :: Qtra, Qrot, Qvib!, Qelec
!!===================================================================================================================================
!Qtra = (2. * Pi * Species(iSpec)%MassIC * BoltzmannConst * Temp / (PlanckConst**2))**(1.5)
!IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
!  IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
!    iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
!    IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec) THEN
!      Qrot = Temp / (SpecDSMC(iSpec)%SymmetryFactor * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))
!    ELSE
!      Qrot = SQRT(Pi) / SpecDSMC(iSpec)%SymmetryFactor * SQRT(Temp**3/( PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)    &
!                                                                      * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(2)    &
!                                                                      * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(3)))
!    END IF
!    Qvib = 1.
!    DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
!      Qvib = Qvib * 1. / (1. - EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / Temp))
!    END DO
!  ELSE
!    Qrot = Temp / (SpecDSMC(iSpec)%SymmetryFactor * SpecDSMC(iSpec)%CharaTRot)
!    Qvib = 1. / (1. - EXP(-SpecDSMC(iSpec)%CharaTVib / Temp))
!  END IF
!ELSE
!  Qrot = 1.
!  Qvib = 1.
!END IF
!!   Qelec = 0.
!!   DO iDOF=1, SpecDSMC(iSpec)%NumElecLevels
!!     Qelec = Qelec + SpecDSMC(iSpec)%ElectronicState(1,iDOF) * EXP(-SpecDSMC(iSpec)%ElectronicState(2,iDOF) / Temp)
!!   END DO
!VarPartitionFuncGas = Qtra * Qrot * Qvib! * Qelec
!
!END SUBROUTINE PartitionFuncGas


SUBROUTINE PartitionFuncAct(iSpec, Temp, VarPartitionFuncAct)!, Surfdensity)
!===================================================================================================================================
!> Calculation of Partitionfunction of activated complex (molecular desorption)
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars ,ONLY: PlanckConst, BoltzmannConst
USE MOD_DSMC_Vars    ,ONLY: SpecDSMC, PolyatomMolDSMC
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
REAL, PARAMETER               :: Pi=3.14159265358979323846_8
INTEGER                       :: iPolyatMole, iDOF
REAL                          :: Qtra, Qrot, Qvib
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
VarPartitionFuncAct = Qtra * Qrot * Qvib

END SUBROUTINE PartitionFuncAct


!SUBROUTINE PartitionFuncAct_dissoc(iSpec,Prod_Spec1,Prod_Spec2, Temp, VarPartitionFuncAct, Surfdensity)
SUBROUTINE PartitionFuncAct_dissoc(Prod_Spec1,Prod_Spec2, Temp, VarPartitionFuncAct)
!===================================================================================================================================
!> Calculation of Partitionfunction of activated complex (dissociation at surface)
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars ,ONLY: PlanckConst, BoltzmannConst
USE MOD_DSMC_Vars    ,ONLY: SpecDSMC, PolyatomMolDSMC
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!INTEGER, INTENT(IN)           :: iSpec
INTEGER, INTENT(IN)           :: Prod_Spec1
INTEGER, INTENT(IN)           :: Prod_Spec2
REAL, INTENT(IN)              :: Temp
!REAL, INTENT(IN)              :: Surfdensity
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: VarPartitionFuncAct
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
REAL, PARAMETER               :: Pi=3.14159265358979323846_8
INTEGER                       :: iPolyatMole, iDOF
REAL                          :: Qtra, Qrot, Qvib
!===================================================================================================================================
Qtra = 1.!((2. * Pi * Species(iSpec)%MassIC * BoltzmannConst * Temp / (PlanckConst**2))/(Surfdensity))
Qrot = 1.
Qvib = 1.
IF(SpecDSMC(Prod_Spec1)%InterID.EQ.2) THEN
  IF(SpecDSMC(Prod_Spec1)%PolyatomicMol) THEN
    iPolyatMole = SpecDSMC(Prod_Spec1)%SpecToPolyArray
    Qrot = 1.
    DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
      Qvib = Qvib * 1. / (1. - EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / Temp))
    END DO
  ELSE
    Qrot = 1.
    Qvib = Qvib * 1. / (1. - EXP(-SpecDSMC(Prod_Spec1)%CharaTVib / Temp))
  END IF
END IF
IF(SpecDSMC(Prod_Spec2)%InterID.EQ.2) THEN
  IF(SpecDSMC(Prod_Spec2)%PolyatomicMol) THEN
    iPolyatMole = SpecDSMC(Prod_Spec2)%SpecToPolyArray
    Qrot = 1.
    DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
      Qvib = Qvib * 1. / (1. - EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / Temp))
    END DO
  ELSE
    Qrot = 1.
    Qvib = Qvib * 1. / (1. - EXP(-SpecDSMC(Prod_Spec2)%CharaTVib / Temp))
  END IF
END IF
VarPartitionFuncAct = Qtra * Qrot * Qvib

END SUBROUTINE PartitionFuncAct_dissoc


SUBROUTINE PartitionFuncAct_recomb(Educt_Spec1, Educt_Spec2, Temp, VarPartitionFuncAct, Surfdensity)
!===================================================================================================================================
!> Calculation of Partitionfunction of activated complex (recombination for desorption)
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars  ,ONLY: PlanckConst, BoltzmannConst
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
REAL, PARAMETER               :: Pi=3.14159265358979323846_8
INTEGER                       :: iPolyatMole, iDOF
REAL                          :: Qtra, Qrot, Qvib
!===================================================================================================================================
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
VarPartitionFuncAct = Qtra * Qrot * Qvib

END SUBROUTINE PartitionFuncAct_recomb


SUBROUTINE PartitionFuncAct_exch(Educt_Spec1, Educt_Spec2, Temp, VarPartitionFuncAct, Surfdensity)
!===================================================================================================================================
!> Calculate Partitionfunction of activated complex (exchange reactions at surface)
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars  ,ONLY: PlanckConst, BoltzmannConst
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
REAL, PARAMETER               :: Pi=3.14159265358979323846_8
INTEGER                       :: iPolyatMole, iDOF
REAL                          :: Qtra, Qrot, Qvib
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
VarPartitionFuncAct = Qtra * Qrot * Qvib

END SUBROUTINE PartitionFuncAct_exch


SUBROUTINE PartitionFuncSurf(iSpec, Temp, VarPartitionFuncSurf, CharaTemp, PartBoundID)
!===================================================================================================================================
!> Calculate partition function of adsorbates on surface for certain species at given partboundID
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars      ,ONLY: PlanckConst, BoltzmannConst
USE MOD_DSMC_Vars         ,ONLY: SpecDSMC, PolyatomMolDSMC
USE MOD_SurfaceModel_Vars ,ONLY: Adsorption
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
REAL, PARAMETER               :: Pi=3.14159265358979323846_8
INTEGER                       :: iPolyatMole, iDOF
REAL                          :: Qtra, Qrot, Qvib
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
IF (Adsorption%Coordination(PartBoundID,iSpec).EQ.1) THEN
Qvib = Qvib * 1. / (1. - EXP(-CharaTemp / Temp))
ELSE IF (Adsorption%Coordination(PartBoundID,iSpec).EQ.2) THEN
  IF ((Adsorption%DiCoord(PartBoundID,iSpec).EQ.1) .OR. (Adsorption%DiCoord(PartBoundID,iSpec).EQ.2)) THEN
    Qvib = Qvib * 1. / (1. - EXP(-CharaTemp / Temp))** 2.
  ELSE
    Qvib = Qvib * 1. / (1. - EXP(-CharaTemp / Temp))
  END IF
ELSE
Qvib = Qvib * 1. / (1. - EXP(-CharaTemp / Temp))
END IF
Qrot = 1.
VarPartitionFuncSurf = Qtra * Qrot * Qvib

END SUBROUTINE PartitionFuncSurf


SUBROUTINE AnalyzePartitionTemp()
!===================================================================================================================================
!> Sampling of variables (part-density, velocity and energy) for Partition function per element
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: BoltzmannConst
USE MOD_SurfaceModel_Vars  ,ONLY: Adsorption
USE MOD_Particle_Vars      ,ONLY: PartState, PDM, PartSpecies, Species, nSpecies, PEM
USE MOD_Particle_Mesh_Vars ,ONLY: IsTracingBCElem
USE MOD_Mesh_Vars          ,ONLY: nElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: ElemID, iElem, i, iSpec
REAL , ALLOCATABLE :: Source(:,:,:)
REAL               :: TempDirec(1:3)
!===================================================================================================================================

ALLOCATE(Source(1:7,1:nElems,1:nSpecies))
Source=0.0

DO i=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(i)) THEN
    ElemID = PEM%Element(i)
    IF(.NOT.IsTracingBCElem(ElemID))CYCLE
    iSpec = PartSpecies(i)
    Source(1:3,ElemID, iSpec) = Source(1:3,ElemID,iSpec) + PartState(i,4:6)
    Source(4:6,ElemID, iSpec) = Source(4:6,ElemID,iSpec) + PartState(i,4:6)**2
    Source(7,ElemID, iSpec) = Source(7,ElemID, iSpec) + 1.0  !density
  END IF
END DO

DO iSpec=1,nSpecies
  DO iElem = 1,nElems
    IF (Source(7,iElem,iSpec).EQ.0.0) THEN
      TempDirec(1:3) = 0.0
    ELSE
      TempDirec(1:3) = Species(iSpec)%MassIC * (Source(1:3,iElem,iSpec)/Source(7,iElem,iSpec) &
                     - Source(4:6,iElem,iSpec)/Source(7,iElem,iSpec)) / BoltzmannConst
    END IF
    Adsorption%PartitionTemp(iElem,iSpec) = (TempDirec(1) + TempDirec(2) + TempDirec(3)) / 3.
  END DO
END DO

END SUBROUTINE AnalyzePartitionTemp

END MODULE MOD_SurfaceModel_PartFunc


