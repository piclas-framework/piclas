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
MODULE MOD_DSMC_SpecXSec
!===================================================================================================================================
! Contains the Argon Ionization
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE

INTERFACE XSec_Argon_DravinLotz
  MODULE PROCEDURE XSec_Argon_DravinLotz
END INTERFACE

PUBLIC :: XSec_Argon_DravinLotz
!===================================================================================================================================

CONTAINS

SUBROUTINE XSec_Argon_DravinLotz(SpecToExec, iPair)
!===================================================================================================================================
! Subroutine computing the collision probability o the Argion ionization
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : Coll_pData, SpecDSMC
  USE MOD_Equation_Vars,          ONLY : eps0
  USE MOD_Globals_Vars,           ONLY : Pi, ElementaryCharge
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: SpecToExec, iPair
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: BohrRad, Rydberg
!===================================================================================================================================

! local constants
BohrRad = 0.5291772109E-10
Rydberg = 13.60569253*ElementaryCharge

!.... Elastic scattering cross section
  Coll_pData(iPair)%Sigma(1) = SQRT(0.5*Pi*SpecDSMC(SpecToExec)%RelPolarizability &
                             * BohrRad**3*ElementaryCharge**2     &   ! AIAA07 Paper
                             / (eps0*Coll_pData(iPair)%Ec))                    ! units checked

!.... Ionization cross section (Lotz)
IF ((Coll_pData(iPair)%Ec/ElementaryCharge).GE.SpecDSMC(SpecToExec)%Eion_eV) THEN
  Coll_pData(iPair)%Sigma(2) = 2.78*SpecDSMC(SpecToExec)%NumEquivElecOutShell*Pi &
             * BohrRad**2*Rydberg**2 &    ! units checked
             / (Coll_pData(iPair)%Ec*SpecDSMC(SpecToExec)%Eion_eV*ElementaryCharge) &
             * LOG(Coll_pData(iPair)%Ec/(ElementaryCharge*SpecDSMC(SpecToExec)%Eion_eV))
ELSE
  Coll_pData(iPair)%Sigma(2) = 0.0
ENDIF
Coll_pData(iPair)%Sigma(0)=Coll_pData(iPair)%Sigma(1)+Coll_pData(iPair)%Sigma(2) ! Calc of Sigma total

END SUBROUTINE XSec_Argon_DravinLotz

END MODULE MOD_DSMC_SpecXSec
