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

MODULE MOD_SEE
!===================================================================================================================================
!> Main Routines of Surface Model: Secondary Electron Emission (SEE) due to particle bombardment
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
PUBLIC :: SecondaryElectronEmission
!===================================================================================================================================

CONTAINS

SUBROUTINE SecondaryElectronEmission(PartSurfaceModel_IN,PartID_IN,Adsorption_prob_OUT,adsorption_case,outSpec)
!----------------------------------------------------------------------------------------------------------------------------------!
! Determine the probability of an electron being emitted due to an impacting particles (ion/electron bombardment)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Vars    ,ONLY: PartState,Species,PartSpecies
USE MOD_Particle_Analyze ,ONLY: PartIsElectron
USE MOD_Globals_Vars     ,ONLY: BoltzmannConst,ElementaryCharge,ElectronMass
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)      :: PartSurfaceModel_IN !< which SEE model?
                                               !< 5: SEE by Levko2015
                                               !< 6: SEE by Pagonakis2016 (originally from Harrower1956)
INTEGER,INTENT(IN)      :: PartID_IN           !< Bombarding Particle ID
REAL   ,INTENT(OUT)     :: Adsorption_prob_OUT !< probability of an electron being emitted due to an impacting particles
                                               !< (ion/electron bombardment)
INTEGER,INTENT(OUT)     :: adsorption_case     !< what happens to the bombarding particle and is a new one created?
INTEGER,INTENT(OUT)     :: outSpec(2)          !< outSpec(1) currently not used
                                               !< outSpec(2) ID of created particle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: eps_e  ! Energy of bombarding electron in eV
REAL              :: v_new  ! Velocity of emitted secondary electron
REAL              :: iRan   ! Random number
REAL              :: k_ee   ! Coefficient of emission of secondary electron
REAL              :: k_refl ! Coefficient for reflection of bombarding electron
!===================================================================================================================================
! Default 0
outSpec = 0
! Select particle surface modeling
SELECT CASE(PartSurfaceModel_IN)
CASE(5) ! 5: SEE by Levko2015 for copper electrodes
  !     ! D. Levko and L. L. Raja, Breakdown of atmospheric pressure microgaps at high excitation, J. Appl. Phys. 117, 173303 (2015)
  ASSOCIATE (&
        phi            => 4.4  ,& ! eV -> cathode work function phi Ref. [20] Y. P. Raizer, Gas Discharge Physics (Springer, 1991)
        I              => 15.6  & ! eV -> ionization threshold of N2
        )
    IF(PARTISELECTRON(PartID_IN))THEN ! Bombarding electron
      ASSOCIATE (&
            delta_star_max => 1.06 ,& !    -> empir. fit. const. copper electrode Ref, [19] R. Cimino et al.,Phys.Rev.Lett. 2004
            s              => 1.35 ,& !    -> empir. fit. const. copper electrode Ref, [19] R. Cimino et al.,Phys.Rev.Lett. 2004
            eps_max        => 262  ,& ! eV -> empir. fit. const. copper electrode Ref, [19] R. Cimino et al.,Phys.Rev.Lett. 2004
            eps_0          => 150  ,& ! eV -> empir. fit. const. copper electrode Ref, [19] R. Cimino et al.,Phys.Rev.Lett. 2004
            velo2          => PartState(PartID_IN,4)**2 + PartState(PartID_IN,5)**2 + PartState(PartID_IN,6)**2 ,&
            mass           => Species(PartSpecies(PartID_IN))%MassIC &! mass of bombarding particle
            )
        ! Electron energy in [eV]
        eps_e = 0.5*mass*velo2/ElementaryCharge
        ! Calculate the electron impact coefficient
        IF(eps_e.LE.5)THEN ! Electron energy <= 5 eV
          k_ee = 0.
        ELSE
          k_ee = delta_star_max*s*( (eps_e/eps_max)/(s-1.+(eps_e/eps_max)**s) )
        END IF
        ! Calculate the elastic electron reflection coefficient
        k_refl = (SQRT(eps_e)-SQRT(eps_e+eps_0))**2/((SQRT(eps_e)+SQRT(eps_e+eps_0))**2)
        ! Decide SEE and/or reflection
        CALL RANDOM_NUMBER(iRan)
        IF(iRan.LT.(k_ee+k_refl))THEN ! Either SEE-E or reflection
          CALL RANDOM_NUMBER(iRan)
          IF(iRan.LT.k_ee/(k_ee+k_refl))THEN ! SEE
            adsorption_case = -2 ! SEE + perfect elastic scattering of the bombarding electron
            outSpec(2)      = 4  ! Species of the injected electron (TODO: make this variable)
            v_new           = SQRT(2.*(eps_e-ElementaryCharge*phi)/ElectronMass) ! Velocity of emitted secondary electron
            eps_e           = 0.5*mass*(v_new**2)/ElementaryCharge               ! Energy of the injected electron
          ELSE
            adsorption_case = -1 ! Only perfect elastic scattering of the bombarding electron
          END IF
          !   IF(k_refl.LT.k_ee)THEN ! -> region with a high chance of SEE
          !     IF(iRan.LT.k_refl/(k_ee+k_refl))THEN ! Reflection
          !       adsorption_case = -1 ! only perfect elastic scattering of the bombarding electron
          !     ELSE ! SEE-E
          !       adsorption_case = -2 ! SEE + perfect elastic scattering of the bombarding electron
          !       outSpec(2)      = 4
          !     END IF
          !   ELSE ! (k_ee.LE.k_refl) -> region with a high chance of reflection
          !     IF(iRan.LT.k_ee/(k_ee+k_refl))THEN ! SEE
          !       adsorption_case = -2 ! SEE + perfect elastic scattering of the bombarding electron
          !       outSpec(2)      = 4
          !     ELSE ! Reflection
          !       adsorption_case = -1 ! only perfect elastic scattering of the bombarding electron
          !     END IF
          !   END IF
        ELSE
          adsorption_case = -4 ! Removal of the bombarding electron
        END IF
      END ASSOCIATE
    ELSEIF(Species(PartSpecies(PartID_IN))%ChargeIC.NE.0.0)THEN ! Positive bombarding ion
      CALL RANDOM_NUMBER(iRan)
      IF(iRan.LT.0.02)THEN ! SEE-I: gamma=0.02 for the N2^+ ions and copper material
        adsorption_case = -2       ! SEE + perfect elastic scattering of the bombarding electron
        outSpec(2)      = 4        ! Species of the injected electron (TODO: make this variable)
        eps_e           = I-2.*phi ! Energy of the injected electron
      ELSE
        adsorption_case = -1 ! Only perfect elastic scattering of the bombarding electron
      END IF
    ELSE ! Neutral bombarding particle
      IF(iRan.LT.0.1)THEN ! SEE-N: from svn-trunk PICLas version
        adsorption_case = -2 ! SEE + perfect elastic scattering of the bombarding electron
        outSpec(2)      = 4  ! Species of the injected electron (TODO: make this variable)
      ELSE
        adsorption_case = -1 ! Only perfect elastic scattering of the bombarding electron
      END IF
    END IF
  END ASSOCIATE
CASE(6) ! 6: SEE by Pagonakis2016 (originally from Harrower1956)
END SELECT

END SUBROUTINE SecondaryElectronEmission


END MODULE MOD_SEE
