!==================================================================================================================================
! Copyright (c) 2010 - 2019 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_Radiation_Molecules
!===================================================================================================================================
! Module for Radiation
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE radiation_molecules
  MODULE PROCEDURE radiation_molecules
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: radiation_molecules
!===================================================================================================================================

CONTAINS


SUBROUTINE radiation_molecules(iElem, em_mol)
!===================================================================================================================================
! Main routine of molecular radiation calculation
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_Globals_Vars,      ONLY   : BoltzmannConst, PlanckConst, AtomicMassUnit, Pi, c
  USE MOD_Radiation_Vars,    ONLY   : RadiationInput, RadiationParameter, SpeciesRadiation,  &
                                      Radiation_Emission_spec, Radiation_Absorption_spec, NumDensElectrons, TElectrons, &
                                      Radiation_ElemEnergy_Species, Radiation_Absorption_SpeciesWave
  USE MOD_Particle_Vars,     ONLY   : nSpecies, Species
  USE MOD_DSMC_Vars,         ONLY   : SpecDSMC

! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT(INOUT)             :: em_mol
  INTEGER, INTENT(IN)             :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

  INTEGER           :: iSpec, iWave, iBand, iTrans, i, iWaveCoarse!, hilf
  REAL              :: rho, ptot                            !in ini !density [kg/m**3], pressure [Pa], number density    
  REAL, ALLOCATABLE :: epsilon_mol(:), epsilon_iSpec(:), abs_iSpec(:)
  INTEGER           :: istart, iend
  REAL              :: v_upper, v_lower                     ! vibrational quantum number of upper/lower state
  REAL              :: ElecTransMoment_squared              ! squared electronic transition moment Re(r_{v_upper,v_lower})^2
  REAL              :: Bv_upper, Bv_lower                   ! rotational constants for the vibrational level v_upper/lower, cm^-1
  REAL              :: Dv_upper, Dv_lower                   ! rotational constants for the vibrational level v_upper/lower, cm^-1
  REAL              :: G_upper, G_lower                     ! G(v) = we(v+0.5) - wexe(v+0.5)**2 + weye(v+0.5)**3 + weze(v+0.5)**4
  REAL              :: cint1, cint2, nubar0, r1             ! TODO: RENAME, nubar0 [cm**-1]
  INTEGER           :: kmax                                 ! TODO: RENAME maximum rotational quantum number (NEQAIR85 determines the maximum rotational quantum number by determining approximately 0.01 of the maximum radiation. This approach is not longer followed since the combination of the radiating band to more transitions had the effect that no rotational lines were computed for higher vibrational levels. Now, the maximum rotational quantum number is determined from the dissociation energy (only bvu is used, since otherwise an iterative procedure would be necessary) mf 2007/07/18)
  REAL              :: DlamG, DlamL, DlamV                  ! Gaussian, Lorentzian and Voigt width_int [m]
  REAL              :: DlamL1, DlamL2, DlamL3, DlamL4       ! 1: collisions with molecules, 2: collisions with all other, 3: natural, 4: Stark [m]
  REAL, PARAMETER   :: Stark = 0.3                          ! Stark width at NumDensElectrons=1.E16 and TElectrons = 10000
  REAL              :: NumDensForeignGas,MolWeightForeignGas! number density and average molecular weight of foreign gas

  INTEGER, PARAMETER:: Voigt_int_nmax = 1001
  REAL, PARAMETER   :: Voigt_mult_discr = 25.
  REAL              :: Voigt_int(-Voigt_int_nmax:Voigt_int_nmax), Voigt_int_step, Voigt
  REAL              :: Voigt_var1, Voigt_const1, Voigt_const2, Voigt_int_distmax, Voigt_arg_dimless, Voigt_int_NormEnergy
  INTEGER           :: iVoigt
  REAL              :: Voigt_dist1, Voigt_dist2, Voigt_dist225

  REAL              :: j_upper, j_lower, j_transition
  REAL              :: k_upper, k_lower
  INTEGER           :: k, iBranch
  REAL              :: HoenlLondon_3p3p(0:500,8), HoenlLondon_2p2p(500,12), HoenlLondon_2s2s(0:500,6) !hoenl-london compution

  INTEGER, PARAMETER:: Deltak_upper_3p3p(8)  = [1,-1,1,-0,-1,1,0,-1]
  INTEGER, PARAMETER:: Deltak_lower_3p3p(8)  = [0, 0,0,0, 0,0,0, 0]
  INTEGER, PARAMETER:: Deltak_upper_2p2p(12) = [-1,0,1, 0,0,1,-1,0,1,-1,0,1]
  INTEGER, PARAMETER:: Deltak_lower_2p2p(12) = [0, 0,0,-1,0,0, 0,0,0, 0,0,0]
  INTEGER, PARAMETER:: Deltak_upper_2s2s(6)  = [1,1,-1,-1,1,-1]
  INTEGER, PARAMETER:: Deltak_lower_2s2s(6)  = [0,0, 0, 0,0, 0]

  REAL              :: Lambda_upper, Lambda_lower           ! Bahndrehimpulsquantenzahl
  REAL              :: Y_upper, Y_lower,Y_transition        ! Y = A/Bv
  REAL              :: F_upper, F_lower
  REAL              :: Z1_upper, Z1_lower, Z2_upper, Z2_lower
  REAL              :: nubar, r2!, lambda                    ! wavenumber [cm**-1], , wavelength [m]

  REAL              :: etot, abstot, eps
  REAL, PARAMETER   :: a0e = 2.5415785E-18
  REAL              :: AlternationFactor                    ! alternation factor for homocuclear molecules
  REAL              :: HoenlLondon, u

  REAL, ALLOCATABLE :: Radiation_Profile(:)
  REAL              :: TempOut_Em, TempOut_Abs
!===================================================================================================================================
  
! --- loop for ptot and rho over all atoms and molecules, additionally electrons
  rho  = 0.0
  ptot = 0.0

  DO iSpec = 1, nSpecies
    rho  = rho + RadiationInput(iSpec)%NumDens * Species(iSpec)%MassIC
    ptot = ptot + RadiationInput(iSpec)%NumDens * BoltzmannConst * RadiationInput(iSpec)%Ttrans(4)
  END DO
  ptot = ptot + NumDensElectrons * BoltzmannConst * TElectrons ! TODO: isp -> should be electrons? Radiation: NumDens(N2), Telec

!  --- initialize emission coefficient
  ALLOCATE(epsilon_mol(RadiationParameter%WaveLenDiscr))
  ALLOCATE(epsilon_iSpec(RadiationParameter%WaveLenDiscr))
  ALLOCATE(abs_iSpec(RadiationParameter%WaveLenDiscr))
  ALLOCATE(Radiation_Profile(RadiationParameter%WaveLenDiscr))
  DO iWave=1, RadiationParameter%WaveLenDiscr
    epsilon_mol(iWave) = 0.0 !Radiation_Emission_spec(iWave,iElem)
  END DO

  istart = RadiationParameter%WaveLenDiscr
  iend   = 1

  DO iSpec = 1, nSpecies
    IF(.NOT.RadiationInput(iSpec)%DoRadiation) CYCLE
    IF((SpecDSMC(iSpec)%InterID .NE. 2) .AND. (SpecDSMC(iSpec)%InterID .NE. 20)) CYCLE
    Radiation_Profile   = 0.0
    IF ((RadiationInput(iSpec)%Telec.LT.10.0).OR.(RadiationInput(iSpec)%Tvib.LT.10.0).OR.(RadiationInput(iSpec)%NumDens.LT.10.0).OR.(RadiationInput(iSpec)%Ttrans(4).LT.10.0))CYCLE 

    epsilon_iSpec = 0.0
    abs_iSpec = 0.0

    IF(SpeciesRadiation(iSpec)%nBands.NE.0) THEN

! --- determine number density and average weight of foreign gas
      NumDensForeignGas = 1.0E10
      IF (ptot .GT. 1.0E-16) &
          NumDensForeignGas = MAX(ptot/(RadiationInput(iSpec)%Ttrans(4)*BoltzmannConst) - RadiationInput(iSpec)%NumDens , &
      1.0D-40*ptot/(RadiationInput(iSpec)%Ttrans(4)*BoltzmannConst))
    
      MolWeightForeignGas = 20. * AtomicMassUnit
      IF(rho .GT. 1.0E-16) MolWeightForeignGas = (rho - Species(iSpec)%MassIC*RadiationInput(iSpec)%NumDens) / NumDensForeignGas
      MolWeightForeignGas = MolWeightForeignGas / (AtomicMassUnit * 1.0E-3)
      IF(MolWeightForeignGas .LT. 10.*1.0E-3*AtomicMassUnit)   MolWeightForeignGas = 10.  * 1.0E-3 * AtomicMassUnit
      IF(MolWeightForeignGas .GT. 28.85*1.0E-3*AtomicMassUnit) MolWeightForeignGas = 28.85* 1.0E-3 * AtomicMassUnit

      DO iBand = 1,SpeciesRadiation(iSpec)%nBands
      ! WRITE(*,*) TRIM(SpeciesRadiation(iSpec)%BandName(iBand))

  !!      TODO : IF ( 2 .EQ. 2) diatomic molecule
  ! -----------------------------------------------------------------------------------------------------------------------------------
  ! ------------------------------------------------- diatomic molecule calculation ---------------------------------------------------
  ! -----------------------------------------------------------------------------------------------------------------------------------

  ! --- cycle over all transitions of electronic bands
        DO iTrans = 1, SpeciesRadiation(iSpec)%NumMolecularTransitions(iBand)
          v_upper = SpeciesRadiation(iSpec)%Bands(iBand)%MolTransLines(iTrans,1)
          v_lower = SpeciesRadiation(iSpec)%Bands(iBand)%MolTransLines(iTrans,2)
          ElecTransMoment_squared = SpeciesRadiation(iSpec)%Bands(iBand)%MolTransLines(iTrans,3) &
            * REAL(SpeciesRadiation(iSpec)%BandProperties(iBand,4))
  ! --- rotational constants for this band
          Bv_upper = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1), 8) &
            - SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1), 9) * (v_upper+0.5)
          Bv_lower = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2), 8) &
            - SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2), 9) * (v_lower+0.5)
          Dv_upper = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),13) &
            + SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),14) * (v_upper+0.5)
          Dv_lower = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),13) &
            + SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),14) * (v_lower+0.5)

  ! --- constants for intensity equation in rotational structure
          G_upper =SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1) ,4) * (v_upper+0.5) &
              - SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1) ,5) * (v_upper+0.5)**2 &
              + SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1) ,6) * (v_upper+0.5)**3 &
              + SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1) ,7) * (v_upper+0.5)**4
          G_lower =SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2) ,4) * (v_lower+0.5) &
              - SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2) ,5) * (v_lower+0.5)**2 &
              + SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2) ,6) * (v_lower+0.5)**3 &
              + SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2) ,7) * (v_lower+0.5)**4

          cint1 = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),2) &
              / RadiationInput(iSpec)%Telec + G_upper / RadiationInput(iSpec)%Tvib
          cint2 = (16.0E-5*c* RadiationInput(iSpec)%NumDens * 1.0E-6 * SpeciesRadiation(iSpec)%Bands(iBand)%MolTransLines(iTrans,4) &
              * ElecTransMoment_squared * Pi**3) / (3.0 * SpeciesRadiation(iSpec)%PartFunc)

          kmax  = INT( SQRT( MAX(0.D0, &
              0.25 + ( SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),3) &
              - G_lower ) / Bv_lower ) + 1. ) )

  ! --- band origin
          nubar0 = (( SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),2) &
              - SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),2) ) &
              + (G_upper - G_lower)) / (100.*PlanckConst*c) ! energy conversion J -> cm**-1

          r1 = EXP( ( SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),2) &
              - SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),2) ) &
              / (BoltzmannConst * RadiationInput(iSpec)%Telec) + (G_upper-G_lower) / (BoltzmannConst * RadiationInput(iSpec)%Tvib) )

  ! --- check if transition is in considered wavelength interval
          IF( (2.0/(100.*nubar0)) .LT. RadiationParameter%WaveLen(1) ) CYCLE
          IF( (0.5/(100.*nubar0)) .GT. RadiationParameter%WaveLen(RadiationParameter%WaveLenDiscr) ) CYCLE

  ! --- widths        
          DlamG = 7.16E-7 / (100.*nubar0) * SQRT(RadiationInput(iSpec)%Ttrans(4)/(Species(iSpec)%MassIC/AtomicMassUnit))

          DlamL1 = 1.33E-29 * SQRT(2.0/(Species(iSpec)%MassIC/AtomicMassUnit)) &
              * (0.01/nubar0)**2 * SQRT(RadiationInput(iSpec)%Ttrans(4)) * RadiationInput(iSpec)%NumDens * 1.0E-6 * 10. * 1.0E10
          DlamL2 = 5.85E-30 * SQRT(1.0/(Species(iSpec)%MassIC/AtomicMassUnit)+1./(MolWeightForeignGas/(AtomicMassUnit*1.0E-3))) &
              * (0.01/nubar0)**2 * SQRT(RadiationInput(iSpec)%Ttrans(4)) * NumDensForeignGas * 1.0E-6 * 10. * 1.0E10
          DlamL3 = 1.18E-4 * 1.0E-10
          DlamL4 = 1.0E-8 * Stark * (1.0E-22 * MAX(NumDensElectrons, 1.E-6))**0.6 * (0.01/nubar0)**2 * 1.0E10
          DlamL  = DlamL1 + DlamL2 + DlamL3 + DlamL4

  ! --- Voigt-line width at half maximum
          DlamV  = Dlaml/2.0 + SQRT( Dlaml**2/4.0 + Dlamg**2 )

  ! --- generate Voigt line profile
          Voigt_int(-Voigt_int_nmax) = 0.                             !initialize Voigt profile
          Voigt_int( Voigt_int_nmax) = 1.

          Voigt_int_distmax = Voigt_mult_discr * DlamV                !determine discretised interval
          Voigt_arg_dimless = DlamL/DlamV
          Voigt_var1        = 1./((1.065 + (0.447 + 0.058 * Voigt_arg_dimless) * Voigt_arg_dimless) * DlamV)
          Voigt_const2      = Voigt_arg_dimless * Voigt_var1
          Voigt_const1      = Voigt_var1 - Voigt_const2
          Voigt_int_step    = Voigt_int_distmax / (Voigt_int_nmax-1) !step width

          DO iVoigt = -(Voigt_int_nmax-1), 0                          !determine lower half of Voigt profile
            Voigt_dist1   = REAL(ABS(iVoigt)) * Voigt_int_step / (DlamV)
            Voigt_dist2   = Voigt_dist1**2
            Voigt_dist225 = Voigt_dist2 * SQRT(SQRT(Voigt_dist1))

            Voigt = Voigt_const1 * EXP(MAX(-1.E8,-2.772*Voigt_dist2)) + Voigt_const2 / (1.+4.*Voigt_dist2) &
                  + 0.016 * Voigt_const2 * (1.-Voigt_arg_dimless) * (EXP(MAX(-1.E5,-0.4*Voigt_dist225) &
                  - 10. / (10. + Voigt_dist225)))
            IF((Voigt_dist225.GT.1E5).OR.(Voigt_dist2.GT.1E8)) CALL abort(&
                __STAMP__&
                ,' ERROR: Voigt_dist225 is too big!')
            Voigt_int(iVoigt) = Voigt_int(iVoigt-1) + Voigt * Voigt_int_step
          END DO

          DO iVoigt = 1, (Voigt_int_nmax-1)                            ! determine upper half of Voigt profile
            Voigt_int(iVoigt) = 2. * Voigt_int(0) - Voigt_int(-iVoigt)
          END DO

          Voigt_int_NormEnergy = 1./ Voigt_int(Voigt_int_nmax-1)

          DO iVoigt = -(Voigt_int_nmax-1), (Voigt_int_nmax-1)         ! normalize profiles that energy below the Voigt profile function is 1 (integrated profiles are normed to 1)
            Voigt_int(iVoigt) = Voigt_int(iVoigt) * Voigt_int_NormEnergy
          END DO


  ! --- different radiating systems
          SELECT CASE (SpeciesRadiation(iSpec)%BandProperties(iBand,3))

            CASE(1) ! parallel transition, Delta Lambda = 0

              etot = 0.0

            ! --- set constants to determine wavelength of center line for triplets (Herzberg p.235)
              Lambda_upper = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),17)
              Lambda_lower = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),17)
              Y_upper = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),16) / Bv_upper
              Y_lower = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),16) / Bv_lower

              DO iBranch = 1,2

                SELECT CASE(iBranch)
                  CASE(1) !set quantum numbers for P branch (+1)
                    k_upper = 0.0
                    k_lower = k_upper + 1.0
                  CASE(2) !set quantum numbers for R branch (-1)
                    k_upper = 0.0
                    k_lower = k_upper - 1.0
                  CASE DEFAULT
                    WRITE(*,*) 'branch not found'
                END SELECT

              ! --- compute and distribute the integrated integrated intensity due to spontaneous emission of all specified rotational lines for the appropriate branch
                DO k = 0, kmax

                  Z1_upper = Lambda_upper**2 * Y_upper * (Y_upper - 4.) + 4./3. + 4. * k_upper * (k_upper + 1.)
                  Z1_lower = Lambda_lower**2 * Y_lower * (Y_lower - 4.) + 4./3. + 4. * k_lower * (k_lower + 1.)

                  IF(Lambda_upper .EQ. 0.0) THEN
                    IF(SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),1) .NE. 3.0) THEN
                      Z2_upper = 0.0
                      Z2_lower = 0.0
                    ELSE
                      Z2_upper = (Lambda_upper**2 * Y_upper * (Y_upper-1.) - 4./9. - 2. * k_upper * (k_upper+1.)) / (3.*Z1_upper) ! TODO: + or - 4/9?
                      Z2_lower = (Lambda_lower**2 * Y_lower * (Y_lower-1.) - 4./9. - 2. * k_lower * (k_lower+1.)) / (3.*Z1_lower) ! TODO: + or - 4/9?
                    END IF
                  ELSE
                    IF((SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),1)/2.0) &
                      .NE. 3.0) THEN
                      Z2_upper = 0.0
                      Z2_lower = 0.0
                    ELSE
                      Z2_upper = (Lambda_upper**2 * Y_upper * (Y_upper-1.) - 4./9. - 2. * k_upper * (k_upper+1.)) / (3.*Z1_upper) ! TODO: + or - 4/9?
                      Z2_lower = (Lambda_lower**2 * Y_lower * (Y_lower-1.) - 4./9. - 2. * k_lower * (k_lower+1.)) / (3.*Z1_lower) ! TODO: + or - 4/9?
                    END IF
                  END IF

                ! --- compute wavelength of line center
                  F_upper = Bv_upper * (k_upper * (k_upper + 1.) + 4. * Z2_upper) - Dv_upper * (k_upper + 0.5)**4
                  F_lower = Bv_lower * (k_lower * (k_lower + 1.) + 4. * Z2_lower) - Dv_lower * (k_lower + 0.5)**4
                  nubar = nubar0 + (F_upper - F_lower) / (100.*PlanckConst*c)
                  r2    = EXP( (F_upper-F_lower) / (BoltzmannConst*RadiationInput(iSpec)%Trot) ) * r1

                  SELECT CASE(iBranch)
                    CASE(1) !set HLF for P branch (+1)
                      HoenlLondon = k_upper + 1.0
                    CASE(2) !set HLF for R branch (-1)
                      HoenlLondon = k_upper
                  END SELECT

                ! --- evaluate the alternation factor for homocuclear molecules
                  IF(SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),12) .EQ. 0.0) THEN
                    AlternationFactor = 1.0
                  ELSE
                    AlternationFactor = 1.0 + (-1.0)**NINT(k_lower &
                      + SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),12)) &
                      / (2.0*SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),11) + 1.0) !Radiation: kexp=INT(jl+altnat(uplev)+0.1)
                  END IF

                ! --- find integrated line intensity due to spontaneous emission (etot) and sum of absorption and stimulated emission (abstot)
                  eps = 1.0E16 * AlternationFactor * HoenlLondon * (nubar**2 * a0e)**2 * cint2 &
                      * EXP(-(cint1 + (Bv_upper * k_upper * (k_upper + 1.)) / RadiationInput(iSpec)%Trot) / BoltzmannConst)
                  abstot = (nubar**(-5)) * (r2 - 1.0) * eps / (2.E10 * PlanckConst * c**2)
!                  abstot = c*c/(8.*Pi*nubar*nubar) &
!                      * 16.*Pi**3*nubar**3/(3.*8.8541878128E-12*PlanckConst*c**3) &
!                      * SpeciesRadiation(iSpec)%Bands(iBand)%MolTransLines(iTrans,4)*ElecTransMoment_squared*HoenlLondon/(2.*J+1.) &
!                      * (RadiationInput(iSpec)%NumDens / SpeciesRadiation(iSpec)%PartFunc &
!                      * SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),2) &
!                      * EXP(-(SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),2) &
!                      / RadiationInput(iSpec)%Telec + G_upper / RadiationInput(iSpec)%Tvib
!                      + (Bv_upper * k_upper * (k_upper + 1.)) / RadiationInput(iSpec)%Trot) / BoltzmannConst) &
!                      * SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),1) &
!                      / SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),2) &
!                      - (RadiationInput(iSpec)%NumDens / SpeciesRadiation(iSpec)%PartFunc &
!                      * SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),1) &
!                      * EXP(-(SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),2) &
!                      / RadiationInput(iSpec)%Telec + G_lower / RadiationInput(iSpec)%Tvib
!                      + (Bv_lower * k_lower * (k_lower + 1.)) / RadiationInput(iSpec)%Trot) / BoltzmannConst) )

!                  abstot = 2.*Pi*Pi*nubar/(3.*8.8541878128E-12*PlanckConst*c) &
!                      * SpeciesRadiation(iSpec)%Bands(iBand)%MolTransLines(iTrans,4)*ElecTransMoment_squared*HoenlLondon/(2.*J+1.) &
!                      * (RadiationInput(iSpec)%NumDens / SpeciesRadiation(iSpec)%PartFunc &
!                      * SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),2) &
!                      * EXP(-(SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),2) &
!                      / RadiationInput(iSpec)%Telec + G_upper / RadiationInput(iSpec)%Tvib
!                      + (Bv_upper * k_upper * (k_upper + 1.)) / RadiationInput(iSpec)%Trot) / BoltzmannConst) &
!                      * SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),1) &
!                      / SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),2) &
!                      - (RadiationInput(iSpec)%NumDens / SpeciesRadiation(iSpec)%PartFunc &
!                      * SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),1) &
!                      * EXP(-(SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),2) &
!                      / RadiationInput(iSpec)%Telec + G_lower / RadiationInput(iSpec)%Tvib
!                      + (Bv_lower * k_lower * (k_lower + 1.)) / RadiationInput(iSpec)%Trot) / BoltzmannConst) )
                  etot = etot + eps

                ! --- calculate line profile and add radiative energy to emission
                  IF(eps .GE. 1.0E-25) THEN
                    CALL Radiation_Molecular_Transition_Line_Profile(Radiation_Profile, nubar, Voigt_int, Voigt_int_distmax, &
                      epsilon_mol, epsilon_iSpec, eps, Radiation_Absorption_spec, abs_iSpec, abstot, iElem)
                  END IF

                  k_lower = k_lower + 1.0
                  k_upper = k_upper + 1.0

                END DO !k
              END DO !iBranch

            CASE(2) ! perpendicular transition, Delta Lambda = +-1

            ! --- adjust factor used in intesity equation to correct partcc for lambda doubling
              IF(SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),17) .GE. 0.5) &
                cint2 = cint2 / 2.0

              etot = 0.0

              DO iBranch = 1,3

                SELECT CASE(iBranch)
                  CASE(1) !set quantum numbers for P branch (+1)
                    k_upper = 1.0
                    k_lower = k_upper + 1.0
                  CASE(2) !set quantum numbers for Q branch ( 0)
                    k_upper = 1.0
                    k_lower = k_upper
                  CASE(3) !set quantum numbers for R branch (-1)
                    k_upper = 1.0
                    k_lower = k_upper - 1.0
                  CASE DEFAULT
                    WRITE(*,*) 'branch not found'
                END SELECT

              ! --- compute and distribute the integrated intensity due to spontaneous emission of all specified rotational lines
                DO k = 1, kmax
                ! --- compute line center
                  F_upper = Bv_upper * (k_upper * (k_upper + 1.0)) - Dv_upper * (k_upper + 0.5)**4
                  F_lower = Bv_lower * (k_lower * (k_lower + 1.0)) - Dv_lower * (k_lower + 0.5)**4
                  nubar = nubar0 + (F_upper - F_lower) / (100.*PlanckConst*c)
    !              lambda = 1./(100.*nubar)
                  r2    = EXP( (F_upper-F_lower) / (BoltzmannConst*RadiationInput(iSpec)%Trot) ) * r1

                ! --- strength factor, set sign of Lambda used in strength equations, appropriate to the sigh of Delta Lambda (Johnson p.150)
                  IF(SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),17) .LT. &
                  SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),17)) THEN
                    SELECT CASE(iBranch)
                      CASE(1)
                        HoenlLondon = (k_upper + 1.0 &
                          + SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),17)) &
                          * (k_upper + 2.0 &
                          + SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),17)) &
                          / (k_upper + 1.0)
                      CASE(2)
                        HoenlLondon = (k_upper &
                          - SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),17)) &
                          * (2.0 * k_upper + 1.0) &
                          * (k_upper + 1.0 &
                          + SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),17)) &
                          / (k_upper * (k_upper + 1.0))
                      CASE(3)
                        HoenlLondon = (k_upper &
                          - SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),17)) &
                          * (k_upper - 1.0 &
                          - SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),17)) &
                          / k_upper
                    END SELECT
                  ELSE
                    SELECT CASE(iBranch)
                      CASE(1)
                        HoenlLondon = (k_upper + 1.0 &
                          - SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),17)) &
                          * (k_upper + 2.0 &
                          - SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),17)) &
                          / (k_upper + 1.0)
                      CASE(2)
                        HoenlLondon = (k_upper &
                          + SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),17)) &
                          * (2.0 * k_upper + 1.0) &
                          * (k_upper + 1.0 &
                          - SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),17)) &
                          / (k_upper * (k_upper + 1.0))
                      CASE(3)
                        HoenlLondon = (k_upper &
                          + SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),17)) &
                          * (k_upper - 1.0 &
                          + SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),17)) &
                          / k_upper
                    END SELECT
                  END IF

                ! --- determine if lines alternate in intesity -> evaluate the alternation factor for homocuclear molecules
                  IF(SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),12) .EQ. 0.0) THEN
                    AlternationFactor = 1.0
                  ELSE
                    AlternationFactor = 1.0 + (-1.0)**NINT(k_lower &
                      + SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),12)) &
                      / (2.0*SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),11) + 1.0) !Radiation: kexp=INT(jl+altnat(uplev)+0.1)
                  END IF

                  eps = 1.0E16 * AlternationFactor * HoenlLondon * (nubar**2 * a0e)**2 * cint2 &
                    * EXP(-(cint1 + (Bv_upper * k_upper * (k_upper + 1.)) / RadiationInput(iSpec)%Trot) / BoltzmannConst)
                  abstot = (nubar**(-5)) * (r2 - 1.0) * eps / (2.E10 * PlanckConst * c**2)
                  etot = etot + eps

                ! --- calculate line profile and add radiative energy to emission
                  IF(eps .GE. 1.0E-25) THEN
                    CALL Radiation_Molecular_Transition_Line_Profile(Radiation_Profile, nubar, Voigt_int, Voigt_int_distmax, &
                      epsilon_mol, epsilon_iSpec, eps, Radiation_Absorption_spec, abs_iSpec, abstot, iElem)
                  END IF

                  k_lower = k_lower + 1.0
                  k_upper = k_upper + 1.0

                END DO !k
              END DO ! iBranch

            CASE(3) ! 2Sigma -> 2Pi

              Lambda_upper = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),17)
              Lambda_lower = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),17)
              Y_upper = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),16) / Bv_upper
              Y_lower = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),16) / Bv_lower

            ! --- set Y_transition depending on direction of transition
              IF(Lambda_upper .GT. Lambda_lower) THEN !Pi -> Sigma transition
                Y_transition = Y_upper
              ELSE                                    !Sigma -> Pi transition
                Y_transition = Y_lower
              END IF

              etot = 0.0

              DO iBranch = 1,8

              ! --- Pi -> Sigma transition
                IF(Lambda_upper .GT. Lambda_lower) THEN
                  SELECT CASE(iBranch)
                    CASE(1) !set quantum numbers for P2 branch
                      j_upper = 1.5
                      j_lower = j_upper + 1.0
                    CASE(2) !set quantum numbers for R1 branch
                      j_upper = 2.5
                      j_lower = j_upper - 1.0
                    CASE(3) !set quantum numbers for R21 branch
                      j_upper = 1.5
                      j_lower = j_upper - 1.0
                    CASE(4) !set quantum numbers for OP12 branch
                      j_upper = 2.5
                      j_lower = j_upper + 1.0
                    CASE(5) !set quantum numbers for Q2/QP21 branch
                      j_upper = 1.5
                      j_lower = j_upper
                    CASE(6) !set quantum numbers for Q1/QR12 branch
                      j_upper = 2.5
                      j_lower = j_upper
                    CASE(7) !set quantum numbers for R2/RQ21 branch
                      j_upper = 1.5
                      j_lower = j_upper - 1.0
                    CASE(8) !set quantum numbers for P1/PQ12 branch
                      j_upper = 2.5
                      j_lower = j_upper + 1.0
                    CASE DEFAULT
                      WRITE(*,*) 'branch not found'
                  END SELECT
              ! --- Sigma -> Pi transition
                ELSE
                  SELECT CASE(iBranch)
                    CASE(1) !set quantum numbers for R2 branch
                      j_upper = 1.5
                      j_lower = j_upper - 1.0
                    CASE(2) !set quantum numbers for P1 branch
                      j_upper = 2.5
                      j_lower = j_upper + 1.0
                    CASE(3) !set quantum numbers for SR21 branch
                      j_upper = 1.5
                      j_lower = j_upper - 1.0
                    CASE(4) !set quantum numbers for OP12 branch
                      j_upper = 2.5
                      j_lower = j_upper + 1.0
                    CASE(5) !set quantum numbers for Q2/QR12 branch
                      j_upper = 1.5
                      j_lower = j_upper
                    CASE(6) !set quantum numbers for Q1/QP21 branch
                      j_upper = 2.5
                      j_lower = j_upper
                    CASE(7) !set quantum numbers for P2/PQ12 branch
                      j_upper = 1.5
                      j_lower = j_upper + 1.0
                    CASE(8) !set quantum numbers for R1/RQ21 branch
                      j_upper = 2.5
                      j_lower = j_upper - 1.0
                    CASE DEFAULT
                      WRITE(*,*) 'branch not found'
                  END SELECT
                END IF

              ! --- set j_transition depending on direction of transition, j is the rotational quatum number of the Pi state (reference by Earls)
                IF(Lambda_upper .GT. Lambda_lower) THEN !Pi -> Sigma transition
                  j_transition = j_upper
                ELSE                                    !Sigma -> Pi transition
                  j_transition = j_lower
                END IF

                DO k = 2, kmax

                  IF((iBranch .EQ. 1) .OR. (iBranch .EQ. 5) .OR. (iBranch .EQ. 7)) THEN
                    F_upper = Bv_upper * ( (j_upper+0.5)**2 - Lambda_upper**2 + .5 * SQRT(4. * (j_upper+0.5)**2 &
                      - 4. * Y_upper * Lambda_upper**2 + (Y_upper * Lambda_upper)**2 ))
                    F_lower = Bv_lower * ( (j_lower+0.5)**2 - Lambda_lower**2 + .5 * SQRT(4. * (j_lower+0.5)**2 &
                      - 4. * Y_lower * Lambda_lower**2 + (Y_lower * Lambda_lower)**2 ))
                  ELSEIF ((iBranch .EQ. 2) .OR. (iBranch .EQ. 6) .OR. (iBranch .EQ. 8)) THEN
                    F_upper = Bv_upper * ( (j_upper+0.5)**2 - Lambda_upper**2 - .5 * SQRT(4. * (j_upper+0.5)**2 &
                      - 4. * Y_upper * Lambda_upper**2 + (Y_upper * Lambda_upper)**2 ))
                    F_lower = Bv_lower * ( (j_lower+0.5)**2 - Lambda_lower**2 - .5 * SQRT(4. * (j_lower+0.5)**2 &
                      - 4. * Y_lower * Lambda_lower**2 + (Y_lower * Lambda_lower)**2 ))
                  ELSEIF (iBranch .EQ. 3) THEN
                    F_upper = Bv_upper * ( (j_upper+0.5)**2 - Lambda_upper**2 + .5 * SQRT(4. * (j_upper+0.5)**2 &
                      - 4. * Y_upper * Lambda_upper**2 + (Y_upper * Lambda_upper)**2 ))
                    F_lower = Bv_lower * ( (j_lower+0.5)**2 - Lambda_lower**2 - .5 * SQRT(4. * (j_lower+0.5)**2 &
                      - 4. * Y_lower * Lambda_lower**2 + (Y_lower * Lambda_lower)**2 ))
                  ELSEIF (iBranch .EQ. 4) THEN
                    F_upper = Bv_upper * ( (j_upper+0.5)**2 - Lambda_upper**2 - .5 * SQRT(4. * (j_upper+0.5)**2 &
                      - 4. * Y_upper * Lambda_upper**2 + (Y_upper * Lambda_upper)**2 ))
                    F_lower = Bv_lower * ( (j_lower+0.5)**2 - Lambda_lower**2 + .5 * SQRT(4. * (j_lower+0.5)**2 &
                      - 4. * Y_lower * Lambda_lower**2 + (Y_lower * Lambda_lower)**2 ))
                  END IF

                  nubar = nubar0 + (F_upper - F_lower) / (100.*PlanckConst*c)
                  r2    = EXP( (F_upper-F_lower) / (BoltzmannConst*RadiationInput(iSpec)%Trot) ) * r1

                ! --- set strength factor for different branches
                  u = 1. / SQRT(Y_transition**2 - 4. * Y_transition + (2.*j_transition+1.)**2)
                  SELECT CASE(iBranch)
                    CASE(1)
                      HoenlLondon = ((2. * j_transition + 1.)**2 + (2.*j_transition+1.) * u &
                          * (4.*j_transition**2 + 4.*j_transition + 1. - 2.*Y_transition) ) / (32. * (j_transition+1.))
                    CASE(2)
                      HoenlLondon = ((2. * j_transition + 1.)**2 + (2.*j_transition+1.) * u &
                          * (4.*j_transition**2 + 4.*j_transition + 1. - 2.*Y_transition) ) / (32. * (j_transition+0.))
                    CASE(3)
                      HoenlLondon = ((2. * j_transition + 1.)**2 - (2.*j_transition+1.) * u &
                          * (4.*j_transition**2 + 4.*j_transition + 1. - 2.*Y_transition) ) / (32. * (j_transition+1.))
                    CASE(4)
                      HoenlLondon = ((2. * j_transition + 1.)**2 - (2.*j_transition+1.) * u &
                          * (4.*j_transition**2 + 4.*j_transition + 1. - 2.*Y_transition) ) / (32. * (j_transition+0.))

                    CASE(5) !double branch
                      HoenlLondon = ((2. * j_transition + 1.)**2 - (2.*j_transition+1.) * u &
                          * (4.*j_transition**2 + 4.*j_transition - 7. + 2.*Y_transition) ) / (32. * (j_transition+1.))
                      HoenlLondon = HoenlLondon &
                          + (2. * j_transition + 1.) * ((4.*j_transition**2 + 4.*j_transition - 1.) + u &
                          * (8.*j_transition**3 + 12.*j_transition**2 - 2.*j_transition + 1. - 2.*Y_transition)) &
                          / (32. * j_transition * (j_transition+1.))
                    CASE(6) !double branch
                      HoenlLondon = ((2. * j_transition + 1.)**2 - (2.*j_transition+1.) * u &
                          * (4.*j_transition**2 + 4.*j_transition - 7. + 2.*Y_transition) ) / (32. * (j_transition+0.))
                      HoenlLondon = HoenlLondon &
                          + (2. * j_transition + 1.) * ((4.*j_transition**2 + 4.*j_transition - 1.) + u &
                          * (8.*j_transition**3 + 12.*j_transition**2 - 2.*j_transition - 7. + 2.*Y_transition)) &
                          / (32. * j_transition * (j_transition+1.))
                    CASE(7) !double branch
                      HoenlLondon = ((2. * j_transition + 1.)**2 + (2.*j_transition+1.) * u &
                          * (4.*j_transition**2 + 4.*j_transition - 7. + 2.*Y_transition) ) / (32. * (j_transition+0.))
                      HoenlLondon = HoenlLondon &
                          + (2. * j_transition + 1.) * ((4.*j_transition**2 + 4.*j_transition - 1.) - u &
                          * (8.*j_transition**3 + 12.*j_transition**2 - 2.*j_transition - 7. + 2.*Y_transition)) &
                          / (32. * j_transition * (j_transition+1.))
                    CASE(8) !double branch
                      HoenlLondon = ((2. * j_transition + 1.)**2 + (2.*j_transition+1.) * u &
                          * (4.*j_transition**2 + 4.*j_transition - 7. + 2.*Y_transition) ) / (32. * (j_transition+1.))
                      HoenlLondon = HoenlLondon &
                          + (2. * j_transition + 1.) * ((4.*j_transition**2 + 4.*j_transition - 1.) - u &
                          * (8.*j_transition**3 + 12.*j_transition**2 - 2.*j_transition + 1. - 2.*Y_transition)) &
                          / (32. * j_transition * (j_transition+1.))
                  END SELECT

                ! --- determine if lines alternate in intesity -> evaluate the alternation factor for homocuclear molecules
                  IF(SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),12) .EQ. 0.0) THEN
                    AlternationFactor = 1.0
                  ELSE
                    AlternationFactor = 1.0 + (-1.0)**NINT(k_lower &
                      + SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),12)) &
                      / (2.0*SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),11) + 1.0) !Radiation: kexp=INT(jl+altnat(uplev)+0.1)
                  END IF

                  eps = 1.0E16 * AlternationFactor * MAX(0.,HoenlLondon) * (nubar**2 * a0e)**2 * cint2 &
                    * EXP(-(cint1 + (Bv_upper * j_upper * (j_upper + 1.)) / RadiationInput(iSpec)%Trot) / BoltzmannConst)
                  abstot = (nubar**(-5)) * (r2 - 1.0) * eps / (2.E10 * PlanckConst * c**2)
                  etot = etot + eps

                ! --- calculate line profile and add radiative energy to emission
                  IF(eps .GE. 1.0E-25) THEN
                    CALL Radiation_Molecular_Transition_Line_Profile(Radiation_Profile, nubar, Voigt_int, Voigt_int_distmax, &
                      epsilon_mol, epsilon_iSpec, eps, Radiation_Absorption_spec, abs_iSpec, abstot, iElem)
                  END IF

                  j_lower      = j_lower + 1.0
                  j_upper      = j_upper + 1.0
                  j_transition = j_transition + 1.0

                END DO !k
              END DO !iBranch

            CASE(4) ! 2Sigma -> 2Sigma

            ! --- computing Hoenl-London factors for 2Sigma->2Sigma transitions for diatomic molecules
              DO k = 0, 500
              ! --- Branch 01: R1 branch
                j_upper = REAL(k + 1.)
                j_lower = REAL(k     )
                HoenlLondon_2s2s(k,1) = (j_lower + 1.) * (j_lower + 2.) / (2. * j_lower + 3.)
              ! --- Branch 02: R2 branch
                j_upper = REAL(k + 1.)
                j_lower = REAL(k     )
                HoenlLondon_2s2s(k,2) = j_lower * (j_lower + 1.) / (2. * j_lower + 1.)
              ! --- Branch 03: P1 branch
                j_upper = REAL(k - 1.)
                j_lower = REAL(k     )
                HoenlLondon_2s2s(k,3) = j_lower * (j_lower - 1.) / (2. * j_lower - 1.)
              ! --- Branch 04: P2 branch
                j_upper = REAL(k - 1.)
                j_lower = REAL(k     )
                HoenlLondon_2s2s(k,4) = j_lower * (j_lower + 1.) / (2. * j_lower + 1.)
              ! --- Branch 05: RG21 branch
                j_upper = REAL(k + 1.)
                j_lower = REAL(k     )
                HoenlLondon_2s2s(k,5) = j_lower / ( (2. * j_lower - 1.) * (2. * j_lower + 1.) )
              ! --- Branch 06: RQ21 Satelliten-Zweig
                j_upper = REAL(k - 1.)
                j_lower = REAL(k     )
                HoenlLondon_2s2s(k,6) = (j_lower + 1.) / ( (2. * j_lower + 1.) * (2. * j_lower + 3.) )
              END DO

              Lambda_upper = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),17)
              Lambda_lower = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),17)
              Y_upper = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),16) / Bv_upper
              Y_lower = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),16) / Bv_lower

              etot = 0.0
  !PRINT*, kmax
              DO iBranch = 1, 6
                DO k = 0, kmax
                  j_upper = REAL(k + Deltak_upper_2s2s(iBranch))
                  j_lower = REAL(k + Deltak_lower_2s2s(iBranch))

                  IF((iBranch .EQ. 1) .OR. (iBranch .EQ. 3)) THEN
                    F_upper = Bv_upper * ( (j_upper+0.5)**2 - Lambda_upper**2 - .5 * SQRT(4. * (j_upper+0.5)**2 &
                      + Y_upper * (Y_upper-4.) * Lambda_upper**2) ) - Dv_upper * (j_upper + 0.)**4
                    F_lower = Bv_lower * ( (j_lower+0.5)**2 - Lambda_lower**2 - .5 * SQRT(4. * (j_lower+0.5)**2 &
                      + Y_lower * (Y_lower-4.) * Lambda_lower**2) ) - Dv_lower * (j_lower + 0.)**4
                  ELSEIF ((iBranch .EQ. 2) .OR. (iBranch .EQ. 4)) THEN
                    F_upper = Bv_upper * ( (j_upper+0.5)**2 - Lambda_upper**2 + .5 * SQRT(4. * (j_upper+0.5)**2 &
                      + Y_upper * (Y_upper-4.) * Lambda_upper**2) ) - Dv_upper * (j_upper + 1.)**4
                    F_lower = Bv_lower * ( (j_lower+0.5)**2 - Lambda_lower**2 + .5 * SQRT(4. * (j_lower+0.5)**2 &
                      + Y_lower * (Y_lower-4.) * Lambda_lower**2) ) - Dv_lower * (j_lower + 1.)**4
                  ELSEIF (iBranch .EQ. 5) THEN
                    F_upper = Bv_upper * ( (j_upper+0.5)**2 - Lambda_upper**2 + .5 * SQRT(4. * (j_upper+0.5)**2 &
                      + Y_upper * (Y_upper-4.) * Lambda_upper**2) ) - Dv_upper * (j_upper + 1.)**4
                    F_lower = Bv_lower * ( (j_lower+0.5)**2 - Lambda_lower**2 - .5 * SQRT(4. * (j_lower+0.5)**2 &
                      + Y_lower * (Y_lower-4.) * Lambda_lower**2) ) - Dv_lower * (j_lower + 0.)**4
                  ELSEIF (iBranch .EQ. 6) THEN
                    F_upper = Bv_upper * ( (j_upper+0.5)**2 - Lambda_upper**2 - .5 * SQRT(4. * (j_upper+0.5)**2 &
                      + Y_upper * (Y_upper-4.) * Lambda_upper**2) ) - Dv_upper * (j_upper + 0.)**4
                    F_lower = Bv_lower * ( (j_lower+0.5)**2 - Lambda_lower**2 + .5 * SQRT(4. * (j_lower+0.5)**2 &
                      + Y_lower * (Y_lower-4.) * Lambda_lower**2) ) - Dv_lower * (j_lower + 1.)**4
                  END IF

                  nubar = nubar0 + (F_upper - F_lower) / (100.*PlanckConst*c)
                  r2    = EXP( (F_upper-F_lower) / (BoltzmannConst*RadiationInput(iSpec)%Trot) ) * r1

                ! --- evaluate the alternation factor for homocuclear molecules
                  IF(SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),12) .EQ. 0.0) THEN
                    AlternationFactor = 1.0
                  ELSE
                    AlternationFactor = 1.0 + (-1.0)**NINT(j_lower &
                      + SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),12)) &
                      / (2.0*SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),11) + 1.0) !Radiation: kexp=INT(jl+altnat(uplev)+0.1)
                  END IF

                  eps = 1.0E16 * AlternationFactor * HoenlLondon_2s2s(k,iBranch) * (nubar**2 * a0e)**2 * cint2 &
                      * EXP(-(cint1 + (Bv_upper * j_upper * (j_upper + 1.)) / RadiationInput(iSpec)%Trot) / BoltzmannConst)
                  abstot = (nubar**(-5)) * (r2 - 1.0) * eps / (2.E10 * PlanckConst * c**2)
                  etot = etot + eps

                ! --- calculate line profile and add radiative energy to emission
                  IF(eps .GE. 1.0E-25) THEN
                    CALL Radiation_Molecular_Transition_Line_Profile(Radiation_Profile, nubar, Voigt_int, Voigt_int_distmax, &
                      epsilon_mol, epsilon_iSpec, eps, Radiation_Absorption_spec, abs_iSpec, abstot, iElem)
                  END IF
                END DO ! k
              END DO ! iBranch

            CASE(5) ! 2Pi -> 2Pi

            ! --- computing Hoenl-London factors for 2Pi->2Pi transitions for diatomic molecules
              DO k = 3, 500
              ! --- Branch 01: P1 branch
                j_upper = REAL(k - 1.) + .5
                j_lower = REAL(k     ) + .5
                HoenlLondon_2p2p(k,1) = (j_lower-1.5) * (j_lower+1.5) / (4.*j_lower*ckovacs(j_upper-1.0,+1.0) &
                  * ckovacs(j_lower,+1.) ) * (ukovacs(j_upper-1.,+1.) * ukovacs(j_lower,+1.) + 4. * (j_lower-.5) * (j_lower+.5 ) )
              ! --- Branch 02: Q1 branch
                j_upper = REAL(k     ) + .5
                j_lower = REAL(k     ) + .5
                HoenlLondon_2p2p(k,2) = (j_lower+.5) / (2. * j_lower*(j_lower+1.) * ckovacs(j_upper,+1.) * ckovacs(j_lower,+1.) ) &
                  * ( 1.5 * ukovacs(j_upper,+1.) * ukovacs(j_lower,+1.) + 2.*(j_lower-0.5 ) * (j_lower+1.5) )**2
              ! --- Branch 03: R1 branch
                j_upper = REAL(k + 1.) + .5
                j_lower = REAL(k     ) + .5
                HoenlLondon_2p2p(k,3) = (j_lower-.5) * (j_lower+2.5) / (4. * (j_lower+1.) * ckovacs(j_upper+1.,+1.) &
                  * ckovacs(j_lower,+1.) ) * (ukovacs(j_upper+1.,+1.) * ukovacs(j_lower,+1.) + 4. * (j_lower+.5) * (j_lower+1.5) )**2
              ! --- Branch 04: QP21 branch
                j_upper = REAL(k     ) + .5
                j_lower = REAL(k - 1.) + .5
                HoenlLondon_2p2p(k,4) = (j_lower-1.5) * (j_lower+1.5) / (4. * j_lower * ckovacs(j_upper-1.,-1.) &
                  * ckovacs(j_lower,+1.) ) * (ukovacs(j_upper-1.,-1.) * ukovacs(j_lower,+1.) - 4. * (j_lower-0.5) * (j_lower+0.5) )**2
              ! --- Branch 05: RQ21 branch
                j_upper = REAL(k     ) + .5
                j_lower = REAL(k     ) + .5
                HoenlLondon_2p2p(k,5) = (j_lower + 0.5) / (2. * j_lower * (j_lower+1.) * ckovacs(j_upper,-1.) &
                  * ckovacs(j_lower,+1.) ) * (1.5*ukovacs(j_upper,-1.) * ukovacs(j_lower,+1.) - 2. * (j_lower-.5) * (j_lower+1.5) )**2
              ! --- Branch 06: SR21 branch
                j_upper = REAL(k + 1.) + .5
                j_lower = REAL(k     ) + .5
                HoenlLondon_2p2p(k,6) = (j_lower-.5) * (j_lower+2.5) / (4. * (j_lower+1.) * ckovacs(j_upper+1.,-1.) &
                  * ckovacs(j_lower,+1.) ) * (ukovacs(j_upper+1.,-1.) * ukovacs(j_lower,+1.) - 4. * (j_lower+0.5) * (j_lower+1.5) )**2
              ! --- Branch 07: OP12 branch
                j_upper = REAL(k - 1.) + .5
                j_lower = REAL(k     ) + .5
                HoenlLondon_2p2p(k,7) = (j_lower-1.5) * (j_lower+1.5) / (4. * j_lower* ckovacs(j_upper-1.,+1.) &
                  * ckovacs(j_lower,-1.) ) * (ukovacs(j_upper-1.,+1.) * ukovacs(j_lower,-1.) - 4. * (j_lower-0.5) * (j_lower+0.5) )**2
              ! --- Branch 08: PQ12 branch
                j_upper = REAL(k     ) + .5
                j_lower = REAL(k     ) + .5
                HoenlLondon_2p2p(k,8) = (j_lower+0.5) / (2. * j_lower * (j_lower+1.) * ckovacs(j_upper,+1.) * ckovacs(j_lower,-1.) ) &
                  * (1.5 * ukovacs(j_upper,+1.) * ukovacs(j_lower,-1.) - 2. * (j_lower-.5) * (j_lower+ 1.5) )**2
              ! --- Branch 09: RQ12 branch
                j_upper = REAL(k + 1.) + .5
                j_lower = REAL(k     ) + .5
                HoenlLondon_2p2p(k,9) = (j_lower-.5) * (j_lower+2.5) / (4. * (j_lower+1.) * ckovacs(j_upper+1.,+1.) &
                  * ckovacs(j_lower,-1.) ) * (ukovacs(j_upper+1.,+1.) * ukovacs(j_lower,-1.)- 4. * (j_lower+.5) * (j_lower+1.5) )**2
              ! --- Branch 10: P2 branch
                j_upper = REAL(k - 1.) + .5
                j_lower = REAL(k     ) + .5
                HoenlLondon_2p2p(k,10) = (j_lower-1.5) * (j_lower+1.5) / (4. * j_lower * ckovacs(j_upper-1.,-1.) &
                  * ckovacs(j_lower,-1.) ) * (ukovacs(j_upper-1.,-1.) * ukovacs(j_lower,-1.) + 4. * (j_lower-.5) * (j_lower+.5) )**2
              ! --- Branch 11: Q2 branch
                j_upper = REAL(k     ) + .5
                j_lower = REAL(k     ) + .5
                HoenlLondon_2p2p(k,11) = (j_lower+.5) / (2. * j_lower * (j_lower+1.) * ckovacs(j_upper,-1.) * ckovacs(j_lower,-1.) ) &
                  * (1.5 * ukovacs(j_upper,-1.) * ukovacs(j_lower,-1.) + 2. * (j_lower-.5) * (j_lower+1.5) )**2
              ! --- Branch 12: R2 branch
                j_upper = REAL(k + 1.) + .5
                j_lower = REAL(k     ) + .5
                HoenlLondon_2p2p(k,12) = (j_lower-.5) * (j_lower+2.5) / (4. * (j_lower+1.) * ckovacs(j_upper+1.,-1.) &
                  * ckovacs(j_lower,-1.) ) * (ukovacs(j_upper+1.,-1.) * ukovacs(j_lower,-1.) + 4. * (j_lower-.5) * (j_lower+1.5) )**2
              END DO

              Lambda_upper = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),17)
              Lambda_lower = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),17)
              Y_upper = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),16) / Bv_upper
              Y_lower = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),16) / Bv_lower

              etot = 0.0

              DO iBranch = 1, 12
                DO k = 3, kmax
                  j_upper = REAL(k + Deltak_upper_2p2p(iBranch)) + .5
                  j_lower = REAL(k + Deltak_lower_2p2p(iBranch)) + .5

                  IF((iBranch .EQ. 1) .OR. (iBranch .EQ. 2) .OR. (iBranch .EQ. 3)) THEN
                    F_upper = Bv_upper * ( (j_upper+0.5)**2 - Lambda_upper**2 - 0.5 * SQRT(4. * (j_upper+0.5)**2 &
                      + Y_upper * (Y_upper-4.) * Lambda_upper**2) ) - Dv_upper * (j_upper + 0.)**4
                    F_lower = Bv_lower * ( (j_lower+0.5)**2 - Lambda_lower**2 - 0.5 * SQRT(4. * (j_lower+0.5)**2 &
                      + Y_lower * (Y_lower-4.) * Lambda_lower**2) ) - Dv_lower * (j_lower + 0.)**4
                  ELSEIF ((iBranch .EQ. 4) .OR. (iBranch .EQ. 5) .OR. (iBranch .EQ. 6)) THEN
                    F_upper = Bv_upper * ( (j_upper+0.5)**2 - Lambda_upper**2 + 0.5 * SQRT(4. * (j_upper+0.5)**2 &
                      + Y_upper * (Y_upper-4.) * Lambda_upper**2) ) - Dv_upper * (j_upper + 1.)**4
                    F_lower = Bv_lower * ( (j_lower+0.5)**2 - Lambda_lower**2 - 0.5 * SQRT(4. * (j_lower+0.5)**2 &
                      + Y_lower * (Y_lower-4.) * Lambda_lower**2) ) - Dv_lower * (j_lower + 0.)**4
                  ELSEIF ((iBranch .EQ. 7) .OR. (iBranch .EQ. 8) .OR. (iBranch .EQ. 9)) THEN
                    F_upper = Bv_upper * ( (j_upper+0.5)**2 - Lambda_upper**2 - 0.5 * SQRT(4. * (j_upper+0.5)**2 &
                      + Y_upper * (Y_upper-4.) * Lambda_upper**2) ) - Dv_upper * (j_upper + 0.)**4
                    F_lower = Bv_lower * ( (j_lower+0.5)**2 - Lambda_lower**2 + 0.5 * SQRT(4. * (j_lower+0.5)**2 &
                      + Y_lower * (Y_lower-4.) * Lambda_lower**2) ) - Dv_lower * (j_lower + 1.)**4
                  ELSEIF ((iBranch .EQ. 10) .OR. (iBranch .EQ. 11) .OR. (iBranch .EQ. 12)) THEN
                    F_upper = Bv_upper * ( (j_upper+0.5)**2 - Lambda_upper**2 + 0.5 * SQRT(4. * (j_upper+0.5)**2 &
                      + Y_upper * (Y_upper-4.) * Lambda_upper**2) ) - Dv_upper * (j_upper + 1.)**4
                    F_lower = Bv_lower * ( (j_lower+0.5)**2 - Lambda_lower**2 + 0.5 * SQRT(4. * (j_lower+0.5)**2 &
                      + Y_lower * (Y_lower-4.) * Lambda_lower**2) ) - Dv_lower * (j_lower + 1.)**4
                  END IF

                  nubar = nubar0 + (F_upper - F_lower) / (100.*PlanckConst*c)
                  r2    = EXP( (F_upper-F_lower) / (BoltzmannConst*RadiationInput(iSpec)%Trot) ) * r1

                ! --- evaluate the alternation factor for homocuclear molecules
                  IF(SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),12) .EQ. 0.0) THEN
                    AlternationFactor = 1.0
                  ELSE
                    AlternationFactor = 1.0 + (-1.0)**NINT(j_lower &
                      + SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),12)) &
                      / (2.0*SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),11) + 1.0) !Radiation: kexp=INT(jl+altnat(uplev)+0.1)
                  END IF

                  eps = 1.0E16 * AlternationFactor * HoenlLondon_2p2p(k,iBranch) * (nubar**2 * a0e)**2 * cint2 &
                      * EXP(-(cint1 + (Bv_upper * j_upper * (j_upper + 1.)) / RadiationInput(iSpec)%Trot) / BoltzmannConst) * 0.66 ! TODO: y 0.66?
                  abstot = (nubar**(-5)) * (r2 - 1.0) * eps / (2.E10 * PlanckConst * c**2)
                  etot = etot + eps

                ! --- calculate line profile and add radiative energy to emission
                  IF(eps .GE. 1.0E-25) THEN
                    CALL Radiation_Molecular_Transition_Line_Profile(Radiation_Profile, nubar, Voigt_int, Voigt_int_distmax, &
                      epsilon_mol, epsilon_iSpec, eps, Radiation_Absorption_spec, abs_iSpec, abstot, iElem)
                  END IF
                END DO ! k
              END DO ! iBranch

            CASE(6) ! 3Pi -> 3Pi

            ! --- computing Hoenl-London factors for 3Pi->3Pi transitions for diatomic molecules
              DO k = 0, 500
              ! --- Branch 01: R1 branch
                j_upper = REAL(k + 1.) + 1.
                j_lower = REAL(k     ) + 1.
                HoenlLondon_3p3p(k,1) = (j_lower + 1.) * (j_lower - 1.) / (3. * j_lower)
              ! --- Branch 02: P1 branch
                j_upper = REAL(k - 1.) + 1.
                j_lower = REAL(k     ) + 1.
                HoenlLondon_3p3p(k,2) = j_lower * (j_lower + 2.) / (3. * (j_lower + 1.))
              ! --- Branch 03: R2 branch
                j_upper = REAL(k + 1.) + 1.
                j_lower = REAL(k     ) + 1.
                HoenlLondon_3p3p(k,3) = (j_lower + 1.) * (j_lower - 1.) / (3. * j_lower)
              ! --- Branch 04: Q2 branch
                j_upper = REAL(k     ) + 1.
                j_lower = REAL(k     ) + 1.
                HoenlLondon_3p3p(k,4) = (2. * j_lower + 1.) / (3. * j_lower * (j_lower + 1.))
              ! --- Branch 05: P2 branch
                j_upper = REAL(k - 1.) + 1.
                j_lower = REAL(k     ) + 1.
                HoenlLondon_3p3p(k,5) = j_lower * (j_lower + 2.) / (3. * (j_lower + 1.))
              ! --- Branch 06: R3 branch
                j_upper = REAL(k + 1.) + 1.
                j_lower = REAL(k     ) + 1.
                HoenlLondon_3p3p(k,6) = (j_lower + 1.) * (j_lower - 1.) / (3. * j_lower)
              ! --- Branch 07: Q3 branch
                j_upper = REAL(k     ) + 1.
                j_lower = REAL(k     ) + 1.
                HoenlLondon_3p3p(k,7) = (2. * j_lower + 1.) / (3. * j_lower * (j_lower + 1.))
              ! --- Branch 08: P3 branch
                j_upper = REAL(k - 1.) + 1.
                j_lower = REAL(k     ) + 1.
                HoenlLondon_3p3p(k,8) = j_lower * (j_lower + 2.) / (3. * (j_lower + 1.))
              END DO

              Lambda_upper = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),17)
              Lambda_lower = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),17)
              Y_upper = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),16) / Bv_upper
              Y_lower = SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),16) / Bv_lower

              etot = 0.0

              DO iBranch = 1, 8
                DO k = 0, kmax
                  j_upper = REAL(k + Deltak_upper_3p3p(iBranch)) + 1.
                  j_lower = REAL(k + Deltak_lower_3p3p(iBranch)) + 1.
                  Z1_upper = Lambda_upper**2 * Y_upper * (Y_upper-4.) + 4./3. + 4. * j_upper * (j_upper + 1.)
                  Z1_lower = Lambda_lower**2 * Y_lower * (Y_lower-4.) + 4./3. + 4. * j_lower * (j_lower + 1.)
                  Z2_upper = (Lambda_upper**2 * Y_upper * (Y_upper-1.) - 4./9. - 2. * j_upper * (j_upper + 1.)) / (3.*Z1_upper) ! TODO: + or - 4/9?
                  Z2_lower = (Lambda_lower**2 * Y_lower * (Y_lower-1.) - 4./9. - 2. * j_lower * (j_lower + 1.)) / (3.*Z1_lower) ! TODO: + or - 4/9?

                  IF((iBranch .EQ. 1) .OR. (iBranch .EQ. 2)) THEN
                    F_upper = Bv_upper * (j_upper * (j_upper + 1.) - SQRT(Z1_upper) - 2. * Z2_upper) - Dv_upper * (j_upper - 0.5)**4
                    F_lower = Bv_lower * (j_lower * (j_lower + 1.) - SQRT(Z1_lower) - 2. * Z2_lower) - Dv_lower * (j_lower - 0.5)**4
                  ELSEIF ((iBranch .EQ. 3) .OR. (iBranch .EQ. 4) .OR. (iBranch .EQ. 5)) THEN
                    F_upper = Bv_upper * (j_upper * (j_upper + 1.)                  + 4. * Z2_upper) - Dv_upper * (j_upper + 0.5)**4
                    F_lower = Bv_lower * (j_lower * (j_lower + 1.)                  + 4. * Z2_lower) - Dv_lower * (j_lower + 0.5)**4
                  ELSEIF ((iBranch .EQ. 6) .OR. (iBranch .EQ. 7) .OR. (iBranch .EQ. 8)) THEN
                    F_upper = Bv_upper * (j_upper * (j_upper + 1.) + SQRT(Z1_upper) - 2. * Z2_upper) - Dv_upper * (j_upper + 1.5)**4
                    F_lower = Bv_lower * (j_lower * (j_lower + 1.) + SQRT(Z1_lower) - 2. * Z2_lower) - Dv_lower * (j_lower + 1.5)**4
                  END IF
                  nubar = nubar0 + (F_upper - F_lower) / (100.*PlanckConst*c)
                  r2    = EXP( (F_upper-F_lower) / (BoltzmannConst*RadiationInput(iSpec)%Trot) ) * r1

                ! --- evaluate the alternation factor for homocuclear molecules
                  IF(SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),12) .EQ. 0.0) THEN
                    AlternationFactor = 1.0
                  ELSE
                    AlternationFactor = 1.0 + (-1.0)**NINT(j_lower &
                      + SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,1),12)) &
                      / (2.0*SpeciesRadiation(iSpec)%EnergyLevelProperties(SpeciesRadiation(iSpec)%BandProperties(iBand,2),11) + 1.0) !Radiation: kexp=INT(jl+altnat(uplev)+0.1)
                  END IF

                  eps = 1.0E16 * AlternationFactor * HoenlLondon_3p3p(k,iBranch) * (nubar**2 * a0e)**2 * cint2 &
                      * EXP(-(cint1 + (Bv_upper * j_upper * (j_upper + 1.)) / RadiationInput(iSpec)%Trot) / BoltzmannConst)
                  abstot = (nubar**(-5)) * (r2 - 1.0) * eps / (2.E10 * PlanckConst * c**2)

                  etot = etot + eps

                ! --- calculate line profile and add radiative energy to emission
                  IF(eps .GE. 1.0E-25) THEN
                    CALL Radiation_Molecular_Transition_Line_Profile(Radiation_Profile, nubar, Voigt_int, Voigt_int_distmax, &
                      epsilon_mol, epsilon_iSpec, eps, Radiation_Absorption_spec, abs_iSpec, abstot, iElem)
                  END IF

                END DO !k

              END DO ! iBranch

            CASE DEFAULT
              CALL abort(&
                __STAMP__&
                ,' ERROR: Radiation transition type is not implemented! (unknown case)')
          END SELECT

        END DO

      END DO

  !    hilf = LEN_TRIM(RadiationInput(iSpec)%RadiationSpectraFileName)
  !    RadiationInput(iSpec)%RadiationSpectraFileName = RadiationInput(iSpec)%RadiationSpectraFileName(1:hilf-4)
  !    WRITE(*,*) '*** CALCULATED ', TRIM(RadiationInput(iSpec)%RadiationSpectraFileName), ' ***'
    END IF
    
    TempOut_Em  = 0.0
    TempOut_Abs = 0.0
    DO iWave=1, RadiationParameter%WaveLenDiscr
      iWaveCoarse = INT((iWave-1)/RadiationParameter%WaveLenReductionFactor) + 1
      IF (iWaveCoarse.GT.RadiationParameter%WaveLenDiscrCoarse) iWaveCoarse = RadiationParameter%WaveLenDiscrCoarse
      TempOut_Em  = TempOut_Em + 4.*Pi*epsilon_iSpec(iWave)*RadiationParameter%WaveLenIncr
      TempOut_Abs = TempOut_Abs + abs_iSpec(iWave)*RadiationParameter%WaveLenIncr
      Radiation_Absorption_SpeciesWave(iWaveCoarse, iSpec) = Radiation_Absorption_SpeciesWave(iWaveCoarse, iSpec) + abs_iSpec(iWave)*RadiationParameter%WaveLenIncr
    END DO
    Radiation_ElemEnergy_Species(iSpec,iElem,1) = TempOut_Em
    Radiation_ElemEnergy_Species(iSpec,iElem,2) = TempOut_Abs

  END DO

  DO iWave = 1, RadiationParameter%WaveLenDiscr
    iWaveCoarse = INT((iWave-1)/RadiationParameter%WaveLenReductionFactor) + 1
    IF (iWaveCoarse.GT.RadiationParameter%WaveLenDiscrCoarse) iWaveCoarse = RadiationParameter%WaveLenDiscrCoarse
    Radiation_Emission_spec(iWaveCoarse,iElem) = Radiation_Emission_spec(iWaveCoarse,iElem) + epsilon_mol(iWave)/RadiationParameter%WaveLenReductionFactor
  END DO

! --- add contribution to total emission
  em_mol = epsilon_mol(1) * RadiationParameter%WaveLenIncr
  DO i=2, RadiationParameter%WaveLenDiscr
    em_mol = em_mol + epsilon_mol(i) * RadiationParameter%WaveLenIncr
  END DO


  ! WRITE(*,*) '*** MOLECULAR BOUND-BOUND RADIATION SUCCESSFULLY DONE***'
  ! WRITE(*,*) ''

END SUBROUTINE radiation_molecules




REAL FUNCTION ukovacs(rj, dlam)
!===================================================================================================================================
! function of kovacs to derive the Hoenl-London factors
!===================================================================================================================================
! MODULES

! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT(IN)                :: rj, dlam
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  
  ukovacs = 2. * (rj + dlam + .5)
  RETURN

END FUNCTION ukovacs




REAL FUNCTION ckovacs(rj, dlam)
!===================================================================================================================================
! function of kovacs to derive the Hoenl-London factors
!===================================================================================================================================
! MODULES

! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT(IN)                :: rj, dlam
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  
  ckovacs = 4. * (rj + .5) * (rj + dlam + .5)
  RETURN

END FUNCTION ckovacs




SUBROUTINE Radiation_Molecular_Transition_Line_Profile(Radiation_Profile, nubar, Voigt_int, Voigt_int_distmax, &
    epsilon_mol, epsilon_iSpec, eps, Radiation_Absorption_spec, abs_iSpec, abstot, iElem)
!===================================================================================================================================
! calculates emission profile functions and adds radiative energy to emission array
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_Radiation_Vars,    ONLY   : RadiationParameter
  USE MOD_Mesh_Tools,            ONLY : GetGlobalElemID
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT(INOUT)             :: Radiation_Profile(:), epsilon_mol(:), epsilon_iSpec(:), Radiation_Absorption_spec(:,:), &
                                     abs_iSpec(:)
  REAL, INTENT(IN)                :: nubar, Voigt_int_distmax,  eps, abstot,Voigt_int(-1001:1001)
  INTEGER, INTENT(IN)             :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER           :: startwavelength_int, endwavelength_int, i, iWave, iWaveCoarse
  LOGICAL           :: add_radline
!===================================================================================================================================

! --- determination of indices for first and last entry of Voigt-profiles on wavelength axis
  startwavelength_int = MAX(1, INT(((1./(100.*nubar)) - Voigt_int_distmax - RadiationParameter%MinWaveLen) &
                      / RadiationParameter%WaveLenIncr) + 1)
  endwavelength_int   = MIN(RadiationParameter%WaveLenDiscr, INT(((1./(100.*nubar)) + Voigt_int_distmax &
      - RadiationParameter%MinWaveLen) / RadiationParameter%WaveLenIncr) + 1)

  add_radline          = .FALSE.

  IF ( (startwavelength_int .LT. endwavelength_int) .AND. &
          (((1./(100.*nubar))+Voigt_int_distmax) .GT. RadiationParameter%MinWaveLen) .AND. &
          (((1./(100.*nubar))-Voigt_int_distmax) .LT. RadiationParameter%MaxWaveLen)) THEN
!                  IF ( (startwavelength_int .LT. endwavelength_int) .AND. &
!                          (((1./(100.*nubar))+Voigt_int_distmax) .GT. RadiationParameter%MinWaveLen) ) THEN

    add_radline = .TRUE.

    ! --- determine transition lines of previous computed Voigt-profiles
      CALL Radiation_Voigt_wavelength_interpolation(Voigt_int, Voigt_int_distmax, (1./(100.*nubar)), &
        Radiation_Profile, startwavelength_int, endwavelength_int)

    IF (add_radline) THEN

!      istart = MIN(istart,startwavelength_int) ! TODO: ???
!      iend   = MAX(iend,  endwavelength_int)

    ! --- add radiative energy to emission
      DO i=0  , (endwavelength_int-startwavelength_int + 0)
        epsilon_iSpec(startwavelength_int+i) &
          = epsilon_iSpec(startwavelength_int+i) + MAX(0.0,eps)/1.E10 * Radiation_Profile(startwavelength_int+i)
        epsilon_mol(startwavelength_int+i) &
          = epsilon_mol(startwavelength_int+i) + MAX(0.0,eps)/1.E10 * Radiation_Profile(startwavelength_int+i)
        abs_iSpec(startwavelength_int+i) &
          = abs_iSpec(startwavelength_int+i)+MAX(0.0,abstot)/1.E10*Radiation_Profile(startwavelength_int+i)
        iWave = startwavelength_int+i
        iWaveCoarse = INT((iWave-1)/RadiationParameter%WaveLenReductionFactor) + 1
        IF (iWaveCoarse.GT.RadiationParameter%WaveLenDiscrCoarse) iWaveCoarse = RadiationParameter%WaveLenDiscrCoarse
        Radiation_Absorption_spec(iWaveCoarse,GetGlobalElemID(iElem)) &
          = Radiation_Absorption_spec(iWaveCoarse,GetGlobalElemID(iElem))+MAX(0.0,abstot)/1.E10*Radiation_Profile(iWave)/RadiationParameter%WaveLenReductionFactor
      END DO

    END IF

  END IF

END SUBROUTINE Radiation_Molecular_Transition_Line_Profile




SUBROUTINE Radiation_Voigt_wavelength_interpolation(Voigt_int, Voigt_int_distmax, centerwavelength, &
    Radiation_Profile, startwavelength_int, endwavelength_int)
!===================================================================================================================================
! distributes transition lines with precomputed integrated Voigt-profiles (Voigt_int)
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_Radiation_Vars,       ONLY   : RadiationParameter
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT(IN)                :: Voigt_int(-1001:1001), centerwavelength, Voigt_int_distmax
                                      ! TODO Voigt_int_nmax in Voigt_int(:)
  REAL, INTENT(INOUT)             :: Radiation_Profile(:)
  INTEGER, INTENT(INOUT)          :: startwavelength_int, endwavelength_int
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER, PARAMETER  :: Voigt_int_nmax = 1001
  INTEGER             :: i!, width_int
  REAL                :: startvalue, Voigt_int_value, Voigt_int_value_old
  REAL                :: weighting_int_left, weighting_int_right
  INTEGER             :: index_left,index_right
  !REAL, ALLOCATABLE   :: help_Radiation_Profile(:)
!===================================================================================================================================

! --- extrapolation of left point
  startvalue          = Voigt_int(-Voigt_int_nmax) &
      - (Voigt_int(-Voigt_int_nmax+1)-Voigt_int(-Voigt_int_nmax))/RadiationParameter%WaveLenIncr &
      * (centerwavelength - Voigt_int_distmax - RadiationParameter%WaveLen(startwavelength_int))
  IF (startvalue .LT. 0.0) startvalue = 0.0


! --- determination of number of wavelength array indices within the corresponding Voigt profile width

!  width_int = endwavelength_int - startwavelength_int

!  IF (width_int .EQ. 0.0) THEN
!    WRITE(*,'(A,F8.4,A)') 'Warning: line at ', centerwavelength*1.E9, &
!    'nm is neglected due to a too small wavelength discretization'
!  END IF

! --- determination of transition lines

  Radiation_Profile(startwavelength_int) = startvalue
!  Radiation_Profile(endwavelength_int+1)   = 1.0

  Voigt_int_value_old = 0.0!startvalue
  Voigt_int_value = 0.0

  DO i=1, (endwavelength_int-startwavelength_int)

    index_left = INT(SIGN(ABS((RadiationParameter%WaveLen(startwavelength_int+i)-centerwavelength)) &
                 / (Voigt_int_distmax/REAL(Voigt_int_nmax)), &
                 RadiationParameter%WaveLen(startwavelength_int+i)-centerwavelength))

    IF ((RadiationParameter%WaveLen(startwavelength_int+i)-centerwavelength).LT.0.0) THEN
      index_right = index_left
      index_left  = index_left - 1
    ELSE
      index_right = index_left + 1
    END IF

    IF(index_right.GT.Voigt_int_nmax) CYCLE
!    IF(index_left.LT.(-Voigt_int_nmax)) CYCLE

    weighting_int_left  = 1. - ABS(REAL(index_left)-((RadiationParameter%WaveLen(startwavelength_int+i)-centerwavelength) &
                 / (Voigt_int_distmax/REAL(Voigt_int_nmax))))
    weighting_int_right = 1. - weighting_int_left

    Voigt_int_value = weighting_int_left * Voigt_int(index_left) + weighting_int_right * Voigt_int(index_right)

    Radiation_Profile(startwavelength_int + i) = Voigt_int_value - Voigt_int_value_old

    Voigt_int_value_old = Voigt_int_value

  END DO

!  DO i=1, width_int+1
  DO i=1, (endwavelength_int-startwavelength_int)
    Radiation_Profile(startwavelength_int+i) = Radiation_Profile(startwavelength_int+i) / RadiationParameter%WaveLenIncr
  END DO

END SUBROUTINE Radiation_Voigt_wavelength_interpolation




END MODULE MOD_Radiation_Molecules
