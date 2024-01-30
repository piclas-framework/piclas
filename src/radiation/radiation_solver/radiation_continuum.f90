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

MODULE MOD_Radiation_Continuum
!===================================================================================================================================
! Module for Radiation
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE radiation_continuum
  MODULE PROCEDURE radiation_continuum
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: radiation_continuum
!===================================================================================================================================

CONTAINS


SUBROUTINE radiation_continuum(iElem, em_cont)
!===================================================================================================================================
! Main routine of continuum radiation calculation
!===================================================================================================================================
! MODULES
  USE MOD_Radiation_Vars,    ONLY   : RadiationInput, RadiationParameter, Radiation_Emission_spec, RadiationSwitches
  USE MOD_PARTICLE_Vars,     ONLY   : nSpecies

! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT(INOUT)             :: em_cont
  INTEGER, INTENT(IN)             :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

  REAL              :: n_atom, n_ion
  INTEGER           :: jSpec, i
  REAL, ALLOCATABLE :: epsilon_cont(:)

!===================================================================================================================================
  
  ALLOCATE(epsilon_cont(RadiationParameter%WaveLenDiscr))

  n_atom = 0.0
  n_ion  = 0.0

  DO jSpec = 1, nSpecies
    IF (RadiationInput(jSpec)%NuclCharge .EQ. 1) n_atom = n_atom + RadiationInput(jSpec)%NumDens
    IF (RadiationInput(jSpec)%NuclCharge .EQ. 2) n_ion  = n_ion  + RadiationInput(jSpec)%NumDens
  END DO

! --- initialize emission coefficient
  DO i=1, RadiationParameter%WaveLenDiscr
    epsilon_cont(i) = Radiation_Emission_spec(i,iElem)  ! TODO: or 0.0?
  END DO

! --- free-free emission
  IF (RadiationSwitches%ff) THEN
! --- free-free emission due to neutrals
    CALL Radiation_continuum_ff(n_atom, 1, iElem)
! --- free-free emission due to ions
    CALL Radiation_continuum_ff(n_ion, 2, iElem)
    ! WRITE(*,*) '*** FREE-FREE CONTINUUM RADIATION SUCCESSFULLY DONE  ***'
    ! WRITE(*,*) ''
  END IF

! --- bound-free emission
  IF (RadiationSwitches%bf) THEN
    CALL Radiation_continuum_bf(iElem)
    ! WRITE(*,*) '*** BOUND-FREE CONTINUUM RADIATION SUCCESSFULLY DONE ***'
    ! WRITE(*,*) ''
  END IF


! --- calculate emission due to continuum radiation
  DO i = 1, RadiationParameter%WaveLenDiscr
    epsilon_cont(i) = Radiation_Emission_spec(i,iElem) - epsilon_cont(i)
  END DO

! --- determine total volumetric emission for continua
  CALL Radiation_continuum_total(epsilon_cont, em_cont)


END SUBROUTINE radiation_continuum









SUBROUTINE Radiation_continuum_ff(n, z, iElem)
!===================================================================================================================================
! determines free-free (ff) continuum radiation (bremsstrahlung)
! T_elec needs to be changed with PIC interface!!!
!===================================================================================================================================
! MODULES
  USE MOD_Globals_Vars,      ONLY   : BoltzmannConst, PlanckConst
  USE MOD_Radiation_Vars,    ONLY   : RadiationInput, RadiationParameter, Radiation_Emission_spec, Radiation_Absorption_spec
  USE MOD_Globals_Vars,      ONLY   : c

! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT(IN)                :: n
  INTEGER, INTENT(IN)             :: z, iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL, PARAMETER    :: gaunt_ff = 1.!0.8                         ! Gaunt factor for free-free radiation
  REAL               :: kbTe, em_brems_const, expon, em_brems  ! kbTe, frequency independent factor, exponent, free-free emission coefficient at particular wavelength
  REAL               :: B_lambda                               ! Planck's function formulated for per wavelength
  INTEGER            :: i
!===================================================================================================================================
  
  kbTe = BoltzmannConst * RadiationInput(1)%Telec
  em_brems_const = 5.443e-52*(z-1)**2*RadiationInput(1)%NumDens*n/sqrt(RadiationInput(1)%Telec)

  DO i=1, RadiationParameter%WaveLenDiscr
    expon = MIN(700., PlanckConst*c/(RadiationParameter%WaveLen(i)*kbTe))
    em_brems = gaunt_ff*em_brems_const/(RadiationParameter%WaveLen(i)**2/c)*EXP(-expon)

!    B_lambda = 2.*PlanckConst*(c/RadiationParameter%WaveLen(i))**3/(c**2)/(RadiationParameter%WaveLen(i)**2/c)/(EXP(expon)-1.)
    B_lambda = 2.*PlanckConst*c**2/(RadiationParameter%WaveLen(i)**5)/(EXP(expon)-1.)
    
    Radiation_Emission_spec(i,iElem)   = Radiation_Emission_spec(i,iElem)   + em_brems
    Radiation_Absorption_spec(i,iElem) = Radiation_Absorption_spec(i,iElem) + em_brems/B_lambda
  END DO


END SUBROUTINE Radiation_continuum_ff









SUBROUTINE Radiation_continuum_bf(iElem)
!===================================================================================================================================
! determines bound-free (bf) continuum radiation due to recombination/photoionization
!===================================================================================================================================
! MODULES
  USE MOD_Globals_Vars,      ONLY   : BoltzmannConst, PlanckConst, ElementaryCharge, ElectronMass, BohrRadius, Pi
  USE MOD_Particle_Vars,     ONLY   : nSpecies, Species
  USE MOD_Radiation_Vars,    ONLY   : RadiationInput, SpeciesRadiation, RadiationParameter, &
                                      Radiation_Emission_spec, Radiation_Absorption_spec, NumDensElectrons
  USE MOD_Globals_Vars,      ONLY   : c
  USE MOD_DSMC_Vars,         ONLY   : SpecDSMC

! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)             :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL               :: kbTe, nue                             ! kbTe, nue
  REAL               :: AdvanceSeriesLimits                   ! Advance of series limits
  REAL               :: NumDensIon, geIon
  INTEGER            :: iAtom, jAtom, iIon, iLevel, iWave     ! counter for loop over all atoms, ions, Levels, wavelength increments
  INTEGER            :: LevelsConsidered                      ! number of considered energy levels
  REAL               :: gaunt_bf                              ! Gaunt factor for bound-free radiation
  REAL               :: low_IonizationPot                     ! lowering of ionization Potenzial [J]
  REAL               :: ActualIonizationEn                    ! lowered ionization energy due to photo-ionization edge shift
  REAL               :: kappaBFFactor, epsilonBFFactor        ! Factor for the calculation of emission and absorption coefficients (constant for each atom)
  REAL               :: epsilonBF, kappaBF
  REAL, PARAMETER    :: Ryd = 13.59983                        ! [eV]
  REAL               :: tau, s, sigmaN, sigma
!===================================================================================================================================

  kbTe = BoltzmannConst * RadiationInput(1)%Telec

! --- loop over all atoms
  DO iAtom = 1, nSpecies
    IF((SpecDSMC(iAtom)%InterID .NE. 1) .AND. (SpecDSMC(iAtom)%InterID .NE. 10)) CYCLE

!    IF (RadiationInput(iAtom)%NuclCharge .GT. 1) CYCLE  !approach only for neutral atoms
    iIon = iAtom
! --- determine ionized species
    DO jAtom = 1, nSpecies ! TODO nAtoms instead of nSpecies
      IF((SpecDSMC(jAtom)%InterID .NE. 1) .AND. (SpecDSMC(jAtom)%InterID .NE. 10)) CYCLE
      IF((Species(iAtom)%MassIC .EQ. Species(jAtom)%MassIC) &
        .AND. (RadiationInput(iAtom)%NuclCharge+1 .EQ. RadiationInput(jAtom)%NuclCharge) ) THEN
        iIon = jAtom
      END IF
    END DO

! --- INITIALIZE LEVEL INDEPENDENT VARIABLES
! --- advance of series limits (Griem 5-47 E_H -> E_Ion)
    AdvanceSeriesLimits = 4. * RadiationInput(iAtom)%NuclCharge**(4./5.) * (BohrRadius**3. * NumDensElectrons)**(4./15.) &
                          * RadiationInput(iAtom)%IonizationEn
    low_IonizationPot   = 2.9E-8*SQRT(NumDensElectrons/1.E6/RadiationInput(iAtom)%Telec)*ElementaryCharge ! TODO correct indices -> To electrons
    ActualIonizationEn  = RadiationInput(iAtom)%IonizationEn - low_IonizationPot
    ActualIonizationEn  = MIN(ActualIonizationEn, RadiationInput(iAtom)%IonizationEn-AdvanceSeriesLimits)

! --- determination of levels to be considered
    DO iLevel = 1, SpeciesRadiation(iAtom)%nLevels
      IF(SpeciesRadiation(iAtom)%Level(iLevel,2) .GT. ActualIonizationEn) THEN
        LevelsConsidered = iLevel - 1
        EXIT
      END IF
    END DO
    

    IF (iIon .NE. iAtom) THEN
      geIon      = SpeciesRadiation(iIon)%Level(1,1) ! values for ion in ground state
      NumDensIon = SpeciesRadiation(iIon)%NumDensExc(1)
      ! PRINT*,'bound-free continuum is calculated using ionized species: ', &
        ! TRIM(RadiationInput(iAtom)%RadiationSpectraFileName)
    ELSE
! --- assumed degeneracy of the ion
      geIon      = 1.0d0
! --- equilibrium ground state ion number density
      NumDensIon = ((Species(iAtom)%MassIC-ElectronMass) * 2. * Pi * BoltzmannConst/PlanckConst**2. * ElectronMass &
        * RadiationInput(1)%Telec / Species(iAtom)%MassIC)**1.5d0 * geIon / SpeciesRadiation(iAtom)%Level(1,1) &
        * EXP(-ActualIonizationEn/kbTe) * SpeciesRadiation(iAtom)%NumDensExc(1) / MAX(1.0d0, NumDensElectrons)
      ! PRINT*,'bound-free continuum is approximated using equilibrium considerations! (missing data file for ionized species): ', &
        ! TRIM(RadiationInput(iAtom)%RadiationSpectraFileName)
    END IF

    
! --- emission coefficient formula is valid for photorecombination from the ion ground state only!
    epsilonBFFactor = 1.719236e-46 * NumDensElectrons * NumDensIon * RadiationInput(iAtom)%NuclCharge**4. / geIon &
      / c * ( Species(iAtom)%MassIC/((Species(iAtom)%MassIC-ElectronMass) * MAX(250., RadiationInput(iAtom)%Telec) ))**1.5
    kappaBFFactor   = 2.815401e25 * RadiationInput(iAtom)%NuclCharge**4.
    
    DO iWave = 1, RadiationParameter%WaveLenDiscr

      nue = c / RadiationParameter%WaveLen(iWave)
      
      epsilonBF = 0.
      kappaBF   = 0.
      
      DO iLevel = 1, LevelsConsidered
        IF(PlanckConst*nue .GE. ActualIonizationEn-SpeciesRadiation(iAtom)%Level(iLevel,2)) THEN
          ! gaunt factors using L.G. D'yachkov - Simple formula for the average Gaunt factor Eq. 9
          IF((SpeciesRadiation(iAtom)%Level(iLevel,5) .LT. 0.0) .OR. (SpeciesRadiation(iAtom)%Level(iLevel,5) .GT. 1.0)) THEN
            tau    = kbTe/(RadiationInput(iAtom)%NuclCharge**2.*Ryd)
            s      = 0.5+0.5*tau**(0.5) 
            sigmaN = 2.*MIN(SpeciesRadiation(iAtom)%Level(iLevel,2), RadiationInput(iAtom)%IonizationEn) &
              / (RadiationInput(iAtom)%NuclCharge**2.*Ryd)
            sigma  = SpeciesRadiation(iAtom)%Level(iLevel,2)/(RadiationInput(iAtom)%NuclCharge**2.*Ryd)
            gaunt_bf = 1. + 0.347/sigma**(2./3.)*(tau+s*sigma-sigmaN/(1.-EXP(-sigmaN/tau))) &
              - 0.0331/sigma**(4./3.)*(2.*tau*(tau+s*sigma)+(s*sigma)**2. &
              - sigmaN*(2.*tau+2.*s*sigma-sigmaN)/(1.-EXP(-sigmaN/tau)))
            SpeciesRadiation(iAtom)%Level(iLevel,5) = gaunt_bf
            PRINT*, 'Gaunt-factor for bound-free radiation calculated'
          END IF
          kappaBF   = kappaBF + SpeciesRadiation(iAtom)%Level(iLevel,5) * SpeciesRadiation(iAtom)%NumDensExc(iLevel) &
            / (nue**3. * SpeciesRadiation(iAtom)%Level(iLevel,4)**5.)
          epsilonBF = epsilonBF + SpeciesRadiation(iAtom)%Level(iLevel,5) * EXP(-SpeciesRadiation(iAtom)%Level(iLevel,2)/kbTe) &
            * nue**2. / SpeciesRadiation(iAtom)%Level(iLevel,4)**5. * SpeciesRadiation(iAtom)%Level(iLevel,1)
        END IF
      END DO

      kappaBF   = kappaBF * kappaBFFactor * (1.d0 - EXP(-PlanckConst*nue/kbTe))
      epsilonBF = epsilonBF   * epsilonBFFactor   * EXP(MIN((ActualIonizationEn-PlanckConst*nue)/kbTe, 7.d2))

!      PRINT*, kappaBF, epsilonBF

      Radiation_Emission_spec(iWave,iElem)   = Radiation_Emission_spec(iWave,iElem)   + epsilonBF
      Radiation_Absorption_spec(iWave,iElem) = Radiation_Absorption_spec(iWave,iElem) + kappaBF
      
    END DO

    ! PRINT*, 'calculated ', LevelsConsidered, 'bound-free levels of ', TRIM(RadiationInput(iAtom)%RadiationSpectraFileName)

  END DO


END SUBROUTINE Radiation_continuum_bf









SUBROUTINE Radiation_continuum_total(epsilon_cont, em_cont)
!===================================================================================================================================
! determines the total emission by numerical 2nd order integration of given epsilon with a Lagrange polynom
!===================================================================================================================================
! MODULES
  USE MOD_Radiation_Vars,       ONLY   : RadiationParameter
  USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT(IN)                :: epsilon_cont(:)
  REAL, INTENT(INOUT)             :: em_cont
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL               :: x1, x2, x3, y1, y2, y3, f1, f2, f3, a, b, c, dr_t, dr_l, dr_l2, dr, x0, y0, left, right
  INTEGER            :: i

!===================================================================================================================================
  !!!!!!!!!TODO!!!!!!!!       
  dr_l2=0.0; y0=0.0; y1=0.0; x0=0.0
  DO i=1, RadiationParameter%WaveLenDiscr-2
 
    x1 = RadiationParameter%WaveLen(i)
    x2 = RadiationParameter%WaveLen(i+1)
    x3 = RadiationParameter%WaveLen(i+2)

    y1 = epsilon_cont(i)
    y2 = epsilon_cont(i+1)
    y3 = epsilon_cont(i+2)

    f1 = y1/(x1-x2)/(x1-x3)
    f2 = y2/(x2-x1)/(x2-x3)
    f3 = y3/(x3-x1)/(x3-x2)

    a  = f1 + f2 + f3
    b  = -(x2+x3)*f1 - (x1+x3)*f2 - (x1+x2)*f3
    c  = x2*x3*f1 + x1*x3*f2 + x1*x2*f3

    left  = x1
    right = x2
    dr_l  = a/3.*(right**3.-left**3.) + b/2.*(right**2.-left**2.) + c*(right-left)
    dr_t  = 0.5 * (y1+y2) * (x2-x1)

    IF( ( dr_l .LT. MIN(y1,y2)*(x2-x1) ) .OR. ( dr_l .GT. MAX(y1,y2)*(x2-x1) ) ) THEN
      dr = dr_t
    ELSE
      dr = dr_l
    END IF

    IF( (i .GT. 1) .AND. (dr_l2 .GT. MIN(y0,y1)*(x1-x0)) .AND. (dr_l2 .LT. MAX(y0,y1)*(x1-x0) ) ) THEN
      dr = 0.5 * (dr + dr_l2)
    END IF

    em_cont = em_cont + dr

    left  = x2
    right = x3
    dr_l2 = a/3.*(right**3.-left**3.) + b/2.*(right**2.-left**2.) + c*(right-left)
    x0    = x1
    y0    = y1

  END DO

  i = RadiationParameter%WaveLenDiscr-1
  dr_t = 0.5 * (epsilon_cont(i)+epsilon_cont(i+1)) * (RadiationParameter%WaveLen(i+1)-RadiationParameter%WaveLen(i))

  em_cont = em_cont + dr_t
 

END SUBROUTINE Radiation_continuum_total



END MODULE MOD_Radiation_Continuum
