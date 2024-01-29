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

MODULE MOD_Radiation_Atoms
!===================================================================================================================================
! Module for Radiation
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE radiation_atoms
  MODULE PROCEDURE radiation_atoms
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: radiation_atoms
!===================================================================================================================================

CONTAINS


SUBROUTINE radiation_atoms(iElem, em_atom)
!===================================================================================================================================
! Main routine of atomic radiation calculation
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_Globals_Vars,      ONLY   : BoltzmannConst, PlanckConst, ElementaryCharge
  USE MOD_Radiation_Vars,    ONLY   : RadiationInput, RadiationParameter, SpeciesRadiation,  &
                                      Radiation_Emission_spec, Radiation_Absorption_spec, &
                                      NumDensElectrons, Radiation_ElemEnergy_Species, Radiation_Absorption_SpeciesWave
  USE MOD_Particle_Vars,     ONLY   : nSpecies, Species
  USE MOD_Globals_Vars,      ONLY   : c, Pi
  USE MOD_DSMC_Vars,         ONLY   : SpecDSMC
!USE MOD_Radiation_Excitation, ONLY   : low_IonizationPot

! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT(INOUT)             :: em_atom
  INTEGER, INTENT(IN)             :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

  REAL              :: c_emi, c_abs, c_dopp
  REAL              :: rho, ntot     !in ini
  REAL              :: sigma_ij, coll_freq_ij, T_mean
  INTEGER           :: i, k, l, iSpec, jSpec, iWave, iWaveCoarse!, hilf, w
  INTEGER           :: nLines_considered             !number of calculated transition lines
  REAL              :: etot, abstot
  REAL, ALLOCATABLE :: epsilon_at(:), epsilon_iSpec(:), abs_iSpec(:)
  REAL              :: cwav
  REAL              :: Dlaml, Dlamd, Dlamn, Dlamp, Dlams, Dlamr, Dlamvw !line broadening mechanisms

  REAL              :: low_IonizationPot    !from excitation
  REAL, ALLOCATABLE :: lamnu(:)
  REAL              :: wrange, wleft, wright
  REAL              :: dpar, beta, rpar, Dlamv
  INTEGER, PARAMETER:: Voigt_int_nmax = 1001
  REAL, PARAMETER   :: Voigt_mult_discr = 25.
  REAL              :: Voigt_int(-Voigt_int_nmax:Voigt_int_nmax), Voigt_int_step, Voigt
  REAL              :: Voigt_var1, Voigt_const1, Voigt_const2, Voigt_int_distmax, Voigt_arg_dimless, Voigt_int_NormEnergy
  INTEGER           :: iVoigt
  REAL              :: Voigt_dist1, Voigt_dist2, Voigt_dist225
  REAL, ALLOCATABLE :: Radiation_Profile(:)
  REAL              :: TempOut_Em, TempOut_Abs
!===================================================================================================================================

!  --- initialize emission coefficient
  ALLOCATE(epsilon_at(RadiationParameter%WaveLenDiscr))
  ALLOCATE(epsilon_iSpec(RadiationParameter%WaveLenDiscr))
  ALLOCATE(abs_iSpec(RadiationParameter%WaveLenDiscr))
  ALLOCATE(Radiation_Profile(RadiationParameter%WaveLenDiscr))
  DO iWave=1, RadiationParameter%WaveLenDiscr
    epsilon_at(iWave) = 0.0 !Radiation_Emission_spec(iWave,iElem)
  END DO

! --- loop for ntot and rho over all atoms
  rho  = 0.0
  ntot = 0.0

  DO jSpec = 1, nSpecies
    rho  = rho + RadiationInput(jSpec)%NumDens * Species(jSpec)%MassIC
    ntot = ntot + RadiationInput(jSpec)%NumDens
  END DO
  
! --- calculation of constants
    c_emi  = PlanckConst * c / (4.*Pi)
    c_abs  = 1. / (8.*Pi*c)
    c_dopp = SQRT(8.*BoltzmannConst*0.69315)/c


  DO iSpec = 1, nSpecies
    IF(.NOT.RadiationInput(iSpec)%DoRadiation) CYCLE
    IF((SpecDSMC(iSpec)%InterID .NE. 1) .AND. (SpecDSMC(iSpec)%InterID .NE. 10)) CYCLE
    IF((RadiationInput(iSpec)%Telec.LT.10.0).OR.(RadiationInput(iSpec)%NumDens.LT.10.0).OR.(RadiationInput(iSpec)%Ttrans(4).LT.10.0)) CYCLE    

    ALLOCATE(lamnu(SpeciesRadiation(iSpec)%nLines))

    lamnu               = 0.0
    Radiation_Profile   = 0.0
    epsilon_iSpec       = 0.0
    abs_iSpec       = 0.0

    IF ((SpeciesRadiation(iSpec)%nLevels.NE.0) .OR. (SpeciesRadiation(iSpec)%nLines.NE.0)) THEN

      DO l=1,SpeciesRadiation(iSpec)%nLines
        lamnu(NINT(SpeciesRadiation(iSpec)%Level(SpeciesRadiation(iSpec)%LinesInt(l,2),3))) &
          = lamnu(NINT(SpeciesRadiation(iSpec)%Level(SpeciesRadiation(iSpec)%LinesInt(l,2),3))) &
          + SpeciesRadiation(iSpec)%LinesReal(l,2)
      END DO

      low_IonizationPot = 2.9E-8*SQRT(NumDensElectrons/1.E6/MAX(1.,RadiationInput(iSpec)%Telec))*ElementaryCharge

  ! --- total collisional frequency for phase broadening
      coll_freq_ij=0.0
      DO jSpec = 1, nSpecies
        sigma_ij = Pi*(RadiationInput(iSpec)%Radius + RadiationInput(jSpec)%Radius)**2
        T_mean = (Species(iSpec)%MassIC*RadiationInput(iSpec)%NumDens*RadiationInput(iSpec)%Ttrans(4)+ &
             Species(jSpec)%MassIC*RadiationInput(jSpec)%NumDens*RadiationInput(jSpec)%Ttrans(4)) / &
            (Species(iSpec)%MassIC*RadiationInput(iSpec)%NumDens+Species(jSpec)%MassIC*RadiationInput(jSpec)%NumDens)
        coll_freq_ij = coll_freq_ij + 2.0 * RadiationInput(iSpec)%NumDens * RadiationInput(jSpec)%NumDens * sigma_ij &
          * SQRT(2*BoltzmannConst * T_mean * (Species(iSpec)%MassIC + Species(jSpec)%MassIC) &
          / (Pi*Species(iSpec)%MassIC*Species(jSpec)%MassIC))
      END DO
      
      nLines_considered = 0

  ! --- loop over all transition lines
      DO k=1,SpeciesRadiation(iSpec)%nLines

        IF (SpeciesRadiation(iSpec)%Level(NINT(SpeciesRadiation(iSpec)%Level(SpeciesRadiation(iSpec)%LinesInt(k,2),3)),2) &
          .LE. (RadiationInput(iSpec)%IonizationEn - low_IonizationPot)) THEN

  ! --- local emission and absorption coefficients
          etot = c_emi * REAL(SpeciesRadiation(iSpec)%LinesInt(k,4)) &
            / SpeciesRadiation(iSpec)%Level(NINT(SpeciesRadiation(iSpec)%Level(SpeciesRadiation(iSpec)%LinesInt(k,2),3)),1) &
            * SpeciesRadiation(iSpec)%NumDensExc(NINT(SpeciesRadiation(iSpec)%Level(SpeciesRadiation(iSpec)%LinesInt(k,2),3))) &
            * SpeciesRadiation(iSpec)%LinesReal(k,2) / (SpeciesRadiation(iSpec)%LinesReal(k,1))
          abstot = c_abs * (SpeciesRadiation(iSpec)%LinesReal(k,1))**4 * SpeciesRadiation(iSpec)%LinesReal(k,2) &
            * (REAL(SpeciesRadiation(iSpec)%LinesInt(k,4)) &
            / SpeciesRadiation(iSpec)%Level(NINT(SpeciesRadiation(iSpec)%Level(SpeciesRadiation(iSpec)%LinesInt(k,1),3)),1) &
            * SpeciesRadiation(iSpec)%NumDensExc(NINT(SpeciesRadiation(iSpec)%Level(SpeciesRadiation(iSpec)%LinesInt(k,1),3))) &
            - REAL(SpeciesRadiation(iSpec)%LinesInt(k,4)) &
            / SpeciesRadiation(iSpec)%Level(NINT(SpeciesRadiation(iSpec)%Level(SpeciesRadiation(iSpec)%LinesInt(k,2),3)),1) &
            * SpeciesRadiation(iSpec)%NumDensExc(NINT(SpeciesRadiation(iSpec)%Level(SpeciesRadiation(iSpec)%LinesInt(k,2),3))) )

  ! --- broadening mechanisms
          cwav = ((SpeciesRadiation(iSpec)%LinesReal(k,1))**2) / c

  ! --- Stark broadening
          Dlams = 2.*RadiationInput(iSpec)%NumDens*SpeciesRadiation(iSpec)%LinesReal(k,3) &
            * (RadiationInput(iSpec)%Telec*1.E-4)**RadiationInput(iSpec)%Starkex*1.E-32

  ! --- Van der Waals broadening
          Dlamvw = 20.0 * cwav * 4.5214D-18 * ntot * (3.* BoltzmannConst* RadiationInput(iSpec)%Ttrans(4)*ntot/rho)**.3

  ! --- natural broadening
          Dlamn = cwav * (lamnu(NINT(SpeciesRadiation(iSpec)%Level(SpeciesRadiation(iSpec)%LinesInt(k,2),3))) &
                  + lamnu(NINT(SpeciesRadiation(iSpec)%Level(SpeciesRadiation(iSpec)%LinesInt(k,1),3))))

  ! --- Lorentz (phase) broadening
          Dlamp = 2.0 * cwav * coll_freq_ij / RadiationInput(iSpec)%NumDens

  ! --- Resonance broadening (NEQAIR)
          Dlamr = 1.03*1.E-11 * SQRT( REAL(SpeciesRadiation(iSpec)%LinesInt(k,4))/REAL(SpeciesRadiation(iSpec)%LinesInt(k,3)) ) &
            * (SpeciesRadiation(iSpec)%LinesReal(k,1))**5 * SpeciesRadiation(iSpec)%LinesReal(k,2) &
            * SpeciesRadiation(iSpec)%NumDensExc(NINT(SpeciesRadiation(iSpec)%Level(SpeciesRadiation(iSpec)%LinesInt(k,1),3)))

  ! --- Total Lorentz width
          Dlaml = Dlams + Dlamvw + Dlamn + Dlamp + Dlamr

  ! --- Doppler broadening
          Dlamd = c_dopp*SpeciesRadiation(iSpec)%LinesReal(k,1)*SQRT(RadiationInput(iSpec)%Ttrans(4)/Species(iSpec)%MassIC )

  ! --- determine actual wavelength range to be calculated
          wrange = 100. * MAX(Dlaml,Dlamd)
          wleft  = MAX(RadiationParameter%WaveLen(1), SpeciesRadiation(iSpec)%LinesReal(k,1) - wrange )
          wright = MIN(RadiationParameter%WaveLen(RadiationParameter%WaveLenDiscr), SpeciesRadiation(iSpec)%LinesReal(k,1) + wrange)

  ! --- check if line is within range
          IF ( wright .GE. RadiationParameter%WaveLen(1) &
            .AND. wleft .LE. RadiationParameter%WaveLen(RadiationParameter%WaveLenDiscr)) THEN

            nLines_considered = nLines_considered + 1

  ! --- parameters for Voigt line profiles
            dpar  = (Dlaml-Dlamd) / (Dlaml+Dlamd)
            beta  = 0.023665 * EXP(0.6*dpar) + 0.00418 * EXP(-1.9*dpar)
            rpar  = 1. - 0.18121 * (1.-dpar**2) - beta * SIN(Pi*dpar)
            Dlamv = rpar*(Dlaml+Dlamd)

  ! --- generate Voigt line profile
            Voigt_int(-Voigt_int_nmax) = 0.                             !initialize Voigt profile
            Voigt_int( Voigt_int_nmax) = 1.

            Voigt_int_distmax = Voigt_mult_discr * Dlamv                !determine discretised interval
            Voigt_arg_dimless = Dlaml/Dlamv
            Voigt_var1        = 1./((1.065 + (0.447 + 0.058 * Voigt_arg_dimless) * Voigt_arg_dimless) * Dlamv)
            Voigt_const2      = Voigt_arg_dimless * Voigt_var1
            Voigt_const1      = Voigt_var1 - Voigt_const2
            Voigt_int_step    = Voigt_int_distmax / (Voigt_int_nmax-1) !step width

            DO iVoigt = -(Voigt_int_nmax-1), 0                          !determine lower half of Voigt profile
              Voigt_dist1   = REAL(ABS(iVoigt)) * Voigt_int_step / (Dlamv)
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

            DO iVoigt = -(Voigt_int_nmax-1), (Voigt_int_nmax-1)          ! normalize profiles that energy below the Voigt profile function is 1 (integrated profiles are normed to 1)
              Voigt_int(iVoigt) = Voigt_int(iVoigt) * Voigt_int_NormEnergy
            END DO

  !          IF(eps .GE. 1.0E-25) THEN
            CALL Radiation_Atomic_Transition_Line_Profile(Radiation_Profile, SpeciesRadiation(iSpec)%LinesReal(k,1), Voigt_int, Voigt_int_distmax, &
                      epsilon_at, epsilon_iSpec, etot, Radiation_Absorption_spec, abs_iSpec, abstot, iElem)
  !          END IF

          END IF

        END IF

      END DO

!    hilf = LEN_TRIM(RadiationInput(iSpec)%RadiationSpectraFileName)
!    RadiationInput(iSpec)%RadiationSpectraFileName = RadiationInput(iSpec)%RadiationSpectraFileName(1:hilf-4)
!    WRITE(*,*) 'calculated ',nLines_considered,' bound-bound lines of ', TRIM(RadiationInput(iSpec)%RadiationSpectraFileName)

    END IF

    TempOut_Em = 0.0
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
    
    DEALLOCATE(lamnu)

  END DO

  DO iWave = 1, RadiationParameter%WaveLenDiscr
    iWaveCoarse = INT((iWave-1)/RadiationParameter%WaveLenReductionFactor) + 1
    IF (iWaveCoarse.GT.RadiationParameter%WaveLenDiscrCoarse) iWaveCoarse = RadiationParameter%WaveLenDiscrCoarse
    Radiation_Emission_spec(iWaveCoarse,iElem) = Radiation_Emission_spec(iWaveCoarse,iElem) + epsilon_at(iWave)/RadiationParameter%WaveLenReductionFactor
  END DO

! --- add contribution to total emission
  em_atom = epsilon_at(1) * RadiationParameter%WaveLenIncr
  DO i=2, RadiationParameter%WaveLenDiscr
    em_atom = em_atom + epsilon_at(i) * RadiationParameter%WaveLenIncr
  END DO

    ! WRITE(*,*) '***  ATOMIC BOUND-BOUND RADIATION SUCCESSFULLY DONE  ***'
    ! WRITE(*,*) ''

END SUBROUTINE radiation_atoms




SUBROUTINE Radiation_Atomic_Transition_Line_Profile(Radiation_Profile, wavelength, Voigt_int, Voigt_int_distmax, &
    epsilon_mol, epsilon_iSpec, eps, Radiation_Absorption_spec, abs_iSpec, abstot, iElem)
!===================================================================================================================================
! calculates emission profile functions and adds radiative energy to emission array
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_Radiation_Vars,    ONLY   : RadiationParameter
  USE MOD_Mesh_Tools,        ONLY : GetGlobalElemID
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT(INOUT)             :: Radiation_Profile(:), epsilon_mol(:), epsilon_iSpec(:), Radiation_Absorption_spec(:,:), &
                                     abs_iSpec(:)
  REAL, INTENT(IN)                :: wavelength, Voigt_int_distmax, eps, abstot, Voigt_int(-1001:1001)
  INTEGER, INTENT(IN)             :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER           :: startwavelength_int, endwavelength_int, i, iWaveCoarse, iWave
  LOGICAL           :: add_radline
!===================================================================================================================================
  ! --- determination of indices for first and last entry of Voigt-profiles on wavelength axis
  startwavelength_int = INT(MAX(1.0 &
    , (((wavelength - Voigt_int_distmax - RadiationParameter%MinWaveLen) / RadiationParameter%WaveLenIncr) + 1.0)))
  endwavelength_int   = INT(MIN(REAL(RadiationParameter%WaveLenDiscr) &
    , (((wavelength + Voigt_int_distmax - RadiationParameter%MinWaveLen) / RadiationParameter%WaveLenIncr) + 1.0)))

!  startwavelength_int = MAX(1, INT((wavelength - Voigt_int_distmax &
!    - RadiationParameter%MinWaveLen) / RadiationParameter%WaveLenIncr) + 1)
!  endwavelength_int   = MIN(RadiationParameter%WaveLenDiscr, INT((wavelength &
!    + Voigt_int_distmax - RadiationParameter%MinWaveLen) / RadiationParameter%WaveLenIncr) + 1)

  add_radline          = .FALSE.

  IF ( (startwavelength_int .LT. endwavelength_int) .AND. &
                      ((wavelength+Voigt_int_distmax) .GT. RadiationParameter%MinWaveLen) .AND. &
                      ((wavelength-Voigt_int_distmax) .LT. RadiationParameter%MaxWaveLen) ) THEN

    add_radline = .TRUE.

  ! --- determine transition lines of previous computed Voigt-profiles
      CALL Radiation_Voigt_wavelength_interpolation(Voigt_int, Voigt_int_distmax, wavelength, &
        Radiation_Profile, startwavelength_int, endwavelength_int)

    IF (add_radline) THEN

      ! --- add radiative energy to emission
      DO i=0  , (endwavelength_int-startwavelength_int + 0)
        epsilon_iSpec(startwavelength_int+i) &
          = epsilon_iSpec(startwavelength_int+i) + MAX(0.0,eps) * Radiation_Profile(startwavelength_int+i)
        epsilon_mol(startwavelength_int+i) &
          = epsilon_mol(startwavelength_int+i) + MAX(0.0,eps) * Radiation_Profile(startwavelength_int+i)
        abs_iSpec(startwavelength_int+i) &
          = abs_iSpec(startwavelength_int+i)+MAX(0.0,abstot)*Radiation_Profile(startwavelength_int+i)
        iWave = startwavelength_int+i
        iWaveCoarse = INT((iWave-1)/RadiationParameter%WaveLenReductionFactor) + 1
        IF (iWaveCoarse.GT.RadiationParameter%WaveLenDiscrCoarse) iWaveCoarse = RadiationParameter%WaveLenDiscrCoarse
        Radiation_Absorption_spec(iWaveCoarse,GetGlobalElemID(iElem)) &
          = Radiation_Absorption_spec(iWaveCoarse,GetGlobalElemID(iElem))+MAX(0.0,abstot)*Radiation_Profile(iWave)/RadiationParameter%WaveLenReductionFactor 
      END DO

    END IF

  END IF

END SUBROUTINE Radiation_Atomic_Transition_Line_Profile




SUBROUTINE Radiation_Voigt_wavelength_interpolation(Voigt_int, Voigt_int_distmax, centerwavelength, &
    Radiation_Profile, startwavelength_int, endwavelength_int)
!===================================================================================================================================
! distributes transition lines with precomputed integrated Voigt-profiles (Voigt_int)
!===================================================================================================================================
! MODULES
  USE MOD_Globals
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

  Voigt_int_value_old = 0.0!startvalue !TODO: check!
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

  DO i=1, (endwavelength_int-startwavelength_int)
    Radiation_Profile(startwavelength_int+i) = Radiation_Profile(startwavelength_int+i) / RadiationParameter%WaveLenIncr
  END DO

!  ALLOCATE(help_Radiation_Profile(RadiationParameter%WaveLenDiscr)) !2nd order integration, but not better for tested lines
!  help_Radiation_Profile = Radiation_Profile
!  Radiation_Profile = 0.0
!  DO i=0, width_int+1
!    IF(i .EQ. 0) THEN
!      Radiation_Profile(startwavelength_int+i) = (help_Radiation_Profile(startwavelength_int+i) &
!        + help_Radiation_Profile(startwavelength_int+i+2)) / (2.*RadiationParameter%WaveLenIncr)
!    ELSE IF(i .EQ. (width_int+1)) THEN
!      Radiation_Profile(startwavelength_int+i) = (help_Radiation_Profile(startwavelength_int+i) &
!        + help_Radiation_Profile(startwavelength_int+i-2)) / (2.*RadiationParameter%WaveLenIncr)
!    ELSE
!      Radiation_Profile(startwavelength_int+i) = (help_Radiation_Profile(startwavelength_int+i-1) &
!        + help_Radiation_Profile(startwavelength_int+i+1)) / (2.*RadiationParameter%WaveLenIncr)
!    END IF
!  END DO


END SUBROUTINE Radiation_Voigt_wavelength_interpolation


END MODULE MOD_Radiation_Atoms
