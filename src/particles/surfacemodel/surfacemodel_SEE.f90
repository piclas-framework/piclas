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

SUBROUTINE SecondaryElectronEmission(PartID_IN,locBCID,ProductSpec,ProductSpecNbr,TempErgy)
!----------------------------------------------------------------------------------------------------------------------------------!
! Determine the probability of an electron being emitted due to an impacting particles (ion/electron bombardment)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals                   ,ONLY: abort,PARTISELECTRON,DOTPRODUCT
USE MOD_Globals_Vars              ,ONLY: c,c2,Joule2eV
USE MOD_Particle_Vars             ,ONLY: PartState,Species,PartSpecies,PartMPF,nSpecies
USE MOD_Globals_Vars              ,ONLY: ElementaryCharge,ElectronMass
USE MOD_SurfaceModel_Vars         ,ONLY: BulkElectronTempSEE
USE MOD_SurfaceModel_Vars         ,ONLY: SurfModResultSpec,SurfModEmissionYield,SurfModEmissionEnergy,SurfModEnergyDistribution
USE MOD_SurfaceModel_Vars         ,ONLY: SurfModSEEPowerFit
USE MOD_Particle_Boundary_Vars    ,ONLY: PartBound
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: CalcElectronSEE,SEE
USE MOD_Particle_Analyze_Tools    ,ONLY: CalcEkinPart2
USE MOD_PARTICLE_Vars             ,ONLY: usevMPF
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER,INTENT(IN)      :: PartID_IN           !< Bombarding Particle ID
INTEGER,INTENT(IN)      :: locBCID
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)     :: ProductSpec(2)      !< ProductSpec(1) new ID of impacting particle (the old one can change)
                                               !< ProductSpec(2) new ID of newly released electron
INTEGER,INTENT(OUT)     :: ProductSpecNbr      !< number of species for ProductSpec(1)
REAL,INTENT(OUT)        :: TempErgy            !< temperature, energy or velocity used for VeloFromDistribution
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: SpecID ! Species index of bombarding particle
REAL              :: eps_e  ! Energy of bombarding electron in eV
REAL              :: iRan   ! Random number
REAL              :: k_ee   ! Coefficient of emission of secondary electron
REAL              :: k_refl ! Coefficient for reflection of bombarding electron
REAL              :: W0,W1,W2
REAL              :: velo2
REAL              :: MPF
REAL              :: SEE_Prob
!===================================================================================================================================
SpecID = PartSpecies(PartID_IN)
! Squared velocity of bombarding particle
velo2=DOTPRODUCT(PartState(4:6,PartID_IN))
! Sanity check: is the impacting particle faster than c
IF(velo2.GT.c2) CALL abort(__STAMP__,'SecondaryElectronEmission: Bombarding particle is faster than the speed of light: ',RealInfoOpt=SQRT(velo2))
! Default 0
ProductSpec    = 0
ProductSpecNbr = 0
TempErgy       = 0.0
! Select particle surface modeling
SELECT CASE(PartBound%SurfaceModel(locBCID))
CASE(4) ! 4: SEE-E by power-law: a*T(eV)^b
  ProductSpecNbr = 0 ! do not create new particle (default value)
  ! Bombarding electron
  IF(PARTISELECTRON(PartID_IN))THEN
    ! Electron energy in [eV]
    eps_e = 0.5*Species(SpecID)%MassIC*velo2*Joule2eV ! Incident electron energy [eV]
    ! Power Fit
    SEE_Prob = SurfModSEEPowerFit(1,locBCID)*eps_e**SurfModSEEPowerFit(2,locBCID)
    ! If the yield is greater than 1.0 (or 2.0 or even higher), set the number of products with the integer and roll the dice for the remainder
    ProductSpecNbr = INT(SEE_Prob)
    SEE_Prob = SEE_Prob - REAL(ProductSpecNbr)

    ! Roll the dice and if the yield is greater than the random number, add an additional electron
    CALL RANDOM_NUMBER(iRan)
    IF(iRan.LT.SEE_Prob) ProductSpecNbr = ProductSpecNbr + 1

    ! If the electron is reflected (ProductSpecNbr=1) or multiple electrons are created (ProductSpecNbr>1)
    IF(ProductSpecNbr.GT.0) ProductSpec(2) = SurfModResultSpec(locBCID,SpecID)

    ! When more than 1 electron is created, give them all part of the impacting energy
    IF(ProductSpecNbr.GT.1) eps_e = eps_e/REAL(ProductSpecNbr) ! [eV]

    ! Velocity of reflected primary or secondary electrons in [m/s]
    TempErgy = SQRT(2.*eps_e*ElementaryCharge/ElectronMass)

  ELSE ! Neutral bombarding particle
    RETURN ! nothing to do
  END IF
CASE(5) ! 5: SEE by Levko2015 for copper electrodes
  !     !    by D. Levko, Breakdown of atmospheric pressure microgaps at high excitation, J. Appl. Phys. 117, 173303 (2015)

  ProductSpec(1)  = SpecID ! old particle

  ASSOCIATE (&
        phi            => 4.4  ,& ! eV -> cathode work function phi Ref. [20] Y. P. Raizer, Gas Discharge Physics (Springer, 1991)
        I              => 15.6  & ! eV -> ionization threshold of N2
        )
    IF(PARTISELECTRON(PartID_IN))THEN ! Bombarding electron
      ASSOCIATE (& ! Empirical fitting constants for a copper electrode Ref. [19] R. Cimino et al., Phys.Rev.Lett. 2004
            delta_star_max => 1.06 ,&
            s              => 1.35 ,&
            eps_max        => 262  ,& ! eV
            eps_0          => 150  ,& ! eV
            mass           => Species(SpecID)%MassIC &! mass of bombarding particle
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
            !ReflectionIndex = 3 ! SEE + perfect elastic scattering of the bombarding electron
            ProductSpec(2)  = SurfModResultSpec(locBCID,SpecID)  ! Species of the injected electron
            ProductSpecNbr = 1
            TempErgy        = SQRT(2.*(eps_e*ElementaryCharge-ElementaryCharge*phi)/ElectronMass) ! Velocity of emitted secondary electron
            eps_e           = 0.5*mass*(TempErgy**2)/ElementaryCharge               ! Energy of the injected electron
!WRITE (*,*) CHAR(27) // "[0;34mPartID_IN                  =", PartID_IN,CHAR(27),"[m"
!WRITE (*,*) CHAR(27) // "[0;34mPartState(1:3,PartID_IN)   =", PartState(1:3,PartID_IN),CHAR(27),"[m"
!WRITE (*,*) CHAR(27) // "[0;34mPartState(4:6,PartID_IN)   =", PartState(4:6,PartID_IN),CHAR(27),"[m"
!WRITE (*,*) CHAR(27) // "[0;34mBombarding electron: TempErgy =", TempErgy,CHAR(27),"[m"
!WRITE (*,*) CHAR(27) // "[0;34m                     eps_e =", eps_e,CHAR(27),"[m"
          ELSE ! Only perfect elastic scattering of the bombarding electron
            ProductSpecNbr = 0 ! do not create new particle
          END IF
          ! Original Code as described in the paper by Levko (2015)
          !   IF(k_refl.LT.k_ee)THEN ! -> region with a high chance of SEE
          !     IF(iRan.LT.k_refl/(k_ee+k_refl))THEN ! Reflection
          !       ReflectionIndex = -1 ! only perfect elastic scattering of the bombarding electron
          !     ELSE ! SEE-E
          !       ReflectionIndex = -2 ! SEE + perfect elastic scattering of the bombarding electron
          !       ProductSpec(2)      = 4
          !     END IF
          !   ELSE ! (k_ee.LE.k_refl) -> region with a high chance of reflection
          !     IF(iRan.LT.k_ee/(k_ee+k_refl))THEN ! SEE
          !       ReflectionIndex = -2 ! SEE + perfect elastic scattering of the bombarding electron
          !       ProductSpec(2)      = 4
          !     ELSE ! Reflection
          !       ReflectionIndex = -1 ! only perfect elastic scattering of the bombarding electron
          !     END IF
          !   END IF
        ELSE ! Removal of the bombarding electron
          !ReflectionIndex = 3 ! Removal of the bombarding electron
          ProductSpec(1) = 0 ! just for sanity check
        END IF
      END ASSOCIATE
    ELSEIF(Species(SpecID)%ChargeIC.GT.0.0)THEN ! Positive bombarding ion
      CALL RANDOM_NUMBER(iRan)
      !IF(iRan.LT.1.)THEN ! SEE-I: gamma=0.02 for the N2^+ ions and copper material
      IF(iRan.LT.0.02)THEN ! SEE-I: gamma=0.02 for the N2^+ ions and copper material
        !ReflectionIndex = -2       ! SEE + perfect elastic scattering of the bombarding electron
        ProductSpec(2) = SurfModResultSpec(locBCID,SpecID)  ! Species of the injected electron
        ProductSpecNbr = 1
        eps_e          = I-2.*phi ! Energy of the injected electron
        TempErgy       = SQRT(2.*(eps_e*ElementaryCharge-ElementaryCharge*phi)/ElectronMass) ! Velocity of emitted secondary electron
!WRITE (*,*) CHAR(27) // "[0;34mPartID_IN                  =", PartID_IN,CHAR(27),"[m"
!WRITE (*,*) CHAR(27) // "[0;31mPartState(1:3,PartID_IN) =", PartState(1:3,PartID_IN),CHAR(27),"[m"
!WRITE (*,*) CHAR(27) // "[0;31mPartState(4:6,PartID_IN) =", PartState(4:6,PartID_IN),CHAR(27),"[m"
!WRITE (*,*) CHAR(27) // "[0;31mPositive bombarding ion: TempErgy =", TempErgy,CHAR(27),"[m"
!WRITE (*,*) CHAR(27) // "[0;31m                         eps_e =", eps_e,CHAR(27),"[m"
      ELSE ! Removal of the bombarding ion
        !ReflectionIndex = -1 ! Only perfect elastic scattering of the bombarding electron
        ProductSpec(1)  = -SpecID ! Negative value: Remove bombarding particle and sample
        ProductSpecNbr = 0 ! do not create new particle
      END IF
    ELSE ! Neutral bombarding particle
    !  IF(iRan.LT.0.1)THEN ! SEE-N: from svn-trunk PICLas version
    !    !ReflectionIndex = -2 ! SEE + perfect elastic scattering of the bombarding electron
    !    ProductSpec(2)  = SurfModelResultSpec(locBCID,SpecID)  ! Species of the injected electron
    !    ProductSpecNbr = 1
    !  ELSE
    !    !ReflectionIndex = -1 ! Only perfect elastic scattering of the bombarding electron
    !    ReflectionIndex = -1 ! Removal of the bombarding neutral
        ProductSpecNbr = 0 ! do not create new particle
    !    WRITE (*,*) "Neutral bombarding particle =", PartID_IN
    !    stop "stop"
    !  END IF
    END IF
  END ASSOCIATE

  ! Sanity check: is the newly created particle faster than c
  IF(TempErgy.GT.c) CALL abort(__STAMP__,'SecondaryElectronEmission: Particle is faster than the speed of light:'&
      ,RealInfoOpt=TempErgy)

CASE(6) ! 6: SEE by Pagonakis2016 (originally from Harrower1956)
  CALL abort(__STAMP__,'Not implemented yet')
CASE(7) ! 7: SEE-I (bombarding electrons are removed, Ar+ on different materials is considered for SEE)
  ProductSpec(1)  = -SpecID ! Negative value: Remove bombarding particle and sample
  ProductSpecNbr = 0 ! do not create new particle (default value)

  IF(PARTISELECTRON(PartID_IN))THEN ! Bombarding electron
    RETURN ! nothing to do
  ELSEIF(Species(SpecID)%ChargeIC.GT.0.0)THEN ! Positive bombarding ion
    ! SEE-I bombarding e- are removed, Ar+ on different materials is considered for secondary e- emission (the default probability
    ! is 0.13 probability, see  D. Depla, Magnetron sputter deposition: Linking discharge voltage with target properties, 2009)

    ! If yield is greater than 1, store the leading integer here
    IF(SurfModEmissionYield(locBCID).GE.1.0) ProductSpecNbr = INT(SurfModEmissionYield(locBCID))
    CALL RANDOM_NUMBER(iRan)
    IF(iRan.LT.MOD(SurfModEmissionYield(locBCID), 1.0) ) ProductSpecNbr = ProductSpecNbr + 1 ! Create one additional new particle

    IF(ProductSpecNbr.GT.0)THEN
      ProductSpec(2) = SurfModResultSpec(locBCID,SpecID)  ! Species of the injected electron
      ! Set TempErgy (velocity in m/s or energy in eV might be required here)
      IF(SurfModEmissionEnergy(locBCID).GE.0.)THEN
        ! Electron energy in [eV] or [m/2] is required here depending on the chosen distribution function
        IF(SurfModEnergyDistribution(locBCID).EQ.'uniform-energy')THEN
          ! TempErgy is in [eV], which will be set to a uniform-energy distribution function
          TempErgy = SurfModEmissionEnergy(locBCID)
        ELSE
          ! TempErgy is in [m/s]
          TempErgy = SQRT(2.0*SurfModEmissionEnergy(locBCID)*ElementaryCharge/Species(ProductSpec(2))%MassIC)
        END IF ! SurfModEnergyDistribution(locBCID).EQ.
      ELSE
        ! Get velocity of new electron (from impacting ion energy)
        TempErgy = CalcEkinPart2(PartState(4:6,PartID_IN),SpecID,1.0) ! [J]
        TempErgy = SQRT(2.0 * Tempergy / ElectronMass) ! [m/s]
      END IF ! SurfModEmissionEnergy
    END IF ! ProductSpecNbr.GT.0
  ELSE ! Neutral bombarding particle
    RETURN ! nothing to do
  END IF

  ! Sanity check: is the newly created particle faster than c
  IF(TempErgy.GT.c) CALL abort(__STAMP__,'SecondaryElectronEmission: Particle is faster than the speed of light:'&
      ,RealInfoOpt=TempErgy)

CASE(8) ! 8: SEE-E (e- on dielectric materials is considered for SEE and three different outcomes)
        !    by A.I. Morozov, "Structure of Steady-State Debye Layers in a Low-Density Plasma near a Dielectric Surface", 2004

    IF(PARTISELECTRON(PartID_IN))THEN ! Bombarding electron
      ASSOCIATE( P0   => 0.9               ,& ! Assumption in paper
                 Te0  => BulkElectronTempSEE  ,& ! Assumed bulk electron temperature [eV] (note this parameter is read as [K])
                 mass => Species(SpecID)%MassIC  ) ! mass of bombarding particle
        eps_e = 0.5*mass*velo2*Joule2eV ! Incident electron energy [eV]
        ASSOCIATE( alpha0 => 1.5*Te0 ,& ! Energy normalization parameter
                   alpha2 => 6.0*Te0  ) ! Energy normalization parameter
          W0 = P0*EXP(-(eps_e/alpha0)**2)
          W2 = 1.0 - EXP(-(eps_e/alpha2)**2)
          W1 = 1.0 - W2 - W0
          CALL RANDOM_NUMBER(iRan) ! 1st random number
          IF(iRan.GT.W0)THEN ! Remove incident electron
            iRan = iRan - W0
            IF(iRan.LT.W1)THEN ! 1 SEE
              !ASSOCIATE( P10 => 1.5*W1/eps_e )
                ProductSpec(2) = SurfModResultSpec(locBCID,SpecID)  ! Species of the injected electron
                ProductSpecNbr = 1 ! Create one new particle
                ProductSpec(1) = SpecID ! Reflect old particle
                !const          = P10 ! Store constant here for usage in VeloFromDistribution()
              !END ASSOCIATE
            ELSE ! 2 SEE
              !ASSOCIATE( P20 => 3.0*W2/(eps_e**2) )
                ProductSpec(2) = SurfModResultSpec(locBCID,SpecID)  ! Species of the injected electron
                ProductSpecNbr = 2 ! Create two new particles
                ProductSpec(1) = SpecID ! Reflect old particle
                !const          = P20 ! Store constant here for usage in VeloFromDistribution()
              !END ASSOCIATE
            END IF
            TempErgy = eps_e ! electron energy
          ELSE
            ProductSpec(1) = -SpecID ! Negative value: Remove bombarding particle and sample
          END IF
        END ASSOCIATE
      END ASSOCIATE
    END IF

CASE(9) ! 9: SEE-I when Ar^+ ion bombards surface with 0.01 probability and fixed SEE electron energy of 6.8 eV
  ProductSpec(1)  = -SpecID ! Negative value: Remove bombarding particle and sample
  IF(Species(SpecID)%ChargeIC.GT.0.0)THEN ! Bombarding positive ion
    CALL RANDOM_NUMBER(iRan) ! 1st random number
    ASSOCIATE( eps_e => 6.8 )! Ejected electron energy [eV]
      IF(iRan.LT.0.01)THEN ! SEE-I: gamma=0.01 for the bombarding Ar^+ ions
        ProductSpec(2) = SurfModResultSpec(locBCID,SpecID) ! Species of the injected electron
        ProductSpecNbr = 1 ! Create one new particle
        TempErgy       = SQRT(2.*eps_e*ElementaryCharge/ElectronMass) ! Velocity of emitted secondary electron in [m/s]
      END IF
    END ASSOCIATE
  END IF

CASE(10) ! 10: SEE-I (bombarding electrons are removed, Ar+ on copper is considered for SEE)
         !     by J.G. Theis "Computing the Paschen curve for argon with speed-limited particle-in-cell simulation", 2021
         !     Plasmas 28, 063513, doi: 10.1063/5.0051095
  ProductSpec(1)  = -SpecID ! Negative value: Remove bombarding particle and sample
  ProductSpecNbr = 0 ! do not create new particle (default value)

  IF(PARTISELECTRON(PartID_IN))THEN ! Bombarding electron
    RETURN ! nothing to do
  ELSEIF(Species(SpecID)%ChargeIC.GT.0.0)THEN ! Positive bombarding ion
    ASSOCIATE (mass  => Species(SpecID)%MassIC )! mass of bombarding particle
      ! Electron energy in [eV]
      eps_e = 0.5*mass*velo2/ElementaryCharge
      IF(eps_e.LT.700) THEN
        SEE_Prob = 0.09*(eps_e/700.0)**0.05
      ELSE
        SEE_Prob = 0.09*(eps_e/700.0)**0.72
      END IF
    END ASSOCIATE
    CALL RANDOM_NUMBER(iRan)
    IF(iRan.LT.SEE_Prob)THEN
      ProductSpec(2) = SurfModResultSpec(locBCID,SpecID)  ! Species of the injected electron
      ProductSpecNbr = 1 ! Create one new particle
      TempErgy       = 0.0 ! emit electrons with zero velocity
    END IF
  ELSE ! Neutral bombarding particle
    RETURN ! nothing to do
  END IF

CASE(11) ! 11: SEE-E by e- on quartz (SiO2) by A. Dunaevsky, "Secondary electron emission from dielectric materials of a Hall
         !     thruster with segmented electrodes", 2003
         !     PHYSICS OF PLASMAS, VOLUME 10, NUMBER 6, DOI: 10.1063/1.1568344
  ProductSpec(1)  = -SpecID ! Negative value: Remove bombarding particle and sample
  ProductSpecNbr = 0 ! do not create new particle (default value)

  IF(PARTISELECTRON(PartID_IN))THEN ! Bombarding electron
    ASSOCIATE (mass  => Species(SpecID)%MassIC )! mass of bombarding particle
      ! Electron energy in [eV]
      eps_e = 0.5*mass*velo2*Joule2eV ! Incident electron energy [eV]

      ! Linear Fit
      SEE_Prob = 0.8 + 0.2 * eps_e/35.0
      ! Power Fit
      !SEE_Prob = (eps_e/30.0)**0.26

    END ASSOCIATE
    ! If the yield is greater than 1.0 (or 2.0 or even higher) store the integer and roll the dice for the remainder
    ProductSpecNbr = INT(SEE_Prob)
    SEE_Prob = SEE_Prob - REAL(ProductSpecNbr)

    ! Roll the dice
    CALL RANDOM_NUMBER(iRan)
    IF(iRan.LT.SEE_Prob) ProductSpecNbr = ProductSpecNbr + 1 ! Create one additional electron

    ! If the electron is reflected (ProductSpecNbr=1) or multiple electrons are created (ProductSpecNbr>1)
    IF(ProductSpecNbr.GT.0) ProductSpec(2) = SurfModResultSpec(locBCID,SpecID)  ! Species of the injected electron

    ! When more than 1 electron is created, give them all part of the impacting energy, otherwise reflect the primary electron
    IF(ProductSpecNbr.GT.1) eps_e = eps_e/REAL(ProductSpecNbr) ! [eV]

    ! Velocity of reflected primary or secondary electrons in [m/s]
    TempErgy = SQRT(2.*eps_e*ElementaryCharge/ElectronMass)

  ELSEIF(Species(SpecID)%ChargeIC.GT.0.0)THEN ! Positive bombarding ion
    RETURN ! nothing to do
  ELSE ! Neutral bombarding particle
    RETURN ! nothing to do
  END IF

END SELECT

IF(ProductSpecNbr.GT.0)THEN
  ! Sanity check
  IF((ProductSpec(2).LE.0).OR.(ProductSpec(2).GT.nSpecies)) CALL abort(__STAMP__,&
      'SEE model trying to create particle with 0, negative or speciesID > nSpecies: ProductSpec(2)=',IntInfoOpt=ProductSpec(2))
  ! Check if SEE counter is active and assign the number of produced electrons to the boundary
  IF(CalcElectronSEE)THEN
    ASSOCIATE( iSEEBC => SEE%BCIDToSEEBCID(locBCID) )
      IF(iSEEBC.EQ.-1) CALL abort(__STAMP__,'SEE%BCIDToSEEBCID(locBCID)) = -1')
      IF(usevMPF)THEN
        ! MPF of impacting particle
        MPF = PartMPF(PartID_IN)
      ELSE
        ! MPF of produced species
        MPF = Species(ProductSpec(2))%MacroParticleFactor
      END IF
      ! Consider the number of produced electrons ProductSpecNbr and their charge q=Species(ProductSpec(2))%ChargeIC
      ! Note that the negative value of the charge -q is used below
      SEE%RealElectronOut(iSEEBC) = SEE%RealElectronOut(iSEEBC) - MPF*ProductSpecNbr*Species(ProductSpec(2))%ChargeIC
    END ASSOCIATE
  END IF ! CalcElectronSEE
END IF ! ProductSpecNbr.GT.0

END SUBROUTINE SecondaryElectronEmission


END MODULE MOD_SEE
