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

MODULE MOD_Equation
!===================================================================================================================================
! Add comments please!
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
INTERFACE InitEquation
  MODULE PROCEDURE InitEquation
END INTERFACE
INTERFACE ExactFunc
  MODULE PROCEDURE ExactFunc
END INTERFACE
INTERFACE CalcSource
  MODULE PROCEDURE CalcSource
END INTERFACE
INTERFACE DivCleaningDamping
  MODULE PROCEDURE DivCleaningDamping
END INTERFACE
PUBLIC::InitEquation,ExactFunc,CalcSource,FinalizeEquation,DivCleaningDamping
!===================================================================================================================================

PUBLIC::DefineParametersEquation
CONTAINS

!==================================================================================================================================
!> Define parameters for equation
!==================================================================================================================================
SUBROUTINE DefineParametersEquation()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Equation")

CALL prms%CreateRealOption(     'c_corr'           , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Multiplied with c0 results in the velocity of '//&
                                                     'introduced artificial correcting waves (HDC)' , '1.')
CALL prms%CreateRealOption(     'c0'               , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Velocity of light (in vacuum)' , '1.')
CALL prms%CreateRealOption(     'eps'              , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Electric constant (vacuum permittivity)' , '1.')
CALL prms%CreateRealOption(     'mu'               , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Magnetic constant (vacuum permeability = 4πE−7H/m)' &
                                                   , '1.')
CALL prms%CreateRealOption(     'fDamping'         , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Apply the damping factor also to PML source terms\n'//&
                                                     'but only to PML variables for Phi_E and Phi_B to prevent charge-related\n'//&
                                                     'instabilities (accumulation of divergence compensation over \n'//&
                                                     'timeU2 = U2 * fDamping' , '0.999')
CALL prms%CreateLogicalOption(  'ParabolicDamping' , 'TODO-DEFINE-PARAMETER' , '.FALSE.')
CALL prms%CreateLogicalOption(  'CentralFlux'      , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Flag for central or upwind flux' , '.FALSE.')
CALL prms%CreateIntOption(      'IniExactFunc'     , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Define exact function necessary for '//&
                                                     'linear scalar advection', '-1')

CALL prms%CreateLogicalOption(  'DoExactFlux'      , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Switch emission to flux superposition at'//&
                                                     ' certain positions' , '.FALSE.')
CALL prms%CreateRealArrayOption('xDipole'          , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Base point of electromagnetic dipole', '0. , 0. , 0.')
CALL prms%CreateRealOption(     'omega'            , 'TODO-DEFINE-PARAMETER\n'//&
                                                     '2*pi*f (f=100 MHz default)' , '6.28318e8')
CALL prms%CreateRealOption(     'tPulse'           , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Half length of pulse' , '30e-9')

CALL prms%CreateRealOption(     'TEFrequency'      , 'Frequency of TE wave')
CALL prms%CreateRealOption(     'TEScale'          , 'Scaling of input TE-wave strength' , '1.')
CALL prms%CreateStringOption(  'TEPolarization'   ,  'Polarization of the TE-mode: x = linear in x-direction\n'//&
                                                     '                             y = linear in y direction\n'//&
                                                     '                             l = left-handed circular\n'//&
                                                     '                             r = right-handed circular', 'x')
CALL prms%CreateLogicalOption(  'TEPulse'          , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Flag for pulsed or continuous wave' , '.FALSE.')
CALL prms%CreateIntArrayOption( 'TEMode'           , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Input of TE_n,m mode', '1 , 1')
CALL prms%CreateRealOption(     'TERadius'         , 'Radius of Input TE wave, if wave is inserted over a plane into a waveguide')
CALL prms%CreateRealOption(     'TEDelay'          , 'TODO-DEFINE-PARAMETER' , '0.5e-9')
CALL prms%CreateRealOption(     'TEPulseSigma'     , 'standard deviation of the Gaussian pulse' , '0.1e-9')
CALL prms%CreateRealOption(     'TEPulseSeriesFrequence', 'if TEPulseSeriesFrequence>0 and TEPulse=T\n'//&
                                                          ' -> a series of gaussian pulses with frequence TEPulseSeriesFrequence' , '0.')
CALL prms%CreateIntOption(      'TEPulseNumber'    , 'number of generated pulses in a pulse-series' , '1')
CALL prms%CreateRealOption(     'TEDirection'      , '+1 for propagation in +z direction, -1 for propagation in -z direction' , '1.0')
CALL prms%CreateStringOption(   'TEPulseShape'     , 'shape of the pulse:\n'//&
                                                     '  gaussian\n'//&
                                                     '  recangular\n'//&
                                                     '  rectangularGaussianEdges', 'gaussian')

CALL prms%CreateRealOption(     'WaveLength'       , 'Wavelength of the electromagnetic wave.' , '1.')
CALL prms%CreateRealArrayOption('WaveVector'       , 'Direction of traveling wave.', '0. , 0. , 1.')
CALL prms%CreateLogicalOption(  'UseWaveVectorE0dir','User defined E_0_vec unit vector in combination with WaveVector' , '.FALSE.')
CALL prms%CreateRealArrayOption('WaveVectorE0dir'  , 'Direction vector for the electric field\n'//&
                                                     ' (is not used on default, only when UseWaveVectorE0dir=T)', '0. , 0. , 1.')
CALL prms%CreateRealArrayOption('WaveBasePoint'    , 'Vector to the position of the beam origin', '0.5 , 0.5 , 0.')
CALL prms%CreateRealOption(     'I_0'              , 'Maximum intensity of the beam' , '1.')
CALL prms%CreateRealOption(     'sigma_t'          , 'First of two definitions for the pulse duration:\n'//&
                                                     'Can be used instead of tFWHM (time For Full '//&
                                                     'Wave Half Maximum)\ntFWHM = 2*sqrt(2*ln(2))*sigma_t' , '0.')
CALL prms%CreateRealOption(     'tFWHM'            , 'Second of two definitions for the pulse duration, '//&
                                                     'Time For Full Wave Half Maximum '//&
                                                     '(pulse duration within which the intensity'//&
                                                     'amplitude is higher than 50% of its maximum)' , '0.')
CALL prms%CreateRealOption(     'BeamEnergy'       , 'Total beam energy [J]' , '-1.0')
CALL prms%CreateRealOption(     'Beam_a0'          , 'Dimensionless beam amplitude \n'//&
                                                     '(value for scaling the max. electric field)' , '-1.0')
CALL prms%CreateRealOption(     'Beam_w0'          , 'Beam spot size (waist radius, where the beam radius has a minimum) \n'//&
                                                     '; the old variable name is "omega_0"' , '1.0')
CALL prms%CreateRealOption(     'Beam_t0'          , 'starting time of the (pulsed) electromagnetic wave. ' , '0.0')
CALL prms%CreateRealOption(     'omega_0'          , 'old variable for "Beam_w0"; remove this variable in 2019', '1.0')
CALL prms%CreateStringOption(   'BCStateFile'      , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Boundary Condition State File', 'no file found')
CALL prms%CreateIntOption(      'AlphaShape'       , 'TODO-DEFINE-PARAMETER', '2')
CALL prms%CreateRealOption(     'r_cutoff'         , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Modified for curved and shape-function influence'//&
                                                     ' (c*dt*SafetyFactor+r_cutoff)' , '1.0')

CALL prms%CreateIntOption(      'FluxDir'          , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Flux direction', '-1')
CALL prms%CreateIntOption(      'ExactFluxDir'     , 'TODO-DEFINE-PARAMETER\n'//&
                                                     'Flux direction for ExactFlux', '3')
CALL prms%CreateRealOption(     'ExactFluxPosition', 'TODO-DEFINE-PARAMETER\n'//&
                                                     'x,y, or z-position of interface')

END SUBROUTINE DefineParametersEquation

SUBROUTINE InitEquation()
!===================================================================================================================================
! Get the constant advection velocity vector from the ini file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: PI,ElectronMass,ElementaryCharge,c,c_inv,c2,mu0,eps0
USE MOD_ReadInTools
#ifdef PARTICLES
USE MOD_Interpolation_Vars ,ONLY: InterpolationInitIsDone
#endif
USE MOD_Equation_Vars
USE MOD_Mesh_Vars          ,ONLY: BoundaryType,nBCs,BC
USE MOD_Globals_Vars       ,ONLY: EpsMach
USE MOD_Mesh_Vars          ,ONLY: xyzMinMax,nSides,nBCSides
USE MOD_Mesh               ,ONLY: GetMeshMinMaxBoundaries
USE MOD_Utils              ,ONLY: RootsOfBesselFunctions
USE MOD_ReadInTools        ,ONLY: PrintOption
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: nRefStates,iBC,ntmp,iRefState
INTEGER,ALLOCATABLE              :: RefStates(:)
LOGICAL                          :: isNew
REAL                             :: PulseCenter
REAL,ALLOCATABLE                 :: nRoots(:)
LOGICAL                          :: DoSide(1:nSides)
INTEGER                          :: locType,locState,iSide
REAL                             :: BeamEnergy_loc,BeamFluency_loc,BeamArea_loc
!===================================================================================================================================
IF(EquationInitIsDone)THEN
#ifdef PARTICLES
  IF(InterpolationInitIsDone)THEN
    LBWRITE(*,*) "InitMaxwell not ready to be called or already called."
    RETURN
  END IF
#else
  LBWRITE(*,*) "InitMaxwell not ready to be called or already called."
  RETURN
#endif /*PARTICLES*/
END IF

LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT MAXWELL ...'

! Read correction velocity
c_corr             = GETREAL('c_corr','1.')
fDamping           = GETREAL('fDamping','0.999')
DoParabolicDamping = GETLOGICAL('ParabolicDamping','.FALSE.')
CentralFlux        = GETLOGICAL('CentralFlux','.FALSE.')
Beam_t0            = GETREAL('Beam_t0','0.0')

c_corr_c  = c_corr*c
c_corr_c2 = c_corr*c2
eta_c     = (c_corr-1.)*c

! default
WaveLength = -1.

! Read in boundary parameters
IniExactFunc = GETINT('IniExactFunc')
nRefStates=nBCs+1
nTmp=0
ALLOCATE(RefStates(nRefStates))
RefStates=0
IF(IniExactFunc.GT.0) THEN
  RefStates(1)=IniExactFunc
  nTmp=1
END IF
DO iBC=1,nBCs
  IF(BoundaryType(iBC,BC_STATE).GT.0)THEN
    isNew=.TRUE.
    ! check if boundarytype already exists
    DO iRefState=1,nTmp
      IF(BoundaryType(iBC,BC_STATE).EQ.RefStates(iRefState)) isNew=.FALSE.
    END DO
    IF(isNew)THEN
      nTmp=nTmp+1
      RefStates(ntmp)=BoundaryType(iBC,BC_STATE)
    END IF
  END IF
END DO
IF(nTmp.GT.0) DoExactFlux = GETLOGICAL('DoExactFlux','.FALSE.')
IF(DoExactFlux) CALL InitExactFlux()
DO iRefState=1,nTmp
  SELECT CASE(RefStates(iRefState))
  CASE(2,22)
    TEFrequency        = GETREAL('TEFrequency')
    TERadius           = GETREAL('TERadius')
    TEScale            = GETREAL('TEScale')
  CASE(4,40,41)
    xDipole(1:3)       = GETREALARRAY('xDipole',3,'0.,0.,0.') ! dipole base point
    DipoleOmega        = GETREAL('omega','6.28318E08')        ! f=100 MHz default
    tPulse             = GETREAL('tPulse','30e-9')            ! half length of pulse
  CASE(5,6)
    TEFrequency        = GETREAL('TEFrequency','35e9')
    TEScale            = GETREAL('TEScale')
    TEPolarization     = TRIM(GETSTR('TEPolarization','x'))
    TEDelay            = GETREAL('TEDelay','0.5e-9')
    TEPulse            = GETLOGICAL('TEPulse','.FALSE.')
    TEPulseSigma       = GETREAL('TEPulseSigma','0.1e-9')
    TEPulseSeriesFrequence = GETREAL('TEPulseSeriesFrequence','0.')
    TEPulseNumber      = GETINT('TEPulseNumber', '1')
    TEMode             = GETINTARRAY('TEMode',2,'1,1')
    TEDirection        = GETREAL('TEDirection','1.0')
    TEPulseShape       = TRIM(GETSTR('TEPulseShape','gaussian'))

    !check the TE-Input-Parameters
    IF (TRIM(TEPolarization).NE.'x'.AND.TRIM(TEPolarization).NE.'y'.AND.TRIM(TEPolarization).NE.'r'.AND.TRIM(TEPolarization).NE.'l') THEN
      SWRITE(UNIT_stdOut,*) 'Wrong value for TEPolarization=', TEPolarization
      CALL abort(__STAMP__,'Wrong value for TEPolarization')
    END IF
    IF ( (2*PI*TEFrequency*c_inv).LT.(TEModeRoot/TERadius)) THEN
      SWRITE(UNIT_stdOut,'(A,E25.14E3)')'(k)**2 = ',((2*PI*TEFrequency*c_inv))**2
      SWRITE(UNIT_stdOut,'(A,E25.14E3)')'kCut**2          = ',TEModeRoot/TERadius**2
      SWRITE(UNIT_stdOut,'(A)')'  Maybe frequency too small?'
      CALL abort(__STAMP__,'kz=SQRT(k**2-kCut**2), but the argument is negative!')
    END IF

    !ExactFluxPosition    = GETREAL('ExactFluxPosition', '0.0') !
    ! compute required roots
    ALLOCATE(nRoots(1:TEMode(2)))
    CALL RootsOfBesselFunctions(TEMode(1),TEMode(2),0,nRoots)
    TEModeRoot=nRoots(TEMode(2))
    DEALLOCATE(nRoots)
    ! check if it is a BC condition
    DO iBC=1,nBCs
      IF(BoundaryType(iBC,BC_STATE).EQ.5)THEN
        DoSide=.FALSE.
        DO iSide=1,nBCSides
          locType =BoundaryType(BC(iSide),BC_TYPE)
          locState=BoundaryType(BC(iSide),BC_STATE)
          IF(locState.EQ.5)THEN
            DoSide(iSide)=.TRUE.
          END IF ! locState.EQ.BCIn
        END DO
        ! call function to get radius
        CALL GetWaveGuideRadius(DoSide)
      END IF
    END DO
    IF((TEPolarization.NE.'x').AND.(TEPolarization.NE.'y').AND.(TEPolarization.NE.'r').AND.(TEPolarization.NE.'l'))THEN
      CALL abort(&
    __STAMP__&
    ,' TEPolarization has to be x,y,l or r.')
    END IF
    IF(TERadius.LT.0.0)THEN ! not set
      TERadius=GETREAL('TERadius','0.0')
      LBWRITE(UNIT_StdOut,*) ' TERadius not determined automatically. Set waveguide radius to ', TERadius
    END IF

    ! display cut-off freequncy for this mode
    LBWRITE(UNIT_stdOut,'(A,I5,A1,I5,A,E25.14E3,A)')&
           '  Cut-off frequency in circular waveguide for TE_[',1,',',1,'] is ',1.8412*c/(2*PI*TERadius),' Hz (lowest mode)'
    LBWRITE(UNIT_stdOut,'(A,I5,A1,I5,A,E25.14E3,A)')&
           '  Cut-off frequency in circular waveguide for TE_[',TEMode(1),',',TEMode(2),'] is ',(TEModeRoot/TERadius)*c/(2*PI),&
           ' Hz (chosen mode)'
  CASE(12,121,14,15,16)
    ! planar wave input: get wavelength and set wave number angular frequency
    WaveLength     = GETREAL('WaveLength','1.')
    BeamWaveNumber = 2.*PI/WaveLength
    BeamOmega      = BeamWaveNumber*c
    CALL PrintOption('BeamOmega [Hz]','CALCUL.',RealOpt=BeamOmega)
    CALL PrintOption('BeamPeriod [s]','CALCUL.',RealOpt=2.*PI/BeamOmega)

    ! set wave vector: direction of traveling wave
    WaveVector(1:3)= GETREALARRAY('WaveVector',3,'0.,0.,1.')
    WaveVector=UNITVECTOR(WaveVector)

    ! construct perpendicular electric field
    IF(ABS(WaveVector(3)).LT.EpsMach)THEN
      E_0_vec=(/ -WaveVector(2)-WaveVector(3)  , WaveVector(1) ,WaveVector(1) /)
    ELSE
      IF(ALMOSTEQUAL(ABS(WaveVector(3)),1.))THEN ! wave vector in z-dir -> E_0_vec in x-dir!
        E_0_vec= (/1., 0., 0./)
      ELSE
        E_0_vec=(/ WaveVector(3) , WaveVector(3) , -WaveVector(1)-WaveVector(2) /)
      END IF
    END IF
    ! normalize E-field
    E_0_vec=UNITVECTOR(E_0_vec)

    ! Use pre-defined E_0_vec vector
    UseWaveVectorE0dir = GETLOGICAL('UseWaveVectorE0dir','.FALSE.')
    IF(UseWaveVectorE0dir)THEN
      WaveVectorE0dir(1:3) = GETREALARRAY('WaveVectorE0dir',3,'1.0,0.0,0.0')
      WaveVectorE0dir      = UNITVECTOR(WaveVectorE0dir)
      IF(.NOT.ALMOSTZERO(DOT_PRODUCT(WaveVector,WaveVectorE0dir))) CALL abort(&
          __STAMP__&
          ,' WaveVector and WaveVectorE0dir must be perpendicular ALMOSTZERO(DOT_PRODUCT(WaveVector,WaveVectorE0dir))=.FALSE. .')
      ! set the user vector
      E_0_vec=WaveVectorE0dir
    END IF

    ! sanity check: perpendicularity
    IF(.NOT.ALMOSTZERO(DOT_PRODUCT(WaveVector,E_0_vec))) CALL abort(&
        __STAMP__&
        ,' WaveVector and E_0_vec must be perpendicular ALMOSTZERO(DOT_PRODUCT(WaveVector,E_0_vec)).')

    IF(RefStates(iRefState).NE.12)THEN
      ! ONLY FOR CASE(121,14,15,16)
      ! -------------------------------------------------------------------
      ! spatial Gaussian beam or plane wave, only in x,y or z direction
      ! additional tFWHM is a temporal Gaussian
      ! note:
      ! 121: Pulsed plane wave (infinite spot size) and temporal Gaussian
      ! 14 : Gaussian pulse is initialized IN the domain
      ! 15 : Gaussian pulse is a boundary condition, HENCE tDelayTime is used
      ! 16 : Gaussian pulse which is initialized in the domain and used as a boundary condition for t>0
      WaveBasePoint = GETREALARRAY('WaveBasePoint',3,'0.5 , 0.5 , 0.')

      ! Pulse duration
      sigma_t       = GETREAL ('sigma_t','0.')
      tFWHM         = GETREAL ('tFWHM','0.')
      IF((sigma_t.GT.0.0).AND.(ABS(tFWHM).LE.0.0))THEN
        tFWHM   = 2.*SQRT(2.*LOG(2.))*sigma_t
      ELSE IF((ABS(sigma_t).LE.0.0).AND.(tFWHM.GT.0.0))THEN
        sigma_t = tFWHM/(2.*SQRT(2.*LOG(2.)))
      ELSE
        CALL abort(&
            __STAMP__&
            ,' Input of pulse length is wrong.')
      END IF

      ! Beam spot size (waist radius, where the beam radius has a minimum)
      omega_0   = GETREAL ('omega_0','-1.')
      IF(omega_0.LE.0.0)THEN
        Beam_w0 = GETREAL ('Beam_w0','1.')
      ELSE
        Beam_w0 = omega_0
        LBWRITE(UNIT_StdOut,'(A)')'Setting Beam_w0=omega_0 (value read from old variable definition)'
      END IF
      Beam_w0_2inv = 2.0/(Beam_w0**2)
      sBeam_w0_2   = 1.0/(Beam_w0**2)

      ! Energy
      BeamEnergy    = GETREAL ('BeamEnergy','-1.0')
      BeamEta       = SQRT(mu0/eps0)

      IF((BeamEnergy.GT.0.0).AND.(RefStates(iRefState).NE.121))THEN ! only for 3D beams
        Beam_a0=-1.0
        I_0 = (BeamEnergy/(eps0*(Beam_w0**2)*c*((PI/2.)**(3./2.))*sigma_t*(EXP(-2*(c**2)*(BeamWaveNumber**2)*(sigma_t**2))+1)))&
            / (2*BeamEta)
        CALL PrintOption('calculated from BeamEnergy: I_0','CALCUL.',RealOpt=I_0)
      ELSE
        ! In 15: scaling by dimensionless laser amplitude Beam_a0 or optical intensity I_0
        Beam_a0 = GETREAL ('Beam_a0','-1.0')
        I_0     = GETREAL ('I_0','1.')
      END IF
      ! Decide if pulse maximum is scaled by intensity or a_0 parameter
      IF(Beam_a0.LE.0.0)THEN ! use I_0 for defining the amplitude
        Beam_a0    = 0.0
        E_0        = SQRT(2.0*BeamEta*I_0)
        CALL PrintOption('calculated from I_0: Beam_a0','CALCUL.',RealOpt=E_0*ElementaryCharge/(c*ElectronMass*BeamOmega))
        CALL PrintOption('calculated from I_0:     E_0','CALCUL.',RealOpt=E_0)
      ELSE ! use Beam_a0 for defining the amplitude
        E_0        = Beam_a0*c*ElectronMass*BeamOmega/ElementaryCharge
        CALL PrintOption('calculated from Beam_a0: I_0','CALCUL.',RealOpt=E_0**2/(2*BeamEta))
        CALL PrintOption('calculated from Beam_a0: E_0','CALCUL.',RealOpt=E_0)
      END IF


      IF(ALMOSTEQUAL(ABS(WaveVector(1)),1.))THEN ! wave in x-direction
        BeamIdir1   = 2
        BeamIdir2   = 3
        BeamMainDir = 1
      ELSE IF(ALMOSTEQUAL(ABS(WaveVector(2)),1.))THEN ! wave in y-direction
        BeamIdir1   = 1
        BeamIdir2   = 3
        BeamMainDir = 2
      ELSE IF(ALMOSTEQUAL(ABS(WaveVector(3)),1.))THEN! wave in z-direction
        BeamIdir1   = 1
        BeamIdir2   = 2
        BeamMainDir = 3
      ELSE
        CALL abort(&
            __STAMP__&
            ,'RefStates CASE(121,14,15,16): wave vector currently only in x,y,z!')
      END IF

      ! determine active time for time-dependent BC: save computational time for BC -> or possible switch to SM BC?
      SELECT CASE(RefStates(iRefState))
      CASE(121,15,16) ! pure BC or mixed IC+BC
        IF(RefStates(iRefState).EQ.16)THEN
          CALL PrintOption('tActive (old for BC=16)','CALCUL.',RealOpt=3*ABS(WaveBasePoint(BeamMainDir))*c_inv)
          ! get xyzMinMax
          CALL GetMeshMinMaxBoundaries()
          PulseCenter = WaveBasePoint(BeamMainDir) - (xyzMinMax(2*BeamMainDir)+xyzMinMax(2*BeamMainDir-1))/2
          IF((PulseCenter*WaveVector(BeamMainDir)).LT.0.0)THEN ! wave vector and base point are pointing in opposite direction
            tActive = (3./2.)*c_inv*(ABS(PulseCenter)+ABS((xyzMinMax(2*BeamMainDir)-xyzMinMax(2*BeamMainDir-1))/2))
          ELSE
            tActive = (1./2.)*c_inv*(ABS(PulseCenter)+ABS((xyzMinMax(2*BeamMainDir)-xyzMinMax(2*BeamMainDir-1))/2))
          END IF
        ELSE
          tActive = 8*sigma_t
        END IF
        CALL PrintOption('tActive (laser pulse time)','CALCUL.',RealOpt=tActive)
        CALL PrintOption('pulse will end at tActive+Beam_t0','CALCUL.',RealOpt=tActive+Beam_t0)
      END SELECT

      ! Determine total pulse energy
      SELECT CASE(RefStates(iRefState))
      CASE(121) ! Pulsed plane wave (pure BC or mixed IC+BC) with infinite spot size
        ! total beam energy in 2D is an energy per area -> [J/m^2]
        BeamEnergy_loc=eps0*(E_0**2)*SQRT(PI/2.0)*sigma_t*(EXP(-2.*(BeamOmega**2)*(sigma_t**2))+1)
        CALL PrintOption(' total beam energy per area [J/m^2]','CALCUL.',RealOpt=BeamEnergy_loc)
        CALL PrintOption('total beam energy per area [J/cm^2]','CALCUL.',RealOpt=BeamEnergy_loc/1.e4)
      CASE(14,15,16) ! 3D pulse with spot size
        ! total beam energy
        BeamEnergy_loc=(E_0**2)*PI*eps0*(Beam_w0**2)*SQRT(PI)*c*(sigma_t/(2.*SQRT(2.)))*&
            (EXP(-2*(c**2)*(BeamWaveNumber**2)*(sigma_t**2))+1)
        CALL PrintOption(' total beam energy [J]','CALCUL.',RealOpt=BeamEnergy_loc)

        ! beam spot area
        BeamArea_loc    = PI*(Beam_w0**2)
        CALL PrintOption('beam spot area (from waist radius) [m^2]','CALCUL.',RealOpt=BeamArea_loc)

        ! beam fluency
        BeamFluency_loc = BeamEnergy_loc / BeamArea_loc
        CALL PrintOption(' beam fluency [J/m^2]','CALCUL.',RealOpt=BeamFluency_loc)
        CALL PrintOption('beam fluency [J/cm^2]','CALCUL.',RealOpt=BeamFluency_loc/1.e4)
      END SELECT
    END IF
  END SELECT
END DO

DEALLOCATE(RefStates)

BCStateFile=GETSTR('BCStateFile','no file found')
!WRITE(DefBCState,'(I3,A,I3,A,I3,A,I3,A,I3,A,I3)') &
!  IniExactFunc,',',IniExactFunc,',',IniExactFunc,',',IniExactFunc,',',IniExactFunc,',',IniExactFunc
!IF(BCType_in(1) .EQ. -999)THEN
!  BCType = GETINTARRAY('BoundaryType',6)
!ELSE
!  BCType=BCType_in
!  SWRITE(UNIT_stdOut,*)'|                   BoundaryType | -> Already read in CreateMPICart!'

!END IF
!BCState   = GETINTARRAY('BoundaryState',6,TRIM(DefBCState))
!BoundaryCondition(:,1) = BCType
!BoundaryCondition(:,2) = BCState
! Read exponent for shape function
alpha_shape = GETINT('AlphaShape','2')
rCutoff     = GETREAL('r_cutoff','1.')
! Compute factor for shape function
ShapeFuncPrefix = 1./(2. * beta(1.5, REAL(alpha_shape) + 1.) * REAL(alpha_shape) + 2. * beta(1.5, REAL(alpha_shape) + 1.)) &
                * (REAL(alpha_shape) + 1.)/(PI*(rCutoff**3))

EquationInitIsDone=.TRUE.
LBWRITE(UNIT_stdOut,'(A)')' INIT MAXWELL DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitEquation



SUBROUTINE ExactFunc(ExactFunction,t_IN,tDeriv,x,resu)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars  ,ONLY: PI,c,c2,eps0,c_inv,c_inv
USE MOD_Equation_Vars ,ONLY: WaveVector,WaveBasePoint,sigma_t,E_0_vec,BeamIdir1,BeamIdir2,BeamMainDir,BeamWaveNumber
USE MOD_Equation_Vars ,ONLY: BeamOmega,E_0,TEScale,TEPulse,TEFrequency,TEPolarization,Beam_w0,TERadius,sBeam_w0_2
USE MOD_Equation_Vars ,ONLY: xDipole,tActive,TEModeRoot,Beam_t0,DoExactFlux,ExactFluxDir,ExactFluxPosition
USE MOD_Equation_Vars ,ONLY: TEMode
USE MOD_TimeDisc_Vars ,ONLY: dt
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: t_IN
INTEGER,INTENT(IN)              :: tDeriv           ! determines the time derivative of the function
REAL,INTENT(IN)                 :: x(3)
INTEGER,INTENT(IN)              :: ExactFunction    ! determines the exact function
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Resu(PP_nVar)    ! state in conservative variables
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: Resu_t(PP_nVar),Resu_tt(PP_nVar) ! state in conservative variables
REAL                            :: Frequency,Amplitude,Omega
REAL                            :: Cent(3),r,r2!,zlen
REAL                            :: a, b, d, l, m, nn, B0            ! aux. Variables for Resonator-Example
REAL                            :: gammaW,Psi,GradPsiX,GradPsiY     !     -"-
REAL                            :: xrel(3), theta, Etheta          ! aux. Variables for Dipole
REAL,PARAMETER                  :: Q=1, dD=1, omegaD=6.28318E8     ! aux. Constants for Dipole
REAL                            :: cos1,sin1,b1,b2                     ! aux. Variables for Gyrotron
REAL                            :: eps,phi,z                       ! aux. Variables for Gyrotron
REAL                            :: Er,Br,Ephi,Bphi,Bz,Ez           ! aux. Variables for Gyrotron
!REAL, PARAMETER                 :: k0=3562.936537
REAL                            :: SqrtN
REAL                            :: omegaG,g,h,B0G
REAL                            :: Bess_mG_R_R_inv,r_inv
REAL                            :: Bess_mG_R,Bess_mGM_R,Bess_mGP_R,costz,sintz,sin2,cos2,costz2,sintz2,dBess_mG_R
INTEGER                         :: MG,nG
REAL                            :: kz
REAL                            :: t ! local time
!===================================================================================================================================
! Apply time shift if needed
t=t_IN - Beam_t0 ! default: Beam_t0 = 0.

Cent=x
SELECT CASE (ExactFunction)
CASE(0) ! Particles
  Resu=0.
CASE(1) ! Constant
  Resu(1:3)=1.
  resu(4:6)=c_inv*resu(1:3)
  Resu(7:8)=0.
  Resu_t=0.
  Resu_tt=0.
CASE(2) ! Coaxial Waveguide
  z = x(3)
  resu      = 0.
  Frequency = TEFrequency
  Amplitude = TEScale
  !zlen     = 2.5
  r         = TERadius
  r2        = (x(1)*x(1)+x(2)*x(2))/r
  omega     = 2.*Pi*Frequency!/zlen ! shift beruecksichtigen
  kz        = (c_inv)
  resu(1)=( x(1))*sin(omega*(kz*z-c*t))/r2
  resu(2)=( x(2))*sin(omega*(kz*z-c*t))/r2
  resu(4)=(-x(2))*sin(omega*(kz*z-c*t))/(r2*c)
  resu(5)=( x(1))*sin(omega*(kz*z-c*t))/(r2*c)
  resu = Amplitude * resu

#if PP_TimeDiscMethod==1
  Resu_t=0.
  resu_t(1)=-omega*c*( x(1))*cos(omega*(z-c*t))/r2
  resu_t(2)=-omega*c*( x(2))*cos(omega*(z-c*t))/r2
  resu_t(4)=-omega*c*(-x(2))*cos(omega*(z-c*t))/(r2*c)
  resu_t(5)=-omega*c*( x(1))*cos(omega*(z-c*t))/(r2*c)
  Resu_tt=0.
  resu_tt(1)=-(omega*c)**2*( x(1))*sin(omega*(z-c*t))/r2
  resu_tt(2)=-(omega*c)**2*( x(2))*sin(omega*(z-c*t))/r2
  resu_tt(4)=-(omega*c)**2*(-x(2))*sin(omega*(z-c*t))/(r2*c)
  resu_tt(5)=-(omega*c)**2*( x(1))*sin(omega*(z-c*t))/(r2*c)
  resu_t = Amplitude * resu_t
  resu_tt = Amplitude * resu_tt
#endif /*PP_TimeDiscMethod==1*/

CASE(22) ! Coaxial Waveguide
  resu      = 0.
  Frequency = TEFrequency
  Amplitude = TEScale
  r         = TERadius
  r2        = (x(1)*x(1)+x(2)*x(2))/r
  omega     = 2.*Pi*Frequency
  resu(1)=( x(1))*sin(-omega*t)/r2
  resu(2)=( x(2))*sin(-omega*t)/r2
  resu(4)=(-x(2))*sin(-omega*t)/(r2*c)
  resu(5)=( x(1))*sin(-omega*t)/(r2*c)
  resu = Amplitude * resu

CASE(3) ! Resonator
  !special initial values
  !geometric parameters
  a=1.5; b=1.0; d=3.0
  !time parameters
  l=5.; m=4.; nn=3.; B0=1.
  IF(a.eq.0) CALL abort(__STAMP__,' Parameter a of resonator is zero!')
  IF(b.eq.0) CALL abort(__STAMP__,' Parameter b of resonator is zero!')
  IF(d.eq.0) CALL abort(__STAMP__,' Parameter d of resonator is zero!')
  omega = Pi*c*sqrt((m/a)**2+(nn/b)**2+(l/d)**2)
  gammaW = sqrt((omega/c)**2-(l*pi/d)**2)
  IF(gammaW.eq.0) CALL abort(__STAMP__,' gammaW is computed to zero!')
  Psi      =   B0          * cos((m*pi/a)*x(1)) * cos((nn*pi/b)*x(2))
  GradPsiX = -(B0*(m*pi/a) * sin((m*pi/a)*x(1)) * cos((nn*pi/b)*x(2)))
  GradPsiY = -(B0*(nn*pi/b) * cos((m*pi/a)*x(1)) * sin((nn*pi/b)*x(2)))

  resu(1)= (-omega/gammaW**2) * sin((l*pi/d)*x(3)) *(-GradPsiY)* sin(omega*t)
  resu(2)= (-omega/gammaW**2) * sin((l*pi/d)*x(3)) *  GradPsiX * sin(omega*t)
  resu(3)= 0.0
  resu(4)=(1/gammaW**2)*(l*pi/d) * cos((l*pi/d)*x(3)) * GradPsiX * cos(omega*t)
  resu(5)=(1/gammaW**2)*(l*pi/d) * cos((l*pi/d)*x(3)) * GradPsiY * cos(omega*t)
  resu(6)= Psi                  * sin((l*pi/d)*x(3))            * cos(omega*t)
  resu(7)=0.
  resu(8)=0.

CASE(4) ! Dipole
  resu(1:8) = 0.
  !RETURN
  eps=1e-10
  xrel    = x - xDipole
  r = SQRT(DOT_PRODUCT(xrel,xrel))
  IF (r.LT.eps) RETURN
  IF (xrel(3).GT.eps) THEN
    theta = ATAN(SQRT(xrel(1)**2+xrel(2)**2)/xrel(3))
  ELSE IF (xrel(3).LT.(-eps)) THEN
    theta = ATAN(SQRT(xrel(1)**2+xrel(2)**2)/xrel(3)) + pi
  ELSE
    theta = 0.5*pi
  END IF
  phi = ATAN2(xrel(2),xrel(1))
  !IF (xrel(1).GT.eps)      THEN  ! <-------------- OLD stuff, simply replaced with ATAN2() ... but not validated
  !  phi = ATAN(xrel(2)/xrel(1))
  !ELSE IF (xrel(1).LT.eps) THEN ! THIS DIVIDES BY ZERO ?!
  !  phi = ATAN(xrel(2)/xrel(1)) + pi
  !ELSE IF (xrel(2).GT.eps) THEN
  !  phi = 0.5*pi
  !ELSE IF (xrel(2).LT.eps) THEN
  !  phi = 1.5*pi
  !ELSE
  !  phi = 0.0                                                                                     ! Vorsicht: phi ist hier undef!
  !END IF

  Er = 2.*cos(theta)*Q*dD/(4.*pi*eps0) * ( 1./r**3*sin(omegaD*t-omegaD*r/c) + (omegaD/(c*r**2)*cos(omegaD*t-omegaD*r/c) ) )
  Etheta = sin(theta)*Q*dD/(4.*pi*eps0) * ( (1./r**3-omegaD**2/(c**2*r))*sin(omegaD*t-omegaD*r/c) &
          + (omegaD/(c*r**2)*cos(omegaD*t-omegaD* r/c) ) )
  Bphi = 1/(c2*eps0)*omegaD*sin(theta)*Q*dD/(4.*pi) &
       * ( - omegaD/(c*r)*sin(omegaD*t-omegaD*r/c) + 1./r**2*cos(omegaD*t-omegaD*r/c) )
  IF (ABS(phi).GT.eps) THEN
    resu(1)= sin(theta)*cos(phi)*Er + cos(theta)*cos(phi)*Etheta
    resu(2)= sin(theta)*sin(phi)*Er + cos(theta)*sin(phi)*Etheta
    resu(3)= cos(theta)         *Er - sin(theta)         *Etheta
    resu(4)=-sin(phi)*Bphi
    resu(5)= cos(phi)*Bphi
    resu(6)= 0.0
  ELSE
    resu(3)= cos(theta)         *Er - sin(theta)         *Etheta
  END IF

CASE(40) ! Dipole without initial condition
  resu(1:8) = 0.

CASE(5) ! Initialization of TE waves in a circular waveguide: Original formulation, for new definition, see CASE(6)
  ! Book: Springer
  ! Elektromagnetische Feldtheorie fuer Ingenieure und Physicker
  ! p. 500ff
  ! polarization:
  ! false - linear polarization
  ! true  - circular polarization
  r=SQRT(x(1)**2+x(2)**2)
  ! if a DOF is located in the origin, prevent division by zero ..
  phi = ATAN2(X(2),X(1))
  z=x(3)
  omegaG=2*PI*TEFrequency ! angular frequency

  ! TE_mG,nG
  mG=TEMode(1) ! azimuthal wave number
  nG=TEMode(2) ! radial wave number

  SqrtN=TEModeRoot/TERadius ! (7.412)
  ! axial wave number
  ! 1/c^2 omegaG^2 - kz^2=mu^2/ro^2
  kz=(omegaG*c_inv)**2-SqrtN**2 ! (7.413)
  IF(kz.LT.0)THEN
    SWRITE(UNIT_stdOut,'(A,ES25.14E3)')'(omegaG*c_inv)**2 = ',(omegaG*c_inv)**2
    SWRITE(UNIT_stdOut,'(A,ES25.14E3)')'SqrtN**2          = ',SqrtN**2
    SWRITE(UNIT_stdOut,'(A)')'  Maybe frequency too small?'
    CALL abort(__STAMP__,'kz=SQRT((omegaG*c_inv)**2-SqrtN**2), but the argument in negative!')
  END IF
  kz=SQRT(kz)
  ! precompute coefficients
  Bess_mG_R  = BESSEL_JN(mG  ,r*SqrtN)
  Bess_mGM_R = BESSEL_JN(mG-1,r*SqrtN)
  Bess_mGP_R = BESSEL_JN(mG+1,r*SqrtN)
  dBess_mG_R = 0.5*(Bess_mGM_R-Bess_mGP_R)
  COSTZ      = COS(kz*z-omegaG*t)
  SINTZ      = SIN(kz*z-omegaG*t)
  sin1       = SIN(REAL(mG)*phi)
  cos1       = COS(REAL(mG)*phi)
  ! barrier for small radii
  IF(r/TERadius.LT.1e-4)THEN
    SELECT CASE(mG)
    CASE(0) ! arbitary
      Bess_mG_R_R_inv=1e6
    CASE(1)
      Bess_mG_R_R_inv=0.5
    CASE DEFAULT
      Bess_mG_R_R_inv=0.
    END SELECT
  ELSE
    r_inv=1./r
    Bess_mG_R_R_inv=Bess_mG_R*r_inv
  END IF
  ! Check polarization
  IF((TRIM(TEPolarization).EQ.'x').OR.(TRIM(TEPolarization).EQ.'y'))THEN ! linear polarization along the x/y-axis
    ! electric field
    Er   =  omegaG*REAL(mG)* Bess_mG_R_R_inv*sin1*SINTZ
    Ephi =  omegaG*SqrtN*dBess_mG_R*cos1*SINTZ
    Ez   =  0.
    ! magnetic field
    Br   = -kz*SqrtN*dBess_mG_R*cos1*SINTZ
    Bphi =  kz*REAL(mG)*Bess_mG_R_R_inv*sin1*SINTZ
    Bz   =  (SqrtN**2)*Bess_mG_R*cos1*COSTZ
  ELSE ! Circular polarization
    ! polarisation if superposition of two fields
    ! circular polarisation requires an additional temporal shift
    ! a) perpendicular shift of TE mode, rotation of 90 degree
    sin2       = SIN(REAL(mG)*phi+0.5*PI)
    cos2       = COS(REAL(mG)*phi+0.5*PI)
    IF(TRIM(TEPolarization).EQ.'l')THEN ! shift for left or right rotating fields
      COSTZ2     = COS(kz*z-omegaG*t-0.5*PI)
      SINTZ2     = SIN(kz*z-omegaG*t-0.5*PI)
    ELSE
      COSTZ2     = COS(kz*z-omegaG*t+0.5*PI)
      SINTZ2     = SIN(kz*z-omegaG*t+0.5*PI)
    END IF
    ! electric field
    Er   =  omegaG*REAL(mG)* Bess_mG_R_R_inv*(sin1*SINTZ+sin2*SINTZ2)
    Ephi =  omegaG*SqrtN*dBess_mG_R*(cos1*SINTZ+cos2*SINTZ2)
    Ez   =  0.
    ! magnetic field
    Br   = -kz*SqrtN*dBess_mG_R*(cos1*SINTZ+cos2*SINTZ2)
    Bphi =  kz*REAL(mG)*Bess_mG_R_R_inv*(sin1*SINTZ+sin2*SINTZ2)
    ! caution: does we have to modify the z entry? yes
    Bz   =  (SqrtN**2)*Bess_mG_R*(cos1*COSTZ+cos2*COSTZ2)
  END IF

  resu(1)= COS(phi)*Er - SIN(phi)*Ephi
  resu(2)= SIN(phi)*Er + COS(phi)*Ephi
  resu(3)= 0.0
  resu(4)= COS(phi)*Br - SIN(phi)*Bphi
  resu(5)= SIN(phi)*Br + COS(phi)*Bphi
  resu(6)= Bz
  resu(1:5)=resu(1:5)
  resu( 6 )=resu( 6 )
  resu(1:6)=TEScale*resu(1:6)
  resu(7)= 0.0
  resu(8)= 0.0
  IF(TEPulse)THEN
    sigma_t=4.*(2.*PI)/omegaG/(2.*SQRT(2.*LOG(2.)))
    IF (t.LE.34*sigma_t) THEN
      ASSOCIATE( t => t-4.*sigma_t )
        ASSOCIATE( temporalWindow => EXP(-0.5*(t/sigma_t)**2) )
          resu(1:8)=resu(1:8)*temporalWindow
        END ASSOCIATE
      END ASSOCIATE
    ELSE
      resu(1:8)=0.
    END IF
  END IF

CASE(6) ! Circular waveguide TEM[n,m]
  CALL ExactFunc_TE_Circular_Waveguide(t,x,resu)

CASE(7) ! Manufactured Solution
  resu(:)=0
  resu(1)=SIN(2*pi*(x(1)-t))
  resu_t(:)=0
  resu_t(1)=-2*pi*COS(2*pi*(x(1)-t))
  resu_tt(:)=0
  resu_tt(1)=-4*pi*pi*resu(1)

CASE(10) !issautier 3D test case with source (Stock et al., divcorr paper), domain [0;1]^3!!!
  resu(:)=0.
  resu(1)=x(1)*SIN(Pi*x(2))*SIN(Pi*x(3)) !*SIN(t)
  resu(2)=x(2)*SIN(Pi*x(3))*SIN(Pi*x(1)) !*SIN(t)
  resu(3)=x(3)*SIN(Pi*x(1))*SIN(Pi*x(2)) !*SIN(t)
  resu(4)=pi*SIN(Pi*x(1))*(x(3)*COS(Pi*x(2))-x(2)*COS(Pi*x(3))) !*(COS(t)-1)
  resu(5)=pi*SIN(Pi*x(2))*(x(1)*COS(Pi*x(3))-x(3)*COS(Pi*x(1))) !*(COS(t)-1)
  resu(6)=pi*SIN(Pi*x(3))*(x(2)*COS(Pi*x(1))-x(1)*COS(Pi*x(2))) !*(COS(t)-1)

  resu_t(:)=0.
  resu_t(1)= COS(t)*resu(1)
  resu_t(2)= COS(t)*resu(2)
  resu_t(3)= COS(t)*resu(3)
  resu_t(4)=-SIN(t)*resu(4)
  resu_t(5)=-SIN(t)*resu(5)
  resu_t(6)=-SIN(t)*resu(6)
  resu_tt=0.
  resu_tt(1)=-SIN(t)*resu(1)
  resu_tt(2)=-SIN(t)*resu(2)
  resu_tt(3)=-SIN(t)*resu(3)
  resu_tt(4)=-COS(t)*resu(4)
  resu_tt(5)=-COS(t)*resu(5)
  resu_tt(6)=-COS(t)*resu(6)

  resu(1)=     SIN(t)*resu(1)
  resu(2)=     SIN(t)*resu(2)
  resu(3)=     SIN(t)*resu(3)
  resu(4)=(COS(t)-1.)*resu(4)
  resu(5)=(COS(t)-1.)*resu(5)
  resu(6)=(COS(t)-1.)*resu(6)

CASE(12) ! planar wave test case
  resu(1:3)=E_0_vec*COS(BeamWaveNumber*DOT_PRODUCT(WaveVector,x)-BeamOmega*t)
  resu(4:6)=c_inv*CROSS(WaveVector,resu(1:3))
  resu(7:8)=0.
CASE(121) ! like CASE(12), but with pulsed wave: boundary condition (BC)
  IF(t.GT.tActive)THEN ! pulse has passesd -> return
    resu(1:8)=0.
  ELSE
    ASSOCIATE( t => t - 4*sigma_t , & ! t: add arbitrary time shift
               k => BeamWaveNumber, & ! k: wave number
               w => BeamOmega)        ! w: angular frequency
      resu(1:3)=E_0*E_0_vec*COS(k*DOT_PRODUCT(WaveVector,x-WaveBasePoint)-w*t)*EXP(-0.25*(t/sigma_t)**2)
      resu(4:6)=c_inv*CROSS(WaveVector,resu(1:3))
      resu(7:8)=0.
    END ASSOCIATE
  END IF

CASE(14) ! 1 of 3: Gauss-shape with perfect focus (w(z)=w_0): initial condition (IC)
         ! spatial gauss beam, still planar wave scaled by intensity spatial and temporal filter are defined according to
         ! Thiele 2016: "Modelling laser matter interaction with tightly focused laser pules in electromagnetic codes"
         ! beam insert is done by a paraxial assumption focus is at basepoint
         ! intensity * Gaussian filter in transversal and longitudinal direction
         ASSOCIATE( spatialWindow  => EXP(    -((x(BeamIdir1)-WaveBasePoint(BeamIdir1))**2+&
                                                (x(BeamIdir2)-WaveBasePoint(BeamIdir2))**2)/((  Beam_w0  )**2)&
                                              -((x(BeamMainDir)-WaveBasePoint(BeamMainDir))**2)/((2*sigma_t*c)**2)  ) , &
                    timeFac        => COS(  BeamWaveNumber*DOT_PRODUCT(WaveVector,x-WaveBasePoint)-BeamOmega*&
                                            (t-ABS(WaveBasePoint(BeamMainDir))/c)  )                                    )
           resu(1:3)=E_0*spatialWindow*E_0_vec*timeFac
           resu(4:6)=c_inv*CROSS( WaveVector,resu(1:3))
           resu(7:8)=0.
         END ASSOCIATE
CASE(15) ! 2 of 3: Gauss-shape with perfect focus (w(z)=w_0): boundary condition (BC)
         ! spatial gauss beam, still planar wave scaled by intensity spatial and temporal filter are defined according to
         ! Thiele 2016: "Modelling laser matter interaction with tightly focused laser pules in electromagnetic codes"
         ! beam insert is done by a paraxial assumption focus is at base point and should be on BC
  IF(t.GT.tActive)THEN ! Pulse has passed -> return
    resu(1:8)=0.
  ELSE
    ASSOCIATE( t => t - 4*sigma_t ) ! t: add (arbitrary) time shift
      ASSOCIATE( spatialWindow  => EXP(-((x(BeamIdir1)-WaveBasePoint(BeamIdir1))**2+&
                                   (x(BeamIdir2)-WaveBasePoint(BeamIdir2))**2)*sBeam_w0_2) , & ! spatial window in (x^2+y^2)/(w_0^2)
                 timeFac        => COS(BeamWaveNumber*DOT_PRODUCT(WaveVector,x-WaveBasePoint)-BeamOmega*t) , & ! COS() function
                 temporalWindow => EXP(-0.25*(t/sigma_t)**2) ) ! temporal Gaussian window
        resu(1:3)=E_0*spatialWindow*E_0_vec*timeFac*temporalWindow
        resu(4:6)=c_inv*CROSS( WaveVector,resu(1:3))
        resu(7:8)=0.
      END ASSOCIATE
    END ASSOCIATE
  END IF
CASE(16) ! 3 of 3: Gauss-shape with perfect focus (w(z)=w_0): initial & boundary condition (BC)
         ! spatial Gauss beam, still planar wave scaled by intensity spatial and temporal filter are defined according to
         ! Thiele 2016: "Modelling laser matter interaction with tightly focused laser pules in electromagnetic codes"
         ! beam insert is done by a paraxial assumption focus is at base point and should be on BC
  IF(t.GT.tActive)THEN ! Pulse has passed -> return
    resu(1:8)=0.
  ELSE ! IC (t=0) or BC (t>0)
    ASSOCIATE( tShift   => t-ABS(WaveBasePoint(BeamMainDir))/c             , & ! substitution: shift to wave base point position
               tShiftBC => t+(WaveBasePoint(BeamMainDir)-x(BeamMainDir))/c   ) ! shift to wave base point position
      IF(t.LT.dt)THEN ! Initial condition: IC
        ASSOCIATE( spatialWindow => EXP(-((x(BeamIdir1)-WaveBasePoint(BeamIdir1))**2+              &
                                          (x(BeamIdir2)-WaveBasePoint(BeamIdir2))**2)*sBeam_w0_2   & ! (x^2+y^2)/(w_0^2)
                                    -((x(BeamMainDir)-WaveBasePoint(BeamMainDir))**2)/((2*sigma_t*c)**2))  , &
                   timeFac       => COS(BeamWaveNumber*DOT_PRODUCT(WaveVector,x-WaveBasePoint)-BeamOmega*tShift) ) ! COS() function
          ! For setting the correct initial condition when using ExactFlux
          IF(DoExactFlux)THEN
            IF(ExactFluxDir.EQ.BeamMainDir)THEN
              ASSOCIATE( plusminus => WaveVector(BeamMainDir)/ABS(WaveVector(BeamMainDir)) )
                ! Set the field "in front of" the ExactFlux plane to the IC condition (compatible with the BC)
                IF( (x(BeamMainDir)*plusminus.GE.ExactFluxPosition*plusminus) .OR.&
                    (ALMOSTEQUALRELATIVE(x(BeamMainDir),ExactFluxPosition,1.0E-4)) )THEN

                  resu(1:3)=E_0*spatialWindow*E_0_vec*timeFac
                ELSE ! Set the field "behind" the ExactFlux plane to zero
                  resu(1:3)=0.
                END IF
              END ASSOCIATE
            ELSE
              CALL abort(__STAMP__,'ExactFunction=16 (laser pulse IC+BC) together with ExactFlux can only '//&
                                   'be used with ExactFluxDir=BeamMainDir (in BeamMainDir-direction).')
            END IF
          ELSE ! Default: no ExactFlux is used
            resu(1:3)=E_0*spatialWindow*E_0_vec*timeFac
          END IF
        END ASSOCIATE
      ELSE ! Boundary condition: BC
        ASSOCIATE( spatialWindow  => EXP(-((x(BeamIdir1)-WaveBasePoint(BeamIdir1))**2+&
                                          (x(BeamIdir2)-WaveBasePoint(BeamIdir2))**2)*sBeam_w0_2) , & ! (x^2+y^2)/(w_0^2)
                   temporalWindow => EXP(-0.25*(tShiftBC/sigma_t)**2)                             , &
                   timeFac        => COS(BeamWaveNumber*DOT_PRODUCT(WaveVector,x-WaveBasePoint)-BeamOmega*tShift) ) ! COS() function
          resu(1:3)=E_0*spatialWindow*E_0_vec*timeFac*temporalWindow
        END ASSOCIATE
      END IF
      resu(4:6)=c_inv*CROSS(WaveVector,resu(1:3))
      resu(7:8)=0.
    END ASSOCIATE
  END IF
CASE(50,51)            ! Initialization and BC Gyrotron - including derivatives
  g      = 3236.706462    ! aux. Constants for Gyrotron
  B0G    = 1.0
  h      = 1489.378411    ! aux. Constants for Gyrotron
  omegaG = 3.562936537e+3 ! aux. Constants for Gyrotron
  eps    = 1e-10
  mG     = 34
  IF ((ExactFunction.EQ.51).AND.(x(3).GT.eps)) RETURN
  r=SQRT(x(1)**2+x(2)**2)
  phi = ATAN2(x(2),x(1))
 !    IF (x(1).GT.eps)      THEN ! <-------------- OLD stuff, simply replaced with ATAN2() ... but not validated
 !      phi = ATAN(x(2)/x(1))
 !    ELSE IF (x(1).LT.(-eps)) THEN
 !      phi = ATAN(x(2)/x(1)) + pi
 !    ELSE IF (x(2).GT.eps) THEN
 !      phi = 0.5*pi
 !    ELSE IF (x(2).LT.(-eps)) THEN
 !      phi = 1.5*pi
 !    ELSE
 !      phi = 0.0         ! Vorsicht: phi ist hier undef!
 !    END IF
  z = x(3)
  a = h*z+mG*phi
  b0 = BESSEL_JN(mG,REAL(g*r))
  b1 = BESSEL_JN(mG-1,REAL(g*r))
  b2 = BESSEL_JN(mG+1,REAL(g*r))
  SELECT CASE(MOD(tDeriv,4))
    CASE(0)
      cos1  =  omegaG**tDeriv * cos(a-omegaG*t)
      sin1  =  omegaG**tDeriv * sin(a-omegaG*t)
    CASE(1)
      cos1  =  omegaG**tDeriv * sin(a-omegaG*t)
      sin1  = -omegaG**tDeriv * cos(a-omegaG*t)
    CASE(2)
      cos1  = -omegaG**tDeriv * cos(a-omegaG*t)
      sin1  = -omegaG**tDeriv * sin(a-omegaG*t)
    CASE(3)
      cos1  = -omegaG**tDeriv * sin(a-omegaG*t)
      sin1  =  omegaG**tDeriv * cos(a-omegaG*t)
    CASE DEFAULT
      cos1  = 0.0
      sin1  = 0.0
      CALL abort(&
          __STAMP__&
          ,'What is that weired tDeriv you gave me?',999,999.)
  END SELECT

  Er  =-B0G*mG*omegaG/(r*g**2)*b0     *cos1
  Ephi= B0G*omegaG/h      *0.5*(b1-b2)*sin1
  Br  =-B0G*h/g           *0.5*(b1-b2)*sin1
  Bphi=-B0G*mG*h/(r*g**2)     *b0     *cos1
  Bz  = B0G                   *b0     *cos1
  resu(1)= cos(phi)*Er - sin(phi)*Ephi
  resu(2)= sin(phi)*Er + cos(phi)*Ephi
  resu(3)= 0.0
  resu(4)= cos(phi)*Br - sin(phi)*Bphi
  resu(5)= sin(phi)*Br + cos(phi)*Bphi
  resu(6)= Bz
  resu(7)= 0.0
  resu(8)= 0.0

CASE(41) ! pulsed Dipole
  resu = 0.0
  RETURN
CASE(100) ! QDS
  resu = 0.0
  RETURN
CASE DEFAULT
  SWRITE(*,*)'Exact function not specified. ExactFunction = ',ExactFunction
END SELECT ! ExactFunction

# if (PP_TimeDiscMethod==1)
! For O3 RK, the boundary condition has to be adjusted
! Works only for O3 RK!!
SELECT CASE(tDeriv)
CASE(0)
  ! resu = g(t)
CASE(1)
  ! resu = g(t) + dt/3*g'(t)
  Resu=Resu + dt/3.*Resu_t
CASE(2)
  ! resu = g(t) + 3/4 dt g'(t) +5/16 dt^2 g''(t)
  Resu=Resu + 0.75*dt*Resu_t+5./16.*dt*dt*Resu_tt
CASE DEFAULT
  ! Stop, works only for 3 Stage O3 LS RK
  CALL abort(__STAMP__,'Exactfuntion works only for 3 Stage O3 LS RK!',999,999.)
END SELECT
#endif
END SUBROUTINE ExactFunc


SUBROUTINE CalcSource(t,coeff,Ut)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals           ,ONLY: abort
USE MOD_Globals_Vars      ,ONLY: PI,eps0
USE MOD_PreProc
USE MOD_Equation_Vars     ,ONLY: c_corr,IniExactFunc, DipoleOmega,tPulse,xDipole
#ifdef PARTICLES
USE MOD_PICDepo_Vars      ,ONLY: PartSource,DoDeposition
USE MOD_Dielectric_Vars   ,ONLY: DoDielectric,isDielectricElem,ElemToDielectric,DielectricEps,ElemToDielectric
#if IMPA
USE MOD_LinearSolver_Vars ,ONLY: ExplicitPartSource
#endif
#endif /*PARTICLES*/
USE MOD_Mesh_Vars         ,ONLY: Elem_xGP
#if defined(LSERK) || defined(IMPA) || defined(ROS)
USE MOD_Equation_Vars     ,ONLY: DoParabolicDamping,fDamping
USE MOD_TimeDisc_Vars     ,ONLY: sdtCFLOne
USE MOD_DG_Vars           ,ONLY: U
#endif /*LSERK*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: t,coeff
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)              :: Ut(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iElem
REAL                            :: eps0inv, x(1:3)
REAL                            :: r           ! for Dipole
REAL,PARAMETER                  :: Q=1, d=1    ! for Dipole
#ifdef PARTICLES
REAL                            :: PartSourceLoc(1:4)
#endif
REAL                            :: coeff_loc
!===================================================================================================================================
eps0inv = 1./eps0
#ifdef PARTICLES
IF(DoDeposition)THEN
  IF(DoDielectric)THEN
    DO iElem=1,PP_nElems
      IF(isDielectricElem(iElem)) THEN ! 1.) PML version - PML element
        DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
#if IMPA
          PartSourceLoc=PartSource(:,i,j,k,iElem)+ExplicitPartSource(:,i,j,k,iElem)
#else
          PartSourceLoc=PartSource(:,i,j,k,iElem)
#endif
          !  Get PartSource from Particles
          !Ut(1:3,i,j,k,iElem) = Ut(1:3,i,j,k,iElem) - eps0inv *coeff* PartSource(1:3,i,j,k,iElem) * DielectricEpsR_inv
          !Ut(  8,i,j,k,iElem) = Ut(  8,i,j,k,iElem) + eps0inv *coeff* PartSource(  4,i,j,k,iElem) * c_corr * DielectricEpsR_inv
          Ut(1:3,i,j,k,iElem) = Ut(1:3,i,j,k,iElem) - eps0inv *coeff* PartSourceloc(1:3) &
                                                      / DielectricEps(i,j,k,ElemToDielectric(iElem)) ! only use x
          Ut(  8,i,j,k,iElem) = Ut(  8,i,j,k,iElem) + eps0inv *coeff* PartSourceloc( 4 ) * c_corr &
                                                      / DielectricEps(i,j,k,ElemToDielectric(iElem)) ! only use x
        END DO; END DO; END DO
      ELSE
        DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
#if IMPA
          PartSourceLoc=PartSource(:,i,j,k,iElem)+ExplicitPartSource(:,i,j,k,iElem)
#else
          PartSourceLoc=PartSource(:,i,j,k,iElem)
#endif
          !  Get PartSource from Particles
          Ut(1:3,i,j,k,iElem) = Ut(1:3,i,j,k,iElem) - eps0inv *coeff* PartSourceloc(1:3)
          Ut(  8,i,j,k,iElem) = Ut(  8,i,j,k,iElem) + eps0inv *coeff* PartSourceloc( 4 ) * c_corr
        END DO; END DO; END DO
      END IF
    END DO
  ELSE
    DO iElem=1,PP_nElems
      DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
#if IMPA
        PartSourceLoc=PartSource(:,i,j,k,iElem)+ExplicitPartSource(:,i,j,k,iElem)
#else
        PartSourceLoc=PartSource(:,i,j,k,iElem)
#endif
        !  Get PartSource from Particles
        Ut(1:3,i,j,k,iElem) = Ut(1:3,i,j,k,iElem) - eps0inv *coeff* PartSourceloc(1:3)
        Ut(  8,i,j,k,iElem) = Ut(  8,i,j,k,iElem) + eps0inv *coeff* PartSourceloc( 4 ) * c_corr
      END DO; END DO; END DO
    END DO
  END IF
END IF
#endif /*PARTICLES*/
SELECT CASE (IniExactFunc)
CASE(0) ! Particles
  ! empty, nothing to do
CASE(1) ! Constant          - no sources
CASE(2,22) ! Coaxial Waveguide - no sources
CASE(3) ! Resonator         - no sources
CASE(4) ! Dipole
  DO iElem=1,PP_nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      r = SQRT(DOT_PRODUCT(Elem_xGP(:,i,j,k,iElem)-xDipole,Elem_xGP(:,i,j,k,iElem)-xDipole))
      IF (shapefunc(r) .GT. 0 ) THEN
        Ut(3,i,j,k,iElem) = Ut(3,i,j,k,iElem) - (shapefunc(r)) *coeff* Q*d*DipoleOmega * COS(DipoleOmega*t) * eps0inv
    ! dipole should be neutral
        Ut(8,i,j,k,iElem) = Ut(8,i,j,k,iElem) + (shapefunc(r)) *coeff* c_corr*Q*d*SIN(DipoleOmega*t) * eps0inv
      END IF
    END DO; END DO; END DO
  END DO
CASE(40) ! Dipole without initial condition
  coeff_loc = 1.0e-11 ! amplitude scaling
  DO iElem=1,PP_nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      r = SQRT(DOT_PRODUCT(Elem_xGP(:,i,j,k,iElem)-xDipole,Elem_xGP(:,i,j,k,iElem)-xDipole))
      IF (shapefunc(r) .GT. 0 ) THEN
        Ut(3,i,j,k,iElem) = Ut(3,i,j,k,iElem) - (shapefunc(r)) *coeff_loc* Q*d*DipoleOmega * COS(DipoleOmega*t) * eps0inv
    ! dipole should be neutral
        Ut(8,i,j,k,iElem) = Ut(8,i,j,k,iElem) + (shapefunc(r)) *coeff_loc* c_corr*Q*d*SIN(DipoleOmega*t) * eps0inv
      END IF
    END DO; END DO; END DO
  END DO
CASE(5) ! TE_34,19 Mode     - no sources
CASE(7) ! Manufactured Solution
  DO iElem=1,PP_nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      Ut(1,i,j,k,iElem) =Ut(1,i,j,k,iElem) - coeff*2*pi*COS(2*pi*(Elem_xGP(1,i,j,k,iElem)-t)) * eps0inv
      Ut(8,i,j,k,iElem) =Ut(8,i,j,k,iElem) + coeff*2*pi*COS(2*pi*(Elem_xGP(1,i,j,k,iElem)-t)) * c_corr * eps0inv
    END DO; END DO; END DO
  END DO
CASE(10) !issautier 3D test case with source (Stock et al., divcorr paper), domain [0;1]^3!!!
  DO iElem=1,PP_nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      x(:)=Elem_xGP(:,i,j,k,iElem)
      Ut(1,i,j,k,iElem) =Ut(1,i,j,k,iElem) + coeff*(COS(t)- (COS(t)-1.)*2*pi*pi)*x(1)*SIN(Pi*x(2))*SIN(Pi*x(3))
      Ut(2,i,j,k,iElem) =Ut(2,i,j,k,iElem) + coeff*(COS(t)- (COS(t)-1.)*2*pi*pi)*x(2)*SIN(Pi*x(3))*SIN(Pi*x(1))
      Ut(3,i,j,k,iElem) =Ut(3,i,j,k,iElem) + coeff*(COS(t)- (COS(t)-1.)*2*pi*pi)*x(3)*SIN(Pi*x(1))*SIN(Pi*x(2))
      Ut(1,i,j,k,iElem) =Ut(1,i,j,k,iElem) - coeff*(COS(t)-1.)*pi*COS(Pi*x(1))*(SIN(Pi*x(2))+SIN(Pi*x(3)))
      Ut(2,i,j,k,iElem) =Ut(2,i,j,k,iElem) - coeff*(COS(t)-1.)*pi*COS(Pi*x(2))*(SIN(Pi*x(3))+SIN(Pi*x(1)))
      Ut(3,i,j,k,iElem) =Ut(3,i,j,k,iElem) - coeff*(COS(t)-1.)*pi*COS(Pi*x(3))*(SIN(Pi*x(1))+SIN(Pi*x(2)))
      Ut(8,i,j,k,iElem) =Ut(8,i,j,k,iElem) + coeff*c_corr*SIN(t)*( SIN(pi*x(2))*SIN(pi*x(3)) &
                                                            +SIN(pi*x(3))*SIN(pi*x(1)) &
                                                            +SIN(pi*x(1))*SIN(pi*x(2)) )
    END DO; END DO; END DO
  END DO

CASE(12) ! plane wave
CASE(14) ! gauss pulse, spatial -> IC
CASE(15) ! gauss pulse, temporal -> BC
CASE(16) ! gauss pulse, temporal -> IC+BC

CASE(41) ! Dipole via temporal Gausspuls
!t0=TEnd/5, w=t0/4 ! for pulsed Dipole (t0=offset and w=width of pulse)
!TEnd=30.E-9 -> short pulse for 100ns runtime
IF(1.EQ.2)THEN ! new formulation with divergence correction considered
  DO iElem=1,PP_nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      Ut(1,i,j,k,iElem) =Ut(1,i,j,k,iElem) - coeff*2*pi*COS(2*pi*(Elem_xGP(1,i,j,k,iElem)-t)) * eps0inv
      Ut(8,i,j,k,iElem) =Ut(8,i,j,k,iElem) + coeff*2*pi*COS(2*pi*(Elem_xGP(1,i,j,k,iElem)-t)) * c_corr * eps0inv
    END DO; END DO; END DO
  END DO
ELSE ! old/original formulation
  DO iElem=1,PP_nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    IF (t.LE.2*tPulse) THEN
      r = SQRT(DOT_PRODUCT(Elem_xGP(:,i,j,k,iElem)-xDipole,Elem_xGP(:,i,j,k,iElem)-xDipole))
      IF (shapefunc(r) .GT. 0 ) THEN
        Ut(3,i,j,k,iElem) = Ut(3,i,j,k,iElem) - ((shapefunc(r))*Q*d*COS(DipoleOmega*t)*eps0inv)*&
                            EXP(-(t-tPulse/5)**2/(2*(tPulse/(4*5))**2))
      END IF
    END IF
  END DO; END DO; END DO; END DO
END IF
CASE(50,51) ! TE_34,19 Mode - no sources
CASE DEFAULT
  CALL abort(__STAMP__,'Exactfunction not specified! IniExactFunc = ',IntInfoOpt=IniExactFunc)
END SELECT ! ExactFunction

#if defined(LSERK) ||  defined(ROS) || defined(IMPA)
IF(DoParabolicDamping)THEN
  !Ut(7:8,:,:,:,:) = Ut(7:8,:,:,:,:) - (1.0-fDamping)*sdtCFLOne/RK_b(iStage)*U(7:8,:,:,:,:)
  Ut(7:8,:,:,:,:) = Ut(7:8,:,:,:,:) - (1.0-fDamping)*sdtCFLOne*U(7:8,:,:,:,:)
END IF
#endif /*LSERK*/

END SUBROUTINE CalcSource


SUBROUTINE DivCleaningDamping()
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars,       ONLY : U
USE MOD_Equation_Vars, ONLY : fDamping,DoParabolicDamping
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iElem
!===================================================================================================================================
IF(DoParabolicDamping) RETURN
DO iElem=1,PP_nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    !  Get source from Particles
    U(7:8,i,j,k,iElem) = U(7:8,i,j,k,iElem) * fDamping
  END DO; END DO; END DO
END DO
END SUBROUTINE DivCleaningDamping

FUNCTION shapefunc(r)
!===================================================================================================================================
! Implementation of (possibly several different) shapefunctions
!===================================================================================================================================
! MODULES
  USE MOD_Equation_Vars, ONLY : shapeFuncPrefix, alpha_shape, rCutoff
! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
    REAL                 :: r         ! radius / distance to center
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
    REAL                 :: shapefunc ! sort of a weight for the source
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
   IF (r.GE.rCutoff) THEN
     shapefunc = 0.0
   ELSE
     shapefunc = ShapeFuncPrefix *(1-(r/rCutoff)**2)**alpha_shape
   END IF
END FUNCTION shapefunc

FUNCTION beta(z,w)
   IMPLICIT NONE
   REAL beta, w, z
   beta = GAMMA(z)*GAMMA(w)/GAMMA(z+w)
END FUNCTION beta


SUBROUTINE GetWaveGuideRadius(DoSide)
!===================================================================================================================================
! routine to find the maximum radius of a  wave-guide at a given BC plane
! radius computation requires interpolation points on the surface, hence
! an additional change-basis is required to map Gauss to Gauss-Lobatto points
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars    ,  ONLY:nSides,Face_xGP
USE MOD_Equation_Vars,  ONLY:TERadius
#if (PP_NodeType==1)
USE MOD_ChangeBasis,    ONLY:ChangeBasis2D
USE MOD_Basis,          ONLY:LegGaussLobNodesAndWeights
USE MOD_Basis,          ONLY:BarycentricWeights,InitializeVandermonde
USE MOD_Interpolation_Vars, ONLY:xGP,wBary
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
LOGICAL,INTENT(IN)      :: DoSide(1:nSides)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Radius
INTEGER                 :: iSide,p,q
#if (PP_NodeType==1)
REAL                    :: xGP_tmp(0:PP_N),wBary_tmp(0:PP_N),wGP_tmp(0:PP_N)
REAL                    :: Vdm_PolN_GL(0:PP_N,0:PP_N)
#endif
REAL                    :: Face_xGL(1:2,0:PP_N,0:PP_N)
!===================================================================================================================================

#if (PP_NodeType==1)
! get Vandermonde, change from Gauss or Gauss-Lobatto Points to Gauss-Lobatto-Points
! radius requires GL-points
CALL LegGaussLobNodesAndWeights(PP_N,xGP_tmp,wGP_tmp)
CALL BarycentricWeights(PP_N,xGP_tmp,wBary_tmp)
!CALL InitializeVandermonde(PP_N,PP_N,wBary_tmp,xGP,xGP_tmp,Vdm_PolN_GL)
CALL InitializeVandermonde(PP_N,PP_N,wBary,xGP,xGP_tmp,Vdm_PolN_GL)
#endif

TERadius=0.
Radius   =0.
DO iSide=1,nSides
  IF(.NOT.DoSide(iSide)) CYCLE
#if (PP_NodeType==1)
  CALL ChangeBasis2D(2,PP_N,PP_N,Vdm_PolN_GL,Face_xGP(1:2,:,:,iSide),Face_xGL)
#else
  Face_xGL(1:2,:,:)=Face_xGP(1:2,:,:,iSide)
#endif
  DO q=0,PP_N
    DO p=0,PP_N
      Radius=SQRT(Face_xGL(1,p,q)**2+Face_xGL(2,p,q)**2)
      TERadius=MAX(Radius,TERadius)
    END DO ! p
  END DO ! q
END DO

#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,TERadius,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/

LBWRITE(UNIT_StdOut,*) ' Found waveguide radius of ', TERadius

END SUBROUTINE GetWaveGuideRadius


SUBROUTINE InitExactFlux()
!===================================================================================================================================
! Get the constant advection velocity vector from the ini file
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals          ,ONLY: abort,UNIT_stdOut,mpiroot,CollectiveStop
#if USE_MPI
USE MOD_Globals          ,ONLY: MPI_COMM_PICLAS,MPI_SUM,MPI_INTEGER,IERROR
#endif
USE MOD_Mesh_Vars        ,ONLY: nElems,ElemToSide,SideToElem,lastMPISide_MINE
USE MOD_Interfaces       ,ONLY: FindElementInRegion,FindInterfacesInRegion,CountAndCreateMappings
USE MOD_Equation_Vars    ,ONLY: ExactFluxDir,ExactFluxPosition,isExactFluxInterFace
USE MOD_ReadInTools      ,ONLY: GETREAL,GETINT
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL,ALLOCATABLE :: isExactFluxElem(:)     ! true if iElem is an element located within the ExactFlux region
LOGICAL,ALLOCATABLE :: isExactFluxFace(:)     ! true if iFace is a Face located wihtin or on the boarder (interface) of the
!                                             ! ExactFlux region
INTEGER,ALLOCATABLE :: ExactFluxToElem(:),ExactFluxToFace(:),ExactFluxInterToFace(:) ! mapping to total element/face list
INTEGER,ALLOCATABLE :: ElemToExactFlux(:),FaceToExactFlux(:),FaceToExactFluxInter(:) ! mapping to ExactFlux element/face list
REAL                :: InterFaceRegion(6)
INTEGER             :: nExactFluxElems,nExactFluxFaces,nExactFluxInterFaces
INTEGER             :: iElem,iSide,SideID,nExactFluxMasterInterFaces,sumExactFluxMasterInterFaces
!===================================================================================================================================
! get x,y, or z-position of interface
ExactFluxDir = GETINT('FluxDir','-1')  ! old name: remove in 2019
IF(ExactFluxDir.LE.0)THEN
  ExactFluxDir = GETINT('ExactFluxDir','3')
END IF
ExactFluxPosition    = GETREAL('ExactFluxPosition') ! initialize empty to force abort when values is not supplied
! set interface region, where one of the bounding box sides coinsides with the ExactFluxPosition in direction of ExactFluxDir
SELECT CASE(ABS(ExactFluxDir))
CASE(1) ! x
  InterFaceRegion(1:6)=(/-HUGE(1.),ExactFluxPosition,-HUGE(1.),HUGE(1.),-HUGE(1.),HUGE(1.)/)
CASE(2) ! y
  InterFaceRegion(1:6)=(/-HUGE(1.),HUGE(1.),-HUGE(1.),ExactFluxPosition,-HUGE(1.),HUGE(1.)/)
CASE(3) ! z
  InterFaceRegion(1:6)=(/-HUGE(1.),HUGE(1.),-HUGE(1.),HUGE(1.),-HUGE(1.),ExactFluxPosition/)
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,' Unknown exact flux direction: ExactFluxDir=',IntInfo=ExactFluxDir)
END SELECT

! set all elements lower/higher than the ExactFluxPosition to True/False for interface determination
CALL FindElementInRegion(isExactFluxElem,InterFaceRegion,ElementIsInside=.FALSE.,DoRadius=.FALSE.,Radius=-1.,DisplayInfo=.FALSE.)

! find all faces in the ExactFlux region
CALL FindInterfacesInRegion(isExactFluxFace,isExactFluxInterFace,isExactFluxElem,info_opt='find all faces in the ExactFlux region')

nExactFluxMasterInterFaces=0
DO iElem=1,nElems ! loop over all local elems
  DO iSide=1,6    ! loop over all local sides
    IF(ElemToSide(E2S_FLIP,iSide,iElem).EQ.0)THEN ! only master sides
      SideID=ElemToSide(E2S_SIDE_ID,iSide,iElem)
      IF(isExactFluxInterFace(SideID))THEN
        nExactFluxMasterInterFaces=nExactFluxMasterInterFaces+1
      END IF
    END IF
  END DO
END DO

#if USE_MPI
  sumExactFluxMasterInterFaces=0
  CALL MPI_REDUCE(nExactFluxMasterInterFaces , sumExactFluxMasterInterFaces , 1 , MPI_INTEGER, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
#else
  sumExactFluxMasterInterFaces=nExactFluxMasterInterFaces
#endif /*USE_MPI*/
LBWRITE(UNIT_StdOut,'(A8,I10,A)') '  Found ',sumExactFluxMasterInterFaces,' interfaces for ExactFlux.'

IF(mpiroot)THEN
  IF(sumExactFluxMasterInterFaces.LE.0)THEN
    CALL abort(__STAMP__&
        ,' [sumExactFluxMasterInterFaces.LE.0]: using ExactFlux but no interfaces found: sumExactFlux=',sumExactFluxMasterInterFaces)
  END IF
END IF


nExactFluxMasterInterFaces=0
DO iSide=1,lastMPISide_MINE ! nSides
  IF(SideToElem(S2E_ELEM_ID,iSide).EQ.-1) CYCLE
  IF(isExactFluxInterFace(SideID))THEN ! if an interface is encountered
    nExactFluxMasterInterFaces=nExactFluxMasterInterFaces+1
  END IF
END DO

#if USE_MPI
  sumExactFluxMasterInterFaces=0
  CALL MPI_REDUCE(nExactFluxMasterInterFaces , sumExactFluxMasterInterFaces , 1 , MPI_INTEGER, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
#else
  sumExactFluxMasterInterFaces=nExactFluxMasterInterFaces
#endif /*USE_MPI*/
LBWRITE(UNIT_StdOut,'(A8,I10,A)') '  Found ',sumExactFluxMasterInterFaces,' interfaces for ExactFlux. <<<<<< DEBUG this'







! Get number of ExactFlux Elems, Faces and Interfaces. Create Mappngs ExactFlux <-> physical region
CALL CountAndCreateMappings('ExactFlux',&
                            isExactFluxElem,isExactFluxFace,isExactFluxInterFace,&
                            nExactFluxElems,nExactFluxFaces, nExactFluxInterFaces,&
                            ElemToExactFlux,ExactFluxToElem,& ! these two are allocated
                            FaceToExactFlux,ExactFluxToFace,& ! these two are allocated
                            FaceToExactFluxInter,ExactFluxInterToFace) ! these two are allocated

! compute the outer radius of the mode in the cylindrical waveguide
CALL GetWaveGuideRadius(isExactFluxInterFace)

! Deallocate the vectors (must be deallocated because the used routine 'CountAndCreateMappings' requires INTENT,IN and ALLOCATABLE)
SDEALLOCATE(isExactFluxElem)
SDEALLOCATE(isExactFluxFace)
SDEALLOCATE(ExactFluxToElem)
SDEALLOCATE(ExactFluxToFace)
SDEALLOCATE(ExactFluxInterToFace)
SDEALLOCATE(ElemToExactFlux)
SDEALLOCATE(FaceToExactFlux)
SDEALLOCATE(FaceToExactFluxInter)
!CALL MPI_BARRIER(MPI_COMM_PICLAS, iError)
!stop
END SUBROUTINE InitExactFlux


SUBROUTINE FinalizeEquation()
!===================================================================================================================================
! Get the constant advection velocity vector from the ini file
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars,ONLY:EquationInitIsDone,isExactFluxInterFace
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
EquationInitIsDone = .FALSE.
SDEALLOCATE(isExactFluxInterFace)
END SUBROUTINE FinalizeEquation


PURE FUNCTION GAUSSIAN_PULSE(t, sigma)
!===================================================================================================================================
! create a gaussian pulse with mean=0, sigma=sigma and amplitude=1
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
REAL, INTENT(IN) :: t, sigma
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: GAUSSIAN_PULSE
!===================================================================================================================================

IF (ABS(t).GT.3.*sigma) THEN
  gaussian_pulse = 0
ELSE
  gaussian_pulse = EXP(-0.5*(t/sigma)**2)
END IF
RETURN
END FUNCTION GAUSSIAN_PULSE


PURE FUNCTION RECTANGULAR_PULSE(t, sigma)
!===================================================================================================================================
! create a rectangular pulse with mean=0, sigma=sigma and amplitude=1
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
REAL, INTENT(IN) :: t, sigma
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: RECTANGULAR_PULSE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: sigma_2
!===================================================================================================================================

sigma_2 = sigma * 0.5

IF ( (t.LT.-sigma_2) .OR. (t.GT.sigma_2) ) THEN
  rectangular_pulse = 0.
ELSE IF ( (t.EQ.-sigma_2) .OR. (t.EQ.sigma_2)) THEN
  rectangular_pulse = 0.5
ELSE
  rectangular_pulse = 1.0
END IF
RETURN

END FUNCTION RECTANGULAR_PULSE


PURE FUNCTION RECTANGULAR_PULSE_WITH_GAUSSIAN_EDGE(t, sigma, riseTime)
!===================================================================================================================================
! create a rectangular pulse with mean=0, sigma=sigma and amplitude=1
! the rising and falling edges of the pulse are gaussian shaped with a gaussian-distribution-sigma = riseTime/3
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
REAL, INTENT(in) :: t, sigma, riseTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: RECTANGULAR_PULSE_WITH_GAUSSIAN_EDGE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: sigma_2
!===================================================================================================================================

sigma_2 = sigma * 0.5

IF ((t.LT.-sigma_2-riseTime) .OR. (t.GT.sigma_2+riseTime)) THEN
  rectangular_pulse_with_gaussian_edge = 0.0
ELSE
  IF (t.LT.-sigma_2) THEN
    !rising edge
    rectangular_pulse_with_gaussian_edge = gaussian_pulse(t+sigma_2, riseTime/3.)
  ELSE IF (t.GT.sigma_2) THEN
    !falling edge
    rectangular_pulse_with_gaussian_edge = gaussian_pulse(t-sigma_2, riseTime/3.)
  ELSE
    !rectangular pulse
    rectangular_pulse_with_gaussian_edge = 1.0
  END IF
END IF
RETURN

END FUNCTION RECTANGULAR_PULSE_WITH_GAUSSIAN_EDGE


PPURE SUBROUTINE ExactFunc_TE_Circular_Waveguide(t,x,resu)
!===================================================================================================================================
! TE waves in a circular waveguide
! Book: Waveguide Handbook - Marcuvitz p. 69 ff.
! TEPolarization: r: circular polarization right handed
!                 l: circular polarization left handed
!                 x: linear polarization x-direction
!                 y: linear polarization y-direction
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars  ,ONLY: c,PI,c_inv,mu0
USE MOD_Equation_Vars ,ONLY: TEScale,TEPulse,TEFrequency,TEPolarization
USE MOD_Equation_Vars ,ONLY: TERadius,TEModeRoot,TEDelay,TEPulseSigma,TEPulseSeriesFrequence
USE MOD_Equation_Vars ,ONLY: TEPulseNumber,TEDirection,TEMode,TEPulseShape
USE MOD_Equation_Vars ,ONLY: ExactFluxPosition
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: t
REAL,INTENT(IN)                 :: x(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Resu(PP_nVar)    ! state in conservative variables
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: loopVar
INTEGER                         :: pG,mG
REAL                            :: omegaG
REAL                            :: k0, kz, kCut
REAL                            :: Z0, Zte
REAL                            :: Er,Br,Ephi,Bphi,Bz,Ez
REAL                            :: Bess_pG, dBess_pG
REAL                            :: cos_m_phi, sin_m_phi, cos_wt_kz, sin_wt_kz
REAL                            :: r,phi,z
REAL                            :: normFactor, V, I
REAL                            :: sigma_t,temporalWindow,tshift
!===================================================================================================================================

!free-space impedance
Z0 = mu0*c

!calculate polar-coordinates
r = SQRT(x(1)**2+x(2)**2)
! if a DOF is located in the origin, prevent division by zero ..
phi = ATAN2(x(2),x(1))
z = x(3) - ExactFluxPosition
omegaG = 2*PI*TEFrequency ! angular frequency

! TE_pG,mG
pG = TEMode(1) ! azimuthal wave number
mG = TEMode(2) ! radial wave number

!free-space wavenumber
k0 = omegaG*c_inv

!cutoff wavenumber
kCut = TEModeRoot/TERadius

!waveguide impedance
Zte = Z0 / SQRT(1.0 - (kCut/k0)**2)

!normalized phase-constant
kz = k0**2 - kCut**2

kz = SQRT(kz)
!kz = TEDirection * SQRT(kz)

! precompute bessel functions
Bess_pG  = BESSEL_JN(pG,r*kCut)
dBess_pG = 0.5*(BESSEL_JN(pG-1,r*kCut)-BESSEL_JN(pG+1,r*kCut))

! precompute the phase-terms
cos_m_phi = COS(REAL(pG)*phi)
sin_m_phi = SIN(REAL(pG)*phi)
cos_wt_kz = COS(kz*z - omegaG*t)
sin_wt_kz = SIN(kz*z - omegaG*t)
!cos_wt_kz = COS(-TEDirection*omegaG*t)
!sin_wt_kz = SIN(-TEDirection*omegaG*t)

! normalization factor
IF (ABS(pG).EQ.0)  THEN
  !special factor for TE_0_X modes
  normFactor = SQRT(1.0 / (PI * (TEModeRoot**2-pG**2))) / (BESSEL_JN(pG,TEModeRoot))
ELSE
  normFactor = SQRT(2.0 / (PI * (TEModeRoot**2-pG**2))) / (BESSEL_JN(pG,TEModeRoot))
END IF

! calculate the normalized amplitudes
V = normFactor * TEScale * SQRT(Zte)
I = normFactor * TEScale / SQRT(Zte)

! calculate the E and B field
IF (r.LT.1e-5) THEN !treat the singularity at r==0
  resu(1) = 0.0
  resu(3) = 0.0
  resu(5) = 0.0

  resu(2) = 0.5 * V * kCut
  resu(4) = mu0 * 0.5 * I * kCut
  resu(6) =  mu0 / k0 / Z0 * V * kCut**2 *  Bess_pG

  SELECT CASE(TEPolarization)
  CASE('x') ! linear polarized in x-direction
  resu(2) = resu(2)          * cos_wt_kz
  resu(4) = resu(4) * (-1.0) * cos_wt_kz
  resu(6) = resu(6) * (-1.0) * sin_wt_kz
  CASE('y') ! linear polarized in y-direction
  resu(2) = resu(2)          * cos_wt_kz
  resu(4) = resu(4) * (-1.0) * cos_wt_kz
  resu(6) = resu(6) * (-1.0) * sin_wt_kz
  CASE('r') ! left-handed circular polarized
  resu(2) = resu(2)          * (cos_wt_kz - sin_wt_kz)
  resu(4) = resu(4) * (-1.0) * (cos_wt_kz + sin_wt_kz)
  resu(6) = resu(6) * (-1.0) * (sin_wt_kz - cos_wt_kz)
  CASE('l') ! right-handed circular polarized
  resu(2) = resu(2)          * (cos_wt_kz + sin_wt_kz)
  resu(4) = resu(4) * (-1.0) * (cos_wt_kz - sin_wt_kz)
  resu(6) = resu(6) * (-1.0) * (sin_wt_kz + cos_wt_kz)
  END SELECT

  IF (ABS(pG).NE.1) THEN
    resu(1) = 0.0
    resu(2) = 0.0
    resu(4) = 0.0
    resu(5) = 0.0
  END IF
ELSE
  Er   = V * REAL(pG) / r *  Bess_pG
  Ephi = V * kCut         * dBess_pG
  Ez   =  0.

  Br   = mu0 * I * kCut              * dBess_pG
  Bphi = mu0 * I * REAL(pG) / r      *  Bess_pG
  Bz   = TEDirection * mu0 / k0 / Z0 * V * kCut**2 *  Bess_pG

  ! create the polarization
  SELECT CASE(TEPolarization)
  CASE('x') ! linear polarized in x-direction
    Er   = Er   * (-1.0) * cos_m_phi*cos_wt_kz
    Ephi = Ephi          * sin_m_phi*cos_wt_kz
    Ez   =  0.

    Br   = Br   * (-1.0) * sin_m_phi*cos_wt_kz
    Bphi = Bphi * (-1.0) * cos_m_phi*cos_wt_kz
    Bz   = Bz   * (-1.0) * sin_m_phi*sin_wt_kz
  CASE('y') ! linear polarized in y-direction
    Er   = Er   * sin_m_phi*cos_wt_kz
    Ephi = Ephi * cos_m_phi*cos_wt_kz
    Ez   =   0.

    Br   = Br   * (-1.0) * cos_m_phi*cos_wt_kz
    Bphi = Bphi          * sin_m_phi*cos_wt_kz
    Bz   = Bz   * (-1.0) * cos_m_phi*sin_wt_kz
  CASE('r') ! left-handed circular polarized
    Er   =   Er   * (sin_m_phi*cos_wt_kz - cos_m_phi*sin_wt_kz)
    Ephi =   Ephi * (cos_m_phi*cos_wt_kz + sin_m_phi*sin_wt_kz)
    Ez   =   0.

    Br   = Br   * (-1.0) * (cos_m_phi*cos_wt_kz + sin_m_phi*sin_wt_kz)
    Bphi = Bphi          * (sin_m_phi*cos_wt_kz - cos_m_phi*sin_wt_kz)
    Bz   = Bz   * (-1.0) * (cos_m_phi*sin_wt_kz - sin_m_phi*cos_wt_kz)
  CASE('l') ! right-handed circular polarized
    Er   = Er   * (sin_m_phi*cos_wt_kz + cos_m_phi*sin_wt_kz)
    Ephi = Ephi * (cos_m_phi*cos_wt_kz - sin_m_phi*sin_wt_kz)
    Ez   =   0.

    Br   = Br   * (-1.0) * (cos_m_phi*cos_wt_kz - sin_m_phi*sin_wt_kz)
    Bphi = Bphi          * (sin_m_phi*cos_wt_kz + cos_m_phi*sin_wt_kz)
    Bz   = Bz   * (-1.0) * (cos_m_phi*sin_wt_kz + sin_m_phi*cos_wt_kz)
  END SELECT

  !map to cartesian vectors
  resu(1) = (COS(phi)*Er - SIN(phi)*Ephi)
  resu(2) = (SIN(phi)*Er + COS(phi)*Ephi)
  resu(3) = Ez
  resu(4) = (COS(phi)*Br - SIN(phi)*Bphi)
  resu(5) = (SIN(phi)*Br + COS(phi)*Bphi)
  resu(6) = Bz
END IF

resu(7) = 0.0
resu(8) = 0.0

!input field as pulse
IF(TEPulse)THEN
  !sigma (standard deviation) of the gaussian pulse
  sigma_t = TEPulseSigma

  !generate series of pulses
  IF (TEPulseSeriesFrequence.GT.0) THEN
    !series of pulses with frequency TEPulseSeriesFrequence
    !add the amplitudes of all pulses
    temporalWindow = 0.
    loopVar = 0
    DO WHILE( loopVar.LT.TEPulseNumber )
      !time shift of the pulse: delay + pos of peak
      tShift = t - TEDelay - loopVar/TEPulseSeriesFrequence
      !temporalWindow = temporalWindow + EXP(-0.5*(tshift/sigma_t)**2)

      SELECT CASE(TRIM(TEPulseShape))
      CASE('gaussian')
        temporalWindow = temporalWindow + gaussian_pulse(tshift, sigma_t)
      CASE('rectangular')
        temporalWindow = temporalWindow + rectangular_pulse(tshift, sigma_t)
      CASE('rectangularGaussianEdges')
        temporalWindow = temporalWindow + rectangular_pulse_with_gaussian_edge(tshift, sigma_t,0.025e-9)
      CASE DEFAULT
      END SELECT

      loopVar = loopVar + 1
    END DO
  ELSE
    !single pulse
    !time shift of the pulse: delay
    tShift = t - TEDelay
    SELECT CASE(TRIM(TEPulseShape))
      CASE('gaussian')
        temporalWindow = gaussian_pulse(tshift, sigma_t)
      CASE('recangular')
        temporalWindow = rectangular_pulse(tshift, sigma_t)
      CASE('rectangularGaussianEdges')
        temporalWindow = rectangular_pulse_with_gaussian_edge(tshift, sigma_t,0.025e-9)
      CASE DEFAULT
      END SELECT
  END IF

  ! multiplicate the electric-field-distribution with the pulse
  ! because the Intensity should have the pulse-shape, the amplitudes have to be multiplicated with the sqrt of the pulse-shape
  resu(1:8) = resu(1:8) * SQRT(temporalWindow)
END IF

END SUBROUTINE ExactFunc_TE_Circular_Waveguide

END MODULE MOD_Equation

