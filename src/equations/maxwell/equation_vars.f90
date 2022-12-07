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
MODULE MOD_Equation_Vars
!===================================================================================================================================
! Contains the constant Advection Velocity Vector used for the linear scalar advection equation
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL           :: CentralFlux                             ! flag for central or upwind flux
REAL              :: c_corr
REAL              :: c_corr_c   !c_corr*c
REAL              :: c_corr_c2  !c_corr*c^2
REAL              :: eta_c      !(c_corr -1 )*c
REAL              :: fDamping
INTEGER           :: IniExactFunc
INTEGER           :: BCType(6)=-999
INTEGER           :: BoundaryCondition(6,2)
LOGICAL           :: EquationInitIsDone=.FALSE.
LOGICAL           :: DoParabolicDamping
REAL              :: DipoleOmega  ! electric dipole angular frequency
REAL              :: xDipole(1:3) ! base point of electromagnetic dipole
REAL              :: tPulse
INTEGER           :: alpha_shape
REAL              :: shapeFuncPrefix
REAL              :: rCutoff
! planar wave and gaussian beam
REAL              :: E_0                                    !> electric field amplitude maximum
REAL              :: E_0_vec(1:3)                           !> electric field unit vector of wave
REAL              :: BeamEta                                !> impedance: BeamEta=SQRT(mu0/eps0)
REAL              :: BeamWaveNumber                         !> wave number: BeamWaveNumber= 2*pi/WaveLength
REAL              :: BeamOmega                              !> angular frequency: BeamOmega  = BeamWaveNumber*c
INTEGER           :: BeamIdir1,BeamIdir2,BeamMainDir        !> wave beam direction aux variables
REAL              :: WaveVector(1:3)                        !> wave vector
REAL              :: WaveVectorE0dir(1:3)                   !> vector in which E_0_vec points (must be perpendicular to WaveVector)
LOGICAL           :: UseWaveVectorE0dir                     !> Use WaveVectorE0dir True/False (Default=.FALSE.)
REAL              :: WaveLength                             !> wave length
REAL,DIMENSION(3) :: WaveBasePoint                          !> wave base point || origin
REAL              :: tFWHM                                  !> time for full wave half maximum
REAL              :: BeamEnergy                             !> Total beam energy [J]
REAL              :: Beam_a0                                !> value to scale max. electric field
REAL              :: Beam_t0                                !> starting time of the (pulsed) electromagnetic wave. Default t0 = 0
REAL              :: tDelay                                 !> delay time filter for gaussian beam
REAL              :: I_0                                    !> max. intensity
REAL              :: sigma_t                                !> sigma_t can be used instead of tFWHM
REAL              :: Beam_w0,Beam_w0_2inv,sBeam_w0_2        !> spot size and inv of spot size (waist radius, smallest radius where
!                                                           !> the optical intensity drops to 1/e^2)
REAL              :: omega_0                                !> old variable name for Beam_w0 (remove in 2019)
REAL              :: tActive                                !> active time for laser pulse
REAL              :: TEScale                                !> scaling of input TE-wave strength
REAL              :: TEFrequency                            !> frequency of TE wave
REAL              :: TERadius                               !> Radius of Input TE wave, if wave is inserted over a plane
INTEGER           :: TEMode(1:2)                            !> input of TE_n,m mode
REAL              :: TEModeRoot                             !> root for the TEMode_n,m (root of derivative of Bessel function)
REAL              :: TEDelay                                !> Delay time for the TE wave (for CW and Pulse)
LOGICAL           :: TEPulse                                !> Flag for pulsed (Gaussian pulse) or continuous wave (CW)
REAL              :: TEPulseSigma                           !> standard deviation of the Gaussian pulse
REAL              :: TEPulseSeriesFrequence                 !> if TEPulseSeriesFrequence>0 and TEPulse=T -> a series of gaussian pulses with frequence TEPulseSeriesFrequence
INTEGER           :: TEPulseNumber                          !> number of generated pulses in a pulse-series
REAL              :: TEDirection                            !> +1 for propagation in +z direction, -1 for propagation in -z direction
CHARACTER(40)     :: TEPulseShape                           !> shape of the pulse: 'gaussian', 'rectangular', 'rectangularGaussianEdges'
!LOGICAL           :: TEPolarization                         !> linear or circular polarized (T=linear, F=circular)
CHARACTER(1)      :: TEPolarization                         !> polarization of the TE-mode: 'x' = linear in x-direction,
                                                            !                               'y' = linear in y direction
                                                            !                               'l' = left-handed circular
                                                            !                               'r' = right-handed circular
LOGICAL           :: DoExactFlux                            !> Flag to switch emission to flux superposition at certain positions
REAL              :: ExactFluxPosition                      !> x,y, or z-position of interface
LOGICAL,ALLOCATABLE::isExactFluxInterFace(:)                !> Flag for each side on which an exact flux is added
INTEGER           :: ExactFluxDir                           !> direction of flux for ExactFlux
! Boundary condition arrays
REAL,ALLOCATABLE     :: BCData(:,:,:,:)
INTEGER,ALLOCATABLE  :: nBCByType(:)
INTEGER,ALLOCATABLE  :: BCSideID(:,:)
! can specify BC state
CHARACTER(LEN=255):: BCStateFile

CHARACTER(LEN=255),DIMENSION(8),PARAMETER :: StrVarNames(8)=(/ CHARACTER(LEN=255) :: 'ElectricFieldX', &
                                                                                     'ElectricFieldY', &
                                                                                     'ElectricFieldZ', &
                                                                                     'MagneticFieldX', &
                                                                                     'MagneticFieldY', &
                                                                                     'MagneticFieldZ', &
                                                                                     'Phi'           , &
                                                                                     'Psi           ' /)
!===================================================================================================================================
END MODULE MOD_Equation_Vars
