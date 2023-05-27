!==================================================================================================================================
! Copyright (c) 2023 - 2023 Marcel Pfeiffer, Stephen Copplestone
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

MODULE MOD_RayTracing_Init
!===================================================================================================================================
! Initialization of Radiation Transport
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

PUBLIC::InitRayTracing, DefineParametersRayTracing
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for FP-Flow
!==================================================================================================================================
SUBROUTINE DefineParametersRayTracing()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Ray Tracing")

CALL prms%CreateIntOption(       'RayTracing-PartBound'      , 'TODO' , '0')
CALL prms%CreateLogicalOption(   'RayTracing-AdaptiveRays'   , 'TODO' , '.FALSE.')
CALL prms%CreateIntOption(       'RayTracing-NumRays'        , 'TODO' , '1')
CALL prms%CreateIntOption(       'RayTracing-RayPosModel'    , 'TODO' , '1')
CALL prms%CreateRealArrayOption( 'RayTracing-RayDirection'   , 'Direction vector for ray emission. Will be normalized after read-in.' , no=3)
CALL prms%CreateIntOption(       'RayTracing-PartBound'      , 'Particle boundary ID where rays are emitted from' , '0')
CALL prms%CreateRealOption(      'RayTracing-PulseDuration'  , 'Pulse duration tau for a Gaussian-type pulse with I~exp(-(t/tau)^2) [s]'                  )
CALL prms%CreateIntOption(       'RayTracing-NbrOfPulses'    , 'Number of pulses [-]'                                                                     , '1')
CALL prms%CreateRealOption(      'RayTracing-WaistRadius'    , 'Beam waist radius (in focal spot) w_b for Gaussian-type pulse with I~exp(-(r/w_b)^2) [m]' , '0.0')
CALL prms%CreateRealOption(      'RayTracing-WaveLength'     , 'Beam wavelength [m]'                                                                      )
CALL prms%CreateRealOption(      'RayTracing-RepetitionRate' , 'Pulse repetition rate (pulses per second) [Hz]'                                           )
CALL prms%CreateRealOption(      'RayTracing-Power'          , 'Average pulse power (energy of a single pulse times repetition rate) [W]'                 )
CALL prms%CreateLogicalOption(   'RayTracing-ForceAbsorption', 'Surface photon sampling is performed independent of the actual absorption/reflection outcome (default=T)', '.TRUE.')

END SUBROUTINE DefineParametersRayTracing


SUBROUTINE InitRayTracing()
!===================================================================================================================================
! Initialization of the radiation transport solver 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_RayTracing_Vars
USE MOD_Globals_Vars           ,ONLY: Pi
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_RadiationTrans_Vars    ,ONLY: RadiationAbsorptionModel,RadObservationPointMethod
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: factor,SurfaceNormal(3),alpha
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT RAY TRACING SOLVER ...'

! Do not absorb rays within the volume!
RadiationAbsorptionModel = 0
RadObservationPointMethod = 0

RayPartBound = GETINT('RayTracing-PartBound')
IF(RayPartBound.EQ.0) RETURN
IF(RayPartBound.LT.0) CALL CollectiveStop(__STAMP__,'RayTracing-PartBound must be > 0 to activate ray tracing on this boundary!')

! Get ray parameters
Ray%PulseDuration  = GETREAL('RayTracing-PulseDuration')
Ray%NbrOfPulses    = GETINT('RayTracing-NbrOfPulses')
Ray%tShift         = SQRT(8.0) * Ray%PulseDuration
Ray%WaistRadius    = GETREAL('RayTracing-WaistRadius')
Ray%WaveLength     = GETREAL('RayTracing-WaveLength')
Ray%RepetitionRate = GETREAL('RayTracing-RepetitionRate')
Ray%Period         = 1./Ray%RepetitionRate
Ray%Power          = GETREAL('RayTracing-Power')
Ray%Direction      = GETREALARRAY('RayTracing-RayDirection',3)
Ray%Direction      = UNITVECTOR(Ray%Direction)

AdaptiveRays       = GETLOGICAL('RayTracing-AdaptiveRays')
NumRays            = GETINT('RayTracing-NumRays')
RayPosModel        = GETINT('RayTracing-RayPosModel')
RayForceAbsorption = GETLOGICAL('RayTracing-ForceAbsorption')

ASSOCIATE( &
      E0      => Ray%Energy             ,&
      wb      => Ray%WaistRadius        ,&
      tau     => Ray%PulseDuration      ,&
      I0      => Ray%IntensityAmplitude ,&
      tShift  => Ray%tShift             ,&
      Period  => Ray%Period             ,&
      tActive => Ray%tActive            ,&
      A       => Ray%Area               )
  ! Derived quantities
  E0 = Ray%Power / Ray%RepetitionRate

  ! Rectangle
  ! Ray emission area
  A = (GEO%xmaxglob-GEO%xminglob) * (GEO%ymaxglob-GEO%yminglob)
  ! Normal vector of the ray emission area
  SurfaceNormal = (/ 0., 0., 1. /)
  ! Angle between emitted rays and emission area
  alpha = (90.-ABS(90.-(180./PI)*ACOS(DOT_PRODUCT(Ray%Direction,SurfaceNormal))))

  ! Calculate the peak intensity (uncorrected)
  I0 = E0 / (SQRT(PI)*tau*A)

  ! Correction factor due to temporal cut-off of the Gaussian pulse
  ! no need for correction in space because the function is not cut-off in space
  ! just consider the temporal cut-off for the rectangle
  factor = ERF(tShift/tau)
  factor = SQRT(PI)*tau*A
  I0 = E0 / factor

  ! Sanity check: overlapping of pulses is not implemented (use multiple emissions for this)
  IF(2.0*tShift.GT.Period) CALL abort(__STAMP__,'Pulse length (2*tShift) is greater than the pulse period. This is not implemented!')
  
  ! Active pulse time
  tActive = REAL(Ray%NbrOfPulses - 1)*Period + 2.0*tShift
END ASSOCIATE

CALL PrintOption('Rectangular ray emission area: A [m2]'                             , 'CALCUL.' , RealOpt=Ray%Area)
CALL PrintOption('Angle between emission area normal and ray direction: alpha [deg]' , 'CALCUL.' , RealOpt=alpha)
CALL PrintOption('Single pulse energy [J]'                                           , 'CALCUL.' , RealOpt=Ray%Energy)
CALL PrintOption('Intensity amplitude: I0 [W/m^2]'                                   , 'CALCUL.' , RealOpt=Ray%IntensityAmplitude)
CALL PrintOption('Corrected Intensity amplitude: I0_corr [W/m^2]'                    , 'CALCUL.' , RealOpt=Ray%IntensityAmplitude)
CALL PrintOption('Pulse period (Time between maximum of two pulses) [s]'             , 'CALCUL.' , RealOpt=Ray%Period)
CALL PrintOption('Temporal pulse width (pulse time 2x tShift) [s]'                   , 'CALCUL.' , RealOpt=2.0*Ray%tShift)
CALL PrintOption('Pulse will end at tActive (pulse final time) [s]'                  , 'CALCUL.' , RealOpt=Ray%tActive)

SWRITE(UNIT_stdOut,'(A)')' INIT RAY TRACING SOLVER DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitRayTracing


END MODULE MOD_RayTracing_Init
