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

CALL prms%CreateIntOption(      'RayTracing-PartBound'    , 'TODO' , '0')
CALL prms%CreateLogicalOption(  'RayTracing-AdaptiveRays' , 'TODO' , '.FALSE.')
CALL prms%CreateIntOption(      'RayTracing-NumRays'      , 'TODO' , '1')
CALL prms%CreateIntOption(      'RayTracing-RayPosModel'  , 'TODO' , '1')
CALL prms%CreateRealArrayOption('RayTracing-RayDirection' , 'Direction vector for ray emission. Will be normalized after read-in.' , no=3)
CALL prms%CreateIntOption(      'RayTracing-PartBound'    , 'Particle boundary ID where rays are emitted from' , '0')

CALL prms%CreateRealOption(     'RayTracing-PulseDuration'  , 'Pulse duration tau for a Gaussian-type pulse with I~exp(-(t/tau)^2) [s]'                  )
CALL prms%CreateRealOption(     'RayTracing-WaistRadius'    , 'Beam waist radius (in focal spot) w_b for Gaussian-type pulse with I~exp(-(r/w_b)^2) [m]' , '0.0')
CALL prms%CreateRealOption(     'RayTracing-WaveLength'     , 'Beam wavelength [m]'                                                                      )
CALL prms%CreateRealOption(     'RayTracing-RepetitionRate' , 'Pulse repetition rate (pulses per second) [Hz]'                                           )
CALL prms%CreateRealOption(     'RayTracing-Power'          , 'Average pulse power (energy of a single pulse times repetition rate) [W]'                 )

END SUBROUTINE DefineParametersRayTracing


SUBROUTINE InitRayTracing()
!===================================================================================================================================
! Initialization of the radiation transport solver 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_RayTracing_Vars
USE MOD_Mesh_Vars              ,ONLY: nGlobalElems
USE MOD_Globals_Vars           ,ONLY: Pi
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_Particle_Boundary_Vars ,ONLY: nComputeNodeSurfTotalSides
#if USE_MPI
USE MOD_MPI_Shared_Vars
USE MOD_MPI_Shared
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT RAY TRACING SOLVER ...'

RayPartBound = GETINT('RayTracing-PartBound')
IF(RayPartBound.EQ.0) RETURN
IF(RayPartBound.LT.0) CALL CollectiveStop(__STAMP__,'RayTracing-PartBound must be > 0!')

! Get ray parameters
Ray%PulseDuration    = GETREAL('RayTracing-PulseDuration')
Ray%WaistRadius      = GETREAL('RayTracing-WaistRadius')
Ray%WaveLength       = GETREAL('RayTracing-WaveLength')
Ray%RepetitionRate   = GETREAL('RayTracing-RepetitionRate')
Ray%Power            = GETREAL('RayTracing-Power')

ASSOCIATE( &
      E0   => Ray%Energy             ,&
      wb   => Ray%WaistRadius        ,&
      tau  => Ray%PulseDuration      ,&
      I0   => Ray%IntensityAmplitude ,&
      A    => Ray%Area               )
  ! Derived quantities
  E0 = Ray%Power / Ray%RepetitionRate
  ! Rectangle
  A  = (GEO%xmaxglob-GEO%xminglob) * (GEO%ymaxglob-GEO%xminglob)
  I0 = E0 / (SQRT(PI)*tau*A)
END ASSOCIATE

ALLOCATE(RayElemPassedEnergy(1:nGlobalElems))
RayElemPassedEnergy=0.0

AdaptiveRays = GETLOGICAL('RayTracing-AdaptiveRays')
NumRays      = GETINT('RayTracing-NumRays')
RayPosModel  = GETINT('RayTracing-RayPosModel')
RayDirection = GETREALARRAY('RayTracing-RayDirection',3)

#if USE_MPI
CALL Allocate_Shared((/nGlobalElems/),RayElemPassedEnergy_Shared_Win,RayElemPassedEnergy_Shared)
CALL MPI_WIN_LOCK_ALL(0,RayElemPassedEnergy_Shared_Win,IERROR)
IF (myComputeNodeRank.EQ.0) RayElemPassedEnergy_Shared = 0.
CALL BARRIER_AND_SYNC(RayElemPassedEnergy_Shared_Win,MPI_COMM_SHARED)  
#endif  /*USE_MPI*/

ALLOCATE(RaySampWall(2,1:nComputeNodeSurfTotalSides))
RaySampWall=0.0

#if USE_MPI
!> Then shared arrays for boundary sampling
CALL Allocate_Shared((/2,nComputeNodeSurfTotalSides/),RaySampWall_Shared_Win,RaySampWall_Shared)
CALL MPI_WIN_LOCK_ALL(0,RaySampWall_Shared_Win,IERROR)

IF (myComputeNodeRank.EQ.0) RaySampWall_Shared = 0.
CALL BARRIER_AND_SYNC(RaySampWall_Shared_Win,MPI_COMM_SHARED)
#endif

SWRITE(UNIT_stdOut,'(A)')' INIT RAY TRACING SOLVER DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitRayTracing


END MODULE MOD_RayTracing_Init
