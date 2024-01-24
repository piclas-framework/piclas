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

MODULE MOD_Particle_SurfaceFlux_Vars
!===================================================================================================================================
!> Variables and types for the surface flux, used directly in MOD_Particle_Vars as types are part of the Species type
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL               :: DoForceFreeSurfaceFlux                              ! switch if the stage reconstruction uses a force

TYPE tSurfFluxSubSideData
  REAL                                   :: projFak                          ! VeloVecIC projected to inwards normal
  REAL                                   :: a_nIn                            ! speed ratio projected to inwards normal
  REAL                                   :: Velo_t1                          ! Velo comp. of first orth. vector
  REAL                                   :: Velo_t2                          ! Velo comp. of second orth. vector
  REAL                                   :: nVFR                             ! normal volume flow rate through subside
  REAL                                   :: Dmax                             ! maximum Jacobian determinant of subside for opt. ARM
  REAL,ALLOCATABLE                       :: BezierControlPoints2D(:,:,:)     ! BCP of SubSide projected to VeloVecIC
                                                                             ! (1:2,0:NGeo,0:NGeo)
END TYPE tSurfFluxSubSideData

TYPE typeSurfaceflux
  INTEGER                                :: BC                               ! PartBound to be emitted from
  INTEGER                                :: Type                             ! Type set based on the read-in parameters
                                                                             ! 0: Standard surface flux
                                                                             ! 1: Adaptive surface flux (includes RadialWeighting)
                                                                             ! 2: Radial weighting in axisymmetric simulations
                                                                             ! 3: DoPoissonRounding .AND. .NOT.DoTimeDepInflow
                                                                             ! 4: DoTimeDepInflow
                                                                             ! 5: Thermionic emission with Schottky effect (requires HDG)
  CHARACTER(30)                          :: velocityDistribution             ! specifying keyword for velocity distribution
  REAL                                   :: VeloIC                           ! velocity for initial Data
  REAL                                   :: VeloVecIC(3)                     ! normalized velocity vector
  REAL                                   :: MWTemperatureIC                  ! Temperature for Maxwell Distribution
  REAL                                   :: PartDensity                      ! PartDensity (real particles per m^3)
  REAL                                   :: EmissionCurrent                  ! Current [A] (if defined replaces PartDensity)
  REAL                                   :: Massflow                         ! Mass flow [kg/s] (if defined replaces PartDensity)
  LOGICAL                                :: UseEmissionCurrent               ! Flag whether the emission current is used
  LOGICAL                                :: UseMassflow                      ! Flag whether the mass flow definition is used
  LOGICAL                                :: VeloIsNormal                     ! VeloIC is in Surf-Normal instead of VeloVecIC
  LOGICAL                                :: ReduceNoise                      ! reduce stat. noise by global calc. of PartIns
  LOGICAL                                :: AcceptReject                     ! perform ARM for skewness of RefMap-positioning
  INTEGER                                :: ARM_DmaxSampleN                  ! number of sample intervals in xi/eta for Dmax-calc.
  REAL                                   :: VFR_total                        ! Total Volumetric flow rate through surface
  REAL                     , ALLOCATABLE :: VFR_total_allProcs(:)            ! -''-, all values for root in ReduceNoise-case
  REAL                                   :: VFR_total_allProcsTotal          !     -''-, total
  REAL                                   :: totalAreaSF                      ! Total area of the respective surface flux
  INTEGER(KIND=8)                        :: InsertedParticle                 ! Number of all already inserted Particles
  INTEGER(KIND=8)                        :: InsertedParticleSurplus          ! accumulated "negative" number of inserted Particles
  TYPE(tSurfFluxSubSideData), ALLOCATABLE :: SurfFluxSubSideData(:,:,:)      ! SF-specific Data of Sides (1:N,1:N,1:SideNumber)
  LOGICAL                                :: CircularInflow                   ! Circular region, which can be used to define small
                                                                             ! geometry features on large boundaries
  INTEGER                                :: dir(3)                           ! axial (1) and orth. coordinates (2,3) of polar system
  REAL                                   :: origin(2)                        ! origin in orth. coordinates of polar system
  REAL                                   :: rmax                             ! max radius of to-be inserted particles
  REAL                                   :: rmin                             ! min radius of to-be inserted particles
  INTEGER, ALLOCATABLE                   :: SurfFluxSideRejectType(:)        ! Type if parts in side can be rejected (1:SideNumber)
  LOGICAL                                :: Adaptive                         ! Is the surface flux an adaptive boundary?
  INTEGER                                :: AdaptiveType                     ! Chose the adaptive type, description in DefineParams
  REAL                                   :: AdaptiveMassflow                 ! Mass flow [kg/s], which is held constant
  REAL                                   :: AdaptivePressure                 ! Static pressure [Pa], which is held constant
  REAL, ALLOCATABLE                      :: ConstMassflowWeight(:,:,:)       ! Adaptive, Type 4: Weighting factor for SF-sides to
                                                                             ! insert the right amount of particles
  REAL, ALLOCATABLE                      :: CircleAreaPerTriaSide(:,:,:)     ! Adaptive, Type 4: Area within a triangle, determined
                                                                             ! through Monte Carlo integration (initially)
  REAL                                   :: SampledMassflow                  ! Actual mass flow rate through a surface flux boundary
  REAL, ALLOCATABLE                      :: nVFRSub(:,:)                     ! normal volume flow rate through subsubside
  LOGICAL                                :: ThermionicEmission               ! Flag for thermionic emission
  LOGICAL                                :: SchottkyEffectTE                 ! Flag for Schottky effect in thermionic emission
  REAL                                   :: WorkFunctionTE                   ! Material-specific work function [Input: eV]
  REAL                                   :: RichardsonConstant               ! Material-specific constant [Input: A/(cm^2*K^2)]
END TYPE

LOGICAL                                 :: UseCircularInflow              ! Flag is set if the circular inflow feature is used:
                                                                          ! Particle insertion only in the defined circular area
                                                                          ! on the surface of a surface flux
INTEGER, ALLOCATABLE                    :: CountCircInflowType(:,:,:)     ! Counter whether cells are inside/partially inside or
                                                                          ! outside of circular region (only with CODE_ANALYZE)

!===================================================================================================================================
END MODULE MOD_Particle_SurfaceFlux_Vars
