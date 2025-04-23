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

MODULE MOD_Symmetry
!===================================================================================================================================
!> Routines for 2D (planar/axisymmetric) simulations
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
PUBLIC :: DefineParametersSymmetry, InitSymmetry
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for particles
!==================================================================================================================================
SUBROUTINE DefineParametersSymmetry()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE

CALL prms%SetSection("Particle Symmetry")
CALL prms%CreateIntOption(    'Particles-Symmetry-Order'        , 'Order of the Simulation 1, 2 or 3 D', '3')
CALL prms%CreateLogicalOption('Particles-Symmetry2DAxisymmetric', 'Activating an axisymmetric simulation with the same mesh requirements as for the 2D case (y is then the radial direction)', '.FALSE.')

END SUBROUTINE DefineParametersSymmetry


SUBROUTINE InitSymmetry()
!===================================================================================================================================
!> Initialize the dimension of the simulation (3D, Axisymmetric, 2D, 1D)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools           ,ONLY: GETLOGICAL,GETINT
USE MOD_Symmetry_Vars         ,ONLY: Symmetry
#if defined(PARTICLES)
USE MOD_Particle_Mesh_Tools   ,ONLY: InitParticleInsideQuad
USE MOD_Particle_TriaTracking ,ONLY: InitSingleParticleTriaTracking
USE MOD_Particle_InterSection ,ONLY: InitParticleThroughSideCheck1D2D
#endif /*defined(PARTICLES)*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

!===================================================================================================================================
Symmetry%Order = GETINT('Particles-Symmetry-Order')

IF((Symmetry%Order.LE.0).OR.(Symmetry%Order.GE.4)) CALL ABORT(__STAMP__&
  ,'Particles-Symmetry-Order (space dimension) has to be in the range of 1 to 3')

#if defined(PARTICLES)
! Initialize the function pointers for triatracking
CALL InitParticleInsideQuad()
CALL InitSingleParticleTriaTracking()
IF(Symmetry%Order.LE.2) CALL InitParticleThroughSideCheck1D2D()
#endif /*defined(PARTICLES)*/

Symmetry%Axisymmetric = GETLOGICAL('Particles-Symmetry2DAxisymmetric')
#if defined(PARTICLES)
! Only abort when particles are active
IF(Symmetry%Axisymmetric.AND.(Symmetry%Order.NE.2)) CALL ABORT(__STAMP__&
  ,'ERROR: Axisymmetric simulations are only implemented for Particles-Symmetry-Order=2 !')
#endif /*defined(PARTICLES)*/

END SUBROUTINE InitSymmetry

END MODULE MOD_Symmetry