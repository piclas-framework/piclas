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
PUBLIC :: DefineParametersSymmetry, Init_Symmetry
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
CALL prms%CreateIntOption(    'Particles-Symmetry-Order',  &
                              'Order of the Simulation 1, 2 or 3 D', '3')
CALL prms%CreateLogicalOption('Particles-Symmetry2D', 'Activating a 2D simulation on a mesh with one cell in z-direction in the '//&
                              'xy-plane (y ranging from 0 to the domain boundaries)', '.FALSE.')
CALL prms%CreateLogicalOption('Particles-Symmetry2DAxisymmetric', 'Activating an axisymmetric simulation with the same mesh '//&
                              'requirements as for the 2D case (y is then the radial direction)', '.FALSE.')

END SUBROUTINE DefineParametersSymmetry


SUBROUTINE Init_Symmetry()
!===================================================================================================================================
!> Initialize if a 2D/1D Simulation is performed and which type
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Symmetry_Vars    ,ONLY: Symmetry
USE MOD_DSMC_Vars        ,ONLY: RadialWeighting
USE MOD_ReadInTools      ,ONLY: GETLOGICAL,GETINT
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
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                 :: Symmetry2D
!===================================================================================================================================
Symmetry%Order = GETINT('Particles-Symmetry-Order')
Symmetry2D = GETLOGICAL('Particles-Symmetry2D')
IF(Symmetry2D.AND.(Symmetry%Order.EQ.3)) THEN
  Symmetry%Order = 2
  LBWRITE(*,*) 'WARNING: Particles-Symmetry-Order is set to 2 because of Particles-Symmetry2D=.TRUE. .'
  LBWRITE(*,*) 'Set Particles-Symmetry-Order=2 and remove Particles-Symmetry2D to avoid this warning'
ELSE IF(Symmetry2D) THEN
  CALL ABORT(__STAMP__&
    ,'ERROR: 2D Simulations either with Particles-Symmetry-Order=2 or (but not recommended) with Symmetry2D=.TRUE.')
END IF

IF((Symmetry%Order.LE.0).OR.(Symmetry%Order.GE.4)) CALL ABORT(__STAMP__&
,'Particles-Symmetry-Order (space dimension) has to be in the range of 1 to 3')

Symmetry%Axisymmetric = GETLOGICAL('Particles-Symmetry2DAxisymmetric')
IF(Symmetry%Axisymmetric.AND.(Symmetry%Order.EQ.3)) CALL ABORT(__STAMP__&
  ,'ERROR: Axisymmetric simulations only for 1D or 2D')
IF(Symmetry%Axisymmetric.AND.(Symmetry%Order.EQ.1))CALL ABORT(__STAMP__&
  ,'ERROR: Axisymmetric simulations are only implemented for Particles-Symmetry-Order=2 !')
IF(Symmetry%Axisymmetric) THEN
  RadialWeighting%DoRadialWeighting = GETLOGICAL('Particles-RadialWeighting')
ELSE
  RadialWeighting%DoRadialWeighting = .FALSE.
  RadialWeighting%PerformCloning = .FALSE.
END IF

END SUBROUTINE Init_Symmetry

END MODULE MOD_Symmetry