!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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
MODULE MOD_PICInterpolation_Vars
!===================================================================================================================================
!> Variables for particle interpolation in PIC:
!> interpolation types, external fields (const. or variable), analytic interpolation
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE        :: FieldAtParticle(:,:)          !< 1st index: Ex,Ey,Ez,Bx,By,Bz
!                                                        !< 2nd index: PIC%maxParticleNumber
CHARACTER(LEN=256)      :: InterpolationType             !< Type of Interpolation-Method
LOGICAL                 :: InterpolationElemLoop         !< Interpolate with outer iElem-loop (not for many elements per processor!)
REAL                    :: externalField(6)              !< ext field is added to the maxwell-solver-field
LOGICAL                 :: DoInterpolation               !< Flag for interpolation
LOGICAL                 :: useBGField                    !< Flag for background field BGField via h5-File
LOGICAL                 :: CalcBField                    !< Calculate the background field BGField from parameters defined in the 
                                                         !< input file

CHARACTER(LEN=256)      :: FileNameVariableExternalField !< filename containing the external field csv table
LOGICAL                 :: useVariableExternalField      !< use given external field. only for Bz variation in z
REAL,ALLOCATABLE        :: VariableExternalField(:,:)    !< z - Pos , Bz
REAL                    :: DeltaExternalField            !< equidistant z-spacing for the VariableExternalField (fast computation)
INTEGER                 :: nIntPoints                    !< number of all interpolation points external field

#ifdef CODE_ANALYZE
LOGICAL                 :: DoInterpolationAnalytic       !< use analytic/algebraic functions for the field at the
!                                                        !< particle position
LOGICAL                 :: DoInitAnalyticalParticleState !< Calculate the initial velocity of the particle from an analytic expression

INTEGER                 :: AnalyticInterpolationType     !< Type of the analytic interpolation method
!                                                        !< 0: const. magnetostatic field: B = B_z = (/ 0 , 0 , 1 T /) = const.
!                                                        !< 1: magnetostatic field: B = B_z = (/ 0 , 0 , B_0 * EXP(x/l) /) = const.
!                                                        !<

INTEGER                 :: AnalyticInterpolationSubType  !< Sub-Type for the analytic interpolation method (in combination with
!                                                        !< AnalyticInterpolationType)

REAL                    :: AnalyticInterpolationP        !< parameter "p" for AnalyticInterpolationType = 1

REAL                    :: AnalyticInterpolationPhase    !< Phase shift angle phi that is used for cos(w*t + phi)

REAL                    :: L_2_Error_Part(1:6)           !< L2 error for the particle state
REAL                    :: L_2_Error_Part_time           !< old time for calculating the time step (when it is variable)
#endif /*CODE_ANALYZE*/
!===================================================================================================================================
END MODULE MOD_PICInterpolation_Vars
