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
REAL              :: c_corr
REAL              :: c_corr_c   !c_corr*c
REAL              :: c_corr_c2  !c_corr*c^2
REAL              :: eta_c      !(c_corr -1 )*c
REAL              :: fDamping, fDamping_pois
INTEGER           :: IniExactFunc
INTEGER           :: BCType(6)=-999
INTEGER           :: BoundaryCondition(6,2)
LOGICAL           :: EquationInitIsDone=.FALSE.
LOGICAL           :: DoParabolicDamping
REAL              :: xDipole(1:3) ! base point of electromagnetic dipole
INTEGER           :: alpha_shape
REAL              :: shapeFuncPrefix
REAL              :: rCutoff
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
                                                                                     'Phi_LMP'       , &
                                                                                     'Potential'     /)

REAL,SAVE,ALLOCATABLE                 :: D(:,:) !D Matrix Kopriva, derivative matrix
REAL,ALLOCATABLE                      :: E(:,:,:,:,:)
REAL,ALLOCATABLE                      :: Phi(:,:,:,:,:)
! DG time derivative
REAL,ALLOCATABLE                      :: Phit(:,:,:,:,:)
! number of array items in U, Ut, gradUx, gradUy, gradUz after allocated
INTEGER                               :: nTotalPhi
! interior face values for all elements
REAL,ALLOCATABLE                      :: Phi_master(:,:,:,:),Phi_slave(:,:,:,:)
REAL,ALLOCATABLE                      :: FluxPhi(:,:,:,:)

!===================================================================================================================================
END MODULE MOD_Equation_Vars
