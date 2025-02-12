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
#if USE_PETSC
#include "petsc/finclude/petsc.h"
#endif
!===================================================================================================================================
!> Contains global variables used by the HDG modules.
!===================================================================================================================================
MODULE MOD_HDG_Vars_PETSc
! MODULES
#if USE_PETSC
USE PETSc
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
#if USE_HDG
#if USE_PETSC
Mat                 :: Smat_petsc
Vec                 :: RHS_petsc
Vec                 :: lambda_petsc
KSP                 :: ksp
Vec                 :: lambda_local_petsc
VecScatter          :: scatter_petsc
IS                  :: idx_local_petsc
IS                  :: idx_global_petsc
Vec                 :: lambda_local_conductors_petsc
VecScatter          :: scatter_conductors_petsc
IS                  :: idx_local_conductors_petsc
IS                  :: idx_global_conductors_petsc
INTEGER,ALLOCATABLE :: PETScGlobal(:)         !< PETScGlobal(SideID) maps the local SideID to global PETScSideID
INTEGER,ALLOCATABLE :: PETScLocalToSideID(:)  !< PETScLocalToSideID(PETScLocalSideID) maps the local PETSc side to SideID
REAL,ALLOCATABLE    :: Smat_BC(:,:,:,:)       !< side to side matrix for dirichlet (D) BCs, (ngpface,ngpface,6Sides,DSides)
INTEGER             :: nPETScSides            !< nSides - nDirichletSides
INTEGER             :: nPETScUniqueSides      !< nPETScSides - nMPISides_YOUR
INTEGER             :: nPETScUniqueSidesGlobal
INTEGER,ALLOCATABLE :: SmallMortarType(:,:)   !< Type of Mortar side ([1] Type, [2] Side, nSides)
                                              !< [1] Type: mortar type this small side belongs to (1-3)
                                              !< [2] Side: Small side number (1-4)
#endif /*USE_PETSC*/
#endif /*USE_HDG*/
END MODULE MOD_HDG_Vars_PETSc