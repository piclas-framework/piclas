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
#include "piclas.h"

MODULE MOD_LD_Vars
!===================================================================================================================================
! Contains the LD variables
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
  LOGICAL                           :: UseLD
  REAL  , ALLOCATABLE               :: LD_RHS(:,:)  ! RHS of the LD Method/ deltaV (npartmax, direction)
  REAL  , ALLOCATABLE               :: LD_DSMC_RHS(:,:)  ! RHS of the LD-DSMC Method/ deltaV (npartmax, direction)
  REAL                              :: LD_RepositionFak
  REAL                              :: LD_DSMC_RelaxationFak_BufferA
  REAL                              :: LD_RelaxationFak
  LOGICAL  , ALLOCATABLE           :: IsDoneLagVelo(:)  ! (nSides)
  REAL  , ALLOCATABLE              :: TempDens(:)
!!!  REAL  , ALLOCATABLE               :: NewNodePosIndx(:,:) ! (1:nDim,1:nNodes)  !!! nur für "Tetra-Methode"
  LOGICAL                           :: LD_CalcDelta_t
  LOGICAL                           :: LD_CalcResidual
  REAL, ALLOCATABLE                :: LD_Residual(:,:) ! Def. for LD Residual (number of Elements, 2nd index:
                                                                                                  ! 1.Velocity ux
                                                                                                  ! 2.Velocity uy
                                                                                                  ! 3.Velocity uz
                                                                                                  ! 4.Velocity |u|
                                                                                                  ! 5.Mass Density
                                                                                                  ! 6.Temperature
TYPE tLD_SecantMethod
  REAL                              :: Guess        ! 2nd guess, plus user defined value, [m/s], (default 10 m/s)
  REAL                              :: MaxIter      ! Max. number of iterations for LAGRANGIAN vell calculation
  REAL                              :: Accuracy     ! accuracy for LAGRANGIAN vell calculation
END TYPE
TYPE (tLD_SecantMethod)             :: LD_SecantMeth

TYPE tBulkValues                                                  ! LD bulk values
  REAL                              :: CellV(3)                   ! Velocity (vx,vy,vz)
  REAL                              :: DegreeOfFreedom            ! Average number of internal degree of freedom
  REAL                              :: Beta                       ! Thermal speed scale
  REAL                              :: MassDens                   ! Mass Density
  REAL                              :: BulkTemperature            ! BulkTemperature
  REAL                              :: DynamicVisc                ! dynamic viscousity
  REAL                              :: ThermalCond                ! thermal conuctivity
  REAL                              :: SpezGasConst               ! specific gas constant
  INTEGER                           :: CellType                   ! CellType (1:DSMC, 2: Bufferzone_A, 3: Bufferzone_B, 4: LD
END TYPE
TYPE(tBulkValues), ALLOCATABLE      :: BulkValues(:)              ! LD bulk values array (number of Elements)

TYPE tBulkValuesOpenBC                                                  ! LD bulk values
  REAL                              :: CellV(3)                   ! Velocity (vx,vy,vz)
  REAL                              :: DegreeOfFreedom            ! Average number of internal degree of freedom
  REAL                              :: Beta                       ! Thermal speed scale
  REAL                              :: MassDens                   ! Mass Density
  REAL                              :: DynamicVisc                ! dynamic viscousity
  REAL                              :: ThermalCond                ! thermal conuctivity
END TYPE
TYPE(tBulkValuesOpenBC), ALLOCATABLE :: BulkValuesOpenBC(:)       ! LD bulk values array for open BCs (number of Elements)

TYPE tSurfLagValues                                               ! LD Lagrangian cellsurface values
  REAL                              :: LagVelo                    ! Lagrangian Velocity
  REAL                              :: Area                       ! area of cellsurface
  REAL                              :: DeltaM(3)                  ! momentumflux (Mx,My,Mz)
  REAL                              :: DeltaE                     ! energyflux
  REAL                              :: LagNormVec(3)              ! corresponding face normal vector, [-]
  REAL                              :: LagTangVec(2,3)            ! cor. face tangential vectors t_1 & t_2, pos=second_indx
END TYPE
TYPE(tSurfLagValues), ALLOCATABLE  :: SurfLagValues(:,:,:)        ! LD Lagrangian cellsurface array (iLocSide, iElem, TriNum)

TYPE tMeanSurfValues                                              ! LD Lagrangian cellsurface values
  REAL                              :: MeanBaseD                  ! Mean Base (Base Vector in coordinat form of Mean Surface)
  REAL                              :: MeanBaseD2                 ! Mean Base for periodoc walls
  REAL                              :: MeanNormVec(3)             ! corresponding face normal vector, [-]
  REAL                              :: MeanLagVelo                ! Lagrangian Velocity for MeanSurf
  REAL                              :: DynamicVisc                ! dynamic viscousity
  REAL                              :: ThermalCond                ! thermal conuctivity
  REAL                              :: MeanBulkVelo(3)            ! mean bulk velocity for viscousity term
  REAL                              :: BulkVeloDiff(3)            ! difference in bulk velocity for viscousity term
  REAL                              :: BulkTempDiff               ! difference in temperature for viscousity term
  REAL                              :: CellCentDist(3)            ! vector difference between cell center for viscousity term
END TYPE
TYPE(tMeanSurfValues), ALLOCATABLE  :: MeanSurfValues(:,:)          ! Mean Surface for LD-Particle push (iLocSide, iElem)

REAL    , ALLOCATABLE               :: PartStateBulkValues(:,:)   ! LD particle values (npartmax, with 2nd index:
                                                                                                  ! 1.Velocity ux
                                                                                                  ! 2.Velocity uy
                                                                                                  ! 3.Velocity uz
                                                                                                  ! 4.Temperature
                                                                                                  ! 5.Degree of freedom
#ifdef MPI
  REAL  , ALLOCATABLE               :: MPINeighborBulkVal(:,:) ! LD values for cells on other procs (SideID, with 2nd index:
                                                                                                  ! 1.Velocity ux
                                                                                                  ! 2.Velocity uy
                                                                                                  ! 3.Velocity uz
                                                                                                  ! 4.Beta
                                                                                                  ! 5.Density
                                                                                                  ! 6.Dynamic Viscousity
                                                                                                  ! 7.Thermal Conuctivity
                                                                                                  ! 8.BulkTemperature
#endif

!===================================================================================================================================

END MODULE MOD_LD_Vars
