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
MODULE MOD_Particle_Sampling_Vars
!===================================================================================================================================
! Contains the variables for sampling of macroscopic properties from the particle information (general for all modules: DSMC, BGK, FP)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! General sampling within the simulation domain

!-----------------------------------------------------------------------------------------------------------------------------------
! Sampling of elements with a boundary for adaptive surface flux and porous BC
LOGICAL                                 :: UseAdaptiveBC                  ! Flag is set if an adaptive boundary is present
LOGICAL                                 :: AdaptBCAverageValBC            ! Flag to enable/disable averaging accross the whole BC
REAL, ALLOCATABLE                       :: AdaptBCAverageMacroVal(:,:,:)  ! Macroscopic values averaged over BC
                                                                          ! (1:3, 1:nSpecies, 1:nSurfaceFluxBCs)
                                                                          !  1: Number density
                                                                          !  2: Translational temperature [K]
                                                                          !  3: Statis pressure [Pa]
REAL                                    :: AdaptBCRelaxFactor             ! weighting factor theta for weighting of average
                                                                          ! instantaneous values with those
                                                                          ! of previous iterations
INTEGER                                 :: AdaptBCSampIter                ! Number of sampling iterations as given by user
INTEGER                                 :: AdaptBCSampIterReadIn          ! Number of sampling iterations as read-in from state
LOGICAL                                 :: AdaptBCTruncAverage            ! Flag to enable/disable a truncated running average
INTEGER                                 :: AdaptBCSampleElemNum           ! Number of elements with an adaptive BC
INTEGER, ALLOCATABLE                    :: AdaptBCMapSampleToElem(:)      ! Mapping from sample element ID to the local element ID
INTEGER, ALLOCATABLE                    :: AdaptBCMapElemToSample(:)      ! Mapping from local element ID to the sample element ID
REAL, ALLOCATABLE                       :: AdaptBCAverage(:,:,:,:)        ! Truncated running average (current value replaces the first)
REAL, ALLOCATABLE                       :: AdaptBCSample(:,:,:)           ! Particle sample near boundaries
REAL, ALLOCATABLE                       :: AdaptBCMacroVal(:,:,:)         ! Macroscopic value near boundaries
                                                                          ! (1:7,1:AdaptBCSampleElemNum,1:nSpecies)
                                                                          !  1:  VELOX
                                                                          !  2:  VELOY
                                                                          !  3:  VELOZ
                                                                          !  4:  NUMBER DENSITY
                                                                          !  5:  Pumping capacity [m3/s]
                                                                          !  6:  Static pressure [Pa]
                                                                          !  7:  Integral pressure difference [Pa]
REAL, ALLOCATABLE                       :: AdaptiveData(:,:)              ! Macroscopic value near boundaries (for output)
REAL, ALLOCATABLE                       :: AdaptBCAreaSurfaceFlux(:,:)    ! UseCircularInflow: Surflux area as the sum of actual elements
REAL, ALLOCATABLE                       :: AdaptBCVolSurfaceFlux(:,:)     ! AdaptiveBC: Volume of adjacent cells
REAL, ALLOCATABLE                       :: AdaptBCBackupVelocity(:,:,:)   ! Velocity is stored as backup for iterations without particles
                                                                          ! in the cell [1:3,1:AdaptBCSampleElemNum,1:nSpecies]
INTEGER, ALLOCATABLE                    :: AdaptBCPartNumOut(:,:)         ! SurfaceFlux, Type 4: Number of particles exiting through
                                                                          ! the adaptive boundary condition
INTEGER                                 :: offSetElemAdaptBCSample
INTEGER                                 :: AdaptBCSampleElemNumGlobal
!-----------------------------------------------------------------------------------------------------------------------------------
END MODULE MOD_Particle_Sampling_Vars
