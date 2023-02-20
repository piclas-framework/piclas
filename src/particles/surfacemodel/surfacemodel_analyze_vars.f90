!==================================================================================================================================
! Copyright (c) 2015 - 2019 Wladimir Reschke
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
MODULE MOD_SurfaceModel_Analyze_Vars
!===================================================================================================================================
! Contains global variables used by the Analyze modules.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                       :: SurfModelAnalyzeInitIsDone = .FALSE.
REAL                          :: SurfModelAnalyzeSampleTime  !< Accumulated simulation time between two outputs to SurfaceAnalyze.csv
INTEGER(KIND=8)               :: SurfaceAnalyzeStep       ! Analyze of surface is performed each Nth time step
! Output flags
LOGICAL                       :: CalcSurfCollCounter      ! Calculate the number of surface collision and number of
                                                          ! adsorbed particles per species
LOGICAL                       :: CalcPorousBCInfo         ! Calculate output for porous BCs (averaged over whole BC)

! Output variables
INTEGER,ALLOCATABLE           :: SurfAnalyzeCount(:)      ! Counter of surface collisions
INTEGER,ALLOCATABLE           :: SurfAnalyzeNumOfAds(:)   ! Number of adsorptions on surfaces
INTEGER,ALLOCATABLE           :: SurfAnalyzeNumOfDes(:)   ! Number of desorptions on surfaces

REAL,ALLOCATABLE              :: PorousBCOutput(:,:)  ! 1: Counter of impinged particles on the BC
                                                      ! 2: Measured pumping speed [m3/s] through # of deleted particles
                                                      ! 3: Pumping speed [m3/s] used to calculate the removal prob.
                                                      ! 4: Removal probability [0-1]
                                                      ! 5: Pressure at the BC normalized with the user-given pressure

! --- BoundaryParticleOutput = BPO
LOGICAL                       :: CalcBoundaryParticleOutput !< Flag for activating this output

TYPE tBoundaryParticleOutput
  REAL,ALLOCATABLE              :: RealPartOut(:,:)         !< Number of particles exiting on boundary X with species X

  INTEGER                       :: NPartBoundaries          !< Total number of boundaries where the particles are counted
  INTEGER,ALLOCATABLE           :: PartBoundaries(:)        !< Part-boundary number on which the particles are counted
                                                            !< Mapping iBPO to BCID: BPO%PartBoundaries(1:BPO%NPartBoundaries)  = 1:nPartBound
  INTEGER,ALLOCATABLE           :: FieldBoundaries(:)       !< Mapping iBPO to iBC : BPO%FieldBoundaries(1:BPO%NPartBoundaries) = 1:nBC
  INTEGER,ALLOCATABLE           :: BCIDToBPOBCID(:)         !< Mapping BCID to iBPO: BPO%BCIDToBPOBCID(1:nPartBound)            = 1:BPO%NPartBoundaries

  INTEGER                       :: NSpecies                 !< Total number of species which are considered for counting
  INTEGER,ALLOCATABLE           :: Species(:)               !< Species IDs which are considered for counting
  INTEGER,ALLOCATABLE           :: SpecIDToBPOSpecID(:)     !< Mapping SpecID to BPOSpecID (1:BpoNSpecies)

  LOGICAL                       :: OutputTotalElectricCurrent !< calculate the sum of all charged particle currents and SEE
END TYPE

TYPE(tBoundaryParticleOutput)   :: BPO

!-- Electron SEE emission counter
LOGICAL :: CalcElectronSEE !< Count the electron emission from BCs where SEE is active

TYPE tSEE
  REAL,ALLOCATABLE    :: RealElectronOut(:) !< Number of electrons emitted on boundary X

  INTEGER             :: NPartBoundaries    !< Total number of boundaries where the particles are counted
  INTEGER,ALLOCATABLE :: PartBoundaries(:)  !< Part-boundary number on which the particles are counted
  INTEGER,ALLOCATABLE :: BCIDToSEEBCID(:)   !< Mapping BCID to iSEE (1:nPartBound)
END TYPE

TYPE(tSEE)   :: SEE

!===================================================================================================================================
END MODULE MOD_SurfaceModel_Analyze_Vars
