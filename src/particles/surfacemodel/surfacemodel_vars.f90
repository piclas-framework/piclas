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
MODULE MOD_SurfaceModel_Vars
!===================================================================================================================================
!> Contains the SurfaceModel variables
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER , ALLOCATABLE            :: SurfModResultSpec(:,:)          ! Resulting species after surface model treatment
                                                                    ! (nPartBound,nSpecies)
CHARACTER(LEN=50) , ALLOCATABLE  :: SurfModEnergyDistribution(:)    ! Energy distribution of the reflected/created particle(s)
REAL , ALLOCATABLE               :: SurfModEmissionEnergy(:)        ! Energy of emitted particle for surface emission model (only available for SurfaceModel=7)
REAL , ALLOCATABLE               :: SurfModEmissionYield(:)         ! Emission yield factor for surface emission model (only changable for SurfaceModel=7)
REAL                             :: BackupVeloABS                   ! Backup of velocity during double-ARMfor 2nd SEE
! === Porous BC ====================================================================================================================
INTEGER                          :: nPorousBC                       ! Number of porous BCs
TYPE tPorousBC
  INTEGER                        :: BC                              ! Number of the reflective BC to be used as a porous BC
  REAL                           :: Pressure                        ! Pressure at the BC [Pa], user-given
  CHARACTER(LEN=50)              :: Type
  REAL                           :: PumpingSpeed                    ! Given/calculated pumping speed [m3/s]
  REAL                           :: DeltaPumpingSpeedKp             ! Proportional factor for the pumping speed controller
  REAL                           :: DeltaPumpingSpeedKi             ! Integral factor for the pumping speed controller
  CHARACTER(LEN=50)              :: Region                          ! Form of the porous BC: 'circular'
  LOGICAL                        :: UsingRegion                     ! Use only a smaller region on the BC as a porous BC (e.g. pump)
  INTEGER                        :: dir(3)                          ! axial (1) and orth. coordinates (2,3) of polar system
  REAL                           :: origin(2)                       ! origin in orth. coordinates of polar system
  REAL                           :: rmax                            ! max radius of to-be inserted particles
  REAL                           :: rmin                            ! min radius of to-be inserted particles
END TYPE
TYPE(tPorousBC), ALLOCATABLE     :: PorousBC(:)                     ! Container for the porous BC, allocated with nPorousBC

! === SEE BC ====================================================================================================================
REAL                             :: BulkElectronTempSEE             ! Bulk electron temperature for SEE model by Morozov2004
                                                                    ! read-in in Kelvin (when using the SEE mode), but is directly
                                                                    ! converted to eV for usage in the code
LOGICAL                          :: SurfModSEEelectronTempAutoamtic ! BulkElectronTempSEE = BulkElectronTemp, which is calculated
                                                                    ! automatically for the first species ID for electrons
REAL, ALLOCATABLE                :: SurfModSEEPowerFit(:,:)         ! Power-fit coefficients (1=a, 2=b) of the form: a*T(ev)^b

! === Sticking coefficient from simple models/interpolation
REAL, ALLOCATABLE                :: StickingCoefficientData(:,:)    ! Data for the model using non-bounce and condensation probability
                                                                    ! [:,1]: Maximum impact angle for model parameters
                                                                    ! [:,2]: Cut-off angle for non-bounce probability
                                                                    ! [:,3:4]: Temperature limits for linear interpolation of condensation probability
!===================================================================================================================================
END MODULE MOD_SurfaceModel_Vars
