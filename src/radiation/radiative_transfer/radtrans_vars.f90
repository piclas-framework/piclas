!==================================================================================================================================
! Copyright (c) 2018 - 2019 Marcel Pfeiffer
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

MODULE MOD_RadiationTrans_Vars
!===================================================================================================================================
! Contains the radiative transfer variables
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                                  :: useParticleRadiationSolver
INTEGER                                  :: RadObservationPointMethod
LOGICAL                                  :: ObservationDoConvolution

TYPE tRadObservationPoint
  REAL                                   :: StartPoint(3)
  REAL                                   :: Area
  REAL                                   :: AngularAperture
  REAL                                   :: ViewDirection(3)
  REAL                                   :: MidPoint(3)
  REAL                                   :: Diameter
  REAL                                   :: OrthoNormBasis(3,3)
  LOGICAL                                :: CalcFullSpectra
  REAL                                   :: SlitFunction(2)
  REAL                                   :: ShockTubeDiameter
END TYPE

TYPE(tRadObservationPoint)               :: RadObservationPoint
REAL,ALLOCATABLE                         :: RadObservation_Emission(:)
REAL,ALLOCATABLE                         :: RadObservation_Emission_Conv(:)  
INTEGER,ALLOCATABLE                      :: RadObservation_EmissionPart(:)  

TYPE tRadTrans
  INTEGER                                :: NumPhotonsPerCell
  REAL                                   :: GlobalRadiationPower
  REAL                                   :: ScaledGlobalRadiationPower
  INTEGER                                :: GlobalPhotonNum
END TYPE

TYPE (tRadTrans)           :: RadTrans

INTEGER               :: RadiationAbsorptionModel
INTEGER               :: RadiationDirectionModel
INTEGER               :: RadiationPhotonPosModel
INTEGER               :: RadiationPhotonWaveLengthModel
LOGICAL               :: RadEmiAdaptPhotonNum

REAL, ALLOCATABLE               :: RadiationElemAbsEnergy(:,:)
REAL, ALLOCATABLE               :: RadiationElemAbsEnergySpec(:,:)
REAL,ALLOCPOINT                 :: Radiation_Emission_Spec_Total(:)
REAL,ALLOCPOINT                 :: Radiation_Emission_Spec_Max(:)
INTEGER,ALLOCPOINT              :: RadTransPhotPerCell(:)     ! (WaveLen(:), number of mesh elements)
INTEGER, ALLOCATABLE            :: RadTransPhotPerCellLoc(:)
REAL, ALLOCPOINT                :: RadTransObsVolumeFrac(:)
REAL, ALLOCPOINT                :: RadObservationPOI(:,:)
#if USE_MPI
INTEGER                         :: RadTransPhotPerCell_Shared_Win
INTEGER,ALLOCPOINT              :: RadTransPhotPerCell_Shared(:)
INTEGER                         :: RadTransObsVolumeFrac_Shared_Win
REAL,ALLOCPOINT                 :: RadTransObsVolumeFrac_Shared(:)
INTEGER                         :: Radiation_Emission_Spec_Total_Shared_Win
REAL,ALLOCPOINT                 :: Radiation_Emission_Spec_Total_Shared(:)
INTEGER                         :: Radiation_Emission_Spec_Max_Shared_Win
REAL,ALLOCPOINT                 :: Radiation_Emission_Spec_Max_Shared(:)
INTEGER                         :: RadiationElemAbsEnergy_Shared_Win
REAL,POINTER                    :: RadiationElemAbsEnergy_Shared(:,:)
INTEGER                         :: RadiationElemAbsEnergySpec_Shared_Win
REAL,POINTER                    :: RadiationElemAbsEnergySpec_Shared(:,:)
INTEGER                         :: RadObservationPOI_Shared_Win
REAL,ALLOCPOINT                 :: RadObservationPOI_Shared(:,:)
#endif
!===================================================================================================================================
END MODULE MOD_RadiationTrans_Vars
