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
MODULE MOD_RadiationTrans_Vars
!===================================================================================================================================
! Contains the tadiation transport variables
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
TYPE tRadTrans
  INTEGER                                :: NumPhotonsPerCell
  INTEGER, ALLOCATABLE                   :: PhotPerCell(:)
  REAL                                   :: GlobalRadiationPower
  REAL                                   :: ScaledGlobalRadiationPower
  INTEGER                                :: GlobalPhotonNum
END TYPE

TYPE (tRadTrans)           :: RadTrans


TYPE tPhotonProps
  REAL                        :: PhotonPos(3)
  REAL                        :: PhotonLastPos(3)
  REAL                        :: PhotonDirection(3)
  REAL                        :: PhotonEnergy
  INTEGER                     :: ElemID
  INTEGER                     :: WaveLength
END TYPE

TYPE (tPhotonProps)           :: PhotonProps

REAL    , ALLOCATABLE :: RadiationElemEnergy(:,:)

REAL    , ALLOCATABLE :: Radiation_Emission_Spec_Total(:)
REAL    , ALLOCATABLE :: Radiation_Absorption_Spec_Total(:)

REAL    , ALLOCATABLE :: PhotonSampWall(:,:)

INTEGER               :: RadiationAbsorptionModel
INTEGER               :: RadiationDirectionModel
INTEGER               :: RadiationPhotonPosModel
LOGICAL               :: RadEmiAdaptPhotonNum
!===================================================================================================================================
END MODULE MOD_RadiationTrans_Vars
