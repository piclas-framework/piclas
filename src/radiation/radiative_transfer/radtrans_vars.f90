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

INTEGER               :: RadiationAbsorptionModel
INTEGER               :: RadiationDirectionModel
INTEGER               :: RadiationPhotonPosModel
INTEGER               :: RadiationPhotonWaveLengthModel
LOGICAL               :: RadEmiAdaptPhotonNum

REAL, ALLOCATABLE               :: RadiationElemAbsEnergy(:)
REAL,ALLOCPOINT                 :: Radiation_Emission_Spec_Total(:)
REAL,ALLOCPOINT                 :: Radiation_Emission_Spec_Max(:)
INTEGER,ALLOCPOINT              :: RadTransPhotPerCell(:)     ! (WaveLen(:), number of mesh elements)
INTEGER, ALLOCATABLE            :: RadTransPhotPerCellLoc(:)
REAL, ALLOCATABLE               :: PhotonSampWall(:,:)
#if USE_MPI
INTEGER                         :: RadTransPhotPerCell_Shared_Win
INTEGER,ALLOCPOINT              :: RadTransPhotPerCell_Shared(:)
INTEGER                         :: Radiation_Emission_Spec_Total_Shared_Win
REAL,ALLOCPOINT                 :: Radiation_Emission_Spec_Total_Shared(:)
INTEGER                         :: Radiation_Emission_Spec_Max_Shared_Win
REAL,ALLOCPOINT                 :: Radiation_Emission_Spec_Max_Shared(:)
INTEGER                         :: PhotonSampWall_Shared_Win
REAL,POINTER                    :: PhotonSampWall_Shared(:,:)
INTEGER                         :: RadiationElemAbsEnergy_Shared_Win
REAL,POINTER                    :: RadiationElemAbsEnergy_Shared(:)
#endif
!===================================================================================================================================
END MODULE MOD_RadiationTrans_Vars
