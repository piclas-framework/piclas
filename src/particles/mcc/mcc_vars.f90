!==================================================================================================================================
! Copyright (c) 2021 boltzplatz - numerical plasma dynamics GmbH
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
MODULE MOD_MCC_Vars
!===================================================================================================================================
! Contains the MCC variables
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

LOGICAL                             :: UseMCC               ! Flag (set automatically) to differentiate between MCC/XSec and regular DSMC
CHARACTER(LEN=256)                  :: XSec_Database        ! Name of the cross-section database
LOGICAL                             :: XSec_NullCollision   ! Flag (read-in) whether null collision method (determining number of pairs based on maximum relaxation frequency)
LOGICAL                             :: XSec_Relaxation      ! Flag (set automatically): usage of XSec data for the total relaxation probability

TYPE tXSecData
  REAL,ALLOCATABLE                  :: XSecData(:,:)        ! Cross-section as read-in from the database
                                                            ! 1: Energy (at read-in in [eV], during simulation in [J])
                                                            ! 2: Cross-section at the respective energy level [m^2]
  REAL                              :: Prob                 ! Event probability
END TYPE tXSecData

TYPE tSpeciesXSec
  LOGICAL                           :: UseCollXSec          ! Flag if the collisions of the species pair should be treated with
                                                            ! read-in collision cross-section (currently only with BGG)
  LOGICAL                           :: CollXSec_Effective   ! Flag whether the given cross-section data is "effective" (complete set
                                                            ! including other processes such as e.g.excitation and ionization) or
                                                            ! "elastic", including only the elastic collision cross-section.
  REAL                              :: CrossSection         ! Current collision cross-section
  REAL,ALLOCATABLE                  :: CollXSecData(:,:)    ! Collision cross-section as read-in from the database
                                                            ! 1: Energy (at read-in in [eV], during simulation in [J])
                                                            ! 2: Cross-section at the respective energy level [m^2]
  REAL,ALLOCATABLE                  :: VibXSecData(:,:)     ! Vibrational cross-section at the same intervals as the effective collision cross-section
                                                            ! 1: Energy (at read-in in [eV], during simulation in [J])
                                                            ! 2: Cross-section at the respective energy level [m^2]
  REAL                              :: ProbNull             ! Collision probability at the maximal collision frequency for the
                                                            ! null collision method of MCC
  LOGICAL                           :: UseVibXSec           ! Flag if cross-section data will be used for the vibrational relaxation
  TYPE(tXSecData),ALLOCATABLE       :: VibMode(:)           ! Vibrational cross-sections (nVib: Number of levels found in database)
  REAL                              :: VibProb              ! Relaxation probability
  REAL                              :: VibCount             ! Event counter
  INTEGER                           :: SpeciesToRelax       ! Save which species shall use the vibrational cross-sections
  TYPE(tXSecData),ALLOCATABLE       :: ReactionPath(:)      ! Reaction cross-sections (nPaths: Number of reactions for that case)
  LOGICAL                           :: UseElecXSec          ! Flag if cross-section data will be used for the electronic relaxation
  TYPE(tXSecData),ALLOCATABLE       :: ElecLevel(:)         ! Electronic cross-sections (nElec: Number of levels found in database)
END TYPE tSpeciesXSec

TYPE(tSpeciesXSec), ALLOCATABLE     :: SpecXSec(:)          ! Species cross-section related data (CollCase)

!===================================================================================================================================
END MODULE MOD_MCC_Vars
