!==================================================================================================================================
! Copyright (c) 2019 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_SuperB
!===================================================================================================================================
!>
!===================================================================================================================================
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: SuperB
!===================================================================================================================================

!===================================================================================================================================

CONTAINS

SUBROUTINE SuperB()
!===================================================================================================================================
!> Routines for the calculation of magnetic fields of different permanent magnets and coils. Possibility to output the geometry
!> of the coil/magnet as a separate VTK file for visualization. Background field is stored in separate HDF5 file and can be utilized
!> as input for future simulations.
!>
!> Workflow:
!>   Step 1: Calculate magnetic fields from permanent magnets
!>   Step 2: Calculate magnetic fields from coils
!>   Step 3: Calculate const. magnetic fields (not time-dependent)
!>   Step 4: Add contribution of permanent magnets
!>   Step 5: Calculate time-dependent magnetic fields
!>   Step 6: Output to HDF5
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_SuperB_PermMag
USE MOD_SuperB_Coil
USE MOD_SuperB_Vars
USE MOD_Preproc
USE MOD_Mesh_Vars             ,ONLY: nElems
USE MOD_Interpolation_Vars    ,ONLY: BGType, BGField, BGFieldAnalytic, BGDataSize, PsiMag
USE MOD_HDF5_Output_Fields    ,ONLY: WriteBGFieldToHDF5,WriteBGFieldAnalyticToHDF5
USE MOD_SuperB_Init           ,ONLY: InitializeSuperB
USE MOD_Globals_Vars          ,ONLY: TimeStampLenStr,TimeStampLenStr2
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE  :: BFieldPermMag(:,:,:,:,:)
INTEGER           :: iMagnet, iCoil, iTimePoint
REAL              :: timestep
INTEGER,PARAMETER :: TimeStampLength=21
LOGICAL           :: BGFieldNotZero              !< flag is set true as soon as BGField is changed due to permanent magnet / coil
!===================================================================================================================================
WRITE(UNIT=TimeStampLenStr ,FMT='(I0)') TimeStampLength
WRITE(UNIT=TimeStampLenStr2,FMT='(I0)') TimeStampLength-4
! Initialization of SuperB
CALL InitializeSuperB()
! Setting the background field type, used in pic_interpolation.f90
BGType = 2
! Datasize not utilized so far but might be required
BGDataSize = 3

! Allocate and nullify the B-Field and the magnetic potential
ALLOCATE(BGField(       1:BGDataSize , 0:PP_N , 0:PP_N , 0:PP_N     , 1:nElems))
BGField    = 0.

IF(DoCalcErrorNormsSuperB)THEN
  ALLOCATE(BGFieldAnalytic(1:BGDataSize    , 0:PP_N , 0:PP_N , 0:PP_N     , 1:nElems))
  BGFieldAnalytic = 0.
END IF ! DoCalcErrorNormsSuperB

BGFieldNotZero = .FALSE.

! ------------------------------------------------------------
! Step 1: Calculate magnetic fields from permanent magnets
! ------------------------------------------------------------
IF(NumOfPermanentMagnets.GT.0) THEN
  BGFieldNotZero = .TRUE.
  ALLOCATE(BFieldPermMag( 1:BGDataSize , 0:PP_N , 0:PP_N , 0:PP_N     , 1:nElems))
  ALLOCATE(MagnetFlag(    0:PP_N       , 0:PP_N , 0:PP_N , 1:nElems))
  ALLOCATE(PsiMag(        0:PP_N       , 0:PP_N , 0:PP_N , 1:nElems))
  MagnetFlag = 0
  PsiMag     = 0
  SWRITE(UNIT_stdOut,'(132("-"))')
  SWRITE(UNIT_stdOUT,'(A)') ' Calculation of the magnetic field from permanent magnets.'
  DO iMagnet=1,NumOfPermanentMagnets
    SWRITE(UNIT_stdOUT,'(A,I2)',ADVANCE='NO') ' Magnet: ', iMagnet
    SELECT CASE(TRIM(PermanentMagnetInfo(iMagnet)%Type))
    CASE('cuboid')
      CALL CalculateCuboidMagneticPotential(iMagnet)
    CASE('sphere')
      CALL CalculateSphericMagneticPotential(iMagnet)
    CASE('cylinder')
      CALL CalculateCylindricMagneticPotential(iMagnet)
    CASE('conic')
      CALL CalculateConicMagneticPotential(iMagnet)
    END SELECT
    SWRITE(UNIT_stdOUT,'(A,I2,A)') '...Magnetic potential of magnet ', iMagnet, ' done!'
  END DO

  SWRITE(UNIT_stdOut,'(A)') ' Calculation of the B-Field'
  CALL CalculateGradient() ! Changes BGField from magnets
  BFieldPermMag = BGField  ! Store in temporary variable
  BGField = 0.
  SWRITE(UNIT_stdOut,'(A)') ' ...Done!'

END IF

! ------------------------------------------------------------
! Step 2: Calculate magnetic fields from coils
! ------------------------------------------------------------
IF (NumOfCoils.GT.0) THEN
  SWRITE(UNIT_stdOut,'(132("-"))')
  SWRITE(UNIT_stdOUT,'(A,I2)') 'Calculation of magnetic field from coils. Total number: ', NumOfCoils

  ! ------------------------------------------------------------
  ! Step 3: Calculate const. magnetic fields (not time-dependent)
  ! ------------------------------------------------------------
  DO iCoil=1,NumOfCoils
    IF(TimeDepCoil(iCoil)) CYCLE ! Skip time-dependent coils here
    BGFieldNotZero = .TRUE.
    CALL SetUpCoils(iCoil)
    CALL BiotSavart(iCoil) ! Changes BGField from coils
    SWRITE(UNIT_stdOut,'(A,I0,A,I0)') '...Done coil #', iCoil," of ",NumOfCoils
  END DO
END IF

! ------------------------------------------------------------
! Step 4: Add contribution of permanent magnets
! ------------------------------------------------------------
IF(NumOfPermanentMagnets.GT.0) BGField = BGField + BFieldPermMag

! ------------------------------------------------------------
! Step 5: Calculate time-dependent magnetic fields
! ------------------------------------------------------------
IF(UseTimeDepCoil) THEN
  IF(DoCalcErrorNormsSuperB) CALL abort(__STAMP__,'DoCalcErrorNormsSuperB=T is not implemented for time-dependent fields')
  ALLOCATE(BGFieldTDep(1:BGDataSize,0:PP_N,0:PP_N,0:PP_N,1:nElems,1:nTimePoints))
  BGFieldTDep = 0.

  DO iCoil=1,NumOfCoils
    IF(.NOT.TimeDepCoil(iCoil)) CYCLE ! Skip time-constant coils here
    CALL SetUpCoils(iCoil)
    ASSOCIATE( f => CurrentInfo(iCoil)%CurrentFreq )
      BGFieldFrequency = f ! Set frequency for output to h5 file as attribute
      BGFieldCurrent   = CoilInfo(iCoil)%Current ! Set current maximum for output to h5 file as attribute
      IF (f.GT.0.) THEN; timestep = 1./(f*REAL(nTimePoints-1))
      ELSE             ; timestep = 0.; END IF
    END ASSOCIATE
    SWRITE(UNIT_stdOut,'(A)') '...Calculation of the B-Field'
    DO iTimePoint = 1, nTimePoints
      CALL Jefimenko(iCoil, timestep, iTimePoint) ! Sets BGFieldTDep(:,:,:,:,:,iTimePoint)
      ! Add contribution of other coils or magnets
      IF(BGFieldNotZero) BGFieldTDep(:,:,:,:,:,iTimePoint) = BGFieldTDep(:,:,:,:,:,iTimePoint) + BGField(:,:,:,:,:)
    END DO ! iTimePoint = 1, nTimePoints
    SDEALLOCATE(CoilNodes)
    SWRITE(UNIT_stdOut,'(A,I0,A,I0)') '...Done coil #', iCoil," of ",NumOfCoils
  END DO

END IF

! ------------------------------------------------------------
! Step 6: Output to HDF5
! ------------------------------------------------------------
! Write BGField (time-constant) or WriteBGFieldToHDF5 (time-dependent) fields to h5
CALL WriteBGFieldToHDF5()
! Output analytic field solution if required
IF(DoCalcErrorNormsSuperB) CALL WriteBGFieldAnalyticToHDF5()

! Deallocate stuff
SDEALLOCATE(PsiMag)
SDEALLOCATE(MagnetFlag)
SDEALLOCATE(PermanentMagnetInfo)
SDEALLOCATE(CoilInfo)
SDEALLOCATE(CurrentInfo)

END SUBROUTINE SuperB

END MODULE MOD_SuperB
