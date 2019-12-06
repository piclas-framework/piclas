!==================================================================================================================================
! Copyright (c) 2019 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

SUBROUTINE SuperB(mode)
!===================================================================================================================================
!> Routines for the calculation of magnetic fields of different permanent magnets and coils. Possibility to output the geometry
!> of the coil/magnet as a separate VTK file for visualization. Background field is stored in separate HDF5 file and can be utilized
!> as input for future simulations.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_SuperB_PermMag
USE MOD_SuperB_Coil
USE MOD_SuperB_Vars
USE MOD_Preproc               ,ONLY: PP_N
USE MOD_TimeDisc_Vars         ,ONLY: TEnd
USE MOD_Mesh_Vars             ,ONLY: nElems
USE MOD_Interpolation_Vars    ,ONLY: BGType, BGField, BGFieldAnalytic, BGFieldVTKOutput, BGDataSize, PsiMag
USE MOD_HDF5_Output_Tools     ,ONLY: WriteBGFieldToHDF5,WriteBGFieldAnalyticToHDF5
USE MOD_SuperB_Init           ,ONLY: InitializeSuperB
#ifdef PARTICLES
USE MOD_PICInterpolation_Vars ,ONLY: InterpolationType
USE MOD_Interpolation_Vars    ,ONLY: BGField_xGP, BGField_wBary
USE MOD_Interpolation_Vars    ,ONLY: xGP, wBary
#endif /*PARTICLES*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: mode ! 1: Standalone
                           ! 2: Called from PICLas
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE              :: BFieldPermMag(:,:,:,:,:)
INTEGER                       :: iMagnet, iCoil, iTimePoint
REAL                          :: timestep
!===================================================================================================================================
! Initialization of SuperB
CALL InitializeSuperB(mode)
! Setting the background field type, used in pic_interpolation.f90
BGType = 2
! Datasize not utilized so far but might be required
BGDataSize = 3

! Allocate and nullify the B-Field and the magnetic potential
ALLOCATE(BGField(       1:BGDataSize , 0:PP_N , 0:PP_N , 0:PP_N     , 1:nElems))
ALLOCATE(BFieldPermMag( 1:BGDataSize , 0:PP_N , 0:PP_N , 0:PP_N     , 1:nElems))
ALLOCATE(PsiMag(        0:PP_N       , 0:PP_N , 0:PP_N , 1:nElems))
ALLOCATE(MagnetFlag(    0:PP_N       , 0:PP_N , 0:PP_N , 1:nElems))
PsiMag     = 0
BGField    = 0
MagnetFlag = 0

IF(DoCalcErrorNormsSuperB)THEN
  ALLOCATE(BGFieldAnalytic(1:BGDataSize    , 0:PP_N , 0:PP_N , 0:PP_N     , 1:nElems))
  BGFieldAnalytic = 0
END IF ! DoCalcErrorNormsSuperB

#ifdef PARTICLES
IF(TRIM(InterpolationType).EQ.'particle_position') THEN
  BGField_xGP = xGP
  BGField_wBary = wBary
END IF
#endif /*PARTICLES*/

IF(NumOfPermanentMagnets.GT.0) THEN
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
END IF

IF (NumOfPermanentMagnets.GT.0) THEN
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') ' Calculation of the B-Field'
  CALL CalculateGradient()
  SWRITE(UNIT_stdOut,'(A)') '...Done!'
END IF

BFieldPermMag = BGField
BGField = 0

IF(ANY(TimeDepCoil)) THEN
  ALLOCATE(BGFieldTDep(1:3,0:PP_N,0:PP_N,0:PP_N,1:nElems,0:nTimePoints))
  timestep = tEnd / (nTimePoints-1)
  DO iTimePoint=0,(nTimePoints-1)
    SWRITE(UNIT_stdOut,'(A)')''
    SWRITE(UNIT_stdOut,'(A33,F10.5)')' Calculating time #', timestep*iTimePoint
    IF (NumOfCoils.GT.0) THEN
      SWRITE(UNIT_stdOut,'(132("-"))')
      SWRITE(UNIT_stdOUT,'(A)') ' Calculation of coils'
      DO iCoil=1,NumOfCoils
        SWRITE(UNIT_stdOut,'(A,I2)') ' Set up coil: ', iCoil
        SELECT CASE(TRIM(CoilInfo(iCoil)%Type))
        CASE('custom')
          CALL SetUpCoil(iCoil)
        CASE('circle')
          CALL SetUpCircleCoil(iCoil)
        CASE('rectangle')
          CALL SetUpRectangleCoil(iCoil)
        CASE('linear')
          CALL SetUpLinearConductor(iCoil)
        CASE DEFAULT
          CALL abort(&
          __STAMP__&
          ,'Unknown time-dependent coil type ['//TRIM(CoilInfo(iCoil)%Type)//']')
        END SELECT
        IF(BGFieldVTKOutput) THEN
          IF(iTimePoint.EQ.0) THEN
            SWRITE(UNIT_stdOut,'(A)') ' Write coil to VTK file'
            CALL WriteCoilVTK(iCoil)
          END IF
        END IF
        SWRITE(UNIT_stdOut,'(A)') ' Calculation of the B-Field'
        IF (TimeDepCoil(iCoil)) THEN
          CALL Jefimenko(iCoil, timestep * iTimePoint)
        ELSE
          CALL BiotSavart(iCoil)
        END IF
        SWRITE(UNIT_stdOut,'(A,I2)') '...Done coil #', iCoil
      END DO
    END IF
    BGField = BGField + BFieldPermMag
    BGFieldTDep(:,:,:,:,:,iTimePoint) = BGField(:,:,:,:,:)
    BGField = 0
  END DO
ELSE
  IF (NumOfCoils.GT.0) THEN
    SWRITE(UNIT_stdOut,'(132("-"))')
    SWRITE(UNIT_stdOUT,'(A,I2)') 'Calculation of the magnetic field from coils. Total number: ', NumOfCoils
    DO iCoil=1,NumOfCoils
      SWRITE(UNIT_stdOut,'(A,I2)',ADVANCE='NO') ' Set up coil: ', iCoil
      SELECT CASE(TRIM(CoilInfo(iCoil)%Type))
      CASE('custom')
        CALL SetUpCoil(iCoil)
      CASE('circle')
        CALL SetUpCircleCoil(iCoil)
      CASE('rectangle')
        CALL SetUpRectangleCoil(iCoil)
      CASE('linear')
        CALL SetUpLinearConductor(iCoil)
      CASE DEFAULT
        CALL abort(&
        __STAMP__&
        ,'Unknown coil type ['//TRIM(CoilInfo(iCoil)%Type)//']')
      END SELECT
      IF(BGFieldVTKOutput) THEN
        SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') '...Write coil to VTK file'
        CALL WriteCoilVTK(iCoil)
      END IF
      SWRITE(UNIT_stdOut,'(A)') '...Calculation of the B-Field'
      CALL BiotSavart(iCoil)
      SWRITE(UNIT_stdOut,'(A,I2)') '...Done coil #', iCoil
    END DO
  END IF
  BGField = BGField + BFieldPermMag
  CALL WriteBGFieldToHDF5()
  IF(DoCalcErrorNormsSuperB)THEN
    CALL WriteBGFieldAnalyticToHDF5()
  END IF ! DoCalcErrorNormsSuperB
END IF

SDEALLOCATE(PsiMag)
SDEALLOCATE(MagnetFlag)
SDEALLOCATE(PermanentMagnetInfo)
SDEALLOCATE(CoilInfo)
SDEALLOCATE(CurrentInfo)

END SUBROUTINE SuperB

END MODULE MOD_SuperB
