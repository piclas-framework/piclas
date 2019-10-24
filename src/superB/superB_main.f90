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

SUBROUTINE SuperB()
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_SuperB_PermMag
USE MOD_SuperB_Coil
USE MOD_SuperB_Vars
USE MOD_Preproc               ,ONLY: PP_N
USE MOD_TimeDisc_Vars         ,ONLY: TEnd
USE MOD_Mesh_Vars             ,ONLY: nElems
USE MOD_PICInterpolation_Vars ,ONLY: InterpolationType, NBG, BGType, BGField, BGFieldVTKOutput
USE MOD_PICInterpolation_Vars ,ONLY: BGField_xGP, BGField_wBary, BGDataSize
USE MOD_Interpolation_Vars    ,ONLY: xGP, wBary
USE MOD_HDF5_Output_Tools     ,ONLY: WriteBFieldToHDF5
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE              :: BFieldPermMag(:,:,:,:,:)
INTEGER                       :: iMagnet, iCoil, iTimePoint
REAL                          :: timestep
!===================================================================================================================================

! Allocate and nullify the B-Field and the magnetic potential
ALLOCATE(BGField(1:3,0:PP_N,0:PP_N,0:PP_N,1:nElems))
ALLOCATE(BFieldPermMag(1:3,0:PP_N,0:PP_N,0:PP_N,1:nElems))
ALLOCATE(PsiMag(0:PP_N,0:PP_N,0:PP_N,1:nElems))
ALLOCATE(MagnetFlag(0:PP_N,0:PP_N,0:PP_N,1:nElems))
PsiMag = 0
BGField = 0
MagnetFlag = 0
! Setting the background field type, used in pic_interpolation.f90
BGType = 2
! Datasize not utilized so far but might be required for other interpolation types (e.g. particle_position)
BGDataSize = 3
! Background field order same as rest
NBG = PP_N

IF((TRIM(InterpolationType).NE.'particle_position').AND.(TRIM(InterpolationType).NE.'nearest_blurrycenter') &
    .AND.(TRIM(InterpolationType).NE.'nearest_gausspoint')) THEN
  CALL abort(__STAMP__,&
    'ERROR: Magnetic DSMC only supports particle_position, nearest_blurrycenter and nearest_gausspoint!')
END IF

IF(TRIM(InterpolationType).EQ.'particle_position') THEN
  ALLOCATE(BGField_xGP(0:NBG), BGField_wBary(0:NBG))
  BGField_xGP = xGP
  BGField_wBary = wBary
END IF

IF(NumOfCuboidMagnets.GT.0) THEN
  SWRITE(UNIT_stdOut,'(132("-"))')
  SWRITE(UNIT_stdOUT,'(A)') ' Calculation of the Magnetic Potential of Cuboid Magnets'
  DO iMagnet=1,NumOfCuboidMagnets
    SWRITE(UNIT_stdOUT,'(A,I2)') ' Magnet: ', iMagnet
    CALL CalculateCuboidMagneticPotential(iMagnet)
    SWRITE(UNIT_stdOUT,'(A,I2)') ' ... Done Magnet #', iMagnet
  ENDDO
ENDIF

IF(NumOfSphericMagnets.GT.0) THEN
  SWRITE(UNIT_stdOut,'(132("-"))')
  SWRITE(UNIT_stdOUT,'(A)') ' Calculation of the Magnetic Potential of Spheric Magnets'
  DO iMagnet=1,NumOfSphericMagnets
    SWRITE(UNIT_stdOUT,'(A,I2)') ' Magnet: ', iMagnet
    CALL CalculateSphericMagneticPotential(iMagnet)
    SWRITE(UNIT_stdOUT,'(A,I2)') ' ... Done Magnet #', iMagnet
  ENDDO
ENDIF

IF(NumOfCylindricMagnets.GT.0) THEN
  SWRITE(UNIT_stdOut,'(132("-"))')
  SWRITE(UNIT_stdOUT,'(A)') ' Calculation of the Magnetic Potential of Cylindric Magnets'
  DO iMagnet=1,NumOfCylindricMagnets
    SWRITE(UNIT_stdOUT,'(A,I2)') ' Magnet: ', iMagnet
    CALL CalculateCylindricMagneticPotential(iMagnet)
    SWRITE(UNIT_stdOUT,'(A,I2)') ' ... Done Magnet #', iMagnet
  ENDDO
ENDIF

IF(NumOfConicMagnets.GT.0) THEN
  SWRITE(UNIT_stdOut,'(132("-"))')
  SWRITE(UNIT_stdOUT,'(A)') ' Calculation of the Magnetic Potential of Conic Magnets'
  DO iMagnet=1,NumOfConicMagnets
    SWRITE(UNIT_stdOUT,'(A,I2)') ' Magnet: ', iMagnet
    CALL CalculateConicMagneticPotential(iMagnet)
    SWRITE(UNIT_stdOUT,'(A,I2)') ' ... Done Magnet #', iMagnet
  ENDDO
ENDIF

IF ((NumOfCuboidMagnets.GT.0).OR.(NumOfSphericMagnets.GT.0).OR.(NumOfCylindricMagnets.GT.0)&
                              .OR.(NumOfConicMagnets.GT.0)) THEN
  SWRITE(UNIT_stdOut,'(132("-"))')
  SWRITE(UNIT_stdOUT,'(A)') ' Calculate the Gradient of the Magnetic Potentail'
  CALL CalculateGradient()
  SWRITE(UNIT_stdOut,'(A)') ' Done Calculating B-Field'
ENDIF

BFieldPermMag = BGField
BGField = 0
! BGField(1,:,:,:,:) = PsiMag(:,:,:,:)

IF(ANY(TimeDepCoil)) THEN

  ALLOCATE(BGFieldTDep(1:3,0:PP_N,0:PP_N,0:PP_N,1:nElems,0:nTimePoints))
  timestep = tEnd / (nTimePoints-1)
  DO iTimePoint=0,(nTimePoints-1)
    SWRITE(UNIT_stdOut,'(A)')''
    SWRITE(UNIT_stdOut,'(A33,F10.5)')' Calculating Time #', timestep*iTimePoint

    IF (NumOfCoils.GT.0) THEN
      SWRITE(UNIT_stdOut,'(132("-"))')
      SWRITE(UNIT_stdOUT,'(A)') ' Calculation of Custom Coils' 
      DO iCoil=1,NumOfCoils
        SWRITE(UNIT_stdOut,'(A,I2)') ' Coil: ', iCoil
        SWRITE(UNIT_stdOut,'(A)') ' Set up coil'
        CALL SetUpCoil(iCoil)

        IF(BGFieldVTKOutput) THEN
          IF (iTimePoint.EQ.0) THEN
            SWRITE(UNIT_stdOut,'(A)') ' Write Coil to VTK File'
            CALL WriteCoilVTk(iCoil, 1)
          ENDIF
        END IF

        SWRITE(UNIT_stdOut,'(A)') ' Calculation of the B-Field'
        IF (TimeDepCoil(iCoil)) THEN
          CALL Jefimenko(iCoil, timestep * iTimePoint, 1)
        ELSE
          CALL BiotSavart(iCoil, 1)
        ENDIF
        
        CALL FinalizeCoil()
        SWRITE(UNIT_stdOut,'(A,I2)') '...Done Coil #', iCoil
      ENDDO
    ENDIF

    IF (NumOfCircleCoils.GT.0) THEN
      SWRITE(UNIT_stdOut,'(132("-"))')
      SWRITE(UNIT_stdOUT,'(A)') ' Calculation of Circle Coils' 
      DO iCoil=1,NumOfCircleCoils
        SWRITE(UNIT_stdOut,'(A,I2)') ' Circle Coil: ', iCoil
        SWRITE(UNIT_stdOut,'(A)') ' Set up coil'
        CALL SetUpCircleCoil(iCoil)

        IF(BGFieldVTKOutput) THEN
          IF (iTimePoint.EQ.0) THEN
            SWRITE(UNIT_stdOut,'(A)') ' Write Coil to VTK File'
            CALL WriteCoilVTk(iCoil, 2)
          ENDIF
        END IF

        SWRITE(UNIT_stdOut,'(A)') ' Calculation of the B-Field'
        IF (TimeDepCoil(iCoil)) THEN
          CALL Jefimenko(iCoil, timestep * iTimePoint, 2)
        ELSE
          CALL BiotSavart(iCoil, 2)
        ENDIF
        
        CALL FinalizeCoil()
        SWRITE(UNIT_stdOut,'(A,I2)') '...Done Coil #', iCoil
      ENDDO
    ENDIF

    IF (NumOfRectangleCoils.GT.0) THEN
      SWRITE(UNIT_stdOut,'(132("-"))')
      SWRITE(UNIT_stdOUT,'(A)') ' Calculation of Rectangle Coils' 
      DO iCoil=1,NumOfRectangleCoils
        SWRITE(UNIT_stdOut,'(A,I2)') ' Rectangle Coil: ', iCoil
        SWRITE(UNIT_stdOut,'(A)') ' Set up coil'
        CALL SetUpRectangleCoil(iCoil)

        IF(BGFieldVTKOutput) THEN
          IF (iTimePoint.EQ.0) THEN
            SWRITE(UNIT_stdOut,'(A)') ' Write Coil to VTK File'
            CALL WriteCoilVTk(iCoil, 3)
          ENDIF
        END IF

        SWRITE(UNIT_stdOut,'(A)') ' Calculation of the B-Field'
        IF (TimeDepCoil(iCoil)) THEN
          CALL Jefimenko(iCoil, timestep * iTimePoint, 3)
        ELSE
          CALL BiotSavart(iCoil, 3)
        ENDIF
        
        CALL FinalizeCoil()
        SWRITE(UNIT_stdOut,'(A,I2)') '...Done Coil #', iCoil
      ENDDO
    ENDIF

    IF (NumOfLinearConductors.GT.0) THEN
      SWRITE(UNIT_stdOut,'(132("-"))')
      SWRITE(UNIT_stdOUT,'(A)') ' Calculation of Linear Conductors' 
      DO iCoil=1,NumOfLinearConductors
        SWRITE(UNIT_stdOut,'(A,I2)') ' Linear Conductor: ', iCoil
        SWRITE(UNIT_stdOut,'(A)') ' Set up conductor'
        CALL SetUpLinearConductor(iCoil)

        IF(BGFieldVTKOutput) THEN
          IF (iTimePoint.EQ.0) THEN
            SWRITE(UNIT_stdOut,'(A)') ' Write Coil to VTK File'
            CALL WriteCoilVTk(iCoil, 4)
          ENDIF
        END IF

        SWRITE(UNIT_stdOut,'(A)') ' Calculation of the B-Field'
        IF (TimeDepCoil(iCoil)) THEN
          CALL Jefimenko(iCoil, timestep * iTimePoint, 4)
        ELSE
          CALL BiotSavart(iCoil, 4)
        ENDIF
        
        CALL FinalizeCoil()
        SWRITE(UNIT_stdOut,'(A,I2)') '...Done Conductor #', iCoil
      ENDDO
    ENDIF
    BGField = BGField + BFieldPermMag
    BGFieldTDep(:,:,:,:,:,iTimePoint) = BGField(:,:,:,:,:)
    BGField = 0
  ENDDO 
ELSE
  IF (NumOfCoils.GT.0) THEN
    SWRITE(UNIT_stdOut,'(132("-"))')
    SWRITE(UNIT_stdOUT,'(A)') ' Calculation of Custom Coils' 
    DO iCoil=1,NumOfCoils
      SWRITE(UNIT_stdOut,'(A,I2)') ' Coil: ', iCoil
      SWRITE(UNIT_stdOut,'(A)') ' Set up coil'
      CALL SetUpCoil(iCoil)

      IF(BGFieldVTKOutput) THEN
        SWRITE(UNIT_stdOut,'(A)') ' Write Coil to VTK File'
        CALL WriteCoilVTk(iCoil, 1)
      END IF
      
      SWRITE(UNIT_stdOut,'(A)') ' Calculation of the B-Field'
      CALL BiotSavart(iCoil, 1)
      
      CALL FinalizeCoil()
      SWRITE(UNIT_stdOut,'(A,I2)') '...Done Coil #', iCoil
    ENDDO
  ENDIF

  IF (NumOfCircleCoils.GT.0) THEN
    SWRITE(UNIT_stdOut,'(132("-"))')
    SWRITE(UNIT_stdOUT,'(A)') ' Calculation of Circle Coils' 
    DO iCoil=1,NumOfCircleCoils
      SWRITE(UNIT_stdOut,'(A,I2)') ' Circle Coil: ', iCoil
      SWRITE(UNIT_stdOut,'(A)') ' Set up coil'
      CALL SetUpCircleCoil(iCoil)

      IF(BGFieldVTKOutput) THEN
        SWRITE(UNIT_stdOut,'(A)') ' Write Coil to VTK File'
        CALL WriteCoilVTk(iCoil, 2)
      END IF

      SWRITE(UNIT_stdOut,'(A)') ' Calculation of the B-Field'
      CALL BiotSavart(iCoil, 2)
      
      CALL FinalizeCoil()
      SWRITE(UNIT_stdOut,'(A,I2)') '...Done Coil #', iCoil
    ENDDO
  ENDIF

  IF (NumOfRectangleCoils.GT.0) THEN
    SWRITE(UNIT_stdOut,'(132("-"))')
    SWRITE(UNIT_stdOUT,'(A)') ' Calculation of Rectangle Coils' 
    DO iCoil=1,NumOfRectangleCoils
      SWRITE(UNIT_stdOut,'(A,I2)') ' Rectangle Coil: ', iCoil
      SWRITE(UNIT_stdOut,'(A)') ' Set up coil'
      CALL SetUpRectangleCoil(iCoil)

      IF(BGFieldVTKOutput) THEN
        SWRITE(UNIT_stdOut,'(A)') ' Write Coil to VTK File'
        CALL WriteCoilVTk(iCoil, 3)
      END IF

      SWRITE(UNIT_stdOut,'(A)') ' Calculation of the B-Field'
      CALL BiotSavart(iCoil, 3)
      
      CALL FinalizeCoil()
      SWRITE(UNIT_stdOut,'(A,I2)') '...Done Coil #', iCoil
    ENDDO
  ENDIF

  IF (NumOfLinearConductors.GT.0) THEN
    SWRITE(UNIT_stdOut,'(132("-"))')
    SWRITE(UNIT_stdOUT,'(A)') ' Calculation of Linear Conductors' 
    DO iCoil=1,NumOfLinearConductors
      SWRITE(UNIT_stdOut,'(A,I2)') ' Linear Conductor: ', iCoil
      SWRITE(UNIT_stdOut,'(A)') ' Set up conductor'
      CALL SetUpLinearConductor(iCoil)

      IF(BGFieldVTKOutput) THEN
        SWRITE(UNIT_stdOut,'(A)') ' Write Coil to VTK File'
        CALL WriteCoilVTk(iCoil, 4)
      END IF
      SWRITE(UNIT_stdOut,'(A)') ' Calculation of the B-Field'
      CALL BiotSavart(iCoil, 4)
      
      CALL FinalizeCoil()
      SWRITE(UNIT_stdOut,'(A,I2)') '...Done Conductor #', iCoil
    ENDDO
  ENDIF

  BGField = BGField + BFieldPermMag

  CALL WriteBFieldToHDF5(0.)
ENDIF

END SUBROUTINE SuperB

END MODULE MOD_SuperB