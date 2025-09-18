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
#if USE_SUPER_B
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
USE MOD_Mesh_Vars             ,ONLY: nElems, offSetElem
USE MOD_Interpolation_Vars    ,ONLY: BGType,BGDataSize, N_BG, BGFieldVTKOutput
USE MOD_HDF5_Output_Fields    ,ONLY: WriteBGFieldToHDF5,WriteBGFieldAnalyticToHDF5
USE MOD_SuperB_Init           ,ONLY: InitializeSuperB
USE MOD_DG_Vars               ,ONLY: N_DG_Mapping
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iMagnet, iCoil, iTimePoint, Nloc, iElem
REAL              :: timestep
LOGICAL           :: BGFieldNotZero              !< flag is set true as soon as BGField is changed due to permanent magnet / coil
TYPE tBPermMag
  REAL,ALLOCATABLE :: Field(:,:,:,:)             !< [1:3,0:NBG,0:NBG,0:NBG,1:PP_nElems,1:nTimePoints]
END TYPE tBPermMag
TYPE(tBPermMag),ALLOCATABLE    :: BPermMag(:)
!===================================================================================================================================
! Initialization of SuperB
CALL InitializeSuperB()
! Setting the background field type, used in pic_interpolation.f90
BGType = 2
! Datasize not utilized so far but might be required
BGDataSize = 3

ALLOCATE(N_BG(1:nElems))
DO iElem = 1, nElems
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  ALLOCATE(N_BG(iElem)%BGField(1:BGDataSize,0:Nloc,0:Nloc,0:Nloc))
  N_BG(iElem)%BGField = 0.
  IF(UseTimeDepCoil) THEN
    ALLOCATE(N_BG(iElem)%BGFieldTDep(1:BGDataSize,0:Nloc,0:Nloc,0:Nloc,1:nTimePoints))
    N_BG(iElem)%BGFieldTDep = 0.
  END IF
  IF(DoCalcErrorNormsSuperB)THEN
    ALLOCATE(N_BG(iElem)%BGFieldAnalytic(1:BGDataSize,0:Nloc,0:Nloc,0:Nloc))
    N_BG(iElem)%BGFieldAnalytic = 0.
  END IF ! DoCalcErrorNormsSuperB
END DO
BGFieldNotZero = .FALSE.
! ------------------------------------------------------------
! Step 1: Calculate magnetic fields from permanent magnets
! ------------------------------------------------------------
IF(NumOfPermanentMagnets.GT.0) THEN
  BGFieldNotZero = .TRUE.
  ALLOCATE(BPermMag(1:nElems))
  ALLOCATE(PermanentMagnets(1:nElems))
  DO iElem = 1, nElems
    Nloc = N_DG_Mapping(2,iElem+offSetElem)
    ALLOCATE(BPermMag(iElem)%Field( 1:BGDataSize , 0:Nloc , 0:NLoc , 0:Nloc))
    ALLOCATE(N_BG(iElem)%PsiMag(0:Nloc , 0:NLoc , 0:Nloc))
    N_BG(iElem)%PsiMag = 0.
    ALLOCATE(PermanentMagnets(iElem)%Flag(0:Nloc , 0:NLoc , 0:Nloc))
    PermanentMagnets(iElem)%Flag = 0
  END DO
  SWRITE(UNIT_stdOut,'(132("-"))')
  SWRITE(UNIT_stdOUT,'(AI0)') ' CALCULATION OF MAGNETIC FIELD FROM PERMANENT MAGNETS - TOTAL NUMBER: ', NumOfPermanentMagnets
  DO iMagnet=1,NumOfPermanentMagnets
    SWRITE(UNIT_stdOUT,'(A,I0)',ADVANCE='NO') ' MAGNETIC POTENTIAL OF MAGNET ', iMagnet
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
    SWRITE(UNIT_stdOUT,'(A)') ' DONE!'
  END DO

  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') ' CALCULATION OF B-FIELD'
  CALL CalculateGradient() ! Changes BGField from magnets
  DO iElem = 1, nElems
    BPermMag(iElem)%Field = N_BG(iElem)%BGField
    N_BG(iElem)%BGField = 0.
  END DO
  SWRITE(UNIT_stdOut,'(A)') ' DONE!'

END IF

! ------------------------------------------------------------
! Step 2: Calculate magnetic fields from coils
! ------------------------------------------------------------
IF (NumOfCoils.GT.0) THEN
  SWRITE(UNIT_stdOut,'(132("-"))')
  SWRITE(UNIT_stdOUT,'(A,I0)') ' CALCULATION OF MAGNETIC FIELD FROM COILS - TOTAL NUMBER: ', NumOfCoils

  ! ------------------------------------------------------------
  ! Step 3: Calculate const. magnetic fields (not time-dependent)
  ! ------------------------------------------------------------
  DO iCoil=1,NumOfCoils
    IF(TimeDepCoil(iCoil)) CYCLE ! Skip time-dependent coils here
    BGFieldNotZero = .TRUE.
    CALL SetUpCoils(iCoil)
    CALL BiotSavart(iCoil) ! Changes BGField from coils
  END DO
END IF

! ------------------------------------------------------------
! Step 4: Add contribution of permanent magnets
! ------------------------------------------------------------
IF(NumOfPermanentMagnets.GT.0) THEN
  DO iElem = 1, nElems
    N_BG(iElem)%BGField = N_BG(iElem)%BGField + BPermMag(iElem)%Field
  END DO
END IF

! ------------------------------------------------------------
! Step 5: Calculate time-dependent magnetic fields
! ------------------------------------------------------------
IF(UseTimeDepCoil) THEN
  IF(DoCalcErrorNormsSuperB) CALL abort(__STAMP__,'DoCalcErrorNormsSuperB=T is not implemented for time-dependent fields')
  DO iCoil=1,NumOfCoils
    IF(.NOT.TimeDepCoil(iCoil)) CYCLE ! Skip time-constant coils here
    CALL SetUpCoils(iCoil)
    ASSOCIATE( f => CurrentInfo(iCoil)%CurrentFreq )
      BGFieldFrequency = f ! Set frequency for output to h5 file as attribute
      BGFieldCurrent   = CoilInfo(iCoil)%Current ! Set current maximum for output to h5 file as attribute
      IF (f.GT.0.) THEN; timestep = 1./(f*REAL(nTimePoints-1))
      ELSE             ; timestep = 0.; END IF
    END ASSOCIATE
    SWRITE(UNIT_stdOut,'(A,I0)',ADVANCE='NO') ' CALCULATION OF TIME-DEPENDENT B-FIELD - COIL ', iCoil
    DO iTimePoint = 1, nTimePoints
      CALL Jefimenko(iCoil, timestep, iTimePoint) ! Sets BGFieldTDep(:,:,:,:,:,iTimePoint)
      ! Add contribution of other coils or magnets
      IF(BGFieldNotZero) THEN
        DO iElem = 1, nElems
          N_BG(iElem)%BGFieldTDep(:,:,:,:,iTimePoint) = N_BG(iElem)%BGFieldTDep(:,:,:,:,iTimePoint) + N_BG(iElem)%BGField
        END DO
      END IF
    END DO ! iTimePoint = 1, nTimePoints
    SDEALLOCATE(CoilNodes)
    SWRITE(UNIT_stdOut,'(A)') ' DONE!'
  END DO

END IF

! Output linear conductors to a single VTK file
IF(BGFieldVTKOutput) CALL WriteLinearConductorVTK()

! ------------------------------------------------------------
! Step 6: Output to HDF5
! ------------------------------------------------------------
! Write BGField (time-constant) or WriteBGFieldToHDF5 (time-dependent) fields to h5
CALL WriteBGFieldToHDF5()
! Output analytic field solution if required
IF(DoCalcErrorNormsSuperB) CALL WriteBGFieldAnalyticToHDF5()

! Deallocate stuff
DO iElem = 1, nElems
  SDEALLOCATE(N_BG(iElem)%PsiMag)
END DO
SDEALLOCATE(PermanentMagnets)
SDEALLOCATE(PermanentMagnetInfo)
DO iCoil=1,NumOfCoils
  SDEALLOCATE(CoilInfo(iCoil)%CoilNodes)
END DO
SDEALLOCATE(CoilInfo)
SDEALLOCATE(CurrentInfo)

END SUBROUTINE SuperB
#endif
END MODULE MOD_SuperB
