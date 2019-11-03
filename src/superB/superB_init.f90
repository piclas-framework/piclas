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

MODULE MOD_SuperB_Init
!===================================================================================================================================
!>
!===================================================================================================================================
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: InitializeSuperB, DefineParametersSuperB, FinalizeSuperB
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for SuperB
!==================================================================================================================================
SUBROUTINE DefineParametersSuperB()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection('SuperB')

CALL prms%CreateLogicalOption('PIC-CalcBField-OutputVTK', 'TO-DO','.FALSE.')

! Input of permanent magnets
CALL prms%CreateIntOption(      'NumOfPermanentMagnets'             , 'TO-DO','0')
CALL prms%CreateStringOption(   'PermanentMagnet[$]-Type'           , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('PermanentMagnet[$]-BasePoint'      , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'PermanentMagnet[$]-NumNodes'       , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('PermanentMagnet[$]-Magnetisation'  , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('PermanentMagnet[$]-BaseVector1'    , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('PermanentMagnet[$]-BaseVector2'    , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('PermanentMagnet[$]-BaseVector3'    , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PermanentMagnet[$]-Radius'         , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PermanentMagnet[$]-Radius2'        , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('PermanentMagnet[$]-HeightVector'   , 'TO-DO', numberedmulti=.TRUE.)

! Input of coils
CALL prms%CreateIntOption(      'NumOfCoils'            , 'TO-DO','0')
CALL prms%CreateStringOption(   'Coil[$]-Type'          , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Coil[$]-BasePoint'     , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Coil[$]-LengthVector'  , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Coil[$]-NumNodes'      , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Coil[$]-LoopNum'       , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Coil[$]-PointsPerLoop' , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Coil[$]-Current'       , 'TO-DO', numberedmulti=.TRUE.)
! Custom coils
CALL prms%CreateRealArrayOption('Coil[$]-AxisVec1'      , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Coil[$]-NumOfSegments' , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Coil[$]-Segment[$]-SegmentType'  , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Coil[$]-Segment[$]-NumOfPoints'  , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Coil[$]-Segment[$]-LineVector'   , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Coil[$]-Segment[$]-Radius'       , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Coil[$]-Segment[$]-Phi1'         , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Coil[$]-Segment[$]-Phi2'         , 'TO-DO', numberedmulti=.TRUE.)
! Circle coils
CALL prms%CreateRealOption(     'Coil[$]-Radius'        , 'TO-DO', numberedmulti=.TRUE.)
! Rectangle coils
CALL prms%CreateRealArrayOption('Coil[$]-RectVec1'      , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Coil[$]-RectVec2'      , 'TO-DO', numberedmulti=.TRUE.)
! Time-dependent coils
CALL prms%CreateLogicalOption(  'Coil[$]-TimeDepCoil'     , 'TO-DO','.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Coil[$]-CurrentAmplitude', 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Coil[$]-CurrentFrequency', 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Coil[$]-CurrentPhase'    , 'TO-DO', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'nTimePoints'             , 'TO-DO')

END SUBROUTINE DefineParametersSuperB


SUBROUTINE InitializeSuperB()
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_SuperB_Vars
USE MOD_Globals_Vars       ,ONLY: PI
USE MOD_Interpolation_Vars ,ONLY: BGFieldVTKOutput
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iMagnet, iCoil, iSegment
CHARACTER(LEN=32)         :: hilf,hilf2
!===================================================================================================================================

! Output of the magnets/coils as separate VTK files
BGFieldVTKOutput     = GETLOGICAL('PIC-CalcBField-OutputVTK','.FALSE.')

! Get the number of magnets
NumOfPermanentMagnets= GETINT('NumOfPermanentMagnets','0')
! Allocate the magnets
ALLOCATE(PermanentMagnetInfo(NumOfPermanentMagnets))
! Read-in of magnet parameters
IF (NumOfPermanentMagnets.GT.0) THEN
  DO iMagnet = 1,NumOfPermanentMagnets
    SWRITE(*,*) "|       Read-in infos of permanent magnet |", iMagnet
    WRITE(UNIT=hilf,FMT='(I0)') iMagnet
    PermanentMagnetInfo(iMagnet)%Type               = GETSTR('PermanentMagnet'//TRIM(hilf)//'-Type')
    PermanentMagnetInfo(iMagnet)%BasePoint(1:3)     = GETREALARRAY('PermanentMagnet'//TRIM(hilf)//'-BasePoint',3)
    PermanentMagnetInfo(iMagnet)%NumNodes           = GETINT('PermanentMagnet'//TRIM(hilf)//'-NumNodes')
    PermanentMagnetInfo(iMagnet)%Magnetisation(1:3) = GETREALARRAY('PermanentMagnet'//TRIM(hilf)//'-Magnetisation',3)
    SELECT CASE(TRIM(PermanentMagnetInfo(iMagnet)%Type))
    CASE('cuboid')
      PermanentMagnetInfo(iMagnet)%BaseVector1(1:3)   = GETREALARRAY('PermanentMagnet'//TRIM(hilf)//'-BaseVector1',3)
      PermanentMagnetInfo(iMagnet)%BaseVector2(1:3)   = GETREALARRAY('PermanentMagnet'//TRIM(hilf)//'-BaseVector2',3)
      PermanentMagnetInfo(iMagnet)%BaseVector3(1:3)   = GETREALARRAY('PermanentMagnet'//TRIM(hilf)//'-BaseVector3',3)
    CASE('sphere')
      PermanentMagnetInfo(iMagnet)%Radius             = GETREAL('PermanentMagnet'//TRIM(hilf)//'-Radius')
    CASE('cylinder')
      PermanentMagnetInfo(iMagnet)%Radius             = GETREAL('PermanentMagnet'//TRIM(hilf)//'-Radius')
      PermanentMagnetInfo(iMagnet)%HeightVector(1:3)  = GETREALARRAY('PermanentMagnet'//TRIM(hilf)//'-HeightVector',3)
    CASE('conic')
      PermanentMagnetInfo(iMagnet)%Radius             = GETREAL('PermanentMagnet'//TRIM(hilf)//'-Radius')
      PermanentMagnetInfo(iMagnet)%Radius2            = GETREAL('PermanentMagnet'//TRIM(hilf)//'-Radius2')
      PermanentMagnetInfo(iMagnet)%HeightVector(1:3)  = GETREALARRAY('PermanentMagnet'//TRIM(hilf)//'-HeightVector',3)
    END SELECT
  END DO
END IF

! Get the number of coils/conductors
NumOfCoils               = GETINT('NumOfCoils','0')
ALLOCATE(CoilInfo(NumOfCoils))
ALLOCATE(TimeDepCoil(NumOfCoils))
ALLOCATE(CurrentInfo(NumOfCoils))
! Read-in of coil/conductor parameters
IF (NumOfCoils.GT.0) THEN
  DO iCoil = 1,NumOfCoils
    SWRITE(*,*) "|       Read-in infos of coil |", iCoil
    WRITE(UNIT=hilf,FMT='(I0)') iCoil
    CoilInfo(iCoil)%Type              = GETSTR('Coil'//TRIM(hilf)//'-Type')
    CoilInfo(iCoil)%BasePoint(1:3)    = GETREALARRAY('Coil'//TRIM(hilf)//'-BasePoint',3)
    CoilInfo(iCoil)%LengthVector(1:3) = GETREALARRAY('Coil'//TRIM(hilf)//'-LengthVector',3)
    CoilInfo(iCoil)%Length = SQRT(CoilInfo(iCoil)%LengthVector(1)**2 + CoilInfo(iCoil)%LengthVector(2)**2 + &
                                  CoilInfo(iCoil)%LengthVector(3)**2)
    ! --------------------- Coil type ---------------------------------------
    SELECT CASE(TRIM(CoilInfo(iCoil)%Type))
    CASE('custom')
      CoilInfo(iCoil)%AxisVec1 = GETREALARRAY('Coil'//TRIM(hilf)//'-AxisVec1',3)
      IF (DOT_PRODUCT(CoilInfo(iCoil)%LengthVector,CoilInfo(iCoil)%AxisVec1).NE.0) THEN
        CALL abort(__STAMP__, &
        'ERROR in pic_interpolation.f90: Length vector and axis vector of coil need to be orthogonal!')
      END IF
      CoilInfo(iCoil)%LoopNum = GETINT('Coil'//TRIM(hilf)//'-LoopNum')
      CoilInfo(iCoil)%NumOfSegments = GETINT('Coil'//TRIM(hilf)//'-NumOfSegments')
      ALLOCATE(CoilInfo(iCoil)%SegmentInfo(CoilInfo(iCoil)%NumOfSegments))
      ! Start with 1 Loop Point as zero
      CoilInfo(iCoil)%PointsPerLoop = 1
      DO iSegment = 1,CoilInfo(iCoil)%NumOfSegments
        WRITE(UNIT=hilf2,FMT='(I0)') iSegment
        CoilInfo(iCoil)%SegmentInfo(iSegment)%SegmentType = GETSTR('Coil'//TRIM(hilf)//'-Segment'//TRIM(hilf2)//'-SegmentType')
        CoilInfo(iCoil)%SegmentInfo(iSegment)%NumOfPoints = GETINT('Coil'//TRIM(hilf)//'-Segment'//TRIM(hilf2)//'-NumOfPoints')
        ! Add the number of segment points to the total loop points
        ! Attention: Add the start/endpoint of two adjacent segments only once
        CoilInfo(iCoil)%PointsPerLoop = CoilInfo(iCoil)%PointsPerLoop + (CoilInfo(iCoil)%SegmentInfo(iSegment)%NumOfPoints - 1)
        SELECT CASE(TRIM(CoilInfo(iCoil)%SegmentInfo(iSegment)%SegmentType))
        CASE('line')
          CoilInfo(iCoil)%SegmentInfo(iSegment)%LineVector = GETREALARRAY('Coil'//TRIM(hilf)//&
                                                              '-Segment'//TRIM(hilf2)//'-LineVector',2)
        CASE('circle')
          CoilInfo(iCoil)%SegmentInfo(iSegment)%Radius = GETREAL('Coil'//TRIM(hilf)//'-Segment'//TRIM(hilf2)//'-Radius')
          CoilInfo(iCoil)%SegmentInfo(iSegment)%Phi1   = GETREAL('Coil'//TRIM(hilf)//'-Segment'//TRIM(hilf2)//'-Phi1')*PI/180.
          CoilInfo(iCoil)%SegmentInfo(iSegment)%Phi2   = GETREAL('Coil'//TRIM(hilf)//'-Segment'//TRIM(hilf2)//'-Phi2')*PI/180.
        CASE DEFAULT
          CALL abort(__STAMP__, &
            'No valid segment type defined! Must be either 1 (Line) or 2 (Circle segment)!')
        END SELECT
      END DO
      ! Multiply the points per loop with the number of loops in the coil
      ! Attention: Add the start/endpoint of two adjacent loops only once and don't forget the starting point
      CoilInfo(iCoil)%NumNodes = (CoilInfo(iCoil)%PointsPerLoop - 1) * CoilInfo(iCoil)%LoopNum + 1
    CASE('circle')
      CoilInfo(iCoil)%Radius            = GETREAL('Coil'//TRIM(hilf)//'-Radius')
      CoilInfo(iCoil)%LoopNum           = GETINT('Coil'//TRIM(hilf)//'-LoopNum')
      CoilInfo(iCoil)%PointsPerLoop     = GETINT('Coil'//TRIM(hilf)//'-PointsPerLoop')
      CoilInfo(iCoil)%NumNodes          = CoilInfo(iCoil)%LoopNum * CoilInfo(iCoil)%PointsPerLoop + 1
    CASE('rectangle')
      CoilInfo(iCoil)%AxisVec1          = GETREALARRAY('Coil'//TRIM(hilf)//'-AxisVec1',3)
      IF (DOT_PRODUCT(CoilInfo(iCoil)%LengthVector,CoilInfo(iCoil)%AxisVec1).NE.0) THEN
        CALL abort(__STAMP__, &
        'ERROR in pic_interpolation.f90: Length vector and axis vector of coil need to be orthogonal!')
      END IF
      CoilInfo(iCoil)%RectVec1          = GETREALARRAY('Coil'//TRIM(hilf)//'-RectVec1',2)
      CoilInfo(iCoil)%RectVec2          = GETREALARRAY('Coil'//TRIM(hilf)//'-RectVec2',2)
      CoilInfo(iCoil)%LoopNum           = GETINT('Coil'//TRIM(hilf)//'-LoopNum')
      CoilInfo(iCoil)%PointsPerLoop     = GETINT('Coil'//TRIM(hilf)//'-PointsPerLoop')
      IF (MOD(CoilInfo(iCoil)%PointsPerLoop - 1,4).NE.0) THEN
        CoilInfo(iCoil)%PointsPerLoop   = CoilInfo(iCoil)%PointsPerLoop + 4 - MOD(CoilInfo(iCoil)%PointsPerLoop - 1,4)
      END IF
      ! Multiply the points per loop with the number of loops in the coil
      ! Attention: Only add the start/endpoint of two adjacent loops only once and don't forget the starting point
      CoilInfo(iCoil)%NumNodes = (CoilInfo(iCoil)%PointsPerLoop - 1) * CoilInfo(iCoil)%LoopNum + 1
    CASE('linear')
      CoilInfo(iCoil)%NumNodes = GETINT('Coil'//TRIM(hilf)//'-NumNodes')
    END SELECT
    ! --------------------- Time-dependent current ---------------------------------------
    TimeDepCoil(iCoil) = GETLOGICAL('Coil'//TRIM(hilf)//'-TimeDepCoil')
    IF(TimeDepCoil(iCoil)) THEN
      CurrentInfo(iCoil)%CurrentAmpl = GETREAL('Coil'//TRIM(hilf)//'-CurrentAmplitude')
      CurrentInfo(iCoil)%CurrentFreq = GETREAL('Coil'//TRIM(hilf)//'-CurrentFrequency')
      CurrentInfo(iCoil)%CurrentPhase = GETREAL('Coil'//TRIM(hilf)//'-CurrentPhase')
    ELSE
      CoilInfo(iCoil)%Current = GETREAL('Coil'//TRIM(hilf)//'-Current')
    END IF
  END DO
END IF

! Discretisation in time
IF (ANY(TimeDepCoil)) THEN
  nTimePoints = GETINT('nTimePoints','0')
END IF

END SUBROUTINE InitializeSuperB


SUBROUTINE FinalizeSuperB()
!----------------------------------------------------------------------------------------------------------------------------------!
! Deallocate the respective arrays used by superB
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_SuperB_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(PsiMag)
SDEALLOCATE(MagnetFlag)
SDEALLOCATE(PermanentMagnetInfo)
SDEALLOCATE(CoilInfo)
SDEALLOCATE(TimeDepCoil)
SDEALLOCATE(CurrentInfo)
END SUBROUTINE FinalizeSuperB


END MODULE MOD_SuperB_Init
