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

CALL prms%CreateLogicalOption('DoCalcErrorNormsSuperB', 'Set true to compute L2 and LInf error norms for magnetic fields.','.FALSE.')
CALL prms%CreateIntOption(    'NAnalyze'              , 'Polynomial degree at which analysis is performed (e.g. for L2 errors).\n'//&
                                                        'Default: 2*N.')
CALL prms%CreateLogicalOption('PIC-CalcBField-OutputVTK', 'Output of the magnets/coils geometry as separate VTK files','.FALSE.')

! Input of permanent magnets
CALL prms%SetSection('Input of permanent magnets')
CALL prms%CreateIntOption(      'NumOfPermanentMagnets'             , 'Number of permanent magnets','0')
CALL prms%CreateStringOption(   'PermanentMagnet[$]-Type'           , 'Permanent magnet type: cuboid, sphere, cylinder, conic', &
                                                                      numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('PermanentMagnet[$]-BasePoint'      , 'Origin (vector) for geometry parametrization', &
                                                                      numberedmulti=.TRUE., no=3)
CALL prms%CreateIntOption(      'PermanentMagnet[$]-NumNodes'       , 'Number of Gauss points for the discretization of the '//&
                                                                      'permanent magnet:\n'//&
                                                                      'Cuboid: N points in each direction (total number: 6N^2)\n'//&
                                                                      'Sphere: N divisions in the zenith direction with 2*N '//&
                                                                      'points in the azimuthal direction\n'//&
                                                                      'Cylinder: N divisions along height vector, 2*N points in '//&
                                                                      'the azimuthal direction, N points in radial direction on '//&
                                                                      'the top and bottom face\n'//&
                                                                      'Conic: see the cylinder NumNodes description', &
                                                                      numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('PermanentMagnet[$]-Magnetisation'  , 'Magnetisation vector in [A/m]', numberedmulti=.TRUE., no=3)
CALL prms%CreateRealArrayOption('PermanentMagnet[$]-BaseVector1'    , 'Vector 1 spanning the cuboid', numberedmulti=.TRUE., no=3)
CALL prms%CreateRealArrayOption('PermanentMagnet[$]-BaseVector2'    , 'Vector 2 spanning the cuboid', numberedmulti=.TRUE., no=3)
CALL prms%CreateRealArrayOption('PermanentMagnet[$]-BaseVector3'    , 'Vector 3 spanning the cuboid', numberedmulti=.TRUE., no=3)
CALL prms%CreateRealOption(     'PermanentMagnet[$]-Radius'         , 'Radius of a spheric, cylindric and conic (first radius) '//&
                                                                      'permanent magnet', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PermanentMagnet[$]-Radius2'        , 'Radius of the second radius of the conic permanent magnet'//&
                                                                      ' or inner radius for hollow cylinders', &
                                                                      numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('PermanentMagnet[$]-HeightVector'   , 'Height vector of cylindric and conic permanent magnet', &
                                                                      numberedmulti=.TRUE., no=3)

! Input of coils
CALL prms%SetSection('Input of coils')
CALL prms%CreateIntOption(      'NumOfCoils'            , 'Number of coils','0')
CALL prms%CreateStringOption(   'Coil[$]-Type'          , 'Coil type: custom, circle, rectangular, linear conductor (straight '//&
                                                          'wire)', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Coil[$]-BasePoint'     , 'Origin vector of the coil/linear conductor', numberedmulti=.TRUE., no=3)
CALL prms%CreateRealArrayOption('Coil[$]-LengthVector'  , 'Length vector of the coil/linear conductor, normal to the cross-'//&
                                                          'sectional area', numberedmulti=.TRUE., no=3)
CALL prms%CreateRealOption(     'Coil[$]-Current'       , 'Electrical coil current [A]', numberedmulti=.TRUE.)

! Linear conductor (calculated from the number of loops and points per loop for coils)
CALL prms%CreateIntOption(      'Coil[$]-NumNodes'      , 'Number of nodes for a linear conductor' &
                                                        , numberedmulti=.TRUE.)
! Coils
CALL prms%CreateIntOption(      'Coil[$]-LoopNum'       , 'Number of coil loops', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Coil[$]-PointsPerLoop' , 'Number of points per loop (azimuthal discretization)', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Coil[$]-AxisVec1'      , 'Axial vector defines the orientation of the cross-section together '//&
                                                          'with the length vector', numberedmulti=.TRUE., no=3)
! Custom coils
CALL prms%SetSection('Custom coils')
CALL prms%CreateIntOption(      'Coil[$]-NumOfSegments' , 'Number of segments for the custom coil definition', numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Coil[$]-Segment[$]-SegmentType'  , 'Possible segment types: line or circle', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Coil[$]-Segment[$]-NumOfPoints'  , 'Number of points to discretize the segment', &
                                                                    numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Coil[$]-Segment[$]-LineVector'   , 'Line segment: Vector (x,y) in the cross-sectional plane '//&
                                                                    'defined by the length and axial vector', numberedmulti=.TRUE., no=2)
CALL prms%CreateRealOption(     'Coil[$]-Segment[$]-Radius'       , 'Circle segment: Radius in the cross-sectional plane '//&
                                                                    'defined by the length and axial vector', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Coil[$]-Segment[$]-Phi1'         , 'Circle segment: Initial angle in [deg]', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Coil[$]-Segment[$]-Phi2'         , 'Circle segment: Final angle in [deg]', numberedmulti=.TRUE.)

! Circle coils
CALL prms%CreateRealOption(     'Coil[$]-Radius'        , 'Radius for circular coils', numberedmulti=.TRUE.)

! Rectangle coils
CALL prms%CreateRealArrayOption('Coil[$]-RectVec1'      , 'Vector 1 (x,y) in the cross-sectional plane defined by the length '//&
                                                          'and axial vector, spanning the rectangular coil', numberedmulti=.TRUE., no=2)
CALL prms%CreateRealArrayOption('Coil[$]-RectVec2'      , 'Vector 2 (x,y) in the cross-sectional plane defined by the length '//&
                                                          'and axial vector, spanning the rectangular coil', numberedmulti=.TRUE., no=2)

! Time-dependent coils
CALL prms%SetSection('Time-dependent coils')
CALL prms%CreateLogicalOption(  'Coil[$]-TimeDepCoil'     , 'Use time-dependent current (sinusoidal curve) for coil', &
                                                            '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Coil[$]-CurrentFrequency', 'Current frequency [1/s]', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Coil[$]-CurrentPhase'    , 'Current phase shift [rad]','0.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'nTimePoints'             , 'Number of points for the discretization of the sinusoidal curve')

END SUBROUTINE DefineParametersSuperB


SUBROUTINE InitializeSuperB()
!===================================================================================================================================
!> Read-in of SuperB parameters for permanent magnets, coils and time-dependent coils
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_SuperB_Vars
USE MOD_Globals_Vars       ,ONLY: PI
USE MOD_Interpolation_Vars ,ONLY: BGFieldVTKOutput
USE MOD_ReadInTools        ,ONLY: PrintOption
USE MOD_Interpolation_Vars ,ONLY: NodeTypeGL,NodeType
USE MOD_Mesh_Vars          ,ONLY: Vdm_GL_N
USE MOD_Interpolation_Vars ,ONLY: NAnalyze
USE MOD_Interpolation      ,ONLY: GetVandermonde
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iMagnet, iCoil, iSegment, NbrOfTimeDepCoils
CHARACTER(LEN=32)  :: hilf,hilf2
REAL               :: FrequencyTmp
!===================================================================================================================================

! Get logical for calculating the error norms L2 and LInf of magnetic field
DoCalcErrorNormsSuperB = GETLOGICAL('DoCalcErrorNormsSuperB')

IF(DoCalcErrorNormsSuperB)THEN
  ! Get Vandermonde from Gauss-Lobatto (GL) nodes to Gauss (G)
  ALLOCATE(Vdm_GL_N(0:PP_N,0:NAnalyze))
  CALL GetVandermonde(NAnalyze , NodeTypeGL   , PP_N , NodeType , Vdm_GL_N , modal=.FALSE.)
END IF ! DoCalcErrorNormsSuperB

! Output of the magnets/coils as separate VTK files
BGFieldVTKOutput     = GETLOGICAL('PIC-CalcBField-OutputVTK')

! Get the number of magnets
NumOfPermanentMagnets= GETINT('NumOfPermanentMagnets')
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
      ! Default value for Radius2 is required because numberedmulti=.TRUE.
      ! For hollow cylinder (ring magnet) choose Radius2 < Radius
      PermanentMagnetInfo(iMagnet)%Radius2            = GETREAL('PermanentMagnet'//TRIM(hilf)//'-Radius2','0.')
      PermanentMagnetInfo(iMagnet)%HeightVector(1:3)  = GETREALARRAY('PermanentMagnet'//TRIM(hilf)//'-HeightVector',3)

      IF(PermanentMagnetInfo(iMagnet)%Radius2.GT.0.0.AND.&
         PermanentMagnetInfo(iMagnet)%Radius2.GT.PermanentMagnetInfo(iMagnet)%Radius)THEN
         CALL abort(&
         __STAMP__&
         ,'Cylindrical magnet: Radius2 cannot be larger than Radius!')
      END IF ! PermanentMagnetInfo(iMagnet)%Radius2.GT.0.0.AND.
    CASE('conic')
      PermanentMagnetInfo(iMagnet)%Radius             = GETREAL('PermanentMagnet'//TRIM(hilf)//'-Radius')
      PermanentMagnetInfo(iMagnet)%Radius2            = GETREAL('PermanentMagnet'//TRIM(hilf)//'-Radius2')
      PermanentMagnetInfo(iMagnet)%HeightVector(1:3)  = GETREALARRAY('PermanentMagnet'//TRIM(hilf)//'-HeightVector',3)
    CASE DEFAULT
      CALL abort(__STAMP__, &
        'ERROR SuperB: Given permanent magnet geometry is not implemented! Permanent magnet number:', iMagnet)
    END SELECT

    ! Sanity Checks
    SELECT CASE(TRIM(PermanentMagnetInfo(iMagnet)%Type))
    CASE('sphere','cylinder','conic')
      IF(PermanentMagnetInfo(iMagnet)%Radius.LE.0.0)THEN
        CALL abort(&
        __STAMP__&
        ,'sphere/cylinder/conic magnet: Radius cannot be <= 0.0')
      END IF ! PermanentMagnetInfo(iMagnet)%Radius.LE.0.0
    END SELECT
  END DO
END IF

! Get the number of coils/conductors
NumOfCoils        = GETINT('NumOfCoils','0')
NbrOfTimeDepCoils = 0. ! Initialize
FrequencyTmp      = 0. ! Initialize
ALLOCATE(CoilInfo(NumOfCoils))
ALLOCATE(TimeDepCoil(NumOfCoils))
ALLOCATE(CurrentInfo(NumOfCoils))
! Read-in of coil/conductor parameters
IF (NumOfCoils.GT.0) THEN
  DO iCoil = 1,NumOfCoils
    !SWRITE(*,*) "|       Read-in infos of coil |", iCoil
    CALL PrintOption('Read-in infos of coil number','superB',IntOpt=iCoil)
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
    CASE DEFAULT
      CALL abort(__STAMP__, &
        'ERROR SuperB: Given coil geometry is not implemented! Coil number:', iCoil)
    END SELECT
    ! --------------------- Time-dependent current ---------------------------------------
    TimeDepCoil(iCoil) = GETLOGICAL('Coil'//TRIM(hilf)//'-TimeDepCoil')
    IF(TimeDepCoil(iCoil)) THEN
      NbrOfTimeDepCoils = NbrOfTimeDepCoils + 1
      CurrentInfo(iCoil)%CurrentFreq = GETREAL('Coil'//TRIM(hilf)//'-CurrentFrequency')
      CurrentInfo(iCoil)%CurrentPhase = GETREAL('Coil'//TRIM(hilf)//'-CurrentPhase')
      ! Check that all time-dependent coils use the same frequency
      IF((NbrOfTimeDepCoils.GT.1).AND.(.NOT.ALMOSTEQUALRELATIVE(CurrentInfo(iCoil)%CurrentFreq,FrequencyTmp,1e-5)))THEN
        CALL abort(__STAMP__,'All time-dependent coils must have the same frequency!')
      END IF
      FrequencyTmp = CurrentInfo(iCoil)%CurrentFreq
    END IF
    CoilInfo(iCoil)%Current = GETREAL('Coil'//TRIM(hilf)//'-Current')
  END DO
END IF

! Discretisation in time
UseTimeDepCoil=.FALSE.
IF (ANY(TimeDepCoil)) THEN
  UseTimeDepCoil = .TRUE.
  nTimePoints    = GETINT('nTimePoints')
  IF(nTimePoints.LT.2) CALL abort(__STAMP__,'nTimePoints cannot be smaller than 2')
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
SDEALLOCATE(TimeDepCoil)
SDEALLOCATE(BGFieldTDep)
END SUBROUTINE FinalizeSuperB

END MODULE MOD_SuperB_Init
