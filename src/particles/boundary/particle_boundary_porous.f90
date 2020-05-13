!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_Particle_Boundary_Porous
!===================================================================================================================================
! Porous Boundary Condition: utilized as a pump with a target pressure, might be extended to other conditions
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitPorousBoundaryCondition
  MODULE PROCEDURE InitPorousBoundaryCondition
END INTERFACE

INTERFACE PorousBoundaryTreatment
  MODULE PROCEDURE PorousBoundaryTreatment
END INTERFACE

INTERFACE PorousBoundaryRemovalProb_Pressure
  MODULE PROCEDURE PorousBoundaryRemovalProb_Pressure
END INTERFACE

PUBLIC::DefineParametersPorousBC
PUBLIC::InitPorousBoundaryCondition, InitPorousCommunication
PUBLIC::PorousBoundaryTreatment, PorousBoundaryRemovalProb_Pressure
PUBLIC::FinalizePorousBoundaryCondition
!===================================================================================================================================

CONTAINS

SUBROUTINE DefineParametersPorousBC()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Porous BC")

CALL prms%CreateIntOption(      'Part-nPorousBC'&
                                , 'Number of porous boundary conditions', '0')
CALL prms%CreateIntOption(      'Part-PorousBC-IterationMacroVal' &
                                , 'Number of iterations the pressure and density will be sampled before updating the value. '//&
                                  'Alternative is to constantly update with a relaxation factor', '0')
CALL prms%CreateIntOption(      'Part-PorousBC[$]-BC' &
                                , 'PartBound to be a porous boundary', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-PorousBC[$]-Pressure' &
                                , 'Pressure [Pa] at the porous boundary', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-PorousBC[$]-Temperature' &
                                , 'Temperature [K] at the porous boundary', numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Part-PorousBC[$]-Type' &
                                , 'Define the type of porous boundary, currently available are sensor and pump.' //&
                                  'Option sensor: Using the defined region/BC to measure the pressure difference between the '//&
                                  'defined pressure and the pressure at the sensor. Option pump: Use the defined region/BC as '//&
                                  'as a pump by supplying a constant -PumpingSpeed or by supplying -DeltaPumpingSpeed-Kp ' //&
                                  'and -DeltaPumpingSpeed-Ki to adapt the pumping speed until the target pressure is reached.' &
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-PorousBC[$]-PumpingSpeed' &
                                , 'PorousBC used as a pump: Removal probability at the BC is determined through the pumping ' //&
                                  'speed, which can be controlled to achieve the given pressure. Initial pumping speed [m3/s] ' //&
                                  'can be zero', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-PorousBC[$]-DeltaPumpingSpeed-Kp' &
                                , 'Proportional factor for the pumping speed controller, Kp=0 -> constant given pumping speed' &
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-PorousBC[$]-DeltaPumpingSpeed-Ki' &
                                , 'Integral factor for the pumping speed controller, Ki=0 -> only proportional part is used' &
                                , numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Part-PorousBC[$]-Region' &
                                , 'Define a region, allowing small geometrical features such as holes w/o meshing. ' //&
                                  'Option: circular, requires the definition of the parameters -normalDir, origin, rmax/rmin' &
                                , 'none', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-PorousBC[$]-normalDir' &
                                , 'Definition of the normal direction (as an integer, 1: x, 2: y, 3: z) of the boundary on ' //&
                                  'which the region shall be created.' &
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-PorousBC[$]-origin' &
                                , 'Coordinates of the middle point of the region, Example: normalDir=1: (/y,z/), ' //&
                                  'normalDir=2: (/z,x/), normalDir=3: (/x,y/)', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-PorousBC[$]-rmax' &
                                , 'Maximum radius [m] of the circular region', '1e21', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-PorousBC[$]-rmin' &
                                , 'Minimal radius [m] of the circular region', '0.', numberedmulti=.TRUE.)
END SUBROUTINE DefineParametersPorousBC


SUBROUTINE InitPorousBoundaryCondition()
!===================================================================================================================================
! 1) Read-in of parameters
! 2) Mapping of the surface side to a porous BC and porous BC
! 3) Allocating the PorousBC arrays per BC, containing only sides with a (partially) porous side
! 4) Mapping of the porous BC sides to the surfaces sides, initialization of the pumping capacity and treatment of regions
!    (determination, which cells are competely and partially inside)
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Mesh_Vars                   ,ONLY: nElems, offsetElem
USE MOD_Particle_Mesh_Vars          ,ONLY: SideInfo_Shared
USE MOD_Particle_Vars               ,ONLY: nSpecies, Adaptive_MacroVal, Symmetry2D, Symmetry2DAxisymmetric
USE MOD_Particle_Boundary_Vars      ,ONLY: PartBound, nPorousBC, nPorousSides, PorousBCMacroVal
USE MOD_Particle_Boundary_Vars      ,ONLY: nComputeNodeSurfTotalSides, SurfSide2GlobalSide
USE MOD_Particle_Boundary_Vars      ,ONLY: PorousBCSampIter, PorousBC
USE MOD_Particle_Boundary_Vars      ,ONLY: MapSurfSideToPorousSide_Shared, MapSurfSideToPorousSide_Shared_Win
USE MOD_Particle_Boundary_Vars      ,ONLY: PorousBCSampWall, PorousBCOutput
USE MOD_Particle_Boundary_Vars      ,ONLY: PorousBCInfo_Shared,PorousBCProperties_Shared,PorousBCSampWall_Shared
USE MOD_Particle_Boundary_Vars      ,ONLY: PorousBCInfo_Shared_Win,PorousBCProperties_Shared_Win,PorousBCSampWall_Shared_Win
USE MOD_Particle_Tracking_Vars      ,ONLY: DoRefMapping
#if USE_MPI
USE MOD_MPI_Shared              ,ONLY: Allocate_Shared
USE MOD_MPI_Shared_Vars              ,ONLY: MPI_COMM_SHARED, myComputeNodeRank
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSurfSide, iPBC, iSurfSideTmp, GlobalSideID, iPorousSide
INTEGER,ALLOCATABLE   :: MapSurfSideToPorousBC_Temp(:,:)
CHARACTER(32)         :: hilf, hilf2
REAL                  :: rmin, rmax
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND)         :: MPISharedSize
#endif /*USE_MPI*/
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)') ' INIT POROUS BOUNDARY CONDITION ...'

IF(DoRefMapping) THEN
  CALL abort(__STAMP__&
      ,'ERROR: Porous boundary conditions are not implemented with DoRefMapping!')
END IF

IF(Symmetry2D.AND.(.NOT.Symmetry2DAxisymmetric)) THEN
  CALL abort(__STAMP__&
      ,'ERROR: Porous boundary conditions are not implemented for 2D simulations!')
END IF

! 1) Read-in of parameters
PorousBCSampIter = GETINT('Part-PorousBC-IterationMacroVal', '0')
IF(PorousBCSampIter.GT.0) THEN
  ALLOCATE(PorousBCMacroVal(1:8,1:nElems,1:nSpecies))
  PorousBCMacroVal = 0.0
END IF

ALLOCATE(PorousBC(1:nPorousBC))
ALLOCATE(MapSurfSideToPorousBC_Temp(1:2,1:nComputeNodeSurfTotalSides))
MapSurfSideToPorousBC_Temp = 0

#if USE_MPI
MPISharedSize = INT((nComputeNodeSurfTotalSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nComputeNodeSurfTotalSides/),MapSurfSideToPorousSide_Shared_Win,MapSurfSideToPorousSide_Shared)
CALL MPI_WIN_LOCK_ALL(0,MapSurfSideToPorousSide_Shared_Win,IERROR)
IF (myComputeNodeRank.EQ.0) THEN
  MapSurfSideToPorousSide_Shared = 0
END IF
! This barrier MIGHT not be required
CALL MPI_WIN_SYNC(MapSurfSideToPorousSide_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#else
ALLOCATE(MapSurfSideToPorousSide_Shared(1:nComputeNodeSurfTotalSides))
MapSurfSideToPorousSide_Shared = 0
#endif /*USE_MPI*/

DO iPBC = 1, nPorousBC
  WRITE(UNIT=hilf,FMT='(I0)') iPBC
  PorousBC(iPBC)%BC = GETINT('Part-PorousBC'//TRIM(hilf)//'-BC')
  IF(PartBound%TargetBoundCond(PorousBC(iPBC)%BC).NE.PartBound%ReflectiveBC) THEN
    CALL abort(__STAMP__&
      ,'ERROR in init of porous BC: given boundary condition must be reflective!')
  END IF
  ! Read-in of the conditions at the porous boundary
  PorousBC(iPBC)%Pressure  = GETREAL('Part-PorousBC'//TRIM(hilf)//'-Pressure')
  PorousBC(iPBC)%Temperature  = GETREAL('Part-PorousBC'//TRIM(hilf)//'-Temperature')
  PorousBC(iPBC)%Type  = TRIM(GETSTR('Part-PorousBC'//TRIM(hilf)//'-Type'))
  SELECT CASE(PorousBC(iPBC)%Type)
    CASE('sensor')
      PorousBC(iPBC)%PumpingSpeed = 0.0
      PorousBC(iPBC)%DeltaPumpingSpeedKp = 0.0
      PorousBC(iPBC)%DeltaPumpingSpeedKi = 0.0
    CASE('pump')
      ! Initial pumping speed at the porous boundary [m3/s]
      PorousBC(iPBC)%PumpingSpeed = GETREAL('Part-PorousBC'//TRIM(hilf)//'-PumpingSpeed')
      ! Proportional and integral factors for the control of the pumping speed
      PorousBC(iPBC)%DeltaPumpingSpeedKp = GETREAL('Part-PorousBC'//TRIM(hilf)//'-DeltaPumpingSpeed-Kp')
      PorousBC(iPBC)%DeltaPumpingSpeedKi = GETREAL('Part-PorousBC'//TRIM(hilf)//'-DeltaPumpingSpeed-Ki')
      ! Determining the order of magnitude
      IF(PorousBC(iPBC)%Pressure.GT.0.0) THEN
        PorousBC(iPBC)%DeltaPumpingSpeedKp = PorousBC(iPBC)%DeltaPumpingSpeedKp / 10.0**(ANINT(LOG10(PorousBC(iPBC)%Pressure)))
        PorousBC(iPBC)%DeltaPumpingSpeedKi = PorousBC(iPBC)%DeltaPumpingSpeedKi / 10.0**(ANINT(LOG10(PorousBC(iPBC)%Pressure)))
      END IF
    CASE DEFAULT
      CALL abort(__STAMP__&
      ,'ERROR in type definition of porous bc:', iPBC)
  END SELECT
  PorousBC(iPBC)%Region = TRIM(GETSTR('Part-PorousBC'//TRIM(hilf)//'-Region','none'))
  IF(PorousBC(iPBC)%Region.EQ.'none') THEN
    PorousBC(iPBC)%UsingRegion = .FALSE.
  ELSE
    PorousBC(iPBC)%UsingRegion = .TRUE.
    ! Normal direction to the surface (could be replaced by the normal vector of the boundary?)
    PorousBC(iPBC)%dir(1) = GETINT('Part-PorousBC'//TRIM(hilf)//'-normalDir')
    IF (PorousBC(iPBC)%dir(1).EQ.1) THEN
        PorousBC(iPBC)%dir(2)=2
        PorousBC(iPBC)%dir(3)=3
    ELSE IF (PorousBC(iPBC)%dir(1).EQ.2) THEN
        PorousBC(iPBC)%dir(2)=3
        PorousBC(iPBC)%dir(3)=1
    ELSE IF (PorousBC(iPBC)%dir(1).EQ.3) THEN
        PorousBC(iPBC)%dir(2)=1
        PorousBC(iPBC)%dir(3)=2
    ELSE
      CALL abort(__STAMP__&
        ,'ERROR in init: normalDir for PorousBC must be between 1 and 3!')
    END IF
    SELECT CASE(PorousBC(iPBC)%Region)
      CASE('circular')
        PorousBC(iPBC)%origin = GETREALARRAY('Part-PorousBC'//TRIM(hilf)//'-origin',2)
        WRITE(UNIT=hilf2,FMT='(E16.8)') HUGE(PorousBC(iPBC)%rmax)
        PorousBC(iPBC)%rmax = GETREAL('Part-PorousBC'//TRIM(hilf)//'-rmax',TRIM(hilf2))
        PorousBC(iPBC)%rmin = GETREAL('Part-PorousBC'//TRIM(hilf)//'-rmin','0.')
        IF(Symmetry2DAxisymmetric) THEN
          IF(PorousBC(iPBC)%dir(1).NE.1) THEN
            CALL abort(__STAMP__&
              ,'ERROR in Porous BC: For axisymmetric simulations, only regions perpendicular to the axis are allowed!', iPBC)
          END IF
          IF(PorousBC(iPBC)%origin(1)*PorousBC(iPBC)%origin(2).NE.0.0) THEN
            CALL abort(__STAMP__&
              ,'ERROR in Porous BC: For axisymmetric simulations, the origin has to be at (0,0)!', iPBC)
          END IF
        END IF
      CASE DEFAULT
        CALL abort(__STAMP__&
        ,'ERROR in region definition of porous bc:', iPBC)
    END SELECT
  END IF    ! Region is given
END DO      ! nPorousBC

! 2) Mapping of the porous BC sides to the respective surface side
IF (myComputeNodeRank.EQ.0) THEN
  nPorousSides = 0
  iSurfSideTmp = 0
  DO iSurfSide=1,nComputeNodeSurfTotalSides
    GlobalSideID = SurfSide2GlobalSide(SURF_SIDEID,iSurfSide)
    ! Skip not reflective BC sides
    IF(PartBound%TargetBoundCond(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,GlobalSideID))).NE.PartBound%ReflectiveBC) CYCLE
    DO iPBC = 1, nPorousBC
      IF(PorousBC(iPBC)%BC.EQ.PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,GlobalSideID))) THEN
        ! Determine which cells are inside/outside/partially inside the defined region (0: inside, 1: partially)
        IF(PorousBC(iPBC)%UsingRegion) THEN
          SELECT CASE(PorousBC(iPBC)%Region)
            CASE('circular')
              CALL GetRadialDistance2D(GlobalSideID,PorousBC(iPBC)%dir,PorousBC(iPBC)%origin,rmin,rmax)
              IF ( (rmin .GT. PorousBC(iPBC)%rmax) .OR. (rmax .LT. PorousBC(iPBC)%rmin) ) CYCLE
          END SELECT
        END IF
        IF(iSurfSide.EQ.iSurfSideTmp) THEN
          CALL abort(__STAMP__&
          ,'ERROR in Porous BC: Side is already defined by another porous BC. Make sure the BCs do not overlap!')
        END IF
        nPorousSides = nPorousSides + 1
        MapSurfSideToPorousBC_Temp(1,nPorousSides) = iPBC
        MapSurfSideToPorousBC_Temp(2,nPorousSides) = iSurfSide
        MapSurfSideToPorousSide_Shared(iSurfSide) = nPorousSides
        iSurfSideTmp = iSurfSide
      END IF
    END DO
  END DO
END IF

#if USE_MPI
CALL MPI_BCAST(nPorousSides,1,MPI_INTEGER,0,MPI_COMM_SHARED,IERROR)

MPISharedSize = INT((3*nPorousSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/3,nPorousSides/),PorousBCInfo_Shared_Win,PorousBCInfo_Shared)
MPISharedSize = INT((2*nPorousSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/2,nPorousSides/),PorousBCProperties_Shared_Win,PorousBCProperties_Shared)
CALL Allocate_Shared(MPISharedSize,(/2,nPorousSides/),PorousBCSampWall_Shared_Win,PorousBCSampWall_Shared)
CALL MPI_WIN_LOCK_ALL(0,PorousBCInfo_Shared_Win,IERROR)
CALL MPI_WIN_LOCK_ALL(0,PorousBCProperties_Shared_Win,IERROR)
CALL MPI_WIN_LOCK_ALL(0,PorousBCSampWall_Shared_Win,IERROR)
MPISharedSize = INT((5*nPorousSides),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
IF (myComputeNodeRank.EQ.0) THEN
  PorousBCInfo_Shared(1,1:nPorousSides) = MapSurfSideToPorousBC_Temp(1,1:nPorousSides)
  PorousBCInfo_Shared(2,1:nPorousSides) = MapSurfSideToPorousBC_Temp(2,1:nPorousSides)
  PorousBCInfo_Shared(3,1:nPorousSides) = -1
  PorousBCProperties_Shared = 0.
  PorousBCSampWall_Shared = 0.
END IF
! This barrier MIGHT not be required
CALL MPI_WIN_SYNC(MapSurfSideToPorousSide_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(PorousBCInfo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(PorousBCProperties_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(PorousBCSampWall_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
#else
ALLOCATE(PorousBCInfo_Shared(1:3,1:nPorousSides))               ! 1: Porous BC,   2: SurfSide ID, 3: SideType
ALLOCATE(PorousBCProperties_Shared(1:2,1:nPorousSides))         ! 1: Remove prob, 2: Pumping speed side
ALLOCATE(PorousBCSampWall_Shared(1:2,1:nPorousSides))           ! 1: Impinging,   2: Deleted particles count
PorousBCInfo_Shared(1,1:nPorousSides) = MapSurfSideToPorousBC_Temp(1,1:nPorousSides)
PorousBCInfo_Shared(2,1:nPorousSides) = MapSurfSideToPorousBC_Temp(2,1:nPorousSides)
PorousBCInfo_Shared(3,1:nPorousSides) = -1
PorousBCProperties_Shared = 0.
PorousBCSampWall_Shared = 0.
#endif /*USE_MPI*/

ALLOCATE(PorousBCSampWall(1:2,1:nPorousSides))
PorousBCSampWall = 0.
ALLOCATE(PorousBCOutput(1:5,1:nPorousBC))
PorousBCOutput = 0.
IF (myComputeNodeRank.EQ.0) THEN
  ! 4) Mapping the porous BC ID and the porous BC side ID to a surface side
  DO iPorousSide=1,nPorousSides
    iPBC = PorousBCInfo_Shared(1,iPorousSide)
    PorousBCProperties_Shared(2,iPorousSide) = PorousBC(iPBC)%PumpingSpeed
    GlobalSideID = SurfSide2GlobalSide(SURF_SIDEID,PorousBCInfo_Shared(2,iPorousSide))
    ! If a pumping speed was read-in during restart, use it (only for local elements and not halo sides)
    ! NOTE: DOES NOT WORK CURRENTLY AS ONLY NODE-LEADER IS PERFORMING THIS LOOP, HOWEVER ADAPTIVEMACROVAL IS STILL PROC LOCAL
    ! IF(Adaptive_MacroVal(11,SideInfo_Shared(SIDE_ELEMID,GlobalSideID)-offsetElem,1).GT.0.0) THEN
    !     PorousBCProperties_Shared(2,iPorousSide) = Adaptive_MacroVal(11,SideInfo_Shared(SIDE_ELEMID,GlobalSideID)-offsetElem,1)
    ! END IF
    ! Determine which cells are inside/outside/partially inside the defined region (0: inside, 1: partially)
    IF(PorousBC(iPBC)%UsingRegion) THEN
      SELECT CASE(PorousBC(iPBC)%Region)
        CASE('circular')
          CALL GetRadialDistance2D(GlobalSideID,PorousBC(iPBC)%dir,PorousBC(iPBC)%origin,rmin,rmax)
          IF ( (rmax .LE. PorousBC(iPBC)%rmax) .AND. (rmin .GE. PorousBC(iPBC)%rmin) ) THEN
            PorousBCInfo_Shared(3,iPorousSide)=0
          ELSE
            PorousBCInfo_Shared(3,iPorousSide)=1
          END IF
      END SELECT
      IF(PorousBCInfo_Shared(3,iPorousSide).LT.0) THEN
          CALL abort(__STAMP__&
          ,'ERROR in Porous BC: Region side type not defined! Porous BC ID:', iPBC)
      END IF
    END IF
  END DO
END IF

#if USE_MPI
CALL MPI_WIN_SYNC(PorousBCInfo_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(PorousBCProperties_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
IF (myComputeNodeRank.EQ.0) CALL InitPorousCommunication()
#endif

SDEALLOCATE(MapSurfSideToPorousBC_Temp)

SWRITE(UNIT_stdOut,'(A)') ' INIT POROUS BOUNDARY CONDITION DONE!'

END SUBROUTINE InitPorousBoundaryCondition

SUBROUTINE PorousBoundaryTreatment(iPart,SideID,alpha,PartTrajectory,ElasticReflectionAtPorousBC)
!===================================================================================================================================
! Treatment of particles impinging on the porous boundary
! 1) (Optional) When using regions on the BC, it is determined whether the particle hit the porous BC region or only the regular BC
! 2) Comparison of the removal probability with a random number to determine whether the particle is deleted. Counting number of
!    impinged and deleted particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars           ,ONLY: LastPartPos
USE MOD_Particle_Boundary_Vars  ,ONLY: PorousBCSampWall, PorousBC, MapSurfSideToPorousSide_Shared
USE MOD_Particle_Boundary_Vars  ,ONLY: PorousBCInfo_Shared, PorousBCProperties_Shared
USE MOD_part_tools              ,ONLY: GetParticleWeight
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound, GlobalSide2SurfSide
USE MOD_part_operations         ,ONLY: RemoveParticle
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPart, SideID
REAL, INTENT(IN)              :: PartTrajectory(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)            :: alpha
LOGICAL,INTENT(INOUT)         :: ElasticReflectionAtPorousBC
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: point(1:2), intersectionPoint(1:3), radius, iRan
INTEGER                       :: SurfSideID, iPBC, pBCSideID
LOGICAL                       :: ParticleHitPorousBC
!===================================================================================================================================
SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)

pBCSideID = MapSurfSideToPorousSide_Shared(SurfSideID)

IF(pBCSideID.GT.0) THEN
  iPBC = PorousBCInfo_Shared(1,pBCSideID)
  ParticleHitPorousBC = .TRUE.
  ! 1) Determination whether the particle hit the porous BC region or only the regular BC
  IF(PorousBC(iPBC)%UsingRegion) THEN
    IF(PorousBCInfo_Shared(3,pBCSideID).EQ.1) THEN
      ! Side is partially inside the porous region (check if its within bounds)
      intersectionPoint(1:3) = LastPartPos(1:3,iPart) + alpha*PartTrajectory(1:3)
      point(1)=intersectionPoint(PorousBC(iPBC)%dir(2))-PorousBC(iPBC)%origin(1)
      point(2)=intersectionPoint(PorousBC(iPBC)%dir(3))-PorousBC(iPBC)%origin(2)
      radius=SQRT( (point(1))**2+(point(2))**2 )
      ! Check if particle hit outside the region
      IF ((radius.GT.PorousBC(iPBC)%rmax).OR.(radius.LT.PorousBC(iPBC)%rmin)) THEN
        ParticleHitPorousBC = .FALSE.
      END IF
    END IF
  END IF
  ! 2) Comparison of the removal probability with a random number to determine whether the particle is deleted
  IF(ParticleHitPorousBC) THEN
    ! Counting particles that are impinging the porous BC (required for the calculation of the removal probability for the next dt)
    PorousBCSampWall(1,pBCSideID)   = PorousBCSampWall(1,pBCSideID) + GetParticleWeight(iPart)
    SELECT CASE(PorousBC(iPBC)%Type)
      CASE('sensor')
        ! Treat the sensor area as a regular boundary condition
        ElasticReflectionAtPorousBC = .FALSE.
      CASE('pump')
        CALL RANDOM_NUMBER(iRan)
        IF(iRan.LE.PorousBCProperties_Shared(1,pBCSideID)) THEN
          ! Counting particles that leave the domain through the porous BC (required for the calculation of the pumping capacity)
          PorousBCSampWall(2,pBCSideID) = PorousBCSampWall(2,pBCSideID) + GetParticleWeight(iPart)
          CALL RemoveParticle(iPart,BCID=PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)),alpha=alpha)
        END IF
        ! Treat the pump as a perfect reflection (particle is not removed and keeps its temperature)
        ElasticReflectionAtPorousBC = .TRUE.
    END SELECT
  END IF
END IF

END SUBROUTINE PorousBoundaryTreatment

SUBROUTINE PorousBoundaryRemovalProb_Pressure()
!===================================================================================================================================
! Determing the removal probability based on a target pressure at the porous boundary condition (can be used as a pump BC)
! 1) MPI communication of particles that impinged on halo sides to corresponding side
! 2) Loop over all porous BCs
!   a) Summing up the number of impinged particles for the whole BC surface
!   2.1) Loop over all sides within each porous BC: a), b) & d) is only performed if the pumping speed is adapted
!       a) Determining the delta between current gas mixture pressure in adjacent cell and target pressure
!       b) Adapting the pumping capacity (m^3/s) according to pressure difference (control through proportional and integral part)
!       c) Calculate the removal probability if any particles hit the pump
!       d) Limit removal probability to values between 0 and 1
!       e) Sampling of the pumping capacity (and other variables for PartAnalyze) for the output
! 3) MPI communication of the removal probability to halo sides
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars,              ONLY:DSMC, RadialWeighting
USE MOD_Mesh_Vars,              ONLY:nElems,offsetElem
USE MOD_Particle_Analyze_Vars,  ONLY:CalcPorousBCInfo
USE MOD_Particle_Boundary_Vars, ONLY:SurfOnNode, SurfSide2GlobalSide
USE MOD_Particle_Boundary_Vars, ONLY:nPorousBC, PorousBC, nPorousSides, PorousBCInfo_Shared, PorousBCOutput, PorousBCSampWall_Shared
USE MOD_Particle_Boundary_Vars, ONLY:PorousBCProperties_Shared, SampWallPumpCapacity
USE MOD_Particle_Vars,          ONLY:Species, nSpecies, Adaptive_MacroVal, usevMPF, VarTimeStep
USE MOD_Particle_VarTimeStep    ,ONLY: CalcVarTimeStep
USE MOD_Timedisc_Vars,          ONLY:dt
USE MOD_Particle_Mesh_Vars
#if USE_MPI
USE MOD_MPI_Shared_Vars         ,ONLY: myComputeNodeRank, MPI_COMM_SHARED,MPI_COMM_LEADERS_SURF
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: LocalElemID,GlobalElemID,SurfSideID,GlobalSideID,iPorousSide,iPBC
REAL                          :: PumpingSpeedTemp,DeltaPressure,partWeight,SumPartPorousBC,dtVar
REAL,ALLOCATABLE              :: SumPartImpinged(:)
!===================================================================================================================================

IF (.NOT.SurfOnNode) RETURN

IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
  partWeight = 1.
ELSE
  partWeight = Species(1)%MacroParticleFactor
END IF

ALLOCATE(SumPartImpinged(nPorousBC))
SumPartImpinged = 0
! 1) MPI communication of particles that impinged on halo sides to corresponding side
#if USE_MPI
CALL ExchangeImpingedPartPorousBC()
#else
PorousBCSampWall_Shared(1:2,1:nPorousSides) = PorousBCSampWall(1:2,1:nPorousSides)
#endif
DO iPorousSide = 1, nPorousSides
  iPBC = PorousBCInfo_Shared(1,iPorousSide)
  SumPartImpinged(iPBC) = SumPartImpinged(iPBC) + PorousBCSampWall_Shared(1,iPorousSide)
END DO

#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,SumPartImpinged(1:nPorousBC),nPorousBC,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEADERS_SURF,iError)
END IF
CALL MPI_BCAST(SumPartImpinged(1:nPorousBC),nPorousBC,MPI_DOUBLE_PRECISION,0,MPI_COMM_SHARED,iError)
#endif

! 2) Loop over all porous sides
DO iPorousSide = 1, nPorousSides
  SurfSideID = PorousBCInfo_Shared(2,iPorousSide)
  GlobalSideID = SurfSide2GlobalSide(SURF_SIDEID,SurfSideID)
  GlobalElemID = SideInfo_Shared(SIDE_ELEMID,GlobalSideID)
  ! Only treat your proc-local elements
  IF ((GlobalElemID.LT.1+offSetElem).OR.(GlobalElemID.GT.nElems+offSetElem)) CYCLE
  LocalElemID = GlobalElemID - offsetElem
  iPBC = PorousBCInfo_Shared(1,iPorousSide)
  SumPartPorousBC = SumPartImpinged(iPBC)
  ! Zero the output variable
  IF(CalcPorousBCInfo) PorousBCOutput = 0.
  ! Skip element if number density is zero
  IF(SUM(Adaptive_MacroVal(DSMC_NUMDENS,LocalElemID,1:nSpecies)).EQ.0.0) CYCLE
  ! Get the correct time step of the cell
  IF(VarTimeStep%UseVariableTimeStep) THEN
    dtVar = dt * CalcVarTimeStep(ElemMidPoint_Shared(1,GlobalElemID), ElemMidPoint_Shared(2,GlobalElemID), LocalElemID)
  ELSE
    dtVar = dt
  END IF
  ! Determine the removal probability based on the pumping speed (adaptive to a target pressure or fixed)
  IF((PorousBC(iPBC)%DeltaPumpingSpeedKp.GT.0.).OR.(PorousBC(iPBC)%DeltaPumpingSpeedKi.GT.0.)) THEN
    ! a) Determining the delta between current gas mixture pressure in adjacent cell and target pressure
    DeltaPressure = SUM(Adaptive_MacroVal(12,LocalElemID,1:nSpecies))-PorousBC(iPBC)%Pressure
    ! Integrating the pressure difference (only utilized later if DeltaPumpingSpeedKi was given)
    IF(PorousBCProperties_Shared(2,iPorousSide).GT.0.0) THEN
      Adaptive_MacroVal(13,LocalElemID,1) = Adaptive_MacroVal(13,LocalElemID,1) + DeltaPressure * dtVar
    ELSE
      Adaptive_MacroVal(13,LocalElemID,1) = 0.0
    END IF
    ! b) Adapting the pumping capacity (m^3/s) according to pressure difference (control through proportional and integral part)
    PumpingSpeedTemp = PorousBCProperties_Shared(2,iPorousSide) + PorousBC(iPBC)%DeltaPumpingSpeedKp * DeltaPressure &
        + PorousBC(iPBC)%DeltaPumpingSpeedKi * Adaptive_MacroVal(13,LocalElemID,1)
    ! c) Calculate the removal probability if any particles hit the pump
    IF(SumPartPorousBC.GT.0) THEN
      PorousBCProperties_Shared(1,iPorousSide) = PumpingSpeedTemp*SUM(Adaptive_MacroVal(DSMC_NUMDENS,LocalElemID,1:nSpecies)) &
                                                      * dtVar / (SumPartPorousBC*partWeight)
    ELSE
      PorousBCProperties_Shared(1,iPorousSide) = 0.0
    END IF
    ! d) Limit removal probability to values between 0 and 1
    IF(PorousBCProperties_Shared(1,iPorousSide).GT.1.0) THEN
      PorousBCProperties_Shared(1,iPorousSide) = 1.0
      ! Setting pumping speed to maximum value (alpha=1)
      PorousBCProperties_Shared(2,iPorousSide) = SumPartPorousBC*partWeight &
                                                    / (SUM(Adaptive_MacroVal(DSMC_NUMDENS,LocalElemID,1:nSpecies))*dtVar)
    ELSE IF(PorousBCProperties_Shared(1,iPorousSide).LE.0.0) THEN
      PorousBCProperties_Shared(1,iPorousSide) = 0.0
      ! Avoiding negative pumping speeds
      PorousBCProperties_Shared(2,iPorousSide) = 0.0
    ELSE
      ! Only adapting the pumping speed if alpha is between zero and one
      PorousBCProperties_Shared(2,iPorousSide) = PumpingSpeedTemp
    END IF
  ELSE IF(PorousBCProperties_Shared(2,iPorousSide).GT.0.0) THEN
    ! Constant given pumping speed
    IF(SumPartPorousBC.GT.0) THEN
      PorousBCProperties_Shared(1,iPorousSide) = PorousBCProperties_Shared(2,iPorousSide) &
                      * SUM(Adaptive_MacroVal(DSMC_NUMDENS,LocalElemID,1:nSpecies)) * dtVar / (SumPartPorousBC*partWeight)
    ELSE
      PorousBCProperties_Shared(1,iPorousSide) = 0.0
    END IF
  END IF
  ! Storing the pumping speed for the restart state file
  Adaptive_MacroVal(11,LocalElemID,1) = PorousBCProperties_Shared(2,iPorousSide)
  ! e) Sampling of the pumping capacity (and other variables for PartAnalyze) for the output
  ! -------- Sampling for output in DSMCSurfState --------------------------------------------------------------------------------
  IF(DSMC%CalcSurfaceVal) THEN
    SampWallPumpCapacity(SurfSideID) = SampWallPumpCapacity(SurfSideID) + PorousBCSampWall_Shared(2,iPorousSide) &
                                                * partWeight / (SUM(Adaptive_MacroVal(DSMC_NUMDENS,LocalElemID,1:nSpecies))*dtVar)
  END IF
  ! -------- Sampling for output in PartAnalyze ----------------------------------------------------------------------------------
  IF(CalcPorousBCInfo) THEN
    ! Sampling the actual instantaneous pumping speed S (m^3/s) through the number of deleted particles  (-PumpSpeed-Measure-)
    PorousBCOutput(2,iPBC) = PorousBCOutput(2,iPBC) + PorousBCSampWall_Shared(2,iPorousSide) * PartWeight &
                                                          / (SUM(Adaptive_MacroVal(DSMC_NUMDENS,LocalElemID,1:nSpecies))*dtVar)
    IF(PorousBCSampWall_Shared(1,iPorousSide).GT.0) THEN
      ! Counting only sides, where a particle hit the pump
      PorousBCOutput(1,iPBC) = PorousBCOutput(1,iPBC) + 1.
      ! Pumping speed S (m^3/s) used at the pump to calculate the removal probability (-PumpSpeed-Control-)
      PorousBCOutput(3,iPBC) = PorousBCOutput(3,iPBC) + PorousBCProperties_Shared(2,iPorousSide)
      ! Removal probability
      PorousBCOutput(4,iPBC) = PorousBCOutput(4,iPBC) + PorousBCProperties_Shared(1,iPorousSide)
      ! Normalized pressure at the pump
      PorousBCOutput(5,iPBC) = PorousBCOutput(5,iPBC) + SUM(Adaptive_MacroVal(12,LocalElemID,1:nSpecies)) / PorousBC(iPBC)%Pressure
    END IF
  END IF
  ! Reset of the sampled particle numbers at the pump
  PorousBCSampWall_Shared = 0
END DO    ! iPorousSide = 1, nPorousSides

! Exchange removal probability for halo cells
#if USE_MPI
CALL ExchangeRemovalProbabilityPorousBC()
#endif

END SUBROUTINE PorousBoundaryRemovalProb_Pressure


#if USE_MPI
SUBROUTINE ExchangeImpingedPartPorousBC()
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED,MPI_COMM_LEADERS_SURF
USE MOD_MPI_Shared_Vars         ,ONLY: nSurfLeaders,myComputeNodeRank,mySurfRank
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfOnNode
USE MOD_Particle_Boundary_Vars  ,ONLY: MapSurfSideToPorousSide_Shared
USE MOD_Particle_Boundary_Vars  ,ONLY: GlobalSide2SurfSide
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfMapping
USE MOD_Particle_Boundary_Vars  ,ONLY: nPorousSides, PorousBCSampWall, PorousBCSampWall_Shared, PorousBCSampWall_Shared_Win
USE MOD_Particle_MPI_Vars       ,ONLY: PorousBCSendBuf,PorousBCRecvBuf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iProc,SideID,PorousSideID,iPos
INTEGER                         :: MessageSize,iSurfSide,SurfSideID
INTEGER                         :: nValues
INTEGER                         :: RecvRequest(0:nSurfLeaders-1),SendRequest(0:nSurfLeaders-1)
!===================================================================================================================================
! nodes without sampling surfaces do not take part in this routine
IF (.NOT.SurfOnNode) RETURN

nValues = 2
! collect the information from the proc-local shadow arrays in the compute-node shared array
MessageSize = nValues*nPorousSides
CALL MPI_REDUCE(PorousBCSampWall,PorousBCSampWall_Shared,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
CALL MPI_WIN_SYNC(PorousBCSampWall_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)

! prepare buffers for surf leader communication
IF (myComputeNodeRank.EQ.0) THEN
  ! open receive buffer
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    ! Only open recv buffer if we are expecting sides from this leader node
    IF (SurfMapping(iProc)%nRecvPorousSides.EQ.0) CYCLE

    ! Message is sent on MPI_COMM_LEADERS_SURF, so rank is indeed iProc
    MessageSize = SurfMapping(iProc)%nRecvPorousSides * nValues
    CALL MPI_IRECV( PorousBCRecvBuf(iProc)%content               &
                  , MessageSize                                  &
                  , MPI_DOUBLE_PRECISION                         &
                  , iProc                                        &
                  , 1209                                         &
                  , MPI_COMM_LEADERS_SURF                        &
                  , RecvRequest(iProc)                           &
                  , IERROR)
  END DO ! iProc

  ! build message
  DO iProc = 0,nSurfLeaders-1
    ! Ignore myself
    IF (iProc .EQ. mySurfRank) CYCLE

    ! Only assemble message if we are expecting sides to send to this leader node
    IF (SurfMapping(iProc)%nSendPorousSides.EQ.0) CYCLE

    ! Nullify everything
    iPos = 0
    PorousBCSendBuf(iProc)%content = 0.

    DO iSurfSide = 1,SurfMapping(iProc)%nSendPorousSides
      SideID     = SurfMapping(iProc)%SendPorousGlobalID(iSurfSide)
      SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)
      PorousSideID = MapSurfSideToPorousSide_Shared(SurfSideID)
      ! Assemble message
      PorousBCSendBuf(iProc)%content(iPos+1:iPos+nValues) = PorousBCSampWall_Shared(:,PorousSideID)
      iPos = iPos + nValues

      PorousBCSampWall_Shared(:,PorousSideID)=0.
    END DO ! iSurfSide=1,nSurfExchange%nSidesSend(iProc)
  END DO

  ! send message
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    ! Only open recv buffer if we are expecting sides from this leader node
    IF (SurfMapping(iProc)%nSendPorousSides.EQ.0) CYCLE

    ! Message is sent on MPI_COMM_LEADERS_SURF, so rank is indeed iProc
    MessageSize = SurfMapping(iProc)%nSendPorousSides * nValues
    CALL MPI_ISEND( PorousBCSendBuf(iProc)%content                   &
                  , MessageSize                                  &
                  , MPI_DOUBLE_PRECISION                         &
                  , iProc                                        &
                  , 1209                                         &
                  , MPI_COMM_LEADERS_SURF                        &
                  , SendRequest(iProc)                           &
                  , IERROR)
  END DO ! iProc

  ! Finish received number of sampling surfaces
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    IF (SurfMapping(iProc)%nSendPorousSides.NE.0) THEN
      CALL MPI_WAIT(SendRequest(iProc),MPIStatus,IERROR)
      IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error',IERROR)
    END IF

    IF (SurfMapping(iProc)%nRecvPorousSides.NE.0) THEN
      CALL MPI_WAIT(RecvRequest(iProc),MPIStatus,IERROR)
      IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error',IERROR)
    END IF
  END DO ! iProc

  ! add data do my list
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    ! Only open recv buffer if we are expecting sides from this leader node
    IF (SurfMapping(iProc)%nRecvPorousSides.EQ.0) CYCLE

    iPos=0
    DO iSurfSide = 1,SurfMapping(iProc)%nRecvPorousSides
      SideID     = SurfMapping(iProc)%RecvPorousGlobalID(iSurfSide)
      SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)
      PorousSideID = MapSurfSideToPorousSide_Shared(SurfSideID)
      PorousBCSampWall_Shared(:,PorousSideID) = PorousBCSampWall_Shared(:,PorousSideID) &
                                              + PorousBCRecvBuf(iProc)%content(iPos+1:iPos+nValues)
      iPos = iPos + nValues
    END DO ! iSurfSide = 1,SurfMapping(iProc)%nRecvPorousSides

     ! Nullify buffer
    PorousBCRecvBuf(iProc)%content = 0.
  END DO ! iProc
END IF

! ensure synchronization on compute node
CALL MPI_WIN_SYNC(PorousBCSampWall_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)

END SUBROUTINE ExchangeImpingedPartPorousBC


SUBROUTINE ExchangeRemovalProbabilityPorousBC
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED,MPI_COMM_LEADERS_SURF
USE MOD_MPI_Shared_Vars         ,ONLY: nSurfLeaders,myComputeNodeRank,mySurfRank
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfOnNode
USE MOD_Particle_Boundary_Vars  ,ONLY: GlobalSide2SurfSide
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfMapping, MapSurfSideToPorousSide_Shared
USE MOD_Particle_Boundary_Vars  ,ONLY: PorousBCProperties_Shared,PorousBCProperties_Shared_Win
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                             :: iProc,SideID,iPos
INTEGER                             :: MessageSize,iSurfSide,SurfSideID,PorousSideID
INTEGER                             :: nValues
INTEGER                             :: RecvRequest(0:nSurfLeaders-1),SendRequest(0:nSurfLeaders-1)
TYPE tTempArrayProc
  REAL, ALLOCATABLE                 :: SendMsg(:)
  REAL, ALLOCATABLE                 :: RecvMsg(:)
  INTEGER                           :: nRecvRemoveProbSides
  INTEGER                           :: nSendRemoveProbSides
END TYPE
TYPE(tTempArrayProc), ALLOCATABLE   :: TempArrayProc(:)
!===================================================================================================================================
! nodes without sampling surfaces do not take part in this routine
IF (.NOT.SurfOnNode) RETURN

! Synchronize the removal probability
CALL MPI_WIN_SYNC(PorousBCProperties_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)

IF (myComputeNodeRank.EQ.0) THEN
  ! 1) Communication of the removal probability to the halo cell of other procs (since the removal probability requires the locally
  ! sampled values of the adjacent elems) for the case, when a particle leaves the local domain and hits the porous BC on a haloside
  ALLOCATE(TempArrayProc(0:nSurfLeaders-1))
  DO iProc=0, nSurfLeaders-1
    TempArrayProc(iProc)%nRecvRemoveProbSides = SurfMapping(iProc)%nSendPorousSides
    TempArrayProc(iProc)%nSendRemoveProbSides = SurfMapping(iProc)%nRecvPorousSides
    ALLOCATE(TempArrayProc(iProc)%SendMsg(1:TempArrayProc(iProc)%nSendRemoveProbSides))
    ALLOCATE(TempArrayProc(iProc)%RecvMsg(1:TempArrayProc(iProc)%nRecvRemoveProbSides))
    TempArrayProc(iProc)%SendMsg(1:TempArrayProc(iProc)%nSendRemoveProbSides) = 0.
    TempArrayProc(iProc)%RecvMsg(1:TempArrayProc(iProc)%nRecvRemoveProbSides) = 0.
  END DO
  nValues = 1
  ! open receive buffer
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    ! Only open recv buffer if we are expecting sides from this leader node
    IF (TempArrayProc(iProc)%nRecvRemoveProbSides.EQ.0) CYCLE

    ! Message is sent on MPI_COMM_LEADERS_SURF, so rank is indeed iProc
    MessageSize = TempArrayProc(iProc)%nRecvRemoveProbSides * nValues
    CALL MPI_IRECV( TempArrayProc(iProc)%RecvMsg                   &
                  , MessageSize                                  &
                  , MPI_DOUBLE_PRECISION                         &
                  , iProc                                        &
                  , 1209                                         &
                  , MPI_COMM_LEADERS_SURF                        &
                  , RecvRequest(iProc)                           &
                  , IERROR)
  END DO ! iProc

  ! build message
  DO iProc = 0,nSurfLeaders-1
    ! Ignore myself
    IF (iProc .EQ. mySurfRank) CYCLE

    ! Only assemble message if we are expecting sides to send to this leader node
    IF (TempArrayProc(iProc)%nSendRemoveProbSides.EQ.0) CYCLE

    ! Nullify everything
    iPos = 0

    DO iSurfSide = 1,TempArrayProc(iProc)%nSendRemoveProbSides
      ! Get the right side id through the receive global id mapping
      SideID     = SurfMapping(iProc)%RecvPorousGlobalID(iSurfSide)
      SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)
      PorousSideID = MapSurfSideToPorousSide_Shared(SurfSideID)
      ! Assemble message
      TempArrayProc(iProc)%SendMsg(iPos+1:iPos+nValues) = PorousBCProperties_Shared(1:nValues,PorousSideID)
      iPos = iPos + nValues
      PorousBCProperties_Shared(1:nValues,PorousSideID)=0.
    END DO ! iSurfSide=1,TempArrayProc(iProc)%nSendRemoveProbSides
  END DO

  ! send message
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    ! Only open recv buffer if we are expecting sides from this leader node
    IF (TempArrayProc(iProc)%nSendRemoveProbSides.EQ.0) CYCLE

    ! Message is sent on MPI_COMM_LEADERS_SURF, so rank is indeed iProc
    MessageSize = TempArrayProc(iProc)%nSendRemoveProbSides * nValues
    CALL MPI_ISEND( TempArrayProc(iProc)%SendMsg                 &
                  , MessageSize                                  &
                  , MPI_DOUBLE_PRECISION                         &
                  , iProc                                        &
                  , 1209                                         &
                  , MPI_COMM_LEADERS_SURF                        &
                  , SendRequest(iProc)                           &
                  , IERROR)
  END DO ! iProc

  ! Finish received number of sampling surfaces
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    IF (TempArrayProc(iProc)%nSendRemoveProbSides.NE.0) THEN
      CALL MPI_WAIT(SendRequest(iProc),MPIStatus,IERROR)
      IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error',IERROR)
    END IF

    IF (TempArrayProc(iProc)%nRecvRemoveProbSides.NE.0) THEN
      CALL MPI_WAIT(RecvRequest(iProc),MPIStatus,IERROR)
      IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error',IERROR)
    END IF
  END DO ! iProc

  ! add data do my list
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    ! Only open recv buffer if we are expecting sides from this leader node
    IF (TempArrayProc(iProc)%nRecvRemoveProbSides.EQ.0) CYCLE

    iPos=0
    DO iSurfSide = 1,TempArrayProc(iProc)%nRecvRemoveProbSides
      SideID     = SurfMapping(iProc)%SendPorousGlobalID(iSurfSide)
      SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)
      PorousSideID = MapSurfSideToPorousSide_Shared(SurfSideID)
      PorousBCProperties_Shared(1:nValues,PorousSideID) = TempArrayProc(iProc)%RecvMsg(iPos+1:iPos+nValues)
      iPos = iPos + nValues
    END DO ! iSurfSide = 1,TempArrayProc(iProc)%nRecvRemoveProbSides

     ! Nullify buffer
    TempArrayProc(iProc)%RecvMsg = 0.
  END DO ! iProc

  DEALLOCATE(TempArrayProc)
END IF

! ensure synchronization on compute node
CALL MPI_WIN_SYNC(PorousBCProperties_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)

END SUBROUTINE ExchangeRemovalProbabilityPorousBC


SUBROUTINE InitPorousCommunication()
!----------------------------------------------------------------------------------------------------------------------------------!
!
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_SURF
USE MOD_MPI_Shared_Vars         ,ONLY: myLeaderGroupRank,nLeaderGroupProcs
USE MOD_MPI_Shared_Vars         ,ONLY: MPIRankSurfLeader
USE MOD_MPI_Shared_Vars         ,ONLY: mySurfRank,nSurfLeaders
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfOnNode,nComputeNodeSurfTotalSides
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfMapping
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfSide2GlobalSide
USE MOD_Particle_Boundary_Vars  ,ONLY: nPorousSides, PorousBCInfo_Shared
USE MOD_Particle_MPI_Vars       ,ONLY: PorousBCSendBuf,PorousBCRecvBuf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: msg_status(1:MPI_STATUS_SIZE)
INTEGER                       :: iProc
INTEGER                       :: LeaderID
INTEGER                       :: iSide, iPorousSide
INTEGER                       :: nSendPorousSidesTmp(0:nSurfLeaders-1)
INTEGER                       :: nRecvPorousSidesTmp(0:nSurfLeaders-1)
!INTEGER                       :: nSurfSidesLeader(1:2,0:nLeaderGroupProcs-1)
INTEGER                       :: RecvRequest(0:nSurfLeaders-1),SendRequest(0:nSurfLeaders-1)
INTEGER                       :: SendPorousGlobalID(0:nSurfLeaders-1,1:nComputeNodeSurfTotalSides)
INTEGER                       :: SampSizeAllocate
!===================================================================================================================================

IF (.NOT.SurfOnNode) RETURN

nRecvPorousSidesTmp = 0

!--- Open receive buffer (number of sampling surfaces in other node's halo region)
DO iProc = 0,nSurfLeaders-1
  IF (iProc.EQ.mySurfRank) CYCLE

  CALL MPI_IRECV( nRecvPorousSidesTmp(iProc)                                    &
                , 1                                                           &
                , MPI_INTEGER                                                 &
                , iProc                                                       &
                , 1211                                                        &
                , MPI_COMM_LEADERS_SURF                                     &
                , RecvRequest(iProc)                                          &
                , IERROR)
END DO

!--- count all surf sides per other compute-node which get sampling data from current leader
nSendPorousSidesTmp = 0

DO iPorousSide = 1, nPorousSides
  ! count surf sides per compute node
  iSide = PorousBCInfo_Shared(2,iPorousSide)
  LeaderID = SurfSide2GlobalSide(SURF_LEADER,iSide)
  nSendPorousSidesTmp(LeaderID) = nSendPorousSidesTmp(LeaderID) + 1
  SendPorousGlobalID(LeaderID,nSendPorousSidesTmp(LeaderID)) = SurfSide2GlobalSide(SURF_SIDEID,iSide)
END DO

!--- send all other leaders the number of sampling sides coming from current node
DO iProc = 0,nSurfLeaders-1
  IF (iProc.EQ.mySurfRank) CYCLE

  CALL MPI_ISEND( nSendPorousSidesTmp(iProc)                                    &
                , 1                                                           &
                , MPI_INTEGER                                                 &
                , iProc                                                       &
                , 1211                                                        &
                , MPI_COMM_LEADERS_SURF                                     &
                , SendRequest(iProc)                                          &
                , IERROR)
END DO

!--- Finish communication
DO iProc = 0,nSurfLeaders-1
  IF (iProc.EQ.mySurfRank) CYCLE

  CALL MPI_WAIT(SendRequest(iProc),MPISTATUS,IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  CALL MPI_WAIT(RecvRequest(iProc),MPISTATUS,IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO

DO iProc = 0,nSurfLeaders-1
  ! Ignore myself
  IF (iProc .EQ. mySurfRank) CYCLE

  ! Save number of send and recv sides
  SurfMapping(iProc)%nRecvPorousSides = nRecvPorousSidesTmp(MPIRankSurfLeader(iProc))
  SurfMapping(iProc)%nSendPorousSides = nSendPorousSidesTmp(MPIRankSurfLeader(iProc))

  ! Only open recv buffer if we are expecting sides from this leader node
  IF (nRecvPorousSidesTmp(MPIRankSurfLeader(iProc)).EQ.0) CYCLE

  ALLOCATE(SurfMapping(iProc)%RecvPorousGlobalID(1:nRecvPorousSidesTmp(MPIRankSurfLeader(iProc))))

  CALL MPI_IRECV( SurfMapping(iProc)%RecvPorousGlobalID                         &
                , nRecvPorousSidesTmp(MPIRankSurfLeader(iProc))                 &
                , MPI_INTEGER                                                 &
                , iProc                                                       &
                , 1211                                                        &
                , MPI_COMM_LEADERS_SURF                                       &
                , RecvRequest(iProc)                                          &
                , IERROR)
END DO

DO iProc = 0,nSurfLeaders-1
  ! Ignore myself
  IF (iProc .EQ. mySurfRank) CYCLE

  ! Only open send buffer if we are expecting sides from this leader node
  IF (nSendPorousSidesTmp(MPIRankSurfLeader(iProc)).EQ.0) CYCLE

  ALLOCATE(SurfMapping(iProc)%SendPorousGlobalID(1:nSendPorousSidesTmp(MPIRankSurfLeader(iProc))))

  SurfMapping(iProc)%SendPorousGlobalID = SendPorousGlobalID(MPIRankSurfLeader(iProc),1:nSendPorousSidesTmp(MPIRankSurfLeader(iProc)))

  CALL MPI_ISEND( SurfMapping(iProc)%SendPorousGlobalID                         &
                , nSendPorousSidesTmp(MPIRankSurfLeader(iProc))                 &
                , MPI_INTEGER                                                 &
                , iProc                                                       &
                , 1211                                                        &
                , MPI_COMM_LEADERS_SURF                                       &
                , SendRequest(iProc)                                          &
                , IERROR)
END DO

!--- Finish communication
DO iProc = 0,nSurfLeaders-1
  ! Ignore myself
  IF (iProc .EQ. mySurfRank) CYCLE

  IF (nSendPorousSidesTmp(MPIRankSurfLeader(iProc)).NE.0) THEN
    CALL MPI_WAIT(SendRequest(iProc),msg_status(:),IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END IF

  IF (nRecvPorousSidesTmp(MPIRankSurfLeader(iProc)).NE.0) THEN
    CALL MPI_WAIT(RecvRequest(iProc),msg_status(:),IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END IF
END DO

!--- Allocate send and recv buffer for each surf leader
ALLOCATE(PorousBCSendBuf(0:nSurfLeaders-1))
ALLOCATE(PorousBCRecvBuf(0:nSurfLeaders-1))

DO iProc = 0,nSurfLeaders-1
  ! Get message size
  SampSizeAllocate = 2

  ! Only allocate send buffer if we are expecting sides from this leader node
  IF (SurfMapping(iProc)%nSendPorousSides.GT.0) THEN
    ALLOCATE(PorousBCSendBuf(iProc)%content(SampSizeAllocate*SurfMapping(iProc)%nSendPorousSides))
    PorousBCSendBuf(iProc)%content = 0.
  END IF

  ! Only allocate recv buffer if we are expecting sides from this leader node
  IF (SurfMapping(iProc)%nRecvPorousSides.GT.0) THEN
    ALLOCATE(PorousBCRecvBuf(iProc)%content(SampSizeAllocate*SurfMapping(iProc)%nRecvPorousSides))
    PorousBCRecvBuf(iProc)%content = 0.
  END IF
END DO ! iProc

END SUBROUTINE InitPorousCommunication
#endif /*USE_MPI*/


SUBROUTINE GetRadialDistance2D(GlobalSideID,dir,origin,rmin,rmax)
!===================================================================================================================================
! Determines the radial distance to a given origin on a surface
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Surfaces,      ONLY: GetSideBoundingBox
USE MOD_Particle_Mesh_Tools       ,ONLY: GetSideBoundingBoxTria
USE MOD_Particle_Tracking_Vars    ,ONLY: TriaTracking
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: GlobalSideID, dir(3)
REAL, INTENT(IN)              :: origin(2)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: rmin,rmax
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iNode
REAL                          :: BoundingBox(1:3,1:8), point(2), vec(2)
REAL                          :: Vector1(3),Vector2(3),Vector3(3),xyzNod(3),corner(3),VecBoundingBox(3),radiusCorner(2,4)
LOGICAL                       :: r0inside
!===================================================================================================================================
! Determine which cells are inside/outside/partially inside the defined region
IF (TriaTracking) THEN
  CALL GetSideBoundingBoxTria(GlobalSideID,BoundingBox)
ELSE
  CALL GetSideBoundingBox(GlobalSideID,BoundingBox)
END IF
r0inside=.FALSE.
Vector1(:)=0.
Vector2(:)=0.
Vector3(:)=0.
xyzNod(1)=MINVAL(BoundingBox(1,:))
xyzNod(2)=MINVAL(BoundingBox(2,:))
xyzNod(3)=MINVAL(BoundingBox(3,:))
VecBoundingBox(1) = MAXVAL(BoundingBox(1,:)) -MINVAL(BoundingBox(1,:))
VecBoundingBox(2) = MAXVAL(BoundingBox(2,:)) -MINVAL(BoundingBox(2,:))
VecBoundingBox(3) = MAXVAL(BoundingBox(3,:)) -MINVAL(BoundingBox(3,:))
Vector1(dir(2)) = VecBoundingBox(dir(2))
Vector2(dir(2)) = VecBoundingBox(dir(2))
Vector2(dir(3)) = VecBoundingBox(dir(3))
Vector3(dir(3)) = VecBoundingBox(dir(3))
!-- determine rmax (and corners)
DO iNode=1,4
  SELECT CASE(iNode)
  CASE(1)
    corner = xyzNod
  CASE(2)
    corner = xyzNod + Vector1
  CASE(3)
    corner = xyzNod + Vector2
  CASE(4)
    corner = xyzNod + Vector3
  END SELECT
  corner(dir(2)) = corner(dir(2)) - origin(1)
  corner(dir(3)) = corner(dir(3)) - origin(2)
  radiusCorner(1,iNode)=SQRT(corner(dir(2))**2+corner(dir(3))**2)
END DO !iNode
rmax=MAXVAL(radiusCorner(1,1:4))
!-- determine rmin
DO iNode=1,4
  SELECT CASE(iNode)
  CASE(1)
    point=(/xyzNod(dir(2)),xyzNod(dir(3))/)-origin
    vec=(/Vector1(dir(2)),Vector1(dir(3))/)
  CASE(2)
    point=(/xyzNod(dir(2)),xyzNod(dir(3))/)-origin
    vec=(/Vector3(dir(2)),Vector3(dir(3))/)
  CASE(3)
    point=(/xyzNod(dir(2)),xyzNod(dir(3))/)+(/Vector2(dir(2)),Vector2(dir(3))/)-origin
    vec=(/-Vector1(dir(2)),-Vector1(dir(3))/)
  CASE(4)
    point=(/xyzNod(dir(2)),xyzNod(dir(3))/)+(/Vector2(dir(2)),Vector2(dir(3))/)-origin
    vec=(/-Vector3(dir(2)),-Vector3(dir(3))/)
  END SELECT
  vec=point + MIN(MAX(-DOT_PRODUCT(point,vec)/DOT_PRODUCT(vec,vec),0.),1.)*vec
  radiusCorner(2,iNode)=SQRT(DOT_PRODUCT(vec,vec)) !rmin
END DO !iNode
!-- determine if r0 is inside of bounding box
IF ((origin(1) .GE. MINVAL(BoundingBox(dir(2),:))) .AND. &
    (origin(1) .LE. MAXVAL(BoundingBox(dir(2),:))) .AND. &
    (origin(2) .GE. MINVAL(BoundingBox(dir(3),:))) .AND. &
    (origin(2) .LE. MAXVAL(BoundingBox(dir(3),:))) ) THEN
    r0inside = .TRUE.
END IF
IF (r0inside) THEN
  rmin = 0.
ELSE
  rmin=MINVAL(radiusCorner(2,1:4))
END IF

END SUBROUTINE GetRadialDistance2D


SUBROUTINE FinalizePorousBoundaryCondition()
!===================================================================================================================================
!> Deallocates
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Boundary_Vars
#if USE_MPI
USE MOD_Particle_MPI_Vars       ,ONLY: PorousBCSendBuf,PorousBCRecvBuf
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(PorousBCMacroVal)
SDEALLOCATE(PorousBC)
SDEALLOCATE(PorousBCSampWall)
SDEALLOCATE(PorousBCOutput)

! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
CALL MPI_WIN_UNLOCK_ALL(MapSurfSideToPorousSide_Shared_Win,iError)
CALL MPI_WIN_FREE(MapSurfSideToPorousSide_Shared_Win,iError)
CALL MPI_WIN_UNLOCK_ALL(PorousBCInfo_Shared_Win,iError)
CALL MPI_WIN_FREE(PorousBCInfo_Shared_Win,iError)
CALL MPI_WIN_UNLOCK_ALL(PorousBCProperties_Shared_Win,iError)
CALL MPI_WIN_FREE(PorousBCProperties_Shared_Win,iError)
CALL MPI_WIN_UNLOCK_ALL(PorousBCSampWall_Shared_Win,iError)
CALL MPI_WIN_FREE(PorousBCSampWall_Shared_Win,iError)

SDEALLOCATE(PorousBCSendBuf)
SDEALLOCATE(PorousBCRecvBuf)
#endif

ADEALLOCATE(MapSurfSideToPorousSide_Shared)
ADEALLOCATE(PorousBCInfo_Shared)
ADEALLOCATE(PorousBCProperties_Shared)
ADEALLOCATE(PorousBCSampWall_Shared)

END SUBROUTINE FinalizePorousBoundaryCondition

END MODULE MOD_Particle_Boundary_Porous