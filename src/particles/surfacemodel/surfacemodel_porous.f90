!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_SurfaceModel_Porous
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
PUBLIC::DefineParametersPorousBC
PUBLIC::InitPorousBoundaryCondition
PUBLIC::PorousBoundaryTreatment
PUBLIC::PorousBoundaryRemovalProb_Pressure
PUBLIC::FinalizePorousBoundaryCondition
#if USE_MPI
PUBLIC::InitPorousCommunication
#endif /*USE_MPI*/
!===================================================================================================================================

CONTAINS

SUBROUTINE DefineParametersPorousBC()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Porous BC")

CALL prms%CreateIntOption(      'Surf-nPorousBC','Number of porous boundary conditions', '0')
CALL prms%CreateIntOption(      'Surf-PorousBC[$]-BC','PartBound to be a porous boundary', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surf-PorousBC[$]-Pressure','Pressure [Pa] at the porous boundary', numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Surf-PorousBC[$]-Type' &
                                , 'Define the type of porous boundary, currently available are sensor and pump.' //&
                                  'Option sensor: Using the defined region/BC to measure the pressure difference between the '//&
                                  'defined pressure and the pressure at the sensor. Option pump: Use the defined region/BC as '//&
                                  'as a pump by supplying a constant -PumpingSpeed or by supplying -DeltaPumpingSpeed-Kp ' //&
                                  'and -DeltaPumpingSpeed-Ki to adapt the pumping speed until the target pressure is reached.' &
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surf-PorousBC[$]-PumpingSpeed' &
                                , 'PorousBC used as a pump: Removal probability at the BC is determined through the pumping ' //&
                                  'speed, which can be controlled to achieve the given pressure. Initial pumping speed [m3/s] ' //&
                                  'can be zero', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surf-PorousBC[$]-DeltaPumpingSpeed-Kp' &
                                , 'Proportional factor for the pumping speed controller, Kp=0 -> constant given pumping speed' &
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surf-PorousBC[$]-DeltaPumpingSpeed-Ki' &
                                , 'Integral factor for the pumping speed controller, Ki=0 -> only proportional part is used' &
                                , numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Surf-PorousBC[$]-Region' &
                                , 'Define a region, allowing small geometrical features such as holes w/o meshing. ' //&
                                  'Option: circular, requires the definition of the parameters -normalDir, origin, rmax/rmin' &
                                , 'none', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Surf-PorousBC[$]-normalDir' &
                                , 'Definition of the normal direction (as an integer, 1: x, 2: y, 3: z) of the boundary on ' //&
                                  'which the region shall be created.' &
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Surf-PorousBC[$]-origin' &
                                , 'Coordinates of the middle point of the region, Example: normalDir=1: (/y,z/), ' //&
                                  'normalDir=2: (/z,x/), normalDir=3: (/x,y/)', numberedmulti=.TRUE., no=2)
CALL prms%CreateRealOption(     'Surf-PorousBC[$]-rmax' &
                                , 'Maximum radius [m] of the circular region', '1e21', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surf-PorousBC[$]-rmin' &
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
USE MOD_Particle_Boundary_Tools     ,ONLY: GetRadialDistance2D
USE MOD_Particle_Mesh_Vars          ,ONLY: SideInfo_Shared
USE MOD_Particle_Vars               ,ONLY: Symmetry, VarTimeStep
USE MOD_SurfaceModel_Vars           ,ONLY: nPorousBC, PorousBC
USE MOD_Particle_Boundary_Vars      ,ONLY: PartBound, nPorousSides, MapSurfSideToPorousSide_Shared, PorousBCSampWall
USE MOD_Particle_Boundary_Vars      ,ONLY: PorousBCInfo_Shared,PorousBCProperties_Shared,PorousBCSampWall_Shared
USE MOD_Particle_Boundary_Vars      ,ONLY: nComputeNodeSurfTotalSides, SurfSide2GlobalSide
USE MOD_Particle_Tracking_Vars      ,ONLY: TrackingMethod
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars             ,ONLY: MPI_COMM_SHARED, myComputeNodeRank
USE MOD_Particle_Boundary_Vars      ,ONLY: MapSurfSideToPorousSide_Shared_Win
USE MOD_Particle_Boundary_Vars      ,ONLY: PorousBCInfo_Shared_Win,PorousBCProperties_Shared_Win,PorousBCSampWall_Shared_Win
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars            ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
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
!===================================================================================================================================

LBWRITE(UNIT_stdOut,'(A)') ' INIT POROUS BOUNDARY CONDITION ...'

IF(TrackingMethod.EQ.REFMAPPING) CALL abort(__STAMP__,'ERROR: Porous boundary conditions are not implemented with RefMapping!')

IF(VarTimeStep%UseSpeciesSpecific) CALL abort(__STAMP__,'ERROR: Porous boundary conditions are not implemented with a species-specific time step!')

IF((Symmetry%Order.LE.2).AND.(.NOT.Symmetry%Axisymmetric)) CALL abort(__STAMP__,'ERROR: Porous boundary conditions are not implemented for 1D/2D simulations!')

! 1) Allocate arrays
ALLOCATE(PorousBC(1:nPorousBC))
ALLOCATE(MapSurfSideToPorousBC_Temp(1:2,1:nComputeNodeSurfTotalSides))
MapSurfSideToPorousBC_Temp = 0

#if USE_MPI
CALL Allocate_Shared((/nComputeNodeSurfTotalSides/),MapSurfSideToPorousSide_Shared_Win,MapSurfSideToPorousSide_Shared)
CALL MPI_WIN_LOCK_ALL(0,MapSurfSideToPorousSide_Shared_Win,IERROR)
IF (myComputeNodeRank.EQ.0) THEN
  MapSurfSideToPorousSide_Shared = 0
END IF
! This barrier MIGHT not be required
CALL BARRIER_AND_SYNC(MapSurfSideToPorousSide_Shared_Win,MPI_COMM_SHARED)
#else
ALLOCATE(MapSurfSideToPorousSide_Shared(1:nComputeNodeSurfTotalSides))
MapSurfSideToPorousSide_Shared = 0
#endif /*USE_MPI*/

DO iPBC = 1, nPorousBC
  WRITE(UNIT=hilf,FMT='(I0)') iPBC
  PorousBC(iPBC)%BC = GETINT('Surf-PorousBC'//TRIM(hilf)//'-BC')
  IF(PartBound%TargetBoundCond(PorousBC(iPBC)%BC).NE.PartBound%ReflectiveBC) THEN
    CALL abort(__STAMP__,'ERROR in init of porous BC: given boundary condition must be reflective!')
  END IF
  ! Read-in of the conditions at the porous boundary
  PorousBC(iPBC)%Type  = TRIM(GETSTR('Surf-PorousBC'//TRIM(hilf)//'-Type'))
  SELECT CASE(PorousBC(iPBC)%Type)
    CASE('sensor')
      PorousBC(iPBC)%PumpingSpeed = 0.0
      PorousBC(iPBC)%DeltaPumpingSpeedKp = 0.0
      PorousBC(iPBC)%DeltaPumpingSpeedKi = 0.0
    CASE('pump')
      ! Define the desired pressure at the pump
      PorousBC(iPBC)%Pressure  = GETREAL('Surf-PorousBC'//TRIM(hilf)//'-Pressure')
      ! Initial pumping speed at the porous boundary [m3/s]
      PorousBC(iPBC)%PumpingSpeed = GETREAL('Surf-PorousBC'//TRIM(hilf)//'-PumpingSpeed')
      ! Proportional and integral factors for the control of the pumping speed
      PorousBC(iPBC)%DeltaPumpingSpeedKp = GETREAL('Surf-PorousBC'//TRIM(hilf)//'-DeltaPumpingSpeed-Kp')
      PorousBC(iPBC)%DeltaPumpingSpeedKi = GETREAL('Surf-PorousBC'//TRIM(hilf)//'-DeltaPumpingSpeed-Ki')
      ! Determining the order of magnitude
      IF(PorousBC(iPBC)%Pressure.GT.0.0) THEN
        PorousBC(iPBC)%DeltaPumpingSpeedKp = PorousBC(iPBC)%DeltaPumpingSpeedKp / 10.0**(ANINT(LOG10(PorousBC(iPBC)%Pressure)))
        PorousBC(iPBC)%DeltaPumpingSpeedKi = PorousBC(iPBC)%DeltaPumpingSpeedKi / 10.0**(ANINT(LOG10(PorousBC(iPBC)%Pressure)))
      END IF
    CASE DEFAULT
      CALL abort(__STAMP__,'ERROR in type definition of porous bc:', iPBC)
  END SELECT
  PorousBC(iPBC)%Region = TRIM(GETSTR('Surf-PorousBC'//TRIM(hilf)//'-Region','none'))
  IF(PorousBC(iPBC)%Region.EQ.'none') THEN
    PorousBC(iPBC)%UsingRegion = .FALSE.
  ELSE
    PorousBC(iPBC)%UsingRegion = .TRUE.
    ! Normal direction to the surface (could be replaced by the normal vector of the boundary?)
    PorousBC(iPBC)%dir(1) = GETINT('Surf-PorousBC'//TRIM(hilf)//'-normalDir')
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
      CALL abort(__STAMP__,'ERROR in init: normalDir for PorousBC must be between 1 and 3!')
    END IF
    SELECT CASE(PorousBC(iPBC)%Region)
      CASE('circular')
        PorousBC(iPBC)%origin = GETREALARRAY('Surf-PorousBC'//TRIM(hilf)//'-origin',2)
        WRITE(UNIT=hilf2,FMT='(E16.8)') HUGE(PorousBC(iPBC)%rmax)
        PorousBC(iPBC)%rmax = GETREAL('Surf-PorousBC'//TRIM(hilf)//'-rmax',TRIM(hilf2))
        PorousBC(iPBC)%rmin = GETREAL('Surf-PorousBC'//TRIM(hilf)//'-rmin','0.')
        IF(Symmetry%Axisymmetric) THEN
          IF(PorousBC(iPBC)%dir(1).NE.1) THEN
            CALL abort(__STAMP__,'ERROR in Porous BC: For axisymmetric simulations, only regions perpendicular to the axis are allowed!', iPBC)
          END IF
          IF(PorousBC(iPBC)%origin(1)*PorousBC(iPBC)%origin(2).NE.0.0) THEN
            CALL abort(__STAMP__,'ERROR in Porous BC: For axisymmetric simulations, the origin has to be at (0,0)!', iPBC)
          END IF
        END IF
      CASE DEFAULT
        CALL abort(__STAMP__&
        ,'ERROR in region definition of porous bc:', iPBC)
    END SELECT
  END IF    ! Region is given
END DO      ! nPorousBC

! 2) Mapping of the porous BC sides to the respective surface side
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
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
#if USE_MPI
END IF

CALL MPI_BCAST(nPorousSides,1,MPI_INTEGER,0,MPI_COMM_SHARED,IERROR)

CALL Allocate_Shared((/3,nPorousSides/),PorousBCInfo_Shared_Win,PorousBCInfo_Shared)
CALL Allocate_Shared((/2,nPorousSides/),PorousBCProperties_Shared_Win,PorousBCProperties_Shared)
CALL Allocate_Shared((/2,nPorousSides/),PorousBCSampWall_Shared_Win,PorousBCSampWall_Shared)
CALL MPI_WIN_LOCK_ALL(0,PorousBCInfo_Shared_Win,IERROR)
CALL MPI_WIN_LOCK_ALL(0,PorousBCProperties_Shared_Win,IERROR)
CALL MPI_WIN_LOCK_ALL(0,PorousBCSampWall_Shared_Win,IERROR)
IF (myComputeNodeRank.EQ.0) THEN
  PorousBCInfo_Shared(1,1:nPorousSides) = MapSurfSideToPorousBC_Temp(1,1:nPorousSides)
  PorousBCInfo_Shared(2,1:nPorousSides) = MapSurfSideToPorousBC_Temp(2,1:nPorousSides)
  PorousBCInfo_Shared(3,1:nPorousSides) = -1
  PorousBCProperties_Shared = 0.
  PorousBCSampWall_Shared = 0.
END IF
! This barrier MIGHT not be required
CALL BARRIER_AND_SYNC(PorousBCInfo_Shared_Win      ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(PorousBCProperties_Shared_Win,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(PorousBCSampWall_Shared_Win  ,MPI_COMM_SHARED)
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
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
  ! 4) Mapping the porous BC ID and the porous BC side ID to a surface side
  DO iPorousSide=1,nPorousSides
    iPBC = PorousBCInfo_Shared(1,iPorousSide)
    PorousBCProperties_Shared(2,iPorousSide) = PorousBC(iPBC)%PumpingSpeed
    GlobalSideID = SurfSide2GlobalSide(SURF_SIDEID,PorousBCInfo_Shared(2,iPorousSide))
    ! If a pumping speed was read-in during restart, use it (only for local elements and not halo sides)
    ! NOTE: DOES NOT WORK CURRENTLY AS ONLY NODE-LEADER IS PERFORMING THIS LOOP, HOWEVER ADAPTIVEMACROVAL IS STILL PROC LOCAL
    ! IF(AdaptBCMacroVal(5,SideInfo_Shared(SIDE_ELEMID,GlobalSideID)-offsetElem,1).GT.0.0) THEN
    !     PorousBCProperties_Shared(2,iPorousSide) = AdaptBCMacroVal(5,SideInfo_Shared(SIDE_ELEMID,GlobalSideID)-offsetElem,1)
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
#if USE_MPI
END IF

CALL BARRIER_AND_SYNC(PorousBCInfo_Shared_Win      ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(PorousBCProperties_Shared_Win,MPI_COMM_SHARED)
IF (myComputeNodeRank.EQ.0) CALL InitPorousCommunication()
#endif

SDEALLOCATE(MapSurfSideToPorousBC_Temp)

LBWRITE(UNIT_stdOut,'(A)') ' INIT POROUS BOUNDARY CONDITION DONE!'

END SUBROUTINE InitPorousBoundaryCondition

SUBROUTINE PorousBoundaryTreatment(iPart,SideID,ElasticReflectionAtPorousBC)
!===================================================================================================================================
! Treatment of particles impinging on the porous boundary
! 1) (Optional) When using regions on the BC, it is determined whether the particle hit the porous BC region or only the regular BC
! 2) Comparison of the removal probability with a random number to determine whether the particle is deleted. Counting number of
!    impinged and deleted particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars           ,ONLY: LastPartPos
USE MOD_SurfaceModel_Vars       ,ONLY: PorousBC
USE MOD_Particle_Boundary_Vars  ,ONLY: PorousBCSampWall, MapSurfSideToPorousSide_Shared
USE MOD_Particle_Boundary_Vars  ,ONLY: PorousBCInfo_Shared, PorousBCProperties_Shared
USE MOD_part_tools              ,ONLY: GetParticleWeight
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound, GlobalSide2SurfSide
USE MOD_part_operations         ,ONLY: RemoveParticle
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackInfo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPart, SideID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
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
      intersectionPoint(1:3) = LastPartPos(1:3,iPart) + TrackInfo%alpha*TrackInfo%PartTrajectory(1:3)
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
          CALL RemoveParticle(iPart,BCID=PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)))
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
! 2) Compute the total number of impinged particles per porous BC on each surface comm leader (=compute node leader)
! 3) Sum-up the number of impinged particles across whole simulation domain
! 4.1) Loop over all sides: a), b) & d) is only performed if the pumping speed is adapted
!     a) Determining the delta between current gas mixture pressure in adjacent cell and target pressure
!     b) Adapting the pumping capacity (m^3/s) according to pressure difference (control through proportional and integral part)
!     c) Calculate the removal probability if any particles hit the pump
!     d) Limit removal probability to values between 0 and 1
!     e) Sampling of the pumping capacity (and other variables for PartAnalyze) for the output
! 3) MPI communication of the removal probability to halo sides
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars                   ,ONLY: DSMC, RadialWeighting
USE MOD_Mesh_Vars                   ,ONLY: nElems,offsetElem
USE MOD_Particle_Boundary_Vars      ,ONLY: SurfTotalSideOnNode, SurfSide2GlobalSide
USE MOD_SurfaceModel_Vars           ,ONLY: nPorousBC, PorousBC
USE MOD_SurfaceModel_Analyze_Vars   ,ONLY: CalcPorousBCInfo, PorousBCOutput
USE MOD_Particle_Boundary_Vars      ,ONLY: nPorousSides, PorousBCProperties_Shared, PorousBCInfo_Shared, SampWallPumpCapacity
USE MOD_Particle_Boundary_Vars      ,ONLY: PorousBCSampWall, PorousBCSampWall_Shared
USE MOD_Particle_Vars               ,ONLY: Species, nSpecies, usevMPF, UseVarTimeStep
USE MOD_Particle_Sampling_Vars      ,ONLY: AdaptBCMacroVal, AdaptBCMapElemToSample
USE MOD_Particle_TimeStep           ,ONLY: GetParticleTimeStep
USE MOD_Timedisc_Vars               ,ONLY: dt
USE MOD_Particle_Mesh_Vars
USE MOD_Mesh_Tools                  ,ONLY: GetCNElemID
#if USE_MPI
USE MOD_MPI_Shared_Vars             ,ONLY: myComputeNodeRank, MPI_COMM_SHARED,MPI_COMM_LEADERS_SURF
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: LocalElemID,CNElemID,GlobalElemID,SurfSideID,GlobalSideID,iPorousSide,iPBC,SampleElemID
REAL                          :: PumpingSpeedTemp,DeltaPressure,partWeight,SumPartPorousBC,dtVar
REAL                          :: SumPartImpinged(nPorousBC)
!===================================================================================================================================

IF (.NOT.SurfTotalSideOnNode) RETURN

IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
  partWeight = 1.
ELSE
  partWeight = Species(1)%MacroParticleFactor
END IF

SumPartImpinged = 0
! 1) MPI communication of particles that impinged on halo sides to corresponding side
#if USE_MPI
CALL ExchangeImpingedPartPorousBC()
#else
PorousBCSampWall_Shared(1:2,1:nPorousSides) = PorousBCSampWall(1:2,1:nPorousSides)
#endif /*USE_MPI*/

! 2) Compute the total number of impinged particles per porous BC on each surface comm leader (=compute node leader)
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
  DO iPorousSide = 1, nPorousSides
    iPBC = PorousBCInfo_Shared(1,iPorousSide)
    SumPartImpinged(iPBC) = SumPartImpinged(iPBC) + PorousBCSampWall_Shared(1,iPorousSide)
  END DO
#if USE_MPI
END IF
#endif /*USE_MPI*/

! 3) Sum-up the number of impinged particles across whole simulation domain by communication of surface comm leaders and distribute
!     the information to the processors on the node
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,SumPartImpinged(1:nPorousBC),nPorousBC,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEADERS_SURF,iError)
END IF
CALL MPI_BCAST(SumPartImpinged(1:nPorousBC),nPorousBC,MPI_DOUBLE_PRECISION,0,MPI_COMM_SHARED,iError)
#endif /*USE_MPI*/

! Zero the output variable
IF(CalcPorousBCInfo) PorousBCOutput = 0.

! 4) Loop over all porous sides
DO iPorousSide = 1, nPorousSides
  SurfSideID = PorousBCInfo_Shared(2,iPorousSide)
  GlobalSideID = SurfSide2GlobalSide(SURF_SIDEID,SurfSideID)
  GlobalElemID = SideInfo_Shared(SIDE_ELEMID,GlobalSideID)
  CNElemID = GetCNElemID(GlobalElemID)
  ! Only treat your proc-local elements
  IF ((GlobalElemID.LT.1+offSetElem).OR.(GlobalElemID.GT.nElems+offSetElem)) CYCLE
  LocalElemID = GlobalElemID - offsetElem
  SampleElemID = AdaptBCMapElemToSample(LocalElemID)
  iPBC = PorousBCInfo_Shared(1,iPorousSide)
  SumPartPorousBC = SumPartImpinged(iPBC)
  ! Skip element if number density is zero
  IF(SUM(AdaptBCMacroVal(4,SampleElemID,1:nSpecies)).EQ.0.0) CYCLE
  ! Get the correct time step of the cell
  IF(UseVarTimeStep) THEN
    dtVar = dt * GetParticleTimeStep(ElemMidPoint_Shared(1,CNElemID), ElemMidPoint_Shared(2,CNElemID), LocalElemID)
  ELSE
    dtVar = dt
  END IF
  ! Determine the removal probability based on the pumping speed (adaptive to a target pressure or fixed)
  IF(TRIM(PorousBC(iPBC)%Type).EQ.'pump') THEN
    IF((PorousBC(iPBC)%DeltaPumpingSpeedKp.GT.0.).OR.(PorousBC(iPBC)%DeltaPumpingSpeedKi.GT.0.)) THEN
      ! a) Determining the delta between current gas mixture pressure in adjacent cell and target pressure
      DeltaPressure = SUM(AdaptBCMacroVal(6,SampleElemID,1:nSpecies))-PorousBC(iPBC)%Pressure
      ! Integrating the pressure difference (only utilized later if DeltaPumpingSpeedKi was given)
      IF(PorousBCProperties_Shared(2,iPorousSide).GT.0.0) THEN
        AdaptBCMacroVal(7,SampleElemID,1) = AdaptBCMacroVal(7,SampleElemID,1) + DeltaPressure * dtVar
      ELSE
        AdaptBCMacroVal(7,SampleElemID,1) = 0.0
      END IF
      ! b) Adapting the pumping capacity (m^3/s) according to pressure difference (control through proportional and integral part)
      PumpingSpeedTemp = PorousBCProperties_Shared(2,iPorousSide) + PorousBC(iPBC)%DeltaPumpingSpeedKp * DeltaPressure &
          + PorousBC(iPBC)%DeltaPumpingSpeedKi * AdaptBCMacroVal(7,SampleElemID,1)
      ! c) Calculate the removal probability if any particles hit the pump
      IF(SumPartPorousBC.GT.0) THEN
        PorousBCProperties_Shared(1,iPorousSide) = PumpingSpeedTemp*SUM(AdaptBCMacroVal(4,SampleElemID,1:nSpecies)) &
                                                        * dtVar / (SumPartPorousBC*partWeight)
      ELSE
        PorousBCProperties_Shared(1,iPorousSide) = 0.0
      END IF
      ! d) Limit removal probability to values between 0 and 1
      IF(PorousBCProperties_Shared(1,iPorousSide).GT.1.0) THEN
        PorousBCProperties_Shared(1,iPorousSide) = 1.0
        ! Setting pumping speed to maximum value (alpha=1)
        PorousBCProperties_Shared(2,iPorousSide) = SumPartPorousBC*partWeight &
                                                      / (SUM(AdaptBCMacroVal(4,SampleElemID,1:nSpecies))*dtVar)
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
                        * SUM(AdaptBCMacroVal(4,SampleElemID,1:nSpecies)) * dtVar / (SumPartPorousBC*partWeight)
      ELSE
        PorousBCProperties_Shared(1,iPorousSide) = 0.0
      END IF
    END IF
  END IF      ! TRIM(PorousBC(iPBC)%Type).EQ.'pump'
  ! Storing the pumping speed for the restart state file
  AdaptBCMacroVal(5,SampleElemID,1) = PorousBCProperties_Shared(2,iPorousSide)
  ! e) Sampling of the pumping capacity (and other variables for PartAnalyze) for the output
  ! -------- Sampling for output in DSMCSurfState --------------------------------------------------------------------------------
  IF(DSMC%CalcSurfaceVal) THEN
    SampWallPumpCapacity(SurfSideID) = SampWallPumpCapacity(SurfSideID) + PorousBCSampWall_Shared(2,iPorousSide) &
                                                * partWeight / (SUM(AdaptBCMacroVal(4,SampleElemID,1:nSpecies))*dtVar)
  END IF
  ! -------- Sampling for output in PartAnalyze ----------------------------------------------------------------------------------
  IF(CalcPorousBCInfo) THEN
    ! Sampling the actual instantaneous pumping speed S (m^3/s) through the number of deleted particles  (-PumpSpeed-Measure-)
    PorousBCOutput(2,iPBC) = PorousBCOutput(2,iPBC) + PorousBCSampWall_Shared(2,iPorousSide) * PartWeight &
                                                          / (SUM(AdaptBCMacroVal(4,SampleElemID,1:nSpecies))*dtVar)
    IF(PorousBCSampWall_Shared(1,iPorousSide).GT.0) THEN
      ! Counting only sides, where a particle hit the pump
      PorousBCOutput(1,iPBC) = PorousBCOutput(1,iPBC) + 1.
      ! Pumping speed S (m^3/s) used at the pump to calculate the removal probability (-PumpSpeed-Control-)
      PorousBCOutput(3,iPBC) = PorousBCOutput(3,iPBC) + PorousBCProperties_Shared(2,iPorousSide)
      ! Removal probability
      PorousBCOutput(4,iPBC) = PorousBCOutput(4,iPBC) + PorousBCProperties_Shared(1,iPorousSide)
      ! Pressure at the pump
      PorousBCOutput(5,iPBC) = PorousBCOutput(5,iPBC) + SUM(AdaptBCMacroVal(6,SampleElemID,1:nSpecies))
    END IF
  END IF
END DO    ! iPorousSide = 1, nPorousSides

! Reset of the sampled particle numbers at the pump
PorousBCSampWall = 0.

! Exchange removal probability for halo cells
#if USE_MPI
CALL ExchangeRemovalProbabilityPorousBC()
#endif

END SUBROUTINE PorousBoundaryRemovalProb_Pressure


#if USE_MPI
SUBROUTINE ExchangeImpingedPartPorousBC()
!===================================================================================================================================
!> Communication of particles impinged on halo sides, required for the calculation of the removal probability
!> 1.) Collect the information from the proc-local shadow arrays in the compute-node shared array
!> 2.) Internode communication between the surface leader processors of each compute node
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_MPI_Shared              ,ONLY: BARRIER_AND_SYNC
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED,MPI_COMM_LEADERS_SURF
USE MOD_MPI_Shared_Vars         ,ONLY: nSurfLeaders,myComputeNodeRank,mySurfRank
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfTotalSideOnNode
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
IF (.NOT.SurfTotalSideOnNode) RETURN

nValues = 2
! 1.) Collect the information from the proc-local shadow arrays in the compute-node shared array
MessageSize = nValues*nPorousSides
IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_REDUCE(PorousBCSampWall,PorousBCSampWall_Shared,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
ELSE
  CALL MPI_REDUCE(PorousBCSampWall,0                      ,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
END IF
CALL BARRIER_AND_SYNC(PorousBCSampWall_Shared_Win,MPI_COMM_SHARED)

! 2.) Internode communication between the surface leader processors of each compute node
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
    END DO ! 1,SurfMapping(iProc)%nSendPorousSides
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
CALL BARRIER_AND_SYNC(PorousBCSampWall_Shared_Win,MPI_COMM_SHARED)

END SUBROUTINE ExchangeImpingedPartPorousBC


SUBROUTINE ExchangeRemovalProbabilityPorousBC
!===================================================================================================================================
!> Routine that communicates the calculated removal probability to the halo sides (since the removal probability requires locally
!> sampled values of the adjacent elements). The halo sides then have the correct removal probability when treating the impinging
!> particles. The sides that communicate the impinged particles, receive the calculated removal probability.
!> Thus, the nSendPorousSides/PorousBCSendBuf variables are utilized for the RECEIVE buffer and nRecvPorousSides/PorousBCRecvBuf
!> are utilized for the SEND buffer.
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
USE MOD_Globals
USE MOD_MPI_Shared              ,ONLY: BARRIER_AND_SYNC
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED,MPI_COMM_LEADERS_SURF
USE MOD_MPI_Shared_Vars         ,ONLY: nSurfLeaders,myComputeNodeRank,mySurfRank
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfTotalSideOnNode
USE MOD_Particle_Boundary_Vars  ,ONLY: GlobalSide2SurfSide
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfMapping, MapSurfSideToPorousSide_Shared
USE MOD_Particle_Boundary_Vars  ,ONLY: PorousBCProperties_Shared,PorousBCProperties_Shared_Win
USE MOD_Particle_MPI_Vars       ,ONLY: PorousBCSendBuf,PorousBCRecvBuf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                             :: iProc,SideID,iPos
INTEGER                             :: MessageSize,iSurfSide,SurfSideID,PorousSideID
INTEGER                             :: nValues
INTEGER                             :: RecvRequest(0:nSurfLeaders-1),SendRequest(0:nSurfLeaders-1)
!===================================================================================================================================
! nodes without sampling surfaces do not take part in this routine
IF (.NOT.SurfTotalSideOnNode) RETURN

! 1) Synchronize the removal probability
CALL BARRIER_AND_SYNC(PorousBCProperties_Shared_Win,MPI_COMM_SHARED)

IF (myComputeNodeRank.EQ.0) THEN
  ! 2) Internode communication
  nValues = 2
  ! open receive buffer
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    ! Only open recv buffer if we are expecting sides from this leader node
    IF (SurfMapping(iProc)%nSendPorousSides.EQ.0) CYCLE

    ! Message is sent on MPI_COMM_LEADERS_SURF, so rank is indeed iProc
    MessageSize = SurfMapping(iProc)%nSendPorousSides * nValues
    CALL MPI_IRECV( PorousBCSendBuf(iProc)%content                   &
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
    IF (SurfMapping(iProc)%nRecvPorousSides.EQ.0) CYCLE

    ! Nullify everything
    iPos = 0

    DO iSurfSide = 1,SurfMapping(iProc)%nRecvPorousSides
      ! Get the right side id through the receive global id mapping
      SideID     = SurfMapping(iProc)%RecvPorousGlobalID(iSurfSide)
      SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)
      PorousSideID = MapSurfSideToPorousSide_Shared(SurfSideID)
      ! Assemble message
      PorousBCRecvBuf(iProc)%content(iPos+1:iPos+nValues) = PorousBCProperties_Shared(1:nValues,PorousSideID)
      iPos = iPos + nValues
    END DO ! iSurfSide=1,SurfMapping(iProc)%nRecvPorousSides
  END DO

  ! send message
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    ! Only open recv buffer if we are expecting sides from this leader node
    IF (SurfMapping(iProc)%nRecvPorousSides.EQ.0) CYCLE

    ! Message is sent on MPI_COMM_LEADERS_SURF, so rank is indeed iProc
    MessageSize = SurfMapping(iProc)%nRecvPorousSides * nValues
    CALL MPI_ISEND( PorousBCRecvBuf(iProc)%content                 &
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

    IF (SurfMapping(iProc)%nRecvPorousSides.NE.0) THEN
      CALL MPI_WAIT(SendRequest(iProc),MPIStatus,IERROR)
      IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error',IERROR)
    END IF

    IF (SurfMapping(iProc)%nSendPorousSides.NE.0) THEN
      CALL MPI_WAIT(RecvRequest(iProc),MPIStatus,IERROR)
      IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error',IERROR)
    END IF
  END DO ! iProc

  ! add data do my list
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    ! Only open recv buffer if we are expecting sides from this leader node
    IF (SurfMapping(iProc)%nSendPorousSides.EQ.0) CYCLE

    iPos=0
    DO iSurfSide = 1,SurfMapping(iProc)%nSendPorousSides
      SideID     = SurfMapping(iProc)%SendPorousGlobalID(iSurfSide)
      SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)
      PorousSideID = MapSurfSideToPorousSide_Shared(SurfSideID)
      PorousBCProperties_Shared(1:nValues,PorousSideID) = PorousBCSendBuf(iProc)%content(iPos+1:iPos+nValues)
      iPos = iPos + nValues
    END DO ! iSurfSide = 1,SurfMapping(iProc)%nSendPorousSides

     ! Nullify buffer
    PorousBCSendBuf(iProc)%content = 0.
  END DO ! iProc
END IF

! ensure synchronization on compute node
CALL BARRIER_AND_SYNC(PorousBCProperties_Shared_Win,MPI_COMM_SHARED)

END SUBROUTINE ExchangeRemovalProbabilityPorousBC


SUBROUTINE InitPorousCommunication()
!===================================================================================================================================
!> Initialize the communication for the porous BCs: Impinging particles and removal probability.
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
USE MOD_Globals
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_LEADERS_SHARED,MPI_COMM_LEADERS_SURF
USE MOD_MPI_Shared_Vars         ,ONLY: MPIRankSurfLeader
USE MOD_MPI_Shared_Vars         ,ONLY: myLeaderGroupRank,nSurfLeaders,nLeaderGroupProcs
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfTotalSideOnNode,nComputeNodeSurfTotalSides
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfMapping
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfSide2GlobalSide
USE MOD_Particle_Boundary_Vars  ,ONLY: nPorousSides, PorousBCInfo_Shared
USE MOD_Particle_MPI_Vars       ,ONLY: PorousBCSendBuf,PorousBCRecvBuf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: msg_status(1:MPI_STATUS_SIZE)
INTEGER                       :: iProc
INTEGER                       :: LeaderID
INTEGER                       :: iSide, iPorousSide
INTEGER                       :: nSendPorousSidesTmp(0:nLeaderGroupProcs-1)
INTEGER                       :: nRecvPorousSidesTmp(0:nLeaderGroupProcs-1)
INTEGER                       :: RecvRequest(0:nLeaderGroupProcs-1),SendRequest(0:nLeaderGroupProcs-1)
INTEGER                       :: SendPorousGlobalID(0:nLeaderGroupProcs-1,1:nComputeNodeSurfTotalSides)
INTEGER                       :: SampSizeAllocate
!===================================================================================================================================

nRecvPorousSidesTmp = 0

!--- Open receive buffer (number of sampling surfaces in other node's halo region)
DO iProc = 0,nLeaderGroupProcs-1
  IF (iProc.EQ.myLeaderGroupRank) CYCLE

  CALL MPI_IRECV( nRecvPorousSidesTmp(iProc)                                    &
                , 1                                                           &
                , MPI_INTEGER                                                 &
                , iProc                                                       &
                , 1213                                                        &
                , MPI_COMM_LEADERS_SHARED                                     &
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
DO iProc = 0,nLeaderGroupProcs-1
  IF (iProc.EQ.myLeaderGroupRank) CYCLE

  CALL MPI_ISEND( nSendPorousSidesTmp(iProc)                                    &
                , 1                                                           &
                , MPI_INTEGER                                                 &
                , iProc                                                       &
                , 1213                                                        &
                , MPI_COMM_LEADERS_SHARED                                     &
                , SendRequest(iProc)                                          &
                , IERROR)
END DO

!--- Finish communication
DO iProc = 0,nLeaderGroupProcs-1
  IF (iProc.EQ.myLeaderGroupRank) CYCLE

  CALL MPI_WAIT(SendRequest(iProc),MPISTATUS,IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  CALL MPI_WAIT(RecvRequest(iProc),MPISTATUS,IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO

IF (.NOT.SurfTotalSideOnNode) RETURN
SurfMapping(:)%nRecvPorousSides = 0
SurfMapping(:)%nSendPorousSides = 0

DO iProc = 0,nLeaderGroupProcs-1
  ! Ignore procs not on surface communicator
  IF (MPIRankSurfLeader(iProc).EQ.MPI_UNDEFINED) CYCLE
  ! Ignore myself
  IF (iProc .EQ. myLeaderGroupRank) CYCLE

  ! Save number of send and recv sides
  SurfMapping(MPIRankSurfLeader(iProc))%nRecvPorousSides = nRecvPorousSidesTmp(iProc)
  SurfMapping(MPIRankSurfLeader(iProc))%nSendPorousSides = nSendPorousSidesTmp(iProc)

  ! Only open recv buffer if we are expecting sides from this leader node
  IF (nRecvPorousSidesTmp(iProc).EQ.0) CYCLE
  ALLOCATE(SurfMapping(MPIRankSurfLeader(iProc))%RecvPorousGlobalID(1:nRecvPorousSidesTmp(iProc)))

  CALL MPI_IRECV( SurfMapping(MPIRankSurfLeader(iProc))%RecvPorousGlobalID    &
                , nRecvPorousSidesTmp(iProc)                                  &
                , MPI_INTEGER                                                 &
                , MPIRankSurfLeader(iProc)                                    &
                , 1213                                                        &
                , MPI_COMM_LEADERS_SURF                                       &
                , RecvRequest(MPIRankSurfLeader(iProc))                       &
                , IERROR)
END DO

DO iProc = 0,nLeaderGroupProcs-1
  ! Ignore procs not on surface communicator
  IF (MPIRankSurfLeader(iProc).EQ.MPI_UNDEFINED) CYCLE
  ! Ignore myself
  IF (iProc .EQ. myLeaderGroupRank) CYCLE

  ! Only open send buffer if we are expecting sides from this leader node
  IF (nSendPorousSidesTmp(iProc).EQ.0) CYCLE

  ALLOCATE(SurfMapping(MPIRankSurfLeader(iProc))%SendPorousGlobalID(1:nSendPorousSidesTmp(iProc)))

  SurfMapping(MPIRankSurfLeader(iProc))%SendPorousGlobalID = SendPorousGlobalID(iProc,1:nSendPorousSidesTmp(iProc))

  CALL MPI_ISEND( SurfMapping(MPIRankSurfLeader(iProc))%SendPorousGlobalID                         &
                , nSendPorousSidesTmp(iProc)                 &
                , MPI_INTEGER                                                 &
                , MPIRankSurfLeader(iProc)                                                       &
                , 1213                                                        &
                , MPI_COMM_LEADERS_SURF                                       &
                , SendRequest(MPIRankSurfLeader(iProc))                                          &
                , IERROR)
END DO

!--- Finish communication
DO iProc = 0,nLeaderGroupProcs-1
  ! Ignore procs not on surface communicator
  IF (MPIRankSurfLeader(iProc).EQ.MPI_UNDEFINED) CYCLE
  ! Ignore myself
  IF (iProc .EQ. myLeaderGroupRank) CYCLE

  IF (nSendPorousSidesTmp(iProc).NE.0) THEN
    CALL MPI_WAIT(SendRequest(MPIRankSurfLeader(iProc)),msg_status(:),IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END IF

  IF (nRecvPorousSidesTmp(iProc).NE.0) THEN
    CALL MPI_WAIT(RecvRequest(MPIRankSurfLeader(iProc)),msg_status(:),IERROR)
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


SUBROUTINE FinalizePorousBoundaryCondition()
!===================================================================================================================================
!> Deallocates
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Boundary_Vars
USE MOD_SurfaceModel_Vars
USE MOD_SurfaceModel_Analyze_Vars
#if USE_MPI
USE MOD_Particle_MPI_Vars       ,ONLY: PorousBCSendBuf,PorousBCRecvBuf
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared
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

IF(nPorousBC.GT.0) THEN
  SDEALLOCATE(PorousBC)
  SDEALLOCATE(PorousBCSampWall)
  SDEALLOCATE(PorousBCOutput)
  ! First, free every shared memory window. This requires MPI_BARRIER as per MPI3.1 specification
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
  CALL UNLOCK_AND_FREE(MapSurfSideToPorousSide_Shared_Win)
  CALL UNLOCK_AND_FREE(PorousBCInfo_Shared_Win)
  CALL UNLOCK_AND_FREE(PorousBCProperties_Shared_Win)
  CALL UNLOCK_AND_FREE(PorousBCSampWall_Shared_Win)
  CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
  SDEALLOCATE(PorousBCSendBuf)
  SDEALLOCATE(PorousBCRecvBuf)
#endif
  ADEALLOCATE(MapSurfSideToPorousSide_Shared)
  ADEALLOCATE(PorousBCInfo_Shared)
  ADEALLOCATE(PorousBCProperties_Shared)
  ADEALLOCATE(PorousBCSampWall_Shared)
END IF

END SUBROUTINE FinalizePorousBoundaryCondition

END MODULE MOD_SurfaceModel_Porous
