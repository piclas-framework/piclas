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

PUBLIC::DefineParametersPorousBC, InitPorousBoundaryCondition, PorousBoundaryTreatment, PorousBoundaryRemovalProb_Pressure
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
USE MOD_Mesh_Vars                   ,ONLY: BC,nElems, SideToElem
USE MOD_Particle_Vars               ,ONLY: nSpecies, Adaptive_MacroVal, Symmetry2D, Symmetry2DAxisymmetric
USE MOD_Particle_Boundary_Vars      ,ONLY: PartBound, nPorousBC, PorousBC, SurfMesh, nPorousBCVars, PorousBCMacroVal
USE MOD_Particle_Boundary_Vars      ,ONLY: MapSurfSideToPorousSide, PorousBCSampIter, MapSurfSideToPorousBC
USE MOD_Particle_Tracking_Vars      ,ONLY: DoRefMapping
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iPorousBC, iSurfSide, BCSideID, SideNumber(nPorousBC), PorousBCID, PorousBCSideID, iSurfSideTmp
CHARACTER(32)         :: hilf, hilf2
REAL                  :: rmin, rmax
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
ALLOCATE(MapSurfSideToPorousSide(1:SurfMesh%nTotalSides))
MapSurfSideToPorousSide = 0
ALLOCATE(MapSurfSideToPorousBC(1:SurfMesh%nTotalSides))
MapSurfSideToPorousBC = 0

DO iPorousBC = 1, nPorousBC
  WRITE(UNIT=hilf,FMT='(I0)') iPorousBC
  PorousBC(iPorousBC)%BC = GETINT('Part-PorousBC'//TRIM(hilf)//'-BC')
  IF(PartBound%TargetBoundCond(PorousBC(iPorousBC)%BC).NE.PartBound%ReflectiveBC) THEN
    CALL abort(__STAMP__&
      ,'ERROR in init of porous BC: given boundary condition must be reflective!')
  END IF
  ! Read-in of the conditions at the porous boundary
  PorousBC(iPorousBC)%Pressure  = GETREAL('Part-PorousBC'//TRIM(hilf)//'-Pressure')
  PorousBC(iPorousBC)%Temperature  = GETREAL('Part-PorousBC'//TRIM(hilf)//'-Temperature')
  PorousBC(iPorousBC)%Type  = TRIM(GETSTR('Part-PorousBC'//TRIM(hilf)//'-Type'))
  SELECT CASE(PorousBC(iPorousBC)%Type)
    CASE('sensor')
      PorousBC(iPorousBC)%PumpingSpeed = 0.0
      PorousBC(iPorousBC)%DeltaPumpingSpeedKp = 0.0
      PorousBC(iPorousBC)%DeltaPumpingSpeedKi = 0.0
    CASE('pump')
      ! Initial pumping speed at the porous boundary [m3/s]
      PorousBC(iPorousBC)%PumpingSpeed = GETREAL('Part-PorousBC'//TRIM(hilf)//'-PumpingSpeed')
      ! Proportional and integral factors for the control of the pumping speed
      PorousBC(iPorousBC)%DeltaPumpingSpeedKp = GETREAL('Part-PorousBC'//TRIM(hilf)//'-DeltaPumpingSpeed-Kp')
      PorousBC(iPorousBC)%DeltaPumpingSpeedKi = GETREAL('Part-PorousBC'//TRIM(hilf)//'-DeltaPumpingSpeed-Ki')
      ! Determining the order of magnitude
      IF(PorousBC(iPorousBC)%Pressure.GT.0.0) THEN
        PorousBC(iPorousBC)%DeltaPumpingSpeedKp = PorousBC(iPorousBC)%DeltaPumpingSpeedKp &
                                                        / 10.0**(ANINT(LOG10(PorousBC(iPorousBC)%Pressure)))
        PorousBC(iPorousBC)%DeltaPumpingSpeedKi = PorousBC(iPorousBC)%DeltaPumpingSpeedKi &
                                                        / 10.0**(ANINT(LOG10(PorousBC(iPorousBC)%Pressure)))
      END IF
    CASE DEFAULT
      CALL abort(__STAMP__&
      ,'ERROR in type definition of porous bc:', iPorousBC)
  END SELECT
  PorousBC(iPorousBC)%Region = TRIM(GETSTR('Part-PorousBC'//TRIM(hilf)//'-Region','none'))
  IF(PorousBC(iPorousBC)%Region.EQ.'none') THEN
    PorousBC(iPorousBC)%UsingRegion = .FALSE.
  ELSE
    PorousBC(iPorousBC)%UsingRegion = .TRUE.
    ! Normal direction to the surface (could be replaced by the normal vector of the boundary?)
    PorousBC(iPorousBC)%dir(1) = GETINT('Part-PorousBC'//TRIM(hilf)//'-normalDir')
    IF (PorousBC(iPorousBC)%dir(1).EQ.1) THEN
        PorousBC(iPorousBC)%dir(2)=2
        PorousBC(iPorousBC)%dir(3)=3
    ELSE IF (PorousBC(iPorousBC)%dir(1).EQ.2) THEN
        PorousBC(iPorousBC)%dir(2)=3
        PorousBC(iPorousBC)%dir(3)=1
    ELSE IF (PorousBC(iPorousBC)%dir(1).EQ.3) THEN
        PorousBC(iPorousBC)%dir(2)=1
        PorousBC(iPorousBC)%dir(3)=2
    ELSE
      CALL abort(__STAMP__&
        ,'ERROR in init: normalDir for PorousBC must be between 1 and 3!')
    END IF
    SELECT CASE(PorousBC(iPorousBC)%Region)
      CASE('circular')
        PorousBC(iPorousBC)%origin = GETREALARRAY('Part-PorousBC'//TRIM(hilf)//'-origin',2)
        WRITE(UNIT=hilf2,FMT='(E16.8)') HUGE(PorousBC(iPorousBC)%rmax)
        PorousBC(iPorousBC)%rmax = GETREAL('Part-PorousBC'//TRIM(hilf)//'-rmax',TRIM(hilf2))
        PorousBC(iPorousBC)%rmin = GETREAL('Part-PorousBC'//TRIM(hilf)//'-rmin','0.')
        IF(Symmetry2DAxisymmetric) THEN
          IF(PorousBC(iPorousBC)%dir(1).NE.1) THEN
            CALL abort(__STAMP__&
              ,'ERROR in Porous BC: For axisymmetric simulations, only regions perpendicular to the axis are allowed!', iPorousBC)
          END IF
          IF(PorousBC(iPorousBC)%origin(1)*PorousBC(iPorousBC)%origin(2).NE.0.0) THEN
            CALL abort(__STAMP__&
              ,'ERROR in Porous BC: For axisymmetric simulations, the origin has to be at (0,0)!', iPorousBC)
          END IF
        END IF
      CASE DEFAULT
        CALL abort(__STAMP__&
        ,'ERROR in region definition of porous bc:', iPorousBC)
    END SELECT
  END IF    ! Region is given
END DO      ! nPorousBC

! 2) Mapping of the porous BC sides to the respective surface side
SideNumber = 0
iSurfSideTmp = 0
DO iSurfSide=1,SurfMesh%nTotalSides
  BCSideID = SurfMesh%SurfIDToSideID(iSurfSide)
  ! Skip not reflective BC sides
  IF(PartBound%TargetBoundCond(PartBound%MapToPartBC(BC(BCSideID))).NE.PartBound%ReflectiveBC) CYCLE
  DO iPorousBC = 1, nPorousBC
    IF(PorousBC(iPorousBC)%BC.EQ.PartBound%MapToPartBC(BC(BCSideID))) THEN
      ! Determine which cells are inside/outside/partially inside the defined region (0: inside, 1: partially)
      IF(PorousBC(iPorousBC)%UsingRegion) THEN
        SELECT CASE(PorousBC(iPorousBC)%Region)
          CASE('circular')
            CALL GetRadialDistance2D(BCSideID,PorousBC(iPorousBC)%dir,PorousBC(iPorousBC)%origin,rmin,rmax)
            IF ( (rmin .GT. PorousBC(iPorousBC)%rmax) .OR. (rmax .LT. PorousBC(iPorousBC)%rmin) ) CYCLE
        END SELECT
      END IF
      IF(iSurfSide.EQ.iSurfSideTmp) THEN
        CALL abort(__STAMP__&
        ,'ERROR in Porous BC: Side is already defined by another porous BC. Make sure the BCs do not overlap!')
      END IF
      SideNumber(iPorousBC) = SideNumber(iPorousBC) + 1
      MapSurfSideToPorousSide(iSurfSide) = SideNumber(iPorousBC)
      MapSurfSideToPorousBC(iSurfSide) = iPorousBC
      iSurfSideTmp = iSurfSide
    END IF
  END DO
END DO

! 3) Allocating the PorousBC arrays per BC, containing only sides with a (partially) porous side
PorousBC(:)%SideNumber = SideNumber(:)
DO iPorousBC = 1, nPorousBC
  ALLOCATE(PorousBC(iPorousBC)%SideList(1:PorousBC(iPorousBC)%SideNumber))
  PorousBC(iPorousBC)%SideList = 0
  ALLOCATE(PorousBC(iPorousBC)%RemovalProbability(1:PorousBC(iPorousBC)%SideNumber))
  PorousBC(iPorousBC)%RemovalProbability = 0.
  ALLOCATE(PorousBC(iPorousBC)%PumpingSpeedSide(1:PorousBC(iPorousBC)%SideNumber))
  PorousBC(iPorousBC)%PumpingSpeedSide(1:PorousBC(iPorousBC)%SideNumber) = PorousBC(iPorousBC)%PumpingSpeed
  ALLOCATE(PorousBC(iPorousBC)%RegionSideType(1:PorousBC(iPorousBC)%SideNumber))
  PorousBC(iPorousBC)%RegionSideType = -1
  ALLOCATE(PorousBC(iPorousBC)%Sample(1:PorousBC(iPorousBC)%SideNumber,1:nPorousBCVars))
  PorousBC(iPorousBC)%Sample = 0
  ! Array for output variables in PartAnalyze.csv
  PorousBC(iPorousBC)%Output = 0.
END DO

! 4) Mapping the porous BC ID and the porous BC side ID to a surface side
DO iSurfSide=1,SurfMesh%nTotalSides
  PorousBCID = MapSurfSideToPorousBC(iSurfSide)
  IF(PorousBCID.EQ.0) CYCLE
  PorousBCSideID = MapSurfSideToPorousSide(iSurfSide)
  PorousBC(PorousBCID)%SideList(PorousBCSideID) = iSurfSide
  BCSideID = SurfMesh%SurfIDToSideID(iSurfSide)
  ! If a pumping speed was read-in during restart, use it (only for local elements and not halo sides)
  IF(iSurfSide.LE.SurfMesh%nSides) THEN
    IF(Adaptive_MacroVal(11,SideToElem(1,BCSideID),1).GT.0.0) THEN
        PorousBC(PorousBCID)%PumpingSpeedSide(PorousBCSideID) = Adaptive_MacroVal(11,SideToElem(1,BCSideID),1)
    END IF
  END IF
  ! Determine which cells are inside/outside/partially inside the defined region (0: inside, 1: partially)
  IF(PorousBC(PorousBCID)%UsingRegion) THEN
    SELECT CASE(PorousBC(PorousBCID)%Region)
      CASE('circular')
        CALL GetRadialDistance2D(BCSideID,PorousBC(PorousBCID)%dir,PorousBC(PorousBCID)%origin,rmin,rmax)
        IF ( (rmax .LE. PorousBC(PorousBCID)%rmax) .AND. (rmin .GE. PorousBC(PorousBCID)%rmin) ) THEN
          PorousBC(PorousBCID)%RegionSideType(PorousBCSideID)=0
        ELSE
          PorousBC(PorousBCID)%RegionSideType(PorousBCSideID)=1
        END IF
    END SELECT
    IF(PorousBC(PorousBCID)%RegionSideType(PorousBCSideID).LT.0) THEN
        CALL abort(__STAMP__&
        ,'ERROR in Porous BC: Region side type not defined! Porous BC ID:', PorousBCID)
    END IF
  END IF
END DO

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
USE MOD_Particle_Vars          ,ONLY: LastPartPos
USE MOD_Particle_Boundary_Vars ,ONLY: SurfMesh, MapSurfSideToPorousBC, PorousBC, MapSurfSideToPorousSide
USE MOD_part_tools             ,ONLY: GetParticleWeight
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_part_operations        ,ONLY: RemoveParticle
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
INTEGER                       :: SurfSideID, PorousBCID, pBCSideID
LOGICAL                       :: ParticleHitPorousBC
!===================================================================================================================================
SurfSideID = SurfMesh%SideIDToSurfID(SideID)
PorousBCID = MapSurfSideToPorousBC(SurfSideID)

IF(PorousBCID.GT.0) THEN
  pBCSideID = MapSurfSideToPorousSide(SurfSideID)
  ParticleHitPorousBC = .TRUE.
  ! 1) Determination whether the particle hit the porous BC region or only the regular BC
  IF(PorousBC(PorousBCID)%UsingRegion) THEN
    IF(PorousBC(PorousBCID)%RegionSideType(pBCSideID).EQ.1) THEN
      ! Side is partially inside the porous region (check if its within bounds)
      intersectionPoint(1:3) = LastPartPos(1:3,iPart) + alpha*PartTrajectory(1:3)
      point(1)=intersectionPoint(PorousBC(PorousBCID)%dir(2))-PorousBC(PorousBCID)%origin(1)
      point(2)=intersectionPoint(PorousBC(PorousBCID)%dir(3))-PorousBC(PorousBCID)%origin(2)
      radius=SQRT( (point(1))**2+(point(2))**2 )
      ! Check if particle hit outside the region
      IF ((radius.GT.PorousBC(PorousBCID)%rmax).OR.(radius.LT.PorousBC(PorousBCID)%rmin)) THEN
        ParticleHitPorousBC = .FALSE.
      END IF
    END IF
  END IF
  ! 2) Comparison of the removal probability with a random number to determine whether the particle is deleted
  IF(ParticleHitPorousBC) THEN
    ! Counting particles that are impinging the porous BC (required for the calculation of the removal probability for the next dt)
    PorousBC(PorousBCID)%Sample(pBCSideID,1)   = PorousBC(PorousBCID)%Sample(pBCSideID,1) + GetParticleWeight(iPart)
    SELECT CASE(PorousBC(PorousBCID)%Type)
      CASE('sensor')
        ! Treat the sensor area as a regular boundary condition
        ElasticReflectionAtPorousBC = .FALSE.
      CASE('pump')
        CALL RANDOM_NUMBER(iRan)
        IF(iRan.LE.PorousBC(PorousBCID)%RemovalProbability(pBCSideID)) THEN
          ! Counting particles that leave the domain through the porous BC (required for the calculation of the pumping capacity)
          PorousBC(PorousBCID)%Sample(pBCSideID,2) = PorousBC(PorousBCID)%Sample(pBCSideID,2) + GetParticleWeight(iPart)
          CALL RemoveParticle(iPart,BCID=PartBound%MapToPartBC(BC(SideID)),alpha=alpha)
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
USE MOD_Particle_Vars,          ONLY:Species, nSpecies, Adaptive_MacroVal, usevMPF, VarTimeStep
USE MOD_Particle_Boundary_Vars, ONLY:SurfMesh, nPorousBC, PorousBC, SampWall
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO
USE MOD_Mesh_Vars,              ONLY:SideToElem
USE MOD_Timedisc_Vars,          ONLY:dt
USE MOD_Particle_Analyze_Vars,  ONLY:CalcPorousBCInfo
USE MOD_DSMC_Vars,              ONLY:DSMC, RadialWeighting
USE MOD_Particle_VarTimeStep    ,ONLY: CalcVarTimeStep
#if USE_MPI
USE MOD_Particle_MPI_Vars,      ONLY:PartMPI
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPBCSideID, iPBC, ElemID, SurfSideID
REAL                          :: PumpingSpeedTemp, DeltaPressure, partWeight, SumPartPorousBC, dtVar
!===================================================================================================================================

IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
  partWeight = 1.
ELSE
  partWeight = Species(1)%MacroParticleFactor
END IF

! 1) MPI communication of particles that impinged on halo sides to corresponding side
#if USE_MPI
CALL ExchangeImpingedPartPorousBC()
#endif

! 2) Loop over all porous BCs
DO iPBC = 1,nPorousBC
  ! a) Summing up the number of impinged particles for the whole BC surface
  SumPartPorousBC = SUM(PorousBC(iPBC)%Sample(1:PorousBC(iPBC)%SideNumber,1))
#if USE_MPI
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,SumPartPorousBC,1,MPI_DOUBLE_PRECISION,MPI_SUM,PartMPI%COMM,iError)
#endif
  ! Zero the output variable
  IF(CalcPorousBCInfo) PorousBC(iPBC)%Output = 0.
  ! 2.1) Loop over all sides within each porous BC
  DO iPBCSideID = 1, PorousBC(iPBC)%SideNumber
    SurfSideID = PorousBC(iPBC)%SideList(iPBCSideID)
    ! Only treat local sides
    IF(SurfSideID.GT.SurfMesh%nSides) CYCLE
    ! Get the adjacent element to the BCSide (<- SurfSide (only reflective) <- PorousBCSide)
    ElemID = SideToElem(1,SurfMesh%SurfIDToSideID(SurfSideID))
    IF (ElemID.LT.1) THEN !not sure if necessary
      ElemID = SideToElem(2,SurfMesh%SurfIDToSideID(SurfSideID))
    END IF
    ! Skip element if number density is zero
    IF(SUM(Adaptive_MacroVal(DSMC_NUMDENS,ElemID,1:nSpecies)).EQ.0.0) CYCLE
    ! Get the correct time step of the cell
    IF(VarTimeStep%UseVariableTimeStep) THEN
      dtVar = dt * CalcVarTimeStep(GEO%ElemMidPoint(1,ElemID), GEO%ElemMidPoint(2,ElemID), ElemID)
    ELSE
      dtVar = dt
    END IF
    ! Determine the removal probability based on the pumping speed (adaptive to a target pressure or fixed)
    IF((PorousBC(iPBC)%DeltaPumpingSpeedKp.GT.0.).OR.(PorousBC(iPBC)%DeltaPumpingSpeedKi.GT.0.)) THEN
      ! a) Determining the delta between current gas mixture pressure in adjacent cell and target pressure
      DeltaPressure = SUM(Adaptive_MacroVal(12,ElemID,1:nSpecies))-PorousBC(iPBC)%Pressure
      ! Integrating the pressure difference (only utilized later if DeltaPumpingSpeedKi was given)
      IF(PorousBC(iPBC)%PumpingSpeedSide(iPBCSideID).GT.0.0) THEN
        Adaptive_MacroVal(13,ElemID,1) = Adaptive_MacroVal(13,ElemID,1) + DeltaPressure * dtVar
      ELSE
        Adaptive_MacroVal(13,ElemID,1) = 0.0
      END IF
      ! b) Adapting the pumping capacity (m^3/s) according to pressure difference (control through proportional and integral part)
      PumpingSpeedTemp = PorousBC(iPBC)%PumpingSpeedSide(iPBCSideID) + PorousBC(iPBC)%DeltaPumpingSpeedKp * DeltaPressure &
          + PorousBC(iPBC)%DeltaPumpingSpeedKi * Adaptive_MacroVal(13,ElemID,1)
      ! c) Calculate the removal probability if any particles hit the pump
      IF(SumPartPorousBC.GT.0) THEN
        PorousBC(iPBC)%RemovalProbability(iPBCSideID) = PumpingSpeedTemp*SUM(Adaptive_MacroVal(DSMC_NUMDENS,ElemID,1:nSpecies)) &
                                                        * dtVar / (SumPartPorousBC*partWeight)
      ELSE
        PorousBC(iPBC)%RemovalProbability(iPBCSideID) = 0.0
      END IF
      ! d) Limit removal probability to values between 0 and 1
      IF(PorousBC(iPBC)%RemovalProbability(iPBCSideID).GT.1.0) THEN
        PorousBC(iPBC)%RemovalProbability(iPBCSideID) = 1.0
        ! Setting pumping speed to maximum value (alpha=1)
        PorousBC(iPBC)%PumpingSpeedSide(iPBCSideID) = SumPartPorousBC*partWeight &
                                                      / (SUM(Adaptive_MacroVal(DSMC_NUMDENS,ElemID,1:nSpecies))*dtVar)
      ELSE IF(PorousBC(iPBC)%RemovalProbability(iPBCSideID).LE.0.0) THEN
        PorousBC(iPBC)%RemovalProbability(iPBCSideID) = 0.0
        ! Avoiding negative pumping speeds
        PorousBC(iPBC)%PumpingSpeedSide(iPBCSideID) = 0.0
      ELSE
        ! Only adapting the pumping speed if alpha is between zero and one
        PorousBC(iPBC)%PumpingSpeedSide(iPBCSideID) = PumpingSpeedTemp
      END IF
    ELSE IF(PorousBC(iPBC)%PumpingSpeedSide(iPBCSideID).GT.0.0) THEN
      ! Constant given pumping speed
      IF(SumPartPorousBC.GT.0) THEN
        PorousBC(iPBC)%RemovalProbability(iPBCSideID) = PorousBC(iPBC)%PumpingSpeedSide(iPBCSideID) &
                        * SUM(Adaptive_MacroVal(DSMC_NUMDENS,ElemID,1:nSpecies)) * dtVar / (SumPartPorousBC*partWeight)
      ELSE
        PorousBC(iPBC)%RemovalProbability(iPBCSideID) = 0.0
      END IF
    END IF
    ! Storing the pumping speed for the restart state file
    Adaptive_MacroVal(11,ElemID,1) = PorousBC(iPBC)%PumpingSpeedSide(iPBCSideID)
    ! e) Sampling of the pumping capacity (and other variables for PartAnalyze) for the output
    ! -------- Sampling for output in DSMCSurfState --------------------------------------------------------------------------------
    IF(DSMC%CalcSurfaceVal) THEN
      SampWall(SurfSideID)%PumpCapacity = SampWall(SurfSideID)%PumpCapacity + PorousBC(iPBC)%Sample(iPBCSideID,2) &
                                                  * partWeight / (SUM(Adaptive_MacroVal(DSMC_NUMDENS,ElemID,1:nSpecies))*dtVar)
    END IF
    ! -------- Sampling for output in PartAnalyze ----------------------------------------------------------------------------------
    IF(CalcPorousBCInfo) THEN
      ! Sampling the actual instantaneous pumping speed S (m^3/s) through the number of deleted particles  (-PumpSpeed-Measure-)
      PorousBC(iPBC)%Output(2) = PorousBC(iPBC)%Output(2) + PorousBC(iPBC)%Sample(iPBCSideID,2) * PartWeight &
                                                            / (SUM(Adaptive_MacroVal(DSMC_NUMDENS,ElemID,1:nSpecies))*dtVar)
      IF(PorousBC(iPBC)%Sample(iPBCSideID,1).GT.0) THEN
        ! Counting only sides, where a particle hit the pump
        PorousBC(iPBC)%Output(1) = PorousBC(iPBC)%Output(1) + 1.
        ! Pumping speed S (m^3/s) used at the pump to calculate the removal probability (-PumpSpeed-Control-)
        PorousBC(iPBC)%Output(3) = PorousBC(iPBC)%Output(3) + PorousBC(iPBC)%PumpingSpeedSide(iPBCSideID)
        ! Removal probability
        PorousBC(iPBC)%Output(4) = PorousBC(iPBC)%Output(4) + PorousBC(iPBC)%RemovalProbability(iPBCSideID)
        ! Normalized pressure at the pump
        PorousBC(iPBC)%Output(5) = PorousBC(iPBC)%Output(5) + SUM(Adaptive_MacroVal(12,ElemID,1:nSpecies)) / PorousBC(iPBC)%Pressure
      END IF
    END IF
    ! Reset of the sampled particle numbers at the pump
    PorousBC(iPBC)%Sample(iPBCSideID,1:2) = 0
  END DO  ! iPBCSideID=1, SideNumber
END DO    ! iPBC=1, nPorousBC

! Exchange removal probability for halo cells
#if USE_MPI
CALL ExchangeRemovalProbabilityPorousBC()
#endif

END SUBROUTINE PorousBoundaryRemovalProb_Pressure


#if USE_MPI
SUBROUTINE ExchangeImpingedPartPorousBC()
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfComm, nPorousBC, nPorousBCVars, MapSurfSideToPorousSide, MapSurfSideToPorousBC
USE MOD_Particle_Boundary_Vars  ,ONLY: PorousBC, SurfMesh
USE MOD_Particle_MPI_Vars       ,ONLY: PorousBCSendBuf, PorousBCRecvBuf, SurfExchange
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iPos,iProc, MessageSize,nValues,iSurfSide,SurfSideID,nVar, PorousBCID, PorousBCSideID
INTEGER                         :: recv_status_list(1:MPI_STATUS_SIZE,1:SurfCOMM%nMPINeighbors)
!===================================================================================================================================

nVar = nPorousBCVars
nValues = nVar*nPorousBC
!
! open receive buffer
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesRecv(iProc).EQ.0) CYCLE
  MessageSize=SurfExchange%nSidesRecv(iProc)*nValues
  CALL MPI_IRECV( PorousBCRecvBuf(iProc)%content                &
                , MessageSize                                   &
                , MPI_DOUBLE_PRECISION                          &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID      &
                , 1009                                          &
                , SurfCOMM%COMM                                 &
                , SurfExchange%RecvRequest(iProc)               &
                , IERROR )
END DO ! iProc

! build message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesSend(iProc).EQ.0) CYCLE
  iPos=0
  PorousBCSendBuf(iProc)%content = 0.
  DO iSurfSide=1,SurfExchange%nSidesSend(iProc)
    SurfSideID = SurfMesh%SideIDToSurfID(SurfCOMM%MPINeighbor(iProc)%SendList(iSurfSide))
    PorousBCID = MapSurfSideToPorousBC(SurfSideID)
    IF(PorousBCID.GT.0) THEN
      PorousBCSideID = MapSurfSideToPorousSide(SurfSideID)
      PorousBCSendBuf(iProc)%content(iPos+1:iPos+nVar) = PorousBC(PorousBCID)%Sample(PorousBCSideID,1:nVar)
      iPos=iPos+nVar
      PorousBC(PorousBCID)%Sample(PorousBCSideID,:)=0.
    END IF
  END DO ! iSurfSide=1,nSurfExchange%nSidesSend(iProc)
END DO

! send message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesSend(iProc).EQ.0) CYCLE
  MessageSize=SurfExchange%nSidesSend(iProc)*nValues
  CALL MPI_ISEND( PorousBCSendBuf(iProc)%content            &
                , MessageSize                               &
                , MPI_DOUBLE_PRECISION                      &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID  &
                , 1009                                      &
                , SurfCOMM%COMM                             &
                , SurfExchange%SendRequest(iProc)           &
                , IERROR )
END DO ! iProc

! 4) Finish Received number of particles
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesSend(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfExchange%SendRequest(iProc),MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error', IERROR)
  END IF
  IF(SurfExchange%nSidesRecv(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfExchange%RecvRequest(iProc),recv_status_list(:,iProc),IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error', IERROR)
  END IF
END DO ! iProc

! add data do my list
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesRecv(iProc).EQ.0) CYCLE
  iPos=0
  DO iSurfSide=1,SurfExchange%nSidesRecv(iProc)
    SurfSideID = SurfMesh%SideIDToSurfID(SurfCOMM%MPINeighbor(iProc)%RecvList(iSurfSide))
    PorousBCID = MapSurfSideToPorousBC(SurfSideID)
    IF(PorousBCID.GT.0) THEN
      PorousBCSideID = MapSurfSideToPorousSide(SurfSideID)
      PorousBC(PorousBCID)%Sample(PorousBCSideID,1:nVar) = PorousBC(PorousBCID)%Sample(PorousBCSideID,1:nVar) &
                                                          + PorousBCRecvBuf(iProc)%content(iPos+1:iPos+nVar)
      iPos=iPos+nVar
    END IF
  END DO ! iSurfSide=1,nSurfExchange%nSidesSend(iProc)
  PorousBCRecvBuf(iProc)%content = 0.
END DO ! iProc

END SUBROUTINE ExchangeImpingedPartPorousBC


SUBROUTINE ExchangeRemovalProbabilityPorousBC
!===================================================================================================================================
! Communication of the removal probability to the halo cell of other procs (since the removal probability requires the locally
! sampled values of the adjacent elems) for the case, when a particle leaves the local domain and hits the porous BC on a haloside
! 1) Building the send message using the receive list (since we usually receive from these halo sides)
! 2) Communication
! 3) Mapping of the received info on our halo elements using the send list (since we usually send away info on our halo sides)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Boundary_Vars      ,ONLY: SurfComm, MapSurfSideToPorousSide, MapSurfSideToPorousBC, PorousBC, SurfMesh
USE MOD_Particle_MPI_Vars           ,ONLY: SurfExchange
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                                 :: iCommSide, iProc, SurfSideID, pBCID
  TYPE tTempArrayProc
    REAL, ALLOCATABLE                     :: SendMsg(:)
    REAL, ALLOCATABLE                     :: RecvMsg(:)
  END TYPE
  TYPE(tTempArrayProc), ALLOCATABLE       :: TempArrayProc(:)
!===================================================================================================================================

! 1) Communication of the removal probability to the halo cell of other procs (since the removal probability requires the locally
! sampled values of the adjacent elems) for the case, when a particle leaves the local domain and hits the porous BC on a haloside
ALLOCATE(TempArrayProc(1:SurfCOMM%nMPINeighbors))
DO iProc=1, SurfCOMM%nMPINeighbors
  ALLOCATE(TempArrayProc(iProc)%SendMsg(1:SurfExchange%nSidesRecv(iProc)))
  ALLOCATE(TempArrayProc(iProc)%RecvMsg(1:SurfExchange%nSidesSend(iProc)))
  TempArrayProc(iProc)%SendMsg(1:SurfExchange%nSidesRecv(iProc)) = 0.
  TempArrayProc(iProc)%RecvMsg(1:SurfExchange%nSidesSend(iProc)) = 0.
END DO

! 2) Building the send message: using the receive list to get the surface sides for which we need to send away the information
! These are the halo sides (from which we usually receive information -> using nSidesRecv and RecvList)
DO iProc=1, SurfCOMM%nMPINeighbors
  DO iCommSide=1, SurfExchange%nSidesRecv(iProc)
    SurfSideID = SurfMesh%SideIDToSurfID(SurfCOMM%MPINeighbor(iProc)%RecvList(iCommSide))
    pBCID = MapSurfSideToPorousBC(SurfSideID)
    IF(pBCID.GT.0) THEN
      TempArrayProc(iProc)%SendMsg(iCommSide) = PorousBC(pBCID)%RemovalProbability(MapSurfSideToPorousSide(SurfSideID))
    END IF
  END DO
END DO

! 3) Communication
DO iProc=1, SurfCOMM%nMPINeighbors
  IF (SurfCOMM%MyRank.LT.SurfCOMM%MPINeighbor(iProc)%NativeProcID) THEN
    IF (SurfExchange%nSidesRecv(iProc).NE.0) &
      CALL MPI_SEND(TempArrayProc(iProc)%SendMsg,SurfExchange%nSidesRecv(iProc), &
                      MPI_DOUBLE_PRECISION,SurfCOMM%MPINeighbor(iProc)%NativeProcID,1101,SurfCOMM%COMM,IERROR)
    IF (SurfExchange%nSidesSend(iProc).NE.0) &
      CALL MPI_RECV(TempArrayProc(iProc)%RecvMsg,SurfExchange%nSidesSend(iProc), &
                      MPI_DOUBLE_PRECISION,SurfCOMM%MPINeighbor(iProc)%NativeProcID,1101,SurfCOMM%COMM,MPISTATUS,IERROR)
  ELSE IF (SurfCOMM%MyRank.GT.SurfCOMM%MPINeighbor(iProc)%NativeProcID) THEN
    IF (SurfExchange%nSidesSend(iProc).NE.0) &
      CALL MPI_RECV(TempArrayProc(iProc)%RecvMsg,SurfExchange%nSidesSend(iProc), &
                      MPI_DOUBLE_PRECISION,SurfCOMM%MPINeighbor(iProc)%NativeProcID,1101,SurfCOMM%COMM,MPISTATUS,IERROR)
    IF (SurfExchange%nSidesRecv(iProc).NE.0) &
      CALL MPI_SEND(TempArrayProc(iProc)%SendMsg,SurfExchange%nSidesRecv(iProc), &
                      MPI_DOUBLE_PRECISION,SurfCOMM%MPINeighbor(iProc)%NativeProcID,1101,SurfCOMM%COMM,IERROR)
  END IF
END DO

! 4) Mapping of the received info on our halo elements (which we usually send away -> using nSidesSend and SendList)
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesSend(iProc).EQ.0) CYCLE
  DO iCommSide=1,SurfExchange%nSidesSend(iProc)
    SurfSideID = SurfMesh%SideIDToSurfID(SurfCOMM%MPINeighbor(iProc)%SendList(iCommSide))
    pBCID = MapSurfSideToPorousBC(SurfSideID)
    IF(pBCID.GT.0) THEN
      PorousBC(pBCID)%RemovalProbability(MapSurfSideToPorousSide(SurfSideID)) = TempArrayProc(iProc)%RecvMsg(iCommSide)
    END IF
  END DO ! iCommSide=1,nSurfExchange%nSidesSend(iProc)
END DO

END SUBROUTINE ExchangeRemovalProbabilityPorousBC
#endif /*USE_MPI*/


SUBROUTINE GetRadialDistance2D(BCSideID,dir,origin,rmin,rmax)
!===================================================================================================================================
! Determines the radial distance to a given origin on a surface
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Surfaces,      ONLY: GetSideBoundingBox
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: BCSideID, dir(3)
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
CALL GetSideBoundingBox(BCSideID,BoundingBox)
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
USE MOD_Particle_Boundary_Vars      ,ONLY: PorousBC, MapSurfSideToPorousBC, MapSurfSideToPorousSide, PorousBCMacroVal
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
SDEALLOCATE(MapSurfSideToPorousBC)
SDEALLOCATE(MapSurfSideToPorousSide)
END SUBROUTINE FinalizePorousBoundaryCondition

END MODULE MOD_Particle_Boundary_Porous
