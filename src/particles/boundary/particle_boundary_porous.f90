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
!! Determines how particles interact with a given boundary condition 
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
CALL prms%CreateIntOption(      'Part-PorousBC-IterationOutput' &
                                , 'Every number of iterations an output (PorousBCInfo.csv) will be given. '//&
                                  'Choose larger numbers on HPC/parallel systems to avoid performance decrease', '1')
CALL prms%CreateIntOption(      'Part-PorousBC[$]-BC' &
                                , 'PartBound to be a porous boundary', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-PorousBC[$]-Pressure' &
                                , 'Pressure [Pa] at the porous boundary', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-PorousBC[$]-Temperature' &
                                , 'Temperature [K] at the porous boundary', numberedmulti=.TRUE.)
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
                                , 'Posibility to define a region, allowing small geometrical features such as holes w/o meshing. ' //&
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
! initialization of porous boundary condition
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Globals_Vars,                   ONLY: BoltzmannConst
USE MOD_ReadInTools
USE MOD_Mesh_Vars,                      ONLY: BC,nElems
USE MOD_Particle_Vars,                  ONLY: nSpecies
USE MOD_Particle_Boundary_Vars,         ONLY: nPartBound, PartBound, nPorousBC, PorousBC, MapBCtoPorousBC, SurfMesh, nPorousBCVars
USE MOD_Particle_Boundary_Vars,         ONLY: MapSurfSideToPorousSide, PorousBCOutputIter, PorousBCSampIter, PorousBCMacroVal
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iPorousBC, iSurfSide, BCSideID, SideNumber, PorousBCID
CHARACTER(32)         :: hilf, hilf2
REAL                  :: rmin, rmax
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)') ' INIT POROUS BOUNDARY CONDITION ...'
PorousBCOutputIter = GETINT('Part-PorousBC-IterationOutput', '1')
PorousBCSampIter = GETINT('Part-PorousBC-IterationMacroVal', '0')

IF(PorousBCSampIter.GT.0) THEN
  ALLOCATE(PorousBCMacroVal(1:7,1:nElems,1:nSpecies))
  PorousBCMacroVal = 0.0
  IF(PorousBCOutputIter.LT.PorousBCSampIter) THEN
    CALL abort(&
    __STAMP__&
    ,'ERROR: IterationForOutput should be larger than or equal to IterationForMacroVal!')
  END IF
END IF

ALLOCATE(PorousBC(1:nPorousBC))
ALLOCATE(MapBCtoPorousBC(1:nPartBound))
MapBCtoPorousBC = 0
ALLOCATE(MapSurfSideToPorousSide(1:SurfMesh%nTotalSides))
MapSurfSideToPorousSide = 0
DO iPorousBC = 1, nPorousBC
  WRITE(UNIT=hilf,FMT='(I0)') iPorousBC
  ! Mapping the porous BC to a boundary
  PorousBC(iPorousBC)%BC = GETINT('Part-PorousBC'//TRIM(hilf)//'-BC')
  IF(PartBound%TargetBoundCond(PorousBC(iPorousBC)%BC).NE.PartBound%ReflectiveBC) THEN
    CALL abort(__STAMP__&
      ,'ERROR in init of porous BC: given boundary condition must be reflective!')
  END IF
  !
  MapBCtoPorousBC(PorousBC(iPorousBC)%BC) = iPorousBC
  ! Read-in of the conditions at the porous boundary
  PorousBC(iPorousBC)%Pressure  = GETREAL('Part-PorousBC'//TRIM(hilf)//'-Pressure')
  PorousBC(iPorousBC)%Temperature  = GETREAL('Part-PorousBC'//TRIM(hilf)//'-Temperature')
  PorousBC(iPorousBC)%NumberDensity  = PorousBC(iPorousBC)%Pressure / (BoltzmannConst * PorousBC(iPorousBC)%Temperature)
  ! Initial pumping speed at the porous boundary
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
        PorousBC(iPorousBC)%origin = GETREALARRAY('Part-PorousBC'//TRIM(hilf)//'-origin',2,'0. , 0.')
        WRITE(UNIT=hilf2,FMT='(E16.8)') HUGE(PorousBC(iPorousBC)%rmax)
        PorousBC(iPorousBC)%rmax = GETREAL('Part-PorousBC'//TRIM(hilf)//'-rmax',TRIM(hilf2))
        PorousBC(iPorousBC)%rmin = GETREAL('Part-PorousBC'//TRIM(hilf)//'-rmin','0.')
      CASE DEFAULT
        CALL abort(__STAMP__&
        ,'ERROR in region definition of porous bc:', iPorousBC)
    END SELECT
  END IF    ! Region is given
END DO      ! nPorousBC

PorousBC(1:nPorousBC)%SideNumber = 0
! Counting the sides at the respective porous boundary
DO iSurfSide=1,SurfMesh%nTotalSides
  PorousBCID = MapBCtoPorousBC(PartBound%MapToPartBC(BC(SurfMesh%SurfIDToSideID(iSurfSide))))
  IF (PorousBCID.GT.0) PorousBC(PorousBCID)%SideNumber = PorousBC(PorousBCID)%SideNumber + 1
END DO

! Allocating the required arrays for each porous boundary condition
DO iPorousBC = 1, nPorousBC
  ALLOCATE(PorousBC(iPorousBC)%SideList(1:PorousBC(iPorousBC)%SideNumber))
  PorousBC(iPorousBC)%SideList = 0
  ALLOCATE(PorousBC(iPorousBC)%RemovalProbability(1:PorousBC(iPorousBC)%SideNumber))
  PorousBC(iPorousBC)%RemovalProbability = 0.
  ALLOCATE(PorousBC(iPorousBC)%PumpingSpeedSide(1:PorousBC(iPorousBC)%SideNumber))
  PorousBC(iPorousBC)%PumpingSpeedSide(1:PorousBC(iPorousBC)%SideNumber) = PorousBC(iPorousBC)%PumpingSpeed
  ALLOCATE(PorousBC(iPorousBC)%RegionSideType(1:PorousBC(iPorousBC)%SideNumber))
  PorousBC(iPorousBC)%RegionSideType = 0
  ALLOCATE(PorousBC(iPorousBC)%Sample(1:PorousBC(iPorousBC)%SideNumber,1:nPorousBCVars))
  PorousBC(iPorousBC)%Sample = 0
END DO

! Mapping of the porous BC sides to the sides
SideNumber = 1
DO iSurfSide=1,SurfMesh%nTotalSides
  BCSideID = SurfMesh%SurfIDToSideID(iSurfSide)
  PorousBCID = MapBCtoPorousBC(PartBound%MapToPartBC(BC(BCSideID)))
  IF (PorousBCID.GT.0) THEN
    PorousBC(PorousBCID)%SideList(SideNumber) = iSurfSide
    MapSurfSideToPorousSide(iSurfSide) = SideNumber
    ! Determine which cells are inside/outside/partially inside the defined region (0: inside, 1: outside, 2: partially)
    IF(PorousBC(PorousBCID)%UsingRegion) THEN
      SELECT CASE(PorousBC(PorousBCID)%Region)
        CASE('circular')
          CALL GetRadialDistance2D(BCSideID,PorousBC(PorousBCID)%dir,PorousBC(PorousBCID)%origin,rmin,rmax)
          IF ( (rmin .GT. PorousBC(PorousBCID)%rmax) .OR. (rmax .LT. PorousBC(PorousBCID)%rmin) ) THEN
            PorousBC(PorousBCID)%RegionSideType(SideNumber)=1
          ELSE IF ( (rmax .LE. PorousBC(PorousBCID)%rmax) .AND. (rmin .GE. PorousBC(PorousBCID)%rmin) ) THEN
            PorousBC(PorousBCID)%RegionSideType(SideNumber)=0
          ELSE
            PorousBC(PorousBCID)%RegionSideType(SideNumber)=2
          END IF
        CASE DEFAULT
          CALL abort(__STAMP__&
          ,'ERROR in region definition of porous bc:', PorousBCID)
      END SELECT
    END IF
    SideNumber = SideNumber + 1
  END IF
END DO

SWRITE(UNIT_stdOut,'(A)') ' INIT POROUS BOUNDARY CONDITION DONE!'

END SUBROUTINE InitPorousBoundaryCondition

SUBROUTINE PorousBoundaryTreatment(iPart,SideID,alpha,PartTrajectory,PorousReflection)
!===================================================================================================================================
! Treatment of particles impinging on the porous boundary
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars,          ONLY:PDM, LastPartPos, PartSpecies
USE MOD_Particle_Boundary_Vars, ONLY:PartBound, SurfMesh, MapBCtoPorousBC, PorousBC, MapSurfSideToPorousSide, HaloImpactCounter
USE MOD_Mesh_Vars,              ONLY:BC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPart, SideID
REAL, INTENT(IN)              :: PartTrajectory(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)            :: alpha
LOGICAL,INTENT(INOUT)         :: PorousReflection
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: point(1:2), intersectionPoint(1:3), radius, iRan
INTEGER                       :: iSpec, SurfSideID, PorousBCID, pBCSideID
!===================================================================================================================================
iSpec = PartSpecies(iPart)
SurfSideID = SurfMesh%SideIDToSurfID(SideID)
PorousBCID = MapBCtoPorousBC(PartBound%MapToPartBC(BC(SideID)))

IF(PorousBCID.GT.0) THEN
  pBCSideID = MapSurfSideToPorousSide(SurfSideID)
  IF(pBCSideID.EQ.0) THEN
    ! Particle will be reflected regurarly
    HaloImpactCounter = HaloImpactCounter + 1
    RETURN
  END IF
  IF(PorousBC(PorousBCID)%UsingRegion) THEN
    IF(PorousBC(PorousBCID)%RegionSideType(pBCSideID).EQ.0) THEN
      ! Side is completey inside the porous region
      PorousReflection = .TRUE.
    ELSEIF(PorousBC(PorousBCID)%RegionSideType(pBCSideID).EQ.2) THEN
      ! Side is partially inside the porous region (check if its within bounds)
      intersectionPoint(1:3) = LastPartPos(iPart,1:3) + alpha*PartTrajectory(1:3)
      point(1)=intersectionPoint(PorousBC(PorousBCID)%dir(2))-PorousBC(PorousBCID)%origin(1)
      point(2)=intersectionPoint(PorousBC(PorousBCID)%dir(3))-PorousBC(PorousBCID)%origin(2)
      radius=SQRT( (point(1))**2+(point(2))**2 )
      IF ((radius.LE.PorousBC(PorousBCID)%rmax).AND.(radius.GE.PorousBC(PorousBCID)%rmin)) THEN
        PorousReflection = .TRUE.
      END IF
    END IF
  ELSE
    PorousReflection = .TRUE.
  END IF
  ! Particle actually hit the porous boundary condition
  IF(PorousReflection) THEN
    PorousBC(PorousBCID)%Sample(pBCSideID,1)   = PorousBC(PorousBCID)%Sample(pBCSideID,1) + 1
    ! check particle for leaving the domain
    CALL RANDOM_NUMBER(iRan)
    IF(iRan.LE.PorousBC(PorousBCID)%RemovalProbability(pBCSideID)) THEN
      PDM%ParticleInside(iPart)=.FALSE.
      alpha=-1.
      PorousBC(PorousBCID)%Sample(pBCSideID,2) = PorousBC(PorousBCID)%Sample(pBCSideID,2) + 1
    END IF
  END IF
END IF

END SUBROUTINE PorousBoundaryTreatment

SUBROUTINE PorousBoundaryRemovalProb_Pressure()
!===================================================================================================================================
! Determing the removal probability and pumping speed for the porous boundary condition
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,           ONLY:BoltzmannConst
USE MOD_Particle_Vars,          ONLY:Species, nSpecies, Adaptive_MacroVal
USE MOD_Particle_Boundary_Vars, ONLY:SurfMesh, nPorousBC, PorousBC, PorousBCOutputIter, HaloImpactCounter
USE MOD_Mesh_Vars,              ONLY:SideToElem
USE MOD_Timedisc_Vars,          ONLY:iter, dt
#ifdef MPI
USE MOD_Particle_MPI_Vars,      ONLY: PartMPI
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPBCSideID, iPBC, ElemID, SumPartPorousBC
REAL                          :: PumpingSpeedTemp, DeltaPressure, PartWeight
!===================================================================================================================================

PartWeight = Species(1)%MacroParticleFactor

#ifdef MPI
CALL ExchangeImpingedPartPorousBC()
#endif

DO iPBC = 1,nPorousBC
  ! Summing up the number of impinged at the porous boundary
  SumPartPorousBC = SUM(PorousBC(iPBC)%Sample(1:PorousBC(iPBC)%SideNumber,1))
#ifdef MPI
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,SumPartPorousBC,1,MPI_INTEGER,MPI_SUM,PartMPI%COMM,iError)
#endif
  DO iPBCSideID = 1, PorousBC(iPBC)%SideNumber
    ! Only treat local sides
    IF(PorousBC(iPBC)%SideList(iPBCSideID).GT.SurfMesh%nSides) CYCLE
    IF(PorousBC(iPBC)%UsingRegion) THEN
      ! Skipping cell which are completely outside specified region
      IF(PorousBC(iPBC)%RegionSideType(iPBCSideID).EQ.1) CYCLE
    END IF
    ! Get the adjacent element to the BCSide (<- SurfSide (only reflective) <- PorousBCSide)
    ElemID = SideToElem(1,SurfMesh%SurfIDToSideID(PorousBC(iPBC)%SideList(iPBCSideID)))
    IF (ElemID.LT.1) THEN !not sure if necessary
      ElemID = SideToElem(2,SurfMesh%SurfIDToSideID(PorousBC(iPBC)%SideList(iPBCSideID)))
    END IF
    ! Determining the delta between current gas mixture pressure in adjacent cell and desired pressure
    DeltaPressure = SUM(Adaptive_MacroVal(12,ElemID,1:nSpecies))-PorousBC(iPBC)%Pressure
    ! Integrating the pressure difference
    IF(PorousBC(iPBC)%PumpingSpeedSide(iPBCSideID).GT.0.0) THEN
      Adaptive_MacroVal(13,ElemID,1) = Adaptive_MacroVal(13,ElemID,1) + DeltaPressure * dt
    ELSE
      Adaptive_MacroVal(13,ElemID,1) = 0.0
    END IF
    ! Adapting the pumping speed S (m^3/s) according to the pressure difference (control through proportional and integral part)
    PumpingSpeedTemp = PorousBC(iPBC)%PumpingSpeedSide(iPBCSideID) + PorousBC(iPBC)%DeltaPumpingSpeedKp * DeltaPressure &
        + PorousBC(iPBC)%DeltaPumpingSpeedKi * Adaptive_MacroVal(13,ElemID,1)
    ! Calculate the removal probability if any particles hit the pump
    IF(SumPartPorousBC.GT.0) THEN
      PorousBC(iPBC)%RemovalProbability(iPBCSideID) = PumpingSpeedTemp*SUM(Adaptive_MacroVal(7,ElemID,1:nSpecies)) * dt &
                                                                  / (REAL(SumPartPorousBC)*PartWeight)
    ELSE
      PorousBC(iPBC)%RemovalProbability(iPBCSideID) = 0.0
    END IF
    ! Making sure that the removal probability is between 0 and 1
    IF(PorousBC(iPBC)%RemovalProbability(iPBCSideID).GT.1.0) THEN
      PorousBC(iPBC)%RemovalProbability(iPBCSideID) = 1.0
      ! Setting pumping speed to maximum value (alpha=1)
      PorousBC(iPBC)%PumpingSpeedSide(iPBCSideID) = REAL(SumPartPorousBC)*PartWeight &
                                                    / (SUM(Adaptive_MacroVal(7,ElemID,1:nSpecies))*dt)
    ELSE IF(PorousBC(iPBC)%RemovalProbability(iPBCSideID).LE.0.0) THEN
      PorousBC(iPBC)%RemovalProbability(iPBCSideID) = 0.0
      ! Avoiding negative pumping speeds
      PorousBC(iPBC)%PumpingSpeedSide(iPBCSideID) = 0.0
    ELSE
      ! Only adapting the pumping speed if alpha is between zero and one
      PorousBC(iPBC)%PumpingSpeedSide(iPBCSideID) = PumpingSpeedTemp
    END IF
    ! Storing the pumping speed for the restart state file
    Adaptive_MacroVal(11,ElemID,1) = PorousBC(iPBC)%PumpingSpeedSide(iPBCSideID)
    ! -------- Sampling for output -------------------------------------------------------------------------------------------------
    ! Sampling the actual instantaneous pumping speed S (m^3/s) through the number of deleted particles  (-PumpSpeed-Measure-)
    PorousBC(iPBC)%Output(2) = PorousBC(iPBC)%Output(2) + REAL(PorousBC(iPBC)%Sample(iPBCSideID,2)) * PartWeight &
                                                                / (SUM(Adaptive_MacroVal(7,ElemID,1:nSpecies))*dt)
    IF(PorousBC(iPBC)%Sample(iPBCSideID,1).GT.0) THEN
      ! Counting only sides, where a particle hit the pump
      PorousBC(iPBC)%Output(1) = PorousBC(iPBC)%Output(1) + 1.
      ! Removal probability
      PorousBC(iPBC)%Output(3) = PorousBC(iPBC)%Output(3) + PorousBC(iPBC)%RemovalProbability(iPBCSideID)
      ! Normalized pressure at the pump (sampled over PumpSampIter number of iterations)
      PorousBC(iPBC)%Output(4) = PorousBC(iPBC)%Output(4) + SUM(Adaptive_MacroVal(12,ElemID,1:nSpecies)) / PorousBC(iPBC)%Pressure
      ! Pumping speed S (m^3/s) used at the pump to calculate the removal probability (-PumpSpeed-Control-)
      PorousBC(iPBC)%Output(5) = PorousBC(iPBC)%Output(5) + PorousBC(iPBC)%PumpingSpeedSide(iPBCSideID)
    END IF
    ! Reset of the sampled particle numbers at the pump
    PorousBC(iPBC)%Sample(iPBCSideID,1:2) = 0
  END DO
END DO

IF(MOD(iter+1,PorousBCOutputIter).EQ.0) THEN
#ifdef MPI
  ! Communicating the values and summation at the root
  DO iPBC = 1, nPorousBC
    IF(MPIRoot) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,PorousBC(iPBC)%Output,5,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,iError)
    ELSE
      CALL MPI_REDUCE(PorousBC(iPBC)%Output,PorousBC(iPBC)%Output,5,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,iError)
    END IF
  END DO
  IF(MPIRoot) THEN
#endif
    ! Averaging of the pump info across the whole pump area
    DO iPBC = 1,nPorousBC
      IF(PorousBC(iPBC)%Output(1).GT.0.0) THEN
        ! Pumping Speed (Output(2)) is the sum of all elements (counter over particles exiting through pump)
        PorousBC(iPBC)%Output(2) = PorousBC(iPBC)%Output(2) / REAL(PorousBCOutputIter)
        ! Other variales are averaged over the elements (sum of elements increasing each PorousBCOutputIter)
        PorousBC(iPBC)%Output(3:5) = PorousBC(iPBC)%Output(3:5) / PorousBC(iPBC)%Output(1)
      END IF
    END DO
    CALL WritePorousBCInfo()
#ifdef MPI
  END IF
#endif
  ! Reset of sampled values after the output
  DO iPBC = 1,nPorousBC
    PorousBC(iPBC)%Output(1:5) = 0.
  END DO
END IF

!! Exchange removal probability for halo cells
! #ifdef MPI
! CALL ExchangeRemovalProbabilityPorousBC()
! #endif

! DEBUG
IF(HaloImpactCounter.GT.0) THEN
  IPWRITE(*,*) 'WARNING: Particles impinging on the halo sides were not treated with the correct removal probability'
  IPWRITE(*,*) 'WARNING: A total of ', HaloImpactCounter, 'out of', SumPartPorousBC, ' particles were mistreated!'
  HaloImpactCounter = 0
END IF

END SUBROUTINE PorousBoundaryRemovalProb_Pressure


#ifdef MPI
SUBROUTINE ExchangeImpingedPartPorousBC()
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Boundary_Vars      ,ONLY:SurfComm, nPorousBC, nPorousBCVars, MapSurfSideToPorousSide, MapBCtoPorousBC
USE MOD_Particle_Boundary_Vars      ,ONLY:PartBound, SurfMesh, PorousBC
USE MOD_Particle_MPI_Vars           ,ONLY:PorousBCSendBuf,PorousBCRecvBuf,SurfExchange
USE MOD_Mesh_Vars                   ,ONLY:BC
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
  CALL MPI_IRECV( PorousBCRecvBuf(iProc)%content_int           &
                , MessageSize                                  &
                , MPI_INTEGER                                  &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID     &
                , 1009                                         &
                , SurfCOMM%COMM                                &
                , SurfExchange%RecvRequest(iProc)              &
                , IERROR )
END DO ! iProc

! build message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesSend(iProc).EQ.0) CYCLE
  iPos=0
  PorousBCSendBuf(iProc)%content_int = 0
  DO iSurfSide=1,SurfExchange%nSidesSend(iProc)
    SurfSideID = SurfCOMM%MPINeighbor(iProc)%SendList(iSurfSide)
    PorousBCSideID = MapSurfSideToPorousSide(SurfSideID)
    PorousBCID = MapBCtoPorousBC(PartBound%MapToPartBC(BC(SurfMesh%SurfIDToSideID(SurfSideID))))
    IF(PorousBCID.GT.0) THEN
      PorousBCSendBuf(iProc)%content_int(iPos+1:iPos+nVar) = PorousBC(PorousBCID)%Sample(PorousBCSideID,1:nVar)
      iPos=iPos+nVar
      PorousBC(PorousBCID)%Sample(PorousBCSideID,:)=0
    END IF
  END DO ! iSurfSide=1,nSurfExchange%nSidesSend(iProc)
END DO

! send message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesSend(iProc).EQ.0) CYCLE
  MessageSize=SurfExchange%nSidesSend(iProc)*nValues
  CALL MPI_ISEND( PorousBCSendBuf(iProc)%content_int       &
                , MessageSize                              &
                , MPI_INTEGER                              &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID &
                , 1009                                     &
                , SurfCOMM%COMM                            &
                , SurfExchange%SendRequest(iProc)          &
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
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%RecvList(iSurfSide)
    PorousBCSideID = MapSurfSideToPorousSide(SurfSideID)
    PorousBCID = MapBCtoPorousBC(PartBound%MapToPartBC(BC(SurfMesh%SurfIDToSideID(SurfSideID))))
    IF(PorousBCID.GT.0) THEN
      PorousBC(PorousBCID)%Sample(PorousBCSideID,1:nVar) = PorousBC(PorousBCID)%Sample(PorousBCSideID,1:nVar) &
                                                          + PorousBCRecvBuf(iProc)%content_int(iPos+1:iPos+nVar)
      iPos=iPos+nVar
    END IF
  END DO ! iSurfSide=1,nSurfExchange%nSidesSend(iProc)
  PorousBCRecvBuf(iProc)%content_int = 0.
END DO ! iProc

END SUBROUTINE ExchangeImpingedPartPorousBC
#endif /*MPI*/


SUBROUTINE WritePorousBCInfo()
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_TimeDisc_Vars          ,ONLY: time, dt
USE MOD_Globals                ,ONLY: FILEEXISTS
USE MOD_Restart_Vars           ,ONLY: DoRestart
USE MOD_Particle_Boundary_Vars ,ONLY: nPorousBC, PorousBC
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: OutputCounter, unit_index, iPBC
CHARACTER(LEN=350)          :: outfile
LOGICAL                     :: isOpen
!===================================================================================================================================
unit_index = 789
outfile = 'PorousBCInfo.csv'
OutputCounter = 2

INQUIRE(UNIT   = unit_index , OPENED = isOpen)
IF (.NOT.isOpen) THEN
  IF (DoRestart.and.FILEEXISTS(outfile)) THEN
    OPEN(unit_index,file=TRIM(outfile),action="WRITE",position="APPEND",status="OLD")
  ELSE
    OPEN(unit_index,file=TRIM(outfile),action="WRITE",status="REPLACE")
    WRITE(unit_index,'(A4)',ADVANCE='NO') 'Time'
    DO iPBC = 1, nPorousBC
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,'(I3.3,A19,I3.3)',ADVANCE='NO') OutputCounter,'-PumpSpeed-Measure-', iPBC
      OutputCounter = OutputCounter + 1
    END DO
    DO iPBC = 1, nPorousBC
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,'(I3.3,A13,I3.3)',ADVANCE='NO') OutputCounter,'-RemovalProb-', iPBC
      OutputCounter = OutputCounter + 1
    END DO
    DO iPBC = 1, nPorousBC
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,'(I3.3,A11,I3.3)',ADVANCE='NO') OutputCounter,'-PressNorm-', iPBC
      OutputCounter = OutputCounter + 1
    END DO
    DO iPBC = 1, nPorousBC
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,'(I3.3,A19,I3.3)',ADVANCE='NO') OutputCounter,'-PumpSpeed-Control-', iPBC
      OutputCounter = OutputCounter + 1
    END DO
    WRITE(unit_index,'(A1)') ' '
  END IF
END IF
WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') Time+dt
DO iPBC=1, nPorousBC
  WRITE(unit_index,'(A1)',ADVANCE='NO') ','
  WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PorousBC(iPBC)%Output(2)
END DO
DO iPBC=1, nPorousBC
  WRITE(unit_index,'(A1)',ADVANCE='NO') ','
  WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PorousBC(iPBC)%Output(3)
END DO
DO iPBC=1, nPorousBC
  WRITE(unit_index,'(A1)',ADVANCE='NO') ','
  WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PorousBC(iPBC)%Output(4)
END DO
DO iPBC=1, nPorousBC
  WRITE(unit_index,'(A1)',ADVANCE='NO') ','
  WRITE(unit_index,WRITEFORMAT,ADVANCE='NO') PorousBC(iPBC)%Output(5)
END DO
WRITE(unit_index,'(A1)') ' '

END SUBROUTINE WritePorousBCInfo


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

END MODULE MOD_Particle_Boundary_Porous