!==================================================================================================================================
! Copyright (c) 2010 - 2019 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_part_emission
!===================================================================================================================================
! module for particle emission
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitializeParticleEmission
  MODULE PROCEDURE InitializeParticleEmission
END INTERFACE

INTERFACE ParticleInserting
  MODULE PROCEDURE ParticleInserting
END INTERFACE

INTERFACE AdaptiveBCAnalyze
  MODULE PROCEDURE AdaptiveBCAnalyze
END INTERFACE

!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC::InitializeParticleEmission,ParticleInserting,AdaptiveBCAnalyze
!===================================================================================================================================
PUBLIC::DefineParametersParticleEmission
CONTAINS

!==================================================================================================================================
!> Define parameters for particle emission (surface flux)
!==================================================================================================================================
SUBROUTINE DefineParametersParticleEmission()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Particle Emission")

CALL prms%CreateIntOption(      'Part-Species[$]-nSurfacefluxBCs'&
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Number of SF emissions', '0', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-Surfaceflux[$]-BC' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'PartBound to be emitted from', '0', numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Part-Species[$]-Surfaceflux[$]-velocityDistribution' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Specifying keyword for velocity distribution' , 'constant'&
, numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-VeloIC' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Velocity for inital Data', '0.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Species[$]-Surfaceflux[$]-VeloIsNormal' &
                                , 'TODO-DEFINE-PARAMETER VeloIC is in Surf-Normal instead of VeloVecIC' &
                                , '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-Surfaceflux[$]-VeloVecIC' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Normalized velocity vector' , '0.0 , 0.0 , 0.0', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Species[$]-Surfaceflux[$]-CircularInflow' &
                                , 'Enables the utilization of a circular region as a surface flux on the selected boundary. '//&
                                  'Only possible on surfaces, which are in xy, xz, and yz-planes.' &
                                , '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-Surfaceflux[$]-axialDir' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Axial direction of coordinates in polar system', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Part-Species[$]-Surfaceflux[$]-origin' &
                                , 'TODO-DEFINE-PARAMETER Origin in orth(ogonal?) coordinates of polar system' , '0.0 , 0.0'&
                                ,  numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-rmax' &
                                , 'TODO-DEFINE-PARAMETER Max radius of to-be inserted particles', '1e21', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-rmin' &
                                , 'TODO-DEFINE-PARAMETER Min radius of to-be inserted particles', '0.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Species[$]-Surfaceflux[$]-SimpleRadialVeloFit' &
                                      , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Fit of veloR/veloTot=-r*(A*exp(B*r)+C)', '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-preFac' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'A , see SimpleRadialVeloFit' &
                                , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-powerFac' &
                                      , 'TODO-DEFINE-PARAMETER\n'//&
                                  'B , see SimpleRadialVeloFit' &
                                , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-shiftFac' &
                                      , 'TODO-DEFINE-PARAMETER\n'//&
                                  'C , see SimpleRadialVeloFit' &
                                , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-MWTemperatureIC' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Temperature for Maxwell Distribution', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-PartDensity' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'PartDensity (real particles per m^3) or  (vpi_)cub./cyl. as alternative  to'//&
                                  ' Part.Emis. in Type1'  , '0.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Species[$]-Surfaceflux[$]-ReduceNoise' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Reduce stat. noise by global calc. of PartIns', '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Part-Species[$]-Surfaceflux[$]-AcceptReject' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  ' Perform ARM for skewness of RefMap-positioning', '.TRUE.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-Surfaceflux[$]-ARM_DmaxSampleN' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Number of sample intervals in xi/eta for Dmax-calc.', '1', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'DoForceFreeSurfaceFlux' &
                                , 'TODO-DEFINE-PARAMETER\n'//&
                                  'Flag if the stage reconstruction uses a force' , '.FALSE.')

CALL prms%CreateLogicalOption(  'OutputSurfaceFluxLinked' &
                                , 'Flag to print the SurfaceFlux-linked Info' , '.FALSE.')
! Parameters for adaptive boundary conditions
CALL prms%CreateLogicalOption(  'Part-Species[$]-Surfaceflux[$]-Adaptive' &
                                      , 'Flag for the definition of adaptive boundary conditions', '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-Surfaceflux[$]-Adaptive-Type' &
                                , 'Define the type of the adaptive boundary condition. Options:\n' //&
                                  '(1) Const. static pressure inlet after Farbar & Boyd 2014 (Type 1)\n' //&
                                  '(2) Const. static pressure outlet after Farbar & Boyd 2014 (Type 1)\n' //&
                                  '(3) Const. massflow inlet after Farbar & Boyd 2014 (Type 2): Number of particles for regular '//&
                                  'surface flux is calculated with velocity and given mass flow. Requires an open BC.' //&
                                  '(4) Const. massflow inlet after Lei 2017 (cf_3): N_in = N_mdot + N_out (counting particles, '//&
                                  'which exist the domain through the adaptive BC).' &
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-Adaptive-Pressure' &
                                , 'Static pressure in [Pa] for the adaptive boundary conditions of type 1 and 2.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surfaceflux[$]-Adaptive-Massflow' &
                                , 'Massflow in [kg/s] for the adaptive boundary conditions of type 3 and 4.', numberedmulti=.TRUE.)

END SUBROUTINE DefineParametersParticleEmission

SUBROUTINE InitializeParticleEmission()
!===================================================================================================================================
! Initialize particles / Insert initial particles
!===================================================================================================================================
! MODULES
#if USE_MPI
USE MOD_Particle_MPI_Vars   ,ONLY: PartMPI
#endif /*USE_MPI*/
USE MOD_Globals
USE MOD_Restart_Vars        ,ONLY: DoRestart
USE MOD_Particle_Vars       ,ONLY: Species,nSpecies,PDM,PEM, usevMPF, SpecReset, Symmetry
USE MOD_Particle_Mesh_Vars  ,ONLY: GEO
USE MOD_part_tools          ,ONLY: UpdateNextFreePosition
USE MOD_ReadInTools
USE MOD_DSMC_Vars           ,ONLY: useDSMC, DSMC, RadialWeighting
USE MOD_part_pressure       ,ONLY: ParticleInsideCheck
USE MOD_Dielectric_Vars     ,ONLY: DoDielectric,isDielectricElem,DielectricNoParticles
USE MOD_part_emission_tools ,ONLY: SetParticleChargeAndMass,SetParticleMPF
USE MOD_part_pos_and_velo   ,ONLY: SetParticlePosition,SetParticleVelocity
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i, NbrOfParticle,iInit,iPart,PositionNbr
INTEGER               :: nPartInside
INTEGER(KIND=8)       :: insertParticles
REAL                  :: EInside,TempInside
LOGICAL               :: EmType6
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)') ' Initial particle inserting... '

CALL UpdateNextFreePosition()
EmType6=.false.
DO i=1, nSpecies
  DO iInit = Species(i)%StartnumberOfInits, Species(i)%NumberOfInits
    IF ((Species(i)%Init(iInit)%ParticleEmissionType.EQ.6)) THEN
      EmType6=.true.
      EXIT
    END IF
  END DO
  IF (EmType6) EXIT
END DO
IF (.NOT.EmType6) DSMC%OutputMeshSamp=.false.
!   CALL Deposition()
!   IF (MESH%t.GE.PIC%DelayTime) PIC%ParticleTreatmentMethod='standard'
  ! for the case of particle insertion per time, the inserted particle number for the current time must
  ! be updated. Otherwise, at the first timestep after restart, these particles will be inserted again
!  DO i=1,nSpecies
!    Species(i)%InsertedParticle = INT(Species(i)%ParticleEmission * Time)
!  END DO
!ELSE
! Do sanity check of max. particle number compared to the number that is to be inserted for certain insertion types
insertParticles = 0
DO i=1,nSpecies
  IF (DoRestart .AND. .NOT.SpecReset(i)) CYCLE
  DO iInit = Species(i)%StartnumberOfInits, Species(i)%NumberOfInits
    IF (TRIM(Species(i)%Init(iInit)%SpaceIC).EQ.'cell_local') THEN
      IF(Symmetry%Order.LE.2) THEN
        ! The correct 2D/axisymmetric LocalVolume could only be calculated after the symmetry axis was defined (through the boundary
        ! conditions). However, the initialParticleNumber was already determined before the 2D volume calculation was performed.
        ! This can lead to initialParticleNumbers of 0, thus skipping the insertion entirely.
        Species(i)%Init(iInit)%initialParticleNumber &
                  = NINT(Species(i)%Init(iInit)%PartDensity / Species(i)%MacroParticleFactor * GEO%LocalVolume)
        ! The radial scaling of the weighting factor has to be considered
        IF(RadialWeighting%DoRadialWeighting) Species(i)%Init(iInit)%initialParticleNumber = &
                                    INT(Species(i)%Init(iInit)%initialParticleNumber * 2. / (RadialWeighting%PartScaleFactor),8)
      END IF
      IF (Species(i)%Init(iInit)%PartDensity.EQ.0) THEN
#if USE_MPI
        insertParticles = insertParticles + INT(REAL(Species(i)%Init(iInit)%initialParticleNumber)/PartMPI%nProcs,8)
#else
        insertParticles = insertParticles + INT(Species(i)%Init(iInit)%initialParticleNumber,8)
#endif
      ELSE
        insertParticles = insertParticles + INT(Species(i)%Init(iInit)%initialParticleNumber,8)
      END IF
    ELSE IF ((TRIM(Species(i)%Init(iInit)%SpaceIC).EQ.'cuboid') &
         .OR.(TRIM(Species(i)%Init(iInit)%SpaceIC).EQ.'cylinder')) THEN
#if USE_MPI
      insertParticles = insertParticles + INT(REAL(Species(i)%Init(iInit)%initialParticleNumber)/PartMPI%nProcs)
#else
      insertParticles = insertParticles + INT(Species(i)%Init(iInit)%initialParticleNumber,8)
#endif
    END IF
  END DO
END DO
IF (insertParticles.GT.PDM%maxParticleNumber) THEN
  IPWRITE(UNIT_stdOut,*)' Maximum particle number : ',PDM%maxParticleNumber
  IPWRITE(UNIT_stdOut,*)' To be inserted particles: ',INT(insertParticles,4)
  CALL abort(&
__STAMP__&
,'Number of to be inserted particles per init-proc exceeds max. particle number! ')
END IF
DO i = 1,nSpecies
  IF (DoRestart .AND. .NOT.SpecReset(i)) CYCLE
  DO iInit = Species(i)%StartnumberOfInits, Species(i)%NumberOfInits
    ! check whether initial particles are defined twice (old and new method) to prevent erroneous doubling
    ! of particles
    !!!Here could be added a check for geometrically overlapping Inits and same Usefor-Flags!!!
    !IF ((Species(i)%initialParticleNumber.NE.0).AND.(Species(i)%NumberOfInits.NE.0)) THEN
    !  WRITE(*,*) 'ERROR in ParticleEmission: Initial emission may only be defined in additional *Init#* blocks'
    !  WRITE(*,*) 'OR the standard initialisation, not both!'
    !  STOP
    !END IF
    IF (((Species(i)%Init(iInit)%ParticleEmissionType .EQ. 4).OR.(Species(i)%Init(iInit)%ParticleEmissionType .EQ. 6)) .AND. &
         (Species(i)%Init(iInit)%UseForInit)) THEN ! Special emission type: constant density in cell, + to be used for init
      CALL abort(&
__STAMP__&
,' particle pressure not moved to picasso!')
      IF (Species(i)%Init(iInit)%ParticleEmissionType .EQ. 4) THEN
        CALL ParticleInsertingCellPressure(i,iInit,NbrofParticle)
        CALL SetParticleVelocity(i,iInit,NbrOfParticle,1)
      ELSE !emission type 6 (constant pressure outflow)
        CALL ParticleInsertingPressureOut(i,iInit,NbrofParticle)
      END IF
      CALL SetParticleChargeAndMass(i,NbrOfParticle)
        IF (usevMPF.AND.(.NOT.RadialWeighting%DoRadialWeighting)) CALL SetParticleMPF(i,NbrOfParticle)
      IF (useDSMC) THEN
        IF(NbrOfParticle.gt.PDM%maxParticleNumber)THEN
          NbrOfParticle = PDM%maxParticleNumber
        END IF
        iPart = 1
        DO WHILE (iPart .le. NbrOfParticle)
          PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
          IF (PositionNbr .ne. 0) THEN
            PDM%PartInit(PositionNbr) = iInit
          END IF
          iPart = iPart + 1
        END DO
      END IF
      !IF (useDSMC) CALL SetParticleIntEnergy(i,NbrOfParticle)
      PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfParticle
      CALL UpdateNextFreePosition()
    ELSE IF (Species(i)%Init(iInit)%UseForInit) THEN ! no special emissiontype to be used
      IF(Species(i)%Init(iInit)%initialParticleNumber.GT.HUGE(1)) CALL abort(&
__STAMP__&
,' Integer of initial particle number larger than max integer size: ',HUGE(1))
      NbrOfParticle = INT(Species(i)%Init(iInit)%initialParticleNumber,4)
      SWRITE(UNIT_stdOut,'(A,I0,A)') ' Set particle position for species ',i,' ... '
#if USE_MPI
      CALL SetParticlePosition(i,iInit,NbrOfParticle,1)
      CALL SetParticlePosition(i,iInit,NbrOfParticle,2)
#else
      CALL SetParticlePosition(i,iInit,NbrOfParticle)
#endif /*USE_MPI*/
      SWRITE(UNIT_stdOut,'(A,I0,A)') ' Set particle velocities for species ',i,' ... '
      CALL SetParticleVelocity(i,iInit,NbrOfParticle,1)
      SWRITE(UNIT_stdOut,'(A,I0,A)') ' Set particle charge and mass for species ',i,' ... '
      CALL SetParticleChargeAndMass(i,NbrOfParticle)
      IF (usevMPF.AND.(.NOT.RadialWeighting%DoRadialWeighting)) CALL SetParticleMPF(i,NbrOfParticle)
      IF (useDSMC) THEN
        IF(NbrOfParticle.gt.PDM%maxParticleNumber)THEN
          NbrOfParticle = PDM%maxParticleNumber
        END IF
        iPart = 1
        DO WHILE (iPart .le. NbrOfParticle)
          PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
          IF (PositionNbr .ne. 0) THEN
            PDM%PartInit(PositionNbr) = iInit
          END IF
          iPart = iPart + 1
        END DO
      END IF
      !IF (useDSMC) CALL SetParticleIntEnergy(i,NbrOfParticle)
      PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfParticle
      CALL UpdateNextFreePosition()
      ! constant pressure condition
      IF ((Species(i)%Init(iInit)%ParticleEmissionType .EQ. 3).OR.(Species(i)%Init(iInit)%ParticleEmissionType .EQ. 5)) THEN
        CALL abort(&
__STAMP__&
,' particle pressure not moved in picasso!')
        CALL ParticleInsideCheck(i, iInit, nPartInside, TempInside, EInside)
        IF (Species(i)%Init(iInit)%ParticleEmission .GT. nPartInside) THEN
          NbrOfParticle = INT(Species(i)%Init(iInit)%ParticleEmission) - nPartInside
          IPWRITE(UNIT_stdOut,*) 'Emission PartNum (Spec ',i,')', NbrOfParticle
#if USE_MPI
          CALL SetParticlePosition(i,iInit,NbrOfParticle,1)
          CALL SetParticlePosition(i,iInit,NbrOfParticle,2)
#else
          CALL SetParticlePosition(i,iInit,NbrOfParticle)
#endif
          CALL SetParticleVelocity(i,iInit,NbrOfParticle,1)
          CALL SetParticleChargeAndMass(i,NbrOfParticle)
            IF (usevMPF.AND.(.NOT.RadialWeighting%DoRadialWeighting)) CALL SetParticleMPF(i,NbrOfParticle)
          !IF (useDSMC) CALL SetParticleIntEnergy(i,NbrOfParticle)
          PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfParticle
          CALL UpdateNextFreePosition()
        END IF
      END IF
    END IF ! not Emissiontype 4
  END DO !inits
END DO ! species

!--- set last element to current element (needed when ParticlePush is not executed, e.g. "delay")
DO i = 1,PDM%ParticleVecLength
  PEM%lastElement(i) = PEM%Element(i)
END DO

!--- Remove particles from dielectric regions if DielectricNoParticles=.TRUE.
IF(DoDielectric)THEN
  IF(DielectricNoParticles)THEN
    DO i = 1,PDM%ParticleVecLength
      ! Remove particles in dielectric elements
      IF(isDielectricElem(PEM%Element(i)))THEN
        PDM%ParticleInside(i) = .FALSE.
      END IF
    END DO
  END IF
END IF

SWRITE(UNIT_stdOut,'(A)') ' ...DONE '

END SUBROUTINE InitializeParticleEmission

#if USE_MPI
SUBROUTINE ParticleInserting(mode_opt)
#else
SUBROUTINE ParticleInserting()
#endif
!===================================================================================================================================
! Particle Inserting
!===================================================================================================================================
! Modules
#if USE_MPI
USE MOD_Particle_MPI_Vars,     ONLY : PartMPI
#endif /*USE_MPI*/
USE MOD_Globals
USE MOD_Timedisc_Vars         , ONLY : dt,time
USE MOD_Timedisc_Vars          ,ONLY: RKdtFrac,RKdtFracTotal
USE MOD_Particle_Vars
USE MOD_PIC_Vars
USE MOD_part_tools             ,ONLY : UpdateNextFreePosition
USE MOD_DSMC_Vars              ,ONLY : useDSMC, CollisMode, SpecDSMC, RadialWeighting
USE MOD_DSMC_Init              ,ONLY : DSMC_SetInternalEnr_LauxVFD
USE MOD_DSMC_PolyAtomicModel   ,ONLY : DSMC_SetInternalEnr_Poly
USE MOD_Particle_Analyze_Vars  ,ONLY: CalcPartBalance,nPartIn,PartEkinIn
USE MOD_Particle_Analyze_Tools ,ONLY: CalcEkinPart
USE MOD_part_pressure          ,ONLY: ParticlePressure, ParticlePressureRem
USE MOD_part_emission_tools    ,ONLY: SetParticleChargeAndMass,SetParticleMPF,SamplePoissonDistri
USE MOD_part_pos_and_velo      ,ONLY: SetParticlePosition,SetParticleVelocity
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
#if USE_MPI
INTEGER, OPTIONAL                :: mode_opt
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Local variable declaration
INTEGER                          :: i , iPart, PositionNbr, iInit, IntSample
INTEGER                , SAVE    :: NbrOfParticle=0
INTEGER(KIND=8)                  :: inserted_Particle_iter,inserted_Particle_time
INTEGER(KIND=8)                  :: inserted_Particle_diff
REAL                             :: PartIns, RandVal1
REAL                             :: RiseFactor, RiseTime
#if USE_MPI
INTEGER                          :: mode
INTEGER                          :: InitGroup
#endif
!===================================================================================================================================
#if USE_MPI
IF (PRESENT(mode_opt)) THEN
  mode=mode_opt
ELSE
  mode=0
END IF
#endif
!---  Emission at time step (initial emission see particle_init.f90: InitializeParticleEmission)
DO i=1,nSpecies
  DO iInit = Species(i)%StartnumberOfInits, Species(i)%NumberOfInits
    IF (((Species(i)%Init(iInit)%ParticleEmissionType .NE. 4).AND.(Species(i)%Init(iInit)%ParticleEmissionType .NE. 6)) .AND. &
         (Species(i)%Init(iInit)%UseForEmission)) THEN ! no constant density in cell type, + to be used for init
#if USE_MPI
      IF (mode.NE.2) THEN
#endif
        SELECT CASE(Species(i)%Init(iInit)%ParticleEmissionType)
        CASE(1) ! Emission Type: Particles per !!!!!SECOND!!!!!!!! (not per ns)
          IF (Species(i)%Init(iInit)%VirtPreInsert .AND. Species(i)%Init(iInit)%PartDensity.GT.0.) THEN
            PartIns=Species(i)%Init(iInit)%ParticleEmission * dt*RKdtFrac  ! emitted particles during time-slab
            NbrOfParticle = 0 ! calculated within SetParticlePosition itself!
          ELSE IF (.NOT.DoPoissonRounding .AND. .NOT.DoTimeDepInflow) THEN
            PartIns=Species(i)%Init(iInit)%ParticleEmission * dt*RKdtFrac  ! emitted particles during time-slab
            inserted_Particle_iter = INT(PartIns,8)                                     ! number of particles to be inserted
            PartIns=Species(i)%Init(iInit)%ParticleEmission * (Time + dt*RKdtFracTotal) ! total number of emitted particle over
                                                                                        ! simulation
            CALL RANDOM_NUMBER(RandVal1)
            !-- random-round the inserted_Particle_time for preventing periodicity
            ! PO & SC: why, sometimes we do not want this add, TB is bad!
            IF (inserted_Particle_iter.GE.1) THEN
              CALL RANDOM_NUMBER(RandVal1)
              inserted_Particle_time = INT(PartIns + RandVal1,8) ! adds up to ONE
            ELSE IF((inserted_Particle_iter.GE.0).AND.(inserted_Particle_iter.LT.1)) THEN
                                                       !needed, since InsertedParticleSurplus can increase
                                                       !and _iter>1 needs to be possible for preventing periodicity
              IF (ALMOSTEQUAL(PartIns,0.)) THEN !dummy
                inserted_Particle_time = INT(PartIns,8)
              ELSE !poisson-distri of PartIns-INT(PartIns)
                CALL SamplePoissonDistri( PartIns-INT(PartIns) , IntSample )
                inserted_Particle_time = INT(INT(PartIns)+IntSample,8) !INT(PartIns) + POISDISTRI( PartIns-INT(PartIns) )
              END IF
            ELSE !dummy
              inserted_Particle_time = INT(PartIns,8)
            END IF
            !-- evaluate inserted_Particle_time and inserted_Particle_iter
            inserted_Particle_diff = inserted_Particle_time - Species(i)%Init(iInit)%InsertedParticle &
              - inserted_Particle_iter - Species(i)%Init(iInit)%InsertedParticleSurplus &
              + Species(i)%Init(iInit)%InsertedParticleMisMatch
            Species(i)%Init(iInit)%InsertedParticleSurplus = ABS(MIN(inserted_Particle_iter + inserted_Particle_diff,0))
            NbrOfParticle = MAX(INT(inserted_Particle_iter + inserted_Particle_diff,4),0)
            !-- if maxwell velo dist and less than 5 parts: skip (to ensure maxwell dist)
            IF (TRIM(Species(i)%Init(iInit)%velocityDistribution).EQ.'maxwell') THEN
              IF (NbrOfParticle.LT.5) NbrOfParticle=0
            END IF
          ELSE IF (DoPoissonRounding .AND. .NOT.DoTimeDepInflow) THEN
            ! linear rise of inflow
            RiseTime=Species(i)%Init(iInit)%InflowRiseTime
            IF(RiseTime.GT.0.)THEN
              IF(Time-DelayTime.LT.RiseTime)THEN
                RiseFactor=(time-DelayTime)/RiseTime
              ELSE
                RiseFactor=1.
              END IF
            ELSE
              RiseFactor=1.
            EnD IF
            PartIns=Species(i)%Init(iInit)%ParticleEmission * dt*RKdtFrac * RiseFactor  ! emitted particles during time-slab
            CALL RANDOM_NUMBER(RandVal1)
            IF (EXP(-PartIns).LE.TINY(PartIns)) THEN
              IPWRITE(*,*)'WARNING: target is too large for poisson sampling: switching now to Random rounding...'
              NbrOfParticle = INT(PartIns + RandVal1)
              DoPoissonRounding = .FALSE.
            ELSE !poisson-sampling instead of random rounding (reduces numerical non-equlibrium effects [Tysanner and Garcia 2004]
              CALL SamplePoissonDistri( PartIns , NbrOfParticle , DoPoissonRounding)
            END IF
          ELSE ! DoTimeDepInflow
            ! linear rise of inflow
            RiseTime=Species(i)%Init(iInit)%InflowRiseTime
            IF(RiseTime.GT.0.)THEN
              IF(Time-DelayTime.LT.RiseTime)THEN
                RiseFactor=(time-DelayTime)/RiseTime
              ELSE
                RiseFactor=1.
              END IF
            ELSE
              RiseFactor=1.
            EnD IF
            ! emitted particles during time-slab
            PartIns=Species(i)%Init(iInit)%ParticleEmission * dt*RKdtFrac * RiseFactor &
                   + Species(i)%Init(iInit)%InsertedParticleMisMatch
            CALL RANDOM_NUMBER(RandVal1)
            NbrOfParticle = INT(PartIns + RandVal1)
          END IF
#if USE_MPI
          InitGroup=Species(i)%Init(iInit)%InitCOMM
          IF(PartMPI%InitGroup(InitGroup)%COMM.NE.MPI_COMM_NULL) THEN
            ! only procs which are part of group take part in the communication
             !NbrOfParticle based on RandVals!
            CALL MPI_BCAST(NbrOfParticle, 1, MPI_INTEGER,0,PartMPI%InitGroup(InitGroup)%COMM,IERROR)
          ELSE
            NbrOfParticle=0
          END IF
          !CALL MPI_BCAST(NbrOfParticle, 1, MPI_INTEGER,0,PartMPI%COMM,IERROR) !NbrOfParticle based on RandVals!
#endif
          Species(i)%Init(iInit)%InsertedParticle = Species(i)%Init(iInit)%InsertedParticle + INT(NbrOfParticle,8)
        CASE(2)    ! Emission Type: Particles per Iteration
          IF (RKdtFracTotal .EQ. 1.) THEN !insert in last stage only, so that no reconstruction is nec. and number/iter matches
            NbrOfParticle = INT(Species(i)%Init(iInit)%ParticleEmission)
          ELSE
            NbrOfParticle = 0
          END IF
        CASE(3)
          CALL abort(&
__STAMP__&
,' particle pressure not moved in picasso!')
          CALL ParticlePressure (i, iInit, NbrOfParticle)
          ! if maxwell velo dist and less than 5 parts: skip (to ensure maxwell dist)
          IF (TRIM(Species(i)%Init(iInit)%velocityDistribution).EQ.'maxwell') THEN
            IF (NbrOfParticle.LT.5) NbrOfParticle=0
          END IF
        CASE(5) ! removal of all parts in pressure area and re-insertion
          CALL abort(&
__STAMP__&
,' particle pressure not moved in picasso!')
          CALL ParticlePressureRem (i, iInit, NbrOfParticle)
        CASE DEFAULT
          NbrOfParticle = 0
        END SELECT
#if USE_MPI
        CALL SetParticlePosition(i,iInit,NbrOfParticle,1)
      END IF
      IF (mode.NE.1) THEN
        CALL SetParticlePosition(i,iInit,NbrOfParticle,2)
#else
        CALL SetParticlePosition(i,iInit,NbrOfParticle)
#endif
       CALL SetParticleVelocity(i,iInit,NbrOfParticle,1)
       CALL SetParticleChargeAndMass(i,NbrOfParticle)
       IF (usevMPF.AND.(.NOT.RadialWeighting%DoRadialWeighting)) CALL SetParticleMPF(i,NbrOfParticle)
       ! define molecule stuff
       IF (useDSMC.AND.(CollisMode.GT.1)) THEN
         iPart = 1
         DO WHILE (iPart .le. NbrOfParticle)
           PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
           IF (PositionNbr .ne. 0) THEN
             IF (SpecDSMC(i)%PolyatomicMol) THEN
               CALL DSMC_SetInternalEnr_Poly(i,iInit,PositionNbr,1)
             ELSE
               CALL DSMC_SetInternalEnr_LauxVFD(i,iInit,PositionNbr,1)
             END IF
           END IF
           iPart = iPart + 1
         END DO
       END IF
       ! instead of UpdateNextfreePosition we update the
       ! particleVecLength only.
       ! and doing it later, after calcpartbalance
       PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + NbrOfParticle
       PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfParticle
       !CALL UpdateNextFreePosition()
#if USE_MPI
      END IF
#endif
    ELSE IF (Species(i)%Init(iInit)%UseForEmission) THEN ! Constant Pressure in Cell Emission (type 4 or 6)
      IF (Species(i)%Init(iInit)%ParticleEmissionType .EQ. 4) THEN
        CALL ParticleInsertingCellPressure(i,iInit,NbrofParticle)
        CALL SetParticleVelocity(i,iInit,NbrOfParticle,1)
      ELSE !emission type 6 (constant pressure outflow)
        CALL abort(&
__STAMP__&
,' particle pressure not moved in picasso!')
        CALL ParticleInsertingPressureOut(i,iInit,NbrofParticle)
      END IF
      CALL SetParticleChargeAndMass(i,NbrOfParticle)
      IF (usevMPF.AND.(.NOT.RadialWeighting%DoRadialWeighting)) CALL SetParticleMPF(i,NbrOfParticle)
      ! define molecule stuff
      IF (useDSMC.AND.(CollisMode.GT.1)) THEN
        iPart = 1
        DO WHILE (iPart .le. NbrOfParticle)
          PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
          IF (PositionNbr .ne. 0) THEN
            IF (SpecDSMC(i)%PolyatomicMol) THEN
              CALL DSMC_SetInternalEnr_Poly(i,iInit,PositionNbr,1)
            ELSE
              CALL DSMC_SetInternalEnr_LauxVFD(i,iInit,PositionNbr,1)
            END IF
        END IF
          iPart = iPart + 1
        END DO
      END IF
      ! instead of UpdateNextfreePosition we update the
      ! particleVecLength only.
      ! and doing it after calcpartbalance
      PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + NbrOfParticle
      PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfParticle
      !CALL UpdateNextFreePosition()
    END IF
    ! compute number of input particles and energy
    IF(CalcPartBalance) THEN
      ! alter history, dirty hack for balance calculation
      PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition - NbrOfParticle
      IF(NbrOfParticle.GT.0)THEN
        nPartIn(i)=nPartIn(i) + NBrofParticle
        DO iPart=1,NbrOfparticle
          PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
          IF (PositionNbr .ne. 0) PartEkinIn(PartSpecies(PositionNbr)) = &
                                  PartEkinIn(PartSpecies(PositionNbr))+CalcEkinPart(PositionNbr)
        END DO ! iPart
      END IF
      ! alter history, dirty hack for balance calculation
      PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + NbrOfParticle
    END IF ! CalcPartBalance
  END DO  ! iInit = Species(i)%StartnumberOfInits, Species(i)%NumberOfInits
END DO  ! i=1,nSpecies

END SUBROUTINE ParticleInserting


SUBROUTINE ParticleInsertingCellPressure(iSpec,iInit,NbrOfParticle)
!===================================================================================================================================
! Insert constant cell pressure particles (and remove additionals)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars
USE MOD_Mesh_Vars,              ONLY:NGeo,XCL_NGeo,XiCL_NGeo,wBaryCL_NGeo
USE MOD_Particle_Mesh_Vars,     ONLY:GEO
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping,TriaTracking
USE MOD_Particle_Mesh,          ONLY:SingleParticleToExactElement,SingleParticleToExactElementNoMap
USE MOD_Eval_xyz,               ONLY:TensorProductInterpolation
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)    :: iSpec, iInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)   :: NbrOfParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iElem, Elem, iPart, i, NbrPartsInCell, NbrNewParts
INTEGER               :: ParticleIndexNbr
REAL                  :: PartDiff, PartDiffRest, RandVal, RandVal3(1:3)
INTEGER, ALLOCATABLE  :: PartsInCell(:)
!===================================================================================================================================

NbrOfParticle = 0
DO iElem = 1,Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside
  Elem = Species(iSpec)%Init(iInit)%ConstPress%ElemTotalInside(iElem)
  ! step 1: count and build array of particles in cell (of current species only)
  ALLOCATE(PartsInCell(1:PEM%pNumber(Elem)))
  NbrPartsInCell = 0
  iPart = PEM%pStart(Elem)
  DO i = 1, PEM%pNumber(Elem)
    IF (PartSpecies(iPart).EQ.iSpec) THEN
      NbrPartsInCell = NbrPartsInCell + 1
      PartsInCell(NbrPartsInCell) = iPart
    END IF
    iPart = PEM%pNext(iPart)
  END DO
  ! step 2: determine number of particles to insert (or remove)
  PartDiff = Species(iSpec)%Init(iInit)%ParticleEmission * GEO%Volume(Elem) - NbrPartsInCell
  PartDiffRest = PartDiff - INT(PartDiff)
  ! step 3: if PartDiff positive (and PartPressAddParts=T), add particles
  IF(PartPressAddParts.AND.PartDiff.GT.0) THEN
    CALL RANDOM_NUMBER(RandVal)
    IF(PartDiffRest.GT.RandVal) PartDiff = PartDiff + 1.0
    NbrNewParts = INT(PartDiff)
    ! insert particles (positions)
    DO i = 1, NbrNewParts
      ! set random position in -1,1 space
      CALL RANDOM_NUMBER(RandVal3)
      RandVal3 = RandVal3 * 2.0 - 1.0
      ParticleIndexNbr = PDM%nextFreePosition(PDM%CurrentNextFreePosition + i + NbrOfParticle)
      IF (ParticleIndexNbr.NE.0) THEN
        CALL TensorProductInterpolation(RandVal3,3,NGeo,XiCL_NGeo,wBaryCL_NGeo,&
                           XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,iElem),PartState(1:3,ParticleIndexNbr))
        !PartState(1:3,ParticleIndexNbr) = MapToGeo(RandVal3,P)
        PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
        IF (.NOT. DoRefMapping) THEN
          IF (TriaTracking) THEN
            CALL SingleParticleToExactElement(ParticleIndexNbr,doHALO=.FALSE.,initFIX=.FALSE.,doRelocate=.FALSE.)
          ELSE
            CALL SingleParticleToExactElementNoMap(ParticleIndexNbr,doHALO=.FALSE.,doRelocate=.FALSE.)
          END IF
        ELSE
          PartPosRef(1:3,ParticleIndexNbr)=RandVal3
        END IF
        IF(.NOT.PDM%ParticleInside(ParticleIndexNbr))THEN
          CALL abort(&
__STAMP__&
,' Particle lost in own MPI region. Need to communicate!')
        END IF
        IF (PDM%ParticleInside(ParticleIndexNbr)) PDM%IsNewPart(ParticleIndexNbr)=.TRUE.
      ELSE
        CALL abort(&
__STAMP__&
,'ERROR in ParticleInsertingCellPressure: ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?')
      END IF
    END DO
    NbrOfParticle = NbrOfParticle + NbrNewParts
  END IF
  ! step 4: if PartDiff negative (and PartPressRemParts=T), remove particles
  IF(PartPressRemParts.AND.PartDiff.LT.0) THEN
    PartDiff = -PartDiff
    CALL RANDOM_NUMBER(RandVal)
    IF(ABS(PartDiffRest).GT.RandVal) PartDiff = PartDiff + 1.0
    NbrNewParts = INT(PartDiff)
    ! remove random part
    DO i = 1, NbrNewParts
      CALL RANDOM_NUMBER(RandVal)
      RandVal = RandVal * REAL(NbrPartsInCell)
      PDM%ParticleInside(PartsInCell(INT(RandVal)+1)) = .FALSE.
    END DO
  END IF
  DEALLOCATE(PartsInCell)
END DO
END SUBROUTINE ParticleInsertingCellPressure

SUBROUTINE ParticleInsertingPressureOut(iSpec,iInit,NbrOfParticle)
!===================================================================================================================================
! Insert constant outflow pressure particles (copied mostly from 'ParticleInsertingCellPressure')
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars
USE MOD_Mesh_Vars,              ONLY:NGeo,XCL_NGeo,XiCL_NGeo,wBaryCL_NGeo
USE MOD_Particle_Mesh,          ONLY:SingleParticleToExactElement,SingleParticleToExactElementNoMap
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping,TriaTracking
USE MOD_Eval_xyz,               ONLY:TensorProductInterpolation
USE MOD_DSMC_Vars,              ONLY:CollisMode
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: iSpec, iInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)           :: NbrOfParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iElem, Elem, iPart, i, NbrPartsInCell,  distnum
INTEGER                       :: ParticleIndexNbr
REAL                          :: RandVal3(1:3)
INTEGER, ALLOCATABLE          :: PartsInCell(:)
REAL                          :: Velo1, Velo2, Velosq, v_sum(3), v2_sum, maxwellfac
REAL                          :: Vec3D(3), RandVal3D(3)
!===================================================================================================================================

NbrOfParticle = 0
IF (CollisMode.EQ.0) THEN
  CALL Abort(&
__STAMP__&
,"Free Molecular Flow (CollisMode=0) is not supported for const pressure outflow BC!")
END IF
IF (TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).NE.'maxwell') THEN
  CALL abort(&
__STAMP__&
,'Only maxwell implemented yet for const pressure outflow BC!')
END IF
DO iElem = 1,Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside
  Elem = Species(iSpec)%Init(iInit)%ConstPress%ElemTotalInside(iElem)
  ! step 1: count and build array of particles in cell (of current species only)
  ALLOCATE(PartsInCell(1:PEM%pNumber(Elem)))
  NbrPartsInCell = 0
  iPart = PEM%pStart(Elem)
  DO i = 1, PEM%pNumber(Elem)
    IF (PartSpecies(iPart).EQ.iSpec) THEN
      NbrPartsInCell = NbrPartsInCell + 1
      PartsInCell(NbrPartsInCell) = iPart
    END IF
    iPart = PEM%pNext(iPart)
  END DO
  ! step 2: sample cell values, remove particles, calculate new NbrPartsInCell
  CALL ParticleInsertingPressureOut_Sampling(iSpec,iInit,Elem,iElem,NbrPartsInCell,PartsInCell)
  DEALLOCATE(PartsInCell)
  ! step 3: add new particles
  IF(NbrPartsInCell.GT.0) THEN
    ! insert particles (positions and velocities)
    v_sum(1:3) = 0.0
    v2_sum = 0.0
    DO i = 1, NbrPartsInCell
      ! set random position in -1,1 space
      CALL RANDOM_NUMBER(RandVal3)
      RandVal3 = RandVal3 * 2.0 - 1.0
      ParticleIndexNbr = PDM%nextFreePosition(PDM%CurrentNextFreePosition + i + NbrOfParticle)
      IF (ParticleIndexNbr.NE.0) THEN
        CALL TensorProductInterpolation(RandVal3,3,NGeo,XiCL_NGeo,wBaryCL_NGeo,&
                           XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,iElem),PartState(1:3,ParticleIndexNbr))
        PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
        IF (.NOT. DoRefMapping) THEN
          IF (TriaTracking) THEN
            CALL SingleParticleToExactElement(ParticleIndexNbr,doHALO=.FALSE.,initFIX=.FALSE.,doRelocate=.FALSE.)
          ELSE
            CALL SingleParticleToExactElementNoMap(ParticleIndexNbr,doHALO=.FALSE.,doRelocate=.FALSE.)
          END IF
        ELSE
          PartPosRef(1:3,ParticleIndexNbr)=RandVal3
        END IF
        IF(.NOT.PDM%ParticleInside(ParticleIndexNbr))THEN
          CALL abort(&
__STAMP__&
,' Particle lost in own MPI region. Need to communicate!')
        END IF
        IF (PDM%ParticleInside(ParticleIndexNbr)) PDM%IsNewPart(ParticleIndexNbr)=.TRUE.
        ! Determine the particle velocity (maxwell, part 1)
        DO distnum = 1, 3
          CALL RANDOM_NUMBER(RandVal3D)
          Velo1 = 2.0*RandVal3D(1)-1.0
          Velo2 = 2.0*RandVal3D(2)-1.0
          Velosq= Velo1**2+Velo2**2
          DO WHILE ((Velosq.LE.0).OR.(Velosq.GE.1))
            CALL RANDOM_NUMBER(RandVal3D)
            Velo1 = 2.0*RandVal3D(1)-1.0
            Velo2 = 2.0*RandVal3D(2)-1.0
            Velosq= Velo1**2+Velo2**2
          END DO
          Vec3D(distnum) = Velo1*SQRT(-2*LOG(Velosq)/Velosq)
        END DO
        PartState(4:6,ParticleIndexNbr) = Vec3D(1:3)
        v_sum(1:3) = v_sum(1:3) + Vec3D(1:3)
        v2_sum = v2_sum + Vec3D(1)**2+Vec3D(2)**2+Vec3D(3)**2
      ELSE
        CALL abort(&
__STAMP__&
,'ERROR in ParticleInsertingCellPressureOut: ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?')
      END IF
    END DO
    v_sum(1:3) = v_sum(1:3) / (NbrPartsInCell+1) !+1 correct?
    v2_sum = v2_sum / (NbrPartsInCell+1)         !+1 correct?
    !maxwellfactor from new calculated values (no vibrational DOF implemented, equilibirium assumed)
    !Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(iElem,:)
    maxwellfac = SQRT(3. * Species(iSpec)%Init(iInit)%ConstantPressure &
                 / (Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(iElem,4)) / & ! velocity of maximum
                 (Species(iSpec)%MassIC*v2_sum))                                                           ! T = p_o / (<n>*k)
    ! particel velocity (maxwell, part 2)
    DO i = 1, NbrPartsInCell
      ParticleIndexNbr = PDM%nextFreePosition(PDM%CurrentNextFreePosition + i + NbrOfParticle)
      IF (ParticleIndexNbr .ne. 0) THEN
        PartState(4:6,ParticleIndexNbr) = (PartState(4:6,ParticleIndexNbr) - v_sum(1:3)) * maxwellfac &  !macro velocity:
                                                                                      !=vi + VeloVecIC*(<p>-p_o)/(SQRT(a**2)*<n>*mt)
             + Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(iElem,1:3) + Species(iSpec)%Init(iInit)%VeloVecIC(1:3) &
             * (Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(iElem,5) - Species(iSpec)%Init(iInit)%ConstantPressure) &
             / (SQRT(Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(iElem,6)) &
                * Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(iElem,4) * Species(iSpec)%MassIC)
      END IF
    END DO

    NbrOfParticle = NbrOfParticle + NbrPartsInCell
  END IF
END DO
END SUBROUTINE ParticleInsertingPressureOut

SUBROUTINE ParticleInsertingPressureOut_Sampling(iSpec, iInit, iElem, ElemSamp, NbrPartsInCell, PartsInCell)
!===================================================================================================================================
! Subroutine to sample current cell values (partly copied from 'LD_DSMC_Mean_Bufferzone_A_Val' and 'dsmc_analyze')
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,          ONLY : BoltzmannConst
USE MOD_Particle_Vars,         ONLY : PartState,usevMPF,Species,PartSpecies,usevMPF,PartMPF,PDM
USE MOD_DSMC_Vars,             ONLY : SpecDSMC
USE MOD_TimeDisc_Vars,         ONLY : iter
USE MOD_Particle_Mesh_Vars,    ONLY : GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)    :: iSpec, iInit, iElem, ElemSamp
INTEGER,INTENT(IN)    :: PartsInCell(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT) :: NbrPartsInCell
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iPart, iPartIndx
REAL                  :: MPFSum, WeightFak, kappa_part, AvogadroConst, RandVal, RealnumberNewParts
REAL                  :: Samp_V2(3), Samp_Temp(4), OldConstPressureSamp(6)
!===================================================================================================================================

IF (NbrPartsInCell .GT. 1) THEN ! Are there more than one particle
  IF(iter.EQ.0) THEN
    OldConstPressureSamp(:) = 0.0
  ELSE
    OldConstPressureSamp(:) = Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(ElemSamp,:)
  END IF
  Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(ElemSamp,:)        = 0.0
  MPFSum                            = 0.0
  Samp_V2(:)                        = 0.0
  Samp_Temp(:)                      = 0.0
  kappa_part                        = 0.0
  AvogadroConst                     = 6.02214129e23 ![1/mol]
  ! Loop over all particles of current species in cell
  DO iPart = 1, NbrPartsInCell
    iPartIndx = PartsInCell(iPart)
    IF (usevMPF) THEN
       WeightFak = PartMPF(iPartIndx)
    ELSE
       WeightFak = Species(iSpec)%MacroParticleFactor
    END IF
    Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(ElemSamp,1:3) &                !vi = vi + vi*w
         = Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(ElemSamp,1:3) &
         + PartState(4:6,iPartIndx) * WeightFak
    Samp_V2(:)                      = Samp_V2(:) + PartState(4:6,iPartIndx)**2 * WeightFak !vi**2 =vi**2 + vi**2*W
    MPFSum                          = MPFSum + WeightFak                                   !MPFsum = MPFsum + W
    PDM%ParticleInside(iPartIndx)=.false. !remove particle
  END DO

  !Calculation of specific heat ratio (no vibrational DOF -> only at low temperatures !!!)
  IF((SpecDSMC(PartSpecies(iPartIndx))%InterID.EQ.2).OR.(SpecDSMC(PartSpecies(iPartIndx))%InterID.EQ.20)) THEN
    kappa_part=1.4
  ELSE IF(SpecDSMC(PartSpecies(iPartIndx))%InterID.EQ.1) THEN
    kappa_part=5.0/3.0
  ELSE
    CALL abort(&
__STAMP__&
,'Wrong PartSpecies for outflow BC!')
  END IF
  ! Calculation of sampling values
  Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(ElemSamp,1:3) &
       = Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(ElemSamp,1:3) / MPFSum              !vi = vi / MPFsum
  Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(ElemSamp,4)   = MPFSum / GEO%Volume(iElem) !n = N / V
  Samp_Temp(1:3) &
       = Species(iSpec)%MassIC / BoltzmannConst * (Samp_V2(:) / MPFSum &                             !Ti = mt/k * (<vi**2>-<vi>**2)
       - Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(ElemSamp,1:3)**2)
  Samp_Temp(4) = (Samp_Temp(1) + Samp_Temp(2) + Samp_Temp(3)) / 3                                    !T = (Tx + Ty + Tz) / 3
  Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(ElemSamp,5) &                              !p = N / V * k * T
       = MPFSum / GEO%Volume(iElem) * BoltzmannConst * Samp_Temp(4)
  Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(ElemSamp,6) &                              !a**2 = kappa * k/mt * T
       = kappa_part * BoltzmannConst/Species(iSpec)%MassIC * Samp_Temp(4)


!----Ralaxationfaktor due to statistical noise in DSMC Results
  IF(iter.NE.0) THEN
    Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(ElemSamp,:) = (1.0 - Species(iSpec)%Init(iInit)%ConstPressureRelaxFac) &
                               * OldConstPressureSamp(:) + Species(iSpec)%Init(iInit)%ConstPressureRelaxFac &
                               * Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(ElemSamp,:)
  END IF
! Calculation of new density and resulting number in cell
  RealnumberNewParts = (Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(ElemSamp,4) & !N=(<n> + (p_o-<p>)/(a**2*mt)) * V/MPF
       + (Species(iSpec)%Init(iInit)%ConstantPressure - Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(ElemSamp,5)) &
       / (Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(ElemSamp,6) * Species(iSpec)%MassIC)) &
       * GEO%Volume(iElem) / Species(iSpec)%MacroParticleFactor !!!not sure if MPF treatment is correct!!!
  IF(RealnumberNewParts.GT.0.) THEN
    CALL RANDOM_NUMBER(RandVal)
    NbrPartsInCell = INT(RealnumberNewParts+RandVal)
  ELSE
    NbrPartsInCell = 0
  END IF

ELSE ! no particles in cell!
  CALL abort(&
__STAMP__&
,'YOU NEED MORE PARTICLES INSIDE THE OUTFLOW REGION!!!')
END IF

END SUBROUTINE ParticleInsertingPressureOut_Sampling


SUBROUTINE AdaptiveBCAnalyze(initSampling_opt)
!===================================================================================================================================
! Sampling of variables (part-density, velocity and energy) for Adaptive BC elements
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars              ,ONLY: PartStateIntEn, DSMC, CollisMode, SpecDSMC, useDSMC, RadialWeighting
USE MOD_Particle_Vars          ,ONLY: PartState, PDM, PartSpecies, Species, nSpecies, PEM, Adaptive_MacroVal, AdaptiveWeightFac
USE MOD_Particle_Vars          ,ONLY: usevMPF
USE MOD_Particle_Boundary_Vars ,ONLY: PorousBCSampIter, PorousBCMacroVal
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO,IsTracingBCElem
USE MOD_DSMC_Analyze           ,ONLY: CalcTVib,CalcTVibPoly,CalcTelec
USE MOD_Timedisc_Vars          ,ONLY: iter
USE MOD_part_tools             ,ONLY: GetParticleWeight
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers     ,ONLY: LBStartTime, LBElemSplitTime, LBPauseTime
USE MOD_LoadBalance_vars       ,ONLY: nPartsPerBCElem
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL, INTENT(IN), OPTIONAL   :: initSampling_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: ElemID, AdaptiveElemID, i, iSpec, SamplingIteration
REAL                            :: TVib_TempFac, TTrans_TempFac, RelaxationFactor
REAL, ALLOCATABLE               :: Source(:,:,:)
LOGICAL                         :: initSampling
#if USE_LOADBALANCE
REAL                            :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
ALLOCATE(Source(1:11,1:nElems,1:nSpecies))
Source=0.0

! Optional flag for the utilization of the routine for an initial sampling of the density and pressure distribution before simstart
IF(PRESENT(initSampling_opt)) THEN
  initSampling = initSampling_opt
ELSE
 initSampling = .FALSE.
END IF

! If no particles are present during the initial sampling for a porous BC, leave the routine, otherwise initial variables for the
! adaptive inlet surface flux will be overwritten by zero's.
IF (PDM%ParticleVecLength.LT.1) RETURN

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
DO i=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(i)) THEN
    ElemID = PEM%Element(i)
    IF(.NOT.IsTracingBCElem(ElemID))CYCLE
#if USE_LOADBALANCE
    nPartsPerBCElem(ElemID) = nPartsPerBCElem(ElemID) + 1
#endif /*USE_LOADBALANCE*/
    !ElemID = BC2AdaptiveElemMap(ElemID)
    iSpec = PartSpecies(i)
    Source(1:3,ElemID, iSpec) = Source(1:3,ElemID,iSpec) + PartState(4:6,i) * GetParticleWeight(i)
    Source(4:6,ElemID, iSpec) = Source(4:6,ElemID,iSpec) + PartState(4:6,i)**2 * GetParticleWeight(i)
    Source(7,ElemID, iSpec) = Source(7,ElemID, iSpec) + 1.0  ! simulation particle number
    IF(useDSMC)THEN
      IF ((CollisMode.EQ.2).OR.(CollisMode.EQ.3)) THEN
        IF ((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
          Source(8:9,ElemID, iSpec) = Source(8:9,ElemID, iSpec) + PartStateIntEn(1:2,i) * GetParticleWeight(i)
        END IF
      END IF
      IF (DSMC%ElectronicModel) THEN
        Source(10,ElemID, iSpec) = Source(10,ElemID, iSpec) + PartStateIntEn(3,i) * GetParticleWeight(i)
      END IF
    END IF
    Source(11,ElemID, iSpec) = Source(11,ElemID, iSpec) + GetParticleWeight(i)
  END IF
END DO
#if USE_LOADBALANCE
CALL LBPauseTime(LB_ADAPTIVE,tLBStart)
#endif /*USE_LOADBALANCE*/

IF(initSampling) THEN
  RelaxationFactor = 1
  SamplingIteration = 1
ELSE
  RelaxationFactor = AdaptiveWeightFac
  SamplingIteration = PorousBCSampIter
END IF

!DO iElem = 1,nElems
!IF(.NOT.IsTracingBCElem(iElem))CYCLE
DO AdaptiveElemID = 1,nElems
IF(.NOT.IsTracingBCElem(AdaptiveElemID))CYCLE
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
DO iSpec = 1,nSpecies
  ! write timesample particle values of bc elements in global macrovalues of bc elements
  IF (Source(7,AdaptiveElemID,iSpec).GT.0.0) THEN
    IF(.NOT.initSampling) THEN
      ! compute flow velocity (during computation, not for the initial distribution, where the velocity from the ini is used)
      Adaptive_MacroVal(1:3,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(1:3,AdaptiveElemID,iSpec) &
          + RelaxationFactor*Source(1:3,AdaptiveElemID, iSpec) / Source(11,AdaptiveElemID,iSpec)
    END IF
    ! compute flow Temperature
    IF (Source(7,AdaptiveElemID,iSpec).GT.1.0) THEN
      Adaptive_MacroVal(4:6,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(4:6,AdaptiveElemID,iSpec) &
        + RelaxationFactor &
          * (Source(7,AdaptiveElemID,iSpec)/(Source(7,AdaptiveElemID,iSpec)-1.0)) &
          * Species(iSpec)%MassIC/ BoltzmannConst &
          * ( Source(4:6,AdaptiveElemID,iSpec) / Source(11,AdaptiveElemID,iSpec) &
          - (Source(1:3,AdaptiveElemID,iSpec)/Source(11,AdaptiveElemID,iSpec))**2)
    ELSE
      Adaptive_MacroVal(4:6,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(4:6,AdaptiveElemID,iSpec)
    END IF
    ! ================================================================
    IF(PorousBCSampIter.GT.0) THEN
      ! Sampling the number density and pressure every given number of iterations and RESETTING it after calculation
      PorousBCMacroVal(1:6,AdaptiveElemID,iSpec) = PorousBCMacroVal(1:6,AdaptiveElemID,iSpec) + Source(1:6,AdaptiveElemID, iSpec)
      ! Sampling the particle weights
      PorousBCMacroVal(7,AdaptiveElemID,iSpec) = PorousBCMacroVal(7,AdaptiveElemID,iSpec) + Source(11,AdaptiveElemID,iSpec)
      ! Sampling the number of simulation particles
      PorousBCMacroVal(8,AdaptiveElemID,iSpec) = PorousBCMacroVal(8,AdaptiveElemID,iSpec) + Source(7,AdaptiveElemID,iSpec)
      IF(MOD(iter,SamplingIteration).EQ.0) THEN
        IF(PorousBCMacroVal(8,AdaptiveElemID,iSpec).GT.1) THEN
          ! number density
          IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
            Adaptive_MacroVal(7,AdaptiveElemID,iSpec)=PorousBCMacroVal(7,AdaptiveElemID,iSpec)/SamplingIteration   &
                                                      /GEO%Volume(AdaptiveElemID)
          ELSE
            Adaptive_MacroVal(7,AdaptiveElemID,iSpec)=PorousBCMacroVal(7,AdaptiveElemID,iSpec)/SamplingIteration   &
                                                      /GEO%Volume(AdaptiveElemID) * Species(iSpec)%MacroParticleFactor
          END IF
          ! instantaneous temperature WITHOUT 1/BoltzmannConst
          PorousBCMacroVal(1:6,AdaptiveElemID,iSpec) = PorousBCMacroVal(1:6,AdaptiveElemID,iSpec) &
                                                        / PorousBCMacroVal(7,AdaptiveElemID,iSpec)
          TTrans_TempFac = (PorousBCMacroVal(8,AdaptiveElemID,iSpec)/(PorousBCMacroVal(8,AdaptiveElemID,iSpec)-1.0)) &
              *Species(iSpec)%MassIC*(PorousBCMacroVal(4,AdaptiveElemID,iSpec) - PorousBCMacroVal(1,AdaptiveElemID,iSpec)**2   &
                                    + PorousBCMacroVal(5,AdaptiveElemID,iSpec) - PorousBCMacroVal(2,AdaptiveElemID,iSpec)**2   &
                                    + PorousBCMacroVal(6,AdaptiveElemID,iSpec) - PorousBCMacroVal(3,AdaptiveElemID,iSpec)**2) / 3.
          ! pressure (BoltzmannConstant canceled out in temperature calculation)
          Adaptive_MacroVal(12,AdaptiveElemID,iSpec)=Adaptive_MacroVal(7,AdaptiveElemID,iSpec)*TTrans_TempFac
          ! Resetting the sampling values
          PorousBCMacroVal(1:8,AdaptiveElemID,iSpec) = 0.0
        END IF
      END IF
    ELSE
      ! Calculation of the number density and pressure with the relaxation factor
      ! compute instantaneous temperature WITHOUT 1/BoltzmannConst
      IF (Source(7,AdaptiveElemID,iSpec).GT.1.0) THEN
        TTrans_TempFac = (Source(7,AdaptiveElemID,iSpec)/(Source(7,AdaptiveElemID,iSpec)-1.0)) &
                        * Species(iSpec)%MassIC * (Source(4,AdaptiveElemID,iSpec) / Source(11,AdaptiveElemID,iSpec)       &
                        - (Source(1,AdaptiveElemID,iSpec)/Source(11,AdaptiveElemID,iSpec))**2   &
                        + Source(5,AdaptiveElemID,iSpec) / Source(11,AdaptiveElemID,iSpec)      &
                        - (Source(2,AdaptiveElemID,iSpec)/Source(11,AdaptiveElemID,iSpec))**2   &
                        + Source(6,AdaptiveElemID,iSpec) / Source(11,AdaptiveElemID,iSpec)      &
                        - (Source(3,AdaptiveElemID,iSpec)/Source(11,AdaptiveElemID,iSpec))**2) / 3.
      ELSE
        TTrans_TempFac = 0.0
      END IF
      IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
        ! compute density
        Adaptive_MacroVal(7,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(7,AdaptiveElemID,iSpec) &
          + RelaxationFactor*Source(11,AdaptiveElemID,iSpec) /GEO%Volume(AdaptiveElemID)
        ! Pressure with relaxation factor
        Adaptive_MacroVal(12,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(12,AdaptiveElemID,iSpec) &
        +RelaxationFactor*Source(11,AdaptiveElemID,iSpec)/GEO%Volume(AdaptiveElemID)*TTrans_TempFac
      ELSE
        ! compute density
        Adaptive_MacroVal(7,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(7,AdaptiveElemID,iSpec) &
          + RelaxationFactor*Source(11,AdaptiveElemID,iSpec)/GEO%Volume(AdaptiveElemID)*Species(iSpec)%MacroParticleFactor
        ! Pressure with relaxation factor
        Adaptive_MacroVal(12,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(12,AdaptiveElemID,iSpec) &
                                                    + RelaxationFactor*Source(11,AdaptiveElemID,iSpec)  &
                                                    / GEO%Volume(AdaptiveElemID)*Species(iSpec)%MacroParticleFactor*TTrans_TempFac
      END IF
    END IF
    ! !==================================================================================================
    IF(useDSMC)THEN
      IF ((CollisMode.EQ.2).OR.(CollisMode.EQ.3))THEN
      IF ((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
          IF (DSMC%VibEnergyModel.EQ.0) THEN              ! SHO-model
            IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
              IF( (Source(8,AdaptiveElemID,iSpec)/Source(11,AdaptiveElemID,iSpec)) .GT. SpecDSMC(iSpec)%EZeroPoint) THEN
                Adaptive_MacroVal(8,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(8,AdaptiveElemID,iSpec) &
                  + RelaxationFactor*CalcTVibPoly(Source(8,AdaptiveElemID,iSpec) / Source(11,AdaptiveElemID,iSpec),iSpec)
              ELSE
                Adaptive_MacroVal(8,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(8,AdaptiveElemID,iSpec)
              END IF
            ELSE
              TVib_TempFac=Source(8,AdaptiveElemID,iSpec)/ (Source(11,AdaptiveElemID,iSpec) &
                *BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)
              IF (TVib_TempFac.LE.DSMC%GammaQuant) THEN
                Adaptive_MacroVal(8,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(8,AdaptiveElemID,iSpec)
              ELSE
                Adaptive_MacroVal(8,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(8,AdaptiveElemID,iSpec) &
                  + RelaxationFactor*SpecDSMC(iSpec)%CharaTVib / LOG(1 + 1/(TVib_TempFac-DSMC%GammaQuant))
              END IF
            END IF
          ELSE                                            ! TSHO-model
            Adaptive_MacroVal(8,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(8,AdaptiveElemID,iSpec) &
              + RelaxationFactor*CalcTVib(SpecDSMC(iSpec)%CharaTVib &
              , Source(8,AdaptiveElemID,iSpec)/Source(11,AdaptiveElemID,iSpec),SpecDSMC(iSpec)%MaxVibQuant)
          END IF
          Adaptive_MacroVal(9,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(9,AdaptiveElemID,iSpec) &
              + RelaxationFactor*Source(9,AdaptiveElemID,iSpec)/(Source(11,AdaptiveElemID,iSpec)*BoltzmannConst)
          IF (DSMC%ElectronicModel) THEN
            Adaptive_MacroVal(10,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(10,AdaptiveElemID,iSpec) &
              + RelaxationFactor*CalcTelec( Source(10,AdaptiveElemID,iSpec)/Source(11,AdaptiveElemID,iSpec),iSpec)
          END IF
        END IF
      END IF
    END IF
  ELSE
    Adaptive_MacroVal(1:10,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(1:10,AdaptiveElemID,iSpec)
  END IF
END DO
#if USE_LOADBALANCE
CALL LBElemSplitTime(AdaptiveElemID,tLBStart)
#endif /*USE_LOADBALANCE*/
END DO

END SUBROUTINE AdaptiveBCAnalyze


END MODULE MOD_part_emission
