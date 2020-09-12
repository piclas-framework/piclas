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
USE MOD_Dielectric_Vars     ,ONLY: DoDielectric,isDielectricElem,DielectricNoParticles
USE MOD_DSMC_Vars           ,ONLY: useDSMC, RadialWeighting
USE MOD_Part_Emission_Tools ,ONLY: SetParticleChargeAndMass,SetParticleMPF,SetParticleTimeStep
USE MOD_Part_Pos_and_Velo   ,ONLY: SetParticlePosition,SetParticleVelocity
USE MOD_Part_Tools          ,ONLY: UpdateNextFreePosition
USE MOD_Particle_Mesh_Vars  ,ONLY: LocalVolume
USE MOD_Particle_Vars       ,ONLY: Species,nSpecies,PDM,PEM, usevMPF, SpecReset, Symmetry, VarTimeStep
USE MOD_ReadInTools
USE MOD_Restart_Vars        ,ONLY: DoRestart
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i, NbrOfParticle,iInit,iPart,PositionNbr
INTEGER(KIND=8)       :: insertParticles
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)') ' Initial particle inserting... '

CALL UpdateNextFreePosition()

! Do sanity check of max. particle number compared to the number that is to be inserted for certain insertion types
insertParticles = 0
DO i=1,nSpecies
  IF (DoRestart .AND. .NOT.SpecReset(i)) CYCLE
  DO iInit = 1, Species(i)%NumberOfInits
    IF (TRIM(Species(i)%Init(iInit)%SpaceIC).EQ.'cell_local') THEN
      IF(Symmetry%Order.LE.2) THEN
        ! The correct 2D/axisymmetric LocalVolume could only be calculated after the symmetry axis was defined (through the boundary
        ! conditions). However, the initialParticleNumber was already determined before the 2D volume calculation was performed.
        ! This can lead to initialParticleNumbers of 0, thus skipping the insertion entirely.
        Species(i)%Init(iInit)%initialParticleNumber &
                  = NINT(Species(i)%Init(iInit)%PartDensity / Species(i)%MacroParticleFactor * LocalVolume)
        ! The radial scaling of the weighting factor has to be considered
        IF(RadialWeighting%DoRadialWeighting) Species(i)%Init(iInit)%initialParticleNumber = &
                                    INT(Species(i)%Init(iInit)%initialParticleNumber * 2. / (RadialWeighting%PartScaleFactor),8)
      END IF
#if USE_MPI
      insertParticles = insertParticles + INT(REAL(Species(i)%Init(iInit)%initialParticleNumber)/PartMPI%nProcs,8)
#else
      insertParticles = insertParticles + INT(Species(i)%Init(iInit)%initialParticleNumber,8)
#endif
    ELSE IF ((TRIM(Species(i)%Init(iInit)%SpaceIC).EQ.'cuboid') &
         .OR.(TRIM(Species(i)%Init(iInit)%SpaceIC).EQ.'cylinder')) THEN
#if USE_MPI
      insertParticles = insertParticles + INT(REAL(Species(i)%Init(iInit)%initialParticleNumber)/PartMPI%nProcs,8)
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
  DO iInit = 1, Species(i)%NumberOfInits
    IF (Species(i)%Init(iInit)%UseForInit) THEN ! no special emissiontype to be used
      IF(Species(i)%Init(iInit)%initialParticleNumber.GT.HUGE(1)) CALL abort(&
__STAMP__&
,' Integer of initial particle number larger than max integer size: ',HUGE(1))
      NbrOfParticle = INT(Species(i)%Init(iInit)%initialParticleNumber,4)
      SWRITE(UNIT_stdOut,'(A,I0,A)') ' Set particle position for species ',i,' ... '
      CALL SetParticlePosition(i,iInit,NbrOfParticle)
      SWRITE(UNIT_stdOut,'(A,I0,A)') ' Set particle velocities for species ',i,' ... '
      CALL SetParticleVelocity(i,iInit,NbrOfParticle)
      SWRITE(UNIT_stdOut,'(A,I0,A)') ' Set particle charge and mass for species ',i,' ... '
      CALL SetParticleChargeAndMass(i,NbrOfParticle)
      IF (usevMPF) CALL SetParticleMPF(i,NbrOfParticle)
      IF (VarTimeStep%UseVariableTimeStep) CALL SetParticleTimeStep(NbrOfParticle)
      IF (useDSMC) THEN
        iPart = 1
        DO WHILE (iPart .le. NbrOfParticle)
          PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
          IF (PositionNbr .ne. 0) THEN
            PDM%PartInit(PositionNbr) = iInit
          ELSE
            CALL abort(&
            __STAMP__&
            ,'ERROR in SetParticlePosition:ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?')    
          END IF
          iPart = iPart + 1
        END DO
      END IF
      !IF (useDSMC) CALL SetParticleIntEnergy(i,NbrOfParticle)
      PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfParticle
      CALL UpdateNextFreePosition()
    END IF ! not Emissiontype 4
  END DO !inits
END DO ! species

!--- set last element to current element (needed when ParticlePush is not executed, e.g. "delay")
DO i = 1,PDM%ParticleVecLength
  PEM%LastGlobalElemID(i) = PEM%GlobalElemID(i)
END DO

!--- Remove particles from dielectric regions if DielectricNoParticles=.TRUE.
IF(DoDielectric)THEN
  IF(DielectricNoParticles)THEN
    DO i = 1,PDM%ParticleVecLength
      ! Remove particles in dielectric elements
      IF(isDielectricElem(PEM%GlobalElemID(i)))THEN
        PDM%ParticleInside(i) = .FALSE.
      END IF
    END DO
  END IF
END IF

SWRITE(UNIT_stdOut,'(A)') ' ...DONE '

END SUBROUTINE InitializeParticleEmission

SUBROUTINE ParticleInserting()
!===================================================================================================================================
! Particle Inserting
!===================================================================================================================================
! Modules
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*USE_MPI*/
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: c
USE MOD_Timedisc_Vars          ,ONLY: dt,time
USE MOD_Timedisc_Vars          ,ONLY: RKdtFrac,RKdtFracTotal
USE MOD_Particle_Vars
USE MOD_PIC_Vars
USE MOD_part_tools             ,ONLY: UpdateNextFreePosition
USE MOD_DSMC_Vars              ,ONLY: useDSMC, CollisMode, SpecDSMC
USE MOD_DSMC_Init              ,ONLY: DSMC_SetInternalEnr_LauxVFD
USE MOD_DSMC_PolyAtomicModel   ,ONLY: DSMC_SetInternalEnr_Poly
USE MOD_Particle_Analyze_Vars  ,ONLY: CalcPartBalance,nPartIn,PartEkinIn
USE MOD_Particle_Analyze_Tools ,ONLY: CalcEkinPart
USE MOD_part_emission_tools    ,ONLY: SetParticleChargeAndMass,SetParticleMPF,SamplePoissonDistri,SetParticleTimeStep,CalcNbrOfPhotons
USE MOD_part_pos_and_velo      ,ONLY: SetParticlePosition,SetParticleVelocity
USE MOD_DSMC_BGGas             ,ONLY: BGGas_PhotoIonization
USE MOD_DSMC_ChemReact         ,ONLY: CalcPhotoIonizationNumber
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
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
REAL                             :: RiseFactor, RiseTime,NbrOfPhotons
#if USE_MPI
INTEGER                          :: InitGroup
#endif
REAL                             :: NbrOfReactions
!===================================================================================================================================

!---  Emission at time step (initial emission see particle_init.f90: InitializeParticleEmission)
DO i=1,nSpecies
  DO iInit = 1, Species(i)%NumberOfInits
    IF (Species(i)%Init(iInit)%UseForEmission) THEN ! no constant density in cell type, + to be used for init
        SELECT CASE(Species(i)%Init(iInit)%ParticleEmissionType)
        CASE(1) ! Emission Type: Particles per !!!!!SECOND!!!!!!!! (not per ns)
          IF (.NOT.DoPoissonRounding .AND. .NOT.DoTimeDepInflow) THEN
            PartIns=Species(i)%Init(iInit)%ParticleEmission * dt*RKdtFrac  ! emitted particles during time-slab
            inserted_Particle_iter = INT(PartIns,8)                                     ! number of particles to be inserted
            PartIns=Species(i)%Init(iInit)%ParticleEmission * (Time + dt*RKdtFracTotal) ! total number of emitted particle over
                                                                                        ! simulation
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
        CASE(7) ! SEE based on photon impact and photo-ionization in the volume
          ASSOCIATE( tShift => Species(i)%Init(iInit)%tShift )
            ! Check if all pulses have terminated
            IF(Time.LE.Species(i)%Init(iInit)%tActive)THEN
              ! Check if pulse is currently active of in between two pulses (in the latter case, do nothing)
              IF(MOD(MERGE(Time-tShift, Time, Time.GE.tShift), Species(i)%Init(iInit)%Period).LE.2.0*tShift)THEN
                ! Calculate the number of currently active photons (both surface SEE and volumetric emission)
                CALL CalcNbrOfPhotons(i, iInit, NbrOfPhotons)

                ! Check if only particles in the first quadrant are to be inserted that is spanned by the vectors 
                ! x=BaseVector1IC and y=BaseVector2IC in the interval x,y in [0,R] and reduce the number of photon accordingly
                IF(Species(i)%Init(iInit)%FirstQuadrantOnly) NbrOfPhotons = NbrOfPhotons / 4.0

                ! Select surface SEE or volumetric emission
                IF(TRIM(Species(i)%Init(iInit)%SpaceIC).EQ.'photon_SEE_disc')THEN
                  ! SEE based on photon impact
                  NbrOfPhotons = Species(i)%Init(iInit)%YieldSEE * NbrOfPhotons / Species(i)%MacroParticleFactor &
                               + Species(i)%Init(iInit)%NINT_Correction 
                  NbrOfParticle = NINT(NbrOfPhotons)
                  Species(i)%Init(iInit)%NINT_Correction = NbrOfPhotons - REAL(NbrOfParticle)
                ELSE
                  ! Photo-ionization in the volume
                  ! Calculation of the number of photons (using actual number and applying the weighting factor on the number of reactions)
                  NbrOfPhotons = Species(i)%Init(iInit)%EffectiveIntensityFactor * NbrOfPhotons
                  ! Calculation of the number of photons depending on the cylinder height (ratio of actual to virtual cylinder height, which
                  ! is spanned by the disk and the length given by c*dt)
                  NbrOfPhotons = NbrOfPhotons * Species(i)%Init(iInit)%CylinderHeightIC / (c*dt)
                  ! Calculation of the number of electron resulting from the chemical reactions in the photoionization region
                  CALL CalcPhotoIonizationNumber(NbrOfPhotons,NbrOfReactions)
                  NbrOfReactions = NbrOfReactions + Species(i)%Init(iInit)%NINT_Correction
                  NbrOfParticle = NINT(NbrOfReactions)
                  Species(i)%Init(iInit)%NINT_Correction = NbrOfReactions - REAL(NbrOfParticle)
                END IF
              ELSE
                NbrOfParticle = 0
              END IF ! MOD(MERGE(Time-T0/2., Time, Time.GE.T0/2.), Period).GT.T0
            ELSE
              NbrOfParticle = 0
            END IF ! Time.LE.Species(i)%Init(iInit)%tActive
          END ASSOCIATE
        CASE DEFAULT
          NbrOfParticle = 0
        END SELECT

       CALL SetParticlePosition(i,iInit,NbrOfParticle)
        ! Pairing of "electrons" with the background species and performing the reaction
        IF(TRIM(Species(i)%Init(iInit)%SpaceIC).EQ.'photon_cylinder') THEN
          CALL BGGas_PhotoIonization(i,iInit,NbrOfParticle)
          CYCLE
        END IF

       CALL SetParticleVelocity(i,iInit,NbrOfParticle)
       CALL SetParticleChargeAndMass(i,NbrOfParticle)
       IF (usevMPF) CALL SetParticleMPF(i,NbrOfParticle)
       IF (VarTimeStep%UseVariableTimeStep) CALL SetParticleTimeStep(NbrOfParticle)
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
  END DO  ! iInit 
END DO  ! i=1,nSpecies

END SUBROUTINE ParticleInserting


SUBROUTINE AdaptiveBCAnalyze(initSampling_opt)
!===================================================================================================================================
! Sampling of variables (part-density, velocity and energy) for Adaptive BC elements
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
USE MOD_DSMC_Analyze           ,ONLY: CalcTVib,CalcTVibPoly,CalcTelec
USE MOD_DSMC_Vars              ,ONLY: PartStateIntEn, DSMC, CollisMode, SpecDSMC, useDSMC, RadialWeighting
USE MOD_Mesh_Vars              ,ONLY: nElems, offsetElem
USE MOD_Part_Tools             ,ONLY: GetParticleWeight
USE MOD_Particle_Boundary_Vars ,ONLY: PorousBCSampIter, PorousBCMacroVal
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemInfo_Shared, SideInfo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemVolume_Shared
USE MOD_Particle_Vars          ,ONLY: PartState, PDM, PartSpecies, Species, nSpecies, PEM, Adaptive_MacroVal, AdaptiveWeightFac
USE MOD_Particle_Vars          ,ONLY: usevMPF
USE MOD_Timedisc_Vars          ,ONLY: iter
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
LOGICAL                         :: initSampling, isBCElem
#if USE_LOADBALANCE
REAL                            :: tLBStart
#endif /*USE_LOADBALANCE*/
INTEGER                         :: nlocSides, GlobalSideID, GlobalElemID, iLocSide
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
    GlobalElemID = PEM%GlobalElemID(i)
    ElemID = GlobalElemID - offsetElem
    IF ((ElemID.LT.1).OR.(ElemID.GT.nElems)) CYCLE
    ! not a BC element
    ! IF (ElemToBCSides(ELEM_NBR_BCSIDES,ElemID).EQ.-1) CYCLE
    ! ======================================
    ! ElemToBCSides is only built for RefMapping, eigenes Array beschreiben statt untere schleife jedes mal abzulaufen
    isBCElem = .FALSE.
    nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,GlobalElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID)
    DO iLocSide = 1,nlocSides
      GlobalSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID) + iLocSide
      IF (SideInfo_Shared(SIDE_BCID,GlobalSideID).GT.0) isBCElem = .TRUE.
    END DO
    IF(.NOT.isBCElem) CYCLE
    ! ======================================
#if USE_LOADBALANCE
    nPartsPerBCElem(ElemID) = nPartsPerBCElem(ElemID) + 1
#endif /*USE_LOADBALANCE*/
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

DO AdaptiveElemID = 1,nElems
  GlobalElemID = AdaptiveElemID + offsetElem
  ! not a BC element
  ! IF (ElemToBCSides(ELEM_NBR_BCSIDES,AdaptiveElemID).EQ.-1) CYCLE
  ! ======================================
  ! ElemToBCSides is only built for RefMapping, eigenes Array beschreiben statt untere schleife jedes mal abzulaufen?
  isBCElem = .FALSE.
  nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,GlobalElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID)
  DO iLocSide = 1,nlocSides
    GlobalSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID) + iLocSide
    IF (SideInfo_Shared(SIDE_BCID,GlobalSideID).GT.0) isBCElem = .TRUE.
  END DO
  IF(.NOT.isBCElem) CYCLE
  ! ======================================
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
                                                        /ElemVolume_Shared(GlobalElemID)
            ELSE
              Adaptive_MacroVal(7,AdaptiveElemID,iSpec)=PorousBCMacroVal(7,AdaptiveElemID,iSpec)/SamplingIteration   &
                                                        /ElemVolume_Shared(GlobalElemID) * Species(iSpec)%MacroParticleFactor
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
            + RelaxationFactor*Source(11,AdaptiveElemID,iSpec) /ElemVolume_Shared(GlobalElemID)
          ! Pressure with relaxation factor
          Adaptive_MacroVal(12,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(12,AdaptiveElemID,iSpec) &
          +RelaxationFactor*Source(11,AdaptiveElemID,iSpec)/ElemVolume_Shared(GlobalElemID)*TTrans_TempFac
        ELSE
          ! compute density
          Adaptive_MacroVal(7,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(7,AdaptiveElemID,iSpec) &
            + RelaxationFactor*Source(11,AdaptiveElemID,iSpec)/ElemVolume_Shared(GlobalElemID)*Species(iSpec)%MacroParticleFactor
          ! Pressure with relaxation factor
          Adaptive_MacroVal(12,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(12,AdaptiveElemID,iSpec) &
                                                      + RelaxationFactor*Source(11,AdaptiveElemID,iSpec)  &
                                                      / ElemVolume_Shared(GlobalElemID)*Species(iSpec)%MacroParticleFactor*TTrans_TempFac
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
