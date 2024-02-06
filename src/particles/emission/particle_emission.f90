!==================================================================================================================================
! Copyright (c) 2010 - 2019 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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
PUBLIC :: ParticleInserting
!===================================================================================================================================
CONTAINS

SUBROUTINE ParticleInserting()
!===================================================================================================================================
!> Particle emission at every iteration (ParticleEmissionType.GT.0)
!===================================================================================================================================
! Modules
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPIInitGroup
#endif /*USE_MPI*/
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: c,PI
USE MOD_Timedisc_Vars          ,ONLY: dt,time
USE MOD_Timedisc_Vars          ,ONLY: RKdtFrac,RKdtFracTotal
USE MOD_Particle_Vars
USE MOD_PIC_Vars
USE MOD_part_tools             ,ONLY: UpdateNextFreePosition, GetNextFreePosition, IncreaseMaxParticleNumber
USE MOD_DSMC_Vars              ,ONLY: useDSMC, CollisMode, SpecDSMC
USE MOD_part_emission_tools    ,ONLY: DSMC_SetInternalEnr_LauxVFD
USE MOD_DSMC_PolyAtomicModel   ,ONLY: DSMC_SetInternalEnr_Poly
USE MOD_Particle_Analyze_Vars  ,ONLY: CalcPartBalance,nPartIn,PartEkinIn
USE MOD_Particle_Analyze_Tools ,ONLY: CalcEkinPart
USE MOD_part_emission_tools    ,ONLY: SetParticleChargeAndMass,SetParticleMPF,SamplePoissonDistri,SetParticleTimeStep,CalcNbrOfPhotons
USE MOD_part_emission_tools    ,ONLY: CountNeutralizationParticles
USE MOD_PICDepo_Tools          ,ONLY: DepositPhotonSEEHoles
USE MOD_part_pos_and_velo      ,ONLY: SetParticlePosition,SetParticleVelocity
USE MOD_DSMC_BGGas             ,ONLY: BGGas_PhotoIonization
USE MOD_DSMC_ChemReact         ,ONLY: CalcPhotoIonizationNumber
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_ReadInTools            ,ONLY: PrintOption
#if defined(MEASURE_MPI_WAIT)
USE MOD_Particle_MPI_Vars      ,ONLY: MPIW8TimePart,MPIW8CountPart
#endif /*defined(MEASURE_MPI_WAIT)*/
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: SEE,CalcElectronSEE
USE MOD_Particle_Photoionization  ,ONLY: PhotoIonization_RayTracing_Volume, PhotoIonization_RayTracing_SEE
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Local variable declaration
INTEGER                          :: i, iPart, PositionNbr, iInit, IntSample
INTEGER                          :: NbrOfParticle,iSEEBC
INTEGER(KIND=8)                  :: inserted_Particle_iter,inserted_Particle_time
INTEGER(KIND=8)                  :: inserted_Particle_diff
REAL                             :: PartIns, RandVal1
REAL                             :: RiseFactor, RiseTime,NbrOfPhotons
REAL                             :: dtVar, TimeVar
#if USE_MPI
INTEGER                          :: InitGroup
#endif
REAL                             :: NbrOfReactions,NbrOfParticlesReal,MPF
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=8)                  :: CounterStart,CounterEnd
REAL(KIND=8)                     :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================

!--- Ray tracing based volume photo-ionization
CALL PhotoIonization_RayTracing_Volume()

!--- Ray tracing based secondary electron emission
CALL PhotoIonization_RayTracing_SEE()

!---  Emission at time step
DO i=1,nSpecies
  ! Species-specific time step
  IF(VarTimeStep%UseSpeciesSpecific) THEN
    dtVar = dt * Species(i)%TimeStepFactor
    TimeVar = Time * Species(i)%TimeStepFactor
  ELSE
    dtVar = dt
    TimeVar = Time
  END IF
  DO iInit = 1, Species(i)%NumberOfInits
    ! Reset the number of particles per species AND init region
    NbrOfParticle = 0
    ! Only use inits defined for emission (every time step)
    IF (Species(i)%Init(iInit)%ParticleEmissionType.LE.0) CYCLE

    SELECT CASE(Species(i)%Init(iInit)%ParticleEmissionType)
      CASE(1) ! Emission Type: Particles per !!!!!SECOND!!!!!!!! (not per ns)
        IF (.NOT.DoPoissonRounding .AND. .NOT.DoTimeDepInflow) THEN
          PartIns=Species(i)%Init(iInit)%ParticleNumber * dtVar*RKdtFrac  ! emitted particles during time-slab
          inserted_Particle_iter = INT(PartIns,8)                                     ! number of particles to be inserted
          PartIns=Species(i)%Init(iInit)%ParticleNumber * (TimeVar + dtVar*RKdtFracTotal) ! total number of emitted particle over
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
          Species(i)%Init(iInit)%InsertedParticleSurplus = ABS(MIN(inserted_Particle_iter + inserted_Particle_diff,0_8))
          NbrOfParticle = MAX(INT(inserted_Particle_iter + inserted_Particle_diff,4),0)
          !-- if maxwell velo dist and less than 5 parts: skip (to ensure maxwell dist)
          IF (TRIM(Species(i)%Init(iInit)%velocityDistribution).EQ.'maxwell') THEN
            IF (NbrOfParticle.LT.5) NbrOfParticle=0
          END IF
        ELSE IF (DoPoissonRounding .AND. .NOT.DoTimeDepInflow) THEN
          ! linear rise of inflow
          RiseTime=Species(i)%Init(iInit)%InflowRiseTime
          IF(RiseTime.GT.0.)THEN
            IF(time-DelayTime.LT.RiseTime)THEN
              RiseFactor=(time-DelayTime)/RiseTime
            ELSE
              RiseFactor=1.
            END IF
          ELSE
            RiseFactor=1.
          EnD IF
          PartIns=Species(i)%Init(iInit)%ParticleNumber * dtVar*RKdtFrac * RiseFactor  ! emitted particles during time-slab
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
            IF(time-DelayTime.LT.RiseTime)THEN
              RiseFactor=(time-DelayTime)/RiseTime
            ELSE
              RiseFactor=1.
            END IF
          ELSE
            RiseFactor=1.
          EnD IF
          ! emitted particles during time-slab
          PartIns=Species(i)%Init(iInit)%ParticleNumber * dtVar*RKdtFrac * RiseFactor &
                  + Species(i)%Init(iInit)%InsertedParticleMisMatch
          CALL RANDOM_NUMBER(RandVal1)
          NbrOfParticle = INT(PartIns + RandVal1)
        END IF
#if USE_MPI
        ! communicate number of particles with all procs in the same init group
        InitGroup=Species(i)%Init(iInit)%InitCOMM
        IF(PartMPIInitGroup(InitGroup)%COMM.NE.MPI_COMM_NULL) THEN
          ! only procs which are part of group take part in the communication
            !NbrOfParticle based on RandVals!
          CALL MPI_BCAST(NbrOfParticle, 1, MPI_INTEGER,0,PartMPIInitGroup(InitGroup)%COMM,IERROR)
        ELSE
          NbrOfParticle=0
        END IF
        !CALL MPI_BCAST(NbrOfParticle, 1, MPI_INTEGER,0,PartMPI%COMM,IERROR) !NbrOfParticle based on RandVals!
#endif
        Species(i)%Init(iInit)%InsertedParticle = Species(i)%Init(iInit)%InsertedParticle + INT(NbrOfParticle,8)
      CASE(2)    ! Emission Type: Particles per Iteration
        IF (RKdtFracTotal .EQ. 1.) THEN !insert in last stage only, so that no reconstruction is nec. and number/iter matches
          NbrOfParticle = INT(Species(i)%Init(iInit)%ParticleNumber)
        ELSE
          NbrOfParticle = 0
        END IF
      CASE(7) ! SEE based on photon impact and photo-ionization in the volume
        ASSOCIATE( tShift => Species(i)%Init(iInit)%tShift )
          ! Check if all pulses have terminated
          IF(time.LE.Species(i)%Init(iInit)%tActive)THEN
            ! Check if pulse is currently active of in between two pulses (in the latter case, do nothing)
            IF(MOD(time, Species(i)%Init(iInit)%Period).LE.2.0*tShift)THEN
              ! Calculate the number of currently active photons (both surface SEE and volumetric emission)
              CALL CalcNbrOfPhotons(i, iInit, NbrOfPhotons)

              ! Check if only particles in the first quadrant are to be inserted that is spanned by the vectors
              ! x=BaseVector1IC and y=BaseVector2IC in the interval x,y in [0,R] and reduce the number of photon accordingly
              IF(Species(i)%Init(iInit)%FirstQuadrantOnly) NbrOfPhotons = NbrOfPhotons / 4.0

              ! Select surface SEE or volumetric emission
              SELECT CASE(TRIM(Species(i)%Init(iInit)%SpaceIC))
              CASE('photon_SEE_disc','photon_SEE_honeycomb','photon_SEE_rectangle')
                ! SEE based on photon impact
                IF(usevMPF)THEN
                  MPF = Species(i)%Init(iInit)%MacroParticleFactor ! Use emission-specific MPF
                ELSE
                  MPF = Species(i)%MacroParticleFactor ! Use species MPF
                END IF ! usevMPF
                NbrOfPhotons = Species(i)%Init(iInit)%YieldSEE * NbrOfPhotons / MPF + Species(i)%Init(iInit)%NINT_Correction
                NbrOfParticle = NINT(NbrOfPhotons)
                Species(i)%Init(iInit)%NINT_Correction = NbrOfPhotons - REAL(NbrOfParticle)
              CASE DEFAULT
                ! Photo-ionization in the volume
                ! Calculation of the number of photons (using actual number and applying the weighting factor on the number of reactions)
                NbrOfPhotons = Species(i)%Init(iInit)%EffectiveIntensityFactor * NbrOfPhotons
                ! Calculation of the number of photons depending on the cylinder height (ratio of actual to virtual cylinder height, which
                ! is spanned by the disk and the length given by c*dt)
                IF(TRIM(Species(i)%Init(iInit)%SpaceIC).EQ.'photon_rectangle')THEN
                  ! Rectangular area -> cuboid: Equally distributed over c*dt
                  NbrOfPhotons = NbrOfPhotons * Species(i)%Init(iInit)%CuboidHeightIC / (c*dt)
                ELSE
                  ! Cylinder and honeycomb: Equally distributed over c*dt
                  NbrOfPhotons = NbrOfPhotons * Species(i)%Init(iInit)%CylinderHeightIC / (c*dt)
                END IF
                ! Calculation of the number of electron resulting from the chemical reactions in the photoionization region
                CALL CalcPhotoIonizationNumber(i,NbrOfPhotons,NbrOfReactions)
                NbrOfReactions = NbrOfReactions + Species(i)%Init(iInit)%NINT_Correction
                NbrOfParticle = NINT(NbrOfReactions)
                Species(i)%Init(iInit)%NINT_Correction = NbrOfReactions - REAL(NbrOfParticle)
              END SELECT
            ELSE
              NbrOfParticle = 0
            END IF ! MOD(time, Period) .LE. 2x tShift
          ELSE
            NbrOfParticle = 0
          END IF ! time.LE.Species(i)%Init(iInit)%tActive
        END ASSOCIATE
      CASE(8) ! SpaceIC='2D_landmark','2D_landmark_copy'
              ! Ionization profile from T. Charoy, 2D axial-azimuthal particle-in-cell benchmark
              ! for low-temperature partially magnetized plasmas (2019)
       ASSOCIATE( x2 => 1.0e-2  ,& ! m
                  x1 => 0.25e-2 ,& ! m
                  Ly => 1.28e-2 ,& ! m
                  S0 => 5.23e23 )  ! m^-3 s^-1
         NbrOfParticlesReal = Ly*dt*S0*2.0*(x2-x1)/PI ! yields 1.60E+08 m^-1 (i.e. per metre in z-direction)
         NbrOfParticlesReal = NbrOfParticlesReal*(GEO%zmaxglob-GEO%zminglob)/Species(i)%MacroParticleFactor &
                              + Species(i)%Init(iInit)%NINT_Correction
         NbrOfParticle = NINT(NbrOfParticlesReal)
         Species(i)%Init(iInit)%NINT_Correction = NbrOfParticlesReal - REAL(NbrOfParticle)
         IF(.NOT.ALLOCATED(PartPosLandmark))THEN
           NbrOfParticleLandmarkMax = MAX(10,NINT(NbrOfParticlesReal*1.2)) ! Add 20% safety
           CALL PrintOption('Landmark volume emission allocation size NbrOfParticleLandmarkMax' , 'CALCUL.' , IntOpt=NbrOfParticleLandmarkMax)
           ALLOCATE(PartPosLandmark(1:3,1:NbrOfParticleLandmarkMax))
           PartPosLandmark=HUGE(1.)
           IF(NbrOfParticleLandmarkMax.LE.0)THEN
             IPWRITE(UNIT_StdOut,*) "NbrOfParticleLandmarkMax =", NbrOfParticleLandmarkMax
             CALL abort(__STAMP__,'NbrOfParticleLandmarkMax.LE.0')
           END IF
         ELSE
           IF(NbrOfParticleLandmarkMax.LT.NbrOfParticle)THEN
             IPWRITE(UNIT_StdOut,*) "NbrOfParticleLandmarkMax,NbrOfParticle =", NbrOfParticleLandmarkMax,NbrOfParticle
             CALL abort(__STAMP__&
             ,'NbrOfParticleLandmarkMax.LT.NbrOfParticle is not allowed! Allocate PartPosLandmark to the appropriate size.')
           END IF ! NbrOfParticleLandmarkMax.LE.NbrOfParticle
         END IF ! .NOT.ALLOCATED()
       END ASSOCIATE
     CASE(9) ! '2D_landmark_neutralization',
             ! '2D_Liu2010_neutralization'      ,'3D_Liu2010_neutralization'
             ! '2D_Liu2010_neutralization_Szabo','3D_Liu2010_neutralization_Szabo'
#if USE_MPI
       ! Communicate number of particles with all procs in the same init group
       InitGroup=Species(i)%Init(iInit)%InitCOMM
       NeutralizationBalanceGlobal=0 ! always nullify
       IF(PartMPIInitGroup(InitGroup)%COMM.NE.MPI_COMM_NULL) THEN
         ! Loop over all elements and count the ion surplus per element if element-local emission is used
         IF(nNeutralizationElems.GT.0) CALL CountNeutralizationParticles()
         ! Only processors which are part of group take part in the communication
         CALL MPI_ALLREDUCE(NeutralizationBalance,NeutralizationBalanceGlobal,1,MPI_INTEGER,MPI_SUM,PartMPIInitGroup(InitGroup)%COMM,IERROR)
       ELSE
         NeutralizationBalanceGlobal=0
       END IF
#else
       NeutralizationBalanceGlobal = NeutralizationBalance
#endif
       ! Insert electrons only when the number is greater than zero
       IF(NeutralizationBalanceGlobal.GT.0)THEN
         ! Insert only when positive
         NbrOfParticle = NeutralizationBalanceGlobal
         ! Reset the counter but only when not using element-local emission, nullify later is this case (in SetParticlePosition)
         IF(nNeutralizationElems.EQ.-1) NeutralizationBalance = 0
       ELSE
         NbrOfParticle = 0
       END IF ! NeutralizationBalance.GT.0

      CASE DEFAULT
        NbrOfParticle = 0
    END SELECT

    ! Create particles by setting their position in space and checking if a host cell can be found
    ! Warning: this routine returns the emitted number of particles for each processor and changes the value of NbrOfParticle here
    CALL SetParticlePosition(i,iInit,NbrOfParticle)

    ! Pairing of "electrons" with the background species and performing the reaction
    SELECT CASE(TRIM(Species(i)%Init(iInit)%SpaceIC))
    CASE('photon_cylinder','photon_honeycomb','photon_rectangle')
      CALL BGGas_PhotoIonization(i,iInit,NbrOfParticle)
      CYCLE
    END SELECT
    ! Check if photon SEE electric current is to be measured
    IF(StringBeginsWith(Species(i)%Init(iInit)%SpaceIC,'photon_SEE').AND.CalcElectronSEE.AND.(NbrOfParticle.GT.0))THEN
      ! Note that the negative value of the charge -q is used below
      iSEEBC = SEE%BCIDToSEEBCID(Species(i)%Init(iInit)%PartBCIndex)
      SEE%RealElectronOut(iSEEBC) = SEE%RealElectronOut(iSEEBC) - MPF*NbrOfParticle*Species(i)%ChargeIC
    END IF

    CALL SetParticleVelocity(i,iInit,NbrOfParticle)
    CALL SetParticleChargeAndMass(i,NbrOfParticle)
    IF (usevMPF) CALL SetParticleMPF(i,iInit,NbrOfParticle)
    ! Check if photon SEE is to be considered in the surface charging
    IF(StringBeginsWith(Species(i)%Init(iInit)%SpaceIC,'photon_SEE').AND.(NbrOfParticle.GT.0)) CALL DepositPhotonSEEHoles(&
        Species(i)%Init(iInit)%PartBCIndex,NbrOfParticle)
    IF (UseVarTimeStep) CALL SetParticleTimeStep(NbrOfParticle)
    ! define molecule stuff
    IF (useDSMC.AND.(CollisMode.GT.1)) THEN
      DO iPart = 1, NbrOfParticle
        PositionNbr = GetNextFreePosition(iPart)
        IF (PositionNbr.NE.0) THEN
          IF (SpecDSMC(i)%PolyatomicMol) THEN
            CALL DSMC_SetInternalEnr_Poly(i,iInit,PositionNbr,1)
          ELSE
            CALL DSMC_SetInternalEnr_LauxVFD(i,iInit,PositionNbr,1)
          END IF
        END IF
      END DO
    END IF
    ! Compute number of input particles and energy
    IF(CalcPartBalance.AND.(NbrOfParticle.GT.0)) THEN
      nPartIn(i)=nPartIn(i) + NbrOfparticle
      DO iPart=1,NbrOfParticle
        PositionNbr = GetNextFreePosition(iPart)
        PartEkinIn(i) = PartEkinIn(i) + CalcEkinPart(PositionNbr)
      END DO ! iPart
    END IF ! CalcPartBalance
    ! Update the current next free position and increase the particle vector length
    IF(NbrOfParticle.GT.0) THEN
      PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + NbrOfParticle
      PDM%ParticleVecLength = MAX(PDM%ParticleVecLength,GetNextFreePosition(0))
    END IF
#ifdef CODE_ANALYZE
    IF(PDM%ParticleVecLength.GT.PDM%maxParticleNumber) CALL Abort(__STAMP__,'PDM%ParticleVeclength exceeds PDM%maxParticleNumber, Difference:',IntInfoOpt=PDM%ParticleVeclength-PDM%maxParticleNumber)
    DO iPart=PDM%ParticleVecLength+1,PDM%maxParticleNumber
      IF (PDM%ParticleInside(iPart)) THEN
        IPWRITE(*,*) iPart,PDM%ParticleVecLength,PDM%maxParticleNumber
        CALL Abort(__STAMP__,'Particle outside PDM%ParticleVeclength',IntInfoOpt=iPart)
      END IF
    END DO
#endif

    ! Complete check if all particles were emitted successfully
#if USE_MPI
    InitGroup=Species(i)%Init(iInit)%InitCOMM
    IF (PartMPIInitGroup(InitGroup)%COMM.NE.MPI_COMM_NULL .AND. Species(i)%Init(iInit)%sumOfRequestedParticles.GT.0) THEN
#if defined(MEASURE_MPI_WAIT)
      CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/
      CALL MPI_WAIT(PartMPIInitGroup(InitGroup)%Request, MPI_STATUS_IGNORE, iError)
#if defined(MEASURE_MPI_WAIT)
      CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
      MPIW8TimePart(5)  = MPIW8TimePart(5) + REAL(CounterEnd-CounterStart,8)/Rate
      MPIW8CountPart(5) = MPIW8CountPart(5) + 1_8
#endif /*defined(MEASURE_MPI_WAIT)*/

      IF(PartMPIInitGroup(InitGroup)%MPIRoot) THEN
#endif
        ! add number of matching error to particle emission to fit
        ! number of added particles
        Species(i)%Init(iInit)%InsertedParticleMisMatch = Species(i)%Init(iInit)%sumOfRequestedParticles - Species(i)%Init(iInit)%sumOfMatchedParticles
        IF (Species(i)%Init(iInit)%sumOfRequestedParticles .GT. Species(i)%Init(iInit)%sumOfMatchedParticles) THEN
          WRITE(UNIT_StdOut,'(A)')      'WARNING in ParticleEmission_parallel:'
          WRITE(UNIT_StdOut,'(A,I0)')   'Fraction Nbr: '  , i
          WRITE(UNIT_StdOut,'(A,I0,A)') 'matched only '   , Species(i)%Init(iInit)%sumOfMatchedParticles  , ' particles'
          WRITE(UNIT_StdOut,'(A,I0,A)') 'when '           , Species(i)%Init(iInit)%sumOfRequestedParticles, ' particles were required!'
        ELSE IF (Species(i)%Init(iInit)%sumOfRequestedParticles .LT. Species(i)%Init(iInit)%sumOfMatchedParticles) THEN
          WRITE(UNIT_StdOut,'(A)')      'ERROR in ParticleEmission_parallel:'
          WRITE(UNIT_StdOut,'(A,I0)')   'Fraction Nbr: '  , i
          WRITE(UNIT_StdOut,'(A,I0,A)') 'matched '        , Species(i)%Init(iInit)%sumOfMatchedParticles  , ' particles'
          WRITE(UNIT_StdOut,'(A,I0,A)') 'when '           , Species(i)%Init(iInit)%sumOfRequestedParticles, ' particles were required!'
        ! ELSE IF (nbrOfParticle .EQ. Species(i)%Init(iInit)%sumOfMatchedParticles) THEN
        !  WRITE(UNIT_stdOut,'(A,I0)')   'Fraction Nbr: '  , FractNbr
        !  WRITE(UNIT_stdOut,'(A,I0,A)') 'ParticleEmission_parallel: matched all (',NbrOfParticle,') particles!'
        END IF
#if USE_MPI
      END IF ! PartMPIInitGroup(InitGroup)%MPIRoot
    END IF
#endif
  END DO  ! iInit
END DO  ! i=1,nSpecies

END SUBROUTINE ParticleInserting

END MODULE MOD_part_emission
