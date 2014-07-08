MODULE MOD_part_emission
!===================================================================================================================================
! module for particle emmision
!===================================================================================================================================
! implicit varible handling
   IMPLICIT NONE                                                                                   !
   PRIVATE                                                                                         !
!----------------------------------------------------------------------------------------------------------------------------------

   PUBLIC         :: InitializeParticleEmission, ParticleInserting &
                , SetParticleChargeAndMass, SetParticleVelocity, SetParticleMPF                              
!----------------------------------------------------------------------------------------------------------------------------------
                                                                                                  
CONTAINS                                                                                           
                                                                                                   
SUBROUTINE InitializeParticleEmission()
!===================================================================================================================================
! Initialize particles / Insert initial particles
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars,  ONLY : Species,nSpecies,PDM,PEM,Time, usevMPF
USE MOD_PIC_Vars,       ONLY : PIC
USE MOD_part_tools,     ONLY : UpdateNextFreePosition
USE MOD_Restart_Vars,   ONLY : DoRestart 
USE MOD_ReadInTools
USE MOD_DSMC_Vars,      ONLY : useDSMC, DSMC
USE MOD_part_pressure,  ONLY : ParticleInsideCheck
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: i, NbrOfParticle,iInit,iPart,PositionNbr
  CHARACTER(32)         :: hilf , hilf2
  !CHARACTER(40)                          :: SpaceIC                          ! specifying Keyword for Particle Space condition
  !CHARACTER(30)                          :: velocityDistribution             ! specifying keyword for velocity distribution
  !INTEGER(8)                             :: initialParticleNumber            ! Number of Particles at time 0.0
  !REAL                                   :: RadiusIC
  !REAL                                   :: Radius2IC
  !REAL                                   :: RadiusICGyro
  !REAL                                   :: NormalIC(3)                      ! Normal / Orientation of circle
  !REAL                                   :: BasePointIC(3)                   ! base point for IC cuboid and IC sphere
  !REAL                                   :: BaseVector1IC(3)                 ! first base vector for IC cuboid
  !REAL                                   :: BaseVector2IC(3)                 ! second base vector for IC cuboid
  !REAL                                   :: CuboidHeightIC                   ! third measure of cuboid
  !REAL                                   :: CylinderHeightIC                 ! third measure of cylinder
  !REAL                                   :: VeloIC                           ! velocity for inital Data
  !REAL                                   :: VeloVecIC(3)                     ! normalized velocity vector
  !REAL                                   :: MWTemperatureIC                  ! Temperature for Maxwell Distribution
  !REAL                                   :: Alpha                            ! WaveNumber for sin-deviation initiation.
  !REAL                                   :: Amplitude                        ! Amplitude for sin-deviation initiation.
  !REAL                                   :: WaveNumber                       ! WaveNumber for sin-deviation initiation.
  !INTEGER(8)                             :: maxParticleNumberX               ! Maximum Number of all Particles in x direction
  !INTEGER(8)                             :: maxParticleNumberY               ! Maximum Number of all Particles in y direction
  !INTEGER(8)                             :: maxParticleNumberZ  
  INTEGER                                :: nPartInside!, nElem_check
  REAL                                   :: EInside,TempInside!,ConstantPressure,ConstPressureRelaxFac,UseForInit,PartDensity
  LOGICAL                                :: EmType6
!==================================================================================================

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
IF (.NOT.DoRestart) THEN
!   CALL Deposition()
!   IF (MESH%t.GE.PIC%DelayTime) PIC%ParticleTreatmentMethod='standard'
  ! for the case of particle insertion per time, the inserted particle number for the current time must
  ! be updated. Otherwise, at the first timestep after restart, these particles will be inserted again
!  DO i=1,nSpecies
!    Species(i)%InsertedParticle = INT(Species(i)%ParticleEmission * Time)
!  END DO
!ELSE
  DO i = 1,nSpecies
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
        IF (Species(i)%Init(iInit)%ParticleEmissionType .EQ. 4) THEN
          CALL ParticleInsertingCellPressure(i,iInit,NbrofParticle)
          CALL SetParticleVelocity(i,iInit,NbrOfParticle)
        ELSE !emission type 6 (constant pressure outflow)
          CALL ParticleInsertingPressureOut(i,iInit,NbrofParticle)
        END IF
        CALL SetParticleChargeAndMass(i,NbrOfParticle)
        IF (usevMPF) CALL SetParticleMPF(i,NbrOfParticle)
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
        NbrOfParticle = Species(i)%Init(iInit)%initialParticleNumber
#ifdef MPI
        CALL SetParticlePosition(i,iInit,NbrOfParticle,1)
        CALL SetParticlePosition(i,iInit,NbrOfParticle,2)
#else
        CALL SetParticlePosition(i,iInit,NbrOfParticle)
#endif /*MPI*/
        CALL SetParticleVelocity(i,iInit,NbrOfParticle)
        CALL SetParticleChargeAndMass(i,NbrOfParticle)
        IF (usevMPF) CALL SetParticleMPF(i,NbrOfParticle)
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
          CALL ParticleInsideCheck(i, iInit, nPartInside, TempInside, EInside)
          IF (Species(i)%Init(iInit)%ParticleEmission .GT. nPartInside) THEN
            NbrOfParticle = Species(i)%Init(iInit)%ParticleEmission - nPartInside
            WRITE(*,*) 'Emission PartNum (Spec ',i,')', NbrOfParticle
#ifdef MPI
            CALL SetParticlePosition(i,iInit,NbrOfParticle,1)
            CALL SetParticlePosition(i,iInit,NbrOfParticle,2)
#else
            CALL SetParticlePosition(i,iInit,NbrOfParticle)
#endif
            CALL SetParticleVelocity(i,iInit,NbrOfParticle)
            CALL SetParticleChargeAndMass(i,NbrOfParticle)
            IF (usevMPF) CALL SetParticleMPF(i,NbrOfParticle)
            !IF (useDSMC) CALL SetParticleIntEnergy(i,NbrOfParticle)
            PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfParticle
            CALL UpdateNextFreePosition()
          END IF
        END IF
      END IF ! not Emissiontype 4
    END DO !inits
  END DO ! species
END IF ! not restart

!--- set last element to current element (needed when ParticlePush is not executed, e.g. "delay")
DO i = 1,PDM%ParticleVecLength
  PEM%lastElement(i) = PEM%Element(i)
END DO

END SUBROUTINE InitializeParticleEmission

#ifdef MPI
SUBROUTINE ParticleInserting(mode_opt)                                                             
#else
SUBROUTINE ParticleInserting()                                                                     
#endif
!===================================================================================================================================
! Particle Inserting
!===================================================================================================================================
! Modules
  USE MOD_Timedisc_Vars         , ONLY : dt, iter
  USE MOD_Particle_Vars       ! , ONLY : time
  USE MOD_PIC_Vars
  USE MOD_part_tools             ,ONLY : UpdateNextFreePosition  
  USE MOD_DSMC_Vars              ,ONLY : useDSMC, CollisMode, DSMC                                   
  USE MOD_DSMC_Init              ,ONLY : DSMC_SetInternalEnr_LauxVFD
#if (PP_TimeDiscMethod==1000)
  USE MOD_LD_Init                ,ONLY : CalcDegreeOfFreedom
  USE MOD_LD_Vars
#endif
  USE MOD_part_pressure          ,ONLY: ParticlePressure, ParticlePressureRem
!===================================================================================================================================
! implicit variable handling
!===================================================================================================================================
   IMPLICIT NONE                                                                                   
!----------------------------------------------------------------------------------------------------------------------------------
! argument list declaration                                                                        
#ifdef MPI
   INTEGER, OPTIONAL                :: mode_opt
#endif
! Local variable declaration                                                                       
   INTEGER                          :: i , iPart, PositionNbr, iInit                                                          
   INTEGER                , SAVE    :: NbrOfParticle=0                                             
   INTEGER(KIND=8)                          :: inserted_Particle_iter,inserted_Particle_time               
   INTEGER(KIND=8)                  :: inserted_Particle_diff  
#ifdef MPI
   INTEGER                          :: mode,iProc,iTag                                             
#endif
!----------------------------------------------------------------------------------------------------------------------------------
!!! VORSICHT: FUNKTIONIERT SO MOMENTAN NUR MIT 1 SPEZIES!!!!
! --- für mehr als eine Spezies gibt es bei der Benutzung des
!     mode_opt Flags Probleme mit den non-blocking communications.
!     Es könnte dann passieren, dass Nachrichten falsch zugeordnet werden.
!     Sicherheitshalber sollte man kein mode_opt Argument bei mehrern
!     Spezies übergeben.
#ifdef MPI
         IF (PRESENT(mode_opt)) THEN
           mode=mode_opt
         ELSE
           mode=0
         END IF
#endif
   !---  Emission at time step (initial emission see pic_init.f90, InitParticles)
   DO i=1,nSpecies
     DO iInit = Species(i)%StartnumberOfInits, Species(i)%NumberOfInits
       IF (((Species(i)%Init(iInit)%ParticleEmissionType .NE. 4).AND.(Species(i)%Init(iInit)%ParticleEmissionType .NE. 6)) .AND. &
            (Species(i)%Init(iInit)%UseForEmission)) THEN ! no constant density in cell type, + to be used for init
#ifdef MPI
         IF (mode.NE.2) THEN
#endif
           SELECT CASE(Species(i)%Init(iInit)%ParticleEmissionType)
           CASE(1) ! Emission Type: Particles per !!!!!SECOND!!!!!!!! (not per ns)
             inserted_Particle_iter = INT(Species(i)%Init(iInit)%ParticleEmission * dt,8)
             inserted_Particle_time = INT(Species(i)%Init(iInit)%ParticleEmission * (Time + dt),8)
             inserted_Particle_diff = INT(inserted_Particle_time - Species(i)%Init(iInit)%InsertedParticle &
                  - inserted_Particle_iter,8)
             NbrOfParticle = inserted_Particle_iter + inserted_Particle_diff
             ! if maxwell velo dist and less than 5 parts: skip (to ensure maxwell dist)
             IF (Species(i)%Init(iInit)%velocityDistribution.EQ.'maxwell') THEN
               IF (NbrOfParticle.LT.5) NbrOfParticle=0
             END IF
             Species(i)%Init(iInit)%InsertedParticle = Species(i)%Init(iInit)%InsertedParticle + NbrOfParticle
           CASE(2)    ! Emission Type: Particles per Iteration
             NbrOfParticle = INT(Species(i)%Init(iInit)%ParticleEmission)
           CASE(3)
             CALL ParticlePressure (i, iInit, NbrOfParticle)
             ! if maxwell velo dist and less than 5 parts: skip (to ensure maxwell dist)
             IF (Species(i)%Init(iInit)%velocityDistribution.EQ.'maxwell') THEN
               IF (NbrOfParticle.LT.5) NbrOfParticle=0
             END IF
           CASE(5) ! removal of all parts in pressure area and re-insertion
             CALL ParticlePressureRem (i, iInit, NbrOfParticle)
           CASE DEFAULT
             NbrOfParticle = 0
           END SELECT
#ifdef MPI
           CALL SetParticlePosition(i,iInit,NbrOfParticle,1)
         END IF
         IF (mode.NE.1) THEN
           CALL SetParticlePosition(i,iInit,NbrOfParticle,2)
#else
           CALL SetParticlePosition(i,iInit,NbrOfParticle)
#endif
           CALL SetParticleVelocity(i,iInit,NbrOfParticle)
           CALL SetParticleChargeAndMass(i,NbrOfParticle)
           IF (usevMPF) CALL SetParticleMPF(i,NbrOfParticle)
           ! define molecule stuff
           IF (useDSMC.AND.(CollisMode.NE.1)) THEN
             iPart = 1
             DO WHILE (iPart .le. NbrOfParticle)
               PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
               IF (PositionNbr .ne. 0) THEN
                 CALL DSMC_SetInternalEnr_LauxVFD(i, iInit, PositionNbr)
               END IF
               iPart = iPart + 1
             END DO
           END IF
#if (PP_TimeDiscMethod==1000)
           iPart = 1
           DO WHILE (iPart .le. NbrOfParticle)
             PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
             IF (PositionNbr .ne. 0) THEN
               PartStateBulkValues(PositionNbr,1) = Species(i)%Init(iInit)%VeloVecIC(1) * Species(i)%Init(iInit)%VeloIC
               PartStateBulkValues(PositionNbr,2) = Species(i)%Init(iInit)%VeloVecIC(2) * Species(i)%Init(iInit)%VeloIC
               PartStateBulkValues(PositionNbr,3) = Species(i)%Init(iInit)%VeloVecIC(3) * Species(i)%Init(iInit)%VeloIC
               PartStateBulkValues(PositionNbr,4) = Species(i)%Init(iInit)%MWTemperatureIC
               PartStateBulkValues(PositionNbr,5) = CalcDegreeOfFreedom(PositionNbr)
             END IF
             iPart = iPart + 1
           END DO
#endif
           ! instead of UpdateNextfreePosition we update the
           ! particleVecLength only.
           PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + NbrOfParticle
           PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfParticle
           !CALL UpdateNextFreePosition()
#ifdef MPI
         END IF
#endif
       ELSE IF (Species(i)%Init(iInit)%UseForEmission) THEN ! Constant Pressure in Cell Emission (type 4 or 6)
         IF (Species(i)%Init(iInit)%ParticleEmissionType .EQ. 4) THEN
           CALL ParticleInsertingCellPressure(i,iInit,NbrofParticle)
           CALL SetParticleVelocity(i,iInit,NbrOfParticle)
         ELSE !emission type 6 (constant pressure outflow)
           CALL ParticleInsertingPressureOut(i,iInit,NbrofParticle)
         END IF
         CALL SetParticleChargeAndMass(i,NbrOfParticle)
         IF (usevMPF) CALL SetParticleMPF(i,NbrOfParticle)
         ! define molecule stuff
         IF (useDSMC.AND.(CollisMode.NE.1)) THEN
           iPart = 1
           DO WHILE (iPart .le. NbrOfParticle)
             PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
             IF (PositionNbr .ne. 0) THEN
               CALL DSMC_SetInternalEnr_LauxVFD(i, iInit, PositionNbr)
             END IF
             iPart = iPart + 1
           END DO
         END IF
#if (PP_TimeDiscMethod==1000)
         iPart = 1
         DO WHILE (iPart .le. NbrOfParticle)
           PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
           IF (PositionNbr .ne. 0) THEN
             PartStateBulkValues(PositionNbr,1) = Species(i)%Init(iInit)%VeloVecIC(1) * Species(i)%Init(iInit)%VeloIC
             PartStateBulkValues(PositionNbr,2) = Species(i)%Init(iInit)%VeloVecIC(2) * Species(i)%Init(iInit)%VeloIC
             PartStateBulkValues(PositionNbr,3) = Species(i)%Init(iInit)%VeloVecIC(3) * Species(i)%Init(iInit)%VeloIC
             PartStateBulkValues(PositionNbr,4) = Species(i)%Init(iInit)%MWTemperatureIC
             PartStateBulkValues(PositionNbr,5) = CalcDegreeOfFreedom(PositionNbr)
           END IF
           iPart = iPart + 1
         END DO
#endif
         ! instead of UpdateNextfreePosition we update the
         ! particleVecLength only.
         PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + NbrOfParticle
         PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfParticle
         !CALL UpdateNextFreePosition()
       END IF
     END DO
   END DO
 RETURN
END SUBROUTINE ParticleInserting
                                                                                                   
#ifdef MPI
SUBROUTINE SetParticlePosition(FractNbr,iInit,NbrOfParticle,mode)
!===================================================================================================================================
  USE MOD_part_MPI_Vars,      ONLY : PMPIVAR,PMPIInsert
!===================================================================================================================================
#else
SUBROUTINE SetParticlePosition(FractNbr,iInit,NbrOfParticle)                                             
#endif
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
  ! modules
  USE MOD_Particle_Vars
  USE MOD_PIC_Vars
  USE MOD_Timedisc_Vars,         ONLY : dt, iter, IterDisplayStep
  USE MOD_BoundaryTools,         ONLY : SingleParticleToExactElement                                  
  USE MOD_PICInterpolation,      ONLY : InterpolateCurvedExternalField
  USE MOD_PICInterpolation_vars, ONLY : CurvedExternalField,useCurvedExternalField
  USE MOD_Equation_vars,         ONLY : c_inv
  USE MOD_LD,                    ONLY : LD_SetParticlePosition
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
   IMPLICIT NONE                                                                                   
!----------------------------------------------------------------------------------------------------------------------------------
#ifdef MPI
    INCLUDE 'mpif.h'                                                                               
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! argument list declaration                                                                        
   INTEGER                          :: FractNbr,iInit,NbrOfParticle                                      
#ifdef MPI
   INTEGER                          :: mode                                                        
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! local variables                                                                                  
#ifdef MPI
   INTEGER                          :: iProc, CellX, CellY, CellZ                                  
   INTEGER                          :: msg_status(1:MPI_STATUS_SIZE)                               
   LOGICAL                          :: InsideMyBGM                                                 
#endif
   REAL,POINTER                     :: particle_positions(:)=>NULL(), gen_particle_positions(:)=>NULL()!
   INTEGER                          :: IERROR, allocStat                                           
   INTEGER                          :: i,j,k,ParticleIndexNbr                                      
   INTEGER                          :: mySumOfMatchedParticles, sumOfMatchedParticles              
   INTEGER                          :: nChunks, chunkSize                                          
   REAL                             :: start_pos(3), end_pos(3),lineVector(3),VectorGap(3)         
   REAL                             :: RandVal(3), Particle_pos(3),lineVector2(3)                  
   REAL                             :: n(3) , radius_vec(3)                                        
   REAL                             :: II(3,3),JJ(3,3),NN(3,3)                                     
   REAL                             :: RandVal1, RandVal2, DeltaTime                              
   REAL                             :: x_0, y_0, z_0, radius, argumentTheta                        
   REAL                             :: Velo(3),q,m,Pos(3),v0(3),s0(3)                              
   REAL                             :: myRealTestValue        ! TestVariable                       
   REAL                             :: rgyrate, Bintpol
   INTEGER                          :: myRealKind             ! 4 oder 8 ?                         
   REAL                             :: x_step, y_step, z_step,  x_pos , y_pos                      
   REAL                             :: xlen, ylen, zlen                                            
   INTEGER                          :: iPart                                                       
   CHARACTER(50)                    :: debugFileName                                               
   REAL, PARAMETER                  :: PI=3.14159265358979323846_8
   REAL,ALLOCATABLE                 :: particle_positions_Temp(:)                                
!----------------------------------------------------------------------------------------------------------------------------------
   INTENT(IN)                       :: FractNbr, iInit                                                    
   INTENT(INOUT)                    :: NbrOfParticle                                               
!----------------------------------------------------------------------------------------------------------------------------------


   IF (nbrOfParticle .LE. 0) RETURN

   DeltaTime = dt

   nChunks = 1                   ! Standard: Nicht-MPI
   sumOfMatchedParticles = 0
   mySumOfMatchedParticles = 0

   chunkSize = nbrOfParticle
   ! process myRank=0 generates the complete list of random positions for all emitted particles
#ifdef MPI
   IF ( (nbrOfParticle.GT.10*PMPIVAR%nProcs                                      ) .AND. &
        (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).NE.'circle_equidistant'                 ) .AND. &
        (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).NE.'sin_deviation'                      ) .AND. &
        (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).NE.'cuboid_with_equidistant_distribution') .AND. &
        (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).NE.'line_with_equidistant_distribution' )) THEN
      nChunks = PMPIVAR%nProcs
   ELSE
      nChunks = 1
   END IF
IF (mode.EQ.1) THEN
   chunkSize = INT(nbrOfParticle/nChunks)
   IF (PMPIVAR%iProc .EQ. 0) chunkSize = chunkSize + ( nbrOfParticle - (nChunks*chunkSize) )
   IF (PMPIVAR%iProc .EQ. 0 .OR. nChunks.GT.1) THEN
#endif
      ALLOCATE( particle_positions(1:chunkSize*3), STAT=allocStat )
      IF (allocStat .NE. 0) THEN
         WRITE(*,*)'ERROR in SetParticlePosition: cannot allocate particle_positions!'
         STOP
      END IF

!      WRITE(*,*)'generating',chunkSize,'*',nChunks,' particles...'

      SELECT CASE(TRIM(Species(FractNbr)%Init(iInit)%SpaceIC))
      CASE ('point')
         Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC
         DO i=1,chunkSize
            particle_positions(i*3-2) = Particle_pos(1)
            particle_positions(i*3-1) = Particle_pos(2)
            particle_positions(i*3  ) = Particle_pos(3)
         END DO
      CASE ('line_with_equidistant_distribution')
        IF(NbrOfParticle.EQ.1)THEN
           Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + 0.5 * Species(FractNbr)%Init(iInit)%BaseVector1IC
        ELSE
           VectorGap = Species(FractNbr)%Init(iInit)%BaseVector1IC/(NbrOfParticle-1)
           DO i=1,chunkSize
              Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + (i-1)*VectorGap
              particle_positions(i*3-2) = Particle_pos(1)
              particle_positions(i*3-1) = Particle_pos(2)
              particle_positions(i*3  ) = Particle_pos(3)
           END DO
        END IF
      CASE ('line')
        DO i=1,chunkSize
           CALL RANDOM_NUMBER(RandVal1)
           Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + Species(FractNbr)%Init(iInit)%BaseVector1IC*RandVal1
           particle_positions(i*3-2) = Particle_pos(1)
           particle_positions(i*3-1) = Particle_pos(2)
           particle_positions(i*3  ) = Particle_pos(3)
        END DO
      CASE('disc')
        IF (Species(FractNbr)%Init(iInit)%NormalIC(3).NE.0) THEN
           lineVector(1) = 1.0
           lineVector(2) = 1.0
           lineVector(3) = -(Species(FractNbr)%Init(iInit)%NormalIC(1)+Species(FractNbr)%Init(iInit)%NormalIC(2))/ &
                             Species(FractNbr)%Init(iInit)%NormalIC(3)
        ELSE
           IF (Species(FractNbr)%Init(iInit)%NormalIC(2).NE.0) THEN
              lineVector(1) = 1.0
              lineVector(3) = 1.0
              lineVector(2) = -(Species(FractNbr)%Init(iInit)%NormalIC(1)+Species(FractNbr)%Init(iInit)%NormalIC(3))/ &
                                Species(FractNbr)%Init(iInit)%NormalIC(2)
           ELSE
              IF (Species(FractNbr)%Init(iInit)%NormalIC(1).NE.0) THEN
                 lineVector(2) = 1.0
                 lineVector(3) = 1.0
                 lineVector(1) = -(Species(FractNbr)%Init(iInit)%NormalIC(2)+Species(FractNbr)%Init(iInit)%NormalIC(3))/ &
                                   Species(FractNbr)%Init(iInit)%NormalIC(1)
              ELSE
                 WRITE(*,*) 'Error in SetParticlePosition, NormalIC Vektor darf nicht Nullvektor sein'
                 STOP
              END IF
           END IF
        END IF
        
        lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * &
             lineVector(2) + lineVector(3) * lineVector(3))
        
        lineVector2(1) = Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(3) - &
             Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(2)
        lineVector2(2) = Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(1) - &
             Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(3)
        lineVector2(3) = Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(2) - &
             Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(1)
        
        lineVector2 = lineVector2 / SQRT(lineVector2(1) * lineVector2(1) + lineVector2(2) * &
             lineVector2(2) + lineVector2(3) * lineVector2(3))

        DO i=1,chunkSize
           radius = Species(FractNbr)%Init(iInit)%RadiusIC + 1
           DO WHILE(radius.GT.Species(FractNbr)%Init(iInit)%RadiusIC)
              CALL RANDOM_NUMBER(RandVal)
              RandVal = RandVal * 2. - 1.
              Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + Species(FractNbr)%Init(iInit)%RadiusIC * &
                       (RandVal(1) * lineVector + RandVal(2) *lineVector2)

              radius = SQRT( (Particle_pos(1)-Species(FractNbr)%Init(iInit)%BasePointIC(1)) * &
                             (Particle_pos(1)-Species(FractNbr)%Init(iInit)%BasePointIC(1)) + &
                             (Particle_pos(2)-Species(FractNbr)%Init(iInit)%BasePointIC(2)) * &
                             (Particle_pos(2)-Species(FractNbr)%Init(iInit)%BasePointIC(2)) + &
                             (Particle_pos(3)-Species(FractNbr)%Init(iInit)%BasePointIC(3)) * &
                             (Particle_pos(3)-Species(FractNbr)%Init(iInit)%BasePointIC(3)) )
           END DO
           particle_positions(i*3-2) = Particle_pos(1)
           particle_positions(i*3-1) = Particle_pos(2)
           particle_positions(i*3  ) = Particle_pos(3)
        END DO
      CASE('circle')
        IF (Species(FractNbr)%Init(iInit)%NormalIC(3).NE.0) THEN
           lineVector(1) = 1.0
           lineVector(2) = 1.0
           lineVector(3) = -(Species(FractNbr)%Init(iInit)%NormalIC(1)+Species(FractNbr)%Init(iInit)%NormalIC(2))/ &
                             Species(FractNbr)%Init(iInit)%NormalIC(3)
        ELSE
           IF (Species(FractNbr)%Init(iInit)%NormalIC(2).NE.0) THEN
              lineVector(1) = 1.0
              lineVector(3) = 1.0
              lineVector(2) = -(Species(FractNbr)%Init(iInit)%NormalIC(1)+Species(FractNbr)%Init(iInit)%NormalIC(3))/ &
                                Species(FractNbr)%Init(iInit)%NormalIC(2)
           ELSE
              IF (Species(FractNbr)%Init(iInit)%NormalIC(1).NE.0) THEN
                 lineVector(2) = 1.0
                 lineVector(3) = 1.0
                 lineVector(1) = -(Species(FractNbr)%Init(iInit)%NormalIC(2)+Species(FractNbr)%Init(iInit)%NormalIC(3))/ &
                                   Species(FractNbr)%Init(iInit)%NormalIC(1)
              ELSE
                 WRITE(*,*) 'Error in SetParticlePosition, NormalIC Vektor darf nicht Nullvektor sein'
                 STOP
              END IF
           END IF
        END IF
        
        lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * &
             lineVector(2) + lineVector(3) * lineVector(3))
        
        lineVector2(1) = Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(3) - &
             Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(2)
        lineVector2(2) = Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(1) - &
             Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(3)
        lineVector2(3) = Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(2) - &
             Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(1)
        
        lineVector2 = lineVector2 / SQRT(lineVector2(1) * lineVector2(1) + lineVector2(2) * &
             lineVector2(2) + lineVector2(3) * lineVector2(3))

        radius = Species(FractNbr)%Init(iInit)%RadiusIC
        DO i=1,chunkSize
           CALL RANDOM_NUMBER(RandVal1)
           argumentTheta = 2.*pi*RandVal1
           Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC +        &
                          linevector * cos(argumentTheta) * radius +  &
                          linevector2 * sin(argumentTheta) * radius
           particle_positions(i*3-2) = Particle_pos(1)
           particle_positions(i*3-1) = Particle_pos(2)
           particle_positions(i*3  ) = Particle_pos(3)
        END DO
      CASE('gyrotron_circle')
        IF (Species(FractNbr)%Init(iInit)%NormalIC(3).NE.0) THEN
           lineVector(1) = 1.0
           lineVector(2) = 1.0
           lineVector(3) = -(Species(FractNbr)%Init(iInit)%NormalIC(1)+Species(FractNbr)%Init(iInit)%NormalIC(2))/ &
                             Species(FractNbr)%Init(iInit)%NormalIC(3)
        ELSE
           IF (Species(FractNbr)%Init(iInit)%NormalIC(2).NE.0) THEN
              lineVector(1) = 1.0
              lineVector(3) = 1.0
              lineVector(2) = -(Species(FractNbr)%Init(iInit)%NormalIC(1)+Species(FractNbr)%Init(iInit)%NormalIC(3))/ &
                                Species(FractNbr)%Init(iInit)%NormalIC(2)
           ELSE
              IF (Species(FractNbr)%Init(iInit)%NormalIC(1).NE.0) THEN
                 lineVector(2) = 1.0
                 lineVector(3) = 1.0
                 lineVector(1) = -(Species(FractNbr)%Init(iInit)%NormalIC(2)+Species(FractNbr)%Init(iInit)%NormalIC(3))/ &
                                   Species(FractNbr)%Init(iInit)%NormalIC(1)
              ELSE
                 WRITE(*,*) 'Error in SetParticlePosition, NormalIC Vektor darf nicht Nullvektor sein'
                 STOP
              END IF
           END IF
        END IF
        
        lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * &
             lineVector(2) + lineVector(3) * lineVector(3))
        
        lineVector2(1) = Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(3) - &
             Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(2)
        lineVector2(2) = Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(1) - &
             Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(3)
        lineVector2(3) = Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(2) - &
             Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(1)
        
        lineVector2 = lineVector2 / SQRT(lineVector2(1) * lineVector2(1) + lineVector2(2) * &
             lineVector2(2) + lineVector2(3) * lineVector2(3))

        radius = Species(FractNbr)%Init(iInit)%RadiusIC
        DO i=1,chunkSize
           CALL RANDOM_NUMBER(RandVal1)
           argumentTheta = 2.*pi*RandVal1
           Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC +        &
                          linevector * cos(argumentTheta) * radius +  &
                          linevector2 * sin(argumentTheta) * radius
           ! Change position of particle on the small gyro circle
           ! take normal vecotr of the circle
           n(1:3) = Species(FractNbr)%Init(iInit)%NormalIC(1:3)
           ! generate radius vector (later it will be multiplied by the length of the
           ! gyro circles. For now we just need the vector)
           radius_vec(1) = Particle_pos(1) - Species(FractNbr)%Init(iInit)%BasePointIC(1)
           radius_vec(2) = Particle_pos(2) - Species(FractNbr)%Init(iInit)%BasePointIC(2)
           radius_vec(3) = Particle_pos(3) - Species(FractNbr)%Init(iInit)%BasePointIC(3)
           !rotate radius vector with random angle
           CALL RANDOM_NUMBER(RandVal1)
           argumentTheta=2.*pi*RandVal1
           JJ(1,1:3) = (/   0.,-n(3), n(2)/)
           JJ(2,1:3) = (/ n(3),   0.,-n(1)/)
           JJ(3,1:3) = (/-n(2), n(1),   0./)
           II(1,1:3) = (/1.,0.,0./)
           II(2,1:3) = (/0.,1.,0./)
           II(3,1:3) = (/0.,0.,1./)
           forall(j=1:3) NN(:,j) = n(:)*n(j)

!          1. determine the z-position in order to get the interpolated curved B-field
           CALL RANDOM_NUMBER(RandVal1)
           IF (NbrOfParticle.EQ.Species(FractNbr)%Init(iInit)%initialParticleNumber) THEN
             particle_positions(i*3  ) = Species(FractNbr)%Init(iInit)%BasePointIC(3) &
                                             + RandVal1 * Species(FractNbr)%Init(iInit)%CuboidHeightIC
           ELSE
             particle_positions(i*3  ) = Species(FractNbr)%Init(iInit)%BasePointIC(3) &
                                             + RandVal1 * dt & 
                                             * Species(FractNbr)%Init(iInit)%VeloIC/Species(FractNbr)%Init(iInit)%alpha 
           END IF

!          2. calculate curved B-field at z-position in order to determin size of gyroradius
           IF (useCurvedExternalField) THEN
              Bintpol = InterpolateCurvedExternalField(particle_positions(i*3))
              rgyrate = 1./ SQRT ( 1 - (Species(FractNbr)%Init(iInit)%VeloIC**2 * (1 + 1./Species(FractNbr)%Init(iInit)%alpha**2)) &
                                  * c_inv * c_inv ) * Species(FractNbr)%MassIC * Species(FractNbr)%Init(iInit)%VeloIC / &
                        ( Bintpol * abs( Species(FractNbr)%ChargeIC) )
           ELSE
             rgyrate =  Species(FractNbr)%Init(iInit)%RadiusICGyro
           END IF

           radius_vec = MATMUL( NN+cos(argumentTheta)*(II-NN)+sin(argumentTheta)*JJ , radius_vec ) 
           radius_vec(1:3) = radius_vec(1:3) / SQRT(radius_vec(1)**2+radius_vec(2)**2+radius_vec(3)**2) &
                         * rgyrate !Species(1)%RadiusICGyro
           ! Set new particles position:
           particle_positions(i*3-2) = Particle_pos(1) + radius_vec(1)
           particle_positions(i*3-1) = Particle_pos(2) + radius_vec(2)
           !particle_positions(i*3  )=0.
        END DO
      CASE('circle_equidistant')
        IF (Species(FractNbr)%Init(iInit)%NormalIC(3).NE.0) THEN
           lineVector(1) = 1.0
           lineVector(2) = 1.0
           lineVector(3) = -(Species(FractNbr)%Init(iInit)%NormalIC(1)+Species(FractNbr)%Init(iInit)%NormalIC(2))/ &
                             Species(FractNbr)%Init(iInit)%NormalIC(3)
        ELSE
           IF (Species(FractNbr)%Init(iInit)%NormalIC(2).NE.0) THEN
              lineVector(1) = 1.0
              lineVector(3) = 1.0
              lineVector(2) = -(Species(FractNbr)%Init(iInit)%NormalIC(1)+Species(FractNbr)%Init(iInit)%NormalIC(3))/ &
                                Species(FractNbr)%Init(iInit)%NormalIC(2)
           ELSE
              IF (Species(FractNbr)%Init(iInit)%NormalIC(1).NE.0) THEN
                 lineVector(2) = 1.0
                 lineVector(3) = 1.0
                 lineVector(1) = -(Species(FractNbr)%Init(iInit)%NormalIC(2)+Species(FractNbr)%Init(iInit)%NormalIC(3))/ &
                                   Species(FractNbr)%Init(iInit)%NormalIC(1)
              ELSE
                 WRITE(*,*) 'Error in SetParticlePosition, NormalIC Vektor darf nicht Nullvektor sein'
                 STOP
              END IF
           END IF
        END IF
                   
        lineVector2(1) = Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(3) - &
             Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(2)
        lineVector2(2) = Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(1) - &
             Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(3)
        lineVector2(3) = Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(2) - &
             Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(1)

        lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * &
             lineVector(2) + lineVector(3) * lineVector(3))

        lineVector2 = lineVector2 / SQRT(lineVector2(1) * lineVector2(1) + lineVector2(2) * &
             lineVector2(2) + lineVector2(3) * lineVector2(3))

        radius = Species(FractNbr)%Init(iInit)%RadiusIC
        DO i=1,chunkSize
           argumentTheta = 2.*pi*i/chunkSize
           Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC +        &
                          linevector * cos(argumentTheta) * radius +  &
                          linevector2 * sin(argumentTheta) * radius
           particle_positions(i*3-2) = Particle_pos(1)
           particle_positions(i*3-1) = Particle_pos(2)
           particle_positions(i*3  ) = Particle_pos(3)
        END DO
      CASE('cuboid')
        DO i=1,chunkSize          
           CALL RANDOM_NUMBER(RandVal)
           Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + Species(FractNbr)%Init(iInit)%BaseVector1IC * RandVal(1)
           Particle_pos = Particle_pos + Species(FractNbr)%Init(iInit)%BaseVector2IC * RandVal(2)
           lineVector(1) = Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * Species(FractNbr)%Init(iInit)%BaseVector2IC(3) - &
                Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * Species(FractNbr)%Init(iInit)%BaseVector2IC(2)
           lineVector(2) = Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * Species(FractNbr)%Init(iInit)%BaseVector2IC(1) - &
                Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * Species(FractNbr)%Init(iInit)%BaseVector2IC(3)
           lineVector(3) = Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * Species(FractNbr)%Init(iInit)%BaseVector2IC(2) - &
                Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * Species(FractNbr)%Init(iInit)%BaseVector2IC(1)
           IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
           ELSE
              lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + &
                   lineVector(3) * lineVector(3))
           END IF          
#if (PP_TimeDiscMethod==201)
           Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%CuboidHeightIC * dt / dt_maxwell * RandVal(3) 
           !scaling due to variable time step
#else
           Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%CuboidHeightIC * RandVal(3) 
#endif
           particle_positions(i*3-2) = Particle_pos(1)
           particle_positions(i*3-1) = Particle_pos(2)
           particle_positions(i*3  ) = Particle_pos(3)
        END DO
      CASE('cylinder')
        DO i=1,chunkSize
           radius = Species(FractNbr)%Init(iInit)%RadiusIC + 1
           DO WHILE((radius.GT.Species(FractNbr)%Init(iInit)%RadiusIC) .OR.(radius.LT.Species(FractNbr)%Init(iInit)%Radius2IC))
              CALL RANDOM_NUMBER(RandVal)
              Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC+Species(FractNbr)%Init(iInit)%BaseVector1IC*(RandVal(1)*2-1)
              Particle_pos = Particle_pos + Species(FractNbr)%Init(iInit)%BaseVector2IC * (RandVal(2)*2-1)

              radius = SQRT( (Particle_pos(1)-Species(FractNbr)%Init(iInit)%BasePointIC(1)) * &
                             (Particle_pos(1)-Species(FractNbr)%Init(iInit)%BasePointIC(1)) + &
                             (Particle_pos(2)-Species(FractNbr)%Init(iInit)%BasePointIC(2)) * &
                             (Particle_pos(2)-Species(FractNbr)%Init(iInit)%BasePointIC(2)) + &
                             (Particle_pos(3)-Species(FractNbr)%Init(iInit)%BasePointIC(3)) * &
                             (Particle_pos(3)-Species(FractNbr)%Init(iInit)%BasePointIC(3)) )
           END DO
           lineVector(1) = Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * Species(FractNbr)%Init(iInit)%BaseVector2IC(3) - &
                Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * Species(FractNbr)%Init(iInit)%BaseVector2IC(2)
           lineVector(2) = Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * Species(FractNbr)%Init(iInit)%BaseVector2IC(1) - &
                Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * Species(FractNbr)%Init(iInit)%BaseVector2IC(3)
           lineVector(3) = Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * Species(FractNbr)%Init(iInit)%BaseVector2IC(2) - &
                Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * Species(FractNbr)%Init(iInit)%BaseVector2IC(1)
           IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
           ELSE
              lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + &
                   lineVector(3) * lineVector(3))
           END IF          
#if (PP_TimeDiscMethod==201)
           Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%CylinderHeightIC * dt / dt_maxwell * RandVal(3) 
           !scaling due to variable time step           
#else
           Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%CylinderHeightIC * RandVal(3)
#endif
           particle_positions(i*3-2) = Particle_pos(1)
           particle_positions(i*3-1) = Particle_pos(2)
           particle_positions(i*3  ) = Particle_pos(3)
        END DO
      CASE('LD_insert')
        CALL LD_SetParticlePosition(chunkSize,particle_positions_Temp,FractNbr)
        DEALLOCATE( particle_positions, STAT=allocStat )
        IF (allocStat .NE. 0) THEN
          WRITE(*,*)'ERROR in ParticleEmission_parallel: cannot deallocate particle_positions!'
          STOP
        END IF
        ALLOCATE(particle_positions(3*chunkSize))
        particle_positions(1:3*chunkSize) = particle_positions_Temp(1:3*chunkSize)
        DEALLOCATE( particle_positions_Temp, STAT=allocStat )
        IF (allocStat .NE. 0) THEN
          WRITE(*,*)'ERROR in ParticleEmission_parallel: cannot deallocate particle_positions!'
          STOP
        END IF
      CASE('cuboid_equal')
#ifdef MPI
         IF (PMPIVAR%nProcs .GT. 1) THEN
           WRITE(*,*)'WARNING in SetParticlePosition:'
           WRITE(*,*)'cannot fully handle Particle Initial Condition \"cuboid equal\"'
           WRITE(*,*)'in parallel mode (with more than one CPU)!'
           WRITE(*,*)'USE WITH CARE!!!'
!           STOP
         END IF
         j=0
         mySumOfMatchedParticles = 0
         DO i=1,PDM%ParticleVecLength
            j=j+1
            ParticleIndexNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
            IF (ParticleIndexNbr .ne. 0) THEN
               PartState(ParticleIndexNbr,1:3) = PartState(j,1:3)
               PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
               CALL SingleParticleToExactElement(ParticleIndexNbr)
               IF (PDM%ParticleInside(ParticleIndexNbr)) THEN
                  mySumOfMatchedParticles = mySumOfMatchedParticles + 1
               ELSE
                  PDM%ParticleInside(ParticleIndexNbr) = .FALSE.
               END IF
            ELSE
               WRITE(*,*)'ERROR in SetParticlePosition:'
               WRITE(*,*)'ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?'
               STOP
            END IF
         END DO
         CALL MPI_ALLREDUCE(mySumOfMatchedParticles, sumOfMatchedParticles, 1, MPI_INTEGER, MPI_SUM, PMPIVAR%COMM, IERROR)
         nbrOfParticle = NbrOfParticle - sumOfMatchedParticles
         IF (nbrOfParticle .NE. 0) THEN
            WRITE(*,*)'ERROR in ParticleEmission_parallel:'
            WRITE(*,'(A,I7,A)')'matched ', sumOfMatchedParticles, ' particles'
            WRITE(*,'(A,I7,A)')'when ', NbrOfParticle+sumOfMatchedParticles, ' particles were required!'
            STOP
!         ELSE IF (nbrOfParticle .EQ. 0) THEN
!            WRITE(*,*)'ParticleEmission_parallel: matched all particles!'
         END IF
         NbrOfParticle = mySumOfMatchedParticles
         DEALLOCATE( particle_positions, STAT=allocStat )
         IF (allocStat .NE. 0) THEN
            WRITE(*,*)'ERROR in ParticleEmission_parallel: cannot deallocate particle_positions!'
            STOP
         END IF
         RETURN
#else
         DO i=1,chunkSize
            ParticleIndexNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
            particle_positions(i*3-2 : i*3) = PartState(ParticleIndexNbr-Species(FractNbr)%Init(iInit)%initialParticleNumber,1:3)
         END DO
#endif
      CASE ('cuboid_with_equidistant_distribution') 
         IF(Species(FractNbr)%Init(iInit)%initialParticleNumber.NE. &
              (Species(FractNbr)%Init(iInit)%maxParticleNumberX * Species(FractNbr)%Init(iInit)%maxParticleNumberY &
              * Species(FractNbr)%Init(iInit)%maxParticleNumberZ)) THEN
           WRITE(*,*) 'ERROR: Number of particles in init / emission region',iInit
           WRITE(*,*) 'of species ',FractNbr,' does not match number of particles in each direction!'
           STOP
         END IF
         xlen = SQRT(Species(FractNbr)%Init(iInit)%BaseVector1IC(1)**2 &
              + Species(FractNbr)%Init(iInit)%BaseVector1IC(2)**2 &
              + Species(FractNbr)%Init(iInit)%BaseVector1IC(3)**2 )
         ylen = SQRT(Species(FractNbr)%Init(iInit)%BaseVector2IC(1)**2 &
              + Species(FractNbr)%Init(iInit)%BaseVector2IC(2)**2 &
              + Species(FractNbr)%Init(iInit)%BaseVector2IC(3)**2 )
         zlen = ABS(Species(FractNbr)%Init(iInit)%CuboidHeightIC)

         ! make sure the vectors correspond to x,y,z-dir
         IF ((xlen.NE.Species(FractNbr)%Init(iInit)%BaseVector1IC(1)).OR. &
            (ylen.NE.Species(FractNbr)%Init(iInit)%BaseVector2IC(2)).OR. &
            (zlen.NE.Species(FractNbr)%Init(iInit)%CuboidHeightIC)) THEN
            WRITE(*,*) 'Basevectors1IC,-2IC and CuboidHeightIC have to be in x,y,z-direction, respectively for emission condition'
            STOP
          END IF
         x_step = xlen/Species(FractNbr)%Init(iInit)%maxParticleNumberX
         y_step = ylen/Species(FractNbr)%Init(iInit)%maxParticleNumberY
         z_step = zlen/Species(FractNbr)%Init(iInit)%maxParticleNumberZ
         iPart = 1
         DO i=1,Species(FractNbr)%Init(iInit)%maxParticleNumberX
           x_pos = (i-0.5) * x_step + Species(FractNbr)%Init(iInit)%BasePointIC(1)
           DO j=1,Species(FractNbr)%Init(iInit)%maxParticleNumberY
             y_pos =  Species(FractNbr)%Init(iInit)%BasePointIC(2) + (j-0.5) * y_step
             DO k=1,Species(FractNbr)%Init(iInit)%maxParticleNumberZ
               particle_positions(iPart*3-2) = x_pos
               particle_positions(iPart*3-1) = y_pos
               particle_positions(iPart*3  ) = Species(FractNbr)%Init(iInit)%BasePointIC(3) &
                    + (k-0.5) * z_step
               iPart = iPart + 1
             END DO
           END DO
         END DO
      CASE('sin_deviation')
         IF(Species(FractNbr)%Init(iInit)%initialParticleNumber.NE. &
              (Species(FractNbr)%Init(iInit)%maxParticleNumberX * Species(FractNbr)%Init(iInit)%maxParticleNumberY &
              * Species(FractNbr)%Init(iInit)%maxParticleNumberZ)) THEN
           WRITE(*,*) 'ERROR: Number of particles in init / emission region',iInit
           WRITE(*,*) 'of species ',FractNbr,' does not match number of particles in each direction!'
           STOP
         END IF
         xlen = abs(GEO%xmaxglob  - GEO%xminglob)  
         ylen = abs(GEO%ymaxglob  - GEO%yminglob)
         zlen = abs(GEO%zmaxglob  - GEO%zminglob)
         x_step = xlen/Species(FractNbr)%Init(iInit)%maxParticleNumberX
         y_step = ylen/Species(FractNbr)%Init(iInit)%maxParticleNumberY
         z_step = zlen/Species(FractNbr)%Init(iInit)%maxParticleNumberZ
         iPart = 1
         DO i=1,Species(FractNbr)%Init(iInit)%maxParticleNumberX
            x_pos = (i * x_step - x_step*0.5)
            x_pos = GEO%xminglob + x_pos + Species(FractNbr)%Init(iInit)%Amplitude &
                    * sin(Species(FractNbr)%Init(iInit)%WaveNumber * x_pos)
            DO j=1,Species(FractNbr)%Init(iInit)%maxParticleNumberY
              y_pos =  GEO%yminglob + j * y_step - y_step * 0.5
              DO k=1,Species(FractNbr)%Init(iInit)%maxParticleNumberZ
                particle_positions(iPart*3-2) = x_pos                                
                particle_positions(iPart*3-1) = y_pos
                particle_positions(iPart*3  ) = GEO%zminglob &
                                          + k * z_step - z_step * 0.5
                iPart = iPart + 1
              END DO
            END DO
         END DO
      END SELECT
#ifdef MPI
   ELSE !no mpi root, nchunks=1
     chunkSize=0
   END IF
  
   IF(nChunks.GT.1) THEN
     ALLOCATE( PMPIInsert%nbrOfSendParticlesEmission(0:PMPIVAR%nProcs-1), STAT=allocStat )
     ALLOCATE( PMPIInsert%nbrOfRecvParticlesEmission(0:PMPIVAR%nProcs-1), STAT=allocStat )
     ALLOCATE( PMPIInsert%send_request(0:PMPIVAR%nProcs-1,1:2), STAT=allocStat )
     ALLOCATE( PMPIInsert%recv_request(0:PMPIVAR%nProcs-1,1:2), STAT=allocStat )
     ALLOCATE( PMPIInsert%send_message(0:PMPIVAR%nProcs-1), STAT=allocStat )
     PMPIInsert%nbrOfSendParticlesEmission(:)=0
     DO i=1,chunkSize
       CellX = INT((particle_positions(i*3-2)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
       CellY = INT((particle_positions(i*3-1)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
       CellZ = INT((particle_positions(i*3  )-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
       InsideMyBGM=.TRUE.
       IF ((CellX.GT.GEO%FIBGMimax).OR.(CellX.LT.GEO%FIBGMimin) .OR. &
           (CellY.GT.GEO%FIBGMkmax).OR.(CellY.LT.GEO%FIBGMkmin) .OR. &
           (CellZ.GT.GEO%FIBGMlmax).OR.(CellZ.LT.GEO%FIBGMlmin)) THEN
         InsideMyBGM=.FALSE.
       END If
       IF (InsideMyBGM) THEN
         IF (.NOT.ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) InsideMyBGM=.FALSE.
       END IF
       IF (InsideMyBGM) THEN
         DO j=2,GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(1)+1
           iProc=GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(j)
           PMPIInsert%nbrOfSendParticlesEmission(iProc)=PMPIInsert%nbrOfSendParticlesEmission(iProc)+1
         END DO
         PMPIInsert%nbrOfSendParticlesEmission(PMPIVAR%iProc)=PMPIInsert%nbrOfSendParticlesEmission(PMPIVAR%iProc)+1
       ELSE
         DO iProc=0,PMPIVAR%nProcs-1
  !         IF (iProc.EQ.PMPIVAR%iProc) CYCLE
           PMPIInsert%nbrOfSendParticlesEmission(iProc)=PMPIInsert%nbrOfSendParticlesEmission(iProc)+1
         END DO
       END IF
     END DO   
   ELSE
      IF(PMPIVAR%iProc .EQ. 0) THEN
        ALLOCATE( PMPIInsert%send_message(0:0), STAT=allocStat )
        ALLOCATE( PMPIInsert%send_message(0)%content(1:3*chunkSize), STAT=allocStat )
        PMPIInsert%send_message(0)%content(:)=particle_positions(:)
        DEALLOCATE(particle_positions, STAT=allocStat)
      END IF
   END IF
   IF (nChunks.GT.1) THEN
      DO iProc=0,PMPIVAR%nProcs-1
        !--- MPI_ISEND lengths of lists of particles leaving local mesh
        CALL MPI_ISEND(PMPIInsert%nbrOfSendParticlesEmission(iProc), 1, MPI_INTEGER, iProc, 1011+FractNbr, PMPIVAR%COMM, &
                       PMPIInsert%send_request(iProc,1), IERROR)
        !--- MPI_IRECV lengths of lists of particles entering local mesh
        CALL MPI_IRECV(PMPIInsert%nbrOfRecvParticlesEmission(iProc), 1, MPI_INTEGER, iProc, 1011+FractNbr, PMPIVAR%COMM, &
                       PMPIInsert%recv_request(iProc,1), IERROR)
        IF (PMPIInsert%nbrOfSendParticlesEmission(iProc).GT.0) THEN
          ALLOCATE( PMPIInsert%send_message(iProc)%content(1:3*PMPIInsert%nbrOfSendParticlesEmission(iProc)), STAT=allocStat )
        END IF
      END DO
      PMPIInsert%nbrOfSendParticlesEmission(:)=0
      DO i=1,chunkSize
        CellX = INT((particle_positions(i*3-2)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
        CellY = INT((particle_positions(i*3-1)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
        CellZ = INT((particle_positions(i*3  )-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
        InsideMyBGM=.TRUE.
        IF ((CellX.GT.GEO%FIBGMimax).OR.(CellX.LT.GEO%FIBGMimin) .OR. &
            (CellY.GT.GEO%FIBGMkmax).OR.(CellY.LT.GEO%FIBGMkmin) .OR. &
            (CellZ.GT.GEO%FIBGMlmax).OR.(CellZ.LT.GEO%FIBGMlmin)) THEN
          InsideMyBGM=.FALSE.
        END If
        IF (InsideMyBGM) THEN
          IF (.NOT.ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) InsideMyBGM=.FALSE.
        END IF
        IF (InsideMyBGM) THEN
          DO j=2,GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(1)+1
            iProc=GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(j)
            PMPIInsert%nbrOfSendParticlesEmission(iProc)=PMPIInsert%nbrOfSendParticlesEmission(iProc)+1
            k=PMPIInsert%nbrOfSendParticlesEmission(iProc)
            PMPIInsert%send_message(iProc)%content(k*3-2:k*3)=particle_positions(i*3-2:i*3)
          END DO
          PMPIInsert%nbrOfSendParticlesEmission(PMPIVAR%iProc)=PMPIInsert%nbrOfSendParticlesEmission(PMPIVAR%iProc)+1
          k=PMPIInsert%nbrOfSendParticlesEmission(PMPIVAR%iProc)
          PMPIInsert%send_message(PMPIVAR%iProc)%content(k*3-2:k*3)=particle_positions(i*3-2:i*3)
        ELSE
          DO iProc=0,PMPIVAR%nProcs-1
   !         IF (iProc.EQ.PMPIVAR%iProc) CYCLE
            PMPIInsert%nbrOfSendParticlesEmission(iProc)=PMPIInsert%nbrOfSendParticlesEmission(iProc)+1
            k=PMPIInsert%nbrOfSendParticlesEmission(iProc)
            PMPIInsert%send_message(iProc)%content(k*3-2:k*3)=particle_positions(i*3-2:i*3)
          END DO
        END IF
      END DO
      DEALLOCATE(particle_positions, STAT=allocStat)
      DO iProc=0,PMPIVAR%nProcs-1
        !--- (non-blocking:) send messages to all procs receiving particles from myself
        IF (PMPIInsert%nbrOfSendParticlesEmission(iProc).GT.0) THEN
          CALL MPI_ISEND(PMPIInsert%send_message(iProc)%content, 3*PMPIInsert%nbrOfSendParticlesEmission(iProc),& 
            MPI_DOUBLE_PRECISION, iProc, 1022+FractNbr, PMPIVAR%COMM, PMPIInsert%send_request(iProc,2), IERROR)
        END IF
      END DO
   END IF
ELSE ! mode.NE.1:
!--- RECEIVE:
   IF(nChunks.EQ.1) THEN
      ALLOCATE(particle_positions(1:chunkSize*3), STAT=allocStat)
      IF(PMPIVAR%iProc .EQ. 0) THEN
         particle_positions(:)=PMPIInsert%send_message(0)%content(:)
         DEALLOCATE( PMPIInsert%send_message(0)%content )      
         DEALLOCATE( PMPIInsert%send_message )
      END IF
      CALL MPI_BCAST(particle_positions, chunkSize*3, MPI_DOUBLE_PRECISION,0,PMPIVAR%COMM,IERROR)
   ELSE   
      DO iProc=0,PMPIVAR%nProcs-1
        CALL MPI_WAIT(PMPIInsert%recv_request(iProc,1),msg_status(:),IERROR)
      END DO
      k=SUM(PMPIInsert%nbrOfRecvParticlesEmission)
      ALLOCATE(particle_positions(1:k*3), STAT=allocStat)
      k=0
      DO iProc=0,PMPIVAR%nProcs-1
        IF (PMPIInsert%nbrOfRecvParticlesEmission(iProc).GT.0) THEN
        !--- MPI_IRECV lengths of lists of particles entering local mesh
          CALL MPI_IRECV(particle_positions(k*3+1), 3*PMPIInsert%nbrOfRecvParticlesEmission(iProc),&
                         MPI_DOUBLE_PRECISION, iProc, 1022+FractNbr, PMPIVAR%COMM, PMPIInsert%recv_request(iProc,2), IERROR)
          CALL MPI_WAIT(PMPIInsert%recv_request(iProc,2),msg_status(:),IERROR)
          k=k+PMPIInsert%nbrOfRecvParticlesEmission(iProc)
        END IF
      END DO
      DEALLOCATE( PMPIInsert%nbrOfRecvParticlesEmission )
      DEALLOCATE( PMPIInsert%recv_request )
      DO iProc=0,PMPIVAR%nProcs-1
        CALL MPI_WAIT(PMPIInsert%send_request(iProc,1),msg_status(:),IERROR)
        IF (PMPIInsert%nbrOfSendParticlesEmission(iProc).GT.0) THEN
          CALL MPI_WAIT(PMPIInsert%send_request(iProc,2),msg_status(:),IERROR)
          DEALLOCATE( PMPIInsert%send_message(iProc)%content )
        END IF
      END DO
      DEALLOCATE( PMPIInsert%nbrOfSendParticlesEmission )
      DEALLOCATE( PMPIInsert%send_message )
      DEALLOCATE( PMPIInsert%send_request )
      chunkSize=k
      nChunks=1
   END IF
#endif
   ! each process checks which particle can be matched to its elements, counting the elements inside (local particles)
!   WRITE(*,*)'locating',chunkSize,'*',nChunks,' particles...'
!   WRITE(UNIT=debugFileName,FMT='(A,I2.2)')'prtcls_',PMPIVAR%iProc
!   OPEN(UNIT=130+PMPIVAR%iProc,FILE=debugFileName)
!   DO i=1,chunkSize*nChunks
!      WRITE(130+PMPIVAR%iProc,'(3(ES15.8))')particle_positions(i*3-2:i*3)
!   END DO
!   CLOSE(130+PMPIVAR%iProc)
   mySumOfMatchedParticles=0
   ParticleIndexNbr = 1
   DO i=1,chunkSize*nChunks
      IF ((i.EQ.1).OR.PDM%ParticleInside(ParticleIndexNbr)) THEN
         ParticleIndexNbr = PDM%nextFreePosition(mySumOfMatchedParticles + 1 &
                                               + PDM%CurrentNextFreePosition)
      END IF
      IF (ParticleIndexNbr .ne. 0) THEN
         PartState(ParticleIndexNbr,1:3) = particle_positions(i*3-2:i*3)
         PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
         CALL SingleParticleToExactElement(ParticleIndexNbr)
         IF (PDM%ParticleInside(ParticleIndexNbr)) THEN
            mySumOfMatchedParticles = mySumOfMatchedParticles + 1
         ELSE
            PDM%ParticleInside(ParticleIndexNbr) = .FALSE.
         END IF
      ELSE
         WRITE(*,*)'ERROR in SetParticlePosition:'
         WRITE(*,*)'ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?'
         STOP
      END IF
   END DO

   IF(MOD(iter,IterDisplayStep).EQ.0) THEN
#ifdef MPI
  !   WRITE(*,*)'mySumOfMatchedParticles=',mySumOfMatchedParticles
     ! check the sum of the matched particles: did each particle find its "home"-CPU?
     CALL MPI_ALLREDUCE(mySumOfMatchedParticles, sumOfMatchedParticles, 1, MPI_INTEGER, MPI_SUM, PMPIVAR%COMM, IERROR)
#else
     ! im seriellen Fall kommen alle Partikel auf einen CPU,
     ! daher ist PIC%maxParticleNumber die harte Grenze
     sumOfMatchedParticles = mySumOfMatchedParticles
#endif

#ifdef MPI
     IF(PMPIVAR%iProc.EQ.0) THEN
#endif
     IF (nbrOfParticle .GT. sumOfMatchedParticles) THEN
        WRITE(*,*)'WARNING in ParticleEmission_parallel:'
        WRITE(*,'(A,I0)')'Fraction Nbr: ', FractNbr
        WRITE(*,'(A,I7,A)')'matched only ', sumOfMatchedParticles, ' particles'
        WRITE(*,'(A,I7,A)')'when ', NbrOfParticle, ' particles were required!'
     ELSE IF (nbrOfParticle .LT. sumOfMatchedParticles) THEN
        WRITE(*,*)'ERROR in ParticleEmission_parallel:'
        WRITE(*,'(A,I0)')'Fraction Nbr: ', FractNbr
        WRITE(*,'(A,I7,A)')'matched ', sumOfMatchedParticles, ' particles'
        WRITE(*,'(A,I7,A)')'when ', NbrOfParticle, ' particles were required!'
#if (PP_TimeDiscMethod==1000)
!      STOP
#else
        STOP
#endif
     ELSE IF (nbrOfParticle .EQ. sumOfMatchedParticles) THEN
!      WRITE(*,'(A,I0)')'Fraction Nbr: ', FractNbr
!      WRITE(*,'(A,I0,A)')'ParticleEmission_parallel: matched all (',NbrOfParticle,') particles!'
     END IF
#ifdef MPI
     END IF ! PMPIVAR%iProc.EQ.0
#endif
  END IF ! IterDisplayStep
   ! Return the *local* NbrOfParticle so that the following Routines only fill in
   ! the values for the local particles
  NbrOfParticle = mySumOfMatchedParticles

  DEALLOCATE( particle_positions, STAT=allocStat )
  IF (allocStat .NE. 0) THEN
     WRITE(*,*)'ERROR in ParticleEmission_parallel: cannot deallocate particle_positions!'
     STOP
  END IF
#ifdef MPI
END IF ! NbrOfParticle.LE.0
#endif
   RETURN
END SUBROUTINE SetParticlePosition

SUBROUTINE SetParticleVelocity(FractNbr,iInit,NbrOfParticle)                                             !
!===================================================================================================================================
! Determine the particle velocity of each inserted particle
!===================================================================================================================================
  USE MOD_Timedisc_Vars , ONLY : dt
  USE MOD_Particle_Vars
  USE MOD_PIC_Vars
  USE MOD_Mesh_Vars,      ONLY : nElems
!----------------------------------------------------------------------------------------------------------------------------------
   IMPLICIT NONE                                                                                   !
!----------------------------------------------------------------------------------------------------------------------------------
! argument list declaration                                                                        !
   INTEGER                          :: FractNbr,iInit,NbrOfParticle                                      !
! Local variable declaration                                                                       !
   INTEGER                          :: i,j,k,PositionNbr,PositionNbr2,iPartStart,iPartEnd          !
   REAL                             :: Radius(3), n_vec(3), tan_vec(3), cdis, Velo1, Angle, Velo2, f, Pos(3) !
   REAL                             :: Vec3D(3), VeloMaxProb, RandVal(3), RanNum, GaussRan(6), RandTmp(6) !
   REAL                             :: II(3,3),JJ(3,3),NN(3,3)                                     !
   INTEGER                          :: elemNbr , distnum                                                    !
   REAL                             :: mu(1:nElems,1:3),sig(1:nElems,1:3),NPart(1:nElems)          !
   REAL                             :: xlen,Wnbr,alpha                                             !
   REAL                             :: r1,r2,x_1,x_2,y_1,y_2,a,b,e,g,x_01,x_02,y_01,y_02, RandVal1
   REAL, PARAMETER                  :: PI=3.14159265358979323846	                                 !
   REAL                             :: Velosq, v_sum(3), v2_sum, maxwellfac
!----------------------------------------------------------------------------------------------------------------------------------
   INTENT(IN)                       :: FractNbr,iInit                                                   !
   INTENT(INOUT)                    :: NbrOfParticle                                               !
!----------------------------------------------------------------------------------------------------------------------------------
   IF(NbrOfParticle.lt.1) RETURN
 
!   IF(NbrOfParticle.gt.PIC%maxParticleNumber)THEN
!      NbrOfParticle = PIC%maxParticleNumber
!   END IF

   SELECT CASE(Species(FractNbr)%Init(iInit)%velocityDistribution)
   CASE('random')
     i = 1
     DO WHILE (i .le. NbrOfParticle)
       PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
       IF (PositionNbr .ne. 0) THEN
         CALL RANDOM_NUMBER(RandVal)
         RandVal(:) = RandVal(:) - 0.5
         RandVal(:) = RandVal(:)/SQRT(RandVal(1)**2+RandVal(2)**2+RandVal(3)**2)
         PartState(PositionNbr,4:6) = RandVal(1:3) * Species(FractNbr)%Init(iInit)%VeloIC
       END IF
       i = i + 1
     END DO
!   CASE('EOC_Test')
!      ! for leapfrog EOC test velo has to be set half dt before ICPos.
!      Radius = Species(1)%BasePointIC 
!      IF (Species(1)%BasePointIC(1) > 0. ) THEN
!        n_vec = (/0.,0.,-1./) 
!        Angle = 0.
!        Angle = Angle - PIC%GyrationFrequency * dt * 0.5
!        Radius(1) = cos(Angle)
!        Radius(2) = sin(Angle)   
!      ELSEIF (Species(1)%BasePointIC(3) > 0. ) THEN  
!        n_vec = (/0.,1.,0./) 
!        Angle = PI
!        Angle = Angle - PIC%GyrationFrequency * dt * 0.5
!        Radius(1) = sin(Angle)
!        Radius(3) = cos(Angle)   
!      END IF
!      tan_vec(1) = Radius(2)*n_vec(3) - Radius(3)*n_vec(2)
!      tan_vec(2) = Radius(3)*n_vec(1) - Radius(1)*n_vec(3)
!      tan_vec(3) = Radius(1)*n_vec(2) - Radius(2)*n_vec(1)
!      PartState(1,4:6) = tan_vec(1:3) * Species(1)%VeloIC
   CASE('constant')
      i = 1
      DO WHILE (i .le. NbrOfParticle)
         PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
         IF (PositionNbr .ne. 0) THEN
            PartState(PositionNbr,4:6) = Species(FractNbr)%Init(iInit)%VeloVecIC(1:3) * Species(FractNbr)%Init(iInit)%VeloIC
         END IF
         i = i + 1
      END DO
   CASE('radial_constant')
      i = 1
      DO WHILE (i .le. NbrOfParticle)
         PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
         IF (PositionNbr .ne. 0) THEN
            Radius(1:3) = PartState(PositionNbr,1:3) - Species(FractNbr)%Init(iInit)%BasePointIC(1:3)
            !  Unity radius
            !Radius(1:3) = Radius(1:3) / Species(FractNbr)%Init(iInit)%RadiusIC
            Radius(1:3) = Radius(1:3) / SQRT(Radius(1)**2+Radius(2)**2+Radius(3)**2) 
            PartState(PositionNbr,4:6) = Radius(1:3) * Species(FractNbr)%Init(iInit)%VeloIC
         END IF
         i = i + 1
         !print*, SQRT(PartState(PositionNbr,4)**2+PartState(PositionNbr,5)**2+PartState(PositionNbr,6)**2)
      END DO
   CASE('tangential_constant')
      i = 1
      DO WHILE (i .le. NbrOfParticle)
         PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
         IF (PositionNbr .ne. 0) THEN
            Radius(1:3) = PartState(PositionNbr,1:3) - Species(FractNbr)%Init(iInit)%BasePointIC(1:3)
            !  Normal Vector of circle
            n_vec(1:3) = Species(FractNbr)%Init(iInit)%NormalIC(1:3)
            ! If we're doing Leapfrog, then use velocities from half-timestep before
            IF (ParticlePushMethod.EQ.'boris_leap_frog_scheme') THEN
              Angle = 0.5 * dt * Species(FractNbr)%Init(iInit)%VeloIC / Species(FractNbr)%Init(iInit)%RadiusIC ! 0.5*dt*(v/r)
              JJ(1,1:3) = (/   0.,-n_vec(3), n_vec(2)/)
              JJ(2,1:3) = (/ n_vec(3),   0.,-n_vec(1)/)
              JJ(3,1:3) = (/-n_vec(2), n_vec(1),   0./)
              II(1,1:3) = (/1.,0.,0./)
              II(2,1:3) = (/0.,1.,0./)
              II(3,1:3) = (/0.,0.,1./)
              forall(j=1:3) NN(:,j) = n_vec(:)*n_vec(j)
              Radius = MATMUL( NN+cos(Angle)*(II-NN)+sin(Angle)*JJ , Radius )
            END IF
            !  Unity radius
            Radius(1:3) = Radius(1:3) / SQRT(Radius(1)**2+Radius(2)**2+Radius(3)**2)
            !  Vector Product rxn
            tan_vec(1) = Radius(2)*n_vec(3) - Radius(3)*n_vec(2)
            tan_vec(2) = Radius(3)*n_vec(1) - Radius(1)*n_vec(3)
            tan_vec(3) = Radius(1)*n_vec(2) - Radius(2)*n_vec(1)
            ! If Gyrotron resonator: Add velocity in normal direction!
            IF (Species(FractNbr)%Init(iInit)%alpha .gt. 0.) THEN 
              n_vec = n_vec * ( 1 / Species(FractNbr)%Init(iInit)%alpha )
            ELSE 
              n_vec = 0
            END IF
            !  And finally the velocities
            PartState(PositionNbr,4:6) = (tan_vec(1:3) + n_vec(1:3)) * Species(FractNbr)%Init(iInit)%VeloIC
         END IF
         i = i + 1
      END DO
   CASE('gyrotron_circle')
      i = 1
      DO WHILE (i .le. NbrOfParticle)
         PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
         IF (PositionNbr .ne. 0) THEN
         !! Position of particle on gyro circle changed in SetParticlePosition.F90: Problem
         !! We don't have the radius-vector any more. Thus transport the radius vector from there to here.
         ! Or do Alternative way: Hack the radius by intersecting two circles (big IC and small gyro circle)
           r1 = Species(FractNbr)%Init(iInit)%RadiusIC
           r2 = Species(FractNbr)%Init(iInit)%RadiusICGyro
           x_1 = 0.
           y_1 = 0.
           x_2 = PartState(PositionNbr,1)
           y_2 = PartState(PositionNbr,2)
           IF (x_1 .eq. x_2) THEN
             a = (x_1 - x_2)/(y_2-y_1)
             b = ((r1**2-r2**2)-(x_1**2-x_2**2)-(y_1**2-y_2**2))&
                 & /(2.*y_2-2.*y_1)
             e = (a**2+1.)
             f = (2.*a*(b-y_1))-2.*x_1
             g = (b-y_1)**2-r1**2+x_1**2
             ! intersection points
             x_01 = (-f + SQRT(ABS(f**2 - 4. * e * g)))/(2.*e) ! the term in SQRT can be -0.0 , therefore the ABS
             x_02 = (-f - SQRT(ABS(f**2 - 4. * e * g)))/(2.*e) ! the term in SQRT can be -0.0 , therefore the ABS
             y_01 = x_01 * a + b
             y_02 = x_02 * a + b
           ELSE
             a = (y_1 - y_2)/(x_2-x_1)
             b = ((r1**2 - r2**2)-(x_1**2-x_2**2)-(y_1**2-y_2**2))&
                  & /(2.*x_2 - 2. * x_1)
             e = (a**2 + 1.)
             f = 2. * a * (b - x_1) -2 *y_1
             g = (b-x_1)**2 - r1**2 + y_1**2
             y_01 = (-f + SQRT(ABS(f**2 - 4. * e * g)))/(2.*e) ! the term in SQRT can be -0.0 , therefore the ABS
             y_02 = (-f - SQRT(ABS(f**2 - 4. * e * g)))/(2.*e) ! the term in SQRT can be -0.0 , therefore the ABS
             x_01 = y_01 * a + b
             x_02 = y_02 * a + b
           END IF
           CALL RANDOM_NUMBER(RandVal1)
           IF (RandVal1 .ge. 0.5) THEN
             Radius(1) = PartState(PositionNbr,1) - x_01 
             Radius(2) = PartState(PositionNbr,2) - y_01 
           ELSE
             Radius(1) = PartState(PositionNbr,1) - x_02
             Radius(2) = PartState(PositionNbr,2) - y_02
           END IF     
         
            Radius(3) = 0.
            !Check if Radius has correct length
            IF ((SQRT(Radius(1)**2+Radius(2)**2)-r1).ge.1E-15) THEN
              WRITE(*,*)"Error in setparticle velocity, gyrotron circle. &
                        & Radius to big after intersection."
            END IF
            !  Normal Vector of circle
            n_vec(1:3) = Species(FractNbr)%Init(iInit)%NormalIC(1:3)

            ! If we're doing Leapfrog, then use velocities from half-timestep before. This only applies in 
            ! x- and y-direction. z has allways same velo. 
!            IF (ParticlePushMethod.EQ.'boris_leap_frog_scheme') THEN
!              ! get angle of part on gyrocircle
!              Angle = ACOS(Radius(1)/Species(1)%RadiusICGyro)
!              IF (Radius(2).LE.0) THEN
!                Angle = 2*PI-Angle
!              END IF
!              ! shift position angle half dt back in time (as particle moves clockwise,
!              ! we add dalpha in ccw direction)
!              Angle = Angle + PIC%GyrationFrequency * dt * 0.5 * PIC%GyroVecDirSIGN  
!              Radius(1) = cos(Angle)
!              Radius(2) = sin(Angle)                 
!            END IF
            !  Unity radius
            Radius(1:3) = Radius(1:3) / SQRT(Radius(1)**2+Radius(2)**2+Radius(3)**2)
            !  Vector Product rxn
            tan_vec(1) = Radius(2)*n_vec(3)*PIC%GyroVecDirSIGN - Radius(3)*n_vec(2)
            tan_vec(2) = Radius(3)*n_vec(1) - Radius(1)*n_vec(3) *PIC%GyroVecDirSIGN
            tan_vec(3) = Radius(1)*n_vec(2) - Radius(2)*n_vec(1)
            ! If Gyrotron resonator: Add velocity in normal direction!
            IF (Species(FractNbr)%Init(iInit)%alpha .gt. 0.) THEN 
              n_vec = n_vec * ( 1. / Species(FractNbr)%Init(iInit)%alpha )
            ELSE 
              n_vec = 0.
            END IF
            !  And finally the velocities
            PartState(PositionNbr,4:6) = (tan_vec(1:3) + n_vec(1:3)) * Species(FractNbr)%Init(iInit)%VeloIC
            IF (ABS(SQRT(PartState(PositionNbr,4)*PartState(PositionNbr,4) &
                       + PartState(PositionNbr,5)*PartState(PositionNbr,5))&
                       - Species(FractNbr)%Init(iInit)%VeloIC) .GT. 10.) THEN
              WRITE(*,'(A,3(E21.15,X))')'ERROR! Velocity=', PartState(PositionNbr,4:6)
              STOP
            END If
            IF (PartState(PositionNbr,4).NE.PartState(PositionNbr,4) .OR. &
                PartState(PositionNbr,5).NE.PartState(PositionNbr,5) .OR. &
                PartState(PositionNbr,6).NE.PartState(PositionNbr,6)     ) THEN
              WRITE(*,'(A,3(E21.15,X))')'WARNING:! NaN: Velocity=', PartState(PositionNbr,4:6)
            END If
         END IF
         i = i + 1
      END DO
   CASE('maxwell')
      v_sum(1:3) = 0.0
      v2_sum = 0.0
   
      i = 1
      DO WHILE (i .le. NbrOfParticle)
         PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
         IF (PositionNbr .ne. 0) THEN
            DO distnum = 1, 3
              CALL RANDOM_NUMBER(RandVal)
              Velo1 = 2.0*RandVal(1)-1.0
              Velo2 = 2.0*RandVal(2)-1.0
              Velosq= Velo1**2+Velo2**2
              DO WHILE ((Velosq.LE.0).OR.(Velosq.GE.1))
                CALL RANDOM_NUMBER(RandVal)
                Velo1 = 2.0*RandVal(1)-1.0
                Velo2 = 2.0*RandVal(2)-1.0
                Velosq= Velo1**2+Velo2**2
              END DO
              Vec3D(distnum) = Velo1*SQRT(-2*LOG(Velosq)/Velosq)
            END DO                    
            PartState(PositionNbr,4:6) = Vec3D(1:3)
            v_sum(1:3) = v_sum(1:3) + Vec3D(1:3)
            v2_sum = v2_sum + Vec3D(1)**2+Vec3D(2)**2+Vec3D(3)**2
         END IF
         i = i + 1
      END DO
      v_sum(1:3) = v_sum(1:3) / NbrOfParticle
      v2_sum = v2_sum / NbrOfParticle
      maxwellfac = SQRT(3. * 1.380650524E-23 * Species(FractNbr)%Init(iInit)%MWTemperatureIC / &              ! velocity of maximum
                     (Species(FractNbr)%MassIC*v2_sum))

      i = 1
      DO WHILE (i .le. NbrOfParticle)
         PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
         IF (PositionNbr .ne. 0) THEN
           PartState(PositionNbr,4:6) = (PartState(PositionNbr,4:6) - v_sum(1:3)) * maxwellfac &
                                        + Species(FractNbr)%Init(iInit)%VeloIC *Species(FractNbr)%Init(iInit)%VeloVecIC(1:3)        
         END IF
         i = i + 1
      END DO
   END SELECT
 RETURN
END SUBROUTINE SetParticleVelocity

SUBROUTINE SetParticleChargeAndMass(FractNbr,NbrOfParticle)                                        
!===================================================================================================================================
! And partilces mass and charge
!===================================================================================================================================
  USE MOD_Particle_Vars,    ONLY : PDM, PartSpecies
!----------------------------------------------------------------------------------------------------------------------------------
  IMPLICIT NONE                                                                                    
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST
  INTEGER                          :: FractNbr,NbrOfParticle                                       
! LOCAL VARIABLES
  INTEGER                          :: i,PositionNbr                                                
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTENT(IN)                       :: FractNbr                                                     
!----------------------------------------------------------------------------------------------------------------------------------
! IN/OUT VARiABLES
  INTENT(INOUT)                    :: NbrOfParticle                                                
!----------------------------------------------------------------------------------------------------------------------------------
  IF(NbrOfParticle.gt.PDM%maxParticleNumber)THEN
    NbrOfParticle = PDM%maxParticleNumber
  END IF
  i = 1
  DO WHILE (i .le. NbrOfParticle)
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr .ne. 0) THEN
      PartSpecies(PositionNbr) = FractNbr
    END IF
    i = i + 1
  END DO
 RETURN
END SUBROUTINE SetParticleChargeAndMass

SUBROUTINE SetParticleMPF(FractNbr,NbrOfParticle)                                        !
!===================================================================================================================================
! finaly, set the MPF
!===================================================================================================================================
  USE MOD_Particle_Vars,    ONLY : PDM, PartMPF, Species
!===================================================================================================================================
  IMPLICIT NONE                                                                                    !
!----------------------------------------------------------------------------------------------------------------------------------
! argument list declaration                                                                        !
  INTEGER                          :: FractNbr,NbrOfParticle                                       !
! Local variable declaration                                                                       !
  INTEGER                          :: i,PositionNbr                                                !
!----------------------------------------------------------------------------------------------------------------------------------
  INTENT(IN)                       :: FractNbr                                                     !
  INTENT(INOUT)                    :: NbrOfParticle                                                !
!----------------------------------------------------------------------------------------------------------------------------------
  IF(NbrOfParticle.gt.PDM%maxParticleNumber)THEN
    NbrOfParticle = PDM%maxParticleNumber
  END IF
  i = 1
  DO WHILE (i .le. NbrOfParticle)
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr .ne. 0) THEN
      PartMPF(PositionNbr) = Species(FractNbr)%MacroParticleFactor
    END IF
    i = i + 1
  END DO
 RETURN
END SUBROUTINE SetParticleMPF

SUBROUTINE ParticleInsertingCellPressure(iSpec,iInit,NbrOfParticle)
!===================================================================================================================================
! Insert constant cell pressure particles (and remove additionals)
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars
USE MOD_part_MPFtools,   ONLY : MapToGEO
USE MOD_BoundaryTools,   ONLY : SingleParticleToExactElement, ParticleInsideQuad3D
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! Argument list variables
  INTEGER               :: iSpec, iInit, NbrOfParticle
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: iElem, Elem, iPart, i, NbrPartsInCell, iNode, NbrNewParts
  INTEGER               :: ParticleIndexNbr
  REAL                  :: PartDiff, PartDiffRest, RandVal, RandVal3(1:3), P(1:3,1:8)
  INTEGER, ALLOCATABLE  :: PartsInCell(:)
  LOGICAL               :: InElementCheck
  REAL                  :: det(16)
!----------------------------------------------------------------------------------------------------------------------------------
! INTENTIONS
  INTENT(IN)            :: iSpec, iInit
  INTENT(OUT)           :: NbrOfParticle
!==================================================================================================

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
    ! Build array for element nodes
    DO iNode = 1,8
      P(1:3,iNode) = GEO%NodeCoords(1:3,GEO%ElemToNodeID(iNode,Elem))
    END DO
    ! insert particles (positions)
    DO i = 1, NbrNewParts
      ! set random position in -1,1 space
      CALL RANDOM_NUMBER(RandVal3)
      RandVal3 = RandVal3 * 2.0 - 1.0 
      ParticleIndexNbr = PDM%nextFreePosition(PDM%CurrentNextFreePosition + i + NbrOfParticle)
      IF (ParticleIndexNbr.NE.0) THEN
        PartState(ParticleIndexNbr, 1:3) = MapToGeo(RandVal3,P)
        PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
        CALL ParticleInsideQuad3D(ParticleIndexNbr,Elem,InElementCheck,det)
        IF (InElementCheck) THEN
          PEM%Element(ParticleIndexNbr) = Elem
        ELSE
          CALL SingleParticleToExactElement(ParticleIndexNbr)
        END IF
      ELSE
        WRITE(*,*)'ERROR in ParticleInsertingCellPressure:'
        WRITE(*,*)'ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?'
        STOP
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
USE MOD_Particle_Vars
!USE MOD_PIC_Vars !hoechstwahrscheinlich nicht benoetigt
USE MOD_part_MPFtools,   ONLY : MapToGEO
USE MOD_BoundaryTools,   ONLY : SingleParticleToExactElement, ParticleInsideQuad3D
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! Argument list variables
  INTEGER               :: iSpec, iInit, NbrOfParticle
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: iElem, Elem, iPart, i, NbrPartsInCell, iNode, distnum
  INTEGER               :: ParticleIndexNbr
  REAL                  :: PartDiff, PartDiffRest, RandVal, RandVal3(1:3), P(1:3,1:8)
  INTEGER, ALLOCATABLE  :: PartsInCell(:)
  LOGICAL               :: InElementCheck
  REAL                  :: det(16)
  REAL                  :: Velo1, Velo2, Velosq, v_sum(3), v2_sum, maxwellfac
  REAL                  :: Vec3D(3), RandVal3D(3)
!----------------------------------------------------------------------------------------------------------------------------------
! INTENTIONS
  INTENT(IN)            :: iSpec, iInit
  INTENT(OUT)           :: NbrOfParticle
!==================================================================================================
NbrOfParticle = 0
IF (Species(iSpec)%Init(iInit)%velocityDistribution.NE.'maxwell') THEN
  PRINT*,'Only maxwell implemented yet!'
  STOP
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
    ! Build array for element nodes
    DO iNode = 1,8
      P(1:3,iNode) = GEO%NodeCoords(1:3,GEO%ElemToNodeID(iNode,Elem))
    END DO
    ! insert particles (positions and velocities)
    v_sum(1:3) = 0.0
    v2_sum = 0.0
    DO i = 1, NbrPartsInCell
      ! set random position in -1,1 space
      CALL RANDOM_NUMBER(RandVal3)
      RandVal3 = RandVal3 * 2.0 - 1.0 
      ParticleIndexNbr = PDM%nextFreePosition(PDM%CurrentNextFreePosition + i + NbrOfParticle)
      IF (ParticleIndexNbr.NE.0) THEN
        PartState(ParticleIndexNbr, 1:3) = MapToGeo(RandVal3,P)
        PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
        CALL ParticleInsideQuad3D(ParticleIndexNbr,Elem,InElementCheck,det)
        IF (InElementCheck) THEN
          PEM%Element(ParticleIndexNbr) = Elem
        ELSE
          CALL SingleParticleToExactElement(ParticleIndexNbr)
        END IF
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
        PartState(ParticleIndexNbr,4:6) = Vec3D(1:3)
        v_sum(1:3) = v_sum(1:3) + Vec3D(1:3)
        v2_sum = v2_sum + Vec3D(1)**2+Vec3D(2)**2+Vec3D(3)**2
      ELSE
        WRITE(*,*)'ERROR in ParticleInsertingCellPressureOut:'
        WRITE(*,*)'ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?'
        STOP
      END IF
    END DO
    v_sum(1:3) = v_sum(1:3) / (NbrPartsInCell+1) !+1 correct?
    v2_sum = v2_sum / (NbrPartsInCell+1)         !+1 correct?
    !maxwellfactor from new calculated values (no vibrational DOF implemented, equilibirium assumed)
    !Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(iElem,:)
    maxwellfac = SQRT(3. * 1.380650524E-23 * Species(iSpec)%Init(iInit)%ConstantPressure &
                 / (Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(iElem,4) * BoltzmannConst) / & ! velocity of maximum
                 (Species(iSpec)%MassIC*v2_sum))                                                           ! T = p_o / (<n>*k)
    ! particel velocity (maxwell, part 2)
    DO i = 1, NbrPartsInCell
      ParticleIndexNbr = PDM%nextFreePosition(PDM%CurrentNextFreePosition + i + NbrOfParticle)
      IF (ParticleIndexNbr .ne. 0) THEN
        PartState(ParticleIndexNbr,4:6) = (PartState(ParticleIndexNbr,4:6) - v_sum(1:3)) * maxwellfac &  !macro velocity:
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
  USE MOD_Particle_Vars,         ONLY : PartState,usevMPF,GEO,Species,PartSpecies,BoltzmannConst,usevMPF,PartMPF,PDM
  USE MOD_DSMC_Vars,             ONLY : SpecDSMC
  USE MOD_TimeDisc_Vars,         ONLY : iter
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER,INTENT(IN) :: iSpec, iInit, iElem, ElemSamp
  INTEGER,INTENT(IN) :: PartsInCell(:)
  INTEGER,INTENT(INOUT) :: NbrPartsInCell
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

  INTEGER             :: iPart, nPart, iPartIndx
  REAL                :: MPFSum, WeightFak, kappa_part, AvogadroConst, RandVal, RealnumberNewParts, RealnumberNewPartsRest
  REAL                :: Samp_V2(3), Samp_Temp(4), OldConstPressureSamp(6)
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
         + PartState(iPartIndx,4:6) * WeightFak
    Samp_V2(:)                      = Samp_V2(:) + PartState(iPartIndx,4:6)**2 * WeightFak !vi**2 =vi**2 + vi**2*W
    MPFSum                          = MPFSum + WeightFak                                   !MPFsum = MPFsum + W
    PDM%ParticleInside(iPartIndx)=.false. !remove particle
  END DO

  !Calculation of specific heat ratio (no vibrational DOF -> only at low temperatures !!!)
  IF(SpecDSMC(PartSpecies(iPartIndx))%InterID.EQ.2) THEN
    kappa_part=1.4
  ELSE IF(SpecDSMC(PartSpecies(iPartIndx))%InterID.EQ.1) THEN
    kappa_part=5.0/3.0
  ELSE
    WRITE(*,*)'Wrong PartSpecies!'
    STOP
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
  RealnumberNewPartsRest = RealnumberNewParts - INT(RealnumberNewParts)
  IF(RealnumberNewParts.GT.0) THEN
    CALL RANDOM_NUMBER(RandVal)
    IF(RealnumberNewPartsRest.GT.RandVal) RealnumberNewParts = RealnumberNewParts + 1.0
    NbrPartsInCell = INT(RealnumberNewParts)
  ELSE
    NbrPartsInCell = 0
  END IF

ELSE ! no particles in cell!
  PRINT*,'YOU NEED MORE PARTICLES INSIDE THE OUTFLOW REGION!!!'
  STOP
END IF

END SUBROUTINE ParticleInsertingPressureOut_Sampling

END MODULE MOD_part_emission
