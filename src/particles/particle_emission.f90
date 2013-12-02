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
USE MOD_DSMC_Vars,      ONLY : useDSMC
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
  CHARACTER(40)                          :: SpaceIC                          ! specifying Keyword for Particle Space condition
  CHARACTER(30)                          :: velocityDistribution             ! specifying keyword for velocity distribution
  INTEGER(8)                             :: initialParticleNumber            ! Number of Particles at time 0.0
  REAL                                   :: RadiusIC
  REAL                                   :: Radius2IC
  REAL                                   :: RadiusICGyro
  REAL                                   :: NormalIC(3)                      ! Normal / Orientation of circle
  REAL                                   :: BasePointIC(3)                   ! base point for IC cuboid and IC sphere
  REAL                                   :: BaseVector1IC(3)                 ! first base vector for IC cuboid
  REAL                                   :: BaseVector2IC(3)                 ! second base vector for IC cuboid
  REAL                                   :: CuboidHeightIC                   ! third measure of cuboid
  REAL                                   :: CylinderHeightIC                 ! third measure of cylinder
  REAL                                   :: VeloIC                           ! velocity for inital Data
  REAL                                   :: VeloVecIC(3)                     ! normalized velocity vector
  REAL                                   :: MWTemperatureIC                  ! Temperature for Maxwell Distribution
  REAL                                   :: alpha                            ! WaveNumber for sin-deviation initiation.
  REAL                                   :: Amplitude                        ! Amplitude for sin-deviation initiation.
  REAL                                   :: WaveNumber                       ! WaveNumber for sin-deviation initiation.
  INTEGER(8)                             :: maxParticleNumberX               ! Maximum Number of all Particles in x direction
  INTEGER(8)                             :: maxParticleNumberY               ! Maximum Number of all Particles in y direction
  INTEGER(8)                             :: maxParticleNumberZ               ! Maximum Number of all Particles in z direction 
!===================================================================================================================================

CALL UpdateNextFreePosition()

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
     ! check whether initial particles are defined twice (old and new method) to prevent erroneous doubling
     ! of particles
     IF ((Species(i)%initialParticleNumber.NE.0).AND.(Species(i)%NumberOfInits.NE.0)) THEN
       WRITE(*,*) 'ERROR in ParticleEmission: Initial emission may only be defined in additional *Init#* blocks'
       WRITE(*,*) 'OR the standard initialisation, not both!'
       STOP
     END IF

     IF (Species(i)%NumberOfInits.EQ.0) THEN
       ! Standard Emission
       NbrOfParticle = Species(i)%initialParticleNumber
#ifdef MPI
       CALL SetParticlePosition(i,NbrOfParticle,1)
       CALL SetParticlePosition(i,NbrOfParticle,2)
#else
       CALL SetParticlePosition(i,NbrOfParticle)
#endif /*MPI*/
       CALL SetParticleVelocity(i,NbrOfParticle)
       CALL SetParticleChargeAndMass(i,NbrOfParticle)
       IF (usevMPF) CALL SetParticleMPF(i,NbrOfParticle)
       !IF (useDSMC) CALL SetParticleIntEnergy(i,NbrOfParticle)
       PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfParticle
       CALL UpdateNextFreePosition()
     ELSE
       ! temporarily store data
       SpaceIC = Species(i)%SpaceIC
       velocityDistribution = Species(i)%velocityDistribution
       initialParticleNumber = Species(i)%initialParticleNumber
       RadiusIC = Species(i)%RadiusIC
       Radius2IC = Species(i)%Radius2IC
       RadiusICGyro = Species(i)%RadiusICGyro
       NormalIC = Species(i)%NormalIC
       BasePointIC = Species(i)%BasePointIC
       BaseVector1IC = Species(i)%BaseVector1IC
       BaseVector2IC = Species(i)%BaseVector2IC
       CuboidHeightIC = Species(i)%CuboidHeightIC
       CylinderHeightIC = Species(i)%CylinderHeightIC
       VeloIC = Species(i)%VeloIC
       VeloVecIC = Species(i)%VeloVecIC
       Amplitude = Species(i)%Amplitude
       WaveNumber = Species(i)%WaveNumber
       maxParticleNumberX = Species(i)%maxParticleNumberX
       maxParticleNumberY = Species(i)%maxParticleNumberY
       maxParticleNumberZ = Species(i)%maxParticleNumberZ
       Alpha = Species(i)%Alpha
       MWTemperatureIC = Species(i)%MWTemperatureIC

       DO iInit = 1, Species(i)%NumberOfInits
         NbrOfParticle = Species(i)%Init(iInit)%initialParticleNumber
         ! to prevent doubling of subroutines and conflicts with earlier versions, the relevant data is 
         ! temporarily copied into the regular species emission data and then re-read from the ini file
         Species(i)%SpaceIC               = Species(i)%Init(iInit)%SpaceIC
         Species(i)%velocityDistribution  = Species(i)%Init(iInit)%velocityDistribution
         Species(i)%initialParticleNumber = Species(i)%Init(iInit)%initialParticleNumber
         Species(i)%RadiusIC              = Species(i)%Init(iInit)%RadiusIC
         Species(i)%Radius2IC             = Species(i)%Init(iInit)%Radius2IC
         Species(i)%RadiusICGyro          = Species(i)%Init(iInit)%RadiusICGyro
         Species(i)%NormalIC              = Species(i)%Init(iInit)%NormalIC
         Species(i)%BasePointIC           = Species(i)%Init(iInit)%BasePointIC
         Species(i)%BaseVector1IC         = Species(i)%Init(iInit)%BaseVector1IC
         Species(i)%BaseVector2IC         = Species(i)%Init(iInit)%BaseVector2IC
         Species(i)%CuboidHeightIC        = Species(i)%Init(iInit)%CuboidHeightIC
         Species(i)%CylinderHeightIC      = Species(i)%Init(iInit)%CylinderHeightIC
         Species(i)%VeloIC                = Species(i)%Init(iInit)%VeloIC
         Species(i)%VeloVecIC             = Species(i)%Init(iInit)%VeloVecIC
         Species(i)%Amplitude             = Species(i)%Init(iInit)%Amplitude
         Species(i)%WaveNumber            = Species(i)%Init(iInit)%WaveNumber
         Species(i)%maxParticleNumberX    = Species(i)%Init(iInit)%maxParticleNumberX
         Species(i)%maxParticleNumberY    = Species(i)%Init(iInit)%maxParticleNumberY
         Species(i)%maxParticleNumberZ    = Species(i)%Init(iInit)%maxParticleNumberZ
         Species(i)%Alpha                 = Species(i)%Init(iInit)%Alpha
         Species(i)%MWTemperatureIC       = Species(i)%Init(iInit)%MWTemperatureIC
#ifdef MPI
         CALL SetParticlePosition(i,NbrOfParticle,1)
         CALL SetParticlePosition(i,NbrOfParticle,2)
#else
         CALL SetParticlePosition(i,NbrOfParticle)
#endif /*MPI*/
         CALL SetParticleVelocity(i,NbrOfParticle)
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
       END DO
       Species(i)%SpaceIC = SpaceIC
       Species(i)%velocityDistribution = velocityDistribution
       Species(i)%initialParticleNumber = initialParticleNumber
       Species(i)%RadiusIC = RadiusIC
       Species(i)%Radius2IC = Radius2IC
       Species(i)%RadiusICGyro = RadiusICGyro
       Species(i)%NormalIC = NormalIC
       Species(i)%BasePointIC = BasePointIC
       Species(i)%BaseVector1IC = BaseVector1IC
       Species(i)%BaseVector2IC = BaseVector2IC
       Species(i)%CuboidHeightIC = CuboidHeightIC
       Species(i)%CylinderHeightIC = CylinderHeightIC
       Species(i)%VeloIC = VeloIC
       Species(i)%VeloVecIC = VeloVecIC
       Species(i)%Amplitude = Amplitude
       Species(i)%WaveNumber = WaveNumber
       Species(i)%maxParticleNumberX = maxParticleNumberX
       Species(i)%maxParticleNumberY = maxParticleNumberY
       Species(i)%maxParticleNumberZ = maxParticleNumberZ
       Species(i)%Alpha = Alpha
       Species(i)%MWTemperatureIC = MWTemperatureIC
     END IF
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
  USE MOD_Timedisc_Vars         , ONLY : dt
  USE MOD_Particle_Vars       ! , ONLY : time
  USE MOD_PIC_Vars
  USE MOD_part_tools             ,ONLY : UpdateNextFreePosition  
  USE MOD_DSMC_Vars              ,ONLY : useDSMC, CollisMode                                   
  USE MOD_DSMC_Init              ,ONLY : DSMC_SetInternalEnr_LauxVFD
#if (PP_TimeDiscMethod==1000)
  USE MOD_LD_Init                ,ONLY : CalcDegreeOfFreedom
  USE MOD_LD_Vars
#endif
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
   INTEGER                          :: i , iPart, PositionNbr                                                          
   INTEGER                , SAVE    :: NbrOfParticle=0                                             
   INTEGER                          :: inserted_Particle_iter,inserted_Particle_time               
   INTEGER                          :: inserted_Particle_diff                                      
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
#ifdef MPI
         IF (mode.NE.2) THEN
#endif
      SELECT CASE(Species(i)%ParticleEmissionType)
      CASE(1) ! Emission Type: Particles per !!!!!SECOND!!!!!!!! (not per ns) 
         inserted_Particle_iter = INT(Species(i)%ParticleEmission * dt)
         inserted_Particle_time = INT(Species(i)%ParticleEmission * (Time + dt))
         inserted_Particle_diff = INT(inserted_Particle_time - Species(i)%InsertedParticle &
                                  - inserted_Particle_iter )
         NbrOfParticle = inserted_Particle_iter + inserted_Particle_diff
         Species(i)%InsertedParticle = Species(i)%InsertedParticle + NbrOfParticle
      CASE(2)    ! Emission Type: Particles per Iteration
         NbrOfParticle = INT(Species(i)%ParticleEmission)
      CASE DEFAULT
         NbrOfParticle = 0
      END SELECT
#ifdef MPI
      CALL SetParticlePosition(i,NbrOfParticle,1)
         END IF
         IF (mode.NE.1) THEN
      CALL SetParticlePosition(i,NbrOfParticle,2)
#else
      CALL SetParticlePosition(i,NbrOfParticle)
#endif
      CALL SetParticleVelocity(i,NbrOfParticle)
      CALL SetParticleChargeAndMass(i,NbrOfParticle)
      IF (usevMPF) CALL SetParticleMPF(i,NbrOfParticle)
      ! define molecule stuff
      IF (useDSMC.AND.(CollisMode.NE.1)) THEN
        iPart = 1
        DO WHILE (iPart .le. NbrOfParticle)
          PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
          IF (PositionNbr .ne. 0) THEN
            CALL DSMC_SetInternalEnr_LauxVFD(i, PositionNbr)
          END IF
          iPart = iPart + 1
        END DO
      END IF
#if (PP_TimeDiscMethod==1000)
      iPart = 1
      DO WHILE (iPart .le. NbrOfParticle)
          PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
          IF (PositionNbr .ne. 0) THEN
            PartStateBulkValues(PositionNbr,1) = Species(i)%VeloVecIC(1) * Species(i)%VeloIC
            PartStateBulkValues(PositionNbr,2) = Species(i)%VeloVecIC(2) * Species(i)%VeloIC
            PartStateBulkValues(PositionNbr,3) = Species(i)%VeloVecIC(3) * Species(i)%VeloIC
            PartStateBulkValues(PositionNbr,4) = Species(i)%MWTemperatureIC
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
   END DO
 RETURN
END SUBROUTINE ParticleInserting
                                                                                                   
#ifdef MPI
SUBROUTINE SetParticlePosition(FractNbr,NbrOfParticle,mode)
!===================================================================================================================================
  USE MOD_part_MPI_Vars,      ONLY : PMPIVAR,PMPIInsert
!===================================================================================================================================
#else
SUBROUTINE SetParticlePosition(FractNbr,NbrOfParticle)                                             
#endif
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
  ! modules
  USE MOD_Particle_Vars
  USE MOD_PIC_Vars
  USE MOD_Timedisc_Vars,         ONLY : dt
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
   INTEGER                          :: FractNbr,NbrOfParticle                                      
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
   INTENT(IN)                       :: FractNbr                                                    
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
        (TRIM(Species(FractNbr)%SpaceIC).NE.'circle_equidistant'                 ) .AND. &
        (TRIM(Species(FractNbr)%SpaceIC).NE.'sin_deviation'                      ) .AND. &
        (TRIM(Species(FractNbr)%SpaceIC).NE.'cuboid_with_equidistant_distribution') .AND. &
        (TRIM(Species(FractNbr)%SpaceIC).NE.'line_with_equidistant_distribution' )) THEN
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

      SELECT CASE(TRIM(Species(FractNbr)%SpaceIC))
      CASE ('point')
         Particle_pos = Species(FractNbr)%BasePointIC
         DO i=1,chunkSize
            particle_positions(i*3-2) = Particle_pos(1)
            particle_positions(i*3-1) = Particle_pos(2)
            particle_positions(i*3  ) = Particle_pos(3)
         END DO
      CASE ('line_with_equidistant_distribution')
        IF(NbrOfParticle.EQ.1)THEN
           Particle_pos = Species(FractNbr)%BasePointIC + 0.5 * Species(FractNbr)%BaseVector1IC
        ELSE
           VectorGap = Species(FractNbr)%BaseVector1IC/(NbrOfParticle-1)
           DO i=1,chunkSize
              Particle_pos = Species(FractNbr)%BasePointIC + (i-1)*VectorGap
              particle_positions(i*3-2) = Particle_pos(1)
              particle_positions(i*3-1) = Particle_pos(2)
              particle_positions(i*3  ) = Particle_pos(3)
           END DO
        END IF
      CASE ('line')
        DO i=1,chunkSize
           CALL RANDOM_NUMBER(RandVal1)
           Particle_pos = Species(FractNbr)%BasePointIC + Species(FractNbr)%BaseVector1IC*RandVal1
           particle_positions(i*3-2) = Particle_pos(1)
           particle_positions(i*3-1) = Particle_pos(2)
           particle_positions(i*3  ) = Particle_pos(3)
        END DO
      CASE('disc')
        IF (Species(FractNbr)%NormalIC(3).NE.0) THEN
           lineVector(1) = 1.0
           lineVector(2) = 1.0
           lineVector(3) = -(Species(FractNbr)%NormalIC(1)+Species(FractNbr)%NormalIC(2))/ &
                             Species(FractNbr)%NormalIC(3)
        ELSE
           IF (Species(FractNbr)%NormalIC(2).NE.0) THEN
              lineVector(1) = 1.0
              lineVector(3) = 1.0
              lineVector(2) = -(Species(FractNbr)%NormalIC(1)+Species(FractNbr)%NormalIC(3))/ &
                                Species(FractNbr)%NormalIC(2)
           ELSE
              IF (Species(FractNbr)%NormalIC(1).NE.0) THEN
                 lineVector(2) = 1.0
                 lineVector(3) = 1.0
                 lineVector(1) = -(Species(FractNbr)%NormalIC(2)+Species(FractNbr)%NormalIC(3))/ &
                                   Species(FractNbr)%NormalIC(1)
              ELSE
                 WRITE(*,*) 'Error in SetParticlePosition, NormalIC Vektor darf nicht Nullvektor sein'
                 STOP
              END IF
           END IF
        END IF
        
        lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * &
             lineVector(2) + lineVector(3) * lineVector(3))
        
        lineVector2(1) = Species(FractNbr)%NormalIC(2) * lineVector(3) - &
             Species(FractNbr)%NormalIC(3) * lineVector(2)
        lineVector2(2) = Species(FractNbr)%NormalIC(3) * lineVector(1) - &
             Species(FractNbr)%NormalIC(1) * lineVector(3)
        lineVector2(3) = Species(FractNbr)%NormalIC(1) * lineVector(2) - &
             Species(FractNbr)%NormalIC(2) * lineVector(1)
        
        lineVector2 = lineVector2 / SQRT(lineVector2(1) * lineVector2(1) + lineVector2(2) * &
             lineVector2(2) + lineVector2(3) * lineVector2(3))

        DO i=1,chunkSize
           radius = Species(FractNbr)%RadiusIC + 1
           DO WHILE(radius.GT.Species(FractNbr)%RadiusIC)
              CALL RANDOM_NUMBER(RandVal)
              RandVal = RandVal * 2. - 1.
              Particle_pos = Species(FractNbr)%BasePointIC + Species(FractNbr)%RadiusIC * &
                       (RandVal(1) * lineVector + RandVal(2) *lineVector2)

              radius = SQRT( (Particle_pos(1)-Species(FractNbr)%BasePointIC(1)) * &
                             (Particle_pos(1)-Species(FractNbr)%BasePointIC(1)) + &
                             (Particle_pos(2)-Species(FractNbr)%BasePointIC(2)) * &
                             (Particle_pos(2)-Species(FractNbr)%BasePointIC(2)) + &
                             (Particle_pos(3)-Species(FractNbr)%BasePointIC(3)) * &
                             (Particle_pos(3)-Species(FractNbr)%BasePointIC(3)) )
           END DO
           particle_positions(i*3-2) = Particle_pos(1)
           particle_positions(i*3-1) = Particle_pos(2)
           particle_positions(i*3  ) = Particle_pos(3)
        END DO
      CASE('circle')
        IF (Species(FractNbr)%NormalIC(3).NE.0) THEN
           lineVector(1) = 1.0
           lineVector(2) = 1.0
           lineVector(3) = -(Species(FractNbr)%NormalIC(1)+Species(FractNbr)%NormalIC(2))/ &
                             Species(FractNbr)%NormalIC(3)
        ELSE
           IF (Species(FractNbr)%NormalIC(2).NE.0) THEN
              lineVector(1) = 1.0
              lineVector(3) = 1.0
              lineVector(2) = -(Species(FractNbr)%NormalIC(1)+Species(FractNbr)%NormalIC(3))/ &
                                Species(FractNbr)%NormalIC(2)
           ELSE
              IF (Species(FractNbr)%NormalIC(1).NE.0) THEN
                 lineVector(2) = 1.0
                 lineVector(3) = 1.0
                 lineVector(1) = -(Species(FractNbr)%NormalIC(2)+Species(FractNbr)%NormalIC(3))/ &
                                   Species(FractNbr)%NormalIC(1)
              ELSE
                 WRITE(*,*) 'Error in SetParticlePosition, NormalIC Vektor darf nicht Nullvektor sein'
                 STOP
              END IF
           END IF
        END IF
        
        lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * &
             lineVector(2) + lineVector(3) * lineVector(3))
        
        lineVector2(1) = Species(FractNbr)%NormalIC(2) * lineVector(3) - &
             Species(FractNbr)%NormalIC(3) * lineVector(2)
        lineVector2(2) = Species(FractNbr)%NormalIC(3) * lineVector(1) - &
             Species(FractNbr)%NormalIC(1) * lineVector(3)
        lineVector2(3) = Species(FractNbr)%NormalIC(1) * lineVector(2) - &
             Species(FractNbr)%NormalIC(2) * lineVector(1)
        
        lineVector2 = lineVector2 / SQRT(lineVector2(1) * lineVector2(1) + lineVector2(2) * &
             lineVector2(2) + lineVector2(3) * lineVector2(3))

        radius = Species(FractNbr)%RadiusIC
        DO i=1,chunkSize
           CALL RANDOM_NUMBER(RandVal1)
           argumentTheta = 2.*pi*RandVal1
           Particle_pos = Species(FractNbr)%BasePointIC +        &
                          linevector * cos(argumentTheta) * radius +  &
                          linevector2 * sin(argumentTheta) * radius
           particle_positions(i*3-2) = Particle_pos(1)
           particle_positions(i*3-1) = Particle_pos(2)
           particle_positions(i*3  ) = Particle_pos(3)
        END DO
      CASE('gyrotron_circle')
        IF (Species(FractNbr)%NormalIC(3).NE.0) THEN
           lineVector(1) = 1.0
           lineVector(2) = 1.0
           lineVector(3) = -(Species(FractNbr)%NormalIC(1)+Species(FractNbr)%NormalIC(2))/ &
                             Species(FractNbr)%NormalIC(3)
        ELSE
           IF (Species(FractNbr)%NormalIC(2).NE.0) THEN
              lineVector(1) = 1.0
              lineVector(3) = 1.0
              lineVector(2) = -(Species(FractNbr)%NormalIC(1)+Species(FractNbr)%NormalIC(3))/ &
                                Species(FractNbr)%NormalIC(2)
           ELSE
              IF (Species(FractNbr)%NormalIC(1).NE.0) THEN
                 lineVector(2) = 1.0
                 lineVector(3) = 1.0
                 lineVector(1) = -(Species(FractNbr)%NormalIC(2)+Species(FractNbr)%NormalIC(3))/ &
                                   Species(FractNbr)%NormalIC(1)
              ELSE
                 WRITE(*,*) 'Error in SetParticlePosition, NormalIC Vektor darf nicht Nullvektor sein'
                 STOP
              END IF
           END IF
        END IF
        
        lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * &
             lineVector(2) + lineVector(3) * lineVector(3))
        
        lineVector2(1) = Species(FractNbr)%NormalIC(2) * lineVector(3) - &
             Species(FractNbr)%NormalIC(3) * lineVector(2)
        lineVector2(2) = Species(FractNbr)%NormalIC(3) * lineVector(1) - &
             Species(FractNbr)%NormalIC(1) * lineVector(3)
        lineVector2(3) = Species(FractNbr)%NormalIC(1) * lineVector(2) - &
             Species(FractNbr)%NormalIC(2) * lineVector(1)
        
        lineVector2 = lineVector2 / SQRT(lineVector2(1) * lineVector2(1) + lineVector2(2) * &
             lineVector2(2) + lineVector2(3) * lineVector2(3))

        radius = Species(FractNbr)%RadiusIC
        DO i=1,chunkSize
           CALL RANDOM_NUMBER(RandVal1)
           argumentTheta = 2.*pi*RandVal1
           Particle_pos = Species(FractNbr)%BasePointIC +        &
                          linevector * cos(argumentTheta) * radius +  &
                          linevector2 * sin(argumentTheta) * radius
           ! Change position of particle on the small gyro circle
           ! take normal vecotr of the circle
           n(1:3) = Species(FractNbr)%NormalIC(1:3)
           ! generate radius vector (later it will be multiplied by the length of the
           ! gyro circles. For now we just need the vector)
           radius_vec(1) = Particle_pos(1) - Species(1)%BasePointIC(1)
           radius_vec(2) = Particle_pos(2) - Species(1)%BasePointIC(2)
           radius_vec(3) = Particle_pos(3) - Species(1)%BasePointIC(3)
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
           IF (NbrOfParticle.EQ.Species(FractNbr)%initialParticleNumber) THEN
             particle_positions(i*3  ) = Species(FractNbr)%BasePointIC(3) &
                                             + RandVal1 * Species(FractNbr)%CuboidHeightIC
           ELSE
             particle_positions(i*3  ) = Species(FractNbr)%BasePointIC(3) &
                                             + RandVal1 * dt & 
                                             * Species(FractNbr)%VeloIC/Species(FractNbr)%alpha 
           END IF

!          2. calculate curved B-field at z-position in order to determin size of gyroradius
           IF (useCurvedExternalField) THEN
              Bintpol = InterpolateCurvedExternalField(particle_positions(i*3))
              rgyrate = 1./ SQRT ( 1 - (Species(FractNbr)%VeloIC**2 * (1 + 1./Species(FractNbr)%alpha**2)) * c_inv * c_inv ) * &
                                  Species(FractNbr)%MassIC * Species(fractNbr)%VeloIC / &
                        ( Bintpol * abs( Species(FractNbr)%ChargeIC) )
           ELSE
             rgyrate =  Species(FractNbr)%RadiusICGyro
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
        IF (Species(FractNbr)%NormalIC(3).NE.0) THEN
           lineVector(1) = 1.0
           lineVector(2) = 1.0
           lineVector(3) = -(Species(FractNbr)%NormalIC(1)+Species(FractNbr)%NormalIC(2))/ &
                             Species(FractNbr)%NormalIC(3)
        ELSE
           IF (Species(FractNbr)%NormalIC(2).NE.0) THEN
              lineVector(1) = 1.0
              lineVector(3) = 1.0
              lineVector(2) = -(Species(FractNbr)%NormalIC(1)+Species(FractNbr)%NormalIC(3))/ &
                                Species(FractNbr)%NormalIC(2)
           ELSE
              IF (Species(FractNbr)%NormalIC(1).NE.0) THEN
                 lineVector(2) = 1.0
                 lineVector(3) = 1.0
                 lineVector(1) = -(Species(FractNbr)%NormalIC(2)+Species(FractNbr)%NormalIC(3))/ &
                                   Species(FractNbr)%NormalIC(1)
              ELSE
                 WRITE(*,*) 'Error in SetParticlePosition, NormalIC Vektor darf nicht Nullvektor sein'
                 STOP
              END IF
           END IF
        END IF
                   
        lineVector2(1) = Species(FractNbr)%NormalIC(2) * lineVector(3) - &
             Species(FractNbr)%NormalIC(3) * lineVector(2)
        lineVector2(2) = Species(FractNbr)%NormalIC(3) * lineVector(1) - &
             Species(FractNbr)%NormalIC(1) * lineVector(3)
        lineVector2(3) = Species(FractNbr)%NormalIC(1) * lineVector(2) - &
             Species(FractNbr)%NormalIC(2) * lineVector(1)

        lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * &
             lineVector(2) + lineVector(3) * lineVector(3))

        lineVector2 = lineVector2 / SQRT(lineVector2(1) * lineVector2(1) + lineVector2(2) * &
             lineVector2(2) + lineVector2(3) * lineVector2(3))

        radius = Species(FractNbr)%RadiusIC
        DO i=1,chunkSize
           argumentTheta = 2.*pi*i/chunkSize
           Particle_pos = Species(FractNbr)%BasePointIC +        &
                          linevector * cos(argumentTheta) * radius +  &
                          linevector2 * sin(argumentTheta) * radius
           particle_positions(i*3-2) = Particle_pos(1)
           particle_positions(i*3-1) = Particle_pos(2)
           particle_positions(i*3  ) = Particle_pos(3)
        END DO
      CASE('cuboid')
        DO i=1,chunkSize          
           CALL RANDOM_NUMBER(RandVal)
           Particle_pos = Species(FractNbr)%BasePointIC + Species(FractNbr)%BaseVector1IC * RandVal(1)
           Particle_pos = Particle_pos + Species(FractNbr)%BaseVector2IC * RandVal(2)
           lineVector(1) = Species(FractNbr)%BaseVector1IC(2) * Species(FractNbr)%BaseVector2IC(3) - &
                Species(FractNbr)%BaseVector1IC(3) * Species(FractNbr)%BaseVector2IC(2)
           lineVector(2) = Species(FractNbr)%BaseVector1IC(3) * Species(FractNbr)%BaseVector2IC(1) - &
                Species(FractNbr)%BaseVector1IC(1) * Species(FractNbr)%BaseVector2IC(3)
           lineVector(3) = Species(FractNbr)%BaseVector1IC(1) * Species(FractNbr)%BaseVector2IC(2) - &
                Species(FractNbr)%BaseVector1IC(2) * Species(FractNbr)%BaseVector2IC(1)
           IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
           ELSE
              lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + &
                   lineVector(3) * lineVector(3))
           END IF          
#if (PP_TimeDiscMethod==201)
           Particle_pos = Particle_pos + lineVector * Species(FractNbr)%CuboidHeightIC * dt / dt_maxwell * RandVal(3) 
           !scaling due to variable time step
#else
           Particle_pos = Particle_pos + lineVector * Species(FractNbr)%CuboidHeightIC * RandVal(3) 
#endif
           particle_positions(i*3-2) = Particle_pos(1)
           particle_positions(i*3-1) = Particle_pos(2)
           particle_positions(i*3  ) = Particle_pos(3)
        END DO
      CASE('cylinder')
        DO i=1,chunkSize
           radius = Species(FractNbr)%RadiusIC + 1
           DO WHILE((radius.GT.Species(FractNbr)%RadiusIC) .OR.(radius.LT.Species(FractNbr)%Radius2IC))          
              CALL RANDOM_NUMBER(RandVal)
              Particle_pos = Species(FractNbr)%BasePointIC + Species(FractNbr)%BaseVector1IC * (RandVal(1)*2-1)
              Particle_pos = Particle_pos + Species(FractNbr)%BaseVector2IC * (RandVal(2)*2-1)

              radius = SQRT( (Particle_pos(1)-Species(FractNbr)%BasePointIC(1)) * &
                             (Particle_pos(1)-Species(FractNbr)%BasePointIC(1)) + &
                             (Particle_pos(2)-Species(FractNbr)%BasePointIC(2)) * &
                             (Particle_pos(2)-Species(FractNbr)%BasePointIC(2)) + &
                             (Particle_pos(3)-Species(FractNbr)%BasePointIC(3)) * &
                             (Particle_pos(3)-Species(FractNbr)%BasePointIC(3)) )
           END DO
           lineVector(1) = Species(FractNbr)%BaseVector1IC(2) * Species(FractNbr)%BaseVector2IC(3) - &
                Species(FractNbr)%BaseVector1IC(3) * Species(FractNbr)%BaseVector2IC(2)
           lineVector(2) = Species(FractNbr)%BaseVector1IC(3) * Species(FractNbr)%BaseVector2IC(1) - &
                Species(FractNbr)%BaseVector1IC(1) * Species(FractNbr)%BaseVector2IC(3)
           lineVector(3) = Species(FractNbr)%BaseVector1IC(1) * Species(FractNbr)%BaseVector2IC(2) - &
                Species(FractNbr)%BaseVector1IC(2) * Species(FractNbr)%BaseVector2IC(1)
           IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
           ELSE
              lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + &
                   lineVector(3) * lineVector(3))
           END IF          
#if (PP_TimeDiscMethod==201)
           Particle_pos = Particle_pos + lineVector * Species(FractNbr)%CylinderHeightIC * dt / dt_maxwell * RandVal(3) 
           !scaling due to variable time step           
#else
           Particle_pos = Particle_pos + lineVector * Species(FractNbr)%CylinderHeightIC * RandVal(3)
#endif
           particle_positions(i*3-2) = Particle_pos(1)
           particle_positions(i*3-1) = Particle_pos(2)
           particle_positions(i*3  ) = Particle_pos(3)
        END DO
      CASE('LD_insert')
        CALL LD_SetParticlePosition(chunkSize,particle_positions_Temp)
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
            particle_positions(i*3-2 : i*3) = PartState(ParticleIndexNbr-Species(1)%initialParticleNumber,1:3)
         END DO
#endif
      CASE ('cuboid_with_equidistant_distribution') 
         IF(Species(FractNbr)%initialParticleNumber.NE. &
              (Species(FractNbr)%maxParticleNumberX*Species(FractNbr)%maxParticleNumberY*Species(FractNbr)%maxParticleNumberZ)) THEN
           WRITE(*,*) 'ERROR: Number of ini particles of species ',FractNbr,' does not match number of particles in each direction!'
           STOP
         END IF
         xlen = SQRT(Species(FractNbr)%BaseVector1IC(1)**2 &
              + Species(FractNbr)%BaseVector1IC(2)**2 &
              + Species(FractNbr)%BaseVector1IC(3)**2 )
         ylen = SQRT(Species(FractNbr)%BaseVector2IC(1)**2 &
              + Species(FractNbr)%BaseVector2IC(2)**2 &
              + Species(FractNbr)%BaseVector2IC(3)**2 )
         zlen = ABS(Species(FractNbr)%CuboidHeightIC)

         ! make sure the vectors correspond to x,y,z-dir
         IF ((xlen.NE.Species(FractNbr)%BaseVector1IC(1)).OR. &
            (ylen.NE.Species(FractNbr)%BaseVector2IC(2)).OR. &
            (zlen.NE.Species(FractNbr)%CuboidHeightIC)) THEN
            WRITE(*,*) 'Basevectors1IC,-2IC and CuboidHeightIC have to be in x,y,z-direction, respectively for emission condition'
            STOP
          END IF
         x_step = xlen/Species(FractNbr)%maxParticleNumberX
         y_step = ylen/Species(FractNbr)%maxParticleNumberY
         z_step = zlen/Species(FractNbr)%maxParticleNumberZ
         iPart = 1
         DO i=1,Species(FractNbr)%maxParticleNumberX
           x_pos = (i-0.5) * x_step + Species(FractNbr)%BasePointIC(1)
           DO j=1,Species(FractNbr)%maxParticleNumberY
             y_pos =  Species(FractNbr)%BasePointIC(2) + (j-0.5) * y_step
             DO k=1,Species(FractNbr)%maxParticleNumberZ
               particle_positions(iPart*3-2) = x_pos
               particle_positions(iPart*3-1) = y_pos
               particle_positions(iPart*3  ) = Species(FractNbr)%BasePointIC(3) &
                    + (k-0.5) * z_step
               iPart = iPart + 1
             END DO
           END DO
         END DO
      CASE('sin_deviation')
         IF(Species(FractNbr)%initialParticleNumber.NE. &
              (Species(FractNbr)%maxParticleNumberX*Species(FractNbr)%maxParticleNumberY*Species(FractNbr)%maxParticleNumberZ)) THEN
           WRITE(*,*) 'ERROR: Number of ini particles of species ',FractNbr,' does not match number of particles in each direction!'
           STOP
         END IF
         xlen = abs(GEO%xmaxglob  - GEO%xminglob)  
         ylen = abs(GEO%ymaxglob  - GEO%yminglob)
         zlen = abs(GEO%zmaxglob  - GEO%zminglob)
         x_step = xlen/Species(FractNbr)%maxParticleNumberX
         y_step = ylen/Species(FractNbr)%maxParticleNumberY
         z_step = zlen/Species(FractNbr)%maxParticleNumberZ
         iPart = 1
         DO i=1,Species(FractNbr)%maxParticleNumberX
            x_pos = (i * x_step - x_step*0.5)
            x_pos = GEO%xminglob + x_pos + Species(FractNbr)%Amplitude &
                    * sin(Species(FractNbr)%WaveNumber * x_pos)
            DO j=1,Species(FractNbr)%maxParticleNumberY
              y_pos =  GEO%yminglob + j * y_step - y_step * 0.5
              DO k=1,Species(FractNbr)%maxParticleNumberZ
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

SUBROUTINE SetParticleVelocity(FractNbr,NbrOfParticle)                                             !
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
   INTEGER                          :: FractNbr,NbrOfParticle                                      !
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
   INTENT(IN)                       :: FractNbr                                                    !
   INTENT(INOUT)                    :: NbrOfParticle                                               !
!----------------------------------------------------------------------------------------------------------------------------------
   IF(NbrOfParticle.lt.1) RETURN
 
!   IF(NbrOfParticle.gt.PIC%maxParticleNumber)THEN
!      NbrOfParticle = PIC%maxParticleNumber
!   END IF

   SELECT CASE(Species(FractNbr)%velocityDistribution)
   CASE('random')
     i = 1
     DO WHILE (i .le. NbrOfParticle)
       PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
       IF (PositionNbr .ne. 0) THEN
         CALL RANDOM_NUMBER(RandVal)
         RandVal(:) = RandVal(:) - 0.5
         RandVal(:) = RandVal(:)/SQRT(RandVal(1)**2+RandVal(2)**2+RandVal(3)**2)
         PartState(PositionNbr,4:6) = RandVal(1:3) * Species(FractNbr)%VeloIC
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
            PartState(PositionNbr,4:6) = Species(FractNbr)%VeloVecIC(1:3) * Species(FractNbr)%VeloIC
         END IF
         i = i + 1
      END DO
   CASE('radial_constant')
      i = 1
      DO WHILE (i .le. NbrOfParticle)
         PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
         IF (PositionNbr .ne. 0) THEN
            Radius(1:3) = PartState(PositionNbr,1:3) - Species(FractNbr)%BasePointIC(1:3)
            !  Unity radius
            !Radius(1:3) = Radius(1:3) / Species(FractNbr)%RadiusIC
            Radius(1:3) = Radius(1:3) / SQRT(Radius(1)**2+Radius(2)**2+Radius(3)**2) 
            PartState(PositionNbr,4:6) = Radius(1:3) * Species(FractNbr)%VeloIC
         END IF
         i = i + 1
         !print*, SQRT(PartState(PositionNbr,4)**2+PartState(PositionNbr,5)**2+PartState(PositionNbr,6)**2)
      END DO
   CASE('tangential_constant')
      i = 1
      DO WHILE (i .le. NbrOfParticle)
         PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
         IF (PositionNbr .ne. 0) THEN
            Radius(1:3) = PartState(PositionNbr,1:3) - Species(FractNbr)%BasePointIC(1:3)
            !  Normal Vector of circle
            n_vec(1:3) = Species(FractNbr)%NormalIC(1:3)
            ! If we're doing Leapfrog, then use velocities from half-timestep before
            IF (ParticlePushMethod.EQ.'boris_leap_frog_scheme') THEN
              Angle = 0.5 * dt * Species(FractNbr)%VeloIC / Species(FractNbr)%RadiusIC ! 0.5*dt*(v/r)
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
            IF (Species(FractNbr)%alpha .gt. 0.) THEN 
              n_vec = n_vec * ( 1 / Species(FractNbr)%alpha )
            ELSE 
              n_vec = 0
            END IF
            !  And finally the velocities
            PartState(PositionNbr,4:6) = (tan_vec(1:3) + n_vec(1:3)) * Species(FractNbr)%VeloIC
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
           r1 = Species(1)%RadiusIC
           r2 = Species(1)%RadiusICGyro
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
            n_vec(1:3) = Species(FractNbr)%NormalIC(1:3)

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
            IF (Species(FractNbr)%alpha .gt. 0.) THEN 
              n_vec = n_vec * ( 1. / Species(FractNbr)%alpha )
            ELSE 
              n_vec = 0.
            END IF
            !  And finally the velocities
            PartState(PositionNbr,4:6) = (tan_vec(1:3) + n_vec(1:3)) * Species(FractNbr)%VeloIC
            IF (ABS(SQRT(PartState(PositionNbr,4)*PartState(PositionNbr,4) &
                       + PartState(PositionNbr,5)*PartState(PositionNbr,5))&
                       - Species(FractNbr)%VeloIC) .GT. 10.) THEN
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
      maxwellfac = SQRT(3. * 1.380650524E-23 * Species(FractNbr)%MWTemperatureIC / &              ! velocity of maximum
                     (Species(FractNbr)%MassIC*v2_sum))

      i = 1
      DO WHILE (i .le. NbrOfParticle)
         PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
         IF (PositionNbr .ne. 0) THEN
           PartState(PositionNbr,4:6) = (PartState(PositionNbr,4:6) - v_sum(1:3)) * maxwellfac &
                                        + Species(FractNbr)%VeloIC *Species(FractNbr)%VeloVecIC(1:3)        
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

END MODULE MOD_part_emission
