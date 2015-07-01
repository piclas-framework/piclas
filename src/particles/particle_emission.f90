#include "boltzplatz.h"

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

INTERFACE SetParticleChargeAndMass
  MODULE PROCEDURE SetParticleChargeAndMass
END INTERFACE

INTERFACE SetParticleVelocity
  MODULE PROCEDURE SetParticleVelocity
END INTERFACE

INTERFACE SetParticleMPF
  MODULE PROCEDURE SetParticleMPF
END INTERFACE

!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC         :: InitializeParticleEmission, ParticleInserting &
                , SetParticleChargeAndMass, SetParticleVelocity, SetParticleMPF                              
!===================================================================================================================================
                                                                                                  
CONTAINS                                                                                           
                                                                                                   
SUBROUTINE InitializeParticleEmission()
!===================================================================================================================================
! Initialize particles / Insert initial particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars,  ONLY : Species,nSpecies,PDM,PEM, usevMPF
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
INTEGER               :: nPartInside
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
          !CALL abort(__STAMP__,&
          !      ' particle pressure not moved to picasso!')
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
!          CALL abort(__STAMP__,&
!                ' particle pressure not moved in picasso!')
          CALL ParticleInsideCheck(i, iInit, nPartInside, TempInside, EInside)
          IF (Species(i)%Init(iInit)%ParticleEmission .GT. nPartInside) THEN
            NbrOfParticle = INT(Species(i)%Init(iInit)%ParticleEmission) - nPartInside
            IPWRITE(UNIT_stdOut,*) 'Emission PartNum (Spec ',i,')', NbrOfParticle
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

SWRITE(UNIT_stdOut,'(A)') ' ...DONE '

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
USE MOD_Globals
USE MOD_Timedisc_Vars         , ONLY : dt
#if (PP_TimeDiscMethod==1) ||  (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6)
USE MOD_Timedisc_Vars         , ONLY : iter
#endif
USE MOD_Particle_Vars
USE MOD_PIC_Vars
USE MOD_part_tools             ,ONLY : UpdateNextFreePosition  
USE MOD_DSMC_Vars              ,ONLY : useDSMC, CollisMode
USE MOD_DSMC_Init              ,ONLY : DSMC_SetInternalEnr_LauxVFD
!#if (PP_TimeDiscMethod==1000) || (PP_TimeDiscMethod==1001)
!USE MOD_LD_Init                ,ONLY : CalcDegreeOfFreedom
!USE MOD_LD_Vars
!#endif
USE MOD_Particle_Analyze_Vars  ,ONLY: CalcPartBalance,nPartIn,PartEkinIn
#if (PP_TimeDiscMethod==1) ||  (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6)
USE MOD_Particle_Analyze_Vars  ,ONLY: nPartInTmp,PartEkinInTmp,PartAnalyzeStep
#endif
USE MOD_Particle_Analyze       ,ONLY: CalcEkinPart
USE MOD_part_pressure          ,ONLY: ParticlePressure, ParticlePressureRem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
#ifdef MPI
INTEGER, OPTIONAL                :: mode_opt
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Local variable declaration                                                                       
INTEGER                          :: i , iPart, PositionNbr, iInit                                                          
INTEGER                , SAVE    :: NbrOfParticle=0                                             
INTEGER(KIND=8)                 :: inserted_Particle_iter,inserted_Particle_time               
INTEGER(KIND=8)                 :: inserted_Particle_diff  
#ifdef MPI
INTEGER                          :: mode                                            
#endif
!===================================================================================================================================

!----------------------------------------------------------------------------------------------------------------------------------
!!! VORSICHT: FUNKTIONIERT SO MOMENTAN NUR MIT 1 SPEZIES!!!!
! --- fuer mehr als eine Spezies gibt es bei der Benutzung des
!     mode_opt Flags Probleme mit den non-blocking communications.
!     Es koennte dann passieren, dass Nachrichten falsch zugeordnet werden.
!     Sicherheitshalber sollte man kein mode_opt Argument bei mehrern
!     Spezies uebergeben.
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
          IF (Species(i)%Init(iInit)%VirtPreInsert .AND. Species(i)%Init(iInit)%PartDensity.GT.0.) THEN
            NbrOfParticle = 0 ! calculated within SetParticlePosition itself!
          ELSE
            inserted_Particle_iter = INT(Species(i)%Init(iInit)%ParticleEmission * dt,8)
            inserted_Particle_time = INT(Species(i)%Init(iInit)%ParticleEmission * (Time + dt),8)
            inserted_Particle_diff = INT(inserted_Particle_time - Species(i)%Init(iInit)%InsertedParticle &
              - inserted_Particle_iter,8)
            NbrOfParticle = INT(inserted_Particle_iter + inserted_Particle_diff,4)
            ! if maxwell velo dist and less than 5 parts: skip (to ensure maxwell dist)
            IF (TRIM(Species(i)%Init(iInit)%velocityDistribution).EQ.'maxwell') THEN
              IF (NbrOfParticle.LT.5) NbrOfParticle=0
            END IF
            Species(i)%Init(iInit)%InsertedParticle = Species(i)%Init(iInit)%InsertedParticle + INT(NbrOfParticle,8)
          END IF
        CASE(2)    ! Emission Type: Particles per Iteration
          NbrOfParticle = INT(Species(i)%Init(iInit)%ParticleEmission)
        CASE(3)
!          CALL abort(__STAMP__,&
!                ' particle pressure not moved in picasso!')
          CALL ParticlePressure (i, iInit, NbrOfParticle)
          ! if maxwell velo dist and less than 5 parts: skip (to ensure maxwell dist)
          IF (TRIM(Species(i)%Init(iInit)%velocityDistribution).EQ.'maxwell') THEN
            IF (NbrOfParticle.LT.5) NbrOfParticle=0
          END IF
        CASE(5) ! removal of all parts in pressure area and re-insertion
         ! CALL abort(__STAMP__,&
         !       ' particle pressure not moved in picasso!')
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
       IF (useDSMC.AND.(CollisMode.GT.1)) THEN
         iPart = 1
         DO WHILE (iPart .le. NbrOfParticle)
           PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
           IF (PositionNbr .ne. 0) THEN
             CALL DSMC_SetInternalEnr_LauxVFD(i, iInit, PositionNbr)
           END IF
           iPart = iPart + 1
         END DO
       END IF
!#if (PP_TimeDiscMethod==1000) || (PP_TimeDiscMethod==1001)
!       iPart = 1
!       DO WHILE (iPart .le. NbrOfParticle)
!         PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
!         IF (PositionNbr .ne. 0) THEN
!           PartStateBulkValues(PositionNbr,1) = Species(i)%Init(iInit)%VeloVecIC(1) * Species(i)%Init(iInit)%VeloIC
!           PartStateBulkValues(PositionNbr,2) = Species(i)%Init(iInit)%VeloVecIC(2) * Species(i)%Init(iInit)%VeloIC
!           PartStateBulkValues(PositionNbr,3) = Species(i)%Init(iInit)%VeloVecIC(3) * Species(i)%Init(iInit)%VeloIC
!           PartStateBulkValues(PositionNbr,4) = Species(i)%Init(iInit)%MWTemperatureIC
!           PartStateBulkValues(PositionNbr,5) = CalcDegreeOfFreedom(PositionNbr)
!         END IF
!         iPart = iPart + 1
!       END DO
!#endif
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
        !  CALL abort(__STAMP__,&
        !        ' particle pressure not moved in picasso!')
        CALL ParticleInsertingPressureOut(i,iInit,NbrofParticle)
      END IF
      CALL SetParticleChargeAndMass(i,NbrOfParticle)
      IF (usevMPF) CALL SetParticleMPF(i,NbrOfParticle)
      ! define molecule stuff
      IF (useDSMC.AND.(CollisMode.GT.1)) THEN
        iPart = 1
        DO WHILE (iPart .le. NbrOfParticle)
          PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
          IF (PositionNbr .ne. 0) THEN
            CALL DSMC_SetInternalEnr_LauxVFD(i, iInit, PositionNbr)
          END IF
          iPart = iPart + 1
        END DO
      END IF
!#if (PP_TimeDiscMethod==1000) || (PP_TimeDiscMethod==1001) !      iPart = 1 !      DO WHILE (iPart .le. NbrOfParticle)
!        PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
!        IF (PositionNbr .ne. 0) THEN
!          PartStateBulkValues(PositionNbr,1) = Species(i)%Init(iInit)%VeloVecIC(1) * Species(i)%Init(iInit)%VeloIC
!          PartStateBulkValues(PositionNbr,2) = Species(i)%Init(iInit)%VeloVecIC(2) * Species(i)%Init(iInit)%VeloIC
!          PartStateBulkValues(PositionNbr,3) = Species(i)%Init(iInit)%VeloVecIC(3) * Species(i)%Init(iInit)%VeloIC
!          PartStateBulkValues(PositionNbr,4) = Species(i)%Init(iInit)%MWTemperatureIC
!          PartStateBulkValues(PositionNbr,5) = CalcDegreeOfFreedom(PositionNbr)
!        END IF
!        iPart = iPart + 1
!      END DO
!#endif
      ! instead of UpdateNextfreePosition we update the
      ! particleVecLength only.
      PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + NbrOfParticle
      PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfParticle
      !CALL UpdateNextFreePosition()
    END IF
    ! compute number of input particles and energy
    IF(CalcPartBalance) THEN
#if (PP_TimeDiscMethod==1) ||  (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6)
      IF((MOD(iter+1,PartAnalyzeStep).EQ.0).AND.(iter.GT.0))THEN ! caution if correct
        nPartInTmp(i)=nPartInTmp(i) + NBrofParticle
        DO iPart=1,NbrOfparticle
          PartEkinInTmp(PartSpecies(iPart))=PartEkinInTmp(PartSpecies(iPart))+CalcEkinPart(iPart)
        END DO ! iPart
      ELSE
        nPartIn(i)=nPartIn(i) + NBrofParticle
        DO iPart=1,NbrOfparticle
          PartEkinIn(PartSpecies(iPart))=PartEkinIn(PartSpecies(iPart))+CalcEkinPart(iPart)
        END DO ! iPart
      END IF
#else
      nPartIn(i)=nPartIn(i) + NBrofParticle
      DO iPart=1,NbrOfparticle
        PartEkinIn(PartSpecies(iPart))=PartEkinIn(PartSpecies(iPart))+CalcEkinPart(iPart)
      END DO ! iPart
#endif
    END IF ! CalcPartBalance
  END DO
END DO

END SUBROUTINE ParticleInserting
                                                                                                   
#ifdef MPI
SUBROUTINE SetParticlePosition(FractNbr,iInit,NbrOfParticle,mode)
#else
SUBROUTINE SetParticlePosition(FractNbr,iInit,NbrOfParticle)                                             
#endif /* MPI*/
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
#ifdef MPI
USE MOD_Particle_MPI_Vars,     ONLY:PartMPI,PartMPIInsert
#endif /* MPI*/
USE MOD_Globals
USE MOD_PIC_Vars
USE MOD_Particle_Vars,         ONLY:Species,BoltzmannConst,PDM,PartState,OutputVpiWarnings
USE MOD_Particle_Mesh_Vars,    ONLY:GEO
USE MOD_Globals_Vars,          ONLY:PI, TwoepsMach
USE MOD_Timedisc_Vars,         ONLY:dt, iter, IterDisplayStep, DoDisplayIter
USE MOD_Particle_Mesh,         ONLY:SingleParticleToExactElement,SingleParticleToExactElementNoMap
USE MOD_Particle_Surfaces_vars,ONLY:DoRefMapping
USE MOD_PICInterpolation,      ONLY:InterpolateCurvedExternalField
USE MOD_PICInterpolation_vars, ONLY:useCurvedExternalField
USE MOD_Equation_vars,         ONLY:c_inv
!USE MOD_LD,                    ONLY : LD_SetParticlePosition
!#ifdef MPI
!! PilleO: to change into use MPi_2003 or so
!INCLUDE 'mpif.h'                                                                               
!#endif
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                       :: FractNbr, iInit                                                    
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                    :: NbrOfParticle                                               
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#ifdef MPI
INTEGER                                  :: mode
INTEGER                                  :: iProc,tProc, CellX, CellY, CellZ                                  
INTEGER                                  :: msg_status(1:MPI_STATUS_SIZE)                               
INTEGER                                  :: MessageSize
LOGICAL                                  :: InsideMyBGM                                                 
#endif
REAL,POINTER                             :: particle_positions(:)=>NULL()
INTEGER                                  :: allocStat                                           
INTEGER                                  :: i,j,k,ParticleIndexNbr                                      
INTEGER                                  :: mySumOfMatchedParticles, sumOfMatchedParticles              
INTEGER                                  :: nChunks, chunkSize, chunkSize2                                        
REAL                                     :: lineVector(3),VectorGap(3)         
REAL                                     :: RandVal(3), Particle_pos(3),lineVector2(3)                  
REAL                                     :: n(3) , radius_vec(3)                                        
REAL                                     :: II(3,3),JJ(3,3),NN(3,3)                                     
REAL                                     :: RandVal1                              
REAL                                     :: radius, argumentTheta                    
REAL                                     :: rgyrate, Bintpol, pilen                    
REAL                                     :: x_step, y_step, z_step,  x_pos , y_pos                      
REAL                                     :: xlen, ylen, zlen                                            
INTEGER                                  :: iPart                                               
REAL,ALLOCATABLE                         :: particle_positions_Temp(:) 
REAL                                     :: Vec3D(3), l_ins, v_line, delta_l, v_drift_line, A_ins, PartIns
REAL                                     :: v_drift_BV(2), lrel_ins_BV(4), BV_lengths(2), v_BV(2), delta_lBV(2)
REAL                                     :: intersecPoint(3), orifice_delta, lPeri, ParaCheck(3)
INTEGER                                  :: DimSend, orificePeriodic
LOGICAL                                  :: orificePeriodicLog(2), insideExcludeRegion
#ifdef MPI
REAL                                     :: PartInsRest
INTEGER                                  :: InitGroup
#endif                        
!===================================================================================================================================

#ifdef MPI
InitGroup=Species(FractNbr)%Init(iInit)%InitCOMM
IF(PartMPI%InitGroup(InitGroup)%COMM.EQ.MPI_COMM_NULL) THEN
  NbrofParticle=0
  RETURN
END IF
#endif /*MPI*/

PartIns=0.
lineVector = 0.0
A_ins = 0.0
l_ins = 0.0
lrel_ins_BV = 0.0
BV_lengths = 0.0
Particle_pos = 0.0
orificePeriodic = 0
orificePeriodicLog = .FALSE.
IF(Species(FractNbr)%Init(iInit)%VirtPreInsert) THEN ! (SpaceIC.EQ.'cuboid_vpi').OR.(SpaceIC.EQ.'cylinder_vpi')
  DimSend=6 !save (and send) velocities and positions
  ! the following is here, and not inside the select case, as it is needed to be excecuted just once and, most importantly, it
  ! calculates the virt. insertion height defining the virt. NbrOfParticle (which is chunked next) when PartDensity is used
  ! (-> could also be moved to particle_init?)
  lineVector(1) = Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * Species(FractNbr)%Init(iInit)%BaseVector2IC(3) - &
    Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * Species(FractNbr)%Init(iInit)%BaseVector2IC(2)
  lineVector(2) = Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * Species(FractNbr)%Init(iInit)%BaseVector2IC(1) - &
    Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * Species(FractNbr)%Init(iInit)%BaseVector2IC(3)
  lineVector(3) = Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * Species(FractNbr)%Init(iInit)%BaseVector2IC(2) - &
    Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * Species(FractNbr)%Init(iInit)%BaseVector2IC(1)
  IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
    CALL abort(__STAMP__,&
                'BaseVectors are parallel!')
  ELSE
    A_ins = SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + lineVector(3) * lineVector(3))
    lineVector = lineVector / A_ins
  END IF
  v_drift_line = Species(FractNbr)%Init(iInit)%VeloIC * &
    ( Species(FractNbr)%Init(iInit)%VeloVecIC(1)*lineVector(1) + Species(FractNbr)%Init(iInit)%VeloVecIC(2)*lineVector(2) &
    + Species(FractNbr)%Init(iInit)%VeloVecIC(3)*lineVector(3) ) !lineVector component of drift-velocity
  l_ins=dt * ( v_drift_line + Species(FractNbr)%Init(iInit)%NSigma & !virt. insertion height
    * SQRT(BoltzmannConst*Species(FractNbr)%Init(iInit)%MWTemperatureIC/Species(FractNbr)%MassIC) )
  IF( (TRIM(Species(FractNbr)%Init(iInit)%vpiDomainType).EQ.'freestream') .OR. &
      (TRIM(Species(FractNbr)%Init(iInit)%vpiDomainType).EQ.'orifice') ) THEN
    BV_lengths(1) = SQRT(Species(FractNbr)%Init(iInit)%BaseVector1IC(1)**2 + &
                         Species(FractNbr)%Init(iInit)%BaseVector1IC(2)**2 + &
                         Species(FractNbr)%Init(iInit)%BaseVector1IC(3)**2)
    v_drift_BV(1) = Species(FractNbr)%Init(iInit)%VeloIC / BV_lengths(1) * &
      ( Species(FractNbr)%Init(iInit)%VeloVecIC(1)*Species(FractNbr)%Init(iInit)%BaseVector1IC(1) &
      + Species(FractNbr)%Init(iInit)%VeloVecIC(2)*Species(FractNbr)%Init(iInit)%BaseVector1IC(2) &
      + Species(FractNbr)%Init(iInit)%VeloVecIC(3)*Species(FractNbr)%Init(iInit)%BaseVector1IC(3) ) !BV1 component of drift-velocity
    BV_lengths(2) = SQRT(Species(FractNbr)%Init(iInit)%BaseVector2IC(1)**2 + &
                         Species(FractNbr)%Init(iInit)%BaseVector2IC(2)**2 + &
                         Species(FractNbr)%Init(iInit)%BaseVector2IC(3)**2)
    v_drift_BV(2) = Species(FractNbr)%Init(iInit)%VeloIC / BV_lengths(2) * &
      ( Species(FractNbr)%Init(iInit)%VeloVecIC(1)*Species(FractNbr)%Init(iInit)%BaseVector2IC(1) &
      + Species(FractNbr)%Init(iInit)%VeloVecIC(2)*Species(FractNbr)%Init(iInit)%BaseVector2IC(2) &
      + Species(FractNbr)%Init(iInit)%VeloVecIC(3)*Species(FractNbr)%Init(iInit)%BaseVector2IC(3) ) !BV2 component of drift-velocity
    IF ( .NOT.AlmostEqual(A_ins,BV_lengths(1)*BV_lengths(2)) ) THEN
      SWRITE(*,'(A72,2(x,I0),A1)') 'cross product and product of theirs lenghts for BaseVectors of Spec/Init',&
        FractNbr, iInit, ':'
      SWRITE(*,*) A_ins, BV_lengths(1)*BV_lengths(2)
      CALL abort(__STAMP__&
        ,' BaseVectors of the current SpaceIC are not parallel?')
    END IF
    IF ( .NOT.AlmostEqual(SQRT(v_drift_BV(1)**2+v_drift_BV(2)**2+v_drift_line**2),ABS(Species(FractNbr)%Init(iInit)%VeloIC)) ) THEN
      SWRITE(*,'(A60,2(x,I0),A1)') 'v_drift_BV1, v_drift_BV2, v_drift_line, VeloIC for Spec/Init',&
        FractNbr, iInit, ':'
      SWRITE(*,*) v_drift_BV(1),v_drift_BV(2),v_drift_line,Species(FractNbr)%Init(iInit)%VeloIC
      CALL abort(__STAMP__&
        ,' Something is wrong with the Basis of the current SpaceIC!')
    END IF
    IF ( TRIM(Species(FractNbr)%Init(iInit)%SpaceIC) .EQ. 'cuboid_vpi' ) THEN
      lrel_ins_BV=dt * ( (/v_drift_BV(1),-v_drift_BV(1),v_drift_BV(2),-v_drift_BV(2)/) + Species(FractNbr)%Init(iInit)%NSigma &
        * SQRT(BoltzmannConst*Species(FractNbr)%Init(iInit)%MWTemperatureIC/Species(FractNbr)%MassIC) )!rel. virt. insertion height:
      lrel_ins_BV(1:2)=lrel_ins_BV(1:2)/BV_lengths(1)                                                        !... in -BV1/+BV1 dir.
      lrel_ins_BV(3:4)=lrel_ins_BV(3:4)/BV_lengths(2)                                                        !... in -BV2/+BV2 dir.
      DO i=1,4
        IF (.NOT.Species(FractNbr)%Init(iInit)%vpiBVBuffer(i)) lrel_ins_BV(i)=0.
      END DO
    ELSE IF ( TRIM(Species(FractNbr)%Init(iInit)%SpaceIC) .EQ. 'orifice' ) THEN !cylinder-orifice
      lrel_ins_BV(1:4)=dt * ( SQRT(v_drift_BV(1)**2+v_drift_BV(2)**2) + Species(FractNbr)%Init(iInit)%NSigma &
        * SQRT(BoltzmannConst*Species(FractNbr)%Init(iInit)%MWTemperatureIC/Species(FractNbr)%MassIC) ) &
        / Species(FractNbr)%Init(iInit)%RadiusIC !rel. virt. insertion height is a single, maximum value for whole circumference
    END IF
    IF( (TRIM(Species(FractNbr)%Init(iInit)%vpiDomainType).EQ.'orifice') .AND. &
        (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC) .EQ. 'cuboid_vpi') ) THEN !needs further devel.!
      SELECT CASE (GEO%nPeriodicVectors)
      CASE (0)
      CASE (1)
        ParaCheck(1) = Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * GEO%PeriodicVectors(3,1) - &
          Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * GEO%PeriodicVectors(2,1)
        ParaCheck(2) = Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * GEO%PeriodicVectors(1,1) - &
          Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * GEO%PeriodicVectors(3,1)
        ParaCheck(3) = Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * GEO%PeriodicVectors(2,1) - &
          Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * GEO%PeriodicVectors(1,1)
        IF ( .NOT.(SQRT(ParaCheck(1)**2+ParaCheck(2)**2+ParaCheck(3)**2).GT.TwoepsMach) ) orificePeriodic=1 !parallel with BV1
        IF (orificePeriodic .EQ. 0) THEN
          ParaCheck(1) = Species(FractNbr)%Init(iInit)%BaseVector2IC(2) * GEO%PeriodicVectors(3,1) - &
            Species(FractNbr)%Init(iInit)%BaseVector2IC(3) * GEO%PeriodicVectors(2,1)
          ParaCheck(2) = Species(FractNbr)%Init(iInit)%BaseVector2IC(3) * GEO%PeriodicVectors(1,1) - &
            Species(FractNbr)%Init(iInit)%BaseVector2IC(1) * GEO%PeriodicVectors(3,1)
          ParaCheck(3) = Species(FractNbr)%Init(iInit)%BaseVector2IC(1) * GEO%PeriodicVectors(2,1) - &
            Species(FractNbr)%Init(iInit)%BaseVector2IC(2) * GEO%PeriodicVectors(1,1)
          IF ( .NOT.(SQRT(ParaCheck(1)**2+ParaCheck(2)**2+ParaCheck(3)**2).GT.TwoepsMach) ) THEN
            orificePeriodic=2 !parallel with BV2
          ELSE
            CALL abort(__STAMP__&
               ,' PeriodicVector is not parallel to any orifice BV -> not implemented yet!')
          END IF
        ELSE IF (orificePeriodic .EQ. 1) THEN
          !PerVec cannot be parallel with BV2, as BV1 already is
        ELSE
          CALL abort(__STAMP__&
             ,' Something is wrong with the PeriodicVector and the orifice BV!')
        END IF
        lPeri=SQRT(GEO%PeriodicVectors(1,1)**2+GEO%PeriodicVectors(2,1)**2+GEO%PeriodicVectors(3,1)**2)
        IF ( .NOT.AlmostEqual(lPeri,BV_lengths(orificePeriodic)) ) THEN
          SWRITE(*,'(A22,I1,x,A1)') 'lPeri and length of BV',orificePeriodic,':'
          SWRITE(*,'(G0,x,G0)') lPeri,BV_lengths(orificePeriodic)
          CALL abort(__STAMP__&
             ,' PeriodicVector and its parallel BV ar not of same length! ')
        END IF
        orificePeriodicLog(1)=(orificePeriodic.EQ.1)
        orificePeriodicLog(2)=(orificePeriodic.EQ.2)
      CASE DEFAULT
        CALL abort(__STAMP__&
           ,' orifice region only implemented for 0 or 1 PeriodicVector!')
      END SELECT
    END IF !cuboid-orifice
  END IF !freestream or orifice
  !--calculation of (virtual) NbrOfParticle from virt. insertion height
  IF(Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) THEN
    SELECT CASE(TRIM(Species(FractNbr)%Init(iInit)%SpaceIC))
    CASE ('cylinder_vpi')
      IF( TRIM(Species(FractNbr)%Init(iInit)%vpiDomainType) .EQ. 'orifice' ) THEN
        A_ins = PI * ( Species(FractNbr)%Init(iInit)%RadiusIC*(1.+lrel_ins_BV(1)) )**2
      ELSE
        A_ins = PI * ( Species(FractNbr)%Init(iInit)%RadiusIC**2 - Species(FractNbr)%Init(iInit)%Radius2IC**2 )
      END IF
    CASE ('cuboid_vpi')
      IF( (TRIM(Species(FractNbr)%Init(iInit)%vpiDomainType).EQ.'freestream') .OR. &
          (TRIM(Species(FractNbr)%Init(iInit)%vpiDomainType).EQ.'orifice') ) THEN
        A_ins = BV_lengths(1)*(lrel_ins_BV(1)+1.0+lrel_ins_BV(2)) * BV_lengths(2)*(lrel_ins_BV(3)+1.0+lrel_ins_BV(4))
      END IF
    CASE DEFAULT
      CALL abort(__STAMP__,&
                'wrong SpaceIC for virtual Pre-Inserting region!')
    END SELECT
    PartIns = Species(FractNbr)%Init(iInit)%PartDensity * l_ins * A_ins / (Species(FractNbr)%MacroParticleFactor)
    IF(PartIns.GT.0.) THEN
      NbrOfParticle = INT(PartIns)
    ELSE
      NbrOfParticle = 0
    END IF
  END IF
ELSE
  DimSend=3 !save (and send) only positions
END IF

IF ( (NbrOfParticle .LE. 0).AND.(PartIns .LE. 0.) ) RETURN !0<Partins<1: statistical handling of exact REAL-INT-conv. below!

nChunks = 1                   ! Standard: Nicht-MPI
sumOfMatchedParticles = 0
mySumOfMatchedParticles = 0

chunkSize = nbrOfParticle
! process myRank=0 generates the complete list of random positions for all emitted particles
#ifdef MPI
!IF ( (nbrOfParticle.GT.10*PMPIVAR%nProcs                                      ) .AND. &
!IF ( (nbrOfParticle.GT.10*PartMPI%nProcs                                      ) .AND. &
IF ( (nbrOfParticle.LE.PartMPI%InitGroup(InitGroup)%nProcs                             ) .AND. &
     (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).NE.'circle_equidistant'                 ) .AND. &
     (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).NE.'sin_deviation'                      ) .AND. &
     (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).NE.'cuboid_with_equidistant_distribution') .AND. &
     (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).NE.'line_with_equidistant_distribution' )) THEN
   nChunks = PartMPI%InitGroup(InitGroup)%nProcs
ELSE
   nChunks = 1
END IF
IF (mode.EQ.1) THEN
  chunkSize = INT(nbrOfParticle/nChunks)
  IF (PartMPI%InitGroup(InitGroup)%MPIROOT) THEN
    IF( Species(FractNbr)%Init(iInit)%VirtPreInsert .AND. (Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) ) THEN
      ! statistical handling of exact REAL-INT-conversion -> values in send(1)- and receive(2)-mode might differ for VPI+PartDens
      ! (NbrOf Particle can differ from root to other procs and, thus, need to be communicated of calculated later again)
      PartInsRest = PartIns - INT(PartIns)
      CALL RANDOM_NUMBER(RandVal1)
      IF(PartInsRest.GT.RandVal1) PartIns = PartIns + 1.0
      NbrOfParticle = INT(PartIns)
    END IF
    chunkSize = chunkSize + ( nbrOfParticle - (nChunks*chunkSize) )
  END IF
  IF (PartMPI%InitGroup(InitGroup)%MPIROOT .OR. nChunks.GT.1) THEN
#endif
    ALLOCATE( particle_positions(1:chunkSize*DimSend), STAT=allocStat )
    IF (allocStat .NE. 0) THEN
      CALL abort(__STAMP__,&
                'ERROR in SetParticlePosition: cannot allocate particle_positions!')
    END IF

    chunkSize2=chunkSize !will be changed during insertion for:
                         !  1.: vpi with PartDensity (orig. chunksize is for buffer region)
                         !  2.: excludeRegions (orig. chunksize is for SpaceIC without taking excludeRegions into account)
    !------------------SpaceIC-cases: start-----------------------------------------------------------!
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
            CALL abort(__STAMP__,&
              'Error in SetParticlePosition, NormalIC Vektor darf nicht Nullvektor sein')
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
              CALL abort(__STAMP__,&
                'Error in SetParticlePosition, NormalIC should not be zero')
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
              CALL abort(__STAMP__,&
                'Error in SetParticlePosition, NormalIC should not be zero')
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

!        1. determine the z-position in order to get the interpolated curved B-field
         CALL RANDOM_NUMBER(RandVal1)
         IF (NbrOfParticle.EQ.Species(FractNbr)%Init(iInit)%initialParticleNumber) THEN
           particle_positions(i*3  ) = Species(FractNbr)%Init(iInit)%BasePointIC(3) &
                                           + RandVal1 * Species(FractNbr)%Init(iInit)%CuboidHeightIC
         ELSE
           particle_positions(i*3  ) = Species(FractNbr)%Init(iInit)%BasePointIC(3) &
                                           + RandVal1 * dt & 
                                           * Species(FractNbr)%Init(iInit)%VeloIC/Species(FractNbr)%Init(iInit)%alpha 
         END IF

!        2. calculate curved B-field at z-position in order to determin size of gyroradius
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
              CALL abort(__STAMP__,&
                'Error in SetParticlePosition, NormalIC should not be zero')
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
      lineVector(1) = Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * Species(FractNbr)%Init(iInit)%BaseVector2IC(3) - &
        Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * Species(FractNbr)%Init(iInit)%BaseVector2IC(2)
      lineVector(2) = Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * Species(FractNbr)%Init(iInit)%BaseVector2IC(1) - &
        Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * Species(FractNbr)%Init(iInit)%BaseVector2IC(3)
      lineVector(3) = Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * Species(FractNbr)%Init(iInit)%BaseVector2IC(2) - &
        Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * Species(FractNbr)%Init(iInit)%BaseVector2IC(1)
      IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
        CALL abort(__STAMP__,&
          'BaseVectors are parallel!')
      ELSE
        lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + &
          lineVector(3) * lineVector(3))
      END IF
      i=1
      chunkSize2=0
      DO WHILE (i .LE. chunkSize)
         CALL RANDOM_NUMBER(RandVal)
         Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + Species(FractNbr)%Init(iInit)%BaseVector1IC * RandVal(1)
         Particle_pos = Particle_pos + Species(FractNbr)%Init(iInit)%BaseVector2IC * RandVal(2)
         IF (Species(FractNbr)%Init(iInit)%CalcHeightFromDt) THEN !directly calculated by timestep
           Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%VeloIC * dt * RandVal(3)
         ELSE
#if (PP_TimeDiscMethod==201)
!           !scaling due to variable time step (for inlet-condition, but already fixed when %CalcHeightFromDt is used!!!)
!           IF (iter.GT.0) THEN
!             Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%CuboidHeightIC * dt / dt_maxwell * RandVal(3)
!           ELSE
             Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%CuboidHeightIC * RandVal(3) 
!           END IF
#else
           Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%CuboidHeightIC * RandVal(3) 
#endif
         END IF
         IF (Species(FractNbr)%Init(iInit)%NumberOfExcludeRegions.GT.0) THEN
           CALL InsideExcludeRegionCheck(FractNbr, iInit, Particle_pos, insideExcludeRegion)
           IF (insideExcludeRegion) THEN
             i=i+1
             CYCLE !particle is in excluded region
           END IF
         END IF
         particle_positions((chunkSize2+1)*3-2) = Particle_pos(1)
         particle_positions((chunkSize2+1)*3-1) = Particle_pos(2)
         particle_positions((chunkSize2+1)*3  ) = Particle_pos(3)
         i=i+1
         chunkSize2=chunkSize2+1
      END DO
    CASE('cylinder')
      lineVector(1) = Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * Species(FractNbr)%Init(iInit)%BaseVector2IC(3) - &
        Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * Species(FractNbr)%Init(iInit)%BaseVector2IC(2)
      lineVector(2) = Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * Species(FractNbr)%Init(iInit)%BaseVector2IC(1) - &
        Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * Species(FractNbr)%Init(iInit)%BaseVector2IC(3)
      lineVector(3) = Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * Species(FractNbr)%Init(iInit)%BaseVector2IC(2) - &
        Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * Species(FractNbr)%Init(iInit)%BaseVector2IC(1)
      IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
        CALL abort(__STAMP__,&
          'BaseVectors are parallel!')
      ELSE
        lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + &
          lineVector(3) * lineVector(3))
      END IF
      i=1
      chunkSize2=0
      DO WHILE (i .LE. chunkSize)
         radius = Species(FractNbr)%Init(iInit)%RadiusIC + 1.
         DO WHILE((radius.GT.Species(FractNbr)%Init(iInit)%RadiusIC) .OR.(radius.LT.Species(FractNbr)%Init(iInit)%Radius2IC))
            CALL RANDOM_NUMBER(RandVal)
            Particle_pos = Species(FractNbr)%Init(iInit)%BaseVector1IC * (RandVal(1)*2-1) &
                         + Species(FractNbr)%Init(iInit)%BaseVector2IC * (RandVal(2)*2-1)
            radius = SQRT( Particle_pos(1) * Particle_pos(1) + &
                           Particle_pos(2) * Particle_pos(2) + &
                           Particle_pos(3) * Particle_pos(3) )
         END DO
         Particle_pos = Particle_pos + Species(FractNbr)%Init(iInit)%BasePointIC
         IF (Species(FractNbr)%Init(iInit)%CalcHeightFromDt) THEN !directly calculated by timestep
           Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%VeloIC * dt * RandVal(3)
         ELSE
#if (PP_TimeDiscMethod==201)
!           !scaling due to variable time step (for inlet-condition, but already fixed when %CalcHeightFromDt is used!!!)
!           IF (iter.GT.0) THEN
!             Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%CylinderHeightIC * dt/dt_maxwell * RandVal(3)
!           ELSE
             Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%CylinderHeightIC * RandVal(3)
!           END IF          
#else
           Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%CylinderHeightIC * RandVal(3)
#endif
         END IF
         IF (Species(FractNbr)%Init(iInit)%NumberOfExcludeRegions.GT.0) THEN
           CALL InsideExcludeRegionCheck(FractNbr, iInit, Particle_pos, insideExcludeRegion)
           IF (insideExcludeRegion) THEN
             i=i+1
             CYCLE !particle is in excluded region
           END IF
         END IF
         particle_positions((chunkSize2+1)*3-2) = Particle_pos(1)
         particle_positions((chunkSize2+1)*3-1) = Particle_pos(2)
         particle_positions((chunkSize2+1)*3  ) = Particle_pos(3)
         i=i+1
         chunkSize2=chunkSize2+1
      END DO
    CASE('cuboid_vpi')
      i=1
      chunkSize2=0
      DO WHILE (i .LE. chunkSize)          
        ! Check if particle would reach comp. domain in one timestep
        CALL CalcVelocity_maxwell_lpn(FractNbr, iInit, Vec3D)
        CALL RANDOM_NUMBER(RandVal)
        v_line = Vec3D(1)*lineVector(1) + Vec3D(2)*lineVector(2) + Vec3D(3)*lineVector(3) !lineVector component of velocity
        delta_l = dt * v_line - l_ins * RandVal(3)
        IF (delta_l .LT. 0.) THEN
          IF (Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) i=i+1
          CYCLE !particle would not reach comp. domain -> try new velo
        END IF
        SELECT CASE(TRIM(Species(FractNbr)%Init(iInit)%vpiDomainType)) 
        CASE('perpendicular_extrusion')
          ! set particle positions depending on SpaceIC
          Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + Species(FractNbr)%Init(iInit)%BaseVector1IC * RandVal(1) &
                                                                   + Species(FractNbr)%Init(iInit)%BaseVector2IC * RandVal(2) &
                                                                   + lineVector * delta_l
        CASE('freestream')
          v_BV(1) = Vec3D(1)*Species(FractNbr)%Init(iInit)%BaseVector1IC(1) &
                  + Vec3D(2)*Species(FractNbr)%Init(iInit)%BaseVector1IC(2) &
                  + Vec3D(3)*Species(FractNbr)%Init(iInit)%BaseVector1IC(3)
          v_BV(1) = v_BV(1) / BV_lengths(1) !BV1 component of velocity
          v_BV(2) = Vec3D(1)*Species(FractNbr)%Init(iInit)%BaseVector2IC(1) &
                  + Vec3D(2)*Species(FractNbr)%Init(iInit)%BaseVector2IC(2) &
                  + Vec3D(3)*Species(FractNbr)%Init(iInit)%BaseVector2IC(3)
          v_BV(2) = v_BV(2) / BV_lengths(2) !BV2 component of velocity
          delta_lBV = dt * v_BV
          delta_lBV(1) = delta_lBV(1) + BV_lengths(1) * ( RandVal(1)*(lrel_ins_BV(1)+1.0+lrel_ins_BV(2)) - lrel_ins_BV(1) )
          delta_lBV(2) = delta_lBV(2) + BV_lengths(2) * ( RandVal(2)*(lrel_ins_BV(3)+1.0+lrel_ins_BV(4)) - lrel_ins_BV(3) )
          IF ( (delta_lBV(1).LT.0.) .OR. (delta_lBV(1).GT.BV_lengths(1)) ) THEN
            IF (Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) i=i+1
            CYCLE !particle would not reach comp. domain in direction of BaseVector1 -> try new velo
          END IF
          IF ( (delta_lBV(2).LT.0.) .OR. (delta_lBV(2).GT.BV_lengths(2)) ) THEN
            IF (Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) i=i+1
            CYCLE !particle would not reach comp. domain in direction of BaseVector2 -> try new velo
          END IF
          ! set particle positions depending on SpaceIC
          Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC &
                       + Species(FractNbr)%Init(iInit)%BaseVector1IC/BV_lengths(1) * delta_lBV(1) &
                       + Species(FractNbr)%Init(iInit)%BaseVector2IC/BV_lengths(2) * delta_lBV(2) &
                       + lineVector * delta_l 
        CASE('orifice')
          ! set particle position (to be tried) depending on SpaceIC
          v_BV(1) = Vec3D(1)*Species(FractNbr)%Init(iInit)%BaseVector1IC(1) &
                  + Vec3D(2)*Species(FractNbr)%Init(iInit)%BaseVector1IC(2) &
                  + Vec3D(3)*Species(FractNbr)%Init(iInit)%BaseVector1IC(3)
          v_BV(1) = v_BV(1) / BV_lengths(1) !BV1 component of velocity
          v_BV(2) = Vec3D(1)*Species(FractNbr)%Init(iInit)%BaseVector2IC(1) &
                  + Vec3D(2)*Species(FractNbr)%Init(iInit)%BaseVector2IC(2) &
                  + Vec3D(3)*Species(FractNbr)%Init(iInit)%BaseVector2IC(3)
          v_BV(2) = v_BV(2) / BV_lengths(2) !BV2 component of velocity
          delta_lBV = dt * v_BV
          delta_lBV(1) = delta_lBV(1) + BV_lengths(1) * ( RandVal(1)*(lrel_ins_BV(1)+1.0+lrel_ins_BV(2)) - lrel_ins_BV(1) )
          delta_lBV(2) = delta_lBV(2) + BV_lengths(2) * ( RandVal(2)*(lrel_ins_BV(3)+1.0+lrel_ins_BV(4)) - lrel_ins_BV(3) )
          Particle_pos = Species(FractNbr)%Init(iInit)%BaseVector1IC/BV_lengths(1) * delta_lBV(1) &
                       + Species(FractNbr)%Init(iInit)%BaseVector2IC/BV_lengths(2) * delta_lBV(2) &
                       + lineVector * delta_l 
          IntersecPoint = Particle_pos - Vec3D * delta_l/v_line !Vector from BP to Intersec point of virt. path with orifice plane
          orifice_delta = (IntersecPoint(1)*Species(FractNbr)%Init(iInit)%BaseVector1IC(1) &
                         + IntersecPoint(2)*Species(FractNbr)%Init(iInit)%BaseVector1IC(2) &
                         + IntersecPoint(3)*Species(FractNbr)%Init(iInit)%BaseVector1IC(3))/BV_lengths(1)
          IF ( orificePeriodicLog(1) ) THEN
            IF ( (delta_lBV(1).LT.0.) .OR. (delta_lBV(1).GT.BV_lengths(1)) ) THEN
              IF (Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) i=i+1
              CYCLE !particle would not reach comp. domain in direction of BaseVector1 -> try new velo
            END IF
          ELSE
            IF ( (orifice_delta.GT.BV_lengths(1)) .OR. (orifice_delta.LT.0.) ) THEN
              IF (Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) i=i+1
              CYCLE !particle would not reach comp. through orifice -> try new velo
            END IF
          END IF
          orifice_delta = (IntersecPoint(1)*Species(FractNbr)%Init(iInit)%BaseVector2IC(1) &
                         + IntersecPoint(2)*Species(FractNbr)%Init(iInit)%BaseVector2IC(2) &
                         + IntersecPoint(3)*Species(FractNbr)%Init(iInit)%BaseVector2IC(3))/BV_lengths(2)
          IF ( orificePeriodicLog(2) ) THEN
            IF ( (delta_lBV(2).LT.0.) .OR. (delta_lBV(2).GT.BV_lengths(2)) ) THEN
              IF (Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) i=i+1
              CYCLE !particle would not reach comp. domain in direction of BaseVector2 -> try new velo
            END IF
          ELSE
            IF ( (orifice_delta.GT.BV_lengths(2)) .OR. (orifice_delta.LT.0.) ) THEN
              IF (Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) i=i+1
              CYCLE !particle would not reach comp. through orifice -> try new velo
            END IF
          END IF
          Particle_pos = Particle_pos + Species(FractNbr)%Init(iInit)%BasePointIC
        CASE DEFAULT
          CALL abort(__STAMP__,&
                'wrong vpiDomainType for virtual Pre-Inserting region!')
        END SELECT
        IF (Species(FractNbr)%Init(iInit)%NumberOfExcludeRegions.GT.0) THEN
          CALL InsideExcludeRegionCheck(FractNbr, iInit, Particle_pos, insideExcludeRegion)
          IF (insideExcludeRegion) THEN
            i=i+1
            CYCLE !particle is in excluded region
          END IF
        END IF
        !store determined values or go to next particle
        particle_positions((chunkSize2+1)*6-5) = Particle_pos(1)
        particle_positions((chunkSize2+1)*6-4) = Particle_pos(2)
        particle_positions((chunkSize2+1)*6-3) = Particle_pos(3)
        particle_positions((chunkSize2+1)*6-2) = Vec3D(1)
        particle_positions((chunkSize2+1)*6-1) = Vec3D(2)
        particle_positions((chunkSize2+1)*6  ) = Vec3D(3)
        i=i+1
        chunkSize2=chunkSize2+1
      END DO
    CASE('cylinder_vpi')
      i=1
      chunkSize2=0
      DO WHILE (i .LE. chunkSize)        
        ! Check if particle would reach comp. domain in one timestep
        CALL CalcVelocity_maxwell_lpn(FractNbr, iInit, Vec3D)
        CALL RANDOM_NUMBER(RandVal)
        v_line = Vec3D(1)*lineVector(1) + Vec3D(2)*lineVector(2) + Vec3D(3)*lineVector(3) !lineVector component of velocity
        delta_l = dt * v_line - l_ins * RandVal(3)
        IF (delta_l .LT. 0.) THEN
          IF (Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) i=i+1
          CYCLE !particle would not reach comp. domain -> try new velo
        END IF
        SELECT CASE(TRIM(Species(FractNbr)%Init(iInit)%vpiDomainType)) 
        CASE('perpendicular_extrusion')
          ! set particle positions depending on SpaceIC
          radius = Species(FractNbr)%Init(iInit)%RadiusIC + 1.
          DO WHILE((radius.GT.Species(FractNbr)%Init(iInit)%RadiusIC) .OR.(radius.LT.Species(FractNbr)%Init(iInit)%Radius2IC))
            CALL RANDOM_NUMBER(RandVal)
            Particle_pos = Species(FractNbr)%Init(iInit)%BaseVector1IC * (RandVal(1)*2-1) &
                         + Species(FractNbr)%Init(iInit)%BaseVector2IC * (RandVal(2)*2-1)
            radius = SQRT( Particle_pos(1) * Particle_pos(1) + &
                           Particle_pos(2) * Particle_pos(2) + &
                           Particle_pos(3) * Particle_pos(3) )
          END DO
          Particle_pos = Particle_pos + Species(FractNbr)%Init(iInit)%BasePointIC + lineVector * delta_l
        CASE('orifice')
          ! set particle position (to be tried) depending on SpaceIC
          radius = Species(FractNbr)%Init(iInit)%RadiusIC * (1.+lrel_ins_BV(1)) + 1. !lrel_ins_BV(1)=lrel_ins_BV(2)=...
          DO WHILE ( radius .GT. Species(FractNbr)%Init(iInit)%RadiusIC * (1.+lrel_ins_BV(1)) )
            CALL RANDOM_NUMBER(RandVal)
            Particle_pos = Species(FractNbr)%Init(iInit)%BaseVector1IC*(1.+lrel_ins_BV(1)) * (RandVal(1)*2-1) &
                         + Species(FractNbr)%Init(iInit)%BaseVector2IC*(1.+lrel_ins_BV(2)) * (RandVal(2)*2-1) !BV_lengths(:)=R
            radius = SQRT( Particle_pos(1) * Particle_pos(1) + &
                           Particle_pos(2) * Particle_pos(2) + &
                           Particle_pos(3) * Particle_pos(3) )
          END DO
          Particle_pos = Particle_pos + dt * Vec3D - lineVector * (dt * v_line - delta_l) !get old RandVal(3) for l_ins * RandVal(3)
          
          IntersecPoint = Particle_pos - Vec3D * delta_l/v_line !Vector from BP to Intersec point of virt. path with orifice plane
          orifice_delta = (IntersecPoint(1)*Species(FractNbr)%Init(iInit)%BaseVector1IC(1) &
                         + IntersecPoint(2)*Species(FractNbr)%Init(iInit)%BaseVector1IC(2) &
                         + IntersecPoint(3)*Species(FractNbr)%Init(iInit)%BaseVector1IC(3))/BV_lengths(1)
          radius = (IntersecPoint(1)*Species(FractNbr)%Init(iInit)%BaseVector2IC(1) &
                  + IntersecPoint(2)*Species(FractNbr)%Init(iInit)%BaseVector2IC(2) &
                  + IntersecPoint(3)*Species(FractNbr)%Init(iInit)%BaseVector2IC(3))/BV_lengths(2)
          radius = SQRT( orifice_delta*orifice_delta + radius*radius )
            IF ( radius.GT.BV_lengths(1) ) THEN
              IF (Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) i=i+1
              CYCLE !particle would not reach comp. through orifice -> try new velo
            END IF
          Particle_pos = Particle_pos + Species(FractNbr)%Init(iInit)%BasePointIC
        CASE DEFAULT
          CALL abort(__STAMP__,&
                'wrong vpiDomainType for virtual Pre-Inserting region!')
        END SELECT
        IF (Species(FractNbr)%Init(iInit)%NumberOfExcludeRegions.GT.0) THEN
          CALL InsideExcludeRegionCheck(FractNbr, iInit, Particle_pos, insideExcludeRegion)
          IF (insideExcludeRegion) THEN
            i=i+1
            CYCLE !particle is in excluded region
          END IF
        END IF
        !store determined values or go to next particle
        particle_positions((chunkSize2+1)*6-5) = Particle_pos(1)
        particle_positions((chunkSize2+1)*6-4) = Particle_pos(2)
        particle_positions((chunkSize2+1)*6-3) = Particle_pos(3)
        particle_positions((chunkSize2+1)*6-2) = Vec3D(1)
        particle_positions((chunkSize2+1)*6-1) = Vec3D(2)
        particle_positions((chunkSize2+1)*6  ) = Vec3D(3)
        i=i+1
        chunkSize2=chunkSize2+1
      END DO
!    CASE('LD_insert')
!      CALL LD_SetParticlePosition(chunkSize,particle_positions_Temp,FractNbr,iInit)
!      DEALLOCATE( particle_positions, STAT=allocStat )
!      IF (allocStat .NE. 0) THEN
!        CALL abort(__STAMP__,&
!          'ERROR in ParticleEmission_parallel: cannot deallocate particle_positions!')
!      END IF
!      ALLOCATE(particle_positions(3*chunkSize))
!      particle_positions(1:3*chunkSize) = particle_positions_Temp(1:3*chunkSize)
!      DEALLOCATE( particle_positions_Temp, STAT=allocStat )
!      IF (allocStat .NE. 0) THEN
!        CALL abort(__STAMP__,&
!          'ERROR in ParticleEmission_parallel: cannot deallocate particle_positions!')
!      END IF
    CASE('cuboid_equal')
#ifdef MPI
      IF (PartMPI%InitGroup(InitGroup)%nProcs.GT. 1) THEN
        SWRITE(UNIT_stdOut,*)'WARNING in SetParticlePosition:'
        SWRITE(UNIT_stdOut,*)'cannot fully handle Particle Initial Condition \"cuboid equal\"'
        SWRITE(UNIT_stdOut,*)'in parallel mode (with more than one CPU)!'
        SWRITE(UNIT_stdOut,*)'USE WITH CARE!!!'
      END IF
      j=0
      mySumOfMatchedParticles = 0
      DO i=1,PDM%ParticleVecLength
         j=j+1
         ParticleIndexNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
         IF (ParticleIndexNbr .ne. 0) THEN
            PartState(ParticleIndexNbr,1:3) = PartState(j,1:3)
            PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
            IF(DoRefMapping)THEN
              CALL SingleParticleToExactElement(ParticleIndexNbr,doHALO=.FALSE.)
            ELSE
              CALL SingleParticleToExactElementNoMap(ParticleIndexNbr,doHALO=.FALSE.)
            END IF
            IF (PDM%ParticleInside(ParticleIndexNbr)) THEN
               mySumOfMatchedParticles = mySumOfMatchedParticles + 1
            ELSE
               PDM%ParticleInside(ParticleIndexNbr) = .FALSE.
            END IF
         ELSE
            CALL abort(__STAMP__,&
                'ERROR in SetParticlePosition: ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?')
         END IF
      END DO
      CALL MPI_ALLREDUCE(mySumOfMatchedParticles, sumOfMatchedParticles, 1, MPI_INTEGER &
          , MPI_SUM, PartMPI%InitGroup(InitGroup)%COMM, IERROR)
      nbrOfParticle = NbrOfParticle - sumOfMatchedParticles
      IF (nbrOfParticle .NE. 0) THEN
        IPWRITE(UNIT_stdOut,*)'ERROR in ParticleEmission_parallel:'
        IPWRITE(UNIT_stdOut,'(A,I8,A)')'matched ', sumOfMatchedParticles, ' particles'
        IPWRITE(UNIT_stdOut,'(A,I8,A)')'when ', NbrOfParticle+sumOfMatchedParticles, ' particles were required!'
        CALL abort(__STAMP__,&
            'ERROR in ParticleEmission_parallel')
      END IF
      NbrOfParticle = mySumOfMatchedParticles
      DEALLOCATE( particle_positions, STAT=allocStat )
      IF (allocStat .NE. 0) THEN
        CALL abort(__STAMP__,&
            'ERROR in ParticleEmission_parallel: cannot deallocate particle_positions!')
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
         SWRITE(*,*) 'for species ',FractNbr,' does not match number of particles in each direction!'
         CALL abort(__STAMP__,&
            'ERROR: Number of particles in init / emission region',iInit)
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
          CALL abort(__STAMP__,&
            'Basevectors1IC,-2IC and CuboidHeightIC have to be in x,y,z-direction, respectively for emission condition')
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
         SWRITE(*,*) 'for species ',FractNbr,' does not match number of particles in each direction!'
         CALL abort(__STAMP__,&
            'ERROR: Number of particles in init / emission region',iInit)
       END IF
       xlen = abs(GEO%xmaxglob  - GEO%xminglob)  
       ylen = abs(GEO%ymaxglob  - GEO%yminglob)
       zlen = abs(GEO%zmaxglob  - GEO%zminglob)
       pilen=2.0*PI/xlen
       x_step = xlen/Species(FractNbr)%Init(iInit)%maxParticleNumberX
       y_step = ylen/Species(FractNbr)%Init(iInit)%maxParticleNumberY
       z_step = zlen/Species(FractNbr)%Init(iInit)%maxParticleNumberZ
       iPart = 1
       DO i=1,Species(FractNbr)%Init(iInit)%maxParticleNumberX
          x_pos = (i * x_step - x_step*0.5)
          x_pos = GEO%xminglob + x_pos + Species(FractNbr)%Init(iInit)%Amplitude &
                  * sin(Species(FractNbr)%Init(iInit)%WaveNumber * pilen * x_pos)
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
    !------------------SpaceIC-cases: end-----------------------------------------------------------!
    chunkSize=chunkSize2

#ifdef MPI
 ELSE !no mpi root, nchunks=1
   chunkSize=0
 END IF

 IF(nChunks.GT.1) THEN
   ALLOCATE( PartMPIInsert%nPartsSend  (0:PartMPI%InitGroup(InitGroup)%nProcs-1), STAT=allocStat )
   ALLOCATE( PartMPIInsert%nPartsRecv  (0:PartMPI%InitGroup(InitGroup)%nProcs-1), STAT=allocStat )
   ALLOCATE( PartMPIInsert%SendRequest (0:PartMPI%InitGroup(InitGroup)%nProcs-1,1:2), STAT=allocStat )
   ALLOCATE( PartMPIInsert%RecvRequest (0:PartMPI%InitGroup(InitGroup)%nProcs-1,1:2), STAT=allocStat )
   ALLOCATE( PartMPIInsert%send_message(0:PartMPI%InitGroup(InitGroup)%nProcs-1), STAT=allocStat )
   PartMPIInsert%nPartsSend(:)=0
   DO i=1,chunkSize
     CellX = INT((particle_positions(DimSend*(i-1)+1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
     CellY = INT((particle_positions(DimSend*(i-1)+2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
     CellZ = INT((particle_positions(DimSend*(i-1)+3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
     InsideMyBGM=.TRUE.
     IF ((CellX.GT.GEO%FIBGMimax).OR.(CellX.LT.GEO%FIBGMimin) .OR. &
         (CellY.GT.GEO%FIBGMjmax).OR.(CellY.LT.GEO%FIBGMjmin) .OR. &
         (CellZ.GT.GEO%FIBGMkmax).OR.(CellZ.LT.GEO%FIBGMkmin)) THEN
       InsideMyBGM=.FALSE.
     END If
     IF (InsideMyBGM) THEN
       IF (.NOT.ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) InsideMyBGM=.FALSE.
     END IF
     IF (InsideMyBGM) THEN
       DO j=2,GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(1)+1
         iProc=GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(j)
         tProc=PartMPI%InitGroup(InitGroup)%CommToGroup(iProc)
         IF(tProc.EQ.-1)CYCLE
         !IF(PartMPI%InitGroup(InitGroup)%COMM.EQ.MPI_COMM_NULL) THEN
         PartMPIInsert%nPartsSend(tProc)=PartMPIInsert%nPartsSend(tProc)+1
       END DO
       PartMPIInsert%nPartsSend(PartMPI%InitGroup(InitGroup)%MyRank)=&
              PartMPIInsert%nPartsSend(PartMPI%InitGroup(InitGroup)%MyRank)+1
     ELSE
       DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
!         IF (iProc.EQ.PMPIVAR%iProc) CYCLE
         PartMPIInsert%nPartsSend(iProc)=PartMPIInsert%nPartsSend(iProc)+1
       END DO
     END IF
   END DO   
 ELSE
    IF(PartMPI%InitGroup(InitGroup)%MPIRoot) THEN
      ALLOCATE( PartMPIInsert%send_message(0:0), STAT=allocStat )
      MessageSize=DimSend*chunkSize
      ALLOCATE( PartMPIInsert%send_message(0)%content(1:MessageSize), STAT=allocStat )
      PartMPIInsert%send_message(0)%content(:)=particle_positions(1:DimSend*chunkSize)
      DEALLOCATE(particle_positions, STAT=allocStat)
    END IF
 END IF
 IF (nChunks.GT.1) THEN
    DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
      !--- MPI_ISEND lengths of lists of particles leaving local mesh
      CALL MPI_ISEND(PartMPIInsert%nPartsSend(iProc), 1, MPI_INTEGER, iProc, 1011+FractNbr, PartMPI%InitGroup(InitGroup)%COMM, &
                     PartMPIInsert%SendRequest(iProc,1), IERROR)
      !--- MPI_IRECV lengths of lists of particles entering local mesh
      CALL MPI_IRECV(PartMPIInsert%nPartsRecv(iProc), 1, MPI_INTEGER, iProc, 1011+FractNbr, PartMPI%InitGroup(InitGroup)%COMM, &
                     PartMPIInsert%RecvRequest(iProc,1), IERROR)
      IF (PartMPIInsert%nPartsSend(iProc).GT.0) THEN
        ALLOCATE( PartMPIInsert%send_message(iProc)%content(1:DimSend*PartMPIInsert%nPartsSend(iProc)), STAT=allocStat )
      END IF
    END DO
    PartMPIInsert%nPartsSend(:)=0
    DO i=1,chunkSize
      CellX = INT((particle_positions(DimSend*(i-1)+1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
      CellY = INT((particle_positions(DimSend*(i-1)+2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
      CellZ = INT((particle_positions(DimSend*(i-1)+3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
      InsideMyBGM=.TRUE.
      IF ((CellX.GT.GEO%FIBGMimax).OR.(CellX.LT.GEO%FIBGMimin) .OR. &
          (CellY.GT.GEO%FIBGMjmax).OR.(CellY.LT.GEO%FIBGMjmin) .OR. &
          (CellZ.GT.GEO%FIBGMkmax).OR.(CellZ.LT.GEO%FIBGMkmin)) THEN
        InsideMyBGM=.FALSE.
      END If
      IF (InsideMyBGM) THEN
        IF (.NOT.ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) InsideMyBGM=.FALSE.
      END IF
      IF (InsideMyBGM) THEN
        DO j=2,GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(1)+1
          iProc=GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(j)
          tProc=PartMPI%InitGroup(InitGroup)%CommToGroup(iProc)
          IF(tProc.EQ.-1)CYCLE
          PartMPIInsert%nPartsSend(tProc)=PartMPIInsert%nPartsSend(tProc)+1
          k=PartMPIInsert%nPartsSend(tProc)
          PartMPIInsert%send_message(tProc)%content(DimSend*(k-1)+1:DimSend*k)=particle_positions(DimSend*(i-1)+1:DimSend*i)
        END DO
        PartMPIInsert%nPartsSend(PartMPI%InitGroup(InitGroup)%MyRank)= &
            PartMPIInsert%nPartsSend(PartMPI%InitGroup(InitGroup)%MyRank)+1
        k=PartMPIInsert%nPartsSend(PartMPI%InitGroup(InitGroup)%MyRank)
        PartMPIInsert%send_message(PartMPI%InitGroup(InitGroup)%MyRank)%content(DimSend*(k-1)+1:DimSend*k)=&
                                                                          particle_positions(DimSend*(i-1)+1:DimSend*i)
      ELSE
        DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
 !         IF (iProc.EQ.PMPIVAR%iProc) CYCLE
          PartMPIInsert%nPartsSend(iProc)=PartMPIInsert%nPartsSend(iProc)+1
          k=PartMPIInsert%nPartsSend(iProc)
          PartMPIInsert%send_message(iProc)%content(DimSend*(k-1)+1:DimSend*k)=particle_positions(DimSend*(i-1)+1:DimSend*i)
        END DO
      END IF
    END DO
    DEALLOCATE(particle_positions, STAT=allocStat)
    DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
      !--- (non-blocking:) send messages to all procs receiving particles from myself
      IF (PartMPIInsert%nPartsSend(iProc).GT.0) THEN
        CALL MPI_ISEND(PartMPIInsert%send_message(iProc)%content, DimSend*PartMPIInsert%nPartsSend(iProc),& 
         MPI_DOUBLE_PRECISION, iProc, 1022+FractNbr, PartMPI%InitGroup(InitGroup)%COMM, PartMPIInsert%SendRequest(iProc,2), IERROR)
      END IF
    END DO
  END IF
ELSE ! mode.NE.1:
!--- RECEIVE:
  IF(nChunks.EQ.1) THEN
    IF(PartMPI%InitGroup(InitGroup)%MPIRoot) THEN !chunkSize can be 1 higher than NbrOfParticle for VPI+PartDens
       chunkSize=INT( REAL(SIZE(PartMPIInsert%send_message(0)%content)) / REAL(DimSend) )
       ALLOCATE(particle_positions(1:chunkSize*DimSend), STAT=allocStat)
       particle_positions(:)=PartMPIInsert%send_message(0)%content(:)
       DEALLOCATE( PartMPIInsert%send_message(0)%content )      
       DEALLOCATE( PartMPIInsert%send_message )
    END IF
    IF( Species(FractNbr)%Init(iInit)%VirtPreInsert .AND. (Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) ) THEN
      CALL MPI_BCAST(chunkSize, 1, MPI_INTEGER,0,PartMPI%InitGroup(InitGroup)%COMM,IERROR)
    ELSE
      chunkSize=NbrOfParticle
    END IF
    IF(.NOT.PartMPI%InitGroup(InitGroup)%MPIROOT) THEN
      ALLOCATE(particle_positions(1:chunkSize*DimSend), STAT=allocStat)
    END IF
    CALL MPI_BCAST(particle_positions, chunkSize*DimSend, MPI_DOUBLE_PRECISION,0,PartMPI%InitGroup(InitGroup)%COMM,IERROR)
  ELSE   
    DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
      CALL MPI_WAIT(PartMPIInsert%RecvRequest(iProc,1),msg_status(:),IERROR)
    END DO
    k=SUM(PartMPIInsert%nPartsRecv)
    ALLOCATE(particle_positions(1:k*DimSend), STAT=allocStat)
    k=0
    DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
      IF (PartMPIInsert%nPartsRecv(iProc).GT.0) THEN
      !--- MPI_IRECV lengths of lists of particles entering local mesh
        CALL MPI_IRECV(particle_positions(k*DimSend+1), DimSend*PartMPIInsert%nPartsRecv(iProc),&
                                                  MPI_DOUBLE_PRECISION, iProc, 1022+FractNbr,   &
                                                  PartMPI%InitGroup(InitGroup)%COMM, PartMPIInsert%RecvRequest(iProc,2), IERROR)
        CALL MPI_WAIT(PartMPIInsert%RecvRequest(iProc,2),msg_status(:),IERROR)
        k=k+PartMPIInsert%nPartsRecv(iProc)
      END IF
    END DO
    DEALLOCATE( PartMPIInsert%nPartsRecv )
    DEALLOCATE( PartMPIInsert%RecvRequest )
    DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
      CALL MPI_WAIT(PartMPIInsert%SendRequest(iProc,1),msg_status(:),IERROR)
      IF (PartMPIInsert%nPartsSend(iProc).GT.0) THEN
        CALL MPI_WAIT(PartMPIInsert%SendRequest(iProc,2),msg_status(:),IERROR)
        DEALLOCATE( PartMPIInsert%send_message(iProc)%content )
      END IF
    END DO
    DEALLOCATE( PartMPIInsert%nPartsSend )
    DEALLOCATE( PartMPIInsert%send_message )
    DEALLOCATE( PartMPIInsert%SendRequest )
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
       PartState(ParticleIndexNbr,1:DimSend) = particle_positions(DimSend*(i-1)+1:DimSend*(i-1)+DimSend)
       PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
       IF(DoRefMapping)THEN
         CALL SingleParticleToExactElement(ParticleIndexNbr,doHALO=.FALSE.)
       ELSE
         CALL SingleParticleToExactElementNoMap(ParticleIndexNbr,doHALO=.FALSE.)
       END IF
       !CALL SingleParticleToExactElement(ParticleIndexNbr)
       IF (PDM%ParticleInside(ParticleIndexNbr)) THEN
          !IPWRITE(UNIT_stdOut,*) 'partpos',PartState(ParticleIndexNbr,1:3)
          mySumOfMatchedParticles = mySumOfMatchedParticles + 1
       ELSE
          PDM%ParticleInside(ParticleIndexNbr) = .FALSE.
       END IF
    ELSE
       CALL abort(__STAMP__,&
            'ERROR in SetParticlePosition:ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?')
    END IF
  END DO
 
  IF(DoDisplayIter)THEN
    IF(MOD(iter,IterDisplayStep).EQ.0) THEN
#ifdef MPI
  !   WRITE(*,*)'mySumOfMatchedParticles=',mySumOfMatchedParticles
     ! check the sum of the matched particles: did each particle find its "home"-CPU?
      CALL MPI_ALLREDUCE(mySumOfMatchedParticles, sumOfMatchedParticles, 1, MPI_INTEGER, MPI_SUM &
                                               , PartMPI%InitGroup(InitGroup)%COMM, IERROR)
#else
      ! im seriellen Fall kommen alle Partikel auf einen CPU,
      ! daher ist PIC%maxParticleNumber die harte Grenze
      sumOfMatchedParticles = mySumOfMatchedParticles
#endif

#ifdef MPI
      IF(PartMPI%InitGroup(InitGroup)%MPIRoot) THEN
#endif
        IF( Species(FractNbr)%Init(iInit)%VirtPreInsert .AND. (Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) ) THEN
          IF ((nbrOfParticle .NE. sumOfMatchedParticles).AND.OutputVpiWarnings) THEN
            SWRITE(*,'(A)')'WARNING in ParticleEmission_parallel:'
            SWRITE(*,'(A,I0)')'Fraction Nbr: ', FractNbr
            SWRITE(*,'(I8,A)') sumOfMatchedParticles, ' particles reached the domain when'
            SWRITE(*,'(I8,A)') NbrOfParticle, '(+1) velocities were calculated with vpi+PartDens'
          END IF
        ELSE
          IF (nbrOfParticle .GT. sumOfMatchedParticles) THEN
            SWRITE(*,*)'WARNING in ParticleEmission_parallel:'
            SWRITE(*,'(A,I0)')'Fraction Nbr: ', FractNbr
            SWRITE(*,'(A,I8,A)')'matched only ', sumOfMatchedParticles, ' particles'
            SWRITE(*,'(A,I8,A)')'when ', NbrOfParticle, ' particles were required!'
          ELSE IF (nbrOfParticle .LT. sumOfMatchedParticles) THEN
            SWRITE(*,*)'ERROR in ParticleEmission_parallel:'
            SWRITE(*,'(A,I0)')'Fraction Nbr: ', FractNbr
            SWRITE(*,'(A,I8,A)')'matched ', sumOfMatchedParticles, ' particles'
            SWRITE(*,'(A,I8,A)')'when ', NbrOfParticle, ' particles were required!'
#if (PP_TimeDiscMethod!=1000) && (PP_TimeDiscMethod!=1001)
            CALL abort(__STAMP__,&
              'selected timedisk does not allow num of inserted part .gt. required')
#endif
          ELSE IF (nbrOfParticle .EQ. sumOfMatchedParticles) THEN
!           SWRITE(UNIT_stdOut,'(A,I0)')'Fraction Nbr: ', FractNbr
!           SWRITE(UNIT_stdOut,'(A,I0,A)')'ParticleEmission_parallel: matched all (',NbrOfParticle,') particles!'
          END IF
        END IF
#ifdef MPI
      END IF ! PMPIVAR%iProc.EQ.0
#endif
    END IF ! IterDisplayStep
  END IF
   ! Return the *local* NbrOfParticle so that the following Routines only fill in
   ! the values for the local particles
  NbrOfParticle = mySumOfMatchedParticles

  DEALLOCATE( particle_positions, STAT=allocStat )
  IF (allocStat .NE. 0) THEN
    CALL abort(__STAMP__,&
              'ERROR in ParticleEmission_parallel: cannot deallocate particle_positions!')
  END IF
#ifdef MPI
END IF ! mode 1/2
#endif

END SUBROUTINE SetParticlePosition

SUBROUTINE SetParticleVelocity(FractNbr,iInit,NbrOfParticle,Is_BGGas_opt)
!===================================================================================================================================
! Determine the particle velocity of each inserted particle
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Timedisc_Vars, ONLY : dt
USE MOD_Equation_Vars, ONLY : c
USE MOD_Particle_Vars
USE MOD_PIC_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: FractNbr,iInit                                                   
LOGICAL, OPTIONAL                :: Is_BGGas_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)            :: NbrOfParticle            
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: i,j,PositionNbr        
REAL                             :: Radius(3), n_vec(3), tan_vec(3), Velo1, Angle, Velo2, f
REAL                             :: Vec3D(3), RandVal(3)
REAL                             :: II(3,3),JJ(3,3),NN(3,3)
INTEGER                          :: distnum
REAL                             :: r1,r2,x_1,x_2,y_1,y_2,a,b,e,g,x_01,x_02,y_01,y_02, RandVal1
REAL                             :: Velosq, v_sum(3), v2_sum, maxwellfac
LOGICAL                          :: Is_BGGas
REAL                             :: sigma(3), ftl, PartVelo 
!===================================================================================================================================

IF (PRESENT(Is_BGGas_opt)) THEN
  Is_BGGas=Is_BGGas_opt
ELSE
  Is_BGGas=.FALSE.
END IF

IF(NbrOfParticle.lt.1) RETURN
IF(Species(FractNbr)%Init(iInit)%VirtPreInsert) RETURN !velocities already set in SetParticlePosition!
 
!   IF(NbrOfParticle.gt.PIC%maxParticleNumber)THEN
!      NbrOfParticle = PIC%maxParticleNumber
!   END IF

SELECT CASE(TRIM(Species(FractNbr)%Init(iInit)%velocityDistribution))
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
!  CASE('EOC_Test')
!     ! for leapfrog EOC test velo has to be set half dt before ICPos.
!     Radius = Species(1)%BasePointIC 
!     IF (Species(1)%BasePointIC(1) > 0. ) THEN
!       n_vec = (/0.,0.,-1./) 
!       Angle = 0.
!       Angle = Angle - PIC%GyrationFrequency * dt * 0.5
!       Radius(1) = cos(Angle)
!       Radius(2) = sin(Angle)   
!     ELSEIF (Species(1)%BasePointIC(3) > 0. ) THEN  
!       n_vec = (/0.,1.,0./) 
!       Angle = PI
!       Angle = Angle - PIC%GyrationFrequency * dt * 0.5
!       Radius(1) = sin(Angle)
!       Radius(3) = cos(Angle)   
!     END IF
!     tan_vec(1) = Radius(2)*n_vec(3) - Radius(3)*n_vec(2)
!     tan_vec(2) = Radius(3)*n_vec(1) - Radius(1)*n_vec(3)
!     tan_vec(3) = Radius(1)*n_vec(2) - Radius(2)*n_vec(1)
!     PartState(1,4:6) = tan_vec(1:3) * Species(1)%VeloIC
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
          IPWRITE(UNIT_stdOut,*)"Error in setparticle velocity, gyrotron circle. &
                    & Radius too big after intersection."
        END IF
        !  Normal Vector of circle
        n_vec(1:3) = Species(FractNbr)%Init(iInit)%NormalIC(1:3)

        ! If we're doing Leapfrog, then use velocities from half-timestep before. This only applies in 
        ! x- and y-direction. z has allways same velo. 
!           IF (ParticlePushMethod.EQ.'boris_leap_frog_scheme') THEN
!             ! get angle of part on gyrocircle
!             Angle = ACOS(Radius(1)/Species(1)%RadiusICGyro)
!             IF (Radius(2).LE.0) THEN
!               Angle = 2*PI-Angle
!             END IF
!             ! shift position angle half dt back in time (as particle moves clockwise,
!             ! we add dalpha in ccw direction)
!             Angle = Angle + PIC%GyrationFrequency * dt * 0.5 * PIC%GyroVecDirSIGN  
!             Radius(1) = cos(Angle)
!             Radius(2) = sin(Angle)                 
!           END IF
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
             SWRITE(*,'(A,3(E21.14,X))') 'Velocity=', PartState(PositionNbr,4:6)
             CALL abort(__STAMP__,&
                'ERROR in gyrotron_circle spaceIC!')
           END If
!           IF(SQRT(DOT_PRODUCT(PartState(PositionNbr,4:6),PartState(PositionNbr,4:6))).GT.c)THEN
!             CALL abort(__STAMP__,&
!                'particle faster than speed of light!')
!           END IF
           IF (PartState(PositionNbr,4).NE.PartState(PositionNbr,4) .OR. &
               PartState(PositionNbr,5).NE.PartState(PositionNbr,5) .OR. &
               PartState(PositionNbr,6).NE.PartState(PositionNbr,6)     ) THEN
             SWRITE(*,'(A,3(E21.14,X))') 'WARNING:! NaN: Velocity=', PartState(PositionNbr,4:6)
           END If
        END IF
        i = i + 1
     END DO
CASE('maxwell_lpn')
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr .NE. 0) THEN
       IF (useVTKFileBGG .AND. Is_BGGas) THEN
         CALL CalcVelocity_maxwell_lpn(FractNbr, iInit, Vec3D, PEM%Element(PositionNbr))
       ELSE
         CALL CalcVelocity_maxwell_lpn(FractNbr, iInit, Vec3D)
       END IF
       PartState(PositionNbr,4:6) = Vec3D(1:3)
    END IF
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
  maxwellfac = SQRT(3. * BoltzmannConst * Species(FractNbr)%Init(iInit)%MWTemperatureIC / &              ! velocity of maximum
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

CASE('weibel')
  v_sum(:)  = 0.0
  sigma(:) = 0.0
  
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr .NE. 0) THEN
      Velosq = 2.
      DO WHILE ((Velosq .GT. 1.) .OR. (Velosq .EQ. 0.))
        CALL RANDOM_NUMBER(RandVal)     
        RandVal(:) = 2*RandVal(:)-1
        Velosq = RandVal(1)**2 + RandVal(2)**2
      END DO
      Velosq = sqrt(-2*log(Velosq)/Velosq)
      Vec3D(1) = RandVal(1)*Velosq
      Vec3D(2) = RandVal(2)*Velosq
      
      Velosq = 2.
      DO WHILE ((Velosq .GT. 1.) .OR. (Velosq .EQ. 0.))
        CALL RANDOM_NUMBER(RandVal)     
        RandVal(:) = 2*RandVal(:)-1
        Velosq = RandVal(1)**2 + RandVal(2)**2
      END DO
      Velosq = sqrt(-2*log(Velosq)/Velosq)
      Vec3D(3) = RandVal(1)*Velosq
      
      v_sum(:) = v_sum(:)  + Vec3D(:)
      sigma(:)   = sigma(:)    + Vec3D(:)**2
      PartState(PositionNbr,4:6) = Vec3D(1:3)
    END IF
  END DO
!  WRITE(*,*) PartVeloX(1), PartVeloY(1), PartVeloZ(1)

!    IF (NbrOfParticle .GT. 0) THEN
!      v_sum(:)  = 0.0
!      sigma(:) = 0.0

!      DO i=1,NbrOfParticle
!  !       elemNbr = PartToElem%Element(i)
!  v_sum(1)  = v_sum(1)  + PartVeloX(i)
!  v_sum(2)  = v_sum(2)  + PartVeloY(i)
!  v_sum(3)  = v_sum(3)  + PartVeloZ(i)
!  sigma(1) = sigma(1) + PartVeloX(i)**2
!  sigma(2) = sigma(2) + PartVeloY(i)**2
!  sigma(3) = sigma(3) + PartVeloZ(i)**2
!  !       NPart(elemNbr) = NPart(elemNbr) + 1
!      END DO

  
  IF (NbrOfParticle .GT. 1) THEN
    v_sum(:)  = v_sum(:)/NbrOfParticle
    sigma(:) = (NbrOfParticle/(NbrOfParticle-1))*(sigma(:)/NbrOfParticle-v_sum(:)**2)  
        ! Verschiebungssatz der korrigierten Stichprobenkovarianz:
  ELSE                                                                            ! s^2(X)=1/(N-1)(N*E(X^2) - N*E(X)^2)   
    v_sum(:) = 0.
    sigma(:) = 1.
  END IF
                                                                                      
  ftl = 0
  DO i=1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr .NE. 0) THEN
      PartState(PositionNbr,4)   = (PartState(PositionNbr,4)  -v_sum(1)) * &
                                    SQRT(Species(FractNbr)%Init(iInit)%WeibelVeloPar**2/sigma(1)) *c
      PartState(PositionNbr,5:6) = (PartState(PositionNbr,5:6)-v_sum(2:3)) * & 
                                    SQRT(Species(FractNbr)%Init(iInit)%WeibelVeloPer**2/sigma(2:3)) *c
      PartVelo = SQRT(PartState(PositionNbr,4)**2+PartState(PositionNbr,5)**2+PartState(PositionNbr,6)**2)
      
      DO WHILE (PartVelo .GE. c)
        ftl = ftl+1
        IPWRITE(UNIT_stdOut,*) 'Number of Particles FTL:', ftl
! 	      PartState(PositionNbr,4:6) = PartState(PositionNbr,4:6) * ((1.-0.00001)/PartVelo)
        Velosq = 2.
        DO WHILE ((Velosq .GT. 1.) .OR. (Velosq .EQ. 0.))
          CALL RANDOM_NUMBER(RandVal)     
          RandVal(:) = 2*RandVal(:)-1
          Velosq = RandVal(1)**2 + RandVal(2)**2
        END DO
        Velosq = sqrt(-2*log(Velosq)/Velosq)
        Vec3D(1) = RandVal(1)*Velosq
        Vec3D(2) = RandVal(2)*Velosq
        
        Velosq = 2.
        DO WHILE ((Velosq .GT. 1.) .OR. (Velosq .EQ. 0.))
          CALL RANDOM_NUMBER(RandVal)     
          RandVal(:) = 2*RandVal(:)-1
          Velosq = RandVal(1)**2 + RandVal(2)**2
        END DO
        Velosq = sqrt(-2*log(Velosq)/Velosq)
        Vec3D(3) = RandVal(1)*Velosq
        PartState(PositionNbr,4:6) = Vec3D(1:3)
        
        PartState(PositionNbr,4)   = (Vec3D(1)  -v_sum(1)) * &
            SQRT(Species(FractNbr)%Init(iInit)%WeibelVeloPar**2/sigma(1)) *c
        PartState(PositionNbr,5:6) = (Vec3D(2:3)-v_sum(2:3)) * & 
            SQRT(Species(FractNbr)%Init(iInit)%WeibelVeloPer**2/sigma(2:3)) *c
        PartVelo = SQRT(PartState(PositionNbr,4)**2+PartState(PositionNbr,5)**2+PartState(PositionNbr,6)**2)
      END DO
    END IF
  END DO

CASE('OneD-twostreaminstabilty')
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr .NE. 0) THEN
      PartState(PositionNbr,4) = Species(FractNbr)%Init(iInit)%OneDTwoStreamVelo
      PartState(PositionNbr,5:6) = Species(FractNbr)%Init(iInit)%OneDTwoStreamTransRatio
      
      Velosq = 2.
      DO WHILE ((Velosq .GT. 1.) .OR. (Velosq .EQ. 0.))
        CALL RANDOM_NUMBER(RandVal)     
        RandVal(:) = 2*RandVal(:)-1
        Velosq = RandVal(1)**2 + RandVal(2)**2
      END DO
      Velosq = sqrt(-2*log(Velosq)/Velosq)
      PartState(PositionNbr,5) = RandVal(1)*Velosq*Species(FractNbr)%Init(iInit)%OneDTwoStreamTransRatio* &
                                                   Species(FractNbr)%Init(iInit)%OneDTwoStreamVelo
      PartState(PositionNbr,6) = RandVal(2)*Velosq*Species(FractNbr)%Init(iInit)%OneDTwoStreamTransRatio* &
                                                   Species(FractNbr)%Init(iInit)%OneDTwoStreamVelo
    END IF  
  END DO

END SELECT

END SUBROUTINE SetParticleVelocity


SUBROUTINE SetParticleChargeAndMass(FractNbr,NbrOfParticle)                                        
!===================================================================================================================================
! And partilces mass and charge
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars,    ONLY : PDM, PartSpecies
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE                                                                                    
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                       :: FractNbr                                                     
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                    :: NbrOfParticle                                                
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                  :: i,PositionNbr                                                
!===================================================================================================================================

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

END SUBROUTINE SetParticleChargeAndMass

SUBROUTINE SetParticleMPF(FractNbr,NbrOfParticle) 
!===================================================================================================================================
! finally, set the MPF
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars,    ONLY : PDM, PartMPF, Species
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE                                          
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                       :: FractNbr    
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                    :: NbrOfParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                  :: i,PositionNbr 
!===================================================================================================================================

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

END SUBROUTINE SetParticleMPF

SUBROUTINE ParticleInsertingCellPressure(iSpec,iInit,NbrOfParticle)
!===================================================================================================================================
! Insert constant cell pressure particles (and remove additionals)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars
USE MOD_Mesh_Vars,          ONLY:NGeo,XCL_NGeo
USE MOD_Particle_Mesh_Vars, ONLY:GEO
USE MOD_Particle_Surfaces_vars,ONLY:DoRefMapping
USE MOD_Particle_Mesh,         ONLY:SingleParticleToExactElement,SingleParticleToExactElementNoMap
USE MOD_Eval_xyz,           ONLY:Eval_XYZ_Poly
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
INTEGER               :: iElem, Elem, iPart, i, NbrPartsInCell, iNode, NbrNewParts
INTEGER               :: ParticleIndexNbr
REAL                  :: PartDiff, PartDiffRest, RandVal, RandVal3(1:3), P(1:3,1:8)
INTEGER, ALLOCATABLE  :: PartsInCell(:)
LOGICAL               :: InElementCheck
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
        CALL Eval_xyz_Poly(RandVal3,3,NGeo,XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,iElem),PartState(ParticleIndexNbr,1:3))
        !PartState(ParticleIndexNbr, 1:3) = MapToGeo(RandVal3,P)
        PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
        IF (.NOT. DoRefMapping) THEN
          CALL SingleParticleToExactElementNoMap(ParticleIndexNbr,doHALO=.FALSE.)
        ELSE
          PartPosRef(1:3,ParticleIndexNbr)=RandVal3
        END IF
        IF(.NOT.PDM%ParticleInside(ParticleIndexNbr))THEN
          CALL abort(__STAMP__,&
           ' Particle lost in own MPI region. Need to communicate!')
        END IF
      ELSE
        CALL abort(__STAMP__,&
             'ERROR in ParticleInsertingCellPressure: ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?')
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
USE MOD_Mesh_Vars,              ONLY:NGeo,XCL_NGeo
USE MOD_Particle_Mesh,          ONLY:SingleParticleToExactElement,SingleParticleToExactElementNoMap
USE MOD_Particle_Surfaces_vars, ONLY:DoRefMapping
USE MOD_Eval_xyz,               ONLY:Eval_XYZ_Poly
USE MOD_DSMC_Vars,              ONLY:useDSMC, CollisMode
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
INTEGER                       :: iElem, Elem, iPart, i, NbrPartsInCell, iNode, distnum
INTEGER                       :: ParticleIndexNbr
REAL                          :: RandVal3(1:3), P(1:3,1:8)
INTEGER, ALLOCATABLE          :: PartsInCell(:)
LOGICAL                       :: InElementCheck
REAL                          :: Velo1, Velo2, Velosq, v_sum(3), v2_sum, maxwellfac
REAL                          :: Vec3D(3), RandVal3D(3)
!===================================================================================================================================

NbrOfParticle = 0
IF (CollisMode.EQ.0) THEN
  CALL Abort(&
    __STAMP__,&
    "Free Molecular Flow (CollisMode=0) is not supported for const pressure outflow BC!")
END IF
IF (TRIM(Species(iSpec)%Init(iInit)%velocityDistribution).NE.'maxwell') THEN
  CALL abort(__STAMP__,&
                'Only maxwell implemented yet for const pressure outflow BC!')
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
    !! Build array for element nodes
    !DO iNode = 1,8
    !  P(1:3,iNode) = GEO%NodeCoords(1:3,GEO%ElemToNodeID(iNode,Elem))
    !END DO
    ! insert particles (positions and velocities)
    v_sum(1:3) = 0.0
    v2_sum = 0.0
    DO i = 1, NbrPartsInCell
      ! set random position in -1,1 space
      CALL RANDOM_NUMBER(RandVal3)
      RandVal3 = RandVal3 * 2.0 - 1.0 
      ParticleIndexNbr = PDM%nextFreePosition(PDM%CurrentNextFreePosition + i + NbrOfParticle)
      IF (ParticleIndexNbr.NE.0) THEN
        CALL Eval_xyz_Poly(RandVal3,3,NGeo,XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,iElem),PartState(ParticleIndexNbr,1:3))
        PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
        IF (.NOT. DoRefMapping) THEN
          CALL SingleParticleToExactElementNoMap(ParticleIndexNbr,doHALO=.FALSE.)
        ELSE
          PartPosRef(1:3,ParticleIndexNbr)=RandVal3
        END IF
        IF(.NOT.PDM%ParticleInside(ParticleIndexNbr))THEN
          CALL abort(__STAMP__,&
           ' Particle lost in own MPI region. Need to communicate!')
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
        CALL abort(__STAMP__,&
                'ERROR in ParticleInsertingCellPressureOut: ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?')
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
USE MOD_Globals
USE MOD_Particle_Vars,         ONLY : PartState,usevMPF,Species,PartSpecies,BoltzmannConst,usevMPF,PartMPF,PDM
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
REAL                  :: MPFSum, WeightFak, kappa_part, AvogadroConst, RandVal, RealnumberNewParts, RealnumberNewPartsRest
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
    CALL abort(__STAMP__,&
            'Wrong PartSpecies for outflow BC!')
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
  CALL abort(__STAMP__,&
      'YOU NEED MORE PARTICLES INSIDE THE OUTFLOW REGION!!!')
END IF

END SUBROUTINE ParticleInsertingPressureOut_Sampling


SUBROUTINE CalcVelocity_maxwell_lpn(FractNbr, iInit, Vec3D, Element)
!===================================================================================================================================
! Subroutine to sample current cell values (partly copied from 'LD_DSMC_Mean_Bufferzone_A_Val' and 'dsmc_analyze')
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars,          ONLY : BoltzmannConst, Species, BGGdataAtElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: FractNbr, iInit
INTEGER, OPTIONAL                :: Element !for BGG from VTK
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                 :: Vec3D(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: RandVal(3), Velo1, Velo2, Velosq, Tx, ty, Tz, v_drift(3)
!===================================================================================================================================

IF (PRESENT(Element)) THEN
  Tx=BGGdataAtElem(1,Element)
  Ty=BGGdataAtElem(2,Element)
  Tz=BGGdataAtElem(3,Element)
  v_drift=BGGdataAtElem(4:6,Element)
ELSE
  Tx=Species(FractNbr)%Init(iInit)%MWTemperatureIC
  Ty=Species(FractNbr)%Init(iInit)%MWTemperatureIC
  Tz=Species(FractNbr)%Init(iInit)%MWTemperatureIC
  v_drift=Species(FractNbr)%Init(iInit)%VeloIC *Species(FractNbr)%Init(iInit)%VeloVecIC(1:3)
END IF

Velosq = 2
DO WHILE ((Velosq .GE. 1.) .OR. (Velosq .EQ. 0.))
  CALL RANDOM_NUMBER(RandVal)
  Velo1 = 2.*RandVal(1) - 1.
  Velo2 = 2.*RandVal(2) - 1.
  Velosq = Velo1**2 + Velo2**2
END DO
Vec3D(1) = Velo1*SQRT(-2*BoltzmannConst*Tx/ &
  Species(FractNbr)%MassIC*LOG(Velosq)/Velosq)                                !x-Komponente
Vec3D(2) = Velo2*SQRT(-2*BoltzmannConst*Ty/ &
  Species(FractNbr)%MassIC*LOG(Velosq)/Velosq)                                !y-Komponente
Velosq = 2
DO WHILE ((Velosq .GE. 1.) .OR. (Velosq .EQ. 0.))
  CALL RANDOM_NUMBER(RandVal)
  Velo1 = 2.*RandVal(1) - 1.
  Velo2 = 2.*RandVal(2) - 1.
  Velosq = Velo1**2 + Velo2**2
END DO
Vec3D(3) = Velo1*SQRT(-2*BoltzmannConst*Tz/ &
  Species(FractNbr)%MassIC*LOG(Velosq)/Velosq)                                !z-Komponente
Vec3D(1:3) = Vec3D(1:3) + v_drift

END SUBROUTINE CalcVelocity_maxwell_lpn


SUBROUTINE InsideExcludeRegionCheck(FractNbr, iInit, Particle_pos, insideExcludeRegion)
!===================================================================================================================================
! Subroutine for checking if calculated particle position would be inside user-defined ExcludeRegion (cuboid or cylinder)
!===================================================================================================================================
! MODULES
USE MOD_Globals,                ONLY : abort
USE MOD_Particle_Vars,          ONLY : Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: FractNbr, iInit
REAL,INTENT(IN)                  :: Particle_pos(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)              :: insideExcludeRegion
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: VecExclude(3), DistExclude
INTEGER                          :: iExclude
!===================================================================================================================================

insideExcludeRegion=.FALSE.
DO iExclude=1,Species(FractNbr)%Init(iInit)%NumberOfExcludeRegions
  VecExclude = Particle_pos - Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BasePointIC
  SELECT CASE (TRIM(Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC))
  CASE ('cuboid')
    !--check normal direction
    DistExclude = VecExclude(1)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1) &
      + VecExclude(2)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2) &
      + VecExclude(3)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3)
    IF ( (DistExclude .LE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%CuboidHeightIC) &
      .AND. (DistExclude .GE. 0.) ) THEN
      insideExcludeRegion = .TRUE.
    ELSE
      insideExcludeRegion = .FALSE.
      CYCLE
    END IF
    !--check BV1 direction
    DistExclude = VecExclude(1)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(1) &
      + VecExclude(2)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(2) &
      + VecExclude(3)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(3)
    IF ( (DistExclude .LE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%ExcludeBV_lenghts(1)**2) &
      .AND. (DistExclude .GE. 0.) ) THEN
      insideExcludeRegion = .TRUE.
    ELSE
      insideExcludeRegion = .FALSE.
      CYCLE
    END IF
    !--check BV2 direction
    DistExclude = VecExclude(1)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(1) &
      + VecExclude(2)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(2) &
      + VecExclude(3)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(3)
    IF ( (DistExclude .LE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%ExcludeBV_lenghts(2)**2) &
      .AND. (DistExclude .GE. 0.) ) THEN
      insideExcludeRegion = .TRUE.
      RETURN !particle is inside current ExcludeRegion based an all dimensions
    ELSE
      insideExcludeRegion = .FALSE.
      CYCLE
    END IF
  CASE ('cylinder')
    !--check normal direction
    DistExclude = VecExclude(1)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1) &
      + VecExclude(2)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2) &
      + VecExclude(3)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3)
    IF ( (DistExclude .LE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%CylinderHeightIC) &
      .AND. (DistExclude .GE. 0.) ) THEN
      insideExcludeRegion = .TRUE.
    ELSE
      insideExcludeRegion = .FALSE.
      CYCLE
    END IF
    !--check radial direction
    DistExclude = SQRT( VecExclude(1)**2 + VecExclude(2)**2 + VecExclude(3)**2 - DistExclude**2 )
    IF ( (DistExclude .LE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%RadiusIC) &
      .AND. (DistExclude .GE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%Radius2IC) ) THEN
      insideExcludeRegion = .TRUE.
      RETURN !particle is inside current ExcludeRegion based an all dimensions
    ELSE
      insideExcludeRegion = .FALSE.
      CYCLE
    END IF
  CASE DEFAULT
    CALL abort(__STAMP__,&
      'wrong SpaceIC for ExcludeRegion!')
  END SELECT
END DO

END SUBROUTINE InsideExcludeRegionCheck

END MODULE MOD_part_emission
