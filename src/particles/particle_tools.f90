MODULE MOD_part_tools                                                                              !
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars, ONLY : useDSMC                                                                !
!--------------------------------------------------------------------------------------------------!
  IMPLICIT NONE                                                                                    !
  PRIVATE                                                                                          !
!--------------------------------------------------------------------------------------------------!
  PUBLIC            :: UpdateNextFreePosition                                                      !
!--------------------------------------------------------------------------------------------------!
                                                                                                   !
CONTAINS                                                                                           !
                                                                                                   !
SUBROUTINE UpdateNextFreePosition()                                                                !
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars, ONLY : PDM,PEM, PartSpecies, doParticleMerge, vMPF_SpecNumElem
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE                                                                                    !
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION                                                                        !
! LOCAL VARIABLE DECLARATION                                                                       !
   INTEGER                          :: counter1,i,n                                                !
   REAL                             :: EndT, StartT
!===================================================================================================================================

   IF(PDM%maxParticleNumber.EQ.0) RETURN
   counter1 = 1
   IF (useDSMC.OR.doParticleMerge) THEN
     PEM%pNumber(:) = 0
     !PDM%nextUsedPosition(:) = 0
   END IF
   !PDM%nextFreePosition(:) = 0
   n = PDM%ParticleVecLength !PDM%maxParticleNumber
   PDM%ParticleVecLength = 0
   PDM%insideParticleNumber = 0
   IF (doParticleMerge) vMPF_SpecNumElem = 0
   !PDM%nextFreePosition(1) = 1
   IF (useDSMC.OR.doParticleMerge) THEN
     DO i=1,n
       IF (.NOT.PDM%ParticleInside(i)) THEN
         PDM%nextFreePosition(counter1) = i
         counter1 = counter1 + 1
       ELSE
         !PDM%nextUsedPosition(i+1-counter1) = i
         IF (PEM%pNumber(PEM%Element(i)).EQ.0) THEN
           PEM%pStart(PEM%Element(i)) = i                    ! Start of Linked List for Particles in Elem
         ELSE
           PEM%pNext(PEM%pEnd(PEM%Element(i))) = i ! Next Particle of same Elem (Linked List)
         END IF
         PEM%pEnd(PEM%Element(i)) = i
         PEM%pNumber(PEM%Element(i)) = &                      ! Number of Particles in Element
         PEM%pNumber(PEM%Element(i)) + 1
         PDM%ParticleVecLength = i
         IF(doParticleMerge) vMPF_SpecNumElem(PEM%Element(i),PartSpecies(i)) = vMPF_SpecNumElem(PEM%Element(i),PartSpecies(i)) + 1
       END IF
     END DO
   ELSE ! no DSMC
     DO i=1,n
       IF (.NOT.PDM%ParticleInside(i)) THEN
         PDM%nextFreePosition(counter1) = i
         counter1 = counter1 + 1
       ELSE
         PDM%ParticleVecLength = i
       END IF
     END DO
   ENDIF
   PDM%insideParticleNumber = PDM%ParticleVecLength - counter1
   PDM%CurrentNextFreePosition = 0
   DO i = n+1,PDM%maxParticleNumber
     PDM%nextFreePosition(counter1) = i
     counter1 = counter1 + 1
   END DO 
   PDM%nextFreePosition(counter1:PDM%MaxParticleNumber)=0

 RETURN
END SUBROUTINE UpdateNextFreePosition

END MODULE MOD_part_tools
