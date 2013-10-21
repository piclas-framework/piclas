MODULE MOD_DSMC_BGGas
!===================================================================================================================================
! module including collisions
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

PUBLIC :: DSMC_InitBGGas, DSMC_pairing_bggas, DSMC_FinalizeBGGas
!-----------------------------------------------------------------------------------------------------------------------------------

CONTAINS

SUBROUTINE DSMC_InitBGGas()

USE MOD_DSMC_Vars,          ONLY : BGGas
USE MOD_PARTICLE_Vars,      ONLY : PDM, PartSpecies, PartState, PEM, usevMPF
USE MOD_part_emission,      ONLY : SetParticleChargeAndMass, SetParticleVelocity, SetParticleMPF
USE MOD_part_tools,         ONLY : UpdateNextFreePosition
!--------------------------------------------------------------------------------------------------!
! init BG Gas and Build Pairs
!--------------------------------------------------------------------------------------------------!
IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !  
INTEGER           :: iNewPart, iPart, PositionNbr
!--------------------------------------------------------------------------------------------------!
  iNewPart=0
  PositionNbr = 0
  DO iPart = 1, PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart).AND.(PartSpecies(iPart).NE.BGGas%BGGasSpecies)) THEN
      iNewPart = iNewPart + 1
      PositionNbr = PDM%nextFreePosition(iNewPart+PDM%CurrentNextFreePosition)
      IF (PositionNbr.EQ.0) THEN
        WRITE(*,*) 'ERROR in BGGas: too many Particles'
        STOP
      END IF
      PartState(PositionNbr,1:3) = PartState(iPart,1:3)
      PartSpecies(PositionNbr) = BGGas%BGGasSpecies
      PEM%Element(PositionNbr) = PEM%Element(iPart)
      PDM%ParticleInside(PositionNbr) = .true.
      PEM%pNext(PEM%pEnd(PEM%Element(PositionNbr))) = PositionNbr ! Next Particle of same Elem (Linked List)
      PEM%pEnd(PEM%Element(PositionNbr)) = PositionNbr
      PEM%pNumber(PEM%Element(PositionNbr)) = &                      ! Number of Particles in Element
      PEM%pNumber(PEM%Element(PositionNbr)) + 1
    END IF
  END DO
  CALL SetParticleVelocity(BGGas%BGGasSpecies, iNewPart)
  !IF (usevMPF) CALL SetParticleMPF(BGGas%BGGasSpecies, iNewPart)
  PDM%ParticleVecLength = MAX(PDM%ParticleVecLength,PositionNbr)
  PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + iNewPart 

  !CALL UpdateNextFreePosition()
END SUBROUTINE DSMC_InitBGGas

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE DSMC_pairing_bggas(iElem)

  USE MOD_DSMC_Vars,              ONLY : Coll_pData, CollInf, BGGas, CollisMode, ChemReac, PartStateIntEn
  USE MOD_Particle_Vars,          ONLY : PEM, PartSpecies, nSpecies, PartState, GEO, Species, usevMPF, PartMPF    
!--------------------------------------------------------------------------------------------------!
! statistical pairing method
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
INTEGER                       :: nPair, iPair, iPart, iLoop, cPart1, cPart2, nPart
INTEGER                       :: cSpec1, cSpec2, iCase
INTEGER, ALLOCATABLE          :: iPartIndx(:) ! List of particles in the cell nec for stat pairing
REAL                          :: iRan
! input variable declaration                                                                       !
INTEGER, INTENT(IN)           :: iElem

!--------------------------------------------------------------------------------------------------!
  nPart = PEM%pNumber(iElem)
  nPair = INT(nPart/2)
  
  CollInf%Coll_SpecPartNum = 0
  CollInf%Coll_CaseNum = 0

  ALLOCATE(Coll_pData(nPair))
  ALLOCATE(iPartIndx(nPart))
  Coll_pData%Ec=0
  iPartIndx = 0
  IF ((CollisMode.EQ.3).AND.ChemReac%MeanEVib_Necc) ChemReac%MeanEVibQua_PerIter(1:nSpecies) = 0.0

  iPart = PEM%pStart(iElem)                         ! create particle index list for pairing
  DO iLoop = 1, nPart
    iPartIndx(iLoop) = iPart
    CollInf%Coll_SpecPartNum(PartSpecies(iPart)) = CollInf%Coll_SpecPartNum(PartSpecies(iPart)) + 1 
          ! counter for part num of spec per cell   
    IF ((CollisMode.EQ.3).AND.ChemReac%MeanEVib_Necc) ChemReac%MeanEVibQua_PerIter(PartSpecies(iPart)) = &
                        ChemReac%MeanEVibQua_PerIter(PartSpecies(iPart)) &
                        + PartStateIntEn(iPart,1) !Calulation of mean evib per cell and iter, necessary for disso prob 
    iPart = PEM%pNext(iPart)    
  END DO
  
  ! Setting NUmber of BGGas Particles per Cell
  BGGas%BGColl_SpecPartNum = BGGas%BGGasDensity * GEO%Volume(iElem)      &
                                               / Species(BGGas%BGGasSpecies)%MacroParticleFactor
  

  DO iPair = 1, nPair 
    DO cPart1 = 1, nPart ! Searching the first non BG particle
      IF (PartSpecies(iPartIndx(cPart1)).NE.BGGas%BGGasSpecies) THEN
        Coll_pData(iPair)%iPart_p1 = iPartIndx(cPart1)
        iPartIndx(cPart1) = iPartIndx(nPart)
        EXIT
      END IF
    END DO
    nPart = nPart - 1

    DO cPart2 = 1, nPart ! Searching the first BG particle
      IF (PartSpecies(iPartIndx(cPart2)).EQ.BGGas%BGGasSpecies) THEN
        Coll_pData(iPair)%iPart_p2 = iPartIndx(cPart2)
        iPartIndx(cPart2) = iPartIndx(nPart)
        EXIT
      END IF
    END DO
    nPart = nPart - 1

    cSpec1 = PartSpecies(Coll_pData(iPair)%iPart_p1) !spec of particle 1   
    cSpec2 = PartSpecies(Coll_pData(iPair)%iPart_p2) !spec of particle 2
    IF (usevMPF) PartMPF(Coll_pData(iPair)%iPart_p2) = PartMPF(Coll_pData(iPair)%iPart_p1) 
    iCase = CollInf%Coll_Case(cSpec1, cSpec2) 
    CollInf%Coll_CaseNum(iCase) = CollInf%Coll_CaseNum(iCase) + 1 !sum of coll case (Sab)
    Coll_pData(iPair)%CRela2 = (PartState(Coll_pData(iPair)%iPart_p1,4) &
                             -  PartState(Coll_pData(iPair)%iPart_p2,4))**2 &
                             + (PartState(Coll_pData(iPair)%iPart_p1,5) &
                             -  PartState(Coll_pData(iPair)%iPart_p2,5))**2 &
                             + (PartState(Coll_pData(iPair)%iPart_p1,6) &
                             -  PartState(Coll_pData(iPair)%iPart_p2,6))**2 
    Coll_pData(iPair)%PairType = iCase
  END DO
  
  DEALLOCATE(iPartIndx)
END SUBROUTINE DSMC_pairing_bggas

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE DSMC_FinalizeBGGas()

USE MOD_DSMC_Vars,          ONLY : BGGas
USE MOD_PARTICLE_Vars,      ONLY : PDM, PartSpecies
USE MOD_part_tools,         ONLY : UpdateNextFreePosition
!--------------------------------------------------------------------------------------------------!
! kills all BG Particles
!--------------------------------------------------------------------------------------------------!
IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !  
!--------------------------------------------------------------------------------------------------!
  WHERE (PartSpecies(1:PDM%ParticleVecLength).EQ.BGGas%BGGasSpecies) &
         PDM%ParticleInside(1:PDM%ParticleVecLength) = .FALSE.
  CALL UpdateNextFreePosition()

END SUBROUTINE DSMC_FinalizeBGGas

!--------------------------------------------------------------------------------------------------!
END MODULE MOD_DSMC_BGGas
