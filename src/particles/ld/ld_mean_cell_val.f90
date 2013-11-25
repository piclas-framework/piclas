MODULE MOD_LD_mean_cell

!===================================================================================================================================
! module for determination of cell quantities
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

PUBLIC :: CalcMacCellLDValues
!===================================================================================================================================

CONTAINS

!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE CalcMacCellLDValues()

USE MOD_LD_Vars
USE MOD_Mesh_Vars,              ONLY : nElems, ElemToSide, SideToElem
USE MOD_Particle_Vars,          ONLY : GEO, PEM, usevMPF, PartMPF, BoltzmannConst, Species, PartSpecies, PDM

!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
  INTEGER             :: iElem, iPart, nPart, iPartIndx, LowPartCount
  REAL                :: MPFSum, WeightFak
  REAL                :: CellMass, CellTemp, CellTempMean, CellVelo2, StanVol
  INTEGER             :: NeiElem, nNei
  INTEGER             :: iLocSide, SideID
!--------------------------------------------------------------------------------------------------!

LowPartCount = 0
DO iElem = 1, nElems ! element=cell main loop
  nPart = PEM%pNumber(iElem)
  IF (nPart.GT. 1) THEN ! Are there more than one particle
    BulkValues(iElem)%CellV           = 0.0
    BulkValues(iElem)%MassDens        = 0.0
    BulkValues(iElem)%DegreeOfFreedom = 0.0
    CellTempMean                      = 0.0
    CellMass                          = 0.0
    MPFSum                            = 0.0
    CellVelo2                         = 0.0
    iPartIndx = PEM%pStart(iElem)
    DO ipart = 1, nPart
      IF (usevMPF) THEN
         WeightFak = PartMPF(iPartIndx)
      ELSE
         WeightFak = Species(PartSpecies(iPartIndx))%MacroParticleFactor
      END IF
      BulkValues(iElem)%CellV(1)        = BulkValues(iElem)%CellV(1) + PartStateBulkValues(iPartIndx,1) &
                                        * WeightFak * Species(PartSpecies(iPartIndx))%MassIC
      BulkValues(iElem)%CellV(2)        = BulkValues(iElem)%CellV(2) + PartStateBulkValues(iPartIndx,2) &
                                        * WeightFak * Species(PartSpecies(iPartIndx))%MassIC
      BulkValues(iElem)%CellV(3)        = BulkValues(iElem)%CellV(3) + PartStateBulkValues(iPartIndx,3) &
                                        * WeightFak * Species(PartSpecies(iPartIndx))%MassIC
      CellVelo2                         = CellVelo2 + ( (PartStateBulkValues(iPartIndx,1))**2 &
                                                      + (PartStateBulkValues(iPartIndx,2))**2 &
                                                      + (PartStateBulkValues(iPartIndx,3))**2 ) &
                                                      * WeightFak * Species(PartSpecies(iPartIndx))%MassIC
      CellTempMean                      = CellTempMean + PartStateBulkValues(iPartIndx,4) * WeightFak
      IF (CellTempMean.lt. 0) THEN
        PRINT*,iElem, iPartIndx, PDM%ParticleVecLength, PartStateBulkValues(iPartIndx,4)
        read*
      END IF
      BulkValues(iElem)%DegreeOfFreedom = BulkValues(iElem)%DegreeOfFreedom + PartStateBulkValues(iPartIndx,5) * WeightFak
      CellMass                          = CellMass + WeightFak * Species(PartSpecies(iPartIndx))%MassIC
      MPFSum = MPFSum + WeightFak
      iPartIndx = PEM%pNext(iPartIndx)
    END DO

    BulkValues(iElem)%CellV(1) = BulkValues(iElem)%CellV(1) / CellMass
    BulkValues(iElem)%CellV(2) = BulkValues(iElem)%CellV(2) / CellMass
    BulkValues(iElem)%CellV(3) = BulkValues(iElem)%CellV(3) / CellMass
    CellVelo2                  = CellVelo2 / CellMass
    CellTempMean               = CellTempMean / MPFSum
    BulkValues(iElem)%DegreeOfFreedom = BulkValues(iElem)%DegreeOfFreedom / MPFSum
    BulkValues(iElem)%MassDens = CellMass / GEO%Volume(iElem)
    CellTemp = CellTempMean &
             + CellMass / MPFSum /(BulkValues(iElem)%DegreeOfFreedom * BoltzmannConst) & 
             * (nPart/(nPart-1)) &
             * (CellVelo2 - ( BulkValues(iElem)%CellV(1)**2 &
                            + BulkValues(iElem)%CellV(2)**2 &
                            + BulkValues(iElem)%CellV(3)**2 ))
    BulkValues(iElem)%Beta = SQRT(CellMass / MPFSum / (2 * CellTemp * BoltzmannConst))
  ELSE ! if there is no particle => keep old LD-state
    LowPartCount = LowPartCount + 1
  END IF
END DO
! Abfangen von Nulldivisionen
IF (LowPartCount.GE. 1) THEN
  DO iElem = 1, nElems
    nPart = PEM%pNumber(iElem)
    IF (nPart.LE. 1) THEN
      BulkValues(iElem)%CellV           = 0.0
      BulkValues(iElem)%DegreeOfFreedom = 0.0
      BulkValues(iElem)%Beta            = 0.0
      BulkValues(iElem)%MassDens        = 0.0
      nNei = 0
      DO iLocSide=1, 6
        SideID = ElemToSide(1,iLocSide,iElem)
        IF (SideToElem(1,SideID) .EQ. iElem) THEN
          IF (SideToElem(2,SideID).LT. 0) CYCLE
          IF (PEM%pNumber(SideToElem(2,SideID)) .LE. 1) CYCLE
          NeiElem = SideToElem(2,SideID)
          BulkValues(iElem)%CellV(1:3) = BulkValues(iElem)%CellV(1:3) &
                                       + BulkValues(NeiElem)%CellV(1:3)
          BulkValues(iElem)%DegreeOfFreedom = BulkValues(iElem)%DegreeOfFreedom &
                                            + BulkValues(NeiElem)%DegreeOfFreedom
          BulkValues(iElem)%Beta = BulkValues(iElem)%Beta &
                                 + BulkValues(NeiElem)%Beta
          BulkValues(iElem)%MassDens = BulkValues(iElem)%MassDens &
                                     + BulkValues(NeiElem)%MassDens
          nNei = nNei + 1
        ELSE
          IF (SideToElem(1,SideID).LT. 0) CYCLE
          IF (PEM%pNumber(SideToElem(1,SideID)) .LE. 1) CYCLE
          NeiElem = SideToElem(1,SideID)
          BulkValues(iElem)%CellV(1:3) = BulkValues(iElem)%CellV(1:3) &
                                       + BulkValues(NeiElem)%CellV(1:3)
          BulkValues(iElem)%DegreeOfFreedom = BulkValues(iElem)%DegreeOfFreedom &
                                            + BulkValues(NeiElem)%DegreeOfFreedom
          BulkValues(iElem)%Beta = BulkValues(iElem)%Beta &
                                 + BulkValues(NeiElem)%Beta
          BulkValues(iElem)%MassDens = BulkValues(iElem)%MassDens &
                                     + BulkValues(NeiElem)%MassDens
          nNei = nNei + 1
        END IF
      END DO
      IF (nNei.EQ.0) THEN
        PRINT*,'YOU NEED MORE PARTCLES FOR LD!!!'
        PRINT*,iElem, nNei, BulkValues(iElem)%Beta
        STOP
      END IF
      BulkValues(iElem)%CellV(1:3) = BulkValues(iElem)%CellV(1:3) / nNei
      BulkValues(iElem)%DegreeOfFreedom = BulkValues(iElem)%DegreeOfFreedom / nNei
      BulkValues(iElem)%Beta = BulkValues(iElem)%Beta / nNei   
    END IF     
  END DO
END IF

END SUBROUTINE CalcMacCellLDValues

!=================================================================================================================================

END MODULE MOD_LD_mean_cell
