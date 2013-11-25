MODULE MOD_LD_reassign_part_prop

!===================================================================================================================================
! module for determination of Lagrangian cell velocity
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

PUBLIC :: LD_reassign_prop, UpdateMacLDValues
!===================================================================================================================================

CONTAINS

!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE LD_reassign_prop(iElem)

USE MOD_LD_Vars
USE MOD_Mesh_Vars,             ONLY : nElems, nSides, SideToElem, ElemToSide
USE MOD_TimeDisc_Vars,         ONLY : dt
USE MOD_Mesh_Vars,             ONLY : ElemToSide
USE MOD_Particle_Vars,         ONLY : GEO, PEM

!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                     !
! Local variable declaration                                                                       !
  INTEGER           :: iLocSide, trinum
  REAL              :: Velo(3)
  REAL              :: Beta, Dens, vLAG
  REAL              :: VeloDiff
  REAL              :: NVec(3)
  REAL              :: kon, VeloDir
  REAL              :: Phi
  REAL, PARAMETER   :: PI=3.14159265358979323846_8
  REAL              :: DeltaM(3) 
  REAL              :: DeltaE
  REAL              :: Area

  INTEGER, INTENT(IN)           :: iElem
!--------------------------------------------------------------------------------------------------!

  Velo = BulkValues(iElem)%CellV
  Beta = BulkValues(iElem)%Beta
  Dens = BulkValues(iElem)%MassDens
  DeltaM(1:3) = 0.0
  DeltaE = 0.0
  DO iLocSide=1, 6
    DO trinum = 1, 2
      Area = SurfLagValues(iLocSide,iElem,trinum)%Area
      NVec = SurfLagValues(iLocSide,iElem,trinum)%LagNormVec
      vLAG = SurfLagValues(iLocSide,iElem,trinum)%LagVelo
      kon = Dens / (2 * SQRT(PI) * Beta**2)
      VeloDir = Velo(1) * NVec(1) &
              + Velo(2) * NVec(2) &  
              + Velo(3) * NVec(3)
      VeloDiff =  Beta * ( VeloDir - vLAG )
      Phi  = kon * ( VeloDiff * EXP(-VeloDiff**2) + SQRT(PI) &
           * ( 1 + ERF(VeloDiff) ) * ( 0.5 + VeloDiff**2 ) ) 
      DeltaM(1) = DeltaM(1) + (-2*Area*dt*Phi*NVec(1))
      DeltaM(2) = DeltaM(2) + (-2*Area*dt*Phi*NVec(2))
      DeltaM(3) = DeltaM(3) + (-2*Area*dt*Phi*NVec(3))
      DeltaE    = DeltaE + (-2*Area*vLAG*dt*Phi)
    END DO      ! END loop over trinum
  END DO        ! END loop over iLocSides

  IF (ABS(DeltaM(1)).LT.1E-14) DeltaM(1) = 0.0
  IF (ABS(DeltaM(2)).LT.1E-14) DeltaM(2) = 0.0
  IF (ABS(DeltaM(3)).LT.1E-14) DeltaM(3) = 0.0
  IF (ABS(DeltaE).LT.1E-14) DeltaE = 0.0
  CALL UpdateMacLDValues(iElem, DeltaM, DeltaE)

END SUBROUTINE LD_reassign_prop


!--------------------------------------------------------------------------------------------------!


SUBROUTINE UpdateMacLDValues(iElem, DeltaM, DeltaE)

USE MOD_LD_Vars
USE MOD_Mesh_Vars,             ONLY : nElems, nSides, SideToElem
USE MOD_Particle_Vars,         ONLY : GEO, BoltzmannConst, PEM
USE MOD_LD_Init,               ONLY : CalcDegreeOfFreedom

!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
  REAL                          :: CellTemp, CellTempNew, CellPartDens
  REAL                          :: CellV_2, CellV_old_2
  INTEGER                       :: iPartIndx, nPart, iPart

  INTEGER, INTENT(IN)           :: iElem
  REAL, INTENT(IN)              :: DeltaE
  REAL, INTENT(IN)              :: DeltaM(3)
!--------------------------------------------------------------------------------------------------!

  CellV_old_2 = BulkValues(iElem)%CellV(1)**2 &
              + BulkValues(iElem)%CellV(2)**2 &
              + BulkValues(iElem)%CellV(3)**2
  BulkValues(iElem)%CellV(1) = BulkValues(iElem)%CellV(1) + DeltaM(1) &
                             / (BulkValues(iElem)%MassDens * GEO%Volume(iElem))
  BulkValues(iElem)%CellV(2) = BulkValues(iElem)%CellV(2) + DeltaM(2) &
                             / (BulkValues(iElem)%MassDens * GEO%Volume(iElem))
  BulkValues(iElem)%CellV(3) = BulkValues(iElem)%CellV(3) + DeltaM(3) &
                             / (BulkValues(iElem)%MassDens * GEO%Volume(iElem))
  IF (ABS(BulkValues(iElem)%CellV(1)).LE. 1E-14) THEN
    BulkValues(iElem)%CellV(1) = 0.0
  END IF
  IF (ABS(BulkValues(iElem)%CellV(2)).LE. 1E-14) THEN
    BulkValues(iElem)%CellV(2) = 0.0
  END IF
  IF (ABS(BulkValues(iElem)%CellV(3)).LE. 1E-14) THEN
    BulkValues(iElem)%CellV(3) = 0.0
  END IF
  CellV_2 = BulkValues(iElem)%CellV(1)**2 &
          + BulkValues(iElem)%CellV(2)**2 &
          + BulkValues(iElem)%CellV(3)**2
  CALL CalcCellTemp_PartDens(iElem, CellTemp, CellPartDens)
  CellTempNew = CellTemp &
              + 2 / (BulkValues(iElem)%DegreeOfFreedom * CellPartDens * BoltzmannConst) &
              * ( DeltaE / GEO%Volume(iElem) - 0.5 * BulkValues(iElem)%MassDens &
              * (CellV_2 - CellV_old_2) )
  IF (CellTempNew.lt. 0) then
    PRINT*,'ERROR in Elem, Temperature is lt zero:', iElem
    STOP
  END IF
  !
  ! reassign particle properties
  !
  nPart     = PEM%pNumber(iElem)
  iPartIndx = PEM%pStart(iElem)
  DO ipart = 1, nPart
    PartStateBulkValues(iPartIndx,1) = BulkValues(iElem)%CellV(1)
    PartStateBulkValues(iPartIndx,2) = BulkValues(iElem)%CellV(2)
    PartStateBulkValues(iPartIndx,3) = BulkValues(iElem)%CellV(3)
    PartStateBulkValues(iPartIndx,4) = CellTempNew
    PartStateBulkValues(iPartIndx,5) = CalcDegreeOfFreedom(iPartIndx)
    iPartIndx = PEM%pNext(iPartIndx)
  END DO
END SUBROUTINE UpdateMacLDValues

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE CalcCellTemp_PartDens(iElem, CellTemp, CellPartDens)

  USE MOD_LD_Vars
  USE MOD_Particle_Vars,      ONLY : Species, PartSpecies, BoltzmannConst, usevMPF, PartMPF, PEM, GEO
!--------------------------------------------------------------------------------------------------!
! calculation of LD-cell temperatur
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE 
! LOCAL VARIABLES
!--------------------------------------------------------------------------------------------------!
  REAL                          :: CellMass, WeightFak, MPFSum
  INTEGER                       :: nPart, iPartIndx, iPart
  CHARACTER(LEN=26)             :: myFileName
!--------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!--------------------------------------------------------------------------------------------------!
  INTEGER, INTENT(IN)           :: iElem
  REAL, INTENT(OUT)             :: CellTemp, CellPartDens
!#ifdef MPI
!#endif
!===================================================================================================

  nPart     = PEM%pNumber(iElem)
  iPartIndx = PEM%pStart(iElem)
  CellMass  = 0.0
  MPFSum    = 0.0
  DO ipart = 1, nPart
    IF (usevMPF) THEN
       WeightFak = PartMPF(iPartIndx)
    ELSE
       WeightFak = Species(PartSpecies(iPartIndx))%MacroParticleFactor
    END IF
    CellMass  = CellMass + WeightFak * Species(PartSpecies(iPartIndx))%MassIC
    MPFSum    = MPFSum + WeightFak
    iPartIndx = PEM%pNext(iPartIndx)
  END DO
  CellPartDens = MPFSum / GEO%Volume(iElem)
  CellTemp = CellMass / MPFSum / (2 * BulkValues(iElem)%Beta**2 * BoltzmannConst)

END SUBROUTINE CalcCellTemp_PartDens

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

END MODULE MOD_LD_reassign_part_prop
