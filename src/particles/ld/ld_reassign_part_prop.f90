!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_LD_reassign_part_prop
!===================================================================================================================================
! module for determination of new ld particle bulk values
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE LD_reassign_prop
  MODULE PROCEDURE LD_reassign_prop
END INTERFACE
INTERFACE UpdateMacLDValues
  MODULE PROCEDURE UpdateMacLDValues
END INTERFACE

PUBLIC :: LD_reassign_prop, UpdateMacLDValues
!===================================================================================================================================

CONTAINS

!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE LD_reassign_prop(iElem)
!===================================================================================================================================
! determination of new ld particle bulk values
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars,          ONLY: PI
USE MOD_LD_Vars
USE MOD_TimeDisc_Vars,         ONLY : dt
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
  INTEGER           :: iLocSide, trinum                                                            !
  REAL              :: Velo(3)                                                                     !
  REAL              :: Beta, Dens, vLAG                                                            !
  REAL              :: VeloDiff                                                                    !
  REAL              :: NVec(3)                                                                     !
  REAL              :: kon, VeloDir                                                                !
  REAL              :: Phi                                                                         !
  !REAL, PARAMETER   :: PI=3.14159265358979323846_8                                                 !
  REAL              :: DeltaM(3)                                                                   !
  REAL              :: DeltaE                                                                      !
  REAL              :: Area                                                                        !
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iElem                                                           !
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
      kon = Dens / (2. * SQRT(PI) * Beta**2)
      VeloDir = Velo(1) * NVec(1) &
              + Velo(2) * NVec(2) &
              + Velo(3) * NVec(3)
      VeloDiff =  Beta * ( VeloDir - vLAG )
      Phi  = kon * ( VeloDiff * EXP(-VeloDiff**2) + SQRT(PI) &
           * ( 1. + ERF(VeloDiff) ) * ( 0.5 + VeloDiff**2 ) )
      DeltaM(1) = DeltaM(1) + (-2.*Area*dt*Phi*NVec(1))
      DeltaM(2) = DeltaM(2) + (-2.*Area*dt*Phi*NVec(2))
      DeltaM(3) = DeltaM(3) + (-2.*Area*dt*Phi*NVec(3))
      DeltaE    = DeltaE + (-2.*Area*vLAG*dt*Phi)
      CALL CalcViscousTerms(iElem,iLocSide,trinum,DeltaM,DeltaE)
    END DO      ! END loop over trinum
  END DO        ! END loop over iLocSides
  CALL UpdateMacLDValues(iElem, DeltaM, DeltaE)

END SUBROUTINE LD_reassign_prop

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE UpdateMacLDValues(iElem, DeltaM, DeltaE)
!===================================================================================================================================
! update new ld particle bulk values
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,          ONLY : BoltzmannConst
USE MOD_LD_Vars
USE MOD_Particle_Vars,         ONLY : PEM
USE MOD_LD_Init,               ONLY : CalcDegreeOfFreedom
USE MOD_Particle_Mesh_Vars,    ONLY : GEO
USE MOD_LD_internal_Temp
USE MOD_DSMC_Vars,             ONLY : CollisMode, LD_MultiTemperaturMod
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
  REAL                          :: CellTemp, CellTempNew, CellPartDens
  REAL                          :: CellV_2, CellV_old_2
  REAL                          :: VX_New, VY_New, VZ_New
  INTEGER                       :: iPartIndx, nPart, iPart
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iElem
  REAL, INTENT(IN)              :: DeltaE
  REAL, INTENT(IN)              :: DeltaM(3)
!--------------------------------------------------------------------------------------------------!

  CellV_old_2 = BulkValues(iElem)%CellV(1)**2 &
              + BulkValues(iElem)%CellV(2)**2 &
              + BulkValues(iElem)%CellV(3)**2
  VX_New = BulkValues(iElem)%CellV(1) + DeltaM(1) &
                             / (BulkValues(iElem)%MassDens * GEO%Volume(iElem))
  VY_New = BulkValues(iElem)%CellV(2) + DeltaM(2) &
                             / (BulkValues(iElem)%MassDens * GEO%Volume(iElem))
  VZ_New = BulkValues(iElem)%CellV(3) + DeltaM(3) &
                             / (BulkValues(iElem)%MassDens * GEO%Volume(iElem))
  IF (ABS(VX_New).LE. 1E-14) THEN
    VX_New = 0.0
  END IF
  IF (ABS(VY_New).LE. 1E-14) THEN
    VY_New = 0.0
  END IF
  IF (ABS(VZ_New).LE. 1E-14) THEN
    VZ_New = 0.0
  END IF
  CellV_2 = VX_New**2 &
          + VY_New**2 &
          + VZ_New**2
  CALL CalcCellTemp_PartDens(iElem, CellTemp, CellPartDens)
  CellTempNew = CellTemp &
              + 2. / (BulkValues(iElem)%DegreeOfFreedom * CellPartDens * BoltzmannConst) &
              * ( DeltaE / GEO%Volume(iElem) - 0.5 * BulkValues(iElem)%MassDens &
              * (CellV_2 - CellV_old_2) )
  BulkValues(iElem)%BulkTemperature = CellTempNew
  IF (CellTempNew.lt. 0) then
    SWRITE(UNIT_StdOut,'(132("-"))')
    SWRITE(UNIT_stdOut,'(A)') 'Element, Temperatur:',iElem, CellTempNew
CALL abort(&
__STAMP__&
,'ERROR LD-DSMC: Temperature is lt zero')
  END IF
  !
  ! reassign particle properties
  !
  IF (CollisMode.GT.1) THEN
    IF(LD_MultiTemperaturMod.NE. 1) THEN
      SELECT CASE(LD_MultiTemperaturMod)
        CASE(0)
          ! Do Nothing
        CASE(1)
CALL abort(&
__STAMP__&
,'ERROR: Wrong MultiTemperature Model here!')
        CASE(2)
          CALL CalcInternalTemp_LD_second(iElem)
        CASE(3)
          CALL CalcInternalTemp_LD_third(iElem)
        CASE DEFAULT
CALL abort(&
__STAMP__&
,'ERROR: Wrong MultiTemperature Model!')
      END SELECT
    END IF
  END IF

  IF (BulkValues(iElem)%BulkTemperature.lt. 0) then
    SWRITE(UNIT_StdOut,'(132("-"))')
    SWRITE(UNIT_stdOut,'(A)') 'Element, Temperatur:',iElem, BulkValues(iElem)%BulkTemperature
CALL abort(&
__STAMP__&
,'ERROR LD-DSMC: Temperature is lt zero after MultiTemperature Model')
  END IF

  nPart     = PEM%pNumber(iElem)
  iPartIndx = PEM%pStart(iElem)
  DO ipart = 1, nPart
    PartStateBulkValues(iPartIndx,1) = VX_New
    PartStateBulkValues(iPartIndx,2) = VY_New
    PartStateBulkValues(iPartIndx,3) = VZ_New
    PartStateBulkValues(iPartIndx,4) = BulkValues(iElem)%BulkTemperature
    PartStateBulkValues(iPartIndx,5) = CalcDegreeOfFreedom(iPartIndx)
    iPartIndx = PEM%pNext(iPartIndx)
  END DO
END SUBROUTINE UpdateMacLDValues

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE CalcCellTemp_PartDens(iElem, CellTemp, CellPartDens)
!===================================================================================================================================
! calculation of cell temperature and density
!===================================================================================================================================
! MODULES
  USE MOD_LD_Vars
  USE MOD_Globals_Vars,       ONLY : BoltzmannConst
  USE MOD_Particle_Vars,      ONLY : Species, PartSpecies, usevMPF, PartMPF, PEM
  USE MOD_Particle_Mesh_Vars, ONLY: GEO
!--------------------------------------------------------------------------------------------------!
! calculation of LD-cell temperatur
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE
! LOCAL VARIABLES
!--------------------------------------------------------------------------------------------------!
  REAL                          :: CellMass, WeightFak, MPFSum
  INTEGER                       :: nPart, iPartIndx, iPart
!--------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL, INTENT(OUT)             :: CellTemp, CellPartDens
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
  CellTemp = CellMass / MPFSum / (2. * BulkValues(iElem)%Beta**2 * BoltzmannConst)

END SUBROUTINE CalcCellTemp_PartDens

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE CalcViscousTerms(iElem,iLocSide,trinum,DeltaM,DeltaE)
!===================================================================================================================================
! calculation viscous transport terms
!===================================================================================================================================
! MODULES
USE MOD_Globals,               ONLY : abort
USE MOD_LD_Vars
USE MOD_TimeDisc_Vars,         ONLY : dt
USE MOD_Mesh_Vars,             ONLY : ElemToSide, nBCSides, BC
USE MOD_Particle_Mesh_Vars,    ONLY : SidePeriodicType
USE MOD_Particle_Boundary_Vars, ONLY:PartBound
!--------------------------------------------------------------------------------------------------!
! calculation of LD-cell temperatur
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE
! LOCAL VARIABLES
!--------------------------------------------------------------------------------------------------!
  REAL              :: Area
  REAL              :: BulkVeloDiff(3),MeanBulkVelo(3)
  REAL              :: NVec(3), T1Vec(3), T2Vec(3), CellCenterDiff(3)
  REAL              :: BulkTempDiff, DynamicVisc, ThermalCond, CellCenterDiffDir
  REAL              :: BulkVeloDiffDirN, BulkVeloDiffDirT1, BulkVeloDiffDirT2
  REAL              :: MeanBulkVeloDirN, MeanBulkVeloDirT1, MeanBulkVeloDirT2
  REAL              :: BulkVeloCellT1, BulkVeloCellT2
  REAL              :: TempCell, TempWall, SpecR
  REAL              :: MomentumACC, TransACC
  REAL              :: WallVelo(3)
  INTEGER           :: SideID
  REAL, PARAMETER   :: PI=3.14159265358979323846_8
!--------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iElem,iLocSide,trinum
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL, INTENT(OUT)             :: DeltaM(3),DeltaE
!===================================================================================================

  Area = SurfLagValues(iLocSide,iElem,trinum)%Area
  NVec = SurfLagValues(iLocSide,iElem,trinum)%LagNormVec
  T1Vec(1)       = SurfLagValues(iLocSide, iElem,trinum)%LagTangVec(1,1)
  T1Vec(2)       = SurfLagValues(iLocSide, iElem,trinum)%LagTangVec(1,2)
  T1Vec(3)       = SurfLagValues(iLocSide, iElem,trinum)%LagTangVec(1,3)
  T2Vec(1)       = SurfLagValues(iLocSide, iElem,trinum)%LagTangVec(2,1)
  T2Vec(2)       = SurfLagValues(iLocSide, iElem,trinum)%LagTangVec(2,2)
  T2Vec(3)       = SurfLagValues(iLocSide, iElem,trinum)%LagTangVec(2,3)

  SideID = ElemToSide(1,iLocSide,iElem)
  IF (SideID.GT.nBCSides) THEN
    IF(SidePeriodicType(SideID).EQ.0)THEN
    !IF(GEO%PeriodicElemSide(iLocSide,iElem).EQ.0) THEN ! only inner side without periodic sides
      MeanBulkVelo   = MeanSurfValues(iLocSide, iElem)%MeanBulkVelo
      BulkVeloDiff   = MeanSurfValues(iLocSide, iElem)%BulkVeloDiff
      BulkTempDiff   = MeanSurfValues(iLocSide, iElem)%BulkTempDiff
      CellCenterDiff = MeanSurfValues(iLocSide, iElem)%CellCentDist
      DynamicVisc    = MeanSurfValues(iLocSide, iElem)%DynamicVisc
      ThermalCond    = MeanSurfValues(iLocSide, iElem)%ThermalCond

      BulkVeloDiffDirN = BulkVeloDiff(1) * NVec(1) &
                       + BulkVeloDiff(2) * NVec(2) &
                       + BulkVeloDiff(3) * NVec(3)
      BulkVeloDiffDirT1 = BulkVeloDiff(1) * T1Vec(1) &
                        + BulkVeloDiff(2) * T1Vec(2) &
                        + BulkVeloDiff(3) * T1Vec(3)
      BulkVeloDiffDirT2 = BulkVeloDiff(1) * T2Vec(1) &
                        + BulkVeloDiff(2) * T2Vec(2) &
                        + BulkVeloDiff(3) * T2Vec(3)
      CellCenterDiffDir = CellCenterDiff(1) * NVec(1) &
                        + CellCenterDiff(2) * NVec(2) &
                        + CellCenterDiff(3) * NVec(3)
      MeanBulkVeloDirN = MeanBulkVelo(1) * NVec(1) &
                       + MeanBulkVelo(2) * NVec(2) &
                       + MeanBulkVelo(3) * NVec(3)
      MeanBulkVeloDirT1 = MeanBulkVelo(1) * T1Vec(1) &
                        + MeanBulkVelo(2) * T1Vec(2) &
                        + MeanBulkVelo(3) * T1Vec(3)
      MeanBulkVeloDirT2 = MeanBulkVelo(1) * T2Vec(1) &
                        + MeanBulkVelo(2) * T2Vec(2) &
                        + MeanBulkVelo(3) * T2Vec(3)
      DeltaM(1) = DeltaM(1) + dt * Area * ( &
                  4. / 3. * DynamicVisc * BulkVeloDiffDirN / CellCenterDiffDir * NVec(1) &
                + DynamicVisc * BulkVeloDiffDirT1 / CellCenterDiffDir * T1Vec(1) &
                + DynamicVisc * BulkVeloDiffDirT2 / CellCenterDiffDir * T2Vec(1) )
      DeltaM(2) = DeltaM(2) + dt * Area * ( &
                  4. / 3. * DynamicVisc * BulkVeloDiffDirN / CellCenterDiffDir * NVec(2) &
                + DynamicVisc * BulkVeloDiffDirT1 / CellCenterDiffDir * T1Vec(2) &
                + DynamicVisc * BulkVeloDiffDirT2 / CellCenterDiffDir * T2Vec(2) )
      DeltaM(3) = DeltaM(3) + dt * Area * ( &
                  4. / 3. * DynamicVisc * BulkVeloDiffDirN / CellCenterDiffDir * NVec(3) &
                + DynamicVisc * BulkVeloDiffDirT1 / CellCenterDiffDir * T1Vec(3) &
                + DynamicVisc * BulkVeloDiffDirT2 / CellCenterDiffDir * T2Vec(3) )
      DeltaE = DeltaE + dt * Area * ( &
               ThermalCond * BulkTempDiff / CellCenterDiffDir &
             + 4. / 3. * DynamicVisc * BulkVeloDiffDirN / CellCenterDiffDir * MeanBulkVeloDirN &
                + DynamicVisc * BulkVeloDiffDirT1 / CellCenterDiffDir * MeanBulkVeloDirT1 &
                + DynamicVisc * BulkVeloDiffDirT2 / CellCenterDiffDir * MeanBulkVeloDirT2 )
    END IF
  ELSE IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(BC(SideID))).EQ.PartBound%ReflectiveBC) THEN
    MomentumACC = PartBound%MomentumACC(PartBound%MapToPartBC(BC(SideID)))
    TransACC = PartBound%TransACC(PartBound%MapToPartBC(BC(SideID)))
    TempWall = PartBound%WallTemp(PartBound%MapToPartBC(BC(SideID)))
    TempCell = BulkValues(iElem)%BulkTemperature
    WallVelo(1) = PartBound%WallVelo(1,PartBound%MapToPartBC(BC(SideID)))
    WallVelo(2) = PartBound%WallVelo(2,PartBound%MapToPartBC(BC(SideID)))
    WallVelo(3) = PartBound%WallVelo(3,PartBound%MapToPartBC(BC(SideID)))
    SpecR = BulkValues(iElem)%SpezGasConst
    BulkVeloCellT1 = (BulkValues(iElem)%CellV(1)-WallVelo(1)) * T1Vec(1) &
                   + (BulkValues(iElem)%CellV(2)-WallVelo(2)) * T1Vec(2) &
                   + (BulkValues(iElem)%CellV(3)-WallVelo(3)) * T1Vec(3)
    BulkVeloCellT2 = (BulkValues(iElem)%CellV(1)-WallVelo(1)) * T2Vec(1) &
                   + (BulkValues(iElem)%CellV(2)-WallVelo(2)) * T2Vec(2) &
                   + (BulkValues(iElem)%CellV(3)-WallVelo(3)) * T2Vec(3)

    DeltaM(1) = DeltaM(1) - MomentumACC * dt * Area * (BulkVeloCellT1 * T1Vec(1) + BulkVeloCellT2 * T2Vec(1)) &
              * BulkValues(iElem)%MassDens * SQRT(SpecR * TempCell /(2 * PI))
    DeltaM(2) = DeltaM(2) - MomentumACC * dt * Area * (BulkVeloCellT1 * T1Vec(2) + BulkVeloCellT2 * T2Vec(2)) &
              * BulkValues(iElem)%MassDens * SQRT(SpecR * TempCell /(2 * PI))
    DeltaM(3) = DeltaM(3) - MomentumACC * dt * Area * (BulkVeloCellT1 * T1Vec(3) + BulkVeloCellT2 * T2Vec(3)) &
              * BulkValues(iElem)%MassDens * SQRT(SpecR * TempCell /(2 * PI))
    DeltaE = DeltaE - TransACC * 0.5 * dt * Area * (BulkValues(iElem)%DegreeOfFreedom + 1.0) &
           * BulkValues(iElem)%MassDens * SpecR * (TempCell - TempWall) * SQRT(SpecR * TempCell /(2 * PI))
  ELSE IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(BC(SideID))).EQ.PartBound%SymmetryBC) THEN
CALL abort(&
__STAMP__&
,'SymmetryBC is not implemented for LD!')
  END IF
END SUBROUTINE CalcViscousTerms

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
END MODULE MOD_LD_reassign_part_prop
