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

#ifdef DONONTCOMPILETHIS
MODULE MOD_LD_DSMC_TOOLS
!==================================================================================================
! module contains tools for coupled DSMC LD method
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE LD_DSMC_Indicate_DSMC_Particles
  MODULE PROCEDURE LD_DSMC_Indicate_DSMC_Particles
END INTERFACE
INTERFACE LD_DSMC_Clone_Particles
  MODULE PROCEDURE LD_DSMC_Clone_Particles
END INTERFACE
INTERFACE LD_DSMC_Clean_Bufferregion
  MODULE PROCEDURE LD_DSMC_Clean_Bufferregion
END INTERFACE
INTERFACE LD_DSMC_data_sampling
  MODULE PROCEDURE LD_DSMC_data_sampling
END INTERFACE
INTERFACE LD_DSMC_output_calc
  MODULE PROCEDURE LD_DSMC_output_calc
END INTERFACE
INTERFACE LD_DSMC_Mean_Bufferzone_A_Val
  MODULE PROCEDURE LD_DSMC_Mean_Bufferzone_A_Val
END INTERFACE

PUBLIC :: LD_DSMC_Indicate_DSMC_Particles, LD_DSMC_Clone_Particles, &
           LD_DSMC_Clean_Bufferregion, LD_DSMC_data_sampling, LD_DSMC_output_calc, LD_DSMC_Mean_Bufferzone_A_Val
!===================================================================================================================================
CONTAINS
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE LD_DSMC_Indicate_DSMC_Particles
!===================================================================================================================================
! Subroutine to reset PartBulkValues as indicator for DSMC-particle in DSMC and Bufferzone_A cells
!===================================================================================================================================
! MODULES
  USE MOD_Mesh_Vars,             ONLY : nElems
  USE MOD_LD_Vars,               ONLY : PartStateBulkValues, BulkValues
  USE MOD_Particle_Vars,         ONLY : PEM
!  USE MOD_TimeDisc_Vars,         ONLY : iter
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

  INTEGER           :: iElem, iPart, iPartIndx, nPart
!===================================================================================================================================

DO iElem = 1, nElems
  IF((BulkValues(iElem)%CellType.EQ.1).OR.(BulkValues(iElem)%CellType.EQ.2)) THEN
    nPart = PEM%pNumber(iElem)
    iPartIndx = PEM%pStart(iElem)
!    IF((iter.EQ.0).AND.(BulkValues(iElem)%CellType.EQ.2)) THEN ! ACHTUNG nur 1 Spezes!!!!!!!!!!!!!!!!!
!        BulkValues(iElem)%CellV(1) = PartStateBulkValues(iPartIndx,1) ! ACHTUNG nur 1 Spezes!!!!!!!!!!!!!!!!!
!        BulkValues(iElem)%CellV(2) = PartStateBulkValues(iPartIndx,2) ! ACHTUNG nur 1 Spezes!!!!!!!!!!!!!!!!!
!        BulkValues(iElem)%CellV(3) = PartStateBulkValues(iPartIndx,3) ! ACHTUNG nur 1 Spezes!!!!!!!!!!!!!!!!!
!        BulkValues(iElem)%BulkTemperature = PartStateBulkValues(iPartIndx,4) ! ACHTUNG nur 1 Spezes!!!!!!!!!!!!!!!!!
!        BulkValues(iElem)%DegreeOfFreedom = PartStateBulkValues(iPartIndx,5) ! ACHTUNG nur 1 Spezes!!!!!!!!!!!!!!!!!
!    END IF

    DO iPart = 1, nPart    ! Reset PartBulkValues as Indicator for DSMC-Particle
      PartStateBulkValues(iPartIndx,1) = 0.0
      PartStateBulkValues(iPartIndx,2) = 0.0
      PartStateBulkValues(iPartIndx,3) = 0.0
      PartStateBulkValues(iPartIndx,4) = 0.0
      iPartIndx = PEM%pNext(iPartIndx)
    END DO
  END IF
END DO

END SUBROUTINE LD_DSMC_Indicate_DSMC_Particles

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE LD_DSMC_Clone_Particles
!===================================================================================================================================
! Subroutine to clone particles in Bufferregion
! Bufferzone_A:  clone all DSMC-Particles to LD-Particles
! Bufferzone_b:  clone all LD-Particles to DSMC-Particles
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_Mesh_Vars,             ONLY : nElems
  USE MOD_LD_Vars,               ONLY : PartStateBulkValues, BulkValues, LD_DSMC_RHS
  USE MOD_Particle_Vars,         ONLY : Species, PartSpecies, usevMPF, PartMPF, PEM, PartState, PDM
  USE MOD_LD_Init,               ONLY : CalcDegreeOfFreedom
  USE MOD_DSMC_Vars,             ONLY : SpecDSMC, PartStateIntEn, CollisMode
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

  INTEGER           :: iElem, iPart, iPartIndx, iCloneIndx, SumOfFormedClones, nPart
  REAL               :: RandVec(2)
  REAL               :: LocalPI=3.14159265359
  REAL               :: CellMeanPartMass, WeightFak, CellMass, MPFSum, iRan

!===================================================================================================================================

SumOfFormedClones = 0
DO iElem = 1, nElems
                ! --------
                ! --------
  IF(BulkValues(iElem)%CellType.EQ.2) THEN
                ! -------- => Bufferzone_A:  clone all DSMC-Particles to LD-Particles
    nPart = PEM%pNumber(iElem)
    iPartIndx = PEM%pStart(iElem)
    DO iPart = 1, nPart
                ! --------get free PartIndex for Clone--------
      SumOfFormedClones = SumOfFormedClones + 1
      iCloneIndx = PDM%nextFreePosition(SumOfFormedClones+PDM%CurrentNextFreePosition)
      IF (iCloneIndx.EQ.0) THEN
        CALL abort(&
        __STAMP__&
        ,'ERROR: New Particle (Clone Region-A) Number greater max Part Num')
      END IF
                ! --------create Clone--------
      PEM%Element(iCloneIndx) = PEM%Element(iPartIndx)
      PDM%ParticleInside(iCloneIndx) = .TRUE.
      PartSpecies(iCloneIndx) = PartSpecies(iPartIndx)
      LD_DSMC_RHS(1,iCloneIndx)=0.0
      LD_DSMC_RHS(2,iCloneIndx)=0.0
      LD_DSMC_RHS(3,iCloneIndx)=0.0
      IF (usevMPF) THEN
        PartMPF(iCloneIndx) = PartMPF(iPartIndx)
      END IF
      PartState(1,iCloneIndx) = PartState(1,iPartIndx)
      PartState(2,iCloneIndx) = PartState(2,iPartIndx)
      PartState(3,iCloneIndx) = PartState(3,iPartIndx)
                ! --------calculate LD-Values for Clone--------
      PartState(4,iCloneIndx) = BulkValues(iElem)%CellV(1)
      PartState(5,iCloneIndx) = BulkValues(iElem)%CellV(2)
      PartState(6,iCloneIndx) = BulkValues(iElem)%CellV(3)
      PartStateBulkValues(iCloneIndx,1) = BulkValues(iElem)%CellV(1)
      PartStateBulkValues(iCloneIndx,2) = BulkValues(iElem)%CellV(2)
      PartStateBulkValues(iCloneIndx,3) = BulkValues(iElem)%CellV(3)
      PartStateBulkValues(iCloneIndx,4) = BulkValues(iElem)%BulkTemperature
      PartStateBulkValues(iCloneIndx,5) = CalcDegreeOfFreedom(iCloneIndx)
      iPartIndx = PEM%pNext(iPartIndx)
    END DO
                ! --------end clone all DSMC-Particles to LD-Particles
                ! --------
                ! --------
  ELSE IF(BulkValues(iElem)%CellType.EQ.3) THEN
                ! -------- => Bufferzone_b:  clone all LD-Particles to DSMC-Particles
                ! --------first calculate mean part mass--------
    CellMass  = 0.0
    MPFSum    = 0.0
    nPart = PEM%pNumber(iElem)
    iPartIndx = PEM%pStart(iElem)
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
    CellMeanPartMass = CellMass / MPFSum
                ! --------end calculate mean part mass--------
    nPart = PEM%pNumber(iElem)
    iPartIndx = PEM%pStart(iElem)
    DO iPart = 1, nPart
                ! --------get free PartIndex for Clone--------
      SumOfFormedClones = SumOfFormedClones + 1
      iCloneIndx = PDM%nextFreePosition(SumOfFormedClones+PDM%CurrentNextFreePosition)
      IF (iCloneIndx.EQ.0) THEN
        CALL abort(&
        __STAMP__&
        ,'ERROR: New Particle (Clone Region-B) Number greater max Part Num')
      END IF
                ! --------create Clone--------
      PEM%Element(iCloneIndx) = PEM%Element(iPartIndx)
      PDM%ParticleInside(iCloneIndx) = .TRUE.
      PartSpecies(iCloneIndx) = PartSpecies(iPartIndx)
      LD_DSMC_RHS(1,iCloneIndx)=0.0
      LD_DSMC_RHS(2,iCloneIndx)=0.0
      LD_DSMC_RHS(3,iCloneIndx)=0.0
      IF (usevMPF) THEN
        PartMPF(iCloneIndx) = PartMPF(iPartIndx)
      END IF
      PartState(1,iCloneIndx) = PartState(1,iPartIndx)
      PartState(2,iCloneIndx) = PartState(2,iPartIndx)
      PartState(3,iCloneIndx) = PartState(3,iPartIndx)
                ! --------calculate DSMC-Values for Clone--------
      CALL RANDOM_NUMBER(RandVec)
      PartState(4,iCloneIndx) = BulkValues(iElem)%CellV(1) &
                              + SIN(2*LocalPI * RandVec(1)) / BulkValues(iElem)%Beta &
                              * SQRT((-1)*CellMeanPartMass / Species(PartSpecies(iCloneIndx))%MassIC * LOG(RandVec(2)))
      CALL RANDOM_NUMBER(RandVec)
      PartState(5,iCloneIndx) = BulkValues(iElem)%CellV(2) &
                              + SIN(2*LocalPI * RandVec(1)) / BulkValues(iElem)%Beta &
                              * SQRT((-1)*CellMeanPartMass / Species(PartSpecies(iCloneIndx))%MassIC * LOG(RandVec(2)))
      CALL RANDOM_NUMBER(RandVec)
      PartState(6,iCloneIndx) = BulkValues(iElem)%CellV(3) &
                              + SIN(2*LocalPI * RandVec(1)) / BulkValues(iElem)%Beta &
                              * SQRT((-1)*CellMeanPartMass / Species(PartSpecies(iCloneIndx))%MassIC * LOG(RandVec(2)))

      IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.&
          (SpecDSMC(PartSpecies(iCloneIndx))%InterID.EQ.2)) THEN
        PartStateIntEn(iCloneIndx,1) = 0.0
        CALL RANDOM_NUMBER(iRan)
        PartStateIntEn(iCloneIndx,2) = - CellMeanPartMass / (2 * BulkValues(iElem)%Beta**2) * LOG(iRan)
      END IF

      PartStateBulkValues(iCloneIndx,1) = 0.0
      PartStateBulkValues(iCloneIndx,2) = 0.0
      PartStateBulkValues(iCloneIndx,3) = 0.0
      PartStateBulkValues(iCloneIndx,4) = PartStateBulkValues(iPartIndx,4)
      PartStateBulkValues(iCloneIndx,5) = CalcDegreeOfFreedom(iCloneIndx)
      PartStateBulkValues(iCloneIndx,4) = 0.0
      iPartIndx = PEM%pNext(iPartIndx)
    END DO
                ! --------end clone all LD-Particles to DSMC-Particles
  END IF
                ! --------
                ! --------
END DO
                ! --------reassign global PartVec and FreePos
PDM%ParticleVecLength = PDM%ParticleVecLength + SumOfFormedClones
PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + SumOfFormedClones

END SUBROUTINE LD_DSMC_Clone_Particles

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE LD_DSMC_Clean_Bufferregion
!===================================================================================================================================
! Subroutine to remove all alien particles in LD or DSMC cells
!===================================================================================================================================
! MODULES
  USE MOD_LD_Vars,               ONLY : PartStateBulkValues, BulkValues
  USE MOD_Particle_Vars,         ONLY : PEM, PDM
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER           :: iPartIndx, iElem
!===================================================================================================================================
DO iPartIndx = 1, PDM%ParticleVecLength
  IF (.NOT. PDM%ParticleInside(iPartIndx)) CYCLE
  iElem = PEM%Element(iPartIndx)
  IF((BulkValues(iElem)%CellType.EQ.1).OR.(BulkValues(iElem)%CellType.EQ.2)) THEN
    IF(PartStateBulkValues(iPartIndx,4).NE. 0.0) THEN ! Delete LD Part in DSMC or Buffer-A Cells
      PDM%ParticleInside(iPartIndx) = .FALSE.
    END IF
  ELSE
    IF(PartStateBulkValues(iPartIndx,4).EQ. 0.0) THEN ! Delete DSMC Part in LD or Buffer-B Cells
      PDM%ParticleInside(iPartIndx) = .FALSE.
    END IF
  END IF
END DO

END SUBROUTINE LD_DSMC_Clean_Bufferregion

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE LD_DSMC_data_sampling()
!===================================================================================================================================
! Sample celldata for output
!===================================================================================================================================
! MODULES
  USE MOD_LD_Vars
  USE MOD_DSMC_Vars,              ONLY : DSMC, SampDSMC, PartStateIntEn, CollisMode, SpecDSMC
  USE MOD_Particle_Vars,          ONLY : PEM, PDM, PartSpecies, PartState, PartMPF, usevMPF
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                 !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
INTEGER                       :: iPart, iElem
!--------------------------------------------------------------------------------------------------!

DSMC%SampNum = DSMC%SampNum + 1
DO iPart = 1, PDM%ParticleVecLength
  IF (PDM%ParticleInside(ipart)) THEN
    iElem = PEM%Element(iPart)
    IF ((BulkValues(iElem)%CellType.EQ.3).OR.(BulkValues(iElem)%CellType.EQ.4)) THEN ! -----------------------------> LD
      IF (usevMPF) THEN
        SampDSMC(iElem,PartSpecies(iPart))%PartV(1:3)  = SampDSMC(iElem,PartSpecies(iPart))%PartV(1:3) &
                                                       + PartStateBulkValues(iPart,1:3) * PartMPF(iPart)
        SampDSMC(iElem,PartSpecies(iPart))%PartV2(1)   = SampDSMC(iElem,PartSpecies(iPart))%PartV2(1) &  ! V2 stands for BulkTemp
                                                       + PartStateBulkValues(iPart,4) * PartMPF(iPart)
        SampDSMC(iElem,PartSpecies(iPart))%PartNum     = SampDSMC(iElem,PartSpecies(iPart))%PartNum + PartMPF(iPart)
        SampDSMC(iElem,PartSpecies(iPart))%SimPartNum  = SampDSMC(iElem,PartSpecies(iPart))%SimPartNum + 1
      ELSE ! normal sampling without weighting
        SampDSMC(iElem,PartSpecies(iPart))%PartV(1:3)  = SampDSMC(iElem,PartSpecies(iPart))%PartV(1:3) &
                                                       + PartStateBulkValues(iPart,1:3)
        SampDSMC(iElem,PartSpecies(iPart))%PartV2(1)   = SampDSMC(iElem,PartSpecies(iPart))%PartV2(1) & !!! V2 stands for BulkTemp
                                                       + PartStateBulkValues(iPart,4)
        SampDSMC(iElem,PartSpecies(iPart))%PartNum     = SampDSMC(iElem,PartSpecies(iPart))%PartNum + 1
      END IF
    ELSE ! -----------------------------> DSMC
      IF (usevMPF) THEN
        SampDSMC(iElem,PartSpecies(iPart))%PartV(1:3)  = SampDSMC(iElem,PartSpecies(iPart))%PartV(1:3) &
                                                       + PartState(4:6,iPart) * PartMPF(iPart)
        SampDSMC(iElem,PartSpecies(iPart))%PartV2(1:3) = SampDSMC(iElem,PartSpecies(iPart))%PartV2(1:3) &
                                                       + PartState(4:6,iPart)**2 * PartMPF(iPart)
        SampDSMC(iElem,PartSpecies(iPart))%PartNum     = SampDSMC(iElem,PartSpecies(iPart))%PartNum + PartMPF(iPart)
        SampDSMC(iElem,PartSpecies(iPart))%SimPartNum  = SampDSMC(iElem,PartSpecies(iPart))%SimPartNum + 1
        ! if usevMPF SampDSMC(iElem,PartSpecies(iPart))%PartNum == real number of particles
        IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.&
            (SpecDSMC(PartSpecies(iPart))%InterID.EQ.2)) THEN
          SampDSMC(iElem,PartSpecies(iPart))%EVib      = SampDSMC(iElem,PartSpecies(iPart))%EVib &
                                                       + PartStateIntEn(iPart,1) * PartMPF(iPart)
          SampDSMC(iElem,PartSpecies(iPart))%ERot      = SampDSMC(iElem,PartSpecies(iPart))%ERot &
                                                       + PartStateIntEn(iPart,2) * PartMPF(iPart)
        END IF
        IF (DSMC%ElectronicModel) THEN
          SampDSMC(iElem,PartSpecies(iPart))%EElec     = SampDSMC(iElem,PartSpecies(iPart))%EElec &
                                                       + PartStateIntEn(iPart,3) * PartMPF(iPart)
        END IF
      ELSE ! normal sampling without weighting
        SampDSMC(iElem,PartSpecies(iPart))%PartV(1:3)  = SampDSMC(iElem,PartSpecies(iPart))%PartV(1:3) &
                                                       + PartState(4:6,iPart)
        SampDSMC(iElem,PartSpecies(iPart))%PartV2(1:3) = SampDSMC(iElem,PartSpecies(iPart))%PartV2(1:3) &
                                                       + PartState(4:6,iPart)**2
        SampDSMC(iElem,PartSpecies(iPart))%PartNum     = SampDSMC(iElem,PartSpecies(iPart))%PartNum + 1
        IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.&
            (SpecDSMC(PartSpecies(iPart))%InterID.EQ.2)) THEN
          SampDSMC(iElem,PartSpecies(iPart))%EVib      = SampDSMC(iElem,PartSpecies(iPart))%EVib &
                                                       + PartStateIntEn(iPart,1)
          SampDSMC(iElem,PartSpecies(iPart))%ERot      = SampDSMC(iElem,PartSpecies(iPart))%ERot &
                                                       + PartStateIntEn(iPart,2)
        END IF
        IF (DSMC%ElectronicModel) THEN
          SampDSMC(iElem,PartSpecies(iPart))%EElec     = SampDSMC(iElem,PartSpecies(iPart))%EElec &
                                                       + PartStateIntEn(iPart,3)
        END IF
      END IF
    END IF
  END IF
END DO

END SUBROUTINE LD_DSMC_data_sampling

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE LD_DSMC_output_calc()
!===================================================================================================================================
! Calculation of outputdata on the basis of sampled values
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : DSMC, SampDSMC, MacroDSMC, CollisMode, SpecDSMC, realtime
  USE MOD_Mesh_Vars,              ONLY : nElems,MeshFile
  USE MOD_Particle_Vars,          ONLY : nSpecies, BoltzmannConst, Species, GEO, usevMPF, Time
  USE MOD_TimeDisc_Vars,          ONLY : TEnd, iter, dt
  USE MOD_Restart_Vars,           ONLY : RestartTime
  USE MOD_LD_Vars,                ONLY : BulkValues
  USE MOD_DSMC_Analyze,           ONLY : CalcTVib, CalcTelec, WriteDSMCToHDF5
!--------------------------------------------------------------------------------------------------!
! statistical pairing method                                                                       !
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                 !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
INTEGER                       :: iElem, iSpec                                              !
REAL                          :: TVib_TempFac, MolecPartNum, HeavyPartNum                                                 !
! input variable declaration                                                      !
!--------------------------------------------------------------------------------------------------!

ALLOCATE(MacroDSMC(nElems,nSpecies + 1))
MacroDSMC(1:nElems,1:nSpecies+1)%PartV(1)  = 0.0
MacroDSMC(1:nElems,1:nSpecies+1)%PartV(2)  = 0.0
MacroDSMC(1:nElems,1:nSpecies+1)%PartV(3)  = 0.0
MacroDSMC(1:nElems,1:nSpecies+1)%PartV(4)  = 0.0
MacroDSMC(1:nElems,1:nSpecies+1)%PartV2(1) = 0.0
MacroDSMC(1:nElems,1:nSpecies+1)%PartV2(2) = 0.0
MacroDSMC(1:nElems,1:nSpecies+1)%PartV2(3) = 0.0
MacroDSMC(1:nElems,1:nSpecies+1)%Temp(1)   = 0.0
MacroDSMC(1:nElems,1:nSpecies+1)%Temp(2)   = 0.0
MacroDSMC(1:nElems,1:nSpecies+1)%Temp(3)   = 0.0
MacroDSMC(1:nElems,1:nSpecies+1)%Temp(4)   = 0.0
MacroDSMC(1:nElems,1:nSpecies+1)%PartNum   = 0.0
MacroDSMC(1:nElems,1:nSpecies+1)%NumDens   = 0.0
MacroDSMC(1:nElems,1:nSpecies+1)%TVib      = 0.0
MacroDSMC(1:nElems,1:nSpecies+1)%TRot      = 0.0
MacroDSMC(1:nElems,1:nSpecies+1)%TElec     = 0.0

DO iSpec = 1, nSpecies
  DO iElem = 1, nElems ! element/cell main loop
    IF(SampDSMC(iElem,iSpec)%PartNum.GT. 0) THEN
! compute flow velocity
      MacroDSMC(iElem,iSpec)%PartV(1) = SampDSMC(iElem,iSpec)%PartV(1) / SampDSMC(iElem,iSpec)%PartNum
      MacroDSMC(iElem,iSpec)%PartV(2) = SampDSMC(iElem,iSpec)%PartV(2) / SampDSMC(iElem,iSpec)%PartNum
      MacroDSMC(iElem,iSpec)%PartV(3) = SampDSMC(iElem,iSpec)%PartV(3) / SampDSMC(iElem,iSpec)%PartNum
      MacroDSMC(iElem,iSpec)%PartV(4) = SQRT( MacroDSMC(iElem,iSpec)%PartV(1)**2 &
                                            + MacroDSMC(iElem,iSpec)%PartV(2)**2 &
                                            + MacroDSMC(iElem,iSpec)%PartV(3)**2 )
! compute flow Temperature
      IF((BulkValues(iElem)%CellType.EQ.1).OR.(BulkValues(iElem)%CellType.EQ.2)) THEN ! -----------------------------> DSMC
        MacroDSMC(iElem,iSpec)%PartV2(1) = SampDSMC(iElem,iSpec)%PartV2(1) / SampDSMC(iElem,iSpec)%PartNum
        MacroDSMC(iElem,iSpec)%PartV2(2) = SampDSMC(iElem,iSpec)%PartV2(2) / SampDSMC(iElem,iSpec)%PartNum
        MacroDSMC(iElem,iSpec)%PartV2(3) = SampDSMC(iElem,iSpec)%PartV2(3) / SampDSMC(iElem,iSpec)%PartNum
        MacroDSMC(iElem,iSpec)%Temp(1:3) = Species(iSpec)%MassIC / BoltzmannConst &
                                             * (MacroDSMC(iElem,iSpec)%PartV2(1:3) - MacroDSMC(iElem,iSpec)%PartV(1:3)**2)
        MacroDSMC(iElem,iSpec)%Temp(4) = (MacroDSMC(iElem,iSpec)%Temp(1) &
                                        + MacroDSMC(iElem,iSpec)%Temp(2) &
                                        + MacroDSMC(iElem,iSpec)%Temp(3)) / 3
      ELSE ! -----------------------------> LD (no direction depended temperature)
        MacroDSMC(iElem,iSpec)%Temp(1)   = SampDSMC(iElem,iSpec)%PartV2(1) / SampDSMC(iElem,iSpec)%PartNum
        MacroDSMC(iElem,iSpec)%Temp(2)   = SampDSMC(iElem,iSpec)%PartV2(1) / SampDSMC(iElem,iSpec)%PartNum
        MacroDSMC(iElem,iSpec)%Temp(3)   = SampDSMC(iElem,iSpec)%PartV2(1) / SampDSMC(iElem,iSpec)%PartNum
        MacroDSMC(iElem,iSpec)%Temp(4)   = SampDSMC(iElem,iSpec)%PartV2(1) / SampDSMC(iElem,iSpec)%PartNum
      END IF
! compute density
      MacroDSMC(iElem,iSpec)%PartNum = SampDSMC(iElem,iSpec)%PartNum / REAL(DSMC%SampNum)
      ! comment: if usevMPF MacroDSMC(iElem,iSpec)%PartNum == real number of particles
      IF (usevMPF) THEN
        MacroDSMC(iElem,iSpec)%NumDens = MacroDSMC(iElem,iSpec)%PartNum / GEO%Volume(iElem)
      ELSE
        MacroDSMC(iElem,iSpec)%NumDens = MacroDSMC(iElem,iSpec)%PartNum * Species(iSpec)%MacroParticleFactor / GEO%Volume(iElem)
      END IF
! compute internal energies / has to be changed for vfd
      IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.&
        (SpecDSMC(iSpec)%InterID.EQ.2)) THEN
        IF((BulkValues(iElem)%CellType.EQ.1).OR.(BulkValues(iElem)%CellType.EQ.2)) THEN ! -----------------------------> DSMC
          IF (DSMC%VibEnergyModel.EQ.0) THEN              ! SHO-model
            TVib_TempFac=SampDSMC(iElem,iSpec)%EVib &
                      /(SampDSMC(iElem,iSpec)%PartNum*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)
            IF (TVib_TempFac.LE.DSMC%GammaQuant) THEN
              MacroDSMC(iElem,iSpec)%TVib = 0.0
            ELSE
              MacroDSMC(iElem,iSpec)%TVib = SpecDSMC(iSpec)%CharaTVib/LOG(1 + 1/(TVib_TempFac-DSMC%GammaQuant))
            END IF
          ELSE                                            ! TSHO-model
            MacroDSMC(iElem,iSpec)%TVib = CalcTVib(SpecDSMC(iSpec)%CharaTVib &
                , SampDSMC(iElem,iSpec)%EVib/SampDSMC(iElem,iSpec)%PartNum, SpecDSMC(iSpec)%MaxVibQuant)
          END IF
          MacroDSMC(iElem,iSpec)%TRot = SampDSMC(iElem, iSpec)%ERot/(BoltzmannConst*SampDSMC(iElem,iSpec)%PartNum)
          IF (DSMC%ElectronicModel) THEN
            MacroDSMC(iElem,iSpec)%TElec= CalcTelec( SampDSMC(iElem,iSpec)%EElec/SampDSMC(iElem,iSpec)%PartNum, iSpec)
          END IF
        ELSE ! -----------------------------> LD (no differentiation in LD)
          MacroDSMC(iElem,iSpec)%TVib = MacroDSMC(iElem,iSpec)%Temp(4)
          MacroDSMC(iElem,iSpec)%TRot = MacroDSMC(iElem,iSpec)%Temp(4)
          IF (DSMC%ElectronicModel) THEN
            MacroDSMC(iElem,iSpec)%TElec= MacroDSMC(iElem,iSpec)%Temp(4)
          END IF
        END IF
      END IF
    END IF
  END DO
END DO
! compute total values
DO iElem = 1, nElems ! element/cell main loop
  MolecPartNum = 0
  HeavyPartNum = 0
  DO iSpec = 1, nSpecies
    IF(SampDSMC(iElem,iSpec)%PartNum.GT. 0) THEN
      MacroDSMC(iElem,nSpecies + 1)%PartV(1) = MacroDSMC(iElem,nSpecies + 1)%PartV(1) &
                                             + MacroDSMC(iElem,iSpec)%PartV(1) * MacroDSMC(iElem,iSpec)%PartNum
      MacroDSMC(iElem,nSpecies + 1)%PartV(2) = MacroDSMC(iElem,nSpecies + 1)%PartV(2) &
                                             + MacroDSMC(iElem,iSpec)%PartV(2) * MacroDSMC(iElem,iSpec)%PartNum
      MacroDSMC(iElem,nSpecies + 1)%PartV(3) = MacroDSMC(iElem,nSpecies + 1)%PartV(3) &
                                             + MacroDSMC(iElem,iSpec)%PartV(3) * MacroDSMC(iElem,iSpec)%PartNum
      MacroDSMC(iElem,nSpecies + 1)%Temp(1)  = MacroDSMC(iElem,nSpecies + 1)%Temp(1) &
                                             + MacroDSMC(iElem,iSpec)%Temp(1) * MacroDSMC(iElem,iSpec)%PartNum
      MacroDSMC(iElem,nSpecies + 1)%Temp(2)  = MacroDSMC(iElem,nSpecies + 1)%Temp(2) &
                                             + MacroDSMC(iElem,iSpec)%Temp(2) * MacroDSMC(iElem,iSpec)%PartNum
      MacroDSMC(iElem,nSpecies + 1)%Temp(3)  = MacroDSMC(iElem,nSpecies + 1)%Temp(3) &
                                             + MacroDSMC(iElem,iSpec)%Temp(3) * MacroDSMC(iElem,iSpec)%PartNum
      MacroDSMC(iElem,nSpecies + 1)%PartNum  = MacroDSMC(iElem,nSpecies + 1)%PartNum + MacroDSMC(iElem,iSpec)%PartNum
      IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.&
        (SpecDSMC(iSpec)%InterID.EQ.2)) THEN
            MacroDSMC(iElem,nSpecies + 1)%TVib = MacroDSMC(iElem,nSpecies + 1)%TVib &
                                               + MacroDSMC(iElem,iSpec)%TVib * MacroDSMC(iElem,iSpec)%PartNum
            MacroDSMC(iElem,nSpecies + 1)%TRot = MacroDSMC(iElem,nSpecies + 1)%TRot &
                                               + MacroDSMC(iElem,iSpec)%TRot * MacroDSMC(iElem,iSpec)%PartNum
            MolecPartNum                       = MolecPartNum + MacroDSMC(iElem,iSpec)%PartNum
      END IF
      IF ( DSMC%ElectronicModel .AND. (SpecDSMC(iSpec)%InterID.NE.4) ) THEN
        MacroDSMC(iElem,nSpecies + 1)%TElec = MacroDSMC(iElem, nSpecies+1)%TElec &
                                            + MacroDSMC(iElem,iSpec)%TElec * MacroDSMC(iElem,iSpec)%PartNum
        HeavyPartNum                        = HeavyPartNum + MacroDSMC(iElem,iSpec)%PartNum
      END IF
    END IF
  END DO
  IF(MacroDSMC(iElem,nSpecies + 1)%PartNum.GT. 0) THEN
    MacroDSMC(iElem,nSpecies + 1)%PartV(1) = MacroDSMC(iElem,nSpecies + 1)%PartV(1) / MacroDSMC(iElem,nSpecies + 1)%PartNum
    MacroDSMC(iElem,nSpecies + 1)%PartV(2) = MacroDSMC(iElem,nSpecies + 1)%PartV(2) / MacroDSMC(iElem,nSpecies + 1)%PartNum
    MacroDSMC(iElem,nSpecies + 1)%PartV(3) = MacroDSMC(iElem,nSpecies + 1)%PartV(3) / MacroDSMC(iElem,nSpecies + 1)%PartNum
    MacroDSMC(iElem,nSpecies + 1)%Temp(1)  = MacroDSMC(iElem,nSpecies + 1)%Temp(1)  / MacroDSMC(iElem,nSpecies + 1)%PartNum
    MacroDSMC(iElem,nSpecies + 1)%Temp(2)  = MacroDSMC(iElem,nSpecies + 1)%Temp(2)  / MacroDSMC(iElem,nSpecies + 1)%PartNum
    MacroDSMC(iElem,nSpecies + 1)%Temp(3)  = MacroDSMC(iElem,nSpecies + 1)%Temp(3)  / MacroDSMC(iElem,nSpecies + 1)%PartNum
    IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.&
          (MolecPartNum.GT. 0)) THEN
              MacroDSMC(iElem,nSpecies + 1)%TVib = MacroDSMC(iElem,nSpecies + 1)%TVib / MolecPartNum
              MacroDSMC(iElem,nSpecies + 1)%TRot = MacroDSMC(iElem,nSpecies + 1)%TRot / MolecPartNum


!!              IF(MolecPartNum .ne. MacroDSMC(iElem,nSpecies)%PartNum) THEN
!!                print*,'MolecPartNum .ne. MacroDSMC(iElem,nSpecies + 1)%PartNum'
!!                stop
!!              END IF
!!              IF(MolecPartNum .ne. MacroDSMC(iElem,nSpecies + 1)%PartNum) THEN
!!                print*,'MolecPartNum .ne. MacroDSMC(iElem,nSpecies + 1)%PartNum'
!!                stop
!!              END IF
!!              IF(MacroDSMC(iElem,nSpecies + 1)%TRot .ne. MacroDSMC(iElem,1)%TRot) THEN

!!                print*,'MolecPartNum .ne. MacroDSMC(iElem,nSpecies + 1)%PartNum'
!!                print*,MacroDSMC(iElem,nSpecies + 1)%TRot,MacroDSMC(iElem,1)%TRot,MacroDSMC(iElem,1)%TRot &
!!                                                                                 - MacroDSMC(iElem,nSpecies + 1)%TRot
!!                print*,MolecPartNum, MacroDSMC(iElem,nSpecies)%PartNum, MolecPartNum - MacroDSMC(iElem,nSpecies)%PartNum
!!!                stop

!!              END IF

    END IF


    IF ( DSMC%ElectronicModel) THEN
      MacroDSMC(iElem,nSpecies + 1)%TElec = MacroDSMC(iElem, nSpecies+1)%TElec / HeavyPartNum
    END IF
  END IF
  MacroDSMC(iElem,nSpecies + 1)%PartV(4) = SQRT(MacroDSMC(iElem,nSpecies + 1)%PartV(1)**2 &
                                          + MacroDSMC(iElem,nSpecies + 1)%PartV(2)**2 &
                                          + MacroDSMC(iElem,nSpecies + 1)%PartV(3)**2)
  MacroDSMC(iElem,nSpecies + 1)%Temp(4)  = (MacroDSMC(iElem,nSpecies + 1)%Temp(1) &
                                          + MacroDSMC(iElem,nSpecies + 1)%Temp(2) &
                                          + MacroDSMC(iElem,nSpecies + 1)%Temp(3)) /3
  MacroDSMC(iElem,nSpecies + 1)%NumDens = MacroDSMC(iElem,nSpecies + 1)%PartNum * Species(1)%MacroParticleFactor &
                                 / GEO%Volume(iElem) ! the calculation is limitied for MPF = const.
  IF (usevMPF) THEN
    MacroDSMC(iElem,nSpecies + 1)%PartNum = 0
    DO iSpec = 1, nSpecies
      MacroDSMC(iElem,iSpec)%PartNum = SampDSMC(iElem,iSpec)%SimPartNum / REAL(DSMC%SampNum)
      MacroDSMC(iElem,nSpecies + 1)%PartNum = MacroDSMC(iElem,nSpecies + 1)%PartNum + MacroDSMC(iElem,iSpec)%PartNum
    END DO
  END IF

!-------------------------------------------DEBUG OUTPUT---------------------------------------------
   MacroDSMC(iElem,nSpecies + 1)%Temp(1)  = CBC_Par(iElem)%Celltype
   MacroDSMC(iElem,nSpecies + 1)%Temp(2)  = CBC_Par(iElem)%TVib
   MacroDSMC(iElem,nSpecies + 1)%Temp(3)  = CBC_Par(iElem)%ZetaRot
   MacroDSMC(iElem,nSpecies + 1)%Temp(4)  = CBC_Par(iElem)%ZetaVib
   MacroDSMC(iElem,nSpecies + 1)%PartV(1) = CBC_Par(iElem)%KnTemp
   MacroDSMC(iElem,nSpecies + 1)%PartV(2) = CBC_Par(iElem)%KnTGLL
   MacroDSMC(iElem,nSpecies + 1)%PartV(4) = REAL(BulkValues(iElem)%CellType)
!   MacroDSMC(iElem,nSpecies + 1)%NumDens  = CBC_Par(iElem)%KnDens
   MacroDSMC(iElem,nSpecies + 1)%PartNum  = REAL(BulkValues(iElem)%CellType)
   MacroDSMC(iElem,nSpecies + 1)%PartNum  = CBC_Par(iElem)%KnMax
!   MacroDSMC(iElem,nSpecies + 1)%PartNum  = GEO%NumNeighborElems(iElem) + MPIGEO%NumNeighborElems(iElem)

!-------------------------------------------DEBUG OUTPUT---------------------------------------------
  IF(DSMC%CalcQualityFactors) THEN
    IF((BulkValues(iElem)%CellType.EQ.1).OR.(BulkValues(iElem)%CellType.EQ.2)) THEN ! -----------------------------> DSMC
      IF ((MacroDSMC(iElem,nSpecies+1)%PartNum.GT.0).AND.(MacroDSMC(iElem,nSpecies + 1)%Temp(4).GT.0)) THEN
        DSMC%QualityFactors(iElem,3) = DSMC%QualityFacSamp(iElem,3) &
                              / CalcMeanFreePath(MacroDSMC(iElem,1:nSpecies)%PartNum, MacroDSMC(iElem,nSpecies+1)%PartNum, &
                                  GEO%Volume(iElem), SpecDSMC(1)%omegaVHS, MacroDSMC(iElem,nSpecies + 1)%Temp(4))
      END IF
      IF(WriteMacroValues) THEN
        DSMC%QualityFactors(iElem,1) = DSMC%QualityFacSamp(iElem,1) / iter_macvalout
        DSMC%QualityFactors(iElem,2) = DSMC%QualityFacSamp(iElem,2) / iter_macvalout
        DSMC%QualityFactors(iElem,3) = DSMC%QualityFactors(iElem,3) / iter_macvalout
      ELSE
        IF (RestartTime.GT.(1-DSMC%TimeFracSamp)*TEnd) THEN
          DSMC%QualityFactors(iElem,1) = DSMC%QualityFacSamp(iElem,1) / iter
          DSMC%QualityFactors(iElem,2) = DSMC%QualityFacSamp(iElem,2) / iter
          DSMC%QualityFactors(iElem,3) = DSMC%QualityFactors(iElem,3) / iter
        ELSE
          DSMC%QualityFactors(iElem,1) = DSMC%QualityFacSamp(iElem,1)*dt / (Time-(1-DSMC%TimeFracSamp)*TEnd)
          DSMC%QualityFactors(iElem,2) = DSMC%QualityFacSamp(iElem,2)*dt / (Time-(1-DSMC%TimeFracSamp)*TEnd)
          DSMC%QualityFactors(iElem,3) = DSMC%QualityFactors(iElem,3)*dt / (Time-(1-DSMC%TimeFracSamp)*TEnd)
        END IF
      END IF
    ELSE
      DSMC%QualityFactors(iElem,:) = 0.0
    END IF
  ELSE
    DSMC%CollProbOut(iElem,:) = 0.0
  END IF
END DO

CALL WriteDSMCToHDF5(TRIM(MeshFile),realtime)

DEALLOCATE(MacroDSMC)

END SUBROUTINE LD_DSMC_output_calc

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE LD_DSMC_Mean_Bufferzone_A_Val(iElem)
!===================================================================================================================================
! Calculation of bulkvalues for bufferregion A
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_LD_Vars,               ONLY : PartStateBulkValues, BulkValues, LD_DSMC_RelaxationFak_BufferA
  USE MOD_Particle_Vars,         ONLY : PEM, PartState, usevMPF, GEO, Species, PartSpecies, BoltzmannConst, usevMPF, PartMPF
  USE MOD_DSMC_Vars,             ONLY : SpecDSMC, PartStateIntEn, CollisMode
!  USE MOD_TimeDisc_Vars,         ONLY : iter
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER,INTENT(IN) :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

  INTEGER             :: iPart, nPart, iPartIndx
  REAL                :: MPFSum, WeightFak, MPFSumOfMolecs
  REAL                :: CellMass, CellTemp, CellVelo2, MeanInternalEnergie
  REAL                :: OldBulkValues(7)
!===================================================================================================================================

nPart = PEM%pNumber(iElem)
IF (nPart.GT. 1) THEN ! Are there more than one particle
  OldBulkValues(1:3) = BulkValues(iElem)%CellV(1:3)
  OldBulkValues(4)   = BulkValues(iElem)%DegreeOfFreedom
  OldBulkValues(5)   = BulkValues(iElem)%MassDens
  OldBulkValues(6)   = BulkValues(iElem)%BulkTemperature
  OldBulkValues(7)   = BulkValues(iElem)%Beta

  BulkValues(iElem)%CellV           = 0.0
  BulkValues(iElem)%MassDens        = 0.0
  BulkValues(iElem)%DegreeOfFreedom = 0.0
  CellMass                          = 0.0
  MPFSum                            = 0.0
  MPFSumOfMolecs                    = 0.0
  CellVelo2                         = 0.0
  MeanInternalEnergie               = 0.0
  iPartIndx = PEM%pStart(iElem)
  DO iPart = 1, nPart
    IF (usevMPF) THEN
       WeightFak = PartMPF(iPartIndx)
    ELSE
       WeightFak = Species(PartSpecies(iPartIndx))%MacroParticleFactor
    END IF
    BulkValues(iElem)%CellV(1)        = BulkValues(iElem)%CellV(1) + PartState(4,iPartIndx) &
                                      * WeightFak * Species(PartSpecies(iPartIndx))%MassIC
    BulkValues(iElem)%CellV(2)        = BulkValues(iElem)%CellV(2) + PartState(5,iPartIndx) &
                                      * WeightFak * Species(PartSpecies(iPartIndx))%MassIC
    BulkValues(iElem)%CellV(3)        = BulkValues(iElem)%CellV(3) + PartState(6,iPartIndx) &
                                      * WeightFak * Species(PartSpecies(iPartIndx))%MassIC
    CellVelo2                         = CellVelo2 + ( (PartState(4,iPartIndx))**2 &
                                                    + (PartState(5,iPartIndx))**2 &
                                                    + (PartState(6,iPartIndx))**2 ) &
                                                    * WeightFak * Species(PartSpecies(iPartIndx))%MassIC
    BulkValues(iElem)%DegreeOfFreedom = BulkValues(iElem)%DegreeOfFreedom + PartStateBulkValues(iPartIndx,5) * WeightFak
    IF(PartStateBulkValues(iPartIndx,5).EQ. 0.0) THEN
        SWRITE(UNIT_StdOut,'(132("-"))')
        SWRITE(UNIT_StdOut,'(A)') 'Element, DOF, Part:',iElem, PartStateBulkValues(iPartIndx,5), iPartIndx
        SWRITE(UNIT_StdOut,'(A)') 'T_part:',PartStateBulkValues(iPartIndx,4)
        SWRITE(UNIT_StdOut,'(A)') 'v_part:',PartStateBulkValues(iPartIndx,1:3)
        CALL abort(&
        __STAMP__&
        ,'ERROR: DOF of part is zero')
    END IF
    CellMass                          = CellMass + WeightFak * Species(PartSpecies(iPartIndx))%MassIC
    MPFSum = MPFSum + WeightFak
    IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.&
        (SpecDSMC(PartSpecies(iPartIndx))%InterID.EQ.2)) THEN
      MeanInternalEnergie = MeanInternalEnergie &
                          + PartStateIntEn(iPartIndx,1)* WeightFak &
                          + PartStateIntEn(iPartIndx,2)* WeightFak
      MPFSumOfMolecs = MPFSumOfMolecs + WeightFak
    END IF
    iPartIndx = PEM%pNext(iPartIndx)
  END DO

  BulkValues(iElem)%CellV(1) = BulkValues(iElem)%CellV(1) / CellMass
  BulkValues(iElem)%CellV(2) = BulkValues(iElem)%CellV(2) / CellMass
  BulkValues(iElem)%CellV(3) = BulkValues(iElem)%CellV(3) / CellMass
  CellVelo2                  = CellVelo2 / CellMass
  BulkValues(iElem)%DegreeOfFreedom = BulkValues(iElem)%DegreeOfFreedom / MPFSum
  IF(MPFSumOfMolecs.NE. 0.0)THEN
    MeanInternalEnergie = MeanInternalEnergie / MPFSumOfMolecs
  END IF
  BulkValues(iElem)%MassDens = CellMass / GEO%Volume(iElem)
  CellTemp = 1.0 /(BulkValues(iElem)%DegreeOfFreedom * BoltzmannConst) &
           * (2.0 * MeanInternalEnergie &
           + (CellMass / MPFSum) * (nPart/(nPart-1)) &
           * (CellVelo2 - ( BulkValues(iElem)%CellV(1)**2 &
                          + BulkValues(iElem)%CellV(2)**2 &
                          + BulkValues(iElem)%CellV(3)**2 )) )
  BulkValues(iElem)%BulkTemperature = CellTemp
  IF (CellTemp.lt. 0) THEN
      SWRITE(UNIT_StdOut,'(132("-"))')
      SWRITE(UNIT_StdOut,'(A)') 'Elem, Temperature:',iElem, PartStateBulkValues(iPartIndx,4)
      CALL abort(&
      __STAMP__&
      ,'ERROR: ERROR LD-DSMC')
  END IF
  BulkValues(iElem)%Beta = SQRT(CellMass / MPFSum / (2 * CellTemp * BoltzmannConst))
!----Ralaxationfaktor due to statistical noise in DSMC Results
!  IF(iter.NE.0) THEN ! now in LD_DSMC_Indicate_DSMC_Particles
    BulkValues(iElem)%CellV(1) =(1.0 - LD_DSMC_RelaxationFak_BufferA) * OldBulkValues(1) &
                               + LD_DSMC_RelaxationFak_BufferA * BulkValues(iElem)%CellV(1)
    BulkValues(iElem)%CellV(2) =(1.0 - LD_DSMC_RelaxationFak_BufferA) * OldBulkValues(2) &
                               + LD_DSMC_RelaxationFak_BufferA * BulkValues(iElem)%CellV(2)
    BulkValues(iElem)%CellV(3) =(1.0 - LD_DSMC_RelaxationFak_BufferA) * OldBulkValues(3) &
                               + LD_DSMC_RelaxationFak_BufferA * BulkValues(iElem)%CellV(3)
    BulkValues(iElem)%DegreeOfFreedom =(1.0 - LD_DSMC_RelaxationFak_BufferA) * OldBulkValues(4) &
                               + LD_DSMC_RelaxationFak_BufferA * BulkValues(iElem)%DegreeOfFreedom
    BulkValues(iElem)%BulkTemperature =(1.0 - LD_DSMC_RelaxationFak_BufferA) * OldBulkValues(6) &
                               + LD_DSMC_RelaxationFak_BufferA * BulkValues(iElem)%BulkTemperature
!  IF(iter.NE.0) THEN
    BulkValues(iElem)%MassDens =(1.0 - LD_DSMC_RelaxationFak_BufferA) * OldBulkValues(5) &
                               + LD_DSMC_RelaxationFak_BufferA * BulkValues(iElem)%MassDens
    BulkValues(iElem)%Beta =(1.0 - LD_DSMC_RelaxationFak_BufferA) * OldBulkValues(7) &
                               + LD_DSMC_RelaxationFak_BufferA * BulkValues(iElem)%Beta
!  ELSE
!    BulkValues(iElem)%MassDens = BulkValues(iElem)%MassDens
!    BulkValues(iElem)%Beta = SQRT(CellMass / MPFSum / (2 * OldBulkValues(6) * BoltzmannConst))
!  END IF

ELSE ! if there is no particle => keep old LD-state
    SWRITE(UNIT_StdOut,'(132("-"))')
    SWRITE(UNIT_StdOut,'(A)') 'Elem, Beta:',iElem, BulkValues(iElem)%Beta
    CALL abort(&
    __STAMP__&
    ,'ERROR: YOU NEED MORE PARTCLES FOR LD (A Buffer Zone)!!!')
END IF

END SUBROUTINE LD_DSMC_Mean_Bufferzone_A_Val


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

END MODULE MOD_LD_DSMC_TOOLS
#endif /* DONONTCOMPILETHIS*/
