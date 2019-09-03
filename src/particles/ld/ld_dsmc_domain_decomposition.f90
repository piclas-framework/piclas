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

MODULE MOD_LD_DSMC_DOMAIN_DEC
!==================================================================================================
! module for LD Sampling and Output Calculation
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE LD_DSMC_DOMAIN_DECOMPOSITION
  MODULE PROCEDURE LD_DSMC_DOMAIN_DECOMPOSITION
END INTERFACE

PUBLIC :: LD_DSMC_DOMAIN_DECOMPOSITION
!===================================================================================================================================
CONTAINS
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE LD_DSMC_DOMAIN_DECOMPOSITION
!===================================================================================================================================
! Subroutine to set celltype
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_Globals_Vars,          ONLY: BoltzmannConst
  USE MOD_Mesh_Vars,             ONLY : nElems, Elem_xGP
  USE MOD_LD_Vars,               ONLY : BulkValues,PartStateBulkValues
  USE MOD_Particle_Vars,         ONLY : PEM, PartState, usevMPF, PartMPF, Species,  PartSpecies
  USE MOD_Particle_Mesh_Vars,    ONLY : GEO
  USE MOD_DSMC_Vars,             ONLY : SpecDSMC, PartStateIntEn, CollisMode
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

  INTEGER           :: iElem, iPart, nPart, iPartIndx
  REAL              :: LocalX_Val, LocalY_Val!, FunctionVal
  REAL              :: F_1, F_2, F_3, F_4, F_5, F_6, F_7, F_8, F_9
  REAL              :: CellMass, MPFSum, CellVelo2, WeightFak, CellTemp
  REAL              :: MeanInternalEnergie, MPFSumOfMolecs
!!            /        /        /               /        /        /
!!           /        /        /               /        /        /
!!  ---LD-  F1  -B-  F2  -A-  F3  ---DSMC---  F4  -A-  F5  -B-  F6  -LD---  ...
!!          /        /        /               /        /        /
!!         /        /        /               /        /        /

!===================================================================================================================================

!CALL LD_DSMC_CBC_Calculation

DO iElem = 1, nElems
  BulkValues(iElem)%CellType = 4
  LocalX_Val = Elem_xGP(1,0,0,0,iElem)
  LocalY_Val = Elem_xGP(2,0,0,0,iElem)


!  IF(LocalY_Val.GT. 0.87) THEN
!    IF((LocalX_Val.GE. 1.15).AND.(LocalX_Val.LE. 1.19)) BulkValues(iElem)%CellType = 1
!    IF( ((LocalX_Val.GE. 1.12).AND.(LocalX_Val.LT. 1.15)) &
!   .OR. ((LocalX_Val.GT. 1.19).AND.(LocalX_Val.LE. 1.22)) ) BulkValues(iElem)%CellType = 2
!    IF( ((LocalX_Val.GE. 1.09).AND.(LocalX_Val.LT. 1.12)) &
!   .OR. ((LocalX_Val.GT. 1.22).AND.(LocalX_Val.LE. 1.25)) ) BulkValues(iElem)%CellType = 3
!  ELSE
!    FunctionVal = (LocalY_Val + 0.6232835822) / 1.298507463
!    F_1 = FunctionVal - 0.06
!    F_2 = FunctionVal - 0.03
!    F_3 = FunctionVal
!    F_4 = FunctionVal + 0.04
!    F_5 = FunctionVal + 0.07
!    F_6 = FunctionVal + 0.1
!    IF((LocalX_Val.GE. F_3).AND.(LocalX_Val.LE. F_4)) BulkValues(iElem)%CellType = 1
!    IF( ((LocalX_Val.GE. F_2).AND.(LocalX_Val.LT. F_3)) &
!   .OR. ((LocalX_Val.GT. F_4).AND.(LocalX_Val.LE. F_5)) ) BulkValues(iElem)%CellType = 2
!    IF( ((LocalX_Val.GE. F_1).AND.(LocalX_Val.LT. F_2)) &
!   .OR. ((LocalX_Val.GT. F_5).AND.(LocalX_Val.LE. F_6)) ) BulkValues(iElem)%CellType = 3
!  END IF

!!! - Tunnel Test Case
!    IF((LocalX_Val.GE. 0.004).AND.(LocalX_Val.LE. 0.006)) BulkValues(iElem)%CellType = 1
!    IF( ((LocalX_Val.GE. 0.003).AND.(LocalX_Val.LT. 0.004)) &
!   .OR. ((LocalX_Val.GT. 0.006).AND.(LocalX_Val.LE. 0.007)) ) BulkValues(iElem)%CellType = 2
!    IF( ((LocalX_Val.GE. 0.002).AND.(LocalX_Val.LT. 0.003)) &
!   .OR. ((LocalX_Val.GT. 0.007).AND.(LocalX_Val.LE. 0.008)) ) BulkValues(iElem)%CellType = 3

  F_1 = 5.340710836*LocalY_Val**2 + 0.003979362*LocalY_Val - 0.03 - 0.0025
  F_2 = 5.340710836*LocalY_Val**2 + 0.003979362*LocalY_Val - 0.0275 - 0.00125
  F_3 = 5.340710836*LocalY_Val**2 + 0.003979362*LocalY_Val - 0.025
  F_4 = 5.816923409*LocalY_Val**2 - 0.015064379*LocalY_Val - 0.015
  F_5 = 5.816923409*LocalY_Val**2 - 0.015064379*LocalY_Val - 0.0125 + 0.00125
  F_6 = 5.816923409*LocalY_Val**2 - 0.015064379*LocalY_Val - 0.01  + 0.0025
  F_7 = 8176.8*LocalY_Val**4 - 846.45*LocalY_Val**3 + 35.228*LocalY_Val**2 - 0.344*LocalY_Val + 0.015
  F_8 = 9628.9*LocalY_Val**4 - 1013.6*LocalY_Val**3 + 42.562*LocalY_Val**2 - 0.4513*LocalY_Val + 0.02
  F_9 = 10087*LocalY_Val**4 - 1001.4*LocalY_Val**3 + 40.246*LocalY_Val**2 - 0.4125*LocalY_Val + 0.025
  IF(LocalX_Val.GE.F_1) THEN
    BulkValues(iElem)%CellType = 3
    IF(LocalX_Val.GE.F_2) THEN
      BulkValues(iElem)%CellType = 2
      IF(LocalX_Val.GE.F_3) THEN
        BulkValues(iElem)%CellType = 1
        IF(LocalX_Val.GE.F_4) THEN
          BulkValues(iElem)%CellType = 2
          IF(LocalX_Val.GE.F_5) THEN
            BulkValues(iElem)%CellType = 3
            IF(LocalX_Val.GE.F_6) THEN
              BulkValues(iElem)%CellType = 4
              IF(LocalX_Val.GE.F_7) THEN
                BulkValues(iElem)%CellType = 3
                IF(LocalX_Val.GE.F_8) THEN
                  BulkValues(iElem)%CellType = 2
                  IF(LocalX_Val.GE.F_9) THEN
                    BulkValues(iElem)%CellType = 1
                  END IF
                END IF
              END IF
            END IF
          END IF
        END IF
      END IF
    END IF
  END IF

! ---adjust particle values in LD-cells and Bufferzone_B
  IF((BulkValues(iElem)%CellType.EQ.3).OR.(BulkValues(iElem)%CellType.EQ.4)) THEN ! overwrite maxwell velocities in case of LD
    nPart = PEM%pNumber(iElem)
    iPartIndx = PEM%pStart(iElem)
    DO iPart = 1, nPart
      PartState(4,iPartIndx) = PartStateBulkValues(iPartIndx,1)
      PartState(5,iPartIndx) = PartStateBulkValues(iPartIndx,2)
      PartState(6,iPartIndx) = PartStateBulkValues(iPartIndx,3)
      iPartIndx = PEM%pNext(iPartIndx)
    END DO
  END IF
! ---adjust cell values for Bufferzone_A
  IF(BulkValues(iElem)%CellType.EQ.2) THEN
    BulkValues(iElem)%CellV           = 0.0
    BulkValues(iElem)%MassDens        = 0.0
    BulkValues(iElem)%DegreeOfFreedom = 0.0
    CellMass                          = 0.0
    MPFSum                            = 0.0
    MPFSumOfMolecs                    = 0.0
    CellVelo2                         = 0.0
    MeanInternalEnergie               = 0.0
    nPart = PEM%pNumber(iElem)
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
    IF (CellTemp.LE. 0) THEN
      SWRITE(UNIT_StdOut,'(132("-"))')
      SWRITE(UNIT_stdOut,'(A)') 'Element, Temperatur:',iElem, CellTemp
CALL abort(&
__STAMP__&
,'ERROR LD-DSMC: Temperature is lt zero')
    END IF
    BulkValues(iElem)%Beta = SQRT(CellMass / MPFSum / (2 * CellTemp * BoltzmannConst))
  END IF
END DO

! MPI COMMUNICATION

!CALL LD_DSMC_SMOOTH_DOMAIN
!CALL LD_DSMC_ASSIGN_BUFFER_REGION

END SUBROUTINE LD_DSMC_DOMAIN_DECOMPOSITION

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

!SUBROUTINE LD_DSMC_CBC_Calculation
!!================================================================================================================================
!! Subroutine to calculate continuum breakdown criteria and assign initial celltyp
!!================================================================================================================================
!! MODULES
!  USE MOD_Mesh_Vars,             ONLY : nElems
!! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
!!--------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES

!  INTEGER           :: iElem
!!================================================================================================================================

!DO iElem = 1, nElems



!END DO

!END SUBROUTINE LD_DSMC_Indicate_DSMC_Particles

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

!SUBROUTINE LD_DSMC_SMOOTH_DOMAIN
!!================================================================================================================================
!! Subroutine to calculate continuum breakdown criteria and assign initial celltyp
!!================================================================================================================================
!! MODULES
!  USE MOD_Mesh_Vars,             ONLY : nElems
!! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
!!--------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES

!  INTEGER           :: iElem
!!================================================================================================================================



!END SUBROUTINE LD_DSMC_SMOOTH_DOMAIN

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

!SUBROUTINE LD_DSMC_ASSIGN_BUFFER_REGION
!!================================================================================================================================
!! Subroutine to assign Bufferregion between DSMC and LD regions
!!================================================================================================================================
!! MODULES
!  USE MOD_Mesh_Vars,             ONLY : nElems
!! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
!!--------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES

!  INTEGER           :: iElem
!!================================================================================================================================

!DO iElem = 1, nElems



!END DO

!END SUBROUTINE LD_DSMC_ASSIGN_BUFFER_REGION

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
END MODULE MOD_LD_DSMC_DOMAIN_DEC
