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

MODULE MOD_DSMC_Collis
!===================================================================================================================================
! Module including collisions, relaxation and reaction decision
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE DSMC_perform_collision
  MODULE PROCEDURE DSMC_perform_collision
END INTERFACE


INTERFACE DSMC_calc_var_P_vib
  MODULE PROCEDURE DSMC_calc_var_P_vib
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: DSMC_perform_collision
PUBLIC :: DSMC_calc_var_P_vib
!===================================================================================================================================

CONTAINS

SUBROUTINE DSMC_perform_collision(iPair, iElem, NodeVolume, NodePartNum)
!===================================================================================================================================
! Collision mode is selected (1: Elastic, 2: Non-elastic, 3: Non-elastic with chemical reactions)
!===================================================================================================================================
! MODULES
  USE MOD_Globals,            ONLY : Abort
  USE MOD_DSMC_Vars,          ONLY : CollisMode, Coll_pData, SelectionProc
  USE MOD_DSMC_Vars,          ONLY : DSMC
  USE MOD_Particle_Vars,      ONLY : PartState, WriteMacroVolumeValues, Symmetry2D
  USE MOD_TimeDisc_Vars,      ONLY : TEnd, Time
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iPair
  INTEGER, INTENT(IN)           :: iElem
  REAL, INTENT(IN), OPTIONAL    :: NodeVolume
  INTEGER, INTENT(IN), OPTIONAL :: NodePartNum
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  LOGICAL                       :: RelaxToDo
  INTEGER                       :: iPart1, iPart2                         ! Colliding particles 1 and 2
  REAL                          :: Distance
!===================================================================================================================================
  iPart1 = Coll_pData(iPair)%iPart_p1
  iPart2 = Coll_pData(iPair)%iPart_p2

  IF(DSMC%CalcQualityFactors) THEN
    IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
      IF(Symmetry2D) THEN
        Distance = SQRT((PartState(iPart1,1) - PartState(iPart2,1))**2 &
                       +(PartState(iPart1,2) - PartState(iPart2,2))**2)
      ELSE
        Distance = SQRT((PartState(iPart1,1) - PartState(iPart2,1))**2 &
                       +(PartState(iPart1,2) - PartState(iPart2,2))**2 &
                       +(PartState(iPart1,3) - PartState(iPart2,3))**2)
      END IF
      DSMC%CollSepDist = DSMC%CollSepDist + Distance
      DSMC%CollSepCount = DSMC%CollSepCount + 1
    END IF
  END IF

  SELECT CASE(CollisMode)
    CASE(1) ! elastic collision
#if (PP_TimeDiscMethod==42)
      ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
      IF (.NOT.DSMC%ReservoirSimuRate) THEN
#endif
        CALL DSMC_Elastic_Col(iPair)
#if (PP_TimeDiscMethod==42)
      END IF
#endif
    CASE(2) ! collision with relaxation
#if (PP_TimeDiscMethod==42)
      ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
      IF (.NOT.DSMC%ReservoirSimuRate) THEN
#endif
        SELECT CASE(SelectionProc)
          CASE(1)
            CALL DSMC_Relax_Col_LauxTSHO(iPair)
          CASE(2)
            CALL DSMC_Relax_Col_Gimelshein(iPair)
          CASE DEFAULT
            CALL Abort(&
__STAMP__&
,'ERROR in DSMC_perform_collision: Wrong Selection Procedure:',SelectionProc)
        END SELECT
#if (PP_TimeDiscMethod==42)
      END IF
#endif
    CASE(3) ! chemical reactions
      RelaxToDo = .TRUE.
      IF (PRESENT(NodeVolume).AND.PRESENT(NodePartNum)) THEN
        CALL ReactionDecision(iPair, RelaxToDo, iElem, NodeVolume, NodePartNum)
      ELSE
        CALL ReactionDecision(iPair, RelaxToDo, iElem)
      END IF
#if (PP_TimeDiscMethod==42)
      ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
      IF (.NOT.DSMC%ReservoirSimuRate) THEN
#endif
        IF (RelaxToDo) THEN
          SELECT CASE(SelectionProc)
            CASE(1)
              CALL DSMC_Relax_Col_LauxTSHO(iPair)
            CASE(2)
              CALL DSMC_Relax_Col_Gimelshein(iPair)
            CASE DEFAULT
              CALL Abort(&
__STAMP__&
,'ERROR in DSMC_perform_collision: Wrong Selection Procedure:',SelectionProc)
          END SELECT
        END IF
#if (PP_TimeDiscMethod==42)
      END IF
#endif
    CASE DEFAULT
      CALL Abort(&
__STAMP__&
,'ERROR in DSMC_perform_collision: Wrong Collision Mode:',CollisMode)
  END SELECT

END SUBROUTINE DSMC_perform_collision


SUBROUTINE DSMC_Elastic_Col(iPair)
!===================================================================================================================================
! Performs simple elastic collision (CollisMode = 1)
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : Coll_pData, CollInf, DSMC_RHS, RadialWeighting
  USE MOD_Particle_Vars,          ONLY : PartSpecies, PartState, VarTimeStep, Species
  USE MOD_part_tools,             ONLY : DiceDeflectedVelocityVector
  USE MOD_part_tools              ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iPair
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                          :: FracMassCent1, FracMassCent2     ! mx/(mx+my)
  REAL                          :: VeloMx, VeloMy, VeloMz           ! center of mass velo
  REAL                          :: cRelaX, cRelaY, cRelaZ           ! pre-collision relative velocities
  REAL                          :: RanVelo(3)                       ! random relative velocity
  INTEGER                       :: iPart1, iPart2, iSpec1, iSpec2   ! Colliding particles 1 and 2, their species
!===================================================================================================================================
  iPart1 = Coll_pData(iPair)%iPart_p1
  iPart2 = Coll_pData(iPair)%iPart_p2

  iSpec1 = PartSpecies(iPart1)
  iSpec2 = PartSpecies(iPart2)

  IF (RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
    FracMassCent1 = Species(iSpec1)%MassIC * GetParticleWeight(iPart1) / (Species(iSpec1)%MassIC &
                  * GetParticleWeight(iPart1) + Species(iSpec2)%MassIC *GetParticleWeight(iPart2))
    FracMassCent2 = Species(iSpec2)%MassIC *GetParticleWeight(iPart2) / (Species(iSpec1)%MassIC  &
                  * GetParticleWeight(iPart1) + Species(iSpec2)%MassIC *GetParticleWeight(iPart2))
  ELSE
    FracMassCent1 = CollInf%FracMassCent(PartSpecies(iPart1), Coll_pData(iPair)%PairType)
    FracMassCent2 = CollInf%FracMassCent(PartSpecies(iPart2), Coll_pData(iPair)%PairType)
  END IF
  !Calculation of center of mass velocity
  VeloMx = FracMassCent1 * PartState(iPart1, 4) + FracMassCent2 * PartState(iPart2, 4)
  VeloMy = FracMassCent1 * PartState(iPart1, 5) + FracMassCent2 * PartState(iPart2, 5)
  VeloMz = FracMassCent1 * PartState(iPart1, 6) + FracMassCent2 * PartState(iPart2, 6)

  IF (CollInf%alphaVSS(iSpec1,iSpec2).GT.1) THEN
    ! Calculation of relative velocities
    cRelax = PartState(iPart1, 4) - PartState(iPart2, 4)
    cRelay = PartState(iPart1, 5) - PartState(iPart2, 5)
    cRelaz = PartState(iPart1, 6) - PartState(iPart2, 6)

    ! Calculation of post collision velocity vector in reference frame and retransformation to center-of-mass frame
    RanVelo(1:3) = DiceDeflectedVelocityVector(Coll_pData(iPair)%cRela2 , CollInf%alphaVSS(iSpec1,iSpec2), cRelaX , cRelaY , cRelaZ)
  ELSE ! alphaVSS .LE. 1
    ! Calculation of post collision velocity vector in reference frame and retransformation to center-of-mass frame
    RanVelo(1:3) = DiceDeflectedVelocityVector(Coll_pData(iPair)%cRela2,CollInf%alphaVSS(iSpec1,iSpec2))
  END IF  ! alphaVSS 

  !    WRITE(*,*) "DiceDeflectedVector in Elastic ",RanVelo/sqrt(Coll_pData(iPair)%cRela2) !to be solved

 ! deltaV particle 1 (post collision particle 1 velocity in laboratory frame)
  DSMC_RHS(iPart1,1) = VeloMx + FracMassCent2 * RanVelo(1) - PartState(iPart1, 4)
  DSMC_RHS(iPart1,2) = VeloMy + FracMassCent2 * RanVelo(2) - PartState(iPart1, 5)
  DSMC_RHS(iPart1,3) = VeloMz + FracMassCent2 * RanVelo(3) - PartState(iPart1, 6)
 ! deltaV particle 2 (post collision particle 2 velocity in laboratory frame)
  DSMC_RHS(iPart2,1) = VeloMx - FracMassCent1 * RanVelo(1) - PartState(iPart2, 4)
  DSMC_RHS(iPart2,2) = VeloMy - FracMassCent1 * RanVelo(2) - PartState(iPart2, 5)
  DSMC_RHS(iPart2,3) = VeloMz - FracMassCent1 * RanVelo(3) - PartState(iPart2, 6)

END SUBROUTINE DSMC_Elastic_Col

SUBROUTINE DSMC_Scat_Col(iPair)
!===================================================================================================================================
! Performs a collision with the possibility of a CEX. In the calculation of the new particle velocities a scattering angle is used,
! which is interpolated from a lookup table.
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : Coll_pData, CollInf, DSMC_RHS, TLU_Data, ChemReac
  USE MOD_Particle_Vars,          ONLY : PartSpecies, PartState
  USE MOD_DSMC_ChemReact,         ONLY : simpleCEX, simpleMEX

! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iPair
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                          :: FracMassCent1, FracMassCent2                ! mx/(mx+my)
  REAL                          :: VeloMx, VeloMy, VeloMz                      ! center of mass velo
  REAL                          :: cRelax, cRelay, cRelaz                      ! pre-collisional relativ velo
  REAL                          :: cRelaxN, cRelayN, cRelazN                   ! post-collisional relativ velo
  REAL                          :: b, bmax                                     ! impact parameters
  REAL                          :: Ekin
  REAL                          :: ScatAngle, RotAngle                         ! scattering and rotational angle
  REAL                          :: sigma_el, sigma_tot                         ! cross-sections
  REAL                          :: P_CEX                                       ! charge exchange probability
  INTEGER                       :: iReac
  REAL                          :: uRan2, uRan3, uRanRot, uRanVHS
  REAL                          :: Pi, aEL, bEL, aCEX, bCEX
  INTEGER                       :: iPart1, iPart2                    ! Colliding particles 1 and 2
!===================================================================================================================================
 iPart1 = Coll_pData(iPair)%iPart_p1
 iPart2 = Coll_pData(iPair)%iPart_p2

  Pi = ACOS(-1.0)
  aCEX = ChemReac%CEXa(ChemReac%ReactNum(PartSpecies(iPart1),PartSpecies(iPart2),1))
  bCEX = ChemReac%CEXb(ChemReac%ReactNum(PartSpecies(iPart1),PartSpecies(iPart2),1))
  aEL  = ChemReac%ELa(ChemReac%ReactNum(PartSpecies(iPart1),PartSpecies(iPart2),1))
  bEL  = ChemReac%ELb(ChemReac%ReactNum(PartSpecies(iPart1),PartSpecies(iPart2),1))
  ! Decision if scattering angle is greater than 1 degree and should be calculated

  sigma_el  = bEL + aEL*0.5 * LOG10(Coll_pData(iPair)%cRela2)

  sigma_tot = ((aCEX+0.5*aEL)*0.5*LOG10(Coll_pData(iPair)%cRela2)+bCEX+0.5*bEL)

  CALL RANDOM_NUMBER(uRan2)

IF ((sigma_el/sigma_tot).GT.uRan2) THEN
    ! Calculation of relative velocities
    cRelax = PartState(iPart1, 4) - PartState(iPart2, 4)
    cRelay = PartState(iPart1, 5) - PartState(iPart2, 5)
    cRelaz = PartState(iPart1, 6) - PartState(iPart2, 6)

    FracMassCent1 = CollInf%FracMassCent(PartSpecies(iPart1), Coll_pData(iPair)%PairType)
    FracMassCent2 = CollInf%FracMassCent(PartSpecies(iPart2), Coll_pData(iPair)%PairType)

    ! Calculation of velo from center of mass
    VeloMx = FracMassCent1 * PartState(iPart1, 4) + FracMassCent2 * PartState(iPart2, 4)
    VeloMy = FracMassCent1 * PartState(iPart1, 5) + FracMassCent2 * PartState(iPart2, 5)
    VeloMz = FracMassCent1 * PartState(iPart1, 6) + FracMassCent2 * PartState(iPart2, 6)

    ! Calculation of impact parameter b
    bmax = SQRT(sigma_el/Pi)
    b = bmax * SQRT(uRan2)
    Ekin = (0.5*CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%cRela2/(1.6021766208E-19))


    ! Determination of scattering angle by interpolation from a lookup table
    ! Check if Collision Energy is below the threshold of table
    IF (Ekin.LT.TLU_Data%Emin) THEN
      ! Isotropic scattering
      CALL RANDOM_NUMBER(uRanVHS)
      ScatAngle = 2*ACOS(SQRT(uRanVHS))
    ELSE
      ! scattering corresponding to table lookup
      CALL TLU_Scat_Interpol(Ekin,b,ScatAngle)
    END IF

    ! Determination of rotational angle by random number
    CALL RANDOM_NUMBER(uRanRot)
    RotAngle = uRanRot * 2 * Pi

    ! Calculation of post-collision relative velocities in center-of-mass frame
    cRelaxN = COS(ScatAngle)*cRelax + SIN(ScatAngle)*SIN(RotAngle)*(cRelay**2+cRelaz**2)**0.5
    cRelayN = COS(ScatAngle)*cRelay &
     +SIN(ScatAngle)*(SQRT(Coll_pData(ipair)%cRela2)*cRelaz*COS(RotAngle)-cRelax*cRelay*SIN(RotAngle))/(cRelay**2+cRelaz**2)**0.5
    cRelazN = COS(ScatAngle)*cRelaz &
     -SIN(ScatAngle)*(SQRT(Coll_pData(ipair)%cRela2)*cRelay*COS(RotAngle)+cRelax*cRelaz*SIN(RotAngle))/(cRelay**2+cRelaz**2)**0.5

    ! Transformation to laboratory frame
    ! deltaV particle 1
    DSMC_RHS(iPart1,1) = VeloMx + FracMassCent2*cRelaxN - PartState(iPart1, 4)
    DSMC_RHS(iPart1,2) = VeloMy + FracMassCent2*cRelayN - PartState(iPart1, 5)
    DSMC_RHS(iPart1,3) = VeloMz + FracMassCent2*cRelazN - PartState(iPart1, 6)
    ! deltaV particle 2
    DSMC_RHS(iPart2,1) = VeloMx - FracMassCent1*cRelaxN - PartState(iPart2, 4)
    DSMC_RHS(iPart2,2) = VeloMy - FracMassCent1*cRelayN - PartState(iPart2, 5)
    DSMC_RHS(iPart2,3) = VeloMz - FracMassCent1*cRelazN - PartState(iPart2, 6)

    ! Decision concerning CEX
    P_CEX = 0.5
    CALL RANDOM_NUMBER(uRan3)
    iReac    = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 1)
    IF (P_CEX.GT.uRan3) THEN
      CALL simpleCEX(iReac, iPair, resetRHS_opt=.FALSE.)
    ELSE
      CALL simpleMEX(iReac, iPair)
    END IF

  ELSE
    ! Perform CEX and leave velocity vectors alone otherwise
    ! CEX
    iReac    = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 1)
    CALL simpleCEX(iReac, iPair)

  END IF

END SUBROUTINE DSMC_Scat_Col

SUBROUTINE TLU_Scat_Interpol(E_p,b_p,ScatAngle)
!===================================================================================================================================
! Interpolates ScatAngle from a lookup table
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_DSMC_Vars,              ONLY :  TLU_Data
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT (IN)              :: E_p, b_p          ! E_p has to have the unit eV
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL, INTENT (OUT)             :: ScatAngle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                           :: i_f_jp1, j_f, i_f_j
  INTEGER                        :: I_j,I_jp1,J
  REAL                           :: w_i_j,w_i_jp1,w_j
  INTEGER                        :: szb,szE
  REAL                           :: chi_b_p_E_j,chi_b_p_E_jp1,chi_b_p_e_p
!===================================================================================================================================
  IF (E_p.GT.TLU_Data%Emax) THEN
    CALL abort(__STAMP__,&
        'Collis_mode - Error in TLU_Scat_Interpol: E_p GT Emax')
  END IF
  !write (*,*) (E_p-TLU_Data%Emin), TLU_Data%deltaE
  j_f = (E_p-TLU_Data%Emin)/TLU_Data%deltaE
  J = FLOOR(j_f)
  w_j = j_f - J
  J = J + 1                                ! Fitting of the indices for the use in FORTRAN matrix
  !write (*,*) j_f, J, w_j
  i_f_j   = ABS((b_p)/TLU_Data%deltabj(J))
  i_f_jp1 = ABS((b_p)/TLU_Data%deltabj(J+1))
  I_j     = FLOOR(i_f_j)
  I_jp1   = FLOOR(i_f_jp1)

  w_i_j = i_f_j - I_j
  w_i_jp1 = i_f_jp1-I_jp1

  I_j     = FLOOR(i_f_j)+1                ! Fitting of the indices for the use in FORTRAN matrix
  I_jp1   = FLOOR(i_f_jp1)+1              !

  szE = SIZE(TLU_Data%Chitable,dim=1)   !SIZE(delta_b_j)
  szB = SIZE(TLU_Data%Chitable,dim=2)



  IF ((I_jp1+1).GE.szB) THEN
    chi_b_p_E_j   = (1 - w_i_j) * TLU_Data%Chitable(J,szB)       !+ w_i_j   * TLU_Data%Chitable(J,szB)
    chi_b_p_E_jp1 = (1-w_i_jp1) * TLU_Data%Chitable((J+1),szB)
    chi_b_p_E_p   = (1-w_j)     * chi_b_p_E_j                    + w_j     * chi_b_p_E_jp1
  ELSE
    chi_b_p_E_j   = (1 - w_i_j) * TLU_Data%Chitable(J,I_j)       + w_i_j   * TLU_Data%Chitable(J,I_jp1)
    chi_b_p_E_jp1 = (1-w_i_jp1) * TLU_Data%Chitable((J+1),I_jp1) + w_i_jp1 * TLU_Data%Chitable((J+1),(I_jp1+1))
    chi_b_p_E_p   = (1-w_j)     * chi_b_p_E_j                    + w_j     * chi_b_p_E_jp1
  END IF
  ScatAngle = chi_b_p_E_p

  !write(*,*) (ScatAngle/ACOS(-1.0)*180), I_jp1, szB
END SUBROUTINE TLU_Scat_Interpol

SUBROUTINE DSMC_Relax_Col_LauxTSHO(iPair)
!===================================================================================================================================
! Performs inelastic collisions with energy exchange (CollisMode = 2/3), allows the relaxation of both collision partners
! Vibrational (of the relaxing molecule), rotational and relative translational energy (of both molecules) is redistributed (V-R-T)
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : Coll_pData, CollInf, DSMC_RHS, DSMC, &
                                         SpecDSMC, PartStateIntEn, RadialWeighting
  USE MOD_Globals_Vars,           ONLY : BoltzmannConst
  USE MOD_Particle_Vars,          ONLY : PartSpecies, PartState, Species, VarTimeStep, PEM
  USE MOD_DSMC_ElectronicModel,   ONLY : ElectronicEnergyExchange, TVEEnergyExchange
  USE MOD_DSMC_PolyAtomicModel,   ONLY : DSMC_RotRelaxPoly, DSMC_VibRelaxPoly
  USE MOD_DSMC_Relaxation,        ONLY : DSMC_VibRelaxDiatomic
  USE MOD_part_tools,             ONLY : DiceDeflectedVelocityVector
USE MOD_part_tools                ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iPair
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                          :: FracMassCent1, FracMassCent2     ! mx/(mx+my)
  REAL                          :: VeloMx, VeloMy, VeloMz           ! center of mass velo
  REAL (KIND=8)                 :: iRan                             ! Random number
  LOGICAL                       :: DoRot1, DoRot2, DoVib1, DoVib2   ! Check whether rot or vib relax is performed
  REAL                          :: cRelaX, cRelaY, cRelaZ           ! pre-collision relative velocities
  REAL (KIND=8)                 :: Xi_rel, Xi, FakXi                ! Factors of DOF
  REAL                          :: RanVelo(3)                       ! random relative velocity
  REAL                          :: ReducedMass
  REAL                          :: ProbRot1, ProbRotMax1, ProbRot2, ProbRotMax2, ProbVib1, ProbVib2
  INTEGER                       :: iSpec1, iSpec2, iPart1, iPart2, iElem ! Colliding particles 1 and 2 and their species
  ! variables for electronic level relaxation and transition
  LOGICAL                       :: DoElec1, DoElec2
!===================================================================================================================================
  iPart1 = Coll_pData(iPair)%iPart_p1
  iPart2 = Coll_pData(iPair)%iPart_p2
  iSpec1 = PartSpecies(iPart1)
  iSpec2 = PartSpecies(iPart2)
  iElem  = PEM%Element(iPart1)

  DoRot1  = .FALSE.
  DoRot2  = .FALSE.
  DoVib1  = .FALSE.
  DoVib2  = .FALSE.
  DoElec1 = .FALSE.
  DoElec2 = .FALSE.

  IF (RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
    ReducedMass = (Species(iSpec1)%MassIC*GetParticleWeight(iPart1) * Species(iSpec2)%MassIC*GetParticleWeight(iPart2)) &
                / (Species(iSpec1)%MassIC*GetParticleWeight(iPart1) + Species(iSpec2)%MassIC*GetParticleWeight(iPart2))
  ELSE
    ReducedMass = CollInf%MassRed(Coll_pData(iPair)%PairType)
  END IF

  Xi_rel = 2*(2. - CollInf%omegaLaux(iSpec1,iSpec2))
    ! DOF of relative motion in VHS model

  Coll_pData(iPair)%Ec = 0.5 * ReducedMass* Coll_pData(iPair)%cRela2

  Xi = Xi_rel !Xi are all DOF in the collision

!--------------------------------------------------------------------------------------------------!
! Decision if Rotation, Vibration and Electronic Relaxation of particles is performed
!--------------------------------------------------------------------------------------------------!

  IF((SpecDSMC(iSpec1)%InterID.EQ.2).OR.(SpecDSMC(iSpec1)%InterID.EQ.20)) THEN
    CALL RANDOM_NUMBER(iRan)
    CALL DSMC_calc_P_rot(iSpec1, iPair, Coll_pData(iPair)%iPart_p1, Xi_rel, ProbRot1, ProbRotMax1)
    IF(ProbRot1.GT.iRan) THEN
      DoRot1 = .TRUE.
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(iPart1,2) * GetParticleWeight(iPart1)
      Xi = Xi + SpecDSMC(iSpec1)%Xi_Rot
    END IF
    IF(DSMC%CalcQualityFactors.AND.(DSMC%RotRelaxProb.GE.2)) THEN
      DSMC%CalcRotProb(iSpec1,2) = MAX(DSMC%CalcRotProb(iSpec1,2),ProbRot1)
      DSMC%CalcRotProb(iSpec1,1) = DSMC%CalcRotProb(iSpec1,1) + ProbRot1
      DSMC%CalcRotProb(iSpec1,3) = DSMC%CalcRotProb(iSpec1,3) + 1
    END IF

    CALL RANDOM_NUMBER(iRan)
    CALL DSMC_calc_P_vib(iSpec1, Xi_rel, iElem, ProbVib1)
    IF(ProbVib1.GT.iRan) DoVib1 = .TRUE.
    IF(DSMC%CalcQualityFactors.AND.(DSMC%VibRelaxProb.EQ.2)) THEN
      DSMC%CalcVibProb(iSpec1,1) = DSMC%CalcVibProb(iSpec1,1) + ProbVib1
      DSMC%CalcVibProb(iSpec1,3) = DSMC%CalcVibProb(iSpec1,3) + 1
    END IF
  END IF
  IF ( DSMC%ElectronicModel ) THEN
    IF((SpecDSMC(iSpec1)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec1)%FullyIonized)) THEN
      CALL RANDOM_NUMBER(iRan)
      IF (SpecDSMC(iSpec1)%ElecRelaxProb.GT.iRan) THEN
        DoElec1 = .TRUE.
      END IF
    END IF
  END IF

  IF((SpecDSMC(iSpec2)%InterID.EQ.2).OR.(SpecDSMC(iSpec2)%InterID.EQ.20)) THEN
    CALL RANDOM_NUMBER(iRan)
    CALL DSMC_calc_P_rot(iSpec2, iPair, Coll_pData(iPair)%iPart_p2, Xi_rel, ProbRot2, ProbRotMax2)
    IF(ProbRot2.GT.iRan) THEN
      DoRot2 = .TRUE.
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(iPart2,2) * GetParticleWeight(iPart2)
      Xi = Xi + SpecDSMC(iSpec2)%Xi_Rot
    END IF
    IF(DSMC%CalcQualityFactors.AND.(DSMC%RotRelaxProb.GE.2)) THEN
      DSMC%CalcRotProb(iSpec2,2) = MAX(DSMC%CalcRotProb(iSpec2,2),ProbRot2)
      DSMC%CalcRotProb(iSpec2,1) = DSMC%CalcRotProb(iSpec2,1) + ProbRot2
      DSMC%CalcRotProb(iSpec2,3) = DSMC%CalcRotProb(iSpec2,3) + 1
    END IF
    CALL RANDOM_NUMBER(iRan)
    CALL DSMC_calc_P_vib(iSpec2, Xi_rel, iElem, ProbVib2)
    IF(ProbVib2.GT.iRan) DoVib2 = .TRUE.
    IF(DSMC%CalcQualityFactors.AND.(DSMC%VibRelaxProb.EQ.2)) THEN
      DSMC%CalcVibProb(iSpec2,1) = DSMC%CalcVibProb(iSpec2,1) + ProbVib2
      DSMC%CalcVibProb(iSpec2,3) = DSMC%CalcVibProb(iSpec2,3) + 1
    END IF
  END IF
  IF ( DSMC%ElectronicModel ) THEN
    IF((SpecDSMC(iSpec2)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec2)%FullyIonized)) THEN
      CALL RANDOM_NUMBER(iRan)
      IF (SpecDSMC(iSpec2)%ElecRelaxProb.GT.iRan) THEN
        DoElec2 = .TRUE.
      END IF
    END IF
  END IF

  FakXi = 0.5*Xi  - 1  ! exponent factor of DOF, substitute of Xi_c - Xi_vib, laux diss page 40


!--------------------------------------------------------------------------------------------------!
! Electronic Relaxation / Transition
!--------------------------------------------------------------------------------------------------!
  IF (DSMC%DoTEVRRelaxation) THEN
    IF(.NOT.SpecDSMC(iSpec1)%PolyatomicMol) THEN
      IF(DoElec1.AND.DoVib1) THEN
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(iPart1,3)  + PartStateIntEn(iPart1,1)
        CALL TVEEnergyExchange(Coll_pData(iPair)%Ec,iPart1,FakXi)
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(iPart1,3)  - PartStateIntEn(iPart1,1)
        DoElec1=.false.
        DoVib1=.false.
      END IF !DoElec1.AND.DoVib1
    END IF ! .NOT.SpecDSMC(iSpec1)%PolyatomicMol

    IF(.NOT.SpecDSMC(iSpec2)%PolyatomicMol) THEN
      IF(DoElec2.AND.DoVib2) THEN
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(Coll_pData(iPair)%iPart_p2,3)  &
                             +    PartStateIntEn(Coll_pData(iPair)%iPart_p2,1)
        CALL TVEEnergyExchange(Coll_pData(iPair)%Ec,Coll_pData(iPair)%iPart_p2,FakXi)
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(Coll_pData(iPair)%iPart_p2,3)  &
                             -    PartStateIntEn(Coll_pData(iPair)%iPart_p2,1)
        DoElec2=.false.
        DoVib2=.false.
      END IF ! .NOT.SpecDSMC(iSpec2)%PolyatomicMol
    END IF !DoElec2.AND.DoVib2
  END IF ! DSMC%DoTEVRRelaxation

  ! Relaxation of first particle
  IF ( DoElec1 ) THEN
    ! calculate energy for electronic relaxation of particle 1
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(iPart1,3) * GetParticleWeight(iPart1)
    CALL ElectronicEnergyExchange(iPair,Coll_pData(iPair)%iPart_p1,FakXi)
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(iPart1,3) * GetParticleWeight(iPart1)
  END IF

  ! Electronic relaxation of second particle
  IF ( DoElec2 ) THEN
    ! calculate energy for electronic relaxation of particle 2
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(iPart2,3) * GetParticleWeight(iPart2)
    CALL ElectronicEnergyExchange(iPair,Coll_pData(iPair)%iPart_p2,FakXi)
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(iPart2,3) * GetParticleWeight(iPart2)
  END IF

#if (PP_TimeDiscMethod==42)
  ! for TimeDisc 42 & only transition counting: prohibit relaxation and energy exchange
  IF (.NOT.DSMC%ReservoirSimuRate) THEN
#endif

!--------------------------------------------------------------------------------------------------!
! Vibrational Relaxation
!--------------------------------------------------------------------------------------------------!

  IF(DoVib1) THEN
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(iPart1,1) * GetParticleWeight(iPart1)
    IF(SpecDSMC(iSpec1)%PolyatomicMol) THEN
      CALL DSMC_VibRelaxPoly(iPair, Coll_pData(iPair)%iPart_p1,FakXi)
    ELSE
      CALL DSMC_VibRelaxDiatomic(iPair, iPart1,FakXi)
    END IF
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(iPart1,1) * GetParticleWeight(iPart1)
  END IF

  IF(DoVib2) THEN
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(iPart2,1) * GetParticleWeight(iPart2)
    IF(SpecDSMC(iSpec2)%PolyatomicMol) THEN
      CALL DSMC_VibRelaxPoly(iPair, Coll_pData(iPair)%iPart_p2,FakXi)
    ELSE
      CALL DSMC_VibRelaxDiatomic(iPair, iPart2,FakXi)
    END IF
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(iPart2,1) * GetParticleWeight(iPart2)
  END IF

!--------------------------------------------------------------------------------------------------!
! Rotational Relaxation
!--------------------------------------------------------------------------------------------------!
  IF(DoRot1) THEN
    IF(SpecDSMC(iSpec1)%PolyatomicMol.AND.(SpecDSMC(iSpec1)%Xi_Rot.EQ.3)) THEN
      FakXi = FakXi - 0.5*SpecDSMC(iSpec1)%Xi_Rot
      CALL DSMC_RotRelaxPoly(iPair, Coll_pData(iPair)%iPart_p1, FakXi)
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(Coll_pData(iPair)%iPart_p1,2)
    ELSE
      CALL RANDOM_NUMBER(iRan)
      PartStateIntEn(Coll_pData(iPair)%iPart_p1,2) = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi))
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(Coll_pData(iPair)%iPart_p1,2)
      FakXi = FakXi - 0.5*SpecDSMC(iSpec1)%Xi_Rot
    END IF
    IF(RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
      PartStateIntEn(iPart1,2) = PartStateIntEn(iPart1,2)/GetParticleWeight(iPart1)
    END IF
  END IF

  IF(DoRot2) THEN
    IF(SpecDSMC(iSpec2)%PolyatomicMol.AND. &
        (SpecDSMC(iSpec2)%Xi_Rot.EQ.3)) THEN
      FakXi = FakXi - 0.5*SpecDSMC(iSpec2)%Xi_Rot
      CALL DSMC_RotRelaxPoly(iPair, Coll_pData(iPair)%iPart_p2, FakXi)
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(Coll_pData(iPair)%iPart_p2,2)
    ELSE
      CALL RANDOM_NUMBER(iRan)
      PartStateIntEn(Coll_pData(iPair)%iPart_p2,2) = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi))
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(Coll_pData(iPair)%iPart_p2,2)
    END IF
    IF(RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
      PartStateIntEn(iPart2,2) = PartStateIntEn(iPart2,2)/GetParticleWeight(iPart2)
    END IF
  END IF

!--------------------------------------------------------------------------------------------------!
! Calculation of new particle velocities
!--------------------------------------------------------------------------------------------------!

  IF (RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
    FracMassCent1 = Species(iSpec1)%MassIC *GetParticleWeight(iPart1)/(Species(iSpec1)%MassIC *GetParticleWeight(iPart1) &
          + Species(iSpec2)%MassIC *GetParticleWeight(iPart2))
    FracMassCent2 = Species(iSpec2)%MassIC *GetParticleWeight(iPart2)/(Species(iSpec1)%MassIC *GetParticleWeight(iPart1) &
          + Species(iSpec2)%MassIC *GetParticleWeight(iPart2))
  ELSE
    FracMassCent1 = CollInf%FracMassCent(iSpec1, Coll_pData(iPair)%PairType)
    FracMassCent2 = CollInf%FracMassCent(iSpec2, Coll_pData(iPair)%PairType)
  END IF

  !Calculation of velo from center of mass
  VeloMx = FracMassCent1 * PartState(iPart1, 4) + FracMassCent2 * PartState(iPart2, 4)
  VeloMy = FracMassCent1 * PartState(iPart1, 5) + FracMassCent2 * PartState(iPart2, 5)
  VeloMz = FracMassCent1 * PartState(iPart1, 6) + FracMassCent2 * PartState(iPart2, 6)

  IF (CollInf%alphaVSS(iSpec1,iSpec2).GT.1) THEN
    !Calculate relative velocities and new squared velocity
    cRelax = PartState(iPart1, 4) - PartState(iPart2, 4)
    cRelay = PartState(iPart1, 5) - PartState(iPart2, 5)
    cRelaz = PartState(iPart1, 6) - PartState(iPart2, 6)

    Coll_pData(iPair)%cRela2 = 2 * Coll_pData(iPair)%Ec/ReducedMass
    
    ! Calculation of post collision velocity vector in reference frame and retransformation to center-of-mass frame
    RanVelo(1:3) = DiceDeflectedVelocityVector(Coll_pData(iPair)%cRela2 , CollInf%alphaVSS(iSpec1,iSpec2), cRelaX , cRelaY , cRelaZ)

  ELSE ! alphaVSS .LE. 1
    ! Calculation of post collision velocity vector in reference frame 
    RanVelo(1:3) = DiceDeflectedVelocityVector(Coll_pData(iPair)%cRela2,CollInf%alphaVSS(iSpec1,iSpec2))
  END IF  ! alphaVSS 
  
!  WRITE(*,*) "DiceDeflectedVector in Relax LauxTHSO",RanVelo/sqrt(Coll_pData(iPair)%cRela2) !to be solved

  ! deltaV particle 1 (post collision particle 1 velocity in laboratory frame)
  DSMC_RHS(iPart1,1) = VeloMx + FracMassCent2*RanVelo(1) - PartState(iPart1, 4)
  DSMC_RHS(iPart1,2) = VeloMy + FracMassCent2*RanVelo(2) - PartState(iPart1, 5)
  DSMC_RHS(iPart1,3) = VeloMz + FracMassCent2*RanVelo(3) - PartState(iPart1, 6)
  ! deltaV particle 2 (post collision particle 2 velocity in laboratory frame)
  DSMC_RHS(iPart2,1) = VeloMx - FracMassCent1*RanVelo(1) - PartState(iPart2, 4)
  DSMC_RHS(iPart2,2) = VeloMy - FracMassCent1*RanVelo(2) - PartState(iPart2, 5)
  DSMC_RHS(iPart2,3) = VeloMz - FracMassCent1*RanVelo(3) - PartState(iPart2, 6)

#if (PP_TimeDiscMethod==42)
  ! for TimeDisc 42 & only transition counting: prohibit relaxation and energy exchange
  END IF
#endif

END SUBROUTINE DSMC_Relax_Col_LauxTSHO


SUBROUTINE DSMC_Relax_Col_Gimelshein(iPair)
!===================================================================================================================================
! Performs inelastic collisions with energy exchange (CollisMode = 2/3)
! Selection procedure according to Gimelshein et al. (Physics of Fluids, V 14, No 12, 2002: 'Vibrational Relaxation rates in the
! DSMC Method') For further understanding see Zhang, Schwarzentruber, Physics of Fluids, V25, 2013: 'inelastic collision selection
! procedures for DSMC calculation of gas mixtures')
!===================================================================================================================================
! MODULES
  USE MOD_Globals,                ONLY : Abort
  USE MOD_Globals_Vars,           ONLY : BoltzmannConst
  USE MOD_DSMC_Vars,              ONLY : Coll_pData, CollInf, DSMC_RHS, DSMC, PolyatomMolDSMC, VibQuantsPar, &
                                         SpecDSMC, PartStateIntEn
  USE MOD_Particle_Vars,          ONLY : PartSpecies, PartState, PEM
  USE MOD_DSMC_PolyAtomicModel,   ONLY : DSMC_RotRelaxPoly, DSMC_VibRelaxPoly
  USE MOD_DSMC_Relaxation,        ONLY : DSMC_VibRelaxDiatomic
  USE MOD_part_tools,             ONLY : DiceDeflectedVelocityVector
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iPair
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                          :: FracMassCent1, FracMassCent2                 ! mx/(mx+my)
  REAL                          :: VeloMx, VeloMy, VeloMz                       ! center of mass velo
  INTEGER                       :: iDOF, iPolyatMole, DOFRelax, iElem
  REAL (KIND=8)                 :: iRan
  LOGICAL                       :: DoRot1, DoRot2, DoVib1, DoVib2               ! Check whether rot or vib relax is performed
  REAL (KIND=8)                 :: FakXi, Xi_rel                                ! Factors of DOF
  INTEGER                       :: iQuaMax, iQua                                ! Quantum Numbers
  REAL                          :: MaxColQua, RanVelo(3)                        ! Max. Quantum Number
  REAL                          :: PartStateIntEnTemp                           ! temp. var for inertial energy (needed for vMPF)
  REAL                          :: ProbFrac1, ProbFrac2, ProbFrac3, ProbFrac4   ! probability-fractions according to Zhang
  REAL                          :: ProbRot1, ProbRot2, ProbVib1, ProbVib2       ! probabilities for rot-/vib-relax for part 1/2
  REAL                          :: BLCorrFact, ProbRotMax1, ProbRotMax2         ! Correction factor for BL-redistribution of energy
  REAL                          :: cRelaX, cRelaY, cRelaZ                       ! pre-collision relative velocities
  INTEGER                       :: iPart1, iPart2, iSpec1, iSpec2               ! Colliding particles 1 and 2 and their species
!===================================================================================================================================
  iPart1 = Coll_pData(iPair)%iPart_p1
  iPart2 = Coll_pData(iPair)%iPart_p2
  iSpec1 = PartSpecies(iPart1)
  iSpec2 = PartSpecies(iPart2)
  iElem  = PEM%Element(iPart1)

  ! set some initial values
  DoRot1  = .FALSE.
  DoRot2  = .FALSE.
  DoVib1  = .FALSE.
  DoVib2  = .FALSE.
  ProbVib1 = 0.
  ProbRot1 = 0.
  ProbVib2 = 0.
  ProbRot2 = 0.

  Xi_rel = 2.*(2. - CollInf%omegaLaux(iSpec1,iSpec2)) ! DOF of relative motion in VHS model
  FakXi  = 0.5*Xi_rel - 1.

  Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType) * Coll_pData(iPair)%cRela2

!--------------------------------------------------------------------------------------------------!
! Decision if Rotation, Vibration and Electronic Relaxation of particles is performed
!--------------------------------------------------------------------------------------------------!

  ! calculate probability for rotational/vibrational relaxation for both particles
  IF ((SpecDSMC(iSpec1)%InterID.EQ.2).OR.(SpecDSMC(iSpec1)%InterID.EQ.20)) THEN
    CALL DSMC_calc_P_vib(iSpec1, Xi_rel, iElem, ProbVib1)
    CALL DSMC_calc_P_rot(iSpec1, iPair, Coll_pData(iPair)%iPart_p1, Xi_rel, ProbRot1, ProbRotMax1)
  ELSE
    ProbVib1 = 0.
    ProbRot1 = 0.
  END IF
  IF ((SpecDSMC(iSpec2)%InterID.EQ.2).OR.(SpecDSMC(iSpec2)%InterID.EQ.20)) THEN
    CALL DSMC_calc_P_vib(iSpec2, Xi_rel, iElem, ProbVib2)
    CALL DSMC_calc_P_rot(iSpec2, iPair, iPart2, Xi_rel, ProbRot2, ProbRotMax2)
  ELSE
    ProbVib2 = 0.
    ProbRot2 = 0.
  END IF

  ! Calculate probability fractions
  IF(SpecDSMC(iSpec1)%PolyatomicMol) THEN
    IF(DSMC%PolySingleMode) THEN
      ! If single-mode relaxation is considered, every mode has its own probability while it accumulates all the previous ones
      ! Here, the last mode of the molecule has the highest probability, later it is found which exact mode is going to be relaxed
      iPolyatMole = SpecDSMC(iSpec1)%SpecToPolyArray
      ProbFrac1 = PolyatomMolDSMC(iPolyatMole)%VibRelaxProb(PolyatomMolDSMC(iPolyatMole)%VibDOF)
    ELSE
      ProbFrac1 = ProbVib1
    END IF
  ELSE
    ProbFrac1 = ProbVib1
  END IF
  IF(SpecDSMC(iSpec2)%PolyatomicMol) THEN
    IF(DSMC%PolySingleMode) THEN
      iPolyatMole = SpecDSMC(iSpec2)%SpecToPolyArray
      ProbFrac2   = ProbFrac1 + PolyatomMolDSMC(iPolyatMole)%VibRelaxProb(PolyatomMolDSMC(iPolyatMole)%VibDOF)
    ELSE
      ProbFrac2   = ProbFrac1 + ProbVib2
    END IF
  ELSE
    ProbFrac2     = ProbFrac1 + ProbVib2
  END IF
  ProbFrac3       = ProbFrac2 + ProbRot1
  ProbFrac4       = ProbFrac3 + ProbRot2

  ! Check if sum of probabilities is less than 1.
  IF (ProbFrac4.GT. 1.0) THEN
    CALL Abort(&
__STAMP__&
,'Error! Sum of internal relaxation probabilities > 1.0 for iPair ',iPair)
  END IF

  ! Select relaxation procedure (vibration, rotation)
  CALL RANDOM_NUMBER(iRan)
  IF(iRan .LT. ProbFrac1) THEN                    !            R1 < A1
    IF (SpecDSMC(iSpec1)%PolyatomicMol.AND.DSMC%PolySingleMode) THEN
      ! Determination which vibrational mode should be relaxed for single-mode polyatomic relaxation
      iPolyatMole = SpecDSMC(iSpec1)%SpecToPolyArray
      DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
        IF(PolyatomMolDSMC(iPolyatMole)%VibRelaxProb(iDOF).GT.iRan) THEN
          DoVib1 = .TRUE.
          DOFRelax = iDOF
          EXIT
        END IF
      END DO
    ELSE
      DoVib1 = .TRUE.
    END IF
  ELSEIF(iRan .LT. ProbFrac2) THEN                !      A1 <= R1 < A2
    IF (SpecDSMC(iSpec2)%PolyatomicMol.AND.DSMC%PolySingleMode) THEN
      ! Determination which vibrational mode should be relaxed for single-mode polyatomic relaxation
      iPolyatMole = SpecDSMC(iSpec2)%SpecToPolyArray
      DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
        IF(ProbFrac1 + PolyatomMolDSMC(iPolyatMole)%VibRelaxProb(iDOF).GT.iRan) THEN
          DoVib2 = .TRUE.
          DOFRelax = iDOF
          EXIT
        END IF
      END DO
    ELSE
      DoVib2 = .TRUE.
    END IF
  ELSEIF(iRan .LT. ProbFrac3) THEN                !      A2 <= R1 < A3
    DoRot1 = .TRUE.
  ELSEIF(iRan .LT. ProbFrac4) THEN                !      A3 <= R1 < A4
    DoRot2 = .TRUE.
  END IF

!--------------------------------------------------------------------------------------------------!
! Vibrational Relaxation
!--------------------------------------------------------------------------------------------------!

  IF(DoVib1) THEN
    ! check if correction term for BL redistribution (depending on relaxation model) is needed
    BLCorrFact = 1.
    ! Adding the interal energy of the particle to be redistributed
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(iPart1,1)
    IF(SpecDSMC(PartSpecies(iPart1))%PolyatomicMol) THEN
      IF (.NOT.DSMC%PolySingleMode) THEN
        ! --------------------------------------------------------------------------------------------------!
        !  Multi-mode relaxation with the Metropolis-Hastings method
        ! --------------------------------------------------------------------------------------------------!
        CALL DSMC_VibRelaxPoly(iPair,iPart1,FakXi)
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(iPart1,1)
      ELSE
        ! --------------------------------------------------------------------------------------------------!
        !  Single-mode relaxation with loop over all vibrational modes
        ! --------------------------------------------------------------------------------------------------!
        !  Not all vibrational energy is redistributed but only the energy of the selected vibrational degree of freedom
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(iPart1,1)
        iPolyatMole          = SpecDSMC(iSpec1)%SpecToPolyArray
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + (VibQuantsPar(iPart1)%Quants(DOFRelax) + DSMC%GammaQuant) &
                             * BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(DOFRelax)
        MaxColQua = Coll_pData(iPair)%Ec/(BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(DOFRelax))- DSMC%GammaQuant
        iQuaMax = MIN(INT(MaxColQua) + 1, PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(DOFRelax))
        CALL RANDOM_NUMBER(iRan)
        iQua = INT(iRan * iQuaMax)
        CALL RANDOM_NUMBER(iRan)
        DO WHILE (iRan.GT.(1 - iQua/MaxColQua)**FakXi)
          !Laux1996 diss page 31
          CALL RANDOM_NUMBER(iRan)
          iQua = INT(iRan * iQuaMax)
          CALL RANDOM_NUMBER(iRan)
        END DO
        PartStateIntEn(iPart1,1) = PartStateIntEn(iPart1,1) - (VibQuantsPar(iPart1)%Quants(DOFRelax) &
                                      + DSMC%GammaQuant) * BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(DOFRelax) &
                                      + (iQua + DSMC%GammaQuant) * BoltzmannConst &
                                      * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(DOFRelax)
        VibQuantsPar(iPart1)%Quants(DOFRelax) = iQua
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec &
                             - (iQua + DSMC%GammaQuant) * BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(DOFRelax)
      END IF
    ELSE
      CALL DSMC_VibRelaxDiatomic(iPair,iPart1,FakXi)
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(iPart1,1)
    END IF

  END IF

  IF(DoVib2) THEN
    ! check if correction term for BL redistribution (depending on relaxation model) is needed
    BLCorrFact = 1.
    ! Adding the interal energy of the particle to be redistributed (not if single-mode polyatomic relaxation is enabled)
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(iPart2,1)
    IF(SpecDSMC(PartSpecies(iPart2))%PolyatomicMol) THEN
      IF (.NOT.DSMC%PolySingleMode) THEN
        ! --------------------------------------------------------------------------------------------------!
        !  Multi-mode relaxation with the Metropolis-Hastings method
        ! --------------------------------------------------------------------------------------------------!
        CALL DSMC_VibRelaxPoly(iPair,iPart2,FakXi)
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(iPart2,1)
      ELSE
        ! --------------------------------------------------------------------------------------------------!
        !  Single-mode relaxation with loop over all vibrational modes
        ! --------------------------------------------------------------------------------------------------!
        !  Not all vibrational energy is redistributed but only the energy of the selected vibrational degree of freedom
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(iPart2,1)
        iPolyatMole = SpecDSMC(iSpec2)%SpecToPolyArray
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + (VibQuantsPar(iPart2)%Quants(DOFRelax) + DSMC%GammaQuant) &
                             * BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(DOFRelax)
        MaxColQua = Coll_pData(iPair)%Ec/(BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(DOFRelax))  &
                  - DSMC%GammaQuant
        iQuaMax = MIN(INT(MaxColQua) + 1, PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(DOFRelax))
        CALL RANDOM_NUMBER(iRan)
        iQua = INT(iRan * iQuaMax)
        CALL RANDOM_NUMBER(iRan)
        DO WHILE (iRan.GT.(1 - iQua/MaxColQua)**FakXi)
         !laux diss page 31
         CALL RANDOM_NUMBER(iRan)
         iQua = INT(iRan * iQuaMax)
         CALL RANDOM_NUMBER(iRan)
        END DO
        PartStateIntEn(iPart2,1) = PartStateIntEn(iPart2,1) - (VibQuantsPar(iPart2)%Quants(DOFRelax) &
                                      + DSMC%GammaQuant) * BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(DOFRelax) &
                                      + (iQua + DSMC%GammaQuant) * BoltzmannConst &
                                      * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(DOFRelax)
        VibQuantsPar(iPart2)%Quants(DOFRelax) = iQua
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - (iQua + DSMC%GammaQuant) * BoltzmannConst &
                             * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(DOFRelax)
      END IF
    ELSE
      CALL DSMC_VibRelaxDiatomic(iPair,iPart2,FakXi)
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(iPart2,1)
    END IF
  END IF

!--------------------------------------------------------------------------------------------------!
! Rotational Relaxation
!--------------------------------------------------------------------------------------------------!

  IF(DoRot1) THEN
    !check if correction term in distribution (depending on relaxation model) is needed
    IF(DSMC%RotRelaxProb.EQ.3.0) THEN
      BLCorrFact = ProbRot1 / ProbRotMax1
    ELSE
      BLCorrFact = 1.
    END IF
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(iPart1,2)    ! adding ro en to collision energy
    ! check for polyatomic treatment
    IF(SpecDSMC(PartSpecies(iPart1))%PolyatomicMol.AND. &
        (SpecDSMC(PartSpecies(iPart1))%Xi_Rot.EQ.3)) THEN
      CALL DSMC_RotRelaxPoly(iPair, iPart1, FakXi)
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(iPart1,2)
    ! no polyatomic treatment
    ELSE
     CALL RANDOM_NUMBER(iRan)
      PartStateIntEnTemp = iRan * Coll_pData(iPair)%Ec
      CALL RANDOM_NUMBER(iRan)
      DO WHILE(iRan.GT.(1. - PartStateIntEnTemp/Coll_pData(iPair)%Ec)**FakXi*BLCorrFact)      ! FakXi hier nur 0.5*Xi_rel - 1 !
        CALL RANDOM_NUMBER(iRan)
        PartStateIntEnTemp = iRan * Coll_pData(iPair)%Ec
        CALL RANDOM_NUMBER(iRan)
      END DO
      PartStateIntEn(iPart1,2) = PartStateIntEnTemp
      Coll_pData(iPair)%Ec     = Coll_pData(iPair)%Ec - PartStateIntEn(iPart1,2)
    END IF
  END IF

  IF(DoRot2) THEN
    !check if correction term in distribution (depending on relaxation model) is needed
    IF(DSMC%RotRelaxProb.EQ.3.0) THEN
      BLCorrFact = ProbRot2 / ProbRotMax2
    ELSE
      BLCorrFact = 1.
    END IF
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(iPart2,2)    ! adding rot en to collision en
    IF(SpecDSMC(PartSpecies(iPart2))%PolyatomicMol.AND.(SpecDSMC(PartSpecies(iPart2))%Xi_Rot.EQ.3)) THEN
      CALL DSMC_RotRelaxPoly(iPair, iPart2, FakXi)
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(iPart2,2)
    ELSE
      CALL RANDOM_NUMBER(iRan)
      PartStateIntEnTemp = iRan * Coll_pData(iPair)%Ec
      CALL RANDOM_NUMBER(iRan)
      DO WHILE(iRan.GT.(1. - PartStateIntEnTemp/Coll_pData(iPair)%Ec)**FakXi*BLCorrFact)       ! FakXi hier nur 0.5*Xi_rel -1 !
        CALL RANDOM_NUMBER(iRan)
        PartStateIntEnTemp = iRan * Coll_pData(iPair)%Ec
        CALL RANDOM_NUMBER(iRan)
      END DO
      PartStateIntEn(iPart2,2) = PartStateIntEnTemp
      Coll_pData(iPair)%Ec     = Coll_pData(iPair)%Ec - PartStateIntEn(iPart2,2)
    END IF
  END IF

!--------------------------------------------------------------------------------------------------!
! Calculation of new particle velocities
!--------------------------------------------------------------------------------------------------!

  FracMassCent1 = CollInf%FracMassCent(PartSpecies(iPart1), Coll_pData(iPair)%PairType)
  FracMassCent2 = CollInf%FracMassCent(PartSpecies(iPart2), Coll_pData(iPair)%PairType)

  ! Calculate center of mass velocity
  VeloMx = FracMassCent1 * PartState(iPart1, 4) + FracMassCent2 * PartState(iPart2, 4)
  VeloMy = FracMassCent1 * PartState(iPart1, 5) + FracMassCent2 * PartState(iPart2, 5)
  VeloMz = FracMassCent1 * PartState(iPart1, 6) + FracMassCent2 * PartState(iPart2, 6)

  IF (CollInf%alphaVSS(iSpec1,iSpec2).GT.1) THEN
    
    ! Calculate relative velocites and the squared velocities
    cRelaX = PartState(iPart1, 4) - PartState(iPart2, 4)
    cRelaY = PartState(iPart1, 5) - PartState(iPart2, 5)
    cRelaZ = PartState(iPart1, 6) - PartState(iPart2, 6)

    Coll_pData(iPair)%cRela2 = 2. * Coll_pData(iPair)%Ec/CollInf%MassRed(Coll_pData(iPair)%PairType)

    ! Calculation of post collision velocity vector in reference frame and retransformation to center-of-mass frame
    RanVelo(1:3) = DiceDeflectedVelocityVector(Coll_pData(iPair)%cRela2 , CollInf%alphaVSS(iSpec1,iSpec2), cRelaX , cRelaY , cRelaZ)

  ELSE ! alphaVSS .LE. 1
    ! Calculation of post collision velocity vector in reference frame 
    RanVelo(1:3) = DiceDeflectedVelocityVector(Coll_pData(iPair)%cRela2, CollInf%alphaVSS(iSpec1,iSpec2))
  END IF  ! alphaVSS 

!WRITE(*,*) "DiceDeflectedVector in Relax Gimelshein",RanVelo/sqrt(Coll_pData(iPair)%cRela2) !to be solved 

  ! deltaV particle 1 (post collision particle 1 velocity in laboratory frame)
  DSMC_RHS(iPart1,1) = VeloMx + FracMassCent2*RanVelo(1) - PartState(iPart1, 4)
  DSMC_RHS(iPart1,2) = VeloMy + FracMassCent2*RanVelo(2) - PartState(iPart1, 5)
  DSMC_RHS(iPart1,3) = VeloMz + FracMassCent2*RanVelo(3) - PartState(iPart1, 6)
  ! deltaV particle 2 (post collision particle 2 velocity in laboratory frame)
  DSMC_RHS(iPart2,1) = VeloMx - FracMassCent1*RanVelo(1) - PartState(iPart2, 4)
  DSMC_RHS(iPart2,2) = VeloMy - FracMassCent1*RanVelo(2) - PartState(iPart2, 5)
  DSMC_RHS(iPart2,3) = VeloMz - FracMassCent1*RanVelo(3) - PartState(iPart2, 6)

END SUBROUTINE DSMC_Relax_Col_Gimelshein


SUBROUTINE ReactionDecision(iPair, RelaxToDo, iElem, NodeVolume, NodePartNum)
!===================================================================================================================================
! Decision of reaction type (recombination, exchange, dissociation, CEX/MEX and multiple combinations of those)
!===================================================================================================================================
! MODULES
  USE MOD_Globals,                ONLY : Abort
  USE MOD_Globals_Vars,           ONLY : BoltzmannConst, ElementaryCharge
  USE MOD_DSMC_Vars,              ONLY : Coll_pData, CollInf, DSMC, SpecDSMC, PartStateIntEn, ChemReac, RadialWeighting
  USE MOD_Particle_Vars,          ONLY : Species, PartSpecies, PEM, VarTimeStep
  USE MOD_DSMC_ChemReact,         ONLY : DSMC_Chemistry, simpleCEX, simpleMEX, CalcReactionProb
  USE MOD_Globals,                ONLY : Unit_stdOut
  USE MOD_Particle_Mesh_Vars,     ONLY : GEO
  USE MOD_DSMC_QK_PROCEDURES,     ONLY : QK_dissociation, QK_recombination, QK_exchange, QK_ImpactIonization, QK_IonRecombination
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iPair
  INTEGER, INTENT(IN)           :: iElem
  LOGICAL, INTENT(INOUT)        :: RelaxToDo
  REAL, INTENT(IN), OPTIONAL    :: NodeVolume
  INTEGER, INTENT(IN), OPTIONAL :: NodePartNum
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: CaseOfReaction, iReac, PartToExec, PartReac2, iPart3, iQuaMax
INTEGER                       :: PartToExecSec, PartReac2Sec, iReac2, iReac3, iReac4, ReacToDo
INTEGER                       :: nPartNode, PairForRec, nPair
REAL                          :: EZeroPoint, Volume, sigmaCEX, sigmaMEX, IonizationEnergy, NumDens
REAL (KIND=8)                 :: ReactionProb, ReactionProb2, ReactionProb3, ReactionProb4
REAL (KIND=8)                 :: iRan, iRan2, iRan3
INTEGER                       :: iPart1, iPart2                         ! Colliding particles 1 and 2
!===================================================================================================================================
 iPart1 = Coll_pData(iPair)%iPart_p1
 iPart2 = Coll_pData(iPair)%iPart_p2

 IF(ChemReac%NumOfReact.EQ.0) THEN
    CaseOfReaction = 0
  ELSE
    CaseOfReaction = ChemReac%ReactCase(PartSpecies(iPart1),PartSpecies(iPart2))
  END IF
  IF (PRESENT(NodeVolume)) THEN
    Volume = NodeVolume
  ELSE
    Volume = GEO%Volume(iElem)
  END IF
  IF (PRESENT(NodePartNum)) THEN
    nPartNode = NodePartNum
  ELSE
    nPartNode = PEM%pNumber(iElem)
  END IF
  nPair = INT(nPartNode/2)
  IF(RadialWeighting%DoRadialWeighting) THEN
    NumDens = SUM(CollInf%Coll_SpecPartNum(:)) / Volume
  ELSE IF (VarTimeStep%UseVariableTimeStep) THEN
    NumDens = SUM(CollInf%Coll_SpecPartNum(:)) / Volume * Species(1)%MacroParticleFactor
  ELSE
    NumDens = nPartNode / Volume * Species(1)%MacroParticleFactor
  END IF
  SELECT CASE(CaseOfReaction)
! ############################################################################################################################### !
    CASE(1) ! Only recombination is possible
! ############################################################################################################################### !
      ! searching third collison partner
      IF(ChemReac%RecombParticle.EQ. 0) THEN
        IF(iPair.LT.(nPair - ChemReac%nPairForRec)) THEN
          PairForRec = nPair - ChemReac%nPairForRec
          iPart3 = Coll_pData(PairForRec)%iPart_p1
        ELSE
          iPart3 = 0
        END IF
      ELSE
        iPart3 = ChemReac%RecombParticle
      END IF
!-----------------------------------------------------------------------------------------------------------------------------------
      IF (iPart3 .GT. 0) THEN
        iReac = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), PartSpecies(iPart3))
        IF(iReac.EQ.0) THEN
          iReac = ChemReac%ReactNumRecomb(PartSpecies(iPart1), PartSpecies(iPart2), PartSpecies(iPart3))
        END IF
        IF (ChemReac%QKProcedure(iReac)) THEN
          ! QK - model
          CALL QK_recombination(iPair,iReac,iPart3,RelaxToDo,NodeVolume,NodePartNum)
        ELSE
!-----------------------------------------------------------------------------------------------------------------------------------
        ! traditional Recombination
          ! Calculation of reaction probability
          CALL CalcReactionProb(iPair,iReac,ReactionProb,iPart3,NumDens)
#if (PP_TimeDiscMethod==42)
          IF (.NOT.DSMC%ReservoirRateStatistic) THEN
            ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + ReactionProb  ! for calculation of reaction rate coefficient
            ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
          END IF
#endif
          CALL RANDOM_NUMBER(iRan)
          IF (ReactionProb.GT.iRan) THEN
#if (PP_TimeDiscMethod==42)
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
#endif
              CALL DSMC_Chemistry(iPair, iReac, iPart3)
              IF(ChemReac%RecombParticle.EQ. 0) THEN
                Coll_pData(PairForRec)%NeedForRec = .TRUE.
                ChemReac%RecombParticle = Coll_pData(PairForRec)%iPart_p2
                ChemReac%nPairForRec    = ChemReac%nPairForRec + 1
              ELSE
                ChemReac%RecombParticle = 0
              END IF
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reaction rate coefficient
            END IF
#endif
            RelaxToDo = .FALSE.
          END IF ! ReactionProb > iRan
        END IF ! Q-K
      END IF ! iPart3 > 0
! ############################################################################################################################### !
    CASE(2) ! Only one dissociation is possible
! ############################################################################################################################### !
      iReac = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 1)
      IF (ChemReac%QKProcedure(iReac)) THEN
        CALL QK_dissociation(iPair,iReac,RelaxToDo)
      ELSE
        ! Arrhenius-based reaction probability
        CALL CalcReactionProb(iPair,iReac,ReactionProb)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac)   = ChemReac%NumReac(iReac)   + ReactionProb  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
        END IF
#endif
        CALL RANDOM_NUMBER(iRan)
        IF (ReactionProb.GT.iRan) THEN
#if (PP_TimeDiscMethod==42)
          ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
          IF (.NOT.DSMC%ReservoirSimuRate) THEN
#endif
            CALL DSMC_Chemistry(iPair, iReac)
#if (PP_TimeDiscMethod==42)
          END IF
          IF (DSMC%ReservoirRateStatistic) THEN
            ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reaction rate coefficient
          END IF
#endif
          RelaxToDo = .FALSE.
        END IF
      END IF
! ############################################################################################################################### !
    CASE(3) ! Only one exchange reaction is possible
! ############################################################################################################################### !
      iReac = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 1)
      IF (ChemReac%QKProcedure(iReac)) THEN
        CALL QK_exchange(iPair,iReac,RelaxToDo)
      ELSE
        ! Arrhenius-based reaction probability
        CALL CalcReactionProb(iPair,iReac,ReactionProb)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac)   = ChemReac%NumReac(iReac)   + ReactionProb  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
        END IF
#endif
        CALL RANDOM_NUMBER(iRan)
        IF (ReactionProb.GT.iRan) THEN
#if (PP_TimeDiscMethod==42)
          ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
          IF (.NOT.DSMC%ReservoirSimuRate) THEN
#endif
            CALL DSMC_Chemistry(iPair, iReac)
#if (PP_TimeDiscMethod==42)
          END IF
          IF (DSMC%ReservoirRateStatistic) THEN
            ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reaction rate coefficient
          END IF
#endif
          RelaxToDo = .FALSE.
        END IF
      END IF
! ############################################################################################################################### !
    CASE(4) ! One dissociation and one exchange reaction are possible
! ############################################################################################################################### !
      iReac  = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 1)
      iReac2 = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 2)
      IF (ChemReac%QKProcedure(iReac)) THEN
        ! first check, if the the molecule dissociate, afterwards, check if an exchange reaction is possible
        CALL QK_dissociation(iPair,iReac,RelaxToDo)
        IF (RelaxToDo) THEN
        ! exchange reactions
          iReac = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 2)
          IF (ChemReac%QKProcedure(iReac)) THEN
            CALL QK_exchange(iPair,iReac,RelaxToDo)
          ELSE
            ! Arrhenius based Exchange Reaction
            Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%cRela2                  &
                                 + PartStateIntEn(iPart1,1) + PartStateIntEn(iPart2,1) &
                                 + PartStateIntEn(iPart1,2) + PartStateIntEn(iPart2,2)
            CALL CalcReactionProb(iPair,iReac,ReactionProb)
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac)   = ChemReac%NumReac(iReac)   + ReactionProb  ! for calculation of reaction rate coefficient
              ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
            END IF
#endif
            CALL RANDOM_NUMBER(iRan)
            IF (ReactionProb.GT.iRan) THEN
#if (PP_TimeDiscMethod==42)
              ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
              IF (.NOT.DSMC%ReservoirSimuRate) THEN
#endif
                CALL DSMC_Chemistry(iPair, iReac)
#if (PP_TimeDiscMethod==42)
              END IF
              IF (DSMC%ReservoirRateStatistic) THEN
                ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reaction rate coefficient
              END IF
#endif
              RelaxToDo = .FALSE.
            END IF
          END IF
        END IF
      ELSE
!-----------------------------------------------------------------------------------------------------------------------------------
        ! Arrhenius-based reaction probability
        ! calculation of dissociation probability
        CALL CalcReactionProb(iPair,iReac,ReactionProb)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac)   = ChemReac%NumReac(iReac)   + ReactionProb  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
        END IF
#endif
        ! calculation of exchange reaction probability
        CALL CalcReactionProb(iPair,iReac2,ReactionProb2)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac2)   = ChemReac%NumReac(iReac2)   + ReactionProb2  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac2) = ChemReac%ReacCount(iReac2) + 1
        END IF
#endif
        CALL RANDOM_NUMBER(iRan)
        IF ((ReactionProb + ReactionProb2).GT.iRan) THEN
          CALL RANDOM_NUMBER(iRan)
          IF((ReactionProb/(ReactionProb + ReactionProb2)).GT.iRan) THEN
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
#endif
              CALL DSMC_Chemistry(iPair, iReac)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reaction rate coefficient
            END IF
#endif
          ELSE
#if (PP_TimeDiscMethod==42)
          ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
#endif
              CALL DSMC_Chemistry(iPair, iReac2)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac2) = ChemReac%NumReac(iReac2) + 1  ! for calculation of reaction rate coefficient
            END IF
#endif
          END IF
          RelaxToDo = .FALSE.
        END IF
      END IF
! ############################################################################################################################### !
    CASE(5) ! Two dissociation reactions are possible
! ############################################################################################################################### !
      iReac  = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 1)
      iReac2 = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 2)
      IF (ChemReac%QKProcedure(iReac).AND.ChemReac%QKProcedure(iReac2)) THEN ! both Reaction QK
        ! collision energy without internal energy
        Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%cRela2
        ! first pseudo reaction probability
        IF (ChemReac%DefinedReact(iReac,1,1).EQ.PartSpecies(iPart1)) THEN
          PartToExec = iPart1
          PartReac2  = iPart2
        ELSE
          PartToExec = iPart2
          PartReac2  = iPart1
        END IF
        ReactionProb = ( Coll_pData(iPair)%Ec + PartStateIntEn(PartToExec,1) - SpecDSMC(PartSpecies(PartToExec))%Ediss_eV &
                     * ElementaryCharge  ) / ( Coll_pData(iPair)%Ec + PartStateIntEn(PartToExec,1) )
        IF (ReactionProb.LE.0.) THEN
          ReactionProb = 0.
        END IF
        IF (ChemReac%DefinedReact(iReac2,1,1).EQ.PartSpecies(iPart1)) THEN
          PartToExecSec = iPart1
          PartReac2Sec  = iPart2
        ELSE
          PartToExecSec = iPart2
          PartReac2Sec  = iPart1
        END IF
        ! pseudo probability for second reaction
        ReactionProb2 = (Coll_pData(iPair)%Ec + PartStateIntEn(PartToExecSec,1 ) - SpecDSMC(PartSpecies(PartToExecSec))%Ediss_eV &
                      * ElementaryCharge ) / ( Coll_pData(iPair)%Ec + PartStateIntEn(PartToExecSec,1) )
        IF (ReactionProb2.LE.0.) THEN
          ReactionProb2 = 0
        END IF
        IF (ReactionProb.GT.0.) THEN
          ! determine if first molecule dissociates
          CALL RANDOM_NUMBER(iRan)
          IF ( ReactionProb / ( ReactionProb + ReactionProb2) .gt. iRan) THEN
            CALL QK_dissociation(iPair,iReac,RelaxToDo)
            ! first molecule does not dissociate. what is the second particle doing?
          ELSE IF ( ReactionProb2 .gt. 0 ) THEN
            ! dissociation second molecule
            CALL QK_dissociation(iPair,iReac2,RelaxToDo)
          END IF
        ! ReactionProb = 0, check if second dissociation is possible
        ELSE IF ( ReactionProb2 .gt. 0 ) THEN
!           ! dissociationof second molecule
            CALL QK_dissociation(iPair,iReac2,RelaxToDo)
        END IF
      ELSE IF ( ChemReac%QKProcedure(iReac) .OR. ChemReac%QKProcedure(iReac2) ) THEN ! only one Reaction QK
        ! collision energy without internal energy
        Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%cRela2
        ! set QK reaction as first tested reaction
        IF (ChemReac%QKProcedure(iReac2)) THEN
          iReac  = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 2)
          iReac2 = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 1)
        END IF
        ! first pseude reaction probability
        IF (ChemReac%DefinedReact(iReac,1,1).EQ.PartSpecies(iPart1)) THEN
          PartToExec = iPart1
          PartReac2  = iPart2
        ELSE
          PartToExec = iPart2
          PartReac2  = iPart1
        END IF
        ReactionProb = ( Coll_pData(iPair)%Ec + PartStateIntEn(PartToExec,1) - SpecDSMC(PartSpecies(PartToExec))%Ediss_eV &
                     * ElementaryCharge  ) / ( Coll_pData(iPair)%Ec + PartStateIntEn(PartToExec,1) )
        IF (ReactionProb.LE.0.) THEN
          ReactionProb = 0.
        END IF
        ! second reaction probability
        IF (ChemReac%DefinedReact(iReac2,1,1).EQ.PartSpecies(iPart1)) THEN
          PartToExecSec = iPart1
          PartReac2Sec  = iPart2
        ELSE
          PartToExecSec = iPart2
          PartReac2Sec  = iPart1
        END IF
        ! pseudo probability for second reaction
        ReactionProb2 = (Coll_pData(iPair)%Ec + PartStateIntEn(PartToExecSec,1 ) - SpecDSMC(PartSpecies(PartToExecSec))%Ediss_eV &
                      * ElementaryCharge ) / ( Coll_pData(iPair)%Ec + PartStateIntEn(PartToExecSec,1) )
        IF (ReactionProb2.LE.0.) THEN
          ReactionProb2 = 0.
        END IF
        IF ( ReactionProb .gt. 0 ) THEN
          ! determine if first molecule dissociates
          CALL RANDOM_NUMBER(iRan)
          IF ( ReactionProb / ( ReactionProb + ReactionProb2) .gt. iRan) THEN
            CALL QK_dissociation(iPair,iReac,RelaxToDo)
          ELSE
          ! perform Arrhenius dissociation
            Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%cRela2 &
                                 + PartStateIntEn(iPart1,1) + PartStateIntEn(iPart2,1) &
                                 + PartStateIntEn(iPart1,2) + PartStateIntEn(iPart2,2)
            EZeroPoint           = DSMC%GammaQuant*BoltzmannConst*SpecDSMC(PartSpecies(PartReac2Sec))%CharaTVib
            IF((Coll_pData(iPair)%Ec-EZeroPoint).GE.ChemReac%EActiv(iReac2)) THEN
              ReactionProb = ChemReac%ReactInfo(iReac2)%Beta_Diss_Arrhenius(                                                     &
                             ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExecSec))                                            &
                             , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2Sec)))                                          &
                             * (Coll_pData(iPair)%Ec - ChemReac%EActiv(iReac2))                                                  &
                             ** (ChemReac%Arrhenius_Powerfactor(iReac2) - 1.5                                                    &
                             + CollInf%omegaLaux(ChemReac%DefinedReact(iReac2,1,1),ChemReac%DefinedReact(iReac2,1,1))            &
                             + ChemReac%ReactInfo(iReac2)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExecSec))      &
                             , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2Sec)))/2)                                       &
                             * Coll_pData(iPair)%Ec                                                                              &
                             ** (1.0 - SpecDSMC(ChemReac%DefinedReact(iReac2,1,1))%VFD_Phi3_Factor                               &
                             - ChemReac%ReactInfo(iReac2)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExecSec))      &
                             , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2Sec)))/2)                                       &
                             * PartStateIntEn(PartToExecSec,1) ** SpecDSMC(ChemReac%DefinedReact(iReac2,1,1))%VFD_Phi3_Factor
            ELSE
              ReactionProb = 0.0
            END IF
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + ReactionProb  ! for calculation of reaction rate coefficient
              ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
            END IF
#endif
            CALL RANDOM_NUMBER(iRan)
            IF (ReactionProb.GT.iRan) THEN
#if (PP_TimeDiscMethod==42)
              ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
              IF (.NOT.DSMC%ReservoirSimuRate) THEN
#endif
                CALL DSMC_Chemistry(iPair, iReac2)
#if (PP_TimeDiscMethod==42)
              END IF
              IF (DSMC%ReservoirRateStatistic) THEN
                ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reaction rate coefficient
              END IF
#endif
              RelaxToDo = .FALSE.
            END IF
          END IF
        END IF
      ELSE ! both reactions Arrhenius
!-----------------------------------------------------------------------------------------------------------------------------------
        ! Arrhenius-based reaction probability
        CALL CalcReactionProb(iPair,iReac,ReactionProb)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac)   = ChemReac%NumReac(iReac)   + ReactionProb  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
        END IF
#endif
        ! calculation of dissociation probability
        CALL CalcReactionProb(iPair,iReac2,ReactionProb2)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac2)   = ChemReac%NumReac(iReac2)   + ReactionProb2  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac2) = ChemReac%ReacCount(iReac2) + 1
        END IF
#endif
        CALL RANDOM_NUMBER(iRan)
        IF ((ReactionProb + ReactionProb2).GT.iRan) THEN
          CALL RANDOM_NUMBER(iRan)
          IF((ReactionProb/(ReactionProb + ReactionProb2)).GT.iRan) THEN
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
#endif
              CALL DSMC_Chemistry(iPair, iReac)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reaction rate coefficient
            END IF
#endif
          ELSE
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
#endif
              CALL DSMC_Chemistry(iPair, iReac2)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac2) = ChemReac%NumReac(iReac2) + 1  ! for calculation of reaction rate coefficient
            END IF
#endif
          END IF
          RelaxToDo = .FALSE.
        END IF
      END IF
! ############################################################################################################################### !
    CASE(6) ! ionization or ion recombination
! ############################################################################################################################### !
      ! searching third collison partner
      IF(ChemReac%RecombParticle.EQ. 0) THEN
        IF(iPair.LT.(nPair - ChemReac%nPairForRec)) THEN
          PairForRec = nPair - ChemReac%nPairForRec
          iPart3   = Coll_pData(PairForRec)%iPart_p1
        ELSE
          iPart3 = 0
        END IF
      ELSE
        iPart3 = ChemReac%RecombParticle
      END IF
!--------------------------------------------------------------------------------------------------!
      iReac  = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 1)
!--------------------------------------------------------------------------------------------------!
#if (PP_TimeDiscMethod==42)
      IF (.NOT.DSMC%ReservoirRateStatistic) THEN
        ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
      END IF
#endif
      ! calculation of recombination probability
      IF (iPart3 .GT. 0) THEN
        iReac2 = ChemReac%ReactNumRecomb(PartSpecies(iPart1), PartSpecies(iPart2), PartSpecies(iPart3))
        CALL CalcReactionProb(iPair,iReac2,ReactionProb2,iPart3,NumDens)
      ELSE
        ReactionProb2 = 0.0
      END IF
      CALL RANDOM_NUMBER(iRan)
      IF(ReactionProb2.GT.iRan) THEN
#if (PP_TimeDiscMethod==42)
      IF (.NOT.DSMC%ReservoirRateStatistic) THEN
        ChemReac%NumReac(iReac2) = ChemReac%NumReac(iReac2) + ReactionProb2  ! for calculation of reaction rate coefficient
        ChemReac%ReacCount(iReac2) = ChemReac%ReacCount(iReac2) + 1
      END IF
#endif
          ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
          CALL DSMC_Chemistry(iPair, iReac2, iPart3)
          IF(ChemReac%RecombParticle.EQ. 0) THEN
            Coll_pData(PairForRec)%NeedForRec = .TRUE.
            ChemReac%RecombParticle = Coll_pData(PairForRec)%iPart_p2
            ChemReac%nPairForRec = ChemReac%nPairForRec + 1
          ELSE
            ChemReac%RecombParticle = 0
          END IF
#if (PP_TimeDiscMethod==42)
        END IF
        IF (DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac2) = ChemReac%NumReac(iReac2) + 1  ! for calculation of reaction rate coefficient
        END IF
        RelaxToDo = .FALSE.
# endif
      ELSE
        ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
          CALL QK_ImpactIonization(iPair,iReac,RelaxToDo)
#if (PP_TimeDiscMethod==42)
        END IF
        IF (DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reaction rate coefficient
        END IF
# endif
      END IF

! ############################################################################################################################### !
    CASE(7) ! three diss reactions possible (at least one molecule is polyatomic)
! ############################################################################################################################### !
      iReac  = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 1)
      iReac2 = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 2)
      iReac3 = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 3)
      IF ( ChemReac%QKProcedure(iReac) .OR. ChemReac%QKProcedure(iReac2) .OR. ChemReac%QKProcedure(iReac3) ) THEN ! all Q-K
          CALL Abort(&
__STAMP__&
,'ERROR! Reaction case not supported with Q-K reactions!')
!--------------------------------------------------------------------------------------------------!
      ELSE ! all reactions Arrhenius
!--------------------------------------------------------------------------------------------------!
        ! calculation of first dissociation probability
        CALL CalcReactionProb(iPair,iReac,ReactionProb)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + ReactionProb  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
        END IF
#endif
        ! calculation of second dissociation probability
        CALL CalcReactionProb(iPair,iReac2,ReactionProb2)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac2) = ChemReac%NumReac(iReac2) + ReactionProb2  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac2) = ChemReac%ReacCount(iReac2) + 1
        END IF
#endif
        ! calculation of third dissociation probability
        CALL CalcReactionProb(iPair,iReac3,ReactionProb3)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac3) = ChemReac%NumReac(iReac3) + ReactionProb3  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac3) = ChemReac%ReacCount(iReac3) + 1
        END IF
#endif
        CALL RANDOM_NUMBER(iRan)
        IF ((ReactionProb + ReactionProb2 + ReactionProb3).GT.iRan) THEN
          CALL RANDOM_NUMBER(iRan)
          CALL RANDOM_NUMBER(iRan2)
          IF((ReactionProb/(ReactionProb + ReactionProb2 + ReactionProb3)).GT.iRan) THEN
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac)
#if (PP_TimeDiscMethod==42)
            END IF
            IF ( DSMC%ReservoirRateStatistic  ) THEN
              ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          ELSEIF(ReactionProb2/(ReactionProb2 + ReactionProb3).GT.iRan2) THEN
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac2)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac2) = ChemReac%NumReac(iReac2) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          ELSE
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac3)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac3) = ChemReac%NumReac(iReac3) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          END IF
          RelaxToDo = .FALSE.
        END IF
      END IF
! ############################################################################################################################### !
    CASE(8) ! four diss reactions possible (at least one polyatomic molecule)
! ############################################################################################################################### !
      iReac  = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 1)
      iReac2 = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 2)
      iReac3 = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 3)
      iReac4 = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 4)
      IF ( ChemReac%QKProcedure(iReac) .OR. ChemReac%QKProcedure(iReac2) &
      .OR. ChemReac%QKProcedure(iReac3) .OR. ChemReac%QKProcedure(iReac4)) THEN ! all Q-K
          CALL Abort(&
__STAMP__&
,'ERROR! Reaction case not supported with Q-K reactions!')
!--------------------------------------------------------------------------------------------------!
      ELSE ! all reactions Arrhenius
!--------------------------------------------------------------------------------------------------!
        ! calculation of first dissociation probability
        CALL CalcReactionProb(iPair,iReac,ReactionProb)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac)   = ChemReac%NumReac(iReac)   + ReactionProb  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
        END IF
#endif
        ! calculation of second dissociation probability
        CALL CalcReactionProb(iPair,iReac2,ReactionProb2)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac2)   = ChemReac%NumReac(iReac2)   + ReactionProb2  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac2) = ChemReac%ReacCount(iReac2) + 1
        END IF
#endif
        ! calculation of third dissociation probability
        CALL CalcReactionProb(iPair,iReac3,ReactionProb3)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac3)   = ChemReac%NumReac(iReac3)   + ReactionProb3  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac3) = ChemReac%ReacCount(iReac3) + 1
        END IF
#endif
        ! calculation of fourth dissociation probability
        CALL CalcReactionProb(iPair,iReac4,ReactionProb4)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac4)   = ChemReac%NumReac(iReac4)   + ReactionProb4  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac4) = ChemReac%ReacCount(iReac4) + 1
        END IF
#endif
        CALL RANDOM_NUMBER(iRan)
        IF ((ReactionProb + ReactionProb2 + ReactionProb3 + ReactionProb4).GT.iRan) THEN
          CALL RANDOM_NUMBER(iRan)
          CALL RANDOM_NUMBER(iRan2)
          CALL RANDOM_NUMBER(iRan3)
          IF((ReactionProb/(ReactionProb + ReactionProb2 + ReactionProb3 + ReactionProb4)).GT.iRan) THEN
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          ELSEIF(ReactionProb2/(ReactionProb2 + ReactionProb3 + ReactionProb4).GT.iRan2) THEN
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac2)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac2) = ChemReac%NumReac(iReac2) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          ELSEIF(ReactionProb3/(ReactionProb3 + ReactionProb4).GT.iRan3) THEN
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac3)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac3) = ChemReac%NumReac(iReac3) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          ELSE
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac4)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac4) = ChemReac%NumReac(iReac4) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          END IF
          RelaxToDo = .FALSE.
        END IF
      END IF
! ############################################################################################################################### !
    CASE(9) ! three diss and one exchange reaction possible (at least one polyatomic molecule)
! ############################################################################################################################### !
      iReac  = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 1)
      iReac2 = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 2)
      iReac3 = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 3)
      iReac4 = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 4)
      IF ( ChemReac%QKProcedure(iReac) .OR. ChemReac%QKProcedure(iReac2) &
      .OR. ChemReac%QKProcedure(iReac3) .OR. ChemReac%QKProcedure(iReac4)) THEN ! all Q-K
          CALL Abort(&
__STAMP__&
,'ERROR! Reaction case not supported with Q-K reactions!')
!--------------------------------------------------------------------------------------------------!
      ELSE ! all reactions Arrhenius
!--------------------------------------------------------------------------------------------------!
        ! calculation of first dissociation probability
        CALL CalcReactionProb(iPair,iReac,ReactionProb)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac)   = ChemReac%NumReac(iReac)   + ReactionProb  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
        END IF
#endif
        ! calculation of second dissociation probability
        CALL CalcReactionProb(iPair,iReac2,ReactionProb2)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac2)   = ChemReac%NumReac(iReac2)   + ReactionProb2  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac2) = ChemReac%ReacCount(iReac2) + 1
        END IF
#endif
        ! calculation of third dissociation probability
        CALL CalcReactionProb(iPair,iReac3,ReactionProb3)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac3)   = ChemReac%NumReac(iReac3)   + ReactionProb3  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac3) = ChemReac%ReacCount(iReac3) + 1
        END IF
#endif
        ! calculation of exchange probability
        CALL CalcReactionProb(iPair,iReac4,ReactionProb4)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac4)   = ChemReac%NumReac(iReac4)   + ReactionProb4  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac4) = ChemReac%ReacCount(iReac4) + 1
        END IF
#endif
        CALL RANDOM_NUMBER(iRan)
        IF ((ReactionProb + ReactionProb2 + ReactionProb3 + ReactionProb4).GT.iRan) THEN
          CALL RANDOM_NUMBER(iRan)
          CALL RANDOM_NUMBER(iRan2)
          CALL RANDOM_NUMBER(iRan3)
          IF((ReactionProb/(ReactionProb + ReactionProb2 + ReactionProb3 + ReactionProb4)).GT.iRan) THEN
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          ELSEIF(ReactionProb2/(ReactionProb2 + ReactionProb3 + ReactionProb4).GT.iRan2) THEN
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac2)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac2) = ChemReac%NumReac(iReac2) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          ELSEIF(ReactionProb3/(ReactionProb3 + ReactionProb4).GT.iRan3) THEN
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac3)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac3) = ChemReac%NumReac(iReac3) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          ELSE
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac4)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac4) = ChemReac%NumReac(iReac4) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          END IF
          RelaxToDo = .FALSE.
        END IF
      END IF
! ############################################################################################################################### !
    CASE(10) ! two diss and one exchange reaction possible
! ############################################################################################################################### !
      iReac  = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 1)
      iReac2 = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 2)
      iReac3 = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 3)
      IF ( ChemReac%QKProcedure(iReac) .OR. ChemReac%QKProcedure(iReac2) .OR. ChemReac%QKProcedure(iReac3) ) THEN ! all Q-K
          CALL Abort(&
__STAMP__&
,'ERROR! Reaction case not supported with Q-K reactions!')
!--------------------------------------------------------------------------------------------------!
      ELSE ! all reactions Arrhenius
!--------------------------------------------------------------------------------------------------!
        ! calculation of first dissociation probability
          CALL CalcReactionProb(iPair,iReac,ReactionProb)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac)   = ChemReac%NumReac(iReac)   + ReactionProb  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
        END IF
#endif
        ! calculation of second dissociation probability
          CALL CalcReactionProb(iPair,iReac2,ReactionProb2)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac2)   = ChemReac%NumReac(iReac2)   + ReactionProb2  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac2) = ChemReac%ReacCount(iReac2) + 1
        END IF
#endif
        ! calculation of exchange probability
        CALL CalcReactionProb(iPair,iReac3,ReactionProb3)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac3)   = ChemReac%NumReac(iReac3)   + ReactionProb3  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac3) = ChemReac%ReacCount(iReac3) + 1
        END IF
#endif
        CALL RANDOM_NUMBER(iRan)
        IF ((ReactionProb + ReactionProb2 + ReactionProb3).GT.iRan) THEN
          CALL RANDOM_NUMBER(iRan)
          CALL RANDOM_NUMBER(iRan2)
          IF((ReactionProb/(ReactionProb + ReactionProb2 + ReactionProb3)).GT.iRan) THEN
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          ELSEIF(ReactionProb2/(ReactionProb2 + ReactionProb3).GT.iRan2) THEN
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac2)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac2) = ChemReac%NumReac(iReac2) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          ELSE
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac3)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac3) = ChemReac%NumReac(iReac3) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          END IF
          RelaxToDo = .FALSE.
        END IF
      END IF
! ############################################################################################################################### !
    CASE(11) ! two diss, one exchange and one recombination reaction possible
! ############################################################################################################################### !
      ! searching third collison partner
      IF(ChemReac%RecombParticle.EQ. 0) THEN
        IF(iPair.LT.(nPair - ChemReac%nPairForRec)) THEN
          PairForRec = nPair - ChemReac%nPairForRec
          iPart3   = Coll_pData(PairForRec)%iPart_p1
        ELSE
          iPart3 = 0
        END IF
      ELSE
        iPart3 = ChemReac%RecombParticle
      END IF
!--------------------------------------------------------------------------------------------------!
      iReac  = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 1)
      iReac2 = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 2)
      iReac3 = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 3)
      IF ( ChemReac%QKProcedure(iReac) .OR. ChemReac%QKProcedure(iReac2) &
      .OR. ChemReac%QKProcedure(iReac3)) THEN ! all Q-K
        CALL Abort(&
__STAMP__&
,'ERROR! Reaction case not supported with Q-K reactions!')
!--------------------------------------------------------------------------------------------------!
      ELSE ! all reactions Arrhenius
!--------------------------------------------------------------------------------------------------!
        ! calculation of first dissociation probability
        CALL CalcReactionProb(iPair,iReac,ReactionProb)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac)   = ChemReac%NumReac(iReac)   + ReactionProb  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
        END IF
# endif
        ! calculation of second dissociation probability
          CALL CalcReactionProb(iPair,iReac2,ReactionProb2)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac2)   = ChemReac%NumReac(iReac2)   + ReactionProb2  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac2) = ChemReac%ReacCount(iReac2) + 1
        END IF
# endif
        ! calculation of exchange probability
        CALL CalcReactionProb(iPair,iReac3,ReactionProb3)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac3)   = ChemReac%NumReac(iReac3)   + ReactionProb3  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac3) = ChemReac%ReacCount(iReac3) + 1
        END IF
#endif
        ! calculation of recombination probability
        IF (iPart3 .GT. 0) THEN
          iReac4 = ChemReac%ReactNumRecomb(PartSpecies(iPart1), PartSpecies(iPart2), PartSpecies(iPart3))
          CALL CalcReactionProb(iPair,iReac4,ReactionProb4,iPart3,NumDens)
        ELSE
          ReactionProb4 = 0.0
        END IF
#if (PP_TimeDiscMethod==42)
        IF ((.NOT.DSMC%ReservoirRateStatistic).AND.(iPart3.GT.0)) THEN
          ChemReac%NumReac(iReac4) = ChemReac%NumReac(iReac4) + ReactionProb4  ! for calculation of reaction rate coefficient
        END IF
# endif
        ! reaction decision
        CALL RANDOM_NUMBER(iRan)
        IF ((ReactionProb + ReactionProb2 + ReactionProb3 + ReactionProb4).GT.iRan) THEN
          CALL RANDOM_NUMBER(iRan)
          CALL RANDOM_NUMBER(iRan2)
          CALL RANDOM_NUMBER(iRan3)
          IF((ReactionProb/(ReactionProb + ReactionProb2 + ReactionProb3 + ReactionProb4)).GT.iRan) THEN
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          ELSEIF(ReactionProb2/(ReactionProb2 + ReactionProb3 + ReactionProb4).GT.iRan2) THEN
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac2)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac2) = ChemReac%NumReac(iReac2) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          ELSEIF(ReactionProb3/(ReactionProb3 + ReactionProb4).GT.iRan3) THEN
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac3)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac3) = ChemReac%NumReac(iReac3) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          ELSEIF(ReactionProb4.GT.0.0) THEN ! Probability is set to zero if no third collision partner is found
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac4, iPart3)
              IF(ChemReac%RecombParticle.EQ. 0) THEN
                Coll_pData(PairForRec)%NeedForRec = .TRUE.
                ChemReac%RecombParticle = Coll_pData(PairForRec)%iPart_p2
                ChemReac%nPairForRec = ChemReac%nPairForRec + 1
              ELSE
                ChemReac%RecombParticle = 0
              END IF
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac4) = ChemReac%NumReac(iReac4) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          END IF
          RelaxToDo = .FALSE.
        END IF
      END IF
! ############################################################################################################################### !
    CASE(12) ! two diss and one recomb reaction possible
! ############################################################################################################################### !
      ! searching third collison partner
      IF(ChemReac%RecombParticle.EQ. 0) THEN
        IF(iPair.LT.(nPair - ChemReac%nPairForRec)) THEN
          PairForRec = nPair - ChemReac%nPairForRec
          iPart3   = Coll_pData(PairForRec)%iPart_p1
        ELSE
          iPart3 = 0
        END IF
      ELSE
        iPart3 = ChemReac%RecombParticle
      END IF
!--------------------------------------------------------------------------------------------------!
      iReac  = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 1)
      iReac2 = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 2)
      IF (ChemReac%QKProcedure(iReac).OR.ChemReac%QKProcedure(iReac2)) THEN ! all Q-K
        CALL Abort(&
         __STAMP__,&
        'ERROR! Reaction case not supported with Q-K reactions!')
!--------------------------------------------------------------------------------------------------!
      ELSE ! all reactions Arrhenius
!--------------------------------------------------------------------------------------------------!
        ! calculation of first dissociation probability
        CALL CalcReactionProb(iPair,iReac,ReactionProb)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + ReactionProb  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
        END IF
#endif
        ! calculation of second dissociation probability
          CALL CalcReactionProb(iPair,iReac2,ReactionProb2)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac2) = ChemReac%NumReac(iReac2) + ReactionProb2  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac2) = ChemReac%ReacCount(iReac2) + 1
        END IF
#endif
        ! calculation of recombination probability
        IF (iPart3 .GT. 0) THEN
          iReac3 = ChemReac%ReactNumRecomb(PartSpecies(iPart1), PartSpecies(iPart2), PartSpecies(iPart3))
          CALL CalcReactionProb(iPair,iReac3,ReactionProb3,iPart3,NumDens)
        ELSE
          ReactionProb3 = 0.0
        END IF
#if (PP_TimeDiscMethod==42)
        IF ((.NOT.DSMC%ReservoirRateStatistic).AND.(iPart3.GT.0)) THEN
          ChemReac%NumReac(iReac3) = ChemReac%NumReac(iReac3) + ReactionProb3  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac3) = ChemReac%ReacCount(iReac3) + 1
        END IF
#endif
        CALL RANDOM_NUMBER(iRan)
        IF ((ReactionProb + ReactionProb2 + ReactionProb3).GT.iRan) THEN
          CALL RANDOM_NUMBER(iRan)
          CALL RANDOM_NUMBER(iRan2)
          IF((ReactionProb/(ReactionProb + ReactionProb2 + ReactionProb3)).GT.iRan) THEN
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          ELSEIF(ReactionProb2/(ReactionProb2 + ReactionProb3).GT.iRan2) THEN
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac2)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac2) = ChemReac%NumReac(iReac2) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          ELSEIF(ReactionProb3.GT.0.0) THEN ! Probability is set to zero if no third collision partner is found
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac3, iPart3)
              IF(ChemReac%RecombParticle.EQ. 0) THEN
                Coll_pData(PairForRec)%NeedForRec = .TRUE.
                ChemReac%RecombParticle = Coll_pData(PairForRec)%iPart_p2
                ChemReac%nPairForRec    = ChemReac%nPairForRec + 1
              ELSE
                ChemReac%RecombParticle = 0
              END IF
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac3) = ChemReac%NumReac(iReac3) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          END IF
          RelaxToDo = .FALSE.
        END IF ! Prob > iRan
      END IF ! Q-K
! ############################################################################################################################### !
    CASE(13) ! one diss, one exchange and one recomb reaction possible
! ############################################################################################################################### !
      ! searching third collison partner
      IF(ChemReac%RecombParticle.EQ. 0) THEN
        IF(iPair.LT.(nPair - ChemReac%nPairForRec)) THEN
          PairForRec = nPair - ChemReac%nPairForRec
          iPart3   = Coll_pData(PairForRec)%iPart_p1
        ELSE
          iPart3   = 0
        END IF
      ELSE
        iPart3 = ChemReac%RecombParticle
      END IF
!--------------------------------------------------------------------------------------------------!
      iReac  = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 1)
      iReac2 = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 2)
      IF (ChemReac%QKProcedure(iReac).OR.ChemReac%QKProcedure(iReac2)) THEN ! all Q-K
        CALL Abort(&
         __STAMP__,&
        'ERROR! Reaction case not supported with Q-K reactions!')
!--------------------------------------------------------------------------------------------------!
      ELSE ! all reactions Arrhenius
!--------------------------------------------------------------------------------------------------!
        ! calculation of first dissociation probability
        CALL CalcReactionProb(iPair,iReac,ReactionProb)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac)   = ChemReac%NumReac(iReac)   + ReactionProb  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
        END IF
#endif
        ! calculation of exchange probability
        CALL CalcReactionProb(iPair,iReac2,ReactionProb2)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac2)   = ChemReac%NumReac(iReac2)   + ReactionProb2  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac2) = ChemReac%ReacCount(iReac2) + 1
        END IF
#endif
        ! calculation of recombination probability
        IF (iPart3 .GT. 0) THEN
          iReac3 = ChemReac%ReactNumRecomb(PartSpecies(iPart1), PartSpecies(iPart2), PartSpecies(iPart3))
          CALL CalcReactionProb(iPair,iReac3,ReactionProb3,iPart3,NumDens)
        ELSE
          ReactionProb3 = 0.0
        END IF
#if (PP_TimeDiscMethod==42)
        IF ((.NOT.DSMC%ReservoirRateStatistic).AND.(iPart3.GT.0)) THEN
          ChemReac%NumReac(iReac3) = ChemReac%NumReac(iReac3) + ReactionProb3  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac3) = ChemReac%ReacCount(iReac3) + 1
        END IF
#endif
        CALL RANDOM_NUMBER(iRan)
        IF ((ReactionProb + ReactionProb2 + ReactionProb3).GT.iRan) THEN
          CALL RANDOM_NUMBER(iRan)
          CALL RANDOM_NUMBER(iRan2)
          IF((ReactionProb/(ReactionProb + ReactionProb2 + ReactionProb3)).GT.iRan) THEN
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          ELSEIF(ReactionProb2/(ReactionProb2 + ReactionProb3).GT.iRan2) THEN
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac2)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac2) = ChemReac%NumReac(iReac2) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          ELSEIF(ReactionProb3.GT.0.0) THEN ! Probability is set to zero if no third collision partner is found
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac3, iPart3)
              IF(ChemReac%RecombParticle.EQ. 0) THEN
                Coll_pData(PairForRec)%NeedForRec = .TRUE.
                ChemReac%RecombParticle           = Coll_pData(PairForRec)%iPart_p2
                ChemReac%nPairForRec              = ChemReac%nPairForRec + 1
              ELSE
                ChemReac%RecombParticle           = 0
              END IF
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac3) = ChemReac%NumReac(iReac3) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          END IF
          RelaxToDo = .FALSE.
        END IF ! Prob > iRan
      END IF ! Q-K
! ############################################################################################################################### !
    CASE(14) ! one diss and one recomb reaction possible
! ############################################################################################################################### !
      ! searching third collison partner
      IF(ChemReac%RecombParticle.EQ. 0) THEN
        IF(iPair.LT.(nPair - ChemReac%nPairForRec)) THEN
          PairForRec = nPair - ChemReac%nPairForRec
          iPart3 = Coll_pData(PairForRec)%iPart_p1
        ELSE
          iPart3 = 0
        END IF
      ELSE
        iPart3 = ChemReac%RecombParticle
      END IF
!--------------------------------------------------------------------------------------------------!
      iReac  = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 1)
      IF(ChemReac%QKProcedure(iReac)) THEN ! all Q-K
        CALL Abort(&
         __STAMP__,&
        'ERROR! Reaction case not supported with Q-K reactions!')
!--------------------------------------------------------------------------------------------------!
      ELSE ! all reactions Arrhenius
!--------------------------------------------------------------------------------------------------!
        ! calculation of first dissociation probability
        CALL CalcReactionProb(iPair,iReac,ReactionProb)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac)   = ChemReac%NumReac(iReac)   + ReactionProb  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
        END IF
#endif
        ! calculation of recombination probability
        IF (iPart3 .GT. 0) THEN
          iReac2 = ChemReac%ReactNumRecomb(PartSpecies(iPart1), PartSpecies(iPart2), PartSpecies(iPart3))
          CALL CalcReactionProb(iPair,iReac2,ReactionProb2,iPart3,NumDens)
        ELSE
          ReactionProb2 = 0.0
        END IF
#if (PP_TimeDiscMethod==42)
        IF ((.NOT.DSMC%ReservoirRateStatistic).AND.(iPart3.GT.0)) THEN
          ChemReac%NumReac(iReac2)   = ChemReac%NumReac(iReac2)   + ReactionProb2  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac2) = ChemReac%ReacCount(iReac2) + 1
        END IF
#endif
        CALL RANDOM_NUMBER(iRan)
        IF ((ReactionProb + ReactionProb2).GT.iRan) THEN
          CALL RANDOM_NUMBER(iRan)
          IF((ReactionProb/(ReactionProb + ReactionProb2)).GT.iRan) THEN
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          ELSEIF(ReactionProb2.GT.0.0) THEN ! Probability is set to zero if no third collision partner is found
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac2, iPart3)
              IF(ChemReac%RecombParticle.EQ. 0) THEN
                Coll_pData(PairForRec)%NeedForRec = .TRUE.
                ChemReac%RecombParticle           = Coll_pData(PairForRec)%iPart_p2
                ChemReac%nPairForRec              = ChemReac%nPairForRec + 1
              ELSE
                ChemReac%RecombParticle           = 0
              END IF
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac2) = ChemReac%NumReac(iReac2) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          END IF
          RelaxToDo = .FALSE.
        END IF
      END IF
! ############################################################################################################################### !
    CASE(15) ! one exchange and one recomb reaction possible
! ############################################################################################################################### !
      ! searching third collison partner
      IF(ChemReac%RecombParticle.EQ. 0) THEN
        IF(iPair.LT.(nPair - ChemReac%nPairForRec)) THEN
          PairForRec = nPair - ChemReac%nPairForRec
          iPart3   = Coll_pData(PairForRec)%iPart_p1
        ELSE
          iPart3   = 0
        END IF
      ELSE
        iPart3 = ChemReac%RecombParticle
      END IF
!--------------------------------------------------------------------------------------------------!
      iReac  = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 1)
      IF(ChemReac%QKProcedure(iReac)) THEN ! all Q-K
        CALL Abort(&
         __STAMP__,&
        'ERROR! Reaction case not supported with Q-K reactions!')
!--------------------------------------------------------------------------------------------------!
      ELSE ! all reactions Arrhenius
!--------------------------------------------------------------------------------------------------!
        ! calculation of exchange probability
        CALL CalcReactionProb(iPair,iReac,ReactionProb)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac)   = ChemReac%NumReac(iReac)   + ReactionProb  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
        END IF
#endif
        ! calculation of recombination probability (only if third collision partner is available)
        IF (iPart3 .GT. 0) THEN
          iReac2 = ChemReac%ReactNumRecomb(PartSpecies(iPart1), PartSpecies(iPart2), PartSpecies(iPart3))
          CALL CalcReactionProb(iPair,iReac2,ReactionProb2,iPart3,NumDens)
        ELSE
          ReactionProb2 = 0.0
        END IF
#if (PP_TimeDiscMethod==42)
        IF ((.NOT.DSMC%ReservoirRateStatistic).AND.(iPart3.GT.0)) THEN
          ChemReac%NumReac(iReac2) = ChemReac%NumReac(iReac2) + ReactionProb2  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac2) = ChemReac%ReacCount(iReac2) + 1
        END IF
#endif
        CALL RANDOM_NUMBER(iRan)
        IF ((ReactionProb + ReactionProb2).GT.iRan) THEN
          CALL RANDOM_NUMBER(iRan)
          IF((ReactionProb/(ReactionProb + ReactionProb2)).GT.iRan) THEN
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
#endif
              CALL DSMC_Chemistry(iPair, iReac)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reaction rate coefficient
            END IF
#endif
          ELSEIF(ReactionProb2.GT.0.0) THEN ! Probability is set to zero if no third collision partner is found
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
#endif
              CALL DSMC_Chemistry(iPair, iReac2, iPart3)
              IF(ChemReac%RecombParticle.EQ. 0) THEN
                Coll_pData(PairForRec)%NeedForRec = .TRUE.
                ChemReac%RecombParticle           = Coll_pData(PairForRec)%iPart_p2
                ChemReac%nPairForRec              = ChemReac%nPairForRec + 1
              ELSE
                ChemReac%RecombParticle           = 0
              END IF
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac2) = ChemReac%NumReac(iReac2) + 1  ! for calculation of reaction rate coefficient
            END IF
#endif
          END IF
          RelaxToDo = .FALSE.
        END IF
      END IF

! ############################################################################################################################### !
    CASE(16) ! simple CEX/MEX
! ############################################################################################################################### !
      iReac    = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 1)
      IF (ChemReac%DoScat(iReac)) THEN! MEX
        CALL DSMC_Scat_Col(iPair)
      ELSE
        sigmaCEX = (ChemReac%CEXa(iReac)*0.5*LOG10(Coll_pData(iPair)%cRela2) + ChemReac%CEXb(iReac))
        sigmaMEX = (ChemReac%MEXa(iReac)*0.5*LOG10(Coll_pData(iPair)%cRela2) + ChemReac%MEXb(iReac))
        ReactionProb=0.
        IF ((sigmaMEX.EQ.0.).AND.(sigmaCEX.GT.0.)) THEN
          ReactionProb=1.
        ELSEIF  ((sigmaMEX.GT.0.).AND.(sigmaCEX.GE.0.)) THEN
          ReactionProb=(sigmaCEX/sigmaMEX)/((sigmaCEX/sigmaMEX)+1)
        ELSE
          CALL Abort(&
            __STAMP__&
            ,'ERROR! CEX/MEX cross sections are both zero or at least one of them is negative.')
        END IF
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac)   = ChemReac%NumReac(iReac)   + ReactionProb  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
        END IF
#endif
        CALL RANDOM_NUMBER(iRan)
        IF (ReactionProb.GT.iRan) THEN !CEX, otherwise MEX
#if (PP_TimeDiscMethod==42)
          ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
          IF (.NOT.DSMC%ReservoirSimuRate) THEN
#endif
            CALL simpleCEX(iReac, iPair)
#if (PP_TimeDiscMethod==42)
          END IF
          IF (DSMC%ReservoirRateStatistic) THEN
            ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reaction rate coefficient
          END IF
#endif
        ELSE
          CALL DSMC_Elastic_Col(iPair)
          CALL simpleMEX(iReac, iPair)
        END IF
      END IF !ChemReac%DoScat(iReac)
      RelaxToDo = .FALSE.
! ############################################################################################################################### !
    CASE(17) ! one dissociation, two exchange possible
! ############################################################################################################################### !
      iReac  = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 1)
      iReac2 = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 2)
      iReac3 = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 3)
      IF ( ChemReac%QKProcedure(iReac) .OR. ChemReac%QKProcedure(iReac2) .OR. ChemReac%QKProcedure(iReac3) ) THEN ! all Q-K
          CALL Abort(&
__STAMP__&
,'ERROR! Reaction case not supported with Q-K reactions!')
!--------------------------------------------------------------------------------------------------!
      ELSE ! all reactions Arrhenius
!--------------------------------------------------------------------------------------------------!
        ! calculation of first dissociation probability
        CALL CalcReactionProb(iPair,iReac,ReactionProb)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac)   = ChemReac%NumReac(iReac)   + ReactionProb  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
        END IF
#endif
        ! calculation of second dissociation probability
        CALL CalcReactionProb(iPair,iReac2,ReactionProb2)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac2)   = ChemReac%NumReac(iReac2)   + ReactionProb2  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac2) = ChemReac%ReacCount(iReac2) + 1
        END IF
#endif
        ! calculation of third dissociation probability
        CALL CalcReactionProb(iPair,iReac3,ReactionProb3)
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac3)   = ChemReac%NumReac(iReac3)   + ReactionProb3  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac3) = ChemReac%ReacCount(iReac3) + 1
        END IF
#endif
        CALL RANDOM_NUMBER(iRan)
        IF ((ReactionProb + ReactionProb2 + ReactionProb3).GT.iRan) THEN
          CALL RANDOM_NUMBER(iRan)
          CALL RANDOM_NUMBER(iRan2)
          IF((ReactionProb/(ReactionProb + ReactionProb2 + ReactionProb3)).GT.iRan) THEN
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac)
#if (PP_TimeDiscMethod==42)
            END IF
            IF ( DSMC%ReservoirRateStatistic  ) THEN
              ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          ELSEIF(ReactionProb2/(ReactionProb2 + ReactionProb3).GT.iRan2) THEN
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac2)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac2) = ChemReac%NumReac(iReac2) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          ELSE
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
#if (PP_TimeDiscMethod==42)
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
# endif
              CALL DSMC_Chemistry(iPair, iReac3)
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac3) = ChemReac%NumReac(iReac3) + 1  ! for calculation of reaction rate coefficient
            END IF
# endif
          END IF
          RelaxToDo = .FALSE.
        END IF
      END IF
!############################################################################################################################### !
    CASE(18) ! only electron impact ionization possible Ar + e -> Ar(+) + e + e
! ############################################################################################################################### !
      iReac = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 1)
      IF (ChemReac%QKProcedure(iReac)) THEN
        IF ( .NOT. DSMC%ElectronicModel ) THEN
          CALL Abort(&
__STAMP__&
,'ERROR! Atomic electron shell has to be initalized.')
        END IF
        CALL QK_ImpactIonization(iPair,iReac,RelaxToDo)
      END IF
!-----------------------------------------------------------------------------------------------------------------------------------
      IF (.NOT.ChemReac%QKProcedure(iReac)) THEN
         CALL Abort(&
__STAMP__&
,'ERROR! Electron impact ionization not implemented without QK')
      END IF
!############################################################################################################################### !
    CASE(19) ! only ion recombination possible Ar(+) + e + e -> Ar + e
! ############################################################################################################################### !
      ! searching third collison partner
      IF(ChemReac%RecombParticle.EQ. 0) THEN
        IF(iPair.LT.(nPair - ChemReac%nPairForRec)) THEN
          PairForRec = nPair - ChemReac%nPairForRec
          iPart3   = Coll_pData(PairForRec)%iPart_p1
        ELSE
          iPart3   = 0
        END IF
      ELSE
        iPart3 = ChemReac%RecombParticle
      END IF
!-----------------------------------------------------------------------------------------------------------------------------------
      IF ( iPart3 .GT. 0 ) THEN
        iReac = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), PartSpecies(iPart3))
        IF(iReac.EQ.0) THEN
          iReac = ChemReac%ReactNumRecomb(PartSpecies(iPart1), PartSpecies(iPart2), PartSpecies(iPart3))
        END IF
!        IF(SpecDSMC(PartSpecies(iPart3))%InterID.NE.4)RETURN
        IF ( ChemReac%QKProcedure(iReac)  ) THEN
          IF ( .NOT. DSMC%ElectronicModel ) THEN
            CALL Abort(&
__STAMP__&
,' ERROR! Atomic electron shell has to be initalized.')
          END IF
          CALL QK_IonRecombination(iPair,iReac,iPart3,RelaxToDo,NodeVolume,NodePartNum)
!          CALL QK_IonRecombination(iPair,iReac,iPart3,RelaxToDo,iElem,NodeVolume,NodePartNum)
!-----------------------------------------------------------------------------------------------------------------------------------
        ELSE
        ! traditional Recombination
          ! Calculation of reaction probability
          CALL CalcReactionProb(iPair,iReac,ReactionProb,iPart3,NumDens)
#if (PP_TimeDiscMethod==42)
          IF (.NOT.DSMC%ReservoirRateStatistic) THEN
            ChemReac%NumReac(iReac)   = ChemReac%NumReac(iReac)   + ReactionProb  ! for calculation of reaction rate coefficient
            ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
          END IF
#endif
          CALL RANDOM_NUMBER(iRan)
          IF (ReactionProb.GT.iRan) THEN
#if (PP_TimeDiscMethod==42)
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
            IF (.NOT.DSMC%ReservoirSimuRate) THEN
#endif
              CALL DSMC_Chemistry(iPair, iReac, iPart3)
              IF(ChemReac%RecombParticle.EQ. 0) THEN
                Coll_pData(PairForRec)%NeedForRec = .TRUE.
                ChemReac%RecombParticle = Coll_pData(PairForRec)%iPart_p2
                ChemReac%nPairForRec = ChemReac%nPairForRec + 1
              ELSE
                ChemReac%RecombParticle = 0
              END IF
#if (PP_TimeDiscMethod==42)
            END IF
            IF (DSMC%ReservoirRateStatistic) THEN
              ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reaction rate coefficient
            END IF
#endif
            RelaxToDo = .FALSE.
          END IF ! ReactionProb > iRan
        END IF ! Q-K
      END IF
! ############################################################################################################################### !
    CASE(20) ! Dissociation and ionization with QK are possible
! ############################################################################################################################### !
      iReac  = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 1)
      iReac2 = ChemReac%ReactNum(PartSpecies(iPart1), PartSpecies(iPart2), 2)
      IF (ChemReac%QKProcedure(iReac).AND.ChemReac%QKProcedure(iReac2)) THEN ! both Reaction QK
        ! first pseudo reaction probability
        IF (ChemReac%DefinedReact(iReac,1,1).EQ.PartSpecies(iPart1)) THEN
          PartToExec = iPart1
          PartReac2  = iPart2
        ELSE
          PartToExec = iPart2
          PartReac2  = iPart1
        END IF
        ! Determine the collision energy (only relative translational)
        Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%cRela2
        IF(DSMC%ElectronicModel) Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(PartToExec,3)
        ! ionization level is last known energy level of species
        iQuaMax          = SpecDSMC(PartSpecies(PartToExec))%MaxElecQuant - 1
        IonizationEnergy = SpecDSMC(PartSpecies(PartToExec))%ElectronicState(2,iQuaMax)*BoltzmannConst
        ! if you have electronic levels above the ionization limit, such limits should be used instead of
        ! the pure energy comparison
        IF(Coll_pData(iPair)%Ec .GT. IonizationEnergy)THEN
          CALL CalcReactionProb(iPair,iReac,ReactionProb)
        ELSE
          ReactionProb = 0.
        END IF
        ! second pseudo reaction probability
        IF (ChemReac%DefinedReact(iReac2,1,1).EQ.PartSpecies(iPart1)) THEN
          PartToExec = iPart1
          PartReac2  = iPart2
        ELSE
          PartToExec = iPart2
          PartReac2  = iPart1
        END IF
        ! Determine the collision energy (relative translational + vibrational energy of dissociating molecule)
        Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%cRela2 &
                             + PartStateIntEn(PartToExec,1)
        ! Correction for second collision partner
        IF ((SpecDSMC(PartSpecies(PartReac2))%InterID.EQ.2).OR.(SpecDSMC(PartSpecies(PartReac2))%InterID.EQ.20)) THEN
          Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - SpecDSMC(PartSpecies(PartReac2))%EZeroPoint
        END IF
        ! Determination of the quantum number corresponding to the collision energy
        iQuaMax   = INT(Coll_pData(iPair)%Ec / ( BoltzmannConst * SpecDSMC(PartSpecies(PartToExec))%CharaTVib ) - DSMC%GammaQuant)
        ! Comparing the collision quantum number with the dissociation quantum number
        IF (iQuaMax.GT.SpecDSMC(PartSpecies(PartToExec))%DissQuant) THEN
          CALL CalcReactionProb(iPair,iReac2,ReactionProb2)
        ELSE
          ReactionProb2 = 0.
        END IF
#if (PP_TimeDiscMethod==42)
        IF (.NOT.DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac)    = ChemReac%NumReac(iReac)    + ReactionProb  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac)  = ChemReac%ReacCount(iReac)  + 1
          ChemReac%NumReac(iReac2)   = ChemReac%NumReac(iReac2)   + ReactionProb2  ! for calculation of reaction rate coefficient
          ChemReac%ReacCount(iReac2) = ChemReac%ReacCount(iReac2) + 1
        END IF
#endif
        ReacToDo = 0
        IF(ReactionProb*ReactionProb2.LE.0.0) THEN
          IF(ReactionProb.GT.0.0) THEN
            ReacToDo = iReac
          END IF
          IF(ReactionProb2.GT.0.0) THEN
            ReacToDo = iReac2
          ENDIF
        ELSE
          CALL RANDOM_NUMBER(iRan)
          IF((ReactionProb/(ReactionProb + ReactionProb2)).GT.iRan) THEN
            ReacToDo = iReac
          ELSE
            ReacToDo = iReac2
          END IF
        END IF
        IF(ReacToDo.NE.0) THEN
#if (PP_TimeDiscMethod==42)
          IF (.NOT.DSMC%ReservoirSimuRate) THEN
#endif
            CALL DSMC_Chemistry(iPair, ReacToDo)
#if (PP_TimeDiscMethod==42)
          END IF
          IF (DSMC%ReservoirRateStatistic) THEN
            ChemReac%NumReac(ReacToDo) = ChemReac%NumReac(ReacToDo) + 1  ! for calculation of reaction rate coefficient
          END IF
#endif
          RelaxToDo = .FALSE.
        END IF
      END IF
!-----------------------------------------------------------------------------------------------------------------------------------
    CASE DEFAULT
      IF(CaseOfReaction.NE.0) THEN
        CALL Abort(&
__STAMP__&
,'Error! Reaction case not defined:',CaseOfReaction)
      END IF
  END SELECT

END SUBROUTINE ReactionDecision


SUBROUTINE DSMC_calc_P_rot(iSpec, iPair, iPart, Xi_rel, ProbRot, ProbRotMax)
!===================================================================================================================================
! Calculation of probability for rotational relaxation. Different Models implemented:
! 0 - Constant Probability
! 1 - No rotational relaxation. RotRelaxProb = 0
! 2 - Boyd
! 3 - Zhang (Nonequilibrium Direction Dependent)
!===================================================================================================================================
! MODULES
  USE MOD_Globals            ,ONLY : Abort
  USE MOD_Globals_Vars       ,ONLY : Pi, BoltzmannConst
  USE MOD_DSMC_Vars          ,ONLY : SpecDSMC, Coll_pData, PartStateIntEn, DSMC, useRelaxProbCorrFactor, CollInf
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)         :: iSpec, iPair, iPart
  REAL, INTENT(IN)            :: Xi_rel
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL, INTENT(OUT)         :: ProbRot, ProbRotMax
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                      :: TransEn, RotEn, RotDOF, CorrFact           ! CorrFact: To correct sample Bias
                                                                          ! (fewer DSMC particles than natural ones)
!===================================================================================================================================

  TransEn    = Coll_pData(iPair)%Ec ! notice that during probability calculation,Collision energy only contains translational part
  RotDOF     = SpecDSMC(iSpec)%Xi_Rot
  RotEn      = PartStateIntEn(iPart,2)
  ProbRot    = 0.
  ProbRotMax = 0.

  ! calculate correction factor according to Lumpkin et al.
  ! - depending on selection procedure. As only one particle undergoes relaxation
  ! - only one RotDOF is needed (of considered species)
  IF(useRelaxProbCorrFactor) THEN
    CorrFact = 1. + RotDOF/Xi_rel
  ELSE
    CorrFact = 1.
  END IF

  ! calculate corrected probability for rotational relaxation
  IF(DSMC%RotRelaxProb.GE.0.0.AND.DSMC%RotRelaxProb.LE.1.0) THEN
    ProbRot = DSMC%RotRelaxProb * CorrFact
  ELSEIF(DSMC%RotRelaxProb.EQ.2.0) THEN ! P_rot according to Boyd (based on Parker's model)

    RotDOF = RotDOF*0.5 ! Only half of the rotational degree of freedom, because the other half is used in the relaxation
                        ! probability of the collision partner, see Boyd (doi:10.1063/1.858531)

    ProbRot = 1./SpecDSMC(iSpec)%CollNumRotInf * (1. + GAMMA(RotDOF+2.-CollInf%omegaLaux(iSpec,iSpec)) &
            / GAMMA(RotDOF+1.5-CollInf%omegaLaux(iSpec,iSpec)) * (PI**(3./2.)/2.)*(BoltzmannConst*SpecDSMC(iSpec)%TempRefRot &
            / (TransEn + RotEn) )**(1./2.) + GAMMA(RotDOF+2.-CollInf%omegaLaux(iSpec,iSpec))  &
            / GAMMA(RotDOF+1.-CollInf%omegaLaux(iSpec,iSpec)) * (BoltzmannConst*SpecDSMC(iSpec)%TempRefRot &
            / (TransEn + RotEn) ) * (PI**2./4. + PI)) &
            * CorrFact

  ELSEIF(DSMC%RotRelaxProb.EQ.3.0) THEN ! P_rot according to Zhang (NDD)
    ! if model is used for further species but N2, it should be checked if factors n = 0.5 and Cn = 1.92 are still valid
    ! (see original eq of Zhang)
    ProbRot = 1.92 * GAMMA(Xi_rel/2.) * GAMMA(RotDOF/2.) / GAMMA(Xi_rel/2.+0.5) / GAMMA(RotDOF/2.-0.5) &
            * (1 + (Xi_rel/2-0.5)*BoltzmannConst*SpecDSMC(iSpec)%TempRefRot/TransEn) * (TransEn/RotEn)**0.5 &
            * CorrFact
    ProbRotMax = MAX(ProbRot, 0.5) ! BL energy redistribution correction factor
    ProbRot    = MIN(ProbRot, 0.5)
  ELSE
    CALL Abort(&
__STAMP__&
,'Error! Model for rotational relaxation undefined:',RealInfoOpt=DSMC%RotRelaxProb)
  END IF

END SUBROUTINE DSMC_calc_P_rot


SUBROUTINE DSMC_calc_P_vib(iSpec, Xi_rel, iElem, ProbVib)
!===================================================================================================================================
! Calculation of probability for vibrational relaxation. Different Models implemented:
! 0 - Constant Probability
! 1 - No vibrational relaxation. VibRelaxProb = 0
! 2 - Boyd with correction of Abe
!===================================================================================================================================
! MODULES
  USE MOD_Globals            ,ONLY : Abort
  USE MOD_Globals_Vars       ,ONLY : BoltzmannConst
  USE MOD_DSMC_Vars          ,ONLY : SpecDSMC, DSMC, VarVibRelaxProb, useRelaxProbCorrFactor
  USE MOD_DSMC_Vars          ,ONLY : PolyatomMolDSMC

! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)       :: iSpec, iElem
REAL, INTENT(IN)          :: Xi_rel
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)         :: ProbVib
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                      :: CorrFact       ! CorrFact: To correct sample Bias
                                            ! (fewer DSMC particles than natural ones)
INTEGER                   :: iPolyatMole, iDOF
!===================================================================================================================================

  ProbVib = 0.

  ! calculate correction factor according to Gimelshein et al.
  ! - depending on selection procedure. As only one particle undergoes relaxation
  ! - only one VibDOF (GammaVib) is needed (of considered species)
  IF(useRelaxProbCorrFactor) THEN
    CorrFact = 1. + SpecDSMC(iSpec)%GammaVib/Xi_rel
  ELSE
    CorrFact = 1.
  END IF

  IF((DSMC%VibRelaxProb.GE.0.0).AND.(DSMC%VibRelaxProb.LE.1.0)) THEN
    IF (SpecDSMC(iSpec)%PolyatomicMol.AND.(DSMC%PolySingleMode)) THEN
      iPolyatMole  = SpecDSMC(iSpec)%SpecToPolyArray
      PolyatomMolDSMC(iPolyatMole)%VibRelaxProb(1) = DSMC%VibRelaxProb * (1. + PolyatomMolDSMC(iPolyatMole)%GammaVib(1)/Xi_rel)
      DO iDOF = 2, PolyatomMolDSMC(iPolyatMole)%VibDOF
        PolyatomMolDSMC(iPolyatMole)%VibRelaxProb(iDOF) = PolyatomMolDSMC(iPolyatMole)%VibRelaxProb(iDOF - 1) + DSMC%VibRelaxProb &
                                                        * (1. + PolyatomMolDSMC(iPolyatMole)%GammaVib(1)/Xi_rel)
      END DO
    ELSE
      ProbVib = DSMC%VibRelaxProb * CorrFact
    END IF
  ELSE IF(DSMC%VibRelaxProb.EQ.2.0) THEN
    ! Calculation of Prob Vib in function DSMC_calc_var_P_vib.
    ! This has to average over all collisions according to Boyd (doi:10.1063/1.858495)
    ! The average value of the cell is only taken from the vector
    ProbVib = VarVibRelaxProb%ProbVibAv(iElem, iSpec) * CorrFact
  ELSE
    CALL Abort(&
    __STAMP__&
    ,'Error! Model for vibrational relaxation undefined:',RealInfoOpt=DSMC%VibRelaxProb)
  END IF

END SUBROUTINE DSMC_calc_P_vib

SUBROUTINE DSMC_calc_var_P_vib(iSpec, jSpec, iPair, ProbVib)
!===================================================================================================================================
  ! Calculation of probability for vibrational relaxation for variable relaxation rates. This has to average over all collisions!
  ! No instantanious variable probability calculateable
!===================================================================================================================================
! MODULES
    USE MOD_Globals            ,ONLY : Abort
    USE MOD_Globals_Vars       ,ONLY : Pi, BoltzmannConst
    USE MOD_DSMC_Vars          ,ONLY : SpecDSMC, Coll_pData, CollInf
! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)       :: iPair, iSpec, jSpec
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL, INTENT(OUT)         :: ProbVib
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                      :: TempCorr, cRela
!===================================================================================================================================
  ! (i) dref changed from dref = 0.5 * (dref_1+dref_2) 
  !                  to   dref(iSpec,jSpec) which is identical to old definition (for averagedCollisionParameters=TRUE (DEFAULT))
  ! in case of averagedCollisionParameter=FALSE dref(iSpec,jSpec) contains collision specific dref see --help for details

  ! P_vib according to Boyd, corrected by Abe, only V-T transfer
  ! determine joint omega and Dref factor and rel velo
  cRela=SQRT(Coll_pData(iPair)%cRela2)
  ! calculate non-corrected probabilities
  ProbVib = 1. /SpecDSMC(iSpec)%CollNumVib(jSpec)* cRela**(3.+2.*CollInf%omegaLaux(iSpec,iSpec)) &
          * EXP(-1.*SpecDSMC(iSpec)%CharaVelo(jSpec)/cRela)
  ! calculate high temperature correction
  TempCorr = SpecDSMC(iSpec)%VibCrossSec / (SQRT(2.)*PI*CollInf%dref(iSpec,jSpec)**2.) &
           * (  CollInf%MassRed(Coll_pData(iPair)%PairType)*cRela & !**2
           / (2.*(2.-CollInf%omegaLaux(iSpec,iSpec))*BoltzmannConst*CollInf%Tref(iSpec,iSpec)))**CollInf%omegaLaux(iSpec,iSpec)
  ! determine corrected probabilities
  ProbVib = ProbVib * TempCorr / (ProbVib + TempCorr)        ! TauVib = TauVibStd + TauTempCorr
  IF(ProbVib.NE.ProbVib) THEN !If is NAN
    ProbVib=0.
    WRITE(*,*) 'WARNING: Vibrational relaxation probability is NAN and is set to zero. cRela:', cRela
    ! CALL Abort(&
    ! __STAMP__&
    ! ,'Error! Vibrational relaxation probability is NAN (cRela);',RealInfoOpt=cRela)!, jSpec, cRela
  END IF

END SUBROUTINE DSMC_calc_var_P_vib

!--------------------------------------------------------------------------------------------------!
END MODULE MOD_DSMC_Collis
