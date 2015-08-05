#include "boltzplatz.h"

MODULE MOD_PIC_Analyze
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE CalcDepositedCharge
  MODULE PROCEDURE CalcDepositedCharge
END INTERFACE

PUBLIC:: CalcDepositedCharge 
!===================================================================================================================================

CONTAINS



SUBROUTINE CalcDepositedCharge() 
!===================================================================================================================================
! calcs the deposited chrages
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,          ONLY : nElems, sJ
USE MOD_Particle_Vars,      ONLY : PDM, Species, PartSpecies 
USE MOD_Interpolation_Vars, ONLY : wGP
USE MOD_PICDepo_Vars,  ONLY : Source
USE MOD_Particle_Analyze_Vars,ONLY:ChargeCalcDone
#ifdef MPI
USE MOD_Particle_MPI_Vars,      ONLY : PartMPI
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iElem
INTEGER           :: i,j,k
REAL              :: J_N(1,0:PP_N,0:PP_N,0:PP_N)
REAL              :: Charge, ChargeLoc, PartCharge
#ifdef MPI
REAL              :: PartCharge_sum, Charge_sum
#endif
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' PERFORMING CHARGE DEPOSITION PLAUSIBILITY CHECK...'

Charge=0.
DO iElem=1,nElems
  !--- Calculate and save volume of element iElem
  ChargeLoc=0. 
  J_N(1,0:PP_N,0:PP_N,0:PP_N)=1./sJ(:,:,:,iElem)
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    ChargeLoc = ChargeLoc + wGP(i)*wGP(j)*wGP(k) * source(4,i,j,k,iElem) * J_N(1,i,j,k)
  END DO; END DO; END DO
  Charge = Charge + ChargeLoc
END DO


PartCharge=0.
DO i=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(i)) THEN
    PartCharge = PartCharge + Species(PartSpecies(i))%ChargeIC * Species(PartSpecies(i))%MacroParticleFactor
  END IF
END DO

#ifdef MPI
   CALL MPI_ALLREDUCE(PartCharge, PartCharge_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, PartMPI%COMM, IERROR)
   CALL MPI_ALLREDUCE(Charge, Charge_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, PartMPI%COMM, IERROR)
   PartCharge = PartCharge_sum
   Charge = Charge_sum
#endif
SWRITE(*,*) "On the grid deposited charge", Charge
SWRITE(*,*) "Charge by the particles:", PartCharge
SWRITE(*,*) "Absolute deposition error:", ABS(PartCharge-Charge)
SWRITE(*,*) "Relative deposition error in percent:", ABS(PartCharge-Charge)/PartCharge*100
SWRITE(UNIT_stdOut,'(A)')' CHARGE DEPOSITION PLAUSIBILITY CHECK DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
ChargeCalcDone = .TRUE.

END SUBROUTINE CalcDepositedCharge



END MODULE MOD_PIC_Analyze
