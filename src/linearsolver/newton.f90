#include "boltzplatz.h"

MODULE MOD_Newton
!===================================================================================================================================
! Contains routines for the fully implicit scheme
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

#if (PP_TimeDiscMethod==120) || (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122) 
INTERFACE ImplicitNorm
  MODULE PROCEDURE ImplicitNorm
END INTERFACE

INTERFACE FullNewton
  MODULE PROCEDURE FullNewton
END INTERFACE

PUBLIC::ImplicitNorm,FullNewton
#endif 
!===================================================================================================================================

CONTAINS

#if (PP_TimeDiscMethod==120) || (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122) 
SUBROUTINE ImplicitNorm(t,coeff,Norm_R) 
!===================================================================================================================================
! The error-norm of the fully implicit scheme is computed
! use same norm as in maxtrix-vector source; initial norm of linearsolver 
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_DG_Vars,                 ONLY:U
USE MOD_LinearSolver_Vars,       ONLY:ImplicitSource,ExplicitSource,LinSolverRHS,mass
#ifndef PP_HDG
USE MOD_DG_Vars,                 ONLY:Ut
USE MOD_DG,                      ONLY:DGTimeDerivative_weakForm
USE MOD_Equation,                ONLY:CalcSource
USE MOD_Equation_Vars,           ONLY:DoParabolicDamping,fDamping
USE MOD_TimeDisc_Vars,           ONLY:sdtCFLOne
#else /* HDG */
USE MOD_Equation,                ONLY:CalcSourceHDG
#endif /*DG*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
REAL,INTENT(IN)        :: t
REAL,INTENT(IN)        :: coeff
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(OUT)       :: Norm_R
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                   :: DeltaX ! difference between electric field and div-correction
INTEGER                :: iElem, i,j,k,iVar
REAL                   :: Norm_e, rTmp(1:8), locMass
!===================================================================================================================================

#ifndef PP_HDG
! compute error-norm-version1, non-optimized
CALL DGTimeDerivative_weakForm(t, t, 0,doSource=.FALSE.)
ImplicitSource=ExplicitSource
CALL CalcSource(t,1.,ImplicitSource)

IF(DoParabolicDamping)THEN
  rTmp(1:6)=1.0
  rTmp( 7 )=1.0-(fDamping-1.0)*coeff*sdTCFLOne
  rTmp( 8 )=1.0-(fDamping-1.0)*coeff*sdTCFLOne
ELSE
  rTmp(1:8)=1.0
END IF

Norm_R=0.
DO iElem=1,PP_nElems
  Norm_e=0.
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        locMass=mass(1,i,j,k,iElem)
        DO iVar=1,8
          DeltaX=locMass*( LinSolverRHS(iVar,i,j,k,iElem)           &
                         - rTmp(iVar)*U(iVar,i,j,k,iElem)           &
                         +     coeff*Ut(iVar,i,j,k,iElem)           &
                         + coeff*ImplicitSource(iVar,i,j,k,iElem)   )
          Norm_e = Norm_e + DeltaX*DeltaX
        END DO ! iVar=1,PP_nVar
      END DO ! i=0,PP_N
    END DO ! j=0,PP_N
  END DO ! k=0,PP_N
  Norm_R=Norm_R+Norm_e
END DO ! iElem=1,PP_nElems
#else /*HDG*/
DO iElem=1,PP_nElems
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    CALL CalcSourceHDG(i,j,k,iElem,ImplicitSource(1:PP_nVar,i,j,k,iElem))
    ImplicitSource(1:PP_nVar,i,j,k,iElem)= ImplicitSource(1:PP_nVar,i,j,k,iElem) + ExplicitSource(1:PP_nVar,i,j,k,iElem)
  END DO; END DO; END DO !i,j,k    
END DO !iElem 
Norm_R=0.
DO iElem=1,PP_nElems
  Norm_e=0.
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        DO iVar=1,PP_nVar
          DeltaX=U(iVar,i,j,k,iElem)+ImplicitSource(iVar,i,j,k,iElem)
          Norm_e = Norm_e + DeltaX*DeltaX
        END DO ! iVar=1,PP_nVar
      END DO ! i=0,PP_N
    END DO ! j=0,PP_N
  END DO ! k=0,PP_N
  Norm_R=Norm_R+Norm_e
END DO ! iElem=1,PP_nElems
#endif /*DG*/

#ifdef MPI
DeltaX=Norm_R
CALL MPI_ALLREDUCE(DeltaX,Norm_R,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,iError)
#endif

END SUBROUTINE ImplicitNorm


SUBROUTINE FullNewton(t,tStage,coeff)
!===================================================================================================================================
! Full Newton with particles and field 
! Newton:
! Init: Implicit particle step and Norm_R0
!       1) Implicit field solver
!       2) ParticleNewton
!       3) Compute Norm_R
! EisenStat-Walker is from 
! Kelly - Iterative Methods for Linear and Nonlinear Equations, p. 105 ff
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Globals_Vars,            ONLY:EpsMach
USE MOD_TimeDisc_Vars,           ONLY:iStage,ESDIRK_a,dt
#ifndef PP_HDG
USE MOD_LinearSolver,            ONLY:LinearSolver
#else
USE MOD_HDG,                     ONLY:HDG
USE MOD_DG_Vars,                 ONLY:U
#endif /*PP_HDG*/
USE MOD_LinearSolver_Vars,       ONLY:ImplicitSource, ExplicitSource,eps_LinearSolver
USE MOD_LinearSolver_Vars,       ONLY:maxFullNewtonIter,totalFullNewtonIter,totalIterLinearSolver
USE MOD_LinearSolver_Vars,       ONLY:Eps2_FullNewton,FullEisenstatWalker,FullgammaEW,DoPrintConvInfo
#ifdef PARTICLES
USE MOD_LinearSolver_Vars,       ONLY:PartRelaxationFac,PartRelaxationFac0,DoPartRelaxation,AdaptIterRelaxation0
USE MOD_Particle_Tracking,       ONLY:ParticleTracing,ParticleRefTracking
USE MOD_Particle_Tracking_vars,  ONLY:DoRefMapping
USE MOD_LinearSolver_Vars,       ONLY:Eps2PartNewton,UpdateInIter
USE MOD_Particle_Vars,           ONLY:PartIsImplicit
USE MOD_Particle_Vars,           ONLY:PartStateN,PartStage
USE MOD_Particle_Vars,           ONLY:PartState, LastPartPos, DelayTime, PEM, PDM
USE MOD_Part_RHS,                ONLY:PartVeloToImp
USE MOD_PICInterpolation,        ONLY:InterpolateFieldToSingleParticle
USE MOD_Part_MPFtools,           ONLY:StartParticleMerge
USE MOD_Particle_Analyze_Vars,   ONLY:DoVerifyCharge
USE MOD_PIC_Analyze,             ONLY:VerifyDepositedCharge
USE MOD_PICDepo,                 ONLY:Deposition
USE MOD_ParticleSolver,          ONLY:ParticleNewton
USE MOD_part_tools,              ONLY:UpdateNextFreePosition
#ifdef MPI
USE MOD_Particle_MPI,            ONLY:IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY:PartMPIExchange
#endif /*MPI*/
#endif /*PARTICLES*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
REAL,INTENT(IN)            :: t
REAL,INTENT(INOUT)         :: tStage
REAL,INTENT(INOUT)         :: coeff
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                       :: Norm_R0,Norm_R,Norm_Rold, Norm_Diff
REAL                       :: etaA,etaB,etaC,etaMax,taut
INTEGER                    :: nFullNewtonIter
#ifdef PARTICLES
INTEGER                    :: iPart,iCounter
REAL                       :: tmpFac
INTEGER                    :: AdaptIterRelaxation
#endif /*PARTICLES*/
REAL                       :: relTolerance,relTolerancePart,Criterion
LOGICAL                    :: IsConverged
#ifdef PP_HDG
INTEGER(KIND=8)            :: iter=0
#endif /*PP_HDG*/
!===================================================================================================================================

#ifdef PARTICLES
IF (t.GE.DelayTime) THEN
  IF(FullEisenstatWalker.GT.1)THEN
    relTolerancePart=0.998
  ELSE
    relTolerancePart=eps2PartNewton
  END IF
  ! now, we have an initial guess for the field  can compute the first particle movement
  CALL ParticleNewton(tstage,coeff,doParticle_In=PartIsImplicit(1:PDM%maxParticleNumber),Opt_In=.TRUE. &
                     ,AbortTol_In=relTolerancePart)
  PartRelaxationFac = PartRelaxationFac0
  ! particle relaxation betweeen old and new position
  IF(DoPartRelaxation)THEN
    AdaptIterRelaxation=AdaptIterRelaxation0
    DO iPart=1,PDM%ParticleVecLength
      IF(PartIsImplicit(iPart))THEN  
        tmpFac=(1.0-PartRelaxationFac)
        PartState(iPart,1)=PartRelaxationFac*PartState(iPart,1)+tmpFac*PartStateN(iPart,1)
        PartState(iPart,2)=PartRelaxationFac*PartState(iPart,2)+tmpFac*PartStateN(iPart,2)
        PartState(iPart,3)=PartRelaxationFac*PartState(iPart,3)+tmpFac*PartStateN(iPart,3)
        PartState(iPart,4)=PartRelaxationFac*PartState(iPart,4)+tmpFac*PartStateN(iPart,4)
        PartState(iPart,5)=PartRelaxationFac*PartState(iPart,5)+tmpFac*PartStateN(iPart,5)
        PartState(iPart,6)=PartRelaxationFac*PartState(iPart,6)+tmpFac*PartStateN(iPart,6)
        DO iCounter=1,iStage-1
          tmpFac=tmpFac*dt*ESDIRK_a(iStage-1,iCounter)
          PartState(iPart,1) = PartState(iPart,1) + tmpFac*PartStage(iPart,1,iCounter)
          PartState(iPart,2) = PartState(iPart,2) + tmpFac*PartStage(iPart,2,iCounter)
          PartState(iPart,3) = PartState(iPart,3) + tmpFac*PartStage(iPart,3,iCounter)
          PartState(iPart,4) = PartState(iPart,4) + tmpFac*PartStage(iPart,4,iCounter)
          PartState(iPart,5) = PartState(iPart,5) + tmpFac*PartStage(iPart,5,iCounter)
          PartState(iPart,6) = PartState(iPart,6) + tmpFac*PartStage(iPart,6,iCounter)
        END DO
      END IF ! ParticleInside
    END DO ! iPart
  END IF ! PartRelaxationFac>0
  ! move particle, if not already done, here, a reduced list could be again used, but a different list...
  ! required to get the correct deposition
#ifdef MPI
  ! open receive buffer for number of particles
  CALL IRecvNbofParticles()
  ! here: could use deposition as hiding, not done yet
  IF(DoPartRelaxation)THEN
    IF(DoRefMapping)THEN
      ! input value: which list:DoPartInNewton or PDM%ParticleInisde?
      CALL ParticleRefTracking(doParticle_In=PartisImplicit(1:PDM%ParticleVecLength)) 
    ELSE
      ! input value: which list:DoPartInNewton or PDM%ParticleInisde?
      CALL ParticleTracing(doParticle_In=PartisImplicit(1:PDM%ParticleVecLength)) 
    END IF
  END IF
  DO iPart=1,PDM%ParticleVecLength
    IF(PartIsImplicit(iPart))THEN
      IF(.NOT.PDM%ParticleInside(iPart)) PartisImplicit(iPart)=.FALSE.
    END IF
  END DO

  ! send number of particles
  CALL SendNbOfParticles(doParticle_In=PartIsImplicit(1:PDM%ParticleVecLength))
  ! finish communication of number of particles and send particles
  CALL MPIParticleSend()
  ! finish communication
  CALL MPIParticleRecv()
  ! ALWAYS require
  PartMPIExchange%nMPIParticles=0
#endif /*MPI*/
  ! map particle from gamma v to v
  CALL PartVeloToImp(VeloToImp=.FALSE.,doParticle_In=PartIsImplicit(1:PDM%ParticleVecLength))
  ! compute particle source terms on field solver of implicit particles :)
  CALL Deposition(doInnerParts=.TRUE.,doParticle_In=PartIsImplicit(1:PDM%ParticleVecLength))
  CALL Deposition(doInnerParts=.FALSE.,doParticle_In=PartIsImplicit(1:PDM%ParticleVecLength))
  ! map particle from v to gamma v
  CALL PartVeloToImp(VeloToImp=.TRUE.,doParticle_In=PartIsImplicit(1:PDM%ParticleVecLength))
END IF
#endif /*PARTICLES*/

CALL ImplicitNorm(tStage,coeff,Norm_R0)
Norm_R=Norm_R0
Norm_Diff=HUGE(1.0)
IF(DoPrintConvInfo.AND.MPIRoot) WRITE(*,*) 'Norm_R0',Norm_R0
IF(FullEisenstatWalker.GT.0)THEN
  etaMax=0.9999
  taut  =epsMach+eps2_FullNewton*Norm_R0
      !SWRITE(*,*) 'taut ', taut
END IF

nFullNewtonIter=0
IsConverged=.FALSE.
DO WHILE ((nFullNewtonIter.LE.maxFullNewtonIter).AND.(.NOT.IsConverged))
  nFullNewtonIter = nFullNewtonIter+1
  IF(FullEisenstatWalker.GT.0)THEN
    IF(nFullNewtonIter.EQ.1)THEN
      !etaA=etaMax
      !etaB=etaMax
      !etaC=etaMax
      relTolerance=etaMax
    ELSE
      etaA=FullgammaEW*Norm_R/Norm_Rold
      !SWRITE(*,*) 'etaA ', etaA
      etaB=MIN(etaMax,etaA)
      !SWRITE(*,*) 'etaB ', etaB
      Criterion  =FullGammaEW*relTolerance*relTolerance
      !SWRITE(*,*) 'criterion ', Criterion
      IF(Criterion.LT.0.1)THEN
        etaC=MIN(etaMax,etaA)
      ELSE
        etaC=MIN(etaMax,MAX(etaA,Criterion))
      END IF
      !SWRITE(*,*) 'etaC ', etaC
      !relTolerance=MIN(MIN(etaMax,MAX(etaC,0.5*taut/Norm_R)),FullgammaEW*relTolerance)
      relTolerance=MIN(etaMax,MAX(etaC,0.5*taut/Norm_R))
    END IF
  ELSE
    relTolerance=eps_LinearSolver
  END IF
  ! solve field to new stage 
  ImplicitSource=ExplicitSource
#ifndef PP_HDG
  CALL LinearSolver(tStage,coeff,relTolerance)
#else
  CALL HDG(tStage,U,iter)
#endif /*HDG*/

#ifdef PARTICLES
  IF (t.GE.DelayTime) THEN
    !DO iPart=1,PDM%ParticleVecLength
    !  IF(PartIsImplicit(iPart))THEN
    !    LastPartPos(iPart,1)=PartState(iPart,1)
    !    LastPartPos(iPart,2)=PartState(iPart,2)
    !    LastPartPos(iPart,3)=PartState(iPart,3)
    !    PEM%lastElement(iPart)=PEM%Element(iPart)
    !  END IF ! PartIsImplicit
    !END DO ! iPart

    ! now, we have an initial guess for the field  can compute the first particle movement
    IF(FullEisenstatWalker.GT.1)THEN
      relTolerancePart=relTolerance*relTolerance
    ELSE
      relTolerancePart=eps2PartNewton
    END IF
    CALL ParticleNewton(tstage,coeff,doParticle_In=PartIsImplicit(1:PDM%maxParticleNumber),Opt_In=.TRUE. &
                       ,AbortTol_In=relTolerancePart)
    ! particle relaxation betweeen old and new position
    IF(DoPartRelaxation)THEN
      DO iPart=1,PDM%ParticleVecLength
        IF(PartIsImplicit(iPart))THEN  
          tmpFac=(1.0-PartRelaxationFac)
          PartState(iPart,1)=PartRelaxationFac*PartState(iPart,1)+tmpFac*PartStateN(iPart,1)
          PartState(iPart,2)=PartRelaxationFac*PartState(iPart,2)+tmpFac*PartStateN(iPart,2)
          PartState(iPart,3)=PartRelaxationFac*PartState(iPart,3)+tmpFac*PartStateN(iPart,3)
          PartState(iPart,4)=PartRelaxationFac*PartState(iPart,4)+tmpFac*PartStateN(iPart,4)
          PartState(iPart,5)=PartRelaxationFac*PartState(iPart,5)+tmpFac*PartStateN(iPart,5)
          PartState(iPart,6)=PartRelaxationFac*PartState(iPart,6)+tmpFac*PartStateN(iPart,6)
          DO iCounter=1,iStage-1
            tmpFac=tmpFac*dt*ESDIRK_a(iStage-1,iCounter)
            PartState(iPart,1) = PartState(iPart,1) + tmpFac*PartStage(iPart,1,iCounter)
            PartState(iPart,2) = PartState(iPart,2) + tmpFac*PartStage(iPart,2,iCounter)
            PartState(iPart,3) = PartState(iPart,3) + tmpFac*PartStage(iPart,3,iCounter)
            PartState(iPart,4) = PartState(iPart,4) + tmpFac*PartStage(iPart,4,iCounter)
            PartState(iPart,5) = PartState(iPart,5) + tmpFac*PartStage(iPart,5,iCounter)
            PartState(iPart,6) = PartState(iPart,6) + tmpFac*PartStage(iPart,6,iCounter)
          END DO
        END IF ! ParticleInside
      END DO ! iPart
    END IF ! PartRelaxationFac>0
    ! move particle, if not already done, here, a reduced list could be again used, but a different list...
#ifdef MPI
    ! open receive buffer for number of particles
    CALL IRecvNbofParticles()
    ! here: could use deposition as hiding, not done yet
    IF(DoPartRelaxation)THEN
      IF(DoRefMapping)THEN
        ! input value: which list:DoPartInNewton or PDM%ParticleInisde?
        CALL ParticleRefTracking(doParticle_In=PartisImplicit(1:PDM%ParticleVecLength)) 
      ELSE
        ! input value: which list:DoPartInNewton or PDM%ParticleInisde?
        CALL ParticleTracing(doParticle_In=PartisImplicit(1:PDM%ParticleVecLength)) 
      END IF
    END IF
    DO iPart=1,PDM%ParticleVecLength
      IF(PartIsImplicit(iPart))THEN
        IF(.NOT.PDM%ParticleInside(iPart)) PartisImplicit(iPart)=.FALSE.
      END IF
    END DO
    ! send number of particles
    CALL SendNbOfParticles(doParticle_In=PartIsImplicit(1:PDM%ParticleVecLength))
    ! finish communication of number of particles and send particles
    CALL MPIParticleSend()
    ! finish communication
    CALL MPIParticleRecv()
    PartMPIExchange%nMPIParticles=0
#endif /*MPI*/

    ! map particle from gamma v to v
    CALL PartVeloToImp(VeloToImp=.FALSE.,doParticle_In=PartIsImplicit(1:PDM%ParticleVecLength))
    ! compute particle source terms on field solver of implicit particles :)
    CALL Deposition(doInnerParts=.TRUE.,doParticle_In=PartIsImplicit(1:PDM%ParticleVecLength))
    CALL Deposition(doInnerParts=.FALSE.,doParticle_In=PartIsImplicit(1:PDM%ParticleVecLength))
    IF(DoVerifyCharge) CALL VerifyDepositedCharge()
    ! and map back
    CALL PartVeloToImp(VeloToImp=.TRUE.,doParticle_In=PartIsImplicit(1:PDM%ParticleVecLength))
  END IF
#endif /*PARTICLES*/

  Norm_Rold=Norm_R
  CALL ImplicitNorm(tStage,coeff,Norm_R)
  IF(DoPrintConvInfo.AND.MPIRoot) WRITE(*,*) 'iter,Norm_R,rel,abort',nFullNewtonIter,Norm_R,Norm_R/Norm_R0,relTolerance

  Norm_Diff=ABS(Norm_Rold-Norm_R)
  IF((Norm_R.LT.Norm_R0*Eps2_FullNewton).OR.(Norm_Diff.LT.Norm_R0*eps2_FullNewton)) IsConverged=.TRUE.
  IF(DoPartRelaxation)THEN
    IF(MOD(nFullNewtonIter,AdaptIterRelaxation).EQ.0)THEN
      IF(Norm_Rold.GT.Norm_R)THEN
        PartRelaxationFac=MAX(PartRelaxationFac*1.55,1.0)
        !PartRelaxationFac=PartRelaxationFac*2
        !IF(PartRelaxationFac.GE.1.0) DoPartRelaxation=.FALSE.
      ELSE
         PartRelaxationFac=MAX(PartRelaxationFac/2,0.001)
      END IF
      AdaptIterRelaxation=MAX(INT(AdaptIterRelaxation*PartRelaxationFac0/PartRelaxationFac),AdaptIterRelaxation0)
    END IF
  END IF ! DoPartRelaxation

#ifdef PARTICLES
  IF((.NOT.IsConverged).AND.(MOD(nFullNewtonIter,UpdateInIter).EQ.0)) CALL UpdateNextFreePosition()
#endif /*PARTICLES*/
END DO ! funny pseudo Newton for all implicit

!IF(PartRelaxationFac0.NE.0) DoPartRelaxation=.TRUE.

totalFullNewtonIter=TotalFullNewtonIter+nFullNewtonIter
IF(nFullNewtonIter.GE.maxFullNewtonIter)THEN
  SWRITE(UNIT_StdOut,'(A)') " Implicit scheme is not converged!"
  SWRITE(UNIT_StdOut,'(A,E20.14,5x,E20.14)') ' NormDiff and NormDiff/Norm_R0: ',Norm_Diff, Norm_Diff/Norm_R0
  SWRITE(UNIT_StdOut,'(A,E20.14,5x,E20.14)') ' Norm_R/Norm_R0               : ',Norm_R/Norm_R0
  IF(MPIRoot) CALL abort(&
 __STAMP__&
   ,' Outer-Newton of semi-fully implicit scheme is running into infinity.',nFullNewtonIter,Norm_R/Norm_R0)
END IF

IF(DoPrintConvInfo.AND.MPIRoot) WRITE(*,*) 'TotalIterlinearsolver',TotalIterlinearSolver

END SUBROUTINE FullNewton
#endif

END MODULE MOD_Newton
