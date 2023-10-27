!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (piclas.boltzplatz.eu/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

MODULE MOD_ParticleSolver
!===================================================================================================================================
! Contains routines to compute the riemann (Advection, Diffusion) for a given Face
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

#if defined(PARTICLES)
#if defined(IMPA)
INTERFACE ParticleNewton
  MODULE PROCEDURE ParticleNewton
END INTERFACE

INTERFACE SelectImplicitParticles
  MODULE PROCEDURE SelectImplicitParticles
END INTERFACE
#endif /*IMPA*/
#if defined(IMPA) || defined(ROS)
INTERFACE InitPartSolver
  MODULE PROCEDURE InitPartSolver
END INTERFACE

INTERFACE FinalizePartSolver
  MODULE PROCEDURE FinalizePartSolver
END INTERFACE

INTERFACE Particle_GMRES
  MODULE PROCEDURE Particle_GMRES
END INTERFACE
#endif /*IMPA or ROS*/

#if defined(IMPA) || defined(ROS)
PUBLIC:: InitPartSolver,FinalizePartSolver
PUBLIC:: Particle_GMRES
#endif /*IMPA or ROS*/
#ifdef IMPA
PUBLIC:: ParticleNewton
PUBLIC:: SelectImplicitParticles
#endif /*IMPA*/
#endif /*PARTICLES*/
!===================================================================================================================================

CONTAINS

#if defined(PARTICLES)
#if defined(IMPA) || defined(ROS)
SUBROUTINE InitPartSolver()
!===================================================================================================================================
! read in and allocation of required global variables for implicit particle treatment
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_PreProc
USE MOD_ReadInTools       ,ONLY: GETINT,GETREAL,GETLOGICAL
USE MOD_Particle_Vars     ,ONLY: PDM
USE MOD_LinearSolver_Vars
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars  ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: allocstat
REAL                        :: scaleps
!===================================================================================================================================

LBWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE SOLVER...'

#if defined(IMPA)
Eps2PartNewton     =GETREAL('EpsPartNewton','0.001')
Eps2PartNewton     =Eps2PartNewton**2
nPartNewtonIter    =GETINT('nPartNewtonIter','20')
FreezePartInNewton =GETINT('FreezePartInNewton','1')
EisenstatWalker    =GETLOGICAL('EisenstatWalker','.FALSE.')
PartgammaEW        =GETREAL('PartgammaEW','0.9')
nPartNewton        =0
PartNewtonLinTolerance  = GETLOGICAL('PartNewtonLinTolerance','.FALSE.')
#elif defined(ROS)
EisenstatWalker = .FALSE.
#endif /*IMPA*/

#if !(USE_HDG)
EpsPartLinSolver   =GETREAL('EpsPartLinSolver','0.')
IF(EpsPartLinSolver.EQ.0.) EpsPartLinSolver=Eps_LinearSolver
#else
EpsPartLinSolver   =GETREAL('EpsPartLinSolver','1e-3')
#endif /*DG*/
nKDIMPart            = GETINT('nKDIMPart','6')

! read in by both
scaleps=GETREAL('scaleps','1.')
! rEps0 = scaleps * 1.E-8
rEps0=scaleps*SQRT(EPSILON(0.0))
srEps0=1./rEps0

ALLOCATE(PartXK(1:6,1:PDM%maxParticleNumber),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
__STAMP__&
,'Cannot allocate PartXK')

ALLOCATE(R_PartXK(1:6,1:PDM%maxParticleNumber),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
__STAMP__&
,'Cannot allocate R_PartXK')

#ifdef IMPA
DoFullNewton = GETLOGICAL('DoFullNewton','.FALSE.')
IF(DoFullNewton)THEN
  LBWRITE(UNIT_stdOut,'(A)') ' Using a full Newton for Particle and Field instead of Piccardi-Iteration.'
  nPartNewtonIter=1
  LBWRITE(UNIT_stdOut,'(A,I0)') ' Setting nPartNewtonIter to: ', nPartNewtonIter
END IF

PartImplicitMethod =GETINT('Part-ImplicitMethod','0')
#endif /*IMPA*/

END SUBROUTINE InitPartSolver
#endif


#if IMPA
SUBROUTINE SelectImplicitParticles()
!===================================================================================================================================
! select if particle is treated implicitly or explicitly, has to be called, after particle are created/emitted
! currently only one criterion is used: the species
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Vars     ,ONLY: Species,PartSpecies,PartIsImplicit,PDM,Pt,PartState
USE MOD_Linearsolver_Vars ,ONLY: PartImplicitMethod
USE MOD_TimeDisc_Vars     ,ONLY: dt,nRKStages,iter
USE MOD_Globals_Vars      ,ONLY: c2_inv
USE MOD_LinearSolver_Vars ,ONLY: DoPrintConvInfo
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER     :: iPart
REAL        :: NewVelo(3),Vabs,PartGamma
INTEGER     :: nImp,nExp
!===================================================================================================================================

PartIsImplicit=.FALSE.
!IF(time.LT.3e-8)THEN
!  RETURN
!END IF
SELECT CASE(PartImplicitMethod)
CASE(0) ! depending on species
  DO iPart=1,PDM%ParticleVecLength
    IF(.NOT.PDM%ParticleInside(iPart)) CYCLE
    IF(Species(PartSpecies(iPart))%IsImplicit) PartIsImplicit(iPart)=.TRUE.
  END DO ! iPart
CASE(1) ! selection after simplified, linear push
  IF(iter.EQ.0)THEN
    DO iPart=1,PDM%ParticleVecLength
      IF(.NOT.PDM%ParticleInside(iPart)) CYCLE
      PartIsImplicit(iPart)=.TRUE.
    END DO ! iPart
  ELSE
    DO iPart=1,PDM%ParticleVecLength
      IF(.NOT.PDM%ParticleInside(iPart)) CYCLE
      NewVelo=PartState(4:6,iPart)+dt/REAL(nRKStages-1)*Pt(1:3,iPart)
      Vabs   =DOT_PRODUCT(NewVelo,NewVelo)
      IF(Vabs*c2_inv.GT.0.9) PartIsImplicit(iPart)=.TRUE.
    END DO ! iPart
  END IF
CASE(2) ! if gamma exceeds a certain treshold
  IF(iter.EQ.0)THEN
    DO iPart=1,PDM%ParticleVecLength
      IF(.NOT.PDM%ParticleInside(iPart)) CYCLE
      PartIsImplicit(iPart)=.TRUE.
    END DO ! iPart
  ELSE
    DO iPart=1,PDM%ParticleVecLength
      IF(.NOT.PDM%ParticleInside(iPart)) CYCLE
      NewVelo=PartState(4:6,iPart)
      Vabs   =DOT_PRODUCT(NewVelo,NewVelo)
      PartGamma=1.0-Vabs*c2_inv
      PartGamma=1.0/SQRT(PartGamma)
      IF(PartGamma.GT.0.3) PartIsImplicit(iPart)=.TRUE.
    END DO ! iPart
  END IF
! CASE(3)
! use the dense output to compute error, if to large, switch to implicit
CASE DEFAULT
  IF(MPIRoot)  CALL abort(&
__STAMP__&
,' Method to select implicit particles is not implemented!')
END SELECT

IF(DoPrintConvInfo)THEN
  nImp=0
  nExp=0
  DO iPart=1,PDM%ParticleVecLength
    IF(.NOT.PDM%ParticleInside(iPart)) CYCLE
    IF(PartIsImplicit(iPart)) nImp=nImp+1
    IF(.NOT.PartIsImplicit(iPart)) nExp=nExp+1
  END DO
#if USE_MPI
  IF(MPIRoot)THEN
    CALL MPI_REDUCE(MPI_IN_PLACE,nExp,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_PICLAS, IERROR)
    CALL MPI_REDUCE(MPI_IN_PLACE,nImp,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_PICLAS, IERROR)
  ELSE
    CALL MPI_REDUCE(nExp       ,iPart,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_PICLAS, IERROR)
    CALL MPI_REDUCE(nImp       ,iPart,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_PICLAS, IERROR)
  END IF
#endif /*USE_MPI*/
  SWRITE(UNIT_StdOut,'(A,I0,x,I0)') ' Particles explicit/implicit ', nExp, nImp
END IF

END SUBROUTINE SelectImplicitParticles


SUBROUTINE ParticleNewton(t,coeff,Mode,doParticle_In,opt_In,AbortTol_In)
!===================================================================================================================================
! Allocate global variable
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_LinearSolver_Vars      ,ONLY: PartXK,R_PartXK,PartNewtonRelaxation
USE MOD_Particle_Vars          ,ONLY: PartQ,F_PartX0,F_PartXk,Norm_F_PartX0,Norm_F_PartXK,Norm_F_PartXK_old,DoPartInNewton    &
                                     ,PartState, Pt, LastPartPos, PEM, PDM, PartLorentzType,PartDeltaX,PartDtFrac,PartStateN  &
                                     ,PartMeshHasReflectiveBCs
USE MOD_LinearOperator         ,ONLY: PartVectorDotProduct
USE MOD_Particle_Tracing       ,ONLY: ParticleTracing
USE MOD_Particle_RefTracking   ,ONLY: ParticleRefTracking
USE MOD_Particle_TriaTracking  ,ONLY: ParticleTriaTracking
USE MOD_Part_RHS               ,ONLY: CalcPartRHS
#if USE_MPI
USE MOD_Particle_MPI           ,ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers     ,ONLY: LBStartTime,LBPauseTime,LBSplitTime
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/
USE MOD_LinearSolver_vars      ,ONLY: Eps2PartNewton,nPartNewton, PartgammaEW,nPartNewtonIter,DoPrintConvInfo
USE MOD_Part_RHS               ,ONLY: PartRHS
USE MOD_Globals_Vars           ,ONLY: c2_inv
USE MOD_PICInterpolation       ,ONLY: InterpolateFieldToSingleParticle
USE MOD_PICInterpolation_Vars  ,ONLY: FieldAtParticle
#ifdef CODE_ANALYZE
USE MOD_Particle_Tracking_Vars ,ONLY: PartOut,MPIRankOut
#endif /*CODE_ANALYZE*/
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: t,coeff
LOGICAL,INTENT(INOUT),OPTIONAL:: doParticle_In(1:PDM%maxParticleNumber)
LOGICAL,INTENT(IN),OPTIONAL   :: opt_In
REAL,INTENT(IN),OPTIONAL      :: AbortTol_In
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)           :: Mode
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                      :: opt
REAL                         :: time
INTEGER                      :: iPart
INTEGER                      :: nInnerPartNewton = 0
REAL                         :: AbortCritLinSolver,gammaA,gammaB
REAL                         :: Pt_tmp(1:6)
LOGICAL                      :: DoNewton,reMap
REAL                         :: AbortTol
REAL                         :: LorentzFacInv
REAL                         :: n_loc(1:3)
INTEGER:: counter
#if USE_LOADBALANCE
REAL                         :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

time = t+coeff
opt=.TRUE.
IF(PRESENT(opt_In)) opt=Opt_in

! quasi-newton:
! hold the system
! real newton:
! update Pt at each iteration
IF(PRESENT(DoParticle_IN))THEN
  DoPartInNewton=DoParticle_In
ELSE
  DoPartInNewton(1:PDM%maxParticleNumber)=PDM%ParticleInside(1:PDM%maxParticleNumber)
END IF

IF(PRESENT(AbortTol_In))THEN
  AbortTol=AbortTol_In
ELSE
  AbortTol=SQRT(Eps2PartNewton)
END IF

IF(opt)THEN ! compute zero state
  ! whole pt array
  DO iPart=1,PDM%ParticleVecLength
    IF(DoPartInNewton(iPart))THEN
      ! compute Lorentz force at particle's position
      CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle(1:6,iPart))
      reMap=.FALSE.
      IF(PartMeshHasReflectiveBCs)THEN
        IF(SUM(ABS(PEM%NormVec(1:3,iPart))).GT.0.)THEN
          n_loc=PEM%NormVec(1:3,iPart)
          ! particle is actually located outside, hence, it moves in the mirror field
          FieldAtParticle(1:3,iPart)=FieldAtParticle(1:3,iPart)-2.*DOT_PRODUCT(FieldAtParticle(1:3,iPart),n_loc)*n_loc
          FieldAtParticle(4:6,iPart)=FieldAtParticle(4:6,iPart)!-2.*DOT_PRODUCT(FieldAtParticle(4:6,iPart),n_loc)*n_loc
          ! and of coarse, the velocity has to be back-rotated, because the particle has not hit the wall
          reMap=.TRUE.
          PEM%NormVec(1:3,iPart)=0.
        END IF
      END IF
      IF(PEM%PeriodicMoved(iPart)) THEN
        reMap=.TRUE.
        PEM%PeriodicMoved(iPart)=.FALSE.
      END IF
      IF(PartMeshHasReflectiveBCs) PEM%NormVec(1:3,iPart)=0.
      PEM%PeriodicMoved(iPart)=.FALSE.
      IF(reMap)THEN
        PartState(1:6,iPart)=PartXK(1:6,iPart)+PartDeltaX(1:6,iPart)
      END IF
      ! update the last part pos and element for particle movement
      LastPartPos(1,iPart)=PartStateN(1,iPart)
      LastPartPos(2,iPart)=PartStateN(2,iPart)
      LastPartPos(3,iPart)=PartStateN(3,iPart)
      PEM%LastGlobalElemID(iPart)=PEM%GlobalElemID(iPart)
      ! HERE: rotate part to partstate back
      IF(PartLorentzType.EQ.5)THEN
        LorentzFacInv=1.0/SQRT(1.0+DOT_PRODUCT(PartState(4:6,iPart),PartState(4:6,iPart))*c2_inv)
        CALL PartRHS(iPart,FieldAtParticle(1:6,iPart),Pt(1:3,iPart),LorentzFacInv)
      ELSE
        LorentzFacInv = 1.0
        CALL PartRHS(iPart,FieldAtParticle(1:6,iPart),Pt(1:3,iPart))
      END IF ! PartLorentzType.EQ.5
      ! PartStateN has to be exchanged by PartQ
      Pt_tmp(1:3) = LorentzFacInv*PartState(4:6,iPart)
      Pt_tmp(4) = Pt(1,iPart)
      Pt_tmp(5) = Pt(2,iPart)
      Pt_tmp(6) = Pt(3,iPart)
      F_PartX0(1:6,iPart) =   PartState(1:6,iPart)-PartQ(1:6,iPart)-PartDtFrac(iPart)*coeff*Pt_tmp(1:6)
      PartXK(1:6,iPart)   =   PartState(1:6,iPart)
      R_PartXK(1:6,iPart) =   Pt_tmp(1:6)
      F_PartXK(1:6,iPart) =   F_PartX0(1:6,iPart)
      CALL PartVectorDotProduct(F_PartX0(:,iPart),F_PartX0(:,iPart),Norm_F_PartX0(iPart))
      Norm_F_PartX0(iPart)=SQRT(Norm_F_PartX0(iPart))
      IF (Norm_F_PartX0(iPart).LT.6E-16) THEN ! do not iterate, as U is already the implicit solution
        Norm_F_PartXk(iPart)=TINY(1.)
        DoPartInNewton(iPart)=.FALSE.
      ELSE ! we need iterations
        Norm_F_PartXk(iPart)=Norm_F_PartX0(iPart)
      END IF
    END IF ! ParticleInside
  END DO ! iPart
ELSE
  DO iPart=1,PDM%ParticleVecLength
    IF(DoPartInNewton(iPart))THEN
      ! update the last part pos and element for particle movement
      !LastPartPos(1,iPart)=StagePartPos(iPart,1)
      !LastPartPos(2,iPart)=StagePartPos(iPart,2)
      !LastPartPos(3,iPart)=StagePartPos(iPart,3)
      !PEM%LastGlobalElemID(iPart)=PEM%StageElement(iPart)
      LastPartPos(1,iPart)=PartStateN(1,iPart)
      LastPartPos(2,iPart)=PartStateN(2,iPart)
      LastPartPos(3,iPart)=PartStateN(3,iPart)
      PEM%LastGlobalElemID(iPart)=PEM%GlobalElemID(iPart)
      reMap=.FALSE.
      IF(PartMeshHasReflectiveBCs)THEN
        IF(SUM(ABS(PEM%NormVec(1:3,iPart))).GT.0.)THEN
          reMap=.TRUE.
          PEM%NormVec(1:3,iPart)=0.
        END IF
      END IF
      IF(PEM%PeriodicMoved(iPart)) reMap=.TRUE.
      IF(reMap)THEN
        PartState(1:6,iPart)=PartXK(1:6,iPart)+PartDeltaX(1:6,iPart)
      END IF
      PEM%PeriodicMoved(iPart)=.FALSE.
    END IF ! ParticleInside
  END DO ! iPart
END IF
#if USE_LOADBALANCE
CALL LBPauseTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/

DoNewton=.FALSE.
Mode=0
IF(ANY(DoPartInNewton)) DoNewton=.TRUE.
#if USE_MPI
!set T if at least 1 proc has to do newton
CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoNewton,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/

IF(DoPrintConvInfo)THEN
  ! newton per particle
  Counter=0
  DO iPart=1,PDM%ParticleVecLength
    IF(DoPartInNewton(iPart))THEN
      Counter=Counter+1
    END IF ! ParticleInside
  END DO ! iPart
#if USE_MPI
  !set T if at least 1 proc has to do newton
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,Counter,1,MPI_INTEGER,MPI_SUM,MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/
  SWRITE(UNIT_StdOut,'(A,I0)') ' Initial particle number in newton: ',Counter
END IF

IF(.NOT.DoNewton)THEN
  Mode=1
  RETURN
END IF

AbortCritLinSolver=0.999
nInnerPartNewton=0
DO WHILE((DoNewton) .AND. (nInnerPartNewton.LT.nPartNewtonIter))  ! maybe change loops, finish particle after particle?
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  nInnerPartNewton=nInnerPartNewton+1
  IF(DoPrintConvInfo)THEN
    SWRITE(UNIT_StdOut,'(A,I0)') ' Particle Newton iteration: ',nInnerPartNewton
  END IF
  DO iPart=1,PDM%ParticleVecLength
    IF(DoPartInNewton(iPart))THEN
      ! set abort crit
      IF (nInnerPartNewton.EQ.1) THEN
        AbortCritLinSolver=0.999
      ELSE
        gammaA = PartgammaEW*(Norm_F_PartXk(iPart)**2)/(Norm_F_PartXk_old(iPart)**2) ! square of norms
        IF (PartgammaEW*AbortCritLinSolver*AbortCritLinSolver < 0.1) THEN
          gammaB = MIN(0.999,gammaA)
        ELSE
          gammaB = MIN(0.999, MAX(gammaA,PartgammaEW*AbortCritLinSolver*AbortCritLinSolver))
        ENDIF
        AbortCritLinSolver = MIN(0.999,MAX(gammaB,0.5*(AbortTol)/(Norm_F_PartXk(iPart))))
      END IF
      Norm_F_PartXk_old(iPart)=Norm_F_PartXk(iPart)
      CALL Particle_GMRES(t,coeff,iPart,-F_PartXK(:,iPart),(Norm_F_PartXk(iPart)),AbortCritLinSolver,PartDeltaX(1:6,iPart))
      PartDeltaX(1:6,iPart)=PartDeltaX(1:6,iPart)*PartNewtonRelaxation
      ! everything else is done in Particle_Armijo
    END IF ! ParticleInside
  END DO ! iPart
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/

  ! DeltaX is going to be global
  CALL Particle_Armijo(t,coeff,AbortTol,nInnerPartNewton)

  ! check if all particles are converged
  DoNewton=.FALSE.
  IF(ANY(DoPartInNewton)) DoNewton=.TRUE.
#if USE_MPI
  !set T if at least 1 proc has to do newton
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoNewton,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/
  IF(DoPrintConvInfo)THEN
    Counter=0
    DO iPart=1,PDM%ParticleVecLength
      IF(DoPartInNewton(iPart))THEN
        Counter=Counter+1
      END IF ! ParticleInside
    END DO ! iPart
#if USE_MPI
    !set T if at least 1 proc has to do newton
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,Counter,1,MPI_INTEGER,MPI_SUM,MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/
  END IF
END DO

IF(DoPrintConvInfo)THEN
  IF (nInnerPartNewton.EQ.nPartNewtonIter) THEN
    SWRITE(UNIT_stdOut,'(A,2x,I10,2x,I10)') ' PartNewton-not done!',nInnerPartNewton,Counter
!    DO iPart=1,PDM%ParticleVecLength
!      IF(DoPartInNewton(iPart))THEN
!        SWRITE(UNIT_stdOut,'(A20,2x,I10)') ' Failed Particle: ',iPart
!        SWRITE(UNIT_stdOut,'(A20,6(2x,E24.12))') ' Failed Position: ',PartState(1:6,iPart)
!        SWRITE(UNIT_stdOut,'(A20,2x,E24.12)') ' Relative Norm:   ', Norm_F_PartXK(iPart)/Norm_F_PartX0(iPart)
!      END IF ! ParticleInside
!    END DO ! iPart
  ELSE
    SWRITE(UNIT_stdOut,'(A20,2x,I10,2x,I10)') ' PartNewton:',nInnerPartNewton,Counter
  END IF
END IF
nPartNewton=nPartNewton+nInnerPartNewton

END SUBROUTINE ParticleNewton
#endif /*IMPA*/

#if defined(ROS) || defined(IMPA)
SUBROUTINE Particle_GMRES(t,coeff,PartID,B,Norm_B,AbortCrit,DeltaX)
!===================================================================================================================================
! Uses matrix free to solve the linear system
! Attention: We use DeltaX=0 as our initial guess   ! why not Un??
!            X0 is allready stored in U
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_LinearSolver_Vars,    ONLY: epsPartlinSolver,TotalPartIterLinearSolver
USE MOD_LinearSolver_Vars,    ONLY: nKDimPart,nRestarts,nPartInnerIter,EisenstatWalker
USE MOD_LinearOperator,       ONLY: PartMatrixVector, PartVectorDotProduct
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: t,coeff,Norm_B
REAL,INTENT(IN)   :: B(1:6)
REAL,INTENT(IN)   :: AbortCrit
REAL,INTENT(OUT)  :: DeltaX(1:6)
INTEGER,INTENT(IN):: PartID
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: AbortCritLoc
REAL              :: V(1:6,1:nKDimPart)
REAL              :: W(1:6)
REAL              :: R0(1:6)
REAL              :: Gam(1:nKDimPart+1),C(1:nKDimPart),S(1:nKDimPart),H(1:nKDimPart+1,1:nKDimPart+1),Alp(1:nKDimPart+1)
REAL              :: Norm_R0,Resu,Temp,Bet
INTEGER           :: Restart
INTEGER           :: m,nn,o
! preconditoner + Vt
#ifdef DLINANALYZE
REAL              :: tS,tE, tS2,tE2,t1,t2
real              :: tstart,tend,tPMV
#endif /* DLINANALYZE */
!===================================================================================================================================

#ifdef DLINANALYZE
! time measurement
CALL CPU_TIME(tS)
! start GMRES
tPMV=0.
#endif /* DLINANALYZE */

Restart=0
nPartInnerIter=0
DeltaX=0.

! ignore particles with zero change
! maybe a large tolerance is feasible, e.g. eps_Mach?
IF(ABS(Norm_B).EQ.0.) RETURN

! select eisenstat-walker
IF (.NOT.EisenstatWalker) THEN
  AbortCritLoc=Norm_B*epsPartlinSolver
ELSE
  AbortCritLoc=Norm_B*AbortCrit
END IF
R0=B
Norm_R0=Norm_B

V(:,1)=R0/Norm_R0
Gam(1)=Norm_R0

DO WHILE (Restart<nRestarts)
  DO m=1,nKDimPart
    nPartInnerIter=nPartInnerIter+1
#ifdef DLINANALYZE
    CALL CPU_TIME(tStart)
#endif /* DLINANALYZE */
    ! matrix vector
    CALL PartMatrixVector(t,coeff,PartID,V(:,m),W)
#ifdef DLINANALYZE
    CALL CPU_TIME(tend)
    tPMV=tPMV+tend-tStart
#endif /* DLINANALYZE */
    ! Gram-Schmidt
    DO nn=1,m
      CALL PartVectorDotProduct(V(:,nn),W,H(nn,m))
      W=W-H(nn,m)*V(:,nn)
    END DO !nn
    CALL PartVectorDotProduct(W,W,Resu)
    H(m+1,m)=SQRT(Resu)
    ! Givens Rotation
    DO nn=1,m-1
      Temp     =   C(nn)*H(nn,m) + S(nn)*H(nn+1,m)
      H(nn+1,m) = - S(nn)*H(nn,m) + C(nn)*H(nn+1,m)
      H(nn,m)   =   Temp
    END DO !nn
    Bet=SQRT(H(m,m)*H(m,m)+H(m+1,m)*H(m+1,m))
    S(m)=H(m+1,m)/Bet
    C(m)=H(m,m)/Bet
    H(m,m)=Bet
    Gam(m+1)=-S(m)*Gam(m)
    Gam(m)=C(m)*Gam(m)
    IF ((ABS(Gam(m+1)).LE.AbortCritloc) .OR. (m.EQ.nKDimPart)) THEN !converge or max Krylov reached
    !IF (m.EQ.nKDimPart) THEN !converge or max Krylov reached
      DO nn=m,1,-1
         Alp(nn)=Gam(nn)
         DO o=nn+1,m
           Alp(nn)=Alp(nn) - H(nn,o)*Alp(o)
         END DO !o
         Alp(nn)=Alp(nn)/H(nn,nn)
      END DO !nn
      DO nn=1,m
        DeltaX=DeltaX+Alp(nn)*V(:,nn)
      END DO !nn
      !IF (ABS(Gam(m+1)).LE.AbortCritloc) THEN !converged
        totalPartIterLinearSolver=totalPartIterLinearSolver+nPartInnerIter
        ! already back transformed,...more storage...but its ok
#ifdef DLINANALYZE
        IF(nPartInnerIter.GT.1)THEN
          IPWRITE(UNIT_stdOut,*) 'nPartInnerIter - in GMRES',nPartInnerIter
        END IF
        CALL CPU_TIME(tE)
        SWRITE(UNIT_stdOut,'(A22,I5)')      ' Part Iter LinSolver: ',nPartInnerIter
        SWRITE(UNIT_stdOut,'(A22,I5)')      ' nRestarts          : ',Restart
        SWRITE(UNIT_stdOut,'(A22,F16.9)')   ' Time in GMRES      : ',tE-tS
        SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R0            : ',Gam(1)
        SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R             : ',Gam(m+1)
#endif /* DLINANALYZE */
        RETURN
      !END IF  ! converged
    ELSE ! no convergence, next iteration   ((ABS(Gam(m+1)).LE.AbortCrit) .OR. (m.EQ.nKDim))
      V(:,m+1)=W/H(m+1,m)
    END IF ! ((ABS(Gam(m+1)).LE.AbortCrit) .OR. (m.EQ.nKDim))
  END DO ! m
  ! Restart needed
#ifdef ROS
CALL abort(&
__STAMP__&
,'No Restart should be required! Computation of wrong RHS! nkDim',nKDimPart)
#endif
  Restart=Restart+1
  ! new settings for source
  !U=DeltaX
! start residuum berrechnen
  CALL PartMatrixVector(t,Coeff,PartID,DeltaX,R0) ! coeff*Ut+Source^n+1 ! only output
  R0=B-R0
  CALL PartVectorDotProduct(R0,R0,Norm_R0)
  Norm_R0=SQRT(Norm_R0)
  ! GMRES(m)  inner loop
  V(:,1)=R0/Norm_R0
  Gam(1)=Norm_R0
END DO ! Restart

IPWRITE(*,*) 'Gam(1+1)',Gam(m),AbortCrit
CALL abort(&
__STAMP__&
,'GMRES_M NOT CONVERGED WITH RESTARTS AND GMRES ITERATIONS:',Restart,REAL(nPartInnerIter))

END SUBROUTINE Particle_GMRES
#endif /*ROS or IMPA*/

#ifdef IMPA
SUBROUTINE Particle_Armijo(t,coeff,AbortTol,nInnerPartNewton)
!===================================================================================================================================
! an intermediate Armijo step to ensure global convergence
! search direction is d = - F'(U)^-1 F(U), e.g. result of Newton-Step
! Step is limited, if no convergence
! See: Algorithm 8.2.1 on p. 130 of: Kelly: Iterative Methods for linear and nonlinear equations
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_LinearOperator         ,ONLY: PartMatrixVector, PartVectorDotProduct
USE MOD_Particle_Vars          ,ONLY: PartState,F_PartXK,Norm_F_PartXK,PartQ,PartLorentzType,DoPartInNewton,PartLambdaAccept &
                                     ,PartDeltaX,PEM,PDM,LastPartPos,Pt,Norm_F_PartX0,PartDtFrac,PartStateN &
                                     ,PartMeshHasReflectiveBCs
USE MOD_LinearSolver_Vars      ,ONLY: PartXK,R_PartXK,DoPrintConvInfo
USE MOD_LinearSolver_Vars      ,ONLY: Part_alpha
USE MOD_Part_RHS               ,ONLY: PartRHS
USE MOD_PICInterpolation       ,ONLY: InterpolateFieldToSingleParticle
USE MOD_PICInterpolation_Vars  ,ONLY: FieldAtParticle
USE MOD_Globals_Vars           ,ONLY: c2_inv
USE MOD_Particle_Tracking_vars ,ONLY: TrackingMethod
USE MOD_Particle_Tracing       ,ONLY: ParticleTracing
USE MOD_Particle_RefTracking   ,ONLY: ParticleRefTracking
USE MOD_Particle_TriaTracking  ,ONLY: ParticleTriaTracking
USE MOD_LinearSolver_Vars      ,ONLY: DoFullNewton!,PartNewtonRelaxation
#if USE_MPI
USE MOD_Particle_MPI           ,ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers     ,ONLY: LBStartTime,LBPauseTime,LBSplitTime
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/
#ifdef CODE_ANALYZE
USE MOD_Particle_Tracking_Vars ,ONLY: PartOut,MPIRankOut
#endif /*CODE_ANALYZE*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
REAL,INTENT(IN)              :: t
REAL,INTENT(IN)              :: coeff
REAL,INTENT(IN)              :: AbortTol
INTEGER,INTENT(IN)           :: nInnerPartNewton
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iPart,iCounter
REAL                         :: lambda, Norm_PartX,DeltaX_Norm
REAL                         :: LorentzFacInv,Xtilde(1:6), DeltaX(1:6)
LOGICAL                      :: DoSetLambda
INTEGER                      :: nLambdaReduce,nMaxLambdaReduce=10
#if USE_MPI
#if USE_LOADBALANCE
REAL                         :: tLBStart
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/
REAL                         :: n_loc(1:3), PartStateTmp(1:6)
LOGICAL                      :: ReMap
!===================================================================================================================================

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
lambda=1. !*PartNewtonRelaxation
DoSetLambda=.TRUE.
PartLambdaAccept=.TRUE.
DO iPart=1,PDM%ParticleVecLength
  IF(DoPartInNewton(iPart))THEN
    ! caution: PartXK has to be used instead of PartState
    LastPartPos(1,iPart)=PartStateN(1,iPart)
    LastPartPos(2,iPart)=PartStateN(2,iPart)
    LastPartPos(3,iPart)=PartStateN(3,iPart)
    PEM%LastGlobalElemID(iPart)=PEM%GlobalElemID(iPart)
    ! and disable periodic movement
    IF(PartMeshHasReflectiveBCs) PEM%NormVec(:,iPart)=0.
    PEM%PeriodicMoved(iPart)=.FALSE.
    ! new part: of Armijo algorithm: check convergence
    ! compute new function value
    CALL PartMatrixVector(t,Coeff,iPart,PartDeltaX(:,iPart),Xtilde) ! coeff*Ut+Source^n+1 ! only output
    XTilde=XTilde+F_PartXK(1:6,iPart)
    CALL PartVectorDotProduct(Xtilde,Xtilde,Norm_PartX)
    Norm_PartX=SQRT(Norm_PartX)
    IF(Norm_PartX.GT.AbortTol*Norm_F_PartXK(iPart))THEN
      ! bad search direction!
      ! new search direction
      DeltaX=PartDeltaX(:,iPart)!*2
      CALL PartMatrixVector(t,Coeff,iPart,DeltaX(:),Xtilde) ! coeff*Ut+Source^n+1 ! only output
      XTilde=XTilde+F_PartXK(1:6,iPart)
      CALL PartVectorDotProduct(Xtilde,Xtilde,Norm_PartX)
      Norm_PartX=SQRT(Norm_PartX)
!      IF(Norm2_PartX.GT.AbortTol*Norm2_F_PartXK(iPart))THEN
!        Norm2_PartX = Norm2_PartX/Norm2_F_PartXk(iPart)
!        IPWRITE(UNIT_stdOut,'(I0,A,6(X,E24.12))') ' found wrong search direction', deltaX
!        IPWRITE(UNIT_stdOut,'(I0,A,6(X,E24.12))') ' found wrong search direction', PartDeltaX(1:6,iPart)
!        CALL abort(&
!  __STAMP__&
!  ,' Found wrond search direction! Particle, Monitored decrease: ', iPart, Norm2_PartX)
!     END IF
    END IF
#ifdef CODE_ANALYZE
    IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
      IF(iPart.EQ.PARTOUT)THEN
        IPWRITE(UNIT_stdOut,'(I0,A,6(X,E24.12))') ' PartDeltaX: ', PartDeltaX(1:6,iPart)
      END IF !(iPart.EQ.PARTOUT)THEN
    END IF !(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
#endif /*CODE_ANALYZE*/
    ! update position
    PartState(1:6,iPart)=PartXK(1:6,iPart)+lambda*PartDeltaX(1:6,iPart)
    PartLambdaAccept(iPart)=.FALSE.
  END IF ! ParticleInside
END DO ! iPart
#if USE_LOADBALANCE
CALL LBSplitTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/

! move particle
#if USE_MPI
! open receive buffer for number of particles
CALL IRecvNbofParticles() ! input value: which list:PartLambdaAccept or PDM%ParticleInisde?
#if USE_LOADBALANCE
CALL LBPauseTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/
SELECT CASE(TrackingMethod)
CASE(REFMAPPING)
  CALL ParticleRefTracking(doParticle_In=.NOT.PartLambdaAccept(1:PDM%ParticleVecLength))
CASE(TRACING)
  CALL ParticleTracing(doParticle_In=.NOT.PartLambdaAccept(1:PDM%ParticleVecLength))
CASE(TRIATRACKING)
  CALL ParticleTriaTracking(doParticle_In=.NOT.PartLambdaAccept(1:PDM%ParticleVecLength))
CASE DEFAULT
  CALL abort(__STAMP__,'TrackingMethod not implemented! TrackingMethod =',IntInfoOpt=TrackingMethod)
END SELECT

DO iPart=1,PDM%ParticleVecLength
  IF(.NOT.PDM%ParticleInside(iPart))THEN
    DoPartInNewton(iPart)=.FALSE.
    PartLambdaAccept(iPart)=.TRUE.
  END IF
  IF(.NOT.PartLambdaAccept(iPart))THEN
    IF(.NOT.PDM%ParticleInside(iPart))THEN
      DoPartInNewton(iPart)=.FALSE.
      PartLambdaAccept(iPart)=.TRUE.
    END IF
  END IF
END DO
#if USE_MPI
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
! send number of particles
CALL SendNbOfParticles(doParticle_In=.NOT.PartLambdaAccept(1:PDM%ParticleVecLength))
! finish communication of number of particles and send particles
CALL MPIParticleSend() ! input value: which list:PartLambdaAccept or PDM%ParticleInisde?
! finish communication
CALL MPIParticleRecv() ! input value: which list:PartLambdaAccept or PDM%ParticleInisde?
#if USE_LOADBALANCE
CALL LBSplitTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif

DO iPart=1,PDM%ParticleVecLength
  IF(.NOT.PartLambdaAccept(iPart))THEN
#if USE_MPI
    IF(.NOT.PDM%ParticleInside(iPart))THEN
      DoPartInNewton(iPart)=.FALSE.
      PartLambdaAccept(iPart)=.TRUE.
      CYCLE
    END IF
#endif /*USE_MPI*/
    ! compute lorentz force at particles position
    CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle(1:6,iPart))
    reMap=.FALSE.
    IF(PartMeshHasReflectiveBCs)THEN
      IF(SUM(ABS(PEM%NormVec(1:3,iPart))).GT.0.)THEN
        n_loc=PEM%NormVec(1:3,iPart)
        ! particle is actually located outside, hence, it moves in the mirror field
        FieldAtParticle(1:3,iPart)=FieldAtParticle(1:3,iPart)-2.*DOT_PRODUCT(FieldAtParticle(1:3,iPart),n_loc)*n_loc
        FieldAtParticle(4:6,iPart)=FieldAtParticle(4:6,iPart)!-2.*DOT_PRODUCT(FieldAtParticle(4:6,iPart),n_loc)*n_loc
        ! reset part state to the not-reflected position
        !PEM%NormVec(1:3,iPart)=0.
        reMap=.TRUE.
      END IF
    END IF
    IF(PEM%PeriodicMoved(iPart)) reMap=.TRUE.
    IF(reMap)THEN
      ! stoare old position within mesh || required for deposition
      PartStateTmp(1:6) = PartState(1:6,iPart)
      PartState(1:6,iPart)=PartXK(1:6,iPart)+lambda*PartDeltaX(1:6,iPart)
    END IF
    IF(PartLorentzType.EQ.5)THEN
      LorentzFacInv=1.0/SQRT(1.0+DOT_PRODUCT(PartState(4:6,iPart),PartState(4:6,iPart))*c2_inv)
      CALL PartRHS(iPart,FieldAtParticle(1:6,iPart),Pt(1:3,iPart),LorentzFacInv)
    ELSE
      LorentzFacInv = 1.0
      CALL PartRHS(iPart,FieldAtParticle(1:6,iPart),Pt(1:3,iPart))
    END IF ! PartLorentzType.EQ.5
    R_PartXK(1:3,iPart)=LorentzFacInv*PartState(4:6,iPart)
    R_PartXK(4,iPart)=Pt(1,iPart)
    R_PartXK(5,iPart)=Pt(2,iPart)
    R_PartXK(6,iPart)=Pt(3,iPart)
    F_PartXK(1:6,iPart)=PartState(1:6,iPart) - PartQ(1:6,iPart) - PartDtFrac(iPart)*coeff*R_PartXK(1:6,iPart)
    ! if check, then here!
    DeltaX_Norm=DOT_PRODUCT(PartDeltaX(1:6,iPart),PartDeltaX(1:6,iPart))
    DeltaX_Norm=SQRT(DeltaX_Norm)
#ifdef CODE_ANALYZE
    IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
      IF(iPart.EQ.PARTOUT)THEN
        IPWRITE(UNIT_stdOut,'(I0,A,G0)') ' DeltaX_Norm', DeltaX_Norm
      END IF !(iPart.EQ.PARTOUT)
    END IF !(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
#endif /*CODE_ANALYZE*/
    IF(DeltaX_Norm.LT.AbortTol*Norm_F_PartX0(iPart)) THEN
       DoPartInNewton(iPart)=.FALSE.
       PartLambdaAccept(iPart)=.TRUE.
       PartXK(1:6,iPart)=PartState(1:6,iPart)
       PartDeltaX(1:6,iPart)=0.
    ELSE
      !IF(nInnerPartNewton.LT.5)THEN
      !  ! accept lambda
      !  PartLambdaAccept(iPart)=.TRUE.
      !  ! set  new position
      !  PartXK(1:6,iPart)=PartState(1:6,iPart)
      !  ! update norm
      !  CALL PartVectorDotProduct(F_PartXK(1:6,iPart),F_PartXK(1:6,iPart),Norm2_PartX)
      !  Norm2_F_PartXK(iPart)=Norm2_PartX
      !  IF((Norm2_F_PartXK(iPart).LT.AbortTol*Norm2_F_PartX0(iPart)).OR.(Norm2_F_PartXK(iPart).LT.1e-12)) &
      !      DoPartInNewton(iPart)=.FALSE.
      !ELSE
        ! check if residual is reduced
        CALL PartVectorDotProduct(F_PartXK(1:6,iPart),F_PartXK(1:6,iPart),Norm_PartX)
        Norm_PartX=SQRT(Norm_PartX)
#ifdef CODE_ANALYZE
        IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
          IF(iPart.EQ.PARTOUT)THEN
            IPWRITE(UNIT_stdOut,'(I0,A,G0,x,G0,x,G0)') ' Norms: ', Norm_PartX, Norm_F_PartXK(iPart), Norm_F_PartX0(iPart)
          END IF !(iPart.EQ.PARTOUT)THEN
        END IF !(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
#endif /*CODE_ANALYZE*/
        IF(Norm_PartX .LT. (1.-Part_alpha*lambda)*Norm_F_PartXK(iPart))THEN
          ! accept lambda
          PartLambdaAccept(iPart)=.TRUE.
          ! set  new position
          PartXK(1:6,iPart)=PartState(1:6,iPart)
          PartDeltaX(1:6,iPart)=0.
          Norm_F_PartXK(iPart)=Norm_PartX
          IF((Norm_F_PartXK(iPart).LT.AbortTol*Norm_F_PartX0(iPart)).OR.(Norm_F_PartXK(iPart).LT.1e-12)) &
              DoPartInNewton(iPart)=.FALSE.
        ELSE
          ! nothing to do, do not accept lambda
        END IF
  !    END IF ! nInnerPartNewton>1
    END IF
    IF(reMap) PartState(1:6,iPart) = PartStateTmp(1:6)
  END IF
END DO ! iPart=1,PDM%ParticleVecLength
#if USE_LOADBALANCE
CALL LBSplitTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/

! disable Armijo iteration and use only one fixed value
!IF(PartNewtonRelaxation.LT.1.)  PartLambdaAccept=.TRUE.

DoSetLambda=.FALSE.
IF(ANY(.NOT.PartLambdaAccept)) DoSetLambda=.TRUE.
#if USE_MPI
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
!set T if at least 1 proc has to do newton
CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoSetLambda,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_PICLAS,iError)
#if USE_LOADBALANCE
CALL LBSplitTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*USE_MPI*/
IF(DoPrintConvInfo)THEN
  SWRITE(UNIT_stdOut,'(A20,2x,L)') ' Lambda-Accept: ', DoSetLambda
END IF

nLambdaReduce=1
DO WHILE((DoSetLambda).AND.(nLambdaReduce.LE.nMaxLambdaReduce))
  nLambdaReduce=nLambdaReduce+1
  lambda=0.5*lambda
  IF(DoPrintConvInfo)THEN
    SWRITE(UNIT_stdOut,'(A20,2x,E24.12)') ' lambda: ', lambda
  END IF
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  DO iPart=1,PDM%ParticleVecLength
    IF(.NOT.PartLambdaAccept(iPart))THEN
#if USE_MPI
    IF(.NOT.PDM%ParticleInside(iPart))THEN
      DoPartInNewton(iPart)=.FALSE.
      PartLambdaAccept(iPart)=.TRUE.
      CYCLE
    END IF
#endif /*USE_MPI*/
      ! update the last part pos and element for particle movement
      LastPartPos(1,iPart)=PartStateN(1,iPart)
      LastPartPos(2,iPart)=PartStateN(2,iPart)
      LastPartPos(3,iPart)=PartStateN(3,iPart)
      PEM%LastGlobalElemID(iPart)=PEM%GlobalElemID(iPart)
      IF(PartMeshHasReflectiveBCs) PEM%NormVec(1:3,iPart)=0.
      PEM%PeriodicMoved(iPart)=.FALSE.
      ! recompute part state
      PartState(1:6,iPart)=PartXK(:,iPart)+lambda*PartDeltaX(:,iPart)
      PartLambdaAccept(iPart)=.FALSE.
    END IF ! ParticleInside
  END DO ! iPart
  ! move particle
#if USE_MPI
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/
  ! open receive buffer for number of particles
  CALL IRecvNbofParticles() ! input value: which list:PartLambdaAccept or PDM%ParticleInisde?
#endif /*USE_MPI*/
  SELECT CASE(TrackingMethod)
  CASE(REFMAPPING)
    CALL ParticleRefTracking(doParticle_In=.NOT.PartLambdaAccept(1:PDM%ParticleVecLength))
  CASE(TRACING)
    CALL ParticleTracing(doParticle_In=.NOT.PartLambdaAccept(1:PDM%ParticleVecLength))
  CASE(TRIATRACKING)
    CALL ParticleTriaTracking(doParticle_In=.NOT.PartLambdaAccept(1:PDM%ParticleVecLength))
  CASE DEFAULT
    CALL abort(__STAMP__,'TrackingMethod not implemented! TrackingMethod =',IntInfoOpt=TrackingMethod)
  END SELECT
  DO iPart=1,PDM%ParticleVecLength
    IF(.NOT.PDM%ParticleInside(iPart))THEN
      DoPartInNewton(iPart)=.FALSE.
      PartLambdaAccept(iPart)=.TRUE.
    END IF
    IF(.NOT.PartLambdaAccept(iPart))THEN
      IF(.NOT.PDM%ParticleInside(iPart))THEN
        DoPartInNewton(iPart)=.FALSE.
        PartLambdaAccept(iPart)=.TRUE.
      END IF
    END IF
  END DO
#if USE_MPI
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  ! send number of particles
  CALL SendNbOfParticles(doParticle_In=.NOT.PartLambdaAccept(1:PDM%ParticleVecLength))
  ! finish communication of number of particles and send particles
  CALL MPIParticleSend() ! input value: which list:PartLambdaAccept or PDM%ParticleInisde?
  ! finish communication
  CALL MPIParticleRecv() ! input value: which list:PartLambdaAccept or PDM%ParticleInisde?
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif

#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  DO iPart=1,PDM%ParticleVecLength
    IF(.NOT.PartLambdaAccept(iPart))THEN
#if USE_MPI
      IF(.NOT.PDM%ParticleInside(iPart))THEN
        DoPartInNewton(iPart)=.FALSE.
        PartLambdaAccept(iPart)=.TRUE.
        CYCLE
      END IF
#endif /*USE_MPI*/
      ! compute lorentz-force at particle's position
      CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle(1:6,iPart))
      reMap=.FALSE.
      IF(PartMeshHasReflectiveBCs)THEN
        IF(SUM(ABS(PEM%NormVec(1:3,iPart))).GT.0.)THEN
          n_loc=PEM%NormVec(1:3,iPart)
          ! particle is actually located outside, hence, it moves in the mirror field
          FieldAtParticle(1:3,iPart)=FieldAtParticle(1:3,iPart)-2.*DOT_PRODUCT(FieldAtParticle(1:3,iPart),n_loc)*n_loc
          FieldAtParticle(4:6,iPart)=FieldAtParticle(4:6,iPart)!-2.*DOT_PRODUCT(FieldAtParticle(4:6,iPart),n_loc)*n_loc
          reMap=.TRUE.
        END IF
      END IF
      IF(PEM%PeriodicMoved(iPart)) reMap=.TRUE.
      IF(reMap)THEn
        PartStateTmp(1:6) = PartState(1:6,iPart)
        PartState(1:6,iPart)=PartXK(1:6,iPart)+lambda*PartDeltaX(1:6,iPart)
      END IF
      IF(PartLorentzType.EQ.5)THEN
        LorentzFacInv=1.0/SQRT(1.0+DOT_PRODUCT(PartState(4:6,iPart),PartState(4:6,iPart))*c2_inv)
        CALL PartRHS(iPart,FieldAtParticle(1:6,iPart),Pt(1:3,iPart),LorentzFacInv)
      ELSE
        LorentzFacInv = 1.0
        CALL PartRHS(iPart,FieldAtParticle(1:6,iPart),Pt(1:3,iPart))
      END IF ! PartLorentzType.EQ.5
      R_PartXK(1:3,iPart)=LorentzFacInv*PartState(4:6,iPart)
      R_PartXK(4,iPart)=Pt(1,iPart)
      R_PartXK(5,iPart)=Pt(2,iPart)
      R_PartXK(6,iPart)=Pt(3,iPart)
      F_PartXK(1:6,iPart)=PartState(1:6,iPart) - PartQ(1:6,iPart) - PartDtFrac(iPart)*coeff*R_PartXK(1:6,iPart)
      ! vector dot product
      CALL PartVectorDotProduct(F_PartXK(:,iPart),F_PartXK(:,iPart),Norm_PartX)
      Norm_PartX=SQRT(Norm_PartX)
      !IF(Norm2_PartX .LT. (1.-Part_alpha*lambda)*Norm2_F_PartXK(iPart))THEN
      IF(DoFullNewton)THEN
        ! accept lambda
        PartLambdaAccept(iPart)=.TRUE.
        ! set  new position
        PartXK(1:6,iPart)=PartState(1:6,iPart)
        PartDeltaX(1:6,iPart)=0.
        Norm_F_PartXK(iPart)=Norm_PartX
        IF((Norm_F_PartXK(iPart).LT.AbortTol*Norm_F_PartX0(iPart)).OR.(Norm_F_PartXK(iPart).LT.1e-12)) &
           DoPartInNewton(iPart)=.FALSE.
      ELSE ! .NOT.DoFullNewton
        IF(Norm_PartX .LE. Norm_F_PartXK(iPart))THEN
          ! accept lambda
          PartLambdaAccept(iPart)=.TRUE.
          ! set  new position
          PartXK(1:6,iPart)=PartState(1:6,iPart)
          PartDeltaX(1:6,iPart)=0.
          Norm_F_PartXK(iPart)=Norm_PartX
          IF((Norm_F_PartXK(iPart).LT.AbortTol*Norm_F_PartX0(iPart)).OR.(Norm_F_PartXK(iPart).LT.1e-12)) &
            DoPartInNewton(iPart)=.FALSE.
        ELSE
          ! test not working
          !IF(Norm2_PartX.GT.Norm2_F_PartX0(iPart))THEN !allow for a local increase in residual
          !  PartXK(1:6,iPart)=PartState(1:6,iPart)
          !  Norm2_F_PartXK(iPart)=Norm2_PartX
          !  Norm2_F_PartX0(iPart)=Norm2_PartX
          !  Norm2_F_PartXk_old(iPart)=Norm2_PartX
          !  PartLambdaAccept(iPart)=.TRUE.
          !END IF
          ! DO not accept lambda, go to next step
          ! scip particle in current iteration and reiterate
          ! DO NOT nullify because Armijo iteration will not work
          ! OR  NO Armijo and particle does is not changed during this step
          PartDeltaX(1:3,iPart)=0.
          DoPartInNewton(iPart)=.FALSE.
          PartLambdaAccept(iPart)=.TRUE.
        END IF
      END IF ! DoFullNewton
      ! remap is performed for deposition
      IF(reMap) PartState(1:6,iPart) = PartStateTmp(1:6)
    END IF
  END DO ! iPart=1,PDM%ParticleVecLength
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_PUSH,tLBStart)
#endif /*USE_LOADBALANCE*/
  ! detect  convergence
  DoSetLambda=.FALSE.
  IF(ANY(.NOT.PartLambdaAccept)) DoSetLambda=.TRUE.
#if USE_MPI
  !set T if at least 1 proc has to do newton
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoSetLambda,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/
  iCounter=0
  DO iPart=1,PDM%ParticleVecLength
    IF(.NOT.PartLambdaAccept(iPart))THEN
      PartLambdaAccept(iPart)=.FALSE.
      iCounter=iCounter+1
    END IF ! ParticleInside
  END DO ! iPart
  IF(DoPrintConvInfo)THEN
#if USE_MPI
    !set T if at least 1 proc has to do newton
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,iCounter,1,MPI_INTEGER,MPI_SUM,MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/
    SWRITE(UNIT_stdOut,'(A20,2x,L,2x,I10)') ' Accept?: ', DoSetLambda,iCounter
  END IF
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_PARTCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
END DO

IF(1.EQ.2)THEN
  iPart=nInnerPartNewton
END IF

END SUBROUTINE Particle_Armijo
#endif /*IMPA*/

#if defined(IMPA) || defined(ROS)
SUBROUTINE FinalizePartSolver()
!===================================================================================================================================
! deallocate global variables
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
USE MOD_LinearSolver_Vars
#ifdef IMPA
USE MOD_Particle_Vars,           ONLY:F_PartX0,F_PartXk,Norm_F_PartX0,Norm_F_PartXK,Norm_F_PartXK_old,DoPartInNewton &
                                     ,PartDeltaX,PartLambdaAccept
#endif /*IMPA*/
USE MOD_Particle_Vars,           ONLY:PartQ

!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE(PartXK)
SDEALLOCATE(R_PartXK)
! variables of particle_vars.f90
SDEALLOCATE(PartQ)
#ifdef IMPA
SDEALLOCATE(F_PartX0)
SDEALLOCATE(F_PartXk)
SDEALLOCATE(PartLambdaAccept)
SDEALLOCATE(PartDeltaX)
SDEALLOCATE(Norm_F_PartX0)
SDEALLOCATE(Norm_F_PartXK)
SDEALLOCATE(Norm_F_PartXK_old)
SDEALLOCATE(DoPartInNewton)
#endif /*IMPA*/
END SUBROUTINE FinalizePartSolver
#endif /*IMPA or ROS*/
#endif /*PARTICLES*/

END MODULE MOD_ParticleSolver
