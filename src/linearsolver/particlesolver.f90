#include "boltzplatz.h"

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

#if defined(PARTICLES) && defined(IMPA)
INTERFACE ParticleNewton
  MODULE PROCEDURE ParticleNewton
END INTERFACE

INTERFACE InitPartSolver
  MODULE PROCEDURE InitPartSolver
END INTERFACE

#if (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122)
INTERFACE SelectImplicitParticles
  MODULE PROCEDURE SelectImplicitParticles
END INTERFACE
#endif

PUBLIC:: InitPartSolver
PUBLIC:: ParticleNewton
#if (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122)
PUBLIC:: SelectImplicitParticles
#endif
#endif /*PARTICLES*/
!===================================================================================================================================

CONTAINS

#if defined(PARTICLES) && defined(IMPA)
SUBROUTINE InitPartSolver() 
!===================================================================================================================================
! read in and allocation of required global variables for implicit particle treatment
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_PreProc
USE MOD_ReadInTools,          ONLY:GETINT,GETREAL,GETLOGICAL
USE MOD_Particle_Vars,        ONLY:PDM
USE MOD_LinearSolver_Vars
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

SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE SOLVER...'

Eps2PartNewton     =GETREAL('EpsPartNewton','0.001')
Eps2PartNewton     =Eps2PartNewton**2
nPartNewtonIter    =GETINT('nPartNewtonIter','20')
FreezePartInNewton =GETINT('FreezePartInNewton','1')
EisenstatWalker    =GETLOGICAL('EisenstatWalker','.FALSE.')
PartgammaEW        =GETREAL('PartgammaEW','0.9')
nPartNewton        =0

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

END SUBROUTINE InitPartSolver


#if (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122)
SUBROUTINE SelectImplicitParticles() 
!===================================================================================================================================
! select if particle is treated implicitly or explicitly, has to be called, after particle are created/emitted
! currently only one criterion is used: the species
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Vars,     ONLY:Species,nSpecies,PartSpecies,PartIsImplicit,PDM
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER     :: iPart
!===================================================================================================================================

PartIsImplicit=.FALSE.
DO iPart=1,PDM%ParticleVecLength
  IF(.NOT.PDM%ParticleInside(iPart)) CYCLE
  IF(Species(PartSpecies(iPart))%IsImplicit) PartIsImplicit(iPart)=.TRUE.
END DO ! iPart
  
END SUBROUTINE SelectImplicitParticles
#endif


SUBROUTINE ParticleNewton(t,coeff,doParticle_In,opt_In)
!===================================================================================================================================
! Allocate global variable 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_LinearSolver_Vars,       ONLY:PartXK,R_PartXK
USE MOD_Particle_Vars,           ONLY:PartQ,F_PartX0,F_PartXk,Norm2_F_PartX0,Norm2_F_PartXK,Norm2_F_PartXK_old,DoPartInNewton
USE MOD_Particle_Vars,           ONLY:PartState, Pt, LastPartPos, DelayTime, PEM, PDM, usevMPF,PartLorentzType
USE MOD_TimeDisc_Vars,           ONLY:dt,iter
USE MOD_PICInterpolation,        ONLY:InterpolateFieldToParticle
USE MOD_LinearOperator,          ONLY:PartVectorDotProduct
USE MOD_Particle_Tracking,       ONLY:ParticleTrackingCurved,ParticleRefTracking
USE MOD_Particle_Tracking_vars,  ONLY:DoRefMapping
USE MOD_Part_RHS,                ONLY:CalcPartRHS
#ifdef MPI
USE MOD_Particle_MPI,            ONLY:IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY:PartMPIExchange,PartMPI
USE MOD_Particle_MPI_Vars,      ONLY:ExtPartState,ExtPartSpecies,ExtPartMPF,ExtPartToFIBGM,NbrOfExtParticles
#endif /*MPI*/
USE MOD_LinearSolver_vars,       ONLY:Eps2PartNewton,nPartNewton, PartgammaEW,nPartNewtonIter,FreezePartInNewton
USE MOD_Part_RHS,                ONLY:SLOW_RELATIVISTIC_PUSH,FAST_RELATIVISTIC_PUSH
USE MOD_PICInterpolation,        ONLY:InterpolateFieldToSingleParticle
USE MOD_PICInterpolation_Vars,   ONLY:FieldAtParticle
#if (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122)
USE MOD_Particle_Vars,           ONLY:PartIsImplicit
#endif
!USE MOD_Equation,       ONLY: CalcImplicitSource
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: t,coeff
LOGICAL,INTENT(IN),OPTIONAL   :: doParticle_In(1:PDM%maxParticleNumber)
LOGICAL,INTENT(IN),OPTIONAL   :: opt_In
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
LOGICAL                      :: opt
REAL                         :: time
INTEGER                      :: iPart
INTEGER                      :: nInnerPartNewton = 0
REAL                         :: AbortCritLinSolver,gammaA,gammaB
!REAL                         :: FieldAtParticle(1:6)
REAL                         :: DeltaX(1:6)
REAL                         :: Pt_tmp(1:6)
LOGICAL                      :: doParticle(1:PDM%maxParticleNumber)
!! maybeeee
!! and thats maybe local??? || global, has to be set false during communication
LOGICAL                      :: DoNewton
INTEGER:: counter
!===================================================================================================================================

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

DoNewton=.FALSE.
IF(ANY(DoPartInNewton)) DoNewton=.TRUE.
#ifdef MPI
!set T if at least 1 proc has to do newton
CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoNewton,1,MPI_LOGICAL,MPI_LOR,PartMPI%COMM,iError)
#endif /*MPI*/

IF(opt)THEN ! compute zero state
  ! whole pt array
  DO iPart=1,PDM%ParticleVecLength
    IF(DoPartInNewton(iPart))THEN
      CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle(iPart,1:6))
      SELECT CASE(PartLorentzType)
      CASE(1)
        Pt(iPart,1:3) = SLOW_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
      CASE(3)
        Pt(iPart,1:3) = FAST_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
      CASE DEFAULT
      END SELECT
      ! PartStateN has to be exchanged by PartQ
      Pt_tmp(1) = PartState(iPart,4) 
      Pt_tmp(2) = PartState(iPart,5) 
      Pt_tmp(3) = PartState(iPart,6) 
      Pt_tmp(4) = Pt(iPart,1) 
      Pt_tmp(5) = Pt(iPart,2) 
      Pt_tmp(6) = Pt(iPart,3)
      F_PartX0(1:6,iPart) =   PartState(iPart,1:6)-PartQ(1:6,iPart)-coeff*Pt_tmp
      PartXK(1:6,iPart)   =   PartState(iPart,1:6)
      R_PartXK(1:6,iPart) =   Pt_tmp(1:6)
      F_PartXK(1:6,iPart) =   F_PartX0(1:6,iPart)
    END IF ! ParticleInside
  END DO ! iPart
  
  ! compute norm for each particle
  DO iPart=1,PDM%ParticleVecLength
    IF(DoPartInNewton(iPart))THEN
      CALL PartVectorDotProduct(F_PartX0(:,iPart),F_PartX0(:,iPart),Norm2_F_PartX0(iPart))
      IF (Norm2_F_PartX0(iPart).LT.6E-16) THEN ! do not iterate, as U is already the implicit solution
        Norm2_F_PartXk(iPart)=TINY(1.)
        DoPartInNewton(iPart)=.FALSE.
      ELSE ! we need iterations
        Norm2_F_PartXk(iPart)=Norm2_F_PartX0(iPart)
      END IF
    END IF ! ParticleInside
  END DO ! iPart
END IF

! newton per particle 
Counter=0
DO iPart=1,PDM%ParticleVecLength
  IF(DoPartInNewton(iPart))THEN
    Counter=Counter+1      
  END IF ! ParticleInside
END DO ! iPart
#ifdef MPI
!set T if at least 1 proc has to do newton
CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoNewton,1,MPI_LOGICAL,MPI_LOR,PartMPI%COMM,iError) 
CALL MPI_ALLREDUCE(MPI_IN_PLACE,Counter,1,MPI_INTEGER,MPI_SUM,PartMPI%COMM,iError) 
#endif /*MPI*/

SWRITE(*,*) 'init part',Counter

nInnerPartNewton=-1
DO WHILE((DoNewton) .AND. (nInnerPartNewton.LT.nPartNewtonIter))  ! maybe change loops, finish particle after particle?
  nInnerPartNewton=nInnerPartNewton+1
  !SWRITE(*,*) 'PartNewton',nInnerPartNewton,Counter
  DO iPart=1,PDM%ParticleVecLength
    IF(DoPartInNewton(iPart))THEN
      ! set abort crit      
      IF (nInnerPartNewton.EQ.0) THEN
        AbortCritLinSolver=0.999
      ELSE
        gammaA = PartgammaEW*(Norm2_F_PartXk(iPart))/(Norm2_F_PartXk_old(iPart))
        IF (PartgammaEW*AbortCritLinSolver*AbortCritLinSolver < 0.1) THEN
          gammaB = min(0.999,gammaA)
        ELSE
          gammaB = min(0.999, max(gammaA,PartgammaEW*AbortCritLinSolver*AbortCritLinSolver))
        ENDIF
        AbortCritLinSolver = min(0.999,max(gammaB,0.5*SQRT(Eps2PartNewton)/SQRT(Norm2_F_PartXk(iPart))))
      END IF 
      Norm2_F_PartXk_old(iPart)=Norm2_F_PartXk(iPart)
      CALL Particle_GMRES(t,coeff,iPart,-F_PartXK(:,iPart),SQRT(Norm2_F_PartXk(iPart)),AbortCritLinSolver,DeltaX)
      ! update to new partstate during Newton iteration
      PartXK(:,iPart)=PartXK(:,iPart)+DeltaX
      PartState(iPart,:)=PartXK(:,iPart)
    END IF ! ParticleInside
  END DO ! iPart
  ! closed form: now move particles
  ! further improvement: add flag for DoPartInNewton, if it has to be considered in tracking or communication
  IF(MOD(nInnerPartNewton,FreezePartInNewton).EQ.0)THEN
#ifdef MPI
    ! open receive buffer for number of particles
    CALL IRecvNbofParticles() ! input value: which list:DoPartInNewton or PDM%ParticleInisde?
#endif /*MPI*/
    IF(DoRefMapping)THEN
      ! input value: which list:DoPartInNewton or PDM%ParticleInisde?
      CALL ParticleRefTracking(doParticle_In=DoPartInNewton(1:PDM%ParticleVecLength)) 
    ELSE
      ! input value: which list:DoPartInNewton or PDM%ParticleInisde?
      CALL ParticleTrackingCurved(doParticle_In=DoPartInNewton(1:PDM%ParticleVecLength)) 
    END IF
#ifdef MPI
    ! send number of particles
    CALL SendNbOfParticles(doParticle_In=DoPartInNewton(1:PDM%ParticleVecLength)) 
    ! finish communication of number of particles and send particles
    CALL MPIParticleSend() ! input value: which list:DoPartInNewton or PDM%ParticleInisde?
    ! finish communication
    CALL MPIParticleRecv() ! input value: which list:DoPartInNewton or PDM%ParticleInisde?
    ! as we do not have the shape function here, we have to deallocate something
    SDEALLOCATE(ExtPartState)
    SDEALLOCATE(ExtPartSpecies)
    SDEALLOCATE(ExtPartToFIBGM)
    SDEALLOCATE(ExtPartMPF)
    ! update the last part pos and element for particle movement
    LastPartPos(iPart,1)=PartState(iPart,1)
    LastPartPos(iPart,2)=PartState(iPart,2)
    LastPartPos(iPart,3)=PartState(iPart,3)
    PEM%lastElement(iPart)=PEM%Element(iPart)
!#if (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122)
!    DO iPart=1,PDM%ParticleVecLength
!      IF(.NOT.PartIsImplicit(iPart)) DoPartInNewton(iPart)=.FALSE.
!    END DO
!#endif
#endif
  END IF

  DO iPart=1,PDM%ParticleVecLength
    IF(DoPartInNewton(iPart))THEN
      IF(MOD(nInnerPartNewton,FreezePartInNewton).EQ.0) CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle(iPart,1:6))
      SELECT CASE(PartLorentzType)
      CASE(1)
        Pt(iPart,1:3) = SLOW_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
      CASE(3)
        Pt(iPart,1:3) = FAST_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
      CASE DEFAULT
      CALL abort(&
__STAMP__&
,' Given PartLorentzType does not exist!',PartLorentzType)
      END SELECT
      R_PartXK(1,iPart)=PartState(iPart,4)
      R_PartXK(2,iPart)=PartState(iPart,5)
      R_PartXK(3,iPart)=PartState(iPart,6)
      R_PartXK(4,iPart)=Pt(iPart,1)
      R_PartXK(5,iPart)=Pt(iPart,2)
      R_PartXK(6,iPart)=Pt(iPart,3)
      F_PartXK(:,iPart)=PartState(iPart,:) - PartQ(:,iPart) - coeff*R_PartXK(:,iPart)
      ! vector dot product 
      CALL PartVectorDotProduct(F_PartXK(:,iPart),F_PartXK(:,iPart),Norm2_F_PartXK(iPart))
      IF(Norm2_F_PartXK(iPart).LT.Eps2PartNewton*Norm2_F_PartX0(iPart)) DoPartInNewton(iPart)=.FALSE.
      !IF(nInnerPartNewton.GT.20)THEN
      !  IPWRITE(*,*) 'blubb',iPart, Norm2_F_PartXK(iPart),Norm2_F_PartX0(iPart)
      !END IF
    END IF
  END DO
  DoNewton=.FALSE.
  IF(ANY(DoPartInNewton)) DoNewton=.TRUE.
  Counter=0
  DO iPart=1,PDM%ParticleVecLength
    IF(DoPartInNewton(iPart))THEN
      Counter=Counter+1      
    END IF ! ParticleInside
  END DO ! iPart
#ifdef MPI
  !set T if at least 1 proc has to do newton
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoNewton,1,MPI_LOGICAL,MPI_LOR,PartMPI%COMM,iError) 
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,Counter,1,MPI_INTEGER,MPI_SUM,PartMPI%COMM,iError) 
#endif /*MPI*/
END DO

SWRITE(*,*) 'PartNewton',nInnerPartNewton,Counter
nPartNewton=nPartNewton+nInnerPartNewton
!IF (nInnerPartNewton.EQ.nPartNewtonIter) THEN
!  IF(PartMPI%MPIRoot)THEN
!  CALL abort(&
!__STAMP__&
!,'NEWTON NOT CONVERGED WITH NEWTON ITERATIONS,EPISLON',nInnerPartNewton,Eps2PartNewton)
!  END IF
!END IF

END SUBROUTINE ParticleNewton

SUBROUTINE Particle_GMRES(t,coeff,PartID,B,Norm_B,AbortCrit,DeltaX)
!===================================================================================================================================
! Uses matrix free to solve the linear system
! Attention: We use DeltaX=0 as our initial guess   ! why not Un??
!            X0 is allready stored in U
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_LinearSolver_Vars,    ONLY: eps_LinearSolver,TotalPartIterLinearSolver
USE MOD_LinearSolver_Vars,    ONLY: nKDim,nRestarts,nPartInnerIter,EisenstatWalker
USE MOD_LinearOperator,       ONLY: PartMatrixVector, PartVectorDotProduct
USE MOD_TimeDisc_Vars,        ONLY: dt,iter
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: t,coeff,Norm_B
REAL,INTENT(IN)   :: B(1:6)
REAL,INTENT(INOUT):: AbortCrit
REAL,INTENT(OUT)  :: DeltaX(1:6)
INTEGER,INTENT(IN):: PartID
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: V(1:6,1:nKDim)
REAL              :: V2P(1:6)
REAL              :: W(1:6)
REAL              :: Z(1:6)
REAL              :: R0(1:6)
REAL              :: Gam(1:nKDim+1),C(1:nKDim),S(1:nKDim),H(1:nKDim+1,1:nKDim+1),Alp(1:nKDim+1)
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
!Un(:)=PartState(PartID,:)
IF(iter.EQ.0) THEN
  AbortCrit=eps_LinearSolver
ELSE
  IF (EisenstatWalker .eqv. .FALSE.) THEN
    AbortCrit=eps_LinearSolver
  END IF
END IF
AbortCrit=Norm_B*AbortCrit
R0=B
Norm_R0=Norm_B
DeltaX=0.

V(:,1)=R0/Norm_R0
Gam(1)=Norm_R0

DO WHILE (Restart<nRestarts)
  DO m=1,nKDim
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
    IF ((ABS(Gam(m+1)).LE.AbortCrit) .OR. (m.EQ.nKDim)) THEN !converge or max Krylov reached
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
      IF (ABS(Gam(m+1)).LE.AbortCrit) THEN !converged
        totalPartIterLinearSolver=totalPartIterLinearSolver+nPartInnerIter
        ! already back transformed,...more storage...but its ok
#ifdef DLINANALYZE
        CALL CPU_TIME(tE)
        SWRITE(UNIT_stdOut,'(A22,I5)')      ' Part Iter LinSolver: ',nPartInnerIter
        SWRITE(UNIT_stdOut,'(A22,I5)')      ' nRestarts          : ',Restart
        SWRITE(UNIT_stdOut,'(A22,F16.9)')   ' Time in GMRES      : ',tE-tS
        SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R0            : ',Gam(1)
        SWRITE(UNIT_stdOut,'(A22,E16.8)')   ' Norm_R             : ',Gam(m+1)
#endif /* DLINANALYZE */
        RETURN
      END IF  ! converged
    ELSE ! no convergence, next iteration   ((ABS(Gam(m+1)).LE.AbortCrit) .OR. (m.EQ.nKDim)) 
      V(:,m+1)=W/H(m+1,m)
    END IF ! ((ABS(Gam(m+1)).LE.AbortCrit) .OR. (m.EQ.nKDim))
  END DO ! m 
  ! Restart needed
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

SUBROUTINE FinalizePartSolver() 
!===================================================================================================================================
! deallocate global variables
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
USE MOD_LinearSolver_Vars
USE MOD_Particle_Vars,           ONLY:PartQ,F_PartX0,F_PartXk,Norm2_F_PartX0,Norm2_F_PartXK,Norm2_F_PartXK_old,DoPartInNewton
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
SDEALLOCATE(F_PartX0)
SDEALLOCATE(F_PartXk)
SDEALLOCATE(Norm2_F_PartX0)
SDEALLOCATE(Norm2_F_PartXK)
SDEALLOCATE(Norm2_F_PartXK_old)
SDEALLOCATE(DoPartInNewton)
END SUBROUTINE FinalizePartSolver
#endif /*PARTICLES*/

END MODULE MOD_ParticleSolver
