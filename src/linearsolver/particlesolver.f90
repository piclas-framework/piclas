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

INTERFACE FinalizePartSolver
  MODULE PROCEDURE FinalizePartSolver
END INTERFACE

#if (PP_TimeDiscMethod==120) || (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122)
INTERFACE SelectImplicitParticles
  MODULE PROCEDURE SelectImplicitParticles
END INTERFACE
#endif

PUBLIC:: InitPartSolver,FinalizePartSolver
PUBLIC:: ParticleNewton
#if (PP_TimeDiscMethod==120) || (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122)
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
EpsPartLinSolver   =GETREAL('EpsPartLinSolver','0.0')
IF(EpsPartLinSolver.EQ.0.) EpsPartLinSolver=Eps_LinearSolver
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

#if (PP_TimeDiscMethod==120) || (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122)
PartImplicitMethod =GETINT('Part-ImplicitMethod','0')
#endif

END SUBROUTINE InitPartSolver


#if (PP_TimeDiscMethod==120) || (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122)
SUBROUTINE SelectImplicitParticles() 
!===================================================================================================================================
! select if particle is treated implicitly or explicitly, has to be called, after particle are created/emitted
! currently only one criterion is used: the species
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Vars,     ONLY:Species,PartSpecies,PartIsImplicit,PDM,Pt,PartState
USE MOD_Linearsolver_Vars, ONLY:PartImplicitMethod
USE MOD_TimeDisc_Vars,     ONLY:dt,nRKStages,iter,time
USE MOD_Equation_Vars,     ONLY:c2_inv
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
      NewVelo=PartState(iPart,4:6)+dt/REAL(nRKStages-1)*Pt(iPart,1:3)
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
      NewVelo=PartState(iPart,4:6)
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
  
END SUBROUTINE SelectImplicitParticles
#endif


SUBROUTINE ParticleNewton(t,coeff,doParticle_In,opt_In,AbortTol_In)
!===================================================================================================================================
! Allocate global variable 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_LinearSolver_Vars,       ONLY:PartXK,R_PartXK
USE MOD_Particle_Vars,           ONLY:PartQ,F_PartX0,F_PartXk,Norm2_F_PartX0,Norm2_F_PartXK,Norm2_F_PartXK_old,DoPartInNewton
USE MOD_Particle_Vars,           ONLY:PartState, Pt, LastPartPos, StagePartPos,PEM, PDM, PartLorentzType,PartDeltaX
USE MOD_PICInterpolation,        ONLY:InterpolateFieldToParticle
USE MOD_LinearOperator,          ONLY:PartVectorDotProduct
USE MOD_Particle_Tracking,       ONLY:ParticleTracing,ParticleRefTracking
USE MOD_Particle_Tracking_vars,  ONLY:DoRefMapping
USE MOD_Part_RHS,                ONLY:CalcPartRHS
#ifdef MPI
USE MOD_Particle_MPI,            ONLY:IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY:PartMPI
USE MOD_Particle_MPI_Vars,      ONLY:ExtPartState,ExtPartSpecies,ExtPartMPF,ExtPartToFIBGM,NbrOfExtParticles
#endif /*MPI*/
USE MOD_LinearSolver_vars,       ONLY:Eps2PartNewton,nPartNewton, PartgammaEW,nPartNewtonIter,FreezePartInNewton,DoPrintConvInfo
USE MOD_Part_RHS,                ONLY:SLOW_RELATIVISTIC_PUSH,FAST_RELATIVISTIC_PUSH &
                                     ,RELATIVISTIC_PUSH,NON_RELATIVISTIC_PUSH
USE MOD_Equation_vars,           ONLY:c2_inv
USE MOD_PICInterpolation,        ONLY:InterpolateFieldToSingleParticle
USE MOD_PICInterpolation_Vars,   ONLY:FieldAtParticle
!USE MOD_Equation,       ONLY: CalcImplicitSource
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: t,coeff
LOGICAL,INTENT(INOUT),OPTIONAL:: doParticle_In(1:PDM%maxParticleNumber)
LOGICAL,INTENT(IN),OPTIONAL   :: opt_In
REAL,INTENT(IN),OPTIONAL      :: AbortTol_In
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
!REAL                         :: DeltaX(1:6), DeltaX_Norm
REAL                         :: Pt_tmp(1:6)
!! maybeeee
!! and thats maybe local??? || global, has to be set false during communication
LOGICAL                      :: DoNewton
REAL                         :: AbortTol
REAL                         :: LorentzFacInv
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

IF(PRESENT(AbortTol_In))THEN
  AbortTol=AbortTol_In
ELSE
  AbortTol=Eps2PartNewton
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
      ! update the last part pos and element for particle movement
      LastPartPos(iPart,1)=StagePartPos(iPart,1)
      LastPartPos(iPart,2)=StagePartPos(iPart,2)
      LastPartPos(iPart,3)=StagePartPos(iPart,3)
      PEM%lastElement(iPart)=PEM%StageElement(iPart)
      !! update the last part pos and element for particle movement
      !LastPartPos(iPart,1)=PartState(iPart,1)
      !LastPartPos(iPart,2)=PartState(iPart,2)
      !LastPartPos(iPart,3)=PartState(iPart,3)
      !PEM%lastElement(iPart)=PEM%Element(iPart)
      CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle(iPart,1:6))
      SELECT CASE(PartLorentzType)
      CASE(0)
        Pt(iPart,1:3) = NON_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
        LorentzFacInv = 1.0
      CASE(1)
        Pt(iPart,1:3) = SLOW_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
        LorentzFacInv = 1.0
      CASE(3)
        Pt(iPart,1:3) = FAST_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
        LorentzFacInv = 1.0
      CASE(5)
        LorentzFacInv=1.0+DOT_PRODUCT(PartState(iPart,4:6),PartState(iPart,4:6))*c2_inv      
        LorentzFacInv=1.0/SQRT(LorentzFacInv)
        Pt(iPart,1:3) = RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6),LorentzFacInvIn=LorentzFacInv)
      CASE DEFAULT
      END SELECT
      ! PartStateN has to be exchanged by PartQ
      Pt_tmp(1) = LorentzFacInv*PartState(iPart,4) 
      Pt_tmp(2) = LorentzFacInv*PartState(iPart,5) 
      Pt_tmp(3) = LorentzFacInv*PartState(iPart,6) 
      Pt_tmp(4) = Pt(iPart,1) 
      Pt_tmp(5) = Pt(iPart,2) 
      Pt_tmp(6) = Pt(iPart,3)
      F_PartX0(1:6,iPart) =   PartState(iPart,1:6)-PartQ(1:6,iPart)-coeff*Pt_tmp
      PartXK(1:6,iPart)   =   PartState(iPart,1:6)
      R_PartXK(1:6,iPart) =   Pt_tmp(1:6)
      F_PartXK(1:6,iPart) =   F_PartX0(1:6,iPart)
      CALL PartVectorDotProduct(F_PartX0(:,iPart),F_PartX0(:,iPart),Norm2_F_PartX0(iPart))
      IF (Norm2_F_PartX0(iPart).LT.6E-16) THEN ! do not iterate, as U is already the implicit solution
        Norm2_F_PartXk(iPart)=TINY(1.)
        DoPartInNewton(iPart)=.FALSE.
      ELSE ! we need iterations
        Norm2_F_PartXk(iPart)=Norm2_F_PartX0(iPart)
      END IF
    END IF ! ParticleInside
  END DO ! iPart
ELSE
  DO iPart=1,PDM%ParticleVecLength
    IF(DoPartInNewton(iPart))THEN
      ! update the last part pos and element for particle movement
      LastPartPos(iPart,1)=StagePartPos(iPart,1)
      LastPartPos(iPart,2)=StagePartPos(iPart,2)
      LastPartPos(iPart,3)=StagePartPos(iPart,3)
      PEM%lastElement(iPart)=PEM%StageElement(iPart)
      !LastPartPos(iPart,1)=PartState(iPart,1)
      !LastPartPos(iPart,2)=PartState(iPart,2)
      !LastPartPos(iPart,3)=PartState(iPart,3)
      !PEM%lastElement(iPart)=PEM%Element(iPart)
    END IF ! ParticleInside
  END DO ! iPart
END IF

#ifdef MPI
!set T if at least 1 proc has to do newton
CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoNewton,1,MPI_LOGICAL,MPI_LOR,PartMPI%COMM,iError) 
#endif /*MPI*/

IF(DoPrintConvInfo)THEN
  ! newton per particle 
  Counter=0
  DO iPart=1,PDM%ParticleVecLength
    IF(DoPartInNewton(iPart))THEN
      Counter=Counter+1      
    END IF ! ParticleInside
  END DO ! iPart
#ifdef MPI
  !set T if at least 1 proc has to do newton
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,Counter,1,MPI_INTEGER,MPI_SUM,PartMPI%COMM,iError) 
#endif /*MPI*/
  SWRITE(*,*) ' Initial particle number in newton: ',Counter
END IF

AbortCritLinSolver=0.999
nInnerPartNewton=0
DO WHILE((DoNewton) .AND. (nInnerPartNewton.LT.nPartNewtonIter))  ! maybe change loops, finish particle after particle?
  nInnerPartNewton=nInnerPartNewton+1
  !SWRITE(*,*) 'PartNewton',nInnerPartNewton,Counter
  DO iPart=1,PDM%ParticleVecLength
    IF(DoPartInNewton(iPart))THEN
      ! set abort crit      
      IF (nInnerPartNewton.EQ.1) THEN
        AbortCritLinSolver=0.999
      ELSE
        gammaA = PartgammaEW*(Norm2_F_PartXk(iPart))/(Norm2_F_PartXk_old(iPart))
        IF (PartgammaEW*AbortCritLinSolver*AbortCritLinSolver < 0.1) THEN
          gammaB = MIN(0.999,gammaA)
        ELSE
          gammaB = MIN(0.999, MAX(gammaA,PartgammaEW*AbortCritLinSolver*AbortCritLinSolver))
        ENDIF
        AbortCritLinSolver = MIN(0.999,MAX(gammaB,0.5*SQRT(AbortTol)/SQRT(Norm2_F_PartXk(iPart))))
      END IF 
      Norm2_F_PartXk_old(iPart)=Norm2_F_PartXk(iPart)
      CALL Particle_GMRES(t,coeff,iPart,-F_PartXK(:,iPart),SQRT(Norm2_F_PartXk(iPart)),AbortCritLinSolver,PartDeltaX(1:6,iPart))
      ! additional Armijo step for global convegence
      ! update to new partstate during Newton iteration
!      PartXK(:,iPart)=PartXK(:,iPart)+DeltaX
!      PartState(iPart,:)=PartXK(:,iPart)
!      ! forbidden, because particle is NOT moved but has to be traced...
!      DeltaX_Norm=DOT_PRODUCT(DeltaX,DeltaX)
!      IF(DeltaX_Norm.LT.AbortTol*Norm2_F_PartX0(iPart)) THEN
!        DoPartInNewton(iPart)=.FALSE.
!      END IF
    END IF ! ParticleInside
  END DO ! iPart

  ! DeltaX is going to be global
  CALL Particle_Armijo(t,coeff,AbortTol,nInnerPartNewton) 

  ! already done in particle armijo

!  ! closed form: now move particles
!  ! further improvement: add flag for DoPartInNewton, if it has to be considered in tracking or communication
!  IF(MOD(nInnerPartNewton,FreezePartInNewton).EQ.0)THEN
!#ifdef MPI
!    ! open receive buffer for number of particles
!    CALL IRecvNbofParticles() ! input value: which list:DoPartInNewton or PDM%ParticleInisde?
!#endif /*MPI*/
!    IF(DoRefMapping)THEN
!      ! input value: which list:DoPartInNewton or PDM%ParticleInisde?
!      CALL ParticleRefTracking(doParticle_In=DoPartInNewton(1:PDM%ParticleVecLength)) 
!    ELSE
!      ! input value: which list:DoPartInNewton or PDM%ParticleInisde?
!      CALL ParticleTracing(doParticle_In=DoPartInNewton(1:PDM%ParticleVecLength)) 
!    END IF
!    DO iPart=1,PDM%ParticleVecLength
!      IF(DoPartInNewton(iPart))THEN
!        IF(.NOT.PDM%ParticleInside(iPart))THEN
!          DoPartInNewton(iPart)=.FALSE.
!        END IF
!      END IF
!    END DO
!#ifdef MPI
!    ! send number of particles
!    CALL SendNbOfParticles(doParticle_In=DoPartInNewton(1:PDM%ParticleVecLength)) 
!    ! finish communication of number of particles and send particles
!    CALL MPIParticleSend() ! input value: which list:DoPartInNewton or PDM%ParticleInisde?
!    ! finish communication
!    CALL MPIParticleRecv() ! input value: which list:DoPartInNewton or PDM%ParticleInisde?
!    ! as we do not have the shape function here, we have to deallocate something
!    SDEALLOCATE(ExtPartState)
!    SDEALLOCATE(ExtPartSpecies)
!    SDEALLOCATE(ExtPartToFIBGM)
!    SDEALLOCATE(ExtPartMPF)
!    NbrOfExtParticles=0
!!#if (PP_TimeDiscMethod==121) || (PP_TimeDiscMethod==122)
!!    DO iPart=1,PDM%ParticleVecLength
!!      IF(.NOT.PartIsImplicit(iPart)) DoPartInNewton(iPart)=.FALSE.
!!    END DO
!!#endif
!#endif
!  END IF

!  DO iPart=1,PDM%ParticleVecLength
!    IF(DoPartInNewton(iPart))THEN
!      ! update the last part pos and element for particle movement
!      LastPartPos(iPart,1)=StagePartPos(iPart,1)
!      LastPartPos(iPart,2)=StagePartPos(iPart,2)
!      LastPartPos(iPart,3)=StagePartPos(iPart,3)
!      PEM%lastElement(iPart)=PEM%StageElement(iPart)
!      !LastPartPos(iPart,1)=PartState(iPart,1)
!      !LastPartPos(iPart,2)=PartState(iPart,2)
!      !LastPartPos(iPart,3)=PartState(iPart,3)
!      !PEM%lastElement(iPart)=PEM%Element(iPart)
!      IF(MOD(nInnerPartNewton,FreezePartInNewton).EQ.0) CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle(iPart,1:6))
!      SELECT CASE(PartLorentzType)
!      CASE(0)
!        Pt(iPart,1:3) = NON_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
!        LorentzFacInv = 1.0
!      CASE(1)
!        Pt(iPart,1:3) = SLOW_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
!        LorentzFacInv = 1.0
!      CASE(3)
!        Pt(iPart,1:3) = FAST_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
!        LorentzFacInv = 1.0
!      CASE(5)
!        LorentzFacInv=1.0+DOT_PRODUCT(PartState(iPart,4:6),PartState(iPart,4:6))*c2_inv      
!        LorentzFacInv=1.0/SQRT(LorentzFacInv)
!        Pt(iPart,1:3) = RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6),LorentzFacInvIn=LorentzFacInv)
!      CASE DEFAULT
!      CALL abort(&
!__STAMP__&
!,' Given PartLorentzType does not exist!',PartLorentzType)
!      END SELECT
!      R_PartXK(1,iPart)=LorentzFacInv*PartState(iPart,4)
!      R_PartXK(2,iPart)=LorentzFacInv*PartState(iPart,5)
!      R_PartXK(3,iPart)=LorentzFacInv*PartState(iPart,6)
!      R_PartXK(4,iPart)=Pt(iPart,1)
!      R_PartXK(5,iPart)=Pt(iPart,2)
!      R_PartXK(6,iPart)=Pt(iPart,3)
!      F_PartXK(:,iPart)=PartState(iPart,:) - PartQ(:,iPart) - coeff*R_PartXK(:,iPart)
!      ! vector dot product 
!      CALL PartVectorDotProduct(F_PartXK(:,iPart),F_PartXK(:,iPart),Norm2_F_PartXK(iPart))
!      IF((Norm2_F_PartXK(iPart).LT.AbortTol*Norm2_F_PartX0(iPart)).OR.(Norm2_F_PartXK(iPart).LT.1e-12)) &
!        DoPartInNewton(iPart)=.FALSE.
!      ELSE
!      !IF(nInnerPartNewton.GT.20)THEN
!      !  IPWRITE(*,*) 'blubb',iPart, Norm2_F_PartXK(iPart),Norm2_F_PartX0(iPart)
!      !END IF
!    END IF
!  END DO
  DoNewton=.FALSE.
  IF(ANY(DoPartInNewton)) DoNewton=.TRUE.
#ifdef MPI
  !set T if at least 1 proc has to do newton
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoNewton,1,MPI_LOGICAL,MPI_LOR,PartMPI%COMM,iError) 
#endif /*MPI*/
  IF(DoPrintConvInfo)THEN
    Counter=0
    DO iPart=1,PDM%ParticleVecLength
      IF(DoPartInNewton(iPart))THEN
        Counter=Counter+1      
      END IF ! ParticleInside
    END DO ! iPart
#ifdef MPI
    !set T if at least 1 proc has to do newton
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,Counter,1,MPI_INTEGER,MPI_SUM,PartMPI%COMM,iError) 
#endif /*MPI*/
    !IF(DoPrintConvInfo)THEN
    !  SWRITE(*,*) 'PartNewton',nInnerPartNewton,Counter
    !END IF
  END IF
END DO

IF(DoPrintConvInfo)THEN
  IF (nInnerPartNewton.EQ.nPartNewtonIter) THEN
    SWRITE(*,*) 'PartNewton-not done!',nInnerPartNewton,Counter
    DO iPart=1,PDM%ParticleVecLength
      IF(DoPartInNewton(iPart))THEN
        SWRITE(*,*) ' Failed Particle: ',iPart
        SWRITE(*,*) ' Failed Position: ',PartState(iPart,1:6)
        SWRITE(*,*) ' relative Norm:   ',Norm2_F_PartXK(iPart)/Norm2_F_PartX0(iPart)
        SWRITE(*,*) ' removing particle !   '
        PDM%ParticleInside(iPart)=.FALSE.
  CALL abort(&
__STAMP__&
,' abagsdagfha')
      END IF ! ParticleInside
    END DO ! iPart
  ELSE
    SWRITE(*,*) 'PartNewton',nInnerPartNewton,Counter
  END IF
END IF
nPartNewton=nPartNewton+nInnerPartNewton
!IF (nInnerPartNewton.EQ.nPartNewtonIter) THEN
!  IF(PartMPI%MPIRoot)THEN
!  CALL abort(&
!__STAMP__&
!,'NEWTON NOT CONVERGED WITH NEWTON ITERATIONS,EPISLON',nInnerPartNewton,AbortTol)
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
USE MOD_LinearSolver_Vars,    ONLY: epsPartlinSolver,TotalPartIterLinearSolver
USE MOD_LinearSolver_Vars,    ONLY: nKDim,nRestarts,nPartInnerIter,EisenstatWalker
USE MOD_LinearOperator,       ONLY: PartMatrixVector, PartVectorDotProduct
USE MOD_TimeDisc_Vars,        ONLY: iter
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
REAL              :: W(1:6)
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
!IF(iter.EQ.0) THEN
!  AbortCrit=epsPartlinSolver
!ELSE
IF (.NOT.EisenstatWalker) THEN
  AbortCrit=epsPartlinSolver
END IF
!END IF
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
        IF(nPartInnerIter.GT.1)THEN
        print*,'nPartInnerIter - in GMRES',nPartInnerIter
        END IF
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
USE MOD_LinearOperator,          ONLY:PartMatrixVector, PartVectorDotProduct
USE MOD_Particle_Vars,           ONLY:PartState,F_PartXK,Norm2_F_PartXK,PartQ,PartLorentzType,DoPartInNewton,PartLambdaAccept &
                                     ,PartDeltaX,PEM,PDM,LastPartPos,StagePartPos,Pt,Norm2_F_PartX0
USE MOD_LinearSolver_Vars,       ONLY:reps0,PartXK,R_PartXK
USE MOD_LinearSolver_Vars,       ONLY:Part_alpha, Part_sigma
USE MOD_Part_RHS,                ONLY:SLOW_RELATIVISTIC_PUSH,FAST_RELATIVISTIC_PUSH &
                                     ,RELATIVISTIC_PUSH,NON_RELATIVISTIC_PUSH
USE MOD_PICInterpolation,        ONLY:InterpolateFieldToSingleParticle
USE MOD_PICInterpolation_Vars,   ONLY:FieldAtParticle
USE MOD_Equation_Vars,           ONLY:c2_inv
USE MOD_Particle_Tracking_vars,  ONLY:DoRefMapping
USE MOD_Particle_Tracking,       ONLY:ParticleTracing,ParticleRefTracking
#ifdef MPI
USE MOD_Particle_MPI,            ONLY:IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
USE MOD_Particle_MPI_Vars,       ONLY:PartMPI
USE MOD_Particle_MPI_Vars,       ONLY:ExtPartState,ExtPartSpecies,ExtPartMPF,ExtPartToFIBGM,NbrOfExtParticles
#endif /*MPI*/

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
INTEGER                      :: iPart,iCounter,iCounter2
REAL                         :: lambda, Norm2_PartX,DeltaX_Norm
REAL                         :: LorentzFacInv,Xtilde(1:6)
LOGICAL                      :: DoSetLambda
INTEGER                      :: nLambdaReduce,nMaxLambdaReduce=50
!===================================================================================================================================

lambda=1.
print*,'lambda and netwon iter', lambda, nInnerPartNewton
DoSetLambda=.TRUE.
iCounter=0
PartLambdaAccept=.TRUE.
DO iPart=1,PDM%ParticleVecLength
  IF(DoPartInNewton(iPart))THEN
    ! update the last part pos and element for particle movement
    LastPartPos(iPart,1)=StagePartPos(iPart,1)
    LastPartPos(iPart,2)=StagePartPos(iPart,2)
    LastPartPos(iPart,3)=StagePartPos(iPart,3)
    PEM%lastElement(iPart)=PEM%StageElement(iPart)
    ! update particle position
    !PartXK(1:6,iPart)=PartXK(1:6,iPart)+lambda*PartDeltaX(1:6,iPart)
    !PartState(iPart,1:6)=PartXK(1:6,iPart)+lambda*PartDeltaX(1:6,iPart)
    ! to not check amout of particle movement
    ! this cannot be done here, because the particle has to be 
    ! moved and communicated, if required, hence, it is commented out
!    ! not YET !!!
!    DeltaX_Norm=DOT_PRODUCT(PartDeltaX(1:6,iPart),PartDeltaX(1:6,iPart))
!    IF(DeltaX_Norm.LT.AbortTol*Norm2_F_PartX0(iPart)) THEN
!      DoPartInNewton(iPart)=.FALSE.
!    ELSE
!      PartLambdaAccept(iPart)=.FALSE.
!      iCounter=iCounter+1
!    END IF
    ! new part: of Armijo algorithm: check convergence
    ! compute new function value
    CALL PartMatrixVector(t,Coeff,iPart,PartDeltaX(:,iPart),Xtilde) ! coeff*Ut+Source^n+1 ! only output
    XTilde=XTilde+F_PartXK(1:6,iPart)
    CALL PartVectorDotProduct(Xtilde,Xtilde,Norm2_PartX)
    IF(Norm2_PartX.GT.AbortTol*Norm2_F_PartXK(iPart))THEN
      Norm2_PartX = Norm2_PartX/Norm2_F_PartXk(iPart)
      CALL abort(&
__STAMP__&
,' Found wrond search direction! Particle, Monitored decrease: ', iPart, Norm2_PartX) 
    END IF
    ! update position
    PartState(iPart,1:6)=PartXK(1:6,iPart)+lambda*PartDeltaX(1:6,iPart)
    PartLambdaAccept(iPart)=.FALSE.
    iCounter=iCounter+1
  END IF ! ParticleInside
END DO ! iPart
print*,'iCounters',iCounter
print*,'part-alpah',Part_alpha

! move particle
#ifdef MPI
! open receive buffer for number of particles
CALL IRecvNbofParticles() ! input value: which list:PartLambdaAccept or PDM%ParticleInisde?
#endif /*MPI*/
IF(DoRefMapping)THEN
  CALL ParticleRefTracking(doParticle_In=.NOT.PartLambdaAccept(1:PDM%ParticleVecLength)) 
ELSE
  CALL ParticleTracing(doParticle_In=.NOT.PartLambdaAccept(1:PDM%ParticleVecLength)) 
END IF
DO iPart=1,PDM%ParticleVecLength
  IF(.NOT.PartLambdaAccept(iPart))THEN
    IF(.NOT.PDM%ParticleInside(iPart))THEN
      DoPartInNewton(iPart)=.FALSE.
      PartLambdaAccept(iPart)=.TRUE.
    END IF
  END IF
END DO
#ifdef MPI
! send number of particles
CALL SendNbOfParticles(doParticle_In=.NOT.PartLambdaAccept(1:PDM%ParticleVecLength)) 
! finish communication of number of particles and send particles
CALL MPIParticleSend() ! input value: which list:PartLambdaAccept or PDM%ParticleInisde?
! finish communication
CALL MPIParticleRecv() ! input value: which list:PartLambdaAccept or PDM%ParticleInisde?
! as we do not have the shape function here, we have to deallocate something
SDEALLOCATE(ExtPartState)
SDEALLOCATE(ExtPartSpecies)
SDEALLOCATE(ExtPartToFIBGM)
SDEALLOCATE(ExtPartMPF)
NbrOfExtParticles=0
#endif

DO iPart=1,PDM%ParticleVecLength
  IF(.NOT.PartLambdaAccept(iPart))THEN
    CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle(iPart,1:6))
    SELECT CASE(PartLorentzType)
    CASE(0)
      Pt(iPart,1:3) = NON_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
      LorentzFacInv = 1.0
    CASE(1)
      Pt(iPart,1:3) = SLOW_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
      LorentzFacInv = 1.0
    CASE(3)
      Pt(iPart,1:3) = FAST_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
      LorentzFacInv = 1.0
    CASE(5)
      LorentzFacInv=1.0+DOT_PRODUCT(PartState(iPart,4:6),PartState(iPart,4:6))*c2_inv      
      LorentzFacInv=1.0/SQRT(LorentzFacInv)
      Pt(iPart,1:3) = RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6),LorentzFacInvIn=LorentzFacInv)
    CASE DEFAULT
    CALL abort(&
__STAMP__&
,' Given PartLorentzType does not exist!',PartLorentzType)
    END SELECT
    R_PartXK(1,iPart)=LorentzFacInv*PartState(iPart,4)
    R_PartXK(2,iPart)=LorentzFacInv*PartState(iPart,5)
    R_PartXK(3,iPart)=LorentzFacInv*PartState(iPart,6)
    R_PartXK(4,iPart)=Pt(iPart,1)
    R_PartXK(5,iPart)=Pt(iPart,2)
    R_PartXK(6,iPart)=Pt(iPart,3)
    F_PartXK(1:6,iPart)=PartState(iPart,1:6) - PartQ(1:6,iPart) - coeff*R_PartXK(1:6,iPart)
    ! if check, then here!
    DeltaX_Norm=DOT_PRODUCT(PartDeltaX(1:6,iPart),PartDeltaX(1:6,iPart))
    IF(DeltaX_Norm.LT.AbortTol*Norm2_F_PartX0(iPart)) THEN
       DoPartInNewton(iPart)=.FALSE.
       PartLambdaAccept(iPart)=.TRUE.
       PartXK(1:6,iPart)=PartState(iPart,1:6)
    ELSE
!      IF(nInnerPartNewton.EQ.1)THEN
!        ! accept lambda
!        PartLambdaAccept(iPart)=.TRUE.
!        ! set  new position
!        PartXK(1:6,iPart)=PartState(iPart,1:6)
!        ! update norm
!        CALL PartVectorDotProduct(F_PartXK(1:6,iPart),F_PartXK(1:6,iPart),Norm2_PartX)
!        Norm2_F_PartXK(iPart)=Norm2_PartX
!        IF((Norm2_F_PartXK(iPart).LT.AbortTol*Norm2_F_PartX0(iPart)).OR.(Norm2_F_PartXK(iPart).LT.1e-12)) &
!            DoPartInNewton(iPart)=.FALSE.
!      ELSE
        ! check if residual is reduced
        CALL PartVectorDotProduct(F_PartXK(1:6,iPart),F_PartXK(1:6,iPart),Norm2_PartX)
        IF(Norm2_PartX .LT. (1.-Part_alpha*lambda)*Norm2_F_PartXK(iPart))THEN
          ! accept lambda
          PartLambdaAccept(iPart)=.TRUE.
          ! set  new position
          PartXK(1:6,iPart)=PartState(iPart,1:6)
          Norm2_F_PartXK(iPart)=Norm2_PartX
          IF((Norm2_F_PartXK(iPart).LT.AbortTol*Norm2_F_PartX0(iPart)).OR.(Norm2_F_PartXK(iPart).LT.1e-12)) &
              DoPartInNewton(iPart)=.FALSE.
        ELSE

        END IF
!     END IF ! nInnerPartNewton>1
    END IF
  END IF
END DO ! iPart=1,PDM%ParticleVecLength

iCounter=0
iCounter2=0
DO iPart=1,PDM%ParticleVecLength
  IF(.NOT.PartLambdaAccept(iPart))THEN
    iCounter=iCounter+1
  ELSE
    iCounter2=iCounter2+1
  END IF ! ParticleInside
END DO ! iPart
print*,'iCounters',iCounter,iCounter2

DoSetLambda=.FALSE.
print*,'DoSetLambda',DoSetLambda
IF(ANY(.NOT.PartLambdaAccept)) DoSetLambda=.TRUE.
print*,'DoSetLambda',DoSetLambda
#ifdef MPI
!set T if at least 1 proc has to do newton
CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoSetLambda,1,MPI_LOGICAL,MPI_LOR,PartMPI%COMM,iError)
#endif /*MPI*/

nLambdaReduce=1
DO WHILE((DoSetLambda).AND.(nLambdaReduce.LE.nMaxLambdaReduce))
  nLambdaReduce=nLambdaReduce+1
  lambda=0.1*lambda
  print*,'lambda', lambda
  DO iPart=1,PDM%ParticleVecLength
    IF(.NOT.PartLambdaAccept(iPart))THEN
      ! update the last part pos and element for particle movement
      LastPartPos(iPart,1)=StagePartPos(iPart,1)
      LastPartPos(iPart,2)=StagePartPos(iPart,2)
      LastPartPos(iPart,3)=StagePartPos(iPart,3)
      PEM%lastElement(iPart)=PEM%StageElement(iPart)
      PartState(iPart,1:6)=PartXK(:,iPart)+lambda*PartDeltaX(:,iPart)
      PartLambdaAccept(iPart)=.FALSE.
    END IF ! ParticleInside
  END DO ! iPart
  
  ! move particle
#ifdef MPI
  ! open receive buffer for number of particles
  CALL IRecvNbofParticles() ! input value: which list:PartLambdaAccept or PDM%ParticleInisde?
#endif /*MPI*/
  IF(DoRefMapping)THEN
    CALL ParticleRefTracking(doParticle_In=.NOT.PartLambdaAccept(1:PDM%ParticleVecLength)) 
  ELSE
    CALL ParticleTracing(doParticle_In=.NOT.PartLambdaAccept(1:PDM%ParticleVecLength)) 
  END IF
  DO iPart=1,PDM%ParticleVecLength
    IF(.NOT.PartLambdaAccept(iPart))THEN
      IF(.NOT.PDM%ParticleInside(iPart))THEN
        DoPartInNewton(iPart)=.FALSE.
        PartLambdaAccept(iPart)=.TRUE.
      END IF
    END IF
  END DO
#ifdef MPI
  ! send number of particles
  CALL SendNbOfParticles(doParticle_In=.NOT.PartLambdaAccept(1:PDM%ParticleVecLength)) 
  ! finish communication of number of particles and send particles
  CALL MPIParticleSend() ! input value: which list:PartLambdaAccept or PDM%ParticleInisde?
  ! finish communication
  CALL MPIParticleRecv() ! input value: which list:PartLambdaAccept or PDM%ParticleInisde?
  ! as we do not have the shape function here, we have to deallocate something
  SDEALLOCATE(ExtPartState)
  SDEALLOCATE(ExtPartSpecies)
  SDEALLOCATE(ExtPartToFIBGM)
  SDEALLOCATE(ExtPartMPF)
  NbrOfExtParticles=0
#endif

  DO iPart=1,PDM%ParticleVecLength
    IF(.NOT.PartLambdaAccept(iPart))THEN
      CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle(iPart,1:6))
      SELECT CASE(PartLorentzType)
      CASE(0)
        Pt(iPart,1:3) = NON_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
        LorentzFacInv = 1.0
      CASE(1)
        Pt(iPart,1:3) = SLOW_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
        LorentzFacInv = 1.0
      CASE(3)
        Pt(iPart,1:3) = FAST_RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6))
        LorentzFacInv = 1.0
      CASE(5)
        LorentzFacInv=1.0+DOT_PRODUCT(PartState(iPart,4:6),PartState(iPart,4:6))*c2_inv      
        LorentzFacInv=1.0/SQRT(LorentzFacInv)
        Pt(iPart,1:3) = RELATIVISTIC_PUSH(iPart,FieldAtParticle(iPart,1:6),LorentzFacInvIn=LorentzFacInv)
      CASE DEFAULT
      CALL abort(&
  __STAMP__&
  ,' Given PartLorentzType does not exist!',PartLorentzType)
      END SELECT
      R_PartXK(1,iPart)=LorentzFacInv*PartState(iPart,4)
      R_PartXK(2,iPart)=LorentzFacInv*PartState(iPart,5)
      R_PartXK(3,iPart)=LorentzFacInv*PartState(iPart,6)
      R_PartXK(4,iPart)=Pt(iPart,1)
      R_PartXK(5,iPart)=Pt(iPart,2)
      R_PartXK(6,iPart)=Pt(iPart,3)
      F_PartXK(:,iPart)=PartState(iPart,:) - PartQ(:,iPart) - coeff*R_PartXK(:,iPart)
      ! vector dot product 
      CALL PartVectorDotProduct(F_PartXK(:,iPart),F_PartXK(:,iPart),Norm2_PartX)
      IF(Norm2_PartX .LT. (1.-Part_alpha*lambda)*Norm2_F_PartXK(iPart))THEN
        ! accept lambda
        PartLambdaAccept(iPart)=.TRUE.
        ! set  new position
        PartXK(1:6,iPart)=PartState(iPart,1:6)
        Norm2_F_PartXK(iPart)=Norm2_PartX
        IF((Norm2_F_PartXK(iPart).LT.AbortTol*Norm2_F_PartX0(iPart)).OR.(Norm2_F_PartXK(iPart).LT.1e-12)) &
          DoPartInNewton(iPart)=.FALSE.
      END IF
    END IF
  END DO ! iPart=1,PDM%ParticleVecLength
  ! detect  convergence
  DoSetLambda=.FALSE.
  IF(ANY(.NOT.PartLambdaAccept)) DoSetLambda=.TRUE.
#ifdef MPI
  !set T if at least 1 proc has to do newton
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoSetLambda,1,MPI_LOGICAL,MPI_LOR,PartMPI%COMM,iError)
#endif /*MPI*/
  iCounter=0
  DO iPart=1,PDM%ParticleVecLength
    IF(.NOT.PartLambdaAccept(iPart))THEN
      PartLambdaAccept(iPart)=.FALSE.
      iCounter=iCounter+1
!      IF(nLambdaReduce.GT.1)THEN
!        CALL PartVectorDotProduct(F_PartXK(:,iPart),F_PartXK(:,iPart),Norm2_PartX)
!        print*,'Norm2_PartX',Norm2_PartX
!        print*,'norm-fcurrent',Norm2_F_PartXK(iPart)
!        print*,'norm-fold',Norm2_F_PartX0(iPart)
!        print*,Norm2_PartX/(1.-Part_alpha*lambda)*Norm2_F_PartXK(iPart)
!      END IF
    END IF ! ParticleInside
  END DO ! iPart
  print*,'iCounters-in setting',iCounter
END DO

END SUBROUTINE Particle_Armijo


SUBROUTINE FinalizePartSolver() 
!===================================================================================================================================
! deallocate global variables
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
USE MOD_LinearSolver_Vars
USE MOD_Particle_Vars,           ONLY:PartQ,F_PartX0,F_PartXk,Norm2_F_PartX0,Norm2_F_PartXK,Norm2_F_PartXK_old,DoPartInNewton &
                                     ,PartDeltaX,PartLambdaAccept
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
SDEALLOCATE(PartLambdaAccept)
SDEALLOCATE(PartDeltaX)
SDEALLOCATE(Norm2_F_PartX0)
SDEALLOCATE(Norm2_F_PartXK)
SDEALLOCATE(Norm2_F_PartXK_old)
SDEALLOCATE(DoPartInNewton)
END SUBROUTINE FinalizePartSolver
#endif /*PARTICLES*/

END MODULE MOD_ParticleSolver
