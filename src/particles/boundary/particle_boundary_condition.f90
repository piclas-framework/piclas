#include "boltzplatz.h"

MODULE MOD_Particle_Boundary_Condition
!===================================================================================================================================
!! Determines how particles interact with a given boundary condition 
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
INTERFACE GetBoundaryInteraction
  MODULE PROCEDURE GetBoundaryInteraction
END INTERFACE

INTERFACE GetBoundaryInteractionSuperSampled
  MODULE PROCEDURE GetBoundaryInteractionSuperSampled
END INTERFACE

PUBLIC::GetBoundaryInteraction,GetBoundaryInteractionSuperSampled
!===================================================================================================================================

CONTAINS

SUBROUTINE GetBoundaryInteractionSuperSampled(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,QuadID,SideID,ElemID)
!===================================================================================================================================
! Computes the post boundary state of a particle that interacts with a boundary condition
!  OpenBC                  = 1  
!  ReflectiveBC            = 2  
!  PeriodicBC              = 3  
!  SimpleAnodeBC           = 4  
!  SimpleCathodeBC         = 5  
!  MPINeighborhoodBC       = 6  
!===================================================================================================================================
! MODULES
USE MOD_Globals,                ONLY:Abort
USE MOD_PreProc
USE MOD_Particle_Surfaces,      ONLY:CalcNormVec
USE MOD_Particle_Vars,          ONLY:PartBound,PDM,PartSpecies,PartState,LastPartPos
USE MOD_Particle_Analyze,       ONLY:CalcEkinPart
!USE MOD_Equation_Vars,          ONLY:c2_inv
USE MOD_Particle_Vars,          ONLY:PartState!, PartSpecies, Species
!USE MOD_PARTICLE_Vars,          ONLY:PartMPF, usevMPF
USE MOD_Particle_Analyze_Vars,  ONLY:CalcPartBalance,nPartOut,PartEkinOut,PartAnalyzeStep
USE MOD_Particle_Surfaces_Vars, ONLY:epsilontol
USE MOD_Particle_Vars,          ONLY:Pt_temp,Pt
USE MOD_TimeDisc_Vars,          ONLY:iter
USE MOD_Mesh_Vars,              ONLY:BC
#if (PP_TimeDiscMethod==1) || (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6)
USE MOD_TimeDisc_Vars,          ONLY: RK4_a,iStage
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: iPart,SideID,QuadID
REAL,INTENT(IN)                      :: xi,eta
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                :: ElemID
REAL,INTENT(INOUT)                   :: alpha,PartTrajectory(1:3),LengthPartTrajectory
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: v_2(1:3),v_aux(1:3),n_loc(1:3)
!REAL                                 :: OldP(1:3),NewP(1:3),gamma1,absOldP,absNewP,MPFMass
#if (PP_TimeDiscMethod==1) || (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6)
REAL                                 :: absPt_temp
#endif
!===================================================================================================================================

IF (.NOT. ASSOCIATED(PartBound%Map)) THEN
  CALL abort(__STAMP__,&
  ' ERROR: PartBound not associated!.',999,999.)
END IF

! Select the corresponding boundary condition and calculate particle treatment
SELECT CASE(PartBound%Map(BC(SideID)))
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(1) !PartBound%OpenBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  IF(CalcPartBalance) THEN
    IF(MOD(iter+1,PartAnalyzeStep).EQ.0)THEN ! caution if correct
      nPartOut(PartSpecies(iPart))=nPartOut(PartSpecies(iPart)) + 1
      PartEkinOut(PartSpecies(iPart))=PartEkinOut(PartSpecies(iPart))+CalcEkinPart(iPart)
    END IF ! iter+1
  END IF ! CalcPartBalance
  PDM%ParticleInside(iPart) = .FALSE.
  alpha=-1.
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(2) !PartBound%ReflectiveBC)
!-----------------------------------------------------------------------------------------------------------------------------------

  n_loc=CalcNormVec(xi,eta,QuadID,SideID)
  ! substract tolerance from length
  LengthPartTrajectory=LengthPartTrajectory-epsilontol
  ! compute old relativistic impulse
!  gamma1=PartState(iPart,4)*PartState(iPart,4)+PartState(iPart,5)*PartState(iPart,5) &
!        +PartState(iPart,6)*PartState(iPart,6)
!  gamma1=gamma1*c2_inv
!  gamma1=1.0/SQRT(1.-gamma1)
!  IF(usevMPF)THEN
!    MPFMass=PartMPF(iPart)*Species(PartSpecies(iPart))%MassIC
!  ELSE
!    MPFMass=Species(PartSpecies(iPart))%MassIC*Species(PartSpecies(iPart))%MacroParticleFactor
!  END IF ! usevMPF
!  oldP=MPFMass*gamma1*PartState(iPart,4:6)
!  absOldP=SQRT(DOT_PRODUCT(oldP,oldP))

!  WRITE(*,'(A)'),'reflective BC'
!  WRITE(*,'(A,4I,4I,4I)') 'ElemId,SideID,QuadiD',ElemID,SideID,QuadID
!  WRITE(*,'(A,E24.15,E24.15,E24.15)') 'alpha,length,percent',alpha,LengthPartTrajectory,alpha/LengthPartTrajectory
!  WRITE(*,'(A,E24.15,E24.15,E24.15)') 'old path',PartTrajectory
!  WRITE(*,'(A,E24.15,E24.15,E24.15)') 'old pos',LastPartPos(iPart,1:3)
!  WRITE(*,'(A,E24.15,E24.15,E24.15)') 'outside state',PartState(iPart,1:3)
!  WRITE(*,'(A,E24.15,E24.15,E24.15)') 'nVec',n_loc
!  WRITE(*,'(A,E24.15)') 'nVec o PartTrajectory',DOT_PRODUCT(n_loc,PartTrajectory)
!  WRITE(*,'(A,E24.15,E24.15,E24.15)') 'BC velo',PartBound%WallVelo(1:3,BC(SideID))
  ! intersection point with surface
  LastPartPos(iPart,1:3) = LastPartPos(iPart,1:3) + PartTrajectory(1:3)*alpha
  ! In vector notation: r_neu = r_alt + T - 2*((1-alpha)*<T,n>)*n
  !v_aux = - 2*((1-alpha)*<T,n>)*n     (auxiliary variable, used twice)
  v_aux                  = -2*((LengthPartTrajectory-alpha)*DOT_PRODUCT(PartTrajectory(1:3),n_loc))*n_loc
  !PartState(iPart,1:3)   = PartState(iPart,1:3)+PartTrajectory(1:3)+v_aux
  PartState(iPart,1:3)   = PartState(iPart,1:3)+v_aux
  !PartState(iPart,1:3)   = LastPartPos(iPart,1:3)+v_aux
  ! new velocity vector 
  v_2=(LengthPartTrajectory-alpha)*PartTrajectory(1:3)+v_aux

  !print*,"PartBound%WallVelo(1:3,BC(SideID))",PartBound%WallVelo(1:3,BC(SideID))
  PartState(iPart,4:6)   = SQRT(DOT_PRODUCT(PartState(iPart,4:6),PartState(iPart,4:6)))*&
                           (1/(SQRT(DOT_PRODUCT(v_2,v_2))))*v_2                         +&
                           PartBound%WallVelo(1:3,BC(SideID))
  !PartState(iPart,4:6)   = 0.
  PartTrajectory=PartState(iPart,1:3) - LastPartPos(iPart,1:3)
  ! new normalization
  lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                           +PartTrajectory(2)*PartTrajectory(2) &
                           +PartTrajectory(3)*PartTrajectory(3) )
  PartTrajectory=PartTrajectory/lengthPartTrajectory
  lengthPartTrajectory=lengthPartTrajectory+epsilontol
!  WRITE(*,'(A,E24.15,E24.15,E24.15)') 'pos at BC',LastPartPos(iPart,1:3)
!  WRITE(*,'(A,E24.15,E24.15,E24.15)') 'new pos',PartState(iPart,1:3)
!  WRITE(*,'(A,E24.15,E24.15,E24.15)') 'new velo1',PartState(iPart,4:6)
!  WRITE(*,'(A,E24.15,E24.15,E24.15)') 'new velo2,n',PartTrajectory
!  WRITE(*,'(A,E24.15,E24.15,E24.15)') 'sanity check'
!  read*

!  gamma1=PartState(iPart,4)*PartState(iPart,4)+PartState(iPart,5)*PartState(iPart,5) &
!        +PartState(iPart,6)*PartState(iPart,6)
!  gamma1=gamma1*c2_inv
!  gamma1=1.0/SQRT(1.-gamma1)
!  newP=MPFMass*gamma1*PartState(iPart,4:6)
!  absNewP=SQRT(DOT_PRODUCT(newP,newP))
!  print*,'ratio impulse',absOldP/absNewP
!  read*
#if (PP_TimeDiscMethod==1) || (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6)
  ! correction for Runge-Kutta (correct position!!)
  !print*,'Pt_temp',Pt_temp(iPart,1:3)
  ! get length of Pt_temp(iPart,1:3) || equals summed velocity change ! only exact for linear movement
!  print*,'acceleration_old',Pt_temp(iPart,4:6)
!  read*
  absPt_temp=SQRT(Pt_temp(iPart,1)*Pt_temp(iPart,1)+Pt_temp(iPart,2)*Pt_temp(iPart,2)+Pt_temp(iPart,3)*Pt_temp(iPart,3))
  ! scale PartTrajectory to new Pt_temp
  Pt_temp(iPart,1:3)=absPt_temp*PartTrajectory(1:3)
  ! deleate force history
  Pt_temp(iPart,4:6)=0.
  ! what happens with force term || acceleration?
#endif 
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(3) !PartBound%PeriodicBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL abort(__STAMP__,&
  ' ERROR: PartBound not associated!. (PartBound%PeriodicBC)',999,999.)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(4) !PartBound%SimpleAnodeBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL abort(__STAMP__,&
  ' ERROR: PartBound not associated!. (PartBound%SimpleAnodeBC)',999,999.)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(5) !PartBound%SimpleCathodeBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL abort(__STAMP__,&
  ' ERROR: PartBound not associated!. (PartBound%SimpleCathodeBC)',999,999.)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(6) !PartBound%MPINeighborhoodBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL abort(__STAMP__,&
  ' ERROR: PartBound not associated!. (PartBound%MPINeighborhoodBC)',999,999.)
CASE DEFAULT
  CALL abort(__STAMP__,&
' ERROR: PartBound not associated!. (unknown case)',999,999.)
END SELECT !PartBound%Map(BC(SideID)

END SUBROUTINE GetBoundaryInteractionSuperSampled

SUBROUTINE GetBoundaryInteraction(PartTrajectory,alpha,xi,eta,iPart,SideID,ElemID)
!===================================================================================================================================
! Computes the post boundary state of a particle that interacts with a boundary condition
!  OpenBC                  = 1  
!  ReflectiveBC            = 2  
!  PeriodicBC              = 3  
!  SimpleAnodeBC           = 4  
!  SimpleCathodeBC         = 5  
!  MPINeighborhoodBC       = 6  
!===================================================================================================================================
! MODULES
USE MOD_Globals,                ONLY:Abort
USE MOD_PreProc
USE MOD_Particle_Surfaces,      ONLY:CalcBiLinearNormVec
USE MOD_Particle_Vars,          ONLY:PartBound,PDM,PartSpecies,PartState,LastPartPos
USE MOD_Particle_Surfaces_vars, ONLY:SideNormVec,SideIsPlanar
USE MOD_Particle_Analyze,       ONLY:CalcEkinPart
USE MOD_Particle_Analyze_Vars,  ONLY:CalcPartBalance,nPartOut,PartEkinOut,PartAnalyzeStep
USE MOD_TimeDisc_Vars,          ONLY:iter
USE MOD_Mesh_Vars,              ONLY:BC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: iPart,SideID
REAL,INTENT(IN)                      :: xi,eta
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                :: ElemID
REAL,INTENT(INOUT)                   :: alpha,PartTrajectory(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: v_2(1:3),v_aux(1:3),n_loc(1:3)
!===================================================================================================================================

IF (.NOT. ASSOCIATED(PartBound%Map)) THEN
  CALL abort(__STAMP__,&
  ' ERROR: PartBound not associated!.',999,999.)
END IF

! Select the corresponding boundary condition and calculate particle treatment
SELECT CASE(PartBound%Map(BC(SideID)))
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(1) !PartBound%OpenBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  IF(CalcPartBalance) THEN
    IF(MOD(iter+1,PartAnalyzeStep).EQ.0)THEN ! caution if correct
      nPartOut(PartSpecies(iPart))=nPartOut(PartSpecies(iPart)) + 1
      PartEkinOut(PartSpecies(iPart))=PartEkinOut(PartSpecies(iPart))+CalcEkinPart(iPart)
    END IF ! iter+1
  END IF ! CalcPartBalance
  PDM%ParticleInside(iPart) = .FALSE.
  alpha=-1.
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(2) !PartBound%ReflectiveBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  IF(SideIsPlanar(SideID))THEN
    n_loc=SideNormVec(1:3,SideID)
  ELSE
    n_loc=CalcBiLinearNormVec(xi,eta,SideID)
  END IF
  ! intersection point with surface
  lastPartPos(iPart,1:3) = PartState(iPart,1:3)+PartTrajectory(1:3)*alpha
  ! In vector notation: r_neu = r_alt + T - 2*((1-alpha)*<T,n>)*n
  !v_aux = - 2*((1-alpha)*<T,n>)*n     (auxiliary variable, used twice)
  v_aux                  = -2*((1-alpha)*DOT_PRODUCT(PartTrajectory(1:3),n_loc))*n_loc
  PartState(iPart,1:3)   = PartState(iPart,1:3)+PartTrajectory(1:3)+v_aux
  ! new velocity vector 
  v_2=(1-alpha)*PartTrajectory(1:3)+v_aux
  PartState(iPart,4:6)   = SQRT(DOT_PRODUCT(PartState(iPart,4:6), PartState(iPart,4:6)))*&
                           (1/(SQRT(DOT_PRODUCT(v_2,v_2))))*v_2                         +&
                           PartBound%WallVelo(1:3,BC(SideID))
  !PartState(iPart,4:6)   = 0.
  PartTrajectory=PartState(iPart,1:3) - LastPartPos(iPart,1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(3) !PartBound%PeriodicBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL abort(__STAMP__,&
  ' ERROR: PartBound not associated!. (PartBound%PeriodicBC)',999,999.)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(4) !PartBound%SimpleAnodeBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL abort(__STAMP__,&
  ' ERROR: PartBound not associated!. (PartBound%SimpleAnodeBC)',999,999.)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(5) !PartBound%SimpleCathodeBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL abort(__STAMP__,&
  ' ERROR: PartBound not associated!. (PartBound%SimpleCathodeBC)',999,999.)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(6) !PartBound%MPINeighborhoodBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL abort(__STAMP__,&
  ' ERROR: PartBound not associated!. (PartBound%MPINeighborhoodBC)',999,999.)
CASE DEFAULT
  CALL abort(__STAMP__,&
' ERROR: PartBound not associated!. (unknown case)',999,999.)
END SELECT !PartBound%Map(BC(SideID)

END SUBROUTINE GetBoundaryInteraction 

END MODULE MOD_Particle_Boundary_Condition
