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

INTERFACE GetBoundaryInteractionRef
  MODULE PROCEDURE GetBoundaryInteractionRef
END INTERFACE

PUBLIC::GetBoundaryInteraction,GetBoundaryInteractionRef
!===================================================================================================================================

CONTAINS

SUBROUTINE GetBoundaryInteraction(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,ElemID)
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
USE MOD_PreProc
USE MOD_Globals,                ONLY:Abort
USE MOD_Particle_Surfaces,      ONLY:CalcBiLinearNormVecBezier,CalcNormVecBezier
USE MOD_Particle_Vars,          ONLY:PDM,PartSpecies,PartState,LastPartPos,PEM
USE MOD_Particle_Mesh_Vars,     ONLY:PartBound
USE MOD_Particle_Surfaces_vars, ONLY:SideNormVec,SideType,epsilontol
!USE MOD_Particle_Surfaces_Vars, ONLY:BoundingBoxIsEmpty
USE MOD_Particle_Analyze,       ONLY:CalcEkinPart
USE MOD_Particle_Analyze_Vars,  ONLY:CalcPartBalance,nPartOut,PartEkinOut!,PartAnalyzeStep
USE MOD_TimeDisc_Vars,          ONLY:iter
USE MOD_Mesh_Vars,              ONLY:BC
!USE MOD_BoundaryTools,          ONLY:SingleParticleToExactElement                                   !
#if (PP_TimeDiscMethod==1) || (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6)
USE MOD_Particle_Vars,          ONLY:Pt_temp,Pt
USE MOD_TimeDisc_Vars,          ONLY:RK_a,iStage
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: iPart,SideID
REAL,INTENT(IN)                      :: xi,eta
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                :: ElemID
REAL,INTENT(INOUT)                   :: alpha,PartTrajectory(1:3),lengthPartTrajectory
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: v_2(1:3),v_aux(1:3),n_loc(1:3)
#if (PP_TimeDiscMethod==1) || (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6)
REAL                                 :: absPt_temp
#endif
!===================================================================================================================================

IF (.NOT. ALLOCATED(PartBound%Map)) THEN
  CALL abort(__STAMP__,&
  ' ERROR: PartBound not allocated!.',999,999.)
END IF

! Select the corresponding boundary condition and calculate particle treatment
SELECT CASE(PartBound%Map(BC(SideID)))
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(1) !PartBound%OpenBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  IF(CalcPartBalance) THEN
    !IF(MOD(iter+1,PartAnalyzeStep).EQ.0)THEN ! caution if correct
      nPartOut(PartSpecies(iPart))=nPartOut(PartSpecies(iPart)) + 1
      PartEkinOut(PartSpecies(iPart))=PartEkinOut(PartSpecies(iPart))+CalcEkinPart(iPart)
    !END IF ! iter+1
  END IF ! CalcPartBalance
  PDM%ParticleInside(iPart) = .FALSE.
  alpha=-1.
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(2) !PartBound%ReflectiveBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  SELECT CASE(SideType(SideID))
  CASE(PLANAR)
    n_loc=SideNormVec(1:3,SideID)
  CASE(BILINEAR)
    n_loc=CalcBiLinearNormVecBezier(xi,eta,SideID)
  CASE(CURVED)
    n_loc=CalcNormVecBezier(xi,eta,SideID)
!    CALL abort(__STAMP__,'nvec for bezier not implemented!',999,999.)
  END SELECT 
!  print*,'n_loc',n_loc
!  print*,'n_loc,partt',DOT_PRODUCT(n_loc,PartTrajectory)
!  ead*
  ! substract tolerance from length
  LengthPartTrajectory=LengthPartTrajectory-epsilontol
  ! intersection point with surface
  LastPartPos(iPart,1:3) = LastPartPos(iPart,1:3) + PartTrajectory(1:3)*alpha

!  print*,'intersection',lastpartpos(ipart,1:3)
!  stop
  ! In vector notation: r_neu = r_alt + T - 2*((1-alpha)*<T,n>)*n
  !v_aux = - 2*((1-alpha)*<T,n>)*n     (auxiliary variable, used twice)
  !v_aux                  = -2*((1-alpha)*DOT_PRODUCT(PartTrajectory(1:3),n_loc))*n_loc
  v_aux                  = -2*((LengthPartTrajectory-alpha)*DOT_PRODUCT(PartTrajectory(1:3),n_loc))*n_loc
  !PartState(iPart,1:3)   = PartState(iPart,1:3)+PartTrajectory(1:3)+v_aux
  PartState(iPart,1:3)   = PartState(iPart,1:3)+v_aux
  ! new velocity vector 
  !v_2=(1-alpha)*PartTrajectory(1:3)+v_aux
  v_2=(LengthPartTrajectory-alpha)*PartTrajectory(1:3)+v_aux
  PartState(iPart,4:6)   = SQRT(DOT_PRODUCT(PartState(iPart,4:6),PartState(iPart,4:6)))*&
                           (1/(SQRT(DOT_PRODUCT(v_2,v_2))))*v_2                         +&
                           PartBound%WallVelo(1:3,BC(SideID))
  PartTrajectory=PartState(iPart,1:3) - LastPartPos(iPart,1:3)
  lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                           +PartTrajectory(2)*PartTrajectory(2) &
                           +PartTrajectory(3)*PartTrajectory(3) )
  PartTrajectory=PartTrajectory/lengthPartTrajectory
  ! move particle eps along line to prevent a detection of alpha=zero
  LastPartPos(iPart,1:3) = LastPartPos(iPart,1:3)+1.0e-8*PartTrajectory
  lengthPartTrajectory=lengthPartTrajectory+epsilontol-1.0e-8
!  print*, ' oldElemID', ElemID
  ! get new element ID
!  CALL SingleParticleToExactElement(iPart)
!  ElemID=PEM%Element(iPart)
!  PEM%lastElement(iPart) = ElemID
!  print*,' newElemID', ElemID

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
  ! new implementation, nothing to due :)
  ! however, never checked
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


SUBROUTINE GetBoundaryInteractionRef(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,ElemID)
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
USE MOD_PreProc
USE MOD_Globals,                ONLY:Abort
USE MOD_Particle_Surfaces,      ONLY:CalcBiLinearNormVecBezier,CalcNormVecBezier
USE MOD_Particle_Vars,          ONLY:PDM,PartSpecies,PartState,LastPartPos,PEM
USE MOD_Particle_Mesh_Vars,     ONLY:PartBound
USE MOD_Particle_Surfaces_vars, ONLY:SideNormVec,SideType,epsilontol
!USE MOD_Particle_Surfaces_Vars, ONLY:BoundingBoxIsEmpty
USE MOD_Particle_Analyze,       ONLY:CalcEkinPart
USE MOD_Particle_Analyze_Vars,  ONLY:CalcPartBalance,nPartOut,PartEkinOut!,PartAnalyzeStep
USE MOD_TimeDisc_Vars,          ONLY:iter
USE MOD_Mesh_Vars,              ONLY:BC,nSides
#if (PP_TimeDiscMethod==1) || (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6)
USE MOD_Particle_Vars,          ONLY:Pt_temp,Pt
USE MOD_TimeDisc_Vars,          ONLY:RK_a,iStage
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: iPart,SideID
REAL,INTENT(IN)                      :: xi,eta
INTEGER,INTENT(IN)                   :: ElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                   :: alpha,PartTrajectory(1:3),lengthPartTrajectory
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: v_2(1:3),v_aux(1:3),n_loc(1:3)
#if (PP_TimeDiscMethod==1) || (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6)
REAL                                 :: absPt_temp
#endif
!===================================================================================================================================

IF (.NOT. ALLOCATED(PartBound%Map)) THEN
  CALL abort(__STAMP__,&
  ' ERROR: PartBound not allocated!.',999,999.)
END IF

! Select the corresponding boundary condition and calculate particle treatment
SELECT CASE(PartBound%Map(BC(SideID)))
!SELECT CASE(PartBound%SideBCType(SideID))
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(1) !PartBound%OpenBC)
!CASE(PartBound%OpenBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  IF(CalcPartBalance) THEN
    !IF(MOD(iter+1,PartAnalyzeStep).EQ.0)THEN ! caution if correct
      nPartOut(PartSpecies(iPart))=nPartOut(PartSpecies(iPart)) + 1
      PartEkinOut(PartSpecies(iPart))=PartEkinOut(PartSpecies(iPart))+CalcEkinPart(iPart)
    !END IF ! iter+1
  END IF ! CalcPartBalance
  PDM%ParticleInside(iPart) = .FALSE.
  alpha=-1.
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(2) !PartBound%ReflectiveBC)
!CASE(PartBound%ReflectiveBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  SELECT CASE(SideType(SideID))
  CASE(PLANAR)
    n_loc=SideNormVec(1:3,SideID)
  CASE(BILINEAR)
    n_loc=CalcBiLinearNormVecBezier(xi,eta,SideID)
  CASE(CURVED)
    n_loc=CalcNormVecBezier(xi,eta,SideID)
!    CALL abort(__STAMP__,'nvec for bezier not implemented!',999,999.)
  END SELECT 
!  print*,'n_loc',n_loc
!  print*,'n_loc,partt',DOT_PRODUCT(n_loc,PartTrajectory)
!  read*
  ! substract tolerance from length
  LengthPartTrajectory=LengthPartTrajectory-epsilontol
  ! intersection point with surface
  LastPartPos(iPart,1:3) = LastPartPos(iPart,1:3) + PartTrajectory(1:3)*alpha

!  print*,'intersection',lastpartpos(ipart,1:3)
!  stop
  ! In vector notation: r_neu = r_alt + T - 2*((1-alpha)*<T,n>)*n
  !v_aux = - 2*((1-alpha)*<T,n>)*n     (auxiliary variable, used twice)
  !v_aux                  = -2*((1-alpha)*DOT_PRODUCT(PartTrajectory(1:3),n_loc))*n_loc
  v_aux                  = -2*((LengthPartTrajectory-alpha)*DOT_PRODUCT(PartTrajectory(1:3),n_loc))*n_loc
  !PartState(iPart,1:3)   = PartState(iPart,1:3)+PartTrajectory(1:3)+v_aux
  PartState(iPart,1:3)   = PartState(iPart,1:3)+v_aux
  ! new velocity vector 
  !v_2=(1-alpha)*PartTrajectory(1:3)+v_aux
  v_2=(LengthPartTrajectory-alpha)*PartTrajectory(1:3)+v_aux
  PartState(iPart,4:6)   = SQRT(DOT_PRODUCT(PartState(iPart,4:6),PartState(iPart,4:6)))*&
                           (1/(SQRT(DOT_PRODUCT(v_2,v_2))))*v_2                         +&
                           PartBound%WallVelo(1:3,BC(SideID))
  PartTrajectory=PartState(iPart,1:3) - LastPartPos(iPart,1:3)
  lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                           +PartTrajectory(2)*PartTrajectory(2) &
                           +PartTrajectory(3)*PartTrajectory(3) )
  PartTrajectory=PartTrajectory/lengthPartTrajectory
  lengthPartTrajectory=lengthPartTrajectory+epsilontol
!  print*, ' oldElemID', ElemID
  ! get new element ID
!  CALL SingleParticleToExactElement(iPart)
!  ElemID=PEM%Element(iPart)
!  PEM%lastElement(iPart) = ElemID
!  print*,' newElemID', ElemID

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
!CASE(PartBound%PeriodicBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  ! new implementation, nothing to due :)
  ! however, never checked
  CALL abort(__STAMP__,&
  ' ERROR: PartBound not associated!. (PartBound%PeriodicBC)',999,999.)
  !compute new bc
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(4) !PartBound%SimpleAnodeBC)
!CASE(PartBound%SimpleAnodeBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL abort(__STAMP__,&
  ' ERROR: PartBound not associated!. (PartBound%SimpleAnodeBC)',999,999.)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(5) !PartBound%SimpleCathodeBC)
!CASE(PartBound%SimpleCathodeBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL abort(__STAMP__,&
  ' ERROR: PartBound not associated!. (PartBound%SimpleCathodeBC)',999,999.)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(6) !PartBound%MPINeighborhoodBC)
!CASE(PartBound%MPINeighborhoodBC)
  CALL abort(__STAMP__,&
  ' ERROR: PartBound not associated!. (PartBound%MPINeighborhoodBC)',999,999.)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE DEFAULT
  CALL abort(__STAMP__,&
' ERROR: PartBound not associated!. BC(SideID)',BC(SideID),REAL(SideID/nSides))
END SELECT !PartBound%Map(BC(SideID)

END SUBROUTINE GetBoundaryInteractionRef

END MODULE MOD_Particle_Boundary_Condition
