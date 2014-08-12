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

SUBROUTINE GetBoundaryInteractionSuperSampled(PartTrajectory,alpha,xi,eta,iPart,QuadID,SideID,ElemID)
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
USE MOD_Particle_Analyze_Vars,  ONLY:CalcPartBalance,nPartOut,PartEkinOut,PartAnalyzeStep
USE MOD_TimeDisc_Vars,          ONLY:iter
USE MOD_Mesh_Vars,              ONLY:BC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: iPart,SideID,QuadID
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

  n_loc=CalcNormVec(xi,eta,QuadID,SideID)
  !print*,'reflective BC'
  !print*,'ElemId,SideID',ElemID,SideID,QuadID
  !print*,PartTrajectory
  !print*,'nVec',n_loc
  !read*
  ! intersection point with surface
  LastPartPos(iPart,1:3) = LastPartPos(iPart,1:3) + PartTrajectory(1:3)*alpha
  !print*,'alpha,inter',alpha,lastPartPos(iPart,1:3)
  ! In vector notation: r_neu = r_alt + T - 2*((1-alpha)*<T,n>)*n
  !v_aux = - 2*((1-alpha)*<T,n>)*n     (auxiliary variable, used twice)
  v_aux                  = -2*((1-alpha)*DOT_PRODUCT(PartTrajectory(1:3),n_loc))*n_loc
  !PartState(iPart,1:3)   = PartState(iPart,1:3)+PartTrajectory(1:3)+v_aux
  PartState(iPart,1:3)   = PartState(iPart,1:3)+v_aux
  !PartState(iPart,1:3)   = LastPartPos(iPart,1:3)+v_aux
  ! new velocity vector 
  v_2=(1-alpha)*PartTrajectory(1:3)+v_aux

  !print*,"PartBound%WallVelo(1:3,BC(SideID))",PartBound%WallVelo(1:3,BC(SideID))
  !read*
  PartState(iPart,4:6)   = SQRT(DOT_PRODUCT(PartState(iPart,4:6),PartState(iPart,4:6)))*&
                           (1/(SQRT(DOT_PRODUCT(v_2,v_2))))*v_2                         +&
                           PartBound%WallVelo(1:3,BC(SideID))
  !PartState(iPart,4:6)   = 0.
  PartTrajectory=PartState(iPart,1:3) - LastPartPos(iPart,1:3)
!  print*,PartTrajectory
!  read*
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
