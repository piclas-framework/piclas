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

IF (.NOT. ALLOCATED(PartBound%MapToPartBC)) THEN
  CALL abort(__STAMP__,&
  ' ERROR: PartBound not allocated!.',999,999.)
END IF

! Select the corresponding boundary condition and calculate particle treatment
!SELECT CASE(PartBound%MapToPartBC(BC(SideID)))
SELECT CASE(PartBound%TargetBoundCond(PartBound%MapToPartBC(BC(SideID))))
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
END SELECT !PartBound%MapToPartBC(BC(SideID)

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
USE MOD_Globals!,                ONLY:Abort
USE MOD_Particle_Surfaces,      ONLY:CalcBiLinearNormVecBezier,CalcNormVecBezier
USE MOD_Particle_Vars,          ONLY:PDM,PartSpecies,PartState,LastPartPos,PEM
USE MOD_Particle_Mesh_Vars,     ONLY:PartBound
USE MOD_Particle_Surfaces_vars, ONLY:SideNormVec,SideType,epsilontol
!USE MOD_Particle_Surfaces_Vars, ONLY:BoundingBoxIsEmpty
USE MOD_Particle_Analyze,       ONLY:CalcEkinPart
USE MOD_Particle_Analyze_Vars,  ONLY:CalcPartBalance,nPartOut,PartEkinOut!,PartAnalyzeStep
USE MOD_TimeDisc_Vars,          ONLY:iter
USE MOD_Mesh_Vars,              ONLY:BC,nSides
USE MOD_Particle_Mesh_Vars,     ONLY:PartBCSideList
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
REAL                                 :: RanNum
INTEGER                              :: BCSideID
!===================================================================================================================================

IF (.NOT. ALLOCATED(PartBound%MapToPartBC)) THEN
  CALL abort(__STAMP__,&
  ' ERROR: PartBound not allocated!.',999,999.)
END IF

! Select the corresponding boundary condition and calculate particle treatment
SELECT CASE(PartBound%TargetBoundCond(PartBound%MapToPartBC(BC(SideID))))
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
  !BCSideID=PartBCSideList(SideID)
  PDM%ParticleInside(iPart) = .FALSE.
  alpha=-1.
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(2) !PartBound%ReflectiveBC)
!CASE(PartBound%ReflectiveBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL RANDOM_NUMBER(RanNum)
  BCSideID=PartBCSideList(SideID)
  IF(RanNum.GE.PartBound%MomentumACC(PartBound%MapToPartBC(BC(SideID)))) THEN
    ! perfectly reflection, specular re-emission
    CALL PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,BCSideID)
  ELSE
    CALL DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,BCSideID)
  END IF

!  SELECT CASE(SideType(SideID))
!  CASE(PLANAR)
!    n_loc=SideNormVec(1:3,SideID)
!  CASE(BILINEAR)
!    n_loc=CalcBiLinearNormVecBezier(xi,eta,SideID)
!  CASE(CURVED)
!    n_loc=CalcNormVecBezier(xi,eta,SideID)
!!    CALL abort(__STAMP__,'nvec for bezier not implemented!',999,999.)
!  END SELECT 
!  ! substract tolerance from length
!  LengthPartTrajectory=LengthPartTrajectory-epsilontol
!  ! intersection point with surface
!  LastPartPos(iPart,1:3) = LastPartPos(iPart,1:3) + PartTrajectory(1:3)*alpha*oneMinus
!
!  ! In vector notation: r_neu = r_alt + T - 2*((1-alpha)*<T,n>)*n
!  !v_aux = - 2*((1-alpha)*<T,n>)*n     (auxiliary variable, used twice)
!  !v_aux                  = -2*((1-alpha)*DOT_PRODUCT(PartTrajectory(1:3),n_loc))*n_loc
!  v_aux                  = -2*((LengthPartTrajectory-alpha)*DOT_PRODUCT(PartTrajectory(1:3),n_loc))*n_loc
!  !PartState(iPart,1:3)   = PartState(iPart,1:3)+PartTrajectory(1:3)+v_aux
!  PartState(iPart,1:3)   = PartState(iPart,1:3)+v_aux
!  ! new velocity vector 
!  !v_2=(1-alpha)*PartTrajectory(1:3)+v_aux
!  v_2=(LengthPartTrajectory-alpha)*PartTrajectory(1:3)+v_aux
!  PartState(iPart,4:6)   = SQRT(DOT_PRODUCT(PartState(iPart,4:6),PartState(iPart,4:6)))*&
!                           (1/(SQRT(DOT_PRODUCT(v_2,v_2))))*v_2                         +&
!                           PartBound%WallVelo(1:3,BC(SideID))
!  PartTrajectory=PartState(iPart,1:3) - LastPartPos(iPart,1:3)
!  lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
!                           +PartTrajectory(2)*PartTrajectory(2) &
!                           +PartTrajectory(3)*PartTrajectory(3) )
!  PartTrajectory=PartTrajectory/lengthPartTrajectory
!  lengthPartTrajectory=lengthPartTrajectory+epsilontol
!
!#if (PP_TimeDiscMethod==1) || (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6)
!  ! correction for Runge-Kutta (correct position!!)
!  !print*,'Pt_temp',Pt_temp(iPart,1:3)
!  ! get length of Pt_temp(iPart,1:3) || equals summed velocity change ! only exact for linear movement
!!  print*,'acceleration_old',Pt_temp(iPart,4:6)
!!  read*
!  absPt_temp=SQRT(Pt_temp(iPart,1)*Pt_temp(iPart,1)+Pt_temp(iPart,2)*Pt_temp(iPart,2)+Pt_temp(iPart,3)*Pt_temp(iPart,3))
!  ! scale PartTrajectory to new Pt_temp
!  Pt_temp(iPart,1:3)=absPt_temp*PartTrajectory(1:3)
!  ! deleate force history
!  Pt_temp(iPart,4:6)=0.
!  ! what happens with force term || acceleration?
!#endif 


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

!-----------------------------------------------------------------------------------------------------------------------------------
!CASE(PartBound%SymmetryBC)
CASE(10)
!-----------------------------------------------------------------------------------------------------------------------------------
  BCSideID=PartBCSideList(SideID)
  CALL PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,BCSideID)


CASE DEFAULT
  CALL abort(__STAMP__,&
' ERROR: PartBound not associated!. BC(SideID)',BC(SideID),REAL(SideID/nSides))
END SELECT !PartBound%MapToPartBC(BC(SideID)

END SUBROUTINE GetBoundaryInteractionRef


SUBROUTINE PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,SideID,BCSideID)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the perfect reflection in 3D
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Mesh_Vars,     ONLY:PartBound
USE MOD_Particle_Surfaces,      ONLY:CalcBiLinearNormVecBezier,CalcNormVecBezier
USE MOD_Particle_Vars,          ONLY:PartState,LastPartPos
USE MOD_Particle_Surfaces_vars, ONLY:SideNormVec,SideType,epsilontol
USE MOD_Particle_Mesh_Vars,     ONLY:epsInCell
USE MOD_TimeDisc_Vars,          ONLY:iter
USE MOD_Mesh_Vars,              ONLY:BC,nSides
USE MOD_Globals_Vars,    ONLY:EpsMach
#if (PP_TimeDiscMethod==1) || (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6)
USE MOD_Particle_Vars,          ONLY:Pt_temp,Pt
USE MOD_TimeDisc_Vars,          ONLY:RK_a,iStage
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES 
REAL,INTENT(INOUT)                :: PartTrajectory(1:3), lengthPartTrajectory, alpha
REAL,INTENT(IN)                   :: xi, eta
INTEGER,INTENT(IN)                :: PartID, SideID
INTEGER,INTENT(IN),OPTIONAL       :: BCSideID
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: v_2(1:3),v_aux(1:3),n_loc(1:3)
#if (PP_TimeDiscMethod==1) || (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6)
REAL                                 :: absPt_temp
#endif
!REAL,PARAMETER                       :: oneMinus=0.99999999
!REAL                                 :: oneMinus!=0.99999999
REAL                                  :: epsLength,epsReflect
!===================================================================================================================================

!OneMinus=1.0-MAX(epsInCell,epsilontol)
epsLength=MAX(epsInCell,epsilontol)*lengthPartTrajectory

IF(PRESENT(BCSideID))THEN
  SELECT CASE(SideType(BCSideID))
  CASE(PLANAR)
    n_loc=SideNormVec(1:3,BCSideID)
  CASE(BILINEAR)
    n_loc=CalcBiLinearNormVecBezier(xi,eta,BCSideID)
  CASE(CURVED)
    n_loc=CalcNormVecBezier(xi,eta,BCSideID)
  !   CALL abort(__STAMP__,'nvec for bezier not implemented!',999,999.)
  END SELECT 
ELSE
  SELECT CASE(SideType(SideID))
  CASE(PLANAR)
    n_loc=SideNormVec(1:3,SideID)
  CASE(BILINEAR)
    n_loc=CalcBiLinearNormVecBezier(xi,eta,SideID)
  CASE(CURVED)
    n_loc=CalcNormVecBezier(xi,eta,SideID)
  !   CALL abort(__STAMP__,'nvec for bezier not implemented!',999,999.)
  END SELECT 
END IF

IF(DOT_PRODUCT(PartTrajectory,n_loc).LE.0.) RETURN

LastPartPos(PartID,1:3) = LastPartPos(PartID,1:3) + PartTrajectory(1:3)*(alpha)

!IF(ALMOSTZERO(alpha))THEN
!  LastPartPos(PartID,1:3) = LastPartPos(PartID,1:3) + PartTrajectory(1:3)*(alpha)
!ELSE
!  LastPartPos(PartID,1:3) = LastPartPos(PartID,1:3) + PartTrajectory(1:3)*(alpha-epsLength)
!END IF

! In vector notation: r_neu = r_alt + T - 2*((1-alpha)*<T,n>)*n
!v_aux = - 2*((1-alpha)*<T,n>)*n     (auxiliary variable, used twice)
!v_aux                  = -2.0*((LengthPartTrajectory-alpha+epsLength)*DOT_PRODUCT(PartTrajectory(1:3),n_loc))*n_loc
v_aux                  = -2.0*((LengthPartTrajectory-alpha)*DOT_PRODUCT(PartTrajectory(1:3),n_loc))*n_loc
!IF(lengthPartTrajectory.LT.epsilontol)THEN
!  epsReflect=epsilontol
!ELSE
epsReflect=epsilontol*lengthPartTrajectory
!END IF

IF((DOT_PRODUCT(v_aux,v_aux)).GT.epsReflect)THEN
  PartState(PartID,1:3)   = PartState(PartID,1:3)+v_aux
  ! new velocity vector 
  !v_2=(1-alpha)*PartTrajectory(1:3)+v_aux
  v_2=(LengthPartTrajectory-alpha)*PartTrajectory(1:3)+v_aux
  PartState(PartID,4:6)   = SQRT(DOT_PRODUCT(PartState(PartID,4:6),PartState(PartID,4:6)))*&
                           (1/(SQRT(DOT_PRODUCT(v_2,v_2))))*v_2                         +&
                           PartBound%WallVelo(1:3,PartBound%MapToPartBC(BC(SideID)))
  
  PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
  lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                           +PartTrajectory(2)*PartTrajectory(2) &
                           +PartTrajectory(3)*PartTrajectory(3) )
  PartTrajectory=PartTrajectory/lengthPartTrajectory
  lengthPartTrajectory=lengthPartTrajectory!+epsilontol
  
#if (PP_TimeDiscMethod==1) || (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6)
  ! correction for Runge-Kutta (correct position!!)
  absPt_temp=SQRT(Pt_temp(PartID,1)*Pt_temp(PartID,1)+Pt_temp(PartID,2)*Pt_temp(PartID,2)+Pt_temp(PartID,3)*Pt_temp(PartID,3))
  ! scale PartTrajectory to new Pt_temp
  Pt_temp(PartID,1:3)=absPt_temp*PartTrajectory(1:3)
  ! deleate force history
  Pt_temp(PartID,4:6)=0.
  ! what happens with force term || acceleration?
#endif 
END IF

END SUBROUTINE PerfectReflection


SUBROUTINE DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,SideID,BCSideID)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the diffuse reflection in 3D
! only implemented for DoRefMapping tracking
! PartBCs are reduced!
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals,                ONLY:CROSS,abort
USE MOD_Globals_Vars,           ONLY:PI
USE MOD_Particle_Mesh_Vars,     ONLY:PartBound
USE MOD_Particle_Surfaces,      ONLY:CalcBiLinearNormAndTang,CalcNormAndTangBezier
USE MOD_Particle_Vars,          ONLY:PartState,LastPartPos,time,Species,BoltzmannConst,PartSpecies
USE MOD_DSMC_Vars,              ONLY:SpecDSMC,CollisMode
USE MOD_Particle_Surfaces_vars, ONLY:SideNormVec,SideType,epsilontol,BezierControlPoints3D
USE MOD_TimeDisc_Vars,          ONLY:iter
USE MOD_Mesh_Vars,              ONLY:BC,nSides,NGEO
USE MOD_DSMC_Vars,              ONLY:PartStateIntEn,SpecDSMC, DSMC, SampWall, SurfMesh, useDSMC, CollisMode
USE MOD_Particle_Mesh_Vars,     ONLY:epsInCell
#if (PP_TimeDiscMethod==1) || (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6)
USE MOD_Particle_Vars,          ONLY:Pt_temp,Pt
USE MOD_TimeDisc_Vars,          ONLY:RK_a,iStage
#endif
USE MOD_TImeDisc_Vars,          ONLY:dt,tend
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES 
REAL,INTENT(INOUT)                :: PartTrajectory(1:3), lengthPartTrajectory, alpha
REAL,INTENT(IN)                   :: xi, eta
INTEGER,INTENT(IN)                :: PartID, SideID
INTEGER,INTENT(IN),OPTIONAL       :: BCSideID
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: locBCID, vibQuant, vibQuantNew, VibQuantWall
REAL                                 :: v_2(1:3),v_aux(1:3)
REAL                                 :: VibQuantNewR                                                !
#if (PP_TimeDiscMethod==1) || (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6)
REAL                                 :: absPt_temp
#endif
!REAL,PARAMETER                       :: oneMinus=0.99999999
REAL                                 :: oneMinus                
REAL                                 :: VeloReal, RanNum, EtraOld, VeloCrad, Fak_D
REAL                                 :: EtraWall, EtraNew
REAL                                 :: WallVelo(1:3), WallTemp, TransACC, VibACC, RotACC
REAL                                 :: n_loc(1:3), tang1(1:3),tang2(1:3), NewVelo(3)
REAL                                 :: ErotNew, ErotWall, EVibNew, Phi, Cmr, VeloCx, VeloCy, VeloCz
REAL                                 :: WallTransACC
!===================================================================================================================================

OneMinus=1.0-epsInCell

! additional states
locBCID=PartBound%MapToPartBC(BC(SideID))
! get BC values
WallVelo     = PartBound%WallVelo(1:3,locBCID)
WallTemp     = PartBound%WallTemp(locBCID)
WallTransACC = PartBound%TransACC(locBCID)
VibACC       = PartBound%VibACC(locBCID)
RotACC       = PartBound%RotACC(locBCID)

IF(PRESENT(BCSideID))THEN
  SELECT CASE(SideType(BCSideID))
  CASE(PLANAR)
    n_loc=SideNormVec(1:3,BCSideID)
    tang1=BezierControlPoints3D(:,NGeo,0,BCSideID)-BezierControlPoints3D(:,0,0,BCSideID)
    tang2=BezierControlPoints3D(:,0,NGeo,BCSideID)-BezierControlPoints3D(:,0,0,BCSideID)
  CASE(BILINEAR)
    CALL CalcBiLinearNormAndTang(n_loc,tang1,tang2,xi,eta,BCSideID)
  CASE(CURVED)
    CALL CalcNormAndTangBezier(n_loc,tang1,tang2,xi,eta,BCSideID)
  !   CALL abort(__STAMP__,'nvec for bezier not implemented!',999,999.)
  END SELECT 
ELSE
  SELECT CASE(SideType(SideID))
  CASE(PLANAR)
    n_loc=SideNormVec(1:3,SideID)
    tang1=BezierControlPoints3D(:,NGeo,0,SideID)-BezierControlPoints3D(:,0,0,SideID)
    tang2=BezierControlPoints3D(:,0,NGeo,SideID)-BezierControlPoints3D(:,0,0,SideID)
  CASE(BILINEAR)
    CALL CalcBiLinearNormAndTang(n_loc,tang1,tang2,xi,eta,SideID)
  CASE(CURVED)
    CALL CalcNormAndTangBezier(n_loc,tang1,tang2,xi,eta,SideID)
  !   CALL abort(__STAMP__,'nvec for bezier not implemented!',999,999.)
  END SELECT 
END IF

!IF(DOT_PRODUCT(n_loc,PartTrajectory).LT.0.)   CALL abort(__STAMP__,&
!  ' parttrajectory inverse to n_loc')

! substract tolerance from length
!LengthPartTrajectory=LengthPartTrajectory!-epsilontol
! intersection point with surface
LastPartPos(PartID,1:3) = LastPartPos(PartID,1:3) + PartTrajectory(1:3)*alpha*oneMinus

! calculate new velocity vector (Extended Maxwellian Model)
VeloReal = SQRT(PartState(PartID,4) * PartState(PartID,4) + &
                PartState(PartID,5) * PartState(PartID,5) + &
                PartState(PartID,6) * PartState(PartID,6))

EtraOld     = 0.5 * Species(PartSpecies(PartID))%MassIC * VeloReal**2 
CALL RANDOM_NUMBER(RanNum)
VeloCrad    = SQRT(-LOG(RanNum))
CALL RANDOM_NUMBER(RanNum)
VeloCz      = SQRT(-LOG(RanNum))
Fak_D       = VeloCrad**2 + VeloCz**2

EtraWall    = BoltzmannConst * WallTemp * Fak_D
EtraNew     = EtraOld + TransACC * (EtraWall - EtraOld)
Cmr         = SQRT(2.0 * EtraNew / (Species(PartSpecies(PartID))%MassIC * Fak_D))
CALL RANDOM_NUMBER(RanNum)
Phi     = 2.0 * PI * RanNum
VeloCx  = Cmr * VeloCrad * COS(Phi)
VeloCy  = Cmr * VeloCrad * SIN(Phi)
VeloCz  = Cmr * VeloCz

! IF ((DSMC%CalcSurfaceVal.AND.(Time.ge.(1-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroValues)) THEN
! !----  Sampling for energy (translation) accommodation at walls
! ! has to be corrected to new scheme, at output a mpi_reduce has to be called
! !  SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(1) = &
! !    SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(1) + &
! !    EtraOld * Species(PartSpecies(i))%MacroParticleFactor  
! !  SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(2) = &
! !    SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(2) + &
! !    EtraWall * Species(PartSpecies(i))%MacroParticleFactor 
! !  SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(3) = &
! !    SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(3) + & 
! !    EtraNew * Species(PartSpecies(i))%MacroParticleFactor  
! END IF
 
!   Transformation local distribution -> global coordinates
! from flux comutaion
! v = nv*u+t1*v+t2*f3

! NewVelo(1) = tang1(1)*VeloCx + (n_loc(3)*tang1(2)-n_loc(2)*tang1(3))*VeloCy - n_loc(1)*VeloCz
! NewVelo(2) = tang1(2)*VeloCx + (n_loc(1)*tang1(3)-n_loc(3)*tang1(1))*VeloCy - n_loc(2)*VeloCz
! NewVelo(3) = tang1(3)*VeloCx + (n_loc(2)*tang1(1)-n_loc(1)*tang1(2))*VeloCy - n_loc(3)*VeloCz
NewVelo = VeloCx*tang1+CROSS(n_loc,tang1)*VeloCy-VeloCz*n_loc

PartState(PartID,1:3)   = LastPartPos(PartID,1:3) + (1. - alpha) * dt * NewVelo(1:3)

! Internal energy accommodation
IF (useDSMC) THEN
  IF (CollisMode.GT.1) THEN
    IF (SpecDSMC(PartSpecies(PartID))%InterID.EQ.2) THEN

    ! Rotational energy accommodation
    
      CALL RANDOM_NUMBER(RanNum)
      ErotWall = - BoltzmannConst * WallTemp * LOG(RanNum)
      ErotNew  = PartStateIntEn(PartID,2) + RotACC *(ErotWall - PartStateIntEn(PartID,2))
    
!       IF ((DSMC%CalcSurfaceVal.AND.(Time.ge.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroValues)) THEN
! !    !----  Sampling for internal energy accommodation at walls
! !    
! !        SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(4) = &
! !          SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(4) &
! !          + PartStateIntEn(i,2) * Species(PartSpecies(i))%MacroParticleFactor   
! !        SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(5) = &
! !          SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(5) &
! !          + ErotWall * Species(PartSpecies(i))%MacroParticleFactor  
! !        SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(6) = &
! !          SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(6) &
! !          + ErotNew * Species(PartSpecies(i))%MacroParticleFactor 
!       END IF 
    
      PartStateIntEn(PartID,2) = ErotNew
    
    !---- Vibrational energy accommodation
    
      VibQuant     = NINT(PartStateIntEn(PartID,1)/(BoltzmannConst*SpecDSMC(PartSpecies(PartID))%CharaTVib) &
                   - DSMC%GammaQuant)
      CALL RANDOM_NUMBER(RanNum)
      VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(PartSpecies(PartID))%CharaTVib)
      DO WHILE (VibQuantWall.GE.SpecDSMC(PartSpecies(PartID))%MaxVibQuant)
        CALL RANDOM_NUMBER(RanNum)
        VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(PartSpecies(PartID))%CharaTVib)
      END DO
      VibQuantNewR = VibQuant + VibACC*(VibQuantWall - VibQuant)
      VibQuantNew = INT(VibQuantNewR)
      CALL RANDOM_NUMBER(RanNum)
      IF (RanNum.LT.(VibQuantNewR - VibQuantNew)) THEN
        EvibNew = (VibQuantNew + DSMC%GammaQuant + 1.0d0)*BoltzmannConst*SpecDSMC(PartSpecies(PartID))%CharaTVib
      ELSE
        EvibNew = (VibQuantNew + DSMC%GammaQuant)*BoltzmannConst*SpecDSMC(PartSpecies(PartID))%CharaTVib
      END IF
    
!       IF ((DSMC%CalcSurfaceVal.AND.(Time.ge.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroValues)) THEN
! !     !----  Sampling for internal energy accommodation at walls
! !     
! !         SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(7) = &
! !           SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(7) & 
! !                             + (VibQuant + DSMC%GammaQuant) &
! !                             * BoltzmannConst * SpecDSMC(PartSpecies(i))%CharaTVib * Species(PartSpecies(i))%MacroParticleFactor
! !         SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(8) = &
! !           SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(8) + VibQuantWall &
! !                             * BoltzmannConst * SpecDSMC(PartSpecies(i))%CharaTVib * Species(PartSpecies(i))%MacroParticleFactor
! !         SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(9) = &
! !           SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Energy(9) &
! !                             + EvibNew * Species(PartSpecies(i))%MacroParticleFactor
!       END IF
    
      PartStateIntEn(PartID,1) = EvibNew
    END IF
  END IF
END IF

! IF ((DSMC%CalcSurfaceVal.AND.(Time.ge.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroValues)) THEN
! ! !----  Sampling force at walls
! !   SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Force(1) = &
! !     SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Force(1) + &
! !     Species(PartSpecies(i))%MassIC * (PartState(i,4) - VelX) * Species(PartSpecies(i))%MacroParticleFactor
! !   SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Force(2) = &
! !     SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Force(2) + &
! !     Species(PartSpecies(i))%MassIC * (PartState(i,5) - VelY) * Species(PartSpecies(i))%MacroParticleFactor
! !   SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Force(3) = &
! !     SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Force(3) + &
! !     Species(PartSpecies(i))%MassIC * (PartState(i,6) - VelZ) * Species(PartSpecies(i))%MacroParticleFactor
! ! !---- Counter for collisions (normal wall collisions - not to count if only SpeciesSwaps to be counted)
! !   IF (.NOT.DSMC%CalcSurfCollis_OnlySwaps) THEN
! !     SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Counter(PartSpecies(i)) = &
! !       SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(MPIGEO%ElemToSide(1,iLocSide,Element)))%Counter(PartSpecies(i)) + 1
! !   END IF
! END IF

!----  saving new particle velocity
PartState(PartID,4:6)   = NewVelo(1:3) + WallVelo(1:3)

END SUBROUTINE DiffuseReflection


END MODULE MOD_Particle_Boundary_Condition
