#include "boltzplatz.h"

MODULE  MOD_PICInterpolation
!===================================================================================================================================
! 
!===================================================================================================================================
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: InterpolateFieldToParticle,InitializeInterpolation,InterpolateFieldToSingleParticle,InterpolateVariableExternalField
!===================================================================================================================================
INTERFACE InitializeInterpolation
  MODULE PROCEDURE InitializeInterpolation
END INTERFACE

INTERFACE InterpolateFieldToParticle
  MODULE PROCEDURE InterpolateFieldToParticle
END INTERFACE

INTERFACE InterpolateVariableExternalField
  MODULE PROCEDURE InterpolateVariableExternalField
END INTERFACE

INTERFACE InterpolateFieldToSingleParticle
  MODULE PROCEDURE InterpolateFieldToSingleParticle
END INTERFACE
!===================================================================================================================================

CONTAINS

SUBROUTINE InitializeInterpolation
!===================================================================================================================================
! Initialize the interpolation variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc,                ONLY:PP_nElems
USE MOD_ReadInTools
USE MOD_Particle_Vars,          ONLY : PDM
USE MOD_PICInterpolation_Vars
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: ALLOCSTAT
!REAL                      :: P(3,8), T(3,3), T_inv(3,3)
!INTEGER                   :: iNode, Elem
REAL                      :: scaleExternalField
!===================================================================================================================================
InterpolationType = GETSTR('PIC-Interpolation-Type','particle_position')
InterpolationElemLoop = GETLOGICAL('PIC-InterpolationElemLoop','.TRUE.')
IF (InterpolationElemLoop) THEN !If user-defined F: F for all procs
  IF (PP_nElems.GT.10) THEN !so far arbitrary threshold...
    InterpolationElemLoop=.FALSE. !switch off for procs with high number of Elems
  END IF
END IF
externalField(1:6)= GETREALARRAY('PIC-externalField',6,'0.,0.,0.,0.,0.,0.')
scaleexternalField= GETREAL('PIC-scaleexternalField','1.0')
externalField=externalField*ScaleExternalField
!SWRITE(*,*) " External fied", externalfield
DoInterpolation   = GETLOGICAL('PIC-DoInterpolation','.TRUE.')
useBGField        = GETLOGICAL('PIC-BG-Field','.FALSE.')

! Variable external field 
useVariableExternalField = .FALSE.
FileNameVariableExternalField=GETSTR('PIC-curvedexternalField','none')
IF (FileNameVariableExternalField.EQ.'none') THEN
  FileNameVariableExternalField=GETSTR('PIC-variableexternalField','none')
END IF
IF (FileNameVariableExternalField.NE.'none') THEN
  useVariableExternalField = .TRUE.
  CALL ReadVariableExternalField()
END IF

!--- Allocate arrays for interpolation of fields to particles
SDEALLOCATE(FieldAtParticle)
ALLOCATE(FieldAtParticle(1:PDM%maxParticleNumber,1:6), STAT=ALLOCSTAT) 
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
  __STAMP__ &
  ,'ERROR in pic_interpolation.f90: Cannot allocate FieldAtParticle array!',ALLOCSTAT)
END IF

SELECT CASE(TRIM(InterpolationType))
CASE('nearest_blurycenter')
   InterpolationType='nearest_blurrycenter'
CASE('nearest_blurrycenter')
CASE('particle_position_slow')
CASE('particle_position')
CASE('nearest_gausspoint')
CASE DEFAULT
  CALL abort(&
  __STAMP__ &
  ,'Unknown InterpolationType in pic_init.f90')
END SELECT
END SUBROUTINE InitializeInterpolation


SUBROUTINE InterpolateFieldToParticle(doInnerParts)
!===================================================================================================================================
! interpolates field to particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Vars,           ONLY:PartPosRef,PDM,PartState,PEM,PartPosGauss
USE MOD_Particle_Tracking_Vars,  ONLY:DoRefMapping
#ifndef PP_HDG
USE MOD_DG_Vars,                 ONLY:U
#endif
USE MOD_PIC_Vars!,      ONLY: 
USE MOD_PICInterpolation_Vars,   ONLY:useVariableExternalField,FieldAtParticle,externalField,DoInterpolation,InterpolationType
USE MOD_PICInterpolation_Vars,   ONLY:InterpolationElemLoop
USE MOD_PICDepo_Vars,            ONLY:DepositionType,GaussBorder
USE MOD_Eval_xyz,                ONLY:Eval_xyz_elemcheck,Eval_XYZ_Curved,Eval_xyz_Part2
#ifdef PP_POIS
USE MOD_Equation_Vars,           ONLY:E
#endif
#ifdef PP_HDG
#if PP_nVar==1
USE MOD_Equation_Vars,        ONLY:E
#elif PP_nVar==3
USE MOD_Equation_Vars,        ONLY:B
#else
USE MOD_Equation_Vars,        ONLY:B,E
#endif /*PP_nVar==1*/
#endif /*PP_HDG*/
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
USE MOD_Particle_Vars,        ONLY:DoSurfaceFlux
#endif /*HDG-LSERK*/
#ifdef MPI
! only required for shape function??
USE MOD_Particle_MPI_Vars,    ONLY:PartMPIExchange
#endif 

!----------------------------------------------------------------------------------------------------------------------------------
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL                          :: doInnerParts
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                                    
INTEGER                          :: firstPart,lastPart
REAL                             :: Pos(3)                                                      
REAL                             :: field(6)                                                    
INTEGER                          :: iPart,iElem
! for Nearest GaussPoint
INTEGER                          :: a,b,k,ii,l,m
#if defined PP_POIS || (defined PP_HDG && PP_nVar==4)
REAL                             :: HelperU(1:6,0:PP_N,0:PP_N,0:PP_N)
#endif /*(PP_POIS||PP_HDG)*/
LOGICAL                          :: NotMappedSurfFluxParts
!===================================================================================================================================
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
NotMappedSurfFluxParts=DoSurfaceFlux !Surfaceflux particles inserted before interpolation and tracking. Field at wall is needed!
#else
NotMappedSurfFluxParts=.FALSE.
#endif /*HDG-LSERK*/
! null field vector
field=0.

IF(doInnerParts)THEN
  firstPart=1
  lastPart =PDM%ParticleVecLength
ELSE
#ifdef MPI
  firstPart=PDM%ParticleVecLength-PartMPIExchange%nMPIParticles+1
  lastPart =PDM%ParticleVecLength
#else
  firstPart=2
  LastPart =1
#endif /*MPI*/
END IF
! thats wrong
IF(firstPart.GT.lastPart) RETURN

IF (.NOT.InterpolationElemLoop) THEN
  DO iPart = firstPart, LastPart
    IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
    CALL InterpolateFieldToSingleParticle(iPart,FieldAtparticle(iPart,1:6))
  END DO
  RETURN
END IF

IF(useVariableExternalField) THEN ! used curved external Bz
  FieldAtParticle(firstPart:lastPart,:) = 0.
  FieldAtParticle(firstPart:lastPart,1) = externalField(1)
  FieldAtParticle(firstPart:lastPart,2) = externalField(2)
  FieldAtParticle(firstPart:lastPart,3) = externalField(3)
#if (PP_nVar==8)
  FieldAtParticle(firstPart:lastPart,4) = externalField(4)
  FieldAtParticle(firstPart:lastPart,5) = externalField(5)
#endif
  ! Bz field strength at particle position
  DO iPart = firstPart, LastPart
    FieldAtparticle(iPart,6) = InterpolateVariableExternalField(PartState(iPart,3))
  END DO
ELSE ! useVariableExternalField
  FieldAtParticle(firstPart:lastPart,:) = 0.
  FieldAtParticle(firstPart:lastPart,1) = externalField(1)
  FieldAtParticle(firstPart:lastPart,2) = externalField(2)
  FieldAtParticle(firstPart:lastPart,3) = externalField(3)
#if (PP_nVar==8)
  FieldAtParticle(firstPart:lastPart,4) = externalField(4)
  FieldAtParticle(firstPart:lastPart,5) = externalField(5)
  FieldAtParticle(firstPart:lastPart,6) = externalField(6)
#endif
END IF ! use constant external field

IF (DoInterpolation) THEN                 ! skip if no self fields are calculated
  SELECT CASE(TRIM(InterpolationType))
  CASE('nearest_blurrycenter')
    ! add fields to fields at particle position
    ! should nearest_blurrycenter not require to
    ! a) evaluate the polynomial at Xi=0??
    ! b) require the mean value of the field
    !m = INT(PP_N/2)+1
    DO iElem=1,PP_nElems
#if (PP_nVar==8)
!
#ifdef PP_POIS
      HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
      HelperU(4:6,:,:,:) = U(4:6,:,:,:,iElem)
      CALL Eval_xyz_Part2((/0.,0.,0./),6,PP_N,HelperU,field(1:6),iElem)
#else
      CALL Eval_xyz_Part2((/0.,0.,0./),6,PP_N,U(1:6,:,:,:,iElem),field(1:6),iElem)
#endif /*PP_POIS*/
!
#else /*PP_nVar not 8*/
!
#ifdef PP_POIS
      CALL Eval_xyz_Part2((/0.,0.,0./),3,PP_N,E(1:3,:,:,:,iElem),field(1:3),iElem)
#elif defined PP_HDG
#if PP_nVar==1
      CALL Eval_xyz_Part2((/0.,0.,0./),3,PP_N,E(1:3,:,:,:,iElem),field(1:3),iElem)
#elif PP_nVar==3
      CALL Eval_xyz_Part2((/0.,0.,0./),3,PP_N,B(1:3,:,:,:,iElem),field(4:6),iElem)
#else
      HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
      HelperU(4:6,:,:,:) = B(1:3,:,:,:,iElem)
      CALL Eval_xyz_Part2((/0.,0.,0./),6,PP_N,HelperU,field(1:6),iElem)
#endif /*PP_nVar==1*/
#else /*not HDG and not POIS and PP_nVar not 8*/
      CALL Eval_xyz_Part2((/0.,0.,0./),3,PP_N,U(1:3,:,:,:,iElem),field(1:3),iElem)
#endif /*PP_POIS*/
!
#endif /*(PP_nVar==8)*/
      DO iPart=firstPart,LastPart
        IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
        IF(PEM%Element(iPart).EQ.iElem)THEN
          FieldAtParticle(iPart,:) = FieldAtParticle(iPart,:) + field(1:6)
        END IF! Element(iPart).EQ.iElem
      END DO ! iPart
    END DO ! iElem
  CASE('particle_position_slow')
    DO iElem=1,PP_nElems
      DO iPart = firstPart, LastPart
        IF(.NOT.PDM%ParticleInside(iPart))CYCLE
        IF(PEM%Element(iPart).EQ.iElem)THEN
          Pos = PartState(iPart,1:3)
          !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
          HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
          HelperU(4:6,:,:,:) = U(4:6,:,:,:,iElem)
          CALL eval_xyz_curved(Pos,6,PP_N,HelperU,field(1:6),iElem)
#else
          CALL eval_xyz_curved(Pos,6,PP_N,U(1:6,:,:,:,iElem),field(1:6),iElem)
#endif /*PP_POIS*/
#else
#ifdef PP_POIS
          CALL eval_xyz_curved(Pos,3,PP_N,E(1:3,:,:,:,iElem),field(1:3),iElem)
#elif defined PP_HDG
#if PP_nVar==1
          CALL eval_xyz_curved(Pos,3,PP_N,E(1:3,:,:,:,iElem),field(1:3),iElem)
#elif PP_nVar==3
          CALL eval_xyz_curved(Pos,3,PP_N,B(1:3,:,:,:,iElem),field(4:6),iElem)
#else
          HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
          HelperU(4:6,:,:,:) = B(1:3,:,:,:,iElem)
          CALL eval_xyz_curved(Pos,6,PP_N,HelperU,field(1:6),iElem)
#endif /*PP_nVar*/
#else
          CALL eval_xyz_curved(Pos,3,PP_N,U(1:3,:,:,:,iElem),field(1:3),iElem)
#endif /*PP_POIS*/
#endif /*(PP_nVar==8)*/
          FieldAtParticle(iPart,:) = FieldAtParticle(iPart,:) + field(1:6)
        END IF ! Element(iPart).EQ.iElem
      END DO ! iPart
    END DO ! iElem=1,PP_nElems
  CASE('particle_position')
    IF(.NOT.NotMappedSurfFluxParts .AND.(DoRefMapping .OR. TRIM(DepositionType).EQ.'nearest_gausspoint'))THEN
      ! particles have already been mapped in deposition, other eval routine used
      DO iElem=1,PP_nElems
        DO iPart=firstPart,LastPart
          IF(.NOT.PDM%ParticleInside(iPart))CYCLE
          IF(PEM%Element(iPart).EQ.iElem)THEN
            IF(.NOT.DoRefMapping)THEN
              CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),iElem)
            END IF
            !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
            HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
            HelperU(4:6,:,:,:) = U(4:6,:,:,:,iElem)
            CALL eval_xyz_part2(PartPosRef(1:3,iPart),6,PP_N,HelperU,field(1:6),iElem)
#else
            CALL eval_xyz_part2(PartPosRef(1:3,iPart),6,PP_N,U(1:6,:,:,:,iElem),field(1:6),iElem)
#endif
#else
#ifdef PP_POIS
            CALL eval_xyz_part2(PartPosRef(1:3,iPart),3,PP_N,E(1:3,:,:,:,iElem),field(1:3),iElem)     
#elif defined PP_HDG
#if PP_nVar==1
            CALL eval_xyz_part2(PartPosRef(1:3,iPart),3,PP_N,E(1:3,:,:,:,iElem),field(1:3),iElem)     
#elif PP_nVar==3
            CALL eval_xyz_part2(PartPosRef(1:3,iPart),3,PP_N,B(1:3,:,:,:,iElem),field(4:6),iElem)     
#else
            HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
            HelperU(4:6,:,:,:) = B(1:3,:,:,:,iElem)
            CALL eval_xyz_part2(PartPosRef(1:3,iPart),6,PP_N,HelperU,field(1:6),iElem)
#endif
#else
            CALL eval_xyz_part2(PartPosRef(1:3,iPart),3,PP_N,U(1:3,:,:,:,iElem),field(1:3),iElem)
#endif
#endif
            FieldAtParticle(iPart,:) = FieldAtParticle(iPart,:) + field(1:6)
          END IF ! Element(iPart).EQ.iElem
        END DO ! iPart
      END DO ! iElem=1,PP_nElems
    ELSE IF(NotMappedSurfFluxParts .AND.(DoRefMapping .OR. TRIM(DepositionType).EQ.'nearest_gausspoint'))THEN
      !some particle are mapped, surfaceflux particles (dtFracPush) are not
      DO iElem=1,PP_nElems
        DO iPart=firstPart,LastPart
          IF(.NOT.PDM%ParticleInside(iPart))CYCLE
          IF(PEM%Element(iPart).EQ.iElem)THEN
            IF(PDM%dtFracPush(iPart))THEN ! same as in "particles are not yet mapped"
              Pos = PartState(iPart,1:3)
              !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
              HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
              HelperU(4:6,:,:,:) = U(4:6,:,:,:,iElem)
              CALL eval_xyz_curved(Pos,6,PP_N,HelperU,field(1:6),iElem,iPart)
#else
              CALL eval_xyz_curved(Pos,6,PP_N,U(1:6,:,:,:,iElem),field(1:6),iElem,iPart)
#endif
#else
#ifdef PP_POIS
              CALL eval_xyz_curved(Pos,3,PP_N,E(1:3,:,:,:,iElem),field(1:3),iElem,iPart)
#elif defined PP_HDG
#if PP_nVar==1
              CALL eval_xyz_curved(Pos,3,PP_N,E(1:3,:,:,:,iElem),field(1:3),iElem,iPart)
#elif PP_nVar==3
              CALL eval_xyz_curved(Pos,3,PP_N,B(1:3,:,:,:,iElem),field(4:6),iElem,iPart)
#else
              HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
              HelperU(4:6,:,:,:) = B(1:3,:,:,:,iElem)
              CALL eval_xyz_curved(Pos,6,PP_N,HelperU,field(1:6),iElem,iPart)
#endif
#else
              CALL eval_xyz_curved(Pos,3,PP_N,U(1:3,:,:,:,iElem),field(1:3),iElem,iPart)
#endif
#endif
            ELSE !.NOT.PDM%dtFracPush(iPart): same as in "particles have already been mapped in deposition, other eval routine used"
              IF(.NOT.DoRefMapping)THEN
                CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),iElem)
              END IF
              !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
              HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
              HelperU(4:6,:,:,:) = U(4:6,:,:,:,iElem)
              CALL eval_xyz_part2(PartPosRef(1:3,iPart),6,PP_N,HelperU,field(1:6),iElem)
#else
              CALL eval_xyz_part2(PartPosRef(1:3,iPart),6,PP_N,U(1:6,:,:,:,iElem),field(1:6),iElem)
#endif
#else
#ifdef PP_POIS
              CALL eval_xyz_part2(PartPosRef(1:3,iPart),3,PP_N,E(1:3,:,:,:,iElem),field(1:3),iElem)
#elif defined PP_HDG
#if PP_nVar==1
              CALL eval_xyz_part2(PartPosRef(1:3,iPart),3,PP_N,E(1:3,:,:,:,iElem),field(1:3),iElem)
#elif PP_nVar==3
              CALL eval_xyz_part2(PartPosRef(1:3,iPart),3,PP_N,B(1:3,:,:,:,iElem),field(4:6),iElem)
#else
              HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
              HelperU(4:6,:,:,:) = B(1:3,:,:,:,iElem)
              CALL eval_xyz_part2(PartPosRef(1:3,iPart),6,PP_N,HelperU,field(1:6),iElem)
#endif
#else
              CALL eval_xyz_part2(PartPosRef(1:3,iPart),3,PP_N,U(1:3,:,:,:,iElem),field(1:3),iElem)
#endif
#endif
            END IF !PDM%dtFracPush(iPart)
            FieldAtParticle(iPart,:) = FieldAtParticle(iPart,:) + field(1:6)
          END IF ! Element(iPart).EQ.iElem
        END DO ! iPart
      END DO ! iElem=1,PP_nElems
    ELSE ! particles are not yet mapped
      DO iElem=1,PP_nElems
        DO iPart=firstPart,LastPart
          IF(.NOT.PDM%ParticleInside(iPart))CYCLE
          IF(PEM%Element(iPart).EQ.iElem)THEN
            Pos = PartState(iPart,1:3)
            !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
            HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
            HelperU(4:6,:,:,:) = U(4:6,:,:,:,iElem)
            CALL eval_xyz_curved(Pos,6,PP_N,HelperU,field(1:6),iElem,iPart)
#else
            CALL eval_xyz_curved(Pos,6,PP_N,U(1:6,:,:,:,iElem),field(1:6),iElem,iPart)
#endif
#else
#ifdef PP_POIS
            CALL eval_xyz_curved(Pos,3,PP_N,E(1:3,:,:,:,iElem),field(1:3),iElem,iPart)
#elif defined PP_HDG
#if PP_nVar==1
            CALL eval_xyz_curved(Pos,3,PP_N,E(1:3,:,:,:,iElem),field(1:3),iElem,iPart)
#elif PP_nVar==3
            CALL eval_xyz_curved(Pos,3,PP_N,B(1:3,:,:,:,iElem),field(4:6),iElem,iPart)
#else
            HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
            HelperU(4:6,:,:,:) = B(1:3,:,:,:,iElem)
            CALL eval_xyz_curved(Pos,6,PP_N,HelperU,field(1:6),iElem,iPart)
#endif
#else
            CALL eval_xyz_curved(Pos,3,PP_N,U(1:3,:,:,:,iElem),field(1:3),iElem,iPart)
#endif         
#endif
            FieldAtParticle(iPart,:) = FieldAtParticle(iPart,:) + field(1:6)
          END IF ! Element(iPart).EQ.iElem
        END DO ! iPart
      END DO ! iElem=1,PP_nElems
    END IF ! DoRefMapping .or. Depositiontype=nearest_gausspoint
  CASE('nearest_gausspoint')
    ! particles have already been mapped in deposition
    IF(MOD(PP_N,2).EQ.0) THEN
      a = PP_N/2
      b = a
    ELSE
      a = (PP_N+1)/2
      b = a-1
    END IF
    DO iElem=1,PP_nElems
      DO iPart=firstPart,LastPart
        IF(.NOT.PDM%ParticleInside(iPart))CYCLE
        IF(PEM%Element(iPart).EQ.iElem)THEN
          IF(.NOT.DoRefMapping .OR. (NotMappedSurfFluxParts .AND. PDM%dtFracPush(iPart)))THEN
            CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),iElem)
          END IF
          ! compute exact k,l,m
          !! x-direction
          k = a
          DO ii = 0,b-1
            IF(ABS(PartPosRef(1,iPart)).GE.GaussBorder(PP_N-ii))THEN
              k = PP_N-ii
              EXIT
            END IF
          END DO
          k = NINT((PP_N+SIGN(2.0*k-PP_N,PartPosRef(1,iPart)))/2)
          !! y-direction
          l = a
          DO ii = 0,b-1
            IF(ABS(PartPosRef(2,iPart)).GE.GaussBorder(PP_N-ii))THEN
              l = PP_N-ii
              EXIT
            END IF
          END DO
          l = NINT((PP_N+SIGN(2.0*l-PP_N,PartPosRef(2,iPart)))/2)
          !! z-direction
          m = a
          DO ii = 0,b-1
            IF(ABS(PartPosRef(3,iPart)).GE.GaussBorder(PP_N-ii))THEN
              m = PP_N-ii
              EXIT
            END IF
          END DO
          m = NINT((PP_N+SIGN(2.0*m-PP_N,PartPosRef(3,iPart)))/2)
          PartPosGauss(iPart,1) = k
          PartPosGauss(iPart,2) = l
          PartPosGauss(iPart,3) = m
          !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
         field(1:3) = E(1:3,PartPosGauss(iPart,1),PartPosGauss(iPart,2),PartPosGauss(iPart,3), iElem)
         field(4:6) = U(4:6,PartPosGauss(iPart,1),PartPosGauss(iPart,2),PartPosGauss(iPart,3), iElem)
         FieldAtParticle(iPart,:) = FieldAtParticle(iPart,:) + field(1:6)
#else
         field = U(1:6,PartPosGauss(iPart,1),PartPosGauss(iPart,2),PartPosGauss(iPart,3), iElem)
         FieldAtParticle(iPart,:) = FieldAtParticle(iPart,:) + field(1:6)
#endif
#else
#ifdef PP_POIS
         field(1:3) = E(1:3,PartPosGauss(iPart,1),PartPosGauss(iPart,2),PartPosGauss(iPart,3), iElem)
         FieldAtParticle(iPart,1:3) = FieldAtParticle(iPart,1:3) + field(1:3)
#elif defined PP_HDG
#if PP_nVar==1
          field(1:3) = E(1:3,PartPosGauss(iPart,1),PartPosGauss(iPart,2),PartPosGauss(iPart,3), iElem)
          FieldAtParticle(iPart,1:3) = FieldAtParticle(iPart,1:3) + field(1:3)
#elif PP_nVar==3
          field(4:6) = B(1:3,PartPosGauss(iPart,1),PartPosGauss(iPart,2),PartPosGauss(iPart,3), iElem)
          FieldAtParticle(iPart,4:6) = FieldAtParticle(iPart,4:6) + field(4:6)
#else
          field(1:3) = E(1:3,PartPosGauss(iPart,1),PartPosGauss(iPart,2),PartPosGauss(iPart,3), iElem)
          field(4:6) = B(1:3,PartPosGauss(iPart,1),PartPosGauss(iPart,2),PartPosGauss(iPart,3), iElem)
          FieldAtParticle(iPart,:) = FieldAtParticle(iPart,:) + field(1:6)
#endif
#else
         field(1:3) = U(1:3,PartPosGauss(iPart,1),PartPosGauss(iPart,2),PartPosGauss(iPart,3), iElem)
         FieldAtParticle(iPart,1:3) = FieldAtParticle(iPart,1:3) + field(1:3)
#endif
#endif
        END IF ! Element(iPart).EQ.iElem
      END DO ! iPart
    END DO ! iElem=1,PP_nElems
  CASE DEFAULT
    CALL abort(&
__STAMP__&
       , 'ERROR: Unknown InterpolationType!')
  END SELECT
END IF
    
RETURN
END SUBROUTINE InterpolateFieldToParticle


SUBROUTINE InterpolateFieldToSingleParticle(PartID,FieldAtParticle)
!===================================================================================================================================
! interpolates field to particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Vars,           ONLY:PartPosRef,PDM,PartState,PEM,PartPosGauss
USE MOD_Particle_Tracking_Vars,  ONLY:DoRefMapping
#ifndef PP_HDG
USE MOD_DG_Vars,                 ONLY:U
#endif
USE MOD_PIC_Vars!,      ONLY: 
USE MOD_PICInterpolation_Vars,   ONLY:useVariableExternalField,externalField,DoInterpolation,InterpolationType
USE MOD_PICDepo_Vars,            ONLY:DepositionType,GaussBorder
USE MOD_Eval_xyz,                ONLY:Eval_xyz_elemcheck,Eval_XYZ_Curved,Eval_xyz_Part2
#ifdef PP_POIS
USE MOD_Equation_Vars,           ONLY:E
#endif
#ifdef PP_HDG
#if PP_nVar==1
USE MOD_Equation_Vars,        ONLY:E
#elif PP_nVar==3
USE MOD_Equation_Vars,        ONLY:B
#else
USE MOD_Equation_Vars,        ONLY:B,E
#endif /*PP_nVar==1*/
#endif /*PP_HDG*/
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
USE MOD_Particle_Vars,        ONLY:DoSurfaceFlux
#endif /*HDG-LSERK*/

!----------------------------------------------------------------------------------------------------------------------------------
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: PartID
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)             :: FieldAtParticle(1:6)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                                    
REAL                         :: Pos(3),Field(1:6)
INTEGER                      :: ElemID
! for Nearest GaussPoint
INTEGER                          :: a,b,k,ii,l,m
#if defined PP_POIS || (defined PP_HDG && PP_nVar==4)
REAL                             :: HelperU(1:6,0:PP_N,0:PP_N,0:PP_N)
#endif /*(PP_POIS||PP_HDG)*/
LOGICAL                          :: NotMappedSurfFluxParts
!===================================================================================================================================
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
NotMappedSurfFluxParts=DoSurfaceFlux !Surfaceflux particles inserted before interpolation and tracking. Field at wall is needed!
#else
NotMappedSurfFluxParts=.FALSE.
#endif /*HDG-LSERK*/
FieldAtParticle=0.
IF(useVariableExternalField) THEN ! used Variable external Bz
  FieldAtParticle(:) = 0.
  FieldAtParticle(1) = externalField(1)
  FieldAtParticle(2) = externalField(2)
  FieldAtParticle(3) = externalField(3)
#if (PP_nVar==8)
  FieldAtParticle(4) = externalField(4)
  FieldAtParticle(5) = externalField(5)
#endif
  ! Bz field strength at particle position
  FieldAtparticle(6) = InterpolateVariableExternalField(PartState(PartID,3))
ELSE ! useVariableExternalField
  FieldAtParticle(:) = 0.
  FieldAtParticle(1) = externalField(1)
  FieldAtParticle(2) = externalField(2)
  FieldAtParticle(3) = externalField(3)
#if (PP_nVar==8)
  FieldAtParticle(4) = externalField(4)
  FieldAtParticle(5) = externalField(5)
  FieldAtParticle(6) = externalField(6)
#endif
END IF ! use constant external field

IF (DoInterpolation) THEN                 ! skip if no self fields are calculated
  ElemID=PEM%Element(PartID)
#ifdef MPI
  IF(ElemID.GT.PP_nElems) RETURN
#endif
  SELECT CASE(TRIM(InterpolationType))
  CASE('nearest_blurrycenter')
    ! add fields to fields at particle position
    ! should nearest_blurrycenter not require to
    ! a) evaluate the polynomial at Xi=0??
    ! b) require the mean value of the field
    !m = INT(PP_N/2)+1
#if (PP_nVar==8)
#ifdef PP_POIS
    HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
    HelperU(4:6,:,:,:) = U(4:6,:,:,:,ElemID)
    CALL Eval_xyz_Part2((/0.,0.,0./),6,PP_N,HelperU,field(1:6),ElemID)
#else
    CALL Eval_xyz_Part2((/0.,0.,0./),6,PP_N,U(1:6,:,:,:,ElemID),field(1:6),ElemID)
#endif /*PP_POIS*/
#else
#ifdef PP_POIS
    CALL Eval_xyz_Part2((/0.,0.,0./),3,PP_N,E(1:3,:,:,:,ElemID),field(1:3),ElemID)
#elif defined PP_HDG
#if PP_nVar==1
    CALL Eval_xyz_Part2((/0.,0.,0./),3,PP_N,E(1:3,:,:,:,ElemID),field(1:3),ElemID)
#elif PP_nVar==3
    CALL Eval_xyz_Part2((/0.,0.,0./),3,PP_N,B(1:3,:,:,:,ElemID),field(4:6),ElemID)
#else
    HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
    HelperU(4:6,:,:,:) = B(1:3,:,:,:,ElemID)
    CALL Eval_xyz_Part2((/0.,0.,0./),6,PP_N,HelperU,field(1:6),ElemID)
#endif /*PP_nVar==1*/
#else
    CALL Eval_xyz_Part2((/0.,0.,0./),3,PP_N,U(1:3,:,:,:,ElemID),field(1:3),ElemID)
#endif /*PP_POIS*/
#endif /*(PP_nVar==8)*/
    FieldAtParticle(:) = FieldAtParticle(:) + field(1:6)
  CASE('particle_position_slow')
    Pos = PartState(PartID,1:3)
    !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
    HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
    HelperU(4:6,:,:,:) = U(4:6,:,:,:,ElemID)
    CALL eval_xyz_curved(Pos,6,PP_N,HelperU,field(1:6),ElemID)
#else
    CALL eval_xyz_curved(Pos,6,PP_N,U(1:6,:,:,:,ElemID),field(1:6),ElemID)
#endif /*PP_POIS*/
#else
#ifdef PP_POIS
    CALL eval_xyz_curved(Pos,3,PP_N,E(1:3,:,:,:,ElemID),field(1:3),ElemID)
#elif defined PP_HDG
#if PP_nVar==1
    CALL eval_xyz_curved(Pos,3,PP_N,E(1:3,:,:,:,ElemID),field(1:3),ElemID)
#elif PP_nVar==3
    CALL eval_xyz_curved(Pos,3,PP_N,B(1:3,:,:,:,ElemID),field(4:6),ElemID)
#else
    HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
    HelperU(4:6,:,:,:) = B(1:3,:,:,:,ElemID)
    CALL eval_xyz_curved(Pos,6,PP_N,HelperU,field(1:6),ElemID)
#endif /*PP_nVar*/
#else
    CALL eval_xyz_curved(Pos,3,PP_N,U(1:3,:,:,:,ElemID),field(1:3),ElemID)
#endif /*PP_POIS*/
#endif /*(PP_nVar==8)*/
    FieldAtParticle(:) = FieldAtParticle(:) + field(1:6)
  CASE('particle_position')
    IF(.NOT.NotMappedSurfFluxParts .AND.(DoRefMapping .OR. TRIM(DepositionType).EQ.'nearest_gausspoint'))THEN
      ! particles have already been mapped in deposition, other eval routine used
      IF(.NOT.DoRefMapping)THEN
        CALL Eval_xyz_ElemCheck(PartState(PartID,1:3),PartPosRef(1:3,PartID),ElemID)
      END IF
      !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
      HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
      HelperU(4:6,:,:,:) = U(4:6,:,:,:,ElemID)
      CALL eval_xyz_part2(PartPosRef(1:3,PartID),6,PP_N,HelperU,field(1:6),ElemID)
#else
      CALL eval_xyz_part2(PartPosRef(1:3,PartID),6,PP_N,U(1:6,:,:,:,ElemID),field(1:6),ElemID)
#endif
#else
#ifdef PP_POIS
      CALL eval_xyz_part2(PartPosRef(1:3,PartID),3,PP_N,E(1:3,:,:,:,ElemID),field(1:3),ElemID)
#elif defined PP_HDG
#if PP_nVar==1
      CALL eval_xyz_part2(PartPosRef(1:3,PartID),3,PP_N,E(1:3,:,:,:,ElemID),field(1:3),ElemID)
#elif PP_nVar==3
      CALL eval_xyz_part2(PartPosRef(1:3,PartID),3,PP_N,B(1:3,:,:,:,ElemID),field(4:6),ElemID)
#else
      HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
      HelperU(4:6,:,:,:) = B(1:3,:,:,:,ElemID)
      CALL eval_xyz_part2(PartPosRef(1:3,PartID),6,PP_N,HelperU,field(1:6),ElemID)
#endif
#else
      CALL eval_xyz_part2(PartPosRef(1:3,PartID),3,PP_N,U(1:3,:,:,:,ElemID),field(1:3),ElemID)
#endif
#endif
      FieldAtParticle(:) = FieldAtParticle(:) + field(1:6)
    ELSE IF(NotMappedSurfFluxParts .AND.(DoRefMapping .OR. TRIM(DepositionType).EQ.'nearest_gausspoint'))THEN
      !some particle are mapped, surfaceflux particles (dtFracPush) are not
      IF(PDM%dtFracPush(PartID))THEN ! same as in "particles are not yet mapped"
        Pos = PartState(PartID,1:3)
        !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
        HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
        HelperU(4:6,:,:,:) = U(4:6,:,:,:,ElemID)
        CALL eval_xyz_curved(Pos,6,PP_N,HelperU,field(1:6),ElemID,PartID)
#else
        CALL eval_xyz_curved(Pos,6,PP_N,U(1:6,:,:,:,ElemID),field(1:6),ElemID,PartID)
#endif
#else
#ifdef PP_POIS
        CALL eval_xyz_curved(Pos,3,PP_N,E(1:3,:,:,:,ElemID),field(1:3),ElemID,PartID)
#elif defined PP_HDG
#if PP_nVar==1
        CALL eval_xyz_curved(Pos,3,PP_N,E(1:3,:,:,:,ElemID),field(1:3),ElemID,PartID)
#elif PP_nVar==3
        CALL eval_xyz_curved(Pos,3,PP_N,B(1:3,:,:,:,ElemID),field(4:6),ElemID,PartID)
#else
        HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
        HelperU(4:6,:,:,:) = B(1:3,:,:,:,ElemID)
        CALL eval_xyz_curved(Pos,6,PP_N,HelperU,field(1:6),ElemID,PartID)
#endif
#else
        CALL eval_xyz_curved(Pos,3,PP_N,U(1:3,:,:,:,ElemID),field(1:3),ElemID,PartID)
#endif
#endif
      ELSE !.NOT.PDM%dtFracPush(PartID): same as in "particles have already been mapped in deposition, other eval routine used"
        IF(.NOT.DoRefMapping)THEN
          CALL Eval_xyz_ElemCheck(PartState(PartID,1:3),PartPosRef(1:3,PartID),ElemID)
        END IF
        !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
        HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
        HelperU(4:6,:,:,:) = U(4:6,:,:,:,ElemID)
        CALL eval_xyz_part2(PartPosRef(1:3,PartID),6,PP_N,HelperU,field(1:6),ElemID)
#else
        CALL eval_xyz_part2(PartPosRef(1:3,PartID),6,PP_N,U(1:6,:,:,:,ElemID),field(1:6),ElemID)
#endif
#else
#ifdef PP_POIS
        CALL eval_xyz_part2(PartPosRef(1:3,PartID),3,PP_N,E(1:3,:,:,:,ElemID),field(1:3),ElemID)
#elif defined PP_HDG
#if PP_nVar==1
        CALL eval_xyz_part2(PartPosRef(1:3,PartID),3,PP_N,E(1:3,:,:,:,ElemID),field(1:3),ElemID)
#elif PP_nVar==3
        CALL eval_xyz_part2(PartPosRef(1:3,PartID),3,PP_N,B(1:3,:,:,:,ElemID),field(4:6),ElemID)
#else
        HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
        HelperU(4:6,:,:,:) = B(1:3,:,:,:,ElemID)
        CALL eval_xyz_part2(PartPosRef(1:3,PartID),6,PP_N,HelperU,field(1:6),ElemID)
#endif
#else
        CALL eval_xyz_part2(PartPosRef(1:3,PartID),3,PP_N,U(1:3,:,:,:,ElemID),field(1:3),ElemID)
#endif
#endif
      END IF !PDM%dtFracPush(PartID)
      FieldAtParticle(:) = FieldAtParticle(:) + field(1:6)
    ELSE ! particles are not yet mapped
      Pos = PartState(PartID,1:3)
      !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
      HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
      HelperU(4:6,:,:,:) = U(4:6,:,:,:,ElemID)
      CALL eval_xyz_curved(Pos,6,PP_N,HelperU,field(1:6),ElemID,PartID)
#else
      CALL eval_xyz_curved(Pos,6,PP_N,U(1:6,:,:,:,ElemID),field(1:6),ElemID,PartID)
#endif
#else
#ifdef PP_POIS
      CALL eval_xyz_curved(Pos,3,PP_N,E(1:3,:,:,:,ElemID),field(1:3),ElemID,PartID)
#elif defined PP_HDG
#if PP_nVar==1
      CALL eval_xyz_curved(Pos,3,PP_N,E(1:3,:,:,:,ElemID),field(1:3),ElemID,PartID)
#elif PP_nVar==3
      CALL eval_xyz_curved(Pos,3,PP_N,B(1:3,:,:,:,ElemID),field(4:6),ElemID,PartID)
#else
      HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
      HelperU(4:6,:,:,:) = B(1:3,:,:,:,ElemID)
      CALL eval_xyz_curved(Pos,6,PP_N,HelperU,field(1:6),ElemID,PartID)
#endif
#else
      CALL eval_xyz_curved(Pos,3,PP_N,U(1:3,:,:,:,ElemID),field(1:3),ElemID,PartID)
#endif
#endif
      FieldAtParticle(:) = FieldAtParticle(:) + field(1:6)
    END IF ! DoRefMapping .or. Depositiontype=nearest_gausspoint
  CASE('nearest_gausspoint')
    ! particles have already been mapped in deposition
    IF(MOD(PP_N,2).EQ.0) THEN
      a = PP_N/2
      b = a
    ELSE
      a = (PP_N+1)/2
      b = a-1
    END IF
    IF(.NOT.DoRefMapping .OR. (NotMappedSurfFluxParts .AND. PDM%dtFracPush(PartID)))THEN
      CALL Eval_xyz_ElemCheck(PartState(PartID,1:3),PartPosRef(1:3,PartID),ElemID)
    END IF
    ! compute exact k,l,m
    !! x-direction
    k = a
    DO ii = 0,b-1
      IF(ABS(PartPosRef(1,PartID)).GE.GaussBorder(PP_N-ii))THEN
        k = PP_N-ii
        EXIT
      END IF
    END DO
    k = NINT((PP_N+SIGN(2.0*k-PP_N,PartPosRef(1,PartID)))/2)
    !! y-direction
    l = a
    DO ii = 0,b-1
      IF(ABS(PartPosRef(2,PartID)).GE.GaussBorder(PP_N-ii))THEN
        l = PP_N-ii
        EXIT
      END IF
    END DO
    l = NINT((PP_N+SIGN(2.0*l-PP_N,PartPosRef(2,PartID)))/2)
    !! z-direction
    m = a
    DO ii = 0,b-1
      IF(ABS(PartPosRef(3,PartID)).GE.GaussBorder(PP_N-ii))THEN
        m = PP_N-ii
        EXIT
      END IF
    END DO
    m = NINT((PP_N+SIGN(2.0*m-PP_N,PartPosRef(3,PartID)))/2)
    PartPosGauss(PartID,1) = k
    PartPosGauss(PartID,2) = l
    PartPosGauss(PartID,3) = m
    !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
    field(1:3) = E(1:3,PartPosGauss(PartID,1),PartPosGauss(PartID,2),PartPosGauss(PartID,3), ElemID)
    field(4:6) = U(4:6,PartPosGauss(PartID,1),PartPosGauss(PartID,2),PartPosGauss(PartID,3), ElemID)
    FieldAtParticle(:) = FieldAtParticle(:) + field(1:6)
#else
    field = U(1:6,PartPosGauss(PartID,1),PartPosGauss(PartID,2),PartPosGauss(PartID,3), ElemID)
    FieldAtParticle(:) = FieldAtParticle(:) + field(1:6)
#endif
#else
#ifdef PP_POIS
    field(1:3) = E(1:3,PartPosGauss(PartID,1),PartPosGauss(PartID,2),PartPosGauss(PartID,3), ElemID)
    FieldAtParticle(1:3) = FieldAtParticle(1:3) + field(1:3)
#elif defined PP_HDG
#if PP_nVar==1
    field(1:3) = E(1:3,PartPosGauss(PartID,1),PartPosGauss(PartID,2),PartPosGauss(PartID,3), ElemID)
    FieldAtParticle(1:3) = FieldAtParticle(1:3) + field(1:3)
#elif PP_nVar==3
    field(4:6) = B(1:3,PartPosGauss(PartID,1),PartPosGauss(PartID,2),PartPosGauss(PartID,3), ElemID)
    FieldAtParticle(4:6) = FieldAtParticle(4:6) + field(4:6)
#else
    field(1:3) = E(1:3,PartPosGauss(PartID,1),PartPosGauss(PartID,2),PartPosGauss(PartID,3), ElemID)
    field(4:6) = B(1:3,PartPosGauss(PartID,1),PartPosGauss(PartID,2),PartPosGauss(PartID,3), ElemID)
    FieldAtParticle(:) = FieldAtParticle(:) + field(1:6)
#endif
#else
    field(1:3) = U(1:3,PartPosGauss(PartID,1),PartPosGauss(PartID,2),PartPosGauss(PartID,3), ElemID)
    FieldAtParticle(1:3) = FieldAtParticle(1:3) + field(1:3)
#endif
#endif
  CASE DEFAULT
    CALL abort(&
__STAMP__&
    , 'ERROR: Unknown InterpolationType!')
  END SELECT
END IF
    
RETURN
END SUBROUTINE InterpolateFieldToSingleParticle


SUBROUTINE ReadVariableExternalField()
!===================================================================================================================================
! ATTENTION: The extrenal field needs to be defined on equidistant data-points
! Usage Information
! The file for the variable Bfield contains only the z coordinates and the static Bz-field
! Use the following format F8.5,1x,F8.5
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PICInterpolation_Vars, ONLY:VariableExternalField,DeltaExternalField,nIntPoints,FileNameVariableExternalField
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: ioUnit, ii, err, ncounts
REAL                  :: dummy, diff_comp, diff_check
CHARACTER(LEN=256)        :: FileName
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A,3X,A,65X,A)') ' INITIALIZATION OF VARIABLE EXTERNAL FIELD FOR PARTICLES '
!OPEN(NEWUNIT=ioUnit,FILE=VariableExternalField,STATUS='OLD',FORM='FORMATTED')
OPEN(NEWUNIT=ioUnit,FILE=FileNameVariableExternalField,STATUS='OLD')
err = 0
ncounts = 0
DO WHILE (err.EQ.0)
  READ(ioUnit,*,IOSTAT = err) dummy
  IF (err.EQ.-1) THEN
    EXIT
  END IF 
  ERR = 0
  ncounts = ncounts + 1 
END DO
print*,"ncounts=",ncounts
REWIND(ioUnit)
nIntPoints = ncounts
! allocate needed space
ALLOCATE(VariableExternalField(1:2,1:nIntPoints))
DO ii = 1, ncounts
  read(ioUnit,*) VariableExternalField(1,ii) , VariableExternalField(2,ii)
  IF (ii.GE.2) THEN
    diff_comp  = VariableExternalField(1,2)  - VariableExternalField(1,1)
    diff_check = VariableExternalField(1,ii) - VariableExternalField(1,ii-1)
    IF( (.NOT.ALMOSTEQUALRELATIVE(diff_comp,diff_check,1E-5)) .AND. ((diff_comp.GT.0.0).AND.(diff_check.GT.0.0)) )THEN
      SWRITE(UNIT_stdOut,'(A)') "ReadVariableExternalField: Non-equidistant OR non-increasing points for variable external field." 
      SWRITE(UNIT_stdOut,OUTPUTFORMAT) diff_comp
      SWRITE(UNIT_stdOut,OUTPUTFORMAT) diff_check
      CALL abort(&
__STAMP__&
        ,' Error in dataset!')
    END IF  
  END IF
END DO
CLOSE (ioUnit)

!IF (VariableExternalField(1,1) .NE.0) THEN
  !CALL abort(&
!__STAMP__&
!,  &
      !"ERROR: Points have to start at 0.")
!END IF
IF(ncounts.GT.1) THEN
  DeltaExternalField = VariableExternalField(1,2)  - VariableExternalField(1,1)
  SWRITE(UNIT_stdOut,'(A,1X,E25.14E3)') ' Delta external field: ',DeltaExternalField
  IF(DeltaExternalField.LE.0) THEN
    SWRITE(*,'(A)') ' ERROR: wrong sign in external field delta-x'
  END IF
ELSE 
  CALL abort(&
__STAMP__&
, &
      " ERROR: not enough data points in variable external field file!")
END IF
SWRITE(UNIT_stdOut,'(A,I4.0,A)')' Found ', ncounts,' data points.'
SWRITE(UNIT_stdOut,'(A)')' ...VARIABLE EXTERNAL FIELD INITIALIZATION DONE'
END SUBROUTINE ReadVariableExternalField


PURE FUNCTION InterpolateVariableExternalField(Pos)
!===================================================================================================================================
! interpolates Variable external field to z position
! NO z-values smaller than VariableExternalField(1,1) are allowed!
!===================================================================================================================================
! MODULES
!USE MOD_Globals
USE MOD_PICInterpolation_Vars   ,ONLY:DeltaExternalField,nIntPoints,VariableExternalField
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)          :: Pos ! partilce z position
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                     :: InterpolateVariableExternalField  ! Bz
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                  :: iPos
!===================================================================================================================================
iPos = INT((Pos-VariableExternalField(1,1))/DeltaExternalField) + 1
IF(iPos.GE.nIntPoints)THEN ! particle outside of range (greater -> use constant value)
  InterpolateVariableExternalField = VariableExternalField(2,nIntPoints)
ELSEIF(iPos.LT.1)THEN ! particle outside of range (lower -> use constant value)
  InterpolateVariableExternalField = VariableExternalField(2,1)
ELSE ! Linear Interpolation between iPos and iPos+1 B point
  InterpolateVariableExternalField = (VariableExternalField(2,iPos+1) - VariableExternalField(2,iPos)) & !  dy
                                   / (VariableExternalField(1,iPos+1) - VariableExternalField(1,iPos)) & ! /dx
                             * (Pos - VariableExternalField(1,iPos) ) + VariableExternalField(2,iPos)    ! *(z - z_i) + z_i
END IF
END FUNCTION InterpolateVariableExternalField 

END MODULE MOD_PICInterpolation
