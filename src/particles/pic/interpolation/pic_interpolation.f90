!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (gitlab.com/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

MODULE  MOD_PICInterpolation
!===================================================================================================================================
!
!===================================================================================================================================
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: InterpolateFieldToParticle
PUBLIC :: InitializeParticleInterpolation
PUBLIC :: InterpolateFieldToSingleParticle
PUBLIC :: InterpolateVariableExternalField
#ifdef CODE_ANALYZE
PUBLIC :: InitAnalyticalParticleState
#endif /*CODE_ANALYZE*/
!===================================================================================================================================
INTERFACE InitializeParticleInterpolation
  MODULE PROCEDURE InitializeParticleInterpolation
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

SUBROUTINE InitializeParticleInterpolation
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
REAL                      :: scaleExternalField
#ifdef CODE_ANALYZE
CHARACTER(LEN=20)         :: tempStr
#endif /*CODE_ANALYZE*/
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE INTERPOLATION...'

InterpolationType = GETSTR('PIC-Interpolation-Type','particle_position')

InterpolationElemLoop = GETLOGICAL('PIC-InterpolationElemLoop')
IF (InterpolationElemLoop) THEN !If user-defined F: F for all procs
  IF (PP_nElems.GT.10) THEN !so far arbitrary threshold...
    InterpolationElemLoop=.FALSE. !switch off for procs with high number of Elems
  END IF
END IF
externalField(1:6) = GETREALARRAY('PIC-externalField',6)
scaleexternalField = GETREAL('PIC-scaleexternalField')
externalField      = externalField*ScaleExternalField

DoInterpolation    = GETLOGICAL('PIC-DoInterpolation')
useBGField         = GETLOGICAL('PIC-BG-Field')

! Variable external field
useVariableExternalField = .FALSE.
FileNameVariableExternalField=GETSTR('PIC-curvedexternalField')     ! old variable name (for backward compatibility)
IF (FileNameVariableExternalField.EQ.'none') THEN                          ! if not supplied, check the new variable name
  FileNameVariableExternalField=GETSTR('PIC-variableexternalField') ! new variable name (overwrites the old)
END IF
IF (FileNameVariableExternalField.NE.'none') THEN ! if supplied, read the data file
  useVariableExternalField = .TRUE.
  CALL ReadVariableExternalField()
END IF

!--- Allocate arrays for interpolation of fields to particles
SDEALLOCATE(FieldAtParticle)
ALLOCATE(FieldAtParticle(1:6,1:PDM%maxParticleNumber), STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
  __STAMP__ &
  ,'ERROR in pic_interpolation.f90: Cannot allocate FieldAtParticle array!',ALLOCSTAT)
END IF

SELECT CASE(TRIM(InterpolationType))
CASE('particle_position')
  ! PASS
CASE DEFAULT
  CALL abort(&
  __STAMP__ &
  ,'Unknown InterpolationType ['//TRIM(ADJUSTL(InterpolationType))//'] in pic_interpolation.f90')
END SELECT

#ifdef CODE_ANALYZE
! Initialize analytic solutions for particle time integration (checking the order of convergence for time discretizations)
DoInterpolationAnalytic   = GETLOGICAL('PIC-DoInterpolationAnalytic')
IF(DoInterpolationAnalytic)THEN
  AnalyticInterpolationType = GETINT('PIC-AnalyticInterpolation-Type')
  AnalyticInterpolationPhase = GETREAL('PIC-AnalyticInterpolationPhase')
  SELECT CASE(AnalyticInterpolationType)
  CASE(0) ! 0: const. magnetostatic field: B = B_z = (/ 0 , 0 , 1 T /) = const.
    ! no special parameters required
  CASE(1) ! 1: magnetostatic field: B = B_z = (/ 0 , 0 , B_0 * EXP(x/l) /) = const.
    AnalyticInterpolationSubType = GETINT('PIC-AnalyticInterpolation-SubType')
    AnalyticInterpolationP       = GETREAL('PIC-AnalyticInterpolationP')
  CASE(2) !2: const. electromagnetic field: B = B_z = (/ 0 , 0 , (x^2+y^2)^0.5 /) = const.
          !                                 E = 1e-2/(x^2+y^2)^(3/2) * (/ x , y , 0. /)
    ! no special parameters required
  CASE DEFAULT
    WRITE(TempStr,'(I5)') AnalyticInterpolationType
    CALL abort(&
        __STAMP__ &
        ,'Unknown PIC-AnalyticInterpolation-Type "'//TRIM(ADJUSTL(TempStr))//'" in pic_interpolation.f90')
  END SELECT

  ! Calculate the initial velocity of the particle from an analytic expression: must be implemented for the different 
  ! AnalyticInterpolationType methods
  ! Note that for time-staggered methods, Leapfrog and Boris, the initial velocity in shifted by -dt/2 into the past
  IF(DoInterpolationAnalytic.AND.ANY((/0,1/).EQ.AnalyticInterpolationType))THEN
    DoInitAnalyticalParticleState = .TRUE.
  ELSE
    DoInitAnalyticalParticleState = .FALSE.
  END IF
END IF
#endif /*CODE_ANALYZE*/

SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE INTERPOLATION DONE!'
END SUBROUTINE InitializeParticleInterpolation


SUBROUTINE InterpolateFieldToParticle(DoInnerParts)
!===================================================================================================================================
! Calculates the electromagnetic fields at all the particle's positions
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Vars          ,ONLY: PartPosRef,PDM,PartState,PEM,PartPosGauss,DoFieldIonization
USE MOD_Particle_Tracking_Vars ,ONLY: DoRefMapping
USE MOD_Part_Tools             ,ONLY: isInterpolateParticle
#if !(USE_HDG)
USE MOD_DG_Vars                ,ONLY: U
#endif
USE MOD_PIC_Vars
USE MOD_PICInterpolation_Vars  ,ONLY: useVariableExternalField,FieldAtParticle,externalField,DoInterpolation,InterpolationType
USE MOD_PICInterpolation_Vars  ,ONLY: InterpolationElemLoop
USE MOD_PICDepo_Vars           ,ONLY: DepositionType,GaussBorder
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem,EvaluateFieldAtPhysPos,EvaluateFieldAtRefPos
#ifdef PP_POIS
USE MOD_Equation_Vars          ,ONLY: E
#endif
#if USE_HDG
#if PP_nVar==1
USE MOD_Equation_Vars          ,ONLY: E
#elif PP_nVar==3
USE MOD_Equation_Vars          ,ONLY: B
#else
USE MOD_Equation_Vars          ,ONLY: B,E
#endif /*PP_nVar==1*/
#endif /*USE_HDG*/
#if (PP_TimeDiscMethod>=500) && (PP_TimeDiscMethod<=509)
USE MOD_Particle_Vars          ,ONLY: DoSurfaceFlux
#endif /*(PP_TimeDiscMethod>=500) && (PP_TimeDiscMethod<=509)*/
#if USE_MPI
! only required for shape function??  only required for shape function??
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPIExchange
#endif
#ifdef CODE_ANALYZE
USE MOD_PICInterpolation_Vars  ,ONLY: DoInterpolationAnalytic,AnalyticInterpolationType
#endif /* CODE_ANALYZE */
USE MOD_PICInterpolation_Vars  ,ONLY: CalcBField
USE MOD_Interpolation_Vars     ,ONLY: BGField
USE MOD_SuperB_Vars            ,ONLY: TimeDepCoil, nTimePoints, BGFieldTDep
USE MOD_TimeDisc_Vars          ,ONLY: Time, TEnd
USE MOD_HDF5_Output_Tools      ,ONLY: WriteBGFieldToHDF5
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL                          :: DoInnerParts
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: firstPart,lastPart
REAL                             :: Pos(3), timestep
REAL                             :: field(6)
INTEGER                          :: iPart,iElem,iTime
! for Nearest GaussPoint
INTEGER                          :: a,b,k,ii,l,m
#if defined PP_POIS || (USE_HDG && PP_nVar==4)
REAL                             :: HelperU(1:6,0:PP_N,0:PP_N,0:PP_N)
#endif /*(PP_POIS||USE_HDG)*/
LOGICAL                          :: NotMappedSurfFluxParts
!===================================================================================================================================

! Calculate the time step of the discretization of the Current
IF (CalcBField) THEN
  IF (ANY(TimeDepCoil)) THEN
    timestep = tEnd / (nTimePoints - 1)
    iTime = FLOOR(Time / timestep)

    ! Interpolate the Background field linear between two timesteps
    BGField(:,:,:,:,:) = BGFieldTDep(:,:,:,:,:,iTime) + (BGFieldTDep(:,:,:,:,:,iTime) - BGFieldTDep(:,:,:,:,:,iTime+1)) &
                         / timestep * (Time - iTime * timestep)
    ! CALL WriteBGFieldToHDF5(Time)
  ENDIF
END IF

#if (PP_TimeDiscMethod>=500) && (PP_TimeDiscMethod<=509)
NotMappedSurfFluxParts=DoSurfaceFlux !Surfaceflux particles inserted before interpolation and tracking. Field at wall is needed!
#else
NotMappedSurfFluxParts=.FALSE.
#endif /*(PP_TimeDiscMethod>=500) && (PP_TimeDiscMethod<=509)*/
! null field vector
field=0.

IF(DoInnerParts)THEN
  firstPart=1
  lastPart =PDM%ParticleVecLength
ELSE
#if USE_MPI
  firstPart=PDM%ParticleVecLength-PartMPIExchange%nMPIParticles+1
  lastPart =PDM%ParticleVecLength
#else
  firstPart=2
  LastPart =1
#endif /*USE_MPI*/
END IF
! thats wrong
IF(firstPart.GT.lastPart) RETURN

! IF PP_nElems.GT.10 (so far arbitrary threshold...) use InterpolateFieldToSingleParticle routine
IF (.NOT.InterpolationElemLoop) THEN
  DO iPart = firstPart, LastPart
    IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
    ! Don't interpolate the field at neutral particles (only when considering field ionization)
    IF(DoFieldIonization.OR.isInterpolateParticle(iPart))THEN
      CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle(1:6,iPart))
    END IF
  END DO
  RETURN
END IF

! If PP_nElems.LE.10 (so far arbitrary threshold...) do NOT use InterpolateFieldToSingleParticle routine
FieldAtParticle(:,firstPart:lastPart) = 0. ! initialize
#ifdef CODE_ANALYZE
IF(DoInterpolationAnalytic)THEN ! use analytic/algebraic functions for the field interpolation
  SELECT CASE(AnalyticInterpolationType)
  CASE(0) ! 0: const. magnetostatic field: B = B_z = (/ 0 , 0 , 1 T /) = const.
    DO iPart = firstPart, LastPart
      FieldAtParticle(6,iPart) = 1.0
    END DO
  CASE(1) ! magnetostatic field: B = B_z = B_0 * EXP(x/l)
    ASSOCIATE( B_0 => 1.0     , &
               l   => 1.0  )
               !l   => 0.2e-3  )
      DO iPart = firstPart, LastPart
        FieldAtParticle(6,iPart) = B_0 * EXP(PartState(1,iPart) / l)
      END DO
    END ASSOCIATE
  ! 2: const. electromagnetic field: B = B_z = (/ 0 , 0 , (x^2+y^2)^0.5 /) = const.
  !                                  E = 1e-2/(x^2+y^2)^(3/2) * (/ x , y , 0. /)
  ! Example from Paper by H. Qin: Why is Boris algorithm so good? (2013) 
  ! http://dx.doi.org/10.1063/1.4818428
  CASE(2)
    DO iPart = firstPart, LastPart
      ASSOCIATE( x => PartState(1,iPart) ,&
                 y => PartState(2,iPart) )
        !WRITE (*,*) "x,y,PartState(4,iPart) =", x,y,PartState(4,iPart)
        ! Ex and Ey
        FieldAtParticle(1,iPart) = 1.0e-2 * (x**2+y**2)**(-1.5) * x
        FieldAtParticle(2,iPart) = 1.0e-2 * (x**2+y**2)**(-1.5) * y
        ! Bz
        FieldAtParticle(6,iPart) = SQRT(x**2+y**2)
      END ASSOCIATE
    END DO
  END SELECT
  ! exit the subroutine after field determination
  RETURN
ELSE ! use variable or fixed external field
#endif /*CODE_ANALYZE*/
  IF(useVariableExternalField) THEN ! used variable external Bz, which is given as 1D function Bz(z)
    FieldAtParticle(1,firstPart:lastPart) = externalField(1)
    FieldAtParticle(2,firstPart:lastPart) = externalField(2)
    FieldAtParticle(3,firstPart:lastPart) = externalField(3)
#if (PP_nVar==8)                       
    FieldAtParticle(4,firstPart:lastPart) = externalField(4)
    FieldAtParticle(5,firstPart:lastPart) = externalField(5)
#endif
    ! Bz field strength at particle position
    DO iPart = firstPart, LastPart
      IF(isInterpolateParticle(iPart))THEN
        FieldAtParticle(6,iPart) = InterpolateVariableExternalField(PartState(3,iPart))
      END IF
    END DO
  ELSE ! useVariableExternalField
    FieldAtParticle(1,firstPart:lastPart) = externalField(1)
    FieldAtParticle(2,firstPart:lastPart) = externalField(2)
    FieldAtParticle(3,firstPart:lastPart) = externalField(3)
!#if (PP_nVar==8)
    FieldAtParticle(4,firstPart:lastPart) = externalField(4)
    FieldAtParticle(5,firstPart:lastPart) = externalField(5)
    FieldAtParticle(6,firstPart:lastPart) = externalField(6)
!#endif
  END IF ! use constant external field
#ifdef CODE_ANALYZE
END IF
#endif /*CODE_ANALYZE*/

IF (DoInterpolation) THEN                 ! skip if no self fields are calculated
  SELECT CASE(TRIM(InterpolationType))
  CASE('particle_position')
    IF(.NOT.NotMappedSurfFluxParts .AND.(DoRefMapping .OR. TRIM(DepositionType).EQ.'nearest_gausspoint'))THEN
      ! particles have already been mapped in deposition, other eval routine used
      DO iElem=1,PP_nElems
        DO iPart=firstPart,LastPart
          IF(.NOT.PDM%ParticleInside(iPart))CYCLE
          ! Don't interpolate the field at neutral particles (only when considering field ionization)
          IF(DoFieldIonization.OR.isInterpolateParticle(iPart))THEN
            IF(PEM%GlobalElemID(iPart).EQ.iElem)THEN
              IF(.NOT.DoRefMapping)THEN
                CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(1:3,iPart),iElem)
              END IF
              !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
              HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
              HelperU(4:6,:,:,:) = U(4:6,:,:,:,iElem)
              CALL EvaluateFieldAtRefPos(PartPosRef(1:3,iPart),6,PP_N,HelperU,field(1:6),iElem)
#else
              CALL EvaluateFieldAtRefPos(PartPosRef(1:3,iPart),6,PP_N,U(1:6,:,:,:,iElem),field(1:6),iElem)
#endif
#else
#ifdef PP_POIS
              CALL EvaluateFieldAtRefPos(PartPosRef(1:3,iPart),3,PP_N,E(1:3,:,:,:,iElem),field(1:3),iElem)
#elif USE_HDG
#if PP_nVar==1
              CALL EvaluateFieldAtRefPos(PartPosRef(1:3,iPart),3,PP_N,E(1:3,:,:,:,iElem),field(1:3),iElem)
#elif PP_nVar==3
              CALL EvaluateFieldAtRefPos(PartPosRef(1:3,iPart),3,PP_N,B(1:3,:,:,:,iElem),field(4:6),iElem)
#else
              HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
              HelperU(4:6,:,:,:) = B(1:3,:,:,:,iElem)
              CALL EvaluateFieldAtRefPos(PartPosRef(1:3,iPart),6,PP_N,HelperU,field(1:6),iElem)
#endif
#else
              CALL EvaluateFieldAtRefPos(PartPosRef(1:3,iPart),3,PP_N,U(1:3,:,:,:,iElem),field(1:3),iElem)
#endif
#endif
              FieldAtParticle(:,iPart) = FieldAtParticle(:,iPart) + field(1:6)
            END IF ! Element(iPart).EQ.iElem
          END IF ! DoFieldIonization.OR.isInterpolateParticle(iPart)
        END DO ! iPart
      END DO ! iElem=1,PP_nElems
    ELSE IF(NotMappedSurfFluxParts .AND.(DoRefMapping .OR. TRIM(DepositionType).EQ.'nearest_gausspoint'))THEN
      !some particle are mapped, surfaceflux particles (dtFracPush) are not
      DO iElem=1,PP_nElems
        DO iPart=firstPart,LastPart
          IF(.NOT.PDM%ParticleInside(iPart))CYCLE
          ! Don't interpolate the field at neutral particles (only when considering field ionization)
          IF(DoFieldIonization.OR.isInterpolateParticle(iPart))THEN
            IF(PEM%GlobalElemID(iPart).EQ.iElem)THEN
              IF(PDM%dtFracPush(iPart))THEN ! same as in "particles are not yet mapped"
                Pos = PartState(1:3,iPart)
                !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
                HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
                HelperU(4:6,:,:,:) = U(4:6,:,:,:,iElem)
                CALL EvaluateFieldAtPhysPos(Pos,6,PP_N,HelperU,field(1:6),iElem,iPart)
#else
                CALL EvaluateFieldAtPhysPos(Pos,6,PP_N,U(1:6,:,:,:,iElem),field(1:6),iElem,iPart)
#endif
#else
#ifdef PP_POIS
                CALL EvaluateFieldAtPhysPos(Pos,3,PP_N,E(1:3,:,:,:,iElem),field(1:3),iElem,iPart)
#elif USE_HDG
#if PP_nVar==1
                CALL EvaluateFieldAtPhysPos(Pos,3,PP_N,E(1:3,:,:,:,iElem),field(1:3),iElem,iPart)
#elif PP_nVar==3
                CALL EvaluateFieldAtPhysPos(Pos,3,PP_N,B(1:3,:,:,:,iElem),field(4:6),iElem,iPart)
#else
                HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
                HelperU(4:6,:,:,:) = B(1:3,:,:,:,iElem)
                CALL EvaluateFieldAtPhysPos(Pos,6,PP_N,HelperU,field(1:6),iElem,iPart)
#endif
#else
                CALL EvaluateFieldAtPhysPos(Pos,3,PP_N,U(1:3,:,:,:,iElem),field(1:3),iElem,iPart)
#endif
#endif
              ELSE !.NOT.PDM%dtFracPush(iPart): same as in "particles have already been mapped in deposition, other eval routine used"
                IF(.NOT.DoRefMapping)THEN
                  CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(1:3,iPart),iElem)
                END IF
                !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
                HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
                HelperU(4:6,:,:,:) = U(4:6,:,:,:,iElem)
                CALL EvaluateFieldAtRefPos(PartPosRef(1:3,iPart),6,PP_N,HelperU,field(1:6),iElem)
#else
                CALL EvaluateFieldAtRefPos(PartPosRef(1:3,iPart),6,PP_N,U(1:6,:,:,:,iElem),field(1:6),iElem)
#endif
#else
#ifdef PP_POIS
                CALL EvaluateFieldAtRefPos(PartPosRef(1:3,iPart),3,PP_N,E(1:3,:,:,:,iElem),field(1:3),iElem)
#elif USE_HDG
#if PP_nVar==1
                CALL EvaluateFieldAtRefPos(PartPosRef(1:3,iPart),3,PP_N,E(1:3,:,:,:,iElem),field(1:3),iElem)
#elif PP_nVar==3
                CALL EvaluateFieldAtRefPos(PartPosRef(1:3,iPart),3,PP_N,B(1:3,:,:,:,iElem),field(4:6),iElem)
#else
                HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
                HelperU(4:6,:,:,:) = B(1:3,:,:,:,iElem)
                CALL EvaluateFieldAtRefPos(PartPosRef(1:3,iPart),6,PP_N,HelperU,field(1:6),iElem)
#endif
#else
                CALL EvaluateFieldAtRefPos(PartPosRef(1:3,iPart),3,PP_N,U(1:3,:,:,:,iElem),field(1:3),iElem)
#endif
#endif
              END IF !PDM%dtFracPush(iPart)
              FieldAtParticle(:,iPart) = FieldAtParticle(:,iPart) + field(1:6)
            END IF ! Element(iPart).EQ.iElem
          END IF ! DoFieldIonization.OR.isInterpolateParticle(iPart)
        END DO ! iPart
      END DO ! iElem=1,PP_nElems
    ELSE ! particles are not yet mapped
      DO iElem=1,PP_nElems
        DO iPart=firstPart,LastPart
          IF(.NOT.PDM%ParticleInside(iPart))CYCLE
          ! Don't interpolate the field at neutral particles (only when considering field ionization)
          IF(DoFieldIonization.OR.isInterpolateParticle(iPart))THEN
            IF(PEM%GlobalElemID(iPart).EQ.iElem)THEN
              Pos = PartState(1:3,iPart)
              !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
              HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
              HelperU(4:6,:,:,:) = U(4:6,:,:,:,iElem)
              CALL EvaluateFieldAtPhysPos(Pos,6,PP_N,HelperU,field(1:6),iElem,iPart)
#else
              CALL EvaluateFieldAtPhysPos(Pos,6,PP_N,U(1:6,:,:,:,iElem),field(1:6),iElem,iPart)
#endif
#else
#ifdef PP_POIS
              CALL EvaluateFieldAtPhysPos(Pos,3,PP_N,E(1:3,:,:,:,iElem),field(1:3),iElem,iPart)
#elif USE_HDG
#if PP_nVar==1
              CALL EvaluateFieldAtPhysPos(Pos,3,PP_N,E(1:3,:,:,:,iElem),field(1:3),iElem,iPart)
#elif PP_nVar==3
              CALL EvaluateFieldAtPhysPos(Pos,3,PP_N,B(1:3,:,:,:,iElem),field(4:6),iElem,iPart)
#else
              HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
              HelperU(4:6,:,:,:) = B(1:3,:,:,:,iElem)
              CALL EvaluateFieldAtPhysPos(Pos,6,PP_N,HelperU,field(1:6),iElem,iPart)
#endif
#else
              CALL EvaluateFieldAtPhysPos(Pos,3,PP_N,U(1:3,:,:,:,iElem),field(1:3),iElem,iPart)
#endif
#endif
              FieldAtParticle(:,iPart) = FieldAtParticle(:,iPart) + field(1:6)
            END IF ! Element(iPart).EQ.iElem
          END IF ! DoFieldIonization.OR.isInterpolateParticle(iPart)
        END DO ! iPart
      END DO ! iElem=1,PP_nElems
    END IF ! DoRefMapping .or. Depositiontype=nearest_gausspoint
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
! Calculates the electromagnetic fields at the particle's position (single particle)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Vars,           ONLY:PartPosRef,PDM,PartState,PEM,PartPosGauss
USE MOD_Particle_Tracking_Vars,  ONLY:DoRefMapping
#if !(USE_HDG)
USE MOD_DG_Vars,                 ONLY:U
#endif
USE MOD_PIC_Vars!,      ONLY:
USE MOD_PICInterpolation_Vars,   ONLY:useVariableExternalField,externalField,DoInterpolation,InterpolationType
USE MOD_PICDepo_Vars,            ONLY:DepositionType,GaussBorder
USE MOD_Eval_xyz,                ONLY:GetPositionInRefElem,EvaluateFieldAtPhysPos,EvaluateFieldAtRefPos
#ifdef PP_POIS
USE MOD_Equation_Vars,           ONLY:E
#endif
#if USE_HDG
#if PP_nVar==1
USE MOD_Equation_Vars,        ONLY:E
#elif PP_nVar==3
USE MOD_Equation_Vars,        ONLY:B
#else
USE MOD_Equation_Vars,        ONLY:B,E
#endif /*PP_nVar==1*/
#endif /*USE_HDG*/
#if (PP_TimeDiscMethod>=500) && (PP_TimeDiscMethod<=509)
USE MOD_Particle_Vars,        ONLY:DoSurfaceFlux
#endif /*(PP_TimeDiscMethod>=500) && (PP_TimeDiscMethod<=509)*/
#ifdef CODE_ANALYZE
USE MOD_PICInterpolation_Vars  ,ONLY: DoInterpolationAnalytic!,AnalyticInterpolationType
#endif /* CODE_ANALYZE */
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
#if defined PP_POIS || (USE_HDG && PP_nVar==4)
REAL                             :: HelperU(1:6,0:PP_N,0:PP_N,0:PP_N)
#endif /*(PP_POIS||USE_HDG)*/
LOGICAL                          :: NotMappedSurfFluxParts
!===================================================================================================================================
#if (PP_TimeDiscMethod>=500) && (PP_TimeDiscMethod<=509)
NotMappedSurfFluxParts=DoSurfaceFlux !Surfaceflux particles inserted before interpolation and tracking. Field at wall is needed!
#else
NotMappedSurfFluxParts=.FALSE.
#endif /*(PP_TimeDiscMethod>=500) && (PP_TimeDiscMethod<=509)*/
FieldAtParticle=0.
#ifdef CODE_ANALYZE
IF(DoInterpolationAnalytic)THEN ! use analytic/algebraic functions for the field interpolation
  CALL abort(&
  __STAMP__&
  ,'DoInterpolationAnalytic: Do not call subroutine InterpolateFieldToSingleParticle()')
ELSE ! use variable or fixed external field
#endif /*CODE_ANALYZE*/
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
    FieldAtParticle(6) = InterpolateVariableExternalField(PartState(3,PartID))
  ELSE ! useVariableExternalField
    FieldAtParticle(:) = 0.
    FieldAtParticle(1) = externalField(1)
    FieldAtParticle(2) = externalField(2)
    FieldAtParticle(3) = externalField(3)
!#if (PP_nVar==8)
    FieldAtParticle(4) = externalField(4)
    FieldAtParticle(5) = externalField(5)
    FieldAtParticle(6) = externalField(6)
!#endif

  END IF ! use constant external field
#ifdef CODE_ANALYZE
END IF
#endif /*CODE_ANALYZE*/

IF (DoInterpolation) THEN                 ! skip if no self fields are calculated
  field(1:6)=0.
  ElemID=PEM%GlobalElemID(PartID)
#if USE_MPI
  IF(ElemID.GT.PP_nElems) RETURN
#endif
  SELECT CASE(TRIM(InterpolationType))
  CASE('particle_position')
    IF(.NOT.NotMappedSurfFluxParts .AND.(DoRefMapping .OR. TRIM(DepositionType).EQ.'nearest_gausspoint'))THEN
      ! particles have already been mapped in deposition, other eval routine used
      IF(.NOT.DoRefMapping)THEN
        CALL GetPositionInRefElem(PartState(1:3,PartID),PartPosRef(1:3,PartID),ElemID)
      END IF
      !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
      HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
      HelperU(4:6,:,:,:) = U(4:6,:,:,:,ElemID)
      CALL EvaluateFieldAtRefPos(PartPosRef(1:3,PartID),6,PP_N,HelperU,field(1:6),ElemID)
#else
      CALL EvaluateFieldAtRefPos(PartPosRef(1:3,PartID),6,PP_N,U(1:6,:,:,:,ElemID),field(1:6),ElemID)
#endif
#else
#ifdef PP_POIS
      CALL EvaluateFieldAtRefPos(PartPosRef(1:3,PartID),3,PP_N,E(1:3,:,:,:,ElemID),field(1:3),ElemID)
#elif USE_HDG
#if PP_nVar==1
      CALL EvaluateFieldAtRefPos(PartPosRef(1:3,PartID),3,PP_N,E(1:3,:,:,:,ElemID),field(1:3),ElemID)
#elif PP_nVar==3
      CALL EvaluateFieldAtRefPos(PartPosRef(1:3,PartID),3,PP_N,B(1:3,:,:,:,ElemID),field(4:6),ElemID)
#else
      HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
      HelperU(4:6,:,:,:) = B(1:3,:,:,:,ElemID)
      CALL EvaluateFieldAtRefPos(PartPosRef(1:3,PartID),6,PP_N,HelperU,field(1:6),ElemID)
#endif
#else
      CALL EvaluateFieldAtRefPos(PartPosRef(1:3,PartID),3,PP_N,U(1:3,:,:,:,ElemID),field(1:3),ElemID)
#endif
#endif
      FieldAtParticle(:) = FieldAtParticle(:) + field(1:6)
    ELSE IF(NotMappedSurfFluxParts .AND.(DoRefMapping .OR. TRIM(DepositionType).EQ.'nearest_gausspoint'))THEN
      !some particle are mapped, surfaceflux particles (dtFracPush) are not
      IF(PDM%dtFracPush(PartID))THEN ! same as in "particles are not yet mapped"
        Pos = PartState(1:3,PartID)
        !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
        HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
        HelperU(4:6,:,:,:) = U(4:6,:,:,:,ElemID)
        CALL EvaluateFieldAtPhysPos(Pos,6,PP_N,HelperU,field(1:6),ElemID,PartID)
#else
        CALL EvaluateFieldAtPhysPos(Pos,6,PP_N,U(1:6,:,:,:,ElemID),field(1:6),ElemID,PartID)
#endif
#else
#ifdef PP_POIS
        CALL EvaluateFieldAtPhysPos(Pos,3,PP_N,E(1:3,:,:,:,ElemID),field(1:3),ElemID,PartID)
#elif USE_HDG
#if PP_nVar==1
        CALL EvaluateFieldAtPhysPos(Pos,3,PP_N,E(1:3,:,:,:,ElemID),field(1:3),ElemID,PartID)
#elif PP_nVar==3
        CALL EvaluateFieldAtPhysPos(Pos,3,PP_N,B(1:3,:,:,:,ElemID),field(4:6),ElemID,PartID)
#else
        HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
        HelperU(4:6,:,:,:) = B(1:3,:,:,:,ElemID)
        CALL EvaluateFieldAtPhysPos(Pos,6,PP_N,HelperU,field(1:6),ElemID,PartID)
#endif
#else
        CALL EvaluateFieldAtPhysPos(Pos,3,PP_N,U(1:3,:,:,:,ElemID),field(1:3),ElemID,PartID)
#endif
#endif
      ELSE !.NOT.PDM%dtFracPush(PartID): same as in "particles have already been mapped in deposition, other eval routine used"
        IF(.NOT.DoRefMapping)THEN
          CALL GetPositionInRefElem(PartState(1:3,PartID),PartPosRef(1:3,PartID),ElemID)
        END IF
        !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
        HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
        HelperU(4:6,:,:,:) = U(4:6,:,:,:,ElemID)
        CALL EvaluateFieldAtRefPos(PartPosRef(1:3,PartID),6,PP_N,HelperU,field(1:6),ElemID)
#else
        CALL EvaluateFieldAtRefPos(PartPosRef(1:3,PartID),6,PP_N,U(1:6,:,:,:,ElemID),field(1:6),ElemID)
#endif
#else
#ifdef PP_POIS
        CALL EvaluateFieldAtRefPos(PartPosRef(1:3,PartID),3,PP_N,E(1:3,:,:,:,ElemID),field(1:3),ElemID)
#elif USE_HDG
#if PP_nVar==1
        CALL EvaluateFieldAtRefPos(PartPosRef(1:3,PartID),3,PP_N,E(1:3,:,:,:,ElemID),field(1:3),ElemID)
#elif PP_nVar==3
        CALL EvaluateFieldAtRefPos(PartPosRef(1:3,PartID),3,PP_N,B(1:3,:,:,:,ElemID),field(4:6),ElemID)
#else
        HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
        HelperU(4:6,:,:,:) = B(1:3,:,:,:,ElemID)
        CALL EvaluateFieldAtRefPos(PartPosRef(1:3,PartID),6,PP_N,HelperU,field(1:6),ElemID)
#endif
#else
        CALL EvaluateFieldAtRefPos(PartPosRef(1:3,PartID),3,PP_N,U(1:3,:,:,:,ElemID),field(1:3),ElemID)
#endif
#endif
      END IF !PDM%dtFracPush(PartID)
      FieldAtParticle(:) = FieldAtParticle(:) + field(1:6)
    ELSE ! particles are not yet mapped
      Pos = PartState(1:3,PartID)
      !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
      HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
      HelperU(4:6,:,:,:) = U(4:6,:,:,:,ElemID)
      CALL EvaluateFieldAtPhysPos(Pos,6,PP_N,HelperU,field(1:6),ElemID,PartID)
#else
      CALL EvaluateFieldAtPhysPos(Pos,6,PP_N,U(1:6,:,:,:,ElemID),field(1:6),ElemID,PartID)
#endif
#else
#ifdef PP_POIS
      CALL EvaluateFieldAtPhysPos(Pos,3,PP_N,E(1:3,:,:,:,ElemID),field(1:3),ElemID,PartID)
#elif USE_HDG
#if PP_nVar==1
      CALL EvaluateFieldAtPhysPos(Pos,3,PP_N,E(1:3,:,:,:,ElemID),field(1:3),ElemID,PartID)
#elif PP_nVar==3
      CALL EvaluateFieldAtPhysPos(Pos,3,PP_N,B(1:3,:,:,:,ElemID),field(4:6),ElemID,PartID)
#else
      HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
      HelperU(4:6,:,:,:) = B(1:3,:,:,:,ElemID)
      CALL EvaluateFieldAtPhysPos(Pos,6,PP_N,HelperU,field(1:6),ElemID,PartID)
#endif
#else
      CALL EvaluateFieldAtPhysPos(Pos,3,PP_N,U(1:3,:,:,:,ElemID),field(1:3),ElemID,PartID)
#endif
#endif
      FieldAtParticle(:) = FieldAtParticle(:) + field(1:6)
    END IF ! DoRefMapping .or. Depositiontype=nearest_gausspoint
  CASE DEFAULT
    CALL abort(&
__STAMP__&
    , 'ERROR: Unknown InterpolationType!')
  END SELECT
END IF

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
      SWRITE(UNIT_stdOut,WRITEFORMAT) diff_comp
      SWRITE(UNIT_stdOut,WRITEFORMAT) diff_check
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
  SWRITE(UNIT_stdOut,'(A,1X,ES25.14E3)') ' Delta external field: ',DeltaExternalField
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
!> Interpolates the variable external field to the z-position
!> NO z-values smaller than VariableExternalField(1,1) are allowed!
!===================================================================================================================================
! MODULES
!USE MOD_Globals
USE MOD_PICInterpolation_Vars   ,ONLY:DeltaExternalField,nIntPoints,VariableExternalField
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)          :: Pos                               !< particle z-position
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                     :: InterpolateVariableExternalField  !< Bz (magnetic field in z-direction)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: iPos                              !< index in array (equidistant subdivision assumed)
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


#ifdef CODE_ANALYZE
SUBROUTINE InitAnalyticalParticleState()
!----------------------------------------------------------------------------------------------------------------------------------!
! Calculates the initial particle position and velocity depending on an analytical expression
! The velocity is time-shifted for staggered-in-time methods (Leapfrog and Boris-Leapfrog)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_PICInterpolation_Vars  ,ONLY: DoInitAnalyticalParticleState
USE MOD_Particle_Analyze       ,ONLY: CalcAnalyticalParticleState
USE MOD_Particle_Vars          ,ONLY: PartState, PDM
USE MOD_TimeDisc_Vars          ,ONLY: dt
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: PartStateAnalytic(6)
INTEGER :: iPart
!===================================================================================================================================
! Return here, if no analytical function can be used
IF(.NOT.DoInitAnalyticalParticleState) RETURN

! Calculate the initial velocity of the particle from an analytic expression
DO iPart=1,PDM%ParticleVecLength
  !-- set analytic position at x(n) from analytic particle solution
  CALL CalcAnalyticalParticleState(0.0, PartStateAnalytic)
  PartState(1:6,iPart) = PartStateAnalytic(1:6)

  !-- Only for time-staggered methods (Leapfrog and Boris-Leapfrog):
  ! Set analytic velocity at v(n-0.5) from analytic particle solution
#if (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)
  CALL CalcAnalyticalParticleState(-dt/2., PartStateAnalytic)
  PartState(4:6,iPart) = PartStateAnalytic(4:6)
#endif /*(PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)*/

  ! Set new part to false to prevent calculation of velocity in timedisc
  PDM%IsNewPart(iPart) = .FALSE.
END DO
END SUBROUTINE InitAnalyticalParticleState
#endif /*CODE_ANALYZE*/


END MODULE MOD_PICInterpolation
