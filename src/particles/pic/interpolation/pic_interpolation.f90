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
#ifdef CODE_ANALYZE
PUBLIC :: InitAnalyticalParticleState
#endif /*CODE_ANALYZE*/
PUBLIC :: DefineParametersPICInterpolation
PUBLIC :: FinalizePICInterpolation
!===================================================================================================================================
INTERFACE InterpolateFieldToSingleParticle
  MODULE PROCEDURE InterpolateFieldToSingleParticle
END INTERFACE
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for PIC Interpolation
!==================================================================================================================================
SUBROUTINE DefineParametersPICInterpolation()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools    ,ONLY: prms
USE MOD_PICDepo_Method ,ONLY: DefineParametersDepositionMethod
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("PIC Interpolation")

#ifdef CODE_ANALYZE
! -- external field 1
CALL prms%CreateLogicalOption('PIC-DoInterpolationAnalytic'      , "Method 1 of 5: Use an analytic/algebraic function for PIC interpolation "//&
                                                                   "(ifdef CODE_ANALYZE)",&
                                                                   '.FALSE.')

CALL prms%CreateIntOption(    'PIC-AnalyticInterpolation-Type'   , "Type of AnalyticInterpolation-Method for calculating the "//&
                                                                   "EM field's value for the particle (ifdef CODE_ANALYZE)",'0')

CALL prms%CreateIntOption(    'PIC-AnalyticInterpolation-SubType', "SubType of AnalyticInterpolation-Method for calculating the "//&
                                                                   "EM field's value for the particle (ifdef CODE_ANALYZE)",'0')

CALL prms%CreateRealOption(   'PIC-AnalyticInterpolationP'       , "parameter 'p' for AnalyticInterpolationType = 1", '1.')
CALL prms%CreateRealOption(   'PIC-AnalyticInterpolationPhase'   , "Phase shift angle phi that is used for cos(w*t + phi)", '0.')
#endif /*CODE_ANALYZE*/

CALL prms%CreateLogicalOption(  'PIC-DoInterpolation'         , "Interpolate electric/magnetic fields at charged particle position"//&
 " for calculating the acting Lorentz forces, which are used in the RHS of the evolution equation (particle time integration). "//&
 "Required flag also for using OPTIONAL external fields, for which 5 methods are available:\n"//&
 "  Method 1: PIC-DoInterpolationAnalytic (convergence tests, CODE_ANALYZE=ON)\n"//&
 "  Method 2: PIC-externalField (const. E/B field)\n"//&
 "  Method 3: PIC-variableExternalField (CSV file for Bz(z) that is interpolated)\n"//&
 "  Method 4: PIC-AlgebraicExternalField (E/B pre-defined algebraic expression)\n"//&
 "  Method 5: PIC-BG-Field (read 3D magnetic field from .h5)\n"&
 , '.TRUE.')
CALL prms%CreateStringOption(   'PIC-Interpolation-Type'      , "Type of Interpolation-Method to calculate the electro(-magnetic)"//&
                                                                " field values for the particle push. Poisson solver (E field) "//&
                                                                "and Maxwell solver (E+B field)", 'particle_position')
CALL prms%CreateLogicalOption(  'PIC-InterpolationElemLoop'   , 'Interpolate with outer iElem-loop (increased comp. performance '//&
                                                                ' when using many elements per processor)', '.TRUE.')
! -- external field 2
CALL prms%CreateRealArrayOption('PIC-externalField'           , 'Method 2 of 5: Vector for applying a const./externally acting E and B field'//&
                                                                ' acting on charged particles. This external field is added to the'//&
                                                                'maxwell/poisson-solver-field', '0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0')
CALL prms%CreateRealOption(     'PIC-scaleexternalField'      , 'Scaling factor for PIC-externalField', '1.0')

! -- external field 3
CALL prms%CreateStringOption(   'PIC-variableExternalField'   , 'Method 3 of 5: CSV file containing the external magnetic field Bz in z-direction '//&
                                                                'for interpolating the variable field at each particle z-position.', 'none')

! -- external field 4
CALL prms%CreateIntOption(      'PIC-AlgebraicExternalField'   , &
     'Method 4 of 5: External E and B field from algebraic expression that is interpolated to the particle position\n'//&
     '[1]: Axial Bz(x) field from T. Charoy "2D axial-azimuthal particle-in-cell benchmark for low-temperature partially magnetized plasmas" (2019)\n'//&
     '[2]: Radial Br(x) field in axial x-direction from H. Liu "Particle-in-cell simulation of a Hall thruster" (2010)\n'//&
     '[3]: same as [2] but 3D case, where Br(z) = (/Bx(z), By(z)/) in axial z-direction'&
     , '0')

CALL prms%CreateIntOption(      'PIC-AlgebraicExternalFieldDelta', 'Delta factor for H. Liu "Particle-in-cell simulation of a Hall thruster" (2010)', '2')
CALL prms%CreateRealArrayOption('PIC-NormVecOfWall'  , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Normal vector for pushTimeStep', '1. , 0. , 0.')
CALL prms%CreateRealArrayOption('PIC-BGMdeltas'      , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Dimensions of PIC background mesh', '0. , 0. , 0.')
CALL prms%CreateRealArrayOption('PIC-FactorBGM'      , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Denominator of PIC-BGMdeltas', '1. , 1. , 1.')
CALL prms%CreateLogicalOption(  'PIC-OutputSource'   , 'Flag for activating the output of particle charge and current density'//&
                                                       'source terms to hdf5', '.FALSE.')
END SUBROUTINE DefineParametersPICInterpolation


SUBROUTINE InitializeParticleInterpolation
!===================================================================================================================================
! Initialize the interpolation variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_ReadInTools
USE MOD_Particle_Vars         ,ONLY: PDM
USE MOD_PICInterpolation_Vars
USE MOD_ReadInTools           ,ONLY: PrintOption
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

IF(.NOT.DoInterpolation) THEN
  ! Fill interpolation type with empty string
  InterpolationType='NONE'
  RETURN
END IF

InterpolationType = GETSTR('PIC-Interpolation-Type','particle_position')

InterpolationElemLoop = GETLOGICAL('PIC-InterpolationElemLoop')
IF (InterpolationElemLoop) THEN !If user-defined F: F for all procs
  IF (PP_nElems.GT.10) THEN !so far arbitrary threshold...
    InterpolationElemLoop=.FALSE. !switch off for procs with high number of Elems
    CALL PrintOption('PP_nElems.GT.10: Changeing PIC-InterpolationElemLoop','OUTPUT',LogOpt=InterpolationElemLoop)
  END IF
END IF
externalField(1:6) = GETREALARRAY('PIC-externalField',6)
scaleexternalField = GETREAL('PIC-scaleexternalField')
externalField      = externalField*ScaleExternalField

useBGField         = GETLOGICAL('PIC-BG-Field')

! Variable external field
useVariableExternalField      = .FALSE.
FileNameVariableExternalField = GETSTR('PIC-variableExternalField')
IF (FileNameVariableExternalField.NE.'none') THEN ! if supplied, read the data file
  useVariableExternalField = .TRUE.
  CALL ReadVariableExternalField()
END IF

! Algebraic external field
useAlgebraicExternalField    = .FALSE.
AlgebraicExternalField = GETINT('PIC-AlgebraicExternalField')
IF(AlgebraicExternalField.GT.0) useAlgebraicExternalField=.TRUE.
! Sanity Check: Add all integer values that are possible to the vector for checking
IF(.NOT.ANY(AlgebraicExternalField.EQ.(/0,1,2,3/))) CALL abort(&
  __STAMP__&
  ,'Value for PIC-AlgebraicExternalField not defined',IntInfoOpt=AlgebraicExternalField)

AlgebraicExternalFieldDelta = GETINT('PIC-AlgebraicExternalFieldDelta')
IF(AlgebraicExternalFieldDelta.LT.0) CALL abort(__STAMP__,'AlgebraicExternalFieldDelta cannot be negative.')

!--- Allocate arrays for interpolation of fields to particles
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


SUBROUTINE InterpolateFieldToParticle()
!===================================================================================================================================
! Calculates the electromagnetic fields at all the particle's positions
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Vars         ,ONLY: PDM,PEM,DoFieldIonization,PartState
USE MOD_Part_Tools            ,ONLY: isInterpolateParticle
USE MOD_PIC_Vars
USE MOD_PICInterpolation_Vars ,ONLY: FieldAtParticle,DoInterpolation,InterpolationType
USE MOD_PICInterpolation_Vars ,ONLY: InterpolationElemLoop
USE MOD_PICInterpolation_Vars ,ONLY: CalcBField
USE MOD_HDF5_Output_Fields    ,ONLY: WriteBGFieldToHDF5
#if USE_HDG
USE MOD_AnalyzeField          ,ONLY: CalculateAverageElectricPotential
USE MOD_Analyze_Vars          ,ONLY: CalcAverageElectricPotential
#endif /*USE_HDG*/
USE MOD_PICInterpolation_tools,ONLY:GetExternalFieldAtParticle,GetInterpolatedFieldPartPos
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: iPart,iElem
!===================================================================================================================================
!0. Return if interpolation is not required
IF(.NOT.DoInterpolation) RETURN

!1.1 Calculate the time step of the discretization of the Current
IF (CalcBField) CALL GetTimeDependentBGField()

#if USE_HDG
!1.2 Calculate external E-field
IF(CalcAverageElectricPotential) CALL CalculateAverageElectricPotential()
#endif /*USE_HDG*/

!2. Select element-particle loop (InterpolationElemLoop) or particle-element (.NOT.InterpolationElemLoop)
IF (InterpolationElemLoop) THEN ! element-particle loop
  ! If PP_nElems.LE.10 (so far arbitrary threshold...) do NOT use InterpolateFieldToSingleParticle routine
  !2.1 InterpolationElemLoop (loop elements and then particles)
  SELECT CASE(TRIM(InterpolationType))
  CASE('particle_position')
    ! particles have already been mapped
    DO iElem=1,PP_nElems
      DO iPart=1,PDM%ParticleVecLength
        IF(.NOT.PDM%ParticleInside(iPart)) CYCLE ! Skip particles outside
        IF(.NOT.(DoFieldIonization.OR.isInterpolateParticle(iPart))) CYCLE ! Skip neutral particles (if field ionization if off)
        IF(PEM%LocalElemID(iPart).NE.iElem) CYCLE ! Skip particles that are not inside the current element
        ! Add the interpolated electro-(magnetic) field
        FieldAtParticle(1:6,iPart) = GetExternalFieldAtParticle(PartState(1:3,iPart))
        FieldAtParticle(:,iPart) = FieldAtParticle(:,iPart) + GetInterpolatedFieldPartPos(PEM%GlobalElemID(iPart),iPart)
      END DO ! iPart
    END DO ! iElem=1,PP_nElems
  CASE DEFAULT
    CALL abort(__STAMP__, 'ERROR: Unknown InterpolationType!')
  END SELECT
ELSE ! .NOT.InterpolationElemLoop -> particle-element loop
  ! 2.2 particle-element loop: Loop particles and select corresponding element
  ! IF PP_nElems.GT.10 (so far arbitrary threshold...) use InterpolateFieldToSingleParticle routine
  DO iPart = 1, PDM%ParticleVecLength
    IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
    ! Don't interpolate the field at neutral particles (only when considering field ionization)
    IF(DoFieldIonization.OR.isInterpolateParticle(iPart))THEN
      CALL InterpolateFieldToSingleParticle(iPart,FieldAtParticle(1:6,iPart))
    END IF
  END DO
END IF ! InterpolationElemLoop

RETURN
END SUBROUTINE InterpolateFieldToParticle


SUBROUTINE InterpolateFieldToSingleParticle(PartID,FieldAtParticle)
!===================================================================================================================================
! Calculates the electromagnetic fields at the particle's position (single particle)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_PIC_Vars
USE MOD_Particle_Vars         ,ONLY: PEM,PartState
USE MOD_PICInterpolation_Vars ,ONLY: DoInterpolation,InterpolationType
USE MOD_PICInterpolation_tools,ONLY:GetExternalFieldAtParticle,GetInterpolatedFieldPartPos
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
!===================================================================================================================================
!0. Return if interpolation is not required
IF(.NOT.DoInterpolation) RETURN

!1. Apply any external fields
FieldAtParticle(1:6) = GetExternalFieldAtParticle(PartState(1:3,PartID))

!2. Calculate fields at particle
#if USE_MPI
IF(PEM%LocalElemID(PartID).GT.PP_nElems)THEN! RETURN
  CALL abort(&
  __STAMP__&
  ,'ERROR: This check used to "RETURN" here but is now set to "ABORT". PEM%LocalElemID(PartID).GT.PP_nElems should not happen here.')
END IF
#endif
SELECT CASE(TRIM(InterpolationType))
CASE('particle_position')
  ! Add the interpolated electro-(magnetic) field
  FieldAtParticle(:) = FieldAtParticle(:) + GetInterpolatedFieldPartPos(PEM%GlobalElemID(PartID),PartID)
CASE DEFAULT
  CALL abort(&
  __STAMP__&
  , 'ERROR: Unknown InterpolationType!')
END SELECT

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


SUBROUTINE GetTimeDependentBGField()
!===================================================================================================================================
! Calculates BGField for at t=Time by interpolation between two already calculated time steps
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Interpolation_Vars ,ONLY: BGField
USE MOD_SuperB_Vars        ,ONLY: TimeDepCoil, nTimePoints, BGFieldTDep
USE MOD_TimeDisc_Vars      ,ONLY: Time, TEnd
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iTime
REAL    :: timestep
!===================================================================================================================================
IF (ANY(TimeDepCoil)) THEN
  timestep = tEnd / (nTimePoints - 1)
  iTime = FLOOR(Time / timestep)

  ! Interpolate the Background field linear between two timesteps
  BGField(:,:,:,:,:) = BGFieldTDep(:,:,:,:,:,iTime) + (BGFieldTDep(:,:,:,:,:,iTime) - BGFieldTDep(:,:,:,:,:,iTime+1)) &
                       / timestep * (Time - iTime * timestep)
  ! CALL WriteBGFieldToHDF5(Time)
ENDIF
END SUBROUTINE GetTimeDependentBGField


#ifdef CODE_ANALYZE
SUBROUTINE InitAnalyticalParticleState()
!----------------------------------------------------------------------------------------------------------------------------------!
! Calculates the initial particle position and velocity depending on an analytical expression
! The velocity is time-shifted for staggered-in-time methods (Leapfrog and Boris-Leapfrog)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_PICInterpolation_Vars  ,ONLY: DoInitAnalyticalParticleState
USE MOD_Particle_Analyze_Code  ,ONLY: CalcAnalyticalParticleState
USE MOD_Particle_Vars          ,ONLY: PartState, PDM
#if (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)
USE MOD_TimeDisc_Vars          ,ONLY: dt
#endif /*(PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)*/
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


SUBROUTINE FinalizePICInterpolation()
!----------------------------------------------------------------------------------------------------------------------------------!
! finalize pic interpolation
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_PICInterpolation_Vars ,ONLY: FieldAtParticle
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(FieldAtParticle)
END SUBROUTINE FinalizePICInterpolation


END MODULE MOD_PICInterpolation
