#include "boltzplatz.h"

MODULE  MOD_PICInterpolation
!===================================================================================================================================
! 
!===================================================================================================================================
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: InterpolateFieldToParticle, InitializeInterpolation !, Calc_inv, &
PUBLIC :: InterpolateCurvedExternalField                                  
!===================================================================================================================================
INTERFACE InitializeInterpolation
  MODULE PROCEDURE InitializeInterpolation
END INTERFACE

INTERFACE InterpolateFieldToParticle
  MODULE PROCEDURE InterpolateFieldToParticle
END INTERFACE

INTERFACE read_curved_external_Field
  MODULE PROCEDURE read_curved_external_Field
END INTERFACE

INTERFACE InterpolateCurvedExternalField
  MODULE PROCEDURE InterpolateCurvedExternalField
END INTERFACE
!===================================================================================================================================

CONTAINS

SUBROUTINE InitializeInterpolation
!===================================================================================================================================
! Initialize the interpolation variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
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
!===================================================================================================================================
InterpolationType = GETSTR('PIC-Interpolation-Type','particle_position')
externalField(1:6)= GETREALARRAY('PIC-externalField',6,'0.,0.,0.,0.,0.,0.')
!SWRITE(*,*) " External fied", externalfield
DoInterpolation   = GETLOGICAL('PIC-DoInterpolation','.TRUE.')
useBGField        = GETLOGICAL('PIC-BG-Field','.FALSE.')

! curved external field 
usecurvedExternalField = .FALSE.
FileNameCurvedExternalField = GETSTR('PIC-curvedexternalField','none')
IF (FileNameCurvedExternalField/='none') THEN
  usecurvedExternalField = .TRUE.
  CALL read_curved_external_Field()
END IF

!--- Allocate arrays for interpolation of fields to particles
ALLOCATE(FieldAtParticle(1:PDM%maxParticleNumber,1:6), STAT=ALLOCSTAT) 
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(__STAMP__, &
      'ERROR in pic_interpolation.f90: Cannot allocate FieldAtParticle array!',ALLOCSTAT)
END IF

SELECT CASE(TRIM(InterpolationType))
CASE('nearest_blurycenter')
   InterpolationType='nearest_blurrycenter'
CASE('nearest_blurrycenter')
CASE('particle_position_slow')
CASE('particle_position')
CASE('nearest_gausspoint')
CASE DEFAULT
  CALL abort(__STAMP__, &
      'Unknown InterpolationType in pic_init.f90')
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
USE MOD_DG_Vars,                 ONLY:U
USE MOD_PIC_Vars!,      ONLY: 
USE MOD_PICInterpolation_Vars,   ONLY:usecurvedExternalField,FieldAtParticle,externalField,DoInterpolation,InterpolationType
USE MOD_PICDepo_Vars,            ONLY:DepositionType,GaussBorder
USE MOD_Eval_xyz,                ONLY:Eval_xyz_elemcheck,Eval_XYZ_Curved,Eval_xyz_Part2
USE MOD_Particle_Mesh_Vars,ONLY:epsOneCell
#ifdef PP_POIS
USE MOD_Equation_Vars,           ONLY:E
#endif
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
#ifdef PP_POIS
REAL                             :: HelperU(1:6,0:PP_N,0:PP_N,0:PP_N)
#endif
!===================================================================================================================================

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


IF(usecurvedExternalField) THEN ! used curved external Bz
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
    FieldAtparticle(iPart,6) = InterpolateCurvedExternalField(PartState(iPart,3))
  END DO
ELSE ! usecurvedExternalField
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
#ifdef PP_POIS
      HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
      HelperU(4:6,:,:,:) = U(4:6,:,:,:,iElem)
      CALL Eval_xyz_Part2((/0.,0.,0./),6,PP_N,HelperU,field,iElem)
#else
      CALL Eval_xyz_Part2((/0.,0.,0./),6,PP_N,U(1:6,:,:,:,iElem),field,iElem)
#endif /*PP_POIS*/
#else 
#ifdef PP_POIS
      CALL Eval_xyz_Part2((/0.,0.,0./),3,PP_N,E(1:3,:,:,:,iElem),field,iElem)
#else
      CALL Eval_xyz_Part2((/0.,0.,0./),3,PP_N,U(1:3,:,:,:,iElem),field,iElem)
#endif /*PP_POIS*/
#endif /*PP_nVar*/
      DO iPart=firstPart,LastPart
        IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
        IF(PEM%Element(iPart).EQ.iElem)THEN
          FieldAtParticle(iPart,:) = FieldAtParticle(iPart,:) + field
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
          CALL eval_xyz_curved(Pos,6,PP_N,HelperU,field,iElem)
#else
          CALL eval_xyz_curved(Pos,6,PP_N,U(1:6,:,:,:,iElem),field,iElem,iPart)
#endif
#else
#ifdef PP_POIS
          CALL eval_xyz_curved(Pos,3,PP_N,E(1:3,:,:,:,iElem),field,iElem)
#else
          CALL eval_xyz_curved(Pos,3,PP_N,U(1:3,:,:,:,iElem),field,iElem)
#endif         
#endif
          FieldAtParticle(iPart,:) = FieldAtParticle(iPart,:) + field
        END IF ! Element(iPart).EQ.iElem
      END DO ! iPart
    END DO ! iElem=1,PP_nElems
  CASE('particle_position')
    IF(DoRefMapping .OR. TRIM(DepositionType).EQ.'nearest_gausspoint')THEN
      ! particles have already been mapped in deposition, other eval routine used
      DO iElem=1,PP_nElems
        DO iPart=firstPart,LastPart
          IF(.NOT.PDM%ParticleInside(iPart))CYCLE
          IF(PEM%Element(iPart).EQ.iElem)THEN
            IF(.NOT.DoRefMapping)THEN
              CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),iElem,iPart)
            END IF
            !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
            HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
            HelperU(4:6,:,:,:) = U(4:6,:,:,:,iElem)
            CALL eval_xyz_part2(PartPosRef(1:3,iPart),6,PP_N,HelperU,field,ielem)
#else
            CALL eval_xyz_part2(PartPosRef(1:3,iPart),6,PP_N,U(1:6,:,:,:,iElem),field,ielem)
#endif
#else
#ifdef PP_POIS
            CALL eval_xyz_part2(PartPosRef(1:3,iPart),3,PP_N,E(1:3,:,:,:,iElem),field,iElem)     
#else
            CALL eval_xyz_part2(PartPosRef(1:3,iPart),3,PP_N,U(1:3,:,:,:,iElem),field,iElem)
#endif
#endif
            FieldAtParticle(iPart,:) = FieldAtParticle(iPart,:) + field
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
            CALL eval_xyz_curved(Pos,6,PP_N,HelperU,field,iElem)
#else
            CALL eval_xyz_curved(Pos,6,PP_N,U(1:6,:,:,:,iElem),field,iElem,iPart)
#endif
#else
#ifdef PP_POIS
            CALL eval_xyz_curved(Pos,3,PP_N,E(1:3,:,:,:,iElem),field,iElem)
#else
            CALL eval_xyz_curved(Pos,3,PP_N,U(1:3,:,:,:,iElem),field,iElem)
#endif         
#endif
            FieldAtParticle(iPart,:) = FieldAtParticle(iPart,:) + field
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
          IF(.NOT.DoRefMapping)THEN
            CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),iElem,iPart)
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
         FieldAtParticle(iPart,:) = FieldAtParticle(iPart,:) + field
#else
         field = U(1:6,PartPosGauss(iPart,1),PartPosGauss(iPart,2),PartPosGauss(iPart,3), iElem)
         FieldAtParticle(iPart,:) = FieldAtParticle(iPart,:) + field
#endif
#else
#ifdef PP_POIS
         field(1:3) = E(1:3,PartPosGauss(iPart,1),PartPosGauss(iPart,2),PartPosGauss(iPart,3), iElem)
         FieldAtParticle(iPart,1:3) = FieldAtParticle(iPart,1:3) + field(1:3)
#else
         field(1:3) = U(1:3,PartPosGauss(iPart,1),PartPosGauss(iPart,2),PartPosGauss(iPart,3), iElem)
         FieldAtParticle(iPart,1:3) = FieldAtParticle(iPart,1:3) + field(1:3)
#endif
#endif
        END IF ! Element(iPart).EQ.iElem
      END DO ! iPart
    END DO ! iElem=1,PP_nElems
  CASE DEFAULT
    CALL abort(__STAMP__, &
        'ERROR: Unknown InterpolationType!')
  END SELECT
END IF
    
RETURN
END SUBROUTINE InterpolateFieldToParticle


SUBROUTINE read_curved_external_Field()
!===================================================================================================================================
! ATTENTION: The extrenal field needs to be defined on equidistant data-points
! Usage Information
! The file for the curved Bfield contains only the z coordinates and the static Bz-field
! Use the following format F8.5,1x,F8.5
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PICInterpolation_Vars    ,ONLY:FileNameCurvedExternalField,CurvedExternalField &
                                      ,DeltaExternalField,nIntPoints
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: unit_index_CEF, ii, err, ncounts
  REAL                  :: dummy, diff_comp, diff_check
!===================================================================================================================================

  SWRITE(UNIT_stdOut,'(A,3X,A,65X,A)') ' INITIALIZATION OF CURVED EXTERNAL FIELD FOR PARTICLES '
  unit_index_CEF=147
  OPEN(UNIT=unit_index_CEF,FILE=FileNameCurvedExternalField,STATUS='OLD',FORM='FORMATTED')
  err = 0
  ncounts = 0
  DO WHILE ( err == 0 )
    READ(unit_index_CEF,*,IOSTAT = err) dummy
    IF (err == -1 ) THEN
      EXIT
    END IF 
    ERR = 0
    ncounts = ncounts + 1 
  END DO
  CLOSE (unit_index_CEF)
  nIntPoints = ncounts
  ! allocate needed space
  ALLOCATE(CurvedExternalField(1:2,1:ncounts))
  OPEN(UNIT=unit_index_CEF,FILE=FileNameCurvedExternalField,STATUS='OLD',FORM='FORMATTED')
  DO ii = 1, ncounts
    read(unit_index_CEF,'(F8.5,x,F8.5)') CurvedExternalField(1,ii) , CurvedExternalField(2,ii)
    IF (ii.GE.2) THEN
      diff_comp  = CurvedExternalField(1,2)  - CurvedExternalField(1,1)
      diff_check = CurvedExternalField(1,ii) - CurvedExternalField(1,ii-1)
      IF(ABS(diff_comp-diff_check).GE.1E-10)THEN
        SWRITE(UNIT_stdOut,'(A)') "ERROR: No equidistant points were used." 
        SWRITE(UNIT_stdOut,'(A)') diff_comp, diff_check
        STOP
      END IF  
    END IF
  END DO
  CLOSE (unit_index_CEF)

  IF (CurvedExternalField(1,1) .NE.0) THEN
    CALL abort(__STAMP__,  &
        "ERROR: Points have to start at 0.")
  END IF
  IF(ncounts.GT.1) THEN
    DeltaExternalField = CurvedExternalField(1,2)  - CurvedExternalField(1,1)
    SWRITE(UNIT_stdOut,'(A,1X,F8.5)') ' Delta external field: ',DeltaExternalField
    IF(DeltaExternalField.LE.0) THEN
      SWRITE(*,'(A)') ' ERROR: wrong sign in external field delta-x'
    END IF
  ELSE 
    CALL abort(__STAMP__, &
        " ERROR: not enough data points in curved external field file!")
  END IF
  SWRITE(UNIT_stdOut,'(A,I4.0,A)')'Found', ncounts,' data points.'
  SWRITE(UNIT_stdOut,'(A)')'...CURVED EXTERNAL FIELD INITIALIZATION DONE'
END SUBROUTINE read_curved_external_Field


FUNCTION InterpolateCurvedExternalField(Pos)
!===================================================================================================================================
! interpolates curved external field to z position
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PICInterpolation_Vars   ,ONLY:DeltaExternalField,nIntPoints,CurvedExternalField
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)          :: Pos ! partilce z position
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                     :: InterpolateCurvedExternalField  ! Bz
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                  :: iPos
!===================================================================================================================================

iPos = 0
iPos = INT(POS/DeltaExternalField) + 1
IF (iPos.GE.nIntPoints) THEN
  IPWRITE(UNIT_stdOut,'(A,F8.5,I10.2)')"Position and Position index, ",POS,iPos
  CALL abort(__STAMP__, &
        "ERROR: particle out of data point region for external curved field interpolation!")
END IF
!  Linear Interpolation between iPos and iPos+1 B point
InterpolateCurvedExternalField = (CurvedExternalField(2,iPos+1) - CurvedExternalField(2,iPos)) &
                               / (CurvedExternalField(1,iPos+1) - CurvedExternalField(1,iPos)) &
                               * (Pos - CurvedExternalField(1,iPos) ) + CurvedExternalField(2,iPos)

END FUNCTION InterpolateCurvedExternalField 

END MODULE MOD_PICInterpolation
