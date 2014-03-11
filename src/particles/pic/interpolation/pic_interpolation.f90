#include "boltzplatz.h"

MODULE  MOD_PICInterpolation                                                                       !
!===================================================================================================================================
! 
!===================================================================================================================================
   IMPLICIT NONE                                                                                   !
   PRIVATE                                                                                         !
!----------------------------------------------------------------------------------------------------------------------------------
   PUBLIC :: InterpolateFieldToParticle, InitializeInterpolation , Calc_inv, &
             InterpolateCurvedExternalField                                  !
!===================================================================================================================================

CONTAINS                                                                                           !
                                                                                                   !
SUBROUTINE InitializeInterpolation
!===================================================================================================================================
! Initialize the interpolation variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals!,                ONLY : UNIT_errOut, UNIT_StdOut
USE MOD_ReadInTools
USE MOD_Particle_Vars,          ONLY : PDM, GEO
USE MOD_PICInterpolation_Vars
USE MOD_Mesh_Vars,              ONLY : nElems
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: ALLOCSTAT
REAL                      :: P(3,8), PT(3,8), T(3,3), T_inv(3,3), DP(3)
INTEGER                   :: iNode, iElem, i, j, k
!===================================================================================================================================
  InterpolationType = GETSTR('PIC-Interpolation-Type','particle_position')
  externalField(1:6)= GETREALARRAY('PIC-externalField',6,'0.,0.,0.,0.,0.,0.')
  DoInterpolation   = GETLOGICAL('PIC-DoInterpolation','.TRUE.')
  useVTKFileEField  = GETLOGICAL('BGEField-UseVTKFile','.FALSE.')
  useVTKFileBField  = GETLOGICAL('BGBField-UseVTKFile','.FALSE.')
  
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
    WRITE(*,*)'ERROR in pic_interpolation.f90: Cannot allocate FieldAtParticle array!',ALLOCSTAT
    STOP
  END IF

  SELECT CASE(TRIM(InterpolationType))
  CASE('nearest_blurycenter')
     InterpolationType='nearest_blurrycenter'
  CASE('nearest_blurrycenter')
  CASE('particle_position_slow')
  CASE('particle_position')
    SWRITE(UNIT_StdOut,'(132("-"))')
    SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE INTERPOLATION...'
    ! prebuild trafo matrices
!    ALLOCATE(ElemPT(1:3,1:8,1:nElems),       &
!             ElemT_inv(1:3,1:3,1:nElems),     STAT=ALLOCSTAT)
    ALLOCATE(ElemT_inv(1:3,1:3,1:nElems),     STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) THEN
      WRITE(*,*) 'ERROR in InterpolationInit: Cannot allocate ElemPT!'
      STOP
    END IF
    ! prepare interpolation: calculate trafo matrices
    DO iElem = 1,nElems
      DO iNode = 1,8
        P(1:3,iNode) = GEO%NodeCoords(1:3,GEO%ElemToNodeID(iNode,iElem))
      END DO
      T(:,1) = 0.5 * (P(:,2)-P(:,1))
      T(:,2) = 0.5 * (P(:,4)-P(:,1))
      T(:,3) = 0.5 * (P(:,5)-P(:,1))
      T_inv = Calc_inv(T)
      ElemT_inv(1:3,1:3,iElem) = T_inv
    END DO

    SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE INTERPOLATION DONE!'
    SWRITE(UNIT_StdOut,'(132("-"))')
  CASE('nearest_gausspoint')
    SWRITE(UNIT_StdOut,'(132("-"))')
    SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE INTERPOLATION...'
    ! prebuild trafo matrices
!    ALLOCATE(ElemPT(1:3,1:8,1:nElems),       &
!             ElemT_inv(1:3,1:3,1:nElems),     STAT=ALLOCSTAT)
    ALLOCATE(ElemT_inv(1:3,1:3,1:nElems),     STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) THEN
      WRITE(*,*) 'ERROR in InterpolationInit: Cannot allocate ElemPT!'
      STOP
    END IF
    ! prepare interpolation: calculate trafo matrices
    DO iElem = 1,nElems
      DO iNode = 1,8
        P(1:3,iNode) = GEO%NodeCoords(1:3,GEO%ElemToNodeID(iNode,iElem))
      END DO
      T(:,1) = 0.5 * (P(:,2)-P(:,1))
      T(:,2) = 0.5 * (P(:,4)-P(:,1))
      T(:,3) = 0.5 * (P(:,5)-P(:,1))
      T_inv = Calc_inv(T)
      ElemT_inv(1:3,1:3,iElem) = T_inv
    END DO

    SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE INTERPOLATION DONE!'
    SWRITE(UNIT_StdOut,'(132("-"))')
  CASE DEFAULT
    !ERRWRITE(*,*) 'Unknown InterpolationType in pic_init.f90'
    WRITE(*,*) 'Unknown InterpolationType in pic_init.f90'
    STOP
  END SELECT
END SUBROUTINE InitializeInterpolation

SUBROUTINE InterpolateFieldToParticle()                                                            !
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Globals,       ONLY : UNIT_errOut
USE MOD_DG_Vars,       ONLY : U
USE MOD_Particle_Vars!, ONLY: 
USE MOD_PIC_Vars!,      ONLY: 
USE MOD_PICInterpolation_Vars
USE MOD_PICDepo_Vars,  ONLY : DepositionType
USE MOD_PreProc
USE MOD_Eval_xyz
#ifdef PP_POIS
USE MOD_Equation_Vars,ONLY: E
#endif
!--------------------------------------------------------------------------------------------------!
  IMPLICIT NONE                                                                                    !
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION                                                                        !
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION                                                                       !
  REAL                             :: Pos(3)                                                       !
  REAL                             :: field(6)                                                     !
  INTEGER                          :: i,j,k,m,iPart,iElem, iNode                                   !
  REAL                             :: dist_node(8), sum_dist
  REAL                             :: HelperU(1:6,0:PP_N,0:PP_N,0:PP_N)
!===================================================================================================================================

IF(usecurvedExternalField) THEN ! used curved external Bz
  FieldAtParticle(1:PDM%ParticleVecLength,:) = 0.
  FieldAtParticle(1:PDM%ParticleVecLength,1) = externalField(1)
  FieldAtParticle(1:PDM%ParticleVecLength,2) = externalField(2)
  FieldAtParticle(1:PDM%ParticleVecLength,3) = externalField(3)
#if (PP_nVar==8)
  FieldAtParticle(1:PDM%ParticleVecLength,4) = externalField(4)
  FieldAtParticle(1:PDM%ParticleVecLength,5) = externalField(5)
#endif
  ! Bz field strength at particle position
  DO iPart = 1, PDM%ParticleVecLength
    FieldAtparticle(iPart,6) = InterpolateCurvedExternalField(PartState(iPart,3))
  END DO
ELSE ! usecurvedExternalField
  FieldAtParticle(1:PDM%ParticleVecLength,:) = 0.
  FieldAtParticle(1:PDM%ParticleVecLength,1) = externalField(1)
  FieldAtParticle(1:PDM%ParticleVecLength,2) = externalField(2)
  FieldAtParticle(1:PDM%ParticleVecLength,3) = externalField(3)
#if (PP_nVar==8)
  FieldAtParticle(1:PDM%ParticleVecLength,4) = externalField(4)
  FieldAtParticle(1:PDM%ParticleVecLength,5) = externalField(5)
  FieldAtParticle(1:PDM%ParticleVecLength,6) = externalField(6)
#endif
END IF ! use constant external field

!temp: fake Efield of a particle at 000
!FieldAtParticle(1:PDM%ParticleVecLength,1)=0.218602365*PartState(1,1)/ &
!(SQRT(PartState(1,1)**2+PartState(1,2)**2+PartState(1,3)**2))**3
!FieldAtParticle(1:PDM%ParticleVecLength,2)=0.218602365*PartState(1,2)/ &
!(SQRT(PartState(1,1)**2+PartState(1,2)**2+PartState(1,3)**2))**3
!FieldAtParticle(1:PDM%ParticleVecLength,3)=0.218602365*PartState(1,3)/ &
!(SQRT(PartState(1,1)**2+PartState(1,2)**2+PartState(1,3)**2))**3


#if (PP_nVar==8)
IF (useVTKFileEField.AND.useVTKFileBField) THEN
  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
      iElem = PEM%Element(iPart)
      DO iNode = 1,8
        dist_node(iNode) = (GEO%NodeCoords(1,GEO%ElemToNodeID(iNode,iElem))-PartState(iPart,1))**2 + &
                           (GEO%NodeCoords(2,GEO%ElemToNodeID(iNode,iElem))-PartState(iPart,2))**2 + &
                           (GEO%NodeCoords(3,GEO%ElemToNodeID(iNode,iElem))-PartState(iPart,3))**2
      END DO
      dist_node = SQRT(dist_node)
      sum_dist  = SUM(dist_node)
      dist_node = 1.0-dist_node/sum_dist
      DO iNode = 1,8
        FieldAtParticle(iPart,1) = FieldAtParticle(iPart,1) + dist_node(iNode) * BGEfieldAtNode(1,GEO%ElemToNodeID(iNode,iElem))
        FieldAtParticle(iPart,2) = FieldAtParticle(iPart,2) + dist_node(iNode) * BGEfieldAtNode(2,GEO%ElemToNodeID(iNode,iElem))
        FieldAtParticle(iPart,3) = FieldAtParticle(iPart,3) + dist_node(iNode) * BGEfieldAtNode(3,GEO%ElemToNodeID(iNode,iElem))
        FieldAtParticle(iPart,4) = FieldAtParticle(iPart,4) + dist_node(iNode) * BGBfieldAtNode(1,GEO%ElemToNodeID(iNode,iElem))
        FieldAtParticle(iPart,5) = FieldAtParticle(iPart,5) + dist_node(iNode) * BGBfieldAtNode(2,GEO%ElemToNodeID(iNode,iElem))
        FieldAtParticle(iPart,6) = FieldAtParticle(iPart,6) + dist_node(iNode) * BGBfieldAtNode(3,GEO%ElemToNodeID(iNode,iElem))
      END DO
    END IF
  END DO
ELSE IF (useVTKFileEField) THEN
  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
      iElem = PEM%Element(iPart)
      DO iNode = 1,8
        dist_node(iNode) = (GEO%NodeCoords(1,GEO%ElemToNodeID(iNode,iElem))-PartState(iPart,1))**2 + &
                           (GEO%NodeCoords(2,GEO%ElemToNodeID(iNode,iElem))-PartState(iPart,2))**2 + &
                           (GEO%NodeCoords(3,GEO%ElemToNodeID(iNode,iElem))-PartState(iPart,3))**2
      END DO
      dist_node = SQRT(dist_node)
      sum_dist  = SUM(dist_node)
      dist_node = 1.0-dist_node/sum_dist
      DO iNode = 1,8
        FieldAtParticle(iPart,1) = FieldAtParticle(iPart,1) + dist_node(iNode) * BGEfieldAtNode(1,GEO%ElemToNodeID(iNode,iElem))
        FieldAtParticle(iPart,2) = FieldAtParticle(iPart,2) + dist_node(iNode) * BGEfieldAtNode(2,GEO%ElemToNodeID(iNode,iElem))
        FieldAtParticle(iPart,3) = FieldAtParticle(iPart,3) + dist_node(iNode) * BGEfieldAtNode(3,GEO%ElemToNodeID(iNode,iElem))
      END DO
    END IF
  END DO
ELSE IF (useVTKFileBField) THEN
  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
      iElem = PEM%Element(iPart)
      DO iNode = 1,8
        dist_node(iNode) = (GEO%NodeCoords(1,GEO%ElemToNodeID(iNode,iElem))-PartState(iPart,1))**2 + &
                           (GEO%NodeCoords(2,GEO%ElemToNodeID(iNode,iElem))-PartState(iPart,2))**2 + &
                           (GEO%NodeCoords(3,GEO%ElemToNodeID(iNode,iElem))-PartState(iPart,3))**2
      END DO
      dist_node = SQRT(dist_node)
      sum_dist  = SUM(dist_node)
      dist_node = 1.0 - dist_node/sum_dist
      DO iNode = 1,8
        FieldAtParticle(iPart,4) = FieldAtParticle(iPart,4) + dist_node(iNode) * BGBfieldAtNode(1,GEO%ElemToNodeID(iNode,iElem))
        FieldAtParticle(iPart,5) = FieldAtParticle(iPart,5) + dist_node(iNode) * BGBfieldAtNode(2,GEO%ElemToNodeID(iNode,iElem))
        FieldAtParticle(iPart,6) = FieldAtParticle(iPart,6) + dist_node(iNode) * BGBfieldAtNode(3,GEO%ElemToNodeID(iNode,iElem))
      END DO
    END IF
  END DO
END IF
#else
IF (useVTKFileEField) THEN
  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
      iElem = PEM%Element(iPart)
      DO iNode = 1,8
        dist_node(iNode) = (GEO%NodeCoords(1,GEO%ElemToNodeID(iNode,iElem))-PartState(iPart,1))**2 + &
                           (GEO%NodeCoords(2,GEO%ElemToNodeID(iNode,iElem))-PartState(iPart,2))**2 + &
                           (GEO%NodeCoords(3,GEO%ElemToNodeID(iNode,iElem))-PartState(iPart,3))**2
      END DO
      dist_node = SQRT(dist_node)
      sum_dist  = SUM(dist_node)
      dist_node = 1.0-dist_node/sum_dist
      DO iNode = 1,8
        FieldAtParticle(iPart,1) = FieldAtParticle(iPart,1) + dist_node(iNode) * BGEfieldAtNode(1,GEO%ElemToNodeID(iNode,iElem))
        FieldAtParticle(iPart,2) = FieldAtParticle(iPart,2) + dist_node(iNode) * BGEfieldAtNode(2,GEO%ElemToNodeID(iNode,iElem))
        FieldAtParticle(iPart,3) = FieldAtParticle(iPart,3) + dist_node(iNode) * BGEfieldAtNode(3,GEO%ElemToNodeID(iNode,iElem))
      END DO
    END IF
  END DO
END IF
#endif
IF (DoInterpolation) THEN                 ! skip if no self fields are calculated
  SELECT CASE(TRIM(InterpolationType))
  CASE('nearest_blurrycenter')
    ! add fields to fields at particle position
    DO iPart=1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        iElem = PEM%Element(iPart)
        m = INT(PP_N/2)+1
#ifdef PP_POIS
        FieldAtParticle(iPart,1) = FieldAtParticle(iPart,1) + E(1,m,m,m,iElem) 
        FieldAtParticle(iPart,2) = FieldAtParticle(iPart,2) + E(2,m,m,m,iElem) 
        FieldAtParticle(iPart,3) = FieldAtParticle(iPart,3) + E(3,m,m,m,iElem) 
#else
        FieldAtParticle(iPart,1) = FieldAtParticle(iPart,1) + U(1,m,m,m,iElem) 
        FieldAtParticle(iPart,2) = FieldAtParticle(iPart,2) + U(2,m,m,m,iElem) 
        FieldAtParticle(iPart,3) = FieldAtParticle(iPart,3) + U(3,m,m,m,iElem) 
#endif
#if (PP_nVar==8)
        FieldAtParticle(iPart,4) = FieldAtParticle(iPart,4) + U(4,m,m,m,iElem) 
        FieldAtParticle(iPart,5) = FieldAtParticle(iPart,5) + U(5,m,m,m,iElem) 
        FieldAtParticle(iPart,6) = FieldAtParticle(iPart,6) + U(6,m,m,m,iElem) 
#endif
      END IF
    END DO
  CASE('particle_position_slow')
    DO iPart = 1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        iElem = PEM%Element(iPart)
        Pos = PartState(iPart,1:3)
        !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
        HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
        HelperU(4:6,:,:,:) = U(4:6,:,:,:,iElem)
        CALL eval_xyz(Pos,6,PP_N,HelperU,field,iElem)
#else
        CALL eval_xyz(Pos,6,PP_N,U(1:6,:,:,:,iElem),field,iElem)
#endif
#else
#ifdef PP_POIS
        CALL eval_xyz(Pos,3,PP_N,E(1:3,:,:,:,iElem),field,iElem)
#else
        CALL eval_xyz(Pos,3,PP_N,U(1:3,:,:,:,iElem),field,iElem)
#endif
#endif
        FieldAtParticle(iPart,:) = FieldAtParticle(iPart,:) + field
      END IF
    END DO
  CASE('particle_position')
    IF(TRIM(DepositionType).EQ.'nearest_gausspoint')THEN
      ! particles have already been mapped in deposition, other eval routine used
      DO iPart = 1,PDM%ParticleVecLength
        IF (PDM%ParticleInside(iPart)) THEN
          iElem = PEM%Element(iPart)
          !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
          HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
          HelperU(4:6,:,:,:) = U(4:6,:,:,:,iElem)
          CALL eval_xyz_part2(PartPosMapped(iPart,1:3),6,PP_N,HelperU,field,iElem)
#else
          CALL eval_xyz_part2(PartPosMapped(iPart,1:3),6,PP_N,U(1:6,:,:,:,iElem),field,iElem)
#endif
#else
#ifdef PP_POIS
          CALL eval_xyz_part2(PartPosMapped(iPart,1:3),3,PP_N,E(1:3,:,:,:,iElem),field,iElem)     
#else
          CALL eval_xyz_part2(PartPosMapped(iPart,1:3),3,PP_N,U(1:3,:,:,:,iElem),field,iElem)
#endif
#endif
          FieldAtParticle(iPart,:) = FieldAtParticle(iPart,:) + field
        END IF
      END DO
    ELSE ! particles are not yet mapped
      DO iPart = 1,PDM%ParticleVecLength
        IF (PDM%ParticleInside(iPart)) THEN
          iElem = PEM%Element(iPart)
          Pos = PartState(iPart,1:3)
          !--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
          HelperU(1:3,:,:,:) = E(1:3,:,:,:,iElem)
          HelperU(4:6,:,:,:) = U(4:6,:,:,:,iElem)
          CALL eval_xyz_fast(Pos,6,PP_N,HelperU,field,iElem)
#else
          CALL eval_xyz_fast(Pos,6,PP_N,U(1:6,:,:,:,iElem),field,iElem)
#endif
#else
#ifdef PP_POIS
         CALL eval_xyz_fast(Pos,3,PP_N,E(1:3,:,:,:,iElem),field,iElem)
#else
         CALL eval_xyz_fast(Pos,3,PP_N,U(1:3,:,:,:,iElem),field,iElem)
#endif         
#endif
          FieldAtParticle(iPart,:) = FieldAtParticle(iPart,:) + field
        END IF
      END DO
    END IF
  CASE('nearest_gausspoint')
      ! particles have already been mapped in deposition
      DO iPart = 1,PDM%ParticleVecLength
        IF (PDM%ParticleInside(iPart)) THEN
          iElem = PEM%Element(iPart)
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
        END IF
      END DO
  CASE DEFAULT
    ERRWRITE(*,*) 'ERROR: Unknown InterpolationType!'
    STOP
  END SELECT
END IF

!! Temp stuff for pseudo recordlines
!IF(MOD(iter,1000).EQ.0) THEN
!  DO iPart = 1,PDM%ParticleVecLength
!    IF (PDM%ParticleInside(iPart)) THEN
!      WRITE(*,'(I4,6(X,E14.6))') ipart, PartState(iPart,1),PartState(iPart,2),PartState(iPart,3), &
!         FieldAtParticle(iPart,4),FieldAtParticle(iPart,5),FieldAtParticle(iPart,6)
!    END IF
!  END DO
!END IF
    
RETURN
END SUBROUTINE InterpolateFieldToParticle


FUNCTION Calc_inv(M)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)          :: M(3,3)     ! 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                     :: Calc_inv(3,3)  !  
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                     :: M_inv(3,3)  !  
REAL                     :: detjb
!===================================================================================================================================
M_inv = 0.

! Determines the determinant of xj and checks for zero values
!
detjb = M (1, 1) * M (2, 2) * M (3, 3) &
      + M (1, 2) * M (2, 3) * M (3, 1) &
      + M (1, 3) * M (2, 1) * M (3, 2) &
      - M (1, 3) * M (2, 2) * M (3, 1) &
      - M (1, 2) * M (2, 1) * M (3, 3) &
      - M (1, 1) * M (2, 3) * M (3, 2)
IF ( detjb <= 0.d0 ) then
  WRITE(*,*)"Negative determinant of Jacobian in M_inv"
  WRITE(*,*)"Determinant is:",detjb
  WRITE(*,*)"KM:",M_inv
  STOP 
END IF
!
! Determines the inverse of xj
!
M_inv (1, 1) = (M (2, 2) * M (3, 3) &
              - M (2, 3) * M (3, 2) ) / detjb
M_inv (1, 2) = (M (1, 3) * M (3, 2) &
              - M (1, 2) * M (3, 3) ) / detjb
M_inv (1, 3) = (M (1, 2) * M (2, 3) &
              - M (1, 3) * M (2, 2) ) / detjb
M_inv (2, 1) = (M (2, 3) * M (3, 1) &
              - M (2, 1) * M (3, 3) ) / detjb
M_inv (2, 2) = (M (1, 1) * M (3, 3) &
              - M (1, 3) * M (3, 1) ) / detjb
M_inv (2, 3) = (M (1, 3) * M (2, 1) &
              - M (1, 1) * M (2, 3) ) / detjb
M_inv (3, 1) = (M (2, 1) * M (3, 2) &
              - M (2, 2) * M (3, 1) ) / detjb
M_inv (3, 2) = (M (1, 2) * M (3, 1) &
              - M (1, 1) * M (3, 2) ) / detjb
M_inv (3, 3) = (M (1, 1) * M (2, 2) &
              - M (1, 2) * M (2, 1) ) / detjb

Calc_inv = M_inv
END FUNCTION Calc_inv 

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
    SWRITE(UNIT_stdOut,'(A)') "ERROR: Points have to start at 0." 
    STOP
  END IF
  IF(ncounts.GT.1) THEN
    DeltaExternalField = CurvedExternalField(1,2)  - CurvedExternalField(1,1)
    SWRITE(UNIT_stdOut,'(A,1X,F8.5)') ' Delta external field: ',DeltaExternalField
    IF(DeltaExternalField.LE.0) THEN
      SWRITE(*,'(A)') ' ERROR: wrong sign in external field delta-x'
    END IF
  ELSE 
    SWRITE(*,'(A)') " ERROR: not enough data points in curved external field file!"
    STOP
  END IF
  SWRITE(UNIT_stdOut,'(A,I4.0,A)')'Found', ncounts,' data points.'
  SWRITE(UNIT_stdOut,'(A)')'...CURVED EXTERNAL FIELD INITIALIZATION DONE'
END SUBROUTINE read_curved_external_Field

FUNCTION InterpolateCurvedExternalField(Pos)
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
INTEGER                  :: intPolPoint, iPos
!===================================================================================================================================

iPos = 0
iPos = INT(POS/DeltaExternalField) + 1
IF (iPos.GE.nIntPoints) THEN
  SWRITE(*,'(A)')"ERROR: particle out of data point region for external field interpolation!"
  SWRITE(*,'(A,F8.5,I10.2)')"Position and Position index, ",POS,iPos
  STOP
END IF
!  Linear Interpolation between iPos and iPos+1 B point
InterpolateCurvedExternalField = (CurvedExternalField(2,iPos+1) - CurvedExternalField(2,iPos)) &
                               / (CurvedExternalField(1,iPos+1) - CurvedExternalField(1,iPos)) &
                               * (Pos - CurvedExternalField(1,iPos) ) + CurvedExternalField(2,iPos)

END FUNCTION InterpolateCurvedExternalField 




END MODULE MOD_PICInterpolation
