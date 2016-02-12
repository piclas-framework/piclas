#include "boltzplatz.h"

MODULE  MOD_PICDepo                                                                                
!===================================================================================================================================
! MOD PIC Depo
!===================================================================================================================================
 IMPLICIT NONE                                                                                   
 PRIVATE                                                                                         
!===================================================================================================================================
INTERFACE Deposition
  MODULE PROCEDURE Deposition
END INTERFACE

INTERFACE InitializeDeposition
  MODULE PROCEDURE InitializeDeposition
END INTERFACE

INTERFACE FinalizeDeposition
  MODULE PROCEDURE FinalizeDeposition
END INTERFACE

PUBLIC:: Deposition, InitializeDeposition, FinalizeDeposition
!INTERFACE DepositionMPF
!  MODULE PROCEDURE DepositionMPF
!END INTERFACE
!===================================================================================================================================

CONTAINS                                                                                           
                                                                                                   
SUBROUTINE InitializeDeposition
!===================================================================================================================================
! Initialize the deposition variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PICDepo_Vars!,  ONLY : DepositionType, source, r_sf, w_sf, r2_sf, r2_sf_inv
USE MOD_Particle_Vars ! crazy??
USE MOD_Globals_Vars,           ONLY:PI
USE MOD_Mesh_Vars,              ONLY:nElems, XCL_NGeo,Elem_xGP, sJ,nGlobalElems
USE MOD_Particle_Mesh_Vars,     ONLY:Geo
USE MOD_Interpolation_Vars,     ONLY:xGP,wBary
USE MOD_Basis,                  ONLY:ComputeBernsteinCoeff
USE MOD_Basis,                  ONLY:BarycentricWeights,InitializeVandermonde
USE MOD_Basis,                  ONLY:LegendreGaussNodesAndWeights,LegGaussLobNodesAndWeights
USE MOD_ChangeBasis,            ONLY:ChangeBasis3D
USE MOD_PreProc,                ONLY:PP_N,PP_nElems
USE MOD_ReadInTools,            ONLY:GETREAL,GETINT,GETLOGICAL,GETSTR,GETREALARRAY
USE MOD_PICInterpolation_Vars,  ONLY:InterpolationType
USE MOD_Eval_xyz,               ONLY:eval_xyz_elemcheck
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping
#ifdef MPI
USE MOD_Particle_MPI_Vars,      ONLY:DoExternalParts
#endif
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE          :: xGP_tmp(:),wBary_tmp(:),wGP_tmp(:),Vdm_GaussN_EquiN(:,:), Vdm_GaussN_NDepo(:,:)
REAL,ALLOCATABLE          :: dummy(:,:,:,:),dummy2(:,:,:,:)
INTEGER                   :: ALLOCSTAT, iElem, i, j, k, m, dir, weightrun, mm, r, s, t
REAL                      :: Temp(3), MappedGauss(1:PP_N+1), xmin, ymin, zmin, xmax, ymax, zmax, x0
REAL                      :: auxiliary(0:3),weight(1:3,0:3), nTotalDOF, VolumeShapeFunction
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE DEPOSITION...'

DepositionType = GETSTR('PIC-Deposition-Type','nearest_blurrycenter')
! check for interpolation type incompatibilities (cannot be done at interpolation_init
! because DepositionType is not known yet)
IF((TRIM(InterpolationType).EQ.'nearest_gausspoint').AND. &
   (TRIM(DepositionType).NE.'nearest_gausspoint')) THEN
  CALL abort(__STAMP__, &
    'ERROR in pic_depo.f90: Interpolation type nearest_gausspoint only allowed with same deposition type!')
END IF
!--- Allocate arrays for charge density collection and initialize
SDEALLOCATE(source)
ALLOCATE(source(1:4,0:PP_N,0:PP_N,0:PP_N,nElems),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(__STAMP__, &
      'ERROR in pic_depo.f90: Cannot allocate source!')
END IF
Source=0.
SELECT CASE(TRIM(DepositionType))
CASE('nearest_blurrycenter')
CASE('nearest_blurycenter')
  DepositionType = 'nearest_blurrycenter'
CASE('nearest_gausspoint')
  ! Allocate array for particle positions in -1|1 space (used for deposition as well as interpolation)
  ! only if NOT DoRefMapping
  IF(.NOT.DoRefMapping)THEN
    ALLOCATE(PartPosRef(1:3,PDM%MaxParticleNumber), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,' Cannot allocate partposref!')
  END IF
  ! compute the borders of the virtual volumes around the gauss points in -1|1 space
  ALLOCATE(GaussBorder(1:PP_N),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(__STAMP__, &
      'ERROR in pic_depo.f90: Cannot allocate Mapped Gauss Border Coords!')
  END IF
  DO i=0,PP_N
    ! bullshit here, use xGP
    !CALL Eval_XYZ_ElemCheck(Elem_xGP(:,i,1,1,1),Temp(:),1)
    !MappedGauss(i+1) = Temp(1)
    MappedGauss(i+1) = xGP(i)
  END DO
  DO i = 1,PP_N
    GaussBorder(i) = (MappedGauss(i+1) + MappedGauss(i))/2
  END DO
  ! allocate array for saving the gauss points of particles for nearest_gausspoint interpolation
  IF(TRIM(InterpolationType).EQ.'nearest_gausspoint')THEN
    SDEALLOCATE(PartPosGauss)
    ALLOCATE(PartPosGauss(1:PDM%maxParticleNumber,1:3),STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) THEN
      CALL abort(__STAMP__, &
      'ERROR in pic_depo.f90: Cannot allocate Part Pos Gauss!')
    END IF
  END IF
CASE('shape_function')
  !IF(.NOT.DoRefMapping) CALL abort(&
  !  __STAMP__&
  !  ,' Shape function has to be used with ref element tracking.')
  !ALLOCATE(PartToFIBGM(1:6,1:PDM%maxParticleNumber),STAT=ALLOCSTAT)
  !IF (ALLOCSTAT.NE.0) CALL abort(&
  !    __STAMP__, &
  !    ' Cannot allocate PartToFIBGM!')
  !ALLOCATE(ExtPartToFIBGM(1:6,1:PDM%ParticleVecLength),STAT=ALLOCSTAT)
  !IF (ALLOCSTAT.NE.0) THEN
  !  CALL abort(__STAMP__, &
  !    ' Cannot allocate ExtPartToFIBGM!')
  r_sf     = GETREAL('PIC-shapefunction-radius','1.')
  alpha_sf = GETINT('PIC-shapefunction-alpha','2')
  DoSFEqui = GETLOGICAL('PIC-shapefunction-equi','F')
  BetaFac = beta(1.5, REAL(alpha_sf) + 1.)
  w_sf = 1./(2. * BetaFac * REAL(alpha_sf) + 2 * BetaFac) &
                        * (REAL(alpha_sf) + 1.)/(PI*(r_sf**3))
  r2_sf = r_sf * r_sf 
  r2_sf_inv = 1./r2_sf
  VolumeShapeFunction=4./3.*PI*r_sf*r2_sf
  nTotalDOF=REAL(nGlobalElems)*REAL((PP_N+1)**3)
  IF(MPIRoot)THEN
    IF(VolumeShapeFunction.GT.GEO%MeshVolume) &
      CALL abort(__STAMP__, &
      'ShapeFunctionVolume > MeshVolume')
  END IF
  
  SWRITE(UNIT_stdOut,'(A,F12.6)') ' | Average DOFs in Shape-Function |                      ' , &
      nTotalDOF*VolumeShapeFunction/GEO%MeshVolume

  ALLOCATE(ElemDepo_xGP(1:3,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__, &
      ' Cannot allocate ElemDepo_xGP!')
  IF(DoSFEqui)THEN
    ALLOCATE( Vdm_EquiN_GaussN(0:PP_N,0:PP_N)     &
            , Vdm_GaussN_EquiN(0:PP_N,0:PP_N)     &
            , wGP_tmp(0:PP_N)                     &
            , xGP_tmp(0:PP_N)                     &
            , wBary_tmp(0:PP_N)                   )
    DO i=0,PP_N
      xGP_tmp(i) = 2./REAL(PP_N) * REAL(i) - 1. 
    END DO
    CALL BarycentricWeights(PP_N,xGP_tmp,wBary_tmp)
    ! to Gauss
    CALL InitializeVandermonde(PP_N,PP_N,wBary_tmp,xGP_tmp,xGP ,Vdm_EquiN_GaussN)
    ! from Gauss
    CALL InitializeVandermonde(PP_N,PP_N,wBary,xGP,xGP_tmp ,Vdm_GaussN_EquiN)
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_GaussN_EquiN,Elem_xGP(:,:,:,:,iElem),ElemDepo_xGP(:,:,:,:,iElem))
    END DO ! iElem=1,PP_nElems
    DEALLOCATE( Vdm_GaussN_EquiN, wGP_tmp, xGP_tmp, wBary_tmp)
  ELSE
    ElemDepo_xGP=Elem_xGP
  END IF

#ifdef MPI
  DoExternalParts=.TRUE.
#endif /*MPI*/

CASE('shape_function_1d')
  r_sf     = GETREAL('PIC-shapefunction-radius','1.')
  alpha_sf = GETINT ('PIC-shapefunction-alpha','2')
  sf1d_dir = GETINT ('PIC-shapefunction1d-direction','1')
  DoSFEqui = GETLOGICAL('PIC-shapefunction-equi','F')
  r2_sf = r_sf * r_sf 
  r2_sf_inv = 1./r2_sf
  SELECT CASE(alpha_sf)
  CASE(2)
    w_sf=16.*r_sf/15.
  CASE(3)
    w_sf=32.*r_sf/35.
  CASE(4)
    w_sf=256.*r_sf/315.
  CASE(5)
    w_sf=512.*r_sf/693.
  CASE(6)
    w_sf=2048.*r_sf/3003.
  CASE(7)
    w_sf=4096.*r_sf/6435.
  CASE(8)
    w_sf=65536.*r_sf/109395.
  CASE DEFAULT
    CALL abort(__STAMP__, &
    ' Correct 1D weight not precomputed!')
  END SELECT

  IF(sf1d_dir.EQ.1)THEN
    w_sf=w_sf*(GEO%ymaxglob-GEO%yminglob)*(GEO%zmaxglob-GEO%zminglob)
  ELSE IF (sf1d_dir.EQ.2)THEN
    w_sf=w_sf*(GEO%xmaxglob-GEO%xminglob)*(GEO%zmaxglob-GEO%zminglob)
  ELSE IF (sf1d_dir.EQ.3)THEN
    w_sf=w_sf*(GEO%xmaxglob-GEO%xminglob)*(GEO%ymaxglob-GEO%yminglob)
  ELSE
    w_sf=2*GEO%MeshVolume
  END IF
  VolumeShapeFunction=w_sf
  w_sf=1.0/w_sf
  nTotalDOF=REAL(nGlobalElems)*REAL((PP_N+1)**3)
  IF(MPIRoot)THEN
    IF(VolumeShapeFunction.GT.GEO%MeshVolume) &
      CALL abort(__STAMP__, &
      'ShapeFunctionVolume > MeshVolume')
  END IF
  
  SWRITE(UNIT_stdOut,'(A,F12.6)') ' | Shape function volume          |                      ' , &
      VolumeShapeFunction
  SWRITE(UNIT_stdOut,'(A,F12.6)') ' | Average DOFs in Shape-Function |                      ' , &
      nTotalDOF*VolumeShapeFunction/GEO%MeshVolume

  ALLOCATE(ElemDepo_xGP(1:3,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__, &
      ' Cannot allocate ElemDepo_xGP!')
  IF(DoSFEqui)THEN
    ALLOCATE( Vdm_EquiN_GaussN(0:PP_N,0:PP_N)     &
            , Vdm_GaussN_EquiN(0:PP_N,0:PP_N)     &
            , wGP_tmp(0:PP_N)                     &
            , xGP_tmp(0:PP_N)                     &
            , wBary_tmp(0:PP_N)                   )
    DO i=0,PP_N
      xGP_tmp(i) = 2./REAL(PP_N) * REAL(i) - 1. 
    END DO
    CALL BarycentricWeights(PP_N,xGP_tmp,wBary_tmp)
    ! to Gauss
    CALL InitializeVandermonde(PP_N,PP_N,wBary_tmp,xGP_tmp,xGP ,Vdm_EquiN_GaussN)
    ! from Gauss
    CALL InitializeVandermonde(PP_N,PP_N,wBary,xGP,xGP_tmp ,Vdm_GaussN_EquiN)
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_GaussN_EquiN,Elem_xGP(:,:,:,:,iElem),ElemDepo_xGP(:,:,:,:,iElem))
    END DO ! iElem=1,PP_nElems
    DEALLOCATE( Vdm_GaussN_EquiN, wGP_tmp, xGP_tmp, wBary_tmp)
  ELSE
    ElemDepo_xGP=Elem_xGP
  END IF

#ifdef MPI
  DoExternalParts=.TRUE.
#endif /*MPI*/

CASE('cylindrical_shape_function')
  !IF(.NOT.DoRefMapping) CALL abort(&
  !  __STAMP__&
  !  ,' Shape function has to be used with ref element tracking.')
  !ALLOCATE(PartToFIBGM(1:6,1:PDM%maxParticleNumber),STAT=ALLOCSTAT)
  !IF (ALLOCSTAT.NE.0) CALL abort(&
  !    __STAMP__, &
  !    ' Cannot allocate PartToFIBGM!')
  !ALLOCATE(ExtPartToFIBGM(1:6,1:PDM%ParticleVecLength),STAT=ALLOCSTAT)
  !IF (ALLOCSTAT.NE.0) THEN
  !  CALL abort(__STAMP__, &
  !    ' Cannot allocate ExtPartToFIBGM!')
  r_sf       = GETREAL('PIC-shapefunction-radius','1.')
  r_sf0      = GETREAL('PIC-shapefunction-radius0','1.')
  r_sf_scale = GETREAL('PIC-shapefunction-scale','0.')
  alpha_sf   = GETINT('PIC-shapefunction-alpha','2')
  DoSFEqui   = GETLOGICAL('PIC-shapefunction-equi','F')
  BetaFac    = beta(1.5, REAL(alpha_sf) + 1.)
  w_sf = 1./(2. * BetaFac * REAL(alpha_sf) + 2 * BetaFac) &
                        * (REAL(alpha_sf) + 1.)!/(PI*(r_sf**3))
  r2_sf = r_sf * r_sf 
  r2_sf_inv = 1./r2_sf

  ALLOCATE(ElemDepo_xGP(1:3,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
      __STAMP__, &
      ' Cannot allocate ElemDepo_xGP!')
  IF(DoSFEqui)THEN
    ALLOCATE( Vdm_EquiN_GaussN(0:PP_N,0:PP_N)     &
            , Vdm_GaussN_EquiN(0:PP_N,0:PP_N)     &
            , wGP_tmp(0:PP_N)                     &
            , xGP_tmp(0:PP_N)                     &
            , wBary_tmp(0:PP_N)                   )
    DO i=0,PP_N
      xGP_tmp(i) = 2./REAL(PP_N) * REAL(i) - 1. 
    END DO
    CALL BarycentricWeights(PP_N,xGP_tmp,wBary_tmp)
    ! to Gauss
    CALL InitializeVandermonde(PP_N,PP_N,wBary_tmp,xGP_tmp,xGP ,Vdm_EquiN_GaussN)
    ! from Gauss
    CALL InitializeVandermonde(PP_N,PP_N,wBary,xGP,xGP_tmp ,Vdm_GaussN_EquiN)
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_GaussN_EquiN,Elem_xGP(:,:,:,:,iElem),ElemDepo_xGP(:,:,:,:,iElem))
    END DO ! iElem=1,PP_nElems
    DEALLOCATE( Vdm_GaussN_EquiN, wGP_tmp, xGP_tmp, wBary_tmp)
  ELSE
    ElemDepo_xGP=Elem_xGP
  END IF

#ifdef MPI
  DoExternalParts=.TRUE.
#endif /*MPI*/

CASE('delta_distri')
  ! Allocate array for particle positions in -1|1 space (used for deposition as well as interpolation)
  IF(.NOT.DoRefMapping)THEN
    ALLOCATE(PartPosRef(1:3,PDM%MaxParticleNumber), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,' Cannot allocate partposref!')
  END IF
  DeltaType = GETINT('PIC-DeltaType','3')
  NDepo     = GETINT('PIC-DeltaType-N','1')
  SELECT CASE(DeltaType)
  CASE(1)
    SWRITE(UNIT_stdOut,'(A)') ' Lagrange-Polynomial'
  CASE(2)
    SWRITE(UNIT_stdOut,'(A)') ' Bernstein-Polynomial'
    !CALL BuildBernSteinVdm(PP_N,xGP)
    ALLOCATE(NDepochooseK(0:NDepo,0:NDepo))
    CALL ComputeBernsteinCoeff(NDepo,NDepoChooseK)
  CASE(3)
    SWRITE(UNIT_stdOut,'(A)') ' Uniform B-Spline '
    NKnots=(NDepo+1)*2-1
    ALLOCATE( Knots(0:NKnots)  )
    ! knots distance is 2 [-1,1)
    X0=-1.-REAL(NDepo)*2.0
    DO i=0,nKnots
      Knots(i) =X0+2.0*REAL(i)
    END DO ! i
  CASE DEFAULT
    CALL abort(__STAMP__, &
        ' Wrong Delta-Type')
  END SELECT
  IF(NDepo.GT.PP_N) CALL abort(&
    __STAMP__&
    ,' NDepo must be smaller than N!')
  DoChangeBasis=.FALSE.
  IF(NDepo.NE.PP_N) DoChangeBasis=.TRUE.
  ALLOCATE( Vdm_NDepo_GaussN(0:PP_N,0:NDepo)             &
          , Vdm_GaussN_NDepo(0:NDepo,0:PP_N)             &
          , dummy(1,0:PP_N,0:PP_N,0:PP_N)                &
          , dummy2(1,0:NDepo,0:NDepo,0:NDepo)            &
          , XiNDepo(0:NDepo)                             &
          , swGPNDepo(0:NDepo)                           &
          , sjNDepo(0:NDepo,0:NDepo,0:NDepo,0:PP_nElems) &
          , wBaryNDepo(0:NDepo)                          )
#if (PP_NodeType==1)
    CALL LegendreGaussNodesAndWeights(NDepo,XiNDepo,swGPNDepo)
#elif (PP_NodeType==2)
    CALL LegGaussLobNodesAndWeights(NDepo,XiNDepo,swGPNDepo)
#endif
  CALL BarycentricWeights(NDepo,XiNDepo,wBaryNDepo)
  IF(DoChangeBasis)THEN
    CALL InitializeVandermonde(NDepo,PP_N,wBaryNDepo,XiNDepo,xGP,Vdm_NDepo_GaussN)
    !CALL InitializeVandermonde(N_Restart_in,N_in,wBary_Restart,xGP_Restart,xGP,Vdm_GaussNRestart_GaussN)
  ELSE
    DEALLOCATE(Vdm_NDepo_GaussN)
  END IF
  ! and inverse
  DO i=0,NDepo
    swGPNDepo(i)=1.0/swGPNDepo(i)
  END DO ! i=0,PP_N
  CALL InitializeVandermonde(PP_N,NDepo,wBary,xGP,xiNDepo,Vdm_GaussN_NDepo)
  DO iElem=1,PP_nElems
    !CALL ChangeBasis3D(1,PP_N,NDepo,Vdm_GaussN_NDepo,sJ(:,:,:,iElem),sjNDepo(:,:,:,iElem))
    dummy(1,:,:,:)=sJ(:,:,:,iElem)
    CALL ChangeBasis3D(1,PP_N,NDepo,Vdm_GaussN_NDepo,dummy,dummy2)
    sJNDepo(:,:,:,iElem) = dummy2(1,:,:,:)
  END DO ! iElem=1,PP_nElems
  DEALLOCATE(Vdm_GaussN_NDepo)
  DEALLOCATE(dummy,dummy2)

CASE('cartmesh_volumeweighting')

  ! read in background mesh size
  BGMdeltas(1:3) = GETREALARRAY('PIC-BGMdeltas',3,'0. , 0. , 0.')
  FactorBGM(1:3) = GETREALARRAY('PIC-FactorBGM',3,'1. , 1. , 1.')
  BGMdeltas(1:3) = 1./FactorBGM(1:3)*BGMdeltas(1:3)
  IF (ANY(BGMdeltas.EQ.0.0)) THEN
    CALL abort(__STAMP__, &
      'ERROR: PIC-BGMdeltas: No size for the cartesian background mesh definded.')
  END IF
  ! calc min and max coordinates for local mesh
  xmin = 1.0E200
  ymin = 1.0E200
  zmin = 1.0E200
  xmax = -1.0E200
  ymax = -1.0E200
  zmax = -1.0E200
  DO iElem = 1, nElems
    xmin=MIN(xmin,MINVAL(XCL_NGeo(1,:,:,:,iElem)))
    xmax=MAX(xmax,MAXVAL(XCL_NGeo(1,:,:,:,iElem)))
    ymin=MIN(ymin,MINVAL(XCL_NGeo(2,:,:,:,iElem)))
    ymax=MAX(ymax,MAXVAL(XCL_NGeo(2,:,:,:,iElem)))
    zmin=MIN(zmin,MINVAL(XCL_NGeo(3,:,:,:,iElem)))
    zmax=MAX(zmax,MAXVAL(XCL_NGeo(3,:,:,:,iElem)))
  END DO
  ! define minimum and maximum backgroundmesh index, compute volume
  BGMVolume = BGMdeltas(1)*BGMdeltas(2)*BGMdeltas(3)
  BGMminX = FLOOR(xmin/BGMdeltas(1)-0.0001)
  BGMminY = FLOOR(ymin/BGMdeltas(2)-0.0001)
  BGMminZ = FLOOR(zmin/BGMdeltas(3)-0.0001)
  BGMmaxX = CEILING(xmax/BGMdeltas(1)+0.0001)
  BGMmaxY = CEILING(ymax/BGMdeltas(2)+0.0001)
  BGMmaxZ = CEILING(zmax/BGMdeltas(3)+0.0001)
  ! mapping from gausspoints to BGM
  ALLOCATE(GaussBGMIndex(1:3,0:PP_N,0:PP_N,0:PP_N,1:nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(__STAMP__, &
      'ERROR in pic_depo.f90: Cannot allocate GaussBGMIndex!') 
  END IF
  ALLOCATE(GaussBGMFactor(1:3,0:PP_N,0:PP_N,0:PP_N,1:nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(__STAMP__, &
      'ERROR in pic_depo.f90: Cannot allocate GaussBGMFactor!')
  END IF
  DO iElem = 1, nElems
    DO j = 0, PP_N
      DO k = 0, PP_N
        DO m = 0, PP_N
          GaussBGMIndex(1,j,k,m,iElem) = FLOOR(Elem_xGP(1,j,k,m,iElem)/BGMdeltas(1))
          GaussBGMIndex(2,j,k,m,iElem) = FLOOR(Elem_xGP(2,j,k,m,iElem)/BGMdeltas(2))
          GaussBGMIndex(3,j,k,m,iElem) = FLOOR(Elem_xGP(3,j,k,m,iElem)/BGMdeltas(3))
          GaussBGMFactor(1,j,k,m,iElem) = (Elem_xGP(1,j,k,m,iElem)/BGMdeltas(1))-REAL(GaussBGMIndex(1,j,k,m,iElem))
          GaussBGMFactor(2,j,k,m,iElem) = (Elem_xGP(2,j,k,m,iElem)/BGMdeltas(2))-REAL(GaussBGMIndex(2,j,k,m,iElem))
          GaussBGMFactor(3,j,k,m,iElem) = (Elem_xGP(3,j,k,m,iElem)/BGMdeltas(3))-REAL(GaussBGMIndex(3,j,k,m,iElem))
        END DO
      END DO
    END DO
  END DO
#ifdef MPI
  CALL MPIBackgroundMeshInit()
#else
  IF(GEO%nPeriodicVectors.GT.0)THEN
    ! Compute PeriodicBGMVectors (from PeriodicVectors and BGMdeltas)
    ALLOCATE(GEO%PeriodicBGMVectors(1:3,1:GEO%nPeriodicVectors),STAT=allocStat)
    IF (allocStat .NE. 0) THEN
      CALL abort(__STAMP__, &
      'ERROR in MPIBackgroundMeshInit: cannot allocate GEO%PeriodicBGMVectors!')
    END IF
    DO i = 1, GEO%nPeriodicVectors
      GEO%PeriodicBGMVectors(1,i) = NINT(GEO%PeriodicVectors(1,i)/BGMdeltas(1))
      IF(ABS(GEO%PeriodicVectors(1,i)/BGMdeltas(1)-REAL(GEO%PeriodicBGMVectors(1,i))).GT.1E-10)THEN
        CALL abort(__STAMP__, &
      'ERROR: Periodic Vector ist not multiple of background mesh delta')
      END IF
      GEO%PeriodicBGMVectors(2,i) = NINT(GEO%PeriodicVectors(2,i)/BGMdeltas(2))
      IF(ABS(GEO%PeriodicVectors(2,i)/BGMdeltas(2)-REAL(GEO%PeriodicBGMVectors(2,i))).GT.1E-10)THEN
        CALL abort(__STAMP__, &
      'ERROR: Periodic Vector ist not multiple of background mesh delta')
      END IF
      GEO%PeriodicBGMVectors(3,i) = NINT(GEO%PeriodicVectors(3,i)/BGMdeltas(3))
      IF(ABS(GEO%PeriodicVectors(3,i)/BGMdeltas(3)-REAL(GEO%PeriodicBGMVectors(3,i))).GT.1E-10)THEN
        CALL abort(__STAMP__, &
      'ERROR: Periodic Vector ist not multiple of background mesh delta')
      END IF
    END DO
  END IF
#endif
  ALLOCATE(BGMSource(BGMminX:BGMmaxX,BGMminY:BGMmaxY,BGMminZ:BGMmaxZ,1:4))
CASE('cartmesh_splines')
  BGMdeltas(1:3) = GETREALARRAY('PIC-BGMdeltas',3,'0. , 0. , 0.')
  FactorBGM(1:3) = GETREALARRAY('PIC-FactorBGM',3,'1. , 1. , 1.')
  BGMdeltas(1:3) = 1./FactorBGM(1:3)*BGMdeltas(1:3)
  IF (ANY(BGMdeltas.EQ.0)) THEN
    CALL abort(__STAMP__, &
      'ERROR: PIC-BGMdeltas: No size for the cartesian background mesh definded.')
  END IF
  ! calc min and max coordinates for local mesh
  xmin = 1.0E200
  ymin = 1.0E200
  zmin = 1.0E200
  xmax = -1.0E200
  ymax = -1.0E200
  zmax = -1.0E200
  DO iElem = 1, nElems
    xmin=MIN(xmin,MINVAL(XCL_NGeo(1,:,:,:,iElem)))
    xmax=MAX(xmax,MAXVAL(XCL_NGeo(1,:,:,:,iElem)))
    ymin=MIN(ymin,MINVAL(XCL_NGeo(2,:,:,:,iElem)))
    ymax=MAX(ymax,MAXVAL(XCL_NGeo(2,:,:,:,iElem)))
    zmin=MIN(zmin,MINVAL(XCL_NGeo(3,:,:,:,iElem)))
    zmax=MAX(zmax,MAXVAL(XCL_NGeo(3,:,:,:,iElem)))
  END DO
  ! define minimum and maximum backgroundmesh index, compute volume
  BGMVolume = BGMdeltas(1)*BGMdeltas(2)*BGMdeltas(3)
  BGMminX = FLOOR(xmin/BGMdeltas(1)-0.0001)-2
  BGMminY = FLOOR(ymin/BGMdeltas(2)-0.0001)-2
  BGMminZ = FLOOR(zmin/BGMdeltas(3)-0.0001)-2
  BGMmaxX = CEILING(xmax/BGMdeltas(1)+0.0001)+2
  BGMmaxY = CEILING(ymax/BGMdeltas(2)+0.0001)+2
  BGMmaxZ = CEILING(zmax/BGMdeltas(3)+0.0001)+2
  ! mapping from gausspoints to BGM
  ALLOCATE(GaussBGMIndex(1:3,0:PP_N,0:PP_N,0:PP_N,1:nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(__STAMP__, &
      'ERROR in pic_depo.f90: Cannot allocate GaussBGMIndex!')
  END IF
  ALLOCATE(GaussBGMFactor(1:3,0:PP_N,0:PP_N,0:PP_N,1:nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(__STAMP__, &
      'ERROR in pic_depo.f90: Cannot allocate GaussBGMFactor!')
  END IF
  DO iElem = 1, nElems
    DO j = 0, PP_N
      DO k = 0, PP_N
        DO m = 0, PP_N
          GaussBGMIndex(1,j,k,m,iElem) = FLOOR(Elem_xGP(1,j,k,m,iElem)/BGMdeltas(1))
          GaussBGMIndex(2,j,k,m,iElem) = FLOOR(Elem_xGP(2,j,k,m,iElem)/BGMdeltas(2))
          GaussBGMIndex(3,j,k,m,iElem) = FLOOR(Elem_xGP(3,j,k,m,iElem)/BGMdeltas(3))
        END DO
      END DO
    END DO
  END DO
  ! pre-build weights for BGM to gausspoint interpolation
  ALLOCATE(GPWeight(1:nElems,0:PP_N,0:PP_N,0:PP_N,1:4,1:4,1:4),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(__STAMP__, &
      'ERROR in pic_depo.f90: Cannot allocate GPWeight!')
  END IF
  DO iElem = 1, nElems
    DO j = 0, PP_N
      DO k = 0, PP_N
        DO m = 0, PP_N
          DO dir = 1,3               ! x,y,z direction
            DO weightrun = 0,3
              DO mm = 0, 3
                IF (mm.EQ.weightrun) then
                  auxiliary(mm) = 1.0
                ELSE
                  auxiliary(mm) = 0.0
                END IF
              END DO
              CALL DeBoor(GaussBGMIndex(dir,j,k,m,iElem),auxiliary,Elem_xGP(dir,j,k,m,iElem),weight(dir,weightrun),dir)
            END DO
          END DO
          DO r = 1,4
            DO s = 1,4
              DO t = 1,4
                GPWeight(iElem,j,k,m,r,s,t) = weight(1,-(r-4)) * weight(2,-(s-4)) * weight(3,-(t-4))
              END DO !t
            END DO !s
          END DO !r
        END DO !m
      END DO !k
    END DO !j
  END DO !iElem
#ifdef MPI
  CALL MPIBackgroundMeshInit()
#else
  IF(GEO%nPeriodicVectors.GT.0)THEN
    ! Compute PeriodicBGMVectors (from PeriodicVectors and BGMdeltas)
    ALLOCATE(GEO%PeriodicBGMVectors(1:3,1:GEO%nPeriodicVectors),STAT=allocStat)
    IF (allocStat .NE. 0) THEN
      CALL abort(__STAMP__, &
      'ERROR in MPIBackgroundMeshInit: cannot allocate GEO%PeriodicBGMVectors!')
    END IF
    DO i = 1, GEO%nPeriodicVectors
      GEO%PeriodicBGMVectors(1,i) = NINT(GEO%PeriodicVectors(1,i)/BGMdeltas(1))
      IF(ABS(GEO%PeriodicVectors(1,i)/BGMdeltas(1)-REAL(GEO%PeriodicBGMVectors(1,i))).GT.1E-10)THEN
        CALL abort(__STAMP__,  &
      'ERROR: Periodic Vector ist not multiple of background mesh delta')
      END IF
      GEO%PeriodicBGMVectors(2,i) = NINT(GEO%PeriodicVectors(2,i)/BGMdeltas(2))
      IF(ABS(GEO%PeriodicVectors(2,i)/BGMdeltas(2)-REAL(GEO%PeriodicBGMVectors(2,i))).GT.1E-10)THEN
        CALL abort(__STAMP__, &
      'ERROR: Periodic Vector ist not multiple of background mesh delta')
      END IF
      GEO%PeriodicBGMVectors(3,i) = NINT(GEO%PeriodicVectors(3,i)/BGMdeltas(3))
      IF(ABS(GEO%PeriodicVectors(3,i)/BGMdeltas(3)-REAL(GEO%PeriodicBGMVectors(3,i))).GT.1E-10)THEN
        CALL abort(__STAMP__, &
      'ERROR: Periodic Vector ist not multiple of background mesh delta')
      END IF
    END DO
  END IF
#endif
  ALLOCATE(BGMSource(BGMminX:BGMmaxX,BGMminY:BGMmaxY,BGMminZ:BGMmaxZ,1:4))
CASE DEFAULT
  CALL abort(__STAMP__, &
      'Unknown DepositionType in pic_depo.f90')
END SELECT

OutputSource = GETLOGICAL('PIC-OutputSource','F')

SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE DEPOSITION DONE!'

END SUBROUTINE InitializeDeposition


SUBROUTINE Deposition(doInnerParts)  
!============================================================================================================================
! This subroutine performes the deposition of the particle charge and current density to the grid
! following list of distribution methods are implemted
! - nearest blurrycenter (barycenter of hexahedra)
! - nearest Gauss Point  (only volome of IP - higher resolution than nearest blurrycenter )
! - shape function       (only one type implemented)
! - delta distributio
! useVMPF added, therefore, this routine contains automatically the use of variable mpfs
!============================================================================================================================
! use MODULES
USE MOD_PICDepo_Vars
USE MOD_Particle_Vars
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars,           ONLY:PI
USE MOD_Mesh_Vars,              ONLY:nElems, Elem_xGP, sJ
USE MOD_ChangeBasis,            ONLY:ChangeBasis3D
USE MOD_Interpolation_Vars,     ONLY:wGP,swGP
USE MOD_PICInterpolation_Vars,  ONLY:InterpolationType
USE MOD_Eval_xyz,               ONLY:eval_xyz_elemcheck
USE MOD_Basis,                  ONLY:LagrangeInterpolationPolys,BernSteinPolynomial
USE MOD_Interpolation_Vars,     ONLY:wBary,xGP
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping
USE MOD_Particle_Mesh_Vars,     ONLY:GEO,ElemBaryNGeo,casematrix, NbrOfCases
#ifdef MPI
! only required for shape function??
USE MOD_Particle_MPI_Vars,      ONLY:ExtPartState,ExtPartSpecies,ExtPartMPF,ExtPartToFIBGM,NbrOfExtParticles
USE MOD_Particle_MPI_Vars,      ONLY:PartMPI,PartMPIExchange
USE MOD_LoadBalance_Vars,       ONLY:nDeposPerElem,tCartMesh,ElemTime
#endif  /*MPI*/
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE                                                                                   
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT variable declaration                                                                       
LOGICAL                          :: doInnerParts
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT variable declaration                                                                       
!-----------------------------------------------------------------------------------------------------------------------------------
! Local variable declaration                                                                       
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                          :: firstPart,lastPart
INTEGER                          :: i,j, k, l, m, iElem, iPart
LOGICAL                          :: chargedone(1:nElems)!, SAVE_GAUSS             
LOGICAL                          :: SAVE_GAUSS
INTEGER                          :: BGMimax,BGMimin,BGMjmax,BGMjmin,BGMkmax,BGMkmin
INTEGER                          :: kmin, kmax, lmin, lmax, mmin, mmax                           
INTEGER                          :: kk, ll, mm, ppp                                              
INTEGER                          :: ElemID, iCase, ind, N_in
REAL                             :: radius2, S, S1, Fac(4)
!REAL                             :: GaussDistance(0:PP_N,0:PP_N,0:PP_N)
REAL                             :: Vec1(1:3), Vec2(1:3), Vec3(1:3), ShiftedPart(1:3)
INTEGER                          :: a,b, ii, expo
REAL                             :: ElemSource(nElems,1:4)
REAL                             :: Charge, TSource(1:4), auxiliary(0:3),weight(1:3,0:3), locweight
REAL                             :: alpha1, alpha2, alpha3
INTEGER                          :: PosInd(3),r,ss,t,u,v,w, dir, weightrun
INTEGER                          :: foundDOF
REAL,DIMENSION(3,0:NDepo)        :: L_xi
REAL                             :: DeltaIntCoeff,prefac
REAL                             :: local_r_sf, local_r2_sf, local_r2_sf_inv
#ifdef MPI
! load balance
REAL                             :: tLBStart,tLBEnd
#endif /*MPI*/
!============================================================================================================================

IF(doInnerParts)THEN
  source=0.0
  firstPart=1
  lastPart =PDM%ParticleVecLength
  !IF(firstPart.GT.lastPart) RETURN
ELSE
#ifdef MPI
  firstPart=PDM%ParticleVecLength-PartMPIExchange%nMPIParticles+1
  lastPart =PDM%ParticleVecLength
#else
  firstPart=1
  lastPart =0
#endif /*MPI*/
END IF
  
!IF((firstPart.GT.lastPart).AND.(DepositionType.NE.'delta_distri').AND.(DepositionType.NE.'shape_function')&
!                          .AND.(DepositionType.NE.'nearest_blurrycenter')) RETURN

SELECT CASE(TRIM(DepositionType))
CASE('nearest_blurrycenter')
  IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
  ElemSource=0.0
  DO iElem=1,PP_nElems
#ifdef MPI
    tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
    DO iPart=firstPart,lastPart
      IF(.NOT.PDM%ParticleInside(iPart))CYCLE
      IF(PEM%Element(iPart).EQ.iElem)THEN
        IF(usevMPF)THEN
!#if (PP_nVar==8)
         ElemSource(iElem,1:3) = ElemSource(iElem,1:3)+ &
                PartState(iPart,4:6)* Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)
!#endif
         ElemSource(iElem,4) = ElemSource(iElem,4) + & 
              Species(PartSpecies(iPart))%ChargeIC* PartMPF(iPart)
        ELSE
!#if (PP_nVar==8)
         ElemSource(iElem,1:3) = ElemSource(iElem,1:3)+ &
                PartState(iPart,4:6)* Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor
!#endif
         ElemSource(iElem,4) = ElemSource(iElem,4) + & 
              Species(PartSpecies(iPart))%ChargeIC* Species(PartSpecies(iPart))%MacroParticleFactor
        END IF ! usevMPF
      END IF ! Element(iPart).EQ.iElem
    END DO ! iPart
!#if (PP_nVar==8)
    source(1,:,:,:,iElem) = source(1,:,:,:,iElem)+ElemSource(iElem,1) 
    source(2,:,:,:,iElem) = source(2,:,:,:,iElem)+ElemSource(iElem,2) 
    source(3,:,:,:,iElem) = source(3,:,:,:,iElem)+ElemSource(iElem,3) 
!#endif                                            
    source(4,:,:,:,iElem) = source(4,:,:,:,iElem)+ElemSource(iElem,4) 
#ifdef MPI
    tLBEnd = LOCALTIME() ! LB Time End
    ElemTime(iElem)=ElemTime(iElem)+tLBEnd-tLBStart
#endif /*MPI*/
  END DO ! iElem=1,PP_nElems
  IF(.NOT.doInnerParts)THEN
    DO iElem=1,PP_nElems
#ifdef MPI
      tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
!#if (PP_nVar==8)
      source(1:4,:,:,:,iElem) = source(1:4,:,:,:,iElem) / GEO%Volume(iElem)
!#else
!      source(4,:,:,:,iElem) = source(4,:,:,:,iElem) / GEO%Volume(iElem)
!#endif
#ifdef MPI
      tLBEnd = LOCALTIME() ! LB Time End
      ElemTime(iElem)=ElemTime(iElem)+tLBEnd-tLBStart
#endif /*MPI*/
    END DO ! iElem=1,PP_nElems
  END IF ! .NOT. doInnerParts
CASE('shape_function')
  IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
  Vec1(1:3) = 0.
  Vec2(1:3) = 0.
  Vec3(1:3) = 0.
  IF (GEO%nPeriodicVectors.EQ.1) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
  END IF
  IF (GEO%nPeriodicVectors.EQ.2) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
    Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
  END IF
  IF (GEO%nPeriodicVectors.EQ.3) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
    Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
    Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
  END IF
  DO iPart=firstPart,LastPart
    IF (PDM%ParticleInside(iPart)) THEN
      IF (usevMPF) THEN
        Fac(4)= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)*w_sf
      ELSE
        Fac(4)= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor*w_sf
      END IF ! usevMPF
      !IF(fac(4).GT.0.) print*,'charge pos'
      Fac(1:3) = PartState(iPart,4:6)*Fac(4)
      !-- determine which background mesh cells (and interpolation points within) need to be considered
      DO iCase = 1, NbrOfCases
        chargedone(:) = .FALSE.
        DO ind = 1,3
          ShiftedPart(ind) = PartState(iPart,ind) + casematrix(iCase,1)*Vec1(ind) + &
               casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
        END DO
        kmax = CEILING((ShiftedPart(1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
        kmax = MIN(kmax,GEO%FIBGMimax)
        kmin = FLOOR((ShiftedPart(1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
        kmin = MAX(kmin,GEO%FIBGMimin)
        lmax = CEILING((ShiftedPart(2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
        lmax = MIN(lmax,GEO%FIBGMjmax)
        lmin = FLOOR((ShiftedPart(2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
        lmin = MAX(lmin,GEO%FIBGMjmin)
        mmax = CEILING((ShiftedPart(3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
        mmax = MIN(mmax,GEO%FIBGMkmax)
        mmin = FLOOR((ShiftedPart(3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
        mmin = MAX(mmin,GEO%FIBGMkmin)
        !-- go through all these cells
        DO kk = kmin,kmax
          DO ll = lmin, lmax
            DO mm = mmin, mmax
              !--- go through all mapped elements not done yet
              DO ppp = 1,GEO%FIBGM(kk,ll,mm)%nElem
                ElemID = GEO%FIBGM(kk,ll,mm)%Element(ppp)
                IF(ElemID.GT.nElems) CYCLE
                IF (.NOT.chargedone(ElemID)) THEN
#ifdef MPI
                  nDeposPerElem(ElemID)=nDeposPerElem(ElemID)+1
#endif /*MPI*/
                  !--- go through all gauss points
                  !CALL ComputeGaussDistance(PP_N,r2_sf_inv,ShiftedPart,ElemDepo_xGP(:,:,:,:,ElemID),GaussDistance)
                  DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
                    !-- calculate distance between gauss and particle
                    radius2 = (ShiftedPart(1) - ElemDepo_xGP(1,k,l,m,ElemID)) * (ShiftedPart(1) - ElemDepo_xGP(1,k,l,m,ElemID)) &
                            + (ShiftedPart(2) - ElemDepo_xGP(2,k,l,m,ElemID)) * (ShiftedPart(2) - ElemDepo_xGP(2,k,l,m,ElemID)) &
                            + (ShiftedPart(3) - ElemDepo_xGP(3,k,l,m,ElemID)) * (ShiftedPart(3) - ElemDepo_xGP(3,k,l,m,ElemID))
                    !-- calculate charge and current density at ip point using a shape function
                    !-- currently only one shapefunction available, more to follow (including structure change)
                    IF (radius2 .LT. r2_sf) THEN
                      S = 1 - r2_sf_inv * radius2
                    !radius2=GaussDistance(k,l,m)
                    !IF (radius2 .LT. 1.0) THEN
                    !  S = 1 -  radius2
                      S1 = S*S
                      DO expo = 3, alpha_sf
                        S1 = S*S1
                      END DO
                      source(1:3,k,l,m,ElemID) = source(1:3,k,l,m,ElemID) + Fac(1:3) * S1
                      source( 4 ,k,l,m,ElemID) = source( 4 ,k,l,m,ElemID) + Fac(4) * S1
                    END IF
                  END DO; END DO; END DO
                  chargedone(ElemID) = .TRUE.
                END IF
              END DO ! ppp
            END DO ! mm
          END DO ! ll
        END DO ! kk
      END DO ! iCase (periodicity)
    END IF ! inside
  END DO ! i
#ifdef MPI           
  IF(.NOT.DoInnerParts)THEN
    Vec1(1:3) = 0.
    Vec2(1:3) = 0.
    Vec3(1:3) = 0.
    IF (NbrOfextParticles .GT. 0) THEN
      IF (GEO%nPeriodicVectors.EQ.1) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
      END IF
      IF (GEO%nPeriodicVectors.EQ.2) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
        Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
      END IF
      IF (GEO%nPeriodicVectors.EQ.3) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
        Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
        Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
      END IF
    END IF
    
    DO iPart=1,NbrOfextParticles  !external Particles
      IF (usevMPF) THEN
        Fac(4)= Species(ExtPartSpecies(iPart))%ChargeIC * ExtPartMPF(iPart)*w_sf
      ELSE
        Fac(4)= Species(ExtPartSpecies(iPart))%ChargeIC * Species(ExtPartSpecies(iPart))%MacroParticleFactor*w_sf
      END IF ! usevMPF
      Fac(1:3) = ExtPartState(iPart,4:6)*Fac(4)
      !-- determine which background mesh cells (and interpolation points within) need to be considered
      DO iCase = 1, NbrOfCases
        chargedone(:) = .FALSE.
        DO ind = 1,3
          ShiftedPart(ind) = ExtPartState(iPart,ind) + casematrix(iCase,1)*Vec1(ind) + &
               casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
        END DO
        kmax = CEILING((ShiftedPart(1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
        kmax = MIN(kmax,GEO%FIBGMimax)
        kmin = FLOOR((ShiftedPart(1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
        kmin = MAX(kmin,GEO%FIBGMimin)
        lmax = CEILING((ShiftedPart(2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
        lmax = MIN(lmax,GEO%FIBGMjmax)
        lmin = FLOOR((ShiftedPart(2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
        lmin = MAX(lmin,GEO%FIBGMjmin)
        mmax = CEILING((ShiftedPart(3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
        mmax = MIN(mmax,GEO%FIBGMkmax)
        mmin = FLOOR((ShiftedPart(3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
        mmin = MAX(mmin,GEO%FIBGMkmin)
        !-- go through all these cells (should go through non if periodic and shiftedpart not in my domain
        DO kk = kmin,kmax
          DO ll = lmin, lmax
            DO mm = mmin, mmax
              !--- go through all mapped elements not done yet
              DO ppp = 1,GEO%FIBGM(kk,ll,mm)%nElem
                ElemID = GEO%FIBGM(kk,ll,mm)%Element(ppp)
                IF(ElemID.GT.nElems) CYCLE
                IF (.NOT.chargedone(ElemID)) THEN
#ifdef MPI
                  nDeposPerElem(ElemID)=nDeposPerElem(ElemID)+1
#endif /*MPI*/
                  !--- go through all gauss points
                  !CALL ComputeGaussDistance(PP_N,r2_sf_inv,ShiftedPart,ElemDepo_xGP(:,:,:,:,ElemID),GaussDistance)
                  DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
                    !-- calculate distance between gauss and particle
                      radius2 = (ShiftedPart(1) - ElemDepo_xGP(1,k,l,m,ElemID)) * (ShiftedPart(1) - ElemDepo_xGP(1,k,l,m,ElemID)) &
                              + (ShiftedPart(2) - ElemDepo_xGP(2,k,l,m,ElemID)) * (ShiftedPart(2) - ElemDepo_xGP(2,k,l,m,ElemID)) &
                              + (ShiftedPart(3) - ElemDepo_xGP(3,k,l,m,ElemID)) * (ShiftedPart(3) - ElemDepo_xGP(3,k,l,m,ElemID))
                      !-- calculate charge and current density at ip point using a shape function
                      !-- currently only one shapefunction available, more to follow (including structure change)
                      IF (radius2 .LT. r2_sf) THEN
                        S = 1. - r2_sf_inv * radius2
                        !IF(S.LT.0.) print*,'dist neg '
                      !radius2=GaussDistance(k,l,m)
                      !IF (radius2 .LT. 1.0) THEN
                      !  S = 1 -  radius2
                        S1 = S*S
                        DO expo = 3, alpha_sf
                          S1 = S*S1
                        END DO
                        source(1:3,k,l,m,ElemID) = source(1:3,k,l,m,ElemID) + Fac(1:3) * S1
                        source( 4 ,k,l,m,ElemID) = source( 4 ,k,l,m,ElemID) + Fac(4) * S1
                      END IF
                  END DO; END DO; END DO
                  chargedone(ElemID) = .TRUE.
                END IF
              END DO ! ppp
            END DO ! mm
          END DO ! ll
        END DO ! kk
      END DO
    END DO
    ! deallocate external state
    SDEALLOCATE(ExtPartState)
    SDEALLOCATE(ExtPartSpecies)
    SDEALLOCATE(ExtPartToFIBGM)
    SDEALLOCATE(ExtPartMPF)
    NbrOfExtParticles=0
  END IF
#endif /*MPI*/
! my work for higher N
!    N_in=PP_N !(PP_N+1)*(PP_N+1)*(PP_N+1)
!    !foundDOF=0
!    ! initialize periodic vectors
!    SELECT CASE(GEO%nPeriodicVectors)
!    CASE(1)
!      Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
!      Vec2(1:3) = 0.
!      Vec3(1:3) = 0.
!    CASE(2)
!      Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
!      Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
!      Vec3(1:3) = 0.
!    CASE(3)
!      Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
!      Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
!      Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
!    CASE DEFAULT
!      Vec1(1:3) = 0.
!      Vec2(1:3) = 0.
!      Vec3(1:3) = 0.
!    END SELECT
!    ! loop over all cases
!    DO iCase = 1, NbrOfCases
!      DO iPart=firstPart,LastPart
!        IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
!        DO ind = 1,3
!          ShiftedPart(ind) = PartState(iPart,ind) + casematrix(iCase,1)*Vec1(ind) + &
!               casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
!        END DO
!        BGMimax = CEILING((ShiftedPart(1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
!        BGMimax = MIN(BGMimax,GEO%FIBGMimax)
!        PartToFIBGM(2,iPart)=BGMimax
!        BGMimin = FLOOR((ShiftedPart(1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
!        BGMimin = MAX(BGMimin,GEO%FIBGMimin)
!        PartToFIBGM(1,iPart)=BGMimin
!        BGMjmax = CEILING((ShiftedPart(2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
!        BGMjmax = MIN(BGMjmax,GEO%FIBGMjmax)
!        PartToFIBGM(4,iPart)=BGMjmax
!        BGMjmin = FLOOR((ShiftedPart(2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
!        BGMjmin = MAX(BGMjmin,GEO%FIBGMjmin)
!        PartToFIBGM(3,iPart)=BGMjmin
!        BGMkmax = CEILING((ShiftedPart(3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
!        BGMkmax = MIN(BGMkmax,GEO%FIBGMkmax)
!        PartToFIBGM(6,iPart)=BGMkmax
!        BGMkmin = FLOOR((ShiftedPart(3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
!        BGMkmin = MAX(BGMkmin,GEO%FIBGMkmin)
!        PartToFIBGM(5,iPart)=BGMkmin
!      END DO ! iPart=firstPart,LastPart
!  #ifdef MPI
!      IF(.NOT.DoInnerParts)THEN
!        ! stuff to add, ExtPartToFIBGM
!        DO iPart=1,NbrOfextParticles  !external Particles
!          DO ind = 1,3
!            ShiftedPart(ind) = ExtPartState(iPart,ind) + casematrix(iCase,1)*Vec1(ind) + &
!                 casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
!          END DO
!          BGMimax = CEILING((ShiftedPart(1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
!          BGMimax = MIN(BGMimax,GEO%FIBGMimax)
!          ExtPartToFIBGM(2,iPart)=BGMimax
!          BGMimin = FLOOR((ShiftedPart(1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
!          BGMimin = MAX(BGMimin,GEO%FIBGMimin)
!          ExtPartToFIBGM(1,iPart)=BGMimin
!          BGMjmax = CEILING((ShiftedPart(2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
!          BGMjmax = MIN(BGMjmax,GEO%FIBGMjmax)
!          ExtPartToFIBGM(4,iPart)=BGMjmax
!          BGMjmin = FLOOR((ShiftedPart(2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
!          BGMjmin = MAX(BGMjmin,GEO%FIBGMjmin)
!          ExtPartToFIBGM(3,iPart)=BGMjmin
!          BGMkmax = CEILING((ShiftedPart(3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
!          BGMkmax = MIN(BGMkmax,GEO%FIBGMkmax)
!          ExtPartToFIBGM(6,iPart)=BGMkmax
!          BGMkmin = FLOOR((ShiftedPart(3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
!          BGMkmin = MAX(BGMkmin,GEO%FIBGMkmin)
!          ExtPartToFIBGM(5,iPart)=BGMkmin
!        END DO ! iPart=1,NbrOfextParticles  !external Particles
!      END IF
!  #endif /*MPI*/
!      ! loop over all elements
!      DO iElem=1,PP_nElems
!        DO iPart=firstPart,LastPart
!          IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
!          ! recomputed
!          DO ind = 1,3
!            ShiftedPart(ind) = PartState(iPart,ind) + casematrix(iCase,1)*Vec1(ind) + &
!                 casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
!          END DO
!          ! check if bounding FIBGM box of particle and element overlap
!          IF(PartToFIBGM(1,iPart).GT.GEO%ElemToFIBGM(2,iElem)) CYCLE
!          IF(PartToFIBGM(2,iPart).LT.GEO%ElemToFIBGM(1,iElem)) CYCLE
!          IF(PartToFIBGM(3,iPart).GT.GEO%ElemToFIBGM(4,iElem)) CYCLE
!          IF(PartToFIBGM(4,iPart).LT.GEO%ElemToFIBGM(3,iElem)) CYCLE
!          IF(PartToFIBGM(5,iPart).GT.GEO%ElemToFIBGM(6,iElem)) CYCLE
!          IF(PartToFIBGM(6,iPart).LT.GEO%ElemToFIBGM(5,iElem)) CYCLE
!          ! check if element circumcircle and shape function can interact
!          radius2= (ShiftedPart(1) - ElemBaryNGeo(1,iElem)) * (ShiftedPart(1) - ElemBaryNGeo(1,iElem) ) &
!                 + (ShiftedPart(2) - ElemBaryNGeo(2,iElem)) * (ShiftedPart(2) - ElemBaryNGeo(2,iElem) ) &
!                 + (ShiftedPart(3) - ElemBaryNGeo(3,iElem)) * (ShiftedPart(3) - ElemBaryNGeo(3,iElem) )
!          IF(radius2.GT.ElemRadius2_sf(iElem)) CYCLE
!          IF (usevMPF) THEN
!            Fac(4)= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)*w_sf
!          ELSE
!            Fac(4)= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor*w_sf
!          END IF ! usevMPF
!          Fac(1:3) = PartState(iPart,4:6)*Fac(4)
!          ! now, compute the actual 
!          CALL ComputeGaussDistance(N_in,r2_sf_inv,ShiftedPart,ElemDepo_xGP(:,:,:,:,iElem),GaussDistance)
!          DO k=0,PP_N
!            DO j=0,PP_N
!              DO i=0,PP_N
!                radius2=GaussDistance(i,j,k)
!                IF(radius2.GT.1.0) CYCLE
!  !              foundDOF=foundDOF+1
!                S = 1. - radius2
!                S1 = S*S
!                DO expo = 3, alpha_sf,1
!                  S1 = S*S1
!                END DO
!                source(1:3,i,j,k,iElem) = source(1:3,i,j,k,iElem) + Fac(1:3) * S1
!                source( 4 ,i,j,k,iElem) = source( 4 ,i,j,k,iElem) + Fac(4) * S1
!              END DO ! i=0,PP_N
!            END DO ! j=0,PP_N
!          END DO ! k=0,PP_N
!        END DO ! iPart=firstPart,LastPart
!        ! and for second step the external particles
!  #ifdef MPI
!        IF(.NOT.DoInnerParts)THEN
!          ! stuff to add, ExtPartToFIBGM
!          DO iPart=1,NbrOfextParticles  !external Particles
!            DO ind = 1,3
!              ShiftedPart(ind) = ExtPartState(iPart,ind) + casematrix(iCase,1)*Vec1(ind) + &
!                   casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
!            END DO
!            ! check if bounding FIBGM box of particle and element overlap
!            IF(ExtPartToFIBGM(1,iPart).GT.GEO%ElemToFIBGM(2,iElem)) CYCLE
!            IF(ExtPartToFIBGM(2,iPart).LT.GEO%ElemToFIBGM(1,iElem)) CYCLE
!            IF(ExtPartToFIBGM(3,iPart).GT.GEO%ElemToFIBGM(4,iElem)) CYCLE
!            IF(ExtPartToFIBGM(4,iPart).LT.GEO%ElemToFIBGM(3,iElem)) CYCLE
!            IF(ExtPartToFIBGM(5,iPart).GT.GEO%ElemToFIBGM(6,iElem)) CYCLE
!            IF(ExtPartToFIBGM(6,iPart).LT.GEO%ElemToFIBGM(5,iElem)) CYCLE
!            ! check if element circumcircle and shape function can interact
!            radius2= (ShiftedPart(1) - ElemBaryNGeo(1,iElem))* (ShiftedPart(1) - ElemBaryNGeo(1,iElem))&
!                   + (ShiftedPart(2) - ElemBaryNGeo(2,iElem))* (ShiftedPart(2) - ElemBaryNGeo(2,iElem))&
!                   + (ShiftedPart(3) - ElemBaryNGeo(3,iElem))* (ShiftedPart(3) - ElemBaryNGeo(3,iElem))
!            IF(radius2.GT.ElemRadius2_sf(iElem)) CYCLE
!            IF (usevMPF) THEN
!              Fac(4)= Species(ExtPartSpecies(iPart))%ChargeIC * ExtPartMPF(iPart)*w_sf
!            ELSE
!              Fac(4)= Species(ExtPartSpecies(iPart))%ChargeIC * Species(ExtPartSpecies(iPart))%MacroParticleFactor*w_sf
!            END IF ! usevMPF
!            Fac(1:3) = ExtPartState(iPart,4:6)*Fac(4)
!            ! now, compute the actual 
!            CALL ComputeGaussDistance(N_in,r2_sf_inv,ShiftedPart,ElemDepo_xGP(:,:,:,:,iElem),GaussDistance)
!            DO k=0,PP_N
!              DO j=0,PP_N
!                DO i=0,PP_N
!                  radius2=GaussDistance(i,j,k)
!                  IF(radius2.GT.1.0) CYCLE
!                  S = 1.0 - radius2
!                  S1 = S*S
!                  DO expo = 3, alpha_sf,1
!                    S1 = S*S1
!                  END DO
!    !              foundDOF=foundDOF+1
!                  source(1:3,i,j,k,iElem) = source(1:3,i,j,k,iElem) + Fac(1:3) * S1
!                  source( 4 ,i,j,k,iElem) = source( 4 ,i,j,k,iElem) + Fac(4) * S1
!                END DO ! i=0,PP_N
!              END DO ! j=0,PP_N
!            END DO ! k=0,PP_N
!          END DO ! iPart=firstPart,LastPart
!        END IF
!  #endif /*MPI*/
!      END DO ! iElem=1,PP_nElems
!    END DO ! iCase
!  #ifdef MPI
!    IF(.NOT.DoInnerParts)THEN
!      ! deallocate external state
!      SDEALLOCATE(ExtPartState)
!      SDEALLOCATE(ExtPartSpecies)
!      SDEALLOCATE(ExtPartToFIBGM)
!      SDEALLOCATE(ExtPartMPF)
!    END IF
!  #endif /*MPI*/
!    !IPWRITE(*,*) 'Number of found DOFs',FoundDOF

  IF( .NOT.DoInnerParts .AND. DoSFEqui) THEN
    ! map source from Equististant points on Gauss-Points
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(4,PP_N,PP_N,Vdm_EquiN_GaussN,source(:,:,:,:,iElem),source(:,:,:,:,iElem))
    END DO ! iElem=1,PP_nElems
  END IF

CASE('shape_function_1d')
  IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
  Vec1(1:3) = 0.
  Vec2(1:3) = 0.
  Vec3(1:3) = 0.
  IF (GEO%nPeriodicVectors.EQ.1) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
  END IF
  IF (GEO%nPeriodicVectors.EQ.2) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
    Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
  END IF
  IF (GEO%nPeriodicVectors.EQ.3) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
    Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
    Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
  END IF
  DO iPart=firstPart,LastPart
    IF (PDM%ParticleInside(iPart)) THEN
      IF (usevMPF) THEN
        Fac(4)= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)*w_sf
      ELSE
        Fac(4)= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor*w_sf
      END IF ! usevMPF
      !IF(fac(4).GT.0.) print*,'charge pos'
      Fac(1:3) = PartState(iPart,4:6)*Fac(4)
      !-- determine which background mesh cells (and interpolation points within) need to be considered
      chargedone(:) = .FALSE.
      DO iCase = 1, NbrOfCases
        DO ind = 1,3
          ShiftedPart(ind) = PartState(iPart,ind) + casematrix(iCase,1)*Vec1(ind) + &
               casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
        END DO
        IF(sf1d_dir.EQ.1)THEN
          kmax = CEILING((ShiftedPart(1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
          kmax = MIN(kmax,GEO%FIBGMimax)
          kmin = FLOOR((ShiftedPart(1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
          kmin = MAX(kmin,GEO%FIBGMimin)
          lmax = GEO%FIBGMjmax
          lmin = GEO%FIBGMjmin
          mmax = GEO%FIBGMkmax
          mmin = GEO%FIBGMkmin
        ELSEIF(sf1d_dir.EQ.2)THEN
          kmax = GEO%FIBGMimax
          kmin = GEO%FIBGMimin
          lmax = CEILING((ShiftedPart(2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
          lmax = MIN(lmax,GEO%FIBGMjmax)
          lmin = FLOOR((ShiftedPart(2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
          lmin = MAX(lmin,GEO%FIBGMjmin)
          mmax = GEO%FIBGMkmax
          mmin = GEO%FIBGMkmin
        ELSE
          kmax = GEO%FIBGMimax
          kmin = GEO%FIBGMimin
          lmax = GEO%FIBGMjmax
          lmin = GEO%FIBGMjmin
          mmax = CEILING((ShiftedPart(3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
          mmax = MIN(mmax,GEO%FIBGMkmax)
          mmin = FLOOR((ShiftedPart(3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
          mmin = MAX(mmin,GEO%FIBGMkmin)
        END IF
        !-- go through all these cells
        DO kk = kmin,kmax
          DO ll = lmin, lmax
            DO mm = mmin, mmax
              !--- go through all mapped elements not done yet
              DO ppp = 1,GEO%FIBGM(kk,ll,mm)%nElem
                ElemID = GEO%FIBGM(kk,ll,mm)%Element(ppp)
                IF(ElemID.GT.nElems) CYCLE
                IF (.NOT.chargedone(ElemID)) THEN
#ifdef MPI
                  nDeposPerElem(ElemID)=nDeposPerElem(ElemID)+1
#endif /*MPI*/
                  !--- go through all gauss points
                  !CALL ComputeGaussDistance(PP_N,r2_sf_inv,ShiftedPart,ElemDepo_xGP(:,:,:,:,ElemID),GaussDistance)
                  DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
                    !-- calculate distance between gauss and particle
                    radius2 = (ShiftedPart(sf1d_dir) - ElemDepo_xGP(sf1d_dir,k,l,m,ElemID)) &
                            * (ShiftedPart(sf1d_dir) - ElemDepo_xGP(sf1d_dir,k,l,m,ElemID))
                    !-- calculate charge and current density at ip point using a shape function
                    !-- currently only one shapefunction available, more to follow (including structure change)
                    IF (radius2 .LT. r2_sf) THEN
                      S = 1 - r2_sf_inv * radius2
                    !radius2=GaussDistance(k,l,m)
                    !IF (radius2 .LT. 1.0) THEN
                    !  S = 1 -  radius2
                      S1 = S*S
                      DO expo = 3, alpha_sf
                        S1 = S*S1
                      END DO
                      source(1:3,k,l,m,ElemID) = source(1:3,k,l,m,ElemID) + Fac(1:3) * S1
                      source( 4 ,k,l,m,ElemID) = source( 4 ,k,l,m,ElemID) + Fac(4) * S1
                    END IF
                  END DO; END DO; END DO
                  chargedone(ElemID) = .TRUE.
                END IF
              END DO ! ppp
            END DO ! mm
          END DO ! ll
        END DO ! kk
      END DO ! iCase (periodicity)
    END IF ! inside
  END DO ! i
#ifdef MPI           
  IF(.NOT.DoInnerParts)THEN
    Vec1(1:3) = 0.
    Vec2(1:3) = 0.
    Vec3(1:3) = 0.
    IF (NbrOfextParticles .GT. 0) THEN
      IF (GEO%nPeriodicVectors.EQ.1) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
      END IF
      IF (GEO%nPeriodicVectors.EQ.2) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
        Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
      END IF
      IF (GEO%nPeriodicVectors.EQ.3) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
        Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
        Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
      END IF
    END IF
    
    DO iPart=1,NbrOfextParticles  !external Particles
      IF (usevMPF) THEN
        Fac(4)= Species(ExtPartSpecies(iPart))%ChargeIC * ExtPartMPF(iPart)*w_sf
      ELSE
        Fac(4)= Species(ExtPartSpecies(iPart))%ChargeIC * Species(ExtPartSpecies(iPart))%MacroParticleFactor*w_sf
      END IF ! usevMPF
      Fac(1:3) = ExtPartState(iPart,4:6)*Fac(4)
      chargedone(:) = .FALSE.
      !-- determine which background mesh cells (and interpolation points within) need to be considered
      DO iCase = 1, NbrOfCases
        DO ind = 1,3
          ShiftedPart(ind) = ExtPartState(iPart,ind) + casematrix(iCase,1)*Vec1(ind) + &
               casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
        END DO
        IF(sf1d_dir.EQ.1)THEN
          kmax = CEILING((ShiftedPart(1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
          kmax = MIN(kmax,GEO%FIBGMimax)
          kmin = FLOOR((ShiftedPart(1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
          kmin = MAX(kmin,GEO%FIBGMimin)
          lmax = GEO%FIBGMjmax
          lmin = GEO%FIBGMjmin
          mmax = GEO%FIBGMkmax
          mmin = GEO%FIBGMkmin
        ELSEIF(sf1d_dir.EQ.2)THEN
          kmax = GEO%FIBGMimax
          kmin = GEO%FIBGMimin
          lmax = CEILING((ShiftedPart(2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
          lmax = MIN(lmax,GEO%FIBGMjmax)
          lmin = FLOOR((ShiftedPart(2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
          lmin = MAX(lmin,GEO%FIBGMjmin)
          mmax = GEO%FIBGMkmax
          mmin = GEO%FIBGMkmin
        ELSE
          kmax = GEO%FIBGMimax
          kmin = GEO%FIBGMimin
          lmax = GEO%FIBGMjmax
          lmin = GEO%FIBGMjmin
          mmax = CEILING((ShiftedPart(3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
          mmax = MIN(mmax,GEO%FIBGMkmax)
          mmin = FLOOR((ShiftedPart(3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
          mmin = MAX(mmin,GEO%FIBGMkmin)
        END IF
        !-- go through all these cells (should go through non if periodic and shiftedpart not in my domain
        DO kk = kmin,kmax
          DO ll = lmin, lmax
            DO mm = mmin, mmax
              !--- go through all mapped elements not done yet
              DO ppp = 1,GEO%FIBGM(kk,ll,mm)%nElem
                ElemID = GEO%FIBGM(kk,ll,mm)%Element(ppp)
                IF(ElemID.GT.nElems) CYCLE
                IF (.NOT.chargedone(ElemID)) THEN
#ifdef MPI
                  nDeposPerElem(ElemID)=nDeposPerElem(ElemID)+1
#endif /*MPI*/
                  !--- go through all gauss points
                  !CALL ComputeGaussDistance(PP_N,r2_sf_inv,ShiftedPart,ElemDepo_xGP(:,:,:,:,ElemID),GaussDistance)
                  DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
                    !-- calculate distance between gauss and particle
                      radius2 = (ShiftedPart(sf1d_dir) - ElemDepo_xGP(sf1d_dir,k,l,m,ElemID)) &
                              * (ShiftedPart(sf1d_dir) - ElemDepo_xGP(sf1d_dir,k,l,m,ElemID))
                      !-- calculate charge and current density at ip point using a shape function
                      !-- currently only one shapefunction available, more to follow (including structure change)
                      IF (radius2 .LT. r2_sf) THEN
                        S = 1. - r2_sf_inv * radius2
                        !IF(S.LT.0.) print*,'dist neg '
                      !radius2=GaussDistance(k,l,m)
                      !IF (radius2 .LT. 1.0) THEN
                      !  S = 1 -  radius2
                        S1 = S*S
                        DO expo = 3, alpha_sf
                          S1 = S*S1
                        END DO
                        source(1:3,k,l,m,ElemID) = source(1:3,k,l,m,ElemID) + Fac(1:3) * S1
                        source( 4 ,k,l,m,ElemID) = source( 4 ,k,l,m,ElemID) + Fac(4) * S1
                      END IF
                  END DO; END DO; END DO
                  chargedone(ElemID) = .TRUE.
                END IF
              END DO ! ppp
            END DO ! mm
          END DO ! ll
        END DO ! kk
      END DO
    END DO
    ! deallocate external state
    SDEALLOCATE(ExtPartState)
    SDEALLOCATE(ExtPartSpecies)
    SDEALLOCATE(ExtPartToFIBGM)
    SDEALLOCATE(ExtPartMPF)
    NbrOfExtParticles=0
  END IF
#endif /*MPI*/
! my work for higher N
!    N_in=PP_N !(PP_N+1)*(PP_N+1)*(PP_N+1)
!    !foundDOF=0
!    ! initialize periodic vectors
!    SELECT CASE(GEO%nPeriodicVectors)
!    CASE(1)
!      Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
!      Vec2(1:3) = 0.
!      Vec3(1:3) = 0.
!    CASE(2)
!      Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
!      Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
!      Vec3(1:3) = 0.
!    CASE(3)
!      Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
!      Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
!      Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
!    CASE DEFAULT
!      Vec1(1:3) = 0.
!      Vec2(1:3) = 0.
!      Vec3(1:3) = 0.
!    END SELECT
!    ! loop over all cases
!    DO iCase = 1, NbrOfCases
!      DO iPart=firstPart,LastPart
!        IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
!        DO ind = 1,3
!          ShiftedPart(ind) = PartState(iPart,ind) + casematrix(iCase,1)*Vec1(ind) + &
!               casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
!        END DO
!        BGMimax = CEILING((ShiftedPart(1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
!        BGMimax = MIN(BGMimax,GEO%FIBGMimax)
!        PartToFIBGM(2,iPart)=BGMimax
!        BGMimin = FLOOR((ShiftedPart(1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
!        BGMimin = MAX(BGMimin,GEO%FIBGMimin)
!        PartToFIBGM(1,iPart)=BGMimin
!        BGMjmax = CEILING((ShiftedPart(2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
!        BGMjmax = MIN(BGMjmax,GEO%FIBGMjmax)
!        PartToFIBGM(4,iPart)=BGMjmax
!        BGMjmin = FLOOR((ShiftedPart(2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
!        BGMjmin = MAX(BGMjmin,GEO%FIBGMjmin)
!        PartToFIBGM(3,iPart)=BGMjmin
!        BGMkmax = CEILING((ShiftedPart(3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
!        BGMkmax = MIN(BGMkmax,GEO%FIBGMkmax)
!        PartToFIBGM(6,iPart)=BGMkmax
!        BGMkmin = FLOOR((ShiftedPart(3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
!        BGMkmin = MAX(BGMkmin,GEO%FIBGMkmin)
!        PartToFIBGM(5,iPart)=BGMkmin
!      END DO ! iPart=firstPart,LastPart
!  #ifdef MPI
!      IF(.NOT.DoInnerParts)THEN
!        ! stuff to add, ExtPartToFIBGM
!        DO iPart=1,NbrOfextParticles  !external Particles
!          DO ind = 1,3
!            ShiftedPart(ind) = ExtPartState(iPart,ind) + casematrix(iCase,1)*Vec1(ind) + &
!                 casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
!          END DO
!          BGMimax = CEILING((ShiftedPart(1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
!          BGMimax = MIN(BGMimax,GEO%FIBGMimax)
!          ExtPartToFIBGM(2,iPart)=BGMimax
!          BGMimin = FLOOR((ShiftedPart(1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
!          BGMimin = MAX(BGMimin,GEO%FIBGMimin)
!          ExtPartToFIBGM(1,iPart)=BGMimin
!          BGMjmax = CEILING((ShiftedPart(2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
!          BGMjmax = MIN(BGMjmax,GEO%FIBGMjmax)
!          ExtPartToFIBGM(4,iPart)=BGMjmax
!          BGMjmin = FLOOR((ShiftedPart(2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
!          BGMjmin = MAX(BGMjmin,GEO%FIBGMjmin)
!          ExtPartToFIBGM(3,iPart)=BGMjmin
!          BGMkmax = CEILING((ShiftedPart(3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
!          BGMkmax = MIN(BGMkmax,GEO%FIBGMkmax)
!          ExtPartToFIBGM(6,iPart)=BGMkmax
!          BGMkmin = FLOOR((ShiftedPart(3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
!          BGMkmin = MAX(BGMkmin,GEO%FIBGMkmin)
!          ExtPartToFIBGM(5,iPart)=BGMkmin
!        END DO ! iPart=1,NbrOfextParticles  !external Particles
!      END IF
!  #endif /*MPI*/
!      ! loop over all elements
!      DO iElem=1,PP_nElems
!        DO iPart=firstPart,LastPart
!          IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
!          ! recomputed
!          DO ind = 1,3
!            ShiftedPart(ind) = PartState(iPart,ind) + casematrix(iCase,1)*Vec1(ind) + &
!                 casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
!          END DO
!          ! check if bounding FIBGM box of particle and element overlap
!          IF(PartToFIBGM(1,iPart).GT.GEO%ElemToFIBGM(2,iElem)) CYCLE
!          IF(PartToFIBGM(2,iPart).LT.GEO%ElemToFIBGM(1,iElem)) CYCLE
!          IF(PartToFIBGM(3,iPart).GT.GEO%ElemToFIBGM(4,iElem)) CYCLE
!          IF(PartToFIBGM(4,iPart).LT.GEO%ElemToFIBGM(3,iElem)) CYCLE
!          IF(PartToFIBGM(5,iPart).GT.GEO%ElemToFIBGM(6,iElem)) CYCLE
!          IF(PartToFIBGM(6,iPart).LT.GEO%ElemToFIBGM(5,iElem)) CYCLE
!          ! check if element circumcircle and shape function can interact
!          radius2= (ShiftedPart(1) - ElemBaryNGeo(1,iElem)) * (ShiftedPart(1) - ElemBaryNGeo(1,iElem) ) &
!                 + (ShiftedPart(2) - ElemBaryNGeo(2,iElem)) * (ShiftedPart(2) - ElemBaryNGeo(2,iElem) ) &
!                 + (ShiftedPart(3) - ElemBaryNGeo(3,iElem)) * (ShiftedPart(3) - ElemBaryNGeo(3,iElem) )
!          IF(radius2.GT.ElemRadius2_sf(iElem)) CYCLE
!          IF (usevMPF) THEN
!            Fac(4)= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)*w_sf
!          ELSE
!            Fac(4)= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor*w_sf
!          END IF ! usevMPF
!          Fac(1:3) = PartState(iPart,4:6)*Fac(4)
!          ! now, compute the actual 
!          CALL ComputeGaussDistance(N_in,r2_sf_inv,ShiftedPart,ElemDepo_xGP(:,:,:,:,iElem),GaussDistance)
!          DO k=0,PP_N
!            DO j=0,PP_N
!              DO i=0,PP_N
!                radius2=GaussDistance(i,j,k)
!                IF(radius2.GT.1.0) CYCLE
!  !              foundDOF=foundDOF+1
!                S = 1. - radius2
!                S1 = S*S
!                DO expo = 3, alpha_sf,1
!                  S1 = S*S1
!                END DO
!                source(1:3,i,j,k,iElem) = source(1:3,i,j,k,iElem) + Fac(1:3) * S1
!                source( 4 ,i,j,k,iElem) = source( 4 ,i,j,k,iElem) + Fac(4) * S1
!              END DO ! i=0,PP_N
!            END DO ! j=0,PP_N
!          END DO ! k=0,PP_N
!        END DO ! iPart=firstPart,LastPart
!        ! and for second step the external particles
!  #ifdef MPI
!        IF(.NOT.DoInnerParts)THEN
!          ! stuff to add, ExtPartToFIBGM
!          DO iPart=1,NbrOfextParticles  !external Particles
!            DO ind = 1,3
!              ShiftedPart(ind) = ExtPartState(iPart,ind) + casematrix(iCase,1)*Vec1(ind) + &
!                   casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
!            END DO
!            ! check if bounding FIBGM box of particle and element overlap
!            IF(ExtPartToFIBGM(1,iPart).GT.GEO%ElemToFIBGM(2,iElem)) CYCLE
!            IF(ExtPartToFIBGM(2,iPart).LT.GEO%ElemToFIBGM(1,iElem)) CYCLE
!            IF(ExtPartToFIBGM(3,iPart).GT.GEO%ElemToFIBGM(4,iElem)) CYCLE
!            IF(ExtPartToFIBGM(4,iPart).LT.GEO%ElemToFIBGM(3,iElem)) CYCLE
!            IF(ExtPartToFIBGM(5,iPart).GT.GEO%ElemToFIBGM(6,iElem)) CYCLE
!            IF(ExtPartToFIBGM(6,iPart).LT.GEO%ElemToFIBGM(5,iElem)) CYCLE
!            ! check if element circumcircle and shape function can interact
!            radius2= (ShiftedPart(1) - ElemBaryNGeo(1,iElem))* (ShiftedPart(1) - ElemBaryNGeo(1,iElem))&
!                   + (ShiftedPart(2) - ElemBaryNGeo(2,iElem))* (ShiftedPart(2) - ElemBaryNGeo(2,iElem))&
!                   + (ShiftedPart(3) - ElemBaryNGeo(3,iElem))* (ShiftedPart(3) - ElemBaryNGeo(3,iElem))
!            IF(radius2.GT.ElemRadius2_sf(iElem)) CYCLE
!            IF (usevMPF) THEN
!              Fac(4)= Species(ExtPartSpecies(iPart))%ChargeIC * ExtPartMPF(iPart)*w_sf
!            ELSE
!              Fac(4)= Species(ExtPartSpecies(iPart))%ChargeIC * Species(ExtPartSpecies(iPart))%MacroParticleFactor*w_sf
!            END IF ! usevMPF
!            Fac(1:3) = ExtPartState(iPart,4:6)*Fac(4)
!            ! now, compute the actual 
!            CALL ComputeGaussDistance(N_in,r2_sf_inv,ShiftedPart,ElemDepo_xGP(:,:,:,:,iElem),GaussDistance)
!            DO k=0,PP_N
!              DO j=0,PP_N
!                DO i=0,PP_N
!                  radius2=GaussDistance(i,j,k)
!                  IF(radius2.GT.1.0) CYCLE
!                  S = 1.0 - radius2
!                  S1 = S*S
!                  DO expo = 3, alpha_sf,1
!                    S1 = S*S1
!                  END DO
!    !              foundDOF=foundDOF+1
!                  source(1:3,i,j,k,iElem) = source(1:3,i,j,k,iElem) + Fac(1:3) * S1
!                  source( 4 ,i,j,k,iElem) = source( 4 ,i,j,k,iElem) + Fac(4) * S1
!                END DO ! i=0,PP_N
!              END DO ! j=0,PP_N
!            END DO ! k=0,PP_N
!          END DO ! iPart=firstPart,LastPart
!        END IF
!  #endif /*MPI*/
!      END DO ! iElem=1,PP_nElems
!    END DO ! iCase
!  #ifdef MPI
!    IF(.NOT.DoInnerParts)THEN
!      ! deallocate external state
!      SDEALLOCATE(ExtPartState)
!      SDEALLOCATE(ExtPartSpecies)
!      SDEALLOCATE(ExtPartToFIBGM)
!      SDEALLOCATE(ExtPartMPF)
!    END IF
!  #endif /*MPI*/
!    !IPWRITE(*,*) 'Number of found DOFs',FoundDOF

  IF( .NOT.DoInnerParts .AND. DoSFEqui) THEN
    ! map source from Equististant points on Gauss-Points
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(4,PP_N,PP_N,Vdm_EquiN_GaussN,source(:,:,:,:,iElem),source(:,:,:,:,iElem))
    END DO ! iElem=1,PP_nElems
  END IF



CASE('cylindrical_shape_function')
  IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
  Vec1(1:3) = 0.
  Vec2(1:3) = 0.
  Vec3(1:3) = 0.
  IF (GEO%nPeriodicVectors.EQ.1) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
  END IF
  IF (GEO%nPeriodicVectors.EQ.2) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
    Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
  END IF
  IF (GEO%nPeriodicVectors.EQ.3) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
    Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
    Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
  END IF
  DO iPart=firstPart,LastPart
    IF (PDM%ParticleInside(iPart)) THEN
      ! compute local radius
      local_r_sf= r_sf0 * (1.0 + r_sf_scale*DOT_PRODUCT(PartState(iPart,1:2),PartState(iPart,1:2)))
      local_r2_sf=local_r_sf*local_r_sf
      local_r2_sf_inv=1./local_r2_sf
      IF (usevMPF) THEN
        Fac(4)= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)*w_sf/(PI*(local_r_sf**3))
      ELSE
        Fac(4)= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor*w_sf/(PI*(local_r_sf**3))
      END IF ! usevMPF
      !IF(fac(4).GT.0.) print*,'charge pos'
      Fac(1:3) = PartState(iPart,4:6)*Fac(4)
      !-- determine which background mesh cells (and interpolation points within) need to be considered
      DO iCase = 1, NbrOfCases
        chargedone(:) = .FALSE.
        DO ind = 1,3
          ShiftedPart(ind) = PartState(iPart,ind) + casematrix(iCase,1)*Vec1(ind) + &
               casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
        END DO
        kmax = CEILING((ShiftedPart(1)+local_r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
        kmax = MIN(kmax,GEO%FIBGMimax)
        kmin = FLOOR((ShiftedPart(1)-local_r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
        kmin = MAX(kmin,GEO%FIBGMimin)
        lmax = CEILING((ShiftedPart(2)+local_r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
        lmax = MIN(lmax,GEO%FIBGMjmax)
        lmin = FLOOR((ShiftedPart(2)-local_r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
        lmin = MAX(lmin,GEO%FIBGMjmin)
        mmax = CEILING((ShiftedPart(3)+local_r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
        mmax = MIN(mmax,GEO%FIBGMkmax)
        mmin = FLOOR((ShiftedPart(3)-local_r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
        mmin = MAX(mmin,GEO%FIBGMkmin)
        !-- go through all these cells
        DO kk = kmin,kmax
          DO ll = lmin, lmax
            DO mm = mmin, mmax
              !--- go through all mapped elements not done yet
              DO ppp = 1,GEO%FIBGM(kk,ll,mm)%nElem
                ElemID = GEO%FIBGM(kk,ll,mm)%Element(ppp)
                IF(ElemID.GT.nElems) CYCLE
                IF (.NOT.chargedone(ElemID)) THEN
#ifdef MPI
                  nDeposPerElem(ElemID)=nDeposPerElem(ElemID)+1
#endif /*MPI*/
                  !--- go through all gauss points
                  !CALL ComputeGaussDistance(PP_N,r2_sf_inv,ShiftedPart,ElemDepo_xGP(:,:,:,:,ElemID),GaussDistance)
                  DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
                    !-- calculate distance between gauss and particle
                    radius2 = (ShiftedPart(1) - ElemDepo_xGP(1,k,l,m,ElemID)) * (ShiftedPart(1) - ElemDepo_xGP(1,k,l,m,ElemID)) &
                            + (ShiftedPart(2) - ElemDepo_xGP(2,k,l,m,ElemID)) * (ShiftedPart(2) - ElemDepo_xGP(2,k,l,m,ElemID)) &
                            + (ShiftedPart(3) - ElemDepo_xGP(3,k,l,m,ElemID)) * (ShiftedPart(3) - ElemDepo_xGP(3,k,l,m,ElemID))
                    !-- calculate charge and current density at ip point using a shape function
                    !-- currently only one shapefunction available, more to follow (including structure change)
                    IF (radius2 .LT. local_r2_sf) THEN
                      S = 1. - local_r2_sf_inv * radius2
                    !radius2=GaussDistance(k,l,m)
                    !IF (radius2 .LT. 1.0) THEN
                    !  S = 1 -  radius2
                      S1 = S*S
                      DO expo = 3, alpha_sf
                        S1 = S*S1
                      END DO
                      source(1:3,k,l,m,ElemID) = source(1:3,k,l,m,ElemID) + Fac(1:3) * S1
                      source( 4 ,k,l,m,ElemID) = source( 4 ,k,l,m,ElemID) + Fac(4) * S1
                    END IF
                  END DO; END DO; END DO
                  chargedone(ElemID) = .TRUE.
                END IF
              END DO ! ppp
            END DO ! mm
          END DO ! ll
        END DO ! kk
      END DO ! iCase (periodicity)
    END IF ! inside
  END DO ! i
#ifdef MPI           
  IF(.NOT.DoInnerParts)THEN
    Vec1(1:3) = 0.
    Vec2(1:3) = 0.
    Vec3(1:3) = 0.
    IF (NbrOfextParticles .GT. 0) THEN
      IF (GEO%nPeriodicVectors.EQ.1) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
      END IF
      IF (GEO%nPeriodicVectors.EQ.2) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
        Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
      END IF
      IF (GEO%nPeriodicVectors.EQ.3) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
        Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
        Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
      END IF
    END IF
    
    DO iPart=1,NbrOfextParticles  !external Particles
      local_r_sf= r_sf0 * (1.0 + r_sf_scale*DOT_PRODUCT(PartState(iPart,1:2),PartState(iPart,1:2)))
      local_r2_sf=local_r_sf*local_r_sf
      local_r2_sf_inv=1./local_r2_sf
      IF (usevMPF) THEN
        Fac(4)= Species(ExtPartSpecies(iPart))%ChargeIC * ExtPartMPF(iPart)*w_sf/(PI*(local_r_sf**3))
      ELSE
        Fac(4)= Species(ExtPartSpecies(iPart))%ChargeIC & 
              * Species(ExtPartSpecies(iPart))%MacroParticleFactor*w_sf/(PI*(local_r_sf**3))
      END IF ! usevMPF
      Fac(1:3) = ExtPartState(iPart,4:6)*Fac(4)
      !-- determine which background mesh cells (and interpolation points within) need to be considered
      DO iCase = 1, NbrOfCases
        chargedone(:) = .FALSE.
        DO ind = 1,3
          ShiftedPart(ind) = ExtPartState(iPart,ind) + casematrix(iCase,1)*Vec1(ind) + &
               casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
        END DO
        kmax = CEILING((ShiftedPart(1)+local_r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
        kmax = MIN(kmax,GEO%FIBGMimax)
        kmin = FLOOR((ShiftedPart(1)-local_r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
        kmin = MAX(kmin,GEO%FIBGMimin)
        lmax = CEILING((ShiftedPart(2)+local_r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
        lmax = MIN(lmax,GEO%FIBGMjmax)
        lmin = FLOOR((ShiftedPart(2)-local_r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
        lmin = MAX(lmin,GEO%FIBGMjmin)
        mmax = CEILING((ShiftedPart(3)+local_r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
        mmax = MIN(mmax,GEO%FIBGMkmax)
        mmin = FLOOR((ShiftedPart(3)-local_r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
        mmin = MAX(mmin,GEO%FIBGMkmin)
        !-- go through all these cells (should go through non if periodic and shiftedpart not in my domain
        DO kk = kmin,kmax
          DO ll = lmin, lmax
            DO mm = mmin, mmax
              !--- go through all mapped elements not done yet
              DO ppp = 1,GEO%FIBGM(kk,ll,mm)%nElem
                ElemID = GEO%FIBGM(kk,ll,mm)%Element(ppp)
                IF(ElemID.GT.nElems) CYCLE
                IF (.NOT.chargedone(ElemID)) THEN
#ifdef MPI
                  nDeposPerElem(ElemID)=nDeposPerElem(ElemID)+1
#endif /*MPI*/
                  !--- go through all gauss points
                  !CALL ComputeGaussDistance(PP_N,r2_sf_inv,ShiftedPart,ElemDepo_xGP(:,:,:,:,ElemID),GaussDistance)
                  DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
                    !-- calculate distance between gauss and particle
                      radius2 = (ShiftedPart(1) - ElemDepo_xGP(1,k,l,m,ElemID)) * (ShiftedPart(1) - ElemDepo_xGP(1,k,l,m,ElemID)) &
                              + (ShiftedPart(2) - ElemDepo_xGP(2,k,l,m,ElemID)) * (ShiftedPart(2) - ElemDepo_xGP(2,k,l,m,ElemID)) &
                              + (ShiftedPart(3) - ElemDepo_xGP(3,k,l,m,ElemID)) * (ShiftedPart(3) - ElemDepo_xGP(3,k,l,m,ElemID))
                      !-- calculate charge and current density at ip point using a shape function
                      !-- currently only one shapefunction available, more to follow (including structure change)
                      IF (radius2 .LT. local_r2_sf) THEN
                        S = 1. - local_r2_sf_inv * radius2
                        !IF(S.LT.0.) print*,'dist neg '
                      !radius2=GaussDistance(k,l,m)
                      !IF (radius2 .LT. 1.0) THEN
                      !  S = 1 -  radius2
                        S1 = S*S
                        DO expo = 3, alpha_sf
                          S1 = S*S1
                        END DO
                        source(1:3,k,l,m,ElemID) = source(1:3,k,l,m,ElemID) + Fac(1:3) * S1
                        source( 4 ,k,l,m,ElemID) = source( 4 ,k,l,m,ElemID) + Fac(4) * S1
                      END IF
                  END DO; END DO; END DO
                  chargedone(ElemID) = .TRUE.
                END IF
              END DO ! ppp
            END DO ! mm
          END DO ! ll
        END DO ! kk
      END DO
    END DO
    ! deallocate external state
    SDEALLOCATE(ExtPartState)
    SDEALLOCATE(ExtPartSpecies)
    SDEALLOCATE(ExtPartToFIBGM)
    SDEALLOCATE(ExtPartMPF)
    NbrOfExtParticles=0
  END IF
#endif /*MPI*/
  IF( .NOT.DoInnerParts .AND. DoSFEqui) THEN
    ! map source from Equististant points on Gauss-Points
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(4,PP_N,PP_N,Vdm_EquiN_GaussN,source(:,:,:,:,iElem),source(:,:,:,:,iElem))
    END DO ! iElem=1,PP_nElems
  END IF

CASE('delta_distri')
  IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
  DO iElem=1,PP_nElems
#ifdef MPI
    tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
    DO iPart=firstPart,LastPart
      IF (PDM%ParticleInside(iPart)) THEN
        IF(PEM%Element(iPart).EQ.iElem)THEN
          IF (usevMPF) THEN
            prefac= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)
          ELSE
            prefac= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor 
          END IF ! usevMPF
          ! Map Particle to -1|1 space (re-used in interpolation)
          IF(.NOT.DoRefMapping)THEN
            CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),iElem,iPart)
          END IF
          ! get value of test function at particle position
          SELECT CASE(DeltaType)
          CASE(1)
            ! xi   -direction
            CALL LagrangeInterpolationPolys(PartPosRef(1,iPart),NDepo,XiNDepo,wBaryNDepo,L_xi(1,:))
            ! eta  -direction                                 
            CALL LagrangeInterpolationPolys(PartPosRef(2,iPart),NDepo,XiNDepo,wBaryNDepo,L_xi(2,:))
            ! zeta -direction                                 
            CALL LagrangeInterpolationPolys(PartPosRef(3,iPart),NDepo,XiNDepo,wBaryNDepo,L_xi(3,:))
          CASE(2)
            DO i=0,NDepo
              ! xi   -direction
               CALL BernsteinPolynomial(NDepo,i,PartPosRef(1,iPart),L_xi(1,i),NDepoChooseK)
              ! eta  -direction                                  
               CALL BernsteinPolynomial(NDepo,i,PartPosRef(2,iPart),L_xi(2,i),NDepoChooseK)
              ! zeta  -direction                                 
               CALL BernsteinPolynomial(NDepo,i,PartPosRef(3,iPart),L_xi(3,i),NDepoChooseK)
            END DO ! i
          CASE(3)
            ! xi - direction
            CALL DeBoorRef(NDepo,NKnots,Knots,PartPosRef(1,iPart),L_xi(1,:))
            ! eta - direction                                  
            CALL DeBoorRef(NDepo,NKnots,Knots,PartPosRef(2,iPart),L_xi(2,:))
            ! zeta - direction                                 
            CALL DeBoorRef(NDepo,NKnots,Knots,PartPosRef(3,iPart),L_xi(3,:))
          END SELECT
          DO k=0,NDepo
            DO j=0,NDepo
              DO i=0,NDepo
           !     print*,'i,j,k,L',i,j,k,L_xi(1,i)* L_xi(2,j)* L_xi(3,k)
                DeltaIntCoeff = L_xi(1,i)* L_xi(2,j)* L_xi(3,k)*prefac 
                source(1:3,i,j,k,iElem) = source(1:3,i,j,k,iElem) + DeltaIntCoeff*PartState(iPart,4:6)
                source( 4 ,i,j,k,iElem) = source( 4 ,i,j,k,iElem) + DeltaIntCoeff
              END DO ! i
            END DO ! j
          END DO ! k
        END IF ! Particle in Element
      END IF ! ParticleInside of domain
    END DO ! ParticleVecLength
#ifdef MPI
    tLBEnd = LOCALTIME() ! LB Time End
    ElemTime(iElem)=ElemTime(iElem)+tLBEnd-tLBStart
#endif /*MPI*/
  END DO ! iElem
  IF(.NOT.DoInnerParts)THEN
    DO iElem=1,PP_nElems
#ifdef MPI
      tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
      DO k=0,NDepo
        DO j=0,NDepo
          DO i=0,NDepo
            source( : ,i,j,k,iElem) = source( : ,i,j,k,iElem) *sJNDepo(i,j,k,iElem)*swGPNDepo(i)*swGPNDepo(j)*swGPNDepo(k)
          END DO ! i
        END DO ! j
      END DO ! k
      IF(DoChangeBasis)THEN
        CALL ChangeBasis3D(4,NDepo,PP_N,Vdm_NDepo_GaussN,source(1:4,0:NDepo,0:NDepo,0:NDepo,iElem)&
                                                        ,source(1:4,0:PP_N ,0:PP_N ,0:PP_N, iElem))
      END IF
!      IF(DoChangeBasis)THEN
!        CALL ChangeBasis3D(4,NDepo,PP_N,Vdm_NDepo_GaussN,source(:,0:NDepo,0:NDepo,0:NDepo,iElem)&
!                                                        ,source(:,0:PP_N ,0:PP_N ,0:PP_N, iElem))
!      END IF
!      DO k=0,PP_N
!        DO j=0,PP_N
!          DO i=0,PP_N
!            source( : ,i,j,k,iElem) = source( : ,i,j,k,iElem) *sJ(i,j,k,iElem)*swGP(i)*swGP(j)*swGP(k)
!          END DO ! i
!        END DO ! j
!      END DO ! k
#ifdef MPI
      tLBEnd = LOCALTIME() ! LB Time End
      ElemTime(iElem)=ElemTime(iElem)+tLBEnd-tLBStart
#endif /*MPI*/
    END DO ! loop over all elems
  END IF ! DoInnerParts
CASE('nearest_gausspoint')
  IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
  SAVE_GAUSS = .FALSE.
  IF(TRIM(InterpolationType).EQ.'nearest_gausspoint') SAVE_GAUSS = .TRUE.
  IF(MOD(PP_N,2).EQ.0) THEN
    a = PP_N/2
    b = a
  ELSE
    a = (PP_N+1)/2
    b = a-1
  END IF
  DO iElem=1,PP_nElems
#ifdef MPI
    tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
    DO iPart=firstPart,LastPart
      IF (PDM%ParticleInside(iPart)) THEN
        IF(PEM%Element(iPart).EQ.iElem)THEN
          IF (usevMPF) THEN
            prefac= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)
          ELSE
            prefac= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor 
          END IF ! usevMPF
          ! Map Particle to -1|1 space (re-used in interpolation)
          !IF(.NOT.DoRefMapping)THEN
          !  CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),iElem,iPart)
          !END IF
          ! Find out which gausspoint is closest and add up charges and currents
          !! x-direction
          IF(.NOT.SAVE_GAUSS) THEN
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
          END IF
!#if (PP_nVar==8)
          source(1:3,k,l,m,iElem) = source(1:3,k,l,m,iElem) + PartState(iPart,4:6) * prefac
!#endif
          source( 4 ,k,l,m,iElem) = source( 4 ,k,l,m,iElem) + prefac
          !IF (SAVE_GAUSS) THEN
          !  PartPosGauss(iPart,1) = k
          !  PartPosGauss(iPart,2) = l
          !  PartPosGauss(iPart,3) = m
          !END IF
        END IF ! Element .EQ. iElem
      END IF ! Particle inside
    END DO ! iPart
#ifdef MPI
    tLBEnd = LOCALTIME() ! LB Time End
    ElemTime(iElem)=ElemTime(iElem)+tLBEnd-tLBStart
#endif /*MPI*/
  END DO ! iElem=1,PP_nElems
  IF(.NOT.DoInnerParts)THEN
    DO iElem=1,PP_nElems
#ifdef MPI
      tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
      DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
      ! get densities by dividing by gauss volume
!#if (PP_nVar==8)
        source(1:4,k,l,m,iElem) = source(1:4,k,l,m,iElem) * sJ(k,l,m,iElem)/(wGP(k)*wGP(l)*wGP(m))
!#else
!        source(4,k,l,m,iElem) = source(4,k,l,m,iElem) * sJ(k,l,m,iElem)/(wGP(k)*wGP(l)*wGP(m))
!#endif
      END DO; END DO; END DO
#ifdef MPI
      tLBEnd = LOCALTIME() ! LB Time End
      ElemTime(iElem)=ElemTime(iElem)+tLBEnd-tLBStart
#endif /*MPI*/
    END DO ! iElem=1,PP_nElems
  END IF
CASE('cartmesh_volumeweighting')
  ! Step 1: Deposition of all particles onto background mesh -> densities
  ! IF(DoInnerParts) BGMSource=0.0 ! not possible due to periodic stuff --> two communications
#ifdef MPI
  tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
  BGMSource(:,:,:,:) = 0.0
  DO iPart = firstPart, lastPart
    IF (PDM%ParticleInside(iPart)) THEN
      IF (usevMPF) THEN
        Charge= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)
      ELSE
        Charge= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor 
      END IF ! usevMPF
      !Charge = Species(PartSpecies(iPart))%ChargeIC*Species(PartSpecies(iPart))%MacroParticleFactor
      k = FLOOR(PartState(iPart,1)/BGMdeltas(1))
      l = FLOOR(PartState(iPart,2)/BGMdeltas(2))
      m = FLOOR(PartState(iPart,3)/BGMdeltas(3))
      alpha1 = (PartState(iPart,1) / BGMdeltas(1)) - k
      alpha2 = (PartState(iPart,2) / BGMdeltas(2)) - l
      alpha3 = (PartState(iPart,3) / BGMdeltas(3)) - m
      TSource(:) = 0.0
!#if (PP_nVar==8)
      TSource(1) = PartState(iPart,4)*Charge
      TSource(2) = PartState(iPart,5)*Charge
      TSource(3) = PartState(iPart,6)*Charge
!#endif
      TSource(4) = Charge

      BGMSource(k,l,m,1:4)       = BGMSource(k,l,m,1:4) + (TSource * (1-alpha1)*(1-alpha2)*(1-alpha3))
      BGMSource(k,l,m+1,1:4)     = BGMSource(k,l,m+1,1:4) + (TSource * (1-alpha1)*(1-alpha2)*(alpha3))
      BGMSource(k,l+1,m,1:4)     = BGMSource(k,l+1,m,1:4) + (TSource * (1-alpha1)*(alpha2)*(1-alpha3))
      BGMSource(k,l+1,m+1,1:4)   = BGMSource(k,l+1,m+1,1:4) + (TSource * (1-alpha1)*(alpha2)*(alpha3))
      BGMSource(k+1,l,m,1:4)     = BGMSource(k+1,l,m,1:4) + (TSource * (alpha1)*(1-alpha2)*(1-alpha3))
      BGMSource(k+1,l,m+1,1:4)   = BGMSource(k+1,l,m+1,1:4) + (TSource * (alpha1)*(1-alpha2)*(alpha3))
      BGMSource(k+1,l+1,m,1:4)   = BGMSource(k+1,l+1,m,1:4) + (TSource * (alpha1)*(alpha2)*(1-alpha3))
      BGMSource(k+1,l+1,m+1,1:4) = BGMSource(k+1,l+1,m+1,1:4) + (TSource * (alpha1)*(alpha2)*(alpha3))
    END IF
  END DO
  BGMSource(:,:,:,:) = BGMSource(:,:,:,:) / BGMVolume

#ifdef MPI
  ! should be treated in this way, unforunately, we would negelct the periodic stuff
  !IF(.NOT.DoInnerParts)
  CALL MPISourceExchangeBGM()
#else
  IF (GEO%nPeriodicVectors.GT.0) CALL PeriodicSourceExchange()
#endif
#ifdef MPI
  tLBEnd = LOCALTIME() ! LB Time End
  tCartMesh=tCartMesh+tLBEnd-tLBStart
#endif /*MPI*/

  ! Step 2: Interpolation of densities onto grid
  DO iElem = 1, nElems
#ifdef MPI
    tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
    DO kk = 0, PP_N
      DO ll = 0, PP_N
        DO mm = 0, PP_N
         k = GaussBGMIndex(1,kk,ll,mm,iElem)
         l = GaussBGMIndex(2,kk,ll,mm,iElem)
         m = GaussBGMIndex(3,kk,ll,mm,iElem)
         alpha1 = GaussBGMFactor(1,kk,ll,mm,iElem)
         alpha2 = GaussBGMFactor(2,kk,ll,mm,iElem)
         alpha3 = GaussBGMFactor(3,kk,ll,mm,iElem)
!#if (PP_nVar==8)
         DO i = 1,3
           source(i,kk,ll,mm,iElem) = source(i,kk,ll,mm,iElem)            + &
                BGMSource(k,l,m,i) * (1-alpha1) * (1-alpha2) * (1-alpha3) + &
                BGMSource(k,l,m+1,i) * (1-alpha1) * (1-alpha2) * (alpha3) + &
                BGMSource(k,l+1,m,i) * (1-alpha1) * (alpha2) * (1-alpha3) + &
                BGMSource(k,l+1,m+1,i) * (1-alpha1) * (alpha2) * (alpha3) + &
                BGMSource(k+1,l,m,i) * (alpha1) * (1-alpha2) * (1-alpha3) + &
                BGMSource(k+1,l,m+1,i) * (alpha1) * (1-alpha2) * (alpha3) + &
                BGMSource(k+1,l+1,m,i) * (alpha1) * (alpha2) * (1-alpha3) + &
                BGMSource(k+1,l+1,m+1,i) * (alpha1) * (alpha2) * (alpha3)
         END DO
!#endif
           source(4,kk,ll,mm,iElem) = source(4,kk,ll,mm,iElem)          + &
              BGMSource(k,l,m,4) * (1-alpha1) * (1-alpha2) * (1-alpha3) + &
              BGMSource(k,l,m+1,4) * (1-alpha1) * (1-alpha2) * (alpha3) + &
              BGMSource(k,l+1,m,4) * (1-alpha1) * (alpha2) * (1-alpha3) + &
              BGMSource(k,l+1,m+1,4) * (1-alpha1) * (alpha2) * (alpha3) + &
              BGMSource(k+1,l,m,4) * (alpha1) * (1-alpha2) * (1-alpha3) + &
              BGMSource(k+1,l,m+1,4) * (alpha1) * (1-alpha2) * (alpha3) + &
              BGMSource(k+1,l+1,m,4) * (alpha1) * (alpha2) * (1-alpha3) + &
              BGMSource(k+1,l+1,m+1,4) * (alpha1) * (alpha2) * (alpha3)
       END DO !mm
     END DO !ll
   END DO !kk
#ifdef MPI
   tLBEnd = LOCALTIME() ! LB Time End
   ElemTime(iElem)=ElemTime(iElem)+tLBEnd-tLBStart
#endif /*MPI*/
 END DO !iElem
 !DEALLOCATE(BGMSource)
CASE('cartmesh_splines')
  ! Step 1: Deposition of all particles onto background mesh -> densities
  !ALLOCATE(BGMSource(BGMminX:BGMmaxX,BGMminY:BGMmaxY,BGMminZ:BGMmaxZ,1:4))
  ! IF(DoInnerParts) BGMSource=0. not possible due to periodic stuff
#ifdef MPI
  tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
  BGMSource(:,:,:,:) = 0.0
  DO iPart = firstPart, lastPart
    IF (PDM%ParticleInside(iPart)) THEN
!      Charge = Species(PartSpecies(iPart))%ChargeIC*Species(PartSpecies(iPart))%MacroParticleFactor
      IF (usevMPF) THEN
        Charge= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)
      ELSE
        Charge= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor 
      END IF ! usevMPF
      PosInd(1) = FLOOR(PartState(iPart,1)/BGMdeltas(1))
      PosInd(2) = FLOOR(PartState(iPart,2)/BGMdeltas(2))
      PosInd(3) = FLOOR(PartState(iPart,3)/BGMdeltas(3))
      !print*,'posind(1:3),charge',posInd,charge
      DO dir = 1,3               ! x,y,z direction
        DO weightrun = 0,3
          DO mm = 0, 3
            IF (mm.EQ.weightrun) then
              auxiliary(mm) = 1.0
            ELSE
              auxiliary(mm) = 0.0
            END IF
          END DO
          CALL DeBoor(PosInd(dir),auxiliary,PartState(iPart,dir),weight(dir,weightrun),dir)
        END DO
      END DO
      DO k = PosInd(1)-1, PosInd(1)+2
        kk = abs(k - PosInd(1) - 2)
        DO l = PosInd(2)-1, PosInd(2)+2
          ll = abs(l - PosInd(2) - 2)
          DO m = PosInd(3)-1, PosInd(3)+2
            mm = abs(m - PosInd(3) - 2)
            locweight = weight(1,kk)*weight(2,ll)*weight(3,mm)*charge
!#if (PP_nVar==8)
            BGMSource(k,l,m,1) = BGMSource(k,l,m,1) + PartState(iPart,4)* locweight
            BGMSource(k,l,m,2) = BGMSource(k,l,m,2) + PartState(iPart,5)* locweight
            BGMSource(k,l,m,3) = BGMSource(k,l,m,3) + PartState(iPart,6)* locweight
!#endif
            BGMSource(k,l,m,4) = BGMSource(k,l,m,4) + locweight
         !   print*,'BMGSOURCE4',BGMSOURCE(k,l,m,4)
          END DO
        END DO
      END DO
    END IF
  END DO
  BGMSource(:,:,:,:) = BGMSource(:,:,:,:) / BGMVolume

#ifdef MPI
  !IF(.NOT.DoInnerParts)THEN has to be communicated each time :(
  CALL MPISourceExchangeBGM()
#else
  IF (GEO%nPeriodicVectors.GT.0) CALL PeriodicSourceExchange()
#endif
#ifdef MPI
  tLBEnd = LOCALTIME() ! LB Time End
  tCartMesh=tCartMesh+tLBEnd-tLBStart
#endif /*MPI*/

  ! Step 2: Interpolation of densities onto grid
  DO iElem = 1, nElems
#ifdef MPI
    tLBStart = LOCALTIME() ! LB Time Start
#endif /*MPI*/
    DO kk = 0, PP_N
      DO ll = 0, PP_N
        DO mm = 0, PP_N
          k = GaussBGMIndex(1,kk,ll,mm,iElem)
          l = GaussBGMIndex(2,kk,ll,mm,iElem)
          m = GaussBGMIndex(3,kk,ll,mm,iElem)
          DO r = k-1,k+2
            u = r-k+2
            DO ss = l-1,l+2
              v = ss-l+2
              DO t = m-1,m+2
                w = t-m+2
!#if (PP_nVar==8)                  
                source(1:4,kk,ll,mm,iElem) = source(1:4,kk,ll,mm,iElem) + BGMSource(r,ss,t,1:4) * GPWeight(iElem,kk,ll,mm,u,v,w)
                !DO i = 1,3
                !  source(i,kk,ll,mm,iElem) = source(i,kk,ll,mm,iElem) + BGMSource(r,ss,t,i) * GPWeight(iElem,kk,ll,mm,u,v,w)
                !END DO
!#endif
                !source(4,kk,ll,mm,iElem) = source(4,kk,ll,mm,iElem) + BGMSource(r,ss,t,4) * GPWeight(iElem,kk,ll,mm,u,v,w)
              END DO !t
            END DO !s
          END DO !r
        END DO !mm
      END DO !ll
    END DO !kk
#ifdef MPI
    tLBEnd = LOCALTIME() ! LB Time End
    ElemTime(iElem)=ElemTime(iElem)+tLBEnd-tLBStart
#endif /*MPI*/
  END DO !iElem
 !DEALLOCATE(BGMSource)
CASE DEFAULT
  CALL abort(__STAMP__, &
       'Unknown DepositionType in pic_depo.f90')
END SELECT

RETURN
END SUBROUTINE Deposition


#ifndef MPI
SUBROUTINE PeriodicSourceExchange()    
!============================================================================================================================
! Exchange sources in periodic case
!============================================================================================================================
! use MODULES                                                    
USE MOD_PICDepo_Vars
USE MOD_Particle_Mesh_Vars,  ONLY: GEO
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE                                                                                  
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!REAL,INTENT(INOUT)         :: BGMSource(BGMminX:BGMmaxX,BGMminY:BGMmaxY,BGMminZ:BGMmaxZ,1:4)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                                           
INTEGER                     :: i,k,l,m,k2,l2,m2
!-----------------------------------------------------------------------------------------------------------------------------------

DO i = 1,GEO%nPeriodicVectors
  DO k = BGMminX, BGMmaxX
    k2 = k + GEO%PeriodicBGMVectors(1,i)
    DO l = BGMminY, BGMmaxY
      l2 = l + GEO%PeriodicBGMVectors(2,i)
      DO m = BGMminZ, BGMmaxZ
        m2 = m + GEO%PeriodicBGMVectors(3,i)
        IF ((k2.GE.BGMminX).AND.(k2.LE.BGMmaxX)) THEN
          IF ((l2.GE.BGMminY).AND.(l2.LE.BGMmaxY)) THEN
            IF ((m2.GE.BGMminZ).AND.(m2.LE.BGMmaxZ)) THEN
              BGMSource(k,l,m,:) = BGMSource(k,l,m,:) + BGMSource(k2,l2,m2,:)
              BGMSource(k2,l2,m2,:) = BGMSource(k,l,m,:)
            END IF
          END IF
        END IF
      END DO
    END DO
  END DO
END DO

END SUBROUTINE PeriodicSourceExchange
#else /*MPI*/
SUBROUTINE MPISourceExchangeBGM()            
!=================================================================================================================================
! Exchange sources in periodic case for MPI
!==================================================================================================================================
! use MODULES                                                         
USE MOD_Particle_MPI_Vars,  ONLY: PartMPI,tMPIMEssage
USE MOD_Particle_Mesh_Vars, ONLY: GEO
USE MOD_PICDepo_Vars
USE MOD_Globals
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE                                                                                  
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!REAL,INTENT(INOUT)        :: BGMSource(BGMminX:BGMmaxX,BGMminY:BGMmaxY,BGMminZ:BGMmaxZ,1:4)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tMPIMessage)           :: send_message(0:PartMPI%nProcs-1)                              
TYPE(tMPIMessage)           :: recv_message(0:PartMPI%nProcs-1)                              
INTEGER                     :: send_request(0:PartMPI%nProcs-1)                              
INTEGER                     :: recv_request(0:PartMPI%nProcs-1)                              
INTEGER                     :: send_status_list(1:MPI_STATUS_SIZE,0:PartMPI%nProcs-1)        
INTEGER                     :: recv_status_list(1:MPI_STATUS_SIZE,0:PartMPI%nProcs-1)        
INTEGER                     :: iProc,k,l,m,n, ppp, Counter, MsgLength(0:PartMPI%nProcs-1), iPer
INTEGER                     :: SourceLength(0:PartMPI%nProcs-1)                              
INTEGER                     :: RecvLength(0:PartMPI%nProcs-1)                                
INTEGER                     :: allocStat, Counter2                                     
INTEGER                     :: messageCounterS, messageCounterR                                
INTEGER                     :: myRealKind, k2,l2,m2                                          
REAL                        :: myRealTestValue                                                
!-----------------------------------------------------------------------------------------------------------------------------------

myRealKind = KIND(myRealTestValue)
IF (myRealKind.EQ.4) THEN
 myRealKind = MPI_REAL
ELSE IF (myRealKind.EQ.8) THEN
 myRealKind = MPI_DOUBLE_PRECISION
ELSE
 myRealKind = MPI_REAL
END IF
     
!--- Determine which Sources actually need to be sent (<> 0) and build corresponding true/false list
!    One list per process
DO iProc = 0,PartMPI%nProcs-1
   MsgLength(iProc) = 0
   IF (PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor) THEN
      MsgLength(iProc) = (PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,1) - PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,1) + 1) * &
                         (PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,2) - PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,2) + 1) * &
                         (PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,3) - PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,3) + 1)
   END IF
   IF(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor) THEN
      DO k = 1, PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount
         MsgLength(iProc) = MsgLength(iProc) + &
              (PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(2,1) -&
               PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(1,1) + 1) * &
              (PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(2,2) -&
               PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(1,2) + 1) * &
              (PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(2,3) -&
               PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(1,3) + 1)
      END DO
   END IF
   IF((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor))THEN
      ALLOCATE(send_message(iProc)%content_log(1:MsgLength(iProc)), STAT=allocStat)
      IF (allocStat .NE. 0) THEN
         CALL abort(__STAMP__, &
            'ERROR in MPISourceExchangeBGM: cannot allocate send_message')
      END IF
      ALLOCATE(recv_message(iProc)%content_log(1:MsgLength(iProc)), STAT=allocStat)
      IF (allocStat .NE. 0) THEN
         CALL abort(__STAMP__, &
            'ERROR in MPISourceExchangeBGM: cannot allocate recv_message')
      END IF
   END IF
   !--- check which sources are <> 0
   Counter = 0
   SourceLength(iProc) = 0
   IF (PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor) THEN
      DO k = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,1), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,1)
        DO l = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,2), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,2)
          DO m = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,3), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,3)
            Counter = Counter + 1
            IF (ANY(BGMSource(k,l,m,:).NE.0.0)) THEN
               send_message(iProc)%content_log(Counter) = .TRUE.
               SourceLength(iProc) = SourceLength(iProc) + 1
            ELSE
               send_message(iProc)%content_log(Counter) = .FALSE.
            END IF
          END DO
        END DO
      END DO
   END IF
   !--- same for periodic
   IF(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor) THEN
      DO n = 1, PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount
         DO k = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,1),&
                PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,1)
          DO l = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,2),&
                 PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,2)
            DO m = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,3),&
                   PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,3)
               Counter = Counter + 1
               IF (ANY(BGMSource(k,l,m,:).NE.0.0)) THEN
                  send_message(iProc)%content_log(Counter) = .TRUE.
                  SourceLength(iProc) = SourceLength(iProc) + 1
               ELSE
                  send_message(iProc)%content_log(Counter) = .FALSE.
               END IF
            END DO
          END DO
         END DO
      END DO
   END IF
END DO
!--- communicate
messageCounterS = 0
DO iProc = 0,PartMPI%nProcs-1
   IF ((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)) THEN
      ! MPI_ISEND true/false list for all border BGM points
      messageCounterS = messageCounterS + 1
      CALL MPI_ISEND(send_message(iProc)%content_log,MsgLength(iProc),MPI_LOGICAL,iProc,1,PartMPI%COMM, &
                     send_request(messageCounterS), IERROR)
   END IF
END DO
messageCounterR = 0
DO iProc = 0,PartMPI%nProcs-1
   IF ((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)) THEN
      ! MPI_IRECV true/false list for all border BGM points from neighbor CPUs
      messageCounterR = messageCounterR + 1
      CALL MPI_IRECV(recv_message(iProc)%Content_log,MsgLength(iProc),MPI_LOGICAL,iProc,1,PartMPI%COMM, &
                     recv_request(messageCounterR), IERROR)
   END IF
END DO
! MPI_WAITALL for the non-blocking MPI-communication to be finished
IF (messageCounterS .GE. 1) THEN
   CALL MPI_WAITALL(messageCounterS,send_request(1:messageCounterS),send_status_list(:,1:messageCounterS),IERROR)
END IF
IF (messageCounterR .GE. 1) THEN
   CALL MPI_WAITALL(messageCounterR,recv_request(1:messageCounterR),recv_status_list(:,1:messageCounterR),IERROR)
END IF

!--- Assemble actual sources to send
 DO iProc = 0,PartMPI%nProcs-1
   IF (((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.&
        (PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)).AND.(SourceLength(iProc).GT.0)) THEN
      ALLOCATE(send_message(iProc)%content(1:SourceLength(iProc)*4), STAT=allocStat)
      IF (allocStat .NE. 0) THEN
         CALL abort(__STAMP__, &
            'ERROR in MPISourceExchangeBGM: cannot allocate send_message')
      END IF
   END IF
   Counter = 0
   Counter2 = 0
   IF (PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor) THEN
      DO k = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,1), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,1)
        DO l = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,2), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,2)
          DO m = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,3), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,3)
             Counter = Counter + 1
             IF (send_message(iProc)%content_log(Counter)) THEN
                Counter2 = Counter2 + 1
                DO n = 1,4
                   send_message(iProc)%content((Counter2-1)*4 +n) = BGMSource(k,l,m,n)
                END DO
             END IF
          END DO
        END DO
      END DO
   END IF
   IF (PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor) THEN
      DO n = 1, PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount
         DO k = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,1),&
                PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,1)
           DO l = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,2),&
                  PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,2)
             DO m = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,3),&
                    PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,3)
                Counter = Counter + 1
                IF (send_message(iProc)%content_log(Counter)) THEN
                   Counter2 = Counter2 + 1
                   DO ppp = 1,4
                      send_message(iProc)%content((Counter2-1)*4 +ppp) = BGMSource(k,l,m,ppp)
                   END DO
                END IF
             END DO
           END DO
         END DO
      END DO
   END IF
END DO

!--- allocate actual source receive buffer
DO iProc = 0,PartMPI%nProcs-1
   IF ((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)) THEN
      Counter = 0
      DO k = 1, MsgLength(iProc)
         IF (recv_message(iProc)%Content_log(k)) THEN
            Counter = Counter + 1
         END IF
      END DO
      RecvLength(iProc) = Counter
      IF (RecvLength(iProc).GT.0) THEN
         ALLOCATE(recv_message(iProc)%content(1:Counter*4), STAT=allocStat)
         IF (allocStat .NE. 0) THEN
            CALL abort(__STAMP__, &
              'ERROR in MPISourceExchangeBGM: cannot allocate recv_message')
         END IF
      END IF
   END IF
END DO
!--- communicate
messageCounterS = 0
DO iProc = 0,PartMPI%nProcs-1
   IF (((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.&
        (PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)).AND.(SourceLength(iProc).GT.0)) THEN
      ! MPI_ISEND true/false list for all border BGM points
      messageCounterS = messageCounterS + 1
      CALL MPI_ISEND(send_message(iProc)%content,SourceLength(iProc)*4,myRealKind,iProc,1,PartMPI%COMM, &
                     send_request(messageCounterS), IERROR)
   END IF
END DO
messageCounterR = 0
DO iProc = 0,PartMPI%nProcs-1
   IF (((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.&
        (PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)).AND.(RecvLength(iProc).GT.0)) THEN
      ! MPI_IRECV true/false list for all border BGM points from neighbor CPUs
      messageCounterR = messageCounterR + 1
      CALL MPI_IRECV(recv_message(iProc)%content,RecvLength(iProc)*4,myRealKind,iProc,1,PartMPI%COMM, &
                     recv_request(messageCounterR), IERROR)
   END IF
END DO
! MPI_WAITALL for the non-blocking MPI-communication to be finished
IF (messageCounterS .GE. 1) THEN
   CALL MPI_WAITALL(messageCounterS,send_request(1:messageCounterS),send_status_list(:,1:messageCounterS),IERROR)
END IF
IF (messageCounterR .GE. 1) THEN
   CALL MPI_WAITALL(messageCounterR,recv_request(1:messageCounterR),recv_status_list(:,1:messageCounterR),IERROR)
END IF
!--- Deallocate Send Message Buffers
DO iProc = 0,PartMPI%nProcs-1
   IF ((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)) THEN
      DEALLOCATE(send_message(iProc)%content_log, STAT=allocStat)
      IF (allocStat .NE. 0) THEN
         CALL abort(__STAMP__, &
            'ERROR in MPISourceExchangeBGM: cannot deallocate send_message')
      END IF
      IF (SourceLength(iProc).GT.0) THEN
         DEALLOCATE(send_message(iProc)%content, STAT=allocStat)
         IF (allocStat .NE. 0) THEN
            CALL abort(__STAMP__, &
              'ERROR in MPISourceExchangeBGM: cannot deallocate send_message')
         END IF
      END IF
   END IF
END DO

!--- add selfperiodic sources, if any (needs to be done after send message is compiled and before
!---           received sources have been added!
IF ((GEO%nPeriodicVectors.GT.0).AND.(GEO%SelfPeriodic)) THEN
   DO iPer = 1, GEO%nPeriodicVectors
      DO k = BGMminX, BGMmaxX
         k2 = k + GEO%PeriodicBGMVectors(1,iPer)
        DO l = BGMminY, BGMmaxY
          l2 = l + GEO%PeriodicBGMVectors(2,iPer)
          DO m = BGMminZ, BGMmaxZ
             m2 = m + GEO%PeriodicBGMVectors(3,iPer)
             IF ((k2.GE.BGMminX).AND.(k2.LE.BGMmaxX)) THEN
             IF ((l2.GE.BGMminY).AND.(l2.LE.BGMmaxY)) THEN
             IF ((m2.GE.BGMminZ).AND.(m2.LE.BGMmaxZ)) THEN
                BGMSource(k,l,m,:) = BGMSource(k,l,m,:) + BGMSource(k2,l2,m2,:)
                BGMSource(k2,l2,m2,:) = BGMSource(k,l,m,:)
             END IF
             END IF
             END IF
          END DO
        END DO
      END DO
   END DO
END IF

!--- Add Sources and Deallocate Receive Message Buffers
DO iProc = 0,PartMPI%nProcs-1
   IF (RecvLength(iProc).GT.0) THEN
      Counter = 0
      Counter2 = 0
      IF (PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor) THEN
         DO k = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,1), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,1)
         DO l = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,2), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,2)
         DO m = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,3), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,3)
            Counter = Counter + 1
            IF(recv_message(iProc)%content_log(Counter))THEN
               Counter2 = Counter2 + 1
               DO n = 1,4
                 BGMSource(k,l,m,n) = BGMSource(k,l,m,n) + recv_message(iProc)%content((Counter2-1)*4+n)
               END DO
            END IF
         END DO
         END DO
         END DO
      END IF
      IF (PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor) THEN
         DO n = 1, PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount
           DO k = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,1),&
                  PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,1)
             DO l = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,2),&
                    PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,2)
               DO m = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,3),&
                      PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,3)
                  Counter = Counter + 1
                  IF(recv_message(iProc)%content_log(Counter))THEN
                     Counter2 = Counter2 + 1
                     DO ppp = 1,4
                        BGMSource(k,l,m,ppp) = BGMSource(k,l,m,ppp) + recv_message(iProc)%content((Counter2-1)*4+ppp)
                     END DO
                  END IF
               END DO
             END DO
           END DO
         END DO
      END IF
      IF ((PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor)) THEN
         DEALLOCATE(recv_message(iProc)%content, STAT=allocStat)
         IF (allocStat .NE. 0) THEN
            CALL abort(__STAMP__, &
              'ERROR in MPISourceExchangeBGM: cannot deallocate recv_message')
         END IF
      END IF
   END IF
   IF ((PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor)) THEN
      DEALLOCATE(recv_message(iProc)%content_log, STAT=allocStat)
      IF (allocStat .NE. 0) THEN
         CALL abort(__STAMP__, &
          'ERROR in MPISourceExchangeBGM: cannot deallocate recv_message')
      END IF
   END IF
END DO

END SUBROUTINE MPISourceExchangeBGM


SUBROUTINE MPIBackgroundMeshInit()  
!==================================================================================================================================
! initialize MPI background mesh
!==================================================================================================================================
! use MODULES                                                                   
USE MOD_PICDepo_Vars
USE MOD_Globals
! USE MOD_Particle_Vars,       ONLY:
USE MOD_Particle_Mesh_Vars,  ONLY:GEO
USE MOD_Particle_MPI_Vars,   ONLY:PartMPI
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT NONE                                                                                  
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iProc,jPorc
INTEGER                     :: k,m,n,m0,n0
INTEGER                     :: localminmax(6), maxofmin, minofmax                              
INTEGER                     :: completeminmax(6*PartMPI%nProcs)                              
INTEGER                     :: allocStat, NeighCount                                                                 
INTEGER                     :: TempBorder(1:2,1:3)            
INTEGER                     :: Periodicminmax(6), coord, PeriodicVec(1:3)                   
INTEGER                     :: TempPeriBord(1:26,1:2,1:3)                                      
LOGICAL                     :: CHECKNEIGHBOR     
!-----------------------------------------------------------------------------------------------------------------------------------

!Periodic Init stuff
IF(GEO%nPeriodicVectors.GT.0)THEN
  ! Compute PeriodicBGMVectors (from PeriodicVectors and BGMdeltas)
  ALLOCATE(GEO%PeriodicBGMVectors(1:3,1:GEO%nPeriodicVectors),STAT=allocStat)
  IF (allocStat .NE. 0) THEN
    CALL abort(__STAMP__, &
       'ERROR in MPIBackgroundMeshInit: cannot allocate GEO%PeriodicBGMVectors!')
  END IF
  DO iProc = 1, GEO%nPeriodicVectors
    GEO%PeriodicBGMVectors(1,iProc) = NINT(GEO%PeriodicVectors(1,iProc)/BGMdeltas(1))
    IF(ABS(GEO%PeriodicVectors(1,iProc)/BGMdeltas(1)-REAL(GEO%PeriodicBGMVectors(1,iProc))).GT.1E-10)THEN
      CALL abort(__STAMP__, &
       'ERROR: Periodic Vector ist not multiple of background mesh delta')
    END IF
    GEO%PeriodicBGMVectors(2,iProc) = NINT(GEO%PeriodicVectors(2,iProc)/BGMdeltas(2))
    IF(ABS(GEO%PeriodicVectors(2,iProc)/BGMdeltas(2)-REAL(GEO%PeriodicBGMVectors(2,iProc))).GT.1E-10)THEN
      CALL abort(__STAMP__, &
       'ERROR: Periodic Vector ist not multiple of background mesh delta')
    END IF
    GEO%PeriodicBGMVectors(3,iProc) = NINT(GEO%PeriodicVectors(3,iProc)/BGMdeltas(3))
    IF(ABS(GEO%PeriodicVectors(3,iProc)/BGMdeltas(3)-REAL(GEO%PeriodicBGMVectors(3,iProc))).GT.1E-10)THEN
      CALL abort(__STAMP__, &
       'ERROR: Periodic Vector ist not multiple of background mesh delta')
    END IF
  END DO
  ! Check whether process is periodic with itself
  GEO%SelfPeriodic = .FALSE.
  !--- virtually move myself according to periodic vectors in order to find overlapping areas
  !--- 26 possibilities,
  localminmax(1) = BGMminX
  localminmax(2) = BGMminY
  localminmax(3) = BGMminZ
  localminmax(4) = BGMmaxX
  localminmax(5) = BGMmaxY
  localminmax(6) = BGMmaxZ
  DO k = -1,1
    DO m = -1,1
      DO n = -1,1
        PeriodicVec = k*GEO%PeriodicBGMVectors(:,1)
        IF (GEO%nPeriodicVectors.GT.1) THEN
          PeriodicVec = PeriodicVec + m*GEO%PeriodicBGMVectors(:,2)
        END IF
        IF (GEO%nPeriodicVectors.GT.2) THEN
          PeriodicVec = PeriodicVec + n*GEO%PeriodicBGMVectors(:,3)
        END IF
        IF (ALL(PeriodicVec(:).EQ.0)) CYCLE
        periodicminmax(1) = localminmax(1) + PeriodicVec(1)
        periodicminmax(2) = localminmax(2) + PeriodicVec(2)
        periodicminmax(3) = localminmax(3) + PeriodicVec(3)
        periodicminmax(4) = localminmax(4) + PeriodicVec(1)
        periodicminmax(5) = localminmax(5) + PeriodicVec(2)
        periodicminmax(6) = localminmax(6) + PeriodicVec(3)
        !--- find overlap
        DO coord = 1,3  !          x y z direction
          maxofmin = MAX(periodicminmax(coord),localminmax(coord))
          minofmax = MIN(periodicminmax(3+coord),localminmax(3+coord))
          IF (maxofmin.LE.minofmax) GEO%SelfPeriodic = .TRUE.     !  overlapping
        END DO
      END DO
    END DO
  END DO
END IF

! --- send and receive min max indices to and from all processes

!--- enter local min max vector (xmin, ymin, zmin, xmax, ymax, zmax)
localminmax(1) = BGMminX
localminmax(2) = BGMminY
localminmax(3) = BGMminZ
localminmax(4) = BGMmaxX
localminmax(5) = BGMmaxY
localminmax(6) = BGMmaxZ
!--- do allgather into complete min max vector
CALL MPI_ALLGATHER(localminmax,6,MPI_INTEGER,completeminmax,6,MPI_INTEGER,PartMPI%COMM,IERROR)
! Allocate MPIConnect
SDEALLOCATE(PartMPI%DepoBGMConnect)
ALLOCATE(PartMPI%DepoBGMConnect(0:PartMPI%nProcs-1),STAT=allocStat)
IF (allocStat .NE. 0) THEN
  CALL abort(__STAMP__, &
    ' Cannot allocate PartMPI%DepoBGMConnect')
END IF

!--- determine borders indices (=overlapping BGM mesh points) with each process
DO iProc=0,PartMPI%nProcs-1
  PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor = .FALSE.
  PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount = 0
   IF (iProc.EQ.PartMPI%MyRank) THEN
      PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor = .FALSE.
   ELSE
      PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor = .TRUE.
      DO k = 1,3          !  x y z direction
         maxofmin = MAX(localminmax(k),completeminmax((iProc*6)+k))
         minofmax = MIN(localminmax(3+k),completeminmax((iProc*6)+3+k))
         IF (maxofmin.LE.minofmax) THEN          !  overlapping
            TempBorder(1,k) = maxofmin
            TempBorder(2,k) = minofmax
         ELSE
            PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor = .FALSE.
         END IF
      END DO          
   END IF
   IF(PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor)THEN
      SDEALLOCATE(PartMPI%DepoBGMConnect(iProc)%BGMBorder)
      ALLOCATE(PartMPI%DepoBGMConnect(iProc)%BGMBorder(1:2,1:3),STAT=allocStat)
      IF (allocStat .NE. 0) THEN
         CALL abort(__STAMP__, &
          ' Cannot allocate PartMPI%DepoMPIConnect%BGMBorder')
      END IF
      PartMPI%DepoBGMConnect(iProc)%BGMBorder(1:2,1:3) = TempBorder(1:2,1:3)
   END IF
END DO ! iProc=0,PartMPI%nProcs-1

!--- determine border indices for periodic meshes  
IF (GEO%nPeriodicVectors.GT.0) THEN   
  !--- m-/-n-loops must be executed just once (with 0) if respective PV is not present
  !    (otherwise the unnec. 0 are sent and the 0-0-0-case occurs two 2 more which are not cycled, see below for more info)!
  IF (GEO%nPeriodicVectors.GT.1) THEN
    m0=1
  ELSE
    m0=0
  END IF
  IF (GEO%nPeriodicVectors.GT.2) THEN
    n0=1
  ELSE
    n0=0
  END IF
  DO iProc = 0,PartMPI%nProcs-1
    IF (iProc.EQ.PartMPI%MyRank) THEN
      PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor = .FALSE.
      PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount = 0
    ELSE
      !--- virtually move myself according to periodic vectors in order to find overlapping areas
      !--- 26 possibilities, processes need to work through them in opposite direction in order
      !--- to get matching areas.
      !--- Example for 2D:  I am process #3, I compare myself with #7
      !--- Periodic Vectors are p1 and p2.
      !--- I check p1, p2, p1+p2, p1-p2, -p1+p2, -p1-p2, -p2, -p1
      !--- #7 has to check -p1, -p2, -p1-p2, -p1+p2, p1-p2, p1+p1, p2, p1
      !--- This is done by doing 3 loops from -1 to 1 (for the higher process number)
      !--- or 1 to -1 (for the lower process number) and multiplying
      !--- these numbers to the periodic vectors
      NeighCount = 0  ! -- counter: how often is the process my periodic neighbor?
      PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor = .FALSE.
      DO k = -SIGN(1,PartMPI%MyRank-iProc),SIGN(1,PartMPI%MyRank-iProc),SIGN(1,PartMPI%MyRank-iProc)
        DO m = -SIGN(m0,PartMPI%MyRank-iProc),SIGN(m0,PartMPI%MyRank-iProc),SIGN(1,PartMPI%MyRank-iProc)
          DO n = -SIGN(n0,PartMPI%MyRank-iProc),SIGN(n0,PartMPI%MyRank-iProc),SIGN(1,PartMPI%MyRank-iProc)
            IF ((k.EQ.0).AND.(m.EQ.0).AND.(n.EQ.0)) CYCLE ! this is not periodic and already done above
            CHECKNEIGHBOR = .TRUE.
            PeriodicVec = k*GEO%PeriodicBGMVectors(:,1)
            IF (GEO%nPeriodicVectors.GT.1) THEN
              PeriodicVec = PeriodicVec + m*GEO%PeriodicBGMVectors(:,2)
            END IF
            IF (GEO%nPeriodicVectors.GT.2) THEN
              PeriodicVec = PeriodicVec + n*GEO%PeriodicBGMVectors(:,3)
            END IF
            periodicminmax(1) = localminmax(1) + PeriodicVec(1)
            periodicminmax(2) = localminmax(2) + PeriodicVec(2)
            periodicminmax(3) = localminmax(3) + PeriodicVec(3)
            periodicminmax(4) = localminmax(4) + PeriodicVec(1)
            periodicminmax(5) = localminmax(5) + PeriodicVec(2)
            periodicminmax(6) = localminmax(6) + PeriodicVec(3)
            !--- find overlap
            DO coord = 1,3           ! x y z direction
              maxofmin = MAX(periodicminmax(coord),completeminmax((iProc*6)+coord))
              minofmax = MIN(periodicminmax(3+coord),completeminmax((iProc*6)+3+coord))
              IF (maxofmin.LE.minofmax) THEN         !   overlapping
                TempBorder(1,coord) = maxofmin
                TempBorder(2,coord) = minofmax
              ELSE
                CHECKNEIGHBOR = .FALSE.
              END IF
            END DO
            IF(CHECKNEIGHBOR)THEN
              NeighCount = NeighCount + 1
              TempBorder(:,1) = TempBorder(:,1) - PeriodicVec(1)
              TempBorder(:,2) = TempBorder(:,2) - PeriodicVec(2)
              TempBorder(:,3) = TempBorder(:,3) - PeriodicVec(3)
              TempPeriBord(NeighCount,1:2,1:3) = TempBorder(1:2,1:3)
              PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor = .TRUE.
            END IF
          END DO
        END DO
      END DO
      PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount = NeighCount
      ALLOCATE(PartMPI%DepoBGMConnect(iProc)%Periodic(1:PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount),STAT=allocStat)
      IF (allocStat .NE. 0) THEN
        CALL abort(__STAMP__,  &
          'ERROR in MPIBackgroundMeshInit: cannot allocate PartMPI%DepoBGMConnect')
      END IF
      DO k = 1,NeighCount
        ALLOCATE(PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(1:2,1:3),STAT=allocStat)
        IF (allocStat .NE. 0) THEN
          CALL abort(__STAMP__, &
            'ERROR in MPIBackgroundMeshInit: cannot allocate PartMPI%DepoBGMConnect')
        END IF
        PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(1:2,1:3) = TempPeriBord(k,1:2,1:3)
      END DO
    END IF
  END DO
ELSE
  !--- initialize to FALSE for completely non-periodic cases
  DO iProc = 0,PartMPI%nProcs-1
    PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor = .FALSE.
  END DO
END IF

END SUBROUTINE MPIBackgroundMeshInit
#endif /*MPI*/

SUBROUTINE DeBoor(PosInd, aux, coord, results, dir)       
!============================================================================================================================
! recursive function for evaluating a b-spline basis function
!============================================================================================================================
! use MODULES 
   USE MOD_PICDepo_Vars                                                          
!-----------------------------------------------------------------------------------------------------------------------------------
   IMPLICIT NONE                                                                                   
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
   INTEGER, INTENT(IN)                          :: PosInd, dir  
   REAL, INTENT(IN)                             :: coord                                              
   REAL, INTENT(INOUT)                          :: aux(0:3), results                     
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES  
   INTEGER                          :: i,k,jL,jR                                                  
   REAL                              :: hlp1,hlp2                                                   
!-----------------------------------------------------------------------------------------------------------------------------------
    DO i = 0, 2
       DO k = 0, 2-i
          jL = PosInd - k
          jR = jL + 3 - i
          
          hlp1 = jR * BGMdeltas(dir) - coord
          hlp2 = coord - jL * BGMdeltas(dir)
         
          aux(k) = (hlp1 * aux(k+1) + hlp2 * aux(k)) / (hlp1+hlp2)
       ENDDO
    ENDDO
    results = aux(0)

    RETURN
END SUBROUTINE DeBoor


SUBROUTINE DeBoorRef(N_in,NKnots,Knots,Xi_in,UBspline)
!===================================================================================================================================
! DeBoor algorithms for uniform B-Splines
! rule of DeBoor algorithm
! N_i^0 (x) = 1 if x in [u_i, u_i+1); else 0     
! N_i^n (x) = (x-u_i) /(u_i+n-u_i) N_i^n-1(x) + (u_i+n+1 - x)/(u_i+n+1 - u_i+1) N_i+1^n-1(x)
! this algorithm evaluates the complete 1D basis fuction, because certain knots can be reused
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
INTEGER,INTENT(IN)      :: N_in,Nknots
REAL,INTENT(IN)         :: Xi_in
REAL,INTENT(OUT)        :: UBspline(0:N_in)
REAL,INTENT(IN)         :: knots(0:Nknots) ! range of parameter
!REAL,INTENT(IN)         :: DXi is 2 for [-1,1)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: i,n, last!,first
REAL                    :: tmpArray(0:NKnots-1,0:NKnots-1)
REAL                    :: sDxiN!,DxiN1
!===================================================================================================================================

 
! init, first layer (constant)
! zero order
n=0
tmpArray=0.
last=nknots-1
DO i=0,last
  IF((knots(i).LE.Xi_in).AND.(Xi_in.LT.knots(i+1)))THEN
    tmpArray(i,n)=1.0
  END IF ! select
END DO ! i

DO n=1,N_in
  last=last-1
  DO i=0,last
    ! standard
!    tmpArray(i,n) = (Xi_in-knots(i))/(knots(i+n)-knots(i))*tmpArray(i,n-1) &
!                  + (knots(i+n+1)-Xi_in)/(knots(i+n+1)-knots(i+1))*tmpArray(i+1,n-1)
    ! optimized
    sDxiN=0.5/REAL(n)
    tmpArray(i,n) = sDxiN*( (Xi_in-knots(i) )*tmparray(i,n-1)+(knots(i+n+1)-Xi_in)*tmpArray(i+1,n-1))
  END DO ! i
END DO ! n

! move back to correct range
UBSpline(0:N_in)=tmpArray(0:N_in,N_in)

END SUBROUTINE DeBoorRef


FUNCTION beta(z,w)  
!============================================================================================================================
! calculates the beta function
!============================================================================================================================
! use MODULES                                                                                               
   USE nr
!-----------------------------------------------------------------------------------------------------------------------------------
   IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
   REAL, INTENT(IN)                             :: w, z                                             
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
   REAL                                          :: beta
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                    
!----------------------------------------------------------------------------------------------------------------------------------
   beta = exp(gammln(z)+gammln(w)-gammln(z+w))                                                                    
END FUNCTION beta 


SUBROUTINE ComputeGaussDistance(N_In,scaleR,X_in,Elem_xGP,GaussDistance) 
!----------------------------------------------------------------------------------------------------------------------------------!
! compute all distance between X_in and given array 
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES 
INTEGER,INTENT(IN)      :: N_In
REAL,INTENT(IN)         :: X_in(1:3)
REAL,INTENT(IN)         :: scaleR
!REAL,INTENT(IN)         :: Elem_xGP(1:3,1:N_In)
REAL,INTENT(IN)         :: Elem_xGP(1:3,0:N_in,0:N_in,0:N_in)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!REAL,INTENT(OUT)        :: GaussDistance(1:N_In)
REAL,INTENT(OUT)        :: GaussDistance(0:N_in,0:N_in,0:N_in)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: i,j,k
!REAL                    :: tmp(3)
REAL                    :: tmp(3)
!===================================================================================================================================

DO k=0,N_in
  DO j=0,N_in
    DO i=0,N_in
      !tmp=X_in-Elem_xGP(:,i,j,k)
      GaussDistance(i,j,k) =((X_in(1)-Elem_xGP(1,i,j,k))*(X_in(1)-Elem_xGP(1,i,j,k)) &
                            +(X_in(2)-Elem_xGP(2,i,j,k))*(X_in(2)-Elem_xGP(2,i,j,k)) &
                            +(X_in(3)-Elem_xGP(3,i,j,k))*(X_in(3)-Elem_xGP(3,i,j,k)))*scaleR
    END DO ! i=0,N_in
  END DO !  j=0,N_in
END DO ! k=0,N_in


!DO i=1,N_in
!  tmp=X_in-Elem_xGP(:,i)
!  GaussDistance(i) = DOT_PRODUCT(tmp,tmp)*scaleR
!END DO ! i = 1,N_in

END SUBROUTINE ComputeGaussDistance


SUBROUTINE FinalizeDeposition() 
!----------------------------------------------------------------------------------------------------------------------------------!
! finalize pic deposition
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_PICDepo_Vars
USE MOD_Particle_Mesh_Vars,   ONLY:Geo
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE(GaussBorder)
SDEALLOCATE(ElemDepo_xGP)
SDEALLOCATE(Vdm_EquiN_GaussN)
SDEALLOCATE(Knots)
SDEALLOCATE(GaussBGMIndex)
SDEALLOCATE(GaussBGMFactor)
SDEALLOCATE(GEO%PeriodicBGMVectors)
SDEALLOCATE(BGMSource)
SDEALLOCATE(GPWeight)
SDEALLOCATE(ElemRadius2_sf)
SDEALLOCATE(Vdm_NDepo_GaussN)
SDEALLOCATE(sJNDepo)
SDEALLOCATE(XiNDepo)
SDEALLOCATE(swGPNDepo)
SDEALLOCATE(wBaryNDepo)
SDEALLOCATE(NDepochooseK)

END SUBROUTINE FinalizeDeposition

END MODULE MOD_PICDepo
