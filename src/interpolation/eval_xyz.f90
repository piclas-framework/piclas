#include "boltzplatz.h"

MODULE MOD_Eval_xyz
!===================================================================================================================================
! Changes a 3D Tensor Product Lagrange Points of Lagrange Basis of degree N_In to  
! Lagrange points of a Lagrange Basis for one point, using two
! arbitrary point disributions xi_In(0:N_In) and xi_Out 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
  INTEGER :: errorflag
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------

#ifdef PARTICLES
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE eval_xyz_curved
  MODULE PROCEDURE eval_xyz_curved
END INTERFACE

INTERFACE eval_xyz_elemcheck
  MODULE PROCEDURE eval_xyz_elemcheck
END INTERFACE

INTERFACE eval_xyz_part2
  MODULE PROCEDURE eval_xyz_part2
END INTERFACE

PUBLIC :: eval_xyz_curved,eval_xyz_elemcheck, eval_xyz_part2
#endif /*PARTICLES*/
!===================================================================================================================================

CONTAINS

#ifdef PARTICLES
SUBROUTINE eval_xyz_curved(x_in,NVar,N_in,X3D_In,X3D_Out,iElem,PartID)
!===================================================================================================================================
! interpolate a 3D tensor product Lagrange basis defined by (N_in+1) 1D interpolation point positions x
! first get xi,eta,zeta from x,y,z...then do tenso product interpolation
! xi is defined in the 1DrefElem xi=[-1,1]
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Basis,                   ONLY:LagrangeInterpolationPolys
USE MOD_Interpolation_Vars,      ONLY:wBary,xGP
USE MOD_Mesh_Vars,               ONLY:dXCL_NGeo,Elem_xGP,XCL_NGeo,NGeo,wBaryCL_NGeo,XiCL_NGeo
USE MOD_Particle_Surfaces_Vars,  ONLY:epsilonOne,MappingGuess,epsMapping
USE MOD_Particle_surfaces_Vars,  ONLY:XiEtaZetaBasis,ElemBaryNGeo,slenXiEtaZetaBasis
USE MOD_PICInterpolation_Vars,   ONLY:NBG,BGField,useBGField,BGDataSize,BGField_wBary, BGField_xGP,BGType
USE MOD_Particle_MPI_Vars,       ONLY:PartMPI
USE MOD_Particle_Surfaces_Vars,  ONLY:ElemRadiusNGeo
!USE MOD_Mesh_Vars,ONLY: X_CP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: NVar                                  ! 6 (Ex, Ey, Ez, Bx, By, Bz) 
INTEGER,INTENT(IN)  :: N_In                                  ! usually PP_N
INTEGER,INTENT(IN)  :: iElem                                 ! elem index
REAL,INTENT(IN)     :: X3D_In(1:NVar,0:N_In,0:N_In,0:N_In)   ! elem state
REAL,INTENT(IN)     :: x_in(3)                                  ! physical position of particle 
INTEGER,INTENT(IN),OPTIONAL :: PartID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: X3D_Out(1:NVar)  ! Interpolated state
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER             :: i,j,k
REAL                :: xi(3)
INTEGER             :: NewTonIter
REAL                :: X3D_Buf1(1:NVar,0:N_In,0:N_In)  ! first intermediate results from 1D interpolations
REAL                :: X3D_Buf2(1:NVar,0:N_In) ! second intermediate results from 1D interpolations
REAL                :: Winner_Dist,Dist,abortcrit
REAL, PARAMETER     :: EPSONE=1.00000001
INTEGER             :: iDir
REAL                :: F(1:3),Lag(1:3,0:NGeo),Lag2(1:3,0:N_In)
REAL                :: Jac(1:3,1:3),sdetJac,sJac(1:3,1:3)
REAL                :: buff,buff2
REAL                :: Ptild(1:3),XiLinear(1:6)
! h5-external e,b field
REAL,ALLOCATABLE    :: L_xi_BGField(:,:), X3D_tmp1(:,:,:), X3D_tmp2(:,:), X3D_tmp3(:)
!===================================================================================================================================

!print*,'iElem',iElem
!print*,'Pos',X_in

! get initial guess by nearest GP search ! simple guess
SELECT CASE(MappingGuess)
CASE(1)
  Ptild=X_in - ElemBaryNGeo(:,iElem)
  ! plus coord system (1-3) and minus coord system (4-6)
  DO iDir=1,6
    XiLinear(iDir)=DOT_PRODUCT(Ptild,XiEtaZetaBasis(:,iDir,iElem))*slenXiEtaZetaBasis(iDir,iElem)
  END DO
  ! compute guess as average value
  DO iDir=1,3
    Xi(iDir)=0.5*(XiLinear(iDir)-XiLinear(iDir+3))
  END DO 
  IF(MAXVAL(ABS(Xi)).GT.epsOne) Xi=0.
!  IF(ANY(ABS(Xi).GT.epsOne)) THEN
!    DO iDir=1,3
!      IF(Xi(iDir).GT.epsOne) Xi(iDir)=1.0
!      IF(Xi(iDir).LT.-epsOne) Xi(iDir)=-1.0
!    END DO ! iDir
!  END IF
CASE(2) 
  ! compute distance on Gauss Points
  Winner_Dist=HUGE(1.)
  DO i=0,N_in; DO j=0,N_in; DO k=0,N_in
    Dist=SUM((x_in(:)-Elem_xGP(:,i,j,k,iElem))*(x_in(:)-Elem_xGP(:,i,j,k,iElem)))
    IF (Dist.LT.Winner_Dist) THEN
      Winner_Dist=Dist
      Xi(:)=(/xGP(i),xGP(j),xGP(k)/) ! start value
    END IF
  END DO; END DO; END DO
CASE(3) 
  ! compute distance on XCL Points
  Winner_Dist=HUGE(1.)
  DO i=0,NGeo; DO j=0,NGeo; DO k=0,NGeo
    Dist=SUM((x_in(:)-XCL_NGeo(:,i,j,k,iElem))*(x_in(:)-XCL_NGeo(:,i,j,k,iElem)))
    IF (Dist.LT.Winner_Dist) THEN
      Winner_Dist=Dist
      Xi(:)=(/XiCL_NGeo(i),XiCL_NGeo(j),XiCL_NGeo(k)/) ! start value
    END IF
  END DO; END DO; END DO
CASE(4)
  ! trival guess 
  xi=0.
END SELECT
!print*,'Winnerdist',Winner_Dist
!print*,'initial guess'
!print*,'xi',xi
!
!print*,'NGeo',NGeo
!print*,'Xi_CLNGeo',XiCL_NGeo
!print*,'wbary_CLNGeo',wBaryCL_NGeo
!read*

!DO k=0,NGeo
!  DO j=0,NGeo
!    DO i=0,NGeo
!     !! Matrix-vector multiplication
!       print*,'dXCL_NGeo',dXCL_NGeo(:,1,i,j,k,iElem)
!       print*,'dXCL_NGeo',dXCL_NGeo(:,2,i,j,k,iElem)
!       print*,'dXCL_NGeo',dXCL_NGeo(:,3,i,j,k,iElem)
!       read*
!    END DO !i=0,NGeo
!  END DO !j=0,NGeo
!END DO !k=0,NGeo


! initial guess
CALL LagrangeInterpolationPolys(Xi(1),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
CALL LagrangeInterpolationPolys(Xi(2),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
CALL LagrangeInterpolationPolys(Xi(3),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))
! F(xi) = x(xi) - x_in
F=-x_in ! xRp
DO k=0,NGeo
  DO j=0,NGeo
    DO i=0,NGeo
      F=F+XCL_NGeo(:,i,j,k,iElem)*Lag(1,i)*Lag(2,j)*Lag(3,k)
      !print*,'XCL',i,j,k,XCL_NGeo(:,i,j,k,iElem)
    END DO !l=0,NGeo
  END DO !i=0,NGeo
END DO !j=0,NGeo
!print*,'F',F
!read*

!DO WHILE ((SUM(F*F).GT.epsMapping).AND.(NewtonIter.LT.100))
!  NewtonIter=NewtonIter+1
NewtonIter=0
abortCrit=ElemRadiusNGeo(iElem)*epsMapping
!DO WHILE ((SUM(F*F).GT.abortCrit).AND.(NewtonIter.LT.50))
DO WHILE ((SUM(F*F).GT.abortCrit).AND.(NewtonIter.LT.100))
  NewtonIter=NewtonIter+1

  ! 
  ! caution, dXCL_NGeo is transposed of required matrix
  Jac=0.
  DO k=0,NGeo
    DO j=0,NGeo
      buff=Lag(2,j)*Lag(3,k)
      DO i=0,NGeo
        buff2=Lag(1,i)*buff
        Jac(1,1:3)=Jac(1,1:3)+dXCL_NGeo(1:3,1,i,j,k,iElem)*buff2
        Jac(2,1:3)=Jac(2,1:3)+dXCL_NGeo(1:3,2,i,j,k,iElem)*buff2
        Jac(3,1:3)=Jac(3,1:3)+dXCL_NGeo(1:3,3,i,j,k,iElem)*buff2
      END DO !i=0,NGeo
    END DO !j=0,NGeo
  END DO !k=0,NGeo
  
  !print*,'print',Jac

  ! Compute inverse of Jacobian
  sdetJac=getDet(Jac)
  IF(sdetJac.GT.0.) THEN
   sdetJac=1./sdetJac
  ELSE !shit
   ! Newton has not converged !?!?
   CALL abort(__STAMP__, &
        'Newton in FindXiForPartPos singular. iter,sdetJac',NewtonIter,sDetJac)
  ENDIF 
  sJac=getInv(Jac,sdetJac)
  
  ! Iterate Xi using Newton step
  ! Use FAIL
  Xi = Xi - MATMUL(sJac,F)
  IF(ANY(ABS(Xi).GT.1.2)) THEN
    IPWRITE(*,*) ' Particle not inside of element, force!!!'
    IPWRITE(*,*) ' xi  ', xi(1)
    IPWRITE(*,*) ' eta ', xi(2)
    IPWRITE(*,*) ' zeta', xi(3)
    IPWRITE(*,*) ' PartPos', X_in
    IPWRITE(*,*) ' ElemID', iElem
    IF(PRESENT(PartID)) IPWRITE(*,*) ' PartID', PartID
    CALL abort(__STAMP__, &
        'Particle Not inSide of Element, iProc, iElem',PartMPI%MyRank,REAL(iElem))
  END IF
  
  ! Compute function value
  CALL LagrangeInterpolationPolys(Xi(1),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
  CALL LagrangeInterpolationPolys(Xi(2),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
  CALL LagrangeInterpolationPolys(Xi(3),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))
  ! F(xi) = x(xi) - x_in
  F=-x_in ! xRp
  DO k=0,NGeo
    DO j=0,NGeo
      buff=Lag(2,j)*Lag(3,k)
      DO i=0,NGeo
        buff2=Lag(1,i)*buff
        F=F+XCL_NGeo(:,i,j,k,iElem)*buff2
      END DO !l=0,NGeo
    END DO !i=0,NGeo
  END DO !j=0,NGeo
END DO !newton

! check if Newton is successful
!IF(ANY(ABS(Xi).GT.epsilonOne)) THEN
!IF(ANY(ABS(Xi).GT.1.0)) THEN
!  IF(PRESENT(PartID)) WRITE(*,*) 'ParticleID', PartID
! ! IF(PRESENT(PartID).AND.PartID.EQ.40) THEN
!!  !   WRITE(*,*) 'ParticleID', PartID
!     IPWRITE(*,*) ' elemid', ielem
!     IPWRITE(*,*) ' Particle outside of parameter range!!!'
!     IPWRITE(*,*) ' xi  ', xi(:)
!!     WRITE(*,*) ' eta ', xi(2)
!!     WRITE(*,*) ' zeta', xi(3)
! ! END IF
!END IF

! 2.1) get "Vandermonde" vectors
DO i=1,3
  CALL LagrangeInterpolationPolys(xi(i),N_in,xGP,wBary,Lag2(i,:))
 !CALL LagrangeInterpolationPolys(xi(i),N_in,xGP,wBary,L_xi(i,:))
END DO

! 2.2) do the tensor product thing 
X3D_buf1=0.
! first direction iN_In
DO k=0,N_In
  DO j=0,N_In
    DO i=0,N_In
      X3D_Buf1(:,j,k)=X3D_Buf1(:,j,k)+Lag2(1,i)*X3D_In(:,i,j,k)
    END DO
  END DO
END DO
X3D_buf2=0.
! second direction jN_In
DO k=0,N_In
  DO j=0,N_In
    X3D_Buf2(:,k)=X3D_Buf2(:,k)+Lag2(2,j)*X3D_Buf1(:,j,k)
  END DO
END DO
X3D_Out=0.
! last direction kN_In
DO k=0,N_In
  X3D_Out(:)=X3D_Out(:)+Lag2(3,k)*X3D_Buf2(:,k)
END DO

IF(useBGField)THEN
  ! use of BG-Field with possible different polynomial order and nodetype
  ALLOCATE( L_xi_BGField(3,0:NBG)            &
          , X3D_tmp1(BGDataSize,0:NBG,0:NBG) &
          , X3D_tmp2(BGDataSize,0:NBG)       &
          , X3D_tmp3(BGDataSize)             )
  DO i=1,3
    CALL LagrangeInterpolationPolys(xi(i),NBG,BGField_xGP,BGField_wBary,L_xi_BGField(i,:))
  END DO
  
  ! 2.2) do the tensor product thing 
  X3D_tmp1=0.
  ! first direction iN_In
  DO k=0,NBG
    DO j=0,NBG
      DO i=0,NBG
        X3D_tmp1(:,j,k)=X3D_tmp1(:,j,k)+L_xi_BGField(1,i)*BGField(:,i,j,k,iElem)
      END DO
    END DO
  END DO
  X3D_tmp2=0.
  ! second direction jN_In
  DO k=0,NBG
    DO j=0,NBG
      X3D_tmp2(:,k)=X3D_tmp2(:,k)+L_xi_BGField(2,j)*X3D_tmp1(:,j,k)
    END DO
  END DO
  X3D_tmp3=0.
  ! last direction kN_In
  DO k=0,NBG
    X3D_tmp3(:)=X3D_tmp3(:)+L_xi_BGField(3,k)*X3D_tmp2(:,k)
  END DO
  SELECT CASE(BGType)
  CASE(1)
    X3d_Out(1:3)=x3d_Out(1:3)+X3d_tmp3
  CASE(2)
    X3d_Out(4:6)=x3d_Out(4:6)+X3d_tmp3
  CASE(3)
    X3d_Out=x3d_Out+X3d_tmp3
  END SELECT
  DEALLOCATE( L_xi_BGField, X3d_tmp1, x3d_tmp2, x3d_tmp3)
END IF ! useBGField

END SUBROUTINE eval_xyz_curved


SUBROUTINE eval_xyz_elemcheck(x_in,xi,iElem,PartID,DoReUseMap)
!===================================================================================================================================
! interpolate a 3D tensor product Lagrange basis defined by (N_in+1) 1D interpolation point positions x
! first get xi,eta,zeta from x,y,z...then do tenso product interpolation
! xi is defined in the 1DrefElem xi=[-1,1]
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Basis,                   ONLY:LagrangeInterpolationPolys
USE MOD_Interpolation_Vars,      ONLY:wBary,xGP
USE MOD_Particle_Surfaces_Vars,  ONLY:MappingGuess,epsMapping
USE MOD_Particle_surfaces_Vars,  ONLY:XiEtaZetaBasis,ElemBaryNGeo,slenXiEtaZetaBasis
USE MOD_Mesh_Vars,               ONLY:dXCL_NGeo,Elem_xGP,XCL_NGeo,NGeo,wBaryCL_NGeo,XiCL_NGeo,NGeo
USE MOD_Particle_Surfaces_Vars,  ONLY:BezierControlPoints3D,SideType
USE MOD_Mesh_Vars,               ONLY:ElemToSide
USE MOD_Particle_Surfaces_Vars,  ONLY:ElemRadiusNGeo
!USE MOD_Mesh_Vars,ONLY: X_CP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)          :: iElem                                 ! elem index
REAL,INTENT(IN)             :: x_in(3)                                  ! physical position of particle 
INTEGER,INTENT(IN),OPTIONAL :: PartID
LOGICAL,INTENT(IN),OPTIONAL :: DoReUseMap
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)          :: xi(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                    :: i,j,k
REAL                       :: epsOne
INTEGER                    :: NewTonIter
REAL                       :: Winner_Dist,Dist
INTEGER                    :: n_Newton,idir
REAL                       :: F(1:3),Lag(1:3,0:NGeo)
REAL                       :: Jac(1:3,1:3),sdetJac,sJac(1:3,1:3)
REAL                       :: buff,buff2,abortcrit
REAL                       :: Ptild(1:3),XiLinear(1:6)
!===================================================================================================================================

epsOne=1.0+epsMapping
IF(.NOT.PRESENT(DoReUseMap))THEN
  SELECT CASE(MappingGuess)
  CASE(1)
    Ptild=X_in - ElemBaryNGeo(:,iElem)
    ! plus coord system (1-3) and minus coord system (4-6)
    DO iDir=1,6
      XiLinear(iDir)=DOT_PRODUCT(Ptild,XiEtaZetaBasis(:,iDir,iElem))*slenXiEtaZetaBasis(iDir,iElem)
    END DO
    ! compute guess as average value
    DO iDir=1,3
      Xi(iDir)=0.5*(XiLinear(iDir)-XiLinear(iDir+3))
    END DO 
    IF(MAXVAL(ABS(Xi)).GT.epsOne) Xi=0.
    !IF(ANY(ABS(Xi).GT.epsOne)) THEN
    !  DO iDir=1,3
    !    IF(Xi(iDir).GT.epsOne) Xi(iDir)=1.0
    !    IF(Xi(iDir).LT.-epsOne) Xi(iDir)=-1.0
    !  END DO ! iDir
    !END IF
  CASE(2)
    Winner_Dist=HUGE(1.)
    DO i=0,PP_N; DO j=0,PP_N; DO k=0,PP_N
      Dist=SUM((x_in(:)-Elem_xGP(:,i,j,k,iElem))*(x_in(:)-Elem_xGP(:,i,j,k,iElem)))
      IF (Dist.LT.Winner_Dist) THEN
        Winner_Dist=Dist
        Xi(:)=(/xGP(i),xGP(j),xGP(k)/) ! start value
      END IF
    END DO; END DO; END DO
  CASE(3) 
    ! compute distance on XCL Points
    Winner_Dist=HUGE(1.)
    DO i=0,NGeo; DO j=0,NGeo; DO k=0,NGeo
      Dist=SUM((x_in(:)-XCL_NGeo(:,i,j,k,iElem))*(x_in(:)-XCL_NGeo(:,i,j,k,iElem)))
      IF (Dist.LT.Winner_Dist) THEN
        Winner_Dist=Dist
        Xi(:)=(/XiCL_NGeo(i),XiCL_NGeo(j),XiCL_NGeo(k)/) ! start value
      END IF
    END DO; END DO; END DO
  CASE(4)
    ! trival guess, cell mean point
    xi=0.
  END SELECT
END IF

! initial guess
CALL LagrangeInterpolationPolys(Xi(1),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
CALL LagrangeInterpolationPolys(Xi(2),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
CALL LagrangeInterpolationPolys(Xi(3),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))
! F(xi) = x(xi) - x_in
F=-x_in ! xRp
DO k=0,NGeo
  DO j=0,NGeo
    DO i=0,NGeo
      F=F+XCL_NGeo(:,i,j,k,iElem)*Lag(1,i)*Lag(2,j)*Lag(3,k)
    END DO !i=0,NGeo
  END DO !j=0,NGeo
END DO !k=0,NGeo

NewtonIter=0
abortCrit=ElemRadiusNGeo(iElem)*epsMapping
!DO WHILE ((SUM(F*F).GT.abortCrit).AND.(NewtonIter.LT.50))
DO WHILE ((SUM(F*F).GT.abortCrit).AND.(NewtonIter.LT.100))
  NewtonIter=NewtonIter+1
  ! 
  ! caution, dXCL_NGeo is transposed of required matrix
  Jac=0.
  DO k=0,NGeo
    DO j=0,NGeo
      buff=Lag(2,j)*Lag(3,k)
      DO i=0,NGeo
        buff2=Lag(1,i)*buff
        Jac(1,1:3)=Jac(1,1:3)+dXCL_NGeo(1:3,1,i,j,k,iElem)*buff2
        Jac(2,1:3)=Jac(2,1:3)+dXCL_NGeo(1:3,2,i,j,k,iElem)*buff2
        Jac(3,1:3)=Jac(3,1:3)+dXCL_NGeo(1:3,3,i,j,k,iElem)*buff2
      END DO !i=0,NGeo
    END DO !j=0,NGeo
  END DO !k=0,NGeo
  
  !print*,'print',Jac

  ! Compute inverse of Jacobian
  sdetJac=getDet(Jac)
  IF(sdetJac.GT.0.) THEN
   sdetJac=1./sdetJac
  ELSE !shit
   ! Newton has not converged !?!?
   CALL abort(__STAMP__, &
        'Newton in FindXiForPartPos singular. iter,sdetJac',NewtonIter,sDetJac)
  ENDIF 
  sJac=getInv(Jac,sdetJac)
  
  ! Iterate Xi using Newton step
  ! Use FAIL
  Xi = Xi - MATMUL(sJac,F)
  IF(ANY(ABS(Xi).GT.1.5)) THEN
    IF(PRESENT(PartID)) THEN
      IPWRITE(*,*) 'ParticleID', PartID
      IPWRITE(*,*) ' Particle not inside of element!!!'
      IPWRITE(*,*) ' Element', iElem
      IPWRITE(*,*) ' xi  ', xi(1)
      IPWRITE(*,*) ' eta ', xi(2)
      IPWRITE(*,*) ' zeta', xi(3)
      IPWRITE(*,*) ' PartPos', X_in
      !CALL abort(__STAMP__, &
      !    'Particle Not inSide of Element')
    ELSE
      EXIT
    END IF
  END IF
  
  ! Compute function value
  CALL LagrangeInterpolationPolys(Xi(1),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(1,:))
  CALL LagrangeInterpolationPolys(Xi(2),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(2,:))
  CALL LagrangeInterpolationPolys(Xi(3),NGeo,XiCL_NGeo,wBaryCL_NGeo,Lag(3,:))
  ! F(xi) = x(xi) - x_in
  F=-x_in ! xRp
  DO k=0,NGeo
    DO j=0,NGeo
      buff=Lag(2,j)*Lag(3,k)
      DO i=0,NGeo
        buff2=Lag(1,i)*buff
        F=F+XCL_NGeo(:,i,j,k,iElem)*buff2
        !F=F+XCL_NGeo(:,i,j,k,iElem)*Lag(1,i)*Lag(2,j)*Lag(3,k)
      END DO !l=0,NGeo
    END DO !i=0,NGeo
  END DO !j=0,NGeo
END DO !newton


!IF(MAXVAL(ABS(Xi)).GT.1.0) THEN
!!  IF(PRESENT(PartID).AND.PartID.EQ.40) THEN
!  IF(PRESENT(PartID)) WRITE(*,*) 'ParticleID', PartID
!    IPWRITE(*,*) ' Particle not inside of element!!!'
!    IPWRITE(*,*) ' Element', iElem
!    IPWRITE(*,*) ' xi  ', xi(:)
!  !  IPWRITE(*,*) ' PartPos', X_in
!!  END IF
!END IF
!

END SUBROUTINE eval_xyz_elemcheck


SUBROUTINE eval_xyz_part2(xi_in,NVar,N_in,X3D_In,X3D_Out,iElem)
!===================================================================================================================================
! interpolate a 3D tensor product Lagrange basis defined by (N_in+1) 1D interpolation point positions x
! hoewver, particle is already mapped to reference space -1|1
! xi instead of x_in
!===================================================================================================================================
! MODULES
USE MOD_Basis,                 ONLY: LagrangeInterpolationPolys
USE MOD_Interpolation_Vars,    ONLY: wBary,xGP
USE MOD_PICInterpolation_Vars, ONLY:NBG,BGField,useBGField,BGDataSize,BGField_xGP,BGField_wBary,BGType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)        :: NVar                                  ! 6 (Ex, Ey, Ez, Bx, By, Bz) 
INTEGER,INTENT(IN)        :: N_In                                  ! usually PP_N
INTEGER,INTENT(IN)        :: iElem                                 ! elem index
REAL,INTENT(IN)           :: X3D_In(1:NVar,0:N_In,0:N_In,0:N_In)   ! elem state
REAL,INTENT(IN)           :: xi_in(3)                              ! reference space position of particle 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)          :: X3D_Out(1:NVar)  ! Interpolated state
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                   :: iN_In,jN_In,kN_In,i
REAL                      :: X3D_Buf1(1:NVar,0:N_In,0:N_In)  ! first intermediate results from 1D interpolations
REAL                      :: X3D_Buf2(1:NVar,0:N_In) ! second intermediate results from 1D interpolations
REAL,DIMENSION(3,0:N_in)  :: L_xi        
REAL,ALLOCATABLE          :: L_xi_BGField(:,:), X3D_tmp1(:,:,:), X3D_tmp2(:,:), X3D_tmp3(:)
!===================================================================================================================================
! --------------------------------------------------
! 2.) Interpolation
! --------------------------------------------------
! 2.1) get "Vandermonde" vectors
DO i=1,3
  CALL LagrangeInterpolationPolys(xi_in(i),N_in,xGP,wBary,L_xi(i,:))
END DO

! 2.2) do the tensor product thing 
X3D_buf1=0.
! first direction iN_In
DO kN_In=0,N_In
  DO jN_In=0,N_In
    DO iN_In=0,N_In
      X3D_Buf1(:,jN_In,kN_In)=X3D_Buf1(:,jN_In,kN_In)+L_xi(1,iN_In)*X3D_In(:,iN_In,jN_In,kN_In)
    END DO
  END DO
END DO
X3D_buf2=0.
! second direction jN_In
DO kN_In=0,N_In
  DO jN_In=0,N_In
    X3D_Buf2(:,kN_In)=X3D_Buf2(:,kN_In)+L_xi(2,jN_In)*X3D_Buf1(:,jN_In,kN_In)
  END DO
END DO
X3D_Out=0.
! last direction kN_In
DO kN_In=0,N_In
  X3D_Out(:)=X3D_Out(:)+L_xi(3,kN_In)*X3D_Buf2(:,kN_In)
END DO

IF(useBGField)THEN
  ! use of BG-Field with possible different polynomial order and nodetype
  ALLOCATE( L_xi_BGField(3,0:NBG)            &
          , X3D_tmp1(BGDataSize,0:NBG,0:NBG) &
          , X3D_tmp2(BGDataSize,0:NBG)       &
          , X3D_tmp3(BGDataSize)             )
  DO i=1,3
    CALL LagrangeInterpolationPolys(xi_in(i),NBG,BGField_xGP,BGField_wBary,L_xi_BGField(i,:))
  END DO
  
  ! 2.2) do the tensor product thing 
  X3D_tmp1=0.
  ! first direction iN_In
  DO kN_In=0,NBG
    DO jN_In=0,NBG
      DO iN_In=0,NBG
        X3D_tmp1(:,jN_In,kN_In)=X3D_tmp1(:,jN_In,kN_In)+L_xi_BGField(1,iN_In)*BGField(:,iN_In,jN_In,kN_In,iElem)
      END DO
    END DO
  END DO
  X3D_tmp2=0.
  ! second direction jN_In
  DO kN_In=0,NBG
    DO jN_In=0,NBG
      X3D_tmp2(:,kN_In)=X3D_tmp2(:,kN_In)+L_xi_BGField(2,jN_In)*X3D_tmp1(:,jN_In,kN_In)
    END DO
  END DO
  X3D_tmp3=0.
  ! last direction kN_In
  DO kN_In=0,NBG
    X3D_tmp3(:)=X3D_tmp3(:)+L_xi_BGField(3,kN_In)*X3D_tmp2(:,kN_In)
  END DO
  SELECT CASE(BGType)
  CASE(1)
    X3d_Out(1:3)=x3d_Out(1:3)+X3d_tmp3
  CASE(2)
    X3d_Out(4:6)=x3d_Out(4:6)+X3d_tmp3
  CASE(3)
    X3d_Out=x3d_Out+X3d_tmp3
  END SELECT

  DEALLOCATE( L_xi_BGField, X3d_tmp1, x3d_tmp2, x3d_tmp3)
END IF ! useBGField
END SUBROUTINE eval_xyz_part2


FUNCTION getDet(Mat)
!=================================================================================================================================
! compute determinant of 3x3 matrix
!=================================================================================================================================
  ! MODULES
  ! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: Mat(3,3)
!---------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL             :: getDet
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!=================================================================================================================================
getDet=   ( Mat(1,1) * Mat(2,2) - Mat(1,2) * Mat(2,1) ) * Mat(3,3) &
        + ( Mat(1,2) * Mat(2,3) - Mat(1,3) * Mat(2,2) ) * Mat(3,1) &
        + ( Mat(1,3) * Mat(2,1) - Mat(1,1) * Mat(2,3) ) * Mat(3,2)
END FUNCTION getDet


FUNCTION getInv(Mat,sdet)
!=================================================================================================================================
! compute inverse of 3x3 matrix, needs sDet=1/det(Mat)
!=================================================================================================================================
  ! MODULES
  ! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: Mat(3,3),sDet
!---------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL             :: getInv(3,3)
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!=================================================================================================================================
getInv(1,1) = ( Mat(2,2) * Mat(3,3) - Mat(2,3) * Mat(3,2) ) * sdet
getInv(1,2) = ( Mat(1,3) * Mat(3,2) - Mat(1,2) * Mat(3,3) ) * sdet
getInv(1,3) = ( Mat(1,2) * Mat(2,3) - Mat(1,3) * Mat(2,2) ) * sdet
getInv(2,1) = ( Mat(2,3) * Mat(3,1) - Mat(2,1) * Mat(3,3) ) * sdet
getInv(2,2) = ( Mat(1,1) * Mat(3,3) - Mat(1,3) * Mat(3,1) ) * sdet
getInv(2,3) = ( Mat(1,3) * Mat(2,1) - Mat(1,1) * Mat(2,3) ) * sdet
getInv(3,1) = ( Mat(2,1) * Mat(3,2) - Mat(2,2) * Mat(3,1) ) * sdet
getInv(3,2) = ( Mat(1,2) * Mat(3,1) - Mat(1,1) * Mat(3,2) ) * sdet
getInv(3,3) = ( Mat(1,1) * Mat(2,2) - Mat(1,2) * Mat(2,1) ) * sdet
END FUNCTION getInv 
#endif /*PARTICLES*/

END MODULE MOD_Eval_xyz
