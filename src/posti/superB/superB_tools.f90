!==================================================================================================================================
! Copyright (c) 2019 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (piclas.boltzplatz.eu/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

MODULE MOD_SuperB_Tools
!===================================================================================================================================
!> Contains the routines and algorithms required for the calculation of magnetic fields with SuperB
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: FindLinIndependentVectors, GramSchmidtAlgo
!===================================================================================================================================
INTERFACE FindLinIndependentVectors
  MODULE PROCEDURE FindLinIndependentVectors
END INTERFACE FindLinIndependentVectors

INTERFACE GramSchmidtAlgo
  MODULE PROCEDURE GramSchmidtAlgo
END INTERFACE GramSchmidtAlgo

INTERFACE CalcErrorSuperB
  MODULE PROCEDURE CalcErrorSuperB
END INTERFACE
!===================================================================================================================================
PUBLIC::CalcErrorSuperB

CONTAINS

SUBROUTINE FindLinIndependentVectors(NormalVector, Vector1, Vector2)
!===================================================================================================================================
!> Finds two linear vectors of a normal vector around a base point
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN) :: NormalVector(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT) :: Vector1(3), Vector2(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! Find the second vector which is in the normal plane
IF (NormalVector(1).NE.0) THEN
  Vector1(1) = (0 - NormalVector(2) - NormalVector(3)) / NormalVector(1)
  Vector1(2) = 1
  Vector1(3) = 1
ELSE IF (NormalVector(2).NE.0) THEN
  Vector1(1) = 1
  Vector1(2) = (0 - NormalVector(1) - NormalVector(3)) / NormalVector(2)
  Vector1(3) = 1
ELSE IF (NormalVector(3).NE.0) THEN
  Vector1(1) = 1
  Vector1(2) = 1
  Vector1(3) = (0 - NormalVector(1) - NormalVector(2)) / NormalVector(3)
ELSE
  CALL abort(__STAMP__&
      ,'The normal direction vector can not be (0,0,0)')
END IF

! Find the third vecord vector with the cross product
Vector2(1) = NormalVector(2)*Vector1(3) - NormalVector(3)*Vector1(2)
Vector2(2) = NormalVector(3)*Vector1(1) - NormalVector(1)*Vector1(3)
Vector2(3) = NormalVector(1)*Vector1(2) - NormalVector(2)*Vector1(1)

END SUBROUTINE FindLinIndependentVectors


PPURE SUBROUTINE GramSchmidtAlgo(Vector1, Vector2, Vector3)
!===================================================================================================================================
!> Contains the Gram Schmidt algorithm for an orthonormal basis
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(INOUT) :: Vector1(3), Vector2(3), Vector3(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

! v1 = w1/||w1||
Vector1(:) = Vector1(:) / SQRT(Vector1(1)**2 + Vector1(2)**2 + Vector1(3)**2)

! v2 = w2 - <v1,w2>*v1
Vector2(:) = Vector2(:) - DOT_PRODUCT(Vector1, Vector2) * Vector1(:)
! v2 = v2/||v2||
Vector2(:) = Vector2(:) / SQRT(Vector2(1)**2 + Vector2(2)**2 + Vector2(3)**2)

! v3 = w3 - <v1,w3>*v1 - <v2,w3>*v2
Vector3(:) = Vector3(:) - DOT_PRODUCT(Vector1, Vector3) * Vector1(:) -&
                          DOT_PRODUCT(Vector2, Vector3) * Vector2(:)
! v3 = v3/||v3||
Vector3(:) = Vector3(:) / SQRT(Vector3(1)**2 + Vector3(2)**2 + Vector3(3)**2)

END SUBROUTINE GramSchmidtAlgo


SUBROUTINE CalcErrorSuperB(L_2_Error,L_Inf_Error,ExactFunctionNumber,iCoilOrMagnet)
!===================================================================================================================================
! Calculates L_infinfity and L_2 norms of state variables using the Analyze Framework (GL points+weights)
! The analyze polynomial does NOT use Gauss nodes, because this would artificially reduce the error (because of the nodal
! interpolation). Therefore, Gauss-Lobatto (GL) nodes are used (the numerical solution is mapped onto GL nodes) because of the better
! integration quality as compared with equidistant nodes
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars          ,ONLY: Elem_xGP,sJ,Vdm_GL_N
USE MOD_Interpolation_Vars ,ONLY: NAnalyze,Vdm_GaussN_NAnalyze,wAnalyze
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_Particle_Mesh_Vars ,ONLY: MeshVolume
USE MOD_Interpolation_Vars ,ONLY: BGField,BGFieldAnalytic
USE MOD_Interpolation      ,ONLY: GetVandermonde
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)   :: ExactFunctionNumber    ! Number of exact function to be used for the calculation of the analytical solution
INTEGER,INTENT(IN)   :: iCoilOrMagnet          ! Number of Coil or Magnet
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)     :: L_2_Error(4)   !< L2 error of the solution
REAL,INTENT(OUT)     :: L_Inf_Error(4) !< LInf error of the solution
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iElem,k,l,m
REAL                 :: U_exact(3)
REAL                 :: U_NAnalyze(1:3,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                 :: U_NAnalyze_tmp(1:3,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                 :: Coords_NAnalyze(3,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                 :: J_NAnalyze(1,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                 :: J_N(1,0:PP_N,0:PP_N,0:PP_N)
REAL                 :: IntegrationWeight
!===================================================================================================================================
! Initialize errors
L_Inf_Error(:)=-1.E10
L_2_Error(:)=0.

! Interpolate values of Error-Grid from GP's
DO iElem=1,PP_nElems
   ! Interpolate the physical position Elem_xGP to the analyze position, needed for exact function
   CALL ChangeBasis3D(3,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,Elem_xGP(1:3,:,:,:,iElem),Coords_NAnalyze(1:3,:,:,:))
   ! Interpolate the Jacobian to the analyze grid: be carefull we interpolate the inverse of the inverse of the jacobian ;-)
   J_N(1,0:PP_N,0:PP_N,0:PP_N)=1./sJ(:,:,:,iElem)
   CALL ChangeBasis3D(1,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,J_N(1:1,0:PP_N,0:PP_N,0:PP_N),J_NAnalyze(1:1,:,:,:))
   ! Interpolate the solution to the analyze grid
   CALL ChangeBasis3D(3,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,BGField(1:3,:,:,:,iElem),U_NAnalyze(1:3,:,:,:))
   DO m=0,NAnalyze
     DO l=0,NAnalyze
       DO k=0,NAnalyze
         CALL ExactFuncSuperB(ExactFunctionNumber,iCoilOrMagnet,Coords_NAnalyze(1:3,k,l,m),U_exact)
         L_Inf_Error(1:3) = MAX(L_Inf_Error(1:3),abs(U_NAnalyze(:,k,l,m) - U_exact))
         ASSOCIATE( U_NAnalyze_Abs => VECNORM(U_NAnalyze(1:3,k,l,m)) ,&
                    U_exact_Abs    => VECNORM(U_exact)             )
           L_Inf_Error(4) = MAX(L_Inf_Error(4),ABS(U_NAnalyze_Abs-U_exact_Abs))
           IntegrationWeight = wAnalyze(k)*wAnalyze(l)*wAnalyze(m)*J_NAnalyze(1,k,l,m)
           ! To sum over the elements, We compute here the square of the L_2 error
           L_2_Error(1:3) = L_2_Error(1:3) + (U_NAnalyze(:,k,l,m) - U_exact)*(U_NAnalyze(:,k,l,m) - U_exact)*IntegrationWeight
           L_2_Error(4)   = L_2_Error(4)   +((U_NAnalyze_Abs - U_exact_Abs)**2)*IntegrationWeight
         END ASSOCIATE
         ! Store exact solution (GL nodes) at current interpolation point
         U_NAnalyze_tmp(1:3,k,l,m) = U_exact(1:3)
       END DO ! k
     END DO ! l
   END DO ! m
   ! Map exact solution from GL to NodeType (usually G) node set for output to .h5 file
   CALL ChangeBasis3D(3,NAnalyze,PP_N,Vdm_GL_N,U_NAnalyze_tmp(1:3,:,:,:),BGFieldAnalytic(1:3,:,:,:,iElem))
END DO ! iElem=1,PP_nElems
#if USE_MPI
IF(MPIroot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE , L_2_Error   , 4 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)
  CALL MPI_REDUCE(MPI_IN_PLACE , L_Inf_Error , 4 , MPI_DOUBLE_PRECISION , MPI_MAX , 0 , MPI_COMM_PICLAS , iError)
ELSE
  CALL MPI_REDUCE(L_2_Error   , 0            , 4 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)
  CALL MPI_REDUCE(L_Inf_Error , 0            , 4 , MPI_DOUBLE_PRECISION , MPI_MAX , 0 , MPI_COMM_PICLAS , iError)
  ! in this case the receive value is not relevant.
END IF
#endif /*USE_MPI*/

! We normalize the L_2 Error with the Volume of the domain and take into account that we have to use the square root
L_2_Error = SQRT(L_2_Error/MeshVolume)

END SUBROUTINE CalcErrorSuperB


SUBROUTINE ExactFuncSuperB(ExactFunctionNumber,iCoilOrMagnet,x,resu)
!===================================================================================================================================
! Calculates the (analytical) solution for a given magnetostatic problem (for subsequent initial conditions or error calculation)
!
! ExactFunctionNumber corresponds to
!   1X : Coils (and linear conductors)
!   2X : Magnets
!===================================================================================================================================
! MODULES
USE MOD_Globals       ,ONLY: Abort,VECNORM,UNITVECTOR,CROSSNORM,DOTPRODUCT
USE MOD_Globals       ,ONLY: SphericalCoordinates,TransformVectorFromSphericalCoordinates
USE MOD_Globals_Vars  ,ONLY: Pi,mu0
USE MOD_SuperB_Vars   ,ONLY: CoilInfo,PermanentMagnetInfo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)   :: ExactFunctionNumber    ! Number of exact function to be used for the calculation of the analytical solution
INTEGER,INTENT(IN)   :: iCoilOrMagnet          ! Number of Coil or Magnet
REAL,INTENT(IN)      :: x(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Resu(1:3)    ! state in conservative variables
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: r,theta,phi
!===================================================================================================================================
SELECT CASE(ExactFunctionNumber)
CASE(10) ! linear conductor
  ! Calculate DOF vector "a" in local coordinate system
  ASSOCIATE( Lvec => CoilInfo(iCoilOrMagnet)%LengthVector(:)         ,&
             a    => x(1:3) - CoilInfo(iCoilOrMagnet)%BasePoint(1:3) ,&
             I    => CoilInfo(iCoilOrMagnet)%Current                  )
    ! project DOF vector in local coordinate system onto Lvec
    ASSOCIATE( b => UNITVECTOR(DOT_PRODUCT(a,Lvec)*Lvec)*VECNORM(a),&
               L => VECNORM(Lvec) )
      ! calculate the radial vector "c" in local coordinate system
      ASSOCIATE( c => a-b          ,&
                 R => VECNORM(a-b) )
        Resu(1:3) = CROSSNORM(b,c) * (mu0*I/(2.*Pi*R)) * (L/2.0) / SQRT((L/2.)**2 + R**2)
      END ASSOCIATE
    END ASSOCIATE
  END ASSOCIATE
CASE(20) ! Spherical hard magnet
  ! Calculate DOF vector "a" in local coordinate system
  ASSOCIATE( c1 => (mu0/3.)*VECNORM(PermanentMagnetInfo(1)%Magnetisation(1:3))*PermanentMagnetInfo(1)%Radius**3 ,&
              P => x(1:3) - PermanentMagnetInfo(1)%BasePoint(1:3)                                               )

    ! Get spherical coordinates
    CALL SphericalCoordinates(P,r,theta,phi)

    ! Transform vector
    ASSOCIATE( XHat => (c1/(r**3))*(/ 2*COS(theta) , SIN(theta) , 0.0 /) )
      CALL TransformVectorFromSphericalCoordinates(XHat,theta,phi,Resu(1:3))
    END ASSOCIATE
  END ASSOCIATE
CASE DEFAULT
  CALL abort(&
  __STAMP__&
  ,'ERROR in ExactFuncSuperB(): Cannot calculate L2/LInf error for case',IntInfoOpt=ExactFunctionNumber)
END SELECT

END SUBROUTINE ExactFuncSuperB


END MODULE MOD_SuperB_Tools
