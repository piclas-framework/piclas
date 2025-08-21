!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_AnalyzeField_Tools
!===================================================================================================================================
! Contains the Poynting Vector Integral part for the power analysis of the field vector
!===================================================================================================================================
! USE MOD_Globals, ONLY:UNIT_stdout
! USE MOD_PreProc
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: CalculateBoundaryFieldOutput
#if !(USE_FV) || (USE_HDG)
PUBLIC :: CalcPotentialEnergy
PUBLIC :: CalcPotentialEnergy_Dielectric
#endif /*no FV alone - !(USE_FV) || (USE_HDG)*/
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Determine the field boundary output (BFO) values for the given iBC
!===================================================================================================================================
SUBROUTINE CalculateBoundaryFieldOutput(iBC,Time,BoundaryFieldOutput)
! MODULES
#if USE_HDG
USE MOD_Mesh_Vars ,ONLY: BoundaryType
USE MOD_Equation  ,ONLY: ExactFunc
#else
USE MOD_Globals   ,ONLY: abort
#endif /*USE_HDG*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: iBC
REAL,INTENT(IN)    :: Time
REAL,INTENT(OUT)   :: BoundaryFieldOutput(1:PP_nVar)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_HDG
INTEGER           :: BCType,BCState
#else
INTEGER           :: dummy
#endif /*USE_HDG*/
!===================================================================================================================================
BoundaryFieldOutput=0.!Initialize
#if USE_HDG
#if (PP_nVar==1)
BCType =BoundaryType(iBC,BC_TYPE)
BCState=BoundaryType(iBC,BC_STATE)
ASSOCIATE( x => (/0., 0., 0./) )
  SELECT CASE(BCType)
  CASE(2) ! exact BC = Dirichlet BC !! ExactFunc via BCState (time is optional)
    CALL ExactFunc(BCState , x , BoundaryFieldOutput , t=Time)
  CASE(4) ! exact BC = Dirichlet BC !! Zero potential
    BoundaryFieldOutput = 0.
  CASE(5,51,52,60) ! exact BC = Dirichlet BC !! ExactFunc via RefState (time is optional)
    CALL ExactFunc(  -1    , x , BoundaryFieldOutput , t=Time  , iRefState=BCState)
  CASE(6) ! exact BC = Dirichlet BC !! ExactFunc via RefState (Time is optional)
    CALL ExactFunc(  -2    , x , BoundaryFieldOutput , t=time  , iRefState=BCState)
  CASE(50) ! exact BC = Dirichlet BC !! ExactFunc via bias voltage DC
    CALL ExactFunc(  -5    , x , BoundaryFieldOutput )
  END SELECT ! BCType
END ASSOCIATE
#else
CALL abort(__STAMP__,'CalculateBoundaryFieldOutput is not implemented for PP_nVar>1')
#endif /*PP_nVar==1*/
#else
CALL abort(__STAMP__,'CalculateBoundaryFieldOutput is not implemented for other equation systems yet (only HDG)')
! Suppress warnings
dummy=iBC
dummy=INT(Time)
#endif /*USE_HDG*/
END SUBROUTINE CalculateBoundaryFieldOutput


#if !(USE_FV) || (USE_HDG)
#if (PP_nVar==8)
SUBROUTINE CalcPotentialEnergy(WEl, WMag, Wphi, Wpsi)
#else
SUBROUTINE CalcPotentialEnergy(WEl)
#endif /*PP_nVar=8*/
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars          ,ONLY: nElems, N_VolMesh, offSetElem
USE MOD_Interpolation_Vars ,ONLY: N_Inter
#if (PP_nVar==8)
USE MOD_Globals_Vars       ,ONLY: smu0
#endif /*PP_nVar=8*/
USE MOD_Globals_Vars       ,ONLY: eps0
USE MOD_DG_Vars            ,ONLY: U_N,N_DG_Mapping
#if !(USE_HDG)
#endif /*HDG*/
#if USE_HDG
#if PP_nVar==1
!USE MOD_Equation_Vars      ,ONLY: E
#elif PP_nVar==3
USE MOD_Equation_Vars      ,ONLY: B
#else
USE MOD_Equation_Vars      ,ONLY: B,E
#endif /*PP_nVar==1*/
#else
USE MOD_PML_Vars           ,ONLY: DoPML,isPMLElem
#endif /*USE_HDG*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: WEl
#if (PP_nVar==8)
REAL,INTENT(OUT)                :: WMag,Wpsi,Wphi
#endif /*PP_nVar=8*/
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iElem,Nloc
INTEGER           :: i,j,k
REAL              :: WEl_tmp, WMag_tmp, E_abs
#if !(USE_HDG)
REAL              :: B_abs , Phi_abs, Psi_abs
#endif
#if USE_MPI
REAL              :: RD
#endif
#if (PP_nVar==8)
REAL              :: Wphi_tmp, Wpsi_tmp
#endif /*PP_nVar=8*/
!===================================================================================================================================

Wel=0.
#if (PP_nVar==8)
WMag=0.
Wphi=0.
Wpsi=0.
#endif /*PP_nVar=8*/
DO iElem=1,nElems
#if !(USE_HDG)
  IF(DoPML)THEN
    IF(isPMLElem(iElem))CYCLE
  END IF
#endif
  !--- Calculate and save volume of element iElem
  WEl_tmp=0.
  WMag_tmp=0.
#if (PP_nVar==8)
    Wphi_tmp = 0.
    Wpsi_tmp = 0.
#endif /*PP_nVar=8*/
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  ASSOCIATE( wGP => N_Inter(Nloc)%wGP )
  DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
#if USE_HDG
  ASSOCIATE( Ex  => U_N(iElem)%E(1,i,j,k) ,&
             Ey  => U_N(iElem)%E(2,i,j,k) ,&
             Ez  => U_N(iElem)%E(3,i,j,k) )
#endif /*USE_HDG*/
! in electromagnetische felder by henke 2011 - springer
! WMag = 1/(2mu) * int_V B^2 dV
#if USE_HDG
#if PP_nVar==1
    E_abs = DOTPRODUCT(U_N(iElem)%E(1:3,i,j,k))
#elif PP_nVar==3
    B_abs = B(1,i,j,k,iElem)*B(1,i,j,k,iElem) + B(2,i,j,k,iElem)*B(2,i,j,k,iElem) + B(3,i,j,k,iElem)*B(3,i,j,k,iElem)
#else /*PP_nVar==4*/
    E_abs = Ex*Ex + Ey*Ey + Ez*Ez
    B_abs = B(1,i,j,k,iElem)*B(1,i,j,k,iElem) + B(2,i,j,k,iElem)*B(2,i,j,k,iElem) + B(3,i,j,k,iElem)*B(3,i,j,k,iElem)
#endif /*PP_nVar==1*/
#else
    E_abs = U_N(iElem)%U(1,i,j,k)**2 + U_N(iElem)%U(2,i,j,k)**2 + U_N(iElem)%U(3,i,j,k)**2
#endif /*USE_HDG*/

#if (PP_nVar==8)
    B_abs = U_N(iElem)%U(4,i,j,k)**2 + U_N(iElem)%U(5,i,j,k)**2 + U_N(iElem)%U(6,i,j,k)**2
    Phi_abs = U_N(iElem)%U(7,i,j,k)*U_N(iElem)%U(7,i,j,k)
    Psi_abs = U_N(iElem)%U(8,i,j,k)*U_N(iElem)%U(8,i,j,k)
#endif /*PP_nVar=8*/
#if USE_HDG
#if PP_nVar==3
      WMag_tmp = WMag_tmp + wGP(i)*wGP(j)*wGP(k) / N_VolMesh(iElem)%sJ(i,j,k) * B_abs
#elif PP_nVar==4
      WMag_tmp = WMag_tmp + wGP(i)*wGP(j)*wGP(k) / N_VolMesh(iElem)%sJ(i,j,k) * B_abs
#endif /*PP_nVar==3*/
#endif /*USE_HDG*/
    WEl_tmp  = WEl_tmp  + wGP(i)*wGP(j)*wGP(k)   / N_VolMesh(iElem)%sJ(i,j,k) * E_abs
#if (PP_nVar==8)
    WMag_tmp = WMag_tmp + wGP(i)*wGP(j)*wGP(k)   / N_VolMesh(iElem)%sJ(i,j,k) * B_abs
    Wphi_tmp = Wphi_tmp + wGP(i)*wGP(j)*wGP(k)   / N_VolMesh(iElem)%sJ(i,j,k) * Phi_abs
    Wpsi_tmp = Wpsi_tmp + wGP(i)*wGP(j)*wGP(k)   / N_VolMesh(iElem)%sJ(i,j,k) * Psi_abs
#endif /*PP_nVar=8*/
#if USE_HDG
  END ASSOCIATE
#endif /*USE_HDG*/
  END DO; END DO; END DO
  END ASSOCIATE
  WEl = WEl + WEl_tmp
#if (PP_nVar==8)
  WMag = WMag + WMag_tmp
  Wphi = Wphi + Wphi_tmp
  Wpsi = Wpsi + Wpsi_tmp
#endif /*PP_nVar=8*/
END DO

WEl = WEl * eps0 * 0.5
#if (PP_nVar==8)
WMag = WMag * smu0 * 0.5
! caution: change of coefficients for divergence energies
Wphi = Wphi * eps0*0.5
Wpsi = Wpsi * smu0*0.5
#endif /*PP_nVar=8*/

#if USE_MPI
! todo: only one reduce with array
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,WEl  , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
#if (PP_nVar==8)
  CALL MPI_REDUCE(MPI_IN_PLACE,WMag , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,Wphi , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,Wpsi , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
#endif /*PP_nVar=8*/
ELSE
  CALL MPI_REDUCE(WEl         ,RD   , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
#if (PP_nVar==8)
  CALL MPI_REDUCE(WMag        ,RD   , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
  CALL MPI_REDUCE(Wphi        ,RD   , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
  CALL MPI_REDUCE(Wpsi        ,RD   , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
#endif /*PP_nVar=8*/
END IF
#endif /*USE_MPI*/

END SUBROUTINE CalcPotentialEnergy


#if (PP_nVar==8)
SUBROUTINE CalcPotentialEnergy_Dielectric(WEl, WMag, Wphi, Wpsi)
#else
SUBROUTINE CalcPotentialEnergy_Dielectric(WEl)
#endif /*PP_nVar=8*/
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
#if USE_HDG
#if PP_nVar==3 || PP_nVar==4
USE MOD_Dielectric_vars    ,ONLY: DielectricVol
#endif /*PP_nVar==3 or 4*/
#endif /*USE_HDG or PP_nVar==8*/
USE MOD_Mesh_Vars          ,ONLY: nElems, N_VolMesh, offSetElem
USE MOD_Interpolation_Vars ,ONLY: N_Inter
#if (PP_nVar==8)
USE MOD_Dielectric_vars    ,ONLY: DielectricVol
USE MOD_Globals_Vars       ,ONLY: smu0
#endif /*PP_nVar=8*/
USE MOD_Globals_Vars       ,ONLY: eps0
USE MOD_Dielectric_vars    ,ONLY: isDielectricElem,DielectricVol,ElemToDielectric
USE MOD_DG_Vars            ,ONLY: U_N,N_DG_Mapping
#if !(USE_HDG)
#endif /*PP_nVar=8*/
#if USE_HDG
#if PP_nVar==1
!USE MOD_Equation_Vars      ,ONLY: E
#elif PP_nVar==3
USE MOD_Equation_Vars      ,ONLY: B
#else
USE MOD_Equation_Vars      ,ONLY: B,E
#endif /*PP_nVar==1*/
#else
USE MOD_PML_Vars           ,ONLY: DoPML,isPMLElem
#endif /*USE_HDG*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: WEl
#if (PP_nVar==8)
REAL,INTENT(OUT)                :: WMag,Wpsi,Wphi
#endif /*PP_nVar=8*/
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iElem,Nloc
INTEGER           :: i,j,k
REAL              :: WEl_tmp, WMag_tmp, E_abs
#if !(USE_HDG)
REAL              :: B_abs , Phi_abs, Psi_abs
#endif
#if USE_MPI
REAL              :: RD
#endif
#if (PP_nVar==8)
REAL              :: Wphi_tmp, Wpsi_tmp
#endif /*PP_nVar=8*/
!===================================================================================================================================

Wel=0.
#if (PP_nVar==8)
WMag=0.
Wphi=0.
Wpsi=0.
#endif /*PP_nVar=8*/

DO iElem=1,nElems
#if !(USE_HDG)
  IF(DoPML)THEN
    IF(isPMLElem(iElem))CYCLE
  END IF
#endif
  !--- Calculate and save volume of element iElem
  WEl_tmp=0.
  WMag_tmp=0.
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  ASSOCIATE( wGP => N_Inter(Nloc)%wGP )
  IF(isDielectricElem(iElem))THEN
    DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
#if USE_HDG
  ASSOCIATE( Ex  => U_N(iElem)%E(1,i,j,k) ,&
             Ey  => U_N(iElem)%E(2,i,j,k) ,&
             Ez  => U_N(iElem)%E(3,i,j,k) )
#else
  ASSOCIATE( Ex  => U_N(iElem)%U(1,i,j,k) ,&
             Ey  => U_N(iElem)%U(2,i,j,k) ,&
             Ez  => U_N(iElem)%U(3,i,j,k) ,&
             Bx  => U_N(iElem)%U(4,i,j,k) ,&
             By  => U_N(iElem)%U(5,i,j,k) ,&
             Bz  => U_N(iElem)%U(6,i,j,k) )
#endif /*USE_HDG*/
      ! in electromagnetische felder by henke 2011 - springer
      ! WMag = 1/(2mu) * int_V B^2 dV
#if USE_HDG
#if PP_nVar==1
      E_abs = Ex*Ex + Ey*Ey + Ez*Ez
#elif PP_nVar==3
      B_abs = B(1,i,j,k,iElem)*B(1,i,j,k,iElem) + B(2,i,j,k,iElem)*B(2,i,j,k,iElem) + B(3,i,j,k,iElem)*B(3,i,j,k,iElem)
#else /*PP_nVar==4*/
      E_abs = E(1,i,j,k,iElem)*E(1,i,j,k,iElem) + E(2,i,j,k,iElem)*E(2,i,j,k,iElem) + E(3,i,j,k,iElem)*E(3,i,j,k,iElem)
      B_abs = B(1,i,j,k,iElem)*B(1,i,j,k,iElem) + B(2,i,j,k,iElem)*B(2,i,j,k,iElem) + B(3,i,j,k,iElem)*B(3,i,j,k,iElem)
#endif /*PP_nVar==1*/
#else
      E_abs = Ex*Ex + Ey*Ey + Ez*Ez
#endif /*USE_HDG*/

#if (PP_nVar==8)
      B_abs = Bx*Bx + By*By + Bz*Bz
      Phi_abs = U_N(iElem)%U(7,i,j,k)*U_N(iElem)%U(7,i,j,k)
      Psi_abs = U_N(iElem)%U(8,i,j,k)*U_N(iElem)%U(8,i,j,k)
#endif /*PP_nVar=8*/
#if USE_HDG
#if PP_nVar==3
      WMag_tmp = WMag_tmp + wGP(i)*wGP(j)*wGP(k) / N_VolMesh(iElem)%sJ(i,j,k) * B_abs / DielectricVol(ElemToDielectric(iElem))%DielectricMu( i,j,k)
#elif PP_nVar==4
      WMag_tmp = WMag_tmp + wGP(i)*wGP(j)*wGP(k) / N_VolMesh(iElem)%sJ(i,j,k) * B_abs / DielectricVol(ElemToDielectric(iElem))%DielectricMu( i,j,k)
#endif /*PP_nVar==3*/
#endif /*USE_HDG*/
      WEl_tmp  = WEl_tmp  + wGP(i)*wGP(j)*wGP(k) / N_VolMesh(iElem)%sJ(i,j,k) * E_abs * DielectricVol(ElemToDielectric(iElem))%DielectricEps(i,j,k)
#if (PP_nVar==8)
      WMag_tmp = WMag_tmp + wGP(i)*wGP(j)*wGP(k) / N_VolMesh(iElem)%sJ(i,j,k) * B_abs / DielectricVol(ElemToDielectric(iElem))%DielectricMu(i,j,k)
      Wphi_tmp = Wphi_tmp + wGP(i)*wGP(j)*wGP(k) / N_VolMesh(iElem)%sJ(i,j,k) * Phi_abs * DielectricVol(ElemToDielectric(iElem))%DielectricEps(i,j,k)
      Wpsi_tmp = Wpsi_tmp + wGP(i)*wGP(j)*wGP(k) / N_VolMesh(iElem)%sJ(i,j,k) * Psi_abs / DielectricVol(ElemToDielectric(iElem))%DielectricMu(i,j,k)
#endif /*PP_nVar=8*/
  END ASSOCIATE
    END DO; END DO; END DO
  ELSE
    DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
#if USE_HDG
  ASSOCIATE( Ex  => U_N(iElem)%E(1,i,j,k) ,&
             Ey  => U_N(iElem)%E(2,i,j,k) ,&
             Ez  => U_N(iElem)%E(3,i,j,k) )
#else
  ASSOCIATE( Ex  => U_N(iElem)%U(1,i,j,k) ,&
             Ey  => U_N(iElem)%U(2,i,j,k) ,&
             Ez  => U_N(iElem)%U(3,i,j,k) ,&
             Bx  => U_N(iElem)%U(4,i,j,k) ,&
             By  => U_N(iElem)%U(5,i,j,k) ,&
             Bz  => U_N(iElem)%U(6,i,j,k) )
#endif /*USE_HDG*/
      ! in electromagnetische felder by henke 2011 - springer
      ! WMag = 1/(2mu) * int_V B^2 dV
#if USE_HDG
#if PP_nVar==1
      E_abs = Ex*Ex + Ey*Ey + Ez*Ez
#elif PP_nVar==3
      B_abs = B(1,i,j,k,iElem)*B(1,i,j,k,iElem) + B(2,i,j,k,iElem)*B(2,i,j,k,iElem) + B(3,i,j,k,iElem)*B(3,i,j,k,iElem)
#else /*PP_nVar==4*/
      E_abs = E(1,i,j,k,iElem)*E(1,i,j,k,iElem) + E(2,i,j,k,iElem)*E(2,i,j,k,iElem) + E(3,i,j,k,iElem)*E(3,i,j,k,iElem)
      B_abs = B(1,i,j,k,iElem)*B(1,i,j,k,iElem) + B(2,i,j,k,iElem)*B(2,i,j,k,iElem) + B(3,i,j,k,iElem)*B(3,i,j,k,iElem)
#endif /*PP_nVar==1*/
#else
      E_abs = Ex*Ex + Ey*Ey + Ez*Ez
#endif /*USE_HDG*/

#if (PP_nVar==8)
      B_abs = Bx*Bx + By*By + Bz*Bz
      Phi_abs = U_N(iElem)%U(7,i,j,k)*U_N(iElem)%U(7,i,j,k)
      Psi_abs = U_N(iElem)%U(8,i,j,k)*U_N(iElem)%U(8,i,j,k)
#endif /*PP_nVar=8*/
#if USE_HDG
#if PP_nVar==3
      WMag_tmp = WMag_tmp + wGP(i)*wGP(j)*wGP(k) / N_VolMesh(iElem)%sJ(i,j,k) * B_abs
#elif PP_nVar==4
      WMag_tmp = WMag_tmp + wGP(i)*wGP(j)*wGP(k) / N_VolMesh(iElem)%sJ(i,j,k) * B_abs
#endif /*PP_nVar==3*/
#endif /*USE_HDG*/
      WEl_tmp  = WEl_tmp  + wGP(i)*wGP(j)*wGP(k) / N_VolMesh(iElem)%sJ(i,j,k) * E_abs
#if (PP_nVar==8)
      WMag_tmp = WMag_tmp + wGP(i)*wGP(j)*wGP(k) / N_VolMesh(iElem)%sJ(i,j,k) * B_abs
      Wphi_tmp = Wphi_tmp + wGP(i)*wGP(j)*wGP(k) / N_VolMesh(iElem)%sJ(i,j,k) * Phi_abs
      Wpsi_tmp = Wpsi_tmp + wGP(i)*wGP(j)*wGP(k) / N_VolMesh(iElem)%sJ(i,j,k) * Psi_abs
#endif /*PP_nVar=8*/
  END ASSOCIATE
    END DO; END DO; END DO
  END IF
  END ASSOCIATE


    WEl = WEl + WEl_tmp
#if (PP_nVar==8)
    WMag = WMag + WMag_tmp
#endif /*PP_nVar=8*/




END DO

WEl = WEl * eps0 * 0.5
#if (PP_nVar==8)
WMag = WMag * smu0 * 0.5
! caution: change of coefficients for divergence energies
Wphi = Wphi * eps0*0.5
Wpsi = Wpsi * smu0*0.5
#endif /*PP_nVar=8*/

#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,WEl  , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
#if (PP_nVar==8)
  CALL MPI_REDUCE(MPI_IN_PLACE,WMag , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
#endif /*PP_nVar=8*/
ELSE
  CALL MPI_REDUCE(WEl         ,RD   , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
#if (PP_nVar==8)
  CALL MPI_REDUCE(WMag        ,RD   , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
#endif /*PP_nVar=8*/
END IF
#endif /*USE_MPI*/

END SUBROUTINE CalcPotentialEnergy_Dielectric
#endif /*no FV alone - !(USE_FV) || (USE_HDG)*/

END MODULE MOD_AnalyzeField_Tools