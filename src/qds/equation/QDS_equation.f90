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
#if USE_QDS_DG
#include "piclas.h"

MODULE MOD_QDS_Equation
!===================================================================================================================================
!> Contains the routines to
!> - determine exact functions for QDS-Ddetermine exact functions for QDS-DG
!===================================================================================================================================
! MODULES
!USE MOD_io_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE QDS_ExactFunc
  MODULE PROCEDURE QDS_ExactFunc
END INTERFACE

INTERFACE QDS_InitEquation
  MODULE PROCEDURE QDS_InitEquation
END INTERFACE

INTERFACE QDS_Q2U
  MODULE PROCEDURE QDS_Q2U
END INTERFACE

PUBLIC::QDS_Q2U
PUBLIC::QDS_ExactFunc
PUBLIC::QDS_InitEquation
!===================================================================================================================================
CONTAINS


SUBROUTINE QDS_InitEquation()
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
!USE MOD_Globals
USE MOD_Preproc
USE MOD_ReadInTools,        ONLY:GETINT
USE MOD_QDS_Equation_vars,  ONLY:QDSIniExactFunc,U_Face_old,QDSnVar
USE MOD_Mesh_Vars,          ONLY:nBCSides
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! set number of variables
!QDSnVar=40

! initial condition: exact function
QDSIniExactFunc = GETINT('QDSIniExactFunc','0')

! auxiliary face array for first order absorbing BC
ALLOCATE(U_Face_old(QDSnVar,0:PP_N,0:PP_N,1:nBCSides))
U_Face_old = 0.0

END SUBROUTINE QDS_InitEquation


SUBROUTINE QDS_ExactFunc(QDSExactFunction,t,tDeriv,x,resu)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_QDS_DG_Vars,        ONLY:QDSnVar_macro
#if (PP_TimeDiscMethod==1)
USE MOD_TimeDisc_Vars,      ONLY: dt
#endif /* (PP_TimeDiscMethod==1) */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: t
INTEGER,INTENT(IN)              :: tDeriv           ! determines the time derivative of the function
REAL,INTENT(IN)                 :: x(3)
INTEGER,INTENT(IN)              :: QDSExactFunction    ! determines the exact function
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Resu(1:QDSnVar_macro)    ! state in conservative variables
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: Resu_t(1:QDSnVar_macro),Resu_tt(1:QDSnVar_macro) ! state in conservative variables
REAL                            :: Cent(3)
REAL                            :: EOS_R
!===================================================================================================================================
Cent=x
SELECT CASE (QDSExactFunction)
#ifdef PARTICLES
CASE(0) ! Particles
  Resu=0.
  !resu(1:3)= x(1:3)!*x(1)
#endif
CASE(1) ! Constant
  Resu=1.
  Resu_t=0.
  Resu_tt=0.
CASE(100) ! SOD test case
  EOS_R=208.1101 ! [J/(kg K)] for Argon
  IF(x(1).GT.0.0)THEN ! right side of interface
    Resu(1)   = 0.125!*EOS_R/BoltzmannConst
    Resu(2:4) = 0.
    Resu(5)   = 0. !?
    Resu(6)   = 0.1E5/(0.125*EOS_R)
  ELSE ! left side of interface
    Resu(1)   = 1.!EOS_R/BoltzmannConst
    Resu(2:4) = 0.
    Resu(5)   = 0. !?
    Resu(6)   = 1E5/(EOS_R)
  END IF
CASE(200) ! SOD box case
  EOS_R=208.1101 ! [J/(kg K)] for Argon
  IF( .NOT. ALL((/ x(1).GT.0.8, x(1).LT.1.2,  x(2).GT.0.8, x(2).LT.1.2, x(3).GT.0.8, x(3).LT.1.2 /)) )THEN
    Resu(1)   = 0.125!*EOS_R/BoltzmannConst
    Resu(2:4) = 0.
    Resu(5)   = 0. !?
    Resu(6)   = 0.1E5/(0.125*EOS_R)
  ELSE ! left side of interface
    Resu(1)   = 2.!EOS_R/BoltzmannConst
    Resu(2:4) = 0.
    Resu(5)   = 0. !?
    Resu(6)   = 1E5/(EOS_R)
  END IF

CASE DEFAULT
  SWRITE(*,*)'Exact function not specified, QDSExactFunction = ',QDSExactFunction
END SELECT ! QDSExactFunction



#if (PP_TimeDiscMethod==1)
! For O3 RK, the boundary condition has to be adjusted
! Works only for O3 RK!!
SELECT CASE(tDeriv)
CASE(0)
  ! resu = g(t)
CASE(1)
  ! resu = g(t) + dt/3*g'(t)
  Resu=Resu + dt/3.*Resu_t
CASE(2)
  ! resu = g(t) + 3/4 dt g'(t) +5/16 dt^2 g''(t)
  Resu=Resu + 0.75*dt*Resu_t+5./16.*dt*dt*Resu_tt
CASE DEFAULT
  ! Stop, works only for 3 Stage O3 LS RK
  CALL abort(&
      __STAMP__&
      ,'Exactfuntion works only for 3 Stage O3 LS RK!',999,999.)
END SELECT
#endif
END SUBROUTINE QDS_ExactFunc




SUBROUTINE QDS_Q2U(Q,U)
!===================================================================================================================================
! Convert the QDS Q varibale vector to the QDS U variable vector
!===================================================================================================================================
! MODULES
USE MOD_QDS_DG_Vars,        ONLY:QDSSpeciesMass,GaussHermitWeiAbs,QDSSpecDOF,QDSnVar_macro
USE MOD_QDS_Equation_vars,  ONLY:QDSnVar
USE MOD_Globals_Vars,       ONLY:PI
USE MOD_Globals_Vars,       ONLY:BoltzmannConst
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: Q(QDSnVar_macro)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: U(QDSnVar)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         ::  L, iPart1, iPart2, iPart3
!===================================================================================================================================
L = 0
DO iPart1=1,2; DO iPart2=1,2; DO iPart3=1,2
  U(1+L) = Q(1)*&
                  GaussHermitWeiAbs(1,iPart1)     *&
                  GaussHermitWeiAbs(1,iPart2)     *&
                  GaussHermitWeiAbs(1,iPart3)/(PI*SQRT(PI))

  U(2+L) = U(1+L) &
           * (Q(2  ) /&
              Q(1  ) &
       + SQRT(2.*BoltzmannConst*Q(6)/QDSSpeciesMass)*GaussHermitWeiAbs(2,iPart1))

  U(3+L) = U(1+L) &
           * (Q(3  ) /&
              Q(1  ) &
       + SQRT(2.*BoltzmannConst*Q(6)/QDSSpeciesMass)*GaussHermitWeiAbs(2,iPart2))

  U(4+L) = U(1+L) &
           * (Q(4  ) /&
              Q(1  ) &
       + SQRT(2.*BoltzmannConst*Q(6)/QDSSpeciesMass)*GaussHermitWeiAbs(2,iPart3))

  U(5+L) =(QDSSpecDOF-3.)*BoltzmannConst*Q(6)/(QDSSpeciesMass*2.)

  L = L + 5
END DO; END DO; END DO
END SUBROUTINE QDS_Q2U


#endif /*USE_QDS_DG*/
END MODULE MOD_QDS_Equation
