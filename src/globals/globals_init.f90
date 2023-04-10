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

MODULE MOD_Globals_Init
!===================================================================================================================================
!> Provides parameters, used globally (please use EXTREMELY carefully!)
!===================================================================================================================================
! MODULES
#if USE_MPI
USE mpi
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
INTERFACE DefineParametersGlobals
  MODULE PROCEDURE DefineParametersGlobals
END INTERFACE

INTERFACE InitGlobals
  MODULE PROCEDURE InitGlobals
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------

PUBLIC :: DefineParametersGlobals,InitGlobals

!===================================================================================================================================
CONTAINS

!==================================================================================================================================
!> Define global parameters
!==================================================================================================================================
SUBROUTINE DefineParametersGlobals()
! MODULES
USE MOD_ReadInTools  ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Globals")
#if USE_READIN_CONSTANTS
CALL prms%CreateRealOption('c0'  , 'Speed of light (in vacuum) [m/s]' , '299792458.0'     )
CALL prms%CreateRealOption('eps' , 'Permittivity of vacuum [F/m])'    , '8.8541878176e-12')
CALL prms%CreateRealOption('mu'  , 'Permeability of vacuum [H/m])'    , '1.2566370614e-6' )
#else
CALL prms%CreateRealOption('c0'  , 'Speed of light (in vacuum) [m/s]. Hard coded (PICLAS_READIN_CONSTANTS=OFF)' , '299792458.0'     )
CALL prms%CreateRealOption('eps' , 'Permittivity of vacuum [F/m]). Hard coded (PICLAS_READIN_CONSTANTS=OFF)'    , '8.8541878176e-12')
CALL prms%CreateRealOption('mu'  , 'Permeability of vacuum [H/m]). Hard coded (PICLAS_READIN_CONSTANTS=OFF)'    , '1.2566370614e-6' )
#endif /*USE_READIN_CONSTANTS*/
END SUBROUTINE DefineParametersGlobals


!===================================================================================================================================
!> Pre-compute required constants
!===================================================================================================================================
SUBROUTINE InitGlobals()
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: LogFile,UNIT_StdOut,UNIT_logOut,Logging,myRank,abort
USE MOD_Globals_Vars ,ONLY: c,eps0,mu0,ProjectName,memory
#if USE_READIN_CONSTANTS
USE MOD_Globals_Vars ,ONLY: c_inv,c2,c2_inv,smu0,RelativisticLimit
USE MOD_ReadInTools  ,ONLY: GETREAL
#else
USE MOD_ReadInTools  ,ONLY: PrintOption
#endif /*USE_READIN_CONSTANTS*/
#if USE_MPI
USE MOD_Globals      ,ONLY: MPIRoot
#endif /*USE_MPI*/
#if defined(PARTICLES)
USE MOD_Globals      ,ONLY: nGlobalNbrOfParticles
#endif /*defined(PARTICLES)*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: OpenStat
CHARACTER(LEN=8)  :: StrDate
CHARACTER(LEN=10) :: StrTime
LOGICAL           :: LogIsOpen
REAL              :: c_test
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT GLOBALS ...'

! RAM monitor
memory = 0.

#if USE_READIN_CONSTANTS
! Natural constants
c      = GETREAL('c0')
eps0   = GETREAL('eps')
mu0    = GETREAL('mu')
c_inv  = 1./c
c2     = c*c
smu0   = 1./mu0
c2_inv = 1./c2
RelativisticLimit = (1e6/299792458.0)**2 * c2 ! corresponds to 0.3% speed of light (which is 1000000 when speed of light is 299792458)
#else
CALL PrintOption('Speed of light (in vacuum) c0 [m/s]. To change, set PICLAS_READIN_CONSTANTS=ON' , 'FIXED' , RealOpt=c)
CALL PrintOption('Permittivity (of vacuum)  eps [F/m]. To change, set PICLAS_READIN_CONSTANTS=ON' , 'FIXED' , RealOpt=eps0)
CALL PrintOption('Permeability (of vacuum)   mu [F/m]. To change, set PICLAS_READIN_CONSTANTS=ON' , 'FIXED' , RealOpt=mu0)
#endif /*USE_READIN_CONSTANTS*/

! Sanity check
c_test = 1./SQRT(eps0*mu0)
IF(.NOT.ALMOSTEQUALRELATIVE(c_test,c,10E-8))THEN
  SWRITE(*,*) "ERROR: c does not equal 1/sqrt(eps*mu)!"
  SWRITE(*,*) "c:", c
  SWRITE(*,*) "mu:", mu0
  SWRITE(*,*) "eps:", eps0
  SWRITE(*,*) "1/sqrt(eps*mu):", c_test
  CALL abort(__STAMP__,' Speed of light coefficients does not match!')
END IF

#if defined(PARTICLES)
nGlobalNbrOfParticles=0
nGlobalNbrOfParticles(4)=HUGE(nGlobalNbrOfParticles(4))
#endif /*defined(PARTICLES)*/

! Open file for logging
IF(Logging)THEN
  INQUIRE(UNIT=UNIT_LogOut,OPENED=LogIsOpen)
  IF(.NOT.LogIsOpen)THEN
    WRITE(LogFile,'(A,A1,I6.6,A4)')TRIM(ProjectName),'_',myRank,'.log'
    OPEN(UNIT=UNIT_logOut,  &
         FILE=LogFile,      &
         STATUS='UNKNOWN',  &
         ACTION='WRITE',    &
         POSITION='APPEND', &
         IOSTAT=OpenStat)
    CALL DATE_AND_TIME(StrDate,StrTime)
    WRITE(UNIT_logOut,*)
    WRITE(UNIT_logOut,'(132("#"))')
    WRITE(UNIT_logOut,*)
    WRITE(UNIT_logOut,*)'STARTED LOGGING FOR PROC',myRank,' ON ',StrDate(7:8),'.',StrDate(5:6),'.',StrDate(1:4),' | ',&
                        StrTime(1:2),':',StrTime(3:4),':',StrTime(5:10)
  END IF !logIsOpen
END IF  ! Logging

SWRITE(UNIT_stdOut,'(A)')' INIT GLOBALS DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitGlobals


END MODULE MOD_Globals_Init
