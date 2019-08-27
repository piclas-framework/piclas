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

MODULE MOD_Filter
!===================================================================================================================================
! Module to handle Strukti's filter(s)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitFilter
  MODULE PROCEDURE InitFilter
END INTERFACE

INTERFACE Filter
  MODULE PROCEDURE Filter
END INTERFACE

INTERFACE FinalizeFilter
  MODULE PROCEDURE FinalizeFilter
END INTERFACE

PUBLIC :: InitFilter
PUBLIC :: Filter
PUBLIC :: FinalizeFilter
!===================================================================================================================================
PUBLIC :: DefineParametersFilter
CONTAINS

!==================================================================================================================================
!> Define parameters for surfaces (particle-sides)
!==================================================================================================================================
SUBROUTINE DefineParametersFilter()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Filter")

CALL prms%CreateIntOption(      'FilterType'      , 'TODO-DEFINE-PARAMETER', '0')
CALL prms%CreateRealArrayOption('HestFilterParam' , 'TODO-DEFINE-PARAMETER', '36. , 12. , 1.')

END SUBROUTINE DefineParametersFilter

SUBROUTINE InitFilter()
!===================================================================================================================================
! Initialize all necessary information to perform filtering
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Filter_Vars
USE MOD_Interpolation_Vars,ONLY:xGP,InterpolationInitIsDone
USE MOD_ReadInTools
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(FilterInitIsDone.OR.(.NOT.InterpolationInitIsDone))THEN
   CALL abort(&
       __STAMP__&
       ,'InitFilter not ready to be called or already called.',999,999.)
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT FILTER...'
CALL InitFilterBasis(PP_N,xGP)
FilterType = GETINT('FilterType','0')

IF(FilterType.GT.0) ALLOCATE(FilterMat(0:PP_N,0:PP_N))
SELECT CASE (FilterType)
CASE (1) ! Modal Filter (e.g. Hesthaven book)
  ! Read in modal filter parameter
  HestFilterParam = GETREALARRAY('HestFilterParam',3,'36.,12.,1.')
  CALL HestFilter()
END SELECT

FilterInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT FILTER DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitFilter



SUBROUTINE InitFilterBasis(N_in,xGP)
!===================================================================================================================================
! Initialize all necessary information to perform filtering
!===================================================================================================================================
! MODULES
USE MOD_Filter_Vars,ONLY:Vdm_Leg,sVdm_Leg
USE MOD_Basis, ONLY :buildLegendreVdm
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                         :: N_in
REAL,INTENT(IN),DIMENSION(0:N_in)          :: xGP
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
!  NODAL <--> MODAL
! Compute the 1D Vandermondematrix, needed to tranform the nodal basis into a modal (Legendre) basis
ALLOCATE(Vdm_Leg(0:N_in,0:N_in),sVdm_Leg(0:N_in,0:N_in))
CALL buildLegendreVdm(N_in,xGP,Vdm_Leg,sVdm_Leg)
END SUBROUTINE InitFilterBasis



SUBROUTINE HestFilter()
!===================================================================================================================================
! Modal Filter Function. Is called after each timestep. Could be used to controll aliasing instabilities. For details see NODALBOOK.
! The magnitude of the modal DOF are reduced, where DOF belonging to higher order are reduced more. The first DOF (=mean value)
! is NOT reduced to keep conservation. It is also possible to combine the Filter with a modal based indicator (Resolution/Persson
! indicator, to keep accuracy in resolved regions.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Filter_Vars,ONLY:HestFilterParam,FilterMat,sVdm_Leg,Vdm_Leg
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: alpha,s,eta,etac ! etac is the modal cutoff. Example:
                                        !   filtering all DOF from x ORDER to highest ORDER: alpha<0,s arbitrary, param(3)=x-1
                                        ! good start set of parameters: alpha=36, s=12, etac=1
                                        !   HestFilterParam=(/36.,12.,1./) for ini file
                                        !   Default is HestFilterParam(1)=0.
INTEGER             :: iDeg
!===================================================================================================================================
alpha = HestFilterParam(1)
s     = HestFilterParam(2)
etac  = HestFilterParam(3)/REAL(PP_N+1)

FilterMat = 0.
DO iDeg=0,MIN(INT(HestFilterParam(3))-1,PP_N)
  FilterMat(iDeg,iDeg) = 1.
END DO
IF(alpha.GE.0.) THEN
  DO iDeg=INT(HestFilterParam(3)),PP_N
    eta = REAL(iDeg+1)/REAL(PP_N+1)
    FilterMat(iDeg,iDeg) = EXP(-alpha*((eta-etac)/(1.-etac))**s)
  END DO
END IF
! Assemble ready-to-use nodal 1D filter matrix
FilterMat=MATMUL(MATMUL(Vdm_Leg,FilterMat),sVdm_Leg)
END SUBROUTINE HestFilter


SUBROUTINE Filter(U_in)
!===================================================================================================================================
! interpolate a 3D tensor product Lagrange basis defined by (N_in+1) 1D interpolation point positions xi_In(0:N_In)
! to another 3D tensor product node positions (number of nodes N_out+1)
! defined by (N_out+1) interpolation point  positions xi_Out(0:N_Out)
!  xi is defined in the 1DrefElem xi=[-1,1]
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Filter_Vars,ONLY:FilterType,FilterMat
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)  :: U_in(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,l,iElem
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N) :: U_Xi,U_Eta
!===================================================================================================================================
IF(filterType.EQ.0) RETURN ! do nothing
! Perform filtering
DO iElem=1,PP_nElems
!  U_Xi = 0.
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        U_Xi(:,i,j,k)       = FilterMat(i,0)*U_in(:,0,j,k,iElem)
        DO l=1,PP_N
          U_Xi(:,i,j,k)       = U_Xi(:,i,j,k)       + FilterMat(i,l)*U_in(:,l,j,k,iElem)
        END DO !l
      END DO !i
    END DO !j
  END DO !k
!  U_Eta= 0.
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        U_Eta(:,i,j,k)      = FilterMat(j,0)*U_Xi(:,i,0,k)
        DO l=1,PP_N
          U_Eta(:,i,j,k)      = U_Eta(:,i,j,k)      + FilterMat(j,l)*U_Xi(:,i,l,k)
        END DO !l
      END DO !i
    END DO !j
  END DO !k
!  U_in(:,:,:,:,iElem)=0.
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        U_in(:,i,j,k,iElem) = FilterMat(k,0)*U_Eta(:,i,j,0)
        DO l=1,PP_N
          U_in(:,i,j,k,iElem) = U_in(:,i,j,k,iElem) + FilterMat(k,l)*U_Eta(:,i,j,l)
        END DO !l
      END DO !i
    END DO !j
  END DO !k
END DO !iElem
END SUBROUTINE Filter

SUBROUTINE FinalizeFilter()
!===================================================================================================================================
! Deallocate global variable U (solution) and Ut (dg time derivative).
!===================================================================================================================================
! MODULES
USE MOD_Filter_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE(Vdm_Leg)
SDEALLOCATE(sVdm_Leg)
SDEALLOCATE(FilterMat)
FilterInitIsDone = .FALSE.
END SUBROUTINE FinalizeFilter

END MODULE MOD_Filter
