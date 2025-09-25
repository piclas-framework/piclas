!=================================================================================================================================
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

MODULE MOD_HDF5_Output_Particles_HDG
#if USE_HDG && defined(PARTICLES)
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MOD_IO_HDF5
USE MOD_HDF5_output
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: AddBRElectronFluidToPartSource
!===================================================================================================================================

CONTAINS


SUBROUTINE AddBRElectronFluidToPartSource()
!===================================================================================================================================
! Add BR electron fluid density to PartSource for output to state.h5
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals            ,ONLY: abort
USE MOD_Mesh_Vars          ,ONLY: nElems, offSetElem
USE MOD_PreProc
USE MOD_HDG_Vars           ,ONLY: ElemToBRRegion,RegionElectronRef
USE MOD_DG_Vars            ,ONLY: U_N
USE MOD_PICDepo_Vars       ,ONLY: PS_N
USE MOD_DG_Vars            ,ONLY: N_DG_Mapping
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iElem,RegionID,Nloc
INTEGER :: i,j,k
REAL    :: source_e
!===================================================================================================================================
! Loop over all elements and all DOF and add the contribution of the BR electron density to PartSource
DO iElem=1,nElems
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  ! BR electron fluid region
  RegionID=ElemToBRRegion(iElem)
  IF (RegionID.GT.0) THEN
    DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
#if ((USE_HDG) && (PP_nVar==1))
      source_e = U_N(iElem)%U(1,i,j,k)-RegionElectronRef(2,RegionID)
#else
      CALL abort(__STAMP__,' CalculateBRElectronsPerCell only implemented for electrostatic HDG!')
#endif
      IF (source_e .LT. 0.) THEN
        source_e = RegionElectronRef(1,RegionID) &         !--- boltzmann relation (electrons as isothermal fluid!)
            * EXP( (source_e) / RegionElectronRef(3,RegionID) )
      ELSE
        source_e = RegionElectronRef(1,RegionID) &         !--- linearized boltzmann relation at positive exponent
            * (1. + ((source_e) / RegionElectronRef(3,RegionID)) )
      END IF
      PS_N(iElem)%PartSource(4,i,j,k) = PS_N(iElem)%PartSource(4,i,j,k) - source_e
    END DO; END DO; END DO
  END IF
END DO ! iElem=1,PP_nElems

END SUBROUTINE AddBRElectronFluidToPartSource

#endif /*USE_HDG && defined(PARTICLES)*/
END MODULE MOD_HDF5_Output_Particles_HDG
