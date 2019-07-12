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

MODULE MOD_Particle_Output
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitParticleOutput
  MODULE PROCEDURE InitParticleOutput
END INTERFACE

INTERFACE Visualize_Particles
  MODULE PROCEDURE Visualize_Particles
END INTERFACE

PUBLIC:: InitParticleOutput,Visualize_Particles
!===================================================================================================================================

CONTAINS

SUBROUTINE InitParticleOutput()
!===================================================================================================================================
! Initialize all output (and analyze) variables.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
!USE MOD_ReadInTools,ONLY:GetStr,GetLogical,GETINT
USE MOD_Particle_Output_Vars,ONLY:ParticleOutputInitIsDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF (ParticleOutputInitIsDone) THEN
  CALL abort(&
  __STAMP__&
  ,'InitOutput not ready to be called or already called.',999,999.)
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE OUTPUT...'

ParticleOutputInitIsDone =.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE OUTPUT DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitParticleOutput




SUBROUTINE Visualize_Particles(OutputTime)
!===================================================================================================================================
! Simple visualization of conservative variables
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Globals_Vars,ONLY:ProjectName
USE MOD_Particle_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: OutputTime
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: nParts, i, index_unit
CHARACTER(LEN=255)            :: FileString
!===================================================================================================================================

FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Particles',OutputTime))//'.dat'
!FileString=TRIM(INTSTAMP(TRIM(FileString),myRank))//'.dat'
! Visualize data
nParts = 0
DO i=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(i)) THEN
    nParts = nParts + 1
  END IF
END DO
!CALL WriteDataToTecplotBinary(NVisu,PP_nElems,PP_nVar,0,VarNames,Coords_NVisu(1:3,:,:,:,:),U_NVisu,TRIM(FileString))
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')" WRITE PARTICLE DATA TO TECPLOT ASCII FILE..."

index_unit = 45
OPEN(index_unit,FILE=TRIM(FileString),Status="REPLACE")
WRITE(index_unit,*)'TITLE = "',TRIM(ProjectName),'"'
WRITE(index_unit,'(A)')'VARIABLES = "x[m]" ,"y[m]" ,"z[m]" ,"v_x[m/s]" ,"v_y[m/s]" ,"v_z[m/s]" ,"Q[As]" ,"m[kg]" ,"Particle_Number"'
WRITE(index_unit,*)'ZONE T= "',TRIM(TIMESTAMP('Particles',OutputTime)),'"'
WRITE(index_unit,*)'I=',nParts,' J=1, K=1, F=POINT'
DO i=1,nParts
  WRITE(index_unit,'(8(1X,e19.12),1X,i0)')PartState(i,1),&
                                 PartState(i,2),&
                                 PartState(i,3),&
                                 PartState(i,4),&
                                 PartState(i,5),&
                                 PartState(i,6),&
                                 Species(PartSpecies(i))%ChargeIC,&
                                 Species(PartSpecies(i))%MassIC,&
                                 i
END DO
CLOSE(index_unit)
SWRITE(UNIT_stdOut,'(A)')"DONE!"

END SUBROUTINE Visualize_Particles



END MODULE MOD_Particle_Output
