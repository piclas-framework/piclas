!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "piclas.h"

!===================================================================================================================================
!> Routines to create several mappings:
!> * Mapping that contains the distribution of DG elements
!> * Mappings from all available variables to the ones that should be calculated and visualized for volume and surface
!> * Mappings from the available boundaries to the ones that should be visualized
!===================================================================================================================================
MODULE MOD_Posti_Mappings
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE Build_mapDepToCalc_mapAllVarsToVisuVars
  MODULE PROCEDURE Build_mapDepToCalc_mapAllVarsToVisuVars
END INTERFACE

INTERFACE Build_mapBCSides
  MODULE PROCEDURE Build_mapBCSides
END INTERFACE

PUBLIC:: Build_mapDepToCalc_mapAllVarsToVisuVars
PUBLIC:: Build_mapBCSides

CONTAINS

!===================================================================================================================================
!> This routine builds the mappings from the total number of variables available for visualization to number of calculation
!> and visualization variables.
!>  1. Read 'VarName' options from the parameter file. This are the quantities that will be visualized.
!>  2. Initialize the dependecy table
!>  3. check wether gradients are needed for any quantity. If this is the case, remove the conservative quantities from the
!>     dependecies of the primitive quantities (the primitive quantities are available directly, since the DGTimeDerivative_weakForm
!>     will be executed.
!>  4. build the 'mapDepToCalc' that holds for each quantity that will be calculated the index in 'UCalc' array (0 if not calculated)
!>  5. build the 'mapAllVarsToVisuVars' that holds for each quantity that will be visualized the index in 'UVisu' array (0 if not visualized)
!>  6. build the 'mapAllBCNamesToVisuBCNames' that holds for each available boundary the visualization index
!===================================================================================================================================
SUBROUTINE Build_mapDepToCalc_mapAllVarsToVisuVars()
USE MOD_Globals
USE MOD_Visu_Vars
USE MOD_ReadInTools     ,ONLY: GETSTR,CountOption
USE MOD_StringTools     ,ONLY: STRICMP
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iVar,iVar2
CHARACTER(LEN=255)  :: VarName
!CHARACTER(LEN=255)  :: BoundaryName
CHARACTER(LEN=20)   :: format
!===================================================================================================================================
! Read Varnames from parameter file and fill
!   mapAllVarsToVisuVars = map, which stores at position x the position/index of the x.th quantity in the UVisu array
!             if a quantity is not visualized it is zero
SDEALLOCATE(mapAllVarsToVisuVars)
SDEALLOCATE(mapAllVarsToSurfVisuVars)
ALLOCATE(mapAllVarsToVisuVars(1:nVarAll))
ALLOCATE(mapAllVarsToSurfVisuVars(1:nVarAll))
mapAllVarsToVisuVars = 0
mapAllVarsToSurfVisuVars = 0
nVarVisu = 0
nVarSurfVisuAll = 0
! Compare varnames that should be visualized with available varnames
DO iVar=1,CountOption("VarName")
  VarName = GETSTR("VarName")
  DO iVar2=1,nVarAll
    IF (STRICMP(VarName, VarnamesAll(iVar2))) THEN
      nVarSurfVisuAll = nVarSurfVisuAll + 1
      mapAllVarsToSurfVisuVars(iVar2) = nVarSurfVisuAll
      IF (iVar2.LE.nVarDep) THEN
        IF(DepSurfaceOnly(iVar2).EQ.1) CYCLE
      END IF
      nVarVisu = nVarVisu + 1
      mapAllVarsToVisuVars(iVar2) = nVarVisu
    END IF
  END DO
END DO

!! check whether gradients are needed for any quantity
!DO iVar=1,nVarDep
!  IF (mapAllVarsToSurfVisuVars(iVar).GT.0) THEN
!    withDGOperator = withDGOperator .OR. (DepTable(iVar,0).GT.0)
!  END IF
!END DO

!! Calculate all dependencies:
!! For each quantity copy from all quantities that this quantity depends on the dependencies.
!DO iVar=1,nVarDep
!  DepTable(iVar,iVar) = 1
!  DO iVar2=1,iVar-1
!    IF (DepTable(iVar,iVar2).EQ.1) &
!      DepTable(iVar,:) = MAX(DepTable(iVar,:), DepTable(iVar2,:))
!  END DO
!END DO

!! print the dependecy table
!SWRITE(*,*) "Dependencies: ", withDGOperator
!WRITE(format,'(I2)') SIZE(DepTable,2)
!DO iVar=1,nVarDep
!  SWRITE (*,'('//format//'I2,A)') DepTable(iVar,:), " "//TRIM(VarnamesAll(iVar))
!END DO

!! Build :
!!   mapDepToCalc = map, which stores at position x the position/index of the x.th quantity in the UCalc array
!!             if a quantity is not calculated it is zero
!SDEALLOCATE(mapDepToCalc)
!ALLOCATE(mapDepToCalc(1:nVarDep))
!mapDepToCalc = 0
!DO iVar=1,nVarDep
!  IF (mapAllVarsToSurfVisuVars(iVar).GT.0) THEN
!    mapDepToCalc = MAX(mapDepToCalc,DepTable(iVar,1:nVarDep))
!  END IF
!END DO
!! enumerate mapDepToCalc
!nVarCalc = 0
!DO iVar=1,nVarDep
!  IF (mapDepToCalc(iVar).GT.0) THEN
!    nVarCalc = nVarCalc + 1
!    mapDepToCalc(iVar) = nVarCalc
!  END IF
!END DO

!! check if any varnames changed
!changedVarNames = .TRUE.
!IF (ALLOCATED(mapAllVarsToSurfVisuVars_old).AND.(SIZE(mapAllVarsToSurfVisuVars).EQ.SIZE(mapAllVarsToSurfVisuVars_old))) THEN
!  changedVarNames = .NOT.ALL(mapAllVarsToSurfVisuVars.EQ.mapAllVarsToSurfVisuVars_old)
!END IF
!SDEALLOCATE(mapAllVarsToSurfVisuVars_old)
!ALLOCATE(mapAllVarsToSurfVisuVars_old(1:nVarAll))
!mapAllVarsToSurfVisuVars_old = mapAllVarsToSurfVisuVars


! print the mappings
WRITE(format,'(I0)') nVarAll
!SWRITE (*,'(A,'//format//'I3)') "mapDepToCalc             ",mapDepToCalc
SWRITE (*,'(A,'//format//'I3)') "mapAllVarsToVisuVars     ",mapAllVarsToVisuVars
!SWRITE (*,'(A,'//format//'I3)') "mapAllVarsToSurfVisuVars ",mapAllVarsToSurfVisuVars


!!---------------------- Surface visualization ----------------------------!

!! Build the mapping for the surface visualization
!! mapAllBCNamesToVisuBCNames(iBC) stores the ascending visualization index of the all boundaries. 0 means no visualization.
!! nBCNamesVisu is the number of boundaries to be visualized.
!SDEALLOCATE(mapAllBCNamesToVisuBCNames)
!ALLOCATE(mapAllBCNamesToVisuBCNames(1:nBCNamesAll))
!mapAllBCNamesToVisuBCNames = 0
!nBCNamesVisu = 0
!! Compare boundary names that should be visualized with available varnames
!DO iVar=1,CountOption("BoundaryName")
!  BoundaryName = GETSTR("BoundaryName")
!  DO iVar2=1,nBCNamesAll
!    IF (STRICMP(BoundaryName, BCNamesAll(iVar2))) THEN
!      mapAllBCNamesToVisuBCNames(iVar2) = nBCNamesVisu+1
!      nBCNamesVisu = nBCNamesVisu + 1
!    END IF
!  END DO
!END DO

!! Set flag indicating if surface visualization is needed
!doSurfVisu = nBCNamesVisu.GT.0
!SWRITE (*,*) 'doSurfVisu: ',doSurfVisu

!! check if any boundary changed
!changedBCnames = .TRUE.
!IF (ALLOCATED(mapAllBCNamesToVisuBCNames_old).AND.(SIZE(mapAllBCNamesToVisuBCNames).EQ.SIZE(mapAllBCNamesToVisuBCNames_old))) THEN
!  changedBCnames = .NOT.ALL(mapAllBCNamesToVisuBCNames.EQ.mapAllBCNamesToVisuBCNames_old)
!END IF
!SDEALLOCATE(mapAllBCNamesToVisuBCNames_old)
!ALLOCATE(mapAllBCNamesToVisuBCNames_old(1:nBCNamesAll))
!mapAllBCNamesToVisuBCNames_old = mapAllBCNamesToVisuBCNames


!SWRITE (*,'(A,'//format//'I3)') "mapAllBCNamesToVisuBCNames ",mapAllBCNamesToVisuBCNames

END SUBROUTINE Build_mapDepToCalc_mapAllVarsToVisuVars

!===================================================================================================================================
!> This routine builds mappings that give for each BC side the index of the visualization side, seperate by FV and DG.
!> We also store the number of visualization sides per BC name.
!===================================================================================================================================
SUBROUTINE Build_mapBCSides()
USE MOD_Visu_Vars
USE MOD_Mesh_Vars   ,ONLY: nBCSides,BC,SideToElem,BoundaryName
USE MOD_StringTools ,ONLY: STRICMP
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iBC,iElem,iSide
!===================================================================================================================================

! Build surface visualization mappings.
! mapAllBCSidesToDGVisuBCSides/FV(iBCSide) contains the ascending index of the visualization boundary sides. They are sorted after the boundary name.
! 0 means no visualization of this boundary side.
! nSidesPerBCNameVisu_DG/FV(iBCNamesVisu) contains how many boundary sides belong to each boundary that should be visualized.
SDEALLOCATE(mapAllBCSidesToDGVisuBCSides)
SDEALLOCATE(nSidesPerBCNameVisu_DG)
ALLOCATE(mapAllBCSidesToDGVisuBCSides(1:nBCSides))
ALLOCATE(nSidesPerBCNameVisu_DG(1:nBCNamesVisu))
mapAllBCSidesToDGVisuBCSides = 0
nSidesPerBCNameVisu_DG = 0
nBCSidesVisu_DG = 0
IF (doSurfVisu) THEN
  DO iBC=1,nBCNamesAll    ! iterate over all bc names
    IF (mapAllBCNamesToVisuBCNames(iBC).GT.0) THEN
      DO iSide=1,nBCSides   ! iterate over all bc sides
        IF (STRICMP(BoundaryName(BC(iSide)),BCNamesAll(iBC))) THEN ! check if side is of specific boundary name
          iElem = SideToElem(S2E_ELEM_ID,iSide)
          nBCSidesVisu_DG = nBCSidesVisu_DG + 1
          mapAllBCSidesToDGVisuBCSides(iSide) = nBCSidesVisu_DG
          nSidesPerBCNameVisu_DG(mapAllBCNamesToVisuBCNames(iBC)) = nSidesPerBCNameVisu_DG(mapAllBCNamesToVisuBCNames(iBC)) + 1
        END IF
      END DO
    END IF
  END DO
END IF ! doSurfVisu

END SUBROUTINE Build_mapBCSides

END MODULE MOD_Posti_Mappings
