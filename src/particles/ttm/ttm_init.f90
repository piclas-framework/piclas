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

#ifdef PARTICLES
MODULE MOD_TTMInit
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitTTM
  MODULE PROCEDURE InitTTM
END INTERFACE

INTERFACE FinalizeTTM
  MODULE PROCEDURE FinalizeTTM
END INTERFACE

INTERFACE InitIMD_TTM_Coupling
  MODULE PROCEDURE InitIMD_TTM_Coupling
END INTERFACE

PUBLIC::InitTTM,FinalizeTTM,InitIMD_TTM_Coupling
!===================================================================================================================================
PUBLIC::DefineParametersTTM

CONTAINS

!==================================================================================================================================
!> Define parameters for Two temperature Model
!==================================================================================================================================
SUBROUTINE DefineParametersTTM()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("TTM")

CALL prms%CreateLogicalOption(  'DoImportTTMFile'      , 'Read IMD Two-Temperature Model (TTM) data (FD grid data with\n electron'&
                                                          //'temperature and other field data)','.FALSE.')

CALL prms%CreateStringOption(   'TTMLogFile'           , 'TTM Log file path','no file specified')
CALL prms%CreateStringOption(   'TTMFile'              , 'TTM Data file path','no file found')
CALL prms%CreateIntArrayOption( 'TTMGridFDdim'         , 'Number of FD grid cells in each direction (/x,y,z/)','0 , 0 , 0')

CALL prms%CreateRealOption(     'TTMElemBaryTolerance' , 'TTM FD bary center tolerance to DG bary center. The tolerance is used '//&
                                                         'for finding the corresponding DG element to which the FD data is saved.'&
                                                       , '1e-6')

END SUBROUTINE DefineParametersTTM

SUBROUTINE InitTTM()
!===================================================================================================================================
! IMD Two Temperature Model (TTM):
! Read TTM data from ASCII file and create FD mesh which is then mapped onto the TTM-DG solution for output
! IMD TTM data format and target arrays for data input (from ASCII data file *.ttm) are shown below
! -------------------------------------------------------------------------------------------------------------------------------
! The properties and array index are chosen as follows
!        1      2    3       4  5      6       7       8       9    10   11       ->  TTM_FD(1:11,ix,iy,iz) -> TTM_Cell_XX
!  1 2 3 4      5    6       7  8      9       10      11      12   13   14 15    ->  tmp_array(1:15)
! #x y z natoms temp md_temp xi source v_com.x v_com.y v_com.z fd_k fd_g Z  proc
! -------------------------------------------------------------------------------------------------------------------------------
! read the following fields ...
!
! 1.)  x:         0 to Nx-1 FD cells in x-direction                          [-]                not stored
! 2.)  y:         0 to Ny-1 FD cells in y-direction                          [-]                not stored
! 3.)  z:         0 to Nz-1 FD cells in z-direction                          [-]                not stored
! 4.)  natoms:    number of atoms in the FD cell                             [-]                stored in TTM_FD(1,.. -> TTM_Cell_1
! 5.)  temp:      temperature of the electrons (fluid model)                 [eV]               stored in TTM_FD(2,.. -> TTM_Cell_2
! 6.)  md_temp:   temperature of the atoms (MD simulated discrete atoms)     [eV]               stored in TTM_FD(3,.. -> TTM_Cell_3
! 7.)  xi: ?                                                                 [?]                stored in TTM_FD(4,.. -> TTM_Cell_4
! 8.)  source:    energy input source, e.g., laser energy density per volume [?]                stored in TTM_FD(5,.. -> TTM_Cell_5
! 9.)  v_com.x:   barycentric velocity of all MD atoms in x-direction        [Ångström/10,18fs] stored in TTM_FD(6,.. -> TTM_Cell_6
! 10.) v_com.y:   barycentric velocity of all MD atoms in y-direction        [Ångström/10,18fs] stored in TTM_FD(7,.. -> TTM_Cell_7
! 11.) v_com.z:   barycentric velocity of all MD atoms in z-direction        [Ångström/10,18fs] stored in TTM_FD(8,.. -> TTM_Cell_8
! 12.) fd_k:      heat conduction coefficient                                [?]                stored in TTM_FD(9,.. -> TTM_Cell_9
! 13.) fd_g:      electron-phonon coupling coefficient                       [?]                stored in TTM_FD(10,. -> TTM_Cell_10
! 14.) Z:         averged charge of the atoms within the FD cell             [e]                stored in TTM_FD(11,. -> TTM_Cell_11
! 15.) proc:      rank number of MPI process                                 [-]                not stored

! Derived quantities are Debye length, warm/cold plasma frequency, PIC time steps estimation
! n_e(ElectronDensity)           = N[natoms]*charge[Z]/TTMCellVolume              -> TTM_Cell_12
! omega_pe_cold(PlasmaFrequency) = w_peTTM=sqrt(neTTM*e^2/(me0*eps0 ))            -> TTM_Cell_13
! omega_pe_warm(PlasmaFrequency) = w_peTTMwarm = w_peTTM + 3 * kB * TeTTM / me0   -> TTM_Cell_14
! dt_PICcold(TimeStep)          = dtPICTM=0.2./w_peTTM                          -> TTM_Cell_15
! dt_PICwarm(TimeStep)'         = dtPICTMwarm=0.2./w_peTTMwarm                  -> TTM_Cell_16
! T_e(ElectronTempInKelvin)      = T[eV]/8.621738E-5 (eV -> K)                    -> TTM_Cell_17
! lambda_D(DebyeLength)          = sqrt(eps0*kB*TeTTM      ./(K_to_eV*neTTM*e^2))
!                                = sqrt(eps0*kB*TeTTM_in_K./(         neTTM*e^2)) -> TTM_Cell_18
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_ReadInTools
USE MOD_TTM_Vars
USE MOD_Restart_Vars      ,ONLY: DoRestart
USE MOD_Mesh_Vars         ,ONLY: ElemBaryNGeo
USE MOD_Globals_Vars      ,ONLY: BoltzmannConst,eps0
USE MOD_IO_HDF5           ,ONLY: AddToElemData,ElementOut
USE MOD_Restart_Vars      ,ONLY: RestartFile
USE MOD_TTM_Vars          ,ONLY: DoImportTTMFile
USE MOD_TTM_Vars          ,ONLY: TTM_Cell_1,TTM_Cell_2,TTM_Cell_3,TTM_Cell_4,TTM_Cell_5,TTM_Cell_6,TTM_Cell_7,TTM_Cell_8
USE MOD_TTM_Vars          ,ONLY: TTM_Cell_9,TTM_Cell_10,TTM_Cell_11,TTM_Cell_12,TTM_Cell_13,TTM_Cell_14,TTM_Cell_15
USE MOD_TTM_Vars          ,ONLY: TTM_Cell_16,TTM_Cell_17,TTM_Cell_18
USE MOD_StringTools       ,ONLY: STRICMP
USE MOD_IO_hdf5
USE MOD_HDF5_Input
USE MOD_ReadInTools       ,ONLY: PrintOption
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars  ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER        :: ioUnit
INTEGER        :: IndNum
INTEGER        :: i,j,ix,iy,iz
INTEGER        :: io_error
INTEGER        :: iLineTTM
INTEGER        :: iElem,iElemFD
CHARACTER(255) :: StrTmp
REAL           :: tmp_array(1:15)
REAL           :: fd_hx,fd_hy,fd_hz
REAL           :: dt_PIC_min(1:2),dt_PIC_max(1:2)
REAL           :: omega_pe_min(1:2),omega_pe_max(1:2)
REAL           :: lambda_D_min,lambda_D_max

INTEGER                        :: nVal(15),iVar,nRestartVars
REAL,ALLOCATABLE               :: ElemData_loc(:,:),tmp(:)
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesElemData_loc(:) ! Element data variable names
!===================================================================================================================================
IF(TTMInitIsDone)THEN
  LBWRITE(UNIT_stdOut,'(A)') "InitTTM already called."
  RETURN
END IF
LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT TTM (IMD Two-Temperature Model) ...'

DoImportTTMFile = GETLOGICAL('DoImportTTMFile','.FALSE.')
TTMLogFile      = TRIM(GETSTR('TTMLogFile','no file specified'))

IF(DoImportTTMFile.EQV..TRUE.)THEN
  IF(DoRestart.EQV..FALSE.)THEN ! if no restart is performed: read the TTM field data from *.ttm file
    TTMFile   =TRIM(GETSTR('TTMFile','no file found'))
    IF((TRIM(TTMFile).NE.'no file found').AND.(TRIM(TTMLogFile).NE.'no file specified'))THEN
      ! get FD grid array dimension in x, y and z direction
      TTMGridFDdim=GETINTARRAY('TTMGridFDdim',3,'0 , 0 , 0')
      ! get FD grid element bary center tolerance for checking distance to elem bary center distance of DG cells
      TTMElemBaryTolerance=GETREAL('TTMElemBaryTolerance','1e-6')
      ! get FD grid cell sizes
      CALL GetFDGridCellSize(fd_hx,fd_hy,fd_hz)
      IF(ANY(TTMGridFDdim(1:3).LE.0).OR.ANY((/fd_hx,fd_hy,fd_hz/).LE.0.0))THEN
        DoImportTTMFile=.FALSE.
      ELSE
        FD_nElems=TTMGridFDdim(1)*TTMGridFDdim(2)*TTMGridFDdim(3)
        CALL PrintOption('FD_nElems','OUTPUT',IntOpt=FD_nElems)
        TTMCellVolume=fd_hx*fd_hy*fd_hz
        CALL PrintOption('TTM mesh cell volume','OUTPUT',RealOpt=TTMCellVolume)
        TTMTotalVolume=REAL(FD_nElems)*TTMCellVolume
        CALL PrintOption('TTM mesh total volume','OUTPUT',RealOpt=TTMTotalVolume)
        ALLOCATE( TTM_FD(1:11,TTMGridFDdim(1),TTMGridFDdim(2),TTMGridFDdim(3)) )
        TTM_FD=0.0
        ALLOCATE( ElemBaryFD(3,FD_nElems) )
        ElemBaryFD=0.0
        ALLOCATE( ElemIndexFD(3,FD_nElems) )
        ElemIndexFD=0
        ALLOCATE( ElemIsDone(FD_nElems) )
        ElemIsDone=.FALSE.
        LBWRITE(UNIT_stdOut,'(A,A)') " Reading TTM data from file (TTMFile): ",TRIM(TTMFile)
#if USE_MPI
        IF(.NOT.MPIRoot)THEN
          CALL abort(&
              __STAMP__&
              ,'ERROR: Cannot SetParticlePosition in multi-core environment for SpaceIC=IMD!')
        END IF
#endif /*USE_MPI*/
        OPEN(NEWUNIT=ioUnit,FILE=TRIM(TTMFile),STATUS='OLD',ACTION='READ',IOSTAT=io_error)
        IF(io_error.NE.0)THEN
          CALL abort(&
              __STAMP__&
              ,'Cannot open specified File (ttm data) from '//TRIM(TTMFile))
        END IF
        ! get index of *.ttm data file, e.g., "5" (corresponds to specific time stamp calculated from timestep)
        IndNum=INDEX(TTMFile, '/',BACK = .TRUE.)
        IF(IndNum.GT.0)THEN
          StrTmp=TRIM(TTMFile(IndNum+1:LEN(TTMFile)))
          IndNum=INDEX(StrTmp,'.',BACK = .TRUE.)
          IF(IndNum.GT.0)THEN
            StrTmp=StrTmp(1:IndNum-1)
            IndNum=INDEX(StrTmp,'.')
            IF(IndNum.GT.0)THEN
              StrTmp=StrTmp(IndNum+1:LEN(StrTmp))
            END IF
          END IF
        END IF

        READ(StrTmp,*,iostat=io_error)  i
        TTMNumber = i
        CALL PrintOption('TTMNumber','OUTPUT',IntOpt=TTMNumber)

        iLineTTM=0
        DO i=1,1 ! header lines: currently 1 -> adjust?
          READ(ioUnit,'(A)',IOSTAT=io_error)StrTmp
        END DO
        DO
          READ(ioUnit,*,IOSTAT=io_error) tmp_array(1:15)
          IF(io_error>0)THEN
            SWRITE(UNIT_stdOut,'(A,I10)') " io_error=",io_error
            SWRITE(UNIT_stdOut,'(A,I10)') " A problem reading the TTM file occured in line: ",iLineTTM+2
            SWRITE(UNIT_stdOut,'(A1)') ' ['
            DO j=1,15
              SWRITE(UNIT_stdOut,'(A1,ES25.14E3)',ADVANCE='NO') " ",tmp_array(j)
            END DO
            SWRITE(UNIT_stdOut,'(A1)') ']'
            CALL abort(__STAMP__,'Error reading specified File (ttm data): '//TRIM(TTMFile))
          ELSE IF(io_error<0)THEN
            LBWRITE(UNIT_stdOut,'(A,I10,A)') " End of file reached: ",iLineTTM+1," -> EXIT"
            EXIT
          ELSE ! io_error = 0
          END IF
          iLineTTM=iLineTTM+1 ! *.ttm data file: line counter

          ! FD index values go from 0 to N_ijk-1
          ix=INT(tmp_array(1),4)+1 ! FD grid index in x
          iy=INT(tmp_array(2),4)+1 ! FD grid index in y
          iz=INT(tmp_array(3),4)+1 ! FD grid index in z

          IF( (ix.LE.TTMGridFDdim(1)) .AND. (iy.LE.TTMGridFDdim(2)) .AND. (iz.LE.TTMGridFDdim(3)) ) THEN
            TTM_FD(1:11,ix,iy,iz) = tmp_array(4:14)
            ElemIndexFD(1,iLineTTM)=ix
            ElemIndexFD(2,iLineTTM)=iy
            ElemIndexFD(3,iLineTTM)=iz
            ElemBaryFD(1,iLineTTM)=(REAL(ix)-0.5)*fd_hx
            ElemBaryFD(2,iLineTTM)=(REAL(iy)-0.5)*fd_hy
            ElemBaryFD(3,iLineTTM)=(REAL(iz)-0.5)*fd_hz
          ELSE
            LBWRITE(UNIT_stdOut,'(A,I10)') " A problem reading the TTM file occured in line: ",iLineTTM+1
            DoImportTTMFile=.FALSE.
            EXIT
          END IF
        END DO
        IF(FD_nElems.NE.ix*iy*iz)THEN
          SWRITE(UNIT_stdOut,'(A)')     ' Error: FD_nElems.NE.ix*iy*iz'
          SWRITE(UNIT_stdOut,'(A,I10)') ' FD_nElems = ',FD_nElems
          SWRITE(UNIT_stdOut,'(A,I10)') ' ix*iy*iz  = ',ix*iy*iz
          CALL abort(__STAMP__&
              ,'Number of FD elements FD_nElems does not comply with the number of FD elements read from *.ttm file ['&
              //TRIM(TTMFile)//']')
        END IF
        CLOSE(ioUnit)
        LBWRITE(UNIT_stdOut,'(A,I10,A)') " Read ",iLineTTM," lines from file (+1 header line)"
        CALL PrintOption('TTM FD grid x'      , 'OUTPUT' , IntOpt=ix)
        CALL PrintOption('TTM FD grid y'      , 'OUTPUT' , IntOpt=iy)
        CALL PrintOption('TTM FD grid z'      , 'OUTPUT' , IntOpt=iz)
        CALL PrintOption('TTM FD grid points' , 'OUTPUT' , IntOpt=iLineTTM)

        IF(FD_nElems.EQ.PP_nElems)THEN ! same number of nodes in FD grid points and DG bary centres -> assume they coincide
          LBWRITE(UNIT_stdOut,'(A)') ' Searching for all FD cells within the DG mesh in order to map the values ...'
          ! the local DG solution in physical and reference space

          ALLOCATE( TTM_Cell_1(1:PP_nElems) )
          TTM_Cell_1=0.0
          CALL AddToElemData(ElementOut,'TTM_N[natoms]'           ,RealArray=TTM_Cell_1(1:PP_nElems))

          ALLOCATE( TTM_Cell_2(1:PP_nElems) )
          TTM_Cell_2=0.0
          CALL AddToElemData(ElementOut,'TTM_T_e[temp]'           ,RealArray=TTM_Cell_2(1:PP_nElems))

          ALLOCATE( TTM_Cell_3(1:PP_nElems) )
          TTM_Cell_3=0.0
          CALL AddToElemData(ElementOut,'TTM_T_i[md_temp]'        ,RealArray=TTM_Cell_3(1:PP_nElems))

          ALLOCATE( TTM_Cell_4(1:PP_nElems) )
          TTM_Cell_4=0.0
          CALL AddToElemData(ElementOut,'TTM_[xi]'                ,RealArray=TTM_Cell_4(1:PP_nElems))

          ALLOCATE( TTM_Cell_5(1:PP_nElems) )
          TTM_Cell_5=0.0
          CALL AddToElemData(ElementOut,'TTM_[source]'            ,RealArray=TTM_Cell_5(1:PP_nElems))

          ALLOCATE( TTM_Cell_6(1:PP_nElems) )
          TTM_Cell_6=0.0
          CALL AddToElemData(ElementOut,'TTM_[v_com.x]'           ,RealArray=TTM_Cell_6(1:PP_nElems))

          ALLOCATE( TTM_Cell_7(1:PP_nElems) )
          TTM_Cell_7=0.0
          CALL AddToElemData(ElementOut,'TTM_[v_com.y]'           ,RealArray=TTM_Cell_7(1:PP_nElems))

          ALLOCATE( TTM_Cell_8(1:PP_nElems) )
          TTM_Cell_8=0.0
          CALL AddToElemData(ElementOut,'TTM_[v_com.z]'           ,RealArray=TTM_Cell_8(1:PP_nElems))

          ALLOCATE( TTM_Cell_9(1:PP_nElems) )
          TTM_Cell_9=0.0
          CALL AddToElemData(ElementOut,'TTM_[fd_k]'              ,RealArray=TTM_Cell_9(1:PP_nElems))

          ALLOCATE( TTM_Cell_10(1:PP_nElems) )
          TTM_Cell_10=0.0
          CALL AddToElemData(ElementOut,'TTM_[fd_g]'              ,RealArray=TTM_Cell_10(1:PP_nElems))

          ALLOCATE( TTM_Cell_11(1:PP_nElems) )
          TTM_Cell_11=0.0
          CALL AddToElemData(ElementOut,'TTM_charge[Z]'           ,RealArray=TTM_Cell_11(1:PP_nElems))

          ALLOCATE( TTM_Cell_12(1:PP_nElems) )
          TTM_Cell_12=0.0
          CALL AddToElemData(ElementOut,'TTM_n_e(ElectronDensity)',RealArray=TTM_Cell_12(1:PP_nElems))

          ALLOCATE( TTM_Cell_13(1:PP_nElems) )
          TTM_Cell_13=0.0
          CALL AddToElemData(ElementOut,'TTM_omega_pe_cold(PlasmaFrequency)',RealArray=TTM_Cell_13(1:PP_nElems))

          ALLOCATE( TTM_Cell_14(1:PP_nElems) )
          TTM_Cell_14=0.0
          CALL AddToElemData(ElementOut,'TTM_omega_pe_warm(PlasmaFrequency)',RealArray=TTM_Cell_14(1:PP_nElems))

          ALLOCATE( TTM_Cell_15(1:PP_nElems) )
          TTM_Cell_15=0.0
          CALL AddToElemData(ElementOut,'TTM_dt_PIC_cold(TimeStep)',RealArray=TTM_Cell_15(1:PP_nElems))

          ALLOCATE( TTM_Cell_16(1:PP_nElems) )
          TTM_Cell_16=0.0
          CALL AddToElemData(ElementOut,'TTM_dt_PIC_warm(TimeStep)',RealArray=TTM_Cell_16(1:PP_nElems))

          ALLOCATE( TTM_Cell_17(1:PP_nElems) )
          TTM_Cell_17=0.0
          CALL AddToElemData(ElementOut,'TTM_T_e(ElectronTempInKelvin)',RealArray=TTM_Cell_17(1:PP_nElems))

          ALLOCATE( TTM_Cell_18(1:PP_nElems) )
          TTM_Cell_18=0.0
          CALL AddToElemData(ElementOut,'TTM_lambda_D(DebyeLength)',RealArray=TTM_Cell_18(1:PP_nElems))

          DO iElem=1,PP_nElems
            DO iELemFD=1,FD_nElems
              IF(ElemIsDone(iElemFD))CYCLE ! the element has already been found
              IF(ALMOSTEQUALRELATIVE(ElemBaryNGeo(1,iElem),ElemBaryFD(1,iElemFD),TTMElemBaryTolerance).EQV..FALSE.) CYCLE
              IF(ALMOSTEQUALRELATIVE(ElemBaryNGeo(2,iElem),ElemBaryFD(2,iElemFD),TTMElemBaryTolerance).EQV..FALSE.) CYCLE
              IF(ALMOSTEQUALRELATIVE(ElemBaryNGeo(3,iElem),ElemBaryFD(3,iElemFD),TTMElemBaryTolerance).EQV..FALSE.) CYCLE
              ElemIsDone(iElemFD)=.TRUE.

              ! Map FD data to DG cell data (cell constant value)
              TTM_Cell_1( iElem)=TTM_FD( 1,ElemIndexFD(1,iElemFD),ElemIndexFD(2,iElemFD),ElemIndexFD(3,iElemFD))
              TTM_Cell_2( iElem)=TTM_FD( 2,ElemIndexFD(1,iElemFD),ElemIndexFD(2,iElemFD),ElemIndexFD(3,iElemFD))
              TTM_Cell_3( iElem)=TTM_FD( 3,ElemIndexFD(1,iElemFD),ElemIndexFD(2,iElemFD),ElemIndexFD(3,iElemFD))
              TTM_Cell_4( iElem)=TTM_FD( 4,ElemIndexFD(1,iElemFD),ElemIndexFD(2,iElemFD),ElemIndexFD(3,iElemFD))
              TTM_Cell_5( iElem)=TTM_FD( 5,ElemIndexFD(1,iElemFD),ElemIndexFD(2,iElemFD),ElemIndexFD(3,iElemFD))
              TTM_Cell_6( iElem)=TTM_FD( 6,ElemIndexFD(1,iElemFD),ElemIndexFD(2,iElemFD),ElemIndexFD(3,iElemFD))
              TTM_Cell_7( iElem)=TTM_FD( 7,ElemIndexFD(1,iElemFD),ElemIndexFD(2,iElemFD),ElemIndexFD(3,iElemFD))
              TTM_Cell_8( iElem)=TTM_FD( 8,ElemIndexFD(1,iElemFD),ElemIndexFD(2,iElemFD),ElemIndexFD(3,iElemFD))
              TTM_Cell_9( iElem)=TTM_FD( 9,ElemIndexFD(1,iElemFD),ElemIndexFD(2,iElemFD),ElemIndexFD(3,iElemFD))
              TTM_Cell_10(iElem)=TTM_FD(10,ElemIndexFD(1,iElemFD),ElemIndexFD(2,iElemFD),ElemIndexFD(3,iElemFD))
              TTM_Cell_11(iElem)=TTM_FD(11,ElemIndexFD(1,iElemFD),ElemIndexFD(2,iElemFD),ElemIndexFD(3,iElemFD))

              ! set derived quantities
              ! 'n_e(ElectronDensity)'
              ! N[natoms]*charge[Z]/TTMCellVolume
              TTM_Cell_12(iElem) = TTM_Cell_11(iElem)*TTM_Cell_1(iElem)/TTMCellVolume

              ! 'omega_pe_cold(PlasmaFrequency)'
              ! w_peTTM=sqrt(neTTM*e^2/(me0*eps0 ))
              TTM_Cell_13(iElem) = SQRT(TTM_Cell_12(iElem)*ElementaryCharge**2/(ElectronMass*eps0))

              ! 'omega_pe_warm(PlasmaFrequency)'
              ! w_peTTMwarm = w_peTTM + 3 * kB * TeTTM_in_K / me0
              TTM_Cell_14(iElem) = TTM_Cell_13(iElem) + 3 * BoltzmannConst * TTM_Cell_2(iElem) / ElectronMass / 8.621738E-5

              ! 'dt_PIC_cold(TimeStep)'
              ! dtPICTTM=0.2./w_peTTM
              IF(INT(TTM_Cell_11(iElem),4).EQ.0.0)THEN ! when no atoms are present, then no electron density is calculated
                TTM_Cell_15(iElem) = 0.0
              ELSE
                TTM_Cell_15(iElem) = 0.2 / TTM_Cell_13(iElem) ! 13 depends on 12 which depends only on 11
              END IF

              ! 'dt_PIC_warm(TimeStep)'
              ! dtPICTTMwarm=0.2./w_peTTMwarm
              IF(TTM_Cell_14(iElem).LE.0.0)THEN
                TTM_Cell_16(iElem) = 0.0
              ELSE
                TTM_Cell_16(iElem) = 0.2 / TTM_Cell_14(iElem)
              END IF

              ! 'T_e(ElectronTempInKelvin)'
              ! eV -> K
              TTM_Cell_17(iElem) = TTM_Cell_2(iElem)/8.621738E-5

              ! 'lambda_D(DebyeLength)'
              ! lambda_D_TTM=sqrt(eps0*kB*TeTTM./(K_to_eV*neTTM*e^2)) or sqrt(eps0*kB*TeTTM_in_K./(neTTM*e^2))
              IF(INT(TTM_Cell_11(iElem),4).EQ.0)THEN ! when no atoms are present, then no electron density is calculated
                TTM_Cell_18(iElem) = 0.0
              ELSE
                TTM_Cell_18(iElem) = SQRT( eps0*BoltzmannConst*TTM_Cell_17(iElem)/&
                    (TTM_Cell_12(iElem)*ElementaryCharge**2)) ! 12 depends only on 11
              END IF
            END DO
          END DO
          LBWRITE(UNIT_stdOut,'(A)') ' Done.'
          IF(ANY(ElemIsDone).EQV..FALSE.)THEN
            LBWRITE(UNIT_stdOut,'(A)') " NOT all IMD TTM FD cells have been located within the DG grid!"
          ELSE
            LBWRITE(UNIT_stdOut,'(A)') " All IMD TTM FD cells have been located within the DG grid!"
          END IF
          SDEALLOCATE(TTM_FD)
          SDEALLOCATE(ElemBaryFD)
          SDEALLOCATE(ElemIndexFD)
          SDEALLOCATE(ElemIsDone)
        ELSE ! use swap-mesh or interpolate the FD grid to the DG solution
          SWRITE(UNIT_stdOut,'(A,I10)') "FD_nElems=",FD_nElems
          SWRITE(UNIT_stdOut,'(A,I10)') "PP_nElems=",PP_nElems
          CALL abort(__STAMP__&
              ,'Error interpolating ttm data (FD grid) to DG solution. FD_nElems.NE.PP_nElems. Different elements not implemented')
        END IF
      END IF
    ELSE
      DoImportTTMFile=.FALSE.
    END IF
  ELSE
    LBWRITE(UNIT_stdOut,'(A)')' INIT TTM: data will be read from restart file!'
    IF(MPIRoot)THEN
      nRestartVars=0
      CALL OpenDataFile(RestartFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
      CALL GetArrayAndName('ElemData','VarNamesAdd',nVal,tmp,VarNamesElemData_loc)
      CALL CloseDataFile()
      IF (ALLOCATED(VarNamesElemData_loc)) THEN
        ALLOCATE(ElemData_loc(nVal(1),nVal(2)))
        ElemData_loc = RESHAPE(tmp,(/nVal(1),nVal(2)/))
        DO iVar=1,nVal(1)
          ! check the variable names
          IF(STRICMP(VarNamesElemData_loc(iVar),"TTM_N[natoms]")) THEN
            nRestartVars=nRestartVars+1
            ALLOCATE( TTM_Cell_1(1:PP_nElems) )
            TTM_Cell_1 = REAL(ElemData_loc(iVar,:))
            CALL AddToElemData(ElementOut,'TTM_N[natoms]'                      ,RealArray=TTM_Cell_1( 1:PP_nElems))

          ELSEIF (STRICMP(VarNamesElemData_loc(iVar),"TTM_T_e[temp]")) THEN
            nRestartVars=nRestartVars+1
            ALLOCATE( TTM_Cell_2(1:PP_nElems) )
            TTM_Cell_2 = REAL(ElemData_loc(iVar,:))
            CALL AddToElemData(ElementOut,'TTM_T_e[temp]'                      ,RealArray=TTM_Cell_2(1:PP_nElems))

          ELSEIF (STRICMP(VarNamesElemData_loc(iVar),"TTM_T_i[md_temp]")) THEN
            nRestartVars=nRestartVars+1
            ALLOCATE( TTM_Cell_3(1:PP_nElems) )
            TTM_Cell_3 = REAL(ElemData_loc(iVar,:))
            CALL AddToElemData(ElementOut,'TTM_T_i[md_temp]'                   ,RealArray=TTM_Cell_3(1:PP_nElems))

          ELSEIF (STRICMP(VarNamesElemData_loc(iVar),"TTM_[xi]")) THEN
            nRestartVars=nRestartVars+1
            ALLOCATE( TTM_Cell_4(1:PP_nElems) )
            TTM_Cell_4 = REAL(ElemData_loc(iVar,:))
            CALL AddToElemData(ElementOut,'TTM_[xi]'                           ,RealArray=TTM_Cell_4(1:PP_nElems))

          ELSEIF (STRICMP(VarNamesElemData_loc(iVar),"TTM_[source]")) THEN
            nRestartVars=nRestartVars+1
            ALLOCATE( TTM_Cell_5(1:PP_nElems) )
            TTM_Cell_5 = REAL(ElemData_loc(iVar,:))
            CALL AddToElemData(ElementOut,'TTM_[source]'                       ,RealArray=TTM_Cell_5(1:PP_nElems))

          ELSEIF (STRICMP(VarNamesElemData_loc(iVar),"TTM_[v_com.x]")) THEN
            nRestartVars=nRestartVars+1
            ALLOCATE( TTM_Cell_6(1:PP_nElems) )
            TTM_Cell_6 = REAL(ElemData_loc(iVar,:))
            CALL AddToElemData(ElementOut,'TTM_[v_com.x]'                       ,RealArray=TTM_Cell_6(1:PP_nElems))

          ELSEIF (STRICMP(VarNamesElemData_loc(iVar),"TTM_[v_com.y]")) THEN
            nRestartVars=nRestartVars+1
            ALLOCATE( TTM_Cell_7(1:PP_nElems) )
            TTM_Cell_7 = REAL(ElemData_loc(iVar,:))
            CALL AddToElemData(ElementOut,'TTM_[v_com.y]'                       ,RealArray=TTM_Cell_7(1:PP_nElems))

          ELSEIF (STRICMP(VarNamesElemData_loc(iVar),"TTM_[v_com.z]")) THEN
            nRestartVars=nRestartVars+1
            ALLOCATE( TTM_Cell_8(1:PP_nElems) )
            TTM_Cell_8 = REAL(ElemData_loc(iVar,:))
            CALL AddToElemData(ElementOut,'TTM_[v_com.z]'                       ,RealArray=TTM_Cell_8(1:PP_nElems))

          ELSEIF (STRICMP(VarNamesElemData_loc(iVar),"TTM_[fd_k]")) THEN
            nRestartVars=nRestartVars+1
            ALLOCATE( TTM_Cell_9(1:PP_nElems) )
            TTM_Cell_9 = REAL(ElemData_loc(iVar,:))
            CALL AddToElemData(ElementOut,'TTM_[fd_k]'                          ,RealArray=TTM_Cell_9(1:PP_nElems))

          ELSEIF (STRICMP(VarNamesElemData_loc(iVar),"TTM_[fd_g]")) THEN
            nRestartVars=nRestartVars+1
            ALLOCATE( TTM_Cell_10(1:PP_nElems) )
            TTM_Cell_10 = REAL(ElemData_loc(iVar,:))
            CALL AddToElemData(ElementOut,'TTM_[fd_g]'                          ,RealArray=TTM_Cell_10(1:PP_nElems))

          ELSEIF (STRICMP(VarNamesElemData_loc(iVar),"TTM_charge[Z]")) THEN
            nRestartVars=nRestartVars+1
            ALLOCATE( TTM_Cell_11(1:PP_nElems) )
            TTM_Cell_11 = REAL(ElemData_loc(iVar,:))
            CALL AddToElemData(ElementOut,'TTM_charge[Z]'                       ,RealArray=TTM_Cell_11(1:PP_nElems))

          ELSEIF (STRICMP(VarNamesElemData_loc(iVar),"TTM_n_e(ElectronDensity)")) THEN
            nRestartVars=nRestartVars+1
            ALLOCATE( TTM_Cell_12(1:PP_nElems) )
            TTM_Cell_12 = REAL(ElemData_loc(iVar,:))
            CALL AddToElemData(ElementOut,'TTM_n_e(ElectronDensity)'            ,RealArray=TTM_Cell_12(1:PP_nElems))

          ELSEIF (STRICMP(VarNamesElemData_loc(iVar),"TTM_omega_pe_cold(PlasmaFrequency)")) THEN
            nRestartVars=nRestartVars+1
            ALLOCATE( TTM_Cell_13(1:PP_nElems) )
            TTM_Cell_13 = REAL(ElemData_loc(iVar,:))
            CALL AddToElemData(ElementOut,'TTM_omega_pe_cold(PlasmaFrequency)'  ,RealArray=TTM_Cell_13(1:PP_nElems))

          ELSEIF (STRICMP(VarNamesElemData_loc(iVar),"TTM_omega_pe_warm(PlasmaFrequency)")) THEN
            nRestartVars=nRestartVars+1
            ALLOCATE( TTM_Cell_14(1:PP_nElems) )
            TTM_Cell_14 = REAL(ElemData_loc(iVar,:))
            CALL AddToElemData(ElementOut,'TTM_omega_pe_warm(PlasmaFrequency)'  ,RealArray=TTM_Cell_14(1:PP_nElems))

          ELSEIF (STRICMP(VarNamesElemData_loc(iVar),"TTM_dt_PIC_cold(TimeStep)")) THEN
            nRestartVars=nRestartVars+1
            ALLOCATE( TTM_Cell_15(1:PP_nElems) )
            TTM_Cell_15 = REAL(ElemData_loc(iVar,:))
            CALL AddToElemData(ElementOut,'TTM_dt_PIC_cold(TimeStep)'           ,RealArray=TTM_Cell_15(1:PP_nElems))

          ELSEIF (STRICMP(VarNamesElemData_loc(iVar),"TTM_dt_PIC_warm(TimeStep)")) THEN
            nRestartVars=nRestartVars+1
            ALLOCATE( TTM_Cell_16(1:PP_nElems) )
            TTM_Cell_16 = REAL(ElemData_loc(iVar,:))
            CALL AddToElemData(ElementOut,'TTM_dt_PIC_warm(TimeStep)'           ,RealArray=TTM_Cell_16(1:PP_nElems))

          ELSEIF (STRICMP(VarNamesElemData_loc(iVar),"TTM_T_e(ElectronTempInKelvin)")) THEN
            nRestartVars=nRestartVars+1
            ALLOCATE( TTM_Cell_17(1:PP_nElems) )
            TTM_Cell_17 = REAL(ElemData_loc(iVar,:))
            CALL AddToElemData(ElementOut,'TTM_T_e(ElectronTempInKelvin)'       ,RealArray=TTM_Cell_17(1:PP_nElems))

          ELSEIF (STRICMP(VarNamesElemData_loc(iVar),"TTM_lambda_D(DebyeLength)")) THEN
            nRestartVars=nRestartVars+1
            ALLOCATE( TTM_Cell_18(1:PP_nElems) )
            TTM_Cell_18 = REAL(ElemData_loc(iVar,:))
            CALL AddToElemData(ElementOut,'TTM_lambda_D(DebyeLength)'           ,RealArray=TTM_Cell_18(1:PP_nElems))
          END IF
        END DO
        DEALLOCATE(ElemData_loc)
        DEALLOCATE(VarNamesElemData_loc)
        DEALLOCATE(tmp)
      END IF
      LBWRITE(UNIT_stdOut,'(A,I5,A,I5,A)')' Found ',nRestartVars,' of ',18,' TTM elem data arrays in restart file.'
    ELSE
      CALL abort(&
          __STAMP__&
          ,' Cannot read TTM cell data in parallel on restart. Use only one processor!')

    END IF ! MPIRoot
  END IF ! DoRestart.EQV..FALSE.

  ! get min/max plasma frequency
  LBWRITE(UNIT_StdOut,'(A)')' TTM - Plasma Frequency'
  omega_pe_min(1:2)=HUGE(1.)
  omega_pe_max(1:2)=0.0
  DO iElem=1,PP_nElems
    IF(TTM_Cell_13(iElem).GT.0.0)THEN!time step for cold electrons
      omega_pe_min(1)=MIN(omega_pe_min(1),TTM_Cell_13(iElem))
      omega_pe_max(1)=MAX(omega_pe_max(1),TTM_Cell_13(iElem))
    END IF
    IF(TTM_Cell_14(iElem).GT.0.0)THEN!time step for warm electrons
      omega_pe_min(2)=MIN(omega_pe_min(2),TTM_Cell_14(iElem))
      omega_pe_max(2)=MAX(omega_pe_max(2),TTM_Cell_14(iElem))
    END IF
  END DO
  CALL PrintOption('TTM cold elec: omega_pe_min','OUTPUT',RealOpt=omega_pe_min(1))
  CALL PrintOption('TTM cold elec: omega_pe_max','OUTPUT',RealOpt=omega_pe_max(1))
  CALL PrintOption('TTM warm elec: omega_pe_min','OUTPUT',RealOpt=omega_pe_min(2))
  CALL PrintOption('TTM warm elec: omega_pe_max','OUTPUT',RealOpt=omega_pe_max(2))

  ! get min/max PIC timestep (0.2 / plasma frequency)
  LBWRITE(UNIT_StdOut,'(A)')' TTM - PIC Time Step approximation'
  dt_PIC_min(1:2)=HUGE(1.)
  dt_PIC_max(1:2)=0.0
  DO iElem=1,PP_nElems
    IF(TTM_Cell_15(iElem).GT.0.0)THEN!time step for cold electrons
      dt_PIC_min(1)=MIN(dt_PIC_min(1),TTM_Cell_15(iElem))
      dt_PIC_max(1)=MAX(dt_PIC_max(1),TTM_Cell_15(iElem))
    END IF
    IF(TTM_Cell_16(iElem).GT.0.0)THEN!time step for warm electrons
      dt_PIC_min(2)=MIN(dt_PIC_min(2),TTM_Cell_16(iElem))
      dt_PIC_max(2)=MAX(dt_PIC_max(2),TTM_Cell_16(iElem))
    END IF
  END DO
  CALL PrintOption('TTM cold elec: dt_PIC_min','OUTPUT',RealOpt=dt_PIC_min(1))
  CALL PrintOption('TTM cold elec: dt_PIC_max','OUTPUT',RealOpt=dt_PIC_max(1))
  CALL PrintOption('TTM warm elec: dt_PIC_min','OUTPUT',RealOpt=dt_PIC_min(2))
  CALL PrintOption('TTM warm elec: dt_PIC_max','OUTPUT',RealOpt=dt_PIC_max(2))

  ! get min/max debye length
  LBWRITE(UNIT_StdOut,'(A)')' TTM - Debye length'
  lambda_D_min=HUGE(1.)
  lambda_D_max=0.0
  DO iElem=1,PP_nElems
    IF(TTM_Cell_18(iElem).GT.0.0)THEN!time step for cold electrons
      lambda_D_min=MIN(lambda_D_min,TTM_Cell_18(iElem))
      lambda_D_max=MAX(lambda_D_max,TTM_Cell_18(iElem))
    END IF
  END DO
  CALL PrintOption('TTM Debye length: lambda_D_min','OUTPUT',RealOpt=lambda_D_min)
  CALL PrintOption('TTM Debye length: lambda_D_max','OUTPUT',RealOpt=lambda_D_max)

  ! Fill .csv file for analysis
  CALL WriteTTMdataToCSV()

END IF !DoImportTTMFile.EQV..TRUE.


TTMInitIsDone=.TRUE.
LBWRITE(UNIT_stdOut,'(A)')' INIT TTM DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitTTM


SUBROUTINE GetFDGridCellSize(fd_hx,fd_hy,fd_hz)
!----------------------------------------------------------------------------------------------------------------------------------!
! search for the following line and extract the values for fd_hx, fd_hy and fd_hz
! the info is extracted from the IMD+TTM output- or log-file
! fd_h.x:24.107143, fd_h.y:22.509134,fd_h.z:22.509134\nnVolume of one FD cell: 1.221415e+04 [cubic Angstroms]
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_TTM_Vars
USE MOD_Particle_Vars ,ONLY: IMDLengthScale
USE MOD_ReadInTools   ,ONLY: PrintOption
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(INOUT)   :: fd_hx,fd_hy,fd_hz
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: ioUnit
INTEGER              :: IndNum
INTEGER              :: iLine
INTEGER              :: io_error
CHARACTER(255)       :: StrTmp
!===================================================================================================================================
fd_hx = -1.0
fd_hy = -1.0
fd_hz = -1.0

OPEN(NEWUNIT=ioUnit,FILE=TRIM(TTMLogFile),STATUS='OLD',ACTION='READ',IOSTAT=io_error)
IF(io_error.NE.0)THEN
  CALL abort(&
  __STAMP__&
  ,'Cannot open specified File (ttm data) from '//TRIM(TTMLogFile))
END IF
SWRITE(UNIT_stdOut,'(A,A)') " Reading TTM log info from file (TTMLogFile): ",TRIM(TTMLogFile)
iLine=1
DO !iLine=1,1 ! header lines: currently 1 -> adjust?
  READ(ioUnit,'(A)',IOSTAT=io_error)StrTmp
  IF(io_error.GT.0)THEN
    SWRITE(UNIT_stdOut,'(A,I10)') "io_error=",io_error
    SWRITE(UNIT_stdOut,'(A,I10)') "A problem reading the TTM file occured in line: ",iLine+1
    SWRITE(UNIT_stdOut,'(A)') StrTmp
    CALL abort(&
    __STAMP__&
    ,'Error reading specified File (ttm data): '//TRIM(TTMFile))
  ELSE IF(io_error.LT.0)THEN
    SWRITE(UNIT_stdOut,'(A,I10)') "End of file reached: ",iLine
    EXIT
  ELSE ! io_error = 0
    iLine=iLine+1
    IndNum=INDEX(TRIM(StrTmp),'fd_h.x:')
    IF(IndNum.GT.0)THEN
      StrTmp=TRIM(StrTmp(IndNum:LEN(StrTmp)))
      IndNum=INDEX(TRIM(StrTmp),',')
      IF(IndNum.GT.0)THEN
        CALL str2real(StrTmp(8:IndNum-1),fd_hx,io_error)
        IF(io_error.NE.0)THEN
          SWRITE(UNIT_stdOut,'(A,I10,A)') "io_error=",io_error,", failed to get fd_hx"
          fd_hx=-99.
          Exit
        END IF
        fd_hx=fd_hx*IMDLengthScale
        CALL PrintOption('TTM: fd_hx','OUTPUT',RealOpt=fd_hx)
        StrTmp=TRIM(StrTmp(IndNum+1:LEN(StrTmp)))
      END IF

      IndNum=INDEX(TRIM(StrTmp),'fd_h.y:')
      IF(IndNum.GT.0)THEN
        StrTmp=TRIM(StrTmp(IndNum:LEN(StrTmp)))
        IndNum=INDEX(TRIM(StrTmp),',')
        IF(IndNum.GT.0)THEN
          CALL str2real(StrTmp(8:IndNum-1),fd_hy,io_error)
          IF(io_error.NE.0)THEN
            SWRITE(UNIT_stdOut,'(A,I10,A)') "io_error=",io_error,", failed to get fd_hy"
            fd_hy=-99.
            Exit
          END IF
          fd_hy=fd_hy*IMDLengthScale
          CALL PrintOption('TTM: fd_hy','OUTPUT',RealOpt=fd_hy)
          StrTmp=TRIM(StrTmp(IndNum+1:LEN(StrTmp)))
        END IF
      END IF

      IndNum=INDEX(TRIM(StrTmp),'fd_h.z:')
      IF(IndNum.GT.0)THEN
        StrTmp=TRIM(StrTmp(IndNum:LEN(StrTmp)))
        IndNum=INDEX(TRIM(StrTmp),'\')
        IF(IndNum.GT.0)THEN
          CALL str2real(StrTmp(8:IndNum-1),fd_hz,io_error)
        ELSE
          CALL str2real(StrTmp(8:LEN(StrTmp)),fd_hz,io_error)
        END IF
        IF(io_error.NE.0)THEN
          SWRITE(UNIT_stdOut,'(A,I10,A)') "io_error=",io_error,", failed to get fd_hz"
          fd_hz=-99.
          Exit
        END IF
        fd_hz=fd_hz*IMDLengthScale
        CALL PrintOption('TTM: fd_hz','OUTPUT',RealOpt=fd_hz)
        StrTmp=TRIM(StrTmp(IndNum+1:LEN(StrTmp)))
      END IF

      EXIT
    END IF
    IF(iLine.GT.100)EXIT ! only read the first 100 lines
  END IF
END DO



CLOSE(ioUnit)
END SUBROUTINE GetFDGridCellSize


SUBROUTINE InitIMD_TTM_Coupling()
!----------------------------------------------------------------------------------------------------------------------------------!
! 1.) assign charges to each atom using the averaged FD cell charge within each FD cell and sum the charge for step 2.)
! 2.) reconstruct the electron phase space using MD and TTM data and the summed charged per FD cell for which an electron is
!     created to achieve an ionization degree of unity (fully ionized plasma)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Globals_Vars        ,ONLY: BoltzmannConst, ElementaryCharge
USE MOD_Preproc
USE MOD_TTM_Vars
USE MOD_Particle_Vars       ,ONLY: PDM,PEM,PartState,nSpecies,Species,PartSpecies,IMDSpeciesCharge,IMDSpeciesID
USE MOD_Mesh_Vars           ,ONLY: NGeo,XCL_NGeo,XiCL_NGeo,wBaryCL_NGeo
USE MOD_DSMC_Vars           ,ONLY: CollisMode,DSMC,PartStateIntEn
USE MOD_part_emission_tools ,ONLY: CalcVelocity_maxwell_lpn
USE MOD_DSMC_Vars           ,ONLY: useDSMC
USE MOD_Eval_xyz            ,ONLY: TensorProductInterpolation
USE MOD_Part_Tools          ,ONLY: GetNextFreePosition
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: ChargeLower,ChargeUpper
INTEGER :: ElemCharge(1:PP_nElems)
INTEGER :: ElecSpecIndx,iSpec,location,iElem,iPart,ParticleIndexNbr
REAL    :: ChargeProbability
REAL    :: iRan
REAL    :: PartPosRef(1:3)
REAL    :: CellElectronTemperature
REAL    :: MaxElectronTemp_eV
!===================================================================================================================================

! ---------------------------------------------------------------------------------------------------------------------------------
! 1.) reconstruct ions and determine charge
! ---------------------------------------------------------------------------------------------------------------------------------
SWRITE(UNIT_stdOut,*)'1.) Reconstructing ions and determining charge'

! Initialize the element charge with zero
ElemCharge(1:PP_nElems)=0

! Loop over all particles in the vector list
DO iPart=1,PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart)) THEN

    ! If the current particle is part of the IMD species list, then a charge can be assigned
    IF(ANY(PartSpecies(iPart).EQ.IMDSpeciesID(:)))THEN !

      ! Check if the average charge (TTM data) for the current cell is zero
      ! (if yes: add a charge to force at least singly charged ions)
      IF(TTM_Cell_11(PEM%GlobalElemID(iPart)).EQ.0)THEN ! this would create uncharged atoms -> force singly charged ions
        ElemCharge(PEM%GlobalElemID(iPart))=ElemCharge(PEM%GlobalElemID(iPart))+1
        location = MINLOC(ABS(IMDSpeciesCharge-1),1) !Determines the location of the element in the array with min value
      ELSE
        ! Get the TTM cell charge avergae value and select and upper and lower charge number
        ChargeLower       = INT(TTM_Cell_11(PEM%GlobalElemID(iPart))) ! FLOOR(): use first DOF (0,0,0) because the data is const. in each cell
        ChargeUpper       = ChargeLower+1                        ! CEILING(): floor + 1
        ChargeProbability = REAL(ChargeUpper)-TTM_Cell_11(PEM%GlobalElemID(iPart)) ! 2-1,4=0.6 -> 60% probability to get lower charge

        ! Distribute the charge using random numbers
        CALL RANDOM_NUMBER(iRan)

        ! Compare the random number with the charge difference
        IF(iRan.LT.ChargeProbability)THEN ! Select the lower charge number
          ! Determines the location of the element in the array with min value: get the index of the corresponding charged ion
          ! species
          location                       = MINLOC(ABS(IMDSpeciesCharge-ChargeLower),1)
          ElemCharge(PEM%GlobalElemID(iPart)) = ElemCharge(PEM%GlobalElemID(iPart))+ChargeLower
        ELSE ! Select the upper charge number
          ! Determines the location of the element in the array with min value: get the index of the corresponding charged ion
          ! species
          location                       = MINLOC(ABS(IMDSpeciesCharge-ChargeUpper),1)
          ElemCharge(PEM%GlobalElemID(iPart)) = ElemCharge(PEM%GlobalElemID(iPart))+ChargeUpper
        END IF
      END IF

      ! Set the species ID to atom/singly charged ion/doubly charged ... and so on
      PartSpecies(iPart)=IMDSpeciesID(location)
    END IF
  END IF
END DO

! ---------------------------------------------------------------------------------------------------------------------------------
! 2.) reconstruct electrons
! ---------------------------------------------------------------------------------------------------------------------------------
SWRITE(UNIT_stdOut,*)'2.) Reconstructing electrons'

! Set the electron temperature (in eV) from the TTM data
MaxElectronTemp_eV=MAXVAL(TTM_Cell_2(:))

! Initialize the species index for the electron species with -1
ElecSpecIndx = -1

! Loop over all species and find the index corresponding to the electron species
DO iSpec = 1, nSpecies
  IF (Species(iSpec)%ChargeIC.GT.0.0) CYCLE
  IF(NINT(Species(iSpec)%ChargeIC/(-ElementaryCharge)).EQ.1) THEN
    ElecSpecIndx = iSpec
    EXIT
  END IF
END DO
IF (ElecSpecIndx.EQ.-1) CALL abort(&
  __STAMP__&
  ,'Electron species not found. Cannot create electrons without the defined species!')

! Loop over all elements and the sum of charges in each element (for each charge assigned in an element, an electron is created)
DO iElem=1,PP_nElems
  DO iPart=1,ElemCharge(iElem) ! 1 electron for each charge of each element

    ParticleIndexNbr                     = GetNextFreePosition()

    !Set new SpeciesID of new particle (electron)
    PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
    PDM%isNewPart(ParticleIndexNbr)      = .TRUE.
    PartSpecies(ParticleIndexNbr)        = ElecSpecIndx

    ! Place the electron randomly in the reference cell
    CALL RANDOM_NUMBER(PartPosRef(1:3)) ! get random reference space
    PartPosRef(1:3)=PartPosRef(1:3)*2. - 1. ! map (0,1) -> (-1,1)
    CALL TensorProductInterpolation(PartPosRef(1:3),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,iElem) &
                      ,PartState(1:3,ParticleIndexNbr)) !Map into phys. space

    ! Set the internal energies (vb, rot and electronic) to zero if needed
    IF ((useDSMC).AND.(CollisMode.GT.1)) THEN
      PartStateIntEn( 1,ParticleIndexNbr) = 0.
      PartStateIntEn( 2,ParticleIndexNbr) = 0.
      IF (DSMC%ElectronicModel.GT.0)  PartStateIntEn( 3,ParticleIndexNbr) = 0.
    END IF

    ! Set the element ID of the electron to the current element ID
    PEM%GlobalElemID(ParticleIndexNbr) = iElem

    ! Determine the electron temperature (for each cell different) in Kelvin for the Maxwellian distribution
    IF(TTM_Cell_2(iElem).LE.0.0)THEN ! not enough atoms in FD cell for averaging a temperature: use max value for electrons
      CellElectronTemperature=(MaxElectronTemp_eV*ElementaryCharge/BoltzmannConst) ! convert eV to K: 1 [eV] = e/kB [K]
    ELSE
      CellElectronTemperature=(TTM_Cell_2(iElem)*ElementaryCharge/BoltzmannConst) ! convert eV to K: 1 [eV] = e/kB [K]
    END IF

    ! Set the electron velocity using the Maxwellian distribution (use the function that is suitable for small numbers)
    CALL CalcVelocity_maxwell_lpn(ElecSpecIndx, PartState(4:6,ParticleIndexNbr),&
                                  Temperature=CellElectronTemperature)
  END DO
END DO


END SUBROUTINE InitIMD_TTM_Coupling


!----------------------------------------------------------------------------------------------------------------------------------!
!> Write TTM info to InitialTTMdata.csv file
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE WriteTTMdataToCSV()
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Preproc
USE MOD_Globals   ,ONLY: MPIRoot,FILEEXISTS,unit_stdout
USE MOD_Globals   ,ONLY: abort
USE MOD_Mesh_Vars ,ONLY: ElemBaryNGeo
USE MOD_TTM_Vars  ,ONLY: TTM_Cell_1,TTM_Cell_2,TTM_Cell_3,TTM_Cell_4,TTM_Cell_5,TTM_Cell_6,TTM_Cell_7,TTM_Cell_8
USE MOD_TTM_Vars  ,ONLY: TTM_Cell_9,TTM_Cell_10,TTM_Cell_11,TTM_Cell_12,TTM_Cell_13,TTM_Cell_14,TTM_Cell_15
USE MOD_TTM_Vars  ,ONLY: TTM_Cell_16,TTM_Cell_17,TTM_Cell_18
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!LOGICAL,INTENT(IN)                  :: WriteHeader
!REAL,INTENT(IN),OPTIONAL            :: time
!INTEGER(KIND=8),INTENT(IN),OPTIONAL :: iter
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=22),PARAMETER              :: outfile='InitialTTMdata.csv'
INTEGER                                  :: ioUnit,I
CHARACTER(LEN=250)                        :: formatStr
INTEGER,PARAMETER                        :: nOutputVar=21

CHARACTER(LEN=255),DIMENSION(nOutputVar) :: StrVarNames(nOutputVar)=(/ CHARACTER(LEN=255) :: &
          'TTM_N[natoms]'           , &
          'TTM_T_e[temp]'           , &
          'TTM_T_i[md_temp]'        , &
          'TTM_[xi]'                , &
          'TTM_[source]'            , &
          'TTM_[v_com.x]'           , &
          'TTM_[v_com.y]'           , &
          'TTM_[v_com.z]'           , &
          'TTM_[fd_k]'              , &
          'TTM_[fd_g]'              , &
          'TTM_charge[Z]'           , &
          'TTM_n_e(ElectronDensity)', &
          'TTM_omega_pe_cold(PlasmaFrequency)', &
          'TTM_omega_pe_warm(PlasmaFrequency)', &
          'TTM_dt_PIC_cold(TimeStep)',          &
          'TTM_dt_PIC_warm(TimeStep)',          &
          'TTM_T_e(ElectronTempInKelvin)',      &
          'TTM_lambda_D(DebyeLength)',          &
           'x',&
           'y',&
           'z' /)

CHARACTER(LEN=255),DIMENSION(nOutputVar) :: tmpStr ! needed because PerformAnalyze is called mutiple times at the beginning
CHARACTER(LEN=1000)                      :: tmpStr2
INTEGER                                  :: iElem
CHARACTER(LEN=1),PARAMETER               :: delimiter=","
!===================================================================================================================================
IF(.NOT.MPIRoot)RETURN

! write the header line
OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfile),STATUS="UNKNOWN")
tmpStr=""
DO I=1,nOutputVar
  WRITE(tmpStr(I),'(A)')delimiter//'"'//TRIM(StrVarNames(I))//'"'
END DO
WRITE(formatStr,'(A1)')'('
DO I=1,nOutputVar
  IF(I.EQ.nOutputVar)THEN ! skip writing "," and the end of the line
    WRITE(formatStr,'(A,A1,I2)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(I))
  ELSE
    WRITE(formatStr,'(A,A1,I2,A1)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(I)),','
  END IF
END DO
WRITE(formatStr,'(A,A1)')TRIM(formatStr),')' ! finish the format
WRITE(tmpStr2,formatStr)tmpStr               ! use the format and write the header names to a temporary string
tmpStr2(1:1) = " "
WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2))    ! clip away the front and rear white spaces of the temporary string
CLOSE(ioUnit)

! write the data lines for each element
IF(FILEEXISTS(outfile))THEN
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfile),POSITION="APPEND",STATUS="OLD")

  WRITE(formatStr,'(A2,I2,A14,A1)')'(',nOutputVar,CSVFORMAT,')'
  DO iElem=1,PP_nElems
    WRITE(tmpStr2,formatStr)&
        " ",TTM_Cell_1(iElem),&
        delimiter,TTM_Cell_2(iElem),&
        delimiter,TTM_Cell_3(iElem),&
        delimiter,TTM_Cell_4(iElem),&
        delimiter,TTM_Cell_5(iElem),&
        delimiter,TTM_Cell_6(iElem),&
        delimiter,TTM_Cell_7(iElem),&
        delimiter,TTM_Cell_8(iElem),&
        delimiter,TTM_Cell_9(iElem),&
        delimiter,TTM_Cell_10(iElem),&
        delimiter,TTM_Cell_11(iElem),&
        delimiter,TTM_Cell_12(iElem),&
        delimiter,TTM_Cell_13(iElem),&
        delimiter,TTM_Cell_14(iElem),&
        delimiter,TTM_Cell_15(iElem),&
        delimiter,TTM_Cell_16(iElem),&
        delimiter,TTM_Cell_17(iElem),&
        delimiter,TTM_Cell_18(iElem),&
        delimiter,ElemBaryNGeo(1,iElem),&
        delimiter,ElemBaryNGeo(2,iElem),&
        delimiter,ElemBaryNGeo(3,iElem)
    WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2)) ! clip away the front and rear white spaces of the data line
  END DO
  CLOSE(ioUnit)
ELSE
  SWRITE(UNIT_StdOut,'(A)')"InitialTTMdata.csv does not exist. Cannot write TTM info!"
END IF
!END IF
END SUBROUTINE WriteTTMdataToCSV


SUBROUTINE FinalizeTTM()
!----------------------------------------------------------------------------------------------------------------------------------!
! finalize TTM variables
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_TTM_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(TTM_FD)
SDEALLOCATE(TTM_Cell_1)
SDEALLOCATE(TTM_Cell_2)
SDEALLOCATE(TTM_Cell_3)
SDEALLOCATE(TTM_Cell_4)
SDEALLOCATE(TTM_Cell_5)
SDEALLOCATE(TTM_Cell_6)
SDEALLOCATE(TTM_Cell_7)
SDEALLOCATE(TTM_Cell_8)
SDEALLOCATE(TTM_Cell_9)
SDEALLOCATE(TTM_Cell_10)
SDEALLOCATE(TTM_Cell_11)
SDEALLOCATE(TTM_Cell_12)
SDEALLOCATE(TTM_Cell_13)
SDEALLOCATE(TTM_Cell_14)
SDEALLOCATE(TTM_Cell_15)
SDEALLOCATE(TTM_Cell_16)
SDEALLOCATE(TTM_Cell_17)
SDEALLOCATE(TTM_Cell_18)
SDEALLOCATE(ElemBaryFD)
SDEALLOCATE(ElemIndexFD)
SDEALLOCATE(ElemIsDone)
END SUBROUTINE FinalizeTTM
#endif /*PARTICLES*/

END MODULE MOD_TTMInit
