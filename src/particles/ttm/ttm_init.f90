#include "boltzplatz.h"

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

CONTAINS

SUBROUTINE InitTTM()
!===================================================================================================================================
! IMD Two Temperature Model (TTM):
! Read TTM data from ASCII file and create FD mesh which is then mapped onto the TTM-DG solution for output
! IMD TTM data format and target arrays for data input (from ASCII data file *.ttm) are shown below
! -------------------------------------------------------------------------------------------------------------------------------
! The properties and array index are chosen as follows
!        1      2    3       4  5      6       7       8       9    10   11       ->  TTM_FD(1:11,ix,iy,iz) -> TTM_DG
!  1 2 3 4      5    6       7  8      9       10      11      12   13   14 15    ->  tmp_array(1:15)
! #x y z natoms temp md_temp xi source v_com.x v_com.y v_com.z fd_k fd_g Z  proc
! -------------------------------------------------------------------------------------------------------------------------------
! read the following fields ...
! 
! 1.)  x:         0 to Nx-1 FD cells in x-direction                          [-]                not stored
! 2.)  y:         0 to Ny-1 FD cells in y-direction                          [-]                not stored
! 3.)  z:         0 to Nz-1 FD cells in z-direction                          [-]                not stored
! 4.)  natoms:    number of atoms in the FD cell                             [-]                stored in TTM_FD(1,.. -> TTM_DG(1,..
! 5.)  temp:      temperature of the electrons (fluid model)                 [eV]               stored in TTM_FD(2,.. -> TTM_DG(2,..
! 6.)  md_temp:   temperature of the atoms (MD simulated discrete atoms)     [eV]               stored in TTM_FD(3,.. -> TTM_DG(3,..
! 7.)  xi: ?                                                                 [?]                stored in TTM_FD(4,.. -> TTM_DG(4,..
! 8.)  source:    energy input source, e.g., laser energy density per volume [?]                stored in TTM_FD(5,.. -> TTM_DG(5,..
! 9.)  v_com.x:   barycentric velocity of all MD atoms in x-direction        [Ångström/10,18fs] stored in TTM_FD(6,.. -> TTM_DG(6,..
! 10.) v_com.y:   barycentric velocity of all MD atoms in y-direction        [Ångström/10,18fs] stored in TTM_FD(7,.. -> TTM_DG(7,..
! 11.) v_com.z:   barycentric velocity of all MD atoms in z-direction        [Ångström/10,18fs] stored in TTM_FD(8,.. -> TTM_DG(8,..
! 12.) fd_k:      heat conduction coefficient                                [?]                stored in TTM_FD(9,.. -> TTM_DG(9,..
! 13.) fd_g:      electron-phonon coupling coefficient                       [?]                stored in TTM_FD(10,. -> TTM_DG(10,.
! 14.) Z:         averged charge of the atoms within the FD cell             [e]                stored in TTM_FD(11,. -> TTM_DG(11,.
! 15.) proc:      rank number of MPI process                                 [-]                not stored
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_TTM_Vars
USE MOD_Restart_Vars,       ONLY:DoRestart
#ifdef MPI
USE MOD_Particle_MPI_Vars,     ONLY:PartMPI
#endif /* MPI*/
USE MOD_Particle_Mesh_Vars,      ONLY:ElemBaryNGeo
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
INTEGER        :: i,j,k,ix,iy,iz
INTEGER        :: io_error
INTEGER        :: iLineTTM
INTEGER        :: iElem,iElemFD
CHARACTER(255) :: StrTmp
REAL           :: tmp_array(1:15)
REAL           :: fd_hx,fd_hy,fd_hz
!===================================================================================================================================
IF(TTMInitIsDone)THEN
   SWRITE(UNIT_stdOut,'(A)') "InitTTM already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT TTM (IMD Two-Temperature Model) ...'

DoImportTTMFile = GETLOGICAL('DoImportTTMFile','.FALSE.')
TTMLogFile      = TRIM(GETSTR('TTMLogFile',''))

IF(DoImportTTMFile.EQV..TRUE.)THEN
  IF(DoRestart.EQV..FALSE.)THEN ! if no restart is performed: read the TTM field data from *.ttm file
    TTMFile   =TRIM(GETSTR('TTMFile',''))
    IF((TRIM(TTMFile).NE.'').AND.(TRIM(TTMLogFile).NE.''))THEN
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
        ALLOCATE( TTM_FD(1:11,TTMGridFDdim(1),TTMGridFDdim(2),TTMGridFDdim(3)) )
        TTM_FD=0.0
        ALLOCATE( ElemBaryFD(3,FD_nElems) )
        ElemBaryFD=0.0
        ALLOCATE( ElemIndexFD(3,FD_nElems) )
        ElemIndexFD=0
        ALLOCATE( ElemIsDone(FD_nElems) )
        ElemIsDone=.FALSE.
        !SWRITE(UNIT_StdOut,'(a3,a30,a3,a33,a3,a7,a3)')' | ',TRIM("Reading TTM data from"),' | ', TRIM(TTMFile),&
                                                      !' | ',TRIM("OUTPUT"),' | '
        SWRITE(UNIT_stdOut,'(A,A)') " Reading from file: ",TRIM(TTMFile)
#ifdef MPI
        IF(.NOT.PartMPI%MPIROOT)THEN
          CALL abort(&
          __STAMP__&
          ,'ERROR: Cannot SetParticlePosition in multi-core environment for SpaceIC=IMD!')
        END IF
#endif /*MPI*/
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
        SWRITE(UNIT_StdOut,'(a3,a30,a3,i33,a3,a7,a3)')' | ',TRIM("TTMNumber")     ,' | ', TTMNumber   ,' | ',TRIM("OUTPUT"),' | '
    
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
              SWRITE(UNIT_stdOut,'(A1,E25.14E3)',ADVANCE='NO') " ",tmp_array(j)
            END DO
            SWRITE(UNIT_stdOut,'(A1)') ']'
            CALL abort(&
            __STAMP__&
            ,'Error reading specified File (ttm data): '//TRIM(TTMFile))
          ELSE IF(io_error<0)THEN
            SWRITE(UNIT_stdOut,'(A,I10,A)') " End of file reached: ",iLineTTM+1," -> EXIT"
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
            SWRITE(UNIT_stdOut,'(A,I10)') " A problem reading the TTM file occured in line: ",iLineTTM+1
            DoImportTTMFile=.FALSE.
            EXIT
          END IF
        END DO
        IF(FD_nElems.NE.ix*iy*iz)THEN
          SWRITE(UNIT_stdOut,'(A)')     ' Error: FD_nElems.NE.ix*iy*iz'
          SWRITE(UNIT_stdOut,'(A,I10)') ' FD_nElems = ',FD_nElems
          SWRITE(UNIT_stdOut,'(A,I10)') ' ix*iy*iz  = ',ix*iy*iz
          CALL abort(&
          __STAMP__&
          ,'Number of FD elements FD_nElems does not comply with the number of FD elements read from *.ttm file ['&
          //TRIM(TTMFile)//']')
        END IF
        CLOSE(ioUnit)
        SWRITE(UNIT_stdOut,'(A,I10,A)') " Read ",iLineTTM," lines from file (+1 header line)"
        SWRITE(UNIT_StdOut,'(a3,a30,a3,I33,a3,a7,a3)')' | ',TRIM("IMD TTM FD grid points"),&
                                                           ' | ',iLineTTM,' | ',TRIM("OUTPUT"),' | '
    
        IF(FD_nElems.EQ.PP_nElems)THEN ! same number of elements in FD grid and DG solution -> assume they coincide
          SWRITE(UNIT_stdOut,'(A)') ' Searching for all FD cells within the DG mesh in order to map the values ...'
          ! the local DG solution in physical and reference space
          ALLOCATE( TTM_DG(1:11,0:PP_N,0:PP_N,0:PP_N,PP_nElems) )
          TTM_DG=0.0 ! initialize
          DO iElem=1,PP_nElems
            DO iELemFD=1,FD_nElems
              IF(ElemIsDone(iElemFD))CYCLE ! the element has already been found
              IF(ALMOSTEQUALRELATIVE(ElemBaryNGeo(1,iElem),ElemBaryFD(1,iElemFD),TTMElemBaryTolerance).EQV..FALSE.) CYCLE
              IF(ALMOSTEQUALRELATIVE(ElemBaryNGeo(2,iElem),ElemBaryFD(2,iElemFD),TTMElemBaryTolerance).EQV..FALSE.) CYCLE
              IF(ALMOSTEQUALRELATIVE(ElemBaryNGeo(3,iElem),ElemBaryFD(3,iElemFD),TTMElemBaryTolerance).EQV..FALSE.) CYCLE
              ElemIsDone(iElemFD)=.TRUE.
              DO i=0,PP_N ! set all DG DOF values equal to FD cell value
                DO j=0,PP_N ! set all DG DOF values equal to FD cell value
                  DO k=0,PP_N ! set all DG DOF values equal to FD cell value
                    TTM_DG(1:11,i,j,k,iElem)=TTM_FD(1:11,ElemIndexFD(1,iElemFD),ElemIndexFD(2,iElemFD),ElemIndexFD(3,iElemFD))
                  END DO
                END DO
              END DO
            END DO
          END DO
          SWRITE(UNIT_stdOut,'(A)') ' Done.'
          IF(ANY(ElemIsDone).EQV..FALSE.)THEN
            SWRITE(UNIT_stdOut,'(A)') " NOT all IMD TTM FD cells have been located within the DG grid!"
          ELSE
            SWRITE(UNIT_stdOut,'(A)') " All IMD TTM FD cells have been located within the DG grid!"
          END IF
          SDEALLOCATE(TTM_FD)
          SDEALLOCATE(ElemBaryFD)
          SDEALLOCATE(ElemIndexFD)
          SDEALLOCATE(ElemIsDone)
        ELSE ! use swap-mesh or interpolate the FD grid to the DG solution
          SWRITE(UNIT_stdOut,'(A,I10)') "FD_nElems=",FD_nElems
          SWRITE(UNIT_stdOut,'(A,I10)') "PP_nElems=",PP_nElems
          CALL abort(&
          __STAMP__&
          ,'Error interpolating ttm data (FD grid) to DG solution. FD_nElems.NE.PP_nElems. Different elements not implemented')
        END IF
      END IF
    ELSE
      DoImportTTMFile=.FALSE.
    END IF
  ELSE
    SWRITE(UNIT_stdOut,'(A)')' INIT TTM: data will be read from restart file!'
  END IF
END IF


TTMInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT TTM DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

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
USE MOD_Particle_Vars, ONLY:IMDLengthScale
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
SWRITE(UNIT_stdOut,'(A,A)') " Reading from file: ",TRIM(TTMLogFile)
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
        SWRITE(UNIT_StdOut,'(a3,a30,a3,E33.14E3,a3,a7,a3)')' | ',TRIM("fd_hx")     ,' | ',fd_hx,' | ',TRIM("OUTPUT"),' | '
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
          SWRITE(UNIT_StdOut,'(a3,a30,a3,E33.14E3,a3,a7,a3)')' | ',TRIM("fd_hy")     ,' | ',fd_hy,' | ',TRIM("OUTPUT"),' | '
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
        SWRITE(UNIT_StdOut,'(a3,a30,a3,E33.14E3,a3,a7,a3)')' | ',TRIM("fd_hz")     ,' | ',fd_hz,' | ',TRIM("OUTPUT"),' | '
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
USE MOD_Globals_Vars,   ONLY:ElectronCharge
USE MOD_PreProc
USE MOD_TTM_Vars
USE MOD_Particle_Vars, ONLY:PDM,PEM,BoltzmannConst,PartState,nSpecies,Species,PartSpecies,IMDSpeciesCharge,IMDSpeciesID
USE MOD_Eval_xyz,      ONLY:eval_xyz_elemcheck
USE MOD_Mesh_Vars,     ONLY:NGeo,XCL_NGeo,XiCL_NGeo,wBaryCL_NGeo
USE MOD_DSMC_Vars,     ONLY:CollisMode,DSMC,PartStateIntEn
USE MOD_part_emission, ONLY:CalcVelocity_maxwell_lpn
USE MOD_DSMC_Vars,     ONLY:useDSMC
USE MOD_Eval_xyz,      ONLY:Eval_XYZ_Poly
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
! 1.) reconstruct ions and determine charge
ElemCharge(1:PP_nElems)=0
DO iPart=1,PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart)) THEN
    IF(ANY(PartSpecies(iPart).EQ.IMDSpeciesID(:)))THEN ! 
      IF(TTM_DG(11,0,0,0,PEM%Element(iPart)).EQ.0)THEN ! this would create uncharged atoms -> force singly charged ions
        ElemCharge(PEM%Element(iPart))=ElemCharge(PEM%Element(iPart))+1
        location = MINLOC(ABS(IMDSpeciesCharge-1),1) !Determines the location of the element in the array with min value
      ELSE
        ! get the TTM cell charge avergae value and select and upper and lower charge number
        ChargeLower       = INT(TTM_DG(11,0,0,0,PEM%Element(iPart))) ! use first DOF (0,0,0) because the data is const. in each cell
        ChargeUpper       = ChargeLower+1
        ChargeProbability = REAL(ChargeUpper)-TTM_DG(11,0,0,0,PEM%Element(iPart)) ! 2-1,4=0.6 -> 60% probability to get lower charge
        ! distribute the charge using random numbers
        CALL RANDOM_NUMBER(iRan)
        IF(iRan.LT.ChargeProbability)THEN ! select the lower charge number
          location = MINLOC(ABS(IMDSpeciesCharge-ChargeLower),1) !Determines the location of the element in the array with min value
          ElemCharge(PEM%Element(iPart))=ElemCharge(PEM%Element(iPart))+ChargeLower
        ELSE ! select the upper charge number
          location = MINLOC(ABS(IMDSpeciesCharge-ChargeUpper),1) !Determines the location of the element in the array with min value
          ElemCharge(PEM%Element(iPart))=ElemCharge(PEM%Element(iPart))+ChargeUpper
        END IF
      END IF
      PartSpecies(iPart)=IMDSpeciesID(location) ! set the species ID to atom/singly charged ion/doubly charged ... and so on
    END IF
  END IF
END DO

! 2.) reconstruct electrons
MaxElectronTemp_eV=MAXVAL(TTM_DG(2,0,0,0,:))
ElecSpecIndx = -1
DO iSpec = 1, nSpecies
  IF (Species(iSpec)%ChargeIC.GT.0.0) CYCLE
  IF(NINT(Species(iSpec)%ChargeIC/(-1.60217653E-19)).EQ.1) THEN
    ElecSpecIndx = iSpec
    EXIT
  END IF
END DO
IF (ElecSpecIndx.EQ.-1) CALL abort(&
  __STAMP__&
  ,'Electron species not found. Cannot create electrons without the defined species!')

DO iElem=1,PP_nElems
  DO iPart=1,ElemCharge(iElem) ! 1 electron for each charge of each element
    PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + 1
    ParticleIndexNbr = PDM%nextFreePosition(PDM%CurrentNextFreePosition)
    PDM%ParticleVecLength = PDM%ParticleVecLength + 1
    !Set new Species of new particle
    PDM%ParticleInside(ParticleIndexNbr) = .true.
    PartSpecies(ParticleIndexNbr) = ElecSpecIndx
     
    CALL RANDOM_NUMBER(PartPosRef(1:3)) ! get random reference space
    PartPosRef(1:3)=PartPosRef(1:3)*2. - 1. ! map (0,1) -> (-1,1)
    CALL Eval_xyz_Poly(PartPosRef(1:3),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,iElem) &
                      ,PartState(ParticleIndexNbr,1:3)) !Map into phys. space

    IF ((useDSMC).AND.(CollisMode.GT.1)) THEN
      PartStateIntEn(ParticleIndexNbr, 1) = 0.
      PartStateIntEn(ParticleIndexNbr, 2) = 0.
      IF ( DSMC%ElectronicModel )  PartStateIntEn(ParticleIndexNbr, 3) = 0.
    END IF
    PEM%Element(ParticleIndexNbr) = iElem
    IF(TTM_DG(2,0,0,0,iElem).LE.0.0)THEN ! not enough atoms in FD cell for averaging a temperature: use max value for electrons
      CellElectronTemperature=(MaxElectronTemp_eV*ElectronCharge/BoltzmannConst) ! convert eV to K: 1 [eV] = e/kB [K]
    ELSE
      CellElectronTemperature=(TTM_DG(2,0,0,0,iElem)*ElectronCharge/BoltzmannConst) ! convert eV to K: 1 [eV] = e/kB [K]
    END IF
    CALL CalcVelocity_maxwell_lpn(ElecSpecIndx, PartState(ParticleIndexNbr,4:6),&
                                  Temperature=CellElectronTemperature)
  END DO
END DO


END SUBROUTINE InitIMD_TTM_Coupling


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
SDEALLOCATE(TTM_DG)
SDEALLOCATE(ElemBaryFD)
SDEALLOCATE(ElemIndexFD)
SDEALLOCATE(ElemIsDone)
END SUBROUTINE FinalizeTTM
#endif /*PARTICLES*/

END MODULE MOD_TTMInit
