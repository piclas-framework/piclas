#include "boltzplatz.h"

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
! IMD Two Temperature Model (TTM)
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
SWRITE(UNIT_stdOut,'(A)') ' INIT TTM ...'

DoImportTTMFile   =GETLOGICAL('DoImportTTMFile','.FALSE.')

IF(DoRestart.EQV..FALSE.)THEN ! if no restart is performed: read the TTM field data from *.ttm file
  TTMFile   =TRIM(GETSTR('TTMFile',''))
  TTMLogFile=TRIM(GETSTR('TTMLogFile',''))
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
      ElemIndexFD=0.0
      ALLOCATE( ElemIsDone(FD_nElems) )
      ElemIsDone=.FALSE.
      ! IMD TTM data format and target arrays for data input (from ASCII data file *.ttm)
      !        1      2    3       4  5      6       7       8       9    10   11       ->  TTM_FD(1:11,ix,iy,iz,iLineTTM) -> TTM_DG
      !  1 2 3 4      5    6       7  8      9       10      11      12   13   14 15    ->  tmp_array(1:15)
      ! #x y z natoms temp md_temp xi source v_com.x v_com.y v_com.z fd_k fd_g Z  proc
      SWRITE(UNIT_stdOut,'(A,A)') "Reading from file: ",TRIM(TTMFile)
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
  
      read(StrTmp,*,iostat=io_error)  i
      TTMNumber = i
      SWRITE(UNIT_stdOut,'(A,A)')  " IMD *.ttm file = ",StrTmp
      SWRITE(UNIT_stdOut,'(A,I5)') " TTMNumber",TTMNumber
  
      iLineTTM=0
      DO i=1,1 ! header lines: currently 1 -> adjust?
        READ(ioUnit,'(A)',IOSTAT=io_error)StrTmp
      END DO
      DO
        READ(ioUnit,*,IOSTAT=io_error) tmp_array(1:15)
        IF(io_error>0)THEN
          SWRITE(UNIT_stdOut,'(A,I5)') " io_error=",io_error
          SWRITE(UNIT_stdOut,'(A,I5)') " A problem reading the TTM file occured in line: ",iLineTTM+2
          SWRITE(UNIT_stdOut,'(A1)') ' ['
          DO j=1,15
            SWRITE(UNIT_stdOut,'(A1,E25.14E3)',ADVANCE='NO') " ",tmp_array(j)
          END DO
          SWRITE(UNIT_stdOut,'(A1)') ']'
          CALL abort(&
          __STAMP__&
          ,'Error reading specified File (ttm data): '//TRIM(TTMFile))
        ELSE IF(io_error<0)THEN
          SWRITE(UNIT_stdOut,'(A,I5)') "End of file reached: ",iLineTTM+1
          EXIT
        ELSE ! io_error = 0
        END IF
        iLineTTM=iLineTTM+1 ! *.ttm data file: line counter
  
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
          SWRITE(UNIT_stdOut,'(A,I5)') "A problem reading the TTM file occured in line: ",iLineTTM+1
          DoImportTTMFile=.FALSE.
          EXIT
        END IF
      END DO
      IF(FD_nElems.NE.ix*iy*iz)THEN
        CALL abort(&
        __STAMP__&
        ,'Number of FD elements FD_nElems does not comply with the number of FD elements read from *.ttm file'//TRIM(TTMFile))
      END IF
      CLOSE(ioUnit)
      SWRITE(UNIT_stdOut,'(A,I5,A)') "Read ",iLineTTM," lines from file (+1 header line)"
  
      IF(FD_nElems.EQ.PP_nElems)THEN ! same number of elements in FD grid and DG solution -> assume they coincide
        ! the local DG solution in physical and reference space
        ALLOCATE( TTM_DG(1:11,0:PP_N,0:PP_N,0:PP_N,PP_nElems) )
        TTM_DG=0.0 ! initialize
        DO iElem=1,PP_nElems
          DO iELemFD=1,FD_nElems
            IF(ElemIsDone(iElemFD))CYCLE ! the element has already been found
            IF(AlmostEqualToTolerance(ElemBaryNGeo(1,iElem),ElemBaryFD(1,iElemFD),TTMElemBaryTolerance).EQV..FALSE.) CYCLE
            IF(AlmostEqualToTolerance(ElemBaryNGeo(2,iElem),ElemBaryFD(2,iElemFD),TTMElemBaryTolerance).EQV..FALSE.) CYCLE
            IF(AlmostEqualToTolerance(ElemBaryNGeo(3,iElem),ElemBaryFD(3,iElemFD),TTMElemBaryTolerance).EQV..FALSE.) CYCLE
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
        SWRITE(UNIT_stdOut,'(A,I5)') "FD_nElems=",FD_nElems
        SWRITE(UNIT_stdOut,'(A,I5)') "PP_nElems=",PP_nElems
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


TTMInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT TTM DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitTTM


SUBROUTINE GetFDGridCellSize(fd_hx,fd_hy,fd_hz)
!----------------------------------------------------------------------------------------------------------------------------------!
! search for the following line and extract "fd_hx,fd_hy,fd_hz"
! fd_h.x:24.107143, fd_h.y:22.509134,fd_h.z:22.509134\nnVolume of one FD cell: 1.221415e+04 [cubic Angstroms] 
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_TTM_Vars
USE MOD_Particle_Vars, ONLY:IMDLengthScale
!USE MOD_RegressionCheck_tools, ONLY: str2real
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(INOUT)   :: fd_hx,fd_hy,fd_hz
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER        :: ioUnit
INTEGER        :: IndNum,IndNum2
INTEGER        :: i
INTEGER        :: io_error
CHARACTER(255) :: StrTmp,StrTmp2
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
SWRITE(UNIT_stdOut,'(A,A)') "Reading from file: ",TRIM(TTMLogFile)
i=1
DO !i=1,1 ! header lines: currently 1 -> adjust?
  READ(ioUnit,'(A)',IOSTAT=io_error)StrTmp
  IF(io_error.GT.0)THEN
    SWRITE(UNIT_stdOut,'(A,I5)') "io_error=",io_error
    SWRITE(UNIT_stdOut,'(A,I5)') "A problem reading the TTM file occured in line: ",i+1
    SWRITE(UNIT_stdOut,'(A)') StrTmp
    CALL abort(&
    __STAMP__&
    ,'Error reading specified File (ttm data): '//TRIM(TTMFile))
  ELSE IF(io_error.LT.0)THEN
    SWRITE(UNIT_stdOut,'(A,I5)') "End of file reached: ",i
    EXIT
  ELSE ! io_error = 0
    i=i+1
    IndNum=INDEX(TRIM(StrTmp),'fd_h.x:')
    IF(IndNum.GT.0)THEN
      StrTmp=TRIM(StrTmp(IndNum:LEN(StrTmp)))
      IndNum=INDEX(TRIM(StrTmp),',')
      IF(IndNum.GT.0)THEN
        CALL str2real(StrTmp(8:IndNum-1),fd_hx,io_error)
        IF(io_error.NE.0)THEN
          SWRITE(UNIT_stdOut,'(A,I5,A)') "io_error=",io_error,", failed to get fd_hx"
          fd_hx=-99.
          Exit
        END IF
        fd_hx=fd_hx*IMDLengthScale
        SWRITE(UNIT_stdOut,'(A6,E25.14E3)') "fd_hx=",fd_hx
        StrTmp=TRIM(StrTmp(IndNum+1:LEN(StrTmp)))
      END IF

      IndNum=INDEX(TRIM(StrTmp),'fd_h.y:')
      IF(IndNum.GT.0)THEN
        StrTmp=TRIM(StrTmp(IndNum:LEN(StrTmp)))
        IndNum=INDEX(TRIM(StrTmp),',')
        IF(IndNum.GT.0)THEN
          CALL str2real(StrTmp(8:IndNum-1),fd_hy,io_error)
          IF(io_error.NE.0)THEN
            SWRITE(UNIT_stdOut,'(A,I5,A)') "io_error=",io_error,", failed to get fd_hy"
            fd_hy=-99.
            Exit
          END IF
          fd_hy=fd_hy*IMDLengthScale
          SWRITE(UNIT_stdOut,'(A6,E25.14E3)') "fd_hy=",fd_hy
          StrTmp=TRIM(StrTmp(IndNum+1:LEN(StrTmp)))
        END IF
      END IF

      IndNum=INDEX(TRIM(StrTmp),'fd_h.z:')
      IF(IndNum.GT.0)THEN
        StrTmp=TRIM(StrTmp(IndNum:LEN(StrTmp)))
          IndNum=INDEX(TRIM(StrTmp),'\')
          IF(IndNum.GT.0)THEN
            CALL str2real(StrTmp(8:IndNum-1),fd_hz,io_error)
            IF(io_error.NE.0)THEN
              SWRITE(UNIT_stdOut,'(A,I5,A)') "io_error=",io_error,", failed to get fd_hz"
              fd_hz=-99.
              Exit
            END IF
            fd_hz=fd_hz*IMDLengthScale
          SWRITE(UNIT_stdOut,'(A6,E25.14E3)') "fd_hz=",fd_hz
            StrTmp=TRIM(StrTmp(IndNum+1:LEN(StrTmp)))
        END IF
      END IF

      EXIT
    END IF
    IF(i.GT.100)EXIT ! only read the first 100 lines
  END IF
END DO



CLOSE(ioUnit)
END SUBROUTINE GetFDGridCellSize


SUBROUTINE InitIMD_TTM_Coupling() 
!----------------------------------------------------------------------------------------------------------------------------------!
! finalize TTM variables
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
INTEGER :: ChargeLower,ChargeUpper    ! 
INTEGER :: ElemCharge(PP_nElems)      !
INTEGER :: ElecSpecIndx,iSpec,location,iElem,iPart,ParticleIndexNbr
REAL    :: ChargeProbability          ! 
REAL    :: iRan, RandVal(2)           !
REAL    :: PartPosRef(1:3)
REAL    :: CellElectronTemperature
!===================================================================================================================================
ElemCharge=0
DO iPart=1,PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart)) THEN
    IF(ANY(PartSpecies(iPart).EQ.IMDSpeciesID(:)))THEN ! 
      ! get the TTM cell charge avergae value and select and upper and lower charge number
      ChargeLower=INT(TTM_DG(11,0,0,0,PEM%Element(iPart)))
      ChargeUpper=ChargeLower+1
      ChargeProbability=REAL(ChargeUpper)-TTM_DG(11,0,0,0,PEM%Element(iPart)) ! 2 - 1,4 = 0.6 -> 60% probability to get lower charge
      ! distribute the charge using random numbers
      CALL RANDOM_NUMBER(iRan)
      IF(iRan.LT.ChargeProbability)THEN ! select the lower charge number
        location = MINLOC(ABS(IMDSpeciesCharge-ChargeLower),1)
        ElemCharge(PEM%Element(iPart))=ElemCharge(PEM%Element(iPart))+ChargeLower
      ELSE ! select the upper charge number
        location = MINLOC(ABS(IMDSpeciesCharge-ChargeUpper),1)
        ElemCharge(PEM%Element(iPart))=ElemCharge(PEM%Element(iPart))+ChargeUpper
      END IF
      PartSpecies(iPart)=IMDSpeciesID(location)
    END IF
  END IF
END DO

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
    CellElectronTemperature=(TTM_DG(2,0,0,0,iElem)*ElectronCharge/BoltzmannConst)
    CALL CalcVelocity_maxwell_lpn(ElecSpecIndx, PartState(ParticleIndexNbr,4:6),&
                                  Temperature=CellElectronTemperature)
  END DO
END DO



!CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3),ElemID)
!IF(MAXVAL(ABS(PartPosRef(1:3))).LT.1.0) THEN ! particle is inside 
!  PEM%Element(iPart)=ElemID
!END IF

!  IF(DoRefMapping)THEN
!    CALL ParticleRefTracking() ! newton
!
!  ELSE
!    CALL ParticleTracing() ! schnittpunkt
!  CALL PartInElemCheck(LastPartPos(iPart,1:3),iPart,ElemID,isHit)
!  END IF

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


END MODULE MOD_TTMInit
