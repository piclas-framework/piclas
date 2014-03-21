#include "boltzplatz.h"

MODULE MOD_DSMC_ElectronicModel
!===================================================================================================================================
! module including qk procedures
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
PUBLIC :: ElectronicEnergyExchange, InitElectronShell, TVEEnergyExchange
PUBLIC :: ReadSpeciesLevel
!===================================================================================================================================
CONTAINS


SUBROUTINE InitElectronShell(iSpecies,iPart)
!===================================================================================================================================
! init electronic shell
!===================================================================================================================================
USE MOD_DSMC_Vars,              ONLY : DSMC, SpecDSMC, PartStateIntEn
USE MOD_Particle_Vars,          ONLY : BoltzmannConst
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE                                                                                    
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPart, iSpecies                                                 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iQua                                                        
REAL                          :: iRan, ElectronicPartition, ElectronicPartitionTemp, iRan2      
!===================================================================================================================================

ElectronicPartition  = 0
ElectronicPartitionTemp = 0
! calculate sum over all energy levels == qartition function for temperature Telec
DO iQua = 0, SpecDSMC(iSpecies)%MaxElecQuant - 1
  ElectronicPartitionTemp = SpecDSMC(iSpecies)%ElectronicState(1,iQua) * &
          exp ( - SpecDSMC(iSpecies)%ElectronicState(2,iQua) / SpecDSMC(iSpecies)%Telec)
  IF ( ElectronicPartitionTemp .gt. ElectronicPartition ) THEN
    ElectronicPartition = ElectronicPartitionTemp
  END IF
END DO
ElectronicPartitionTemp = 0
! select level
CALL RANDOM_NUMBER(iRan2)
DO WHILE ( iRan2 .ge. ElectronicPartitionTemp / ElectronicPartition )
  CALL RANDOM_NUMBER(iRan)
  iQua = int( ( SpecDSMC(iSpecies)%MaxElecQuant ) * iRan)
  ElectronicPartitionTemp  = SpecDSMC(iSpecies)%ElectronicState(1,iQua) * &
           exp ( - SpecDSMC(iSpecies)%ElectronicState(2,iQua) / SpecDSMC(iSpecies)%Telec)
  CALL RANDOM_NUMBER(iRan2)
END DO
#if ( PP_TimeDiscMethod == 42 )
  SpecDSMC(iSpecies)%levelcounter(iQua) = SpecDSMC(iSpecies)%levelcounter(iQua) + 1
#endif
PartStateIntEn(iPart,3) = BoltzmannConst * SpecDSMC(iSpecies)%ElectronicState(2,iQua)

END SUBROUTINE InitElectronShell


SUBROUTINE ElectronicEnergyExchange(CollisionEnergy,iPart1,FakXi,iPart2,iElem)
!===================================================================================================================================
! electronic energy exchange
!===================================================================================================================================
USE MOD_DSMC_Vars,              ONLY : Coll_pData, DSMC, SpecDSMC, PartStateIntEn
USE MOD_Particle_Vars,          ONLY : PartSpecies, BoltzmannConst, GEO, usevMPF,PartMPF
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE                                                                                    
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPart1
INTEGER, INTENT(IN), OPTIONAL :: iPart2,iElem
REAL, INTENT(IN)              :: FakXi
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(INOUT)           :: CollisionEnergy                                                !
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iQuaMax, MaxElecQuant, iQua, iQuaold, iQuaMax2, iQuaMax3
REAL                          :: iRan, iRan2, gmax, gtemp, PartStateTemp
! vMPF
REAL                          :: DeltaPartStateIntEn, Phi, PartStateIntEnTemp
!===================================================================================================================================

! Determine max electronic quant
MaxElecQuant = SpecDSMC(PartSpecies(iPart1))%MaxElecQuant - 1
#if ( PP_TimeDiscMethod==42 )
! determine old Quant
DO iQua = 0, MaxElecQuant
  IF ( PartStateIntEn(iPart1,3) / BoltzmannConst .ge. &
    SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua) ) THEN
    iQuaold = iQua
  ELSE
  ! exit loop
    EXIT
  END IF
END DO
#endif
! determine maximal Quant and term according to Eq (7) of Liechty
gmax = 0
PartStateTemp = CollisionEnergy / BoltzmannConst
DO iQua = 0, MaxElecQuant
  IF ( PartStateTemp - &
           SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua) .ge. 0 ) THEN
    gtemp = SpecDSMC(PartSpecies(iPart1))%ElectronicState(1,iQua) * &
            ( CollisionEnergy - BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua))**FakXi
    ! maximal possible Quant before term goes negative
    iQuaMax = iQua
    IF ( gtemp .gt. gmax ) THEN
    ! Quant of largest value of Eq (7)
      gmax = gtemp
      iQuaMax2 = iQua
    END IF
  ELSE
  ! exit loop
!    EXIT
    CYCLE
  END IF
END DO
! max value for denominator == max. propability
!  iQuaMax3 = min(iQuaMax,iQuaMax2)
!  gmax = SpecDSMC(PartSpecies(iPart1))%ElectronicState(1,iQuaMax3) * &
!       ( CollisionEnergy - BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQuaMax3))**FakXi
! max iQuant for dicing
iQuaMax  = max(iQuaMax,iQuaMax2)
CALL RANDOM_NUMBER(iRan)
iQua = int( ( iQuaMax +1 ) * iRan)
gtemp = SpecDSMC(PartSpecies(iPart1))%ElectronicState(1,iQua) * &
        ( CollisionEnergy - BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua))**FakXi
CALL RANDOM_NUMBER(iRan2)
! acceptance-rejection for iQuaElec
DO WHILE ( iRan2 .ge. gtemp / gmax )
  CALL RANDOM_NUMBER(iRan)
  iQua = int( ( iQuaMax +1 ) * iRan)
  gtemp = SpecDSMC(PartSpecies(iPart1))%ElectronicState(1,iQua) * &
          ( CollisionEnergy - BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua))**FakXi
  CALL RANDOM_NUMBER(iRan2)
END DO

IF (usevMPF) THEN
  IF (PartMPF( iPart1).GT.PartMPF( iPart2)) THEN
!      DeltaPartStateIntEn = 0.0
    Phi = PartMPF( iPart2) / PartMPF( iPart1)
    PartStateIntEnTemp = BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua)
    CollisionEnergy = CollisionEnergy - PartStateIntEnTemp
    PartStateIntEnTemp = (DBLE(1)-Phi) * PartStateIntEn( iPart1,3) + Phi * PartStateIntEnTemp
    PartStateTemp = PartStateIntEnTemp / BoltzmannConst
    ! searche for new vib quant
    iQuaMax = 0
    DO iQua = 0, MaxElecQuant
      IF ( PartStateTemp .ge. &
        SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua) ) THEN
        iQuaMax = iQua
      ELSE
      ! exit loop
        EXIT
      END IF
    END DO
    iQua = iQuaMax
    PartStateIntEn( iPart1,3) = BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua)
    DeltaPartStateIntEn = PartMPF( iPart1) &
                        * (PartStateIntEnTemp - PartStateIntEn( iPart1,3))
!      CollisionEnergy = CollisionEnergy + DeltaPartStateIntEn / PartMPF(PairE_vMPF(2)) 
    GEO%DeltaEvMPF(iElem) = GEO%DeltaEvMPF(iElem) + DeltaPartStateIntEn
  END IF
ELSE
#if (PP_TimeDiscMethod==42)
! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
  IF ( DSMC%ReservoirSimuRate .EQV. .FALSE. ) THEN
# endif
   PartStateIntEn(iPart1,3) = BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua)
   CollisionEnergy = CollisionEnergy - PartStateIntEn(iPart1,3)
#if (PP_TimeDiscMethod==42)
  END IF
# endif
END IF

#if ( PP_TimeDiscMethod ==42 )
    ! list of number of particles in each energy level
IF ( DSMC%ReservoirSimuRate .EQV. .FALSE. ) THEN
  SpecDSMC(PartSpecies(iPart1))%levelcounter(iQuaold) = SpecDSMC(PartSpecies(iPart1))%levelcounter(iQuaold) - 1
  SpecDSMC(PartSpecies(iPart1))%levelcounter(iQua)    = SpecDSMC(PartSpecies(iPart1))%levelcounter(iQua)    + 1
  SpecDSMC(PartSpecies(iPart1))%dtlevelcounter(iQua)  = SpecDSMC(PartSpecies(iPart1))%dtlevelcounter(iQua)  + 1
END IF
! collision with X resulting in a transition from i to j
IF ( present(iPart2) .AND. (.NOT.usevMPF) ) THEN
SpecDSMC(PartSpecies(iPart1))%ElectronicTransition(PartSpecies(iPart2),iQuaold,iQua) = &
                              SpecDSMC(PartSpecies(iPart1))%ElectronicTransition(PartSpecies(iPart2),iQuaold,iQua) + 1
END IF
#endif
END SUBROUTINE ElectronicEnergyExchange

SUBROUTINE TVEEnergyExchange(CollisionEnergy,iPart1,FakXi,iPart2,iElem)
!===================================================================================================================================
! electronic energy exchange
!===================================================================================================================================
USE MOD_DSMC_Vars,              ONLY : Coll_pData, DSMC, SpecDSMC, PartStateIntEn
USE MOD_Particle_Vars,          ONLY : PartSpecies, BoltzmannConst, GEO, usevMPF,PartMPF
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE                                                                                    
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPart1
INTEGER, INTENT(IN), OPTIONAL :: iPart2,iElem
REAL, INTENT(IN)              :: FakXi
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(INOUT)           :: CollisionEnergy                                                !
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iQuaMax, MaxElecQuant, iQua, iQuaold, iQuaMax2, iQuaMax3
INTEGER                       :: jQVib, QMaxVib
REAL                          :: iRan, iRan2, gmax, gtemp, PartStateTemp, iRanVib
! vMPF
REAL                          :: DeltaPartStateIntEn, Phi, PartStateIntEnTemp
!===================================================================================================================================

! Determine max electronic quant
MaxElecQuant = SpecDSMC(PartSpecies(iPart1))%MaxElecQuant - 1
#if ( PP_TimeDiscMethod==42 )
! determine old Quant
DO iQua = 0, MaxElecQuant
  IF ( PartStateIntEn(iPart1,3) / BoltzmannConst .ge. &
    SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua) ) THEN
    iQuaold = iQua
  ELSE
  ! exit loop
    EXIT
  END IF
END DO
#endif
! determine maximal Quant and term according to Eq (7) of Liechty
gmax = 0
PartStateTemp = CollisionEnergy / BoltzmannConst
DO iQua = 0, MaxElecQuant
  IF ( (PartStateTemp  &
           - SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua) &
           - DSMC%GammaQuant * SpecDSMC(PartSpecies(iPart1))%CharaTVib) &
      .ge. 0 ) THEN
    gtemp = SpecDSMC(PartSpecies(iPart1))%ElectronicState(1,iQua) * &
            ( CollisionEnergy - BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua) &
            -DSMC%GammaQuant * SpecDSMC(PartSpecies(iPart1))%CharaTVib * BoltzmannConst)**FakXi
    ! maximal possible Quant before term goes negative
    iQuaMax = iQua
    IF ( gtemp .gt. gmax ) THEN
    ! Quant of largest value of Eq (7)
      gmax = gtemp
      iQuaMax2 = iQua
    END IF
  ELSE
  ! exit loop
  !  EXIT
    CYCLE
  END IF
END DO
! max value for denominator == max. propability
! iQuaMax3 = min(iQuaMax,iQuaMax2)
! gmax = SpecDSMC(PartSpecies(iPart1))%ElectronicState(1,iQuaMax3) * &
!       ( CollisionEnergy - BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQuaMax3) &
!         -DSMC%GammaQuant * SpecDSMC(PartSpecies(iPart1))%CharaTVib * BoltzmannConst)**FakXi
! max iQuant for dicing
iQuaMax  = max(iQuaMax,iQuaMax2)
QMaxVib = CollisionEnergy/(BoltzmannConst*SpecDSMC(PartSpecies(iPart1))%CharaTVib)  &
            - DSMC%GammaQuant
QMaxVib = MIN(INT(QMaxVib) + 1, SpecDSMC(PartSpecies(iPart1))%MaxVibQuant)
CALL RANDOM_NUMBER(iRan)
CALL RANDOM_NUMBER(iRanVib)
iQua = int( ( iQuaMax +1 ) * iRan)
jQVib =  INT(iRanVib * QMaxVib)
!  gtemp = SpecDSMC(PartSpecies(iPart1))%ElectronicState(1,iQua) * &
!          ( CollisionEnergy - BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua) &
!          -(DSMC%GammaQuant + jQVib) * SpecDSMC(PartSpecies(iPart1))%CharaTVib * BoltzmannConst)**FakXi
gtemp =( CollisionEnergy - BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua) &
      -(DSMC%GammaQuant + jQVib) * SpecDSMC(PartSpecies(iPart1))%CharaTVib * BoltzmannConst)
IF (gtemp.LE.0.0) THEN
  gtemp = 0.0
ELSE
  gtemp = SpecDSMC(PartSpecies(iPart1))%ElectronicState(1,iQua) *(gtemp)**FakXi
END IF
CALL RANDOM_NUMBER(iRan2)
! acceptance-rejection for iQuaElec
DO WHILE ( iRan2 .ge. gtemp / gmax )
  CALL RANDOM_NUMBER(iRan)
  CALL RANDOM_NUMBER(iRanVib)
  iQua = int( ( iQuaMax +1 ) * iRan)
  jQVib =  INT(iRanVib * QMaxVib)
  gtemp =( CollisionEnergy - BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua) &
        -(DSMC%GammaQuant + jQVib) * SpecDSMC(PartSpecies(iPart1))%CharaTVib * BoltzmannConst)
  IF (gtemp.LE.0.0) THEN
    gtemp = 0.0
  ELSE
    gtemp = SpecDSMC(PartSpecies(iPart1))%ElectronicState(1,iQua) *(gtemp)**FakXi
  END IF
  CALL RANDOM_NUMBER(iRan2)
END DO

!vmpf muss noch gemacht werden !!!!
IF (usevMPF) THEN
  IF (PartMPF( iPart1).GT.PartMPF( iPart2)) THEN
!      DeltaPartStateIntEn = 0.0
    Phi = PartMPF( iPart2) / PartMPF( iPart1)
    PartStateIntEnTemp = BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua)
    CollisionEnergy = CollisionEnergy - PartStateIntEnTemp
    PartStateIntEnTemp = (DBLE(1)-Phi) * PartStateIntEn( iPart1,3) + Phi * PartStateIntEnTemp
    PartStateTemp = PartStateIntEnTemp / BoltzmannConst
    ! searche for new vib quant
    iQuaMax = 0
    DO iQua = 0, MaxElecQuant
      IF ( PartStateTemp .ge. &
        SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua) ) THEN
        iQuaMax = iQua
      ELSE
      ! exit loop
        EXIT
      END IF
    END DO
    iQua = iQuaMax
    PartStateIntEn( iPart1,3) = BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua)
    DeltaPartStateIntEn = PartMPF( iPart1) &
                        * (PartStateIntEnTemp - PartStateIntEn( iPart1,3))
!      CollisionEnergy = CollisionEnergy + DeltaPartStateIntEn / PartMPF(PairE_vMPF(2)) 
    GEO%DeltaEvMPF(iElem) = GEO%DeltaEvMPF(iElem) + DeltaPartStateIntEn
  END IF
ELSE
#if (PP_TimeDiscMethod==42)
! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
  IF ( DSMC%ReservoirSimuRate .EQV. .FALSE. ) THEN
# endif
   PartStateIntEn(iPart1,3) = BoltzmannConst * SpecDSMC(PartSpecies(iPart1))%ElectronicState(2,iQua)
   PartStateIntEn(iPart1,1) = (jQVib + DSMC%GammaQuant) * BoltzmannConst &
                  * SpecDSMC(PartSpecies(iPart1))%CharaTVib
   CollisionEnergy = CollisionEnergy - PartStateIntEn(iPart1,3) - PartStateIntEn(iPart1,1)
#if (PP_TimeDiscMethod==42)
  END IF
# endif
END IF

#if ( PP_TimeDiscMethod ==42 )
    ! list of number of particles in each energy level
IF ( DSMC%ReservoirSimuRate .EQV. .FALSE. ) THEN
  SpecDSMC(PartSpecies(iPart1))%levelcounter(iQuaold) = SpecDSMC(PartSpecies(iPart1))%levelcounter(iQuaold) - 1
  SpecDSMC(PartSpecies(iPart1))%levelcounter(iQua)    = SpecDSMC(PartSpecies(iPart1))%levelcounter(iQua)    + 1
  SpecDSMC(PartSpecies(iPart1))%dtlevelcounter(iQua)  = SpecDSMC(PartSpecies(iPart1))%dtlevelcounter(iQua)  + 1
END IF
! collision with X resulting in a transition from i to j
IF ( present(iPart2) .AND. (.NOT.usevMPF) ) THEN
SpecDSMC(PartSpecies(iPart1))%ElectronicTransition(PartSpecies(iPart2),iQuaold,iQua) = &
                              SpecDSMC(PartSpecies(iPart1))%ElectronicTransition(PartSpecies(iPart2),iQuaold,iQua) + 1
END IF
#endif
END SUBROUTINE TVEEnergyExchange

SUBROUTINE ReadSpeciesLevel ( Dsetname, iSpec )
!===================================================================================================================================
! Subroutine to read the electronic levels from DSMCSpeciesElectronicState.h5
!===================================================================================================================================
! use module
USE MOD_io_hdf5
USE MOD_Globals
USE MOD_DSMC_Vars,                    ONLY: DSMC, SpecDSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                    :: iSpec
CHARACTER(LEN=64),INTENT(IN)                          :: dsetname
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                               :: err
! HDF5 specifier taken from extractParticles
INTEGER(HSIZE_T), DIMENSION(2)                        :: dims,sizeMax
INTEGER(HID_T)                                        :: file_id_dsmc                       ! File identifier
INTEGER(HID_T)                                        :: dset_id_dsmc                       ! Dataset identifier
INTEGER(HID_T)                                        :: dataspace                          ! Dataspace identifier
INTEGER(HID_T)                                        :: filespace                          ! filespace identifier
!===================================================================================================================================

! Initialize FORTRAN interface.
CALL h5open_f(err)
! Open the file.
CALL h5fopen_f (DSMC%ElectronicStateDatabase, H5F_ACC_RDONLY_F, file_id_dsmc, err)
! Open the  dataset.
CALL h5dopen_f(file_id_dsmc, dsetname, dset_id_dsmc, err)
! Get the file space of the dataset.
CALL H5DGET_SPACE_F(dset_id_dsmc, FileSpace, err)
! get size
CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, dims, SizeMax, err)
! Allocate electronic_state
ALLOCATE ( SpecDSMC(iSpec)%ElectronicState( 1:dims(1), 0:dims(2)-1 ) )
! read data
CALL H5dread_f(dset_id_dsmc, H5T_NATIVE_DOUBLE, SpecDSMC(iSpec)%ElectronicState, dims, err)
! Close the file.
CALL h5fclose_f(file_id_dsmc, err)
! Close FORTRAN interface.
CALL h5close_f(err)

SWRITE(*,*) 'Read entry ',dsetname,' from ',DSMC%ElectronicStateDatabase

END SUBROUTINE ReadSpeciesLevel

END MODULE MOD_DSMC_ElectronicModel
