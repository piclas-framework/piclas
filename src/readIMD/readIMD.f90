#ifdef MPI
#ifdef PARTICLES
#include "piclas.h"
module mod_readIMD

implicit none
private

public :: read_IMD_results
public :: DefineParametersReadIMDdata, initReadIMDdata
contains

! ==============================================================================

subroutine DefineParametersReadIMDdata()

  use MOD_ReadInTools ,only: prms
  ! -----------------------------------------------------------------------------
  call prms%CreateLogicalOption('useIMDresults', 'Use Data from IMD', '.FALSE.')
  call prms%CreateStringOption('filenameIMDresults', 'Filename of the IMD result file', 'None')
  call prms%CreateLogicalOption('killPIClas', 'Debugging Variable to kill IMD, only for testing', '.FALSE.')

end subroutine DefineParametersReadIMDdata

! ==============================================================================

subroutine initReadIMDdata()

  use mod_readIMD_vars
  use mod_ReadInTools, only: getLogical, getStr
  ! -----------------------------------------------------------------------------

  useIMDresults = getLogical( 'useIMDresults' )
  filenameIMDresults = getStr( 'filenameIMDresults' )
  killPIClas = getLogical( 'killPIClas' )

end subroutine initReadIMDdata

! ==============================================================================

subroutine read_IMD_results()

  use mod_readIMD_vars
  use mod_globals
  use MOD_Particle_Vars,only:PartState,PDM
  use MOD_DSMC_Vars,only:PartStateIntEn
  use MOD_Particle_Mesh,only:SingleParticleToExactElementNoMap
  use mod_particle_mpi,only:IRecvNbofParticles, SendNbOfParticles, MPIParticleSend, MPIParticleRecv
  use mod_hdf5_output,only:WriteStateToHDF5
  use mod_mesh_vars,only:meshfile
  use mod_part_tools,only:UpdateNextFreePosition

  implicit none
  ! --------------------------------------------------------
  ! local variable
  integer                                   :: mpiInfo
  integer                                   :: filehandle
  integer                                   :: mpiFileError
  real,dimension(3)                         :: Box_X=.0_8, Box_Y=.0_8, Box_Z=.0_8
  integer(kind=MPI_OFFSET_KIND)             :: disp, myFileOffset
  integer(kind=2)                           :: disp_tmp
  integer(kind=2)                           :: observables_tmp
  integer                                   :: observables
  integer(kind=8)                           :: nGlobalAtoms
  integer                                   :: ioStatus(MPI_STATUS_SIZE)
  integer(kind=8)                           :: nAtoms, iAtom, iProc
  integer(kind=8),dimension(0:nProcessors)  :: FileOffsets
  integer(kind=MPI_OFFSET_KIND)             :: myOffset
  character(len=1),dimension(:),allocatable :: AtomsBuffer
  integer(kind=4)                           :: atomBufferSize
  integer                                   :: atomsBufferPos = 0
  integer                                   :: errorLen
  character(len=254)                        :: errorString
  integer(kind=4)                           :: iPart
  real                                      :: MaxX,MaxX_glob
  real                                      :: MinX,MinX_glob
  real                                      :: MaxY,MaxY_glob
  real                                      :: MinY,MinY_glob
  real                                      :: MaxZ,MaxZ_glob
  real                                      :: MinZ,MinZ_glob
  ! -----------------------------------------------------------------------------

  SWRITE(UNIT_stdOut,'(A,A)')'Read IMD-results from file: ',trim(filenameIMDresults)
  call MPI_INFO_CREATE(mpiInfo, iError)

  if( myRank == 0 )then
    call MPI_FILE_OPEN(MPI_COMM_SELF, trim(filenameIMDresults), MPI_MODE_RDONLY,&
                       MPI_INFO_NULL, filehandle, mpiFileError)

    if( mpiFileError /= 0)then
      call MPI_Error_string(mpiFileError, errorString, errorLen, iError)
      write(UNIT_errOut,'(A34,I8,6A)')'Error during MPI_File_open! Rank: ', myRank, '\n',&
                                      'Error-Message: ', trim(errorString), ' "',trim(filenameIMDresults), '"',&
                                      'The program will be terminated!!!'
      call MPI_Abort(MPI_COMM_WORLD, mpiFileError, iError)
    end if

    call MPI_FILE_READ_AT(filehandle,3_8,disp_tmp,2,MPI_BYTE,ioStatus,iError)
    call MPI_FILE_READ_AT(filehandle,5_8,nGlobalAtoms,8,MPI_BYTE,ioStatus,iError)
    call MPI_FILE_READ_AT(filehandle,13_8,observables_tmp,2,MPI_BYTE,ioStatus,iError)

    call MPI_FILE_READ_AT(filehandle,15_8,Box_X,3,MPI_DOUBLE_PRECISION,ioStatus,iError)
    call MPI_FILE_READ_AT(filehandle,39_8,Box_Y,3,MPI_DOUBLE_PRECISION,ioStatus,iError)
    call MPI_FILE_READ_AT(filehandle,63_8,Box_Z,3,MPI_DOUBLE_PRECISION,ioStatus,iError)

    disp = int(disp_tmp,8)
    observables = int(observables_tmp,4)

    SWRITE(UNIT_stdOut,'(3A,I15)')'Number of atmos in ',trim(filenameIMDresults),': ',nGlobalAtoms
    call MPI_FILE_CLOSE(filehandle, mpiFileError)

  end if

  call MPI_BCAST(Box_X, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, iError)
  call MPI_BCAST(Box_Y, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, iError)
  call MPI_BCAST(Box_Z, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, iError)
  call MPI_BCAST(nGlobalAtoms, 8, MPI_BYTE, 0, MPI_COMM_WORLD, iError)
  call MPI_BCAST(observables, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, iError)
  call MPI_BCAST(disp, 8, MPI_BYTE, 0, MPI_COMM_WORLD, iError)

  nAtoms = nGlobalAtoms/nProcessors
  iAtom = nGlobalAtoms - nAtoms * nProcessors
  FileOffsets(0) = 0
  DO iProc=0,nProcessors-1
    FileOffsets(iProc) = nAtoms * iProc + MIN( iProc,iAtom )
  END DO
  FileOffsets(nProcessors) = nGlobalAtoms
  nAtoms = FileOffsets(myRank+1) - FileOffsets(myRank)
  PDM%ParticleVecLength = int ( nAtoms, 4 )
  myOffset = FileOffsets(myRank)

  myFileOffset = disp + myOffset * observables * 8_8
  atomBufferSize = 8 * observables * int ( nAtoms, 4 )
  allocate(AtomsBuffer(atomBufferSize))

  call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(filenameIMDresults), MPI_MODE_RDONLY,&
                      mpiInfo, filehandle, mpiFileError)

  if( mpiFileError /= 0)then
    call MPI_Error_string(mpiFileError, errorString, errorLen, iError)
    write(UNIT_errOut,'(A34,I8,6A)')'Error during MPI_File_open! Rank: ', myRank, '\n',&
                                    'Error-Message: ', trim(errorString), ' "',trim(filenameIMDresults), '"',&
                                    'The program will be terminated!!!'
    call MPI_Abort(MPI_COMM_WORLD, mpiFileError, iError)
  end if

  call MPI_FILE_SET_VIEW(filehandle, myFileOffset, MPI_PACKED, MPI_PACKED, 'native', mpiInfo, iError)

  call MPI_FILE_READ_ALL(filehandle, AtomsBuffer, atomBufferSize, MPI_PACKED, ioStatus, iError)

  call MPI_FILE_CLOSE(filehandle, iError)

  do iPart=1,PDM%ParticleVecLength
    call MPI_UNPACK(AtomsBuffer, atomBufferSize, atomsBufferPos, PartState(iPart,1:6),&
                    6_4, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, iError)

    call MPI_UNPACK(AtomsBuffer, atomBufferSize, atomsBufferPos, PartStateIntEn(iPart,1:2),&
                    int( observables-6_8, 4 ), MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, iError)
  end do

  deallocate(AtomsBuffer)

  PartState = PartState * 1e-10_8
  PartState(:,4:6) = PartState(:,4:6) * 10.18e15_8

  call MPI_Info_free(mpiInfo, iError)

  PDM%ParticleInside(:) = .False.

  do iPart=1,PDM%ParticleVecLength
    PDM%ParticleInside(iPart) = .True.
    CALL SingleParticleToExactElementNoMap(iPart,doHALO=.TRUE.,doRelocate=.TRUE.)
    if( .not. PDM%ParticleInside(iPart) )then
      WRITE (*,*) "Particle Lost: iPart=", iPart," position=",PartState(iPart,1),PartState(iPart,2),PartState(iPart,3)
    end if
  end do


  CALL MPI_REDUCE(MaxX , MaxX_glob , 1 , MPI_DOUBLE_PRECISION , MPI_MAX , 0 , MPI_COMM_WORLD , iError)
  CALL MPI_REDUCE(MinX , MinX_glob , 1 , MPI_DOUBLE_PRECISION , MPI_MIN , 0 , MPI_COMM_WORLD , iError)

  CALL MPI_REDUCE(MaxY , MaxY_glob , 1 , MPI_DOUBLE_PRECISION , MPI_MAX , 0 , MPI_COMM_WORLD , iError)
  CALL MPI_REDUCE(MinY , MinY_glob , 1 , MPI_DOUBLE_PRECISION , MPI_MIN , 0 , MPI_COMM_WORLD , iError)

  CALL MPI_REDUCE(MaxZ , MaxZ_glob , 1 , MPI_DOUBLE_PRECISION , MPI_MAX , 0 , MPI_COMM_WORLD , iError)
  CALL MPI_REDUCE(MinZ , MinZ_glob , 1 , MPI_DOUBLE_PRECISION , MPI_MIN , 0 , MPI_COMM_WORLD , iError)

  if( MPIroot )then
    write(*,*) "Global particle information"
    write(*,*) "MinX_glob,MaxX_glob: ", MinX_glob,MaxX_glob
    write(*,*) "MinY_glob,MaxY_glob: ", MinY_glob,MaxY_glob
    write(*,*) "MinZ_glob,MaxZ_glob: ", MinZ_glob,MaxZ_glob
  end if



  call IRecvNbofParticles()
  call SendNbOfParticles()
  call MPIParticleSend()

  ! UpdateNextFreePosition must be called after particles are sent (because halo particles crash this routine) and
  ! before particles are received (where the next free position is required)
  call UpdateNextFreePosition()

  call MPIParticleRecv()

  !call WriteStateToHDF5( trim(meshfile), t, tFuture )
  call WriteStateToHDF5( trim(meshfile), 0.0_8, 0.0_8 )

  SWRITE(UNIT_stdOut,'(A)')'Read IMD data done'
  if( killPIClas )then
    call MPI_FINALIZE( iError )
    stop
  end if

end subroutine read_IMD_results

end module mod_readIMD
#endif /*MPI*/
#endif /*PARTICLES*/
