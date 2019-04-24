#include "piclas.h"
module mod_readIMD

implicit none
private

public :: read_IMD_results, DefineParametersReadIMDdata, initReadIMDdata
contains

! ==============================================================================

subroutine DefineParametersReadIMDdata()

  use MOD_ReadInTools ,only: prms
  ! -----------------------------------------------------------------------------
  call prms%CreateLogicalOption('useIMDresults', 'Use Data from IMD', '.FALSE.')
  call prms%CreateStringOption('filenameIMDresults', 'Filename of the IMD result file')!, '')
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
  implicit none
  ! --------------------------------------------------------
  ! local variable
  integer                                   :: mpiInfo
  integer                                   :: filehandle
  integer                                   :: mpiFileError
  real,dimension(3)                         :: Box_X=.0_8, Box_Y=.0_8, Box_Z=.0_8
  integer(kind=MPI_OFFSET_KIND)             :: disp, myFileOffset
  integer(kind=2)                           :: disp_tmp
  integer(kind=2)                           :: observables
  integer(kind=8)                           :: nGlobalAtoms
  integer                                   :: ioStatus(MPI_STATUS_SIZE)
  integer                                   :: nAtoms, iAtom, iProc
  integer(kind=8),dimension(0:nProcessors)  :: FileOffsets
  integer(kind=MPI_OFFSET_KIND)             :: myOffset
  character(len=1),dimension(:),allocatable :: AtomsBuffer
  integer(kind=4)                           :: atomBufferSize
  integer(kind=8)                           :: ii
  real,dimension(:,:),allocatable           :: Atoms
  integer                                   :: atomsBufferPos = 0
  integer                                   :: errorLen
  character(len=254)                        :: errorString
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
    call MPI_FILE_READ_AT(filehandle,13_8,observables,2,MPI_BYTE,ioStatus,iError)

    call MPI_FILE_READ_AT(filehandle,15_8,Box_X,3,MPI_DOUBLE_PRECISION,ioStatus,iError)
    call MPI_FILE_READ_AT(filehandle,39_8,Box_Y,3,MPI_DOUBLE_PRECISION,ioStatus,iError)
    call MPI_FILE_READ_AT(filehandle,63_8,Box_Z,3,MPI_DOUBLE_PRECISION,ioStatus,iError)

    disp = int(disp_tmp,8)

    call MPI_FILE_CLOSE(filehandle, mpiFileError)

  end if

  call MPI_BCAST(Box_X, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, iError)
  call MPI_BCAST(Box_Y, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, iError)
  call MPI_BCAST(Box_Z, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, iError)
  call MPI_BCAST(nGlobalAtoms, 8, MPI_BYTE, 0, MPI_COMM_WORLD, iError)
  call MPI_BCAST(observables, 2, MPI_BYTE, 0, MPI_COMM_WORLD, iError)
  call MPI_BCAST(disp, 8, MPI_BYTE, 0, MPI_COMM_WORLD, iError)

  nAtoms = nGlobalAtoms/nProcessors
  iAtom = nGlobalAtoms - nAtoms * nProcessors
  FileOffsets(0) = 0
  DO iProc=0,nProcessors-1
    FileOffsets(iProc) = nAtoms * iProc + MIN( iProc,iAtom )
  END DO
  FileOffsets(nProcessors) = nGlobalAtoms
  nAtoms = FileOffsets(myRank+1) - FileOffsets(myRank)
  myOffset = FileOffsets(myRank)

  myFileOffset = disp + myOffset * observables * 8_8
  atomBufferSize = 8 * observables * nAtoms
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

  allocate(Atoms(observables,nAtoms))

  do ii=1,nAtoms
    call MPI_UNPACK(AtomsBuffer, atomBufferSize, atomsBufferPos, Atoms(:,ii),&
                    int(observables,4), MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, iError)
  end do

  deallocate(AtomsBuffer)
  call MPI_Info_free(mpiInfo, iError)

  if( killPIClas )then
    call MPI_FINALIZE( iError )
    stop
  end if
end subroutine read_IMD_results

end module mod_readIMD
