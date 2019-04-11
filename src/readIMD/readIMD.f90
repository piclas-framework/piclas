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

  ! -----------------------------------------------------------------------------

  if ( myRank == 0 )then
    write(*,*)'FilenameIMD: ', trim(filenameIMDresults)
    write(*,*)'killPIClas:  ', killPIClas
  end if

  if( killPIClas )then
    call MPI_FINALIZE( iError )
    stop
  end if
end subroutine read_IMD_results

end module mod_readIMD
