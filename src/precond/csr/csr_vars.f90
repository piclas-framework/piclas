MODULE MOD_CSR_Vars
!===================================================================================================================================
! Contains global variables used by the DG modules.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! sparse test for calculation
INTEGER,ALLOCATABLE,DIMENSION(:,:)    :: nNonZeros
!REAL,ALLOCATABLE,DIMENSION(:,:)       :: Ucsr,Utcsr
TYPE tSparse
 REAL,ALLOCATABLE,DIMENSION(:)        :: Entry
 INTEGER,ALLOCATABLE,DIMENSION(:)     :: IEntry,JEntry
END TYPE                                                                     
REAL,ALLOCATABLE                      :: DebugMat(:,:)
TYPE(tSparse), ALLOCATABLE            :: SparseMatrix(:,:)
REAL,ALLOCATABLE                      :: L_HatPlusMinus(:,:)
REAL,ALLOCATABLE                      :: GlobalAA(:)
REAL,ALLOCATABLE                      :: GlobalBAA(:,:,:)
INTEGER,ALLOCATABLE                   :: GlobalIA(:),GlobalJA(:), GlobalDiag(:)
INTEGER                               :: nonZerosGlobal
REAL                                  :: epsZero
INTEGER,ALLOCATABLE,DIMENSION(:)      :: nUNonZeros,nLNonZeros
INTEGER                               :: nMTriangle
REAL,ALLOCATABLE,DIMENSION(:,:)       :: DE
TYPE tIL                                                                     ! ILU for each element
 REAL,ALLOCATABLE,DIMENSION(:)        :: Entry
 INTEGER,ALLOCATABLE,DIMENSION(:)     :: IEntry,JEntry
END TYPE                                                                     
TYPE(tIL), ALLOCATABLE                :: IL(:)             
TYPE(tIL), ALLOCATABLE                :: IU(:)             
!TYPE tAA                                                                     ! ILU for each element
! REAL,ALLOCATABLE,DIMENSION(:)        :: AA
! INTEGER,ALLOCATABLE,DIMENSION(:)     :: IA,JA,Diag
!END TYPE                                                                     
!TYPE(taa),ALLOCATABLE,DIMENSION(:)    :: EILU
!===================================================================================================================================
END MODULE MOD_CSR_Vars
