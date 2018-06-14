MODULE MOD_PICInterpolation_Vars
!===================================================================================================================================
!> Variables for particle interpolation in PIC:
!> interpolation types, external fields (const. or variable), analytic interpolation
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE        :: FieldAtParticle(:,:)          !< (PIC%maxParticleNumber,6) 2nd index: Ex,Ey,Ez,Bx,By,Bz
CHARACTER(LEN=256)      :: InterpolationType             !< Type of Interpolation-Method
LOGICAL                 :: InterpolationElemLoop         !< Interpolate with outer iElem-loop (not for many Elems per proc!)
REAL                    :: externalField(6)              !< ext field is added to the maxwell-solver-field
LOGICAL                 :: DoInterpolation               !< Flag for interpolation
LOGICAL                 :: useBGField                    !< Flag for BGField via h5-File
INTEGER                 :: NBG                           !< Polynomial degree of BG-Field
INTEGER                 :: BGType                        !< Type of BG-Field (Electric,Magnetic,Both)
INTEGER                 :: BGDataSize                    !< Type of BG-Field (Electric,Magnetic,Both)
REAL, ALLOCATABLE       :: BGField(:,:,:,:,:)            !< BGField data
                                                         !< (1:x,0:NBG,0:NBG,0:NBG,1:PP_nElems)
REAL,ALLOCATABLE        :: BGField_xGP(:)                !< Gauss point coordinates
REAL,ALLOCATABLE        :: BGField_wGP(:)                !< GP integration weights
REAL,ALLOCATABLE        :: BGField_wBary(:)              !< barycentric weights


CHARACTER(LEN=256)      :: FileNameVariableExternalField !< filename containing the externanl field csv table
LOGICAL                 :: useVariableExternalField      !< use given external field. only for Bz variation in z
REAL,ALLOCATABLE        :: VariableExternalField(:,:)    !< z - Pos , Bz
REAL                    :: DeltaExternalField            !< equidistant z-spacing for the VariableExternalField (fast computation)
INTEGER                 :: nIntPoints                    !< number of all interpolation points external field

#ifdef CODE_ANALYZE
LOGICAL                 :: DoInterpolationAnalytic       !< use analytic/algebraic functions for the field at the
!                                                                      !< particle position
#endif /*CODE_ANALYZE*/
!===================================================================================================================================
END MODULE MOD_PICInterpolation_Vars
