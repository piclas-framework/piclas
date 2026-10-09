# Style Guide

* Why do we need a style guide?
    * It creates a unified appearance and coding structure
    * It makes the code more understandable and therefore important information is understood more
        easily
    * It forces the developers to think more actively about their work
* General rules
    * Coding language: English
    * A maximum of 132 characters are allowed per line (incl. Comments)
    * Indentation: 2 spaces (no tabs!)
    * Line breaks in comments -> the following line must be indented appropriately
    * Comments of modules and input-/output variables: Doxygen style
    * Comments of preprocessor directives in C-Style

## File Header

Every new source file must start with the following license header. When editing existing files, feel free to change the copyright statement to the current line using `PICLas Developers, see CONTRIBUTORS.md` and only keeping the first year.

```Fortran
!==================================================================================================================================
! Copyright (c) 2010 PICLas Developers, see CONTRIBUTORS.md
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
```

## Header of Functions and Subroutines

Function calls must always supply the variable name of optional arguments. Always use `USE` statements with `ONLY`

    USE MODULE, ONLY: ...

this accounts for variables and function/subroutines. An exception are the initialization and finalization routines.

    !==============================================================
    !> \brief Fills the solution array U with a initial solution.
    !>
    !> Fills the solution array U with a initial solution provided by the ExactFunc subroutine through interpolation. Function is
    !> specified with the IniExactFunc paramter.
    !==============================================================
    SUBROUTINE FillIni(NLoc,xGP,U)
    !--------------------------------------------------------------
    ! MODULES
    USE MOD_PreProc
    USE MOD_Equation_Vars ,ONLY: IniExactFunc
    USE MOD_Exactfunc     ,ONLY: ExactFunc
    USE MOD_Mesh_Vars     ,ONLY: nElems
    IMPLICIT NONE
    !--------------------------------------------------------------
    ! INPUT/OUTPUT VARIABLES
    INTEGER,INTENT(IN)              :: NLoc                                      !< Polynomial degree of solution
    REAL,INTENT(IN)                 :: xGP(3,    0:NLoc,0:NLoc,0:NLoc,nElems)    !< Coordinates of Gauss-points
    REAL,INTENT(OUT)                :: U(1:PP_nVar,0:NLoc,0:NLoc,0:NLoc,nElems)  !< Solution array
    !--------------------------------------------------------------
    ! LOCAL VARIABLES
    INTEGER                         :: i,j,k,iElem
    !==============================================================

    ! Evaluate the initial solution at the nodes and fill the solution vector U.
    DO iElem=1,nElems
      DO k=0,NLoc; DO j=0,NLoc; DO i=0,NLoc
        CALL ExactFunc(IniExactFunc,0.,xGP(1:3,i,j,k,iElem),U(:,i,j,k,iElem))
      END DO; END DO; END DO
    END DO
    END SUBROUTINE FillIni

The separators `!====` and `!----` are exactly 132 characters long (here they have been shortened for visualization purposes).

## Variables

* Preprocessor variables: `PP_$var`
    ```
    PP_nVar
    ```

    Note that

    * `USE MOD_Preproc` cannot be used with `ONLY` because the pro-processor flags sometimes result in constants and not variables
    * `PP_N` and other pre-processor variables that may be constants cannot be assigned in `ASSOCIATE` constructs

* Counters: the counting variable (lower case) + description (the first character is capital case)
    ```
    DO iVar=1,PP_nVar
    ```

* Variables generally begin with a capital letter (composite words also)
    ```
    ActualElem
    ```

* Dimension allocations must be specified by giving both a lower and an upper boundary
    ```
    ALLOCATE(U(1:PP_nVar,0:NLoc,0:NLoc,0:NLoc,nElems))
    ```

* When using single characters: small at the beginning when using composite words otherwise in
  capital letters. Both is possible when purely single characters are used. Exceptions are allowed in
  special cases, but they are not recommended.
    ```
    hTilde, TildeH, (Elem%U)
    ```

## Functions and Control Structures
* User-defined functions and subroutines should carry meaning in their name. If their name is
    composed of multiple words, they are to be fused together without underscores (`_`) and the
    first letter of each words should be capital.
    ```
    GetParticleWeight(), isChargedParticle()
    ```
    An exception to this rule is the usage of underscores (`_`) for shared memory arrays, where
    `_Shared` indicates that the property is available to all processors on the same shared-memory
    domain (generally a node).
    Furthermore, the words `Get`, `get`, `Is`, `is`, `Do`, `do` indicate the intention of the
    function/subroutine in supplying or acquiring specific properties or values.
    User-defined functions and subroutines should not, in general, be named in all-capital letters.
* FORTRAN intrinsics generally in capital letters
    ```
    ALLOCATE(), DO, MAX(), SQRT(), INT(), etc.
    ```
* END-X is to be separated by a space
    ```
    END DO, END IF, END SUBROUTINE
    ```
* For loops and `IF` statements etc. comments are to be inserted at the end (and in-between, e.g. when
`ELSE IF` is used)

    ```
    DO iVar=1,PP_nVar
      IF (a.EQ.b) THEN
    ...
      ELSE ! a.NE.b
    ...
      END IF ! a.EQ.b
    ...
    END DO ! PP_nVar
    ```

## Workflow Description
Additionally to the header description, a short workflow table of contents at the beginning of the
subroutine or function  is required for longer subroutines in which multiple tasks are completed.
Example:

    ! -----------------------------------------------------------------------------
    ! MAIN STEPS        []=FV only
    ! -----------------------------------------------------------------------------
    ! 1.  Filter solution vector
    ! 2.  Convert volume solution to primitive
    ! 3.  Prolong to face (fill U_master/slave)
    ! 4.  ConsToPrim of face data (U_master/slave)
    ![5.] Second order reconstruction for FV
    ! 6.  Lifting
    ! 7.  Volume integral (DG only)
    ![8.] FV volume integral
    ! 9.  IF EDDYVISCOSITY: Prolong muSGS to face and send from slave to master
    ! 10. Fill flux (Riemann solver) + surface integral
    ! 11. Ut = -Ut
    ! 12. Sponge and source terms
    ! 13. Perform overintegration and apply Jacobian
    ! -----------------------------------------------------------------------------

Furthermore, the steps are required to be found at the appropriate position within the code. It is
not allowed to just incorporate the corresponding number of the step within the code.

    ! (0. Nullify arrays)
    ! NOTE: UT and U are nullified in DGInit, and Ut is set directly

    ! 1. Filter the solution vector if applicable, filter_pointer points to cut-off
    IF(FilterType.GT.0) CALL Filter_Pointer(U,FilterMat)

    ! 2. Convert Volume solution to primitive
    CALL ConsToPrim(PP_N,UPrim,U)

    ! X. Update mortar operators and neighbour connectivity for the sliding mesh
    CALL PrepareSM()

    ! 3. Prolong the solution to the face integration points for flux computation
    ! --------------------------------------------------------------
    ! General idea: The slave sends its surface data to the master
    ! where the flux is computed and sent back to the slaves.
    ! Steps:
    ! (these steps are done for all slave MPI sides and then for all remaining sides):
    ! 3.1)  Prolong solution to faces and store in U_master/slave.
    !       Use them to build mortar data (split into 2/4 smaller sides).
    ![3.2)] The information which element is a DG or FV subcells element is stored
    !       in FV_Elems per element.
    ![3.3)] The reconstruction of slopes over element interfaces requires,
    !       besides U_slave and FV_Elems_slave, some more information that
    !       has to be transmitted from the slave to the master MPI side.
    ! 3.4)  Finish all started MPI communications (after step 2. due to latency hiding)

    #if USE_MPI
    ! Step 3 for all slave MPI sides
    ! 3.1) Prolong solution to faces and store in U_master/slave.
    !      Use them to build mortar data (split into 2/4 smaller sides).
    CALL StartReceiveMPIData(U_slave,DataSizeSide,1,nSides,MPIRequest_U(:,SEND),SendID=2) ! Receive MINE / U_slave: slave -> master
    CALL StartReceiveSM_MPIData(PP_nVar,U_MorRot,MPIRequestSM_U,SendID=2) ! Receive MINE / U_slave: slave -> master
    CALL ProlongToFaceCons(PP_N,U,U_master,U_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
    CALL U_MortarCons(U_master,U_slave,doMPISides=.TRUE.)
    CALL U_MortarConsSM(U_master,U_slave,U_MorStat,U_MorRot,doMPISides=.TRUE.)
    CALL StartSendMPIData(   U_slave,DataSizeSide,1,nSides,MPIRequest_U(:,RECV),SendID=2) ! SEND YOUR / U_slave: slave -> master
    CALL StartSendSM_MPIData(   PP_nVar,U_MorRot,MPIRequestSM_U,SendID=2) ! SEND YOUR / U_slave: slave -> master
    #if FV_ENABLED
    ! 3.2) The information which element is a DG or FV subcells element is stored
    !      in FV_Elems per element.
    CALL FV_Elems_Mortar(FV_Elems_master,FV_Elems_slave,doMPISides=.TRUE.)
    CALL StartExchange_FV_Elems(FV_Elems_slave,1,nSides,MPIRequest_FV_Elems(:,SEND),MPIRequest_FV_Elems(:,RECV),SendID=2)
    #endif /* FV_ENABLED */

## Special Rules
### CALL Allocate_Shared()
Subroutine calls for *Allocate_Shared()* must be placed in a single line, i.e., no line break is allowed due to the compile flag
*PICLAS_DEBUG_MEMORY*, which changes the number of arguments of the call for this routine by adding a debugging parameter.

### USE MOD_Preproc
Using variables via *USE MOD_Preproc, ONLY: PP_N* is not possible due to the compile flag *PICLAS_POLYNOMIAL_DEGREE*, which switches
the polynomial degree between a constant values, e.g., 1, 2, 3 .. or the default value of *N*. Therefore, the *ONLY* statement is
not allowed here.

## Code Annotations
In the context of the tags (`TODO`, `FIXME`, `HACK`, `NOTE`, `OPTIMIZE`), code annotations are standardized
keywords used within comments to label specific tasks, issues, or insights directly in the source code.
They offer the two main functionalities:

1. Searchability: They act as "bookmarks." Instead of searching for "to do," you can search for `TODO:`
   across thousands of files to see exactly what’s left to finish.
1. Tool Integration: Modern IDEs (e.g. VS Code, IntelliJ, Xcode) scan your project for these keywords and
   compile them into a Task List or `TODO` Board, so you don't have to remember where you left a note.

|Keyword    | Description                                                                                                                                                   |
| --------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|`TODO`     | Keep track of the tasks that have to completed in the feature branch. For optimization and improvement, instead use the keyword `OPTIMIZE`.                   |
|`FIXME`    | Broken or buggy code that requires immediate attention. It marks places where a bug is suspected or confirmed but hasn't been fixed yet.                      |
|`NOTE`     | Provide general information, explanations, or observations about a specific section of code and explains "the why" behind a non-obvious implementation.       |
|`HACK`     | Unorthodox, non-ideal, or "quick and dirty" solutions. Also used when bypassing a bug in a third-party framework or handle an edge case in an inelegant way.  |
|`OPTIMIZE` | Mark areas of code that are fully functional but inefficient: Single-core performance optimization, improved multi-core scaling or reduced memory consumption.|
|`FEATURE`  | Track a specific missing functionality in the code. Marks a new feature, which is not implemented yet, but intended to be in the future.                      |