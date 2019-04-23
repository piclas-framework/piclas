\hypertarget{guidelines}{}

# Development guidelines \label{chap:guidelines}

Code development is performed on the GitLab platform, with the protected `master` and `master.dev` branches. The actual development is performed on feature branches, which can be merged to `master.dev` following a merge request. After a successful pass of the nightly and weekly regression test, the `master.dev` can be merged into the `master`. A merge of the `master.dev` to the `master` should be associated with a release tag, where the changes to previous version are summarized.

In the following the envisioned development process, code style guidelines, the release & deploy procedure as well as other developer relevant issues are discussed.

## Development process
Issues are created for bugs, improvements, features, regression testing and documentation. The issues should be named with a few keywords. Try to avoid describing the complete issue already in the title. The issue can be assigned to a certain milestone (if appropriate).

Milestones are created based on planned releases (e.g. Release 1.2.1) or as a grouping of multiple related issues (e.g. Documentation Version 1, Clean-up Emission Routines). A deadline can be given if applicable. The milestone should contain a short summary of the work performed (bullet-points) as its contents will be added to the description of the releases. Generally, merge requests should be associated with a milestone containing a release tag, while issues should be associated with the grouping milestones.

As soon as a developer wants to start working on an issue, she/he shall assign himself to the issue and a branch and merge request denoted as work in progress (`WIP: ...`) should be created to allow others to contribute and track the progress. For this purpose, it should be created directly from the web interface within the issue (`Create merge request`). This creates a branch, naming it automatically, leading with the issue number (e.g. `60-fix-boundary-condition`) and associates the branch and merge request to the issue (visible in the web interface below the description). To start working on the issue, the branch can be checked out as usually.

Ideally, issues should be created for every code development for documentation purposes. Branches without an issue should be avoided to reduce the number of orphaned/stale branches. However, if branches are created outside of the issue context, they should be named with a prefix indicating the purpose of the branch, according to the existing labels in GitLab. Examples are given below:

    feature.chemistry.polyatomic
    improvement.tracking.curved
    bug.compiler.warnings
    reggie.chemistry.reservoir
    documentation.pic.maxwell

Progress tracking, documentation and collaboration on the online platform can be enabled through creating a merge request with the WIP prefix for this branch instead of an issue. An issues created afterwards cannot be associated with an already created branch, without renaming the branch to include the issue number at the beginning. However, this should be avoided.

## Release and deploy

After the successful completion of all regression checks (check-in, nightly, weekly), the master.dev branch can be merged into the master. This merge request should be associated with a milestone (e.g. Release 1.2.1)

### Create a Release Tag

**WORK IN PROGRESS**

### Collaborative Numerics Group

The master branch of development group can be merged after the successful regression check with the master of the collaborative group. For this purpose, the collaborative repository can be added as a remote

    git remote add remote_name git@gitlab.com:collaborative-numerics-group/piclas/piclas.git

Now you can checkout the most recent version of the master branch of the collaborative-numerics-group and create a local branch with that version (a simple checkout will create a detached HEAD state)

    git fetch
    git checkout -b branch_name remote_name/master

The master branch of the development repository can now be merged into the newly created branch. Make sure to have the most recent version of the master branch (of the development repository) as well.

    git merge origin/master

Finally, the changes can be pushed from the *branch_name* to the master of collaborative-numerics-group

    git push remote_name master

If a tag has also been created, it should be pushed separately.

    git push remote_name tag_name

### GitHub

Upon completion of a milestone leading to tagged version, the tag should be deployed to GitHub.

**WORK IN PROGRESS**

## Style Guide

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

### Header of Functions and Subroutines

Always use `USE` statements with `ONLY`

    USE MODULE, ONLY: ...
    
this accounts for variables and function/subroutines. An exception are the initilization and finalization routines.

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
    INTEGER,INTENT(IN)              :: NLoc                                    !< Polynomial degree of solution 
    REAL,INTENT(IN)                 :: xGP(3,    0:NLoc,0:NLoc,0:NLoc,nElems)  !< Coordinates of Gauss-points
    REAL,INTENT(OUT)                :: U(PP_nVar,0:NLoc,0:NLoc,0:NLoc,nElems)  !< Solution array
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

### Variables

* Preprocessor variables: `PP_$var`
    ```
    PP_nVar
    ```

* Counters: the counting variable (lower case) + description (the first character is capital case)
    ```
    DO iVar=1,PP_nVar
    ```

* Variables generally begin with a capital letter (composite words also)
    ```
    ALLOCATE(ActualElem)
    ```

* When using single characters: small at the beginning when using composite words otherwise in
  capital letters. Both is possible when purely single characters are used. Exceptions are allowed in
  special cases, but they are not recommened.
    ```
    hTilde, TildeH, (Elem%U)
    ```

### Functions and Control Structures
* FORTRAN intrinsics generally in capital letters
    ```
    ALLOCATE(), DO, MAX(), SQRT(), INT(), etc.
    ```
* END-X is to be separated by a space
    ```
    END DO, END IF, END SUBROUTINE
    ```
* For loops and `IF` statements etc. comments are to be inserted at the end (and inbetween, e.g. when
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
## Compiler flags
  * RELEASE: optimized with -O3 for execution runs
  * DEBUG: debugger options
  * SANI: GNU sanitizer for further debugging
  
    | Compiler-Flag           | Options,List  | What does it do?                                                                                                                                                                                                                                                                                                                                                                                                   |
    |-------------------------|---------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
    |     --ffpe-trap=list    |  *invalid*    | invalid floating point operation, such as SQRT(-1.0)                                                                                                                                                                                                                                                                                                                                                               |
    |                         |     *zero*    | division by zero                                                                                                                                                                                                                                                                                                                                                                                                   |
    |                         |   *overflow*  | overflow in a floating point operation                                                                                                                                                                                                                                                                                                                                                                             |
    |                         |  *underflow*  | underflow in a floating point. **DO NOT USE**. Because a small value can occure, such as exp(-766.2). operation                                                                                                                                                                                                                                                                                                      |
    |                         |  *precision*  | loss of precision during operation                                                                                                                                                                                                                                                                                                                                                                                 |
    |                         |               | Some of the routines in the Fortran runtime library, like **CPU_TIME**, are likely to trigger floating point exceptions when ffpe-trap=precision is used. For this reason, the use of ffpe-trap=precision is not recommended.                                                                                                                                                                                        |
    |                         |   *denormal*  | operation produced a denormal value                                                                                                                                                                                                                                                                                                                                                                                |
    | -fbacktrace             |               | runtime error should lead to a backtrace of the error                                                                                                                                                                                                                                                                                                                                                              |
    | -fcheck=keyword         |     *all*     | enable all run-time check                                                                                                                                                                                                                                                                                                                                                                                          |
    |                         | *array-temps* | Warns at run time when for passing an actual argument a temporary array had to be generated. The information generated by this warning is sometimes useful in optimization, in order to avoid such temporaries.                                                                                                                                                                                                    |
    |                         |    *bounds*   | Enable generation of run-time checks for array subscripts and against the declared minimum and maximum values. It also checks array indices for assumed and deferred shape arrays against the actual allocated bounds and ensures that all string lengths are equal for character array constructors without an explicit typespec.                                                                                 |
    |                         |      *do*     | Enable generation of run-time checks for invalid modification of loop iteration variables                                                                                                                                                                                                                                                                                                                          |
    |                         |     *mem*     | Enable generation of run-time checks for memory allocation. Note: This option does not affect explicit allocations using theALLOCATE statement, which will be always checked.                                                                                                                                                                                                                                      |
    |                         |   *pointer*   | Enable generation of run-time checks for pointers and allocatables.                                                                                                                                                                                                                                                                                                                                                |
    |                         | *recursion*   | Enable generation of run-time checks for recursively called subroutines and functions which are not marked as recursive. See also -frecursive. Note: This check does not work for OpenMP programs and is disabled if used together with -frecursiveand -fopenmp.                                                                                                                                                   |
    | -fdump-core             |               | Request that a core-dump file is written to disk when a runtime error is encountered on systems that support core dumps. This option is only effective for the compilation of the Fortran main program                                                                                                                                                                                                             |
    | -fstack-arrays          |               | Adding this option will make the Fortran compiler put all local arrays, even those of unknown size onto stack memory. If your program uses very large local arrays it is possible that you will have to extend your runtime limits for stack memory on some operating systems. This flag is enabled by default at optimization level -Ofast.                                                                       |
    | -frepack-arrays         |               | In some circumstances GNU Fortran may pass assumed shape array sections via a descriptor describing a noncontiguous area of memory. This option adds code to the function prologue to repack the data into a contiguous block at runtime.This should result in faster accesses to the array. However it can introduce significant overhead to the function call, especially when the passed data is noncontiguous. |
    | -finline-matmul-limit=n |               |                                                                                                                                                                                                                                                                                                                                                                                                                    |
    | -finit-local-zero       |               | The -finit-local-zero option instructs the compiler to initialize local INTEGER, REAL, and COMPLEX variables to zero, LOGICALvariables to false, and CHARACTER variables to a string of null bytes                                                                                                                                                                                                                 |
