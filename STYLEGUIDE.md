# Style Guide

* Warum überhaupt ein Styleguide?
    * Schafft ein einheitliches Erscheinungsbild und damit eine einheitliche Code-Struktur
    * Macht das Programm verständlicher, da wichtige Informationen sofort ersichtlich sind
    * Zwingt zu überlegter und aktiver Mitarbeit
* Generelle Regeln:
    * Codesprache: Englisch
    * Maximal 132 Zeichen pro Zeile (incl. Kommentare)
    * Einrückungen: Generell 2 Leerzeichen (keine Tabs)
    * Zeilenumbrüche in Kommentaren -> nächste Zeile entsprechend einrücken
    * Kommentare von Modulen bzw. Input-/Output-Variablen: Doxygen-Style
    * Kommentare von Präprozessor-Anweisungen im C-StyleA comparison of Rosenbrock and ESDIRK methods

## Header von Funktionen und Subroutinen

    !==================================================================================================================================
    !> \brief Fills the solution array U with a initial solution.
    !>
    !> Fills the solution array U with a initial solution provided by the ExactFunc subroutine through interpolation. Function is
    !> specified with the IniExactFunc paramter.
    !==================================================================================================================================
    SUBROUTINE FillIni(NLoc,xGP,U)
    !----------------------------------------------------------------------------------------------------------------------------------
    ! MODULES
    USE MOD_PreProc
    USE MOD_Equation_Vars ,ONLY: IniExactFunc
    USE MOD_Exactfunc     ,ONLY: ExactFunc
    USE MOD_Mesh_Vars     ,ONLY: nElems
    IMPLICIT NONE
    !----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT/OUTPUT VARIABLES
    INTEGER,INTENT(IN)              :: NLoc                                    !< Polynomial degree of solution 
    REAL,INTENT(IN)                 :: xGP(3,    0:NLoc,0:NLoc,0:NLoc,nElems)  !< Coordinates of Gauss-points
    REAL,INTENT(OUT)                :: U(PP_nVar,0:NLoc,0:NLoc,0:NLoc,nElems)  !< Solution array
    !----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    INTEGER                         :: i,j,k,iElem
    !==================================================================================================================================
    
    ! Evaluate the initial solution at the nodes and fill the solution vector U. 
    DO iElem=1,nElems
      DO k=0,NLoc; DO j=0,NLoc; DO i=0,NLoc
        CALL ExactFunc(IniExactFunc,0.,xGP(1:3,i,j,k,iElem),U(:,i,j,k,iElem))
      END DO; END DO; END DO
    END DO
    END SUBROUTINE FillIni

## Variablen

* Präprozessor-Variablen: ```PP_$var```
     ```
     PP_nVar
     ```
     
* Zählvariablen: Zählindex (klein) + Beschreibung (erster Buchstabe groß)
    ```
    DO iVar=1,PP_nVar
    ```
    
* Variablen beginnen generell mit einem Großbuchstaben, zusammengesetzte Wörter beginnen generell groß
    ```
    ALLOCATE(ActualElem)
    ```
    
* Bei einzelnen Buchstaben: Klein am Anfang zusammengesetzer Variablennamen, sonst groß. Beides möglich bei Einzelbuchstaben. Ausnahmen möglich, aber nicht erwünscht
    ```
    hTilde, TildeH, (Elem%U)
    ```



## Funktionen und Kontrollstrukturen
* FORTRAN intrinsics generell groß
    ```
    ALLOCATE(), DO, MAX(), SQRT(), INT(), etc.
    ```
* END-Befehl abtrennen
    ```
    END DO, END IF, END SUBROUTINE
    ```
* Besonders bei längeren Schleifen, If-Abfragen, etc.: Kommentare hinzufügen

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

## USE
* Wichtig, immer USE Statements mit

    ```USE, ONLY: ...```
    verwenden
* gilt für Variablen und für Funktion/Prozedur
* Aussnahme: Init und Finalize Routinen der entsprechenden Module
