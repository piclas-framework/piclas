This directory contains scripts used to profile PICLas
Currently only working with no mpi

CONTAINS:
===================================================================================================================================
callgrind.sh:
===================================================================================================================================
  This script uses valgrind-function callgrind
  Callgrind is a profiling tool that records the call history among functions in a program's run as a call-graph 
  with the percentage weighted runtime necessery for each function
------------------------------------------------------------------------------------------------------------------------------------
  NEDDED APPLICATIONS:
    - Valgrind    [ http://valgrind.org/docs/manual/cl-manual.html ]
    - gprof2dot.py    [ https://github.com/jrfonseca/gprof2dot ]
------------------------------------------------------------------------------------------------------------------------------------
  Neccessary EDITS inside the script:
    PICLAS_PATH="./piclas"
    PICLAS_FLAGS="parameter.ini DSMCSpecies_N2.ini"
    PROFTODOT_PATH="./gprof2dot.py"
------------------------------------------------------------------------------------------------------------------------------------
  OUTPUT:
    flowchart-init-###.png
      - graph showing callhistory in PICLas BEFORE timedisc starts
    flowchart-timedisc-###.png
      - graph showing callhistory in PICLas AFTER timedisc starts
===================================================================================================================================

massif.sh:
===================================================================================================================================
  This script uses valgrind-function massif
  Massif is a heap profiler. It measures how much heap memory your program uses.
------------------------------------------------------------------------------------------------------------------------------------
  NEDDED APPLICATIONS:
    - Valgrind    [ http://valgrind.org/docs/manual/ms-manual.html ]
------------------------------------------------------------------------------------------------------------------------------------
  Neccessary EDITS inside the script:
    PICLAS_PATH="./piclas"
    PICLAS_FLAGS="parameter.ini DSMCSpecies_N2.ini"
    
    not many iterations needed
------------------------------------------------------------------------------------------------------------------------------------
  OUTPUT:
    heap-allocation-###.out
      - file showing allocation history in PICLas
===================================================================================================================================

