#!/bin/bash 
PICLAS_PATH="./piclas"
PICLAS_FLAGS="parameter.ini DSMCSpecies_N2.ini"

STD_INDEX=$(ls heap-allocation-*.out | tail -n 1 | sed s/heap-allocation-/\\n/ | tail -n 1 | sed s/.out/\\n/ | head -n 1 )
valgrind --tool=massif $PICLAS_PATH $PICLAS_FLAGS 2>&1 | tee massif-calc-`printf "%03d" $((${STD_INDEX##*0}+1))`.out

OUT_INDEX=$(ls massif.out.* | tail -n 1 | sed s/massif./\\n/ | grep out. | sed s/out./\\n/ | grep [[:digit:]] | head -n 1)
ms_print massif.out.`printf "%d" $(($OUT_INDEX))` | tee heap-allocation-`printf "%03d" $((${STD_INDEX##*0}+1))`.out
rm massif.out.`printf "%d" $(($OUT_INDEX))`
