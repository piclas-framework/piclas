#!/bin/bash 
PICLAS_PATH="./piclas"
PICLAS_FLAGS="parameter.ini DSMCSpecies_N2.ini"
FUNCTION="__mod_timedisc_MOD_timedisc"
PROFTODOT_PATH="./gprof2dot.py"

STD_INDEX=$(ls flowchart-init-*.png | tail -n 1 | sed s/flowchart-init-/\\n/ | tail -n 1 | sed s/.png/\\n/ | head -n 1 )
valgrind --tool=callgrind --dump-before=$FUNCTION $PICLAS_PATH $PICLAS_FLAGS 2>&1 | tee flowchart-calc-`printf "%03d" $((${STD_INDEX##*0}+1))`.out

OUT_INDEX=$(ls callgrind.out* | tail -n 2 | sed s/callgrind./\\n/ | grep out. | sed s/out./\\n/ | grep [[:digit:]] | head -n 1)
$PROFTODOT_PATH --format=callgrind -n1 callgrind.out.`printf "%d" $(($OUT_INDEX))`.1 | dot -Tpng -o flowchart-init-`printf "%03d" $((${STD_INDEX##*0}+1))`.png 
$PROFTODOT_PATH --format=callgrind -n1 callgrind.out.`printf "%d" $(($OUT_INDEX))` | dot -Tpng -o flowchart-timedisc-`printf "%03d" $((${STD_INDEX##*0}+1))`.png
rm callgrind.out.`printf "%d" $(($OUT_INDEX))`*
