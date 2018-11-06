#!/bin/bash 
STD_INDEX=$(ls std-*.out | tail -n 1 | sed s/std-/\\n/ | tail -n 1 | sed s/.out/\\n/ | head -n 1 )

mpirun -np 32 piclas parameter_pressurebc.ini DSMCSpecies.ini 1>std-`printf "%03d" $((${STD_INDEX##*0}+1))`.out 2>err-`printf "%03d" $((${STD_INDEX##*0}+1))`.out
wait
vim std-`printf "%03d" $((${STD_INDEX##*0}+1))`.out
