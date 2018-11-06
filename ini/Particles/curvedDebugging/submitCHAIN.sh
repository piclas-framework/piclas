#!/bin/bash
#MSUB -N onePartDEBUG6-N4
#MSUB -M binder@irs.uni-stuttgart.de
#MSUB -m bea
#MSUB -l nodes=2:ppn=16
#MSUB -l walltime=0:04:00:00
#MSUB -q multinode


#--------------------------------------------------------------------------------------------------
# ForHLR I
#--------------------------------------------------------------------------------------------------

# switch to submit dir 
cd $MOAB_SUBMITDIR

## State and std-out
STATE='' #$(ls *_State_*.h5 | tail -n 1 | sed s/:/\\n/ | sed s/\ /\\n/ | tail -n 1)
STD_INDEX=$(ls std-*.out | tail -n 1 | sed s/std-/\\n/ | tail -n 1 | sed s/.out/\\n/ | head -n 1 ) 
#cp -p $STATE ${STATE}_ORIG 

## --------------------------------impi--------------------------------
module swap compiler/intel compiler/intel/16.0
module load mpi/impi/5.1

mpiexec.hydra -bootstrap slurm ./piclas pic.ini DSMCSpecies_CEXMEX.ini $STATE 1>std-`printf "%03d" $((${STD_INDEX##*0}+1))`.out 2>err-`printf "%03d" $((${STD_INDEX##*0}+1))`.out
