#!/bin/bash
#PBS -N ZK3
#PBS -M stock@iag.uni-stuttgart.de
#PBS -m abe
#PBS -l mppwidth=1024
#PBS -l mppnppn=32
#PBS -l walltime=12:00:00             

cd $PBS_O_WORKDIR

aprun -n 1024 ./flexi SparkPlug_BGG_fine.ini DSMCSpecies.ini 1>std-1.out 2>err-1.out
#aprun -n 1024 ./flexi SparkPlug_BGG_fine.ini DSMCSpecies.ini [Restart-File] 1>std-X.out 2>err-X.out
