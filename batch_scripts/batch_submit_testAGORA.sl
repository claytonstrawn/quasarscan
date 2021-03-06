#!/bin/bash -l

#SBATCH -N 1         #Use 1 node
#SBATCH -t 00:30:00  #Set 3 hr time limit
#SBATCH -q debug     #Submit to the regular QOS
#SBATCH -C haswell   #Use Haswell nodes
#SBATCH -A agora     #Charge job to mp363 project
#SBATCH -L SCRATCH   #Job requires $SCRATCH file system

quasarscan/batch_scripts/./run_some_agora_files.sh