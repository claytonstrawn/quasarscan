#!/bin/bash -l

#SBATCH -N 1         #Use 1 node
#SBATCH -t 03:00:00  #Set 3 hr time limit
#SBATCH -q regular   #Submit to the regular QOS
#SBATCH -C haswell   #Use Haswell nodes
#SBATCH -A mp363     #Charge job to mp363 project
#SBATCH -L SCRATCH   #Job requires $SCRATCH file system

quasarscan/./run_continuously.sh