#!/bin/bash -l

#SBATCH -N 1         #Use 1 node
#SBATCH -t 00:20:00  #Set 30 minute time limit
#SBATCH -q debug   #Submit to the regular QOS
#SBATCH -C haswell   #Use Haswell nodes
#SBATCH -A mp363     #Charge job to mp363 project

export HDF5_USE_FILE_LOCKING=FALSE
quasarscan/./run_one_cont_snapshot_nersc.sh 2 17 4.0 srun