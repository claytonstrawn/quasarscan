#!/bin/bash -l

#SBATCH -N 2         #Use 1 node
#SBATCH -t 03:00:00  #Set 3 hr time limit
#SBATCH -q regular   #Submit to the regular QOS
#SBATCH -C haswell   #Use Haswell nodes
#SBATCH -A mp363     #Charge job to mp363 project

export HDF5_USE_FILE_LOCKING=FALSE
quasarscan/./run_one_new_snapshot_nersc.sh 2 17 200