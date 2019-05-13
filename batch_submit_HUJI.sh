#!/bin/bash -l

#SBATCH -n 1         #Use 1 cpu
#SBATCH -t 02:00:00  #Set 6 hr time limit
#SBATCH --mem 400m   #400 MB of memory per node
#SBATCH -c 1        #perform 10 tasks at once


export HDF5_USE_FILE_LOCKING=FALSE
source myenv/bin/activate.csh
quasarscan/./test_on_huji.sh
