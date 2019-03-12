#!/bin/bash -l

#SBATCH -c 1         #Use 1 cpu
#SBATCH -t 02:00:00  #Set 6 hr time limit
#SBATCH --mem 400m   #400 MB of memory per node
#SBATCH -n 10        #perform 10 tasks at once


export HDF5_USE_FILE_LOCKING=FALSE
quasarscan/./test_on_huji.sh