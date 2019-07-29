#!/bin/bash -l

#SBATCH -n 1         #Use 1 cpu
#SBATCH -t 00:30:00  #Set 6 hr time limit
#SBATCH --mem 400m   #400 MB of memory per node
#SBATCH -c 1        #perform 10 tasks at once

source myenv3/bin/activate.csh
quasarscan/batch_scripts/./run_continuously.sh