#!/bin/bash -l

#SBATCH -N 1         #Use 1 node
#SBATCH -t 03:30:00  #Set 6 hr time limit
#SBATCH -q regular   #Submit to the regular QOS
#SBATCH -C haswell   #Use Haswell nodes
#SBATCH -A mp363     #Charge job to mp363 project
#SBATCH -L SCRATCH   #Job requires $SCRATCH file system

v=$1
num=$2
a=$3
z=$4
n=$5

module load python
conda activate myenv3v2
export HDF5_USE_FILE_LOCKING=FALSE
quasarscan/batch_scripts/./run_one_new_snapshot_nersc.sh $v $num $a $z $n
