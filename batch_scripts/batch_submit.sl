#!/bin/bash -l

#SBATCH -N 1         #Use 1 node
#SBATCH -t 00:30:00  #Set 6 hr time limit
#SBATCH -q debug   #Submit to the regular QOS
#SBATCH -C haswell   #Use Haswell nodes
#SBATCH -A mp363     #Charge job to mp363 project
#SBATCH -L SCRATCH   #Job requires $SCRATCH file system

simname=$1
filename=$2
z=$3
n=$4
procs=$5

module load python
conda activate myenv3v2
export HDF5_USE_FILE_LOCKING=FALSE
quasarscan/batch_scripts/./run_one_new_snapshot_nersc.sh $simname $filename $z $n $procs
