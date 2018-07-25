#!/bin/bash -l

#SBATCH -N 1         #Use 1 node
#SBATCH -t 00:20:00  #Set 50 minute time limit
#SBATCH -q debug     #Submit to the regular QOS
#SBATCH -L SCRATCH   #Job requires $SCRATCH file system (not sure what this means)
#SBATCH -C haswell   #Use Haswell nodes

module load python
source activate myenv

export HDF5_USE_FILE_LOCKING=FALSE
srun -n 10 python quasarscan/./quasar_scan.py 'n' 'VELA_v2_17/10MpcBox_csf512_a0.200.d' 'VELA_v2_17' '-p' 
