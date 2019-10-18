#!/bin/bash -l

#SBATCH -N 4         #Use 1 node
#SBATCH -t 05:00:00  #Set 6 hr time limit
#SBATCH -q regular   #Submit to the regular QOS
#SBATCH -C haswell   #Use Haswell nodes
#SBATCH -A mp363     #Charge job to mp363 project
#SBATCH -L SCRATCH   #Job requires $SCRATCH file system
export HDF5_USE_FILE_LOCKING=FALSE

module load python
source activate myenv3
#python quasarscan/create_qso_endpoints.py VELA_v2_art_08 /global/cscratch1/sd/cstrawn/VELA_v2_08/10MpcBox_csf512_a0.500.d  
srun -n 128 -u python -u quasarscan/get_coldens.py quasarscan/output/VELA_v2_art_08coldensinfo/0_of_448-agoraions_z1.0.txt 128 p
