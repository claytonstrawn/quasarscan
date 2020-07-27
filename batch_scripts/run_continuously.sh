#!/bin/bash

n=$1
module load python
conda activate myenv3v2
chmod u+x quasarscan/batch_scripts/run_one_new_snapshot_nersc.sh

export HDF5_USE_FILE_LOCKING=FALSE
/global/homes/c/cstrawn/.conda/envs/myenv3v2/bin/python quasarscan/decide_next_onlyVELA.py $n
until ! [ -s quasarscan/nextfile.sh ]
do
    chmod u+x quasarscan/nextfile.sh
    quasarscan/./nextfile.sh
    /global/homes/c/cstrawn/.conda/envs/myenv3v2/bin/python quasarscan/decide_next_onlyVELA.py $n
done