#!/bin/bash
module load python
source activate myenv
export HDF5_USE_FILE_LOCKING=FALSE
while :
do
    python quasarscan/decide_next.py
    chmod u+x quasarscan/nextfile.sh
    quasarscan/./nextfile.sh
done