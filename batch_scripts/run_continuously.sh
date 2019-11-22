#!/bin/bash

n=$1

source activate myenv3
export HDF5_USE_FILE_LOCKING=FALSE
python quasarscan/decide_next_onlyVELA.py
until ! [ -s quasarscan/nextfile.sh ]
do
    chmod u+x quasarscan/nextfile.sh
    quasarscan/./nextfile.sh
    python quasarscan/decide_next_onlyVELA.py $n
done