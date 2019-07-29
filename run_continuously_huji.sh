#!/bin/bash

source myenv3/bin/activate.csh
export HDF5_USE_FILE_LOCKING=FALSE
python quasarscan/decide_next_nihao.py
until ! [ -s quasarscan/nextfile.sh ]
do
    chmod u+x quasarscan/nextfile.sh
    quasarscan/./nextfile.sh
    python quasarscan/decide_next_nihao.py
done