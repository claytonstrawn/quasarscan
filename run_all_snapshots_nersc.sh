#!/bin/bash

#from VELA01 to VELA07, then VELA16 to VELA35
#10MpcBox_csf512_a0.020.d 
#from a0.330 to a0.500
for version in 1 2; do
    for i in 01 02 03 04 05 06 07 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35; do
        for a in 200 250 330 400 500 670; do
           echo $version $i $a
           quasarscan/./run_one_snapshot_nersc.sh $version $i $a
        done
    done
done 
