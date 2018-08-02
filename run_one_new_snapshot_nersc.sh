#!/bin/bash

#from VELA01 to VELA07, then VELA16 to VELA35
#10MpcBox_csf512_a0.020.d 
#from a0.330 to a0.500

SOURCEDIR1="/home/c/ceverino/VELA"
SOURCEDIR2="/project/projectdirs/mp363/VELA_ANALYSIS"

i=$2
a=$3
z=$4

if [[ $1 = "1" ]]
then
    v=""
    DIRNAME=$SOURCEDIR2
    mkdir VELA$v$i
    cd VELA$v$i
    scp $DIRNAME/VELA$i/*a0.$a* .
else
    v="_v"$1"_"
    DIRNAME=$SOURCEDIR1"_v"$1
    mkdir VELA$v$i
    cd VELA$v$i
    hsi "cd $DIRNAME/VELA$i; get *a0.$a*"
fi
cd ..



module load python
source activate myenv

quasarscan/quasar_scan.py n /global/homes/c/cstrawn/VELA$v$i/10MpcBox_csf512_a0.$a.d VELA$v$i
srun -n 32 -u python -u quasarscan/get_coldens.py -fn quasarscan/output/VELA$v${i}coldensinfo/0_of_400-_AlII_AlIII_ArI_ArII_ArVII_CI_CII_CIII_CIV_CaX_FeII_HI_MgII_MgX_NI_NII_NIII_NIV_NV_NaIX_NeV_NeVI_NeVII_NeVIII_OI_OII_OIII_OIV_OV_OVI_PIV_PV_SII_SIII_SIV_SV_SVI_SXIV_SiII_SiIII_SiIV_SiXII_z$z.txt -s 96

echo "remember to rm -rf VELA$v$i"
source deactivate
