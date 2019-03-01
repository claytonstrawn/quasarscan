#!/bin/bash
module load python
source activate myenv
export HDF5_USE_FILE_LOCKING=FALSE
#simname, filename
#z=20
quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_art_01 /global/cscratch1/sd/cstrawn/AGORAfiles/art/10MpcBox_csf512_a0.040.d 20
quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_ramses_01 /global/cscratch1/sd/cstrawn/AGORAfiles/ramses/output_00007/info_00007.txt 20
quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_gizmo_01 /global/cscratch1/sd/cstrawn/AGORAfiles/gizmo/snapshot.006.hdf5 20
quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_gadget_01 /global/cscratch1/sd/cstrawn/AGORAfiles/gadget/z20/snapshot_000.0.hdf5 20
quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_gear_01 /global/cscratch1/sd/cstrawn/AGORAfiles/gear/snapshot_0047.hdf5 20
quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_enzo_01 /global/cscratch1/sd/cstrawn/AGORAfiles/enzo/DD0008/DD0008 20

##z=15
#quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_art_01 /global/cscratch1/sd/cstrawn/AGORAfiles/art/10MpcBox_csf512_a0.062.d 15
#quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_ramses_01 /global/cscratch1/sd/cstrawn/AGORAfiles/ramses/output_00009/info_00009.txt 15
#quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_gizmo_01 /global/cscratch1/sd/cstrawn/AGORAfiles/gizmo/snapshot.009.hdf5 15
#quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_gadget_01 /global/cscratch1/sd/cstrawn/AGORAfiles/gadget/z15/snapshot_006.0.hdf5 15
#quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_gear_01 /global/cscratch1/sd/cstrawn/AGORAfiles/gear/snapshot_0056.hdf5 15
#quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_enzo_01 /global/cscratch1/sd/cstrawn/AGORAfiles/enzo/RD0009/RD0009 15

##z=10
##quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_art_01 /global/cscratch1/sd/cstrawn/AGORAfiles/art/10MpcBox_csf512_a0.062.d 10
#quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_ramses_01 /global/cscratch1/sd/cstrawn/AGORAfiles/ramses/output_00011/info_00011.txt 10
#quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_gizmo_01 /global/cscratch1/sd/cstrawn/AGORAfiles/gizmo/snapshot.017.hdf5 10
#quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_gadget_01 /global/cscratch1/sd/cstrawn/AGORAfiles/gadget/z10/snapshot_021.0.hdf5 10
#quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_gear_01 /global/cscratch1/sd/cstrawn/AGORAfiles/gear/snapshot_0067.hdf5 10
#quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_enzo_01 /global/cscratch1/sd/cstrawn/AGORAfiles/enzo/RD0017/RD0017 10

##z=8
##quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_art_01 /global/cscratch1/sd/cstrawn/AGORAfiles/art/10MpcBox_csf512_a0.062.d 8
#quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_ramses_01 /global/cscratch1/sd/cstrawn/AGORAfiles/ramses/output_00014/info_00014.txt 8
#quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_gizmo_01 /global/cscratch1/sd/cstrawn/AGORAfiles/gizmo/snapshot.022.hdf5 8
#quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_gadget_01 /global/cscratch1/sd/cstrawn/AGORAfiles/gadget/z8/snapshot_034.0.hdf5 8
#quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_gear_01 /global/cscratch1/sd/cstrawn/AGORAfiles/gear/snapshot_0074.hdf5 8
##quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_enzo_01 /global/cscratch1/sd/cstrawn/AGORAfiles/enzo/RD0017/RD0017 08

##z=6
##quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_art_01 /global/cscratch1/sd/cstrawn/AGORAfiles/art/10MpcBox_csf512_a0.062.d 6
#quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_ramses_01 /global/cscratch1/sd/cstrawn/AGORAfiles/ramses/output_00019/info_00019.txt 6
#quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_gizmo_01 /global/cscratch1/sd/cstrawn/AGORAfiles/gizmo/snapshot.029.hdf5 6
#quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_gadget_01 /global/cscratch1/sd/cstrawn/AGORAfiles/gadget/z6/snapshot_055.0.hdf5 6
#quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_gear_01 /global/cscratch1/sd/cstrawn/AGORAfiles/gear/snapshot_0083.hdf5 6
##quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_enzo_01 /global/cscratch1/sd/cstrawn/AGORAfiles/enzo/RD0017/RD0017 6

##z=5
##quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_art_01 /global/cscratch1/sd/cstrawn/AGORAfiles/art/10MpcBox_csf512_a0.062.d 5
##quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_ramses_01 /global/cscratch1/sd/cstrawn/AGORAfiles/ramses/output_00011/info_00011.txt 5
#quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_gizmo_01 /global/cscratch1/sd/cstrawn/AGORAfiles/gizmo/snapshot.035.hdf5 5
##quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_gadget_01 /global/cscratch1/sd/cstrawn/AGORAfiles/gadget/z10/snapshot_021.0.hdf5 5
#quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_gear_01 /global/cscratch1/sd/cstrawn/AGORAfiles/gear/snapshot_0088.hdf5 5
##quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_enzo_01 /global/cscratch1/sd/cstrawn/AGORAfiles/enzo/RD0017/RD0017 5

##z=4
##quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_art_01 /global/cscratch1/sd/cstrawn/AGORAfiles/art/10MpcBox_csf512_a0.062.d 4
##quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_ramses_01 /global/cscratch1/sd/cstrawn/AGORAfiles/ramses/output_00011/info_00011.txt 4
##quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_gizmo_01 /global/cscratch1/sd/cstrawn/AGORAfiles/gizmo/snapshot.017.hdf5 4
##quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_gadget_01 /global/cscratch1/sd/cstrawn/AGORAfiles/gadget/z10/snapshot_021.0.hdf5 4
#quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_gear_01 /global/cscratch1/sd/cstrawn/AGORAfiles/gear/snapshot_0095.hdf5 4
##quasarscan/./run_one_new_snapshot_nersc.sh AGORA_v1_enzo_01 /global/cscratch1/sd/cstrawn/AGORAfiles/enzo/RD0017/RD0017 4
