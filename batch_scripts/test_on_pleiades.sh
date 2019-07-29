#!/bin/bash

python quasarscan/create_qso_endpoints.py VELA_v2_art_07 /nobackupp2/cstrawn/mydir/VELA07v2/10MpcBox_csf512_a0.200.d
python quasarscan/get_coldens.py quasarscan/output/VELA_v2_art_07coldensinfo/0_of_448-agoraions_z4.0.txt 96 n
python quasarscan/create_qso_endpoints.py VELA_v2_art_07 /nobackupp2/cstrawn/mydir/VELA07v2/10MpcBox_csf512_a0.250.d
python quasarscan/get_coldens.py quasarscan/output/VELA_v2_art_07coldensinfo/0_of_448-agoraions_z4.0.txt 96 n
python quasarscan/create_qso_endpoints.py VELA_v2_art_07 /nobackupp2/cstrawn/mydir/VELA07v2/10MpcBox_csf512_a0.330.d
python quasarscan/get_coldens.py quasarscan/output/VELA_v2_art_07coldensinfo/0_of_448-agoraions_z4.0.txt 96 n
python quasarscan/create_qso_endpoints.py VELA_v2_art_07 /nobackupp2/cstrawn/mydir/VELA07v2/10MpcBox_csf512_a0.400.d
python quasarscan/get_coldens.py quasarscan/output/VELA_v2_art_07coldensinfo/0_of_448-agoraions_z4.0.txt 96 n
python quasarscan/create_qso_endpoints.py VELA_v2_art_07 /nobackupp2/cstrawn/mydir/VELA07v2/10MpcBox_csf512_a0.500.d
python quasarscan/get_coldens.py quasarscan/output/VELA_v2_art_07coldensinfo/0_of_448-agoraions_z4.0.txt 96 n
