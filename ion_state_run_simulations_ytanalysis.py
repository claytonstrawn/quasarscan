import sys
import yt
import trident
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import os
import shutil
import parse_metadata
import os.path
from os import path
from ionization_state_ytanalysis import *

def run_functions:
#     simulations = ["VELA_v2_art_01", "VELA_v2_art_02", "VELA_v2_art_03", "VELA_v2_art_04", "VELA_v2_art_05", "VELA_v2_art_16",
#                    "VELA_v2_art_17", "VELA_v2_art_18", "VELA_v2_art_19", "VELA_v2_art_20", "VELA_v2_art_21", "VELA_v2_art_22",
#                    "VELA_v2_art_23", "VELA_v2_art_24", "VELA_v2_art_25", "VELA_v2_art_26", "VELA_v2_art_27", "VELA_v2_art_28",
#                    "VELA_v2_art_29", "VELA_v2_art_30", "VELA_v2_art_31", "VELA_v2_art_32", "VELA_v2_art_33", "VELA_v2_art_34",
#                    "VELA_v2_art_35"]
    simulations = ["VELA_v2_art_01", "VELA_v2_art_02"]
    files = []
    file_exists = []
    already_done = []

    for i in range(len(simulations)):
        files.append("/global/cfs/cdirs/mp363/SIP_INTERNS_2020/" + simulations[i] + "/10MpcBox_csf512_a0.500.d")
        file_exists.append(path.exists(files[i]))
        filename = "/global/cfs/cdirs/mp363/SIP_INTERNS_2020/" + simulations[i]
        filename = filename + "/" + simulations[i] + "_" + str(redshift) + "/o_ion_fraction_plot.png"
        already_done.append(path.exists(filename))

    for i in range(len(simulations)):
        if(file_exists[i] == True && already_done[i] == False):
            rahuls_function(simulations[i],files[i],1.0)
            sallys_function(simulations[i],files[i],1.0)

if __name__ == '__main__':
    run_functions()