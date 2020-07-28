import sys
# import yt
# import trident
# import numpy as np
# import matplotlib.pyplot as plt
# import PIL
# from PIL import Image
import os
# import shutil
# import parse_metadata
import os.path
from os import path
from ionization_state_ytanalysis import *

def run_functions():
#     simulations = ["VELA_v2_art_01", "VELA_v2_art_02", "VELA_v2_art_03", "VELA_v2_art_04", "VELA_v2_art_05", "VELA_v2_art_16",
#                    "VELA_v2_art_17", "VELA_v2_art_18", "VELA_v2_art_19", "VELA_v2_art_20", "VELA_v2_art_21", "VELA_v2_art_22",
#                    "VELA_v2_art_23", "VELA_v2_art_24", "VELA_v2_art_25", "VELA_v2_art_26", "VELA_v2_art_27", "VELA_v2_art_28",
#                    "VELA_v2_art_29", "VELA_v2_art_30", "VELA_v2_art_31", "VELA_v2_art_32", "VELA_v2_art_33", "VELA_v2_art_34",
#                    "VELA_v2_art_35"]
    simulations = ["VELA_v2_art_01", "VELA_v2_art_02", "VELA_v2_art_03", "VELA_v2_art_04"]
    files = []
    file_exists = []
    already_done = []

    for i in range(len(simulations)):
        # edit this if redshift is changed!!
        files.append("/global/cfs/cdirs/mp363/SIP_INTERNS_2020/" + simulations[i] + "/10MpcBox_csf512_a0.500.d")
        file_exists.append(path.exists(files[i]))
        # edit this if redshift is changed!!
        filename = "ion_state_ytanalysis/" + simulations[i] + "_" + str(1.0) + "/o_ion_fraction_plot.png"
        already_done.append(path.exists(filename))

    for i in range(len(simulations)):
        if(file_exists[i] == True and already_done[i] == False):
            # edit this if red shift is changed!!
            rahuls_function(simulations[i],files[i],1.0)
            for j in range(9):
                # edit this if redshift is changed!!
                new_file = '10MpcBox_csf512_a0.500.d_Projection_x_O_p' + str(j) + '_number_density.png'
                if(os.path.exists(new_file)):
                    os.remove(new_file)
            sallys_function(simulations[i],files[i],1.0)
        elif(already_done[i] == True): 
            print(simulations[i] + " was already done")
        else:
            print(simulations[i] + " file did not exist")

if __name__ == '__main__':
    run_functions()