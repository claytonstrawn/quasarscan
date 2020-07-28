import sys
# import yt
# import trident
# import numpy as np
# import matplotlib.pyplot as plt
# import PIL
# from PIL import Image
import os
import datetime
import pytz
# import shutil
# import parse_metadata
import os.path
from os import path
from ionization_state_ytanalysis import *

def run_functions():
    simulations = ["VELA_v2_art_01", "VELA_v2_art_02", "VELA_v2_art_03", "VELA_v2_art_04", "VELA_v2_art_05", "VELA_v2_art_16",
                   "VELA_v2_art_17", "VELA_v2_art_18", "VELA_v2_art_19", "VELA_v2_art_20", "VELA_v2_art_21", "VELA_v2_art_22",
                   "VELA_v2_art_23", "VELA_v2_art_24", "VELA_v2_art_25", "VELA_v2_art_26", "VELA_v2_art_27", "VELA_v2_art_28",
                   "VELA_v2_art_29", "VELA_v2_art_30", "VELA_v2_art_31", "VELA_v2_art_32", "VELA_v2_art_33", "VELA_v2_art_34",
                   "VELA_v2_art_35"]
    files = []
    file_exists = []
    already_done = []
    redshifts = [1.0, 2.0]
    names = ["/10MpcBox_csf512_a0.500.d", "/10MpcBox_csf512_a0.400.d"]

    for z in range(len(redshifts)):
        print("Starting redshift = " + str(redshifts[z]))
        for i in range(len(simulations)):
            # edit this if redshift is changed!!
            files.append("/global/cfs/cdirs/mp363/SIP_INTERNS_2020/" + simulations[i] + names[i])
            file_exists.append(path.exists(files[i]))
            # edit this if redshift is changed!!
            filename = "ion_state_ytanalysis/" + simulations[i] + "_" + str(redshifts[z]) + "/o_ion_fraction_plot.png"
            already_done.append(path.exists(filename))

        for i in range(len(simulations)):
            current_time = datetime.datetime.now(pytz.timezone('America/Los_Angeles'))
            print("Time for simulation " + simulations[i] + " before Rahul's function: " + str(current_time))
            if(file_exists[i] == True and already_done[i] == False):
                # edit this if red shift is changed!!
                rahuls_function(simulations[i],files[i],redshifts[z])
                for j in range(9):
                    # edit this if redshift is changed!!
                    new_file = names + '_Projection_x_O_p' + str(j) + '_number_density.png'
                    if(os.path.exists(new_file)):
                        os.remove(new_file)
                current_time = datetime.datetime.now(pytz.timezone('America/Los_Angeles'))
                print("Time for simulation " + simulations[i] + " after Rahul's function: " + str(current_time))
                # edit this if redshift is changed!!
                sallys_function(simulations[i],files[i],redshifts[z])
                current_time = datetime.datetime.now(pytz.timezone('America/Los_Angeles'))
                print("Time for simulation " + simulations[i] + " after Sally's function: " + str(current_time))
            elif(already_done[i] == True): 
                print(simulations[i] + " was already done")
            else:
                print(simulations[i] + " file did not exist")

if __name__ == '__main__':
    run_functions()