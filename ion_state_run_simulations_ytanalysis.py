import sys
import os
import datetime
import pytz
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
    names = ["10MpcBox_csf512_a0.500.d", "10MpcBox_csf512_a0.400.d"]

    for z in range(len(redshifts)):
        print("Starting redshift = " + str(redshifts[z]))
        files = []
        file_exists = []
        already_done = []
        for i in range(len(simulations)):
            # edit this if redshift is changed!!
            files.append("/global/cfs/cdirs/mp363/SIP_INTERNS_2020/" + simulations[i] + '/' + names[z])
            file_exists.append(path.exists(files[i]))
            # edit this if redshift is changed!!
            filename = "quasarscan/ion_state_ytanalysis/" + simulations[i] + "_" + str(redshifts[z]) + "/o_ion_fraction_plot.png"
            already_done.append(path.exists(filename))

        for i in range(len(simulations)):
            current_time = datetime.datetime.now(pytz.timezone('America/Los_Angeles'))
            print("Time for simulation " + simulations[i] + " before Rahul's function: " + str(current_time))
            if(file_exists[i] == True and already_done[i] == False):
                ds, redshift, Rvir, Mvir, center = setup(simulations[i],files[i],redshifts[z])
                # edit this if red shift is changed!!
                rahuls_function(simulations[i],files[i], ds, redshift, Rvir, Mvir, center)
                for j in range(9):
                    # edit this if redshift is changed!!
                    new_file = names[z] + '_Projection_x_O_p' + str(j) + '_number_density.png'
                    if(os.path.exists(new_file)):
                        os.remove(new_file)
                    new_file = names[z] + '_Projection_x_O_p' + str(j) + '_PI_number_density.png'
                    if(os.path.exists(new_file)):
                        os.remove(new_file)
                    new_file = names[z] + '_Projection_x_O_p' + str(j) + '_CI_number_density.png'
                    if(os.path.exists(new_file)):
                        os.remove(new_file)
                current_time = datetime.datetime.now(pytz.timezone('America/Los_Angeles'))
                print("Time for simulation " + simulations[i] + " after Rahul's function: " + str(current_time))
                # edit this if redshift is changed!!
                sallys_function(simulations[i],files[i],ds, redshift, Rvir, Mvir, center)
                current_time = datetime.datetime.now(pytz.timezone('America/Los_Angeles'))
                print("Time for simulation " + simulations[i] + " after Sally's function: " + str(current_time))
            elif(already_done[i] == True): 
                print(simulations[i] + "for redshift z = " + str(redshifts[z]) + " was already done")
            else:
                print(simulations[i] + "for redshift z = " + str(redshifts[z]) + " file did not exist")

if __name__ == '__main__':
    run_functions()