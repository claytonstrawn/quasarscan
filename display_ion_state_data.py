#***import statements***
import numpy as np
import matplotlib.pyplot as plt
import os
# from quasarscan import parse_metadata
import parse_metadata

#define analysis functions

def get_sorted_array(criteria):
    #skip this one for now
    pass
    
def read_ion_masses(list_of_galaxies):
    array = np.empty((len(list_of_galaxies),27)
    num = len(list_of_galaxies)
    print(num)
    for i in range(num): 
        filename = 'ion_state_ytanalysis/' + list_of_galaxies[i] + '/o_ion_fraction_data.txt'
        f = open(filename, "r")
        contents = f.readline()
        curindex = 9
        for j in range(8):
            array[i,j] = int(contents[curindex:contents.index(', ',curindex)])
            curindex = contents.index(',',curindex + 1) + 2
        array[i,8] = int(contents[curindex:contents.index(']',curindex)])
        f.readline()
        contents = f.readline()
        curindex = 10
        for j in range(8):
            array[i,j+9] = int(contents[curindex:contents.index(', ',curindex)])
            curindex = contents.index(',',curindex + 1) + 2
        array[i,17] = int(contents[curindex:contents.index(']',curindex)])
        f.readline()
        contents = f.readline()
        curindex = 13
        for j in range(8):
            array[i,j+18] = int(contents[curindex:contents.index(', ',curindex)])
            curindex = contents.index(',',curindex + 1) + 2
        array[i,26] = int(contents[curindex:contents.index(']',curindex)])
    print(array)
    return array

def read_ion_number_densities(list_of_galaxies):
    return array

def combined_mass_plots(list_of_galaxies,data_array,weighting = 'mass',logy=True, plot_data = 'all'):
    #possible weightings = 'mass' (add all the stacked masses) 
    #and 'snapshot' (average all the fractions for individual galaxies)
    #possible plot_datas = 'all' (plot PI, CI, and total)
    #and 'PI', 'CI', or 'total' if you just want one 
    #processing of the data array (list_of_galaxies should generate captions or something)
    plt.plot()

def covering_fraction_plots(list_of_galaxies,data_array,thresholds = [13,14,15],logy=False, plot_data = 'all'):
    #possible plot_datas = 'all' (plot PI, CI, and total)
    #and 'PI', 'CI', or 'total' if you just want one 
    #and 'snapshot' (average all the fractions for individual galaxies)
    #processing of the data array (list_of_galaxies should generate captions or something)
    plt.plot()

if __name__ == '__main__':
	# this one doesn't really make sense as a script
    # I'd probably just import these functions from a 
    # jupyter notebook or something instead
    array = ["VELA_v2_art_01_1.0", "VELA_v2_art_01_2.0"]
    read_ion_masses(array)
	print('no script, just for imports')
