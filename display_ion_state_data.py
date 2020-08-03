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
    array = np.empty((len(list_of_galaxies),27))
    for i in range(len(list_of_galaxies)):
        filename = 'quasarscan/ion_state_ytanalysis/' + list_of_galaxies[i] + '/o_ion_mass_data.txt'
        f = open(filename, "r")
        contents = f.readline()
        curindex = contents.index('[') + 1
        for j in range(8):
            array[i,j] = float(contents[curindex:contents.index(', ',curindex)])
            curindex = contents.index(',',curindex + 1) + 2
        array[i,8] = float(contents[curindex:contents.index(']',curindex)])
        f.readline()
        contents = f.readline()
        curindex = contents.index('[') + 1
        for j in range(8):
            array[i,j+9] = float(contents[curindex:contents.index(', ',curindex)])
            curindex = contents.index(',',curindex + 1) + 2
        array[i,17] = float(contents[curindex:contents.index(']',curindex)])
        f.readline()
        contents = f.readline()
        curindex = contents.index('[') + 1
        for j in range(8):
            array[i,j+18] = float(contents[curindex:contents.index(', ',curindex)])
            curindex = contents.index(',',curindex + 1) + 2
        array[i,26] = float(contents[curindex:contents.index(']',curindex)])
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
    avg = np.average(data_array,axis=0)
    x=['O I', 'O II', 'O III', 'O IV', 'O V', 'O VI', 'O VII', 'O VIII', 'O IX']
    plt.clf()
    plt.plot(x, np.log10(avg[:9]), label="PI oxygen mass",color='b',linewidth=1,marker='o',linestyle='-')
    plt.plot(x, np.log10(avg[9:18]), label="CI oxygen mass",color='r',linewidth=1,marker='o',linestyle='-')
    plt.plot(x, np.log10(avg[18:]), label="total oxygen mass",color='g',linewidth=1,marker='o',linestyle='-')
    plt.legend(loc=0, fontsize=10)
    plt.xlabel("Oxygen ion")
    plt.ylabel("ion fraction of total oxygen mass in log10")
    plt.plot()
    plt.show()

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
    data_array = read_ion_masses(array)
    combined_mass_plots(array,data_array)
	print('no script, just for imports')
