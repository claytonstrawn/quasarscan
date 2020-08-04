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
    #list_of_galaxies = [('VELA_v2_art_10',1.0),('VELA_v2_art_10',2.0),('VELA_v2_art_11',1.0),('VELA_v2_art_11',2.0)]
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

def read_ion_number_densities(list_of_galaxies):
    return array

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
