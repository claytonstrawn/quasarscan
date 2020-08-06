#***import statements***
import numpy as np
import matplotlib.pyplot as plt
import os
# from quasarscan import parse_metadata
import parse_metadata
import os.path
from os import path

#define analysis functions

def get_sorted_array(criteria):
    #skip this one for now
    pass
    
def read_ion_masses(list_of_galaxies):
    #list_of_galaxies = [('VELA_v2_art_10',1.0),('VELA_v2_art_10',2.0),('VELA_v2_art_11',1.0),('VELA_v2_art_11',2.0)]
    num_exist = 0
    gal = list_of_galaxies
    list_of_galaxies = []
    for i in range(len(gal)):
        filename = 'quasarscan/ion_state_ytanalysis/' + gal[i] + '/o_ion_mass_data.txt'
        if (path.exists(filename)):
            list_of_galaxies.append(gal[i])
        else: 
            print(gal[i] + " file did not exist")
    PI = np.empty((len(list_of_galaxies),9))
    CI = np.empty((len(list_of_galaxies),9))
    total = np.empty((len(list_of_galaxies),9))
    mass = np.empty(len(list_of_galaxies))
    for i in range(len(list_of_galaxies)):
        filename = 'quasarscan/ion_state_ytanalysis/' + list_of_galaxies[i] + '/o_ion_mass_data.txt'
        f = open(filename, "r")
        contents = f.readline()
        arr = np.fromstring(contents[9:-1],dtype=float,sep=', ')
        for j in range(9): 
            PI[i,j] = arr[j]
        f.readline()
        contents = f.readline()
        arr = np.fromstring(contents[10:-1],dtype=float,sep=', ')
        for j in range(9): 
            CI[i,j] = arr[j]
#         if(PI[i,7] < CI[i,7]):
#             print(list_of_galaxies[i])
        f.readline()
        contents = f.readline()
        arr = np.fromstring(contents[13:-1],dtype=float,sep=', ')
        for j in range(9): 
            total[i,j] = arr[j]
        f.readline()
        f.readline()
        f.readline()
        contents = f.readline()
        mass[i] = float(contents[7:])
    return PI, CI, total, mass


def combined_mass_plots(list_of_galaxies,weighting = 'mass',logy=True, plot_data = 'all'):
    #possible weightings = 'mass' (add all the stacked masses) 
    #and 'snapshot' (average all the fractions for individual galaxies)
    #possible plot_datas = 'all' (plot PI, CI, and total)
    #and 'PI', 'CI', or 'total' if you just want one 
    #processing of the data array (list_of_galaxies should generate captions or something)
    PI, CI, total, mass = read_ion_masses(list_of_galaxies)
    x=['O I', 'O II', 'O III', 'O IV', 'O V', 'O VI', 'O VII', 'O VIII', 'O IX']
    x1=['O I', 'O II', 'O III', 'O IV', 'O V', 'O VI', 'O VII', 'O VIII', 'O IX']
    PI_avg = np.empty(9)
    CI_avg = np.empty(9)
    total_avg = np.empty(9)
    mass_sum = np.sum(mass)
    if (weighting == 'mass'):
        for i in range(len(mass)):
            for j in range(9):
                PI_avg[j] += PI[i,j]*mass[i] / mass_sum
                CI_avg[j] += CI[i,j]*mass[i] / mass_sum
                total_avg[j] += total[i,j]*mass[i] / mass_sum
        x1 = x1[1:]
        PI_avg = PI_avg[1:]
    elif (weighting == 'snapshot'):
        PI_avg = np.average(PI, axis=0)
        CI_avg = np.average(CI, axis=0)
        total_avg = np.average(total, axis=0)
    plt.clf()
    if (plot_data == 'all' or plot_data == 'CI'):
        plt.plot(x, np.log10(CI_avg), label="CI oxygen mass",color='r',linewidth=1,marker='o',linestyle='-')
    if (plot_data == 'all' or plot_data == 'PI'):
        plt.plot(x1, np.log10(PI_avg), label="PI oxygen mass",color='b',linewidth=1,marker='o',linestyle='-')
    if (plot_data == 'all' or plot_data == 'total'):
        plt.plot(x, np.log10(total_avg), label="total oxygen mass",color='g',linewidth=1,marker='o',linestyle='-')
    plt.legend(loc=0, fontsize=10)
    plt.xlabel("Oxygen ion")
    plt.ylabel("ion fraction of total oxygen mass in log10")
    plt.plot()
    plt.show()

def read_ion_number_densities(list_of_galaxies):
    for galaxy in range(len(list_of_galaxies)):
        filename = 'quasarscan/ion_state_ytanalysis/' + list_of_galaxies[i] + '/o_number_density_data.txt'
        f = open(filename, "r")
        lines = f.readlines()
        current_ion_state = None
        all_data = np.zeros((3,9,16))
        for line in lines:
            if line == ' \n':
                continue
            if 'O' in line:
                current_ion_state = line.strip('\n')
                continue
            arys = line.split('], ')

            for k,ary in enumerate(arys):
                stripped_array = ary.replace('\n', '').replace('[','').replace(']', '')
                key = stripped_array.split(', ')[0]
                value = float(stripped_array.split(', ')[1])
                if current_ion_state.startswith('O'):
                    i = 0
                elif current_ion_state.startswith('P'):
                    i = 1
                else:
                    i = 2
                
                numbers = {'I': 0, 'II':1, 'III':2, 'IV':3, 'V':4, 'VI':5, 'VII':6, 'VIII':7, 'IX':8}
                j = numbers[current_ion_state.split('O')[1]]
                
                all_data[i,j,k] = (all_data[i,j,k]*galaxy+value)/(galaxy+1)

    all_data = all_data / 100
    return all_data

def covering_fraction_plots(list_of_galaxies,thresholds,plot_data):
    #list_of_galaxies in form (ex) ['VELA_v2_art_10_1.0','VELA_v2_art_1_1.0',...]
    #thresholds ex) [13,14,15]
    #possible plot_data = 'all' (plot PI, CI, and total) and 'PI', 'CI', or 'total' if you just want one 
    #and 'snapshot' (average all the fractions for individual galaxies)
    #processing of the data array (list_of_galaxies should generate captions or something)

    # 1st dimension: types of oxygen
    # 2nd dimension: oxygen state (ionization)
    # 3rd dimension: density value

    all_data = read_ion_number_densities(list_of_galaxies)
    labels = ['OI','OII','OIII','OIV','OV','OVI','OVII','OVIII','OIX']
    output = []
    for threshold in thresholds:
        values = []
        for i in range(all_data.shape[0]):
            oxygen_type = []
            for j in range(all_data.shape[1]):
                oxygen_type.append(sum(all_data[i][j][threshold-3:all_data.shape[2]-1]))
            values.append(oxygen_type)
        output.append(values)

    plt.figure(figsize=(20,20))
    colors = ['blue','red','green','black','purple','gray','yellow','cyan','magenta']
    for i in range(len(thresholds)):
        if plot_data == "all" or plot_data == "total":
            plt.plot(labels, output[i][0], linestyle = 'solid', color = colors[i], linewidth = 3, label='Total Oxygen_'+str(thresholds[i]))
            plt.scatter(labels,output[i][0], color = colors[i])
        if plot_data == "all" or plot_data == "PI":
            plt.plot(labels, output[i][1], linestyle = 'dashed', color = colors[i], linewidth = 3, label='PI_'+str(thresholds[i]))
            plt.scatter(labels,output[i][1], color = colors[i])
        if plot_data == "all" or plot_data == "CI":
            plt.plot(labels, output[i][2], linestyle = 'dotted', color = colors[i], linewidth = 3, label='CI_'+str(thresholds[i]))
            plt.scatter(labels,output[i][2], color = colors[i])

    plt.legend(prop={'size': 25})
    plt.xlabel('Oxygen States', fontsize=15)
    plt.ylabel('Covering Fraction', fontsize=15)
    plt.show() 

if __name__ == '__main__':
	# this one doesn't really make sense as a script
    # I'd probably just import these functions from a 
    # jupyter notebook or something instead
    array = ["VELA_v2_art_01_1.0", "VELA_v2_art_01_2.0"]
    data_array = read_ion_masses(array)
    combined_mass_plots(array,data_array)
    print('no script, just for imports')
