import numpy as np
import matplotlib.pyplot as plt
import os
matplotlib_default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
default_color_assignments = {('O VI',1031.912000): matplotlib_default_colors[0],\
                             ('O VI',1037.613000): matplotlib_default_colors[1],\
                             ('C IV',1548.187000): matplotlib_default_colors[2],\
                             ('C IV',1550.772000): matplotlib_default_colors[3]}
from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
#/project/projectdirs/agora/paper_CGM/spectra_from_trident/Ion_Spectra/AGORA_art_CR/1.998998976095549/Line_0/C I.txt
def load_file(filename):
    folders = filename.split('/')
    index_of_enclosing = folders.index('spectra_from_trident')
    folders = folders[index_of_enclosing:]
    simname = folders[2]
    redshift = float(folders[3])
    line_num = int(folders[4].split('_')[1])
    if filename[-4:] == '.txt':
        file = open(filename)
    else:
        file = open(filename + '.txt')
    string_file = file.read()
    array_file = string_file.split('\n')
    array_file.pop(0)
    array_file.pop(len(array_file)-1)
    
    xs = []
    ys = []
    
    for i in range(len(array_file)):
        
        line = array_file[i]
        linelist = line.split(' ')
        xs.append(float(linelist[0]))
        ys.append(float(linelist[2]))
        
    wl = np.array(xs)
    flux = np.array(ys)
    return wl, flux, redshift

def trident_file_reader(filename = 'default'):
    if filename == 'default':
        filename = os.path.expanduser('~/trident/trident/data/line_lists/lines.txt')
    try:
        with open(filename) as f:
            lines = f.readlines()
    except FileNotFoundError as e:
        print(f'File "{filename}" not found. Please specify the location of the trident file ".../lines.txt"')
        raise e
    lines_dict = {}
    for line in lines[1:]:
        elements = line.split()
        ion = f'{elements[0]} {elements[1]}'
        wl = float(elements[2])
        gamma = float(elements[3])
        f = float(elements[4])
        if ion not in lines_dict:
            lines_dict[ion] = [(wl,gamma,f)]
        else:
            lines_dict[ion] += [(wl,gamma,f)]
        lines_dict['fakeion'] = [(1000,1e7,1e-1)]
    return lines_dict

def trident_lines_starting_points(ions,lines_dict,**kwargs):
    to_return = {}
    if isinstance(ions,str):
        ions = [ions]
    for ion in ions:
        ion_name_nospace = ion.replace(' ','')
        ion_parameters = {'name':ion_name_nospace,
                     'maxN':1e17,'minN':1e11,
                     'maxb':300, 'minb':1,
                     'maxz':6, 'minz':0,
                     'init_b':20,
                     'init_N':1e12}
        for kwarg in kwargs:
            ion_parameters[kwarg] = kwargs[kwarg]
        fs,gammas,wavelengths,numlines = [],[],[],0
        try:
            lines_to_use = lines_dict[ion]
        except KeyError:
            continue
        for line in lines_to_use:
            numlines+=1
            wavelengths+=[line[0]]
            gammas+=[line[1]]
            fs+=[line[2]]
        ion_parameters['f'] = fs
        ion_parameters['Gamma'] = gammas
        ion_parameters['wavelength'] = wavelengths
        ion_parameters['numLines'] = numlines
        to_return[ion_name_nospace] = ion_parameters
    return list(to_return.keys()),to_return

nolinesdict = {'N': np.array([], dtype=np.float64),\
               'b': np.array([], dtype=np.float64),\
               'z': np.array([], dtype=np.float64),\
               'group#': np.array([], dtype=np.float64)}

def call_trident_fitter(wavelength,flux,ions,filename = 'default',**kwargs):
    lines_dict = trident_file_reader(filename=filename)
    orderFits,speciesDicts = trident_lines_starting_points(ions,lines_dict,**kwargs)
    fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits, speciesDicts)
    dict_to_return = {}
    for ion in ions:
        if ion.replace(' ','') not in fitted_lines:
            dict_to_return[ion] = nolinesdict
        else: 
            dict_to_return[ion] = fitted_lines[ion.replace(' ','')]
    return dict_to_return, fitted_flux

def plot_wl_around_line(wl,flux,line,redshift,noise = 0,color = 'default',
                        left_distance = 20,right_distance = "default",
                        label = 'default',ax = None,plot_chunk = 'all'):
    if right_distance == 'default':
        right_distance = left_distance
    if ax is None:
        _,ax = plt.subplots()
    
    center = line[1] * (1+redshift)
    
    ax.plot([center,center],[0,1],color = 'red',linestyle = ':')
    
    add_noise = np.random.normal(0,noise,len(flux))
    noise_flux = flux + add_noise
    
    if color =='default':
        if line in default_color_assignments:
            color = default_color_assignments[line]
        else:
            color = None
    if label == 'default':
        label = str(line)
    if plot_chunk == 'all':
        chunk = np.ones(len(wl),dtype = bool)
    elif plot_chunk == 'lims':
        chunk = np.logical_and(wl>=center-left_distance,wl<=center+right_distance)
    elif isinstance(plot_chunk,(int,float)):
        chunk = np.logical_and(wl>=center-left_distance*plot_chunk,wl<=center+right_distance*plot_chunk)
    ax.plot(wl[chunk],noise_flux[chunk],label = label,color = color)
    ax.legend()
    ax.set_ylabel("Relative Flux")
    ax.set_xlabel("Wavelength (Å)")
    ax.set_xlim(center-left_distance,center+right_distance)

def plot_vel_around_line(wl,flux,line,redshift,noise = 0,color = 'default',
                        left_distance = 200,right_distance = "default",\
                        label = 'default',ax = None,plot_chunk = 'lims'):
    
    if right_distance == 'default':
        right_distance = left_distance
    if ax is None:
        _,ax = plt.subplots()
    speedoflight = 299792458/1000
    v_ion = speedoflight*(wl/(line[1]*(1+redshift))-1)

    add_noise = np.random.normal(0,noise,len(flux))
    noise_flux = flux + add_noise
    
    if color =='default':
        if line in default_color_assignments:
            color = default_color_assignments[line]
        else:
            color = None
    if label == 'default':
        label = str(line)
    if plot_chunk == 'all':
        chunk = np.ones(len(wl),dtype = bool)
    elif plot_chunk == 'lims':
        chunk = np.logical_and(v_ion>=0-left_distance,v_ion<=0+right_distance)
    elif isinstance(plot_chunk,(int,float)):
        chunk = np.logical_and(v_ion>=0-left_distance*plot_chunk,v_ion<=0+right_distance*plot_chunk)
    ax.plot(v_ion[chunk],noise_flux[chunk],label = label,color = color)
    ax.legend()
    ax.set_ylabel("Relative Flux")
    ax.set_xlabel("Velocity (km/s)")
    ax.set_xlim(0-left_distance,0+right_distance)

    

class AbsorptionLine(object):
    
    def __init__(self, line, velocity, min_flux):
        self.line = line
        self.ion = line[0]
        self.rest_wavelength = line[1]
        self.velocity = velocity
        self.min_flux = min_flux
    
    def plot_data(self, ax, color = 'default'): 
        if color == 'default':
            if self.line in default_color_assignments:
                color = default_color_assignments[self.line]
            else:
                color = None
        if self.min_flux is None:
            flux_to_plot = 1.1
        else:
            flux_to_plot = self.min_flux - 0.05
        ax.plot(self.velocity, flux_to_plot, "o", color = color)
        
class Component(object):
    def __init__(self, list_of_lines):
        self.list_of_lines = list_of_lines
       # ion_list = []
        self.min_flux = np.inf
        total_velocity = 0
        for line in list_of_lines:
          #  ion_list = ion_list + [(line.ion,line.rest_wavelength)]
            total_velocity = total_velocity + line.velocity
            new_flux = line.min_flux
            if(new_flux < self.min_flux):
                self.min_flux = new_flux
        self.velocity = total_velocity/len(list_of_lines) 
        
    def alignment_printer(self):
        print('There is a component at ' + str(self.velocity) + ' km/s with ions: ' , self.list_of_lines)
        
    def plot_data(self, ax, color = 'black'): 
        flux_to_plot = self.min_flux - 0.1
        ax.plot([self.velocity - 1, self.velocity + 1], [flux_to_plot, flux_to_plot], linewidth = 4, color = color)
        
    
            
def trident_component_interpreter(dict_of_trident_components):
    return list_of_AbsorptionLine_objects
        
def automatic_component_detector_v2(wl,flux,line,redshift,
                        left_distance = 200,right_distance = "default", extreme_boundary = 0.01):
    if right_distance == 'default':
        right_distance = left_distance
    #line = (ion,line_wavelength)
    speedoflight = 299792458/1000
    v_ion = speedoflight*(wl/(line[1]*(1+redshift))-1)
    
    diff_flux = np.gradient(flux)
    
    minimums = []
    maximums = []
    
    for i in range(1, len(wl[1:])): 
        if v_ion[i] > (0 - left_distance) and v_ion[i] < (0 + right_distance):
            if diff_flux[i]>0 and diff_flux[i-1]<0:
                #minimum
                minimums = minimums + [i]
                
            if diff_flux[i]<0 and diff_flux[i-1]>0:
                #maximum
                maximums = maximums + [i]
                
    #print(len(minimums))
    #print(len(maximums))
    min_vion = []
    max_vion = []
    min_flux = []
    max_flux = []
    
   
    for i in range(len(maximums)):

        if (flux[maximums[i]] - flux[minimums[i+1]])> 0.0 and \
         (flux[maximums[i]] - flux[minimums[i]]) > extreme_boundary :
            
            max_vion_add = (line[0] + ' ' + str(line[1]), v_ion[maximums[i]])
            max_vion = max_vion + [max_vion_add]
            
            max_flux_add = (line[0] + ' ' + str(line[1]), flux[maximums[i]])
            max_flux = max_flux  + [max_flux_add]
            
            min_vion_add = (line[0] + ' ' + str(line[1]), v_ion[minimums[i]])
            min_vion = min_vion + [min_vion_add]
            
            min_flux_add = (line[0] + ' ' + str(line[1]), flux[minimums[i]])
            min_flux = min_flux  + [min_flux_add]
    
    if flux[minimums[len(maximums)]] < 0.98: 
        last_v_min = ((line[0] + ' ' + str(line[1])), v_ion[minimums[len(maximums)]])
        min_vion = min_vion + [last_v_min]
        last_f_min = ((line[0] + ' ' + str(line[1])), flux[minimums[len(maximums)]])
        min_flux = min_flux + [last_f_min]

   # print('\n')
    #print('min_vions ' + line[0] + ' ' + str(line[1]) + ': ')
   # print(min_vion)
    #print('\n')
    #print('min_flux ' + line[0] + ' ' + str(line[1]) + ': ')
    #print(min_flux)
    #print('\n')
    #print('max_vions ' + line[0] + ' ' + str(line[1]) + ': ')
    #print(max_vion)     
    #print('\n')
    #print('max_flux '  + line[0] + ' ' + str(line[1]) + ': ')
    #print(max_flux)   
    #print('\n')

    list_of_lines = []
    for i in range (len(min_vion)):
        list_of_lines_add = AbsorptionLine(line, min_vion[i][1], min_flux[i][1])
        list_of_lines = list_of_lines + [list_of_lines_add]
        
    print(list_of_lines)
    return list_of_lines

def alignment_checker(v_min, v_max, bin_step, list_of_lines):
    bins = list(np.arange(v_min, v_max, bin_step))
    counter = [[]] * len(bins)
    list_of_comps = []
    total_velocity = 0
    for i in range(1, len(bins)):
        for j in range(len(list_of_lines)):
            if(list_of_lines[j].velocity > bins[i-1] and list_of_lines[j].velocity < bins[i]):
                counter[i] = counter[i] + [list_of_lines[j]]
        if(len(counter[i]) != 0):
            add_comp = Component(counter[i])
            list_of_comps = list_of_comps + [add_comp]
    return list_of_comps