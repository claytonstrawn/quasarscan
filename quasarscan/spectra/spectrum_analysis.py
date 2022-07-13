import numpy as np
import matplotlib.pyplot as plt
matplotlib_default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
default_color_assignments = {('O VI',1031.912000): matplotlib_default_colors[0],\
                             ('O VI',1037.613000): matplotlib_default_colors[1],\
                             ('C IV',1548.187000): matplotlib_default_colors[2],\
                             ('C IV',1550.772000): matplotlib_default_colors[3]}
from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit

def load_file(filename):
    z_redshift = filename.split('+')[2]
    redshift = float(z_redshift.split('z')[1])
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
    if filename = 'default':
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
    for ion in ions:
        if ion.replace(' ','') not in fitted_lines:
            fitted_lines[ion.replace(' ','')] = nolinesdict
    return fitted_lines, fitted_flux

def plot_wl_around_line(wl,flux,ion,line_wavelength,redshift,noise = 0,color = 'default',
                        left_distance = 20,right_distance = "default",ax = None):
    if right_distance == 'default':
        right_distance = left_distance
    if ax is None:
        _,ax = plt.subplots()
    
    center = line_wavelength * (1+redshift)
    
    ax.plot([center,center],[0,1],color = 'red',linestyle = ':')
    
    add_noise = np.random.normal(0,noise,len(flux))
    noise_flux = flux + add_noise
    
    if color =='default':
        if (ion,line_wavelength) in default_color_assignments:
            color = default_color_assignments[(ion,line_wavelength)]
        else:
            color = None
    ax.plot(wl,noise_flux,label = ion,color = color)
    ax.legend()
    ax.set_ylabel("Relative Flux")
    ax.set_xlabel("Wavelength (Å)")
    ax.set_xlim(center-left_distance,center+right_distance)

def plot_vel_around_line(wl,flux,ion,line_wavelength,redshift,noise = 0,color = 'default',
                        left_distance = 200,right_distance = "default",ax = None):
    
    if right_distance == 'default':
        right_distance = left_distance
    if ax is None:
        _,ax = plt.subplots()
    speedoflight = 299792458/1000
    v_ion = speedoflight*(wl/(line_wavelength*(1+redshift))-1)

    add_noise = np.random.normal(0,noise,len(flux))
    noise_flux = flux + add_noise
    
    if color =='default':
        if (ion,line_wavelength) in default_color_assignments:
            color = default_color_assignments[(ion,line_wavelength)]
        else:
            color = None
    ax.plot(v_ion,noise_flux,label = ion,color = color)
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
        ax.plot(self.velocity, self.min_flux - 0.05, "o", color = color)
        
class Component(object):
    def __init__(self, velocity, list_of_lines):
        pass
    def print_out_lines(self):
        pass
        
def automatic_component_detector_v2(wl,flux,ion,line_wavelength,redshift,
                        left_distance = 200,right_distance = "default"):
    if right_distance == 'default':
        right_distance = left_distance
    #line = (ion,line_wavelength)
    speedoflight = 299792458/1000
    v_ion = speedoflight*(wl/(line_wavelength*(1+redshift))-1)
    
    diff_flux = np.gradient(flux)
    
    absorption_cutoff = 0.99
    extreme_boundary = 0.05
    
    minimums = []
    maximums = []
    
    for i in range(1, len(wl[1:])): 
        if v_ion[i] > (0 - left_distance) and v_ion[i] < (0 + right_distance):
            if diff_flux[i]>0 and diff_flux[i-1]<0 and flux[i]<absorption_cutoff:
                #minimum
                minimums = minimums + [i]
                
            if diff_flux[i]<0 and diff_flux[i-1]>0 and flux[i]<absorption_cutoff:
                #maximum
                maximums = maximums + [i]
                
    
    min_vion = []
    max_vion = []
    min_flux = []
    max_flux = []
    
    for i in range(len(maximums)):

        if  (flux[maximums[i]] - flux[minimums[i+1]])> 0.0 and \
            (flux[maximums[i]] - flux[minimums[i]]) > extreme_boundary :
            
            max_vion_add = (ion + ' ' + str(line_wavelength), v_ion[maximums[i]])
            max_vion = max_vion + [max_vion_add]
            
            max_flux_add = (ion + ' ' + str(line_wavelength), flux[maximums[i]])
            max_flux = max_flux  + [max_flux_add]
            
            min_vion_add = (ion + ' ' + str(line_wavelength), v_ion[minimums[i]])
            min_vion = min_vion + [min_vion_add]
            
            min_flux_add = (ion + ' ' + str(line_wavelength), flux[minimums[i]])
            min_flux = min_flux  + [min_flux_add]
            
        last_v_min = ((ion + ' ' + str(line_wavelength)), v_ion[minimums[len(maximums)]])
        min_vion = min_vion + [last_v_min]
        last_f_min = ((ion + ' ' + str(line_wavelength)), flux[minimums[len(maximums)]])
        min_flux = min_flux + [last_f_min]

    list_of_lines = []
    for i in range (len(min_vion)):
        list_of_lines_add = AbsorptionLine((ion,line_wavelength), min_vion[i][1], min_flux[i][1])
        list_of_lines = list_of_lines + [list_of_lines_add]
    return list_of_lines
