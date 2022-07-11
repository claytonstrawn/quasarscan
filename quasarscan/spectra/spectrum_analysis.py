import numpy as np
import matplotlib.pyplot as plt

def load_file(filename):
    print('running')
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

def plot_wl_around_line(wl,flux,ion,line_wavelength,redshift,noise = 0,
                        left_distance = 20,right_distance = "default"):
    if right_distance == 'default':
        right_distance = left_distance
    
    center = line_wavelength * (1+redshift)
    
    plt.plot([center,center],[0,1],color = 'red',linestyle = ':')
    
    add_noise = np.random.normal(0,noise,len(flux))
    noise_flux = flux + add_noise
    
    plt.plot(wl,noise_flux,label = ion)
    plt.legend()
    plt.ylabel("Relative Flux")
    plt.xlabel("Wavelength (Å)")
    plt.xlim(center-left_distance,center+right_distance)
    
    return

def plot_vel_around_line(wl,flux,ion,line_wavelength,redshift,noise = 0,
                        left_distance = 200,right_distance = "default"):
    
    if right_distance == 'default':
        right_distance = left_distance
    
    speedoflight = 299792458/1000
    v_ion = speedoflight*(wl/(line_wavelength*(1+redshift))-1)

    add_noise = np.random.normal(0,noise,len(flux))
    noise_flux = flux + add_noise
        
    plt.plot(v_ion,noise_flux,label = ion)
    plt.legend()
    plt.ylabel("Relative Flux")
    plt.xlabel("Velocity (km/s)")
    plt.xlim(0-left_distance,0+right_distance)
    
    return

class AbsorptionLine(object):
    #pass
    def __init__(self, line, velocity, min_flux):
        self.line = line
        self.ion = line[0]
        self.rest_wavelength = line[1]
        self.velocity = velocity
        self.min_flux = min_flux
    
    def plot_data(self, ax, line_color_assignments):      
        ax.plot(self.velocity, self.min_flux - 0.05, "o", color = line_color_assignments[self.line])
        
def automatic_component_detector_v2(wl,flux,ion,line_wavelength,redshift,
                        left_distance = 200,right_distance = "default"):
    if right_distance == 'default':
        right_distance = left_distance
    #line = (ion,line_wavelength)
    speedoflight = 299792458/1000
    v_ion = speedoflight*(wl/(line[1]*(1+redshift))-1)
    
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

        if  (flux[maximums[i]] - flux[minimums[i+1]])> 0.0 and (flux[maximums[i]] - flux[minimums[i]]) > extreme_boundary :
            
            max_vion_add = (ion + ' ' + str(line_wavelength), v_ion[maximums[i]])
            max_vion = max_vion + [max_vion_add]
            
            max_flux_add = (ion + ' ' + str(line_wavelength), flux[maximums[i]])
            max_flux = max_flux  + [max_flux_add]
            
            min_vion_add = (ion + ' ' + str(line_wavelength), v_ion[minimums[i]])
            min_vion = min_vion + [min_vion_add]
            
            min_flux_add = (ion + ' ' + str(line_wavelength), flux[minimums[i]])
            min_flux = min_flux  + [min_flux_add]
            
        last_v_min = ((ion + ' ' + line_wavelength), v_ion[minimums[len(maximums)]])
        min_vion = min_vion + [last_v_min]
        last_f_min = ((ion + ' ' + str(line_wavelength)), flux[minimums[len(maximums)]])
        min_flux = min_flux + [last_f_min]

    list_of_lines = []
    for i in range (len(min_vion)):
        list_of_lines_add = AbsorptionLine(line, min_vion[i][1], min_flux[i][1])
        list_of_lines = list_of_lines + [list_of_lines_add]
    return list_of_lines
