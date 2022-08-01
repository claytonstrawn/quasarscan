import numpy as np
from quasarscan.spectra.spectrum_analysis import default_color_assignments,\
                                                    matplotlib_default_colors,\
                                                    speedoflight

class AbsorptionLine(object):
    def __init__(self, line, z_cos, z_tot = None, velocity = None,\
                 wavelength_detected = None, min_flux = None, N = None, \
                 b = None, z = None):
        #need at least one of these
        self.line = line
        self.ion = line[0]
        self.wavelength = line[1]
        self.wl_nat = self.wavelength
        self.z_cos = z_cos
        self.wl_rest = self.wavelength*(1+self.z_cos)
        assert (z_tot is not None) or (velocity is not None) or (wavelength_detected is not None)
        if z_tot is not None:
            self.z = z_tot
            self.wl_det = self.wl_nat*(1+self.z)
            self.velocity = speedoflight*(self.wl_det/self.wl_rest-1)
        elif velocity is not None:
            self.velocity = velocity
            self.wl_det = self.wl_rest*(1+self.velocity/speedoflight)
            self.z = self.wl_det/self.wl_nat-1
        elif wavelength_detected is not None:
            self.wl_det = wavelength_detected
            self.z = self.wl_det/self.wl_nat-1
            self.velocity = speedoflight*(self.wl_det/self.wl_rest-1)
        
        #these four should be well defined, but if they aren't, its ok
        self.min_flux = min_flux
        self.N = N
        self.b = b
        self.z = z
    
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
        ax.plot(self.velocity, flux_to_plot, "^", color = color)
        
class Component(object):
    def __init__(self, list_of_lines,component_threshold = None):
        self.list_of_lines = list_of_lines
        self.min_flux = np.inf
        total_velocity = 0
        for line in list_of_lines:
            total_velocity = total_velocity + line.velocity
            new_flux = line.min_flux
            if new_flux is None:
                self.min_flux = None
                break
            if(new_flux < self.min_flux):
                self.min_flux = new_flux
        self.velocity = total_velocity/len(list_of_lines) 
        self.z_cos = None
        for line in list_of_lines:
            if self.z_cos is None:
                self.z_cos = line.z_cos
            else:
                assert self.z_cos == line.z_cos
        self.component_threshold = component_threshold
    
    def __contains__(self, item):
        ion, wl = item
        for line in self.list_of_lines:
            if ion == line.ion and wl == line.wavelength:
                if self.component_threshold is None:
                    return True
                elif self.component_threshold < line.N:
                    return True
        return False
        
    def alignment_printer(self):
        print('There is a component at ' + str(self.velocity) + \
              ' km/s with ions: ', print_ion(self.list_of_lines))
        
    def plot_data(self, ax, color = 'black'): 
        if self.min_flux is None:
            flux_to_plot = 1.2
        else:
            flux_to_plot = self.min_flux - 0.1
        ax.plot([self.velocity - 2,self.velocity - 1,self.velocity, self.velocity + 1,self.velocity + 2], 
                [flux_to_plot]*5, 's',markersize = 2, color = color)
        
def print_ion(list_of_lines):
    lines = []
    for i in range(len(list_of_lines)):
        add_line = str(list_of_lines[i].line)
        lines = lines + [add_line]
    return str(lines)

def alignment_checker(list_of_lines, acceptable_width = 4):
    x = len(list_of_lines)
    swapped = False
    for i in range(x-1):
        for j in range(0, x-i-1):
            if (list_of_lines[j].velocity > list_of_lines[j+1].velocity):
                swapped = True
                list_of_lines[j], list_of_lines[j+1] = list_of_lines[j+1], list_of_lines[j]
        if not swapped:
            break
    counter = [[]] * len(list_of_lines)
    list_of_comps = []
    total_velocity = 0
    c = 0
    n = 0
    for i in range(len(list_of_lines)):
        i = n 
        j = n
        while (j < len(list_of_lines)):
            if(abs(list_of_lines[i].velocity - list_of_lines[j].velocity) <= acceptable_width):
                counter[c] = counter[c] + [list_of_lines[j]]
                n = j
            j = j+1
        if(len(counter[c]) != 0):
            add_comp = Component(counter[c])
            list_of_comps = list_of_comps + [add_comp]
        c = c + 1
        n = n + 1
    return list_of_comps

def automatic_component_detector_v2(wl,flux,line,redshift,
                        left_distance = 200,right_distance = "default", extreme_boundary = 0.01):
    if right_distance == 'default':
        right_distance = left_distance
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

    list_of_lines = []
    for i in range (len(min_vion)):
        list_of_lines_add = AbsorptionLine(line, min_vion[i][1], min_flux[i][1])
        list_of_lines = list_of_lines + [list_of_lines_add]
        
    return list_of_lines

def load_components(loc):
    read_comp_list = []
    with open(loc,'r') as f:
        lines = f.readlines()
    x = 0
    for i in range(len(lines)):
        i = x
        if(i >= len(lines)):
            break
        line = lines[i]
        element = line.split('_')
        if(element[0] == 'Component'):
            comp_line_list = []
            x = i + 1
            line_elements = lines[x].split('_')   
            while(line_elements[0] != 'Component'):
                comp_line_list_add = AbsorptionLine((line_elements[1], float(line_elements[2])), float(line_elements[3]),
                                velocity = float(line_elements[4]), min_flux = float(line_elements[5]), N = float(line_elements[6]),
                                                    b = float(line_elements[7]), z = float(line_elements[8]))
                comp_line_list = comp_line_list + [comp_line_list_add]
                x = x+1
                if(x >= len(lines)):
                    break
                line_elements = lines[x].split('_')
        make_comp_add = Component(comp_line_list)
        make_comp_add.velocity = float(element[1])
        make_comp_add.min_flux = float(element[2])
        make_comp_add.z_cos = float(element[3])
        read_comp_list = read_comp_list + [make_comp_add]
        #create the Components back from the text in the lines as a list
    return read_comp_list
