import numpy as np
import matplotlib.pyplot as plt
import os
from quasarscan.utils.utils import data_path

speedoflight = 299792458/1000  #units are km/s (so 1/1000 of the m/s value)

matplotlib_default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

default_color_assignments = {
                             ('O VI',1031.912000): matplotlib_default_colors[0],\
                             ('O VI',1037.613000): matplotlib_default_colors[5],\
                             ('C IV',1548.187000): matplotlib_default_colors[1],\
                             ('C IV',1550.772000): matplotlib_default_colors[6],\
                             ('Ne VIII', 780.324000): matplotlib_default_colors[7],\
                             ('Ne VIII', 770.409000): matplotlib_default_colors[2],\
                             ('Mg X', 624.941000): matplotlib_default_colors[8],\
                             ('Mg X', 609.793000): matplotlib_default_colors[3],\
                             ('Si IV',1402.770000): matplotlib_default_colors[4],\
                            }

#/project/projectdirs/agora/paper_CGM/spectra_from_trident/Ion_Spectra/AGORA_art_CR/1.998998976095549/Line_0/C I.txt

def load_sightline_from_qsdata(fullname,approx_redshift,sightline_num,line=None,tolerance = 0.1):
    PATH = data_path()
    root = os.path.join(PATH,"output","spectra_from_trident","all")
    for z in os.listdir(os.path.join(root,fullname)):
        if np.abs(float(z)-approx_redshift)<tolerance:
            redshift = z
    full_path = os.path.join(f"{root}",f"{fullname}",f"{redshift}",f"Line_{sightline_num}")
    if line is not None:
        ion,wavelength = line
        atom = ion.split(' ')[0]
        line_name = f'{ion}_{wavelength}'
        filename  = os.path.join(full_path,f"{atom}",f"{line_name}.txt")
        if not os.path.exists(filename):
            one_level_up = os.path.join(full_path,f"{atom}")
            print(f"file '{filename}' does not exist. Lines for this ion are:")
            options = os.listdir(one_level_up)
            for option in options:
                print('\t',option)
    else:
        filename  = os.path.join(full_path,"full.txt")
    return load_file(filename)

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
    list_file = string_file.split('\n')
    list_file.pop(0)
    list_file.pop(len(list_file)-1)
    xs = []
    ys = []
    for i in range(len(list_file)):
        line = list_file[i]
        linelist = line.split(' ')
        xs.append(float(linelist[0]))
        ys.append(float(linelist[2]))
    wl = np.array(xs)
    flux = np.array(ys)
    return wl, flux, redshift

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
    return ax

def plot_vel_around_line(wl,fl,line,cosmo_redshift,noise = 0,color = 'default',
                        left_distance = 200,right_distance = "default",\
                        label = 'default',ax = None,plot_chunk = 'lims',
                        bv_adjust = None,fitted = False):
    
    if right_distance == 'default':
        right_distance = left_distance
    if ax is None:
        _,ax = plt.subplots()
    v_ion = speedoflight*(wl/(line[1]*(1+cosmo_redshift))-1)
    if bv_adjust is not None:
        v_ion = v_ion+bv_adjust

    add_noise = np.random.normal(0,noise,len(fl))
    noise_flux = fl + add_noise
    
    if color =='default':
        if line in default_color_assignments:
            color = default_color_assignments[line]
        else:
            color = None
    if label == 'default':
        label = f'{line[0]} ({line[1]:.3f})'
    if plot_chunk == 'all':
        chunk = np.ones(len(wl),dtype = bool)
    elif plot_chunk == 'lims':
        chunk = np.logical_and(v_ion>=0-left_distance,v_ion<=0+right_distance)
    elif isinstance(plot_chunk,(int,float)):
        chunk = np.logical_and(v_ion>=0-left_distance*plot_chunk,v_ion<=0+right_distance*plot_chunk)
    if fitted:
        linewidth = 4
        linestyle = '-'
        alpha = 0.5
        label = None
    else:
        linewidth = 2
        linestyle = '-'
        alpha = 1
    ax.plot(v_ion[chunk],noise_flux[chunk],label = label,color = color,
            linestyle = linestyle, linewidth = linewidth, alpha = alpha)
    ax.legend()
    ax.set_ylabel("Relative Flux")
    ax.set_xlabel("Velocity (km/s)")
    ax.set_xlim(0-left_distance,0+right_distance)
    return ax



