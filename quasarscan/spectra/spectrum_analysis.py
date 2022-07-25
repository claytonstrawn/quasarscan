import numpy as np
import matplotlib.pyplot as plt
import os

speedoflight = 299792458/1000  #units are km/s (so 1/1000 of the m/s value)

matplotlib_default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
default_color_assignments = {('O VI',1031.912000): matplotlib_default_colors[0],\
                             ('O VI',1037.613000): matplotlib_default_colors[1],\
                             ('C IV',1548.187000): matplotlib_default_colors[2],\
                             ('C IV',1550.772000): matplotlib_default_colors[3]}
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



