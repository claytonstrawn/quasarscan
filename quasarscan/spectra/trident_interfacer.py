import os
from quasarscan.spectra.data_objects import AbsorptionLine
from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
import numpy as np

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

def ion_strongest_line(ion):
    lines_dict = trident_file_reader()
    lines_to_use = lines_dict[ion]
    gammas = []
    for line in lines_to_use:
        #numlines+=1
        #wavelengths+=[line[0]]
        gammas+=[line[1]]
        #fs+=[line[2]]
    max_gamma = np.argmax(gammas)
    return lines_to_use[max_gamma]

def convert_equivalent_width_to_coldens(ew,ew_err,ion,approx_wavelength):
    lines_dict = trident_file_reader()
    lines_to_check = lines_dict[ion]
    for line in lines_to_check:
        wavelengths+=[line[0]]
        gammas+=[line[1]]
        fs+=[line[2]]
    closest_wl_index = np.argmax(np.abs(np.array(wavelengths)-approx_wavelength))
    f = fs[closest_wl_index]
    wl = wavelengths[closest_wl_index]
    if ew < 0.2:
        #formula from https://www.aanda.org/articles/aa/full/2004/04/aa0003/aa0003.right.html
        cdens = 1.13e20*ew/(wl**2*f)
        cdens_low = 1.13e20*(ew-ew_err)/(wl**2*f)
        cdens_high = 1.13e20*(ew+ew_err)/(wl**2*f)
        cdens_eb = np.mean(np.abs(cdens_low-cdens),np.abs(cdens_high-cdens))
    else:
        print('equivalent width not in linear regime, conversion unknown.')
        cdens,cdens_eb = np.nan,np.nan
    return cdens, cdens_eb
        
def get_line_list(atom_list = None, ion_list = None,only_strongest = False):
    if atom_list is None:
        atom_list = []
    if ion_list is None:
        ion_list = []
    dictionary = trident_file_reader()
    for key in dictionary.keys():
        if key not in ion_list and key.split(' ')[0] in atom_list:
            add_ion = key
            ion_list = ion_list + [add_ion]
    line_list = []
    for ion in ion_list:
        for i in range(len(dictionary[ion])):
            if only_strongest:
                if dictionary[ion][i][0] != ion_strongest_line(ion)[0]:
                    continue
            line_list = line_list + [(ion,dictionary[ion][i][0])]
    return atom_list, ion_list, line_list

def trident_lines_starting_points(ions,lines_dict,**kwargs):
    to_return = {}
    if isinstance(ions,str):
        ions = [ions]
    for ion in ions:
        ion_name_nospace = ion.replace(' ','')
        ion_parameters = {'name':ion_name_nospace,
                     'maxN':1e20,'minN':1e8,
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

def call_trident_fitter(wavelength,flux,trident_filename = 'default',
                        line = None,maxNumComps = 8,**kwargs):
    lines_dict = trident_file_reader(filename=trident_filename)
    new_lines_dict = {}
    if line != None:
        ions = [line[0]]
        for ion in ions:
            ion_lines = lines_dict[ion]
            line_parameters = []
            for line_parameter in ion_lines:
                if line[1] == line_parameter[0]:
                    line_parameters = line_parameters + [line_parameter]
            new_lines_dict[ion] = line_parameters
    orderFits,speciesDicts = trident_lines_starting_points(ions,new_lines_dict,**kwargs)
    fitted_lines, fitted_flux = generate_total_fit(wavelength, flux, orderFits,
                                                   speciesDicts,maxLength = 5000,
                                                  )
                                                   #maxNumComps = maxNumComps)
    dict_to_return = {}
    for ion in ions:
        if ion.replace(' ','') not in fitted_lines:
            dict_to_return[ion] = nolinesdict
        else: 
            dict_to_return[ion] = fitted_lines[ion.replace(' ','')]
    return dict_to_return, fitted_flux

def one_line_interpreter(wl,fl,fitted_lines,cosmo_redshift,line,
                         fluxthreshold = 0.99,bv_adjust = None,loud = False):
    ion,nat_wavelength = line
    list_of_AbsorptionLine_objects = []
    lines_dict = fitted_lines[ion]
    for i in range(len(lines_dict['N'])):
        N = fitted_lines[ion]['N'][i]
        b = fitted_lines[ion]['b'][i]
        z = fitted_lines[ion]['z'][i]
        speedoflight = 299792458/1000
        wavelength_detected = nat_wavelength * (1+z)
        rest_wavelength = nat_wavelength * (1+cosmo_redshift)
        doppler_redshift = wavelength_detected/rest_wavelength-1
        velocity = speedoflight*doppler_redshift
        minflux = np.interp(wavelength_detected,wl,fl)
        minflux = minflux if minflux > .01 else 0
        if minflux > fluxthreshold:
            continue
        if loud:
            print(f'best fit for line #{i} at wavelength {wavelength_detected}: '+\
                  f'\n{line}, with wavelength {rest_wavelength}, velocity {velocity}'+\
                  f'column density {N:.2e}')  
        absorption_line = AbsorptionLine(line, cosmo_redshift,bv_adjust = bv_adjust, \
                                         velocity = velocity, N=N, b=b, z=z,min_flux = minflux)
        list_of_AbsorptionLine_objects.append(absorption_line)
    return list_of_AbsorptionLine_objects


