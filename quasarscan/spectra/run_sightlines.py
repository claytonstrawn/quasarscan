import os
from quasarscan.utils.utils import data_path
from quasarscan.data_objects import simulation_quasar_sphere
from quasarscan.preprocessing.code_specific_setup import load_and_setup
from quasarscan.spectra.trident_interfacer import get_line_list
import trident
import numpy as np

def find_endpoints(fullname,redshift,sightline_num,tolerance = .1):
    PATH = data_path()
    foldername = os.path.join(PATH,'output',f'{fullname}coldensinfo')
    if os.path.exists(foldername):
        contents = os.listdir(foldername)
        filename = None
        for filename in contents:
            try:
                z = float(filename.split('z')[1].split('.txt')[0])
            except:
                continue
            if np.abs(z-redshift)<tolerance:
                break
    else:
        print(f'{foldername} does not exist. Try running '\
              '"quasarscan.create_qso_endpoints(fullname, filename, redshift, ions)"')
    full_file_name = os.path.join(foldername,filename)
    readvalsoutput = simulation_quasar_sphere.read_values(full_file_name)
    q = simulation_quasar_sphere.SimQuasarSphere(start_up_info_packet = readvalsoutput)
    assert sightline_num < q.length
    start = q.info[sightline_num][5:8]
    end = q.info[sightline_num][8:11]
    path = q.simparams[6]
    return path,start,end

def create_dirs_in_path(path):
    dirs_in_path = path.split('/')
    last_is_file = '.' in path
    if last_is_file:
        iterate_over = dirs_in_path[:-1]
    else:
        iterate_over = dirs_in_path
    path_so_far = ''
    for d in iterate_over:
        current_name_to_use = os.path.join(path_so_far,d)
        if not os.path.exists(current_name_to_use):
            os.mkdir(current_name_to_use)
        path_so_far = current_name_to_use

def run_line(fullname,approx_redshift,atom_list,sightline_num,wl_distance = 15):
    path,raw_start,raw_end = find_endpoints(fullname,approx_redshift,sightline_num)
    ds,_ = load_and_setup(path,fullname)
    start = ds.arr(raw_start,'unitary')
    end = ds.arr(raw_end,'unitary')
    redshift = ds.current_redshift
    ray = trident.make_simple_ray(ds,start_position = start, end_position = end, 
                              data_filename =  f"ray{sightline_num}.h5", lines = atom_list)
    _, _, line_list = get_line_list(atom_list)
    PATH = data_path()
    root = os.path.join(PATH,"output","spectra_from_trident","all")
    full_path = os.path.join(f"{root}",f"{fullname}",f"{redshift}",f"Line_{sightline_num}")
    create_dirs_in_path(full_path)
    for line in line_list:
        ion,wavelength = line
        atom = ion.split(' ')[0]
        line_name = f'{ion}_{wavelength}'
        lambda_min = wavelength*(1+redshift) - wl_distance
        lambda_max = wavelength*(1+redshift) + wl_distance
        sg = trident.SpectrumGenerator(lambda_min = np.floor(lambda_min),
                                    lambda_max = np.ceil(lambda_max), dlambda = 0.01)
        sg.make_spectrum(ray, lines = [line_name])
        filename  = os.path.join(full_path,f"{atom}",f"{line_name}.txt")
        create_dirs_in_path(filename)
        sg.save_spectrum(filename)
    sg = trident.SpectrumGenerator(lambda_min = 1500,
                                   lambda_max = 5000,
                                   dlambda = 0.01)
    sg.make_spectrum(ray, lines = atom_list)
    filename  = os.path.join(full_path,"full.txt")
    sg.save_spectrum(filename)
