import os
from quasarscan.utils.utils import data_path
from quasarscan.data_objects import simulation_quasar_sphere
from quasarscan.preprocessing.code_specific_setup import load_and_setup
from quasarscan.spectra.trident_interfacer import get_line_list
from quasarscan.utils.utils import agora_custom_metals
import trident
import numpy as np

def set_up_folders(path):
    segments = path.split('/')
    partial_path = '/'
    for segment in segments[:-1]:
        partial_path = os.path.join(partial_path,segment)
        if not os.path.exists(partial_path):
            os.mkdir(partial_path)

def run_line(sim_name,sightline_num,ds,start,end,lines,loc = 'default',
             metal_function = None,line_spec_width = 15,dlambda = 0.1):
    redshift = ds.current_redshift
    if loc == 'default':
        root = "~/spectra_from_trident/Ion_Spectra"
    all_done = True
    for line in line_list:
        atom = line[0].split(' ')[0]
        line_fname  = f'{root}/{sim_name}/{redshift}/'+\
                      f'Line_{sightline_num}/{atom}/'+\
                      f'{line[0]}_{line[1]}.txt'
        if os.path.exists(line_fname):
            print(f'line {sightline_num} ({line}) is already done! skipping...')
            continue
        else:
            all_done = False
            break
    if all_done:
        continue
    ray = trident.make_simple_ray(ds,start_position = start, end_position = end, 
                                      data_filename =  "ray.h5",fields = fields_to_keep)
    if metal_function == 'agora':
        trident.add_ion_fields(ray,ion_list,
                            metal_source = 'custom',
                            custom_metal_function = agora_custom_metals,
                            H_source = 'density',
                            custom_H_function = None)
    else:
        trident.add_ion_fields(ray,ion_list)
    for line in line_list:
        atom = line[0].split(' ')[0]
        trident_line = line[0] + " " + str(line[1])
        lambda_min = (line[1]*(1+redshift)) - line_spec_width
        lambda_max = (line[1]*(1+redshift)) + line_spec_width

        sg = trident.SpectrumGenerator(lambda_min = np.floor(lambda_min),
            lambda_max = np.ceil(lambda_max), dlambda = dlambda) 

        sg.make_spectrum(ray, lines= trident_line)

        filename  = f"{root}/{sim_name}/{redshift}/Line_{sightline_num}/{atom}/{line[0]}_{line[1]}.txt"
        set_up_folders(filename)
        sg.save_spectrum(filename)
