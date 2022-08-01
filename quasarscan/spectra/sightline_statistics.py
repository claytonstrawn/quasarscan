import numpy as np
from quasarscan.spectra.data_objects import load_components

#for one code:
#iterate over all the lines
#if any lines have the ion
#count how many components within that line
#but only counting ones where this ion is found
#return the average number

def comps_fraction(ion,comp_list):
    line_comp_tot = 0
    comp_fraction = 0
    for comp in comp_list:
        for line in comp.list_of_lines:
            if ion == line.ion:
                line_comp_tot = line_comp_tot + 1
                break
    if len(comp_list) != 0:
        comp_fraction = line_comp_tot/len(comp_list)
    return comp_fraction

def ncomps(ion,comp_list):
    line_comp_tot = 0
    for comp in comp_list:
        for line in comp.list_of_lines:
            if ion == line.ion:
                line_comp_tot = line_comp_tot + 1
                break
    return line_comp_tot

#if any lines have the ion
#collect b for each component within that line
#but only counting ones where this ion is found
#return the average b
def bparam(ion, comp_list):
    total_b = 0
    line_b_avg = 0
    counter = 0
    for comp in comp_list:
        for line in comp.list_of_lines:
            if ion == line.ion:
                counter += 1
                total_b = total_b + line.b
    if counter != 0:
        line_b_avg = total_b / counter
    return line_b_avg

def sightline_looping(code,sightline_range,ions):
    root = '/global/project/projectdirs/agora/paper_CGM/spectra_from_trident/Ion_Spectra/'
    redshifts = {'art':1.998998976095549, 'enzo':1.9999988698963, 'ramses':2.0005652511032306, 
             'gadget':2.003003004441382, 'gear':2.0000000000018687, 'gizmo':2.0012562123472377}
    redshift = redshifts[code]
    avg_num_list = []
    avg_b_list = []
    all_comp_tots = []
    for ion in ions:
        total_ion_num = 0
        total_ion_b = 0
        for sightline_num in range(sightline_range):
            f_name = f"{root}/AGORA_{code}_CR/{redshift}/Line_{sightline_num}/components.txt"
            data = load_components(f_name)
            line_comp_tot = ncomps(ion, data)
            all_comp_tots =  all_comp_tots + [line_comp_tot]
            total_ion_num = total_ion_num + line_comp_tot
            line_b_avg = bparam(ion, data)
            if line_b_avg != 0:
                total_ion_b = total_ion_b + line_b_avg
            #print("total_ion_num: ", total_ion_num)
        avg_num = total_ion_num/20
        avg_b = total_ion_b/20
        avg_num_list = avg_num_list + [avg_num]
        avg_b_list = avg_b_list + [avg_b]
    print("avg num: ", avg_num_list)
    print("avg b: ", avg_b_list)
    print("comp_tots ", all_comp_tots)
    return avg_num_list,avg_b_list,all_comp_tots