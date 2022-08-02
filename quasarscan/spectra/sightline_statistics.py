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

def nparam(ion, comp_list):
    total_n = 0
    line_n_avg = 0
    counter = 0
    for comp in comp_list:
        for line in comp.list_of_lines:
            if ion == line.ion:
                counter += 1
                total_n = total_n + np.log10(line.N)
    if counter != 0:
        line_n_avg = total_n / counter
    return line_n_avg

def sightline_looping(code,sightline_range,ions,include_zero = False):
    root = '/global/project/projectdirs/agora/paper_CGM/spectra_from_trident/Ion_Spectra/'
    redshifts = {'art':1.998998976095549, 'enzo':1.9999988698963, 'ramses':2.0005652511032306, 
             'gadget':2.003003004441382, 'gear':2.0000000000018687, 'gizmo':2.0012562123472377}
    redshift = redshifts[code]
    avg_num_list = []
    avg_b_list = []
    avg_n_list = []
    comp_tots = {}
    line_bs = {}
    line_ns = {}
    for ion in ions:
        total_ion_num = 0
        total_ion_b = 0
        total_ion_n = 0
        all_comp_tots = []
        all_bs = []
        all_ns = []
        for sightline_num in range(sightline_range):
            f_name = f"{root}/AGORA_{code}_CR/{redshift}/Line_{sightline_num}/components.txt"
            data = load_components(f_name)
            line_comp_tot = ncomps(ion, data)
            if include_zero == True:
                all_comp_tots =  all_comp_tots + [line_comp_tot]
            elif line_comp_tot != 0:
                all_comp_tots =  all_comp_tots + [line_comp_tot]
            line_b_avg = bparam(ion, data)
            if line_b_avg != 0:
                all_bs =  all_bs + [line_b_avg]
            line_n_avg = nparam(ion, data)
            if line_n_avg != 0:
                all_ns =  all_ns + [line_n_avg]
                
        avg_num = np.average(all_comp_tots)
        avg_b = np.average(all_bs)
        avg_n = np.average(all_ns)
        avg_num_list = avg_num_list + [avg_num]
        avg_b_list = avg_b_list + [avg_b]
        avg_n_list = avg_n_list + [avg_n]
        comp_tots[ion] = all_comp_tots
        line_bs[ion] = all_bs
        line_ns[ion] = all_ns
    return avg_num_list,avg_b_list,avg_n_list,comp_tots,line_bs,line_ns