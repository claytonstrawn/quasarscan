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
    return [line_comp_tot]

#if any lines have the ion
#collect b for each component within that line
#but only counting ones where this ion is found
#return the average b
def bparam(ion, comp_list):
    all_bs = []
    for comp in comp_list:
        for line in comp.list_of_lines:
            if ion == line.ion:
                all_bs.append(line.b)
    return all_bs

def covering_fraction(ion,comp_list):
    for comp in comp_list:
        for line in comp.list_of_lines:
            if ion == line.ion:
                return [1]
    return [0]

def nparam(ion, comp_list):
    all_ns = []
    for comp in comp_list:
        for line in comp.list_of_lines:
            if ion == line.ion:
                all_ns.append(np.log10(line.N))
    return all_ns

def sightline_looping(code,sightline_range,ions,vars_to_return = [],include_zero = False):
    root = '/global/project/projectdirs/agora/paper_CGM/spectra_from_trident/Ion_Spectra/'
    redshifts = {'art':1.998998976095549, 'enzo':1.9999988698963, 'ramses':2.0005652511032306, 
             'gadget':2.003003004441382, 'gear':2.0000000000018687, 'gizmo':2.0012562123472377}
    redshift = redshifts[code]
    vars_funcs_dict = {'num':ncomps,'b':bparam,'n':nparam,'covering_fraction':covering_fraction}
    to_return = np.zeros((len(vars_to_return),len(ions)))+np.nan
    for i,var in enumerate(vars_to_return):
        try:
            func = vars_funcs_dict[var]
        except KeyError as e:
            print(f'key {var} is not a known function. Known functions are {list(vars_funcs_dict.keys())}')
            raise e
        for j,ion in enumerate(ions):
            all_var = []
            for sightline_num in range(sightline_range):
                f_name = f"{root}/AGORA_{code}_CR/{redshift}/Line_{sightline_num}/components.txt"
                data = load_components(f_name)
                one_ion_results = func(ion,data)
                if var in ['num']:
                    if include_zero:
                        all_var =  all_var + one_ion_results
                    else:
                        if one_ion_results[0] != 0:
                            all_var =  all_var + one_ion_results
                else:
                    all_var =  all_var + one_ion_results
            avg_var = np.average(all_var)
            to_return[i][j] = avg_var
    return to_return