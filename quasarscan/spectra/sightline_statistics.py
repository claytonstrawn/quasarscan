import numpy as np
from quasarscan.spectra.data_objects import load_components
import math

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

#summary: returns total number of components for an ion in a sightline as an array 
#
#inputs: ion: ion to count
#        comp_list: list of components for the sightline
#
#output: array of the total number of components for an ion in a sightline
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

#summary: returns total b for components of the ion in a line as an array 
#
#inputs: ion: ion to count
#        comp_list: list of components for the sightline
#
#output: array of the total b for components of the ion in a line as an array 
def bparam(ion, comp_list):
    all_bs = []
    for comp in comp_list:
        for line in comp.list_of_lines:
            if ion == line.ion:
                all_bs.append(line.b)
                break
    return all_bs

#summary: returns total column density components of the ion in a line as an array 
#
#inputs: ion: ion to count
#        comp_list: list of components for the sightline
#
#output: array of the total column density for an ion of a component in a sightline
def nparam(ion, comp_list):
    all_ns = []
    for comp in comp_list:
        for line in comp.list_of_lines:
            if ion == line.ion:
                all_ns.append(np.log10(line.N))
                break
    return all_ns

def covering_fraction(ion,comp_list):
    for comp in comp_list:
        for line in comp.list_of_lines:
            if ion == line.ion:
                return [1]
    return [0]

def remove_outliers_nans(code,redshift,ion,var,all_var):
    all_var_array = np.array(all_var)
    all_var_array = all_var_array[~np.isnan(all_var_array)]
    if code == 'gear' and redshift == 3.00000000000869 and var == 'b':
        #gear has very few detections (<5) at z=3 but 
        #one of them has extremely wide (>5x as much as all others)
        #components, giving a huge spread
        print('gear',all_var_array,'removing b outliers...')
        all_var_array = all_var_array[all_var_array<100]
    return list(all_var_array)

#summary: returns average and standard error of chosen variables(total components, b, column density, covering fraction) for each ion chosen for the code chosen
#
#inputs: code: galaxy code
#       sightline_range: sightlines used to find average
#       ions: list of ions to find averages for
#       vars_to_return: list of what to find the averages for, can choose from total components, b, column density, or covering fraction
#       include_zero: choose to include even if there's 0 components for an ion, default is false
#       min_data_points: minimum number of components needed to be included
#
#output: array of averages for each chosen variable for each chosen ion and an array of standard errors for each chosen variable for each chosen ion
def sightline_looping(code,redshift,sightline_range,ions,
                      vars_to_return = [],include_zero = False, 
                      min_data_points = 2,root = 'Ion_Spectra_new',
                      errors = 'warn',snr = None):
    root = '/global/cfs/projectdirs/agora/paper_CGM/spectra_from_trident/'+root
    vars_funcs_dict = {'num':ncomps,'b':bparam, 
                       'n':nparam,'covering_fraction':covering_fraction}
    #add a second thing, in addition to 
    #to_return values, return the errors on each
    to_return = np.zeros((len(vars_to_return),len(ions)))+np.nan        
    return_SE = np.zeros((len(vars_to_return),len(ions)))+np.nan 
    for i,var in enumerate(vars_to_return):
        try:
            func = vars_funcs_dict[var]
        except KeyError as e:
            print(f'key {var} is not a known function. '+\
                  f'Known functions are {list(vars_funcs_dict.keys())}')
            raise e
        for j,ion in enumerate(ions):
            all_var = []
            for sightline_num in range(sightline_range):
                add_snr = '' if snr is None else f"_snr{snr}"
                f_name = f"{root}/AGORA_{code}_CR/{redshift}/"+\
                        f"Line_{sightline_num}/components{add_snr}.txt"
                try:
                    data = load_components(f_name)
                except FileNotFoundError as e:
                    if errors == 'warn':
                        print(e)
                    elif errors is True:
                        raise(e)
                    continue
                one_ion_results = func(ion,data)
                if var in ['num']:
                    if include_zero:
                        all_var =  all_var + one_ion_results
                    else:
                        if one_ion_results[0] != 0:
                            all_var =  all_var + one_ion_results
                else:
                    all_var =  all_var + one_ion_results
            all_var = remove_outliers_nans(code,redshift,ion,var,all_var)
            if var != 'num' and len(all_var) < min_data_points:
                avg_var = np.nan
            else:
                avg_var = np.average(all_var)
            sd_var = np.std(all_var, ddof = 1)
            se_var = sd_var/math.sqrt((len(all_var)))
            return_SE[i][j] = se_var             
            to_return[i][j] = avg_var
                
    return to_return, return_SE