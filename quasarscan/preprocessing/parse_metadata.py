import os
import numpy as np

class NoMetadataError(Exception):
    def __init__(self, message):
        self.message = message
        print(self.message)

def dict_of_all_info(fullname):
    #This is a guide for what info can be found in what column, starting from column 0.
    try:
        simname,version,code,simnum = fullname.split('_')
    except:
        print(fullname)
        stop
  
    #Setting up the dictionary of the desired property, where the dictionary key is the name of the galaxy folder : the redshift, and the dictionary value is the respective desired property quantity
    ret_dict = {}

    float_quantities = ['a','center_x','center_y','center_z','L_x','L_y','L_z',\
        'bulk_velocity_x','bulk_velocity_y','bulk_velocity_z',\
        'Rvir','Mvir','Lmag','sfr','Mstar','Mgas','Mgas_galaxy','Mdm']

    string_quantities = ['compaction_stage']

    if os.path.isdir("quasarscan_data/galaxy_catalogs"):
        basepath = "quasarscan_data/galaxy_catalogs/%s"%(fullname)
    else: 
        raise NoMetadataError('folder "quasarscan_data/galaxy_catalogs" not found!')

    filename = os.path.join(basepath,fullname+'_metadata.txt')
    if os.path.exists(filename):
        f=open(filename)
    else:
        raise NoMetadataError('file "%s" not found!'%filename)
    lines = f.readlines()
    f.close()

    # the following lines including the while loop retrieve the desired property value by iterating through the textfile line by line 
    firstline = lines[0]
    quantities = firstline.replace('\n','').split(', ')
    column_numbers = {}
    assert firstline[0] == 'a',"Metadata Formatted Incorrectly, \nexpansion factor must be first column"
    for i,quantity in enumerate(quantities):
        if quantity not in float_quantities+string_quantities:
            print('quantity %s not recognized. Skipping.'%quantity)
            continue
        ret_dict[quantity] = {}
        column_numbers[i] = quantity

    for line in lines[1:]:
        values = line.replace('\n','').split(', ')
        a = float(values[0])
        for i,value in enumerate(values):
            if column_numbers[i] in float_quantities:
                ret_dict[column_numbers[i]][a] = float(value)
            elif column_numbers[i] in string_quantities:
                ret_dict[column_numbers[i]][a] = value
    return ret_dict
    
def dict_of_quantity(quantity,fullname,all_quantities_dict = None):
    if not all_quantities_dict:
        all_quantities_dict = dict_of_all_info(fullname)
    if quantity in all_quantities_dict.keys():
        return dict_of_all_info(fullname)[quantity]
    else:
        return None

def all_avals(fullname,all_quantities_dict = None):
    avals = np.array(list(dict_of_quantity('a',fullname,all_quantities_dict).keys()))
    avals.sort()
    return avals

def get_closest_value_for_a(fullname,a,all_quantities_dict = None):
    avals = np.array(list(dict_of_quantity('a',fullname,all_quantities_dict).keys()))
    avals.sort()
    best = -1.0
    for aval in avals:
        if np.abs(a-aval)<np.abs(best-a):
            best = aval
    if np.abs(a-best)>.2:
        return -1.0
    return best

def get_last_a_for_sim(fullname,all_quantities_dict = None):
    avals = np.array(list(dict_of_quantity('a',fullname,all_quantities_dict).keys()))
    avals.sort()
    return avals[-1]

def get_value(quantity, fullname, redshift = None,a = None,z=None, check_exists = False,loud = False):
    assert a or z or redshift, 'no timestep specified'
    if a is None:
        redshift = z if redshift is None else redshift
        a = 1./(redshift+1)
    if check_exists:
        return not np.isnan(get_value(quantity, fullname, a = a))
    simname,version,code,simnum = fullname.split("_")
    all_quantities_dict = dict_of_all_info(fullname)
    a0 = get_closest_value_for_a(fullname,a=a,all_quantities_dict=all_quantities_dict)
    one_quantity_dict = dict_of_quantity(quantity,fullname,all_quantities_dict=all_quantities_dict)
    if one_quantity_dict is None:
        if loud:
            print('quantity "%s" not in ~/quasarscan_info/galaxy_catalogs/%s_%s_%s/%s.txt'%(quantity,simname,version,code,simnum))
        return np.nan
    if a0 == -1.0:
        if loud:
            print("simulation %s does not reach expansion factor %s"%(fullname,a))
        return np.nan
    return one_quantity_dict[a0]


