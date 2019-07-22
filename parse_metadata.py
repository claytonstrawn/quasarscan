try:
    from quasarscan import parse_vela_metadata
    from quasarscan import parse_nihao_metadata
    from quasarscan import parse_agora_metadata
    level = 0
except: 
    import parse_vela_metadata
    import parse_nihao_metadata
    import parse_agora_metadata
import numpy as np
import os

functions = {'VELA':parse_vela_metadata.dict_of_vela_info,'NIHAO':parse_nihao_metadata.dict_of_nihao_info,'AGORA':parse_agora_metadata.dict_of_agora_info}
available_quantities = {'VELA':list(parse_vela_metadata.quantity_dict.keys()),'NIHAO':list(parse_nihao_metadata.quantity_dict.keys()),'AGORA':list(parse_agora_metadata.quantity_dict.keys())}
avalsdict = {'VELA':parse_vela_metadata.adict,'NIHAO':parse_nihao_metadata.adict,'AGORA':parse_agora_metadata.adict}

#the problem is that NIHAO is not organized nice round expansion factors (a)
#this will probably be more common, so I will use it as the default behavior
#so I need this function to find the best a value in a given system, starting from redshift 
#(which is what people are more likely to be interested in asking for)

def get_closest_value_for_a(simname,name,redshift=None,a0=None,loud = False):
    if simname not in functions.keys():
        if loud:
            print("that simulation %s does not have metadata in parse_metadata yet"%simname)
        return None
    if not redshift is None:
        a0 = 1./(redshift+1)
    elif a0:
        a0 = a0
    else:
        print("called without specifying a time")
        return None
    avals = avalsdict[simname]
    best = -1.0
    for a in avals[name].keys():
        if np.abs(a0-a)<np.abs(best-a0):
            best = a
    if np.abs(a0-best)>.2:
        return -1.0
    return best

def get_compaction_times(quantity, name):
    try:
        if os.path.isdir('galaxy_catalogs'):
            f = open('galaxy_catalogs/vela_compaction_times.txt')
        else:
            f = open('quasarscan/galaxy_catalogs/vela_compaction_times.txt')
    except:
        print('Couldn\'t open file.')
        return np.nan
    line = f.readline().split(' ')
    simname = name.split('_')[0] + name.split('_')[3]
    if name.split('_')[1] != 'v2' or name.split('_')[2] != 'art':
        return np.nan
    while line[0] != simname:
        line = f.readline()
        if line == '':
            #print('Reached end of file. Couldn\'t find desired simulation.')
            f.close()
            return np.nan
        line = line.split(' ')
    f.close()
    if quantity == 'compaction_start':
        return float(line[1])
    elif quantity == 'compaction_end':
        return float(line[2])

def get_value(quantity, name, redshift = None,a0 = None, check_exists = False):
    if check_exists:
        return not np.isnan(get_value(quantity, name, redshift = redshift,a0 = a0))
    if quantity == 'compaction_start' or quantity == 'compaction_end':
        return get_compaction_times(quantity, name)
    simname = name.split("_")[0]
    a0 = get_closest_value_for_a(simname,name,redshift=redshift,a0=a0)
    if a0 == -1.0:
        print("simulation %s does not reach redshift %s"%(name,redshift))
        return np.nan
    if simname in functions.keys() and quantity in available_quantities[simname]:
        return functions[simname](quantity)[name][a0]
    else:
        return np.nan
