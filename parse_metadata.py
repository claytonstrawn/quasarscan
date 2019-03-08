import parse_vela_metadata
import parse_nihao_metadata
import numpy as np

functions = {'VELA':parse_vela_metadata.dict_of_vela_info,'NIHAO':parse_nihao_metadata.dict_of_nihao_info}
avalsdict = {'VELA':parse_vela_metadata.adict,'NIHAO':parse_nihao_metadata.adict}

#the problem is that NIHAO is not organized nice round expansion factors (a)
#this will probably be more common, so I will use it as the default behavior
#so I need this function to find the best a value in a given system, starting from redshift 
#(which is what people are more likely to be interested in asking for)

def get_closest_value_for_a(simname,name,redshift=None,a0=None,loud = False):
    if simname not in functions.keys():
        if loud:
            print "that simulation %s does not have metadata in parse_metadata yet"%simname
        return None
    if redshift:
        a0 = 1./(redshift+1)
    elif a0:
        a0 = a0
    else:
        print "called without specifying a time"
        return None
    avals = avalsdict[simname]
    best = -1.0
    for a in avals[name].keys():
        if np.abs(a0-a)<np.abs(best-a0):
            best = a
    if np.abs(a0-best)>.2:
        return -1.0
    return best

def get_value(quantity, name, redshift = None,a = None):
    simname = name.split("_")[0]
    a0 = get_closest_value_for_a(simname,name,redshift=redshift,a0=a)
    if a0 == -1.0:
        print "simulation %s does not reach redshift %s"%(name,redshift)
        return np.nan
    if simname in functions.keys():
        return functions[simname](quantity)[name][a0]
    else:
        return np.nan
