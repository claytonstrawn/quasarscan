import parse_vela_metadata
import parse_nihao_metadata
import numpy as np

functions = {'VELA':parse_vela_metadata.dict_of_vela_info,'NIHAO':parse_nihao_metadata.dict_of_nihao_info}
avalsdict = {'VELA':parse_vela_metadata.adict,'NIHAO':parse_nihao_metadata.adict}

#the problem is that NIHAO is not organized nice round expansion factors (a)
#this will probably be more common, so I will use it as the default behavior
#so I need this function to find the best a value in a given system, starting from redshift 
#(which is what people are more likely to be interested in asking for)
def get_closest_value_for_a(redshift,simname,name):
    if simname not in functions.keys():
        return None
    a0 = 1./(redshift+1)
    avals = avalsdict[simname]
    best = '-1.0'
    for a in avals[name].keys():
        if np.abs(a0-float(a))<np.abs(float(best)-a0):
            best = a
    if a0-float(a)>.2:
        return None
    return float(best)

def get_value(quantity, name, redshift = None,a = None):
    simname = name.split("_")[0]
    if a:
        a0 = a
    elif redshift:
        a0 = get_closest_value_for_a(redshift,simname,name)
    else:
        print "called without specifying a time"
        raise ValueError
    if a0 is None:
        print "simulation %s does not reach redshift %s"%(name,redshift)
    try: 
        return float(functions[simname](quantity)[name][a0])
    except:
        return None
