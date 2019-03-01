import parse_vela_metadata
import parse_nihao_metadata
import numpy as np

functions = {'VELA':parse_vela_metadata.dict_of_vela_info,'NIHAO':parse_nihao_metadata.dict_of_nihao_info}
avalsdict = {'VELA':parse_vela_metadata.adict,'NIHAO':parse_nihao_metadata.adict}


def get_closest_value_for_a(redshift,simname,name):
    try:
        a0 = 1./(redshift+1)
        avals = avalsdict[simname]
        best = '-1.0'
        for a in avals.keys()[name]:
            if np.abs(a0-float(a))<np.abs(float(best)-a0):
                best = a
        if a0-float(a)>.2:
            return None
        return best
    except:
        return None

def get_value(quantity, name, redshift):
    simname = name.split("_")[0]
    a0 = get_closest_value_for_a(redshift,simname,name)
    if a0 is None:
        print "simulation %s does not go down to redshift %s"%(name,redshift)
    try: 
        return float(functions[simname](quantity)[name][a0])
    except:
        return None
