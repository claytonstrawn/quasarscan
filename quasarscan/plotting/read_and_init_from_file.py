from quasarscan.data_objects import observation_quasar_sphere,\
                                     quasar_sphere,\
                                     simulation_quasar_sphere
from quasarscan.preprocessing import parse_metadata
from quasarscan.utils.utils import data_path
import os
import numpy as np

class NoFilesError(Exception):
    def __init__(self):
        pass

#presumes you pass it a string, or a list or tuple of strings
#will return all files/folders whose names contain that string
def one_is_in_name(name,loadonly):
    if loadonly == 'all':
        return True
    elif loadonly == 'none':
        return False
    #if a string is passed, 
    if isinstance(loadonly,str):
        loadonly = [loadonly]        
    for fragment in loadonly:
        if fragment in name:
            return True
    return False

#summary: search directory for textfiles
#
#inputs: inquasarscan: if True, only look down one level. If false, look two.
#        loadonly: if not 'all' only load certain simulations (e.g. 'VELA')
#
#outputs: textfiles: list of names of textfiles
def get_all_textfiles(loadonly,qtype,throw_errors = False):
    PATH = data_path()
    assert os.path.exists(PATH), "folder 'quasarscan_data' not found, so nothing available to plot!"
    if qtype == 'sim':
        path = os.path.join(PATH,"output")
        if not os.path.exists(path):
            print("folder 'quasarscan_data/output' not found, so no simulations available to plot!")
            if throw_errors:
                raise NoFilesError()
        dirs = os.listdir(path)
    elif qtype == 'obs':
        path = os.path.join(PATH,"observations")
        if not os.path.exists(path):
            print("folder 'quasarscan_data/observations' not found, so no observations available to plot!")
            if throw_errors:
                raise NoFilesError()
    elif qtype == 'empty':
        path = os.path.join(PATH,"galaxy_catalogs")
        if not os.path.exists(path):
            print("folder 'quasarscan_data/galaxy_catalogs' not found, so no metadata available to plot!")
            if throw_errors:
                raise NoFilesError()
    dirs = os.listdir(path)
    textfiles = []
    #gets all folders in output
    for folder_name in dirs:
        if (not folder_name.startswith(".")) and one_is_in_name(folder_name,loadonly):
            folder_path = os.path.join(path,folder_name)
            folder_dirs = os.listdir(folder_path)
            for file_name in folder_dirs:
                if not file_name.startswith("."):
                    textfiles.append(os.path.join(folder_path,file_name))
    return textfiles


def read_and_init(loadonly,qtype):
    quasar_array = []
    textfiles = get_all_textfiles(loadonly,qtype)
    if qtype == 'sim':
        for file in textfiles:
            readvalsoutput = simulation_quasar_sphere.read_values(file)
            q = simulation_quasar_sphere.SimQuasarSphere(start_up_info_packet = readvalsoutput)
            quasar_array.append(q)
    elif qtype == 'obs':
        for file in textfiles:
            fullname_nonum,header,lines = observation_quasar_sphere.read_obs_textfile(file)
            for i,line in enumerate(lines):
                fullname = fullname_nonum+'_'+str(i+1)
                oq = observation_quasar_sphere.ObsQuasarSphere(fullname, header, line)
                quasar_array.append(oq)
    elif qtype == 'empty':
        for file_name in textfiles:
            fullname = file_name.split('/')[2]
            all_avals = parse_metadata.all_avals(fullname)
            for a in all_avals:
                eq = quasar_sphere.EmptyQuasarSphere(fullname, redshift=1/a - 1)
                quasar_array.append(eq)
    return np.array(quasar_array)


