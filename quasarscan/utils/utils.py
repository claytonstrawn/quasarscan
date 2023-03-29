from quasarscan.utils import roman
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import os

def round_redshift(redshift):
    if abs(redshift - 0.0) <= .05: return 0.00
    elif abs(redshift - 0.25) <= .05: return 0.25
    elif abs(redshift - 0.5) <= .05: return 0.50
    elif abs(redshift - 1) <= .05: return 1.00
    elif abs(redshift - 1.5) <= .1: return 1.50
    elif abs(redshift - 2) <= .1: return 2.00
    elif abs(redshift - 2.5) <= .1: return 2.50
    elif abs(redshift - 3) <= .1: return 3.00
    elif abs(redshift - 4) <= .1: return 4.00
    elif abs(redshift - 5) <= .5: return 5.00
    elif abs(redshift - 6) <= .5: return 6.00
    elif abs(redshift - 8) <= 1: return 8.00
    elif abs(redshift - 10) <= 2: return 10.00
    elif abs(redshift - 15) <= 2: return 15.00
    elif abs(redshift - 20) <= 4: return 20.00
    else: return float(str(redshift)[:3])

def split_by_ops(s):
    #NOTE: This will have a bug if a variable legitimately ends in e and adds/subtracts
    #the next term! this is because e-XX is how python floats are written as strings
    s = s.replace("e-","_temp_minus_char_")
    s = s.replace("e+","_temp_plus_char_")
    s = s.replace("(","_splitchar_")
    s = s.replace(")","_splitchar_")
    s = s.replace("/","_splitchar_")
    s = s.replace("*","_splitchar_")
    s = s.replace("+","_splitchar_")
    s = s.replace("-","_splitchar_")
    s = s.replace("_temp_minus_char_","e-")
    s = s.replace("_temp_plus_char_","e+")
    return list(filter(None,s.split("_splitchar_")))

def sort_ions(ions,flat = True):
    def sort_ions_one_element(ions,element):
        nums = [None]*len(ions)
        toreturn = []
        for i in range(len(ions)):
            nums[i] = roman.from_roman(ions[i].split(" ")[1])
        nums.sort()
        for val in nums:
            toreturn.append("%s %s"%(element,roman.to_roman(val)))
        return toreturn
    ions = list(ions)
    ions.sort()
    index = 0
    element = ions[index].split(" ")[0]
    tosort = []
    toreturn = []
    while index < len(ions):
        if ions[index].split(" ")[0] == element:
            tosort.append(ions[index])
            index += 1
        else:
            toreturn.append(sort_ions_one_element(tosort,element))
            element = ions[index].split(" ")[0]
            tosort = []
    toreturn.append(sort_ions_one_element(tosort,element))
    if flat:
        toreturn = [item for sublist in toreturn for item in sublist]
    return toreturn

#summary: reverse an array
#
#inputs: ary: an array
#
#outputs: the reversed array
def reversearray(ary):
    ary = list(ary)
    ary.reverse()
    try:
        return np.array(ary)
    except ValueError:
        return np.array(ary,dtype = object)

def string_represents_ion(string):
    if ' ' not in string:
        return False
    atom,ionization = string.split(' ')
    if not (len(atom)==1 or len(atom) == 2):
        return False
    try:
        roman.from_roman(ionization)
    except:
        return False
    return True

#summary: SPH codes use a "smoothing length" instead of "cell_size"
#
#inputs: code: which kind of AGORA code is under analysis here
#
#outputs: whether this is sph or amr code in question
sphcodes = ['gizmo','gadget','gear','tipsy','changa']
def get_gasbins_arg(code):
    if code in sphcodes:
        return 'all_sph'
    else:
        return 'all_amr'

def definecolorbar(bar_type = 'HotCustom',**kwargs):
    f= 256.0
    if bar_type not in ('HotCustom','RainbowCustom','BlackandWhite'):
        raise Exception("Not a ColorMap. Please try another one.")
    if bar_type == 'HotCustom':
        cdict = {'red':   ((0.0,  255/f, 255/f),
                           (1e-9, 255/f, 255/f),
                           (0.1,  255/f, 255/f),
                           (0.3,  255/f, 255/f),
                           (0.5,  139/f, 139/f),
                           (1.0,  0/f, 0/f)),

                 'green': ((0.0,  255/f, 255/f),
                           (1e-9, 255/f, 255/f),
                           (0.1,  165/f, 165/f),
                           (0.3,  0/f, 0/f),
                           (0.5,  0/f, 0/f),
                           (1.0,  0/f, 0/f)),

                 'blue':  ((0.0,  255/f, 255/f),
                           (1e-9, 255/f, 0/f),
                           (0.1,  0/f, 0/f),
                           (0.3,  0/f, 0/f),
                           (0.5,  0/f, 0/f),
                           (1.0,  0/f, 0/f))}

        custom = LinearSegmentedColormap('HotCustom', cdict)
    elif bar_type == 'RainbowCustom':
        bowdict = {'red': ((0.0, 1.0, 1.0),
                             (1e-9, 0.0, 0.0),
                             (0.3, 0.5, 0.5),
                             (0.6, 0.7, 0.7),
                             (0.9, 0.8, 0.8),
                             (1.0, 0.0, 0.0)),
                     'green': ((0.0, 1.0, 1.0),
                               (1e-9, 0.0, 0.0),
                               (0.3, 0.8, 0.8),
                               (0.6, 0.7, 0.7),
                               (0.9, 0.0, 0.0),
                               (1.0, 0.0, 0.0)),
                     'blue': ((0.0, 1.0, 1.0),
                              (1e-9, 1.0, 1.0),
                              (0.3, 1.0, 1.0),
                              (0.6, 0.0, 0.0),
                              (0.9, 0.0, 0.0), 
                              (1.0, 0.0, 0.0))}
        custom = LinearSegmentedColormap('RainbowCustom', bowdict)
    elif bar_type == 'BlackandWhite':
        bwdict = {'red': ((0.0, 1.0, 1.0),
                             (1e-9, 0.8, 0.8),
                             (0.3, 0.6, 0.6),
                             (0.6, 0.4, 0.4),
                             (0.9, 0.2, 0.2),
                             (1.0, 0.0, 0.0)),
                     'green': ((0.0, 1.0, 1.0),
                             (1e-9, 0.8, 0.8),
                             (0.3, 0.6, 0.6),
                             (0.6, 0.4, 0.4),
                             (0.9, 0.2, 0.2),
                             (1.0, 0.0, 0.0)),
                     'blue': ((0.0, 1.0, 1.0),
                             (1e-9, 0.8, 0.8),
                             (0.3, 0.6, 0.6),
                             (0.6, 0.4, 0.4),
                             (0.9, 0.2, 0.2),
                             (1.0, 0.0, 0.0))}
        custom = LinearSegmentedColormap('BlackandWhite', bwdict)
    return custom

def data_path():
    HOME = os.path.expanduser('~')
    if os.path.exists("quasarscan_data"):
        PATH = 'quasarscan_data'
    elif os.path.exists(os.path.join(HOME,"quasarscan_data")):
        PATH = os.path.join(HOME,"quasarscan_data")
    else:
        PATH = None
    return PATH

#taken from trident
# Taken from Cloudy documentation.
solar_abundance = {
    'H' : 1.00e+00, 'He': 1.00e-01, 'Li': 2.04e-09,
    'Be': 2.63e-11, 'B' : 6.17e-10, 'C' : 2.45e-04,
    'N' : 8.51e-05, 'O' : 4.90e-04, 'F' : 3.02e-08,
    'Ne': 1.00e-04, 'Na': 2.14e-06, 'Mg': 3.47e-05,
    'Al': 2.95e-06, 'Si': 3.47e-05, 'P' : 3.20e-07,
    'S' : 1.84e-05, 'Cl': 1.91e-07, 'Ar': 2.51e-06,
    'K' : 1.32e-07, 'Ca': 2.29e-06, 'Sc': 1.48e-09,
    'Ti': 1.05e-07, 'V' : 1.00e-08, 'Cr': 4.68e-07,
    'Mn': 2.88e-07, 'Fe': 2.82e-05, 'Co': 8.32e-08,
    'Ni': 1.78e-06, 'Cu': 1.62e-08, 'Zn': 3.98e-08}

atomic_mass = {
    'H' : 1.00794,   'He': 4.002602,  'Li': 6.941,
    'Be': 9.012182,  'B' : 10.811,    'C' : 12.0107,
    'N' : 14.0067,   'O' : 15.9994,   'F' : 18.9984032,
    'Ne': 20.1797,   'Na': 22.989770, 'Mg': 24.3050,
    'Al': 26.981538, 'Si': 28.0855,   'P' : 30.973761,
    'S' : 32.065,    'Cl': 35.453,    'Ar': 39.948,
    'K' : 39.0983,   'Ca': 40.078,    'Sc': 44.955910,
    'Ti': 47.867,    'V' : 50.9415,   'Cr': 51.9961,
    'Mn': 54.938049, 'Fe': 55.845,    'Co': 58.933200,
    'Ni': 58.6934,   'Cu': 63.546,    'Zn': 65.409}

from yt.utilities.physical_constants import mh

def agora_custom_metals(field,data):
    if isinstance(field.name, tuple):
        ftype = field.name[0]
        field_name = field.name[1]
    else:
        ftype = "gas"
        field_name = field.name
    atom = field_name.split("_")[1]
    
    H_mass_fraction = 0.76
    to_nH = H_mass_fraction / mh
    return data['gas', "agora_metallicity"]*\
            data['gas', 'density']*\
            data.ds.quan(solar_abundance[atom],'1/Zsun')*\
            to_nH
