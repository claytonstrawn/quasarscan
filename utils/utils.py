from quasarscan.utils import roman
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

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
    return np.array(ary)

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

#summary: SPH codes do not have a simple "resolution" so for now we need to 
#         not use that information
#
#inputs: code: which kind of AGORA code is under analysis here
#
#outputs: 'noresolution' if code is SPH, 'all' if code is AMR
sphcodes = ['gizmo','gadget','gear','tipsy']
def get_gasbins_arg(code):
    if code in sphcodes:
        return 'noresolution'
    else:
        return 'all'

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