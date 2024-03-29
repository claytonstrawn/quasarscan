import numpy as np
from quasarscan.data_objects.quasar_sphere import QuasarSphere
from quasarscan.utils.utils import string_represents_ion,round_redshift
from quasarscan.data_objects.gasbinning import GasBinsHolder
from quasarscan.utils.variable_lists import param_xVars
from quasarscan.spectra.trident_interfacer import convert_equivalent_width_to_coldens

def get_header_columns_dict(header):
    to_return = {}
    for i,name in enumerate(header.split(',')):
        if name.startswith('log '):
            islog = True
            retname = name.split('log ')[1]
        else:
            islog = False
            retname = name
        to_return[retname] = i,islog
    return to_return

def access_obs_data_values(variable,header_dict,parsed_line):
    strvars = ['sl','cgm','sightline']
    numvars = ['z','Rvir','Mstar','sfr','sfr:eb','r']
    if variable in strvars:
        try:
            i,islog = header_dict[variable]
            assert not islog
            return parsed_line[i]
        except Exception as e:
            return None
    elif variable in numvars or 'cdens' in variable or 'eb' in variable:
        try:
            i,islog = header_dict[variable]
            if parsed_line[i] == '<':
                value = -2.0
            elif parsed_line[i] == '>':
                value = -3.0
            elif parsed_line[i] == 'nan':
                value = np.nan
            else:
                value = float(parsed_line[i])
            if islog and value > 0.0:
                return 10.**float(value)
            else:
                return value
        except Exception as e:
            return np.nan
    else:
        print("%s should not be stored in the sightline!"%variable)
        assert False

#input: filename to go into
#output: a line of headers representing columns, and the lines of the file
def read_obs_textfile(filename):
    with open(filename) as f:
        fullname_nonum = f.readline()[:-1]
        f.readline()
        f.readline()
        header = f.readline()[:-1]
        f.readline()
        lines = f.read().splitlines()
        lines = list(filter(None, lines))
    return fullname_nonum,header, lines

class ObsQuasarSphere(QuasarSphere):
    def __init__(self, fullname, header, line_from_file):
        self.type = 'obs'
        self.code = 'obs'
        self.number = 1
        self.length = 1
        self.length_reached = 1
        self.gasbins = GasBinsHolder(bins = None)
        parsed_line = line_from_file.split(",") 
        header_columns_dict = get_header_columns_dict(header)
        temp = header_columns_dict.keys()
        ions = []
        for key in temp:
            if ':' in key:
                ion,modifier = key.split(':')
                if string_represents_ion(ion) and ion not in ions:
                    ions.append(ion)
        self.ions = ions
        self.intensives = []
        iondata = []
        for ion in self.ions:
            if f'{ion}:cdens' in temp:
                iondata.append(access_obs_data_values('%s:cdens'%ion,header_columns_dict,parsed_line))
                iondata.append(access_obs_data_values('%s:eb'%ion,header_columns_dict,parsed_line))
            elif f'{ion}:EW' in temp:
                ew = access_obs_data_values('%s:EW'%ion,header_columns_dict,parsed_line)
                ew_eb = access_obs_data_values('%s:EW_eb'%ion,header_columns_dict,parsed_line)
                cdens,cdens_eb = convert_equivalent_width_to_coldens(ew,ew_eb,ion)
                iondata.append(cdens)
                iondata.append(cdens_eb)
        self.info = np.array([[-1, -1, -1]+
                             [access_obs_data_values('r',header_columns_dict,parsed_line)]+
                             [-1, -1, -1,-1,-1, -1, -1] + 
                             iondata + 
                             [-1, -1, -1]])

        redshift = access_obs_data_values('z',header_columns_dict,parsed_line)
        self.Rvir = access_obs_data_values('Rvir',header_columns_dict,parsed_line)
        self.a = 1./(redshift+1)
        self.Mvir = 4./3*np.pi*self.Rvir**3*(200*1.5e-7) #crit_dens = 1.5e-7 M⊙ pc-3 
        self.Mstar = access_obs_data_values('Mstar',header_columns_dict,parsed_line)
        self.Mstar_eb = access_obs_data_values('Mstar:eb',header_columns_dict,parsed_line)
        self.sfr = access_obs_data_values('sfr',header_columns_dict,parsed_line)
        self.sfr_eb = access_obs_data_values('sfr:eb',header_columns_dict,parsed_line)
        self.ssfr = self.sfr/self.Mstar
        self.ssfr_eb = self.sfr_eb/self.Mstar
        super().__init__(fullname,redshift)
        self.author = self.version
        self.survey = self.simname
        
        for param in param_xVars:
            if param not in self.__dict__.keys():
                self.__dict__[param] = np.nan
        
