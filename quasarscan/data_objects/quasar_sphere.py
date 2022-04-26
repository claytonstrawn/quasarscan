import numpy as np
from quasarscan.preprocessing.parse_metadata import get_value,get_last_a_for_sim,dict_of_all_info
from quasarscan.utils.utils import string_represents_ion,round_redshift
from quasarscan.data_objects.gasbinning import GasBinsHolder

class NoTypeError(Exception):
    def __init__(self, message):
        self.message = message
        print(self.message)

class IonNotFoundError(Exception):
    def __init__(self, message):
        self.message = message
        print(self.message)

class QuasarSphere(object):
    def __init__(self, fullname, redshift):
        try: 
            self.type
        except AttributeError:
            raise NoTypeError('Cannot instantiate QuasarSphere superclass')
        self.fullname = fullname
        self.redshift = redshift
        self.add_id_variables()
        self.add_metadata_dict()
        self.lookup_fields_by_name()
        self.create_ion_list_w_bins()

    def add_id_variables(self):
        if 'multi' in self.type:
            return
        name_fields = self.fullname.split('_')
        # if reading an observation, the fullname will be 
        # survey (COS-Halos, CASBAH, etc), author for cite, 'obs', sightline #
        self.simname  = name_fields[0]
        self.version  = name_fields[1]
        self.code     = name_fields[2]
        self.simnum   = name_fields[3]
        self.rounded_redshift = round_redshift(self.redshift)
        self.a = 1./(1+self.redshift)
        
    def add_metadata_dict(self):
        if 'multi' in self.type or 'obs' in self.type:
            return
        self.metadata_dict = dict_of_all_info(self.fullname)

    def lookup_fields_by_name(self):
        if 'multi' in self.type or 'obs' in self.type:
            return
        #look for this info in metadata files
        #any not found will be nan
        for key in self.metadata_dict.keys():
            if key[-2:] == '_x':
                x = get_value(key,self.fullname,redshift = self.redshift,\
                              all_quantities_dict = self.metadata_dict)
                y = get_value(key.replace('_x','_y'),self.fullname,redshift = self.redshift,\
                              all_quantities_dict = self.metadata_dict)
                z = get_value(key.replace('_x','_z'),self.fullname,redshift = self.redshift,\
                              all_quantities_dict = self.metadata_dict)
                self.__setattr__(key.replace('_x',''),np.array([x,y,z]))
            elif key[-2:] in ['_y','_z']:
                continue
            else:
                v = get_value(key,self.fullname,redshift = self.redshift,\
                              all_quantities_dict = self.metadata_dict)
                self.__setattr__(key,v)
        self.final_a0 = get_last_a_for_sim(self.fullname,all_quantities_dict = self.metadata_dict)

    def create_ion_list_w_bins(self):
        if 'empty' in self.type:
            self.ion_list_w_bins = []
            return
        toreturn = ['i', 'theta', 'phi', 'r', 'alpha', \
                    'startx', 'starty', 'startz', 'endx', 'endy', 'endz']
        for ion in self.ions:
            toreturn.append('%s:cdens'%ion)
            if 'sim' in self.type:
                toreturn.append('%s:fraction'%ion)
                for key in self.gasbins.get_all_keys():
                    toreturn.append('%s:%s'%(ion,key))
            elif 'obs' in self.type:
                toreturn.append('%s:eb'%ion)
        if 'sim' in self.type:
            for intensive in self.intensives:
                toreturn.append(intensive)
        self.ion_list_w_bins = toreturn

    def get_ion_column_num(self,ion):
        if ':' not in ion and string_represents_ion(ion):
            ion = ion+':cdens'
        if ion in self.ion_list_w_bins:
            return self.ion_list_w_bins.index(ion)
        else:
            print(self.ion_list_w_bins[-10:])
            raise IonNotFoundError('Ion %s not found in this Sphere, which is %s. Available ions are: %s and sorted into bins: %s. It could also be a metadata error. Known quantities are %s'%(ion,self.fullname,self.ions,self.gasbins.get_all_keys(),self.__dict__.keys()))

    def get_ion_values(self,ion):
        if ion == 'zeros':
            return self.info[:,0]*0.0
        elif isinstance(ion,list):
            return np.array([self.get_ion_values(i) for i in ion])
        column_num = self.get_ion_column_num(ion)
        return self.info[:,column_num]

    def get_sightline_values(self,xVar):
        vardict = {"theta":1,"phi":2,"r":3}
        column_num = vardict[xVar]
        return self.info[:,column_num]

class EmptyQuasarSphere(QuasarSphere):
    def __init__(self,fullname,redshift):
        self.type = 'empty'
        self.number = 1
        self.length = 0
        self.ions = []
        self.info = []
        self.gasbins = GasBinsHolder(bins = None)
        super().__init__(fullname,redshift)

