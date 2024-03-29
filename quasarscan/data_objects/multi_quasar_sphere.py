import numpy as np
from quasarscan.data_objects.quasar_sphere import QuasarSphere
from quasarscan.data_objects.gasbinning import GasBinsHolder
from functools import reduce

class NoRvirError(Exception):
    def __init__(self, message):
        self.message = message
        print(self.message)

class MultiQuasarSphere(QuasarSphere):
    def __init__(self, list_of_quasar_spheres, distance = "kpc"):
        self.number = len(list_of_quasar_spheres)
        self.items = list_of_quasar_spheres
        self.distance = distance
        self.gasbins = GasBinsHolder(bins = None)
        self.fullname = []
        mytype = 'unknown'
        
        ions_lists = []
        intensives_lists = []
        sum_of_lengths = 0
        
        for i in range(self.number):
            q = list_of_quasar_spheres[i]
            self.fullname.append(q.fullname)
            if mytype == 'unknown':
                mytype = q.type
            else:
                assert mytype == q.type, "QuasarSphere types must be the same."
            ions_lists.append(q.ions) 
            intensives_lists.append(q.intensives) 
            sum_of_lengths += q.length_reached
            self.gasbins.combine_holders(q.gasbins)
        self.type = 'multi_'+mytype
        if self.number > 0:
            ions_in_all = list(reduce(set.intersection, map(set, ions_lists)))
            intensives_in_all = list(reduce(set.intersection, map(set, intensives_lists)))
        else:
            ions_in_all = []
            intensives_in_all = []
        self.ions = ions_in_all
        self.intensives = intensives_in_all
        self.length = sum_of_lengths
        num_extra_columns = self.gasbins.get_length()
        self.info = np.zeros((self.length,11+len(self.ions)*(num_extra_columns+2)+len(self.intensives)))
        currentpos = 0
        for i in range(self.number):
            q = list_of_quasar_spheres[i]
            size = q.length_reached
            self.info[currentpos:currentpos+size,:11] = q.info[:size,:11]
            self.create_ion_list_w_bins()
            for ion in self.ion_list_w_bins:
                pos_in_q = q.get_ion_column_num(ion)
                pos_in_self = self.get_ion_column_num(ion)
                self.info[currentpos:currentpos+size,pos_in_self] = q.info[:size,pos_in_q]
            convert = 1.0
            if distance == "Rvir":
                if not np.isnan(q.Rvir):
                    convert/=q.Rvir
                else:
                    raise NoRvirError('Rvir is not known for all simulations, cannot combine this way!')
            for i in range(len(self.intensives)):
                self.info[currentpos:currentpos+size,-i] = q.info[:size,-i]
            self.info[currentpos:currentpos+size,3]*=convert
            currentpos += size