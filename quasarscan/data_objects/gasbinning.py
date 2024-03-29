import numpy as np
from quasarscan import mode
if 'write' in mode:
    from pi_or_ci import PI_field_defs

class BadBinNameError(Exception):    
    def __init__(self, message):
        self.message = message
        print(self.message)

#read a dictionary of a string and return a dictionary object
def parse_string_for_dict(s):
    s = s.strip('{}')
    toreturn = {}
    list_of_terms = s.split('),')
    for word in list_of_terms:
        key,value = word.split(':')
        value = value.strip('()')
        value0,value1 = value.split(',')
        toreturn[key] = (value0,value1)
    return toreturn

# GasBin object contains a field, units, and edges for several 
# bins to put data in
class GasBin(object):
    #summary: creates GasBin object
    #
    #inputs: name: how the gasbin is referred to
    #        binnames: what the bins are "called" (cold, cool, PI, etc)
    #        binvalstr: what the bin edges are (numerically), in 'units'
    #        field: the yt field you will check against the binval. Default 
    #               is 'name'
    #        units: the units of the field. Defaults to whatever yt outputs as default.
    #
    #outputs: GasBin object
    def __init__(self,name,binnames,binvalstr,field = None,units = None):
        self.name = name
        self.units = units
        if not field:
            self.field = ('gas',name)
        elif field == "NotImplemented":
            print("cannot instantiate field %s!"%name)
            return
        else:
            self.field = field
        assert len(binnames) == len(binvalstr)-1
        self.binvalstr = binvalstr
        self.binvals = [None]*len(binvalstr)
        for i,item in enumerate(self.binvalstr):
            self.binvals[i] = eval(item)
        #TODO - convert bin names to be stored in all caps
        for binname in binnames:
            if binname[-1] == 'e':
                raise BadBinNameError('Bin names cannot end in "e" (it looks too much like Python float syntax)')
        self.binnames = binnames

    def __eq__(self,other):
        if self.name == other.name:
            if self.binvalstr == other.binvalstr and self.binnames == other.binnames and self.units == other.units:
                return True
            else:
                return False
        return False

    def __ne__(self,other):
        return not self.__eq__(other)

    def get_length(self):
        return len(self.binnames)

    def get_keys(self):
        newlist = []
        for val in self.binnames:
            assert val[-1] != 'e', "variable names cannot end with 'e', it looks too much like scientific notation in a formula"
            newlist.append(self.name+":"+val)
        return newlist

if 'write' in mode:
    # below are all the known gasbins. More could easily be made
    density_bin = GasBin("density",["low","between6_4","between4_2","between2_0","selfshielding","starforming"],["0.0",'1.674e-30','1.674e-28','1.674e-26',"1.6737e-25","1.6737e-23","np.inf"],units = 'g/cm**3')#[0,1e-6,1e-4,1e-2,1e-1,1e1,inf]*mh
    temperature_bin = GasBin("temperature",["cold","cool","warm_hot","hot"],["0.0","10**3.8","10**4.5","10**6.5","np.inf"],units = "K")
    radial_velocity_bin = GasBin("radial_velocity",["inflow",'inflow_slow',"tangential",'outflow_slow',"outflow"],['-np.inf','-500','-10','10','500','np.inf'],units = "km/s")
    resolution_bin_amr = GasBin("resolution",["high","between1_5","between5_15","low"],['0.00e+00', '1e9', '1.250e+11','3.375e+12', 'np.inf'],field = ('gas','cell_volume'),units = "pc**3")# (['0.00e+00', '1000', '5000','15000', 'np.inf']pc)**3
    resolution_bin_sph = GasBin("resolution",["high","between1_5","between5_15","low"],['0.00e+00', '1000', '5000','15000', 'np.inf'],field = ('gas','smoothing_length'),units = "pc")# (['0.00e+00', '1000', '5000','15000', 'np.inf']pc)**3
    all_ions_w_known_PI_defs = PI_field_defs.make_funcs()[0]
    allionizationfields = {}
    for ion in all_ions_w_known_PI_defs:
        allionizationfields[ion]=('gas','PI_%s'%(ion.replace(' ','')))
    pi_bin = GasBin('ionization_mechanism',['PI'],['0.9','1.1'], field = allionizationfields)
    possible_bin_types = ["density","temperature","radial_velocity","resolution","ionization_mechanism"]

def ds_has_field(ds,gb):
    if ds is None:
        return True
    if isinstance(gb.field,str):
        return ('gas',gb.field) in ds.derived_field_list
    elif isinstance(gb.field,tuple):
        return gb.field in ds.derived_field_list
    elif isinstance(gb.field,dict):
        key = list(gb.field.keys())[0]
        return gb.field[key] in ds.derived_field_list
    else:
        return False


# GasBinsHolder object contains multiple GasBin objects
# and can query them
class GasBinsHolder(object):
    def __init__(self,ds=None,bins = None,string = None):
        self.bin_types = []
        if string:
            fields_with_data = string.strip("[]").split(", ")
            field_names,bins,edges,ytfields,units,unique_names = [],[],[],[],[],[]
            for item in fields_with_data:
                split_by_colon = item.split(":")
                field_name,current_bin,current_edges = split_by_colon[0:3]
                field_names.append(field_name)
                unique_names.append(field_name) if field_name not in unique_names else None
                bins.append(current_bin)
                edges.append(current_edges)
                if len(split_by_colon)>3 and "field" in split_by_colon[3]:
                    if split_by_colon[3].split("-")[1][0] == '{':
                        ytfields.append(parse_string_for_dict(item.split('{')[1].split('}')[0]))
                    else:
                        fieldlist = split_by_colon[3].split("-")[1].strip("()").split(",")
                        ytfields.append((fieldlist[0],fieldlist[1]))
                else:
                    ytfields.append(None)
                if len(split_by_colon) > 3 and "units" in split_by_colon[-1]:
                    units.append(split_by_colon[-1].split("-")[1])
                else:
                    units.append(None)
            for field_name in unique_names:
                current_bins = []
                current_edges = []
                current_ytfield = None
                current_units = None
                for i in range(len(field_names)):
                    if field_name == field_names[i]:
                        current_bins.append(bins[i])
                        begin,end = edges[i].split('_')
                        current_edges.append(begin)
                        current_ytfield = ytfields[i]
                        current_units = units[i]
                newBin = GasBin(field_name,current_bins,current_edges+[end],
                                field = current_ytfield,units = current_units)
                self.bin_types.append(newBin)
            return
        if bins == "all_amr":
            bins = possible_bin_types
            res_bin = resolution_bin_amr
        elif bins == "all_sph":
            bins = possible_bin_types
            res_bin = resolution_bin_sph
        elif bins is None:
            bins = []
        if "density" in bins:
            assert ds_has_field(ds,density_bin), f'{density_bin.field} field not found'
            self.bin_types.append(density_bin)
        if "temperature" in bins:
            assert ds_has_field(ds,temperature_bin), f'{temperature_bin.field} field not found'
            self.bin_types.append(temperature_bin)
        if "radial_velocity" in bins:
            assert ds_has_field(ds,radial_velocity_bin), f'{radial_velocity_bin.field} field not found'
            self.bin_types.append(radial_velocity_bin)
        if "resolution" in bins:
            assert ds_has_field(ds,res_bin), f'{res_bin.field} field not found'
            self.bin_types.append(res_bin)
        if "ionization_mechanism" in bins:
            assert ds_has_field(ds,pi_bin), f'{pi_bin.field} field not found'
            self.bin_types.append(pi_bin)

    def __contains__(self, key):
        for b in self.bin_types:
            if key == b.name:
                return True
        return False

    def combine_holders(self,other):
        if len(self.bin_types) == 0:
            self.bin_types = other.bin_types
        elif len(other.bin_types) == 0:
            pass
        else:
            bins = []
            for obj in self.bin_types:
                if obj in other.bin_types:
                    bins.append(obj)
            self.bin_types = bins

    def get_all_keys(self):
        mylist = []
        for binvar in self.bin_types:
            mylist += binvar.get_keys()
        return mylist

    def get_field_binedges_for_key(self,key,ion):
        var_name,bin_name = key.split(":")
        for gb in self.bin_types:
            if gb.name == var_name:
                if isinstance(gb.field, dict):
                    try:
                        fld = gb.field[ion]
                    except KeyError:
                        fld = None
                elif isinstance(gb.field,tuple):
                    fld = gb.field
                i = gb.binnames.index(bin_name)
                return fld,(gb.binvals[i],gb.binvals[i+1]),gb.units
        print("that field does not exist!")
        assert 0 == 1

    def get_field_binedges_for_num(self,num,ion):
        mykeys = self.get_all_keys()
        return self.get_field_binedges_for_key(mykeys[num],ion)

    def get_all_bintypes(self):
        return self.bin_types

    def get_length(self):
        mysum = 0
        for binvar in self.bin_types:
            mysum += binvar.get_length()
        return mysum

    def get_bin_str(self,numbers = True,append = ""):
        if self.bin_types == []:
            return "No bins, GasBinsHolder is empty"
        to_return = "filler"
        current_to_add = append
        for gb in self.bin_types:
            current_to_add+=gb.name
            for i,val in enumerate(gb.binnames):
                if numbers:
                    specific_val=":%s:%s_%s"%(val,gb.binvalstr[i],gb.binvalstr[i+1])
                else:
                    specific_val=":%s"%val
                if gb.field != ('gas',gb.name):
                    specific_val += ":field-%s"%(str(gb.field).replace(', ',',').replace(': ',':').replace("'",""))
                if gb.units:
                    specific_val += ":units-%s"%gb.units
                to_return+=", %s%s"%(current_to_add,specific_val)
            current_to_add = append
        return to_return.replace("filler, ","")
