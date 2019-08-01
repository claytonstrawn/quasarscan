import numpy as np

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

class GasBin(object):

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
        self.binnames = binnames

    def __eq__(self,other):
        if self.name == other.name:
            if self.field == other.field and self.binvalstr == other.binvalstr\
             and self.binnames == other.binnames and self.units == other.units:
                return True
            print("using two different definitions of %s!"%self.name)
            return False
        return False

    def __ne__(self,other):
        return not self.__eq__(other)

    def get_length(self):
        return len(self.binnames)

    def get_keys(self):
        newlist = []
        for val in self.binnames:
            #can't have that, it looks too much like scientific notation
            assert val[-1] != 'e'
            newlist.append(self.name+":"+val)
        return newlist

density_bin = GasBin("density",["low","between6_4","between4_2","between2_0","selfshielding","starforming"],["0.0",'1.674e-30','1.674e-28','1.674e-26',"1.6737e-25","1.6737e-23","np.inf"],units = 'g/cm**3')#[0,1e-6,1e-4,1e-2,1e-1,1e1,inf]*mh
temperature_bin = GasBin("temperature",["cold","cool","warm_hot","hot"],["0.0","10**3.8","10**4.5","10**6.5","np.inf"],units = "K")
radial_velocity_bin = GasBin("radial_velocity",["inflow",'inflow_slow',"tangential",'outflow_slow',"outflow"],['-np.inf','-500','-10','10','500','np.inf'],units = "km/s")
resolution_bin = GasBin("resolution",["high","between1_5","between5_15","low"],['0.00e+00', '1e9', '1.250e+11','3.375e+12', 'np.inf'],field = ('gas','cell_volume'),units = "pc**3")# (['0.00e+00', '1000', '5000','15000', 'np.inf']pc)**3
pi_bin = GasBin('ionization_mechanism',['PI'],['0.9','1.1'], field = {'O IV':('gas','PI_OIV'), 'O V':('gas','PI_OV'), 'O VI':('gas','PI_OVI'), 'O VII':('gas','PI_OVII'), 'O VIII':('gas','PI_OVIII')})
#todo: I have to think about this one a little more
possible_bin_types = ["density","temperature","radial_velocity","resolution"]


class GasBinsHolder(object):

    def __init__(self,bins = None,string = None):
        self.bin_types = []
        if string:
            fields_with_data = string.strip("[]").split(", ")
            field_names,bins,edges,fields,units = [],[],[],[],[]
            for item in fields_with_data:
                split_by_colon = item.split(":")
                field_name,current_bin,current_edges = split_by_colon[0:3]
                field_names.append(field_name)
                bins.append(current_bin)
                edges.append(current_edges)
                if len(split_by_colon)>3 and "field" in split_by_colon[3]:
                    if split_by_colon[3].split("-")[1][0] == '{':
                        fields.append(parse_string_for_dict(item.split('{')[1].split('}')[0]))
                    else:
                        fieldlist = split_by_colon[3].split("-")[1].strip("()").split(",")
                        fields.append((fieldlist[0],fieldlist[1]))
                else:
                    fields.append(None)
                if len(split_by_colon) > 3 and "units" in split_by_colon[-1]:
                    units.append(split_by_colon[-1].split("-")[1])
                else:
                    units.append(None)
            for field_name in np.unique(field_names):
                current_bins = []
                current_edges = []
                current_field = None
                current_units = None
                for i in range(len(field_names)):
                    if field_name == field_names[i]:
                        current_bins.append(bins[i])
                        begin,end = edges[i].split('_')
                        current_edges.append(begin)
                        current_field = fields[i]
                        current_units = units[i]
                newBin = GasBin(field_name,current_bins,current_edges+[end],field = current_field,units = current_units)
                self.bin_types.append(newBin)
            return
        if bins == "all":
            bins = possible_bin_types
        elif bins == "noresolution":
            bins = list(possible_bin_types)
            bins.remove('resolution')
        elif bins is None:
            bins = []
        if "density" in bins:
            self.bin_types.append(density_bin)
        if "temperature" in bins:
            self.bin_types.append(temperature_bin)
        if "radial_velocity" in bins:
            self.bin_types.append(radial_velocity_bin)
        if "resolution" in bins:
            self.bin_types.append(resolution_bin)
            #not sure how this works in demeshed
        if "ionization_mechanism" in bins:
            self.bin_types.append(pi_bin)

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
                    specific_val += ":field-%s"%str(gb.field).replace(" ","").replace("'","")
                if gb.units:
                    specific_val += ":units-%s"%gb.units
                to_return+=", %s%s"%(current_to_add,specific_val)
            current_to_add = append
        return to_return.replace("filler, ","")
