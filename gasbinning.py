import numpy as np
from yt import YTArray

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
            print "using two different definitions of %s!"%name
            return False
        return False

    def __ne__(self,other):
        return not self.__eq__(other)

    def get_length(self):
        return len(self.binnames)

    def get_keys(self):
        newlist = []
        for val in self.binnames:
            newlist.append(self.name+":"+val)
        return newlist

density_bin = GasBin("density",["low","medium","selfshielding","starforming"],["0.0","1e-2","1e-1","1e1","np.inf"])
temperature_bin = GasBin("temperature",["cold","cool","warm-hot","hot"],["0.0","10**3.8","10**4.5","10**6.5","np.inf"])
radial_velocity_bin = GasBin("radial_velocity",["inflow","tangential","outflow"],['-np.inf','-10','10','np.inf'],units = "km/s")
resolution_bin = GasBin("resolution",["high","medium","low"],['0','200','1000','np.inf'],field = ('gas','dx'),units = "pc")


class GasBinsHolder(object):

    def __init__(self,bins = None,string = None):
        assert not (bins is None and string is None)
        self.bin_types = []
        if string:
            allnames = string.strip("[]").split(", ")
            i = 0 
            while i < len(allnames):
                currentBin = allnames[i].split(":")[0]
                allinfo = allnames[i].split(":")
                currentBinNames = [allinfo[1]]
                currentBinVals = [allinfo[2].split("_")[0],allinfo[2].split("_")[1]]
                j = 1
                field = None
                while i+j < len(allnames) and allnames[i+j].startswith(currentBin):
                    allinfo = allnames[i+j].split(":")
                    currentBinNames += [allinfo[1]]
                    currentBinVals += [allinfo[2].split("_")[1]]
                    if len(allinfo) >= 4:
                        if "field" in allinfo[3]:
                            fieldlist = allinfo[3].split("-")[1].strip("()").split(",")
                            field = (fieldlist[0],fieldlist[1])
                        if len(allinfo) > 4 and "units" in allinfo[4]:
                            units = allinfo[4].split("-")[1]
                        elif "units" in allinfo[3]:
                            units = allinfo[3].split("-")[1]
                    else:
                        field = None
                    j+=1
                i+=j
                newBin = GasBin(currentBin,currentBinNames,currentBinVals,field = field)
                self.bin_types.append(newBin)
            return

        self.possible_bin_types = ["density","temperature","radial_velocity","resolution"]
        if bins == "all":
            bins = self.possible_bin_types
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

    def combine_holders(self,other):
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

    def get_field_binedges_for_key(self,key):
        var_name,bin_name = key.split(":")
        for gb in self.bin_types:
            if gb.name == var_name:
                i = gb.binnames.index(bin_name)
                return gb.field,(gb.binvals[i],gb.binvals[i+1]),gb.units
        print "that field does not exist!"
        assert 0 == 1

    def get_field_binedges_for_num(self,num):
        mykeys = self.get_all_keys()
        return self.get_field_binedges_for_key(mykeys[num])

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

    def add_temperature(self):
        if "temperature" in self.bin_types:
            return
        self.bin_types += ["temperature"]
        self.bins_dict["temperature:cold"] = YTArray([0.0,10.0**3.8],"K")
        self.bins_dict["temperature:cool"] = YTArray([10.0**3.8,10.0**4.5],"K")

    def add_resolution(self):
        if "resolution" in self.bin_types:
            return
        print "add_resolution is not yet implemented" 