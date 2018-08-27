import numpy as np
from yt import YTArray

class GasBinsHolder(Object):
    self.bins_dict = {}
    self.bin_types = []
    self.possible_bin_types = ["density","temperature","radial_velocity","resolution"]

    def __init__(self,bins):
        if bins == "all":
            bins = self.possible_bin_types
        elif bins is None:
            bins = []
        if "density" in bins:
            self.add_density()
        if "temperature" in bins:
            self.add_temperature()
        if "radial_velocity" in bins:
            self.add_radial_velocity()
        if "resolution" in bins:
            self.add_resolution()

    def combine_holders(self,others):
        bins = self.bin_types
        for other in others:
            bins = list(set(bins).union(set(other.bin_types)))
        if "density" in bins:
            self.add_density()
        if "temperature" in bins:
            self.add_temperature()
        if "radial_velocity" in bins:
            self.add_radial_velocity()
        if "resolution" in bins:
            self.add_resolution()

    def get_bins_for_key(self,key):
        try: 
            if not self.bins_dict[key] is None:
                return self.bins_dict[key]
        except:
            print("key '%s' is not recognized. Add it to gasbinning.py?"%key)
        raise KeyError

    def get_all_keys(self):
        mykeys = self.bins_dict.keys()
        mykeys.sort()
        return mykeys

    def get_all_bintypes(self):
        mytypes = self.bin_types
        return mytypes

    def get_length(self):
        mykeys = self.get_all_keys()
        return len(mykeys)

    def add_density(self):
        if "density" in self.bin_types:
            return
        self.bin_types += ["density"]
        self.bins_dict["density:low"] = YTArray([1e-4,1e-2],"g/cm**3")
        self.bins_dict["density:high"] = YTArray([1e-2,1e+0],"g/cm**3")

    def add_temperature(self):
        if "temperature" in self.bin_types:
            return
        self.bin_types += ["temperature"]
        self.bins_dict["temperature:cold"] = YTArray([0.0,10.0**3.8],"K")
        self.bins_dict["temperature:cool"] = YTArray([10.0**3.8,10.0**4.5],"K")
        self.bins_dict["temperature:warm_hot"] = YTArray([10.0**4.5,10.0**6.5],"K")

    def add_radial_velocity(self):
        if "radial_velocity" in self.bin_types:
            return
        self.bin_types += ["radial_velocity"]      
        self.bins_dict["radial_velocity:inflow"] = YTArray([-np.inf,-10],"km/s")
        self.bins_dict["radial_velocity:tangential"] = YTArray([-10,10],"km/s")
        self.bins_dict["radial_velocity:outflow"] = YTArray([10,np.inf],"km/s")

    def add_resolution(self):
        if "resolution" in self.bin_types:
            return
        print "add_resolution is not yet implemented" 
