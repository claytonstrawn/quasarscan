#!/usr/bin/env python

import numpy as np
import os
import datetime
from functools import reduce

try:
    from quasarscan import parse_metadata
    from quasarscan import ion_lists
    from quasarscan import gasbinning
    from quasarscan import roman
    #from quasarscan.observational_quasar_sphere import Observation
    level = 0
except:
    import parse_metadata
    import ion_lists
    import gasbinning
    import roman
    #from observational_quasar_sphere import Observation
    level = 1
 

def ions_to_field_name(ions):
    lst = []
    for ion in ions:
        lst += [('gas',ion_to_field_name(ion))]
    return lst

def ion_to_field_name(ion,field_type = "number_density"):
    atom = ion.split(" ")[0]
    ionization = roman.from_roman(ion.split(" ")[1])-1
    return "%s_p%s_%s"%(atom,ionization,field_type)

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

class GeneralizedQuasarSphere(object):
    def __init__(self, list_of_quasar_spheres, distance = "kpc"):
        if isinstance(list_of_quasar_spheres, GeneralizedQuasarSphere):
            list_of_quasar_spheres = [list_of_quasar_spheres]
        self.number = len(list_of_quasar_spheres)
        self.distance = distance
        self.gasbins = gasbinning.GasBinsHolder(bins = None)
        
        ions_lists = []
        sum_of_lengths = 0
        
        for i in range(self.number):
            q = list_of_quasar_spheres[i]
            if i == 0:
                self.type = q.type
            else:
                assert self.type == q.type, "quasar_sphere types must be the same."
            ions_lists.append(q.ions) 
            sum_of_lengths += q.length
            self.gasbins.combine_holders(q.gasbins)
        if self.number > 0:
            ions_in_all = list(reduce(set.intersection, map(set, ions_lists)))
        else:
            ions_in_all = []
        self.ions = ions_in_all
        self.length = sum_of_lengths
        num_extra_columns = self.gasbins.get_length()
        self.info = np.zeros((self.length,11+len(self.ions)*(num_extra_columns+2)+3))
        currentpos = 0
        for i in range(self.number):
            q = list_of_quasar_spheres[i]
            size = min(q.length_reached,q.length)
            self.info[currentpos:currentpos+size,:11] = q.info[:size,:11]
            for ion in self.get_ion_list_w_bins():
                pos_in_q = q.get_ion_column_num(ion)
                pos_in_self = self.get_ion_column_num(ion)
                self.debug = q.info
                self.info[currentpos:currentpos+size,pos_in_self] = q.info[:size,pos_in_q]
                convert = 1.0
                if distance == "Rvir":
                    convert/=q.Rvir
            self.info[currentpos:currentpos+size,-1] = q.info[:size,-1]
            self.info[currentpos:currentpos+size,-2] = q.info[:size,-2] 
            self.info[currentpos:currentpos+size,-3] = q.info[:size,-3] 
            self.info[currentpos:currentpos+size,3]*=convert
            currentpos += size
                   
    def get_ion_column_num(self,ion):
        intensivesdict = {'Z':-1,'rho':-2,'T':-3}
        if not ":" in ion:
            bintype = "cdens"
        else:
            i = ion.index(":")
            bintype = ion[i+1:]
            ion = ion[:i]
        if bintype == "cdens":
            plus = 0
        elif bintype == "fraction":
            plus = 1
        elif bintype == "eb":
            plus = 1
        else:
            plus = self.gasbins.get_all_keys().index(bintype)+2
        if ion in intensivesdict.keys():
            return len(self.info[0])+intensivesdict[ion]
        try:
            return 11 + self.ions.index(ion)*(self.gasbins.get_length()+2)+plus
        except:
            print("Ion %s not found in this Sphere. Please try restricting to ions = ['%s']."%(ion,ion))

    def get_ion_list_w_bins(self):
        
        toreturn = []
        for ion in self.ions:
            toreturn.append("%s:cdens"%ion)
            if self.type == "Simulation":
                toreturn.append("%s:fraction"%ion)
            if self.type == "Observation":
                toreturn.append("%s:eb"%ion)
            for key in self.gasbins.get_all_keys():
                toreturn.append("%s:%s"%(ion,key))
        return toreturn

class QuasarSphere(GeneralizedQuasarSphere):
    def __init__(self,ions=None,data = None,\
                 simparams = None,scanparams = None,gasbins = None,readvalsoutput = None):
        self.number = 1
        self.type = "Simulation"
        if readvalsoutput:
            simparams  = readvalsoutput[0]
            scanparams = readvalsoutput[1]
            ions       = readvalsoutput[2]
            data       = readvalsoutput[3]
            gasbins    = readvalsoutput[4]
        self.simparams = simparams
        self.scanparams = scanparams
        self.add_extra_scanparam_fields()
        if type(ions) is list:
            self.ions = ions
        elif type(ions) is str:
            self.ions = ions[1:-1].split(", ")
        else:
            self.ions = []
        self.info = data
        self.gasbins = gasbins
        self.add_extra_simparam_fields()

   
    #renames basic simparams data into new instance variables,
    #finds and saves stellar mass, total mass (virial mass), and the star formation rate
    def add_extra_simparam_fields(self):
        self.name     = self.simparams[0]
        name_fields   = self.name.split("_")
        self.simname  = name_fields[0]
        self.version  = name_fields[1]
        self.code     = name_fields[2]
        self.simnum   = name_fields[3]
        self.redshift = self.simparams[1]
        self.rounded_redshift = round_redshift(self.redshift)
        self.center = np.array([self.simparams[2], self.simparams[3], self.simparams[4]])
        self.Rvir = self.simparams[5]
        self.Rvir_is_real = str(parse_metadata.get_value("Rvir",self.name,redshift = self.redshift)==self.Rvir)
        self.dspath = self.simparams[6]
        self.a0 = 1./(1+self.redshift)
        self.L = np.array([self.simparams[7], self.simparams[8], self.simparams[9]])
        self.L_mag = np.sqrt(np.sum(self.L**2))
        self.conversion_factor = self.simparams[10]
        self.code_unit = self.simparams[11]
        #start looking for metadata files
        self.Mvir = parse_metadata.get_value("Mvir",self.name,redshift = self.redshift)
        self.Mvir_cgm = self.Mvir
        self.gas_Rvir = parse_metadata.get_value("gas_Rvir",self.name,redshift = self.redshift)
        self.Mgas = self.gas_Rvir
        self.Mgas_cgm = self.gas_Rvir - parse_metadata.get_value("gas_.1Rvir",self.name,redshift = self.redshift)
        self.star_Rvir = parse_metadata.get_value("star_Rvir",self.name,redshift = self.redshift)
        self.Mstar = self.star_Rvir
        self.Mstar_cgm = self.star_Rvir - parse_metadata.get_value("star_.1Rvir",self.name,redshift = self.redshift)
        self.dm_Rvir = parse_metadata.get_value("dm_Rvir",self.name,redshift = self.redshift)
        self.Mdm = self.dm_Rvir
        self.sfr = parse_metadata.get_value("SFR",self.name,redshift = self.redshift)
        if self.sfr and self.star_Rvir:
            self.ssfr = self.sfr / self.star_Rvir
        else:
            self.ssfr = None
        compaction_start = parse_metadata.get_value("compaction_start",self.name)
        compaction_end = parse_metadata.get_value("compaction_end",self.name)
        if ~np.isnan(compaction_start):
            if self.redshift > compaction_start:
                self.compaction_stage = "pre"
            elif self.redshift < compaction_end:
                self.compaction_stage = "post"
            else:
                self.compaction_stage = "during"
        else:
            self.compaction_stage = "unknown"
        try:
            aDict = parse_metadata.avalsdict[self.simname][self.name]
            aDict_list = sorted(aDict.keys())
            self.final_a0 = aDict_list[-1]
        except:
            self.final_a0 = None

    def get_criteria_at_a(self, a0, criteria):
        if self.final_a0 < a0:
            print("Inputted a value that exceeds the greatest 'a' value in %s" %(self.simname))
        if criteria == "ssfr":
            sfr = parse_metadata.get_value("SFR",self.name,a0=a0)
            star_Rvir = parse_metadata.get_value("star_Rvir",self.name,a0=a0)
            return sfr/star_Rvir
        if criteria == "sfr":
            criteria = "SFR"
        return parse_metadata.get_value(criteria,self.name, a0=a0)
    #renames basic scanparams data into new instance variables
    def add_extra_scanparam_fields(self):
        self.R = self.scanparams[0]
        self.len_th_arr = self.scanparams[1]
        self.len_phi_arr = self.scanparams[2]
        self.len_r_arr = self.scanparams[3]
        self.rmax = self.scanparams[4]
        self.length = self.scanparams[5] 
        self.length_reached = self.scanparams[6]

    def save_values(self,dest = None,at_level = 1,test = False):
        if len(self.info[0]) <= 11:
            print("No ions!")
        linesfinished = self.length_reached
        numlines = self.length
        redshift = self.rounded_redshift
        simname = self.name
        ionsstr = ion_lists.filename_stringform(self.ions)
        if ionsstr in ion_lists.dict_of_ionstrsfilename.keys():
            ionsstr = ion_lists.dict_of_ionstrsfilename[ionsstr]
        if dest:
            filename = dest
        else:
            if level == 0 or at_level == 0:
                levelvar = "quasarscan/"
            elif level == 1:
                levelvar = ""
            foldername = levelvar+"output/"+simname+"coldensinfo"
            if not os.path.exists(foldername):
                os.makedirs(foldername)
            specificfilename = "%s_of_%s-"%(str(linesfinished),str(numlines)) +ionsstr+"_z"+str(redshift)[:4]+".txt"
            filename = foldername+"/"+specificfilename
            if test:
                print(filename)
                return
            prev = os.listdir(foldername)
            for item in prev:
                if item.endswith("of_%s-"%str(numlines) +ionsstr+"_z"+str(redshift)[:4]+".txt"):
                    os.remove(foldername+"/"+item)
        f = open(filename,"w+")
        firstfew = [None]*7
        firstfew[0] = "Simparams: [dsname, z, center[0], center[1], center[2], Rvir, pathname, L[0], L[1], L[2], convert]\n"
        firstfew[1] = "Simparams: "+str(self.simparams)+"\n"
        firstfew[2] = "Scanparams: [R, n_th, n_phi, n_r, r_max, num_lines, line_reached]\n"
        firstfew[3] = "Scanparams: "+str(self.scanparams)+"\n"
        firstfew[4] = "Ions: %s\n"%str(self.ions)
        firstfew[5] = "Bins: [%s]\n"%self.gasbins.get_bin_str()
        firstfew[6] = "Data: [i, theta, phi, r, alpha, startx, starty, startz, endx, endy, endz, "
        for ion in self.ions:
            toAdd = "%s:cdens, %s:fraction, %s, "%(ion,ion,self.gasbins.get_bin_str(numbers = False,append = ion+":"))
            firstfew[6] += toAdd
        firstfew[6] += "T, n, Z]\n"
        f.write(firstfew[0])
        f.write(firstfew[1])
        f.write(firstfew[2])
        f.write(firstfew[3])
        f.write(firstfew[4])
        f.write(firstfew[5])
        f.write(firstfew[6])
        for vector in self.info:
            f.write(str(vector).replace("\n",""))
            f.write("\n")
        f.close()
        print("saved file %s"%filename)
        return filename
    
def read_values(filename):
    f = open(filename)
    firstfew = [None]*7
    firstfew[0] = f.readline()
    firstfew[1] = f.readline()[:-1].split(": ")[1]
    firstfew[2] = f.readline()
    firstfew[3] = f.readline()[:-1].split(": ")[1]
    firstfew[4] = f.readline()[:-1].split(": ")[1]
    firstfew[5] = f.readline()[:-1].split(": ")[1]
    firstfew[6] = f.readline()[:-1].split(": ")[1]
    simparams = eval(firstfew[1])
    scanparams = eval(firstfew[3])
    ions = eval(firstfew[4])
    gasbins = gasbinning.GasBinsHolder(string = firstfew[5])
    length = scanparams[5]
    num_extra_columns = gasbins.get_length()
    data = np.zeros((int(length),11+len(ions)*(num_extra_columns+2)+3))
    for i in range(length):
        myline = f.readline()[1:-1]
        data[i] = np.fromstring(myline,sep = " ")
    return simparams,scanparams,ions,data,gasbins


