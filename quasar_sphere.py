#!/usr/bin/env python

import numpy as np
import os
import datetime
try:
    from quasarscan import parse_metadata
    from quasarscan import ion_lists
    from quasarscan import gasbinning
    from quasarscan import roman
    level = 0
except:
    import parse_metadata
    import ion_lists
    import gasbinning
    import roman
    level = 1

def get_aux_files_art(dspath):
    projectdir = dspath.split("10MpcBox")[0]
    a0 = dspath.split("a0.")[1][:3]
    file_particle_header = projectdir+"PMcrda0.%s.DAT"%a0
    file_particle_data = projectdir+"PMcrs0a0.%s.DAT"%a0
    file_particle_stars = projectdir+"stars_a0.%s.dat"%a0
    return file_particle_header,file_particle_data,file_particle_stars    

def ions_to_field_name(ions):
    lst = []
    for ion in ions:
        lst += [('gas',ion_to_field_name(ion))]
    return lst

def ion_to_field_name(ion):
    atom = ion.split(" ")[0]
    ionization = roman.from_roman(ion.split(" ")[1])-1
    return "%s_p%s_number_density"%(atom,ionization)


class GeneralizedQuasarSphere(object):
    def __init__(self, list_of_quasar_spheres, distance = "kpc"):
        if isinstance(list_of_quasar_spheres, GeneralizedQuasarSphere):
            list_of_quasar_spheres = [list_of_quasar_spheres]
        self.number = len(list_of_quasar_spheres)
        self.distance = distance
        if self.number > 0:
            self.gasbins = list_of_quasar_spheres[0].gasbins
        else:
            self.gasbins = gasbinning.GasBinsHolder(bins = None)
        
        ions_lists = []
        sum_of_lengths = 0
        
        for i in range(self.number):
            q = list_of_quasar_spheres[i]
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
                self.info[currentpos:currentpos+size,pos_in_self] = q.info[:size,pos_in_q]
                if distance == "kpc":
                    convert = q.code_unit_in_kpc
                elif distance == "Rvir":
                    convert = q.code_unit_in_kpc/q.Rvir
                else:
                    print "not sure what distance = %s means..."%distance
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
            toreturn.append("%s:fraction"%ion)
            for key in self.gasbins.get_all_keys():
                toreturn.append("%s:%s"%(ion,key))
        return toreturn

class QuasarSphere(GeneralizedQuasarSphere):
    def __init__(self,ions=None,data = None,\
                 simparams = None,scanparams = None,gasbins = None,readvalsoutput = None):
        self.number = 1
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
        print self.name
        name_fields   = self.name.split("_")
        self.simname  = name_fields[0]
        self.version  = name_fields[1]
        self.code     = name_fields[2]
        self.simnum   = name_fields[3]
        self.redshift = self.simparams[1]
        self.rounded_redshift = self.redshift
        if abs(self.redshift - 1) <= .05: self.rounded_redshift = 1.00
        elif abs(self.redshift - 1.5) <= .05: self.rounded_redshift = 1.50
        elif abs(self.redshift - 2) <= .05: self.rounded_redshift = 2.00
        elif abs(self.redshift - 3) <= .05: self.rounded_redshift = 3.00
        elif abs(self.redshift - 4) <= .05: self.rounded_redshift = 4.00
        elif abs(self.redshift - 5) <= .5: self.rounded_redshift = 5.00
        elif abs(self.redshift - 6) <= .5: self.rounded_redshift = 6.00
        elif abs(self.redshift - 8) <= 1: self.rounded_redshift = 8.00
        elif abs(self.redshift - 10) <= 2: self.rounded_redshift = 10.00
        elif abs(self.redshift - 15) <= 2: self.rounded_redshift = 15.00
        elif abs(self.redshift - 20) <= 4: self.rounded_redshift = 20.00
        else: self.rounded_redshift = self.redshift
        self.center = np.array([self.simparams[2], self.simparams[3], self.simparams[4]])
        self.Rvir = self.simparams[5]
        self.Rvir_is_real = str(parse_metadata.get_value("Rvir",self.name,redshift = self.redshift)==self.Rvir)
        self.dspath = self.simparams[6]
        self.a0 = 1./(1+self.redshift)
        self.L = np.array([self.simparams[7], self.simparams[8], self.simparams[9]])
        self.L_mag = np.sqrt(np.sum(self.L**2))

        #start looking for metadata files
        self.code_unit_in_kpc = self.simparams[10]
        self.Mvir = parse_metadata.get_value("Mvir",self.name,redshift = self.redshift)
        self.gas_Rvir = parse_metadata.get_value("Mvir",self.name,redshift = self.redshift)
        self.star_Rvir = parse_metadata.get_value("Mvir",self.name,redshift = self.redshift)
        self.dm_Rvir = parse_metadata.get_value("Mvir",self.name,redshift = self.redshift)
        self.sfr = parse_metadata.get_value("Mvir",self.name,redshift = self.redshift)
        if self.sfr and self.star_Rvir:
            self.ssfr = self.sfr / self.star_Rvir
        else:
            self.ssfr = None
        try:
            aDict = parse_metadata.avalsdict[self.simname]
            aDict.sort()
            self.final_a0 = float(aDict[-1])
        except:
            self.final_a0 = None

    def get_criteria_at_a(self, a, criteria):
        if float(self.final_a0) < float(a):
            print "Inputted a value that exceeds the greatest 'a' value in %s" %(self.simname)
        if criteria == "ssfr":
            sfr_dict = parse_vela_metadata.dict_of_vela_info("SFR")
            star_Rvir_dict = parse_vela_metadata.dict_of_vela_info("star_Rvir")
            return float(sfr_dict[self.simname][a])/float(star_Rvir_dict[self.simname][a])
        if criteria == "sfr":
            criteria = "SFR"
        criteria_dict = parse_vela_metadata.dict_of_vela_info(criteria)
        return float(criteria_dict[self.simname][a])
    #renames basic scanparams data into new instance variables
    def add_extra_scanparam_fields(self):
        self.R = self.scanparams[0]
        self.len_th_arr = self.scanparams[1]
        self.len_phi_arr = self.scanparams[2]
        self.len_r_arr = self.scanparams[3]
        self.rmax = self.scanparams[4]
        self.length = self.scanparams[5] 
        self.length_reached = self.scanparams[6]
    
    def save_values(self,dest = None,at_level = 1):
        if len(self.info[0]) <= 11:
            print "No ions!"
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
        print "saved file %s"%filename
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


