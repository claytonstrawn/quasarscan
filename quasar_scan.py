#!/usr/bin/env python

import numpy as np
import trident
import yt
import os
import sys
import itertools
import logging

from mpi4py import MPI

try:
    from quasarscan import parse_vela_metadata
except:
    import parse_vela_metadata

yt.funcs.mylog.setLevel(50)

def convert_to_xyz(r, theta, phi):
    return np.array([r*np.sin(theta)*np.cos(phi),r*np.sin(theta)*np.sin(phi),r*np.cos(theta)])


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2.0)
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def get_rotation_matrix(L):
    zhat = np.array([0,0,1])
    if np.array_equal(L,zhat):
        return np.diag([1,1,1])
    theta = np.arccos(L[2]/np.linalg.norm(L))
    axis = np.cross(zhat,L)
    return rotation_matrix(axis,theta)
        

def ray_endpoints_spherical(R,r,theta,phi,alpha,endonsph):
    start = convert_to_xyz(R,theta,phi)
    xhat = convert_to_xyz(1,np.pi/2,np.pi/2+phi)
    yhat = convert_to_xyz(1,np.pi/2-theta,np.pi+phi)
    mid = r*(np.cos(alpha)*xhat+np.sin(alpha)*yhat)
    diff = start-mid
    if endonsph:
        t = 2*np.dot(start,diff)/np.dot(diff,diff)
    else:
        t = 2*R/np.linalg.norm(diff)
    end = start*(1-t)+mid*t
    return np.array([start,end])

def weights(array,function):
    if function == "sin":
        probs = np.sin(array)/2
        probs[0] = probs[-1]
    elif function == "lin":
        probs = np.linspace(0,1,len(array)+1)[1:]
    probs /= np.sum(probs)
    return probs

def ions_to_field_name(ions):
    lst = []
    for ion in ions:
        lst += [('gas',ion_to_field_name(ion))]
    return lst

def ion_to_field_name(ion):
    atom = ion.split(" ")[0]
    ionization = trident.roman.from_roman(ion.split(" ")[1])-1
    return "%s_p%s_number_density"%(atom,ionization)

class GeneralizedQuasarSphere(object):
    def __init__(self, list_of_quasar_spheres, name, distance = "kpc"):
        if isinstance(list_of_quasar_spheres, GeneralizedQuasarSphere):
            list_of_quasar_spheres = [list_of_quasar_spheres]
        self.number = len(list_of_quasar_spheres)
        self.distance = distance
        self.simname = name
        
        ions_lists = []
        sum_of_lengths = 0
        
        for i in range(self.number):
            q = list_of_quasar_spheres[i]
            ions_lists.append(q.ions) 
            sum_of_lengths += q.length
        if self.number > 0:
            ions_in_all = list(reduce(set.intersection, map(set, ions_lists)))
        else:
            ions_in_all = []
        self.ions = ions_in_all
        self.length = sum_of_lengths
        self.info = np.zeros((self.length,11+len(self.ions)+1))
        self.simname_arr = []
        self.redshift_arr = []
        self.center_arr = []
        self.Rvir_arr = []
        self.a0_arr = []
        self.L_arr = []
        self.conversion_arr = []
        self.Mvir_arr = []
        self.gas_Rvir_arr = []
        self.star_Rvir_arr = []
        self.dm_Rvir_arr = []
        self.sfr_arr = []

        
        currentpos = 0
        for i in range(self.number):
            q = list_of_quasar_spheres[i]
            size = min(q.length_reached,q.length)
            self.info[currentpos:currentpos+size,:11] = q.info[:size,:11]
            for ion in self.ions:
                pos_in_q = q.get_ion_column_num(ion)
                pos_in_self = self.get_ion_column_num(ion)
                self.info[currentpos:currentpos+size,pos_in_self] = q.info[:size,pos_in_q]
                if distance == "kpc":
                    convert = q.code_unit_in_kpc
                elif distance == "Rvir":
                    convert = q.code_unit_in_kpc/q.Rvir
                else:
                    print("not sure what distance = %s means..."%distance)
            self.info[currentpos:currentpos+size,3]*=convert
            currentpos += size
                
            self.simname_arr+=q.simname_arr
            self.redshift_arr+=q.redshift_arr
            self.center_arr+=q.center_arr
            self.Rvir_arr+=q.Rvir_arr
            self.a0_arr+=q.a0_arr
            self.L_arr+=q.L_arr
            self.conversion_arr+=q.conversion_arr
            self.Mvir_arr+=q.Mvir_arr
            self.gas_Rvir_arr+=q.gas_Rvir_arr
            self.star_Rvir_arr+=q.star_Rvir_arr
            self.dm_Rvir_arr+=q.dm_Rvir_arr
            self.sfr_arr+=q.sfr_arr
            
    def get_ion_column_num(self,ion):
        return 11 + self.ions.index(ion)

class QuasarSphere(GeneralizedQuasarSphere):
    def __init__(self,ions=None,simname=None,dspath=None,data = None,\
                 simparams = None,scanparams = None,Rvir = None,L=np.array([0,0,1]),\
                 ytlevel = "quiet",readonly = False):
        self.number = 1
        if ytlevel == "loud":
            yt.funcs.mylog.setLevel(1)
        elif ytlevel == "medium":
            yt.funcs.mylog.setLevel(10)
        if simparams == None:
            #need to load simulation from filename
            if dspath:
                projectdir = dspath.split("10MpcBox")[0]
                a0 = dspath.split("a0.")[1][:3]
                self.a0 = "0." + a0
                file_particle_header = projectdir+"PMcrda0.%s.DAT"%a0
                file_particle_data = projectdir+"PMcrs0a0.%s.DAT"%a0
                file_particle_stars = projectdir+"stars_a0.%s.dat"%a0
                self.ds = yt.load(dspath,file_particle_header=file_particle_header,\
                                  file_particle_data=file_particle_data,\
                                  file_particle_stars=file_particle_stars)
                z = self.ds.current_redshift
                #cstrings = parse_vela_metadata.dict_of_vela_info("cm")[simname][self.a0]
                #c = np.array([float(cstrings[0]),float(cstrings[1]),float(cstrings[2])])
                c = self.ds.find_max("density")[1]

                convert = self.ds.length_unit.in_units('kpc').value.item()
            else:
                #for testing without loading real sim
                self.ds = None
                z = -1.0
                c = np.zeros(3)
            self.simparams = [None]*11
            self.simparams[0]  = simname
            self.simparams[1]  = z
            self.simparams[2]  = c[0]
            self.simparams[3]  = c[1]
            self.simparams[4]  = c[2]
            self.simparams[5]  = Rvir
            self.simparams[6]  = dspath
            self.simparams[7]  = L[0]
            self.simparams[8]  = L[1]
            self.simparams[9]  = L[2]
            self.simparams[10] = convert
            self.scanparams = None
        else:
            self.simparams = simparams
            self.scanparams = scanparams
            self.add_extra_scanparam_fields()
            if not readonly:
                self.ds = yt.load(simparams[6])
        if type(ions) is list:
            self.ions = ions
        elif type(ions) is str:
            self.ions = ions[1:-1].split(", ")
        else:
            self.ions = []
        self.info = data
        self.add_extra_simparam_fields()
        
    def get_ion_column_num(self,ion):
        return 11 + self.ions.index(ion)
   
    #renames basic simparams data into new instance variables,
    #finds and saves stellar mass, total mass (virial mass), and the star formation rate
    def add_extra_simparam_fields(self):
        self.simname = self.simparams[0]
        self.simnum = self.simname[-2:]
        self.simname_arr = [self.simname]
        self.simnum_arr = [self.simnum]
        self.version = "v1.0"
        if self.simname.find("v2") != -1:
            self.version = "v2.0"
        elif self.simname.find("v3") != -1:
            self.version = "v3.0"
        elif self.simname.find("v3.1") != -1:
            self.version = "v3.1"
        self.version_arr = [self.version]
        self.redshift = self.simparams[1]
        self.redshift_arr = [self.redshift]
        self.rounded_redshift = self.redshift
        if abs(self.redshift - 1) <= .05: self.rounded_redshift = 1.00
        if abs(self.redshift - 1.5) <= .05: self.rounded_redshift = 1.50
        if abs(self.redshift - 2) <= .05: self.rounded_redshift = 2.00
        if abs(self.redshift - 3) <= .05: self.rounded_redshift = 3.00
        if abs(self.redshift - 4) <= .05: self.rounded_redshift = 4.00
        self.rounded_redshift_arr = [self.rounded_redshift]
        self.center = np.array([self.simparams[2], self.simparams[3], self.simparams[4]])
        self.center_arr = [self.center]
        self.Rvir = self.simparams[5]
        self.Rvir_arr = [self.Rvir]
        self.dspath = self.simparams[6]
        self.dspath_arr = [self.dspath]
        self.a0 = "0."+self.dspath.split("a0.")[1][:3]
        self.a0_arr = [self.a0]
        self.L = np.array([self.simparams[7], self.simparams[8], self.simparams[9]])
        self.L_arr = [self.L]
        self.code_unit_in_kpc = self.simparams[10]
        self.conversion_arr = [self.code_unit_in_kpc]
        try:
            self.Mvir = float(parse_vela_metadata.dict_of_vela_info("Mvir")[self.simname][self.a0])
            self.Mvir_arr = [self.Mvir]
            self.gas_Rvir = float(parse_vela_metadata.dict_of_vela_info("gas_Rvir")[self.simname][self.a0])
            self.gas_Rvir_arr = [self.gas_Rvir]
            self.star_Rvir = float(parse_vela_metadata.dict_of_vela_info("star_Rvir")[self.simname][self.a0])
            self.star_Rvir_arr = [self.star_Rvir]
            self.dm_Rvir = float(parse_vela_metadata.dict_of_vela_info("dm_Rvir")[self.simname][self.a0])
            self.dm_Rvir_arr = [self.dm_Rvir]
            self.sfr = float(parse_vela_metadata.dict_of_vela_info("SFR")[self.simname][self.a0])
            self.sfr_arr = [self.sfr]
            
            
            aDict = parse_vela_metadata.dict_of_vela_info("a")[self.simname].keys()
            aDict.sort()
            self.final_a0 = float(aDict[-1])
            self.final_a0_arr = [self.final_a0]
        except:
            "test"

    def get_criteria_at_a(self, a, criteria):
        if float(self.final_a0) < float(a):
            print "Inputted a value that exceeds the greatest 'a' value in %s" %(self.simname)
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

    def create_QSO_endpoints(self, R, n_th, n_phi, n_r, rmax, length,\
                             distances = "kpc", overwrite = False, endonsph = False):
        if not overwrite and self.scanparams:
            print("overwrite is FALSE, set to TRUE to create new scan.")
            return None
        r_arr = np.linspace(0,rmax,n_r)
        th_arr = np.linspace(0,np.pi,n_th,endpoint = False)
        phi_arr = np.linspace(0,2*np.pi,n_phi,endpoint = False)
        if distances == "kpc":
            convert = self.ds.length_unit.in_units('kpc').value
        elif distances == "Rvir":
            convert = self.ds.length_unit.in_units('kpc').value
            convert /= self.Rvir
        else:
            convert = 1
        R /= convert
        r_arr /= convert
        self.scanparams = [None]*7
        self.scanparams[0] = R
        self.scanparams[1] = len(th_arr)
        self.scanparams[2] = len(phi_arr)
        self.scanparams[3] = len(r_arr)
        self.scanparams[4] = rmax
        self.scanparams[5] = length
        self.scanparams[6] = 0
        self.add_extra_scanparam_fields()
        
        self.info = np.zeros((int(length),11+len(self.ions)+1))-1.0
        weightth = weights(th_arr, "sin")
        weightr = weights(r_arr, "lin")
        L = self.L
        rot_matrix = get_rotation_matrix(L)
        for i in range(int(length)):
            theta = np.random.choice(th_arr,p = weightth)
            r = np.random.choice(r_arr,p = weightr)
            phi= np.random.choice(phi_arr)
            alpha = 2*np.pi*np.random.random()
            self.info[i][:5] = np.array([i,theta,phi,r,alpha])
            self.info[i][5:8] = np.matmul(rot_matrix, ray_endpoints_spherical(R,r,theta,phi,alpha,endonsph)[0]) + self.center
            self.info[i][8:11] = np.matmul(rot_matrix, ray_endpoints_spherical(R,r,theta,phi,alpha,endonsph)[1]) + self.center 
        print(str(length)+" LOSs to scan.")
        return length

    def get_coldens(self, save = 10, parallel = False, test = False,comm = None):
        tosave = save
        starting_point = self.length_reached
        if not parallel:
            for vector in self.info[starting_point:]:
                self.scanparams[6]+=1
                self.length_reached = self.scanparams[6]
                print("%s/%s"%(self.length_reached,self.length))
                vector = _get_coldens_helper(self.ds,vector,self.ions)
                tosave -= 1
                if tosave == 0 and not test :
                    output = self.save_values()
                    print("file saved to "+output+".")
                    tosave = save
        if parallel:
            bins = np.append(np.arange(starting_point,self.length,save),self.length)
            for i in range(0, len(bins)-1):
                current_info = self.info[bins[i]:bins[i+1]]
                print("%s-%s /%s"%(bins[i],bins[i+1],len(self.info)))
                for process in range(0,comm.Get_size()-1):
                    try:
                        infoToSend = current_info[process]
                    except:
                        infoToSend = None
                    toSend = (infoToSend,self.ions)
                    comm.send(toSend, dest = process+1)
                for process in range(0,comm.Get_size()-1):
                    new_info = comm.recv(source = process+1)
                    if bins[i] + process >= self.length:
                        continue
                    self.info[bins[i]+process] = new_info
                self.scanparams[6]+=(bins[i+1]-bins[i])
                self.length_reached = self.scanparams[6]
                if not test:
                    output = self.save_values()
                    print("file saved to "+output+".")
            if not test:
                output = self.save_values()
        if not test:
            print("file saved to "+output+".")
        return self.info
    
    def save_values(self,dest = None):
        if len(self.info[0]) <= 11:
            print("No ions!")
        linesfinished = self.length_reached
        numlines = self.length
        redshift = self.rounded_redshift
        simname = self.simname
        ionsstr = ""
        for ion in self.ions:
            ionsstr += "_"+ion.replace(" ","")
        if dest:
            filename = dest
        else:
            foldername = "quasarscan/output/"+simname+"coldensinfo"
            if not os.path.exists(foldername):
                os.makedirs(foldername)
            specificfilename = "%s_of_%s-"%(str(linesfinished),str(numlines)) +ionsstr+"_z"+str(redshift)[:4]+".txt"
            filename = foldername+"/"+specificfilename
            prev = os.listdir(foldername)
            for item in prev:
                if item.endswith("of_%s-"%str(numlines) +ionsstr+"_z"+str(redshift)[:4]+".txt"):
                    os.remove(foldername+"/"+item)
        f = open(filename,"w+")
        firstline = "[dsname, z, center[0], center[1], center[2], Rvir, pathname]\n"
        secondline = str(self.simparams)+"\n"
        thirdline = "[R, n_th, n_phi, n_r, r_max, num_lines, line_reached]\n"
        fourthline = str(self.scanparams)+"\n"
        fifthline = "ions\n"
        sixthline = "["+str(self.ions[0])
        for ion in self.ions[1:]:
            sixthline += ", "+ion
        f.write(firstline)
        f.write(secondline)
        f.write(thirdline)
        f.write(fourthline)
        f.write(fifthline)
        f.write(sixthline+"]\n")
        for vector in self.info:
            f.write(str(vector).replace("\n",""))
            f.write("\n")
        f.close()
        return filename
    
def read_values(filename):
    """ firstline = "[dsname, z, center[0], center[1], center[2], Rvir, pathname]\n"
        secondline = str(self.simparams)+"\n"
        thirdline = "[R, n_th, n_phi, n_r, r_max, num_lines, line_reached]\n"
        fourthline = str(self.scanparams)+"\n"
        fifthline = "ions"
        sixthline = "["+str(self.ions[0])+", "+...+"]"
    """
    f = open(filename)
    firstline = f.readline()
    secondline = f.readline()[:-1]
    thirdline = f.readline()
    fourthline = f.readline()[:-1]
    fifthline = f.readline()
    sixthline = f.readline()[:-1]
    simparams = eval(secondline)
    scanparams = eval(fourthline)
    ions = sixthline[1:-1].split(", ")
    length = scanparams[5]
    data = np.zeros((int(length),11+len(ions)+1))
    for i in range(length):
        myline = f.readline()[1:-1]
        data[i] = np.fromstring(myline,sep = " ")
    return simparams,scanparams,ions,data

def _get_coldens_helper(ds,vector,ions,rank = 0):
    try:
        print("<rank %d, starting process> "%rank)
        ident = str(rank)
        start = vector[5:8]
        end = vector[8:11]  
        print("saving ray to %s"%"ray"+ident+".h5") 
        ray = trident.make_simple_ray(ds,
                    start_position=start,
                    end_position=end,
                    data_filename="ray"+ident+".h5",
                    fields = [('gas',"metallicity")],
                    ftype='gas')
        trident.add_ion_fields(ray,ions)
        field_data = ray.all_data()
        for i in range(len(ions)):
            ion = ions[i]
            cdens = np.sum(field_data[("gas",ion_to_field_name(ion))] * field_data['dl'])
            #outcdens = np.sum((field_data['radial_velocity']>0)*field_data[ion_to_field_name(ion)]*field_data['dl'])
            #incdens = np.sum((field_data['radial_velocity']<0)*field_data[ion_to_field_name(ion)]*field_data['dl'])
            vector[11+i] = cdens
            #vector[12+3*i+1] = outcdens
            #vector[12+3*i+2] = incdens
        Z = np.average(field_data[('gas',"metallicity")],weights=field_data['dl'])
        vector[-1] = Z
    except Exception:
        logging.exception("failed")
    try:
        os.remove("ray"+ident+".h5")
    except:
        pass 
    print("vector = "+str(vector))
    return vector


def convert_a0_to_redshift(a0):
    return 1.0/float(a0)-1

def read_Choi_metadata():
    files = [948,908,858,763,721,664,616,549,501,408,380,329,305,290,259,227,224,220,215,209,204,190,189,175,163,162,125,53]
    all_data_dict = {}
    for num in files:
        filename = "/Users/claytonstrawn/Downloads/Choi17_z1set/m0{0:03d}_info_044.txt".format(num)
        f = open(filename)
        lines = f.readlines()
        data_dict = {}
        for line in lines:
            if len(line.split(":"))>1:
                key,value = line.split(":")
                data_dict[key] = value
        all_data_dict[num] = data_dict
    return files,all_data_dict

def get_vela_metadata(simname,a0):
    try:
        Rvir = float(parse_vela_metadata.Rdict[simname][a0])
        Lstrings = parse_vela_metadata.Ldict[simname][a0]
    except KeyError:
        Rvir = -1.0
        Lstrings = [0,0,1]
    L = np.zeros(3)
    L[0] = float(Lstrings[0])
    L[1] = float(Lstrings[1])
    L[2] = float(Lstrings[2])
    if np.isnan(L).any():
        print("VELA L data not found, returning None.")
        L = np.array([0,0,1])
    elif Rvir == 0.0:
        print("VELA Rvir data not found, returning None.")
        Rvir = -1.0
    return Rvir,L

def read_command_line_args(args, shortform,longform, tograb, defaults = None):
    if shortform in args or longform in args:
        if tograb == 0:
            return 1
        if shortform in args:
            i = args.index(shortform)
        else:
            i = args.index(longform)
        paramsstr = args[i+1:i+1+tograb]
        params = [None]*len(paramsstr)
        for i in range(len(defaults)):
            if type(defaults[i]) is str:
                params[i] = paramsstr[i]
            else:
                toinsert = eval(paramsstr[i])
                if type(defaults[i]) is type(toinsert):
                    params[i] = toinsert
                else:
                    print("Called with incorrect arguments: \
                        The %dth argument after '%s' or '%s' should be type: %s"%
                        (i, shortform, longform, type(defaults[i])))
                    return None
        return params
    else: 
        return defaults

def get_filename_from_simname(simname,redshift):
    dirname = 'quasarscan/output/%scoldensinfo/'%simname
    for _, _, filenames in os.walk(dirname):
        for file in filenames:
            z = file[-8:-4]
            if z[0] == "z":
                z = z[1:]
            if redshift == -1:
                progress = file[:10]
                print z, progress
            elif abs(redshift - float(z))< .05:
                filename = dirname+file
    return filename


def main_for_other_rank(comm, rank, dspath):
    projectdir = dspath.split("10MpcBox")[0]
    a0 = dspath.split("a0.")[1][:3]
    file_particle_header = projectdir+"PMcrda0.%s.DAT"%a0
    file_particle_data = projectdir+"PMcrs0a0.%s.DAT"%a0
    file_particle_stars = projectdir+"stars_a0.%s.dat"%a0
    ds = yt.load(dspath,file_particle_header=file_particle_header,\
                      file_particle_data=file_particle_data,\
                      file_particle_stars=file_particle_stars)
    numlines = 40
    print("Loaded additional dataset for process %d"%rank)
    size = comm.Get_size()
    cycles = numlines//(size-1) if numlines%(size-1) == 0 else numlines//(size-1)
    for i in range(cycles):
        vectorions = comm.recv(source = 0)
        vector = vectorions[0]
        ions = vectorions[1]
        myvector = _get_coldens_helper(ds,vector,ions,rank = rank)
        comm.send(myvector, dest = 0)
        if vector is None:
            comm.send(None, dest = 0)
        else:
            myvector = _get_coldens_helper(ds,vector,ions,rank = rank)
            comm.send(myvector, dest = 0)
def main_for_rank_0(simname = None, dspath = None, ions = None, Rvir = None,L = None,\
                    R = None,n_r = None,n_th = None,n_phi = None,rmax = None,length = None,distances = None,\
                    simparams = None,scanparams = None,data = None,\
                    comm = None,parallel = None,save = None):
    
    q = QuasarSphere(simname = simname, dspath = dspath, ions = ions, Rvir = Rvir,L = L,\
                     simparams = simparams, scanparams = scanparams,data = data)
    if not simparams:
        q.create_QSO_endpoints(R, n_th, n_phi, n_r, rmax, length, distances=distances)
    q.get_coldens(parallel=parallel,comm=comm,save=save)


if __name__ == "__main__":
    simname = None
    dspath = None
    ions = None
    Rvir = None
    L = None
    R = None
    n_r = None
    n_th = None
    n_phi = None
    rmax = None
    length = None
    simparams = None
    scanparams = None
    data = None
    comm = None
    parallel = None
    save = None
    
    if len(sys.argv) == 1:
        print "Using defaults (testing)"
        sys.argv = ["","n",'VELA_v2_17/10MpcBox_csf512_a0.200.d', 'VELA_v2_17',\
                    "-qp","6", "12", "12", "12", "1.5", "40","-i","[H I]","-p"]
    new = sys.argv[1]
    if new == "n":
        dspath = sys.argv[2]
        a0 = "0."+dspath.split("a0.")[1][:3]
        simname = sys.argv[3]
        if simname.lower().startswith("vela"):
            Rvir, L = get_vela_metadata(simname,a0)
            if Rvir > 0:
                distances = "Rvir"
            else:
                print("metadata not found in table! Assuming r_arr, L")
                distances = "kpc"
        elif simname.lower().startswith("choi"):
            Rvir = read_Choi_metadata()[1]["Rvir"]
            L = read_Choi_metadata()[1]["L"]
            distances = "Rvir"
        else: 
            print("simulation not found! Rvir, L not set")
            Rvir = None
            L = None
            distances = "kpc"

        if distances == "Rvir":
            defaultsphere = 6,12,12,12,1.5,400
        else:
            defaultsphere = 1000,12,12,12,250,400
        defaultions = ["[O VI, Ne VIII, H I, C III, O IV, N III, Mg II, O V, "+\
                        "O III, N IV, Mg X, N V, S IV, O II, S III, S II, S V, S VI, N II]"]
        defaultsave = [10]

        R,n_r,n_th,n_phi,rmax,length = read_command_line_args(sys.argv, "-qp","--sphereparams", 6, defaultsphere)
        save = read_command_line_args(sys.argv, "-s","--save", 1, defaultsave)[0]
        ions = read_command_line_args(sys.argv, "-i","--ions", 1, defaultions)[0]
    elif new == "c":
        filename = read_command_line_args(sys.argv, "-fn","--filename", 1, ["None"])[0]
        simname, redshift = read_command_line_args(sys.argv, "-sz","--simnameredshift", 2, ["None",-1.0])
        if filename == "None" and simname == "None":
            print "no file to load"
        elif filename == "None":
            filename = get_filename_from_simname(simname, redshift)
        simparams,scanparams,ions,data = read_values(filename)
        dspath = simparams[0]
    else: 
        print("Run this program with argument 'n' (new) or 'c' (continue).")
        
    parallelint = read_command_line_args(sys.argv, "-p","--parallel", 0)
    parallel = (parallelint == 1)
    if parallel:
        comm = MPI.COMM_WORLD
        save = comm.Get_size()-1        
        rank = MPI.COMM_WORLD.Get_rank()
        if rank == 0:
            print("Starting quasar_scan.py main script (rank 0)")
            main_for_rank_0(simname = simname, dspath = dspath, ions = ions, Rvir = Rvir,L = L,\
                    R = R,n_r = n_r,n_th = n_th,n_phi = n_phi,rmax = rmax,length = length,distances = distances,\
                    simparams = simparams,scanparams = scanparams,data = data,\
                    comm = comm,parallel = parallel,save = save)
        else:
            print("rank "+str(rank)+": activating")
            main_for_other_rank(comm, rank,dspath)
    else:
        main_for_rank_0(simname = simname, dspath = dspath, ions = ions, Rvir = Rvir,L = L,\
                    R = R,n_r = n_r,n_th = n_th,n_phi = n_phi,rmax = rmax,length = length,distances = distances,\
                    simparams = simparams,scanparams = scanparams,data = data,\
                    comm = comm,parallel = parallel,save = save)



