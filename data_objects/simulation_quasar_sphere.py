import numpy as np
from quasarscan3.data_objects.quasar_sphere import QuasarSphere
from quasarscan3.data_objects.gasbinning import GasBinsHolder
from quasarscan3.utils import ion_lists
import os

class SimQuasarSphere(QuasarSphere):
    def __init__(self,start_up_info_packet):
        self.number = 1
        self.type = "sim"
        simparams  = start_up_info_packet[0]
        scanparams = start_up_info_packet[1]
        ions       = start_up_info_packet[2]
        data       = start_up_info_packet[3]
        gasbins    = start_up_info_packet[4]
        self.simparams = simparams
        self.scanparams = scanparams
        self.add_extra_scanparam_fields()
        self.ions = ions
        self.info = data
        self.gasbins = gasbins
        super().__init__(simparams[0],simparams[1])

        #renames basic scanparams data into new instance variables
    def add_extra_scanparam_fields(self):
        self.R = self.scanparams[0]
        self.len_th_arr = self.scanparams[1]
        self.len_phi_arr = self.scanparams[2]
        self.len_r_arr = self.scanparams[3]
        self.rmax = self.scanparams[4]
        self.length = self.scanparams[5] 
        self.length_reached = self.scanparams[6]

    def save_values(self,test = False,oldfilename=None):
        linesfinished = self.length_reached
        numlines = self.length
        rounded_redshift = self.rounded_redshift
        simname = self.fullname
        ionsstr = ion_lists.filename_stringform(self.ions)
        if ionsstr in ion_lists.dict_of_ionstrsfilename.keys():
            ionsstr = ion_lists.dict_of_ionstrsfilename[ionsstr]
        foldername1 = "quasarscan_data"
        foldername = foldername1+'/output/'+simname+"coldensinfo"
        assert os.path.exists(foldername1), 'need a quasarscan_data folder'
        if not os.path.exists(foldername):
            os.makedirs(foldername)
        specificfilename = "%s_of_%s-"%(str(linesfinished),str(numlines)) +ionsstr+"_z"+str(rounded_redshift)[:4]+".txt"
        filename = os.path.join(foldername,specificfilename)
        if test:
            print(filename)
            return
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
        for i in range(7):
	        f.write(firstfew[i])
        for vector in self.info:
            f.write(str(vector).replace("\n",""))
            f.write("\n")
        f.close()
        print("saved file %s"%filename)
        if oldfilename is not None:
            os.remove(oldfilename)
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
    gasbins = GasBinsHolder(string = firstfew[5])
    length = scanparams[5]
    num_extra_columns = gasbins.get_length()
    data = np.zeros((int(length),11+len(ions)*(num_extra_columns+2)+3))
    for i in range(length):
        myline = f.readline()[1:-1]
        if myline.strip('\n \t') == '':
            continue
        data[i] = np.fromstring(myline,sep = " ")
    f.close()
    return simparams,scanparams,ions,data,gasbins

