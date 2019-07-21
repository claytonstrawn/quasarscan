import numpy as np
import os
import sys
try:
    from quasarscan import quasar_sphere
    from quasarscan import gasbinning
    from quasarscan.quasar_sphere import GeneralizedQuasarSphere
    level = 0
except:
    import quasar_sphere
    from quasar_sphere import GeneralizedQuasarSphere
    import gasbinning
    
    level = 1

#input: level (0 meaning working directory is '~', 1 meaning, working directory is '~/quasarscan')
#output: list of strings, names of the files where observations are being stored
def get_all_textfiles(level):
    return files

#input: filename to go into
#output: a line of headers representing columns, and the lines of the file
def read_textfile(filename):
    f = open(filename)
    headers = f.readline()[:-1]
    f.readline()
    lines = f.read().splitlines()
    lines = list(filter(None, lines))
    return headers, lines

def create_observational_quasarArray(headers,lines):
    quasarArray = []
    for line in lines:
        quasarArray.append(Observation(headers,line))
    return quasarArray

class Observation(GeneralizedQuasarSphere):
    # input: the info from read_textfile
    # output: None
    # action: read the file, add all of the relevant attributes listed here to the object
    def __init__(self, header, line_from_file):
        parsed_line = line_from_file.split(",") 
        header = header.replace('"','').split(",");
        self.number = 1 # 1 galaxy
        self.type = "Observation"
        self.length = 1 # 1 sightline
        self.length_reached = 1
        self.gasbins = gasbinning.GasBinsHolder()
        self.ions = [header[11][4:header[11].find(":")]]
        
        if header[11][0:3] == "log":
            self.info = np.array([[-1, -1, -1, float(parsed_line[5]), -1, -1, -1,-1,-1, -1, -1,
                                   10**float(parsed_line[11]), float(parsed_line[12]), -1, -1]])
        else:        
            self.info = np.array([[-1, -1, -1, float(parsed_line[5]), -1, -1, -1,-1,-1, -1, -1,
                                   float(parsed_line[11]), float(parsed_line[12]), -1, -1]])
        
        
        #self.info needs to have a particular form:
        """
        np.array([[i, theta, phi, r, alpha, x1, y1, z1, x2, y2, z2, ion1:cdens, ion1:eb, ion2:cdens, ion2:eb ...]])
        """
        #this form can be changed up *except* info[:,1], info[:,2] and info[:,3] must be theta, phi, and r respectively 
        #(and we probably won't know theta/phi)
        self.name     = parsed_line[0]
        self.simname  = None
        self.version  = None
        self.code     = None
        self.simnum   = None
        self.Rvir_is_real = None
        self.code_unit = None
        self.cgm_name= parsed_line[0]
        self.sightline = parsed_line[1]
        #above are strings
        
        #below are numbers
        self.redshift = float(parsed_line[4])
        self.rounded_redshift = None
        self.center = None
        self.Rvir = float(parsed_line[8])
        self.dspath = None
        self.a0 = None
        self.L = None
        self.L_mag = None
        self.conversion_factor = None
        self.Mvir = None
        self.gas_Rvir = None
        self.Mgas = self.gas_Rvir
        self.star_Rvir = float(parsed_line[6])
        self.Mstar = self.star_Rvir
        self.dm_Rvir = None
        self.Mdm = self.dm_Rvir
        self.sfr = float(parsed_line[9])
        self.ssfr = None
        self.final_a0 = None
        return
    
    #input: the ion you are looking for as a string, either like "Ne VIII" or "Ne VIII:cdens" which should both return the column density
    #output: the position i, i.e. you'll find that value in info[0,i]
    def get_ion_column_num(self,ion):
        if not ":" in ion:
            bintype = "cdens"
        else:
            i = ion.index(":")
            bintype = ion[i+1:]
            ion = ion[:i]
        if bintype == "cdens":
            plus = 0
        elif bintype == "eb":
            plus = 1
        try:
            return 11 + self.ions.index(ion)*(2)+plus
        except:
            print("Ion %s not found in this Sphere. Please try restricting to ions = ['%s']."%(ion,ion))
    
#input: somewhere data is
#output: None
#action: save a text file which Observation can read later
def read_from_table(pdf_file):
    #might be interactive? or require pdf crawling
    #which are tricky things
    return 