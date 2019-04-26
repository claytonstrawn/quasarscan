import numpy as np
import os
import sys
try:
    from quasarscan import quasar_sphere
    level = 0
except:
    import quasar_sphere
    level = 1

#input: level (0 meaning working directory is '~', 1 meaning, working directory is '~/quasarscan')
#output: list of strings, names of the files where observations are being stored
def get_all_textfiles(level):
    return files

#input: filename to go into
#output: a line of headers representing columns, and the lines of the file
def read_textfile(filename):
    return headers, lines

class Observation(quasar_sphere.GeneralizedQuasarSphere):
    # input: the info from read_textfile
    # output: None
    # action: read the file, add all of the relevant attributes listed here to the object
    def __init__(self, header, line_from_file):
        self.number = 1 # 1 galaxy
        self.length = 1 # 1 sightline
        self.gasbins = None
        self.ions = None
        self.info = None #self.info needs to have a particular form:
        """
        np.array([[i, theta, phi, r, alpha, x1, y1, z1, x2, y2, z2, ion1:cdens, ion1:eb, ion2:cdens, ion2:eb ...]])
        """
        #this form can be changed up *except* info[:,1], info[:,2] and info[:,3] must be theta, phi, and r respectively 
        #(and we probably won't know theta/phi)
        self.name     = None
        self.simname  = None
        self.version  = None
        self.code     = None
        self.simnum   = None
        self.Rvir_is_real = None
        self.code_unit = None
        self.cgm_name=None
        self.sightline =None
        #above are strings
        
        #below are numbers
        self.redshift = None
        self.rounded_redshift = None
        self.center = None
        self.Rvir = None
        self.dspath = None
        self.a0 = None
        self.L = None
        self.L_mag = None
        self.conversion_factor = None
        self.Mvir = None
        self.gas_Rvir = None
        self.Mgas = self.gas_Rvir
        self.star_Rvir = None
        self.Mstar = self.star_Rvir
        self.dm_Rvir = None
        self.Mdm = self.dm_Rvir
        self.sfr = None
        self.ssfr = None
        self.final_a0 = None
        return
    
    #input: the ion you are looking for as a string, either like "Ne VIII" or "Ne VIII:cdens" which should both return the column density
    #output: the position i, i.e. you'll find that value in info[0,i]
    def get_ion_column_num(self,ion):
        return num
    
#input: somewhere data is
#output: None
#action: save a text file which Observation can read later
def read_from_table(pdf_file):
    #might be interactive? or require pdf crawling
    #which are tricky things
    return 