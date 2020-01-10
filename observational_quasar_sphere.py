import numpy as np
import os
import sys
try:
    from quasarscan import quasar_sphere
    from quasarscan import gasbinning
    from quasarscan.quasar_sphere import GeneralizedQuasarSphere
    from quasarscan.multi_quasar_sphere_plotter import sort_ions
    level = 0
except:
    import quasar_sphere
    from quasar_sphere import GeneralizedQuasarSphere
    import gasbinning
    from multi_quasar_sphere_plotter import sort_ions
    level = 1

#input: level (0 meaning working directory is '~', 1 meaning, working directory is '~/quasarscan')
#output: list of strings, names of the files where observations are being stored
def get_all_textfiles(level):
    return files

#input: filename to go into
#output: a line of headers representing columns, and the lines of the file
def read_textfile(filename):
    f = open(filename)
    papername = f.readline()
    f.readline()
    f.readline()
    headers = f.readline()[:-1]
    f.readline()
    lines = f.read().splitlines()
    lines = list(filter(None, lines))
    return papername.replace('\n',''),headers, lines

def access_data_values(variable,header_dict,parsed_line):
    strvars = ['sl','cgm','sightline']
    numvars = ['z','Rvir','Mstar','sfr','sfr:eb','r']
    if variable in strvars:
        try:
            i,islog = header_dict[variable]
            assert not islog
            return parsed_line[i]
        except:
            return None
    elif variable in numvars or 'cdens' in variable or 'eb' in variable:
        try:
            i,islog = header_dict[variable]
            if parsed_line[i] == '<':
                value = -2.0
            elif parsed_line[i] == '>':
                value = -3.0
            elif parsed_line[i] == 'nan':
                value = np.nan
            else:
                value = float(parsed_line[i])
            if islog:
                return 10.**float(value)
            else:
                return value
        except:
            return np.nan
    else:
        print("%s should not be stored in the sightline!"%variable)
        assert False

def create_observational_quasarArray(papername,headers,lines):
    quasarArray = []
    for line in lines:
        newObs = Observation(papername,headers,line)
        add = True
        #check we don't already have this sightline from a different survey
        for obs in quasarArray:
            if obs.cgm_name == newObs.cgm_name:
                add = False
                obs.simname += ', %s'%newObs.simname
                for ion in newObs.ions:
                    if ion not in obs.ions:
                        obs.ions.append(ion)
                        ion_cdens = newObs.info[:,newObs.get_ion_column_num("%s:cdens"%ion)]
                        ion_eb = newObs.info[:,newObs.get_ion_column_num("%s:eb"%ion)]
                        new_info = np.array([[ion_cdens,ion_eb]])
                        obs.info = np.concatenate((obs.info[:,:-3],new_info,obs.info[:,-3:]))
        if add:
            quasarArray.append(Observation(papername,headers,line))
    return np.array(quasarArray)

    
#input: the ion you are looking for as a string, either like "Ne VIII" or "Ne VIII:cdens" which should both return the column density
#output: the position i, i.e. you'll find that value in info[0,i]

def get_header_columns_dict(header):
    to_return = {}
    for i,name in enumerate(header):
        if name.startswith('log '):
            islog = True
            retname = name.split('log ')[1]
        else:
            islog = False
            retname = name
        to_return[retname] = i,islog
    return to_return

class Observation(GeneralizedQuasarSphere):
    # input: the info from read_textfile
    # output: None
    # action: read the file, add all of the relevant attributes listed here to the object
    def __init__(self, papername, header, line_from_file):
        parsed_line = line_from_file.split(",") 
        header = header.split(",")
        self.ions = []
        header_columns_dict = get_header_columns_dict(header)
        temp = header_columns_dict.keys()
        keylist = []
        for key in temp:
            if (key.split(':')[0] not in keylist) and (' ' in key.split(':')[0]):
                keylist.append(key.split(':')[0])
        self.ions = sort_ions(keylist)
        self.number = 1 # 1 galaxy
        self.type = "Observation"
        self.length = 1 # 1 sightline
        self.length_reached = 1
        self.gasbins = gasbinning.GasBinsHolder()
        iondata = []
        for ion in self.ions:
            iondata.append(access_data_values('%s:cdens'%ion,header_columns_dict,parsed_line))
            iondata.append(access_data_values('%s:eb'%ion,header_columns_dict,parsed_line))
        self.info = np.array([[-1, -1, -1]+
                             [access_data_values('r',header_columns_dict,parsed_line)]+
                             [-1, -1, -1,-1,-1, -1, -1] + 
                             iondata + 
                             [-1, -1, -1]])
        
        
        #self.info needs to have a particular form:
        """
        np.array([[i, theta, phi, r, alpha, x1, y1, z1, x2, y2, z2, ion1:cdens, ion1:eb, ion2:cdens, ion2:eb ...]])
        """
        #this form can be changed up *except* info[:,1], info[:,2] and info[:,3] must be theta, phi, and r respectively 
        #(and we probably won't know theta/phi)
        
        #stringparams - knowable
        self.name     = access_data_values('cgm',header_columns_dict,parsed_line)
        self.simname  = papername
        self.Rvir_is_real = 'True'
        self.cgm_name = self.name
        self.sightline = access_data_values('sl',header_columns_dict,parsed_line)
        self.compaction_stage = "unknown"

        #stringparams - unknowable (should all be None)
        self.version  = None
        self.code     = None
        self.simnum   = None
        self.code_unit = None
        self.dspath = None
        
        #numparams - knowable
        self.redshift = access_data_values('z',header_columns_dict,parsed_line)
        self.rounded_redshift = quasar_sphere.round_redshift(self.redshift)
        self.Rvir = access_data_values('Rvir',header_columns_dict,parsed_line)
        self.a0 = 1./(self.redshift+1)
        self.Mvir = 4./3*np.pi*self.Rvir**3*(200*1.5e-7) #crit_dens = 1.5e-7 M⊙ pc-3 
        self.Mstar = access_data_values('Mstar',header_columns_dict,parsed_line)
        self.star_Rvir = self.Mstar
        self.sfr = access_data_values('sfr',header_columns_dict,parsed_line)
        self.ssfr = self.sfr/self.Mstar
        
        #numparams - unknowable (should all be nan)
        self.center = np.array([np.nan]*3)
        self.L = np.array([np.nan]*3)
        self.L_mag = np.nan
        self.conversion_factor = np.nan
        self.gas_Rvir = np.nan
        self.Mgas = self.gas_Rvir
        self.dm_Rvir = np.nan
        self.Mdm = self.dm_Rvir
        self.final_a0 = np.nan
        return

    

            
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