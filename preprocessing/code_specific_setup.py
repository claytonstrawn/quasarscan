import yt
import numpy as np
from quasarscan.utils import PI_field_defs
from yt.utilities.physical_constants import mh
from yt.utilities.cosmology import Cosmology
from yt.data_objects.particle_filters import add_particle_filter

codes = ['art','ramses','gizmo','gadget','gear','enzo','tipsy','changa']
sphcodes = ['gizmo','gadget','gear','tipsy','changa']
yt_dstype_names = {'art':'art','ramses':'ramses','gizmo':'gadget_hdf5','gadget':'gadget_hdf5','gear':'gadget_hdf5','enzo':None,'tipsy':'tipsy','changa':'tipsy'}

atoms = ['C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', \
        'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', \
        'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']

all_ions_w_known_PI_defs = PI_field_defs.make_funcs()[0]

#summary: ART does not store all of the data in a single snapshot file, 
#         but several files with similar names. At least some of these are
#         required in order to know it is a cosmological simulation
#
#inputs: dspath: file location of main ('.d') simulation file
#
#outputs: tuple of names of all other necessary files
def get_aux_files_art(dspath):
    projectdir = dspath.split("10MpcBox")[0]
    if "a0." in dspath:
        a0 = dspath.split("a0.")[1][:3]
        file_particle_header = projectdir+"PMcrda0.%s.DAT"%a0
        file_particle_data = projectdir+"PMcrs0a0.%s.DAT"%a0
        file_particle_stars = projectdir+"stars_a0.%s.dat"%a0
    else:
        timestep = dspath.split("_")[-1][:-2]
        file_particle_header = projectdir+"PMcrd_%s.DAT"%timestep
        file_particle_data = projectdir+"PMcrs0_%s.DAT"%timestep
        file_particle_stars = projectdir+"stars_%s.dat"%timestep 
    return file_particle_header,file_particle_data,file_particle_stars   

#summary: wrapper for yt.load that will make sure to put in ART
#         filenames correctly
#
#inputs: path: file location of main file, can be passed to yt.load
#        code: which kind of AGORA code is under analysis here
#
#outputs: ds: yt DataSet object
def ytload(path,code):
    if code == 'art':
        h,d,s = get_aux_files_art(path)
        ds = yt.load(path,file_particle_header=h,\
                                  file_particle_data=d,\
                                  file_particle_stars=s)
        ds.length_unit.in_units('unitary')
    else:
        ds = yt.load(path)
    return ds



#summary: goes through each nonnative field necessary for quasarscan
#         analysis and adds it to the dataset. This includes
#         X_nuclei_density for several atoms X and a general
#         metal_density (as opposed to metallicity)
#
#inputs: code: which kind of AGORA code is under analysis here
#        ds: yt DataSet object
#        add_pi_fracs: whether we want to add PI/CI field definitions 
#                     as defined in Strawn et al. (2020) as binary 
#                     fields (1 if PI/CI, 0 otherwise). Default False
#
#outputs: None
def add_necessary_fields_to_ds(code,ds,add_pi_fracs=True):
    if code not in codes:
        print("add_necessary_fields_to_ds was not prepared for the code %s!"%code)
        print("please edit that file first.")
        raise KeyError
    assert ds.dataset_type == yt_dstype_names[code]
    if add_pi_fracs:
        PI_field_defs.make_funcs(ds=ds,add_fields=True)
    if code == 'art':
        set_up_art(ds)
    elif code == 'ramses':
        set_up_ramses(ds)
    elif code == 'gizmo':
        set_up_gizmo(ds)
    elif code == 'gadget':
        pass
    elif code == 'changa':
        set_up_changa(ds)
    elif code == 'gear':
        set_up_gear(ds)
    elif code == 'enzo':
        set_up_enzo(ds)
    elif code == 'tipsy':
        pass
    check_necessary_fields_exist(ds)


def add_radial_distance_fields(ds,center):
    for i,ax in enumerate('xyz'):
        def radial_distance_ax(field,data):
            return data['gas',ax]-center[i]
        ds.add_field(('gas','radial_distance_%s'%ax),
                   sampling_type='cell',
                   function=radial_distance_ax,
                   units='cm',force_override = True)
        
def check_necessary_fields_exist(ds):
    necessary = [('stars','particle_mass'),\
                 ('darkmatter','particle_mass'),\
                 ('stars','particle_creation_time'),\
                 ('gas','mass'),\
                 ('gas','velocity_x'),
                 ('gas','x')]
    for f in necessary:
        assert f in ds.derived_field_list

def set_up_art(ds):
    #art is used as the archetype for all codes
    #its default field names are what all other codes will
    #implement
    pass
    
def set_up_enzo(ds):
    def star_filter(pfilter, data):
        return ((data[(pfilter.filtered_type, 'creation_time')] > 0))
    add_particle_filter('stars', function=star_filter, filtered_type='all', requires=['creation_time'])
    ds.add_particle_filter('stars')
    def dm_filter(pfilter, data):
        return ((data[(pfilter.filtered_type, 'creation_time')] == 0))
    add_particle_filter('darkmatter', function=dm_filter, filtered_type='all', requires=['creation_time'])
    ds.add_particle_filter('darkmatter')
    
    def particle_creation_time(field,data):
        return data['stars','creation_time']
    ds.add_field(('stars','particle_creation_time'),
                sampling_type = 'particle',
                function = particle_creation_time,
                units = 's')
    
def set_up_ramses(ds):
    def star_mass_rename(field,data):
        return data['star','particle_mass']
    ds.add_field(('stars','particle_mass'),
                sampling_type = 'particle',
                function = star_mass_rename,
                units = 'g')
    def dm_mass_rename(field,data):
        return data['DM','mass']
    ds.add_field(('darkmatter','particle_mass'),
                sampling_type = 'particle',
                function = star_mass_rename,
                units = 'g')
    def particle_creation_time(field,data):
        return ds.current_time.in_units('s') - data['star','age'].in_units('s')
    ds.add_field(('stars','particle_creation_time'),
                sampling_type = 'particle',
                function = particle_creation_time,
                units = 's')
    
    def metal_density(field,data):
        return data['gas','metal_mass']/data['gas','cell_volume']
    ds.add_field(('gas','metal_density'),
                sampling_type = 'cell',
                function = metal_density,
                units = 'g/cm**3')
    
def set_up_changa(ds):
    def star_mass_rename(field,data):
        return data['Stars','particle_mass']
    ds.add_field(('stars','particle_mass'),
                sampling_type = 'particle',
                function = star_mass_rename,
                units = 'g')
    def dm_mass_rename(field,data):
        return data['DarkMatter','mass']
    ds.add_field(('darkmatter','particle_mass'),
                sampling_type = 'particle',
                function = star_mass_rename,
                units = 'g')
    def particle_creation_time_rename(field,data):
        return data['Stars','creation_time']
    ds.add_field(('stars','particle_creation_time'),
                sampling_type = 'particle',
                function = particle_creation_time_rename,
                units = 's')
        
    def metal_density(field,data):
        return data['gas','metallicity']*data['gas','density']
    ds.add_field(('gas','metal_density'),
                sampling_type = 'cell',
                function = metal_density,
                units = 'g/cm**3')
    
def set_up_gear(ds):
    def star_mass_rename(field,data):
        return data['PartType1','particle_mass']
    ds.add_field(('stars','particle_mass'),
                sampling_type = 'particle',
                function = star_mass_rename,
                units = 'g')
    
    ad = ds.all_data()
    m2 = ad['PartType2','Masses']
    low,high = np.amin(m2),np.amax(m2)
    m5 = ad['PartType5','Masses']
    unique_dm5_masses = np.unique(m5)
    def dm_filter(pfilter, data):
        abovemin = data[(pfilter.filtered_type, 'Masses')]>=low
        belowmax = data[(pfilter.filtered_type, 'Masses')]<=high
        allowed = np.logical_and(abovemin,belowmax)
        for v in unique_dm5_masses:
            allowed = np.logical_or(allowed,data[(pfilter.filtered_type, 'Masses')]==v)
        return allowed
    add_particle_filter('darkmatter', function=dm_filter, filtered_type='all', requires=['Masses'])
    ds.add_particle_filter('darkmatter')
    
    co = Cosmology(hubble_constant=ds.hubble_constant, omega_matter=ds.omega_matter,
               omega_lambda=ds.omega_lambda, omega_curvature=0.0)
    def particle_creation_time_rename(field,data):
        return co.t_from_a(data['PartType1','StarFormationTime'])
    ds.add_field(('stars','particle_creation_time'),
                sampling_type = 'particle',
                function = particle_creation_time_rename,
                units = 's')
    
def set_up_gizmo(ds):
    def star_mass_rename(field,data):
        return data['PartType4','particle_mass']
    ds.add_field(('stars','particle_mass'),
                sampling_type = 'particle',
                function = star_mass_rename,
                units = 'g')
    
    ad = ds.all_data()
    m1 = ad['PartType1','Masses']
    unique_dm1_masses = np.unique(m1)
    m2 = ad['PartType2','Masses']
    unique_dm2_masses = np.unique(m2)    
    def dm_filter(pfilter, data):
        allowed = np.zeros(data[(pfilter.filtered_type, 'Masses')].shape,dtype=bool)
        for v in ds.arr(np.concatenate([unique_dm1_masses,unique_dm2_masses]).v,m2.units):
            allowed = np.logical_or(allowed,data[(pfilter.filtered_type, 'Masses')]==v)
        return allowed
    add_particle_filter('darkmatter', function=dm_filter, filtered_type='all', requires=['Masses'])
    ds.add_particle_filter('darkmatter')
    
    co = Cosmology(hubble_constant=ds.hubble_constant, omega_matter=ds.omega_matter,
               omega_lambda=ds.omega_lambda, omega_curvature=0.0)
    def particle_creation_time_rename(field,data):
        return co.t_from_a(data['PartType4','StellarFormationTime'])
    ds.add_field(('stars','particle_creation_time'),
                sampling_type = 'particle',
                function = particle_creation_time_rename,
                units = 's')
    
    def metal_density(field,data):
        return data['gas','metallicity']*data['gas','density']
    ds.add_field(('gas','metal_density'),
                sampling_type = 'cell',
                function = metal_density,
                units = 'g/cm**3')
    
    
#summary: get list of fields which need to be passed to TRIDENT
#         and told to forward to the sightlines created by 
#         trident.make_simple_ray
#
#inputs: code: which kind of AGORA code is under analysis here
#        ions: which ions we will attempt to track
#        add_pi_fracs: whether we want to add PI/CI field definitions 
#                     as defined in Strawn et al. (2020) as binary 
#                     fields (1 if PI/CI, 0 otherwise). Default False
#
#outputs: fields_to_keep: list of yt field names.
def fields_to_keep_in_sightline(code,ions,add_pi_fracs=True):
    fields_to_keep = [('gas',"density"),('gas',"mass"),('gas',"temperature"),('gas',"radial_velocity")]
    def atoms_from_ions(ions):
        toret = []
        for ion in ions:
            toret.append(ion.split(" ")[0])
        return list(set(toret))
    if add_pi_fracs:
        for ion in ions:
            if ion not in all_ions_w_known_PI_defs:
                continue
            PI_field_name = ('gas','PI_%s'%(ion.replace(' ','')))
            CI_field_name = ('gas','CI_%s'%(ion.replace(' ','')))
            fields_to_keep += [PI_field_name,CI_field_name]
    if code == 'art':
        fields_to_keep.append(('gas',"metal_density"))
        fields_to_keep.append(('gas', 'cell_volume'))
        for atom in atoms:
            if atom in atoms_from_ions(ions):
                fields_to_keep.append(('gas','%s_nuclei_mass_density'%atom))
        fields_to_keep.append(('gas',"H_nuclei_density"))
    elif code == 'enzo':
        fields_to_keep.append(('gas', 'cell_volume'))
        fields_to_keep.append(('gas','metal_density'))
        fields_to_keep.append(('gas','metallicity'))
    elif code == 'ramses':
        fields_to_keep.append(('gas', 'cell_volume'))
        fields_to_keep.append(('gas',"metal_density"))
        fields_to_keep.append(('gas','metallicity'))
    elif code == 'changa':
        fields_to_keep.append(('gas','smoothing_length'))
        fields_to_keep.append(('gas','metal_density'))
        fields_to_keep.append(('gas','metallicity'))
    elif code == 'gizmo':
        fields_to_keep.append(('gas','smoothing_length'))
        fields_to_keep.append(('gas',"metal_density"))
        fields_to_keep.append(('gas','metallicity'))
    elif code == 'gadget':
        fields_to_keep.append(('gas',"metal_density"))
        fields_to_keep.append(('gas','metallicity'))
    elif code == 'gear':
        fields_to_keep.append(('gas','smoothing_length'))
        fields_to_keep.append(('gas',"metal_density"))
        fields_to_keep.append(('gas','metallicity'))
    elif code == 'tipsy':
        fields_to_keep.append(('gas',"metal_density"))
        fields_to_keep.append(('gas',"metallicity"))
        fields_to_keep.append(('gas',"H_nuclei_density"))
        if 'O' in atoms_from_ions(ions):
            fields_to_keep.append(('gas',"O_nuclei_mass_density"))
    else:
        print("set_up_fields_for_sims was not prepared for the code %s!"%code)
        print("please edit that file first.")
        raise KeyError
    return fields_to_keep

#summary: wrapper for ytload which also sets up all necessary fields
#         
#inputs: path: file location of main file, can be passed to yt.load
#        code: which kind of AGORA code is under analysis here
#        ions: which ions we will attempt to track
#        add_pi_fracs: whether we want to add PI/CI field definitions 
#                     as defined in Strawn et al. (2020) as binary 
#                     fields (1 if PI/CI, 0 otherwise). Default False
#
#outputs: ds: yt DataSet object of simulation
#         fields_to_keep: list of yt field names.
def load_and_setup(path,code,ions=None,add_pi_fracs=True):
    if "_" in code:
        code = code.split("_")[2]
    if code not in codes:
        print("load_and_setup was not prepared for the code %s!"%dstype_name)
        print("Please edit that file first.")
        raise NotImplementedError
    ds = ytload(path,code)
    try:
        assert yt_dstype_names[code] == ds.dataset_type
    except:
        print("the code stored at: %s is not of type %s, but of type %s"%(path,code,ds.dataset_type))
        raise AssertionError
    add_necessary_fields_to_ds(code,ds,add_pi_fracs=add_pi_fracs)
    if ions is not None:
        fields_to_keep = fields_to_keep_in_sightline(code,ions,add_pi_fracs=add_pi_fracs)
    else:
        fields_to_keep = []
    return ds, fields_to_keep

