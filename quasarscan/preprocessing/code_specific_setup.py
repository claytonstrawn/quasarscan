import yt
import numpy as np
from quasarscan.utils import PI_field_defs
from quasarscan.preprocessing.add_common_fields import set_up_funcs


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
    set_up_funcs[code](ds)
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
        fields_to_keep.append(('gas','smoothing_length'))
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

