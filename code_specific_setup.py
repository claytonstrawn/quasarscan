import yt

codes = ['art','ramses','gizmo','gadget','gear','enzo','tipsy']
yt_dstype_names = {'art':'art','ramses':'ramses','gizmo':'gizmo','gadget':'gadget','gear':'gear','enzo':'enzo','tipsy':'tipsy'}


def get_aux_files_art(dspath):
    projectdir = dspath.split("10MpcBox")[0]
    a0 = dspath.split("a0.")[1][:3]
    file_particle_header = projectdir+"PMcrda0.%s.DAT"%a0
    file_particle_data = projectdir+"PMcrs0a0.%s.DAT"%a0
    file_particle_stars = projectdir+"stars_a0.%s.dat"%a0
    return file_particle_header,file_particle_data,file_particle_stars   

def ytload(path,code):
    if code == 'art':
        h,d,s = get_aux_files_art(path)
        ds = yt.load(path,file_particle_header=h,\
                                  file_particle_data=d,\
                                  file_particle_stars=s)
    else:
        ds = yt.load(path)
    return ds

def add_necessary_fields_to_ds(ds):
    dstype_name = ds.dataset_type
    if dstype_name not in yt_dstype_names.keys():
        print "set_up_fields_for_sims was not prepared for the code %s!"%dstype_name
        print "please edit that file first."
        raise KeyError
    else:
        code = yt_dstype_names[dstype_name]
    if code == 'art':
        #yt has been updated to do ART pre-processing already by me
        pass
    elif code == 'ramses':
        print "code %s not implemented yet!"%code
    elif code == 'gizmo':
        print "code %s not implemented yet!"%code
    elif code == 'gadget':
        print "code %s not implemented yet!"%code
    elif code == 'gear':
        print "code %s not implemented yet!"%code
    elif code == 'enzo':
        print "code %s not implemented yet!"%code
    elif code == 'tipsy':
        def gas_mass(field, data):
            return data['deposit','Gas_mass']
        ds.add_field(('gas','mass'),units = 'g', function = gas_mass, sampling_type = 'cell')


def fields_to_keep_in_sightline(code):
    fields_to_keep = [('gas',"H_nuclei_density"),('gas',"density"),('gas',"mass"),('gas',"temperature"),('gas',"radial_velocity")]
    if code not in codes:
        print "set_up_fields_for_sims was not prepared for the code %s!"%code
        print "please edit that file first."
        raise KeyError
    elif code == 'art':
        fields_to_keep.append(('gas',"metal_density"))
        fields_to_keep.append(('gas', 'cell_volume'))
        def atoms_from_ions(ions):
            toret = []
            for ion in ions:
                toret.append(ion.split(" ")[0])
            return list(set(toret))
        for atom in atoms:
            if atom in atoms_from_ions(q.ions):
                fields_to_keep.append(('gas','%s_nuclei_mass_density'%atom))
    elif code == 'ramses':
        print "code %s not implemented yet!"%code
    elif code == 'gizmo':
        print "code %s not implemented yet!"%code
    elif code == 'gadget':
        print "code %s not implemented yet!"%code
    elif code == 'gear':
        print "code %s not implemented yet!"%code
    elif code == 'enzo':
        print "code %s not implemented yet!"%code
    elif code == 'tipsy':
        fields_to_keep.append(('gas','metallicity'))       
    return fields_to_keep

def load_and_setup(path,code):
    if code not in codes:
        print "set_up_fields_for_sims was not prepared for the code %s!"%dstype_name
        print "Please edit that file first."
        raise NotImplementedError    
    ds = ytload(path,code)
    try:
        assert yt_dstype_names[ds.dataset_type] == code
    except:
        print "the code stored at: %s is not of type %s, but of type %s"%(path,code,ds.dataset_type) 
        raise AssertionError
    add_necessary_fields_to_ds(ds)
    fields_to_keep = fields_to_keep_in_sightline(code)
    return ds, fields_to_keep

