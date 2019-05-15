import yt

codes = ['art','ramses','gizmo','gadget','gear','enzo','tipsy']
yt_dstype_names = {'art':'art','ramses':'ramses','gizmo':'gadget_hdf5','gadget':'gadget_hdf5','gear':'gadget_hdf5','enzo':None,'tipsy':'tipsy'}

atoms = ['C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', \
        'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', \
        'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']

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
        ds.length_unit.in_units('unitary')
    else:
        ds = yt.load(path)
    return ds

def add_necessary_fields_to_ds(code,ds):
    dstype_name = ds.dataset_type
    if dstype_name not in yt_dstype_names.values():
        print "add_necessary_fields_to_ds was not prepared for the code %s!"%dstype_name
        print "please edit that file first."
        raise KeyError
    assert dstype_name == yt_dstype_names[code]
    if code == 'art':
        #yt has been updated to do ART pre-processing already by me
        pass
    elif code == 'ramses':
        def _metal_density(field, data):
            tr = data['gas','H_nuclei_density']*yt.utilities.physical_constants.mh
            tr /= data['gas','metallicity'].in_units('dimensionless')
            return tr
        ds.add_field(('gas','metal_density'),
                       sampling_type="cell",
                       function=_metal_density,
                       units='g/cm**3')
    elif code == 'gizmo':
        def _metal_density(field, data):
            tr = data['gas','H_nuclei_density']*yt.utilities.physical_constants.mh
            tr /= data['gas','metallicity'].in_units('dimensionless')
            return tr
        ds.add_field(('gas','metal_density'),
                       sampling_type="cell",
                       function=_metal_density,
                       units='g/cm**3')
    elif code == 'gadget':
        def _metal_density(field, data):
            tr = data['gas','H_nuclei_density']*yt.utilities.physical_constants.mh
            tr /= data['gas','metallicity'].in_units('dimensionless')
            return tr
        ds.add_field(('gas','metal_density'),
                       sampling_type="cell",
                       function=_metal_density,
                       units='g/cm**3')
    elif code == 'gear':
        print "code %s not implemented yet!"%code
    elif code == 'enzo':
        # no new fields needed for ENZO
        pass
    elif code == 'tipsy':
        def _gas_mass(field, data):
            return data['deposit','Gas_mass']
        def _metal_density(field, data):
            tr = data['gas','H_nuclei_density']*yt.utilities.physical_constants.mh
            tr /= data['gas','metallicity'].in_units('dimensionless')
            return tr
        try:
            ds.add_field(('gas','mass'),units = 'g', function = _gas_mass, sampling_type = 'cell')
        except:
            pass
        try:
            ds.add_field(('gas','metal_density'),
                       sampling_type="cell",
                       function=_metal_density,
                       units='g/cm**3')
        except:
            pass

def fields_to_keep_in_sightline(code,ions):
    fields_to_keep = [('gas',"density"),('gas',"mass"),('gas',"temperature"),('gas',"radial_velocity")]
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
            if atom in atoms_from_ions(ions):
                fields_to_keep.append(('gas','%s_nuclei_mass_density'%atom))
        fields_to_keep.append(('gas',"H_nuclei_density"))
    elif code == 'ramses':
        fields_to_keep.append(('gas',"metal_density"))
        fields_to_keep.append(('gas','metallicity'))
    elif code == 'gizmo':
        fields_to_keep.append(('gas',"metal_density"))
        fields_to_keep.append(('gas','metallicity'))
    elif code == 'gadget':
        fields_to_keep.append(('gas',"metal_density"))
        fields_to_keep.append(('gas','metallicity'))
    elif code == 'gear':
        print "code %s not implemented yet!"%code
    elif code == 'enzo':
        fields_to_keep.append(('gas','metal_density'))
        fields_to_keep.append(('gas','metallicity'))
    elif code == 'tipsy':
        fields_to_keep.append(('gas',"metal_density"))
        fields_to_keep.append(('gas',"metallicity"))
        fields_to_keep.append(('gas',"H_nuclei_density"))
    return fields_to_keep

def load_and_setup(path,code,ions):
    if "_" in code:
        code = code.split("_")[2]
    if code not in codes:
        print "load_and_setup was not prepared for the code %s!"%dstype_name
        print "Please edit that file first."
        raise NotImplementedError
    ds = ytload(path,code)
    try:
        assert yt_dstype_names[code] == ds.dataset_type
    except:
        print "the code stored at: %s is not of type %s, but of type %s"%(path,code,ds.dataset_type) 
        raise AssertionError
    add_necessary_fields_to_ds(code,ds)
    fields_to_keep = fields_to_keep_in_sightline(code,ions)
    return ds, fields_to_keep

