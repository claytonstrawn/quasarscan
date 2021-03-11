from yt_astro_analysis.halo_analysis.api  import HaloCatalog
from yt.utilities.cosmology import Cosmology
from quasarscan.utils.utils import sphcodes
from quasarscan.preprocessing import code_specific_setup
import numpy as np

class NotEnoughHalosError(Exception):
    def __init__(self, message):
        self.message = message
        print(self.message)

def initialize_halo_catalog(ds,hc=None):
    if hc is None:
        hc = HaloCatalog(data_ds=ds, finder_method='hop')
        hc.create()
    return hc

def get_halo_center(code,ds,hc=None,ith_largest = 1):
    hc = initialize_halo_catalog(ds,hc=hc)
        #check that the resolution is ok

    ad = ds.all_data()
    if code not in sphcodes:
        best_resolution = np.min(ad['gas','cell_volume'])
    else:
        best_resolution = None
    sorted_halos = sorted(hc.catalog,key=lambda x:x['particle_mass'])
    catalog_of_acceptable_halos = []
    for i,h in enumerate(sorted_halos):
        if len(catalog_of_acceptable_halos)==ith_largest:
            return catalog_of_acceptable_halos[-1]
        center = np.array([h['particle_position_x'].in_units('unitary').v,
                            h['particle_position_y'].in_units('unitary').v,
                            h['particle_position_z'].in_units('unitary').v])
        approx_Rvir = h['virial_radius'].in_units('kpc')*1.126
        approx_halo = ds.sphere(center,approx_Rvir)
        if best_resolution is not None and np.min(approx_halo['gas','cell_volume']) == best_resolution:
            catalog_of_acceptable_halos.append(center)
    raise NotEnoughHalosError('Ran out of halos that reach max refinement!')

def find_virial_radius(ds,center):
    z = ds.current_redshift

    co = Cosmology(hubble_constant=ds.hubble_constant, omega_matter=ds.omega_matter,
                   omega_lambda=ds.omega_lambda, omega_curvature=0.0)
    rho_c = co.critical_density(z=z)
    
    r = 0
    rs = []
    densities = []
    density = np.inf
    while density > 200*rho_c:
        r+=5
        rs.append(r)
        sp = ds.sphere(center,(r,'kpc'))
        cell_mass,particle_mass = sp.quantities.total_mass()
        volume = 4/3*np.pi*ds.arr(r,'kpc')**3
        new_density = (cell_mass+particle_mass)/volume
        if new_density > density:
            print('density not decreasing!')
            break
        density = new_density.in_units('g/cm**3')
        densities.append(density)
        toprint = (r,cell_mass+particle_mass,(cell_mass+particle_mass).units,particle_mass,particle_mass.units,cell_mass,cell_mass.units)
    rs = np.array(rs)
    densities = np.array(densities)
    Rvir = np.interp(178*rho_c,np.flip(densities),np.flip(rs))
    Mvir = cell_mass+particle_mass
    return ds.arr(Rvir,'kpc'),ds.arr(Mvir,'g')


def add_all_fields(ds,center,Rvir,Mvir):
    def all_cic_mass(field,data):
        return data[('deposit', 'all_cic')]*data[('gas','cell_volume')]
    ds.add_field(('deposit','all_cic_mass'),
               sampling_type="cell",
               function=all_cic_mass,
               units='g',force_override = True)

    def particle_momentum_x(field,data):
        return data[('all', 'particle_velocity_x')]*data[('all','particle_mass')]
    ds.add_field(('deposit','particle_momentum_x'),
               sampling_type="particle",
               function=particle_momentum_x,
               units='g*cm/s',force_override = True)
    def particle_momentum_y(field,data):
        return data[('all', 'particle_velocity_y')]*data[('all','particle_mass')]
    ds.add_field(('deposit','particle_momentum_y'),
               sampling_type="particle",
               function=particle_momentum_y,
               units='g*cm/s',force_override = True)
    def particle_momentum_z(field,data):
        return data[('all', 'particle_velocity_z')]*data[('all','particle_mass')]
    ds.add_field(('deposit','particle_momentum_z'),
               sampling_type="particle",
               function=particle_momentum_z,
               units='g*cm/s',force_override = True)


    def complete_mass(field,data):
        return data[('deposit', 'all_cic_mass')]+data['gas','mass']
    ds.add_field(('deposit','complete_mass'),
               sampling_type="cell",
               function=complete_mass,
               units='g',force_override = True)

    def radial_distance_x(field,data):
        return data['index','x']-center[0]
    ds.add_field(('gas','radial_distance_x'),
               sampling_type="cell",
               function=radial_distance_x,
               units='cm',force_override = True)
    def radial_distance_y(field,data):
        return data['index','y']-center[1]
    ds.add_field(('gas','radial_distance_y'),
               sampling_type="cell",
               function=radial_distance_y,
               units='cm',force_override = True)
    def radial_distance_z(field,data):
        return data['index','z']-center[2]
    ds.add_field(('gas','radial_distance_z'),
               sampling_type="cell",
               function=radial_distance_z,
               units='cm',force_override = True)

    def momentum_x(field,data):
        return data['deposit','all_cic_mass']*data[('deposit', 'all_cic_velocity_x')]+data['gas','cell_mass']*data[('gas', 'velocity_x')]
    ds.add_field(('deposit','momentum_x'),
               sampling_type="cell",
               function=momentum_x,
               units='cm*g/s',force_override = True)
    def momentum_y(field,data):
        return data['deposit','all_cic_mass']*data[('deposit', 'all_cic_velocity_y')]+data['gas','cell_mass']*data[('gas', 'velocity_y')]
    ds.add_field(('deposit','momentum_y'),
               sampling_type="cell",
               function=momentum_y,
               units='cm*g/s',force_override = True)
    def momentum_z(field,data):
        return data['deposit','all_cic_mass']*data[('deposit', 'all_cic_velocity_z')]+data['gas','cell_mass']*data[('gas', 'velocity_z')]
    ds.add_field(('deposit','momentum_z'),
               sampling_type="cell",
               function=momentum_z,
               units='cm*g/s',force_override = True)

    def bulk_velocity_x(field,data):
        sp = ds.sphere(center,Rvir)
        total_x_mom = sp.quantities.total_quantity([('deposit','momentum_x')])
        total_x_vel = total_x_mom/Mvir
        tr = data['gas','velocity_x']*0.0+total_x_vel
    ds.add_field(('index','bulk_velocity_x'),
               sampling_type="cell",
               function=bulk_velocity_x,
               units='cm/s',force_override = True)
    def bulk_velocity_y(field,data):
        sp = ds.sphere(center,Rvir)
        total_y_mom = sp.quantities.total_quantity(('deposit','momentum_y'))
        total_y_vel = total_y_mom/Mvir
        tr = data['gas','velocity_y']*0.0+total_y_vel
    ds.add_field(('index','bulk_velocity_y'),
               sampling_type="cell",
               function=bulk_velocity_y,
               units='cm/s',force_override = True)
    def bulk_velocity_z(field,data):
        sp = ds.sphere(center,Rvir)
        total_z_mom = sp.quantities.total_quantity(('deposit','momentum_z'))
        total_z_vel = total_z_mom/Mvir
        tr = data['gas','velocity_z']*0.0+total_z_vel
    ds.add_field(('index','bulk_velocity_z'),
               sampling_type="cell",
               function=bulk_velocity_z,
               units='cm/s',force_override = True)

    def relative_gas_velocity_x(field,data):
        return data['gas','velocity_x']-data['index','bulk_velocity_x']
    ds.add_field(('gas','relative_velocity_x'),
           sampling_type="cell",
           function=relative_gas_velocity_x,
           units='cm/s',force_override = True)
    def relative_gas_velocity_y(field,data):
        return data['gas','velocity_y']-data['index','bulk_velocity_y']
    ds.add_field(('gas','relative_velocity_y'),
           sampling_type="cell",
           function=relative_gas_velocity_y,
           units='cm/s',force_override = True)
    def relative_gas_velocity_z(field,data):
        return data['gas','velocity_z']-data['index','bulk_velocity_z']
    ds.add_field(('gas','relative_velocity_z'),
           sampling_type="cell",
           function=relative_gas_velocity_z,
           units='cm/s',force_override = True)

    def relative_particle_velocity_x(field,data):
        return data['deposit','all_cic_velocity_x']-data['index','bulk_velocity_x']
    ds.add_field(('deposit','relative_velocity_x'),
           sampling_type="cell",
           function=relative_particle_velocity_x,
           units='cm/s',force_override = True)
    def relative_particle_velocity_y(field,data):
        return data['deposit','all_cic_velocity_y']-data['index','bulk_velocity_y']
    ds.add_field(('deposit','relative_velocity_y'),
           sampling_type="cell",
           function=relative_particle_velocity_y,
           units='cm/s',force_override = True)
    def relative_particle_velocity_z(field,data):
        return data['deposit','all_cic_velocity_z']-data['index','bulk_velocity_z']
    ds.add_field(('deposit','relative_velocity_z'),
           sampling_type="cell",
           function=relative_particle_velocity_z,
           units='cm/s',force_override = True)

    def relative_momentum_x(field,data):
        return data['gas','relative_velocity_x']*data['gas','mass']+data['deposit','relative_velocity_x']*data['deposit','all_cic_mass']
    ds.add_field(('deposit','relative_momentum_x'),
           sampling_type="cell",
           function=relative_momentum_x,
           units='cm/s',force_override = True)
    def relative_momentum_y(field,data):
        return data['gas','relative_velocity_y']*data['gas','mass']+data['deposit','relative_velocity_y']*data['deposit','all_cic_mass']
    ds.add_field(('deposit','relative_momentum_y'),
           sampling_type="cell",
           function=relative_momentum_y,
           units='cm/s',force_override = True)
    def relative_momentum_z(field,data):
        return data['gas','relative_velocity_z']*data['gas','mass']+data['deposit','relative_velocity_z']*data['deposit','all_cic_mass']
    ds.add_field(('deposit','relative_momentum_z'),
           sampling_type="cell",
           function=relative_momentum_z,
           units='cm/s',force_override = True)

    def angular_momentum_x(field,data):
        return data['index','radial_distance_y']*data[('deposit', 'relative_momentum_z')] - data['index','radial_distance_z']*data[('deposit', 'relative_momentum_y')]
    ds.add_field(('deposit','angular_momentum_x'),
               sampling_type="cell",
               function=angular_momentum_x,
               units='cm**2*g/s',force_override = True)
    def angular_momentum_y(field,data):
        return data['index','radial_distance_z']*data[('deposit', 'relative_momentum_x')] - data['index','radial_distance_x']*data[('deposit', 'relative_momentum_z')]
    ds.add_field(('deposit','angular_momentum_y'),
               sampling_type="cell",
               function=angular_momentum_y,
               units='cm**2*g/s',force_override = True)
    def angular_momentum_z(field,data):
        return data['index','radial_distance_x']*data[('deposit', 'relative_momentum_y')] - data['index','radial_distance_y']*data[('deposit', 'relative_momentum_x')]
    ds.add_field(('deposit','angular_momentum_z'),
               sampling_type="cell",
               function=angular_momentum_z,
               units='cm**2*g/s',force_override = True)

def calculate_sfr(ds,center,Rvir):
    # I have no idea how to do this
    return np.nan

def get_required_quantities(ds,center,Rvir,Mvir):
    sphere_Rvir = ds.sphere(center,(Rvir,'kpc'))
    sphere_gal = ds.sphere(center,(0.1*Rvir,'kpc'))
    ad = ds.all_data()
    z = ds.current_redshift
    a = 1./(z+1.)
    dict_of_quantities = {'a':a,'center_x':center[0],'center_y':center[1],'center_z':center[2],'Rvir':Rvir,'Mvir':Mvir}
    Mgal = sphere_gal.quantities.total_quantity(["complete_mass"])
    Mstar = sphere_gal.quantities.total_quantity([('stars','particle_mass')])
    dict_of_quantities['Mstar'] = Mstar
    Mgas = sphere_Rvir.quantities.total_quantity([('gas','mass')])
    dict_of_quantities['Mgas'] = Mgas
    Mdm = sphere_Rvir.quantities.total_quantity([('darkmatter','particle_mass')])
    dict_of_quantities['Mdm'] = Mdm
    L_x = sphere_gal.quantities.total_quantity([('deposit','angular_momentum_x')])
    L_y = sphere_gal.quantities.total_quantity([('deposit','angular_momentum_y')])
    L_z = sphere_gal.quantities.total_quantity([('deposit','angular_momentum_z')])
    Lmag = np.sqrt(L_x**2+L_y**2+L_z**2)
    dict_of_quantities['Lmag'] = Lmag
    L_x,L_y,L_z = L_x/Lmag,L_y/Lmag,L_z/Lmag
    dict_of_quantities['L_x'] = L_x
    dict_of_quantities['L_y'] = L_y
    dict_of_quantities['L_z'] = L_z
    bulk_velocity_x = ad[('index','bulk_velocity_x')][0]
    bulk_velocity_y = ad[('index','bulk_velocity_y')][0]
    bulk_velocity_z = ad[('index','bulk_velocity_z')][0]
    dict_of_quantities['bulk_velocity_x'] = bulk_velocity_x
    dict_of_quantities['bulk_velocity_y'] = bulk_velocity_y
    dict_of_quantities['bulk_velocity_z'] = bulk_velocity_z
    sfr = calculate_sfr(ds,center,Rvir)
    dict_of_quantities['sfr'] = sfr
    return dict_of_quantities

def write_quantities_to_file(fullname,dict_of_quantities):
    pathname = './quasarscan_data/galaxy_catalogs/%s'%fullname
    a = dict_of_quantities['a']
    try:
        dict_of_all_info = parse_metadata.dict_of_all_info(fullname)
        all_avals = parse_metadata.all_avals(fullname,dict_of_all_info)
        index = len(all_avals[a<all_avals])
        if a not in all_avals:
            all_avals = np.array(list(all_avals[:index])+[a]+list(all_avals[index:]))
        path_exists = True
    except NoMetadataError:
        dict_of_all_info = {}
        index = 0
        path_exists = False

    all_keys = ['a']
    new_keys = dict_of_quantities.keys()
    existing_keys = dict_of_all_info.keys()
    for key in new_keys+existing_keys:
        if key not in all_keys:
            all_keys.append(key)
    all_lines = [all_keys]
    for aval in all_avals:
        current_line = ''
        if aval == a:
            for quantity in all_keys:
                try:
                    current_line+='%s, '%dict_of_quantities[quantity]
                except KeyError:
                    current_line+='np.nan, '
        else:
            for quantity in all_keys:
                try:
                    current_line+='%s, '%dict_of_all_info[quantity][aval]
                except KeyError:
                    current_line+='np.nan, '
                    continue
        all_lines.append(current_line[:-2])

    with open(pathname,'w') as f:
        for line in all_lines:
            f.write(line)


def create_metadata_table(fullname,filepath,hc = None,ith_largest = 1):
    code = fullname.split('_')[2]
    ds = code_specific_setup.ytload(filepath,code)
    center = get_halo_center(code,ds,hc=None,ith_largest = 1)
    Rvir,Mvir = find_virial_radius(ds,center)
    add_all_fields(ds,center,Rvir,Mvir)
    dict_of_quantities = get_required_quantities(ds,center,Rvir)
    write_quantities_to_file(fullname,dict_of_quantities)







