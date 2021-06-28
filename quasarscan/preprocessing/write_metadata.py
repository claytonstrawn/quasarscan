from yt_astro_analysis.halo_analysis import HaloCatalog
from yt.utilities.cosmology import Cosmology
from yt import YTArray
from quasarscan.utils.utils import sphcodes
from quasarscan.preprocessing import code_specific_setup,parse_metadata,add_common_fields
from quasarscan import __version__

import numpy as np
import os
import datetime

class NotEnoughHalosError(Exception):
    def __init__(self, message):
        self.message = message
        print(self.message)

        
class UnitNotUnderstoodError(Exception):
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
        best_resolution = np.min(ad['index','cell_volume'])
    else:
        best_resolution = None
    sorted_halos = sorted(hc.catalog,key=lambda x:x['particle_mass'])
    catalog_of_acceptable_halos = []
    for i,h in enumerate(sorted_halos):
        if len(catalog_of_acceptable_halos)==ith_largest:
            return ds.arr(catalog_of_acceptable_halos[-1],'code_length')
        center = np.array([h['particle_position_x'].in_units('unitary').v,
                            h['particle_position_y'].in_units('unitary').v,
                            h['particle_position_z'].in_units('unitary').v])
        approx_Rvir = h['virial_radius'].in_units('kpc')*1.126
        approx_halo = ds.sphere(center,approx_Rvir)
        if best_resolution is not None and np.min(approx_halo['index','cell_volume']) == best_resolution:
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
    Rvir = np.interp(200*rho_c,np.flip(densities),np.flip(rs))
    Mvir = cell_mass+particle_mass
    return ds.arr(Rvir,'kpc'),ds.arr(Mvir,'g')

def calculate_bulk_velocity(ds,sphere_Rvir):
    return sphere_Rvir.quantities.bulk_velocity()

def calculate_ang_mom(ds,sphere_gal,v):
    L = ds.arr([0.,0.,0.],'cm**2*g/s')
    cp_ax_num = {0:(1,2),1:(2,0),2:(0,1)}
    num2str = {0:'x',1:'y',2:'z'}
    for i,ax in enumerate('xyz'):
        ax0 = num2str[cp_ax_num[i][0]]
        ax1 = num2str[cp_ax_num[i][1]]
        term_1 = sphere_gal['gas','radial_distance_%s'%ax0]*(sphere_gal['gas','velocity_%s'%ax1]-v[cp_ax_num[i][1]])*sphere_gal['gas','mass']
        term_2 = sphere_gal['gas','radial_distance_%s'%ax1]*(sphere_gal['gas','velocity_%s'%ax0]-v[cp_ax_num[i][0]])*sphere_gal['gas','mass']
        L[i] = np.sum(term_1 - term_2)
    Lmag = np.sqrt(np.sum(L**2))
    L_x,L_y,L_z = L[0]/Lmag,L[1]/Lmag,L[2]/Lmag
    return Lmag,(L_x,L_y,L_z)

def calculate_sfr(ds,sphere_star,throw_errors = False):
    masses = sphere_star['stars','particle_mass']
    if len(masses) == 0:
        print('no stars in sphere!')
        return ds.arr(0,'Msun/yr')
    creation_times = sphere_star['stars','particle_creation_time']
    last_5_times = np.unique(creation_times)[-5:]
    first = last_5_times[0]
    mass = np.sum(masses[creation_times>=first])
    time = ds.current_time - first
    return (mass/time).in_units('Msun/yr')


def get_required_quantities(ds,center,Rvir,Mvir,stars_boundary,gal_edge):
    z = ds.current_redshift
    a = 1./(z+1.)
    dict_of_quantities = {'a':a,'center_x':center[0],'center_y':center[1],'center_z':center[2],'Rvir':Rvir,'Mvir':Mvir}

    sphere_Rvir = ds.sphere(center,Rvir)
    bulk_velocity = calculate_bulk_velocity(ds,sphere_Rvir)
    dict_of_quantities['bulk_velocity_x'] = bulk_velocity[0]
    dict_of_quantities['bulk_velocity_y'] = bulk_velocity[1]
    dict_of_quantities['bulk_velocity_z'] = bulk_velocity[2]
        
    Mgas = np.sum(sphere_Rvir['gas','mass'])
    dict_of_quantities['Mgas'] = Mgas
    Mdm = np.sum(sphere_Rvir['darkmatter','particle_mass'])
    dict_of_quantities['Mdm'] = Mdm
    
    sphere_gal = ds.sphere(center,gal_edge*Rvir)
    Lmag,L_norm = calculate_ang_mom(ds,sphere_gal,bulk_velocity)
    dict_of_quantities['Lmag'] = Lmag
    dict_of_quantities['L_x'] = L_norm[0]
    dict_of_quantities['L_y'] = L_norm[1]
    dict_of_quantities['L_z'] = L_norm[2]

    sphere_star = ds.sphere(center,stars_boundary*Rvir)
    Mstar = np.sum(sphere_star['stars','particle_mass'])
    dict_of_quantities['Mstar'] = Mstar
    sfr = calculate_sfr(ds,sphere_star)
    dict_of_quantities['sfr'] = sfr
    return dict_of_quantities


unit_conversion = {'a':'', 'center_x':'unitary', 'center_y':'unitary','center_z':'unitary', 'Rvir':'kpc', 'Mvir':'Msun', 'Mstar':'Msun', 'Mgas':'Msun', 'Mdm':'Msun', 'Lmag':'cm**2*g/s', 'L_x':'', 'L_y':'', 'L_z':'', 'bulk_velocity_x':'km/s', 'bulk_velocity_y':'km/s', 'bulk_velocity_z':'km/s', 'sfr':'Msun/yr'}
    
def write_quantities_to_file(fullname,dspath,dict_of_quantities,tolerance = .001):
    pathname = os.path.expanduser('~/quasarscan_data/galaxy_catalogs/%s'%fullname)
    if not os.path.exists(pathname):
        os.mkdir(pathname)
    a = dict_of_quantities['a']
    try:
        dict_of_all_info = parse_metadata.dict_of_all_info(fullname)
        all_avals = parse_metadata.all_avals(fullname,dict_of_all_info)
        if a not in all_avals:
            all_avals = sorted([a]+list(all_avals))
    except parse_metadata.NoMetadataError:
        dict_of_all_info = {}
        all_avals = [a]

    dict_of_quantities['a'] = YTArray(a,'')
    all_keys = ['a']
    new_keys = list(dict_of_quantities.keys())
    existing_keys = list(dict_of_all_info.keys())
    for key in new_keys+existing_keys:
        if key not in all_keys:
            all_keys.append(key)
    all_lines = ["Metadata recorded on file %s with quasarscan version %s on date %s"%\
                 (dspath,__version__,str(datetime.datetime.now()))]
    all_lines.append(str(all_keys)[1:-1].replace("'",""))
    for aval in all_avals:
        current_line = ''

        for quantity in all_keys:
            try:
                if np.abs(aval-a)<tolerance:
                    to_write = dict_of_quantities[quantity]
                    try:
                        convert_to = unit_conversion[quantity]
                        to_write = float(to_write.in_units(convert_to).v)
                    except:
                        raise UnitNotUnderstoodError(quantity)
                else:
                    to_write = dict_of_all_info[quantity][aval]
                current_line += '%s, '%to_write
            except KeyError:
                current_line+='np.nan, '
                continue
        all_lines.append(current_line[:-2])

    with open(os.path.join(pathname,fullname+'_metadata.txt'),'w') as f:
        for line in all_lines:
            f.write(line+'\n')


def create_metadata_table(fullname,filepath,hc = None,ith_largest = 1,Rvir=None,
                          center=None,stars_boundary = 0.1, gal_edge = 0.1):
    code = fullname.split('_')[2]
    ds,_ = code_specific_setup.load_and_setup(filepath,code)
    if center is None:
        center = get_halo_center(code,ds,hc=hc,ith_largest = 1)
    else:
        center = ds.arr(center,'unitary')
    print(center)
    if Rvir is None:
        Rvir,Mvir = find_virial_radius(ds,center)
    else:
        Rvir = ds.arr(Rvir,'kpc')
        sp = ds.sphere(center,Rvir)
        cell_mass,particle_mass = sp.quantities.total_mass()
        Mvir = cell_mass+particle_mass
    add_common_fields.add_radial_distance_fields(ds,center)
    dict_of_quantities = get_required_quantities(ds,center,Rvir,Mvir,stars_boundary,gal_edge)
    dspath = ds.fullpath
    write_quantities_to_file(fullname,dspath,dict_of_quantities)
