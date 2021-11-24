from yt.utilities.cosmology import Cosmology
from yt.data_objects.particle_filters import add_particle_filter
from yt.utilities.physical_constants import mh
import numpy as np 

def set_up_art(ds):
    def H_nuclei_density(field,data):
        return data['gas','density']*data['gas','H_mass_fraction']/mh
    ds.add_field(('gas','H_nuclei_density'),
                sampling_type = 'cell',
                function = H_nuclei_density,
                units = 'cm**-3')
    
def set_up_enzo(ds):
    def star_filter(pfilter, data):
        return (data[pfilter.filtered_type, 'creation_time'] > 0)
    add_particle_filter('stars', function=star_filter, filtered_type='all', requires=['creation_time'])
    ds.add_particle_filter('stars')
    def dm_filter(pfilter, data):
        return (data[pfilter.filtered_type, 'creation_time'] == 0)
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
    
def set_up_gadget(ds):
    def star_mass_rename(field,data):
        return data['PartType4','particle_mass']
    ds.add_field(('stars','particle_mass'),
                sampling_type = 'particle',
                function = star_mass_rename,
                units = 'g')
    
    ad = ds.all_data()
    m1 = ad['PartType1','Masses']
    assert np.amin(m1)==np.amax(m1)
    unique_dm1_mass = np.amin(m1)
    m5 = ad['PartType5','Masses']
    unique_dm5_masses = np.unique(m5)
    def dm_filter(pfilter, data):
        allowed = data[pfilter.filtered_type, 'Masses']==unique_dm1_mass
        for v in unique_dm5_masses:
            allowed = np.logical_or(allowed,data[pfilter.filtered_type, 'Masses']==v)
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
        abovemin = data[pfilter.filtered_type, 'Masses']>=low
        belowmax = data[pfilter.filtered_type, 'Masses']<=high
        allowed = np.logical_and(abovemin,belowmax)
        for v in unique_dm5_masses:
            allowed = np.logical_or(allowed,data[pfilter.filtered_type, 'Masses']==v)
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

    Z_Solar= 0.02041
    def new_metallicity(field, data):
        if len(data['PartType0', "Metals"].shape) == 1:
            return data['PartType0', "Metals"].in_units("")/Z_Solar
        else:
            return data['PartType0', "Metals"][:,9].in_units("")/Z_Solar
    ds.add_field(('gas','metallicity'),
                sampling_type = 'cell',
                function = new_metallicity,
                units = '',
                force_override = True)
    def metal_density(field,data):
        return data['gas','metallicity']*data['gas','density']
    ds.add_field(('gas','metal_density'),
                sampling_type = 'cell',
                function = metal_density,
                units = 'g/cm**3')

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
        allowed = np.zeros(data[pfilter.filtered_type, 'Masses'].shape,dtype=bool)
        for v in ds.arr(np.concatenate([unique_dm1_masses,unique_dm2_masses]).v,m2.units):
            allowed = np.logical_or(allowed,data[pfilter.filtered_type, 'Masses']==v)
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
    
def set_up_mockstreams(ds):
    def metal_density(field,data):
        return data['gas','metallicity']*data['gas','density']
    ds.add_field(('gas','metal_density'),
                sampling_type = 'cell',
                function = metal_density,
                units = 'g/cm**3')
    
set_up_funcs = {'art':set_up_art,'enzo':set_up_enzo,'ramses':set_up_ramses,'changa':set_up_changa,'gear':set_up_gear,'gizmo':set_up_gizmo,'gadget':set_up_gadget,'mockstreams':set_up_mockstreams}    

def add_radial_distance_fields(ds,center):
    for i,ax in enumerate('xyz'):
        def radial_distance_ax(field,data):
            return data['gas',ax]-center[i]
        ds.add_field(('gas','radial_distance_%s'%ax),
                   sampling_type='cell',
                   function=radial_distance_ax,
                   units='cm',force_override = True)
    
def add_relative_velocity_fields(ds,v):
    for i,ax in enumerate('xyz'):
        def relative_velocity_ax(field,data):
            return data['gas','velocity_%s'%ax]-v[i]
        ds.add_field(('gas','relative_velocity_%s'%ax),
                   sampling_type='cell',
                   function=relative_velocity_ax,
                   units='cm/s',force_override = True)
        def relative_momentum_ax(field,data):
            return data['gas','relative_velocity_%s'%ax]*data['gas','mass']
        ds.add_field(('gas','relative_momentum_%s'%ax),
                   sampling_type='cell',
                   function=relative_momentum_ax,
                   units='g*cm/s',force_override = True)
        
def add_radial_velocity_fields(ds,v):
    def radial_velocity(field,data):
        tot_distance = np.sqrt(data['gas','radial_distance_x']**2+data['gas','radial_distance_y']**2+data['gas','radial_distance_z']**2)
        x_comp = data['gas','radial_distance_x']*data['gas','relative_velocity_x']/tot_distance
        y_comp = data['gas','radial_distance_y']*data['gas','relative_velocity_y']/tot_distance
        z_comp = data['gas','radial_distance_z']*data['gas','relative_velocity_z']/tot_distance
        return x_comp+y_comp+z_comp
    ds.add_field(('gas','radial_velocity'),
               sampling_type='cell',
               function=radial_velocity,
               units='km/s',force_override = True)