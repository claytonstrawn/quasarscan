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