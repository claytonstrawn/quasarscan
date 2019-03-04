import yt
import numpy as np
import trident
try:
    from quasarscan import quasar_sphere
except:
    import quasar_sphere
import sys
import os
import datetime
from yt.utilities.physical_constants import mh

use_tprint = True
def tprint(*args,**kwargs):
    if use_tprint:
        print(datetime.datetime.now())
    print(args)

yt.funcs.mylog.setLevel(50)

def get_cmd_args():
    filename = sys.argv[1]
    save = int(sys.argv[2])
    parallel = sys.argv[3]=='p'
    return filename,save,parallel


atoms = ['C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', \
        'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', \
        'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']
fields_to_keep_art = [('gas',"H_nuclei_density"),('gas',"metal_density"),('gas',"density"),('gas',"temperature"),('gas',"radial_velocity"),('gas', 'cell_volume')]
fields_to_keep_ramses = [('gas',"H_nuclei_density"),('gas',"metal_density"),('gas',"density"),('gas',"temperature"),('gas',"radial_velocity"),('gas', 'cell_volume')]
fields_to_keep_gizmo = [('gas',"H_nuclei_density"),('gas',"metal_density"),('gas',"density"),('gas',"temperature"),('gas',"radial_velocity"),('gas', 'cell_volume')]
fields_to_keep_gadget = [('gas',"H_nuclei_density"),('gas',"metal_density"),('gas',"density"),('gas',"temperature"),('gas',"radial_velocity"),('gas', 'cell_volume')]
fields_to_keep_gear = [('gas',"H_nuclei_density"),('gas',"metal_density"),('gas',"density"),('gas',"temperature"),('gas',"radial_velocity"),('gas', 'cell_volume')]
fields_to_keep_enzo = [('gas',"H_nuclei_density"),('gas',"metal_density"),('gas',"density"),('gas',"temperature"),('gas',"radial_velocity"),('gas', 'cell_volume')]
fields_to_keep_dict = {'art':fields_to_keep_art,'ramses':fields_to_keep_ramses,'gizmo':fields_to_keep_gizmo,'gadget':fields_to_keep_gadget,'gear':fields_to_keep_gear,'enzo':fields_to_keep_enzo}

def atoms_from_ions(ions):
    toret = []
    for ion in ions:
        toret.append(ion.split(" ")[0])
    return list(set(toret))

filename,save,parallel = get_cmd_args()
if parallel:
    yt.enable_parallelism()
try:
    readvalsoutput = quasar_sphere.read_values(filename)
except:
    readvalsoutput = quasar_sphere.read_values(filename.split('quasarscan/')[1])
test = False
q = quasar_sphere.QuasarSphere(readvalsoutput=readvalsoutput)
fields_to_keep = fields_to_keep_dict[q.code]

g = q.gasbins
for atom in atoms:
    if q.code == 'art' and atom in atoms_from_ions(q.ions):
        fields_to_keep.append(('gas','%s_nuclei_mass_density'%atom))
if q.code == 'art':
    h,d,s = quasar_sphere.get_aux_files_art(q.dspath)
    ds = yt.load(q.dspath,file_particle_header=h,\
                                  file_particle_data=d,\
                                  file_particle_stars=s)
else:
    ds = yt.load(q.dspath)
if ('deposit','Gas_mass') in ds.derived_field_list:
    def gas_mass(field, data):
        return data['deposit','Gas_mass']
    ds.add_field(('gas','mass'),units = 'g', function = gas_mass, sampling_type = 'cell')
num_bin_vars = q.gasbins.get_length()
starting_point = q.length_reached 
bins = np.append(np.arange(starting_point,q.length,save)[:-1],q.length)
for i in range(0, len(bins)-1):
    current_info = q.info[bins[i]:bins[i+1]]
    if yt.is_root():
        tprint("%s-%s /%s"%(bins[i],bins[i+1],len(q.info)))
    my_storage = {}
    for sto,in_vec in yt.parallel_objects(current_info, storage = my_storage):
        vector = np.copy(in_vec)
        index = vector[0]
        tprint("<line %d, starting process> "%index)
        ident = str(index)
        start = vector[5:8]
        end = vector[8:11]
        ray = trident.make_simple_ray(ds,
            start_position=start,
            end_position=end,
            data_filename="ray"+ident+".h5",
            fields = fields_to_keep,
            ftype='gas')
        trident.add_ion_fields(ray,q.ions)
        field_data = ray.all_data()
        for j in range(len(q.ions)):
            ion = q.ions[j]
            dl = field_data['dl']
            ionfield = field_data[("gas",quasar_sphere.ion_to_field_name(ion))]
            cdens = np.sum(ionfield * dl)
            vector[11+j*(num_bin_vars+2)] = cdens
            atom = ion.split(" ")[0]
            total_nucleus = np.sum(ionfield[ionfield>0]/\
                                       field_data[("gas",quasar_sphere.ion_to_field_name(ion).replace("number_density","ion_fraction"))][ionfield>0]\
                                        * dl[ionfield>0])
            vector[11+j*(num_bin_vars+2)+1] = cdens / total_nucleus
            for k in range(num_bin_vars):
                try:
                    variable_name,edges,units = q.gasbins.get_field_binedges_for_num(k)
                    if units:
                        data = field_data[variable_name].in_units(units)
                    else:
                        data = field_data[variable_name]
                    abovelowerbound = data>edges[0]
                    belowupperbound = data<edges[1]
                    withinbounds = np.logical_and(abovelowerbound,belowupperbound)
                    coldens_in_line = (ionfield[withinbounds])*(dl[withinbounds])
                    coldens_in_bin = np.sum(coldens_in_line)
                    vector[11+j*(num_bin_vars+2)+k+2] = coldens_in_bin/cdens
                except Exception as e:
                    print "Could not bin into %s with edges %s because of error %s"%(variable_name,edges,e)
        try:
            Z = np.sum(field_data[('gas',"metal_density")]*dl)/ \
                np.sum(field_data[('gas',"H_nuclei_density")]*mh*dl)
            vector[-1] = Z
        except Exception as e:
            print "Could not get average metallicity because of error %s"%(e)
        try:
            n = np.average(field_data['density'],weights=field_data['density']*dl)
            vector[-2] = n
        except Exception as e:
            print "Could not get average density because of error %s"%(e)
        try:
            T = np.average(field_data['temperature'],weights=field_data['density']*dl)
            vector[-3] = T
        except Exception as e:
            print "Could not get average temperature because of error %s"%(e)
        try:
            os.remove("ray"+ident+".h5")
        except:
            pass 
        tprint("vector = "+str(vector))
        sto.result_id = index
        sto.result = vector
    if yt.is_root():
        keys = my_storage.keys()
        for key in keys:
            q.info[int(key)] = my_storage[key]
        q.scanparams[6]+=(bins[i+1]-bins[i])
        q.length_reached = q.scanparams[6]
        if not test:
            output = q.save_values(at_level = 0)
            tprint("file saved to "+output+".")
