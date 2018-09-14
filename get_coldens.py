import yt
import numpy as np
import trident
import quasar_scan
import sys
import os
from quasar_scan import tprint
import gasbinning
from yt.utilities.physical_constants import mh

yt.funcs.mylog.setLevel(50)

def get_filename_and_save():
    filename = quasar_scan.read_command_line_args(sys.argv, "-fn","--filename", 1, ["None"])[0]
    simname, redshift = quasar_scan.read_command_line_args(sys.argv, "-sz","--simnameredshift", 2, ["None",-1.0])
    if filename == "None" and simname == "None":
        tprint("no file to load")
    elif filename == "None":
        filename = quasar_scan.get_filename_from_simname(simname, redshift)
    save = quasar_scan.read_command_line_args(sys.argv, "-s","--save", 1, [12])[0]
    print filename
    return filename,save

atoms = ['C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', \
        'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', \
        'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']
fields_to_keep = [('gas',"H_nuclei_density"),('gas',"metal_density"),('gas',"density"),('gas',"temperature"),('gas',"radial_velocity"),('gas', 'dx')]
for atom in atoms:
    fields_to_keep.append(('gas','%s_nuclei_mass_density'%atom))

yt.enable_parallelism()
filename,save = get_filename_and_save()
simparams,scanparams,ions,data,gasbins = quasar_scan.read_values(filename)
test = False
q = quasar_scan.QuasarSphere(simparams=simparams,scanparams=scanparams,ions=ions,data=data,gasbins = gasbins)
num_bin_vars = gasbins.get_length()
starting_point = q.length_reached 
bins = np.append(np.arange(starting_point,q.length,save)[:-1],q.length)
try:
    os.remove("ray0.0.h5")
except:
    pass 
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
        ray = trident.make_simple_ray(q.ds,
            start_position=start,
            end_position=end,
            data_filename="ray"+ident+".h5",
            fields = fields_to_keep,
            ftype='gas')
        trident.add_ion_fields(ray,ions)
        field_data = ray.all_data()
        for j in range(len(ions)):
            ion = ions[j]
            cdens = np.sum(field_data[("gas",quasar_scan.ion_to_field_name(ion))] * field_data['dl'])
            vector[11+j*(num_bin_vars+2)] = cdens
            atom = ion.split(" ")[0]
            total_nucleus = np.sum(field_data[("gas","%s_nuclei_mass_density"%atom)]/(mh*trident.ion_balance.atomic_mass[atom]) * field_data['dl'])
            vector[11+j*(num_bin_vars+2)+1] = cdens / total_nucleus
            for k in range(num_bin_vars):
                variable_name,edges,units = gasbins.get_field_binedges_for_num(k)
                if units:
                    data = field_data[variable_name].in_units(units)
                else:
                    data = field_data[variable_name]
                abovelowerbound = data>edges[0]
                belowupperbound = data<edges[1]
                withinbounds = np.logical_and(abovelowerbound,belowupperbound)
                coldens_in_line = (field_data[("gas",quasar_scan.ion_to_field_name(ion))][withinbounds])*(field_data['dl'][withinbounds])
                coldens_in_bin = np.sum(coldens_in_line)
                vector[11+j*(num_bin_vars+2)+k+2] = coldens_in_bin/cdens
        Z = np.sum(field_data[('gas',"metal_density")]*field_data['dl'])/ \
            np.sum(field_data[('gas',"H_nuclei_mass_density")]*field_data['dl'])
        vector[-1] = Z
        n = np.average(field_data['density'],weights=field_data['density']*field_data['dl'])
        vector[-2] = n
        T = np.average(field_data['temperature'],weights=field_data['density']*field_data['dl'])
        vector[-3] = T
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
