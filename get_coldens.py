import yt
import numpy as np
import trident
import quasar_scan
import sys
import os
from quasar_scan import tprint

yt.funcs.mylog.setLevel(50)

yt.enable_parallelism()

filename = quasar_scan.read_command_line_args(sys.argv, "-fn","--filename", 1, ["None"])[0]
simname, redshift = quasar_scan.read_command_line_args(sys.argv, "-sz","--simnameredshift", 2, ["None",-1.0])
if filename == "None" and simname == "None":
    tprint("no file to load")
elif filename == "None":
    filename = quasar_scan.get_filename_from_simname(simname, redshift)
simparams,scanparams,ions,data = quasar_scan.read_values(filename)

test = False
save = quasar_scan.read_command_line_args(sys.argv, "-s","--save", 1, [12])[0]

q = quasar_scan.QuasarSphere(simparams=simparams,scanparams=scanparams,ions=ions,data=data)
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
        
        ray = trident.make_simple_ray(q.ds,
            start_position=start,
            end_position=end,
            data_filename="ray"+ident+".h5",
            fields = [('gas',"metallicity")],
            ftype='gas')
        trident.add_ion_fields(ray,ions)
        field_data = ray.all_data()
        for j in range(len(ions)):
            ion = ions[j]
            cdens = np.sum(field_data[("gas",quasar_scan.ion_to_field_name(ion))] * field_data['dl'])
            #outcdens = np.sum((field_data['radial_velocity']>0)*field_data[ion_to_field_name(ion)]*field_data['dl'])
            #incdens = np.sum((field_data['radial_velocity']<0)*field_data[ion_to_field_name(ion)]*field_data['dl'])
            vector[11+j] = cdens
            #vector[12+3*i+1] = outcdens
            #vector[12+3*i+2] = incdens
        Z = np.average(field_data[('gas',"metallicity")],weights=field_data['dl'])
        vector[-1] = Z
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
            output = q.save_values()
            tprint("file saved to "+output+".")
