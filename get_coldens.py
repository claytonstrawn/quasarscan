import yt
import numpy as np
import trident
try:
    from quasarscan import quasar_sphere
    from quasarscan import code_specific_setup
except:
    import quasar_sphere
    import code_specific_setup
import sys
import os
import datetime

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


filename,save,parallel = get_cmd_args()
if parallel:
    yt.enable_parallelism()
try:
    readvalsoutput = quasar_sphere.read_values(filename)
except:
    readvalsoutput = quasar_sphere.read_values(filename.split('quasarscan/')[1])

test = False
q = quasar_sphere.QuasarSphere(readvalsoutput=readvalsoutput)

ds,fields_to_keep = code_specific_setup.load_and_setup(q.dspath,q.code)
convert_unit = ds.length_unit.units
print convert_unit
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
        start = yt.YTArray(vector[5:8],convert_unit)
        end = yt.YTArray(vector[8:11],convert_unit)
        print start,end, q.center
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
