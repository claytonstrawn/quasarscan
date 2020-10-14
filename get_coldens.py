#this is the main "processing" file of quasarscan. It runs
#as a script, following running 'create_qso_endpoints.py'.
#loading those endpoints, it then fills in the data one line at 
#a time (but parallelizable)

import yt
from yt.utilities.physical_constants import mh
import numpy as np
import trident
try:
    from quasarscan import quasar_sphere
    from quasarscan import code_specific_setup
    level = 0
except:
    import quasar_sphere
    import code_specific_setup
    level = 1
import sys
import os
import datetime

#can print out the time to see how fast the script is running
use_tprint = True
def tprint(*args,**kwargs):
    if use_tprint:
        print(args,end=str(datetime.datetime.now())+'\n')
    else:
        print(args)

#do not print out anything from yt
yt.funcs.mylog.setLevel(50)

#filename = output/simtype_version_code_simnum_coldensinfo.txt
#it has details of simulation location inside it
#save = how many lines to process before saving (if parallel,
#       also how many processors to use at once)
#parallel = 'p' if True. If false run on one processor, if true use 
#           'save' processors. (saving is pretty fast, so 'save' = 1 
#           would not change much)
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
except IOError:
    print("unable to read from %s, checking after 'quasarscan'"%filename)
    readvalsoutput = quasar_sphere.read_values(filename.split('quasarscan/')[1])

test = False
#by creating a QuasarSphere, it knows all its metadata and other
#information from simparams and scanparams (first lines of file at 
#'filename')
q = quasar_sphere.QuasarSphere(readvalsoutput=readvalsoutput)

ds,fields_to_keep = code_specific_setup.load_and_setup(q.dspath,q.code,q.ions,add_pi_fracs = True)
convert_unit = ds.length_unit.units
num_bin_vars = q.gasbins.get_length()
#Can start at a position further than 0 if reached
starting_point = q.length_reached 
bins = np.append(np.arange(starting_point,q.length,save)[:-1],q.length)
#first for loop is non-parallel. If 32 processors available, it will break up
#into bins of size 32 at a time for example. At end, saves data from all 32.
#this is (~12 bins) in usual circumstances
for i in range(0, len(bins)-1):
    current_info = q.info[bins[i]:bins[i+1]]
    if yt.is_root():
        tprint("%s-%s /%s"%(bins[i],bins[i+1],len(q.info)))
    my_storage = {}
    #2nd for loop is parallel. Each vector goes to a different processor, and creates 
    #a separate trident sightline (~32 sightlines [in a bin]).
    #the longest processing step is ray = trident.make_simple_ray, and it's 
    #the only step which actually takes any time (below two for loops go by fast)
    for sto,in_vec in yt.parallel_objects(current_info, storage = my_storage):
        vector = np.copy(in_vec)
        index = vector[0]
        toprint = "line %s, (r = %.0f) densities "%(str(int(index)),vector[3])
        tprint("<line %d, starting process> "%index)
        ident = str(index)
        start = yt.YTArray(np.copy(vector[5:8]),convert_unit)
        end = yt.YTArray(np.copy(vector[8:11]),convert_unit)
        
        try:
            ray = trident.make_simple_ray(ds,
            start_position=start,
            end_position=end,
            data_filename="ray"+ident+".h5",
            fields = fields_to_keep,
            ftype='gas')
        except RuntimeError as e:
            print("there was a problem in making the ray!")
            print(e)
            continue
        trident.add_ion_fields(ray,q.ions)
        field_data = ray.all_data()
        dl = field_data['dl']
        #3rd for loop is for processing each piece of info about each ion
        #including how much that ion is in each bin according to gasbinning
        #here just process topline data (column densities and ion fractions)
        #(~10 ions)
        for j in range(len(q.ions)):
            ion = q.ions[j]
            ionfield = field_data["gas",quasar_sphere.ion_to_field_name(ion)]
            cdens = np.sum(ionfield * dl)
            vector[11+j*(num_bin_vars+2)] = cdens
            total_nucleus = np.sum(ionfield[ionfield>0]/\
                                       field_data[("gas",quasar_sphere.ion_to_field_name(ion,'ion_fraction'))][ionfield>0]\
                                        * dl[ionfield>0])
            vector[11+j*(num_bin_vars+2)+1] = cdens / total_nucleus
            #4th for loop is processing each gasbin for the current ion
            #(~20 bins)
            for k in range(num_bin_vars):
                try:
                    variable_name,edges,units = q.gasbins.get_field_binedges_for_num(k, ion)
                    if variable_name is None:
                        vector[11+j*(num_bin_vars+2)+k+2] = np.nan
                    elif variable_name in ray.derived_field_list:
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
                    else:
                        print(str(variable_name)+" not in ray.derived_field_list")
                except Exception as e:
                    print("Could not bin into %s with edges %s because of error %s"%(variable_name,edges,e))
            toprint+="%s:%e "%(ion,cdens)
        #gets some more information from the general sightline.
        #metallicity, average density (over the whole sightline)
        #mass-weighted temperature
        try:
            Z = np.sum(field_data[('gas',"metal_density")]*dl)/ \
                np.sum(field_data[('gas',"H_nuclei_density")]*mh*dl)
            vector[-1] = Z
        except Exception as e:
            print("Could not get average metallicity because of error %s"%(e))
        try:
            n = np.average(field_data['density'],weights=field_data['density']*dl)
            vector[-2] = n
        except Exception as e:
            print("Could not get average density because of error %s"%(e))
        try:
            T = np.average(field_data['temperature'],weights=field_data['density']*dl)
            vector[-3] = T
        except Exception as e:
            print("Could not get average temperature because of error %s"%(e))
        try:
            os.remove("ray"+ident+".h5")
        except:
            pass 
        tprint(toprint)
        #'vector' now contains real data, not just '-1's
        sto.result_id = index
        sto.result = vector
    #save all parallel sightlines after they finish (every 32 lines are saved at once)
    if yt.is_root():
        keys = my_storage.keys()
        for key in keys:
            q.info[int(key)] = my_storage[key]
        q.scanparams[6]+=(bins[i+1]-bins[i])
        q.length_reached = q.scanparams[6]
        if not test:
            output = q.save_values(at_level = 0)
            tprint("file saved to "+output+".")
