import yt
from yt.utilities.physical_constants import mh
import numpy as np
import trident
from quasarscan.data_objects import simulation_quasar_sphere
from quasarscan.preprocessing import code_specific_setup
from quasarscan.utils import roman
import sys
import os
import datetime
import time

class NoSimulationError(Exception):
    def __init__(self, message):
        self.message = message
        print(self.message)
    
class IllegalSightlineError(Exception):
    def __init__(self, message):
        self.message = message
        print(self.message)

#can print out the time to see how fast the script is running
use_tprint = True
def tprint(*args,**kwargs):
    if use_tprint is True:
        print(args,end=str(datetime.datetime.now())+'\n')
    elif use_tprint == 'no_time':
        print(args)

def ion_to_field_name(ion,field_type = "number_density"):
    atom = ion.split(" ")[0]
    ionization = roman.from_roman(ion.split(" ")[1])-1
    return "%s_p%s_%s"%(atom,ionization,field_type)

def throw_errors_if_allowed(e,throwerrors,message=None):
    if throwerrors == True:
        print(message)
        raise e
    elif throwerrors == 'warn':
        print(message)
        print(e)
    elif throwerrors == False:
        return

def run_sightlines(outputfilename,save_after_num,parallel,simulation_dest = None,run = 'default',throwerrors = 'warn'):
    if run not in ['default','test']:
        print('unknown option for "run" %s. Please restart with "run = default" or "run = test".'%run)
    #do not print out anything from yt (it prints plenty)
    yt.funcs.mylog.setLevel(50)
    if parallel:
        yt.enable_parallelism()
    readvalsoutput = simulation_quasar_sphere.read_values(outputfilename)
    #by creating a QuasarSphere, it knows all its metadata and other
    #information from simparams and scanparams (first lines of file at 
    #'filename')
    q = simulation_quasar_sphere.SimQuasarSphere(start_up_info_packet=readvalsoutput)
    if q.simparams[6] is None:
        if simulation_dest:
            q.simparams[6] = simulation_dest
        else:
            raise NoSimulationError('Simulation file location unknown, run with "simulation_dest" to process')
    else:
        simulation_dest = q.simparams[6]
    ds,fields_to_keep = code_specific_setup.load_and_setup(simulation_dest,q.code,q.ions)
    code_specific_setup.check_redshift(ds,outputfilename=outputfilename)
    num_bin_vars = q.gasbins.get_length()
    #Can start at a position further than 0 if reached
    starting_point = q.length_reached 
    bins = np.append(np.arange(starting_point,q.length,save_after_num)[:-1],q.length)
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
            start = ds.arr(tuple(vector[5:8]),'unitary')
            end = ds.arr(tuple(vector[8:11]),'unitary')
            try:
                ray = trident.make_simple_ray(ds,
                                            start_position=start,
                                            end_position=end,
                                            data_filename="ray"+ident+".h5",
                                            fields = fields_to_keep,
                                            ftype='gas')
            except KeyboardInterrupt:
                print('skipping sightline %s ...'%index)
                print('Interrupt again within 5 seconds to *actually* end')
                time.sleep(5)
                continue
            except Exception as e:
                throw_errors_if_allowed(e,throwerrors,'problem with making ray')
                continue
            trident.add_ion_fields(ray,q.ions)
            field_data = ray.all_data()
            dl = field_data['gas','dl']
            #3rd for loop is for processing each piece of info about each ion
            #including how much that ion is in each bin according to gasbinning
            #here just process topline data (column densities and ion fractions)
            #(~10 ions)
            for j in range(len(q.ions)):
                ion = q.ions[j]
                ionfield = field_data["gas",ion_to_field_name(ion)]
                cdens = np.sum((ionfield * dl).in_units('cm**-2')).value
                vector[11+j*(num_bin_vars+2)] = cdens
                total_nucleus = np.sum(ionfield[ionfield>0]/\
                                           field_data["gas",ion_to_field_name(ion,'ion_fraction')][ionfield>0]\
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
                        throw_errors_if_allowed(e,throwerrors,'Could not bin into %s with edges %s'%(variable_name,edges))
                toprint+="%s:%e "%(ion,cdens)
            #gets some more information from the general sightline.
            #metallicity, average density (over the whole sightline)
            #mass-weighted temperature
            try:
                if ('gas',"H_nuclei_density") in ray.derived_field_list:
                    Z = np.sum(field_data['gas',"metal_density"]*dl)/ \
                        np.sum(field_data['gas',"H_nuclei_density"]*mh*dl)
                else:
                    Z = np.sum(field_data['gas',"metal_density"]*dl)/ \
                        np.sum(field_data['gas',"number_density"]*mh*dl)
                vector[-1] = Z
            except Exception as e:
                throw_errors_if_allowed(e,throwerrors,'problem with average metallicity')
            try:
                n = np.average(field_data['gas','number_density'],weights=field_data['gas','number_density']*dl)
                vector[-2] = n
            except Exception as e:
                throw_errors_if_allowed(e,throwerrors,'problem with average density')
            try:
                T = np.average(field_data['gas','temperature'],\
                               weights=field_data['gas','density']*dl)
                vector[-3] = T
            except Exception as e:
                throw_errors_if_allowed(e,throwerrors,'problem with average temperature')
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
            if run != 'test':
                outputfilename = q.save_values(oldfilename=outputfilename)
                tprint("file saved to "+outputfilename+".")
