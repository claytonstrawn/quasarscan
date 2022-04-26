import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import datetime
from functools import reduce

from quasarscan.preprocessing import parse_metadata
from quasarscan.utils import roman,ion_lists
from quasarscan.utils.utils import sort_ions,reversearray,split_by_ops
from quasarscan.plotting.sorter import MultiSphereSorter
from quasarscan.plotting.read_and_init_from_file import read_and_init
from quasarscan.plotting.var_labels_interpreter import decide_plot_type
from quasarscan.data_objects import gasbinning,\
                                     observation_quasar_sphere,\
                                     quasar_sphere,\
                                     simulation_quasar_sphere,\
                                     multi_quasar_sphere
from quasarscan.utils.variable_lists import stringcriteria,intensives,intensiveslabels,\
                                            sightline_xVars,param_xVars,\
                                            sightline_unit_labels,param_unit_labels,\
                                            all_known_variables

#summary: Fetches y values from saved data. First, has to get relevant info, then has 
#         to perform math formula if any
#
#inputs: gq: GeneralizedQuasarSphere object (one or more quasarSpheres squished into one)
#        stringVar: formula to show. Can use basic '+-*/()' characters. Variables have the names
#                   'O VI' or 'O VI:temperature:hot' or things like this. Look at gasbinning.py for
#                   more info
#        
#outputs: ys: unprocessed y data            
def get_yVar_from_str(gq,stringVar):
    strings_to_find = split_by_ops(stringVar)
    new_str_to_eval = stringVar
    strings_to_replace_with = {}
    for i,s in enumerate(strings_to_find):
        if len(s)>0:
            try:
                eval(s)
                string_to_replace_with = s
            except:
                string_to_replace_with = "gq.info[:,%d]"%gq.get_ion_column_num(s)
            strings_to_replace_with[s] = string_to_replace_with
    for s in sorted(strings_to_replace_with.keys(),key=len,reverse=True):
        new_str_to_eval = new_str_to_eval.replace(s,strings_to_replace_with[s])
    to_return = eval(new_str_to_eval)
    return to_return


#summary: In a type 0 plot, get sightline values for y and either discrete sightline values for x
#         or metadata values for x. (it really just calls 'get_xy_type1' or 'get_xy_type2')
#
#inputs: xVar: xVar or xVar formula being displayed
#        yVars: list of ions or formulas being displayed
#        rlims: what impact parameters to accept. (from 'get_sightline_xy_vals')
#        
#outputs: xs: raw x data for processing
#         ys: raw y data for processing 
def get_xy_type0(xVar,yVars,quasar_array,rlims):
    xs = np.empty(len(yVars),dtype = object)
    ys = np.empty(len(yVars),dtype = object)
    for i in range(len(yVars)):
        yVar = yVars[i]
        plot_type = decide_plot_type(xVar,yVar)
        if plot_type == 1:
            x,y = get_xy_type1(xVar,yVar,quasar_array,rlims)
        elif plot_type == 2:
            x,y = get_xy_type2(xVar,yVar,quasar_array,rlims)
        elif plot_type == 3:
            x,y = get_xy_type3(xVar,yVar,quasar_array,rlims)
        elif plot_type == 4:
            x,y = get_xy_type4(xVar,yVar,quasar_array,rlims)
        xs[i] = x[0]
        ys[i] = y[0]
    return xs,ys

#summary: In a type 1 plot, get sightline values for y and discrete sightline values for x
#
#inputs: xVar: xVar or xVar formula being displayed
#        yVar: ion or formula being displayed
#        quasarArray: array to use (from 'get_sightline_xy_vals')
#        rlims: what impact parameters to accept. (from 'get_sightline_xy_vals')
#        
#outputs: xs: raw x data for processing
#         ys: raw y data for processing  
def get_xy_type1(xVar,yVar,quasar_array,rlims):
    to_look_up = 'r' if xVar not in ["theta","phi"] else xVar
    distances = "kpc" if xVar == "r" or xVar == "rMpc" else "Rvir"
    gqary = []
    xs = np.empty(len(quasar_array),dtype = object)
    ys = np.empty(len(quasar_array),dtype = object)
    for i,array in enumerate(quasar_array):
        gq = multi_quasar_sphere.MultiQuasarSphere(array,distance=distances)
        if gq.number == 0:
            xs[i] = np.empty(0)
            ys[i] = np.empty(0)
            continue
        xs[i] = gq.get_sightline_values(to_look_up)
        if xVar == "rMpc":
            xs[i]/=1000
        ys[i] = get_yVar_from_str(gq,yVar)
        rs = gq.get_sightline_values('r')
        acceptedLines = np.logical_and(rlims[0]<=rs,rs<=rlims[1])
        xs[i] = xs[i][acceptedLines]
        ys[i] = ys[i][acceptedLines]
    return xs,ys
    
#summary: In a type 2 plot, get sightline values for y and metadata values for x
#
#inputs: xVar: xVar or xVar formula being displayed
#        yVar: ion or formula being displayed
#        quasarArray: array to use (from 'get_sightline_xy_vals')
#        rlims: what impact parameters to accept. (from 'get_sightline_xy_vals')
#        
#outputs: xs: raw x data for processing
#         ys: raw y data for processing    
def get_xy_type2(xVar,yVar,quasar_array,rlims):
    xs = np.empty(len(quasar_array),dtype = object)
    ys = np.empty(len(quasar_array),dtype = object)
    for i,ary in enumerate(quasar_array):
        xs[i] = np.empty(0)
        ys[i] = np.empty(0)
        for q in ary:
            x = q.__getattribute__(xVar)
            y = get_yVar_from_str(q,yVar)
            rs = q.get_sightline_values('r')
            acceptedLines = np.logical_and(rlims[0]<=rs/q.Rvir,\
                                           rs/q.Rvir<=rlims[1])
            y = y[acceptedLines]
            xs[i] = np.append(xs[i],np.tile(x,len(y)))
            ys[i] = np.append(ys[i],y)
    return xs,ys
    
#summary: In a type 3 plot, get values in the same way for both x and y. And they're just metadata table lookups
#
#inputs: xVar: xVar or xVar formula being displayed
#        yVar: ion or formula being displayed
#        quasarArray: array to use (from 'get_sightline_xy_vals')
#        
#outputs: xs: raw x data for processing
#         ys: raw y data for processing
def get_xy_type3(xVar,yVar,quasar_array,rlims):
    xs = np.empty(len(quasar_array),dtype = object)
    ys = np.empty(len(quasar_array),dtype = object)
    for i,ary in enumerate(quasar_array):
        xs[i] = np.empty(0)
        ys[i] = np.empty(0)
        for q in ary:
            x = q.__getattribute__(xVar)
            y = q.__getattribute__(yVar)
            xs[i] = np.append(xs[i],x)
            ys[i] = np.append(ys[i],y)
    return xs,ys
    
#summary: In a type 4 plot, get values in the same way for both x and y
#
#inputs: xVar: xVar or xVar formula being displayed
#        yVar: ion or formula being displayed
#        quasarArray: array to use (from 'get_sightline_xy_vals')
#        rlims: what impact parameters to accept. (from 'get_sightline_xy_vals')
#        
#outputs: xs: raw x data for processing
#         ys: raw y data for processing
def get_xy_type4(xVar, yVar, quasar_array, rlims):
    xs = np.empty(len(quasar_array),dtype = object)
    ys = np.empty(len(quasar_array),dtype = object)
    for i,ary in enumerate(quasar_array):
        xs[i] = np.empty(len(quasar_array[i]),dtype = object)
        ys[i] = np.empty(len(quasar_array[i]),dtype = object)
        for j,q in enumerate(ary):
            x = get_yVar_from_str(q,xVar)
            y = get_yVar_from_str(q,yVar)
            rs = q.get_sightline_values('r')
            acceptedLines = np.logical_and(rlims[0]<=rs/q.Rvir,\
                                           rs/q.Rvir<=rlims[1])
            x = x[acceptedLines]
            y = y[acceptedLines]
            xs[i][j] = x
            ys[i][j] = y
    return xs,ys


def mask_positive_vals(array_to_mask,masking_array,include_0):
    if include_0:
        m = np.logical_and(masking_array>=0,masking_array<np.inf)
    else:
        m = np.logical_and(masking_array>0,masking_array<np.inf)
    return array_to_mask[m]

#summary: Get raw data for processing from saved data. Basically this just 
#         calls one of four different functions depending on plot_type
#
#inputs: plot_type: 0,1,2,3,4 depending on kinds of variables (see 'decide_plot_type')
#        ion: ion or formula being displayed
#        xVar: xVar or xVar formula being displayed
#        lq: tuple of splitting, from 'sort_by'
#        quasarArray: array to use if not 'currentQuasarArray'
#        rlims: what impact parameters to accept. Default is .1 Rvir to infinity
#        
#outputs: xarys: raw x data for processing
#         yarys: raw y data for processing
def get_xy_vals(plot_type,xVar_packet,yVar_packet,quasar_array,mask_positive = True,rlims=None,**kwargs):
    xVar,yVar = xVar_packet[0],yVar_packet[0]
    if rlims is None:
        rlims = np.array([0.1,np.inf])
    elif rlims == "all":
        rlims = [0.0,np.inf]
    if plot_type==0:
        xarys,yarys = get_xy_type0(xVar,yVar,quasar_array,rlims)
    elif plot_type==1:
        xarys,yarys = get_xy_type1(xVar,yVar,quasar_array,rlims)
    elif plot_type==2:
        xarys,yarys = get_xy_type2(xVar,yVar,quasar_array,rlims)
    elif plot_type==3:
        xarys,yarys = get_xy_type3(xVar,yVar,quasar_array,rlims)
    elif plot_type==4:
        xarys,yarys = get_xy_type4(xVar,yVar,quasar_array,rlims)
    if isinstance(mask_positive,tuple):
        mask_positive,include_0 = mask_positive
    else:
        mask_positive,include_0 = mask_positive,False
    if mask_positive:
        if plot_type in [0,1,2,3]:
            for i in range(len(yarys)):
                xarys[i] = mask_positive_vals(xarys[i],yarys[i],include_0)
                yarys[i] = mask_positive_vals(yarys[i],yarys[i],include_0)
        elif plot_type in [4]:
            for i in range(len(yarys)):
                for j in range(len(yarys[i])):
                    xarys[i][j] = mask_positive_vals(xarys[i][j],yarys[i][j],include_0)
                    yarys[i][j] = mask_positive_vals(yarys[i][j],yarys[i][j],include_0)
                    yarys[i][j] = mask_positive_vals(yarys[i][j],xarys[i][j],include_0)
                    xarys[i][j] = mask_positive_vals(xarys[i][j],xarys[i][j],include_0)
    return xarys,yarys
    
