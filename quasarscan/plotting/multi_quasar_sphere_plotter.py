import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import datetime
from functools import reduce

from quasarscan.preprocessing import parse_metadata
from quasarscan.utils import roman,ion_lists
from quasarscan.utils.utils import sort_ions,reversearray
from quasarscan.utils.variable_lists import stringcriteria,intensives,intensiveslabels,\
                                            sightline_xVars,param_xVars,\
                                            sightline_unit_labels,param_unit_labels,\
                                            all_known_variables
from quasarscan.plotting.quasar_array_handler import QuasarArrayHandler
from quasarscan.plotting import var_labels_interpreter,\
                                plot_data_processor,\
                                errorbar_processor,\
                                matplotlib_interfacer

#these are a number of global lists and dictionaries which are checkedvariable_lists against in various places

class BadPlotError(Exception):
    def __init__(self,message):
        self.message = message
        print(self.message)

class MultiQuasarSpherePlotter():
    #summary: initialize mq and load all data
    #
    #inputs: loadonly: if not 'all' only load certain simulations (e.g. 'VELA')
    #        loadobs: if not 'all' only load certain simulations (e.g. 'COS-Halos') [I'm not sure this works]
    #        textfiles: usually None. Could give a specific list of textfiles to use if you don't want to search here
    #        cleanup: delete textfiles that fail safety check (will ask user permission first)
    #        plots: default plots value to use. Will default to mean if not given
    #        throwErrors: if true, throw errors when reading textfiles if broken. If false, skip ones that create errors
    #        safetycheck: if False, skip safetychecks and just use whatever you load
    #
    #outputs: MultiQuasarSpherePlotter object, usually called 'mq'
    def __init__(self, loadsim = "all",loadobs = 'all',loadempty = 'none',average = 'median'):
        self.defaultaverage = average
        self.quasar_array_handler = QuasarArrayHandler(loadsim,loadobs,loadempty)

    #summary: print length of list of QuasarSpheres
    #
    #inputs: None
    #
    #outputs: length of currentQuasarArray
    def length(self,include_nonsims=False):
        return self.quasar_array_handler.length(include_nonsims)

    def list_all_quasar_spheres(self,*criteria,qtype='sim',log=False):
        self.quasar_array_handler.list_all_quasar_spheres(*criteria,qtype=qtype,log=log)

    #summary: cancel all constraints
    #
    #inputs: None
    #        
    #outputs: None, changes state of mq.currentQuasarArray and mq.currentQuasarArrayName
    def reset_current_quasar_array(self):
        self.quasar_array_handler.reset_current_quasar_array()
    
    #summary: splits currentQuasarArray into particular bins, either calculated on the fly or given
    #
    #inputs: criteria: what criteria to constrain by (simname, simnum, Mvir, Rvir, SFR, 
    #                           etc. See 'quasar_sphere.py' for full list)
    #        bins: a list of n numbers, which determine n-1 bins in between them. If string param, n strings which each
    #              constitute a bin
    #        atEnd: if True, compare values by their final (z=1 or z=minimum among remaining values) value, 
    #               not the current one
    #        splitEven: a number of bins to split into. The bins will be chosen so each has the same number of members.
    #        **kwargs: onlyNonempty,reverse ['postprocess_sorted']
    #        
    #outputs: labels: list of strings for labelling points in legend
    #         bins: the bins to compare to 
    #         quasarBins: list of lists of quasarSphere objects that fit in the bins
    #         obsBins: list of lists of observationalQuasarSphere objects that fit in the bins
    #          
    #         NOTE: These are usually combined together and considered an 'lq' object, passed directly 
    #            into most plots (any except type 0) general use case is e.g.
    #            >>>lq = mq.sort_by('Mvir',[0,10**11,np.inf])
    #            >>>mq.plot_err('O VI',lq=lq)
    def sort_by(self, criteria, bins = [0,np.inf],at_end = False,split_even = False,reverse=False,**kwargs):
        return self.quasar_array_handler.sort_by(criteria,bins,at_end,split_even,reverse,**kwargs)

    def sort_by_2D(self, criteria_x,criteria_y, bins_x = [0,np.inf],bins_y = [0,np.inf],\
                    at_end_x = False,at_end_y = False,split_even_x = False,split_even_y = False,\
                    reverse_x = False,reverse_y = False,**kwargs):
        return self.quasar_array_handler.sort_by_2D(criteria_x,criteria_y, bins_x,bins_y,at_end_x,at_end_y,\
                                                split_even_x,split_even_y,reverse_x,reverse_y,**kwargs)

    #summary: restricts to only quasarspheres with galaxy parameters within certain limits
    #
    #inputs: constrainCriteria: what criteria to constrain by (simname, simnum, Mvir, Rvir, SFR, 
    #                           etc. See 'quasar_sphere.py' for full list)
    #        bins: either a list of two numbers, if a numerical criteria, or several strings if string criteria
    #              can leave as None if splitEven is used
    #        **kwargs: changeArrayName, exclude ['change_array_name']
    #                  splitEven,atEnd,set_main_array,sortobs ['constrain_array_helper']
    #        
    #outputs: return the bins used (in the case of 'low' or 'high' for example it'll tell you the cutoff)
    def constrain_current_quasar_array(self, constrain_criteria,bins=None,qtype='all',**kwargs):
        return self.quasar_array_handler.constrain_current_quasar_array(constrain_criteria,bins,qtype,**kwargs)
    
    def plot_err(self,yVar,xVar='rdivR',qtype = 'sim',average = 'default',force_averaging = False,**kwargs):
        average = self.defaultaverage if average == 'default' else average
        average = 'scatter' if qtype == 'obs' and not force_averaging else average
        plot_type,xVar_packet,yVar_packet,labels,filter_for = var_labels_interpreter.configure_variables(xVar,yVar,average,**kwargs)
        unfiltered_qlist = self.quasar_array_handler.impose_requirements(filter_for,qtype)
        xlabel,ylabel,title_final = var_labels_interpreter.get_labels_and_titles(plot_type,xVar_packet,yVar_packet,average,**kwargs)
        quasar_array = var_labels_interpreter.decide_quasar_array(qtype,self.quasar_array_handler.get_qlist(qtype),**kwargs)
        xarys,yarys = plot_data_processor.get_xy_vals(plot_type,xVar_packet,yVar_packet,quasar_array,**kwargs)
        if qtype == 'sim' or force_averaging:
            xs,ys,xerrs,yerrs,empty = errorbar_processor.get_sim_errs(plot_type,xVar_packet,yVar_packet,xarys,yarys,average = average,**kwargs)
            if not empty:
                to_return = matplotlib_interfacer.plot_sim_on_ax(plot_type, xs, ys, xerrs, yerrs, xlabel, ylabel, labels, title_final, **kwargs)        
        elif qtype == 'obs':
            xs,ys,empty = errorbar_processor.process_scatter_points(xVar_packet,yVar_packet,xarys,yarys,**kwargs)
            xerrs,yerrs = errorbar_processor.handle_scatter_errs(xVar_packet,yVar_packet,quasar_array)
            if not empty:
                to_return = matplotlib_interfacer.plot_obs_on_ax(plot_type, xs, ys, xerrs, yerrs, xlabel, ylabel, labels, title_final, quasar_array, **kwargs)
        elif qtype == 'empty':
            assert plot_type == 3
            xs,ys,empty = errorbar_processor.process_scatter_points(xVar_packet,yVar_packet,xarys,yarys,**kwargs)
            xerrs,yerrs = None,None
            if not empty:
                to_return = matplotlib_interfacer.plot_sim_on_ax(plot_type, xs, ys, xerrs, yerrs, xlabel, ylabel, labels, title_final, **kwargs)
        if empty:
            to_return = None,None
        self.quasar_array_handler.update_qlist(qtype,unfiltered_qlist)
        return to_return

    def plot_scatter(self,yVar,xVar='rdivR',qtype = 'sim',**kwargs):
        plot_type,xVar_packet,yVar_packet,labels,filter_for = var_labels_interpreter.configure_variables(xVar,yVar,'scatter',**kwargs)
        unfiltered_qlist = self.quasar_array_handler.impose_requirements(filter_for,qtype)
        xlabel,ylabel,title = var_labels_interpreter.get_labels_and_titles(plot_type,xVar_packet,yVar_packet,'scatter',**kwargs)
        quasar_array = var_labels_interpreter.decide_quasar_array(qtype,self.quasar_array_handler.get_qlist(qtype),**kwargs)
        xarys,yarys = plot_data_processor.get_xy_vals(plot_type,xVar_packet,yVar_packet,quasar_array,**kwargs)
        xs,ys,empty = errorbar_processor.process_scatter_points(xVar_packet,yVar_packet,xarys,yarys,**kwargs)
        if qtype == 'obs':
            xerrs,yerrs = errorbar_processor.handle_scatter_errs(xVar_packet,yVar_packet,quasar_array)
            if not empty:
                to_return = matplotlib_interfacer.plot_obs_on_ax(plot_type, xs, ys, xerrs, yerrs, xlabel, ylabel, labels, title, quasar_array, **kwargs)
        elif qtype in ['sim','empty']:
            xerrs,yerrs = None,None
            if not empty:
                to_return = matplotlib_interfacer.plot_scatter_on_ax(plot_type, xs, ys, xlabel, ylabel, labels, title, **kwargs)
        if empty:
            to_return = None,None
        self.quasar_array_handler.update_qlist(qtype,unfiltered_qlist)
        return to_return

    def plot_hist(self,yVar,xVar='rdivR',qtype = 'sim',**kwargs):
        if qtype != 'sim':
            raise BadPlotError('can only plot simulation sightlines for plot_hist')
        plot_type,xVar_packet,yVar_packet,labels,filter_for = var_labels_interpreter.configure_variables(xVar,yVar,'scatter',**kwargs)
        assert plot_type in [1,2], "'plot_hist' can only plot one continuous variable against one discrete variable"
        unfiltered_qlist = self.quasar_array_handler.impose_requirements(filter_for,qtype)
        xlabel,ylabel,title = var_labels_interpreter.get_labels_and_titles(plot_type,xVar_packet,yVar_packet,'hist',**kwargs)
        quasar_array = var_labels_interpreter.decide_quasar_array('sim',self.quasar_array_handler.get_qlist('sim'),**kwargs)
        xarys,yarys = plot_data_processor.get_xy_vals(plot_type,xVar_packet,yVar_packet,quasar_array,**kwargs)
        xs,ys,weight,cbarlabel,empty = errorbar_processor.process_xy_vals_hist(xVar_packet,yVar_packet,xarys,yarys,**kwargs)
        if not empty:
            to_return = matplotlib_interfacer.plot_hist_on_ax(plot_type, xs, ys, xlabel, ylabel, title, weight, cbarlabel, **kwargs)
        else:
            to_return = None,None
        self.quasar_array_handler.update_qlist(qtype,unfiltered_qlist)
        return to_return

    def faberplot(self,yVar,xVar='rdivR',plot_kind='err',lq2=None,qtype='sim',lq=None,fig = None, axes = None,figsize='guess',sharex=True,sharey=True,\
                  **kwargs):
        #after using sort_by_2d to get a set of labels and a 2d array of quasarspheres,
        #ask plot_err or plot_hist for completed plots of type given, for 
        #quasars in that cell of quasarArray, put them in subplots of a n by m subplots object
        #and show that plot
        if lq2 is None:
            raise BadPlotError('required to use sort_by_2D before plotting a faberplot')
        fig,axes = matplotlib_interfacer.setup_faberplot_subplots(lq2,fig,axes,figsize,sharex,sharey)
        old_quasar_array = self.quasar_array_handler.get_qlist(qtype)
        quasar_array = lq2[4][self.quasar_array_handler.get_qtype_index(qtype)]
        for i,axlist in enumerate(axes):
            for j,ax in enumerate(axlist):
                self.quasar_array_handler.update_qlist(qtype,quasar_array[i][j])
                if lq is not None:
                    lq = self.sort_by(lq[3],lq[1],**kwargs)
                    lq = [None]*len(lq[0]) if i>0 or j>0 else lq[0],lq[1],lq[2],lq[3]
                if plot_kind=='err':
                    self.plot_err(yVar, xVar=xVar, qtype=qtype, fig = fig,ax = ax, lq=lq,  **kwargs)
                elif plot_kind=='scatter':
                    self.plot_scatter(yVar, xVar=xVar, qtype=qtype, fig = fig,ax = ax, lq=lq, **kwargs)
                elif plot_kind=='hist':
                    self.plot_hist(yVar, xVar=xVar, qtype=qtype, fig = fig,ax = ax, **kwargs)
                matplotlib_interfacer.handle_faberplot_titles(i,j,axes,lq2)
        self.quasar_array_handler.update_qlist(qtype,old_quasar_array)













