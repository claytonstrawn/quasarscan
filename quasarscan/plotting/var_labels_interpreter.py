import numpy as np
from quasarscan.utils import utils
from quasarscan.utils.variable_lists import stringcriteria,intensives,intensiveslabels,\
                                            sightline_xVars,param_xVars,\
                                            sightline_unit_labels,param_unit_labels,\
                                            all_known_variables


class UnknownPlotError(Exception):
    def __init__(self,message):
        self.message = message
        print(self.message)
#summary: given ion and xVar, figure out which kind of plot it is and 
#         check the other information given is consistent with that type
#
#inputs: ion: ion or formula being displayed. Can be a tuple (formula, name)
#        xVar: xVar or xVar formula being displayed. Can be a tuple (formula, name)
#        labels: list of labels for all lines to plot. if labels, return them as is
#        lq: this will include labels, if labels return as is
#        
#outputs: plot_type: 0,1,2,3,4 depending on kinds of variables
#             - Type 0 compares several different y variables for the same set of sightlines. Each curve is a different 
#               variable, and all variables are applied to the whole sample. Any xVar works.
#             - Type 1 compares a single y variable for the same set of sightlines. Each curve is a different 
#               subsample (i.e. a different galaxy snapshot/subset). Only sightline xVars accepted
#             - Type 2 compares a single y variable for the same set of sightlines. Each curve is a different 
#               subsample (i.e. a different galaxy subset). Only param xVars accepted (so probably only one datapoint per 
#               galaxy)
#             - Type 3 compares 2 params. Doesn't use our data really so not important. 
#               Just lets us fit in Nir's data catalogue into our plotting system
#             - Type 4 compares 2 ions or formulas against each other. Very useful for scatterplots.
#               colored by subset or snapshot.
def decide_plot_type(xVar,yVar,labels=None,quasar_array=None,lq=None):
    if lq:
        labels = lq[0]
    if isinstance(yVar,tuple):
        yVar = yVar[0]
    try:
        if isinstance(yVar,list):
            assert labels is None and quasar_array is None and xVar in sightline_xVars+param_xVars
            return 0
        elif isinstance(yVar,str) and xVar in sightline_xVars:
            return 1
        elif xVar in param_xVars and (not yVar in param_xVars):
            assert isinstance(yVar,str)
            return 2
        elif yVar in param_xVars and xVar in param_xVars:
            return 3
        elif xVar not in sightline_xVars + param_xVars:
            return 4
        else:
            assert False
    except AssertionError:
        raise UnknownPlotError('Plot type could not be detected with the variables chosen. x axis = "%s", y axis = "%s"'%(xVar,yVar))

def decide_quasar_array(qtype,default_quasar_array,quasar_array=None,lq=None,**kwargs):
    if quasar_array:
        quasar_array = quasar_array
    elif (quasar_array is None) and (lq is not None):
        index = {'sim':0,'obs':1,'empty':2}[qtype]
        quasar_array = lq[2][index]
    else:
        return np.array([default_quasar_array])
    final_quasar_array = np.empty(len(quasar_array),dtype=object)
    for i in range(len(quasar_array)):
        final_quasar_array[i] = []
        for q in quasar_array[i]:
            if q in default_quasar_array:
                final_quasar_array[i].append(q)
    final_quasar_array = np.array(final_quasar_array)
    return final_quasar_array





#summary: In a type 0 plot, or a plot with no labels specified, guess what the labels should be
#         also disaggregates ion formula and ion name if given
#
#inputs: plot_type: 0,1,2,3,4 depending on kinds of variables (see 'decide_plot_type')
#        ion: ion or formula being displayed. Can be a tuple (formula, name)
#        xVar: xVar or xVar formula being displayed. Can be a tuple (formula, name)
#        labels: list of labels for all lines to plot. if labels, return them as is
#        lq: this will include labels, if labels return as is
#        
#outputs: ion: formula to use to calculate y vals
#         ion_name: name to use for those y vals
#         xVar: formula to use to calculate x vals
#         xVar_name: name to use for those y vals
#         labels: list of list of labels for all lines to plot
def get_labels_from_vars(plot_type,xVar,yVar,labels,lq):
    if isinstance(yVar,tuple):
        yVar_name=yVar[1]
        yVar = yVar[0]
    else:
        yVar_name=yVar
    if isinstance(xVar,tuple):
        xVar_name = xVar[1]
        xVar = xVar[0]
    else:
        xVar_name = xVar
    if plot_type==0:
        yVars = yVar
        labels = []
        yVars_notuple = []
        for yVar in yVars:
            if isinstance(yVar,tuple):
                yVar_name = yVar[1]
                yVar = yVar[0]
            else:
                yVar_name = yVar
            labels.append(yVar_name)
            yVars_notuple.append(yVar)
            yVar_name = None
        yVar = yVars_notuple
        yVar_name = labels
    if labels is not None: 
        labels = labels
    elif labels is None and lq is not None:
        labels=lq[0]
    elif labels is None and lq is None:
        labels = ['all data']
    return xVar,xVar_name,yVar,yVar_name,labels

#summary: Detect whether x and y variables are naturally plotted as logs
#
#inputs: ion: ion or formula to use
#        xVar: xVar or formula to use
#        logx: whether to plot x logarithmically. If not 'guess', return whatever was given
#        logy: whether to plot y logarithmically. If not 'guess', return whatever was given
#        average: basically checking if you're plotting 'covering_fraction'
#        
#outputs: logx: whether to plot x logarithmically
#         logy: whether to plot y logarithmically
def should_take_logs_xy(xVar,yVar,logx,logy,average):
    probablylinear = ['redshift','rounded_redshift','a0','Rvir']
    if logx=='guess':
        if xVar in sightline_xVars:
            logx=False
        elif xVar in probablylinear:
            logx=False
        elif xVar in param_xVars:
            logx=True
        else:
            logx=True
    if logy=='guess':
        if yVar in sightline_xVars or average=='covering_fraction':
            logy=False
        elif yVar in param_xVars and xVar not in probablylinear:
            logy=True
        else:
            logy=True
    return logx,logy    


def what_to_filter_for(xVar,yVar):
    if not isinstance(yVar,list):
        yVars = [yVar]
    else:
        yVars = yVar
    filter_for = []
    for var in yVars+[xVar]:
        all_quantities = utils.split_by_ops(var)
        for s in all_quantities:
            if s == 'rdivR':
                continue
            elif s in all_known_variables:
                filter_for.append(s)
            elif len(s)>0:
                try:
                    eval(s)
                except:
                    filter_for.append(s)
    if xVar == 'rdivR':
        filter_for.append('Rvir')
    return filter_for

def configure_variables(xVar,yVar,average,logx='guess',logy='guess',labels=None,quasar_array=None,lq=None,**kwargs):
    plot_type = decide_plot_type(xVar,yVar,labels,quasar_array,lq)
    logx,logy = should_take_logs_xy(xVar,yVar,logx,logy,average)
    xVar,xVar_name,yVar,yVar_name,labels = get_labels_from_vars(plot_type,xVar,yVar,labels,lq)
    filter_for = what_to_filter_for(xVar,yVar)
    xVar_packet,yVar_packet = (xVar,xVar_name,logx),(yVar,yVar_name,logy)
    return plot_type,xVar_packet,yVar_packet,labels,filter_for

#summary: helper function for plot_err that gets the labels of things
#
#inputs: plot_type: 0,1,2,3,4 depending on kinds of variables (see 'decide_plot_type')
#        ion: ion or formula being displayed
#        ion_name: display name for ion. Especially useful for formulas
#        xVar: xVar or xVar formula being displayed
#        xVar_name: display name for xVar. Especially useful for formulas and for Type 4 plots
#        replacement_title: clear title and replace with this
#        extra_title: add this to the autogenerated title
#        average: what kind of average to use. See 'setPlots' for details. default is take mean and stderr
#        logx: whether to plot x logarithmically. Default is "guess", it'll ask 'should_take_logs_xy'
#        logy: whether to plot y logarithmically. Default is "guess", it'll ask 'should_take_logs_xy'
#        
#outputs:  xlabel: label for x axis
#          ylabel: label for y axis
#          title: label for top of figure
def get_labels_and_titles(plot_type,xVar_packet,yVar_packet,average,title=None,add_averaging_type=False,**kwargs):
    xVar,xVar_name,logx = xVar_packet
    yVar,yVar_name,logy = yVar_packet
    if isinstance(title,str):
        title = title
    else:
        if plot_type == 0:
            title = "CGM Sightline Data"
        elif plot_type == 1:
            title = yVar_name+" CGM Sightline Data"
        elif plot_type == 2:
            title = yVar_name+" CGM Parameter Dependence"
        elif plot_type == 3:
            title = "Galaxy 2 Parameter Dependence"
        elif plot_type == 4:
            title = "Sightline 2 Variable Correlation"
    if add_averaging_type:
        title+=" ("+str(average)+") "

    if plot_type == 0 and ":" in yVar[0] and yVar[0].split(":")[1] != "cdens":
        ylabel = "fraction of ion in state"
    elif plot_type == 0 and (":" not in yVar[0] or yVar.split(":")[1] == "cdens"):
        ylabel = "column density"
    elif yVar in intensives:
        ylabel = intensiveslabels[yVar]
    elif yVar in param_xVars:
        ylabel = param_unit_labels[yVar]
    elif average == 'covering_fraction':
        ylabel = yVar_name + ' covering fraction'
    else:
        ylabel = yVar_name
    if logy and average != 'covering_fraction':
        ylabel = 'log '+ylabel

    if xVar in intensives:
        xlabel = intensiveslabels[xVar]
    elif xVar in param_xVars:
        xlabel = param_unit_labels[xVar]
    elif xVar in sightline_xVars:
        xlabel = sightline_unit_labels[xVar]
    else:
        xlabel = xVar_name
    if logx:
        xlabel = 'log '+xlabel

    return xlabel,ylabel,title



