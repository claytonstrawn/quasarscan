import numpy as np
from quasarscan.utils.utils import split_by_ops,string_represents_ion
from quasarscan.utils.variable_lists import stringcriteria,intensives,intensiveslabels,\
                                            intensivespositions,sightline_xVars,param_xVars,\
                                            sightline_unit_labels,param_unit_labels,\
                                            all_known_variables
class BadVariableError(Exception):
    def __init__(self,message):
        self.message = message
        print(self.message)

class AveragesAndErrors(object):
    def __init__(self,plots,quartiles=None):
        if isinstance(plots,tuple):
            if len(plots) == 3:
                quartiles = (plots[1],plots[2])
            elif len(plots) == 2:
                if isinstance(plots[1],tuple):
                    quartiles = plots[1]
                elif plots[0] == 'median' or plots[0] == 'med':
                    quartiles = (50-plots[1],50+plots[1])
            plots = plots[0]
        if not quartiles:
            quartiles = (40,60)
        def getquartiles(data):
            return np.array([np.nanmedian(data) - np.nanpercentile(data,quartiles[0]),\
                             np.nanpercentile(data,quartiles[1])-np.nanmedian(data)])
        def getstderr(data):
            l = len(data[~np.isnan(data)])
            return np.array([np.nanstd(data)/np.sqrt(l),np.nanstd(data)/np.sqrt(l)])
        def getstddev(data):
            return np.array([np.nanstd(data),np.nanstd(data)])
        if plots in ["mean","stderr"]:
            self.plots = "mean"
            self.avgfn = np.nanmean
            self.errfn = getstderr
        if plots in ["stddev"]:
            self.plots = "stddev"
            self.avgfn = np.nanmean
            self.errfn = getstddev
        elif plots in ["median","med"]:
            self.plots = "median"
            self.avgfn = np.nanmedian
            self.errfn = getquartiles
        elif plots in ["med_noquartiles","median_std"]:
            self.plots = "median_std"
            self.avgfn = np.nanmedian
            self.errfn = getstderr
        elif plots in ["scatter"]:
            self.plots = "scatter"
            self.avgfn = None
            self.errfn = None
        elif plots in ['covering_fraction','cvf']:
            self.plots = "covering_fraction"
            self.avgfn = np.nanmean
            self.errfn = getstderr


#summary: helper for combine_xs. Checks if values are within tolerance
#
#inputs: x_in_list: value we're checking proximity to
#        x_comp: value we're checking for proximity of
#        tolerance: If x values are within tolerance, then count them as the same. 
#        symmetric: if true, check values above and below. If false, only check above
#                   (want to return from combine_xs the bottom end of the bins then)
#        
#outputs: True if close enough, False if not
def values_within_tolerance(x_in_list,x_comp,tolerance,symmetric = True):
    tocompare = x_in_list-x_comp
    if symmetric:
        tocompare = np.abs(tocompare)
    return np.logical_and(-tolerance<=tocompare,tocompare<=tolerance)

#summary: In a type 0, 1, or 2 plot, this is used to combine together a number of sightlines into 
#         single data points affiliated with certain xVar values. This function only creates the 
#         xVar values to use, they are applied in 'process_errbars_onlyvertical'
#
#inputs: x_variable: the set of all xVar numerical values
#        tolerance: If x values are within tolerance, then count them as the same. 
#        
#outputs: unique_xs (for use in 'process_errbars_onlyvertical')
def combine_xs(x_variable,tolerance):
    if tolerance == 0:
        return x_variable
    x_values_not_averaged = np.unique(x_variable)
    x_values = []
    i = 0
    while i < len(x_values_not_averaged)-1:
        currentList = [x_values_not_averaged[i]]
        numtocombine=1
        while i+numtocombine < len(x_values_not_averaged) and \
                values_within_tolerance(x_values_not_averaged[i+numtocombine],x_values_not_averaged[i],tolerance):
            currentList.append(x_values_not_averaged[i+numtocombine])
            numtocombine+=1
        i+=numtocombine
        x_values.append(np.min([currentList]))
    if i == len(x_values_not_averaged)-1:
        x_values.append(x_values_not_averaged[-1]) 
    return np.array(x_values)

#summary: if not plotting errorbars but instead plotting random individual sightlines 
#         data (average = 'scatter'), this will return the values needed
#
#inputs: xarys: raw x data (from 'get_sightline_xy_vals')
#        yarys: raw y data (from 'get_sightline_xy_vals')
#        logx: whether to plot x logarithmically. Default is "guess", it'll ask 'should_take_logs_xy'
#        logy: whether to plot y logarithmically. Default is "guess", it'll ask 'should_take_logs_xy'
#        subsample: only plot this fraction of points. (half of points, 1/3 of points etc)
#        offsetx: if True, calculate a small offset for x values and add it. (so errorbars don't overlap)
#        tolerance: If x values are within tolerance, then count them as the same. The only effect on this is 
#                   if offsetx is true
#        
#outputs: xs: processed x data to plot
#         ys: processed y data to plot
#         empty: if True, no data points to plot
def process_scatter_points(xVar_packet,yVar_packet,xarys,yarys,subsample=1.0,offsetx=False,tolerance=1e-5,**kwargs):
    logx,logy = xVar_packet[2],yVar_packet[2]
    def flatten_if_needed(ary):
        try:
            for item in ary:
                len(item)
            return np.concatenate(ary)
        except:
            return ary
    def getsubsample(data,frac):
        l = int(len(data)*frac)
        mask = np.zeros(len(data))
        mask[:l] = 1
        np.random.shuffle(mask)
        return mask.astype(bool)
    xs = np.zeros(len(xarys),dtype=object)
    ys = np.zeros(len(xarys),dtype=object)
    empty = True
    for i in range(len(xarys)):
        xs[i] = flatten_if_needed(xarys[i])
        ys[i] = flatten_if_needed(yarys[i])
        mask = getsubsample(xs[i],subsample)
        xs[i]=xs[i][mask]
        ys[i]=ys[i][mask]
        if np.any(ys[i]>0) and np.any(xs[i]>0):
            empty=False
        if logx:
            xs[i] = np.log10(xs[i])
        if logy:
            ys[i] = np.log10(ys[i])
        if offsetx and np.any(ys[i]>0) and np.any(xs[i]>0):
            decimals = int(-np.log10(tolerance))
            try:
                randomscatterwidth = np.min(np.diff(np.unique(np.round(xs[i],decimals = decimals))))
            except ValueError:
                randomscatterwidth = 0.
            add = (np.random.random(len(xs[i]))-.5)*randomscatterwidth*.5
            xs[i] += add
    return xs,ys,empty
    
#summary: calculate the averages and errorbars in any plot except type 4
#
#inputs: xarys: raw x data (from 'get_sightline_xy_vals')
#        yarys: raw y data (from 'get_sightline_xy_vals')
#        logx: whether to plot x logarithmically. Default is "guess", it'll ask 'should_take_logs_xy'
#        logy: whether to plot y logarithmically. Default is "guess", it'll ask 'should_take_logs_xy'
#        offsetx: if True, calculate a small offset for x values and add it. (so errorbars don't overlap)
#        tolerance: If x values are within tolerance, then count them as the same. Might need to change for ssfr
#        
#outputs: xs: processed x data to plot
#         ys: processed y data to plot
#         yerrs: processed y error data to plot
#         empty: if True, no data points to plot
def process_errbars_onlyvertical(xVar_packet,yVar_packet,xarys,yarys,averager,offsetx=False,tolerance=1e-5,**kwargs):
    logx,logy = xVar_packet[2],yVar_packet[2]
    xs = np.empty(len(yarys),dtype = object)
    ys = np.empty(len(yarys),dtype = object)
    yerrs = np.empty(len(yarys),dtype = object)
    empty = True
    for i in range(len(yarys)):
        unique_xs = combine_xs(xarys[i],tolerance)
        xs[i] = unique_xs
        ys[i] = np.zeros(len(unique_xs))
        yerrs[i] = np.zeros((len(unique_xs),2))
        mask = np.ones(len(unique_xs),dtype = bool)
        for j in range(len(unique_xs)):
            x_value = unique_xs[j]
            yvals = yarys[i][values_within_tolerance(xarys[i],x_value,tolerance,symmetric=False)]
            if len(yvals) == 0:
                mask[j] = False
                continue
            ys[i][j] = averager.avgfn(yvals)
            yerrs[i][j] = averager.errfn(yvals)
            if ys[i][j] > 0:
                empty = False
            if logy:
                yerrs[i][j][0] = np.log10(ys[i][j])-np.log10(ys[i][j]-yerrs[i][j][0])
                if yerrs[i][j][0] == np.inf or np.isnan(yerrs[i][j][0]):
                    yerrs[i][j][0] = max(np.log10(ys[i][j]),0)
                yerrs[i][j][1] = np.log10(ys[i][j]+yerrs[i][j][1])-np.log10(ys[i][j])
                ys[i][j] = np.log10(ys[i][j])
        xs[i] = xs[i][mask]
        ys[i] = ys[i][mask]
        yerrs[i] = yerrs[i][mask]
        if isinstance(offsetx,tuple):
            xs[i]+=offsetx[1]
            offsetx_bool=offsetx[0]
        elif isinstance(offsetx,bool):
            offsetx_bool = offsetx
        elif isinstance(offsetx,float) or isinstance(offsetx,int):
            xs[i]+=offsetx
            offsetx_bool = False
        if offsetx_bool and len(xs[i])>0:
            scatterwidth = np.min(np.diff(np.unique(xs[i])))
            if logx:
                xs[i][xs[i]>0]*=10**((float(i)/len(xs)-.5)*scatterwidth*.2)
            else:
                xs[i]+=(float(i)/len(xs)-.5)*scatterwidth*.5

        if logx:
            xs[i] = np.log10(xs[i])
    return xs,ys,yerrs,empty
    
#summary: calculate the averages and errorbars in type 4 plot
#
#inputs: xarys: raw x data (from 'get_sightline_xy_vals')
#        yarys: raw y data (from 'get_sightline_xy_vals')
#        logx: whether to plot x logarithmically. Default is "guess", it'll ask 'should_take_logs_xy'
#        logy: whether to plot y logarithmically. Default is "guess", it'll ask 'should_take_logs_xy'
#        
#outputs: xs: processed x data to plot
#         ys: processed y data to plot
#         xerrs: processed x error data to plot (if type 4, otherwise None)
#         yerrs: processed y error data to plot
#         empty: if True, no data points to plot
def process_errbars_vertandhoriz(xVar_packet,yVar_packet,xarys,yarys,averager,**kwargs):
    logx,logy = xVar_packet[2],yVar_packet[2]
    xs = np.empty(len(yarys),dtype = object)
    ys = np.empty(len(yarys),dtype = object)
    xerrs = np.empty(len(yarys),dtype = object)
    yerrs = np.empty(len(yarys),dtype = object)
    empty = True
    for i in range(len(xs)):
        xs[i] = np.zeros(len(xarys[i]))
        ys[i] = np.zeros(len(xarys[i]))
        xerrs[i] = np.zeros((len(xarys[i]),2))
        yerrs[i] = np.zeros((len(xarys[i]),2))
        for j,xlist in enumerate(xarys[i]):
            xs[i][j] = averager.avgfn(xlist)
            xerrs[i][j] = averager.errfn(xlist)
            if logx:
                xerrs[i][j][0] = np.log10(xs[i][j])-np.log10(xs[i][j]-xerrs[i][j][0])
                xerrs[i][j][1] = np.log10(xs[i][j]+xerrs[i][j][1])-np.log10(xs[i][j])
                xs[i][j] = np.log10(xs[i][j])
        for j,ylist in enumerate(yarys[i]):
            ys[i][j] = averager.avgfn(ylist)
            yerrs[i][j] = averager.errfn(ylist)
            if (xs[i][j]>0 and ys[i][j]>0) or (logx and ys[i][j]>0):
                empty = False
            if logy:
                yerrs[i][j][0] = np.log10(ys[i][j])-np.log10(ys[i][j]-yerrs[i][j][0])
                yerrs[i][j][1] = np.log10(ys[i][j]+yerrs[i][j][1])-np.log10(ys[i][j])
                ys[i][j] = np.log10(ys[i][j])
    return xs,ys,xerrs,yerrs,empty

def get_sim_errs(plot_type,xVar_packet,yVar_packet,xarys,yarys,average = 'median',quartiles = None,**kwargs):
    averager = AveragesAndErrors(average,quartiles)
    if plot_type in [0,1,2,3]:
        xs,ys,yerrs,empty = process_errbars_onlyvertical(xVar_packet,yVar_packet,xarys,yarys,averager,**kwargs)
        xerrs = None
    elif plot_type in [4]:
        xs,ys,xerrs,yerrs,empty = process_errbars_vertandhoriz(xVar_packet,yVar_packet,xarys,yarys,averager,**kwargs)
    return xs,ys,xerrs,yerrs,empty

def get_xVar_err(xVar):
    if len(split_by_ops(xVar))>1:
        return None
    if ':' in xVar:
        if len(xVar.split(':')) == 2 and xVar.split(':')[1] == 'cdens' and \
                        string_represents_ion(xVar.split(':')[0]):
            return xVar.split(':')[0]+':eb'
        else:
            return None
    if string_represents_ion(xVar):
        return xVar+':eb'
    elif xVar in param_xVars:
        return xVar+'_eb'
    else:
        return None


def handle_scatter_errs(xVar_packet,yVar_packet,quasar_array):
    xVar,yVar = xVar_packet[0],yVar_packet[0]
    logx,logy = xVar_packet[2],yVar_packet[2]
    xVar_err = get_xVar_err(xVar)
    if xVar_err is None:
        xerrs = None
    yVar_err = get_xVar_err(yVar)
    if yVar_err is None:
        yerrs = None

    if xVar_err is not None:
        xerrs = np.empty(len(quasar_array),dtype=object)
        for i,l in enumerate(quasar_array):
            xerrs[i] = np.zeros(len(l))
            for j,q in enumerate(l):
                if xVar in param_xVars:
                    xerrs[i][j] = eval('q.%s'%xVar_err)
                elif string_represents_ion(xVar):
                    xerrs[i][j] = q.get_ion_values(xVar_err)[0]
                xerrs[i][j] = np.log10(xerrs[i][j]) if logx and xerrs[i][j]>0 else xerrs[i][j]
    if yVar_err is not None:
        yerrs = np.empty(len(quasar_array),dtype=object)
        for i,l in enumerate(quasar_array):
            yerrs[i] = np.zeros(len(l))
            for j,q in enumerate(l):
                if yVar in param_xVars:
                    yerrs[i][j] = eval('q.%s'%yVar_err)
                elif string_represents_ion(yVar):
                    yerrs[i][j] = q.get_ion_values(yVar_err)[0]
                yerrs[i][j] = np.log10(yerrs[i][j]) if logy and yerrs[i][j]>0 else yerrs[i][j]
    return xerrs,yerrs

def process_xy_vals_hist(xVar_packet,yVar_packet,xarys,yarys,tolerance=1e-5,weights=True,**kwargs):
    assert len(xarys) == 1 and len(yarys) == 1
    xs,ys = xarys[0],yarys[0]
    xVar,_,logx = xVar_packet
    yVar,_,logy = yVar_packet
    unique_xs = combine_xs(xs,tolerance)
    retxs = xs*0.0
    for x_value in unique_xs:
        mask = values_within_tolerance(xs,x_value,tolerance,symmetric=False).astype(float)
        retxs+= (mask*x_value)
    if weights:
        weight = xs*0.0
        for i,x in enumerate(xs):
            weight[i] = 1/np.sum(values_within_tolerance(xs,x,tolerance,symmetric=False).astype(float))
        cbarlabel = "Fraction of lines for fixed %s"%(xVar)
    else:
        weight = xs*0.0+1.0
        cbarlabel = "Total number of lines"
    xs=retxs
    if logx:
        xs = np.log10(xs)
    if logy:
        ys = np.log10(ys)
    return xs,ys,weight,cbarlabel,len(xs)<=0