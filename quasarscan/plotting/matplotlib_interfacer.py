import matplotlib.pyplot as plt
matplotlib_default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
import numpy as np
from quasarscan.utils.utils import definecolorbar
matplotlib_default_symbols = np.array(['s', 'v', '^', '<', '>', '8', 'p', 'P', '*', 'h', 'H', '+', 'x', 'X', 'D', 'd', '|', '_', '1', '2', '3', '4'])
import logging
logging.getLogger().setLevel(logging.CRITICAL)

def transpose_err_array(errs):
    if errs is None:
        return None
    else:
        return errs.transpose()

#summary: helper function for plot_err that handles simulation data
#
#inputs: ax: axis (from 'plot_err')
#        plot_type: 0,1,2,3,4 depending on kinds of variables (see 'decide_plot_type')
#        xs: x values (calculated in 'process_datapoints')
#        ys: y values (calculated in 'process_datapoints')
#        xerrs: x errors (calculated in 'process_datapoints'). Should be None unless type 4 plot
#        yerrs: y values (calculated in 'process_datapoints')
#        xlabel: label to write along x axis (from 'get_title_and_axislabels')
#        ylabel: label to write along y axis (from 'get_title_and_axislabels')
#        title: label to write along top (from 'get_title_and_axislabels')
#        ylabel: label to write along y axis (from 'get_title_and_axislabels')
#        labels: list of labels for curves
#        average: what kind of average to use. See 'setPlots' for details. default is take mean and stderr
#        dots: If True, plot points instead of errorbars
#        grid: whether to plot grid underneath values
#        linestyle: how to connect points, default is no connection
#        linewidth: thickness of connecting line
#        fmt: style of points ('o','s','x','v', etc)
#        coloration: a list of colors. If none go through the default matplotlib colors
#        xlims: limits of x axis of plot. Default is matplotlib default
#        ylims: limits of y axis of plot. Default is matplotlib default
#        markersize: size of datapoints on plot
#        alpha: transparency of datapoints
#        elinewidth: thickness of errorbars. Default is linewidth
#        
#outputs:  future_colors: list of colors plotted to replicate in other plots if desired
def plot_sim_on_ax(plot_type,xs,ys,xerrs,yerrs,xlabel,ylabel,labels,title_final,ax=None,fig=None,
               average='default',dots=False,grid=True,linestyle='',ls='',linewidth = 1.5,
               fmt=None,coloration=None,xlims='default',ylims='default',markersize='default',
               alpha = 1.0,elinewidth=None,capsize=3,zorder = 100,legend = True,**kwargs):
    if ax is None:
        assert fig is None
        fig,ax = plt.subplots(1)
    lencolorlist = len(xs)//len(matplotlib_default_colors)+len(xs)%len(matplotlib_default_colors)
    coloration = coloration or np.tile(matplotlib_default_colors,lencolorlist)
    if xerrs is None:
        xerrs=[None]*len(xs)
    if yerrs is None:
        yerrs=[None]*len(xs)

    linestyle = linestyle if linestyle!='' else ls
    if markersize=='default':
        markersize=6
    if dots or plot_type in [3]:
        for i in range(len(xs)):
            fmt=fmt or 'o'
            ax.plot(xs[i],ys[i],marker = fmt,linestyle=linestyle,label=labels[i],zorder = zorder,
                                  color=coloration[i],markersize=markersize,alpha = alpha,linewidth=linewidth)
    else:
        fmtdict = {"mean":'.',"median_std":',',"covering_fraction":',',\
                   "stddev":'.',"median":".",'default':'.'}
        if not fmt:
            fmt = fmtdict[average]
        if fmt=='fill':
            for i in range(len(xs)):
                    ax.plot(xs[i],ys[i],label=labels[i],ls=linestyle,zorder = 100,markersize = markersize,\
                                 color = coloration[i],marker = fmtdict[average],alpha = alpha, linewidth=linewidth)
                    ax.fill_between(xs[i],ys[i]-yerrs[i][0,:],ys[i]+yerrs[i][1,:],color = future_colors[-1],\
                                    alpha = alpha/2,linewidth = linewidth/2,zorder = zorder/2)
        else:
            for i in range(len(xs)):
                    ax.errorbar(xs[i],ys[i],xerr=transpose_err_array(xerrs[i]),markersize = markersize,\
                                            yerr=transpose_err_array(yerrs[i]),label=labels[i],\
                                            ls=linestyle,color = coloration[i],fmt = fmt,capsize = capsize,\
                                            alpha = alpha, linewidth=linewidth, elinewidth=elinewidth,\
                                            zorder = zorder)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title_final)
    if xlims != 'default':
        ax.set_xlim(xlims)        
    if ylims != 'default':
        ax.set_ylim(ylims)
    if grid:
        if grid is True:
            grid = 'major'
        ax.grid(which = grid,visible = True)
    if legend:
        ax.legend().set_zorder(2)
    return fig,ax
#plot_type, xs, ys, xerrs, yerrs, xlabel, ylabel, labels, title, 

def plot_scatter_on_ax(plot_type,xs,ys,xlabel,ylabel,labels,title_final,ax=None,fig=None,\
                        grid=False,fmt=None,coloration=None,xlims='default',ylims='default',\
                        markersize='default',alpha = 1.0,legend = True,zorder = -100,stars = None,**kwargs):
    if ax is None:
        assert fig is None
        fig,ax = plt.subplots(1)
    lencolorlist = len(xs)//len(matplotlib_default_colors)+len(xs)%len(matplotlib_default_colors)
    coloration = coloration or np.tile(matplotlib_default_colors,lencolorlist)

    if markersize=='default':
        markersize=6
    if fmt is None:
        fmt = '.'
    if zorder == 'random':
        #this slows down plotting quite a bit -- instead of say six calls to ax.plot
        #which each plots around 200 points, it now does 1200 calls to ax.plot which 
        #each plot 1 point. I could not figure out any way to give points individual 
        #zorders without doing this however. This does successfully eliminate z-bias
        for i in range(len(xs)):
            zordering_bot = np.random.choice(np.arange(-100,0),len(xs[i]))
            for j in range(len(xs[i])):
                if j==0:
                    label=labels[i]
                else:
                    label = None
                ax.plot(xs[i][j],ys[i][j],linestyle='',label=label,color=coloration[i],
                        markersize=markersize,alpha = alpha,zorder=zordering_bot[j],marker=fmt)
    elif isinstance(zorder,str):
        print(f"I don't understand the command: zorder = {zorder}. Known commands: ['random'].")
    else:
        for i in range(len(xs)):
            ax.plot(xs[i],ys[i],linestyle='',label=labels[i],color=coloration[i],
                    markersize=markersize,alpha = alpha,zorder=zorder,marker=fmt)
            
    if stars is not None:
        used = False
        for i in range(len(xs)):
            label = labels[i]
            try:
                js_to_use = stars[label]
                used = True
            except KeyError:
                continue
            for j in js_to_use:
                ax.plot(xs[i][j],ys[i][j],linestyle='',label=None,color=coloration[i],
                        markersize=markersize*2,alpha = alpha,zorder=200,marker="*",markeredgecolor = 'k')
        if not used:
            print(f'stars was given as {stars} but no stars added.'+\
                  f'Did you have the correct keys? {labels}')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title_final)
    if xlims != 'default':
        ax.set_xlim(xlims)        
    if ylims != 'default':
        ax.set_ylim(ylims)
    if grid:
        if grid is True:
            grid = 'major'
        ax.grid(which = grid,visible = True)
    if legend:
        ax.legend().set_zorder(2)
    return fig,ax

def figure_out_limits(xerrs,i,j,default=.15,nan_err_behavior = 'zero'):
    if xerrs is None:
        return None,False,False
    xerr = xerrs[i][j]
    if np.isnan(xerr) or np.isinf(xerr):
        if nan_err_behavior == 'zero':
            xerr,xlolims,xuplims = 0,False,False
        elif nan_err_behavior == 'nan':
            xerr,xlolims,xuplims = np.nan,False,False
    elif xerr>=0:
        xerr,xlolims,xuplims = xerr,False,False
    elif xerr == -2:
        xerr,xlolims,xuplims = default,False,True
    elif xerr == -3:
        xerr,xlolims,xuplims = default,True,False
    return xerr,xlolims,xuplims

#summary: helper function for plot_err that handles observational data
#
#inputs: ion: the ion or formula to plot(can only use single ion b/c shows detailed ion info as color)
#        xVar: the xVar to plot. Can use any variable, will figure out what kind in 'decide_plot_type'
#        ax: axis (from 'plot_err')
#        fig: fig (from 'plot_err')
#        plot_type: 0,1,2,3,4 depending on kinds of variables (see 'decide_plot_type')
#        include_observations: 'only' or 'both' (if False, this func won't be called)
#        obs_data: package of data (from 'get_observational_data')
#        xlabel: label to write along x axis (from 'get_title_and_axislabels')
#        ylabel: label to write along y axis (from 'get_title_and_axislabels')
#        title: label to write along top (from 'get_title_and_axislabels')
#        labels: list of labels for curves
#        grid: whether to plot grid underneath values
#        linestyle: how to connect points, default is no connection (which is probably best in all cases for observations)
#        obs_coloration: if observations should be colored different from simulation, this is a list of colors
#        alpha: transparency of observations
#        fmt: style of points ('o','s','x','v', etc)
#        coloration: a list of colors, inherited from simulation colors. If obs_coloration, ignore this
#        xlims: limits of x axis of plot. Default is matplotlib default
#        ylims: limits of y axis of plot. efault is matplotlib default
#
#outputs:  None. This actually does the plot operation (ax.errorbar)


def plot_obs_on_ax(plot_type,xs,ys,xerrs,yerrs,xlabel,ylabel,labels,title,quasar_array,ax=None,fig=None,\
                        grid=False,symbols=None,coloration=None,xlims='default',ylims='default',\
                        markersize='default',alpha = 1.0,capsize=3,legend = True,\
                       **kwargs):
    if ax is None:
        assert fig is None
        fig,ax = plt.subplots(1)
    coloration = coloration or matplotlib_default_colors
    symbols = symbols or matplotlib_default_symbols

    author_symbol_dict = {}
    survey_color_dict = {}
    index_label_dict = {}
    for i in range(len(xs)):
        survey = labels[i]
        if survey not in survey_color_dict:
            survey_color_dict[survey] = coloration[len(survey_color_dict)]
        for j in range(len(xs[i])):
            author = quasar_array[i][j].author
            if author not in author_symbol_dict:
                author_symbol_dict[author] = symbols[len(author_symbol_dict)]
            index_label_dict[(i,j)] = (survey,author)


    used_labels = []
    for i in range(len(xs)):
        for j in range(len(xs[i])):
            (survey,author) = index_label_dict[(i,j)]
            color = survey_color_dict[survey]
            symbol = author_symbol_dict[author]
            label = survey+' - '+author
            if label not in used_labels:
                ax.errorbar([],[],color = color,fmt = symbol,capsize = capsize,\
                    alpha = alpha,label = label)
                used_labels.append(label)
            xerr,xlolims,xuplims = figure_out_limits(xerrs,i,j)
            yerr,ylolims,yuplims = figure_out_limits(yerrs,i,j)
            if not (xlolims or ylolims or xuplims or yuplims):
                ax.errorbar(xs[i][j],ys[i][j],xerr=xerr,yerr=yerr,color = color,\
                            fmt = symbol,capsize = 3,alpha = alpha)
            else:
                ax.errorbar(xs[i][j],ys[i][j],xerr=xerr,yerr=yerr,xlolims=xlolims,xuplims=xuplims,\
                            lolims=ylolims,uplims=yuplims,\
                             mec = color,ecolor = color,fmt = symbol,\
                             capsize = 3,mfc='w',alpha = alpha)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if xlims != 'default':
        ax.set_xlim(xlims)        
    if ylims != 'default':
        ax.set_ylim(ylims)
    if grid:
        if grid is True:
            grid = 'major'
        ax.grid(which = grid,visible = True)
    if legend:
        ax.legend().set_zorder(2)
    return fig,ax

def plot_hist_on_ax(plot_type,xs,ys,xlabel,ylabel,title,weight,cbarlabel,ax=None,fig=None,ns = (42,15),\
                    xlims='default',ylims='default',cbarlims = None,show_cbar = True,**kwargs):
    if ax is None:
        assert fig is None
        fig,ax = plt.subplots(1)
    cmap = definecolorbar(**kwargs)
    H, xedges, yedges = np.histogram2d(xs, ys, bins=ns,weights = weight)
    H = H.T
    X, Y = np.meshgrid(xedges, yedges)
    if cbarlims is None:
        cms=ax.pcolormesh(X,Y, H, cmap=cmap)
    else:
        cms=ax.pcolormesh(X,Y, H, cmap=cmap,vmin = cbarlims[0],vmax = cbarlims[1])
    if show_cbar:
        fig.colorbar(cms,label = cbarlabel,ax=ax)
    def get_new_axislim(current,divideBy = 40):
        x1 = current[0]
        x2 = current[1]
        return x1-(x2-x1)/divideBy,x2+(x2-x1)/divideBy
    if xlims=='default':
        ax.set_xlim(get_new_axislim(ax.get_xlim()))
    else:
        ax.set_xlim(xlims)
    if ylims=='default':
        ax.set_ylim(get_new_axislim(ax.get_ylim()))
    else:
        ax.set_ylim(ylims)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    return fig,ax

def setup_faberplot_subplots(lq2,fig,axes,figsize,sharex,sharey):
    ncols,nrows = len(lq2[0]),len(lq2[1])
    if fig:
        assert axes is not None
        assert axes.shape == (nrows,ncols)
    else:
        if figsize=='guess':
            figsize=(min(15,5*ncols),3.5*nrows)
        fig, axes = plt.subplots(nrows,ncols,figsize = figsize,sharex=sharex,sharey=sharey, squeeze=False)
    return fig,axes

def handle_faberplot_titles(i,j,axes,lq2):
    ax = axes[i,j]
    labels_x,labels_y = lq2[0],lq2[1]
    nrows,ncols = len(labels_x),len(labels_y)
    if i == 0:
        ax.set_xlabel('')
        ax.set_title(labels_x[j])
    if i == nrows-1:
        ax.set_title('')
    if i!=0 and i!= nrows-1:
        ax.set_xlabel('')
        ax.set_title('')
    if j == ncols-1:
        ax.set_ylabel('')
        right_ax = ax.twinx()
        right_ax.set_ylabel(labels_y[i])
        right_ax.yaxis.set_ticks_position('none') 
        right_ax.yaxis.set_ticklabels('') 
    if j!=0 and j!=ncols-1:
        ax.set_ylabel('')


