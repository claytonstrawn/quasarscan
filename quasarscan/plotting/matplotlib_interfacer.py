import matplotlib.pyplot as plt
matplotlib_default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
import numpy as np
from quasarscan.utils.utils import definecolorbar
matplotlib_default_symbols = np.array(['s', 'v', '^', '<', '>', '1', '2', '3', '4', '8', 'p', 'P', '*', 'h', 'H', '+', 'x', 'X', 'D', 'd', '|', '_'])
import logging
logging.getLogger().setLevel(logging.CRITICAL)

def transpose_err_array(errs):
    if errs is None:
        return None
    else:
        return errs.transpose()
    
    
class ZBiasFreePlotter(object):
    def __init__(self):
        self.plot_calls = []

    def add_plot(self, f, xs, ys, *args, **kwargs):
        self.plot_calls.append((f, xs, ys, args, kwargs))

    def draw_plots(self, chunk_size=512):
        scheduled_calls = []
        for f, xs, ys, args, kwargs in self.plot_calls:
            assert(len(xs) == len(ys))
            index = np.arange(len(xs))
            np.random.shuffle(index)
            index_blocks = [index[i:i+chunk_size] for i in np.arange(len(index))[::chunk_size]]
            for i, index_block in enumerate(index_blocks):
                # Only attach a label for one of the chunks
                if i != 0 and kwargs.get("label") is not None:
                    kwargs = kwargs.copy()
                    kwargs["label"] = None
                scheduled_calls.append((f, xs[index_block], ys[index_block], args, kwargs))

        np.random.shuffle(scheduled_calls)

        for f, xs, ys, args, kwargs in scheduled_calls:
            f(xs, ys, *args, **kwargs)

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
               alpha = 1.0,elinewidth=None,capsize=3,zorder = 100,**kwargs):
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

    if dots or plot_type in [3]:
        for i in range(len(xs)):
            fmt=fmt or 'o'
            if markersize=='default':
                markersize=6
            ax.plot(xs[i],ys[i],marker = fmt,linestyle=linestyle,label=labels[i],zorder = zorder,
                                  color=coloration[i],markersize=markersize,alpha = alpha,linewidth=linewidth)
    else:
        fmtdict = {"mean":'.',"median_std":',',"covering_fraction":',',"stddev":'.',"median":".",'default':'.'}
        if not fmt:
            fmt = fmtdict[average]
        if fmt=='fill':
            for i in range(len(xs)):
                    ax.plot(xs[i],ys[i],label=labels[i],ls=linestyle,zorder = 100,\
                                 color = coloration[i],marker = fmtdict[average],alpha = alpha, linewidth=linewidth)
                    ax.fill_between(xs[i],ys[i]-yerrs[i][0,:],ys[i]+yerrs[i][1,:],color = future_colors[-1],\
                                    alpha = alpha/2,linewidth = linewidth/2,zorder = zorder/2)
        else:
            for i in range(len(xs)):
                    ax.errorbar(xs[i],ys[i],xerr=transpose_err_array(xerrs[i]),\
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
        ax.grid()
    ax.legend().set_zorder(2)
    return fig,ax
#plot_type, xs, ys, xerrs, yerrs, xlabel, ylabel, labels, title, 

def plot_scatter_on_ax(plot_type,xs,ys,xlabel,ylabel,labels,title_final,ax=None,fig=None,\
                        grid=False,fmt=None,coloration=None,xlims='default',ylims='default',\
                        markersize='default',alpha = 1.0,zorder = -100,**kwargs):
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
        bias_free_plotter = ZBiasFreePlotter()
        for i in range(len(xs)):
            bias_free_plotter.add_plot(ax.plot, xs[i],ys[i],'o',label=labels[i],
                                       color=coloration[i],markersize=markersize,
                                       alpha = alpha,marker=fmt)
        bias_free_plotter.draw_plots()
    else:
        for i in range(len(xs)):
            ax.plot(xs[i],ys[i],'o',label=labels[i],color=coloration[i],markersize=markersize,
                                    alpha = alpha,zorder=zorder,marker=fmt)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title_final)
    if xlims != 'default':
        ax.set_xlim(xlims)        
    if ylims != 'default':
        ax.set_ylim(ylims)
    if grid:
        ax.grid()
    ax.legend().set_zorder(2)
    return fig,ax

def figure_out_limits(xerrs,i,j,default=.15):
    if xerrs is None:
        return None,False,False
    xerr = xerrs[i][j]
    if np.isnan(xerr) or np.isinf(xerr):
        xerr,xlolims,xuplims = 0,False,False
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
                        markersize='default',alpha = 1.0,capsize=3,**kwargs):
    if ax is None:
        assert fig is None
        fig,ax = plt.subplots(1)
    lencolorlist = len(xs)//len(matplotlib_default_colors)+len(xs)%len(matplotlib_default_colors)
    coloration = coloration or np.tile(matplotlib_default_colors,lencolorlist)
    symbols = symbols or np.tile(matplotlib_default_symbols,lencolorlist)
    citations = []
    used_labels = []
    for i in range(len(xs)):
        color = coloration[i]
        for j in range(len(xs[i])):
            current_label = quasar_array[i][j].author+'-'+labels[i]
            if current_label[0] not in citations:
                citations.append(current_label[0])
            if current_label not in used_labels:
                used_labels.append(current_label)
                fmt = symbols[len(citations)-1]
                ax.errorbar([],[],color = color,fmt = fmt,capsize = capsize,\
                    alpha = alpha,label = current_label)
            xerr,xlolims,xuplims = figure_out_limits(xerrs,i,j)
            yerr,ylolims,yuplims = figure_out_limits(yerrs,i,j)
            if not (xlolims or ylolims or xuplims or yuplims):
                ax.errorbar(xs[i][j],ys[i][j],xerr=xerr,yerr=yerr,color = color,\
                            fmt = fmt,capsize = 3,alpha = alpha)
            else:
                ax.errorbar(xs[i][j],ys[i][j],xerr=xerr,yerr=yerr,xlolims=xlolims,xuplims=xuplims,\
                            lolims=ylolims,uplims=yuplims,\
                             mec = color,ecolor = color,fmt = fmt,\
                             capsize = 3,mfc='w',alpha = alpha)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if xlims != 'default':
        ax.set_xlim(xlims)        
    if ylims != 'default':
        ax.set_ylim(ylims)
    if grid:
        ax.grid()
    ax.legend().set_zorder(2)
    return fig,ax

def plot_hist_on_ax(plot_type,xs,ys,xlabel,ylabel,title,weight,cbarlabel,ax=None,fig=None,ns = (42,15),\
                    xlims='default',ylims='default',**kwargs):
    if ax is None:
        assert fig is None
        fig,ax = plt.subplots(1)
    cmap = definecolorbar(**kwargs)
    H, xedges, yedges = np.histogram2d(xs, ys, bins=ns,weights = weight)
    H = H.T
    X, Y = np.meshgrid(xedges, yedges)
    cms=ax.pcolormesh(X,Y, H, cmap=cmap)
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


