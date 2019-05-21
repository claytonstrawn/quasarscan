import yt
import trident
import numpy as np
import matplotlib.pyplot as plt
from quasarscan import quasar_sphere
from quasarscan import parse_metadata
from quasarscan.quasar_sphere import ion_to_field_name
from quasarscan.code_specific_setup import ytload
from trident.ion_balance import atomic_number
from quasarscan.roman import to_roman

def loadfile(name,path):
    ds=ytload(path,name)
    return ds

def make_full_ionlist(atom):
    num_states = atomic_number[atom]+1
    mylist = []
    for i in range(1,num_states+1):
        mylist.append(atom + ' ' + to_roman(i))
    return mylist

def addions_get_ad(atom,ds,name,CGM = True):
    ions = make_full_ionlist(atom)
    trident.add_ion_fields(ds,ions)
    Rvir = parse_metadata.get_value('Rvir',name,redshift=ds.current_redshift)
    if CGM:
        ad = ds.sphere("c", (Rvir, "kpc"))-ds.sphere("c", (.1*Rvir, "kpc"))
    else:
        ad = ds.all_data()
    return ad,ions

def convert_to_cbar_scale(field,min_val,max_minus_min,log = True):
    if log:
        return (np.log10(field)-min_val)/max_minus_min
    else:
        return (field-min_val)/max_minus_min

def get_original_lists(ad,ions):
    length = len(ad[('gas',quasar_sphere.ion_to_field_name(ions[0],"mass"))])
    myionmasses_old = np.zeros((len(ions),length))
    for i,ion in enumerate(ions):
        print "working on ion %s..."%(ion,)
        myionmasses_old[i,:] = ad[('gas',ion_to_field_name(ion,"mass"))]
        print np.sum(myionmasses_old[i])

    intensives = {0:('gas','mass'),1:('gas','temperature'),2:('gas','density'),3:('gas','cell_volume')}
    myintensives_old = np.zeros((4,length))
    for i in [0,1,2,3]:
        print "working on value %s..."%(intensives[i],)
        myintensives_old[i,:] = ad[intensives[i]]
        print np.sum(ad[intensives[i]])
    return myionmasses_old,myintensives_old

def split_into_bins(low_T,high_T,low_rho,high_rho,n_Tbins,n_rhobins,myintensives_old,myionmasses_old,ions):
    t_bins = np.linspace(high_T,low_T,n_Tbins+1)
    rho_bins = np.linspace(low_rho,high_rho,n_rhobins+1)+np.log10(1.6737236e-24)
    mask_lowbound_rho = np.log10(myintensives_old[2])>np.min(rho_bins)
    mask_highbound_rho = np.log10(myintensives_old[2])<np.max(rho_bins)
    mask_lowbound_T = np.log10(myintensives_old[1])>np.min(t_bins)
    mask_highbound_T = np.log10(myintensives_old[1])<np.max(t_bins)
    mask = np.logical_and(mask_lowbound_rho,mask_highbound_rho)
    mask = np.logical_and(mask,mask_lowbound_T)
    mask = np.logical_and(mask,mask_highbound_T)
    newlen =len(myionmasses_old[0][mask])
    myionmasses = np.zeros((len(ions),newlen))
    myintensives = np.zeros((len(ions),newlen))
    for i in range(9):
        myionmasses[i] = myionmasses_old[i][mask]
    for i in range(3):
        myintensives[i] = myintensives_old[i][mask]
    return myionmasses,myintensives,t_bins,rho_bins

def get_overall_plot(myionmasses,myintensives,t_bins,rho_bins,ions,log = True,savefigname = None,cbar = "jet",specialion = 5,colorwith = ['edge','line','dot']):
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    max_val = 1e-10
    min_val = 1e10

    ion_mass_in_each_state = np.sum(myionmasses,axis=1)
    all_O_mass = np.sum(ion_mass_in_each_state)
    total_gas_mass = np.sum(myintensives[0])

    for i,t in enumerate(t_bins[:-1]):
        t_high = t
        t_low = t_bins[i+1]
        t_mask = np.logical_and(np.log10(myintensives[1])>t_low,np.log10(myintensives[1])<t_high)
        for j,rho in enumerate(rho_bins[:-1]):
            rho_low = rho
            rho_high = rho_bins[j+1]
            rho_mask = np.logical_and(np.log10(myintensives[2])>rho_low,np.log10(myintensives[2])<rho_high)
            my_mask = np.logical_and(t_mask,rho_mask)

            ion_mass_in_current_bin = np.sum(myionmasses[:,my_mask],axis=1)
            all_O_mass_in_current_bin = np.sum(ion_mass_in_current_bin)
            all_gas_mass_in_current_bin = np.sum(myintensives[0][my_mask])
            to_compare = []
            if 'edge' in colorwith:
                to_compare.append(all_gas_mass_in_current_bin/total_gas_mass)
            if 'line' in colorwith:
                to_compare.append(all_O_mass_in_current_bin/all_O_mass)
            if 'dot' in colorwith and specialion is not None:
                to_compare.append(ion_mass_in_current_bin[specialion]/ion_mass_in_each_state[specialion])
            myarray = np.array(to_compare)
            if log:
                myarray = myarray[myarray>1e-8]
            max_val = np.max(np.append(myarray,max_val))
            min_val = np.min(np.append(myarray,min_val))
    print min_val,max_val
    if log:
        min_val = np.log10(min_val)
        max_val = np.log10(max_val)
    max_minus_min = max_val-min_val
    min_cb, max_cb = (min_val, max_val)
    if len(colorwith)==0:
        min_cb, max_cb = 0.,1.
    step = 100

    # Setting up a colormap that's a simple transtion
    mymap = eval("mpl.cm.%s"%cbar)

    # Using contourf to provide my colorbar info, then clearing the figure
    Z = [[0,0],[0,0]]
    levels = np.linspace(min_cb,max_cb,step)
    CS3 = plt.contourf(Z, levels, cmap=mymap)
    plt.clf()

    plt.colorbar(CS3) # using the colorbar info I got from contourf

    plt.plot(np.arange(1,len(ions)+1),np.log10(ion_mass_in_each_state/all_O_mass),'k')
    if specialion:
        plt.plot([specialion+1],[np.log10(ion_mass_in_each_state[specialion]/all_O_mass)],'ok')
    plt.xlabel("ionization state")
    plt.ylabel("log fraction in state")
    if savefigname:
        plt.savefig("quasarscan/plots/%s"%savefigname)
    plt.show()
    return min_val,max_minus_min,ion_mass_in_each_state,all_O_mass,total_gas_mass,CS3,mymap

def get_phaseplot(myionmasses,myintensives,t_bins,rho_bins,min_val,max_minus_min,\
                  ion_mass_in_each_state,all_O_mass,total_gas_mass,CS3,mymap,ions,\
                  dsname=None,ds=None,log = True,savefigname = None,specialion=5,colorwith = ['edge','line','dot']):
    import matplotlib.pylab as pl
    fig = plt.figure(1,figsize=(15,11))
    if dsname:
        Mstar = parse_metadata.get_value('star_Rvir',dsname,redshift=ds.current_redshift)
        fig.suptitle("%s at z=%1.3f (Mstar = %1.2e)"%(dsname,ds.current_redshift,Mstar))
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.9, 0.15, 0.02, 0.7])
    fig.colorbar(CS3, cax=cbar_ax)
    size_t = len(t_bins)-1
    size_rho = len(rho_bins)-1
    for i,t in enumerate(t_bins[:-1]):
        t_high = t
        t_low = t_bins[i+1]
        t_mask = np.logical_and(np.log10(myintensives[1])>t_low,np.log10(myintensives[1])<t_high)
        for j,rho in enumerate(rho_bins[:-1]):
            rho_low = rho
            rho_high = rho_bins[j+1]
            rho_mask = np.logical_and(np.log10(myintensives[2])>rho_low,np.log10(myintensives[2])<rho_high)
            my_mask = np.logical_and(t_mask,rho_mask)

            ion_mass_in_current_bin = np.sum(myionmasses[:,my_mask],axis=1)
            all_O_mass_in_current_bin = np.sum(ion_mass_in_current_bin)
            all_gas_mass_in_current_bin = np.sum(myintensives[0][my_mask])
            if i == 0 and j == 0:
                ax0 = fig.add_subplot(size_t,size_rho,size_t*i+(j+1))
                ax = ax0
            else:
                ax = fig.add_subplot(size_t,size_rho,size_t*i+(j+1),sharex=ax0)

            edge_number = all_gas_mass_in_current_bin/total_gas_mass
            if edge_number==0 or 'edge' not in colorwith:
                edge_color = 'k'
                edge_width = 1
            else:
                edge_color = mymap(convert_to_cbar_scale(edge_number,min_val,max_minus_min,log=log))
                edge_width = 2
            for spine in ax.spines.values():
                spine.set_edgecolor(edge_color)
                spine.set_linewidth(edge_width)

            line_number = all_O_mass_in_current_bin/all_O_mass
            if line_number==0 or 'line' not in colorwith:
                line_color = 'k'
                line_width = 1
            else:
                line_color = mymap(convert_to_cbar_scale(line_number,min_val,max_minus_min,log=log))
                line_width = 2
            ax.plot(np.arange(1,len(ions)+1),np.log10(ion_mass_in_current_bin/all_O_mass_in_current_bin),color=line_color,linewidth=line_width)

            if specialion is not None:
                dot_number = ion_mass_in_current_bin[specialion]/ion_mass_in_each_state[specialion]
                if dot_number==0 or 'dot' not in colorwith:
                    dot_color = 'k'
                    dot_size = 3.
                else:
                    dot_color = mymap(convert_to_cbar_scale(dot_number,min_val,max_minus_min,log=log))
                    dot_size = 7.
                ax.plot([specialion+1],[np.log10(ion_mass_in_current_bin[specialion]/all_O_mass_in_current_bin)],'o',mec = 'k',mfc=dot_color,ms=dot_size)

            if i==0:
                ax.set_title('[%.2f,%.2f]'%(rho_low-np.log10(1.6737236e-24),rho_high-np.log10(1.6737236e-24)))
            if j==0:
                ax.set_ylabel("[%.2f,%.2f]"%(t_low,t_high))

    if savefigname:
        plt.savefig("quasarscan/plots/%s"%savefigname)
    plt.show()

def make_ionmasses_and_intensives(name,path,atom = 'O',CGM=True):
    ds = loadfile(name,path)
    ad,ions = addions_get_ad(atom,ds,name,CGM)
    myionmasses_old,myintensives_old = get_original_lists(ad,ions)
    return myionmasses_old,myintensives_old,ds,ions
    
def script(name,path,low_T,high_T,low_rho,high_rho,n_Tbins,n_rhobins,atom = 'O',saveoverallname=None,savephaseplotname=None,\
           CGM=True,myionmasses_old=None,myintensives_old=None,ds = None,log=True,cbar='jet',specialion=5,colorwith = ['edge','line','dot']):
    if myionmasses_old is None and myintensives_old is None and ds is None:
        myionmasses_old,myintensives_old,ds,ions = make_ionmasses_and_intensives(name,path,atom,CGM)
    else:
        ions = make_full_ionlist(atom)
    try:
        myionmasses,myintensives,t_bins,rho_bins = split_into_bins(low_T,high_T,low_rho,high_rho,n_Tbins,n_rhobins,myintensives_old,myionmasses_old,ions)
        min_val,max_minus_min,ion_mass_in_each_state,all_O_mass,total_gas_mass,CS3,mymap = get_overall_plot(myionmasses,myintensives,t_bins,\
                                                                                                            rho_bins,ions,log,saveoverallname,cbar,specialion=specialion,\
                                                                                                               colorwith=colorwith)
        get_phaseplot(myionmasses,myintensives,t_bins,rho_bins,min_val,max_minus_min,\
                          ion_mass_in_each_state,all_O_mass,total_gas_mass,CS3,mymap,ions,name,ds,log,savephaseplotname,specialion=specialion,colorwith=colorwith)
    except:
        pass
    return myionmasses_old,myintensives_old,ds