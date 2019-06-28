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

def add_mass_fields(ds):
    atoms = ['C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', \
            'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', \
            'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']
    def _specific_metal_mass_function(atom):
        def _specific_metal_mass(field, data):
            nucleus_densityIa = data['gas','metal_ia_density']*\
                                SNIa_abundance[atom]
            nucleus_densityII = data['gas','metal_ii_density']*\
                                SNII_abundance[atom]
            return (nucleus_densityIa + nucleus_densityII)*data['gas','cell_volume']
        return _specific_metal_mass
    for atom in atoms:
        ds.add_field(('gas','%s_nuclei_mass'%atom),
                        sampling_type="cell",
                        function=_specific_metal_mass_function(atom),
                        units=ds.unit_system["mass"])


# based on Iwamoto et al 1999
# number of grams atom per gram of SNIa metal
SNIa_abundance = {
    'H'  : 0.00E+00, 'He' : 0.00E+00, 'C'  : 3.52E-02, 
    'N'  : 8.47E-07, 'O'  : 1.04E-01, 'F'  : 4.14E-10, 
    'Ne' : 3.30E-03, 'Na' : 4.61E-05, 'Mg' : 6.25E-03, 
    'Al' : 7.19E-04, 'Si' : 1.14E-01, 'P'  : 2.60E-04, 
    'S'  : 6.35E-02, 'Cl' : 1.27E-04, 'Ar' : 1.14E-02, 
    'K'  : 5.72E-05, 'Ca' : 8.71E-03, 'Sc' : 1.61E-07, 
    'Ti' : 2.50E-04, 'V'  : 5.46E-05, 'Cr' : 6.19E-03, 
    'Mn' : 6.47E-03, 'Fe' : 5.46E-01, 'Co' : 7.59E-04, 
    'Ni' : 9.17E-02, 'Cu' : 2.19E-06, 'Zn' : 2.06E-05}

# number of grams atom per gram of SNII metal
SNII_abundance = {
    'H'  : 0.00E+00, 'He' : 0.00E+00, 'C'  : 3.12E-02, 
    'N'  : 6.15E-04, 'O'  : 7.11E-01, 'F'  : 4.57E-10, 
    'Ne' : 9.12E-02, 'Na' : 2.56E-03, 'Mg' : 4.84E-02, 
    'Al' : 5.83E-03, 'Si' : 4.81E-02, 'P'  : 4.77E-04, 
    'S'  : 1.62E-02, 'Cl' : 4.72E-05, 'Ar' : 3.15E-03, 
    'K'  : 2.65E-05, 'Ca' : 2.31E-03, 'Sc' : 9.02E-08, 
    'Ti' : 5.18E-05, 'V'  : 3.94E-06, 'Cr' : 5.18E-04, 
    'Mn' : 1.52E-04, 'Fe' : 3.58E-02, 'Co' : 2.86E-05, 
    'Ni' : 2.35E-03, 'Cu' : 4.90E-07, 'Zn' : 7.46E-06}
        
def loadfile(name,path):
    ds=ytload(path,name)
    add_mass_fields(ds)
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


def convert_to_cbar_scale(field,min_val,max_minus_min,log = True,floor=-5):
    if log and floor is not None:
        min_val = floor
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

    intensives = {0:('gas','mass'),1:('gas','temperature'),2:('gas','density')}
    myintensives_old = np.zeros((3,length))
    for i in [0,1,2]:
        print "working on value %s..."%(intensives[i],)
        myintensives_old[i,:] = ad[intensives[i]]
        print np.sum(ad[intensives[i]])
    return myionmasses_old,myintensives_old

def split_into_bins(low_T,high_T,low_rho,high_rho,n,myintensives_old,myionmasses_old,ions):
    t_bins = np.linspace(high_T,low_T,n+1)
    rho_bins = np.linspace(low_rho,high_rho,n+1)+np.log10(1.6737236e-24)
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

def get_overall_plot(myionmasses,myintensives,t_bins,rho_bins,ions,log = True,savefigname = None,cbar = "jet",specialion = 5,
                     colorwith = ['edge','line','dot'],dsname = None,ds=None,show=True,makeoverall=True):

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
        min_val = -5#np.log10(min_val)
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

    #plt.colorbar(CS3) # using the colorbar info I got from contourf
    if makeoverall:
        plt.plot(np.arange(1,len(ions)+1),np.log10(ion_mass_in_each_state/all_O_mass),'k')
        if specialion:
            plt.plot([specialion+1],[np.log10(ion_mass_in_each_state[specialion]/all_O_mass)],'ok')
        plt.xlabel("ionization state")
        plt.ylabel("log fraction in state")
        if dsname and ds:
            Mstar = parse_metadata.get_value('star_Rvir',dsname,redshift=ds.current_redshift)
            plt.title("%s at z=%1.3f (Mstar = %1.2e)"%(dsname,ds.current_redshift,Mstar))
        if savefigname:
            plt.savefig("quasarscan/plots/%s"%savefigname)
        if show:
            plt.show()
        else:
            plt.clf()
    return min_val,max_minus_min,ion_mass_in_each_state,all_O_mass,total_gas_mass,CS3,mymap

def get_phaseplot(myionmasses,myintensives,t_bins,rho_bins,min_val,max_minus_min,\
                  ion_mass_in_each_state,all_O_mass,total_gas_mass,CS3,mymap,ions,\
                  dsname=None,ds=None,log = True,savefigname = None,specialion=5,colorwith = ['edge','line','dot'],show=True):
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
                ax = fig.add_subplot(size_t,size_rho,size_t*i+(j+1),sharex=ax0,sharey=ax0)

            edge_number = all_gas_mass_in_current_bin/total_gas_mass
            if edge_number==0 or 'edge' not in colorwith or np.log10(edge_number)<-5:
                edge_color = 'k'
                edge_width = 1
            else:
                edge_color = mymap(convert_to_cbar_scale(edge_number,min_val,max_minus_min,log=log))
                edge_width = 2
            for spine in ax.spines.values():
                spine.set_edgecolor(edge_color)
                spine.set_linewidth(edge_width)

            line_number = all_O_mass_in_current_bin/all_O_mass
            if line_number==0 or 'line' not in colorwith or np.log10(line_number)<-5:
                line_color = 'k'
                line_width = 1
                ax.tick_params(axis='y',          
                                which='both',    
                                left=False,   
                                right=False,    
                                labelleft=False)
            else:
                line_color = mymap(convert_to_cbar_scale(line_number,min_val,max_minus_min,log=log))
                line_width = 2
            ax.plot(np.arange(1,len(ions)+1), np.log10(ion_mass_in_current_bin/all_O_mass_in_current_bin), color=line_color, linewidth=line_width)

            if specialion is not None:
                dot_number = ion_mass_in_current_bin[specialion]/ion_mass_in_each_state[specialion]
                if dot_number==0 or 'dot' not in colorwith or np.log10(dot_number)<-5:
                    dot_color = 'k'
                    dot_size = 3.
                else:
                    dot_color = mymap(convert_to_cbar_scale(dot_number,min_val,max_minus_min,log=log))
                    dot_size = 7.
                ax.plot([specialion+1],[np.log10(ion_mass_in_current_bin[specialion]/all_O_mass_in_current_bin)],'o',mec = 'k',mfc=dot_color,ms=dot_size)
                ax.set_ylim(-5,0)
            
            if len(t_bins)>9:
                ax.tick_params(axis='x',          
                                which='both',    
                                bottom=False,   
                                top=False,    
                                labelbottom=False)
            if len(t_bins)>11:
                ax.tick_params(axis='y',          
                                which='both',    
                                left=False,   
                                right=False,    
                                labelleft=False)
            if i==0:
                if len(t_bins)>11:
                    ax.set_title('%.1f'%((rho_low-np.log10(1.6737236e-24)+rho_high-np.log10(1.6737236e-24))/2))
                else:
                    ax.set_title('[%.2f,%.2f]'%(rho_low-np.log10(1.6737236e-24),rho_high-np.log10(1.6737236e-24)))
            if j==0:
                if len(t_bins)>11:
                    ax.set_ylabel('%.1f'%((t_low+t_high)/2))
                else:
                    ax.set_ylabel("[%.2f,%.2f]"%(t_low,t_high))
    if savefigname:
        plt.savefig("quasarscan/plots/%s"%savefigname)
    if show:
        plt.show()
    else:
        plt.clf()

def make_ionmasses_and_intensives(name,path,atom = 'O',CGM=True):
    ds = loadfile(name,path)
    ad,ions = addions_get_ad(atom,ds,name,CGM)
    myionmasses_old,myintensives_old = get_original_lists(ad,ions)
    return myionmasses_old,myintensives_old,ds,ions
    
def script(name,path,low_T,high_T,low_rho,high_rho,n,atom = 'O',saveoverallname=None,savephaseplotname=None,\
           CGM=True,myionmasses_old=None,myintensives_old=None,ds = None,log=True,cbar='jet',specialion=5,colorwith = ['edge','line','dot'],\
          show=True, makeoverall=True):
    if myionmasses_old is None and myintensives_old is None and ds is None:
        myionmasses_old,myintensives_old,ds,ions = make_ionmasses_and_intensives(name,path,atom,CGM)
    else:
        ions = make_full_ionlist(atom)
    try:
        myionmasses,myintensives,t_bins,rho_bins = split_into_bins(low_T,high_T,low_rho,high_rho,n,myintensives_old,myionmasses_old,ions)
        min_val,max_minus_min,ion_mass_in_each_state,all_O_mass,total_gas_mass,CS3,mymap = get_overall_plot(myionmasses,myintensives,t_bins,\
                                                                                                            rho_bins,ions,log,saveoverallname,cbar,specialion=specialion,\
                                                                                                               colorwith=colorwith,dsname=name,ds=ds,\
                                                                                                            show=show,makeoverall=makeoverall)
        get_phaseplot(myionmasses,myintensives,t_bins,rho_bins,min_val,max_minus_min,\
                          ion_mass_in_each_state,all_O_mass,total_gas_mass,CS3,mymap,ions,name,ds,log,savephaseplotname,specialion=specialion,colorwith=colorwith,show=show)
    except:
        pass
    return myionmasses_old,myintensives_old,ds