import numpy as np

def read_table(filename = "CI_PI_cutoff_tables.txt"):
    try:
        f=open(filename,'r')
    except:
        f=open('quasarscan/'+filename)
    lines = f.readlines()[6:]
    f.close()
    ions = []
    redshifts = []
    rhos = {}
    ts = {}
    for line in lines:
        if len(line.split())==2:
            current_ion = line.replace('\n','')
            ions.append(current_ion)
        elif len(line.split())>2:
            redshift = line.split(',')[0]
            if not redshift in redshifts:
                redshifts.append(redshift)
            rho_or_t = line.split()[1].replace(':','')
            start_of_list = line.index('[')
            if rho_or_t == 'rho':
                rhos[(current_ion,redshift)] = np.fromstring(line[start_of_list+1:-2],sep = ' ')
            elif rho_or_t == 't':
                ts[(current_ion,redshift)] = np.fromstring(line[start_of_list+1:-2],sep = ' ')
    return ions,redshifts,rhos, ts
table = read_table()

def cutoffs_for_ion_at_redshift(ion,redshift):
    if isinstance(redshift,str):
        redshift = float(redshift)
    ions,redshifts,rhos,ts = table
    if ion in ions and str(redshift) in redshifts:
        return rhos[(ion,str(redshift))],ts[(ion,str(redshift))]
    else:
        for i,z in enumerate(redshifts):
            if float(z)>redshift:
                break
            if i == len(redshifts)-1:
                return rhos[(ion,str(z))],ts[(ion,str(z))]
        z_below = redshifts[i-1]
        z_above = redshifts[i]
        fraction_z_below = (float(z_above)-redshift)/(float(z_above)-float(z_below))
        fraction_z_above = 1-fraction_z_below
        if len(rhos[(ion,z_below)]) > len(rhos[(ion,z_above)]):
            adjusted_rho_below = np.append(rhos[(ion,z_below)][:len(rhos[(ion,z_above)])-1],rhos[(ion,z_below)][-1])
            adjusted_rho_above = rhos[(ion,z_above)]
            adjusted_t_below = np.append(ts[(ion,z_below)][:len(ts[(ion,z_above)])-1],ts[(ion,z_below)][-1])
            adjusted_t_above = ts[(ion,z_above)]
        if len(rhos[(ion,z_below)]) < len(rhos[(ion,z_above)]):
            adjusted_rho_below = rhos[(ion,z_below)]
            adjusted_rho_above = np.append(rhos[(ion,z_above)][:len(rhos[(ion,z_below)])-1],rhos[(ion,z_above)][-1])
            adjusted_t_below = ts[(ion,z_below)]
            adjusted_t_above = np.append(ts[(ion,z_above)][:len(ts[(ion,z_below)])-1],ts[(ion,z_above)][-1])
        if len(rhos[(ion,z_below)]) == len(rhos[(ion,z_above)]):
            adjusted_rho_below = rhos[(ion,z_below)]
            adjusted_rho_above = rhos[(ion,z_above)]
            adjusted_t_below = ts[(ion,z_below)]
            adjusted_t_above = ts[(ion,z_above)]
        rho_final = adjusted_rho_below*fraction_z_below+adjusted_rho_above*fraction_z_above
        t_final = adjusted_t_below*fraction_z_below+adjusted_t_above*fraction_z_above
        return rho_final,t_final

def make_PI_CI_funcs(ion,redshift):
    rhos,ts = cutoffs_for_ion_at_redshift(ion,redshift)
    def PI_ion(field, data):
        #0 if CI, 1 if PI
        temps = data['gas','temperature']
        comp = np.zeros(temps.shape)
        for i in range(1,len(ts)):
            if i==1:
                mask = temps<ts[i]
                comp[mask] = np.inf
            elif i < len(ts)-1:
                mask = np.logical_and(temps<ts[i],temps>ts[i-1])
                comp[mask] = np.sqrt(rhos[i]*rhos[i-1])
            else: 
                mask = temps>ts[i]
                comp[mask] = -np.inf
        tr = yt.YTArray((data['gas','density']/mh)<comp).astype(float)
        return tr
    def CI_ion(field,data):
        PI_field_name = ('gas','PI_%s'%ion.replace(' ',''))
        return 1-(data[PI_field_name])
    return PI_ion,CI_ion

ions = table[0]
def make_funcs(ds=None,z=None,ions = ions,add_fields=False,loud = False):
    if not z and not ds:
        if loud:
            print('assuming z=1')
        z=1.0
    elif z is None:
        try:
            z = ds.current_redshift
        except:
            z = 1/ds.current_time.value.item()-1
    all_PI_funcs = {}
    all_CI_funcs = {}
    for ion in ions:
        PI_field_name = 'PI_%s'%(ion.replace(' ',''))
        CI_field_name = 'CI_%s'%(ion.replace(' ',''))
        PI_ion,CI_ion = make_PI_CI_funcs(ion,z)
        all_PI_funcs[ion]=PI_ion
        all_CI_funcs[ion]=CI_ion
        if add_fields and ds:                     
            ds.add_field(('gas',PI_field_name),
                           sampling_type="cell",
                           function=PI_ion,
                           units='')
            ds.add_field(('gas',CI_field_name),
                           sampling_type="cell",
                           function=CI_ion,
                           units='')
        elif add_fields:
            print("no ds loaded!")
    return ions,all_PI_funcs,all_CI_funcs