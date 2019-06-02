#!/usr/bin/env python

import numpy as np
import yt
import sys

try:
    from quasarscan import parse_metadata
    from quasarscan import quasar_sphere
    from quasarscan import ion_lists
    from quasarscan import gasbinning
    from quasarscan import code_specific_setup

except:
    import parse_metadata
    import quasar_sphere
    import ion_lists
    import gasbinning
    import code_specific_setup

def convert_to_xyz(r, theta, phi):
    return np.array([r*np.sin(theta)*np.cos(phi),r*np.sin(theta)*np.sin(phi),r*np.cos(theta)])


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians. (taken from stack exchange question)
    """
    axis = np.asarray(axis)
    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2.0)
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def get_rotation_matrix(L):
    zhat = np.array([0,0,1])
    if np.array_equal(L,zhat):
        return np.diag([1,1,1])
    theta = np.arccos(L[2]/np.linalg.norm(L))
    axis = np.cross(zhat,L)
    return rotation_matrix(axis,theta)
        

def ray_endpoints_spherical(R,r,theta,phi,alpha,endonsph):
    start = convert_to_xyz(R,theta,phi)
    xhat = convert_to_xyz(1,np.pi/2,np.pi/2+phi)
    yhat = convert_to_xyz(1,np.pi/2-theta,np.pi+phi)
    mid = r*(np.cos(alpha)*xhat+np.sin(alpha)*yhat)
    diff = start-mid
    if endonsph:
        t = 2*np.dot(start,diff)/np.dot(diff,diff)
    else:
        t = 2*R/np.linalg.norm(diff)
    end = start*(1-t)+mid*t
    return np.array([start,end])

def weights(array,function):
    if function == "sin":
        probs = np.sin(array)/2
        probs[0] = probs[-1]
    elif function == "lin":
        probs = np.linspace(0,1,len(array)+1)[1:]
    probs /= np.sum(probs)
    return probs

def create_QSO_endpoints(sphere, ions,code_unit=None,gasbins=None,include_0 = True,\
                         L = None, center=None, endonsph = False, dsbounds = None):
    R=sphere[0]
    n_th=sphere[1]
    n_phi=sphere[2]
    n_r=sphere[3]
    rmax=sphere[4]
    length=sphere[5]
    if dsbounds is not None and np.any(2*R>dsbounds):
        ratio = float(rmax) / R
        #this is the minimum ratio needed to ensure a sightline within a box of edge 2*R
        #stays within the box
        assert ratio < .57
        R = np.min(dsbounds-center.in_units('kpc').value)
        rmax = R*ratio
        print "R > dsbounds: adjusting to R=%2.2f kpc, rmax = %2.2f kpc"%(R,rmax)
    r_arr = np.linspace(0,rmax,n_r)
    if not include_0:
        r_arr = r_arr[1:]
    th_arr = np.linspace(0,np.pi,n_th,endpoint = False)
    phi_arr = np.linspace(0,2*np.pi,n_phi,endpoint = False)

    scanparams = [None]*7
    scanparams[0] = R
    scanparams[1] = len(th_arr)
    scanparams[2] = len(phi_arr)
    scanparams[3] = len(r_arr)
    scanparams[4] = rmax
    scanparams[5] = length
    scanparams[6] = 0
    weightth = weights(th_arr, "sin")
    weightr = weights(r_arr, "lin")        
    num_extra_columns = gasbins.get_length()
    info = np.zeros((int(length),11+len(ions)*(num_extra_columns+2)+3))-1.0
    if L is None:
        L = np.array([0,0,1.])
    if center is None:
        center = yt.YTArray([0,0,0],'kpc')
        code_unit = 'kpc'
    rot_matrix = get_rotation_matrix(L)
    for i in range(int(length)):
        theta = np.random.choice(th_arr,p = weightth)
        r = np.random.choice(r_arr,p = weightr)
        phi= np.random.choice(phi_arr)
        alpha = 2*np.pi*np.random.random()
        info[i][:5] = np.array([i,theta,phi,r,alpha])
        info[i][5:8] = (yt.YTArray(np.matmul(rot_matrix, ray_endpoints_spherical(R,r,theta,phi,alpha,endonsph)[0]),'kpc') + center).in_units(code_unit).value
        info[i][8:11] = (yt.YTArray(np.matmul(rot_matrix, ray_endpoints_spherical(R,r,theta,phi,alpha,endonsph)[1]),'kpc') + center).in_units(code_unit).value
    return scanparams,info

if __name__ == "__main__":
    name = sys.argv[1]
    path = sys.argv[2]
    if len(sys.argv)==4:
        ionlist = sys.argv[3]
    else:
        ionlist = None
    code = name.split("_")[2]
    ds = code_specific_setup.ytload(path,code)
    try:
        z = ds.current_redshift
    except:
        a = ds.current_time.value
        z = 1./a-1.
    Rvir = parse_metadata.get_value("Rvir",name,z)
    if np.isnan(Rvir):
        Rvir = 100#kpc
    center = ds.find_max(('gas','density'))[1]
    L = parse_metadata.get_value("L",name,z)
    if np.isnan(L).all():
        L = np.array([0,0,1.])
    convert_unit = ds.length_unit.units
    include_0 = True
    if ds.dataset_type == 'tipsy':
        include_0 = False
    convert = ds.length_unit.in_units('kpc').value.item()
    defaultsphere = 6*Rvir,12,12,12,1.5*Rvir,448
    testsphere = 6*Rvir,12,12,12,2*Rvir,10
    dsbounds = ds.domain_width.to('kpc').value
    defaultions = ion_lists.agoraions
    if ionlist:
        try:
            ions = ion_lists.dict_of_ionlists[ionlist]
        except:
            ions = ionlist.replace("[","").replace("]","").split(",")
    else:
        ions = defaultions
    gasbins = gasbinning.GasBinsHolder(code_specific_setup.get_gasbins_arg(code))
    scanparams, info = create_QSO_endpoints(defaultsphere,ions,code_unit = convert_unit,L=L,center=center,gasbins = gasbins,dsbounds=dsbounds,include_0 = include_0)
    center = center.in_units(convert_unit).value
    simparams = [name,z,center[0],center[1],center[2],Rvir,path,L[0],L[1],L[2],convert, str(convert_unit)]
    q = quasar_sphere.QuasarSphere(simparams=simparams,scanparams= scanparams,ions= ions,data=info,gasbins= gasbins)
    q.save_values(at_level = 0)

