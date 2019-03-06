#!/usr/bin/env python

import numpy as np
import yt
import sys

try:
    from quasarscan import parse_metadata
    from quasarscan import quasar_sphere
    from quasarscan import ion_lists
    from quasarscan import gasbinning

except:
    import parse_metadata
    import quasar_sphere
    import ion_lists
    import gasbinning

def convert_to_xyz(r, theta, phi):
    return np.array([r*np.sin(theta)*np.cos(phi),r*np.sin(theta)*np.sin(phi),r*np.cos(theta)])


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
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

def create_QSO_endpoints(sphere, convert_code_unit_to_kpc,ions,gasbins=None,\
                         L = None, center=None, endonsph = False):
    R=sphere[0]
    n_th=sphere[1]
    n_phi=sphere[2]
    n_r=sphere[3]
    rmax=sphere[4]
    length=sphere[5]
    r_arr = np.linspace(0,rmax,n_r)
    th_arr = np.linspace(0,np.pi,n_th,endpoint = False)
    phi_arr = np.linspace(0,2*np.pi,n_phi,endpoint = False)

    R /= convert_code_unit_to_kpc
    r_arr /= convert_code_unit_to_kpc
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
        center = np.array([0,0,0.])
    rot_matrix = get_rotation_matrix(L)
    for i in range(int(length)):
        theta = np.random.choice(th_arr,p = weightth)
        r = np.random.choice(r_arr,p = weightr)
        phi= np.random.choice(phi_arr)
        alpha = 2*np.pi*np.random.random()
        info[i][:5] = np.array([i,theta,phi,r,alpha])
        info[i][5:8] = np.matmul(rot_matrix, ray_endpoints_spherical(R,r,theta,phi,alpha,endonsph)[0]) + center
        info[i][8:11] = np.matmul(rot_matrix, ray_endpoints_spherical(R,r,theta,phi,alpha,endonsph)[1]) + center 
    return scanparams,info

if __name__ == "__main__":
    name = sys.argv[1]
    path = sys.argv[2]
    if 'art' in name:
        h,d,s = quasar_sphere.get_aux_files_art(path)
        ds = yt.load(path,file_particle_header=h,\
                                  file_particle_data=d,\
                                  file_particle_stars=s)
    else:
        ds = yt.load(path)
    z = ds.current_redshift
    Rvir = parse_metadata.get_value("Rvir",name,z)
    if Rvir is None:
        Rvir = 100#kpc
    center = ds.find_max(('gas','density'))[1]
    L = parse_metadata.get_value("L",name,z)
    if L is None:
        L = np.array([0,0,1.])
    convert = float(ds.length_unit.in_units('kpc').value)
    defaultsphere = 6*Rvir,12,12,12,2*Rvir,448
    testsphere = 6*Rvir,12,12,12,2*Rvir,10
    defaultions = ion_lists.agoraions
    gasbins = gasbinning.GasBinsHolder("all")
    scanparams, info = create_QSO_endpoints(defaultsphere,convert,defaultions,L=L,center=center,gasbins = gasbins)
    simparams = [name,z,center[0],center[1],center[2],Rvir,path,L[0],L[1],L[2],convert]
    q = quasar_sphere.QuasarSphere(simparams=simparams,scanparams= scanparams,ions= defaultions,data=info,gasbins= gasbins)
    q.save_values(at_level = 0)

