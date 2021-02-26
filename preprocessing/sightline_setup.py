#skeleton for sightline_setup
#designed to create and save QuasarSphere with realistic sightlines but with no data filled in
#relies on yt (for making unit conversions)
#relies on quasarsphere (makes QuasarSphere object)
#relies on code_specific_setup (to load any kind of code, though with our 
#              metadata in a prescribed format this doesn't matter as much)

#!/usr/bin/env python

import numpy as np
import sys
from quasarscan3.preprocessing import parse_metadata
from quasarscan3.data_objects import quasar_sphere,simulation_quasar_sphere,gasbinning
from quasarscan3.utils import ion_lists,utils

class BadMetadataError(Exception):
    def __init__(self, message):
        self.message = message
        print(self.message)
class BadListError(Exception):
    def __init__(self, message):
        self.message = message
        print(self.message)

#summary: Return the rotation matrix associated with counterclockwise rotation about
#         the given axis by theta radians. (taken from stack exchange question)
#
#inputs: axis: 3d vector for the axis around which to rotate
#        theta: angle in radians to rotate around this axis
#
#outputs: rotation matrix to rotate all vectors by this angle around this axis
def rotation_matrix(axis, theta):
    axis = np.asarray(axis)
    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2.0)
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

#summary: Return the rotation matrix which identifies the z axis
#         with the angular momentum of the simulation L
#
#inputs: L: 3d vector for the angular momentum of simulation
#
#outputs: rotation matrix to rotate all vectors by to make them 
#         'face-on' = theta<pi/2 or theta>3pi/2
#         'edge-on' = pi/2<theta<3pi/2
def get_rotation_matrix(L):
    zhat = np.array([0,0,1])
    if np.array_equal(L,zhat):
        return np.diag([1,1,1])
    theta = np.arccos(L[2]/np.linalg.norm(L))
    axis = np.cross(zhat,L)
    return rotation_matrix(axis,theta)
        
#summary: given the position in terms of 5 parameters (and a 
#         optional requirement to end the sightline on the sphere 
#         on the other side), return the xyz coordinates of the two
#         endpoints
#
#inputs: R: "outer" radius to use, the origin of the sightline. If 
#           too far, resolution will be bad (dep. on simulation). 
#           If too close, geometry will be bad. In Strawn et al. 2020
#           this was 6Rvir
#        r: impact parameter of sightline compared to galaxy center
#        theta: polar angle of starting point
#        phi: azimuthal angle of starting point
#        alpha: angle of rotation along a circle of radius r. (circle 
#               is defined with x' direction horizontal w.r.t. angular momentum
#               y' at a right angle to x' and at right angle to the line from
#               galaxy center to starting point
#        endonsph: If True, require the sightlines to stop when they hit the sphere.
#                  if False, require the sightlines to have a length of 2x the
#                  distance from starting point to midpoint
#
#outputs: array: 2 3-d arrays of the xyz coordinates of the starting and ending
#                points of the sightline
def ray_endpoints_spherical(R,r,theta,phi,alpha,endonsph):

    start = R*np.array([np.sin(theta)*np.cos(phi),
                      np.sin(theta)*np.sin(phi),
                      np.cos(theta)])

    xhat = np.array([-np.cos(theta)*np.cos(phi),
                     -np.cos(theta)*np.sin(phi),
                     np.sin(theta)])

    yhat = np.array([np.sin(phi),
                     -np.cos(phi),
                     0.])

    mid = r*(np.cos(alpha)*xhat+np.sin(alpha)*yhat)
    diff = start-mid
    if endonsph:
        t = 2*np.dot(start,diff)/np.dot(diff,diff)
    else:
        t = 2*R/np.linalg.norm(diff)
    end = start*(1-t)+mid*t
    return np.array([start,end])

#summary: to evenly distribute startpoints on surface of outer sphere, and 
#         midpoints in 3D sphere of radius "max-r" a weighting will be 
#         added to random choice. NOTE: since this will naively reject all chance
#         of selection of "0" in both cases, both probabilities are adjusted to be
#         nonzero to allow theta=0 and r=0 a chance of selection
#
#inputs: array: list of possible points to be chosen
#        function: if "sin" assign probabilities according to sin of array [theta]
#                  if "lin" assign probabilities linearly (later in list = more likely) [r]
#
#outputs: probs: probability of selection of each point.
def weights(array,function):
    if function == "sin":
        probs = np.sin(array)/2
        probs[0] = probs[-1]
    elif function == "lin":
        probs = np.linspace(0,1,len(array)+1)[1:]
    probs /= np.sum(probs)
    return probs

#summary: Create the raw data array for this simulation, each sightline is one 
#         line of the array
# 
#inputs: sphere: tuple of (R,n_th,n_phi,n_r,rmax, and length), where R is 
#                outer sphere radius, n_th is number of possible polar angle 
#                values between 0 and pi, n_phi is the number of possible 
#                azimuthal angle values between 0 and 2pi, n_r is the number
#                of possible impact parameters between 0 and rmax, rmax is the 
#                maximum impact parameter probed by sightline, length is the 
#                number of such sightlines to generate
#        ions: which ions will be tracked
#        code_unit: if the requires a special unit (i.e. "code_length") which 
#                   needs to be understood in advance, add it
#        gasbins: what extra information will be stored in addition to column 
#                 densities of the different ions (see quasarscan.gasbinning)
#        include_0: whether to allow for 0 impact parameters. Default True
#        L: angular momentum axis of galaxy. if None assume L = z axis
#        center: center of galaxy. if None assume center = 0,0,0
#        endonsph: If True, require the sightlines to stop when they hit the sphere.
#                  if False, require the sightlines to have a length of 2x the
#                  distance from starting point to midpoint
#        dsbounds: furthest allowed distance from center in any direction (otherwise 
#                  you would leave the simulaton)
#
#outputs: scanparams: a bunch of information about the simulation (mostly the 
#                     arguments to this function)
#         info: 2d array of one line = one sightline, each line is 11 setup numbers
#               and len(ions)*len(gasbins) '-1's, which will be filled in in 
#               'get_coldens.py'
def create_QSO_endpoints_helper(R, n_th,n_phi,n_r,rmax,length, ions,gasbins,L, center,\
                         Rvir,include_0=True,endonsph=False):
    assert isinstance(R,tuple) and isinstance(rmax,tuple) and R[1] in ['Rvir','kpc'] and rmax[1] in ['Rvir','kpc'], \
            "Please provide R and rmax as tuples with float size and unit 'Rvir' or 'kpc'"

    if R[1] == 'Rvir':
        R = R[0]*Rvir
    elif R[1] == 'kpc':
        R = R[0]
    if rmax[1] == 'Rvir':
        rmax = rmax[0]*Rvir
    elif rmax[1] == 'kpc':
        rmax = rmax[0]
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
    rot_matrix = get_rotation_matrix(L)
    for i in range(int(length)):
        theta = np.random.choice(th_arr,p = weightth)
        r = np.random.choice(r_arr,p = weightr)
        phi= np.random.choice(phi_arr)
        alpha = 2*np.pi*np.random.random()
        info[i][:5] = np.array([i,theta,phi,r,alpha])
        start_needrotation, end_needrotation = ray_endpoints_spherical(R,r,theta,phi,alpha,endonsph)
        info[i][5:8] = np.matmul(rot_matrix, start_needrotation) + center
        info[i][8:11] = np.matmul(rot_matrix, end_needrotation) + center
    return scanparams,info

def create_QSO_endpoints(fullname,redshift,ions,gasbins='default',R=(6,'Rvir'),\
                        n_th=12,n_phi=12,n_r=12,rmax=(2,'Rvir'),length=384,\
                        run='default',filename=None,**kwargs):
    code = fullname.split('_')[2]
    if isinstance(ions,str):
        try:
            ions = ion_lists.dict_of_ionlists[ions]
        except KeyError:
            raise BadListError('list nickname not recognized, enter ions as list of strings')
    Rvir = parse_metadata.get_value("Rvir",fullname,redshift=redshift)
    center_x = parse_metadata.get_value("center_x",fullname,redshift=redshift)
    center_y = parse_metadata.get_value("center_y",fullname,redshift=redshift)
    center_z = parse_metadata.get_value("center_z",fullname,redshift=redshift)
    center = np.array([center_x,center_y,center_z])
    L_x = parse_metadata.get_value("L_x",fullname,redshift=redshift)
    L_y = parse_metadata.get_value("L_y",fullname,redshift=redshift)
    L_z = parse_metadata.get_value("L_z",fullname,redshift=redshift)
    L = np.array([L_x,L_y,L_z])
    if np.isnan(Rvir) or np.any(np.isnan(center)) or np.any(np.isnan(L)):
        if run != 'force':
            raise BadMetadataError("metadata file for %s is missing one of Rvir, center, L")
        elif run == 'force':
            Rvir = 100 #assume 100 kpc
            L = np.array([0,0,1])
            assert filename is not None, "Need source file to force center finding"
            # these dependencies are only present if you will be forcing this, 
            # it is recommended that you run 'create_metadata' first
            import yt
            from quasarscan3.preprocessing import code_specific_setup
            ds = code_specific_setup.ytload(filename,code)
            center = ds.find_max('gas','density')[1].in_units('kpc').value

    gasbins = gasbinning.GasBinsHolder(utils.get_gasbins_arg(code))
    scanparams,info = create_QSO_endpoints_helper(R, n_th,n_phi,n_r,rmax,length, ions,gasbins,L, center,Rvir)
    simparams = [fullname,redshift,center[0],center[1],center[2],Rvir,filename,L[0],L[1],L[2]]
    q = simulation_quasar_sphere.SimQuasarSphere((simparams,scanparams,ions,info,gasbins))
    if run == 'test':
        return q
    outputfilename = q.save_values()
    return outputfilename

