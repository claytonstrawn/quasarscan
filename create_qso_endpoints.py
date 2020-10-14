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

#summary: change 3D spherical coordinates to cartesian
#
#inputs: r, theta, phi: radial distance, polar angle, azimuthal angle
#
#outputs: x,y,z in a single array
def convert_to_xyz(r, theta, phi):
    return np.array([r*np.sin(theta)*np.cos(phi),r*np.sin(theta)*np.sin(phi),r*np.cos(theta)])

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
        print("R > dsbounds: adjusting to R=%2.2f kpc, rmax = %2.2f kpc"%(R,rmax))
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

#when running this as a script, it will use the simulation name 
#and path to the simulation snapshot to load it. It reads some information
#(redshift, center) from yt's snapshot and gets other information (Rvir,L)
#from pre-loaded data tables through 'parse-metadata'. Everything else is 
#simply a default, change if you want to collect different information
if __name__ == "__main__":
    name = sys.argv[1]
    path = sys.argv[2]
    if len(sys.argv)>=4:
        ionlist = sys.argv[3]
    else:
        ionlist = None
    if 'test' in sys.argv:
        test = True
    else:
        test = False
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
    if name == 'VELA_v2_art_29' and np.abs(z-1)<.1:
        print('Using VELA29 special center....')
        center = ds.arr([0.49281657, 0.49691165, 0.49609502], 'code_length')
    else:
        center = ds.find_max(('gas','density'))[1]
    L = parse_metadata.get_value("L",name,z)
    if np.isnan(L).all():
        L = np.array([0,0,1.])
    convert_unit = ds.length_unit.units
    include_0 = True
    if ds.dataset_type == 'tipsy':
        include_0 = False
    convert = ds.length_unit.in_units('kpc').value.item()
    defaultsphere = 6*Rvir,12,12,12,2*Rvir,384
    if test:
        defaultsphere = .1*Rvir,12,12,12,.01*Rvir,10
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

