#!/bin/evn python
import numdifftools
import numpy as np
from scipy import interpolate
from scipy import spatial

class fkr_wrp(object):
    def __init__(self, kr):
        self.kr = kr

    def calc_free_energy(self, pt):
        return -0.6*np.log(self.kr.evaluate(pt))

def get_gradient(string, kernel):
    return_vec = np.zeros(string.shape)
    kr_wrp = fkr_wrp(kernel)
    dfun = numdifftools.Gradient(kr_wrp.calc_free_energy)
    for ind, v in enumerate(string):
        return_vec[ind] = dfun(v)[0]
    return return_vec


def finite_t_string_method(prt, kernel, start,end,N=30,
                           MaxIter=50,interp_type="linear",lxyt=None,starting_string=None):
    dt=500*(prt.tic_max[0]-prt.tic_min[0])
    d = np.shape(start)[0]
    string_h = dt
    mu = 0.00001
    kap = 0.1
    ha = 0.1
    nstep2 = 0
    r1 = np.ones((N,2))
    r1[0,:] = 0
    r1[-1,:] = 0
    alpha_eq = np.linspace(0,1,N)
    #actual string
    string = np.zeros((N,d))
    #start and end points
    string[0,:] = start
    string[-1,:] = end

    if lxyt == None:
          #start by linearly interpolatin
        f = interpolate.interp1d([0,1], [start,end],kind=interp_type ,axis=0)
        string = f(alpha_eq)
        initial_string  = np.copy(string)

        string_av = np.copy(string)
        kdtree = spatial.kdtree.KDTree(string)
    else:
        f = interpolate.interp1d(lxyt, starting_string,kind=interp_type ,axis=0)
        string = f(alpha_eq)
        initial_string  = np.copy(string)

        string_av = np.copy(string)
        ##get the triangluation of each initial string.can use triangulation but i think kdtree
        #should work too.
        kdtree = spatial.kdtree.KDTree(initial_string)

    for i in np.arange(MaxIter):
        string_copy = np.copy(string)
        grd = get_gradient(string, kernel)
        string = string -dt*grd +np.sqrt(2*string_h*mu)*np.random.normal(loc=0.0,scale=1.0,size=(N,d))

        distances,vertex_indices = kdtree.query(string)

        neighbor_based_update =  vertex_indices == np.arange(N)
        vertex_indices_copy = vertex_indices
        #either keep the string if in the same region or move it to an new region
        for i,image in enumerate(string):
            dis_mat = np.zeros(N)
            for j,original in enumerate(initial_string):
                dis_mat[j] = np.linalg.norm(original - image)
            vertex_indices_copy[i] = np.argmin(dis_mat)
        string[:,0] = string[:,0]*neighbor_based_update +string_copy[:,0]*(1-neighbor_based_update)
        string[:,1] = string[:,1]*neighbor_based_update +string_copy[:,1]*(1-neighbor_based_update)
        ##getting the new rolling averages
        string_av = (string_av*(i+nstep2)+string)/(i+nstep2+1)
        string_shift = np.roll(initial_string,-1,axis=0) + np.roll(initial_string,1,axis=0)\
        - 2*initial_string
        #move initial string towards the average
        #print string_shift
        initial_string = initial_string -(initial_string-string_av)*ha+kap*ha*N*r1*string_shift

        string_shift = initial_string - np.roll(initial_string,1,axis=0)
        string_shift[0,:] = 0

        #we need to associate every point on the string with number
        lxyt = np.cumsum(np.sqrt(string_shift[:,0] * string_shift[:,0]\
        + string_shift[:,1] * string_shift[:,1]))
        #normalize the number
        lxyt = lxyt/lxyt[-1]
        #set up a function to take in lxyt
        my_func =  interpolate.interp1d(lxyt,initial_string,kind='linear',axis=0)
        #extrapolate equally spaced points to get the initial string.
        initial_string = my_func(alpha_eq)

        kdtree = spatial.kdtree.KDTree(initial_string)
        distances,vertex_indices = kdtree.query(string)
        neighbor_based_update =  vertex_indices == np.arange(N)
        #either keep the string if in the same region or move it to an new region
        string[:,0] = string[:,0]*neighbor_based_update +\
        initial_string[:,0]*(1-neighbor_based_update)
        string[:,1] = string[:,1]*neighbor_based_update +\
        initial_string[:,1]*(1-neighbor_based_update)
    return lxyt, string ,initial_string
