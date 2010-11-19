# Author: Anand Patil
# Date: 6 Feb 2009
# License: Creative Commons BY-NC-SA
####################################

import pymc as pm
import numpy as np
import os
from copy import copy
# from scipy import interpolate as interp
from fst_cov_fun import gtf
from pymc.gp.cov_funs import imul, symmetrize, nsst
from pymc.gp.cov_funs import aniso_geo_rad, euclidean, default_h
from pymc import get_threadpool_size, map_noreturn
#import MAPdata

__all__ = ['nonstationary_spatiotemporal', 'gtf']

t_gam_fun = gtf

# TODO: Do this using the thread pool. There should be a version of the code around that does.

def nonstationary_spatiotemporal(x,y,amp,scale,diff_degree,h=default_h,symm=None,geometry='aniso_geo_rad',**kwds):
    """
    Spatiotemporal covariance function. Converts x and y
    to a matrix of covariances. x and y are assumed to have
    columns (long,lat,t). Parameters are:
    - t_gam_fun: A function returning a matrix of variogram values.
      Inputs will be the 't' columns of x and y, as well as kwds.
    - amp: The MS amplitude of realizations.
    - scale: Scales distance.
    - diff_degree: A function that returns local degree of differentiability at x.
    - h: A function that returns local relative amplitude at x.
    - inc, ecc: Anisotropy parameters. Needed if geometry=='aniso_geo_rad'.
    - n_threads: Maximum number of threads available to function.
    - symm: Flag indicating whether matrix will be symmetric (optional).
    - geometry: Must be 'aniso_geo_rad' or 'euclidean'.
    - kwds: Passed to t_gam_fun.
    
    Output value should never drop below -1. This happens when:
    -1 > -sf*c+k
    
    References:
    
    Stein, 2005. "Space-Time Covariance Functions". Journal of the American Statistical 
        Association 100(469).
    
    Pintore and Holmes, 2010, "Spatially adaptive non-stationary covariance functions
        via spatially adaptive spectra". Journal of the American Statistical Association.
        Forthcoming.
    
    """
    # Allocate 
    nx = x.shape[0]
    ny = y.shape[0]
    
    k=kwds['tlc']/kwds['sd']
    c=1./kwds['sd']-k
    sf=kwds['sf']
    tlc=kwds['tlc']
    sd=kwds['sd']
    
    if kwds.has_key('n_threads'):
        kwds.pop('n_threads')
    
    if geometry=='aniso_geo_rad':
        inc = kwds.pop('inc')
        ecc = kwds.pop('ecc')
    
    if geometry not in ['aniso_geo_rad','euclidean']:
        raise ValueError, 'Geometry %s unknown, must be aniso_geo_rad or euclidean.'%geometry
        
    
    # If parameter values are illegal, just return zeros.
    # This case will be caught by the Potential.
    if -sd >= 1./(-sf*(1-tlc)+tlc):
        return np.zeros((nx,ny),order='F')
    
    D = np.asmatrix(np.empty((nx,ny),order='F'))
    GT = np.asmatrix(np.empty((nx,ny),order='F'))
    
    # Figure out symmetry and threading
    if symm is None:
        symm = (x is y)

    n_threads = min(get_threadpool_size(), nx*ny / 10000)    
    if n_threads > 1:
        if not symm:
            bounds = np.linspace(0,ny,n_threads+1)
        else:
            bounds = np.array(np.sqrt(np.linspace(0,ny*ny,n_threads+1)),dtype=int)

    # Target function for threads
    def targ(D,GT,x,y,cmin,cmax,symm,inc=inc,ecc=ecc,amp=amp,scale=scale,diff_degree=diff_degree,h=h,geometry=geometry,kwds=kwds):
        # Spatial distance
        if geometry=='aniso_geo_rad':
            aniso_geo_rad(D, x[:,:-1], y[:,:-1], inc, ecc,cmin=cmin,cmax=cmax,symm=symm)    
        else:
            euclidean(D, x[:,:-1], y[:,:-1], cmin=cmin,cmax=cmax,symm=symm)    
        imul(D,1./scale,cmin=cmin,cmax=cmax,symm=symm)            
        # Temporal variogram
        origin_val = t_gam_fun(GT, x[:,-1], y[:,-1],cmin=cmin,cmax=cmax,symm=symm,**kwds)
        # Local properties
        ddx, ddy = diff_degree(x), diff_degree(y)
        hx, hy = h(x), h(y)
        # Covariance
        stein_spatiotemporal(D,GT,origin_val,ddx,ddy,hx,hy,cmin=cmin,cmax=cmax,symm=symm)                        
        imul(D,amp*amp,cmin=cmin,cmax=cmax,symm=symm)            
        # if symm:
        #     symmetrize(D, cmin=cmin, cmax=cmax)
    
    # Serial version
    if n_threads <= 1:
        targ(D,GT,x,y,0,-1,symm)
    
    # Parallel version
    else:   
        thread_args = [(D,GT,x,y,bounds[i],bounds[i+1],symm) for i in xrange(n_threads)]
        map_noreturn(targ, thread_args)

    if symm:
        symmetrize(D)
    
    return D
    
def nsstd(x,amp,scale,inc,ecc,diff_degree,h=default_h):
    return (h(x)*amp)**2
        
nonstationary_spatiotemporal.diag_call = nsstd