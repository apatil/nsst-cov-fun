import pymc as pm
import numpy as np
import os
from copy import copy
# from scipy import interpolate as interp
from fst_cov_fun import my_gt_fun
from threading import Thread, Lock
import threading
from pymc.gp.cov_funs import imul, symmetrize, stein_spatiotemporal
from pymc.gp.cov_funs import aniso_geo_rad
#import MAPdata
from IPython import Debugger
from IPython.Debugger import Pdb

__all__ = ['my_st', 'my_gt_fun']

t_gam_fun = my_gt_fun

def my_st(x,y,amp,scale,inc,ecc,n_threads=1,symm=None,**kwds):
    """
    Spatiotemporal covariance function. Converts x and y
    to a matrix of covariances. x and y are assumed to have
    columns (long,lat,t). Parameters are:
    - t_gam_fun: A function returning a matrix of variogram values.
      Inputs will be the 't' columns of x and y, as well as kwds.
    - amp: The MS amplitude of realizations.
    - scale: Scales distance.
    - inc, ecc: Anisotropy parameters.
    - n_threads: Maximum number of threads available to function.
    - symm: Flag indicating whether matrix will be symmetric (optional).
    - kwds: Passed to t_gam_fun.
    
    Output value should never drop below -1. This happens when:
    -1 > -sf*c+k
    
    """
    # Allocate 
    nx = x.shape[0]
    ny = y.shape[0]
    
    k=kwds['tlc']/kwds['sd']
    c=1./kwds['sd']-k
    sf=kwds['sf']
    tlc=kwds['tlc']
    sd=kwds['sd']
    
    # If parameter values are illegal, just return zeros.
    # This case will be caught by the Potential.
    if -sd >= 1./(-sf*(1-tlc)+tlc):
        return np.zeros((nx,ny))
    
    D = np.asmatrix(np.empty((nx,ny),order='F'))
    GT = np.asmatrix(np.empty((nx,ny),order='F'))
    
    # Figure out symmetry and threading
    if symm is None:
        symm = (x is y)

    n_threads = min(n_threads, nx*ny / 10000)        
    
    if n_threads > 1:
        if not symm:
            bounds = np.linspace(0,ny,n_threads+1)
        else:
            bounds = np.array(np.sqrt(np.linspace(0,ny*ny,n_threads+1)),dtype=int)
    
    # Serial version
    if n_threads <= 1:
        # Spatial distance
        aniso_geo_rad(D, x[:,:-1], y[:,:-1], inc, ecc, symm=symm)
        imul(D,1./scale,symm=symm)
        # Temporal variogram
        origin_val = t_gam_fun(GT, x[:,-1], y[:,-1],symm=symm,**kwds)
        # Covariance
        stein_spatiotemporal(D,GT,origin_val,symm=symm)
        imul(D,amp*amp,symm=symm)
        if symm:
            symmetrize(D)
                    
    
    # Parallel version
    else:
        
        # Target function for threads
        def targ(cmin,cmax):
            # Spatial distance
            aniso_geo_rad(D, x[:,:-1], y[:,:-1], inc, ecc,cmin=cmin,cmax=cmax,symm=symm)    
            imul(D,1./scale,cmin=cmin,cmax=cmax,symm=symm)            
            # Temporal variogram
            origin_val = t_gam_fun(GT, x[:,-1], y[:,-1],cmin=cmin,cmax=cmax,symm=symm,**kwds)
            # Covariance
            stein_spatiotemporal(D,GT,origin_val,cmin=cmin,cmax=cmax,symm=symm)                        
            imul(D,amp*amp,cmin=cmin,cmax=cmax,symm=symm)            
            if symm:
                symmetrize(D, cmin=cmin, cmax=cmax)        
        
        # Dispatch threads        
        threads = []
        for i in xrange(n_threads):
            new_thread= Thread(target=targ, args=(bounds[i],bounds[i+1]))
            threads.append(new_thread)
            try:
                new_thread.start()
            except:
                print i,threading.activeCount()
                Pdb(color_scheme='Linux').set_trace()
        [thread.join() for thread in threads]       
    
    # return D
    return D

    # def my_GT_fun(tx,ty,scal_t,t_lim_corr,sin_frac,space_diff):
    #     """
    #     Converts two vectors of times, tx and ty, into a 
    #     matrix whose i,j'th entry is gamma(abs(t[i]-t[j])),
    #     gamma being Stein's 'valid variogram'. Parameters of
    #     this variogram are:
    #     - scal_t: Scales time.
    #     - t_lim_corr: The limiting correlation between two points
    #       as the distance between them grows. Note that this will be
    #       multiplied by the overall 'amp' parameter.
    #     - space_diff: The desired degree of differentiability of
    #       the spatial margin.
    #     - sin_frac: The fraction of the partial sill taken up by the first harmonic.
    #     """
    # 
    #     k = t_lim_corr/space_diff
    #     c = 1./space_diff-k
    #     dt = np.asarray(abs(np.asmatrix(tx).T-ty))
    #     GT = 1./((np.exp(-dt/scal_t)*(1.-sin_frac) + sin_frac*np.cos(2.*np.pi*dt))*c+k)
    #     return GT, 1./(k+c)


# def my_cov_fun(x, y, amp, scale, scale_t, inc, ecc, symm=False):
#     """
#     Exponential model, common metric.
#     """
#     tx = np.asmatrix(x[:,2])
#     ty = np.asmatrix(y[:,2])
#     
#     Dx = pm.gp.cov_funs.aniso_geo_rad(x[:,:2],y[:,:2],inc,ecc,symm)
#     Dt = tx.T - ty
#     
#     C = np.asmatrix(np.exp(-np.sqrt((np.asarray(Dx)/scale)**2 + (np.asarray(Dt)/scale_t)**2))) * amp**2
#     return C
