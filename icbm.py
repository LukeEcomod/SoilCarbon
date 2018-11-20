# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 17:01:23 2018

@author: slauniai

ICBM (Introductory Carbon Balance Model) by Andren & Kätterer, 1997 Ecol. Appl.
7(4), 1226 - 1236.

A two-pool model for soil C parameterized for Swedish agricultural soils.

Samuli Launiainen 20.11.2018

TODO:
    Include environmental effects
    Make forward Euler solution for vectors.
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

EPS  = np.finfo(float).eps  # machine epsilon
NT = 273.15  # 0 degC in Kelvin
NP = 101300.0  # Pa, sea level normal pressure
R = 8.314462175  # J mol-1 K-1, universal gas constant
O2_IN_AIR = 0.209 # volume fraction of O2 in free air

class model():
    def __init__(self, para, ini, gridded=False):
        """
        ICBM -model for soil C balance consists of two pools, young (Y) and
        old (O). External inputs to the system (I) are through Y only.
        """
        self.ky = para['ky']    # decomposition rate of Y-pool
        self.ko = para['ko']    # decomposition rate of O-pool
        self.h = para['h']      # transfer coefficient from Y to O (-), f = 1 - h
        
        self.Y = ini['Y']       # initial Y-pool
        self.O = ini['O']       # initial O-pool
        #self.fT = 1.0           # decomposition modifier (-) for temperature
        #self.fW = 1.0           # decomposition modifier (-) for moisture
        
        #self.gridded=gridded    # boolean for point vs. gridded model
    
    def compute(self, t, I=0.0, fenv=1.0):
        """
        runs IBDM from initial state to time t. State and fluxes at each time
        point will be returned.
        Args:
            t - time (floar or array)
            I - inputs to Y pool (scalar or array)
            fenv - environmental modifier (-) for temperature and moisture
        """
        # massage inputs to arrays if they are not
        t = np.array(t, ndmin=1)
        I = I*np.ones(np.shape(t))
        
        # environmental modifier (-)        
        fenv = fenv*np.ones(np.shape(t))
        
        nsteps = len(t)
        # results
        C = np.zeros((2, nsteps)) # pools, nsteps
        C[0,0] = self.Y
        C[1,0] = self.O

        for m in range(1,nsteps):
            k = [self.ky * fenv[m], self.ko * fenv[m]]
            x0 = [self.Y, self.O]
            tspan = [t[m-1], t[m]] # this must be list of successive time points!
            x = odeint(time_derivative, x0, tspan, args=(I[m], k, self.h)) 
            #x, out = odeint(time_derivative, x0, tspan, args=(I[m], k, self.h), full_output=True) 
 
            # update output and model state
            C[0,m] = x[1,0]
            C[1,m] = x[1,1]
            
            self.Y = x[1,0]
            self.O = x[1,1]
        return C


def time_derivative(x, t, I, k, h):
    """
    returns time derivative of C pools in format suitable for scipy.odeint
    Args:
        x - initial pools, 
        t - time
        I - input to pool Y
        k - rate coefficients (adjusted with env. conditions)
        h - humification coefficient 
    Returs:
        derivatives
    """
    # x[0] = Y, x[1] = O
    dydt = I -k[0]*x[0]
    dodt = h*k[0]*x[0] - k[1]*x[1]
    return [dydt, dodt] 

def test():
    # test ICBM for scenarios 'fallow' and '+N+straw'
    para = {
            'fallow': 
                {'ky': 0.8,      # yr-1
                 'ko': 6.05e-3,  # yr-1
                 'h': 0.13       # -
                 },
            'N_straw':
                {'ky': 0.8,      # yr-1
                 'ko': 6.05e-3,  # yr-1
                 'h': 0.125       # -
                 }
            }
                
    ini = {'fallow': {'Y': 0.3, 'O': 3.96}, # kg C m-3
           'N_straw': {'Y': 0.3, 'O': 4.11}
           }
    litter = {'fallow': 0.0, 'N_straw': 0.285} # kg C a-1
    
    # create model instances
    run1 = model(para['fallow'], ini['fallow'])
    run2 = model(para['N_straw'], ini['N_straw'])
    
    n = 40 # yr
    t = np.linspace(0, n)
    
    # run models and return results
    C1 = run1.compute(t, I=litter['fallow'], fenv=1.32)
    C2 = run2.compute(t, I=litter['N_straw'], fenv=1.0)
    
    # plot figure()
    plt.figure()
    
    plt.subplot(121)
    plt.plot(t, C1[0], label='Y')
    plt.plot(t, C1[1], label='O')
    plt.plot(t, C1[0]+C1[1], label='Y+O')
    plt.plot(t[[0,-1]], [litter['fallow'], litter['fallow']], 'r--', label='L (kgCa-1)')
    # plt.legend()
    plt.xlabel('t (yr)')
    plt.ylabel('C pools (kg C m-3)')
    plt.title('fallow')

    plt.subplot(122)
    plt.plot(t, C2[0], label='Y')
    plt.plot(t, C2[1], label='O')
    plt.plot(t, C2[0]+C2[1], label='Y+O')
    plt.plot(t[[0,-1]], [litter['N_straw'], litter['N_straw']], 'r--', label='L (kgCa-1)')
    plt.legend()
    plt.xlabel('t (yr)')
    plt.ylabel('C pools (kgC m-3)')
    plt.title('+N+Straw')
    plt.savefig('ICBM_test.png', dpi=300)
    
