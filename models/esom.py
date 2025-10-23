# -*- coding: utf-8 -*-
"""
ESOM - Extendable soil organic matter model

Created on Tue Feb 15 12:01:18 2022

@author: Samuli Launiainen & Annamari Laurén
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

EPS  = np.finfo(float).eps  # machine epsilon
NT = 273.15  # 0 degC in Kelvin

def temperature_functions():
    t2 = interp1d([-40.,-5., -1, 25., 35., 60.],         # effect of temperature on the decomposition rate
                [0.,   0.,  0.2, 1.53, 1.53, 0.])
    t3 = interp1d([-40., -3., 0., 7., 60.],
                   [0., 0., 1.3, 1.3, 0.])
    t4 = interp1d([-40., -5., 1., 20., 40., 80.],
                   [0., 0., 0.2, 1., 1., 0.])
    t5 = interp1d([-40., -5., 1., 13., 25., 50.],
                   [0., 0., 0.2, 1., 1., 0.])
    t6 = interp1d([-40., -5., 1., 27.5, 35., 60.],
                  [0., 0., 0.2, 1.95, 1.95, 0.])
    #t6 = interp1d([-40., -30., -20. ,-10.,   0.,  10.,  20.,  30.,  40.,  50.],             #Q10 = 2
    #               [ 0.03125,0.0625, 0.125, 0.25, 0.5, 1., 2., 4., 4., 4.]) 
    
    #t7 = interp1d([-40., -5., 1., 27.5, 35., 60.],
    #               [0., 0., 0.2, 1.95, 1.95, 0.])
    t7 = interp1d([-40., -30., -20. ,-10.,   0.,  10.,  20.,  30.,  40.,  50.],             #Q10 = 2
                   [ 0.03125,0.0625, 0.125, 0.25, 0.5, 1., 2., 4., 4., 4.]) 
    return (t2, t3, t4, t5, t6, t7)

def pH_from_sfc(sfc):
    # Mese data Raija Laiho
    # input site fertility class, output pH in water solution
    dph = {1:3.5, 2:3.4, 3:3.3, 4:3.2, 5:3.1, 6:3.0}
    pH = np.zeros(np.shape(sfc))
    for key in dph.keys():
        ix = np.where(sfc==key)
        pH[ix] = dph[key]
    return pH

def moisture_functions():
    # Description of Romul model Table 2
    phi1236 = interp1d([0.02, 0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4, 0.417, 1.333, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8,4], 
                      [0.0, 0.004, 0.026, 0.074, 0.154, 0.271, 0.432, 0.64,  0.899, 1.0, 1.0, 0.844, 0.508, 0.305, 0.184, 0.111, 0.067, 0.04, 0.024, 0])
    phi4 = interp1d([0.0, 0.133, 1.333, 2.333,4.0],
                    [0.0, 1.0, 1.0 ,0.0, 0.0])
    phi5 = interp1d([0.0, 0.067, 0.5, 2.333, 4.0, 10.0],
                     [0.0, 0.0, 1.0, 1.0 ,0.0, 0.0])
    return phi1236, phi4, phi5

    
def lignin_corrections(nitrogen = 0.7, lignin = 25.0):
    #N contents, branches Skonieczna,et al 2014 unit %
    #Lignin contents Kilpeläinen et al. 2003 unit %
    adjust = 2. #4.
    mu_k1 = 0.092 * (lignin / nitrogen)**-0.7396 * adjust
    mu_k2 = 0.0027 * (lignin / nitrogen)**-0.3917 * adjust
    mu_k3 = 0.062 * (lignin / nitrogen)**-0.3972 * adjust
    return  mu_k1, mu_k2, mu_k3


def get_rates(ash, N, pH, tair, tp_top, tp_middle, tp_bottom, t2, t3, t4, t5, t6, t7, wn, phi1236, phi4, phi5, peat_w1,peat_w2, peat_w3, H_w):
    """
    Args:
        ash - ash content in gravimetric %
        N - N content in gravimetric %
        Ta  - air temperature in deg C
        t2...t7 temperature functions
        wn normalaized water content w/wfc
    phi1236, phi4, phi5 moisture functions
    """
    nu = np.clip(0, 0.701*pH -1.6018 - 0.038*pH**2, 1)   # ph Romul documentation Table 1

    k1= (0.002 + 0.00009*ash + 0.003*N)*min(0.1754*np.exp(0.0871*tair), 1.)*phi1236(wn)*nu # adjusted decomposition rates
    k2= np.clip((0.00114 -0.00028*N)*t2(tair)*phi1236(wn)*nu, 0., 1.)     #
    k3= np.clip((0.04 - 0.003*N)*t3(tair)*phi1236(wn), 0., 1.) 
    k4= 0.005*N*t4(tair)*phi4(wn) 
    k5= 0.007*t5(tair)*phi5(wn)
    k6= 0.0006*t6(tp_top)*phi1236(wn) #* 0.5
    #k6= 0.0006*t6(tp_top)*H_w 
    
    #####  THESE can be modified by you
    k7c = 2.0
    k8c = 1.0
    k9c = 1.0
    ####   UNTIL HERE

    k7= 0.0001*t7(tp_top)*peat_w1 * k7c       #Change this                                # Lappalainen et al 2018, gamma/VfAir slightly decomposed peat
    k8 = 0.0001*t7(tp_middle)*peat_w2 * k8c                                               # Lappalainen et al. 2018 gamma/VfAir highly decomposed
    k9 = 0.0001*t7(tp_bottom)*peat_w3 * k9c

    
    return (k1, k2, k3, k4, k5, k6, k7, k8, k9)
