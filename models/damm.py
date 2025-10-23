# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 20:22:18 2018

@author: slauniai

DAMM (Dual Arhenius Michelis-Menten model), Davidson et al. 2012 Global Change
Biol. 18, 371-384.

Formulation as in original paper, parameters in test-function from the field-
scale parameterization.
"""
import numpy as np
import matplotlib.pyplot as plt

EPS  = np.finfo(float).eps  # machine epsilon
NT = 273.15  # 0 degC in Kelvin
NP = 101300.0  # Pa, sea level normal pressure
R = 8.314462175  # J mol-1 K-1, universal gas constant
O2_IN_AIR = 0.209 # volume fraction of O2 in free air


class Damm():
    def __init__(self, para, St, poros):
        """
        sets up DAMM model from parameter dict
        para: with keys
            alpha - pre-exponential term
            Ea - activation energy
            kMs - michaelis constant for S
            kMo2 - michaelis constant for O2
            p - fraction of soluble to total substrate
            Dliq - liquid phase diffusion coefficent (-)
            Dair - air-phase diffusion coefficient for O2
        St - total carbon substrate concentration in soil
        poros - porosity (m3m-3)
        """
        self.alpha = para['alpha']
        self.Ea = para['Ea']
        self.kMs = para['kMs']
        self.kMo2 = para['kMo2']
        self.p = para['p']
        self.Dliq = para['Dliq']
        self.Dgas = para['Dgas']
        #self.para = para
        
        self.St = St # total substrate concentration in soil
        self.poros = poros
        
        
    def reaction_velocity(self, T, W):
        """
        reaction velocity (eq. 1-6)
        Args:
            para - parameter dict with keys:
                    alpha - pre-exponential term (mg C cm-3 soil)
                    Ea - activation energy (kJ mol-1)
                    kMs - michaelis constant for S (g C cm-3 soil)
                    kMo2 - michaelis constant for O2 (cm3 O2 cm-3 air)
            T - temperature (K)
            W - liquid water content (m3m-3)
        Returns:
            v - reaction velocity (units s-1)
            c - dict of component terms (Vmax, fs, fo2)
        """
        # maximum reaction velocity (1e-3 converts R to kJ mol-1 K-1)
        Vmax = self.alpha * np.exp( -self.Ea / (1e-3*R*T) )
        
        # substrate and oxygen concentrations
        S = self.p * self.St * self.Dliq * W**3.0 
      
        a = np.maximum(0.0, self.poros - W)  # air-filled porosity
        O2 = self.Dgas * O2_IN_AIR * a**(4./3.)
            
        # michaelis-menten equations
        fs = S / (self.kMs + S)
        fo2 = O2 / (self.kMo2 + O2)
        
        # reaction velocity
        v = Vmax * fs * fo2
        
        return v, {'Vmax': Vmax, 'fs': fs, 'fo2': fo2}


def test_model():
    # set up DAMM and reproduce fig. 5 in Davidson et al. 2012.
    
    para = {'alpha': 5.38e10,  # mg C cm-3 soil h-1
            'Ea': 72.26,     # kJ mol-1
            'kMs': 9.95e-7, # g C cm-3 soil
            'kMo2': 0.121,  # cm3 O2 cm-3 air
            'p': 4.14e-4,   # -
            'Dliq': 3.17,   # -
            'Dgas': 1.67    # -
            }
    St = 0.048 # gC cm-3 soil
    poros = 0.68 # cm3 cm-3
    
    model = Damm(para, St, poros)
    
    T = np.arange(5.0, 25.0, 0.1)  + NT
    W = np.arange(0.1, 0.68, 0.01)
    
    T, W = np.meshgrid(T, W)
    
    R, c = model.reaction_velocity(T, W)
    
    f = 100. * 100. *10.  # unit conversion: mg C cm-3 to mg C m-2 from 10cm layer
    
    
    '''
    ======================
    3D surface (color map)
    ======================
    
    Demonstrates plotting a 3D surface colored with the coolwarm color map.
    The surface is made opaque by using antialiased=False.
    
    Also demonstrates using the LinearLocator and custom formatting for the
    z axis tick labels.
    '''
    
    from mpl_toolkits.mplot3d import Axes3D
    # import matplotlib.pyplot as plt
    from matplotlib import cm
    # from matplotlib.ticker import LinearLocator, FormatStrFormatter
    
    ele = 35.0
    azm = 145.0
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.view_init(elev=ele, azim=azm)
    
    x = T - NT
    y = W
    z = f*R
    
    # Plot the surface.
    surf = ax.plot_surface(x, y, z, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False, alpha=0.6)
    ax.set_xlabel('T (deg C)'); ax.set_ylabel('W (-)'); ax.set_zlabel(r'R')
    ax.invert_yaxis()
    # Add a color bar which maps values to colors.
    # fig.colorbar(surf, shrink=0.5, aspect=5)
    
    ax.set_title('Damm sensitivity')
    plt.show()
    plt.savefig('damm_sensitivities.png')
    
    
#def carbon_concentration(p, Dliq, St, W):
#    """
#    soluble carbon substrate concentration at the reaction site
#    Args:
#        p - soluble fraction of total substrate (-)
#        Dliq - non-dimensional diffusion (liquid phase) parameter (-)
#        St - total substrate concentration in soil (units)
#        W - vol. liquid water content (m3m-3)
#    Returns:
#        s - concentration at reaction site (soluble) (units)
#    """
#    s = p * Dliq * W**3.0
#    
#    return s
#
#def oxygen_concentration(Dgas, a):
#    """
#    Oxygen concentration at the reaction site
#    Args:
#        Dgas - non-dimensional (gas-phase) diffusion parameter
#        a - air-filled porosity (m3m-3)
#    Returns:
#        oxygen concentration at the reaction site
#    """
#    return Dgas * O2_IN_AIR * a**(4./3.)

