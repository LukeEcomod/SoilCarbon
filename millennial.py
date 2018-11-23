# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 21:55:17 2018

@author: slauniai

Millennial - model (Abramoff et al. 2017 Biogeochemistry)

see also: https://github.com/email-clm/Millennial/blob/master/main.F90
"""

import numpy as np
# from scipy.integrate import odeint
import matplotlib.pyplot as plt
from collections import namedtuple

from millennial_parameters import param   # get model default parameters

EPS  = np.finfo(float).eps  # machine epsilon
NT = 273.15  # 0 degC in Kelvin


# define namedtuple constructor for inputting model parameters to odeint
millennial_param = namedtuple('millennial_param', ' '.join(sorted(param.keys())))


class Millennial():
    def __init__(self, p, soilp, C0, results=False):   
        """
        Implementation of single-layer Millennial model.
        Args:
            p - Millennial parameters (dict)
            soilp - soil type related parameters (dict)
            C0 - initial pools (g C m-2)
            fT, fW - temperature and moisture functions to use
            results - boolean
            
        Note:
            x[0] = POM (particulate organic matter)
            x[1] = LMWC (doc)
            x[2] = B (microbial biomass)
            x[3] = A (aggregate C)
            x[4] = MAOM (mineral associated C)
        Note:
            Qmax and Langmuir adsorption equation is likely wrong!
        """
        
        Qmax = soilp['bd']*10**(p['c'][0] * np.log(soilp['clay'] + p['c'][1]))
        print(Qmax)
        p['Qmax'] = 4550.0; # Qmax # g C m-2
        self.dt = p['dt']  # d-1
        self.para = p
        self.soilpara = soilp
        
        self.temperature_response = fT_century
        self.moisture_response = fW_century
        
        # carbon pools
        self.Cpools = C0

        # results dict for testing outputs
        if results:
            self.results = {'Cpools': None, 'flx': None}

    def decompose(self, T, W, F_in, F_adv=0.0):
        """
        Computes decomposition and pool transitions during timestep dt using
        Eulerian method. Updates state variable self.Cpools and returns fluxes.
        Args:
            T - temperature (degC)
            W - vol. moisture (m3 m-3)
            F_in - array of litter input (g C m-2)
            F_adv - net outflow of LMWC (advection kg C m-2 timestep-1)
    
        Returns:
            flx (dict), all in (g C m-2 timestep-1)
                Fpl - decomposition of POM
                Fpa - aggregate C formation from POM
                Fa  - aggregate C breakdown
                Flb - microbial uptake of LMWC
                Fbm - adsorption of microbial necromass to MAOM
                Fl -  leaching loss of LMWC
                Fma - MAOM to aggregate
                Flm - adsorption of LMWC to minerals
                Fmr - microbial maintenance respiration
                Fgr - microbial growth respiration
        Updates state variable self.Cpools
        """
        
        x = self.Cpools.copy()
        dt = self.dt
        p = self.para
        p['CUE'] = p['CUEp'][0] - p['CUEp'][2] * (T - p['CUEp'][1])
        #print(p['CUE'])
        
        # parameters into named tuple (immutable)
        p = millennial_param(**p)

        # environmental modifiers
        fT = self.temperature_response(T)
        fW = self.moisture_response(W / self.soilpara['fc'])
        
        # fT = 1.0; fW=1.0
        
        # compute fluxes
        F = fluxes(x, p, env_f=fT*fW) 
        x0 = x.copy()
        
        """ integrate in time and update new pools """
        
        x[0] += F_in[0] + dt * (p.pa*F['Fa'] - F['Fpa'] - F['Fpl'])
        x[1] += F_in[1] + dt * (F['Fpl'] - F['Flb'] - F['Flm'] - F['Fl'])
        x[2] += dt * (F['Flb'] - F['Fbm'] - F['Fmr'])
        x[3] += dt * (F['Fpa'] + F['Fma'] - F['Fa'])
        x[4] += dt * (F['Flm'] + F['Fbm'] + (1.0 - p.pa)*F['Fa'] - F['Fma'])
        
        self.Cpools = x.copy()
        
        # delta C = F_in - Fmr since growth respiration Fgr is bypass
        mbe = sum(self.Cpools) - sum(x0)  - sum(F_in) + dt*F['Fmr']
        
        return F, mbe

def fluxes(x, p, dt=1.0, env_f=1.0):
    """
    Computes fluxes between C pools.
    Args:
        x - C pools (list or array)
            x[0] = POM (particulate organic matter)
            x[1] = LMWC (doc)
            x[2] = B (microbial biomass)
            x[3] = A (aggregate C)
            x[4] = MAOM (mineral associated C)
        p - parameters (namedtuple millennial_param)
        dt - timestep (d-1)
        env_f - environmental effects modifier (-)
    Returns:
        flx (dict), all in (g C m-2 d-1)
            Fpl - decomposition of POM
            Fpa - aggregate C formation from POM
            Fa  - aggregate C breakdown
            Flb - microbial uptake of LMWC
            Fbm - adsorption of microbial necromass to MAOM
            Fl -  leaching loss of LMWC
            Fma - MAOM to aggregate
            Flm - adsorption of LMWC to minerals
            Fmr - microbial maintenance respiration
            Fgr - microbial growth respiration
    """
    
    """ fluxes (g C m-2 d-1). Constraints as in Fortran-code """
    Fpl = env_f * p.V_pl * (x[0] / (p.K_pl + x[0])) * (x[2] / (p.K_pe + x[2])) 
    Fpl = np.minimum(dt*Fpl, 0.9*x[0]) / dt # constrains flux to be <= 0.9*poolsize/dt
    
    Fpa = env_f * p.V_pa * (x[0] / (p.K_pa + x[0])) * (1.0 - x[3] / p.A_max)
    Fpa = np.minimum(dt*Fpa, 0.9*x[0]) / dt
    
    Fa = env_f * p.k_b * x[3]

    
    # Note! additional M-M -term from fortran code aa:    
    #aa = (x[1] / (20.0 + x[1]))
    aa = 1.0
    Flb = env_f * p.V_lm * x[1] * aa * p.CUE
    Flb = np.minimum(dt*Flb, 0.9*x[1]) / dt
    
    Fgr = env_f * Flb *(1.0 - p.CUE) / p.CUE
        
    Fmr = env_f * p.k_m * x[2]

    Fbm = env_f * p.k_mm * x[2]
    Fbm = np.minimum(dt*Fbm, 0.9*x[2]) / dt
    
    """ leaching: parameterize p.k_l through net water flow? """
    Fl = env_f * p.k_l * x[1]
    # Fl = 0.0
    
    Fma = env_f * p.V_ma * (x[4] / (p.K_ma + x[4])) * (1.0 - x[3] / p.A_max)
    Fma = np.minimum(dt*Fma, 0.9*x[4]) / dt
    
    """ check Langmuir isotherm parameters and eq. 9-11 and Mayes et al. 2012
    Adsorption of L to M depends on %clay and pH:
    p.Qmax = BD*10**(c1*log(%clay + c2)) but given here as pre-computed input
    binding affinity:
    p.K_lm = 10**(-186*pH - 0.216). for pH=7 this gives 0.03
    but Table 3 gives p.K_lm = 0.25 kg C m-2. Units should be m2 kgC-1.
    
    Also rate-coefficient [d-1] is missing in the eq.9 t convert the units to units
    of flux (kg C m-2 d-1).
    
    See Millenium Fortran code L ~701-->
    """
    
    # added rate-constant p.k_s (d-1) to parameters
    Flm = env_f * p.k_s * x[1] * ( (p.K_lm * p.Qmax * x[1]) / 
                                   (1.0 + p.K_lm * x[1]) - x[4]) / p.Qmax

    #Flm = np.minimum(dt*Flm, 0.9*x[1]) / dt
    
    flx = {'Fpl': Fpl, 'Fpa': Fpa, 'Fa': Fa, 'Flb': Flb, 'Fbm': Fbm, 'Fl': Fl,
           'Fma': Fma, 'Flm': Flm, 'Fmr': Fmr, 'Fgr': Fgr}
    
    return flx

def fT_century(T):
    """ 
    century temperature function (-)
    """
    c = [15.4, 11.75, 29.7, 0.031]
    
    nom = c[1] + (c[2] / np.pi) * np.arctan(np.pi * c[3] * (T - c[0]))
    den = c[1] + (c[2] / np.pi)* np.arctan(np.pi * c[3] * (30.0 - c[0]))

    f = nom / den
        
    return f

def fW_century(rwc):	
    """
    century soil moisture function (-)
    Arg:
        rwc - relative moisture, rwc = w / fc
    """
    f = 1.0 / (1.0 + 30.0 * np.exp(-9.0 * rwc))
    
    return f


""" *** testing scripts *** """

def test_millennial():
    """
    tests millennial using forcing data from Abramoff et al. 2017.
    Global average soil temperature, vol moisture and example litter input.
    Loop data over M years
    """
    # import parameters
    from millennial_parameters import param
    
    # load forcing file
    forc = np.loadtxt(r'c:\repositories\soilcarbon\data\millennial_globalaverage_data.txt', skiprows=1)
    T = forc[:,0] # degC
    W = forc[:,1] # m3m-3
    F_litter = forc[:,2] # g C d-1  

    M = 200 # yrs
    N = 365 * M # days
    
    soilp = {'clay': 40.0, 'bd': 1350.0, 'poros': 0.5, 'fc': 0.3}
    C0 = 1.0 * np.ones(5) # g C m-2 initial pools
    
    # create instance
    model = Millennial(param, soilp, C0, results=True)
    
    # create holders for daily data
    res = np.zeros((5, N))*np.NaN
    F = {'Fpl': np.zeros(N), 'Fpa': np.zeros(N), 'Fa': np.zeros(N), 'Flb': np.zeros(N),
           'Fbm': np.zeros(N), 'Fl': np.zeros(N), 'Fma': np.zeros(N), 'Flm': np.zeros(N),
           'Fmr': np.zeros(N), 'Fgr': np.zeros(N)}
    
    mbe = np.zeros(N)* np.NaN
    
    j = 0
    for yr in range(M):
        print('Run year: ', yr)
        for k in range(365):   
            F_in = [0.66*F_litter[k], 0.34*F_litter[k]]
            flx, err = model.decompose(T[k], W[k], F_in, F_adv=0.0)
            res[:,j] = model.Cpools
            for m in F.keys():
                F[m][j] = flx[m] 
            #F['Fgr'][k] = flx['Fgr']
            #F['Fmr'][k] = flx['Fmr']
            mbe[j] = err
            j +=1
    
    # plot figs
    tt = np.arange(N) / 365.0
    poolname = ['POM', 'LMWC', 'MIC', 'AGG', 'MAOM']
    
    plt.figure(200)
    for n in range(5):
        plt.plot(tt, res[n,:]/1000, label=poolname[n])
    plt.legend()
    plt.ylabel('kg C m-2'); plt.xlabel('yr')
    plt.savefig('millennial_pools.png')

    plt.figure(300)
    for m in ['Fbm','Fma','Flm','Fa']:
        plt.plot(tt, F[m], label=m)
        plt.legend()
    for m in ['Fpl','Fpa','Flb']:
        plt.plot(tt, F[m], ':', label=m)
        plt.legend()    
    for m in ['Fgr','Fmr']:
        plt.plot(tt, F[m], '-', label=m)
        plt.legend()      
    plt.ylabel('Flux g C d-1')
    
    return model, res, F
    
#def dCdt(x, t, F_in, F_adv, env_f, para):
#    """
#    Derivatives of C pools as required by scipy.odeint
#    Args:
#        x - array of initial C pools:
#            x[0] = POM (particulate organic matter)
#            x[1] = LMWC (doc)
#            x[2] = B (microbial biomass)
#            x[3] = A (aggregate C)
#            x[4] = MAOM (mineral associated C)
#        t - array of timepoints [timeunit]
#        F_in - array of litter input
#        F_adv - net outflow of LMWC (advection)
#        env_ef - environmental modifier
#        para - namedtuple of parameters
#    
#    Returns:
#        derivatives in array
#    """
#    
#    # parameters into named tuple
#    p = millennial_param(*para)
#    #print(p)
#    
#    """ fluxes (g C m-2 d-1); later multiplied with fS """
#    Fpl = p.V_pl * (x[0] / (p.K_pl + x[0])) * (x[2] / (p.K_pl + x[2]))  # POM --> L
#    
#    Fpa = p.V_pa * (x[0] / (p.K_pa + x[0])) * (1.0 - x[3] / p.A_max)  # POM --> A
#    
#    Fa = p.k_b * x[3]  # A breakdown
#    
#    Flb = p.V_lm * x[1] * (x[0] / (p.K_lb + x[0])) * p.CUE  # L --> B
#    # f_gr = Flb *(1.0 - CUE) / CUE  # growth respiration, output from system 
#    
#    Fbm = p.k_mm * x[2]  # B --> M, necromass absorption
#    
#    Fmr = p.k_m * x[2]  # microbial maintenance respration
#
#    """ parameterize p.k_l through net water flow? """
#    Fl = p.k_l * x[1]  # leaching of L
#    
#    Fma = p.V_ma * (x[4] / (p.K_ma + x[4])) * (1.0 - x[3] / p.A_max)  # MAOM --> A
#    
#    """ check Langmuir isotherm parameters and eq. 9-11 and Mayes et al. 2012
#    Adsorption of L to M depends on %clay and pH:
#    p.Qmax = BD*10**(c1*log(%clay + c2)) but given here as pre-computed input
#    binding affinity:
#    p.K_lm = 10**(-186*pH - 0.216). for pH=7 this gives 0.03
#    but Table 3 gives p.K_lm = 0.25 kg C m-2. Units should be m2 kgC-1.
#    Also rate-coefficient [d-1] is missing in the eq.9 t convert the units to units
#    of flux (kg C m-2 d-1).
#    
#    See Millenium Fortran code L ~700-->
#    """
#    
#    # added rate-constant p.k_s (d-1)
#    Flm = p.k_s * x[1] * ( (p.K_lm * p.Qmax * x[1]) / 
#                         (1.0 + p.K_lm * x[1]) - x[4]) / p.Qmax  # L --> M 
#    
#    """ derivatives y = dx / dt """
#    y = [None]*5
#    
#    y[0] = F_in[0] + env_f * (p.pa*Fa - Fpa - Fpl)    # dP/dt
#    y[1] = F_in[1] + env_f * (Fpl - Flb - Flm - Fl)   # dL/dt
#    y[2] = env_f * (Flb - Fbm - Fmr)                  # dB/dt
#    y[3] = env_f *(Fpa + Fma - Fa)                    # dA/dt
#    y[4] = env_f *(Flm + Fbm + (1.0 - p.pa)*Fa - Fma) # dM/dt
#    
#    return y