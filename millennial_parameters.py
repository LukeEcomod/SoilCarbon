# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 09:25:35 2018

@author: slauniai

Parameter dictionary for Millennial -model
"""

"""
Parameter values from Abramoff et al. (2017) Table 3.
Environmental modifiers are given separately.
"""
param = {
        'dt': 1.0,      # timestep (d)
        'pa': 0.66,     # fraction of A breakdown allocated to P
        'V_pl': 10.0,   # max rate of P decomosition (g C m-2 d-1)'
        'K_pl': 150.0,  # half-sat. const. P --> LMWC (g C m-2)
        'K_pe': 12.0,   # half-sat. const. for mibrob. control of P mineralization (g C m-2)
        'V_pa': 2e-3,   # max rate of P --> A formation (g C m-2 d-1)
        'K_pa': 50.0,   # half-sat. const. P --> A formation (g C m-2)   
        'A_max': 500.0, # max. capacity of A carbon (g C m-2)
        'k_b': 2e-4,    # rate of A breakdown (-)
        'k_l': 1.5e-3,  # leaching rate of L (g C m-2 d-1) !!! NOTE: relate to water fow in MLM !!!
        'K_lm': 0.25,   # binding affinity for L sorption (g C m-2): NOTE: wrong units!
        'pH': 7.0,      # (-)
        'c': [0.297, 3.355, 0.5], # coeffs of max. sorption capacity (-)
        'Qmax': None,   # max sorption capacity (g C m-2): COMPUTED IN initialization
        'k_s': 2e-1,    # rate of L sorption (from FORTRAN-code of Abramoff et al. 2017)
        'V_lm': 0.35,   # pot. L turnover rate (g C m-2 d-1)
        'K_lb': 7.2,    # halt-sat. const. for microbial activity (g C m-2)
        'CUE': None,     # microbial carbon use efficiency (-) 
        'CUEp': [0.6, 15.0, -1.2e-2], # CUE params [CUEref, Tref, sensitivity]               
        'V_ma': 7e-2,   # max rate of M --> A (g C m-2 d-1)
        'K_ma': 200.0,  # half-sat. const M --> A (g Cm-2)
        'k_mm': 2.5e-2, # microbial biomass adsorption rate (-)
        'k_m': 3.6e-2,  # microbial turn-over rate (d-1)
        }

