#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 18:44:56 2021

@author: u92876da
"""

import config as cf

class sigma_8_gal:
    
    def __init__(self, r_0, r_0_err, delta): # delta err here too
    
        numerator = 72 * ((r_0.value / 8) ** (1 + delta))
        denominator = (3 - (1 + delta)) * (4 - (1 + delta)) * (6 - (1 + delta)) \
            * (2 ** (1 + delta))
        
        self.sigma_8_gal = cf.np.sqrt(numerator / denominator)
        self.sigma_8_gal_err = (1 + delta) * self.sigma_8_gal * r_0_err / (2 * r_0)
        
class sigma_8_gal_evolution:
    pass