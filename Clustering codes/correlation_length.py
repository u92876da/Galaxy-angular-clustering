#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 10:44:24 2021

@author: u92876da
"""

import config as cf

class corr_length_evolution:
    
    def __init__(self, gal_field_arr, fitting_params = "A", method = "Lee06"):
        
        z_array = []
        z_err_array = []
        r0_array = []
        r0_err_array = []
        two_pt_angular_arr = [gal_field.w_theta_obj for gal_field in gal_field_arr]
        for loc_object in two_pt_angular_arr:
            z_array = cf.np.append(z_array, loc_object.z)
            z_err_array = cf.np.append(z_err_array, loc_object.z_std)
            r0 = loc_object.calculate_r0(fitting_params, method)
            r0_array = cf.np.append(r0_array, r0[0])
            r0_err_array = cf.np.append(r0_err_array, r0[1])
        
        self.z_data = z_array
        self.z_data_errs = z_err_array
        self.r0_data = r0_array
        self.r0_data_errs = r0_err_array
        
    @property
    def min_z(self):
        return cf.np.min(self.z_data)
    
    @property
    def max_z(self):
        return cf.np.max(self.z_data)
        