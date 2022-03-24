#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 10:52:41 2021

@author: u92876da
"""

import config as cf
from halo_model import halo_model

class HOD_abstract(cf.ABC): # abstract HOD class

    constrained_params : dict = {}

    def define_model(self, model):
        model_type = {"pure power law": self.pure_power_law}
        return model_type[model]
    
    @cf.abstractmethod
    def pure_power_law(self):
        pass
    
    def param_constraints(self, number_density_obj_data, w_theta_obj_data):
        pass
    
    def calc_w_theta(self): # in two_point_angular_dm derived class
        
        pass
    
    # put these two functions in the halo model class
    def calc_P_k_1h(self):
        pass
    
    def calc_P_k_2h(self): # y ~ 1 on large scales much larger than halo virial radius
        # assuming linear bias model of Mo, White 1996
        
        
        pass

    
    def calc_mean_gal_pairs(self, M, params):
        
        cen_params, sat_params = params
        N_g = self.full_HOD_model(M, cen_params, sat_params)
        
        if N_g > 1:
            return N_g ** 2
        elif N_g > 0.25:
            return (N_g ** 2) * cf.np.log10(4 * N_g) / cf.np.log10(4)
        else:
            return 0

# -----------------------------------------------------------------------------
    
class HOD_centrals(HOD_abstract):
    
    def pure_power_law(self, M, params):
        return 0
   
# -----------------------------------------------------------------------------

class HOD_satellites(HOD_abstract):

    def pure_power_law(self, M, params):
        if M >= params["M_min"]:
            return (M / params["M_1"]) ** params["alpha"]
        else:
            return 0

# -----------------------------------------------------------------------------

class HOD(HOD_centrals, HOD_satellites):
    
    def __init__(self, cen_model, sat_model):
        self.central_model = HOD_centrals().define_model(cen_model)
        self.satellite_model = HOD_satellites().define_model(sat_model)

    def full_HOD_model(self, M, params):
        cen_params, sat_params = params
        return self.central_model(M, cen_params) + self.satellite_model(M, sat_params)

# -----------------------------------------------------------------------------    

def main():
    r = HOD("pure power law", "pure power law")
    print(r.central_model(1e10, [0]))
    print(r.satellite_model(1e10, {"M_min": 1e10, "M_1": 1e2, "alpha": 1}))
    r.calc_mean_gal_pairs(1e10, [{}, {"M_min": 1e10, "M_1": 1e2, "alpha": 1}])
 
if __name__ == "__main__":   
    main()

    
