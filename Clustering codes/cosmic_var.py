#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 15:25:38 2021

@author: u92876da
"""

import config as cf
import field_geometry

class cosmic_var:
    
    def __init__(self, fieldname, z, gal_bias = 1):
        self.fieldname = fieldname
        self.z = z
        self.gal_bias = gal_bias
    
    @property
    def geometry(self):
        return field_geometry.field_geometry(self.fieldname)
    
    # calculated over the unmasked area, dz = 0.1
    @property # this function can be improved to remove s1 and s2 arguments
    def unmasked_dm_cosmic_var(self, s1 = 0.92858, s2 = 0.92858): #, gal_bias = 2.1):
        dm_cosmic_var_data = cf.pd.read_csv(cf.dm_cosmic_var_file)
        dm_cosmic_var_data = dm_cosmic_var_data[lambda x: x["s1"] == s1]
        dm_cosmic_var_data = dm_cosmic_var_data[lambda x: x["s2"] == s2]
        dm_cosmic_var_data = dm_cosmic_var_data[lambda x: x["sigma8"] == cf.cosmo.sigma8]
        dm_cosmic_var_data = dm_cosmic_var_data[lambda x: x["z_mid"] == self.z]
        dm_cosmic_var = dm_cosmic_var_data["cosmic_var"]
        # dm_cosmic_var = dm_cosmic_var * (0.1 / self.bin_width)
        # gal_cosmic_var = (gal_bias ** 2) * dm_cosmic_var
        return dm_cosmic_var
    
    def unmasked_gal_cosmic_var(self, z_width):
        return cf.np.array(self.unmasked_dm_cosmic_var)[0] * (self.gal_bias ** 2) * 0.1 / z_width
    
def main():
    var = cosmic_var("UDS", 1.0, 3.0)
    print(cf.np.sqrt(var.unmasked_gal_cosmic_var(0.1)))
    
if __name__ == "__main__":
    main()
        