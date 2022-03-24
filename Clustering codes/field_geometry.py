#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 13:32:54 2021

@author: u92876da
"""

import config as cf

field_masked_areas = {"UDS": 0.58202} # UDS best galaxies
field_corners = {"UDS": [[34.9159924, -4.5994879], [34.9167096, -5.5917635], \
                         [33.9841102, -5.5917636], [33.9848297, -4.5994866]]}

class field_geometry:
    
    def __init__(self, fieldname):
        
        self.fieldname = fieldname
    
    @property
    def unmasked_geometry(self): # this geometry is NOT rectangular, it is a PARALLELOGRAM
        # data = cf.fits.open(cf.filenames["UDS"])[1].data
        # data = cf.rewrite_matched_table_columns(data)
        min_ra = cf.np.min([field_corners[self.fieldname][i][0] for i in range(len(field_corners[self.fieldname]))]) * cf.u.deg
        max_ra = cf.np.max([field_corners[self.fieldname][i][0] for i in range(len(field_corners[self.fieldname]))]) * cf.u.deg
        min_dec = cf.np.min([field_corners[self.fieldname][i][1] for i in range(len(field_corners[self.fieldname]))]) * cf.u.deg
        max_dec = cf.np.max([field_corners[self.fieldname][i][1] for i in range(len(field_corners[self.fieldname]))]) * cf.u.deg
        return {"min_ra": min_ra, "max_ra": max_ra, "min_dec": min_dec, "max_dec": max_dec}
    
    @property
    def unmasked_s1(self): # parallelogram base
        delta_ra = self.mask.shape[1] * self.mask_pixel_scaling
        print("delta_ra = " + str(delta_ra))
        unmasked_s1 = cf.np.mean([field_corners[self.fieldname][0][0] - field_corners[self.fieldname][3][0], \
                                  field_corners[self.fieldname][1][0] - field_corners[self.fieldname][2][0]]) * cf.u.deg
        print(field_corners[self.fieldname][0][0] - field_corners[self.fieldname][3][0])
        print(field_corners[self.fieldname][1][0] - field_corners[self.fieldname][2][0])
        #print("gals ra: " + str(self.unmasked_geometry["min_ra"]) + ", " + str(self.unmasked_geometry["max_ra"]))
        #print(self.unmasked_geometry["max_ra"] - self.unmasked_geometry["min_ra"])
        return delta_ra
    
    @property
    def unmasked_s2(self): # parallelogram height
        delta_dec = self.mask.shape[0] * self.mask_pixel_scaling
        print("delta_dec = " + str(delta_dec))
        #print("gals dec: " + str(self.unmasked_geometry["min_dec"]) + ", " + str(self.unmasked_geometry["max_dec"]))
        #print(self.unmasked_geometry["max_dec"] - self.unmasked_geometry["min_dec"])
        return delta_dec
        
    @property
    def unmasked_area(self):
        area = self.mask.shape[0] * self.mask.shape[1] * self.mask_pixel_scaling ** 2
        print("area = " + str(area))
        return self.unmasked_s1 * self.unmasked_s2
    
    @property # can now calculate this directly from the mask
    def masked_area(self):
        # should count number of 0s in masked region
        return field_masked_areas[self.fieldname] * cf.u.deg ** 2
    
    @property
    def mask(self):
        mask_dir = cf.os.environ["MASK_DIR"]
        if mask_dir == "local":
            return cf.fits.open(cf.masks[self.fieldname])[0].data
        else:
            return cf.fits.open(mask_dir + "/" + cf.UDS_mask_filename)[0].data
    
    @property
    def mask_pixel_scaling(self):
        # for UDS only
        pixel_scale = (0.1342 * cf.u.arcsec).to(cf.u.deg) #deg/pixel
        return pixel_scale
    
    @property
    def mask_bin_edges(self):
        dec_bin_edges = cf.np.linspace(self.unmasked_geometry["min_dec"].value, self.unmasked_geometry["max_dec"].value, \
                                      self.mask.shape[0])
        ra_bin_edges = cf.np.linspace(self.unmasked_geometry["max_ra"].value, self.unmasked_geometry["min_ra"].value, \
                                      self.mask.shape[1])
        return {"ra": ra_bin_edges, "dec": dec_bin_edges}

    def plot_field(self, data_gals_ra_dec, random_gals_ra_dec, data_label = "data", save_fig = False):
        extent = (self.unmasked_geometry["max_ra"].value, self.unmasked_geometry["min_ra"].value, \
                  self.unmasked_geometry["max_dec"].value, self.unmasked_geometry["min_dec"].value)
        # more accurate image without stretching (extent argument)
        cf.plt.imshow(self.mask, cmap = cf.binary, extent = extent)
        cf.plt.scatter(data_gals_ra_dec["ra"], data_gals_ra_dec["dec"], s = 2.0, c = "red", \
                       label = data_label, zorder = 2)
        cf.plt.scatter(random_gals_ra_dec["ra"], random_gals_ra_dec["dec"], s = 2.0, c = "limegreen", \
                       label = "random", zorder = 1)
        cf.plt.title(self.fieldname + " field", fontsize = 14)
        cf.plt.xlabel("RA /deg", fontsize = 12)
        cf.plt.ylim(self.unmasked_geometry["min_dec"].value, self.unmasked_geometry["max_dec"].value)
        cf.plt.ylabel("DEC /deg", fontsize = 12)
        lgd = cf.plt.legend(loc="center right", frameon=True, fontsize=10, bbox_to_anchor = (1.2, 1.0), fancybox=True, shadow=True)
        if save_fig == True:
            cf.os.chdir("/Users/user/Documents/PGR/UDS field/Plots")
            cf.plt.savefig(self.fieldname + " field " + data_label, dpi = 600, \
                           box_extra_artists = (lgd,), bbox_inches = "tight")
        cf.plt.show()
    
def test_mask():
    geom1 = field_geometry("UDS")
    #print(g1.unmasked_s1)
    #print(g1.unmasked_s2)
    #print(g1.unmasked_area)
    geom1.plot_field({"ra": [None], "dec": [None]}, {"ra": [None], "dec": [None]})    

def main():
    pass
    # g1 = field_geometry("UDS")
    # print(cf.np.round(g1.unmasked_s1, 5))
    # print(cf.np.round(g1.unmasked_s2, 5))
    
if __name__ == "__main__":
    test_mask()