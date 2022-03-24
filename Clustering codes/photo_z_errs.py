#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 19 19:29:53 2021

@author: u92876da
"""

import config as cf

class photo_z_errs:
    
    def __init__(self, cat_1_filename, cat_2_filename = None, z_1 = "z_spec", z_2 = "z_p",
                 catastrophic_err_percentage = 15):

        # data = cf.fits.open(cf.filenames[fieldname])[1].data
        # data = cf.rewrite_matched_table_columns(data)
        # data = data[data["Best galaxies"] == True]
        # data = data[data["z_spec"] > 0]
        # data = data[data["z_p"] > 0]
        # self.photo_z = data['z_p']
        # self.photo_z_u1 = data['z_u1']
        # self.photo_z_l1 = data['z_l1']
        # self.spec_z = data['z_spec']

        if cat_2_filename == None:
            cat_2_filename = cat_1_filename
            
        catalogues = []
        for filename in [cat_1_filename, cat_2_filename]:
            data = cf.fits.open(filename)[1].data
            if filename == "/Users/user/Documents/PGR/UDS field/DR11-2arcsec-Jun-30-2019.fits":
                data = data[data["Best galaxies"] == True]
                data = data[data["z_spec"] > 0]
                data = data[data["z_p"] > 0]
            catalogues.append(data)
            
        self.z_1 = catalogues[0][z_1]
        self.z_2 = catalogues[1][z_2]
        self.z_1_name = z_1
        self.z_2_name = z_2
        self.catastrophic_err_percentage = catastrophic_err_percentage
    
    @property
    def delta_z_over_1_plus_z_spec(self):
        return (self.z_1 - self.z_2) / (1 + self.z_1)
    
    @property
    def percentage_outliers(self, round_places = 1):
        high_z_p_outliers = len(self.z_1[self.z_2 > self.z_1 + \
                        self.catastrophic_err_percentage * (1 + self.z_1) / 100])
        low_z_p_outliers = len(self.z_1[self.z_2 < self.z_1 - \
                        self.catastrophic_err_percentage * (1 + self.z_1) / 100])
        outliers = high_z_p_outliers + low_z_p_outliers
        return cf.np.round(outliers * 100 / len(self.z_1), round_places)
    
    @property
    def NMAD(self): # normalized median absolute deviation
        return 1.48 * cf.np.median(abs(self.z_2 - self.z_1) / (1 + self.z_1))
    
    def plot_photo_vs_spec_z(self, save_fig = False):
        cf.plt.scatter(self.z_1, self.z_2, s = 1, c = 'black')
        z_1 = cf.np.linspace(-0.1, 10, 2) # shorter spec_z for plotting use
        cf.plt.plot(z_1, z_1, c = 'blue') # photo_z = spec_z line
        cf.plt.plot(z_1, z_1 - self.catastrophic_err_percentage * (1 + z_1) / 100, \
                    c = 'red')
        cf.plt.plot(z_1, z_1 + self.catastrophic_err_percentage * (1 + z_1) / 100, \
                    c = 'red')
        cf.plt.text(1, 6, "η = %1.2f" % self.percentage_outliers + "%", \
                    ha = "center", va = "center", c = "red", fontsize = 10)
        
        # plot parameters
        cf.plt.xlabel(self.z_1_name, fontsize = 12)
        cf.plt.ylabel(self.z_2_name, fontsize = 12)
        cf.plt.xlim(cf.np.min(z_1), cf.np.max(self.z_1) + 0.1)
        cf.plt.ylim(cf.np.min(z_1), cf.np.max(self.z_1) + 0.5)
        cf.plt.xticks(fontsize = 10)
        cf.plt.yticks(fontsize = 10)
        cf.plt.show()
        
    def plot_delta_z_over_1_plus_z_spec(self, Nbins = 50, min_bin = -0.1, max_bin = 0.1, \
                                        save_fig = False):
        
        bins = cf.np.linspace(min_bin, max_bin, Nbins + 1)
        mid_bins = 0.5 * (bins[1:] + bins[:-1])
        def gauss_function(x, a, mu, sigma):
            return a * cf.np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))
        popt, pcov = cf.curve_fit(gauss_function, mid_bins, \
                    cf.np.histogram(self.delta_z_over_1_plus_z_spec, bins)[0])
        cf.plt.hist(self.delta_z_over_1_plus_z_spec, bins = bins, color = "gray")
        cf.plt.plot(mid_bins, gauss_function(mid_bins, *popt), c = "black")
        cf.plt.plot(mid_bins, gauss_function(mid_bins, popt[0], 0, self.NMAD), c = "red")
        cf.plt.text(-0.06, popt[0] * 0.95, "μ = %1.4f, σ = %1.4f" % (popt[1], popt[2]), \
                    ha = "center", va = "center", c = "black", fontsize = 10)
        cf.plt.text(-0.06, popt[0] * 0.85, "NMAD = %1.4f" % self.NMAD, ha = "center", \
                    va = "center", c = "red", fontsize = 10)
        
        # plot parameters
        cf.plt.xlim(-0.1, 0.1)
        cf.plt.xlabel("Δz / (1+z$_s$)", fontsize = 12)
        cf.plt.ylabel("N", fontsize = 12)
        cf.plt.xticks(fontsize = 10)
        cf.plt.yticks(fontsize = 10)
        cf.plt.show()

# -----------------------------------------------------------------------------    

def main(cat_1_filename, cat_2_filename = None, z_1 = "z_spec", z_2 = "z_p"):
    photo_z_err_obj = photo_z_errs(cat_1_filename, cat_2_filename, z_1, z_2)
    photo_z_err_obj.plot_photo_vs_spec_z()
    photo_z_err_obj.plot_delta_z_over_1_plus_z_spec()

if __name__ == "__main__":
    filename = "UDS_def_test_EM_lines_ZP_corrected_small"
    cf.os.chdir("/Users/user/Documents/PGR/LePhare/output_cat/" + filename)
    cat_1_filename = filename + ".fits"
    
    # Rachanas_UDS_DR11 = "/Users/user/Documents/PGR/UDS field/DR11-2arcsec-Jun-30-2019.fits"
    # cat_1_filename = Rachanas_UDS_DR11
    # cat_2_filename = ""
    z_1 = "zSPEC"
    z_2 = "z_BEST"
    main(cat_1_filename, z_1 = z_1, z_2 = z_2)
        