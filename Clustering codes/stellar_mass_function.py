#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 15:46:25 2022

@author: u92876da
"""

import config as cf

class stellar_mass_function:
    
    def __init__(self, photo_z, smf_bins, smf_values):
        self.photo_z = photo_z
        self.smf_bins = smf_bins
        self.smf_values = smf_values
        
    @classmethod
    def calculate_smf(cls, stellar_masses, photo_z, field_geometry, Nbins = 10):
        #stellar_mass_bin_edges = cf.np.linspace(min_bin, max_bin, num = Nbins)
        def calculate_comoving_volume(z1, z2, solid_angle): # z2 > z1; input solid angle in deg**2
            def integrand(u):
                integrand = cf.astropy_cosmo.differential_comoving_volume(u).value # potentially shouldn't be comoving volume
                return integrand
            comoving_volume_element, err = cf.scipy.integrate.quad(integrand, z1, z2)
            comoving_volume = comoving_volume_element * (cf.u.Mpc ** 3 / cf.u.sr) * (solid_angle * cf.u.deg ** 2).to(cf.u.sr) # in Mpc**3
            return comoving_volume
        # scale the histogram y-axis to form a number density
        # comoving_vol = cf.astropy_cosmo.comovingDistance(photo_z["min_z"], photo_z["max_z"]) \
        #     * field_geometry.masked_area
        comoving_vol = calculate_comoving_volume(photo_z["min_z"], photo_z["max_z"], field_geometry.masked_area.value)
        
        # calculate stellar mass function
        stellar_masses = stellar_masses[stellar_masses > 0] # remove galaxies with negative stellar masses
        weights = cf.np.full(len(stellar_masses), Nbins / comoving_vol)
        stellar_mass_func, stellar_mass_bin_edges = cf.np.histogram(cf.np.log10(stellar_masses), \
                                                                    bins = Nbins, weights = weights)
        stellar_mass_mid_bins = 0.5 * (stellar_mass_bin_edges[1:] + stellar_mass_bin_edges[:-1])
        
        return cls(photo_z, stellar_mass_mid_bins, stellar_mass_func)
    
    def plot_smf(self, colour = "black", stellar_mass_cut = 10):
        
        label = "z = %1.1f" % self.photo_z["mid_z"]
        cf.plt.plot(self.smf_bins, cf.np.log10(self.smf_values), label = label)
        cf.plt.axvline(x = stellar_mass_cut, color = "black", linestyle = "--")
        
    @staticmethod
    def plot_stellar_mass_functions(gal_field_obj_arr, plot_title, plot_literature = ["Mortlock_2015_z1_25"], save_fig = False):
        
        # define appropriate colours
        colours = []
        # plot the models
        for i in range(len(gal_field_obj_arr)):
            gal_field_obj_arr[i].stellar_mass_function.plot_smf()
            
        for author_year_z in plot_literature:
            if author_year_z == "Mortlock_2015_z1_25":
                stellar_mass_function.plot_literature_data("Mortlock15", 1.25, ["black"])
        
        # plot parameters
        #plot_title = "z = %1.1f-%1.1f galaxy bias" % (self.min_z, self.max_z)
        cf.plt.title(plot_title, fontsize = 14)
        cf.plt.xlabel('log$_{10}$(M$_{\star}$/M$_{\odot}$)', fontsize = 12)
        #cf.plt.ylabel("Arbitrary units (Havn't checked the units yet)")
        cf.plt.ylabel('log$_{10}$($\Phi$/Mpc$^{-3}$dex$^{-1}$)', fontsize = 12)
        cf.plt.xticks(fontsize = 10)
        cf.plt.yticks(fontsize = 10)
        cf.plt.xlim(7, 12)
        #cf.plt.yscale("log")
        #cf.plt.ylim()
        # legend handles
        cf.plt.legend(loc = "lower left", bbox_to_anchor = (0, 0))
        if save_fig == True:
            cf.os.chdir("/Users/user/Documents/PGR/UDS field/Plots/Galaxy bias")
            cf.plt.savefig(plot_title + ".jpeg", dpi=600)
        cf.plt.show()
        
    @staticmethod
    def plot_literature_data(author_year, redshift, colours, additional_info = None, ls = "dashed"):
        
        if additional_info == None:
            label = author_year + " z=" + str(redshift)
        else:
            label = author_year + " z=" + str(redshift) + " " + additional_info
        
        # open the .json file
        cf.os.chdir("/Users/user/Documents/PGR/Literature/Data/Mass function")
        with open(label + ".json") as jsonFile:
            jsonObject = cf.json.load(jsonFile)
            jsonFile.close()
        #print(jsonObject)
        # print out all data sets from the specific .json file
        for i in range(len(jsonObject['datasetColl'])):
            name = jsonObject['datasetColl'][i]['name']
            data = jsonObject['datasetColl'][i]['data']
            x_y_values = [data[index]['value'] for index in range(len(data))]
            x_values = cf.np.array(list(x_y_values[index][0] for index in range(len(x_y_values))))
            y_values = cf.np.log10(cf.np.array(list(x_y_values[index][1] for index in range(len(x_y_values)))))
            #print(x_values)
            #print(x_y_values)
            if author_year == "MM99":
                y_values = y_values - 1
            cf.plt.plot(x_values, y_values, c = colours[i], ls = ls, label = label + " " + name)  
        
if __name__ == "__main__":
    pass