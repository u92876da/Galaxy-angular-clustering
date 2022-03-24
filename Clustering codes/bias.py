#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 22:12:31 2021

@author: u92876da
"""

# bias.py

import config as cf

# dicts for colours and shapes for different author_year and gal_types
author_year_colours = {"This work": "black", "Foucaud 2010": "red"}
stellar_mass_bin_shapes_sizes_fill = {"M$_{\star}$ > 10$^{10.0}$ M$_{\odot}$": ["o", 7.5, "full"], \
                      "10$^{10.0}$ < M$_{\star}$ / M$_{\odot}$ < 10$^{10.5}$": ["s", 5, "none"], \
                      "10$^{10.5}$ < M$_{\star}$ / M$_{\odot}$ < 10$^{11.0}$": ["o", 7.5, "none"], \
                      "10$^{11.0}$ < M$_{\star}$ / M$_{\odot}$ < 10$^{11.5}$": ["o", 10, "full"], \
                      "10$^{11.0}$ < M$_{\star}$ / M$_{\odot}$ < 10$^{12.0}$": ["^", 10, "none"]}


class bias_evolution:
    
    def __init__(self, z_data, z_data_errs, b_data, b_data_errs, stellar_mass_min, \
                 stellar_mass_max, author_year):
        
        self.z_data = z_data
        self.z_data_errs = z_data_errs
        self.b_data = b_data
        self.b_data_errs = b_data_errs
        #self.sample_type = sample_type
        self.stellar_mass_min = cf.np.round(stellar_mass_min, 1)
        self.stellar_mass_max = cf.np.round(stellar_mass_max, 1)
        self.author_year = author_year
     
    @classmethod
    def from_gal_field_arr(cls, gal_field_arr, author_year = "This work"):
        
        z_array = []
        z_err_array = []
        b_array = []
        b_err_array = []
        two_pt_angular_arr = [gal_field.w_theta_ACF for gal_field in gal_field_arr]
        for loc_object in two_pt_angular_arr:
            z_array = cf.np.append(z_array, loc_object.photo_z["mean_z"])
            z_err_array = cf.np.append(z_err_array, loc_object.photo_z["std_z"])
            b_array = cf.np.append(b_array, loc_object.b)
            b_err_array = cf.np.append(b_err_array, loc_object.b_err)
            
        # calculate min and max stellar masses    
        stellar_mass_min = cf.np.log10(two_pt_angular_arr[0].stellar_masses["min_Mstar"])
        if two_pt_angular_arr[0].sample_type.name == "LOWER_STELLAR_MASS_LIM":
            stellar_mass_max = 99
        else:
            stellar_mass_max = cf.np.log10(two_pt_angular_arr[0].stellar_masses["max_Mstar"])
        
        return cls(z_array, z_err_array, b_array, b_err_array, stellar_mass_min, stellar_mass_max, author_year)
    
    @classmethod
    def from_literature(cls, stellar_mass_min, stellar_mass_max, author_year = "Foucaud 2010", \
                        sample_type = "Stellar mass bin", filename = "gal_bias_literature"):
        
        cf.os.chdir("/Users/user/Documents/PGR/Literature")
        data = cf.pd.read_csv(filename + ".csv")
        data = data.loc[lambda x: x["gal_type"] == sample_type]
        data = data.loc[lambda x: x["author_year"] == author_year]
        
        # if stellar_mass_min != None and stellar_mass_max != None:
        data = data.loc[lambda x: x["Mstar min"] == stellar_mass_min]
        data = data.loc[lambda x: x["Mstar max"] == stellar_mass_max]
        
        return cls(cf.np.array(data["z_mid"]), cf.np.array(data["z_std"]), cf.np.array(data["bias"]), \
                   cf.np.array(data["bias_u1"]), stellar_mass_min, stellar_mass_max, author_year = author_year)
    
    @property
    def min_z(self):
        return cf.np.min(self.z_data)
    
    @property
    def max_z(self):
        return cf.np.max(self.z_data)
     
    @property
    def stellar_mass_label(self):
        if self.stellar_mass_max == 99:
            stellar_mass_label = "M$_{\star}$ > 10$^{%1.1f}$ M$_{\odot}$" % self.stellar_mass_min
        else:
            stellar_mass_label = "10$^{%1.1f}$ < M$_{\star}$ / M$_{\odot}$ < 10$^{%1.1f}$" \
                 % (self.stellar_mass_min, self.stellar_mass_max)
        return stellar_mass_label
        
    def calc_no_evolution_model(self): # B0 in MM2000, give options to fit or give parameters and determine red_chi_sq
    
        def const_fit(z, b_0):
            return b_0
        popt, pcov = cf.scipy.optimize.curve_fit(const_fit, self.z_data, self.b_data, \
                                                 sigma = self.b_data_errs, absolute_sigma = True)
        # compute reduced chi squared
        model = cf.np.full(len(self.z_data), const_fit(self.z_data, *popt))
        chi_sq = 0
        for i in range(len(self.z_data)):
            chi_sq += ((self.b_data[i] - model[i]) / self.b_data_errs[i]) ** 2
        red_chi_sq = chi_sq / (len(self.z_data) - len(popt))
        
        # recalculate model with more bins
        long_z_array = cf.np.linspace(self.min_z - 1, self.max_z + 1, 1000)
        model = cf.np.full(len(long_z_array), const_fit(long_z_array, *popt))
        
        # print output parameters
        print("No evolution model: [b_0, \u03C7^2 red.] = [%1.2f \u00B1 %1.2f, %1.2f]" \
              % (popt[0], cf.np.sqrt(pcov[0][0]), red_chi_sq))
        
        return model, popt, pcov, red_chi_sq
    
    def calc_test_particle_bias_model(self): # B1 in MM2000, galaxy conserving model
        
        def test_particle_bias_fit(z, b_0):
            return 1 + (b_0 - 1) / cf.cosmo.growthFactor(z)
        popt, pcov = cf.scipy.optimize.curve_fit(test_particle_bias_fit, self.z_data, self.b_data, \
                                                 sigma = self.b_data_errs, absolute_sigma = True)
        # compute reduced chi squared
        model = test_particle_bias_fit(self.z_data, *popt)
        chi_sq = 0
        for i in range(len(self.z_data)):
            chi_sq += ((self.b_data[i] - model[i]) / self.b_data_errs[i]) ** 2
        red_chi_sq = chi_sq / (len(self.z_data) - len(popt))
        
        # recalculate model with more bins
        long_z_array = cf.np.linspace(self.min_z - 1, self.max_z + 1, 1000)
        model = test_particle_bias_fit(long_z_array, *popt)
        
        # print output parameters
        print("Test particle bias model: [b_0, \u03C7^2 red.] = [%1.2f \u00B1 %1.2f, %1.2f]" \
              % (popt[0], cf.np.sqrt(pcov[0][0]), red_chi_sq))
        
        return model, popt, pcov, red_chi_sq
    
    @staticmethod
    def calc_transient_model(M_min): # B2ii in MM2000
        M_min_params = cf.pd.read_csv(cf.M_min_params_file)
        M_min_params = M_min_params[cf.np.log10(M_min_params["M_min"]) == M_min]
        b0_eff_z0 = M_min_params["b0_eff_z0"]
        beta = M_min_params["beta"]
        
        # calculate model with more bins
        long_z_array = cf.np.linspace(0, 10, 1000)
        D_z = cf.np.array(cf.cosmo.growthFactor(long_z_array))
        gal_bias = 0.41 + (cf.np.full(len(D_z), b0_eff_z0) - 0.41) / cf.np.power(D_z, cf.np.full(len(D_z), beta))
        return gal_bias
    
    def calc_merging_model(self): # B2i in MM2000
        pass
    
    @staticmethod
    def plot_models(models_to_plot = [], plot_labels = False, linestyles = None, colours = None):
        
        long_z_array = cf.np.linspace(0, 10, 1000)
        if colours == None:
            colours = []
            for i in range(len(models_to_plot)):
                colours.append(cf.plasma(i * 0.9 / len(models_to_plot[:-1])))
        
        for i in range(len(models_to_plot)):
            name = models_to_plot[i]
            
            if name == "No evolution":
                pass
                #model, popt, pcov, red_chi_sq = self.calc_no_evolution_model()

            elif name == "Test particle bias":
                pass
                #model, popt, pcov, red_chi_sq = self.calc_test_particle_bias_model()

            elif name.split(" ")[1] == "Transient": # based on Press, Schechter 1974 dn/dM and Mo, White 1996 b(z)
                M_min = float(name.split(" ")[0])
                model = bias_evolution.calc_transient_model(M_min)
                plot_name = "M$_{\mathrm{min}}$=10$^{%1.1f}h^{-1}$M$_{\odot}$" %(M_min)
                #cf.plt.text(long_z_array[850] + 0.3, model[850] - 0.8, "10$^{%2.0f}$" %M_min, c = colour)
            
            if plot_labels == True:
                if name.split(" ")[1] == "Transient":
                    label = plot_name
                else:
                    label = name + " model"   
            else:
                label = None
                
            cf.plt.plot(long_z_array, model, c = colours[i], label = label)
        
    def plot_data(self):
        
        # plot the data
        cf.plt.errorbar(self.z_data, self.b_data, xerr = self.z_data_errs, yerr = self.b_data_errs, \
                        ls = "None", fmt = stellar_mass_bin_shapes_sizes_fill[self.stellar_mass_label][0], \
                        ms = stellar_mass_bin_shapes_sizes_fill[self.stellar_mass_label][1], \
                 fillstyle = stellar_mass_bin_shapes_sizes_fill[self.stellar_mass_label][2], \
                         c = author_year_colours[self.author_year], capsize = 4, capthick = 1)

        # legend
        legend_label = self.author_year + ", " + self.stellar_mass_label 
        legend_handle = cf.Line2D([0], [0], color = author_year_colours[self.author_year], \
                                  marker = stellar_mass_bin_shapes_sizes_fill[self.stellar_mass_label][0], \
                                      ms = stellar_mass_bin_shapes_sizes_fill[self.stellar_mass_label][1], \
                               fillstyle = stellar_mass_bin_shapes_sizes_fill[self.stellar_mass_label][2], \
                                   label = legend_label)
        return legend_handle
        
    @staticmethod
    def plot_bias_evolution(bias_evolution_obj_arr, plot_title, save_fig = False, \
                            models_to_plot = ["10 Transient", "11 Transient", "12 Transient", "13 Transient"]):
        
        # plot the models
        if models_to_plot != None:
            bias_evolution.plot_models(models_to_plot)
            
        z_list = []
        b_list = []
        b_errs_list = []
        legend_handles = []
        for bias_data in bias_evolution_obj_arr:
            # plot the data
            handle = bias_data.plot_data()
            legend_handles = cf.np.append(legend_handles, handle)
            for i in range(len(bias_data.z_data)):
                z_list = cf.np.append(z_list, bias_data.z_data[i])
                b_list = cf.np.append(b_list, bias_data.b_data[i])
                b_errs_list = cf.np.append(b_errs_list, bias_data.b_data_errs[i])
        # plot parameters
        #plot_title = "z = %1.1f-%1.1f galaxy bias" % (self.min_z, self.max_z)
        cf.plt.title(plot_title, fontsize = 14)
        cf.plt.xlabel("z", fontsize = 12)
        cf.plt.ylabel("b(z)", fontsize = 12)
        cf.plt.xticks(fontsize = 10)
        cf.plt.yticks(fontsize = 10)
        cf.plt.xlim(0, cf.np.max(z_list) + 1)
        cf.plt.ylim(0, 1.2 * cf.np.max(b_list + b_errs_list))
        # legend handles
        cf.plt.legend(handles = list(legend_handles), loc = "upper center", \
                      bbox_to_anchor = (0.5, -0.15))
        if save_fig == True:
            cf.os.chdir("/Users/user/Documents/PGR/UDS field/Plots/Galaxy bias")
            cf.plt.savefig(plot_title + ".jpeg", dpi=600)
        cf.plt.show()
        
    @staticmethod
    def calc_bias(z, b0, b1, b2):
        b = b0 * cf.np.power((1 + z), b1) + b2
        return b
     
def main():
    bias_evolution.from_literature()
     
if __name__ == "__main__":
    main()
        
            