#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 17:02:44 2022

@author: u92876da
"""

import config as cf

author_year_colours = {"This work": "black", "Foucaud 2010": "red"}
stellar_mass_bin_shapes_sizes_fill = {"M$_{\star}$ > 10$^{10.0}$ M$_{\odot}$": ["o", 7.5, "full"], \
                      "10$^{10.0}$ < M$_{\star}$ / M$_{\odot}$ < 10$^{10.5}$": ["s", 5, "none"], \
                      "10$^{10.5}$ < M$_{\star}$ / M$_{\odot}$ < 10$^{11.0}$": ["o", 7.5, "none"], \
                      "10$^{11.0}$ < M$_{\star}$ / M$_{\odot}$ < 10$^{11.5}$": ["o", 10, "full"], \
                      "10$^{11.0}$ < M$_{\star}$ / M$_{\odot}$ < 10$^{12.0}$": ["^", 10, "none"]}

# could be derived class (along with bias and SHMR, r_0 etc) from "resulting graph" base class    
    
class SHM_evolution:
    
    def __init__(self, z_data, z_data_errs, SHM_ratios, SHM_ratios_u1, SHM_ratios_l1, \
                 stellar_mass_min, stellar_mass_max, author_year):
        self.z_data = z_data
        self.z_data_errs = z_data_errs
        self.SHM_ratios = SHM_ratios
        self.SHM_ratios_u1 = SHM_ratios_u1
        self.SHM_ratios_l1 = SHM_ratios_l1
        self.stellar_mass_min = cf.np.round(stellar_mass_min, 1)
        self.stellar_mass_max = cf.np.round(stellar_mass_max, 1)
        self.author_year = author_year
        
    @classmethod
    def from_gal_field_arr(cls, gal_field_arr, author_year = "This work"):
        
        z_array = []
        z_err_array = []
        SHM_ratio_array = []
        SHM_ratio_err_array = []
        two_pt_angular_arr = [gal_field.w_theta_ACF for gal_field in gal_field_arr]
        for loc_object in two_pt_angular_arr:
            z_array = cf.np.append(z_array, loc_object.photo_z["mean_z"])
            z_err_array = cf.np.append(z_err_array, loc_object.photo_z["std_z"])
            SHM_ratio_array = cf.np.append(SHM_ratio_array, loc_object.SHM_ratio)
            SHM_ratio_err_array = cf.np.append(SHM_ratio_err_array, loc_object.SHM_ratio_err)
            
        # calculate min and max stellar masses    
        stellar_mass_min = cf.np.log10(two_pt_angular_arr[0].stellar_masses["min_Mstar"])
        if two_pt_angular_arr[0].sample_type.name == "LOWER_STELLAR_MASS_LIM":
            stellar_mass_max = 99
        else:
            stellar_mass_max = cf.np.log10(two_pt_angular_arr[0].stellar_masses["max_Mstar"])
        
        return cls(z_array, z_err_array, SHM_ratio_array, SHM_ratio_err_array, stellar_mass_min, \
                   stellar_mass_max, author_year)
        
    @classmethod
    def from_literature(cls, stellar_mass_min, stellar_mass_max, author_year = "Foucaud 2010", \
                        sample_type = "Stellar mass bin", filename = "Clustering data"):
        
        cf.os.chdir("/Users/user/Documents/PGR/Literature/Data")
        data = cf.pd.read_csv(filename + ".csv")
        # data = data.loc[lambda x: x["gal_type"] == sample_type]
        data = data.loc[lambda x: x["author_year"] == author_year]
        
        # if stellar_mass_min != None and stellar_mass_max != None:
        data = data.loc[lambda x: x["Mstar_min"] == stellar_mass_min]
        data = data.loc[lambda x: x["Mstar_max"] == stellar_mass_max]
        
        return cls(cf.np.array(data["z"]), cf.np.array(data["z_u1"]), cf.np.array(data["Mstar/M_DM"]), \
                   cf.np.array(data["Mstar/M_DM err"]), cf.np.array(data["Mstar/M_DM err"]), \
                       stellar_mass_min, stellar_mass_max, author_year = author_year)
            
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
    
    def plot_data(self):
        
        # plot the data
        # upper error bars
        cf.plt.errorbar(self.z_data, self.SHM_ratios, xerr = self.z_data_errs, yerr = self.SHM_ratios_u1, \
                        ls = "None", lolims = cf.np.full(len(self.z_data), True), \
                        fmt = stellar_mass_bin_shapes_sizes_fill[self.stellar_mass_label][0], \
                        ms = stellar_mass_bin_shapes_sizes_fill[self.stellar_mass_label][1], \
                  fillstyle = stellar_mass_bin_shapes_sizes_fill[self.stellar_mass_label][2], \
                          c = author_year_colours[self.author_year], capsize = 4, capthick = 1)
        # lower error bars    
        cf.plt.errorbar(self.z_data, self.SHM_ratios, xerr = self.z_data_errs, yerr = self.SHM_ratios_l1, \
                        ls = "None", uplims = cf.np.full(len(self.z_data), True), \
                        fmt = stellar_mass_bin_shapes_sizes_fill[self.stellar_mass_label][0], \
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
    def plot_horizontal_line(value, colour = "black", ls = "-"):
        cf.plt.axhline(y = value, color = colour, linestyle = ls)
    
    @staticmethod
    def plot_bias_evolution(SHM_evolution_obj_arr, plot_title, save_fig = False):
        
        # # plot the models
        # if models_to_plot != None:
        #     SHM_evolution.plot_models(models_to_plot)
        
        z_list = []
        SHM_ratio_list = []
        SHM_ratio_u1_list = []
        SHM_ratio_l1_list = []
        legend_handles = []
        for data in SHM_evolution_obj_arr:
            # plot the data
            handle = data.plot_data()
            legend_handles = cf.np.append(legend_handles, handle)
            for i in range(len(data.z_data)):
                z_list = cf.np.append(z_list, data.z_data[i])
                SHM_ratio_list = cf.np.append(SHM_ratio_list, data.SHM_ratios[i])
                SHM_ratio_u1_list = cf.np.append(SHM_ratio_u1_list, data.SHM_ratios_u1[i])
                SHM_ratio_l1_list = cf.np.append(SHM_ratio_l1_list, data.SHM_ratios_l1[i])

        # plot f_b * h
        SHM_evolution.plot_horizontal_line(cf.cosmo.Ob0 * cf.cosmo.H0 / (100 * cf.cosmo.Om0))
        # plot f_* * h, from Cole 2001 via Foucaud 2010
        SHM_evolution.plot_horizontal_line(2.03 * 1e-2 * cf.cosmo.H0 / 100, ls = "--")

        # plot parameters
        # plot_title = "z = %1.1f-%1.1f galaxy bias" % (self.min_z, self.max_z)
        cf.plt.title(plot_title, fontsize = 14)
        cf.plt.xlabel("z", fontsize = 12)
        cf.plt.ylabel("M$_{\star}$/M$_{DM}$ $h$", fontsize = 12)
        cf.plt.xticks(fontsize = 10)
        cf.plt.yticks(fontsize = 10)
        cf.plt.xlim(0, cf.np.max(z_list) + 1)
        cf.plt.ylim(0.5 * cf.np.min(SHM_ratio_list - SHM_ratio_l1_list), 1)#\
                    #1.5 * cf.np.max(SHM_ratio_list + SHM_ratio_u1_list))
        cf.plt.yscale("log")
        # legend handles
        cf.plt.legend(handles = list(legend_handles), loc = "upper center", \
                      bbox_to_anchor = (0.5, -0.15))
        if save_fig == True:
            cf.os.chdir("/Users/user/Documents/PGR/UDS field/Plots/SHM evolution")
            cf.plt.savefig(plot_title + ".jpeg", dpi=600)
        cf.plt.show()
   

def main():
    obj = SHM_evolution.from_literature(11, 12)
    SHM_evolution.plot_bias_evolution([obj], "test")
     
if __name__ == "__main__":
    main()