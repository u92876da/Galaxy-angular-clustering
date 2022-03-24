# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 19:31:52 2021

@author: dunca
"""

# two_pt_angular_correlation

import config as cf
from redshift_distribution import redshift_distribution
from field_geometry import field_geometry
from cosmic_var import cosmic_var

class two_pt_angular_corr:
    
    def __init__(self, fieldname, theta_mid_bins, w_theta, w_theta_err, bootstraps, \
                 galaxy_redshifts, galaxy_stellar_masses, sample_type, random_gals_ra_dec, delta = 0.8):
        
        # construct class member data
        self.fieldname = fieldname
        self.bootstraps = bootstraps
        self.galaxy_redshifts = galaxy_redshifts
        self.galaxy_stellar_masses = galaxy_stellar_masses
        self.sample_type = sample_type
        self.random_gals_ra_dec = random_gals_ra_dec
        
        # ensure entered data is not nan
        indices_to_remove = list([i for i in range(len(w_theta)) if cf.np.isnan(w_theta[i])])
        if len(indices_to_remove) != 0:
            print("z = " + str(self.photo_z["mid_z"]) + " indices removed = " \
                  + str(indices_to_remove) + ", Ngals = " + str(self.Ngals))
                
        # save numerical data that is not nan
        self.theta_mid_bins = cf.np.array([theta_mid_bins[i] for i in range(len(w_theta)) if i not in indices_to_remove])
        self.w_theta = cf.np.array([w_theta[i] for i in range(len(w_theta)) if i not in indices_to_remove])
        self.w_theta_err = cf.np.array([w_theta_err[i] for i in range(len(w_theta)) if i not in indices_to_remove])
        
        self.delta = delta
        self.delta_err = 0
        
        # integral constraint
        C = self.C_as_function_of_Nrandom(plot_graph = False)
        self.C = C[0]
        self.C_err = C[1]
        self.C_fit_red_chi_sq = C[2]
        
        self.calculate_A()
        self.calculate_b(flat_N_z = True)
        #self.cosmic_variance()
        # self.calculate_b(flat_N_z = False) # repeat this calculation many times
        
        # include cosmic variance and errors from adjacent bin scattering here
        self.A_err_cv = cf.np.sqrt(self.cosmic_var.unmasked_gal_cosmic_var \
                                   (self.photo_z["z_width"])) * self.A
        print(self.A_err * 1e3)
        print(self.A_err_cv * 1e3)
        self.A_err = cf.np.sqrt(self.A_err ** 2 + self.A_err_cv ** 2)
        self.calculate_b(flat_N_z = True)
        
        # remove integral constraint
        self.w_theta_true = self.w_theta + self.C * self.A
        self.w_theta_true_err = cf.np.sqrt(self.w_theta_err ** 2 + (self.C * self.A_err) ** 2 + \
                                           (self.C_err * self.A) ** 2)
        
    @classmethod # include galaxy redshifts here for now
    def from_pkl(cls, w_theta_pkl_name):
        cf.os.chdir("/Users/user/Documents/PGR/UDS field/pkl")
        with open(w_theta_pkl_name + '.pkl', 'rb') as f:
            fieldname, theta, w_theta, w_theta_err, bootstraps, galaxy_redshifts, \
                galaxy_stellar_masses, sample_type, random_gals_ra_dec = cf.pickle.load(f)
        return cls(fieldname, theta, w_theta, w_theta_err, bootstraps, galaxy_redshifts, \
                   galaxy_stellar_masses, sample_type, random_gals_ra_dec)
            
    # -----------------------------------------------------------------------------------
        
    @property
    def Ngals(self):
        return len(self.galaxy_redshifts)
    
    @property
    def photo_z(self):
        min_z = cf.np.min(self.galaxy_redshifts)
        max_z = cf.np.max(self.galaxy_redshifts)
        photo_z = {"mean_z": cf.np.mean(self.galaxy_redshifts), "std_z": \
                cf.np.std(self.galaxy_redshifts), "min_z": min_z, "max_z": max_z, \
                "z_width": max_z - min_z, "mid_z": cf.np.round((max_z - min_z) / 2 + min_z, 1)} 
        return photo_z    
    
    @property
    def stellar_masses(self):
        min_stellar_mass = cf.np.round(cf.np.min(self.galaxy_stellar_masses), 1)
        max_stellar_mass = cf.np.round(cf.np.max(self.galaxy_stellar_masses), 1)
        return {"mean_Mstar": cf.np.mean(self.galaxy_stellar_masses), "std_Mstar": \
                cf.np.std(self.galaxy_stellar_masses), "min_Mstar": min_stellar_mass, \
                "max_Mstar": max_stellar_mass, "mid_Mstar": \
                cf.np.round((max_stellar_mass - min_stellar_mass) / 2 + min_stellar_mass, 1)}  
    
    @property
    def Nbootstraps(self):
        return len(self.bootstraps)
    
    @property
    def theta_params(self):
        delta_theta = self.theta_mid_bins[1] - self.theta_mid_bins[0]
        theta_min_bin = cf.np.min(self.theta_mid_bins)
        theta_max_bin = cf.np.max(self.theta_mid_bins)
        return {"theta_min_bin": theta_min_bin, "theta_max_bin": theta_max_bin, "delta_theta": delta_theta, \
                "theta_min": theta_min_bin - delta_theta / 2, "theta_max": theta_max_bin + delta_theta / 2}
    
    @property
    def fiducial_hm(self):
        return cf.AngularCF(z = self.z, zmin = self.z_min, zmax = self.z_max, \
                    theta_min = (self.theta_min * cf.u.deg).to(cf.u.rad).value, \
                    theta_max = (self.theta_max * cf.u.deg).to(cf.u.rad).value, \
                    theta_log = True, theta_num = len(self.theta_mid_bins)) #, hod_model = HOD_model)

    @property
    def Nrandom_pairs(self):
        cf.os.chdir("/Users/user/Documents/PGR/UDS field")
        return cf.pd.read_csv(self.fieldname + " N random pairs.csv")
    
    @property
    def cosmic_var(self):
        return cosmic_var(self.fieldname, self.photo_z["mid_z"], self.b)

    def calc_z_distribution(self, Nbins = 10):
        #self.z_distribution = redshift_distribution(self.fieldname, self.galaxy_redshifts, Nbins, gal_bias)
        return redshift_distribution(self.fieldname, self.galaxy_redshifts, Nbins, self.b)

    # plot the data and models ----------------------------------------------------------

    @staticmethod
    def plot_two_pt_angular_corr(w_theta_obj_array, fit_power_law_models = True, fit_halo_models = False, \
        label_plot = False, colours = None, shapes = None, vline = 1, save_fig = False, plot_literature = []):
        # define appropriate colours from plasma colourmap
        if colours == None:
            colours = []
            for i in range(len(w_theta_obj_array)):
                colours.append(cf.plasma(i * 0.9 / len(w_theta_obj_array)))
        # default shapes are squares
        if shapes == None:
            shapes = ["s" for x in range(len(w_theta_obj_array))]
        
        for i in range(len(w_theta_obj_array)):
            w_theta_obj = w_theta_obj_array[i]
            w_theta_obj.plot_data(colours[i], shape = shapes[i], vline = vline)
            if fit_power_law_models == True:
                w_theta_obj.plot_power_law_model(colours[i], label_plot = label_plot)
            elif fit_halo_models == True:
                w_theta_obj.plot_halo_model(colours[i], label_plot = label_plot, params = "fiducial")
        
        # determine plot title
        if len(w_theta_obj_array) == 1:
            redshift_label = "z = %1.f, " % w_theta_obj_array[0].photo_z["mid_z"]
        else:
            z_list = list([w_theta_obj.photo_z["mid_z"] for w_theta_obj in w_theta_obj_array])
            redshift_label = "z = %1.1f-%1.1f, " % (cf.np.min(z_list), cf.np.max(z_list))
        if w_theta_obj_array[0].sample_type.name == "ALL":
            mass_label = "All galaxies, "
            saving_mass_label = mass_label
        elif w_theta_obj_array[0].sample_type.name == "LOWER_STELLAR_MASS_LIM": # stellar masses here
            mass_label = "M$_{\star}$ > 10$^{%1.1f}$M$_{\odot}$, " % cf.np.log10(w_theta_obj_array[0].stellar_masses["min_Mstar"])
            saving_mass_label = "log(Mstar) > %1.1f, " % cf.np.log10(w_theta_obj_array[0].stellar_masses["min_Mstar"])
        else:
            mass_label = "10$^{%1.1f}$ < M$_{\star}$/M$_{\odot}$ < 10$^{%1.1f}$, " \
                % (cf.np.log10(w_theta_obj_array[0].stellar_masses["min_Mstar"]), \
                   cf.np.log10(w_theta_obj_array[0].stellar_masses["max_Mstar"]))
            saving_mass_label = "%1.1f < log(Mstar) < %1.1f, " \
                            % (cf.np.log10(w_theta_obj_array[0].stellar_masses["min_Mstar"]), \
                               cf.np.log10(w_theta_obj_array[0].stellar_masses["max_Mstar"]))
        if fit_power_law_models == True:
            delta = w_theta_obj_array[0].delta
            delta_label = "δ = {}".format(delta)
        else:
            delta_label = None
            
        for author_year_z in plot_literature:
            if author_year_z == "Foucaud_2010_z1":
                two_pt_angular_corr.plot_literature_data("Foucaud10", 1, ["red", "purple", "green"])
            
        # plot parameters    
        plot_title = redshift_label + mass_label + delta_label
        cf.plt.title(plot_title, fontsize = 14)
        cf.plt.xlabel("θ /deg", fontsize = 12)
        cf.plt.ylabel("w(θ)", fontsize = 12)
        lgd = cf.plt.legend(loc="upper center", frameon=True, fontsize=10, bbox_to_anchor = (0.5, -0.2), fancybox=True, shadow=True)
        cf.plt.xticks(fontsize = 10)
        cf.plt.yticks(fontsize = 10)
        cf.plt.xscale("log")
        cf.plt.yscale("log")
        cf.plt.xlim(w_theta_obj_array[0].theta_params["theta_min"] * 0.95, \
                    w_theta_obj_array[0].theta_params["theta_max"] * 1.1)
        if save_fig == True:
            cf.os.chdir("/Users/user/Documents/PGR/UDS field/Plots/2pt angular correlation")
            save_title = redshift_label + saving_mass_label + delta_label + " w(θ).jpeg"
            cf.plt.savefig(save_title, dpi = 600, box_extra_artists = (lgd,), bbox_inches = "tight")
        cf.plt.show()

    def plot_data(self, colour, shape = "s", vline = 1):
        
        # plot w_theta including errors
        cf.plt.errorbar(self.theta_mid_bins, self.w_theta_true, yerr = self.w_theta_true_err, \
                      label = ("z = %1.1f, N = %5.0f" \
                    % (self.photo_z["mid_z"], self.Ngals)), c = colour, fmt = shape,\
                      capsize = 5, capthick = 1, ls = "none") # , M$_{\star, min}$ = $10^{%1.1f}$
        
        # plot vertical n Mpc/h line
        if vline != None: 
            theta_Mpc_h = (cf.astropy_cosmo.arcsec_per_kpc_comoving(self.photo_z["mid_z"])).to(cf.u.deg / cf.u.Mpc) * vline
            cf.plt.axvline(x = theta_Mpc_h.value, color = colour, linestyle = "--")

    def plot_power_law_model(self, colour, theta_bins = 1000, label_plot = False, model_name = None):

        #theta = 10 ** cf.np.linspace(cf.np.log10(cf.np.min(self.theta_edge_bins)), cf.np.log10(cf.np.max(self.theta_edge_bins)), 1000) 
        theta = 10 ** cf.np.linspace(cf.np.log10(self.theta_params["theta_min"]), \
                                     cf.np.log10(self.theta_params["theta_max"]), theta_bins)
        model = self.A * (theta ** (- self.delta))
        # model name here
        if label_plot == True:
            plot_label = ("z = {} power law".format(self.z))
        else:
            plot_label = None
        
        cf.plt.plot(theta, model, c = colour, label = plot_label)

    def plot_halo_model(self, colour, HOD_model = "Zehavi05", params = "fiducial", label_plot = True):
        
        hm = cf.AngularCF(z = self.z, zmin = self.z_min, zmax = self.z_max, theta_min = (self.theta_min * cf.u.deg).to(cf.u.rad).value,\
                                  theta_max = (self.theta_max * cf.u.deg).to(cf.u.rad).value, theta_log = True, hod_model = HOD_model)
            
        if params != "fiducial":
            hm.hod_params = params
        
        plot_label = "z = %1.1f " %self.z + HOD_model
        cf.plt.plot((hm.theta * cf.u.rad).to(cf.u.deg).value, hm.angular_corr_gal, c = colour, label = plot_label)

    @staticmethod
    def plot_literature_data(author_year, redshift, colours, additional_info = None, ls = "dashed"):
        
        if additional_info == None:
            label = author_year + " z=" + str(redshift)
        else:
            label = author_year + " z=" + str(redshift) + " " + additional_info
        
        # open the .json file
        cf.os.chdir("/Users/user/Documents/PGR/Literature/Data/2pt angular correlation")
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
            y_values = cf.np.array(list(x_y_values[index][1] for index in range(len(x_y_values))))
            #print(x_values)
            #print(x_y_values)
            if author_year == "MM99":
                y_values = y_values - 1
            cf.plt.plot(x_values, y_values, c = colours[i], ls = ls, label = label + " " + name)  

    # put these in a new integral_constraint.py file ------------------------------------

    #@classmethod
    def calc_N_random_pairs(cls, min_max_ra_dec, Nrandom = 1000, Nrepeats = 5, re_initialize = False):
        
        if re_initialize == True:
            cls.initialize_N_random_pairs()
        
        for repeat in range(Nrepeats): # accomodate different random fields
            # setup random field
            cf.np.random.seed() # set random number seed
            random_ra = cf.np.random.uniform(min_max_ra_dec["min_ra"], min_max_ra_dec["max_ra"], size = Nrandom)
            cf.np.random.seed() # set different random number seed
            random_dec = cf.np.random.uniform(min_max_ra_dec["min_dec"], min_max_ra_dec["max_dec"], size = Nrandom)
            random_field = cf.SkyCoord(ra = random_ra * cf.u.deg, dec = random_dec * cf.u.deg)
            
            # calculate the number of pairs of random galaxies in each theta bin
            sep2d = random_field.search_around_sky(random_field, 100 * cf.u.deg)[2]
            sep2d = [x for x in sep2d.value if x != 0] # remove all occurrences of 0
            # sum extends to largest separations in the field
            bins = 10 ** cf.np.linspace(cf.np.log10(1e-5), cf.np.log10(1e1), 1000)
            Npairs, bin_edges = cf.np.histogram(sep2d, bins = bins)
            Npairs = Npairs / 2 # avoid double counting
            
            cls.N_random_pairs.insert(1, "run={}".format(repeat + 1),\
                                        cf.np.append(Nrandom, Npairs), allow_duplicates = True)
        print("Nrandom={} calculated".format(Nrandom))
        return cls.N_random_pairs
        
    #@classmethod
    def initialize_N_random_pairs(cls):
        
        if hasattr(cls, "__N_random_pairs"):
            del cls.N_random_pairs

        # initialize class attribute
        bins = 10 ** cf.np.linspace(cf.np.log10(1e-5), cf.np.log10(1e1), 1000)
        theta_bins = 0.5 * (bins[1:] + bins[:-1])
        N_random_random_bins = cf.pd.DataFrame(data = [theta_bins], index = ["theta_bins"])
        N_random_random_bins.insert(0, "Nrandom", [None])
        cls.N_random_pairs = (N_random_random_bins).T

    def calc_integral_constraint(self, Nrandom = 1000):
        
        theta_bins = self.Nrandom_pairs["theta_bins"][1:]

        total_pairs = Nrandom * (Nrandom - 1) / 2
        
        loc_random_pairs = ((self.Nrandom_pairs.T).loc[lambda x: x[0] == Nrandom]).T
        C_array = []
        
        for repeat in range(len(loc_random_pairs.T)): 
            #self.Nrandom_pairs.T
            Npairs = loc_random_pairs.iloc[1:, repeat]
            
            # calculate C using Roche and Eales formula
            C = cf.np.sum(Npairs * theta_bins ** (- self.delta)) / total_pairs
            C_array = cf.np.append(C_array, C)
                        
        C_mean = cf.np.round(cf.np.mean(C_array), 3)
        C_std = cf.np.round(cf.np.std(C_array), 3)

        return [C_mean, C_std]
    
    def C_as_function_of_Nrandom(self, Nrandom_arr = [200, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 10000],\
                                 shape = "s", plot_graph = True, save_fig = False):
        C_mean = []
        C_std = []
        for Nrandom in Nrandom_arr:
            #cls.calc_N_random_random(Nrepeats, Nrandom)#, Nbins, min_bin)
            C_loc = self.calc_integral_constraint(Nrandom)
            C_mean = cf.np.append(C_mean, C_loc[0])
            C_std = cf.np.append(C_std, C_loc[1])
        
        # least squared fit to data
        def C_fitting_func(theta, const):
            return const
        
        # least squares fitting programme
        # perform least squared fit
        popt, pcov = cf.curve_fit(C_fitting_func, Nrandom_arr, C_mean, sigma = C_std, absolute_sigma = True, method = "lm")
        model = cf.np.full(len(Nrandom_arr), C_fitting_func(Nrandom_arr, *popt))
        # compute the reduced chi squared value
        chi_sq = 0
        for i in range(len(Nrandom_arr)):
            chi_sq += ((C_mean[i] - model[i]) / C_std[i]) ** 2
        red_chi_sq = chi_sq / (len(Nrandom_arr) - len(popt))
        
        if plot_graph == True:
            # plot data
            cf.plt.errorbar(Nrandom_arr, C_mean, yerr = C_std, c = "black", fmt = shape, capsize = 5, capthick = 1, ls = "none")
            # plot model
            cf.plt.plot(Nrandom_arr, model, color = "red", label = "C = %1.3f \u00B1 %1.3f, \u03C7^2 red. = %1.2f"\
                        % (popt[0], cf.np.sqrt(pcov[0][0]), red_chi_sq))
            
            # plot parameters
            plot_title = "C (N$_{random}$), \u03B4 = %1.1f" % self.delta
            cf.plt.title(plot_title, fontsize = 14)
            cf.plt.xlabel("N$_{random}$", fontsize = 12)
            cf.plt.ylabel("C", fontsize = 12)
            cf.plt.legend(loc="upper right", frameon=True, fontsize=10, bbox_to_anchor=(1, 1), fancybox=True, shadow=True)
            cf.plt.xticks(fontsize = 10)
            cf.plt.yticks(fontsize = 10)
            #cf.plt.xscale("log")
            #cf.plt.yscale("log")
            if save_fig == True:
                cf.os.chdir("/Users/user/Documents/PGR/UDS field/Plots/2pt angular correlation")
                cf.plt.savefig("C (N_random), \u03B4 = %1.1f" % self.delta + ".jpeg", dpi=600) 
            cf.plt.show()
        
        # self.C = popt[0]
        # self.C_err = cf.np.sqrt(pcov[0][0])
        # self.C_fit_red_chi_sq = red_chi_sq
        
        return popt[0], cf.np.sqrt(pcov[0][0]), red_chi_sq
        
    # calculate best fitting power law parameters ---------------------------------------
    
    def calculate_A(self, print_params = False): # A_err = None, delta_err = None, C_err = None):

        def two_pt_angular_fitting_func(theta, A_loc):
            return A_loc * (theta ** (- self.delta) - self.C)
        
        # least squares fitting programme
        # perform least squared fit
        popt, pcov = cf.curve_fit(two_pt_angular_fitting_func, self.theta_mid_bins,\
                self.w_theta, sigma = self.w_theta_err, absolute_sigma = True, method = "lm")
        model = two_pt_angular_fitting_func(self.theta_mid_bins, *popt)
            
        # compute the reduced chi squared value
        chi_sq = 0
        for i in range(len(self.theta_mid_bins)):
            chi_sq += ((self.w_theta[i] - model[i]) / self.w_theta_err[i]) ** 2
        red_chi_sq = chi_sq / (len(self.theta_mid_bins) - len(popt))
        
        # store the fitting parameters in the object
        self.A = popt[0]
        self.A_err = cf.np.sqrt(pcov[0][0])
        self.A_fit_red_chi_sq = red_chi_sq
        # params = [popt[0], delta, C[0]]
        # param_errs = [cf.np.sqrt(pcov[0][0]), 0, C[1]]
        # self.A_fit = list(zip(params, param_errs)), red_chi_sq
        if print_params == True:
            self.print_model_params(print_A_delta_C = True)
            
    def calculate_A_delta(self): # iterative MCMC
        self.A_delta_fit = 0
        self.print_model_params("A, delta")
    
    # functions to define the parameters to be used and print these out -----------------
    
    # def define_fitting_params(self, fitting_params):
    #     if fitting_params == "A":
    #         params = self.A_fit
    #     elif fitting_params == "A, delta":
    #         params = self.A_delta_fit
    #     return params
    
    # this function could do with a bit of neatening
    def print_model_params(self, print_A_delta_C = True, print_r0_s8gal_b = True):
        
        if print_A_delta_C == True:
            print('z = %1.1f: [A, \u03B4, C, \u03C7^2 red.] = [%1.2e \u00B1 %1.2e, %5.3f \u00B1 %5.3f, ' \
                '%5.3f \u00B1 %5.3f, %5.3f]' % (self.photo_z["mid_z"], self.A, self.A_err, self.delta, \
                self.delta_err, self.C, self.C_err, self.A_fit_red_chi_sq))
            
        if print_r0_s8gal_b == True:
            print('z = %1.1f: [r_0 / Mpc/h, σ_8_gal, b] = [%1.2f \u00B1 %1.2f, %1.2f \u00B1 %1.2f, %1.2f \u00B1 %1.2f]' \
                  % (self.photo_z["mid_z"], self.r_0.value, self.r_0_err.value, self.sigma_8_gal, self.sigma_8_gal_err, \
                     self.b, self.b_err))

    # def write_two_pt_data(self):
        
    #     filename = "z = {}, w_theta data".format(self.z)
        
    #     column_labels = ["theta", "w_theta", "w_theta_err"]
    #     data = cf.np.transpose([self.theta_mid_bins, self.w_theta, self.w_theta_err])
    #     print(data)
        
    #     cf.os.chdir("/Users/user/Documents/PGR/UDS field/Plots/2pt angular correlation/Data/")
    #     write_data = cf.pd.DataFrame(data = data, columns = column_labels)
    #     write_data.to_csv(filename + ".csv", mode='w')
    #     print("Written to \'" + filename + "\'")

    # calculate clustering lengths, biases and halo masses ------------------------------ 
    
    def calculate_r_0(self, method = "Adel05", flat_N_z = False, Nbins = 10): # for a flat Universe with non-zero cosmological const
        
        # if A == None and A_err == None and delta == None and z == None and z_width == None:
        #     A = self.A
        #     A_err = self.A_err
        #     delta = self.delta
        #     z = self.photo_z["mid_z"]
        #     z_width = self.photo_z["max_z"] - self.photo_z["min_z"]
        #     z_min = self.photo_z["min_z"]
        #     z_max = self.photo_z["max_z"]
        # elif A != None and A_err != None and delta != None and z != None and z_width != None:
        #     z_min = z - (z_width / 2)
        #     z_max = z + (z_width / 2)
        # else:
        #     raise ValueError("Must specify either all or neither of (A, A_err, delta, z, z_width)!")
    
        #print("r_0 method: " + method)
        
        H_gamma = cf.math.gamma(1/2) * cf.math.gamma(self.delta / 2) / cf.math.gamma((self.delta + 1) / 2)
        #print("H_gamma = %1.2f" %H_gamma)
        F_z = 1 # curvature correction factor, from MM99
        
        # N(z) linear fit function parameters
        if flat_N_z == False:
            model, popt, pcov, red_chi_sq = self.calc_z_distribution(self.b, Nbins).lin_fit_params
        
        # N_z_dz_integral
        if flat_N_z == True:
            N_z_dz_integral = 1
        else:
            def N_z_dz_integrand(u):
                return popt[0] + (u * popt[1])
            N_z_dz_integral, numerical_err = cf.scipy.integrate.quad(N_z_dz_integrand, \
                                            self.photo_z["min_z"], self.photo_z["max_z"])
            # error propagation
            # N_z_dz_A_err = (z_max - z_min) * cf.np.sqrt(pcov[0][0])
            # N_z_dz_B_err = 0.5 * ((z_max ** 2) - (z_min ** 2)) * cf.np.sqrt(pcov[1][1])
            # #print("N_z_dz_A_err = {}, N_z_dz_B_err = {}".format(N_z_dz_A_err * 100 / N_z_dz_integral, \
            # #                                N_z_dz_B_err * 100 / N_z_dz_integral))
            # N_z_dz_err = cf.np.sqrt((N_z_dz_A_err ** 2) + (N_z_dz_B_err ** 2))
        
        if method == "MM99":
            def integrand(u):
                x_z = cf.cosmo.comovingDistance(z_min = 0, z_max = u) * cf.u.Mpc / cf.u.littleh # comoving distance
                P_Om0_z = cf.np.sqrt(cf.cosmo.Om0 * ((1 + u) ** 3 + (1 / cf.cosmo.Om0) - 1))
                return (x_z ** (- self.delta) * P_Om0_z * F_z).value
            integral = cf.scipy.integrate.quad(integrand, self.photo_z["min_z"], self.photo_z["max_z"])[0] \
                * ((cf.u.Mpc / cf.u.littleh) ** (- self.delta))
            r_0 = (cf.const.c * self.A * (self.photo_z["z_width"] ** 2) / (100 * \
                    (cf.u.littleh * cf.u.km / cf.u.s / cf.u.Mpc) * H_gamma * integral)) ** (1 / (1 + self.delta))
            #print(P_Om0_z) # checked!
            #rms_z_err = 0.1 * (1 + z)
            #z_bin_width_increase = cf.np.sqrt(12 * cf.np.power(rms_z_err / (z_width), 2) + 1) # assuming Gaussian errors
            #z_width = z_width / z_bin_width_increase # NOT SURE IF THIS IS 100% CORRECT
            # REQUIRES APPENDING!!!
            #r_0 = r_0 / (1 + z)
            r_0_err = (r_0 / (1 + self.delta)) * self.A_err / self.A
            
        elif method == "Efstathiou91":
            def integrand(u):
                epsilon = -1.2
                z_dependence = (1 + u) ** (- (3 + epsilon - (1 + self.delta)))
                x_z = cf.cosmo.comovingDistance(z_min = 0, z_max = u) * cf.u.Mpc / cf.u.littleh # comoving distance
                dz_dx = 1 / cf.np.sqrt(cf.cosmo.Om0 * ((1 + u) ** 3 + (1 / cf.cosmo.Om0) - 1))
                print((x_z ** (- self.delta) * dz_dx * F_z * z_dependence).unit)
                return (x_z ** (- self.delta) * dz_dx * F_z * z_dependence).value
            integral = cf.scipy.integrate.quad(integrand, self.photo_z["min_z"], self.photo_z["max_z"])[0] \
                * ((cf.u.Mpc / cf.u.littleh) ** (- self.delta))
            print(integral.unit)
            r_0 = (self.A / (H_gamma * integral)) ** (1 / (1 + self.delta))
            print(r_0.unit)
            r_0_err = (r_0 / (1 + self.delta)) * self.A_err / self.A
        
        elif method == "Adel05": # NOT YET IMPLEMENTED NON-FLAT N(z) YET
            
            def full_integrand(u):
                if flat_N_z == True:
                    N_z = 1 / self.photo_z["z_width"]
                else:
                    N_z = (popt[0] + (popt[1] * u)) / N_z_dz_integral # unitless
                x_z = cf.cosmo.comovingDistance(0, u) * cf.u.Mpc / cf.u.littleh # comoving distance
                H_z = cf.cosmo.Ez(u) * 100 * cf.u.littleh * cf.u.km / (cf.u.s * cf.u.Mpc) # Hubble constant
                integrand = (H_z * (x_z ** (- self.delta)) * (N_z ** 2) / cf.const.c).to((cf.u.littleh / cf.u.Mpc) ** (1 + self.delta), \
                                                cf.u.with_H0(cf.cosmo.H0 * cf.u.km / (cf.u.s * cf.u.Mpc)))
                return integrand.value
            full_integral = cf.scipy.integrate.quad(full_integrand, self.photo_z["min_z"], self.photo_z["max_z"])[0] \
                * ((cf.u.littleh / cf.u.Mpc) ** (1 + self.delta))
            r_0 = (((self.A * (cf.u.deg ** self.delta)).to(cf.u.rad ** self.delta)).value \
                   / (H_gamma * full_integral)) ** (1 / (1 + self.delta))
            r_0_err = (r_0 / (1 + self.delta)) * self.A_err / self.A
            
        elif method == "Lee06" or method == "Ouchi04":
            
            def full_integrand(u):
                if flat_N_z == True:
                    N_z = 1 / self.photo_z["z_width"]
                else:
                    N_z = (popt[0] + (popt[1] * u)) / N_z_dz_integral # unitless
                if method == "Ouchi04":
                    epsilon = -1.2
                    F_z = (1 + u) / ((1 + self.photo_z["mean_z"]) ** (-(3 + epsilon)))
                else:
                    F_z = 1
                d_A = cf.cosmo.angularDiameterDistance(u) * cf.u.Mpc / cf.u.littleh # angular diameter distance (astropy CHECKED)
                #print((d_A**(-delta)).unit)
                g_z = ((cf.cosmo.H0 * cf.u.km / (cf.u.s * cf.u.Mpc)) / cf.const.c) * ((1 + u) ** 2 * \
                    cf.np.sqrt(1 + cf.cosmo.Om0 * u + cf.cosmo.Ode0 * ((1 + u) ** (-2) - 1)))
                integrand = (cf.np.power(d_A, (- self.delta)) * g_z * F_z * (N_z ** 2)).to(1 / (cf.u.Mpc ** (1 + self.delta)), \
                                                cf.u.with_H0(cf.cosmo.H0 * cf.u.km / (cf.u.s * cf.u.Mpc)))
                return integrand.value
            full_integral, numerical_err = cf.scipy.integrate.quad(full_integrand, self.photo_z["min_z"], self.photo_z["max_z"])
            
            r_0 = (self.A / (H_gamma * full_integral)) ** (1 / (1 + self.delta)) * cf.u.Mpc # * (N_z_dz_integral ** 2)
            
            # error propagation
            
            # include Ouchi04 F_z here too
            if flat_N_z == False:
                def dI_by_dA_integrand(u):
                    N_z = (popt[0] + (popt[1] * u)) / N_z_dz_integral # unitless
                    d_A = cf.cosmo.angularDiameterDistance(u) * cf.u.Mpc / cf.u.littleh # angular diameter distance (astropy CHECKED)
                    g_z = ((cf.cosmo.H0 * cf.u.km / (cf.u.s * cf.u.Mpc)) / cf.const.c) * (1 + u) * \
                        cf.np.sqrt(cf.cosmo.Om0 * ((1 + u) ** 3 + (1 / cf.cosmo.Om0) - 1))
                    integrand = ((2 * N_z / N_z_dz_integral) * cf.np.power(d_A, (- self.delta)) * g_z * F_z).to(1 / (cf.u.Mpc ** (1 + self.delta)),\
                                                    cf.u.with_H0(cf.cosmo.H0 * cf.u.km / (cf.u.s * cf.u.Mpc)))
                    return integrand.value
                dI_by_dA, numerical_err = cf.scipy.integrate.quad(dI_by_dA_integrand, self.photo_z["min_z"], self.photo_z["max_z"])
                
                def dI_by_dB_integrand(u):
                    N_z = (popt[0] + (popt[1] * u)) / N_z_dz_integral # unitless
                    d_A = cf.cosmo.angularDiameterDistance(u) * cf.u.Mpc / cf.u.littleh # angular diameter distance (astropy CHECKED)
                    g_z = ((cf.cosmo.H0 * cf.u.km / (cf.u.s * cf.u.Mpc)) / cf.const.c) * (1 + u) * \
                        cf.np.sqrt(cf.cosmo.Om0 * ((1 + u) ** 3 + (1 / cf.cosmo.Om0) - 1))
                    integrand = ((2 * u * N_z / N_z_dz_integral) * cf.np.power(d_A, (- self.delta)) * g_z * F_z).to(1 / (cf.u.Mpc ** (1 + self.delta)),\
                                                    cf.u.with_H0(cf.cosmo.H0 * cf.u.km / (cf.u.s * cf.u.Mpc)))
                    return integrand.value
                dI_by_dB, numerical_err = cf.scipy.integrate.quad(dI_by_dB_integrand, self.photo_z["min_z"], self.photo_z["max_z"])
                
                full_integral_err = cf.np.sqrt((dI_by_dA * cf.np.sqrt(pcov[0][0]) / N_z_dz_integral) ** 2 + (dI_by_dB * cf.np.sqrt(pcov[1][1]) / N_z_dz_integral) ** 2)
            
            else:
                full_integral_err = 0  
            
            #print("A_err = {}, N_z_dz_err = {}, full_integral_err = {}".format(A_err * 100 / A, \
            #                            N_z_dz_err * 100 / N_z_dz_integral, full_integral_err * 100 / full_integral))
            #print("N_z_dz_integral = {}".format(N_z_dz_integral))
            r_0_err = (r_0 / (1 + self.delta)) * cf.np.sqrt((self.A_err / self.A) ** 2 + (full_integral_err / full_integral) ** 2) #+ (2 * N_z_dz_err / N_z_dz_integral) ** 2 \
                
        self.r_0  = r_0.to((cf.u.Mpc / cf.u.littleh), cf.u.with_H0(cf.cosmo.H0 * cf.u.km / (cf.u.s * cf.u.Mpc))) # unit conversion
        self.r_0_err = r_0_err.to((cf.u.Mpc / cf.u.littleh), cf.u.with_H0(cf.cosmo.H0 * cf.u.km / (cf.u.s * cf.u.Mpc))) # unit conversion

    def calculate_sigma_8_gal(self, method = "Adel05", flat_N_z = False, Nbins = 10, calc_r_0 = True):
        
        # if r_0 == None and delta == None:
        #     delta = self.delta
        #     r_0 = self.calculate_r0(gal_bias, method = r_0_method, flat_N_z = flat_N_z)
        # elif r_0 != None and delta != None:
        #     pass
        # else:
        #     raise ValueError("Must specify either all or neither of (r_0, delta)!")
        
        if calc_r_0 == True:
            self.calculate_r_0(method, flat_N_z, Nbins)
        
        numerator = 72 * ((self.r_0 / (8 * cf.u.Mpc / cf.u.littleh)) ** (1 + self.delta))
        denominator = (3 - (1 + self.delta)) * (4 - (1 + self.delta)) * (6 - (1 + self.delta)) * (2 ** (1 + self.delta))
        
        self.sigma_8_gal = cf.np.sqrt(numerator / denominator) # this square root makes minimal difference
        self.sigma_8_gal_err = (1 + self.delta) * self.sigma_8_gal * self.r_0_err / (2 * self.r_0)
    
    def calculate_b(self, method = "Adel05", flat_N_z = False, Nbins = 10, calc_r_0 = True, calc_s8 = True):
        
        # if r_0 == None and delta == None and z == None:
        #     z = self.photo_z["mid_z"]
        #     sigma8_gal = self.calculate_sigma8_gal(fitting_params, gal_bias, r_0_method = r_0_method)[0]
        #     sigma8_gal_err = self.calculate_sigma8_gal(fitting_params, gal_bias, r_0_method = r_0_method)[1]
        # elif r_0 != None and delta != None and z != None:
        #     sigma8_gal_params = self.calculate_sigma8_gal(fitting_params, gal_bias, r_0 = r_0, delta = delta)
        #     sigma8_gal = sigma8_gal_params[0]
        #     sigma8_gal_err = sigma8_gal_params[1]
        # else:
        #     raise ValueError("Must specify either all or neither of (r_0, delta, z)!")
        
        # calculate linear growth factor and it's normalization
        # def integrand(u):
        #     E_u = cf.cosmo.H((1 / u) - 1) / ((cf.cosmo.H0 * cf.u.km / (cf.u.s * cf.u.Mpc)) * cf.u.km / (cf.u.s * cf.u.Mpc))
        #     E_z = cf.cosmo.H(z) / (cf.cosmo.H0 * cf.u.km / (cf.u.s * cf.u.Mpc))
        #     return 2.5 * cf.cosmo.Om0 * E_z * cf.np.power(u * E_u, -3)
        # D_z, err = cf.scipy.integrate.quad(integrand, 0., 1 / (1 + z))
        # D_0, err = cf.scipy.integrate.quad(integrand, 0., 1.) # D(a=1)=1
        
        if calc_s8 == True:
            self.calculate_sigma_8_gal(method, flat_N_z, Nbins, calc_r_0)
        
        # correctly normalized
        sigma_8_m_z = cf.cosmo.sigma8 * cf.cosmo.growthFactor(self.photo_z["mean_z"])
        
        self.b = self.sigma_8_gal / sigma_8_m_z
        self.b_err = self.sigma_8_gal_err / sigma_8_m_z
        #print("z = %1.1f: b = %2.2f \u00B1 %2.2f" % (z, b, b_err))
        
    @staticmethod
    def create_random_gals(geometry, Ndata_gals, save_gals = False):
        # time this
        start_time_randoms = cf.time.time()
        
        Nrandom = Ndata_gals
        random_ra = []
        random_dec = []
        # generate random ra, dec over masked survey area
        while(Nrandom != 0):
            # choose random ra, dec over entire survey
            cf.np.random.seed() # set random number seed
            random_ra_loc = list(cf.np.random.uniform(geometry.unmasked_geometry["min_ra"].value, \
                                geometry.unmasked_geometry["max_ra"].value, size = Nrandom))
            cf.np.random.seed() # set different random number seed
            random_dec_loc = list(cf.np.random.uniform(geometry.unmasked_geometry["min_dec"].value, \
                                geometry.unmasked_geometry["max_dec"].value, size = Nrandom))
            
            # remove (ra, dec) co-ordinates that are within the masked region
            ra_index = cf.np.digitize(random_ra_loc, geometry.mask_bin_edges["ra"])
            dec_index = cf.np.digitize(random_dec_loc, geometry.mask_bin_edges["dec"])
            masked_ra_dec_indices = sorted([i for i in range(len(ra_index)) \
                             if geometry.mask[dec_index[i]][ra_index[i]] == 1], reverse = True)
            for index in masked_ra_dec_indices:
                random_ra_loc.pop(index)
                random_dec_loc.pop(index)
                
            # add to global ra, dec arrays
            random_ra = cf.np.append(random_ra, random_ra_loc)
            random_dec = cf.np.append(random_dec, random_dec_loc)
            
            # repeat until appropriate number of galaxies has been obtained
            Nrandom = Nrandom - len(random_ra_loc)
            print("Nrandom = " + str(Ndata_gals - Nrandom))
        
        # flatten random ra and random dec arrays
        random_ra = cf.np.reshape(random_ra, Ndata_gals)
        random_dec = cf.np.reshape(random_dec, Ndata_gals)
        
        if save_gals == True:
            # write random galaxies to csv, and append to the same file if it already exists
            #cf.os.chdir("/Users/user/Documents/PGR/UDS field")
            #filename = "UDS random gals ra dec"
            filename = cf.os.environ["RANDOM_GALS_FILE"]
            new_data = cf.pd.DataFrame({"ra": random_ra, "dec": random_dec})
            try:
                orig_data = cf.pd.read_csv(filename + ".csv")
                write_data = cf.pd.concat([orig_data, new_data])
            except FileNotFoundError:
                write_data = new_data
            print(write_data.shape)
            write_data.to_csv(filename + ".csv", mode = "w", index = False)
            print("Written to \'" + filename + "\'")
            
        end_time_randoms = cf.time.time()
        print("{} random galaxies took {} seconds to simulate!".format(str(Ndata_gals), \
                                cf.np.round(end_time_randoms - start_time_randoms, 2)))
        
        return random_ra, random_dec
    
    # HOD MCMC fitting with emcee and halomod -------------------------------------------
    # currently only works for Zehavi05 HOD model
    def calc_w_theta_hod(self, params):
        hm = self.fiducial_hm
        hm.update(hod_params = params)
        #pk = cf.scipy.interpolate.interp1d(hm.k_hm, hm.power_auto_tracer)
        return hm.angular_corr_gal
    @staticmethod
    def lnprior_flat(params):
        if params["M_min"] > 9.00 and params["M_min"] < 18.00 and params["M_1"] < 18.00 and params["M_1"] > 7.00 \
        and params["alpha"] < 3 and params["alpha"] > 0:
            return 0.0
        else:
            return -cf.np.inf
    def log_likelihood_flat(self, params, w_theta, w_theta_err):
        params_loc = dict({"M_min": params[0], "M_1": params[1], "alpha": params[2]})
        #print(params_loc)
        lp = self.lnprior_flat(params_loc)
        if not cf.np.isfinite(lp):
            return -cf.np.inf
        w_theta_model = self.calc_w_theta_hod(params_loc)
        return lp - 0.5 * cf.np.sum(((w_theta - w_theta_model) / w_theta_err) ** 2)
    
    def constrain_hod_params(self, ndim = 3, nwalkers = 50, nsteps = 100, HOD_model = "Zehavi05",\
                             backend_filename = "HOD_test_4"):
        
        # make dicts of parameters and their fiducial values for each of the different models

        # test against fiducial parameters
        fid_pars = [11.6222, 12.851, 1.049] # Zehavi05 model
        #print(log_likelihood_flat(fid_pars, hm_fid.angular_corr_gal, self.w_theta_err))
        
        # define initial positions
        pos = [fid_pars + 1e-2 * cf.np.random.uniform(0, 1, ndim) * fid_pars for i in range(nwalkers)]

        with cf.Pool() as pool:
            cf.os.chdir("/Users/user/Documents/PGR/UDS field/MCMC backends/")
            backend = cf.emcee.backends.HDFBackend(backend_filename + ".h5")
            sampler = cf.emcee.EnsembleSampler(nwalkers, ndim, self.log_likelihood_flat, \
                    args = (self.w_theta, self.w_theta_err), backend = backend, pool = pool)
            sampler.run_mcmc(pos, nsteps, progress = True)
            print("Final size: {0}".format(backend.iteration))  
            
    @staticmethod
    def plot_hod_params_corner(backend_filename, print_autocorr_time = False):
        
        fid_pars = [11.6222, 12.851, 1.049]
        cf.os.chdir("/Users/user/Documents/PGR/UDS field/MCMC backends/")
        reader = cf.emcee.backends.HDFBackend(backend_filename + ".h5")
        print(reader.get_last_sample())
        if print_autocorr_time == True:
            tau = reader.get_autocorr_time()
            print("autocorr time = " + str(tau))
        flat_samples = reader.get_chain(flat = True) # discard=100, thin=15,
        
        labels = ["M_min", "M_1", "alpha"]
        fig = cf.corner.corner(flat_samples, labels = labels, truths = fid_pars, \
                quantiles = [0.16, 0.5, 0.84], show_titles = True, title_kwargs = {"fontsize": 12}) 
    def add_to_backend(self, backend_filename, nsamples = 100, ndim = 3, nwalkers = 50):
        
        with cf.Pool() as pool:
            cf.os.chdir("/Users/user/Documents/PGR/UDS field/MCMC backends/")
            new_backend = cf.emcee.backends.HDFBackend(backend_filename + ".h5")
            print("Initial size: {0}".format(new_backend.iteration))
            new_sampler = cf.emcee.EnsembleSampler(nwalkers, ndim, self.log_likelihood_flat, \
                                    args = (self.w_theta, self.w_theta_err), backend = new_backend, pool = pool)
            new_sampler.run_mcmc(None, nsamples, progress = True)
            print("Final size: {0}".format(new_backend.iteration))
            
    # create an array of objects ---------------------------------------------------------
    @staticmethod
    def w_theta_obj_array(redshifts, sample_type, corr_type = "ACF", stellar_masses = None, fieldname = "UDS", w_theta_method = "halotools"):
        
        if corr_type not in ['ACF', 'CCF']:
            raise ValueError("corr_type must be 'ACF' or 'CCF'")
        
        sample_type_enum = cf.galaxy_sample_type[sample_type]
        w_theta_obj_array = []
        for z in redshifts:
            if sample_type_enum.name == "ALL":
                if stellar_masses != None:
                    raise TypeError("Stellar masses should not be given for ALL galaxies!")
                filename = "z = %1.0f, ALL Mstar %s w(theta) " % (z, fieldname) + corr_type
            elif sample_type_enum.name == "LOWER_STELLAR_MASS_LIM":
                try:
                    filename = "z = %1.0f, log(Mstar) > %1.1f, %s w(theta) " \
                              % (z, stellar_masses, fieldname) + corr_type
                except TypeError:
                    print("Stellar masses should be single valued (min. stellar mass)!")
            elif sample_type_enum.name == "STELLAR_MASS_BIN":
                try:
                    filename = "z = %1.0f, %1.1f < log(Mstar) < %1.1f, %s w(theta) " \
                        % (z, stellar_masses[0], stellar_masses[1], fieldname) + corr_type + ", " + w_theta_method
                except TypeError:
                    print("Stellar masses should have length of 2 (min. and max. stellar mass)!")
            w_theta_obj_loc = two_pt_angular_corr.from_pkl(filename)
            w_theta_obj_array = cf.np.append(w_theta_obj_array, w_theta_obj_loc)
            
        return w_theta_obj_array
    
# -----------------------------------------------------------------------------------------------------------

class two_pt_angular_auto_corr(two_pt_angular_corr):
    
    def __init__(self, fieldname, theta_mid_bins, w_theta, w_theta_err, bootstraps, \
                 galaxy_redshifts, galaxy_stellar_masses, sample_type, delta = 0.8):
        super().__init__(fieldname, theta_mid_bins, w_theta, w_theta_err, bootstraps, \
                     galaxy_redshifts, galaxy_stellar_masses, sample_type, delta)
    
    @classmethod # computationally expensive alternative constructor
    def calc_two_pt_angular_corr(cls, fieldname, data, sample_type, method = "landy-szalay", \
                        min_bin = 1e-10, max_bin = 2e2, Nbins = 20, Nbootstraps = 100, Ngals = None, \
                            w_theta_method = "halotools"):
        
        # time this
        start_time = cf.time.time()
        
        sample_type_enum = cf.galaxy_sample_type[sample_type]
        geometry = field_geometry(fieldname)
        
        # potential stumbling block in code here
        if Ngals == None:
            Ngals = len(data)
        else: # take a random sample of galaxies of length n_gals
            cf.np.random.seed()
            row_numbers_to_keep = cf.np.random.choice(len(data), Ngals, replace = False)
            data = data[row_numbers_to_keep]
            #print(len(data))
        
        # create an array of bins evenly spaced in log-space
        theta_edge_bins = 10 ** cf.np.linspace(cf.np.log10(min_bin), cf.np.log10(max_bin), Nbins + 1)
        theta = 0.5 * (theta_edge_bins[1:] + theta_edge_bins[:-1])
        
        # calculate w(theta) from data using bootstrapping to generate errors
        if w_theta_method == "astroML":
            results = cf.astroML.bootstrap_two_point_angular(data["ALPHA_J2000"], data["DELTA_J2000"], \
                            theta_edge_bins, geometry, method = method, Nbootstraps = Nbootstraps)
            w_theta, w_theta_err, bootstraps, random_gals_ra_dec = results
            
        elif w_theta_method == "halotools":        
            
            # bootstrap samples
            w_theta_repeats = []
            for repeat in range(Nbootstraps):
                # time this
                start_time_repeat = cf.time.time()   
                
                print("bootstrap = " + str(repeat))
                # resample distribution
                bootstrap_sample_data = cf.sklearn.utils.resample(data, replace = True)
                bootstrap_sample_ra_dec_data = cf.np.asarray([bootstrap_sample_data["ALPHA_J2000"], \
                                                         bootstrap_sample_data["DELTA_J2000"]]).T
                
                # create random galaxies in the masked field outside of the masked region
                random_ra, random_dec = two_pt_angular_corr.create_random_gals(geometry, Ngals)
                random_ra_dec = cf.np.vstack((random_ra, random_dec)).T
                
                # calculate w(theta) from data using bootstrapping to generate errors
                w_theta_loc = cf.angular_tpcf(bootstrap_sample_ra_dec_data, theta_edge_bins, randoms = random_ra_dec, \
                                          num_threads = "max", estimator = method)
                w_theta_repeats = cf.np.append(w_theta_repeats, w_theta_loc)
                
                end_time_repeat = cf.time.time()
                print("Bootstrap took {} seconds!".format(cf.np.round(end_time_repeat - start_time_repeat, 2)))
            
            w_theta_repeats = (cf.np.reshape(w_theta_repeats, (Nbootstraps, Nbins))).T
            bootstraps = w_theta_repeats
            
            w_theta = []
            w_theta_err = []
            for i in range(Nbins):
                w_theta_bin = cf.np.mean(w_theta_repeats[i])
                w_theta_err_bin = cf.np.std(w_theta_repeats[i])
                w_theta = cf.np.append(w_theta, w_theta_bin)
                w_theta_err = cf.np.append(w_theta_err, w_theta_err_bin)
                
            random_gals_ra_dec = {"ra": random_ra, "dec": random_dec}

        elif w_theta_method == "manual":
            
            def galaxy_pairs(data, correlation_type, random_ra, random_dec):
                
                # catch the error if correlation type is incorrectly inputted
                if correlation_type == "data-data" or correlation_type == "data-random" or correlation_type == "random-random":
                    # define correlation fields
                    if correlation_type == "data-data":
                        field_1 = cf.SkyCoord(ra = data["ALPHA_J2000"] * cf.u.deg, dec = data["DELTA_J2000"] * cf.u.deg)
                        field_2 = field_1
                        total_pairs = Ngals * (Ngals - 1) / 2 # total data-data pairs for normalization
                    else:
                        field_1 = cf.SkyCoord(ra = random_ra * cf.u.deg, dec = random_dec * cf.u.deg)
                        if correlation_type == "data-random":
                            field_2 = cf.SkyCoord(ra = data["ALPHA_J2000"] * cf.u.deg, dec = data["DELTA_J2000"] * cf.u.deg)
                            total_pairs = Ngals * len(random_ra) # total data-random pairs for normalization
                        elif correlation_type == "random-random":
                            field_2 = cf.SkyCoord(ra = random_ra * cf.u.deg, dec = random_dec * cf.u.deg)
                            total_pairs = len(random_ra) * (len(random_ra) - 1) / 2 # total random-random pairs for normalization
                else:
                    print("Incorrect input of the variable 'correlation_type'")
                
                # calculate the number of pairs between the data and/or random fields in each theta bin
                sep2d = field_2.search_around_sky(field_1, 100 * cf.u.deg)[2]
                sep2d = [x for x in sep2d.value if x != 0] # remove all occurrences of 0
                # add total number of pairs of galaxies for each angular separation bin to 'number_of_pairs'
                Npairs, bin_edges = cf.np.histogram(sep2d, bins = theta_edge_bins)
                if correlation_type != "data-random":
                    Npairs = Npairs / 2 # avoid double counting
                    
                if correlation_type == "data-data":
                    return theta, Npairs, total_pairs
                else:
                    return theta, Npairs, total_pairs
            
            # save repeated data to then calculate bootstrapped mean and std from
            w_theta_repeats = []
            for repeat in range(Nbootstraps):
                # time this
                start_time_repeat = cf.time.time()  
                
                # resample distribution
                bootstrap_sample = cf.sklearn.utils.resample(data, replace = True) # change sample size here
            
                # create random galaxies in the masked field outside of the masked region
                random_ra, random_dec = two_pt_angular_corr.create_random_gals(geometry, Ngals)
            
                theta_bins, data_data, total_data_data = \
                    galaxy_pairs(bootstrap_sample, "data-data", random_ra, random_dec)
                # # seed for random field
                # seed = cf.np.random.randint(0, 2 ** 32, 2, dtype = cf.np.uint32)
                # print(seed)
                print("finished DD")
                theta_bins, data_random, total_data_random = \
                    galaxy_pairs(bootstrap_sample, "data-random", random_ra, random_dec)
                print("finished DR")
                theta_bins, random_random, total_random_random = \
                    galaxy_pairs(bootstrap_sample, "random-random", random_ra, random_dec)
                print("finished RR")
                # normalize the data, random counts
                data_data_normed = data_data / total_data_data
                data_random_normed = data_random / total_data_random
                random_random_normed = random_random / total_random_random
                w_theta_loc = (data_data_normed - 2 * data_random_normed + random_random_normed) / random_random_normed
                w_theta_repeats = cf.np.append(w_theta_repeats, w_theta_loc)
                
                end_time_repeat = cf.time.time()
                print("Bootstrap took {} seconds!".format(cf.np.round(end_time_repeat - start_time_repeat, 2)))
            
            w_theta_repeats = (cf.np.reshape(w_theta_repeats, (Nbootstraps, Nbins))).T
            bootstraps = w_theta_repeats
            
            w_theta = []
            w_theta_err = []
            for i in range(Nbins):
                w_theta_bin = cf.np.mean(w_theta_repeats[i])
                w_theta_err_bin = cf.np.std(w_theta_repeats[i])
                w_theta = cf.np.append(w_theta, w_theta_bin)
                w_theta_err = cf.np.append(w_theta_err, w_theta_err_bin)
                
            random_gals_ra_dec = {"ra": random_ra, "dec": random_dec}
            
        else:
            raise SystemError("Invalid 'w_theta_method' input!")
            
        # plot the field containing data and random galaxies
        geometry.plot_field({"ra": data["ALPHA_J2000"], "dec": data["DELTA_J2000"]}, random_gals_ra_dec)
        
        end_time = cf.time.time()
        print("This took {} seconds!".format(cf.np.round(end_time - start_time, 2)))
        
        # calculate mid-photo-z and mid-stellar mass for naming
        min_z = cf.np.min(data["z_p"])
        max_z = cf.np.max(data["z_p"])
        mid_z = (max_z - min_z) / 2 + min_z
        min_stellar_mass = cf.np.log10(cf.np.min(data["Mstar_z_p"]))
        max_stellar_mass = cf.np.log10(cf.np.max(data["Mstar_z_p"]))
        
        # save data as pickle
        #cf.os.chdir("/Users/user/Documents/PGR/UDS field/pkl")
        # include sample type here
        if sample_type_enum.name == "ALL":
            save_name = "z = %1.0f, ALL Mstar %s w(theta) ACF" \
                      % (mid_z, fieldname)
        elif sample_type_enum.name == "LOWER_STELLAR_MASS_LIM":
            save_name = "z = %1.0f, log(Mstar) > %1.1f, %s w(theta) ACF" \
                      % (mid_z, min_stellar_mass, fieldname)
        elif sample_type_enum.name == "STELLAR_MASS_BIN":
            save_name = "z = %1.0f, %1.1f < log(Mstar) < %1.1f, %s w(theta) ACF, " \
                % (mid_z, min_stellar_mass, max_stellar_mass, fieldname) + w_theta_method
                      
        with open(save_name + ".pkl", 'wb') as f:
            cf.pickle.dump([fieldname, theta, w_theta, w_theta_err, bootstraps, \
                            data["z_p"], data["Mstar_z_p"], sample_type_enum, random_gals_ra_dec], f)
        return cls(fieldname, theta, w_theta, w_theta_err, bootstraps, data["z_p"], \
                   data["Mstar_z_p"], sample_type_enum, random_gals_ra_dec)     
    
    
class two_pt_angular_cross_corr(two_pt_angular_corr):
    
    def __init__(self, fieldname, theta_mid_bins, w_theta, w_theta_err, bootstraps, \
                 galaxy_redshifts, galaxy_stellar_masses, sample_type, delta = 0.6):
        super().__init__(fieldname, theta_mid_bins, w_theta, w_theta_err, bootstraps, \
                     galaxy_redshifts, galaxy_stellar_masses, sample_type, delta)
    
    @classmethod # computationally expensive alternative constructor
    def calc_two_pt_angular_corr(cls, fieldname, data, sample_type, cross_data, \
                                 method = "Landy-Szalay", min_bin = 1e-3, max_bin = 1e-1, \
                                     Nbins = 10, Nbootstraps = 100, Ncrossgals = None):
        
        # time this
        start_time = cf.time.time()
        
        sample_type_enum = cf.galaxy_sample_type[sample_type]
        
        if Ncrossgals == None:
            Ncrossgals = len(cross_data)
            
        # create an array of bins evenly spaced in log-space
        theta_edge_bins = 10 ** cf.np.linspace(cf.np.log10(min_bin), cf.np.log10(max_bin), Nbins + 1)
        theta = 0.5 * (theta_edge_bins[1:] + theta_edge_bins[:-1])    
        
        w_theta_repeats = []
        for repeat in range(Nbootstraps):
            # resample distributions
            bootstrap_sample_data = cf.sklearn.utils.resample(data, replace = True)
            bootstrap_sample_ra_dec_data = cf.np.asarray([bootstrap_sample_data["ALPHA_J2000"], \
                                                     bootstrap_sample_data["DELTA_J2000"]]).T
            bootstrap_sample_data_2 = cf.sklearn.utils.resample(cross_data, replace = True, n_samples = Ncrossgals)
            bootstrap_sample_ra_dec_data_2 = cf.np.asarray([bootstrap_sample_data_2["ALPHA_J2000"], \
                                                     bootstrap_sample_data_2["DELTA_J2000"]]).T
            # calculate w(theta) from data using bootstrapping to generate errors
            w_theta_loc = cf.angular_tpcf(bootstrap_sample_ra_dec_data, theta_edge_bins, \
                        sample2 = bootstrap_sample_ra_dec_data_2, do_auto = False, estimator = method)
            w_theta_repeats = cf.np.append(w_theta_repeats, w_theta_loc)
        
        w_theta_repeats = (cf.np.reshape(w_theta_repeats, (Nbootstraps, Nbins))).T
        bootstraps = w_theta_repeats
        print(w_theta_repeats)
        
        w_theta = []
        w_theta_err = []
        for i in range(Nbins):
            
            w_theta_bin = cf.np.mean(w_theta_repeats[i])
            w_theta_err_bin = cf.np.std(w_theta_repeats[i])
            w_theta = cf.np.append(w_theta, w_theta_bin)
            w_theta_err = cf.np.append(w_theta_err, w_theta_err_bin)
        
        end_time = cf.time.time()
        print("This took {} seconds!".format(cf.np.round(end_time - start_time, 2)))
        
        # calculate mid-photo-z and mid-stellar mass for naming
        min_z = cf.np.min(data["z_p"])
        max_z = cf.np.max(data["z_p"])
        mid_z = (max_z - min_z) / 2 + min_z
        min_stellar_mass = cf.np.log10(cf.np.min(data["Mstar_z_p"]))
        max_stellar_mass = cf.np.log10(cf.np.max(data["Mstar_z_p"]))
        
        # save data as pickle
        cf.os.chdir("/Users/user/Documents/PGR/UDS field/pkl")
        # include sample type here
        if sample_type_enum.name == "ALL":
            save_name = "z = %1.0f, ALL Mstar %s w(theta) ACF" \
                      % (mid_z, fieldname)
        elif sample_type_enum.name == "LOWER_STELLAR_MASS_LIM":
            save_name = "z = %1.0f, log(Mstar) > %1.1f, %s w(theta) CCF" \
                      % (mid_z, min_stellar_mass, fieldname)
        elif sample_type_enum.name == "STELLAR_MASS_BIN":
            save_name = "z = %1.0f, %1.1f < log(Mstar) < %1.1f, %s w(theta) CCF" \
                % (mid_z, min_stellar_mass, max_stellar_mass, fieldname)
                      
        with open(save_name + ".pkl", 'wb') as f:
            cf.pickle.dump([fieldname, theta, w_theta, w_theta_err, bootstraps, \
                            data["z_p"], data["Mstar_z_p"], sample_type_enum], f)
        return cls(fieldname, theta, w_theta, w_theta_err, bootstraps, data["z_p"], \
                   data["Mstar_z_p"], sample_type_enum)
    

# class two_pt_angular_dm_corr(two_pt_angular_corr):
#     def __init__(self, theta_mid_bins, w_theta, w_theta_err, bootstraps,\
#              Ngals, redshift_distribution, HOD, theta_edge_bins = None):
#         pass
    
# -----------------------------------------------------------------------------
    
def test_literature_paper(z_min, z_max, A, A_err, r_0 = None, r_0_err = None, s8 = None, s8_err = None, \
                          delta = 0.8, flat_N_z = True, method = "Wake11"):
    # create blank object
    w_theta_obj = two_pt_angular_auto_corr("UDS", [1e-3, 2e-1], [1.01, 1.01], [0.01, 0.01], 2, \
                                           [z_min, z_max], [1e10, 1e12], "STELLAR_MASS_BIN")
    
    # set relevant object parameters for calculating bias
    w_theta_obj.A = A
    w_theta_obj.A_err = A_err
    w_theta_obj.delta = delta
    w_theta_obj.C = 0
    w_theta_obj.C_err = 0
    w_theta_obj.A_fit_red_chi_sq = 0
    
    if r_0 != None and r_0_err != None:
        w_theta_obj.r_0 = r_0 * cf.u.Mpc / cf.u.littleh
        w_theta_obj.r_0_err = r_0_err * cf.u.Mpc / cf.u.littleh
        w_theta_obj.calculate_b(flat_N_z = flat_N_z, calc_r_0 = False, calc_s8 = True, method = method)
    elif s8 != None and s8_err != None:
        w_theta_obj.sigma_8_gal = s8
        w_theta_obj.sigma_8_gal_err = s8_err
        w_theta_obj.calculate_b(flat_N_z = flat_N_z, calc_r_0 = False, calc_s8 = False, method = method)
    else:
        w_theta_obj.calculate_b(flat_N_z = flat_N_z, calc_r_0 = True, calc_s8 = True, method = method)
    print(w_theta_obj.photo_z["mean_z"])
    w_theta_obj.print_model_params()


def main():
    redshifts = [1, 2, 3, 4]
    #gal_biases = [2.1, 2.7, 2.7, 3.0]
    w_theta_obj_array = two_pt_angular_corr.w_theta_obj_array(redshifts, "ALL", "ACF")
    for i in range(len(w_theta_obj_array)):
        #w_theta_obj_array[i].calculate_A(print_params = True)
        w_theta_obj_array[i].print_model_params(print_r0_s8gal_b = True)
        print(cf.np.log10(w_theta_obj_array[i].stellar_masses["mean_Mstar"]))
    #two_pt_angular_corr.plot_two_pt_angular_corr(w_theta_obj_array, "A")
    
if __name__ == "__main__":
    z_mid = int(cf.os.environ["Z_MID"])
    print(z_mid)
    #two_pt_angular_corr.create_random_gals(field_geometry("UDS"), Ngals, save_gals = True)
    #main()
    #test_literature_paper(0.9, 1.3, 70.55e-3, 15.07e-3, delta = 0.6, method = "Adel05")
    #print((5.7 * cf.u.Mpc).to(cf.u.Mpc / cf.u.littleh, cf.u.with_H0(cf.cosmo.H0 * cf.u.km / (cf.u.s * cf.u.Mpc))))
    

