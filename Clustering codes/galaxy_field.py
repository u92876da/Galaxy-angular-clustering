# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 21:56:28 2021

@author: dunca
"""

# galaxy_field.py

import config as cf
from two_pt_angular_correlation import two_pt_angular_corr
from two_pt_angular_correlation import two_pt_angular_auto_corr
from two_pt_angular_correlation import two_pt_angular_cross_corr
from redshift_distribution import redshift_distribution
from field_geometry import field_geometry
from bias import bias_evolution
from cosmic_var import cosmic_var
from photo_z_errs import photo_z_errs
from stellar_mass_function import stellar_mass_function as smf
import MCMC

class galaxy_field:
    
    def __init__(self, fieldname = "Rachana UDS DR11", z = None, z_width = None, \
        min_stellar_mass = None, max_stellar_mass = None, w_theta_ACF_pkl = None, \
        N_z_bins = 10):
        
        self.fieldname = fieldname
        #self.HDU_list = cf.fits.open(datafile) # CANNOT PICKLE THIS
        if cf.os.getenv("FILENAME", False):
            data = cf.fits.open(cf.os.environ["FILENAME"])
            data = cf.fits.open(cf.os.environ("FIELDNAME"))[1].data
        else:
            data = cf.fits.open(cf.filenames[fieldname])[1].data
        #data = cf.rewrite_matched_table_columns(data)
        data["Mstar_z_p"] = data["Mstar_m62_z_p"] # solar metallicity Z = 0.02
        data = data[data["Best galaxies"] == True]
        self.data = data
        print("data length =", len(self.data))
        
        if z != None and z_width != None:
            self.photo_z_cut(z, z_width)
            
        self.stellar_mass_cut(min_stellar_mass, max_stellar_mass)
        
        if w_theta_ACF_pkl != None:
            self.w_theta_ACF = two_pt_angular_auto_corr.from_pkl(w_theta_ACF_pkl)
            
        # self.z_distribution = redshift_distribution(self.fieldname, self.data["z_p"], \
        #                                             N_z_bins, gal_bias)
    
    @property
    def geometry(self):
        return field_geometry(self.fieldname)
    
    @property
    def survey_volume(self):
        def integrand(u):
            integrand = cf.astropy_cosmo.differential_comoving_volume(u).value # potentially shouldn't be comoving volume
            return integrand
        comoving_volume_element = cf.scipy.integrate.quad(integrand, \
                                    self.photo_z["min_z"], self.photo_z["max_z"])[0]
        comoving_volume = comoving_volume_element * (cf.u.Mpc ** 3 / cf.u.sr) \
            * self.geometry.masked_area.to(cf.u.sr) # in Mpc**3
        return comoving_volume.to((cf.u.Mpc / cf.u.littleh) ** 3, \
                cf.u.with_H0(cf.cosmo.H0 * cf.u.km / (cf.u.s * cf.u.Mpc)))
    
    @property
    def Ngals(self):
        return len(self.data['z_p'])
    
    def Ngals_err(self):
        # make this over the masked area
        gal_cosmic_var = self.cosmic_variance().unmasked_gal_cosmic_var(self.photo_z["z_width"])
        
        return cf.np.sqrt(self.Ngals + gal_cosmic_var * self.Ngals ** 2)
    
    @property
    def photo_z(self):
        min_z = cf.np.min(self.data['z_p'])
        max_z = cf.np.max(self.data['z_p'])
        return {"mean_z": cf.np.mean(self.data['z_p']), "std_z": \
                cf.np.std(self.data['z_p']), "min_z": min_z, "max_z": max_z, \
                "mid_z": cf.np.round((max_z - min_z) / 2 + min_z, 1), \
                    "z_width": max_z - min_z}        
            
    @property
    def photo_z_errs(self):
        return photo_z_errs(self.fieldname)
    
    @property
    def stellar_mass_function(self):
        return smf.calculate_smf(self.data["Mstar_z_p"], self.photo_z, self.geometry)
    
    def gal_bias(self):
        return [self.w_theta_ACF.b, self.w_theta_ACF.b_err]
        
    # def M_halo(self):
    #     return self.w_theta_ACF.calculate_M_h()
    
    def galaxy_density(self): # [arcmin^-2]
        density = (self.Ngals / (self.geometry.masked_area)).to(1 / (cf.u.arcmin ** 2))
        density_err = (self.Ngals_err() / (self.geometry.masked_area)).to(1 / (cf.u.arcmin ** 2))
        print("z = %1.1f: N = (%1.2f ± %1.2f) arcmin^-2" \
              % (self.photo_z["mid_z"], density.value, density_err.value))
        self.density = density
        self.density_err = density_err
        return density, density_err
    
    def galaxy_abundance(self): # [h^3 Mpc^-3]
        abundance = self.Ngals / self.survey_volume
        abundance_err = self.Ngals_err() / self.survey_volume
        print("z = %1.1f: n = (%1.1f ± %1.1f) x 10^-4 h^3 Mpc^-3" \
              % (self.photo_z["mid_z"], abundance.value * 1e4, abundance_err.value * 1e4))
        self.abundance = abundance
        self.abundance_err = abundance_err
        return abundance, abundance_err

    def cosmic_variance(self):
        # gal_bias = self.w_theta_ACF.b
        return cosmic_var(self.fieldname, self.photo_z["mid_z"], self.gal_bias()[0])
    
    def stellar_mass(self, metallicity = None):
        if metallicity == None:
            mean_mass = cf.np.mean(self.data["Mstar_z_p"])
            min_mass = cf.np.min(self.data["Mstar_z_p"])
            max_mass = cf.np.max(self.data["Mstar_z_p"])
        else:
            mean_mass = cf.np.mean(self.data["Mstar_m" + metallicity + "_z_p"])
            min_mass = cf.np.min(self.data["Mstar_m" + metallicity + "_z_p"])
            max_mass = cf.np.max(self.data["Mstar_m" + metallicity + "_z_p"])   
        return {"mean_stellar_mass": mean_mass, "min_stellar_mass": min_mass, \
                "max_stellar_mass": max_mass, "mid_stellar_mass": (max_mass - min_mass) \
                    / 2 + min_mass}
    
    def z_distribution(self, Nbins = 10):
        # gal_bias = self.w_theta_ACF.b
        return redshift_distribution(self.fieldname, self.data["z_p"], Nbins, self.gal_bias())
    
    def photo_z_cut(self, z, z_width):
        self.data = self.data[self.data["z_p"] >= z - z_width / 2]
        self.data = self.data[self.data["z_p"] <= z + z_width / 2]
    
    def photo_band_cut(self, band, min_mag, max_mag):
        pass
    
    def gal_colour_cut(self, method, colour): # BzK, UVJ, etc..
        pass
    
    def stellar_mass_cut(self, min_mass, max_mass):
        if min_mass != None:
            print("min_mass =", min_mass)
            self.data = self.data[self.data['Mstar_z_p'] >= 10 ** min_mass]
        if max_mass != None:
            print("max_mass =", max_mass)
            self.data = self.data[self.data['Mstar_z_p'] <= 10 ** max_mass]
    
    def plot_galaxies(self):
        data_gals = {"ra": self.data["ALPHA_J2000"], "dec": self.data["DELTA_J2000"]}
        random_gals = self.w_theta_ACF.random_gals_ra_dec
        self.geometry.plot_field(data_gals, random_gals)
    
    # computationally expensive
    def calc_w_theta_ACF(self, method = "landy-szalay", min_bin = 1e-3, max_bin = 2e-1, Nbins = 20, 
                         Nbootstraps = 100, Ngals = None, w_theta_method = "astroML"):
        
        self.w_theta_ACF = two_pt_angular_auto_corr.calc_two_pt_angular_corr(self.fieldname, \
            self.data, method, min_bin, max_bin, Nbins, Nbootstraps, Ngals, w_theta_method)
    
    def calc_w_theta_CCF(self, sample_type, cross_data, method = "landy-szalay", min_bin = 1e-3,\
                    max_bin = 2e-1, Nbins = 20, Nbootstraps = 100, Ngals = None):
        
        self.w_theta_CCF = two_pt_angular_cross_corr.calc_two_pt_angular_corr(self.fieldname, \
                    self.data, cross_data, method, min_bin, max_bin, \
                        Nbins, Nbootstraps, Ngals)
    
    def M_h_MCMC(self, fit_type, nwalkers, backend_filename, \
                 method = "1 parameter, Foucaud 2010"): # blobs here
        if method == "1 parameter, Foucaud 2010":
            if fit_type == "n":
                fit_data = [self.abundance]
                fit_data_err = [self.abundance_err]
            self.MCMC_HOD_fit = MCMC.Foucaud_2010_HOD_MCMC(fit_type, fit_data, fit_data_err, backend_filename, nwalkers)
    
    # @classmethod # may need to make this a class function
    # def from_pkl(pkl_name, fieldname = "UDS"):
    #     cf.os.chdir("/Users/user/Documents/PGR/UDS field/pkl")
    #     with open(pkl_name + '.pkl', 'rb') as f:
    #         gal_field_obj = cf.pickle.load(f)
    #         gal_field_obj.fieldname = fieldname
    #         #gal_field_obj.define_post_mask_data()
    #         return gal_field_obj
        
    @staticmethod # in two_point_angular_corr function
    def write_two_pt_params(gal_field_array, fitting_params, include_r0 = False, include_b = False):
        
        parameters = []
        row_labels = []
        for gal_field in gal_field_array:
            params_loc = gal_field.w_theta_obj.define_fitting_params(fitting_params)
            A, delta, C = list(zip(*params_loc[0]))[0]
            A_err, delta_err, C_err = list(zip(*params_loc[0]))[1]
            parameters = cf.np.append(parameters, (A, delta, C, A_err, delta_err, C_err, params_loc[1]))
            row_labels = cf.np.append(row_labels, "z={}".format(gal_field.mean_photo_z))
        
        parameters = cf.np.reshape(parameters, (len(gal_field_array), 7))
        column_labels = ["A", "delta", "C", "A_err", "delta_err", "C_err", "red_chi_sq"]
        
        # determine filename
        z_list = list([gal_field.mean_photo_z for gal_field in gal_field_array])
        delta = gal_field_array[0].w_theta_obj.define_fitting_params(fitting_params)[0][1][0]
        if len(z_list) == 1:
            filename = "z = {}, \u03B4 = {} w(theta) params".format(z_list[0], delta)
        else:
            filename = "z = {}-{}, \u03B4 = {} w(theta) params".format(cf.np.min(z_list),\
                                                                       cf.np.max(z_list), delta)
        # write data to csv
        cf.os.chdir("/Users/user/Documents/PGR/UDS field/Plots/2pt angular correlation/Fit parameters")
        write_data = cf.pd.DataFrame(data = parameters, columns = column_labels, index = row_labels)
        write_data.to_csv(filename + ".csv", mode='w')
        print("Written to \'" + filename + "\'")    
        
    # @staticmethod
    # def gal_field_obj_array(redshifts):
    #     gal_field_array = []
    #     for z in redshifts:
    #         gal_field_loc = galaxy_field.from_pkl("z={} UDS".format(z))
    #         gal_field_array = cf.np.append(gal_field_array, gal_field_loc)
    #     return gal_field_array
    
    @staticmethod
    def gal_field_obj_array(redshifts, z_width, sample_type, stellar_masses = None, \
            include_w_theta = True, corr_type = "ACF", fieldname = "UDS", w_theta_method = "astroML"):
        
        if corr_type not in ['None', 'ACF', 'CCF']:
            raise ValueError("corr_type must be 'None', 'ACF' or 'CCF'")
        
        sample_type_enum = cf.galaxy_sample_type[sample_type]
        gal_field_obj_array = []
        
        cf.os.chdir("/Users/user/Documents/PGR/UDS field/pkl")
        for i in range(len(redshifts)):

            if sample_type_enum.name == "ALL":
                if stellar_masses != None:
                    raise TypeError("Stellar masses should not be given for ALL galaxies!")
                filename = "z = %1.0f, All Mstar %s w(theta) " % (redshifts[i], fieldname) + corr_type
                if corr_type == "None":
                    gal_field_obj_loc = galaxy_field(fieldname, redshifts[i], z_width)
                elif corr_type == "ACF":
                    gal_field_obj_loc = galaxy_field(fieldname, redshifts[i], z_width, \
                                            w_theta_ACF_pkl = filename)
            elif sample_type_enum.name == "LOWER_STELLAR_MASS_LIM":
                try:
                    filename = "z = %1.0f, log(Mstar) > %1.1f, %s w(theta) " \
                              % (redshifts[i], stellar_masses, fieldname) + corr_type
                except TypeError:
                    print("Stellar masses should be single valued (min. stellar mass)!")
                if corr_type == "None":
                    gal_field_obj_loc = galaxy_field(fieldname, redshifts[i], z_width, stellar_masses, 99)
                elif corr_type == "ACF":
                    gal_field_obj_loc = galaxy_field(fieldname, redshifts[i], z_width, stellar_masses, 99, \
                                            w_theta_ACF_pkl = filename)
            elif sample_type_enum.name == "STELLAR_MASS_BIN":
                try:
                    filename = "z = %1.0f, %1.1f < log(Mstar) < %1.1f, %s w(theta) " \
                        % (redshifts[i], stellar_masses[0], stellar_masses[1], fieldname) + corr_type + ", " + w_theta_method
                except TypeError:
                    print("Stellar masses should have length of 2 (min. and max. stellar mass)!")
                if corr_type == "None":
                    gal_field_obj_loc = galaxy_field(fieldname, redshifts[i], z_width, stellar_masses[0], stellar_masses[1])
                elif corr_type == "ACF":
                    gal_field_obj_loc = galaxy_field(fieldname, redshifts[i], z_width, stellar_masses[0], \
                            stellar_masses[1], w_theta_ACF_pkl = filename)
            
            gal_field_obj_array = cf.np.append(gal_field_obj_array, gal_field_obj_loc)
        
        return gal_field_obj_array
    
# main -----------------------------------

def plot_w_theta(w_theta_methods, redshifts = [1, 2, 3, 4], z_width = 0.6, sample_type = "STELLAR_MASS_BIN", \
                 stellar_masses = [10.0, 10.5], corr_type = "ACF", save_fig = False):
    #UDS_field = galaxy_field(cf.UDS_filename) #, 1, 0.6, two_pt_angular_corr_pkl = "z=1 w_theta")
    
    gal_field_arr = []
    w_theta_obj_arr = []
    for i in range(len(w_theta_methods)):
        gal_field_arr_loc = galaxy_field.gal_field_obj_array(redshifts, z_width, sample_type, stellar_masses, \
                                                         w_theta_method = w_theta_methods[i])
        w_theta_obj_arr_loc = two_pt_angular_corr.w_theta_obj_array(redshifts, sample_type, corr_type, stellar_masses, \
                                                                w_theta_method = w_theta_methods[i])
        if i == 0:
            gal_field_arr = cf.np.append(gal_field_arr, gal_field_arr_loc)
            w_theta_obj_arr = cf.np.append(w_theta_obj_arr, w_theta_obj_arr_loc)
        else:
            gal_field_arr = cf.np.concatenate((gal_field_arr, gal_field_arr_loc))
            w_theta_obj_arr = cf.np.concatenate((w_theta_obj_arr, w_theta_obj_arr_loc))

    for i in range(len(gal_field_arr)):
        
        gal_field_arr[i].galaxy_density()
        gal_field_arr[i].galaxy_abundance()
        #gal_field_arr[i].plot_galaxies()
        #gal_field_array[i].z_distribution().plot_hist(save_fig = False)
        # plot w_theta here
        gal_field_arr[i].w_theta_ACF.print_model_params(print_r0_s8gal_b = True) # power law parameters
        #print(gal_field_arr[i].w_theta_ACF.A_err_cv_corrected)
        
    two_pt_angular_corr.plot_two_pt_angular_corr(w_theta_obj_arr, plot_literature = ["Foucaud_2010_z1"], save_fig = save_fig)


def plot_bias(redshifts = [1, 2, 3, 4], z_width = 0.6, save_fig = False):
    
    # Foucaud 2010 data
    Foucaud2010_10_10_5 = bias_evolution.from_literature(stellar_mass_min = 10, stellar_mass_max = 10.5)
    Foucaud2010_10_5_11 = bias_evolution.from_literature(stellar_mass_min = 10.5, stellar_mass_max = 11)
    Foucaud2010_11_12 = bias_evolution.from_literature(stellar_mass_min = 11, stellar_mass_max = 12)
    
    # This work
    # gal_field_arr_10 = galaxy_field.gal_field_obj_array(redshifts, z_width, "LOWER_STELLAR_MASS_LIM", 10)
    gal_field_arr_10_10_5 = galaxy_field.gal_field_obj_array(redshifts, z_width, "STELLAR_MASS_BIN", [10, 10.5])
    # gal_field_arr_10_5_11 = galaxy_field.gal_field_obj_array(redshifts, z_width, "STELLAR_MASS_BIN", [10.5, 11])
    # b_gal_10 = bias_evolution.from_gal_field_arr(gal_field_arr_10)
    b_gal_10_10_5 = bias_evolution.from_gal_field_arr(gal_field_arr_10_10_5)
    # b_gal_10_5_11 = bias_evolution.from_gal_field_arr(gal_field_arr_10_5_11)
    
    bias_evolution.plot_bias_evolution([b_gal_10_10_5, \
                                Foucaud2010_10_10_5, Foucaud2010_10_5_11, Foucaud2010_11_12], \
                                "z = 1.0-4.0 galaxy bias", save_fig = save_fig)

def plot_stellar_mass_functions():
    redshifts = [1, 2, 3, 4]
    gal_field_array = galaxy_field.gal_field_obj_array(redshifts, 0.6, \
                            "ALL", corr_type = "None")
    smf.plot_stellar_mass_functions(gal_field_array, "smf")

def calculate_w_theta(z, z_width, min_stellar_mass, max_stellar_mass, Nbootstraps = 100, w_theta_method = "astroML"):
    
    gal_field = galaxy_field("Rachana UDS DR11", z, z_width, min_stellar_mass, max_stellar_mass)
    #cross_field = galaxy_field("UDS", 4, 0.6, 10, 99)
    print("Ngals = %1.0f" %gal_field.Ngals)
    # if gal_field.Ngals > 10000:
    #     Ngals = 10000
    # else:
    #     Ngals = gal_field.Ngals
    Ngals = gal_field.Ngals
    #print("Ngals cross = %1.0f" %cross_field.Ngals)
    # gal_field.galaxy_density(4.9)
    # gal_field.galaxy_abundance(4.9)
    gal_field.calc_w_theta_ACF(Ngals = Ngals, Nbootstraps = Nbootstraps, w_theta_method = w_theta_method)
    # gal_field.galaxy_density()
    # gal_field.galaxy_abundance()
    # two_pt_angular_corr.plot_two_pt_angular_corr([gal_field.w_theta_ACF])
    # gal_field.w_theta_ACF.print_model_params(print_r0_s8gal_b = True)
    
def plot_completeness(save_fig = False):
    cf.os.chdir("/Users/user/Documents/PGR/UDS field/")
    comp_curve = cf.np.genfromtxt("K_band_source_completeness.txt")
    #print(comp_curve[:, 1])
    cf.plt.plot(comp_curve[:34, 0], 100 * comp_curve[:34, 1], c = "red", label = "UKIDSS UDS")
    cf.plt.title("UDS completeness curve", fontsize = 14)
    cf.plt.xlabel("K$_{AB}$", fontsize = 12)
    cf.plt.ylabel("Completeness / %", fontsize = 12)
    lgd = cf.plt.legend(loc="lower left", frameon=True, fontsize=10, bbox_to_anchor = (0, 0), fancybox=True, shadow=True)
    if save_fig == True:
        cf.os.chdir("/Users/user/Documents/PGR/UDS field/Plots")
        save_title = "UDS completeness curve.jpeg"
        cf.plt.savefig(save_title, dpi = 600, box_extra_artists = (lgd,), bbox_inches = "tight")
    cf.plt.show()
    
def test_M_h_from_Foucaud_2010(z_min, z_max, n, n_err, nwalkers, nsteps, \
                    backend_filename = "1 parameter HOD, Foucaud 2010, MCMC test 2"):
    gal_field = galaxy_field("UDS", z_min, z_max)
    gal_field.abundance = n
    gal_field.abundance_err = n_err
    gal_field.M_h_MCMC("n", nwalkers, backend_filename = backend_filename)
    gal_field.MCMC_HOD_fit.run_chain(nsteps)
    gal_field.MCMC_HOD_fit.plot_corner(True)

def from_shell():
    z_mid = float(cf.os.environ["Z_MID"])
    z_width = float(cf.os.environ["Z_WIDTH"])
    if cf.os.environ["STELLAR_MASS_MIN"] == "":
        min_stellar_mass = None
    else:
        min_stellar_mass = float(cf.os.environ["STELLAR_MASS_MIN"])
    if cf.os.environ["STELLAR_MASS_MAX"] == "":
        max_stellar_mass = None
    else:
        max_stellar_mass = float(cf.os.environ["STELLAR_MASS_MAX"])
    if bool(cf.os.environ["CALCULATE_ACF"]):
        calculate_w_theta(z_mid, z_width, min_stellar_mass, max_stellar_mass)

if __name__ == "__main__":
    if bool(cf.os.environ["FROM_SHELL"]):
        from_shell()
    #calculate_w_theta(4, 0.6, 11.0, 12.0, Nbootstraps = 10, w_theta_method = "astroML")
    #plot_w_theta(["astroML", "halotools", "manual"], save_fig = False)
    #plot_w_theta(["astroML"])
    # plot_completeness()
    #plot_bias()
    
    


    
        
