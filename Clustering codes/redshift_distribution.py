#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 08:36:57 2021

@author: u92876da
"""

# redshift_distribution.py

import config as cf
from cosmic_var import cosmic_var

class redshift_distribution:
    
    def __init__(self, fieldname, galaxy_redshifts, Nbins, gal_bias):
        
        self.fieldname = fieldname
        self.galaxy_redshifts = galaxy_redshifts
        self.Nbins = Nbins
        self.gal_bias = gal_bias
        
        
        #self.P_z_distributions
    
    @property
    def z_edge_bins(self):
        return cf.np.linspace(cf.np.min(self.galaxy_redshifts), cf.np.max(self.galaxy_redshifts), self.Nbins + 1)
    
    @property
    def z_mid_bins(self):
        return 0.5 * (self.z_edge_bins[1:] + self.z_edge_bins[:-1])
    
    @property 
    def N_z_over_dz(self):
        return cf.np.histogram(self.galaxy_redshifts, bins = self.z_edge_bins,\
                                           weights = cf.np.full(len(self.galaxy_redshifts), 1 / self.bin_width))[0]
    @property
    def N_z(self):
        return self.N_z_over_dz * self.bin_width
       
    @property
    def N_z_over_dz_errs(self):
        return cf.np.sqrt(self.N_z + (self.N_z ** 2) * cf.np.full(len(self.N_z), \
            self.cosmic_variance.unmasked_gal_cosmic_var * (0.1 / self.bin_width))) \
            / self.bin_width
        
    @property
    def bin_width(self):
        return self.z_mid_bins[1] - self.z_mid_bins[0]
    
    @property
    def min_z(self):
        return cf.np.round(cf.np.min(self.z_mid_bins) - (self.bin_width / 2), 1)
    
    @property
    def max_z(self):
        return cf.np.round(cf.np.max(self.z_mid_bins) + (self.bin_width / 2), 1)
    
    @property
    def mid_z(self):
        return cf.np.round(self.min_z + (self.z_width / 2), 1)
    
    @property
    def z_std(self):
        return cf.np.std(self.galaxy_redshifts)
        
    @property
    def z_width(self):
        return self.max_z - self.min_z
    
    @property
    def Ngals(self):
        return len(self.galaxy_redshifts)
    
    @property
    def lin_fit_params(self):
        return self.calc_linear_fit()
    
    @property
    def cosmic_variance(self):
        return cosmic_var(self.fieldname, self.mid_z, self.gal_bias)
    
    @staticmethod
    def lin_fit(z, A, B):
        return A + B * z
    
    def calc_linear_fit(self):
        
        popt, pcov = cf.scipy.optimize.curve_fit(self.lin_fit, self.z_mid_bins, self.N_z_over_dz, \
                                        sigma = self.N_z_over_dz_errs, absolute_sigma = True)
        model = self.lin_fit(self.z_mid_bins, *popt)
        chi_sq = 0
        for i in range(len(self.z_mid_bins)):
            chi_sq += ((self.N_z_over_dz[i] - model[i]) / self.N_z_over_dz_errs[i]) ** 2
        red_chi_sq = chi_sq / (len(self.z_mid_bins) - len(popt))
        
        return model, popt, pcov, red_chi_sq
    
    # MCMC fitting ---------------------------
    @staticmethod
    def lnprior_flat(params):
        if params["N_bins"] >= 2 and params["N_bins"] <= 300 and params["A"] < 500000 and params["A"] > 0 \
        and params["B"] < 150000 and params["B"] > -300000:
            return 0.0
        else:
            return -cf.np.inf
    def log_likelihood_flat(self, params):
        params_loc = dict({"N_bins": int(params[0]), "A": params[1], "B": params[2]})
        #print(params_loc)
        lp = self.lnprior_flat(params_loc)
        if not cf.np.isfinite(lp):
            return -cf.np.inf
        self.Nbins = params_loc["N_bins"]
        #print(self.Nbins)
        lin_hist_model = self.lin_fit(self.z_mid_bins, params_loc["A"], params_loc["B"])
        return lp - 0.5 * cf.np.sum(((self.N_z_over_dz - lin_hist_model) / self.N_z_over_dz_errs) ** 2)
    
    def lin_fit_monte_carlo(self, ndim = 3, nwalkers = 50, nsteps = 500, backend_filename = "Hist_MC_test", fid_Nbins = 100):
        
        fid_pars = [fid_Nbins, 175000, -100000]
        pos = [fid_pars + 1e-2 * cf.np.random.uniform(0, 1, ndim) * fid_pars for i in range(nwalkers)]
        
        with cf.Pool() as pool:
            cf.os.chdir("/Users/user/Documents/PGR/UDS field/MCMC backends/")
            backend = cf.emcee.backends.HDFBackend(backend_filename + ".h5")
            sampler = cf.emcee.EnsembleSampler(nwalkers, ndim, self.log_likelihood_flat, backend = backend, pool = pool)
            sampler.run_mcmc(pos, nsteps, progress = True)
            print("Final size: {0}".format(backend.iteration))
    
    def plot_corner(self, backend_filename = "Hist_MC_test", print_autocorr_time = True, save_fig = True, fid_Nbins = -1000):
        #fid_pars = [11.6222, 12.851, 1.049]
        cf.os.chdir("/Users/user/Documents/PGR/UDS field/MCMC backends/")
        reader = cf.emcee.backends.HDFBackend(backend_filename + ".h5")
        fid_pars = [fid_Nbins, 175000, -100000]
        #print(reader.get_last_sample())
        if print_autocorr_time == True:
            tau = reader.get_autocorr_time()
            print("autocorr time = " + str(tau))
        flat_samples = reader.get_chain(flat = True) # discard=100, thin=15,
        
        labels = ["N_bins", "A", "B"]
        fig = cf.corner.corner(flat_samples, labels = labels, truths = fid_pars, \
                quantiles = [0.16, 0.5, 0.84], show_titles = True, title_kwargs = {"fontsize": 12})  # truths = fid_pars, \
        
        if save_fig == True:
            cf.os.chdir("/Users/user/Documents/PGR/UDS field/Plots/Redshift distributions/")
            cf.plt.savefig("z = %1.1f " %self.mid_z + "redshift dist. MCMC.jpeg")
            
    def add_to_backend(self, backend_filename, nsamples = 1000, ndim = 3, nwalkers = 50):
        
        with cf.Pool() as pool:
            cf.os.chdir("/Users/user/Documents/PGR/UDS field/MCMC backends/")
            new_backend = cf.emcee.backends.HDFBackend(backend_filename + ".h5")
            print("Initial size: {0}".format(new_backend.iteration))
            new_sampler = cf.emcee.EnsembleSampler(nwalkers, ndim, self.log_likelihood_flat, backend = new_backend, pool = pool)
            new_sampler.run_mcmc(None, nsamples, progress = True)
            print("Final size: {0}".format(new_backend.iteration))
            
    def plot_hist(self, fit_type = "linear", save_fig = False):
        
        # plot histogram
        # cf.plt.hist(self.galaxy_redshifts, bins = self.z_edge_bins,# stacked = True, \
        #                                    weights = cf.np.full(len(self.galaxy_redshifts), 1 / self.bin_width))
        cf.plt.bar(self.z_mid_bins, self.N_z_over_dz, width = self.bin_width, fill = True,\
                label = "N = %5.0f" % (self.Ngals), alpha = 0.5)
        # plot histogram data
        cf.plt.errorbar(self.z_mid_bins, self.N_z_over_dz, self.N_z_over_dz_errs,\
                        ls = "None", fmt = "s", c = "black", capsize = 5, capthick = 1)
            
        if fit_type == "linear":
            model, popt, pcov, red_chi_sq = self.lin_fit_params
            cf.plt.plot(self.z_mid_bins, model, c = "red")
        if fit_type != None:    
            self.print_fit_params(fit_type)

        # plot parameters
        plot_title = "z = %1.1f \u00B1 %1.1f redshift distribution" % (self.mid_z, self.z_width / 2)
        cf.plt.title(plot_title, fontsize = 14)
        cf.plt.xlabel("z", fontsize = 12)
        cf.plt.ylabel("N(z) / \u0394z", fontsize = 12)
        cf.plt.xticks(fontsize = 10)
        cf.plt.yticks(fontsize = 10)
        cf.plt.xlim(self.min_z, self.max_z)
        cf.plt.legend(loc="upper right", frameon=True, fontsize=10, bbox_to_anchor=(1, 1), fancybox=True, shadow=True)
        if save_fig == True:
            cf.plt.savefig(plot_title + ".jpeg", dpi=600)
        cf.plt.show() 
    
    def print_fit_params(self, fit_type):
        
        if fit_type == "linear":
            model, popt, pcov, red_chi_sq = self.lin_fit_params
            popt = cf.np.round(popt, 1)
            pcov = cf.np.round(pcov, 1)
            red_chi_sq = cf.np.round(red_chi_sq, 1)
            print("f(z) = A + Bz \n [A, B, \u03C7^2 red.] = [%1.1f \u00B1 %1.1f, %1.1f \u00B1 %1.1f, %1.1f]" \
                  % (popt[0], cf.np.sqrt(pcov[0][0]), popt[1], cf.np.sqrt(pcov[1][1]), red_chi_sq))

if __name__ == "__main__":
    z = 4.0
    z_width = 0.6
    
    HDU_list = cf.fits.open("/Users/user/Documents/PGR/UDS field/DR11-2arcsec-Jun-30-2019.fits")
    data = HDU_list[1].data
    data_indices = cf.np.linspace(0, len(data) - 1, len(data))
    indice_col = cf.fits.ColDefs([
            cf.fits.Column(name = "index", format='J', array = data_indices)])
    data = (cf.fits.BinTableHDU.from_columns(data.columns + indice_col)).data
    
    data = data[data["z_p"] >= z - z_width / 2]
    data = data[data["z_p"] <= z + z_width / 2]
    print(data["index"].shape)
    
    P_z_distributions = cf.np.load("/Users/user/Documents/PGR/UDS field/P_z distribution/pz_z10_32bit.npy")
    P_z = cf.np.array([P_z_distributions[index] for index in data["index"]])
    z_grid = cf.np.loadtxt("/Users/user/Documents/PGR/UDS field/P_z distribution/DR11pz_zgrid.txt")
    print(P_z.shape)
    # P_z_distributions = cf.pd.read_csv("/Users/user/Documents/PGR/UDS field/P_z distribution/DR11pz.txt", \
    #                                   delim_whitespace = True, skipfooter = 280000)
    cf.plt.plot(z_grid, P_z[3600])
    