#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 13:19:13 2022

@author: u92876da
"""

import config as cf

class MCMC(cf.ABC):
    
    def __init__(self, data, data_errs, param_keys, fid_pars, flat_priors, backend_filename, blob_keys, nwalkers):
        
        if not len(data) == len(data_errs):
            raise IndexError("Length of data, data errors and model functions should be the same")
        
        if not len(param_keys) == len(fid_pars) == len(flat_priors):
            raise IndexError("Length of parameter keys and fiducial and prior parameters should be the same")
        
        self.data: list = data
        self.data_errs: list = data_errs
        self.param_keys: list = param_keys
        self.fid_pars: list = fid_pars
        self.flat_priors: list = flat_priors
        self.blob_keys: list = blob_keys
        self.backend_filename: str = backend_filename
        self.nwalkers: int = nwalkers
        
    @property
    def ndim(self):
        return len(self.param_keys)
    
    @cf.abstractmethod
    def model_funcs(self):
        pass
    
    @cf.abstractmethod
    def blob_funcs(self):
        pass
    
    
    def lnprior_flat(self, params):
        
        for i in range(len(self.param_keys)):
            if params[self.param_keys[i]] > self.flat_priors[i][0] and \
               params[self.param_keys[i]] < self.flat_priors[i][1]:
                pass
            else:
                return -cf.np.inf
        return 0.0
    
    
    def log_likelihood_flat(self, params, data, data_errs):
        
        params_loc = {}
        for i in range(len(self.param_keys)):
            params_loc[self.param_keys[i]] = params[i]
        #print(params_loc)
        lp = self.lnprior_flat(params_loc)
        
        if not cf.np.isfinite(lp):
            return cf.np.full(1 + len(self.blob_keys), -cf.np.inf)

        return_value = lp
        for i in range(len(self.data)):
            return_value = return_value - 0.5 * cf.np.sum(((data[i] - \
                                self.model_funcs(params_loc)[i]) / data_errs[i]) ** 2)
       
        if len(self.blob_keys) == 0:
            return return_value
        else:
            return return_value, self.blob_funcs(params_loc)
    
    
    def run_chain(self, nsteps):
        
        # change fid_pars into a list
        pos = [self.fid_pars + 1e-2 * cf.np.random.uniform(0, 1, self.ndim) * \
               self.fid_pars for i in range(self.nwalkers)]
        
        with cf.Pool() as pool:
            cf.os.chdir("/Users/user/Documents/PGR/UDS field/MCMC backends/")
            backend = cf.emcee.backends.HDFBackend(self.backend_filename + ".h5")
            blobs_dtype = [(key, float) for key in self.blob_keys]
            sampler = cf.emcee.EnsembleSampler(self.nwalkers, self.ndim, self.log_likelihood_flat, \
                    args = (self.data, self.data_errs), blobs_dtype = blobs_dtype, \
                        pool = pool, backend = backend)
            sampler.run_mcmc(pos, nsteps, progress = True)
            print("Final size: {0}".format(backend.iteration))
        pool.close()
        
            
    def plot_corner(self, print_autocorr_time = True):
        
        cf.os.chdir("/Users/user/Documents/PGR/UDS field/MCMC backends/")
        reader = cf.emcee.backends.HDFBackend(self.backend_filename + ".h5")
        # print(reader.get_last_sample())
        if print_autocorr_time == True:
            tau = reader.get_autocorr_time()
            print("autocorr time = " + str(tau))
        flat_samples = reader.get_chain(flat = True) # discard=100, thin=15,
        flat_blobs = reader.get_blobs(flat = True).reshape([self.nwalkers * reader.iteration, self.ndim])
        
        corner_data = cf.np.hstack([flat_samples, flat_blobs[self.blob_keys[0]]])
        corner_labels = cf.np.append([self.param_keys[i] for i in range(len(self.param_keys))], \
                                     [self.blob_keys[j] for j in range(len(self.blob_keys))])
        for i in range(len(self.blob_keys) - 1): # this part is untested
            corner_data = cf.np.hstack([corner_data, flat_blobs[self.blob_keys[i]]])
            
        fig = cf.corner.corner(corner_data, labels = corner_labels, quantiles = [0.16, 0.5, 0.84], \
                                   show_titles = True, title_kwargs = {"fontsize": 12}) 
          
# -----------------------------------------------------------------------------        
        
class Foucaud_2010_HOD_MCMC(MCMC):
    
    def __init__(self, fit_type, data, data_errs, backend_filename, nwalkers, \
                 hmf_model = "ST", bias_model = "SMT01"):
        self.hm = cf.TracerHaloModel(hmf_model = hmf_model, bias_model = bias_model, \
                            cosmo_model = cf.astropy_cosmo, hod_model = "Zehavi05")
        self.fit_type = fit_type
        super().__init__(data, data_errs, ["M_min"], [12.], [[8., 16.]], backend_filename, ["M_eff"], nwalkers)
      
    def model_funcs(self, params):
        if self.fit_type == "n":
            self.hm.update(hod_params = {"M_min": params["M_min"], "M_1": 20 * params["M_min"], "alpha": 1})
            return [self.hm.mean_tracer_den]
            
    def blob_funcs(self, params):
        if "M_eff" in self.blob_keys:
            self.hm.update(hod_params = {"M_min": params["M_min"], "M_1": 20 * params["M_min"], "alpha": 1})
            return self.hm.mass_effective
    
    