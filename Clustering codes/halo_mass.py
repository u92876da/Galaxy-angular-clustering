#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 22 15:26:33 2022

@author: u92876da
"""

import config as cf

def calc_halo_bias(params, z):
    M = (10 ** params["M"]) * cf.u.solMass / cf.u.littleh # M in M_solar / h
    rho_m0 = cf.cosmo.rho_m(0) * cf.u.solMass * (cf.u.littleh ** 2) * (cf.u.kpc ** (-3)) # check this!
    R = ((3 * M) / (4 * cf.np.pi * rho_m0)) ** (1/3) # factor of 200?
    R = R.to((cf.u.Mpc / cf.u.littleh), cf.u.with_H0(cf.cosmo.H0 * cf.u.km / (cf.u.s * cf.u.Mpc))) # unit conversion
    sigma = cf.cosmo.sigma(R.value, z)
    delta_c = 1.686 #47 
    nu = delta_c / sigma
    bias = (cf.hm.bias.Tinker10(nu, n = cf.cosmo.ns, sigma_8 = cf.cosmo.sigma8, \
                               cosmo = cf.astropy_cosmo)).bias()
    return bias

def lnprior_flat(params):
    if params["M"] > 6.0 and params["M"] < 18.0:
        return 0.0
    else:
        return -cf.np.inf
    
def log_likelihood_flat(params, z, bias, bias_err):
    params_loc = dict({"M": params[0]})
    #print(params_loc)
    lp = lnprior_flat(params_loc)
    if not cf.np.isfinite(lp):
        return -cf.np.inf
    bias_model = calc_halo_bias(params_loc, z)
    return lp - 0.5 * cf.np.sum(((bias - bias_model) / bias_err) ** 2)

def plot_corner(backend_filename = "M_halo_test", print_autocorr_time = True):
    
    fid_pars = [12.8]
    cf.os.chdir("/Users/user/Documents/PGR/UDS field/MCMC backends/")
    reader = cf.emcee.backends.HDFBackend(backend_filename + ".h5")
    # print(reader.get_last_sample())
    if print_autocorr_time == True:
        tau = reader.get_autocorr_time()
        print("autocorr time = " + str(tau))
    flat_samples = reader.get_chain(flat = True) # discard=100, thin=15,
    
    labels = ["M", "z"]
    fig = cf.corner.corner(flat_samples, labels = labels,  \
            quantiles = [0.16, 0.5, 0.84], show_titles = True, title_kwargs = {"fontsize": 12}) 

def calculate_M_h(z, gal_bias, ndim = 1, nwalkers = 100, nsteps = 5000, \
                  backend_filename = "M_halo_test", fid_pars = [12.8]):
    
    # M = {}
    # M["M"] = 12
    # print(calc_halo_bias(M, z))
    # print(log_likelihood_flat(fid_pars, z, gal_bias[0], gal_bias[1]))
    
    # define initial positions
    pos = [fid_pars + 1e-2 * cf.np.random.uniform(0, 1, ndim) * fid_pars for i in range(nwalkers)]

    with cf.Pool() as pool:
        cf.os.chdir("/Users/user/Documents/PGR/UDS field/MCMC backends/")
        backend = cf.emcee.backends.HDFBackend(backend_filename + ".h5")
        sampler = cf.emcee.EnsembleSampler(nwalkers, ndim, log_likelihood_flat, \
                args = (z, gal_bias[0], gal_bias[1]), pool = pool, backend = backend)
        sampler.run_mcmc(pos, nsteps, progress = True)
        print("Final size: {0}".format(backend.iteration))

if __name__ == "__main__":
    calculate_M_h(0.98, [1.8, 0.4])
    plot_corner()
    
    