# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 13:03:42 2021

@author: dunca
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy
import pandas as pd
import scipy.integrate
import scipy.stats
# astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=65, Om0=0.4)
Omega_0 = cosmo.Om0 #+ cosmo.Ogamma0 + cosmo.Onu0 # generic density parameter
from astropy import units as u
from astropy import constants as const
from astropy.io import fits
from astroML.correlation import bootstrap_two_point_angular
from astroML.correlation import two_point_angular
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
#from astroML.utils.decorators import pickle_results
from scipy.optimize import curve_fit
import pickle
import time
import math
import sklearn.utils
import matplotlib as mpl
from matplotlib import cm
mpl.colors.Normalize(vmin=0, vmax=1)
plasma = cm.get_cmap("plasma")

# global variables
UDS_field = 'file:///C:/Users/dunca/Documents/PGR/HSCDR3/Catalogues/UDS_HSC3_Cut.fits'
UDS_solid_angle = 0.63 # deg^2
redshift_width = 0.3

def read_in_fits_table(fits_filename):
    HDU_list = fits.open(fits_filename)
    #HDU_list.info()
    data = HDU_list[1].data
    #print(UDS_data.columns)
    return(data)

def calculate_comoving_volume(z1, z2, solid_angle): # z2 > z1; input solid angle in deg**2
    def integrand(u):
        integrand = cosmo.differential_comoving_volume(u).value # potentially shouldn't be comoving volume
        return integrand
    comoving_volume_element, err = scipy.integrate.quad(integrand, z1, z2)
    comoving_volume = comoving_volume_element * (u.Mpc ** 3 / u.sr) * (solid_angle * u.deg ** 2).to(u.sr) # in Mpc**3
    return comoving_volume

def photo_z_group(orig_table, redshift, width): # make new column of photo-z grouped galaxies
    orig_cols = orig_table.columns
    # loop through photoZ data
    # store boolean array of whether galaxy is at required photo-z
    is_correct_photo_z = []
    for photo_z in orig_table['z_p']:
        if photo_z <= (redshift + width) and photo_z >= (redshift - width):
            is_correct_photo_z.append(1)
        else:
            is_correct_photo_z.append(0)
    new_cols = fits.ColDefs([
            fits.Column(name='z_p=%1.1f' %redshift, format='L', array = is_correct_photo_z)])
    new_table = fits.BinTableHDU.from_columns(orig_cols + new_cols)
    return new_table.data

def stellar_mass_function(data, redshift, obs_comoving_volume, stellar_mass_of_specific_metallicity, n_bins):
    # create pandas dataframe
    df = pd.DataFrame({"is_correct_z": data['z_p=%1.1f'%redshift],\
                           "log10_stellar_mass": np.log10(data[stellar_mass_of_specific_metallicity])})
    # retain only galaxies with a suitable stellar mass
    df = df.loc[lambda temp_df: temp_df["log10_stellar_mass"] > 0]
    # retain only galaxies with correct redshift
    correct_z = df.loc[lambda temp_df: temp_df["is_correct_z"] == True]
    
    # array of n_bins equally spaced in log space between the min. and max. stellar masses
    stellar_mass_bins = np.linspace(correct_z["log10_stellar_mass"].min(),\
                                    correct_z["log10_stellar_mass"].max(), num = n_bins)
    # scale the histogram y-axis to form a number density
    weights = np.full(len(correct_z), n_bins / obs_comoving_volume.value)
    
    stellar_mass_func, bin_edges = np.histogram(correct_z["log10_stellar_mass"], bins = n_bins, weights = weights)
    #stellar_mass_func = stellar_mass_func * n_bins / obs_comoving_volume
    #print(bin_edges)
    
#    # split up the data into n_bins based on stellar mass of specific metallicity stars
#    edge_to_mid_bin = (correct_z["log10_stellar_mass"].max() - correct_z["log10_stellar_mass"].min()) / (2 * n_bins)
#    correct_z["mean_bin_log10_stellar_mass"] = pd.cut(correct_z["log10_stellar_mass"], bins = n_bins,\
#                                   labels = np.linspace(correct_z["log10_stellar_mass"].min() + edge_to_mid_bin, \
#                                   correct_z["log10_stellar_mass"].max() - edge_to_mid_bin, num = n_bins))
#    #print(correct_z)
#    # create an ordered array containing the mid-bin points
#    stellar_mass_bins = []
#    for i in correct_z["mean_bin_log10_stellar_mass"]:
#        if not i in stellar_mass_bins:
#            stellar_mass_bins = np.append(stellar_mass_bins, i)
#    stellar_mass_bins = np.sort(stellar_mass_bins)
#    #print(stellar_mass_bins)
    return (bin_edges[:-1], stellar_mass_func)

def plot_stellar_mass_function(data, redshift, z_width, solid_angle, hist_bins):
    
    # create group of galaxies at required redshift
    data = photo_z_group(data, redshift, z_width)
    
    # calculate co-moving volume
    comoving_volume = calculate_comoving_volume(redshift - z_width, redshift + z_width, solid_angle)
    
    # calculate stellar mass functions for different redshift, stellar metallicity
    z_22_bins, z_22_hist = stellar_mass_function(data, redshift, comoving_volume, 'Mstar_m22_z_p', hist_bins)
    z_32_bins, z_32_hist = stellar_mass_function(data, redshift, comoving_volume, 'Mstar_m32_z_p', hist_bins)
    z_42_bins, z_42_hist = stellar_mass_function(data, redshift, comoving_volume, 'Mstar_m42_z_p', hist_bins)
    z_52_bins, z_52_hist = stellar_mass_function(data, redshift, comoving_volume, 'Mstar_m52_z_p', hist_bins)
    z_62_bins, z_62_hist = stellar_mass_function(data, redshift, comoving_volume, 'Mstar_m62_z_p', hist_bins)
    z_72_bins, z_72_hist = stellar_mass_function(data, redshift, comoving_volume, 'Mstar_m72_z_p', hist_bins)
    
    # plot the graphs
    plt.plot(z_22_bins, np.log10(z_22_hist), label = 'Metallicity = 22', lw = 2)
    plt.plot(z_32_bins, np.log10(z_32_hist), label = 'Metallicity = 32', lw = 2)
    plt.plot(z_42_bins, np.log10(z_42_hist), label = 'Metallicity = 42', lw = 2)
    plt.plot(z_52_bins, np.log10(z_52_hist), label = 'Metallicity = 52', lw = 2)
    plt.plot(z_62_bins, np.log10(z_62_hist), label = 'Metallicity = 62', lw = 2)
    plt.plot(z_72_bins, np.log10(z_72_hist), label = 'Metallicity = 72', lw = 2)
    
    # plot parameters
    plt.title('z=%1.1f stellar mass functions' %redshift, fontsize = 14)
    plt.xlabel('log$_{10}$(M/M$_{\odot}$)', fontsize = 12)
    cf.plt.ylabel("Arbitrary units (Havn't checked the units yet)")
    #plt.ylabel('log$_{10}$($\Phi$/Mpc$^{-3}$dex$^{-1}$)', fontsize = 12)
    plt.xlim(10, 12)
    plt.legend(loc="upper right", frameon=True, fontsize=10, bbox_to_anchor=(1, 1), fancybox=True, shadow=True)
    plt.xticks(fontsize = 10)
    plt.yticks(fontsize = 10)
    plt.yscale("log")
    plt.savefig('z=%1.1f stellar mass functions.jpeg' %redshift, dpi=600)  # save a copy of the plot
    plt.show()
    return 0

def calculate_number_counts(data, redshift, z_width):
    
    # remove this later on and put in a separate function
    data = photo_z_group(data, redshift, z_width)
    
    # create pandas dataframe
    df = pd.DataFrame({"is_correct_z": data['z_p=%1.1f'%redshift],\
                           "K_AB": data['MAG_BEST_K']})
    # retain only galaxies with a suitable stellar mass
    #df = df.loc[lambda temp_df: temp_df["log10_stellar_mass"] > 0]
    # retain only galaxies with correct redshift
    correct_z = df.loc[lambda temp_df: temp_df["is_correct_z"] == True]
    
    # define histogram bin edges with half magnitude width
    half_mag_bin_edges = np.linspace(18, 27, 19)
    #print(half_mag_bin_edges)
    number_counts, bin_edges = np.histogram(correct_z["K_AB"], bins = half_mag_bin_edges)
    
    # calculate mid-bin points
    mid_bin_pts = bin_edges[:-1] + 0.25
    
    # calculate Poisson errors
    y_errs = np.sqrt(number_counts)
    
    return mid_bin_pts, number_counts, y_errs

def plot_number_counts(data, redshifts, z_width, colors, shapes):
    if len(redshifts) == len(colors) == len(shapes):
        for i in range(len(redshifts)):
            app_K_mag, counts, counts_err = calculate_number_counts(data, redshifts[i], z_width)
            plt.errorbar(app_K_mag, counts, yerr = counts_err, color = colors[i], label = "z = %1.1f" %redshifts[i],\
                         fmt = shapes[i], capsize = 5, capthick = 1)
        # plot parameters
        plt.xlabel("K$_{AB}$", fontsize = 12)
        plt.ylabel("N / 0.5 mag / deg$^2$", fontsize = 12)
        plt.legend(loc="upper right", frameon=True, fontsize=10, bbox_to_anchor=(1, 1), fancybox=True, shadow=True)
        plt.xticks(fontsize = 10)
        plt.yticks(fontsize = 10)
        plt.yscale("log")
        if len(redshifts) == 1:
            plt.title('z=%1.1f differential number counts.jpeg' %redshifts[0], fontsize = 14)
            plt.savefig('z=%1.1f differential number counts.jpeg' %redshifts[0], dpi=600)
        else:
            plt.title('Differential number counts for z=%1.1f-%1.1f.jpeg' %(np.min(redshifts), np.max(redshifts)), fontsize = 14)
            plt.savefig('Differential number counts for z=%1.1f-%1.1f.jpeg' %(np.min(redshifts), np.max(redshifts)), dpi=600)
        plt.show()
    else:
        print("Your redshift array must be the same length as your colors and shapes arrays!")
    return 0

def calculate_2_pt_angular_correlation(data, redshift, z_width, min_bin, max_bin, Nbins, Nbootstraps, pickle_name, n_gals = None):
    # select only galaxies in required redshift range
    data = photo_z_group(data, redshift, z_width)
    data = data[data["z_p=%1.1f" %redshift] == 1]
    if n_gals == None:
        n_gals = len(data)
    else: # take a random sample of galaxies of length n_gals
        np.random.seed()
        row_numbers_to_keep = np.random.choice(len(data), n_gals, replace=False)
        data = data[row_numbers_to_keep]
        print(len(data))
    # compute results from data     
    # create an array of bins evenly spaced in log-space
    bins = 10 ** np.linspace(np.log10(min_bin), np.log10(max_bin), Nbins)
    # calculate w(theta) from data using bootstrapping to generate errors
    if Nbootstraps != 0:
        results = bootstrap_two_point_angular(data["ALPHA_J2000"], data["DELTA_J2000"], bins = bins,\
                                              method = "landy-szalay", Nbootstraps = Nbootstraps)
        corr, corr_err, bootstraps = results
    else:
        results = two_point_angular(data["ALPHA_J2000"], data["DELTA_J2000"], bins = bins, method = "landy-szalay")
        corr = results
        corr_err = np.zeros(len(results))
        bootstraps = 0
        
    mid_bins = 0.5 * (bins[1:] + bins[:-1])
    
    # save object as pickle
    with open(pickle_name, 'wb') as f:
        pickle.dump([mid_bins, corr, corr_err, bootstraps, redshift, n_gals], f)
    
    return mid_bins, corr, corr_err, bootstraps, n_gals

def plot_2_pt_angular_correlation(pickle_names, colors, shapes, plot_title, fit_models,\
                                  fitting_params, int_constraints, int_constraint_errs = 0):
    # empty arrays storing model parameters
    parameters = []
    # open data file to plot
    for i in range(len(pickle_names)):
        with open(pickle_names[i], 'rb') as f:
            theta, corr, corr_err, bootstraps, redshift, galaxy_number = pickle.load(f)
        # plot 2-pt angular correlation function data
        plt.errorbar(theta, corr, yerr = corr_err, label = ("z = %1.1f, N = %5.0f" % (redshift, galaxy_number)),\
                     c = colors[i], fmt = shapes[i], capsize = 5, capthick = 1, ls = "none")
        # plot vertical 1 Mpc/h line
        theta_1Mpc_h = (cosmo.arcsec_per_kpc_comoving(redshift)).to(u.deg / u.Mpc)
        plt.axvline(x = theta_1Mpc_h.value, color = colors[i], linestyle = "--")

        if fit_models == True:
            theta, model, popt, pcov, red_chi_sq = \
            model_2_pt_angular_correlation(fitting_params, pickle_names[i], int_constraints[i])
            # store popt, pcov, red_chi_sq for each data set
            if fitting_params == "A, delta, C":
                params_loc = np.append(np.concatenate((popt, np.diag(pcov))), red_chi_sq)
            elif fitting_params == "A, C":
                params_loc = (popt[0], 0.8, popt[1], pcov[0][0], 0, pcov[1][1], red_chi_sq)
            elif fitting_params == "A, delta":
                params_loc = (popt[0], popt[1], int_constraints[i], pcov[0][0], pcov[1][1], int_constraint_errs[i], red_chi_sq)
            elif fitting_params == "A":
                params_loc = (popt[0], 0.8, int_constraints[i], pcov[0][0], 0, 0, red_chi_sq)
            parameters = np.append(parameters, params_loc)
            # print model and fitting values and errors
            if len(pickle_names) == 1:
                if fitting_params == "A, delta, C":
                    fit_label = 'fit: A=%1.2e\u00B1%1.2e,\n     \u03B4=%5.3f\u00B1%5.3f,\n     C=%5.3f\u00B1%5.3f, \n     \u03C7$^2$ red. = %5.3f'\
                         % (popt[0], np.sqrt(pcov[0][0]), popt[1], np.sqrt(pcov[1][1]), popt[2], np.sqrt(pcov[2][2]), red_chi_sq)
                elif fitting_params == "A, C":
                    fit_label = 'fit: A=%1.2e\u00B1%1.2e,\n     \u03B4=0.8,\n     C=%5.3f\u00B1%5.3f, \n     \u03C7$^2$ red. = %5.3f'\
                         % (popt[0], np.sqrt(pcov[0][0]), popt[1], np.sqrt(pcov[1][1]), red_chi_sq)
                elif fitting_params == "A, delta":
                    fit_label = 'fit: A=%1.2e\u00B1%1.2e,\n     \u03B4=%5.3f\u00B1%5.3f,\n     C=%5.3f\u00B1%5.3f, \n     \u03C7$^2$ red. = %5.3f'\
                         % (popt[0], np.sqrt(pcov[0][0]), popt[1], np.sqrt(pcov[1][1]), int_constraints[i], int_constraint_errs[i], red_chi_sq)
                elif fitting_params == "A":
                    fit_label = 'fit: A=%1.2e\u00B1%1.2e,\n     \u03B4=0.8,\n     C=%5.3f\u00B1%5.3f, \n     \u03C7$^2$ red. = %5.3f'\
                         % (popt[0], np.sqrt(pcov[0][0]), int_constraints[i], int_constraint_errs[i], red_chi_sq)
            else:
                fit_label = None
            # plot the fit
            plt.plot(theta, model, color = colors[i], label = fit_label)
    # plot parameters
    plt.title(plot_title, fontsize = 14)
    plt.xlabel("$\u03B8$ /deg", fontsize = 12)
    plt.ylabel("\u03C9($\u03B8$)", fontsize = 12)
    plt.legend(loc="upper right", frameon=True, fontsize=10, bbox_to_anchor=(1, 1), fancybox=True, shadow=True)
    plt.xticks(fontsize = 10)
    plt.yticks(fontsize = 10)
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig(plot_title + ".jpeg", dpi=600)
    plt.show()
    # return based on input
    if fit_models == True:
        parameters = np.reshape(parameters, (len(pickle_names), 7))
        return parameters
    else:
        return 0

def model_2_pt_angular_correlation(fitting_params, pickle_names, int_constraint = 0):
    # 2-pt angular correlation function model
    if fitting_params == "A, delta, C":
        def two_pt_angular_fitting_func(theta, A, delta, C):
            return A * (theta ** (- delta) - C)
    elif fitting_params == "A, C": # delta = 0.8
        def two_pt_angular_fitting_func(theta, A, C):
            return A * (theta ** (- 0.8) - C)
    elif fitting_params == "A, delta": # C constrained
        def two_pt_angular_fitting_func(theta, A, delta):
            return A * (theta ** (- delta) - int_constraint)
    elif fitting_params == "A": # C constrained and delta = 0.8
        def two_pt_angular_fitting_func(theta, A):
            return A * (theta ** (- 0.8) - int_constraint)
    else:
        print("Fitting parameters incorrectly specified!")
        raise SystemExit
    # fit the model to the data
    with open(pickle_names, 'rb') as f:
        theta, corr, corr_err, bootstraps, redshift, galaxy_number = pickle.load(f)
    popt, pcov = curve_fit(two_pt_angular_fitting_func, theta, corr, sigma = corr_err,\
                           absolute_sigma = True, method = "lm")
    # compute the reduced chi squared value
    model = two_pt_angular_fitting_func(theta, *popt)
    chi_sq = 0
    for i in range(len(theta)):
        chi_sq += ((corr[i] - model[i]) / corr_err[i]) ** 2
    red_chi_sq = chi_sq / (len(theta) - len(popt))
    #print("\u03C7$^2$ red. = %5.3f" % red_chi_sq)
    return theta, model, popt, pcov, red_chi_sq

def galaxy_pairs(data, correlation_type, redshift, z_width, min_bin, max_bin, Nbins, pickle_name, seed = [None, None], n_random = None):
    
    # store the min and max ra and dec for the data field before the redshift cut
    min_ra = np.min(data["ALPHA_J2000"])
    max_ra = np.max(data["ALPHA_J2000"])
    min_dec = np.min(data["DELTA_J2000"])
    max_dec = np.max(data["DELTA_J2000"])
    
    # store the number of galaxies
    number_of_gals = len(data)
    #print(number_of_gals)
    # create an array of theta bin edges evenly spaced logarithmically
    bins = 10 ** np.linspace(np.log10(min_bin), np.log10(max_bin), Nbins)   
    
    # catch the error if correlation type is incorrectly inputted
    if correlation_type == "data-data" or correlation_type == "data-random" or correlation_type == "random-random":
        # define correlation fields
        if correlation_type == "data-data":
            field_1 = SkyCoord(ra = data["ALPHA_J2000"] * u.deg, dec = data["DELTA_J2000"] * u.deg)
            field_2 = field_1
            total_pairs = number_of_gals * (number_of_gals - 1) / 2 # total data-data pairs for normalization
        else: # randomly place the same number of galaxies within the RA, Dec bounds of the survey
            if n_random == None:
                n_random = number_of_gals
            np.random.seed(seed[0]) # set random number seed
            random_ra = np.random.uniform(min_ra, max_ra, size = n_random)
            np.random.seed(seed[1]) # set different random number seed
            random_dec = np.random.uniform(min_dec, max_dec, size = n_random)
            field_1 = SkyCoord(ra = random_ra * u.deg, dec = random_dec * u.deg)
            if correlation_type == "data-random":
                field_2 = SkyCoord(ra = data["ALPHA_J2000"] * u.deg, dec = data["DELTA_J2000"] * u.deg)
                total_pairs = number_of_gals * n_random # total data-random pairs for normalization
            elif correlation_type == "random-random":
                field_2 = SkyCoord(ra = random_ra * u.deg, dec = random_dec * u.deg)
                total_pairs = n_random * (n_random - 1) / 2 # total random-random pairs for normalization
    else:
        print("Incorrect input of the variable 'correlation_type'")
    
    # calculate the number of pairs between the data and/or random fields in each theta bin
    sep2d = field_2.search_around_sky(field_1, 100 * u.deg)[2]
    sep2d = [x for x in sep2d.value if x != 0] # remove all occurrences of 0
    # add total number of pairs of galaxies for each angular separation bin to 'number_of_pairs'
    number_of_pairs, bin_edges = np.histogram(sep2d, bins = bins)
    if correlation_type != "data-random":
        number_of_pairs = number_of_pairs / 2 # avoid double counting
    
    # calculate the theta values at the bin centres
    mid_bins = 0.5 * (bins[1:] + bins[:-1])
    
    return mid_bins, number_of_pairs, total_pairs

def overplot_2_pt_angular(parameters, min_bin, max_bin, Nbins, label):
    
    # unpack parameters
    A = parameters[0]
    delta = parameters[1]
    C = parameters[2]
    
    bins = 10 ** np.linspace(np.log10(min_bin), np.log10(max_bin), Nbins)
    # calculate the theta values at the bin centres
    mid_bins = 0.5 * (bins[1:] + bins[:-1])
    
    # calculate values of the model at these bins
    model = A * (np.power(mid_bins, -delta) - C)
    
    plt.plot(mid_bins, model, color = "blue", linestyle = "--", label = label)
    return 0

def calculate_r0(amplitude, delta, redshift, delta_z): # for a flat Universe with non-zero cosmological const
    
    H_gamma = math.gamma(1/2) * math.gamma(delta / 2) / math.gamma((delta + 1) / 2)
    print("H_gamma = %1.2f" %H_gamma)
    F_z = 1 # curvature correction factor
    x = cosmo.comoving_distance(redshift) # comoving distance
    print(x)
    P_Om0_z = np.sqrt(Omega_0 * ((1 + redshift) ** 3 + (1 / Omega_0) - 1))
    print(P_Om0_z) # checked!
    rms_z_err = 0.1 * (1 + redshift)
    z_bin_width_increase = np.sqrt(12 * np.power(rms_z_err / (delta_z), 2) + 1) # assuming Gaussian errors
    delta_z = delta_z / z_bin_width_increase
    r_0 = ((const.c * amplitude * delta_z) / (cosmo.H0 * H_gamma * (x ** (-delta)) * P_Om0_z * F_z)) ** (1 / (delta + 1))
    r_0 = r_0.to(u.Mpc / u.littleh, u.with_H0(cosmo.H0)) # unit conversion
    print("r_0 = %2.2f Mpc/h" % r_0.value)
    r_0_err = 0
    return r_0, r_0_err

# manually check that x in the "calculate_r0" function is the co-moving distance
def comoving_distance_check(redshift):
    def integrand(u):
        integrand = 1 / np.sqrt(np.power(1 + u, 3) + (1 / Omega_0) - 1)
        return integrand
    integral, err = scipy.integrate.quad(integrand, 0, redshift)
    comoving_distance = integral * const.c / (cosmo.H0 * np.sqrt(Omega_0))
    comoving_distance = comoving_distance.to(u.Mpc / u.littleh, u.with_H0(cosmo.H0))
    return comoving_distance # CHECK COMPLETE

def redshift_distribution(data, redshift, z_width, Nbins, plot_results):
    # calculate histogram of the data
    # select only galaxies in required redshift range
    data = photo_z_group(data, redshift, z_width)
    data = data[data["z_p=%1.1f" %redshift] == 1]
    # bin the data by redshift
    z_bins = np.linspace(redshift - z_width, redshift + z_width, Nbins)
    z_mid_bins = 0.5 * (z_bins[1:] + z_bins[:-1])
    N_z_data, z_bin_edges = np.histogram(data["z_p"], bins = z_bins)
    
    # fit polynomial to histogram
    
    # plot the histogram and fit
    if plot_results == True:
        plt.bar(z_mid_bins, N_z_data, width = z_bins[1] - z_bins[0], fill = True,\
                label = "z = %1.1f, N = %5.0f" % (redshift, len(data)))
        # plot parameters
        plt.title("z = %1.1f \u00B1 %1.1f redshift distribution" % (redshift, z_width), fontsize = 14)
        plt.xlabel("z", fontsize = 12)
        plt.ylabel("N(z)", fontsize = 12)
        plt.xticks(fontsize = 10)
        plt.yticks(fontsize = 10)
        plt.xlim(redshift - z_width, redshift + z_width)
        plt.legend(loc="upper right", frameon=True, fontsize=10, bbox_to_anchor=(1, 1), fancybox=True, shadow=True)
        plt.savefig("z = %1.1f \u00B1 %1.1f redshift distribution.jpeg" % (redshift, z_width), dpi=600)
        plt.show()
        
    return z_bins, N_z_data

def manual_2_pt_correlation(data, redshift, z_width, min_bin, max_bin, Nbins, method, Nbootstraps, pickle_name): # APPEND TO CATCH DIVIDE BY ZERO ERRORS AND CORRECT  
        
    # select only galaxies in required redshift range
    data = photo_z_group(data, redshift, z_width)
    data = data[data["z_p=%1.1f" %redshift] == 1]
    
    # save repeated data to then calculate bootstrapped mean and std from
    w_theta_repeats = []
    w_theta_err_repeats = []
    
    for repeat in range(Nbootstraps):
        # resample distribution
        bootstrap_sample = sklearn.utils.resample(data, replace = True) # change sample size here
        
        # calculate DD, DR, and RR galaxy pairs in each theta bin for the sample
        while True: # repeat if there is a divide by zero error
            try:
                theta_bins, data_data, total_data_data = galaxy_pairs(bootstrap_sample, "data-data", redshift,\
                        z_width, min_bin, max_bin, Nbins, "z = %1.1f data_data_pickle.pkl" % redshift)
                # seed for random field
                seed = np.random.randint(0, 2 ** 32, 2, dtype = np.uint32)
                print(seed)
                theta_bins, data_random, total_data_random = galaxy_pairs(bootstrap_sample, "data-random", redshift,\
                        z_width, min_bin, max_bin, Nbins, "z = %1.1f data_random_pickle.pkl" % redshift, seed = seed)
                
                theta_bins, random_random, total_random_random = galaxy_pairs(bootstrap_sample, "random-random", redshift,\
                        z_width, min_bin, max_bin, Nbins, "z = %1.1f random_random_pickle.pkl" % redshift, seed = seed)
                # normalize the data, random counts
                data_data_normed = data_data / total_data_data
                data_random_normed = data_random / total_data_random
                random_random_normed = random_random / total_random_random
                if method == "landy-szalay":
                    w_theta = (data_data_normed - 2 * data_random_normed + random_random_normed) \
                        / random_random_normed
                    # Poisson errors
                    w_theta_err = np.sqrt(((data_data / (total_data_data ** 2)) + (data_random / (total_data_random ** 2)) \
                                        + (random_random / (total_random_random ** 2)) / np.power(data_data_normed \
                                        - 2 * data_random_normed + random_random_normed, 2)) + (1 / random_random)) * w_theta
                break
            except ZeroDivisionError:
                print("Division by zero!")
        
        # add repeat to 'global' w_theta, w_theta_err
        w_theta_repeats = np.append(w_theta_repeats, w_theta)
        w_theta_err_repeats = np.append(w_theta_err_repeats, w_theta_err)
        
    # reshape 'global' w_theta, w_theta_err
    w_theta_repeats = (np.reshape(w_theta_repeats, (Nbootstraps, Nbins - 1))).T
    w_theta_err_repeats = (np.reshape(w_theta_err_repeats, (Nbootstraps, Nbins - 1))).T
    
    # calculate a weighted mean and error of w_theta for each theta bin
    mean = []
    mean_err = []
    for i in range(Nbins - 1):
        # perform weighted mean including Poisson errors
        #weights = 1 / np.power(w_theta_err_repeats[i], 2)
        mean_loc = np.mean(w_theta_repeats[i]) #np.average(w_theta_repeats[i], weights = weights)
        mean_err_loc = np.std(w_theta_repeats[i]) #np.sqrt(1 / np.sum(weights))
        mean = np.append(mean, mean_loc)
        mean_err = np.append(mean_err, mean_err_loc)
    
    # save object as pickle
    bootstraps = [w_theta_repeats, w_theta_err_repeats]
    with open(pickle_name, 'wb') as f:
        pickle.dump([theta_bins, mean, mean_err, bootstraps, redshift, len(bootstrap_sample)], f)
        
    return theta_bins, mean, mean_err

def repeat_2_pt_angular_correlation(data, redshift, z_width, min_bin, max_bin, Nbins, method,\
                                                  Nbootstraps, Nrepeats, calc_type): # IMPLEMENT WEIGHTED MEAN HERE
    w_theta_repeats = []
    w_theta_err_repeats = []
    # repeat manual 2-pt correlation calculations and save results
    for repeat in range(Nrepeats):
        print(repeat)
        if calc_type == "manual":
            theta_bins, w_theta, w_theta_err = manual_2_pt_correlation(data, redshift, z_width, min_bin, max_bin, Nbins, method, Nbootstraps)
        elif calc_type == "function":
            theta_bins, w_theta, w_theta_err, bootstraps, n_gals = calculate_2_pt_angular_correlation(data, redshift, z_width,\
                                                                   min_bin, max_bin, Nbins, Nbootstraps, "repeat.pkl")
        # add repeat to 'global' w_theta, w_theta_err
        w_theta_repeats = np.append(w_theta_repeats, w_theta)
        w_theta_err_repeats = np.append(w_theta_err_repeats, w_theta_err)
    
    # reshape 'global' w_theta, w_theta_err
    w_theta_repeats = (np.reshape(w_theta_repeats, (Nrepeats, Nbins - 1))).T
    w_theta_err_repeats = (np.reshape(w_theta_err_repeats, (Nrepeats, Nbins - 1))).T
    # calculate a weighted mean and error of w_theta for each theta bin
    mean = []
    mean_err = []
    for i in range(Nbins - 1):
        if calc_type == "manual": # perform weighted mean including Poisson errors
            weights = 1 / np.power(w_theta_err_repeats[i], 2)
            mean_loc = np.average(w_theta_repeats[i], weights = weights)
            mean_err_loc = np.sqrt(1 / np.sum(weights))
        elif calc_type == "function":
            mean_loc = np.mean(w_theta_repeats[i])
            mean_err_loc = np.std(w_theta_repeats[i])
        mean = np.append(mean, mean_loc)
        mean_err = np.append(mean_err, mean_err_loc)
    
    # save object as pickle
    with open("z = %1.0f " %(redshift) + calc_type + " 2pt angular correlation.pkl", 'wb') as f:
        pickle.dump([theta_bins, mean, mean_err, w_theta_repeats, w_theta_err_repeats], f)
            
    return theta_bins, mean, mean_err, w_theta_repeats, w_theta_err_repeats

def plot_repeated_2_pt_angular_residuals(data, redshift, z_width, min_bin, max_bin, Nbins, method,\
                                                  Nbootstraps, Nrepeats):
    theta_bins, func_mean, func_std, func_repeats, func_repeats_err = repeat_2_pt_angular_correlation(data, redshift, z_width, min_bin, max_bin, Nbins,\
                                                                      "landy-szalay", Nbootstraps, Nrepeats, "function")
    theta_bins, manual_mean, manual_std, manual_repeats, manual_repeats_err = repeat_2_pt_angular_correlation(data, redshift, z_width, min_bin, max_bin, Nbins,\
                                                                          "landy-szalay", Nbootstraps, Nrepeats, "manual")
    # plot w(theta) comparison
    plt.figure()
    plt.errorbar(theta_bins, func_mean, func_std, label = "z = %1.1f function, %1.0f repeats" % (redshift, Nrepeats), ls = "none",\
                 fmt = "s", capsize = 5, capthick = 1)
    plt.errorbar(theta_bins, manual_mean, manual_std, label = "z = %1.1f manual, %1.0f repeats" % (redshift, Nrepeats), ls = "none",\
                 fmt = "s", capsize = 5, capthick = 1)
    plt.legend(loc="upper right", frameon=True, fontsize=10, bbox_to_anchor=(1, 1), fancybox=True, shadow=True)
    plt.xscale("log")
    plt.yscale("log")
    plt.show()
    
    # calculate residuals
    residuals = func_mean - manual_mean
    residuals_err = np.sqrt((func_std)**2 + (manual_std)**2)

    plt.figure()
    plt.errorbar(theta_bins, residuals, residuals_err, label = "z = %1.1f residuals, %1.0f repeats" % (redshift, Nrepeats), ls = "none",\
                 fmt = "s", capsize = 5, capthick = 1, c = "black")
    plt.legend(loc="upper right", frameon=True, fontsize=10, bbox_to_anchor=(1, 1), fancybox=True, shadow=True)
    
    chi_sq = 0
    for i in range(len(residuals)):
        if np.isnan(residuals_err[i]) == False:
            chi_sq += (residuals[i] / residuals_err[i]) ** 2
    red_chi_sq = chi_sq / len(residuals)
    print(red_chi_sq)
    
    plt.xscale("log")
    #plt.yscale("log")
    #plt.ylim(1e-100, 1e100)
    plt.show()
    
    return theta_bins, func_mean, func_std, manual_mean, manual_std

def calculate_const_integral_constraint(data, delta, redshift, z_width, min_bin, Nrepeats, n_random = None):
    
    # select only galaxies in required redshift range
    data = photo_z_group(data, redshift, z_width)
    data = data[data["z_p=%1.1f" %redshift] == 1]
    
    # if unspecified number of random galaxies
    if n_random == None:
        n_random = len(data)
    
    # sum extends to largest separations in the field
    C_array = []
    for repeat in range(Nrepeats): # accomodate different random fields
        while True:
            try:
                theta_bins, random_random, total_random_random = galaxy_pairs(data, "random-random", redshift,\
                                    z_width, min_bin, 1e1, 1000, "z = %1.1f random_random_pickle.pkl" % redshift, n_random = n_random)
                
                C = np.sum(random_random * theta_bins ** (- delta)) / total_random_random
                print(C)
                C_array = np.append(C_array, C)
                break
            except ZeroDivisionError:
                print("Division by zero!")
    
    C_mean = np.round(np.mean(C_array), 3)
    C_std = np.round(np.std(C_array), 3)
    
    return C_mean, C_std, C_array

def write_data_to_csv(filename, variables, column_labels, row_labels):
    write_data = pd.DataFrame(data = variables, columns = column_labels, index = row_labels)
    #print(write_data)
    write_data.to_csv(filename + ".csv", mode='w')
    print("Written to..." + filename)
    return write_data

# main code starts here ---------------------------------------------------------------------------

# start time
start_time=time.time()

# read in UDS data
UDS_data = read_in_fits_table(UDS_field)

# plot stellar mass function
#plot_stellar_mass_function(UDS_data, 1, redshift_width, UDS_solid_angle, 20)
#plot_number_counts(UDS_data, [1, 2, 3], redshift_width, ["black", "blue", "red"], ["s", "s", "s"])
#calculate_2_pt_angular_correlation(UDS_data, 4, 0.3, 1e-3, 1e-1, 10, 100, "function z=1 10000 gals.pkl", n_gals = 10000)
#parameters = plot_2_pt_angular_correlation(["z=4 w_theta.pkl"], ["black"], ["s"], "z=4 2-pt angular correlation function", True, "A, delta, C")

#C_mean_z1, C_std_z1, C_array_z1 = calculate_const_integral_constraint(UDS_data, 0.8, 1, 0.3, 1e-3, 10, n_random = 5000)
#C_mean_z2, C_std_z2, C_array_z2 = calculate_const_integral_constraint(UDS_data, 0.8, 2, 0.3, 1e-3, 10, n_random = 5000)
#C_mean_z3, C_std_z3, C_array_z3 = calculate_const_integral_constraint(UDS_data, 0.8, 3, 0.3, 1e-3, 10, n_random = 5000)
#C_mean_z4, C_std_z4, C_array_z4 = calculate_const_integral_constraint(UDS_data, 0.8, 4, 0.3, 1e-3, 10)

#theta_bins, mean, std = repeat_2_pt_angular_correlation(UDS_data, 4, 0.3, 1e-3, 1e-1, 10, "landy-szalay", 0, 5, "function")
#theta_bins, func_mean, func_std, manual_mean, manual_std = plot_repeated_2_pt_angular_residuals(UDS_data, 4, 0.3, 1e-3, 1e-1, 10, "landy-szalay", 0, 100)

#theta_bins, w_theta, w_theta_err = manual_2_pt_correlation(UDS_data, 4, 0.3, 1e-3, 1e-1, 10, "landy-szalay", 100, "manual z=4.pkl")
parameters = plot_2_pt_angular_correlation(["function z=1 10000 gals.pkl", "function z=2 10000 gals.pkl", "function z=3 10000 gals.pkl", "function z=4.pkl"],\
                                           [plasma(0), plasma(0.333333), plasma(0.666666), plasma(0.999999)],\
                                           ["s", "s", "s", "s"], "w(theta) redshift evolution, 100 boots, C constrained, 10000 gals",\
            True, "A", [2.393, 2.397, 2.421, 2.450])#,\
            #int_constraint_err = [0.004, 0.005, 0.007, C_std_z4])

#print(parameters)

# PUT THIS IN THE MODEL_2PT_ANGULAR_CORRELATION FUNCTION
#write_data_to_csv("z=1-4 2-pt angular correlation parameters", parameters,\
#                  ["A", "delta", "C", "A_err", "delta_err", "C_err", "red_chi_sq"], ["z=1", "z=2", "z=3", "z=4"])

#r_0, r_0_err = calculate_r0(5.9e-3, 0.8, 4.4, 0.8)
#redshift_distribution(UDS_data, 2.15, 2.15, 50, True)

# end time
end_time=time.time()
print("This took {} seconds!".format(np.round(end_time-start_time, 2)))
