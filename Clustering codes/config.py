#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 14:38:53 2021

@author: u92876da
"""

# global variables and imports
    
# imports
import os
from abc import ABC, abstractmethod
import camb
from halomod import TracerHaloModel
from halomod.integrate_corr import AngularCF
import halomod as hm
import time
import emcee
import corner
import inspect
from types import FunctionType
# import timeit
import dill as pickle
import json
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from multiprocessing import Pool
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astroML_v2 as astroML
#from astroML.correlation import bootstrap_two_point_angular
import scipy.integrate
import math
from astropy import units as u
from astropy import constants as const
import matplotlib as mpl
from scipy.optimize import curve_fit
from numpy.polynomial import Polynomial
from matplotlib import cm
from matplotlib.lines import Line2D
from dataclasses import dataclass
from enum import Enum
import sklearn.utils
#import halotools
#from halotools.mock_observables import angular_tpcf
mpl.colors.Normalize(vmin=0, vmax=1)
plasma = cm.get_cmap("plasma")
tab10 = cm.get_cmap("tab10")
binary = cm.get_cmap("binary")
from astropy.cosmology import FlatLambdaCDM # does not include sigma8 or n_s
from colossus.cosmology import cosmology
astropy_cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3, Ob0 = 0.05, Tcmb0=2.725) # set astropy cosmology
cosmo = cosmology.setCosmology("cosmo", dict(H0 = 70, Om0 = 0.3, Ob0 = 0.05, sigma8 = 0.9, ns = 1.0)) # set colossus cosmology
# define global variables
filenames = {"Rachana UDS DR11": "/Users/user/Documents/PGR/UDS field/DR11-2arcsec-Jun-30-2019.fits"}
UDS_mask_filename = "DR11.multiband-mask.final.binary-best.fits"
masks = {"UDS": "/Users/user/Documents/PGR/UDS field/" + UDS_mask_filename}

def rewrite_matched_table_columns(data, save_name = "DR11 UDS with Nathans cuts"): # remove _1 from columns
    col_names = data.columns.names
    for i in range(len(col_names)):
        if col_names[i][-2:] == "_1":
            data.columns[i].name = col_names[i][:-2]
    # save table
    return data

class galaxy_sample_type(Enum): # applies to data in a given z-bin
    ALL, LOWER_STELLAR_MASS_LIM, STELLAR_MASS_BIN = range(3)
    
def linear_scale_factor(z):
    
    print("From colossus cosmology: D(z = %1.1f) = %1.5f" % (z, cosmo.growthFactor(z)))
    
    # Carroll 1992
    a = 1 / (1 + z)
    def D(a_prime):
        def da_dτ(a_prime):
            return np.sqrt(1 + cosmo.Om0 * (1 / a_prime - 1) + (1 - cosmo.Om0) * (a_prime ** 2 - 1))
        def integrand(a_prime):
            return np.power(da_dτ(a_prime), -3)
        integral = scipy.integrate.quad(integrand, 0, a_prime)[0]
        return integral * da_dτ(a_prime) * 5 * cosmo.Om0 / (2 * a_prime)
    D_a = D(a) / D(1)
    print("From Carroll 1992: D(z = %1.1f) = %1.5f" % (z, D_a))


#UDS_solid_angle = 0.63 # deg^2 FROM THE PRE-MASK DATA = 0.8622561357749562
dm_cosmic_var_file = "/Users/user/Documents/PGR/UDS field/UDS DM cosmic variance.csv"
M_min_params_file = "/Users/user/Documents/PGR/Literature/Data/M_min params.csv" # 1997 paper

if __name__ == "__main__":
    linear_scale_factor(200)
