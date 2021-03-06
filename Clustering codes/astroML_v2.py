#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tools for computing two-point correlation functions.
"""

import numpy as np

from sklearn.neighbors import KDTree
from sklearn.utils import check_random_state
import sklearn.utils
from two_pt_angular_correlation import two_pt_angular_corr
import time


def uniform_sphere(RAlim, DEClim, size=1):
    """Draw a uniform sample on a sphere

    Parameters
    ----------
    RAlim : tuple
        select Right Ascension between RAlim[0] and RAlim[1]
        units are degrees
    DEClim : tuple
        select Declination between DEClim[0] and DEClim[1]
    size : int (optional)
        the size of the random arrays to return (default = 1)

    Returns
    -------
    RA, DEC : ndarray
        the random sample on the sphere within the given limits.
        arrays have shape equal to size.
    """
    zlim = np.sin(np.pi * np.asarray(DEClim) / 180.)

    z = zlim[0] + (zlim[1] - zlim[0]) * np.random.random(size)
    DEC = (180. / np.pi) * np.arcsin(z)
    RA = RAlim[0] + (RAlim[1] - RAlim[0]) * np.random.random(size)

    return RA, DEC


def ra_dec_to_xyz(ra, dec):
    """Convert ra & dec to Euclidean points

    Parameters
    ----------
    ra, dec : ndarrays

    Returns
    x, y, z : ndarrays
    """
    sin_ra = np.sin(ra * np.pi / 180.)
    cos_ra = np.cos(ra * np.pi / 180.)

    sin_dec = np.sin(np.pi / 2 - dec * np.pi / 180.)
    cos_dec = np.cos(np.pi / 2 - dec * np.pi / 180.)

    return (cos_ra * sin_dec,
            sin_ra * sin_dec,
            cos_dec)


def angular_dist_to_euclidean_dist(D, r=1):
    """convert angular distances to euclidean distances"""
    return 2 * r * np.sin(0.5 * D * np.pi / 180.)


def two_point(data, bins, method = "standard", data_R = None, random_state = None):
    """Two-point correlation function

    Parameters
    ----------
    data : array_like
        input data, shape = [n_samples, n_features]
    bins : array_like
        bins within which to compute the 2-point correlation.
        shape = Nbins + 1
    method : string
        "standard" or "landy-szalay".
    data_R : array_like (optional)
        if specified, use this as the random comparison sample
    random_state : integer, np.random.RandomState, or None
        specify the random state to use for generating background

    Returns
    -------
    corr : ndarray
        the estimate of the correlation function within each bin
        shape = Nbins
    """
    data = np.asarray(data)
    bins = np.asarray(bins)
    rng = check_random_state(random_state)

    if method not in ['standard', 'landy-szalay']:
        raise ValueError("method must be 'standard' or 'landy-szalay'")

    if bins.ndim != 1:
        raise ValueError("bins must be a 1D array")

    if data.ndim == 1:
        data = data[:, np.newaxis]
    elif data.ndim != 2:
        raise ValueError("data should be 1D or 2D")

    n_samples, n_features = data.shape

    # shuffle all but one axis to get background distribution
    if data_R is None:
        data_R = data.copy()
        for i in range(n_features - 1):
            rng.shuffle(data_R[:, i])
    else:
        data_R = np.asarray(data_R)
        if (data_R.ndim != 2) or (data_R.shape[-1] != n_features):
            raise ValueError('data_R must have same n_features as data')

    #print(bins)
    # print("len(DD) = %1.1f" % (len(data) * (len(data) - 1) / 2))
    # print("len(RR) = %1.1f" % (len(data_R) * (len(data_R) - 1) / 2))

    # Fast two-point correlation functions added in scikit-learn v. 0.14
    KDT_D = KDTree(data)
    KDT_R = KDTree(data_R)

    counts_DD = KDT_D.two_point_correlation(data, bins)
    counts_RR = KDT_R.two_point_correlation(data_R, bins)

    DD = np.diff(counts_DD)
    RR = np.diff(counts_RR)
    
    factor = len(data_R) * 1. / len(data)
    
    #print(DD)
    #print(sum(DD))
    #print(sum(RR))
    #n_DD = 0.5 + np.sqrt(1 + 8 * sum(DD)) / 2
    #print(n_DD)
    # print(sum(RR))
    #n_RR = 0.5 + np.sqrt(1 + 8 * sum(RR)) / 2

    # check for zero in the denominator
    RR_zero = (RR == 0)
    RR[RR_zero] = 1

    if method == 'standard':
        corr = factor ** 2 * DD / RR - 1
    elif method == 'landy-szalay':
        counts_DR = KDT_R.two_point_correlation(data, bins)

        DR = np.diff(counts_DR)
        #print(DR)
        #print("len(DR) = %1.1f" % (len(data) * len(data_R)))
        #print(sum(DR))

        corr = (factor ** 2 * DD - 2 * factor * DR + RR) / RR

    corr[RR_zero] = np.nan

    return corr



def bootstrap_two_point(data, bins, Nbootstrap=10, method = 'standard',
                        return_bootstraps = False, random_state = None):
    """Bootstrapped two-point correlation function

    Parameters
    ----------
    data : array_like
        input data, shape = [n_samples, n_features]
    bins : array_like
        bins within which to compute the 2-point correlation.
        shape = Nbins + 1
    Nbootstrap : integer
        number of bootstrap resamples to perform (default = 10)
    method : string
        "standard" or "landy-szalay".
    return_bootstraps: bool
        if True, return full bootstrapped samples
    random_state : integer, np.random.RandomState, or None
        specify the random state to use for generating background

    Returns
    -------
    corr, corr_err : ndarrays
        the estimate of the correlation function and the bootstrap
        error within each bin. shape = Nbins
    """
    data = np.asarray(data)
    bins = np.asarray(bins)
    rng = check_random_state(random_state)

    if method not in ['standard', 'landy-szalay']:
        raise ValueError("method must be 'standard' or 'landy-szalay'")

    if bins.ndim != 1:
        raise ValueError("bins must be a 1D array")

    if data.ndim == 1:
        data = data[:, np.newaxis]
    elif data.ndim != 2:
        raise ValueError("data should be 1D or 2D")

    if Nbootstrap < 2:
        raise ValueError("Nbootstrap must be greater than 1")

    n_samples, n_features = data.shape

    # get the baseline estimate
    corr = two_point(data, bins, method=method, random_state=rng)

    bootstraps = np.zeros((Nbootstrap, len(corr)))

    for i in range(Nbootstrap):
        indices = rng.randint(0, n_samples, n_samples)
        bootstraps[i] = two_point(data[indices, :], bins, method=method,
                                  random_state=rng)

    # use masked std dev in case of NaNs
    corr_err = np.asarray(np.ma.masked_invalid(bootstraps).std(0, ddof=1))

    if return_bootstraps:
        return corr, corr_err, bootstraps
    else:
        return corr, corr_err



def two_point_angular(ra, dec, ra_R, dec_R, bins, geometry, method='standard', random_state=None):
    """Angular two-point correlation function

    A separate function is needed because angular distances are not
    euclidean, and random sampling needs to take into account the
    spherical volume element.

    Parameters
    ----------
    ra : array_like
        input right ascention, shape = (n_samples,)
    dec : array_like
        input declination
    bins : array_like
        bins within which to compute the 2-point correlation.
        shape = Nbins + 1
    method : string
        "standard" or "landy-szalay".
    random_state : integer, np.random.RandomState, or None
        specify the random state to use for generating background

    Returns
    -------
    corr : ndarray
        the estimate of the correlation function within each bin
        shape = Nbins
    """
    ra = np.asarray(ra)
    dec = np.asarray(dec)
    rng = check_random_state(random_state)

    if method not in ['standard', 'landy-szalay']:
        raise ValueError("method must be 'standard' or 'landy-szalay'")

    if bins.ndim != 1:
        raise ValueError("bins must be a 1D array")

    if (ra.ndim != 1) or (dec.ndim != 1) or (ra.shape != dec.shape):
        raise ValueError('ra and dec must be 1-dimensional '
                         'arrays of the same length')

    # draw a random sample with N points
    ra_R, dec_R = two_pt_angular_corr.create_random_gals(geometry, len(ra))
    # ra_R, dec_R = uniform_sphere((min(ra), max(ra)),
    #                              (min(dec), max(dec)),
    #                              2 * len(ra))

    data = np.asarray(ra_dec_to_xyz(ra, dec), order='F').T
    data_R = np.asarray(ra_dec_to_xyz(ra_R, dec_R), order='F').T

    # convert spherical bins to cartesian bins
    bins_transform = angular_dist_to_euclidean_dist(bins)

    return two_point(data, bins_transform, method=method,
                     data_R = data_R, random_state=rng)



def bootstrap_two_point_angular(ra, dec, bins, ra_rand = [], dec_rand = [], geometry = None,
                                method = "standard", Nbootstraps = 10, random_state = None):
    """Angular two-point correlation function

    A separate function is needed because angular distances are not
    euclidean, and random sampling needs to take into account the
    spherical volume element.

    Parameters
    ----------
    ra : array_like
        input right ascention, shape = (n_samples,)
    dec : array_like
        input declination
    bins : array_like
        bins within which to compute the 2-point correlation.
        shape = Nbins + 1
    method : string
        "standard" or "landy-szalay".
    Nbootstraps : int
        number of bootstrap resamples
    random_state : integer, np.random.RandomState, or None
        specify the random state to use for generating background

    Returns
    -------
    corr : ndarray
        the estimate of the correlation function within each bin
        shape = Nbins
    dcorr : ndarray
        error estimate on dcorr (sample standard deviation of
        bootstrap resamples)
    bootstraps : ndarray
        The full sample of bootstraps used to compute corr and dcorr
    """
    ra = np.asarray(ra)
    dec = np.asarray(dec)
    rng = check_random_state(random_state)

    if method not in ['standard', 'landy-szalay']:
        raise ValueError("method must be 'standard' or 'landy-szalay'")

    if bins.ndim != 1:
        raise ValueError("bins must be a 1D array")

    if (ra.ndim != 1) or (dec.ndim != 1) or (ra.shape != dec.shape):
        raise ValueError('ra and dec must be 1-dimensional '
                         'arrays of the same length')

    data = np.asarray(ra_dec_to_xyz(ra, dec), order='F').T

    # convert spherical bins to cartesian bins
    bins_transform = angular_dist_to_euclidean_dist(bins)

    bootstraps = []

    for i in range(Nbootstraps):
        # time this
        start_time_repeat = time.time() 
        
        # draw a random sample with N points if galaxy positions not already given
        if len(ra_rand) < 1000 * Nbootstraps and len(dec_rand) < 1000 * Nbootstraps:
            ra_R, dec_R = two_pt_angular_corr.create_random_gals(geometry, len(ra))
        # use the inputted random galaxy positions
        elif len(ra_rand) > 1000 * Nbootstraps and len(dec_rand) > 1000 * Nbootstraps:
            ra_R = ra_rand[int(i * len(ra_rand) / Nbootstraps) : int((i + 1) * len(ra_rand) / Nbootstraps)]
            dec_R = dec_rand[int(i * len(dec_rand) / Nbootstraps) : int((i + 1) * len(dec_rand) / Nbootstraps)]
            
        # ra_R, dec_R = uniform_sphere((min(ra), max(ra)),
        #                              (min(dec), max(dec)),
        #                              2 * len(ra))

        data_R = np.asarray(ra_dec_to_xyz(ra_R, dec_R), order='F').T
        
        #print(len(data_R))

        if i > 0:
            # random sample of the data
            data_b = sklearn.utils.resample(data, replace = True)
            #ind = np.random.randint(0, data.shape[0], data.shape[0])
            # print(len(ind))
            # print(len(set(list(ind))))
            #data_b = data[ind]
            # print("len(data_b) = %1.1f" %len(data_b))
        else:
            data_b = data

        bootstraps.append(two_point(data_b, bins_transform, method = method,
                                    data_R = data_R, random_state=rng))
        
        end_time_repeat = time.time()
        print("Bootstrap {} took {} seconds!".format(i+1, np.round(end_time_repeat - start_time_repeat, 2)))
     
    print("Nrandom =", len(ra_R))

    bootstraps = np.asarray(bootstraps)
    corr = np.mean(bootstraps, 0)
    corr_err = np.std(bootstraps, 0, ddof=1)

    return corr, corr_err, bootstraps #{"ra": ra_R, "dec": dec_R}


