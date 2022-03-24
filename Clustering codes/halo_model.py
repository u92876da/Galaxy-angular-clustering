
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 14:38:52 2021

@author: u92876da
"""

import config as cf
import HOD

@cf.dataclass
class halo_model: # no need for this to be a class, can be a data struct instead
    
    delta_c0: float = 1.68647 # / D(z) # overdensity of spherical collapse at z = 0
    #sigma_M_z0_interp = 0 # interpolated 2d spline of sigma8, saved in a file somewhere
    
    @staticmethod
    def calc_halo_mass_func(M, z, dM_frac = 0.1, author_year = "ST99"): # dn(M,z)/dM
        
        # from Jenkins, 2001; Sheth-Tormen 1999
        
        M = M * cf.u.solMass / cf.u.littleh # M in M_solar / h
        # rho_m0 = cf.cosmo.rho_m(0) * cf.u.solMass * (cf.u.littleh ** 2) * (cf.u.kpc ** (-3)) # check this!
        # R = ((3 * M) / (4 * cf.np.pi * rho_m0)) ** (1/3)
        # R = R.to((cf.u.Mpc / cf.u.littleh), cf.u.with_H0(cf.cosmo.H0 * cf.u.km / (cf.u.s * cf.u.Mpc))) # unit conversion
        # dlnR_dM = 1 / (3 * M)
        # print(R)
        # dlnsigma_dM = dlnR_dM * cf.cosmo.sigma(cf.np.round(R.value,2), z, derivative = True, kmin = 0, kmax = cf.np.inf)
        # sigma = cf.cosmo.sigma(R.value, z, kmin = 0, kmax = cf.np.inf)
        # #print(sigma)
        # f = halo_model.halo_mass_func_f(sigma, author_year)
        
        sigma_M_z = halo_model.calc_sigma_M_z(M.value, z)
        sigma_M_dM_z = halo_model.calc_sigma_M_z((1 + dM_frac) * M.value, z)
        dlnsigma_dM = (cf.np.log(sigma_M_dM_z) - cf.np.log(sigma_M_z)) / (dM_frac * M)
        f = halo_model.halo_mass_func_f(halo_model.calc_sigma_M_z(M.value, z), author_year)
        
        # check this, makes no difference for b_h as it falls out in normalization!
        rho_m_z = cf.cosmo.rho_m(z) * cf.u.solMass * (cf.u.littleh ** 2) * (cf.u.kpc ** (-3))

        dn_dM = - dlnsigma_dM * f * rho_m_z / M
        dn_dM = dn_dM.to((cf.u.littleh ** 4 / (cf.u.solMass * cf.u.Mpc ** 3)), \
                         cf.u.with_H0(cf.cosmo.H0 * cf.u.km / (cf.u.s * cf.u.Mpc)))
        return dn_dM
    
    @staticmethod
    def halo_mass_func_f(sigma, author_year = "ST99"):
        
        nu = halo_model.delta_c0 / sigma
        
        if author_year == "ST99":
            A = 0.3222
            a = 0.707
            p = 0.3
            f = A * cf.np.sqrt(2 * a / cf.np.pi) * 1 + cf.np.power(1 / (a * nu ** 2), p) \
                * nu * cf.np.exp(-0.5 * a * nu ** 2)
        elif author_year == "PS74":
            f = cf.np.sqrt(2 / cf.np.pi) * nu * cf.np.exp(- 0.5 * nu ** 2)
            
        return f
    
    @staticmethod
    def calc_sigma_M_z(M, z, method = "manual"): # correctly returns input sigma8
        M = M * cf.u.solMass / cf.u.littleh # M in M_solar / h
        # calculate interpolated 2d spline of sigma8, mass
        # from Jenkins 2001
        #cls.sigma_M_z0_interp = 0
        rho_m0 = cf.cosmo.rho_m(0) * cf.u.solMass * (cf.u.littleh ** 2) * (cf.u.kpc ** (-3)) # check this!
        R = ((3 * M) / (4 * cf.np.pi * rho_m0)) ** (1/3) # factor of 200?
        R = R.to((cf.u.Mpc / cf.u.littleh), cf.u.with_H0(cf.cosmo.H0 * cf.u.km / (cf.u.s * cf.u.Mpc))) # unit conversion
        print(R)
        
        if method == "manual":
            def sigma_squared_integrand(k, R, z):
                P_k = cf.cosmo.matterPowerSpectrum(k, z) # k in h/Mpc comoving
                W_kR = cf.cosmo.filterFunction("tophat", k, R) #3 * (cf.np.sin(k * R) - k * R * cf.np.cos(k * R)) / ((k * R) ** 3)
                return (P_k * (k * W_kR) ** 2) / (2 * cf.np.pi ** 2)
            sigma_squared = cf.scipy.integrate.quad(sigma_squared_integrand, 0, cf.np.inf,\
                                                    args = (R.value, z,), epsrel=1e-6)[0]
            sigma = sigma_squared ** (1 / 2)
        
        elif method == "colossus":
            sigma = cf.cosmo.sigma(R.value, z)
            
        return sigma
    
    @staticmethod
    def calc_halo_bias(M, z, author_year = "SMT01"): # Eulerian bias
        
        sigma_M_z = halo_model.calc_sigma_M_z(M, z)
        delta_c_z = halo_model.delta_c0 # evolve delta_c0 over z, sufficiently weak to ignore
        nu = delta_c_z / sigma_M_z
        
        if author_year == "SMT01": # Sheth, Mo and Tormen 2001
            a = 0.707
            b = 0.5
            c = 0.6
            term1 = cf.np.sqrt(a) * delta_c_z
            term2 = (a ** (3/2)) * (nu ** 2)
            term3 = cf.np.sqrt(a) * b * ((a * (nu ** 2)) ** (1 - c))
            term4 = (a * (nu ** 2)) ** c
            term5 = b * (1 - c) * (1 - (c / 2))
            b_h = 1 + (1 / term1) * (term2 + term3 - (term4 / (term4 + term5)))
        
        if author_year == "MW96": # Mo and White 1996
            b_h = 1 + ((nu ** 2) - 1) / delta_c_z
        return b_h
    
    # put galaxy method in bias class
    @staticmethod
    def calc_halo_number_density(z, M_min, author_year):
        
        def number_density_integrand(u):
            return halo_model.calc_halo_mass_func(u, z, author_year = author_year).value
            
        number_density = cf.scipy.integrate.quad(number_density_integrand, M_min, cf.np.inf)[0]
        return number_density
    
    @staticmethod
    def calc_gal_number_density(z, M_min, author_year = "ST99", HOD_model = None, HOD_params = None):
        
        def number_density_integrand(u):
            return (halo_model.calc_halo_mass_func(u, z, author_year = author_year) * HOD_model(u, HOD_params)).value
            
        number_density = cf.scipy.integrate.quad(number_density_integrand, M_min, cf.np.inf)[0]
        return number_density
    
    # put galaxy method in HOD class
    @staticmethod
    def calc_average_halo_bias(z, M_min, b_author_year = "SMT01", f_author_year = "ST99"):
        
        def average_bias_numerator_integrand(u):
            return (halo_model.calc_halo_mass_func(u, z, author_year = f_author_year) \
                    * halo_model.calc_halo_bias(u, z, author_year = b_author_year)).value
            
        average_bias_numerator = cf.scipy.integrate.quad(average_bias_numerator_integrand, M_min, cf.np.inf)[0]
        number_density = halo_model.calc_number_density(z, M_min, b_author_year)
        average_bias = average_bias_numerator / number_density
        print(average_bias)
        return average_bias
    
    def calc_average_gal_bias(z, M_min, b_author_year = "SMT01", f_author_year = "ST99", HOD_model = None, HOD_params = None):
        def average_bias_numerator_integrand(u):
            return (halo_model.calc_halo_mass_func(u, z, author_year = f_author_year) \
                    * halo_model.calc_halo_bias(u, z, author_year = b_author_year) * HOD_model(u, HOD_params)).value

# -----------------------------------------------------------------------------
# halomod example functions

def plot_halo_mass_function(hm):
    cf.plt.plot(hm.m, hm.dndm)
    cf.plt.xscale('log')
    cf.plt.yscale('log')
    
    cf.plt.xlabel("Halo Mass [$h^{-1} M_\odot$]")
    cf.plt.ylabel(r"dn/dm [$h^2 M_\odot^{-1} {\rm Mpc}^{-3}$]");
    cf.plt.show()

def plot_galaxy_power_spectrum(hm):
    cf.plt.plot(hm.k_hm, hm.power_auto_tracer, label='Galaxy-Galaxy Power')
    cf.plt.plot(hm.k_hm, hm.power_1h_auto_tracer, ls='--', label='1-halo term')
    cf.plt.plot(hm.k_hm, hm.power_2h_auto_tracer, ls='--', label='2-halo term')
    
    cf.plt.xscale('log')
    cf.plt.yscale('log')
    
    cf.plt.ylim(1e-5,1e6)
    cf.plt.legend()
    cf.plt.xlabel("Wavenumber [h/Mpc]")
    cf.plt.ylabel(r"Galaxy Power Spectrum [${\rm Mpc^3} h^{-3}$]");
    cf.plt.show()
    
def plot_halo_profile(hm, halo_masses, r = cf.np.logspace(-3, 1, 20)):
    for m in halo_masses:
        cf.plt.plot(r, hm.halo_profile.rho(r=r, m=m), label=f'm={m:1.2e}')
    
    cf.plt.legend()
    cf.plt.yscale('log')
    cf.plt.xscale('log')
    
    cf.plt.xlabel("Distance from Centre [Mpc/h]")
    cf.plt.ylabel(r"Halo Density [$h^2 M_\odot {\rm Mpc}^{-3}$]")
    cf.plt.show()
    
def plot_correlation_function(hm):
    cf.plt.plot(hm.r, hm.corr_auto_tracer, label='Tinker at z=0')

    cf.plt.xscale('log')
    cf.plt.yscale('log')
    
    cf.plt.xlabel("r [Mpc/h]")
    cf.plt.ylabel("Correlation Function")
    cf.plt.show()
    
def plot_HOD(hm):
    #hm.hod_model = hod
    print(hm.hod_params)
    cf.plt.plot(hm.m, hm.total_occupation) #, label=hod)

    cf.plt.xscale('log')
    cf.plt.yscale('log')
    cf.plt.ylim(1e-4, 1e4)
    cf.plt.xlim(1e9, 1e16)
    #cf.plt.legend()
    
    cf.plt.xlabel("m [$h^{-1}M_\odot$]")
    cf.plt.ylabel("N(m)")
    cf.plt.show()
    
def plot_two_pt_angular_corr(hm):
    
    cf.plt.plot((hm.theta * cf.u.rad).to(cf.u.deg).value, hm.angular_corr_gal)
    cf.plt.xscale('log')
    cf.plt.yscale('log')
    
# -----------------------------------------------------------------------------

def main():
    hm = cf.TracerHaloModel()
    two_pt = cf.AngularCF(z = 2.0, zmin = 1.7, zmax = 2.3, theta_min = (1e-3 * cf.u.deg).to(cf.u.rad).value,\
                          theta_max = (1e-1 * cf.u.deg).to(cf.u.rad).value, theta_log = True)
    plot_two_pt_angular_corr(two_pt)
    print(two_pt.hod_params)
    #cf.inspect.getmembers(hm, predicate=cf.inspect.ismethod)
    #print(dir(two_pt))
    #print([x for x, y in hm.__dict__.items() if type(y) == cf.FunctionType])

if __name__ == "__main__":
    #main()
    redshift = 0
    rho_m0 = cf.cosmo.rho_m(0) * cf.u.solMass * (cf.u.littleh ** 2) * (cf.u.kpc ** (-3))
    rho_m0 = rho_m0.to(cf.u.solMass * (cf.u.littleh ** 2) * (cf.u.Mpc ** (-3)), cf.u.with_H0(cf.cosmo.H0 * cf.u.km / (cf.u.s * cf.u.Mpc)))
    
    print(halo_model.calc_sigma_M_z((8 ** 3) * 4 * cf.np.pi * rho_m0.value / 3, 0))
    #halo_model.calc_average_halo_bias(2, 1e13, b_author_year = "MW96", f_author_year = "PS74")
        