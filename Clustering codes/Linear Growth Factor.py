# Linear Growth Factor in the matter and dark energy domination epoch

# Imports
import numpy as np
import scipy.integrate

# Initial parameters
Omega_m = 0.25
Omega_Lambda = 1 - Omega_m
Omega_k = 0
z = 1
a = 1/(1+z)

def calculate_E(u):
    E_u = np.power((Omega_m*np.power(u,-3)+Omega_Lambda+Omega_k*np.power(u, -2)), 0.5)
    return(E_u)

def integrand(u):
    integrand = np.power(u*calculate_E(u),-3)
    return(integrand)

# Calculate Linear Growth Factor D(a)
integral, err = scipy.integrate.quad(integrand, 0., a)
linear_growth_factor = integral*5*Omega_m*calculate_E(a)/2
#print(linear_growth_factor)

# Normalise using D(a=1)=1
integral, err = scipy.integrate.quad(integrand, 0., 1.)
normalisation_factor = 1/(integral*5*Omega_m*calculate_E(1.)/2)
#print(normalisation_factor)

# Normalised Linear Growth Factor
norm_linear_growth_factor = normalisation_factor*linear_growth_factor
print(norm_linear_growth_factor)

# f = d(ln(D))/d(ln(a))

