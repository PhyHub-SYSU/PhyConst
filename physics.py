#!/usr/bin/env python
# coding: utf-8

# In[1]:


# import numpy as np
# import scipy
# from scipy.optimize import fsolve, leastsq
# from scipy.interpolate import interp1d
# from scipy.fft import fft
# from scipy.optimize import curve_fit

import math

# In[2]:


"""
Basic physical constant
"""
c = 299792458 # speed of light
h = 6.62607015e-34 # Plank's constant
hbar = h/(2 * math.pi) # Planck constant, reduced
e = 1.602176634e-19 # charge of a proton
hbar_c = hbar * c # conversion constant
hbar_c_square = (hbar * c) ** 2 # conversion constant


eV = e # electronvolt
hc = h*c
MeV_per_c_square = 1e6 * eV / c**2 # M eV/c^2
GeV_per_c_square = 1e9 * eV / c**2


u =  1.66053906660e-27 # unified atomic mass unit
m_e = 9.1093837015e-31 # electron mass
m_p = 1.007276466622 * u # proton mass
m_n =  1.00866491595 * u # neutron mass
m_D = 1875.61294258 * MeV_per_c_square # deuteron mass 
m_alpha = 6.644657e-27 # alpha

mu_0 = 4 * math.pi * 1e-7 # permeability of free space
epsilon_0 = 1/(mu_0 * c **2) # permittivity of free space

e_square_per_4pi_epsilon_0 = e**2/(4 * math.pi * epsilon_0)

alpha = e_square_per_4pi_epsilon_0/hbar_c # fine-structure constant
r_e = e_square_per_4pi_epsilon_0 /(m_e * c **2) # classical electron radius
lambdabar_e = hbar/(m_e * c) # electron Canpton wavelength / 2π
a_infty = r_e / alpha **2
R_infty = m_e*(e_square_per_4pi_epsilon_0)**2 / (2 * hc * hbar**2) # Rydberg constant
sigma_T = 0.66524587321 * 1e-28 # Thomson cross section

mu_B = e * hbar / (2 * m_e) # Bohr magneton
mu_N = e * hbar / (2 * m_p) # nuclear magneton
omega_cycl_e_per_B = e / m_e # electron cyclotron freq./field
omega_cycl_p_per_B = e / m_e # proton cyclotron freq./field

G_N = 6.67430e-11 # gravitational constant
g_N = 9.80665 # standard gravitational accel.


N_A = 6.02214076e23 # Avogadro constant
k = 1.380649e-23 # Boltzmann's constant
V_m = 22.41396954 * 1e-3 # molar volume, ideal gas at STP
b = 2.897771955e-3 # Wien displacement law constant
sigma = math.pi**2 * k **4 / (60 * hbar **3 * c**2) # Stefan-Boltzmann constant 

G_F_per_cubic_hbar_c = 1.1663788e-5 * (1e9 * eV) **(-2) # Fermi coupling constant
G_F = G_F_per_cubic_hbar_c * hbar_c **3 # Fermi coupling constant * (\hbar c)^3
square_sin_hat_theta = 0.23121 # weak-mixing angle
m_W = 80.379 * GeV_per_c_square # W± boson mass
m_Z = 91.1876 * GeV_per_c_square # Z^0 boson mass 
alpha_s = 0.1179 # strong coupling constant

"""
Useful data of some particle 
"""

m_pi = 139.570 * MeV_per_c_square # π^±
m_pi_neu = 134.977 * MeV_per_c_square # π^0
m_K = 493.7 * MeV_per_c_square # K^±
m_Sigma_neg = 1197.4 * MeV_per_c_square # \Sigma^-
m_Lambda =  115.68 * MeV_per_c_square # \Lambda^0
m_rho = 775.11 * MeV_per_c_square # ρ^±
m_rho_neu = 775.26 * MeV_per_c_square # ρ^0
m_sigma = 500 * MeV_per_c_square # f_0(500)
m_omega = 782.65 * MeV_per_c_square # omega 

tau_Lambda = 2.6e-10 # Mean lifetime of \Lambda^0


m_H = 1.00784 * u
R_H = R_infty*(1-m_e/m_p)

"""
Other constant
"""

E_0 = hc * R_infty # Hydrogen-like ionization energy
atm = 101325 # atmospheric pressure


# In[3]:


"""
Astronomical constant
"""
cm = 1e-2
g = 1e-3

M_P = math.sqrt(hbar_c / G_N) # Planck mass
l_P = math.sqrt(hbar * G_N / c **3) # Planck length

ys = 31556925.1 # tropical year (equinox to equinox, 2020)
ys_sidereal = 31558149.8 # sidereal year (period of Earth around Sun relative to stars) 

au = 149597870700 # astronomical unit
pc = 3.08567758149e16 # parsec
ly = 0.946073e16 # light year

square_deg = (math.pi / 180) **2 # solid angle

M_sun = 1.98841e30
R_Schwarzschild_sun = 2 * G_N * M_sun / c **2 # Schwarzschild radius of the Sun
R_sun = 6.597e8 # nominal Solar equatorial radius
S_sun = 1361 # nominal Solar constant
T_sun = 5772 # nominal Solar photosphere temperature 
L_sun = 3.828e26 # nominal Solar luminosity 
R_Schwarzschild_earth = 8.870055940e-3 # Schwarzschild radius of the Earth
M_earth = 5.97217e24 # earth mass
R_earth = 6.3781e6 # nominal Earth equatorial radius

M_Ch_mu = 3.097972 * M_P **3 / m_H **2  # Chandrasekhar mass * average molecular weight per electron
L_Ed_per_M =  1.2570651798e31 / M_sun # Eddington luminosity / mass of the central object

Jy = 1e-26 # jansky (flux density)
# f_0 = # luminosity conversion 

# TODO


# In[8]:


def gamma(v):
    return 1/math.sqrt(1 - v**2 / c**2)


# In[4]:


def scale(L):
    return math.sqrt(L*(L+1))

def g_J(S,L,J):
    return 1+1/2*(scale(J)**2-scale(L)**2+scale(S)**2)/scale(J)**2


# In[5]:


def bind_energy(t, p, n):
    return (n*m_n + p*m_p - t*u) * c**2/(p+n)

def nuclear_mass(DeltaE_per_keV, Z, A):
    return (A-Z)*m_n + Z*m_p - A*DeltaE_per_keV*1e3*eV/c**2


# In[6]:


# def FFT(func, x, N=300):
#    func_fft = math.abs(fft(func)[0:N])
#    k = math.linspace(0, 1/( (x.max() - x.min()) / x.shape[0]), N)
#    return (k, func_fft)



