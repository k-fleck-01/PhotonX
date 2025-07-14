import numpy as np
from scipy.constants import physical_constants
from scipy.special import spence

###############################################################################
### Physical constants
###############################################################################
ALPHA = physical_constants['fine-structure constant'][0]
RADIUS_E = physical_constants['classical electron radius'][0]
COMPTON_E = physical_constants['Compton wavelength'][0]

###############################################################################
### Spence's function ∫ln(1-z)/z from scipy convention
###############################################################################
def phi(z):
    return - spence(1.0 - z)

###############################################################################
### Auxilliary functions
###############################################################################
def a(r, s):
    return np.real((1.0 - (r + s) / (r * s))**0.5)

def b(r):
    return (1.0 - 1.0 / r)**0.5

###############################################################################
### Mandelstam kinematic integrals
### Defined in De Tollis 1964 and Karplus and Neumann 1951
###############################################################################
def B(r):
    b_real = np.real(-1.0 + 0.5 * b(r) * np.log((b(r) + 1.0) / (b(r) - 1.0)))
    if r.real >= 1.0:
        return b_real - 0.5j * np.pi * b(r)
    else:
        return b_real + 0j

def T(r):
    t_real = np.real((0.5 * np.log((b(r) + 1.0) / (b(r) - 1.0)))**2.0)
    if r.real >= 1.0:
        return t_real - 1.0j * np.pi * np.arccosh(r**0.5)
    else:
        return t_real + 0j

def F(r, s):
    f_real = np.real(0.5 / a(r, s) * (np.log(r * (a(r, s)**2.0 - b(r)**2.0))
        * np.log((a(r, s) + 1.0) / (a(r, s) - 1.0)) 
        + phi((a(r, s) + 1.0) / (a(r, s) + b(r)))
        - phi((a(r, s) - 1.0) / (a(r, s) + b(r)))
        + phi((a(r, s) + 1.0) / (a(r, s) - b(r)))
        - phi((a(r, s) - 1.0) / (a(r, s) - b(r)))))
    if r.real >= 1.0:
        return f_real + 0.5j * np.pi / a(r, s) * np.log((a(r, s) - b(r)) / (a(r, s) + b(r)))
    else:
        return f_real + 0j

def I(r, s):
    return F(r, s) + F(s, r)

###############################################################################
### Polarization resolved matrix elements
###############################################################################
def Mpppp(r, s, t):
    return 1.0 + (2.0 + 4.0 * s / r) * B(s) + (2.0 + 4.0 * t / r) * B(t) \
        + (2.0 * (s**2.0 + t**2.0) / r**2.0 - 2.0 / r) * (T(s) + T(t)) \
        + (-1.0 / s + 0.5 / (r * s)) * I(r, s) \
        + (-1.0 / t + 0.5 / (r * t)) * I(r, t) \
        + (- 2.0 * (s**2.0 + t**2.0) / r**2.0 + 4.0 / r + 1.0 / s
            + 1.0 / t + 0.5 / (s * t)) * I(s, t)

def Mpppn(r, s, t):
    return -1.0 + (-1.0 / r - 1.0 / s - 1.0 / t) * (T(r) + T(s) + T(t)) \
        + (1.0 / t + 0.5 / (r * s)) * I(r, s) \
        + (1.0 / s + 0.5 / (r * t)) * I(r, t) \
        + (1.0 / r + 0.5 / (s * t)) * I(s, t)

def Mppnn(r, s ,t):
    return -1.0 + 0.5 / (r * s) * I(r, s) + 0.5 / (r * t) * I(r, t) \
        + 0.5 / (s * t) * I(s, t)

###############################################################################
### Squared matrix elements (real-valued)
###############################################################################
def amplitudes_polarized(r, s, t):
    return {
        "++++" : np.abs(Mpppp(r, s, t))**2,
        "+++-" : np.abs(Mpppn(r, s, t))**2,
        "++--" : np.abs(Mppnn(r, s, t))**2,
    }

def amplitude(r, s, t):
    return (0.5 * (Mpppp(r, s, t) * Mpppp(r, s, t).conjugate()
                + Mppnn(r, s, t) * Mppnn(r, s, t).conjugate()
                + Mpppp(t, s, r) * Mpppp(t, s, r).conjugate()
                + Mpppp(s, r, t) * Mpppp(s, r, t).conjugate()
                + 4.0 * Mpppn(r, s, t) * Mpppn(r, s, t).conjugate())).real

###############################################################################
### Differential cross section
###############################################################################
@np.vectorize
def diff_cross(k, theta):
    r = complex(k**2.0)
    s = complex(-k**2.0 * (1.0 - np.cos(theta)) / 2.0)
    t = complex(-k**2.0 * (1.0 + np.cos(theta)) / 2.0)
    return ALPHA**4.0 * amplitude(r, s, t) / (k**2.0 * 4.0 * np.pi**2.0)

def diff_cross_polarized(k, theta):
    k = np.asarray(k)
    theta = np.asarray(theta)

    shape = np.broadcast(k, theta).shape
    result = {key: np.zeros(shape) for key in ["++++", "+++-", "++--"]}
    
    ### Flatten arrays for vectorization
    k_flat = k.flatten()
    teta_flat = theta.flatten()

    for idx in range(k_flat.size):
        r = complex(k_flat[idx]**2.0)
        s = complex(-k_flat[idx]**2.0 * (1.0 - np.cos(teta_flat[idx])) / 2.0)
        t = complex(-k_flat[idx]**2.0 * (1.0 + np.cos(teta_flat[idx])) / 2.0)
        squared_amplitudes = amplitudes_polarized(r, s, t)
        prefractor = ALPHA**4.0 / (k_flat[idx]**2.0 * 4.0 * np.pi**2.0)

        result["++++"].flat[idx] = prefractor * squared_amplitudes["++++"]
        result["+++-"].flat[idx] = prefractor * squared_amplitudes["+++-"]
        result["++--"].flat[idx] = prefractor * squared_amplitudes["++--"]

    return {key: value.reshape(shape) for key, value in result.items()}


