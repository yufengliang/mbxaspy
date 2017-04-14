""" Compute spectra """

from __future__ import print_function


import sys
import os
import math
from bisect import bisect_left, bisect_right

from constants import *
from utils import *
from init import *


# *** you don't need to use sp here ! Just import init !!!
def gaussian(x, sigma, mu):
    return 1.0 / sigma / sp.sqrt(2 * sp.pi) * sp.exp(- 0.5 * ((x - mu) / sigma) ** 2)

def gaussian_slice(x, sigma, mu):
    """ 
    return a slice of x in which the gaussian is significant 
    exp(-0.5 * ((x - mu) / sigma) ** 2) < given_threshold
    """
    r = sp.sqrt(-2.0 * sigma * sp.log(small_thr))
    x_lo = bisect_left(x, mu - r)
    x_hi = bisect_right(x, mu + r)
    return slice(x_lo, x_hi)

def lorentzian(x, sigma, mu):
    return 1.0 / sp.pi * 0.5 * sigma / ((x - mu) ** 2 + (0.5 * sigma) ** 2)


def stick_to_spectrum(stick, spec_info, smear_func = gaussian):
    """
    Given a stick list: [energy, intensity] * #stick, produce
    the corresponding spectrum

    spec_info: information of the spectrum (energy range, broadening ...). 
    You may use a user_input_class for this
    
    Assume a uniform energy grid (which may not be ideal)
    """
    ener_axis = sp.linspace(spec_info.ELOW, spec_info.EHIGH, spec_info.NENER + 1)
    spectrum = sp.zeros(len(ener_axis))
    for i in range(len(stick)):
        xslice = slice(0, len(ener_axis))
        # Narrow down the range if doing guassian
        if smear_func == gaussian:
            xslice = gaussian_slice(ener_axis, spec_info.SIGMA, stick[i][0])
        spectrum[xslice] += stick[i][1] * smear_func(ener_axis[xslice], spec_info.SIGMA, stick[i][0])
    return ener_axis, spectrum

def Af_to_stick(Af):
    """
    Given the final-state amplitudy Af, return a stick list
    Af: a dictionary of array([energy, amplitude])
    """
    return [ [Af[conf][0].real, float(abs(Af[conf][1]) ** 2)] for conf in Af]
    
def eff_nocc(nelec, nspin, ispin):
    """ 
    Find the number of effective electrons to indicate the Fermi level

    example:
            spin =      1               2
    nelec
    13(odd)            6.5    (up = 7, down = 6)
    12(even)            6     (up = 6, down = 6)
    """
    if nspin == 1:
        return nelec / 2.0
    else:
        return nelec / 2 + ispin * (nelec % 2)


def spectrum0(scf, ixyz, nocc = 0, smearing = 'gauss'):
    """
    Convert the xmat of scf into a non-interacting spectrum
    """
    # *** multiple core levels unsupported now (ncp > 1)
    empty = slice(int(nocc), scf.nbnd)
    size = scf.nbnd - int(nocc)
    # print(empty, size, scf.eigval[empty].shape, scf.xmat.shape) # debug
    # stick: energy v.s. intensity
    stick = sp.concatenate((scf.eigval[empty].reshape(size, 1), 
            (abs(scf.xmat[empty, 0, ixyz]) ** 2).reshape(size, 1)),
            axis = 1)
    # partial occupancy of LUMO
    if nocc % 1 > small_thr:
        stick[0, 1] *= nocc % 1 # adjust intensity
    smear_func = gaussian
    if smearing == 'lor': smear_func = lorentzian
    return stick_to_spectrum(stick, scf.userin, smear_func)


def convolute_spec(spec, spec_xps):
    """
    Convolute some spectra with a given XPS spectrum

    continuous version:

    A_total (E) = int dE' A_xps (E - E') A(E')

    discrete version:
    A_{total, i} = sum^i_{j = j_lo} A_{i - j} B_j
    
    Assume they have the same energy axes
    """
    nener = min(spec.shape[0], spec_xps.shape[0])
    spec_ = spec.copy()
    for ie in range(nener):
        # e.g., c3 = a0 b3 + a1 b2 + a2 b1 + a3 b0
        spec_[ie, 1 :: ] = sp.matrix(spec_xps[:ie, 1]).T * sp.matrix(spec[ie :: -1, 1 ::])
    return spec_
