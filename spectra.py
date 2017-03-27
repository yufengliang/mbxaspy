""" Compute spectra """

from __future__ import print_function


import sys
import os

from constants import *
from utils import *
import math
from bisect import bisect_left, bisect_right


def gaussian(x, sigma, mu, sp):
    return 1.0 / sigma / sp.sqrt(2 * sp.pi) * sp.exp(- 0.5 * ((x - mu) / sigma) ** 2)

def gaussian_slice(x, sigma, mu, sp):
    """ 
    return a slice of x in which the gaussian is significant 
    exp(-0.5 * ((x - mu) / sigma) ** 2) < given_threshold
    """
    r = sp.sqrt(-2.0 * sigma * sp.log(small_thr))
    x_lo = bisect_left(x, mu - r)
    x_hi = bisect_right(x, mu + r)
    return slice(x_lo, x_hi)

def lorentzian(x, sigma, mu, sp):
    return 1.0 / sp.pi * 0.5 * sigma / ((x - mu) ** 2 + (0.5 * sigma) ** 2)


def stick_to_spectrum(stick, spec_info, smear_func = gaussian):
    """
    Given a stick list: [energy, intensity] * #stick, produce
    the corresponding spectrum

    spec_info: information of the spectrum (energy range, broadening ...). 
    You may use a user_input_class for this
    
    Assume a uniform energy grid (which may not be ideal)
    """
    sp = spec_info.sp
    ener_axis = sp.linspace(spec_info.ELOW, spec_info.EHIGH, spec_info.NENER + 1)
    spectrum = sp.zeros(len(ener_axis))
    for i in range(len(stick)):
        xslice = slice(0, len(ener_axis))
        # Narrow down the range if doing guassian
        if smear_func == gaussian:
            xslice = gaussian_slice(ener_axis, spec_info.SIGMA, stick[i][0], sp)
        spectrum[xslice] += stick[i][1] * smear_func(ener_axis[xslice], spec_info.SIGMA, stick[i][0], sp)
    return ener_axis, spectrum


def eff_nelec(nelec, nspin, ispin):
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


def spectrum0(scf, ixyz, eff_nelec = 0, smearing = 'gauss'):
    """
    Convert the xmat of scf into a non-interacting spectrum
    """
    sp = scf.sp
    # *** multiple core levels unsupported now (ncp > 1)
    empty = slice(int(eff_nelec), scf.nbnd)
    size = scf.nbnd - int(eff_nelec)
    # print(empty, size, scf.eigval[empty].shape, scf.xmat.shape) # debug
    # stick: energy v.s. intensity
    stick = sp.concatenate((scf.eigval[empty].reshape(size, 1), 
            (abs(scf.xmat[empty, 0, ixyz]) ** 2).reshape(size, 1)),
            axis = 1)
    # partial occupancy of LUMO
    if eff_nelec % 1 > small_thr:
        stick[0, 1] *= eff_nelec % 1 # adjust intensity
    smear_func = gaussian
    if smearing == 'lor': smear_func = lorentzian
    return stick_to_spectrum(stick, scf.userin, smear_func)
    
        
