""" Compute spectra """

from __future__ import print_function


import sys
import os
import math
from bisect import bisect_left, bisect_right

from constants import *
from utils import *
from init import *


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

# This is obsolete and could be incorrect for nspin = 2.   
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
        return nelec / 2 + (1 - ispin) * (nelec % 2)


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
        spec_[ie, 1 :: ] = sp.matrix(spec_xps[:ie + 1, 1]) * sp.matrix(spec[ie :: -1, 1 ::])
    return spec_

# ========================================================================================================================


def xmat_ixyz(xmat, ixyz, evec):
    """
    ixyz > 0:   normal xmat
    ixyz = -1:  angular average of x, y, z
    ixyz = -2:  weight with the electric field direction: evec
    """
    if ixyz >= 0:
        return xmat[ixyz]
    if ixyz == -1:
        return sp.sqrt(sum([abs(xmat[i]) ** 2 for i in range(3)]) / 3.0)
    if ixyz == -2:
        evec_ = sp.array(evec) / la.norm(evec)
        return sum([xmat[i] * evec_[i] for i in range(3)])

def xmat_to_sticks(scf, ixyz_list, nocc = 0, offset = 0.0, evec = None):
    """
    Convert the empty states and their oscillator strengths stored in scf.xmat
    into a stick array

    sticks: [[energy, "info", os_1, os_2, os_3, ...], ...]
    ixyz_list: allowed values of each element: -2, -1, 0, 1, 2, ...
    nocc: effective occupation number
    offset: energy offset in eV
    evec: electric field vector

    """
    # *** multiple core levels unsupported now (ncp > 1)
    sticks = [[scf.eigval[ib] + offset, ''] + 
              [abs(xmat_ixyz(scf.xmat[ib, 0, :], ixyz, evec)) ** 2 for ixyz in ixyz_list] 
              for ib in range(int(nocc), scf.nbnd_use)]
    # partial occupancy of LUMO
    if nocc % 1 > small_thr:
        sticks[0][2 : ] *= nocc % 1 # adjust intensity
    return sticks


def same_axis(self, other):
    if len(self.ener_axis) != len(other.ener_axis) or any(abs(self.ener_axis - other.ener_axis) > zero_thr): return False
    else: return True

def add_I(I1, I2):
    """ 
    Add two 2D arrays.
    Non-overlap region takes the values from one of the two arrays
    """
    if I1 is None and I2 is None: return None
    if I1 is None: return I2
    if I2 is None: return I1
    row, col = max(I1.shape[0], I2.shape[0]), max(I1.shape[1], I2.shape[1])
    I = sp.zeros((row, col))
    I[: I1.shape[0], : I1.shape[1]] += I1
    I[: I2.shape[0], : I2.shape[1]] += I2
    return I

def sticks_filter(sticks, sticks_thr = 1e-1, max_sticks = 10):
    pass


class spec_class(object):
    """
    A class for storing spectra and relevant data
    """
    def __init__(self, spec_info = None, ener_axis = None):
        if spec_info:
            lo, hi = spec_info.ELOW, spec_info.EHIGH # Definitions of ELOW and EHIGH in shirley_xas are slightly different !
            nener = spec_info.NENER
        else:
            lo, hi = -2.0, 8.0 # eV
            nener = 100
        if ener_axis is not None:
            self.dE = (max(ener_axis) - min(ener_axis)) / (len(ener_axis) - 1)
            self.ener_axis = ener_axis.copy()
        else:
            self.dE = (hi - lo) / nener
            self.ener_axis = sp.array([e * self.dE for e in range(int(lo / self.dE) - 1, int(hi / self.dE) + 2)])
        self.lener = len(self.ener_axis)
        self.nener = self.lener - 1
        # indicate where is zero. Ideally, self.lo < 0 and self.hi > 0
        self.zero_ind = min(range(len(self.ener_axis)), key = lambda i : abs(self.ener_axis[i])) 
        self.I = None
        self.ncol = 0
        self.sticks_all = []
            
    def add_sticks(self, sticks, spec_info = None, prefac = 1.0, mode = 'append'):
        #append = True, save_sticks = False, sticks_thr = 1e-1, max_sticks = 10):
        """
        Convert sticks into a spectrum

        Args:
        sticks      : [[energy, "info", os_1, os_2, os_3, ...], ...]
        append      : append the newly obtained spectra to previous results
        save_sticks : whether we save the sticks
        thr_sticks  : threshold for saving sticks
        max_sticks  : maximum no. sticks to save
        """

        ncol = len(sticks[0]) - 2
        new_spec = sp.zeros([self.lener, ncol])
        for i, stick in enumerate(sticks):
            xslice = slice(0, self.lener)
            # Narrow down the range if doing guassian
            if spec_info.smearing == 'gauss':
                smear_func = gaussian
                xslice = gaussian_slice(self.ener_axis, spec_info.SIGMA, stick[0])
            else:
                smear_func = lorentzian
            for col in range(2, 2 + ncol):
                # convolute with the chosen smearing function
                new_spec[xslice, col - 2] += prefac * stick[col] * smear_func(self.ener_axis[xslice], spec_info.SIGMA, stick[0])

        # ensure no. of stick arrays <= no. of spec (self.ncol)
        # save_sticks = save_sticks and self.ncol == len(sticks_all)

        # add new_spec to self.I
        if mode == 'append':
            self.I = sp.concatenate((self.I, new_spec), axis = 1) if self.I is not None else new_spec
        else: # 'additive'
            self.I = add_I(self.I, new_spec)

        self.ncol = self.I.shape[1]

        #if save_sticks:
        #    # filter sticks
        #    new_sticks = sticks_filter(sticks, sticks_thr, max_sticks)
        #    self.sticks_all += new_sticks

    def savetxt(self, fname, offset = 0.0):
         sp.savetxt(fname, sp.concatenate( ( sp.matrix(self.ener_axis + offset).T, self.I ), axis = 1), 
                    delimiter = ' ', fmt = '%.6e')

    def __add__(self, other):
        if not same_axis(self, other):
            raise IndexError('cannot add spectra with different energy axes.')
        spec = spec_class(ener_axis = self.ener_axis)
        spec.I = add_I(self.I, other.I)
        return spec

    def __or__(self, other):
        if not same_axis(self, other):
            raise IndexError('cannot alternate spectra with different energy axes.')
        if self.ncol != other.ncol:
            raise IndexError('cannot alternate spectra with different no. of cols ({0}, {1}). '.format(self.ncol, other.ncol))
        spec = spec_class(ener_axis = self.ener_axis)
        spec.I = sp.zeros((self.lener, self.ncol * 2))
        spec.I[:, 0 :: 2] = self.I
        spec.I[:, 1 :: 2] = other.I
        return spec

    def mp_sum(self, comm = None):
        if comm and comm != MPI.COMM_NULL:
            self.I = comm.allreduce(self.I, op = MPI.SUM) # check if allreduce works for older mpi4py
