""" Compute spectra """

from __future__ import print_function


import sys
import os
import math
from bisect import bisect_left, bisect_right
import numbers

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

## ========================================================================================================================
# rewritten with OOP

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

## defs related to sticks

# a stick array is:
# [[energy, "info", os_1, os_2, os_3, ...], ...]
# info contains information of the stick
# os's are oscillator strengths
# *** Maybe it will be become an object one day.

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

def Af_to_sticks(Af, offset = 0.0):
    """
    Given the final-state amplitudy Af, return a stick array:
    sticks: [[energy, "info", os_1, os_2, os_3, ...], ...]
    
    Af: a dictionary of array([energy, amplitude])
    """
    return [ [ complex(Af[conf][0]).real + offset, conf, float(abs(Af[conf][1]) ** 2) ] for conf in Af]

def calc_occ_pdos(scf, ixyz_list, nocc = 0, evec = None):
    """
    Calculate the integral of the occupied part of the spectra
    Q = int_{E < E_f} d E sigma(E)
    """
    # *** multiple core levels unsupported now (ncp > 1)
    return [sum([abs(xmat_ixyz(scf.xmat[ib, 0, :], ixyz, evec)) ** 2 for ib in range(int(nocc))]) for ixyz in ixyz_list]

def os_sum(sticks):
    """
    sum up the oscillator strengths of all sticks
    """
    return sp.sum(sp.array([s[2 : ] for s in sticks]), axis = 0)

def sticks_filter(sticks, os_min = 1e-1, maxn = 10):
    pass

##

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
    # I[: I1.shape[0], : I1.shape[1]] += I1 # left aligned
    I[: I1.shape[0], col - I1.shape[1] :] += I1 # right aligned
    # I[: I2.shape[0], : I2.shape[1]] += I2 # left aligned
    I[: I2.shape[0], col - I2.shape[1] :] += I2 # right aligned
    return I


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

        ncol = len(sticks[0]) - 2 if len(sticks) > 0 else 1 # 1 column for no stick
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
            #  1 2   5   -\   1 2 5 
            #  3 4   6   -/   3 4 6
            self.I = sp.concatenate((self.I, new_spec), axis = 1) if self.I is not None else new_spec
        elif mode == 'additive':
            #  1 2   5   -\   1 7 
            #  3 4   6   -/   3 10
            self.I = add_I(self.I, new_spec)
        else:
            raise ValueError('spec_class().add_sticks: unknown mode {}'.format(mode))

        self.ncol = self.I.shape[1]

        #if save_sticks:
        #    # filter sticks
        #    new_sticks = sticks_filter(sticks, sticks_thr, max_sticks)
        #    self.sticks_all += new_sticks

    def savetxt(self, fname, offset = 0.0):
         sp.savetxt(fname, sp.concatenate( ( sp.matrix(self.ener_axis + offset).T, self.I ), axis = 1), 
                    delimiter = ' ', fmt = '%.6e')

    def __add__(self, other):
        """ define spec1 + spec2 """
        if not same_axis(self, other):
            raise IndexError('cannot add spectra with different energy axes.')
        spec = spec_class(ener_axis = self.ener_axis)
        spec.I = add_I(self.I, other.I)
        spec.ncol = max(self.ncol, other.ncol)
        return spec

    def __iadd__(self, other):
        return self + other

    def __or__(self, other):
        """ 
        define spec1 | spec2
        Alternate columns of spec1.I and spec2.I
        """
        if not same_axis(self, other):
            raise IndexError('cannot alternate spectra with different energy axes.')
        if self.ncol != other.ncol:
            raise IndexError('cannot alternate spectra with different no. of cols ({0}, {1}). '.format(self.ncol, other.ncol))
        spec = spec_class(ener_axis = self.ener_axis)
        spec.I = sp.zeros((self.lener, self.ncol * 2))
        spec.I[:, 0 :: 2] = self.I
        spec.I[:, 1 :: 2] = other.I
        spec.ncol = spec.I.shape[1]
        return spec

    def __mul__(self, other):
        """
        define spec1 * spec2

        usage 1: multiple intensity I by a constant

        usage 2: spectral convolution

        spec(E) = integrate dE' spec1(E - E') spec2(E')

        Let's zero_ind = 4 for I1 and I2, then the convoluted intensity at ind = 7 is:

        I[7] = ( ... + I1[8] * I2[3] +  I1[7] * I2[4] + I1[6] * I2[5] + ... ) dE
        
        index relations:
        ie + zero_ind = ie1 + ie2

        0 <= ie1, ie2 < lener, then
        0 <= ie + zero_ind - ie1 < lener, so
        ie + zero_ind - lener < ie1 <= ie + zero_ind, so
        
        max(ie + zero_ind - lener + 1, 0) <= ie1 < min(ie + zero_ind + 1, lener)

        usage 3: convolute a spec_class with a stick array

        spec(E) = \sum_j I_j spec1(E - E_j) 

        I_j's are stick heights and E_j's stick energies.

        I_j's are expected to be dimensionless.
        
        """
    
        # usage 1: a spec_class and a constant
        if isinstance(other, numbers.Number):
            spec = spec_class(ener_axis = self.ener_axis)
            spec.I = self.I * other
            return spec
        # usage 3: a spec_class and a stick array
        if isinstance(other, list):
            if len(other) > 0 and len(other[0]) > 2: # then it is a stick array
                spec = spec_class(ener_axis = self.ener_axis)
                spec.I = sp.zeros(self.I.shape)
                spec.ncol = self.ncol
                for stick in other: # *** to be parallelized
                    e_ind = int(stick[0] / self.dE)
                    l = stick[0] % self.dE
                    if e_ind > 0:
                        spec.I[e_ind + 1 : , :] += ( l * self.I[: -e_ind - 1, :] + (1 - l) * self.I[1 : -e_ind, :] ) * stick[2]
                    else:
                        spec.I[: e_ind - 1, :] += ( (1 - l) * self.I[-e_ind : -1, :] + l * self.I[-e_ind + 1 :, :] ) * stick[2]
                return spec
        if not isinstance(other, spec_class):
            raise TypeError('unsupported operand type(s) for *')
        # usage 2
        if not same_axis(self, other):
            raise IndexError('cannot convolute spectra with different energy axes.')
        if abs(self.ener_axis[self.zero_ind]) > small_thr:
            raise IndexError('Energy axis does not contain E=0. Cannot convolute spectra.')
        spec = spec_class(ener_axis = self.ener_axis)
        spec.I = sp.zeros(self.I.shape)
        spec.ncol = self.ncol
        for ie in range(self.lener):
            e1_lo = max(ie + other.zero_ind + 1 - self.lener, 0)
            e1_hi = min(ie + other.zero_ind + 1, self.lener)
            e2_lo = ie + other.zero_ind - e1_hi
            e2_hi = ie + other.zero_ind - e1_lo 
            # only use the 1st col of the second spec
            spec.I[ie, :] = sp.dot(other.I[range(e2_hi, e2_lo, -1), 0],
                                    self.I[range(e1_lo, e1_hi,  1), :]) * other.dE
        return spec

    def __imul__(self, other):
        """ is this meaningful  ? """
        if isinstance(other, numbers.Number):
            self.I *= other
            return self
        return self * other

    def mp_sum(self, comm = None):
        """ sum and bcast I within given comm """
        if comm and comm != MPI.COMM_NULL:
            self.I = comm.allreduce(self.I, op = MPI.SUM) # check if allreduce works for older mpi4py

    def average(self, cols = [], target_col = 0):
        """
        Take the average of the given cols and put it into target_col
        """
        self.I[:, target_col] = sp.average(self.I[:, cols], axis = 1)

    def os_sum(self):
        """
        return the spectral area. 
        Use this with caution if the energy windonw is not wide enough.
        """
        return sp.sum(self.I, axis = 0) * self.dE

def afi(xi, iscf, fscf, nocc, ixyz_list, offset = 0.0, evec = None):
    """
    Calculate the afi (final-initial projection) spectrum defined as:

    sum_f | sum_{i in empty} < f~ | i > < i | r | h > | ^ 2 delta(E - e~_f)

    The amplitude 

    sum_{i in empty} < f~ | i > < i | r | h >

    is calculated as:

    < f~ | r | h > - sum_{i in occ} < f~ | i > < i | r | h >

    to make connections with the final-state rule:

    < f~ | r | h > 

    """
    
    nocc = int(nocc) # do this for now ***
    nbnd_f = min(xi.shape[0], fscf.xmat.shape[0], len(fscf.eigval))

    ixmat = sp.matrix([[xmat_ixyz(iscf.xmat[ib, 0, :], ixyz, evec) for ixyz in ixyz_list] 
            for ib in range(nocc)])
    # sum_{i in occ} < f~ | i > < i | r | h >
    occ_proj = sp.matrix(xi[nocc : nbnd_f, : nocc]).conjugate() * ixmat[:, :]
    #  < f~ | r | h >
    fxmat = sp.matrix([[xmat_ixyz(fscf.xmat[ib, 0, :], ixyz, evec) for ixyz in ixyz_list] 
            for ib in range(nocc : nbnd_f)])
    sticks = fxmat - occ_proj
    sticks = abs(sp.array(sticks)) ** 2
    sticks = sp.column_stack(fscf.eigval[nocc : nbnd_f] + offset, sticks)
    return [list[stick] for stick in sticks]
    
    
    

    
