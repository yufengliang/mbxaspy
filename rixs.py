"""
Calculate RIXS spectra using the determinant formalism
"""

from __future__ import print_function


import sys
import os
import bisect
import numbers

from constants import *
from utils import *
from init import *
from spectra import *

class rixs_map_class(object):
    """
    a class for storing and manipulating the RIXS 2D map
    """
    def __init__(self, rixs_map = None):
        if isinstance(rixs_map, rixs_map_class):
            self.I = sp.matrix(sp.zeros(rixs_map.I.shape))
            # need improvements ***
            self.omega_in = rixs_map.omega_in.copy()
            self.omega_out = rixs_map.omega_out.copy()
            self.loss_mode = rixs_map.loss_mode
            return
        # other constructions ***

    def __imul__(self, other):
        if isinstance(other, numbers.Number):
            self.I *= other
            return self

    def __iadd__(self, other):
        if isinstance(other, rixs_map_class):
            # There also should be a checker for energy axis ***
            self.I += other.I
            return self

    def mp_sum(self, comm = None):
        if valid_comm(pool.rootcomm):
            self.I = comm.allreduce(self.I, op = MPI.SUM)

    def savetxt(self, fname):
        sp.savetxt(fname, self.I, delimiter = ' ', fmt = '%.6e')

    def imshow(self, fname):
        if self.loss_mode:
            extent = [min(self.omega_in), max(self.omega_in), min(self.omega_out), max(self.omega_out)]
        else:
            extent = [min(self.omega_out), max(self.omega_out), min(self.omega_in), max(self.omega_in)]        
        plt.imshow(self.I, extent = extent, origin = 'lower')
        #plt.imshow(self.I)
        plt.axes().set_aspect('equal')
        plt.savefig(fname, format = 'png', bbox_inches = 'tight', dpi = 100)
        plt.close()      

def rixs_f1(xi, nelec, xmat_in, xmat_out, ener_i, ener_f,
            nbnd_i = -1, nbnd_f = -1, 
            e_in_lo = -2.0, e_in_hi = 8.0, nener_in = 100,
            loss_mode = False, eloss_range = -10, nener_out = 100,
            Gamma_h = 0.2, Gamma_f = 0.2,
            I_thr = 1e-3, 
            comm = None):
    """
    
    Calculate the RIXS at the f(1)-level using brute-force enumeration
    This is just an ad-hoc solution.

    xi:         given xi matrix (nbnd_f * nbnd_i)
    nelec:      number of electrons
    xmat_in:    matrix element for the incoming photon (1d array) with specific polarization
    xmat_out:   matrix element for the outgoing photon (1d array) with specific polarization
    ener_i:     initial-state orbital energies
    ener_f:     final-state orbital energies

    e_in_lo:    lower energy bound of the incoming photon (eV)
    e_in_hi:    higher ***
    nener_in:   #energy points
    loss_mode:  use out-going photon energy (False) or energy loss (True)
    e_out_lo:   lower energy bound of the outgoing photon (eV)
    e_out_hi:   higher ***
    nener_out:  #energy points

    Gamma_h:    core-hole broadening (eV)
    Gamma_f:    final-state broadening (eV)

    I_thr:      intensity filter (useful for plotting spectra)

    comm:       mpi communicator
    
    """

    # in case they will be specified by users
    if nbnd_i < 0: nbnd_i = xi.shape[1]
    if nbnd_f < 0: nbnd_f = xi.shape[0]
    
    # make sure xmat_in and xmat_out are 1d row vectors
    xmat_in = sp.matrix(xmat_in).reshape(len(xmat_in), 1)
    xmat_out = sp.matrix(xmat_out).reshape(len(xmat_out), 1)

    # sum up all initial-state channels after absorbing the incoming photon
    xi_in = sp.matrix(xi[:, nelec : nbnd_i]) * sp.matrix(xmat_in[nelec : nbnd_i, 0])
    
    # sum up amplitudes for all emission channels
    xi_out = sp.matrix(xi[:, nelec : nbnd_i]) * sp.matrix(xmat_out[nelec : nbnd_i, 0])    


    ## auxiliary matrix for absorption
    A_mat = sp.matrix(sp.zeros((nbnd_f, nelec + 1), dtype = sp.complex128))
    # initialize the common part
    A_mat[ : nbnd_f, : nelec] = xi[ : nbnd_f, : nelec]
    A_mat[ : nbnd_f, -1] = xi_in[: nbnd_f, 0]
    A_det = la.det(A_mat[: nelec + 1, : nelec + 1])
    A_inv = la.inv(A_mat[: nelec + 1, : nelec + 1])
    # get the expansion coefficient as for the absorption case
    A_zeta = A_mat[nelec : nbnd_f, :] * A_inv

    ## auxiliary matrices for emission
    # emission from conduction bands
    xi_c_mat = sp.matrix(sp.zeros((nelec + 1, nbnd_i), dtype = sp.complex128))
    xi_c_mat[ : nelec, : nelec] = xi[ : nelec, : nelec]
    # xi_c_mat[ : nelec, -1] = xi_out[ : nelec, 0]
    # emission from valence bands
    xi_v_mat = sp.matrix(sp.zeros((nelec + 1, nbnd_i), dtype = sp.complex128))
    xi_v_mat[ : nelec, : nbnd_i] = xi[ : nelec, : nbnd_i]

    ## RIXS matrix elements - depending on the omega_in
    
    # energy axis
    omega_in = sp.linspace(e_in_lo, e_in_hi, nener_in + 1)

    # distribute jobs over the omega_in axis
    if valid_comm(comm):
        rank, size = comm.Get_rank(), comm.Get_size()
    else:
        rank, size = 0, 1
    
    iw_local = range(rank, nener_in + 1, size) # indices of omega_in that will be processed in this rank

    # M_v1c1 (omega_in): there are nener_in frequencies (excluding the point at e_hi_in)
    # *** Be careful of memory issue here ***
    Mv1c1 = [sp.matrix(sp.zeros((nelec, nbnd_i - nelec), dtype = sp.complex128)) for iw in range(len(iw_local))]
    
    # compute Mv1c1
    for c1p in range(nelec, nbnd_f):

        # transition amplitude to the final state c1p: just a complex number
        Ac1p = A_zeta[c1p - nelec, -1] * A_det

        # update the auxiliary matrices with the new c1p row

        xi_c_mat[ -1, : nbnd_i] = xi[c1p, : nbnd_i]
        xi_c_square_mat = sp.concatenate((xi_c_mat[:, : nelec], xi_out[: nelec + 1, 0]), axis = 1)
        xi_c_square_mat[-1, -1] = xi_out[c1p, 0]
        xi_v_mat[ -1, : nbnd_i] = xi[c1p, : nbnd_i]

        # I am gonna save the sherman-morrison formula for later
        xi_c_det = la.det(xi_c_square_mat)
        xi_v_det = la.det(xi_v_mat[:, : nelec + 1])
        xi_c_inv = la.inv(xi_c_square_mat)
        xi_v_inv = la.inv(xi_v_mat[:, : nelec + 1])
        xi_c_zeta = xi_c_inv * xi_c_mat[:, nelec : nbnd_i]
        xi_v_zeta = xi_v_inv * xi_v_mat[:, nelec : nbnd_i]

        # emission amplitude: indexing: (v1, c1)
        Ev1c1 = sp.matrix(sp.zeros((nelec, nbnd_i - nelec), dtype = sp.complex128))
        Ev1c1 += xi_c_zeta[: nelec, :] * xi_c_det # emission from conduction bands
        Ev1c1 -= la.kron(xmat_out[: nelec, 0], xi_v_zeta[nelec, : nbnd_i - nelec]) # emission from valence bands

        for iw, w_ind in enumerate(iw_local):
            # Now this is the RIXS matrix element
            Mv1c1[iw] += Ev1c1 * Ac1p / (omega_in[w_ind] - ener_f[c1p] + 1j * Gamma_h) 


    ## Construct the RIXS map for the given range of omega_in

    rixs_intensity = sp.matrix(sp.zeros((nener_in + 1, nener_out + 1), dtype = sp.float64))

    class spec_info_class: pass
    spec_info = spec_info_class()
    if loss_mode:
        e_out_lo, e_out_hi = -5.0, eloss_range  
    else:
        e_out_lo, e_out_hi = e_in_lo - eloss_range, e_in_hi
    spec_info.ELOW, spec_info.EHIGH, spec_info.NENER, spec_info.SIGMA = e_out_lo, e_out_hi, nener_out, Gamma_f

    # omega_out energy axis
    omega_out = sp.linspace(e_out_lo, e_out_hi, nener_out + 1)

    # find the brightest transition
    Mv1c1_max = 0
    for iw in range(len(iw_local)):
        Mv1c1_max = max(Mv1c1_max, abs(Mv1c1[iw]).max())
    if size > 1:
        Mv1c1_max = comm.allreduce(sp.array([Mv1c1_max]), op = MPI.MAX)[0]

    # threshold for plotting
    M_thr = Mv1c1_max * sp.sqrt(I_thr)

    # broaden the spectral value for a fixed omega_in into a spectrum of omega_out
    for iw, w_ind in enumerate(iw_local):
        # look for significant matrix elements
        coords = sp.where(abs(Mv1c1[iw]) > M_thr)
        if loss_mode:
            stick = [[ener_i[c1 + nelec] - ener_i[v1], abs(Mv1c1[iw][v1, c1]) ** 2] for v1, c1 in zip(coords[0], coords[1])]
        else:
            stick = [[omega_in[iw] - (ener_i[c1 + nelec] - ener_i[v1]), abs(Mv1c1[iw][v1, c1]) ** 2] for v1, c1 in zip(coords[0], coords[1])]
        omega_out, rixs_intensity[w_ind, :] = stick_to_spectrum(stick, spec_info, smear_func = gaussian)
    
    # if more than one MPI task
    if size > 1:
        rixs_intensity = comm.allreduce(rixs_intensity, op = MPI.SUM)

    # store intensity and axes under rixs_map
    rixs_map = rixs_map_class()
    rixs_map.I = rixs_intensity
    rixs_map.omega_in = omega_in
    rixs_map.omega_out = omega_out
    rixs_map.loss_mode = loss_mode

    return rixs_map

if __name__ == "main":
    """
    Below is for test purpose
    """
    pass
