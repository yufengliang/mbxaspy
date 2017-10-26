"""
Calculate RIXS spectra using the determinant formalism
"""

from __future__ import print_function


import sys
import os
import bisect

from constants import *
from utils import *
from init import *


def rixs_f1(xi, nelec, xmat_in, xmat_out, ener_i, ener_f, 
            e_in_lo = -2.0, e_in_hi = 8.0, nener_in = 100,
            e_out_lo = -10.0, e_out_hi = 8.0, nener_out = 100,
            Gamma_h = 0.2, Gamma_f = 0.2,
            I_thr = 1e-3, 
            comm = None)
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
    e_out_lo:   lower energy bound of the outgoing photon (eV)
    e_out_hi:   higher ***
    nener_out:  #energy points

    Gamma_h:    core-hole broadening (eV)
    Gamma_f:    final-state broadening (eV)

    I_thr:      intensity filter (useful for plotting spectra)

    comm:       mpi communicator
    
    """

    # in case they will be specified by users
    nbnd_i = xi.shape[1]
    nbnd_f = xi.shape[0]
    
    # make sure xmat_in and xmat_out are 1d row vectors
    xmat_in = sp.matrix(xmat_in).reshape(1, len(xmat_in))
    xmat_out = sp.matrix(xmat_out).reshape(1, len(xmat_out))

    # sum up all initial-state channels after absorbing the incoming photon
    xi_in = sp.matrix(xi[:, nelec : nbnd_i]) * sp.matrix(xmat_in[nelec : nbnd_i]).H
    
    # sum up amplitudes for all emission channels
    xi_out = sp.matrix(xi[:, nelec : nbnd_i]) * sp.matrix(xmat_out[nelec : nbnd_i]).H    


    ## auxiliary matrix for absorption
    A_mat = sp.matrix(sp.zeros((nelec + 1, nelec + 1), dtype = sp.complex128))
    # initialize the common part
    A_mat[ : nelec + 1, : nelelc] = xi[ : nelec + 1, : nelec]
    A_mat[ :, -1] = xi_in
    A_det = la.det(A_det)
    A_inv = la.inv(A_mat)
    # get the expansion coefficient as for the absorption case
    A_zeta = A_mat[nelec : nbnd_f, :] * A_inv


    ## auxiliary matrices for emission
    # emission from conduction bands
    xi_c_mat = sp.matrix(sp.zeros((nelec + 1, nelec + 1), dtype = sp.complex128))
    xi_c_mat[ : nelec, : nelec] = xi[ : nelec, : nelec]
    xi_c_mat[ : nelec, -1] = xi_out[ : nelec, 0]
    # emission from valence bands
    xi_v_mat = sp.matrix(sp.zeros((nelec + 1, nelec + 1), dtype = sp.complex128))
    xi_v_mat[ : nelec, : nelec + 1] = xi[ : nelec, : nelec + 1]

    ## RIXS matrix elements - depending on the omega_in
    
    # energy axis
    omega_in = sp.linspace(e_lo_in, e_hi_in, nener_in + 1)
    
    # M_v1c1 (omega_in): there are nener_in frequencies (excluding the point at e_hi_in)
    # *** Be careful of memory issue here ***
    Mv1c1 = [sp.matrix(sp.zeros((nelec, nbnd_i - nelec), dtype = sp.complex128)) for iw in range(nener_in)]

    # compute Mv1c1
    for c1p in range(nelec, nbnd_f):

        # transition amplitude to the final state c1p: just a complex number
        Ac1p = A_zeta[c1p - nelec, -1] * A_det

        # update the auxiliary matrices with the new c1p row

        xi_c_mat[ -1, : nelec] = xi[c1p, : nelec]
        xi_c_mat[ -1, -1] = xi_out[c1p, 0]
        
        xi_v_mat[ -1, : nelec + 1] = xi[c1p, : nelec + 1]

        # I am gonna save the sherman-morrison formula for later
        xi_c_det = la.det(xi_c_mat)
        xi_v_det = la.det(xi_v_mat)
        xi_c_inv = la.inv(xi_c_mat)
        xi_v_inv = la.inv(xi_v_mat)
        xi_c_zeta = xi_c_inv * xi_c_mat[:, nelec : nbnd_i]
        xi_v_zeta = xi_v_inv * xi_v_mat[:, nelec : nbnd_i]

        # emission amplitude: indexing: (v1, c1)
        Ev1c1 = sp.matrix(sp.zeros((nelec, nbnd_i - nelec), dtype = sp.complex128))
        Ev1c1 += xi_c_zeta[: nelec, :] * xi_c_det # emission from conduction bands
        Ev1c1 -= la.kron(xmat_out[: nelec, 0], xi_v_zeta[nelec, nelec : nbnd_i]) # emission from valence bands

        for iw in range(nener_in):

            Mv1c1[iw] += Ev1c1 * Ac1p / (omega_in[iw] - ener_f[c1p] + 1j * Gamma_h) 


    ## Construct the RIXS map for the given range of omega_in

    rixs_map = sp.matrix(sp.zeros((nener_in, nener_out + 1), dtype = sp.complex128)
    omega_out = sp.linspace(e_lo_out, e_hi_out, nener_out + 1)

    # find the brightest transition
    Mv1c1_max = 0
    for iw in range(nener_in):
        Mv1c1_max = max(Mv1c1_max, abs(Mv1c1[iw]).max())

    # threshold for plotting
    M_thr = Mv1c1_max * sp.sqrt(I_thr)

    # broaden the spectral value for a fixed omega_in into a spectrum of omega_out
    for iw in range(nener_in):
        # look for coordinates
        pass
        
    
        
