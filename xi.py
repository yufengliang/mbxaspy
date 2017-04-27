""" Compute the xi matrix based on i and f scf """

from __future__ import print_function


import sys
import os

from constants import *
from utils import *
from init import *


def compute_full_sij(fproj):
    """
    Given the final proj_class, calculate the full S_ij matrix
    for each kind of atom

    ground-state atom: sij is just qij
    excited sij stored under proj of fscf
    """
    fproj.full_sij = []
    for kind in range(fproj.nspecies):
        full_sij = sp.matrix(sp.zeros([fproj.nprojs[kind], fproj.nprojs[kind]]))
        # Extract sij for this kind of atom
        if fproj.xs == kind:
            # if it is an excited atom
            sij = fproj.sij
        else:
            sij = fproj.qij[kind]
        # fproj.para.print(sij) # debug
        l_offset1 = il1 = 0
        for l1 in fproj.l[kind]:
            l_offset2 = il2 = 0
            for l2 in fproj.l[kind]:
                if l2 == l1:
                    full_sij[l_offset1 : l_offset1 + 2 * l1 + 1, l_offset2 : l_offset2 + 2 * l2 + 1] = sp.identity(2 * l1 + 1) * sij[il1][il2]
                l_offset2 += 2 * l2 + 1
                il2 += 1
            l_offset1 += 2 * l1 + 1
            il1 += 1
        # fproj.para.print(full_sij) # debug
        fproj.full_sij.append(full_sij)


def compute_xi(iscf, fscf):
    """ 
    compute the xi matrix using two given scfs
    """
    userin  = user_input

    if userin.scf_type == 'shirley_xas':

        # The pseudo part: xi_{mn}^PS = < nk | B_j > < B_j | ~ B_i > < ~ B_i | ~mk >
        xi      = iscf.obf.eigvec.H * fscf.obf.overlap * fscf.obf.eigvec # All 3 must be sp.matrix

        # PAW corrections:
        # \sum_{I, l, l'} < nk | beta_Il > S_ll' < ~beta_Il' | ~mk >
        proj_offset = 0
        proj = fscf.proj

        if userin.do_paw_correction:
            # loop over species first, then loop over atoms within each species
            # just to respect the order as in Qespresso/shirley_xas
            for I in range(proj.natom):
                atom_name   = proj.atomic_pos[I][0]
                s           = proj.ind[atom_name]
                nprojs      = proj.nprojs[s]
                full_sij    = proj.full_sij[s]
                proj_range  = slice(proj_offset, proj_offset + nprojs)
                xi          += iscf.proj.beta_nk[proj_range, :].H \
                            * full_sij \
                            * fscf.proj.beta_nk[proj_range, :]
                proj_offset += nprojs

        xi      = xi.transpose() # fscf.nbnd x iscf.nbnd. This is important !

        return xi

    return None

def plot_xi(xi):
    """
    plot a heap map for a complex matrix xi
    Take the abs of each element and plot
    
    sp: scipy or numpy
    """
    heatmap = plt.pcolor(abs(sp.array(xi)), cmap = 'seismic')
    plt.axis([0, xi.shape[1], 0, xi.shape[0]])
    plt.axes().set_aspect('equal')
    plt.savefig('test_xi.eps', format = 'eps', dpi = 1000)
    plt.close()

def eig_analysis_xi(xi, postfix = ''):
    """
    Analyze the eigenvalues of the transformation matrix xi

    sp: scipy or numpy
    la: linalg
    """
    size = min(xi.shape)
    xi_eigval, xi_eigvec = la.eig(xi[0 : size, 0 : size])
    # Now I plot the abs of eigenvalues
    out_eigval = sorted(abs(xi_eigval), reverse = True)
    plt.stem(out_eigval)
    plt.xlim([-1, size])
    #plt.savefig('test_xi_eig.eps', format = 'eps', dpi = 1000)
    plt.savefig('test_xi_eig{0}.png'.format(postfix), format = 'png')
    plt.close()
    sp.savetxt('test_xi_eig{0}.dat'.format(postfix), out_eigval, delimiter = ' ')
    # Analyze eigenvalues and return a message
    msg = ''
    return msg


def compute_xi_c(xi, xmat_c, nocc, nbnd_i_use = nbnd_max):
    """
    Compute
        sum_c xi_{i, c} < phi_c | O | phi_h > ^ *
    =   sum_c < phi_c | ~phi_i > < phi_h | O | phi_c >
    =   (sum_c < ~phi_i | phi_c > < phi_c | O | phi_h >)^* 
    for all final-state orbital i.

    c sums over all EMPTY initial-state orbitals

    arguments:
    xi:     xi matrix that is already calculated
    xmat_c: *1d* array, < phi_c | O | phi_h > for user-chosen O with c ranging from LUMO to nbnd_i
    nocc: effective occupation number
    """
    xmat_c_ = xmat_c.copy()
    if nocc % 1 > small_thr: # partial occupation
        xmat_c_[0] *= nocc % 1
    nocc = int(nocc)
    nbnd_i = min(xi.shape[1], len(xmat_c_), nbnd_i_use)
    xmat_c_ = sp.matrix(xmat_c_)
    if xmat_c_.shape[1] > 1: xmat_c_ = xmat_c_.H
    return xi[:, nocc : nbnd_i] * xmat_c_[nocc : nbnd_i, 0]

    

