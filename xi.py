""" Compute the xi matrix based on i and f scf """

from __future__ import print_function


import sys
import os

from constants import *
from utils import *


def compute_full_sij(fproj):
    """
    Given the final proj_class, calculate the full S_ij matrix
    for each kind of atom

    ground-state atom: sij is just qij
    excited sij stored under proj of fscf
    """
    sp = fproj.sp
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
    sp      = iscf.sp
    para    = iscf.para

    if iscf.userin.scf_type == 'shirley_xas':

        # The pseudo part: xi_{mn}^PS = < nk | B_j > < B_j | ~ B_i > < ~ B_i | ~mk >
        xi      = iscf.obf.eigvec.H * fscf.obf.overlap * fscf.obf.eigvec # All 3 must be sp.matrix

        # PAW corrections:
        # \sum_{I, l, l'} < nk | beta_Il > S_ll' < ~beta_Il' | ~mk >
        proj_offset = 0
        proj = fscf.proj

        do_paw_correction = True # Just a flag for turning on PAW core correction
        if do_paw_correction:
            for I in range(proj.natom):
                atom_name   = proj.atomic_pos[I][0]
                kind        = proj.ind[atom_name]
                nprojs      = proj.nprojs[kind]
                full_sij    = proj.full_sij[kind]
                xi          += iscf.proj.beta_nk[proj_offset : proj_offset + nprojs, :].H \
                            * full_sij \
                            * fscf.proj.beta_nk[proj_offset : proj_offset + nprojs, :]
                proj_offset += nprojs

        xi      = xi.transpose() # fscf.nbnd x iscf.nbnd. This is important !

        return xi

    return None

def plot_xi(xi, sp):
    """
    plot a heap map for a complex matrix xi
    Take the abs of each element and plot
    
    sp: scipy or numpy
    """
    from matplotlib import pyplot as plt # This import is temporarily here ***
    heatmap = plt.pcolor(abs(sp.array(xi)), cmap = 'seismic')
    plt.axis([0, xi.shape[1], 0, xi.shape[0]])
    plt.axes().set_aspect('equal')
    plt.savefig('test_xi.eps', format = 'eps', dpi = 1000)
    plt.close()

def eig_analysis_xi(xi, sp, la):
    """
    Analyze the eigenvalues of the transformation matrix xi

    sp: scipy or numpy
    la: linalg
    """
    size = min(xi.shape)
    xi_eigval, xi_eigvec = la.eig(xi[0 : size, 0 : size])
    from matplotlib import pyplot as plt # This import is temporarily here ***
    # Now I plot the abs of eigenvalues
    plt.stem(abs(xi_eigval))
    plt.xlim([-1, size])
    plt.savefig('test_xi_eig.eps', format = 'eps', dpi = 1000)
    plt.close()
    # Analyze eigenvalues and return a message
    msg = ''
    return msg


