""" Compute BSE spectra based on the DSCF core-hole potential """

from __future__ import print_function

from constants import *

from init import *
from spectra import *

def bse(xi, iscf, fscf, nocc, ixyz_list, offset = 0.0, evec = None):
    """
    Construct the BSE Hamiltonian with only empty initial-state orbitals

    H_{ii'} = < i | f > < f | H | f' > < f' | i' >

    i, i' are empty-orbital indices.

    < f | H | f' > = delta_{ff'} e_f

    < i | f > is the *transpose* of xi

    """
    
    # BSE Hamiltonian
    H_BSE = xi[:, int(nocc) : ].T * sp.diag(fscf.eigval[ : fscf.nbnd_use]) * xi[:, int(nocc) : ].conjugate()
    #H_BSE = xi[:, int(nocc) : ].H * sp.diag(fscf.eigval[ : fscf.nbnd_use]) * xi[:, int(nocc) : ]

    # Solve BSE
    E, Af = la.eig(H_BSE)

    ind = sp.array(sorted(range(len(E)), key = lambda i : E[i]))
    # Make sticks
    sticks = sp.zeros((len(E), len(ixyz_list) + 2))
    sticks[:, 0] = E[ind].real + offset

    for i, ixyz in enumerate(ixyz_list):

        ixmat = sp.array([ xmat_ixyz( iscf.xmat[ib, 0, :], ixyz, evec = evec ) for ib in range(int(nocc), iscf.nbnd_use) ])
        
        # transition amplitude: A_op = sum_i Af*_i < i | O | h >
        A_op = sp.dot(ixmat, sp.matrix(Af).conjugate())[0, ind]
        #A_op = sp.dot(ixmat, sp.matrix(Af))
        
        sticks[:, i + 2] = sp.square(abs(A_op))

    return [list(stick) for stick in list(sticks)]

        
