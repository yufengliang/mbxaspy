""" Analyze the origins of major transitions """

from __future__ import print_function


import sys
import os
import bisect

from constants import *
from utils import *
from init import *

def xatom(proj, xmat):
    """
    Given < beta | nk >  and < nk | r | phi_c >, obtain
    sum p = rx, ry, rz
    | sum_{nk} < beta | nk > < nk | x | phi_c > | ^ 2

    This should reflect which atom is excited.
    """
    beta_nk = proj.beta_nk
    nbnd = min(beta_nk.shape[1], xmat.shape[0])
    beta_c = sp.array([0.0] * proj.nproj)
    for ixyz in range(3):
        # why so ugly ?! *** 
        beta_c += sp.array( abs( sp.matrix(beta_nk[:, : nbnd]) * sp.matrix(xmat[: nbnd, 0, ixyz]).T ) ) [:, 0] ** 2
    atom_proj = [0.0] * proj.natom
    for p in range(len(beta_c)):
        atom_proj[proj.iproj2atom[p]] += beta_c[p]
    return max(range(proj.natom), key = lambda x : atom_proj[x])

