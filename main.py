""" mbxaspy main """

from __future__ import print_function

import sys
import os

from utils import *
from defs import *
from init import *
from spectra import *

# input user-defined arguments from stdin
user_input.read()

# Check initial and final state and perform sanity checks
# Check the initial-state scf calculations
para.sep_line()
para.print(' Checking initial-state scf from: \n ' + user_input.path_i + '\n')
iscf.input(isk = -1)

# Check the final-state scf calculations
para.sep_line()
para.print(' Checking final-state scf from: \n ' + user_input.path_f + '\n')
fscf.input(is_initial = False, isk = -1)

# Input < B_i | \tilde{B}_j >
fscf.obf.input_overlap(user_input.path_f, iscf.nbnd, fscf.nbnd)

from xi import *
# Compute full atomic overlap sij
if user_input.scf_type == 'shirley_xas': 
    compute_full_sij(fscf.proj)

nspin = iscf.nspin
# loop over spin and kpoints
for isk in range(para.pool.nsk):

    ispin, ik  = para.pool.sk_list[isk] # acquire current spin
    # weight the sticks according to the k-grid
    weight =  iscf.kpt.weight[isk]; 
    # para.print('weight = {0}'.format(weight)) # debug

    # Import the initial-state scf calculation
    para.sep_line()
    para.print(' Import initial-state scf for (ispin, ik) = ({0},{1}) \n'.format(ispin, ik))
    iscf.input(isk = isk)

    # Import the final-state scf calculation
    para.sep_line()
    para.print(' Import final-state scf for (ispin, ik) = ({0},{1}) \n'.format(ispin, ik), flush = True)
    fscf.input(is_initial = False, isk = isk)

    # Obtain the effective occupation number: respect the initial-state #electrons
    nocc = eff_nocc(iscf.nelec, nspin, ispin)

    # Compute non-interacting spectra *** should I put it in a def ?
    for ixyz in range(3):
        # *** currently only x, y, z
        ener_axis, spec0 = spectrum0(iscf, ixyz, nocc, user_input.smearing)
        ener_axis += user_input.ESHIFT_FINAL
        # print(ener_axis.shape, spec0.shape) # debug
        if isk == 0 and ixyz == 0:
            # initialize spectra and energy axis
            spec0_i = sp.zeros([len(ener_axis), nspin * 4 + 1]) # ener (total x y z) * nspin
            spec0_i[:, 0] = ener_axis
            if user_input.final_1p: spec0_f = spec0_i.copy()
        col_offset = (1 + ixyz) * iscf.nspin + (ispin + 1)
        spec0 *= weight
        spec0_i[:, col_offset] += spec0
        spec0_i[:, ispin + 1] += spec0 # angular average
        if user_input.final_1p:
            ener_axis, spec0 \
            = spectrum0(fscf, ixyz, nocc, user_input.smearing)
            spec0 *= weight
            spec0_f[:, col_offset] += spec0
            spec0_f[:, ispin + 1] += spec0 # angular average

    if not user_input.spec0_only:

        # Compute the transformation matrix xi
        xi = compute_xi(iscf, fscf)

        #para.print(xi) # debug
        if user_input.xi_analysis and para.isroot() and ik == 0:
            # plot_xi(xi) # debug
            if nspin > 1:
                msg = eig_analysis_xi(xi, '_spin_{0}'.format(ispin)) # debug
            else:
                msg = eig_analysis_xi(xi) # debug

        for ixyz in range(3):
            # Compute xi_c
            xi_c = compute_xi_c(xi, iscf.xmat[:, 0, ixyz], nocc)
            # para.print('xi_c.shape = {0}'.format(str(xi_c.shape))) # debug

# Output Spectra

# intial-state one-body
spec0_i[:, 1 : 1 + nspin] /= 3.0
if ismpi() and para.pool.isroot():
    spec0_i[:, 1 : nspin * 4 + 1] = para.pool.rootcomm.reduce(spec0_i[:, 1 : nspin * 4 + 1], op = MPI.SUM)

if para.isroot(): sp.savetxt(spec0_i_fname, spec0_i, delimiter = ' ')

# final-state one-body
if user_input.final_1p:
    spec0_f[:, 1 : 1 + nspin] /= 3.0
    if ismpi() and para.pool.isroot():
        spec0_f[:, 1 : nspin * 4 + 1] = para.pool.rootcomm.reduce(spec0_f[:, 1 : nspin * 4 + 1], op = MPI.SUM) 
    if para.isroot(): sp.savetxt(spec0_f_fname, spec0_f, delimiter = ' ')




