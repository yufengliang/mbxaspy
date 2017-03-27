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

    ispin = para.pool.sk_list[isk][0] # acquire current spin
    # weight the sticks according to the k-grid
    weight =  iscf.kpt.weight[isk]; 
    para.print('weight = {0}'.format(weight))

    # Import the initial-state scf calculations
    para.sep_line()
    para.print(' Import initial-state scf for (ispin, ik) = ({0},{1}) \n'.format(ispin, para.pool.sk_list[isk][1]))
    iscf.input(isk = isk)

    # Check the final-state scf calculations
    para.sep_line()
    para.print(' Import final-state scf for (ispin, ik) = ({0},{1}) \n'.format(ispin, para.pool.sk_list[isk][1]))
    fscf.input(is_initial = False, isk = isk)

    # Compute \xi
    xi = compute_xi(iscf, fscf)
    #para.print(xi) # debug
    if user_input.xi_analysis and para.pool.isroot():
        plot_xi(xi, sp) # debug
        msg = eig_analysis_xi(xi, sp, la) # debug

    # Compute non-interacting spectra *** should I put it in a def ?
    for ixyz in range(3):
        # *** currently only x, y, z
        ener_axis, spec0 = spectrum0(iscf, ixyz, eff_nelec(iscf.nelec, nspin, ispin), user_input.smearing)
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
            = spectrum0(fscf, ixyz, eff_nelec(iscf.nelec, nspin, ispin), user_input.smearing) # yes, electron number from iscf
            spec0 *= weight
            spec0_f[:, col_offset] += spec0
            spec0_f[:, ispin + 1] += spec0 # angular average

# Output Spectra
# mpi_sum spectrum ***
spec0_i[:, 1 : 1 + nspin] /= 3.0
sp.savetxt(spec0_i_fname, spec0_i, delimiter = ' ')
if user_input.final_1p:
    spec0_f[:, 1 : 1 + nspin] /= 3.0
    sp.savetxt(spec0_f_fname, spec0_f, delimiter = ' ')




