""" mbxaspy main """

from __future__ import print_function

import sys
import os

from utils import *
from defs import *
from init import *

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

# loop over spin and kpoints
for isk in range(para.pool.nsk):

    # Import the initial-state scf calculations
    para.sep_line()
    para.print(' Import initial-state scf for (ispin, ik) = ({0},{1}) \n'.format(para.pool.sk_list[isk][0], para.pool.sk_list[isk][1]))
    iscf.input(isk = isk)

    # Check the final-state scf calculations
    para.sep_line()
    para.print(' Import final-state scf for (ispin, ik) = ({0},{1}) \n'.format(para.pool.sk_list[isk][0], para.pool.sk_list[isk][1]))
    fscf.input(is_initial = False, isk = isk)

    # Compute \xi
    xi = compute_xi(iscf, fscf)
    para.print(xi) # debug
    plot_xi(xi, sp) # debug
    msg = eig_analysis_xi(xi, sp, la) # debug



