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
para.print(' Checking initial-state scf from: \n ' + user_input.path_i + '\n')
iscf.input(user_input = user_input,
           is_initial = True,
           isk = -1)

# Check the final-state scf calculations
para.print(' Checking final-state scf from: \n ' + user_input.path_f + '\n')
fscf.input(user_input = user_input,
           is_initial = False,
           isk = -1)

# loop over spin and kpoints
for isk in range(para.pool.nsk):

    # Import the initial-state scf calculations
    para.print(' Import initial-state scf for (ispin, ik) = ({0},{1}) \n'.format(para.pool.sk_list[isk][0], para.pool.sk_list[isk][1]))
    iscf.input(user_input = user_input,
               is_initial = True,
               isk = isk)

    # Check the final-state scf calculations
    para.print(' Import final-state scf for (ispin, ik) = ({0},{1}) \n'.format(para.pool.sk_list[isk][0], para.pool.sk_list[isk][1]))
    fscf.input(user_input = user_input,
               is_initial = False,
               isk = isk)


