""" mbxaspy main """

from __future__ import print_function

import sys
import os

from utils import *
from defs import *
from init import *

# input user-defined arguments from stdin
user_input.read()

# import initial-state scf calculations
para.print(' Checking initial-state scf from: \n ' + user_input.path_i + '\n')
iscf.input(user_input = user_input,
           is_initial = True,
           isk = -1)

# import final-state scf calculations
para.print(' Checking final-state scf from: \n ' + user_input.path_f + '\n')
fscf.input(user_input = user_input,
           is_initial = False,
           isk = -1)

# import final-state scf calculations
pass

# check consistency between i and f scf
pass


