""" mbxaspy main """

from __future__ import print_function

import sys
import os

from utils import *
from defs import *
from init import *

user_input.read()

para.print(' Reading initial-state scf from: \n ' + user_input.ipath)
iscf.input(user_input = user_input,
           path = user_input.ipath, 
           mol_name = user_input.mol_name_i,
           is_initial = True)
