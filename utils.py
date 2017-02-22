""" Utilities """

from __future__ import print_function

import sys
import os
from ast import parse

def is_valid_variable_name(name):
    """test if name is a valid python variable name"""
    try:
        parse('{} = None'.format(name))
        return True
    except (SyntaxError, ValueError, TypeError) as err:
        return False


def ispython3x():
    if sys.version_info.major == 3:
        return True
    else:
        return False


def ismpi():
    mpi_var = {'SLURM_JOB_ID', 'PBS_O_WORKDIR'}
    if set(mpi_var) & set(os.environ):
        return True
    else:
        return False


