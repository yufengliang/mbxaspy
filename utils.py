""" Utilities """

from __future__ import print_function

import sys
import os

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

