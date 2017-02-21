""" Initialize mbxaspy """

from __future__ import print_function

import sys
import os

from defs import *
from utils import *

__all__ = ['sys', 'os']

# version
print("Running mbxaspy with \n" + sys.version)

# initialize mpi environment
para = para_class()
if ismpi():
    try:
        from mpi4py import MPI
        __all__ += ['MPI']
        para = para_class(MPI)
    except ImportError:
        raise ImportError(' Please install or load mpi4py first. Running with one processor. ')
else:
    print(" Not an MPI environment. Running the serial version. ")

if __name__ == '__main__':
    print(__file__ + ": initialization module for mbxaspy")
