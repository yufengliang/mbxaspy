""" Initialize mbxaspy """

from __future__ import print_function

import sys
import os

from utils import *
from defs import *

__all__ = ['sys', 'os']

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
    print(' Not an MPI environment. Running the serial version. ')
__all__ += ['para']
para.print()

# python version
para.print(' Running mbxaspy with: \n Python' + sys.version)
para.print()


# numpy / scipy version: define sp and la
try:
    import scipy as sp
    from scipy import linalg as la
    para.print(' Using Python scientific computation module: \n SciPy ' + sp.__version__)
except ImportError:
    para.print(' Fail to import SciPy. Turn to import NumPy... ')
    import numpy as sp # use sp for scipy or numpy
    from numpy import linalg as la
    para.print(' Using Python scientific computation module: \n NumPy ' + sp.__version__) 
except ImportError:
    para.print(' Fail to import NumPy. Please install or load it first. Halt. ')
    para.stop()
__all__ += ['sp', 'la']
para.print()


# user input
user_input = user_input_class()
__all__ += ['user_input']


# initial- and final-state scf
scf_class.sp                    = sp    # There could be a better way of doing this
scf_class.para                  = para
optimal_basis_set_class.sp      = sp
optimal_basis_set_class.para    = para
iscf = scf_class()
fscf = scf_class()
__all__ += ['iscf', 'fscf']

if __name__ == '__main__':
    print(__file__ + ': initialization module for mbxaspy')
