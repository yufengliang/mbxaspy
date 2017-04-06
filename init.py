""" Initialize mbxaspy """

from __future__ import print_function

import sys
import os

from utils import *
from defs import *
from para_defs import *

__all__ = ['sys', 'os']

# initialize mpi environment
para = para_class()
if ismpi():
    try:
        from mpi4py import MPI
        __all__ += ['MPI']
        para = para_class(MPI)
    except ImportError:
        raise ImportError(' Please install or load mpi4py first. MPI init fails. Halt. ')
else:
    print(' Not an MPI environment. Running the serial version. ')
# define pools
para.pool = pool_class(para)
__all__ += ['para']
para.print()

# python version
para.print(' Running mbxaspy with: \n Python' + sys.version)
para.print()
if 'anaconda' in sys.version.lower():
    para.print(' Detect Python Anaconda. Standard input is blocked. Read user input from the 1st arg. ')
    if len(sys.argv) < 2:
        para.print(' Please enter the file name of user input as the 1st arg.', flush = True)
        para.stop()

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

# Pass attributes
for obj in ['scf_class', 'optimal_basis_set_class', 'proj_class', 'user_input_class']:
    for attr in ['sp', 'para']:
        setattr(eval(obj), attr, eval(attr))

scf_class.userin = user_input

# initial- and final-state scf
iscf = scf_class()
fscf = scf_class()
__all__ += ['iscf', 'fscf']

# flush all the output
para.print('', flush = True)

if __name__ == '__main__':
    print(__file__ + ': initialization module for mbxaspy')
