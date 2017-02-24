""" Definitions of classes """

from __future__ import print_function


import sys
import os

from utils import *
from io_mod import *


class para_class(object):
    """ para class: wrap up user-defined mpi variables for parallelization """

    def __init__(self, MPI = None):
        if MPI is not None:
            self.comm = MPI.COMM_WORLD
            self.size = self.comm.Get_size()
            self.rank = self.comm.Get_rank()

        else:
            self.comm = None
            self.size = 1
            self.rank = 0

        self.npool = 1


    def print(self, msg = '', rank = 0):
        """ print at given rank """
        if self.rank == rank:
            print(msg) # recursively ?


    def stop(self):
        """ stop the code """
        if self.comm is not None:
            # in mpi
            self.comm.Abort(0)
        else:
            # non-mpi
            sys.exit(0)


class user_input_class(object):
    """ input and record user-defined argument from stdin """

    def __init__(self):
        # user-defined variables and default values
        self.ipath      = '.'
        self.fpath      = '.'
        self.nbnd       = 0
        self.nelect     = 0
        self.gamma_only = False
        self.scf_type   = 'shirley_xas'

    def read(self):
        """ input from stdin """
        lines = sys.stdin.read()
        var_input = input_arguments(lines)
        for var in set(vars(self)) & set(var_input): # This can be improved
            try:
                setattr(self, var, convert_val(var_input[var], type(getattr(self, var))))
            except:
                pass


class optimal_basis_set_class(object):
    """"""
    pass

class scf_class(object):
    """ pipeline data from self-consistent-field calculations """


    def __init__(self):
        # variable list and default values
        self.nbnd   = 0                         # number of bands    
        self.nk     = 0                         # number of k-points
        self.nelect = 0                         # number of electrons
        # should I use numpy array ?
        self.eigval = []                        # eigenvalues (band energies)
        self.eigvec = []                        # eigenvectors (wavefunctions)
        self.xmat   = []                        # single-particle matrix elements
        self.obf = optimal_basis_set_class()


    def input_shirley(self):
        pass
        
    def input(self, scf_type):
        """ input from one shirley run """
        if scf_type == 'shirley_xas':
            self.input_shirley()
        else:
            pass
            


__all__ = ['para_class', 'user_input_class']


if __name__ == '__main__':
    print(__file__ + ": definitions of classes for mbxaspy")
    user_input = user_input_class()
    user_input.read()
    print(vars(user_input))

    initial_scf = scf_class()
    initial_scf.input('shirley_xas')
