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


class user_input(object):
    """ input and record user-defined argument from stdin """

    def __init__(self):
        # user-defined variables and default values
        self.ipath = '.'
        self.fpath = '.'
        self.nbnd = 0
        self.nelect = 0

    def read(self):
        """ input from stdin """
        lines = sys.stdin.read()
        var_dict = input_arguments(lines)
    
__all__ = ['para_class', 'user_input']
