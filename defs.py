""" Definitions of classes """

from __future__ import print_function


import sys


class para_class(object):
    """ para class: wrap up user-defined mpi variables for parallelization """

    def __init__(self, MPI = None):
        if MPI:
            self.comm = MPI.COMM_WORLD
            self.size = self.comm.Get_size()
            self.rank = self.comm.Get_rank()

        else:
            self.comm = None
            self.size = 1
            self.rank = 0

        self.npool = 1


    def print(self, msg = '', rank = 0):
        if self.rank == rank:
            print(msg) # recursively ?


    def stop(self):
        if self.comm:
            self.comm.Abort(0)
        else:
            sys.exit(0)

