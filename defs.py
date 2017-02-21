""" Definitions of classes """

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
        

