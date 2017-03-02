""" Definitions of classes """

from __future__ import print_function


import sys
import os
import struct

from constants import *
from utils import *
from io_mod import *

class pool_class(object):
    """ 
    pool class: define parallelization over k-points

    Each pool has a given number of procs (nproc_per_pool) and 
    takes care of ONE k-point of ONE spin at a time.

    """

    def set_pool(self, nproc_per_pool = 1):
        """ 
        set up pools so that each pool has at least nproc_per_pool procs 
        
        examples:
        para.size       = 10    # total number of processors
        nproc_per_pool  = 3     # least number of processors per pool
        10 / 3 = 3 pools

        pool    size    proc
           0       4    0, 3, 6, 9   
           1       3    1, 4, 7
           2       3    2, 5, 8
        """
        para = self.para
        if nproc_per_pool > 0:
            self.npp    = nproc_per_pool
            if para.size < self.npp:
                para.print(' Insufficient number of procs for nproc_per_pool = ' + str(self.npp))
                para.print(' Reduce nproc_per_pool to ' + str(para.size))
                self.npp = para.size
            self.n      = para.size / self.npp  # number of pools
            self.i      = para.rank % self.n    # the index of the pool
            self.rank   = para.rank / self.n    # rank within the pool
            # actual pool size (npp plus residue)
            self.size   = self.npp + (self.i < para.size % self.npp)
    
    def set_sklist(self, nspin = 1, nk = 1):
        """ 
        set up a list of spin and kpoint tuples to be processed on this proc

        example:
        nspin = 2, nk = 5
        self.n  =   3   # number of pools
        2 * 5 / 3 = 3   # each pool at least deal with 3 tuples

        pool    tuple
           0    (0, 0), (0, 3), (1, 1), (1, 4)
           1    (0, 1), (0, 4), (1, 2)
           2    (0, 2), (1, 0), (1, 3)
        """
        self.sklist = []
        for s in range(nspin):
            for k in range(nk):
                if (s * nk + k) % self.n == self.i:
                    self.sklist.append((s, k))

        
    def __init__(self, para):
        self.para = para
        self.set_pool(nproc_per_pool = 1)

        
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

        # initialize pool for k-points
        self.pool = pool_class(self)

    def print(self, msg = '', rank = 0, flush = False):
        """ print at given rank """
        if self.rank == rank:
            print(msg) # recursively ?
        if flush:
            sys.stdout.flush()


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
        self.nelec      = 0
        self.gamma_only = False
        self.scf_type   = 'shirley_xas'
        self.xas_arg    = 5
        self.mol_name_i = 'mol_name'
        self.mol_name_f = 'xatom'
        self.nproc_per_pool = 1

    def read(self):
        """ input from stdin or userin"""
        if isanaconda():
            try:
                userin = open(sys.argv[1], 'r')
            except IOError:
                para.print(" Can't open user-defined input " + sys.arg[1] + " Halt. ", flush = True)
                para.stop()
        else: userin = sys.stdin
        lines = userin.read()
        var_input = input_arguments(lines)
        for var in set(vars(self)) & set(var_input): # This can be improved
            try:
                # convert var into correct data type as implied in __init__ and set attributes
                setattr(self, var, convert_val(var_input[var], type(getattr(self, var))))
            except:
                pass
        self.ipath = os.path.abspath(self.ipath)
        self.fpath = os.path.abspath(self.fpath)
        userin.close()

class kpoints_class(object):
    """ store information related to kpoints """

    def __init__(self, nk = 1):
        # variable list and default values
        self.nk         = nk
        # in future there may be a list of k-vectors (not needed now)

    def set_kpool(self, para):
        """ distribute k-points over pools """
        pass
        

class optimal_basis_set_class(object):
    """ store information related to shirley optimal basis functions """
    

    def __init__(self, nbasis = 0, nbnd = 0):
        sp = self.sp
        # variable list and default values
        self.nbasis     = nbasis
        self.nbnd       = nbnd
        self.eigval     = sp.array([])    # eigenvalues (band energies)
        self.eigvec     = sp.array([])    # eigenvectors (wavefunctions)


    # Have you considered k-points and spins ?!
    def input_eigval(self, fh):
        sp = self.sp
        try:
            self.eigval = input_from_binary(fh, 'double', self.nbnd, 0)
        except struct.error:
            pass
        self.eigval = sp.array(self.eigval)
        self.para.print(self.eigval) # debug

    def input_eigvec(self, fh):
        pass
        

class paw_class(object):
    """ store information related to projector-argumentated wave (PAW) method """


    def __init__(self):
        # variable list and default values
        self.natom      = 0             # number of PAW atoms


class scf_class(object):
    """ pipeline data from self-consistent-field calculations """

    def __init__(self):
        sp = self.sp
        # variable list and default values
        self.nbnd   = 0                         # number of bands    
        self.nk     = 0                         # number of k-points
        self.nelec  = 0                         # number of electrons
        self.ncp    = 1                         # number of core levels
        self.nspin  = 1                         # number of spins
        self.xmat   = sp.array([])         # single-particle matrix elements


    def input_shirley(self, user_input, path, mol_name, is_initial):
        """ input from shirley xas """
        para = self.para

        # construct file names
        xas_prefix = mol_name + '.xas'
        xas_data_prefix = xas_prefix + '.' + str(user_input.xas_arg)
        ftype = ['info', 'eigval', 'eigvec', 'proj']

        if is_initial: ftype += ['xmat'] # need pos matrix element for the initial state
        
        for f in ftype:
            fname = os.path.abspath(path + '/' + xas_data_prefix + '.' + f)
            if f != 'info': binary = '' # Does it matter if not use b ?
            else: binary = 'b'
            try:
                fh = open(fname, 'r' + binary)
            except:
                para.print(" Can't open " + fname + '. Check if shirley_xas finishes properly. Halt. ', flush = True)
                para.stop()
            if f == 'info':
                lines = fh.read()
                var_input = input_arguments(lines, lower = True)
                # para.print(var_input) # debug
                for var in ['nbnd', 'nk', 'nbasis', 'nelec', 'ncp', 'nspin']:
                    if var in var_input:
                        try:
                        # convert var into correct data type as implied in __init__ and set attributes
                            setattr(self, var, convert_val(var_input[var], type(getattr(self, var))))
                        except:
                            pass
                    else:
                        para.print(' Variable "' + var + '" missed in ' + fname, flush = True)
                        para.stop()

                # initialize k-points
                if user_input.gamma_only: self.nk = 1
                self.kpt = kpoints_class(nk = self.nk)

                # get spin and k-point index processed by this pool
                para.pool.set_sklist(nspin = self.nspin, nk = self.nk)

                # initialize obf
                self.obf = optimal_basis_set_class(nbnd = self.nbnd) # there're more: nbasis, ...

            if f == 'eigval':
                self.obf.input_eigval(fh)
            if f == 'eigvec':
                        pass
                        
        
    def input(self, user_input, path, mol_name, is_initial):
        """ input from one shirley run """
        para = self.para
        if user_input.scf_type == 'shirley_xas':
            para.print(' Import wavefunctions and energies from shirley_xas calculation .. ')
            self.input_shirley(user_input, path, mol_name, is_initial)
        else:
            para.print(' Unsupported scf input: ' + scf_type + ' Halt. ')
            para.stop()
            


__all__ = [c for c in dir() if c.endswith('_class')]


if __name__ == '__main__':
    print(__file__ + ": definitions of classes for mbxaspy")

    # test user_input
    #user_input = user_input_class()
    #user_input.read()
    #print(vars(user_input))

    # test para and pool
    para = para_class()
    para.size = 50
    for r in range(para.size):
        para.rank = r
        para.pool.set_pool(para, 10)
        print(para.pool.i)
