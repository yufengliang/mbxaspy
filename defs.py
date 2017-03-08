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
        comm = para.comm
        if nproc_per_pool > 0:
            self.npp    = int(nproc_per_pool)
            if para.size < self.npp:
                para.print(' Insufficient number of procs for nproc_per_pool = ' + str(self.npp))
                para.print(' Reduce nproc_per_pool to ' + str(para.size))
                self.npp = para.size
            self.n      = int(para.size / self.npp)  # number of pools
            self.i      = int(para.rank % self.n)    # the index of the pool
            if comm:
                # set up pool communicators
                self.comm     = comm.Split(para.rank % self.n, para.rank) # intrapool comm
                self.rank     = self.comm.Get_rank()      # rank within the pool
                self.roots    = range(self.n)             # a list of all the pool roots
                roots_group   = para.MPI.Group.Incl(comm.Get_group(), self.roots)
                self.rootcomm = comm.Create_group(roots_group) # communication among the roots of all pools
                # actual pool size (npp plus residue)
                self.size     = self.comm.Get_size()
            else:
                self.comm   = None
                self.roots  = [0]
                self.rank   = 0
                self.size   = 1
    
    def info(self):
        """ Collect pool information and print """
        para = self.para
        comm = para.comm
        if not self.pool_list:
            if comm and self.comm:
                mypool = self.comm.gather(para.rank, root = 0)
                if self.rootcomm != para.MPI.COMM_NULL:
                    self.pool_list = self.rootcomm.gather(mypool, root = 0)
                self.pool_list = comm.bcast(self.pool_list, root = 0)
            else:
                mypool = [0]
                self.pool_list = [mypool]
        para.print(' Setting up pools ...')
        para.print(' {0:<4}{1:>6}{2:<6}'.format('pool', 'size', '  proc'))
        for i in range(self.n):
            para.print(' {0:>4}{1:>6}{2:<6}'.format(i, len(self.pool_list[i]), '  ' + str(self.pool_list[i]).strip('[]')))
        para.print()
        
    def set_sklist(self, nspin = 1, nk = 1):
        """ 
        set up a list of spin and kpoint tuples to be processed on this proc

        example:
        nspin = 2, nk = 5
        self.n  =   3   # number of pools
        2 * 5 / 3 = 3   # each pool at least deal with 3 tuples

        pool    tuple                           offset
           0    (0, 0), (0, 3), (1, 1), (1, 4)  0, 3, 6, 9
           1    (0, 1), (0, 4), (1, 2)          1, 4, 7
           2    (0, 2), (1, 0), (1, 3)          2, 5, 8
        """
        self.nspin = nspin
        self.nk = nk
        self.sk_list = []
        self.sk_offset = []
        for s in range(nspin):
            for k in range(nk):
                offset = s * nk + k
                if offset % self.n == self.i:
                    self.sk_list.append((s, k))
                    self.sk_offset.append(offset)
        self.nsk = len(self.sk_list)

    def sk_info(self):
        """ collect and print out the spin-kpoint tuples on each pool """

        para = self.para
        if para.comm:
            self.sk_list_all = para.comm.gather(self.sk_list, root = 0)
        else:
            self.sk_list_all = self.sk_list
        
    def __init__(self, para = None):
        self.para = para
        self.pool_list = []
        self.set_pool(nproc_per_pool = 1)

        
class para_class(object):
    """ para class: wrap up user-defined mpi variables for parallelization """


    def __init__(self, MPI = None):
        if MPI is not None:
            self.MPI  = MPI
            self.comm = MPI.COMM_WORLD
            self.size = self.comm.Get_size()
            self.rank = self.comm.Get_rank()

        else:
            self.MPI  = None
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
        para = self.para
        if isanaconda():
            try:
                userin = open(sys.argv[1], 'r')
                para.print(' Reading user input from ' + sys.argv[1] + ' ...')
            except IOError:
                para.print(" Can't open user-defined input " + sys.argv[1] + " Halt. ", flush = True)
                para.stop()
        else: userin = sys.stdin
        lines = userin.read()
        var_input = input_arguments(lines)
        # para.print(var_input) # debug
        for var in set(vars(self)) & set(var_input): # This can be improved
            setattr(self, var, convert_val(var_input[var], type(getattr(self, var))))
            try:
                # convert var into correct data type as implied in __init__ and set attributes
                setattr(self, var, convert_val(var_input[var], type(getattr(self, var))))
            except:
                pass
        self.ipath = os.path.abspath(self.ipath)
        self.fpath = os.path.abspath(self.fpath)

        if not para.comm: self.nproc_per_pool = 1
        # para.print(vars(self)) # debug
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
        para = self.para
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
            if f == 'info': binary = '' # Need this because python3 will try to decode the binary file erroneously
            else: binary = 'b'
            try:
                fh = open(fname, 'r' + binary)
            except:
                para.print(" Can't open " + fname + '. Check if shirley_xas finishes properly. Halt. ', flush = True)
                para.stop()
            if f == 'info':
                lines = fh.read()
                var_input = input_arguments(lines, lower = True)
                #para.print(var_input) # debug
                for var in ['nbnd', 'nk', 'nelec', 'ncp', 'nspin']:
                    if var in var_input:
                        try:
                        # convert var into correct data type as implied in __init__ and set attributes
                            setattr(self, var, convert_val(var_input[var], type(getattr(self, var))))
                            #para.print(var) # debug
                            #para.print(convert_val(var_input[var], type(getattr(self, var)))) # debug
                        except:
                            para.print(' Convert ' + var + ' = ' + var_input[var] + ' failed.')
                    else:
                        para.print(' Variable "' + var + '" missed in ' + fname, flush = True)
                        para.stop()

                # initialize k-points
                # print(self.nspin, self.nk, user_input.gamma_only) # debug
                if user_input.gamma_only: self.nk = 1
                self.kpt = kpoints_class(nk = self.nk) # do we need this ?

                # Check k-grid consistency between the initial and final state
                if not is_initial:
                    if para.pool.nk != self.nk: # if the initial-state # of kpoints is not the same with the final-state one
                        para.print(' Inconsistent k-point number nk_i = ' + str(para.pool.nk) + ' v.s. nk_f = ' + str(self.nk) + ' Halt.')
                        para.stop()

                if is_initial:
                    # set up pools
                    if self.nk < para.size / user_input.nproc_per_pool:
                        para.print(' Too few k-points (' + str(self.nk) + ') for ' + str(int(para.size / user_input.nproc_per_pool)) + ' pools.')
                        para.print(' Increat nproc_per_pool to ' + str(int(para.size / self.nk)))
                        para.print()
                        user_input.nproc_per_pool = int(para.size / self.nk)

                    para.pool.set_pool(user_input.nproc_per_pool)
                    para.pool.info()

                # get spin and k-point index processed by this pool
                para.pool.set_sklist(nspin = self.nspin, nk = self.nk)

                # initialize obf
                self.obf = optimal_basis_set_class(nbnd = self.nbnd) # there're more: nbasis, ...

            if f == 'eigval':
                self.obf.input_eigval(fh)
            if f == 'eigvec':
                        pass
            fh.close()           
        
    def input(self, user_input, path, mol_name, is_initial):
        """ input from one shirley run """
        para = self.para
        if user_input.scf_type == 'shirley_xas':
            para.print(' Import wavefunctions and energies from shirley_xas calculation ...\n ')
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
