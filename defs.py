""" Definitions of classes """

from __future__ import print_function


import sys
import os
import struct

from constants import *
from utils import *
from io_mod import *

class user_input_class(object):
    """ input and record user-defined argument from stdin """

    def __init__(self):
        # user-defined variables and default values
        self.path_i     = '.'
        self.path_f     = '.'
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
                para.print(' Reading user input from ' + sys.argv[1] + ' ...\n')
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
        self.path_i = os.path.abspath(self.path_i)
        self.path_f = os.path.abspath(self.path_f)

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
    

    def __init__(self, nbasis = 0, nbnd = 0, nk = 1, nspin = 1, nelec = 0):
        sp = self.sp
        # variable list and default values
        self.nbasis     = nbasis
        self.nbnd       = nbnd
        self.nk         = nk
        self.nspin      = nspin
        self.nelec      = nelec
        self.eigval     = sp.array([])    # eigenvalues (band energies)
        self.eigvec     = sp.array([])    # eigenvectors (wavefunctions)

    # Should I input them at the pool root and then broadcast them ?
    def input_eigval(self, fh, sk_offset):
        sp = self.sp
        para = self.para
        pool = para.pool
        try:
            self.eigval = input_from_binary(fh, 'double', self.nbnd, sk_offset * self.nbnd)
        except struct.error:
            pass
        # Output a part of eigenvalues
        para.print('  ' + list2str_1d(self.eigval)) 
        self.eigval = sp.array(self.eigval)

    def input_eigvec(self, fh, sk_offset):
        sp = self.sp
        para = self.para
        pool = para.pool
        size = self.nbnd * self.nbasis
        try:
            """
            eigvec is given as a 1D array (transpose(eigvec) in column-major order)
            [ <B_1|1k>, <B_1|2k>, <B_2|1k>, <B_2|2k>]

            which corresponds to such a matrix:
            <B_1|1k>  <B_1|2k>
            <B_2|1k>  <B_2|2k>
            """
            # Note that the input is the transpose of < B_i | nk >
            self.eigvec = input_from_binary(fh, 'complex', size, sk_offset * size)
        except struct.error:
            pass
        # Output some major components of a part of eigenvalues near the fermi level
        para.print(eigvec2str(self.eigvec, self.nbasis, self.nbnd, int(self.nelec / (3 - self.nspin))))
        # Note that reshape works in row-major order
        self.eigvec = sp.array(self.eigvec).reshape(self.nbasis, self.nbnd)
        

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
        self.nbasis = 1                         # number of basis
        self.xmat   = sp.array([])              # single-particle matrix elements


    def input_shirley(self, user_input, is_initial, isk = 0):
        """ 
        input from shirley xas

        Arguments:
        is_initial: is this an initial-state scf
        isk: the index of the spin-k-point block to be input
             isk == 0 indicates this is the first time reading the data and 
             we need to extract the basic scf information from the *.info file.
        """
        para = self.para

        # stripe the i/f postfix
        if is_initial: postfix = '_i'
        else: postfix = '_f'
        
        path = getattr(user_input, 'path' + postfix)
        mol_name = getattr(user_input, 'mol_name' + postfix)

        # construct file names
        xas_prefix = mol_name + '.xas'
        xas_data_prefix = xas_prefix + '.' + str(user_input.xas_arg)
        ftype = []
        if isk < 0: ftype.append('info')
        if isk >= 0: 
            ftype += ['eigval', 'eigvec', 'proj']
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
                for var in ['nbnd', 'nk', 'nelec', 'ncp', 'nspin', 'nbasis']:
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

                # print out basis information
                info_str = ('  number of bands (nbnd)                    = {0}\n'\
                         +  '  number of spins (nspin)                   = {1}\n'\
                         +  '  number of k-points (nk)                   = {2}\n'\
                         +  '  number of electrons (nelec)               = {3} \n'\
                         +  '  number of optimal-basis function (nbasis) = {4} \n'\
                           ).format(self.nbnd, self.nspin, self.nk, self.nelec, self.nbasis)
                para.print(info_str)

                # initialize k-points
                # print(self.nspin, self.nk, user_input.gamma_only) # debug
                if user_input.gamma_only: self.nk = self.nspin
                self.kpt = kpoints_class(nk = self.nk) # do we need this ?

                # Check k-grid consistency between the initial and final state
                if not is_initial:
                    if para.pool.nk != self.nk: # if the initial-state # of kpoints is not the same with the final-state one
                        para.print(' Inconsistent k-point number nk_i = ' + str(para.pool.nk) + ' v.s. nk_f = ' + str(self.nk) + ' Halt.')
                        para.stop()
                    para.print(' Consistency check OK between the initial and final scf. \n')

                if is_initial and not para.pool.up:
                    # set up pools
                    if self.nk < para.size / user_input.nproc_per_pool:
                        para.print(' Too few k-points (' + str(self.nk) + ') for ' + str(int(para.size / user_input.nproc_per_pool)) + ' pools.')
                        para.print(' Increat nproc_per_pool to ' + str(int(para.size / self.nk)))
                        para.print()
                        user_input.nproc_per_pool = int(para.size / self.nk)

                    para.pool.set_pool(user_input.nproc_per_pool)
                    para.pool.info()

                # get spin and k-point index processed by this pool
                para.pool.set_sk_list(nspin = self.nspin, nk = self.nk)
                # para.pool.print(str(para.pool.sk_list) + ' ' + str(para.pool.sk_offset)) # debug

                # initialize obf
                # Typing this list of parameter again seems to be redundant. Use inheritance ?
                self.obf = optimal_basis_set_class(nbnd   = self.nbnd, 
                                                   nbasis = self.nbasis,
                                                   nk     = self.nk,
                                                   nspin  = self.nspin,
                                                   nelec  = self.nelec) # there're more: nbasis, ...

            if f == 'eigval':
                para.print('  Band energies: ')
                self.obf.input_eigval(fh, para.pool.sk_offset[isk])
                para.print()

            if f == 'eigvec':
                para.print('  Obf Wavefunctions : ')
                self.obf.input_eigvec(fh, para.pool.sk_offset[isk])
                para.print()

            fh.close()           
        
    def input(self, user_input, is_initial, isk):
        """ input from one shirley run """
        para = self.para
        if user_input.scf_type == 'shirley_xas':
            if isk < 0: para.print(' Wavefunctions and energies will be imported from shirley_xas calculation. \n ')
            self.input_shirley(user_input, is_initial, isk)
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
