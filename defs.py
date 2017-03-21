""" Definitions of classes """

from __future__ import print_function


import sys
import os
import struct
import glob

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
        # *** You can do this via import scf only 
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
            para.error('Problem converting eigval file.')
        # Output a part of eigenvalues
        para.print('  ' + list2str_1d(self.eigval)) 
        self.eigval = sp.array(self.eigval)

    def input_eigvec(self, fh, sk_offset):
        sp = self.sp
        para = self.para
        pool = para.pool
        size = self.nbasis * self.nbnd
        try:
            """
            eigvec is given as a 1D array (transpose(eigvec) in column-major order)
            [ <B_1|1k>, <B_1|2k>, <B_2|1k>, <B_2|2k>]

            which corresponds to such a matrix:
            <B_1|1k>  <B_1|2k>
            <B_2|1k>  <B_2|2k>
            """
            # Note that the input is the transpose of < B_i | nk > ***
            self.eigvec = input_from_binary(fh, 'complex', size, sk_offset * size)
        except struct.error:
            para.error('Problem converting eigvec file.')
        # Output some major components of a part of eigenvalues near the fermi level
        para.print(eigvec2str(self.eigvec, self.nbasis, self.nbnd, int(self.nelec / (3 - self.nspin))))
        # Note that reshape works in row-major order
        self.eigvec = sp.array(self.eigvec).reshape(self.nbasis, self.nbnd)

    def input_overlap(self, nbnd1, nbnd2):
        """ 
        Input overlap < B_i | \tilde{B}_j > from file

        It is of the specified size nbnd1 x nbnd2
        """
        sp = self.sp
        para = self.para
        try:
            fh = open(overlap_fname, 'rb')
        except IOError:
            para.error('cannot open {0}'.format(overlap_fname))
        self.overlap = input_from_binary(fh, 'complex', nbnd1 * nbnd2, 0)
        # Note that in fortran this is in column-major order
        self.overlap = sp.array(self.overlap).reshape(nbnd2, nbnd1).transpose()
        fh.close()


class proj_class(object):
    """ store information related to the projectors used in the projector-argumentated wave (PAW) method """


    def __init__(self, scf = None):
        sp = self.sp
        # variable list and default values
        self.nspecies   = 0             # number of species of PAW atoms
        self.natom      = 0             # number of atoms
        self.scf        = scf           # so as to take the variables from the scf
        self.xs, self.x = -1, -1        # Assume no excited atoms
        self.l          = []            # angular momenta of all atom species
        self.qij        = []            # Q_int of all atom species
        self.nproj      = 0             # total number of projectors
        self.beta_nk    = sp.array([])  # < beta | nk > as in shirley_xas
        self.sij        = []            # PAW atomic overlap S_int between i and f
        if scf:
            self.import_from_iptblk(scf.tmp_iptblk)
            self.import_l_qij()


    def import_from_iptblk(self, iptblk):
        """ import atomic species and positions (names) from TMP_INPUT by shirley_xas """
        para = self.para
        self.atomic_species = atomic_species_to_list(iptblk['TMP_ATOMIC_SPECIES']) # element pseudopotential_file
        self.atomic_pos     = atomic_positions_to_list(iptblk['TMP_ATOMIC_POSITIONS']) # element x y z
        # para.print(self.atomic_species) # debug
        # para.print(self.atomic_pos) # debug
        # identify the excited atom if any
        self.nspecies   = len(self.atomic_species)
        self.natom      = len(self.atomic_pos)
        self.ind        = {}
        for i in range(self.nspecies):
            if self.atomic_species[i][0][-1] == 'X': self.xs = i
            self.ind[self.atomic_species[i][0]] = i
        for i in range(self.natom):
            if self.atomic_pos[i][0][-1] == 'X': self.x = i
        para.print('  number of atom species                    = {0}'.format(self.nspecies))
        para.print('  number of atoms                           = {0}'.format(self.natom))
        if self.x >= 0: para.print(' The excited atom is ' + str(self.x)) # debug


    def import_l_qij(self):
        """ import Q_int for all ground-state atoms in the supercell """
        para    = self.para
        scf     = self.scf
        para.print('  {0:6}{1:<30}{2:<20}'.format('Kind', 'UPF', 'Beta L'))
        # Import l and qij for each species of atom
        for i in range(self.nspecies): 
            pseudo_fname = scf.tmp_iptblk['TMP_PSEUDO_DIR'] + '/' + self.atomic_species[i][1]
            l, qij, errmsg = read_qij_from_upf(pseudo_fname)
            # para.print(l) # debug
            # para.print(qij) # debug
            if errmsg:
                para.error('read_qij_from_upf: {0}'.format(errmsg))
            # print out PAW atom information
            para.print('  {0:6}{1:<30}{2:<20}'.format(self.atomic_species[i][0], self.atomic_species[i][1], str(l).strip('[]')))
            self.l.append(l)
            self.qij.append(qij)
        # Calculate total # of projectors in the system
        for i in range(self.natom):
            self.nproj += sum([2 * l + 1 for l in self.l[self.ind[self.atomic_pos[i][0]]]])
        para.print('  number of projectors = {0}'.format(self.nproj))


    def input_sij(self):
        """ 
        input the atomic overlap term sij

        There should be one sij file for one UPF file of each excited atom.
        For example:
        O.pbe-van-yufengl-1s1.sij => O.pbe-van-yufengl-1s1.UPF
        So we would need to look for sij for a given UPF
        """
        scf = self.scf
        para = self.para
        if self.x >= 0: # if there exists an excited atom
            fname = self.atomic_species[self.xs][1]
            fname = scf.tmp_iptblk['TMP_PSEUDO_DIR'] + '/' + os.path.splitext(fname)[0] + '.sij'
            try:
                fh = open(fname, 'r')
            except IOError:
                para.error('cannot open sij file: {0}'.format(fname))
            for line in fh:
                self.sij.append([float(_) for _ in line.split()])
            fh.close()
            para.print('  Sij imported from {0}'.format(fname))
            # para.print('sij = ' + str(self.sij)) # debug


    def input_proj(self, fh, sk_offset):
        """ input projectors from shirley_xas calculation """
        para    = self.para
        scf     = self.scf
        sp      = self.sp
        size    = self.nproj * scf.nbnd
        try:
            self.beta_nk = input_from_binary(fh, 'complex', size, sk_offset * size)
        except struct.error:
            para.error('Problem converting proj file.')
        self.beta_nk = sp.array(self.beta_nk).reshape(self.nproj, scf.nbnd)


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

        user_input: user_input_class that stores user input arguments

        is_initial: is this an initial-state scf ?

        isk: the index of the spin-k-point block to be input
             isk == 0 indicates this is the first time reading the data and 
             we need to extract the basic scf information from the *.info file.
        """
        para = self.para

        # stripe the i/f postfix
        if is_initial: postfix = '_i'
        else: postfix = '_f'
        
        # construct the path to the scf calculation and the file prefix
        path = getattr(user_input, 'path' + postfix)
        mol_name = getattr(user_input, 'mol_name' + postfix)

        # import Input_Block.in from shirley_xas calculation if the first time to read
        if isk < 0:
            # fname = os.path.abspath(path + '/' + iptblk_fname) # Input_Block.in
            fname = sorted(glob.glob(path + '/' + tmp_iptblk_fname + '*'))[-1]
            with open(fname, 'r') as fh:
                lines = fh.read()
            self.tmp_iptblk = input_arguments(lines) # store variables in iptblk
            # para.print(self.tmp_iptblk['TMP_ATOMIC_SPECIES']) # debug
            # para.print(self.tmp_iptblk['TMP_ATOMIC_POSITIONS']) # debug

        # construct file names
        xas_prefix = mol_name + '.xas'
        xas_data_prefix = xas_prefix + '.' + str(user_input.xas_arg)

        # Figure out which types of files to input
        ftype = []
        if isk < 0: ftype += ['info', 'proj'] # if this is the first time to read, then read the information first
        else: 
            ftype += ['eigval', 'eigvec', 'proj']
            if is_initial: ftype += ['xmat'] # need pos matrix element for the initial state
        
        # Open and read the relevant files
        for f in ftype:
            fname = os.path.abspath(path + '/' + xas_data_prefix + '.' + f)
            if f == 'info': binary = '' # Need this because python3 will try to decode the binary file erroneously
            else: binary = 'b'
            try:
                fh = open(fname, 'r' + binary)
            except:
                para.print(" Can't open " + fname + '. Check if shirley_xas finishes properly. Halt. ', flush = True)
                para.stop()

            # information file
            if f == 'info':
                lines = fh.read()
                var_input = input_arguments(lines, lower = True)
                #para.print(var_input) # debug
                # take specific variables of interest from the info file
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
                         +  '  number of electrons (nelec)               = {3}\n'\
                         +  '  number of optimal-basis function (nbasis) = {4}\n'\
                           ).format(self.nbnd, self.nspin, self.nk, self.nelec, self.nbasis)
                para.print(info_str)

                # initialize k-points
                # print(self.nspin, self.nk, user_input.gamma_only) # debug
                if user_input.gamma_only: self.nk = self.nspin
                self.kpt = kpoints_class(nk = self.nk) # do we need this ?

                # Check k-grid consistency between the initial and final state *** We should move this check outside 
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

            # eigenvalue file
            if f == 'eigval':
                para.print('  Band energies: ')
                self.obf.input_eigval(fh, para.pool.sk_offset[isk])
                para.print()

            # eigenvector file
            if f == 'eigvec':
                para.print('  Obf Wavefunctions : ')
                self.obf.input_eigvec(fh, para.pool.sk_offset[isk])
                para.print()

            # projector file
            if f == 'proj':
                if isk < 0:
                    # initialize proj
                    para.print('  Reading information on PAW atoms ... ')
                    self.proj = proj_class(self)
                    if not is_initial:
                        self.proj.input_sij()
                else:
                    para.print('  Reading projectors ... ')
                    self.proj.input_proj(fh, para.pool.sk_offset[isk])
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
