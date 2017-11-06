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
        self.path_i         = '.'
        self.path_f         = '.'
        self.mol_name_i     = 'mol_name'
        self.mol_name_f     = 'xatom'
        self.xas_arg        = 5
        self.nelec          = -1
        self.scf_type       = 'shirley_xas'
        self.nproc_per_pool = 1

        # spectral information
        self.ELOW           = -2.0      # eV
        self.EHIGH          = 8.0       # eV
        self.ESHIFT_FINAL   = 0.0       # eV
        self.NENER          = 100       # energy grid size
        self.SIGMA          = 0.2       # eV
        self.EVEC           = None      # electric field vector E
        self.smearing       = 'gauss'   # gauss or lor (lorentzian)

        self.nbnd_i         = nbnd_max
        self.nbnd_f         = nbnd_max
        self.nk_use         = 0         # no. of k-points used

        self.maxfn          = 2         # final-state shakeup order
        self.I_thr          = 1e-3      # intensity cutoff
        
        # RIXS parameters
        self.RIXS_I_thr     = 1e-3          # rixs intensity cutoff
        self.eloss_range    = 10            # the range of energy loss (below ELOW) in eV
        self.loss_mode      = False         # True: output RIXS(x = omega_out = omega_in - eloss, y = omega_in)
                                            # False: output RIXS(x = omega_in, y = energy loss)
        self.inout_pol      = 'xx yy zz'    # polarization of incoming and outgoing photon. can also be xy, yz, zx ...
                                            # xy means the in-photon polarized along x, out-photon along y
        self.NENER_out      = 100           # omega_out or eloss grid size
        self.SIGMA_out      = 0.2           # lifetime broadening of the RIXS final state
        self.draw_rixs      = False         # output RIXS spectra with matplotlib
        self.rixs_nmajor    = 10            # number of major transitions to output
       
        # control flags
        self.gamma_only         = False         # Using Gamma-point only
        self.final_1p           = False         # Calculate one-body final-state spectra
        self.xi_analysis        = False         # perform full analysis on the xi matrix
        self.zeta_analysis      = False         # perform full analysis on the zeta matrix
        self.do_paw_correction  = True          # perform PAW corrections
        self.use_pos            = True          # Use the .pos files instead of the .xmat files (recommended)
        self.spec0_only         = False         # Calculate one-body / non-interacing spectra only
        self.xps_only           = False         # Calculate one-body and XPS spectra only
        self.want_bse           = False         # Want to calculate BSE spectra
        self.want_spec_o        = False         # Convolute the spec0_i with spec_xps
        self.spec_analysis      = False         # perform analysis on spectra

    def read(self):
        """ input from stdin or userin"""
        para = self.para
        if isanaconda():
            try:
                userin = open(sys.argv[1], 'r')
                para.print(' Reading user input from {0}\n'.format(sys.argv[1]))
            except IOError:
                para.error(" Can't open user-defined input {0}" .format(sys.argv[1]))
        else: userin = sys.stdin
        lines = userin.read()
        var_input = input_arguments(lines)
        # para.print(var_input) # debug
        for var in set(vars(self)) & set(var_input): # This can be improved
            try:
                # convert var into correct data type as implied in __init__ and set attributes
                setattr(self, var, convert_val(var_input[var], type(getattr(self, var))))
            except:
                pass
        self.path_i = os.path.abspath(self.path_i)
        self.path_f = os.path.abspath(self.path_f)

        # convert str EVEC into a list
        if self.EVEC is not None:
            try:
                evec = [float(e) for e in self.EVEC.split()]
                self.EVEC = evec  
            except:
                para.error(" E-field vector not correct. ")

        # test polarization format of inout_pol
        for pol_str in self.inout_pol.split():
            if len(pol_str) != 2:
                para.error("The length of '{}' is not 2. Only accept two polarization directions.".format(pol_str))
            for s in pol_str:
                if s not in pol_index:
                    para.error("Don't recognize polarization symbol {} in {}".format(s, pol_str))

        if not para.comm: self.nproc_per_pool = 1
        # para.print(vars(self)) # debug
        userin.close()


class kpoints_class(object):
    """ store information related to kpoints """

    def __init__(self, nk = 1):
        # variable list and default values
        self.nk         = nk
        self.weight     = [2.0 / nk] * nk
        # in future there may be a list of k-vectors and nontrivial weights (not needed now)


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
        self.eigvec     = sp.matrix([])    # eigenvectors (wavefunctions)

    def input_eigval(self, fh, sk_offset, output_msg = True, mid = -1):
        sp = self.sp
        para = self.para
        try:
            self.eigval = input_from_binary(fh, 'double', self.nbnd, sk_offset * self.nbnd)
        except struct.error:
            para.error('Problem converting eigval file.')
        self.eigval = [ e * Ryd for e in self.eigval ] # Don't forget Rydberg
        # Output a part of eigenvalues
        if output_msg: para.print('  ' + list2str_1d(self.eigval, mid)) 
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
        except: # struct.error:
            para.error('Problem converting eigvec file.')
        # Output some major components of a part of eigenvalues near the fermi level
        para.print(eigvec2str(self.eigvec, self.nbasis, self.nbnd, int(self.nelec / 2)), flush = True)
        # Note that reshape works in row-major order 
        # eigvec = < B_i | nk > 
        self.eigvec = sp.matrix(self.eigvec).reshape(self.nbasis, self.nbnd)

    def input_overlap(self, path, nbnd1, nbnd2):
        """ 
        Input overlap < B_i | \tilde{B}_j > from file

        It is of the specified size nbnd1 x nbnd2
        """
        fname = path + '/' + overlap_fname
        try:
            # fh = open(overlap_fname, 'rb') # Fortran unformatted output may vary depending on compilers
            fh = open(fname, 'r')
        except IOError:
            self.para.error('cannot open {0}'.format(overlap_fname))
        # self.overlap = input_from_binary(fh, 'complex', nbnd1 * nbnd2, 0) # binary
        self.overlap = []
        for line in fh:
            self.overlap.append(float(line.split()[0]) + 1j * float(line.split()[1]))
        # Note that in fortran this is in column-major order
        try:
            self.overlap = self.sp.matrix(self.overlap).reshape(nbnd2, nbnd1).transpose()
        except ValueError as vle:
            self.para.error(str(vle) + '\n Insufficient data in {}({})for nbnd1 = {} and nbnd2 = {}'.format(fname, len(self.overlap), nbnd1, nbnd2))
        # self.para.print(self.overlap[0:2, 0:5]) # debug
        fh.close()
        self.para.print('  Overlap matrix < B_i | ~B_j > imported from {0}'.format(fname))


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
        self.nprojs     = []            # number of projectors for each atom
        self.beta_nk    = sp.matrix([]) # < beta | nk > as in shirley_xas
        self.sij        = []            # PAW atomic overlap S_int between i and f
        self.icore      = 0             # index of this excited atom among all excited atoms
        if scf:
            self.import_from_iptblk(scf.tmp_iptblk)
            self.import_l_qij()
            self.find_icore()


    def import_from_iptblk(self, tmp_iptblk):
        """ import atomic species and positions (names) from TMP_INPUT by shirley_xas """
        para = self.para
        self.atomic_species = atomic_species_to_list(tmp_iptblk['TMP_ATOMIC_SPECIES']) # element pseudopotential_file
        self.atomic_pos     = atomic_positions_to_list(tmp_iptblk['TMP_ATOMIC_POSITIONS']) # element x y z
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
        if self.x >= 0: 
            para.print('  The {0}th atom in ATOMIC_POSITIONS ({1}) is excited.'.format(self.x + 1, self.atomic_pos[self.x][0][:-1])) # debug


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
            # Calculate # projectors for each kind of atom
            self.nprojs.append(sum([2 * _ + 1 for _ in l]))
        # Calculate total # of projectors in the system
        self.iproj2atom = []
        for i in range(self.natom):
            nproj = self.nprojs[self.ind[self.atomic_pos[i][0]]]
            self.nproj += nproj
            # Find which atom each projector belongs to
            self.iproj2atom += [i] * nproj 
        para.print('  number of projectors = {0}'.format(self.nproj))


    def find_icore(self):
        """ Calculate the index of this excited atom among all the excited atoms """
        self.ind_excitation = [0] * self.natom
        for key in self.scf.iptblk:
            if 'IND_EXCITATION' in key:
                self.ind_excitation[get_index(key) - 1] = int(self.scf.iptblk[key])
        # self.para.print(self.ind_excitation) # debug
        self.ncore = sum(self.ind_excitation)
        if self.x > 0: 
            self.icore = sum(self.ind_excitation[0 : self.x])
            self.para.print('  This is the {0}th atom among {1} core-excited atoms.\n'.format(self.icore + 1, self.ncore))

    def get_s(self, ith):
        """ Get the species of ith atom in the ATOMIC_POSITION list"""
        return  self.ind[ self.atomic_pos[ith][0] ] if ith < self.natom else None

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
        scf     = self.scf
        size    = self.nproj * scf.nbnd
        try:
            self.beta_nk = input_from_binary(fh, 'complex', size, sk_offset * size)
        except struct.error:
            self.para.error('Problem converting proj file.')
        self.beta_nk = self.sp.matrix(self.beta_nk).reshape(scf.nbnd, self.nproj).transpose()


class scf_class(object):
    """ pipeline data from self-consistent-field calculations """

    def __init__(self):
        sp = self.sp
        # variable list and default values
        self.nbnd       = 0                         # number of bands    
        self.nbnd_use   = 0                         # number of bands used actually  
        self.nk         = 0                         # number of k-points
        self.nk_use     = 0                         # number of k-points to be used
        self.nelec      = 0                         # number of electrons
        self.ncp        = 1                         # number of core levels
        self.nspin      = 1                         # number of spins
        self.nbasis     = 1                         # number of basis
        self.xmat       = sp.array([])              # single-particle matrix elements
        self.e_lowest   = None                      # the energy of the lowest-lying empty/partiall occupied state in this scf


    def input_xmat(self, fh, offset, sk_offset, is_initial = True):
        """ input matrix elements from fh """
        size = self.nbnd * self.ncp * nxyz
        try:
            self.xmat = input_from_binary(fh, 'complex', size, offset + sk_offset * size)
        except struct.error:
            if self.userin.final_1p and not is_initial:
                self.para.print('Problem converting xmat file. Skip one-body final-state spectrum. ')
                self.userin.final_1p = False
                return
            else:
                self.para.error('Problem converting xmat file.')
        # Note that the indexing in fortran is reversed
        # In shirley_xas, posn is nbnd x ncp x nxyz 3d array; xmat is the same here
        # xmat is calculated as < nk | O | phi_h > (phi_h being the core levels)
        self.xmat = self.sp.array(self.xmat).reshape(nxyz, self.ncp, self.nbnd).T

    def calc_xmat(self):
        """ 
        Calculate xmat from < beta | nk > and the *.pos file

        < nk | r_i | h_c > = sum_{l} < nk | beta_l > < beta_l | r_i | h_c > 
        
        l loops over the beta functions of the excited atom.
        < nk | beta_l > is from proj.bet_nk.
        < beta_l | r_i | h > is extracted from the pos file. 
        """
        para    = self.para
        proj    = self.proj
        sp      = self.sp
        # find the pos file
        xs = proj.get_s(proj.x) # Note that you can't use proj.xs for GS because there isn't any excited-atom species.
        pseudo_fname = proj.atomic_species[xs][1]
        pos_fname = self.tmp_iptblk['TMP_PSEUDO_DIR'] + '/' + os.path.splitext(pseudo_fname)[0] + '.pos'
        try:
            fh = open(pos_fname, 'r')
        except IOError:
            para.error('cannot open the pos file {}'.format(pos_fname))
            
        para.print('  Reading the pos file {}'.format(pos_fname))
        lwfc1, lwfc2, elem = import_from_pos(fh) # lwfc1: an l array, lwfc2: an l number, elem: a list
        fh.close()

        if elem is None:
            para.error('Problem reading the pos file.')
            
        # consistency check between the projectors in the pos file and the proj type
        if lwfc1 != proj.l[xs]:
            para.error('lwfc1 in the pos file {} not consistent with the pseudopotential {}. '.format(lwfc1, proj.l[xs]))
        
        # find the projector offset for the excited atom
        # proj.x for iscf has been assigned by proj.x from fscf in main.py in the pre-input step (isk < 0)
        proj_offset = 0
        for I in range(proj.x):
            proj_offset +=  proj.nprojs[proj.get_s(I)]

        # calculate < nk | r_i | h_c >
        self.xmat = sp.zeros((self.nbnd, 2 * lwfc2 + 1, nxyz), dtype = sp.complex128)
        for pos in elem:
            lm_valence, m_core, ixyz, pos_val = pos[0] - 1, pos[1] - 1, pos[2] - 1, pos[3]
            # Note that beta_nk has been converted into a matrix type. You can't broadcast it into an array directly
            for b in range(self.nbnd):
                # In future I should rewrite this using matrix multiplication
                self.xmat[b, m_core, ixyz] += proj.beta_nk[proj_offset + lm_valence, b].conjugate() * pos_val.conjugate()

    def input_shirley(self, is_initial = True, isk = 0, nelec = -1):
        """ 
        input from shirley xas

        Arguments:

        user_input: user_input_class that stores user input arguments

        is_initial: is this an initial-state scf ?

        isk: the index of the spin-k-point block to be input
             isk == 0 indicates this is the first time reading the data and 
             we need to extract the basic scf information from the *.info file.
        
        nelec: overwrite no. of electrons
        """
        para    = self.para
        userin  = self.userin

        # stripe the i/f postfix
        if is_initial: postfix = '_i'
        else: postfix = '_f'
        
        # construct the path to the scf calculation and the file prefix
        path, mol_name, nbnd_use = tuple([getattr(userin, _ + postfix) for _ in ['path', 'mol_name', 'nbnd']])

        # import Input_Block.in and TMP_INPUT* from shirley_xas calculation if the first time to read
        if isk < 0:

            fname = path + '/' + iptblk_fname # Input_Block.in
            try:
                with open(fname, 'r') as fh:
                    lines = fh.read()
            except IOError:
                para.error('Problem reading {0}'.format(fname))
            self.iptblk = input_arguments(lines) # store variables in iptblk

            # Use the latest TMP_INPUT file
            tmp_file_list = sorted(glob.glob(path + '/' + tmp_iptblk_fname + '*'))
            if not tmp_file_list:
                para.error('cannot find any {0} file in {1}'.format(tmp_iptblk_fname, path))
            fname = tmp_file_list[-1]
            try:
                with open(fname, 'r') as fh:
                    lines = fh.read()
            except IOError:
                para.error('Problem reading {0}'.format(fname))
            self.tmp_iptblk = input_arguments(lines) # store variables in iptblk

            # para.print(self.iptblk['IND_EXCITATION[0]']) # debug
            # para.print(self.tmp_iptblk['TMP_ATOMIC_SPECIES']) # debug
            # para.print(self.tmp_iptblk['TMP_ATOMIC_POSITIONS']) # debug

        # construct file names
        xas_prefix = mol_name + '.xas'
        xas_data_prefix = xas_prefix + '.' + str(userin.xas_arg)

        # Figure out which types of files to input
        if isk < 0: ftype = ['info', 'eigval', 'proj'] # if this is the first time to read, then read the information first
        else: 
            ftype = ['eigval', 'eigvec', 'proj']
            if is_initial or userin.final_1p: ftype += ['xmat'] # need pos matrix element for the initial state
        
        # Open and read the relevant files
        for f in ftype:
            fname = os.path.abspath(path + '/' + xas_data_prefix + '.' + f)
            if f == 'info': binary = '' # Need this because python3 will try to decode the binary file erroneously
            else: binary = 'b'
            try:
                fh = open(fname, 'r' + binary)
            except:
                if f != 'xmat' or not userin.use_pos:
                    para.error(" Can't open " + fname + '. Check if shirley_xas finishes properly. Halt. ')

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
                        para.error(' Variable "' + var + '" missed in ' + fname)

                # adjust the no. of bands actually used
                self.nbnd_use = nbnd_use if 0 < nbnd_use < self.nbnd else self.nbnd

                # overwrite no. of electrons
                if nelec >= 0: self.nelec = nelec

                # print out basis information
                info_str = ('  number of bands (nbnd)                    = {0}\n'\
                         +  '  number of bands used (nbnd_use)           = {1}\n'\
                         +  '  number of spins (nspin)                   = {2}\n'\
                         +  '  number of k-points (nk)                   = {3}\n'\
                         +  '  number of electrons (nelec)               = {4}\n'\
                         +  '  number of optimal-basis function (nbasis) = {5}\n'\
                           ).format(self.nbnd, self.nbnd_use, self.nspin, self.nk, self.nelec, self.nbasis)
                para.print(info_str, flush = True)

                # check no. of electrons
                if self.nelec > self.nbnd * 2:
                    para.error('Too many electrons ({0}) for {1} bands !'.format(self.nelec, self.nbnd))

                # initialize k-points
                # print(self.nspin, self.nk, userin.gamma_only) # debug
                # self.nk_use = 1 if userin.gamma_only else self.nk
                self.nk_use = userin.nk_use if userin.nk_use > 0 else self.nk
                if userin.gamma_only:
                    self.nk_use = 1

                if self.nk_use > self.nk: para.error('Number of kpoints to be used ({0}) larger than kpoints provided ({1})'.format(self.nk_use, self.nk))
                
                self.kpt = kpoints_class(nk = self.nk_use)

                # Check k-grid consistency between the initial and final state *** We should move this check outside 
                if not is_initial:
                    if para.pool.nk != self.nk: # if the initial-state # of kpoints is not the same with the final-state one
                        para.error(' Inconsistent k-point number nk_i = ' + str(para.pool.nk) + ' v.s. nk_f = ' + str(self.nk) + ' Halt.')
                    para.print(' Consistency check OK between the initial and final scf. \n')

                if is_initial and not para.pool.up:
                    # set up pools
                    nsk = self.nk_use * self.nspin
                    # if (0, k) and (1, k) can be treated on different pools
                    # if nsk < para.size / userin.nproc_per_pool:
                    #    para.print(' Too few (spin, k) tuples ({0}) to calculate for {1} pools.'.format(nsk, int(para.size / userin.nproc_per_pool)))
                    #    userin.nproc_per_pool = int(para.size / nsk)
                    #    para.print(' Increase nproc_per_pool to {0}\n'.format(int(userin.nproc_per_pool))

                    # if (0, k) and (1, k) can only be treated on the same pool
                    if self.nk_use < para.size / userin.nproc_per_pool:
                        para.print(' Too few k-points ({0}) to calculate for {1} pools.'.format(self.nk_use, int(para.size / userin.nproc_per_pool)))
                        userin.nproc_per_pool = int(para.size / self.nk_use) # for the contiguous mode
                        para.print(' Increase nproc_per_pool to {0}\n'.format(int(userin.nproc_per_pool)))

                    para.pool.set_pool(userin.nproc_per_pool)
                    para.pool.info()

                    # get spin and k-point index processed by this pool
                    para.pool.set_sk_list(nspin = self.nspin, nk = self.nk, nk_use = self.nk_use)
                    # para.pool.print(str(para.pool.sk_list) + ' ' + str(para.pool.sk_offset)) # debug
                    para.pool.sk_info()

                # initialize obf
                # Typing this list of parameter again seems to be redundant. Use inheritance ?
                self.obf = optimal_basis_set_class(nbnd   = self.nbnd, 
                                                   nbasis = self.nbasis,
                                                   nk     = self.nk,
                                                   nspin  = self.nspin,
                                                   nelec  = self.nelec) # there're more: nbasis, ...

            # eigenvalue file
            if f == 'eigval':
                if isk < 0:
                    # determine the occupation number for each k-point
                    if self.nspin == 1:
                        nocc = self.nocc = self.nelec / 2.0
                        for k in range(self.nk_use):
                            offset = k
                            self.obf.input_eigval(fh, offset, output_msg = False)
                            if not self.e_lowest or self.obf.eigval[int(nocc)] < self.e_lowest:
                                self.e_lowest = self.obf.eigval[int(nocc)]
                    else:
                        self.nocc = []
                        for k in range(self.nk_use):
                            eigval = [None] * 2
                            for s in range(self.nspin):
                                offset = s * self.nk + k
                                self.obf.input_eigval(fh, offset, output_msg = False)
                                eigval[s] = list(self.obf.eigval)
                            # occupation numbers for spin up and down channels for this k
                            nocc = find_nocc(eigval, self.nelec)
                            self.nocc.append(nocc)
                            emin = min(eigval[0][int(nocc[0])], eigval[1][int(nocc[1])])
                            if not self.e_lowest or emin < self.e_lowest:
                                self.e_lowest = emin
                    para.print('  Energy of LUMO: {0} eV '.format(self.e_lowest), flush = True)

                if isk >= 0:
                    para.print('  Band energies (eV): ')
                    if self.nspin == 1: nocc = self.nocc
                    else:
                        ispin, ik  = para.pool.sk_list[isk] 
                        nocc = self.nocc[ik][ispin]
                    self.obf.input_eigval(fh, para.pool.sk_offset[isk], mid = int(nocc))
                    self.eigval = self.obf.eigval
                    para.print('  occupation number: {0}'.format(nocc))
                    para.print('  local efermi = {0:.4f}'.format(self.eigval[int(nocc) - 1]))
                    para.print(flush = True)

            # eigenvector file
            if f == 'eigvec':
                para.print('  Obf Wavefunctions : ')
                self.obf.input_eigvec(fh, para.pool.sk_offset[isk])
                self.obf.eigvec = self.obf.eigvec[:, : self.nbnd_use] # resize to save memory
                self.eigvec = self.obf.eigvec
                para.print(flush = True)

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
                para.print(flush = True)

            # xmat file: matrix elements
            if f == 'xmat':
                para.print('  Reading single-body matrix elements ... ')

                if userin.use_pos:
                    # calculate xmat from < beta | nk > and the *.pos file in the pseudo library
                    self.calc_xmat()
                else:
                    # Notes:
                    # There is a problem with the xmat file produced by shirley_xas.overlap
                    # It's also easy to make mistakes to use the xmat files !!!

                    # This is important: for the ground-state system, xmat for all excited atoms
                    # are stored in the same xmat file. You need to set up an offset to locate the 
                    # right block for this excited atom being processed
                    if is_initial: 
                        size = self.nk * self.nbnd * self.ncp * nxyz
                        offset = size * self.proj.icore
                    else: offset = 0
                    self.input_xmat(fh, offset, para.pool.sk_offset[isk], is_initial)

                para.print(flush = True)
                
            fh.close()
        
        
    def input(self, is_initial = True, isk = -1, nelec = -1):
        """ input from one shirley run """
        if self.userin.scf_type == 'shirley_xas':
            if isk < 0: self.para.print(' Wavefunctions and energies will be imported from shirley_xas calculation. \n ')
            self.input_shirley(is_initial, isk, nelec)
        else:
            self.para.error(' Unsupported scf input: ' + scf_type)
            


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
