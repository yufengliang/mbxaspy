""" mbxaspy main """

from __future__ import print_function

import sys
import os

from utils import *
from defs import *
from init import *
from spectra import *
from analysis import *

# input user-defined arguments from stdin
user_input.read()

# Check initial and final state and perform sanity checks
# Check the initial-state scf calculations
para.sep_line()
para.print(' Checking initial-state scf from: \n ' + user_input.path_i + '\n')
iscf.input(isk = -1)

# Check the final-state scf calculations
para.sep_line()
para.print(' Checking final-state scf from: \n ' + user_input.path_f + '\n')
fscf.input(is_initial = False, isk = -1, nelec = iscf.nelec)

# Important: Need to tell iscf the index of core
if user_input.scf_type == 'shirley_xas':
    iscf.proj.icore = fscf.proj.icore

# Input < B_i | \tilde{B}_j >
fscf.obf.input_overlap(user_input.path_f, iscf.nbnd, fscf.nbnd)

from xi import *
# Compute full atomic overlap sij
if user_input.scf_type == 'shirley_xas': 
    compute_full_sij(fscf.proj)

if not user_input.spec0_only:
    from determinants import *
    spec_xps_all = []
    spec_xas_all = []

nspin = iscf.nspin

# initialize the energy axis for spectra
global_ener_axis = spec_class(user_input).ener_axis
# if para.isroot(): sp.savetxt('ener_axis.dat', global_ener_axis) # debug

# global energy shift
global_offset = user_input.ESHIFT_FINAL + fscf.e_lowest

# initialize i and f spec0
spec0_i = [spec_class(ener_axis = global_ener_axis) for s in range(nspin)]
if user_input.final_1p: spec0_f = [spec_class(ener_axis = global_ener_axis) for s in range(nspin)]

ixyz_list = [-1, 0, 1, 2] # user_input.ixyz_list
# for many-body
ixyz_list_mb = ixyz_list[:]
if -1 in ixyz_list:
    for ixyz in [0, 1, 2]: if ixyz not in ixyz_list: ixyz_list_mb.append(ixyz)

## loop over spin and kpoints
for isk in range(pool.nsk):

    ispin, ik  = pool.sk_list[isk] # acquire current spin

    para.sep_line()
    para.print(' Processing (ispin, ik) = ({0},{1}) \n'.format(ispin, ik))

    # weight the sticks according to the k-grid
    weight = iscf.kpt.weight[ik]; 
    prefac = weight * Ryd
    # para.print('weight = {0}'.format(weight)) # debug

    # Import the initial-state scf calculation
    para.sep_line(second_sepl)
    para.print(' Importing initial-state scf\n')
    iscf.input(isk = isk)
    para.print('  xmat: the {0}th atom in ATOMIC_POSITIONS is excited.'.format(xatom(iscf.proj, iscf.xmat) + 1))

    # Import the final-state scf calculation
    para.sep_line(second_sepl)
    para.print(' Importing final-state scf\n', flush = True)
    fscf.input(is_initial = False, isk = isk)
    if user_input.final_1p: para.print('  xmat: the {0}th atom in ATOMIC_POSITIONS card is excited.'.format(xatom(fscf.proj, fscf.xmat) + 1))

    # Obtain the effective occupation number: respect the initial-state #electrons
    if nspin == 1: nocc = iscf.nocc
    else: nocc = iscf.nocc[ik][ispin]

    ## Compute non-interacting spectra *** should I put it in a def ?
    para.print('  Calculating one-body spectra ...\n', flush = True)

    # sticks = xmat_to_sticks(iscf, [-2], nocc, evec = [1.0, 0.0, 0.0]) # debug
    # print(sticks[0]) debug

    # initial-state
    sticks = xmat_to_sticks(iscf, ixyz_list, nocc, offset = -fscf.e_lowest, evec = user_input.EVEC)
    spec0_i[ispin].add_sticks(sticks, user_input, prefac, mode = 'additive')

    # final-state
    if user_input.final_1p:
        sticks = xmat_to_sticks(fscf, ixyz_list, nocc, offset = -fscf.e_lowest, evec = user_input.EVEC)
        spec0_f[ispin].add_sticks(sticks, user_input, prefac, mode = 'additive')
    
    para.print('  One-body spectra finished.', flush = True)

    ## Compute many-body spectra
    if not user_input.spec0_only:

        para.print('  Calculating many-body spectra ... ')

        ## Compute the transformation matrix xi
        para.print('  Calculating transformation matrix xi ... ')
        xi = compute_xi(iscf, fscf)

        if user_input.xi_analysis and para.isroot() and ik == 0:
            # plot_xi(xi) # debug
            if nspin > 1:
                msg = eig_analysis_xi(xi, '_spin_{0}'.format(ispin)) # debug
            else:
                msg = eig_analysis_xi(xi) # debug

        para.print('  Matrix xi finished. ', flush = True)

        ## XPS spectra (N X N)
        para.sep_line(second_sepl)
        para.print('  Calculating many-body XPS spectra ... ')

        Af_list, msg = quick_det(xi[:, 0 : int(nocc)], ener = fscf.obf.eigval,
                                 fix_v1 = False, maxfn = user_input.maxfn - 1,
                                 I_thr = user_input.I_thr,
                                 e_lo_thr = user_input.ELOW, e_hi_thr = user_input.EHIGH, 
                                 comm = pool.comm, 
                                 zeta_analysis = user_input.zeta_analysis and ik == 0)

        spec_xps_isk = spec_class(ener_axis = global_ener_axis)

        for order, Af in enumerate(Af_list):

            sticks = Af_to_sticks(Af)

            # important information for understanding shakeup effects and convergence 
            para.print("order {0:>2}: no. of sticks = {1:>7}, max stick = {2} ".
                        format( order, len(sticks), max([s[2] for s in sticks] + [0.0]) ))

            spec_xps_isk.add_sticks(sticks, user_input, mode = 'additive')

        spec_xps_all.append(spec_xps_isk)

        # output for debug
        postfix = '_ik{0}'.format(ik)
        if nspin == 2:
            postfix += '_ispin{0}'.format(ispin)
        postfix += '.dat'

        spec_xps_isk.savetxt(spec_xps_fname + postfix, offset = 0.0)
        
        para.print('  XPS spectra finished. ', flush = True)
        para.print()

        ## XAS spectra ( (N + 1) x (N + 1) )
        para.sep_line(second_sepl)
        para.print('  Calculating many-body XAS spectra ... ')

        spec_xas_isk = spec_class(ener_axis = global_ener_axis)

        for ixyz in ixyz_list_mb:

            # ixyz == -1 (non-polarized) is special. Can only be obtained by superposition of sticks
            if ixyz == -1:
                # deal with this at the end
                spec_xas_isk.add_sticks([], mode = 'append')
                continue

            para.print('  ixyz = {0}'.format(ixyz))
            # Compute xi_c
            ixmat = [ xmat_ixyz( iscf.xmat[ib, 0, :], ixyz, evec = user_input.EVEC ) for ib in range(iscf.nbnd) ]
            # xi_c = compute_xi_c(xi, iscf.xmat[:, 0, ixyz], nocc, iscf.nbnd_use)
            xi_c = compute_xi_c(xi, ixmat, nocc, iscf.nbnd_use)
            # para.print('xi_c.shape = {0}'.format(str(xi_c.shape))) # debug

            # Add the last column
            xi_c_ = sp.concatenate((xi[:, 0 : int(nocc)], xi_c), axis = 1)

            Af_list, msg = quick_det(xi_c_, ener = fscf.obf.eigval,
                                     fix_v1 = True, maxfn = user_input.maxfn,
                                     I_thr = user_input.I_thr,
                                     e_lo_thr = user_input.ELOW, e_hi_thr = user_input.EHIGH, 
                                     comm = pool.comm, 
                                     zeta_analysis = user_input.zeta_analysis and ik == 0)

            for order, Af in enumerate(Af_list):

                sticks = Af_to_sticks(Af)

                # important information for understanding shakeup effects and convergence 
                para.print("order {0:>2}: no. of sticks = {1:>7}, max stick = {2} ".
                            format( order + 1, len(sticks), max([s[2] for s in sticks] + [0.0]) ))

                spec_xas_isk.add_sticks(sticks, user_input, prefac, mode = 'append')

            para.print()
        # end of ixyz
        # go back and deal with average
        if -1 in ixyz_list_mb:
            spec_xas_isk.average([ixyz for ixyz in ixyz_list_mb if ixyz in [0, 1, 2]], ixyz_list_mb.index(-1))

        spec_xas_isk.savetxt(spec_xas_fname, offset = global_offset)
        spec_xas_all.append(spec_xas_isk)

        para.print('  Many-body XAS spectra finished. ')
        # end of ixyz
    # end if spec0_only
# end of isk

## Output one-body spectra

# intial-state one-body
for ispin in range(nspin): spec0_i[ispin].mp_sum(pool.rootcomm) 

if nspin == 1: spec0_i = spec0_i[0]
else:   spec0_i = spec0_i[0] | spec0_i[1] # mix spin up and down

if para.isroot(): spec0_i.savetxt(spec0_i_fname, offset = global_offset)

# final-state one-body
if user_input.final_1p:
    for ispin in range(nspin): spec0_f[ispin].mp_sum(pool.rootcomm) 

    if nspin == 1: spec0_f = spec0_f[0]
    else:   spec0_f = spec0_f[0] | spec0_f[1] # mix spin up and down

    if para.isroot(): spec0_f.savetxt(spec0_f_fname, offset = global_offset)

# spec0_sum = spec0_i[0] + spec0_f[0] # test operator overload
# spec0_sum.savetxt('spec0_sum.dat')

if user_input.spec0_only:
    para.done() # debug

## Calculate total many-body spectra 

# convolute spin-up and -down spectra if nspin == 2
if nspin == 2:

    # convolute xas spectra: do this before xps
    for isk, sk in enumerate(pool.sk_list):
        ispin, ik = sk
        if (1 - ispin, ik) not in pool.sk_list:
            pool.log('The twin tuple ({0}, {1}) for ({2},{3}) is not on this pool with {5}'
                    .format(ispin, ik, 1 - ispin, ik, str(pool.sk_list)))
            pool.log('Unsupported distribution of sk-tuples', flush = True)
            para.exit()
        ind = pool.sk_list.index((1 - ispin, ik))
        spec_xps_twin = spec_xps[ind]
        spec_xas_all[isk] *= spec_xps_twin

    # convolute xps spectra
    isk_done = []
    # Add spectra from each k-point
    for isk, sk in enumerate(pool.sk_list):
        if isk in isk_done: continue
        ispin, ik = sk
        ind = pool.sk_list.index((1 - ispin, ik)) # don't need to check existence again
        spec_xps_twin = spec_xps[ind]
        spec_xps_all[isk] *= spec_xps_twin
        # the two spin channels are the same for xps
        spec_xps_all[ind] = spec_xps_all[isk] # note that the twins are correlated now (use deepcopy to uncorrelate them)
        isk_done +=[isk, ind]

# local summation

spec_xps = spec_class(ener_axis = global_ener_axis)
spec_xas = [spec_class(ener_axis = global_ener_axis) for s in range(nspin)]

for isk, sk in enumerate(pool.sk_list):
    ispin, ik = sk
    weight = iscf.kpt.weight[ik]    
    if ispin == 0: spec_xps += spec_xps_all[isk] * weight # two spin channels are the same
    spec_xas[ispin] += spec_xas_all[isk]

# mpi reduce
spec_xps.mp_sum(pool.rootcomm)
for ispin in range(nspin):
    spec_xas[ispin].mp_sum(pool.rootcomm)

if nspin == 1: spec_xas = spec_xas[0]
else:   spec_xas = spec_xas[0] | spec_xas[1] # mix spin up and down

# This requires the world root is also one of the pool roots: can be made more robust
if para.isroot():
    postfix = '.dat'
    spec_xps.savetxt(spec_xps_fname + postfix)
    spec_xas.savetxt(spec_xas_fname + postfix, offset = global_offset)

para.done()
# Bye ! ~
