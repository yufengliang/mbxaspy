""" mbxaspy main """

from __future__ import print_function

import sys
import os

from utils import *
from defs import *
from init import *
from spectra import *

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
fscf.input(is_initial = False, isk = -1)

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
    spec_xps = []
    spec_xas = []

nspin = iscf.nspin
# loop over spin and kpoints
for isk in range(para.pool.nsk):

    ispin, ik  = para.pool.sk_list[isk] # acquire current spin
    para.print(' Processing (ispin, ik) = ({0},{1}) \n'.format(ispin, ik))
    # weight the sticks according to the k-grid
    weight =  iscf.kpt.weight[ik]; 
    # para.print('weight = {0}'.format(weight)) # debug

    # Import the initial-state scf calculation
    para.sep_line()
    para.print(' Importing initial-state scf\n')
    iscf.input(isk = isk)

    # Import the final-state scf calculation
    para.sep_line(second_sepl)
    para.print(' Importing final-state scf\n', flush = True)
    fscf.input(is_initial = False, isk = isk)

    # Obtain the effective occupation number: respect the initial-state #electrons
    nocc = eff_nocc(iscf.nelec, nspin, ispin)

    # Compute non-interacting spectra *** should I put it in a def ?
    for ixyz in range(3):
        # *** currently only x, y, z
        ener_axis, spec0 = spectrum0(iscf, ixyz, nocc, user_input.smearing)
        ener_axis += user_input.ESHIFT_FINAL
        # print(ener_axis.shape, spec0.shape) # debug
        if isk == 0 and ixyz == 0:
            # initialize spectra and energy axis
            spec0_i = sp.zeros([len(ener_axis), nspin * 4 + 1]) # ener (total x y z) * nspin
            spec0_i[:, 0] = ener_axis
            if user_input.final_1p: spec0_f = spec0_i.copy()
        col = (1 + ixyz) * iscf.nspin + (ispin + 1)
        spec0 *= weight
        spec0_i[:, col] += spec0
        spec0_i[:, ispin + 1] += spec0 # angular average
        if user_input.final_1p:
            ener_axis, spec0 \
            = spectrum0(fscf, ixyz, nocc, user_input.smearing)
            spec0 *= weight
            spec0_f[:, col] += spec0
            spec0_f[:, ispin + 1] += spec0 # angular average

    if not user_input.spec0_only:

        para.print('  Calculating many-body spectra ... ')

        # Compute the transformation matrix xi
        xi = compute_xi(iscf, fscf)

        #para.print(xi) # debug
        if user_input.xi_analysis and para.isroot() and ik == 0:
            # plot_xi(xi) # debug
            if nspin > 1:
                msg = eig_analysis_xi(xi, '_spin_{0}'.format(ispin)) # debug
            else:
                msg = eig_analysis_xi(xi) # debug

        ## XPS spectra (N X N)
        para.sep_line(second_sepl)
        para.print('  Calculating many-body XPS spectra ... ')

        Af_list, msg = quick_det(xi[:, 0 : int(nocc)], ener = fscf.obf.eigval,
                                 fix_v1 = False, maxfn = user_input.maxfn - 1,
                                 I_thr = user_input.I_thr,
                                 e_lo_thr = user_input.ELOW, e_hi_thr = user_input.EHIGH, 
                                 comm = para.pool.comm, 
                                 zeta_analysis = user_input.zeta_analysis and ik == 0)

        first = True
        for order, Af in enumerate(Af_list):

            stick = Af_to_stick(Af)
            ener_axis, spec = stick_to_spectrum(stick, user_input)
            ener_axis += user_input.ESHIFT_FINAL + fscf.obf.eigval[int(nocc)]

            # important information for understanding shakeup effects and convergence 
            para.print("order {0:>2}: no. of sticks = {1:>7}, max stick = {2} ".
                        format( order, len(stick), max([s[1] for s in stick] + [0.0]) ))

            if first:
                spec_xps_ = sp.zeros([len(ener_axis), 2])
                spec_xps_[:, 0] = ener_axis
                first = False

            spec_xps_[:, 1] += spec

        para.print()

        ## XAS spectra ( (N + 1) x (N + 1) )
        para.sep_line(second_sepl)
        para.print('  Calculating many-body XAS spectra ... ')

        first = True
        for ixyz in range(3):

            para.print('  ixyz = {0}'.format(ixyz))
            # Compute xi_c
            xi_c = compute_xi_c(xi, iscf.xmat[:, 0, ixyz], nocc, user_input.nbnd_i)
            # xi_c = compute_xi_c(xi, iscf.xmat[:, 0, ixyz], nocc)
            # para.print('xi_c.shape = {0}'.format(str(xi_c.shape))) # debug

            # Add the last column
            xi_c_ = sp.concatenate((xi[:, 0 : int(nocc)], xi_c), axis = 1)

            Af_list, msg = quick_det(xi_c_, ener = fscf.obf.eigval,
                                     fix_v1 = True, maxfn = user_input.maxfn,
                                     I_thr = user_input.I_thr,
                                     e_lo_thr = user_input.ELOW, e_hi_thr = user_input.EHIGH, 
                                     comm = para.pool.comm, 
                                     zeta_analysis = user_input.zeta_analysis and ik == 0)

            col = 2 + ixyz

            for order, Af in enumerate(Af_list):

                stick = Af_to_stick(Af)
                ener_axis, spec = stick_to_spectrum(stick, user_input)

                # important information for understanding shakeup effects and convergence 
                para.print("order {0:>2}: no. of sticks = {1:>7}, max stick = {2} ".
                            format( order + 1, len(stick), max([s[1] for s in stick] + [0.0]) ))

                ener_axis += user_input.ESHIFT_FINAL + fscf.obf.eigval[int(nocc)]

                if first:
                    spec_xas_ = sp.zeros([len(ener_axis), 4 + 1])
                    spec_xas_[:, 0] = ener_axis
                    first = False

                spec_xas_[:, col] += spec

            para.print()
        # end of ixyz

        spec_xas_[:, 1] = spec_xas_[:, 2] + spec_xas_[:, 3] + spec_xas_[:, 4]

        # output for debug
        postfix = ''
        if not user_input.gamma_only:
            postfix += '_ik{0}'.format(ik)
        if nspin == 2:
            postfix += '_ispin{0}'.format(ispin)
        postfix += '.dat'
        
        sp.savetxt(spec_xps_fname + postfix, spec_xps_, delimiter = ' ')
        spec_xps.append(spec_xps_)

        sp.savetxt(spec_xas_fname + postfix, spec_xas_, delimiter = ' ')
        spec_xas.append(spec_xas_)

    # end if spec0_only
# end of isk

## Output Spectra

# intial-state one-body
spec0_i[:, 1 : 1 + nspin] /= 3.0
if ismpi() and para.pool.isroot():
    spec0_i[:, 1 : ] = para.pool.rootcomm.reduce(spec0_i[:, 1 : ], op = MPI.SUM)

if para.isroot(): sp.savetxt(spec0_i_fname, spec0_i, delimiter = ' ')

# final-state one-body
if user_input.final_1p:
    spec0_f[:, 1 : 1 + nspin] /= 3.0
    if ismpi() and para.pool.isroot():
        spec0_f[:, 1 : ] = para.pool.rootcomm.reduce(spec0_f[:, 1 : ], op = MPI.SUM) 
    if para.isroot(): sp.savetxt(spec0_f_fname, spec0_f, delimiter = ' ')

para.stop() # debug

## Calculate total many-body spectra 

# convolute spin-up and -down spectra if nspin == 2
# *** Is this too cumbersome ?
if pool.isroot():

    spec_xas_cvlt = [None] * nspin
    spec_xps_cvlt = [None] * nspin

    first = True
    # go over all the isk th elements
    for isk in range(pool.sk_list_maxl):

        if nspin == 2:
            # find the xps that needs to be sent for the isk tuples overall all pools
            for pool_i_recv, skl in enumerate(pool.sk_list_all):
                # if wanted xps (of opposite spin) is on this pool
                if isk < len(skl): # if skl has the isk th element
                    twin_sk = (1 - skl[isk][0], skl[isk][1]) # find its twin
                    if twin_sk in pool.sk_list:
                        ind = pool.sk_list.index(twin_sk)
                        if pool.rootcomm and pool_i_recv != pool.i:
                            pool.rootcomm.send(spec_xps[ind], dest = pool_i_recv)
                        else:
                            spec_xps_twin = spec_xps[ind]
            # receive xps_twin if it is not on the same pool
            if isk < len(pool.sk_list) and (1 - pool.sk_list[isk][0], pool.sk_list[isk][1]) not in pool.sk_list[isk]:
                spec_xps_twin = pool.rootcomm.recv(source = MPI.ANY_SOURCE)

        if isk < len(pool.sk_list):

            # convolute xas with xps_twin
            spec_cvlt = convolute_spec(spec_xas[isk], spec_xps_twin) if nspin == 2 else spec_xas[isk].copy()
            ispin = pool.sk_list[isk][0]
            if first:
                spec_xas_cvlt[ispin] = spec_cvlt
            else:
                spec_xas_cvlt[ispin][:, 1 :: ] += spec_cvlt[:, 1 :: ]

            # convolute xps with xps_twin
            spec_cvlt = convolute_spec(spec_xps[isk], spec_xps_twin) if nspin == 2 else spec_xps[isk].copy()
            if first:
                spec_xps_cvlt[ispin] = spec_cvlt
            else:
                spec_xps_cvlt[ispin][:, 1 :: ] += spec_cvlt[:, 1 :: ]
                first = False

            # *** convolute major sticks (gamma_only)

        if pool.rootcomm: pool.rootcomm.barrier()

# Add spectra from each k-point
spec_xas_final = sp.zeros([spec_xas_cvlt[0].shape[0], 4 * nspin + 1])
spec_xas_final[:, 0] = spec_xas_cvlt[0][:, 0]
for ispin in range(nspin):
    if pool.rootcomm:
        spec_xas_final[:, 1 + ispin :: nspin] = pool.rootcomm.reduce(spec_xas_cvlt[ispin][:, 1 :: ], op = MPI.SUM)
    if not ismpi():
        spec_xas_final[:, 1 + ispin :: nspin] = spec_xas_cvlt[ispin][:, 1 ::]

spec_xps_final = sp.zeros([spec_xps_cvlt[0].shape[0], 2])
spec_xps_final[:, 0] = spec_xps_cvlt[0][:, 0]
# the two spin channels are the same for xps
if pool.rootcomm:
    spec_xps_final[:, 1] = pool.rootcomm.reduce(spec_xps_cvlt[ispin][:, 1], op = MPI.SUM)
if not ismpi():
    spec_xps_final[:, 1] = spec_xps_cvlt[ispin][:, 1]

# the end
postfix = '.dat'
sp.savetxt(spec_xas_fname + postfix, spec_xas_final, delimiter = ' ')
sp.savetxt(spec_xps_fname + postfix, spec_xps_final, delimiter = ' ')

# Bye ! ~
