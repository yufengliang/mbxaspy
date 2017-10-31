""" 

rixs f1 main 

Calculate 2D RIXS map to the f(1) order
Modified from mbxaspy (minimally)

"""

from __future__ import print_function

import sys
import os

from utils import *
from defs import *
from init import *
from spectra import *
from analysis import *
from bse import *
from rixs import *

## input user-defined arguments from stdin
para.print(' Running RIXS calculatoin at f(1) order ...')
userin.read()
userin.maxfn = 1 # hardcoded as 1

## Check initial and final state and perform sanity checks
# Check the initial-state scf calculations
para.sep_line()
para.print(' Checking initial-state scf from: \n ' + userin.path_i + '\n')
iscf.input(isk = -1, nelec = userin.nelec)

# Check the final-state scf calculations
para.sep_line()
para.print(' Checking final-state scf from: \n ' + userin.path_f + '\n')
fscf.input(is_initial = False, isk = -1, nelec = iscf.nelec)

if userin.scf_type == 'shirley_xas':
    # Important: Need to tell iscf the index of the core in the xmat file with multiple atoms
    iscf.proj.icore = fscf.proj.icore
    # or tell iscf which atom is excited (for use_pos = True)
    iscf.proj.x     = fscf.proj.x
    # Input < B_i | \tilde{B}_j >
    fscf.obf.input_overlap(userin.path_f, iscf.nbnd, fscf.nbnd)

nspin = iscf.nspin


## initialize spectral information
if not userin.spec0_only:

    from xi import *
    # Compute full atomic overlap sij
    if userin.scf_type == 'shirley_xas': 
        compute_full_sij(fscf.proj)

    from determinants import *
    spec_xps_all = []   
    sticks_xps_all = []
    # spec_xas_all = [] 
    rixs_maps_all = []

global_ener_axis = spec_class(userin).ener_axis
# if para.isroot(): sp.savetxt('ener_axis.dat', global_ener_axis) # debug

# global energy shift
global_offset = userin.ESHIFT_FINAL + fscf.e_lowest

# initialize i and f spec0
def init_spec(nspin = 1):
    return [spec_class(ener_axis = global_ener_axis) for s in range(nspin)]
spec0_i = init_spec(nspin)
if userin.final_1p: spec0_f = init_spec(nspin)
if userin.want_bse: spec_bse = init_spec(nspin)

if userin.spec_analysis:
    # perform analysis on spectra at Gamma-point only
    def init_order(): return [[None] * 2 for order in range(userin.maxfn)]
    spec_xps_g,     spec_xas_g      = init_order(), init_order()
    sticks_xps_g,   sticks_xas_g    = init_order(), init_order()


## setup ixyz list *** need more works

ixyz_list = [-1, 0, 1, 2] # userin.ixyz_list
# if not completed
ixyz_list_ = ixyz_list[:]
if -1 in ixyz_list:  # need x, y, z to do "-1"
    ixyz_list_ += [ixyz for ixyz in [0, 1, 2] if ixyz not in ixyz_list ]
if userin.xps_only: ixyz_list_ = []

# RIXS polarization
inout_pols = userin.inout_pol.split()

# output format
fn_fmt      = 'f^(n)/statistics\n!{:>8} {:>10}  {:>12}  {:>12}'
fn_num_fmt  =                   '!{:>8} {:>10}  {:>12.7}  {:>12.7}' 

# charge transfer
qi = [0] * len(ixyz_list_)
qf = [0] * len(ixyz_list_)

## loop over spin and kpoints
for isk in range(pool.nsk):

    ispin, ik  = pool.sk_list[isk] # acquire current spin

    para.sep_line()
    para.print(' Processing (ispin, ik) = ({0},{1}) \n'.format(ispin, ik), flush = True)

    # weight the sticks according to the k-grid
    weight = iscf.kpt.weight[ik]; 
    prefac = weight * Ryd
    # para.print('weight = {0}'.format(weight)) # debug

    # Import the initial-state scf calculation
    para.sep_line(second_sepl)
    para.print(' Importing initial-state scf\n', flush = True)
    iscf.input(isk = isk)
    para.print('  xmat: the {0}th atom in ATOMIC_POSITIONS is excited.'.format(xatom(iscf.proj, iscf.xmat) + 1))

    # Import the final-state scf calculation
    para.sep_line(second_sepl)
    para.print(' Importing final-state scf\n', flush = True)
    fscf.input(is_initial = False, isk = isk)
    if userin.final_1p: para.print('  xmat: the {0}th atom in ATOMIC_POSITIONS card is excited.'.format(xatom(fscf.proj, fscf.xmat) + 1))

    # Obtain the effective occupation number: respect the initial-state #electrons
    if nspin == 1: nocc = iscf.nocc
    else: nocc = iscf.nocc[ik][ispin]

    ## Compute non-interacting spectra *** should I put it in a def ?
    para.print('  Calculating one-body spectra ...\n')

    # sticks = xmat_to_sticks(iscf, [-2], nocc, evec = [1.0, 0.0, 0.0]) # debug
    # print(sticks[0]) debug

    # initial-states: pec0_i.dat
    sticks = xmat_to_sticks(iscf, ixyz_list_, nocc, offset = -fscf.e_lowest, evec = userin.EVEC)
    spec0_i[ispin].add_sticks(sticks, userin, prefac, mode = 'additive')
    spec0_i_os_sum = os_sum(sticks)
    # sp.savetxt('spec0_sticks.dat', sp.array(sticks), delimiter = ' ', fmt = '%s')# debug

    # final-state: spec0_f.dat
    if userin.final_1p:
        sticks = xmat_to_sticks(fscf, ixyz_list_, nocc, offset = -fscf.e_lowest, evec = userin.EVEC)
        spec0_f[ispin].add_sticks(sticks, userin, prefac, mode = 'additive')
    
    para.print('  One-body spectra finished.', flush = True)

    ## Calculate charge transfer
    q_isk = calc_occ_pdos(iscf, ixyz_list_, nocc, evec = userin.EVEC)
    qi = [qi[_] + q_isk[_] for _ in range(len(ixyz_list_))]
    q_isk = calc_occ_pdos(fscf, ixyz_list_, nocc, evec = userin.EVEC)
    qf = [qf[_] + q_isk[_] for _ in range(len(ixyz_list_))]

    ## Compute many-body spectra
    if not userin.spec0_only:

        para.print('  Calculating many-body spectra ... ')

        ## Compute the transformation matrix xi
        para.print('  Calculating transformation matrix xi ... ')
        xi = compute_xi(iscf, fscf)

        size = min(xi.shape)
        xi_eigvals = la.eigvals(xi[: size, : size])
        para.print('  Average of the eigenvalues of xi: {}'.format( sp.sum(abs(xi_eigvals)) / size ))

        if userin.xi_analysis and para.isroot() and ik == 0:
            plot_xi(xi) # debug
            if nspin > 1:
                msg = eig_analysis_xi(xi, '_spin_{0}'.format(ispin)) # debug
            else:
                msg = eig_analysis_xi(xi) # debug

        para.print('  Matrix xi finished. ', flush = True)

        ## BSE spectra
        if userin.want_bse:
            sticks = bse(xi, iscf, fscf, nocc, ixyz_list, offset = -fscf.e_lowest, evec = userin.EVEC)
            spec_bse[ispin].add_sticks(sticks, userin, prefac, mode = 'additive')
            # sp.savetxt('bse_sticks.dat', sp.array(sticks), delimiter = ' ', fmt = '%s')# debug

        ## XPS spectra (N X N)
        para.sep_line(second_sepl)
        para.print('  Calculating many-body XPS spectra ... ')

        Af_list, msg = quick_det(xi[:, 0 : int(nocc)], ener = fscf.obf.eigval,
                                 fix_v1 = False, maxfn = userin.maxfn - 1,
                                 I_thr = userin.I_thr,
                                 e_lo_thr = userin.ELOW, e_hi_thr = userin.EHIGH, 
                                 comm = pool.comm, 
                                 zeta_analysis = userin.zeta_analysis and ik == 0)

        spec_xps_isk = init_spec()[0]

        sticks_all_order = []
        para.print(fn_fmt.format('order', '#sticks', 'stick max', 'os sum'))
        for order, Af in enumerate(Af_list):

            sticks = Af_to_sticks(Af)
            #sticks_all_order += sticks
            sticks_all_order.append(sticks)

            # important information for understanding shakeup effects and convergence
            if len(sticks) > 0:
                para.print(fn_num_fmt.format( order, len(sticks), max([s[2] for s in sticks]), os_sum(sticks)[0]))

            # don't use prefac here
            spec_xps_isk.add_sticks(sticks, userin, mode = 'additive')

            if userin.spec_analysis and ik == 0:
                spec_xps_g[order][ispin] = copy.deepcopy(spec_xps_isk)
                sticks_xps_g[order][ispin] = sticks

        spec_xps_all.append(spec_xps_isk)
        sticks_xps_all.append(sticks_all_order)

        # output for debug
        postfix = '_ik{0}'.format(ik)
        if nspin == 2:
            postfix += '_ispin{0}'.format(ispin)
        postfix += '.dat'

        spec_xps_isk.savetxt(spec_xps_fname + postfix, offset = 0.0)
        
        para.print('  XPS spectra finished. ', flush = True)
        para.print()

        if userin.xps_only: continue

        ## 2D RIXS map
        para.sep_line(second_sepl)
        para.print('  Calculating RIXS map ... ')

        rixs_maps = [] # for different combinations of polarizations
        
        for pol in inout_pols:
            
            ixyz_in = pol_index[pol[0]]
            ixyz_out = pol_index[pol[1]]

            rixs_maps.append(

            rixs_f1(xi = xi, nelec = int(nocc), xmat_in = iscf.xmat[:, 0, ixyz_in], xmat_out = iscf.xmat[:, 0, ixyz_out], 
                    ener_i = iscf.obf.eigval, ener_f = fscf.obf.eigval + fscf.e_lowest,
                    nbnd_i = iscf.nbnd_use, nbnd_f = fscf.nbnd_use,
                    e_in_lo = userin.ELOW, e_in_hi = userin.EHIGH, nener_in = userin.NENER,
                    loss_mode = userin.loss_mode, eloss_range = userin.eloss_range, nener_out = userin.NENER_out,
                    Gamma_h = userin.SIGMA, Gamma_f = userin.SIGMA_out,
                    I_thr = userin.RIXS_I_thr, 
                    comm = None)
            )

        rixs_maps_all.append(rixs_maps)

        para.print('  RIXS map finished. \n', flush = True)
    # end if spec0_only
# end of isk

## Output one-body spectra

para.print('Output one-body spectra ...', flush = True)

def mix_spin(spec):
    # mix spin up and down
    return spec[0] if nspin == 1 else spec[0] | spec[1]

# intial-state one-body
para.print('mp_summing ...', flush = True)
for ispin in range(nspin):
    print('ispin = ', ispin, 'rank = ', para.rank, 'rootcomm = ', pool.rootcomm)#debug
    sys.stdout.flush()
    spec0_i[ispin].mp_sum(pool.rootcomm) 
para.print('mp_summed ...', flush = True)
spec0_i = mix_spin(spec0_i)
if para.isroot(): spec0_i.savetxt(spec0_i_fname, offset = global_offset)

# final-state one-body
if userin.final_1p:
    for ispin in range(nspin): spec0_f[ispin].mp_sum(pool.rootcomm) 
    spec0_f = mix_spin(spec0_f)
    if para.isroot(): spec0_f.savetxt(spec0_f_fname, offset = global_offset)

# spec0_sum = spec0_i[0] + spec0_f[0] # test operator overload
# spec0_sum.savetxt('spec0_sum.dat')

para.print('one-body spectra output to files.\n', flush = True)

# output charge-transfer
qi, qf = sp.array(qi), sp.array(qf)
if pool.rootcomm and pool.rootcomm != MPI.COMM_NULL:
    qi = pool.rootcomm.allreduce(qi, op = MPI.SUM)
    qf = pool.rootcomm.allreduce(qf, op = MPI.SUM)
para.print(' Charge transfer analysis: ')
para.print(' initial-state')
para.print([pol_label[ixyz_list_[_]] + ': {:12.7}'.format(qi[_]) for _ in range(len(ixyz_list_))])
para.print(' final-state')
para.print([pol_label[ixyz_list_[_]] + ': {:12.7}'.format(qf[_]) for _ in range(len(ixyz_list_))])

if userin.spec0_only:
    para.done() # debug

## Output BSE spectra
if userin.want_bse:
    for ispin in range(nspin): spec_bse[ispin].mp_sum(pool.rootcomm) 
    spec_bse = mix_spin(spec_bse)
    if para.isroot(): spec_bse.savetxt(spec_bse_fname, offset = global_offset)

## Calculate total many-body spectra 

# convolute spin-up and -down spectra if nspin == 2
if nspin == 2:

    if not userin.xps_only:
        para.print('Calculating total RIXS map for nspin = 2.', flush = True)

        # convolute RIXS spectra: do this before xps
        for isk, sk in enumerate(pool.sk_list):
            ispin, ik = sk
            if (1 - ispin, ik) not in pool.sk_list:
                pool.log('The twin tuple ({0}, {1}) for ({2},{3}) is not on this pool with {5}'
                        .format(ispin, ik, 1 - ispin, ik, str(pool.sk_list)))
                pool.log('Unsupported distribution of sk-tuples', flush = True)
                para.exit()
            ind = pool.sk_list.index((1 - ispin, ik))
            xps_weight = sticks_xps_all[ind][0][1] # [0]: zero-order, [1]: for intensity (0 for energy)
            for rixs_map in rixs_maps:
                rixs_map *= xps_weight
        para.print(' RIXS map convoluted with xps.', flush = True)

    # convolute xps spectra
    isk_done = []
    # Add spectra from each k-point
    for isk, sk in enumerate(pool.sk_list):
        if isk in isk_done: continue
        ispin, ik = sk
        ind = pool.sk_list.index((1 - ispin, ik)) # don't need to check existence again
        sticks_xps_twin = []
        for order in range(userin.maxfn): sticks_xps_twin += sticks_xps_all[ind][order]
        spec_xps_all[isk] *= sticks_xps_twin
        # the two spin channels are the same for xps
        spec_xps_all[ind] = spec_xps_all[isk] # note that the twins are correlated now (use deepcopy to uncorrelate them)
        isk_done += [isk, ind]
    para.print(' many-body XPS convoluted.', flush = True)

    # convolute gamma-point spectra for analysis
    if userin.spec_analysis:
        for ispin in range(nspin):
            ind = pool.sk_list.index((1 - ispin, 0))
            sticks_xps_twin = []
            for order in range(userin.maxfn):
                sticks_xps_twin += sticks_xps_all[ind][order]
                #spec_xas_g[order][ispin] *= sticks_xps_g[order][1 - ispin]
                if ispin == 0: spec_xps_g[order][ispin] *= sticks_xps_g[order][1 - ispin]
    para.print(' many-body XPS at gamma convoluted.\n', flush = True)

# intrapool summation

spec_xps = init_spec()[0]
if not userin.xps_only:
    # rixs_maps[spin][polarization]
    rixs_maps = [[rixs_map_class(rixs_maps_all[0][0]) for pol in inout_pols] for s in range(nspin)]

for isk, sk in enumerate(pool.sk_list):
    ispin, ik = sk
    weight = iscf.kpt.weight[ik]    
    if ispin == 0: spec_xps += spec_xps_all[isk] # two spin channels are the same
    if not userin.xps_only:
        for ip in range(len(inout_pols)):
            rixs_maps[ispin][ip] += rixs_maps_all[isk][ip]
spec_xps *= weight

# mpi reduce
para.print('Collecting spectra at all k-points...', flush = True)
spec_xps.mp_sum(pool.rootcomm)
if not userin.xps_only:
    for ispin in range(nspin):
        for ip in range(len(inout_pols)):
            rixs_maps[ispin][ip].mp_sum(pool.rootcomm)
para.print('Spectra at all k-points collected.\n', flush = True)

# This requires the world root is also one of the pool roots: can be made more robust
if para.isroot():
    postfix = '.dat'
    spec_xps.savetxt(spec_xps_fname + postfix)
    if not userin.xps_only:
        for ispin in range(nspin):
            for ip, pol in enumerate(inout_pols):
                if userin.loss_mode: rixs_maps[ispin][ip].I = rixs_maps[ispin][ip].I.T
                rixs_fname_tmp = rixs_fname + '.ispin{}.{}'.format(ispin, pol) + postfix
                rixs_maps[ispin][ip].savetxt(rixs_fname_tmp)

## Convolute the initial-state spectrum with XPS: test orthogonality
if userin.want_spec_o:
    spec_o = spec0_i * spec_xps
    if para.isroot():
        spec_o.savetxt(spec_o_fname, offset = global_offset) 

if userin.spec_analysis and para.isroot():
    for order in range(userin.maxfn):
        if not userin.xps_only:
            pass
            # mix_spin(spec_xas_g[order]).savetxt(spec_xas_fname + '_maxfn_{}.dat'.format(order + 1), offset = global_offset)
        spec_xps_g[order][0].savetxt(spec_xps_fname + '_maxfn_{}.dat'.format(order))

# Bye ! ~
para.done()
