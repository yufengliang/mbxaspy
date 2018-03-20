""" Find all nontrivial determinants """

from __future__ import print_function


import sys
import os
import bisect

from constants import *
from utils import *
from init import *


def quick_det(xi_mat, ener, fix_ind1 = True, I_thr = 1e-3, maxfn = 2, 
              e_lo_thr = -8.0, e_hi_thr = 2.0,
              det_scaling_fac = 1.0,
              xes = False,
              comm = None,
              zeta_analysis = False):
    """
    Given a m x n (m > n) matrix, find all of its significant n x n minors
    This procedure employs an efficient algorithm that makes use of the multi-linearity
    of determinants.

    Examples:
    A is 5 x 3 matrix and A contains 5 row vectors (a1, a2, a3, a4, a5)
    Then its 3 x 3 minors can be expressed as ai x aj x ak, where i, j, k
    are 3 distinct integers chosen from 1 .. 5. 'x' stands for wedge (outer) product 
    of two vectors.

    Assume (a1, a2, a3) forms a full-rank matrix. Then a4 and a5 can be expressed as linear
    superposition of the 3 vectors:

    a4 = sum_{i=1..3} zeta_{4, i} ai
    a5 = sum_{i=1..3} zeta_{5, i} ai

    If a1 x a2 x a3 is known, then a1 x a2 x a4 can be obtained simply by using the expansion
    coefficient above:

    a1 x a2 x a4 = a1 x a2 x (zeta_{41} a1 + zeta_{42} a2 + zeta_{43} a3) = zeta_{43}

    This is simply the last coefficient.

    This code first obtains this zeta-matrix and then search for all its minors
    with one column being the last one, starting from the size 1 x 1, 2 x 2, 3 x 3, and so forth

    xes: if doing XES calculations, then the above column / row convention will be different compared to XAS/XPS
    
    """

    Af_list = []
    msg = ''
    if xes and fix_ind1:
        n, m = xi_mat.shape # n < m in XES
    else:
        m, n = xi_mat.shape

    if m < n:
        msg = 'not enough no. of orbitals '
        return Af_list, msg

    ## Calculate the mother determinant: topmost n x n minor

    det_ref = la.det(xi_mat[ : n,  : n])
    para.print('det_ref = {0}\n'.format(det_ref)) # debug
    xi_mat_tmp = xi_mat[ : n,  : n]

    if not fix_ind1: # if doing xps then, then add the f^(0) term
        Af_list.append({'' : sp.array([0.0, det_ref])})

    """
    If the reference determinant is too small, then replace the last row vector.
    A small reference determinant may be caused by a weak first transition. The
    idea is to replace it with a higher-energy and bright transition.

    This works if we believe the f^(1) group always has very bright transitions. 
    """
    # *** THIS MAY NOT ALWAYS WORK *** / not for XES now
    if not xes and abs(det_ref) < small_thr:

        xi_mat_q, xi_mat_r = la.qr(xi_mat.T)
        xi_mat_tmp[n - 1] = xi_mat_q[:, n - 1].T
        det_ref = la.det(xi_mat_tmp)

    # Construct the zeta-matrix
    xi_mat_inv = la.inv(xi_mat_tmp)

    ## Now the zeta matrix !

    if xes and fix_ind1:
        
        """
        zeta matrix for XES:
        
        | < 1~ | o | h >   xi_{1, N - 1}  ... xi_{1, 1} |       |       |       xi_{1, nbnd_use}    ... xi_{1, N}
        |                                               |       |       |
        | < 2~ | o | h >   xi_{2, N - 1}  ... xi_{2, 1} |       |       |       xi_{2, nbnd_use}    ... xi_{2, N}
        |                                               |   x   | zeta  | =
        |    ...                                        |       |       |
        |                                               |       |       |
        | < N~ | o | h >   xi_{N, N - 1}  ... xi_{N, 1} |       |       |       xi_{N, nbnd_use}    ... xi_{N, N}
        """
        
        zeta_mat = sp.matrix(xi_mat_inv) * sp.matrix(xi_mat[:, n : ])

        # overwrite the first line because transformation coeff. related to < i~ | o | h > are not needed.
        zeta_mat[0, :] = 0.0
        zeta_mat[0, -1] = 1.0
        
    else:
        zeta_mat = sp.matrix(xi_mat[n - 1 : m, :]) * sp.matrix(xi_mat_inv)

    # para.print(zeta_mat) # debug

    # print out a corner of zeta_mat
    ms, ns = min(ms_const, m - n + 1), min(ns_const, n)
    para.print('| zeta matrix | (left-upper corner): ')
    para.print('{0:5}'.format('') + ' '.join(['{0:>12}'.format(_) for _ in range(zeta_mat.shape[1] - ns, zeta_mat.shape[1])]))
    for im in range(ms):
        para.print('{0:5}'.format(im) + ' '.join(['{0:>12.3e}'.format(abs(_)) for _ in sp.array(zeta_mat[im, - ns : ])[0]]))
    para.print()
    
    ## Check the sparsity of the zeta-matrix
    
    if zeta_analysis:
        plot_zeta(zeta_mat)

    """
    Carry out a breath-first search for nontrivial determinants.
     
    For each configuration that has been obtained and stored in the queue, the search procedure 
    matches one significant matrix element in the zeta matrix and calculates the contributions of 
    this contribution to a child configuration. The matrix element is matched in an ascending row order
    so as to avoid double-counting, which means the new matrix element can only be chosen
    from the row vectors below the rows that form the parent configuration.

    """

    ## Now starting the breadth-first search process
	 
    # Format of Af: {configuration string: [energy, amplitude/intensity]}

    Af = {}
    max_conf = ''
    Af[max_conf] = sp.array([0.0, det_ref]) # This is considered as the f^(0) configuration

    # estimate the overall determinant theshold according to the maximum of the lowest-order contribution
    max_zeta = abs(zeta_mat[:, -1]).max() if fix_ind1 else abs(zeta_mat).max()
    para.print('max_zeta = {0}'.format(max_zeta))
    f_lowest_max = max_zeta * abs(det_ref)
    para.print('filter out transitions with intensity < {0:>6.3}% of the maximal intensity {1}'.format(I_thr * 100, abs(f_lowest_max) ** 2))
    det_thr = abs(I_thr) * f_lowest_max # you should not take the square root of I_thr because | dI / I | ~ | 2 dA / A |
    para.print('determinant threshold = {0}'.format(det_thr))
    para.print()

    # parallelism
    if comm:
        rank, size = comm.Get_rank(), comm.Get_size()
    else:
        rank, size = 0, 1

    for ndepth in range(1, maxfn + 1):

        """
        Find out major matrix elements of the zeta matrix.
        Record the coordinates and values of all significant matrix elements
        in such an order
        0 0 4 0 1 
        0 6 5 3 2
        0 7 0 0 0
        which means we need to flip our zeta_mat from left to right, and then tranpose it:
        1 2 0
        0 3 0
        4 5 0
        0 6 7
        0 0 0
        """

        # estimate sparse threshold 'sparse_thr' dynamically according to previous max Af
        if len(Af) == 0 or abs(Af[max_conf][1]) < zero_thr:
            Af = {}
            Af_list.append(Af)
            para.print('f^({0}) contributions: '.format(ndepth))
            para.print('None')
            continue

        sparse_thr = det_thr / abs(Af[max_conf][1])

        zeta_coord = sp.where(abs(sp.fliplr(zeta_mat).T) > sparse_thr)
        zeta_coord = zip(zeta_coord[0], zeta_coord[1])

        # Now convert the coordinates back
        # i = j', j = n_col - 1 - i'
        zeta_coord = [(j, zeta_mat.shape[1] - 1 - i) for i, j in zeta_coord]

        # para.print(zeta_coord) # debug
        para.print('No. of signiciant matrix elements of zeta: {0}'.format(len(zeta_coord)))
        
	# Record Af from all possible final-state indices f
	# 'If' is an array of dictionaries (len = nspin)
        Af_new = {}

        """
        configuration string:

        XPS
        'v1 c1 v2 c2 v3 c3 ... '

        XAS
        'c1 v2 c2 v3 c3 ... ' (v1 = n - 1)

        XES
        'v1 c2 v2 c3 v3 ... ' (c1 = n - 1)

        that satisfies:
        v1 > v2 > v3 > ..., and c1 < c2 < c3 < ...

	When fix_ind1 = True, v1(c1) is defaulted to (n - 1) and is omitted from conf.
        
        fix_ind1 = True,      for XAS / XES
                   False,     for XPS

        """

        for conf in Af:

            # Basic information of the parent configuration

            # para.print(conf)
				
            conf_ = [int(_) for _ in conf.split()]

            """
            indexing of zeta for XES (n rows x m cols )

                    ... N + 2   N + 1   N
                    ______________________
            N       |
            N - 1   |
            N - 2   |
            ...     |
            2       |
            1       |
                
            """
            if xes and fix_ind1:
                # indices of empty (c) states
                conf_c = conf_[1 :: 2]
                # the maximum of the chosen c states
                conf_maxc = max(conf_c) if ndepth > 1 else n - 2
                # indices of occupied (v) states; conf_v is guaranteed to be sorted
                conf_v = conf_[ :: 2]
                conf_v_set = set(conf_v)
            else:
                # indices of occupied (v) states
                conf_v = conf_[int(fix_ind1) :: 2]
                # the minimum of the chosen v states
                if fix_ind1 and ndepth > 1 : conf_minv = min(conf_v + [n - 1])
                else: conf_minv = min(conf_v + [n])
                # indices of empty (c) states ; conf_c is guaranteed to be sorted
                conf_c = conf_[1 - int(fix_ind1) :: 2]
                conf_c_set = set(conf_c)

            # energy of the parent configuration
            f_ener = Af[conf][0]

            ## peform the breath-first search in descending column order

            # Make sure it is in descending column order
            if xes and fix_ind1:
                low_izeta = bisect.bisect_left([-j for i, j in zeta_coord], -(m - 2 - (conf_maxc + 1)))
            else:
                low_izeta = bisect.bisect_left([-j for i, j in zeta_coord], -(conf_minv - 1))

            # Loop over all nontrivial matrix elements of zeta
            # Distribute over the given communicator comm
            # *** You may also consider distribute conf
            for izeta in range(low_izeta + rank, len(zeta_coord), size):


                if xes and fix_ind1:

                    new_c = m - 2 - zeta_coord[izeta][1]
                    new_v = n - 1 - zeta_coord[izeta][0]

		    # ndepth = 1 is special when fix_ind1 = True
                    if ndepth == 1 and fix_ind1 and new_c > n - 1: break

                    # Make sure v doesn't appear twice in a configuration
                    if new_v in conf_v_set: continue

		    # Energy filter
	
		    # don't go too high into the conduction band
                    if ener[n - 1] - ener[new_c] + f_ener < e_low_thr: break
                    
                else:
                    
                    new_c, new_v = zeta_coord[izeta]
                    # zeta_coord[0] = 0 corresponds to c = n - 1
                    new_c += n - 1
                    
		    # ndepth = 1 is special when fix_ind1 = True
                    if ndepth == 1 and fix_ind1 and new_v < n - 1: break

		    # Make sure c doesn't appear twice in a configuration
                    if new_c in conf_c_set: continue

		    # Energy filter
	
		    # don't go too deep into the valence band
                    if ener[n - 1] - ener[new_v] + f_ener > e_hi_thr: break
                    
                if xes:
                    # the sign of XES energy is reverse !
                    enew = ener[new_v] - ener[new_c] + f_ener
                    if enew < e_low_thr: continue
                else:
                    enew = ener[new_c] - ener[new_v] + f_ener
                    if enew > e_hi_thr: continue

                """
                Find out the sign for the child configuration
                {     c1, new_v}   ...  {     c1, v_{n-1}}   ...  {     c1, v1}
                ...............................................................
                {  new_c, new_v}   ...  {  new_c, v_{n-1}}   ...  {  new_c, v1}
                ...............................................................
                {c_{n-1}, new_v}   ...  {c_{n-1},   new_v}   ...  {c_{n-1}, v1}
                sign = (-1) ** i, i is the order new_c locates in c1, c2, ..., c_{n-1}

                Example:

                ......  ......  c1, v2  c1, v1
                new_c, new_v    ......  ......
                ......  ......  c2, v2  c2, v1

                c1 < new_c < c2, so it contributes a (-1).
                
                """

                if xes and fix_ind1:

                    insert_pos =  bisect.bisect_left(conf_v, new_v) # conf_v is sorted
                    v_sgn = (-1) ** insert_pos

                    ## Construct the new configuration
                
                    # add new_c and new_v
                    conf_v_ = conf_v[slice(0, insert_pos)] + [new_v] + conf_v[slice(insert_pos, ndepth - 1)]
                    conf_c_ = list(conf_c)
                    if not fix_ind1 or ndepth > 1: conf_c_ += [new_c]
                    conf_new = [None] * (2 * ndepth - int(fix_ind1))
                    conf_new[1 - int(fix_ind1) :: 2] = conf_v_
                    conf_new[int(fix_ind1) :: 2] = conf_c_

                else:

                    insert_pos =  bisect.bisect_left(conf_c, new_c) # conf_c is sorted
                    v_sgn = (-1) ** insert_pos

                    ## Construct the new configuration
                
                    # add new_c and new_v
                    conf_c_ = conf_c[slice(0, insert_pos)] + [new_c] + conf_c[slice(insert_pos, ndepth - 1)]
                    conf_v_ = list(conf_v)
                    if not fix_ind1 or ndepth > 1: conf_v_ += [new_v]
                    conf_new = [None] * (2 * ndepth - int(fix_ind1))
                    conf_new[1 - int(fix_ind1) :: 2] = conf_c_
                    conf_new[int(fix_ind1) :: 2] = conf_v_

                conf_new = ' '.join([str(_) for _ in conf_new])

                ## Calculate the determinantal contribution to the child configuration
                new_contribution = v_sgn * Af[conf][1] * zeta_mat[zeta_coord[izeta]]

                # Too small ? throw it away !
                if abs(new_contribution) < det_thr: continue
                
                # Check if this configuration already exists
                if conf_new in Af_new:
                    Af_new[conf_new][1] += new_contribution
                else:
                    Af_new[conf_new] = sp.array([enew, new_contribution])
            # end for izeta
        # end for conf

        # reduce det_thr by a factor of the number of orbitals
        # The reason we reduce det_thr is because there are more states in higher-order f^(n)
        # We may need to lower the threshold in case they may form a background of weak transitions 
        # use this with cautions
        det_thr *= det_scaling_fac

        para.print('Breadth-first search finished.', flush = True)
        if comm:
            para.print('Gathering Af ...', flush = True)
            # Gather and combine Af_new
            Af_gather = comm.gather(Af_new, root = 0)
            # *** doing this on one core may be inefficient
            # A divide and conquer algorithm should be much faster: learn from mpi_reduce
            if rank == 0:
                Af = {}
                for iter_Af in Af_gather:
                    for conf in iter_Af:
                        if conf in Af:
                            Af[conf][1] += iter_Af[conf][1]
                        else:
                            Af[conf] = iter_Af[conf].copy()
            para.print('Finish gathering Af ...', flush = True)
        else:
            Af = Af_new.copy()
            del Af_new

        # Filter out the small terms in Af
        # *** In future, do this when combine Af
        if rank == 0:
            Af = {conf : Af[conf] for conf in Af if abs(Af[conf][1]) > det_thr}

        # This is important: avoid double-counting the f^(0) term for xps
        if rank == 0 and not fix_ind1:
            conf0 = ' '.join([str(n - 1), str(n - 1)])
            if conf0 in Af: Af.pop(conf0)

        if comm:
            para.print('Broadcasting Af ...', flush = True)
            Af = comm.bcast(Af, root = 0)
            para.print('Finish broadcasting Af.', flush = True)

        # Find out the max amplitude
        try:
            max_conf = max(Af, key = lambda conf : abs(Af[conf][1]))
        except ValueError:
            max_conf = None

        # print out major distributions
        if max_conf:
            para.print('f^({0}) major contributions: '.format(ndepth))
            first = True
            for conf in Af:
                conf_ = ' '.join(['{0:>5}'.format(_) for _ in conf.split()])
                if first:
                    para.print(' ' * (len(conf_) + 2) + '{0:>12}  {1:>12}'.format('energy (eV)', 'amplitude'))
                    first = False
                if abs(Af[conf][1]) > abs(Af[max_conf][1] * det_thr_print):
                    label = ''
                    if abs(Af[conf][1]) > abs(Af[max_conf][1] * det_thr_label): label = '*'
                    para.print('{0}: {1:>12.5}  {2:>12.5e} {3}'.format(conf_, abs(Af[conf][0]), abs(Af[conf][1]), label))

            para.print('max_conf: {0}, amp = {1}'.format(max_conf, abs(Af[max_conf][1])))
            para.print(flush = True)

        Af_list.append(Af)
    # end for ndepth

    return Af_list, msg


def plot_zeta(zeta, postfix = ''):
    """ Plot a heatmap of zeta """
    max_zeta = abs(zeta).max()
    plt.imshow(abs(zeta) / max_zeta, cmap='gray', interpolation='none')
    #plt.pcolor(sp.array(abs(zeta) / max_zeta), cmap='gray')
    plt.savefig('test_zeta{0}.png'.format(postfix), format = 'png')
    plt.close()
