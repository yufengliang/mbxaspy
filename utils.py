""" Utilities """

from __future__ import print_function

import sys
import os
from ast import parse

def is_valid_variable_name(name):
    """test if name is a valid python variable name"""
    try:
        parse('{} = None'.format(name))
        return True
    except (SyntaxError, ValueError, TypeError) as err:
        return False

def ispython3x():
    if sys.version_info.major == 3: return True
    else: return False

def isanaconda():
    if 'anaconda' in sys.version.lower(): return True
    else: return False

def ismpi():
    mpi_var = {'SLURM_JOB_ID', 'PBS_O_WORKDIR'}
    if set(mpi_var) & set(os.environ): return True
    else: return False

def valid_comm(comm):
    """ test if comm has been initialized and is a valid communicator"""
    if 'MPI' not in dir(): return False
    return comm and comm != MPI.COMM_NULL

def find_nocc(two_arr, n):
    """
    Given two sorted arrays of the SAME lengths and a number,
    find the nth smallest number a_n and use two indices to indicate
    the numbers that are no larger than a_n.

    n can be real. Take the floor.
    """
    l = len(two_arr[0])
    if n >= 2 * l: return l, l 
    if n == 0: return 0, 0
    res, n = n % 1, int(n)
    lo, hi = max(0, n - l - 1), min(l - 1, n - 1)
    while lo <= hi:
        mid = int((lo + hi) / 2) # image mid is the right answer
        if mid + 1 < l and n - mid - 2 >= 0:
            if two_arr[0][mid + 1] < two_arr[1][n - mid - 2]:
                lo = mid + 1
                continue
        if n - mid - 1 < l:
            if two_arr[1][n - mid - 1] < two_arr[0][mid]:
                hi = mid
                continue
        break
    if n - mid - 1 >= l or mid + 1 < l and two_arr[0][mid + 1] < two_arr[1][n - mid - 1]:
        return mid + res + 1, n - mid - 1
    else:
        return mid + 1, n - mid - 1 + res
          
