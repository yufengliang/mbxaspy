""" a module for input and output """

from __future__ import print_function


from struct import pack, unpack
import re
import sys
import inspect
import heapq


from constants import *
from utils import *


_quote = {'"', "'"}
_delimiter = {';', ',', ' ', '\t', '\n'}


def input_from_binary(fhandle, data_type, ndata, offset):
    """ input data from a binary file 

    Args:
        fhandle: file handle. The file needs to be opened first.
        data_type: 'float', 'double', 'complex'.
        ndata: length of the data measured in data_type
        offset: start to read at the offset measured in data_type.
                count from the head of fhandle.
        
    Returns:
        a list of data of specified data_type
    """

    if not data_type in data_set:
        raise TypeError(' data_type must be in ' + str(set(data_set)) + '.' )

    # current file position
    pos = fhandle.tell() 

    # seek and read
    fhandle.seek(offset * data_set[data_type][0])
    #data = fhandle.read(ndata * data_set[data_type][0])
    data = fhandle.read(ndata * data_set[data_type][0])

    # convert
    data_len = data_set[data_type][0] * ndata
    reslist = list(unpack(data_set[data_type][1] * ndata, data[0 : data_len]))
    if data_type == 'complex':
        # Add odd and even-index elements into complexes
        reslist = [re + 1j * im for re, im in zip(reslist[::2], reslist[1::2])]
    # return to previous position
    fhandle.seek(pos)

    return reslist


def input_arguments(lines, lower = False):
    """ input arguments from a user-defined file

    Given lines = 

        "
    	nbnd_f = 300 # number of final-state orbitals
	# is_gamma = 
	
	nbnd_i=400; job_done = .true., is_absorption =FALSE
	ecut = 30.0
	"

    Should return a dictionary like:
	{'nbnd_f' : '300', 'nbnd_i' : '400', 'job_done' : '.true.', 'is_absorption' : 'FALSE', 'ecut' : '30.0'}

    Args:
        fname: file name

    Returns:
        a dictionary as above
    """

    var_dict = {}
    if len(lines) == 0: return var_dict
    if lines[-1] != '\n': lines += '\n'
    lines = re.sub('#.*\n', '#', lines) # convert all comments into _delimiters
    for block in lines.split('#'):
        name = None
        # look for assignment
        for item in block.split('='):
            value = None
            new_name = item
            # if value is string
            for s in _quote:
                item_str = item.split(s)
                if len(item_str) > 2: # found quotation marks
                    value = item_str[1] # the string in the first _quote
                    new_name = item_str[-1].strip() # last term
            # value not a string
            if value is None:
                value = item
                for s in _delimiter:
                    try:
                        value = list(filter(None, value.split(s)))[0] # always take the first meaningful string
                    except IndexError:
                        value = ''
                        break
                for s in _delimiter:
                    try:
                        new_name = list(filter(None, new_name.split(s)))[-1] # always take the last meaningful string
                    except IndexError:
                        new_name = ''
                        break
            if is_valid_variable_name(name) and value is not None:
                if lower: name = name.lower()
                var_dict.update({name : value})
            name = new_name
    return var_dict


def convert_val(val_str, val):
    """ Given a string, convert into the correct data type """
    if val is bool:
        if 'true' in val_str.lower(): val_str = 'true' 
        else: val_str = '' # otherwise set to false
    val_type = val
    try:
        return val_type(val_str)
    except ValueError:
        # Can it be a float ?
        return val_type(float(val_str))


def list2str_1d(nums, mid = -1):
    """ 
    Give a list of nums, output the head, the middle, and the tail of it
    with nice format. Return the formatted string.

    Args:
    mid: define the middle point you are interested in
    """

    nvis = 3 # numbers printed out in each part
    l = len(nums) 
    mid = mid if mid > 0 else l / 2
    fmtstr = '{0:.4f} '
    resstr = ''
    for i in range(min(nvis, l)):
        resstr += fmtstr.format(nums[i])
    irange = range(max(nvis + 1, int(mid - nvis / 2)), min(l, int(mid + nvis / 2 + 1)))
    if len(irange) > 0: resstr += ' ... '
    for i in irange:
        resstr += fmtstr.format(nums[i])
    irange = range(max(int(mid + nvis / 2 + 1), l - nvis), l)
    if len(irange) > 0: resstr += ' ... '
    for i in irange:
        resstr += fmtstr.format(nums[i])
    return resstr


def eigvec2str(eigvec, m, n, nctr, nvis = 6, npc = 6, iws = '  '):
    """
    Output some prominent matrix elements for an eigenvector matrix
    
    eigvec is given as a 1D array:
    [ <B_1|1k>, <B_1|2k>, <B_2|1k>, <B_2|2k>]

    which corresponds to such a matrix (m rows x n cols):
    <B_1|1k>  <B_1|2k>
    <B_2|1k>  <B_2|2k>

    nctr: list states around the center nctr
    nvis: number of printed out states
    npc:  number of principal components
    iws:  initial white spaces for indentation

    """

    resstr = iws + '{0:<10}{1}\n'.format('norm', 'principal components')
    # guarantee len(eigvec) = n * m
    for j in range(max(0, nctr - int(nvis / 2) + 1), min(n, nctr + int(nvis / 2) + 1)):
        eabs = [ abs(eigvec[i * n + j]) ** 2 for i in range(m) ]
        norm = sum(eabs)
        # print out the norm
        resstr += iws + '{0:<10.5f}'.format(norm)
        # the wavefunction to print: |nk> = ...
        resstr += '|n={0:>4},k> = '.format(j)
        # Find the npc elements with the largest norm
        indices = heapq.nlargest(npc, range(len(eabs)), key = lambda i : eabs[i])
        for i in indices:
            resstr += '({0:11.3f})|B_{1}> + '.format(eigvec[i * n + j], i)
        # resstr = resstr[:-3] # delete the last +
        resstr += '... \n'
    return resstr


def atomic_species_to_list(asp_str):
    """
    Convert a atomic_species block (as in Qespresso) into a list like:

    [['Br', 'Br.pbe-van_mit.UPF'], ['C', 'C.pbe-van_bm.UPF'], ... ]

    """
    res = []
    for l in asp_str.split('\n'):
        words = l.split()
        if len(words) == 3 and len(words[0]) < 3 and words[2].split('.')[-1] == 'UPF': 
            # *** There should be more robust sanity checks: check elements
            res.append([words[0], words[2]])
    return res


def atomic_positions_to_list(apos_str):
    """
     Convert a atomic_positions block (as in Qespresso) into a list like:

    [['Pb', '0.0', '0.0', '0.0'], ['Br', '0.0', '0.0', '0.5'], ...]

    Most interested in the atoms' names rather than their positions 
    """
    res = []
    for l in apos_str.split('\n'):
        words = l.split()
        if len(words) >= 4 and len(words[0]) < 3:
            # *** There should be more robust sanity checks: check elements
            res.append(words)
    return res


def read_qij_from_upf(upf_fname):
    """
    Given a PAW/ultrasoft pseudopotential in UPF format, find the projectors' angular momenta
    and the corresponding Q_int matrices in the file.

    """
    l   = [] # angular momentum number
    qij = [] # Q_int matrix
    i, j = 0, 0
    errmsg = ''
    fh = []
    try:
        fh = open(upf_fname, 'r')
    except IOError:
        errmsg = 'cannot open UPF file: ' + str(upf_fname)
    for line in fh:
        words = line.split()
        if len(words) >= 4 and words[2 : 4] == ['Beta', 'L']:
            l.append(int(words[1]))
        if len(words) >= 2 and words[1] == 'Q_int':
            if (i, j) == (0, 0):
                # if first time, initialize the qij matrix
                qij = [[0.0] * len(l) for _ in l]
            try:
                qij[i][j] = float(words[0])
            except IndexError:
                errmsg = 'too many Q_int for given l = ' + str(len(l))
                break
            qij[j][i] = qij[i][j]
            j += 1
            if j > len(l) - 1: i, j = i + 1, i + 1
    if fh: fh.close()
    return l, qij, errmsg

def get_index(s):
    """ 
    get the index in the string like 'a[4]'.
    should return 4
    """
    return int(s[s.find("[")+1:s.find("]")])

# export function only
__all__ = [s for s in dir() if not s.startswith('_') and inspect.isfunction(getattr(sys.modules[__name__],s))]


if __name__ == '__main__':
    print(__file__ + ": the i/o module for mbxaspy")
    # debug
    var_dict = input_arguments(sys.stdin.read())
    print(var_dict)
