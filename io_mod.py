""" a module for input and output """

from __future__ import print_function


from struct import pack, unpack
import re
import sys
import inspect


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


def list2str_1d(nums):
    """ 
    Give a list of nums, output the head, the middle, and the tail of it
    with nice format. Return the formatted string.
    """

    nvis = 3 # numbers printed out in each part
    l = len(nums)
    fmtstr = '{0:.4f} '
    resstr = ''
    for i in range(min(nvis, l)):
        resstr += fmtstr.format(nums[i])
    irange = range(max(nvis + 1, int(l / 2 - nvis / 2)), min(l, int(l / 2 + nvis / 2 + 1)))
    if len(irange) > 0: resstr += ' ... '
    for i in irange:
        resstr += fmtstr.format(nums[i])
    irange = range(max(int(l / 2 + nvis / 2 + 1), l - nvis), l)
    if len(irange) > 0: resstr += ' ... '
    for i in irange:
        resstr += fmtstr.format(nums[i])
    return resstr

# export function only
__all__ = [s for s in dir() if not s.startswith('_') and inspect.isfunction(getattr(sys.modules[__name__],s))]


if __name__ == '__main__':
    print(__file__ + ": the i/o module for mbxaspy")
    # debug
    var_dict = input_arguments(sys.stdin.read())
    print(var_dict)
