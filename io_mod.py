""" a module for input and output """

from __future__ import print_function
from ast import parse
from struct import pack, unpack

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

    # acceptable data: name, size (bytes), and format in pack/unpack
    data_set = {'integer' : (4, 'i'), 'float' : (4, 'f') , 'double' : (8, 'd'), 'complex' : (16, 'dd')}

    if not data_type in data_set:
        raise TypeError(' data_type must be in ' + str(set(data_set)) + '.' )

    # current file position
    pos = fhandle.tell() 

    # seek and read
    fhandle.seek(offset * data_set[data_type][0])
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


def is_valid_variable_name(name):
    """test if name is a valid python variable name"""
    try:
        parse('{} = None'.format(name))
        return True
    except (SyntaxError, ValueError, TypeError) as err:
        return False


def input_arguments(fname):
    """ input arguments from a user-defined file

    Given:

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
    with open(fname, 'r') as fh:
        # read lines and only keep the non-comment part
        for line in fh:
            for block in line.split('#')[0].split(';'):
                for item in block.split(','):
                    # look for assignment
                    items = item.split('=')
                    if len(items) > 1:
                        name, value = items[0].strip(), items[1].strip()
                        if is_valid_variable_name(name):
                            var_dict.update({name : value})
    return var_dict

if __name__ == '__main__':
    print(__file__ + ": the i/o module for mbxaspy")
    # debug
    var_dict = input_arguments('test/test_io/sample.in')
    print(var_dict)
