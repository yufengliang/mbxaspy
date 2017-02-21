""" a module for input and output """

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
        reslist = [re + 1j * im for re, im in zip(reslist[::2], reslist[1::2])]
    # return to previous position
    fhandle.seek(pos)

    return reslist

