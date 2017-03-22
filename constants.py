
# input of shirley_xas calculation
iptblk_fname        = 'Input_Block.in'
tmp_iptblk_fname    = 'TMP_INPUT.in'
overlap_fname       = 'overlap.dat'

# acceptable data: name, size (bytes), and format in pack/unpack
data_set = {'integer' : (4, 'i'), 'float' : (4, 'f') , 'double' : (8, 'd'), 'complex' : (16, 'dd')}

# number of transition operators (dipole + quadrupole)
nxyz = 9

# separation line
separation_line = "===================================================================================================="
