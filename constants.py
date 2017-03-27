
# File names
# input of shirley_xas calculation
iptblk_fname        = 'Input_Block.in'
tmp_iptblk_fname    = 'TMP_INPUT.in'
overlap_fname       = 'overlap.dat'

spec0_i_fname       = 'spec0_i.dat'
spec0_f_fname       = 'spec0_f.dat'

# acceptable data: name, size (bytes), and format in pack/unpack
data_set = {'integer' : (4, 'i'), 'float' : (4, 'f') , 'double' : (8, 'd'), 'complex' : (16, 'dd')}

# number of transition operators (dipole + quadrupole)
nxyz = 9

# threshold
small_thr = 1e-8
zero_thr = 1e-16

# separation line *** you can use [=] * length
separation_line = "===================================================================================================="

# physics constants
Ryd = 13.605698065894
