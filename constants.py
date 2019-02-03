
# File names
# input of shirley_xas calculation
iptblk_fname        = 'Input_Block.in'
tmp_iptblk_fname    = 'TMP_INPUT.in'
overlap_fname       = 'overlap.dat'

spec0_i_fname       = 'spec0_i.dat'
spec0_f_fname       = 'spec0_f.dat'
spec_xps_fname      = 'spec_xps'
spec_xas_fname      = 'spec_xas'
spec_bse_fname      = 'spec_bse.dat'
spec_o_fname        = 'spec_o.dat'
spec_afi_fname       = 'spec_afi.dat'

# acceptable data: name, size (bytes), and format in pack/unpack
data_set = {'integer' : (4, 'i'), 'float' : (4, 'f') , 'double' : (8, 'd'), 'complex' : (16, 'dd')}

# number of transition operators (dipole + quadrupole)
nxyz = 9

# labels of polarization directions
pol_label = {-2 : 'User-defined vector (EVEC)', -1 : 'Angular average', 0 : 'x', 1 : 'y', 2 : 'z'}
pol2num = {'evec': -2, 'all' : -1, 'x' : 0, 'y' : 1, 'z' : 2}
# pol2num is for mapping input polarization string into polarization code:
# e.g., ixyz_list = 'evec all x y z' --> [-2, -1, 0, 1, 2]

# the reverse function
def num2pol(n):
    for p in pol2num:
        if pol2num[p] == n: return p
    return 'undefined'

nbnd_max = 10000

# threshold
small_thr = 1e-8
zero_thr = 1e-16

# separation line
main_sepl = '=' * 90
second_sepl = '-' * 90

# constants related to spectra
ms_const, ns_const      = 30, 6     # define the corner of the zeta matrix for printing
det_thr_print           = 0.3       # % threshold for printing out determinants  
det_thr_label           = 0.8       # % threshold for labeling determinants  

# physics constants
Ryd = 13.605698065894

# other
elem_maxl = 5
