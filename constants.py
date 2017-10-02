
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

# acceptable data: name, size (bytes), and format in pack/unpack
data_set = {'integer' : (4, 'i'), 'float' : (4, 'f') , 'double' : (8, 'd'), 'complex' : (16, 'dd')}

# number of transition operators (dipole + quadrupole)
nxyz = 9

# labels of polarization directions
pol_label = {-2 : 'user evec', -1 : 'angular average', 0 : 'x', 1 : 'y', 2 : 'z'}

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
