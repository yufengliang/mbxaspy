"""
a standalone script for:
combining calculated spectra from mbxaspy for all excited atoms

usage:
__file__ spec1.dat spec2.dat ... 
__file__ means this file.
output as spec_all.dat

(1) number of spec?.dat can vary.

(2) The first column of spec?.dat is the energy axis and the rest are spectral column.
The energy axis and number of columns of all files must match the first one, otherwise
the file will be skipped.


"""

import sys, os
import scipy as sp
from scipy.interpolate import interp1d

if len(sys.argv) <= 1:
    print('usage: {} spec1.dat spec2.dat ...'.format(os.path.basename(__file__)))
    sys.exit(0)

for i, f in enumerate(sys.argv[1 : ]):
    spec = sp.loadtxt(f)
    if i == 0:
        spec_all = spec
        continue

    # energy axis check
    #if not sp.allclose(spec_all[:, 0], spec[:, 0], atol = 1e-6):
    #    print('skipping {} ... (energy mismatch)'.format(f))
    #    continue
    
    if spec.shape[1] != spec_all.shape[1]:
        print('skipping {} ... (column mismatch)'.format(f))
        continue

    ## doing interpolation for different energy axes

    # to cover both axes
    emin = min(spec[0, 0], spec_all[0, 0])
    emax = max(spec[-1, 0], spec_all[-1, 0])


    # find out the size of the energy grid, assuming it is uniform
    de = min(spec[1, 0] - spec[0, 0], spec_all[1, 0] - spec_all[0, 0])
    ne = int((emax - emin) / de) + 2
    de = (emax - emin) / ne
    spec_new = sp.zeros((ne, spec_all.shape[1]))
    spec_new[:, 0] = [_ * de + emin for _ in range(ne)]

    for j in range(1, spec_all.shape[1]):
        f_all = interp1d(spec_all[:, 0], spec_all[:, j], bounds_error = False, fill_value = 0.0)
        f = interp1d(spec[:, 0], spec[:, j], bounds_error = False, fill_value = 0.0)
        spec_new[:, j] = f_all(spec_new[:, 0]) + f(spec_new[:, 0])

    spec_all = spec_new

# take average
spec_all[:, 1 : ] /= len(sys.argv) - 1

sp.savetxt('spec_all.dat', spec_all)
        
        
    
