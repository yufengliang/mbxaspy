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
import numbers

if len(sys.argv) <= 1:
    print('usage: {} spec1.dat [w1] spec2.dat [w2] ...'.format(os.path.basename(__file__)))
    print('If provided, wi (equal to 1.0 by default) will be recognized as the weight of the proceeding spectrum. ')
    print('All the weights will be normalized in the end.')
    sys.exit(0)

w = 0 # total weight
i = 1
nargv = len(sys.argv)
spec_all = None
while i < nargv:
    arg = sys.argv[i]
    if arg.split('.')[-1] == 'dat' : # then this is recognized as a spectrum file
        f = arg
        wi = 1.0
        if i + 1 < nargv:
            i += 1
            try:
                wi = float(arg[i])
            except ValueError:
                i += 1
                continue
        if wi < 1:
            print('non-positive weight {} reset to 1.0.'.format(wi))
            wi = 1.0
        w += wi
    else:
        i += 1
        continue
    print('{} {}'.format(f, wi))
    i += 1
    spec = sp.loadtxt(f)
    if spec_all is None:
        spec_all = spec * wi
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
        spec_new[:, j] = f_all(spec_new[:, 0]) + f(spec_new[:, 0]) * wi

    spec_all = spec_new

# take average
spec_all[:, 1 : ] /= w

sp.savetxt('spec_all.dat', spec_all)
        
        
    
