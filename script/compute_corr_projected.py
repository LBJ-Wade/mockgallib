"""
This script computes projected correlation function from given files of
galalxy and random catalgues

python3 compute_corr_projected.py

Args:
    n: index of lightcone

Options:
    --igalaxies=1:1
    --irandoms=1:1

Input:
    mocks/mock_<n>.txt
    randoms/random_<n>.txt
    
    File format is:
    x y z
    in each line

Output:
    Standard output
    Column 1: rp
    Column 2: wp
    Column 3: rms(wp)  (wp if number of galaxy catalogues in 1)
"""

import os
import argparse
import json
import signal
import numpy as np
import h5py
import mockgallib as mock

signal.signal(signal.SIGINT, signal.SIG_DFL) # stop with ctrl-c


#
# Command-line options
#
parser = argparse.ArgumentParser()
parser.add_argument('reg', help='region w1 or w4')
parser.add_argument('--igalaxies', default='1:1',
                    help='index range of galaxy catalogues')
parser.add_argument('--irandoms', default='1:1',
                    help='index range of random catalogues')
parser.add_argument('--param', default='param.json',
                    help='parameter json file')
parser.add_argument('--rr-only', default=False, action='store_true',
                    help='compute RR only')
parser.add_argument('--rr', default='', help='RR filename')



parser.add_argument('--zmin', type=float, default=0.5, help='minimum redshift')
parser.add_argument('--zmax', type=float,  default=1.2, help='minimum redshift')

arg = parser.parse_args()

igalaxies = arg.igalaxies.split(':')
irandoms = arg.irandoms.split(':')

#
# Read parameter file
#
print('# Parameter file: %s' % arg.param)

with open(arg.param, 'r') as f:
    param = json.load(f)

omega_m = param['omega_m']
print('# Setting cosmology: omega_m= %.4f' % omega_m)
print('# redshift-range %f %f' % (arg.zmin, arg.zmax))

#
# Initilise
#
mock.set_loglevel(0)
mock.cosmology.set(omega_m)
mock.distance.init(1.2)

def read_catalogues(filebase, irange):
     cats = mock.Catalogues()
     for i in range(int(irange[0]), int(irange[1]) + 1):
         filename = '%s%05d.txt' % (filebase, i)
         a = np.loadtxt(filename, delimiter=' ', usecols=[1,2,3,6,7])
         cats.append(a, z_min=arg.zmin, z_max= arg.zmax)
     return cats

#galaxies = read_catalogues('../mock_%s_' % arg.reg, igalaxies)
#


corr = mock.CorrelationFunction(rp_min=0.1, rp_max=60.0, nbin=24,
                                pi_max=60.0, pi_nbin=20,
                                ra_min=0.001388889, dec_min=0.0375)

if arg.rr_only:
    randoms  = read_catalogues('../rand_%s_' % arg.reg, irandoms)
    
    rr = mock.corr.Hist2D(rp_min=0.1, rp_max=60.0, rp_nbin=24, pi_max=60.0, pi_nbin=20)

    npairs = corr.compute_corr_projected_rr(randoms, rr)

    f = h5py.File(arg.rr, 'w')
    f['npairs'] = npairs
    f['rr'] = rr[:]
    f.close()

#corr.compute_corr_projected_with_rr(galaxies, randoms, rr)
#rp = corr.rp_i(0)
#wp = corr.wp_i(0)

#print(rp)
#print(wp)

#n= len(rp)
#with open('corr.txt', 'w') as f:
#    for i in range(n):
#        f.write('%e %e\n' % (rp[i], wp[i]))

#print('corr.txt written')
