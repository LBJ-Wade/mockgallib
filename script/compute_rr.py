"""
This script computes projected correlation function from given files of
galalxy and random catalgues

python3 compute_corr_projected.py

Args:
    n: index of lightcone

Options:
    -i=1

Input:
   ../rands/rand_<reg>_<n>.txt
    
    File format is:
    Column 1-3: x y z (starting from 0)
    Column 6-7: RA Dec
    in each line
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
parser.add_argument('-i', default='1',
                    help='index range of random catalogues')
parser.add_argument('-o', default='', help='output RR filename')
parser.add_argument('--param', default='param.json',
                    help='parameter json file')
parser.add_argument('--zmin', type=float, default=0.5, help='minimum redshift')
parser.add_argument('--zmax', type=float,  default=1.2, help='minimum redshift')

arg = parser.parse_args()

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

#
# Read catalogues
#
randoms = mock.Catalogues()
filename = '../rands/%s/rand_%s_%05d.txt' % (arg.reg, arg.reg, int(arg.i))
a = np.loadtxt(filename, delimiter=' ', usecols=[1,2,3,6,7])
randoms.append(a, z_min=arg.zmin, z_max= arg.zmax)


#
# Compute RR paris
#
corr = mock.CorrelationFunction(rp_min=0.1, rp_max=60.0, nbin=24,
                                pi_max=60.0, pi_nbin=20,
                                ra_min=0.001388889, dec_min=0.0375)


rr = mock.corr.Hist2D(rp_min=0.1, rp_max=60.0, rp_nbin=24, 
                      pi_max=60.0, pi_nbin=20)


npairs = corr.compute_corr_projected_rr(randoms, rr)

#
# Write
#
f = h5py.File(arg.o, 'w')
f['npairs'] = npairs
f['rr'] = rr[:]
f.close()
