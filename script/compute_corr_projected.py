"""
This script computes projected correlation function from given files of
galalxy and random catalgues

python3 compute_corr_projected.py

Args:
    n: index of lightcone

Options:
    --igalaxies=1:1
    --irandoms=1.1

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
import numpy as np
import mockgallib as mock


#
# Command-line options
#
parser = argparse.ArgumentParser()
parser.add_argument('--igalaxies', default='1:1',
                    help='index range of galaxy catalogues')
parser.add_argument('--irandoms', default='1:1',
                    help='index range of random catalogues')
parser.add_argument('--param', default='param.json',
                    help='parameter json file')
arg = parser.parse_args()

igalaxies = arg.igalaxies.split(':')
irandoms = arg.irandoms.split(':')

#
# Read parameter file
#
print('Parameter file: %s' % arg.param)

with open(arg.param, 'r') as f:
    param = json.load(f)

omega_m = param['omega_m']
print('Setting cosmology: omega_m= %.4f' % omega_m)


#
# Initilise
#
mock.set_loglevel(0)
mock.cosmology.set(omega_m)

def read_catalogues(filebase, irange):
    cats = mock.Catalogues()
    for i in range(int(irange[0]), int(irange[1]) + 1):
        filename = '%s%05d.txt' % (filebase, i)
        a = np.loadtxt(filename, delimiter=' ')[:,0:3]
        cats.append(a)
    return cats

galaxies = read_catalogues('mocks/mock_', igalaxies)
randoms  = read_catalogues('randoms/random_', igalaxies)

corr = mock.CorrelationFunction()
a = corr.compute_corr_projected(galaxies, randoms)

#n = a.shape[0]
#for i in range(n):
#    print('%e %e %e' % (a[i,0], a[i,1], a[i,2]))


