"""
This script generates a mock/random catalogue from a lightcone 

python3 generate_mock.py [--random] <n>

Args:
    n: index of lightcone

Options:
    --param [=param.json]: parameter file
    --random:              generate random catalogue

Input: 
    halo_lightcone/lightcone_<n>.h5
    rand_lightcone/lightcone_<n>.h5

Output:
    mocks/mock_<n>.txt
    randoms/random_<n>.txt
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
parser.add_argument('n', help='index of lightcone')
parser.add_argument('--param', default='param.json',
                    help='parameter json file')
parser.add_argument('--random', help='generate random catalogue',
                    action="store_true")
arg = parser.parse_args()


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


#
# Set HOD parameters
#
hod = mock.Hod()
hod.set_coef([12.0, 0.5, 0, 0, 1.0, 0.0, 13.0, 0.0, 1.5, 0.0])

lightcones = mock.LightCones()
cats = mock.Catalogues()

n = int(arg.n)

if arg.random:
    filename = 'rand_lightcone/lightcone_%05d.h5' % n
    lightcones.load_h5([filename])
    cats.generate_randoms(hod, lightcones, 0.4, 1.2)
    
    outdir = 'randoms'
    ofilename = '%s/random_%05d.txt' % (outdir, n)
else:
    filename = 'halo_lightcone/lightcone_%05d.h5' % n
    lightcones.load_h5([filename])
    cats.generate_galaxies(hod, lightcones, 0.4, 1.2)

    outdir = 'mocks'
    ofilename = '%s/mock_%05d.txt' % (outdir, n)


if not os.path.exists(outdir):
    os.mkdir(outdir)

cat = cats[0]

with open(ofilename, 'w') as f:
    for i in range(cat.shape[0]):
        f.write('%e %e %e %e\n' % (cat[i,0], cat[i,1], cat[i,2], cat[i,3]))

# x,y,z,redshift
        
print(ofilename, 'written')
