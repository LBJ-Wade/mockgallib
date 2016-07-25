""" lightcone generator

This script generates lightcone

python3 lightcone.py [--random] --ibegin=6001 --iend=6001 w1

Options:
    --param [=param.json]: parameter file
    --random:              generate lightcone with randomised points

Input:
    <fof_dir>/<abc>/fof<isnp><abc>.b

Output:
    halo_lightcone/lightcone_<n>.h5
    rand_lightcone/lightcone_<n>.h5 for --random
"""

import sys
import os
import argparse
import json
import numpy as np
import mockgallib as mock


def redshift_from_a(a):
    return 1.0/a - 1.0


#
# Command-line options
#
parser = argparse.ArgumentParser()
parser.add_argument('reg', help='region w1 or w4')
parser.add_argument('--param', default='param.json',
                    help='parameter json file')
parser.add_argument('--ibegin', type=int, default=1,
                    help='first index of output')
parser.add_argument('--iend', type=int, default=1, help='last index of output')
parser.add_argument('--random', help='generate random lightcone',
                    action="store_true")
arg = parser.parse_args()


#
# Read parameter file
#

print('Parameter file: %s' % arg.param)

with open(arg.param, 'r') as f:
    param = json.load(f)

# find region name w1 or w4
for reg in param['reg']:
    if reg['name'] == arg.reg:
        break

if reg['name'] != arg.reg:
    raise Exception('region %s not found in parameter file %s' %
                    (arg.reg, arg.param))

print('Region: %s' % reg['name'])

mock.set_loglevel(0)

# omega_m
omega_m = param['omega_m']
print('Setting cosmology: omega_m= %.4f' % omega_m)
mock.cosmology.set(omega_m)

# power spectrum
mock.power.init(param['power_spectrum'])
mock.sigma.init()

# Random
if arg.random:
    random = True
    out_dir = 'rand_lightcone'
    print('random')
else:
    random = False
    out_dir = 'halo_lightcone'
    print('halo lightcone')

# index range
isnp0 = reg['index'][0]

# redshift range
z_mins = []
z_maxs = []

for snp in param['snapshots']:
    z_mins.append(redshift_from_a(snp['a'][2]))
    z_maxs.append(redshift_from_a(snp['a'][1]))

z_min = min(z_mins)
z_max = max(z_maxs)

print('Lightcone redshift range: z= (%.5f, %.5f)' % (z_min, z_max))

# sky
print('Sky ra=(%.5f, %.5f), dec=(%.5f, %.5f)' % tuple(reg['ra'] + reg['dec']))
sky = mock.Sky(reg['ra'], reg['dec'], [z_min, z_max])

# box remapping
print('Remap u=' + reg['remap'].__repr__())
remap = mock.Remap(reg['remap'], param['boxsize'])

# slice
slice = mock.Slice(remap, sky)

nslice = len(slice)

# lightcones
lightcones = mock.LightCones()

# snapshots
snapshots = mock.Snapshots()
fof_dir = param['fof']
halomass_dir = param['halomass']

# output directory
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

for isnp in range(arg.ibegin, arg.iend+1):
    print(isnp)

    snapshots.clear()
    lightcones.clear()

    for snp in param['snapshots']:
        abc = snp['abc']
        filename_fof = '%s/%s/fof%05d%s.b' % (fof_dir, abc, isnp, abc)
        filename_halo_mass = '%s/halomass_%s.txt' % (halomass_dir, abc)

        snapshots.insert(filename_fof, filename_halo_mass, snp['a'])

    lightcones.create_from_snapshots(snapshots, sky, remap, slice, random)

    print('Write lightcones')
    for islice, lightcone in enumerate(lightcones):
        iout = nslice*(isnp - isnp0) + islice + 1
        filename = '%s/lightcone_%05d.h5' % (out_dir, iout)
        print(filename)
        lightcone.save_h5(filename)
