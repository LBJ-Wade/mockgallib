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
    rands/random_<n>.txt
"""

import os
import argparse
import json
import signal
import numpy as np
import mockgallib as mock

signal.signal(signal.SIGINT, signal.SIG_DFL) # stop with ctrl-c

#
# Command-line options
#
parser = argparse.ArgumentParser()
parser.add_argument('n', help='index of lightcone')
parser.add_argument('--reg', default='w1', help='region w1/w4')
parser.add_argument('--dir', default='.', help='base data directory') 
parser.add_argument('--param', default='param.json',
                    help='parameter json file')
parser.add_argument('--mock', help='generate mock catalogue',
                    action="store_true")
parser.add_argument('--rand', help='generate random catalogue',
                    action="store_true")

arg = parser.parse_args()

data_dir = '/workplace/wp2e/como5/data'

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
mock.power.init(arg.dir + '/' + param['power_spectrum'])

#
# redshift range
#
z_min = 0.39
z_max = 1.21
print('redshift range %f %f' % (z_min, z_max))

# nz
# nbar_obs=  mock.array.loadtxt(arg.dir + '/' + param['nz'])

# sky
sky = {}
for reg in param['reg']:
    sky[reg['name']] = mock.Sky(reg['ra'], reg['dec'], [z_min, z_max])


#
# Set HOD parameters
#
hod = mock.Hod()
hod_param = [11.632682100874081, -0.5706390738948128, 4.904043697780981, -1.0126352684312565, 0.45, 0.9, 1.05, 0.0, 0.9, 0.0, 4.0, 2.0]

hod.set_coef(hod_param)

#nbar= mock.NbarFitting(hod, nbar_obs, 0.6, 1.2)
#x0 = [1.0, 0.15, 1.0]
#nbar.fit()
    
lightcones = mock.LightCones()
cats = mock.Catalogues()

n = int(arg.n)

def write_catalogue(filename, a):
    with open(filename, 'w') as f:
        for i in range(a.shape[0]):
            f.write('%d %e %e %e %e %e %e %e %e %e %e\n' % (
                    i,
                    a[i,0], a[i,1], a[i,2], 
                    a[i,4], a[i,3],
                    a[i,5], a[i,6],
                    a[i,7], a[i,10], a[i,11]))

reg = arg.reg

# mock
if arg.mock:
    halo_lightcones = mock.LightCones()
    halo_lightcones.load_h5(
        ['%s/halo_lightcone/%s/lightcone_%05d.h5' % (arg.dir, reg, n)])
    galaxy_catalogues = mock.Catalogues()
    galaxy_catalogues.generate_galaxies(hod, halo_lightcones, sky[reg],
                                        z_min, z_max)

    write_catalogue('mocks/%s/mock_%s_%05d.txt' % (reg, reg, n), 
                    galaxy_catalogues[0])

if arg.rand:
    rand_lightcones = mock.LightCones()
    rand_lightcones.load_h5(
        ['%s/rand_lightcone/%s/lightcone_%05d.h5' % (arg.dir, reg, n)])

    random_catalogues = mock.Catalogues()
    random_catalogues.generate_randoms(hod,  rand_lightcones, sky[reg],
                                       z_min, z_max)

    write_catalogue('rands/%s/rand_%s_%05d.txt' % (reg, reg, n), 
                    random_catalogues[0])

# Column 0: index
# Column 1: x realspace [1/h Mpc]
# Column 2: y
# Column 3: z
# Column 4: vr [km/s]
# Column 5: redshift
# Column 6: RA
# Column 7: Dec
# Column 8: M_host_halo
# Column 9: r_satellite
# Column 10: vr_satellite
