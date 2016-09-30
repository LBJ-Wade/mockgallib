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
parser.add_argument('--random', help='generate random catalogue',
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

# nz
nbar_obs=  mock.array.loadtxt(arg.dir + '/' + param['nz'])

#
# Set HOD parameters
#
hod = mock.Hod()
hod_param = [11.5659763286151, -1.3653860293894984, 7.372119497149173, -3.905197669231047, 0.4114530081890546, 0.0, 4.480542558267717, 0.0, 1.1909781896845204, 0.0]
hod.set_coef(hod_param)
nbar= mock.NbarFitting(hod, nbar_obs, 0.6, 1.2)

x0 = [1.0, 0.15, 1.0]
nbar.fit()
    
lightcones = mock.LightCones()
cats = mock.Catalogues()

n = int(arg.n)
domain = ('w1', '0.8', '1.0')
z_min = float(domain[1])
z_max = float(domain[2])

halo_lightcones = mock.LightCones()
halo_lightcones.load_h5(
    ['%s/halo_lightcone/%s/lightcone_%05d.h5' % (arg.dir, domain[0], n)])

rand_lightcones = mock.LightCones()
rand_lightcones.load_h5(
    ['%s/rand_lightcone/%s/lightcone_%05d.h5' % (arg.dir, domain[0], n)])

galaxy_catalogues = mock.Catalogues()
random_catalogues = mock.Catalogues()

wp_obs = mock.array.loadtxt(
    '%s/vipers/run4/0.1/corr_projected_%s_%s_%s.txt' % ((data_dir,) + domain))
        
covinv = np.load('%s/vipers/run4/0.1/covinv_%s_%s_%s.npy' %
                 ((data_dir,) + domain))

galaxy_catalogues.generate_galaxies(hod, halo_lightcones, z_min, z_max)
random_catalogues.generate_randoms(hod,  rand_lightcones, z_min, z_max)



with open('mock.txt', 'w') as f:
    a = galaxy_catalogues[0]
    print(a.shape)
    for i in range(a.shape[0]):
        f.write('%e %e %e %e %e %e\n' % (a[i,0], a[i,1], a[i,2], a[i,3],
                                   a[i,7], a[i,10]))

    
# Column 1: x
# Column 2: y
# Column 3: z
# Column 4: redshift
# Column 5: M
# Column 6: r satellite
                
        
print('mock.txt written')

with open('nz.txt', 'w') as f:
    for i in range(len(nbar)):
        f.write('%e %e %e\n' % (nbar.z[i], nbar.nbar_obs[i], nbar.nbar_hod[i]))
