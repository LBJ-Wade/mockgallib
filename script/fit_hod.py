"""
This script find best fitting Halo Occcupation Distribution (HOD) parameters
that simultaneously fit the number density n(z) and projected correlation
function.

python3 fit_hod.py param.json

Options:
    param [=param.json]: parameter file
    nmocks [=1]: number of mock catalogues used for corr_projected
    nrands [=1]: number of random catalogues used for corr_projected

Input: 
    halo_lightcone/lightcone_<imock>.h5;  imock = 1 - nmocks
    rand_lightcone/lightcone_<irand>.h5;  irand = 1 - nrands   

Output:
    mocks/mock_<n>.txt
    randoms/random_<n>.txt
"""

import os
import sys
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
parser.add_argument('--param', default='param.json',
                    help='parameter json file')
parser.add_argument('--nmocks', default='1', help='number of mock catalogues')
parser.add_argument('--nrands', default='1', help='number of random catalogues')

arg = parser.parse_args()


mock.set_loglevel(0)

#
# Read parameter file and initialise modules
#
#
print('Parameter file: %s' % arg.param)

with open(arg.param, 'r') as f:
    param = json.load(f)

# omega_m
omega_m = param['omega_m']
print('Setting cosmology: omega_m= %.4f' % omega_m)
mock.cosmology.set(omega_m)

# power_spectrum
print('Using linear power spectrum: ', param['power_spectrum'])
mock.power.init(param['power_spectrum'])

# redshift_bins
def read_redshift_bins(redshift_bins):
    arr = []
    for zbin in redshift_bins:
        arr.append((float(zbin['zmin']), float(zbin['zmax'])))

    return arr

redshift_bins = read_redshift_bins(param['redshift_bins'])
print(redshift_bins)

# nz
nbar_obs= np.loadtxt(param['nz'], delimiter=' ')



#
# Load lightcones
#
regs = ['w1', 'w4']

halo_lightcones = mock.LightCones()
halo_lightcones.load_h5(['halo_lightcone/lightcone_%05d.h5' % (n + 1)
                         for n in range(int(arg.nmocks))])

rand_lightcones = mock.LightCones()
rand_lightcones.load_h5(['rand_lightcone/lightcone_%05d.h5' % (n + 1)
                         for n in range(int(arg.nrands))])


#
# Set HOD parameters (initial guess)
#
hod = mock.Hod()
hod.set_coef([12, 0.0, 0.0, 0, 0.1, 0.0, 15.0, 0.0, 1.5, 0.0])


#
# Setup HOD parameter fitting
#
nbar= mock.NbarFitting(hod, nbar_obs, 0.6, 1.2)
corr = mock.CorrelationFunction()


#
# Setup output
#
outdir = 'fit'
if not os.path.exists(outdir):
    os.mkdir(outdir)

def write_nbar_fitting(nz, iter):
    """Write nbar_fitting to an ascii file
    Args:
        nz: NbarFitting object
        iter: index of output filename

    Output:
        outdir/nbar_<iter>.txt file
        Column 1: z
        Column 2: nbar_obs
        Column 3: nbar_HOD
    """
    filename = '%s/nz_%05d.txt' % (outdir, iter)
    with open(filename, 'w') as f:
        for i in range(len(nz)):
            f.write('%e %e %e\n' % (nz.z[i], nz.nbar_obs[i], nz.nbar_hod[i]))


def write_corr_projected(corr, zbin, iter):
    """Write projected correlation function
    """
    filename = '%s/corr_%.2f_%.2f_%05d.txt' % (odir, zbin[0], zbin[1], iter)
    rp = corr.rp
    wp = corr.wp
    with open(filename, 'w') as f:
        for i in range(len(rp)):
            f.write('%e %e\n' % (rp[i], wp[i]))


mock_catalogues = mock.Catalogues()
random_catalogues = mock.Catalogues()

iter = 1
nbar.fit()
write_nbar_fitting(nbar, iter)

sys.exit()

for zbin in redshift_bins:
    random_catalogues.generate_randoms(hod, rand_lightcones, zbin[0], zbin[1])
    galaxy_catalogues.generate_galaxies(hod, halo_lightcones, zbin[0], zbin[1])
    corr.compute_corr_projected(galaxy_catalogues, random_catalogues)


    write_corr_projected(corr, zbin, iter)

