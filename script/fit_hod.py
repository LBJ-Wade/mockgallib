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

Notation:
    *domain* is a tuble of (region, z_min, z_max) such as ('w1', '1.0', '1.2')
    *data* is defined for each domain which contains:
       halo_lightcones
       random_lightcones
       galaxy_catalogues
       random_catalogues
"""



import os
import sys
import argparse
import json
import signal
import numpy as np
import scipy.optimize
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


mock.set_loglevel(2)

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
# The main data structure defined for each region x redshift bin
#

class Data:
    """Data is defined on each domain (region, z_min, z_max)
    """
    def __init__(self, domain):
        self.domain = domain
        self.reg = domain[0]
        self.z_min = float(domain[1])
        self.z_max = float(domain[2])

        # Load lightcones
        self.halo_lightcones = mock.LightCones()
        self.halo_lightcones.load_h5(
                           ['halo_lightcone/lightcone_%05d.h5' % (n + 1)
                            for n in range(int(arg.nmocks))])

        self.rand_lightcones = mock.LightCones()
        self.rand_lightcones.load_h5(
                           ['rand_lightcone/lightcone_%05d.h5' % (n + 1)
                            for n in range(int(arg.nrands))])

        # Catalogues will be generated from lightcones for given hod
        self.galaxy_catalogues = mock.Catalogues()
        self.random_catalogues = mock.Catalogues()

        self.corr = mock.CorrelationFunction(
            rp_min=0.5, rp_max=60.0, nbin=20, pi_max= 60.0, pi_nbin= 20)


        # VIPERS projected correlation function
        self.wp_obs = np.loadtxt('data/vipers/try1/corr_projected_%s_%s_%s.txt'
                                 % domain, delimiter=' ')
        

    def generate_catalogues(self, hod):
        """Generate galaxy and random catalogues from given HOD"""
        self.galaxy_catalogues.generate_galaxies(hod, self.halo_lightcones,
                                                 self.z_min, self.z_max)
        
        self.random_catalogues.generate_randoms(hod, self.rand_lightcones,
                                                self.z_min, self.z_max)


    def chi2(self, hod):
        """Compute chi2 between HOD and observation 
        projected correlation functions wp's.
        Assumed that both wps are computed in the same coditions.
        """
        self.generate_catalogues(hod)
        
        self.corr.compute_corr_projected(self.galaxy_catalogues,
                                         self.random_catalogues)

        wp = self.corr.wp
        diff = (self.corr.wp - self.wp_obs[:,1])/self.wp_obs[:,2]
        chi2 = np.sum(diff**2)

        return chi2


    def write_corr_projected(self, index):
        """Write projected correlation function"""
        arg = self.domain + (index,)
        print(arg)
        filename = 'log/corr_%s_%s_%s_%05d.txt' % arg
        rp = self.corr.rp
        wp = self.corr.wp
        with open(filename, 'w') as f:
            for i in range(len(rp)):
                f.write('%e %e\n' % (rp[i], wp[i]))


        
regs = ['w1', 'w4']


domains = []
for zbin in redshift_bins:
    domains.append(('w1', zbin[0], zbin[1]))
data = {}

for domain in domains:
    data[domain] = Data(domain)



#
# Set HOD parameters (initial guess)
#
hod = mock.Hod()
hod.set_coef([12, 0.0, 0.0, 0, 0.1, 0.0, 15.0, 0.0, 1.5, 0.0])


#
# Setup HOD parameter fitting
#
nbar= mock.NbarFitting(hod, nbar_obs, 0.6, 1.2)



#
# Setup output
#
outdir = 'log'
if not os.path.exists(outdir):
    os.mkdir(outdir)

def write_nbar_fitting(nz, index):
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
    filename = '%s/nz_%05d.txt' % (outdir, index)
    with open(filename, 'w') as f:
        for i in range(len(nz)):
                f.write('%e %e %e\n' % (nz.z[i], nz.nbar_obs[i], nz.nbar_hod[i]))


flog = open('%s/fit_hod.log' % outdir, 'w')
iter = 0

def test(x):
    return (x[0] - 13.0)**2 + (x[1] - 0.2)**2 + (x[2] - 1.0)**2
    
def cost_function(x):
    """Compute the chi^2
    Args:
        domains: An arrain of domains
        domain is a tuble containing the survey region x redshift bin info
          domain[0] (str): 'w1' or 'w4'
          domain[1] (str): redshift_min
          domain[2] (str): redshift_max
    """

    hod[4] = x[0]
    hod[6] = x[1]
    hod[8] = x[2]

    # Find best fitting logMmin(z) function
    nbar.fit()

    chi2 = 0
    for domain, d in data.items():
        print('computing chi2 from domain', domain)
        chi2 += d.chi2(hod)

    #print('chi2 = %.3f | %.3f %.3f %.3f' % (chi2, x[0], x[1], x[2]))


    return chi2

def logging_minimization(x):
    global iter
    iter = iter + 1

    print(x)
    print('callback called')
    hod[4] = x[0]
    hod[6] = x[1]
    hod[8] = x[2]

    # Find best fitting logMmin(z) function
    nbar.fit()
    write_nbar_fitting(nbar, iter)

    chi2 = 0
    for domain, d in data.items():
        chi2 += d.chi2(hod)
        d.write_corr_projected(iter)

    print('chi2 = %.3f | %.3f %.3f %.3f' % (chi2, x[0], x[1], x[2]))
    #flog.write('%.3 %.4f %.4f %.4f\n' % (chi2, x[0], x[1], x[2]))

    return None

#
# Chi2 minimization
#
# x = [logM1, sigma, alpha]  HOD parameters to be optimised
#
x0 = [13.0, 0.1, 1.5] # starting point
ss = [1.5, 0.05, 0.5] 

#opt = scipy.optimize.minimize(cost_function, x0, method='Nelder-Mead',
#                              tol=0.01,
#                              callback=logging_minimization)

#x = mock.minimise(cost_function, None, x0, ss)
x = mock.minimise(cost_function, logging_minimization, x0, ss)
#x = mock.minimise(test, logging_minimization, x0, ss)
print('minimum', x)

flog.close()
#write_nbar_fitting(nbar, iter)
#write_corr_projected(corr, zbin, iter)

