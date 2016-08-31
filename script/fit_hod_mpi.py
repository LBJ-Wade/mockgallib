import os
import sys
import argparse
import json
import signal
import numpy as np
import mockgallib as mock

signal.signal(signal.SIGINT, signal.SIG_DFL) # stop with ctrl-c

def print0(str):
    if mock.comm.rank == 0:
        print(str)

#
# Command-line options
#
parser = argparse.ArgumentParser()
parser.add_argument('--param', default='param.json',
                    help='parameter json file')
parser.add_argument('--nmocks', default='1', help='number of mock catalogues')
parser.add_argument('--nrands', default='1', help='number of random catalogues')
parser.add_argument('--dir', default='.', help='base directory')
parser.add_argument('--outdir', default='log', help='output directory')
parser.add_argument('--loglevel', default='0', help='loglevel')

arg = parser.parse_args()
outdir = arg.outdir

mock.set_loglevel(0)

#
# Read parameter file and initialise modules
#
#
print0('Parameter file: %s' % arg.param)

param_str = ''
if mock.comm.rank == 0:
    with open(arg.param, 'r') as f:
        param_str = f.read()

param_str = mock.comm.bcast_str(param_str)

param = json.loads(param_str)
    
omega_m = param['omega_m']

# omega_m
omega_m = param['omega_m']
print0('Setting cosmology: omega_m= %.4f' % omega_m)
mock.cosmology.set(omega_m)

# power_spectrum
print0('Using linear power spectrum: ' + param['power_spectrum'])
mock.power.init(arg.dir + '/' + param['power_spectrum'])

# redshift_bins
def read_redshift_bins(redshift_bins):
    arr = []
    for zbin in redshift_bins:
        arr.append((float(zbin['zmin']), float(zbin['zmax'])))

    return arr

redshift_bins = read_redshift_bins(param['redshift_bins'])
print(redshift_bins)

# nz
nbar_obs=  mock.array.loadtxt(arg.dir + '/' + param['nz'])

flog = None
iter = 0
if mock.comm.rank == 0:
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    flog = open('%s/fit_hod.log' % outdir, 'w', 1)



#
# The main data structure defined for each region x redshift bin
#

class Data:
    """Data is defined on each domain = (region, z_min, z_max)
    """
    def __init__(self, domain):
        self.domain = domain
        self.reg = domain[0]
        self.z_min = float(domain[1])
        self.z_max = float(domain[2])

        # Load lightcones
        self.halo_lightcones = mock.LightCones()
        self.halo_lightcones.load_h5(
                           ['%s/halo_lightcone/%s/lightcone_%05d.h5' %
                            (arg.dir, domain[0], n + 1)
                            for n in range(int(arg.nmocks))
                            if n % mock.comm.n_nodes == mock.comm.rank
                           ])
        print("len lightcone", len(self.halo_lightcones))


        self.rand_lightcones = mock.LightCones()
        self.rand_lightcones.load_h5(
                           ['%s/rand_lightcone/%s/lightcone_%05d.h5' %
                            (arg.dir, domain[0], n + 1)
                            for n in range(int(arg.nrands))
                            if n % mock.comm.n_nodes == mock.comm.rank
                           ])

        # Catalogues will be generated from lightcones for given hod
        self.galaxy_catalogues = mock.Catalogues()
        self.random_catalogues = mock.Catalogues()

        self.corr = mock.CorrelationFunction(
            rp_min=0.5, rp_max=60.0, nbin=20, pi_max= 60.0, pi_nbin= 20,
            ra_min=0.001388889, dec_min=0.0375)


        # VIPERS projected correlation function
        self.wp_obs = mock.array.loadtxt(
            '%s/data/vipers/corr_projected_%s_%s_%s.txt' % ((arg.dir,) +domain))
        

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
        chi2 = 0.0;
        wp = self.corr.wp
        n = len(wp)
        np= 0

        for i in range(n):
            if self.wp_obs[i,2] > 0.0:
                diff = (wp[i] - self.wp_obs[i,1])/self.wp_obs[i,2]
                chi2 += diff*diff
                np += 1

        #diff = (self.corr.wp[:] - self.wp_obs[:,1])/self.wp_obs[:,2]
        #chi2 = np.sum(diff**2)
        print('chi2', chi2, np)

        return chi2


    def write_corr_projected(self, index):
        """Write projected correlation function"""
        arg = self.domain + (index,)
        filename = '%s/corr_%s_%s_%s_%05d.txt' % (outdir, arg)
        rp = self.corr.rp
        wp = self.corr.wp
        with open(filename, 'w') as f:
            for i in range(len(rp)):
                f.write('%e %e\n' % (rp[i], wp[i]))


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
hod.set_coef([12, 0.0, 0.0, 0, 0.1, 0.0, 1.5, 0.0, 1.5, 0.0])


#
# Setup HOD parameter fitting
#
nbar= mock.NbarFitting(hod, nbar_obs, 0.6, 1.2)

#
# Setup output
#

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


def write_hod_params(h, index):
    """Write HOD parameters as a function of z
    Column 1: z
    Column 2: log10 M
    Column 3: sigma
    Column 4: log10 M1
    Column 5: alpha
    """
    filename = '%s/hod_%05d.txt' % (arg.outdir, index)
    f = open(filename, 'w')
    f.write('# c= ' + hod.get_coef().__repr__() + '\n')
    
    for z in np.arange(0.4, 1.2, 0.01):
        logMmin = h.logMmin(z)
        sigma   = h.sigma(z)
        logM1   = h.logM1(z)
        alpha   = h.alpha(z)
        f.write('%.4f %.6f %.6f %.6f %.6f\n' %
                (z, logMmin, sigma, logM1, alpha))
    f.close()
                
#flog = open('%s/fit_hod.log' % outdir, 'w', 1)
iter = 0

def cost_function(x):
    """Compute the chi^2
    Args:
        x[0]: M1/M_min
        x[1]: sigma
        x[2]: alpha

        domains: An arrain of domains
        domain is a tuble containing the survey region x redshift bin info
          domain[0] (str): 'w1' or 'w4'
          domain[1] (str): redshift_min
          domain[2] (str): redshift_max
    
    """

    hod[6] = x[0]  #log10 M1 = log10 (M_min + c_6 + c_7*(z - z_0)
    hod[4] = x[1] # sigma
    hod[8] = x[2] # alpha

    #print('eval %.3f %.3f' % (x[0], x[1]))
    #print(hod.coef)
    
    # Find best fitting logMmin(z) function
    nbar.fit()

    chi2 = 0
    for domain, d in data.items():
        chi2 += d.chi2(hod)


    return chi2

def logging_minimization(x):
    global iter
    iter += 1

    hod[6] = x[0]  #log10 M1 = log10 (M_min + c_6 + c_7*(z - z_0)
    hod[4] = x[1] # sigma
    hod[8] = x[2] # alpha

    # Find best fitting logMmin(z) function
    nbar.fit()

    chi2 = 0
    for domain, d in data.items():
        chi2 += d.chi2(hod)

    if flog:
        flog.write('%.3f %.4f %.4f %.4f\n' % (chi2, x[0], x[1], x[2]))
        for domain, d in data.items():
            d.write_corr_projected(iter)
        write_nbar_fitting(nbar, iter)
        write_hod_params(hod, iter)


x0 = [1.5, 0.1, 1.0]
ss = [0.2, 0.05, 0.1]
x = mock.minimise(cost_function, logging_minimization, x0, ss)

if mock.comm.rank == 0:
    print('minimum', x)

if flog:
    flog.close()
