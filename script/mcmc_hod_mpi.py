import os
import sys
import argparse
import json
import signal
import numpy as np
import mockgallib as mock
import emcee

signal.signal(signal.SIGINT, signal.SIG_DFL) # stop with ctrl-c

def print0(*a):
    if mock.comm.rank == 0:
        print(*a)

#
# Command-line options
#
parser = argparse.ArgumentParser()
parser.add_argument('--param', default='param.json',
                    help='parameter json file')
parser.add_argument('--nmocks', default='1', help='number of mock catalogues')
parser.add_argument('--nrands', default='1', help='number of random catalogues')
parser.add_argument('--subsample', default='10000', help='number of data points')
parser.add_argument('--dir', default='.', help='base directory')
parser.add_argument('--outdir', default='log', help='output directory')
parser.add_argument('--loglevel', default='0', help='loglevel')
parser.add_argument('--test', default=False, action='store_true', help='compute one chi2 for given parameter')
parser.add_argument('--x0', default="1.0,0.15,1.0", help='initial parameters')
parser.add_argument('--nchain', default="10", help='length of MCMC chain')

arg = parser.parse_args()
outdir = arg.outdir
nsample = int(arg.subsample)

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
z_mins = []
z_maxs = []

def read_redshift_bins(redshift_bins):
    arr = []
    for zbin in redshift_bins:
        arr.append((float(zbin['zmin']), float(zbin['zmax'])))
        z_mins.append(float(zbin['zmin']))
        z_maxs.append(float(zbin['zmax']))

    return arr

redshift_bins = read_redshift_bins(param['redshift_bins'])
print0(redshift_bins)

z_min_all = min(z_mins)
z_max_all = max(z_maxs)

# nz
nbar_obs=  mock.array.loadtxt(arg.dir + '/' + param['nz'])

# sky
sky = {}
for reg in param['reg']:
    sky[reg['name']] = mock.Sky(reg['ra'], reg['dec'], [z_min_all, z_max_all])

flog = None
iter = 0
if mock.comm.rank == 0:
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if arg.test:
        for i in range(1000):
            filename = '%s/nz_%05d.txt' % (outdir, i)
            if not os.path.isfile(filename):
                break
        iter = i
        flog = open('%s/fit_hod.log' % outdir, 'a', 1)
    else:
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
        self.sky = sky[self.reg]
        
        print0('Create Data for ', domain) 
        
        # Load lightcones
        self.halo_lightcones = mock.LightCones()
        self.halo_lightcones.load_h5(
                           ['%s/halo_lightcone/%s/lightcone_%05d.h5' %
                            (arg.dir, domain[0], n + 1)
                            for n in range(int(arg.nmocks))
                            if n % mock.comm.n_nodes == mock.comm.rank
                           ])
        print0("len halo lightcone %d" % len(self.halo_lightcones))


        self.rand_lightcones = mock.LightCones()
        self.rand_lightcones.load_h5(
                           ['%s/rand_lightcone/%s/lightcone_%05d.h5' %
                            (arg.dir, domain[0], n + 1)
                            for n in range(int(arg.nrands))
                            if n % mock.comm.n_nodes == mock.comm.rank
                           ])
        print0("len rand lightcone %d" % len(self.rand_lightcones))

        # Catalogues will be generated from lightcones for given hod
        self.galaxy_catalogues = mock.Catalogues()
        self.random_catalogues = mock.Catalogues()

        self.corr = mock.CorrelationFunction(
            rp_min=0.1, rp_max=60.0, nbin=24, pi_max= 60.0, pi_nbin= 20,
            ra_min=0.001388889, dec_min=0.0375)


        # VIPERS projected correlation function
        self.wp_obs = mock.array.loadtxt(
            '%s/data/vipers/run4/0.1/corr_projected_%s_%s_%s.txt' % ((arg.dir,) +domain))[:,1]
        
        
        if mock.comm.rank == 0:
            self.covinv = np.load('%s/data/vipers/run4/0.1/covinv_%s_%s_%s.npy'\
 % ((arg.dir,) + domain))
            print0('Inverse covariance matrix', self.covinv.shape)
        else:
            self.covinv = None
        

    def generate_catalogues(self, hod):
        """Generate galaxy and random catalogues from given HOD"""
        print0('generate catalogues for %.2f - %.2f' % (self.z_min, self.z_max))
        self.galaxy_catalogues.generate_galaxies(hod, self.halo_lightcones,
                                                 self.sky,
                                                 self.z_min, self.z_max)
        
        self.random_catalogues.generate_randoms(hod, self.rand_lightcones,
                                                self.sky,
                                                self.z_min, self.z_max)

        if nsample > 0:
            self.galaxy_catalogues.subsample(nsample)
            self.random_catalogues.subsample(nsample)


    def chi2(self, hod):
        """Compute chi2 between HOD and observation 
        projected correlation functions wp's.
        Assumed that both wps are computed in the same coditions.
        """
        self.generate_catalogues(hod)

        print('compute corr projected %.2f - %.2f' % (self.z_min, self.z_max))
        self.corr.compute_corr_projected(self.galaxy_catalogues,
                                         self.random_catalogues)

        wp = self.corr.wp
        chi2 = 0.0;

        if mock.comm.rank == 0:
             chi2 = np.dot(wp - self.wp_obs, self.covinv.dot(wp - self.wp_obs))

        print0('chi2 %f / %d' % (chi2, len(wp)))
        self.x2 = chi2
        
        return chi2


    def write_corr_projected(self, index):
        """Write projected correlation function"""
        print0('write projected correlation function')
        arg = (outdir,) + self.domain + (index,)
        filename = '%s/corr_%s_%s_%s_%05d.txt' % arg
        rp = self.corr.rp
        wp = self.corr.wp
        with open(filename, 'w') as f:
            for i in range(len(rp)):
                f.write('%e %e\n' % (rp[i], wp[i]))


domains = []
for zbin in redshift_bins:
    domains.append(('w1', zbin[0], zbin[1]))
    domains.append(('w4', zbin[0], zbin[1]))
data = {}

for domain in domains:
    data[domain] = Data(domain)

#
# Set HOD parameters (initial guess)
#
z0 = 0.5
hod = mock.Hod(z0=z0)
hod.set_coef([12, 0.0, 0.0, 0, 0.1, 0.0, 1.5, 0.0, 1.5, 0.0, 3.0])


#
# Setup HOD parameter fitting
#
nbar= mock.NbarFitting(hod, nbar_obs, 0.5, 1.2)
nbar1 = mock.Nbar(hod)

#
# Setup output
#

def write_nbar_fitting(nz, nbar1, index):
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
        print('write_nbar_fitting', hod)
        for i in range(len(nz)):
            ncen = nbar1.ncen(nz.z[i])
            nsat = nbar1.nsat(nz.z[i])
            f.write('%e %e %e %e %e\n' %
                    (nz.z[i], nz.nbar_obs[i], nz.nbar_hod[i], ncen, nsat))


def write_hod_params(h, index):
    """Write HOD parameters as a function of z
    Column 1: z
    Column 2: log10 M
    Column 3: sigma
    Column 4: log10 M1
    Column 5: alpha
    """
    filename = '%s/hod_%05d.txt' % (outdir, index)
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

def lnprob_function(x):
    """Compute log probability - 0.5*chi^2
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

    #sys.stderr.write('This is lnprof_function %d\n' % mock.comm.rank)
    print('This is lnprof_function %d' % mock.comm.rank)
    print('x= ', x, mock.comm.rank)

    print(type(x).__module__)
    if type(x).__module__ == 'numpy':
        x = x.tolist()
    mock.callback.sync(x) 
    print('callback sync done %d' % mock.comm.rank)

    if x[0] < 0 or x[0] > 4 or x[1] < 0 or x[2] < 0:
        print('out of range return', x)
        return -np.inf

    hod[6] = x[0]  #log10 M1 = log10 (M_min + c_6 + c_7*(z - z_0)
    #hod[7] = x[1] 
    hod[4] = x[1] # sigma
    hod[8] = x[2] # alpha

    print0('eval', x)
    
    # Find best fitting logMmin(z) function
    nbar.fit()
    nbar_fitting_iter = nbar.iter
    nbar_fitting_iter = mock.comm.bcast_int(nbar_fitting_iter)

    print0('nbar.fit', nbar.iter)
    if nbar_fitting_iter > 900:
        return -1000.0

    chi2 = 0

    for domain in domains:
        d = data[domain]
        chi2 += d.chi2(hod)

    print0('cost_function ', chi2)

    #print('returning from lnprob')
    return -0.5*chi2

def print_array(format, v):
    str=''
    for x in v:
        str += (format % x) + ' '
    return str
    
def logging_minimization(x):
    global iter
    iter += 1

    hod[6] = x[0] #log10 M1 = log10 (M_min + c_6 + c_7*(z - z_0)
    hod[7] = x[1] #log10 M1 = log10 (M_min + c_6 + c_7*(z - z_0)
    hod[4] = x[2] # sigma
    hod[8] = x[3] # alpha

    # Find best fitting logMmin(z) function
    nbar.fit()
    nbar_fitting_iter = nbar.iter
    nbar_fitting_iter = mock.comm.bcast_int(nbar_fitting_iter)

    print('logging nbar.fit %d' % nbar_fitting_iter, hod)
    if nbar_iter > 900:
        return -1000

    chi2_total = 0.0
    chi2_each = []

    for domain in domains:
        d = data[domain]
        chi2 = d.chi2(hod)
        chi2_each.append(chi2)
        chi2_total += chi2

    if flog:
        log_str = ('%.3f' % chi2_total) + ' '
        log_str += print_array('%.4f', x)
        log_str += print_array('%.4f', chi2_each)
        flog.write(log_str + '\n')

        for domain, d in data.items():
            d.write_corr_projected(iter)
        write_nbar_fitting(nbar, nbar1, iter)
        write_hod_params(hod, iter)

    return chi2_total
#
#x0 = [ float(x) for x in arg.x0.split(',') ]
#print(x0)

x0 = [1.5, 0.15, 1.0]
ss = [0.5, 0.05, 0.5]

#if arg.test:
#    x = x0
#    print('test run at x= ', x0)
#else:
#    x = mock.minimise(cost_function, logging_minimization, x0, ss)
#    print0('minimum', x)

#chi2_final= logging_minimization(x)
#print0(hod)

ndim = len(x0)
print0('ndim= ', ndim)
nwalkers = 2*ndim

def init_x0(n):
    x = [np.random.normal(loc=x0[i], scale=ss[i]) for i in range(n)]
    return np.array(x)

# nwalker initial points
p0 = [ init_x0(ndim) for i in range(nwalkers) ]

if mock.comm.rank == 0:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_function)
    print('run mcmc', arg.nchain)
    
    sampler.run_mcmc(p0, int(arg.nchain))

    print('mcmc done')
    mock.callback.release()
    
    #print(sampler.flatchain[:])
    chain = sampler.chain[:, :, :].reshape((-1, ndim))
    np.savetxt('chain.txt', chain) #sampler.flatchain[:])
else:
    print('standby')
    mock.callback.standby(lnprob_function)


#
# Write hod.json
#

# if mock.comm.rank == 0:
#     domains_json = []
#     for domain in domains:
#         d = data[domain]
#         domains_json.append({
#             'reg': domain[0], 'zmin': domain[1], 'zmax': domain[2],
#             'chi2': d.x2})

#     hod_json = {
#         'z0': z0,
#         'c': hod.coef,
#         'domains': domains_json,
#         'chi2': chi2_final
#     }

#     with open('fit_hod.json', 'w') as f:
#         json.dump(hod_json, f)

# if flog:
#     flog.close()
