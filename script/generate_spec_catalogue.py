"""python3 spec.py

This script create a spectroscopic catalogue from:
  Target catalogue:     venice/mocks/w1/mock_w1_00001.txt
  Target Sampling Rate: tsr/out/w1/tsr_w1_00001.txt

Applies:
  Spectroscopic Sampling Rate (SSR)
    Redshift and quadrant dependent subsampling
  Redshift space distortion (RSD)
  Redshift error

Input:
  arg.data/vipers_ssr.txt
  ../venice/mocks/<reg>/mock_<reg>_<n>.txt
  ../tsr/out/<reg>/tsr_<reg>_<n>.txt

Output:
  std out

  Column  1:   ID       unique nonnegative integer
  Column  2-4: x,y,z    [1/h Mpc] comoving coordinate originated at observer
  Column  5:   vr       line-of-sight velocity [km/s]
  Column  6:   redshift
  Column  7:   z_cosmological 
  Column  8:   RA       [deg]
  Colum   9:   Dec      [deg]
  Column 10:   M_halo   host halo mass [1/h M_solar]
  Column 11:   r_sat    distance of galaxy from halo centre
  Column 12:   v_sat    line-of-sight virial velocity
  Column 13:   TSR
  Column 14:   SSR
  Column 15:   quadrant e.g W1P033Q1
"""

import math
import sys
import argparse
import json
import signal
import numpy as np
import mockgallib as mock

#
# Command-line options
#
parser = argparse.ArgumentParser()
parser.add_argument('n', help='index of mock')
parser.add_argument('--reg', default='w1', help='region w1/w4')
parser.add_argument('--param', default='param.json',
                    help='parameter json file')
parser.add_argument('--zmin', default='0.5', help='region w1/w4')
parser.add_argument('--zmax', default='1.2', help='region w1/w4')
#parser.add_argument('--tsr', help='generate mock catalogue',
#                    action="store_true")
#parser.add_argument('--ssr-quad', default='vipers_quad_ssr.txt',
#                    help='ssr quad file')

arg = parser.parse_args()

#
# Read parameter file
#
# print('Parameter file: %s' % arg.param)

with open(arg.param, 'r') as f:
    param = json.load(f)

omega_m = param['omega_m']
omega_l = 1.0 - omega_m
mock.set_loglevel(7)
mock.cosmology.set(omega_m)
mock.distance.init(1.21)
sys.stderr.write('Setting cosmology: omega_m= %.4f\n' % omega_m)

zmin = float(arg.zmin)
zmax = float(arg.zmax)
#rmin = mock.cosmology.compute_comoving_distance(1.0/(1.0 + arg.zmin))
#rmax = mock.cosmology.compute_comoving_distance(1.0/(1.0 + arg.zmax)))
#sys.stderr.write('rmin, max= %.4f %.4f\n' % (rmin, rmax))

def read_ssr_quad(filename):
    """read ssr file vipers_quad_ssr.txt and set to ssr_quad dictionary"""
    ssr_quad = {}
    
    for line in open(filename, 'r'):
        # line = 'W1P033 Q1 0.934210526316'
        v = line.rstrip().split(' ')
        ssr_quad[v[0] + v[1]] = float(v[2])

    return ssr_quad


def read_pointings(filename, reg):
    """
    Arg:
      filename: pointings_w1.txt
      reg: w1 or w4
    """
    reg = reg.upper()
    pointings = []

    for line in open(filename, 'r'):
        if line[0] == '#':
            continue
        if line[:2] == reg and line[3] != '9':
            v = line.split()
            pointings.append(v[0])

    return pointings

class SSRQuad:
    def __init__(self, reg):
        self.reg = reg.upper()
        self._ssr_quad = read_ssr_quad('vipers_quad_ssr.txt')
          ## dictionary 'W1P033Q1' -> 0.934210526316
        self._pointings = read_pointings('pointings_%s.txt' % reg, reg)
          ## index -> W1P033 

    def __call__(self, quad):
        """
          Arg:
              quad (str): W1-083-1
          Returns:
              SSR (float)
              name (str): P033Q1
        """
        v = quad.split('-')
        assert(len(v) == 3 and self.reg == v[0])
        n = int(v[1]) - 1
        name = self._pointings[n] + 'Q' + v[2]

        if name in self._ssr_quad:
            return (self._ssr_quad[name], name)
        
        ssr_default = 1
        return (ssr_default, name)
        

def main(reg, n):
    ssr_quad = SSRQuad(arg.reg)
    
    # TSR file
    ftsr = open('../tsr/out/%s/tsr_%s_%05d.txt' % (reg, reg, n), 'r')
    
    filename = '../venice/mocks/%s/mock_%s_%05d.txt' % (reg, reg, n)
    for line in open(filename):
        v = line.rstrip().split(' ')

        # Read TSR
        sline = ftsr.readline()
        while sline[0] == '#':
            sline = ftsr.readline()
        vtsr = sline.rstrip().split(' ')

        assert(v[0] == vtsr[0]) # check ID
        tsr = vtsr[1]
        
        # Compute RSD + redshift error
        x = float(v[1])
        y = float(v[2])
        z = float(v[3])
        vr = float(v[4])

        redshift = float(v[5])

        r = math.sqrt(x*x + y*y + z*z)
        a = 1.0/(1.0 + redshift)

        # aH(a) = a H_0 sqrt(Omega_m*a^-3 + Omega_lambda)
        aHinv = 1.0/(100.0*a*math.sqrt(omega_m/a**3 + omega_l))

        zsigma = 141.0*(1.0 + redshift)
        # s = x + 1/aH vr x^
        verr = zsigma*np.random.normal() # added redshift error
        #print(vr, verr)
        vr += verr
        #print(vr)
        x += vr*aHinv*x/r
        y += vr*aHinv*y/r
        z += vr*aHinv*z/r

        d = math.sqrt(x*x + y*y + z*z)
        zobs = mock.distance.redshift(d)
        if zobs < zmin or zobs > zmax:
            continue


        quad = v[11]

        ssr_z = 0.89910836
        if zobs > 0.5:
            ssr_z = 0.89910836/(1.0 + 1.02856938*(zobs - 0.5)**3.05048812)
        
        ssr_q, quad_name = ssr_quad(quad)
        ssr = ssr_q*ssr_z
        #print(ssr)
        #assert(ssr > 0.0)


        if np.random.uniform() < ssr:
            print('%s %e %e %e %s %e %s %s %s %s %s %s %s %e %s'
                % (v[0], x, y, z, v[4], zobs, v[5], v[6], v[7], v[8], v[9], v[10],
                   tsr, ssr, quad_name))
            
    ftsr.close()


main(arg.reg, int(arg.n))
