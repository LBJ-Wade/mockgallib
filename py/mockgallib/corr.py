import mockgallib._mockgallib as c
import numpy
import h5py

class Hist2D:
    def __init__(self, *,
                 rp_min=0.1, rp_max=60.0, rp_nbin=24,
                 pi_max=60.0, pi_nbin=20):
        self.rp_min = rp_min
        self.rp_max = rp_max
        self.rp_nbin= rp_nbin
        self.pi_max= pi_max
        self.pi_nbin= pi_nbin
        
        self._hist2d = c._corr_projected_hist2d_alloc(
            self.rp_min, self.rp_max, self.rp_nbin,
            self.pi_max, self.pi_nbin)

        

    def __getitem__(self, key):
        """
        Returns [ix,iy] component of the 2D histogram
        """
        a = c._corr_projected_hist2d_as_array(self._hist2d)
        return a[key]

    def load(self, filename):
        f = h5py.File(filename, 'r')
        npairs = f['npairs'][()]
        a = f['rr'][:]
        f.close()

        c._corr_projected_hist2d_set(self._hist2d, a, npairs)
    

class CorrelationFunction:
    """A collection of correlation function

    Keyword arguments:
        rp_min (=0.01): minimum rp
        rp_max (=60.0): maximum rp
        nbin   (=100):  number of bins in rp
        pi_max (=60.0): integration pi range
        pi_nbin (=10):  number of bins in pi
        ra_min (= 0.0): minimum pairwise ra separation
        dec_min (=0.0):                  dec
        pair_correction (=None): pair wise correction filename

    Attributes:
        rp_min, rp_max, nbin, pi_max, pi_nbin (same as above)
        rp: array of rp
        xi: array of xi
        dxi: array of rms(xi) if more than 1 galaxy catalogue is given
    """


    def __init__(self, *,
                 rp_min=0.1, rp_max=60.0, nbin=24,
                 pi_max=60.0, pi_nbin=10,
                 ra_min=0.0, dec_min=0.0,
                 pair_correction= None):
        self.rp_min= rp_min
        self.rp_max= rp_max
        self.nbin = nbin
        self.pi_max= pi_max
        self.pi_nbin= pi_nbin
        self.ra_min= ra_min
        self.dec_min= dec_min
        self._array = None
        
        self._corr = c._corr_projected_alloc(
                          self.rp_min, self.rp_max, self.nbin,
                          self.pi_max, self.pi_nbin)

        c._corr_set_radec_min(self.ra_min, self.dec_min)

        if pair_correction:
            c._corr_set_pair_correction(pair_correction)


    def __getitem__(self, key):
        self._array[key, :]


    def __len__(self):
        return self.nbin

    def __repr__(self):
        return ("Projected correlation function\n"
                + "%d logarithmic bins in rp (%e, %e)\n"
                + "%d linear bins in pi (0, %e)") % (self.nbin,
                self.rp_min, self.rp_max, self.pi_nbin, self.pi_max)

    def compute_corr_projected(self, cats_galaxies, cats_randoms, *,
                               direct=False):
        """Computes projected correlation function."""

        if direct:
            c._corr_projected_compute_direct(cats_galaxies._cats,
                                             cats_randoms._cats,
                                             self._corr)
        else:
            c._corr_projected_compute(cats_galaxies._cats,
                                      cats_randoms._cats,
                                      self._corr)

        self._array = numpy.transpose(c._corr_as_array(self._corr))
        return self._array

    def compute_corr_projected_rr(self, cats_randoms, rr, *,
                                  direct=False):
        """Compute RR Hist2D
        Returns:
          rr->npairs
        """

        if direct:
            return c._corr_projected_compute_rr_direct(cats_randoms._cats, rr._hist2d)

        return c._corr_projected_compute_rr(cats_randoms._cats, rr._hist2d)

    def compute_corr_projected_pairs(self, cat_galaxies, cat_randoms,
                                     dd, dr, rr, *, direct=False):
        return c._corr_projected_compute_all(cat_galaxies._cats,
                                    cat_randoms._cats,
                                    dd._hist2d, dr._hist2d, rr._hist2d,
                                    int(direct))

    def compute_corr_projected_with_rr(self, cats_galaxies, cats_randoms, rr):
        """Compute projected correlation function with given RR
        """
        return c._corr_projected_compute_with_rr(cats_galaxies._cats,
                                            cats_randoms._cats, rr._hist2d)


    @property
    def rp(self):
        return c._corr_rp(self._corr)


    @property
    def wp(self):
        return c._corr_wp(self._corr)


    @property
    def dwp(self):
        return c._corr_dwp(self._corr)

    def rp_i(self, i):
        return c._corr_rp_i(i)

    def wp_i(self, i):
        return c._corr_wp_i(i)

    
