import mockgallib._mockgallib as c
import numpy

class CorrelationFunction:
    """A collection of correlation function

    Keyword arguments:
        rp_min (=0.01): minimum rp
        rp_max (=60.0): maximum rp
        nbin   (=100):  number of bins in rp
        pi_max (=60.0): integration pi range
        pi_nbin (=10):  number of bins in pi

    Attributes:
        rp_min, rp_max, nbin, pi_max, pi_nbin (same as above)
        rp: array of rp
        xi: array of xi
        dxi: array of rms(xi) if more than 1 galaxy catalogue is given
    """


    def __init__(self, *args, **kwargs):
        self.rp_min = kwargs.get('rp_min', 0.01)
        self.rp_max = kwargs.get('rp_max', 60.0)
        self.nbin = kwargs.get('nbin', 101)
        self.pi_max = kwargs.get('pi_max', 60.0)
        self.pi_nbin = kwargs.get('pi_nbin', 10)
        self._array = None
        
        self._corr = c._corr_projected_alloc(
                          self.rp_min, self.rp_max, self.nbin,
                          self.pi_max, self.pi_nbin)


    def __getitem__(self, key):
        self._array[key, :]


    def __len__(self):
        return self.nbin


    def __repr__(self):
        return ("Projected correlation function\n"
                + "%d logarithmic bins in rp (%e, %e)\n"
                + "%d linear bins in pi (0, %e)") % (self.nbin,
                self.rp_min, self.rp_max, self.pi_nbin, self.pi_max)


    def compute_corr_projected(self, cats_galaxies, cats_randoms):
        """Computes projected correlation function."""
        c._corr_projected_compute(cats_galaxies._cats,
                                  cats_randoms._cats,
                                  self._corr)

        self._array = numpy.transpose(c._corr_as_array(self._corr))
        return self._array

    @property
    def rp(self):
        return c._corr_rp(self._corr)


    @property
    def wp(self):
        return c._corr_wp(self._corr)


    @property
    def dwp(self):
        return c._corr_dwp(self._corr)
