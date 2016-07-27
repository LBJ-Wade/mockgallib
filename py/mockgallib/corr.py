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
    """
    def __init__(self, *args, **kwargs):
        rp_min = kwargs.get('rp_min', 0.01)
        rp_max = kwargs.get('rp_max', 60.0)
        nbin = kwargs.get('nbin', 101)
        pi_max = kwargs.get('pi_max', 60.0)
        pi_nbin = kwargs.get('pi_nbin', 10) 
        
        self._corr = c._corr_projected_alloc(rp_min, rp_max, nbin,
                                             pi_max, pi_nbin)

    #def as_array(self):
        #return numpy.transpose(c._corr_as_array(self._corr))

    def __repr__(self):
        return "A correlation function"


    def compute_corr_projected(self, cats_galaxies, cats_randoms):
        """Computes projected correlation function."""
        c._corr_projected_compute(cats_galaxies._cats,
                                  cats_randoms._cats,
                                  self._corr)

        self._asarray = numpy.transpose(c._corr_as_array(self._corr))
        return self._asarray
