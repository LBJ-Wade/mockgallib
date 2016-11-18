import mockgallib._mockgallib as c
import numpy
import h5py

class Hist2D:
    def __init__(self, *,
                 r_min=0.1, p_max=60.0, rp_nbin=24, mu_nbin=20):
        self.r_min = rp_min
        self.r_max = rp_max
        self.r_nbin= rp_nbin
        self.mu_nbin= mu_nbin
        
        self._hist2d = c._corr_multipole_hist2d_alloc(
            self.r_min, self.r_max, self.r_nbin, self.mu_nbin)


    def __getitem__(self, key):
        """
        Returns [ix,iy] component of the 2D histogram
        """
        a = c._corr_multipole_hist2d_as_array(self._hist2d)
        return a[key]

    def load(self, filename):
        f = h5py.File(filename, 'r')
        npairs = f['npairs'][()]
        a = f['rr'][:]
        f.close()

        c._corr_multipole_hist2d_set(self._hist2d, a, npairs)
    

class CorrelationFunctionMultipole:
    """A collection of correlation function

    Keyword arguments:
        r_min (=0.01): minimum r
        r_max (=60.0): maximum r
        nbin   (=100):  number of bins in r
        mu_nbin (=10):  number of bins in mu
        ra_min (= 0.0): minimum pairwise ra separation
        dec_min (=0.0):                  dec
        pair_correction (=None): pair wise correction filename

    Attributes:
        r_min, r_max, nbin, mu_nbin (same as above)
    """


    def __init__(self, *,
                 r_min=0.1, r_max=60.0, nbin=24,
                 mu_nbin=20,
                 ra_min=0.0, dec_min=0.0,
                 pair_correction= None):
        self._array = None
        
        self._corr = c._corr_multipole_alloc(r_min, r_max, nbin, mu_nbin)

        c._corr_multipole_set_radec_min(ra_min, dec_min)

        if pair_correction:
            c._corr_multipole_set_pair_correction(pair_correction)


    def __getitem__(self, key):
        self._array[key, :]


    def __len__(self):
        return self.nbin

    def __repr__(self):
        return ("Correlation function multipole\n")

    def compute_corr_multipole(self, cats_galaxies, cats_randoms, *,
                               direct=False):
        """Computes correlation function multipoles."""

        if direct:
            c._corr_multipole_compute_direct(cats_galaxies._cats,
                                             cats_randoms._cats,
                                             self._corr)
        else:
            c._corr_multipole_compute(cats_galaxies._cats,
                                      cats_randoms._cats,
                                      self._corr)

        self._array = numpy.transpose(c._corr_multipole_as_array(self._corr))
        return self._array

    def compute_corr_multipole_rr(self, cats_randoms, rr, *,
                                  direct=False):
        """Compute RR Hist2D
        Returns:
          rr->npairs
        """

        return c._corr_multipole_compute_rr(cats_randoms._cats, rr._hist2d)

    def compute_corr_multipole_pairs(self, cat_galaxies, cat_randoms,
                                     dd, dr, rr, *, direct=False):
        return c._corr_multipole_compute_all(cat_galaxies._cats,
                                    cat_randoms._cats,
                                    dd._hist2d, dr._hist2d, rr._hist2d,
                                    int(direct))

    def compute_corr_multipole_with_rr(self, cats_galaxies, cats_randoms, rr):
        """Compute correlation function multipoles with given RR
        """
        return c._corr_multipole_compute_with_rr(cats_galaxies._cats,
                                            cats_randoms._cats, rr._hist2d)


    def r_i(self, i):
        return c._corr_multipole_r_i(i)

    def xi0_i(self, i):
        return c._corr_multipole_xi0_i(i)

    def xi2_i(self, i):
        return c._corr_multipole_xi2_i(i)

    
