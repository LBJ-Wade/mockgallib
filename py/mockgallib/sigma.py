import mockgallib._mockgallib as c


class Sigma:
    """Sigma0(M): Density fluctuation at z = 0 on mass scale M.

    s = Sigma(s, M_min=1.0e10, M_max=1.0e16, n=1001)

    Args:
        ps (PowerSpectrum).
        [M_min (float)]: minimum mass for sigma(M).
        [M_max (float)]: maximum mass for sigma(M).
        [n (int)]: Number of data points for interpolation.

    Attributes:
        s.M:  array of mass [1/h Solar mass]
        s.sinv: array of 1/sigma(M)
        s(M): sigma0(M)
        s.inv(M): inverse function sigma0 = M(sigma0)
        s.M_range: (M_min, M_max)
        s.sigma0_range: (sigma0(M_max), sigma0(M_min))
    """
    
    def __init__(self, ps, *args, **kwargs):
        M_min = kwargs.get('M_min', 1.0e10)
        M_max = kwargs.get('M_max', 1.0e16)
        n     = kwargs.get('n', 1001)
    
        self._s = c._sigma_alloc(ps._ps, M_min, M_max, n)
        self.n= c._sigma_n(self._s)
        self.M= c._sigma_M_array(self._s)
        self.sinv= c._sigma_sinv_array(self._s)
        self.M_range= c._sigma_M_range(self._s)
        self.sigma0_range= (self(self.M_range[0]),
                            self(self.M_range[1]))

    def __len__(self):
        """number of precomputed data points; length of M and sinv"""
        return self.n

    def __getitem__(self, i):
        return (self.M[i], 1.0/self.sinv[i])

    def __repr__(self):
        str = "M Sigma\n"
        for i in range(0,min(5, self.n)):
            str += "% .3e % .3e\n" % (self.M[i], 1.0/self.sinv[i])
        if self.n > 5:
            str += "......\n"
            str += "% .3e % .3e\n" % (self.M[self.n-1], 1.0/self.sinv[self.n-1])
        return str

    def __call__(self, M):
        return 1.0/c._sigma_0inv(self._s, M)

    def inv(self, sigma0):
        return c._sigma_M(self._s, sigma0)

