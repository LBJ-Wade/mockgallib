import mockgallib._mockgallib as c


class PowerSpectrum:
    """Power spectrum function P(k) read from file

    ps = PowerSpectrum(filename)

    Args:
        filename (str): file name of k P(k) file; CAMB matterpower.dat format.
                        Assume that the aplitude is for redshift z = 0.

    Attributes:
        n: number of (k P) data pairs.
        k: a vector of k
        P: a vector of P
    """
    def __init__(self, filename):
        self._ps = c._power_alloc(filename)
        if self._ps is None:
            raise Exception("Unable to open file " + filename)
        self.n = c._power_n(self._ps)
        self.k = c._power_k(self._ps)
        self.P = c._power_P(self._ps)

    def __len__(self):
        return self.n

    def __getitem__(self, i):
        return (c._power_ki(self._ps, i), c._power_Pi(self._ps, i))

    def __repr__(self):
        str = "PowerSpectrum\n"
        for i in range(0, min(5, self.n)):
            str += "% .3e % .3e\n" % (self.k[i], self.P[i])
        if self.n > 5:
            str += "......\n"
            str += "% .3e % .3e\n" % (self.k[self.n-1], self.P[self.n-1])
        return str

    def sigma(self, r):
        """ps.sigma(R): Rms fluctuation of the density field; sigma8 for R=8.

        Args:
            R: Top hat smoothing length [1/h Mpc].

        Returns:
            sigma(R)
        """

        return c._power_sigma(self._ps, r)
