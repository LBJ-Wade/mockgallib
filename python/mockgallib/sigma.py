import mockgallib._mockgallib as c

class Sigma:
    def __init__(self, ps):
        self._s = c._sigma_alloc(ps._ps)
        self.n= c._sigma_n(self._s)
        self.M= c._sigma_M_array(self._s)
        self.sinv= c._sigma_sinv_array(self._s)
        self.M_range= c._sigma_M_range(self._s)
        self.sigma0_range= (self.sigma0(self.M_range[0]),
                          self.sigma0(self.M_range[1]))
        
    def __len__(self):
        return self.n

    def __getitem__(self, i):
        if i < 0 or i > self.n:
            return None
        return (self.M[i], 1.0/self.sinv[i])

    def __repr__(self):
        str = "M Sigma\n"
        for i in range(0,min(5, self.n)):
            str += "% .3e % .3e\n" % (self.M[i], 1.0/self.sinv[i])
        if self.n > 5:
            str += "......\n"
            str += "% .3e % .3e\n" % (self.M[self.n-1], 1.0/self.sinv[self.n-1])
        return str

    def mass(self, sigma0):
        return c._sigma_M(self._s, sigma0)

    def sigma0(self, M):
        return 1.0/c._sigma_0inv(self._s, M)
