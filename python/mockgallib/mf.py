import mockgallib._mockgallib as c

# ps= PowerSpectrum("filename")
# s= Sigma(ps)
# mf= MassFunction(s)

class MassFunction:
    def __init__(self, sigma):
        self.s = sigma
        self._mf = c._mf_alloc()
        self.delta_c = c._const_deltac();
        self.rho_m= c.cosmology_rhom()
        self.a = -1.0
        
    def __repr__(self):
        return "MassFunction"

    def dndlnM(self, M, z):
        sigma0 = self.s(M)
        a =1.0/(1.0 + z)
        if a != self.a:
            c._mf_set_redshift(self._mf, a)
            self.a = a
        
        nu= self.delta_c/sigma0  # ToDo: and growth factor D(z)
        return c._mf_f(self._mf, nu)*self.rho_m/M
