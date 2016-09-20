import mockgallib
import mockgallib._mockgallib as c

# ps= PowerSpectrum("filename")
# s= Sigma(ps)
# mf= MassFunction(s)

class MassFunction:
    def __init__(self):
        self._mf = c._mf_alloc()
        self.delta_c = mockgallib.const.delta_c;
        self.rho_m= mockgallib.cosmology.rhom()
        self.a = -1.0
        self._mfc = None
        self._mfc_a= None
        
    def __repr__(self):
        return "MassFunction"

    def dndlnM(self, M, z):
        sigma0 = self.s(M)
        a =1.0/(1.0 + z)
        if a != self.a:
            c._mf_set_redshift(self._mf, a)
            self.a = a
            self._D = growth_D(a)
        
        nu= self.delta_c/(self._D*sigma0)
        return c._mf_f(self._mf, nu)*self.rho_m/M

    def M(self, n, z):
        """
        args
            n  (float): halo number density M(>n)
            z  (float): redshift
        Returns
            M such that M(>n) = n
        """
        a = 1.0/(1.0 + z)
        if self._mfc_a != a:
            self._mfc = c._mf_cumulative_alloc(a)
        self._mfc_a = a

        return c._mf_cumulative_M(self._mfc, n)
