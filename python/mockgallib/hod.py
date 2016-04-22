import mockgallib._mockgallib as c

class Hod:
    """HOD parameters as a function of z parameterised by coefficients c"""
    def __init__(self, _hod=None):
        if _hod:
            self._hod= _hod
        else:
            self._hod= c._hod_alloc()

    def __repr__(self):
        return "HOD coef= " + self.coef().__repr__()

    def __len__(self):
        """number of parameterisation coefficients c[i]"""
        return len(self.coef())
    
    def __getitem__(self, i):
        """parameterisation coefficent c[i]"""
        return c._hod_get_coef(self._hod, i)

    def __setitem__(self, i, ci):
        """set c[i]"""
        c._hod_set_coef(self._hod, i, ci)

    def set(self, cs):
        """set all coefficients c"""
        for i, ci in enumerate(cs):
            c._hod_set_coef(self._hod, i, ci)
            
    def coef(self):
        return c._hod_get_coef_all(self._hod)

    def logMmin(self, z):
        """logMmin(z)"""
        x = z - 0.5
        co = self.coef()
        return co[0] + co[1]*x + co[2]*x**2 + co[3]*x**3

    def sigma(self, z):
        """sigma(z)"""
        x = z - 0.5
        co = self.coef()
        return co[4] + co[5]*x

    def logM1(self, z):
        """logM1(z)"""
        x = z - 0.5
        co = self.coef()
        return co[6] + co[7]*x

    def alpha(self, z):
        """alpha(z)"""
        x = z - 0.5
        co = self.coef()
        return co[8] + co[9]*x

    def ncen(self, M, z):
        """ncen(M, z): probability of having a central galaxy"""
        return c._hod_ncen(self._hod, M, z)

    def nsat(self, M, z):
        """nsat(M, z): mean number of satellite galaxies if the halo has a central"""
        return c._hod_nsat(self._hod, M, z)
