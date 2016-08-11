import mockgallib._mockgallib as c


"""HOD class

HOD coeffficients:
        log10 M_min = c[0] + c[1]*x + c[2]*x**2 + c[3]*x**3
        sigma       = c[4] + c[5]*x
        M1    = (c[6] + c[7]*x)*M_min
        alpha       = c[8] + c[9]*x
"""


class Hod:
    """HOD parameters as a function of z parameterised by coefficients c
    Attributes:
        hod[i]    : HOD coeffcient c[i]
        logMmin(z): HOD parameter log_10 Mmin
        sigma(z)  : HOD parameter sigma
        logM1(z)  : HOD parameter log 10 M1
        alpha(z)  : HOD parameter alpha
        ncen(M, z): mean number of central galaxy per halo
        nsat(M, z): mean number of satellite galaxies per halo
    """
    def __init__(self, _hod=None):
        if _hod:
            self._hod = _hod
        else:
            self._hod = c._hod_alloc()

    def __repr__(self):
        return "HOD coef = " + self.get_coef().__repr__()

    def __len__(self):
        """number of parameterisation coefficients c[i]"""
        return len(self.get_coef())

    def __getitem__(self, i):
        """parameterisation coefficent c[i]"""
        return c._hod_get_coef(self._hod, i)

    def __setitem__(self, i, ci):
        """set c[i]"""
        c._hod_set_coef(self._hod, i, ci)

    @property
    def coef(self):
        return c._hod_get_coef_all(self._hod)
        
    def set_coef(self, cs):
        """set all coefficients c"""
        for i, ci in enumerate(cs):
            c._hod_set_coef(self._hod, i, ci)

    def get_coef(self):
        return c._hod_get_coef_all(self._hod)

    def logMmin(self, z):
        """logMmin(z)"""
        x = z - 0.5
        co = self.get_coef()
        return co[0] + co[1]*x + co[2]*x**2 + co[3]*x**3

    def sigma(self, z):
        """sigma(z)"""
        x = z - 0.5
        co = self.get_coef()
        return co[4] + co[5]*x

    def logM1(self, z):
        """M1(z)/M_min"""
        x = z - 0.5
        co = self.get_coef()
        return co[6] + co[7]*x

    def alpha(self, z):
        """alpha(z)"""
        x = z - 0.5
        co = self.get_coef()
        return co[8] + co[9]*x

    def ncen(self, M, z):
        """ncen(M, z): probability of having a central galaxy"""
        return c._hod_ncen(self._hod, M, z)

    def nsat(self, M, z):
        """nsat(M, z): mean number of satellite galaxies if the halo has
        a central
        """
        return c._hod_nsat(self._hod, M, z)
