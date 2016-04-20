import mockgallib._mockgallib as c

class Hod:
    """HOD parameters as a function of z parameterised by coefficients c"""
    def __init__(self):
        self._hod= c._hod_alloc()
        self.coef = c._hod_get_coef(self._hod)
        
    def __repr__(self):
        self._update()
        return "HOD coef= " + self.coef.__repr__()

    def __len__(self):
        """number of parameterisation coefficients c[i]"""
        return len(self.coef)
    
    def __getitem__(self, i):
        """parameterisation coefficent c[i]"""
        self._update()
        return self.coef[i]

    def __setitem__(self, i, ci):
        """set c[i]"""
        c._hod_set_coef(self._hod, i, ci)
        self._update()

    def _update(self):
        self.coef = c._hod_get_coef(self._hod)

    def logMmin(self, z):
        """logMmin(z)"""
        x = z - 0.5
        self._update()
        return self.coef[0] + self.coef[1]*x + self.coef[2]*x**2 + self.coef[3]*x**3

    def sigma(self, z):
        """sigma(z)"""
        x = z - 0.5
        self._update()
        return self.coef[4] + self.coef[5]*x

    def logM1(self, z):
        """logM1(z)"""
        x = z - 0.5
        self._update()
        return self.coef[6] + self.coef[7]*x

    def alpha(self, z):
        """alpha(z)"""
        x = z - 0.5
        self._update()
        return self.coef[8] + self.coef[9]*x

    def ncen(self, M, z):
        """ncen(M, z): probability of having a central galaxy"""
        self._update()
        return c._hod_ncen(self._hod, M, z)

    def nsat(self, M, z):
        """nsat(M, z): mean number of satellite galaxies if the halo has a central"""
        self._update()
        return c._hod_nsat(self._hod, M, z)
