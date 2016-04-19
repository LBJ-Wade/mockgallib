import mockgallib._mockgallib as c

class Hod:
    def __init__(self):
        self._hod= c._hod_alloc()
        self.coef = c._hod_get_coef(self._hod)
    
    def __repr__(self):
        return "HOD coef= " + self.coef.__repr__()

    def logMmin(self, z):
        x = z - 0.5
        return self.coef[0] + self.coef[1]*x + self.coef[2]*x**2 + self.coef[3]*x**3

    def sigma(self, z):
        x = z - 0.5
        return self.coef[4] + self.coef[5]*x

    def logM1(self, z):
        x = z - 0.5
        return self.coef[6] + self.coef[7]*x

    def alpha(self, z):
        x = z - 0.5
        return self.coef[8] + self.coef[9]*x

    def ncen(self, M, z):
        return c._hod_ncen(self._hod, M, z)

    def nsat(self, M, z):        
        return c._hod_nsat(self._hod, M, z)
