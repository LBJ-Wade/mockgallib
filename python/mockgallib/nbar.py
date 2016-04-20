import mockgallib._mockgallib as c

class Nbar:
    """Nbar(ps, hod): compute nbar(z) from power spectrum ps and HOD hod"""
    def __init__(self, ps, hod):
        self._ni= c._nbar_alloc(ps._ps, hod._hod)

    def __repr__(self):
        return "nbar integration object"

    def __call__(self, z):        
        """compute nbar(z)"""
        return c._nbar_compute(self._ni, z)

class NbarFitting:
    def __init__(self, ps, hod, array_obs, z_min, z_max):
        self._f= c._nbar_fitting_alloc(ps._ps, hod._hod, array_obs, z_min, z_max)
        
    
