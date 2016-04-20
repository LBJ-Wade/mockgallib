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
        self.z_min= z_min
        self.z_max= z_max
        self._z = [ z for z in array_obs[:,0] if z_min <= z <= z_max ]
        self._nbar_hod= Nbar(ps, hod)
        self._nbar_obs= [ nbar for z,nbar in array_obs[:,[0,1]] if z_min <= z <= z_max ]
        self._fitting= c._nbar_fitting_alloc(ps._ps, hod._hod, array_obs,
                                             z_min, z_max)
        self.n = c._nbar_fitting_len(self._fitting)

    def fit(self):
        c._nbar_fitting_compute(self._fitting)

    def __repr__(self):
        return "NbarFitting: %d data in range %.3f <= z <= %.3f" % (self.n, self.z_min, self.z_max)

    def __len__(self):
        return self.n

    def __getitem__(self, i):
        return (self._z[i], self._nbar_obs[i], self._nbar_hod(self._z[i]))

    def nbar_hod(self):
        return [ self._nbar_hod(zz) for zz in self.z ]


    
        
    
