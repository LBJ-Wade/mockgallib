import mockgallib._mockgallib as c

class HaloMass:
    def __init__(self, filename):
        self._halo_mass = c._halo_mass_alloc(filename)
