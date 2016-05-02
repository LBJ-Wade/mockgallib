import mockgallib._mockgallib as c

class Sky:
    def __init__(self, ra, dec, z):
        
        self.ra  = ra
        self.dec = dec
        self.z   = z
        
        self._sky= c._sky_alloc(ra[0], ra[1], dec[0], dec[1], z[0], z[1])

        self.bounding_box= c._sky_box(self._sky)
