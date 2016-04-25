import mockgallib._mockgallib as c

class Sky:
    def __init__(ra_min, ra_max, dec_min, dec_max, z_min, z_max):
        self.ra  = [ra_min, ra_max]
        self.dec = [dec_min, dec_max]
        self.z   = [z_min, z_max]
        
        self._sky= c._sky_alloc(ra_min, ra_max, dec_min, dec_max, z_min, z_max)

        self.bounding_box= c._sky_box(self._sky)
