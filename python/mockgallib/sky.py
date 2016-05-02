import mockgallib._mockgallib as c

class Sky:
    def __init__(self, ra, dec, z):
        
        self.ra  = ra
        self.dec = dec
        self.z   = z
        
        self._sky= c._sky_alloc(ra[0], ra[1], dec[0], dec[1], z[0], z[1])

        self.boxsize= c._sky_boxsize(self._sky)
        self.left= c._sky_left(self._sky)
        self.right= c._sky_right(self._sky)
        self.r = c._sky_r_range(self._sky)

    def __repr__(self):
        return ('Sky in ra= (%.5f, %.5f), dec=(%.5f, %.5f), z=(%.5f, %.5f)'
                % tuple(self.ra + self.dec + self.z))

    def compute_x(self, r, ra, dec):
        return c._sky_compute_x(self._sky, r, ra, dec)

    def compute_radec(self, x):
        return c._sky_compute_radec(self._sky, x[0], x[1], x[2])
