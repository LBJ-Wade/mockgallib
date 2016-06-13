import mockgallib._mockgallib as c

class Distance:
    def __init__(self, z_max):
        self.z_max = z_max
        self.d_max = c._distance_init(z_max)

    def __repr__(self):
        return "distance redshift lookup z_max= %.4f" % self.z_max

    def redshift(self, d):
        z = c._distance_redshift(d)
        if z == -1.0:
            raise LookupError()
        return z
