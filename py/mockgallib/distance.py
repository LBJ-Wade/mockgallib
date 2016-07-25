import mockgallib._mockgallib as c

def init(z_max):
    c._distance_init(z_max)


def redshift(d):
    return c._distance_redshift(d)
