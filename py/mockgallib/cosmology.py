import mockgallib._mockgallib as c


def set(omega_m):
    """Set omega_m
    Args:
        omega_m: cosmological matter density Omega_m(z=0)
    """
    
    c._cosmology_set(omega_m)


def compute_comoving_distance(a):
    """Computes comoving distance correcspons to scale factor a"""
    
    return c._cosmology_compute_comoving_distance(a)


def rhom():
    """Mean matter density at z=0 in [M_solar (1/h Mpc)^-3]"""
    return c._cosmology_rhom()
