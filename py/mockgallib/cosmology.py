import mockgallib._mockgallib as c


def set(omega_m):
    """Set omega_m
    Args:
        omega_m: cosmological matter density Omega_m(z=0)
    """
    
    c.cosmology_set(omega_m)

