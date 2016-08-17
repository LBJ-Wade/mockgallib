import mockgallib._mockgallib as c

"""Power Spectrum function P(k)

This module handles linear power spectrum read from file
"""

n = 0
k = None
P = None


def init(filename):
    """Initialise power module with filename

    Args:
        filename (str): filename of tabulated power spectrum k P

    Raises:
        FileNotFoundError: If unable to open file
    """

    c._power_init(filename)
    global n, k, P
    
    n = c._power_n()
    k = c._power_k()
    P = c._power_P()


def sigma(R):
    """Computes sigma, rms amplitude at spherical scale R"""

    return c._power_sigma(R)
