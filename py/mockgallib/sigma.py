import mockgallib._mockgallib as c

"""Sigma function

This module handles sigma(M), rms fluctuation of density field

Attributes:
    n: number of interpolation points
    M: array of M
    sinv: array of 1/sigma0(M)
"""

n = 0
M = None
sinv = None


def init(*args, **kwargs):
    """Initilises sigma module

    Optional keyword arguments:
       M_min = 1.0e10: interpolation range minimum of sigma(M)
       M_max = 1.0e16: interpolation range maximum of sigma(M)
       n = 1001: number of interpolation points
    """

    print("sigma.init")
    
    M_min = kwargs.get('M_min', 1.0e10)
    M_max = kwargs.get('M_max', 1.0e16)
    n = kwargs.get('n', 1001)
    
    c._sigma_init(M_min, M_max, n)

    n = c._sigma_n()
    M = c._sigma_M_array()
    sinv = c._sigma_sinv_array()


