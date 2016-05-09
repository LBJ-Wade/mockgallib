import mockgallib._mockgallib as c

class Random:
    def __init__(self):
        c._rng= py_random_alloc()

    
