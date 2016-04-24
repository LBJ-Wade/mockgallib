import mockgallib._mockgallib as c

class Catalogues:
    """A collection of mock catalogues"""
    def __init__(self):
        self._cats = c._catalogues_alloc()

    def __len__(self):
        return c._catalogues_len(self._cats)

    def __repr__(self):
        return "A collection of %d catalogues" % len(self)

    def __getitem__(self, i):
        """ith catalogue"""
        n = len(self)
        if i < 0 or i >= n:
            raise LookupError()
        return c._catalogues_catalogue(self._cats, i)

    def generate(self, hod, lightcones, z_min, z_max):
        c._catalogues_generate(self._cats, hod._hod, lightcones._lt,
                               z_min, z_max)
        
