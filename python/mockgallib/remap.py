import mockgallib._mockgallib as c

class Remap:
    def __init__(self, u, boxsize):
        self._remap = c._remap_alloc(u, boxsize)

    def __repr__(self):
        return "remaping to box " + self.boxsize().__repr__()
    
    def boxsize(self):
        return c._remap_boxsize(self._remap)
