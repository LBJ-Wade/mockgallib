import mockgallib._mockgallib as c

class Remap:
    def __init__(self, u, boxsize):
        self._remap = c._remap_alloc(u, boxsize)
        self.boxsize = c._remap_boxsize(self._remap)

    def __repr__(self):
        return "remaping to box " + self.boxsize().__repr__()
    


