import mockgallib._mockgallib as c

class Slice:
    def __init__(self, remap, sky):
        self._slice = c._slice_alloc(remap._remap, sky._sky)
        self.boxsize = c._slice_boxsize(self._slice)

    def __len__(self):
        return c._slice_len(self._slice)


