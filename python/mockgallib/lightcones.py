import mockgallib._mockgallib as c

class LightCones:
    """A collection of halo light cones"""
    def __init__(self):
        self._lt = c._lightcones_alloc()

    def __len__(self):
        return c._lightcones_len(self._lt)

    def __repr__(self):
        return "A collection of %d lightcones" % len(self)
        
    def load(self, filenames):
        """load([filename1, filename2]); load lightcones from files"""
        for filename in filenames:
            c._lightcones_load(self._lt, filename)

    def create_from_snapshots(self, snapshots, sky, remap):
        """create_from_snapshots(snapshots, sky, remap); create lightcones from halo snapshots"""
        cola_lightcones_create(snapshots._snps, sky._sky, remap._remap,
                               self._lt);
