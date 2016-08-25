import mockgallib._mockgallib as c


class Catalogues:
    """A collection of mock catalogues"""
    def __init__(self):
        self._cats = c._catalogues_alloc()

    def __len__(self):
        """Returns number of catalogues stored"""
        return c._catalogues_len(self._cats)

    def __repr__(self):
        return "A collection of %d catalogues" % len(self)

    def __getitem__(self, i):
        """ith catalogue as an np.array
        Column 0-2: x;
        Column 3:   redshift;
        Column 4:   vr;
        Column 5-6: ra, dec
        Column 7:   M;
        Column 8:   w;
        Column 9:   flag;
        """
        n = len(self)
        if i < 0 or i >= n:
            raise IndexError()
        return c._catalogues_catalogue(self._cats, i)

    def load_h5(self, filename):
        c._catalogues_load_h5(self._cats, filename)

    def generate_galaxies(self, hod, lightcones, z_min, z_max):
        """Generates n mock catalogues from n lightcones.
        Existing catalogues will be overwritten
        """
        c._catalogues_generate(self._cats, hod._hod, lightcones._lt,
                               z_min, z_max, 0)

    def generate_randoms(self, hod, lightcones, z_min, z_max):
        c._catalogues_generate(self._cats, hod._hod, lightcones._lt,
                               z_min, z_max, 1)

    def append(self, arr):
        """Append a catalogue
        
        Input: arr is an array with 3 or 4 columns
           Column 0 (a[:,0]): x
           Column 1 (a[:,1]): y
           Column 2 (a[:,2]): z
           Column 3 (a[:,3], optional): weight
           Column 4 (a[:,4], optional): ra
           Column 5 (a[:,5], optional): dec
        """

        c._catalogues_append(self._cats, arr)

    def ngal(self):
        """Return number of galaxies
        """
        return c._catalogues_ngal(self._cats)
