import unittest
import mockgallib as mock

class TestSlice(unittest.TestCase):
    def setUp(self):
        mock.set_loglevel(3)
        mock.cosmology_set(0.31)

        u =  [2, 2, 1, 1, 0, 0, 0, 1, 0]
        boxsize = 600.0
        self.remap = mock.Remap(u, boxsize)

        # sky
        ra_range  = [30.1, 38.9]
        dec_range = [-6.0, -4.1]
        z_range   = [0.38889, 1.21239]
        self.sky = mock.Sky(ra_range, dec_range, z_range)

        # slice
        self.slice = mock.Slice(self.remap, self.sky)

    def test_len(self):
        self.assertEqual(len(self.slice), 3)

    def test_boxsize(self):
        boxsize = self.slice.boxsize
        self.assertAlmostEqual(boxsize[0], 1800.0)
        self.assertAlmostEqual(boxsize[1], 447.213623046875)
        self.assertAlmostEqual(boxsize[2], 89.44271850585938)

if __name__ == '__main__':
    unittest.main()
