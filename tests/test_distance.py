import unittest
import numpy as np
import mockgallib as mock

class TestPower(unittest.TestCase):
    def setUp(self):
        self.z_max = 1.5
        mock.set_loglevel(3)
        mock.cosmology_set(0.31)
        self.distance = mock.Distance(self.z_max)

    def test_distance_redshift(self):
        d = np.arange(0, 3000, 100)
        z = [self.distance.redshift(dd) for dd in d]

        for dd,zz in zip(d,z):
            aa = 1.0/(1.0 + zz)
            dcheck = mock.cosmology_compute_comoving_distance(aa)
            self.assertAlmostEqual(dcheck, dd)


if __name__ == '__main__':
    unittest.main()
